'''
General multi-level, multi-laser Lindblad master-equation engine.

This module builds the equations of motion

    d rho/dt = -i [H(t), rho] + sum_k D[C_k] rho - 1/2 {K, rho}

for any set of atomic sub-states coupled by any number of lasers, and solves
them with either an exact matrix-exponential propagator (time-independent
problems), the Rust backend (time-dependent problems and ensembles, run in
parallel), or a SciPy fallback.

Units: angular frequencies and rates in Grad/s (= rad/ns), times in ns,
wavelengths in m, intensities in mW/mm^2, velocities in m/s.

The Hamiltonian is assembled as a Hermitian matrix in the rotating frame of
the lasers, H = H0 + sum_k (V_k + V_k^dagger), where V_k[u, l] couples lower
state l to upper state u of laser k. Hermiticity (and therefore real
populations and conservation of probability) holds for any complex
polarisation, envelope or phase. The superoperator is built in row-major
vectorisation, vec(rho)[i*n + j] = rho_ij, matching ``index(i, j, n)``.

Example:
    >>> import LASED as las, numpy as np
    >>> g = las.generateSubStates(label_from=1, w=0, L=0, S=0)
    >>> e = las.generateSubStates(label_from=2, w=las.angularFreq(422.8e-9), L=1, S=0)
    >>> laser = las.Laser(lower=g, upper=e, wavelength=422.8e-9, intensity=100, tau=4.6,
    ...                   polarisation=las.elliptical(angle=0, ellipticity=0.3))
    >>> system = las.AtomicSystem(g + e, lasers=[laser], decays=[las.Decay(upper=e, lower=g, tau=4.6)])
    >>> rho = system.evolve(np.linspace(0, 50, 501))   # (501, 4, 4) complex array
'''

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Union
from collections import deque
import warnings

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

from LASED.constants import PI, C
from LASED.half_rabi_freq import coupling, halfRabiFreq
from LASED.polarisation import asPolarisation
from LASED import _backend

ComplexFn = Callable[[np.ndarray], np.ndarray]


# --------------------------------------------------------------------------------------
# Physical components
# --------------------------------------------------------------------------------------

@dataclass
class Laser:
    """A laser coupling a set of lower sub-states to a set of upper sub-states.

    The coupling between sub-states u (upper) and l (lower) is
    ``V[u, l] = (rabi/2) * A(t) * sum_q a_q * coupling(u, l, q)`` with V[l, u] = conj(V[u, l]).

    Attributes:
        lower (list of State): Lower sub-states coupled by this laser.
        upper (list of State): Upper sub-states coupled by this laser.
        wavelength (float): Vacuum wavelength of the laser in metres. Sets the laser frequency
            (and the Doppler shift); detunings are measured from it.
        intensity (float): Intensity in mW/mm^2. Requires ``tau`` to obtain the dipole strength.
        tau (float): Lifetime (ns) of the upper -> lower transition used to convert intensity to Rabi frequency.
        rabi (float): Alternatively, the reduced Rabi frequency Omega_0 in Grad/s. Each sub-transition
            has Rabi frequency Omega_0 * sum_q a_q coupling(u, l, q).
        polarisation (dict): Spherical amplitudes ``{q: a_q}`` (see ``LASED.polarisation``). May be complex.
        detuning (float): Static detuning of the laser in Grad/s (added to the laser frequency).
        envelope (callable or Pulse): Complex amplitude A(t) (t in ns) multiplying the coupling, e.g. a
            pulse shape times exp(i phase(t)). ``None`` means A = 1 (continuous wave).
        detuning_t (callable): Time-dependent detuning delta(t) in Grad/s added to ``detuning`` (chirps, ARP).
        direction (sequence of 3 floats): Unit propagation direction, used for Doppler shifts.
        linewidth (float): Lorentzian laser linewidth (FWHM) in Grad/s, modelled as white frequency
            noise (dephasing of the coherences this laser drives).
        name (str): Optional label.
    """
    lower: Sequence
    upper: Sequence
    wavelength: float
    intensity: Optional[float] = None
    tau: Optional[float] = None
    rabi: Optional[float] = None
    polarisation: Union[dict, int, Sequence] = field(default_factory=lambda: {0: 1+0j})
    detuning: float = 0.0
    envelope: Optional[ComplexFn] = None
    detuning_t: Optional[Callable] = None
    direction: Sequence = (0.0, 0.0, 1.0)
    linewidth: float = 0.0
    name: str = ""

    @property
    def angular_frequency(self):
        """Laser angular frequency in Grad/s (without detuning)."""
        return 2*PI*C/self.wavelength*1e-9

    @property
    def wavenumber(self):
        """Wavenumber k = 2 pi / lambda in rad/m."""
        return 2*PI/self.wavelength

    def halfRabi(self):
        """Reduced half-Rabi frequency Omega_0/2 in Grad/s."""
        if self.rabi is not None:
            return self.rabi/2
        if self.intensity is None or self.tau is None:
            raise ValueError(f"Laser {self.name!r}: give either rabi, or intensity together with tau.")
        return halfRabiFreq(self.intensity, self.tau, self.wavelength)

    def dopplerShift(self, velocity):
        """Detuning (Grad/s) seen by an atom moving with ``velocity`` (m/s, scalar along z or 3-vector)."""
        v = np.asarray(velocity, dtype=float)
        k_hat = np.asarray(self.direction, dtype=float)
        k_hat = k_hat/np.linalg.norm(k_hat)
        v_par = float(v)*k_hat[2] if v.ndim == 0 else float(np.dot(k_hat, v))
        return -self.wavenumber*v_par*1e-9

    @property
    def time_dependent(self):
        return self.envelope is not None or self.detuning_t is not None


@dataclass
class Decay:
    """Spontaneous emission from ``upper`` to ``lower`` sub-states.

    Branching between sub-states follows the squared dipole couplings, normalised for each upper
    sub-state so that its total decay rate into ``lower`` is ``branching/tau``. Coherences between upper
    sub-states are transferred to the lower states (including vertical coherences between different F
    levels). Implemented with Lindblad collapse operators, so the evolution stays physical.

    Attributes:
        upper (list of State): Decaying sub-states.
        lower (list of State): Sub-states decayed into.
        tau (float): Lifetime in ns.
        branching (float): Fraction of the decay going into ``lower`` (default 1).
    """
    upper: Sequence
    lower: Sequence
    tau: float
    branching: float = 1.0


@dataclass
class Loss:
    """Population loss out of the modelled system (e.g. decay to unmodelled levels, ionisation).

    Attributes:
        states (list of State): Sub-states that lose population.
        tau (float): Lifetime in ns of the loss channel.
    """
    states: Sequence
    tau: float


@dataclass
class Dephasing:
    """Pure dephasing of the given sub-states at ``rate`` (Grad/s): collapse operator sqrt(rate)*P_states."""
    states: Sequence
    rate: float


@dataclass
class Ensemble:
    """A weighted set of atoms that differ in velocity and/or laser intensity.

    The result of an evolution is ``sum_m weights[m] * rho_m(t)``.

    Attributes:
        weights (array, shape (M,)): Weight of each member.
        detuning (array, shape (M, n_lasers)): Extra detuning (Grad/s) of each laser for each member.
        rabi_scale (array, shape (M, n_lasers)): Factor multiplying each laser's Rabi frequency.
    """
    weights: np.ndarray
    detuning: np.ndarray
    rabi_scale: np.ndarray

    def __post_init__(self):
        self.weights = np.asarray(self.weights, dtype=float)
        self.detuning = np.atleast_2d(np.asarray(self.detuning, dtype=float))
        self.rabi_scale = np.atleast_2d(np.asarray(self.rabi_scale, dtype=float))
        if self.detuning.shape != self.rabi_scale.shape or self.detuning.shape[0] != len(self.weights):
            raise ValueError("Ensemble arrays must have shapes (M,), (M, n_lasers), (M, n_lasers).")

    @property
    def size(self):
        return len(self.weights)

    @staticmethod
    def single(n_lasers):
        return Ensemble(np.ones(1), np.zeros((1, n_lasers)), np.ones((1, n_lasers)))

    def __mul__(self, other):
        """Outer product of two ensembles (e.g. Doppler x beam profile)."""
        M1, M2 = self.size, other.size
        w = np.outer(self.weights, other.weights).ravel()
        det = (self.detuning[:, None, :] + other.detuning[None, :, :]).reshape(M1*M2, -1)
        scale = (self.rabi_scale[:, None, :]*other.rabi_scale[None, :, :]).reshape(M1*M2, -1)
        return Ensemble(w, det, scale)


def dopplerEnsemble(lasers, temperature=None, mass=None, sigma_v=None, n=41, n_sigma=4.0, direction=(0, 0, 1)):
    """Maxwell-Boltzmann velocity ensemble along ``direction`` (1D), discretised on n points.

    Parameters:
        lasers (list of Laser): The lasers of the system (for their Doppler shifts).
        temperature (float): Temperature in kelvin (with ``mass`` in kg), or give ``sigma_v`` in m/s.
        mass (float): Atomic mass in kg.
        sigma_v (float): 1D velocity standard deviation in m/s.
        n (int): Number of velocity classes.
        n_sigma (float): Velocity range is +- n_sigma * sigma_v.
        direction (sequence): Direction of the atomic motion.

    Returns:
        Ensemble
    """
    if sigma_v is None:
        if temperature is None or mass is None:
            raise ValueError("Give sigma_v or (temperature, mass).")
        sigma_v = np.sqrt(1.380649e-23*temperature/mass)
    v = np.linspace(-n_sigma*sigma_v, n_sigma*sigma_v, n)
    w = np.exp(-0.5*(v/sigma_v)**2)
    w /= w.sum()
    d = np.asarray(direction, dtype=float)
    d = d/np.linalg.norm(d)
    det = np.array([[laser.dopplerShift(vi*d) for laser in lasers] for vi in v])
    return Ensemble(w, det, np.ones_like(det))


def gaussianBeamEnsemble(n_lasers, waist_ratio=None, n=20, r_max=2.0, lasers_affected=None):
    """Average over the transverse intensity profile of Gaussian beams.

    Atoms are uniformly spread over a disc of radius ``r_max`` (in units of the 1/e^2 intensity radius w).
    The Rabi frequency of each affected laser is scaled by exp(-r^2/w^2) (i.e. intensity exp(-2r^2/w^2)),
    where the beam's peak Rabi frequency is the one set on the Laser.

    Parameters:
        n_lasers (int): Number of lasers in the system.
        waist_ratio (sequence): Optional per-laser ratio of beam waist to the reference waist (default all 1).
        n (int): Number of radial rings.
        r_max (float): Outer radius in units of the waist.
        lasers_affected (sequence of int): Indices of lasers with a Gaussian profile (default all).

    Returns:
        Ensemble
    """
    edges = np.linspace(0, r_max, n+1)
    r = 0.5*(edges[1:] + edges[:-1])
    w = edges[1:]**2 - edges[:-1]**2
    w /= w.sum()
    ratio = np.ones(n_lasers) if waist_ratio is None else np.asarray(waist_ratio, dtype=float)
    affected = range(n_lasers) if lasers_affected is None else lasers_affected
    scale = np.ones((n, n_lasers))
    for k in affected:
        scale[:, k] = np.exp(-(r/ratio[k])**2)
    return Ensemble(w, np.zeros((n, n_lasers)), scale)


# --------------------------------------------------------------------------------------
# Superoperator helpers (row-major vectorisation)
# --------------------------------------------------------------------------------------

def _spre(A, I):
    return sp.kron(A, I, format="csr")


def _spost(A, I):
    return sp.kron(I, A.T, format="csr")


def commutatorSuperoperator(H):
    """Returns the sparse superoperator of rho -> -i[H, rho] (row-major vec)."""
    H = sp.csr_matrix(H)
    I = sp.identity(H.shape[0], dtype=complex, format="csr")
    return (-1j*(_spre(H, I) - _spost(H, I))).tocsr()


def dissipatorSuperoperator(Cop):
    """Returns the sparse superoperator of rho -> C rho C^dag - 1/2 {C^dag C, rho} (row-major vec)."""
    Cop = sp.csr_matrix(Cop)
    I = sp.identity(Cop.shape[0], dtype=complex, format="csr")
    CdC = (Cop.conj().T @ Cop).tocsr()
    return (sp.kron(Cop, Cop.conj(), format="csr") - 0.5*_spre(CdC, I) - 0.5*_spost(CdC, I)).tocsr()


def anticommutatorSuperoperator(K):
    """Returns the sparse superoperator of rho -> -1/2 {K, rho} (row-major vec)."""
    K = sp.csr_matrix(K)
    I = sp.identity(K.shape[0], dtype=complex, format="csr")
    return (-0.5*(_spre(K, I) + _spost(K, I))).tocsr()


def reachableSubspace(L0, ops, y0):
    """Indices of vec(rho) that can become non-zero when evolving y0 under L0 and ops.

    The remaining elements are exactly zero for all times, so the problem can be solved in the
    reachable subspace only (e.g. with pi light and a diagonal rho0 only Delta m = 0 elements couple).

    Returns:
        ndarray: Sorted integer indices.
    """
    pattern = abs(sp.csr_matrix(L0))
    for op in ops:
        pattern = pattern + abs(sp.csr_matrix(op))
    pattern = (pattern != 0).astype(np.float64).tocsr()
    reached = np.abs(np.asarray(y0)) > 0
    while True:
        new = reached | ((pattern @ reached.astype(np.float64)) > 0)
        if np.array_equal(new, reached):
            return np.nonzero(reached)[0]
        reached = new


# --------------------------------------------------------------------------------------
# The system
# --------------------------------------------------------------------------------------

@dataclass
class Model:
    """The assembled equations of motion d vec(rho)/dt = [L0 + sum_j c_j(t) M_j] vec(rho).

    Attributes:
        n (int): Number of sub-states.
        L0 (scipy.sparse.csr_matrix): Static Liouvillian superoperator (n^2 x n^2).
        ops (list): Superoperators M_j multiplying time-dependent or ensemble-dependent coefficients.
        op_info (list of tuple): (laser index, kind) for each op, kind in {"re", "im", "detuning"}.
        H0 (ndarray): Static diagonal Hamiltonian in the rotating frame (n x n, Grad/s).
        V (list of ndarray): Coupling matrix V_k for each laser (n x n, upper-lower block only).
        D (list of ndarray): Diagonal detuning generator for each laser: dH/d(delta_k).
        collapse (list of ndarray): Collapse operators.
        K (ndarray): Loss operator (anti-Hermitian part of the effective Hamiltonian is -i K/2).
    """
    n: int
    L0: sp.csr_matrix
    ops: list
    op_info: list
    H0: np.ndarray
    V: list
    D: list
    collapse: list
    K: np.ndarray


class AtomicSystem:
    """A set of atomic sub-states driven by lasers, with spontaneous decay, losses and dephasing.

    Parameters:
        states (list of State): All sub-states. Labels must be 1..n (the matrix index is label-1).
        lasers (list of Laser): Lasers coupling the states.
        decays (list of Decay): Spontaneous decay channels.
        losses (list of Loss): Loss channels out of the system.
        dephasings (list of Dephasing): Pure dephasing channels.

    Example:
        >>> system = las.AtomicSystem(G + E, lasers=[laser], decays=[las.Decay(E, G, tau=26.2)])
        >>> rho_t = system.evolve(time)               # (len(time), n, n)
        >>> pops = system.populations(rho_t)          # {label: array over time}
    """

    def __init__(self, states, lasers=(), decays=(), losses=(), dephasings=()):
        self.states = sorted(states, key=lambda s: s.label)
        labels = [s.label for s in self.states]
        if labels != list(range(1, len(labels)+1)):
            raise ValueError(f"State labels must be exactly 1..n with no gaps or duplicates, got {labels}.")
        self.lasers = list(lasers)
        self.decays = list(decays)
        self.losses = list(losses)
        self.dephasings = list(dephasings)

    @property
    def n(self):
        return len(self.states)

    def _idx(self, state):
        return state.label - 1

    # ---------------------------------------------------------------- rotating frame
    def _frame(self):
        """Finds the rotating frame. Returns (f0, d): frame frequencies at zero detuning (n,)
        and d[k, i] = d f_i / d omega_k (the number of photons of laser k "absorbed" to reach state i)."""
        n, nl = self.n, len(self.lasers)
        adj = [[] for _ in range(n)]
        for k, laser in enumerate(self.lasers):
            for l in laser.lower:
                for u in laser.upper:
                    adj[self._idx(l)].append((self._idx(u), k, +1))
                    adj[self._idx(u)].append((self._idx(l), k, -1))
        d = np.zeros((nl, n), dtype=int)
        visited = np.zeros(n, dtype=bool)
        root_of = np.arange(n)
        for start in range(n):
            if visited[start]:
                continue
            visited[start] = True
            queue = deque([start])
            while queue:
                i = queue.popleft()
                for j, k, sgn in adj[i]:
                    expected = d[:, i].copy()
                    expected[k] += sgn
                    if not visited[j]:
                        visited[j] = True
                        root_of[j] = start
                        d[:, j] = expected
                        queue.append(j)
                    elif not np.array_equal(d[:, j], expected):
                        raise ValueError(
                            "The lasers form a closed loop, so no rotating frame removes all time dependence. "
                            "If two Laser objects describe the same beam, merge them into one Laser with a "
                            "combined polarisation.")
        omegas = np.array([laser.angular_frequency for laser in self.lasers])
        w = np.array([s.w for s in self.states], dtype=float)
        f0 = w[root_of] + (d.T @ omegas if nl else 0.0)
        return f0, d

    # ---------------------------------------------------------------- assembly
    def build(self):
        """Assembles the Hamiltonian, collapse operators and Liouvillian superoperators.

        Returns:
            Model
        """
        n = self.n
        f0, d = self._frame()
        w = np.array([s.w for s in self.states], dtype=float)
        H0 = np.diag(w - f0).astype(complex)

        V, D = [], []
        for k, laser in enumerate(self.lasers):
            pol = asPolarisation(laser.polarisation)
            half = laser.halfRabi()
            Vk = np.zeros((n, n), dtype=complex)
            for u in laser.upper:
                for l in laser.lower:
                    amp = sum(a*coupling(u, l, q) for q, a in pol.items() if a != 0)
                    Vk[self._idx(u), self._idx(l)] += half*amp
            V.append(Vk)
            # H_ii = w_i - f_i and f_i depends on (omega_k + delta_k) through d[k, i]
            D.append(np.diag(-d[k].astype(float)).astype(complex))

        collapse = []
        K = np.zeros((n, n), dtype=complex)
        for decay in self.decays:
            rate = decay.branching/decay.tau
            norms = {}
            for e in decay.upper:
                norms[e.label] = sum(coupling(e, g, q)**2 for g in decay.lower for q in (-1, 0, 1))
            for q in (-1, 0, 1):
                Cq = np.zeros((n, n), dtype=complex)
                for e in decay.upper:
                    if norms[e.label] == 0:
                        continue
                    for g in decay.lower:
                        c = coupling(e, g, q)
                        if c != 0:
                            Cq[self._idx(g), self._idx(e)] = np.sqrt(rate/norms[e.label])*c
                if np.any(Cq):
                    collapse.append(Cq)
            for e in decay.upper:  # No allowed channel into `lower`: population simply leaves
                if norms[e.label] == 0:
                    K[self._idx(e), self._idx(e)] += rate
        for loss in self.losses:
            for s in loss.states:
                K[self._idx(s), self._idx(s)] += 1/loss.tau
        for deph in self.dephasings:
            P = np.zeros((n, n), dtype=complex)
            for s in deph.states:
                P[self._idx(s), self._idx(s)] = 1
            collapse.append(np.sqrt(deph.rate)*P)
        for k, laser in enumerate(self.lasers):
            if laser.linewidth:
                # White frequency noise on laser k: collapse operator sqrt(Gamma) dH/d(delta_k)
                collapse.append(np.sqrt(laser.linewidth)*D[k])

        L0 = commutatorSuperoperator(H0)
        for Cop in collapse:
            L0 = L0 + dissipatorSuperoperator(Cop)
        if np.any(K):
            L0 = L0 + anticommutatorSuperoperator(K)

        ops, op_info = [], []
        for k in range(len(self.lasers)):
            Vk = sp.csr_matrix(V[k])
            Vd = Vk.conj().T
            ops.append(commutatorSuperoperator(Vk + Vd))
            op_info.append((k, "re"))
            ops.append(commutatorSuperoperator(1j*(Vk - Vd)))
            op_info.append((k, "im"))
            ops.append(commutatorSuperoperator(D[k]))
            op_info.append((k, "detuning"))
        return Model(n, sp.csr_matrix(L0), ops, op_info, H0, V, D, collapse, K)

    def hamiltonian(self, t=None, velocity=0.0):
        """The Hermitian Hamiltonian (Grad/s) in the rotating frame at time ``t`` (ns) for an atom with
        ``velocity`` (m/s). Static parts only when ``t`` is None."""
        model = self.build()
        H = model.H0.copy()
        for k, laser in enumerate(self.lasers):
            A = 1.0 if (laser.envelope is None or t is None) else complex(np.asarray(laser.envelope(np.asarray([t])))[0])
            delta = laser.detuning + laser.dopplerShift(velocity)
            if laser.detuning_t is not None and t is not None:
                delta += float(np.asarray(laser.detuning_t(np.asarray([t])))[0])
            H += A*model.V[k] + np.conj(A)*model.V[k].conj().T + delta*model.D[k]
        return H

    # ---------------------------------------------------------------- coefficients
    def _coefficients(self, model, time, ensemble, velocity, n_ctrl):
        """Returns (ctrl_t, ctrl_u, alpha, beta) so that op j of member m has coefficient
        alpha[m, j]*u_j(t) + beta[m, j]."""
        nl = len(self.lasers)
        M = ensemble.size
        nops = len(model.ops)
        tdep = any(laser.time_dependent for laser in self.lasers)
        if tdep:
            t0, t1 = float(np.min(time)), float(np.max(time))
            ctrl_t = np.linspace(t0, t1, max(int(n_ctrl), 2))
        else:
            ctrl_t = np.zeros(0)
        ctrl_u = np.ones((nops, len(ctrl_t)))
        alpha = np.zeros((M, nops))
        beta = np.zeros((M, nops))
        for j, (k, kind) in enumerate(model.op_info):
            laser = self.lasers[k]
            scale = ensemble.rabi_scale[:, k]
            if kind in ("re", "im"):
                if laser.envelope is not None:
                    A = np.asarray(laser.envelope(ctrl_t), dtype=complex)*np.ones(len(ctrl_t))
                    ctrl_u[j] = A.real if kind == "re" else A.imag
                    alpha[:, j] = scale
                else:
                    if kind == "re":
                        beta[:, j] = scale  # A = 1
            else:
                static = laser.detuning + laser.dopplerShift(velocity)
                beta[:, j] = static + ensemble.detuning[:, k]
                if laser.detuning_t is not None:
                    ctrl_u[j] = np.asarray(laser.detuning_t(ctrl_t), dtype=float)*np.ones(len(ctrl_t))
                    alpha[:, j] = 1.0
        return ctrl_t, ctrl_u, alpha, beta

    def liouvillian(self, velocity=0.0, detuning=None, rabi_scale=None, sparse=True):
        """The full (time-independent) Liouvillian superoperator for a continuous-wave problem.

        Parameters:
            velocity (float or 3-vector): Atomic velocity in m/s.
            detuning (sequence): Extra detuning per laser (Grad/s).
            rabi_scale (sequence): Rabi-frequency factor per laser.
            sparse (bool): Return a scipy.sparse matrix (else a dense ndarray).

        Returns:
            (n^2 x n^2) matrix acting on row-major vec(rho).
        """
        nl = len(self.lasers)
        ens = Ensemble(np.ones(1),
                       np.zeros((1, nl)) if detuning is None else np.atleast_2d(detuning),
                       np.ones((1, nl)) if rabi_scale is None else np.atleast_2d(rabi_scale))
        model = self.build()
        ctrl_t, ctrl_u, alpha, beta = self._coefficients(model, np.array([0.0, 1.0]), ens, velocity, 2)
        if len(ctrl_t):
            raise ValueError("The system is time dependent; use evolve() instead.")
        L = model.L0.copy()
        for j, op in enumerate(model.ops):
            c = alpha[0, j] + beta[0, j]
            if c != 0:
                L = L + c*op
        return L.tocsr() if sparse else L.toarray()

    # ---------------------------------------------------------------- dynamics
    def groundStateRho0(self, states=None):
        """Density matrix with population shared equally between ``states`` (default: lower states of
        the first laser, or all states)."""
        if states is None:
            states = self.lasers[0].lower if self.lasers else self.states
        rho0 = np.zeros((self.n, self.n), dtype=complex)
        for s in states:
            rho0[self._idx(s), self._idx(s)] = 1/len(states)
        return rho0

    def evolve(self, time, rho0=None, ensemble=None, velocity=0.0, backend="auto",
               rtol=1e-8, atol=1e-10, max_step=None, n_ctrl=4001):
        """Time evolves the density matrix.

        Parameters:
            time (array): Increasing output times in ns. rho0 is the state at time[0].
            rho0 (ndarray): Initial density matrix (n x n). Default: equal ground-state populations.
            ensemble (Ensemble): Optional average over velocity classes / beam intensities.
            velocity (float or 3-vector): Velocity (m/s) of the atom (added to any ensemble).
            backend (str): "auto", "rust", "scipy" or "expm" (exact propagator; time-independent only).
            rtol, atol (float): Tolerances of the adaptive integrator (rust/scipy).
            max_step (float): Maximum step in ns for adaptive integration.
            n_ctrl (int): Number of samples used to tabulate time-dependent envelopes/detunings.

        Returns:
            ndarray: complex array of shape (len(time), n, n) with rho(t).
        """
        time = np.asarray(time, dtype=float)
        if time.ndim != 1 or len(time) == 0:
            raise ValueError("time must be a non-empty 1D array")
        if np.any(np.diff(time) < 0):
            raise ValueError("time must be increasing")
        n = self.n
        rho0 = self.groundStateRho0() if rho0 is None else np.asarray(rho0, dtype=complex)
        if rho0.shape != (n, n):
            raise ValueError(f"rho0 must have shape ({n}, {n})")
        if not np.allclose(rho0, rho0.conj().T, atol=1e-12):
            warnings.warn("rho0 is not Hermitian.")
        ensemble = Ensemble.single(len(self.lasers)) if ensemble is None else ensemble
        model = self.build()
        ctrl_t, ctrl_u, alpha, beta = self._coefficients(model, time, ensemble, velocity, n_ctrl)
        time_dependent = len(ctrl_t) > 0

        if backend == "auto":
            if not time_dependent and (ensemble.size == 1 or not _backend.HAVE_RUST):
                backend = "expm"
            else:
                backend = "rust" if _backend.HAVE_RUST else "scipy"
        if backend == "expm" and time_dependent:
            raise ValueError("backend='expm' only supports time-independent problems.")

        y0 = rho0.reshape(-1)
        # Restrict to the elements of rho reachable from rho0 (exact: the others stay zero)
        active = [j for j in range(len(model.ops)) if np.any(alpha[:, j]) or np.any(beta[:, j])]
        idx = reachableSubspace(model.L0, [model.ops[j] for j in active], y0)
        L0 = model.L0[idx][:, idx].tocsr()
        ops = [op[idx][:, idx].tocsr() for op in model.ops]
        ys = y0[idx]

        if backend == "expm":
            sub = np.zeros((len(time), len(idx)), dtype=complex)
            for m in range(ensemble.size):
                L = L0.copy()
                for j in active:
                    c = alpha[m, j] + beta[m, j]
                    if c != 0:
                        L = L + c*ops[j]
                sub += ensemble.weights[m]*_backend.propagate_expm(L.toarray(), ys, time)
        else:
            sub = _backend.evolve_ensemble(L0, ops, ctrl_t, ctrl_u, alpha, beta,
                                           ensemble.weights, ys, time, rtol=rtol, atol=atol,
                                           max_step=max_step, backend=backend)
        out = np.zeros((len(time), n*n), dtype=complex)
        out[:, idx] = sub
        return out.reshape(len(time), n, n)

    def steadyState(self, velocity=0.0, ensemble=None):
        """Steady-state density matrix of a continuous-wave (time-independent) system.

        Returns:
            ndarray: (n x n) density matrix (ensemble-weighted if an ensemble is given).
        """
        ensemble = Ensemble.single(len(self.lasers)) if ensemble is None else ensemble
        n = self.n
        rho = np.zeros((n, n), dtype=complex)
        for m in range(ensemble.size):
            L = self.liouvillian(velocity, ensemble.detuning[m], ensemble.rabi_scale[m], sparse=False)
            A = L.copy()
            b = np.zeros(n*n, dtype=complex)
            A[0, :] = np.eye(n).reshape(-1)  # replace one equation by tr(rho) = 1
            b[0] = 1
            x, *_ = la.lstsq(A, b)
            rho += ensemble.weights[m]*x.reshape(n, n)
        return rho

    # ---------------------------------------------------------------- outputs
    def populations(self, rho_t):
        """Returns {label: population(t)} from an (T, n, n) array."""
        return {s.label: rho_t[:, i, i].real for i, s in enumerate(self.states)}

    def describe(self):
        """A JSON-serialisable summary of the system (useful for logging and agents)."""
        model = self.build()
        lasers = []
        for k, laser in enumerate(self.lasers):
            Vk = model.V[k]
            couplings = [{"upper": int(i+1), "lower": int(j+1), "half_rabi": [float(Vk[i, j].real), float(Vk[i, j].imag)]}
                         for i, j in zip(*np.nonzero(Vk))]
            lasers.append({"name": laser.name, "wavelength_m": laser.wavelength,
                           "angular_frequency_Grad_s": laser.angular_frequency,
                           "half_rabi_Grad_s": laser.halfRabi(),
                           "polarisation": {int(q): [complex(a).real, complex(a).imag] for q, a in asPolarisation(laser.polarisation).items()},
                           "detuning_Grad_s": laser.detuning, "time_dependent": laser.time_dependent,
                           "linewidth_Grad_s": laser.linewidth, "couplings": couplings})
        return {
            "n_states": self.n,
            "states": [{"label": int(s.label), "w": float(s.w), "L": float(s.L), "S": float(s.S),
                        "J": float(s.J), "I": float(s.I), "F": float(s.F), "m": float(s.m)} for s in self.states],
            "rotating_frame_energies_Grad_s": [float(x) for x in np.diag(model.H0).real],
            "lasers": lasers,
            "decays": [{"upper": [int(s.label) for s in d.upper], "lower": [int(s.label) for s in d.lower],
                        "tau_ns": d.tau, "branching": d.branching} for d in self.decays],
            "losses": [{"states": [int(s.label) for s in l.states], "tau_ns": l.tau} for l in self.losses],
            "dephasings": [{"states": [int(s.label) for s in d.states], "rate_Grad_s": d.rate} for d in self.dephasings],
            "n_collapse_operators": len(model.collapse),
        }
