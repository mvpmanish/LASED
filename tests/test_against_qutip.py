"""Independent check of LASED against a QuTiP model built from Cartesian dipole operators.

The reference model uses only standard Clebsch-Gordan coefficients (Condon-Shortley phases) and the
Cartesian dot product d.eps, so it shares no polarisation-convention code with LASED. Agreement for
elliptical light and coherent initial states tests the spherical-basis convention, the Hermitian
coupling (Omega vs Omega*), and the coherence transfer in spontaneous decay.
"""
import numpy as np
import pytest
from sympy.physics.quantum.cg import CG
from sympy import S

import LASED as las

qutip = pytest.importorskip("qutip")


def reference_operators(Fg, Fe):
    """Dipole operators d_x, d_y, d_z on the space [ground m=-Fg..Fg, excited m=-Fe..Fe]."""
    ng, ne = int(2*Fg+1), int(2*Fe+1)
    n = ng + ne
    mg = np.arange(-Fg, Fg+1)
    me = np.arange(-Fe, Fe+1)
    d = {q: np.zeros((n, n), dtype=complex) for q in (-1, 0, 1)}
    for q in (-1, 0, 1):
        for i, m in enumerate(mg):
            for j, mp in enumerate(me):
                c = float(CG(S(Fg), S(m), 1, q, S(Fe), S(mp)).doit())
                d[q][ng+j, i] = c          # <e|d_q|g> (raising part only)
    # Hermitian Cartesian components from the spherical (raising) parts:
    # d_{+-1} = -+(d_x +- i d_y)/sqrt2  =>  d_x = (d_-1 - d_+1)/sqrt2, d_y = i(d_-1 + d_+1)/sqrt2
    up_x = (d[-1] - d[1])/np.sqrt(2)
    up_y = 1j*(d[-1] + d[1])/np.sqrt(2)
    up_z = d[0]
    return ng, ne, up_x, up_y, up_z


def reference_evolution(Fg, Fe, jones_vec, Omega, Gamma, delta, psi0, t):
    ng, ne, ux, uy, uz = reference_operators(Fg, Fe)
    Ex, Ey, Ez = jones_vec
    up = Ex*ux + Ey*uy + Ez*uz                     # excited <- ground block of d.eps
    H = -(Omega/2)*(up + up.conj().T)
    H[ng:, ng:] -= delta*np.eye(ne)                # rotating frame, laser detuning delta
    # Spontaneous emission: jump operators from the Cartesian lowering parts (unitary mix of d_q)
    c_ops = [np.sqrt(Gamma)*u.conj().T for u in (ux, uy, uz)]
    norm = sum((c.conj().T @ c)[ng, ng] for c in c_ops).real / Gamma
    c_ops = [qutip.Qobj(c/np.sqrt(norm)) for c in c_ops]
    rho0 = qutip.Qobj(np.outer(psi0, psi0.conj()))
    res = qutip.mesolve(qutip.Qobj(H), rho0, t, c_ops=c_ops, options={"atol": 1e-12, "rtol": 1e-10})
    return np.array([s.full() for s in res.states])


def lased_evolution(Fg, Fe, jones_vec, Omega, Gamma, delta, psi0, t, backend):
    w0 = las.angularFreq(780e-9)
    G = las.generateSubStates(label_from=1, w=0, L=0, S=1/2, I=3/2, F=Fg)
    E = las.generateSubStates(label_from=len(G)+1, w=w0, L=1, S=1/2, J=3/2, I=3/2, F=Fe)
    # Scale LASED's reduced Rabi frequency so the pi (m=0 -> m=0) element matches the CG reference
    g0 = [g for g in G if g.m == 0][0]
    e0 = [e for e in E if e.m == 0][0]
    cg = float(CG(S(Fg), 0, 1, 0, S(Fe), 0).doit())
    rabi = Omega*abs(cg/las.coupling(e0, g0, 0))
    laser = las.Laser(lower=G, upper=E, wavelength=780e-9, rabi=rabi, detuning=delta,
                      polarisation=las.jones(*jones_vec, normalise=False))
    system = las.AtomicSystem(G + E, lasers=[laser], decays=[las.Decay(E, G, tau=1/Gamma)])
    return system.evolve(t, np.outer(psi0, psi0.conj()), backend=backend)


CASES = [
    # (Fg, Fe, Jones vector)
    (1, 2, (np.cos(0.3), 1j*np.sin(0.3), 0)),                         # elliptical
    (1, 2, (np.cos(np.pi/6), np.sin(np.pi/6), 0)),                    # linear at 30 deg
    (1, 2, (0.6, 0.48*np.exp(0.7j), 0.64)),                            # general 3D, complex
    (2, 1, (np.cos(0.4), np.exp(1j*1.1)*np.sin(0.4), 0)),
]


_REFERENCE_CACHE = {}


@pytest.mark.parametrize("case", range(len(CASES)))
@pytest.mark.parametrize("backend", ["expm", "scipy", "rust"])
def test_matches_qutip_reference(case, backend):
    if backend == "rust" and not las.HAVE_RUST:
        pytest.skip("Rust backend not built")
    Fg, Fe, jones_vec = CASES[case]
    rng = np.random.default_rng(1)
    ng = 2*Fg + 1
    psi0 = np.zeros(ng + 2*Fe + 1, dtype=complex)
    psi0[:ng] = rng.normal(size=ng) + 1j*rng.normal(size=ng)    # coherent ground superposition
    psi0 /= np.linalg.norm(psi0)
    t = np.linspace(0, 60, 31)
    Omega, Gamma, delta = 2*np.pi*0.03, 1/26.2, 2*np.pi*0.004
    if case not in _REFERENCE_CACHE:
        _REFERENCE_CACHE[case] = reference_evolution(Fg, Fe, jones_vec, Omega, Gamma, delta, psi0, t)
    ref = _REFERENCE_CACHE[case]
    got = lased_evolution(Fg, Fe, jones_vec, Omega, Gamma, delta, psi0, t, backend)
    # Compare full density matrices up to the (convention-dependent) sign of each excited sub-state
    # by comparing populations and |coherences|, and the ground-state block exactly.
    ng = 2*Fg + 1
    assert np.allclose(np.diagonal(got, axis1=1, axis2=2), np.diagonal(ref, axis1=1, axis2=2), atol=2e-6)
    assert np.allclose(got[:, :ng, :ng], ref[:, :ng, :ng], atol=2e-6)
    assert np.allclose(np.abs(got), np.abs(ref), atol=2e-6)
