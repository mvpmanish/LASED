"""Physics and consistency tests for LASED."""
import numpy as np
import pytest

import LASED as las

RUST = pytest.mark.skipif(not las.HAVE_RUST, reason="Rust backend not built")


# ----------------------------------------------------------------------------- fixtures

def calcium():
    """Ca 4s2 1S0 -> 4s4p 1P1 (no hyperfine structure)."""
    wl = 422.8e-9
    G = las.generateSubStates(label_from=1, w=0, L=0, S=0)
    E = las.generateSubStates(label_from=2, w=las.angularFreq(wl), L=1, S=0)
    return G, E, wl, 4.6


def rb87_d2():
    """Rb-87 5S1/2 (F=1,2) -> 5P3/2 (F'=0..3), n = 24."""
    wl = 780.241209686e-9
    w0 = las.angularFreq(wl)
    twopi = 2*np.pi
    G = (las.generateSubStates(label_from=1, w=-twopi*6.834682611, L=0, S=1/2, I=3/2, F=1)
         + las.generateSubStates(label_from=4, w=0, L=0, S=1/2, I=3/2, F=2))
    hf = {0: -twopi*(0.2666500+0.1569470+0.0722180), 1: -twopi*(0.2666500+0.1569470),
          2: -twopi*0.2666500, 3: 0.0}
    E, label = [], 9
    for F in range(4):
        E += las.generateSubStates(label_from=label, w=w0 + hf[F], L=1, S=1/2, J=3/2, I=3/2, F=F)
        label += 2*F + 1
    return G, E, wl, 26.2348


def is_physical(rho_t, atol=1e-8):
    herm = np.abs(rho_t - rho_t.conj().transpose(0, 2, 1)).max()
    trace = np.abs(np.trace(rho_t, axis1=1, axis2=2) - 1).max()
    min_eig = min(np.linalg.eigvalsh(0.5*(r + r.conj().T)).min() for r in rho_t)
    return herm < atol and trace < 1e-6 and min_eig > -1e-7


# ----------------------------------------------------------------------------- the elliptical bug

@pytest.mark.parametrize("factors", [[1, 1j], [np.cos(0.3), 1j*np.sin(0.3)], [0.4, 0.3*np.exp(2.1j)]])
def test_complex_rabi_factors_stay_physical(factors):
    """Regression test: complex rabi_factors used to make the generator non-Hermitian (populations < 0)."""
    G, E, wl, tau = calcium()
    system = las.LaserAtomSystem(E, G, tau, [1, -1], wl, laser_intensity=100, rabi_factors=factors)
    system.setRho_0(1, 1, 1)
    rho = system.timeEvolution(np.linspace(0, 200, 201))
    assert is_physical(rho)


def test_relative_phase_is_gauge_for_j0_to_j1():
    """For J=0 -> J=1 a phase between the sigma+ and sigma- amplitudes only rephases one excited state."""
    G, E, wl, tau = calcium()
    t = np.linspace(0, 100, 51)
    pops = []
    for phase in (0, np.pi/2, 2.0):
        system = las.LaserAtomSystem(E, G, tau, [1, -1], wl, laser_intensity=50,
                                     rabi_factors=[1/np.sqrt(2), np.exp(1j*phase)/np.sqrt(2)])
        system.setRho_0(1, 1, 1)
        rho = system.timeEvolution(t)
        pops.append(np.real(np.diagonal(rho, axis1=1, axis2=2)))
    assert np.allclose(pops[0], pops[1], atol=1e-10)
    assert np.allclose(pops[0], pops[2], atol=1e-10)


def test_liouvillian_structure_for_random_elliptical_light():
    G, E, wl, tau = rb87_d2()
    rng = np.random.default_rng(0)
    pol = las.jones(*(rng.normal(size=3) + 1j*rng.normal(size=3)))
    laser = las.Laser(G, E, wl, intensity=10, tau=tau, polarisation=pol, detuning=0.3)
    system = las.AtomicSystem(G + E, lasers=[laser], decays=[las.Decay(E, G, tau)])
    L = system.liouvillian(sparse=False)
    n = system.n
    # Trace preservation: d tr(rho)/dt = vec(I)^T L vec(rho) = 0 for all rho
    assert np.abs(np.eye(n).reshape(-1) @ L).max() < 1e-10
    # Hermiticity preservation
    X = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    X = X + X.conj().T
    Y = (L @ X.reshape(-1)).reshape(n, n)
    assert np.abs(Y - Y.conj().T).max() < 1e-10
    # Stability: no eigenvalue with positive real part (on the F=2 -> F'=3 cycling sub-system, for speed)
    G2, E3 = [g for g in G if g.F == 2], [e for e in E if e.F == 3]
    for i, s in enumerate(G2 + E3):
        s.label = i + 1
    small = las.AtomicSystem(G2 + E3, [las.Laser(G2, E3, wl, intensity=10, tau=tau, polarisation=pol)],
                             [las.Decay(E3, G2, tau)])
    assert np.linalg.eigvals(small.liouvillian(sparse=False)).real.max() < 1e-8
    # The Hamiltonian is Hermitian
    H = system.hamiltonian()
    assert np.allclose(H, H.conj().T)


def test_rb87_elliptical_dynamics_physical():
    G, E, wl, tau = rb87_d2()
    system = las.LaserAtomSystem(E, G, tau, laser_wavelength=wl, laser_intensity=5,
                                 polarisation=las.elliptical(angle=0.4, ellipticity=0.2))
    rho = system.timeEvolution(np.linspace(0, 200, 41), detuning=2*np.pi*0.1)
    assert is_physical(rho)


# ----------------------------------------------------------------------------- polarisation

def test_polarisation_helpers():
    s2 = 1/np.sqrt(2)
    x = las.jones(1, 0)
    assert np.allclose([x[-1], x[0], x[1]], [s2, 0, s2])
    sp_ = las.elliptical(0, np.pi/4)
    assert abs(sp_[1]) == pytest.approx(1) and abs(sp_[-1]) < 1e-12
    sm = las.elliptical(0.7, -np.pi/4)
    assert abs(sm[-1]) == pytest.approx(1) and abs(sm[1]) < 1e-12
    assert las.linear(polar=0)[0] == pytest.approx(1)
    assert las.fromRabiFactors([1, -1], [1, -1], 1/np.sqrt(2)) == {1: s2, -1: -s2}


def test_rotating_linear_polarisation_is_a_rotation_about_z():
    """Rotating the polarisation by psi about z multiplies coherences rho_mm' by exp(-i(m - m')psi)."""
    G, E, wl, tau = calcium()
    t = np.linspace(0, 30, 7)
    psi = 0.6
    results = []
    for angle in (0.0, psi):
        laser = las.Laser(G, E, wl, intensity=50, tau=tau, polarisation=las.linear(angle))
        sys_ = las.AtomicSystem(G + E, lasers=[laser], decays=[las.Decay(E, G, tau)])
        results.append(sys_.evolve(t))
    m = np.array([0, -1, 0, 1])
    phase = np.exp(-1j*(m[:, None] - m[None, :])*psi)
    assert np.allclose(results[1], results[0]*phase, atol=1e-10)


# ----------------------------------------------------------------------------- analytic limits

def two_level(rabi_eff, tau=None, detuning=0.0, envelope=None, detuning_t=None):
    """J=0 -> J=1 with pi light: an exact two-level system (m=0 <-> m=0)."""
    G, E, wl, tau_ca = calcium()
    c0 = abs(las.coupling(E[1], G[0], 0))
    laser = las.Laser(G, E, wl, rabi=rabi_eff/c0, polarisation=las.pi_polarised(), detuning=detuning,
                      envelope=envelope, detuning_t=detuning_t)
    decays = [] if tau is None else [las.Decay(E, G, tau)]
    return las.AtomicSystem(G + E, lasers=[laser], decays=decays), G, E


def test_two_level_rabi_oscillation():
    rabi = 2*np.pi*0.05
    system, G, E = two_level(rabi)
    t = np.linspace(0, 60, 121)
    rho = system.evolve(t)
    assert np.allclose(rho[:, 2, 2].real, np.sin(rabi*t/2)**2, atol=1e-10)


def test_two_level_steady_state():
    rabi, tau, delta = 2*np.pi*0.02, 4.6, 2*np.pi*0.015
    system, G, E = two_level(rabi, tau=tau, detuning=delta)
    gamma = 1/tau
    expected = (rabi**2/4)/(delta**2 + rabi**2/2 + gamma**2/4)
    assert system.steadyState()[2, 2].real == pytest.approx(expected, rel=1e-9)
    rho = system.evolve(np.linspace(0, 200, 3))
    assert rho[-1, 2, 2].real == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("backend", ["scipy", pytest.param("rust", marks=RUST)])
def test_gaussian_pi_pulse_inverts(backend):
    rabi = 2*np.pi*0.05
    sigma = np.pi/(rabi*np.sqrt(2*np.pi))     # pulse area pi
    pulse = las.pulses.gaussian(t0=8*sigma, sigma=sigma)
    system, G, E = two_level(rabi, envelope=pulse)
    rho = system.evolve(np.linspace(0, 16*sigma, 5), backend=backend, n_ctrl=20001)
    assert rho[-1, 2, 2].real == pytest.approx(1.0, abs=1e-5)
    assert pulse.area(rabi, 0, 16*sigma) == pytest.approx(np.pi, rel=1e-6)


@pytest.mark.parametrize("backend", ["scipy", pytest.param("rust", marks=RUST)])
def test_pulse_phase_controls_rotation_axis(backend):
    """pi/2 then pi/2 with phase pi returns the atom to the ground state (Ramsey echo)."""
    rabi = 2*np.pi*0.05
    T = (np.pi/2)/rabi
    pulse = (las.pulses.square(0, T) + las.pulses.square(T, 2*T).withPhase(np.pi))
    system, G, E = two_level(rabi, envelope=pulse)
    rho = system.evolve(np.linspace(0, 2*T, 3), backend=backend, n_ctrl=40001, max_step=T/50)
    assert rho[1, 2, 2].real == pytest.approx(0.5, abs=2e-3)
    assert rho[-1, 2, 2].real == pytest.approx(0.0, abs=2e-3)


@pytest.mark.parametrize("backend", ["scipy", pytest.param("rust", marks=RUST)])
def test_adiabatic_rapid_passage(backend):
    rabi, T = 2*np.pi*0.05, 3000.0
    chirp = las.pulses.linearChirp(rate=6/T, t0=T/2)    # sweep -3 -> +3 Grad/s
    system, G, E = two_level(rabi, envelope=las.pulses.sinSquared(0, T), detuning_t=chirp)
    rho = system.evolve(np.linspace(0, T, 3), backend=backend)
    assert rho[-1, 2, 2].real > 0.99


def test_doppler_shift_equals_detuning():
    G, E, wl, tau = calcium()
    t = np.linspace(0, 50, 11)
    v = 30.0
    laser = las.Laser(G, E, wl, intensity=20, tau=tau)
    sys_v = las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)])
    shift = laser.dopplerShift(v)
    assert shift == pytest.approx(-2*np.pi*v/wl*1e-9)
    laser_d = las.Laser(G, E, wl, intensity=20, tau=tau, detuning=shift)
    sys_d = las.AtomicSystem(G + E, [laser_d], [las.Decay(E, G, tau)])
    assert np.allclose(sys_v.evolve(t, velocity=v), sys_d.evolve(t), atol=1e-12)


# ----------------------------------------------------------------------------- multi-laser

def lambda_system(delta1, delta2):
    """J=1 (m=-1, +1) -> J'=0 Lambda system driven by two lasers (sigma+ and sigma-)."""
    wl = 795e-9
    gm = las.State(label=1, w=0, m=-1, L=0, S=1, J=1)
    gp = las.State(label=2, w=0, m=1, L=0, S=1, J=1)
    e = las.State(label=3, w=las.angularFreq(wl), m=0, L=1, S=1, J=0)
    rabi = 2*np.pi*0.01
    l1 = las.Laser([gm], [e], wl, rabi=rabi, polarisation=las.circular(+1), detuning=delta1)
    l2 = las.Laser([gp], [e], wl, rabi=rabi, polarisation=las.circular(-1), detuning=delta2)
    return las.AtomicSystem([gm, gp, e], lasers=[l1, l2], decays=[las.Decay([e], [gm, gp], tau=26.0)])


def test_coherent_population_trapping():
    dark = lambda_system(2*np.pi*0.01, 2*np.pi*0.01).steadyState()
    bright = lambda_system(2*np.pi*0.01, 2*np.pi*0.013).steadyState()
    assert dark[2, 2].real < 1e-8
    assert bright[2, 2].real > 1e-3
    # The dark state is the antisymmetric superposition (up to the coupling signs)
    assert abs(dark[0, 1]) == pytest.approx(0.5, abs=1e-6)


def test_closed_laser_loop_rejected():
    wl = 795e-9
    a = las.State(label=1, w=0, m=0, L=0, S=0)
    b = las.State(label=2, w=las.angularFreq(wl), m=0, L=1, S=0)
    l1 = las.Laser([a], [b], wl, rabi=0.1)
    l2 = las.Laser([a], [b], wl, rabi=0.1, detuning=0.2)
    with pytest.raises(ValueError, match="closed loop"):
        las.AtomicSystem([a, b], lasers=[l1, l2]).build()


def test_laser_linewidth_dephases_optical_coherence():
    rabi = 1e-6     # negligible drive: the coherence just decays
    system, G, E = two_level(rabi)
    system.lasers[0].linewidth = 0.2
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[0, 0] = rho0[2, 2] = rho0[0, 2] = rho0[2, 0] = 0.5
    t = np.linspace(0, 10, 6)
    rho = system.evolve(t, rho0)
    assert np.allclose(np.abs(rho[:, 0, 2]), 0.5*np.exp(-0.1*t), atol=1e-6)


# ----------------------------------------------------------------------------- backends & ensembles

@RUST
def test_backends_agree_for_doppler_ensemble():
    G, E, wl, tau = calcium()
    laser = las.Laser(G, E, wl, intensity=50, tau=tau, polarisation=las.elliptical(0.2, 0.3))
    system = las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)])
    ens = las.dopplerEnsemble([laser], sigma_v=40.0, n=15)*las.gaussianBeamEnsemble(1, n=4)
    t = np.linspace(0, 40, 21)
    a = system.evolve(t, ensemble=ens, backend="expm")
    b = system.evolve(t, ensemble=ens, backend="rust", rtol=1e-10, atol=1e-12)
    c = system.evolve(t, ensemble=ens, backend="scipy", rtol=1e-10, atol=1e-12)
    assert np.abs(a - b).max() < 1e-7
    assert np.abs(a - c).max() < 1e-7
    assert is_physical(a)


@RUST
def test_rust_and_scipy_agree_time_dependent_multilevel():
    G, E, wl, tau = rb87_d2()
    pulse = las.pulses.flatTop(2, 12, rise=4).withPhase(lambda t: 0.05*t)
    laser = las.Laser(G, E, wl, intensity=30, tau=tau, polarisation=las.linear(0.3),
                      envelope=pulse, detuning_t=las.pulses.tanhChirp(2*np.pi*0.05, 8, 3))
    system = las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)])
    t = np.linspace(0, 15, 4)
    a = system.evolve(t, backend="rust", rtol=1e-8, atol=1e-10)
    b = system.evolve(t, backend="scipy", rtol=1e-8, atol=1e-10)
    assert np.abs(a - b).max() < 1e-6
    assert is_physical(a)


def test_moving_atom_through_beam():
    G, E, wl, tau = calcium()
    laser = las.Laser(G, E, wl, rabi=2*np.pi*0.2, tau=tau, direction=(0, 0, 1))
    w0, v = 50e-6, 100.0       # 50 um waist crossed at 100 m/s transverse (500 ns transit)
    env, det = las.pulses.movingAtom(laser, lambda t: np.stack([v*(t - 500)*1e-9, 0*t, 0*t], axis=1), waist=w0)
    t = np.array([0.0, 500.0, 1000.0])
    assert np.allclose(env(t).real, np.exp(-(v*(t - 500)*1e-9/w0)**2))
    assert np.allclose(det(t), 0, atol=1e-9)    # transverse motion: no first-order Doppler shift
    env2, det2 = las.pulses.movingAtom(laser, lambda t: np.stack([0*t, 0*t, 20.0*t*1e-9], axis=1), waist=w0)
    assert np.allclose(det2(t), laser.dopplerShift(20.0))


# ----------------------------------------------------------------------------- legacy API

def test_legacy_api_features_run():
    G, E, wl, tau = calcium()
    t = np.linspace(0, 20, 11)
    s = las.LaserAtomSystem(E, G, tau, [0], wl, laser_intensity=100, laser_power=100, tau_f=50, tau_b=80)
    s.timeEvolution(t, detuning=0.5, atomic_velocity=10)
    assert np.asarray(s.rho_t).shape == (16, 11)
    assert len(s.Rho_t(E[1], E[1])) == 11
    assert s.rho_t_matrix.shape == (11, 4, 4)
    assert len(s.rho_et) == 11 and s.rho_et[0].shape == (3, 3)
    s.clearRho_0()
    dets = np.linspace(-3, 3, 7)
    s.timeEvolution(t, doppler_averaging=True, doppler_width=1.0, doppler_detunings=dets,
                    beam_profile_averaging=True, r_sigma=0.75, n_beam_averaging=3)
    rho = s.rho_t_matrix
    assert np.abs(rho - rho.conj().transpose(0, 2, 1)).max() < 1e-12


def test_legacy_matrix_matches_engine():
    G, E, wl, tau = calcium()
    A = las.timeEvolutionMatrix(4, E, G, [1, -1], [1, 0, -1], tau, wl, 10, rabi_factors=[1, 1j])
    assert A.shape == (16, 16)
    rho = np.diag([1, 0, 0, 0]).astype(complex).reshape(-1)
    d = (A @ rho).reshape(4, 4)
    assert np.allclose(d, d.conj().T)


def test_zeeman_g_factors():
    assert las.landeGF(L=0, S=1/2, J=1/2, I=3/2, F=2) == pytest.approx(0.5006, abs=1e-4)
    assert las.landeGF(L=0, S=1/2, J=1/2, I=3/2, F=1) == pytest.approx(-0.5006, abs=1e-4)
    s = las.State(label=1, w=0, m=2, L=0, S=1/2, I=3/2, F=2)
    las.applyZeemanShift([s], B=1.0)
    assert s.w == pytest.approx(2*np.pi*1.39962449361e-3*0.5*2*1.00116, rel=1e-4)


def test_describe_is_json_serialisable():
    import json
    G, E, wl, tau = calcium()
    laser = las.Laser(G, E, wl, intensity=10, tau=tau, polarisation=las.elliptical(0, 0.2))
    info = las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)]).describe()
    json.dumps(info)
    assert info["n_states"] == 4 and len(info["lasers"][0]["couplings"]) == 2


def test_symbolic_equations_match_numeric_matrix():
    """The printed (sympy) equations are exactly the equations that are solved."""
    import sympy
    from LASED.symbolic_print import symbolicEquations
    G, E, wl, tau = calcium()
    Q, factors, detuning, I0 = [1, -1], [np.cos(0.3), 1j*np.sin(0.3)], 0.7, 10
    eqs = symbolicEquations(4, E, G, Q, [1, 0, -1], 50.0, None, detuning, wl, 0, 1, factors)
    A = las.timeEvolutionMatrix(4, E, G, Q, [1, 0, -1], tau, wl, I0, tau_f=50.0, detuning=detuning,
                                rabi_factors=factors)
    half = las.halfRabiFreq(I0, tau, wl)
    subs = {sympy.Symbol("tau", positive=True): tau, sympy.Symbol("tau_f", positive=True): 50.0,
            sympy.Symbol("omega_L"): las.angularFreq(wl), sympy.Symbol("delta"): detuning}
    subs.update({sympy.Symbol(f"w_{s.label}"): s.w for s in G + E})
    Omega = sympy.Function("Omega", real=True)
    subs.update({Omega(e.label, g.label, q): half*las.coupling(e, g, q) for e in E for g in G for q in Q})
    rho = [sympy.Symbol(f"rho_{i}{j}") for i in range(1, 5) for j in range(1, 5)]
    M, _ = sympy.linear_eq_to_matrix([eq.rhs.subs(subs) for eq in eqs], rho)
    M = np.array(M.evalf(), dtype=complex)
    # The two rotating frames differ by a constant energy shift, which drops out of the commutator
    assert np.abs(M - A).max() < 1e-6*np.abs(A).max()


@pytest.mark.parametrize("angles,pol", [((0, np.pi/2, 0), las.linear(0)),                       # z -> x
                                        ((np.pi/2, np.pi/2, -np.pi/2), las.linear(np.pi/2))])  # z -> y
def test_rotation_matches_polarisation_convention(angles, pol):
    """Rotate rho0, evolve with pi light, rotate back == evolve with the rotated (linear) polarisation."""
    wl = 900e-9
    rng = np.random.default_rng(3)
    psi = rng.normal(size=5) + 1j*rng.normal(size=5)
    psi /= np.linalg.norm(psi)
    t = np.linspace(0, 2000, 5)

    def make(**kw):
        G = las.generateSubStates(1, 0, 2, 0)
        E = las.generateSubStates(6, las.angularFreq(wl), 1, 0)
        s = las.LaserAtomSystem(E, G, 60e3, laser_wavelength=wl, laser_intensity=100, **kw)
        s.appendDensityMatrixToRho_0(np.outer(psi, psi.conj()), "g")
        return s
    direct = make(polarisation=pol)
    direct.timeEvolution(t)
    a, b, c = angles
    rotated = make(Q=[0])
    rotated.rotateRho_0(a, b, c)
    rotated.timeEvolution(t)
    rotated.rotateRho_t(-c, -b, -a)
    assert np.abs(direct.rho_t_matrix[:, :5, :5] - rotated.rho_t_matrix[:, :5, :5]).max() < 1e-12


def test_half_integer_sub_states():
    assert [s.m for s in las.generateSubStates(1, 0, 1, 1/2, J=3/2)] == [-1.5, -0.5, 0.5, 1.5]
    assert [s.m for s in las.generateSubStates(1, 0, 1, 0)] == [-1, 0, 1]
