"""Tests of the array-based backend: numeric Wigner symbols, the Basis storage and decay operators."""
import numpy as np
import pytest
from sympy import Rational
from sympy.physics.wigner import wigner_3j, wigner_6j

import LASED as las
from LASED.basis import Basis
from LASED.half_rabi_freq import couplingExact
from LASED.wigner import wigner3j, wigner6j
from LASED.engine import commutatorSuperoperator, dissipatorSuperoperator
from LASED._backend import propagate_expm


def test_wigner3j_matches_sympy():
    rng = np.random.default_rng(0)
    for _ in range(2000):
        j1, j2 = rng.integers(0, 10, 2)
        j3 = rng.choice(np.arange(abs(j1 - j2), j1 + j2 + 1, 2))
        m1 = rng.choice(np.arange(-j1, j1 + 1, 2))
        m2 = rng.choice(np.arange(-j2, j2 + 1, 2))
        m3 = -(m1 + m2)
        exact = 0.0 if abs(m3) > j3 else float(wigner_3j(*(Rational(int(x), 2) for x in (j1, j2, j3, m1, m2, m3))))
        assert wigner3j(j1, j2, j3, m1, m2, m3) == pytest.approx(exact, abs=1e-13)


def test_wigner6j_matches_sympy():
    rng = np.random.default_rng(1)
    checked = 0
    while checked < 500:
        a, b = rng.integers(0, 9, 2)
        c = rng.choice(np.arange(abs(a - b), a + b + 1, 2))
        d = rng.integers(0, 9)
        e = rng.choice(np.arange(abs(d - c), d + c + 1, 2))      # (d e c) triangle
        f_choices = [f for f in range(0, 17) if abs(a - e) <= f <= a + e and (a + e + f) % 2 == 0
                     and abs(d - b) <= f <= d + b and (d + b + f) % 2 == 0]
        if not f_choices:
            continue
        f = rng.choice(f_choices)
        exact = float(wigner_6j(*(Rational(int(x), 2) for x in (a, b, c, d, e, f))))
        assert wigner6j(a, b, c, d, e, f) == pytest.approx(exact, abs=1e-13)
        checked += 1
    assert wigner6j(2, 2, 8, 2, 2, 2) == 0.0          # triangle violated


def test_wigner_vectorises():
    m = np.array([-2, 0, 2])
    vals = wigner3j(2, 2, 0, m, -m, 0)
    assert np.allclose(vals, [1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3)])


def cs_states():
    I = 7/2
    G = las.generateSubStates(1, 0, 0, 1/2, I=I, F=3) + las.generateSubStates(8, 0, 0, 1/2, I=I, F=4)
    E = (las.generateSubStates(17, 1, 1, 1/2, J=3/2, I=I, F=2) + las.generateSubStates(22, 1, 1, 1/2, J=3/2, I=I, F=3)
         + las.generateSubStates(29, 1, 1, 1/2, J=3/2, I=I, F=4) + las.generateSubStates(38, 1, 1, 1/2, J=3/2, I=I, F=5))
    return G, E


def test_coupling_matches_exact_sympy():
    G, E = cs_states()
    for e in E[::3]:
        for g in G:
            for q in (-1, 0, 1):
                assert las.coupling(e, g, q) == pytest.approx(couplingExact(e, g, q), abs=1e-14)


def test_basis_dipole_block_matches_coupling():
    G, E = cs_states()
    basis = Basis.fromStates(G + E)
    assert basis.n == 48 and basis.qn2["m"].dtype.kind == "i"
    for q in (-1, 0, 1):
        ref = np.array([[las.coupling(e, g, q) for g in G] for e in E])
        assert np.abs(basis.dipoleBlock(E, G, q) - ref).max() < 1e-15


def test_basis_rejects_bad_labels():
    a = las.State(label=1, w=0, m=0, L=0, S=0)
    b = las.State(label=3, w=0, m=0, L=1, S=0)
    with pytest.raises(ValueError, match="labels"):
        Basis.fromStates([a, b])


def test_state_table():
    G, E = cs_states()
    table = las.AtomicSystem(G + E).stateTable()
    F = table["F"] if isinstance(table, dict) else table["F"].to_numpy()
    assert len(F) == 48 and F[0] == 3 and F[-1] == 5


def test_decay_needs_one_jump_operator_per_polarisation():
    """Spontaneous emission must be described by one jump operator per photon polarisation q,
    C_q = sum_{e,g} c(e,g,q)|g><e|, not one operator per e -> g path. Only the former transfers
    the coherence of an excited superposition to the ground states (here |rho_g-,g+| -> 1/6)."""
    w0 = las.angularFreq(795e-9)
    G = las.generateSubStates(1, 0, L=0, S=1/2)               # S1/2
    E = las.generateSubStates(3, w0, L=1, S=1/2, J=1/2)       # P1/2
    tau = 27.7
    system = las.AtomicSystem(G + E, decays=[las.Decay(E, G, tau)])
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[2:, 2:] = 0.5                                        # (|e,-1/2> + |e,+1/2>)/sqrt(2)
    t = np.linspace(0, 20*tau, 5)
    per_q = system.evolve(t, rho0)
    assert abs(per_q[-1, 0, 1]) == pytest.approx(1/6, abs=1e-6)

    # The same decay rates, but with one jump operator per e -> g path
    model = system.build()
    L = commutatorSuperoperator(model.H0)
    for e in E:
        norm = sum(las.coupling(e, g, q)**2 for g in G for q in (-1, 0, 1))
        for g in G:
            rate = sum(las.coupling(e, g, q)**2 for q in (-1, 0, 1))/norm/tau
            if rate:
                C = np.zeros((4, 4))
                C[g.label - 1, e.label - 1] = np.sqrt(rate)
                L = L + dissipatorSuperoperator(C)
    per_path = propagate_expm(L.toarray(), rho0.reshape(-1), t).reshape(-1, 4, 4)
    assert np.allclose(np.diagonal(per_path[-1]), np.diagonal(per_q[-1]), atol=1e-6)   # same populations
    assert abs(per_path[-1, 0, 1]) < 1e-10                                             # coherence lost
