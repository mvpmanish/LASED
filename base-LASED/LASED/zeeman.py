'''
Linear (weak-field) Zeeman shifts of atomic sub-states.

The shift of the sub-state ``|F, m_F>`` in a magnetic field B along the quantisation axis is
``Delta omega = g_F * m_F * mu_B * B / hbar``, valid while the shift is small compared with the
hyperfine splitting (use a full Breit-Rabi or
Paschen-Back treatment beyond that).

Example:
    >>> G = las.generateSubStates(label_from=1, w=0, L=0, S=1/2, I=3/2, F=2)
    >>> las.applyZeemanShift(G, B=1.0)     # 1 gauss, modifies State.w in place
'''

import numpy as np

MU_B_OVER_H_GHZ_PER_GAUSS = 1.39962449361e-3  # mu_B/h in GHz/G (CODATA 2018)
G_S = 2.00231930436256  # electron spin g-factor


def landeGJ(L, S, J, g_s=G_S, g_l=1.0):
    """Lande g-factor g_J (LS coupling)."""
    if J == 0:
        return 0.0
    return (g_l*(J*(J+1) - S*(S+1) + L*(L+1)) + g_s*(J*(J+1) + S*(S+1) - L*(L+1)))/(2*J*(J+1))


def landeGF(L, S, J, I, F, g_I=0.0):
    """Lande g-factor g_F, neglecting the nuclear g-factor by default."""
    if F == 0:
        return 0.0
    gJ = landeGJ(L, S, J)
    return (gJ*(F*(F+1) - I*(I+1) + J*(J+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1)))/(2*F*(F+1))


def zeemanShift(state, B, g_F=None):
    """Linear Zeeman shift of a sub-state in Grad/s for a field B in gauss."""
    if g_F is None:
        g_F = landeGF(state.L, state.S, state.J, state.I, state.F)
    return 2*np.pi*MU_B_OVER_H_GHZ_PER_GAUSS*g_F*state.m*B


def applyZeemanShift(states, B, g_F=None):
    """Adds the linear Zeeman shift for field B (gauss) to ``State.w`` of each state, in place.

    Returns:
        list of State: the same states (for chaining).
    """
    for s in states:
        s.w = s.w + zeemanShift(s, B, g_F)
    return states
