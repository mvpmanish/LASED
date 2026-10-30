'''
Polarisation of the laser field in the spherical basis used by LASED.

Convention
----------
The quantisation axis is z. LASED expands the (complex, unit) polarisation
vector of the laser as

    eps = a_{-1} u_{-1} + a_0 u_0 + a_{+1} u_{+1},
    u_{+1} = (x + i y)/sqrt(2),   u_0 = z,   u_{-1} = (x - i y)/sqrt(2),

so a_q = u_q^* . eps. The amplitude a_q drives Delta m = m_upper - m_lower = q
transitions (q = +1 is sigma+, q = 0 is pi, q = -1 is sigma-). This is the
convention of ``coupling(e, g, q)`` and of the legacy ``Q``/``rabi_factors``
arguments (``rabi_factors[i]`` is a_{Q[i]}).

A polarisation is represented as a plain ``dict`` ``{q: a_q}``. Amplitudes may
be complex: a relative phase between the sigma+ and sigma- amplitudes gives
linear polarisation at an angle, and unequal magnitudes give elliptical light.

Example:
    >>> import LASED as las
    >>> las.linear(angle=0)           # x-polarised light travelling along z
    {-1: (0.7071067811865475+0j), 0: 0j, 1: (0.7071067811865475+0j)}
    >>> las.elliptical(angle=0, ellipticity=0.3)  # tilt 0, ellipticity angle 0.3 rad
'''

import numpy as np

SIGMA_PLUS = 1
PI_POL = 0
SIGMA_MINUS = -1


def _normalise(pol, normalise):
    pol = {q: complex(a) for q, a in pol.items()}
    if normalise:
        norm = np.sqrt(sum(abs(a)**2 for a in pol.values()))
        if norm == 0:
            raise ValueError("Polarisation vector has zero length.")
        pol = {q: a/norm for q, a in pol.items()}
    return pol


def jones(Ex, Ey, Ez=0, normalise=True):
    """Converts a Cartesian (Jones) polarisation vector to spherical amplitudes.

    Parameters:
        Ex (complex): x component of the polarisation vector.
        Ey (complex): y component of the polarisation vector.
        Ez (complex): z component (along the quantisation axis).
        normalise (bool): Rescale so that ``sum |a_q|^2 = 1``.

    Returns:
        dict: ``{-1: a_-1, 0: a_0, 1: a_+1}`` with a_{+-1} = (Ex -+ i Ey)/sqrt(2), a_0 = Ez.
    """
    s2 = np.sqrt(2)
    pol = {-1: (Ex + 1j*Ey)/s2, 0: Ez, 1: (Ex - 1j*Ey)/s2}
    return _normalise(pol, normalise)


def linear(angle=0.0, polar=np.pi/2):
    """Linearly polarised light.

    Parameters:
        angle (float): Azimuthal angle of the polarisation vector from the x-axis, in radians.
        polar (float): Polar angle of the polarisation vector from the quantisation (z) axis, in
            radians. The default pi/2 puts the polarisation in the xy-plane (light travelling along z);
            polar = 0 gives pi light.

    Returns:
        dict: Spherical amplitudes ``{q: a_q}``.
    """
    return jones(np.sin(polar)*np.cos(angle), np.sin(polar)*np.sin(angle), np.cos(polar))


def pi_polarised():
    """Linear polarisation along the quantisation axis (drives Delta m = 0).

    Returns:
        dict: ``{0: 1}``
    """
    return {0: 1+0j}


def circular(handedness=+1):
    """Circularly polarised light travelling along the quantisation axis.

    Parameters:
        handedness (int): +1 for sigma+ (Delta m = +1), -1 for sigma- (Delta m = -1).

    Returns:
        dict: Spherical amplitudes ``{q: a_q}``.
    """
    if handedness not in (1, -1):
        raise ValueError("handedness must be +1 (sigma+) or -1 (sigma-)")
    return {handedness: 1+0j}


def elliptical(angle=0.0, ellipticity=0.0):
    """Elliptically polarised light travelling along the quantisation axis (z).

    The polarisation ellipse has its major axis at ``angle`` from the x-axis and
    ellipticity angle ``ellipticity`` (tan(ellipticity) = minor/major axis ratio).
    ellipticity = 0 is linear, +pi/4 is sigma+ and -pi/4 is sigma-.

    Parameters:
        angle (float): Orientation of the major axis from the x-axis in radians.
        ellipticity (float): Ellipticity angle chi in [-pi/4, pi/4] radians.

    Returns:
        dict: Spherical amplitudes ``{q: a_q}``, normalised.
    """
    # Jones vector of an ellipse rotated by `angle`
    Ex = np.cos(angle)*np.cos(ellipticity) - 1j*np.sin(angle)*np.sin(ellipticity)
    Ey = np.sin(angle)*np.cos(ellipticity) + 1j*np.cos(angle)*np.sin(ellipticity)
    return jones(Ex, Ey, 0)


def fromRabiFactors(Q, rabi_factors=None, rabi_scaling=None):
    """Builds a polarisation dict from the legacy ``Q``/``rabi_factors``/``rabi_scaling`` arguments.

    Parameters:
        Q (list): Polarisations q present in the laser.
        rabi_factors (list): Complex amplitude multiplying each q in Q (default 1 each).
        rabi_scaling (float): Overall factor multiplying every amplitude (default 1).

    Returns:
        dict: ``{q: a_q}`` (not normalised, to reproduce the legacy behaviour exactly).
    """
    if rabi_factors is None:
        rabi_factors = [1]*len(Q)
    if len(rabi_factors) != len(Q):
        raise ValueError("rabi_factors must be the same length as Q: each element of Q is multiplied by the corresponding rabi_factor.")
    scale = 1 if rabi_scaling is None else rabi_scaling
    pol = {}
    for q, f in zip(Q, rabi_factors):
        if q not in (-1, 0, 1):
            raise ValueError(f"Polarisation q must be -1, 0 or +1, got {q}")
        pol[q] = pol.get(q, 0) + complex(f)*scale
    return pol


def asPolarisation(pol):
    """Coerces a polarisation given as a dict, an int q, or a length-3 sequence (a_-1, a_0, a_+1) into a dict."""
    if isinstance(pol, dict):
        return {int(q): complex(a) for q, a in pol.items()}
    if isinstance(pol, (int, np.integer)):
        return {int(pol): 1+0j}
    arr = np.asarray(pol, dtype=complex).ravel()
    if arr.shape != (3,):
        raise ValueError("Polarisation must be a dict {q: a_q}, an int q, or a sequence (a_-1, a_0, a_+1).")
    return {-1: arr[0], 0: arr[1], 1: arr[2]}
