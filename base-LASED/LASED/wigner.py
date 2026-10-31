'''
Fast, vectorised Wigner 3j and 6j symbols (floating point, NumPy broadcasting).

All angular momenta are passed as **twice their value** (integers), so half-integers are exact:
``wigner3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3)``. Every argument may be an array; the arrays
broadcast against each other and the result has the broadcast shape. Symbols that violate the
triangle, projection or parity rules are returned as 0.

The Racah formulas are evaluated with a table of factorials, which is exact to double precision
for the angular momenta found in atoms (j up to ~40). The test suite checks them against sympy.

Example:
    >>> from LASED.wigner import wigner3j, wigner6j
    >>> wigner3j(2, 2, 0, 2, -2, 0)        # (1 1 0; 1 -1 0) = 1/sqrt(3)
    0.5773502691896258
'''

import math
import numpy as np

_FACT = np.array([math.factorial(k) for k in range(171)], dtype=float)


def _fact(x):
    """Factorial of integer array x (x >= 0 assumed where used)."""
    return _FACT[np.clip(x, 0, 170)]


def _triangle_ok(a2, b2, c2):
    """Triangle and integer-sum condition for doubled angular momenta."""
    return (c2 >= np.abs(a2 - b2)) & (c2 <= a2 + b2) & ((a2 + b2 + c2) % 2 == 0)


def _delta(a2, b2, c2):
    """Triangle coefficient Delta(abc) for doubled arguments (valid triangles only)."""
    return np.sqrt(_fact((a2 + b2 - c2)//2)*_fact((a2 - b2 + c2)//2)*_fact((-a2 + b2 + c2)//2)
                   / _fact((a2 + b2 + c2)//2 + 1))


def wigner3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (j1 j2 j3; m1 m2 m3) with all arguments doubled. Vectorised."""
    j1, j2, j3, m1, m2, m3 = np.broadcast_arrays(*(np.asarray(x, dtype=np.int64) for x in (j1, j2, j3, m1, m2, m3)))
    ok = (m1 + m2 + m3 == 0) & _triangle_ok(j1, j2, j3)
    for j, m in ((j1, m1), (j2, m2), (j3, m3)):
        ok &= (np.abs(m) <= j) & ((j + m) % 2 == 0)
    # Racah formula, in terms of integer combinations
    a = (j1 + j2 - j3)//2
    b = (j1 - m1)//2
    c = (j2 + m2)//2
    d = (j3 - j2 + m1)//2
    e = (j3 - j1 - m2)//2
    kmin = np.maximum(0, np.maximum(-d, -e))
    kmax = np.minimum(a, np.minimum(b, c))
    ok &= kmax >= kmin
    total = np.zeros(j1.shape)
    if np.any(ok):
        k_lo = int(kmin[ok].min())
        k_hi = int(kmax[ok].max())
        for k in range(k_lo, k_hi + 1):
            valid = ok & (k >= kmin) & (k <= kmax)
            term = ((-1.0)**k)/(_fact(k)*_fact(a - k)*_fact(b - k)*_fact(c - k)*_fact(d + k)*_fact(e + k))
            total += np.where(valid, term, 0.0)
    pref = _delta(j1, j2, j3)*np.sqrt(_fact((j1 + m1)//2)*_fact((j1 - m1)//2)*_fact((j2 + m2)//2)
                                      * _fact((j2 - m2)//2)*_fact((j3 + m3)//2)*_fact((j3 - m3)//2))
    phase = np.where(((j1 - j2 - m3)//2) % 2 == 0, 1.0, -1.0)
    out = np.where(ok, phase*pref*total, 0.0)
    return out[()] if out.ndim == 0 else out


def wigner6j(j1, j2, j3, j4, j5, j6):
    """Wigner 6j symbol {j1 j2 j3; j4 j5 j6} with all arguments doubled. Vectorised."""
    j1, j2, j3, j4, j5, j6 = np.broadcast_arrays(*(np.asarray(x, dtype=np.int64) for x in (j1, j2, j3, j4, j5, j6)))
    ok = (_triangle_ok(j1, j2, j3) & _triangle_ok(j1, j5, j6) & _triangle_ok(j4, j2, j6)
          & _triangle_ok(j4, j5, j3))
    t1 = (j1 + j2 + j3)//2
    t2 = (j1 + j5 + j6)//2
    t3 = (j4 + j2 + j6)//2
    t4 = (j4 + j5 + j3)//2
    p1 = (j1 + j2 + j4 + j5)//2
    p2 = (j2 + j3 + j5 + j6)//2
    p3 = (j3 + j1 + j6 + j4)//2
    kmin = np.maximum(np.maximum(t1, t2), np.maximum(t3, t4))
    kmax = np.minimum(p1, np.minimum(p2, p3))
    ok &= kmax >= kmin
    total = np.zeros(j1.shape)
    if np.any(ok):
        for k in range(int(kmin[ok].min()), int(kmax[ok].max()) + 1):
            valid = ok & (k >= kmin) & (k <= kmax)
            term = ((-1.0)**k)*_fact(k + 1)/(_fact(k - t1)*_fact(k - t2)*_fact(k - t3)*_fact(k - t4)
                                            * _fact(p1 - k)*_fact(p2 - k)*_fact(p3 - k))
            total += np.where(valid, term, 0.0)
    pref = _delta(j1, j2, j3)*_delta(j1, j5, j6)*_delta(j4, j2, j6)*_delta(j4, j5, j3)
    out = np.where(ok, pref*total, 0.0)
    return out[()] if out.ndim == 0 else out


def dipoleCoupling(upper, lower, q):
    """Dipole coupling coefficients ``<u| d.u_q |l>`` (in units of the reduced matrix element) between
    every upper and lower sub-state, as an array of shape (len(upper), len(lower)).

    ``upper`` and ``lower`` are dicts of doubled quantum numbers (integer arrays) with keys
    "L", "S", "J", "I", "F", "m", as produced by ``LASED.basis.Basis``. This is the vectorised
    equivalent of ``coupling(u, l, q)`` for all pairs at once.
    """
    u = {k: np.asarray(v, dtype=np.int64)[:, None] for k, v in upper.items()}
    l = {k: np.asarray(v, dtype=np.int64)[None, :] for k, v in lower.items()}
    # phase (-1)^(q(1+q)/2 + Fu + Fl + Ju + Jl + Iu + Lu + Su - mu + 1), exponent in doubled units
    exponent2 = (q*(1 + q)) + u["F"] + l["F"] + u["J"] + l["J"] + u["I"] + u["L"] + u["S"] - u["m"] + 2
    sign = np.where((exponent2//2) % 2 == 0, 1.0, -1.0)
    factor = np.sqrt((u["F"] + 1.0)*(l["F"] + 1.0)*(u["J"] + 1.0)*(l["J"] + 1.0)*(u["L"] + 1.0))
    sixj_1 = wigner6j(u["J"], u["F"], u["I"], l["F"], l["J"], 2)
    sixj_2 = wigner6j(u["L"], u["J"], u["S"], l["J"], l["L"], 2)
    threej = wigner3j(u["F"], 2, l["F"], -u["m"], 2*q, l["m"])
    return sign*factor*sixj_1*sixj_2*threej
