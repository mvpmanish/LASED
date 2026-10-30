'''
Numerical backends for LASED.

* ``propagate_expm``: exact propagation of a time-independent Liouvillian with matrix exponentials
  (LAPACK via SciPy). Used for continuous-wave problems.
* ``evolve_ensemble``: adaptive integration of d vec(rho)/dt = [L0 + sum_j c_{m,j}(t) M_j] vec(rho)
  for an ensemble of members m, with c_{m,j}(t) = alpha[m, j]*u_j(t) + beta[m, j].
  Runs in the compiled Rust extension ``LASED._lased_rs`` (parallel over members with rayon) when
  available, otherwise falls back to ``scipy.integrate.solve_ivp``.

``HAVE_RUST`` tells you whether the Rust extension was found.
'''

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la

try:
    from LASED import _lased_rs
    HAVE_RUST = True
except ImportError:  # pragma: no cover - depends on the build
    _lased_rs = None
    HAVE_RUST = False


def backendInfo():
    """Returns a dict describing the available numerical backends."""
    info = {"rust": HAVE_RUST}
    if HAVE_RUST:
        info["rust_version"] = _lased_rs.__version__
        info["rust_threads"] = _lased_rs.num_threads()
    return info


def propagate_expm(L, y0, time):
    """Propagates vec(rho) with the exact propagator exp(L dt).

    Parameters:
        L (ndarray): Dense (N x N) Liouvillian.
        y0 (ndarray): vec(rho) at time[0].
        time (ndarray): Increasing times.

    Returns:
        ndarray: (len(time), N) array of vec(rho(t)).
    """
    out = np.empty((len(time), len(y0)), dtype=complex)
    out[0] = y0
    cache = {}
    y = np.asarray(y0, dtype=complex)
    for i in range(1, len(time)):
        dt = time[i] - time[i-1]
        key = round(dt, 12)
        P = cache.get(key)
        if P is None:
            P = la.expm(L*dt)
            if len(cache) < 64:
                cache[key] = P
        y = P @ y
        out[i] = y
    return out


def _csr_triple(A):
    A = sp.csr_matrix(A, dtype=complex)
    A.sort_indices()
    return (np.ascontiguousarray(A.data, dtype=np.complex128),
            np.ascontiguousarray(A.indices, dtype=np.int64),
            np.ascontiguousarray(A.indptr, dtype=np.int64))


def evolve_ensemble(L0, ops, ctrl_t, ctrl_u, alpha, beta, weights, y0, time,
                    rtol=1e-8, atol=1e-10, max_step=None, backend="rust"):
    """Weighted ensemble evolution; see the module docstring.

    Returns:
        ndarray: (len(time), N) array, sum_m weights[m] * vec(rho_m(t)).
    """
    ctrl_t = np.ascontiguousarray(ctrl_t, dtype=float)
    ctrl_u = np.ascontiguousarray(ctrl_u, dtype=float).reshape(len(ops), len(ctrl_t))
    alpha = np.ascontiguousarray(alpha, dtype=float)
    beta = np.ascontiguousarray(beta, dtype=float)
    weights = np.ascontiguousarray(weights, dtype=float)
    y0 = np.ascontiguousarray(y0, dtype=np.complex128)
    time = np.ascontiguousarray(time, dtype=float)
    N = len(y0)

    # Drop operators that never contribute
    keep = [j for j in range(len(ops)) if np.any(alpha[:, j]) or np.any(beta[:, j])]
    ops = [ops[j] for j in keep]
    ctrl_u = ctrl_u[keep]
    alpha = alpha[:, keep]
    beta = beta[:, keep]

    if backend == "rust":
        if not HAVE_RUST:
            raise RuntimeError("The Rust backend is not available; build it with `maturin develop --release` "
                               "or use backend='scipy'.")
        return _lased_rs.evolve_ensemble(
            N, _csr_triple(L0), [_csr_triple(M) for M in ops], ctrl_t, ctrl_u, alpha, beta, weights,
            y0, time, rtol=rtol, atol=atol, max_step=np.inf if max_step is None else float(max_step))
    if backend != "scipy":
        raise ValueError(f"Unknown backend {backend!r}")

    from scipy.integrate import solve_ivp
    L0 = sp.csr_matrix(L0)
    ops = [sp.csr_matrix(M) for M in ops]
    out = np.zeros((len(time), N), dtype=complex)

    def u(j, t):
        if len(ctrl_t) == 0:
            return 1.0
        return np.interp(t, ctrl_t, ctrl_u[j])

    for m in range(len(weights)):
        a, b = alpha[m], beta[m]
        static = [j for j in range(len(ops)) if a[j] == 0]
        Ls = L0.copy()
        for j in static:
            if b[j] != 0:
                Ls = Ls + b[j]*ops[j]
        Ls = Ls.tocsr()
        dyn = [j for j in range(len(ops)) if a[j] != 0]

        def rhs(t, y):
            dy = Ls @ y
            for j in dyn:
                dy += (a[j]*u(j, t) + b[j])*(ops[j] @ y)
            return dy

        if len(time) == 1:
            out[0] += weights[m]*y0
            continue
        sol = solve_ivp(rhs, (time[0], time[-1]), y0, method="DOP853", t_eval=time,
                        rtol=rtol, atol=atol, max_step=np.inf if max_step is None else max_step)
        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")
        out += weights[m]*sol.y.T
    return out
