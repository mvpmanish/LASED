Backends and performance
========================

``AtomicSystem.evolve`` / ``LaserAtomSystem.timeEvolution`` take ``backend=``:

``"expm"``
    Exact propagation with matrix exponentials (LAPACK through SciPy). This is the default for
    continuous-wave problems with a single atom.
``"rust"``
    Compiled Rust extension (PyO3 + maturin). Adaptive Dormand–Prince 5(4) on sparse superoperators,
    with ensemble members run in parallel by rayon. This is the default for pulses and ensembles.
``"scipy"``
    Pure-Python fallback (``scipy.integrate.solve_ivp``), used when the extension is not built.

Before solving, LASED drops the density-matrix elements that cannot be reached from the initial state.
This is exact. For π light and a diagonal initial state, only the :math:`\Delta m = 0` elements remain.

Check what is available with ``LASED.backendInfo()``. Benchmarks are in ``benchmarks/``. On the
Cs D2 tutorial (48 sub-states, 2304 equations), v1 did not finish after 35 minutes; v2 takes a few
seconds. For a shaped, chirped pulse on Rb-87 D2, the Rust backend is about 160 times faster than the
SciPy fallback.
