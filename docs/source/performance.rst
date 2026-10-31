Backends and performance
========================

``AtomicSystem.evolve`` / ``LaserAtomSystem.timeEvolution`` take ``backend=``:

``"expm"``
    Exact propagation with matrix exponentials (LAPACK through SciPy). This is the default for
    every continuous-wave problem, including ensembles.
``"rust"``
    Compiled Rust extension (PyO3 + maturin). Adaptive Dormand–Prince 5(4) on sparse superoperators,
    with ensemble members run in parallel by rayon. This is the default for time-dependent problems
    (pulses, chirps, moving atoms).
``"scipy"``
    ``scipy.integrate.solve_ivp`` on the same sparse operators, used for time-dependent problems when
    the extension is not built.

Before solving, LASED drops the density-matrix elements that cannot be reached from the initial state.
This is exact. For π light and a diagonal initial state, only the :math:`\Delta m = 0` elements remain.

Check what is available with ``LASED.backendInfo()``.

Where the time goes
-------------------

In v1 almost all of the run time went into *building* the equations: Python loops over matrix elements,
each calling exact (sympy) Wigner symbols. For the He tutorial, building the matrix took 4.1 s and the
LAPACK solve took 0.15 s. v2 stores the sub-states as quantum-number arrays
(:class:`LASED.basis.Basis`) and builds every operator with vectorised Wigner symbols
(:mod:`LASED.wigner`) and sparse matrices. The 48-state Cs D2 system assembles in about 13 ms.

Which backend to use
--------------------

Timings from ``benchmarks/bench_backends.py`` on an idle 4-core machine:

=============================================  =========  ==========
Problem                                        Backend    Time (s)
=============================================  =========  ==========
Cs D2 (48 states), CW π light, 500 ns          expm       0.16
\                                              rust       2.0
Cs D2, CW elliptical light                     expm       1.1
\                                              rust       10.7
Cs D2, CW, 16 Doppler classes                  expm       2.2
\                                              rust       8.6
Cs D2, Gaussian pulse + chirp                  rust       0.99
\                                              scipy      4.5
Rb-87 D2, shaped + chirped pulse               rust       0.35
\                                              scipy      1.1
Rb-87 D2 pulse, 16 Doppler classes             rust       1.5
\                                              scipy      16.7
=============================================  =========  ==========

For continuous-wave problems, SciPy's matrix exponential is the right tool and Rust adds nothing. For
time-dependent problems the Rust stepper avoids calling Python at every step and runs ensemble members
in parallel. Very stiff problems with large hyperfine splittings and long pulses are where explicit
stepping is slowest; a piecewise-constant propagator for these is on the roadmap.
