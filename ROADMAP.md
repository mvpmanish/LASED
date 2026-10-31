# LASED roadmap

This roadmap comes from the survey in [docs/research/neutral_atom_frontier.md](docs/research/neutral_atom_frontier.md)
(what neutral-atom QC and atomic-physics groups need from a multi-level master-equation code, 2023–2026)
and [docs/research/elliptical_and_simulators.md](docs/research/elliptical_and_simulators.md)
(numerical methods, other simulators, Rust tooling).

**Positioning:** LASED is the pulse-level, multi-level master-equation layer. It sits between atomic-data
tools (ARC, pairinteraction) and many-body tools (Pulser, Bloqade), with full fine, hyperfine and
Zeeman structure and error budgets.

## Done in v2.0

- [x] Hermitian Hamiltonian and Lindblad decay. This fixes elliptical light and decay coherence transfer.
- [x] **Time-dependent multi-laser engine** (`AtomicSystem`, `Laser`): Ω(t), φ(t) and Δ(t) per laser,
  any level graph without closed loops, automatic rotating frame. This was roadmap item #1, which
  everything else depends on.
- [x] Pulse library: Gaussian, square, sin², Blackman, flat-top, sampled, phase, chirps, ARP.
- [x] Ensembles (Doppler, Gaussian beam) run in parallel in Rust. Moving-atom trajectories become
  Ω(t) and Δ(t).
- [x] Laser linewidth (white frequency noise), dephasing, loss channels, steady state.
- [x] Linear Zeeman shifts and Landé g-factors.
- [x] Rust backend, packaging with `pyproject.toml` and maturin, tests (including against QuTiP),
  `AGENTS.md`.
- [x] Array-based backend: `Basis` (quantum-number arrays), vectorised Wigner 3j/6j and dipole blocks,
  no loops over matrix elements; exact SciPy propagator for all continuous-wave problems; sparse
  steady-state solve; reproducible backend benchmark.

## Next (ranked by value ÷ effort)

1. **Two-photon Rydberg presets.** Rb 420/1013 nm and 780/480 nm, Cs 459/1039 nm, Sr 3P0→Rydberg.
   Keep the intermediate hyperfine structure explicit, and report light shifts and scattering against
   the adiabatically eliminated model. *(The engine already supports this; see
   `examples/rydberg_two_photon.py`. What's missing is data presets and convenience functions.)*
2. **Leakage sinks and `budget()`.** Labelled sink levels for radiative decay, blackbody transfer,
   photoionisation and trap loss, each tagged as erasure-detectable or not. A structured error-budget
   dict, plus fidelity metrics.
3. **Noise module.** Turn a phase or intensity noise spectrum (servo bumps) into stochastic φ(t)/I(t)
   traces averaged over an ensemble, plus a linear-response infidelity kernel (Day 2022, Jiang 2023,
   Tsai 2025).
4. **Two-atom blockade gates.** Tensor-product systems with a C6 interaction or pair-state channels,
   CZ fidelity, and time-optimal gate pulse templates (Jandura & Pupillo 2022, Evered 2023).
5. **Atomic data layer.** ARC and pairinteraction adapters, and curated presets with sources for Rb87,
   Cs, Sr, Yb, Ca and He.
6. **Fields and motion.** Breit–Rabi / Paschen–Back treatment, tweezer vector and tensor light shifts,
   and thermal motion in a trap.
7. **Raman, optical pumping, readout and imaging.** Mid-circuit-measurement budgets, plus Rabi/Ramsey
   fitting helpers.
8. **Numerics.** A piecewise-constant (Magnus) propagator for stiff pulsed problems with large hyperfine
   splittings, a differentiable JAX backend for pulse optimisation, exact real (Hermitian-packed)
   storage, and Krylov `expm_multiply` for very large systems.
