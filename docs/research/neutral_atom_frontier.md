> **Provenance:** AI-assisted literature survey written while modernising LASED (October 2026). Some sources could only be checked through search metadata, so verify citations before relying on specific numbers.

# LASED and neutral-atom quantum computing: where the value is (2023–2026)

## TL;DR

Leading groups now report CZ fidelities of about 99.4–99.7%: Harvard/QuEra 99.5% across 60 parallel pairs, Caltech 99.62–99.71% on Sr, Atom Computing 99.40% (99.72% post-selected) on Yb-171, and Infleqtion 99.35% (99.73% post-selected) on Cs. At that level the remaining error comes from many small effects in a **multi-level** atom. These include intermediate-state scattering, Rydberg decay and leakage into other states, hyperfine and Förster structure in the Rydberg manifold, Doppler and position fluctuations, and laser phase and intensity noise. Every group builds its own error-budget code to account for these. Pulser (Pasqal) and Bloqade (QuEra) work at the level of *effective two-level, many-atom* models. ARC and pairinteraction provide atomic *data*. No open tool takes a real level structure, applies time-dependent pulses to one or two atoms, and returns a calibrated error budget. LASED already has the hard part, the angular-momentum algebra, so it can fill this gap. That gap is "the alpha".

**LASED today.** I inspected the code in `LASED/time_evolution.py`, `time_evolution_matrix.py` and `laser_atom_system.py`. It has:
- exactly two manifolds (`E`, `G`) and one laser wavelength;
- constant intensity and detuning, solved by eigen-decomposing a dense n²×n² Liouvillian built element by element with sympy;
- no magnetic field;
- no multi-atom support.

Almost every frontier feature needs a **time-dependent, multi-laser, N-manifold engine** first, so that is roadmap item #1.

## Landscape: what exists and where LASED fits

| Tool | What it does | Gap LASED can fill |
|---|---|---|
| Pulser (Pasqal) [docs](https://pulser.readthedocs.io/en/stable/noise_model.html) | Pulse-level, many-atom; effective 2–3 level; phenomenological noise (SPAM, dephasing, amplitude, Doppler, leakage "error state") | No hyperfine/Zeeman structure, no intermediate state, no microscopic decay branching |
| Bloqade (QuEra) [github](https://github.com/QuEraComputing/bloqade-python) | Analog Rydberg Hamiltonians, 2-level | Same |
| RydIQule v1/v2 [CPC 2024](https://doi.org/10.1016/j.cpc.2023.108952), [v2 arXiv:2510.21628](https://arxiv.org/pdf/2510.21628) | Graph-based multi-level OBEs, Doppler averaging, aimed at Rydberg *vapor-cell sensing* | Not aimed at tweezer gates, multi-atom blockade or gate-fidelity metrics. Closest competitor architecturally, so copy its graph/rotating-frame idea |
| pylcp (NIST) [CPC 270, 108166](https://github.com/jqiAMO/pylcp) | OBE/rate equations with B-fields and atomic motion, for laser cooling | No Rydberg, no gates |
| ARC 3.0 [CPC 261, 107814](https://arc-alkali-rydberg-calculator.readthedocs.io/) / pairinteraction v2 [PyPI](https://pypi.org/project/pairinteraction/) | Matrix elements, lifetimes, C6, pair potentials, polarisabilities | Data only, no dynamics; use as data back-ends |

**Positioning:** "the pulse-level master-equation layer between ARC/pairinteraction (data) and Pulser/Bloqade (many-body), with full fine/hyperfine/Zeeman structure and error budgets".

---

## Top capabilities

### 1. Time-dependent pulses on a general multi-laser level graph (foundation)

**Why it matters.** Every gate is a shaped pulse.
- Time-optimal CZ gates use a phase-modulated global pulse φ(t) at fixed Ω ([Jandura & Pupillo, Quantum 6, 712 (2022)](https://quantum-journal.org/papers/q-2022-05-13-712/)). These were used for the 99.5% parallel CZ of [Evered et al., Nature 622, 268 (2023)](https://doi.org/10.1038/s41586-023-06481-y) and the 99.71% Sr CZ of [Tsai et al., PRX Quantum 6, 010331 (2025)](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.6.010331).
- Amplitude- and Doppler-robust variants appear in [Fromonteil, Bluvstein & Pichler, PRX Quantum 4, 020335 (2023)](https://arxiv.org/pdf/2210.08824).
- Adiabatic rapid passage, STIRAP and quasi-adiabatic Δ(t) sweeps are standard for analog simulation ([Scholl et al., Nature 622, 273 (2023)](https://doi.org/10.1038/s41586-023-06516-4)).

**Physics.** In a rotating frame chosen by a spanning tree over the laser couplings (RydIQule's trick):

H(t) = Σ_i Δ_i(t)|i⟩⟨i| + Σ_L Σ_(a,b) [Ω_L(t)/2] e^{iφ_L(t)} c^{q_L}_{ab}|a⟩⟨b| + h.c.

- c^q_{ab} are LASED's existing 3j/6j coupling factors.
- Δ_i(t) includes Zeeman shifts, light shifts and the Doppler term −Σ k_L·v(t).

The solver needs to:
- integrate dρ/dt = −i[H(t),ρ] + Σ_k D[L_k]ρ in **matrix form**, not as an n²×n² superoperator, using scipy `solve_ivp` (DOP853) or a QuTiP/dynamiqs backend;
- handle piecewise-constant segments with `expm_multiply`;
- keep the existing eigen-solution as the fast path when the Hamiltonian is constant.

**API sketch.**
```python
sys = lased.System(levels=[g_manifold, p_manifold, r_manifold],
                   lasers=[Laser("420", couples=("5S1/2","6P3/2"), pol="sigma+", k=+z,
                                 rabi=Blackman(peak=2*pi*160e6, T=250e-9), detuning=Const(2*pi*-2e9)),
                           Laser("1013", couples=("6P3/2","70S1/2"), pol="sigma-", k=-z,
                                 rabi=FlatTop(rise=50e-9,...), phase=Spline(t_knots, phi_knots))],
                   B=Field(Bz=8e-4))
res = sys.evolve(rho0, t=np.linspace(0, 400e-9, 2001))
```
- Pulse primitives: `Gaussian`, `Blackman`, `Sin2`, `FlatTop(ramp)`, `Linear/TanhChirp`, `Spline`, `Composite([...])` (BB1, CORPSE), and `Callable`.
- Every pulse is JSON-serialisable so agents can round-trip configurations.

**Difficulty:** Medium. This is mainly a refactor of the coupling and decay builders into sparse operators. The rest of the roadmap depends on it.

### 2. Two-photon Rydberg excitation with the intermediate state kept explicitly

**Why it matters.** Rb and Cs gates go through 5P/6P (420+1013 nm for Rb, [Levine et al., PRL 123, 170503 (2019)](https://doi.org/10.1103/PhysRevLett.123.170503)), and the intermediate state is a leading error.
- Evered et al. credit part of their improvement to "atomic dark states to reduce scattering".
- [de Léséleuc et al., PRA 97, 053803 (2018)](https://arxiv.org/abs/1802.10424) showed that Rabi-oscillation damping is a sum of several small effects: intermediate scattering, Doppler, laser noise and SPAM.
- For alkaline-earth atoms, single-photon UV excitation from ³P₀ (Sr, Yb) still needs the clock and metastable structure. Hyperfine mixing in ¹⁷¹Yb Rydberg states plus an anomalous Förster resonance limited earlier gates ([Peper et al., PRX 15 (2025)](https://arxiv.org/abs/2406.01482v2)).

**Physics.**
- Ladder g→e→r with Ω₁, Ω₂ and a large Δ_e.
- Full hyperfine structure of 5P₃/₂ or 6P₃/₂ F′ levels, so polarisation-dependent dark and bright paths appear automatically.
- Hyperfine structure of the Rydberg state where it matters (Yb-171, low-n Rb).
- Differential AC Stark shifts Ω₁²/4Δ − Ω₂²/4Δ, plus shifts from off-resonant F′ levels.
- Scattering probability ≈ Γ_e·∫(Ω₁²+Ω₂²)/(4Δ²)dt.
- An `adiabatic_eliminate()` option and a built-in comparison of the eliminated model against the full model.

**API.**
- `lased.presets.Rb87.rydberg_ladder(n=53, L="S", J=0.5, via="6P3/2")`.
- Reports from a run: `res.scatter_events("6P3/2")`, `res.light_shift(qubit_levels)`.

**Difficulty:** Small once item 1 exists. This is very high value because LASED's hyperfine engine is the differentiator.

### 3. Two-atom (and small-N) blockade gates with real level structure

**Why it matters.** The CZ gate is the bottleneck in every platform:
- [Muniz et al. (Atom Computing), PRX Quantum 2025](https://arxiv.org/pdf/2411.11708): 99.40%.
- [Radnaev et al. (Infleqtion), arXiv:2408.08288](https://arxiv.org/pdf/2408.08288): 99.35%.
- [Finkelstein et al., Nature 634, 321 (2024)](https://doi.org/10.1038/s41586-024-08005-8): 99.62%.

The 2026 theory shows that a single effective pair state is the wrong model near Förster resonances; keeping the exchange channels changes predicted fidelities by up to 100× ([Norrell, Shen, Saffman & Otten, arXiv:2605.19245](https://arxiv.org/pdf/2605.19245)). Dual-species Rb–Cs gates for syndrome extraction also appear ([Miles et al., arXiv:2603.13492 (2026)](https://arxiv.org/pdf/2603.13492)).

**Physics.**
- Tensor-product Hilbert space H₁⊗I + I⊗H₂ + V, where V = Σ V_{rr′}(R)|rr′⟩⟨rr′| for van der Waals (C₆/R⁶).
- Optional pair-state matrix V_{(r1r2),(r3r4)} for Förster/dipolar channels, imported from pairinteraction.
- Per-atom collapse operators, and support for different species and different local Ω or Δ per atom.
- Symmetric-subspace reduction for global pulses.
- Gate metrics: CZ fidelity after optimising single-qubit Z phases, Bell-state fidelity, and the average gate fidelity of [Pedersen et al. 2007](https://doi.org/10.1016/j.physleta.2007.02.069).

**API.**
- `pair = lased.Pair(sys_a, sys_b, R=4e-6, interaction=C6(…)|PairInteraction(...))`
- `pair.evolve(...)`
- `lased.gates.cz_fidelity(pair, pulse, qubit_map)`
- Pulse library `lased.pulses.time_optimal_cz(Omega, V)` returns the Jandura–Pupillo phase profile.

**Difficulty:** Medium–Large. Dimension grows to n² (ρ has n⁴ entries), so this needs sparse or matrix-form propagation and optional quantum trajectories (MCWF). Keep atoms to about 10–20 levels each after truncation.

### 4. Laser noise: from a PSD to infidelity

**Why it matters.**
- Frequency noise, especially PDH **servo bumps** near the Rabi frequency, dominates over intensity noise for well-stabilised lasers ([Day et al., npj QI 8, 72 (2022)](https://doi.org/10.1038/s41534-022-00586-4)).
- A perturbative master-equation theory matches the numerics ([Jiang, Scott, Friesen & Saffman, PRA 107, 042611 (2023)](https://link.aps.org/doi/10.1103/PhysRevA.107.042611)).
- Tsai et al. 2025 introduce "fidelity response theory": a linear-response kernel I(f) such that 1−F ≈ ∫S(f)I(f)df. They use it to *learn laser noise from gate data* and to choose between protocols.
- Feed-forward phase-noise cancellation is now a hardware trend ([Denecker et al., arXiv:2411.10021](https://arxiv.org/pdf/2411.10021)), which creates demand for a "what if the servo bump is 30 dB lower?" tool.

**Physics.** Three levels of treatment:
- (a) **Quasi-static ensembles**: shot-to-shot δΩ, δΔ, position and velocity. This generalises LASED's existing Doppler and Gaussian averaging to time-dependent runs.
- (b) **Stochastic trajectories**: φ(t) or δI(t) realisations synthesised from a user PSD (white + 1/f + Lorentzian/Gaussian servo bump) by the FFT method, then ensemble-averaged.
- (c) **Linear response**: compute the sensitivity kernel I(f) of any pulse by Fourier-transforming the first-order perturbation. This is cheap and gives an interpretable budget.

**API.**
- `noise = PhaseNoisePSD(h0=…, bump=(f=1.2e6, width=…, height=…))`
- `res = sys.evolve(..., noise=[noise, IntensityNoise(rel_rms=0.005)], n_traj=200, seed=1)`
- `lased.response.kernel(sys, pulse, observable=fidelity)`

**Difficulty:** Medium. (a) Small, (b) Medium, (c) Medium.

### 5. Leakage, decay channels and error-budget bookkeeping

**Why it matters.**
- Erasure conversion is now central to fault-tolerance roadmaps. Detected leakage raises the threshold from about 0.94% to about 4.15% ([Wu, Kolkowitz, Puri & Thompson, Nat. Commun. 13, 4657 (2022)](https://doi.org/10.1038/s41467-022-32094-6)).
- Experiments: [Ma et al., Nature 622, 279 (2023)](https://doi.org/10.1038/s41586-023-06438-1) (Yb-171 metastable qubit; 33% of two-qubit errors converted to erasures), Scholl et al. 2023 (Sr), and leakage detection units in Cs ([Chow et al., PRX Quantum 2024](https://arxiv.org/abs/2405.10434v3)).
- Peper et al. attribute the remaining Yb gate error to Rydberg decay and Doppler shifts.
- Experiments need to know where population ends up: back in the qubit, in other hyperfine states, in non-qubit ground states, in other Rydberg states via BBR, ionised, or lost from the trap.

**Physics.**
- Explicit auxiliary "sink" levels with labelled jump operators:
  - radiative Rydberg decay branched through intermediate levels;
  - BBR transfer at temperature T (rates from ARC);
  - photoionisation by the tweezer;
  - an anti-trapped-Rydberg loss probability.
- Each sink is tagged `detectable_as_erasure=True/False`.
- Error budget = integrated jump-rate counts per channel, plus coherent infidelity from the residual.

**API.**
- `res.budget()` returns a dict, e.g. `{"intermediate_scatter": 1.2e-3, "rydberg_decay_radiative": 8e-4, "bbr": 3e-4, "leak_nonqubit_hf": …, "doppler": …, "phase_noise": …}` together with erasure versus Pauli classification.

This is ideal for agents: it is structured, comparable and diffable.

**Difficulty:** Small–Medium, and the value is very high.

### 6. Atom motion, thermal effects and position dependence

**Why it matters.**
- Doppler shifts and finite temperature are named limits in Peper 2025, Tsai 2025 and Miles 2026. Miles et al. use an adiabatic trap-depth ramp that trades motional error against worse localisation.
- Photon recoil and focused beams entangle internal states with motion ([Robicheaux, Graham & Saffman, PRA 103, 022424 (2021)](https://doi.org/10.1103/PhysRevA.103.022424)).
- Transport during "zoned" architectures moves atoms over hundreds of µm ([Bluvstein et al., Nature 626, 58 (2024)](https://doi.org/10.1038/s41586-023-06927-3); [arXiv:2506.20661](https://arxiv.org/abs/2506.20661)). In continuously reloaded arrays the qubits must stay coherent while others are loaded nearby ([Chiu et al., Nature 646, 1075 (2025)](https://arxiv.org/abs/2506.20660)).

**Physics (moving frame).** Sample (x₀, v₀) from a thermal distribution in a harmonic trap. A trajectory x(t) = x₀cos ωt + (v₀/ω)sin ωt, or ballistic motion when the trap is off, enters the model as:
- Δ(t) ← Δ − k_eff·v(t), with counter-propagating beams reducing k_eff;
- Ω(t) ← Ω·u(x(t)) from the Gaussian beam profile;
- δ_LS(t) ∝ Δα·U(x(t), t), the differential light shift along a transport trajectory U(t).

The first stage is a classical-trajectory ensemble, which extends LASED's existing averaging. A later stage adds an optional quantised motional mode in the Lamb–Dicke expansion for recoil, Raman sideband cooling and motional heating.

**API.**
- `Thermal(T=10e-6, trap=Harmonic(omega_r, omega_z))`
- `Trajectory.from_tweezer_move(path, U(t))`
- `sys.evolve(..., motion=…)`

**Difficulty:** Small–Medium for the classical version; Large for the quantised mode.

### 7. Static fields, Landé factors and trap light shifts

**Why it matters.**
- Qubit frequencies and gate phases depend on Zeeman shifts. Atom Computing works at 500 G ([Norcia et al., PRX 13, 041034 (2023)](https://doi.org/10.1103/PhysRevX.13.041034)), which is beyond the linear-Zeeman regime for many states.
- Tweezer **vector and tensor** light shifts cause inhomogeneous dephasing.
- Magic wavelengths decide the clock and omg qubits ([Lis et al. (Kaufman), PRX 13, 041035 (2023)](https://doi.org/10.1103/PhysRevX.13.041035)).

**Physics.**
- Full Hamiltonian H_hfs + μ_B(g_J J + g_I I)·B, diagonalised per manifold (Breit–Rabi to Paschen–Back), with lasers then coupling the eigenstates.
- Tweezer shift −U·[α_s + α_v (C) m_F/(2F) + α_t (…)], with polarisabilities supplied as data or computed by ARC.

**API.**
- `System(..., B=Field(…), traps=[Tweezer(λ=813e-9, U0=…, pol=…)])`
- `lased.atomic.landé_g(state)`

**Difficulty:** Medium.

### 8. Λ/Raman qubits, optical pumping, readout and erasure imaging

**Why it matters.**
- Hyperfine qubit gates are made with Raman Λ systems.
- State preparation by optical pumping limits SPAM.
- Fast, low-loss state-selective imaging and mid-circuit measurement need scattering-rate and heating budgets ([Norcia 2023](https://doi.org/10.1103/PhysRevX.13.041034); [Lis 2023](https://doi.org/10.1103/PhysRevX.13.041035); Finkelstein 2024 ancilla readout).
- Dark states, CPT and Λ-enhanced gray molasses also come from multi-level coherent dynamics.

**Physics.** These are already within reach of LASED's formalism once multiple lasers are allowed:
- two-frequency Λ (effective Raman Ω, spontaneous Raman scattering rates);
- repumpers and the steady state;
- photons-per-time and recoil-heating estimates.

**API.**
- `res.photon_scatter_rate()`, `res.steady_state()`
- `lased.protocols.optical_pumping_fidelity(...)`
- `imaging_budget(exposure, NA, collection_eff)`

**Difficulty:** Small–Medium.

### 9. Atomic-data layer

Agents need trustworthy numbers with provenance. Proposed design:
- `lased.data` adapters to [ARC 3.0](https://arc-alkali-rydberg-calculator.readthedocs.io/): alkali Rydberg matrix elements, lifetimes including BBR, C₆; Sr-88, Ca-40 and Yb-174 single-electron.
- Adapters to [pairinteraction](https://pypi.org/project/pairinteraction/) for pair potentials and Förster channels.
- Curated presets with citations for D lines and hyperfine constants: [Steck alkali data](https://steck.us/alkalidata/), [NIST ASD](https://physics.nist.gov/asd), and the [UD high-precision portal](https://www1.udel.edu/atom/) for polarisabilities.
- Yb-171 Rydberg data from the MQDT models of Peper et al.
- He (LASED's historical use case) and Ca.

Priority presets: Rb87 (QuEra, Pasqal, Lukin, Browaeys), Cs133 (Infleqtion, Saffman), Sr88/Sr87 (Endres, Kaufman, planqc), Yb171/Yb174 (Atom Computing, Thompson, Kaufman).

**Difficulty:** Medium, mostly curation.

### 10. Benchmarks, fitting and calibration (agent-facing)

Provide:
- fits for damped Rabi oscillations and Ramsey/echo data;
- metrics: Bell, CZ and average gate fidelity, plus an RB-like decay estimate;
- calibration by gradient-free or gradient-based fits of (Ω, Δ, T) to data.

Muniz et al. highlight "a simple and optimized method for calibrating multi-parameter quantum gates". Ship a test suite that reproduces published curves:
- de Léséleuc 2018 Rabi damping;
- Jandura–Pupillo CZ fidelity versus Ω/V;
- Day 2022 infidelity versus servo bump.

These double as worked examples for agents. A JAX/dynamiqs backend would add autodiff for optimal control (GRAPE) and fidelity-response kernels. **Difficulty:** Small for fitting; Large for differentiability.

---

## Agent-oriented design notes

- **Declarative, serialisable configs** (pydantic/JSON) with explicit units (rad/s vs Hz) and named presets carrying citations.
- **Structured outputs**: labelled arrays (xarray), the `budget()` dict, and validity checks (trace = 1, positivity, convergence of averages, warnings when the RWA or adiabatic elimination breaks down).
- **Deterministic seeds**, cheap default resolutions, and `estimate_cost()` before large two-atom runs.
- **Interop**: `to_qutip()` and `to_pulser_sequence()` exporters so agents can hand off to many-body tools.

---

## Prioritised roadmap (value ÷ effort)

| Rank | Feature | Value | Effort | Notes |
|---|---|---|---|---|
| 1 | **Time-dependent engine + N-manifold, multi-laser level graph + pulse library** (Ω(t), Δ(t), φ(t), chirps, composite pulses; matrix-form ODE, expm fast path) | ★★★★★ | M | Unlocks everything; fixes dense n²×n² sympy bottleneck |
| 2 | **Two-photon Rydberg ladder presets** with explicit intermediate hyperfine structure, light shifts, scattering, adiabatic-elimination comparison | ★★★★★ | S | Unique selling point versus Pulser/Bloqade |
| 3 | **Leakage/decay sinks + `budget()`** (radiative branching, BBR, erasure tagging) and fidelity metrics | ★★★★★ | S–M | Agent-friendly; directly answers the erasure-conversion design question |
| 4 | **Noise module**: quasi-static ensembles (time-dependent Doppler/intensity/position), PSD → stochastic φ(t)/I(t), linear-response kernel | ★★★★☆ | M | Day 2022 / Jiang 2023 / Tsai 2025 workflows |
| 5 | **Two-atom blockade**: tensor product, C₆ or pairinteraction pair matrix (Förster channels), CZ fidelity, time-optimal pulse templates, dual-species | ★★★★★ | M–L | Most-wanted capability; needs #1 and #3 |
| 6 | **Atomic data layer** (ARC/pairinteraction adapters; Rb87, Cs, Sr88/87, Yb171/174, Ca, He presets with provenance) | ★★★★☆ | M | Makes agents reliable |
| 7 | **B-field (Breit–Rabi/Paschen–Back) + tweezer vector/tensor light shifts + thermal-trajectory motion/transport** | ★★★★☆ | M | Clock/omg qubits, 500 G operation, transport dephasing |
| 8 | **Λ/Raman, optical pumping, readout/imaging and mid-circuit-measurement budgets; fitting/calibration utilities** | ★★★☆☆ | S–M | Extends the existing steady-state OBE strengths |
| (next) | Differentiable JAX/dynamiqs backend for GRAPE; quantised motion (Lamb–Dicke) for recoil and sideband cooling; MCWF for 3+ atoms | ★★★★☆ | L | Longer-term |

**Suggested sequencing:**
1. Items 1 and 3 together, as a "LASED 2.0 core".
2. Item 2 with Rb87/Cs presets, validated against de Léséleuc 2018.
3. Items 4 and 6.
4. Item 5, validated against Jandura–Pupillo and Evered/Tsai.
5. Items 7 and 8.

Items 1–3 alone would let LASED answer this question with citations: "What is the intermediate-state-limited fidelity of a 420/1013 nm Rb π-pulse to 53S at Δ = 2 GHz, and how much of the error is erasure-detectable?" Today no open package can.

*Note: arxiv.org could not be fetched directly in this session. Citations were verified through search-result metadata, and a few numbers (e.g. 0.9971 for Tsai et al.) come from abstracts or conference posters.*
