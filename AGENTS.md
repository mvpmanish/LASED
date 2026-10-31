# AGENTS.md: guide for AI agents (and humans in a hurry)

LASED simulates lasers driving multi-level atoms (fine, hyperfine and Zeeman structure) by solving the
Lindblad master equation for the density matrix. This file has what you need to use it correctly
without reading the source. Everything below is tested in `tests/`.

## Install and check

```bash
pip install LASED                     # wheels include the compiled Rust backend
# from a source checkout (needs a Rust toolchain: https://rustup.rs)
pip install maturin && maturin develop --release
python -c "import LASED; print(LASED.backendInfo())"   # {'rust': True, 'rust_threads': ...}
```

Without the Rust extension, everything still works through the SciPy fallback. Continuous-wave
problems are unaffected (they always use SciPy). Pulsed problems are about 3× slower, and pulsed
ensembles about 10× slower.

## Units (always)

| Quantity | Unit |
|---|---|
| time | ns |
| angular frequency, detuning, Rabi frequency, rate, linewidth | Grad/s = rad/ns (so `2*np.pi*0.1` = 100 MHz) |
| wavelength | m |
| intensity | mW/mm² |
| power | mW |
| lifetime `tau` | ns |
| velocity | m/s |
| magnetic field | gauss |

## Conventions

* **States.** Every `State` has `label` (1..n, contiguous, unique), `w` (angular frequency in Grad/s),
  quantum numbers `L, S, J, I, F, m`. The matrix index of a state is `label - 1`.
  `generateSubStates(label_from, w, L, S, J=None, I=None, F=None)` makes all `m` sub-states of one
  level (half-integers allowed).
* **Density matrix layout.** `rho_t_matrix[k, i-1, j-1]` = ρ_ij at `time[k]`. The legacy flat layout
  `rho_t[index(i, j, n), k]` uses row-major order, `index = (i-1)*n + (j-1)`.
* **Polarisation.** The quantisation axis is z. A polarisation is a dict `{q: a_q}` of spherical
  amplitudes, where `a_q` drives Δm = m_upper − m_lower = q (+1 is σ+, 0 is π, −1 is σ−). It is
  defined by ε = Σ a_q u_q with u_±1 = (x ± iy)/√2 and u_0 = z, so `a_q = u_q* · ε`. Complex values
  are fine (elliptical light). Use the helpers instead of writing them by hand:
  `las.linear(angle, polar=pi/2)`, `las.circular(±1)`, `las.elliptical(angle, ellipticity)`,
  `las.jones(Ex, Ey, Ez)`, `las.pi_polarised()`.
  The legacy `Q=[...]` plus `rabi_factors=[...]` means `{Q[i]: rabi_factors[i]*rabi_scaling}`.
* **Rabi frequency.** `Laser(rabi=Ω0)` sets the *reduced* Rabi frequency. Sub-transition u←l has
  coupling `V[u,l] = (Ω0/2)·A(t)·Σ_q a_q·coupling(u, l, q)`, and `V[l,u] = conj(V[u,l])`, so the
  Hamiltonian is always Hermitian. To set a *specific* sub-transition's Rabi frequency Ω, use
  `rabi = Ω/abs(las.coupling(u, l, q))`. Alternatively give `intensity` and `tau`; then
  Ω0/2 = `halfRabiFreq(intensity, tau, wavelength)`.
* **Detuning sign.** `detuning > 0` means the laser frequency is above `angularFreq(wavelength)`
  (blue). Resonance is set by the states' `w`: the laser is resonant when
  `angularFreq(wavelength) + detuning == w_upper - w_lower`.
* **Doppler.** An atom with velocity v sees the laser shifted by −k·v, where k = (2π/λ)·`direction`.
  A scalar velocity points along +z.
* **Rotating frame.** The frame is chosen automatically, one laser frequency per edge of the laser
  graph. Lasers must not form a closed loop (you get `ValueError("... closed loop ...")`). If two
  `Laser` objects describe the same beam, merge them into one with a combined polarisation.
* **Decay.** `Decay(upper, lower, tau, branching=1)` uses Lindblad jump operators. Branching follows
  the squared couplings, normalised for each upper state, so each upper state decays at
  `branching/tau` into `lower`. Coherences are transferred (including Zeeman and hyperfine ones).
  `Loss(states, tau)` removes population from the model, so the trace decreases.
  `Dephasing(states, rate)` and `Laser(linewidth=Γ)` (FWHM, white frequency noise) add dephasing.

## API cheat sheet

```python
import numpy as np, LASED as las

# --- Classic API: one laser, ground states G -> excited states E, decay E -> G
sys = las.LaserAtomSystem(E, G, tau=26.2, laser_wavelength=780.24e-9, laser_intensity=1.0,
                          polarisation=las.elliptical(angle=0.3, ellipticity=0.2),
                          tau_f=None, tau_b=None)            # optional extra losses (ns)
sys.setRho_0(g, g, 1.0)                    # or appendDensityMatrixToRho_0(rho_g, "g"); default: equal ground pops
rho = sys.timeEvolution(time,              # returns (T, n, n); also stored in sys.rho_t
        detuning=2*np.pi*0.01, atomic_velocity=0.0,
        pulse=las.pulses.gaussian(t0=50, fwhm=20),          # complex envelope A(t) (phase allowed)
        detuning_t=las.pulses.linearChirp(rate=1e-3, t0=50),# delta(t), Grad/s
        linewidth=0.0, backend="auto")
sys.population(state); sys.excited_population; sys.rho_t_matrix; sys.steadyState()
sys.toAtomicSystem()                        # -> general engine object (add lasers, inspect H)

# --- General engine: any states, any number of lasers
red  = las.Laser(lower=S, upper=P, wavelength=780.24e-9, rabi=2*np.pi*0.1, polarisation=las.circular(+1),
                 detuning=2*np.pi*2.0, envelope=None, detuning_t=None, direction=(0, 0, 1), linewidth=0.0)
blue = las.Laser(lower=P, upper=R, wavelength=479.6e-9, rabi=2*np.pi*0.1, polarisation=las.circular(-1),
                 detuning=-2*np.pi*2.0, direction=(0, 0, -1))
atom = las.AtomicSystem(S + P + R, lasers=[red, blue],
                        decays=[las.Decay(P, S, tau=26.24)], losses=[las.Loss(R, tau=150e3)],
                        dephasings=[])
rho = atom.evolve(time, rho0=None, ensemble=None, velocity=0.0, backend="auto",
                  rtol=1e-8, atol=1e-10, max_step=None, n_ctrl=4001)    # (T, n, n)
atom.steadyState(); atom.hamiltonian(t=None, velocity=0.0); atom.liouvillian(); atom.describe()  # JSON-able dict
atom.stateTable()        # the sub-states as a table (pandas DataFrame if installed)
atom.build()             # Model: H0, couplings V, detuning generators D, collapse ops, superoperators

# --- Ensembles (run in parallel in Rust); combine with *
ens = las.dopplerEnsemble(atom.lasers, temperature=10e-6, mass=87*1.6605e-27, n=41) \
      * las.gaussianBeamEnsemble(n_lasers=2, n=10, r_max=2.0)
rho = atom.evolve(time, ensemble=ens)

# --- Pulses (las.pulses): gaussian, square, sinSquared, blackman, flatTop, fromSamples, constant;
#     p1*p2, p1+p2, p.withPhase(phi or phi(t)), p.delayed(dt), p.area(rabi, t0, t1)
#     detunings: linearChirp, tanhChirp; phaseFromDetuning; movingAtom(laser, trajectory, waist)
# --- Zeeman: las.applyZeemanShift(states, B_gauss); las.landeGF(L, S, J, I, F)
```

## Recipes

Each of these is a runnable script in `examples/`:

* `examples/elliptical_light.py`: elliptical light on Rb-87 D2 (the case that used to blow up).
* `examples/pulses_and_chirps.py`: Gaussian π pulse, Ramsey phase scan, adiabatic rapid passage.
* `examples/rydberg_two_photon.py`: 5S→5P→70S ladder, effective Rabi frequency, scattering error budget.
* `examples/doppler_and_moving_atoms.py`: thermal dephasing, atom flying through a focused beam.

## Choosing a backend

`backend="auto"` picks:
- `"expm"` (SciPy's exact matrix-exponential propagator) whenever the lasers are constant, including
  for ensembles. This is the fastest option.
- `"rust"` (compiled adaptive Dormand–Prince 5(4), parallel over ensemble members) when anything is
  time-dependent (envelope, `detuning_t`, moving atom).
- `"scipy"` for time-dependent problems if the Rust extension is missing.

Timings for each backend: `python benchmarks/bench_backends.py`.

Before solving, LASED restricts the problem to the elements of ρ that are reachable from ρ0. This is
exact, and for π light it is often 10× smaller.

Accuracy knobs:
- `rtol`/`atol` set the adaptive integrator tolerances.
- `max_step` caps the step size. Set it smaller than any feature of your pulse if the pulse has sharp
  edges.
- `n_ctrl` is the number of samples used to tabulate envelopes and detunings. They are linearly
  interpolated, so raise it for long sequences with fast features.

## Validate every result

```python
herm = np.abs(rho - rho.conj().transpose(0, 2, 1)).max()       # ~1e-15
trace = np.trace(rho, axis1=1, axis2=2).real                    # 1 (less if you used Loss / tau_f / tau_b)
min_eig = min(np.linalg.eigvalsh(r).min() for r in rho)        # >= -1e-10
```

Also check convergence: halve `rtol`, double `n_ctrl` or the ensemble size, and confirm the
observable doesn't change.

## Gotchas

* Labels must be exactly 1..n. When you build a subset of states, relabel them.
* `State.w` sets the energy levels. Apply Zeeman shifts to the states *before* building the system.
* The laser wavelength sets the frame and the Doppler shift. The *states'* `w` set the resonance.
  When you build a system from transition wavelengths, use the same number for both (or use `detuning`).
* `tau_f`/`tau_b` (and `Loss`) are losses, so the trace decreases and `steadyState()` is trivially zero.
* Legacy `doppler_averaging`/`beam_profile_averaging` keep the v1 weighting. For new work use
  `dopplerEnsemble` / `gaussianBeamEnsemble`, which are properly normalised.
* `LaserAtomSystem.timeEvolution` prints "Populating ground states equally..." when `rho_0` is empty.
* Square pulses are tabulated, so their edges have a rise time of about (duration of the time grid)/`n_ctrl`.

## Repository map

```
base-LASED/LASED/       Python package
  engine.py             AtomicSystem, Laser, Decay, Loss, Dephasing, Ensemble (the core)
  basis.py              Basis: the states stored as quantum-number arrays (doubled integers)
  wigner.py             vectorised Wigner 3j/6j symbols and dipole coupling blocks
  laser_atom_system.py  LaserAtomSystem (classic API, wraps the engine)
  polarisation.py       spherical-basis polarisation helpers and convention
  pulses.py             envelopes, chirps, moving atoms
  half_rabi_freq.py     coupling() for one pair of states, couplingExact() (sympy reference), halfRabiFreq
  _backend.py           expm / Rust / SciPy solvers
  symbolic_print.py     sympy equations of motion (same equations as the numerics)
  zeeman.py, rotation.py, angular_shape.py, save_to_csv.py, ...
base-LASED/src/lib.rs   Rust backend (PyO3): parallel adaptive sparse Lindblad integrator
tests/                  pytest suite, including an independent QuTiP cross-check
examples/               runnable scripts
benchmarks/             performance scripts
docs/                   Sphinx docs (readthedocs) and research notes (docs/research)
```

Developer loop: `maturin develop --release && pytest -q`. The Rust code is checked with `cargo clippy`
(run it from `base-LASED/`).
