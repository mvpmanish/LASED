> **Provenance:** AI-assisted literature survey written while modernising LASED (October 2026). Some sources could only be checked through search metadata, so verify citations before relying on specific numbers.

# LASED: elliptical-polarisation bug, robust Lindblad numerics, simulator survey, Rust tooling

How this was researched: web search, plus reading the source of three wheels downloaded from PyPI: pylcp 1.0.2, rydiqule 2.1.3 and centrex-tlf 0.3.1. The bugs described below were reproduced numerically on LASED v1.0.1 against QuTiP 5.3.1; the equivalent checks now live in `tests/test_against_qutip.py`. arxiv.org, steck.us and readthedocs.io could not be fetched directly, so claims about those sources rest on search snippets and standard textbook results.

**Status in LASED v2.0:** Bug A and Bug B below are both fixed. The equations are now built from a Hermitian Hamiltonian and Lindblad jump operators (`LASED/engine.py`).

---

## 0. Findings checked against the actual LASED code

### Bug A (the elliptical bug): confirmed, and the exact lines are known

In `LASED/time_evolution_matrix.py` the coupling is `rabi*rabi_factors[i]*coupling(e,g,q)`. That is H_eg = Ω_eg = rabi·Σ_q ε_q c(e,g,q). It is used unconjugated everywhere. Starting from dρ/dt = −i[H,ρ] with H_ge = Ω_eg*, the correct terms and their status in the code are:

| Row | Correct term | Code at HEAD |
|---|---|---|
| ρ_gg'' | `+i Ω_eg'' ρ_ge` | OK |
| ρ_gg'' | `−i Ω_eg* ρ_eg''` | **missing conj** (`rho_ggpp`, 2nd `for e in E` loop) |
| ρ_ee'' | `−i Ω_eg ρ_ge''` | OK |
| ρ_ee'' | `+i Ω_e''g* ρ_eg` | **missing conj** (`rho_eepp`, 1st loop) |
| ρ_ge | `−i Ω_e'g* ρ_e'e` and `+i Ω_eg'* ρ_gg'` | **both missing conj** (`rho_ge`, both loops) |
| ρ_eg | `−i Ω_eg' ρ_g'g` and `+i Ω_e'g ρ_ee'` | OK |

The minimal fix is `rabi_factors[i]` → `np.conj(rabi_factors[i])` at those 4 sites. The `symbolic_print.py` equations need the same change.

Measured result (J=0→J=1, and an F=1→F'=2 system like ⁸⁷Rb D2; elliptical light with ψ=30°, χ=20°, so ε₊₁ = −0.785+0.453i and ε₋₁ = 0.366+0.211i; propagated with expm of A over 200 ns):

| System | Version | max\|ρ−ρ†\| | min eig(ρ) | max Re eig(L) | max\|ρ−ρ_QuTiP\| |
|---|---|---|---|---|---|
| J=0→1 | HEAD | 4.3e3 | −2.5e3 | +4.5e−2 | 2.4e3 |
| J=0→1 | conj fix | 1e−15 | −2e−16 | 1e−18 | **4.6e−8** |
| F=1→2 | HEAD | 1.9e2 | −1.1e2 | +3.5e−2 | 1.1e2 |
| F=1→2 | conj fix | 1e−15 | −9e−17 | 5e−17 | 5.6e−2 (see Bug B) |

An important point for the tests: **the trace stays exactly 1 even with the bug**. The population rows of A still sum to zero over each column. What blows up is ρ−ρ†: L has eigenvalues with positive real part, so populations and coherences grow without bound, with opposite signs, while their sum stays 1. A trace-only test will not catch this bug. Test Hermiticity preservation (`P L* P == L`, where P swaps (i,j)↔(j,i)), max Re eig(L) ≤ 0, and the eigenvalues of ρ(t).

Real factors (linear x/y light with the axes along x/y, or pure σ±) hide the bug, because then conj(ε)=ε. It shows up as soon as a relative phase is present, such as a rotated linear axis or ellipticity.

### Bug B (separate issue): Zeeman coherence transfer by spontaneous emission is missing

Even after the conj fix, the F=1→F'=2 result disagrees with QuTiP by about 0.1 for **real** x-linear light. Pure π or pure σ⁺ light agrees to 1e−8. Inspecting A shows why. The full Lindblad dissipator with jump operators L_q = √(Γ/N) Σ c(e,g,q)|g⟩⟨e| feeds ρ_{g(−1)g(+1)} from ρ_{e(−2)e(0)}, ρ_{e(−1)e(1)} and ρ_{e(0)e(2)}, with coefficients 1.6e−2, 1.9e−2 and 1.6e−2 ns⁻¹. LASED has **0** for all three.

The cause is in `rho_ggpp`. The `sum_decay_channels` gate needs both excited states to couple to the same g' with the same q, which only happens when m_e' = m_e''. On top of that, `generalisedDecayConstant(ep, epp, gpp, …)` ignores g. So coherence transfer with Δm ≠ 0 is dropped, and the hyperfine (same-m) branch adds a term independent of g. This matters for any configuration that mixes σ⁺ and σ⁻, including every linear polarisation not along the quantisation axis. It is the mechanism behind EIA and alignment-to-orientation effects.

Fix: build the dissipator from jump operators (section 5). The coherence transfer then comes out automatically, with the correct signs.

---

## 1. Formalism for arbitrary (elliptical) polarisation

**Spherical basis.** This is the convention of Steck (*Quantum and Atom Optics*; [Rb-87 D-line data](https://steck.us/alkalidata/)), Auzinsh–Budker–Rochester (*Optically Polarized Atoms*, OUP 2010) and Metcalf & van der Straten:
ê₊₁ = −(x̂ + iŷ)/√2, ê₀ = ẑ, ê₋₁ = (x̂ − iŷ)/√2. These satisfy ê_q* = (−1)^q ê₋q and ê_q*·ê_q' = δ_qq'.
For any vector, A_q ≡ ê_q·A and A = Σ_q (−1)^q A_q ê₋q = Σ_q A_q ê_q*. The dot product of two vectors is A·B = Σ_q (−1)^q A_q B₋q.

**Field and Rabi frequency.** Write E(t) = (E₀/2)[ε e^{−iωt} + c.c.] and expand the *polarisation vector* as ε = Σ_q ε_q ê_q, so ε_q = ê_q*·ε. Then d·ε = Σ_q ε_q d_q, where d_q = ê_q·d. In the RWA:

H_eg = −(E₀/2)⟨e|d·ε|g⟩ = −(E₀/2) Σ_q ε_q ⟨e|d_q|g⟩, nonzero only for m_e = m_g + q (Wigner–Eckart: ⟨F_e m_e|d_q|F_g m_g⟩ ∝ (−1)^{F_e−m_e} (F_e 1 F_g; −m_e q m_g)).
**H_ge = H_eg\*** holds by definition. Equivalently, ⟨e|d_q|g⟩* = (−1)^q ⟨g|d₋q|e⟩.
H = ħ Σ_{e,g} (Ω_eg|e⟩⟨g| + Ω_eg*|g⟩⟨e|) + Σ Δ-terms is Hermitian by construction. Global phases of ε drop out. **Relative phases between the ε_q set the phases of the Zeeman coherences.** This is how the orientation of the linear-polarisation axis and the ellipticity enter.

LASED's `coupling(e,g,q)` uses (F_e 1 F_g; −m_e q m_g), so q=+1 means Δm=+1 (σ⁺ relative to the quantisation axis) and `rabi_factors` play the role of ε_q. The docstring labels "+1 right-hand circular", which is ambiguous because optics and helicity conventions differ. Rename to σ±.

**Elliptical light along z** (major axis at angle ψ, ellipticity angle χ, ε = R_z(ψ)(cosχ x̂ + i sinχ ŷ)):
ε₊₁ = −e^{−iψ}(cosχ + sinχ)/√2, ε₋₁ = e^{+iψ}(cosχ − sinχ)/√2, ε₀ = 0.
Checks: x-linear gives (−1/√2, +1/√2); y-linear gives (i/√2, i/√2); χ=π/4 gives pure σ⁺.
For k not along z, compute ε in Cartesian form and take ε_q = ê_q*·ε directly. Alternatively, rotate the spherical components with the rank-1 Wigner D-matrix, ε'_q = Σ_q' D¹_{q q'}(α,β,γ) ε_q', which LASED already has in `rotation.py`. A rotation by ψ about z multiplies ε_q by e^{−iqψ}.

**References and code that do this correctly:**
- [pylcp](https://github.com/JQIamo/pylcp) ([paper](https://arxiv.org/abs/2011.07979)). Accepts `pol` in a Cartesian or spherical basis. In `obe.py` it builds *separate* superoperators for d_q and d_q*, multiplied by E_q and E_q*: `drhodt -= (-1)**q*Eq[2-ii]*(ev_mat['d_q'][ii]@rho)` and `drhodt -= (-1)**q*conj(Eq[2-ii])*(ev_mat['d_q*'][ii]@rho)`.
- [RydIQule](https://github.com/QTC-UMD/rydiqule). `hamiltonian[idx] = cc*Ω*exp(1j*phase)/2; hamiltonian[conj_idx] = np.conj(hamiltonian[idx])`.
- [CeNTREX-TlF](https://github.com/ograsdijk/CeNTREX-TlF). Builds the symbolic H with `H[i,j] += val*…; H[j,i] += conjugate(val)*…`.
- ADM ([AtomicDensityMatrix](http://rochesterscientific.com/ADM/), Mathematica). Builds H from a field vector in the spherical basis, with polarisation parametrised by ellipticity and angle as in the Auzinsh–Budker–Rochester book.

**Lesson:** never write the h.c. half of the equations by hand. Build H (n×n complex) and set H_ge = conj(H_eg) once, or build H then do `H = H + H.conj().T`. Then derive the superoperator from H mechanically.

---

## 2. Numerical techniques that preserve structure

- **Vectorised Liouvillian.** For row-major vec, which is what LASED's `index(i,j,n)=(i−1)n+(j−1)` uses: vec(AρB) = (A⊗Bᵀ)vec(ρ). So
  L_H = −i(H⊗I − I⊗Hᵀ), and D[L] = L⊗L* − ½(L†L⊗I) − ½(I⊗(L†L)ᵀ).
  This form is Hermiticity-preserving and trace-preserving by construction. Use sparse storage (CSR) for n ≳ 20. pylcp switches to sparse automatically; RydIQule builds batched dense arrays with `einsum`.
- **Real parametrisation.** Use n real diagonal entries plus Re/Im of the upper triangle, which is n² reals. The ODE becomes real (half the flops of complex n²·2) and Hermiticity holds exactly.
  - pylcp: `transform_into_re_im=True`. It forms U⁻¹LU, asserts the result is real, and raises otherwise. That assertion is itself a Hermiticity test.
  - RydIQule: `make_real`, `get_basis_transform`.
  - CeNTREX-TlF 0.3: `PackedHermitianLayout`, laid out as [diag | re₀₁, im₀₁, …].
  - Gell-Mann / generalised Bloch-vector bases are an equivalent choice.
- **Trace constraint and steady state.** Replace one population equation with Σρ_ii = 1, or eliminate ρ₀₀ (RydIQule `remove_ground`), and solve the resulting non-singular linear system. RydIQule batches `np.linalg.solve` over "stacks" of detunings and Doppler classes. Iterative alternatives: QuantumOptics.jl `steadystate.iterative` (a matrix-free LinearMap), and QuTiP `steadystate(method="direct"|"eigen"|"iterative-gmres"|"power")`.
- **Propagators for piecewise-constant H.** P = expm(L·Δt) is exact and CPTP for a valid Lindbladian. Apply it repeatedly on a uniform grid, which costs one expm.
  - Use `scipy.sparse.linalg.expm_multiply` (Al-Mohy–Higham) for sparse L and time grids. Krylov `expv` is the other option; QuTiP 5 has Krylov; dynamiqs has an `Expm` method; qiskit-dynamics has `jax_expm`/`jax_expm_parallel`.
  - **LASED currently uses `la.eig(A)` plus `la.inv(V)`.** L is non-normal and can be defective or ill-conditioned near degeneracies or critical damping. Replace this with `scipy.linalg.expm` / `expm_multiply`, or at least `la.solve(V, rho0)`.
- **Adaptive RK** (Tsit5, DOPRI5, DOP853) for time-dependent pulses. Positivity is preserved only up to tolerance.
- **Positivity-preserving schemes:** Rouchon/Kraus-map integrators (dynamiqs `Rouchon1/2`). Strang splitting e^{L_D dt/2} e^{L_H dt} e^{L_D dt/2} keeps each factor CPTP.
- **Stiff problems** (large fields and small Γ, or many hyperfine levels): BDF/ESDIRK, for example diffsol, or SciPy `BDF`/`LSODA` with the sparse Jacobian = L.

---

## 3. Simulator survey (what each does well, and what to borrow)

| Tool | Language/backend | Polarisation | Speed/robustness tricks | Lesson for LASED |
|---|---|---|---|---|
| [QuTiP 5](https://github.com/qutip/qutip) ([paper](https://arxiv.org/abs/2412.04705)) | Python/Cython; data layer (Dense/CSR/Dia); [qutip-jax](https://github.com/qutip/qutip-jax) | User builds H | `mesolve` (adams/bdf/vern7/9/dop853), `liouvillian`, `steadystate`, `propagator`, Krylov | **Reference for tests**; the `H + H.dag()` idiom |
| [dynamiqs](https://github.com/dynamiqs/dynamiqs) | JAX + Diffrax, CPU/GPU | User builds H | Batching over H/ρ₀/jump ops, autodiff, Tsit5/Dopri/Expm/Rouchon | Doppler classes as a batch axis; gradients for fitting |
| ADM (Rochester) | Mathematica | Ellipticity/angle vectors, spherical basis | Symbolic + numeric, steady state, Doppler | Gold-standard formalism (Auzinsh–Budker–Rochester) |
| [pylcp](https://github.com/JQIamo/pylcp) (NIST/JQI) | NumPy/SciPy `solve_ivp` | Per-beam Cartesian/spherical; separate d_q and d_q\* operators | Real Re/Im basis, auto-sparse, rate-equation reduction | Closest analogue; copy its d_q/d_q\* split |
| [RydIQule](https://github.com/QTC-UMD/rydiqule) (ARL), [JOSS v2](https://www.theoj.org/joss-papers/joss.08539/10.21105.joss.08539.pdf) | NumPy; networkx graph; scipy/CyRK time solvers | Per-coupling `phase`; explicit conj | `make_real`, `remove_ground`, batched stack solves, exact Doppler | Graph API of levels and couplings; parameter "stacks" |
| [CeNTREX-TlF](https://github.com/ograsdijk/CeNTREX-TlF) 0.3 | Python + **Rust (PyO3)** `centrex_tlf_rust`, optional CUDA (CuPy) | Symbolic polarisation, conj h.c. | Packed Hermitian real layout, native dopri5/tsit5, `solve_lindblad_batch`, AVX2/AVX-512 kernels | **Direct precedent for a Rust LASED backend** |
| [ARC](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator) | Python/C | `getRabiFrequency(..., q)`, q=−1,0,1 ↔ σ⁻,π,σ⁺ | — | Source of dipole matrix elements and lifetimes |
| [ElecSus](https://github.com/durham-qlm/ElecSus) | Python | Arbitrary B-field geometry and polarisation | Weak-probe susceptibility only | Validates LASED's weak-field limit |
| [pairinteraction](https://github.com/pairinteraction/pairinteraction) | C++ backend + Python | — | Sparse diagonalisation | Pattern: C++/Rust core, thin Python API |
| [AtomECS](https://github.com/TeamAtomECS/AtomECS) | **Rust**, ECS (specs→bevy), rayon | Beam polarisation in scattering rates | Rate equations, millions of atoms | Rust architecture; no coherences |
| [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) | Julia | User-built | Lazy operators, `master_dynamic`, matrix-free steady state | Osite et al. alkali OBE toolkit (CPC 314, 109678 (2025)) is built on it |
| `OpticalBlochEquations.jl` | not found; closest: [CoOMBE](https://arxiv.org/abs/2406.19144) (Fortran OBE/MBE), [OBE-Solver](https://github.com/QiSun97/OBE-Solver) (Py+Julia) | — | — | — |
| [Bloqade](https://github.com/QuEraComputing/bloqade-analog) / [Pulser](https://github.com/pasqal-io/Pulser) | Julia→Python / Python (QutipEmulator) | Effective two-level drives (amplitude, phase, detuning) | Builder APIs; Pulser separates sequence from backend | API ergonomics; not multilevel OBE |
| [PennyLane pulse](https://docs.pennylane.ai/en/stable/code/qml_pulse.html), [qiskit-dynamics](https://github.com/qiskit-community/qiskit-dynamics) | JAX | Complex drive amplitude with phase | `LindbladModel(vectorized=True)`, rotating frames, `jax_expm_parallel` | Rotating-frame and evaluation-mode abstractions |
| [MaxwellBloch](https://github.com/tommyogden/maxwellbloch) | Python/QuTiP | Hyperfine structure helper | Propagation (Maxwell–Bloch) | Future propagation extension |

---

## 4. Rust and Python tooling

- **Bindings:** [PyO3](https://pyo3.rs) plus [maturin](https://github.com/PyO3/maturin). Use `pyproject.toml` with `build-backend = "maturin"` and the mixed Python/Rust layout `python-source`. [rust-numpy](https://github.com/PyO3/rust-numpy) (`numpy` crate) gives zero-copy `PyReadonlyArray2<Complex64>` ↔ `ndarray`.
- **ABI:** `pyo3 = { features = ["extension-module", "abi3-py39"] }` produces one wheel per OS/arch.
- **CI:** [maturin-action](https://github.com/PyO3/maturin-action). `maturin generate-ci github` writes the manylinux/musllinux/macOS (universal2)/Windows matrix plus PyPI trusted publishing. Keep a pure-Python fallback so `pip install` works on platforms without wheels.
- **Linear algebra:**
  - `ndarray` is NumPy-like and interoperates with rust-numpy.
  - `nalgebra` is used by `ode_solvers`.
  - [faer](https://docs.rs/faer) offers fast dense LU/QR/eigen for complex (`c64`) and real matrices, sparse LU/Cholesky, and rayon-parallel kernels.
  - None of them ships a production `expm`. Implement Higham's scaling-and-squaring Padé-13, roughly 100 lines on faer matmul and LU solve, or use RK.
- **ODE crates:**
  - [diffsol](https://docs.rs/diffsol) offers BDF, ESDIRK34, TR-BDF2 and Tsit45, sparse KLU, forward sensitivities, and Python bindings [pydiffsol](https://pypi.org/project/pydiffsol/).
  - [ode_solvers](https://docs.rs/ode_solvers) offers Dopri5/Dop853/RK4 for f64 only, which suits the real-packed ρ.
  - [russell_ode](https://docs.rs/russell_ode) offers DoPri5/8 and Radau5, but needs external BLAS/SuiteSparse.
- **Parallel sweeps:** `rayon::par_iter` over Doppler velocity classes, beam-profile radii and detuning scans. Each job owns its own ρ, and the shared L₀ is reused. Release the GIL with `py.allow_threads`.
- **Alternatives:**
  - **numba**: an `@njit` sparse mat-vec RHS with SciPy `solve_ivp`, or numbalsoda. Least effort, but no wheels to build; 3–10× gain on the RHS.
  - **JAX**: GPU, `vmap` over velocity classes, and autodiff for fitting intensities and lifetimes to data. The cost is a heavy dependency and compile latency. Could reuse dynamiqs directly.
  - **C++ with pybind11/nanobind**: same performance as Rust, plus Eigen's `MatrixExponential`. Packaging (scikit-build-core, cibuildwheel) and memory safety are worse than with maturin.
- **Expected sizes:** realistic LASED systems have n ≲ 30, for example the ⁸⁷Rb D2 F=1,2 → F'=0..3 manifold has n=24, so n²=576. A dense expm of 576×576 takes tens of ms in NumPy/LAPACK. The real speedups therefore come from (a) building L vectorised rather than with Python quadruple loops and sympy Wigner calls (the dominant cost now: about 2 min for n=8 in this test), (b) one expm per Δt reused across the time grid, and (c) parallel Doppler and beam averaging.

---

## 5. Recommendations for LASED

**1. Fix pattern.** Replace the hand-expanded equations with a generator built from operators:
```python
# couplings: real reduced CG factors c[e,g,q]; eps: complex spherical amplitudes eps_q
Om = rabi * np.einsum('q,egq->eg', eps, C)           # H_eg block
H = np.diag(energies_rot_frame).astype(complex)
H[np.ix_(E_idx, G_idx)] += Om
H[np.ix_(G_idx, E_idx)] += Om.conj().T               # Hermitian by construction
L = -1j*(np.kron(H, I) - np.kron(I, H.T))            # row-major vec, matches index()
for iq, q in enumerate((-1, 0, 1)):                   # decay incl. coherence transfer
    Lq = np.zeros((n, n)); Lq[np.ix_(G_idx, E_idx)] = np.sqrt(Gamma/N) * C[:, :, iq].T  # N = sum_{g,q} c^2 per e
    LdL = Lq.conj().T @ Lq
    L += np.kron(Lq, Lq.conj()) - 0.5*np.kron(LdL, I) - 0.5*np.kron(I, LdL.T)
# plus 1/tau_f, 1/tau_b, Doppler terms on the diagonal of H / as extra decay
```
For a minimal patch release, the 4-site `np.conj(rabi_factors[i])` change from section 0 is enough. Also:
- Add `polarisation.py` helpers: `from_jones`, `from_cartesian(ε)`, `elliptical(psi, chi, k_dir)` → ε_q, normalised.
- Replace `np.complex`, which was removed in NumPy 1.24.
- Replace `eig`/`inv` with `expm`/`expm_multiply`.
- Cache the Wigner 3j/6j symbols, for example `functools.lru_cache` on float arguments, or `sympy.physics.wigner` once per tuple.

**2. Validation tests** (pytest, run for real, x-linear, rotated-linear, elliptical, and σ± ε):
- **Generator checks:** `P @ L.conj() @ P == L` (Hermiticity preserving); `ones_pop @ L == 0` (trace preserving); `max(eig(L).real) <= 1e-12`.
- **Trajectory checks:** ‖ρ−ρ†‖ < 1e−12; |Tr ρ − 1| < 1e−12; min eig ρ > −1e−9.
- **QuTiP comparison** (`mesolve` with Hermitian H and the c_ops above): J=0→J=1 and F=1→F'=2 with ψ=30°, χ=20°, agreement ≤ 1e−6. The current code fails this at the 10² level (Bug A), and the conj-only fix still fails at 0.1 for F=1→2 (Bug B).
- **Invariance checks:**
  - Rotating ψ about the quantisation axis must leave all populations unchanged and multiply the Zeeman coherences ρ_{m m'} by e^{−i(m−m')ψ} (equivalently, compare with `rotateRho_t`).
  - A global phase on ε must change nothing.
  - χ → −χ must swap σ⁺↔σ⁻ populations under m → −m.
- **Limits:** steady state for J=0→1 σ⁺ against the analytic two-level result ρ_ee = s/2/(1+s+(2Δ/Γ)²); weak-probe absorption against ElecSus.

**3. Architecture:**
- **Python front end** (keep the `LaserAtomSystem` API):
  - Builds the states and computes c[e,g,q] with cached Wigner symbols.
  - Applies the polarisation helpers to get ε_q.
  - Emits a **coupling list** (i, j, complex Ω_ij) for i<j only, the diagonal energies, and the **jump operators** as sparse triplets.
  - The Hermitian half is generated by the backend, never by the user. Alternatively it emits the sparse real-packed Liouvillian (CSR, n²×n²) directly.
- **Backend interface:**
  - `solve(L_csr, rho0_packed, t_grid, method="expm"|"dopri5"|"bdf")`
  - `steady_state(L)`
  - `batch_solve(L0, diag_shifts[k], weights[k])` for Doppler and beam averaging, where only the diagonal detuning changes between classes. This maps directly onto rayon.
- **Implementation order:**
  1. NumPy/SciPy reference: dense expm, `expm_multiply`, `solve_ivp`.
  2. Optional Rust extension (PyO3/maturin, abi3, faer + ode_solvers/diffsol, rayon), following CeNTREX-TlF 0.3's packed-Hermitian layout and batch API.
  3. Optionally a JAX path through dynamiqs for GPU and gradients.
- **Keep the structural assertions** (Hermiticity-preserving L, real transform residual ≈ 0, pylcp style) as runtime debug checks in every backend.
