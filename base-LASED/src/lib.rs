//! Rust backend for LASED.
//!
//! Integrates the vectorised Lindblad master equation
//!
//!     d vec(rho)/dt = [ L0 + sum_j c_{m,j}(t) M_j ] vec(rho)
//!
//! for an ensemble of members m (e.g. Doppler classes, beam-profile rings),
//! where every member shares the sparse superoperators L0 and M_j but has its
//! own scalar coefficients
//!
//!     c_{m,j}(t) = alpha[m, j] * u_j(t) + beta[m, j]
//!
//! and u_j(t) are control waveforms (pulse envelopes, chirps, ...) given as
//! samples that are linearly interpolated. Each member is integrated with an
//! adaptive Dormand-Prince 5(4) scheme; members run in parallel with rayon and
//! their weighted results are summed.

use num_complex::Complex64 as C64;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Compressed sparse row matrix (square, n x n).
struct Csr {
    data: Vec<C64>,
    indices: Vec<usize>,
    indptr: Vec<usize>,
}

impl Csr {
    fn from_py(
        n: usize,
        data: PyReadonlyArray1<C64>,
        indices: PyReadonlyArray1<i64>,
        indptr: PyReadonlyArray1<i64>,
    ) -> PyResult<Self> {
        let data = data.as_slice()?.to_vec();
        let indices: Vec<usize> = indices.as_slice()?.iter().map(|&i| i as usize).collect();
        let indptr: Vec<usize> = indptr.as_slice()?.iter().map(|&i| i as usize).collect();
        if indptr.len() != n + 1 {
            return Err(PyValueError::new_err("CSR indptr has the wrong length"));
        }
        if data.len() != indices.len() || indptr[n] != data.len() {
            return Err(PyValueError::new_err("CSR data/indices are inconsistent"));
        }
        if indices.iter().any(|&i| i >= n) {
            return Err(PyValueError::new_err("CSR column index out of range"));
        }
        Ok(Csr { data, indices, indptr })
    }

    /// sum_k c_k A_k for matrices sharing the same shape (sparse accumulator per row).
    fn linear_combination(n: usize, terms: &[(&Csr, C64)]) -> Csr {
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = Vec::with_capacity(n + 1);
        indptr.push(0);
        let mut acc = vec![C64::new(0.0, 0.0); n];
        let mut used = vec![false; n];
        let mut cols: Vec<usize> = Vec::new();
        for row in 0..n {
            for (a, c) in terms {
                for k in a.indptr[row]..a.indptr[row + 1] {
                    let col = a.indices[k];
                    if !used[col] {
                        used[col] = true;
                        cols.push(col);
                    }
                    acc[col] += a.data[k] * c;
                }
            }
            cols.sort_unstable();
            for &col in &cols {
                if acc[col] != C64::new(0.0, 0.0) {
                    data.push(acc[col]);
                    indices.push(col);
                }
                acc[col] = C64::new(0.0, 0.0);
                used[col] = false;
            }
            cols.clear();
            indptr.push(data.len());
        }
        Csr { data, indices, indptr }
    }

    /// y += coef * A x
    #[inline]
    fn matvec_add(&self, coef: C64, x: &[C64], y: &mut [C64]) {
        for (row, yr) in y.iter_mut().enumerate() {
            let mut acc = C64::new(0.0, 0.0);
            for k in self.indptr[row]..self.indptr[row + 1] {
                acc += self.data[k] * x[self.indices[k]];
            }
            *yr += coef * acc;
        }
    }
}

/// Piecewise-linear control waveforms u_j(t), held constant outside the grid.
struct Controls {
    t: Vec<f64>,
    u: Vec<Vec<f64>>, // [op][sample]
}

impl Controls {
    fn value(&self, j: usize, t: f64) -> f64 {
        let ts = &self.t;
        if ts.is_empty() {
            return 1.0;
        }
        let u = &self.u[j];
        if t <= ts[0] {
            return u[0];
        }
        let last = ts.len() - 1;
        if t >= ts[last] {
            return u[last];
        }
        // Binary search for the interval containing t
        let k = ts.partition_point(|&s| s <= t) - 1;
        let w = (t - ts[k]) / (ts[k + 1] - ts[k]);
        u[k] * (1.0 - w) + u[k + 1] * w
    }
}

struct Problem {
    n: usize,
    l0: Csr,
    ops: Vec<Csr>,
    controls: Controls,
}

/// The generator of one ensemble member: a merged static matrix plus time-dependent terms.
struct MemberGenerator<'a> {
    stat: Csr,
    dynamic: Vec<(usize, f64, f64)>, // (op index, alpha, beta)
    problem: &'a Problem,
}

impl<'a> MemberGenerator<'a> {
    fn new(problem: &'a Problem, alpha: &[f64], beta: &[f64]) -> Self {
        let mut terms: Vec<(&Csr, C64)> = vec![(&problem.l0, C64::new(1.0, 0.0))];
        let mut dynamic = Vec::new();
        for (j, op) in problem.ops.iter().enumerate() {
            if alpha[j] != 0.0 {
                dynamic.push((j, alpha[j], beta[j]));
            } else if beta[j] != 0.0 {
                terms.push((op, C64::new(beta[j], 0.0)));
            }
        }
        let stat = Csr::linear_combination(problem.n, &terms);
        MemberGenerator { stat, dynamic, problem }
    }

    fn rhs(&self, t: f64, x: &[C64], dx: &mut [C64]) {
        dx.iter_mut().for_each(|v| *v = C64::new(0.0, 0.0));
        self.stat.matvec_add(C64::new(1.0, 0.0), x, dx);
        for &(j, a, b) in &self.dynamic {
            let c = a * self.problem.controls.value(j, t) + b;
            if c != 0.0 {
                self.problem.ops[j].matvec_add(C64::new(c, 0.0), x, dx);
            }
        }
    }
}

struct Tolerances {
    rtol: f64,
    atol: f64,
    max_step: f64,
}

// Dormand-Prince 5(4) Butcher tableau
const C2: f64 = 1.0 / 5.0;
const C3: f64 = 3.0 / 10.0;
const C4: f64 = 4.0 / 5.0;
const C5: f64 = 8.0 / 9.0;
const A21: f64 = 1.0 / 5.0;
const A31: f64 = 3.0 / 40.0;
const A32: f64 = 9.0 / 40.0;
const A41: f64 = 44.0 / 45.0;
const A42: f64 = -56.0 / 15.0;
const A43: f64 = 32.0 / 9.0;
const A51: f64 = 19372.0 / 6561.0;
const A52: f64 = -25360.0 / 2187.0;
const A53: f64 = 64448.0 / 6561.0;
const A54: f64 = -212.0 / 729.0;
const A61: f64 = 9017.0 / 3168.0;
const A62: f64 = -355.0 / 33.0;
const A63: f64 = 46732.0 / 5247.0;
const A64: f64 = 49.0 / 176.0;
const A65: f64 = -5103.0 / 18656.0;
const B1: f64 = 35.0 / 384.0;
const B3: f64 = 500.0 / 1113.0;
const B4: f64 = 125.0 / 192.0;
const B5: f64 = -2187.0 / 6784.0;
const B6: f64 = 11.0 / 84.0;
const E1: f64 = 71.0 / 57600.0;
const E3: f64 = -71.0 / 16695.0;
const E4: f64 = 71.0 / 1920.0;
const E5: f64 = -17253.0 / 339200.0;
const E6: f64 = 22.0 / 525.0;
const E7: f64 = -1.0 / 40.0;

/// Integrate one ensemble member; returns the state at every output time
/// (row-major, n_out x n2).
fn integrate_member(
    p: &MemberGenerator,
    rho0: &[C64],
    t_out: &[f64],
    tol: &Tolerances,
) -> Result<Vec<C64>, String> {
    let n = rho0.len();
    let mut out = Vec::with_capacity(n * t_out.len());
    out.extend_from_slice(rho0);
    if t_out.len() < 2 {
        return Ok(out);
    }

    let mut y = rho0.to_vec();
    let mut t = t_out[0];
    let mut k1 = vec![C64::new(0.0, 0.0); n];
    let mut k2 = k1.clone();
    let mut k3 = k1.clone();
    let mut k4 = k1.clone();
    let mut k5 = k1.clone();
    let mut k6 = k1.clone();
    let mut k7 = k1.clone();
    let mut tmp = k1.clone();
    let mut ynew = k1.clone();

    p.rhs(t, &y, &mut k1);

    // Initial step guess from the size of the derivative
    let norm = |v: &[C64], scale: &[C64]| -> f64 {
        let s: f64 = v
            .iter()
            .zip(scale)
            .map(|(a, b)| {
                let sc = tol.atol + tol.rtol * b.norm();
                (a.norm() / sc).powi(2)
            })
            .sum();
        (s / n as f64).sqrt()
    };
    let d1 = norm(&k1, &y);
    let span = t_out[t_out.len() - 1] - t_out[0];
    let mut h = if d1 > 1e-12 { 0.01 / d1 } else { span * 1e-3 };
    h = h.min(tol.max_step).min(span).max(1e-12 * span.abs().max(1.0));

    let max_steps: usize = 50_000_000;
    let mut steps = 0usize;

    for &t_target in &t_out[1..] {
        while t < t_target {
            steps += 1;
            if steps > max_steps {
                return Err("maximum number of integration steps exceeded".into());
            }
            let h_unclipped = h;
            let mut last = false;
            if t + h >= t_target {
                h = t_target - t;
                last = true;
            }

            macro_rules! stage {
                ($k:ident, $c:expr, [$(($a:expr, $kk:ident)),*]) => {
                    for i in 0..n {
                        tmp[i] = y[i] $(+ $kk[i] * (h * $a))*;
                    }
                    p.rhs(t + $c * h, &tmp, &mut $k);
                };
            }
            stage!(k2, C2, [(A21, k1)]);
            stage!(k3, C3, [(A31, k1), (A32, k2)]);
            stage!(k4, C4, [(A41, k1), (A42, k2), (A43, k3)]);
            stage!(k5, C5, [(A51, k1), (A52, k2), (A53, k3), (A54, k4)]);
            stage!(k6, 1.0, [(A61, k1), (A62, k2), (A63, k3), (A64, k4), (A65, k5)]);
            for i in 0..n {
                ynew[i] = y[i] + (k1[i] * B1 + k3[i] * B3 + k4[i] * B4 + k5[i] * B5 + k6[i] * B6) * h;
            }
            p.rhs(t + h, &ynew, &mut k7);

            // Error estimate (difference between the 5th and 4th order solutions)
            let mut err_sq = 0.0;
            for i in 0..n {
                let e = (k1[i] * E1 + k3[i] * E3 + k4[i] * E4 + k5[i] * E5 + k6[i] * E6 + k7[i] * E7) * h;
                let sc = tol.atol + tol.rtol * y[i].norm().max(ynew[i].norm());
                err_sq += (e.norm() / sc).powi(2);
            }
            let err = (err_sq / n as f64).sqrt();

            if err <= 1.0 || h <= 1e-14 * t.abs().max(1.0) {
                t = if last { t_target } else { t + h };
                std::mem::swap(&mut y, &mut ynew);
                std::mem::swap(&mut k1, &mut k7); // first-same-as-last
                if !y.iter().all(|v| v.re.is_finite() && v.im.is_finite()) {
                    return Err("integration produced non-finite values".into());
                }
                let fac = if err == 0.0 { 5.0 } else { (0.9 * err.powf(-0.2)).clamp(0.2, 5.0) };
                // After a step clipped to an output time, resume from the unclipped size
                h = if last { h_unclipped.max(h * fac) } else { h * fac }.min(tol.max_step);
            } else {
                let fac = (0.9 * err.powf(-0.2)).clamp(0.1, 1.0);
                h *= fac;
            }
        }
        out.extend_from_slice(&y);
    }
    Ok(out)
}

/// Integrate an ensemble of Lindblad master equations sharing sparse superoperators.
///
/// Arguments (n2 = dimension of vec(rho), n_ops = number of time-dependent terms):
///     n2: size of the vectorised density matrix.
///     l0: (data, indices, indptr) CSR arrays of the static superoperator L0.
///     ops: list of CSR triples for the superoperators M_j.
///     ctrl_t: sample times of the control waveforms (may be empty: then u_j = 1).
///     ctrl_u: (n_ops, len(ctrl_t)) control samples u_j(t).
///     alpha, beta: (n_members, n_ops) coefficients, c_{m,j}(t) = alpha*u_j(t) + beta.
///     weights: (n_members,) weights of the members in the returned average.
///     rho0: initial vec(rho) (complex, length n2), given at t_out[0].
///     t_out: increasing output times.
///     rtol, atol, max_step: adaptive step control.
///
/// Returns:
///     complex array of shape (len(t_out), n2): the weighted sum over members.
#[pyfunction]
#[pyo3(signature = (n2, l0, ops, ctrl_t, ctrl_u, alpha, beta, weights, rho0, t_out, rtol=1e-8, atol=1e-10, max_step=f64::INFINITY))]
#[allow(clippy::too_many_arguments)]
fn evolve_ensemble<'py>(
    py: Python<'py>,
    n2: usize,
    l0: (PyReadonlyArray1<'py, C64>, PyReadonlyArray1<'py, i64>, PyReadonlyArray1<'py, i64>),
    ops: Vec<(PyReadonlyArray1<'py, C64>, PyReadonlyArray1<'py, i64>, PyReadonlyArray1<'py, i64>)>,
    ctrl_t: PyReadonlyArray1<'py, f64>,
    ctrl_u: PyReadonlyArray2<'py, f64>,
    alpha: PyReadonlyArray2<'py, f64>,
    beta: PyReadonlyArray2<'py, f64>,
    weights: PyReadonlyArray1<'py, f64>,
    rho0: PyReadonlyArray1<'py, C64>,
    t_out: PyReadonlyArray1<'py, f64>,
    rtol: f64,
    atol: f64,
    max_step: f64,
) -> PyResult<Bound<'py, PyArray2<C64>>> {
    let l0 = Csr::from_py(n2, l0.0, l0.1, l0.2)?;
    let ops = ops
        .into_iter()
        .map(|(d, i, p)| Csr::from_py(n2, d, i, p))
        .collect::<PyResult<Vec<_>>>()?;
    let n_ops = ops.len();

    let ctrl_t = ctrl_t.as_slice()?.to_vec();
    let ctrl_u = ctrl_u.as_array();
    if !ctrl_t.is_empty() && (ctrl_u.nrows() != n_ops || ctrl_u.ncols() != ctrl_t.len()) {
        return Err(PyValueError::new_err("ctrl_u must have shape (len(ops), len(ctrl_t))"));
    }
    if ctrl_t.windows(2).any(|w| w[1] <= w[0]) {
        return Err(PyValueError::new_err("ctrl_t must be strictly increasing"));
    }
    let u: Vec<Vec<f64>> = if ctrl_t.is_empty() {
        vec![Vec::new(); n_ops]
    } else {
        ctrl_u.rows().into_iter().map(|r| r.to_vec()).collect()
    };

    let alpha = alpha.as_array().to_owned();
    let beta = beta.as_array().to_owned();
    let weights = weights.as_slice()?.to_vec();
    let n_members = weights.len();
    if alpha.dim() != (n_members, n_ops) || beta.dim() != (n_members, n_ops) {
        return Err(PyValueError::new_err("alpha and beta must have shape (len(weights), len(ops))"));
    }
    let rho0 = rho0.as_slice()?.to_vec();
    if rho0.len() != n2 {
        return Err(PyValueError::new_err("rho0 must have length n2"));
    }
    let t_out = t_out.as_slice()?.to_vec();
    if t_out.windows(2).any(|w| w[1] < w[0]) {
        return Err(PyValueError::new_err("t_out must be non-decreasing"));
    }
    let tol = Tolerances { rtol, atol, max_step: if max_step > 0.0 { max_step } else { f64::INFINITY } };
    let problem = Problem { n: n2, l0, ops, controls: Controls { t: ctrl_t, u } };
    let n_out = t_out.len();

    let result: Result<Vec<C64>, String> = py.detach(|| {
        (0..n_members)
            .into_par_iter()
            .map(|m| {
                let a: Vec<f64> = alpha.row(m).to_vec();
                let b: Vec<f64> = beta.row(m).to_vec();
                let generator = MemberGenerator::new(&problem, &a, &b);
                let mut traj = integrate_member(&generator, &rho0, &t_out, &tol)?;
                let w = weights[m];
                traj.iter_mut().for_each(|v| *v *= w);
                Ok(traj)
            })
            .try_reduce(
                || vec![C64::new(0.0, 0.0); n_out * n2],
                |mut acc, x| {
                    acc.iter_mut().zip(x).for_each(|(a, b)| *a += b);
                    Ok(acc)
                },
            )
    });
    let flat = result.map_err(PyValueError::new_err)?;
    let arr = numpy::ndarray::Array2::from_shape_vec((n_out, n2), flat)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(arr.into_pyarray(py))
}

/// Number of threads rayon will use for ensemble parallelism.
#[pyfunction]
fn num_threads() -> usize {
    rayon::current_num_threads()
}

#[pymodule]
fn _lased_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(evolve_ensemble, m)?)?;
    m.add_function(wrap_pyfunction!(num_threads, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
