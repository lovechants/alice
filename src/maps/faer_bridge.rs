use alloc::vec::Vec;
use faer::Mat;
use faer::prelude::Solve;
use num_complex::Complex;

type C64 = Complex<f64>;

pub fn matrix_exp(data: &[f64], n: usize) -> Vec<f64> {
    let a = to_faer(data, n);
    let result = pade_exp(&a, n);
    from_faer(&result)
}

pub fn matrix_log(data: &[f64], n: usize) -> Option<Vec<f64>> {
    let a = to_faer(data, n);
    let result = eigen_log(&a, n)?;
    let out = from_faer(&result);
    if out.iter().any(|v| v.is_nan() || v.is_infinite()) {
        None
    } else {
        Some(out)
    }
}

fn to_faer(data: &[f64], n: usize) -> Mat<f64> {
    assert_eq!(data.len(), n * n);
    Mat::from_fn(n, n, |i, j| data[i * n + j])
}

fn from_faer(m: &Mat<f64>) -> Vec<f64> {
    let n = m.nrows();
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            out.push(m[(i, j)]);
        }
    }
    out
}

fn pade_exp(a: &Mat<f64>, n: usize) -> Mat<f64> {
    let norm = mat_inf_norm(a);
    let mut s = 0i32;
    let mut scaled = a.clone();
    if norm > 0.5 {
        s = (libm::ceil(libm::log2(norm / 0.5)) as i32).max(1);
        let factor = libm::pow(2.0, -(s as f64));
        scaled = faer::Scale(factor) * a;
    }
    let p = horner(&scaled, &PADE_NUM, n);
    let q = horner(&scaled, &PADE_DEN, n);
    let lu = q.partial_piv_lu();
    let mut r: Mat<f64> = lu.solve(p);
    for _ in 0..s {
        r = &r * &r;
    }
    r
}

fn eigen_log(a: &Mat<f64>, n: usize) -> Option<Mat<f64>> {
    let evd = a.eigen().ok()?;
    let u = evd.U();
    let s_col = evd.S();

    let mut log_diag: Vec<C64> = Vec::with_capacity(n);
    for i in 0..n {
        let lambda: C64 = s_col.column_vector()[i];
        let r = libm::sqrt(lambda.re * lambda.re + lambda.im * lambda.im);
        if r <= 0.0 {
            return None;
        }
        let theta = libm::atan2(lambda.im, lambda.re);
        log_diag.push(Complex::new(libm::log(r), theta));
    }

    let log_d: Mat<C64> = Mat::from_fn(n, n, |i, j| {
        if i == j {
            log_diag[i]
        } else {
            Complex::new(0.0, 0.0)
        }
    });

    let u_owned: Mat<C64> = u.to_owned();
    let u_lu = u_owned.partial_piv_lu();
    let identity_c: Mat<C64> = Mat::identity(n, n);
    let u_inv: Mat<C64> = u_lu.solve(identity_c);
    let result_c: Mat<C64> = &u_owned * &log_d * &u_inv;

    Some(Mat::from_fn(n, n, |i, j| result_c[(i, j)].re))
}

fn horner(a: &Mat<f64>, coeffs: &[f64], n: usize) -> Mat<f64> {
    let m = coeffs.len();
    let mut result: Mat<f64> = faer::Scale(coeffs[m - 1]) * Mat::<f64>::identity(n, n);
    for k in (0..m - 1).rev() {
        result = &result * a + faer::Scale(coeffs[k]) * Mat::<f64>::identity(n, n);
    }
    result
}

fn mat_inf_norm(a: &Mat<f64>) -> f64 {
    let n = a.nrows();
    let mut max_row = 0.0f64;
    for i in 0..n {
        let row_sum: f64 = (0..n).map(|j| libm::fabs(a[(i, j)])).sum();
        if row_sum > max_row {
            max_row = row_sum;
        }
    }
    max_row
}

const PADE_NUM: [f64; 14] = [
    1.0,
    0.5,
    0.12,
    1.833333333333333e-2,
    1.992063492063492e-3,
    1.630434782608696e-4,
    1.035196687370600e-5,
    5.175983436853003e-7,
    2.043151389366194e-8,
    6.306659613335511e-10,
    1.483524052786891e-11,
    2.529153491597966e-13,
    2.810170546428327e-15,
    1.544049750670308e-17,
];

const PADE_DEN: [f64; 14] = [
    1.0,
    -0.5,
    0.12,
    -1.833333333333333e-2,
    1.992063492063492e-3,
    -1.630434782608696e-4,
    1.035196687370600e-5,
    -5.175983436853003e-7,
    2.043151389366194e-8,
    -6.306659613335511e-10,
    1.483524052786891e-11,
    -2.529153491597966e-13,
    2.810170546428327e-15,
    -1.544049750670308e-17,
];
