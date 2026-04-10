use crate::maps::faer_bridge;

pub trait HasExpMap: Sized {
    type Algebra: Clone + PartialEq;
    fn exp(x: &Self::Algebra) -> Self;
    fn log(&self) -> Option<Self::Algebra>;
}

pub trait MatrixExpMap: HasExpMap {
    fn matrix_data(x: &Self::Algebra) -> (alloc::vec::Vec<f64>, usize);
    fn from_matrix_data(data: alloc::vec::Vec<f64>, n: usize) -> Self;
    fn algebra_from_matrix_data(data: alloc::vec::Vec<f64>, n: usize) -> Self::Algebra;
    fn to_matrix_data(&self) -> (alloc::vec::Vec<f64>, usize);
}

impl<T: MatrixExpMap> HasExpMap for T {
    type Algebra = T::Algebra;
    fn exp(x: &Self::Algebra) -> Self {
        let (data, n) = T::matrix_data(x);
        let result = faer_bridge::matrix_exp(&data, n);
        T::from_matrix_data(result, n)
    }
    fn log(&self) -> Option<Self::Algebra> {
        let (data, n) = self.to_matrix_data();
        let result = faer_bridge::matrix_log(&data, n)?;
        Some(T::algebra_from_matrix_data(result, n))
    }
}

pub fn exp_log_roundtrip_error(
    data: &[f64],
    n: usize,
) -> Option<f64> {
    let exp_data = faer_bridge::matrix_exp(data, n);
    let log_data = faer_bridge::matrix_log(&exp_data, n)?;
    let max_err = data.iter()
        .zip(log_data.iter())
        .map(|(a, b)| libm::fabs(a - b))
        .fold(0.0f64, f64::max);
    Some(max_err)
}

pub fn log_exp_roundtrip_error(
    data: &[f64],
    n: usize,
) -> Option<f64> {
    let log_data = faer_bridge::matrix_log(data, n)?;
    let exp_data = faer_bridge::matrix_exp(&log_data, n);
    let max_err = data.iter()
        .zip(exp_data.iter())
        .map(|(a, b)| libm::fabs(a - b))
        .fold(0.0f64, f64::max);
    Some(max_err)
}
