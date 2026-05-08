use alloc::vec::Vec;
use faer::Mat;
use faer::prelude::Solve;
use crate::maps::faer_bridge;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixGroup {
    data: Vec<OrderedF64>,
    n: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrderedF64(pub f64);

impl Eq for OrderedF64 {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupError {
    DimensionMismatch,
    NotInvertible,
    ConstraintViolated(&'static str),
}

impl core::fmt::Display for GroupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GroupError::DimensionMismatch => write!(f, "matrix dimension mismatch"),
            GroupError::NotInvertible => write!(f, "matrix is not invertible"),
            GroupError::ConstraintViolated(msg) => write!(f, "group constraint violated: {}", msg),
        }
    }
}

impl MatrixGroup {
    pub fn new(data: Vec<f64>, n: usize) -> Result<Self, GroupError> {
        if data.len() != n * n {
            return Err(GroupError::DimensionMismatch);
        }
        Ok(Self {
            data: data.into_iter().map(OrderedF64).collect(),
            n,
        })
    }

    pub fn identity(n: usize) -> Self {
        let data = (0..n * n)
            .map(|k| OrderedF64(if k / n == k % n { 1.0 } else { 0.0 }))
            .collect();
        Self { data, n }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn data(&self) -> Vec<f64> {
        self.data.iter().map(|x| x.0).collect()
    }

    pub fn data_ref(&self) -> &[OrderedF64] {
        &self.data
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j].0
    }

    pub fn mul(&self, other: &Self) -> Result<Self, GroupError> {
        if self.n != other.n {
            return Err(GroupError::DimensionMismatch);
        }
        let a = to_faer_inner(&self.data(), self.n);
        let b = to_faer_inner(&other.data(), self.n);
        let c = a * b;
        Ok(Self {
            data: from_faer_inner(&c).into_iter().map(OrderedF64).collect(),
            n: self.n,
        })
    }

    pub fn inverse(&self) -> Result<Self, GroupError> {
        let a = to_faer_inner(&self.data(), self.n);
        let det = a.determinant();
        if libm::fabs(det) < 1e-14 {
            return Err(GroupError::NotInvertible);
        }
        let lu = a.partial_piv_lu();
        let identity = Mat::<f64>::identity(self.n, self.n);
        let inv = lu.solve(identity);
        Ok(Self {
            data: from_faer_inner(&inv).into_iter().map(OrderedF64).collect(),
            n: self.n,
        })
    }

    pub fn det(&self) -> f64 {
        let a = to_faer_inner(&self.data(), self.n);
        a.determinant()
    }

    pub fn transpose(&self) -> Self {
        let mut data = alloc::vec![OrderedF64(0.0); self.n * self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                data[j * self.n + i] = OrderedF64(self.data[i * self.n + j].0);
            }
        }
        Self { data, n: self.n }
    }

    pub fn exp(&self) -> Self {
        let result = faer_bridge::matrix_exp(&self.data(), self.n);
        Self {
            data: result.into_iter().map(OrderedF64).collect(),
            n: self.n,
        }
    }

    pub fn log(&self) -> Option<Self> {
        let result = faer_bridge::matrix_log(&self.data(), self.n)?;
        Some(Self {
            data: result.into_iter().map(OrderedF64).collect(),
            n: self.n,
        })
    }
}

pub(crate) fn to_faer_inner(data: &[f64], n: usize) -> Mat<f64> {
    Mat::from_fn(n, n, |i, j| data[i * n + j])
}

pub(crate) fn from_faer_inner(m: &Mat<f64>) -> Vec<f64> {
    let n = m.nrows();
    let mut out = alloc::vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            out[i * n + j] = m[(i, j)];
        }
    }
    out
}
