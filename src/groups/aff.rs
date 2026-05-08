use alloc::vec::Vec;
use core::marker::PhantomData;
use crate::core::ops::{Group, Magma, Monoid, Semigroup, TopologicalGroup};
use crate::groups::matrix_group::{GroupError, MatrixGroup};
use crate::topology::manifold::{Atlas, Dim, Manifold};
use crate::core::scalar::FiniteF64;
use crate::maps::exp_log::HasExpMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Aff<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffAlgebra<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

impl<D: Dim> Aff<D> {
    pub fn from_homogeneous(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let nh = n + 1;
        if data.len() != nh * nh {
            return Err(GroupError::DimensionMismatch);
        }
        for j in 0..n {
            if libm::fabs(data[n * nh + j]) > 1e-10 {
                return Err(GroupError::ConstraintViolated("bottom row must be [0..0,1]"));
            }
        }
        if libm::fabs(data[n * nh + n] - 1.0) > 1e-10 {
            return Err(GroupError::ConstraintViolated("bottom right must be 1"));
        }
        let inner = MatrixGroup::new(data, nh)?;
        let linear_det: f64 = {
            let linear: Vec<f64> = (0..n)
                .flat_map(|i| (0..n).map(move |j| (i, j)))
                .map(|(i, j)| inner.get(i, j))
                .collect();
            MatrixGroup::new(linear, n)
                .map(|m| m.det())
                .unwrap_or(0.0)
        };
        if libm::fabs(linear_det) < 1e-10 {
            return Err(GroupError::ConstraintViolated("linear part must be invertible"));
        }
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn from_parts(linear: Vec<f64>, translation: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        if linear.len() != n * n || translation.len() != n {
            return Err(GroupError::DimensionMismatch);
        }
        let nh = n + 1;
        let mut data = alloc::vec![0.0f64; nh * nh];
        for i in 0..n {
            for j in 0..n {
                data[i * nh + j] = linear[i * n + j];
            }
            data[i * nh + n] = translation[i];
        }
        data[n * nh + n] = 1.0;
        Self::from_homogeneous(data, dim)
    }

    pub fn identity(dim: D) -> Self {
        Self { inner: MatrixGroup::identity(dim.value() + 1), _dim: PhantomData }
    }

    pub fn linear_part(&self, dim: D) -> Vec<f64> {
        let n = dim.value();
        (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| self.inner.get(i, j))
            .collect()
    }

    pub fn translation_part(&self, dim: D) -> Vec<f64> {
        let n = dim.value();
        (0..n).map(|i| self.inner.get(i, n)).collect()
    }

    pub fn to_homogeneous(&self) -> Vec<f64> { self.inner.data() }
    pub fn n(&self) -> usize { self.inner.n() - 1 }
}

impl<D: Dim> AffAlgebra<D> {
    pub fn from_homogeneous(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let nh = dim.value() + 1;
        if data.len() != nh * nh {
            return Err(GroupError::DimensionMismatch);
        }
        let inner = MatrixGroup::new(data, nh)?;
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn zero(dim: D) -> Self {
        let nh = dim.value() + 1;
        Self {
            inner: MatrixGroup::new(alloc::vec![0.0; nh * nh], nh).unwrap(),
            _dim: PhantomData,
        }
    }

    pub fn data(&self) -> Vec<f64> { self.inner.data() }
}

impl<D: Dim> Magma for Aff<D> {
    fn op(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.mul(&other.inner).expect("Aff multiplication must succeed"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Semigroup for Aff<D> {}

impl<D: Dim> Monoid for Aff<D> {
    fn identity() -> Self {
        panic!("Aff identity requires dimension; use Aff::identity(dim)")
    }
}

impl<D: Dim> Group for Aff<D> {
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse().expect("Aff element must be invertible"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> TopologicalGroup for Aff<D> {}

impl<D: Dim> Manifold for Aff<D> {
    type Scalar = FiniteF64;
    fn dim(&self) -> usize { let n = self.inner.n() - 1; n * n + n }
    fn atlas(&self) -> &Atlas<FiniteF64> { unimplemented!("Aff atlas not yet constructed") }
}

impl<D: Dim> HasExpMap for Aff<D> {
    type Algebra = AffAlgebra<D>;
    fn exp(x: &AffAlgebra<D>) -> Self {
        Self { inner: x.inner.exp(), _dim: PhantomData }
    }
    fn log(&self) -> Option<AffAlgebra<D>> {
        Some(AffAlgebra { inner: self.inner.log()?, _dim: PhantomData })
    }
}
