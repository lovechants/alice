use alloc::vec::Vec;
use core::marker::PhantomData;
use crate::core::ops::{Group, Magma, Monoid, Semigroup, TopologicalGroup};
use crate::groups::matrix_group::{GroupError, MatrixGroup};
use crate::topology::manifold::{Atlas, Dim, Manifold};
use crate::core::scalar::FiniteF64;
use crate::maps::exp_log::HasExpMap;

const DET_TOLERANCE: f64 = 1e-10;
const ORTHO_TOLERANCE: f64 = 1e-8;
const SKEW_TOLERANCE: f64 = 1e-10;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct So<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SoAlgebra<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

impl<D: Dim> So<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let inner = MatrixGroup::new(data, n)?;
        if libm::fabs(inner.det() - 1.0) > DET_TOLERANCE {
            return Err(GroupError::ConstraintViolated("det must equal 1"));
        }
        let t = inner.transpose();
        let rrt = inner.mul(&t).expect("multiplication must succeed");
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                if libm::fabs(rrt.get(i, j) - expected) > ORTHO_TOLERANCE {
                    return Err(GroupError::ConstraintViolated("R*R^T must equal identity"));
                }
            }
        }
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn identity(dim: D) -> Self {
        Self { inner: MatrixGroup::identity(dim.value()), _dim: PhantomData }
    }

    pub fn n(&self) -> usize { self.inner.n() }
    pub fn data(&self) -> Vec<f64> { self.inner.data() }

    pub fn transpose(&self) -> Self {
        Self { inner: self.inner.transpose(), _dim: PhantomData }
    }
}

impl<D: Dim> SoAlgebra<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let inner = MatrixGroup::new(data, n)?;
        for i in 0..n {
            for j in 0..n {
                if libm::fabs(inner.get(i, j) + inner.get(j, i)) > SKEW_TOLERANCE {
                    return Err(GroupError::ConstraintViolated("matrix must be skew-symmetric"));
                }
            }
        }
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn zero(dim: D) -> Self {
        Self {
            inner: MatrixGroup::new(alloc::vec![0.0; dim.value() * dim.value()], dim.value()).unwrap(),
            _dim: PhantomData,
        }
    }

    pub fn data(&self) -> Vec<f64> { self.inner.data() }
    pub fn n(&self) -> usize { self.inner.n() }

    pub fn bracket(&self, other: &Self) -> Self {
        let n = self.inner.n();
        let mut result = alloc::vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut val = 0.0f64;
                for k in 0..n {
                    val += self.inner.get(i, k) * other.inner.get(k, j)
                         - other.inner.get(i, k) * self.inner.get(k, j);
                }
                result[i * n + j] = val;
            }
        }
        Self {
            inner: MatrixGroup::new(result, n).unwrap(),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Magma for So<D> {
    fn op(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.mul(&other.inner).expect("SO multiplication must succeed"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Semigroup for So<D> {}

impl<D: Dim> Monoid for So<D> {
    fn identity() -> Self {
        panic!("SO identity requires dimension; use So::identity(dim)")
    }
}

impl<D: Dim> Group for So<D> {
    fn inverse(&self) -> Self {
        Self { inner: self.inner.transpose(), _dim: PhantomData }
    }
}

impl<D: Dim> TopologicalGroup for So<D> {}

impl<D: Dim> Manifold for So<D> {
    type Scalar = FiniteF64;
    fn dim(&self) -> usize { let n = self.inner.n(); n * (n - 1) / 2 }
    fn atlas(&self) -> &Atlas<FiniteF64> { unimplemented!("SO atlas not yet constructed") }
}

impl<D: Dim> HasExpMap for So<D> {
    type Algebra = SoAlgebra<D>;
    fn exp(x: &SoAlgebra<D>) -> Self {
        Self { inner: x.inner.exp(), _dim: PhantomData }
    }
    fn log(&self) -> Option<SoAlgebra<D>> {
        Some(SoAlgebra { inner: self.inner.log()?, _dim: PhantomData })
    }
}
