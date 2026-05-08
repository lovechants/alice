use alloc::vec::Vec;
use core::marker::PhantomData;
use crate::core::ops::{Group, Magma, Monoid, Semigroup, TopologicalGroup};
use crate::groups::matrix_group::{GroupError, MatrixGroup};
use crate::topology::manifold::{Atlas, Dim, Manifold};
use crate::core::scalar::FiniteF64;
use crate::maps::exp_log::HasExpMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gl<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlAlgebra<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

impl<D: Dim> Gl<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let inner = MatrixGroup::new(data, n)?;
        if libm::fabs(inner.det()) < 1e-14 {
            return Err(GroupError::ConstraintViolated("matrix must be invertible"));
        }
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn identity(dim: D) -> Self {
        Self { inner: MatrixGroup::identity(dim.value()), _dim: PhantomData }
    }

    pub fn n(&self) -> usize { self.inner.n() }
    pub fn data(&self) -> Vec<f64> { self.inner.data() }
    pub fn det(&self) -> f64 { self.inner.det() }
}

impl<D: Dim> GlAlgebra<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let inner = MatrixGroup::new(data, dim.value())?;
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
}

impl<D: Dim> Magma for Gl<D> {
    fn op(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.mul(&other.inner).expect("GL multiplication must succeed"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Semigroup for Gl<D> {}

impl<D: Dim> Monoid for Gl<D> {
    fn identity() -> Self {
        panic!("GL identity requires dimension; use Gl::identity(dim)")
    }
}

impl<D: Dim> Group for Gl<D> {
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse().expect("GL element must be invertible"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> TopologicalGroup for Gl<D> {}

impl<D: Dim> Manifold for Gl<D> {
    type Scalar = FiniteF64;
    fn dim(&self) -> usize { self.inner.n() * self.inner.n() }
    fn atlas(&self) -> &Atlas<FiniteF64> { unimplemented!("GL atlas not yet constructed") }
}

impl<D: Dim> HasExpMap for Gl<D> {
    type Algebra = GlAlgebra<D>;
    fn exp(x: &GlAlgebra<D>) -> Self {
        Self { inner: x.inner.exp(), _dim: PhantomData }
    }
    fn log(&self) -> Option<GlAlgebra<D>> {
        Some(GlAlgebra { inner: self.inner.log()?, _dim: PhantomData })
    }
}
