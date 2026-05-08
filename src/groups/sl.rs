use alloc::vec::Vec;
use core::marker::PhantomData;
use crate::core::ops::{Group, Magma, Monoid, Semigroup, TopologicalGroup};
use crate::groups::matrix_group::{GroupError, MatrixGroup};
use crate::topology::manifold::{Atlas, Dim, Manifold};
use crate::core::scalar::FiniteF64;
use crate::maps::exp_log::HasExpMap;

const DET_TOLERANCE: f64 = 1e-10;
const TRACE_TOLERANCE: f64 = 1e-10;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sl<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SlAlgebra<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

impl<D: Dim> Sl<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let inner = MatrixGroup::new(data, n)?;
        if libm::fabs(inner.det() - 1.0) > DET_TOLERANCE {
            return Err(GroupError::ConstraintViolated("det must equal 1"));
        }
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn identity(dim: D) -> Self {
        Self { inner: MatrixGroup::identity(dim.value()), _dim: PhantomData }
    }

    pub fn n(&self) -> usize { self.inner.n() }
    pub fn data(&self) -> Vec<f64> { self.inner.data() }
}

impl<D: Dim> SlAlgebra<D> {
    pub fn new(data: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        let inner = MatrixGroup::new(data, n)?;
        let trace: f64 = (0..n).map(|i| inner.get(i, i)).sum();
        if libm::fabs(trace) > TRACE_TOLERANCE {
            return Err(GroupError::ConstraintViolated("trace must equal 0"));
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
}

impl<D: Dim> Magma for Sl<D> {
    fn op(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.mul(&other.inner).expect("SL multiplication must succeed"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Semigroup for Sl<D> {}

impl<D: Dim> Monoid for Sl<D> {
    fn identity() -> Self {
        panic!("SL identity requires dimension; use Sl::identity(dim)")
    }
}

impl<D: Dim> Group for Sl<D> {
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse().expect("SL element must be invertible"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> TopologicalGroup for Sl<D> {}

impl<D: Dim> Manifold for Sl<D> {
    type Scalar = FiniteF64;
    fn dim(&self) -> usize { let n = self.inner.n(); n * n - 1 }
    fn atlas(&self) -> &Atlas<FiniteF64> { unimplemented!("SL atlas not yet constructed") }
}

impl<D: Dim> HasExpMap for Sl<D> {
    type Algebra = SlAlgebra<D>;
    fn exp(x: &SlAlgebra<D>) -> Self {
        Self { inner: x.inner.exp(), _dim: PhantomData }
    }
    fn log(&self) -> Option<SlAlgebra<D>> {
        Some(SlAlgebra { inner: self.inner.log()?, _dim: PhantomData })
    }
}
