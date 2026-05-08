use alloc::vec::Vec;
use core::marker::PhantomData;
use crate::core::ops::{Group, Magma, Monoid, Semigroup, TopologicalGroup};
use crate::groups::matrix_group::{GroupError, MatrixGroup};
use crate::groups::so::So;
use crate::topology::manifold::{Atlas, Dim, Manifold};
use crate::core::scalar::FiniteF64;
use crate::maps::exp_log::HasExpMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Se<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeAlgebra<D: Dim> {
    inner: MatrixGroup,
    _dim: PhantomData<D>,
}

impl<D: Dim> Se<D> {
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
            return Err(GroupError::ConstraintViolated("bottom right element must be 1"));
        }
        let rotation_data: Vec<f64> = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| data[i * nh + j])
            .collect();
        So::new(rotation_data, dim)?;
        let inner = MatrixGroup::new(data, nh)?;
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn from_parts(rotation: So<D>, translation: Vec<f64>, dim: D) -> Result<Self, GroupError> {
        let n = dim.value();
        if translation.len() != n {
            return Err(GroupError::DimensionMismatch);
        }
        let nh = n + 1;
        let mut data = alloc::vec![0.0f64; nh * nh];
        let rot_data = rotation.data();
        for i in 0..n {
            for j in 0..n {
                data[i * nh + j] = rot_data[i * n + j];
            }
            data[i * nh + n] = translation[i];
        }
        data[n * nh + n] = 1.0;
        let inner = MatrixGroup::new(data, nh)?;
        Ok(Self { inner, _dim: PhantomData })
    }

    pub fn identity(dim: D) -> Self {
        Self { inner: MatrixGroup::identity(dim.value() + 1), _dim: PhantomData }
    }

    pub fn rotation(&self, dim: D) -> So<D> {
        let n = dim.value();
        let _nh = n + 1;
        let data: Vec<f64> = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| self.inner.get(i, j))
            .collect();
        So::new(data, dim).expect("rotation block must be valid SO element")
    }

    pub fn translation(&self, dim: D) -> Vec<f64> {
        let n = dim.value();
        (0..n).map(|i| self.inner.get(i, n)).collect()
    }

    pub fn to_homogeneous(&self) -> Vec<f64> { self.inner.data() }
    pub fn n(&self) -> usize { self.inner.n() - 1 }
}

impl<D: Dim> SeAlgebra<D> {
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

impl<D: Dim> Magma for Se<D> {
    fn op(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.mul(&other.inner).expect("SE multiplication must succeed"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> Semigroup for Se<D> {}

impl<D: Dim> Monoid for Se<D> {
    fn identity() -> Self {
        panic!("SE identity requires dimension; use Se::identity(dim)")
    }
}

impl<D: Dim> Group for Se<D> {
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse().expect("SE element must be invertible"),
            _dim: PhantomData,
        }
    }
}

impl<D: Dim> TopologicalGroup for Se<D> {}

impl<D: Dim> Manifold for Se<D> {
    type Scalar = FiniteF64;
    fn dim(&self) -> usize {
        let n = self.inner.n() - 1;
        n * (n - 1) / 2 + n
    }
    fn atlas(&self) -> &Atlas<FiniteF64> { unimplemented!("SE atlas not yet constructed") }
}

impl<D: Dim> HasExpMap for Se<D> {
    type Algebra = SeAlgebra<D>;
    fn exp(x: &SeAlgebra<D>) -> Self {
        Self { inner: x.inner.exp(), _dim: PhantomData }
    }
    fn log(&self) -> Option<SeAlgebra<D>> {
        Some(SeAlgebra { inner: self.inner.log()?, _dim: PhantomData })
    }
}
