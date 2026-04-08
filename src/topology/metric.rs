use crate::core::scalar::Scalar;
use crate::core::ops::AbelianGroup;
use crate::core::ring::Field;

pub trait MetricSpace {
    type Point: Clone + Eq;
    type Distance: Scalar;
    fn distance(&self, a: &Self::Point, b: &Self::Point) -> Self::Distance;
    fn is_positive(&self, a: &Self::Point, b: &Self::Point) -> bool {
        if a == b {
            self.distance(a, b) == Self::Distance::zero()
        } else {
            self.distance(a, b) > Self::Distance::zero()
        }
    }
    fn is_symmetric(&self, a: &Self::Point, b: &Self::Point) -> bool {
        self.distance(a, b) == self.distance(b, a)
    }
    fn satisfies_triangle_inequality(
        &self,
        a: &Self::Point,
        b: &Self::Point,
        c: &Self::Point,
    ) -> bool {
        let ab = self.distance(a, b);
        let bc = self.distance(b, c);
        let ac = self.distance(a, c);
        ac <= ab.add(&bc)
    }
}

pub trait NormedSpace: Sized {
    type Scalar: Scalar;
    fn norm(&self) -> Self::Scalar;
    fn norm_positive(&self) -> bool
    where
        Self: Eq + AbelianGroup,
    {
        if *self == Self::zero() {
            self.norm() == Self::Scalar::zero()
        } else {
            self.norm() > Self::Scalar::zero()
        }
    }
    fn norm_homogeneous(&self, scalar: &Self::Scalar) -> Self::Scalar
    where
        Self: crate::core::module::VectorSpace<Self::Scalar>,
    {
        let scaled = self.scale(scalar);
        scaled.norm()
    }
}

pub trait InnerProductSpace: NormedSpace {
    fn inner_product(&self, other: &Self) -> Self::Scalar;
    fn norm_from_inner(&self) -> Option<Self::Scalar> {
        self.inner_product(self).sqrt()
    }
    fn is_orthogonal(&self, other: &Self) -> bool {
        self.inner_product(other) == Self::Scalar::zero()
    }
    fn is_symmetric_inner(&self, other: &Self) -> bool {
        self.inner_product(other) == other.inner_product(self)
    }
}
