// Ring through field 

use crate::core::ops::AbelianGroup;

pub trait Ring: AbelianGroup {
    fn mul(&self, other: &Self) -> Self;
    fn one() -> Self;
}

pub trait CommutativeRing: Ring {}

pub trait IntegralDomain: CommutativeRing {
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
}

pub trait Field: IntegralDomain {
    fn mul_inverse(&self) -> Option<Self>;
    fn div(&self, other: &Self) -> Option<Self> {
        other.mul_inverse().map(|inv| self.mul(&inv))
    }
}
