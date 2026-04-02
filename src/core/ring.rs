// Ring through field 

use crate::core::ops::AlbeianGroup;

pub trait Ring: AlbeianGroup {
    fn mul(&self, other: &Self) -> Self;
    fn one() -> Self;
}

pub trait CommunicativeRing: Ring {}

pub trait IntegralDomain: CommunicativeRing {
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
