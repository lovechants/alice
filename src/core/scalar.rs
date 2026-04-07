// Finite f64 & Rational with full field

use crate::core::ops::{AbelianGroup, Group, Magma, Monoid, Semigroup};
use crate::core::ring::{CommutativeRing, Field, IntegralDomain, Ring};

#[derive(Debug, Clone, Copy)]
pub struct FiniteF64(f64);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarError {
    NonFinite,
    DivisionByZero,
    NegativeSqrt,
}

impl core::fmt::Display for ScalarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ScalarError::NonFinite => write!(f, "value is not finite"),
            ScalarError::DivisionByZero => write!(f, "division by zero"),
            ScalarError::NegativeSqrt => write!(f, "square root of negative number"),
        }
    }
}

impl FiniteF64 {
    pub fn new(v: f64) -> Result<Self, ScalarError> {
        if v.is_finite() {
            Ok(Self(v))
        } else {
            Err(ScalarError::NonFinite)
        }
    }

    pub fn get(&self) -> f64 {
        self.0
    }

    pub fn checked_add(self, other: Self) -> Result<Self, ScalarError> {
        Self::new(self.0 + other.0)
    }

    pub fn checked_sub(self, other: Self) -> Result<Self, ScalarError> {
        Self::new(self.0 - other.0)
    }

    pub fn checked_mul(self, other: Self) -> Result<Self, ScalarError> {
        Self::new(self.0 * other.0)
    }

    pub fn checked_div(self, other: Self) -> Result<Self, ScalarError> {
        if other.0 == 0.0 {
            Err(ScalarError::DivisionByZero)
        } else {
            Self::new(self.0 / other.0)
        }
    }

    pub fn checked_sqrt(self) -> Result<Self, ScalarError> {
        if self.0 < 0.0 {
            Err(ScalarError::NegativeSqrt)
        } else {
            Self::new(libm::sqrt(self.0))
        }
    }

    pub fn abs(self) -> Self {
        Self(libm::fabs(self.0))
    }
}

impl PartialEq for FiniteF64 {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for FiniteF64 {}

impl PartialOrd for FiniteF64 {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FiniteF64 {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.0.partial_cmp(&other.0).expect("FiniteF64 values are always comparable")
    }
}

impl Magma for FiniteF64 {
    fn op(&self, other: &Self) -> Self {
        Self::new(self.0 + other.0)
            .expect("addition of two finite f64 values must remain finite")
    }
}

impl Semigroup for FiniteF64 {}

impl Monoid for FiniteF64 {
    fn identity() -> Self {
        Self(0.0)
    }
}

impl Group for FiniteF64 {
    fn inverse(&self) -> Self {
        Self(-self.0)
    }
}

impl AbelianGroup for FiniteF64 {}

impl Ring for FiniteF64 {
    fn mul(&self, other: &Self) -> Self {
        Self::new(self.0 * other.0)
            .expect("multiplication of two finite f64 values must remain finite")
    }

    fn one() -> Self {
        Self(1.0)
    }
}

impl CommutativeRing for FiniteF64 {}

impl IntegralDomain for FiniteF64 {}

impl Field for FiniteF64 {
    fn mul_inverse(&self) -> Option<Self> {
        if self.0 == 0.0 {
            None
        } else {
            Self::new(1.0 / self.0).ok()
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Rational {
    numer: i64,
    denom: u64,
}

impl Rational {
    pub fn new(numer: i64, denom: u64) -> Self {
        assert!(denom != 0, "denominator must be nonzero");
        let g = gcd(numer.unsigned_abs(), denom);
        let sign = if numer < 0 { -1i64 } else { 1i64 };
        Self {
            numer: sign * (numer.unsigned_abs() / g) as i64,
            denom: denom / g,
        }
    }

    pub fn numer(&self) -> i64 {
        self.numer
    }

    pub fn denom(&self) -> u64 {
        self.denom
    }

    pub fn is_zero(&self) -> bool {
        self.numer == 0
    }

    fn checked_add_inner(self, other: Self) -> Option<Self> {
        let g = gcd(self.denom, other.denom);
        let self_scale = (other.denom / g) as i64;
        let other_scale = (self.denom / g) as i64;
        let lhs = self.numer.checked_mul(self_scale)?;
        let rhs = other.numer.checked_mul(other_scale)?;
        let n = lhs.checked_add(rhs)?;
        let d = self.denom.checked_mul(other.denom / g)?;
        Some(Self::new(n, d))
    }

    fn checked_mul_inner(self, other: Self) -> Option<Self> {
        let g1 = gcd(self.numer.unsigned_abs(), other.denom);
        let g2 = gcd(other.numer.unsigned_abs(), self.denom);
        let n1 = self.numer / g1 as i64;
        let n2 = other.numer / g2 as i64;
        let d1 = self.denom / g2;
        let d2 = other.denom / g1;
        let n = n1.checked_mul(n2)?;
        let d = d1.checked_mul(d2)?;
        Some(Self::new(n, d))
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        a %= b;
        core::mem::swap(&mut a, &mut b);
    }
    if a == 0 { 1 } else { a }
}

//fn lcm(a: u64, b: u64) -> Option<u64> {
//    a.checked_mul(b / gcd(a, b))
//}

impl Magma for Rational {
    fn op(&self, other: &Self) -> Self {
        self.checked_add_inner(*other)
            .expect("rational addition overflowed i64/u64 bounds")
    }
}

impl Semigroup for Rational {}

impl Monoid for Rational {
    fn identity() -> Self {
        Self::new(0, 1)
    }
}

impl Group for Rational {
    fn inverse(&self) -> Self {
        Self::new(-self.numer, self.denom)
    }
}

impl AbelianGroup for Rational {}

impl Ring for Rational {
    fn mul(&self, other: &Self) -> Self {
        self.checked_mul_inner(*other)
            .expect("rational multiplication overflowed i64/u64 bounds")
    }

    fn one() -> Self {
        Self::new(1, 1)
    }
}

impl CommutativeRing for Rational {}

impl IntegralDomain for Rational {}

impl Field for Rational {
    fn mul_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        if self.numer < 0 {
            Some(Self::new(-(self.denom as i64), self.numer.unsigned_abs()))
        } else {
            Some(Self::new(self.denom as i64, self.numer.unsigned_abs()))
        }
    }
}

pub trait Scalar: Field + Ord + core::fmt::Debug {
    fn abs(&self) -> Self;
    fn sqrt(&self) -> Option<Self>;
    fn from_f64(v: f64) -> Option<Self>;
    fn to_f64(&self) -> f64;
}

impl Scalar for FiniteF64 {
    fn abs(&self) -> Self {
        Self(libm::fabs(self.0))
    }

    fn sqrt(&self) -> Option<Self> {
        self.checked_sqrt().ok()
    }

    fn from_f64(v: f64) -> Option<Self> {
        Self::new(v).ok()
    }

    fn to_f64(&self) -> f64 {
        self.0
    }
}

impl Scalar for Rational {
    fn abs(&self) -> Self {
        Self::new(self.numer.abs(), self.denom)
    }

    fn sqrt(&self) -> Option<Self> {
        if self.numer < 0 {
            return None;
        }
        let v = libm::sqrt(self.numer as f64 / self.denom as f64);
        if v.is_finite() {
            let denom = 1_000_000u64;
            let numer = libm::round(v * denom as f64) as i64;
            Some(Self::new(numer, denom))
        } else {
            None
        }
    }

    fn from_f64(v: f64) -> Option<Self> {
        if !v.is_finite() {
            return None;
        }
        let denom = 1_000_000u64;
        let numer = libm::round(v * denom as f64) as i64;
        Some(Self::new(numer, denom))
    }

    fn to_f64(&self) -> f64 {
        self.numer as f64 / self.denom as f64
    }
}
