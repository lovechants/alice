// modules + vector spaces 

use crate::core::ops::AbelianGroup;
use crate::core::ring::{Field, Ring};

pub trait Module<R: Ring>: AbelianGroup {
    fn scale(&self, r: &R) -> Self;
}

pub trait VectorSpace<F: Field>: Module<F> {}
