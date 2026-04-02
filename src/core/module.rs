// modules + vector spaces 

use crate::core::ops::AlbeianGroup;
use crate::core::ring::{Field, Ring};

pub trait Module<R: Ring>: AlbeianGroup {
    fn scale(&self, r: &R) -> Self;
}

pub trait VectorSpace<F: Field>: Module<F> {}
