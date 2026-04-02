// Algebraic hierachy from Magma -> AlbeianGroup


pub trait Magma: Clone + Eq + Sized {
    fn op(&self, other: &Self) -> Self;
}

pub trait Semigroup: Magma {

}

pub trait Monoid: Magma {
    fn identity() -> Self;
}

pub trait Group: Monoid {
    fn inverse(&self) -> Self;
}

pub trait AlbeianGroup: Group {
   fn add(&self, other: &Self) -> Self{
       self.op(other)
   }

   fn zero() -> Self {
       Self::identity()
   }

   fn neg(&self) -> Self {
       self.inverse()
   }
   
   fn sub(&self, other: &Self) -> Self {
       self.add(&other.neg())
   }
}

pub trait TopologicalGroup: Group {}

pub trait FiniteGroup: Group {
    fn order() -> usize;
    fn elements() -> alloc::vec::Vec<Self>;
}

pub trait SimpleGroup: FiniteGroup {}
