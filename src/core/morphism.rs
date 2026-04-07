use alloc::boxed::Box;

pub trait Morphism {
    type Domain: Clone;
    type Codomain: Clone;

    fn apply(&self, x: &Self::Domain) -> Self::Codomain;
}

pub trait Composable<M>: Morphism
where
    M: Morphism<Codomain = Self::Domain>,
{
    type Composed: Morphism<Domain = M::Domain, Codomain = Self::Codomain>;
    fn compose(self, other: M) -> Self::Composed;
}

pub trait Invertible: Morphism {
    type Inverse: Morphism<Domain = Self::Codomain, Codomain = Self::Domain>;
    fn inverse(self) -> Self::Inverse;
}

pub struct MapMorphism<A, B> {
    f: Box<dyn Fn(&A) -> B + Send + Sync>,
}

impl<A: Clone, B: Clone> MapMorphism<A, B> {
    pub fn new<F>(f: F) -> Self
    where
        F: Fn(&A) -> B + Send + Sync + 'static,
    {
        Self { f: Box::new(f) }
    }
}

impl<A: Clone, B: Clone> Morphism for MapMorphism<A, B> {
    type Domain = A;
    type Codomain = B;

    fn apply(&self, x: &A) -> B {
        (self.f)(x)
    }
}

pub struct ComposedMorphism<F, G> {
    outer: F,
    inner: G,
}

impl<F, G> Morphism for ComposedMorphism<F, G>
where
    F: Morphism,
    G: Morphism<Codomain = F::Domain>,
{
    type Domain = G::Domain;
    type Codomain = F::Codomain;

    fn apply(&self, x: &Self::Domain) -> Self::Codomain {
        let intermediate = self.inner.apply(x);
        self.outer.apply(&intermediate)
    }
}

impl<F, G> Composable<G> for F
where
    F: Morphism,
    G: Morphism<Codomain = F::Domain>,
{
    type Composed = ComposedMorphism<F, G>;

    fn compose(self, other: G) -> Self::Composed {
        ComposedMorphism {
            outer: self,
            inner: other,
        }
    }
}

pub struct IdentityMorphism<A> {
    _marker: core::marker::PhantomData<A>,
}

impl<A> IdentityMorphism<A> {
    pub fn new() -> Self {
        Self {
            _marker: core::marker::PhantomData,
        }
    }
}

impl<A> Default for IdentityMorphism<A> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: Clone> Morphism for IdentityMorphism<A> {
    type Domain = A;
    type Codomain = A;

    fn apply(&self, x: &A) -> A {
        x.clone()
    }
}

impl<A: Clone> Invertible for IdentityMorphism<A> {
    type Inverse = IdentityMorphism<A>;

    fn inverse(self) -> Self::Inverse {
        Self::new()
    }
}
