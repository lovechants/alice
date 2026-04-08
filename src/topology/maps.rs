use alloc::boxed::Box;
use alloc::collections::BTreeSet;
use crate::core::morphism::Morphism;
use crate::topology::space::TopologicalSpace;

pub struct ContinuousMap<X, Y>
where
    X: TopologicalSpace,
    Y: TopologicalSpace,
{
    f: Box<dyn Fn(&X::Point) -> Y::Point + Send + Sync>,
}

impl<X, Y> ContinuousMap<X, Y>
where
    X: TopologicalSpace,
    Y: TopologicalSpace,
{
    pub fn new<F>(f: F) -> Self
    where
        F: Fn(&X::Point) -> Y::Point + Send + Sync + 'static,
    {
        Self { f: Box::new(f) }
    }

    pub fn apply(&self, p: &X::Point) -> Y::Point {
        (self.f)(p)
    }

    pub fn preimage(
        &self,
        domain_points: &BTreeSet<X::Point>,
        open_set: &BTreeSet<Y::Point>,
    ) -> BTreeSet<X::Point>
    where
        X::Point: Ord,
        Y::Point: Ord,
    {
        domain_points
            .iter()
            .filter(|p| open_set.contains(&self.apply(p)))
            .cloned()
            .collect()
    }

    pub fn is_continuous_on(
        &self,
        domain: &X,
        codomain: &Y,
        domain_points: &BTreeSet<X::Point>,
    ) -> bool
    where
        X::Point: Ord,
        Y::Point: Ord,
        Y: OpenSetsIter,
    {
        for open in codomain.open_sets_iter() {
            let pre = self.preimage(domain_points, &open);
            if !domain.is_open(&pre) {
                return false;
            }
        }
        true
    }
}

impl<X, Y> Morphism for ContinuousMap<X, Y>
where
    X: TopologicalSpace,
    Y: TopologicalSpace,
    X::Point: Clone,
    Y::Point: Clone,
{
    type Domain = X::Point;
    type Codomain = Y::Point;

    fn apply(&self, x: &X::Point) -> Y::Point {
        (self.f)(x)
    }
}

pub trait OpenSetsIter: TopologicalSpace {
    fn open_sets_iter(&self) -> alloc::vec::IntoIter<BTreeSet<Self::Point>>;
}

use crate::topology::space::FiniteTopologicalSpace;

impl<P: Clone + Eq + Ord> OpenSetsIter for FiniteTopologicalSpace<P> {
    fn open_sets_iter(&self) -> alloc::vec::IntoIter<BTreeSet<P>> {
        self.open_sets().iter().cloned().collect::<alloc::vec::Vec<_>>().into_iter()
    }
}
