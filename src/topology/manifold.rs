use alloc::boxed::Box;
use alloc::vec::Vec;
use crate::core::scalar::Scalar;

pub trait Dim: Copy + Eq + core::fmt::Debug {
    fn value(&self) -> usize;
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Const<const N: usize>;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Dynamic(usize);

impl<const N: usize> Dim for Const<N> {
    fn value(&self) -> usize {
        N
    }
}

impl Dynamic {
    pub fn new(n: usize) -> Self {
        Self(n)
    }
}

impl Dim for Dynamic {
    fn value(&self) -> usize {
        self.0
    }
}

pub struct Point<S: Scalar> {
    coords: Vec<S>,
}

impl<S: Scalar> Point<S> {
    pub fn new(coords: Vec<S>) -> Self {
        Self { coords }
    }

    pub fn dim(&self) -> usize {
        self.coords.len()
    }

    pub fn coords(&self) -> &[S] {
        &self.coords
    }

    pub fn get(&self, i: usize) -> Option<&S> {
        self.coords.get(i)
    }
}

impl<S: Scalar + PartialEq> PartialEq for Point<S> {
    fn eq(&self, other: &Self) -> bool {
        self.coords == other.coords
    }
}

impl<S: Scalar + Eq> Eq for Point<S> {}

impl<S: Scalar + Clone> Clone for Point<S> {
    fn clone(&self) -> Self {
        Self {
            coords: self.coords.clone(),
        }
    }
}

pub struct TangentVector<S: Scalar> {
    components: Vec<S>,
}

impl<S: Scalar> TangentVector<S> {
    pub fn new(components: Vec<S>) -> Self {
        Self { components }
    }

    pub fn dim(&self) -> usize {
        self.components.len()
    }

    pub fn components(&self) -> &[S] {
        &self.components
    }
}

impl<S: Scalar + Clone> Clone for TangentVector<S> {
    fn clone(&self) -> Self {
        Self {
            components: self.components.clone(),
        }
    }
}

pub struct Chart<S: Scalar> {
    dim: usize,
    to_euclidean: Box<dyn Fn(&Point<S>) -> Point<S> + Send + Sync>,
    from_euclidean: Box<dyn Fn(&Point<S>) -> Point<S> + Send + Sync>,
    domain_check: Box<dyn Fn(&Point<S>) -> bool + Send + Sync>,
}

impl<S: Scalar> Chart<S> {
    pub fn new<F, G, D>(dim: usize, to_euclidean: F, from_euclidean: G, domain_check: D) -> Self
    where
        F: Fn(&Point<S>) -> Point<S> + Send + Sync + 'static,
        G: Fn(&Point<S>) -> Point<S> + Send + Sync + 'static,
        D: Fn(&Point<S>) -> bool + Send + Sync + 'static,
    {
        Self {
            dim,
            to_euclidean: Box::new(to_euclidean),
            from_euclidean: Box::new(from_euclidean),
            domain_check: Box::new(domain_check),
        }
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn contains(&self, p: &Point<S>) -> bool {
        (self.domain_check)(p)
    }

    pub fn to_euclidean(&self, p: &Point<S>) -> Point<S> {
        (self.to_euclidean)(p)
    }

    pub fn from_euclidean(&self, p: &Point<S>) -> Point<S> {
        (self.from_euclidean)(p)
    }
}

pub struct Atlas<S: Scalar> {
    charts: Vec<Chart<S>>,
    dim: usize,
}

impl<S: Scalar> Atlas<S> {
    pub fn new(dim: usize) -> Self {
        Self {
            charts: Vec::new(),
            dim,
        }
    }

    pub fn add_chart(&mut self, chart: Chart<S>) {
        assert_eq!(
            chart.dim(),
            self.dim,
            "chart dimension must match atlas dimension"
        );
        self.charts.push(chart);
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn covers(&self, p: &Point<S>) -> bool {
        self.charts.iter().any(|c| c.contains(p))
    }

    pub fn chart_for(&self, p: &Point<S>) -> Option<&Chart<S>> {
        self.charts.iter().find(|c| c.contains(p))
    }

    pub fn charts(&self) -> &[Chart<S>] {
        &self.charts
    }
}

pub trait Manifold {
    type Scalar: Scalar;
    fn dim(&self) -> usize;
    fn atlas(&self) -> &Atlas<Self::Scalar>;
    fn tangent_space_dim(&self) -> usize {
        self.dim()
    }
    fn contains(&self, p: &Point<Self::Scalar>) -> bool {
        self.atlas().covers(p)
    }
    fn chart_for(&self, p: &Point<Self::Scalar>) -> Option<&Chart<Self::Scalar>> {
        self.atlas().chart_for(p)
    }
}
