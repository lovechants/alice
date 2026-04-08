use alloc::boxed::Box;
use crate::core::scalar::Scalar;
use crate::topology::manifold::{Atlas, Chart, Manifold, Point, TangentVector};

pub struct SmoothMap<S: Scalar> {
    domain_dim: usize,
    codomain_dim: usize,
    f: Box<dyn Fn(&Point<S>) -> Point<S> + Send + Sync>,
    differential: Box<dyn Fn(&Point<S>, &TangentVector<S>) -> TangentVector<S> + Send + Sync>,
}

impl<S: Scalar> SmoothMap<S> {
    pub fn new<F, D>(domain_dim: usize, codomain_dim: usize, f: F, differential: D) -> Self
    where
        F: Fn(&Point<S>) -> Point<S> + Send + Sync + 'static,
        D: Fn(&Point<S>, &TangentVector<S>) -> TangentVector<S> + Send + Sync + 'static,
    {
        Self {
            domain_dim,
            codomain_dim,
            f: Box::new(f),
            differential: Box::new(differential),
        }
    }

    pub fn apply(&self, p: &Point<S>) -> Point<S> {
        (self.f)(p)
    }

    pub fn differential(&self, p: &Point<S>, v: &TangentVector<S>) -> TangentVector<S> {
        (self.differential)(p, v)
    }

    pub fn domain_dim(&self) -> usize {
        self.domain_dim
    }

    pub fn codomain_dim(&self) -> usize {
        self.codomain_dim
    }
}

pub trait SmoothManifold: Manifold {
    fn transition_map(
        &self,
        from_chart: &Chart<Self::Scalar>,
        to_chart: &Chart<Self::Scalar>,
    ) -> SmoothMap<Self::Scalar>;

    fn pushforward(
        &self,
        map: &SmoothMap<Self::Scalar>,
        p: &Point<Self::Scalar>,
        v: &TangentVector<Self::Scalar>,
    ) -> TangentVector<Self::Scalar> {
        map.differential(p, v)
    }
}

pub struct EuclideanManifold<S: Scalar> {
    dim: usize,
    atlas: Atlas<S>,
}

impl<S: Scalar + 'static> EuclideanManifold<S> {
    pub fn new(dim: usize) -> Self {
        let mut atlas = Atlas::new(dim);
        let d = dim;
        atlas.add_chart(Chart::new(
            d,
            |p| p.clone(),
            |p| p.clone(),
            move |p| p.dim() == d,
        ));
        Self { dim, atlas }
    }
}

impl<S: Scalar> Manifold for EuclideanManifold<S> {
    type Scalar = S;

    fn dim(&self) -> usize {
        self.dim
    }

    fn atlas(&self) -> &Atlas<S> {
        &self.atlas
    }
}

impl<S: Scalar + 'static> SmoothManifold for EuclideanManifold<S> {
    fn transition_map(
        &self,
        _from_chart: &Chart<S>,
        _to_chart: &Chart<S>,
    ) -> SmoothMap<S> {
        let d = self.dim;
        SmoothMap::new(
            d,
            d,
            |p| p.clone(),
            |_p, v| v.clone(),
        )
    }
}

pub fn verify_transition_smoothness<M: SmoothManifold>(
    manifold: &M,
    test_points: &[Point<M::Scalar>],
) -> bool
where
    M::Scalar: PartialEq,
{
    let charts = manifold.atlas().charts();
    for i in 0..charts.len() {
        for j in 0..charts.len() {
            let t = manifold.transition_map(&charts[i], &charts[j]);
            for p in test_points {
                if charts[i].contains(p) && charts[j].contains(p) {
                    let q = t.apply(p);
                    if q.dim() != manifold.dim() {
                        return false;
                    }
                }
            }
        }
    }
    true
}

pub fn compose_smooth_maps<S: Scalar + 'static>(
    f: SmoothMap<S>,
    g: SmoothMap<S>,
) -> SmoothMap<S> {
    assert_eq!(
        g.codomain_dim(),
        f.domain_dim(),
        "map dimensions must be compatible for composition"
    );

    let domain_dim = g.domain_dim();
    let codomain_dim = f.codomain_dim();

    let f = alloc::sync::Arc::new(f);
    let g = alloc::sync::Arc::new(g);

    let f1 = f.clone();
    let g1 = g.clone();
    let f2 = f.clone();
    let g2 = g.clone();

    SmoothMap::new(
        domain_dim,
        codomain_dim,
        move |p| f1.apply(&g1.apply(p)),
        move |p, v| {
            let gp = g2.apply(p);
            let dg_v = g2.differential(p, v);
            f2.differential(&gp, &dg_v)
        },
    )
}
