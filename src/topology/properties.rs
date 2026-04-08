use alloc::collections::BTreeSet;
use alloc::vec::Vec;
use crate::topology::space::{FiniteTopologicalSpace, TopologicalSpace};

pub trait Hausdorff: TopologicalSpace {
    fn are_separable(&self, p: &Self::Point, q: &Self::Point) -> bool;
}

pub trait Compact: TopologicalSpace {
    fn has_finite_subcover(
        &self,
        cover: &[BTreeSet<Self::Point>],
    ) -> bool;
}

pub trait Connected: TopologicalSpace {
    fn is_connected(&self) -> bool;
}

pub trait SecondCountable: TopologicalSpace {
    type Basis: IntoIterator<Item = BTreeSet<Self::Point>>;
    fn countable_basis(&self) -> Self::Basis;
}

impl<P: Clone + Eq + Ord> Hausdorff for FiniteTopologicalSpace<P> {
    fn are_separable(&self, p: &P, q: &P) -> bool {
        if p == q {
            return false;
        }
        for u in self.open_sets() {
            if u.contains(p) {
                for v in self.open_sets() {
                    if v.contains(q) {
                        let inter: BTreeSet<P> = u.intersection(v).cloned().collect();
                        if inter.is_empty() {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }
}

impl<P: Clone + Eq + Ord> FiniteTopologicalSpace<P> {
    pub fn is_hausdorff(&self) -> bool {
        let points: Vec<&P> = self.points().iter().collect();
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                if !self.are_separable(points[i], points[j]) {
                    return false;
                }
            }
        }
        true
    }
}

impl<P: Clone + Eq + Ord> Compact for FiniteTopologicalSpace<P> {
    fn has_finite_subcover(&self, cover: &[BTreeSet<P>]) -> bool {
        let covered: BTreeSet<P> = cover
            .iter()
            .flat_map(|s| s.iter().cloned())
            .collect();
        self.points().is_subset(&covered)
    }
}

impl<P: Clone + Eq + Ord> FiniteTopologicalSpace<P> {
    pub fn is_compact(&self) -> bool {
        let open_vec: Vec<BTreeSet<P>> = self.open_sets().iter().cloned().collect();
        let total_cover: BTreeSet<P> = open_vec
            .iter()
            .flat_map(|s| s.iter().cloned())
            .collect();
        if !self.points().is_subset(&total_cover) {
            return false;
        }
        self.has_finite_subcover(&open_vec)
    }
    pub fn is_connected(&self) -> bool {
        for u in self.open_sets() {
            if u.is_empty() || u == self.points() {
                continue;
            }
            let complement: BTreeSet<P> =
                self.points().difference(u).cloned().collect();
            if self.open_sets().contains(&complement) {
                return false;
            }
        }
        true
    }
}

impl<P: Clone + Eq + Ord> Connected for FiniteTopologicalSpace<P> {
    fn is_connected(&self) -> bool {
        FiniteTopologicalSpace::is_connected(self)
    }
}

impl<P: Clone + Eq + Ord> SecondCountable for FiniteTopologicalSpace<P> {
    type Basis = alloc::vec::Vec<BTreeSet<P>>;
    fn countable_basis(&self) -> Self::Basis {
        self.open_sets().iter().cloned().collect()
    }
}
