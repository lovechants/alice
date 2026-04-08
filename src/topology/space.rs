use alloc::collections::BTreeSet;
use alloc::vec::Vec;

pub trait TopologicalSpace {
    type Point: Clone + Eq;
    fn is_open(&self, subset: &BTreeSet<Self::Point>) -> bool;
    fn is_closed(&self, subset: &BTreeSet<Self::Point>) -> bool;
    fn closure(&self, subset: &BTreeSet<Self::Point>) -> BTreeSet<Self::Point>;
    fn interior(&self, subset: &BTreeSet<Self::Point>) -> BTreeSet<Self::Point>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteTopologicalSpace<P: Clone + Eq + Ord> {
    points: BTreeSet<P>,
    open_sets: BTreeSet<BTreeSet<P>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopologyError {
    MissingEmpty,
    MissingTotal,
    NotClosedUnderUnion,
    NotClosedUnderIntersection,
    OpenSetContainsExternalPoint,
}

impl<P: Clone + Eq + Ord> FiniteTopologicalSpace<P> {
    pub fn new(
        points: BTreeSet<P>,
        open_sets: BTreeSet<BTreeSet<P>>,
    ) -> Result<Self, TopologyError> {
        for open in &open_sets {
            for p in open {
                if !points.contains(p) {
                    return Err(TopologyError::OpenSetContainsExternalPoint);
                }
            }
        }

        if !open_sets.contains(&BTreeSet::new()) {
            return Err(TopologyError::MissingEmpty);
        }

        if !open_sets.contains(&points) {
            return Err(TopologyError::MissingTotal);
        }

        let open_vec: Vec<&BTreeSet<P>> = open_sets.iter().collect();
        let n = open_vec.len();

        for i in 0..n {
            for j in 0..n {
                let union: BTreeSet<P> = open_vec[i]
                    .union(open_vec[j])
                    .cloned()
                    .collect();
                if !open_sets.contains(&union) {
                    return Err(TopologyError::NotClosedUnderUnion);
                }

                let inter: BTreeSet<P> = open_vec[i]
                    .intersection(open_vec[j])
                    .cloned()
                    .collect();
                if !open_sets.contains(&inter) {
                    return Err(TopologyError::NotClosedUnderIntersection);
                }
            }
        }

        Ok(Self { points, open_sets })
    }

    pub fn indiscrete(points: BTreeSet<P>) -> Self {
        let mut open_sets = BTreeSet::new();
        open_sets.insert(BTreeSet::new());
        open_sets.insert(points.clone());
        Self { points, open_sets }
    }

    pub fn discrete(points: BTreeSet<P>) -> Self {
        let mut open_sets = BTreeSet::new();
        Self::all_subsets(&points, &mut open_sets);
        Self { points, open_sets }
    }

    pub fn points(&self) -> &BTreeSet<P> {
        &self.points
    }

    pub fn open_sets(&self) -> &BTreeSet<BTreeSet<P>> {
        &self.open_sets
    }

    fn all_subsets(set: &BTreeSet<P>, out: &mut BTreeSet<BTreeSet<P>>) {
        let items: Vec<&P> = set.iter().collect();
        let n = items.len();
        for mask in 0u64..(1u64 << n) {
            let subset: BTreeSet<P> = (0..n)
                .filter(|&i| mask & (1 << i) != 0)
                .map(|i| items[i].clone())
                .collect();
            out.insert(subset);
        }
    }

    pub fn verify_axioms(&self) -> Result<(), TopologyError> {
        Self::new(self.points.clone(), self.open_sets.clone()).map(|_| ())
    }
}

impl<P: Clone + Eq + Ord> TopologicalSpace for FiniteTopologicalSpace<P> {
    type Point = P;

    fn is_open(&self, subset: &BTreeSet<P>) -> bool {
        self.open_sets.contains(subset)
    }

    fn is_closed(&self, subset: &BTreeSet<P>) -> bool {
        let complement: BTreeSet<P> = self
            .points
            .difference(subset)
            .cloned()
            .collect();
        self.is_open(&complement)
    }

    fn closure(&self, subset: &BTreeSet<P>) -> BTreeSet<P> {
        self.open_sets
            .iter()
            .filter(|closed| {
                let complement: BTreeSet<P> = self
                    .points
                    .difference(*closed)
                    .cloned()
                    .collect();
                self.is_open(closed) && subset.is_subset(&complement)
            })
            .fold(self.points.clone(), |acc, closed| {
                let complement: BTreeSet<P> = self
                    .points
                    .difference(closed)
                    .cloned()
                    .collect();
                acc.intersection(&complement).cloned().collect()
            })
    }

    fn interior(&self, subset: &BTreeSet<P>) -> BTreeSet<P> {
        self.open_sets
            .iter()
            .filter(|open| open.is_subset(subset))
            .fold(BTreeSet::new(), |acc, open| {
                acc.union(open).cloned().collect()
            })
    }
}
