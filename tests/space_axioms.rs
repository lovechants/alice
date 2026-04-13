use alice::topology::space::{FiniteTopologicalSpace, TopologyError, TopologicalSpace};
use alice::topology::properties::{Hausdorff, SecondCountable};
use alice::topology::metric::MetricSpace;
use alice::topology::manifold::{
    Const, Dim, Dynamic, Manifold, Point, TangentVector,
};
use alice::topology::smooth::{
    EuclideanManifold, SmoothManifold, compose_smooth_maps, verify_transition_smoothness,
    SmoothMap,
};
use alice::core::scalar::FiniteF64;

use alloc::collections::BTreeSet;
use alloc::vec;
use alloc::vec::Vec;

extern crate alloc;

fn singleton(x: u32) -> BTreeSet<u32> {
    let mut s = BTreeSet::new();
    s.insert(x);
    s
}

fn set(xs: &[u32]) -> BTreeSet<u32> {
    xs.iter().cloned().collect()
}

fn open_sets(sets: &[&[u32]]) -> BTreeSet<BTreeSet<u32>> {
    sets.iter().map(|s| set(s)).collect()
}

mod finite_topology {
    use super::*;

    #[test]
    fn indiscrete_topology_is_valid() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::indiscrete(points.clone());
        assert!(space.verify_axioms().is_ok());
    }

    #[test]
    fn discrete_topology_is_valid() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points.clone());
        assert!(space.verify_axioms().is_ok());
    }

    #[test]
    fn discrete_contains_all_subsets() {
        let points = set(&[1, 2]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.is_open(&set(&[])));
        assert!(space.is_open(&set(&[1])));
        assert!(space.is_open(&set(&[2])));
        assert!(space.is_open(&set(&[1, 2])));
    }

    #[test]
    fn indiscrete_contains_only_empty_and_total() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::indiscrete(points.clone());
        assert!(space.is_open(&set(&[])));
        assert!(space.is_open(&points));
        assert!(!space.is_open(&set(&[1])));
        assert!(!space.is_open(&set(&[1, 2])));
    }

    #[test]
    fn custom_topology_valid() {
        let points = set(&[1, 2, 3]);
        let opens = open_sets(&[&[], &[1, 2, 3], &[1], &[2, 3]]);
        let space = FiniteTopologicalSpace::new(points, opens);
        assert!(space.is_ok());
    }

    #[test]
    fn rejects_missing_empty_set() {
        let points = set(&[1, 2]);
        let mut opens = BTreeSet::new();
        opens.insert(points.clone());
        let result = FiniteTopologicalSpace::new(points, opens);
        assert_eq!(result, Err(TopologyError::MissingEmpty));
    }

    #[test]
    fn rejects_missing_total_set() {
        let points = set(&[1, 2]);
        let mut opens = BTreeSet::new();
        opens.insert(set(&[]));
        let result = FiniteTopologicalSpace::new(points, opens);
        assert_eq!(result, Err(TopologyError::MissingTotal));
    }

    #[test]
    fn rejects_not_closed_under_union() {
        let points = set(&[1, 2, 3]);
        let opens = open_sets(&[&[], &[1, 2, 3], &[1], &[2]]);
        let result = FiniteTopologicalSpace::new(points, opens);
        assert_eq!(result, Err(TopologyError::NotClosedUnderUnion));
    }

    #[test]
    fn rejects_not_closed_under_intersection() {
        let points = set(&[1, 2, 3]);
        let opens = open_sets(&[&[], &[1, 2, 3], &[1, 2], &[2, 3]]);
        let result = FiniteTopologicalSpace::new(points, opens);
        assert_eq!(result, Err(TopologyError::NotClosedUnderIntersection));
    }

    #[test]
    fn rejects_external_points_in_open_set() {
        let points = set(&[1, 2]);
        let opens = open_sets(&[&[], &[1, 2], &[3]]);
        let result = FiniteTopologicalSpace::new(points, opens);
        assert_eq!(result, Err(TopologyError::OpenSetContainsExternalPoint));
    }

    #[test]
    fn is_open_empty_set() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.is_open(&set(&[])));
    }

    #[test]
    fn is_closed_complement_of_open() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.is_closed(&set(&[1])));
        assert!(space.is_closed(&set(&[2, 3])));
    }

    #[test]
    fn interior_of_open_set_is_itself() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points);
        let open = set(&[1, 2]);
        assert_eq!(space.interior(&open), open);
    }

    #[test]
    fn interior_of_non_open_is_smaller() {
        let points = set(&[1, 2, 3]);
        let opens = open_sets(&[&[], &[1, 2, 3], &[1]]);
        let space = FiniteTopologicalSpace::new(points, opens).unwrap();
        let subset = set(&[1, 2]);
        let interior = space.interior(&subset);
        assert!(interior.is_subset(&subset));
        assert_eq!(interior, set(&[1]));
    }
}

mod topology_axioms_proptest {
    use super::*;
    use proptest::prelude::*;

    fn small_points_strategy() -> impl Strategy<Value = BTreeSet<u32>> {
        prop::collection::btree_set(0u32..5, 2..5)
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(256))]

        #[test]
        fn discrete_always_valid(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::discrete(points);
            prop_assert!(space.verify_axioms().is_ok());
        }

        #[test]
        fn indiscrete_always_valid(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::indiscrete(points);
            prop_assert!(space.verify_axioms().is_ok());
        }

        #[test]
        fn discrete_open_sets_closed_under_union(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::discrete(points.clone());
            let opens: Vec<BTreeSet<u32>> = space.open_sets().iter().cloned().collect();
            for u in &opens {
                for v in &opens {
                    let union: BTreeSet<u32> = u.union(v).cloned().collect();
                    prop_assert!(space.is_open(&union));
                }
            }
        }

        #[test]
        fn discrete_open_sets_closed_under_intersection(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::discrete(points.clone());
            let opens: Vec<BTreeSet<u32>> = space.open_sets().iter().cloned().collect();
            for u in &opens {
                for v in &opens {
                    let inter: BTreeSet<u32> = u.intersection(v).cloned().collect();
                    prop_assert!(space.is_open(&inter));
                }
            }
        }

        #[test]
        fn interior_always_subset(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::discrete(points.clone());
            for subset in space.open_sets() {
                let interior = space.interior(subset);
                prop_assert!(interior.is_subset(subset));
            }
        }

        #[test]
        fn closed_complement_of_open(points in small_points_strategy()) {
            let space = FiniteTopologicalSpace::discrete(points.clone());
            for open in space.open_sets() {
                let complement: BTreeSet<u32> = points.difference(open).cloned().collect();
                prop_assert!(space.is_closed(&complement));
            }
        }
    }
}

mod hausdorff_tests {
    use super::*;

    #[test]
    fn discrete_is_hausdorff() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.is_hausdorff());
    }

    #[test]
    fn indiscrete_two_points_not_hausdorff() {
        let points = set(&[1, 2]);
        let space = FiniteTopologicalSpace::indiscrete(points);
        assert!(!space.is_hausdorff());
    }

    #[test]
    fn single_point_indiscrete_vacuously_hausdorff() {
        let points = singleton(1);
        let space = FiniteTopologicalSpace::indiscrete(points);
        assert!(space.is_hausdorff());
    }

    #[test]
    fn hausdorff_separates_distinct_points() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.are_separable(&1, &2));
        assert!(space.are_separable(&1, &3));
        assert!(space.are_separable(&2, &3));
    }
}

mod compact_connected_tests {
    use super::*;

    #[test]
    fn finite_space_is_compact() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points.clone());
        assert!(space.is_compact());
    }

    #[test]
    fn indiscrete_is_connected() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::indiscrete(points);
        assert!(space.is_connected());
    }

    #[test]
    fn discrete_two_points_not_connected() {
        let points = set(&[1, 2]);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(!space.is_connected());
    }

    #[test]
    fn single_point_discrete_is_connected() {
        let points = singleton(1);
        let space = FiniteTopologicalSpace::discrete(points);
        assert!(space.is_connected());
    }

    #[test]
    fn second_countable_basis_covers_space() {
        let points = set(&[1, 2, 3]);
        let space = FiniteTopologicalSpace::discrete(points.clone());
        let basis: Vec<BTreeSet<u32>> = space.countable_basis();
        let covered: BTreeSet<u32> = basis
            .iter()
            .flat_map(|s| s.iter().cloned())
            .collect();
        assert!(points.is_subset(&covered));
    }
}

mod dim_tests {
    use super::*;

    #[test]
    fn const_dim_value() {
        let d: Const<3> = Const;
        assert_eq!(d.value(), 3);
    }

    #[test]
    fn dynamic_dim_value() {
        let d = Dynamic::new(5);
        assert_eq!(d.value(), 5);
    }

    #[test]
    fn const_dims_equal() {
        let a: Const<4> = Const;
        let b: Const<4> = Const;
        assert_eq!(a, b);
    }

    #[test]
    fn dynamic_dims_equal() {
        let a = Dynamic::new(3);
        let b = Dynamic::new(3);
        assert_eq!(a, b);
    }

    #[test]
    fn dynamic_dims_unequal() {
        let a = Dynamic::new(3);
        let b = Dynamic::new(4);
        assert_ne!(a, b);
    }
}

mod manifold_tests {
    use super::*;

    //fn f64(v: f64) -> FiniteF64 {
    //    FiniteF64::new(v).unwrap()
    //}

    fn point(coords: &[f64]) -> Point<FiniteF64> {
        Point::new(coords.iter().map(|&v| FiniteF64::new(v).unwrap()).collect())
    }

    fn tangent(comps: &[f64]) -> TangentVector<FiniteF64> {
        TangentVector::new(comps.iter().map(|&v| FiniteF64::new(v).unwrap()).collect())
    }

    #[test]
    fn point_dim_correct() {
        let p = point(&[1.0, 2.0, 3.0]);
        assert_eq!(p.dim(), 3);
    }

    #[test]
    fn point_coords_accessible() {
        let p = point(&[1.0, 2.0]);
        assert_eq!(p.get(0).unwrap().get(), 1.0);
        assert_eq!(p.get(1).unwrap().get(), 2.0);
    }

    #[test]
    fn euclidean_manifold_dim() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(3);
        assert_eq!(m.dim(), 3);
    }

    #[test]
    fn euclidean_manifold_contains_correct_dim_points() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(3);
        let p = point(&[1.0, 2.0, 3.0]);
        assert!(m.contains(&p));
    }

    #[test]
    fn euclidean_manifold_rejects_wrong_dim() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(3);
        let p = point(&[1.0, 2.0]);
        assert!(!m.contains(&p));
    }

    #[test]
    fn atlas_covers_correct_dim_point() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(2);
        let p = point(&[0.5, 0.5]);
        assert!(m.atlas().covers(&p));
    }

    #[test]
    fn chart_roundtrip() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(2);
        let p = point(&[1.0, 2.0]);
        let chart = m.chart_for(&p).unwrap();
        let q = chart.to_euclidean(&p);
        let r = chart.from_euclidean(&q);
        assert_eq!(r.dim(), p.dim());
        for i in 0..p.dim() {
            assert!((r.get(i).unwrap().get() - p.get(i).unwrap().get()).abs() < 1e-10);
        }
    }

    #[test]
    fn transition_map_identity_on_euclidean() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(2);
        let p = point(&[1.0, 2.0]);
        let charts = m.atlas().charts();
        let t = m.transition_map(&charts[0], &charts[0]);
        let q = t.apply(&p);
        assert_eq!(q.dim(), 2);
    }

    #[test]
    fn smooth_map_applies_correctly() {
        let double: SmoothMap<FiniteF64> = SmoothMap::new(
            1,
            1,
            |p: &Point<FiniteF64>| Point::new(vec![FiniteF64::new(p.get(0).unwrap().get() * 2.0).unwrap()]),
            |_p: &Point<FiniteF64>, v: &TangentVector<FiniteF64>| TangentVector::new(vec![FiniteF64::new(v.components()[0].get() * 2.0).unwrap()]),
        );
        let p = point(&[3.0]);
        let q = double.apply(&p);
        assert!((q.get(0).unwrap().get() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn smooth_map_differential_correct() {
        let double: SmoothMap<FiniteF64> = SmoothMap::new(
            1,
            1,
            |p: &Point<FiniteF64>| Point::new(vec![FiniteF64::new(p.get(0).unwrap().get() * 2.0).unwrap()]),
            |_p: &Point<FiniteF64>, v: &TangentVector<FiniteF64>| TangentVector::new(vec![FiniteF64::new(v.components()[0].get() * 2.0).unwrap()]),
        );
        let p = point(&[1.0]);
        let v = tangent(&[1.0]);
        let dv = double.differential(&p, &v);
        assert!((dv.components()[0].get() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn compose_smooth_maps_applies_correctly() {
        let double: SmoothMap<FiniteF64> = SmoothMap::new(
            1, 1,
            |p: &Point<FiniteF64>| Point::new(vec![FiniteF64::new(p.get(0).unwrap().get() * 2.0).unwrap()]),
            |_p: &Point<FiniteF64>, v: &TangentVector<FiniteF64>| TangentVector::new(vec![FiniteF64::new(v.components()[0].get() * 2.0).unwrap()]),
        );
        let add_one: SmoothMap<FiniteF64> = SmoothMap::new(
            1, 1,
            |p: &Point<FiniteF64>| Point::new(vec![FiniteF64::new(p.get(0).unwrap().get() + 1.0).unwrap()]),
            |_p: &Point<FiniteF64>, v: &TangentVector<FiniteF64>| v.clone(),
        );
        let composed = compose_smooth_maps(double, add_one);
        let p = point(&[3.0]);
        let q = composed.apply(&p);
        assert!((q.get(0).unwrap().get() - 8.0).abs() < 1e-10);
    }

    #[test]
    fn verify_transition_smoothness_euclidean() {
        let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(2);
        let test_points = vec![
            point(&[0.0, 0.0]),
            point(&[1.0, 0.0]),
            point(&[0.0, 1.0]),
        ];
        assert!(verify_transition_smoothness(&m, &test_points));
    }

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(256))]

        #[test]
        fn euclidean_chart_roundtrip_proptest(
            x in -1e6f64..1e6f64,
            y in -1e6f64..1e6f64
        ) {
            let m: EuclideanManifold<FiniteF64> = EuclideanManifold::new(2);
            let p = Point::new(vec![
                FiniteF64::new(x).unwrap(),
                FiniteF64::new(y).unwrap(),
            ]);
            let chart = m.chart_for(&p).unwrap();
            let q = chart.from_euclidean(&chart.to_euclidean(&p));
            for i in 0..2 {
                proptest::prop_assert!(
                    libm::fabs(q.get(i).unwrap().get() - p.get(i).unwrap().get()) < 1e-10
                );
            }
        }

        #[test]
        fn smooth_map_differential_linearity(
            ax in -1e3f64..1e3f64,
            bx in -1e3f64..1e3f64,
            _s in -1e3f64..1e3f64,
        ) {
            let double: SmoothMap<FiniteF64> = SmoothMap::new(
                1, 1,
                |p: &Point<FiniteF64>| Point::new(vec![FiniteF64::new(p.get(0).unwrap().get() * 2.0).unwrap()]),
                |_p: &Point<FiniteF64>, v: &TangentVector<FiniteF64>| TangentVector::new(vec![FiniteF64::new(v.components()[0].get() * 2.0).unwrap()]),
            );
            let p = Point::new(vec![FiniteF64::new(0.0).unwrap()]);
            let va = TangentVector::new(vec![FiniteF64::new(ax).unwrap()]);
            let vb = TangentVector::new(vec![FiniteF64::new(bx).unwrap()]);

            let d_a = double.differential(&p, &va).components()[0].get();
            let d_b = double.differential(&p, &vb).components()[0].get();

            let vab = TangentVector::new(vec![FiniteF64::new(ax + bx).unwrap()]);
            let d_ab = double.differential(&p, &vab).components()[0].get();

            proptest::prop_assert!(libm::fabs(d_ab - (d_a + d_b)) < 1e-9);
        }
    }
}

mod metric_tests {
    use super::*;
    use proptest::prelude::*;

    struct EuclideanMetric;

    impl MetricSpace for EuclideanMetric {
        type Point = (FiniteF64, FiniteF64);
        type Distance = FiniteF64;

        fn distance(&self, a: &Self::Point, b: &Self::Point) -> FiniteF64 {
            let dx = a.0.get() - b.0.get();
            let dy = a.1.get() - b.1.get();
            FiniteF64::new(libm::sqrt(dx * dx + dy * dy)).unwrap()
        }
    }

    fn fv(v: f64) -> FiniteF64 { FiniteF64::new(v).unwrap() }
    fn pt(x: f64, y: f64) -> (FiniteF64, FiniteF64) { (fv(x), fv(y)) }

    #[test]
    fn distance_self_is_zero() {
        let m = EuclideanMetric;
        let p = pt(1.0, 2.0);
        assert!(m.is_positive(&p, &p));
    }

    #[test]
    fn distance_is_symmetric() {
        let m = EuclideanMetric;
        let p = pt(1.0, 2.0);
        let q = pt(3.0, 4.0);
        assert!(m.is_symmetric(&p, &q));
    }

    #[test]
    fn triangle_inequality_holds() {
        let m = EuclideanMetric;
        let a = pt(0.0, 0.0);
        let b = pt(3.0, 0.0);
        let c = pt(0.0, 4.0);
        assert!(m.satisfies_triangle_inequality(&a, &b, &c));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(256))]

        #[test]
        fn metric_symmetry(
            ax in -1e6f64..1e6f64, ay in -1e6f64..1e6f64,
            bx in -1e6f64..1e6f64, by in -1e6f64..1e6f64,
        ) {
            let m = EuclideanMetric;
            let p = pt(ax, ay);
            let q = pt(bx, by);
            let dpq = m.distance(&p, &q).get();
            let dqp = m.distance(&q, &p).get();
            prop_assert!(libm::fabs(dpq - dqp) < 1e-10);
        }

        #[test]
        fn metric_non_negative(
            ax in -1e6f64..1e6f64, ay in -1e6f64..1e6f64,
            bx in -1e6f64..1e6f64, by in -1e6f64..1e6f64,
        ) {
            let m = EuclideanMetric;
            let p = pt(ax, ay);
            let q = pt(bx, by);
            prop_assert!(m.distance(&p, &q).get() >= 0.0);
        }

        #[test]
        fn metric_triangle_inequality(
            ax in -1e3f64..1e3f64, ay in -1e3f64..1e3f64,
            bx in -1e3f64..1e3f64, by in -1e3f64..1e3f64,
            cx in -1e3f64..1e3f64, cy in -1e3f64..1e3f64,
        ) {
            let m = EuclideanMetric;
            let a = pt(ax, ay);
            let b = pt(bx, by);
            let c = pt(cx, cy);
            prop_assert!(m.satisfies_triangle_inequality(&a, &b, &c));
        }
    }
}
