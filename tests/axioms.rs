use alice::core::ops::{AbelianGroup, Group, Magma, Monoid, Semigroup};
use alice::core::ring::{Field, Ring};
use alice::core::scalar::{FiniteF64, Rational};

fn check_closure<T: Magma>(a: &T, b: &T) -> bool {
    let _ = a.op(b);
    true
}

fn check_associativity<T: Semigroup>(a: &T, b: &T, c: &T) -> bool {
    a.op(&b.op(c)) == a.op(b).op(c)
}

fn check_identity_left<T: Monoid>(a: &T) -> bool {
    T::identity().op(a) == *a
}

fn check_identity_right<T: Monoid>(a: &T) -> bool {
    a.op(&T::identity()) == *a
}

fn check_inverse_left<T: Group>(a: &T) -> bool {
    a.inverse().op(a) == T::identity()
}

fn check_inverse_right<T: Group>(a: &T) -> bool {
    a.op(&a.inverse()) == T::identity()
}

fn check_commutativity<T: AbelianGroup>(a: &T, b: &T) -> bool {
    a.op(b) == b.op(a)
}

fn check_add_sub_roundtrip<T: AbelianGroup>(a: &T, b: &T) -> bool {
    a.add(b).sub(b) == *a
}

fn check_distributivity_left<T: Ring>(a: &T, b: &T, c: &T) -> bool {
    a.mul(&b.add(c)) == a.mul(b).add(&a.mul(c))
}

fn check_distributivity_right<T: Ring>(a: &T, b: &T, c: &T) -> bool {
    b.add(c).mul(a) == b.mul(a).add(&c.mul(a))
}

fn check_mul_identity_left<T: Ring>(a: &T) -> bool {
    T::one().mul(a) == *a
}

fn check_mul_identity_right<T: Ring>(a: &T) -> bool {
    a.mul(&T::one()) == *a
}

fn check_mul_associativity<T: Ring>(a: &T, b: &T, c: &T) -> bool {
    a.mul(&b.mul(c)) == a.mul(b).mul(c)
}

fn check_mul_inverse<T: Field>(a: &T) -> bool {
    if *a == T::zero() {
        a.mul_inverse().is_none()
    } else {
        match a.mul_inverse() {
            None => false,
            Some(inv) => a.mul(&inv) == T::one() && inv.mul(a) == T::one(),
        }
    }
}

fn check_div_exact<T: Field>(a: &T, b: &T) -> bool {
    if *b == T::zero() {
        a.div(b).is_none()
    } else {
        match a.div(b) {
            None => false,
            Some(q) => q.mul(b) == *a,
        }
    }
}

mod finite_f64_axioms {
    use super::*;
    use proptest::prelude::*;

    fn finite_f64_strategy() -> impl Strategy<Value = FiniteF64> {
        prop::num::f64::NORMAL
            .prop_filter("must be finite and within safe range for all ops", |v| {
                v.is_finite() && (*v < 1e15) && (*v > -1e15)
            })
            .prop_map(|v| FiniteF64::new(v).unwrap())
    }

    fn same_order_pair() -> impl Strategy<Value = (FiniteF64, FiniteF64)> {
        finite_f64_strategy().prop_flat_map(|a| {
            let a_abs = libm::fabs(a.get()).max(1.0);
            let lo = -(a_abs * 1e4);
            let hi = a_abs * 1e4;
            (Just(a), lo..=hi)
                .prop_filter("b must be finite", |(_, b)| b.is_finite())
                .prop_map(|(a, b)| (a, FiniteF64::new(b).unwrap_or(FiniteF64::new(0.0).unwrap())))
        })
    }

    fn safe_div_pair() -> impl Strategy<Value = (FiniteF64, FiniteF64)> {
        finite_f64_strategy().prop_flat_map(|a| {
            let a_abs = libm::fabs(a.get()).max(1.0);
            let lo = a_abs * 1e-4;
            let hi = a_abs * 1e4;
            (Just(a), lo..=hi)
                .prop_filter("b must be finite and nonzero", |(_, b)| {
                    b.is_finite() && *b != 0.0
                })
                .prop_map(|(a, b)| (a, FiniteF64::new(b).unwrap()))
        })
    }

    fn f64_near_relative(a: f64, b: f64, rel_eps: f64) -> bool {
        let diff = libm::fabs(a - b);
        let scale = libm::fabs(a).max(libm::fabs(b)).max(1.0);
        diff / scale < rel_eps
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn closure(a in finite_f64_strategy(), b in finite_f64_strategy()) {
            prop_assert!(check_closure(&a, &b));
        }

        #[test]
        fn associativity(
            a in finite_f64_strategy(),
            b in finite_f64_strategy(),
            c in finite_f64_strategy()
        ) {
            let lhs = a.op(&b.op(&c));
            let rhs = a.op(&b).op(&c);
            prop_assert!(
                f64_near_relative(lhs.get(), rhs.get(), 1e-9),
                "associativity violated: lhs={} rhs={}",
                lhs.get(),
                rhs.get()
            );
        }

        #[test]
        fn identity_left(a in finite_f64_strategy()) {
            prop_assert!(check_identity_left(&a));
        }

        #[test]
        fn identity_right(a in finite_f64_strategy()) {
            prop_assert!(check_identity_right(&a));
        }

        #[test]
        fn inverse_left(a in finite_f64_strategy()) {
            prop_assert!(check_inverse_left(&a));
        }

        #[test]
        fn inverse_right(a in finite_f64_strategy()) {
            prop_assert!(check_inverse_right(&a));
        }

        #[test]
        fn commutativity(a in finite_f64_strategy(), b in finite_f64_strategy()) {
            prop_assert!(check_commutativity(&a, &b));
        }

        #[test]
        fn add_sub_roundtrip(
            (a, b) in same_order_pair()
        ) {
            let result = a.add(&b).sub(&b);
            prop_assert!(
                f64_near_relative(result.get(), a.get(), 1e-9),
                "roundtrip failed: got {} expected {}",
                result.get(),
                a.get()
            );
        }

        #[test]
        fn distributivity_left(
            a in finite_f64_strategy(),
            b in finite_f64_strategy(),
            c in finite_f64_strategy()
        ) {
            let lhs = a.mul(&b.add(&c));
            let rhs = a.mul(&b).add(&a.mul(&c));
            prop_assert!(
                f64_near_relative(lhs.get(), rhs.get(), 1e-6),
                "left distributivity violated: lhs={} rhs={}",
                lhs.get(),
                rhs.get()
            );
        }

        #[test]
        fn distributivity_right(
            a in finite_f64_strategy(),
            b in finite_f64_strategy(),
            c in finite_f64_strategy()
        ) {
            let lhs = b.add(&c).mul(&a);
            let rhs = b.mul(&a).add(&c.mul(&a));
            prop_assert!(
                f64_near_relative(lhs.get(), rhs.get(), 1e-6),
                "right distributivity violated: lhs={} rhs={}",
                lhs.get(),
                rhs.get()
            );
        }

        #[test]
        fn mul_identity_left(a in finite_f64_strategy()) {
            prop_assert!(check_mul_identity_left(&a));
        }

        #[test]
        fn mul_identity_right(a in finite_f64_strategy()) {
            prop_assert!(check_mul_identity_right(&a));
        }

        #[test]
        fn mul_associativity(
            a in finite_f64_strategy(),
            b in finite_f64_strategy(),
            c in finite_f64_strategy()
        ) {
            let lhs = a.mul(&b.mul(&c));
            let rhs = a.mul(&b).mul(&c);
            prop_assert!(
                f64_near_relative(lhs.get(), rhs.get(), 1e-6),
                "mul associativity violated: lhs={} rhs={}",
                lhs.get(),
                rhs.get()
            );
        }

        #[test]
        fn mul_inverse_zero_returns_none(a in finite_f64_strategy()) {
            if a == FiniteF64::zero() {
                prop_assert!(a.mul_inverse().is_none());
            }
        }

        #[test]
        fn mul_inverse_nonzero_roundtrip(a in finite_f64_strategy()) {
            if a != FiniteF64::zero() {
                let inv = a.mul_inverse().unwrap();
                let product = a.mul(&inv);
                prop_assert!(
                    f64_near_relative(product.get(), FiniteF64::one().get(), 1e-9),
                    "a * a_inv != 1: got {}",
                    product.get()
                );
            }
        }

        #[test]
        fn div_nonzero_roundtrip(
            (a, b) in safe_div_pair()
        ) {
            let q = a.div(&b).unwrap();
            let recovered = q.mul(&b);
            prop_assert!(
                f64_near_relative(recovered.get(), a.get(), 1e-6),
                "div roundtrip failed: q*b={} a={}",
                recovered.get(),
                a.get()
            );
        }

        #[test]
        fn finitef64_is_always_finite(a in finite_f64_strategy()) {
            prop_assert!(a.get().is_finite());
        }

        #[test]
        fn finitef64_op_stays_finite(
            a in finite_f64_strategy(),
            b in finite_f64_strategy()
        ) {
            prop_assert!(a.op(&b).get().is_finite());
        }

        #[test]
        fn finitef64_mul_stays_finite(
            a in finite_f64_strategy(),
            b in finite_f64_strategy()
        ) {
            prop_assert!(a.mul(&b).get().is_finite());
        }
    }
}

mod rational_axioms {
    use super::*;
    use proptest::prelude::*;

    fn rational_strategy() -> impl Strategy<Value = Rational> {
        (i16::MIN as i32..=i16::MAX as i32, 1u32..=u16::MAX as u32)
            .prop_map(|(n, d)| Rational::new(n as i64, d as u64))
    }

    fn nonzero_rational_strategy() -> impl Strategy<Value = Rational> {
        (1i32..=i16::MAX as i32, 1u32..=u16::MAX as u32)
            .prop_map(|(n, d)| Rational::new(n as i64, d as u64))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn associativity(
            a in rational_strategy(),
            b in rational_strategy(),
            c in rational_strategy()
        ) {
            prop_assert!(check_associativity(&a, &b, &c));
        }

        #[test]
        fn identity_left(a in rational_strategy()) {
            prop_assert!(check_identity_left(&a));
        }

        #[test]
        fn identity_right(a in rational_strategy()) {
            prop_assert!(check_identity_right(&a));
        }

        #[test]
        fn inverse_left(a in rational_strategy()) {
            prop_assert!(check_inverse_left(&a));
        }

        #[test]
        fn inverse_right(a in rational_strategy()) {
            prop_assert!(check_inverse_right(&a));
        }

        #[test]
        fn commutativity(a in rational_strategy(), b in rational_strategy()) {
            prop_assert!(check_commutativity(&a, &b));
        }

        #[test]
        fn add_sub_roundtrip(a in rational_strategy(), b in rational_strategy()) {
            prop_assert!(check_add_sub_roundtrip(&a, &b));
        }

        #[test]
        fn distributivity_left(
            a in rational_strategy(),
            b in rational_strategy(),
            c in rational_strategy()
        ) {
            prop_assert!(check_distributivity_left(&a, &b, &c));
        }

        #[test]
        fn distributivity_right(
            a in rational_strategy(),
            b in rational_strategy(),
            c in rational_strategy()
        ) {
            prop_assert!(check_distributivity_right(&a, &b, &c));
        }

        #[test]
        fn mul_identity_left(a in rational_strategy()) {
            prop_assert!(check_mul_identity_left(&a));
        }

        #[test]
        fn mul_identity_right(a in rational_strategy()) {
            prop_assert!(check_mul_identity_right(&a));
        }

        #[test]
        fn mul_associativity(
            a in rational_strategy(),
            b in rational_strategy(),
            c in rational_strategy()
        ) {
            prop_assert!(check_mul_associativity(&a, &b, &c));
        }

        #[test]
        fn mul_inverse_field(a in rational_strategy()) {
            prop_assert!(check_mul_inverse(&a));
        }

        #[test]
        fn div_exact(a in rational_strategy(), b in rational_strategy()) {
            prop_assert!(check_div_exact(&a, &b));
        }

        #[test]
        fn mul_inverse_nonzero(a in nonzero_rational_strategy()) {
            let inv = a.mul_inverse().expect("nonzero rational must have inverse");
            let product = a.mul(&inv);
            prop_assert_eq!(product, Rational::one());
        }

        #[test]
        fn canonical_form_denom_positive(a in rational_strategy()) {
            prop_assert!(a.denom() > 0);
        }
    }
}

mod scalar_unit {
    use super::*;
    use alice::core::scalar::{Scalar, ScalarError};

    #[test]
    fn finitef64_rejects_nan() {
        assert_eq!(FiniteF64::new(f64::NAN), Err(ScalarError::NonFinite));
    }

    #[test]
    fn finitef64_rejects_pos_inf() {
        assert_eq!(FiniteF64::new(f64::INFINITY), Err(ScalarError::NonFinite));
    }

    #[test]
    fn finitef64_rejects_neg_inf() {
        assert_eq!(FiniteF64::new(f64::NEG_INFINITY), Err(ScalarError::NonFinite));
    }

    #[test]
    fn finitef64_accepts_zero() {
        assert!(FiniteF64::new(0.0).is_ok());
    }

    #[test]
    fn finitef64_accepts_negative() {
        assert!(FiniteF64::new(-1.0).is_ok());
    }

    #[test]
    fn finitef64_div_by_zero_returns_none() {
        let a = FiniteF64::new(1.0).unwrap();
        let b = FiniteF64::new(0.0).unwrap();
        assert_eq!(b.mul_inverse(), None);
        assert_eq!(a.div(&b), None);
    }

    #[test]
    fn finitef64_sqrt_negative_returns_none() {
        let a = FiniteF64::new(-1.0).unwrap();
        assert!(a.sqrt().is_none());
    }

    #[test]
    fn finitef64_sqrt_positive() {
        let a = FiniteF64::new(4.0).unwrap();
        let s = a.sqrt().unwrap();
        assert!(libm::fabs(s.get() - 2.0) < 1e-10);
    }

    #[test]
    fn rational_addition() {
        let a = Rational::new(1, 2);
        let b = Rational::new(1, 3);
        assert_eq!(a.add(&b), Rational::new(5, 6));
    }

    #[test]
    fn rational_subtraction() {
        let a = Rational::new(3, 4);
        let b = Rational::new(1, 4);
        assert_eq!(a.sub(&b), Rational::new(1, 2));
    }

    #[test]
    fn rational_multiplication() {
        let a = Rational::new(2, 3);
        let b = Rational::new(3, 4);
        assert_eq!(a.mul(&b), Rational::new(1, 2));
    }

    #[test]
    fn rational_inverse() {
        let a = Rational::new(3, 4);
        let inv = a.mul_inverse().unwrap();
        assert_eq!(inv, Rational::new(4, 3));
        assert_eq!(a.mul(&inv), Rational::one());
    }

    #[test]
    fn rational_zero_has_no_inverse() {
        assert!(Rational::new(0, 1).mul_inverse().is_none());
    }

    #[test]
    fn rational_negative_inverse() {
        let a = Rational::new(-2, 3);
        let inv = a.mul_inverse().unwrap();
        assert_eq!(inv, Rational::new(-3, 2));
    }

    #[test]
    fn rational_canonical_reduction() {
        assert_eq!(Rational::new(4, 6), Rational::new(2, 3));
    }

    #[test]
    fn rational_zero_identity() {
        let z = Rational::zero();
        let a = Rational::new(5, 7);
        assert_eq!(z.add(&a), a);
        assert_eq!(a.add(&z), a);
    }

    #[test]
    fn rational_one_identity() {
        let o = Rational::one();
        let a = Rational::new(5, 7);
        assert_eq!(o.mul(&a), a);
        assert_eq!(a.mul(&o), a);
    }
}

mod morphism_unit {
    use alice::core::morphism::{Composable, IdentityMorphism, MapMorphism, Morphism};

    #[test]
    fn identity_morphism() {
        let id: IdentityMorphism<i32> = IdentityMorphism::new();
        assert_eq!(id.apply(&42), 42);
    }

    #[test]
    fn map_morphism_apply() {
        let double = MapMorphism::new(|x: &i32| x * 2);
        assert_eq!(double.apply(&5), 10);
    }

    #[test]
    fn composition_applies_inner_then_outer() {
        let double = MapMorphism::new(|x: &i32| x * 2);
        let add_one = MapMorphism::new(|x: &i32| x + 1);
        let composed = double.compose(add_one);
        assert_eq!(composed.apply(&3), 8);
    }

    #[test]
    fn composition_associativity() {
        let f = MapMorphism::new(|x: &i32| x * 2);
        let g = MapMorphism::new(|x: &i32| x + 1);
        let h = MapMorphism::new(|x: &i32| x - 1);
        let fg_h = f.compose(g).compose(h);

        let f2 = MapMorphism::new(|x: &i32| x * 2);
        let g2 = MapMorphism::new(|x: &i32| x + 1);
        let h2 = MapMorphism::new(|x: &i32| x - 1);
        let f_gh = f2.compose(g2.compose(h2));

        for i in -10i32..=10 {
            assert_eq!(fg_h.apply(&i), f_gh.apply(&i));
        }
    }
}
