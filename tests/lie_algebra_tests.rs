use alice::core::scalar::FiniteF64;
use alice::core::ops::{AbelianGroup, Group, Magma, Monoid, Semigroup};
use alice::core::ring::Ring;
use alice::core::module::{Module, VectorSpace};
use alice::algebra::lie_algebra::LieAlgebra;
use alice::algebra::bch::{
    bch_commuting_check, bch_first_order, bch_inverse_check, bch_second_order, bch_third_order,
};

fn f(v: f64) -> FiniteF64 {
    FiniteF64::new(v).unwrap()
}

/// Element of so(2): real skew-symmetric 2x2 matrix.
/// Parameterized by a single value a: [[0, -a], [a, 0]].
/// so(2) is one-dimensional and abelian: bracket is always zero.
#[derive(Clone, Debug, PartialEq, Eq)]
struct So2(FiniteF64);

impl Magma for So2 {
    fn op(&self, other: &Self) -> Self {
        So2(self.0.add(&other.0))
    }
}

impl Semigroup for So2 {}

impl Monoid for So2 {
    fn identity() -> Self {
        So2(FiniteF64::zero())
    }
}

impl Group for So2 {
    fn inverse(&self) -> Self {
        So2(self.0.neg())
    }
}

impl AbelianGroup for So2 {}

impl Module<FiniteF64> for So2 {
    fn scale(&self, r: &FiniteF64) -> Self {
        So2(self.0.mul(r))
    }
}

impl VectorSpace<FiniteF64> for So2 {}

impl LieAlgebra<FiniteF64> for So2 {
    fn bracket(&self, _other: &Self) -> Self {
        So2(FiniteF64::zero())
    }
}

/// Element of so(3): real skew-symmetric 3x3 matrix.
/// Represented as (a, b, c) corresponding to:
/// [[0, -c, b], [c, 0, -a], [-b, a, 0]]
/// Bracket is the cross product: [X,Y] = X x Y.
#[derive(Clone, Debug, PartialEq, Eq)]
struct So3(FiniteF64, FiniteF64, FiniteF64);

impl So3 {
    fn new(a: f64, b: f64, c: f64) -> Self {
        So3(f(a), f(b), f(c))
    }
}

impl Magma for So3 {
    fn op(&self, other: &Self) -> Self {
        So3(
            self.0.add(&other.0),
            self.1.add(&other.1),
            self.2.add(&other.2),
        )
    }
}

impl Semigroup for So3 {}

impl Monoid for So3 {
    fn identity() -> Self {
        So3(FiniteF64::zero(), FiniteF64::zero(), FiniteF64::zero())
    }
}

impl Group for So3 {
    fn inverse(&self) -> Self {
        So3(self.0.neg(), self.1.neg(), self.2.neg())
    }
}

impl AbelianGroup for So3 {}

impl Module<FiniteF64> for So3 {
    fn scale(&self, r: &FiniteF64) -> Self {
        So3(self.0.mul(r), self.1.mul(r), self.2.mul(r))
    }
}

impl VectorSpace<FiniteF64> for So3 {}

impl LieAlgebra<FiniteF64> for So3 {
    fn bracket(&self, other: &Self) -> Self {
        let a1 = &self.0;
        let b1 = &self.1;
        let c1 = &self.2;
        let a2 = &other.0;
        let b2 = &other.1;
        let c2 = &other.2;
        So3(
            b1.mul(c2).add(&b2.mul(c1).neg()),
            c1.mul(a2).add(&c2.mul(a1).neg()),
            a1.mul(b2).add(&a2.mul(b1).neg()),
        )
    }
}

mod so2_axioms {
    use super::*;
    use proptest::prelude::*;

    fn so2_strategy() -> impl Strategy<Value = So2> {
        (-1e6f64..1e6f64)
            .prop_filter("must be finite", |v| v.is_finite())
            .prop_map(|v| So2(f(v)))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn antisymmetry(a in so2_strategy(), b in so2_strategy()) {
            prop_assert!(a.antisymmetry_holds(&b));
        }

        #[test]
        fn jacobi(a in so2_strategy(), b in so2_strategy(), c in so2_strategy()) {
            prop_assert!(a.jacobi_holds(&b, &c));
        }

        #[test]
        fn bracket_bilinear_left(
            a in so2_strategy(),
            b in so2_strategy(),
            s in (-1e3f64..1e3f64).prop_filter("finite", |v| v.is_finite()).prop_map(|v| f(v))
        ) {
            prop_assert!(a.bracket_bilinear_left(&b, &s));
        }

        #[test]
        fn bracket_bilinear_right(
            a in so2_strategy(),
            b in so2_strategy(),
            s in (-1e3f64..1e3f64).prop_filter("finite", |v| v.is_finite()).prop_map(|v| f(v))
        ) {
            prop_assert!(a.bracket_bilinear_right(&b, &s));
        }

        #[test]
        fn bracket_additive_left(
            a in so2_strategy(),
            b in so2_strategy(),
            c in so2_strategy()
        ) {
            prop_assert!(a.bracket_additive_left(&b, &c));
        }

        #[test]
        fn bracket_additive_right(
            a in so2_strategy(),
            b in so2_strategy(),
            c in so2_strategy()
        ) {
            prop_assert!(a.bracket_additive_right(&b, &c));
        }
    }

    #[test]
    fn so2_bracket_is_always_zero() {
        let x = So2(f(1.0));
        let y = So2(f(2.0));
        assert_eq!(x.bracket(&y), So2(FiniteF64::zero()));
    }

    #[test]
    fn so2_bracket_self_is_zero() {
        let x = So2(f(3.14));
        assert_eq!(x.bracket(&x), So2(FiniteF64::zero()));
    }
}

mod so3_axioms {
    use super::*;
    use proptest::prelude::*;

    fn bounded_f64() -> impl Strategy<Value = f64> {
        (-1e3f64..1e3f64).prop_filter("finite", |v| v.is_finite())
    }

    fn so3_strategy() -> impl Strategy<Value = So3> {
        (bounded_f64(), bounded_f64(), bounded_f64())
            .prop_map(|(a, b, c)| So3::new(a, b, c))
    }

    fn f64_near(a: f64, b: f64) -> bool {
        libm::fabs(a - b) < 1e-6
    }

    fn so3_near(x: &So3, y: &So3) -> bool {
        f64_near(x.0.get(), y.0.get())
            && f64_near(x.1.get(), y.1.get())
            && f64_near(x.2.get(), y.2.get())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn antisymmetry(a in so3_strategy(), b in so3_strategy()) {
            let xy = a.bracket(&b);
            let yx = b.bracket(&a);
            prop_assert!(so3_near(&xy, &yx.inverse()));
        }

        #[test]
        fn jacobi(
            a in so3_strategy(),
            b in so3_strategy(),
            c in so3_strategy()
        ) {
            let xy_z = a.bracket(&b).bracket(&c);
            let yz_x = b.bracket(&c).bracket(&a);
            let zx_y = c.bracket(&a).bracket(&b);
            let sum = xy_z.add(&yz_x).add(&zx_y);
            prop_assert!(
                so3_near(&sum, &So3::zero()),
                "Jacobi violated: sum = ({}, {}, {})",
                sum.0.get(), sum.1.get(), sum.2.get()
            );
        }

        #[test]
        fn bracket_bilinear_left(
            a in so3_strategy(),
            b in so3_strategy(),
            s in bounded_f64().prop_map(|v| f(v))
        ) {
            let lhs = a.scale(&s).bracket(&b);
            let rhs = a.bracket(&b).scale(&s);
            prop_assert!(so3_near(&lhs, &rhs));
        }

        #[test]
        fn bracket_bilinear_right(
            a in so3_strategy(),
            b in so3_strategy(),
            s in bounded_f64().prop_map(|v| f(v))
        ) {
            let lhs = a.bracket(&b.scale(&s));
            let rhs = a.bracket(&b).scale(&s);
            prop_assert!(so3_near(&lhs, &rhs));
        }

        #[test]
        fn bracket_additive_left(
            a in so3_strategy(),
            b in so3_strategy(),
            c in so3_strategy()
        ) {
            let lhs = a.add(&b).bracket(&c);
            let rhs = a.bracket(&c).add(&b.bracket(&c));
            prop_assert!(so3_near(&lhs, &rhs));
        }

        #[test]
        fn bracket_additive_right(
            a in so3_strategy(),
            b in so3_strategy(),
            c in so3_strategy()
        ) {
            let lhs = a.bracket(&b.add(&c));
            let rhs = a.bracket(&b).add(&a.bracket(&c));
            prop_assert!(so3_near(&lhs, &rhs));
        }

        #[test]
        fn bracket_self_is_zero(a in so3_strategy()) {
            let result = a.bracket(&a);
            prop_assert!(so3_near(&result, &So3::zero()));
        }
    }

    #[test]
    fn so3_basis_bracket_e1_e2() {
        let e1 = So3::new(1.0, 0.0, 0.0);
        let e2 = So3::new(0.0, 1.0, 0.0);
        let e3 = So3::new(0.0, 0.0, 1.0);
        assert_eq!(e1.bracket(&e2), e3);
    }

    #[test]
    fn so3_basis_bracket_e2_e3() {
        let e1 = So3::new(1.0, 0.0, 0.0);
        let e2 = So3::new(0.0, 1.0, 0.0);
        let e3 = So3::new(0.0, 0.0, 1.0);
        assert_eq!(e2.bracket(&e3), e1);
    }

    #[test]
    fn so3_basis_bracket_e3_e1() {
        let e1 = So3::new(1.0, 0.0, 0.0);
        let e2 = So3::new(0.0, 1.0, 0.0);
        let e3 = So3::new(0.0, 0.0, 1.0);
        assert_eq!(e3.bracket(&e1), e2);
    }

    #[test]
    fn so3_bracket_antisymmetric_on_basis() {
        let e1 = So3::new(1.0, 0.0, 0.0);
        let e2 = So3::new(0.0, 1.0, 0.0);
        assert_eq!(e1.bracket(&e2), e2.bracket(&e1).inverse());
    }

    #[test]
    fn so3_is_non_abelian() {
        let e1 = So3::new(1.0, 0.0, 0.0);
        let e2 = So3::new(0.0, 1.0, 0.0);
        assert_ne!(e1.bracket(&e2), So3::zero());
    }
}

mod bch_tests {
    use super::*;
    use proptest::prelude::*;

    fn small_so3() -> impl Strategy<Value = So3> {
        (-0.1f64..0.1f64, -0.1f64..0.1f64, -0.1f64..0.1f64)
            .prop_filter("all finite", |(a, b, c)| {
                a.is_finite() && b.is_finite() && c.is_finite()
            })
            .prop_map(|(a, b, c)| So3::new(a, b, c))
    }

    fn so3_near(x: &So3, y: &So3, eps: f64) -> bool {
        libm::fabs(x.0.get() - y.0.get()) < eps
            && libm::fabs(x.1.get() - y.1.get()) < eps
            && libm::fabs(x.2.get() - y.2.get()) < eps
    }

    #[test]
    fn bch_zero_zero_is_zero() {
        let z = So3::zero();
        let result = bch_third_order::<FiniteF64, So3>(&z, &z);
        assert_eq!(result, So3::zero());
    }

    #[test]
    fn bch_x_zero_is_x() {
        let x = So3::new(0.05, 0.0, 0.0);
        let z = So3::zero();
        let result = bch_third_order::<FiniteF64, So3>(&x, &z);
        assert!(so3_near(&result, &x, 1e-10));
    }

    #[test]
    fn bch_zero_y_is_y() {
        let z = So3::zero();
        let y = So3::new(0.0, 0.05, 0.0);
        let result = bch_third_order::<FiniteF64, So3>(&z, &y);
        assert!(so3_near(&result, &y, 1e-10));
    }

    #[test]
    fn bch_inverse_gives_zero() {
        let x = So3::new(0.05, 0.02, 0.01);
        assert!(bch_inverse_check::<FiniteF64, So3>(&x));
    }

    #[test]
    fn bch_commuting_elements_give_sum() {
        let x = So3::new(0.1, 0.0, 0.0);
        let y = So3::new(0.2, 0.0, 0.0);
        assert!(bch_commuting_check::<FiniteF64, So3>(&x, &y));
    }

    #[test]
    fn bch_orders_agree_for_commuting_elements() {
        let x = So3::new(0.05, 0.0, 0.0);
        let y = So3::new(0.03, 0.0, 0.0);
        let first = bch_first_order::<FiniteF64, So3>(&x, &y);
        let second = bch_second_order::<FiniteF64, So3>(&x, &y);
        let third = bch_third_order::<FiniteF64, So3>(&x, &y);
        assert!(so3_near(&first, &second, 1e-10));
        assert!(so3_near(&second, &third, 1e-10));
    }

    #[test]
    fn bch_noncommuting_second_differs_from_first() {
        let x = So3::new(0.1, 0.0, 0.0);
        let y = So3::new(0.0, 0.1, 0.0);
        let first = bch_first_order::<FiniteF64, So3>(&x, &y);
        let second = bch_second_order::<FiniteF64, So3>(&x, &y);
        assert!(!so3_near(&first, &second, 1e-10));
    }

    #[test]
    fn bch_bracket_correction_has_correct_sign() {
        let x = So3::new(0.1, 0.0, 0.0);
        let y = So3::new(0.0, 0.1, 0.0);
        let result = bch_second_order::<FiniteF64, So3>(&x, &y);
        let xy_bracket = x.bracket(&y);
        let expected = x.add(&y).add(&xy_bracket.scale(&f(0.5)));
        assert!(so3_near(&result, &expected, 1e-10));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn bch_first_order_is_sum(x in small_so3(), y in small_so3()) {
            let result = bch_first_order::<FiniteF64, So3>(&x, &y);
            let expected = x.add(&y);
            prop_assert!(so3_near(&result, &expected, 1e-12));
        }

        #[test]
        fn bch_second_contains_first(x in small_so3(), y in small_so3()) {
            let first = bch_first_order::<FiniteF64, So3>(&x, &y);
            let second = bch_second_order::<FiniteF64, So3>(&x, &y);
            let correction = second.add(&first.inverse());
            let half_bracket = x.bracket(&y).scale(&f(0.5));
            prop_assert!(so3_near(&correction, &half_bracket, 1e-9));
        }

        #[test]
        fn bch_third_order_x_self_is_double(x in small_so3()) {
            let result = bch_third_order::<FiniteF64, So3>(&x, &x);
            let expected = x.scale(&f(2.0));
            prop_assert!(
                so3_near(&result, &expected, 1e-9),
                "bch(x,x) != 2x: result=({},{},{}), expected=({},{},{})",
                result.0.get(), result.1.get(), result.2.get(),
                expected.0.get(), expected.1.get(), expected.2.get()
            );
        }

        #[test]
        fn bch_symmetry_abelian(
            a in (-0.05f64..0.05f64).prop_filter("finite", |v| v.is_finite())
        ) {
            let x = So3::new(a, 0.0, 0.0);
            let y = So3::new(a * 0.5, 0.0, 0.0);
            let xy = bch_third_order::<FiniteF64, So3>(&x, &y);
            let yx = bch_third_order::<FiniteF64, So3>(&y, &x);
            prop_assert!(so3_near(&xy, &yx, 1e-10));
        }
    }
}
