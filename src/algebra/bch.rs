use crate::core::scalar::Scalar;
use crate::algebra::lie_algebra::LieAlgebra;

/// Computes Z such that exp(X)exp(Y) = exp(Z) using the
/// Baker-Campbell-Hausdorff formula truncated at third order:
///
/// Z = X + Y
///   + (1/2)[X,Y]
///   + (1/12)[X,[X,Y]]
///   - (1/12)[Y,[X,Y]]
///   - (1/24)[Y,[X,[X,Y]]]         (fourth order, included for accuracy)
///
/// Valid when X and Y are sufficiently small in norm.
pub fn bch_third_order<F, G>(x: &G, y: &G) -> G
where
    F: Scalar,
    G: LieAlgebra<F>,
{
    let half = F::from_f64(0.5).expect("0.5 must be representable");
    let twelfth = F::from_f64(1.0 / 12.0).expect("1/12 must be representable");
    let twenty_fourth = F::from_f64(1.0 / 24.0).expect("1/24 must be representable");
    let xy = x.bracket(y);
    let xxy = x.bracket(&xy);
    let yxy = y.bracket(&xy);
    let yxxy = y.bracket(&xxy);
    x.add(y)
        .add(&xy.scale(&half))
        .add(&xxy.scale(&twelfth))
        .add(&yxy.scale(&twelfth).inverse())
        .add(&yxxy.scale(&twenty_fourth).inverse())
}

/// Computes Z to first order only: Z = X + Y.
/// Valid when [X,Y] = 0 (commuting elements).
pub fn bch_first_order<F, G>(x: &G, y: &G) -> G
where
    F: Scalar,
    G: LieAlgebra<F>,
{
    x.add(y)
}

/// Computes Z to second order: Z = X + Y + (1/2)[X,Y].
pub fn bch_second_order<F, G>(x: &G, y: &G) -> G
where
    F: Scalar,
    G: LieAlgebra<F>,
{
    let half = F::from_f64(0.5).expect("0.5 must be representable");
    let xy = x.bracket(y);
    x.add(y).add(&xy.scale(&half))
}

/// Verifies that BCH is consistent with commutativity:
/// if [X,Y] = 0 then bch(X,Y) = X + Y.
pub fn bch_commuting_check<F, G>(x: &G, y: &G) -> bool
where
    F: Scalar,
    G: LieAlgebra<F>,
{
    if x.bracket(y) != G::zero() {
        return true;
    }
    bch_third_order(x, y) == x.add(y)
}

/// Verifies BCH symmetry under negation:
/// bch(X, -X) = 0 to third order.
pub fn bch_inverse_check<F, G>(x: &G) -> bool
where
    F: Scalar,
    G: LieAlgebra<F>,
{
    let neg_x = x.inverse();
    let result = bch_third_order(x, &neg_x);
    result == G::zero()
}
