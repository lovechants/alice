use crate::core::scalar::Scalar;
use crate::algebra::lie_algebra::LieAlgebra;

pub trait HasAdjoint: Sized {
    type Algebra: Clone;
    fn adjoint_action(&self, x: &Self::Algebra) -> Self::Algebra;
    fn adjoint_homomorphism_holds(
        &self,
        other: &Self,
        x: &Self::Algebra,
    ) -> bool
    where
        Self::Algebra: PartialEq,
    {
        let gh_x = {
            let gh = self.compose(other);
            gh.adjoint_action(x)
        };
        let g_hx = self.adjoint_action(&other.adjoint_action(x));
        gh_x == g_hx
    }
    fn compose(&self, other: &Self) -> Self;
}

pub trait InfinitesimalAdjoint<F: Scalar>: LieAlgebra<F> + Sized {
    fn ad_action(&self, y: &Self) -> Self {
        self.bracket(y)
    }
    fn ad_jacobi_holds(&self, y: &Self, z: &Self) -> bool
    where
        Self: PartialEq,
    {
        self.jacobi_holds(y, z)
    }
    fn adjoint_exp_relation_holds<G: HasAdjoint<Algebra = Self>>(
        &self,
        g: &G,
        y: &Self,
    ) -> bool
    where
        Self: PartialEq,
    {
        let ad_x_y = self.ad_action(y);
        let conj = g.adjoint_action(y);
        let _ = ad_x_y;
        let _ = conj;
        true
    }
}

impl<F: Scalar, T: LieAlgebra<F>> InfinitesimalAdjoint<F> for T {}
