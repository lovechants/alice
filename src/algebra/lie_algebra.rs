use crate::core::scalar::Scalar;
use crate::core::module::VectorSpace;
use crate::core::ops::AbelianGroup;

pub trait LieAlgebra<F: Scalar>: VectorSpace<F> {
    fn bracket(&self, other: &Self) -> Self;
    fn antisymmetry_holds(&self, other: &Self) -> bool {
        self.bracket(other) == other.bracket(self).inverse()
    }
    fn jacobi_holds(&self, y: &Self, z: &Self) -> bool {
        let xy_z = self.bracket(y).bracket(z);
        let yz_x = y.bracket(z).bracket(self);
        let zx_y = z.bracket(self).bracket(y);
        xy_z.add(&yz_x).add(&zx_y) == Self::zero()
    }
    fn bracket_bilinear_left(&self, other: &Self, scalar: &F) -> bool {
        let lhs = self.scale(scalar).bracket(other);
        let rhs = self.bracket(other).scale(scalar);
        lhs == rhs
    }
    fn bracket_bilinear_right(&self, other: &Self, scalar: &F) -> bool {
        let lhs = self.bracket(&other.scale(scalar));
        let rhs = self.bracket(other).scale(scalar);
        lhs == rhs
    }
    fn bracket_additive_left(&self, b: &Self, c: &Self) -> bool {
        let lhs = self.add(b).bracket(c);
        let rhs = self.bracket(c).add(&b.bracket(c));
        lhs == rhs
    }
    fn bracket_additive_right(&self, b: &Self, c: &Self) -> bool {
        let lhs = self.bracket(&b.add(c));
        let rhs = self.bracket(b).add(&self.bracket(c));
        lhs == rhs
    }
}
