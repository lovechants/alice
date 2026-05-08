use crate::core::ops::TopologicalGroup;
use crate::maps::exp_log::HasExpMap;
use crate::topology::manifold::Manifold;

pub trait LieGroup: TopologicalGroup + Manifold + HasExpMap {
    type Scalar;
    fn identity_element() -> Self;
    fn lie_algebra_dim(&self) -> usize {
        self.dim()
    }
}
