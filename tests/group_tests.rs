use alice::groups::gl::Gl;
use alice::groups::sl::Sl;
use alice::groups::so::So;
use alice::groups::se::Se;
use alice::groups::aff::Aff;
use alice::groups::matrix_group::GroupError;
use alice::maps::exp_log::HasExpMap;
use alice::topology::manifold::Dynamic;
use alice::core::ops::{Group, Magma};

fn near(a: f64, b: f64, eps: f64) -> bool {
    libm::fabs(a - b) < eps
}

fn near_vec(a: &[f64], b: &[f64], eps: f64) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| near(*x, *y, eps))
}

fn identity2() -> Vec<f64> { vec![1.0, 0.0, 0.0, 1.0] }
fn identity3() -> Vec<f64> { vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0] }

fn rotation2(theta: f64) -> Vec<f64> {
    vec![libm::cos(theta), -libm::sin(theta), libm::sin(theta), libm::cos(theta)]
}

fn rotation3_z(theta: f64) -> Vec<f64> {
    vec![
        libm::cos(theta), -libm::sin(theta), 0.0,
        libm::sin(theta),  libm::cos(theta), 0.0,
        0.0,               0.0,              1.0,
    ]
}

fn se2_homogeneous(theta: f64, tx: f64, ty: f64) -> Vec<f64> {
    vec![
        libm::cos(theta), -libm::sin(theta), tx,
        libm::sin(theta),  libm::cos(theta), ty,
        0.0,               0.0,              1.0,
    ]
}

mod gl_tests {
    use super::*;

    #[test]
    fn gl_identity_construction() {
        let g = Gl::identity(Dynamic::new(2));
        assert!(near_vec(&g.data(), &identity2(), 1e-14));
    }

    #[test]
    fn gl_valid_construction() {
        assert!(Gl::new(vec![2.0, 0.0, 0.0, 3.0], Dynamic::new(2)).is_ok());
    }

    #[test]
    fn gl_rejects_singular() {
        assert_eq!(
            Gl::new(vec![1.0, 0.0, 0.0, 0.0], Dynamic::new(2)),
            Err(GroupError::ConstraintViolated("matrix must be invertible"))
        );
    }

    #[test]
    fn gl_composition() {
        let a = Gl::new(vec![2.0, 0.0, 0.0, 3.0], Dynamic::new(2)).unwrap();
        let b = Gl::new(vec![4.0, 0.0, 0.0, 5.0], Dynamic::new(2)).unwrap();
        let c = a.op(&b);
        assert!(near(c.data()[0], 8.0, 1e-12));
        assert!(near(c.data()[3], 15.0, 1e-12));
    }

    #[test]
    fn gl_inverse() {
        let a = Gl::new(vec![2.0, 0.0, 0.0, 4.0], Dynamic::new(2)).unwrap();
        let prod = a.op(&a.inverse());
        assert!(near_vec(&prod.data(), &identity2(), 1e-10));
    }

    #[test]
    fn gl_exp_log_roundtrip() {
        let alg = alice::groups::gl::GlAlgebra::new(
            vec![0.0, -0.3, 0.3, 0.0], Dynamic::new(2),
        ).unwrap();
        let g = Gl::exp(&alg);
        let recovered = g.log().unwrap();
        assert!(near_vec(&recovered.data(), &alg.data(), 1e-7));
    }
}

mod sl_tests {
    use super::*;

    #[test]
    fn sl_identity_construction() {
        assert!(near_vec(&Sl::identity(Dynamic::new(2)).data(), &identity2(), 1e-14));
    }

    #[test]
    fn sl_valid_construction() {
        assert!(Sl::new(vec![2.0, 0.0, 0.0, 0.5], Dynamic::new(2)).is_ok());
    }

    #[test]
    fn sl_rejects_non_unit_det() {
        assert_eq!(
            Sl::new(vec![2.0, 0.0, 0.0, 2.0], Dynamic::new(2)),
            Err(GroupError::ConstraintViolated("det must equal 1"))
        );
    }

    #[test]
    fn sl_composition_preserves_det() {
        let a = Sl::new(vec![2.0, 1.0, 0.0, 0.5], Dynamic::new(2)).unwrap();
        let b = Sl::new(vec![3.0, 0.0, 1.0, 1.0 / 3.0], Dynamic::new(2)).unwrap();
        let c = a.op(&b);
        let det = alice::groups::matrix_group::MatrixGroup::new(c.data(), 2).unwrap().det();
        assert!(near(det, 1.0, 1e-10));
    }

    #[test]
    fn sl_algebra_rejects_nonzero_trace() {
        assert_eq!(
            alice::groups::sl::SlAlgebra::new(vec![1.0, 0.0, 0.0, 0.0], Dynamic::new(2)),
            Err(GroupError::ConstraintViolated("trace must equal 0"))
        );
    }

    #[test]
    fn sl_algebra_accepts_traceless() {
        assert!(alice::groups::sl::SlAlgebra::new(vec![1.0, 0.0, 0.0, -1.0], Dynamic::new(2)).is_ok());
    }
}

mod so_tests {
    use super::*;

    #[test]
    fn so_identity_construction() {
        assert!(near_vec(&So::identity(Dynamic::new(3)).data(), &identity3(), 1e-14));
    }

    #[test]
    fn so2_rotation_construction() {
        assert!(So::new(rotation2(0.7), Dynamic::new(2)).is_ok());
    }

    #[test]
    fn so_rejects_non_orthogonal() {
        assert!(So::new(vec![2.0, 0.0, 0.0, 1.0], Dynamic::new(2)).is_err());
    }

    #[test]
    fn so_rejects_det_minus_one() {
        assert!(So::new(vec![-1.0, 0.0, 0.0, 1.0], Dynamic::new(2)).is_err());
    }

    #[test]
    fn so_composition_stays_in_group() {
        let a = So::new(rotation2(0.3), Dynamic::new(2)).unwrap();
        let b = So::new(rotation2(0.5), Dynamic::new(2)).unwrap();
        let c = a.op(&b);
        assert!(near_vec(&c.data(), &rotation2(0.8), 1e-10));
    }

    #[test]
    fn so_inverse_is_transpose() {
        let r = So::new(rotation2(1.1), Dynamic::new(2)).unwrap();
        assert!(near_vec(&r.inverse().data(), &rotation2(-1.1), 1e-10));
    }

    #[test]
    fn so_group_axiom_inverse() {
        let r = So::new(rotation2(0.5), Dynamic::new(2)).unwrap();
        assert!(near_vec(&r.op(&r.inverse()).data(), &identity2(), 1e-10));
    }

    #[test]
    fn so3_rotation_z() {
        assert!(So::new(rotation3_z(0.4), Dynamic::new(3)).is_ok());
    }

    #[test]
    fn so_algebra_rejects_non_skew() {
        assert!(alice::groups::so::SoAlgebra::new(vec![1.0, 0.0, 0.0, 1.0], Dynamic::new(2)).is_err());
    }

    #[test]
    fn so_algebra_accepts_skew() {
        assert!(alice::groups::so::SoAlgebra::new(vec![0.0, -0.5, 0.5, 0.0], Dynamic::new(2)).is_ok());
    }

    #[test]
    fn so_exp_log_roundtrip_2d() {
        let alg = alice::groups::so::SoAlgebra::new(
            vec![0.0, -0.4, 0.4, 0.0], Dynamic::new(2),
        ).unwrap();
        let g = So::exp(&alg);
        assert!(near_vec(&g.log().unwrap().data(), &alg.data(), 1e-7));
    }

    #[test]
    fn so_exp_log_roundtrip_3d() {
        let alg = alice::groups::so::SoAlgebra::new(
            vec![0.0,-0.1,0.2, 0.1,0.0,-0.3, -0.2,0.3,0.0], Dynamic::new(3),
        ).unwrap();
        let g = So::exp(&alg);
        assert!(near_vec(&g.log().unwrap().data(), &alg.data(), 1e-7));
    }

    #[test]
    fn so_bracket_e1_e2_gives_e3() {
        let e1 = alice::groups::so::SoAlgebra::new(
            vec![0.0,0.0,0.0, 0.0,0.0,-1.0, 0.0,1.0,0.0], Dynamic::new(3),
        ).unwrap();
        let e2 = alice::groups::so::SoAlgebra::new(
            vec![0.0,0.0,1.0, 0.0,0.0,0.0, -1.0,0.0,0.0], Dynamic::new(3),
        ).unwrap();
        let e3 = alice::groups::so::SoAlgebra::new(
            vec![0.0,-1.0,0.0, 1.0,0.0,0.0, 0.0,0.0,0.0], Dynamic::new(3),
        ).unwrap();
        assert!(near_vec(&e1.bracket(&e2).data(), &e3.data(), 1e-12));
    }

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(256))]

        #[test]
        fn so2_composition_det_preserved(a in -3.14f64..3.14f64, b in -3.14f64..3.14f64) {
            let rc = So::new(rotation2(a), Dynamic::new(2)).unwrap()
                .op(&So::new(rotation2(b), Dynamic::new(2)).unwrap());
            let det = alice::groups::matrix_group::MatrixGroup::new(rc.data(), 2).unwrap().det();
            proptest::prop_assert!(near(det, 1.0, 1e-10));
        }

        #[test]
        fn so2_inverse_gives_identity(a in -3.14f64..3.14f64) {
            let r = So::new(rotation2(a), Dynamic::new(2)).unwrap();
            proptest::prop_assert!(near_vec(&r.op(&r.inverse()).data(), &identity2(), 1e-10));
        }

        #[test]
        fn so3_composition_orthogonal(a in -1.0f64..1.0f64, b in -1.0f64..1.0f64) {
            let rc = So::new(rotation3_z(a), Dynamic::new(3)).unwrap()
                .op(&So::new(rotation3_z(b), Dynamic::new(3)).unwrap());
            let rt = rc.transpose();
            let prod = alice::groups::matrix_group::MatrixGroup::new(rc.data(), 3).unwrap()
                .mul(&alice::groups::matrix_group::MatrixGroup::new(rt.data(), 3).unwrap())
                .unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    let expected = if i == j { 1.0 } else { 0.0 };
                    proptest::prop_assert!(near(prod.get(i, j), expected, 1e-9));
                }
            }
        }
    }
}

mod se_tests {
    use super::*;

    #[test]
    fn se_identity_construction() {
        let expected = vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0];
        assert!(near_vec(&Se::identity(Dynamic::new(2)).to_homogeneous(), &expected, 1e-14));
    }

    #[test]
    fn se_from_homogeneous_valid() {
        assert!(Se::from_homogeneous(se2_homogeneous(0.5, 1.0, 2.0), Dynamic::new(2)).is_ok());
    }

    #[test]
    fn se_from_parts_roundtrip() {
        let theta = 0.6;
        let t = vec![1.5, -0.5];
        let r = So::new(rotation2(theta), Dynamic::new(2)).unwrap();
        let se = Se::from_parts(r, t.clone(), Dynamic::new(2)).unwrap();
        assert!(near_vec(&se.translation(Dynamic::new(2)), &t, 1e-12));
        assert!(near_vec(&se.rotation(Dynamic::new(2)).data(), &rotation2(theta), 1e-10));
    }

    #[test]
    fn se_rejects_invalid_bottom_row() {
        let mut data = se2_homogeneous(0.5, 1.0, 2.0);
        data[6] = 0.1;
        assert!(Se::from_homogeneous(data, Dynamic::new(2)).is_err());
    }

    #[test]
    fn se_composition_rotation_valid() {
        let a = Se::from_homogeneous(se2_homogeneous(0.3, 1.0, 0.0), Dynamic::new(2)).unwrap();
        let b = Se::from_homogeneous(se2_homogeneous(0.2, 0.0, 1.0), Dynamic::new(2)).unwrap();
        let c = a.op(&b);
        let det = alice::groups::matrix_group::MatrixGroup::new(
            c.rotation(Dynamic::new(2)).data(), 2
        ).unwrap().det();
        assert!(near(det, 1.0, 1e-10));
    }

    #[test]
    fn se_inverse() {
        let se = Se::from_homogeneous(se2_homogeneous(0.5, 1.0, 2.0), Dynamic::new(2)).unwrap();
        let prod = se.op(&se.inverse());
        let expected = vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0];
        assert!(near_vec(&prod.to_homogeneous(), &expected, 1e-10));
    }
}

mod aff_tests {
    use super::*;

    #[test]
    fn aff_identity_construction() {
        let expected = vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0];
        assert!(near_vec(&Aff::identity(Dynamic::new(2)).to_homogeneous(), &expected, 1e-14));
    }

    #[test]
    fn aff_from_parts_valid() {
        assert!(Aff::from_parts(vec![2.0,0.0,0.0,3.0], vec![1.0,-1.0], Dynamic::new(2)).is_ok());
    }

    #[test]
    fn aff_rejects_singular_linear() {
        assert!(Aff::from_parts(vec![1.0,0.0,0.0,0.0], vec![0.0,0.0], Dynamic::new(2)).is_err());
    }

    #[test]
    fn aff_composition() {
        let a = Aff::from_parts(vec![2.0,0.0,0.0,1.0], vec![1.0,0.0], Dynamic::new(2)).unwrap();
        let b = Aff::from_parts(vec![1.0,0.0,0.0,3.0], vec![0.0,1.0], Dynamic::new(2)).unwrap();
        let c = a.op(&b);
        let linear = c.linear_part(Dynamic::new(2));
        assert!(near(linear[0], 2.0, 1e-12));
        assert!(near(linear[3], 3.0, 1e-12));
    }

    #[test]
    fn aff_inverse() {
        let a = Aff::from_parts(vec![2.0,0.0,0.0,2.0], vec![1.0,1.0], Dynamic::new(2)).unwrap();
        let prod = a.op(&a.inverse());
        let expected = vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0];
        assert!(near_vec(&prod.to_homogeneous(), &expected, 1e-10));
    }

    #[test]
    fn aff_parts_roundtrip() {
        let linear = vec![3.0, 1.0, 0.0, 2.0];
        let translation = vec![4.0, -2.0];
        let a = Aff::from_parts(linear.clone(), translation.clone(), Dynamic::new(2)).unwrap();
        assert!(near_vec(&a.linear_part(Dynamic::new(2)), &linear, 1e-12));
        assert!(near_vec(&a.translation_part(Dynamic::new(2)), &translation, 1e-12));
    }
}
