use alice::maps::exp_log::{exp_log_roundtrip_error, log_exp_roundtrip_error};

fn mat2(a: f64, b: f64, c: f64, d: f64) -> [f64; 4] {
    [a, b, c, d]
}

fn mat3(data: [f64; 9]) -> [f64; 9] {
    data
}

mod exp_known_values {
    use super::*;

    #[test]
    fn exp_zero_matrix_is_identity_2x2() {
        let z = mat2(0.0, 0.0, 0.0, 0.0);
        let result = alice::maps::faer_bridge_pub::matrix_exp(&z, 2);
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 0.0).abs() < 1e-10);
        assert!((result[2] - 0.0).abs() < 1e-10);
        assert!((result[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn exp_zero_matrix_is_identity_3x3() {
        let z = [0.0f64; 9];
        let result = alice::maps::faer_bridge_pub::matrix_exp(&z, 3);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((result[i * 3 + j] - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn exp_so2_element() {
        let theta = core::f64::consts::PI / 4.0;
        let x = mat2(0.0, -theta, theta, 0.0);
        let result = alice::maps::faer_bridge_pub::matrix_exp(&x, 2);
        let cos = libm::cos(theta);
        let sin = libm::sin(theta);
        assert!((result[0] - cos).abs() < 1e-8, "r[0,0]={} expected={}", result[0], cos);
        assert!((result[1] - (-sin)).abs() < 1e-8, "r[0,1]={} expected={}", result[1], -sin);
        assert!((result[2] - sin).abs() < 1e-8, "r[1,0]={} expected={}", result[2], sin);
        assert!((result[3] - cos).abs() < 1e-8, "r[1,1]={} expected={}", result[3], cos);
    }
    #[test]
    fn exp_so3_rotation_x_axis() {
        let theta = 0.3;
        let x = mat3([
            0.0,   0.0,    0.0,
            0.0,   0.0,    -theta,
            0.0,   theta,  0.0,
        ]);
        let result = alice::maps::faer_bridge_pub::matrix_exp(&x, 3);
        let cos = libm::cos(theta);
        let sin = libm::sin(theta);
        assert!((result[0] - 1.0).abs() < 1e-8);
        assert!((result[4] - cos).abs() < 1e-8);
        assert!((result[5] - (-sin)).abs() < 1e-8);
        assert!((result[7] - sin).abs() < 1e-8);
        assert!((result[8] - cos).abs() < 1e-8);
    }

    #[test]
    fn log_identity_is_zero_2x2() {
        let id = mat2(1.0, 0.0, 0.0, 1.0);
        let result = alice::maps::faer_bridge_pub::matrix_log(&id, 2);
        assert!(result.is_some());
        let log = result.unwrap();
        for v in &log {
            assert!(v.abs() < 1e-10, "expected zero, got {}", v);
        }
    }

    #[test]
    fn log_identity_is_zero_3x3() {
        let id: Vec<f64> = (0..9).map(|i| if i % 4 == 0 { 1.0 } else { 0.0 }).collect();
        let result = alice::maps::faer_bridge_pub::matrix_log(&id, 3);
        assert!(result.is_some());
        for v in result.unwrap() {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn exp_log_roundtrip_zero() {
        let z = [0.0f64; 4];
        let err = exp_log_roundtrip_error(&z, 2);
        assert!(err.is_some());
        assert!(err.unwrap() < 1e-10);
    }

    #[test]
    fn exp_log_roundtrip_so2() {
        let theta = 0.5;
        let x = mat2(0.0, -theta, theta, 0.0);
        let err = exp_log_roundtrip_error(&x, 2);
        assert!(err.is_some());
        assert!(err.unwrap() < 1e-7, "roundtrip error: {}", err.unwrap());
    }

    #[test]
    fn log_exp_roundtrip_so2_rotation() {
        let theta = 0.5;
        let x = mat2(0.0, -theta, theta, 0.0);
        let rot = alice::maps::faer_bridge_pub::matrix_exp(&x, 2);
        let err = log_exp_roundtrip_error(&rot, 2);
        assert!(err.is_some());
        assert!(err.unwrap() < 1e-8, "roundtrip error: {}", err.unwrap());
    }

    #[test]
    fn exp_log_roundtrip_so3_small() {
        let x = mat3([
            0.0,  -0.1,  0.2,
            0.1,   0.0, -0.3,
           -0.2,   0.3,  0.0,
        ]);
        let err = exp_log_roundtrip_error(&x, 3);
        assert!(err.is_some());
        assert!(err.unwrap() < 1e-8, "roundtrip error: {}", err.unwrap());
    }
}

mod exp_properties {
    use super::*;
    use proptest::prelude::*;

    fn skew2_strategy() -> impl Strategy<Value = [f64; 4]> {
        (-0.5f64..0.5f64).prop_map(|t| mat2(0.0, -t, t, 0.0))
    }

    fn skew3_strategy() -> impl Strategy<Value = [f64; 9]> {
        (-0.3f64..0.3f64, -0.3f64..0.3f64, -0.3f64..0.3f64)
            .prop_map(|(a, b, c)| mat3([
                0.0, -c,  b,
                c,   0.0, -a,
               -b,   a,   0.0,
            ]))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(256))]

        #[test]
        fn exp_log_roundtrip_so2_proptest(x in skew2_strategy()) {
            let err = exp_log_roundtrip_error(&x, 2);
            prop_assert!(err.is_some());
            prop_assert!(
                err.unwrap() < 1e-7,
                "roundtrip error too large: {}",
                err.unwrap()
            );
        }

        #[test]
        fn exp_log_roundtrip_so3_proptest(x in skew3_strategy()) {
            let err = exp_log_roundtrip_error(&x, 3);
            prop_assert!(err.is_some());
            prop_assert!(
                err.unwrap() < 1e-7,
                "roundtrip error too large: {}",
                err.unwrap()
            );
        }

        #[test]
        fn exp_so2_preserves_det(t in -1.0f64..1.0f64) {
            let x = mat2(0.0, -t, t, 0.0);
            let r = alice::maps::faer_bridge_pub::matrix_exp(&x, 2);
            let det = r[0] * r[3] - r[1] * r[2];
            prop_assert!(
                libm::fabs(det - 1.0) < 1e-8,
                "det={} expected 1.0",
                det
            );
        }

        #[test]
        fn exp_so3_preserves_det(
            a in -0.3f64..0.3f64,
            b in -0.3f64..0.3f64,
            c in -0.3f64..0.3f64
        ) {
            let x = mat3([0.0, -c, b, c, 0.0, -a, -b, a, 0.0]);
            let r = alice::maps::faer_bridge_pub::matrix_exp(&x, 3);
            let det = r[0]*(r[4]*r[8]-r[5]*r[7])
                    - r[1]*(r[3]*r[8]-r[5]*r[6])
                    + r[2]*(r[3]*r[7]-r[4]*r[6]);
            prop_assert!(
                libm::fabs(det - 1.0) < 1e-7,
                "det={} expected 1.0",
                det
            );
        }

        #[test]
        fn exp_skew3_orthogonal(
            a in -0.3f64..0.3f64,
            b in -0.3f64..0.3f64,
            c in -0.3f64..0.3f64
        ) {
            let x = mat3([0.0, -c, b, c, 0.0, -a, -b, a, 0.0]);
            let r = alice::maps::faer_bridge_pub::matrix_exp(&x, 3);
            for i in 0..3 {
                for j in 0..3 {
                    let dot: f64 = (0..3)
                        .map(|k| r[i * 3 + k] * r[j * 3 + k])
                        .sum();
                    let expected = if i == j { 1.0 } else { 0.0 };
                    prop_assert!(
                        libm::fabs(dot - expected) < 1e-7,
                        "R*R^T[{},{}]={} expected {}",
                        i, j, dot, expected
                    );
                }
            }
        }
    }
}
