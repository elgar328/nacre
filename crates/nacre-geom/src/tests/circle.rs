use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn unit_circle() -> Circle {
    Circle::from_center_normal(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
    )
    .unwrap()
}

// --- golden ---

/// `angle_of` is `point_at`'s inverse, on all four quadrants and on the seam — and it
/// answers for a point off the circle by its in-plane direction (a scaled point maps to the
/// same angle), which is how the arc consumers ask it about realized endpoints.
#[test]
fn angle_of_inverts_point_at() {
    let c = Circle::from_center_normal(
        Point3::from_array([1.0, -2.0, 0.5]),
        Vector3::from_array([0.0, 3.0, 4.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.5,
    )
    .unwrap();
    for k in 0..8 {
        let theta = std::f64::consts::TAU * (k as f64) / 8.0;
        let back = c.angle_of(c.point_at(theta));
        let diff = (back - theta).abs();
        assert!(
            diff < 1e-12 || (std::f64::consts::TAU - diff) < 1e-12,
            "theta {theta} came back as {back}"
        );
    }
    // Off-circle: the direction speaks, not the distance.
    let p = c.center() + (c.point_at(1.0) - c.center()) * 0.25;
    assert!((c.angle_of(p) - 1.0).abs() < 1e-12);
    // The seam itself is angle zero, not 2π.
    assert!(c.angle_of(c.point_at(0.0)) < 1e-15);
}

#[test]
fn unit_circle_basics() {
    let c = unit_circle();
    // sin_cos(0) = (0, 1) exactly, so these are exact.
    assert_eq!(c.point_at(0.0).as_array(), [1.0, 0.0, 0.0]);
    assert_eq!(c.derivative(0.0).as_array(), [0.0, 1.0, 0.0]);
    // Distances are exact (same bits as `f64::sqrt`).
    assert_eq!(c.distance(Point3::origin()), 1.0); // the axis point, = radius
    assert_eq!(c.distance(Point3::from_array([1.0, 0.0, 0.0])), 0.0); // on the circle
    assert_eq!(c.distance(Point3::from_array([2.0, 0.0, 0.0])), 1.0); // radially outside
    assert_eq!(
        c.distance(Point3::from_array([0.0, 0.0, 5.0])),
        26.0_f64.sqrt()
    ); // up the axis
}

#[test]
fn ref_dir_is_orthogonalized_and_normalized() {
    // A ref_dir with a normal component and non-unit length is projected and
    // normalized: normal +Z, ref_dir (2, 0, 7) → angle-0 axis +X.
    let c = Circle::from_center_normal(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([2.0, 0.0, 7.0]),
        3.0,
    )
    .unwrap();
    assert_eq!(c.ref_dir().as_array(), [1.0, 0.0, 0.0]);
    assert_eq!(c.point_at(0.0).as_array(), [3.0, 0.0, 0.0]);
}

#[test]
fn degenerate_constructions_return_none() {
    let c = Point3::origin();
    let z = Vector3::from_array([0.0, 0.0, 1.0]);
    let x = Vector3::from_array([1.0, 0.0, 0.0]);
    assert!(Circle::from_center_normal(c, Vector3::zero(), x, 1.0).is_none()); // zero normal
    assert!(Circle::from_center_normal(c, z, z, 1.0).is_none()); // ref_dir ∥ normal
    assert!(Circle::from_center_normal(c, z, x, 0.0).is_none()); // zero radius
    assert!(Circle::from_center_normal(c, z, x, -2.0).is_none()); // negative radius
}

// --- proptest ---

fn pt3() -> impl Strategy<Value = Point3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point3::from_array)
}
fn vec3() -> impl Strategy<Value = Vector3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Vector3::from_array)
}

/// A well-conditioned circle: normal and ref_dir are kept clearly non-parallel
/// so the in-plane frame is stable (mirrors `plane_and_points`).
fn circle() -> impl Strategy<Value = Circle> {
    (pt3(), vec3(), vec3(), 0.5f64..10.0).prop_filter_map(
        "degenerate circle frame",
        |(center, normal, ref_dir, radius)| {
            let separated = matches!(
                (normal.normalize(), ref_dir.normalize()),
                (Some(n), Some(r)) if n.cross(r).norm() >= 1e-3
            );
            if separated {
                Circle::from_center_normal(center, normal, ref_dir, radius)
            } else {
                None
            }
        },
    )
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn derivative_matches_central_difference(c in circle(), theta in -10.0f64..10.0) {
        // The design-mandated numerical-differentiation check. The derivative is
        // center-independent (impl uses only ref_dir/binormal/radius), so evaluate
        // it on an origin-centered copy — differencing `point_at` at a large center
        // would just add catastrophic cancellation from the huge shared term.
        let c0 =
            Circle::from_center_normal(Point3::origin(), c.normal(), c.ref_dir(), c.radius())
                .unwrap();
        let h = 1e-6;
        let central = (c0.point_at(theta + h) - c0.point_at(theta - h)) / (2.0 * h);
        let tol = 1e-6 * (c0.radius() + 1.0);
        prop_assert!((c0.derivative(theta) - central).norm() <= tol);
    }

    #[test]
    fn point_on_circle_has_zero_distance(c in circle(), theta in -10.0f64..10.0) {
        let p = c.point_at(theta);
        let scale = 1e-6 * (p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(c.distance(p) <= scale);
    }

    #[test]
    fn derivative_norm_is_radius(c in circle(), theta in -10.0f64..10.0) {
        prop_assert!(approx_eq(c.derivative(theta).norm(), c.radius(), EPS, EPS));
    }

    #[test]
    fn derivative_is_perpendicular_to_radius(c in circle(), theta in -10.0f64..10.0) {
        // Velocity ⟂ radius vector — an analytic cross-check independent of the
        // finite difference above.
        let radial = c.point_at(theta) - c.center();
        let scale = 1e-6 * (c.radius() * c.radius() + 1.0);
        prop_assert!(radial.dot(c.derivative(theta)).abs() <= scale);
    }
}
