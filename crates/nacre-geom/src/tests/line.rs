use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn x_axis() -> Line {
    Line::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
    )
    .unwrap()
}

// --- golden ---

#[test]
fn x_axis_basics() {
    let l = x_axis();
    assert_eq!(l.direction().as_array(), [1.0, 0.0, 0.0]);
    assert_eq!(l.point_at(0.0), l.origin()); // exact at t=0
    assert_eq!(l.point_at(3.0).as_array(), [3.0, 0.0, 0.0]);
    assert_eq!(l.distance(Point3::from_array([0.0, 3.0, 0.0])), 3.0);
    assert_eq!(l.distance(Point3::from_array([5.0, 0.0, 0.0])), 0.0);
    assert!(l.contains(Point3::from_array([5.0, 0.0, 0.0]), EPS));
    assert!(!l.contains(Point3::from_array([0.0, 3.0, 0.0]), EPS));
}

#[test]
fn from_point_direction_normalizes() {
    let l =
        Line::from_point_direction(Point3::origin(), Vector3::from_array([2.0, 0.0, 0.0])).unwrap();
    assert_eq!(l.direction().as_array(), [1.0, 0.0, 0.0]);
}

#[test]
fn degenerate_constructions_return_none() {
    let q = Point3::from_array([1.0, 2.0, 3.0]);
    assert!(Line::through_points(q, q).is_none());
    assert!(Line::from_point_direction(q, Vector3::zero()).is_none());
}

// --- proptest ---

fn pt3() -> impl Strategy<Value = Point3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point3::from_array)
}

/// A well-conditioned line plus its two defining points.
fn line_and_points() -> impl Strategy<Value = (Line, Point3, Point3)> {
    (pt3(), pt3()).prop_filter_map("coincident pair", |(a, b)| {
        if (b - a).norm() >= 1e-3 {
            Line::through_points(a, b).map(|l| (l, a, b))
        } else {
            None
        }
    })
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn contains_defining_points_and_unit_direction((l, a, b) in line_and_points()) {
        prop_assert!(approx_eq(l.direction().norm(), 1.0, EPS, EPS));
        let scale = 1e-6 * (a.as_array().iter().chain(&b.as_array())
            .map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(l.distance(a) <= scale);
        prop_assert!(l.distance(b) <= scale);
    }

    #[test]
    fn point_at_zero_is_origin((l, ..) in line_and_points()) {
        prop_assert_eq!(l.point_at(0.0), l.origin());
    }

    #[test]
    fn on_line_points_have_zero_distance((l, ..) in line_and_points(), t in -1e3f64..1e3) {
        let on = l.point_at(t);
        let scale = 1e-6 * (on.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(l.distance(on) <= scale);
    }

    #[test]
    fn distance_matches_projection_form((l, ..) in line_and_points(), p in pt3()) {
        // Cross-check the cross formula against the projection formula.
        let w = p - l.origin();
        let d = l.direction();
        let proj = (w - d * w.dot(d)).norm();
        let scale = 1e-6 * (w.norm() + 1.0);
        prop_assert!((l.distance(p) - proj).abs() <= scale);
    }

    #[test]
    fn distance_invariant_along_line((l, ..) in line_and_points(), p in pt3(), s in -1e3f64..1e3) {
        // Moving p along the line direction doesn't change perpendicular distance.
        let moved = p + s * l.direction();
        let scale = 1e-6 * (p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max)
            + s.abs() + 1.0);
        prop_assert!((l.distance(moved) - l.distance(p)).abs() <= scale);
    }
}
