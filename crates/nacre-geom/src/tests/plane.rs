use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn z0() -> Plane {
    Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap()
}

// --- golden ---

/// **A unit normal handed in is kept bit for bit** — normalizing it again is not the identity.
/// The vector below is itself `normalize` of a random direction; one more `normalize` moves all
/// three components by an ulp, which is why a caller holding the realization it wants kept goes
/// through `from_point_unit_normal`.
#[test]
fn a_unit_normal_handed_in_is_kept_bit_for_bit() {
    let n = Vector3::from_array([
        -0.9792291343701904,
        -0.018250430468946345,
        -0.20193371236201613,
    ]);
    let kept = Plane::from_point_unit_normal(Point3::origin(), n);
    assert_eq!(kept.normal().as_array(), n.as_array());
    let renormalized = Plane::from_point_normal(Point3::origin(), n).unwrap();
    assert_ne!(
        renormalized.normal().as_array(),
        n.as_array(),
        "normalizing this unit vector again is the identity — the fixture shows nothing"
    );
}

/// Mirroring twice about the origin plane returns the plane bit-for-bit: a sign flip is
/// exact, so nothing accumulates.
#[test]
fn mirroring_twice_about_the_origin_is_bit_identical() {
    let p = Plane::through_points(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.5, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 3.25]),
    )
    .unwrap();
    for axis in 0..3 {
        let m = crate::AxisMirror::new(axis, 0.0).unwrap();
        assert_eq!(p.mirrored(m).mirrored(m), p, "axis {axis}");
    }
}

/// An offset mirror places the plane where the geometry says, and its normal flips on the
/// mirrored axis only.
#[test]
fn mirroring_about_an_offset_plane() {
    let p = Plane::from_point_normal(
        Point3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
    )
    .unwrap();
    // x = 1 mirrored about x = 3 lands on x = 5, facing −x.
    let q = p.mirrored(crate::AxisMirror::new(0, 3.0).unwrap());
    assert_eq!(q.origin().as_array(), [5.0, 0.0, 0.0]);
    assert_eq!(q.normal().as_array(), [-1.0, 0.0, 0.0]);
}

#[test]
fn z0_plane_distances() {
    let p = z0();
    assert_eq!(p.signed_distance(Point3::from_array([1.0, 2.0, 0.0])), 0.0);
    assert!(p.contains(Point3::from_array([1.0, 2.0, 0.0]), EPS));
    assert_eq!(p.signed_distance(Point3::from_array([0.0, 0.0, 5.0])), 5.0);
    assert_eq!(
        p.signed_distance(Point3::from_array([0.0, 0.0, -5.0])),
        -5.0
    );
    assert_eq!(p.distance(Point3::from_array([0.0, 0.0, -5.0])), 5.0);
    assert!(!p.contains(Point3::from_array([0.0, 0.0, 5.0]), EPS));
}

#[test]
fn z0_plane_project() {
    let p = z0();
    assert_eq!(
        p.project(Point3::from_array([3.0, 4.0, 5.0])).as_array(),
        [3.0, 4.0, 0.0]
    );
}

#[test]
fn through_axis_points_gives_unit_z_normal() {
    let p = Plane::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([0.0, 1.0, 0.0]),
    )
    .unwrap();
    assert_eq!(p.normal().as_array(), [0.0, 0.0, 1.0]);
    assert_eq!(p.origin().as_array(), [0.0, 0.0, 0.0]);
}

#[test]
fn from_point_normal_normalizes() {
    let p =
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 2.0])).unwrap();
    assert_eq!(p.normal().as_array(), [0.0, 0.0, 1.0]);
}

#[test]
fn degenerate_constructions_return_none() {
    assert!(Plane::from_point_normal(Point3::origin(), Vector3::zero()).is_none());
    // collinear
    assert!(
        Plane::through_points(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 0.0, 0.0]),
        )
        .is_none()
    );
    // coincident
    let q = Point3::from_array([5.0, 6.0, 7.0]);
    assert!(Plane::through_points(q, q, Point3::from_array([1.0, 1.0, 1.0])).is_none());
}

// --- proptest ---

fn pt3() -> impl Strategy<Value = Point3> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point3::from_array)
}

/// A well-conditioned plane plus its three defining points.
fn plane_and_points() -> impl Strategy<Value = (Plane, Point3, Point3, Point3)> {
    (pt3(), pt3(), pt3()).prop_filter_map("degenerate triple", |(a, b, c)| {
        // Reject near-collinear triples so normalize is well-conditioned.
        if (b - a).cross(c - a).norm() >= 1e-3 {
            Plane::through_points(a, b, c).map(|pl| (pl, a, b, c))
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
    fn contains_its_defining_points((pl, a, b, c) in plane_and_points()) {
        let scale = 1e-6 * (a.as_array().iter().chain(&b.as_array()).chain(&c.as_array())
            .map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(pl.distance(a) <= scale);
        prop_assert!(pl.distance(b) <= scale);
        prop_assert!(pl.distance(c) <= scale);
    }

    #[test]
    fn normal_is_unit_and_perpendicular((pl, a, b, c) in plane_and_points()) {
        prop_assert!(approx_eq(pl.normal().norm(), 1.0, EPS, EPS));
        let scale = 1e-6 * ((b - a).norm() + (c - a).norm() + 1.0);
        prop_assert!(pl.normal().dot(b - a).abs() <= scale);
        prop_assert!(pl.normal().dot(c - a).abs() <= scale);
    }

    #[test]
    fn projection_lands_on_plane((pl, ..) in plane_and_points(), p in pt3()) {
        let foot = pl.project(p);
        // foot lies on the plane
        let scale = 1e-6 * (p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!(pl.distance(foot) <= scale);
        // residual (p - foot) is parallel to the normal
        prop_assert!((p - foot).cross(pl.normal()).norm() <= scale);
    }

    #[test]
    fn signed_distance_flips_across_plane((pl, ..) in plane_and_points(), p in pt3()) {
        let q = pl.project(p);            // on the plane
        let n = pl.normal();
        // h = 1.0 dominates the ~1e-10 projection residual, avoiding flakiness.
        prop_assert!(pl.signed_distance(q + n) > 0.0);
        prop_assert!(pl.signed_distance(q - n) < 0.0);
    }
}
