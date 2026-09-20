use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn z0() -> Plane {
    Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap()
}

// --- golden ---

/// A mirror keeps the exact `raw` coefficients — the reason this lives here rather than in a
/// caller: rebuilding through `from_point_normal` would substitute the rounded unit normal.
/// The plane is built from three integer points so `raw` is an exact cross product that a
/// unit normal cannot represent.
#[test]
fn mirroring_keeps_the_exact_raw_normal() {
    let p = Plane::through_points(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 0.0]),
        Point3::from_array([0.0, 1.0, 3.0]),
    )
    .unwrap();
    let m = crate::AxisMirror::new(0, 0.0).unwrap();
    let q = p.mirrored(m);

    // `raw` is mirrored, not renormalised: exactly the source `raw` with x negated.
    let [rx, ry, rz] = p.raw.as_array();
    assert_eq!(q.raw.as_array(), [-rx, ry, rz]);
    // …and it is *not* the unit normal, which is what the naive rebuild would have stored.
    assert_ne!(q.raw.as_array(), q.normal.as_array());
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
fn coefficients_of_z5_plane() {
    // z = 5: normal +z, origin (0,0,5) ⇒ d = signed_distance(0) = −5.
    let p = Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    assert_eq!(p.coefficients(), [0.0, 0.0, 1.0, -5.0]);
    // A point on the plane evaluates the implicit form to exactly zero.
    let [a, b, c, d] = p.coefficients();
    assert_eq!(a * 0.0 + b * 0.0 + c * 5.0 + d, 0.0);
}

#[test]
fn coefficients_are_exact_on_defining_points() {
    // Un-normalized coefficients from integer points evaluate the implicit form to
    // **exactly** zero on the plane — the sqrt-rounded unit normal could not. The
    // tilted plane through these three is `3x − 6y = 0` (normal (3,−6,0), un-normalized).
    let a = Point3::from_array([0.0, 0.0, 0.0]);
    let b = Point3::from_array([2.0, 1.0, 0.0]);
    let c = Point3::from_array([0.0, 0.0, 3.0]);
    let pl = Plane::through_points(a, b, c).unwrap();
    assert_eq!(pl.coefficients(), [3.0, -6.0, 0.0, 0.0]);
    let [ca, cb, cc, cd] = pl.coefficients();
    let eval = |p: Point3| {
        let [x, y, z] = p.as_array();
        ca * x + cb * y + cc * z + cd
    };
    assert_eq!(eval(a), 0.0);
    assert_eq!(eval(b), 0.0);
    assert_eq!(eval(c), 0.0);
    assert_eq!(eval(Point3::from_array([2.0, 1.0, 5.0])), 0.0); // a 4th exactly-coplanar point
    // `normal()` is still the unit normal, for magnitude consumers.
    assert!((pl.normal().norm() - 1.0).abs() < 1e-15);
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

    /// The implicit form `a·x + b·y + c·z + d` is `signed_distance` scaled by the
    /// un-normalized `raw`'s length: `raw = |raw|·n̂`, so `eval = |raw|·signed_distance`.
    /// Its **sign** matches (scale-invariant), which is all the predicates read.
    #[test]
    fn coefficients_evaluate_to_scaled_signed_distance((pl, ..) in plane_and_points(), p in pt3()) {
        let [a, b, c, d] = pl.coefficients();
        let [x, y, z] = p.as_array();
        let eval = a * x + b * y + c * z + d;
        let raw_len = (a * a + b * b + c * c).sqrt();
        let scale = 1e-9 * raw_len * (p.as_array().iter().map(|v| v.abs()).fold(0.0, f64::max) + 1.0);
        prop_assert!((eval - raw_len * pl.signed_distance(p)).abs() <= scale);
    }
}
