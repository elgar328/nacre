//! `Plane` — an unbounded analytic surface: a point and a unit normal.

use nacre_math::{Point3, Vector3};

/// An unbounded plane, stored as an origin point and a **unit** normal.
///
/// Invariant: `normal` is unit length (to machine precision) — every
/// constructor normalizes and rejects a zero normal, so downstream code may
/// assume unit length without re-checking. `origin` is any point on the plane
/// and is not canonicalized (two `Plane`s describing the same geometric plane
/// but with different origins, or opposite normals, are distinct under
/// `PartialEq`).
///
/// Minimal by design (M1): no uv-frame / parametric `evaluate(u, v)` yet. A
/// parametric frame (two in-plane basis vectors) arrives in M3, when tess
/// uv-tagging (design §5) and NURBS need surface parameters.
///
/// `PartialEq` is exact `f64` comparison — for tests and literal coincidence
/// only. "Is this point on the plane?" goes through [`Plane::distance`] /
/// [`Plane::contains`] with a caller-supplied tolerance, never `==` (overview
/// 절대원칙 2 & 4). `Eq`/`Hash` are deliberately not implemented.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    origin: Point3,
    normal: Vector3,
}

impl Plane {
    /// From an origin and a normal of any nonzero length. Normalizes the
    /// normal; returns `None` iff `normal` is the zero vector.
    #[inline]
    pub fn from_point_normal(origin: Point3, normal: Vector3) -> Option<Plane> {
        normal.normalize().map(|normal| Plane { origin, normal })
    }

    /// The plane through three points: `normal = (b − a) × (c − a)` normalized,
    /// `origin = a`. Returns `None` if the points are collinear or coincident
    /// (the cross product is the zero vector).
    ///
    /// The normal follows the right-hand rule of `(a, b, c)`; choosing the
    /// outward orientation of a face is the caller's (topology's) job.
    #[inline]
    pub fn through_points(a: Point3, b: Point3, c: Point3) -> Option<Plane> {
        Self::from_point_normal(a, (b - a).cross(c - a))
    }

    /// A point on the plane (the stored origin).
    #[inline]
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    /// The unit normal.
    #[inline]
    pub fn normal(&self) -> Vector3 {
        self.normal
    }

    /// Signed distance from `p` to the plane, `(p − origin) · normal`: positive
    /// on the normal's side, negative on the other, zero on the plane.
    #[inline]
    pub fn signed_distance(self, p: Point3) -> f64 {
        (p - self.origin).dot(self.normal)
    }

    /// Unsigned distance from `p` to the plane. The tolerance-free residual;
    /// callers apply their own epsilon (see [`Plane::contains`]).
    #[inline]
    pub fn distance(self, p: Point3) -> f64 {
        self.signed_distance(p).abs()
    }

    /// Whether `p` lies on the plane within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }

    /// The closest point on the plane to `p` (`p − signed_distance(p)·normal`).
    #[inline]
    pub fn project(self, p: Point3) -> Point3 {
        p - self.signed_distance(p) * self.normal
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx_eq;
    use proptest::prelude::*;

    const EPS: f64 = 1e-9;

    fn z0() -> Plane {
        Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0])).unwrap()
    }

    // --- golden ---

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
        let p = Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 2.0]))
            .unwrap();
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
}
