//! `Cylinder` — an unbounded analytic surface: a circle swept along an axis.

use crate::Line;
use nacre_math::{Point3, Vector3};

/// An unbounded (infinite) cylinder, stored as its **axis** [`Line`], a **unit**
/// `ref_dir` perpendicular to the axis (the angle-0 direction), and a positive
/// `radius`.
///
/// Invariant: `ref_dir` is unit and orthogonal to `axis.direction()`, and
/// `radius > 0` — the constructor normalizes, orthogonalizes, and rejects a zero
/// axis or non-positive radius. The second in-plane axis (`axis.direction() ×
/// ref_dir`) is derived on demand.
///
/// Parameterized by `(u, v)`: `u` is the angle (radians) from `ref_dir` around
/// the axis, `v` is the signed height along the axis. Storing the axis as a
/// [`Line`] lets [`Cylinder::distance`] reuse [`Line::distance`] (the same
/// cancellation-safe cross form) and [`Cylinder::point_at`] reuse
/// [`Line::point_at`]. This matches the STEP `CYLINDRICAL_SURFACE`, whose
/// `axis2_placement_3d` is exactly `axis.origin()` / `axis.direction()` /
/// `ref_dir`.
///
/// [`normal_at`](Cylinder::normal_at) returns the natural (radial, outward)
/// normal; a face's actual sense is topology's job (as with [`Plane`](crate::Plane)).
///
/// Same `PartialEq` / no-`Eq` / no-`Hash` rationale as [`Line`]: exact `==` is
/// for tests only; on-surface queries use [`Cylinder::distance`] /
/// [`Cylinder::contains`] with a caller-supplied tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cylinder {
    axis: Line,
    ref_dir: Vector3,
    radius: f64,
}

impl Cylinder {
    /// A cylinder about the axis through `origin` along `direction`, with the
    /// angle-0 direction taken from `ref_dir` (its component perpendicular to the
    /// axis, via Gram–Schmidt) and the given `radius`.
    ///
    /// Returns `None` if `direction` is zero, `ref_dir` is parallel to the axis
    /// (no perpendicular component), or `radius` is not positive (rejects `0`,
    /// negatives, and `NaN`). Fallible `from_*` constructor, mirroring
    /// [`Circle::from_center_normal`](crate::Circle).
    pub fn from_axis(
        origin: Point3,
        direction: Vector3,
        ref_dir: Vector3,
        radius: f64,
    ) -> Option<Cylinder> {
        let axis = Line::from_point_direction(origin, direction)?;
        // Strip the axial component so the angle-0 axis is perpendicular to it.
        let d = axis.direction();
        let ref_dir = (ref_dir - d * ref_dir.dot(d)).normalize()?;
        (radius > 0.0).then_some(Cylinder {
            axis,
            ref_dir,
            radius,
        })
    }

    /// The axis line.
    #[inline]
    pub fn axis(&self) -> Line {
        self.axis
    }

    /// The cylinder translated by `offset` — the axis shifts, the direction,
    /// ref_dir, and radius are unchanged.
    #[inline]
    pub fn translated(self, offset: Vector3) -> Cylinder {
        Cylinder {
            axis: self.axis.translated(offset),
            ref_dir: self.ref_dir,
            radius: self.radius,
        }
    }

    /// The unit angle-0 direction (perpendicular to the axis).
    #[inline]
    pub fn ref_dir(&self) -> Vector3 {
        self.ref_dir
    }

    /// The radius.
    #[inline]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The second in-plane axis (angle 90°): `axis.direction() × ref_dir` (unit).
    #[inline]
    fn binormal(self) -> Vector3 {
        self.axis.direction().cross(self.ref_dir)
    }

    /// The point at parameters `(u, v)`:
    /// `axis.point_at(v) + r·(cos u · ref_dir + sin u · (axis × ref_dir))`.
    #[inline]
    pub fn point_at(self, u: f64, v: f64) -> Point3 {
        let (sin, cos) = u.sin_cos();
        self.axis.point_at(v)
            + self.ref_dir * (self.radius * cos)
            + self.binormal() * (self.radius * sin)
    }

    /// The natural outward unit normal at angle `u` (the radial direction, away
    /// from the axis): `cos u · ref_dir + sin u · (axis × ref_dir)`. Independent
    /// of `v`.
    #[inline]
    pub fn normal_at(self, u: f64) -> Vector3 {
        let (sin, cos) = u.sin_cos();
        self.ref_dir * cos + self.binormal() * sin
    }

    /// The natural outward unit normal **at a point** — [`normal_at`](Cylinder::normal_at)'s
    /// sibling for a caller holding a position rather than an angle: the direction from the
    /// axis to `p`, with the axial part removed.
    ///
    /// `None` when `p` sits on the axis, where there is no radial direction to name. `p` is
    /// assumed to lie on the surface; nothing here checks it, because the answer — the radial
    /// direction — is well defined for any point off the axis, and a mesh corner is on the
    /// surface by construction.
    ///
    /// ★ Why a point version exists at all: a mesh corner *does* carry its angle in the
    /// tessellation's provenance tag, but a **seam** vertex is tagged as a model vertex and has
    /// no angle. Reading the position asks no question about how the point was made.
    #[inline]
    pub fn normal_toward(self, p: Point3) -> Option<Vector3> {
        let d = self.axis.direction();
        let from_axis = p - self.axis.origin();
        (from_axis - d * from_axis.dot(d)).normalize()
    }

    /// The circumferential partial derivative `∂P/∂u` at angle `u`:
    /// `r·(−sin u · ref_dir + cos u · (axis × ref_dir))`.
    ///
    /// Magnitude is `radius` (**not** unit) and independent of `v`. The axial
    /// partial `∂P/∂v` is the constant `axis().direction()`, so it has no method.
    /// (Term follows *The NURBS Book*, the M3 reference.)
    #[inline]
    pub fn du(self, u: f64) -> Vector3 {
        let (sin, cos) = u.sin_cos();
        self.ref_dir * (-self.radius * sin) + self.binormal() * (self.radius * cos)
    }

    /// Distance from `p` to the cylinder surface — the tolerance-free residual:
    /// `|perpendicular distance to the axis − r|`, reusing [`Line::distance`].
    #[inline]
    pub fn distance(self, p: Point3) -> f64 {
        (self.axis.distance(p) - self.radius).abs()
    }

    /// What [`distance`](Cylinder::distance) can report for a point exactly on the surface — see
    /// `Plane::distance_eps` for why the term exists at all.
    ///
    /// A cross product, its norm (a `sqrt`), and the subtraction of the radius. Bounding the lot by
    /// `4ε` of the magnitudes that enter — the axis distance and the radius — is loose for a
    /// `sqrt`, which is correctly rounded, and loose is the safe direction here.
    #[inline]
    pub fn distance_eps(self, p: Point3) -> f64 {
        4.0 * f64::EPSILON * (self.axis.distance(p).abs() + self.radius.abs())
    }

    /// Whether `p` lies on the surface within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx_eq;
    use proptest::prelude::*;

    const EPS: f64 = 1e-9;

    fn unit_cylinder() -> Cylinder {
        Cylinder::from_axis(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            1.0,
        )
        .unwrap()
    }

    // --- golden ---

    #[test]
    fn unit_cylinder_basics() {
        let c = unit_cylinder();
        // sin_cos(0) = (0, 1) exactly, so u=0 evaluations are exact.
        assert_eq!(c.point_at(0.0, 0.0).as_array(), [1.0, 0.0, 0.0]);
        assert_eq!(c.point_at(0.0, 3.0).as_array(), [1.0, 0.0, 3.0]);
        assert_eq!(c.normal_at(0.0).as_array(), [1.0, 0.0, 0.0]);
        assert_eq!(c.du(0.0).as_array(), [0.0, 1.0, 0.0]);
        // Distances are exact.
        assert_eq!(c.distance(Point3::from_array([2.0, 0.0, 0.0])), 1.0); // radially outside
        assert_eq!(c.distance(Point3::from_array([0.5, 0.0, 0.0])), 0.5); // inside
        assert_eq!(c.distance(Point3::from_array([0.0, 0.0, 5.0])), 1.0); // on the axis, = radius
        assert_eq!(c.distance(Point3::from_array([1.0, 0.0, 7.0])), 0.0); // on the surface
    }

    #[test]
    fn ref_dir_is_orthogonalized_and_normalized() {
        // ref_dir with an axial component and non-unit length → angle-0 axis +X.
        let c = Cylinder::from_axis(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([2.0, 0.0, 7.0]),
            3.0,
        )
        .unwrap();
        assert_eq!(c.ref_dir().as_array(), [1.0, 0.0, 0.0]);
        assert_eq!(c.point_at(0.0, 0.0).as_array(), [3.0, 0.0, 0.0]);
    }

    #[test]
    fn degenerate_constructions_return_none() {
        let o = Point3::origin();
        let z = Vector3::from_array([0.0, 0.0, 1.0]);
        let x = Vector3::from_array([1.0, 0.0, 0.0]);
        assert!(Cylinder::from_axis(o, Vector3::zero(), x, 1.0).is_none()); // zero axis
        assert!(Cylinder::from_axis(o, z, z, 1.0).is_none()); // ref_dir ∥ axis
        assert!(Cylinder::from_axis(o, z, x, 0.0).is_none()); // zero radius
        assert!(Cylinder::from_axis(o, z, x, -2.0).is_none()); // negative radius
    }

    // --- proptest ---

    fn pt3() -> impl Strategy<Value = Point3> {
        prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point3::from_array)
    }
    fn vec3() -> impl Strategy<Value = Vector3> {
        prop::array::uniform3(-1e6f64..1e6f64).prop_map(Vector3::from_array)
    }

    /// A well-conditioned cylinder: axis direction and ref_dir kept clearly
    /// non-parallel so the in-plane frame is stable (mirrors `circle`).
    fn cylinder() -> impl Strategy<Value = Cylinder> {
        (pt3(), vec3(), vec3(), 0.5f64..10.0).prop_filter_map(
            "degenerate cylinder frame",
            |(origin, direction, ref_dir, radius)| {
                let separated = matches!(
                    (direction.normalize(), ref_dir.normalize()),
                    (Some(a), Some(r)) if a.cross(r).norm() >= 1e-3
                );
                if separated {
                    Cylinder::from_axis(origin, direction, ref_dir, radius)
                } else {
                    None
                }
            },
        )
    }

    proptest! {
        #[test]
        fn point_on_surface_has_zero_distance(c in cylinder(), u in -10.0f64..10.0, v in -1e3f64..1e3) {
            let p = c.point_at(u, v);
            let scale = 1e-6 * (p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
            prop_assert!(c.distance(p) <= scale);
        }

        #[test]
        fn normal_is_unit_and_perpendicular(c in cylinder(), u in -10.0f64..10.0) {
            prop_assert!(approx_eq(c.normal_at(u).norm(), 1.0, EPS, EPS));
            prop_assert!(c.normal_at(u).dot(c.axis().direction()).abs() <= EPS);
            let scale = 1e-6 * (c.radius() + 1.0);
            prop_assert!(c.normal_at(u).dot(c.du(u)).abs() <= scale);
        }

        #[test]
        fn du_norm_is_radius(c in cylinder(), u in -10.0f64..10.0) {
            prop_assert!(approx_eq(c.du(u).norm(), c.radius(), EPS, EPS));
        }

        #[test]
        fn normal_matches_du_cross_axis(c in cylinder(), u in -10.0f64..10.0) {
            // normalize(∂u × axis) is the outward normal — an independent cross-check.
            let n = c.du(u).cross(c.axis().direction()).normalize().unwrap();
            prop_assert!((n - c.normal_at(u)).norm() <= 1e-6);
        }

        #[test]
        fn partials_match_central_differences(c in cylinder(), u in -10.0f64..10.0, v in -10.0f64..10.0) {
            // Design-mandated numerical differentiation. The partials are origin-
            // independent, so evaluate on an origin-centered copy (a large center
            // would just add catastrophic cancellation to the differenced points).
            let c0 = Cylinder::from_axis(
                Point3::origin(),
                c.axis().direction(),
                c.ref_dir(),
                c.radius(),
            )
            .unwrap();
            let h = 1e-6;
            let tol = 1e-6 * (c0.radius() + 1.0);
            let du_num = (c0.point_at(u + h, v) - c0.point_at(u - h, v)) / (2.0 * h);
            prop_assert!((c0.du(u) - du_num).norm() <= tol);
            let dv_num = (c0.point_at(u, v + h) - c0.point_at(u, v - h)) / (2.0 * h);
            prop_assert!((c0.axis().direction() - dv_num).norm() <= tol);
        }
    }
    /// The point form answers what the angle form does, at the point the angle names —
    /// two roads to one direction, which is the only way to know either is right.
    #[test]
    fn the_point_normal_agrees_with_the_angle_normal() {
        let c = unit_cylinder();
        for k in 0..16 {
            let u = std::f64::consts::TAU * f64::from(k) / 16.0;
            for v in [-3.0, 0.0, 7.5] {
                let n = c.normal_toward(c.point_at(u, v)).expect("off the axis");
                assert!((n - c.normal_at(u)).norm() < EPS, "u={u} v={v}");
            }
        }
    }

    /// Radial means: unit, perpendicular to the axis, and pointing away from it. Measured on a
    /// **tilted** cylinder so that no coordinate axis is doing the work by accident.
    #[test]
    fn a_point_normal_is_radial_whatever_the_axis() {
        let c = Cylinder::from_axis(
            Point3::from_array([1.0, -2.0, 0.5]),
            Vector3::from_array([1.0, 2.0, 3.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.5,
        )
        .expect("a tilted cylinder");
        let d = c.axis().direction();
        for k in 0..8 {
            let u = std::f64::consts::TAU * f64::from(k) / 8.0;
            let p = c.point_at(u, 4.0);
            let n = c.normal_toward(p).expect("off the axis");
            assert!((n.norm() - 1.0).abs() < EPS, "unit");
            assert!(n.dot(d).abs() < EPS, "perpendicular to the axis");
            let from_axis = p - c.axis().origin();
            assert!(
                n.dot(from_axis - d * from_axis.dot(d)) > 0.0,
                "away from the axis"
            );
        }
    }

    /// On the axis there is no radial direction, and the answer says so rather than picking one.
    #[test]
    fn a_point_on_the_axis_has_no_radial_direction() {
        let c = unit_cylinder();
        assert!(c.normal_toward(Point3::origin()).is_none());
        assert!(c.normal_toward(c.axis().point_at(9.0)).is_none());
    }
}
