//! Exact analytic geometry — the truth layer of the nacre kernel (design.md §3).
//!
//! Surfaces and curves are kept in analytic form forever; meshes are derived
//! (overview 절대원칙 1). Coordinates here are the *cache* side of the
//! truth/cache split — this crate does exact geometry and stores **no**
//! tolerance: every containment query takes the caller's epsilon (overview
//! 절대원칙 4).
//!
//! M1 defined [`Plane`] and [`Line`]; M3 adds [`Circle`] (the full-circle
//! carrier), [`Cylinder`], and [`NurbsCurve`]/[`NurbsSurface`] (rational B-spline
//! evaluation), with `Sphere` and the `Curve::Intersection` variant to follow.
//! ([`NurbsCurve`]/[`NurbsSurface`] are standalone evaluators for now — the
//! `Curve::Nurbs`/`Surface::Nurbs` variants are wired with a producer later.)
//! That intersection variant (holding `Handle<Surface>`) is the sole reason geom
//! will later depend on `nacre-store`; today it uses no `Handle` and has no store
//! dependency.

mod circle;
mod cylinder;
mod line;
mod nurbs;
mod plane;

/// Closed-form intersections (design §3 isolates robustness-sensitive
/// intersection/classification code in one module). Kept as `pub mod` — callers
/// write `intersect::plane_plane(..)`, keeping that isolation visible.
pub mod intersect;

mod region;

pub use circle::Circle;
pub use cylinder::Cylinder;
pub use line::Line;
pub use nurbs::{NurbsCurve, NurbsSurface};
pub use plane::Plane;
pub use region::planar_region_area_centroid;

use nacre_math::{Point3, Vector3};

/// A reflection across the coordinate plane `p[axis] = offset`.
///
/// Axis-aligned only, which is the same restriction rigid motion already has (the kernel rotates
/// about the coordinate axes only): a reflection in a general plane needs an irrational unit
/// normal, so its image could not be reproduced from an exact definition.
///
/// The axis is an index rather than an enum because this crate does not depend on `nacre-scalar`
/// (whose `Axis` names a *rotation* axis); the caller maps its typed axis once, at the boundary.
///
/// **Exactness:** `dir` only flips a sign, so it is always exact. `point` computes `2·offset − p`
/// on one coordinate, which is exact when `offset == 0` and otherwise rounds once — the same
/// rounding a translation by a rational offset already has.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisMirror {
    axis: usize,
    offset: f64,
}

impl AxisMirror {
    /// `None` if `axis` is not 0, 1 or 2.
    #[inline]
    pub fn new(axis: usize, offset: f64) -> Option<AxisMirror> {
        (axis < 3).then_some(AxisMirror { axis, offset })
    }

    /// `p` with its mirrored coordinate replaced by `2·offset − p[axis]`.
    #[inline]
    pub fn point(self, p: Point3) -> Point3 {
        let mut q = p.as_array();
        q[self.axis] = 2.0 * self.offset - q[self.axis];
        Point3::from_array(q)
    }

    /// `v` with its mirrored component negated. A direction is offset-free, so this is a pure
    /// sign flip and stays exact even for a mirror plane away from the origin.
    #[inline]
    pub fn dir(self, v: Vector3) -> Vector3 {
        let mut w = v.as_array();
        w[self.axis] = -w[self.axis];
        Vector3::from_array(w)
    }
}

/// A surface — the exact truth of a face's geometry.
///
/// `Plane` and `Cylinder` are wired; a `Nurbs(NurbsSurface)` variant (the
/// [`NurbsSurface`] evaluator already exists) and `Sphere` arrive when wired with
/// a producer. **Not `Copy`**: the coming `Nurbs` variant owns heap-allocated
/// control points, so this type is non-`Copy` from the start to match its
/// eventual nature.
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    Plane(Plane),
    Cylinder(Cylinder),
}

impl Surface {
    /// Unsigned distance from `p` to the surface (the tolerance-free residual).
    ///
    /// Note: closed-form and always finite for analytic surfaces (planes,
    /// cylinders, and the M6 quadrics). When iterative surfaces (NURBS) arrive
    /// this may gain fallibility/cost (a `Result` or a separate `try_distance`).
    #[inline]
    pub fn distance(&self, p: Point3) -> f64 {
        match self {
            Surface::Plane(s) => s.distance(p),
            Surface::Cylinder(s) => s.distance(p),
        }
    }

    /// **The rounding [`distance`](Surface::distance) itself contributes**, for a caller comparing
    /// that residual against a tolerance that describes only where `p` came from. See
    /// `Plane::distance_eps`.
    #[inline]
    pub fn distance_eps(&self, p: Point3) -> f64 {
        match self {
            Surface::Plane(s) => s.distance_eps(p),
            Surface::Cylinder(s) => s.distance_eps(p),
        }
    }

    /// Whether `p` lies on the surface within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(&self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }

    /// The surface translated by `offset` (a rigid translation).
    #[inline]
    pub fn translated(&self, offset: Vector3) -> Surface {
        match self {
            Surface::Plane(s) => Surface::Plane(s.translated(offset)),
            Surface::Cylinder(s) => Surface::Cylinder(s.translated(offset)),
        }
    }

    /// The surface reflected in `m`, or `None` for a cylinder.
    ///
    /// A reflection reverses the sense of a circle: `R(n × ref) = −(R(n) × R(ref))`, so the image
    /// of a cylinder's parametrisation runs the opposite way round its axis. Which convention a
    /// mirrored quadric should take is decided with the curved-geometry milestone rather than
    /// guessed here, so the function does not exist for that variant instead of existing and
    /// being wrong.
    #[inline]
    pub fn mirrored(&self, m: AxisMirror) -> Option<Surface> {
        match self {
            Surface::Plane(s) => Some(Surface::Plane(s.mirrored(m))),
            Surface::Cylinder(_) => None,
        }
    }
}

/// A curve — the exact truth of an edge's geometry.
///
/// `Line` and `Circle` (the full-circle carrier) are wired; `Nurbs` and the
/// `Intersection { surfaces: [Handle<Surface>; 2], .. }` variant (design §3)
/// arrive later — that intersection variant is the sole reason geom will later
/// depend on `nacre-store`. **Not `Copy`** (future heap-backed variants), same
/// as [`Surface`].
#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line(Line),
    Circle(Circle),
}

impl Curve {
    /// Unsigned distance from `p` to the curve (the tolerance-free residual).
    /// See [`Surface::distance`] for the M3 fallibility note.
    #[inline]
    pub fn distance(&self, p: Point3) -> f64 {
        match self {
            Curve::Line(c) => c.distance(p),
            Curve::Circle(c) => c.distance(p),
        }
    }

    /// The curve translated by `offset` (a rigid translation).
    #[inline]
    pub fn translated(&self, offset: Vector3) -> Curve {
        match self {
            Curve::Line(c) => Curve::Line(c.translated(offset)),
            Curve::Circle(c) => Curve::Circle(c.translated(offset)),
        }
    }

    /// The curve reflected in `m`, or `None` for a circle — see [`Surface::mirrored`] for why the
    /// curved variant is left undefined rather than guessed.
    #[inline]
    pub fn mirrored(&self, m: AxisMirror) -> Option<Curve> {
        match self {
            Curve::Line(c) => Some(Curve::Line(c.mirrored(m))),
            Curve::Circle(_) => None,
        }
    }

    /// Whether `p` lies on the curve within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(&self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

/// Combined relative-or-absolute float comparison for tests.
///
/// Test-only; production geometry never compares coordinates with a bare `==`.
/// (Duplicated from nacre-math, whose copy is `#[cfg(test)] pub(crate)` and thus
/// unreachable here; promote to a shared util only on a third consumer.)
#[cfg(test)]
pub(crate) fn approx_eq(a: f64, b: f64, rel: f64, abs: f64) -> bool {
    let diff = (a - b).abs();
    diff <= abs || diff <= rel * a.abs().max(b.abs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;

    #[test]
    fn surface_delegates_to_plane() {
        let s = Surface::Plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                .unwrap(),
        );
        assert_eq!(s.distance(Point3::from_array([0.0, 0.0, 5.0])), 5.0);
        assert!(s.contains(Point3::from_array([1.0, 2.0, 0.0]), 1e-9));
        assert!(!s.contains(Point3::from_array([0.0, 0.0, 5.0]), 1e-9));
    }

    #[test]
    fn curve_delegates_to_line() {
        let c = Curve::Line(
            Line::through_points(
                Point3::from_array([0.0, 0.0, 0.0]),
                Point3::from_array([1.0, 0.0, 0.0]),
            )
            .unwrap(),
        );
        assert_eq!(c.distance(Point3::from_array([0.0, 3.0, 0.0])), 3.0);
        assert!(c.contains(Point3::from_array([5.0, 0.0, 0.0]), 1e-9));
        assert!(!c.contains(Point3::from_array([0.0, 3.0, 0.0]), 1e-9));
    }

    #[test]
    fn curve_delegates_to_circle() {
        // Unit circle in the XY plane, centre origin.
        let c = Curve::Circle(
            Circle::from_center_normal(
                Point3::origin(),
                Vector3::from_array([0.0, 0.0, 1.0]),
                Vector3::from_array([1.0, 0.0, 0.0]),
                1.0,
            )
            .unwrap(),
        );
        assert_eq!(c.distance(Point3::from_array([2.0, 0.0, 0.0])), 1.0);
        assert!(c.contains(Point3::from_array([1.0, 0.0, 0.0]), 1e-9));
        assert!(!c.contains(Point3::from_array([2.0, 0.0, 0.0]), 1e-9));
    }

    #[test]
    fn surface_delegates_to_cylinder() {
        // Unit-radius cylinder about the Z axis.
        let s = Surface::Cylinder(
            Cylinder::from_axis(
                Point3::origin(),
                Vector3::from_array([0.0, 0.0, 1.0]),
                Vector3::from_array([1.0, 0.0, 0.0]),
                1.0,
            )
            .unwrap(),
        );
        assert_eq!(s.distance(Point3::from_array([2.0, 0.0, 5.0])), 1.0);
        assert!(s.contains(Point3::from_array([1.0, 0.0, 3.0]), 1e-9));
        assert!(!s.contains(Point3::from_array([2.0, 0.0, 5.0]), 1e-9));
    }
}
