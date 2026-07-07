//! Exact analytic geometry — the truth layer of the nacre kernel (design.md §3).
//!
//! Surfaces and curves are kept in analytic form forever; meshes are derived
//! (overview 절대원칙 1). Coordinates here are the *cache* side of the
//! truth/cache split — this crate does exact geometry and stores **no**
//! tolerance: every containment query takes the caller's epsilon (overview
//! 절대원칙 4).
//!
//! M1 defined [`Plane`] and [`Line`]; M3 adds [`Circle`] (the full-circle
//! carrier), with `Cylinder`/`Sphere`/`Nurbs` and the `Curve::Intersection`
//! variant to follow. That intersection variant (holding `Handle<Surface>`) is
//! the sole reason geom will later depend on `nacre-store`; today it uses no
//! `Handle` and has no store dependency.

mod circle;
mod line;
mod plane;

pub use circle::Circle;
pub use line::Line;
pub use plane::Plane;

use nacre_math::Point3;

/// A surface — the exact truth of a face's geometry.
///
/// Single-variant for now so topology can hold `Handle<Surface>`; variants grow
/// in M3. **Not `Copy`**: the coming `Nurbs(NurbsSurface)` variant owns
/// heap-allocated control points, so this type is non-`Copy` from the start to
/// match its eventual nature.
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    Plane(Plane),
}

impl Surface {
    /// Unsigned distance from `p` to the surface (the tolerance-free residual).
    ///
    /// Note: closed-form and always finite for analytic surfaces (planes, and
    /// the M6 quadrics). When iterative surfaces (NURBS) arrive in M3 this may
    /// gain fallibility/cost (a `Result` or a separate `try_distance`).
    #[inline]
    pub fn distance(&self, p: Point3) -> f64 {
        match self {
            Surface::Plane(s) => s.distance(p),
        }
    }

    /// Whether `p` lies on the surface within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(&self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

/// A curve — the exact truth of an edge's geometry.
///
/// Single-variant for now so topology can hold `Handle<Curve>`; a `Circle(Circle)`
/// variant (the [`Circle`] carrier already exists), `Nurbs`, and the
/// `Intersection { surfaces: [Handle<Surface>; 2], .. }` variant (design §3)
/// arrive in M3 — that intersection variant is the sole reason geom will later
/// depend on `nacre-store`. **Not `Copy`** (future heap-backed variants), same
/// as [`Surface`].
#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line(Line),
}

impl Curve {
    /// Unsigned distance from `p` to the curve (the tolerance-free residual).
    /// See [`Surface::distance`] for the M3 fallibility note.
    #[inline]
    pub fn distance(&self, p: Point3) -> f64 {
        match self {
            Curve::Line(c) => c.distance(p),
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
}
