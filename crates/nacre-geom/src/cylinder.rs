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

    /// A cylinder whose axis direction and `ref_dir` **are already** the unit, mutually
    /// perpendicular vectors to keep, stored bit for bit — the cache a realization rounded
    /// correctly, which [`Cylinder::from_axis`] would round again (it normalizes both and strips
    /// the axial part of `ref_dir` in `f64`). `None` for a radius that is not positive.
    #[inline]
    pub fn from_unit_frame(
        origin: Point3,
        direction: Vector3,
        ref_dir: Vector3,
        radius: f64,
    ) -> Option<Cylinder> {
        debug_assert!(
            (ref_dir.dot(ref_dir) - 1.0).abs() < 1e-12 && direction.dot(ref_dir).abs() < 1e-12,
            "a unit ref_dir across the axis: {direction:?} {ref_dir:?}"
        );
        (radius > 0.0).then_some(Cylinder {
            axis: Line::from_point_unit_direction(origin, direction),
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
#[path = "tests/cylinder.rs"]
mod tests;
