//! `Line` — an unbounded analytic curve: a point and a unit direction.

use nacre_math::{Point3, Vector3};

/// An unbounded line, stored as an origin point and a **unit** direction.
///
/// Invariant: `direction` is unit length (to machine precision) — constructors
/// normalize and reject a zero direction. `origin` is any point on the line and
/// is not canonicalized. Same `PartialEq` / no-`Eq` / no-`Hash` rationale as
/// [`Plane`](crate::Plane): exact `==` is for tests only; on-line queries use
/// [`Line::distance`] / [`Line::contains`] with a caller-supplied tolerance.
///
/// A `Line` is the *unbounded* curve. An edge trims it to a segment via its
/// endpoint vertices — the bounds live on the edge, not here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line {
    origin: Point3,
    direction: Vector3,
}

impl Line {
    /// From an origin and a direction of any nonzero length. Normalizes;
    /// returns `None` iff `direction` is the zero vector.
    #[inline]
    pub fn from_point_direction(origin: Point3, direction: Vector3) -> Option<Line> {
        direction
            .normalize()
            .map(|direction| Line { origin, direction })
    }

    /// The line through two points: `direction = (b − a)` normalized,
    /// `origin = a`. Returns `None` if `a == b` (coincident).
    #[inline]
    pub fn through_points(a: Point3, b: Point3) -> Option<Line> {
        Self::from_point_direction(a, b - a)
    }

    /// A point on the line (the stored origin).
    #[inline]
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    /// The line translated by `offset` — the origin shifts, the direction is
    /// unchanged (a rigid translation does not rotate a line).
    #[inline]
    pub fn translated(self, offset: Vector3) -> Line {
        Line {
            origin: self.origin + offset,
            direction: self.direction,
        }
    }

    /// The line reflected in `m`.
    ///
    /// The direction is mirrored component-wise rather than rebuilt through
    /// [`Line::from_point_direction`]: a reflection only flips one component's sign, so the
    /// mirrored direction is still exactly unit, and re-normalising an already-unit vector would
    /// only round it.
    #[inline]
    pub fn mirrored(self, m: crate::AxisMirror) -> Line {
        Line {
            origin: m.point(self.origin),
            direction: m.dir(self.direction),
        }
    }

    /// The unit direction.
    #[inline]
    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    /// The point at parameter `t`: `origin + t · direction`. Exact at `t == 0.0`
    /// (returns `origin`); `t` is unclamped.
    #[inline]
    pub fn point_at(self, t: f64) -> Point3 {
        self.origin + t * self.direction
    }

    /// Perpendicular distance from `p` to the line, `‖(p − origin) × direction‖`
    /// (`direction` is unit). The tolerance-free residual (the primitive).
    ///
    /// The cross form avoids the catastrophic cancellation that the projection
    /// form `‖w − (w·d)d‖` suffers when `p` is nearly on the line.
    #[inline]
    pub fn distance(self, p: Point3) -> f64 {
        (p - self.origin).cross(self.direction).norm()
    }

    /// Whether `p` lies on the line within `tol` (a caller-supplied epsilon).
    #[inline]
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        self.distance(p) <= tol
    }
}

#[cfg(test)]
#[path = "tests/line.rs"]
mod tests;
