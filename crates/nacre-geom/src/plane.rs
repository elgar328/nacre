//! `Plane` — an unbounded analytic surface: a point and a unit normal.

use nacre_math::{Point3, Vector3};

/// An unbounded plane, stored as an origin point and a **unit** normal.
///
/// Invariant: `normal` is unit length (to machine precision) — every
/// constructor normalizes or is handed a unit vector, and rejects a zero normal, so magnitude
/// consumers (`signed_distance`, `project`, the conditioning gates in [`crate::intersect`])
/// may assume unit length without re-checking. `origin` is any point on the plane and
/// is not canonicalized.
///
/// ★★ **This is a rounded image of a plane, never a description to decide on.** Its own
/// defining points need not satisfy it exactly (`z = 0.1` is not `fl(0.1)`), and two different
/// planes can share one image. Judging reads the plane's name, derived from its defining points
/// without rounding — the exact shortcuts and the plane-class merge alike.
///
/// Minimal by design (M1): no uv-frame / parametric `evaluate(u, v)` yet. A
/// parametric frame (two in-plane basis vectors) arrives in M3, when tess
/// uv-tagging and NURBS need surface parameters.
///
/// `PartialEq` is exact `f64` comparison — for tests and
/// literal coincidence only. "Is this point on the plane?" goes through
/// [`Plane::distance`] / [`Plane::contains`] with a caller-supplied tolerance,
/// never `==`. `Eq`/`Hash` are deliberately not
/// implemented.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    origin: Point3,
    normal: Vector3,
}

impl Plane {
    /// From an origin and a normal of any nonzero length, normalized here; returns `None` iff
    /// `normal` is the zero vector.
    #[inline]
    pub fn from_point_normal(origin: Point3, normal: Vector3) -> Option<Plane> {
        normal
            .normalize()
            .map(|unit| Plane::from_point_unit_normal(origin, unit))
    }

    /// From an origin and a normal that **is already** the unit normal, stored bit for bit.
    ///
    /// Normalizing a unit `f64` vector again is not the identity — the `sqrt` and the division
    /// each round — so a caller holding the realization it wants kept (a correctly rounded
    /// normal, or one read off another plane) hands it here rather than to
    /// [`Plane::from_point_normal`].
    #[inline]
    pub fn from_point_unit_normal(origin: Point3, normal: Vector3) -> Plane {
        debug_assert!(
            (normal.dot(normal) - 1.0).abs() < 1e-12,
            "a unit normal: {normal:?}"
        );
        Plane { origin, normal }
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

    /// The plane translated by `offset` — the origin shifts, the normal is unchanged (a rigid
    /// translation does not rotate a plane).
    #[inline]
    pub fn translated(self, offset: Vector3) -> Plane {
        Plane {
            origin: self.origin + offset,
            normal: self.normal,
        }
    }

    /// The plane reflected in `m`. A reflection only flips the sign of one component, so the
    /// mirrored normal is still unit bit for bit and is not normalized again.
    #[inline]
    pub fn mirrored(self, m: crate::AxisMirror) -> Plane {
        Plane {
            origin: m.point(self.origin),
            normal: m.dir(self.normal),
        }
    }

    /// The same plane facing the other way: the normal negated, which is exact, and the origin
    /// kept.
    #[inline]
    pub fn reversed(self) -> Plane {
        Plane {
            origin: self.origin,
            normal: -self.normal,
        }
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

    /// **What [`distance`](Plane::distance) can report for a point that is exactly on the plane** —
    /// the f64 rounding of its own eight operations, and nothing about where `p` came from.
    ///
    /// A caller asking "is this vertex on this face?" compares a residual against a tolerance, and
    /// the tolerance it has to hand is the *vertex's* — how far the point may sit from where it
    /// should be. That says nothing about the arithmetic performed here, so a point that is exactly
    /// on the plane can still produce a nonzero residual, and the comparison then fails for a
    /// reason neither operand is responsible for. This is that missing term, derived where the
    /// operations are rather than re-spelled at the call site.
    ///
    /// Three differences, three products and two sums, each round-to-nearest at `≤ ε/2` of its own
    /// magnitude; with `|normal| = 1` an over-estimate of the whole is `3ε · Σ|pᵢ − oᵢ|`.
    ///
    /// ★ A loose vertex tolerance absorbs an under-estimate here. The one place it shows — a
    /// four-plane concurrency whose vertex is genuinely *on* the plane — puts the residual *exactly
    /// equal* to the claimed tolerance, saved only by the comparison being strict.
    #[inline]
    pub fn distance_eps(self, p: Point3) -> f64 {
        let d = p - self.origin;
        let spread: f64 = (0..3).map(|i| d[i].abs()).sum();
        3.0 * f64::EPSILON * spread
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
#[path = "tests/plane.rs"]
mod tests;
