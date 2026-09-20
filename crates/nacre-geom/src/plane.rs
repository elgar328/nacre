//! `Plane` — an unbounded analytic surface: a point and a unit normal.

use nacre_math::{Point3, Vector3};

/// An unbounded plane, stored as an origin point, a **unit** normal, and the
/// **un-normalized** normal it was built from.
///
/// Invariant: `normal` is unit length (to machine precision) — every
/// constructor normalizes and rejects a zero normal, so magnitude consumers
/// (`signed_distance`, `project`, the conditioning gates in [`crate::intersect`])
/// may assume unit length without re-checking. `raw` is a positive multiple of
/// `normal` (the pre-normalization normal), kept because [`Plane::coefficients`]
/// — the handoff to the exact predicates — needs the cross product the plane was
/// built from, not the `sqrt`-rounded unit normal. The predicates are
/// scale-invariant, so `raw`'s length does not affect their sign
/// (`prop_scaling_a_plane_is_invariant`). `origin` is any point on the plane and
/// is not canonicalized.
///
/// ★★ **Keeping `raw` does not make the plane's own points satisfy its form exactly.** That
/// holds when the defining vertices are integers and fails as soon as they are not, because `d`
/// is an `f64` product — see [`Plane::coefficients`]. A caller that needs the two descriptions
/// to be one plane asks [`Plane::spans_exactly`] rather than assuming.
///
/// Minimal by design (M1): no uv-frame / parametric `evaluate(u, v)` yet. A
/// parametric frame (two in-plane basis vectors) arrives in M3, when tess
/// uv-tagging and NURBS need surface parameters.
///
/// `PartialEq` is exact `f64` comparison (including `raw`) — for tests and
/// literal coincidence only. "Is this point on the plane?" goes through
/// [`Plane::distance`] / [`Plane::contains`] with a caller-supplied tolerance,
/// never `==` (overview, principle 2 & 4). `Eq`/`Hash` are deliberately not
/// implemented.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    origin: Point3,
    normal: Vector3,
    raw: Vector3,
}

impl Plane {
    /// From an origin and a normal of any nonzero length. Stores the unit normal
    /// for magnitude use and the given normal verbatim as `raw` for exact
    /// coefficients; returns `None` iff `normal` is the zero vector.
    #[inline]
    pub fn from_point_normal(origin: Point3, normal: Vector3) -> Option<Plane> {
        normal.normalize().map(|unit| Plane {
            origin,
            normal: unit,
            raw: normal,
        })
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

    /// The plane translated by `offset` — the origin shifts, the normal (and its
    /// exact `raw`) are unchanged, so `raw` exactness is preserved (a rigid
    /// translation does not rotate a plane).
    #[inline]
    pub fn translated(self, offset: Vector3) -> Plane {
        Plane {
            origin: self.origin + offset,
            normal: self.normal,
            raw: self.raw,
        }
    }

    /// The plane reflected in `m`.
    ///
    /// Origin **and `raw`** are mirrored, so the exact coefficients survive: rebuilding through
    /// [`Plane::from_point_normal`] would store the *unit* normal as `raw` and lose the exact
    /// (often integer-arithmetic) one. A reflection only flips the sign of one component, so a
    /// mirrored `raw` is as exact as the original — and the unit normal stays unit, so no
    /// re-normalisation is needed either.
    #[inline]
    pub fn mirrored(self, m: crate::AxisMirror) -> Plane {
        Plane {
            origin: m.point(self.origin),
            normal: m.dir(self.normal),
            raw: m.dir(self.raw),
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
    /// ★ It went unnoticed while `Discovered` tolerances were loose enough to absorb it. The one
    /// place it surfaced — a four-plane concurrency whose vertex is genuinely *on* the plane — had
    /// been passing with the residual *exactly equal* to the claimed tolerance, saved only by the
    /// comparison being strict.
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

    /// The coefficients `[a, b, c, d]` of the implicit form `a·X + b·Y + c·Z + d = 0`
    /// — `[raw, −(raw·origin)]`, using the **un-normalized** `raw` normal.
    ///
    /// This is the handoff to `nacre-predicates`: the exact indirect
    /// predicates take plane coefficients as plain arrays, never kernel types, and
    /// are scale-invariant — so `raw`'s length does not affect their sign, only
    /// its exactness matters.
    ///
    /// ★★★ **"Exact" holds of the coefficients, not of any particular point on the plane.**
    /// This doc used to claim that a plane built through exact vertices satisfies the form at
    /// those vertices to *exactly zero*, by exact integer arithmetic. That is true when the
    /// vertices are integers and false as soon as they are not: `d` is the `f64` product
    /// `raw·origin`, and a face at `y = −0.2` with `raw = [0, −3.5, 0]` gets
    /// `d = 0.7000000000000001` — a plane `2⁻⁵⁴` from the one its own points span.
    ///
    /// The claim went unchecked and a caller relied on it, describing one plane by these
    /// coefficients in one question and by its witness triangle in the next; answers composed
    /// across the two were not even transitive. [`contains_exactly`](Plane::contains_exactly) is
    /// that claim made checkable.
    #[inline]
    pub fn coefficients(&self) -> [f64; 4] {
        let [a, b, c] = self.raw.as_array();
        [a, b, c, -self.raw.dot(self.origin - Point3::origin())]
    }

    /// **Does `p` satisfy this plane's [`coefficients`](Plane::coefficients) *exactly*?**
    ///
    /// Not "within a tolerance" and not "in `f64`": the sum is accumulated in exact expansion
    /// arithmetic, so a `true` means the point is on the coefficient plane and nothing else does.
    ///
    /// The question a caller is really asking is *"may I describe this plane by its coefficients
    /// where someone else describes it by these points?"* — and for three non-collinear points
    /// that satisfy the form, the answer is yes, because they span exactly this plane.
    /// [`spans_exactly`](Plane::spans_exactly) is that question in one call.
    pub fn contains_exactly(&self, p: Point3) -> bool {
        nacre_predicates::plane_contains(self.coefficients(), p.as_array())
    }

    /// **Do these three points span exactly this plane?** — the licence to describe one plane two
    /// ways and expect the same answers.
    ///
    /// Both halves are needed. Satisfying the form is not enough on its own: three *collinear*
    /// points satisfy infinitely many planes, so they would license a description that is not this
    /// one. Non-collinearity is therefore tested exactly too, as the direction the three span
    /// being non-zero.
    pub fn spans_exactly(&self, tri: [Point3; 3]) -> bool {
        nacre_predicates::plane_spanned_by(self.coefficients(), tri.map(|p| p.as_array()))
    }

    /// **Is this plane's stored normal parallel to what `tri` spans?** — the weaker licence, for a
    /// caller that reads only the direction.
    ///
    /// ★ **Parallel, not co-directed.** A class's stored normal is allowed to *oppose* its witness
    /// triangle's; `nacre_ops`' `WorkingPlane::frame_sign` records exactly that, and the predicates
    /// that care carry the convention. Demanding agreement of direction here would refuse planes
    /// that agree perfectly about *where* they are.
    ///
    /// This is the half that survives `d` — and `d` is where the two descriptions actually part,
    /// so a normals-only predicate keeps its exact route on a plane the full test rejects.
    pub fn normal_spans(&self, tri: [Point3; 3]) -> bool {
        nacre_predicates::plane_normal_spanned_by(self.coefficients(), tri.map(|p| p.as_array()))
    }
}

#[cfg(test)]
#[path = "tests/plane.rs"]
mod tests;
