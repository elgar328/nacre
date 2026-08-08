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
/// scale-invariant, so `raw`'s length does not affect their sign (design §9;
/// `prop_scaling_a_plane_is_invariant`). `origin` is any point on the plane and
/// is not canonicalized.
///
/// ★★ **Keeping `raw` does not make the plane's own points satisfy its form exactly.** That
/// holds when the defining vertices are integers and fails as soon as they are not, because `d`
/// is an `f64` product — see [`Plane::coefficients`]. A caller that needs the two descriptions
/// to be one plane asks [`Plane::spans_exactly`] rather than assuming.
///
/// Minimal by design (M1): no uv-frame / parametric `evaluate(u, v)` yet. A
/// parametric frame (two in-plane basis vectors) arrives in M3, when tess
/// uv-tagging (design §5) and NURBS need surface parameters.
///
/// `PartialEq` is exact `f64` comparison (including `raw`) — for tests and
/// literal coincidence only. "Is this point on the plane?" goes through
/// [`Plane::distance`] / [`Plane::contains`] with a caller-supplied tolerance,
/// never `==` (overview 절대원칙 2 & 4). `Eq`/`Hash` are deliberately not
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
    /// This is the handoff to `nacre-predicates` (design §9): the exact indirect
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
mod tests {
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
        let p = Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 2.0]))
            .unwrap();
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
}
