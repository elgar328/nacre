//! Closed-form intersections among analytic surfaces (design §3, §8 M5).
//!
//! Robustness-sensitive intersection code is isolated in this one module
//! (design §3). Everything here is the **construction (coordinate)** side of the
//! precision split (design §3): the returned `Line`/`Point3` are f64 *caches* —
//! the defining surfaces are the truth. So near-degenerate inputs are gated by a
//! conditioning threshold (there is no exact answer to cache), whereas exact
//! **sign** decisions live in `nacre-predicates`.

use crate::{Line, Plane};
use nacre_math::{Point2, Point3, Vector3};

/// The exact ray/segment vs triangle crossing outcomes, re-exported so callers
/// (`nacre-ops`) reach them through geom without depending on `nacre-predicates`
/// directly (the same `topo → geom → predicates` layering as the surface
/// handoffs).
pub use nacre_predicates::{RayCross, SegCross, orient2d};

/// `sin²θ` below which two plane normals count as parallel. Unit normals make
/// `‖n1 × n2‖² = sin²θ ∈ [0, 1]`, so this absolute cutoff is scale-free.
const PARALLEL_EPS: f64 = 1e-16;

/// The line where two planes meet, or `None` if they are parallel (or
/// coincident).
///
/// Closed form — no SSI march or spline cache (design §8 M5). Direction is
/// `n1 × n2`; the base point is the point of the line closest to the origin.
/// Returns a bare [`Line`] (the closed-form truth), not a cached
/// `Curve::Intersection` (that variant, holding `Handle<Surface>`, is for the
/// marched intersections of M7).
///
/// The base point is `p0 = (h1·(n2×d) + h2·(d×n1)) / (d·d)` with `d = n1×n2` and
/// `hi = ni·originᵢ`; one verifies `n1·p0 = h1` and `n2·p0 = h2` (using
/// `n1·(n2×d) = ‖d‖²` and `n1·(d×n1) = 0`), so `p0` lies on both planes, and
/// `p0 ⟂ d` makes it the closest to the origin.
pub fn plane_plane(a: &Plane, b: &Plane) -> Option<Line> {
    let n1 = a.normal();
    let n2 = b.normal();
    let d = n1.cross(n2);
    let dd = d.norm_squared();
    if dd <= PARALLEL_EPS {
        return None; // parallel or coincident — no unique line
    }
    // Plane i is `nᵢ·X = hᵢ` with `hᵢ = nᵢ·originᵢ`.
    let h1 = n1.dot(a.origin() - Point3::origin());
    let h2 = n2.dot(b.origin() - Point3::origin());
    let base = Point3::origin() + (h1 * n2.cross(d) + h2 * d.cross(n1)) / dd;
    // `dd > 0` guarantees `d` is nonzero, so this is always `Some`; keep the
    // `Option` rather than unwrapping (defensive).
    Line::from_point_direction(base, d)
}

/// `|det|` below which three unit normals are too close to a common plane to
/// yield a usable vertex coordinate. `det = n1·(n2×n3) ∈ [−1, 1]` for unit
/// normals, so this is scale-free; and since the coordinate error is roughly
/// `ε_f64 / |det|`, the cutoff bounds that error (~`1e-7` here) rather than
/// being an arbitrary constant.
const COPLANAR_DET_EPS: f64 = 1e-9;

/// The point where three planes meet, or `None` if they do not meet in a
/// well-conditioned point (two parallel, or the three normals near-coplanar).
///
/// Closed-form Cramer: `P = (h1·(n2×n3) + h2·(n3×n1) + h3·(n1×n2)) / det`,
/// `det = n1·(n2×n3)` (the scalar triple product = determinant of rows
/// `[n1, n2, n3]`), `hᵢ = nᵢ·originᵢ`.
///
/// The f64 coordinate is a **cache** — the three planes are the truth (design
/// §3/§4). Its residual to the planes is the `Origin::Discovered` tolerance,
/// which the caller (M5-c) measures when it forms the vertex (closed-form
/// residual is recomputable, unlike an iterative solve's, so it is not returned
/// here). The gate is the conditioning threshold [`COPLANAR_DET_EPS`], not the
/// exact `det3_sign`: an exact-nonzero determinant can still round its f64
/// counterpart to ~0 and blow the coordinate up, and this function's output is
/// an inexact coordinate, so an exact existence test would be the wrong gate.
/// (The exact predicate is used where a *decision* must be exact — M5-c's
/// combinatorial "do these meet".)
pub fn three_planes(a: &Plane, b: &Plane, c: &Plane) -> Option<Point3> {
    let n1 = a.normal();
    let n2 = b.normal();
    let n3 = c.normal();
    let det = n1.dot(n2.cross(n3));
    if det.abs() <= COPLANAR_DET_EPS {
        return None; // parallel or near-coplanar normals — no usable vertex
    }
    let h1 = n1.dot(a.origin() - Point3::origin());
    let h2 = n2.dot(b.origin() - Point3::origin());
    let h3 = n3.dot(c.origin() - Point3::origin());
    let num = h1 * n2.cross(n3) + h2 * n3.cross(n1) + h3 * n1.cross(n2);
    Some(Point3::origin() + num / det)
}

/// The exact sign of `orient3d(V, q, r, s)`, where `V` is the implicit point at
/// which planes `a, b, c` meet and `q, r, s` are explicit points — the
/// geom→predicates handoff (design §9).
///
/// **Sign convention.** `orient3d(a,b,c,d) = det[a−d, b−d, c−d]`, so with `V`
/// first this is `(V − s) · ((q − s) × (r − s))`, and `(q−s)×(r−s) = (r−q)×(s−q)`
/// is the right-hand normal of the triangle `(q, r, s)`. Therefore
///
/// > `+1` ⇔ `V` lies on the **right-hand-normal side** of the triangle `(q, r, s)`,
/// > `-1` on the other side, `0` on its plane.
///
/// Callers that pass a face's outer-CCW triple get "positive = outside" for free,
/// because an outer-CCW loop's RH normal is the outward normal. Do not substitute
/// a `Plane`'s stored normal for the triangle's — the two need not agree in sign.
///
/// The three planes' [`coefficients`](Plane::coefficients) fill the rows of a
/// [`nacre_predicates::ThreePlane`], keeping the predicate dependency inside
/// geom (so `topo → geom → predicates`, and `topo` need not know predicates).
/// `V` is never materialized — the sign is exact even though `V`'s coordinates
/// are generally irrational (design §3, Attene 2020).
///
/// **Precondition (inherited from [`nacre_predicates::indirect_orient3d`]):**
/// `a, b, c` must meet in a single point (`det` of their normals ≠ 0), and
/// `q, r, s` must be non-collinear. A degenerate input yields `0`, which in a
/// release build is indistinguishable from a true coplanar `0` (the debug
/// assertion is compiled out) — callers must pass a valid vertex, as M5-c does.
pub fn three_plane_orient3d(
    a: &Plane,
    b: &Plane,
    c: &Plane,
    q: Point3,
    r: Point3,
    s: Point3,
) -> i8 {
    let tp = nacre_predicates::ThreePlane([a.coefficients(), b.coefficients(), c.coefficients()]);
    nacre_predicates::indirect_orient3d(&tp, q.as_array(), r.as_array(), s.as_array())
}

/// The exact side of triangle `tri`'s plane that the **explicit** point `p` lies on:
/// `+1` on the right-hand-normal side, `-1` on the other, `0` on the plane.
///
/// The all-explicit twin of [`three_plane_orient3d`], and it shares that function's
/// convention exactly — `p` takes `V`'s slot, so `orient3d(p, tri…) = (p − tri[2]) ·
/// (RH normal of tri)`. A face's outward-oriented `tri` therefore reads `+1` for
/// "outside the face's plane".
///
/// This is what decides whether a segment straddles a face's plane, and it is the only
/// place a coordinate enters that decision. For a vertex the operations built, the
/// coordinate *is* the truth; for a `Origin::Discovered` vertex it is a rounded cache of
/// a plane triple, and an exact answer would come from [`three_plane_orient3d`] on that
/// triple instead.
pub fn plane_side(tri: [Point3; 3], p: Point3) -> i8 {
    let d = nacre_predicates::orient3d(
        p.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    );
    match d.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => 1,
        Some(std::cmp::Ordering::Less) => -1,
        _ => 0,
    }
}

/// The exact sign of `a[axis] − b[axis]`, where `a` and `b` are the implicit points at
/// which each plane triple meets — the two-implicit handoff to `nacre-predicates`.
///
/// `0` means the coordinates are exactly equal. Neither point is materialized, and no
/// tolerance is involved: this is the sign of a decision, so it belongs on the predicate
/// side of the precision split (design §3), exactly like [`plane_pair_dir_sign`].
///
/// **Precondition (inherited):** each triple meets in a single point (`det` of its normals
/// ≠ 0), as for [`three_plane_orient3d`].
pub fn three_plane_cmp_coord(a: [&Plane; 3], b: [&Plane; 3], axis: usize) -> i8 {
    let tp = |t: [&Plane; 3]| {
        nacre_predicates::ThreePlane([
            t[0].coefficients(),
            t[1].coefficients(),
            t[2].coefficients(),
        ])
    };
    nacre_predicates::indirect_cmp_coord(&tp(a), &tp(b), axis)
}

/// The exact sign of `(nA × nB) · nC` — how the line `A ∩ B` runs relative to
/// plane `C`'s normal.
///
/// The line is oriented the way [`plane_plane`] orients it (direction `nA × nB`),
/// so a caller may use this to reason about order along that line without ever
/// building it. `0` iff the three normals are coplanar, i.e. the planes have no
/// well-conditioned common point — the same condition that makes
/// [`three_planes`] return `None`, but decided exactly.
///
/// Exact: `det3_sign` of the three normals as rows, since
/// `det[nA; nB; nC] = nA · (nB × nC) = (nA × nB) · nC`. Not an f64 dot product —
/// the sign of a decision, so it belongs to the predicate side of the precision
/// split (design §3).
///
/// **Reads the rows from [`Plane::coefficients`] (the un-normalized `raw`), not [`Plane::normal`].**
/// `raw` is a positive multiple of `normal`, and a determinant's sign is invariant under scaling
/// each row by a positive factor, so this is the *same* function either way — except that `normal`
/// is `raw / ‖raw‖`, and that division rounds each component differently, which can destroy an exact
/// cancellation. Two exactly-parallel walls of a slanted prism have `raw` rows that are exact
/// negatives (`det = 0`), but their unit rows round to `det = ±1`; the guard then admitted a triple
/// its own consumer ([`three_plane_orient3d`], also on `coefficients`) rejects as `D = 0`, aborting.
/// Feeding the guard the same `raw` the consumer reads keeps the two in lockstep. (Measured
/// 2026-07-22; the invariant is already stated on `Plane`: exact predicates take `coefficients`.)
pub fn plane_pair_dir_sign(a: &Plane, b: &Plane, c: &Plane) -> i8 {
    let row = |p: &Plane| {
        let [x, y, z, _] = p.coefficients();
        [x, y, z]
    };
    nacre_predicates::det3_sign([row(a), row(b), row(c)])
}

/// Whether `a` and `b` are the **same plane** — coplanar, exactly.
///
/// The geom→predicates handoff for [`nacre_predicates::planes_coplanar`]: it hands
/// over each plane's exact (un-normalized) [`coefficients`](Plane::coefficients) and
/// asks whether the two `[a, b, c, d]` rows are proportional (rank ≤ 1). Coplanarity
/// is a topological decision, so it lives on the predicate side of the precision
/// split (design §3), like [`plane_side`] and [`plane_pair_dir_sign`] — not a
/// `Plane` method with a length tolerance. Scale-invariant and direction-agnostic:
/// opposite normals still name the same plane.
pub fn planes_coplanar(a: &Plane, b: &Plane) -> bool {
    nacre_predicates::planes_coplanar(a.coefficients(), b.coefficients())
}

/// Exact forward-ray/triangle crossing for kernel types — the geom→predicates
/// handoff for [`nacre_predicates::ray_triangle_cross`] (design §8 M5-d
/// point-in-polyhedron). `d` is the ray direction; `tri` is a single triangle
/// (a caller triangulates a face into these). See [`RayCross`] for the outcome.
pub fn ray_face_cross(p: Point3, d: Vector3, tri: [Point3; 3]) -> RayCross {
    nacre_predicates::ray_triangle_cross(
        p.as_array(),
        d.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    )
}

/// Exact segment/triangle crossing for kernel types — the handoff for
/// [`nacre_predicates::segment_triangle_cross`]. See [`SegCross`].
pub fn segment_face_cross(a: Point3, b: Point3, tri: [Point3; 3]) -> SegCross {
    nacre_predicates::segment_triangle_cross(
        a.as_array(),
        b.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    )
}

/// Where a point sits relative to a closed 2-D ring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RingSide {
    Inside,
    Outside,
    /// On the ring itself — a vertex or a point along an edge.
    OnBoundary,
}

/// Whether `p` is inside `ring`, decided by **exact** signs only (`orient2d`), never by a
/// tolerance.
///
/// Crossing parity along a ray: an edge that straddles `p`'s horizontal line contributes a
/// crossing when `p` lies on the correct side of it, which `orient2d` answers exactly. The
/// half-open straddle test (`>` on one end, `<=` on the other) is what keeps a vertex exactly on
/// the ray from being counted twice.
///
/// The ring is assumed **simple** (no self-intersections); for such a ring parity and winding
/// agree, so this needs no fill rule. A self-intersecting ring is out of contract, and the
/// callers enforce that rather than assume it: `sketch::from_rings` runs
/// [`ring_self_intersection`] over every ring before it classifies containment, and
/// `Profile2d::check` does the same before an operation consumes a profile.
pub fn point_in_ring_2d(p: Point2, ring: &[Point2]) -> RingSide {
    let n = ring.len();
    if n < 3 {
        return RingSide::Outside;
    }
    let pa = p.as_array();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i].as_array(), ring[(i + 1) % n].as_array());
        let side = orient2d(a, b, pa);
        if side == 0.0 && on_segment_2d(a, b, pa) {
            return RingSide::OnBoundary;
        }
        // Half-open in y so a vertex on the ray counts for exactly one of its two edges.
        if (a[1] > pa[1]) != (b[1] > pa[1]) {
            // Upward edge with `p` to its left, or downward edge with `p` to its right.
            let upward = b[1] > a[1];
            if (side > 0.0) == upward {
                inside = !inside;
            }
        }
    }
    if inside {
        RingSide::Inside
    } else {
        RingSide::Outside
    }
}

/// Whether two closed rings meet at all — a proper crossing or a mere touch. Exact.
///
/// A profile's rings must be disjoint (a hole lies strictly inside the outer ring and strictly
/// outside its siblings), so "touching" is a failure just as much as "crossing"; the two are not
/// worth distinguishing here. `O(n·m)`, which is the right complexity for sketch-sized rings.
pub fn rings_cross(a: &[Point2], b: &[Point2]) -> bool {
    let (n, m) = (a.len(), b.len());
    if n < 2 || m < 2 {
        return false;
    }
    for i in 0..n {
        let (p1, p2) = (a[i].as_array(), a[(i + 1) % n].as_array());
        for j in 0..m {
            let (q1, q2) = (b[j].as_array(), b[(j + 1) % m].as_array());
            if segments_meet_2d(p1, p2, q1, q2) {
                return true;
            }
        }
    }
    false
}

/// The first pair of edge indices where a closed ring meets itself, or `None` if the ring is a
/// **simple polygon**. Exact — every decision is an `orient2d` sign.
///
/// [`rings_cross`] cannot answer this for a ring against itself: consecutive edges share an
/// endpoint, so it would report every ring. The two cases are therefore split here.
///
/// - **Non-adjacent edges**: any contact at all is a self-intersection. A mere touch counts —
///   a ring that pinches through a point has no strict inside there, so even-odd is not the
///   author's meaning of it.
/// - **Adjacent edges** share one endpoint legitimately. They fail only when *collinear and
///   overlapping* — a spike that doubles back. Both directions must be tested: the doubling-back
///   edge may be longer than the one it retraces, in which case it is the *first* edge's far
///   endpoint that lies inside the second.
/// - A **zero-length edge** (a repeated consecutive point) is a failure of its own, reported as
///   `(i, i)`. It leaves a degenerate edge in the topology, and it breaks the adjacency rule's
///   premise that neighbours share exactly one point.
///
/// A merely *collinear* vertex (a flat corner in the middle of a straight run) is **not** an
/// error: the neighbours are collinear but do not overlap.
///
/// A ring that passes is a simple polygon, and a simple polygon has nonzero area — which is the
/// unstated precondition of every winding decision taken from a signed area.
///
/// `O(n²)`, the right complexity for sketch-sized rings.
pub fn ring_self_intersection(ring: &[Point2]) -> Option<(usize, usize)> {
    let n = ring.len();
    if n < 3 {
        return None;
    }
    let pt = |i: usize| ring[i % n].as_array();
    for i in 0..n {
        if pt(i) == pt(i + 1) {
            return Some((i, i));
        }
    }
    // Adjacency is cyclic: edge `n-1` ends where edge `0` begins, so that pair shares a point too.
    // Missing this would report every ring as self-intersecting.
    for i in 0..n {
        for j in (i + 1)..n {
            let hit = if j == i + 1 {
                spike(pt(i), pt(i + 1), pt(i + 2))
            } else if i == 0 && j == n - 1 {
                spike(pt(1), pt(0), pt(n - 1))
            } else {
                segments_meet_2d(pt(i), pt(i + 1), pt(j), pt(j + 1))
            };
            if hit {
                return Some((i, j));
            }
        }
    }
    None
}

/// Whether the edges `[u, s]` and `[s, v]`, which share `s`, double back over one another.
fn spike(u: [f64; 2], s: [f64; 2], v: [f64; 2]) -> bool {
    orient2d(u, s, v) == 0.0 && (on_segment_2d(u, s, v) || on_segment_2d(s, v, u))
}

/// Exact segment/segment test: `true` for a proper crossing **or** any touching contact.
fn segments_meet_2d(p1: [f64; 2], p2: [f64; 2], q1: [f64; 2], q2: [f64; 2]) -> bool {
    let d1 = orient2d(q1, q2, p1);
    let d2 = orient2d(q1, q2, p2);
    let d3 = orient2d(p1, p2, q1);
    let d4 = orient2d(p1, p2, q2);
    // Proper crossing: each segment separates the other's endpoints.
    if ((d1 > 0.0) != (d2 > 0.0))
        && ((d3 > 0.0) != (d4 > 0.0))
        && d1 != 0.0
        && d2 != 0.0
        && d3 != 0.0
        && d4 != 0.0
    {
        return true;
    }
    // Collinear or endpoint contact.
    (d1 == 0.0 && on_segment_2d(q1, q2, p1))
        || (d2 == 0.0 && on_segment_2d(q1, q2, p2))
        || (d3 == 0.0 && on_segment_2d(p1, p2, q1))
        || (d4 == 0.0 && on_segment_2d(p1, p2, q2))
}

/// Whether the **collinear** point `p` lies within segment `ab`'s extent (callers check
/// collinearity with `orient2d` first, so this is a bounding-box question only).
fn on_segment_2d(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    p[0] >= a[0].min(b[0])
        && p[0] <= a[0].max(b[0])
        && p[1] >= a[1].min(b[1])
        && p[1] <= a[1].max(b[1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;
    use proptest::prelude::*;

    fn plane(origin: [f64; 3], normal: [f64; 3]) -> Plane {
        Plane::from_point_normal(Point3::from_array(origin), Vector3::from_array(normal)).unwrap()
    }

    // --- golden ---

    #[test]
    fn plane_plane_axis_planes_give_z_axis() {
        // x = 0 (normal +x) ∩ y = 0 (normal +y) = the z-axis.
        let line = plane_plane(
            &plane([0.0; 3], [1.0, 0.0, 0.0]),
            &plane([0.0; 3], [0.0, 1.0, 0.0]),
        )
        .unwrap();
        let dir = line.direction().as_array();
        assert!(dir[0].abs() < 1e-15 && dir[1].abs() < 1e-15 && dir[2].abs() > 1.0 - 1e-15);
        let o = line.origin().as_array();
        assert!(o[0].abs() < 1e-15 && o[1].abs() < 1e-15);
    }

    #[test]
    fn plane_plane_parallel_is_none() {
        // Same normal, different offset — parallel, never meet.
        assert!(
            plane_plane(
                &plane([0.0; 3], [0.0, 0.0, 1.0]),
                &plane([0.0, 0.0, 3.0], [0.0, 0.0, 1.0])
            )
            .is_none()
        );
        // Anti-parallel normals are parallel too.
        assert!(
            plane_plane(
                &plane([0.0; 3], [0.0, 0.0, 1.0]),
                &plane([0.0, 0.0, 3.0], [0.0, 0.0, -1.0])
            )
            .is_none()
        );
    }

    // --- proptest ---

    fn coord() -> impl Strategy<Value = f64> {
        -100.0f64..100.0
    }

    fn vec3() -> impl Strategy<Value = Vector3> {
        prop::array::uniform3(-1.0f64..1.0).prop_map(Vector3::from_array)
    }

    /// Two planes whose normals are well-separated (so their line is
    /// well-conditioned).
    fn two_planes() -> impl Strategy<Value = (Plane, Plane)> {
        (
            prop::array::uniform3(coord()),
            vec3(),
            prop::array::uniform3(coord()),
            vec3(),
        )
            .prop_filter_map("zero/near-parallel normals", |(o1, n1, o2, n2)| {
                let a = Plane::from_point_normal(Point3::from_array(o1), n1)?;
                let b = Plane::from_point_normal(Point3::from_array(o2), n2)?;
                (a.normal().cross(b.normal()).norm() >= 0.1).then_some((a, b))
            })
    }

    proptest! {
        /// The returned line lies on both planes and runs perpendicular to both
        /// normals.
        #[test]
        fn plane_plane_line_lies_on_both((a, b) in two_planes(), t in -100.0f64..100.0) {
            let line = plane_plane(&a, &b).unwrap();
            let dir = line.direction();
            prop_assert!(a.normal().dot(dir).abs() <= 1e-9);
            prop_assert!(b.normal().dot(dir).abs() <= 1e-9);
            let p = line.point_at(t);
            let mag = p.as_array().iter().map(|v| v.abs()).fold(0.0, f64::max);
            let scale = 1e-9 * (mag + 1.0);
            prop_assert!(a.distance(p) <= scale);
            prop_assert!(b.distance(p) <= scale);
        }
    }

    // --- three_planes / three_plane_orient3d ---

    /// The sign of an f64 (`+1`/`-1`/`0`) — not `f64::signum`, which maps `0.0`
    /// to `+1.0`; a coplanar `orient3d` (exactly `0.0`) must read as `0`.
    fn sign_f64(x: f64) -> i8 {
        if x > 0.0 {
            1
        } else if x < 0.0 {
            -1
        } else {
            0
        }
    }

    #[test]
    fn three_planes_axis_gives_unit_vertex() {
        // x = 1, y = 1, z = 1 ⇒ (1, 1, 1).
        let v = three_planes(
            &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            &plane([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
            &plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
        )
        .unwrap();
        assert_eq!(v.as_array(), [1.0, 1.0, 1.0]);
    }

    #[test]
    fn three_planes_parallel_pair_is_none() {
        // x = 0 and x = 1 are parallel ⇒ no vertex, whatever the third plane.
        assert!(
            three_planes(
                &plane([0.0; 3], [1.0, 0.0, 0.0]),
                &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
                &plane([0.0; 3], [0.0, 1.0, 0.0]),
            )
            .is_none()
        );
    }

    /// Well-conditioned unit normals: three vectors whose triple product is well
    /// clear of zero (so the vertex is well-conditioned).
    fn three_unit_normals() -> impl Strategy<Value = (Vector3, Vector3, Vector3)> {
        (vec3(), vec3(), vec3()).prop_filter_map("zero/near-coplanar normals", |(a, b, c)| {
            let n1 = a.normalize()?;
            let n2 = b.normalize()?;
            let n3 = c.normalize()?;
            (n1.dot(n2.cross(n3)).abs() >= 0.1).then_some((n1, n2, n3))
        })
    }

    fn ivec3() -> impl Strategy<Value = [i64; 3]> {
        prop::array::uniform3(-30i64..=30)
    }

    /// Pinned against `plane_plane` itself: the sign must agree with the dot of
    /// that function's actual line direction and `c`'s normal. The seam ordering
    /// in `nacre-ops` assumes exactly this coupling.
    #[test]
    fn three_plane_cmp_coord_orders_two_meets() {
        // The unit cube's corners `(0,0,0)` and `(1,1,0)`, each as a triple of its faces.
        let px = |d: f64| {
            Plane::from_point_normal(
                Point3::from_array([d, 0.0, 0.0]),
                Vector3::from_array([1.0, 0.0, 0.0]),
            )
            .unwrap()
        };
        let py = |d: f64| {
            Plane::from_point_normal(
                Point3::from_array([0.0, d, 0.0]),
                Vector3::from_array([0.0, 1.0, 0.0]),
            )
            .unwrap()
        };
        let pz = Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
            .unwrap();
        let (x0, x1, y0, y1) = (px(0.0), px(1.0), py(0.0), py(1.0));
        let a = [&x0, &y0, &pz];
        let b = [&x1, &y1, &pz];
        assert_eq!(three_plane_cmp_coord(a, b, 0), -1);
        assert_eq!(three_plane_cmp_coord(b, a, 0), 1);
        assert_eq!(three_plane_cmp_coord(a, b, 1), -1);
        assert_eq!(three_plane_cmp_coord(a, b, 2), 0); // both on z = 0

        // A `Plane`'s normal is unit, so the coefficients carry a `d` that is not an
        // integer here — the predicate is scale-invariant and never divides.
        let tilt = Plane::through_points(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([0.0, 0.0, 1.0]),
        )
        .unwrap();
        let c = [&x0, &y0, &tilt];
        assert_eq!(three_plane_cmp_coord(c, c, 2), 0);
        assert_eq!(three_plane_cmp_coord(a, c, 2), -1); // (0,0,0) below (0,0,1)
    }

    #[test]
    fn plane_pair_dir_sign_agrees_with_plane_plane() {
        let cases = [
            // x=0 ∩ y=0 is the z axis, direction (1,0,0)×(0,1,0) = +z.
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 1),
            // Flip c's normal ⇒ flip the sign.
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0], -1),
            // Swap a and b ⇒ the line reverses ⇒ flip the sign.
            ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], -1),
            // c parallel to the line ⇒ normals coplanar ⇒ 0.
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], 0),
        ];
        for (na, nb, nc, expect) in cases {
            let (a, b, c) = (
                plane([0.0; 3], na),
                plane([0.0; 3], nb),
                plane([0.0; 3], nc),
            );
            assert_eq!(
                plane_pair_dir_sign(&a, &b, &c),
                expect,
                "{na:?} {nb:?} {nc:?}"
            );
            let l = plane_plane(&a, &b).expect("a and b are not parallel here");
            let dot = l.direction().dot(c.normal());
            let dot_sign = if dot > 1e-12 {
                1
            } else if dot < -1e-12 {
                -1
            } else {
                0
            };
            assert_eq!(dot_sign, expect, "plane_plane's direction disagrees");
        }
    }

    /// The two anti-parallel side walls of a `(1,1,1)`-slanted prism, as they actually came off
    /// `build_prism` (measured 2026-07-22). Their `raw` rows sum to an exact zero in the 2×2 minors,
    /// so `det = 0` against any third plane — but the `sqrt`-rounded **unit** rows do not, and the
    /// old `normal()`-based predicate returned `±1`, admitting a triple its consumer rejects as
    /// `D = 0`. Pinned with the real coefficients so a regression to `normal()` fails with the
    /// reason attached. `raw` is set verbatim via `from_point_normal` (which stores its argument as
    /// `raw`); the origin is irrelevant here (the predicate reads only the normal rows).
    #[test]
    fn parallel_planes_stay_degenerate_after_normalization() {
        let p = |raw: [f64; 3]| {
            Plane::from_point_normal(Point3::origin(), Vector3::from_array(raw)).unwrap()
        };
        let a = p([
            -2.220446049250313e-16,
            -2.828427124746191,
            2.8284271247461907,
        ]);
        let b = p([0.0, 2.8284271247461907, -2.8284271247461907]);
        // A third, independent plane (the prism's tilted cap direction).
        let c = p([
            -0.5773502691896258,
            -0.5773502691896258,
            -0.5773502691896258,
        ]);
        assert_eq!(
            plane_pair_dir_sign(&a, &b, &c),
            0,
            "anti-parallel walls have no well-conditioned common line — det is exactly 0"
        );
        // The bug this pins: the unit-normal determinant does NOT vanish.
        let unit_det = nacre_predicates::det3_sign([
            a.normal().as_array(),
            b.normal().as_array(),
            c.normal().as_array(),
        ]);
        assert_ne!(
            unit_det, 0,
            "sanity: normalized normals round the exact zero away — the reason we read coefficients"
        );
    }

    #[test]
    fn planes_coplanar_names_the_same_plane_regardless_of_scale_or_direction() {
        // z = 0, built three ways: two positive normals of different magnitude and
        // one opposite normal. All name the same plane.
        let a = Plane::through_points(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        )
        .unwrap();
        let bigger = Plane::through_points(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([0.0, 3.0, 0.0]),
        )
        .unwrap();
        let opposite = Plane::through_points(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
        )
        .unwrap();
        assert!(planes_coplanar(&a, &bigger));
        assert!(planes_coplanar(&a, &opposite));
        // Parallel but offset (z = 1) and non-parallel (y = 0): distinct planes.
        assert!(!planes_coplanar(
            &a,
            &plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0])
        ));
        assert!(!planes_coplanar(&a, &plane([0.0; 3], [0.0, 1.0, 0.0])));
    }

    /// `+1` means `V` is on the triangle's right-hand-normal side — pinned here
    /// because the whole seam-ordering algebra in `nacre-ops` hangs off it, and
    /// the doc comment said the opposite until this test existed.
    #[test]
    fn three_plane_orient3d_is_positive_on_the_rh_normal_side() {
        let at = |n: [f64; 3], d: f64| {
            Plane::from_point_normal(
                Point3::from_array([n[0] * d, n[1] * d, n[2] * d]),
                Vector3::from_array(n),
            )
            .unwrap()
        };
        // x=1, y=1, z=1 meet at V = (1,1,1).
        let (px, py, pz) = (
            at([1.0, 0.0, 0.0], 1.0),
            at([0.0, 1.0, 0.0], 1.0),
            at([0.0, 0.0, 1.0], 1.0),
        );
        // Triangle in the z=0 plane, CCW seen from +z ⇒ RH normal is +z. V is above it.
        let (q, r, s) = (
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        );
        assert_eq!(three_plane_orient3d(&px, &py, &pz, q, r, s), 1);
        // Swapping two triangle points flips its RH normal, hence the sign.
        assert_eq!(three_plane_orient3d(&px, &py, &pz, r, q, s), -1);
        // A plane's *stored* normal is not the triangle's RH normal: `pz` above has
        // normal +z, but the triangle (q, s, r) spans the same plane with RH normal
        // −z. Reading the convention off `Plane::normal()` would invert the answer.
        assert_eq!(three_plane_orient3d(&px, &py, &pz, q, s, r), -1);
    }

    /// `plane_side` is the all-explicit twin, and shares the convention exactly: put the
    /// implicit point's coordinates in and the two agree, sign for sign. Pinning that
    /// here is what lets a caller mix them without thinking.
    #[test]
    fn plane_side_shares_three_plane_orient3d_s_convention() {
        // The same triangle in `z = 0`, RH normal `+z`.
        let (q, r, s) = (
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        );
        let v = Point3::from_array([1.0, 1.0, 1.0]); // above it
        assert_eq!(plane_side([q, r, s], v), 1);
        assert_eq!(plane_side([r, q, s], v), -1); // flip the triangle, flip the sign
        assert_eq!(
            plane_side([q, r, s], Point3::from_array([1.0, 1.0, -1.0])),
            -1
        );
        // On the plane is exactly zero, however far from the triangle itself.
        assert_eq!(
            plane_side([q, r, s], Point3::from_array([9.0, -4.0, 0.0])),
            0
        );
        assert_eq!(plane_side([q, r, s], q), 0);
    }

    proptest! {
        /// The two orient3d handoffs are one predicate seen from two sides: a plane
        /// triple's meet, fed to `plane_side` as coordinates, gives the same sign the
        /// implicit form gives without ever building it.
        #[test]
        fn prop_plane_side_agrees_with_three_plane_orient3d(
            v in prop::array::uniform3(-20.0f64..20.0),
            t in prop::array::uniform3(prop::array::uniform3(-20.0f64..20.0)),
        ) {
            let tri = t.map(Point3::from_array);
            let e1 = tri[1] - tri[0];
            let e2 = tri[2] - tri[0];
            prop_assume!(e1.cross(e2).norm() > 1e-6);
            // Three axis planes meeting exactly at `v`.
            let at = |n: [f64; 3]| {
                Plane::from_point_normal(Point3::from_array(v), Vector3::from_array(n)).unwrap()
            };
            let (px, py, pz) = (
                at([1.0, 0.0, 0.0]),
                at([0.0, 1.0, 0.0]),
                at([0.0, 0.0, 1.0]),
            );
            prop_assert_eq!(
                plane_side(tri, Point3::from_array(v)),
                three_plane_orient3d(&px, &py, &pz, tri[0], tri[1], tri[2])
            );
        }
    }

    proptest! {
        /// Three well-conditioned planes through a target point recover it.
        #[test]
        fn three_planes_recovers_constructed_vertex(
            p in prop::array::uniform3(-100.0f64..100.0),
            (n1, n2, n3) in three_unit_normals(),
        ) {
            let p = Point3::from_array(p);
            let v = three_planes(
                &Plane::from_point_normal(p, n1).unwrap(),
                &Plane::from_point_normal(p, n2).unwrap(),
                &Plane::from_point_normal(p, n3).unwrap(),
            )
            .unwrap();
            let mag = p.as_array().iter().map(|x| x.abs()).fold(0.0, f64::max);
            prop_assert!(v.distance(p) <= 1e-6 * (mag + 1.0));
        }

        /// End-to-end wiring: `three_plane_orient3d` (via `coefficients()`) agrees
        /// with `orient3d` evaluated at the constructed vertex. Because `Plane`
        /// normalizes, the implicit vertex is `p + O(1e-16)`, not exactly `p`; so
        /// the config is restricted to well-conditioned normals and a non-coplanar
        /// (p, q, r, s) — there the true `orient3d` is a nonzero integer, which the
        /// tiny vertex perturbation cannot flip. This exercises the new geom code
        /// (coefficient extraction + assembly), not `indirect_orient3d` itself.
        #[test]
        fn three_plane_orient3d_matches_materialized(
            p in ivec3(),
            normals in prop::array::uniform3(ivec3()),
            q in ivec3(),
            r in ivec3(),
            s in ivec3(),
        ) {
            let pf = Point3::from_array(p.map(|v| v as f64));
            let planes: Option<Vec<Plane>> = normals
                .iter()
                .map(|n| {
                    Plane::from_point_normal(pf, Vector3::from_array(n.map(|v| v as f64)))
                })
                .collect();
            prop_assume!(planes.is_some()); // reject a zero integer normal
            let planes = planes.unwrap();
            let u: Vec<_> = planes.iter().map(|pl| pl.normal()).collect();
            prop_assume!(u[0].dot(u[1].cross(u[2])).abs() >= 0.1); // well-conditioned

            let qf = Point3::from_array(q.map(|v| v as f64));
            let rf = Point3::from_array(r.map(|v| v as f64));
            let sf = Point3::from_array(s.map(|v| v as f64));
            let expected = sign_f64(nacre_predicates::orient3d(
                pf.as_array(),
                qf.as_array(),
                rf.as_array(),
                sf.as_array(),
            ));
            prop_assume!(expected != 0); // coplanar: perturbation could flip it

            prop_assert_eq!(
                three_plane_orient3d(&planes[0], &planes[1], &planes[2], qf, rf, sf),
                expected
            );
        }
    }

    #[test]
    fn ray_and_segment_face_cross_handoff() {
        // CCW triangle in z = 0, right-hand normal +z.
        let tri = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        ];
        // Ray from below, straight up through the interior ⇒ forward Cross(+1).
        assert_eq!(
            ray_face_cross(
                Point3::from_array([0.25, 0.25, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                tri,
            ),
            RayCross::Cross(1)
        );
        // Ray pointing away from the triangle ⇒ Miss.
        assert_eq!(
            ray_face_cross(
                Point3::from_array([0.25, 0.25, 1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                tri,
            ),
            RayCross::Miss
        );
        // Segment straddling the plane through the interior ⇒ Cross(+1).
        assert_eq!(
            segment_face_cross(
                Point3::from_array([0.25, 0.25, -1.0]),
                Point3::from_array([0.25, 0.25, 1.0]),
                tri,
            ),
            SegCross::Cross(1)
        );
    }

    // --- ring_self_intersection ---

    fn ring(pts: &[[f64; 2]]) -> Vec<Point2> {
        pts.iter().map(|p| Point2::from_array(*p)).collect()
    }

    /// The net that catches the likeliest way to get this wrong: edges `0` and `n-1` share a point
    /// too, so a non-cyclic adjacency test reports *every* ring. Convex, non-convex, and a flat
    /// (collinear) corner all have to pass.
    #[test]
    fn simple_polygons_have_no_self_intersection() {
        let square = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
        assert_eq!(ring_self_intersection(&square), None);
        let triangle = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
        assert_eq!(ring_self_intersection(&triangle), None);
        // Reflex corner (an L).
        let l = ring(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 2.0],
            [2.0, 2.0],
            [2.0, 4.0],
            [0.0, 4.0],
        ]);
        assert_eq!(ring_self_intersection(&l), None);
        // A flat corner: `(2,0)` sits mid-run on a straight edge. Collinear but not overlapping.
        let flat = ring(&[[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]);
        assert_eq!(ring_self_intersection(&flat), None);
    }

    #[test]
    fn a_bowtie_crosses_itself() {
        let bowtie = ring(&[[0.0, 0.0], [4.0, 4.0], [4.0, 0.0], [0.0, 4.0]]);
        // Edges 0 (`(0,0)→(4,4)`) and 2 (`(4,0)→(0,4)`) are the crossing pair.
        assert_eq!(ring_self_intersection(&bowtie), Some((0, 2)));
    }

    /// A touch is as fatal as a crossing: the ring has no strict inside at the pinch point.
    #[test]
    fn a_pinch_touching_without_crossing_is_rejected() {
        let pinch = ring(&[
            [0.0, 0.0],
            [2.0, 2.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [2.0, 2.0],
            [0.0, 4.0],
        ]);
        assert!(ring_self_intersection(&pinch).is_some());
    }

    /// A spike shorter than the edge it retraces, and one longer. Neither *isolates* the two-sided
    /// adjacency test — with four or more points a long spike puts the ring's start point in the
    /// interior of the doubling-back edge, so the non-adjacent rule catches it first. The case that
    /// needs both directions is the collinear triple below.
    #[test]
    fn a_spike_doubling_back_is_rejected_either_length() {
        let short = ring(&[[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 4.0]]);
        assert!(ring_self_intersection(&short).is_some());
        let long = ring(&[[0.0, 0.0], [4.0, 0.0], [-2.0, 0.0], [0.0, 4.0]]);
        assert!(ring_self_intersection(&long).is_some());
    }

    #[test]
    fn a_repeated_consecutive_point_is_a_zero_length_edge() {
        let dup = ring(&[[0.0, 0.0], [4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]);
        assert_eq!(ring_self_intersection(&dup), Some((1, 1)));
        // Also across the wrap: the last point repeats the first.
        let closed = ring(&[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0], [0.0, 0.0]]);
        assert_eq!(ring_self_intersection(&closed), Some((3, 3)));
    }

    /// A ring of collinear points encloses nothing, and a signed area cannot orient it. It falls
    /// out of the same rule — the return leg always overlaps the outbound one.
    ///
    /// **This is the net for the two-sided adjacency test.** A brute force over every ring of 3–5
    /// points on a 4×4 grid found the one-sided and two-sided predicates disagreeing on 88 rings,
    /// *all* of them collinear triples like this one: at each corner the shared point sits at an
    /// end of the retraced span, so only the `u ∈ [s, v]` direction fires.
    #[test]
    fn a_zero_area_ring_is_rejected() {
        let flat = ring(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
        assert!(ring_self_intersection(&flat).is_some());
    }
}
