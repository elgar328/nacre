//! Closed-form intersections among analytic surfaces.
//!
//! Robustness-sensitive intersection code is isolated in this one module
//! Everything here is the **construction (coordinate)** side of the
//! precision split: the returned `Line`/`Point3` are f64 *caches* —
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
// Re-exported for the same reason: `nacre-ops` consumes the rational ring predicates below and
// states their coordinates in `Rat` without needing its own view of the sign primitive.
pub use nacre_exact::{Orient, Rat, orient2d_rat};

/// `sin²θ` below which two plane normals count as parallel. Unit normals make
/// `‖n1 × n2‖² = sin²θ ∈ [0, 1]`, so this absolute cutoff is scale-free.
const PARALLEL_EPS: f64 = 1e-16;

/// The line where two planes meet, or `None` if they are parallel (or
/// coincident).
///
/// Closed form — no SSI march or spline cache. Direction is
/// `n1 × n2`; the base point is the point of the line closest to the origin.
/// Returns a bare [`Line`] (the closed form), not a marched approximation — the M7
/// shape for those is an open question this crate cannot spell (see the crate doc).
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

/// The point where a line crosses a plane — `p = l(t)` with
/// `t = n·(o_p − o_l) / (n·d)`. `None` when `n·d` is exactly zero (the line runs in or
/// parallel to the plane — no unique crossing).
///
/// ★ The parallel test is **exact**, not toleranced: the production caller (the rim-circle
/// derivation, cylinder axis × cap plane) is perpendicular by construction, so `n·d ≈ ±|n|` is
/// far from zero; a near-parallel pair would give a far-away but well-defined crossing, which
/// is the honest answer to the question asked.
pub fn line_plane(l: &Line, p: &Plane) -> Option<Point3> {
    let n = p.normal();
    let denom = n.dot(l.direction());
    if denom == 0.0 {
        return None;
    }
    let t = n.dot(p.origin() - l.origin()) / denom;
    Some(l.point_at(t))
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
/// The f64 coordinate is a **cache** — the three planes are the truth.
/// Its residual to the planes is the vertex cache's measured tolerance,
/// which the caller measures when it forms the vertex (closed-form
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
/// geom→predicates handoff.
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
/// are generally irrational (Attene 2020).
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
/// coordinate *is* the truth; for a measured (boolean-made) vertex it is a rounded cache of
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
/// side of the precision split, exactly like [`plane_pair_dir_sign`].
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
/// split.
///
/// **Reads the rows from [`Plane::coefficients`] (the un-normalized `raw`), not [`Plane::normal`].**
/// `raw` is a positive multiple of `normal`, and a determinant's sign is invariant under scaling
/// each row by a positive factor, so this is the *same* function either way — except that `normal`
/// is `raw / ‖raw‖`, and that division rounds each component differently, which can destroy an exact
/// cancellation. Two exactly-parallel walls of a slanted prism have `raw` rows that are exact
/// negatives (`det = 0`), but their unit rows round to `det = ±1`; the guard then admitted a triple
/// its own consumer ([`three_plane_orient3d`], also on `coefficients`) rejects as `D = 0`, aborting.
/// Feeding the guard the same `raw` the consumer reads keeps the two in lockstep. (The
/// invariant is already stated on `Plane`: exact predicates take `coefficients`.)
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
/// split, like [`plane_side`] and [`plane_pair_dir_sign`] — not a
/// `Plane` method with a length tolerance. Scale-invariant and direction-agnostic:
/// opposite normals still name the same plane.
pub fn planes_coplanar(a: &Plane, b: &Plane) -> bool {
    nacre_predicates::planes_coplanar(a.coefficients(), b.coefficients())
}

/// Exact forward-ray/triangle crossing for kernel types — the geom→predicates
/// handoff for [`nacre_predicates::ray_triangle_cross`] (the
/// point-in-polyhedron test). `d` is the ray direction; `tri` is a single triangle
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
    let sgn = |v: f64, at: f64| {
        if v > at {
            Orient::Positive
        } else if v < at {
            Orient::Negative
        } else {
            Orient::Zero
        }
    };
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i].as_array(), ring[(i + 1) % n].as_array());
        let side = orient2d(a, b, pa);
        if side == 0.0 && on_segment_2d(a, b, pa) {
            return RingSide::OnBoundary;
        }
        match ray_step_crossing(sgn(a[1], pa[1]), sgn(b[1], pa[1]), sgn(side, 0.0)) {
            Some(true) => inside = !inside,
            Some(false) => {}
            // A straddling step through the probe is on the segment — already answered above;
            // the rule's own word for it is the same.
            None => return RingSide::OnBoundary,
        }
    }
    if inside {
        RingSide::Inside
    } else {
        RingSide::Outside
    }
}

/// **The rightward ray's crossing rule, spelled once** — for [`point_in_ring_2d`],
/// [`point_in_ring_2d_rat`], and the kernel's mixed (arc-bearing) ring parity above this crate,
/// which reads it lazily through [`ray_straddle`].
///
/// A ring step `a → b` against the probe's ray `{y = p.y, x > p.x}`: `ya` and `yb` are the signs
/// of `a.y − p.y` and `b.y − p.y`. The rule is **half-open in y**: an end *on* the ray is "not
/// above", so the step straddles the ray iff exactly one end is strictly above. A corner on the
/// ray is therefore counted by exactly one of its two steps — the one that leaves it upward — and
/// a step lying along the ray by neither: two steps leaving a corner to opposite sides count
/// once (a crossing), two leaving to the same side count twice or not at all (a touch). That is
/// what makes a corner on the ray a decision and not a tie; only the probe *on* the boundary is
/// left to the caller.
///
/// `None`: the step does not straddle. `Some(upward)`: it does, going up (`b` strictly above) or
/// down.
pub fn ray_straddle(ya: Orient, yb: Orient) -> Option<bool> {
    let above = |y: Orient| y == Orient::Positive;
    (above(ya) != above(yb)).then(|| above(yb))
}

/// Whether the step crosses the ray: a straddling step crosses iff the probe is to its left when
/// it goes up and to its right when it goes down — `side` is the sign of `orient2d(a, b, p)`,
/// which a caller may compute only once [`ray_straddle`] says the step straddles. `Some(false)`
/// for a step that does not straddle; `None` when a straddling step has `side == Zero`: the probe
/// is on the step, a boundary and not a crossing.
pub fn ray_step_crossing(ya: Orient, yb: Orient, side: Orient) -> Option<bool> {
    match ray_straddle(ya, yb) {
        None => Some(false),
        Some(upward) => match side {
            Orient::Zero => None,
            s => Some((s == Orient::Positive) == upward),
        },
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

// ---------------------------------------------------------------------------
// Rational twins of the ring predicates.
//
// A profile's truth is the rational the author's decimal spelled (`Rat::from_decimal`), and the
// f64 realization can carry a *different sign*: three points collinear in decimal land a hair off
// the line in binary. `Profile2d::check` and `sketch::from_rings` therefore judge the truth, and
// these are the predicates they do it with — verbatim ports of the f64 versions above (which stay:
// tessellation and the f64 fallback path still consume them), with `orient2d` swapped for the
// total `nacre_exact::orient2d_rat` and f64 comparisons for `Rat`'s exact `Ord`. Everything
// else in the originals is comparison, min/max, and array equality, so nothing changes meaning.
// ---------------------------------------------------------------------------

/// [`point_in_ring_2d`] over the ring's rational truth. Same crossing-parity contract, including
/// the assumption that the ring is simple.
pub fn point_in_ring_2d_rat(p: [Rat; 2], ring: &[[Rat; 2]]) -> RingSide {
    let n = ring.len();
    if n < 3 {
        return RingSide::Outside;
    }
    let sgn = |o: core::cmp::Ordering| match o {
        core::cmp::Ordering::Greater => Orient::Positive,
        core::cmp::Ordering::Less => Orient::Negative,
        core::cmp::Ordering::Equal => Orient::Zero,
    };
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        let side = orient2d_rat(a, b, p);
        if side == 0 && on_segment_2d_rat(a, b, p) {
            return RingSide::OnBoundary;
        }
        match ray_step_crossing(
            sgn(a[1].cmp(&p[1])),
            sgn(b[1].cmp(&p[1])),
            sgn(side.cmp(&0)),
        ) {
            Some(true) => inside = !inside,
            Some(false) => {}
            None => return RingSide::OnBoundary,
        }
    }
    if inside {
        RingSide::Inside
    } else {
        RingSide::Outside
    }
}

/// [`rings_cross`] over the rings' rational truth — any contact at all, proper or touching.
pub fn rings_cross_rat(a: &[[Rat; 2]], b: &[[Rat; 2]]) -> bool {
    let (n, m) = (a.len(), b.len());
    if n < 2 || m < 2 {
        return false;
    }
    for i in 0..n {
        let (p1, p2) = (a[i], a[(i + 1) % n]);
        for j in 0..m {
            if segments_meet_2d_rat(p1, p2, b[j], b[(j + 1) % m]) {
                return true;
            }
        }
    }
    false
}

/// [`ring_self_intersection`] over the ring's rational truth — same case split (zero-length
/// edge first, then cyclic-adjacent spikes, then any contact between non-adjacent edges), same
/// "a flat corner is not an error" stance.
pub fn ring_self_intersection_rat(ring: &[[Rat; 2]]) -> Option<(usize, usize)> {
    let n = ring.len();
    if n < 3 {
        return None;
    }
    let pt = |i: usize| ring[i % n];
    for i in 0..n {
        if pt(i) == pt(i + 1) {
            return Some((i, i));
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let hit = if j == i + 1 {
                spike_rat(pt(i), pt(i + 1), pt(i + 2))
            } else if i == 0 && j == n - 1 {
                spike_rat(pt(1), pt(0), pt(n - 1))
            } else {
                segments_meet_2d_rat(pt(i), pt(i + 1), pt(j), pt(j + 1))
            };
            if hit {
                return Some((i, j));
            }
        }
    }
    None
}

/// The ring with every flat corner dissolved — the profile constructor's lossless
/// normalization pass (a collinear vertex's two walls are one plane,
/// so the vertex has no three-plane definition; deleting it changes no geometry).
///
/// A vertex dissolves only when it is **strictly interior** to the segment its neighbours span:
/// collinear, and equal to neither neighbour. That strictness is load-bearing —
/// - a **repeated point** (`p == prev`) must survive, so `check` can still name it
///   `ZeroLengthProfileEdge` ("you typed the same point twice" is an author's mistake to report,
///   not to erase);
/// - a **spike** (collinear but past the far neighbour) must survive, so `check` still reports
///   the self-intersection.
///
/// Removal can make the two ex-neighbours' own corners newly flat (four points on one line), so
/// the scan repeats to a fixpoint. A ring that collapses below three points is returned as-is
/// for `check` to reject as degenerate — that a fully-collinear "ring" encloses nothing is the
/// honest report.
pub fn drop_collinear_midpoints(mut ring: Vec<[Rat; 2]>) -> Vec<[Rat; 2]> {
    loop {
        let n = ring.len();
        if n < 3 {
            return ring;
        }
        let flat = (0..n).find(|&i| {
            let (prev, p, next) = (ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]);
            // Collinear + inside the neighbours' box + distinct from both = strictly between.
            // Reusing `on_segment_2d_rat` alone would be wrong: it includes the endpoints, and
            // an endpoint hit here is a zero-length edge that must survive to be reported.
            p != prev
                && p != next
                && orient2d_rat(prev, p, next) == 0
                && on_segment_2d_rat(prev, next, p)
        });
        match flat {
            Some(i) => {
                ring.remove(i);
            }
            None => return ring,
        }
    }
}

/// [`spike`]'s rational twin.
pub(crate) fn spike_rat(u: [Rat; 2], s: [Rat; 2], v: [Rat; 2]) -> bool {
    orient2d_rat(u, s, v) == 0 && (on_segment_2d_rat(u, s, v) || on_segment_2d_rat(s, v, u))
}

/// [`segments_meet_2d`]'s rational twin: proper crossing or any touching contact.
pub(crate) fn segments_meet_2d_rat(p1: [Rat; 2], p2: [Rat; 2], q1: [Rat; 2], q2: [Rat; 2]) -> bool {
    let d1 = orient2d_rat(q1, q2, p1);
    let d2 = orient2d_rat(q1, q2, p2);
    let d3 = orient2d_rat(p1, p2, q1);
    let d4 = orient2d_rat(p1, p2, q2);
    if ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0
    {
        return true;
    }
    (d1 == 0 && on_segment_2d_rat(q1, q2, p1))
        || (d2 == 0 && on_segment_2d_rat(q1, q2, p2))
        || (d3 == 0 && on_segment_2d_rat(p1, p2, q1))
        || (d4 == 0 && on_segment_2d_rat(p1, p2, q2))
}

/// [`on_segment_2d`]'s rational twin — the **collinear** point `p` within `ab`'s box.
pub(crate) fn on_segment_2d_rat(a: [Rat; 2], b: [Rat; 2], p: [Rat; 2]) -> bool {
    p[0] >= a[0].min(b[0])
        && p[0] <= a[0].max(b[0])
        && p[1] >= a[1].min(b[1])
        && p[1] <= a[1].max(b[1])
}

#[cfg(test)]
#[path = "tests/intersect.rs"]
mod tests;
