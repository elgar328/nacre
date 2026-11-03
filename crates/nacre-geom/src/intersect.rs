//! Closed-form intersections among analytic surfaces.
//!
//! This module holds the closed-form intersections and the ring predicates over straight steps,
//! kept as `pub mod` so callers write `intersect::three_planes(..)`. It is not the one home of
//! robustness-sensitive code: rings with arcs are [`crate::mixed`]'s. The intersections are the
//! **construction (coordinate)** side of the precision split: the returned `Point3` is an f64
//! *cache* — the defining surfaces are the truth. So near-degenerate inputs are gated by a
//! conditioning threshold (there is no exact answer to cache). The ring predicates decide on the
//! rational truth, by `nacre_exact::orient2d_rat`'s exact signs.

use crate::{Line, Plane};
use nacre_math::Point3;

// The types and the sign the rational ring predicates below are stated in.
pub use nacre_exact::{Orient, Rat, orient2d_rat};

/// The point where a line crosses a plane — `p = l(t)` with
/// `t = n·(o_p − o_l) / (n·d)`. `None` when `n·d` is exactly zero (the line runs in or
/// parallel to the plane — no unique crossing).
///
/// The parallel test is **exact**, not toleranced: the production caller (the rim-circle
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
/// (The exact predicate is used where a *decision* must be exact — the
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

/// Where a point sits relative to a closed 2-D ring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RingSide {
    Inside,
    Outside,
    /// On the ring itself — a vertex or a point along an edge.
    OnBoundary,
}

/// **The rightward ray's crossing rule, spelled once** — for [`point_in_ring_2d_rat`] and the
/// kernel's mixed (arc-bearing) ring parity above this crate,
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

// ---------------------------------------------------------------------------
// The ring predicates, on the rational truth.
//
// A profile's truth is the rational the author's decimal spelled (`Rat::from_decimal`), and the
// f64 realization can carry a *different sign*: three points collinear in decimal land a hair off
// the line in binary. `Profile2d::check` and `sketch::from_rings` therefore judge the truth, through
// `mixed`'s arc-bearing ports of these. `point_in_ring_2d_rat` serves `combinatorics::in_faces`;
// all three are what `mixed` is held to (`tests/mixed.rs`). Every sign is the total
// `nacre_exact::orient2d_rat` and every comparison `Rat`'s exact `Ord`.
// ---------------------------------------------------------------------------

/// Whether `p` is inside `ring`, decided by **exact** signs only, never by a tolerance.
///
/// Crossing parity along a ray: an edge that straddles `p`'s horizontal line contributes a
/// crossing when `p` lies on the correct side of it, which `orient2d_rat` answers exactly. The
/// half-open straddle test ([`ray_straddle`]) is what keeps a vertex exactly on the ray from being
/// counted twice.
///
/// The ring is assumed **simple** (no self-intersections); for such a ring parity and winding
/// agree, so this needs no fill rule. A self-intersecting ring is out of contract and nothing here
/// checks it — run [`ring_self_intersection_rat`] first. A sketch or profile is judged through
/// [`crate::mixed::point_in_mixed_ring`], which carries arcs.
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

/// Whether two closed rings meet at all — a proper crossing or a mere touch. Exact.
///
/// A profile's rings must be disjoint (a hole lies strictly inside the outer ring and strictly
/// outside its siblings), so "touching" is a failure just as much as "crossing"; the two are not
/// worth distinguishing here. `O(n·m)`, which is the right complexity for sketch-sized rings.
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

/// The first pair of edge indices where a closed ring meets itself, or `None` if the ring is a
/// **simple polygon**. Exact — every decision is an `orient2d_rat` sign.
///
/// [`rings_cross_rat`] cannot answer this for a ring against itself: consecutive edges share an
/// endpoint, so it would report every ring. The cases are therefore split here, in this order:
///
/// - A **zero-length edge** (a repeated consecutive point) is a failure of its own, reported as
///   `(i, i)`. It leaves a degenerate edge in the topology, and it breaks the adjacency rule's
///   premise that neighbours share exactly one point.
/// - **Adjacent edges** share one endpoint legitimately. They fail only when *collinear and
///   overlapping* — a spike that doubles back. Both directions must be tested: the doubling-back
///   edge may be longer than the one it retraces, in which case it is the *first* edge's far
///   endpoint that lies inside the second.
/// - **Non-adjacent edges**: any contact at all is a self-intersection. A mere touch counts —
///   a ring that pinches through a point has no strict inside there, so even-odd is not the
///   author's meaning of it.
///
/// A merely *collinear* vertex (a flat corner in the middle of a straight run) is **not** an
/// error: the neighbours are collinear but do not overlap.
///
/// A ring that passes is a simple polygon, and a simple polygon has nonzero area — which is the
/// unstated precondition of every winding decision taken from a signed area.
///
/// `O(n²)`, the right complexity for sketch-sized rings.
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

/// **Is `p` a flat corner** between its ring neighbours `prev` and `next` — strictly mid-run on
/// the straight segment they span: collinear, and equal to neither. Such a vertex's two walls are
/// one plane, so it has no three-plane definition, and dissolving it changes no geometry.
///
/// The strictness is load-bearing — what fails it must survive to be reported:
/// - a **repeated point** (`p == prev`) is a zero-length edge, an author's mistake to name rather
///   than to erase (a bounding-box test alone would accept it: it includes the ends);
/// - a **spike** (collinear but past the far neighbour) is a self-intersection, which
///   [`ring_self_intersection_rat`] names.
///
/// Exact, on the rational truth: three points collinear in decimal are often a hair off the line
/// in binary, and the decimal is the author's meaning.
pub fn flat_corner(prev: [Rat; 2], p: [Rat; 2], next: [Rat; 2]) -> bool {
    p != prev && p != next && orient2d_rat(prev, p, next) == 0 && on_segment_2d_rat(prev, next, p)
}

/// Whether the edges `[u, s]` and `[s, v]`, which share `s`, double back over one another.
pub(crate) fn spike_rat(u: [Rat; 2], s: [Rat; 2], v: [Rat; 2]) -> bool {
    orient2d_rat(u, s, v) == 0 && (on_segment_2d_rat(u, s, v) || on_segment_2d_rat(s, v, u))
}

/// Exact segment/segment test: `true` for a proper crossing **or** any touching contact.
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

/// Whether the **collinear** point `p` lies within segment `ab`'s extent (callers check
/// collinearity with `orient2d_rat` first, so this is a bounding-box question only).
pub(crate) fn on_segment_2d_rat(a: [Rat; 2], b: [Rat; 2], p: [Rat; 2]) -> bool {
    p[0] >= a[0].min(b[0])
        && p[0] <= a[0].max(b[0])
        && p[1] >= a[1].min(b[1])
        && p[1] <= a[1].max(b[1])
}

#[cfg(test)]
#[path = "tests/intersect.rs"]
mod tests;
