//! Planar polygon triangulation with holes.
//!
//! The mesh is a **cache**, so this works in `f64` — exact geometry is the truth
//! and lives above. What it must never do is return a *wrong* mesh: a ring that
//! cannot be triangulated yields [`TessError`], never a silent fan.
//!
//! Nothing here knows `Model` or `Tessellation`. It takes points and index rings,
//! which is why it can be tested on hand-built polygons, and why both the
//! provenance tessellator and the bootstrap OBJ writer can share it.

mod delaunay;
mod monotone;
mod sos;

pub(crate) use monotone::{Meets, Touch, TouchKind, Witness};
use sos::Twins;

/// Two rings that share a vertex — the same point at two indices, one in each — that the caller
/// wants bridged into one ring there. Which of the two is the ring whose straight segment was
/// split is decided here, by geometry, not by the caller: it is the one whose neighbours at the
/// shared point are collinear with it.
/// What [`triangulate_uv`] hands back: the triangles, and the rings it actually triangulated —
/// the input rings, or with a bridged pair spliced into one.
pub(crate) type Triangulated = (Vec<[usize; 3]>, Vec<Vec<usize>>);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Bridge {
    pub(crate) ring_x: usize,
    pub(crate) x: usize,
    pub(crate) ring_y: usize,
    pub(crate) y: usize,
}

use crate::TessError;
use nacre_math::{Point3, Vector3};
use nacre_predicates::orient2d;
use std::collections::HashSet;

/// A point in the plane's projected frame.
pub(crate) type P2 = [f64; 2];

/// The ring's own normal (Newell): robust for a non-planar-ish polygon and, more
/// to the point, **independent of the surface**. A `Reversed` face's loop still
/// winds CCW about its outward normal, so asking the ring rather than the surface
/// removes the "`Orientation::Forward` only" assumption entirely.
pub(crate) fn newell(pts: &[Point3], ring: &[usize]) -> Vector3 {
    let mut n = [0.0f64; 3];
    for w in 0..ring.len() {
        let a = pts[ring[w]].as_array();
        let b = pts[ring[(w + 1) % ring.len()]].as_array();
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    Vector3::from_array(n)
}

/// Which two coordinates to keep: drop the axis the normal leans on hardest.
///
/// **Only the magnitude is read.** Any axis the polygon does not lie edge-on to gives
/// a non-degenerate projection, and the normal's *sign* — which decides whether that
/// projection preserves or reverses orientation — is not consulted at all. Consulting it
/// (`if a[k] >= 0.0 { .. } else { swap }`) would put the whole frame on an `f64` sign
/// computed by summing `n` cross products. [`ring_orientation`] settles the same
/// question afterwards with one exact predicate, on the ring itself.
pub(crate) fn drop_axis(n: Vector3) -> (usize, usize) {
    let a = n.as_array();
    let k = (0..3)
        .max_by(|&i, &j| a[i].abs().total_cmp(&a[j].abs()))
        .unwrap();
    // (y,z) / (z,x) / (x,y) — the cyclic choices, up to the handedness settled later.
    [(1, 2), (2, 0), (0, 1)][k]
}

/// `a` before `b` in the sweep's total order: **`v` descending, then `u` ascending**.
///
/// Every vertex comparison in this module goes through here, and it is exact — the
/// projected coordinates are copied `f64`s, not computed ones, so `<` and `==` on them
/// answer about the real numbers. Making the order *total* is what lets a face full of
/// shared `v` values (which is what an axis-aligned CAD face is) be handled without a
/// single tie.
fn lex_less(a: P2, b: P2) -> bool {
    a[1] > b[1] || (a[1] == b[1] && a[0] < b[0])
}

/// `+1` CCW, `-1` CW, `0` the ring is degenerate — **exactly**.
///
/// At the ring's *last* vertex in [`lex_less`] order (lowest `v`, then greatest `u`)
/// a simple polygon's interior angle is strictly convex: both neighbours lie above, so
/// they cannot be collinear with it unless the ring doubles back on itself. One
/// `orient2d` there therefore decides the whole ring — where a shoelace sum would need
/// the exact addition of `n` products to say the same thing, and `f64` addition does
/// not give it.
pub(crate) fn ring_orientation(ring: &[usize], uv: &[P2]) -> i8 {
    let n = ring.len();
    let mut m = 0;
    for i in 1..n {
        if lex_less(uv[ring[m]], uv[ring[i]]) {
            m = i;
        }
    }
    let (p, c, q) = (
        uv[ring[(m + n - 1) % n]],
        uv[ring[m]],
        uv[ring[(m + 1) % n]],
    );
    match orient2d(p, c, q) {
        d if d > 0.0 => 1,
        d if d < 0.0 => -1,
        _ => 0,
    }
}

/// Triangulate rings that are **already in a chart** — the chart-free core of this layer's one
/// triangulation road: monotone decomposition, then the Lawson pass.
///
/// ★★ **Precondition: the outer ring is CCW and every hole is CW, in `uv`.** The repair is the
/// *caller's* because charts differ in what a repair costs. A planar face's two axes are
/// interchangeable, so swapping them mirrors the frame for free — that is what
/// `planar_chart` does. A curved face's are not: `v` is the sweep axis, and swapping it
/// would send the sweep the wrong way across the surface, which is how a triangulation grows
/// diagonals that leave the surface. Such a caller negates one axis instead. Both repairs are
/// checked here rather than trusted.
///
/// ★★ **`candidates` are points the caller offers the interior**, in its own order; `uv` grows by
/// the ones that were taken, in that same order, so a caller can map the tail back through its
/// chart. A candidate with nowhere to go is **dropped, not forced** — see [`insert_interior`].
/// Every way a ring set's boundary meets itself, without triangulating it — the chart layer's
/// question before it decides to bridge a touch. Rings are linked exactly as [`triangulate_uv`]
/// would link them, so the answer is the one the sweep would refuse on.
pub(crate) fn meets(uv: &[P2], rings: &[&[usize]]) -> Result<Meets, TessError> {
    let (prev, next) = monotone::link(uv.len(), rings)?;
    Ok(monotone::self_touch(uv, &prev, &next))
}

///
/// ★★★★★ **`bridges` — a hole touching another ring at one point is spliced into it.** The
/// touching point sits at two indices (the chart layer put the curved ring's sample into the
/// straight edge it touches, so both rings carry it); the two rings are joined there into one
/// ring, `a[..=i] ++ b[j+1..] ++ b[..=j] ++ a[i+1..]`, the textbook bridge with a cut of length
/// zero. Both rings are traversed in their stored direction, so same-winding rings merge into
/// a ring of that winding (two holes stay a hole). The winding gate runs on the rings as
/// given, *before* the splice, so a mis-wound hole is still named. The returned rings are the
/// ones actually triangulated, and a caller that reads boundary edges off rings must use them.
///
/// The one pair of coincident indices this leaves in the ring is handed to the sweep as
/// [`Twins`], with the order [`sos`] defines. At most one bridge per face is bridged; a second
/// would need a second symbolic pair, and that population has not been seen — it comes back
/// [`TessError::SelfTouchingBoundary`], as it always did.
pub(crate) fn triangulate_uv(
    uv: &mut Vec<P2>,
    rings: &[&[usize]],
    candidates: &[P2],
    bridges: &[Bridge],
) -> Result<Triangulated, TessError> {
    if ring_orientation(rings[0], uv) != 1 {
        return Err(TessError::DegenerateRing);
    }
    if rings[1..].iter().any(|h| ring_orientation(h, uv) != -1) {
        return Err(TessError::HoleWinding);
    }
    let mut rings_used: Vec<Vec<usize>> = rings.iter().map(|r| r.to_vec()).collect();
    let mut twins: Option<Twins> = None;
    if let Some(&bridge) = bridges.first() {
        if bridges.len() > 1 {
            return Err(TessError::SelfTouchingBoundary);
        }
        twins = Some(splice(uv, &mut rings_used, bridge)?);
    }
    let refs: Vec<&[usize]> = rings_used.iter().map(|r| r.as_slice()).collect();
    let mut out = Vec::new();
    for piece in monotone::decompose(uv, &refs, twins.as_ref())? {
        monotone::triangulate_monotone(uv, &piece, &mut out, twins.as_ref())?;
    }
    // The decomposition answers *whether* the face meshes; this answers *how well*.
    // It moves diagonals only — never the rings — so the count, the area and the
    // boundary are the same on both sides of it.
    let constrained: HashSet<(usize, usize)> = refs
        .iter()
        .flat_map(|r| {
            (0..r.len()).map(move |k| {
                let (a, b) = (r[k], r[(k + 1) % r.len()]);
                (a.min(b), a.max(b))
            })
        })
        .collect();
    delaunay::refine(uv, &mut out, &constrained);
    if !candidates.is_empty() {
        insert_interior(uv, &mut out, &constrained, candidates);
        // The insertions are what make the flip pass able to do anything here: a long edge stays
        // Delaunay while nothing sits inside its circumcircle, and these points are what sit there.
        delaunay::refine(uv, &mut out, &constrained);
    }
    Ok((out, rings_used))
}

/// Join two rings at their shared point into one, and name the twins that leaves.
///
/// The split ring is the one whose two neighbours of the shared point are collinear with it —
/// that is what a straight edge split at a sample looks like — and the other ring's neighbours
/// are strictly off that line (the tangency witness certified as much). Both collinear or
/// neither is not the shape this bridges and is refused by the sweep's own name for it.
fn splice(uv: &[P2], rings: &mut Vec<Vec<usize>>, b: Bridge) -> Result<Twins, TessError> {
    let straight = |ring: &[usize], at: usize| -> bool {
        let n = ring.len();
        let (p, c, q) = (ring[(at + n - 1) % n], ring[at], ring[(at + 1) % n]);
        orient2d(uv[p], uv[c], uv[q]) == 0.0
    };
    let (sx, sy) = (
        straight(&rings[b.ring_x], b.x),
        straight(&rings[b.ring_y], b.y),
    );
    let (ra, i, rb, j) = match (sx, sy) {
        (true, false) => (b.ring_x, b.x, b.ring_y, b.y),
        (false, true) => (b.ring_y, b.y, b.ring_x, b.x),
        _ => return Err(TessError::SelfTouchingBoundary),
    };
    let (a, hb) = (rings[ra].clone(), rings[rb].clone());
    let (o, h) = (a[i], hb[j]);
    let mut merged: Vec<usize> = Vec::with_capacity(a.len() + hb.len());
    merged.extend_from_slice(&a[..=i]);
    merged.extend_from_slice(&hb[j + 1..]);
    merged.extend_from_slice(&hb[..=j]);
    merged.extend_from_slice(&a[i + 1..]);
    // `o` keeps its on-line predecessor, `h` its on-line successor: the two ends of the split
    // segment, and the directions the twins slide toward.
    let o_along = a[(i + a.len() - 1) % a.len()];
    let h_along = a[(i + 1) % a.len()];
    rings[ra] = merged;
    rings.remove(rb);
    Ok(Twins {
        o,
        h,
        o_along,
        h_along,
    })
}

/// Whether `p` is **strictly** inside the counter-clockwise triangle `t` — exactly.
///
/// Strict on purpose: a point on an edge or at a corner would make a zero-area triangle if fanned
/// from there, so it is handled as an edge split or dropped instead.
fn strictly_inside(uv: &[P2], t: [usize; 3], p: P2) -> bool {
    (0..3).all(|k| orient2d(uv[t[k]], uv[t[(k + 1) % 3]], p) > 0.0)
}

/// Whether `p`, already known collinear with `a`–`b`, lies strictly between them — exactly.
///
/// Compares one coordinate, the one the segment actually varies along, so the answer is a pair of
/// `<` on stored `f64`s rather than a rounded dot product that a near-endpoint could flip.
fn strictly_between(a: P2, b: P2, p: P2) -> bool {
    let k = usize::from(a[0] == b[0]);
    (a[k] < p[k] && p[k] < b[k]) || (b[k] < p[k] && p[k] < a[k])
}

/// Place each candidate the sweep will accept, and drop the rest.
///
/// ★★★★★ **Two places, and which one is used depends on how the face is spelled.** A candidate may
/// fall **strictly inside** a triangle, which splits it in three; or **on an interior edge**, which
/// splits that edge — and with it both triangles sharing it — in two. The second is not the exotic
/// one: the caller's points sit on the lines where the face's shape changes, and a triangulation
/// naturally already has edges running along those lines.
///
/// **Both roads are live in production**, and the split is not gradual — measured over the five
/// merged lateral faces (a boss on each of four congruent walls, and on the corner):
///
/// | face | inside | on an edge | dropped |
/// |---|---|---|---|
/// | `−x`, `+x`, `+y` walls | **0** | 180 | 180 |
/// | corner | **0** | 270 | 90 |
/// | `−y` wall | **102** | 76 | 182 |
///
/// The odd one out is the wall where `ref_dir` puts θ = 0 *inside* the notch, so the notch stays an
/// honest **hole** instead of being bridged into the outer walk — and a hole's neighbourhood has no
/// long diagonal for the points to land on. So neither branch may be deleted as unreachable: which
/// one runs is decided by the seam's accident, not by the shape the user drew.
///
/// ★★ **A constrained edge is never split.** A point the neighbouring face does not know about is
/// a T-vertex, and a T-vertex is a crack — the one thing this layer may not produce. A candidate
/// that lands on the boundary, or outside the rings entirely, is simply dropped: the caller offers
/// points, it does not get to place them.
fn insert_interior(
    uv: &mut Vec<P2>,
    tris: &mut Vec<[usize; 3]>,
    constrained: &HashSet<(usize, usize)>,
    candidates: &[P2],
) {
    for &p in candidates {
        // ★ **One pass answers both questions.** The two places are mutually exclusive across the
        // whole triangulation — a point strictly inside one triangle is on no triangle's edge, and
        // vice versa — so there is no need to exhaust the "inside" search before starting the
        // "on an edge" one. Measured: scanning twice cost **8× per triangle** what an ordinary
        // band costs (9.48 µs vs 1.17 µs), almost all of it in the scan that was going to fail.
        let mut inside = None;
        let mut edge = None;
        'search: for (ti, t) in tris.iter().enumerate() {
            // A point outside the triangle's own bounding box is neither inside it nor on any of
            // its edges — four comparisons on stored coordinates that skip six exact predicates.
            let (mut lo, mut hi) = (uv[t[0]], uv[t[0]]);
            for &i in &t[1..] {
                lo = [lo[0].min(uv[i][0]), lo[1].min(uv[i][1])];
                hi = [hi[0].max(uv[i][0]), hi[1].max(uv[i][1])];
            }
            if p[0] < lo[0] || p[0] > hi[0] || p[1] < lo[1] || p[1] > hi[1] {
                continue;
            }
            if strictly_inside(uv, *t, p) {
                inside = Some(ti);
                break;
            }
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                if constrained.contains(&(a.min(b), a.max(b))) {
                    continue;
                }
                if orient2d(uv[a], uv[b], p) == 0.0 && strictly_between(uv[a], uv[b], p) {
                    edge = Some((a, b));
                    break 'search;
                }
            }
        }
        if let Some(k) = inside {
            let t = tris.remove(k);
            let m = uv.len();
            uv.push(p);
            for e in 0..3 {
                tris.push([t[e], t[(e + 1) % 3], m]);
            }
            continue;
        }
        // On an interior edge: split it, and with it every triangle that owns it.
        let Some((a, b)) = edge else { continue };
        let m = uv.len();
        uv.push(p);
        let mut fresh: Vec<[usize; 3]> = Vec::new();
        tris.retain(|t| {
            for e in 0..3 {
                // ★ Each owner meets the edge in its **own** direction — the two sharing it run
                // opposite ways — so the split is written from that rotation, never from `(a, b)`.
                let (x, y, c) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
                if (x == a && y == b) || (x == b && y == a) {
                    fresh.push([x, m, c]);
                    fresh.push([m, y, c]);
                    return false;
                }
            }
            true
        });
        tris.extend(fresh);
    }
}

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;
