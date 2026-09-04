//! Planar polygon triangulation with holes (design §5).
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
/// projection preserves or reverses orientation — is not consulted at all. It used to
/// be (`if a[k] >= 0.0 { .. } else { swap }`), which put the whole frame on an `f64`
/// sign computed by summing `n` cross products. [`ring_orientation`] settles the same
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
/// triangulation road (design §5): monotone decomposition, then the Lawson pass.
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
pub(crate) fn triangulate_uv(
    uv: &mut Vec<P2>,
    rings: &[&[usize]],
    candidates: &[P2],
) -> Result<Vec<[usize; 3]>, TessError> {
    if ring_orientation(rings[0], uv) != 1 {
        return Err(TessError::DegenerateRing);
    }
    if rings[1..].iter().any(|h| ring_orientation(h, uv) != -1) {
        return Err(TessError::HoleWinding);
    }
    let mut out = Vec::new();
    for piece in monotone::decompose(uv, rings)? {
        monotone::triangulate_monotone(uv, &piece, &mut out)?;
    }
    // The decomposition answers *whether* the face meshes; this answers *how well*.
    // It moves diagonals only — never the rings — so the count, the area and the
    // boundary are the same on both sides of it.
    let constrained: HashSet<(usize, usize)> = rings
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
    Ok(out)
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
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// The goldens below are written **as a chart** — the numbers are `uv` already.
    ///
    /// ★ They used to be lifted to `z = 0` and handed to `triangulate_polygon`, the projection
    /// wrapper the bootstrap OBJ writer needed. That wrapper is gone (the projection rule lives
    /// once, in `planar_chart`), and for a `z = 0` ring its projection was the identity — so no
    /// coordinate here changes.
    fn pts(v: &[[f64; 2]]) -> Vec<P2> {
        v.to_vec()
    }

    /// The one road, given a golden's rings — a slice adapter, not a second rule: no projection,
    /// no handedness repair, no candidates.
    fn tri(uv: &[P2], outer: &[usize], holes: &[&[usize]]) -> Result<Vec<[usize; 3]>, TessError> {
        let mut v = uv.to_vec();
        let refs: Vec<&[usize]> = std::iter::once(outer)
            .chain(holes.iter().copied())
            .collect();
        triangulate_uv(&mut v, &refs, &[])
    }

    /// The checks every golden gets, each catching a different fault: the count pins
    /// the triangulation's identity, the **unsigned** area sum pins that no triangle
    /// escaped the polygon (a signed sum cancels and shoelace agrees), and all-CCW
    /// pins that none is folded.
    fn check(p: &[P2], tris: &[[usize; 3]], want_n: usize, want_area: f64) {
        assert_eq!(tris.len(), want_n, "triangle count");
        let mut sum = 0.0;
        for t in tris {
            let (a, b, c) = (p[t[0]], p[t[1]], p[t[2]]);
            let cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            assert!(cr > 0.0, "triangle {t:?} is not CCW");
            sum += 0.5 * cr.abs();
        }
        assert!((sum - want_area).abs() < 1e-12, "area {sum} vs {want_area}");
    }

    /// **Is this a partition, or merely the right amount of area?**
    ///
    /// Count + area + winding can all agree while triangles overlap and leave a
    /// matching hole elsewhere. What cannot survive that is the edge bookkeeping of a
    /// simplicial complex, so this checks it directly:
    ///
    /// - every **directed** edge is used at most once — two triangles covering the
    ///   same ground the same way would repeat one;
    /// - an undirected edge is used **twice** (interior, once each way) or **once**
    ///   (boundary) — never more;
    /// - and the once-used directed edges are **exactly the input rings**, which is
    ///   what says the triangles cover the polygon that was asked for rather than
    ///   some other region of the same area.
    ///
    /// Ear clipping over a bridged ring was measured clipping a triangle straight
    /// across another hole's slit; the ring self-intersected and two more triangles
    /// were emitted before it stalled. Had it not stalled, that mesh would have
    /// passed the three checks above. This is the net for it.
    fn check_partition(tris: &[[usize; 3]], outer: &[usize], holes: &[&[usize]]) {
        let mut once: HashSet<(usize, usize)> = HashSet::new();
        for t in tris {
            for k in 0..3 {
                let e = (t[k], t[(k + 1) % 3]);
                assert!(once.insert(e), "directed edge {e:?} used twice");
            }
        }
        let mut boundary: HashSet<(usize, usize)> = HashSet::new();
        for ring in std::iter::once(outer).chain(holes.iter().copied()) {
            for w in 0..ring.len() {
                boundary.insert((ring[w], ring[(w + 1) % ring.len()]));
            }
        }
        let unmatched: HashSet<(usize, usize)> = once
            .iter()
            .filter(|&&(a, b)| !once.contains(&(b, a)))
            .copied()
            .collect();
        assert_eq!(
            unmatched, boundary,
            "the once-used directed edges are not the input rings"
        );
    }

    #[test]
    fn a_square_is_two_triangles() {
        let p = pts(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let t = tri(&p, &[0, 1, 2, 3], &[]).unwrap();
        check(&p, &t, 2, 1.0);
        check_partition(&t, &[0, 1, 2, 3], &[]);
    }

    #[test]
    fn a_u_polygon_is_not_star_shaped_from_its_first_vertex() {
        // The `u_prism` cap. Fanning it from (0,0) drags a triangle across the notch:
        // that triangle comes out clockwise and the unsigned area overshoots. Ear
        // clipping owes nothing to the first vertex.
        let p = pts(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 2.3],
            [2.0, 2.3],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let ring: Vec<usize> = (0..8).collect();
        let t = tri(&p, &ring, &[]).unwrap();
        // Bar 3×1, left prong 1×1, right prong 1×1.3.
        check(&p, &t, 6, 3.0 + 1.0 + 1.3);
        check_partition(&t, &ring, &[]);
    }

    #[test]
    fn a_reflex_corner_polygon() {
        // The L-prism cap: star-shaped from (0,0), so even the old fan was right —
        // which is exactly why the fan bug hid for so long.
        let p = pts(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let ring: Vec<usize> = (0..6).collect();
        let t = tri(&p, &ring, &[]).unwrap();
        check(&p, &t, 4, 3.0);
        check_partition(&t, &ring, &[]);
    }

    #[test]
    fn a_square_with_a_square_hole() {
        // The pocketed cube's lid: outer 4 + hole 4 = 8 vertices, one hole, so
        // 8 + 2·1 − 2 = 8 triangles. The hole winds CW about the outer's normal.
        let p = pts(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.3, 0.3],
            [0.3, 0.7],
            [0.7, 0.7],
            [0.7, 0.3],
        ]);
        let hole = [4, 5, 6, 7];
        let t = tri(&p, &[0, 1, 2, 3], &[&hole]).unwrap();
        check(&p, &t, 8, 1.0 - 0.16);
        check_partition(&t, &[0, 1, 2, 3], &[&hole]);
    }

    #[test]
    fn two_holes() {
        let p = pts(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [0.0, 1.0],
            [0.2, 0.2],
            [0.2, 0.8],
            [0.8, 0.8],
            [0.8, 0.2],
            [2.2, 0.2],
            [2.2, 0.8],
            [2.8, 0.8],
            [2.8, 0.2],
        ]);
        let (h0, h1) = ([4, 5, 6, 7], [8, 9, 10, 11]);
        let t = tri(&p, &[0, 1, 2, 3], &[&h0, &h1]).unwrap();
        check(&p, &t, 12 + 4 - 2, 3.0 - 2.0 * 0.36);
        check_partition(&t, &[0, 1, 2, 3], &[&h0, &h1]);
    }

    /// **The wall with two windows — the case that used to have no ear.**
    ///
    /// Both holes' first mutually visible ring vertex was the *same* outer corner, so
    /// the bridged ring visited it three times and the polygon was pinched there. A
    /// pinched ring is not simple, Meisters' two-ears theorem does not apply to it, and
    /// measurement agreed: all six convex vertices had a diagonal that genuinely crossed
    /// the other hole. The sweep never merges the rings, so the pinch is not built.
    ///
    /// Everything here is transcribed from the failing model (a hub wall with two fins
    /// through it, projected to its own plane) — coordinates, index layout, and the
    /// order the two inner loops arrive in, which is what decided where each bridge
    /// landed. Handing the same two holes in the other order used to triangulate fine.
    #[test]
    fn a_wall_with_two_windows() {
        let p = pts(&[
            // Outer: the 3 × 2 wall, CCW.
            [0.0, 1.0],
            [0.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            // Upper hole, CW.
            [1.0, 0.25959136597258015],
            [1.0, 0.7035578716424774],
            [2.0, 0.7035578716424774],
            [2.0, 0.25959136597258015],
            // Lower hole, CW.
            [1.0, -0.703557871642477],
            [1.0, -0.25959136597258003],
            [2.0, -0.25959136597258003],
            [2.0, -0.703557871642477],
        ]);
        let (upper, lower) = ([4, 5, 6, 7], [8, 9, 10, 11]);
        let t = tri(&p, &[0, 1, 2, 3], &[&upper, &lower]).unwrap();
        let window = |a: usize, b: usize| p[b][1] - p[a][1];
        check(
            &p,
            &t,
            12 + 2 * 2 - 2,
            3.0 * 2.0 - window(4, 5) - window(8, 9),
        );
        check_partition(&t, &[0, 1, 2, 3], &[&upper, &lower]);
    }

    /// **The corpus that measured the old triangulator's 27.6% failure rate**, now a
    /// standing gate on the whole function: random axis-aligned rectangles with one to
    /// four non-overlapping rectangular holes on an integer grid.
    ///
    /// Hand-picked fixtures only cover the degeneracies someone thought of. This one
    /// produces shared sweep coordinates, horizontal edges and collinear vertices by
    /// the thousand — which is what an axis-aligned CAD face is — and every sample must
    /// come back a genuine triangulation: the topological count, all CCW, the right
    /// area, and the edge bookkeeping of a simplicial complex.
    #[test]
    fn random_rectilinear_faces_triangulate() {
        let mut st = 0xC0FF_EE00_1234_5678u64;
        let mut lcg = move || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            st >> 33
        };
        let mut rng = move |lo: i64, hi: i64| lo + (lcg() % ((hi - lo + 1) as u64)) as i64;
        let (w, h) = (24i64, 16i64);
        let mut samples = 0;
        for _ in 0..2000 {
            let k = rng(1, 4);
            let mut boxes: Vec<[i64; 4]> = Vec::new();
            for _ in 0..k {
                for _try in 0..20 {
                    let x0 = rng(1, w - 3);
                    let y0 = rng(1, h - 3);
                    let x1 = x0 + rng(1, (w - 1 - x0).min(5));
                    let y1 = y0 + rng(1, (h - 1 - y0).min(5));
                    if boxes
                        .iter()
                        .any(|b| x0 <= b[2] && b[0] <= x1 && y0 <= b[3] && b[1] <= y1)
                    {
                        continue;
                    }
                    boxes.push([x0, y0, x1, y1]);
                    break;
                }
            }
            if boxes.is_empty() {
                continue;
            }
            let mut flat = vec![
                [0.0, 0.0],
                [w as f64, 0.0],
                [w as f64, h as f64],
                [0.0, h as f64],
            ];
            let mut holes: Vec<Vec<usize>> = Vec::new();
            let mut hole_area = 0.0;
            for b in &boxes {
                let base = flat.len();
                for c in [[b[0], b[1]], [b[0], b[3]], [b[2], b[3]], [b[2], b[1]]] {
                    flat.push([c[0] as f64, c[1] as f64]);
                }
                holes.push((base..base + 4).collect());
                hole_area += ((b[2] - b[0]) * (b[3] - b[1])) as f64;
            }
            let p = pts(&flat);
            let hs: Vec<&[usize]> = holes.iter().map(|v| v.as_slice()).collect();
            let t = tri(&p, &[0, 1, 2, 3], &hs).expect("meshes");
            check(
                &p,
                &t,
                flat.len() + 2 * holes.len() - 2,
                (w * h) as f64 - hole_area,
            );
            check_partition(&t, &[0, 1, 2, 3], &hs);
            samples += 1;
        }
        assert!(samples > 1500, "only {samples} samples reached the gate");
    }

    /// ★★★★ **Offered points, and the two places they can land.**
    ///
    /// A lattice over a square: the ones on the diagonal the sweep already drew arrive by the
    /// **edge split** (this is the road the real face uses — every placed candidate on a merged
    /// lateral face came in this way), the rest by the three-way split. The triangle count is the
    /// arithmetic of both roads at once — `2I + V − 2` holds however each point got in — and
    /// `check_partition` is what says nothing cracked: the once-used directed edges must still be
    /// exactly the input ring, so no interior point leaked onto the boundary.
    #[test]
    fn offered_points_land_inside_and_on_interior_edges() {
        let ring = [0, 1, 2, 3];
        let mut uv: Vec<P2> = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
        let lattice: Vec<P2> = (1..4)
            .flat_map(|i| (1..4).map(move |j| [f64::from(i), f64::from(j)]))
            .collect();
        let t = triangulate_uv(&mut uv, &[&ring], &lattice).unwrap();
        assert_eq!(uv.len(), 4 + 9, "every offered point was placed");
        assert_eq!(t.len(), 2 * 9 + 4 - 2);
        check_partition(&t, &ring, &[]);
        let area: f64 = t
            .iter()
            .map(|t| {
                let (a, b, c) = (uv[t[0]], uv[t[1]], uv[t[2]]);
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            })
            .sum();
        assert!((area - 16.0).abs() < 1e-12, "area {area}");
        assert!(t.iter().all(|t| {
            let (a, b, c) = (uv[t[0]], uv[t[1]], uv[t[2]]);
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0.0
        }));
    }

    /// ★★★★★ **A candidate on the boundary, or outside it, is dropped — never forced.**
    ///
    /// The boundary polyline is *shared* with the neighbouring face, so a vertex added to it that
    /// the neighbour does not know about is a T-vertex, and a T-vertex is a crack. That is the one
    /// thing this layer may not produce, so the sweep declines the point instead: the caller offers
    /// candidates, it does not place them.
    #[test]
    fn a_candidate_on_the_boundary_or_outside_is_declined() {
        let ring = [0, 1, 2, 3];
        let base: Vec<P2> = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
        for p in [
            [2.0, 0.0],  // on a boundary edge
            [0.0, 0.0],  // on a boundary vertex
            [5.0, 5.0],  // outside the ring
            [-1.0, 2.0], // outside, and level with the interior
        ] {
            let mut uv = base.clone();
            let t = triangulate_uv(&mut uv, &[&ring], &[p]).unwrap();
            assert_eq!(uv.len(), 4, "{p:?} was placed");
            assert_eq!(t.len(), 2, "{p:?} changed the mesh");
            check_partition(&t, &ring, &[]);
        }
    }

    #[test]
    fn degenerate_rings_are_errors() {
        let two = pts(&[[0.0, 0.0], [1.0, 0.0]]);
        assert!(matches!(
            tri(&two, &[0, 1], &[]),
            Err(TessError::DegenerateRing)
        ));
        let line = pts(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
        assert!(matches!(
            tri(&line, &[0, 1, 2], &[]),
            Err(TessError::DegenerateRing)
        ));
    }

    /// ★★★★★ **A boundary that meets itself, and the two different sentences that says.**
    ///
    /// A **touch** — a vertex on the boundary elsewhere — is a boundary this decomposition has
    /// no answer for (an exact tangency, which the sampled boundary reproduces when a sample lands
    /// on it). A **crossing** is a ring set that has no triangulation at all: any one would be
    /// wrong. Different claims about *the mesh*, different names — neither is a verdict on the
    /// solid, which is `validate`'s to give.
    ///
    /// ☑ **The crossing branch exists because the gap was measured.** A bow-tie already came back
    /// `DegenerateRing` from the sweep, but the hole below — crossing the outer ring with no
    /// vertex at either crossing — came back **`Ok` with eight confident, wrong triangles**. That
    /// is the one thing this layer may not do.
    #[test]
    fn a_boundary_that_meets_itself_is_named_by_how() {
        // A hole's apex exactly on an outer edge: a touch.
        let touch = pts(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [2.0, 0.0],
            [1.0, 2.0],
            [3.0, 2.0],
        ]);
        assert!(matches!(
            tri(&touch, &[0, 1, 2, 3], &[&[4, 5, 6]]),
            Err(TessError::SelfTouchingBoundary)
        ));
        // A hole hanging out through the outer ring's right edge: two crossings, no vertex at
        // either. This used to mesh.
        let cross = pts(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [3.0, 1.0],
            [3.0, 3.0],
            [5.0, 3.0],
            [5.0, 1.0],
        ]);
        assert!(matches!(
            tri(&cross, &[0, 1, 2, 3], &[&[4, 5, 6, 7]]),
            Err(TessError::DegenerateRing)
        ));
        // And a single ring crossing itself — the sweep caught this one all along.
        let bowtie = pts(&[[0.0, 0.0], [2.0, 2.0], [2.0, 0.0], [0.0, 2.0]]);
        assert!(matches!(
            tri(&bowtie, &[0, 1, 2, 3], &[]),
            Err(TessError::DegenerateRing)
        ));
    }

    /// ★★★★★ **A touch vouches for nothing — two ways the old check let a crossing through.**
    ///
    /// `self_touch` used to return on the *first* touch it found and only then, if it found none,
    /// scan for crossings. Two defects follow; this fixture is the first, the next test the second.
    ///
    /// **A touch hides a crossing elsewhere.** A hole whose apex sits exactly on the outer ring's
    /// bottom edge (a genuine tangency — its two neighbours are both inside) *and* whose far side
    /// pokes out through the right edge with no vertex at the crossing. The touch was found first
    /// and the crossing never looked at, so the answer was `SelfTouchingBoundary` — "a boundary
    /// this decomposition cannot draw" — when the true answer is `DegenerateRing`: any
    /// triangulation of this ring set would be wrong.
    ///
    /// ☑ Measured red first (`SelfTouchingBoundary`), then fixed.
    #[test]
    fn a_touch_must_not_hide_a_crossing_elsewhere() {
        // A real tangency at (2,0) plus a real crossing of x = 4 at y = 2.25.
        let p = pts(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [2.0, 0.0],
            [1.0, 2.0],
            [3.0, 3.0],
            [5.0, 1.5],
        ]);
        let r = tri(&p, &[0, 1, 2, 3], &[&[4, 5, 6, 7]]);
        assert!(matches!(r, Err(TessError::DegenerateRing)), "{r:?}");
    }

    /// **A touch *is* a crossing.** A diamond hole poking out through the right edge whose two
    /// crossings land exactly on its own vertices `(4,3)` and `(4,1)`. Both were recorded as
    /// touches and the crossing scan never ran. The one `orient2d` pair that tells them apart:
    /// at a vertex touching a segment, its two neighbours lie on the **same** side of that
    /// segment's line for a tangency (the boundary turns back) and on **opposite** sides for a
    /// crossing (the boundary passes through). Here `(3,2)` and `(5,2)` straddle `x = 4`.
    ///
    /// ☑ Measured red first (`SelfTouchingBoundary`), then fixed.
    #[test]
    fn a_touch_that_is_really_a_crossing_is_named_so() {
        let p = pts(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [3.0, 2.0],
            [4.0, 3.0],
            [5.0, 2.0],
            [4.0, 1.0],
        ]);
        let r = tri(&p, &[0, 1, 2, 3], &[&[4, 5, 6, 7]]);
        assert!(matches!(r, Err(TessError::DegenerateRing)), "{r:?}");
    }

    #[test]
    fn a_hole_wound_the_wrong_way_is_an_error() {
        // A CCW hole means the b-rep does not keep material on the left. That is a
        // broken solid — `validate`'s business — not something to silently reverse.
        let p = pts(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.3, 0.3],
            [0.7, 0.3],
            [0.7, 0.7],
            [0.3, 0.7],
        ]);
        let ccw_hole = [4, 5, 6, 7];
        assert!(matches!(
            tri(&p, &[0, 1, 2, 3], &[&ccw_hole]),
            Err(TessError::HoleWinding)
        ));
    }
}
