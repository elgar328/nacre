//! Planar polygon triangulation with holes (design §5).
//!
//! The mesh is a **cache**, so this works in `f64` — exact geometry is the truth
//! and lives above. What it must never do is return a *wrong* mesh: a ring that
//! cannot be triangulated yields [`TessError`], never a silent fan.
//!
//! Nothing here knows `Model` or `Tessellation`. It takes points and index rings,
//! which is why it can be tested on hand-built polygons, and why both the
//! provenance tessellator and the bootstrap OBJ writer can share it.

mod monotone;

use crate::TessError;
use nacre_math::{Point3, Vector3};
use std::collections::HashMap;

/// A point in the plane's projected frame.
type P2 = [f64; 2];

/// The ring's own normal (Newell): robust for a non-planar-ish polygon and, more
/// to the point, **independent of the surface**. A `Reversed` face's loop still
/// winds CCW about its outward normal, so asking the ring rather than the surface
/// removes the "`Orientation::Forward` only" assumption entirely.
fn newell(pts: &[Point3], ring: &[usize]) -> Vector3 {
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
fn drop_axis(n: Vector3) -> (usize, usize) {
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
fn ring_orientation(ring: &[usize], uv: &[P2]) -> i8 {
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
    match nacre_predicates::orient2d(p, c, q) {
        d if d > 0.0 => 1,
        d if d < 0.0 => -1,
        _ => 0,
    }
}

/// Triangulate a planar polygon with holes, indices into `pts`.
///
/// Each hole must wind **opposite** to `outer` — the b-rep invariant that every loop
/// keeps material on its left. A ring that violates it is a broken solid, not
/// something to quietly repair.
///
/// Which side is "out" is the caller's to know, and it is never asked: the frame comes
/// from the ring's own Newell normal, which lets a `Reversed` face mesh without anyone
/// consulting its surface. Returns `V + 2H − 2` triangles for `V` ring vertices and `H`
/// holes — a topological count, so it does not depend on how they were found.
///
/// **This used to bridge each hole into the outer ring and clip ears.** The bridge
/// repeats a vertex, which makes the ring non-simple, which is the hypothesis Meisters'
/// two-ears theorem needs — so a bridged ring can have no ear at all, and 27.6% of
/// random rectilinear faces with two to four holes hit that. The sweep in
/// [`monotone`] never merges rings, so the question does not arise.
pub(crate) fn triangulate_polygon(
    pts: &[Point3],
    outer: &[usize],
    holes: &[&[usize]],
) -> Result<Vec<[usize; 3]>, TessError> {
    if outer.len() < 3 || holes.iter().any(|h| h.len() < 3) {
        return Err(TessError::DegenerateRing);
    }
    let n = newell(pts, outer);
    if n.norm() <= 0.0 {
        return Err(TessError::DegenerateRing);
    }
    let (iu, iv) = drop_axis(n);

    // **Only the ring vertices are projected.** `to_obj` hands over the whole vertex
    // store — indices there are handles, so it has to — and projecting all of it once
    // per face was already waste. It would have become worse here, where the sweep
    // allocates per-vertex state: `faces × store` instead of `faces × ring`.
    let mut back: Vec<usize> = Vec::new();
    let mut local: HashMap<usize, usize> = HashMap::new();
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for r in std::iter::once(outer).chain(holes.iter().copied()) {
        rings.push(
            r.iter()
                .map(|&i| {
                    *local.entry(i).or_insert_with(|| {
                        back.push(i);
                        back.len() - 1
                    })
                })
                .collect(),
        );
    }
    let mut uv: Vec<P2> = back
        .iter()
        .map(|&i| {
            let c = pts[i].as_array();
            [c[iu], c[iv]]
        })
        .collect();

    // **The frame's handedness is measured, not asserted.** The outer ring is CCW about
    // its own Newell normal by construction, so if it comes out CW here the projection
    // reversed orientation, and swapping `u` and `v` mirrors it back. This used to be a
    // `debug_assert` over a shoelace sum, which meant that in release nobody checked and
    // a near-edge-on face could run the whole triangulation in a flipped frame.
    match ring_orientation(&rings[0], &uv) {
        1 => {}
        -1 => {
            for p in &mut uv {
                p.swap(0, 1);
            }
        }
        _ => return Err(TessError::DegenerateRing),
    }
    if rings[1..].iter().any(|h| ring_orientation(h, &uv) != -1) {
        return Err(TessError::HoleWinding);
    }

    let refs: Vec<&[usize]> = rings.iter().map(|r| r.as_slice()).collect();
    let mut out = Vec::new();
    for piece in monotone::decompose(&uv, &refs)? {
        monotone::triangulate_monotone(&uv, &piece, &mut out)?;
    }
    Ok(out.into_iter().map(|t| t.map(|i| back[i])).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn pts(v: &[[f64; 2]]) -> Vec<Point3> {
        v.iter()
            .map(|p| Point3::from_array([p[0], p[1], 0.0]))
            .collect()
    }

    /// The checks every golden gets, each catching a different fault: the count pins
    /// the triangulation's identity, the **unsigned** area sum pins that no triangle
    /// escaped the polygon (a signed sum cancels and shoelace agrees), and all-CCW
    /// pins that none is folded.
    fn check(p: &[Point3], tris: &[[usize; 3]], want_n: usize, want_area: f64) {
        assert_eq!(tris.len(), want_n, "triangle count");
        let mut sum = 0.0;
        for t in tris {
            let (a, b, c) = (p[t[0]], p[t[1]], p[t[2]]);
            let cr = (b - a).cross(c - a);
            assert!(cr[2] > 0.0, "triangle {t:?} is not CCW");
            sum += 0.5 * cr.norm();
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
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[]).unwrap();
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
        let t = triangulate_polygon(&p, &ring, &[]).unwrap();
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
        let t = triangulate_polygon(&p, &ring, &[]).unwrap();
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
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[&hole]).unwrap();
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
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[&h0, &h1]).unwrap();
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
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[&upper, &lower]).unwrap();
        let window = |a: usize, b: usize| p[b].as_array()[1] - p[a].as_array()[1];
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
            let t = triangulate_polygon(&p, &[0, 1, 2, 3], &hs).expect("meshes");
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

    #[test]
    fn a_plane_other_than_xy() {
        // Same square, on x = 5, wound CCW about −x. The projector must follow the
        // ring's own normal, not a coordinate convention.
        let p: Vec<Point3> = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]
            .iter()
            .map(|c| Point3::from_array([5.0, c[0], c[1]]))
            .collect();
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[]).unwrap();
        assert_eq!(t.len(), 2);
        let area: f64 = t
            .iter()
            .map(|t| 0.5 * (p[t[1]] - p[t[0]]).cross(p[t[2]] - p[t[0]]).norm())
            .sum();
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_rings_are_errors() {
        let p = pts(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(matches!(
            triangulate_polygon(&p, &[0, 1], &[]),
            Err(TessError::DegenerateRing)
        ));
        // Collinear points have no normal.
        let line = pts(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]);
        assert!(matches!(
            triangulate_polygon(&line, &[0, 1, 2], &[]),
            Err(TessError::DegenerateRing)
        ));
    }

    #[test]
    fn a_ring_is_always_ccw_about_its_own_normal() {
        // Reversing the outer ring does not make it "clockwise" — it flips the Newell
        // normal, and the same triangles come out with the opposite winding. This
        // function cannot know which side is out, and does not pretend to.
        let p = pts(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let fwd = triangulate_polygon(&p, &[0, 1, 2, 3], &[]).unwrap();
        let rev = triangulate_polygon(&p, &[3, 2, 1, 0], &[]).unwrap();
        assert_eq!(fwd.len(), rev.len());
        let z = |t: &[usize; 3]| (p[t[1]] - p[t[0]]).cross(p[t[2]] - p[t[0]])[2];
        assert!(fwd.iter().all(|t| z(t) > 0.0));
        assert!(rev.iter().all(|t| z(t) < 0.0));
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
            triangulate_polygon(&p, &[0, 1, 2, 3], &[&ccw_hole]),
            Err(TessError::HoleWinding)
        ));
    }
}
