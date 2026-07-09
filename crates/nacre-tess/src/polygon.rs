//! Planar polygon triangulation with holes (design §5).
//!
//! The mesh is a **cache**, so this works in `f64` — exact geometry is the truth
//! and lives above. What it must never do is return a *wrong* mesh: a ring that
//! cannot be triangulated yields [`TessError`], never a silent fan.
//!
//! Nothing here knows `Model` or `Tessellation`. It takes points and index rings,
//! which is why it can be tested on hand-built polygons, and why both the
//! provenance tessellator and the bootstrap OBJ writer can share it.

use crate::TessError;
use nacre_math::{Point3, Vector3};

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

/// Drop the axis the normal leans on hardest, keeping a right-handed `(u, v)` so
/// that a ring wound CCW about `n` stays CCW in 2D. Choosing an axis is exact —
/// no arithmetic, no tolerance.
fn projector(n: Vector3) -> impl Fn(Point3) -> P2 {
    let a = n.as_array();
    let k = (0..3)
        .max_by(|&i, &j| a[i].abs().total_cmp(&a[j].abs()))
        .unwrap();
    // (y,z) / (z,x) / (x,y) are the cyclic choices that keep `(u, v, e_k)` right
    // handed; a negative component flips handedness, so swap `u` and `v` back.
    let (iu, iv) = [(1, 2), (2, 0), (0, 1)][k];
    let (iu, iv) = if a[k] >= 0.0 { (iu, iv) } else { (iv, iu) };
    move |p: Point3| {
        let c = p.as_array();
        [c[iu], c[iv]]
    }
}

/// Twice the signed area of a ring. Positive ⇔ CCW.
fn area2(ring: &[usize], uv: &[P2]) -> f64 {
    let mut s = 0.0;
    for w in 0..ring.len() {
        let a = uv[ring[w]];
        let b = uv[ring[(w + 1) % ring.len()]];
        s += a[0] * b[1] - b[0] * a[1];
    }
    s
}

/// `> 0` when `c` is to the left of `a → b`.
fn cross2(a: P2, b: P2, c: P2) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Inside **or on the boundary**. A vertex merely touching the candidate diagonal
/// must veto the ear: the L-prism's cap has `(0,2)–(2,0)` running exactly through
/// its reflex corner `(1,1)`, and a strict test would clip an ear that folds the
/// next triangle backwards.
///
/// Corners are excluded by index at the call site, which also covers the two
/// vertices a bridge repeats.
fn inside_or_on(t: [P2; 3], p: P2) -> bool {
    cross2(t[0], t[1], p) >= 0.0 && cross2(t[1], t[2], p) >= 0.0 && cross2(t[2], t[0], p) >= 0.0
}

/// Do the open segments `a→b` and `c→d` properly cross, or touch other than at a
/// shared endpoint? Endpoint sharing is fine — a bridge starts on a ring vertex.
fn segments_conflict(a: P2, b: P2, c: P2, d: P2) -> bool {
    let shares =
        |p: P2, q: P2| (p[0] - q[0]).abs() < f64::EPSILON && (p[1] - q[1]).abs() < f64::EPSILON;
    if shares(a, c) || shares(a, d) || shares(b, c) || shares(b, d) {
        return false;
    }
    let (d1, d2) = (cross2(c, d, a), cross2(c, d, b));
    let (d3, d4) = (cross2(a, b, c), cross2(a, b, d));
    if (d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0) {
        return true;
    }
    // Collinear touching counts as a conflict: a bridge lying along an edge would
    // be traversed by three triangles, and the mesh would stop being watertight.
    let on = |p: P2, q: P2, r: P2| {
        cross2(p, q, r) == 0.0
            && r[0] >= p[0].min(q[0])
            && r[0] <= p[0].max(q[0])
            && r[1] >= p[1].min(q[1])
            && r[1] <= p[1].max(q[1])
    };
    on(c, d, a) || on(c, d, b) || on(a, b, c) || on(a, b, d)
}

/// Merge each hole into the outer ring with a bridge to a mutually visible vertex.
///
/// Brute force over vertex pairs: rings here have tens of vertices, not thousands,
/// and a visible pair always exists for a hole strictly inside a simple polygon.
/// Eberly's ray cast would be `O(n)` instead of `O(n·m·E)` and much easier to get
/// subtly wrong.
fn bridge_holes(outer: &[usize], holes: &[&[usize]], uv: &[P2]) -> Result<Vec<usize>, TessError> {
    let mut ring: Vec<usize> = outer.to_vec();
    // Merge the rightmost hole first: its bridge cannot be blocked by a hole that
    // is still further left.
    let mut order: Vec<usize> = (0..holes.len()).collect();
    let max_u = |h: &[usize]| {
        h.iter()
            .map(|&i| uv[i][0])
            .fold(f64::NEG_INFINITY, f64::max)
    };
    order.sort_by(|&i, &j| max_u(holes[j]).total_cmp(&max_u(holes[i])));

    for &hi in &order {
        let hole = holes[hi];
        let rest: Vec<&[usize]> = order
            .iter()
            .filter(|&&j| j != hi)
            .map(|&j| holes[j])
            .collect();
        let mut bridge = None;
        'search: for (bi, &b) in ring.iter().enumerate() {
            for (mi, &m) in hole.iter().enumerate() {
                let (pb, pm) = (uv[b], uv[m]);
                let blocked = |r: &[usize]| {
                    (0..r.len()).any(|w| {
                        let (c, d) = (uv[r[w]], uv[r[(w + 1) % r.len()]]);
                        segments_conflict(pb, pm, c, d)
                    })
                };
                if blocked(&ring) || blocked(hole) || rest.iter().any(|r| blocked(r)) {
                    continue;
                }
                bridge = Some((bi, mi));
                break 'search;
            }
        }
        let (bi, mi) = bridge.ok_or(TessError::NoBridge)?;

        // `outer[..=bi]` + the hole walked from `mi` all the way round + `outer[bi]`
        // again. Both bridge endpoints appear twice; no new point is created, so the
        // shared edge polylines stay shared and the mesh stays crack-free.
        let mut merged: Vec<usize> = ring[..=bi].to_vec();
        merged.extend(hole[mi..].iter().chain(hole[..=mi].iter()));
        merged.push(ring[bi]);
        merged.extend_from_slice(&ring[bi + 1..]);
        ring = merged;
    }
    Ok(ring)
}

/// Triangulate a planar polygon, indices into `pts`.
///
/// Each hole must wind **opposite** to `outer` — the b-rep invariant that every loop
/// keeps material on its left. A ring that violates it is a broken solid, not
/// something to quietly repair.
///
/// `outer`'s own winding cannot be checked here, and is not: the Newell normal comes
/// *from* the ring, so a ring is always CCW about it. Which side is "out" is the
/// caller's to know. That is the point — it lets a `Reversed` face mesh without
/// anyone consulting its surface normal.
///
/// Returns `outer.len() + Σ hole.len() + 2·holes.len() − 2` triangles: ear clipping
/// on the bridged ring, which gains two vertices per hole.
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
    let project = projector(n);
    let uv: Vec<P2> = pts.iter().map(|&p| project(p)).collect();

    debug_assert!(
        area2(outer, &uv) > 0.0,
        "a ring is CCW about its own normal"
    );
    if holes.iter().any(|h| area2(h, &uv) >= 0.0) {
        return Err(TessError::HoleWinding);
    }

    let mut ring = bridge_holes(outer, holes, &uv)?;
    let mut out = Vec::with_capacity(ring.len().saturating_sub(2));

    // Ear clipping. An ear is a convex corner whose triangle touches no *reflex*
    // vertex. Only reflex ones can veto: a convex vertex inside would mean the ring
    // self-intersects, and testing them too over-blocks valid ears. `O(n²)`; rings
    // here are tens of vertices.
    while ring.len() > 3 {
        let k = ring.len();
        let reflex = |j: usize| {
            let (a, b, c) = (ring[(j + k - 1) % k], ring[j], ring[(j + 1) % k]);
            cross2(uv[a], uv[b], uv[c]) <= 0.0
        };
        let ear = (0..k).find(|&i| {
            let (a, b, c) = (ring[(i + k - 1) % k], ring[i], ring[(i + 1) % k]);
            let tri = [uv[a], uv[b], uv[c]];
            if cross2(tri[0], tri[1], tri[2]) <= 0.0 {
                return false; // reflex or collinear
            }
            !(0..k).any(|j| {
                let v = ring[j];
                v != a && v != b && v != c && reflex(j) && inside_or_on(tri, uv[v])
            })
        });
        let Some(i) = ear else {
            // Self-intersecting or otherwise not a simple polygon. No correct mesh
            // exists; do not invent one.
            return Err(TessError::NoEar);
        };
        out.push([ring[(i + k - 1) % k], ring[i], ring[(i + 1) % k]]);
        ring.remove(i);
    }
    out.push([ring[0], ring[1], ring[2]]);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pts(v: &[[f64; 2]]) -> Vec<Point3> {
        v.iter()
            .map(|p| Point3::from_array([p[0], p[1], 0.0]))
            .collect()
    }

    /// The three checks every golden gets, and each catches a different fault:
    /// the count pins ear clipping's identity, the **unsigned** area sum pins that
    /// no triangle escaped the polygon (a signed sum cancels and shoelace agrees),
    /// and all-CCW pins that none is folded.
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

    #[test]
    fn a_square_is_two_triangles() {
        let p = pts(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let t = triangulate_polygon(&p, &[0, 1, 2, 3], &[]).unwrap();
        check(&p, &t, 2, 1.0);
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
        let t = triangulate_polygon(&p, &(0..8).collect::<Vec<_>>(), &[]).unwrap();
        // Bar 3×1, left prong 1×1, right prong 1×1.3.
        check(&p, &t, 6, 3.0 + 1.0 + 1.3);
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
        let t = triangulate_polygon(&p, &(0..6).collect::<Vec<_>>(), &[]).unwrap();
        check(&p, &t, 4, 3.0);
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
