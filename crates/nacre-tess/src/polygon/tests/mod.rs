use super::*;
use std::collections::HashSet;

/// The goldens below are written **as a chart** — the numbers are `uv` already.
///
/// ★ No projection is involved: the projection rule lives once, in `planar_chart`, and these
/// rings are already in its output space.
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
    triangulate_uv(&mut v, &refs, &[], &[]).map(|(t, _)| t)
}

/// The same road with one bridge named — what the chart layer hands over for a pinched face.
fn tri_bridged(
    uv: &[P2],
    outer: &[usize],
    holes: &[&[usize]],
    bridge: Bridge,
) -> Result<Triangulated, TessError> {
    let mut v = uv.to_vec();
    let refs: Vec<&[usize]> = std::iter::once(outer)
        .chain(holes.iter().copied())
        .collect();
    triangulate_uv(&mut v, &refs, &[], &[bridge])
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
    // The L-prism cap: star-shaped from (0,0), so even a fan from vertex 0 is right —
    // which is why a fan hides its defect on this shape.
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

/// **The wall with two windows — a case ear clipping finds no ear in.**
///
/// Both holes' first mutually visible ring vertex was the *same* outer corner, so
/// the bridged ring visited it three times and the polygon was pinched there. A
/// pinched ring is not simple, Meisters' two-ears theorem does not apply to it, and
/// measurement agreed: all six convex vertices had a diagonal that genuinely crossed
/// the other hole. The sweep never merges the rings, so the pinch is not built.
///
/// Everything here is transcribed from the failing model (a hub wall with two fins
/// through it, projected to its own plane) — coordinates, index layout, and the
/// order the two inner loops arrive in, which decides where each bridge lands.
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

/// **The corpus on which ear clipping fails 27.6% of the time** (measured), as a
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
    let t = triangulate_uv(&mut uv, &[&ring], &lattice, &[]).unwrap().0;
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
        let t = triangulate_uv(&mut uv, &[&ring], &[p], &[]).unwrap().0;
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
/// no answer for *as it stands* (an exact tangency, which the sampled boundary reproduces when
/// a sample lands on it); given a bridge it is drawn (see `a_bridged_touch_meshes_the_pinched_face`),
/// and without one, as here, it is named. A **crossing** is a ring set that has no
/// triangulation at all: any one would be wrong. Different claims about *the mesh*, different
/// names — neither is a verdict on the solid, which is `validate`'s to give.
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

/// ★★★★★ **A touch vouches for nothing — two ways a check that trusts one lets a crossing
/// through.**
///
/// Returning on the *first* touch found, and scanning for crossings only when there is none, has
/// two defects; this fixture is the first, the next test the second.
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

/// ★★★★★ **A pinched face meshes once its touch is bridged.** The square's bottom edge
/// carries the hole's apex as a vertex of its own (index 1 — what the chart layer's pre-pass
/// does to the shared edge), the hole's apex is index 5 at the same point, and the bridge
/// names the two. The sweep then sees one ring, orders the twins by where they slide, and
/// produces a triangulation that covers exactly the square minus the hole — six triangles
/// for the eight-vertex ring, all CCW, whose once-used edges are the merged ring.
#[test]
fn a_bridged_touch_meshes_the_pinched_face() {
    let p = pts(&[
        [0.0, 0.0],
        [2.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [2.0, 0.0],
        [1.0, 2.0],
        [3.0, 2.0],
    ]);
    let (tris, rings) = tri_bridged(
        &p,
        &[0, 1, 2, 3, 4],
        &[&[5, 6, 7]],
        Bridge {
            ring_x: 0,
            x: 1,
            ring_y: 1,
            y: 0,
        },
    )
    .expect("the bridged face meshes");
    check(&p, &tris, 6, 16.0 - 2.0);
    assert_eq!(rings.len(), 1, "the hole was spliced into the outer ring");
    check_partition(&tris, &rings[0], &[]);
    // Naming the rings the other way round makes no difference: which is the split ring is
    // read off the geometry.
    let (tris2, _) = tri_bridged(
        &p,
        &[0, 1, 2, 3, 4],
        &[&[5, 6, 7]],
        Bridge {
            ring_x: 1,
            x: 0,
            ring_y: 0,
            y: 1,
        },
    )
    .expect("the bridged face meshes");
    check(&p, &tris2, 6, 16.0 - 2.0);
    // And without the bridge the same rings are what they always were: a touch.
    let r = tri(&p, &[0, 1, 2, 3, 4], &[&[5, 6, 7]]);
    assert!(matches!(r, Err(TessError::SelfTouchingBoundary)), "{r:?}");
}

/// Two holes touching: the same bridge, and the merged ring is still a hole — CW, so the
/// winding gate on the rings as given still runs, and the outer ring is untouched.
#[test]
fn two_holes_that_touch_are_bridged_into_one_hole() {
    let p = pts(&[
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        // hole A, CW, with the touching point (3,2) as a vertex of its right edge
        [1.0, 1.0],
        [1.0, 3.0],
        [3.0, 3.0],
        [3.0, 2.0],
        [3.0, 1.0],
        // hole B, CW, its apex on A's right edge
        [3.0, 2.0],
        [3.5, 2.5],
        [3.5, 1.5],
    ]);
    let (tris, rings) = tri_bridged(
        &p,
        &[0, 1, 2, 3],
        &[&[4, 5, 6, 7, 8], &[9, 10, 11]],
        Bridge {
            ring_x: 1,
            x: 3,
            ring_y: 2,
            y: 0,
        },
    )
    .expect("the bridged holes mesh");
    check(&p, &tris, 12, 16.0 - 4.0 - 0.25);
    assert_eq!(rings.len(), 2, "outer, and one merged hole");
    assert_eq!(rings[0], vec![0, 1, 2, 3], "the outer ring is untouched");
    assert_eq!(
        ring_orientation(&rings[1], &p),
        -1,
        "the merged hole is still a hole"
    );
    check_partition(&tris, &rings[0], &[&rings[1]]);
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
