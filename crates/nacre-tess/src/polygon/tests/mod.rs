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
    triangulate_uv(&mut v, &refs, &[], &[Revisit::Bridge(bridge)])
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

/// The same road with one ring passing a point twice, at positions `x < y`.
fn tri_pinched(uv: &[P2], ring: &[usize], x: usize, y: usize) -> Result<Triangulated, TessError> {
    let mut v = uv.to_vec();
    triangulate_uv(&mut v, &[ring], &[], &[Revisit::Pinch { ring: 0, x, y }])
}

/// ★★★★★ **A pinched ring meshes without a bridge.** The square's bottom is touched at `(2, 0)`
/// by the tip of a triangular groove, so the one ring passes that point twice (indices 0 and 5).
/// Six triangles for the eight-vertex ring, all CCW, covering exactly the square less the groove,
/// whose once-used edges are the ring.
///
/// Both spellings of the ring are asked, because they slide the twins differently. Starting at
/// the point, the earlier visit slides up the groove's right side and the later up its left — a
/// slide with a `v` component while `(0, 0)` and `(4, 0)` sit at the twins' own `v`, the tie the
/// sweep breaks by `u` (`sos`: the sweep's tie-break is a turn `δ`, the twins' slide `ε ≪ δ`).
/// Starting one vertex on, both slide along the bottom.
#[test]
fn a_pinched_ring_meshes_without_a_bridge() {
    let p = pts(&[
        [2.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [0.0, 0.0],
        [2.0, 0.0],
        [1.5, 1.0],
        [2.5, 1.0],
    ]);
    let ring: Vec<usize> = (0..8).collect();
    let (tris, rings) = tri_pinched(&p, &ring, 0, 5).expect("the pinched ring meshes");
    check(&p, &tris, 6, 16.0 - 0.5);
    assert_eq!(rings, vec![ring.clone()], "nothing spliced");
    check_partition(&tris, &ring, &[]);
    let mut turned = ring.clone();
    turned.rotate_left(1);
    let (tris, _) = tri_pinched(&p, &turned, 4, 7).expect("the turned ring meshes");
    check(&p, &tris, 6, 16.0 - 0.5);
    check_partition(&tris, &turned, &[]);
}

/// ★★ **A ring pinched at two points meshes — one pair of twins each.** A `6 × 4` rectangle whose
/// bottom is touched by the tips of two triangular grooves, at `(2, 0)` and `(4, 0)`: the one ring
/// passes each point twice (indices 1 and 4, 5 and 8). Ten triangles for the twelve-vertex ring,
/// all CCW, covering exactly the rectangle less the grooves, whose once-used edges are the ring.
#[test]
fn a_ring_pinched_at_two_points_meshes() {
    let p = pts(&[
        [0.0, 0.0],
        [2.0, 0.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [2.0, 0.0],
        [4.0, 0.0],
        [3.5, 1.0],
        [4.5, 1.0],
        [4.0, 0.0],
        [6.0, 0.0],
        [6.0, 4.0],
        [0.0, 4.0],
    ]);
    let ring: Vec<usize> = (0..12).collect();
    let mut v = p.clone();
    let revisits = [
        Revisit::Pinch {
            ring: 0,
            x: 1,
            y: 4,
        },
        Revisit::Pinch {
            ring: 0,
            x: 5,
            y: 8,
        },
    ];
    let (tris, rings) =
        triangulate_uv(&mut v, &[&ring], &[], &revisits).expect("the twice-pinched ring meshes");
    check(&p, &tris, 10, 24.0 - 0.5 - 0.5);
    assert_eq!(rings, vec![ring.clone()], "nothing spliced");
    check_partition(&tris, &ring, &[]);
}

/// ★★★★★ **A pinch whose sector is not convex meshes — that visit stays, the other slides.** The
/// top of a `6 × 6` plate with a triangular hole `(2,2) (4,3) (4,4)` and the footprint of a
/// triangular boss `(2,2) (4,0.5) (4,1)` meeting it at `(2, 2)`: one hole ring passing that point
/// twice (indices 5 and 8), as the kernel builds it. One visit's sector is the 53° between the two
/// triangles; the other's is the 278° around them, and it stays put. Ten triangles, the plate
/// less the two triangles, whose once-used edges are the rings.
#[test]
fn a_pinch_with_a_reflex_sector_meshes() {
    let p = pts(&[
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, 6.0],
        [0.0, 6.0],
        [4.0, 3.0],
        [2.0, 2.0],
        [4.0, 1.0],
        [4.0, 0.5],
        [2.0, 2.0],
        [4.0, 4.0],
    ]);
    let (outer, hole) = ([0, 1, 2, 3], [4, 5, 6, 7, 8, 9]);
    let mut v = p.clone();
    let revisits = [Revisit::Pinch {
        ring: 1,
        x: 1,
        y: 4,
    }];
    let (tris, rings) =
        triangulate_uv(&mut v, &[&outer, &hole], &[], &revisits).expect("the pinched hole meshes");
    check(&p, &tris, 10, 36.0 - 1.0 - 0.5);
    check_partition(&tris, &rings[0], &[&rings[1]]);
}

/// ★★ **Two pinches on one ring, one of them with a reflex sector.** The same plate with the
/// boss's corner on the plate's edge, `(2,2) (4,0) (4,1)`: the outer ring passes `(4, 0)` twice
/// (both sectors convex, both slide) and `(2, 2)` twice (one sector 270°, which stays) — the top
/// face the kernel builds for that boss. Nine triangles for the eleven-vertex ring.
#[test]
fn a_ring_pinched_twice_with_a_reflex_sector_meshes() {
    let p = pts(&[
        [0.0, 0.0],
        [4.0, 0.0],
        [2.0, 2.0],
        [4.0, 4.0],
        [4.0, 3.0],
        [2.0, 2.0],
        [4.0, 1.0],
        [4.0, 0.0],
        [6.0, 0.0],
        [6.0, 6.0],
        [0.0, 6.0],
    ]);
    let ring: Vec<usize> = (0..11).collect();
    let mut v = p.clone();
    let revisits = [
        Revisit::Pinch {
            ring: 0,
            x: 2,
            y: 5,
        },
        Revisit::Pinch {
            ring: 0,
            x: 1,
            y: 7,
        },
    ];
    let (tris, rings) =
        triangulate_uv(&mut v, &[&ring], &[], &revisits).expect("the twice-pinched ring meshes");
    check(&p, &tris, 9, 36.0 - 1.0 - 1.0);
    assert_eq!(rings, vec![ring.clone()], "nothing spliced");
    check_partition(&tris, &ring, &[]);
}

/// ★★ **Two lobes meeting only at a point are not a face.** Two counter-clockwise triangles above
/// `T = (0, 0)` (rays at 11°/22° and 39°/50°), one ring through `T` twice: each visit's wedge holds
/// both of the other's rays. The visit crossing from the upper lobe to the lower turns left — a
/// convex wedge — so «slide the convex one» alone would move it into the other lobe. That ring is
/// not a face, and `pinch` says so; read as a valid pinch instead, the sweep refuses it as
/// `DegenerateRing` (measured) — the shape's name is what this answer decides, not whether a
/// wrong mesh comes out.
#[test]
fn two_lobes_meeting_at_a_point_are_not_a_face() {
    let p = pts(&[
        [0.0, 0.0],
        [10.0, 2.0],
        [10.0, 4.0],
        [0.0, 0.0],
        [10.0, 8.0],
        [10.0, 12.0],
    ]);
    let ring: Vec<usize> = (0..6).collect();
    let r = tri_pinched(&p, &ring, 0, 3);
    assert!(matches!(r, Err(TessError::SelfTouchingBoundary)), "{r:?}");
}

/// What a pinch the sweep cannot order comes back as. The boundary passing *through* the point —
/// one visit's ray inside the other's wedge and the other outside — is a crossing every
/// triangulation gets wrong (`DegenerateRing`); without `pinch`'s check the sweep is handed it and
/// builds a piece whose chains do not descend. Two lobes meeting only at the point — each visit's
/// wedge holding both of the other's rays — are not a face's interior, and a ray of one visit
/// along the other's is a touch: both `SelfTouchingBoundary`.
#[test]
fn a_pinch_the_sweep_cannot_order_is_named() {
    // The later visit leaves the point into the earlier one's wedge.
    let through = pts(&[
        [2.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [0.0, 0.0],
        [2.0, 0.0],
        [3.7, 1.5],
        [1.5, 0.1],
    ]);
    let ring: Vec<usize> = (0..8).collect();
    let r = tri_pinched(&through, &ring, 0, 5);
    assert!(matches!(r, Err(TessError::DegenerateRing)), "{r:?}");
    // An L whose reflex corner a triangle touches from outside: the L and the triangle meet only
    // there, and each visit's wedge holds both of the other's rays — two lobes.
    let reflex = pts(&[
        [0.0, 0.0],
        [2.0, 1.0],
        [1.0, 2.0],
        [0.0, 0.0],
        [0.0, 2.0],
        [-2.0, 2.0],
        [-2.0, -2.0],
        [2.0, -2.0],
        [2.0, 0.0],
    ]);
    let ring: Vec<usize> = (0..9).collect();
    let r = tri_pinched(&reflex, &ring, 0, 3);
    assert!(matches!(r, Err(TessError::SelfTouchingBoundary)), "{r:?}");
    // The later visit leaves along the line the earlier one arrived on.
    let along = pts(&[
        [2.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 2.0],
        [3.0, 3.0],
        [1.5, 1.0],
    ]);
    let ring: Vec<usize> = (0..9).collect();
    let r = tri_pinched(&along, &ring, 0, 5);
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

/// ★★ **A ring pinched at its extreme point winds by the first ray past the outside gap.** Three
/// rings, each passing `T = (0, 0)` — its lowest point — twice, each asked from every starting
/// index and reversed, so neither where the ring starts nor which pass comes first decides it:
///
/// - a clockwise figure-eight hole, two lobes above `T`: one pass turns left there (its own sector
///   is the 37° between the lobes), the other right — `−1`;
/// - a counter-clockwise outer ring holding a triangular hole that touches it at `T`: `+1`;
/// - two counter-clockwise lobes meeting only at `T` (rays at 11°/22° and 39°/50°): the pass that
///   crosses from one lobe to the other turns right, and the ring still winds left — `+1`.
///
/// Reading only the first pass gets the hole wrong from half its starting indices; reading «any
/// pass turning right means clockwise» gets the two lobes wrong.
#[test]
fn a_ring_pinched_at_its_extreme_point_winds_by_its_first_ray() {
    let figure_eight_hole = [
        [0.0, 0.0],
        [-2.0, 3.0],
        [-1.0, 3.0],
        [0.0, 0.0],
        [1.0, 3.0],
        [2.0, 3.0],
    ];
    let outer_holding_a_hole = [
        [0.0, 0.0],
        [4.0, 4.0],
        [-4.0, 4.0],
        [0.0, 0.0],
        [-1.0, 2.0],
        [1.0, 2.0],
    ];
    let two_lobes = [
        [0.0, 0.0],
        [10.0, 2.0],
        [10.0, 4.0],
        [0.0, 0.0],
        [10.0, 8.0],
        [10.0, 12.0],
    ];
    for (name, ring, want) in [
        ("figure-eight hole", figure_eight_hole, -1i8),
        ("outer ring holding a hole", outer_holding_a_hole, 1),
        ("two lobes", two_lobes, 1),
    ] {
        let p = pts(&ring);
        for reversed in [false, true] {
            for start in 0..ring.len() {
                let mut r: Vec<usize> = (0..ring.len()).collect();
                if reversed {
                    r.reverse();
                }
                r.rotate_left(start);
                let want = if reversed { -want } else { want };
                assert_eq!(
                    ring_orientation(&r, &p),
                    want,
                    "{name}, reversed = {reversed}, start = {start}"
                );
            }
        }
    }
}

/// A deterministic draw in `lo..=hi` — the generator `random_rectilinear_faces_triangulate` uses.
fn lcg(seed: u64) -> impl FnMut(i64, i64) -> i64 {
    let mut st = seed;
    move |lo: i64, hi: i64| {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + ((st >> 33) % ((hi - lo + 1) as u64)) as i64
    }
}

/// Twice the signed area of a ring — exact here: the corpora below are small integers.
fn area2(p: &[P2], ring: &[usize]) -> f64 {
    (0..ring.len())
        .map(|k| {
            let (a, b) = (p[ring[k]], p[ring[(k + 1) % ring.len()]]);
            a[0] * b[1] - a[1] * b[0]
        })
        .sum()
}

/// ★★★★ **Figure-eight holes at random angles — every one meshes.** A `20 × 20` square with one
/// hole made of two triangular lobes meeting at `T`, their four rays drawn at random integer
/// directions and the lobes taken as two runs of adjacent rays, so the two visits never
/// interleave. Each visit bounds the sector of material between the lobes on its side; in about
/// a quarter of the samples one of those is reflex (27.5% of 119,016 in the measurement this
/// harness was drawn from), and it stays while the other slides. Every sample must come
/// back a genuine triangulation of the square less the two lobes.
#[test]
fn random_figure_eight_holes_triangulate() {
    let mut rng = lcg(0xD1B5_4A32_D192_ED03);
    let cross = |a: (i64, i64), b: (i64, i64)| a.0 * b.1 - a.1 * b.0;
    let (mut samples, mut reflex) = (0, 0);
    for _ in 0..20_000 {
        let t = (rng(9, 11), rng(9, 11));
        let mut dirs: Vec<(i64, i64)> = (0..4).map(|_| (rng(-4, 4), rng(-4, 4))).collect();
        // No zero direction, and no two rays along one direction (their edges would overlap).
        let same_way = |a: (i64, i64), b: (i64, i64)| cross(a, b) == 0 && a.0 * b.0 + a.1 * b.1 > 0;
        if dirs.contains(&(0, 0)) || (0..4).any(|i| (i + 1..4).any(|j| same_way(dirs[i], dirs[j])))
        {
            continue;
        }
        let half = |d: &(i64, i64)| u8::from(!(d.1 > 0 || (d.1 == 0 && d.0 > 0)));
        dirs.sort_by(|a, b| half(a).cmp(&half(b)).then(cross(*b, *a).cmp(&0)));
        let rot = rng(0, 1) as usize;
        let (l1, l2) = (
            (dirs[rot], dirs[rot + 1]),
            (dirs[rot + 2], dirs[(rot + 3) % 4]),
        );
        // Each lobe strictly less than a half-turn, so it is a triangle.
        if cross(l1.0, l1.1) <= 0 || cross(l2.0, l2.1) <= 0 {
            continue;
        }
        let at = |d: (i64, i64), k: i64| [(t.0 + k * d.0) as f64, (t.1 + k * d.1) as f64];
        let tp = [t.0 as f64, t.1 as f64];
        let (k1, k2, k3, k4) = (rng(1, 2), rng(1, 2), rng(1, 2), rng(1, 2));
        // The outer ring, then the hole clockwise: `T`, one lobe, `T`, the other.
        let p = pts(&[
            [0.0, 0.0],
            [20.0, 0.0],
            [20.0, 20.0],
            [0.0, 20.0],
            tp,
            at(l1.1, k2),
            at(l1.0, k1),
            tp,
            at(l2.1, k4),
            at(l2.0, k3),
        ]);
        let (outer, hole) = ([0, 1, 2, 3], [4, 5, 6, 7, 8, 9]);
        let lobes = (cross(l1.0, l1.1) * k1 * k2 + cross(l2.0, l2.1) * k3 * k4) as f64 / 2.0;
        let mut v = p.clone();
        let (tris, rings) = triangulate_uv(
            &mut v,
            &[&outer, &hole],
            &[],
            &[Revisit::Pinch {
                ring: 1,
                x: 0,
                y: 3,
            }],
        )
        .unwrap_or_else(|e| panic!("{e:?}: {p:?}"));
        check(&p, &tris, 10, 400.0 - lobes);
        check_partition(&tris, &rings[0], &[&rings[1]]);
        samples += 1;
        if orient2d(p[9], tp, p[5]) < 0.0 || orient2d(p[6], tp, p[8]) < 0.0 {
            reflex += 1;
        }
    }
    assert!(samples > 10_000, "only {samples} samples reached the gate");
    assert!(reflex > 2_500, "only {reflex} samples had a reflex sector");
}

/// The boundary rings of a union of unit cells on an `n × n` grid — material on the left of every
/// edge — with the pairing at each corner where two cells touch only diagonally chosen at random
/// (turn left or right there), and the number of edge-connected components of the cells.
fn grid_rings(rng: &mut impl FnMut(i64, i64) -> i64, n: i64) -> (Vec<Vec<(i64, i64)>>, usize) {
    use std::collections::{BTreeMap, BTreeSet};
    let fill = rng(30, 80);
    let cells: BTreeSet<(i64, i64)> = (0..n)
        .flat_map(|x| (0..n).map(move |y| (x, y)))
        .filter(|_| rng(1, 100) <= fill)
        .collect();
    let mut sides = BTreeSet::new();
    for &(x, y) in &cells {
        let c = [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)];
        for k in 0..4 {
            sides.insert((c[k], c[(k + 1) % 4]));
        }
    }
    let kept: Vec<_> = sides
        .iter()
        .filter(|&&(a, b)| !sides.contains(&(b, a)))
        .copied()
        .collect();
    let mut out: BTreeMap<(i64, i64), Vec<(i64, i64)>> = BTreeMap::new();
    for &(a, b) in &kept {
        out.entry(a).or_default().push(b);
    }
    let mut turn_left: BTreeMap<(i64, i64), bool> = BTreeMap::new();
    let mut used = BTreeSet::new();
    let mut rings = Vec::new();
    for &(a0, b0) in &kept {
        if used.contains(&(a0, b0)) {
            continue;
        }
        let (mut a, mut b) = (a0, b0);
        let mut ring = vec![a0];
        loop {
            used.insert((a, b));
            let outs = &out[&b];
            let next = if outs.len() == 1 {
                outs[0]
            } else {
                let left = *turn_left.entry(b).or_insert_with(|| rng(0, 1) == 0);
                let d = (b.0 - a.0, b.1 - a.1);
                let turn = |o: &(i64, i64)| d.0 * (o.1 - b.1) - d.1 * (o.0 - b.0);
                *outs
                    .iter()
                    .find(|o| if left { turn(o) > 0 } else { turn(o) < 0 })
                    .expect("a diagonal touch leaves one way to each side")
            };
            if (b, next) == (a0, b0) {
                break;
            }
            ring.push(b);
            (a, b) = (b, next);
        }
        // Straight corners dissolved, as a b-rep face's loops come.
        let m = ring.len();
        rings.push(
            (0..m)
                .filter(|&i| {
                    let (p, c, q) = (ring[(i + m - 1) % m], ring[i], ring[(i + 1) % m]);
                    (c.0 - p.0) * (q.1 - c.1) - (c.1 - p.1) * (q.0 - c.0) != 0
                })
                .map(|i| ring[i])
                .collect(),
        );
    }
    // Edge-connected components of the cells.
    let mut seen = BTreeSet::new();
    let mut components = 0;
    for &c in &cells {
        if !seen.insert(c) {
            continue;
        }
        components += 1;
        let mut stack = vec![c];
        while let Some((x, y)) = stack.pop() {
            for nb in [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)] {
                if cells.contains(&nb) && seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
    }
    (rings, components)
}

/// How many random grid unions [`random_pinched_grid_faces_triangulate`] draws: at 200,000 the
/// cell oracle agreed with every answer (3,598 pinched faces meshed, 19,393 lobe rings and 6,884
/// corner bridges refused); a tenth keeps the commit hook quick.
const GRID_SAMPLES: usize = 20_000;

/// ★★★★ **Unions of grid cells, sheared — every pinched face meshes, and only a face is meshed.**
/// Random unions of unit cells on a small grid, mapped by two integer shears in turn (determinant
/// one, so orientation and every pinch survive, at angles other than right ones). Where two cells
/// touch only at a corner the boundary passes that corner twice, and the pairing there is drawn
/// at random. The oracle is the cells, not the rings: one outer ring over **edge-connected** cells
/// is one face pinched where cells touch at corners — it must mesh, a genuine triangulation of
/// the cells; one ring over cells in several pieces is lobes meeting only at points —
/// `SelfTouchingBoundary`; and two rings sharing a corner (`AtEnd`, no straight edge to splice
/// into) stay what they were, `SelfTouchingBoundary`.
#[test]
fn random_pinched_grid_faces_triangulate() {
    let mut rng = lcg(0x9E37_79B9_7F4A_7C15);
    let (mut pinched, mut lobes, mut corner_bridges) = (0, 0, 0);
    for _ in 0..GRID_SAMPLES {
        let n = rng(2, 5);
        let (rings, components) = grid_rings(&mut rng, n);
        let (k, m2) = (rng(-2, 2), rng(-1, 1));
        let map = |q: (i64, i64)| {
            let x = q.0 + k * q.1;
            [x as f64, (q.1 + m2 * x) as f64]
        };
        let mut p: Vec<P2> = Vec::new();
        let mut idx: Vec<Vec<usize>> = Vec::new();
        for r in &rings {
            idx.push((p.len()..p.len() + r.len()).collect());
            p.extend(r.iter().map(|&q| map(q)));
        }
        let outer: Vec<usize> = (0..idx.len())
            .filter(|&i| area2(&p, &idx[i]) > 0.0)
            .collect();
        if outer.len() != 1 {
            continue;
        }
        let order: Vec<usize> = std::iter::once(outer[0])
            .chain((0..idx.len()).filter(|&i| i != outer[0]))
            .collect();
        let idx: Vec<Vec<usize>> = order.iter().map(|&i| idx[i].clone()).collect();
        // Revisits by coordinate, as `shared_vertices` finds them by handle.
        let mut first: std::collections::BTreeMap<(u64, u64), (usize, usize)> =
            std::collections::BTreeMap::new();
        let mut revisits = Vec::new();
        for (ri, ring) in idx.iter().enumerate() {
            for (kk, &i) in ring.iter().enumerate() {
                match first.get(&(p[i][0].to_bits(), p[i][1].to_bits())) {
                    Some(&(rj, kj)) if rj != ri => revisits.push(Revisit::Bridge(Bridge {
                        ring_x: rj,
                        x: kj,
                        ring_y: ri,
                        y: kk,
                    })),
                    Some(&(_, kj)) => revisits.push(Revisit::Pinch {
                        ring: ri,
                        x: kj,
                        y: kk,
                    }),
                    None => {
                        first.insert((p[i][0].to_bits(), p[i][1].to_bits()), (ri, kk));
                    }
                }
            }
        }
        if revisits.is_empty() {
            continue;
        }
        let refs: Vec<&[usize]> = idx.iter().map(|r| r.as_slice()).collect();
        let mut v = p.clone();
        let r = triangulate_uv(&mut v, &refs, &[], &revisits);
        if revisits.iter().any(|r| matches!(r, Revisit::Bridge(_))) {
            assert!(
                matches!(r, Err(TessError::SelfTouchingBoundary)),
                "{r:?}: {p:?} {idx:?}"
            );
            corner_bridges += 1;
        } else if components > 1 {
            assert!(
                matches!(r, Err(TessError::SelfTouchingBoundary)),
                "{r:?}: {p:?} {idx:?}"
            );
            lobes += 1;
        } else {
            let (tris, used) = r.unwrap_or_else(|e| panic!("{e:?}: {p:?} {idx:?}"));
            let area: f64 = idx.iter().map(|r| area2(&p, r)).sum::<f64>() / 2.0;
            let holes: Vec<&[usize]> = idx[1..].iter().map(|r| r.as_slice()).collect();
            check(&p, &tris, p.len() + 2 * holes.len() - 2, area);
            let used_holes: Vec<&[usize]> = used[1..].iter().map(|r| r.as_slice()).collect();
            check_partition(&tris, &used[0], &used_holes);
            pinched += 1;
        }
    }
    assert!(corner_bridges > 0, "no corner bridge reached the gate");
    assert!(
        pinched > GRID_SAMPLES / 100,
        "only {pinched} pinched faces reached the gate"
    );
    assert!(
        lobes > GRID_SAMPLES / 100,
        "only {lobes} lobe rings reached the gate"
    );
}
