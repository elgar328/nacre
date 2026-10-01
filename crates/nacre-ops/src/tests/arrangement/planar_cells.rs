//! The planar arrangement: edge masks, chords, angular order, splitting, cell labels, keep rules,
//! end-to-end planar booleans.

use super::*;

/// ★ **A fixture with no cylinders, said as a fact rather than left as a hole.** The table is
/// how a `NodeId::Pierce` reaches its definition, so an all-plane fixture has nothing to put
/// in it — and a bare `&[]` at a call site reads like something forgotten.
const NO_CYLS: &[crate::planes::WorkingCyl] = &[];

/// The graze rows of `edge_mask`'s algebra, one by one. The missing row was the
/// same-side pair: a knife edge's two faces both leave `W` upward, the material
/// between them pinches to measure zero at the plane, and flipping anything there
/// hands the label propagation a contradiction it reports two layers later.
#[test]
fn edge_mask_graze_algebra() {
    let g = |above: bool| (SolidSide::A, SegKind::Graze { body_above: above });
    let flips = |m: &[(SolidSide, SegKind)]| edge_mask(m).unwrap();

    // One graze: a step — the face fills that arc, its bit flips.
    assert_eq!(flips(&[g(true)]), [true, false, false, false]);
    assert_eq!(flips(&[g(false)]), [false, true, false, false]);
    // An opposite pair: the solid's edge lies in W with material above on one
    // in-plane side and below on the other — both flip (the [T,T] picture).
    assert_eq!(flips(&[g(true), g(false)]), [true, true, false, false]);
    // ★ A same-side pair: knife edge (or its reflex complement) — nothing flips,
    // whichever side the pair is on.
    assert_eq!(flips(&[g(true), g(true)]), [false; 4]);
    assert_eq!(flips(&[g(false), g(false)]), [false; 4]);
    // The two solids are independent lanes.
    let gb = (SolidSide::B, SegKind::Graze { body_above: true });
    assert_eq!(flips(&[g(true), g(true), gb]), [false, false, true, false]);

    // More than two grazes from one solid is degenerate input, not a case.
    assert!(edge_mask(&[g(true), g(true), g(true)]).is_err());
    // Unchanged refusals: a graze coincident with a same-solid true crossing.
    let t = (SolidSide::A, SegKind::Transversal { mat: 1 });
    assert!(edge_mask(&[g(true), t]).is_err());
}

/// Point of a named vertex, for asserting geometry by hand.
fn pt(n: NodeId, jd: &Judge<'_, WorkingPlane>) -> [f64; 3] {
    let t = three_plane_name(n).expect("a three-plane node");
    let planes = jd.planes;
    three_planes(
        &planes[t[0]].plane,
        &planes[t[1]].plane,
        &planes[t[2]].plane,
    )
    .unwrap()
    .as_array()
}

/// A seated face's whole boundary is emitted, with `body_above` matching the geometry.
///
/// Two stacked cubes sharing the plane z=1: `a=[0,1]³`, `b=[0,1]²×[1,2]`. On the shared
/// class, `a`'s top cap is seated with its body **below** (body_above=false) and `b`'s bottom
/// cap is seated with its body **above** (body_above=true). Each emits its 4 boundary edges.
#[test]
fn a_seated_cap_emits_its_boundary_with_the_right_body_side() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);

    // The shared class z=1: the class both solids seat a cap on.
    let wc = (0..planes.len())
        .find(|&c| {
            let seats = |s: Handle<Solid>| {
                solid_shell_handles(&m, s).into_iter().any(|sh| {
                    m.shell(sh).faces.iter().any(|fh| {
                        plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, &surf_ix, &faces_tab)
                    })
                })
            };
            seats(a) && seats(b)
        })
        .expect("a shared cap class");

    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    // Two seated caps × 4 edges each. (The side walls only *graze* z=1 — each cube's body is
    // entirely on one side — so they trace as `Graze` chords, not `Transversal`; the two far
    // caps are parallel to z=1 and contribute nothing — misses, not declines.)
    let seated: Vec<&Seg> = tr
        .segs
        .iter()
        .filter(|s| matches!(s.kind, SegKind::Seated { .. }))
        .collect();
    assert_eq!(seated.len(), 8, "two seated caps, four edges each: {tr:?}");
    // Four side walls per cube each graze z=1 in one chord, body on the cube's side.
    let grazes: Vec<&Seg> = tr
        .segs
        .iter()
        .filter(|s| matches!(s.kind, SegKind::Graze { .. }))
        .collect();
    assert_eq!(grazes.len(), 8, "eight side-wall graze chords: {tr:?}");
    assert_eq!(
        tr.segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .count(),
        0,
        "no wall crosses z=1: {tr:?}"
    );
    for s in &grazes {
        match (s.solid, s.kind) {
            (SolidSide::A, SegKind::Graze { body_above }) => assert!(!body_above, "a below"),
            (SolidSide::B, SegKind::Graze { body_above }) => assert!(body_above, "b above"),
            _ => unreachable!(),
        }
    }
    assert!(
        tr.declined.is_empty(),
        "axis-aligned cubes decline nothing: {tr:?}"
    );
    assert!(tr.touches.is_empty());
    for s in &seated {
        // On z=1, riding a side wall (not W itself), body on the geometrically correct side.
        for e in &s.end {
            assert!((pt(*e, &jd)[2] - 1.0).abs() < 1e-12, "endpoint on z=1");
        }
        assert_ne!(s.wall, wc, "a seated edge rides a side wall, not W");
        match (s.solid, s.kind) {
            (SolidSide::A, SegKind::Seated { body_above }) => assert!(!body_above, "a below"),
            (SolidSide::B, SegKind::Seated { body_above }) => assert!(body_above, "b above"),
            _ => unreachable!(),
        }
    }
}

/// Round a point for set membership by hand.
fn rp(p: [f64; 3]) -> [i64; 3] {
    [
        (p[0] * 1e6).round() as i64,
        (p[1] * 1e6).round() as i64,
        (p[2] * 1e6).round() as i64,
    ]
}

/// The chord across a non-convex reflex plane — the cap chord — which the seated brick provably
/// cannot produce (it declines this face). l_prism (profile (0,0),(2,0),(2,1),(1,1),(1,2),(0,2),
/// extruded z∈[0,1]) cut at y=1: the z=0 cap covers x∈[0,2], as a **transversal** stretch
/// x∈[0,1] where `y=1` runs through the cap's interior plus a **graze** x∈[1,2] where the cap's
/// own boundary edge rides the line and the cap lies on one side of it.
#[test]
fn the_cap_chord_stops_where_the_on_line_edge_begins() {
    let profile = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 0.0]),
        Point2::from_array([2.0, 1.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([1.0, 2.0]),
        Point2::from_array([0.0, 2.0]),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile,
        dist: 1.0,
    }])
    .unwrap();
    // A far cube so plane_index_setup has two solids; it never touches the y=1 class.
    crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([10.0, 10.0, 10.0]),
        Point3::from_array([11.0, 11.0, 11.0]),
    );
    m.rebuild_adjacency();
    let a2 = m.live_solids()[0];
    let b2 = m.live_solids()[1];
    let PlaneSetup {
        cyls,
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a2, b2).unwrap();
    let jd = Judge::new(&planes, standard, &notes);

    // The y=1 class: a plane through all-y=1 points that a's reflex face sits on.
    let wc = (0..planes.len())
        .find(|&c| {
            planes[c]
                .witness_coords()
                .iter()
                .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("a y=1 class");

    let mut out = Trace::default();
    trace_one_of(
        &m,
        a2,
        SolidSide::A,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_a,
        &plane_ix,
        Default::default(),
        &mut out,
    );

    // Collect transversal segments as rounded endpoint-pairs.
    let chords: Vec<([i64; 3], [i64; 3])> = out
        .segs
        .iter()
        .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
        .map(|s| {
            let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            (u, v)
        })
        .collect();

    // The cap chords: the z=0 (and z=1) cap ∩ (y=1) covers x∈[0,2]. A cap is a horizontal
    // segment (both endpoints share z), which distinguishes it from the reflex wall x=1's own
    // vertical chord (1,1,0)-(1,1,1), a legitimate different face's trace.
    let horiz_z0: Vec<_> = chords
        .iter()
        .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
        .collect();
    // ★ Falsifiable core: the on-line edge does not merge into the chord that runs through the
    // cap's interior. x∈[0,1] is interior (the cap straddles y=1) and stays transversal;
    // x∈[1,2] is the cap's own boundary edge, so the cap is on one side there and it leaves as
    // a graze — see `run_body_above`. A single chord spanning [0,2] would be occupancy-blind
    // bridging.
    assert_eq!(
        horiz_z0.len(),
        1,
        "one transversal chord, the interior stretch: {chords:?}"
    );
    assert_eq!(
        *horiz_z0[0],
        (rp([0.0, 1.0, 0.0]), rp([1.0, 1.0, 0.0])),
        "the interior chord is x∈[0,1]: {chords:?}"
    );
    // The seated brick provably cannot produce this: it declines the z=0 cap as `not-seated`
    // (its class is z=0, not y=1). So a transversal chord riding fp=z=0 is new capability.
    // mat is a per-face constant, never 0.
    for s in &out.segs {
        if let SegKind::Transversal { mat } = s.kind {
            assert!(mat == 1 || mat == -1, "mat is a clean side: {mat}");
        }
    }
}

/// The flanks-equal (tangential on-line edge) branch: u_prism's notch bottom on y=1. The run
/// x∈[1,2] is flanked by material on **both** sides (the two prongs), so it does not toggle
/// parity yet must still appear inside the trace.
///
/// The cap ∩ y=1 covers x∈[0,3], but **not as one kind**: over the notch bottom the line is on
/// the cap's *boundary* and the cap lies below it, while on either side the line runs through
/// the cap's *interior* and the cap straddles. So the trace is `T[0,1] · G[1,2] · T[2,3]`, and
/// the graze's side is **below** — note the run's flanks (2,2.3) and (1,2) are both *above*,
/// which is why reading the side off a flank gets a notch backwards.
#[test]
fn a_tangential_on_line_edge_spans_as_transversal_graze_transversal() {
    let profile = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([3.0, 0.0]),
        Point2::from_array([3.0, 2.3]),
        Point2::from_array([2.0, 2.3]),
        Point2::from_array([2.0, 1.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([1.0, 2.0]),
        Point2::from_array([0.0, 2.0]),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile,
        dist: 1.0,
    }])
    .unwrap();
    crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([10.0, 10.0, 10.0]),
        Point3::from_array([11.0, 11.0, 11.0]),
    );
    m.rebuild_adjacency();
    let a2 = m.live_solids()[0];
    let b2 = m.live_solids()[1];
    let PlaneSetup {
        cyls,
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a2, b2).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = (0..planes.len())
        .find(|&c| {
            planes[c]
                .witness_coords()
                .iter()
                .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("a y=1 class");
    let mut out = Trace::default();
    trace_one_of(
        &m,
        a2,
        SolidSide::A,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_a,
        &plane_ix,
        Default::default(),
        &mut out,
    );

    let mut horiz_z0: Vec<_> = out
        .segs
        .iter()
        .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
        .map(|s| {
            let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            (u, v)
        })
        .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
        .collect();
    horiz_z0.sort();
    assert_eq!(
        horiz_z0,
        vec![
            (rp([0.0, 1.0, 0.0]), rp([1.0, 1.0, 0.0])),
            (rp([2.0, 1.0, 0.0]), rp([3.0, 1.0, 0.0])),
        ],
        "the interior stretches straddle; the notch bottom is not one of them"
    );

    // ★ The run itself is covered, as a one-sided graze — coverage is not lost, the kind
    //   differs. Its body is **below** (y < 1), the side the notch's material is on.
    let graze_z0: Vec<_> = out
        .segs
        .iter()
        .filter_map(|s| match s.kind {
            SegKind::Graze { body_above } => {
                let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
                    .then_some((u, v, body_above))
            }
            _ => None,
        })
        .collect();
    assert_eq!(
        graze_z0,
        vec![(rp([1.0, 1.0, 0.0]), rp([2.0, 1.0, 0.0]), false)],
        "the notch bottom x∈[1,2], body below"
    );
}

/// ★ Spike: a coordinate-free exact cyclic order of edges around an arrangement vertex is
/// buildable — the "DNA question" the winding-based engine design flagged as a possible death
/// condition for the per-plane arrangement. A degree-5 vertex with directions +x, +y, −x, −y and one **oblique**
/// (the oblique is non-optional: without it every open half-plane bucket holds one element and
/// transitivity never fires). `angular_order` returns the CCW order reading no coordinate; we
/// check it against the CCW order computed *with* coordinates (atan2), which is the oracle.
#[test]
fn angular_order_around_a_vertex_is_coordinate_free_and_ccw() {
    // A pentagon prism giving y-, x-, and diagonal-normal side faces plus z caps.
    let profile = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]), // (0,0)-(3,0): y=0
        Point2::from_array([3.0, 0.0]), // (3,0)-(3,3): x=3
        Point2::from_array([3.0, 3.0]), // (3,3)-(2,3): y=3
        Point2::from_array([2.0, 3.0]), // (2,3)-(0,1): diagonal y=x+1
        Point2::from_array([0.0, 1.0]), // (0,1)-(0,0): x=0
    ])
    .unwrap();
    let m = replay(&[Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile,
        dist: 1.0,
    }])
    .unwrap();
    let a = m.live_solids()[0];
    let faces_tab = collect_planes(&m, a).unwrap();
    // One prism: no two faces are coplanar, so `plane_ix` is the identity and a face index and
    // its plane id coincide. Built through the real path anyway, so the test cannot drift.
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, _plane_ix, _cyls) = crate::planes::dense_planes(&faces_tab, &canon);
    assert_eq!(planes.len(), faces_tab.len(), "no coplanar pair in a prism");

    // Find a face by its outward normal direction (z cap, y-wall, x-wall, diagonal wall).
    let axis = |i: usize| {
        let n = planes[i].plane.normal();
        [n[0], n[1], n[2]]
    };
    let find = |f: &dyn Fn([f64; 3]) -> bool| (0..planes.len()).find(|&i| f(axis(i)));
    let w = find(&|n| n[0].abs() < 1e-9 && n[1].abs() < 1e-9).expect("z cap"); // z-normal
    let fpy = find(&|n| n[0].abs() < 1e-9 && n[2].abs() < 1e-9).expect("y wall"); // y-normal
    let fpx = find(&|n| n[1].abs() < 1e-9 && n[2].abs() < 1e-9).expect("x wall"); // x-normal
    let fpd =
        find(&|n| n[2].abs() < 1e-9 && n[0].abs() > 1e-6 && (n[0].abs() - n[1].abs()).abs() < 1e-6)
            .expect("diagonal wall");

    // The direction an edge (fp, s) actually runs, WITH coordinates — the oracle only.
    let n_w = planes[w].plane.normal();
    let dir = |fp: usize, s: i8| {
        let d = n_w.cross(planes[fp].plane.normal());
        [d[0] * s as f64, d[1] * s as f64, d[2] * s as f64]
    };
    // Five edges: both directions on the y-wall and x-wall lines, one on the diagonal.
    //
    // ★ **The pairs stay, and the oracle below reads *them*, not the `EdgeDir`s.** A direction's
    // representation is private to `combinatorics` on purpose (see `EdgeDir`); a test that needs
    // coordinates is an oracle, and an oracle reads its own inputs — reading them back out of
    // the value under test would derive the oracle from the answer.
    let raw = [(fpy, 1i8), (fpy, -1), (fpx, 1), (fpx, -1), (fpd, 1)];
    let edges: Vec<combinatorics::EdgeDir> = raw
        .iter()
        .map(|&(c, s)| combinatorics::EdgeDir::new(c, s))
        .collect();

    let order = angular_order(&crate::planes::test_judge(&planes), w, &edges)
        .expect("five distinct directions leave one vertex — nothing to refuse here");
    assert_eq!(
        order.len(),
        edges.len(),
        "every edge placed exactly once: {order:?}"
    );
    assert_eq!(order[0], 0, "order starts at the reference edge");

    // Oracle: the angular order of the actual directions (atan2), starting from edge 0.
    let ang = |e: (usize, i8)| {
        let d = dir(e.0, e.1);
        d[1].atan2(d[0])
    };
    let a0 = ang(raw[0]);
    let mut want: Vec<usize> = (0..edges.len()).collect();
    want.sort_by(|&i, &j| {
        let (ci, cj) = (
            (ang(raw[i]) - a0).rem_euclid(std::f64::consts::TAU),
            (ang(raw[j]) - a0).rem_euclid(std::f64::consts::TAU),
        );
        ci.partial_cmp(&cj).unwrap()
    });
    // The order is a consistent cyclic order — CW or CCW depending on orient_sign(w)'s
    // convention; both are correct rotations. So it matches the oracle, or the oracle with its
    // tail reversed (the same cycle traversed the other way, reference fixed).
    let mut rev_tail = order.clone();
    rev_tail[1..].reverse();
    assert!(
        order == want || rev_tail == want,
        "coordinate-free order is the atan2 cycle (either direction): got {order:?}, want {want:?}"
    );

    // ★ Transitivity fired: the oblique lands strictly between the +x and +y axis directions in
    // the cyclic order — a bucket-sort no-op could not place it. Checked as a cyclic adjacency
    // so it holds regardless of traversal direction.
    let cyc = if order == want { &order } else { &rev_tail };
    let pos = |e: usize| cyc.iter().position(|&x| x == e).unwrap();
    // edges: 0=(y,+) 1=(y,-) 2=(x,+) 3=(x,-) 4=(diag,+). The oblique (4) is 45°, between the
    // two axis directions flanking it. Identify its neighbours are axis edges, not each other.
    let obl = pos(4);
    let lo = cyc[(obl + cyc.len() - 1) % cyc.len()];
    let hi = cyc[(obl + 1) % cyc.len()];
    assert!(
        [0, 1, 2, 3].contains(&lo) && [0, 1, 2, 3].contains(&hi),
        "oblique is flanked by axis edges (transitivity placed it): {cyc:?}"
    );
}

/// On a cap plane, the cap rim and each side wall's top edge are one geometric edge traced
/// twice (seated + transversal). Merge collapses coincident edges by (wall, endpoint set) while
/// preserving contributions. Two footprint cases pin the count by hand.
#[test]
fn coincident_cap_edges_merge_by_footprint() {
    // (a) Same footprint: stacked cubes share the z=1 plane AND the same [0,1]² rim. Each of
    // the 4 rim edges is produced 4× (a-seated, a-wall, b-seated, b-wall) → 4 merged edges.
    {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        assert_eq!(merged.len(), 4, "one merged edge per rim edge: {merged:#?}");
        for e in &merged {
            assert_eq!(
                e.merged.len(),
                4,
                "a·b × seated·transversal: {:?}",
                e.merged
            );
        }
    }
    // (b) Different footprint: a cross. a's cap [0,3]×[1,2] and b's cap [1,2]×[0,3] are
    // different rims, so no a↔b coincidence: a's 4 pairs → 4, b's 4 pairs → 4 = 8.
    {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        assert_eq!(
            merged.len(),
            8,
            "4 per box, no a↔b coincidence: {merged:#?}"
        );
        for e in &merged {
            assert_eq!(
                e.merged.len(),
                2,
                "one solid × seated·transversal: {:?}",
                e.merged
            );
        }
    }
}

/// Split at crossings: a cross fixture makes four strict interior crossings; the chords cut
/// into the right pieces, ending at the right triples, and the crossing vertex is a clean
/// degree-4 the DCEL can order.
#[test]
fn split_cuts_the_cross_into_the_right_pieces() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    assert_eq!(merged.len(), 8, "8 merged edges before split");

    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
    // 4 chords crossed twice → 3 pieces each = 12; 4 outer walls uncrossed = 4; total 16.
    assert_eq!(split.len(), 16, "16 sub-segments: {}", split.len());

    // a's y=1 chord splits into exactly 3, ending at x=0,1,2,3 — checked by coordinate.
    let seg_pts = |s: &MergedSeg| {
        let mut u = pt(s.end[0], &jd);
        let mut v = pt(s.end[1], &jd);
        if u[0] > v[0] {
            std::mem::swap(&mut u, &mut v);
        }
        (
            [(u[0] * 1e6).round() as i64, (u[1] * 1e6).round() as i64],
            [(v[0] * 1e6).round() as i64, (v[1] * 1e6).round() as i64],
        )
    };
    // a's y=1 chord pieces: horizontal (both y=1), z=1, spanning consecutive x's.
    let y1: Vec<_> = split
        .iter()
        .filter(|s| {
            let (u, v) = seg_pts(s);
            u[1] == 1_000_000 && v[1] == 1_000_000 // both endpoints at y=1
                    && s.merged.iter().any(|(sd, _)| *sd == SolidSide::A)
        })
        .map(&seg_pts)
        .collect();
    assert_eq!(y1.len(), 3, "a's y=1 chord → 3 pieces: {y1:?}");
    let mut xs: Vec<[i64; 2]> = y1.iter().map(|(u, v)| [u[0], v[0]]).collect();
    xs.sort();
    assert_eq!(
        xs,
        vec![
            [0, 1_000_000],
            [1_000_000, 2_000_000],
            [2_000_000, 3_000_000]
        ],
        "pieces are [0,1][1,2][2,3], adjacent and endpoint-exact"
    );

    // The crossing vertex (1,1,1) is a clean degree-4 the DCEL can order: 4 incident pieces,
    // and their (wall, far-R) pairs — i.e. (fp, direction) — are all distinct (no coincidence
    // that angular_order's zero bucket would collapse). This is what the merge earned.
    let is_v = |n: NodeId| {
        let p = pt(n, &jd);
        (p[0] - 1.0).abs() < 1e-9 && (p[1] - 1.0).abs() < 1e-9 && (p[2] - 1.0).abs() < 1e-9
    };
    let incident: Vec<&MergedSeg> = split
        .iter()
        .filter(|s| is_v(s.end[0]) || is_v(s.end[1]))
        .collect();
    assert_eq!(incident.len(), 4, "degree-4 at (1,1,1): {}", incident.len());
    let mut keys = std::collections::HashSet::new();
    for s in &incident {
        // The far end as a handle on this edge's own line — carried by the segment now.
        let far_h = if is_v(s.end[0]) {
            s.end_h[1]
        } else {
            s.end_h[0]
        };
        keys.insert((s.wall, far_h));
    }
    assert_eq!(keys.len(), 4, "no coincident (fp,direction) — DCEL-ready");
}

/// Same-wall partial overlap is **resolved** by the per-wall overlay: a=[0,2], b=[1,3]
/// share y=1, their chords overlap on x∈[1,2]. The overlay splits the y=1 wall into three
/// non-overlapping pieces `[0,1] [1,2] [2,3]`, and the shared middle `[1,2]` carries a
/// contribution from **both** solids (which the label brick then reads per solid), while the
/// flanks are single-solid.
#[test]
fn partial_overlap_is_resolved() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();

    // The shared y=1 wall (a face at y=1).
    let y1 = planes
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("a y=1 face");
    // x-extent of a y=1 sub-segment, plus which solids contribute.
    let piece = |s: &MergedSeg| -> ([i64; 2], bool, bool) {
        let mut u = pt(s.end[0], &jd)[0];
        let mut v = pt(s.end[1], &jd)[0];
        if u > v {
            std::mem::swap(&mut u, &mut v);
        }
        let has = |sd: SolidSide| s.merged.iter().any(|(x, _)| *x == sd);
        (
            [(u * 1e6).round() as i64, (v * 1e6).round() as i64],
            has(SolidSide::A),
            has(SolidSide::B),
        )
    };
    let mut pieces: Vec<([i64; 2], bool, bool)> =
        split.iter().filter(|s| s.wall == y1).map(piece).collect();
    pieces.sort();
    assert_eq!(
        pieces,
        vec![
            ([0, 1_000_000], true, false),         // [0,1] A only
            ([1_000_000, 2_000_000], true, true),  // [1,2] BOTH solids
            ([2_000_000, 3_000_000], false, true), // [2,3] B only
        ],
        "y=1 overlap resolved into three pieces, middle carries both solids: {pieces:?}"
    );
}

/// ★★★★★ **The circle-against-ruling net is a check, and this is it being walked.**
///
/// A class carrying a circle (⊥ one cylinder) and rulings (∥ another) is fine until the two
/// **meet**, and then the crossing is `plane ∩ cylinder ∩ cylinder` — a point no road mints.
/// The
/// population reaches `ClassEdges::of` through a boolean, but the check is the first thing
/// that function does, so it can be walked here directly — which is the point: it must not be
/// a device behind a wall.
#[test]
fn a_circle_meeting_a_ruling_is_refused_by_name() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
    let q = |n: i128, d: i128| nacre_exact::Rat::new(n, d).unwrap();
    let z = |v: i128| nacre_exact::Rat::from_int(v);
    let def = |o: [nacre_exact::Rat; 3],
               dir: [nacre_exact::Rat; 3],
               e: [nacre_exact::Rat; 3],
               r: nacre_exact::Rat| {
        // `r` is the radius; the truth takes its square.
        nacre_topo::CylinderDef::new(o, dir, e, nacre_exact::BigRat::square_of(r)).unwrap()
    };
    // ⊥ the shared cap: a circle of radius 6/5 about `(3/2, 3/2)` on that plane.
    //
    // ★ The contribution is stated **whole** rather than left empty: the check reads
    // the arc a contribution covers, and an empty list is its fallback path. Carrying a real
    // contribution keeps this test on the road a model takes.
    let whole = vec![(SolidSide::A, SegKind::Seated { body_above: true }, None)];
    let circle = MergedCircle {
        cyl: 0,
        def: def(
            [q(3, 2), q(3, 2), z(0)],
            [z(0), z(0), z(1)],
            [z(1), z(0), z(0)],
            q(6, 5),
        ),
        merged: whole.clone(),
    };
    // ∥ it, and near enough (`1/2 < 1`) that both caps carry its rulings — at
    // `y = 3/2 ± √(3)/2`, which the circle above reaches past on either side.
    let ruling = MergedRuling {
        cyl: 1,
        def: def(
            [z(0), q(3, 2), q(1, 2)],
            [z(1), z(0), z(0)],
            [z(0), z(0), z(1)],
            z(1),
        ),
        side: 1,
        end: [
            crate::combinatorics::NodeId::three_planes(crate::combinatorics::Canon3::three([
                0, 1, 2,
            ])),
            crate::combinatorics::NodeId::three_planes(crate::combinatorics::Canon3::three([
                0, 1, 3,
            ])),
        ],
        merged: Vec::new(),
        orient: 1,
    };
    let out = ClassEdges::of(
        &jd,
        NO_CYLS,
        wc,
        &split,
        std::slice::from_ref(&circle),
        std::slice::from_ref(&ruling),
        &Aliases::default(),
    );
    assert!(
        matches!(
            out,
            Err(BoolError::Rejected {
                reason: RejectReason::CircleCrossesRuling,
                ..
            })
        ),
        "a circle across a ruling must be refused by name"
    );
    // ★ And a circle that stays clear of both rulings is not refused — without this the check
    // could be a constant and still be green. It sits **outside** the strip rather than inside
    // it, so it adds nothing to `MIXED_CLASS_AUDIT.inside`, which another test reads from the
    // same process-wide counters.
    let clear = MergedCircle {
        def: def(
            [q(3, 2), z(10), z(0)],
            [z(0), z(0), z(1)],
            [z(1), z(0), z(0)],
            q(1, 10),
        ),
        merged: whole,
        ..circle
    };
    // ★ Not «builds» — these edges are hand-made and the rest of the walk has no class table
    // for them — but **not refused by this name**, which is what says the check reads its
    // inputs rather than the mere presence of a circle beside a ruling.
    assert!(
        !matches!(
            ClassEdges::of(
                &jd,
                NO_CYLS,
                wc,
                &split,
                std::slice::from_ref(&clear),
                std::slice::from_ref(&ruling),
                &Aliases::default(),
            ),
            Err(BoolError::Rejected {
                reason: RejectReason::CircleCrossesRuling,
                ..
            })
        ),
        "a circle well inside the rulings shares no point with them"
    );
}

/// DCEL face-walk on the cross: 6 cells (1 outer + 5 bounded), exactly one winding -1, and the
/// center square is the unique cell with 2 A-edges + 2 B-edges. All checkable without labels.
#[test]
fn the_cross_arrangement_has_six_cells() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();

    let edges = ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();

    // ★★★ **An independent reading of every winding: the shoelace sign.**
    //
    // `loop_winding` reads the turn at one extreme node, which is correct only if the
    // comparison it picks that node with is an order. When it was not, exactly one ring of
    // five came back with the wrong sign and stayed *confident* about it — the engine only
    // noticed two layers later, as "no outer contour", and named the symptom. Realizing the
    // nodes and summing cross products says the same thing directly, so a disagreement points
    // at the ring rather than at the count.
    {
        let pl = &jd.planes[wc].plane;
        let n = pl.normal();
        let ax = if n[0].abs() < 0.9 {
            Vector3::from_array([1.0, 0.0, 0.0])
        } else {
            Vector3::from_array([0.0, 1.0, 0.0])
        };
        let u = n.cross(ax).normalize().expect("in-plane axis");
        let v = n.cross(u);
        let o = Point3::from_array([0.0; 3]);
        for c in &cells {
            let pts: Vec<[f64; 2]> = c
                .half_edges
                .iter()
                .filter_map(|&h| {
                    let t = three_plane_name(split[h / 2].end[h % 2]).expect("a three-plane node");
                    three_planes(
                        &jd.planes[t[0]].plane,
                        &jd.planes[t[1]].plane,
                        &jd.planes[t[2]].plane,
                    )
                    .map(|q| [u.dot(q - o), v.dot(q - o)])
                })
                .collect();
            assert_eq!(pts.len(), c.half_edges.len(), "every node realizes");
            let area: f64 = (0..pts.len())
                .map(|i| {
                    let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                    a[0] * b[1] - b[0] * a[1]
                })
                .sum();
            assert!(
                area.abs() > 1e-9,
                "a degenerate ring gives the oracle nothing to say"
            );
            assert_eq!(
                if area > 0.0 { 1i8 } else { -1 },
                c.winding,
                "shoelace {area:+.6e} disagrees with the reported winding {}",
                c.winding
            );
        }
    }
    assert_eq!(cells.len(), 6, "1 outer + 5 bounded: {}", cells.len());
    assert_eq!(
        cells.iter().filter(|c| c.winding == -1).count(),
        1,
        "exactly one outer (winding -1)"
    );
    assert_eq!(
        cells.iter().filter(|c| c.winding == 1).count(),
        5,
        "five bounded (winding +1)"
    );
    // Every half-edge belongs to exactly one cell.
    assert_eq!(face_of.len(), 2 * split.len(), "2E half-edges all placed");

    // The center square [1,2]²: the unique bounded cell with 2 A-edges and 2 B-edges. (a and b
    // ride disjoint walls, so there are no both-solid edges — provenance, not coincidence.)
    let solids_of =
        |he: usize| -> Vec<SolidSide> { split[he / 2].merged.iter().map(|(sd, _)| *sd).collect() };
    let mut center = 0;
    for c in cells.iter().filter(|c| c.winding == 1) {
        let (mut na, mut nb) = (0, 0);
        for &he in &c.half_edges {
            let s = solids_of(he);
            if s.contains(&SolidSide::A) {
                na += 1;
            }
            if s.contains(&SolidSide::B) {
                nb += 1;
            }
        }
        if na == 2 && nb == 2 {
            center += 1;
        }
    }
    assert_eq!(
        center, 1,
        "exactly one center square (2 A-edges + 2 B-edges)"
    );
}

/// Cell labels propagate from the void, and match footprint containment (an independent
/// oracle). Cross: W=z=1, a=[0,3]×[1,2], b=[1,2]×[0,3]. Both bodies below z=1 → aboves all F;
/// each cell's X_below = "cell centroid inside X's footprint".
#[test]
fn cell_labels_propagate_and_match_footprints() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
    let edges = ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
    let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

    for (i, c) in cells.iter().enumerate() {
        if c.winding == -1 {
            assert_eq!(labels[i], [false; 4], "unbounded is void");
            continue;
        }
        // Centroid of the cell (convex here) via the average of its vertex points.
        let verts: Vec<[f64; 3]> = {
            let mut vs: Vec<NodeId> = c
                .half_edges
                .iter()
                .map(|&h| split[h / 2].end[h % 2])
                .collect();
            vs.dedup();
            vs.iter().map(|&t| pt(t, &jd)).collect()
        };
        let cx = verts.iter().map(|p| p[0]).sum::<f64>() / verts.len() as f64;
        let cy = verts.iter().map(|p| p[1]).sum::<f64>() / verts.len() as f64;
        let in_a = (0.0..=3.0).contains(&cx) && (1.0..=2.0).contains(&cy);
        let in_b = (1.0..=2.0).contains(&cx) && (0.0..=3.0).contains(&cy);
        // aboves F (bodies below the cap); belows = footprint containment.
        assert_eq!(
            labels[i],
            [false, in_a, false, in_b],
            "cell centroid ({cx},{cy}) label vs footprint"
        );
    }
}

/// A straddle: a=[0,1]²×[0,2] passes THROUGH z=1 (transversal, no face there), b=[0,1]²×[1,2]
/// caps at z=1 from above. Same [0,1]² footprint → one square, 2 cells. The inner cell exercises
/// transversal-flips-both (a's above reaches T) and seated-wins on b in the same edge.
#[test]
fn a_straddle_drives_the_above_bit_and_seated_wins() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    // a passes through z=1 (no seated face); find the class b caps at z=1.
    let wc = (0..planes.len())
        .find(|&c| {
            solid_shell_handles(&m, b).into_iter().any(|sh| {
                m.shell(sh).faces.iter().any(|fh| {
                    plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, &surf_ix, &faces_tab)
                })
            })
        })
        .expect("b's z=1 cap class");
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
    let edges = ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
    let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

    assert_eq!(cells.len(), 2, "one square: inner + unbounded");
    let inner = cells.iter().position(|c| c.winding == 1).unwrap();
    // "above/below" is relative to +n_w (the class plane's stored normal), whose sign is a
    // convention the boolean's keep-rule treats symmetrically. Here n_w points -z, so +n_w
    // ("above") is physically below. A straddles z=1 → material on BOTH sides (T,T); b caps
    // from z>1 → material only on -n_w ("below") → (B_above=F, B_below=T). The point of the
    // fixture stands: the transversal drives A's above-bit to T (impossible on a pure cap),
    // and B is seated-wins (a single bit) on the very same rim edge.
    assert_eq!(
        labels[inner],
        [true, true, false, true],
        "A straddles both sides; B is one-sided (seated-wins); transversal reached the above-bit"
    );
}

/// Boolean keep-decision + result-face emission on the cross. Fuse keeps 5 cells (the plus
/// cap), Cut a−b keeps the 2 a-arms, Common keeps the 1 center — hand-verified — and the
/// emitted loops share interior edges in opposite directions with a consistent winding and a
/// flip that points the normal out of the kept solid.
#[test]
fn boolean_keep_and_result_faces_on_the_cross() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
    let edges = ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
    let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

    // Centroid of a face's ring (convex cells here).
    let centroid = |f: &LocalFace| -> [f64; 2] {
        let ps: Vec<[f64; 3]> = f.outer.expect_ring().iter().map(|&n| pt(n, &jd)).collect();
        [
            ps.iter().map(|p| p[0]).sum::<f64>() / ps.len() as f64,
            ps.iter().map(|p| p[1]).sum::<f64>() / ps.len() as f64,
        ]
    };
    let near = |c: [f64; 2], x: f64, y: f64| (c[0] - x).abs() < 1e-9 && (c[1] - y).abs() < 1e-9;

    // Fuse: all 5 bounded cells (the plus cap).
    let (fuse, _, _) = emit_faces(
        BoolKind::Fuse,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(fuse.len(), 5, "Fuse keeps the whole plus cap");
    // Cut a−b: exactly the two a-arms.
    let (cut, _, _) = emit_faces(
        BoolKind::Cut,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(cut.len(), 2, "Cut keeps the a-arms");
    let cut_c: Vec<[f64; 2]> = cut.iter().map(centroid).collect();
    assert!(
        cut_c.iter().any(|c| near(*c, 0.5, 1.5)) && cut_c.iter().any(|c| near(*c, 2.5, 1.5)),
        "the two survivors are the a-arms at (0.5,1.5),(2.5,1.5): {cut_c:?}"
    );
    // Common: exactly the center square.
    let (common, _, _) = emit_faces(
        BoolKind::Common,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(common.len(), 1, "Common keeps the center");
    assert!(near(centroid(&common[0]), 1.5, 1.5), "center at (1.5,1.5)");

    // Each emitted loop winds +1 (CCW about n_out(wc)) and has ≥3 distinct nodes.
    for f in fuse.iter().chain(&cut).chain(&common) {
        let ring: Vec<[usize; 3]> = f
            .outer
            .expect_ring()
            .iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect();
        assert!(ring.len() >= 3);
        assert_eq!(
            combinatorics::loop_winding(
                &jd,
                NO_CYLS,
                wc,
                &combinatorics::ring_from_names(wc, &ring).unwrap()
            )
            .unwrap(),
            1,
            "CCW about n_out"
        );
    }

    // flip oracle (coordinate, test-only): the result normal points away from the kept
    // chamber. Fuse keeps below (bodies below the cap), so n_result·n_w > 0.
    let n_w = planes[wc].plane.normal();
    let os = planes[wc].frame_sign as f64;
    for f in &fuse {
        let n_result = os * if f.flip { -1.0 } else { 1.0 };
        let dot = n_result * n_w.dot(n_w); // n_result·n_w, |n_w|²>0
        assert!(
            dot > 0.0,
            "Fuse (keep below) normal points +n_w side: flip={}",
            f.flip
        );
    }

    // Edge-parity: count how many times each undirected edge is emitted, and the net
    // direction. An edge shared by two survivors (count 2) must net to 0 (opposite directions
    // — assemble's "each edge twice"); a plus-cap boundary edge (count 1) is excluded (it
    // closes against a wall face only at assembly). At least one interior edge must exist.
    let mut edges: HashMap<([usize; 3], [usize; 3]), (i32, i32)> = HashMap::new();
    for f in &fuse {
        let ns: Vec<[usize; 3]> = f
            .outer
            .expect_ring()
            .iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect();
        for w in ns
            .windows(2)
            .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
        {
            let (key, sign) = if w[0] < w[1] {
                ((w[0], w[1]), 1)
            } else {
                ((w[1], w[0]), -1)
            };
            let e = edges.entry(key).or_insert((0, 0));
            e.0 += 1;
            e.1 += sign;
        }
    }
    let mut interior = 0;
    for (e, (count, net)) in &edges {
        if *count == 2 {
            interior += 1;
            assert_eq!(
                *net, 0,
                "interior edge {e:?} emitted in opposite directions"
            );
        }
    }
    assert!(
        interior > 0,
        "the plus cap has interior edges between arms and center"
    );
}

/// Intermediate lock (before assembly, A/B isolation): the driver emits exactly 10 faces for
/// stacked Fuse (z=0 cap + z=2 cap + 4 walls × 2 z-split), no z=1 face, and every undirected
/// edge appears exactly twice across all faces — so a later NON_MANIFOLD reject is a weld gap,
/// not an emit gap.
#[test]
fn stacked_fuse_emits_ten_closed_faces() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        plane_ix,
        class_owner,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        Default::default(),
    );
    let (faces, _, _) = trace_result_faces(
        &m,
        BoolKind::Fuse,
        a,
        b,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .unwrap();
    assert_eq!(faces.len(), 10, "z=0 + z=2 + 4 walls×2");

    // Every undirected edge (sorted triple pair) is used exactly twice — a closed shell.
    let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
    for f in &faces {
        let ns: Vec<[usize; 3]> = f
            .outer
            .expect_ring()
            .iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect();
        for w in ns
            .windows(2)
            .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
        {
            let key = if w[0] < w[1] {
                (w[0], w[1])
            } else {
                (w[1], w[0])
            };
            *count.entry(key).or_insert(0) += 1;
        }
    }
    assert!(
        count.values().all(|&c| c == 2),
        "every edge used exactly twice (closed shell): {:?}",
        count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
    );
}

/// First end-to-end: drive every plane class, assemble, check volume + manifold. Stacked cubes
/// Fuse = a 1×1×2 box, volume 2.0; the z=1 interface face vanishes (kept both sides).
#[test]
fn end_to_end_stacked_fuse_is_a_tall_box() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
    assert_eq!(solids.len(), 1, "one connected solid");
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "manifold: {vs:?}");
    let vol = nacre_props::mass_props(&m, solids[0]).unwrap().volume;
    assert!(
        (vol - 2.0).abs() < 1e-9,
        "stacked Fuse volume 2.0, got {vol}"
    );
}

/// Non-degenerate all-three: overlapping cubes a=[0,1]³, b=[0.5,1.5]³.
/// Fuse = 2 − 0.5³ = 1.875, Common = 0.5³ = 0.125, Cut = 1 − 0.125 = 0.875.
#[test]
fn end_to_end_overlapping_cubes_all_three() {
    let vol_of = |kind: BoolKind| -> f64 {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        m.rebuild_adjacency();
        let solids = boolean(&mut m, kind, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum()
    };
    assert!((vol_of(BoolKind::Fuse) - 1.875).abs() < 1e-9, "Fuse 1.875");
    assert!((vol_of(BoolKind::Cut) - 0.875).abs() < 1e-9, "Cut 0.875");
    assert!(
        (vol_of(BoolKind::Common) - 0.125).abs() < 1e-9,
        "Common 0.125"
    );
}

/// cube[0,3]³ and a square prism [1,2]²×[−1,4] piercing it in z; run the boolean, validate,
/// and sum (volume, area) over the result solids.
fn tunnel_vol_area(kind: BoolKind) -> (f64, f64) {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    m.rebuild_adjacency();
    let solids = boolean(&mut m, kind, a, b).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
    solids.iter().fold((0.0, 0.0), |acc, &s| {
        let p = nacre_props::mass_props(&m, s).unwrap();
        (acc.0 + p.volume, acc.1 + p.area)
    })
}

/// The first result with a hole. On the z=0/z=3 caps the arrangement is two disjoint loops
/// (cube rim `O`, tunnel mouth `I`); the mouth is a hole of the annular cap, which cannot fall
/// out as an absent cell (tiling the annulus needs a rim node the adjacent wall lacks) — it
/// must emit as an inner ring. Cut volume 24 (27−3); AREA 64 (not 54) proves a bore, not a
/// filled cap — volume alone cannot tell a through-hole from a blind dent.
#[test]
fn tunnel_cut_is_bored_not_dented() {
    let (v, area) = tunnel_vol_area(BoolKind::Cut);
    assert!((v - 24.0).abs() < 1e-9, "Cut vol 24 (27−1·1·3), got {v}");
    assert!(
        (area - 64.0).abs() < 1e-9,
        "Cut area 64 (36 walls + 8+8 punched caps + 12 tunnel walls; a filled mouth is 54): {area}"
    );
}

/// Fuse: the bar protrudes z∈[−1,0] and z∈[3,4], so each cap is still annular (the peg passes
/// through the mouth, kept on both sides) — volume 29 (27+5−3), area 62.
#[test]
fn tunnel_fuse_has_annular_caps() {
    let (v, area) = tunnel_vol_area(BoolKind::Fuse);
    assert!((v - 29.0).abs() < 1e-9, "Fuse vol 29 (27+5−3), got {v}");
    assert!(
        (area - 62.0).abs() < 1e-9,
        "Fuse area 62 (cube 52 + two 1×1×1 pegs at 5 each): {area}"
    );
}

/// Common: `nest_cells` still finds the same annulus grouping, but `keep` leaves only the mouth
/// (inside `I`); the annular host is skipped, its inner ring never leaks, and the mouth emits a
/// plain `[1,2]²` cap. Net a plain 1×1×3 box — vol 3, area 14. A stronger control than a
/// hole-free fixture: it exercises the nesting path and proves it does not leak a spurious hole.
#[test]
fn tunnel_common_is_the_bar_box_control() {
    let (v, area) = tunnel_vol_area(BoolKind::Common);
    assert!(
        (v - 3.0).abs() < 1e-9,
        "Common vol 3 ([1,2]²×[0,3]), got {v}"
    );
    assert!(
        (area - 14.0).abs() < 1e-9,
        "Common area 14 (2·1 caps + 4·3 walls): {area}"
    );
}

/// Intermediate lock (before assembly, A/B isolation): the Cut tunnel emits exactly 10 faces
/// (2 annular caps + 4 cube walls + 4 tunnel walls), exactly 2 carry a non-empty inner ring
/// (the caps), and every undirected edge across all rings (outer + inner) is used exactly twice
/// — a closed shell, so a later NON_MANIFOLD is a weld gap, not an emit gap.
#[test]
fn tunnel_cut_emits_ten_faces_two_annular() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        plane_ix,
        class_owner,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        Default::default(),
    );
    let (faces, _, _) = trace_result_faces(
        &m,
        BoolKind::Cut,
        a,
        b,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .unwrap();
    assert_eq!(
        faces.len(),
        10,
        "2 annular caps + 4 cube walls + 4 tunnel walls"
    );
    assert_eq!(
        faces.iter().filter(|f| !f.inner.is_empty()).count(),
        2,
        "exactly the two annular caps carry a hole"
    );

    // Every undirected edge across outer + inner rings is used exactly twice (closed shell).
    let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
        ns.iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect()
    };
    let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
    for f in &faces {
        for ring in f.poly_rings() {
            let ns = triples(ring);
            for w in ns
                .windows(2)
                .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
            {
                let key = if w[0] < w[1] {
                    (w[0], w[1])
                } else {
                    (w[1], w[0])
                };
                *count.entry(key).or_insert(0) += 1;
            }
        }
    }
    assert!(
        count.values().all(|&c| c == 2),
        "every edge used exactly twice (closed shell): {:?}",
        count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
    );
}

/// Multi-hole lock: a slab fused with a U (its two prongs pierce the slab) emits exactly one
/// face with **two** inner rings — the y=2 slab annulus, holed by both prongs — and every
/// undirected edge across all rings is used exactly twice (closed shell). This exercises
/// `nest_cells`' multi-hole support.
#[test]
fn u_slab_fuse_emits_a_two_hole_face() {
    let u_profile = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([3.0, 0.0]),
        Point2::from_array([3.0, 2.3]),
        Point2::from_array([2.0, 2.3]),
        Point2::from_array([2.0, 1.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([1.0, 2.0]),
        Point2::from_array([0.0, 2.0]),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        frame: SketchFrame::world(&Model::new(), Axis::Z),
        profile: u_profile,
        dist: 1.0,
    }])
    .unwrap();
    let u = m.live_solids()[0];
    let slab = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-0.5, 1.5, -0.5]),
        Point3::from_array([3.5, 2.5, 1.5]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        plane_ix,
        class_owner,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, u, slab).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(u, &inc_a), (slab, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        Default::default(),
    );
    let (faces, _, _) = trace_result_faces(
        &m,
        BoolKind::Fuse,
        u,
        slab,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .unwrap();

    // Exactly one face carries two inner rings: the y=2 slab annulus, holed by both prongs.
    assert_eq!(
        faces.iter().filter(|f| f.inner.len() == 2).count(),
        1,
        "the slab face on y=2 has two holes (the U's two prongs)"
    );

    // Every undirected edge across outer + inner rings is used exactly twice (closed shell).
    let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
        ns.iter()
            .map(|&n| three_plane_name(n).expect("a three-plane node"))
            .collect()
    };
    let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
    for f in &faces {
        for ring in f.poly_rings() {
            let ns = triples(ring);
            for w in ns
                .windows(2)
                .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
            {
                let key = if w[0] < w[1] {
                    (w[0], w[1])
                } else {
                    (w[1], w[0])
                };
                *count.entry(key).or_insert(0) += 1;
            }
        }
    }
    assert!(
        count.values().all(|&c| c == 2),
        "every edge used exactly twice (closed shell): {:?}",
        count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
    );
}
