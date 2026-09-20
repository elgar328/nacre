use super::*;
use crate::SketchFrame;
use nacre_scalar::Axis;

/// ★ **A fixture with no cylinders, said as a fact rather than left as a hole.** The table is
/// how a `NodeId::Pierce` reaches its definition, so an all-plane fixture has nothing to put
/// in it — and a bare `&[]` at a call site reads like something forgotten.
const NO_CYLS: &[crate::planes::WorkingCyl] = &[];

/// ★★★ **Where a named curved ring stops, and under which name.**
///
/// ★★★★★ **The vessel carries a curved ring through, and only a *collapsed* name is refused.**
///
/// This used to assert the opposite: a pierce corner was `RingFail::Pierce` and a curved
/// carrier `RingFail::CurvedWall`, because the ring came out as plane ids and could not hold
/// either. Both refusals were the type running out rather than a decision, and they are gone
/// with the projections. What still stops here is a name that is degenerate *as a name* — two
/// of its three planes equal — which is a fact about the triple and not about cylinders.
///
/// ★ Rewritten rather than patched: the old proposition ("the plane road cannot use such a
/// ring") became false, and a lock whose sentence is false is worse than no lock.
#[test]
fn a_named_curved_ring_rides_through_and_only_a_collapsed_name_stops() {
    let three =
        |a, b, c| combinatorics::NodeId::three_planes(combinatorics::Canon3::three([a, b, c]));
    let pierce = combinatorics::NodeId::pierce(0, 1, 0, nacre_topo::QuadRoot::Lo);
    let plane = |c| crate::boolean::Wall::Plane(c);
    let ruling = crate::boolean::Wall::Ruling {
        cyl: 0,
        side: 1,
        up: true,
    };
    // A curved carrier rides through, carried as itself.
    let curved_carrier = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
        walls: vec![plane(1), ruling, plane(3)],
        arc_ccw: vec![None; 3],
        concurrencies: vec![],
    };
    let (_, walls) = plane_ring(&curved_carrier).expect("a curved carrier is describable");
    assert_eq!(
        walls[1], ruling,
        "the carrier is the producer's, unflattened"
    );
    // So does a pierce corner, as its own name.
    let curved_corner = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), pierce, three(0, 1, 3)],
        walls: vec![plane(1), ruling, plane(3)],
        arc_ccw: vec![None; 3],
        concurrencies: vec![],
    };
    let (ts, _) = plane_ring(&curved_corner).expect("a pierce corner is describable");
    assert_eq!(ts[1], pierce, "the corner is the producer's, unflattened");
    // ★ What is still refused, and for a reason that has nothing to do with cylinders.
    let collapsed = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 0, 3), three(0, 1, 3)],
        walls: vec![plane(1), plane(2), plane(3)],
        arc_ccw: vec![None; 3],
        concurrencies: vec![],
    };
    assert!(matches!(plane_ring(&collapsed), Err(RingFail::Collapsed)));
    // And the plane-only ring still comes back with both halves.
    let plain = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
        walls: vec![plane(1), plane(2), plane(3)],
        arc_ccw: vec![None; 3],
        concurrencies: vec![],
    };
    let (ts, walls) = plane_ring(&plain).expect("a plane ring");
    assert_eq!(ts.len(), 3);
    assert_eq!(walls, vec![plane(1), plane(2), plane(3)]);
}

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

/// The engine entry with the evidence dropped — these tests assert geometry, and the report
/// has its own tests. Shadows [`super::boolean`] so the call sites read as they always did.
fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    super::boolean(model, kind, a, b).map(|(solids, ..)| solids)
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
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    m.add_cuboid(
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
                .tri
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
    // a graze — see `run_body_above`. A single chord spanning [0,2] would be the old,
    // occupancy-blind bridging.
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
    m.add_cuboid(
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
                .tri
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
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
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
        let a = m.add_cuboid(
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = m.add_cuboid(
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

/// Same-wall partial overlap (E5) is **resolved** by the per-wall overlay: a=[0,2], b=[1,3]
/// share y=1, their chords overlap on x∈[1,2]. The overlay splits the y=1 wall into three
/// non-overlapping pieces `[0,1] [1,2] [2,3]`, and the shared middle `[1,2]` carries a
/// contribution from **both** solids (which the label brick then reads per solid), while the
/// flanks are single-solid.
#[test]
fn partial_overlap_is_resolved() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
        .position(|p| p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12))
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let q = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
    let z = |v: i128| nacre_scalar::Rat::from_int(v);
    let def = |o: [nacre_scalar::Rat; 3],
               dir: [nacre_scalar::Rat; 3],
               e: [nacre_scalar::Rat; 3],
               r: nacre_scalar::Rat| {
        // `r` is the radius; the truth takes its square.
        nacre_topo::CylinderDef::new(o, dir, e, nacre_scalar::BigRat::square_of(r)).unwrap()
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
            crate::combinatorics::NodeId::ThreePlane([0, 1, 2]),
            crate::combinatorics::NodeId::ThreePlane([0, 1, 3]),
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 1.0, 0.0]),
        Point3::from_array([3.0, 2.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
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
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
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
/// `nest_cells`' multi-hole support (the dropped `HOLE_MULTI` reject).
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
    let slab = m.add_cuboid(
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

// --- Rotation-generality: the trace engine's decisions are coordinate-free. Rigidly rotating
// both operands by the same isometry must leave the result invariant. A non-90° angle flags the
// solid, whose faces record the motion, routing every predicate to the exact frame3 backend. `rot30` (lib.rs)
// is in a sibling test module and unreachable here, so the isometries are built inline.

/// 30° about `axis` through (1,1,0) — non-90°, so the motion is recorded (exact frame3 path).
fn rot_iso(axis: nacre_scalar::Axis) -> nacre_scalar::Isometry {
    use nacre_scalar::{Angle, Isometry, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// Rotate a solid by each axis in turn. A compound tilt needs `rebuild_adjacency` BETWEEN the
/// transforms (matching the oracle's `rotated_boolean_matches_occt`), or the second reads a
/// stale topology.
fn tilt(m: &mut Model, mut s: Handle<Solid>, axes: &[nacre_scalar::Axis]) -> Handle<Solid> {
    for &ax in axes {
        s = transform(m, s, &rot_iso(ax)).unwrap();
        m.rebuild_adjacency();
    }
    s
}

/// The boolean's error and the class audit's `failed_at` are the **same** reject.
///
/// They are two consumers of one `DeclineKind → RejectReason` mapping, and before
/// `decline_to_reject` they were two copies of it. A copy that drifts makes the audit — the
/// tool used to debug a reject — disagree with the reject being debugged, which is the worst
/// possible time to be lying.
///
/// The model is **chosen by measurement, not by taste**: a sweep over the fixture corpus found
/// no boolean that declines inside the arrangement at all (every fixture builds), so the
/// agreement had to be pinned on an input that still stops after the classes have run. If a
/// later capability makes it build, the fix is to take whatever still declines — not to weaken
/// the assertion.
#[test]
fn the_audit_does_not_invent_failures() {
    // The audit's scope is the per-class pipeline, and its duty is to run **the pipeline the
    // boolean runs** — with the alias fixpoint. Audited against an empty alias table it
    // reported `UnorderedEdges` for three classes of this input (names two classes discover
    // for each other were missing), failures the boolean never had: an instrument that
    // invents readings. The boolean's own reject here (`StraightAngle`) comes from the
    // assembly's vertex naming, after every class pipeline has run and outside the audit's
    // scope — so the audit's honest answer for this input is "no class failed".
    //
    // ★ **The fixture has moved twice, exactly as the note above prescribes**, and each move
    // is a capability the kernel gained. First it was the same pair at 30° with `Cut`,
    // rejecting `DegenerateWitness` — a reject from the component outwardness test, which the
    // nesting-parity label replaced. Then it was that pair at 60°, rejecting `StraightAngle`:
    // two unit cubes that only *touch*, which now come back as the two bodies they are
    // (measured: `Common` empty, fused volume 2.0, `validate` clean).
    //
    // So the fixture is now a pinch that **cannot** part: A and B meet only along the line
    // `x = 2, y = 2`, and a bridge overlapping both runs the material around the contact, so
    // cutting there leaves one piece. The closed-shell guard in `boolean::reconstruct` rejects
    // it — after every class pipeline has run, which is the property this test needs.
    //
    // (This test once asserted the audit reports the boolean's *class-level* reject, on a
    // fixture chosen as "some input that rejects" — a bar rotated through an L-shaped
    // target. The knife-edge fix turned that family, and every class-level-rejecting valid
    // input we could construct, into answers; the shared mapping the old test guarded,
    // `decline_to_reject`, is one function called by both consumers, so it cannot drift.)
    let build = || -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
            let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            s
        };
        let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
        let bridge = cub(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
        let b = cub(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
        let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
        m.rebuild_adjacency();
        (m, ab[0], b)
    };
    let (mut m, a, b) = build();
    let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::NonManifoldResultEdge,
                ..
            }
        ),
        "the fixture's premise: a reject from outside the class pipeline (got {err:?})"
    );
    let (m, a, b) = build();
    let audits = frame_audit(&m, BoolKind::Fuse, a, b).unwrap();
    let failed: Vec<RejectReason> = audits.iter().filter_map(|x| x.failed_at).collect();
    assert_eq!(
        failed,
        vec![],
        "every class runs clean under the boolean's own alias table — a failure invented \
             here is an artifact of running a different pipeline than the boolean runs"
    );
}

/// Rigidly rotating both operands (same single-Z tilt) leaves all three booleans' volumes
/// invariant — the axis values (hand-anchored by `end_to_end_overlapping_cubes_all_three`)
/// transfer to the rotated case. A coordinate-dependent decision that flipped under rotation
/// would add/drop a cell and move the volume by O(0.1), far past 1e-9.
#[test]
fn rotated_overlapping_cubes_all_three() {
    use nacre_scalar::Axis;
    let vol_of = |kind: BoolKind| -> f64 {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let a = tilt(&mut m, a, &[Axis::Z]);
        let b = tilt(&mut m, b, &[Axis::Z]);
        let solids = boolean(&mut m, kind, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum()
    };
    assert!(
        (vol_of(BoolKind::Fuse) - 1.875).abs() < 1e-9,
        "rotated Fuse invariant 1.875"
    );
    assert!(
        (vol_of(BoolKind::Cut) - 0.875).abs() < 1e-9,
        "rotated Cut invariant 0.875"
    );
    assert!(
        (vol_of(BoolKind::Common) - 0.125).abs() < 1e-9,
        "rotated Common invariant 0.125"
    );
}

/// The through-tunnel Cut is invariant under every orientation: single Z, X, Y, and compound
/// Z∘X. The axis-DEPENDENCE was the tell of the bug (Y worked; Z/X declined before the
/// `Judge::orient3d` on-plane fix), so all four orientations returning 24 is the fix's direct
/// regression lock. The cube's own z=0/z=3 caps exercise the seated path under rotation.
#[test]
fn rotated_tunnel_cut_all_orientations() {
    use nacre_scalar::Axis;
    let vol_of = |axes: &[Axis]| -> f64 {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let a = tilt(&mut m, a, axes);
        let b = tilt(&mut m, b, axes);
        let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "manifold: {vs:?}");
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum()
    };
    for axes in [
        &[Axis::Z][..],
        &[Axis::X][..],
        &[Axis::Y][..],
        &[Axis::Z, Axis::X][..],
    ] {
        let v = vol_of(axes);
        assert!(
            (v - 24.0).abs() < 1e-9,
            "tunnel Cut vol 24 for {axes:?}, got {v}"
        );
    }
}

/// Compound-tilted (no face normal axis-aligned) tunnel Cut: AREA 64 is invariant (the bore
/// discriminator volume cannot see), and the pre-assembly face set is combinatorially identical
/// to the axis case (`tunnel_cut_emits_ten_faces_two_annular`) — 10 faces, 2 with an inner ring,
/// every undirected edge twice — proving the arrangement itself survived rotation.
#[test]
fn rotated_tunnel_area_and_faces() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
    let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);

    // Combinatorial invariant (pre-assembly, A/B isolation).
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
    assert_eq!(faces.len(), 10, "rotated arrangement keeps 10 faces");
    assert_eq!(
        faces.iter().filter(|f| !f.inner.is_empty()).count(),
        2,
        "the two annular caps survive rotation"
    );
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
        "every edge used exactly twice: {:?}",
        count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
    );

    // Area (assembled).
    let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let area: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().area)
        .sum();
    assert!(
        (area - 64.0).abs() < 1e-9,
        "rotated Cut area 64 (bore survives rotation), got {area}"
    );
}

/// The sharp tripwire: a compound-tilted tunnel must decline nothing, exactly like the axis
/// baseline `axis_aligned_cubes_decline_nothing`. This is the EXACT site the bug broke
/// (`coincident-features` from a `wall == wc` degenerate), so it fails immediately if the
/// `Judge::orient3d` on-plane fix is reverted.
#[test]
fn rotated_tunnel_declines_nothing() {
    use nacre_scalar::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
    let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);
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
    for wc in 0..planes.len() {
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
        assert!(
            tr.declined.is_empty(),
            "rotated class {wc} declined: {tr:?}"
        );
    }
}

/// The z=1 class both solids seat a cap on.
fn shared_cap_class(
    m: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
) -> usize {
    let n_planes = plane_ix
        .iter()
        .map(|c| c.plane())
        .max()
        .map_or(0, |m| m + 1);
    (0..n_planes)
        .find(|&c| {
            let seats = |s: Handle<Solid>| {
                solid_shell_handles(m, s).into_iter().any(|sh| {
                    m.shell(sh).faces.iter().any(|fh| {
                        plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, surf_ix, faces)
                    })
                })
            };
            seats(a) && seats(b)
        })
        .expect("a shared z=1 cap class")
}

/// ★★ **The production path, run on a cylinder-bearing model.**
///
/// A `.plane()` call on the production road with cylinder rows flowing into it —
/// the `work` class list, the decline witness, the report's `class_of` — aborts the
/// kernel on the first drill, and reading the ~30 call sites does not find one. **Running
/// the road is what finds them**, which is why this test exists.
///
/// What it asserts is deliberately weak on geometry and strong on survival: the tracer
/// completes, every face it emits is still a plane face, and the
/// drilled cap carries its circular hole.
#[test]
fn the_production_road_survives_a_cylinder_past_the_stopper() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // a through-hole: circles on both box caps
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        class_owner,
        n_a,
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
        // `Proved` on purpose: the reuse guard must be what turns the shortcut off, not the
        // caller. Without it `pass_through` would carry planes across and drop the lateral.
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .expect("the drill population traces");
    assert!(
        faces.iter().all(|f| matches!(f.surf, ClassIx::Plane(_))),
        "the plane arrangement emits plane faces only; bands are C4b-2"
    );
    // The box: 4 walls + 2 caps, and each cap carries the drill's circular hole.
    assert_eq!(faces.len(), 6, "{faces:?}");
    let holed = faces
        .iter()
        .filter(|f| {
            f.inner
                .iter()
                .any(|b| matches!(b, crate::boolean::Bound::Circle { .. }))
        })
        .count();
    assert_eq!(holed, 2, "both caps are drilled: {faces:?}");
}

/// The rulings-road harness: the through-boss geometry
/// — a plate and a cylinder overlapping its full height, axis exactly
/// on the `x = 40` wall plane — assembled **past the standing gate** from production parts:
/// `plane_index_setup_inner` (the gate-free half) plus the cylinder table built the way the
/// gate builds it. Returns everything the two locks below read, and the record the
/// gate produces: `(wall class, cylinder)` marked not-proven-clear.
#[allow(clippy::type_complexity)]
fn armed_through_boss() -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::planes::PlaneSetup,
    usize,
    std::collections::HashSet<(usize, usize)>,
) {
    armed_through_boss_z(-10.0, 50.0)
}

/// [`armed_through_boss`] with the boss's axial extent chosen: `z_lo` and height `h`. The
/// default runs through the plate (`−10`, `50`); a boss whose lower cap sits **inside** the
/// plate (`10`, `30`) has a lateral whose lower boundary is a **chain** — arcs at z = 10
/// (outside the plate) and z = 20 (inside) joined by rulings — the staircase the corner
/// rule is watched on.
#[allow(clippy::type_complexity)]
fn armed_through_boss_z(
    z_lo: f64,
    h: f64,
) -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::planes::PlaneSetup,
    usize,
    std::collections::HashSet<(usize, usize)>,
) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 20.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array([40.0, 20.0, z_lo]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        5.0,
        h,
    );
    m.rebuild_adjacency();
    let (mut setup, cyl_surfs) = crate::planes::plane_index_setup_inner(&m, plate, boss).unwrap();
    for &surf in &cyl_surfs {
        let nacre_topo::Surface::Cylinder { def, .. } = m.surface(surf) else {
            unreachable!("a cylinder row carries a cylinder truth")
        };
        let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(surf) else {
            unreachable!("push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache")
        };
        setup.cyls.push(crate::planes::WorkingCyl {
            surf,
            def: def.clone(),
            realized: *cache,
            owner: crate::planes::SolidSide::A,
        });
    }
    let wc = setup
        .geom
        .iter()
        .position(|p| p.tri.iter().all(|q| (q.as_array()[0] - 40.0).abs() < 1e-12))
        .expect("the x = 40 wall class");
    let crossings: std::collections::HashSet<(usize, usize)> = [(wc, 0)].into_iter().collect();
    (m, plate, boss, setup, wc, crossings)
}

/// ★ **The gate's record arms the rulings road — and only the record.** With the
/// through-boss pair listed, the boss's trace on the wall class is the rectangle: two
/// rulings (one per side, Pierce-named ends, no plain segments) and two cap chords
/// (`Lo`/`Hi` roots on the cap classes). With the record empty — every production call
/// today — the same trace is empty: the negative control that pins "an empty record
/// changes nothing", which is the very population an unconditionally-firing arm broke
/// (14 arc-family tests red, measured).
#[test]
fn the_gates_record_arms_the_ruling_trace() {
    let (m, _plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut tr = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        crossings.clone(),
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(
        tr.segs.len(),
        2,
        "one chord per cap — a segment of the sweep since cell ⑩"
    );
    let mut sides: Vec<i8> = tr.rulings.iter().map(|r| r.side).collect();
    sides.sort_unstable();
    assert_eq!(sides, [-1, 1], "one ruling per side");
    for r in &tr.rulings {
        for nd in r.end {
            let (_, cyl, _) = combinatorics::pierce_name(nd).expect("a Pierce end");
            assert_eq!(cyl, 0);
        }
        assert!(matches!(r.kind, SegKind::Transversal { .. }));
    }
    for c in &tr.segs {
        assert_eq!(
            c.end_h,
            [combinatorics::EndPin::Cylinder; 2],
            "a chord's ends are the two roots"
        );
        // The chord rides the cap's own class — a ⊥ plane at z = −10 or z = 40.
        let z = setup.geom[c.wall].tri[0].as_array()[2];
        assert!(
            setup.geom[c.wall]
                .tri
                .iter()
                .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
                && ((z + 10.0).abs() < 1e-12 || (z - 40.0).abs() < 1e-12),
            "cap class at z = {z}"
        );
    }
    // The negative control: today's record.
    let mut tr0 = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        Default::default(),
        &mut tr0,
    );
    assert!(
        tr0.rulings.is_empty() && tr0.segs.is_empty(),
        "empty record, empty road"
    );
}

/// ★ **The armed arrangement digests the rectangle** — production bricks end to end on the
/// wall class. Each ruling is cut at its T-junctions with the plate's `z = 0` and `z = 20`
/// lines (three pieces each), those lines split at the same Pierce nodes, the chords close
/// the far ends, and the walk closes the subdivision: one unbounded contour and five
/// bounded cells — plate-left, plate-right, the overlap band, and the rectangle's two
/// overhangs.
#[test]
fn the_armed_wall_class_walks_to_closed_cells() {
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let tr = trace_on_class_of(
        &m,
        plate,
        boss,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_a,
        &setup.inc_b,
        &setup.plane_ix,
        crossings,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
    let split = drop_newsless(split).unwrap();
    let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
    // This fixture's class is parallel to the only cylinder in it, so it traces rulings and
    // no circle. Not a general law: a class may carry both, when the circle
    // belongs to a *different*, perpendicular cylinder.
    assert!(circles.is_empty(), "this class meets no ⊥ cylinder");
    let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
    assert_eq!(rulings.len(), 2);
    let edges = ClassEdges::of(
        &jd,
        &setup.cyls,
        wc,
        &split,
        &circles,
        &rulings,
        &Aliases::default(),
    )
    .unwrap();
    assert_eq!(
        edges.rulings.len(),
        6,
        "each ruling cut at z = 0 and z = 20"
    );
    let (cells, _face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 6, "five bounded cells and the outer");
    assert_eq!(
        cells.iter().filter(|c| c.winding == -1).count(),
        1,
        "one connected skeleton, one outer contour"
    );
}

/// One armed class's arrangement, through the production bricks (the chain
/// `the_armed_wall_class_walks_to_closed_cells` spells out) — for the locks that need
/// several classes' products at once.
fn armed_class_edges<'a>(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &'a crate::planes::PlaneSetup,
    jd: &Judge<'a, WorkingPlane>,
    wc: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
) -> ClassEdges<'static> {
    let tr = trace_on_class_of(
        m,
        plate,
        boss,
        wc,
        jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_a,
        &setup.inc_b,
        &setup.plane_ix,
        crossings.clone(),
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    let merged = merge_coincident(jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
    let split = drop_newsless(split).unwrap();
    let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
    let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
    let e = ClassEdges::of(
        jd,
        &setup.cyls,
        wc,
        &split,
        &circles,
        &rulings,
        &Aliases::default(),
    )
    .unwrap();
    // Owned copies so the borrows above may end — a test convenience, not a production shape.
    ClassEdges {
        segs: std::borrow::Cow::Owned(e.segs.into_owned()),
        arcs: std::borrow::Cow::Owned(e.arcs.into_owned()),
        rulings: std::borrow::Cow::Owned(e.rulings.into_owned()),
        circles: std::borrow::Cow::Owned(e.circles.into_owned()),
        cut_rims: e.cut_rims,
    }
}

/// **A lateral face answers the crossing question by the same parity, on its own chart**
/// — the digon oracle's lateral twin, on the through-boss Fuse's lateral. The
/// truth is known in world coordinates, so every crossing a rational ray makes with the
/// cylinder is checked: rays `{y = y₀, z = z₀}` over a half-step lattice (two roots each,
/// at x = 40 ± √(25 − (y₀ − 20)²) — irrational θ) and the **station column**
/// `{x = 40, z = z₀}`, whose roots are exactly the rulings' points (40, 15) and (40, 25):
/// on the ruling within the notch's z (a corner at its ends), on the face beyond it.
///
/// Two builds. **Through** (caps at −10 and 40): a band between the caps' whole circles
/// with the plate's notch as its one hole (the plate-side half-circle × z ∈ (0, 20), two
/// arcs on the cut rims and two rulings on the wall x = 40). **Staircase** (caps at 10 and
/// 40, the lower cap inside the plate): the lower boundary is a chain — the arc at z = 10
/// outside the plate, the arc at z = 20 inside, the two rulings between — so at the
/// station (40, 25) one arc ends `hi` and the other starts `lo`: the corner rule's one
/// discriminating population. ☑ In a notch both arcs share their ends, and «both ends
/// count» is invisible there — measured; the staircase is why the second build exists.
/// The other station, (40, 15), is the seam (`add_cylinder`'s `ref_dir` is −y): a
/// seam-incident root with an arc above it is the one tie the loops road keeps.
///
/// ☑ Measured on the through build's rays: on the two whole-circle rims alone, «opposite
/// sides of the two rim planes» and the loops road's «exactly one rim above» agree on all
/// 2,180 rays — every root, the graze on a rim included — which is why no banded arm
/// exists.
#[test]
fn the_lateral_parity_agrees_with_the_notch_it_bounds() {
    for staircase in [false, true] {
        lateral_lattice(staircase);
    }
}

fn lateral_lattice(staircase: bool) {
    use nacre_scalar::quad::{CylinderMeet, plane_plane_cylinder, plane_side};
    use nacre_scalar::{Orient, Rat};
    let caps = if staircase {
        [10.0, 40.0]
    } else {
        [-10.0, 40.0]
    };
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss_z(caps[0], caps[1] - caps[0]);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, plane_faces, rows, _) =
        armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, caps);
    let fuse = crate::cyl_chart::emit_lateral(
        BoolKind::Fuse,
        &jd,
        &setup.cyls,
        &plane_faces,
        &curved,
        &rows,
    )
    .unwrap();
    assert_eq!(fuse.len(), 1);
    let cf = crate::boolean::comp_face(&jd, &setup.cyls, &fuse[0]).unwrap();
    let combinatorics::CompSurf::Cylinder(def) = &cf.surf else {
        panic!("a lateral")
    };
    let combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
        panic!("loops, got {:?}", cf.outer)
    };
    assert!(cf.inner.is_empty(), "a lateral's holes are among its loops");
    let rims: Vec<usize> = loops
        .iter()
        .filter_map(|l| match l {
            combinatorics::LateralLoop::Circle(c) => Some(*c),
            combinatorics::LateralLoop::Ring(_) => None,
        })
        .collect();
    if staircase {
        assert_eq!(rims.len(), 1, "the upper cap's whole circle");
        assert_eq!(loops.len(), 2, "and the chain below");
    } else {
        assert_eq!(rims.len(), 2, "the caps' two whole circles");
        assert_eq!(loops.len(), 3, "and the notch as one ring");
    }
    let (o, mm, r2) = (def.origin(), def.dir(), def.r2());
    let rat = |k: i128| Rat::new(k, 2).unwrap();
    let ri = Rat::from_int;
    let zero = ri(0);
    let neg = |v: Rat| zero.checked_sub(v).unwrap();
    let x40 = [ri(1), zero, zero, neg(ri(40))];
    let y20 = [zero, ri(1), zero, neg(ri(20))];
    let (cap_lo, cap_hi) = (ri(caps[0] as i128), ri(caps[1] as i128));
    // The notch's z range: the plate's own (0, 20) through the boss, or — with the lower
    // cap inside the plate — the chain's step from the cap (10) to the plate's top (20).
    let (n_lo, n_hi) = (if staircase { cap_lo } else { zero }, ri(20));
    // The truth on the cylinder, by the root's side of the wall and the ray's z. `seam`:
    // the root is the station (40, 15), the seam generator; a seam-incident root with an
    // arc above it is the tie the loops road keeps (`arc_span`'s `SeamRoot`; z is asked
    // first, so with no arc above the rims still decide).
    let truth = |x_side: Orient, z0: Rat, seam: bool| -> Option<bool> {
        if seam && z0 < n_hi {
            return None;
        }
        if z0 >= cap_hi {
            return if z0 == cap_hi { None } else { Some(false) };
        }
        if z0 < cap_lo {
            return Some(false);
        }
        if z0 == cap_lo {
            // The lower cap: through the plate a whole rim (a tie everywhere); in the
            // staircase an arc outside the plate only — a tie there, a corner on the wall,
            // and nothing inside, where the face starts at the plate's top.
            return match (staircase, x_side) {
                (true, Orient::Negative) => Some(false),
                _ => None,
            };
        }
        match x_side {
            // Outside the plate: the band, whole.
            Orient::Positive => Some(true),
            // Inside the plate's footprint: the notch (a hole, or the chain's step), its
            // arcs the boundary — both arcs through the plate, the upper one alone in the
            // staircase (its lower boundary there is the cap, handled above).
            Orient::Negative => {
                if (!staircase && z0 == n_lo) || z0 == n_hi {
                    None
                } else {
                    Some(!(z0 > n_lo && z0 < n_hi))
                }
            }
            // On the wall: a ruling for z within the notch (corners at its ends), the
            // face beyond.
            Orient::Zero => {
                if z0 >= n_lo && z0 <= n_hi {
                    None
                } else {
                    Some(true)
                }
            }
        }
    };
    // on face, off (hole), boundary, tangent, station on-ruling, station on-face, seam ties
    let mut n = [0usize; 7];
    let mut ask = |pa: [Rat; 4], pb: [Rat; 4], z0: Rat, station: bool| {
        let roots = match plane_plane_cylinder(&pa, &pb, &o, &mm, r2).unwrap() {
            CylinderMeet::Pair { line, s } => (line, s),
            CylinderMeet::Tangent { .. } => {
                n[3] += 1;
                return;
            }
            other => panic!("an ordinary ray: {other:?}"),
        };
        let (line, s) = roots;
        for root in &s {
            let x_side = plane_side(&x40, &line, root);
            let seam = station && plane_side(&y20, &line, root) == Orient::Negative;
            let want = truth(x_side, z0, seam);
            let got = combinatorics::loop_parity(&jd, def, loops, &line, root);
            assert_eq!(
                got, want,
                "staircase {staircase} z0 {z0:?} side {x_side:?} station {station} seam {seam}"
            );
            match (got, station, seam) {
                (Some(true), false, _) => n[0] += 1,
                (Some(false), false, _) => n[1] += 1,
                (None, false, _) => n[2] += 1,
                (None, true, true) if z0 < n_lo || z0 > n_hi => n[6] += 1,
                (None, true, _) => n[4] += 1,
                (Some(_), true, _) => n[5] += 1,
            }
        }
    };
    for k in 0..=108 {
        let z0 = ri(-12).checked_add(rat(k)).unwrap();
        let pz = [zero, zero, ri(1), neg(z0)];
        for j in 0..=20 {
            let y0 = ri(15).checked_add(rat(j)).unwrap();
            ask([zero, ri(1), zero, neg(y0)], pz, z0, false);
        }
        ask(x40, pz, z0, true);
    }
    eprintln!(
        "lateral lattice (staircase {staircase}): on {} hole {} boundary {} tangent {} \
             station on-ruling {} station on-face {} seam ties {}",
        n[0], n[1], n[2], n[3], n[4], n[5], n[6]
    );
    assert!(
        n.iter().all(|&c| c > 0),
        "every arm has a population: {n:?}"
    );
}

/// ★ **A cut circle is a band boundary**. On the armed
/// through-boss, the per-class products of the four ⊥ classes are collected the production
/// way (`per_class` on each), and the chart must break the lateral at the two **cut**
/// circles (z = 0, z = 20) as well as the rims: three intervals, not one full-height band
/// (`boundary_lines` and `chart_of`'s z-lines agree on the four). The middle interval is
/// both-cut and is emitted as panel rings; blinding **both** of the chart's axes (the arc
/// labels and the rulings') leaves the reader no answer there and `emit_lateral` refuses by
/// name (`CylinderGateUndecided`) — blinding only the arcs does not, because the vertical
/// lines answer in their place.
/// The through-boss's curved carriers and plane faces, collected the production way
/// (`per_class` on each of the four ⊥ classes and the wall class), with the cylinder's rows
/// and the four ⊥ classes `[z0, z20, cap_lo, cap_hi]`. ★ This mirrors production's fold by
/// hand (`trace_result_faces`' accumulation). A drift between them is not caught by
/// anything: a test would simply start measuring a map the boolean never builds.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn armed_curved(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &crate::planes::PlaneSetup,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
    caps: [f64; 2],
) -> (
    Curved,
    Vec<LocalFace>,
    Vec<crate::bands::CylRow>,
    [usize; 4],
) {
    let z_class = |z: f64| -> usize {
        setup
            .geom
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[2] - z).abs() < 1e-12))
            .unwrap_or_else(|| panic!("a class at z = {z}"))
    };
    let (z0, z20, cap_lo, cap_hi) = (
        z_class(0.0),
        z_class(20.0),
        z_class(caps[0]),
        z_class(caps[1]),
    );
    // The production carriers: disk labels from every ⊥ class, cut rims from the cut ones.
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut plane_faces: Vec<LocalFace> = Vec::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    // ★ The wall class too: the chart's rulings are the wall's pieces, and without them the
    // both-cut middle is one whole-circle cell the emitter cannot read per sector.
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in [z0, z20, cap_lo, cap_hi, wc] {
        let edges = armed_class_edges(m, plate, boss, setup, jd, c, crossings);
        let staged = per_class(jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in &edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        plane_faces.extend(staged.faces);
    }
    let curved = Curved {
        aliases: Aliases::default(),
        disk_labels,
        arc_labels,
        cut_rims,
        rulings,
    };
    let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
    (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi])
}

#[test]
fn a_cut_circle_bounds_the_bands() {
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi]) =
        armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, [-10.0, 40.0]);
    let (disk_labels, arc_labels, cut_rims) =
        (&curved.disk_labels, &curved.arc_labels, &curved.cut_rims);
    assert!(cut_rims.contains_key(&(0, z0)) && cut_rims.contains_key(&(0, z20)));
    assert!(disk_labels.contains_key(&(0, cap_lo)) && disk_labels.contains_key(&(0, cap_hi)));
    assert_eq!(arc_labels[&(0, z0)].len(), 2, "two arcs, two sector labels");
    assert_eq!(arc_labels[&(0, z20)].len(), 2);
    assert_eq!(rows.len(), 1);
    // ★ A cut circle is a band boundary (the rulings ladder): the chart's boundary rule names
    // the two cut rims beside the caps, and the chart's own lines are exactly those four.
    let def = &setup.cyls[0].def;
    let t_of = |c: usize| crate::bands::axis_param(&jd, c, def).unwrap();
    let expect: Vec<Rat> = [cap_lo, z0, z20, cap_hi].map(t_of).to_vec();
    assert_eq!(
        crate::cyl_chart::boundary_lines(&jd, 0, def, &plane_faces, &curved, &rows).unwrap(),
        expect,
        "three intervals, cut circles included"
    );
    let chart = crate::cyl_chart::chart_of(&jd, &setup.cyls, 0, &plane_faces, &curved).unwrap();
    assert_eq!(
        chart.z_lines.iter().map(|l| l.t).collect::<Vec<_>>(),
        expect,
        "the chart's lines are the boundaries and nothing else here"
    );
    // ★ **The emitter answers with regions.** Under `keep` the wall splits the middle
    // interval into two sectors, one kept and one not, and the kept one joins the whole
    // bands above and below it: **fuse** emits one face — a `Band` between the two cap
    // circles with the unkept sector as its one **hole** (a 4-node ring: two arcs on the
    // cut rims' own nodes, two rulings); **cut** keeps the other sector alone — one 4-node
    // `Ring` face. Which sector is the outer one is the assembly's and the volume oracle's to
    // measure (the gate-opening cell), not this harness's to re-derive.
    let faces_for = |kind: BoolKind| -> Vec<LocalFace> {
        crate::cyl_chart::emit_lateral(kind, &jd, &setup.cyls, &plane_faces, &curved, &rows)
            .unwrap()
    };
    let (fuse, cut) = (faces_for(BoolKind::Fuse), faces_for(BoolKind::Cut));
    assert_eq!(fuse.len(), 1, "fuse: one lateral face, a band with a hole");
    assert_eq!(cut.len(), 1, "cut: one lateral face, the kept sector");
    let rim = &cut_rims[&(0, z0)];
    // The ccw arc a ring carries on the lower cut rim, as an ordered node pair.
    let arc_on_z0 = |r: &crate::boolean::Ring| -> [NodeId; 2] {
        let n = r.nodes.len();
        assert_eq!(n, 4, "two arcs and two rulings");
        let arcs = r
            .walls
            .iter()
            .filter(|w| matches!(w, crate::boolean::Wall::Arc { .. }))
            .count();
        assert_eq!(arcs, 2);
        for i in 0..n {
            let (a, b) = (r.nodes[i], r.nodes[(i + 1) % n]);
            if let crate::boolean::Wall::Arc { ccw, .. } = r.walls[i]
                && rim.nodes.contains(&a)
                && rim.nodes.contains(&b)
            {
                return if ccw { [a, b] } else { [b, a] };
            }
        }
        panic!("no arc on the lower cut rim");
    };
    let crate::boolean::Bound::Band { lo, hi } = &fuse[0].outer else {
        panic!("fuse: a band between the caps, got {:?}", fuse[0].outer);
    };
    assert!(
        matches!(
            (lo, hi),
            (
                crate::boolean::Rim::Circle(_),
                crate::boolean::Rim::Circle(_)
            )
        ),
        "the caps' whole circles are the band's rims"
    );
    assert_eq!(
        fuse[0].inner.len(),
        1,
        "the unkept sector is the band's one hole"
    );
    let crate::boolean::Bound::Ring(hole) = &fuse[0].inner[0] else {
        panic!("a hole is a ring");
    };
    let crate::boolean::Bound::Ring(panel) = &cut[0].outer else {
        panic!("cut: the kept sector is a ring, got {:?}", cut[0].outer);
    };
    assert!(cut[0].inner.is_empty());
    // ★ **Fuse's hole is cut's panel**: the sector fuse drops (inside the plate) is exactly
    // the sector cut keeps (the groove's wall), so the two rings carry the same ccw arc on
    // the lower cut rim — the old «complementary panels» claim, restated for regions.
    let (pf, pc) = (arc_on_z0(hole), arc_on_z0(panel));
    assert_eq!(pf, pc, "fuse's hole is cut's panel");
    // Negative control: without the sector labels the both-cut interval's cells have no
    // speaking end, and the emitter refuses the class rather than guess a chamber. Called
    // directly, not through the census — this hand-broken input violates the very premise
    // (`src0_present == 0`) the census asserts where the fact is made.
    let mut blind = Curved {
        aliases: Aliases::default(),
        disk_labels: curved.disk_labels.clone(),
        arc_labels: HashMap::new(),
        cut_rims: curved.cut_rims.clone(),
        rulings: curved.rulings.clone(),
    };
    blind.arc_labels.clear();
    // ★ **Blinding one axis is no longer blinding the reader** — the chart has two, and the
    // rulings answer where the rims are silent (the vertical read's own cell). So this now
    // says the weaker, truer thing: with the arcs gone the emitter still gets an answer, and
    // it is only when **both** axes are blinded that it refuses by name.
    assert!(
        crate::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &plane_faces,
            &blind,
            &rows
        )
        .is_ok(),
        "the vertical lines answer what the blinded rims cannot"
    );
    for v in blind.rulings.values_mut() {
        for r in v.iter_mut() {
            r.label = None;
        }
    }
    assert!(matches!(
        crate::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &plane_faces,
            &blind,
            &rows
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::CylinderGateUndecided,
            ..
        })
    ));
}

/// The armed through-boss fuse **assembled to the end of the road** — the pipeline the
/// cell-3 lock walks, extracted so the consumer locks (props' θ-range integrals, tess's
/// open-rim merge) measure the very same solid through one spelling: every class's
/// arrangement, the coplanar unify, the band/panel pass, the seam table, `reconstruct`.
/// Returns the model with the one welded solid live, plus the setup and the wall class.
fn armed_assembled_through_boss() -> (Model, crate::planes::PlaneSetup, usize) {
    let (mut m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut faces: Vec<LocalFace> = Vec::new();
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in 0..setup.geom.len() {
        let edges = armed_class_edges(&m, plate, boss, &setup, &jd, c, &crossings);
        let staged = per_class(&jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        // ★ This mirrors production's fold by hand (`trace_result_faces`' accumulation). A
        // drift between them is not caught by anything: the test would simply start measuring
        // a map the boolean never builds.
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in &edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        faces.extend(staged.faces);
    }
    let curved = Curved {
        aliases: Aliases::default(),
        disk_labels,
        arc_labels,
        cut_rims,
        rulings,
    };
    let faces = crate::boolean::unify_coplanar_faces(faces, &jd, &setup.cyls).unwrap();
    let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
    let mut faces = faces;
    faces.extend(
        crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, &setup.cyls, &faces, &curved, &rows)
            .unwrap(),
    );
    let seam = seam_table(&faces, &setup.cyls, &jd).unwrap();
    let out = crate::boolean::reconstruct(
        &mut m,
        &jd,
        &seam,
        &faces,
        &setup.cyls,
        &curved.cut_rims,
        None,
        crate::boolean::Tangencies::none(),
    )
    .unwrap();
    assert_eq!(out.len(), 1, "one welded solid");
    m.restore_live(out);
    m.rebuild_adjacency();
    (m, setup, wc)
}

/// ★ **The armed assembly welds the panels** (the lock): the whole
/// through-boss fuse, from production parts past the standing gate — every class's
/// arrangement, the coplanar unify, the band/panel pass, the seam table, and
/// `reconstruct` — comes back one solid with a clean `validate`. The ruling edges are
/// pinned structurally: exactly two straight lateral edges (the kept outer panel's), each
/// used exactly twice, carriers stated as **the edge's own fact** — the cylinder and the
/// wall plane — and the outer panel's seam-holding arc is split at an `OnSeam` vertex
/// (the panel road runs the wrap-arc split the Band road already had).
#[test]
fn the_armed_assembly_welds_the_panels() {
    let (m, setup, wc) = armed_assembled_through_boss();
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    // The ruling edges, structurally: straight lateral edges = carriers {cylinder, a plane}
    // with a Line curve (the seam's [cyl, cyl] spelling is excluded by the mixed pair).
    let cyl_surf = setup.cyls[0].surf;
    let wall_surf = setup.geom[wc].surf;
    let reach = m.reachable();
    let mut rulings = 0;
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        let mixed = (e.surfaces[0] == cyl_surf) != (e.surfaces[1] == cyl_surf);
        if !mixed {
            continue;
        }
        let curve = m.derive_edge_curve(e.surfaces, e.vertices).unwrap();
        if matches!(curve, nacre_geom::Curve::Line(_)) {
            rulings += 1;
            assert!(
                e.surfaces.contains(&wall_surf),
                "a ruling's plane carrier is the wall: {:?}",
                e.surfaces
            );
        }
    }
    assert_eq!(rulings, 2, "the kept outer panel's two rulings");
    // The seam split ran on the panel: an OnSeam vertex is reachable.
    let on_seam = reach
        .vertices
        .iter()
        .filter(|&&v| matches!(*m.vertex(v), nacre_topo::Vertex::OnSeam(_)))
        .count();
    assert!(on_seam >= 1, "the outer panel's wrap arc split at the seam");
}

/// ★ **The armed solid's mass properties are exact** (the props
/// instrument): `mass_props` on the assembled through-boss fuse, with the θ-range lateral
/// integrals live, answers the derived closed forms — volume `32000 + 1000π` (plate plus
/// the boss outside it), area `6200 + 425π` (walls 3000 + split x = 40 wall 600, caps
/// bitten `3200 − 25π`, boss disks `50π`, full bands `300π`, the outer half-panel `100π`).
/// The lifted probe's full-2π mis-answer (7849.34 = truth + the panel counted whole,
/// `+100π`) cross-checks the area derivation. This is the volume oracle standing while the
/// gate still refuses production input.
#[test]
fn the_armed_solids_mass_properties_are_exact() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let props = nacre_props::mass_props(&m, m.live_solids()[0]).unwrap();
    let volume = 32000.0 + 1000.0 * std::f64::consts::PI;
    let area = 6200.0 + 425.0 * std::f64::consts::PI;
    assert!(
        (props.volume - volume).abs() <= 1e-9 * volume,
        "volume {} != {volume}",
        props.volume
    );
    assert!(
        (props.area - area).abs() <= 1e-9 * area,
        "area {} != {area}",
        props.area
    );
}

/// ★ **The armed solid tessellates watertight** (the tess
/// instrument): the open-rim merge triangulates the θ-panels of the assembled through-boss
/// fuse — the outer panel's rims are two arcs split at the seam vertex, so this one solid
/// exercises the multi-polyline open chain *and* the seam-crossing unwrap — and every
/// undirected triangle edge is used exactly twice.
#[test]
fn the_armed_solid_tessellates_watertight() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    let open = uses.values().filter(|&&n| n != 2).count();
    assert_eq!(open, 0, "the mesh is watertight");
}

/// **The mixed parity reads a bitten ring — exactly, on the production pieces.**
///
/// The straddling boss cuts the plate-top ring at `(40, 15)` and `(40, 25)`, so that ring
/// carries two pierce corners and one arc step; the overhang digon is chord + outer arc.
/// Probes are derived, not read back: the chart for `n = +z` picks `e1 = [0, −1, 0]`, so
/// the ray runs toward −y at fixed x. `[37, 30]` is inside and its ray crosses the **arc
/// twice** (`(37−40)² + (y−20)² = 25` → y = 16, 24) before the bottom edge — the arc arm is
/// what that probe measures, and a chord-minded parity would answer it wrong. `[37, 18]`
/// sits inside the bite (outside the face), `[50, 20]` outside everything; `[43, 20]` is
/// inside the overhang (one arc crossing at y = 16). Cells are identified structurally,
/// not by size: the outside cell (winding −1) owns the outer arc piece's twin, so that
/// piece's forward cell is the overhang digon; the other digon is the bite.
#[test]
fn the_mixed_parity_reads_a_bitten_ring() {
    with_bitten_rings(|jd, cyls, wc, big, overhang, bite| {
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
        let rat = |x: i128, y: i128| {
            [
                nacre_scalar::Rat::from_int(x),
                nacre_scalar::Rat::from_int(y),
                nacre_scalar::Rat::from_int(20),
            ]
        };
        let ask = |ring: &[combinatorics::RingEdge], p: [nacre_scalar::Rat; 3]| {
            combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring)
        };
        assert_eq!(ask(big, rat(12, 12)), Some(true), "plain interior");
        assert_eq!(
            ask(big, rat(37, 30)),
            Some(true),
            "interior whose ray crosses the inner arc twice"
        );
        assert_eq!(ask(big, rat(37, 18)), Some(false), "inside the bite");
        assert_eq!(ask(big, rat(50, 20)), Some(false), "outside everything");
        assert_eq!(
            ask(overhang, rat(43, 20)),
            Some(true),
            "inside the overhang"
        );
        assert_eq!(
            ask(overhang, rat(37, 18)),
            Some(false),
            "the bite is not the overhang"
        );
        assert_eq!(ask(bite, rat(37, 18)), Some(true), "inside the bite digon");
        assert_eq!(
            ask(bite, rat(43, 20)),
            Some(false),
            "the overhang is not the bite"
        );
    });
}

/// ★★★★★ **The digon's truth is known independently, so the parity can be swept rather than
/// spot-checked.** A digon of a chord and an arc is `disk ∩ half-space`, and both halves are
/// exact rational predicates the arrangement already owns —
/// [`nacre_scalar::quad::cylinder_radial_side`] and the sign of the wall's plane equation. A
/// grid over the circle's neighbourhood therefore checks **every** answer, and it is the only
/// control here that crosses the straight arm, the arc arm **and their junction**: the
/// non-mixed road's oracle cannot see an arc, and a whole circle's arm collapses to
/// `Ordering::Equal` on a single step.
///
/// ☑ **What it does and does not reach, measured.** On this fixture the chord's ends are the
/// seam and the tangent columns, and a ray along the chord meets **both** arc ends at once —
/// so this sweep is green under a flipped arc-end departure sign (two ends on one ray flip
/// together and keep the parity) and green under the straight arm's old corner abstention
/// (negative controls, both measured). It is the oracle for the mixed road as a
/// whole — the straight arm, the arc arm and their junction on ordinary roots — and the
/// corner-bitten plate (`a_root_at_an_arc_end_is_a_corner_on_the_ray`) is the watch on the
/// arc-end arm, where a ray meets one end alone.
#[test]
fn the_mixed_parity_agrees_with_the_digon_it_bounds() {
    use nacre_scalar::{Orient, Rat};
    with_bitten_rings(|jd, cyls, wc, _big, overhang, bite| {
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
        let def = cyls
            .iter()
            .map(|c| &c.def)
            .find(|d| *d.r2() == nacre_scalar::BigRat::from(Rat::from_int(25)))
            .expect("the bitten circle");
        let (mut swept, mut on_boundary) = (0usize, 0usize);
        let (mut abstained, mut inside_seen) = (0usize, 0usize);
        let mut on_chord = 0usize;
        for i in 60..=100 {
            for j in 20..=60 {
                let p = [
                    Rat::new(i.into(), 2).unwrap(),
                    Rat::new(j.into(), 2).unwrap(),
                    Rat::from_int(20),
                ];
                let radial = nacre_scalar::quad::cylinder_radial_side(
                    &p,
                    &def.origin(),
                    &def.dir(),
                    def.r2(),
                );
                // The chord is the plate's wall `x = 40`; the bite keeps `x < 40`.
                let wall = p[0].checked_sub(Rat::from_int(40)).unwrap();
                if radial == Orient::Zero || wall == Rat::from_int(0) {
                    on_boundary += 1;
                    // ★★★★★ **A point strictly inside the circle and *on* the chord is on
                    // both digons' boundary, and "inside" has no answer there.** These are the
                    // grid's sharpest points: the ray along the chart's first axis leaves such
                    // a point straight through **both arc endpoints** — the chord is a step
                    // lying along the ray, and the straight arm names the probe between its
                    // ends as the boundary before any arc is asked.
                    if radial == Orient::Negative && wall == Rat::from_int(0) {
                        for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                            assert_eq!(
                                combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                None,
                                "{who} must abstain on its own chord at {p:?}"
                            );
                        }
                        on_chord += 1;
                    }
                    continue;
                }
                let inside_bite = radial == Orient::Negative && wall < Rat::from_int(0);
                let inside_over = radial == Orient::Negative && wall > Rat::from_int(0);
                for (ring, want, who) in [
                    (bite, inside_bite, "bite"),
                    (overhang, inside_over, "overhang"),
                ] {
                    // ★ **An abstention is allowed and a wrong answer is not.** The ray runs
                    // along the chart's first axis, so a grid point whose ray meets a ring
                    // corner has no parity to report — the caller's remedy is another point,
                    // which is exactly what `coord_probes` does with its candidate list. What
                    // the sweep locks is that every answer it *does* give is the truth.
                    match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                        Some(got) => {
                            assert_eq!(got, want, "{who} at {p:?}");
                            swept += 1;
                            inside_seen += usize::from(want);
                        }
                        None => abstained += 1,
                    }
                }
            }
        }
        // Neither vacuous nor all-outside: the grid straddles the circle and the chord, and
        // the two digons between them own a real interior.
        assert!(swept > 1_000, "answers swept: {swept}");
        assert!(inside_seen > 100, "points inside a digon: {inside_seen}");
        assert!(
            on_boundary > 0,
            "the grid meets the boundary: {on_boundary}"
        );
        // Recorded rather than bounded: every abstention here is the tangent ray (the rows
        // `y = 15, 25` graze the circle — measured: 160 of 160, no corner among
        // them), and its count is a property of this grid, not of the rule.
        assert!(abstained > 0, "abstentions: {abstained}");
        // What those abstentions were, by kind.
        {
            let rows = combinatorics::tie_probe::ROWS
                .lock()
                .expect("the probe's lock is never held across a panic");
            let me = std::thread::current().name().unwrap_or("?").to_string();
            let mut hist: Vec<(combinatorics::tie_probe::Tie, usize)> = Vec::new();
            for (_, t) in rows.iter().filter(|(n, _)| *n == me) {
                match hist.iter_mut().find(|(k, _)| k == t) {
                    Some((_, c)) => *c += 1,
                    None => hist.push((*t, 1)),
                }
            }
            eprintln!("P2 digon abstained {abstained} on_chord {on_chord} kinds {hist:?}");
        }
        assert!(
            on_chord > 0,
            "points on the chord inside the circle: {on_chord}"
        );
    });
}

/// The bitten-plate fixture, handed to `f` as `(judge, cyls, class, big, overhang, bite)`.
/// The bitten fixture: one wall (`x = 40`) cuts the circle into two digons.
fn with_bitten_rings(
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    with_cut_rings([40.0, 40.0, 20.0], false, f);
}

/// A plate `[0, plate]` with a z-cylinder (r = 5, h = 10) standing on its top at
/// `(40, 20)`: the class carrying the cut circle, its judge, and the three bounded rings —
/// the bitten top, the overhang (disk ∖ plate) and the bite (disk ∩ plate). The plate's
/// extent picks the cut: `[40, 40, 20]` bites with one wall (two digons); `[40, 24, 20]`
/// with two — the corner `(40, 24)` sits inside the circle and both rings are trigons
/// whose arc ends at the rational `(37, 24)`, off the seam and off the tangent columns.
/// `seam_off` states the cylinder exactly with `ref_dir = x`, so the seam sits at
/// `(45, 20)` — on no ring corner — instead of `add_cylinder`'s `−y`, which puts it on the
/// chord end `(40, 15)`.
fn with_cut_rings(
    plate: [f64; 3],
    seam_off: bool,
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(plate));
    let b = if seam_off {
        use nacre_scalar::Rat;
        let r = Rat::from_int;
        m.add_cylinder_exact(
            [r(40), r(20), r(20)],
            [r(0), r(0), r(1)],
            [r(1), r(0), r(0)],
            r(5),
            r(10),
            None,
        )
        .expect("an exact cylinder on an axis frame")
        .0
    } else {
        m.add_cylinder(
            Point3::from_array([40.0, 20.0, 20.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            10.0,
        )
    };
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
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    // The one class whose circle is cut: find it by running the split everywhere.
    let mut found = None;
    for wc in 0..planes.len() {
        if !matches!(plane_ix.get(wc), Some(ClassIx::Plane(_)) | None) && wc < plane_ix.len() {
            continue;
        }
        if combinatorics::class_coeffs_rat(&jd, wc).is_none() {
            continue;
        }
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
        let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()) else {
            continue;
        };
        let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
            continue;
        };
        let Ok(edges) = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
        else {
            continue;
        };
        if edges.arcs.is_empty() {
            continue;
        }
        let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
        let ns = edges.segs.len();
        let na = edges.arcs.len();
        // Four cells share this class: the bitten plate-top face (+1, many steps), the
        // outside (−1 — it also borders the arc, via the outer bulge), and two digons
        // (chord + arc): the overhang (disk minus plate) and the bite (disk ∩ plate).
        // Tell them apart structurally: the outside owns the outer piece's twin, so the
        // outer piece's forward cell is the overhang; the remaining digon is the bite.
        let outside = cells
            .iter()
            .position(|c| c.winding == -1)
            .expect("the outside cell");
        let bounded: Vec<usize> = (0..cells.len())
            .filter(|&i| cells[i].winding == 1)
            .collect();
        assert_eq!(
            bounded.len(),
            3,
            "the bitten top, the overhang and the bite"
        );
        let big_ix = *bounded
            .iter()
            .max_by_key(|&&i| cells[i].half_edges.len())
            .expect("the bitten top ring");
        let outer_ai = (0..na)
            .find(|ai| face_of[&(2 * (ns + ai) + 1)] == outside)
            .expect("the outer piece borders the outside");
        let overhang_ix = face_of[&(2 * (ns + outer_ai))];
        let bite_ix = *bounded
            .iter()
            .find(|&&i| i != big_ix && i != overhang_ix)
            .expect("the bite");
        assert_eq!(
            cells[overhang_ix].half_edges.len(),
            cells[bite_ix].half_edges.len(),
            "the same walls cut the overhang and the bite"
        );
        let ring = |ix: usize| -> Vec<combinatorics::RingEdge> {
            cells[ix]
                .half_edges
                .iter()
                .map(|&he| edges.edge_at(he))
                .collect()
        };
        found = Some((wc, ring(big_ix), ring(overhang_ix), ring(bite_ix)));
        break;
    }
    let (wc, big, overhang, bite) = found.expect("one class carries the cut circle");
    f(&jd, &cyls, wc, &big, &overhang, &bite);
}

/// **A root at an arc's own end is a corner on the ray, and the arc's tangent there says
/// whether the arc counts it** — the half-open rule in the arc arm, on the one
/// fixture that reaches it.
///
/// The bitten fixture cannot: its chord's ends are the seam and the tangent columns, and
/// each abstains first (measured over the whole suite before this cell: `ArcRootAtEnd` 0
/// while the straight arm's corner tie spoke for the same corners). The plate's second wall
/// `y = 24` puts the corner `(40, 24)` inside the circle, so the bite's arc ends at
/// `(37, 24)` — a 3-4-5 point, rational, off the seam `(40, 15)` and off the tangent
/// columns `x = 35, 45` — and the lattice column `x = 37` sends its rays through that end
/// **alone**, beside an ordinary second root at `(37, 16)`. That is the single-end crossing
/// a global sign error cannot hide: two ends on one ray flip together and keep the parity
/// (the bitten fixture's chord), one end on a ray does not. Every answer the sweep gives
/// is the truth, the boundary abstains, and the arm decided the end exactly once per point
/// on the shooting side of that column, for both rings — under both seam placements, so
/// the seam-incident arms and the `(false, false)` arm each decide an end.
#[test]
fn a_root_at_an_arc_end_is_a_corner_on_the_ray() {
    use nacre_scalar::{Orient, Rat};
    for seam_off in [false, true] {
        with_cut_rings(
            [40.0, 24.0, 20.0],
            seam_off,
            |jd, cyls, wc, _big, overhang, bite| {
                let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
                let def = cyls
                    .iter()
                    .map(|c| &c.def)
                    .find(|d| *d.r2() == nacre_scalar::BigRat::from(Rat::from_int(25)))
                    .expect("the cut circle");
                let decided0 = combinatorics::tie_probe::arc_end_decisions_here();
                let (mut swept, mut abstained, mut inside_seen, mut boundary) = (0usize, 0, 0, 0);
                let mut column = [0usize; 2];
                let rat = |k: i128| Rat::new(k, 2).unwrap();
                for i in 60..=100 {
                    for j in 20..=60 {
                        let p = [rat(i), rat(j), Rat::from_int(20)];
                        let radial = nacre_scalar::quad::cylinder_radial_side(
                            &p,
                            &def.origin(),
                            &def.dir(),
                            def.r2(),
                        );
                        let (x, y) = (p[0], p[1]);
                        let (wx, wy) = (x == Rat::from_int(40), y == Rat::from_int(24));
                        // Both rings' boundary: the circle and the two chords, `x = 40` for
                        // `15 ≤ y ≤ 24` and `y = 24` for `37 ≤ x ≤ 40`, ends included.
                        let on_chord = (wx && y >= Rat::from_int(15) && y <= Rat::from_int(24))
                            || (wy && x >= Rat::from_int(37) && x <= Rat::from_int(40));
                        if radial == Orient::Zero || on_chord {
                            boundary += 1;
                            for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                                assert_eq!(
                                    combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                    None,
                                    "{who} must abstain on its boundary at {p:?}"
                                );
                            }
                            continue;
                        }
                        let in_disk = radial == Orient::Negative;
                        let on_plate = x < Rat::from_int(40) && y < Rat::from_int(24);
                        for (k, (ring, want, who)) in [
                            (bite, in_disk && on_plate, "bite"),
                            (overhang, in_disk && !on_plate, "overhang"),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                                Some(got) => {
                                    assert_eq!(got, want, "{who} at {p:?}");
                                    swept += 1;
                                    inside_seen += usize::from(want);
                                    if x == Rat::from_int(37) && radial == Orient::Positive {
                                        column[k] += 1;
                                    }
                                }
                                None => abstained += 1,
                            }
                        }
                    }
                }
                assert!(swept > 1_000, "answers swept: {swept}");
                assert!(inside_seen > 100, "points inside a ring: {inside_seen}");
                assert!(boundary > 0, "the grid meets the boundary: {boundary}");
                // The `x = 37` column outside the circle: 24 lattice points, every one answered by
                // both rings — the shooting side through the end, the other side by a miss — and
                // the end decided exactly once per point on the shooting side: 12 × 2 rings.
                assert_eq!(column, [24, 24], "the single-end column is answered");
                let decided = combinatorics::tie_probe::arc_end_decisions_here() - decided0;
                eprintln!(
                    "arc-end sweep (seam_off {seam_off}): swept {swept} abstained {abstained} \
                 boundary {boundary} arc-end decisions {decided}"
                );
                // With the seam on the chord end `(40, 15)` (the `add_cylinder` build) the
                // `(true, false)` / `(false, true)` arms decide the `x = 37` column's shooting side:
                // 12 points × 2 rings. With the seam off every corner (`ref_dir = x`) both ends are
                // ordinary and the `(false, false)` arm decides — and the `x = 40` column's rays now
                // meet `(40, 15)` alone as well (its other root `(40, 25)` is the overhang's
                // interior), 12 more per ring, one of them a call that then abstains at the corner
                // `(40, 24)`; the column's seam-root abstentions are gone (182 → 160).
                assert_eq!(
                    decided,
                    if seam_off { 48 } else { 24 },
                    "the arc-end arm decided the single-end rays (seam_off {seam_off})"
                );
            },
        );
    }
}

/// **On a ring whose every corner is rational, the mixed road and the chart road are one
/// rule**: `point_in_ring_2d_rat` says Inside/Outside/OnBoundary and
/// `point_in_mixed_ring` `Some(true)`/`Some(false)`/`None`, and the pairs match at every
/// lattice point of every bounded ring of every class of two overlapping plates — squares
/// and L-shapes with reflex corners. A lattice at half steps shares a row with every corner
/// and every horizontal step, so every corner-on-the-ray configuration the half-open rule
/// distinguishes (both steps up, both down, one each, a step along the ray, the probe at
/// the corner) is on the sweep. Before this cell the mixed road abstained at every such
/// corner (5 abstentions where the chart road answered, over the whole suite's rational
/// rings; 0 disagreements — P3); a whole-suite shadow of the two roads cost 61 % of the
/// serial sweep and is not kept — this lattice is the lock.
#[test]
fn the_two_roads_agree_on_every_rational_ring() {
    use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
    use nacre_scalar::Rat;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([2.0, 2.0, 0.0]),
        Point3::from_array([6.0, 6.0, 2.0]),
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
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let (mut rings_swept, mut points, mut inside, mut boundary, mut reflex) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    for wc in 0..planes.len() {
        let Some(coeffs) = combinatorics::class_coeffs_rat(&jd, wc) else {
            continue;
        };
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
        let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()) else {
            continue;
        };
        let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
            continue;
        };
        let Ok(edges) = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
        else {
            continue;
        };
        let Ok((cells, _)) = walk_cells(&jd, &cyls, wc, &edges) else {
            continue;
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let chart = combinatorics::Chart2dRat::of_normal(&n).unwrap();
        let (e1, e2) = chart.axes();
        // Axis-aligned classes: unit chart axes and a unit normal, so a chart point
        // `(X, Y)` lifts to `X·e1 + Y·e2 − d·n`.
        for e in [e1, e2, &n] {
            assert_eq!(combinatorics::dot3_rat(e, e), Some(Rat::from_int(1)));
        }
        let lift = |xy: [Rat; 2]| -> [Rat; 3] {
            let mut p = [Rat::from_int(0); 3];
            for k in 0..3 {
                p[k] = xy[0]
                    .checked_mul(e1[k])
                    .unwrap()
                    .checked_add(xy[1].checked_mul(e2[k]).unwrap())
                    .unwrap()
                    .checked_sub(coeffs[3].checked_mul(n[k]).unwrap())
                    .unwrap();
            }
            p
        };
        for c in cells.iter().filter(|c| c.winding == 1) {
            let ring: Vec<combinatorics::RingEdge> =
                c.half_edges.iter().map(|&he| edges.edge_at(he)).collect();
            assert!(!combinatorics::ring_is_mixed(&ring), "a planar fixture");
            let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
            let ring2 = chart.ring(&jd, &nodes).expect("rational corners");
            if ring2.len() > 4 {
                reflex += 1;
            }
            let bound = |k: usize, lo: bool| -> i128 {
                let it = ring2.iter().map(|q| {
                    assert_eq!(q[k].denom(), 1, "integer corners");
                    q[k].numer()
                });
                if lo {
                    it.min().unwrap() - 1
                } else {
                    it.max().unwrap() + 1
                }
            };
            for xi in 2 * bound(0, true)..=2 * bound(0, false) {
                for yi in 2 * bound(1, true)..=2 * bound(1, false) {
                    let q = [Rat::new(xi, 2).unwrap(), Rat::new(yi, 2).unwrap()];
                    let p = lift(q);
                    assert_eq!(
                        chart.project(&p),
                        Some(q),
                        "the lift is the chart's inverse"
                    );
                    let plane = point_in_ring_2d_rat(q, &ring2);
                    let mixed = combinatorics::point_in_mixed_ring(&jd, &cyls, &coeffs, &p, &ring);
                    match (plane, mixed) {
                        (RingSide::Inside, Some(true)) => inside += 1,
                        (RingSide::Outside, Some(false)) => {}
                        (RingSide::OnBoundary, None) => boundary += 1,
                        other => panic!("class {wc} ring {ring2:?} at {q:?}: {other:?}"),
                    }
                    points += 1;
                }
            }
            rings_swept += 1;
        }
    }
    assert!(rings_swept >= 3, "rings swept: {rings_swept}");
    assert!(reflex > 0, "an L-shaped ring: {reflex}");
    assert!(
        inside > 0 && boundary > 0,
        "points {points} inside {inside} boundary {boundary}"
    );
    eprintln!("two roads: rings {rings_swept} points {points} inside {inside} boundary {boundary}");
}

/// The gated drill population's fixture: a `[0,2]³` box and an axis-aligned cylinder at
/// `(1,1)`, r=0.5 — every wall is a full unit from the axis, so the population gate passes
/// and [`plane_index_setup`] hands the arrangement bricks a cylinder-bearing table.
fn drilled(m: &mut Model, z0: f64, h: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = m.add_cylinder(
        Point3::from_array([1.0, 1.0, z0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        h,
    );
    m.rebuild_adjacency();
    (a, b)
}

/// The plane class whose defining triangle lies wholly at `z` — a z-cap.
fn class_at_z(planes: &[WorkingPlane], z: f64) -> usize {
    planes
        .iter()
        .position(|p| p.tri.iter().all(|q| (q.as_array()[2] - z).abs() < 1e-12))
        .expect("a z-cap class")
}

/// A cylinder cap alone on its class: the disk's circular outer traces as a **seated
/// circle** — and the lateral, whose rim lies on this very plane, adds its **graze** beside it.
/// The cell bricks turn the pair into a disk `+1` / contour `−1` pair whose labels say "body on
/// the cap's inside", with no segments anywhere. The emitted face's outer is the circle itself
/// ([`crate::boolean::Bound::Circle`]), the vocabulary the assembly consumes.
#[test]
fn a_cylinder_cap_is_a_seated_circle_and_a_disk_face() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // caps at z=-1 and z=3, clear of the box
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 3.0);

    let mut tr = Trace::default();
    trace_one_of(
        &m,
        b,
        SolidSide::B,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_b,
        &plane_ix,
        Default::default(),
        &mut tr,
    );
    assert!(
        tr.segs.is_empty(),
        "a circle owes the segment machinery nothing"
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    // The class's normal is the cap's own outward +z (the cap is its only member), so the
    // body — below z=3 — is `body_above: false`. ★ **Two contributions, not one**: the cap is
    // seated here and the lateral's rim grazes here, and on this convex cap they agree.
    assert!(planes[wc].plane.normal().as_array()[2] > 0.0);
    let mut kinds: Vec<SegKind> = tr
        .circles
        .iter()
        .inspect(|c| assert!(c.cyl == 0 && c.solid == SolidSide::B, "{:?}", tr.circles))
        .map(|c| c.kind)
        .collect();
    kinds.sort_by_key(|k| matches!(k, SegKind::Graze { .. }));
    assert!(
        matches!(
            kinds[..],
            [
                SegKind::Seated { body_above: false },
                SegKind::Graze { body_above: false },
            ]
        ),
        "{:?}",
        tr.circles
    );

    let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
    let edges = ClassEdges::of(&jd, &cyls, wc, &[], &circles, &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
    // Pseudo-half-edges 0 and 1 (no segments): the disk (+1) and its contour (−1).
    assert_eq!(cells.len(), 2, "{cells:?}");
    assert_eq!(
        (cells[0].half_edges.as_slice(), cells[0].winding),
        (&[0][..], 1)
    );
    assert_eq!(
        (cells[1].half_edges.as_slice(), cells[1].winding),
        (&[1][..], -1)
    );
    let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.root_groups,
        vec![1],
        "the contour bounds the unbounded region"
    );
    assert!(nesting.holes.is_empty());
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    assert_eq!(labels[1], [false; 4]);
    assert_eq!(
        labels[0],
        [false, false, false, true],
        "inside the circle, B below the plane only"
    );
    // Fuse keeps below and not above across the disk → the disk is a result face, and its
    // outer boundary is the circle — no ring, no nodes.
    let (out, _, _) = emit_faces(
        BoolKind::Fuse,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "{out:?}");
    assert!(matches!(
        out[0].outer,
        crate::boolean::Bound::Circle { cyl: 0 }
    ));
    assert!(out[0].inner.is_empty());
}

/// The through-drill's cap arrangement, brick by brick: the box cap is a seated 4-ring, the
/// lateral crosses the cap plane in a **transversal circle**, and the bricks nest the circle
/// as the cap cell's hole, label the disk with the cylinder straddling `W`, and emit — for
/// `Cut` — exactly one face: the cap with a circular hole. The per-kind `edge_mask`
/// difference is what the two labels measure (seated flips one side, transversal flips
/// both).
#[test]
fn a_drill_circle_is_a_hole_of_the_cap_ring() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // z∈[-1,3]: through both caps of the box
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 0.0);

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
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::B,
                kind: SegKind::Transversal { .. },
                // No hole in this band, so the mark is the whole circle.
                arc: None,
            }]
        ),
        "{:?}",
        tr.circles
    );

    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
    let split = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()).unwrap();
    assert_eq!(
        split.len(),
        4,
        "the cap ring alone — the circle owes the splitter nothing"
    );

    let edges = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 4, "cap ±1 and circle ±1: {cells:?}");
    let at = |he: usize| cells.iter().position(|c| c.half_edges == [he]).unwrap();
    let (disk, contour) = (at(2 * split.len()), at(2 * split.len() + 1));
    let cap = cells
        .iter()
        .position(|c| c.winding == 1 && c.half_edges.len() == 4)
        .unwrap();

    let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.holes.get(&cap).map(Vec::as_slice),
        Some(&[contour][..]),
        "the circle contour is the cap cell's hole"
    );

    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    // The class is the box's **bottom** cap, so which of `W`'s two sides carries the box is
    // read off the class normal rather than assumed: the box occupies z>0.
    let up = planes[wc].plane.normal().as_array()[2] > 0.0;
    let box_side = usize::from(!up); // 0 = above `n_W`, 1 = below
    let mut cap_label = [false; 4];
    cap_label[box_side] = true;
    assert_eq!(
        labels[cap], cap_label,
        "outside the circle: box material on the z>0 side only"
    );
    assert_eq!(
        labels[contour], labels[cap],
        "a hole bounds its host's region"
    );
    let mut disk_label = cap_label;
    disk_label[2] = true;
    disk_label[3] = true;
    assert_eq!(
        labels[disk], disk_label,
        "inside the circle the cylinder straddles W, the box side is unchanged"
    );

    let (out, _, _) = emit_faces(
        BoolKind::Cut,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "one face: the drilled cap — {out:?}");
    assert_eq!(out[0].outer.expect_ring().len(), 4);
    assert!(
        matches!(out[0].inner[..], [crate::boolean::Bound::Circle { cyl: 0 }]),
        "{:?}",
        out[0].inner
    );
}

/// The hole arm of the seated tracer: a face whose inner loop is a circle (a two-hole
/// plate's input shape — a bore `Cut` builds one in production now, so the hand-doctored
/// loops this arm was written against have company) emits its polygon segments **and** a
/// seated circle that inherits the face's own
/// body side.
#[test]
fn a_circular_hole_ring_traces_as_a_seated_circle() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0);
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
    let wc = class_at_z(&planes, 0.0);

    let input = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &[],
        Default::default(),
    );
    let (fp, fl) = input.faces[0]
        .iter()
        .find(|(fp, _)| matches!(plane_ix[*fp], ClassIx::Plane(c) if c == wc))
        .expect("the box's bottom cap sits on wc");
    let doctored = combinatorics::FaceLoops {
        outer: fl.outer.clone(),
        holes: Some(vec![combinatorics::LoopRing::Circle { cyl: 0 }]),
        cycles: None,
    };
    let mut tr = Trace::default();
    trace_one(
        &[(*fp, doctored)],
        SolidSide::A,
        wc,
        &jd,
        &[],
        &faces_tab,
        &plane_ix,
        &Default::default(),
        &Aliases::default(),
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(tr.segs.len(), 4, "the polygon outer still emits its edges");
    let SegKind::Seated { body_above } = tr.segs[0].kind else {
        panic!("a face on wc traces seated: {:?}", tr.segs[0]);
    };
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::A,
                kind: SegKind::Seated { body_above: ba },
                arc: None
            }]
            if ba == body_above
        ),
        "the hole circle inherits the face's body side: {:?}",
        tr.circles
    );
}

/// The existence condition, negatively: a cap plane **outside** the lateral's rim span gets
/// no circle (and no decline — a miss, like a parallel plane), while the rim-interior cap
/// still does, and a rim-**coincident** plane carries the cap's seated circle rather than a
/// transversal one. A ghost circle here would corrupt every label on the class.
#[test]
fn no_ghost_circle_outside_the_rim_span() {
    let mut m = Model::new();
    // z∈[-1,1.5]: through the box's bottom cap, short of its top cap at z=2.
    let (a, b) = drilled(&mut m, -1.0, 2.5);
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

    // The box's top cap at z=2 is past the rim span [−1, 1.5]: no circle, no decline.
    let top = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 2.0),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert!(top.circles.is_empty(), "{:?}", top.circles);
    assert!(top.declined.is_empty(), "{:?}", top.declined);

    // The bottom cap at z=0 is strictly inside the span: the transversal circle is there.
    let bottom = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 0.0),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert_eq!(
        bottom
            .circles
            .iter()
            .filter(|c| matches!(c.kind, SegKind::Transversal { .. }))
            .count(),
        1,
        "{:?}",
        bottom.circles
    );

    // The cylinder's own top cap at z=1.5 (inside the box): the rim-coincident plane carries
    // **two** contributions — the cap's `Seated` circle and the lateral's `Graze`, because a
    // rim is a touch, not a miss. It adds no *transversal* twin (t = span end, not strictly
    // inside). ★ On a **convex** cap the two agree about which side the body is on; the whole
    // point of carrying both is the reflex corner (a blind bore's ceiling) where they do not,
    // and `edge_mask`'s `Graze > Seated` then picks the right one.
    let cap = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 1.5),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert!(cap.declined.is_empty(), "{:?}", cap.declined);
    let sides: Vec<(bool, bool)> = cap
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some((false, body_above)),
            SegKind::Graze { body_above } => Some((true, body_above)),
            SegKind::Transversal { .. } | SegKind::Tangent { .. } => None,
        })
        .collect();
    assert_eq!(sides.len(), 2, "seated and graze: {:?}", cap.circles);
    assert!(
        sides.iter().any(|(g, _)| *g) && sides.iter().any(|(g, _)| !*g),
        "one of each kind: {:?}",
        cap.circles
    );
    assert_eq!(
        sides[0].1, sides[1].1,
        "a convex cap's two contributions agree on the body's side: {:?}",
        cap.circles
    );
}

/// ★★ **The corner the graze exists for.** On a convex cap the lateral's rim-graze and the
/// cap's seated circle say the same thing, so carrying both changes nothing — the test above
/// locks that. Here they **disagree**: at a blind bore's ceiling the plate's material is above
/// the cap while the bore's wall hangs below it, a reflex dihedral in the plane. `edge_mask`'s
/// `Graze > Seated` precedence then picks the wall's side, which is the whole mechanism of the
/// repair; before the graze existed, the seated rule flipped the wrong label bit and the next
/// boolean on that solid came back `CylinderGateUndecided`.
#[test]
fn at_a_blind_bores_ceiling_the_graze_and_the_seated_circle_disagree() {
    let mut m = Model::new();
    let plate = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let hole = m.add_cylinder(
        Point3::from_array([5.0, 5.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        3.0,
    );
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // A second operand so the arrangement runs on the bored solid as an operand.
    let boss = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 10.0]),
        Point3::from_array([2.0, 2.0, 12.0]),
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
    } = plane_index_setup(&m, bored, boss).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let tr = trace_on_class_of(
        &m,
        bored,
        boss,
        class_at_z(&planes, 3.0), // the bore's ceiling
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let seated: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    let grazes: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Graze { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    assert_eq!(
        seated.len(),
        1,
        "the ceiling is seated here: {:?}",
        tr.circles
    );
    assert_eq!(grazes.len(), 1, "the wall grazes here: {:?}", tr.circles);
    assert_ne!(
        seated[0], grazes[0],
        "a reflex dihedral: the cap and the wall put the body on opposite sides — this is the \
             only place the precedence matters, and the only reason the graze is emitted at all"
    );
    // The wall hangs below the ceiling, and that is the side `edge_mask` believes.
    assert!(!grazes[0], "the bore's wall is below its ceiling");
}

/// Partial overlap (E5) is NOT merged — different endpoints mean different edges. a and b share
/// the y=1 plane; a's y=1 chord is x∈[0,2], b's is x∈[1,3] — overlapping on x∈[1,2] but not
/// coincident. They must stay separate.
#[test]
fn partial_overlap_is_not_merged() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
    // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
    // seated≡transversal coincidence).
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    // The shared y=1 wall class (a face at y=1).
    let y1 = planes
        .iter()
        .position(|p| p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12))
        .expect("a y=1 face");
    // a's chord x∈[0,2] and b's chord x∈[1,3] ride y=1 but have different endpoints, so they
    // stay as two distinct MergedSegs. A merge that ignored extent would collapse them to one.
    let on_y1 = merged.iter().filter(|e| e.wall == y1).count();
    assert!(
        on_y1 >= 2,
        "a's and b's y=1 chords stay distinct (partial overlap not merged): {on_y1}"
    );
}

fn face_on_z1(fh: Handle<Face>, surf_ix: &HashMap<Handle<Face>, usize>, faces: &[FaceRow]) -> bool {
    let p = &faces[surf_ix[&fh]];
    // A cap in the z=1 plane: all three defining points at z=1.
    p.plane()
        .tri
        .iter()
        .all(|q| (q.as_array()[2] - 1.0).abs() < 1e-12)
}

/// Axis-aligned cubes never decline: every face is seated, a clean transversal chord, or a
/// parallel miss. This is falsifiable — a `declined` entry would mean the producer hit a
/// degeneracy it should not on this input.
#[test]
fn axis_aligned_cubes_decline_nothing() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
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
    for wc in 0..planes.len() {
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
        assert!(tr.declined.is_empty(), "class {wc} declined: {tr:?}");
    }
}

fn tilt_by(m: &mut Model, s: Handle<Solid>, deg: nacre_scalar::Rat) -> Handle<Solid> {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let out = transform(
        m,
        s,
        &Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(deg).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    out
}

/// **The crossing collector's direction families really are equivalence classes.**
///
/// ★ `split_at_crossings` replaced a predicate per wall pair with "same family?", which is sound
/// only because `wc ∩ w ∥ wc ∩ r` is transitive. The argument is in that function; this is the
/// **check**, over every wall pair of every class of a rotated and an axis-aligned fold: the
/// predicate's own answer must agree with the partition, everywhere.
///
/// It matters that this is a test and not a spike. A violation would not reject or panic — it
/// would invent a crossing point where the lines never meet, or lose one where they do, and the
/// census would drift by one vertex somewhere. **A silent wrong answer is exactly what a corpus
/// this size can hide**, so the invariant is asserted rather than inspected.
///
/// Separate from the timing spike on purpose: this asks the predicate about every pair, which
/// fills `WitnessPoint`'s realization cells and would make the phase timers read the collector 1.76x
/// cheaper than it is.
#[test]
#[ignore = "slow: every wall pair of every class, two folds"]
fn direction_families_partition_the_walls() {
    for rotated in [true, false] {
        let n = 24i128;
        let mut m = Model::new();
        let mut acc = if rotated {
            m.add_cuboid(
                Point3::from_array([-3.0, -3.0, 0.0]),
                Point3::from_array([3.0, 3.0, 2.0]),
            )
        } else {
            m.add_cuboid(
                Point3::from_array([-1.0, -1.0, 0.0]),
                Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
            )
        };
        m.rebuild_adjacency();
        let mut pairs = 0usize;
        for i in 0..n {
            let fin = if rotated {
                let f = m.add_cuboid(
                    Point3::from_array([2.0, -0.4, 0.0]),
                    Point3::from_array([8.0, 0.4, 1.0]),
                );
                m.rebuild_adjacency();
                tilt_by(&mut m, f, nacre_scalar::Rat::new(360 * i, n).unwrap())
            } else {
                let x = i as f64 * 0.5;
                let f = m.add_cuboid(
                    Point3::from_array([x, 0.5, 0.0]),
                    Point3::from_array([x + 0.2, 4.0, 2.0]),
                );
                m.rebuild_adjacency();
                f
            };
            // The audit needs the *same* judging context the boolean uses, and the walls of each
            // class as the tracer finds them — so it re-runs the setup and the trace, then checks
            // the partition the collector would build.
            let setup = plane_index_setup(&m, acc, fin).expect("setup");
            let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
            let trace_in = combinatorics::trace_input(
                &m,
                [(acc, &setup.inc_a), (fin, &setup.inc_b)],
                &setup.surf_ix,
                setup.planes.len(),
                &jd,
                &setup.plane_ix,
                &setup.cyls,
                Default::default(),
            );
            for wc in 0..setup.geom.len() {
                let tr = trace_on_class(
                    &trace_in,
                    wc,
                    &jd,
                    &setup.cyls,
                    &setup.planes,
                    &setup.plane_ix,
                    &Aliases::default(),
                );
                let mut walls: Vec<usize> = Vec::new();
                for s in &tr.segs {
                    if !walls.contains(&s.wall) {
                        walls.push(s.wall);
                    }
                }
                let par = |a: usize, b: usize| jd.plane_pair_dir_sign(wc, a, b) == 0;
                // The partition, exactly as `split_at_crossings` builds it.
                let mut reps: Vec<usize> = Vec::new();
                let dir: Vec<usize> = walls
                    .iter()
                    .map(|&w| {
                        reps.iter().position(|&rep| par(rep, w)).unwrap_or_else(|| {
                            reps.push(w);
                            reps.len() - 1
                        })
                    })
                    .collect();
                for (i, &a) in walls.iter().enumerate() {
                    for (j, &b) in walls.iter().enumerate().skip(i + 1) {
                        pairs += 1;
                        assert_eq!(
                            par(a, b),
                            dir[i] == dir[j],
                            "class {wc}: walls {a} and {b} — the predicate and the partition \
                                 disagree, so parallelism is not transitive here"
                        );
                    }
                }
            }
            acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
                .expect("fuse")
                .0[0];
            m.rebuild_adjacency();
        }
        println!(
            "  {} fold: {pairs} wall pairs audited, all agree",
            if rotated { "rotated" } else { "axis-aligned" }
        );
        assert!(
            pairs > 10_000,
            "{} fold audited only {pairs} pairs — too few to mean anything",
            if rotated { "rotated" } else { "axis-aligned" }
        );
    }
}

/// **S2a: where does a whole boolean's time go?** Every earlier profile answered a share *of a
/// phase* — and the one that mattered was never taken.
///
/// ★ Two ways to be wrong that this is built against:
///
/// 1. **Mixed clocks.** The last breakdown summed per-thread timers (CPU) and read the result
///    against the fold's wall clock. The part came out larger than the whole — 3,719ms of a
///    2,540ms fold — and the ratio it implied was meaningless. **Run this with
///    `--no-default-features`**, where every phase and the total are the same clock.
/// 2. **Phases that do not add up.** The sum is printed against the measured whole, so anything
///    unaccounted for shows as a gap rather than hiding inside a phase's share.
/// 3. ★★★ **Another test adding to the same counters.** `phase::` are process-global atomics
///    and `cargo test` runs this binary's tests on parallel threads, so `reset()` here does not
///    fence anything: every concurrently running test that calls `boolean` — and in an
///    `--ignored` pass that is the *other* spikes and the two rotation stress tests — lands in
///    the same buckets, while `whole` below is this thread's wall clock alone. The share
///    percentages then exceed 100 for reasons that have nothing to do with the code.
///    **Run it alone**: `cargo test -p nacre-ops --no-default-features spike_where -- --ignored
///    --nocapture --test-threads=1`.
///
/// Timers live on the production path (`phase::` in this module), not in a replica of it.
#[test]
#[ignore = "spike"]
fn spike_where_the_boolean_spends_it() {
    for rotated in [true, false] {
        spend(60, rotated);
    }
}

/// One fold's phase breakdown. `rotated` picks which predicate routes the judgements take:
/// an axis-aligned fold answers on the exact path, a rotated one mostly on the certified one,
/// and the difference between the two breakdowns is what the certification actually costs.
fn spend(n: i128, rotated: bool) {
    let mut m = Model::new();
    // ★ **The two folds must be the same *shape* of work, not the same model minus a rotation.**
    // Dropping the tilt stacks all 60 fins on one another — 7.6x fewer segments, a degenerate
    // model, and a comparison that says nothing. The axis-aligned arm places each fin at its
    // own x instead, so both arms fuse `n` distinct blades onto a growing solid and the only
    // difference is which predicate route the judgements take.
    let mut acc = if rotated {
        m.add_cuboid(
            Point3::from_array([-3.0, -3.0, 0.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        )
    } else {
        m.add_cuboid(
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
        )
    };
    m.rebuild_adjacency();
    let blade = |m: &mut Model, i: i128| {
        let f = if rotated {
            m.add_cuboid(
                Point3::from_array([2.0, -0.4, 0.0]),
                Point3::from_array([8.0, 0.4, 1.0]),
            )
        } else {
            let x = i as f64 * 0.5;
            m.add_cuboid(
                Point3::from_array([x, 0.5, 0.0]),
                Point3::from_array([x + 0.2, 4.0, 2.0]),
            )
        };
        m.rebuild_adjacency();
        if rotated {
            tilt_by(m, f, nacre_scalar::Rat::new(360 * i, n).unwrap())
        } else {
            f
        }
    };
    // Warm the code paths, then zero the counters: the first boolean pays for lazily-built
    // caches that the other n do not.
    {
        let f = blade(&mut m, if rotated { 1 } else { n });
        acc = super::boolean(&mut m, BoolKind::Fuse, acc, f)
            .expect("warm")
            .0[0];
        m.rebuild_adjacency();
    }
    phase::reset();
    phase::scale::reset();

    let mut whole = std::time::Duration::ZERO;
    for i in 0..n {
        let fin = blade(&mut m, i);
        let t = std::time::Instant::now();
        acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
            .expect("fuse")
            .0[0];
        whole += t.elapsed();
        m.rebuild_adjacency();
    }

    let rows = phase::all();
    // ★ Indentation is nesting: depth 0 are `boolean`'s own phases, depth 1 the inside of
    // `trace_result_faces`, depth 2 the inside of `split_at_crossings`. Depth 1 partitions the
    // work depth 0 does not name, so 0 and 1 add to the whole — and adding depth 2 on top would
    // count `split_at_crossings` twice.
    let depth = |l: &str| (l.len() - l.trim_start().len()) / 2;
    let sum: u64 = rows
        .iter()
        .filter(|(l, _)| depth(l) <= 1)
        .map(|(_, ns)| ns)
        .sum();
    let total = whole.as_nanos() as u64;
    println!(
        "\n{n}-fin fold ({}), whole boolean, serial build:",
        if rotated { "rotated" } else { "axis-aligned" }
    );
    for (label, ns) in &rows {
        println!(
            "  {label:<32} {:>8.1?}  {:>5.1}%",
            std::time::Duration::from_nanos(*ns),
            100.0 * *ns as f64 / total as f64
        );
    }
    println!(
        "  {:<32} {:>8.1?}  {:>5.1}%   ← accounted",
        "sum of the above",
        std::time::Duration::from_nanos(sum),
        100.0 * sum as f64 / total as f64
    );
    println!("  {:<32} {whole:>8.1?}  100.0%   ← measured", "the fold");
    println!(
        "\n  ★ split_at_crossings is {:.1}% of the whole boolean. The plan continues at 25%.",
        100.0
            * rows
                .iter()
                .find(|(l, _)| l.trim() == "split_at_crossings")
                .map(|(_, ns)| *ns)
                .unwrap_or(0) as f64
            / total as f64
    );

    // ★ **Scale, which decides whether the structural fixes are worth building.** A hull test
    // over wall *pairs* replaces `2 × |segs on r|` predicate calls with `2` — worth nothing when
    // a wall carries one segment, and worth the ratio when it carries many.
    use phase::scale as sc;
    let (segs, walls) = (sc::get(&sc::SEGS), sc::get(&sc::WALLS));
    println!("\n  scale, summed over every class of every boolean:");
    println!(
        "    segments / walls        {segs:>10} / {walls:<10} = {:.2}   ← S3b lives or dies here",
        segs as f64 / walls.max(1) as f64
    );
    println!("    split points (reps)     {:>10}", sc::get(&sc::PTS));
    println!(
        "    (1) collect trips       {:>10}",
        sc::get(&sc::COLLECT_TRIPS)
    );
    println!(
        "    (3) cover trips         {:>10}   = {:.2}x the collecting loop",
        sc::get(&sc::COVER_TRIPS),
        sc::get(&sc::COVER_TRIPS) as f64 / sc::get(&sc::COLLECT_TRIPS).max(1) as f64
    );
}

/// S3a gate: **what did hoisting the loop invariant do to the evidence?** `plane_pair_dir_sign`
/// records on the rotated path, and `BoolReport::coincidences` is a *count*, so asking the same
/// question fewer times moves it. The answer must not move; the count may.
#[test]
#[ignore = "spike"]
fn spike_report_after_hoisting() {
    let n = 24i128;
    let mut m = Model::new();
    let mut acc = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let (mut total, mut loosest) = (0usize, String::new());
    for i in 0..n {
        let fin = m.add_cuboid(
            Point3::from_array([2.0, -0.4, 0.0]),
            Point3::from_array([8.0, 0.4, 1.0]),
        );
        m.rebuild_adjacency();
        let fin = tilt_by(&mut m, fin, nacre_scalar::Rat::new(360 * i, n).unwrap());
        let (solids, report) =
            crate::boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
        total += report.coincidences;
        if let Some(e) = &report.loosest {
            loosest = format!("{e:?}");
        }
        acc = solids[0];
        m.rebuild_adjacency();
    }
    let v = nacre_props::mass_props(&m, acc).unwrap().volume;
    println!("\n  coincidences over {n} booleans : {total}");
    println!("  loosest (last)                : {loosest}");
    println!("  volume                        : {v:.9}");
}

/// B0: what is actually left to save, **in production's configuration**?
///
/// The earlier phase timing used `ClassReuse::Off`, so it counted work production never does —
/// `reuse.rs` skips most classes outright. This counts what survives that, how much of it a
/// per-class bounding box would cull, and how often the two things that would make the cull
/// unsound actually occur. Counting only; nothing is built and nothing is timed (the machine is
/// not quiet).
#[test]
#[ignore = "spike"]
fn spike_cull_potential() {
    let n = 60i128;
    let mut m = Model::new();
    let mut acc = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let mut tot = Counts::default();
    for i in 0..n {
        let fin = m.add_cuboid(
            Point3::from_array([2.0, -0.4, 0.0]),
            Point3::from_array([8.0, 0.4, 1.0]),
        );
        m.rebuild_adjacency();
        let fin = tilt_by(&mut m, fin, nacre_scalar::Rat::new(360 * i, n).unwrap());
        tot.add(&count_cull(&m, acc, fin));
        acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
            .expect("fuse")
            .0[0];
        m.rebuild_adjacency();
    }
    println!("over the fold, production config (ClassReuse::Proved):");
    println!("  classes                 {:>9}", tot.classes);
    println!(
        "  arranged after reuse    {:>9}  ({:.1}%)",
        tot.arranged,
        100.0 * tot.arranged as f64 / tot.classes as f64
    );
    println!(
        "  (face, class) pairs     {:>9}   in arranged classes",
        tot.pairs
    );
    println!(
        "  ★ cullable by box       {:>9}  ({:.1}%)",
        tot.cullable,
        100.0 * tot.cullable as f64 / tot.pairs.max(1) as f64
    );
    println!(
        "  seated footprint is one piece: {} of {} arranged classes",
        tot.one_piece, tot.arranged
    );
    // ★ Recorded, and **not** the operative number. The plan restricted the cull to classes
    // whose seated footprint is one connected piece, fearing that several pieces would need
    // several seeds. They do not: what the correction needs is that the culled chords' parity
    // be *constant over the footprint*, and the cull criterion (`box(F)` disjoint from the
    // footprint box) already guarantees no culled chord enters that box. A box is connected, so
    // the parity is one constant however many pieces the seated faces form.
    println!(
        "  (one-piece classes only:  {:>6}  of {} — not the limit, see the note)",
        tot.cullable_1p, tot.pairs_1p
    );
    println!(
        "  (class, solid) needing a seed: {}  of {}   — of those, {} are free (solid box misses)",
        tot.needs_seed,
        2 * tot.arranged,
        tot.seed_free
    );
}

#[derive(Default)]
struct Counts {
    classes: usize,
    arranged: usize,
    pairs: usize,
    cullable: usize,
    one_piece: usize,
    needs_seed: usize,
    /// Pairs in classes whose footprint is one piece — the realistic cull, since the
    /// multi-piece ones cannot take a single seed.
    pairs_1p: usize,
    cullable_1p: usize,
    /// Of `needs_seed`, the ones where the solid's whole box misses the footprint, so it
    /// **cannot** enclose it and the seed is false without any work.
    seed_free: usize,
}

impl Counts {
    fn add(&mut self, o: &Counts) {
        self.classes += o.classes;
        self.arranged += o.arranged;
        self.pairs += o.pairs;
        self.cullable += o.cullable;
        self.one_piece += o.one_piece;
        self.needs_seed += o.needs_seed;
        self.pairs_1p += o.pairs_1p;
        self.cullable_1p += o.cullable_1p;
        self.seed_free += o.seed_free;
    }
}

fn count_cull(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> Counts {
    let mut c = Counts::default();
    let Ok(setup) = plane_index_setup(m, a, b) else {
        return c;
    };
    let PlaneSetup {
        planes: faces_tab,
        n_a,
        geom,
        plane_ix,
        class_owner,
        standard,
        notes,
        ..
    } = setup;
    let jd = Judge::new(&geom, standard, &notes);
    let _ = &jd;
    let plans = crate::reuse::class_plans(
        m,
        crate::reuse::ClassReuse::Proved,
        BoolKind::Fuse,
        a,
        b,
        &geom,
        &class_owner,
    );
    let boxes: Vec<[[f64; 2]; 3]> = faces_tab
        .iter()
        .map(|f| face_box(m, f.face().expect("real")))
        .collect();
    let overlap = |x: &[[f64; 2]; 3], y: &[[f64; 2]; 3]| {
        (0..3).all(|k| x[k][0] <= y[k][1] && y[k][0] <= x[k][1])
    };
    c.classes = geom.len();
    for (wc, plan) in plans.iter().enumerate() {
        if *plan != crate::reuse::ClassPlan::Arrange {
            continue;
        }
        c.arranged += 1;
        let seated: Vec<usize> = (0..faces_tab.len())
            .filter(|&f| plane_ix[f].plane() == wc)
            .collect();
        if seated.is_empty() {
            continue;
        }
        // The class's region of interest, and whether it is one connected piece (box graph).
        let mut foot = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
        for &s in &seated {
            for k in 0..3 {
                foot[k][0] = foot[k][0].min(boxes[s][k][0]);
                foot[k][1] = foot[k][1].max(boxes[s][k][1]);
            }
        }
        let mut parent: Vec<usize> = (0..seated.len()).collect();
        for i in 0..seated.len() {
            for j in (i + 1)..seated.len() {
                if overlap(&boxes[seated[i]], &boxes[seated[j]]) {
                    let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                    parent[ri] = rj;
                }
            }
        }
        let pieces = (0..seated.len())
            .map(|i| uf_find(&mut parent, i))
            .collect::<std::collections::HashSet<_>>()
            .len();
        if pieces == 1 {
            c.one_piece += 1;
        }
        // How many of this class's traced faces the box would cull, and whether a solid that
        // only crosses `W` has any culled face (that is the one needing a seed).
        let mut culled_side = [false; 2];
        let mut seated_side = [false; 2];
        for &s in &seated {
            seated_side[usize::from(s >= n_a)] = true;
        }
        for f in 0..faces_tab.len() {
            c.pairs += 1;
            if pieces == 1 {
                c.pairs_1p += 1;
            }
            if plane_ix[f].plane() == wc {
                continue; // seated: never culled
            }
            if !overlap(&boxes[f], &foot) {
                c.cullable += 1;
                if pieces == 1 {
                    c.cullable_1p += 1;
                }
                culled_side[usize::from(f >= n_a)] = true;
            }
        }
        // Each operand's whole box: if it misses the footprint entirely it cannot enclose it,
        // so its parity is false for free.
        let mut solid_box = [[[f64::INFINITY, f64::NEG_INFINITY]; 3]; 2];
        for (f, fb) in boxes.iter().enumerate() {
            let side = usize::from(f >= n_a);
            for k in 0..3 {
                solid_box[side][k][0] = solid_box[side][k][0].min(fb[k][0]);
                solid_box[side][k][1] = solid_box[side][k][1].max(fb[k][1]);
            }
        }
        for side in 0..2 {
            if culled_side[side] && !seated_side[side] {
                c.needs_seed += 1;
                if !overlap(&solid_box[side], &foot) {
                    c.seed_free += 1;
                }
            }
        }
    }
    c
}
