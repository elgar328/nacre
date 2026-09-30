//! Coplanar contact — flush faces, shared caps, hollow parts, coplanar seams.

#![allow(unused_imports)]
use crate::common::*;
use nacre_exact::Axis;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean,
    replay,
};
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Orientation, Shell, Solid, Vertex};

/// A holed operand already survives the arrangement — it names holes as inner rings rather
/// than guarding against them, so the pocket rides through untouched. Nothing tested it.
/// Pin it before the guard comes down.
#[test]
fn containment_boolean_already_keeps_a_pocket() {
    for (kind, want) in [(BoolKind::Cut, 0.919), (BoolKind::Fuse, 0.92)] {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.05, 0.05, 0.05]),
            Point3::from_array([0.15, 0.15, 0.15]),
        );
        let r = boolean_one(&mut m, kind, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} {vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - want).abs() < 1e-9, "{kind:?} volume {vol}");
    }
}

#[test]
fn cut_with_a_hollow_operand_far_from_the_void() {
    // A cavitied operand whose seam misses the void: the corner cut is far from
    // the [1,2]³ void, so the void is carried through and preserved.
    // Every shell is walked, so the void does not vanish (dropping it reads as convex and
    // returns vol 26.875, cavities 0): 25.875 (27 − 1 void − 0.125 corner) with the cavity
    // intact.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    assert_eq!(m.solid(hollow).cavities.len(), 1);
    m.rebuild_adjacency();
    let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, cutter).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solid(r).cavities.len(), 1);
}

#[test]
fn a_slab_splits_a_hollow_box_into_two() {
    // A slab cut through the whole box (and its void) severs it into two solids.
    // The slab spans the full cross-section, so it opens the void — both pieces are
    // cavity-free.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let slab = m.add_cuboid(
        Point3::from_array([-0.5, 1.4, -0.5]),
        Point3::from_array([3.5, 1.6, 3.5]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    for &s in &solids {
        assert_eq!(m.solid(s).cavities.len(), 0);
    }
    // Hollow 26 (= 27 − 1 void); the slab removes 1.6 of material (8 area × 0.2 thick).
    let vol: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((vol - 24.4).abs() < 1e-9, "total volume {vol}");
}

// A hollow operand in a COPLANAR contact. The cavity goldens above take transversal cuts and
// every coplanar golden uses solid operands, so this is the intersection of the two. Emitting
// only outer-shell faces loses the void: the fuse reads 27.0625 — the *un-hollowed* cube plus the
// boss — with `cavities: 0` and a clean `validate`, because what remains is still a closed
// shell. Silent-wrong, invisible to every guard; every shell is walked.
#[test]
fn a_hollow_part_takes_a_coplanar_boss() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    // Top-flush boss on z=3: a genuine coplanar contact, clear of the void's planes.
    let boss = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 3.0]),
        Point3::from_array([0.75, 0.75, 4.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    // 27 − 1 void + 0.25·0.25·1 boss.
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 26.0625).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solid(r).cavities.len(), 1, "the void survives");
}

#[test]
fn a_hollow_part_takes_a_coplanar_pocket() {
    // The Cut twin of the boss case: a top-flush pocket sunk into a hollow part. Same blind
    // spot, same silent-wrong with outer shells only (26.96875 with the void gone).
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let tool = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 2.5]),
        Point3::from_array([0.75, 0.75, 3.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, tool).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    // 27 − 1 void − 0.25·0.25·0.5 pocket.
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.96875).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solid(r).cavities.len(), 1, "the void survives");
}

#[test]
fn a_coplanar_boss_over_a_void_plane_is_solved() {
    // The boss straddles the void's x=1 and y=1 planes (it spans [0.9,1.1]²), so classifying
    // the void's walls against it is not the clean whole-face case. This was pinned as an
    // honest reject with the standing instruction that a future change may "turn it into a
    // correct 26.04, never a silent answer" — and family #3 did: once one geometric plane is
    // one class whatever the two faces' sizes, the arrangement solves it. 26 (hollow) + 0.04
    // (boss); area 60 + 0.8 (boss sides) + 0.04 (its top) − 0.04 (its footprint); void intact.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let boss = m.add_cuboid(
        Point3::from_array([0.9, 0.9, 3.0]),
        Point3::from_array([1.1, 1.1, 4.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - 26.04).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 60.8).abs() < 1e-9, "area {}", props.area);
    assert_eq!(m.solid(r).cavities.len(), 1, "the void survives");
}

#[test]
fn a_hollow_part_takes_a_second_far_cut() {
    // The headline: keep cutting a part after it is hollow. A bore at the corner
    // opposite the void — the void survives, and the result is fed back as an
    // operand (chaining past the first cavity-producing op).
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    let bore = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, bore).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solid(r).cavities.len(), 1);
}

#[test]
fn a_non_convex_pad_cantilevers_and_runs_flush() {
    // **Three hard properties at once**, which is what makes this footprint worth keeping. The
    // frame maps local `[px, py]` to world `(0.5 + py, 0.5 − px)`, so the L below lands on
    // `(0.25,0.75) (0.25,−0.25) (0.75,−0.25) (0.75,0.25) (1.0,0.25) (1.0,0.75)`:
    //   1. **non-convex** — the L has a reflex corner at `(0.75, 0.25)`;
    //   2. **overhanging** — `y < 0` cantilevers past the cube's `y = 0` edge;
    //   3. **flush** — the edge `x = 1.0, y∈[0.25,0.75]` lies *exactly* on the face's `x = 1`
    //      boundary, the "profile rim shares the face rim" case.
    //
    // Hand-checked shape: footprint `0.25 + 0.375 = 0.625`, prism wholly above `z = 1`, so
    //   volume 1 + 0.625 = 1.625
    //   area   5 (cube minus its top) + 0.5 (top left uncovered) + 3.5 (prism sides)
    //          + 0.625 (prism cap) + 0.125 (the cantilever's underside) = 9.75
    // The underside term is the cantilever: a contained pad would not have one.
    //
    // OCCT cannot score this directly — `pad` builds its tool prism internally, and rebuilding
    // it here would lean on the same frame mapping the assertion is testing.
    let (mut m, top) = cube_with_top();
    // On the lid, frame coordinates are world x and y: the origin is the world origin projected
    // onto `z = 1` and the axes are `u = +x̂`, `v = +ŷ`.
    let l_over = Profile2d::polygon(vec![
        p2(0.25, 0.75),
        p2(0.25, -0.25),
        p2(0.75, -0.25),
        p2(0.75, 0.25),
        p2(1.0, 0.25),
        p2(1.0, 0.75),
    ])
    .unwrap();
    let OpOutput::PadOnFace { solid, top_face } =
        apply(&mut m, &pad_op(top, l_over, 1.0)).expect("the cantilevered L pad")
    else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let p = nacre_props::mass_props(&m, solid).unwrap();
    assert!((p.volume - 1.625).abs() < 1e-12, "volume {}", p.volume);
    assert!((p.area - 9.75).abs() < 1e-12, "area {}", p.area);
    assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
}

#[test]
fn a_corner_flush_common_keeps_the_non_convex_overlap() {
    // A **corner-flush** `Common`: the L-prism and the box both start at the origin, so
    // **three** of their face planes coincide — `z = 0` (both floors), `x = 0`, `y = 0`. Every
    // vertex of the shared corner lies exactly on the other solid's face planes. The arrangement
    // names such a point by its plane triple like any other, so the configuration is not
    // special.
    //
    // Not covered by the other two non-convex `Common` locks:
    // `common_non_convex_overlap_is_their_intersection` (l_and_corner_box) and
    // `common_non_convex_containment_is_inner` both meet transversally, with no coplanar pair.
    //
    // Hand-checked shape, not just volume. The overlap is the L
    // `x∈[0,1.5]×y∈[0,1]` (1.5) plus `x∈[0,1]×y∈[1,1.5]` (0.5) = 2.0, over `z∈[0,0.5]`:
    //   volume 2.0 · 0.5 = 1.0
    //   area   2 · 2.0 (caps) + 6.0 (the L's perimeter) · 0.5 = 7.0
    // and the L has six sides, so eight faces.
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let mut m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    let lsolid = *m.live_solids().first().unwrap();
    let b = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.5, 1.5, 0.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Common, lsolid, b).expect("corner-flush Common");
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let p = nacre_props::mass_props(&m, r).unwrap();
    assert!((p.volume - 1.0).abs() < 1e-12, "volume {}", p.volume);
    assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
    assert_eq!(m.solid(r).cavities.len(), 0);
}

/// Two boxes that touch only along a single edge share a plane (`x = 4`) whose cross-section is
/// two rectangles meeting at one point; the unbounded contour of that plane pinches through the
/// point, tracing a figure-8. Their `Common` is empty (a measure-zero intersection), and it must
/// be empty in **both** operand orders — the figure-8's winding is read at a lex-extreme corner,
/// not at the pinch, so the verdict does not depend on where the ring happens to start.
#[test]
fn edge_contact_common_is_empty_in_both_orders() {
    let build = || {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([4.0, 2.0, 1.0]),
            Point3::from_array([5.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 0.0, 2.0]),
            Point3::from_array([4.0, 2.0, 3.0]),
        );
        (m, a, b)
    };
    for (label, swap) in [("a,b", false), ("b,a", true)] {
        let (mut m, a, b) = build();
        let (x, y) = if swap { (b, a) } else { (a, b) };
        let out = boolean(&mut m, BoolKind::Common, x, y).unwrap_or_else(|e| {
            panic!("edge-contact Common ({label}) should be empty, not reject: {e:?}")
        });
        assert!(
            out.is_empty(),
            "edge-contact Common ({label}) should be empty, got {} solid(s)",
            out.len()
        );
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
    }
}

/// Every vertex of a solid, cavities and inner loops included, as sorted coordinate bits — the
/// identity [`a_box_resting_across_a_rim_leaves_what_it_cuts_whole`] holds a result to.
fn vertex_bits(m: &Model, s: Handle<Solid>) -> Vec<[u64; 3]> {
    let sol = m.solid(s);
    let mut seen = std::collections::HashSet::new();
    for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
        for &fh in &m.shell(sh).faces {
            let f = m.face(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    seen.extend(m.edge(he.edge).vertices);
                }
            }
        }
    }
    let mut bits: Vec<[u64; 3]> = seen
        .into_iter()
        .map(|v| m.vertex_point(v).as_array().map(f64::to_bits))
        .collect();
    bits.sort_unstable();
    bits
}

fn face_count(m: &Model, s: Handle<Solid>) -> usize {
    let sol = m.solid(s);
    std::iter::once(sol.outer)
        .chain(sol.cavities.iter().copied())
        .map(|sh| m.shell(sh).faces.len())
        .sum()
}

/// A cylinder of radius 2 from `z = 0` to `z = 1` about the `z` axis, through the product's own road.
fn round_pin(m: &mut Model) -> Handle<Solid> {
    nacre_ops::fixtures::cylinder(
        m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        1.0,
    )
    .solid
}

/// An 8 × 8 plate from `z = −1` to `z = 1` with a bore of radius 2 through it about the `z` axis.
fn bored_plate(m: &mut Model) -> Handle<Solid> {
    let plate = m.add_cuboid(
        Point3::from_array([-4.0, -4.0, -1.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    let drill = nacre_ops::fixtures::cylinder(
        m,
        Point3::from_array([0.0, 0.0, -2.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    boolean_one(m, BoolKind::Cut, plate, drill).expect("the plate is bored")
}

/// [`plate_with_a_fillet`] one unit thick.
fn fillet_plate(m: &mut Model) -> Handle<Solid> {
    let profile = plate_with_a_fillet();
    let frame = nacre_ops::SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("the plate extrudes") else {
        unreachable!()
    };
    solid
}

/// ★ **A box resting on a cap across its rim leaves what it cuts whole.** `B` touches `A`'s cap
/// face to face and nowhere else, and the contact's boundary crosses the cap's rim, so the cap
/// class splits the rim at `B`'s planes and the cut merges the cap back. `A − B` is `A` itself —
/// the same vertices to the bit, the same faces — `B − A` is `B`, the Fuse is one body of both
/// volumes and Common is empty; every result valid, in both operand orders.
///
/// What it holds, two rules. The lateral and the assembly meet on the rim as the cleaned faces
/// hold it (`draft::HeldRims`), not as the arrangement split it — the split cut the lateral's rim
/// at nodes the merged cap had dissolved (an edge used once, `OpenResultShell`; the fillet
/// cells). And a whole rim a tool crossed dissolves back into its circle: crossed twice, the cap
/// ring is two arcs whose nodes have one neighbour each (they stayed, naming the box's planes —
/// `VertexNamesAbsentSurface`); crossed four times it dissolved to no node at all (an empty loop).
/// The fillet slab runs flush with the plate's bottom wall on one row and stops short of it on
/// the other: the wall plane is not what the rim meets.
#[test]
fn a_box_resting_across_a_rim_leaves_what_it_cuts_whole() {
    type Build = fn(&mut Model) -> Handle<Solid>;
    let cases: &[(&str, Build, [f64; 3], [f64; 3])] = &[
        (
            "cylinder, a box inside the top cap (no crossing)",
            round_pin,
            [-0.5, -0.5, 1.0],
            [0.5, 0.5, 2.0],
        ),
        (
            "cylinder, a box across the top rim",
            round_pin,
            [1.5, -0.5, 1.0],
            [2.5, 0.5, 2.0],
        ),
        (
            "cylinder, a box across the base rim",
            round_pin,
            [-0.5, 1.5, -1.0],
            [0.5, 2.5, 0.0],
        ),
        (
            "cylinder, a bar across the top rim twice",
            round_pin,
            [-3.0, -0.5, 1.0],
            [3.0, 0.5, 2.0],
        ),
        (
            "bored plate, a box across the bore's rim",
            bored_plate,
            [1.5, -0.5, 1.0],
            [2.5, 0.5, 2.0],
        ),
        (
            "bored plate, a bar across the bore's rim twice",
            bored_plate,
            [-3.0, -0.5, 1.0],
            [3.0, 0.5, 2.0],
        ),
        (
            "fillet plate, a slab flush with the wall plane",
            fillet_plate,
            [1.6, -4.0, 1.0],
            [2.6, 0.0, 2.0],
        ),
        (
            "fillet plate, a slab short of the wall plane",
            fillet_plate,
            [1.6, -3.9, 1.0],
            [2.6, -2.9, 2.0],
        ),
    ];
    for &(name, build, lo, hi) in cases {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            for swapped in [false, true] {
                let at = format!("{name}: {kind:?}, swapped {swapped}");
                let mut m = Model::new();
                let a = build(&mut m);
                let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
                m.rebuild_adjacency();
                let (va, vb) = (volume(&m, a), volume(&m, b));
                let (bits_a, faces_a) = (vertex_bits(&m, a), face_count(&m, a));
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let out = boolean(&mut m, kind, x, y).unwrap_or_else(|e| panic!("{at}: {e:?}"));
                m.rebuild_adjacency();
                assert!(nacre_validate::validate(&m).is_empty(), "{at}: validate");
                match (kind, swapped) {
                    (BoolKind::Common, _) => assert!(out.is_empty(), "{at}: {out:?}"),
                    (BoolKind::Fuse, _) => {
                        assert_eq!(out.len(), 1, "{at}: one body");
                        let v = volume(&m, out[0]);
                        assert!((v - (va + vb)).abs() < 1e-9, "{at}: volume {v}");
                    }
                    (BoolKind::Cut, true) => {
                        assert_eq!(out.len(), 1, "{at}: one body");
                        let v = volume(&m, out[0]);
                        assert!((v - vb).abs() < 1e-9, "{at}: volume {v}");
                    }
                    (BoolKind::Cut, false) => {
                        assert_eq!(out.len(), 1, "{at}: one body");
                        assert_eq!(vertex_bits(&m, out[0]), bits_a, "{at}: A's vertices");
                        assert_eq!(face_count(&m, out[0]), faces_a, "{at}: A's faces");
                        let v = volume(&m, out[0]);
                        assert!((v - va).abs() < 1e-9, "{at}: volume {v}");
                    }
                }
            }
        }
    }
}
