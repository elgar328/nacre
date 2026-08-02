//! Coplanar contact — flush faces, shared caps, hollow parts, coplanar seams.

#![allow(unused_imports)]
use crate::common::*;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean,
    replay,
};
use nacre_scalar::Axis;
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
    // the [1,2]³ void, so the void is carried through and preserved (cell (5c-in)).
    // The seam front-end walks all shells, so the void no longer silently vanishes
    // (it used to read as convex and return vol 26.875, cavities 0). Now: correct
    // 25.875 (27 − 1 void − 0.125 corner) with the cavity intact.
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);
    m.rebuild_adjacency();
    let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, hollow, cutter).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
    assert_eq!(m.solids.get(r).cavities.len(), 1);
}

#[test]
fn a_slab_splits_a_hollow_box_into_two() {
    // A slab cut through the whole box (and its void) severs it into two solids (cell 0.4).
    // The slab spans the full cross-section, so it opens the void — both pieces are
    // cavity-free. (Before cell 0.4 this was rejected as `DISCONNECTED_RESULT`.)
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
        assert_eq!(m.solids.get(s).cavities.len(), 0);
    }
    // Hollow 26 (= 27 − 1 void); the slab removes 1.6 of material (8 area × 0.2 thick).
    let vol: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((vol - 24.4).abs() < 1e-9, "total volume {vol}");
}

// A hollow operand in a COPLANAR contact — the combination nothing covered until now. The
// cavity goldens above all take the seam path (transversal cuts) and every coplanar golden uses
// solid operands, so the intersection of the two was a blind spot, and the coplanar driver had
// never received the all-shell patch the seam front-end got in cell (5c-in). It emitted only
// outer-shell faces, so the void vanished: the fuse read 27.0625 — the *un-hollowed* cube plus
// the boss — with `cavities: 0` and a clean `validate`, because what remained was still a
// closed shell. Silent-wrong, invisible to every guard. Now the driver walks all shells.
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
    assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
}

#[test]
fn a_hollow_part_takes_a_coplanar_pocket() {
    // The Cut twin of the boss case: a top-flush pocket sunk into a hollow part. Same blind
    // spot, same silent-wrong before the fix (26.96875 with the void gone).
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
    assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
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
    assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
}

#[test]
fn a_hollow_part_takes_a_second_far_cut() {
    // The headline: keep cutting a part after it is hollow. A bore at the corner
    // opposite the void — the void survives, and the result is fed back as an
    // operand (chaining past the first cavity-producing op, which the door used to
    // block).
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
    assert_eq!(m.solids.get(r).cavities.len(), 1);
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
    // It used to reject because the overhang sidecars gated on convexity; that gate is gone.
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
    // In the lid's own frame, whose origin is the world origin projected onto `z = 1` and whose
    // axes are `u = −ŷ`, `v = +x̂`: a frame point `(a, b)` is world `(b, −a, 1)`.
    let l_over = Profile2d::polygon(vec![
        p2(-0.75, 0.25),
        p2(0.25, 0.25),
        p2(0.25, 0.75),
        p2(-0.25, 0.75),
        p2(-0.25, 1.0),
        p2(-0.75, 1.0),
    ]);
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
    // vertex of the shared corner lies exactly on the other solid's face planes, which is what
    // the old reject tag said: `VERTEX_ON_FACE_PLANE`. The arrangement engine names such a
    // point by its plane triple like any other, so the configuration is no longer special.
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
    ]);
    let mut m = replay(&[extrude_op(l, 1.0)]).unwrap();
    let lsolid = *m.live_solids.first().unwrap();
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
    assert_eq!(m.solids.get(r).cavities.len(), 0);
}

/// Two boxes that touch only along a single edge share a plane (`x = 4`) whose cross-section is
/// two rectangles meeting at one point; the unbounded contour of that plane pinches through the
/// point, tracing a figure-8. Their `Common` is empty (a measure-zero intersection), and it must
/// be empty in **both** operand orders — the figure-8's winding is read at a lex-extreme corner,
/// not at the pinch, so the verdict no longer depends on where the ring happens to start.
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
