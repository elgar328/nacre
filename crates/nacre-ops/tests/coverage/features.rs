//! Feature operations — pad (boss) and pocket + honest rejections.

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

#[test]
fn pad_boss_on_cube_top() {
    let (mut m, top) = cube_with_top();
    let out = apply(&mut m, &pad_op(top, small_square(), 0.5)).unwrap();
    let OpOutput::PadOnFace { top_face, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let v = nacre_validate::validate(&m);
    assert!(v.is_empty(), "{v:?}");
    let reach = m.reachable();
    // 6 cube faces − top + (outer' + 4 walls + cap) = 11.
    assert_eq!(reach.faces.len(), 11);
    assert_eq!(reach.vertices.len(), 16); // 8 cube + 4 base + 4 top
    assert_eq!(reach.edges.len(), 24); // 12 cube + 4 base + 4 top + 4 vertical
    let inner: usize = reach.faces.iter().map(|fh| m.face(*fh).inner.len()).sum();
    assert_eq!(inner, 1);
    assert!(reach.faces.contains(&top_face));
}

#[test]
fn pocket_on_cube_top() {
    let (mut m, top) = cube_with_top();
    let out = apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap();
    let OpOutput::PocketOnFace { bottom_face, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let v = nacre_validate::validate(&m);
    assert!(v.is_empty(), "{v:?}");
    // Same topology as a boss (a downward prism instead of upward).
    let reach = m.reachable();
    assert_eq!(reach.faces.len(), 11);
    assert_eq!(reach.vertices.len(), 16);
    assert_eq!(reach.edges.len(), 24);
    let inner: usize = reach.faces.iter().map(|fh| m.face(*fh).inner.len()).sum();
    assert_eq!(inner, 1);
    assert!(reach.faces.contains(&bottom_face));
}

/// The lateral face of a cylinder — a non-planar target both pad and pocket reject.
fn cylinder_lateral(m: &mut Model) -> nacre_store::Handle<nacre_topo::Face> {
    nacre_ops::fixtures::cylinder_with_seam(
        m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    );
    let shell = m.solid(m.live_solids()[0]).outer;
    *m.shell(shell)
        .faces
        .iter()
        .find(|&&fh| matches!(m.surface_cache(m.face(fh).surface), Surface::Cylinder(_)))
        .unwrap()
}

#[test]
fn pad_rejects_nonplanar_face() {
    let mut m = Model::new();
    let lateral = cylinder_lateral(&mut m);
    assert!(matches!(
        apply(&mut m, &pad_op(lateral, small_square(), 0.5)),
        Err(OpError::NonPlanarFace)
    ));
}

#[test]
fn pocket_rejects_nonplanar_face() {
    let mut m = Model::new();
    let lateral = cylinder_lateral(&mut m);
    assert!(matches!(
        apply(&mut m, &pocket_op(lateral, small_square(), 0.5)),
        Err(OpError::NonPlanarFace)
    ));
}

#[test]
fn pad_rejects_nonpositive_dist() {
    let (mut m, top) = cube_with_top();
    assert!(matches!(
        apply(&mut m, &pad_op(top, small_square(), 0.0)),
        Err(OpError::NonPositiveDistance)
    ));
}

#[test]
fn pocket_rejects_nonpositive_dist() {
    let (mut m, top) = cube_with_top();
    assert!(matches!(
        apply(&mut m, &pocket_op(top, small_square(), 0.0)),
        Err(OpError::NonPositiveDistance)
    ));
}

#[test]
fn pad_rejects_degenerate_profile() {
    let (mut m, top) = cube_with_top();
    let two = Profile2d::polygon(vec![p2(0.0, 0.0), p2(0.1, 0.0)]).unwrap();
    assert!(matches!(
        apply(&mut m, &pad_op(top, two, 0.5)),
        Err(OpError::DegenerateProfile)
    ));
}

#[test]
fn pocket_rejects_degenerate_profile() {
    let (mut m, top) = cube_with_top();
    let two = Profile2d::polygon(vec![p2(0.0, 0.0), p2(0.1, 0.0)]).unwrap();
    assert!(matches!(
        apply(&mut m, &pocket_op(top, two, 0.5)),
        Err(OpError::DegenerateProfile)
    ));
}

#[test]
fn pocket_through_the_solid_is_rejected() {
    let (mut m, top) = cube_with_top(); // 1.0-thick cube
    let got = apply(&mut m, &pocket_op(top, small_square(), 1.5));
    assert!(
        matches!(got, Err(OpError::Boolean(_)) | Err(OpError::PocketNotBlind)),
        "through-pocket must reject honestly, got {got:?}"
    );
}

#[test]
fn cut_by_a_box_inside_the_pocket_is_a_no_op() {
    // The box sits wholly in the void, so the solids are disjoint and `A − B = A`.
    // Reading the lid as filled instead classified the box's eight corners five
    // Inside and three Outside, and the seam-free path's own `debug_assert`
    // ("classification must be consistent per solid") caught it.
    let (mut m, pc) = pocketed_cube();
    let bx = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.4, 0.4, 0.6]),
        Point3::from_array([0.6, 0.6, 0.9]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    assert!(m.solid(r).cavities.is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.92).abs() < 1e-9, "volume {vol}");
}

#[test]
fn pad_an_overhanging_boss() {
    // The profile reaches past one face edge: part of the boss sits on the face, part
    // cantilevers into the air. Relaxing the containment gate routes it to the overhang Fuse
    // sidecar (Ok here proves the routing — a contained-only pad would reject). The boss lives
    // wholly above z=1, so vol = cube 1 + footprint 0.5 · dist 1 = 1.5.
    let (mut m, top) = cube_with_top();
    let OpOutput::PadOnFace { solid, top_face } =
        apply(&mut m, &pad_op(top, edge_overhang_profile(), 1.0)).unwrap()
    else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.5).abs() < 1e-12);
    assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
}

#[test]
fn pad_a_spanning_slab_boss() {
    // A slab crossing the whole face (overhangs two opposite edges). vol = 1 + 0.75 · 1 = 1.75.
    let (mut m, top) = cube_with_top();
    let out = apply(&mut m, &pad_op(top, spanning_slab_profile(), 1.0)).unwrap();
    let OpOutput::PadOnFace { solid, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.75).abs() < 1e-12);
}

#[test]
fn pocket_an_edge_slot() {
    // A blind pocket whose footprint overhangs one edge — an edge slot open to the side.
    // Only the on-face part (world x[0.25,0.75]×y[0,0.75] = 0.375) carves: 1 − 0.375·0.5 = 0.8125.
    let (mut m, top) = cube_with_top();
    let OpOutput::PocketOnFace { solid, bottom_face } =
        apply(&mut m, &pocket_op(top, edge_overhang_profile(), 0.5)).unwrap()
    else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.8125).abs() < 1e-12);
    assert!(m.reachable().faces.contains(&bottom_face)); // slot floor recovered
}

#[test]
fn pocket_a_slab_channel() {
    // A blind channel crossing the whole face (breaches two opposite walls). On-face carve
    // world x[0.25,0.75]×y[0,1] = 0.5: 1 − 0.5·0.5 = 0.75.
    let (mut m, top) = cube_with_top();
    let out = apply(&mut m, &pocket_op(top, spanning_slab_profile(), 0.5)).unwrap();
    let OpOutput::PocketOnFace { solid, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.75).abs() < 1e-12);
}

#[test]
fn pad_overhang_off_the_face_is_rejected() {
    // A footprint that does not touch the face at all. The boolean is not what fails here — it
    // fuses the two into a base plus a detached boss, which is the right answer (see
    // `a_touchless_boss_fuses_into_two_solids`). What breaks is the *pad's* premise, so the
    // error names that, and the model the caller is left holding is the one it started with.
    let (mut m, top) = cube_with_top();
    let far =
        Profile2d::polygon(vec![p2(1.8, 1.8), p2(2.2, 1.8), p2(2.2, 2.2), p2(1.8, 2.2)]).unwrap();
    let before = m.live_solids().to_vec();
    assert_eq!(
        apply(&mut m, &pad_op(top, far, 0.3)),
        Err(OpError::PadMissesFace)
    );
    let (mut a, mut b) = (before, m.live_solids().to_vec());
    a.sort_by_key(|h| h.index());
    b.sort_by_key(|h| h.index());
    assert_eq!(a, b, "a rejected pad must leave the live model untouched");
}

/// The kernel's answer for a boss that misses the face, stated on its own so nobody "fixes" the
/// boolean to reject it: fusing two solids that do not touch **is** two solids, and both are
/// whole. Only `pad` refuses that outcome, because a pad is defined as material joined to a face
/// (`pad_overhang_off_the_face_is_rejected`).
#[test]
fn a_touchless_boss_fuses_into_two_solids() {
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    // Shares the z = 1 plane class with the base's top, but sits far away in x/y — so the plane
    // carries two separate bodies.
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.8, 1.8, 1.0]),
        Point3::from_array([2.2, 2.2, 1.3]),
    );
    let solids = boolean(&mut m, BoolKind::Fuse, base, boss).unwrap();
    assert_eq!(solids.len(), 2, "disjoint operands stay two solids");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let mut vols: Vec<f64> = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .collect();
    vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!(
        (vols[0] - 0.048).abs() < 1e-12 && (vols[1] - 1.0).abs() < 1e-12,
        "both pieces whole: {vols:?}"
    );
}

#[test]
fn pocket_through_overhang_is_rejected() {
    // An overhang pocket deep enough to pierce the far side is not blind — no single floor.
    // Honest reject via whichever path fires (the overhang detector declines, the seam path
    // rejects the mixed contact), mirroring `pocket_through_the_solid_is_rejected`.
    let (mut m, top) = cube_with_top();
    let got = apply(&mut m, &pocket_op(top, edge_overhang_profile(), 1.5));
    assert!(
        matches!(got, Err(OpError::Boolean(_)) | Err(OpError::PocketNotBlind)),
        "through overhang must reject honestly, got {got:?}"
    );
}

#[test]
fn pad_step_exports() {
    // The boss (holed outer face + walls + cap) exports without error.
    let (mut m, top) = cube_with_top();
    apply(&mut m, &pad_op(top, small_square(), 0.5)).unwrap();
    let step = nacre_step::to_step(&m).expect("boss exports");
    assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
}

/// A prism raised on a `(1,1,1)`-slanted sketch plane, its far cap holding a blind pocket. The
/// cap and its two anti-parallel side walls meet in a triple whose `raw` coefficients are
/// exactly dependent (`det = 0`); before family #3's dir-sign fix the guard read `sqrt`-rounded
/// unit normals, called that triple non-degenerate, and the consumer aborted on `D = 0`. Now
/// the guard reads the same coefficients the consumer does, so the arrangement runs. Volume:
/// a `2×2` base × `2` deep block is `8`, less the `0.4²×0.5` pocket.
#[test]
fn a_pocket_on_a_slanted_face() {
    let plane =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    let mut m = Model::new();
    let big = Profile2d::polygon(vec![
        p2(-1.0, -1.0),
        p2(1.0, -1.0),
        p2(1.0, 1.0),
        p2(-1.0, 1.0),
    ])
    .unwrap();
    let __f101 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f101,
            profile: big,
            dist: 2.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let out = apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap();
    let OpOutput::PocketOnFace { solid, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((vol - (8.0 - 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
}

/// The same slanted cap, but a boss (pad, Fuse) instead of a pocket — the sweep runs the other
/// way, a different code path. Volume: the `8` block plus a `0.4²×0.5` stub.
#[test]
fn a_pad_on_a_slanted_face() {
    let plane =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    let mut m = Model::new();
    let big = Profile2d::polygon(vec![
        p2(-1.0, -1.0),
        p2(1.0, -1.0),
        p2(1.0, 1.0),
        p2(-1.0, 1.0),
    ])
    .unwrap();
    let __f100 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f100,
            profile: big,
            dist: 2.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let out = apply(
        &mut m,
        &Operation::PadOnFace {
            face: faces[1],
            profile: small_square(),
            dist: 0.5,
        },
    )
    .unwrap();
    let OpOutput::PadOnFace { solid, .. } = out else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((vol - (8.0 + 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn pocket_step_exports() {
    let (mut m, top) = cube_with_top();
    apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap();
    let step = nacre_step::to_step(&m).expect("pocket exports");
    assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
}

/// `pocket_corner_cut` by hand, so the pocket family keeps a regression net that runs without
/// OCCT: the unit cube less a `0.4²×0.5` pocket is `0.92`, and the corner box `[0.85,1.15]³`
/// bites `0.15³` of solid (it clears the pocket, whose footprint stops at `x = 0.7`).
#[test]
fn a_corner_cut_off_a_pocketed_cube() {
    let (mut m, pc) = pocketed_cube();
    let bx = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.85, 0.85, 0.85]),
        Point3::from_array([1.15, 1.15, 1.15]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let props = nacre_props::mass_props(&m, r).unwrap();
    assert!(
        (props.volume - (0.92 - 0.15 * 0.15 * 0.15)).abs() < 1e-12,
        "volume {}",
        props.volume
    );
    // A corner bite replaces three 0.15² squares with three more: the area is unchanged at
    // 6 − 0.16 (the lid's hole) + 0.8 (four pocket walls) + 0.16 (its floor).
    assert!((props.area - 6.8).abs() < 1e-12, "area {}", props.area);
}

#[test]
fn an_edge_slot_through_the_bottom() {
    // The prism pokes out the base's bottom too, and the boolean carves the slot exactly:
    // base 1.0 − (x∈[0.5,1] · y∈[0.25,0.75] · z∈[0,1]) = 1 − 0.25 = 0.75.
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let prism = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.25, -0.5]),
        Point3::from_array([1.5, 0.75, 1.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
}

#[test]
fn a_boss_that_pierces_the_base_is_not_an_overhang() {
    // The boss dips below the base's top (its walls cross the base) — a transversal seam cut,
    // not a coplanar overhang, which the arrangement handles as a normal crossing.
    // Union = 1.0 + boss 0.75 − overlap 0.125 = 1.625.
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let through = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.25, 0.5]),
        Point3::from_array([1.5, 0.75, 2.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.625).abs() < 1e-12, "volume {vol}");
}

// The Cut and Common twins of the corner-overhanging boss fuse below — the same two solids,
// the same shared z=1 plane. The boss sits entirely above it, so it removes nothing and shares
// nothing: Cut is the base untouched and Common is empty.
#[test]
fn cut_by_a_corner_overhanging_boss_removes_nothing() {
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let corner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
    // Structure, not just volume: the base comes through as itself. Six faces means the cap was
    // not split along ∂Q and the boss contributed nothing.
    let s = m.solid(r);
    assert!(s.cavities.is_empty());
    assert_eq!(m.shell(s.outer).faces.len(), 6, "a clean cube");
}

/// A coplanar-contact undercut: the tool's boss annulus sits flush on the base's top (z=3, a
/// contact that removes nothing) while its pin reaches back below that plane into the base and
/// carves a notch. The arrangement cuts the notch (27 − 1 = 26) and keeps the flush contact a
/// no-op.
#[test]
fn a_coplanar_boss_with_a_pin_undercuts_the_base() {
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 3.0]),
        Point3::from_array([2.5, 2.5, 4.0]),
    );
    let pin = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 2.0]),
        Point3::from_array([2.0, 2.0, 3.0]),
    );
    let tool = boolean_one(&mut m, BoolKind::Fuse, boss, pin).unwrap();
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, base, tool).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 26.0).abs() < 1e-12, "volume {vol}");
    assert!(
        m.solid(r).cavities.is_empty(),
        "the notch is a surface indentation, not a void"
    );
}

#[test]
fn common_with_a_corner_overhanging_boss_is_empty() {
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let corner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    m.rebuild_adjacency();
    // They meet only along the base's top face — a contact of zero volume.
    assert!(
        boolean(&mut m, BoolKind::Common, base, corner)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn fuse_an_overhanging_boss_onto_a_non_convex_solid() {
    // A boss cantilevers off the +x side face of a top-pocketed cube (non-convex solid),
    // overhanging the bottom edge. Volume: pocketed 0.92 + boss 0.25 = 1.17.
    let (mut m, pc) = top_pocketed_cube();
    let boss = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 0.25, -0.25]),
        Point3::from_array([1.5, 0.75, 0.75]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, pc, boss).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    assert!((nacre_props::mass_props(&m, r).unwrap().volume - 1.17).abs() < 1e-12);
}

#[test]
fn cut_a_blind_pocket_into_a_non_convex_solid() {
    // A blind pocket carved into an already-pocketed (non-convex) cube: a second contained
    // top-flush prism in a corner away from the first pocket. The kept solid `a` is non-convex,
    // which the pocket contact now admits (the gates are convexity-agnostic). Removed
    // 0.2·0.1·0.4 = 0.008 on top of the first pocket's 0.08 → 1 − 0.08 − 0.008 = 0.912.
    let (mut m, pc) = pocketed_cube();
    let corner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.05, 0.1, 0.6]),
        Point3::from_array([0.25, 0.2, 1.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, pc, corner).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.912).abs() < 1e-12, "volume {vol}");
}

#[test]
fn cut_a_non_convex_blind_pocket() {
    // A blind pocket with a non-convex (L-shaped) footprint: the cutter prism is non-convex,
    // which the pocket contact now admits. The L extrudes to z∈[0,0.5], top-flush on the base's
    // z=0.5 face, blind. L area = 0.6² − 0.3² = 0.27, depth 0.5 → removed 0.135; base 3²·1.5 =
    // 13.5 → 13.365.
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, -1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 0.5]),
    );
    let l = Profile2d::polygon(vec![
        p2(-0.3, -0.3),
        p2(0.3, -0.3),
        p2(0.3, 0.0),
        p2(0.0, 0.0),
        p2(0.0, 0.3),
        p2(-0.3, 0.3),
    ])
    .unwrap();
    let __op = extrude_op(&m, l, 0.5);
    let OpOutput::Extrude { solid: lp, .. } = apply(&mut m, &__op).unwrap() else {
        unreachable!()
    };
    let r = boolean_one(&mut m, BoolKind::Cut, base, lp).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 13.365).abs() < 1e-12, "volume {vol}");
}

#[test]
fn a_pocket_that_punches_through_drills_a_bore() {
    // A top-flush tool that pokes out the base's bottom: the pocket becomes a through hole.
    // The exit face has to come out annular, and until the coplanar reconstruct learned to
    // emit a hole it came out whole instead, leaving the bore's walls nothing to close
    // against — an open shell the assembly guard rejected. (Honest reject, never a wrong
    // answer; the previous cell pinned it as such.)
    //
    // Area is the assertion that matters here: volume alone cannot tell a bore from a shape
    // that merely displaces the same material. 0.75 (top) + 0.75 (bottom) + 4 (sides) +
    // 2.0 (the bore's four inner walls) = 7.5, against 6.0 for the cube.
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let through = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.25, 0.25, -0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let p = nacre_props::mass_props(&m, r).unwrap();
    assert!((p.volume - 0.75).abs() < 1e-12, "volume {}", p.volume);
    assert!((p.area - 7.5).abs() < 1e-12, "area {}", p.area);
    // A bore, not a void: no cavity shell, and both caps carry the hole (the top from the
    // coincident contact, the bottom from the section the tool cuts through it).
    let s = m.solid(r);
    assert!(s.cavities.is_empty(), "a through hole is not a cavity");
    let faces = &m.shell(s.outer).faces;
    assert_eq!(faces.len(), 10, "6 base faces + the bore's 4 walls");
    assert_eq!(
        faces
            .iter()
            .filter(|&&fh| !m.face(fh).inner.is_empty())
            .count(),
        2,
        "both caps are annular"
    );
}

#[test]
fn a_boss_that_punches_through_keeps_the_stub() {
    // The Fuse twin, and the same emission path: the tool's a-side face keeps `P∖Q`, so the
    // base's bottom needs the same hole for the stub below it to join on. Volume
    // 1 + 0.5·0.5·0.5 = 1.125; area 1.0 (top, the flush tool cap dissolves into it) + 0.75
    // (bottom) + 4 (sides) + 1.0 (stub walls) + 0.25 (stub floor) = 7.0.
    let mut m = Model::new();
    let base = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let through = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.25, 0.25, -0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let p = nacre_props::mass_props(&m, r).unwrap();
    assert!((p.volume - 1.125).abs() < 1e-12, "volume {}", p.volume);
    assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
    let s = m.solid(r);
    assert_eq!(
        m.shell(s.outer)
            .faces
            .iter()
            .filter(|&&fh| !m.face(fh).inner.is_empty())
            .count(),
        1,
        "only the bottom is annular — the flush top merges away"
    );
}

// --- face_plane: the frame pad/pocket actually use ---

/// **The contract that would rot silently.** An app asks `face_plane` where the
/// sketch origin is, then calls `pad` — and the kernel recomputes the frame. If the
/// two ever diverge, the boss lands somewhere the app did not predict and nothing
/// complains. So: place an **off-centre** profile through `face_plane` by hand and
/// require `pad` to put the boss exactly there. Centred profiles would pass even if
/// the frames disagreed only in origin, which is the divergence most likely to happen.
#[test]
fn face_plane_is_the_frame_pad_places_profiles_in() {
    let (mut m, top) = cube_with_top();
    let plane = nacre_ops::face_plane(&m, top).unwrap();

    // A square from (0.7, 0.2) to (0.9, 0.4) in face coordinates — asymmetric in both axes, so any
    // origin or axis mismatch moves its centre. On a lid these are world x and y, which is what
    // puts the boss inside the cube.
    let profile =
        Profile2d::polygon(vec![p2(0.7, 0.4), p2(0.7, 0.2), p2(0.9, 0.2), p2(0.9, 0.4)]).unwrap();
    let dist = 0.5;
    let OpOutput::PadOnFace { top_face, .. } = apply(&mut m, &pad_op(top, profile, dist)).unwrap()
    else {
        unreachable!()
    };
    m.rebuild_adjacency();

    // Where the caller predicts the boss's cap centre is, from `face_plane` alone.
    let n = plane.x_axis().cross(plane.y_axis());
    let want = plane.origin() + plane.x_axis() * 0.8 + plane.y_axis() * 0.3 + n * dist;

    let got = nacre_props::face_props(&m, top_face).unwrap().centroid;
    assert!(
        (got - want).norm() < 1e-12,
        "pad placed the boss at {got:?}, face_plane predicted {want:?}"
    );
}

/// **The same contract, on a turned solid — the population the test above cannot reach.**
///
/// An axis-aligned face never goes near the flip/sketch-frame branch, so the pin above cannot see
/// a report that spells `flip` apart from the realization. On a turned face it can: a report that
/// half-turned about one axis while the realization (`frame_chain`) half-turned about the other
/// left the two point-symmetric in the plane, a footprint centred on the face through
/// `face_plane`'s own coordinates was built on the opposite side — outside the face — and a legal
/// pad came back `PadMissesFace` (measured: 2 of 6 faces of a 30°-turned block, every axis).
///
/// All six faces, deliberately: the two flip=true walls are where the two spellings part, the two
/// flip=false walls are the contrast, and the two rotation-invariant planes take the world-axis
/// branch — three populations under one assertion, so none can drift out of the contract
/// unnoticed.
#[test]
fn face_plane_is_the_frame_pad_places_profiles_in_on_a_turned_face() {
    use nacre_exact::{Angle, Isometry, Rat, Rotation};
    for axis in [Axis::X, Axis::Y, Axis::Z] {
        for fi in 0..6 {
            let mut m = Model::new();
            let __op = extrude_op(&m, square(), 1.0);
            let OpOutput::Extrude { solid, .. } = apply(&mut m, &__op).unwrap() else {
                unreachable!()
            };
            m.rebuild_adjacency();
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid,
                    isometry: Isometry::rotation(Rotation {
                        axis,
                        pivot: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
                    }),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            m.rebuild_adjacency();
            let face = m.shell(m.solid(solid).outer).faces[fi];

            // Off-centre about the face's own centroid (the frame origin is a property of the
            // plane, not of the face, so absolute sketch coordinates can lie off the face
            // entirely on a turned solid). Asymmetric in both axes, as in the pin above.
            let plane = nacre_ops::face_plane(&m, face).unwrap();
            let d = nacre_props::face_props(&m, face).unwrap().centroid - plane.origin();
            let (cu, cv) = (d.dot(plane.x_axis()), d.dot(plane.y_axis()));
            let profile = Profile2d::polygon(vec![
                p2(cu + 0.05, cv - 0.15),
                p2(cu + 0.25, cv - 0.15),
                p2(cu + 0.25, cv + 0.05),
                p2(cu + 0.05, cv + 0.05),
            ])
            .unwrap();
            let dist = 0.5;
            let OpOutput::PadOnFace { top_face, .. } = apply(&mut m, &pad_op(face, profile, dist))
                .unwrap_or_else(|e| {
                    panic!("{axis:?} face[{fi}]: the pad missed the frame it was promised: {e:?}")
                })
            else {
                unreachable!()
            };
            m.rebuild_adjacency();

            let n = plane.x_axis().cross(plane.y_axis());
            let want = plane.origin()
                + plane.x_axis() * (cu + 0.15)
                + plane.y_axis() * (cv - 0.05)
                + n * dist;
            let got = nacre_props::face_props(&m, top_face).unwrap().centroid;
            assert!(
                (got - want).norm() < 1e-9,
                "{axis:?} face[{fi}]: pad placed the boss at {got:?}, face_plane predicted {want:?}"
            );
        }
    }
}

/// A rotated face's sketch frame is deterministic and right-handed about its own
/// outward normal — pinned so a change to `any_perpendicular` cannot quietly
/// rotate every sketch on such a face.
///
/// This face is **not** the discontinuous case, though the first draft of this
/// test claimed it was: `any_perpendicular` picks the *least*-aligned world axis,
/// and for a normal of (0.707, 0.707, 0) that is z by a mile. The genuine tie is a
/// normal whose two **smallest** components are equal, which axis-aligned rotation
/// by a rational angle cannot even produce here — so it is pinned where the rule
/// lives, in `nacre_math`'s own tests, not through a face.
#[test]
fn a_rotated_face_has_a_pinned_sketch_axis() {
    let (mut m, _) = cube_with_top();
    let s = m.live_solids()[0];
    let s = xf(&mut m, s, rot_iso(Axis::Z, 45));
    m.rebuild_adjacency();

    // The side face whose outward normal points into +x+y.
    let diag = Vector3::from_array([1.0, 1.0, 0.0]).normalize().unwrap();
    let face = *m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .find(|&&f| {
            nacre_props::face_props(&m, f)
                .unwrap()
                .normal
                .is_some_and(|v| v.dot(diag) > 0.99)
        })
        .expect("a rotated cube has a face facing the diagonal");

    let plane = nacre_ops::face_plane(&m, face).unwrap();
    // |n| = (0.707, 0.707, 0) — z is the least-aligned axis, so x = ẑ × n̂.
    let want = Vector3::from_array([0.0, 0.0, 1.0]).cross(diag);
    assert!(
        (plane.x_axis() - want).norm() < 1e-12,
        "x axis {:?}, expected {want:?}",
        plane.x_axis()
    );
    // And the frame stays right-handed about the face's outward normal.
    assert!((plane.x_axis().cross(plane.y_axis()) - diag).norm() < 1e-12);
}

/// **Why the sketch origin does not come from the face at all.**
///
/// A mean of the outer loop's corners drags sideways when a vertex is added along a straight edge
/// — the same face would then seat a boss somewhere else — and the face region's area centroid,
/// which does not move under subdivision, still reads the face's `f64` vertices: `construct.rs`
/// lifts the frame origin with `Rat::from_decimal`, so a rounded cache becomes the truth, and
/// padding one footprint twice leaves faces of area `2.2e-16`. The origin is the **world origin
/// projected onto the plane**, which does not read the face at all.
///
/// That is strictly stronger, and this test says so: subdividing an edge cannot move it, and
/// neither can replacing the outline with a different shape on the same plane. The old rule fails
/// the second of those, and a corner mean fails both.
#[test]
fn the_sketch_origin_does_not_depend_on_the_outline_at_all() {
    // An L, the same L with two extra points along its bottom edge, and a wholly different
    // outline — all extruded to a cap on the plane `z = 1`.
    let l = vec![
        p2(0.0, 0.0),
        p2(4.0, 0.0),
        p2(4.0, 2.0),
        p2(2.0, 2.0),
        p2(2.0, 4.0),
        p2(0.0, 4.0),
    ];
    let subdivided = vec![
        p2(0.0, 0.0),
        p2(1.0, 0.0), // mid-run on a straight edge — shape unchanged
        p2(3.0, 0.0),
        p2(4.0, 0.0),
        p2(4.0, 2.0),
        p2(2.0, 2.0),
        p2(2.0, 4.0),
        p2(0.0, 4.0),
    ];
    let elsewhere = vec![p2(7.0, 7.0), p2(9.0, 7.0), p2(9.0, 9.5), p2(7.0, 9.5)];

    let top_origin = |pts: Vec<Point2>| {
        let mut m = Model::new();
        let __op = extrude_op(&m, Profile2d::polygon(pts).unwrap(), 1.0);
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &__op).unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        nacre_ops::face_plane(&m, faces[1]).unwrap().origin()
    };

    let a = top_origin(l.clone());
    for (other, what) in [
        (subdivided, "subdividing an edge"),
        (elsewhere, "a different outline"),
    ] {
        assert_eq!(
            a.as_array(),
            top_origin(other).as_array(),
            "{what} moved the sketch origin"
        );
    }
    // It is the world origin projected onto the cap's plane — exactly, not nearly.
    assert_eq!(a.as_array(), [0.0, 0.0, 1.0], "{a:?}");

    // An outline-reading rule fails above, either one: the area centroid of
    // the L is (5/3, 5/3) and of the far rectangle (8, 8.25); the corner mean sits at (2, 2) on
    // the plain L and moves to (1.75, 1.5) once subdivided.
    let mean = |pts: &[Point2]| {
        let (sx, sy) = pts
            .iter()
            .fold((0.0, 0.0), |(x, y), p| (x + p[0], y + p[1]));
        [sx / pts.len() as f64, sy / pts.len() as f64]
    };
    let (ma, mb) = (
        mean(&l),
        mean(&[
            p2(0.0, 0.0),
            p2(1.0, 0.0),
            p2(3.0, 0.0),
            p2(4.0, 0.0),
            p2(4.0, 2.0),
            p2(2.0, 2.0),
            p2(2.0, 4.0),
            p2(0.0, 4.0),
        ]),
    );
    let moved = ((ma[0] - mb[0]).powi(2) + (ma[1] - mb[1]).powi(2)).sqrt();
    assert!(
        moved > 0.2,
        "the corner-mean rule must visibly move: {ma:?} vs {mb:?}"
    );
}

// ---- the signed extrude builds what pad and pocket build ----

/// One face to work on, built afresh for each road — `Model` does not clone, and its
/// construction is deterministic, so two builds hold the same handles.
type FaceBuilder = fn() -> (Model, Handle<Face>);

/// The profile a case draws on its face, read in that face's frame.
type ProfileOn = fn(&Model, Handle<Face>) -> Profile2d;

/// A `2×2×2` prism on the `(1,1,1)`-slanted plane, and its far cap.
fn slanted_cap() -> (Model, Handle<Face>) {
    let plane =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    let mut m = Model::new();
    let big = Profile2d::polygon(vec![
        p2(-1.0, -1.0),
        p2(1.0, -1.0),
        p2(1.0, 1.0),
        p2(-1.0, 1.0),
    ])
    .unwrap();
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: big,
            dist: 2.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    (m, faces[1])
}

/// A `2`-cube turned 15° about `y`, optionally mirrored across `x = 0`, and its turned `+x` face —
/// a face frame carried by a chain with a reflection in it.
fn turned_face(mirror: bool) -> (Model, Handle<Face>) {
    use nacre_exact::{Angle, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let frame = nacre_ops::SketchFrame::world(&m, Axis::Z);
    let square = Profile2d::polygon(vec![
        p2(-1.0, -1.0),
        p2(1.0, -1.0),
        p2(1.0, 1.0),
        p2(-1.0, 1.0),
    ])
    .unwrap();
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square,
            dist: 2.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let down = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(-1)]);
    let solid = xf(&mut m, solid, down);
    let turn = Isometry::rotation(Rotation {
        axis: Axis::Y,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(15)).unwrap(),
    });
    let mut solid = xf(&mut m, solid, turn);
    if mirror {
        let OpOutput::Mirror { solid: mirrored } = apply(
            &mut m,
            &Operation::Mirror {
                solid,
                axis: Axis::X,
                offset: Rat::from_int(0),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        solid = mirrored;
    }
    let (s15, c15) = 15f64.to_radians().sin_cos();
    let n = [if mirror { -c15 } else { c15 }, 0.0, -s15];
    let face = *m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find(|&&fh| {
            let f = m.face(fh);
            let Surface::Plane(pl) = m.surface_cache(f.surface) else {
                return false;
            };
            let out = pl.normal() * f64::from(f.orientation.sign());
            (0..3).all(|k| (out.as_array()[k] - n[k]).abs() < 1e-9)
        })
        .expect("the turned +x face");
    (m, face)
}

/// The census `wf` base: a prism on a tilted decimal orthonormal frame, and the wall whose
/// canonical name has no rational default frame — the kernel's natural producer of wide names.
fn wf_wall() -> (Model, Handle<Face>) {
    let mut m = Model::new();
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(0.1111111111111111, 0.1234567890123456),
                p2(4.123456789012345, 0.2345678901234567),
                p2(3.9876543210987654, 3.1234567890123459),
                p2(0.2222222222222222, 2.765432109876543),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let wall = *m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            let s = m.face(f).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| nacre_exact::plane_frame_default(*c).is_none())
        })
        .expect("the wf population");
    (m, wall)
}

/// A `0.6` square centred on the face's outer-loop vertex average, in the face's own frame.
fn centred_square(m: &Model, face: Handle<Face>) -> Profile2d {
    let sp = nacre_ops::face_plane(m, face).expect("planar");
    let pts: Vec<Point3> = m
        .face(face)
        .outer
        .half_edges
        .iter()
        .map(|&he| m.vertex_point(m.he_start(he)))
        .collect();
    let n = pts.len() as f64;
    let c = Point3::from_array(std::array::from_fn(|k| {
        pts.iter().map(|p| p.as_array()[k]).sum::<f64>() / n
    }));
    let d = c - sp.origin();
    let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
    Profile2d::polygon(vec![
        p2(cu - 0.3, cv - 0.3),
        p2(cu + 0.3, cv - 0.3),
        p2(cu + 0.3, cv + 0.3),
        p2(cu - 0.3, cv + 0.3),
    ])
    .unwrap()
}

/// What a road left behind, bit for bit: the arena's size, then every live solid — volume, area
/// and centroid bits, its vertices' coordinate bits, and every face's surface handle, orientation
/// and surface cache (`Debug` of an `f64` is its shortest round-trip, so equal text is equal bits).
fn road_digest(m: &Model) -> Vec<String> {
    let mut out = vec![format!(
        "arena s{} v{} e{} f{}",
        m.surface_count(),
        m.vertex_count(),
        m.edge_count(),
        m.face_count()
    )];
    for &s in m.live_solids() {
        let p = nacre_props::mass_props(m, s).expect("props");
        let c = nacre_props::centroid(m, s).expect("a planar result has a centroid");
        let mut coords: Vec<[u64; 3]> = Vec::new();
        let mut faces = Vec::new();
        let solid = m.solid(s);
        for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let f = m.face(fh);
                faces.push(format!(
                    "{} {:?} {:?}",
                    f.surface.index(),
                    f.orientation,
                    m.surface_cache(f.surface)
                ));
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        coords.push(m.vertex_point(m.he_start(he)).as_array().map(f64::to_bits));
                    }
                }
            }
        }
        coords.sort_unstable();
        coords.dedup();
        out.push(format!(
            "solid {} {:016x}/{:016x}/{:016x},{:016x},{:016x} {coords:?}",
            s.index(),
            p.volume.to_bits(),
            p.area.to_bits(),
            c[0].to_bits(),
            c[1].to_bits(),
            c[2].to_bits(),
        ));
        out.extend(faces);
    }
    out
}

/// ★★ **Pad and pocket are a face's frame, a signed extrude, and a boolean — bit for bit.**
///
/// The face operations build their tool by the extrude's own road (`frame_rings` + `build_prism`,
/// the face's surface as the base cap) and sweep a pocket against `ŵ`; the public `Extrude` now
/// takes that sign. Each case runs both roads on the same face, built twice: the old operation,
/// and `face_sketch_frame` + `Extrude(±d)` + `Boolean` with the operands in the operation's order
/// (the face's solid first). What either road leaves behind — arena, live solids, coordinates,
/// surface caches — must agree. The cases are the shapes the face operations are tested on: a lid,
/// an overhang, a channel, a profile with a hole, the slanted cap, a turned face with and without
/// a reflection in its chain, and the wide-name wall.
#[test]
fn the_signed_extrude_builds_what_pad_and_pocket_build() {
    let lid_ring = |lo: f64, hi: f64| vec![p2(lo, lo), p2(hi, lo), p2(hi, hi), p2(lo, hi)];
    let holed = || Profile2d::with_holes(lid_ring(0.2, 0.8), vec![lid_ring(0.4, 0.6)]).unwrap();
    let cases: Vec<(&str, FaceBuilder, ProfileOn, f64, f64)> = vec![
        ("lid", cube_with_top, |_, _| small_square(), 0.5, 0.5),
        (
            "overhang",
            cube_with_top,
            |_, _| edge_overhang_profile(),
            1.0,
            0.5,
        ),
        (
            "channel",
            cube_with_top,
            |_, _| spanning_slab_profile(),
            1.0,
            0.5,
        ),
        ("slanted", slanted_cap, |_, _| small_square(), 0.5, 0.5),
        ("turned", || turned_face(false), centred_square, 1.0, 0.5),
        (
            "turned mirrored",
            || turned_face(true),
            centred_square,
            1.0,
            0.5,
        ),
        ("wf wall", wf_wall, centred_square, 0.4, 0.4),
    ];
    let run = |name: &str,
               build: &dyn Fn() -> (Model, Handle<Face>, Profile2d),
               pad: f64,
               pocket: f64| {
        for (kind, dist) in [(BoolKind::Fuse, pad), (BoolKind::Cut, pocket)] {
            let (mut old, face, profile) = build();
            let op = match kind {
                BoolKind::Fuse => pad_op(face, profile.clone(), dist),
                _ => pocket_op(face, profile.clone(), dist),
            };
            apply(&mut old, &op)
                .unwrap_or_else(|e| panic!("{name} {kind:?}: the operation: {e:?}"));
            let (mut new, face, _) = build();
            let solid = *new
                .live_solids()
                .iter()
                .find(|&&s| new.shell(new.solid(s).outer).faces.contains(&face))
                .expect("the face's solid");
            let frame = nacre_ops::face_sketch_frame(&new, face).expect("a planar live face");
            let signed = if matches!(kind, BoolKind::Fuse) {
                dist
            } else {
                -dist
            };
            let OpOutput::Extrude { solid: tool, .. } = apply(
                &mut new,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist: signed,
                },
            )
            .unwrap_or_else(|e| panic!("{name} {kind:?}: the extrude: {e:?}")) else {
                unreachable!()
            };
            apply(
                &mut new,
                &Operation::Boolean {
                    kind,
                    a: solid,
                    b: tool,
                },
            )
            .unwrap_or_else(|e| panic!("{name} {kind:?}: the boolean: {e:?}"));
            let (a, b) = (road_digest(&old), road_digest(&new));
            assert_eq!(a.len(), b.len(), "{name} {kind:?}: {a:#?}\nvs\n{b:#?}");
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(x, y, "{name} {kind:?}");
            }
        }
    };
    for (name, build, profile, pad, pocket) in cases {
        run(
            name,
            &|| {
                let (m, f) = build();
                let p = profile(&m, f);
                (m, f, p)
            },
            pad,
            pocket,
        );
    }
    run(
        "holed",
        &|| {
            let (m, f) = cube_with_top();
            (m, f, holed())
        },
        0.5,
        0.5,
    );
}
