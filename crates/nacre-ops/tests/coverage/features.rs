//! Feature operations — pad (boss) and pocket, each a tool prism + boolean sugar
//! over `apply(Operation::Pad/PocketOnFace)`. Covers the happy path (topology of
//! a boss/pocket on a cube top) and the honest rejections (non-planar face,
//! non-positive distance, degenerate profile, through-pocket).

use crate::common::{cube_with_top, p2, pad_op, pocket_op, small_square};
use nacre_geom::Surface;
use nacre_math::{Point3, Vector3};
use nacre_ops::{OpError, OpOutput, Profile2d, apply};
use nacre_topo::Model;

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
    let inner: usize = reach.faces.iter().map(|fh| m.faces.get(*fh).inner.len()).sum();
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
    let inner: usize = reach.faces.iter().map(|fh| m.faces.get(*fh).inner.len()).sum();
    assert_eq!(inner, 1);
    assert!(reach.faces.contains(&bottom_face));
}

/// The lateral face of a cylinder — a non-planar target both pad and pocket reject.
fn cylinder_lateral(m: &mut Model) -> nacre_store::Handle<nacre_topo::Face> {
    m.add_cylinder(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        5.0,
    );
    let shell = m.solids.get(m.live_solids[0]).outer;
    *m.shells
        .get(shell)
        .faces
        .iter()
        .find(|&&fh| matches!(m.surfaces.get(m.faces.get(fh).surface), Surface::Cylinder(_)))
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
    let two = Profile2d {
        points: vec![p2(0.0, 0.0), p2(0.1, 0.0)],
    };
    assert!(matches!(
        apply(&mut m, &pad_op(top, two, 0.5)),
        Err(OpError::DegenerateProfile)
    ));
}

#[test]
fn pocket_rejects_degenerate_profile() {
    let (mut m, top) = cube_with_top();
    let two = Profile2d {
        points: vec![p2(0.0, 0.0), p2(0.1, 0.0)],
    };
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
