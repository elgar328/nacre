//! Multi-solid output — disjoint operands and face-contact stacks.

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

#[test]
fn fuse_of_disjoint_boxes_is_two_solids() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    for &s in &solids {
        let v = nacre_props::mass_props(&m, s).unwrap().volume;
        assert!((v - 1.0).abs() < 1e-12, "each unit box survives whole: {v}");
    }
}

#[test]
fn cut_of_disjoint_is_a() {
    // A − B with B disjoint from A removes nothing ⇒ the result is A.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - vol_a).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.live_solids().to_vec(), vec![r]);
}

#[test]
fn common_stacked_cubes_is_empty() {
    // The stack shares only its interface plane, so the intersection has no volume.
    let (mut m, a, b) = stacked_cubes();
    assert!(boolean(&mut m, BoolKind::Common, a, b).unwrap().is_empty());
    assert!(m.live_solids().is_empty(), "both operands are consumed");
}

#[test]
fn cut_stacked_cubes_is_a() {
    let (mut m, a, b) = stacked_cubes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
}

/// The same two solids the other way round: the bar severs the rod, and `Cut` answers with
/// two solids (cell 0.4). Each severed stub is its own genus-0 box, so `validate` is clean —
/// the pieces share no vertices or edges. (Before cell 0.4 this was `DISCONNECTED_RESULT`:
/// one handle could not name two solids, and forcing both into one shell read as
/// `NegativeGenus { genus: -1 }`. `pierced_multi` had been hiding it: severing A takes an
/// edge of A through B.)
#[test]
fn cut_rod_by_l_severs_into_two() {
    let (mut m, l, rod) = l_and_rod();
    let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
}

/// An edge of one convex solid, threading the other, severs it. The bar runs through the
/// cube and out both ends, so `Cut(bar, cube)` leaves the bar in two 1×1×1 stubs — two solids
/// (cell 0.4), each a clean genus-0 box (`validate` clean, pieces share nothing). The convex
/// path rejected this as `poke_through`; the seam path returns both pieces. (Before cell 0.4
/// this was `disconnected_result`; cell 3e-3's non-convex sibling was `Cut(rod, L)`.)
#[test]
fn a_convex_cut_severs_its_operand() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let bar = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let solids = boolean(&mut m, BoolKind::Cut, bar, a).unwrap();
    assert_eq!(solids.len(), 2);
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    for &s in &solids {
        assert!((nacre_props::mass_props(&m, s).unwrap().volume - 1.0).abs() < 1e-9);
    }
}

#[test]
fn cut_a_seated_block_by_the_part_below_it() {
    // The operands swapped: now the canonical contact face is the upper block's *lower* cap, so
    // the separation test runs with the plane's normal the other way round. Same answer — the
    // block keeps its volume.
    let mut m = Model::new();
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 1.0]),
        Point3::from_array([1.5, 1.5, 2.0]),
    );
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, block, base).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
}

/// **An empty result is an answer.** Two solids that miss each other have no intersection, and
/// that is what `Common` reports: `Ok` with no solids, both operands consumed like any other
/// successful boolean. Stated on its own because the name is the contract — if someone makes
/// this an error again, the failure points straight at what was decided (2026-07-22), and the
/// `live_solids` assertion pins the retire that an early return would otherwise skip.
#[test]
fn a_disjoint_common_is_empty_not_an_error() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    let solids = boolean(&mut m, BoolKind::Common, a, b).expect("empty is not a failure");
    assert!(solids.is_empty());
    assert!(
        m.live_solids().is_empty(),
        "a successful boolean consumes its operands"
    );
}
