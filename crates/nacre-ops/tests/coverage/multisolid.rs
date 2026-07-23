//! Multi-solid output — disjoint operands and face-contact stacks, where a boolean
//! yields several solids (or none), each a clean manifold body.

use crate::common::{boolean_one, stacked_cubes};
use nacre_math::Point3;
use nacre_ops::{BoolKind, boolean};
use nacre_topo::Model;

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
    assert_eq!(m.live_solids, vec![r]);
}

#[test]
fn common_stacked_cubes_is_empty() {
    // The stack shares only its interface plane, so the intersection has no volume.
    let (mut m, a, b) = stacked_cubes();
    assert!(boolean(&mut m, BoolKind::Common, a, b).unwrap().is_empty());
    assert!(m.live_solids.is_empty(), "both operands are consumed");
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
