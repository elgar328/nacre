//! Convex operands — basic cut/fuse/common of axis-aligned boxes. Pins the
//! in/out sign convention and the volume/manifold invariants end to end.

use crate::common::{boolean_one, two_boxes};
use nacre_math::Point3;
use nacre_ops::{BoolKind, boolean};
use nacre_topo::Model;

#[test]
fn cut_of_two_cubes() {
    // A − B where A = [0,1]³, B = [0.5,1.5]³ ⇒ 1 − 0.125 = 0.875.
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.875).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.live_solids, vec![r]);
}

#[test]
fn fuse_of_two_cubes() {
    // A = [0,1]³, B = [0.5,1.5]³ ⇒ A∪B volume 1+1−0.125 = 1.875.
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 1.875).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.live_solids, vec![r]);
}

#[test]
fn common_of_two_cubes_is_their_overlap() {
    // A = [0,1]³, B = [0.5,1.5]³ ⇒ A∩B = [0.5,1]³ (volume 0.125). This also
    // pins the in/out sign convention end to end: a flipped parity would take
    // the complement and give a wrong (or non-closed) result.
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Common, a, b).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let reach = m.reachable();
    assert_eq!(reach.vertices.len(), 8);
    assert_eq!(reach.edges.len(), 12);
    assert_eq!(reach.faces.len(), 6);
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 0.125).abs() < 1e-12, "volume {vol}");
    assert_eq!(m.live_solids, vec![r]); // A and B superseded
}

#[test]
fn common_with_enclosing_box_is_the_inner_solid() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let c = m.add_cuboid(Point3::from_array([-5.0; 3]), Point3::from_array([5.0; 3]));
    let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
    let r = boolean_one(&mut m, BoolKind::Common, a, c).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let reach = m.reachable();
    assert_eq!(reach.vertices.len(), 8);
    assert_eq!(reach.faces.len(), 6);
    let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol_r - vol_a).abs() < 1e-9, "{vol_r} vs {vol_a}");
}

#[test]
fn common_of_disjoint_boxes_is_empty() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    // Offset in all axes so no faces are coplanar with A.
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    assert!(boolean(&mut m, BoolKind::Common, a, b).unwrap().is_empty());
}

#[test]
fn fuse_common_inclusion_exclusion() {
    // vol(A∪B) + vol(A∩B) == vol(A) + vol(B). boolean supersedes inputs, so
    // fuse and common run on independent copies.
    let corner = |o: f64| {
        (
            Point3::from_array([o; 3]),
            Point3::from_array([o + 1.0, o + 1.0, o + 1.0]),
        )
    };
    let (amin, amax) = corner(0.0);
    let (bmin, bmax) = corner(0.5);
    let vol_of = |mn, mx| {
        let mut m = Model::new();
        let s = m.add_cuboid(mn, mx);
        nacre_props::mass_props(&m, s).unwrap().volume
    };
    let (va, vb) = (vol_of(amin, amax), vol_of(bmin, bmax));

    let mut m1 = Model::new();
    let (a1, b1) = (m1.add_cuboid(amin, amax), m1.add_cuboid(bmin, bmax));
    let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1).unwrap();
    let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;

    let mut m2 = Model::new();
    let (a2, b2) = (m2.add_cuboid(amin, amax), m2.add_cuboid(bmin, bmax));
    let rc = boolean_one(&mut m2, BoolKind::Common, a2, b2).unwrap();
    let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;

    assert!((vf + vc - va - vb).abs() < 1e-9, "{vf}+{vc} vs {va}+{vb}");
}
