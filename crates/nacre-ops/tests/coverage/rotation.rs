//! Rotated operands — rotation-invariance of the boolean.

#![allow(unused_imports)]
use crate::common::*;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean, replay,
};
use nacre_scalar::Axis;
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Orientation, Shell, Solid, Vertex};

#[test]
fn rotated_corner_bite_cut_is_rotation_invariant() {
    let (mut m, l, bx) = l_and_corner_box();
    let l = xf(&mut m, l, rot30());
    let bx = xf(&mut m, bx, rot30());
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_corner_bite_fuse_is_rotation_invariant() {
    let (mut m, l, bx) = l_and_corner_box();
    let l = xf(&mut m, l, rot30());
    let bx = xf(&mut m, bx, rot30());
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 0.924 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_sever_cut_severs_into_two() {
    let (mut m, l, rod) = l_and_rod();
    let l = xf(&mut m, l, rot30());
    let rod = xf(&mut m, rod, rot30());
    let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
    m.rebuild_adjacency();
    assert_eq!(solids.len(), 2, "rotated sever still yields two solids");
    assert!(nacre_validate::validate(&m).is_empty());
    let vol: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
}

#[test]
fn rotated_containment_cut_makes_cavity() {
    let (mut m, a, b) = l_and_inner_box();
    let a = xf(&mut m, a, rot30());
    let b = xf(&mut m, b, rot30());
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(
        m.solids.get(r).cavities.len(),
        1,
        "the inner box becomes a cavity"
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn a_rotated_profile_is_the_same_solid_to_the_boolean() {
    // Rotating a profile's winding cannot change the solid. The cap here carries the
    // dimple's hole, so an inward `n_out` would reverse it — only `validate` and the
    // signed volume would ever say so.
    let (mut m, l) = rotated_l_prism();
    let stub = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_coplanar_contact_is_never_silently_wrong() {
    let tilt = |m: &mut Model, s| xf(m, s, rot_iso(Axis::X, 30));
    // boss fuse (correct fused volume 1.25 if solved).
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let boss = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 1.0]),
        Point3::from_array([0.75, 0.75, 2.0]),
    );
    let (base, boss) = (tilt(&mut m, base), tilt(&mut m, boss));
    if let Ok(r) = boolean_one(&mut m, BoolKind::Fuse, base, boss) {
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "solved boss is valid");
        let v = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((v - 1.25).abs() < 1e-9, "solved boss fuse is correct: {v}");
    }
    // pocket cut (correct carved volume 0.875 if solved).
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let prism = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let (base, prism) = (tilt(&mut m, base), tilt(&mut m, prism));
    if let Ok(r) = boolean_one(&mut m, BoolKind::Cut, base, prism) {
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "solved pocket is valid");
        let v = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((v - 0.875).abs() < 1e-9, "solved pocket cut is correct: {v}");
    }
}
