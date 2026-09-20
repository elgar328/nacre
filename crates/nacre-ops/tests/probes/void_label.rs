//! **What makes a component material, and what makes it a void.**
//!
//! The label used to be read off the normals at the component's lexicographically-minimal vertex:
//! outward iff *some* face there faced −x. That existential is a shortcut for "the face with the
//! largest `|n_x|` faces −x", and the shortcut only holds while the normals are **axis-aligned**.
//! A slanted sketch breaks it with no rotation in sight — and then a void came back as its own
//! *material* solid of negative volume, with the operand unchanged beside it. `validate` had
//! nothing to say about it, which is why the propositions here are about **what the model is**,
//! not about whether the boolean returned `Ok`.
//!
//! The label is containment depth now — even is material, odd is a void — the rule
//! `sketch::from_rings` already uses one dimension down.

use nacre_exact::{Axis, Isometry, Rat};
use nacre_math::{Point2, Point3};
use nacre_ops::{BoolKind, Operation, Profile2d, SketchFrame, apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn prism(m: &mut Model, pts: &[[f64; 2]], h: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile,
            dist: h,
        },
    )
    .expect("extrude") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// Raise a solid by `num/den` along z — the void fixtures need the cutter to stop short of both
/// caps, which is what makes it an enclosed cavity rather than a through pocket.
fn lift(m: &mut Model, s: Handle<Solid>, num: i128, den: i128) -> Handle<Solid> {
    let nacre_ops::OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::translation([
                Rat::from_int(0),
                Rat::from_int(0),
                Rat::new(num, den).unwrap(),
            ]),
        },
    )
    .expect("translate") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// One solid, one cavity, the volume the set difference has — and **no body wound inside out**.
/// The last one is the proposition the old label broke: it returned the cavity as a second
/// "solid" whose volume was negative, and every other check was happy with that.
fn one_solid_with_a_cavity(m: &mut Model, got: Vec<Handle<Solid>>, volume: f64) {
    m.rebuild_adjacency();
    assert_eq!(got.len(), 1, "one body, not the operand plus its own void");
    let s = m.solid(got[0]);
    assert_eq!(s.cavities.len(), 1, "the void is this solid's cavity");
    let v = nacre_props::mass_props(m, got[0]).expect("props").volume;
    assert!(
        v > 0.0,
        "a body wound inside out came back as a solid: volume {v}"
    );
    assert!((v - volume).abs() < 1e-9, "volume {v}, expected {volume}");
    assert!(nacre_validate::validate(m).is_empty());
}

/// ★ The defect: a **slanted** wedge cut wholly inside a cube, touching nothing at all. Before the
/// parity label this returned two bodies — the untouched cube (volume 1) and the wedge wound
/// inside out (volume −0.108) — whose volumes still summed to the right answer.
#[test]
fn a_slanted_void_is_a_cavity_not_a_second_body() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let b = prism(&mut m, &[[0.8, 0.5], [0.2, 0.2], [0.2, 0.8]], 0.6);
    let b = lift(&mut m, b, 2, 10);
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the cut builds");
    one_solid_with_a_cavity(&mut m, got, 1.0 - 0.18 * 0.6);
}

/// The negative control: the **axis-aligned** void, which is the corpus the old rule was written
/// for and which it got right. If this moves, the new label is wrong in the easy direction.
#[test]
fn an_axis_aligned_void_is_unchanged() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.2, 0.2, 0.2]),
        Point3::from_array([0.8, 0.8, 0.8]),
    );
    m.rebuild_adjacency();
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the cut builds");
    one_solid_with_a_cavity(&mut m, got, 1.0 - 0.6 * 0.6 * 0.6);
}

/// The wedge whose sharp tip stops just **short** of the cube's wall — a plain modelling move, and
/// the shape that showed the defect first.
///
/// ★ The tip was at `1.0`, *on* the wall, until the kernel learned to refuse a solid whose surface
/// touches itself: that version is now `SelfTouchingResult`, and it lives in `self_touch.rs` as the
/// case that reject is for. Moving it to `0.999` keeps this file asking what it was always asking —
/// that a cut wholly inside a body comes back as one body with a cavity, whatever the sketch's
/// angles — without also asserting the self-contact rule, which belongs next door.
#[test]
fn a_wedge_whose_tip_nearly_reaches_the_wall_is_still_one_body() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let b = prism(&mut m, &[[0.999, 0.5], [0.1, 0.1], [0.1, 0.9]], 0.6);
    let b = lift(&mut m, b, 2, 10);
    let got = boolean(&mut m, BoolKind::Cut, a, b).expect("the cut builds");
    // Area of the triangle (0.999,0.5), (0.1,0.1), (0.1,0.9): base 0.8 on x = 0.1, height 0.899.
    one_solid_with_a_cavity(&mut m, got, 1.0 - 0.5 * 0.8 * 0.899 * 0.6);
}
