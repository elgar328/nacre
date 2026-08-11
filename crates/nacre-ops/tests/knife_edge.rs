//! A knife edge on a foreign face plane — the tangency that broke the arrangement.
//!
//! The failing population is not "rotation": it is an operand **edge** lying exactly in the
//! other operand's face plane with both incident faces leaving that plane to one side. An
//! axis-aligned box pair can never produce it (a box edge on a plane always brings a face with
//! it — a seated contact), which is why the corpus missed it; a rotation's fixed line and a
//! slanted sketch both produce it routinely.
//!
//! Before the fix the same root wore two names, two layers from its cause: the tangency's two
//! same-side `Graze` segments were combined by union into a flip that changes nothing (caught by
//! label propagation as `LabelConflict`), and where the tangency was far from everything it
//! dangled in the skeleton and broke the face walk (`RingOrientation`).

use nacre_math::Point3;
use nacre_ops::{BoolKind, Operation, apply, boolean};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn corner_box(m: &mut Model) -> Handle<Solid> {
    let s = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    m.rebuild_adjacency();
    s
}

fn rot_z(m: &mut Model, s: Handle<Solid>, deg: i128) -> Handle<Solid> {
    let nacre_ops::OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            }),
        },
    )
    .expect("rotate") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

fn translate(m: &mut Model, s: Handle<Solid>, d: [i128; 3]) -> Handle<Solid> {
    let nacre_ops::OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::translation(d.map(Rat::from_int)),
        },
    )
    .expect("translate") else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// ① The phantom: the rotated copy is 300 away — no contact of any kind — yet its fixed edge
/// still lies in the *plane* of A's wall, and the graze it traces there used to dangle in the
/// skeleton and break the face walk. Two separate bodies is the whole answer.
#[test]
fn disjoint_phantom_graze_fuses_to_two_bodies() {
    let mut m = Model::new();
    let a = corner_box(&mut m);
    let b = corner_box(&mut m);
    let b = rot_z(&mut m, b, 40);
    let b = translate(&mut m, b, [300, 0, 0]);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a phantom tangency is not contact");
    assert_eq!(out.len(), 2, "disjoint solids fuse to two bodies");
}

/// ② The pivot on the boundary: the rotation's fixed line is A's own corner edge, so B's knife
/// edge coincides with a true boundary. The union is a 130° material sector around that line —
/// manifold, buildable, and previously `LabelConflict`. Volume is scored against inclusion-
/// exclusion via Cut: |A| + |B| − |A∩B| must equal |A∪B|, all four measured.
#[test]
fn overlapping_rotated_copy_fuses() {
    let volume = |m: &Model, solids: &[Handle<Solid>]| -> f64 {
        solids
            .iter()
            .map(|&s| nacre_props::mass_props(m, s).expect("props").volume)
            .sum()
    };
    let build = || {
        let mut m = Model::new();
        let a = corner_box(&mut m);
        let b = corner_box(&mut m);
        let b = rot_z(&mut m, b, 40);
        (m, a, b)
    };
    let (mut m, a, b) = build();
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("the 130-degree union builds");
    let v_union = volume(&m, &fused);
    let (mut m2, a2, b2) = build();
    let common = boolean(&mut m2, BoolKind::Common, a2, b2).expect("the overlap builds");
    let v_common = volume(&m2, &common);
    let box_v = 1000.0;
    assert!(
        (v_union - (2.0 * box_v - v_common)).abs() < 1e-6,
        "inclusion-exclusion: |A∪B| = |A|+|B|−|A∩B| ({v_union} vs {})",
        2.0 * box_v - v_common
    );
    assert!(
        v_union > box_v + 1.0,
        "the union is strictly bigger than one box"
    );
}

/// A triangular prism whose sharp vertical edge is at `(x, y)`, body on the +x side of it.
fn prism(m: &mut Model, x: f64, y: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    use nacre_ops::{Profile2d, SketchFrame};
    let profile = Profile2d::polygon(vec![
        Point2::from_array([x, y]),
        Point2::from_array([x + 10.0, y + 5.0]),
        Point2::from_array([x + 10.0, y - 5.0]),
    ])
    .unwrap();
    let op = Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile,
        dist: 10.0,
    };
    let out = apply(m, &op).expect("prism");
    let nacre_ops::OpOutput::Extrude { solid, .. } = out else {
        panic!()
    };
    m.rebuild_adjacency();
    solid
}

/// ③ No rotation anywhere: a prism whose sharp edge lies in the *plane* of a box's wall but far
/// outside its face — the phantom again, produced by a plain axis-aligned model. The lock that
/// this never was a rotation bug: the knife edge is the whole story.
#[test]
fn axis_aligned_knife_edge_fuses_to_two_bodies() {
    let mut m = Model::new();
    let b = m.add_cuboid(
        Point3::from_array([-30.0, -20.0, 0.0]),
        Point3::from_array([-20.0, 20.0, 10.0]),
    );
    m.rebuild_adjacency();
    // In the plane x = −20, 80 above the face's y-extent: no contact of any kind.
    let w = prism(&mut m, -20.0, 100.0);
    let out =
        boolean(&mut m, BoolKind::Fuse, w, b).expect("an axis-aligned knife is the same case");
    assert_eq!(out.len(), 2);
}

/// ④ The knife edge ON the face itself: the prism touches the box along that line and nowhere
/// else, and the material sectors on the two sides of the plane are not adjacent — the union is
/// genuinely non-manifold along the line. ★ The answer is already in the kernel: the
/// non-manifold-vertex detection sees the edge's endpoints and rejects by name. This is the
/// vertex-touch precedent (corner-coincidence) extended to a line for free — locked here so a
/// future change cannot silently start emitting the non-manifold result.
#[test]
fn edge_only_contact_is_a_named_reject() {
    let mut m = Model::new();
    let b = m.add_cuboid(
        Point3::from_array([-30.0, -20.0, 0.0]),
        Point3::from_array([-20.0, 20.0, 10.0]),
    );
    m.rebuild_adjacency();
    let w = prism(&mut m, -20.0, 0.0); // the sharp edge lies inside the wall's face
    let err = boolean(&mut m, BoolKind::Fuse, w, b).expect_err("a non-manifold union");
    assert!(
        format!("{err:?}").contains("NonManifoldVertex"),
        "named, not silent: {err:?}"
    );
}

/// ⑤ The negative control: the centred pair worked before the fix (the pivot is interior, no
/// edge meets a foreign plane) and its volume must not move.
#[test]
fn centred_rotated_copy_is_untouched() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([-5.0, -5.0, -5.0]),
        Point3::from_array([5.0, 5.0, 5.0]),
    );
    m.rebuild_adjacency();
    let b = m.add_cuboid(
        Point3::from_array([-5.0, -5.0, -5.0]),
        Point3::from_array([5.0, 5.0, 5.0]),
    );
    m.rebuild_adjacency();
    let b = rot_z(&mut m, b, 40);
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("was ok before, stays ok");
    let v: f64 = fused
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    assert!(
        v > 1000.0 && v < 2000.0,
        "a union strictly between one and two boxes: {v}"
    );
}
