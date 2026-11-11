//! A knife edge on a foreign face plane — a tangency the arrangement must survive.
//!
//! The population is not "rotation": it is an operand **edge** lying exactly in the other operand's
//! face plane with both incident faces leaving that plane to one side. An axis-aligned box pair can
//! never produce it (a box edge on a plane always brings a face with it — a seated contact), which
//! is why a corpus of box pairs misses it; a rotation's fixed line and a slanted sketch both
//! produce it routinely.
//!
//! A same-side pair must flip nothing. Combining the tangency's two same-side `Graze` segments by
//! union would make a flip that changes nothing (which label propagation would catch as
//! `LabelConflict`), and where the tangency is far from everything the segment would dangle in the
//! skeleton and break the face walk (`RingOrientation`) — one root wearing two names, two layers
//! from its cause (measured).

use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_ops::{BoolKind, Operation, apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn corner_box(m: &mut Model) -> Handle<Solid> {
    nacre_ops::fixtures::cuboid(m, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
}

fn rot_z(m: &mut Model, s: Handle<Solid>, deg: i128) -> Handle<Solid> {
    let nacre_ops::OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            }),
        },
    )
    .expect("rotate") else {
        panic!()
    };
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
    solid
}

/// ① The phantom: the rotated copy is 300 away — no contact of any kind — yet its fixed edge
/// still lies in the *plane* of A's wall, and the graze it traces there must not dangle in the
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
/// manifold and buildable. Volume is scored against inclusion-
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
    solid
}

/// ③ No rotation anywhere: a prism whose sharp edge lies in the *plane* of a box's wall but far
/// outside its face — the phantom again, produced by a plain axis-aligned model. The lock that the
/// knife edge, not rotation, is the whole story.
#[test]
fn axis_aligned_knife_edge_fuses_to_two_bodies() {
    let mut m = Model::new();
    let b = nacre_ops::fixtures::cuboid(&mut m, [-30.0, -20.0, 0.0], [-20.0, 20.0, 10.0]);
    // In the plane x = −20, 80 above the face's y-extent: no contact of any kind.
    let w = prism(&mut m, -20.0, 100.0);
    let out =
        boolean(&mut m, BoolKind::Fuse, w, b).expect("an axis-aligned knife is the same case");
    assert_eq!(out.len(), 2);
}

/// ④ The knife edge ON the face itself: the prism touches the box along that line and nowhere
/// else. Not a shared ring edge — the prism's sharp edge lies in the *interior* of the box's wall
/// face — so it is a different shape of contact from two cubes sharing a corner edge, and it is
/// checked the way the criterion has to be checked: **by volume**. 4000 + 500 with nothing
/// removed means the interiors are disjoint, which is what "two bodies" claims.
///
/// Not a `NonManifoldVertex` reject: "the material sectors on the two sides of the plane are
/// not adjacent — the union is genuinely non-manifold along the line" is true of a *single* union
/// body, and two bodies is the answer.
#[test]
fn edge_only_contact_separates_into_two_bodies() {
    let mut m = Model::new();
    let b = nacre_ops::fixtures::cuboid(&mut m, [-30.0, -20.0, 0.0], [-20.0, 20.0, 10.0]);
    let w = prism(&mut m, -20.0, 0.0); // the sharp edge lies inside the wall's face
    let out = boolean(&mut m, BoolKind::Fuse, w, b).expect("a line contact separates");
    assert_eq!(out.len(), 2, "a contact along a line is two bodies");
    let total: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    assert!(
        (total - 4500.0).abs() < 1e-9,
        "box 4000 + prism 500, nothing removed — got {total}"
    );
    assert_eq!(
        nacre_validate::validate(&m),
        Vec::new(),
        "each body is a solid on its own"
    );
}

/// ⑤ The negative control: the centred pair (the pivot is interior, no edge meets a foreign plane),
/// whose union lies strictly between one box and two.
#[test]
fn centred_rotated_copy_fuses() {
    let mut m = Model::new();
    let a = nacre_ops::fixtures::cuboid(&mut m, [-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]);
    let b = nacre_ops::fixtures::cuboid(&mut m, [-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]);
    let b = rot_z(&mut m, b, 40);
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("the centred pair fuses");
    let v: f64 = fused
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    assert!(
        v > 1000.0 && v < 2000.0,
        "a union strictly between one and two boxes: {v}"
    );
}

/// A block `x ∈ [−2, 2]`, `z ∈ [0, depth]`, whose top is `y = 5 + |x|/2` — a V notch whose crease
/// is the line `x = 0, y = 5`.
fn notched(m: &mut Model, depth: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    use nacre_ops::{Profile2d, SketchFrame};
    let profile = Profile2d::polygon(vec![
        Point2::from_array([-2.0, 0.0]),
        Point2::from_array([2.0, 0.0]),
        Point2::from_array([2.0, 6.0]),
        Point2::from_array([0.0, 5.0]),
        Point2::from_array([-2.0, 6.0]),
    ])
    .unwrap();
    let op = Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile,
        dist: depth,
    };
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(m, &op).expect("notched block") else {
        panic!()
    };
    solid
}

/// The tool `y ≥ 5`, holding the notch's crease in its floor.
fn above_the_crease(m: &mut Model) -> Handle<Solid> {
    nacre_ops::fixtures::cuboid(m, [-6.0, 5.0, -2.0], [6.0, 12.0, 6.0])
}

/// ⑥ **A knife that parts the kept material is where the face ends.** The notch's crease lies in the
/// tool's floor `y = 5` with both notch faces leaving it upward, so no label changes across it —
/// yet the part of the block above `y = 5` is two wings that meet along the crease alone. Their
/// common part is two bodies (each `∫₀² x/2 dx · 4 = 4`), the floor under them two faces; the other
/// kinds build whole: block 88, tool 672, common 8.
///
/// red: knives never walked (or the coplanar merge joining across the crease) → the floor is one
/// face under both wings and the common part is refused `SelfTouchingResult`.
#[test]
fn a_knife_parting_the_kept_material_splits_the_face_under_it() {
    for swap in [false, true] {
        for (kind, want) in [
            (BoolKind::Fuse, vec![752.0]),
            (BoolKind::Cut, if swap { vec![664.0] } else { vec![80.0] }),
            (BoolKind::Common, vec![4.0, 4.0]),
        ] {
            let mut m = Model::new();
            let b = notched(&mut m, 4.0);
            let t = above_the_crease(&mut m);
            let (x, y) = if swap { (t, b) } else { (b, t) };
            let out = boolean(&mut m, kind, x, y)
                .unwrap_or_else(|e| panic!("{kind:?} swap={swap}: {e:?}"));
            let got: Vec<f64> = out
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                .collect();
            assert_eq!(got.len(), want.len(), "{kind:?} swap={swap}: {got:?}");
            for (g, w) in got.iter().zip(&want) {
                assert!((g - w).abs() < 1e-9, "{kind:?} swap={swap}: {got:?}");
            }
            assert_eq!(
                nacre_validate::validate(&m),
                Vec::new(),
                "{kind:?} swap={swap}"
            );
        }
    }
}

/// ⑦ **A knife that ends inside the face is the body touching itself.** Notched over `z ∈ [0, 2]`
/// only, the block is whole above `y = 5` over `z ∈ [2, 4]`, so the two wings join there and the
/// crease, ending at `z = 2` inside the floor, cannot cut it in two. The common part is one body
/// whose two wings meet along the crease — refused `SelfTouchingResult` at the crease, in both
/// orders.
#[test]
fn a_knife_ending_inside_the_face_is_a_self_touch() {
    for swap in [false, true] {
        let mut m = Model::new();
        let n = notched(&mut m, 2.0);
        let top = nacre_ops::fixtures::cuboid(&mut m, [-2.0, 0.0, 2.0], [2.0, 6.0, 4.0]);
        let b = boolean(&mut m, BoolKind::Fuse, n, top).expect("the half-notched block")[0];
        let t = above_the_crease(&mut m);
        let (x, y) = if swap { (t, b) } else { (b, t) };
        match boolean(&mut m, BoolKind::Common, x, y) {
            Err(nacre_ops::BoolError::Rejected {
                reason: nacre_ops::RejectReason::SelfTouchingResult,
                at: Some(nacre_ops::RejectWhere::Segment([p, q])),
            }) => {
                let mut ends = [p.as_array(), q.as_array()];
                ends.sort_by(|a, b| a[2].total_cmp(&b[2]));
                assert_eq!(ends, [[0.0, 5.0, 0.0], [0.0, 5.0, 2.0]], "swap={swap}");
            }
            other => panic!("swap={swap}: {other:?}"),
        }
    }
}
