//! `Operation::Mirror` — reflection in a coordinate plane.
//!
//! Driven through the public front door (`apply`). A reflection is the one motion the kernel
//! could not express: `Isometry` holds proper motions only, so no amount of rotating and moving
//! produces a mirror image. What has to be right is the **orientation algebra** — reflecting the
//! coordinates alone turns a solid inside out, because the winding that defines "outward" flips
//! with `det = −1` — and `validate` cannot see that mistake, so the nets here are the *signed*
//! volume and the boolean.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolKind, OpError, OpOutput, Operation, apply};
use nacre_scalar::{Axis, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn mirror_solid(m: &mut Model, s: Handle<Solid>, axis: Axis, offset: i128) -> Handle<Solid> {
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid: s,
            axis,
            offset: Rat::from_int(offset),
        },
    )
    .unwrap() else {
        unreachable!("Mirror yields a Mirror output")
    };
    m.rebuild_adjacency();
    solid
}

fn copy_solid(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let OpOutput::Copy { solid } = apply(m, &Operation::Copy { solid: s }).unwrap() else {
        unreachable!("Copy yields a Copy output")
    };
    m.rebuild_adjacency();
    solid
}

/// An L-shaped prism in `x ∈ [0, 2]`, deliberately **not** symmetric about any coordinate plane —
/// a symmetric fixture would let a broken mirror pass by doing nothing.
fn ell(m: &mut Model) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            plane: nacre_ops::SketchPlane::world_xy(),
            profile: nacre_ops::Profile2d::polygon(vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 3.0),
                p2(0.0, 3.0),
            ]),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("Extrude yields an Extrude output")
    };
    m.rebuild_adjacency();
    solid
}

/// Volume and area are reflection-invariant, and the volume's **sign** is the check that matters:
/// `validate` has no inside-out test (its violations are all per-cell), so a solid turned inside
/// out would pass validation while the divergence-theorem volume goes negative.
#[test]
fn a_mirror_preserves_volume_and_stays_outward() {
    let mut m = Model::new();
    let a = ell(&mut m);
    let before = nacre_props::mass_props(&m, a).unwrap();

    let b = mirror_solid(&mut m, a, Axis::X, 0);

    let after = nacre_props::mass_props(&m, b).unwrap();
    assert!(after.volume > 0.0, "not inside out: {}", after.volume);
    assert!((after.volume - before.volume).abs() < 1e-12);
    assert!((after.area - before.area).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());
    assert!(!m.live_solids.contains(&a), "mirror supersedes its input");
}

/// The mirror actually moves the shape: an asymmetric solid's coordinates must change. Without
/// this, a mirror that quietly did nothing would pass every other test here.
#[test]
fn a_mirror_is_not_a_no_op() {
    let mut m = Model::new();
    let a = ell(&mut m);
    let before = outer_points(&m, a);

    let b = mirror_solid(&mut m, a, Axis::X, 0);

    let after = outer_points(&m, b);
    assert_ne!(before, after, "an asymmetric solid must change");
    // Every mirrored point is the source point with x negated.
    let mut want: Vec<[f64; 3]> = before.iter().map(|p| [-p[0], p[1], p[2]]).collect();
    let mut got = after;
    want.sort_by(|a, b| a.partial_cmp(b).unwrap());
    got.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(got, want);
}

/// Mirroring twice about the origin plane returns the coordinates **bit for bit**: reflecting
/// across `x = 0` is a sign flip, which f64 does exactly. (An offset plane computes `2c − x` and
/// rounds once, so it gets the tolerance comparison in `an_offset_mirror_plane` instead.)
#[test]
fn mirroring_twice_about_the_origin_is_exact() {
    let mut m = Model::new();
    let a = ell(&mut m);
    let before = outer_points(&m, a);

    let b = mirror_solid(&mut m, a, Axis::Y, 0);
    let c = mirror_solid(&mut m, b, Axis::Y, 0);

    assert_eq!(outer_points(&m, c), before, "bit-identical round trip");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A mirrored solid is a proper operand: if the orientation algebra were wrong the result would
/// be inside out, and the boolean's own outward test would reject or produce nonsense.
#[test]
fn a_mirrored_solid_is_a_usable_operand() {
    let mut m = Model::new();
    let a = ell(&mut m);
    let b = mirror_solid(&mut m, a, Axis::X, 0);
    // A knife across the mirrored L's arm.
    let knife = m.add_cuboid(
        Point3::from_array([-3.0, -1.0, 0.5]),
        Point3::from_array([3.0, 4.0, 2.0]),
    );
    m.rebuild_adjacency();

    let r = boolean_one(&mut m, BoolKind::Cut, b, knife).unwrap();
    m.rebuild_adjacency();

    // The L is 4 units of area, 1 tall; the knife takes everything above z = 0.5.
    assert!((volume(&m, r) - 2.0).abs() < 1e-12, "{}", volume(&m, r));
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **Why this operation exists.** Design half a symmetric part, mirror it, join the two — and
/// because the half touches the mirror plane, the two halves meet there as a *coplanar contact*,
/// the hardest path in the boolean. The interface must dissolve rather than survive as an
/// interior face.
#[test]
fn a_half_and_its_mirror_fuse_into_a_symmetric_part() {
    let mut m = Model::new();
    let half = ell(&mut m); // x ∈ [0, 2], touching x = 0 but not crossing it
    let v_half = volume(&m, half);

    let twin = copy_solid(&mut m, half);
    let other = mirror_solid(&mut m, twin, Axis::X, 0); // x ∈ [−2, 0]
    let whole = boolean_one(&mut m, BoolKind::Fuse, half, other).unwrap();
    m.rebuild_adjacency();

    assert!(
        (volume(&m, whole) - 2.0 * v_half).abs() < 1e-12,
        "symmetric part is twice the half: {} vs {}",
        volume(&m, whole),
        2.0 * v_half
    );
    assert!(
        !has_face_on_plane(
            &m,
            whole,
            Point3::from_array([0.0, 0.5, 0.5]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        ),
        "the joining plane must not survive as a face"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A boolean result carries `Discovered` vertices whose definitions name three surfaces; mirroring
/// has to remap those onto the reflected surfaces, and the result must still be usable.
#[test]
fn a_boolean_result_mirrors() {
    let (mut m, a, b) = two_boxes();
    m.rebuild_adjacency();
    let fused = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    let before = volume(&m, fused);

    let mirrored = mirror_solid(&mut m, fused, Axis::Z, 0);

    assert!((volume(&m, mirrored) - before).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());

    let knife = m.add_cuboid(
        Point3::from_array([-2.0, -2.0, -2.0]),
        Point3::from_array([2.0, 2.0, -0.75]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, mirrored, knife).unwrap();
    m.rebuild_adjacency();
    assert!(volume(&m, r) > 0.0);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A cavity is a second shell; the walk has to reflect and rewind it like any other.
#[test]
fn a_hollow_solid_mirrors_with_its_cavity() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);

    let mirrored = mirror_solid(&mut m, hollow, Axis::X, 0);

    assert_eq!(m.solids.get(mirrored).cavities.len(), 1, "void survives");
    assert!((volume(&m, mirrored) - 26.0).abs() < 1e-9);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// An offset plane places the image where the geometry says. The coordinates go through
/// `2c − x`, which rounds once, so this compares within a tolerance rather than bit-for-bit.
#[test]
fn an_offset_mirror_plane() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    m.rebuild_adjacency();

    let b = mirror_solid(&mut m, a, Axis::X, 3); // x ∈ [0,1] ↦ x ∈ [5,6]

    let xs: Vec<f64> = outer_points(&m, b).iter().map(|p| p[0]).collect();
    assert!(
        xs.iter()
            .all(|&x| (x - 5.0).abs() < 1e-12 || (x - 6.0).abs() < 1e-12),
        "{xs:?}"
    );
    assert!((volume(&m, b) - 1.0).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Curved geometry is declined rather than guessed: a reflection reverses a circle's
/// parametrisation, and that convention is settled with the curved-geometry milestone.
#[test]
fn mirroring_a_cylinder_is_declined() {
    let mut m = Model::new();
    let c = m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    m.rebuild_adjacency();

    let err = apply(
        &mut m,
        &Operation::Mirror {
            solid: c,
            axis: Axis::X,
            offset: Rat::from_int(0),
        },
    )
    .unwrap_err();
    assert!(matches!(err, OpError::MirrorNotPlanar), "{err:?}");
}

/// A rotated solid mirrors: the reflection extends the chain with a `Mirror` node, so the image
/// keeps an exact `Origin::Moved` definition rather than degrading to bare coordinates.
#[test]
fn a_rotated_solid_mirrors() {
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    m.rebuild_adjacency();
    let r = xf(&mut m, c, rot30());
    m.rebuild_adjacency();
    let before = volume(&m, r);

    let mirrored = mirror_solid(&mut m, r, Axis::X, 0);

    assert!((volume(&m, mirrored) - before).abs() < 1e-12);
    let rotated_verts = m
        .shells
        .get(m.solids.get(mirrored).outer)
        .faces
        .iter()
        .flat_map(|&fh| m.faces.get(fh).outer.half_edges.clone())
        .filter_map(|he| m.edges.get(he.edge).bounds)
        .flatten()
        .filter(|&vh| matches!(m.vertices.get(vh).origin, nacre_topo::Origin::Moved { .. }))
        .count();
    assert!(rotated_verts > 0, "the image keeps its rotation provenance");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// **The improper chain's actual test.** A mirrored rotated vertex is described by a chain that
/// ends in a reflection, and it must still land exactly where a plain reflection would put it. A
/// chain that is wrong in a *rigid* way (a sign on an angle, a pivot that did not move) yields a
/// solid merely turned away from the truth: same volume, valid topology, wrong place. Only
/// comparing against the independent reflection catches that.
#[test]
fn a_mirrored_rotated_solid_lands_where_reflection_says() {
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let r = xf(&mut m, c, rot30());
    m.rebuild_adjacency();
    let source = outer_points(&m, r);

    let mirrored = mirror_solid(&mut m, r, Axis::X, 0);

    let mut want: Vec<[f64; 3]> = source.iter().map(|p| [-p[0], p[1], p[2]]).collect();
    let mut got = outer_points(&m, mirrored);
    let key = |a: &[f64; 3], b: &[f64; 3]| a.partial_cmp(b).unwrap();
    want.sort_by(key);
    got.sort_by(key);
    assert_eq!(want.len(), got.len());
    for (w, g) in want.iter().zip(&got) {
        for k in 0..3 {
            assert!(
                (w[k] - g[k]).abs() < 1e-12,
                "replayed {g:?} vs reflected {w:?}"
            );
        }
    }
}

/// A mirrored *rotated* solid is still a boolean operand. Moved operands are judged on their
/// exact `Pt3` definitions, and this one's chain has **odd parity** — so a handedness correction
/// that is wrong (or missing) shows up as a wrong or refused decision here, where the coordinates
/// alone would look perfectly fine.
#[test]
fn a_mirrored_rotated_solid_is_a_usable_operand() {
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();
    let r = xf(&mut m, c, rot30());
    m.rebuild_adjacency();
    let mirrored = mirror_solid(&mut m, r, Axis::X, 0);

    let knife = m.add_cuboid(
        Point3::from_array([-10.0, -10.0, 1.0]),
        Point3::from_array([10.0, 10.0, 10.0]),
    );
    m.rebuild_adjacency();
    let cut = boolean_one(&mut m, BoolKind::Cut, mirrored, knife).unwrap();
    m.rebuild_adjacency();

    // The knife takes the top half of a 2³ box, whatever its orientation about z.
    assert!((volume(&m, cut) - 4.0).abs() < 1e-9, "{}", volume(&m, cut));
    assert!(nacre_validate::validate(&m).is_empty());
}
