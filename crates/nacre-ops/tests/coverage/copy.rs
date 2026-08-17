//! `Operation::Copy` — the one operation that adds to the live set without removing anything.
//!
//! Driven through the public front door (`apply`), the way a consumer uses it: the kernel could
//! *move* a solid but never *copy* one, so "cut with the same tool twice" and pattern/mirror sugar
//! had no expression. What matters here is that the twin is genuinely independent (no shared
//! cells) and that a vertex's exact definition survives the duplication.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolKind, OpError, OpOutput, Operation, apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, VertexDef};

/// `apply(Copy)` through the public API, asserting the output shape.
fn copy_solid(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let OpOutput::Copy { solid } = apply(m, &Operation::Copy { solid: s }).unwrap() else {
        unreachable!("Copy yields a Copy output")
    };
    m.rebuild_adjacency();
    solid
}

fn unit_cube(m: &mut Model) -> Handle<Solid> {
    m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]))
}

/// Measured (boolean-made) vertices of a solid — the count is how the `ThreePlane` remap is
/// observed (S7: the measured tolerance is what "discovered" now means).
fn discovered_count(m: &Model, s: Handle<Solid>) -> usize {
    let mut n = 0;
    for &fh in &m.shells.get(m.solids.get(s).outer).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            for vh in m.edges.get(he.edge).vertices {
                if m.vertex_tol(vh).is_some() {
                    n += 1;
                }
            }
        }
    }
    n
}

/// The original stays live beside its twin, and two coincident solids are a valid model: the Euler
/// count is per reachable cell (`χ = 4`, `S = 2`, `G = 0`), and nothing checks solids for overlap.
#[test]
fn a_copy_leaves_the_original_live() {
    let mut m = Model::new();
    let a = unit_cube(&mut m);
    m.rebuild_adjacency();

    let b = copy_solid(&mut m, a);

    assert_ne!(a, b, "the twin is a new solid");
    assert_eq!(m.live_solids.len(), 2, "both live: {:?}", m.live_solids);
    assert!(m.live_solids.contains(&a) && m.live_solids.contains(&b));
    assert!((volume(&m, b) - volume(&m, a)).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Independence, measured rather than assumed: cutting the twin cannot touch the original, which
/// is only true if the copy duplicated every cell instead of sharing them.
#[test]
fn cutting_the_copy_leaves_the_original_untouched() {
    let mut m = Model::new();
    let a = unit_cube(&mut m);
    m.rebuild_adjacency();
    let b = copy_solid(&mut m, a);

    let knife = m.add_cuboid(
        Point3::from_array([0.5, -0.5, -0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    m.rebuild_adjacency();
    let cut = boolean_one(&mut m, BoolKind::Cut, b, knife).unwrap();
    m.rebuild_adjacency();

    assert!(
        (volume(&m, cut) - 0.5).abs() < 1e-12,
        "half the twin remains"
    );
    assert!((volume(&m, a) - 1.0).abs() < 1e-12, "the original is whole");
    assert!(nacre_validate::validate(&m).is_empty());
}

/// The reason this operation exists: one tool, two cuts. Without `Copy` the first boolean consumes
/// the tool and the second cannot be written.
#[test]
fn one_tool_cuts_two_parts() {
    let mut m = Model::new();
    let part1 = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let part2 = m.add_cuboid(
        Point3::from_array([5.0, 0.0, 0.0]),
        Point3::from_array([7.0, 2.0, 2.0]),
    );
    let tool = m.add_cuboid(
        Point3::from_array([-0.5, -0.5, 1.5]),
        Point3::from_array([2.5, 2.5, 2.5]),
    );
    m.rebuild_adjacency();

    let tool2 = copy_solid(&mut m, tool);
    let r1 = boolean_one(&mut m, BoolKind::Cut, part1, tool).unwrap();
    m.rebuild_adjacency();
    // The copy is placed over part2 and cuts it the same way.
    let shift = nacre_scalar::Isometry::translation([
        nacre_scalar::Rat::from_int(5),
        nacre_scalar::Rat::from_int(0),
        nacre_scalar::Rat::from_int(0),
    ]);
    let tool2 = xf(&mut m, tool2, shift);
    let r2 = boolean_one(&mut m, BoolKind::Cut, part2, tool2).unwrap();
    m.rebuild_adjacency();

    assert!((volume(&m, r1) - 6.0).abs() < 1e-12, "{}", volume(&m, r1));
    assert!((volume(&m, r2) - 6.0).abs() < 1e-12, "{}", volume(&m, r2));
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A boolean *result* carries `Discovered` vertices whose definition names three surfaces; copying
/// has to remap those onto the twin's own surfaces. The copy is then fed back into a boolean,
/// which only works if the remap produced usable definitions.
#[test]
fn a_boolean_result_copies_with_its_discovered_vertices() {
    let (mut m, a, b) = two_boxes();
    m.rebuild_adjacency();
    let fused = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    let before = discovered_count(&m, fused);
    assert!(before > 0, "a fuse produces Discovered corners");

    let twin = copy_solid(&mut m, fused);

    assert_eq!(discovered_count(&m, twin), before, "definitions carried");
    assert!((volume(&m, twin) - volume(&m, fused)).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());

    // Reusable: the twin still behaves as an operand.
    let knife = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 1.2]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, twin, knife).unwrap();
    m.rebuild_adjacency();
    assert!(volume(&m, r) > 0.0);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A rotated solid's vertices keep their definitions, and those definitions name the **moved**
/// surfaces — so the twin's corners are defined against the twin's own planes (S7: the motion
/// lives on the faces, and a vertex follows the planes it names).
#[test]
fn a_rotated_solid_copies_with_its_rotation_origin() {
    let mut m = Model::new();
    let c = unit_cube(&mut m);
    m.rebuild_adjacency();
    let r = xf(&mut m, c, rot30());
    m.rebuild_adjacency();

    let twin = copy_solid(&mut m, r);

    let rotated = |s: Handle<Solid>| {
        m.shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .flat_map(|&fh| m.faces.get(fh).outer.half_edges.clone())
            .flat_map(|he| m.edges.get(he.edge).vertices)
            .filter(|&vh| {
                let VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
                    return false;
                };
                // ★ `any`, not `all`: since the invariant-plane restatement the caps a
                // rotation fixes are world-stated, so every corner names two moved walls
                // and one restated cap — the moved provenance lives on the walls.
                tri.iter().any(|&h| {
                    !matches!(
                        m.surface_truth(h),
                        nacre_topo::SurfaceTruth::Plane { motion: None, .. }
                    )
                })
            })
            .count()
    };
    assert!(
        rotated(r) > 0,
        "the rotated solid's corners name moved planes"
    );
    assert_eq!(rotated(twin), rotated(r), "provenance carried");
    assert!((volume(&m, twin) - volume(&m, r)).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Rotation and boolean-result reuse are each covered on their own; their intersection is not.
/// The expectation is *not* preserved `Rotated` provenance — a result is uniformly `Discovered` —
/// but that the twin survives the remap and is still usable as an operand.
#[test]
fn a_rotated_boolean_result_copies_and_stays_usable() {
    let (mut m, a, b) = two_boxes();
    m.rebuild_adjacency();
    let fused = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    let turned = xf(&mut m, fused, rot30());
    m.rebuild_adjacency();

    let twin = copy_solid(&mut m, turned);

    assert!((volume(&m, twin) - volume(&m, turned)).abs() < 1e-12);
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A cavity is a second shell hanging off the solid; the walk has to carry it, and nothing tested
/// that before (not even for `transform`).
#[test]
fn a_hollow_solid_copies_with_its_cavity() {
    let mut m = Model::new();
    let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();
    let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(m.solids.get(hollow).cavities.len(), 1);

    let twin = copy_solid(&mut m, hollow);

    assert_eq!(m.solids.get(twin).cavities.len(), 1, "the void survives");
    assert!(
        (volume(&m, twin) - 26.0).abs() < 1e-9,
        "{}",
        volume(&m, twin)
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Fully coincident operands: the value semantics `nacre-kit` settled on turn `fuse(a, a)` into a
/// legitimate script, so this input stops being exotic. Idempotence (`A ∪ A = A`) is the answer if
/// the engine takes it; an honest reject would be a fact worth knowing instead.
#[test]
fn fusing_a_solid_with_its_own_copy() {
    let mut m = Model::new();
    let a = unit_cube(&mut m);
    m.rebuild_adjacency();
    let twin = copy_solid(&mut m, a);

    match boolean(&mut m, BoolKind::Fuse, a, twin) {
        Ok(solids) => {
            assert_eq!(solids.len(), 1, "one body");
            m.rebuild_adjacency();
            assert!(
                (volume(&m, solids[0]) - 1.0).abs() < 1e-12,
                "A ∪ A = A, got {}",
                volume(&m, solids[0])
            );
            assert!(nacre_validate::validate(&m).is_empty());
        }
        Err(e) => panic!("coincident fuse rejected: {e:?} — record this as a known limit"),
    }
}

/// A superseded handle is a caller mistake; copying it would resurrect a shape the log says is
/// gone, and hide the bug.
#[test]
fn copying_a_consumed_solid_is_rejected() {
    let (mut m, a, b) = two_boxes();
    m.rebuild_adjacency();
    let _ = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();

    let err = apply(&mut m, &Operation::Copy { solid: a }).unwrap_err();
    assert!(matches!(err, OpError::SolidNotLive), "{err:?}");
}
