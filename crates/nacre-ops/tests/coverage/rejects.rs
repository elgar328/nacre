//! What a **consumer** sees when the engine declines.
//!
//! The reason rides in the returned value, so an application can say *why* instead of "it
//! failed". This module is that application's job written as a test — the shape `nacre-kit`
//! will take — and it exists to keep the error surface usable from outside the crate, which
//! is the only place that can prove it.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, RejectClass, RejectReason, boolean};
use nacre_scalar::Rat;
use nacre_topo::Model;

/// The message an application would show.
///
/// Branching is on [`RejectClass`], never on the variant: reason names are engine vocabulary
/// (plane triples, seam runs), kept for logs and bug reports rather than for control flow.
/// Localization is the application's business — the kernel supplies a stable identifier.
fn user_message(err: BoolError) -> String {
    match err {
        BoolError::InputNotLive => "an operand is no longer part of the model".to_string(),
        BoolError::Unsupported { reason } => {
            let what = match reason.class() {
                RejectClass::NotSupportedYet => "not supported yet",
                RejectClass::Impossible => "cannot produce a valid solid",
                RejectClass::SuspectedDefect => "could not produce a valid result (please report)",
            };
            format!("{what} [{reason}]")
        }
    }
}

/// Two cubes meeting at exactly one point: the fuse would pinch two solids at a single vertex,
/// which no valid 2-manifold can be. The caller learns that much — the reason, its class, and a
/// message — from the returned error alone.
#[test]
fn a_corner_touching_fuse_tells_the_caller_why() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    m.rebuild_adjacency();

    let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();

    let BoolError::Unsupported { reason } = err else {
        panic!("expected an Unsupported rejection, got {err:?}");
    };
    assert_eq!(reason, RejectReason::NonManifoldVertex);
    // Not a coverage limit and not a defect: no milestone will make this input buildable.
    assert_eq!(reason.class(), RejectClass::Impossible);
    // The identifier is stable across the tag→enum move, so logs and reports keep their meaning.
    assert_eq!(reason.as_str(), "non_manifold_vertex");
    assert_eq!(
        user_message(err),
        "cannot produce a valid solid [non_manifold_vertex]"
    );
}

/// Two cubes meeting along one whole edge — the edge twin of the corner touch above, and the
/// reject the grid proptest lands on most often. Four faces would meet along that edge, which no
/// 2-manifold boundary allows, so it is `Impossible` rather than a coverage limit.
#[test]
fn an_edge_touching_fuse_is_impossible_not_unsupported() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();

    let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();

    assert_eq!(
        err,
        BoolError::Unsupported {
            reason: RejectReason::NonManifoldResultEdge
        }
    );
    assert_eq!(
        RejectReason::NonManifoldResultEdge.class(),
        RejectClass::Impossible
    );
}

/// A handle that no longer names a live solid is a *caller* mistake, not a kernel limit, so it
/// stays its own variant rather than becoming a reason.
#[test]
fn a_stale_operand_is_not_a_reject_reason() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    m.rebuild_adjacency();
    let fused = boolean(&mut m, BoolKind::Fuse, a, b).expect("overlapping cubes fuse");
    m.rebuild_adjacency();

    // `a` was consumed by the fuse above (supersede semantics), so reusing it is stale.
    let err = boolean(&mut m, BoolKind::Fuse, a, fused[0]).unwrap_err();
    assert_eq!(err, BoolError::InputNotLive);
    assert_eq!(
        user_message(err),
        "an operand is no longer part of the model"
    );
}

/// A tool edge lying inside one of the target's planes: three planes share that line, so a fourth
/// crossing it makes a vertex no plane *triple* can name. This is the model reported from the
/// playground — a bar spun 45° next to a step — which is an easy thing to draw by accident.
///
/// The caller must hear the cause — "four planes meet at a point, not supported yet" — and not the
/// symptom the tracer happens to trip over ("could not name a run's bounding node"). This is the
/// firing test for [`RejectReason::FourPlane`], and the switch that flips the day the substrate can
/// express such a vertex: then this boolean succeeds and this test says so.
///
/// **The fused block is load-bearing** — and for a more specific reason than "it adds a face": it
/// supplies the plane **x = 0.5**, and *that* is the plane the bar's bottom corner edge lands in.
/// The bar is 0.4 wide centred on x = 0.5 and 0.4 tall centred on z = 1, so its half-width and its
/// pivot-to-bottom offset are equal (both `0.2`, and bit-identical in `f64`: `0.7 − 0.5` and
/// `1 − 0.8` round the same way) — spun 45°, the corner travels to x = 0.5 exactly. Without the
/// block that plane does not exist and the same cut succeeds; [`crate::rotation`] pins that along
/// with the other near misses.
#[test]
fn a_tool_edge_in_the_target_plane_is_a_four_plane_concurrency() {
    let (mut m, cube, bar) = cube_and_spun_bar(0.2, 45, Rat::from_int(1));
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let target =
        boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses onto the cube")[0];
    m.rebuild_adjacency();

    let err = boolean(&mut m, BoolKind::Cut, target, bar).unwrap_err();

    let BoolError::Unsupported { reason } = err else {
        panic!("expected an Unsupported rejection, got {err:?}");
    };
    assert_eq!(reason, RejectReason::FourPlane);
    // A coverage limit, not a defect and not an impossible input: the cut is perfectly well
    // defined, the arrangement just cannot name its vertices yet.
    assert_eq!(reason.class(), RejectClass::NotSupportedYet);
    assert_eq!(reason.as_str(), "fourplane");
    assert_eq!(user_message(err), "not supported yet [fourplane]");
}
