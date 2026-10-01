//! What a **consumer** sees when the engine declines.
//!
//! The reason rides in the returned value, so an application can say *why* instead of "it
//! failed". This module is that application's job written as a test — the shape `nacre-kit`
//! will take — and it exists to keep the error surface usable from outside the crate, which
//! is the only place that can prove it.

use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, RejectClass, RejectReason, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// The message an application would show.
///
/// Branching is on [`RejectClass`], never on the variant: reason names are engine vocabulary
/// (plane triples, seam runs), kept for logs and bug reports rather than for control flow.
/// Localization is the application's business — the kernel supplies a stable identifier.
fn user_message(err: BoolError) -> String {
    match err {
        BoolError::InputNotLive => "an operand is no longer part of the model".to_string(),
        BoolError::Rejected { reason, .. } => {
            // Facts only, never advice or prediction — the standing wording rule (a
            // diagnosis, not a prompt; no smuggled tense). These mirror what the kit ships.
            let what = match reason.class() {
                RejectClass::NotSupported => "the kernel does not build this",
                RejectClass::Impossible => "no valid solid exists for this input",
                RejectClass::SuspectedDefect => "an engine invariant broke",
            };
            format!("{what} [{reason}]")
        }
    }
}

/// A fuse that would pinch one solid at a single vertex: A and B meet only at `(2,2,1)` while two
/// bridges run around the contact and join them elsewhere, so the material loops and there is no
/// pair of solids to hand back. The caller learns that much — the reason, its class, and a message
/// — from the returned error alone.
///
/// ★ The fixture is a shape that cannot part — two cubes touching at a corner are two bodies —
/// because what this test is about is the error *surface*.
#[test]
fn a_pinched_fuse_tells_the_caller_why() {
    let mut m = Model::new();
    let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        let s = nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        s
    };
    let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let g1 = cub(&mut m, [1.0, 0.3, 0.2], [4.5, 1.3, 0.8]);
    let g2 = cub(&mut m, [3.5, 0.3, 0.2], [4.5, 3.8, 1.6]);
    let b = cub(&mut m, [2.0, 2.0, 1.0], [4.0, 4.0, 2.0]);
    let t1 = boolean(&mut m, BoolKind::Fuse, a, g1).expect("a and g1 overlap");
    m.rebuild_adjacency();
    let t2 = boolean(&mut m, BoolKind::Fuse, t1[0], g2).expect("g1 and g2 overlap");
    m.rebuild_adjacency();

    let err = boolean(&mut m, BoolKind::Fuse, t2[0], b).unwrap_err();

    let BoolError::Rejected { reason, .. } = err else {
        panic!("expected a named rejection, got {err:?}");
    };
    assert_eq!(reason, RejectReason::NonManifoldVertex);
    // Not a coverage limit and not a defect: no milestone will make this input buildable.
    assert_eq!(reason.class(), RejectClass::Impossible);
    // The identifier is stable across the tag→enum move, so logs and reports keep their meaning.
    assert_eq!(reason.as_str(), "non_manifold_vertex");
    assert_eq!(
        user_message(err),
        "no valid solid exists for this input [non_manifold_vertex]"
    );
}

/// The edge twin of the pinch above: A and B meet only along the line `x = 2, y = 2`, and a bridge
/// overlapping both takes the material around the contact. Four faces use that segment and the body
/// cannot be parted there, which no 2-manifold boundary allows — `Impossible`, not a coverage limit.
///
/// ★ Not two cubes sharing one whole edge: those are two bodies, and the grid proptest scores
/// them.
#[test]
fn a_pinched_fuse_is_impossible_not_unsupported() {
    let mut m = Model::new();
    let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
        let s = nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        s
    };
    let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let bridge = cub(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
    let b = cub(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
    let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
    m.rebuild_adjacency();

    let err = boolean(&mut m, BoolKind::Fuse, ab[0], b).unwrap_err();

    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::NonManifoldResultEdge,
                ..
            }
        ),
        "got {err:?}"
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5; 3]),
        Point3::from_array([1.5; 3]),
    );
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

/// **A rejected boolean leaves the live model exactly as it found it.**
///
/// The operands are retired when the result is accepted, and a reject can be raised on either side
/// of that — `check_result_topology` runs after the whole engine has. So "reject" has to mean the
/// caller still holds the two solids it passed in, whichever step said no; otherwise a refused
/// operation silently costs the user their model.
///
/// Only the *live set* is restored. Result cells the refused attempt appended to the arena stay
/// there, unreachable — the arena is append-only by design and superseded solids leave the same
/// residue.
///
/// ★ **What this cannot show today, stated rather than implied.** Every reject the kernel currently
/// raises happens *before* the retire, so the restore this locks is not what keeps these particular
/// cases whole — they were never at risk. The assertion is here because the property is the
/// contract, and because the next reject to be added is the first one that sits after the retire.
/// The `!= 0` count below is what stops the loop from passing on an empty set of rejects.
#[test]
fn a_rejected_boolean_leaves_the_operands_live() {
    // A pinched pair: A and B meet only along `x = 2, y = 2`, and a bridge overlapping both runs
    // the material around the contact. `Fuse` refuses it; `Cut` and `Common` of the same pair
    // build, so one fixture covers a reject and a success against the same assertion.
    let pinched = || -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
            let s = nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            s
        };
        let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
        let bridge = cub(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
        let b = cub(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
        let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
        m.rebuild_adjacency();
        (m, ab[0], b)
    };
    let mut rejected = 0;
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        let (mut m, a, b) = pinched();
        let before = m.live_solids().to_vec();
        match boolean(&mut m, kind, a, b) {
            Err(e) => {
                rejected += 1;
                assert_eq!(
                    m.live_solids().to_vec(),
                    before,
                    "{kind:?} rejected with {e:?} and changed the live set"
                );
            }
            // A success *must* move the live set on: the operands are consumed.
            Ok(_) => assert_ne!(
                m.live_solids().to_vec(),
                before,
                "{kind:?} succeeded without retiring"
            ),
        }
    }
    assert_ne!(rejected, 0, "no kind rejected — the assertion never ran");
}
