//! What a **consumer** sees when the engine declines.
//!
//! The reason rides in the returned value, so an application can say *why* instead of "it
//! failed". This module is that application's job written as a test — the shape `nacre-kit`
//! will take — and it exists to keep the error surface usable from outside the crate, which
//! is the only place that can prove it.

use crate::common::{rot_iso, two_boxes, xf};
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, RejectClass, RejectReason, boolean};
use nacre_scalar::Axis;
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

/// **The cause is named, not whatever broke downstream.**
///
/// The original subject was `DegenerateWitness`: rotating one operand of a face-touching pair
/// left a judgement with no distance to measure and no precision that would create one, and
/// before the evidence had a channel out the arrangement read that as "the same point" and failed
/// several steps later as `LoopOrientMismatch` — a symptom of the guess, reported as the problem.
///
/// ★ **That subject went unfired**, and the reason is worth keeping: the only judgement in the
/// corpus that came back degenerate was the component **outwardness** test, and the nesting-parity
/// label replaced it (it asks containment, which the substrate answers, instead of the sign of a
/// rotated plane's normal, which it sometimes cannot). A sweep over three rotation axes, ten
/// angles, two operand shapes, four overlaps and all three kinds — 720 booleans — found no
/// `DegenerateWitness` and no `JudgeExhausted` left. `undecided_reject` is still wired; it has no
/// fixture. See `RejectReason::DegenerateWitness`.
///
/// So the proposition is locked on what the same family still produces. At 60° the pair's `Fuse`
/// stops in the assembly's vertex naming: a corner with no turn, named `StraightAngle` — the
/// cause — rather than the ring that will not close afterwards. `Cut` and `Common` of the same
/// pair still build, which is what makes this a statement about the *judgement* rather than about
/// the geometry as a whole.
#[test]
fn the_cause_is_named_instead_of_its_downstream_symptom() {
    let (mut m, a, b) = two_boxes();
    let a = xf(&mut m, a, rot_iso(Axis::Z, 60));
    let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();
    assert_eq!(
        err,
        BoolError::Unsupported {
            reason: RejectReason::StraightAngle
        },
        "named a symptom instead of the cause"
    );
    for kind in [BoolKind::Cut, BoolKind::Common] {
        let (mut m, a, b) = two_boxes();
        let a = xf(&mut m, a, rot_iso(Axis::Z, 60));
        boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?} of the same pair: {e:?}"));
    }
}

/// **A rejected boolean leaves the live model exactly as it found it.**
///
/// The operands are retired when the result is accepted, and a reject can be raised on either side
/// of that — `check_result_topology` runs after the whole engine has. So "reject" has to mean the
/// caller still holds the two solids it passed in, whichever step said no; otherwise a refused
/// operation silently costs the user their model. The kernel already restored the live set for the
/// topology reject and for nothing else, and nothing asserted the rule at all.
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
    let mut rejected = 0;
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        let (mut m, a, b) = two_boxes();
        // The rotation the arrangement refuses at 60° for `Fuse`; `Cut`/`Common` of the same pair
        // build, so one fixture covers a reject and a success against the same assertion.
        let a = xf(&mut m, a, rot_iso(Axis::Z, 60));
        let before = m.live_solids.clone();
        match boolean(&mut m, kind, a, b) {
            Err(e) => {
                rejected += 1;
                assert_eq!(
                    m.live_solids, before,
                    "{kind:?} rejected with {e:?} and changed the live set"
                );
            }
            // A success *must* move the live set on: the operands are consumed.
            Ok(_) => assert_ne!(m.live_solids, before, "{kind:?} succeeded without retiring"),
        }
    }
    assert_ne!(rejected, 0, "no kind rejected — the assertion never ran");
}
