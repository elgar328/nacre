//! What the kernel had to assume, reported — `boolean_with_report`.
//!
//! Rotated geometry has no exact zero, so "these two faces are one plane" is proved *to within a
//! distance*, never outright. These fix what the report says about that: that it says nothing when
//! nothing was assumed, that it names the plane merges first, and that the closest call it quotes
//! is far below the coincidence limit rather than merely under it.

use crate::common::*;
use nacre_exact::Axis;
use nacre_ops::{BoolKind, Decision, Site, boolean_with_report};

/// **An exact model assumes nothing, and the report says so.**
///
/// Axis-aligned coordinates are exact `f64`, so every predicate answers exactly — including the
/// zeros. If this list were non-empty the channel would be recording proved answers, which would
/// bury the assumed ones it exists to surface.
#[test]
fn an_exact_boolean_reports_nothing_assumed() {
    let (mut m, a, b) = stacked_cubes();
    let (_, r) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).expect("axis-aligned fuse");
    assert_eq!(r.coincidences, 0, "an exact model assumed a coincidence");
    assert!(r.merges.is_empty(), "an exact merge was reported as judged");
    assert!(r.loosest.is_none());
}

/// …and neither does a rotation both operands **share**: the motion cancels out of every
/// determinant, so the answers are exact again. Rotating a model must not turn its exact
/// questions into assumed ones.
#[test]
fn a_shared_rotation_still_assumes_nothing() {
    let (mut m, a, b) = stacked_cubes();
    let a = xf(&mut m, a, rot_iso(Axis::X, 30));
    let b = xf(&mut m, b, rot_iso(Axis::X, 30));
    let (_, r) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).expect("shared-rotation fuse");
    assert_eq!(r.coincidences, 0, "a shared rotation assumed a coincidence");
    assert!(r.merges.is_empty());
}

/// **A plane-class merge decided on toleranced evidence is reported, with the evidence.**
///
/// This is the most consequential judgement the kernel makes — a merge changes which planes exist
/// before a single vertex is computed — and one `Judge::planes_coplanar`'s bare `bool` would
/// hide.
///
/// The fixture is **one motion spelled two ways**: `30°` about X, against the same `30°` followed
/// by a translation and its exact inverse. The net motion is identical, so the planes coincide —
/// but the second chain carries three nodes instead of one, its coordinates are realized through
/// three steps rather than one, and the composition shortcut declines outright on a translation
/// (a plane's `d` is a statement about position, so nothing cancels). The toleranced path decides.
///
/// ★ A pivot moved *along* the rotation axis does not work here, and the reason is worth keeping:
/// it names the same rotation and realizes to the same bits, so the coefficients agree and an
/// earlier, cheaper route answers before this one is reached. The spelling has to differ in the
/// arithmetic, not only in the record.
///
/// ★ **The fixture is a difference the composition cannot cancel.** Same-axis rotations compose
/// (`coplanar_by_composed_rotation` sums same-axis angles), so `30°` against `10°` then `20°` is
/// exact — that pair is pinned in [`the_same_motion_two_ways_is_now_proved`].
#[test]
fn a_toleranced_plane_merge_is_reported_with_its_evidence() {
    use nacre_exact::{Isometry, Rat};
    // Stacked cubes share the plane `z = 1`, and the turn is about X so that plane is *tilted* by
    // it. A plane perpendicular to the rotation axis would keep its coordinate exactly and the
    // merge would be proved rather than judged — a real outcome, but not this one.
    let shift =
        |k: i128| Isometry::translation([Rat::from_int(0), Rat::from_int(k), Rat::from_int(0)]);
    let (mut m, a, b) = stacked_cubes();
    let a = xf(&mut m, a, rot_iso(Axis::X, 30));
    let b = xf(&mut m, b, rot_iso(Axis::X, 30));
    let b = xf(&mut m, b, shift(7));
    let b = xf(&mut m, b, shift(-7));
    let (solids, r) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).expect("fuse");

    // The answer is still right: two unit cubes fused are one solid of volume 2.
    assert_eq!(solids.len(), 1);
    assert!((volume(&m, solids[0]) - 2.0).abs() < 1e-9);

    assert!(
        !r.merges.is_empty(),
        "the shared plane was merged on toleranced evidence, and the report did not say so"
    );
    for e in &r.merges {
        assert!(matches!(e.site, Site::PlanesCoplanar { .. }));
        let Decision::Coincident { within } = e.outcome else {
            panic!("a merge rested on {:?}, which is not evidence", e.outcome);
        };
        // Proved far inside the limit — the point of quoting it is that a reader can check.
        assert!(
            within.exp2().expect("a bound") < -200,
            "merge evidence {:?} is loose enough to want a second look",
            within.exp2()
        );
    }
}

/// ★ **One motion spelled two ways is now *proved*, not judged.**
///
/// `30°` about X against `10°` then `20°` about the same axis and pivot: the chains differ
/// structurally, but they amount to the same rotation, and summing same-axis angles is exact
/// (`Angle` is rational degrees). So the merge rests on `==` between two canonical rational
/// planes, and there is nothing to report.
///
/// This is the other half of [`a_toleranced_plane_merge_is_reported_with_its_evidence`]: that one
/// proves the report speaks when the kernel assumed something, this one proves it stays quiet when
/// the kernel proved it instead. Neither is worth much alone.
#[test]
fn the_same_motion_two_ways_is_now_proved() {
    let (mut m, a, b) = stacked_cubes();
    let a = xf(&mut m, a, rot_iso(Axis::X, 30));
    let b = xf(&mut m, b, rot_iso(Axis::X, 10));
    let b = xf(&mut m, b, rot_iso(Axis::X, 20));
    let (solids, r) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).expect("fuse");
    assert_eq!(solids.len(), 1);
    assert!((volume(&m, solids[0]) - 2.0).abs() < 1e-9);
    assert!(
        r.merges.is_empty(),
        "the shared plane was proved, so nothing should be reported: {:?}",
        r.merges
    );
}

/// **The closest call is quoted, and it is nowhere near the limit.**
///
/// Two operands turned by *different* angles share no cancellation, so hundreds of judgements come
/// down to proved coincidences. What matters is not that they happened but how close the closest
/// one came: the coincidence limit sits at `scale · 2⁻¹⁸⁰`, and the widest bound the kernel had to
/// lean on is some sixty bits below that. A report that ever quoted something near the limit would
/// be telling the model's author to go look at their dimensions.
#[test]
fn the_closest_call_is_reported_and_is_far_inside_the_limit() {
    let (mut m, a, b) = two_boxes();
    let a = xf(&mut m, a, rot_iso(Axis::Z, 30));
    let b = xf(&mut m, b, rot_iso(Axis::Z, 37));
    let (_, r) = boolean_with_report(&mut m, BoolKind::Cut, a, b).expect("mixed-rotation cut");
    assert!(r.coincidences > 0, "a mixed rotation proved no coincidence");
    let e = r.loosest.expect("a widest bound");
    let Decision::Coincident { within } = e.outcome else {
        panic!("the loosest evidence is not a coincidence: {e:?}");
    };
    let exp = within.exp2().expect("a bound");
    assert!(
        (-260..-200).contains(&exp),
        "the closest call was 2^{exp}; the limit for this model is near 2^-178"
    );
}
