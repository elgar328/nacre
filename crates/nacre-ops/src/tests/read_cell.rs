//! The cell reader's own rules, asked of hand-written ends.

use super::*;

/// **A silent cell is named by why its ends are silent, and a defect wins.** The two refusals
/// that read [`silence`] — the emitter's present cell with no chamber and the reader's two silent
/// ends — have no input in the suite or the census, so the rule is held here on ends written by
/// hand.
#[test]
fn a_silent_cell_takes_its_ends_cause_and_a_defect_outranks_the_rest() {
    let (side, above) = (SolidSide::A, [true, true]);
    // The lower end is a chart gap, the upper one has no circle: the defect wins.
    assert_eq!(
        silence(
            &[End::Other(RejectReason::RulingBoundNotYet), End::NoCircle],
            side,
            above
        ),
        RejectReason::CylinderStagesDisagree
    );
    // Neither is a defect: the lower end's.
    assert_eq!(
        silence(
            &[
                End::Other(RejectReason::WitnessNotRational),
                End::Other(RejectReason::RulingBoundNotYet)
            ],
            side,
            above
        ),
        RejectReason::WitnessNotRational
    );
    // An end that is not silent for a reason says nothing here.
    assert_eq!(
        silence(
            &[End::Uncovered, End::Other(RejectReason::ArcBoundNotYet)],
            side,
            above
        ),
        RejectReason::ArcBoundNotYet
    );
    // A disk always speaks: with the other end silent for a reason, that reason.
    assert_eq!(
        silence(
            &[
                End::Disk([false; 4]),
                End::Other(RejectReason::CoincidentNodes)
            ],
            side,
            above
        ),
        RejectReason::CoincidentNodes
    );
}
