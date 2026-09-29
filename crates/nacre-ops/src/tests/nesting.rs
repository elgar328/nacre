//! Tests for `nesting`'s private reading of a road's refusal.

use super::*;

fn refused(reason: RejectReason) -> BoolError {
    // Built directly: `reject()` would record a raise in the reject census.
    BoolError::Rejected { reason, at: None }
}

/// **The ray road's refusals, read one way each** — the witness's circumstance, the road's, or an
/// error that is not a witness's to absorb. Every name the road raises, under every combination of
/// the two ring facts the reading depends on.
#[test]
fn the_ray_roads_refusals_are_read_by_whose_they_are() {
    use RejectReason::*;
    for has_pierce in [false, true] {
        for mixed in [false, true] {
            let at = format!("has_pierce = {has_pierce}, mixed = {mixed}");
            let said = |r| Said::of_ray_error(refused(r), has_pierce, mixed);
            // A witness on the ring, or one whose rays all graze: the next witness is the remedy.
            assert!(matches!(said(PointOnRing), Ok(Said::Abstain)), "{at}");
            assert!(matches!(said(NoClearRay), Ok(Said::Abstain)), "{at}");
            // A pierce corner no ray plane can side: the ring's, so the road declines. Without one
            // the name came from an arc on the ray's own plane, which the witness chose.
            if has_pierce {
                assert!(
                    matches!(said(PierceVertexUnnamed), Ok(Said::Declines)),
                    "{at}"
                );
            } else {
                assert!(
                    matches!(said(PierceVertexUnnamed), Ok(Said::Abstain)),
                    "{at}"
                );
            }
            // Fewer than three edges: a lens when mixed, a producer defect when straight.
            if mixed {
                assert!(matches!(said(DegenerateRing), Ok(Said::Declines)), "{at}");
            } else {
                assert!(
                    matches!(
                        said(DegenerateRing),
                        Err(BoolError::Rejected {
                            reason: DegenerateRing,
                            ..
                        })
                    ),
                    "{at}"
                );
            }
            // Names the road cannot read are not a witness's circumstance.
            assert!(
                matches!(
                    said(RingNaming),
                    Err(BoolError::Rejected {
                        reason: RingNaming,
                        ..
                    })
                ),
                "{at}"
            );
        }
    }
}
