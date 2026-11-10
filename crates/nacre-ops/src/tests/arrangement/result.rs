//! The alias rounds' fixed point.

use super::*;

/// **A round's refusal is read only on the settled table** ([`fixed_point`]). Class 0 refuses
/// until the table knows two names as one point; class 1 learns it, in the same round. The first
/// round's table grew, so class 0 is asked again and stands — red when the rounds raise the first
/// refusal they meet. A refusal that no round undoes is still the answer: class 0 alone refuses
/// in a round that learnt nothing.
#[test]
fn a_round_refusal_is_read_on_the_settled_table() {
    let name = |k: usize| NodeId::three_planes(Canon3::three([0, 1, k]));
    let (a, b) = (name(2), name(3));
    // The one union this round can learn, absorbed where the class's table holds it.
    let absorb = |t: &mut Aliases, l: &Aliases| {
        if l.one_point(a, b) {
            t.union_for_test(a, b)
                .map_err(|_| crate::reject(RejectReason::CylinderPairContact))?;
        }
        Ok(())
    };
    let refused = || crate::reject(RejectReason::CoincidentNodes);
    let class0 = |s: &Aliases| {
        if s.one_point(a, b) {
            Ok(((), s.clone()))
        } else {
            Err(refused())
        }
    };
    let learns = |s: &Aliases| {
        let mut l = s.clone();
        l.union_for_test(a, b).expect("one cylinder");
        Ok(((), l))
    };
    let (table, classes) =
        fixed_point(Aliases::default(), absorb, |s| vec![class0(s), learns(s)]).expect("absorbed");
    assert!(table.one_point(a, b));
    assert!(
        classes.iter().all(Result::is_ok),
        "class 0 stands on the settled table: {classes:?}"
    );
    let (_, alone) =
        fixed_point(Aliases::default(), absorb, |s| vec![class0(s)]).expect("absorbed");
    assert!(
        matches!(
            alone[..],
            [Err(BoolError::Rejected {
                reason: RejectReason::CoincidentNodes,
                ..
            })]
        ),
        "nothing learnt, so the refusal stands: {alone:?}"
    );
}
