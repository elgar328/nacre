//! The indirect cmp_coord: its sign combination and its soundness.

use super::*;
use crate::kernel::frame3::tests::{
    deg, fixture, rand_base, rand_triple, ri, rng, tr, triple_pts, truth,
};

/// `cmp_combine` maps the parity of negative signs to the ordering.
#[test]
fn cmp_combine_counts_negatives() {
    assert_eq!(
        cmp_combine(Some(true), Some(true), Some(true)),
        Some(Orient::Positive)
    );
    assert_eq!(
        cmp_combine(Some(false), Some(true), Some(true)),
        Some(Orient::Negative)
    );
    assert_eq!(
        cmp_combine(Some(false), Some(false), Some(true)),
        Some(Orient::Positive)
    );
    assert_eq!(cmp_combine(None, Some(true), Some(true)), None);
}

/// Indirect cmp_coord soundness over heterogeneous provenance (corpus A) and a
/// near-tie corpus (corpus B: two points sharing an axis coordinate up to ε, rotated
/// about that same axis so the near-tie survives). The judge must never disagree with
/// a GT-stable 512-bit truth; both the fast filter and the escalation are exercised.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn indirect_cmp_coord_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let (mut zero_asserted, mut zero_signed) = (0usize, Vec::<Orient>::new());
    let mut check =
        |a: &[[WitnessPoint; 3]; 3], b: &[[WitnessPoint; 3]; 3], axis: usize, zero: bool| {
            let (ta, tb) = (tr(a), tr(b));
            let judged = indirect_cmp_coord_judge(ta, tb, axis, fixture()).orient();
            let iv = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
                [
                    plane_iv(t[0].0, t[0].1, t[0].2),
                    plane_iv(t[1].0, t[1].1, t[1].2),
                    plane_iv(t[2].0, t[2].1, t[2].2),
                ]
            };
            if cmp_filter(iv(ta), iv(tb), axis).is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            if zero {
                zero_asserted += 1;
                if judged != Orient::Zero {
                    zero_signed.push(judged);
                }
                return;
            }
            match truth::cmp(ta, tb, axis, GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        };

    // Corpus A — heterogeneous provenance (each triple its own rotation and pivot).
    let mut st = 0x6A11_C0DE_5151_2323u64;
    for _ in 0..2000 {
        let a = rand_triple(&mut st);
        let b = rand_triple(&mut st);
        check(&a, &b, rng(&mut st, 0, 2) as usize, false);
    }

    // Corpus B — near-tie in z, shared Z-rotation (a Z-rotation leaves z unchanged, so
    // the pre-rotation z near-equality survives → M near 0 → escalation forced).
    for _ in 0..2000 {
        let z0 = rng(&mut st, -200, 200);
        let eps = match rng(&mut st, 0, 2) {
            0 => ri(0, 1),
            1 => ri(1, rng(&mut st, 5_000, 200_000)),
            _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
        };
        let va = [
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
            ri(z0, 1),
        ];
        let vb = [
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
            ri(z0, 1).checked_add(eps).unwrap(),
        ];
        let angle = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let a = triple_pts(va, &mut st, Axis::Z, angle, piv);
        let b = triple_pts(vb, &mut st, Axis::Z, angle, piv);
        check(&a, &b, 2, eps == ri(0, 1));
    }

    eprintln!(
        "[indirect_cmp_coord] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}; zero_asserted: {zero_asserted}"
    );
    assert!(
        zero_asserted > 0,
        "corpus must exercise the constructed-zero path"
    );
    assert!(
        zero_signed.is_empty(),
        "indirect_cmp_coord: a constructed zero was judged with a sign — {} of {zero_asserted}, first {:?}",
        zero_signed.len(),
        &zero_signed[..zero_signed.len().min(5)]
    );

    assert_eq!(
        wrong, 0,
        "cmp judge must never disagree with GT (soundness)"
    );
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path"
    );
}
