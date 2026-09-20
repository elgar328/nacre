//! The retired plane-fence extent test, kept for one commit so the rule that replaced it can be
//! **differenced against it** rather than argued equal to it.
//!
//! The two are the same predicate wherever the fence is well posed: the fence plane is the class
//! that pins one end, so it crosses the line exactly *at* that end, and "the side the far endpoint
//! is on" is the half-line from there. Where it is **not** well posed — the far endpoint sitting on
//! the fence plane, so `want == 0` and every nonzero side is read as outside — the two disagree, and
//! that is what the counters below are for.
//!
//! ☑ It expires with the gate: it needs both endpoints' rational coordinates, which is the very
//! demand this cell is removing.

use super::{Judge, MergedSeg, WorkingPlane, combinatorics};
use core::sync::atomic::{AtomicUsize, Ordering};

pub(crate) static ASKED: AtomicUsize = AtomicUsize::new(0);
pub(crate) static DISAGREED: AtomicUsize = AtomicUsize::new(0);

fn by_fences(
    jd: &Judge<'_, WorkingPlane>,
    sg: &MergedSeg,
    ends: [&[nacre_scalar::Rat; 3]; 2],
    meet: &(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal),
) -> Option<bool> {
    for k in 0..2 {
        let e = combinatorics::class_coeffs_rat(jd, sg.end_h[k].class()?)?;
        let far = ends[1 - k];
        let mut at_far = e[3];
        for i in 0..3 {
            at_far = at_far.checked_add(e[i].checked_mul(far[i])?)?;
        }
        let want = at_far.numer().signum();
        let got = match nacre_scalar::quad::plane_side(&e, &meet.0, &meet.1) {
            nacre_scalar::Orient::Positive => 1,
            nacre_scalar::Orient::Negative => -1,
            nacre_scalar::Orient::Zero => 0,
        };
        if got != 0 && got != want {
            return Some(false);
        }
    }
    Some(true)
}

/// Both verdicts for one crossing, counted. `None` from either side is not a disagreement —
/// it is a description that could not be formed, and the two roads decline for different causes.
#[allow(clippy::too_many_arguments)]
pub(super) fn against_the_fences(
    jd: &Judge<'_, WorkingPlane>,
    sg: &MergedSeg,
    ends: [&[nacre_scalar::Rat; 3]; 2],
    n: combinatorics::NodeId,
    def: &nacre_topo::CylinderDef,
    at: &combinatorics::Located<'_>,
    seg_ends: &[combinatorics::OnLine; 2],
    wc: usize,
) {
    let Some((_, cyl, _)) = combinatorics::pierce_name(n) else {
        return;
    };
    let Some(meet) = combinatorics::pierce_meet(jd, cyl, def, n) else {
        return;
    };
    let (a, b) = (
        by_fences(jd, sg, ends, &meet),
        combinatorics::closed_contains(jd, wc, sg.wall, at, seg_ends),
    );
    if let (Some(a), Some(b)) = (a, b) {
        ASKED.fetch_add(1, Ordering::Relaxed);
        if a != b {
            DISAGREED.fetch_add(1, Ordering::Relaxed);
        }
    }
}
