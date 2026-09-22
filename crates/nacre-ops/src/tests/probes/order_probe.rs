//! A ruler production does not lay, kept so the order rule can be **differenced against it**.
//!
//! ★★★★★ **The two do not agree pointwise, and that is not a defect — it is measured, and it is
//! the reason the difference has to be stated as a *shape* rather than as equality.** The ruler's
//! direction is `n₁ × n₂` built from the classes' **rational coefficients**;
//! [`combinatorics::order_located`] takes its axis sign from the judge's **stored** planes. Those
//! two spellings name the same plane but not the same *side*, so on a class where one of them is
//! stored negated the whole comparison flips. ☑ Measured over the suite: 9,630 of 20,453 segments
//! come out read the other way.
//!
//! So what must hold is not "same answer" but **"same or exactly opposite, per segment"** — a
//! wholesale reversal cancels, because the pieces are emitted in the comparator's own ascending
//! order and `forward` is taken with that same comparator, while a *partial* disagreement would be
//! a genuine reshuffle and would put the sub-segments in the wrong places. That is what is counted.
//!
//! ★ The pair handed to the rule is the **sorted** one, matching `NodeId::Pierce`'s own convention
//! (`planes` ascending, and `root` defined against that order). ☑ Measured: **9,946 of 20,453**
//! segments have `wc > wall`, so the corpus reaches the case where the two orders are opposite —
//! and by the paragraph above it does not matter which is taken, because swapping the pair can only
//! turn "same" into "reversed" for a whole segment at once.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `order_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use super::{Judge, MergedSeg, NodeId, WorkingPlane, along, cmp_along, combinatorics};
use crate::ledger::Ledger;

/// One row per segment the rule ordered, owned by the test whose work ordered it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Row {
    /// Reached with `wc > wall` — where the sorted pair and the call-order pair are opposite
    /// calls. The relation that decides whether the corpus can tell the two apart at all.
    pub(crate) wc_above_wall: bool,
    /// The ruler could be laid and **every** pair compared. The verdicts below describe a
    /// segment only when this holds; the counts still carry what was read before it stopped.
    pub(crate) compared: bool,
    /// The rule read the ruler's order **exactly backwards** — benign, and the common case.
    pub(crate) reversed: bool,
    /// The rule agreed with the ruler on some pairs and not others. **This is the defect.**
    pub(crate) scrambled: bool,
    /// Comparisons both roads called a tie — no direction in them, so they are counted apart.
    pub(crate) eq_both: usize,
    /// Comparisons where one road called two points **the same place** and the other did not — a
    /// claim about coincidence, not about sequence, so it is kept out of the verdicts above.
    pub(crate) equality_disagreed: usize,
}

pub(crate) static ROWS: Ledger<Row> = Ledger::new();

/// The old key: a crossing's parameter on the canonical meet line, an endpoint's `along` on the
/// same line. `None` where the ruler could not be laid — which is the very shape this cell is
/// removing, so it is skipped rather than counted.
fn ruler_keys(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    keyed: &[(combinatorics::PointOn, NodeId, combinatorics::EndPin)],
) -> Option<Vec<nacre_exact::quad::QuadVal>> {
    let mut line = None;
    let mut out = vec![None; keyed.len()];
    for (i, k) in keyed.iter().enumerate() {
        // The endpoints are three-plane named and take the second pass; only a crossing lays
        // the ruler. (Written with `?` at first, which made every call return `None` — the
        // aliveness check below is what said so.)
        let Some((_, cyl, _)) = combinatorics::pierce_name(k.1) else {
            continue;
        };
        let def = &cyls.get(cyl)?.def;
        let (l, s) = combinatorics::pierce_meet(jd, cyl, def, k.1)?;
        line = Some(l);
        out[i] = Some(s);
    }
    let line = line?;
    for (i, k) in keyed.iter().enumerate() {
        if out[i].is_none() {
            out[i] = Some(along(&line, &combinatorics::node_coords_rat(jd, k.1)?)?);
        }
    }
    out.into_iter().collect()
}

pub(crate) fn against_the_ruler(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    sg: &MergedSeg,
    keyed: &[(combinatorics::PointOn, NodeId, combinatorics::EndPin)],
    pair: [usize; 2],
) {
    let mut row = Row {
        wc_above_wall: wc >= sg.wall,
        compared: false,
        reversed: false,
        scrambled: false,
        eq_both: 0,
        equality_disagreed: 0,
    };
    // One row per call, whatever it reaches — a segment the ruler could not be laid on says so
    // with `compared: false` rather than by being absent.
    compare(jd, cyls, keyed, pair, &mut row);
    ROWS.push(row);
}

fn compare(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    keyed: &[(combinatorics::PointOn, NodeId, combinatorics::EndPin)],
    pair: [usize; 2],
    row: &mut Row,
) {
    let Some(keys) = ruler_keys(jd, cyls, keyed) else {
        return;
    };
    let rule = |p: [usize; 2], i: usize, j: usize| -> Option<core::cmp::Ordering> {
        let a = combinatorics::locate(jd, cyls, p[0], p[1], keyed[i].0)?;
        let b = combinatorics::on_line(jd, cyls, p[0], p[1], keyed[j].0)?;
        Some(
            match combinatorics::order_located(jd, p[0], p[1], &a, &b)? {
                -1 => core::cmp::Ordering::Less,
                1 => core::cmp::Ordering::Greater,
                _ => core::cmp::Ordering::Equal,
            },
        )
    };
    use core::cmp::Ordering::Equal;
    let (mut same, mut opposite) = (0usize, 0usize);
    for i in 0..keyed.len() {
        for j in 0..keyed.len() {
            if i == j {
                continue;
            }
            let (Some(want), Some(got)) = (cmp_along(&keys[i], &keys[j]), rule(pair, i, j)) else {
                return;
            };
            // ★ **Coincidence is a different disagreement from order, and is counted apart.**
            // One side calling two points the same place while the other separates them says
            // nothing about the sequence; it says the `CoincidentNodes` refusal is reachable
            // from one road and not the other, which is its own fact.
            if want == Equal && got == Equal {
                // ★ Both roads put the two points in the same place. That is agreement about
                // *coincidence*, and it carries no direction, so counting it as "same order"
                // would make a wholly reversed segment look scrambled. ☑ Measured 0 over the
                // suite: the class is kept because it is a different claim, not because a
                // fixture reaches it.
                row.eq_both += 1;
            } else if (want == Equal) != (got == Equal) {
                row.equality_disagreed += 1;
            } else if got == want {
                same += 1;
            } else {
                opposite += 1;
            }
        }
    }
    row.compared = true;
    match (same, opposite) {
        (0, n) if n > 0 => row.reversed = true,
        (_, 0) => {}
        _ => row.scrambled = true,
    }
}
