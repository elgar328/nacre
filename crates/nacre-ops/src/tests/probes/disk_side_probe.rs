//! **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `disk_side_probe`, which is what keeps
//! `super::` here meaning that module — the file's place in the tree does not change what it
//! belongs to.

use crate::ledger::Ledger;

/// The sentence the disk-side check panics with — one spelling, so a commuting oracle `KNOWN`
/// entry can name the site by the constant the check prints.
pub(crate) const NOT_OWN_SOLID: &str = "a cut circle's disk side does not carry its own solid";

/// One entry per arc.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Row {
    /// The cell's corners agreed on a side, so the rule was checked against the geometry here.
    pub(crate) checked: bool,
    /// Both sides spoke — the strongest form of the check.
    pub(crate) both_spoke: bool,
    /// The class's `frame_sign` was `-1`, so the factor decided this arc.
    pub(crate) frame_negative: bool,
}

pub(crate) static ROWS: Ledger<Row> = Ledger::new();

/// Check the disk-side rule against the geometry, and record the row. Each side arrives as
/// `cell_side`'s answer for the arc's cell on that half-edge (`None` = no corner spoke, the corners
/// disagree, or the cell was not found).
///
/// Where the geometry speaks it must **agree with the derived rule** — the check that keeps the
/// derivation honest. Two cells claiming the same side is the witness failing rather than the arc
/// (below), and then both abstain.
pub(crate) fn record(even: Option<bool>, odd: Option<bool>, rule_says_even: bool, frame_sign: i8) {
    let (a, b) = (even, odd);
    // **Two cells on one side is the witness failing, not the arc**: a **convex** arc — a fillet's
    // quarter, a slot's end — bounds a cell that lies inside the circle at the arc and reaches far
    // outside it, so every rational corner of that cell is outside while the arc bounds it from the
    // disk side. The corner witness cannot see a side there; where the two cells' corners agree,
    // both abstain and the rule goes unchecked for that arc (the volume oracles of the tangent
    // fixtures are what measure it). The bite population — arcs concave into a plate — keeps its
    // witness.
    let (a, b) = match (a, b) {
        (Some(x), Some(y)) if x == y => (None, None),
        other => other,
    };
    // The geometry's verdict on which half-edge borders the disk, where it has one.
    let witness = a.or_else(|| b.map(|inside| !inside));
    if let Some(even_is_disk) = witness {
        assert_eq!(
            even_is_disk, rule_says_even,
            "the disk-side rule and the cell's own corners disagree (frame_sign {frame_sign}, even {a:?}, odd {b:?})"
        );
    }
    ROWS.push(Row {
        checked: witness.is_some(),
        both_spoke: a.is_some() && b.is_some(),
        frame_negative: frame_sign < 0,
    });
}
