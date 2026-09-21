//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `disk_side_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use std::sync::Mutex;

/// The sentence the disk-side check panics with — one spelling, shared with the commuting
/// oracle's `KNOWN` list.
pub(crate) const NOT_OWN_SOLID: &str = "a cut circle's disk side does not carry its own solid";

/// One entry per arc.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Row {
    /// A corner named the side, so the rule was checked against the geometry here.
    pub(crate) checked: bool,
    /// Both sides spoke — the strongest form of the check.
    pub(crate) both_spoke: bool,
    /// The class's `frame_sign` was `-1`, so the factor decided this arc.
    pub(crate) frame_negative: bool,
}

pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

/// Assert the premises **and** the rule, and record the row. Each side arrives as its
/// **first** corner that names a side (`None` = no corner spoke, or the cell was not found).
///
/// Two propositions, where the fact is made: an arc's two cells are never on the same side of
/// its circle, and where the geometry speaks it **agrees with the derived rule**. The second
/// is what keeps the derivation honest — the rule was a guess once, and the guess was wrong
/// by one factor.
///
/// ★ A third premise — *no cell has corners on both sides of a circle it cannot straddle* —
/// is what licenses taking the **first** corner as the cell's answer. Scanning every corner
/// to assert it costs a rational three-plane solve per corner per arc: the lib suite measured
/// **526 s** that way against **58 s** taking the first (one session, warm builds). It was
/// measured **0 violations** over the whole suite twice — when the geometric road was built
/// and again here — so it is recorded rather than re-proved on every run.
pub(crate) fn record(even: Option<bool>, odd: Option<bool>, rule_says_even: bool, frame_sign: i8) {
    let (a, b) = (even, odd);
    // ★ **Two cells on one side is the witness failing, not the arc**: a **convex**
    // arc — a fillet's quarter, a slot's end — bounds a cell that lies inside the circle at
    // the arc and reaches far outside it, so every rational corner of that cell is outside
    // while the arc bounds it from the disk side. The corner witness cannot see a side
    // there; where the two cells' corners agree, both abstain and the rule goes unchecked
    // for that arc (the volume oracles of the tangent fixtures are what measure it). The
    // bite population — arcs concave into a plate — keeps its witness exactly as before.
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
    ROWS.lock()
        .expect("the probe's lock is never held across a panic")
        .push(Row {
            checked: witness.is_some(),
            both_spoke: a.is_some() && b.is_some(),
            frame_negative: frame_sign < 0,
        });
}
