//! **Why the mixed road abstains**, one row per `None`, by the site that said it: the
//! probe *at* a ring corner, the probe *on* a step (along the ray or across it), a tangent ray,
//! a root at the probe on its arc, a horizontal tangent at an arc end the root lands on, a
//! zero-span arc — or `Other` for the silent `?` arms (no chart, overflow).
//! Read by the crossing census (per cell) and the ledger (whole suite). A corner on the ray is
//! not among them: the half-open rule decides it, and `ARC_END` counts the arc-end arm's
//! decisions so a fixture can say it ran.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::combinatorics`] mounts it with `#[path]` as `tie_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use crate::ledger::Ledger;
use std::cell::Cell;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Tie {
    ProbeAtCorner,
    ProbeOnStep,
    TangentRay,
    ArcRootAtProbe,
    TangentAtEnd,
    ZeroSpanArc,
    Other,
    /// The lateral road's abstentions: the crossing on a whole-circle rim, on an
    /// arc, on a ruling piece; an arc whose ends name no ⊥ class (a tilted cut); a
    /// loop edge the road cannot read (a plane carrier on a lateral, a corner without a
    /// wall class).
    OnRim,
    OnArc,
    OnRuling,
    TiltedArc,
    Producer,
}

pub(crate) static ROWS: Ledger<Tie> = Ledger::new();

/// One entry per arc-end departure the arc arm *decided*.
pub(crate) static ARC_END: Ledger<()> = Ledger::new();

pub(crate) fn arc_end_decided() {
    ARC_END.push(());
}

/// This test's arc-end decisions so far.
pub(crate) fn arc_end_decisions_mine() -> usize {
    ARC_END.mine_len()
}

thread_local! {
    static LAST: Cell<Option<Tie>> = const { Cell::new(None) };
}

pub(crate) fn begin() {
    LAST.with(|c| c.set(None));
}

pub(crate) fn mark(t: Tie) {
    LAST.with(|c| c.set(Some(t)));
}

/// The lateral road's reading of a mark: `arc_span` marks its ties into `LAST` for the
/// mixed road's wrapper to flush, and the lateral road — which has no wrapper — records
/// the mark as a row at once (or `or` when the abstention was arithmetic, unmarked).
pub(crate) fn flush_or(or: Tie) {
    push(LAST.with(|c| c.take()).unwrap_or(or));
}

/// A row recorded at once — the lateral road has no wrapper to flush `LAST`.
pub(crate) fn push(t: Tie) {
    ROWS.push(t);
}

/// **This test's** rows recorded since `from` (a snapshot of [`len`]), as a histogram by kind —
/// a window alone would carry whatever a parallel test recorded in the same stretch.
pub(crate) fn since(from: usize) -> Vec<(Tie, usize)> {
    let mut hist: Vec<(Tie, usize)> = Vec::new();
    for t in &ROWS.mine_since(from) {
        match hist.iter_mut().find(|(k, _)| k == t) {
            Some((_, c)) => *c += 1,
            None => hist.push((*t, 1)),
        }
    }
    hist
}

pub(crate) fn len() -> usize {
    ROWS.len()
}

pub(crate) fn abstained() {
    let t = LAST.with(|c| c.take()).unwrap_or(Tie::Other);
    ROWS.push(t);
}
