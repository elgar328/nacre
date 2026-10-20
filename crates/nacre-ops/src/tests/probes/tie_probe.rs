//! **How often the mixed road's arc-end arm decided** — one row per arc-end departure it read, so
//! a fixture can say the arm ran. A corner on the ray is decided by the half-open rule, and an arc
//! whose root lands on its own end by the departure there; this counts the second.
//!
//! **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::combinatorics`] mounts it with `#[path]` as `tie_probe`, which is what keeps `super::`
//! here meaning that module — the file's place in the tree does not change what it belongs to.

use crate::ledger::Ledger;

/// One entry per arc-end departure the arc arm *decided*.
pub(crate) static ARC_END: Ledger<()> = Ledger::new();

pub(crate) fn arc_end_decided() {
    ARC_END.push(());
}

/// This test's arc-end decisions so far.
pub(crate) fn arc_end_decisions_mine() -> usize {
    ARC_END.mine_len()
}
