//! How many probes a component's depth classification tried before one
//! decided — `(tried, offered, decided)`, each row owned by the test whose work made it. Expected 1 under the corner rule, more only
//! where a tangent or boundary tie remains.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::assembly`] mounts it with `#[path]` as `probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

pub(crate) mod deciding {
    use crate::combinatorics::tie_probe::Tie;
    use crate::ledger::Ledger;
    /// `(tried, offered, decided, ties of this call when exhausted)`.
    pub(crate) type Row = (usize, usize, bool, Vec<(Tie, usize)>);
    pub(crate) static ROWS: Ledger<Row> = Ledger::new();
    pub(crate) fn record(tried: usize, offered: usize, decided: bool, ties: Vec<(Tie, usize)>) {
        ROWS.push((tried, offered, decided, ties));
    }
}
