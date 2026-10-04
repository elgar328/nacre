//! **What one lateral face's cycles laid on one plane class** — the shape of the answer
//! `cycle_on_class` produced, station by station.
//!
//! ★ **A ledger, not a check.** Each row states the class's axis station, whether the face's
//! outer answer grazed or crossed, how the face's loops broke down, and how many
//! carved runs became spans. A fixture reads the distribution to see which arms of the lateral
//! road its models actually exercise — a zero here means *unexercised*, never *verified*.
//!
//! ★★ **Nothing downstream reads [`HITS`]**, so the engine's answer is the same with this module
//! compiled out; `cfg(test)` at the mount states that.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `cycle_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use super::{CylOnClass, combinatorics};
use crate::ledger::Ledger;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Hit {
    /// The class's axis station, realized.
    pub t: f64,
    pub outer: Option<CylOnClass>,
    /// `[whole rims, other outer loops, other inner loops]` among the face's loops (whole rims are
    /// not carved).
    pub kinds: [usize; 3],
    pub carved: usize,
    pub spans: usize,
}

pub(crate) static HITS: Ledger<Hit> = Ledger::new();

pub(crate) fn record(
    t: nacre_exact::Rat,
    outer: Option<CylOnClass>,
    rims: usize,
    cycles: &[(bool, combinatorics::LoopRing)],
    carved: usize,
    spans: usize,
) {
    let mut kinds = [rims, 0, 0];
    for &(inner, _) in cycles {
        kinds[if inner { 2 } else { 1 }] += 1;
    }
    HITS.push(Hit {
        t: t.to_f64(),
        outer,
        kinds,
        carved,
        spans,
    });
}
