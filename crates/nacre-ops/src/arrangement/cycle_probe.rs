//! **What one lateral face's cycles laid on one plane class** — the shape of the answer
//! `cycle_on_class` produced, station by station.
//!
//! ★ **A ledger, not a check.** Each row states the class's axis station, whether the face's
//! outer answer grazed or crossed, how the face's cycles broke down by kind, and how many
//! carved runs became spans. A fixture reads the distribution to see which arms of the lateral
//! road its models actually exercise — a zero here means *unexercised*, never *verified*.
//!
//! ★★ **Nothing downstream reads [`HITS`]**, so the engine's answer is the same with this module
//! compiled out; `cfg(test)` at the mount states that.

use super::{CylOnClass, combinatorics};
use std::sync::Mutex;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Hit {
    /// The class's axis station, realized.
    pub t: f64,
    pub outer: Option<CylOnClass>,
    /// `[rims, chains, panels, holes]` among the face's cycles (rims are not carved).
    pub kinds: [usize; 4],
    pub carved: usize,
    pub spans: usize,
}

pub(crate) static HITS: Mutex<Vec<Hit>> = Mutex::new(Vec::new());

pub(crate) fn record(
    t: nacre_exact::Rat,
    outer: Option<CylOnClass>,
    rims: usize,
    cycles: &[(combinatorics::CycleKind, combinatorics::LoopRing)],
    carved: usize,
    spans: usize,
) {
    let mut kinds = [rims, 0, 0, 0];
    for (k, _) in cycles {
        kinds[match k {
            combinatorics::CycleKind::Rim => 0,
            combinatorics::CycleKind::Chain => 1,
            combinatorics::CycleKind::Panel => 2,
            combinatorics::CycleKind::Hole => 3,
        }] += 1;
    }
    HITS.lock()
        .expect("the probe's lock is never held across a panic")
        .push(Hit {
            t: t.to_f64(),
            outer,
            kinds,
            carved,
            spans,
        });
}
