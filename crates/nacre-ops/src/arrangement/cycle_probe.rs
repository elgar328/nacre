use super::{CylOnClass, combinatorics};
use std::sync::Mutex;

#[derive(Clone, Copy, Debug)]
#[cfg_attr(not(test), allow(dead_code))]
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
    t: nacre_scalar::Rat,
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
