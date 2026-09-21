//! How many probes a component's depth classification tried before one
//! decided — `(test, tried, offered, decided)`. Expected 1 under the corner rule, more only
//! where a tangent or boundary tie remains.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::assembly`] mounts it with `#[path]` as `probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

pub(crate) mod deciding {
    use crate::combinatorics::tie_probe::Tie;
    use std::sync::Mutex;
    /// `(test, tried, offered, decided, ties of this call when exhausted)`.
    pub(crate) type Row = (String, usize, usize, bool, Vec<(Tie, usize)>);
    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());
    pub(crate) fn record(tried: usize, offered: usize, decided: bool, ties: Vec<(Tie, usize)>) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push((
                std::thread::current().name().unwrap_or("?").to_string(),
                tried,
                offered,
                decided,
                ties,
            ));
    }
}
