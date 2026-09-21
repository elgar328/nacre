//! **Does the ring's lexicographic minimum node really support the ring?**
//! `loop_winding` reads the turn there and its doc argues the node is a hull vertex — true for a
//! polygon, and an open question the moment an edge is an arc. Counts only.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::combinatorics`] mounts it with `#[path]` as `hull_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use std::sync::Mutex;

/// `(rings, arc_rings, circle_below_lo, undecided, PREMISE_BROKEN)` — the last is the
/// true count: `E` is lexicographically below `lo` **and** in an arc's interior.
pub(crate) static ROWS: Mutex<(usize, usize, usize, usize, usize)> = Mutex::new((0, 0, 0, 0, 0));

/// Arcs whose circle's minimum is irrational — the axis is not ⊥ to the first world axis the
/// circle spans — so this instrument says nothing about them.
pub(crate) static TILTED: Mutex<usize> = Mutex::new(0);

pub(crate) fn note(arcs: usize, below: usize, undecided: usize, inside: usize, tilted: usize) {
    *TILTED
        .lock()
        .expect("the probe's lock is never held across a panic") += tilted;
    let mut g = ROWS
        .lock()
        .expect("the probe's lock is never held across a panic");
    g.0 += 1;
    g.1 += usize::from(arcs > 0);
    g.2 += usize::from(below > 0);
    g.3 += undecided;
    g.4 += usize::from(inside > 0);
}
