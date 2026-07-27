//! Counters for **measuring** the judge ladder, off by default.
//!
//! These answer a question the crate cannot answer from the outside: *of the judgements that miss
//! the f64 filter and escalate, how many did not need the toleranced path at all?* A rigid rotation
//! preserves every determinant these judges take, so when a judgement's inputs all carry the **same
//! rotation chain**, its answer equals the answer on their pre-rotation rational bases — exactly,
//! with no tolerance. The share of escalations in that position decides whether routing them to the
//! exact predicate is worth building.
//!
//! **Why a runtime switch rather than `#[cfg(test)]`.** The corpus that produces real escalations
//! lives in `nacre-ops`, and a `cfg(test)` counter here is not compiled when *that* crate's tests
//! run. The switch costs one relaxed load per escalation, on a path that is about to do
//! hundreds of astro-float operations.

use crate::Pt3;
use core::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};

static ON: AtomicBool = AtomicBool::new(false);
static ESCALATED: AtomicU64 = AtomicU64::new(0);
/// Per-judge tallies, in `KINDS` order: (escalated, cancellable).
static BY_KIND: [(AtomicU64, AtomicU64); 5] = [
    (AtomicU64::new(0), AtomicU64::new(0)),
    (AtomicU64::new(0), AtomicU64::new(0)),
    (AtomicU64::new(0), AtomicU64::new(0)),
    (AtomicU64::new(0), AtomicU64::new(0)),
    (AtomicU64::new(0), AtomicU64::new(0)),
];
/// Judge names, indexed by the `kind` passed to [`escalation`].
pub const KINDS: [&str; 5] = [
    "orient3d",
    "indirect_orient3d",
    "cmp_coord",
    "dir_sign",
    "dir_orient3d",
];
static SAME_CHAIN: AtomicU64 = AtomicU64::new(0);
static SAME_CHAIN_EXACT_BASE: AtomicU64 = AtomicU64::new(0);

/// What the ladder did, since [`reset`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Counts {
    /// Judgements that missed the f64 filter and escalated.
    pub escalated: u64,
    /// …of those, every input carries the same rotation chain, so a rigid motion could be
    /// cancelled out.
    pub same_chain: u64,
    /// …of those, every base is also f64-representable, which is what the exact predicate needs
    /// (it takes `[f64; 3]`, and `Rat`'s `i128` cannot carry a 3×3 determinant).
    pub same_chain_exact_base: u64,
    /// `(escalated, cancellable)` per judge, indexed by [`KINDS`].
    pub by_kind: [(u64, u64); 5],
}

/// Start counting from zero. Not synchronized against concurrent judging — a measurement run
/// should be single-threaded.
pub fn reset() {
    ESCALATED.store(0, Relaxed);
    for (a, b) in &BY_KIND {
        a.store(0, Relaxed);
        b.store(0, Relaxed);
    }
    SAME_CHAIN.store(0, Relaxed);
    SAME_CHAIN_EXACT_BASE.store(0, Relaxed);
    ON.store(true, Relaxed);
}

/// Stop counting and read the tally.
pub fn take() -> Counts {
    ON.store(false, Relaxed);
    Counts {
        escalated: ESCALATED.load(Relaxed),
        same_chain: SAME_CHAIN.load(Relaxed),
        same_chain_exact_base: SAME_CHAIN_EXACT_BASE.load(Relaxed),
        by_kind: std::array::from_fn(|i| (BY_KIND[i].0.load(Relaxed), BY_KIND[i].1.load(Relaxed))),
    }
}

/// Is this base point representable in `f64` without loss? The exact predicate takes `f64`
/// coordinates, so a base that does not round-trip cannot be handed to it.
fn base_is_f64_exact(p: &Pt3) -> bool {
    p.base
        .iter()
        .all(|&r| nacre_scalar::Rat::try_from_f64(r.to_f64()) == Some(r))
}

/// Record one escalation and whether its inputs were in a position to avoid it.
///
/// Chains are compared **structurally** — same nodes, same order, same pivot. Two chains that are
/// different spellings of one motion (`30° + 30°` vs `60°`) read as different, which under-counts
/// rather than over-counts: a missed cancellation is slower, never wrong. Pivots must match too:
/// a rotation about `c` is `Rx + (c − Rc)`, and two different pivots leave two different
/// translations that do not cancel in a difference.
pub fn escalation(kind: usize, points: &[&Pt3]) {
    if !ON.load(Relaxed) {
        return;
    }
    ESCALATED.fetch_add(1, Relaxed);
    BY_KIND[kind].0.fetch_add(1, Relaxed);
    let Some(first) = points.first() else { return };
    let same = points.iter().all(|p| {
        p.chain.len() == first.chain.len()
            && p.chain
                .iter()
                .zip(first.chain.iter())
                .all(|(a, b)| a.axis == b.axis && a.angle == b.angle && a.point == b.point)
    });
    if same {
        SAME_CHAIN.fetch_add(1, Relaxed);
        if points.iter().all(|p| base_is_f64_exact(p)) {
            SAME_CHAIN_EXACT_BASE.fetch_add(1, Relaxed);
            BY_KIND[kind].1.fetch_add(1, Relaxed);
        }
    }
}
