//! **Which guard rang, where — and which of those the caller actually saw.**
//!
//! Every [`crate::BoolError::Rejected`] in this crate is built by [`crate::reject`], so one
//! `#[track_caller]` there records the whole population: the reason, the source line of the guard
//! that raised it, and separately the reasons that actually left a public entry point.
//!
//! ★★★ **The two columns are not the same population, and reading one for the other is how this
//! measurement goes wrong.** Measured over the whole workspace suite while `point_in_component`
//! still raised `no_clear_ray` and its caller retried other nodes: 119 raises, 20 surfaced — 97
//! of them that one swallowed guard. A census keyed on raises alone would
//! have reported it as 82% of all rejects, which was false about everything a caller ever saw.
//! ★ **That population is gone on purpose** — the probe's "this node cannot
//! decide" is a typed abstention now, not an error to swallow — so the columns' gap reads the
//! other way: it is the count of *remaining* swallowed raises, near-zero is the goal state, and a
//! gap reopening is news of a new swallow.
//!
//! ★★ **What this does *not* answer: which *site* surfaced.** A reject is raised on a rayon worker
//! and returned from the main thread, so nothing cheap connects the two (a thread-local cannot
//! cross the boundary, so a tag-based census cannot see it). It is also not
//! needed: "is this guard swallowed by its caller?" is a **static** question, answered by reading
//! the caller, not by counting.
//!
//! ★ **Unconditional, like `nacre_judge::climb_census`** — and for its reason: the cost is a mutex
//! on a path that has already decided to fail, noise beside the work that decision ends. (A reader
//! in another crate is not the reason: a `test-util` gate reaches integration tests, as
//! `nacre_topo`'s push counters show — they sit on the ordinary push path, so they are gated.)
//!
//! # Reading it
//!
//! - **The standing gate** is `tests/reject_census.rs`: a frozen corpus of rejecting shapes whose
//!   `(reason, detail, file)` sets are pinned. Its population is those shapes and nothing else.
//! - **The whole-suite sweep** is the `reject-trace` feature, which prints each raise and surface
//!   as it happens (30-odd test binaries are 30-odd processes, so an in-memory table cannot
//!   outlive them):
//!
//!   ```text
//!   cargo test --workspace --features nacre-ops/reject-trace -- --nocapture 2>&1 \
//!     | grep '^RAISE ' | sort | uniq -c | sort -rn
//!   ```
//!
//! ★ **Raise *counts* are not stable; the *sets* are.** `debug` runs every traced boolean twice
//! (`arrangement`'s `#[cfg(debug_assertions)]` reference run), and `parallel` evaluates the
//! remaining plane classes of a failing input while `serial` stops at the first
//! (`crate::par::try_map_range`). Surfacing is unaffected: that same function returns the
//! **lowest-index** error by construction, so *which* reason comes back is schedule-independent.
//! Assert on sets and on distinct-site counts, never on raise counts.

use crate::RejectReason;
use std::collections::{BTreeMap, BTreeSet};
use std::panic::Location;
use std::sync::{Mutex, MutexGuard, PoisonError};

/// A reason as the census keys it: the stable identifier, plus the **categorical** part of its
/// payload.
///
/// ★★★ `detail` exists because [`RejectReason::TraceDeclined`]'s eleven
/// [`crate::DeclineKind`]s are all raised from **one line**. Keyed by `as_str()` alone they would
/// collapse into a single row — the census reproducing the very defect it is here to find (ten
/// causes reported under one label is what `TraceDeclined` was split out of). The rule: **the key
/// must be as fine as the reason's own vocabulary.** `PrecisionBudget`'s `needed`/`cap` are
/// measurements rather than a category, so they stay out of the key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReasonId {
    /// [`RejectReason::as_str`].
    pub reason: &'static str,
    /// [`crate::DeclineKind::as_str`] for `TraceDeclined`, `None` otherwise.
    pub detail: Option<&'static str>,
}

impl ReasonId {
    fn of(reason: RejectReason) -> Self {
        Self {
            reason: reason.as_str(),
            detail: match reason {
                RejectReason::TraceDeclined { kind, .. } => Some(kind.as_str()),
                _ => None,
            },
        }
    }
}

impl std::fmt::Display for ReasonId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.detail {
            Some(d) => write!(f, "{}({d})", self.reason),
            None => f.write_str(self.reason),
        }
    }
}

/// One place a guard built a rejection: the reason it named and the source line it named it at.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RaiseSite {
    pub id: ReasonId,
    pub file: &'static str,
    pub line: u32,
}

/// Guards that rang, and how often. Includes raises a caller swallowed.
static RAISED: Mutex<BTreeMap<RaiseSite, u64>> = Mutex::new(BTreeMap::new());
/// Reasons that left a public entry point in an `Err`.
static SURFACED: Mutex<BTreeMap<ReasonId, u64>> = Mutex::new(BTreeMap::new());

/// An instrument must never turn someone else's panic into a second one. Nothing here can poison
/// a lock (the guard is held across an increment and nothing else), but recovering costs a line.
fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Record a guard ringing. Called only from [`crate::reject`], which is `#[track_caller]`, so
/// `loc` is the guard's own line rather than `reject`'s.
pub(crate) fn raised(reason: RejectReason, loc: &'static Location<'static>) {
    let site = RaiseSite {
        id: ReasonId::of(reason),
        file: loc.file(),
        line: loc.line(),
    };
    trace_raise(reason, &site);
    *lock(&RAISED).entry(site).or_insert(0) += 1;
}

/// Record a rejection leaving the kernel. Called at the public boolean entry points — the
/// boundary the earlier one-shot census settled on, for the reason above.
pub(crate) fn surfaced(reason: RejectReason) {
    trace_surface(reason);
    *lock(&SURFACED).entry(ReasonId::of(reason)).or_insert(0) += 1;
}

/// What the guards did, **and reset** so the next window is its own.
///
/// ★ The reset is process-global: calling this steals the counts of anything running concurrently.
/// `nacre-ops`' unit-test binary raises three reasons of its own from tests that run in parallel,
/// so this is called from **one dedicated integration-test binary** (`tests/reject_census.rs`) and
/// nowhere else.
pub fn take() -> Census {
    Census {
        raised: std::mem::take(&mut *lock(&RAISED)).into_iter().collect(),
        surfaced: std::mem::take(&mut *lock(&SURFACED)).into_iter().collect(),
    }
}

/// One window of [`take`]. Sorted, because the maps are ordered.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Census {
    pub raised: Vec<(RaiseSite, u64)>,
    pub surfaced: Vec<(ReasonId, u64)>,
}

impl Census {
    /// Nothing rang and nothing surfaced.
    pub fn is_empty(&self) -> bool {
        self.raised.is_empty() && self.surfaced.is_empty()
    }

    /// **What the gate pins on the raise side**: which reason is raised in which file.
    ///
    /// The *line* is deliberately not part of it — a gate keyed on line numbers goes red on every
    /// unrelated edit above a guard, which is how a lock stops being read. Lines are printed by
    /// [`Self::report`] for diagnosis instead.
    pub fn raised_shape(&self) -> BTreeSet<(ReasonId, &'static str)> {
        self.raised.iter().map(|(s, _)| (s.id, s.file)).collect()
    }

    /// **What the gate pins on the surfaced side** — the strong half: this set is the same in
    /// every build (see the module doc on `try_map_range`).
    pub fn surfaced_shape(&self) -> BTreeSet<ReasonId> {
        self.surfaced.iter().map(|(id, _)| *id).collect()
    }

    /// How many **distinct sites** raised `reason` — the question "how many of a reason's guards
    /// actually fire", which raise counts cannot answer stably.
    pub fn distinct_sites(&self, reason: &str) -> usize {
        self.raised
            .iter()
            .filter(|(s, _)| s.id.reason == reason)
            .count()
    }

    /// The table, for `--nocapture`.
    pub fn report(&self) -> String {
        use std::fmt::Write;
        let mut out = String::new();
        let _ = writeln!(out, "raised ({} sites):", self.raised.len());
        for (s, n) in &self.raised {
            // Rendered first: a hand-written `Display` ignores the formatter's width, so
            // `{:<28}` on `s.id` directly would silently not pad.
            let id = s.id.to_string();
            let _ = writeln!(out, "  {n:>5}  {id:<26} {}:{}", s.file, s.line);
        }
        let _ = writeln!(out, "surfaced ({} reasons):", self.surfaced.len());
        for (id, n) in &self.surfaced {
            let _ = writeln!(out, "  {n:>5}  {id}");
        }
        out
    }
}

// The crate denies printing outside tests. This is the one exception, and it is the whole point of
// the `reject-trace` feature: an in-memory table dies with its process, and the population worth
// sweeping is spread over ~30 test binaries. Off by default, so a product build compiles no
// printing at all.
#[cfg(feature = "reject-trace")]
#[allow(clippy::print_stderr)]
fn trace_raise(reason: RejectReason, site: &RaiseSite) {
    // `Display` rather than the key, so `trace_declined(hole-ring)` reads as itself.
    eprintln!("RAISE {reason} {}:{}", site.file, site.line);
}

#[cfg(not(feature = "reject-trace"))]
fn trace_raise(_reason: RejectReason, _site: &RaiseSite) {}

#[cfg(feature = "reject-trace")]
#[allow(clippy::print_stderr)]
fn trace_surface(reason: RejectReason) {
    eprintln!("SURFACED {reason}");
}

#[cfg(not(feature = "reject-trace"))]
fn trace_surface(_reason: RejectReason) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DeclineKind;
    use nacre_math::Point3;
    use nacre_topo::Model;

    /// **The eleven [`DeclineKind`]s all ring from one line**, so `file:line` cannot tell them
    /// apart — `detail` is the only thing that can, and without it the census would report ten
    /// different causes as one row, which is the defect `TraceDeclined` was split out of.
    ///
    /// Asserted on the key function rather than on a fixture because `trace_declined` has never
    /// fired anywhere in the suite: there is no shape to run. Pure — it reads no static, so it is
    /// safe beside the parallel tests of this binary (which do ring guards of their own).
    #[test]
    fn the_key_is_as_fine_as_the_reasons_own_vocabulary() {
        let mut m = Model::new();
        let s = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0; 3]),
        );
        let face = m.shell(m.solid(s).outer).faces[0];
        let of = |kind| {
            ReasonId::of(RejectReason::TraceDeclined {
                kind,
                face: Some(face),
            })
        };

        let hole = of(DeclineKind::HoleRing);
        let outer = of(DeclineKind::OuterRing);
        assert_ne!(
            hole, outer,
            "two decline kinds collapsed into one census row"
        );
        assert_eq!(
            (hole.reason, outer.reason),
            ("trace_declined", "trace_declined"),
            "the stable identifier is shared — only `detail` separates them"
        );
        assert_eq!(hole.to_string(), "trace_declined(hole-ring)");

        // A reason whose payload is a measurement, not a category, keys on its name alone.
        assert_eq!(ReasonId::of(RejectReason::NoClearRay).detail, None);
        assert_eq!(
            ReasonId::of(RejectReason::PrecisionBudget {
                needed: 300,
                cap: 256
            })
            .detail,
            None
        );
    }
}
