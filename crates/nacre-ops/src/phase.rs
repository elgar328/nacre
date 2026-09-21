//! **Phase timers, on the production path.**
//!
//! ★ They live *inside* `boolean` rather than in a spike that replays it, because a replica measures
//! the proposition next to the one that matters — the last profile of these phases was taken with
//! `ClassReuse::Off` and read as if it were production's.
//!
//! ★ **One table, so it sits under everything that charges it.** Two stages write here — the
//! arrangement (`timed!`/`watch!`, and the scale counters in its split) and the plane setup
//! (`planes::standard::Watch`) — and `all()` renders every counter in one ordered report, which is
//! the spike test's whole output. A table split per stage could not produce that report, so the
//! table is one and lives below both rather than inside either.
//!
//! Read them with `--no-default-features`: parallel accumulation would sum CPU across threads, and a
//! share of a wall-clock whole computed from a CPU sum is not a share of anything. (That mistake
//! once made a part measure larger than its whole.) `cfg(test)` so release binaries carry nothing.

use std::sync::atomic::{AtomicU64, Ordering};

macro_rules! counters {
        ($($name:ident = $label:literal),+ $(,)?) => {
            $(pub(crate) static $name: AtomicU64 = AtomicU64::new(0);)+
            /// Every counter with its label, in report order.
            pub(crate) fn all() -> Vec<(&'static str, u64)> {
                vec![$(($label, $name.load(Ordering::Relaxed))),+]
            }
            pub(crate) fn reset() { $($name.store(0, Ordering::Relaxed);)+ }
        };
    }

counters! {
    SETUP      = "plane_index_setup",
    S_TRIPT3   = "    collect: tri_pt3 (chain replay)",
    S_COLLECT  = "    collect: the per-face loop",
    S_STD      = "    standard_for",
    S_EDGES    = "    edge_faces x2",
    S_CLASSES  = "    plane_classes (pairwise)",
    S_DENSE    = "    dense_planes + owners",
    TRACE_IN   = "trace_input (face table)",
    TRACE_ON   = "  trace_on_class",
    MERGE      = "  merge_coincident",
    SPLIT      = "  split_at_crossings",
    S_PART     = "    build the direction partition",
    COLLECT    = "    (1) collect crossings",
    SORT       = "    (2) sort + flush groups",
    COVER      = "    (3) cover sub-intervals",
    CELLS      = "  split + cells + nest + label + emit",
    C_SPLIT    = "    split_circles (ClassEdges::of)",
    C_EXTRACT  = "    walk_cells",
    E_ORDER    = "      per-half-edge order_along",
    E_ANGULAR  = "      angular_order",
    E_WALK     = "      the rest (walk + cells)",
    C_NEST     = "    nest_cells",
    C_LABEL    = "    label_cells",
    C_EMIT     = "    emit_faces",
    REUSE      = "  reuse pass-through",
    UNIFY      = "unify_coplanar_faces",
    SEAM       = "seam table + alias scan",
    ASSEMBLE   = "assemble_fuse_cut",
}

/// **Scale, not time.** The two ratios that decide whether a hull test over wall pairs or a
/// tighter covering loop is worth building: how many segments a wall carries (the hull test buys
/// nothing at 1), and how big the covering loop is against the collecting one.
pub(crate) mod scale {
    use std::sync::atomic::{AtomicU64, Ordering};
    pub(crate) static SEGS: AtomicU64 = AtomicU64::new(0);
    pub(crate) static WALLS: AtomicU64 = AtomicU64::new(0);
    pub(crate) static PTS: AtomicU64 = AtomicU64::new(0);
    /// `Σ_w |pts on w| × |segs on w|` — the covering loop's actual trip count.
    pub(crate) static COVER_TRIPS: AtomicU64 = AtomicU64::new(0);
    /// `Σ_w |segs|` — the collecting loop's actual trip count.
    pub(crate) static COLLECT_TRIPS: AtomicU64 = AtomicU64::new(0);
    pub(crate) fn add(c: &AtomicU64, n: usize) {
        c.fetch_add(n as u64, Ordering::Relaxed);
    }
    pub(crate) fn get(c: &AtomicU64) -> u64 {
        c.load(Ordering::Relaxed)
    }
    pub(crate) fn reset() {
        for c in [&SEGS, &WALLS, &PTS, &COVER_TRIPS, &COLLECT_TRIPS] {
            c.store(0, Ordering::Relaxed);
        }
    }
}

/// Run `f`, adding its wall time to `c`. Returns what `f` returned.
pub(crate) fn timed<R>(c: &AtomicU64, f: impl FnOnce() -> R) -> R {
    let t = std::time::Instant::now();
    let r = f();
    c.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
    r
}

/// Charges the enclosing scope to `c` **on drop** — for a block whose `?` or early `return`
/// would jump past the end of a closure.
pub(crate) struct Watch(&'static AtomicU64, std::time::Instant);

impl Watch {
    pub(crate) fn new(c: &'static AtomicU64) -> Self {
        Self(c, std::time::Instant::now())
    }
}

impl Drop for Watch {
    fn drop(&mut self) {
        self.0
            .fetch_add(self.1.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
}
