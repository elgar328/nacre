//! **The one place the parallel switch lives.**
//!
//! The boolean's expensive phases are independent per plane class, and evaluating them on
//! several cores is the only lever that touches every judgement at once. What makes that
//! safe is not the absence of shared state — there is none to speak of — but **order**:
//! `assemble_fuse_cut` assigns vertex handles by first appearance across the faces it is
//! given, so `Store::push` order *is* handle identity, and handle identity is what replay
//! determinism rests on. A result that depended on which worker finished first
//! would not be wrong once; it would be a different model every run.
//!
//! So every helper here **collects in index order**. Nothing reduces across items, and no
//! item's computation is itself split, so no floating-point or `BigFloat` sum is ever
//! reassociated. A caller reads the same code in both builds and cannot reintroduce the
//! hazard by accident, because there is no call site where the ordering decision is made.
//!
//! The `parallel` feature (default on) is also the switch for `WitnessPoint`'s high-precision cache
//! (`Arc<OnceLock>` vs `Rc<OnceCell>`) — which is why the two bodies below need different
//! bounds, and why **both feature combinations have to be built**. `cargo` will not check
//! the one you are not using.

/// Below this many items, run on one thread. **Measured, not chosen** — and only for
/// [`map_range`], whose items are single-point realizations, where the smallest models
/// (a dozen faces, ~36 points) were paying more to dispatch than to compute.
///
/// The crossover, 25-fin fold against 200 small axis-aligned fuses:
///
/// ```text
///   threshold      25 fins    small boolean
///        none      176.6ms          575.6µs   <- the small case regressed
///          64      177.8ms          530.7µs
///         128      184.6ms          514.8µs
///         256      214.4ms          522.6µs
/// ```
///
/// 64 removes the regression without costing the large case anything; past it the fold's
/// mid-sized booleans start falling back to one thread and the large case pays for it.
///
/// It is a count where the quantity that matters is *work*, so a small model made entirely
/// of deeply-rotated points would go serial when it need not have. That case is bounded by
/// being small: it is microseconds either way. Changing the answer is not among the risks —
/// both branches compute the same values in the same order.
#[cfg(feature = "parallel")]
const PARALLEL_FLOOR: usize = 64;

/// Map `f` over `0..n`, collecting **in index order**.
#[cfg(feature = "parallel")]
pub(crate) fn map_range<R: Send>(n: usize, f: impl Fn(usize) -> R + Sync + Send) -> Vec<R> {
    if n < PARALLEL_FLOOR {
        return (0..n).map(f).collect();
    }
    use rayon::prelude::*;
    // `Range<usize>` is an `IndexedParallelIterator`, so `collect` restores index order
    // however the work was scheduled.
    (0..n).into_par_iter().map(f).collect()
}

/// Map `f` over `0..n`, collecting **in index order**.
#[cfg(not(feature = "parallel"))]
pub(crate) fn map_range<R>(n: usize, f: impl Fn(usize) -> R) -> Vec<R> {
    (0..n).map(f).collect()
}

/// Map a fallible `f` over `0..n`, collecting **in index order**, and return the
/// **lowest-index** error if any item fails.
///
/// The error rule is the subtle part, and the reason this exists rather than being written out at
/// each call site: the sequential loops it replaces return at the first failing class, so
/// anything else would change *which* rejection a model reports. **Never use rayon's
/// `collect::<Result<_, _>>()` here** — it short-circuits, and which error survives then
/// depends on the schedule.
///
/// The two builds differ in how much work they do on a failing input, not in the answer:
/// the serial one stops at the first error (exactly as before), while the parallel one has
/// already evaluated the later items. That matters only for cost, and for the fact that a
/// latent panic in a later item becomes reachable — which the reject corpus covers.
#[cfg(feature = "parallel")]
pub(crate) fn try_map_range<R: Send, E: Send>(
    n: usize,
    f: impl Fn(usize) -> Result<R, E> + Sync + Send,
) -> Result<Vec<R>, E> {
    use rayon::prelude::*;
    let all: Vec<Result<R, E>> = (0..n).into_par_iter().map(f).collect();
    // A *sequential* collect over an ordered vec: the first `Err` by index wins.
    all.into_iter().collect()
}

/// Map a fallible `f` over `0..n`, returning the **lowest-index** error.
#[cfg(not(feature = "parallel"))]
pub(crate) fn try_map_range<R, E>(
    n: usize,
    f: impl Fn(usize) -> Result<R, E>,
) -> Result<Vec<R>, E> {
    (0..n).map(f).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Order, under both builds. A shuffled workload would show up here as a permutation.
    #[test]
    fn results_arrive_in_index_order() {
        let v: Result<Vec<usize>, ()> = try_map_range(1000, |i| Ok(i * 7));
        assert_eq!(v.unwrap(), (0..1000).map(|i| i * 7).collect::<Vec<_>>());
    }

    /// **The proposition the sequential loops relied on**: the error a caller sees is the
    /// one belonging to the lowest index, not whichever finished first. Two failing items
    /// far apart, so a schedule that got there in the other order would be visible.
    #[test]
    fn the_lowest_index_error_is_the_one_returned() {
        let r: Result<Vec<usize>, usize> =
            try_map_range(500, |i| if i == 17 || i == 400 { Err(i) } else { Ok(i) });
        assert_eq!(r, Err(17));
        // And with no failure, the values still come back in order.
        let ok: Result<Vec<usize>, ()> = try_map_range(500, Ok);
        assert_eq!(ok.unwrap(), (0..500).collect::<Vec<_>>());
    }

    /// **The floor is a cost knob, never an answer knob.** Sizes on both sides of it must
    /// give the same values in the same order, or a model would compute one thing when small
    /// and another when large.
    #[test]
    fn both_sides_of_the_parallel_floor_agree() {
        let f = |i: usize| (i * 2654435761) % 1_000_003;
        for n in [1, 63, 64, 65, 200] {
            assert_eq!(
                map_range(n, f),
                (0..n).map(f).collect::<Vec<_>>(),
                "n = {n}"
            );
        }
    }

    /// **rayon is named here and nowhere else in the crate.**
    ///
    /// The serial build is what the wasm playground takes, and **nothing builds it**: the
    /// pre-commit hook runs the default features, and a workspace-wide
    /// `--no-default-features` unifies `parallel` back on through `nacre-oracle`'s and
    /// `nacre-props`'s dev-dependencies. So it has to be checked deliberately
    /// (`cargo test -p nacre-ops --no-default-features`), and between checks the thing that
    /// silently breaks it is a second place where parallelism lives — a `par_iter` reached
    /// for at a call site, without the `#[cfg]` pair that keeps the serial build compiling.
    ///
    /// This does not replace running that command; it catches the drift that makes running
    /// it necessary, at no cost. Same shape as the facade's test that *reads* `Cargo.toml`
    /// rather than building a feature combination it cannot build.
    ///
    /// What is forbidden is the **parallel-iteration decision**, not the word rayon: the
    /// thread-order test builds a one-thread pool on purpose, and that is the opposite of a
    /// hazard. So this looks for the call syntax, which is also why it does not trip over
    /// prose — it caught this file's own doc comment when it matched the bare name.
    ///
    /// Being textual, it can be talked around (write `.par_iter()` in a comment and it
    /// fires; write the call some other way and it does not). It is a tripwire on the
    /// ordinary way to get this wrong, not a proof.
    #[test]
    fn the_parallel_switch_lives_only_here() {
        let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        // Every file under `src`, folders included: a module that became a folder is still
        // this crate's code, and a scan that stopped at the top level would go quiet about it.
        let mut files = Vec::new();
        let mut dirs = vec![src];
        while let Some(dir) = dirs.pop() {
            for entry in std::fs::read_dir(&dir).expect("src") {
                let path = entry.expect("entry").path();
                if path.is_dir() {
                    dirs.push(path);
                } else {
                    files.push(path);
                }
            }
        }
        let mut offenders = Vec::new();
        for path in files {
            if path.extension().is_some_and(|e| e == "rs")
                && path.file_name().is_some_and(|f| f != "par.rs")
                && {
                    let text = std::fs::read_to_string(&path).expect("read");
                    [
                        ".par_iter(",
                        ".into_par_iter(",
                        ".par_chunks(",
                        ".par_bridge(",
                        "use rayon::prelude",
                    ]
                    .iter()
                    .any(|api| text.contains(api))
                }
            {
                offenders.push(path);
            }
        }
        assert!(
            offenders.is_empty(),
            "a parallel iterator outside par.rs, so the serial build has a second way to \
             break -- and a second place where index order could be lost: {offenders:?}"
        );
    }

    /// Nothing to do is not an error, and not a panic — `0..0` is a real case here (a
    /// boolean whose operands share no plane class at all).
    #[test]
    fn an_empty_range_is_fine() {
        let r: Result<Vec<usize>, ()> = try_map_range(0, Ok);
        assert!(r.unwrap().is_empty());
    }
}
