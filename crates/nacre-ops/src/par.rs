//! **The one place the parallel switch lives.**
//!
//! The boolean's expensive phases are independent per plane class, and evaluating them on
//! several cores is the only lever that touches every judgement at once. What makes that
//! safe is not the absence of shared state — there is none to speak of — but **order**:
//! `assemble_fuse_cut` assigns vertex handles by first appearance across the faces it is
//! given, so `Store::push` order *is* handle identity, and handle identity is what replay
//! determinism rests on (design §2). A result that depended on which worker finished first
//! would not be wrong once; it would be a different model every run.
//!
//! So every helper here **collects in index order**. Nothing reduces across items, and no
//! item's computation is itself split, so no floating-point or `BigFloat` sum is ever
//! reassociated. A caller reads the same code in both builds and cannot reintroduce the
//! hazard by accident, because there is no call site where the ordering decision is made.
//!
//! The `parallel` feature (default on) is also the switch for `Pt3`'s high-precision cache
//! (`Arc<OnceLock>` vs `Rc<OnceCell>`) — which is why the two bodies below need different
//! bounds, and why **both feature combinations have to be built**. `cargo` will not check
//! the one you are not using.

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

    /// Nothing to do is not an error, and not a panic — `0..0` is a real case here (a
    /// boolean whose operands share no plane class at all).
    #[test]
    fn an_empty_range_is_fine() {
        let r: Result<Vec<usize>, ()> = try_map_range(0, Ok);
        assert!(r.unwrap().is_empty());
    }
}
