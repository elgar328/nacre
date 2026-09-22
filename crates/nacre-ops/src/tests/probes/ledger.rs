//! **Every probe row says whose work made it**, so a test can count its own.
//!
//! The probes are process-global: one list per question, written by every test in this binary at
//! once. A reader that takes the whole list, or a window of it, sees other tests' rows too — which
//! is why a count could not be asserted, only existence or a universal. A row carries its
//! [`Owner`] instead: the test the work belongs to. [`Ledger::all`] is still the binary's, for the
//! universals and the whole-suite ledger; [`Ledger::mine`] is this test's alone, and a count over
//! it is a fact about the fixture.
//!
//! **The owner is the thread's name, up to `#`.** A boolean's expensive phases run on rayon
//! workers (`crate::par`), so the pushing thread is usually not the test's. [`owned`] runs the
//! test's work in a pool whose workers are named `{test}#{i}`, and every row they push is the
//! test's. Work done on the test thread itself is already named after the test, by libtest, in
//! parallel and `--test-threads=1` runs alike.
//!
//! ★ A test that reads `mine()` must build its population inside [`owned`]: in the parallel build a
//! boolean called outside it hands the work to the **global** pool, whose workers are nameless, and
//! the rows go to `"?"`. That shows up as a count of zero, not as a quiet pass.
//!
//! ★ **An instrument, so the file lives in the test tree**; `crate::ledger` mounts it with
//! `#[path]`.

use std::sync::{Arc, Mutex, MutexGuard};

/// Whose work a row belongs to — a test's name, shared by its pool's workers.
pub(crate) type Owner = Arc<str>;

/// This thread's owner, made once per thread: pushing a row is then a refcount bump, not a
/// string copy, which matters because the biggest ledgers take thousands of rows per suite.
fn make_owner() -> Owner {
    let t = std::thread::current();
    let name = t.name().unwrap_or("?");
    // `{test}#{i}` is one of `owned`'s workers; a test's own thread has no `#`.
    Arc::from(name.split('#').next().unwrap_or(name))
}

thread_local! {
    static OWNER: Owner = make_owner();
}

/// The owner of the current thread's rows.
pub(crate) fn owner() -> Owner {
    OWNER.with(Clone::clone)
}

/// One probe's rows, each with its owner.
pub(crate) struct Ledger<T>(Mutex<Vec<(Owner, T)>>);

impl<T> Ledger<T> {
    pub(crate) const fn new() -> Self {
        Ledger(Mutex::new(Vec::new()))
    }

    /// A reader that panics while holding this lock (an assertion inside the loop) poisons it;
    /// the next test still gets its rows, because a probe's list has no invariant a panic could
    /// have left half-written.
    fn rows(&self) -> MutexGuard<'_, Vec<(Owner, T)>> {
        self.0.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Record a row for this thread's owner.
    pub(crate) fn push(&self, row: T) {
        let who = owner();
        self.rows().push((who, row));
    }

    /// How many rows this test's work made.
    pub(crate) fn mine_len(&self) -> usize {
        let me = owner();
        self.rows().iter().filter(|(who, _)| *who == me).count()
    }

    /// How many rows the whole binary has made.
    pub(crate) fn len(&self) -> usize {
        self.rows().len()
    }
}

impl<T: Clone> Ledger<T> {
    /// Every row in the binary — what a universal assertion and the whole-suite ledger read.
    pub(crate) fn all(&self) -> Vec<T> {
        self.rows().iter().map(|(_, r)| r.clone()).collect()
    }

    /// The rows this test's own work made.
    pub(crate) fn mine(&self) -> Vec<T> {
        self.owned_by(&owner())
    }

    /// Every row with its owner — for the whole-suite ledger, which names the test a row came
    /// from.
    pub(crate) fn all_owned(&self) -> Vec<(Owner, T)> {
        self.rows().clone()
    }

    /// This test's rows recorded since `from`, a snapshot of [`len`](Ledger::len) — for a reader
    /// that wants one stretch of its own work rather than all of it.
    pub(crate) fn mine_since(&self, from: usize) -> Vec<T> {
        let me = owner();
        self.rows()
            .iter()
            .skip(from)
            .filter(|(who, _)| *who == me)
            .map(|(_, r)| r.clone())
            .collect()
    }

    /// The rows owned by `who` — [`mine`](Ledger::mine) names the caller. Equality, never
    /// `contains`: one test's name is a substring of another's often enough (`a`, `a_b`).
    pub(crate) fn owned_by(&self, who: &str) -> Vec<T> {
        self.rows()
            .iter()
            .filter(|(owner, _)| &**owner == who)
            .map(|(_, r)| r.clone())
            .collect()
    }
}

/// Run `f` so that every probe row it makes is **this test's**.
///
/// Under `parallel` that means a pool of its own whose workers carry the test's name; the pool's
/// size is rayon's default, so a fixture that runs a whole corpus is not squeezed onto fewer cores
/// than the suite gives it. The serial build has no workers to name and runs `f` where it stands.
#[cfg(feature = "parallel")]
pub(crate) fn owned<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    let me = {
        let t = std::thread::current();
        t.name().unwrap_or("?").to_string()
    };
    rayon::ThreadPoolBuilder::new()
        .thread_name(move |i| format!("{me}#{i}"))
        .build()
        .expect("a probe pool")
        .install(f)
}

/// The serial build's [`owned`]: `f` runs on the test's own thread, which libtest has already
/// named after the test. **No `Send` bound** — the serial judge's caches are `Rc`-backed, so a
/// bound copied from the parallel arm would refuse code that never leaves this thread.
#[cfg(not(feature = "parallel"))]
pub(crate) fn owned<R>(f: impl FnOnce() -> R) -> R {
    struct Mark;
    impl Drop for Mark {
        fn drop(&mut self) {
            OWNED.with(|c| c.set(false));
        }
    }
    OWNED.with(|c| c.set(true));
    let _mark = Mark;
    f()
}

#[cfg(not(feature = "parallel"))]
thread_local! {
    static OWNED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether the work running here belongs to a test that asked for its rows — the switch for a
/// probe too expensive to record unconditionally (`nesting::nesting_probe`).
///
/// The two builds answer differently because [`owned`] does: a pool worker says so by its name,
/// and the serial arm, which has no worker, leaves a mark on the thread it runs on.
#[cfg(feature = "parallel")]
pub(crate) fn owned_thread() -> bool {
    std::thread::current()
        .name()
        .is_some_and(|n| n.contains('#'))
}

#[cfg(not(feature = "parallel"))]
pub(crate) fn owned_thread() -> bool {
    OWNED.with(std::cell::Cell::get)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **A row another test's thread pushed is not mine, and a row my pool's worker pushed is.**
    ///
    /// Both halves matter: the first is the contamination this type exists to stop, and the second
    /// is the whole reason a test may count at all — its rows are made on threads it never names
    /// itself.
    #[test]
    fn a_row_belongs_to_the_test_whose_thread_made_it() {
        static L: Ledger<u32> = Ledger::new();
        let me = owner().to_string();
        L.push(1);
        // A stranger, and a stranger's worker.
        for (name, row) in [("someone_else", 2u32), ("someone_else#3", 3)] {
            std::thread::Builder::new()
                .name(name.to_string())
                .spawn(move || L.push(row))
                .expect("a probe thread")
                .join()
                .expect("the probe thread");
        }
        // My own worker, named the way `owned` names them.
        std::thread::Builder::new()
            .name(format!("{me}#1"))
            .spawn(|| L.push(4))
            .expect("a probe thread")
            .join()
            .expect("the probe thread");
        assert_eq!(L.mine(), vec![1, 4], "rows: {:?}", L.all());
        assert_eq!(L.mine_len(), 2);
        assert_eq!(L.len(), 4);
        assert_eq!(L.owned_by("someone_else"), vec![2, 3]);
    }

    /// **One owner's name may be a prefix of another's**, which is how a `contains` filter would
    /// hand `a` the rows of `a_b`. ☑ Rewriting the filter as `contains` turns this red.
    #[test]
    fn an_owner_that_is_a_prefix_of_another_does_not_take_its_rows() {
        static L: Ledger<u32> = Ledger::new();
        for (name, row) in [("a", 1u32), ("a_b", 2), ("a#0", 3)] {
            std::thread::Builder::new()
                .name(name.to_string())
                .spawn(move || L.push(row))
                .expect("a probe thread")
                .join()
                .expect("the probe thread");
        }
        assert_eq!(L.owned_by("a"), vec![1, 3]);
        assert_eq!(L.owned_by("a_b"), vec![2]);
    }

    /// **`owned` puts the work where the rows come back as the test's.** The body runs on a pool
    /// worker in the parallel build and on the test's thread in the serial one; either way it says
    /// it is owned, and what it pushes is `mine`.
    #[test]
    fn owned_work_pushes_rows_this_test_can_count() {
        static L: Ledger<u32> = Ledger::new();
        assert!(!owned_thread(), "a test's own thread is not owned work");
        let n = owned(|| {
            assert!(owned_thread(), "the body of `owned` is owned work");
            L.push(7);
            L.mine_len()
        });
        assert_eq!(n, 1);
        assert_eq!(L.mine(), vec![7]);
        assert!(!owned_thread(), "the mark does not outlive the call");
    }
}
