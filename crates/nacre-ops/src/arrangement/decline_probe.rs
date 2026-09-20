//! The ⊥ road's answers, recorded: one entry per lateral face × ⊥ class it meets — the station,
//! the outer answer, what kinds of cycles the face has, and how many extents were carved and how
//! many spans came out. What says the road ran on a panel or a chain, and what it said there.
//! **Where an `OuterRing` decline was produced**. Three sites push that one kind and
//! the census key stops at the kind, so a fixture could say *that* it was declined but not by
//! which sentence. Test-only tally by producing site, keyed by thread name like `tie_probe` — a
//! test reads its own rows and never a parallel test's.

use std::sync::Mutex;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Site {
    /// The transversal walk named no outer loop for the face.
    NoOuter,
    /// A seated face's polygon ring collapsed in `plane_ring` — a carrier it could not spell.
    PolyCollapsed,
    /// A seated face whose outer loop is a lateral rim, or unnamed.
    RimOrNone,
}

pub(crate) static ROWS: Mutex<Vec<(String, Site)>> = Mutex::new(Vec::new());

pub(crate) fn mark(site: Site) {
    let name = std::thread::current().name().unwrap_or("?").to_string();
    ROWS.lock()
        .expect("the probe's lock is never held across a panic")
        .push((name, site));
}

/// The sites produced on threads whose name contains `tag`, in order. The trace runs on
/// rayon workers under `parallel`, so a test attributes rows by running its boolean in a pool
/// named after itself; a sequential build's test thread carries the test's name already.
pub(crate) fn named_like(tag: &str) -> Vec<Site> {
    ROWS.lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        .filter(|(n, _)| n.contains(tag))
        .map(|(_, s)| *s)
        .collect()
}
