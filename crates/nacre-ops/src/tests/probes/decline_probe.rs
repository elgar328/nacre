//! The ⊥ road's answers, recorded: one entry per lateral face × ⊥ class it meets — the station,
//! the outer answer, what kinds of cycles the face has, and how many extents were carved and how
//! many spans came out. What says the road ran on a panel or a chain, and what it said there.
//! **Where an `OuterRing` decline was produced**. Three sites push that one kind and
//! the census key stops at the kind, so a fixture could say *that* it was declined but not by
//! which sentence. Test-only tally by producing site: the rows carry their owner, so a test reads
//! its own and never a parallel test's.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `decline_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use crate::ledger::Ledger;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Site {
    /// The transversal walk named no outer loop for the face.
    NoOuter,
    /// A seated face's polygon ring collapsed in `plane_ring` — a carrier it could not spell.
    PolyCollapsed,
    /// A seated face whose outer loop is a lateral rim, or unnamed.
    RimOrNone,
}

pub(crate) static ROWS: Ledger<Site> = Ledger::new();

pub(crate) fn mark(site: Site) {
    ROWS.push(site);
}
