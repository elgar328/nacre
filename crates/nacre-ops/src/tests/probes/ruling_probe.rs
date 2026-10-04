//! **What a holed lateral's ruling actually came out as** — the lock for
//! `a_holed_laterals_ruling_grazes_where_the_hole_is`.
//!
//! ★ This cell's result lives *inside* an operation that still refuses further down, so there is no
//! solid to open and count faces on. The trace's own answer is the thing to hold, the way
//! `arrangement`'s `sides == [-1, 1]` lock already holds the hole-free one.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `ruling_probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use super::SegKind;
use crate::ledger::Ledger;

/// The sentence the ruling label's postcondition panics with — one spelling, so the
/// commuting oracle's `KNOWN` list names the site by the same constant the
/// `assert` prints.
pub(crate) const WRONG_SIDE: &str = "the ruling label took the wrong side of the wall";

/// One entry per ruling a cycle grazed, in emission order: the kinds of its pieces, beside
/// the cylinder's origin and the ruling's first and last station — so a reader can tell one of
/// its own fixtures from another (a panel's ruling is one graze, a chain's a transversal then a
/// graze).
pub(crate) static CARVED: Ledger<Carved> = Ledger::new();

#[derive(Clone, Debug)]
pub(crate) struct Carved {
    pub origin: [f64; 3],
    pub span: [f64; 2],
    pub kinds: Vec<SegKind>,
}

/// **The ruling label's postcondition, one entry per ruling piece**.
///
/// ★★★★★ **A second, independent description of the side the derivation picked.**
/// [`super::ruling_interior_is_even`] derives it from `side · κ · frame_sign`; the check asks the *content*
/// instead — on the side of the surface the lateral's material lies (inside for a boss,
/// outside for a bore or a notch: `MergedRuling::orient`) its own solid has material, and on
/// the other side it does not.
/// The two share no step, so a disagreement is real. (`ArcLabels`' doc set its own side by
/// this same content rule, measured.)
///
/// `Some(true)` the two agree · `Some(false)` they contradict · `None` the content does not
/// distinguish, so the check is **blind** there and only the derivation speaks. A boss whose
/// own plate surrounds it is exactly such a case, which is why this is recorded rather than
/// asserted.
pub(crate) static SIDE_CHECK: Ledger<Option<bool>> = Ledger::new();

/// One entry per ruling piece: whether it got a label at all (a class without a world sense declines).
pub(crate) static LABELLED: Ledger<bool> = Ledger::new();

/// One entry per graze: the stated `body_above`, beside the **realized** stored normal of the
/// wall class it is stated against.
///
/// ★★★★★ **The sign has no oracle downstream yet** — the operation this exercises refuses at
/// the labelling for an unrelated incompleteness, so flipping any factor of `body_above`
/// changes nothing a test can see (☑ measured: all four controls green). So the claim is
/// checked against the fixture's own geometry instead: the buried half of the boss lies at
/// `y < 0` of the wall, so the face is on the stored-normal side exactly when that normal
/// points at `−y`.
/// Beside them the cylinder's origin and the run's stations, the fixture's identity.
pub(crate) static GRAZE_SIDE: Ledger<GrazeSide> = Ledger::new();

#[derive(Clone, Copy, Debug)]
pub(crate) struct GrazeSide {
    /// Whether the cycle whose run this is is an inner loop — a hole's face lies on one side of
    /// the wall, a panel's on the other, so a reader must say which it is asking about.
    pub inner: bool,
    pub body_above: bool,
    /// The realized stored normal's `y` component of the wall class.
    pub ny: f64,
    pub origin: [f64; 3],
    pub span: [f64; 2],
}
