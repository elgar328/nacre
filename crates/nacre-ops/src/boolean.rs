//! **The front door of the boolean.** Two live solids in, result solids out.
//!
//! This is the first stage, not the last: it validates, hands the operands to the arrangement
//! engine ([`crate::arrangement::boolean`]), and reports. The engine calls the **assembly**
//! ([`crate::assembly`]) when it has faces to build from, so the three read in the order they
//! run -- `boolean` -> `arrangement` -> `assembly` -- rather than the round trip that showed
//! while the front door and the assembly shared one module name.

use crate::assembly::check_result_topology;
#[cfg(test)]
use crate::assembly::tess_census;
use crate::draft::BoolKind;
use crate::{BoolError, reject};
use nacre_judge::Decision;
use nacre_judge::predicate::{Evidence, Notes, Site};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};
/// Boolean of two live solids (M5; the honest-reject strategy).
///
/// **Coverage:** planar solids. All three kinds go through the single per-plane-class arrangement
/// engine ([`crate::arrangement::boolean`]), which handles transverse, coplanar-contact,
/// coincident, contained and disjoint cases in one path and cleans its own output (coplanar-face
/// merge) so results are chainable. Anything it cannot resolve is rejected with [`BoolError`] —
/// never a silent wrong answer (DNA). Transactional: it computes the result in local structures
/// and pushes only after every degeneracy check passes, so a rejected boolean leaves the model
/// untouched.
pub fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    boolean_with_report(model, kind, a, b).map(|(solids, _)| solids)
}

/// [`boolean`], and **what the kernel had to assume to get there** — see [`BoolReport`].
///
/// A parallel entry point rather than a wider return type: the report is wanted by roughly one
/// caller in thirty, and changing `boolean`'s signature would rewrite every other one for nothing.
pub fn boolean_with_report(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(Vec<Handle<Solid>>, BoolReport), BoolError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids().contains(&a) || !model.live_solids().contains(&b) {
        return Err(BoolError::InputNotLive);
    }
    // The arrangement engine (`arrangement`) is the sole boolean path: one per-plane-class 2D
    // arrangement handles transverse, coplanar-contact, coincident, contained and disjoint cases,
    // and cleans its own output (coplanar-face merge) so results are chainable.
    // ★ **A rejected boolean leaves the live set as it found it** — every reject, not just the
    // topology one below. The engine retires the operands once the result is accepted, so a reject
    // raised after that point would otherwise hand back a model whose inputs had vanished. Restoring
    // here makes the rule hold no matter where inside the engine a reject is raised or added later.
    // (Orphaned result cells stay in the append-only arena, unreachable, as any superseded solid's
    // do — the arena is not restored and is not meant to be.)
    let snapshot = model.live_solids().to_vec();
    let (result, notes) = match crate::arrangement::boolean(model, kind, a, b) {
        Ok(v) => v,
        Err(e) => {
            model.restore_live(snapshot);
            return Err(surfacing(e));
        }
    };
    // Topological self-check on the assembled result (DNA: never return a malformed solid). Reject
    // rather than return. Valid results always pass, so this never false-rejects; the traversal
    // reads the topology stores directly (no adjacency rebuild, no coordinates).
    if let Some((t, at)) = check_result_topology(model, &result) {
        model.restore_live(snapshot);
        return Err(surfacing(match at {
            Some(w) => crate::reject_at(t, w),
            None => reject(t),
        }));
    }
    #[cfg(test)]
    tess_census::record(model);
    Ok((result, BoolReport::of(&notes)))
}

/// Record that `e` is leaving the kernel — the *surfaced* column of a test build's
/// [`crate::reject_census`]; a product build passes `e` through untouched.
///
/// Sits at the public entry points rather than at [`reject`] because those are two different
/// populations: guards are raised and swallowed (a retry tries the ring's next node), and only
/// what comes back here is something a caller ever sees. `InputNotLive` is not a
/// [`RejectReason`] — it is a caller mistake, not a guard — so it stays out of the census.
fn surfacing(e: BoolError) -> BoolError {
    #[cfg(any(test, feature = "test-util"))]
    if let BoolError::Rejected { reason, .. } = e {
        crate::reject_census::surfaced(reason);
    }
    e
}

/// **What the boolean had to take on faith**, so a user can see it and act on it.
///
/// Rotated geometry has no exact zero: a plane through turned points meets another at coordinates
/// no finite precision writes down, so "these two faces are the same plane" is *proved to within a
/// distance*, never proved outright. The kernel decides such a question by proving the separation
/// is below the coincidence limit — and then this says which questions those were and how close
/// the closest call came.
///
/// **It is a diagnosis, not a prompt.** Nothing here asks the user to choose; the choice was made
/// on evidence and the evidence is here. What it is *for* is noticing an unintended coincidence —
/// two features that met because a dimension made them meet — and going back to fix the design,
/// which is a thing only the author of the model can do.
#[derive(Clone, Debug, Default)]
pub struct BoolReport {
    /// **Faces merged into one plane class on toleranced evidence — first, because everything
    /// else follows from them.** A merge decided here changes which planes exist before a single
    /// vertex is computed, so a surprise in this list explains surprises everywhere else.
    pub merges: Vec<Evidence>,
    /// How many judgements were answered by a proved coincidence rather than a proved sign.
    pub coincidences: usize,
    /// The closest call: the coincidence with the **widest** bound, the one nearest to having
    /// been wrong. `None` when nothing was assumed at all — an axis-aligned model, typically,
    /// where every question has an exact answer.
    pub loosest: Option<Evidence>,
}

impl BoolReport {
    fn of(notes: &Notes) -> BoolReport {
        let all = notes.sorted();
        let widest = |e: &Evidence| match e.outcome {
            Decision::Coincident { within } => within.exp2(),
            _ => None,
        };
        BoolReport {
            merges: all
                .iter()
                .filter(|e| matches!(e.site, Site::PlanesCoplanar { .. }))
                .copied()
                .collect(),
            coincidences: all
                .iter()
                .filter(|e| matches!(e.outcome, Decision::Coincident { .. }))
                .count(),
            // Ties keep the first in sorted order, so the answer does not depend on the schedule.
            loosest: all
                .iter()
                .filter(|e| matches!(e.outcome, Decision::Coincident { .. }))
                .max_by_key(|e| widest(e))
                .copied(),
        }
    }
}
