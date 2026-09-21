//! Operations for the nacre kernel, plus a replayable operation log.
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::PadOnFace`]/[`Operation::PocketOnFace`] (M4) consume a prior op's face by
//! `Handle` (exposed via [`OpOutput`]) and supersede a solid (live-solid
//! semantics) — each is a tool prism plus a boolean, not a direct face-split. Ops are
//! applied by [`apply`] and folded by [`replay`]; every result is a **closed** solid, so
//! `nacre-validate` applies fully.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_math::Point3;
// The test modules read the root's names through their globs.
use nacre_exact::Rat;
#[cfg(test)]
use nacre_math::{Point2, Vector3};
use nacre_store::Handle;
use nacre_topo::{Face, HalfEdge, Model, Solid, Vertex};

mod arrangement;
mod assembly;
/// The cylinder-band pass — see the module docs.
mod bands;
mod boolean;
mod combinatorics;
mod construct;
mod draft;
mod error;
mod nesting;
mod ops;
mod par;
/// Phase timers — see the module docs. `cfg(test)`: a release binary carries nothing.
#[cfg(test)]
mod phase;
mod planes;
mod realize;
/// Public because the measurements that read it live in other crates — see the module doc.
pub mod reject_census;
mod reuse;
mod rotated_vertex;
mod sketch;
mod tolerant;
mod transform;

pub use boolean::{BoolReport, boolean, boolean_with_report};
pub use error::{BoolError, DeclineKind, RejectClass, RejectReason, RejectWhere};
pub(crate) use error::{reject, reject_at};
// The report's vocabulary: what a judgement was asked about, and what it established. Re-exported
// so a consumer reads one crate, not two.
pub use draft::BoolKind;
pub use nacre_geom::mixed::Edge2d;
pub use nacre_judge::Decision;
pub use nacre_judge::predicate::{Evidence, Site};
pub use ops::{
    DatumDef, LogCell, OpError, OpOutput, Operation, PlaneDef, Profile2d, ProfileRing, Ring2d,
    SketchFrame, SketchPlane, apply, face_plane, face_sketch_frame, frame_plane, replay,
};
pub use realize::{
    CacheDecline, Precision, RealizeError, Realized, RefineReport, realize_cache, realize_vertex,
    realize_vertex_decimal, refine_vertex_cache,
};
pub use sketch::{SketchError, arc_to_rat, arc_turns, arc_turns_rat, from_paths, from_rings};

/// Assert that `f` rejects *through the intended guard* — a reject test whose fixture drifts
/// onto a different guard then fails instead of silently passing. Compares the projected
/// reason only: the location payload is a measurement, asserted (approximately) by the
/// per-reason payload tests, not here.
#[cfg(test)]
fn assert_rejects<T: std::fmt::Debug + PartialEq>(
    f: impl FnOnce() -> Result<T, BoolError>,
    expect: RejectReason,
) {
    match f() {
        Err(BoolError::Rejected { reason, .. }) => assert_eq!(reason, expect),
        other => panic!("expected a reject with {expect:?}, got {other:?}"),
    }
}

/// The start vertex of a half-edge (`vertices[0]` if forward, else `vertices[1]`).
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    model.he_start(he)
}

use std::collections::HashMap;

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

// The lib's own tests share fixture files with the integration tests, which name this crate
// `nacre_ops`; under `cfg(test)` the crate answers to that name too.
#[cfg(test)]
extern crate self as nacre_ops;

#[cfg(test)]
#[path = "tests/suite/mod.rs"]
pub mod tests;

/// **Ledger totals, printed** — a measurement, not an assertion. Run last, single-threaded, so
/// the sums cover the whole lib suite:
/// `cargo test -p nacre-ops --lib -- --include-ignored --test-threads=1 --nocapture zzz_ledger`
/// is *not* it (a filter runs only this); the whole-suite run is
/// `cargo test -p nacre-ops --lib -- --include-ignored --test-threads=1 --nocapture --skip stress --skip spike --skip direction_families`.
#[cfg(test)]
#[path = "tests/zzz_ledger.rs"]
mod zzz_ledger;
