//! Combinatorial queries over the per-face planar arrangement.
//!
//! This module answers combinatorial questions about how one solid's boundary cuts
//! a face of the other. It lives in `nacre-ops` and not in `nacre-geom` because it
//! needs `Model`/`Face`/`Edge`, and geom sits below topo (dependencies
//! flow upward only).
//!
//! Everything decided here is decided by an exact predicate. Coordinates that
//! appear (`three_planes`' cache) are never the basis of a decision — the truth of
//! a seam point is its plane triple, as it is for a measured vertex.
//!
//! # One `usize`, two meanings — now two tables
//!
//! An arrangement reasons **per plane**, but a solid gives you **faces**: an earlier boolean can
//! split one geometric plane between two faces (a base's exposed top and the cantilever underside
//! above it) whose outward normals **oppose**. Both were rows of one `planes` table, so the same
//! `usize` meant "face" here and "plane" there, and when the two meanings met in one comparison the
//! answer was wrong *silently* — four times on this branch, most recently as 117748 predicate calls
//! that read "different plane" for two faces of one plane and answered from rounding noise.
//!
//! That used to be held by a naming convention (`fp` / `fc`) and a debug-time net. It is now the
//! type: the predicates here take [`crate::planes::WorkingPlane`], which has no face geometry to offer, and
//! `plane_ix` is the one place a face index becomes a plane index (in [`loop_triples`]).
//!
//! **Exception:** the code that *defines* the classes (`crate::fill_classes` →
//! `crate::shares_or_coplanar` → `Judge::planes_coplanar`) necessarily runs before a
//! plane table exists, so it takes face indices — hence that predicate's generic `Witness` bound.

use crate::planes::{ClassIx, WorkingPlane, edge_incidence};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Model, Solid, Vertex};
use std::collections::HashMap;

mod component;
mod direction;
#[cfg(test)]
pub(crate) mod hull_probe;
mod in_faces;
mod incidence_order;
mod loops;
mod mixed_ring;
mod names;
mod rat;
mod rays;
mod ring_walk;
#[cfg(test)]
#[path = "../tests/combinatorics.rs"]
mod tests;
#[cfg(test)]
pub(crate) mod tie_probe;
mod winding;
mod witnesses;

pub(crate) use component::*;
pub(crate) use direction::*;
pub(crate) use in_faces::*;
pub(crate) use incidence_order::*;
pub(crate) use loops::*;
pub(crate) use mixed_ring::*;
pub(crate) use names::*;
pub(crate) use rat::*;
pub(crate) use rays::*;
pub(crate) use ring_walk::*;
pub(crate) use winding::*;
pub(crate) use witnesses::*;

/// **A direction on a working plane**: the carrier that supplies it, and which way along it.
///
/// ★★★ **The fields are private, and that is the whole point.** Every rule that reads this
/// representation — [`turn`], [`antiparallel`], [`parallel_carriers`] via [`continuation`] — lives
/// in this module, and Rust's module privacy is what keeps it that way: a sibling module cannot
/// reach in even by accident. The node identity needs a `rg` gate (it caught fifteen bypasses)
/// because it is an enum whose variants are nameable everywhere; here
/// the boundary is a compile error instead of a grep.
///
/// ★★ **A test that must read a direction is not a violation — it is an oracle.** It should read
/// **its own** inputs (the `(carrier, sense)` it built) rather than this type, both because that
/// keeps the boundary and because an oracle read back out of the value under test is an oracle
/// derived from its own answer. Do not add accessors for it.
///
/// ★★★ **`carrier`, not `wall` — and it grew when arcs arrived.** A straight edge's direction is a plane carrier and a sign;
/// an arc's is a **tangent at one end**, and the two ends differ. So the type is a sum, and the arc
/// arm carries *everything the comparison needs* — the end as `(line, s)`, the circle's centre, and
/// which way the axis points against the class's stored normal — rather than an index a reader
/// would have to resolve against a table. That is what keeps [`turn`] from needing one.
#[derive(Clone, Debug)]
pub(crate) enum EdgeDir {
    /// A straight edge: the plane whose meet with `P` carries it, `+1` when travel runs along
    /// `n_P × n_carrier`.
    Line { carrier: usize, sense: i8 },
    /// An arc **at one of its ends**. Boxed: the exact end is a `MeetLine` and a `QuadVal`, an
    /// order of magnitude wider than a line's two words, and every direction on the hot path is a
    /// line.
    Arc(Box<ArcDir>),
    /// A ruling: travel along `±m`, the cylinder's axis — the same at both ends, like a line, but
    /// with no plane-class carrier to name a `(carrier, sense)` pair by. Boxed for the axis
    /// vector's width.
    Ruling(Box<RulingDir>),
}

/// The payload of [`EdgeDir::Ruling`] — private fields like [`ArcDir`]'s, for the same reason.
#[derive(Clone, Debug)]
pub(crate) struct RulingDir {
    cyl: usize,
    /// Which of the two parallel rulings ([`RulingCarrier::side`]) — read by [`antiparallel`] and
    /// [`continuation`], where "one carrier" means one line, not one cylinder.
    side: i8,
    /// The axis direction `m`, rational — the one geometric fact [`turn`] needs.
    axis: [nacre_exact::Rat; 3],
    /// Travel runs along `+m`.
    up: bool,
}

/// The payload of [`EdgeDir::Arc`] — private fields for the same reason the enum has them: every
/// rule that reads a direction lives in this module.
#[derive(Clone, Debug)]
pub(crate) struct ArcDir {
    cyl: usize,
    /// The end this direction is taken at, exactly.
    at: (nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    /// The circle's centre on `P` — rational, because a circle bound's plane is ⊥ to the axis.
    centre: [nacre_exact::Rat; 3],
    /// The cylinder's axis direction — what `ccw` is about, and what the travel tangent at `at`
    /// runs along as `m × (N − c)` ([`tangent_travel_agrees`]).
    axis: [nacre_exact::Rat; 3],
    /// `n_P · m > 0`: the class's **stored** normal against the cylinder's axis.
    ///
    /// ★★ **It used to say "measured unexercised, `true` on every class the corpus reaches", and
    /// that went stale the moment a new population arrived.** ☑ Re-measured with the hole trace in:
    /// [`smooth_extremum_winding`] reads it `false` on half its firings, and dropping it there moves
    /// the chained fixtures' wall. What the old note still gets right is that no *fixture* forced
    /// the factor into being — the algebra did — and [`arc_side`] records which factors its own
    /// population locks.
    axis_up: bool,
    /// Travel runs counter-clockwise about the circle's own normal (the axis).
    ///
    /// ★ This one **is** locked: dropping it makes the two arcs at a crossing share a bucket in
    /// `arrangement::angular_order`, and the walk comes back `UnorderedEdges`.
    ccw: bool,
}

/// How often the two roads ask a **cylinder** face (a ledger line). A lateral's holes are
/// among its loops, so holed laterals are not counted apart.
#[cfg(test)]
pub(crate) mod cylinder_asks {
    use std::sync::Mutex;
    pub(crate) static COUNT: Mutex<usize> = Mutex::new(0);
    pub(crate) fn asked(_f: &super::CompFace) {
        *COUNT
            .lock()
            .expect("the probe's lock is never held across a panic") += 1;
    }
}

/// How often a cut cap's candidate list holds **no** point the ring says is
/// inside — the fall-through `ring_interior_candidates`' doc calls a guard without a population.
#[cfg(test)]
pub(crate) mod witness_probe {
    use std::sync::Mutex;
    pub(crate) static NO_CANDIDATE: Mutex<Vec<String>> = Mutex::new(Vec::new());
    /// Which candidate the ring accepted: `[centre, an axis step, a chord point]` — the axis
    /// steps are the eight after the centre, the chord points follow.
    pub(crate) static ANSWERED: Mutex<[usize; 3]> = Mutex::new([0; 3]);
    pub(crate) fn answered(i: usize) {
        let k = match i {
            0 => 0,
            1..=8 => 1,
            _ => 2,
        };
        ANSWERED
            .lock()
            .expect("the probe's lock is never held across a panic")[k] += 1;
    }
    pub(crate) fn no_candidate() {
        NO_CANDIDATE
            .lock()
            .expect("the probe's lock is never held across a panic")
            .push(std::thread::current().name().unwrap_or("?").to_string());
    }
}

/// How many failed judgements the ring-vs-ring retry swallowed (predicted 0).
/// ★ That retry lives in `nesting::cell_inside`, and this counter's push with it; it feeds a
/// ledger line.
#[cfg(test)]
pub(crate) mod swallowed_probe {
    use std::sync::Mutex;
    pub(crate) static COUNT: Mutex<usize> = Mutex::new(0);
}
