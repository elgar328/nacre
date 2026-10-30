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
//! # One `usize` would mean two things — so there are two tables
//!
//! An arrangement reasons **per plane**, but a solid gives you **faces**: an earlier boolean can
//! split one geometric plane between two faces (a base's exposed top and the cantilever underside
//! above it) whose outward normals **oppose**. In one table the same `usize` would mean "face" here
//! and "plane" there, and where the two meanings meet in one comparison the answer would be wrong
//! *silently* — measured: 117748 predicate calls read "different plane" for two faces of one plane
//! and answered from rounding noise.
//!
//! That is held by the type, not by a naming convention (`fp` / `fc`): the predicates here take
//! [`crate::planes::WorkingPlane`], which has no face geometry to offer, and
//! `plane_ix` is the one place a face index becomes a plane index (in `loops::loop_triples`).
//!
//! **Exception:** the code that *defines* the classes (`planes::plane_classes` →
//! `Judge::planes_coplanar`) necessarily runs before a
//! plane table exists, so it takes face indices — hence that predicate's generic `Witness` bound.

use crate::planes::{ClassIx, WorkingPlane, edge_incidence};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Model, Solid, Vertex};
use std::collections::HashMap;

mod class_line;
mod component;
mod cyl_crossings;
mod direction;
mod in_faces;
mod incidence_order;
mod loops;
mod mixed_ring;
mod names;
mod rat;
mod rays;
mod ring_walk;
mod sort;
#[cfg(test)]
#[path = "../tests/combinatorics.rs"]
mod tests;
#[cfg(test)]
#[path = "../tests/probes/tie_probe.rs"]
pub(crate) mod tie_probe;
mod winding;
mod witnesses;

pub(crate) use class_line::*;
pub(crate) use component::*;
pub(crate) use cyl_crossings::*;
pub(crate) use direction::*;
pub(crate) use in_faces::*;
pub(crate) use incidence_order::*;
pub(crate) use loops::*;
pub(crate) use mixed_ring::*;
// A pierce node's root (`NodeKind::Pierce`) is class-space vocabulary: a result vertex names the same
// point by its rank about the cylinder (`nacre_topo::Vertex::Pierce`), and `Def::vertex` reads one as
// the other. The derived order `Lo < Hi` is read by comparison keys here — the arrangement's reject
// witness is the smallest name — so reordering the variants would silently reverse those choices.
pub(crate) use nacre_exact::quad::{QuadRoot, ruling_side_signed};
pub(crate) use nacre_judge::cylinder::PiercePoint;
pub(crate) use names::*;
pub(crate) use rat::*;
pub(crate) use rays::*;
pub(crate) use ring_walk::*;
pub(crate) use sort::*;
pub(crate) use winding::*;
pub(crate) use witnesses::*;

/// **What one witness said about a target** — a cell against a cell ([`crate::nesting`]), a
/// component's probe against another component (`assembly::grouping`): the same question one
/// dimension apart, so one vocabulary and one name for running out ([`WitnessTally`]).
pub(crate) enum Said {
    In,
    Out,
    /// This witness cannot rank itself against the target — it sits on the target's boundary, or
    /// every ray from it grazes. **The next witness is the remedy**, which is why this is not an
    /// error.
    Abstain,
    /// The value needed to ask could not be formed. The next witness may still answer, but if
    /// none does, the honest name is the arithmetic one and not "no clear ray".
    Unformed,
    /// The question reads a plane class the gate proved silent about the cylinders
    /// ([`crate::planes::WorkingCyl::is_silent`]) — a class with no world statement — and nothing
    /// here could answer it without one. Not [`Said::Unformed`]: the same arithmetic also fails on
    /// a `Wide` name or an overflow, and those are width, where this is the gate's own
    /// undecidedness carried to the place that met it.
    Silent,
    /// **The road cannot read the target** — the ray road against a ring with a pierce corner,
    /// cast from a plane with a world statement (no cylinder table there), or a two-edge lens. Not
    /// [`Said::Abstain`]: that is one witness's circumstance and the next witness of the same road
    /// is its remedy; here every witness whose name holds no class the gate passed silent would say
    /// the same, so the rest of those are skipped and the other road (coordinates) answers. A
    /// witness whose name holds one never says this: `nesting::ask` answers it by the ring's box.
    /// For the refusal's name it counts as an abstention.
    Declines,
}

/// **What a list of witnesses that all failed to decide is called** — the one rule for the name,
/// read by both roads that try witnesses in turn. The cause that a later witness could not have
/// cured wins: the gate's undecidedness ([`Said::Silent`]) over a value that could not be formed
/// ([`Said::Unformed`]) over witnesses that only grazed; with none asked at all, the cell had no
/// witness. It names the reason and does not raise it — the raise stays where the list runs out.
#[derive(Default)]
pub(crate) struct WitnessTally {
    asked: bool,
    unformed: bool,
    silent: bool,
}

impl WitnessTally {
    /// Count one witness that did not decide. A deciding answer (`In`, `Out`) is the caller's to
    /// act on and counts for nothing here.
    pub(crate) fn note(&mut self, said: &Said) {
        match said {
            Said::Abstain | Said::Declines => self.asked = true,
            Said::Unformed => self.unformed = true,
            Said::Silent => self.silent = true,
            Said::In | Said::Out => {}
        }
    }

    pub(crate) fn reason(&self) -> RejectReason {
        if self.silent {
            RejectReason::CylinderGateUndecided
        } else if self.unformed {
            RejectReason::WitnessNotRational
        } else if self.asked {
            RejectReason::NoClearRay
        } else {
            RejectReason::RingHasNoWitness
        }
    }
}

/// **A direction on a working plane**: the carrier that supplies it, and which way along it.
///
/// **The arc and ruling payloads' fields are private, and that is the point.** Every rule that
/// reads them — [`turn`], [`antiparallel`] and [`ruling_pole`] — lives in this module, and Rust's
/// module privacy is what keeps it that way: a sibling module cannot reach in even by accident. A
/// straight edge's `Line { carrier, sense }` is readable anywhere: it is the pair itself.
///
/// **A test that must read a direction is not a violation — it is an oracle.** It should read
/// **its own** inputs (the `(carrier, sense)` it built) rather than this type, both because that
/// keeps the boundary and because an oracle read back out of the value under test is an oracle
/// derived from its own answer. Do not add accessors for it.
///
/// **`carrier`, not `wall`.** A straight edge's direction is a plane carrier and a sign; an arc's
/// is a **tangent at one end**, and the two ends differ. So the type is a sum, and the arc arm
/// carries *everything the comparison needs* — the end as `(line, s)`, the circle's centre, and
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

impl EdgeDir {
    /// The same direction travelled the other way — a ring's arriving edge seen from the node it
    /// arrives at, as a departure: a line's `sense` negated, an arc's `ccw`, a ruling's `up`.
    pub(crate) fn reversed(&self) -> EdgeDir {
        match self {
            EdgeDir::Line { carrier, sense } => EdgeDir::Line {
                carrier: *carrier,
                sense: -sense,
            },
            EdgeDir::Arc(a) => EdgeDir::Arc(Box::new(ArcDir {
                ccw: !a.ccw,
                ..(**a).clone()
            })),
            EdgeDir::Ruling(r) => EdgeDir::Ruling(Box::new(RulingDir {
                up: !r.up,
                ..(**r).clone()
            })),
        }
    }
}

/// The payload of [`EdgeDir::Ruling`] — private fields like [`ArcDir`]'s, for the same reason.
#[derive(Clone, Debug)]
pub(crate) struct RulingDir {
    cyl: usize,
    /// Which of the two parallel rulings ([`RulingCarrier::side`]) — read by [`antiparallel`],
    /// where "one carrier" means one line, not one cylinder.
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
    /// Measured: `winding::smooth_extremum_winding` reads it `false` on half its firings, and
    /// dropping it there moves the chained fixtures' wall. No *fixture* forced the factor into
    /// being — the algebra did — and `direction::arc_side` records which factors its own population
    /// locks.
    axis_up: bool,
    /// Travel runs counter-clockwise about the circle's own normal (the axis).
    ///
    /// This one **is** locked: dropping it makes the two arcs at a crossing share a bucket in
    /// [`angular_order`], and the walk comes back `UnorderedEdges`.
    ccw: bool,
}
