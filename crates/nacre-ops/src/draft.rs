//! **What the engine says a result face is, before anything is built.**
//!
//! The arrangement decides a face's boundary long before a `Handle` exists for any of it, and the
//! assembly turns that description into b-rep. Both halves speak this vocabulary, so it belongs to
//! neither: [`LocalFace`] with its [`Bound`]s, a [`Ring`] of [`NodeId`]s and the [`Wall`] each edge
//! rides, a [`Rim`], a [`SeamVertex`], and the record a cut circle carries out of the split
//! ([`CutRim`]). Alongside them the two terms that decide whether a draft face survives at all --
//! [`BoolKind`] and [`keep`], one predicate on one chamber's `(inA, inB)`.
//!
//! ★ **This module is under the engine, not beside it.** It names `combinatorics` (a ring is a ring
//! of names) and `planes` (a face lives in a class table) and nothing else in this crate -- no
//! `arrangement`, no `assembly`, no `ops`. That is what a shared vocabulary means here, and the
//! module-graph gate asserts exactly it.

use crate::combinatorics::{NodeId, Wall};
use crate::planes::{ClassIx, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, combinatorics, reject};
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Surface, Vertex};
use std::collections::HashMap;
/// A seam vertex — a three-plane point on both `∂A` and `∂B` (2 A-planes + 1
/// B-plane, or 1 A + 2 B). Shared (one `Handle`) by every incident result piece.
pub(crate) struct SeamVertex {
    pub(crate) point: Point3,
    pub(crate) triple: NodeId,
    pub(crate) tol: f64,
}

/// One ring of a result face: its nodes, and **the plane each edge rides**.
///
/// ★ **The walls are carried, not derived.** Reading an edge's supporting plane back out of its two
/// endpoint names is sound only while every vertex lies on exactly three planes — see
/// [`Ring::edges`]. Every producer here knows the wall (the arrangement's
/// half-edge was told it; a pass-through face reads it off the edge's other face), so it hands it
/// over instead of leaving it to be guessed.
///
/// Derefs to its nodes, so the many places that only walk the ring read unchanged.
#[derive(Clone, Debug, Default)]
pub(crate) struct Ring {
    pub(crate) nodes: Vec<NodeId>,
    /// `walls[i]` is the carrier of the edge `nodes[i] -> nodes[i+1]`.
    pub(crate) walls: Vec<Wall>,
}

/// **The edge-welding key** — the key half of "a line is unordered, a circle is ordered"
/// (`Wall` is the type half).
///
/// Keys live in *handle* space, the space `edge_of` always keyed. A line is its unordered
/// endpoint pair, as before. Between one pair of pierce vertices a circle offers **two**
/// complementary pieces, so the endpoints alone cannot name an arc — the key carries them **in
/// CCW order about the axis** (`from → to`), and the two complementary arcs get the two orders.
/// The minted edge stores its vertices in that same order, which is what the `[A, B]`-CCW
/// convention means downstream (`Model::derive_edge_curve`'s circle arm).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum EdgeKey {
    Line((usize, usize)),
    Arc {
        cyl: usize,
        from: usize,
        to: usize,
    },
    /// A ruling piece: straight, so the pair is unordered like a line's — keyed apart from
    /// plane edges by `(cyl, side)`, mirroring [`Wall::Ruling`]'s identity (a plane edge
    /// collinear with a ruling is a refused degeneracy, not a legal share).
    Ruling {
        cyl: usize,
        side: i8,
        pair: (usize, usize),
    },
}

impl Ring {
    /// A ring whose walls are **derived from the node names** — for hand-built fixtures, where
    /// every vertex is a clean three-plane point so "the class the two names share besides `p`" is
    /// well defined. Production never derives; see the type's note.
    #[cfg(test)]
    pub(crate) fn from_clean_names(p: usize, nodes: Vec<NodeId>) -> Ring {
        let k = nodes.len();
        let walls = (0..k)
            .map(|i| {
                // Its own contract is "clean fixture rings only", so a pierce node here is a
                // fixture bug, not an input the kernel must survive.
                let name = |n| {
                    combinatorics::three_plane_name(n).expect("a clean fixture ring names triples")
                };
                let (a, b) = (name(nodes[i]), name(nodes[(i + 1) % k]));
                Wall::Plane(
                    a.iter()
                        .copied()
                        .find(|&c| c != p && b.contains(&c))
                        .expect("a clean fixture ring edge rides one wall"),
                )
            })
            .collect();
        Ring { nodes, walls }
    }
}

impl std::ops::Deref for Ring {
    type Target = [NodeId];
    fn deref(&self) -> &[NodeId] {
        &self.nodes
    }
}

impl Ring {
    pub(crate) fn new(nodes: Vec<NodeId>, walls: Vec<Wall>) -> Ring {
        debug_assert_eq!(nodes.len(), walls.len(), "one wall per edge");
        Ring { nodes, walls }
    }

    /// The carrier-level reading of [`combinatorics::ring_is_mixed`], asked of the
    /// **walls and nodes** rather than of derived edges - the legacy `edges()` below dies
    /// on a mixed ring before any `RingEdge` exists, so a consumer that must pass such a
    /// ring over has to ask the `Ring` itself.
    pub(crate) fn is_mixed(&self) -> bool {
        self.walls.iter().any(|w| !matches!(w, Wall::Plane(_)))
            || self
                .nodes
                .iter()
                .any(|&n| combinatorics::pierce_name(n).is_some())
    }

    /// This ring's edges, ready for the exact predicates - built from the walls the ring
    /// carries, in the engine's own vocabulary.
    ///
    /// A wall becomes its carrier (`Wall::Arc` an [`combinatorics::ArcCarrier`] with the class
    /// table's def - the same clone convention every producer follows), a pierce end is pinned by
    /// its cylinder, and a three-plane end by `pin_on_line` - so a mixed ring yields edges its
    /// consumers fork on (`ring_is_mixed`). Flattening walls to plane indices for a plane-only
    /// derivation would refuse every mixed ring here, before any consumer could abstain.
    ///
    /// `p` is the plane class the ring lies in, which pins a **three-plane** end on its edge's
    /// line; a lateral ring (on a cylinder) has none — its corners are pierce names — and passes
    /// `None`, so a three-plane end there is refused by name rather than pinned by a guess.
    pub(crate) fn edges(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        p: Option<usize>,
    ) -> Result<Vec<combinatorics::RingEdge>, BoolError> {
        if self.nodes.len() != self.walls.len() {
            return Err(reject(RejectReason::RingNaming));
        }
        let pin = |n: NodeId, wall: &Wall| -> Result<combinatorics::EndPin, BoolError> {
            if combinatorics::pierce_name(n).is_some() {
                return Ok(combinatorics::EndPin::Cylinder);
            }
            let Some(t) = combinatorics::three_plane_name(n) else {
                return Err(reject(RejectReason::RingNaming));
            };
            let (Wall::Plane(c), Some(p)) = (wall, p) else {
                // A three-plane point on a curved carrier, or on a lateral ring, has no third
                // plane to pin it with - no producer builds the shape today, and it is refused
                // by name rather than guessed at.
                return Err(reject(RejectReason::RingNaming));
            };
            combinatorics::pin_on_line(jd, p, *c, t)
                .map(combinatorics::EndPin::Class)
                .ok_or_else(|| reject(RejectReason::RingNaming))
        };
        let def_of = |cyl: usize| -> Result<nacre_topo::CylinderDef, BoolError> {
            Ok(cyls
                .get(cyl)
                .ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?
                .def
                .clone())
        };
        (0..self.nodes.len())
            .map(|i| {
                let jn = (i + 1) % self.nodes.len();
                let wall = &self.walls[i];
                let carrier = match wall {
                    Wall::Plane(c) => combinatorics::Carrier::plane(*c),
                    Wall::Arc { cyl, ccw } => {
                        combinatorics::Carrier::Arc(Box::new(combinatorics::ArcCarrier {
                            cyl: *cyl,
                            def: def_of(*cyl)?,
                            ccw: *ccw,
                        }))
                    }
                    Wall::Ruling { cyl, side, up } => {
                        combinatorics::Carrier::Ruling(Box::new(combinatorics::RulingCarrier {
                            cyl: *cyl,
                            def: def_of(*cyl)?,
                            side: *side,
                            up: *up,
                        }))
                    }
                };
                Ok(combinatorics::RingEdge {
                    node: self.nodes[i],
                    to: self.nodes[jn],
                    carrier,
                    from_h: pin(self.nodes[i], wall)?,
                    to_h: pin(self.nodes[jn], wall)?,
                })
            })
            .collect()
    }
}

/// One boundary of a result face: a polygon of seam nodes, a **full circle** of a
/// cylinder class, or a lateral **band** between two rims. A circle has no nodes or walls and is
/// assembled through the rim machinery (`push_edge([lateral, plane], [v, v])` + an `OnSeam`
/// vertex), never through the seam-vertex table; so is a band's rim while it is a whole circle —
/// a band's **chain** rim is a node ring like a polygon, and goes through both.
#[derive(Clone, Debug)]
pub(crate) enum Bound {
    Ring(Ring),
    Circle {
        cyl: usize,
    },
    /// A lateral band's whole boundary: its two rims, `lo` the one walked forward (winding `+1`
    /// about the axis) and `hi` the one walked backward (`−1`) — for two whole circles, the one
    /// with the smaller axis parameter and the larger. The face it bounds is the cylinder
    /// itself, so the class is on [`LocalFace::surf`] rather than repeated here.
    Band {
        lo: Rim,
        hi: Rim,
    },
}

/// One rim of a lateral band: the **whole circle** of a plane class, or a **wrapping chain** —
/// a closed ring of arcs and rulings going once around the cylinder, the shape the cleaning
/// pass leaves when a band and a panel merge (a boss whose cap sits inside the other body). A
/// chain is stored in the direction it is walked (a `lo` chain forward, a `hi` chain backward)
/// and unflipped, like every bound; its nodes are polygon nodes ([`Bound::rings`]), so the seam
/// table, the vertex minting and the node join all read it as they read a panel ring.
#[derive(Clone, Debug)]
pub(crate) enum Rim {
    Circle(usize),
    Chain(Ring),
}

impl Rim {
    /// The plane class, `None` for a chain.
    pub(crate) fn circle(&self) -> Option<usize> {
        match self {
            Rim::Circle(c) => Some(*c),
            Rim::Chain(_) => None,
        }
    }

    /// The chain's ring, `None` for a whole circle.
    pub(crate) fn ring(&self) -> Option<&Ring> {
        match self {
            Rim::Circle(_) => None,
            Rim::Chain(r) => Some(r),
        }
    }

    fn ring_mut(&mut self) -> Option<&mut Ring> {
        match self {
            Rim::Circle(_) => None,
            Rim::Chain(r) => Some(r),
        }
    }
}

impl Bound {
    /// The polygon ring, `None` for a curved bound — the consumers that want a *polygon
    /// boundary* (not merely nodes) filter on this; a band's chain rims are not one.
    pub(crate) fn ring(&self) -> Option<&Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } | Bound::Band { .. } => None,
        }
    }

    /// Every **node ring** this bound carries: the polygon itself, or a band's chain rims (a
    /// whole circle carries none). The node-level consumers — the seam table, vertex minting,
    /// the node join, the naming tables — read this, so a chain's nodes exist wherever a panel
    /// ring's do.
    pub(crate) fn rings(&self) -> impl Iterator<Item = &Ring> {
        let (a, b) = match self {
            Bound::Ring(r) => (Some(r), None),
            Bound::Band { lo, hi } => (lo.ring(), hi.ring()),
            Bound::Circle { .. } => (None, None),
        };
        a.into_iter().chain(b)
    }

    pub(crate) fn rings_mut(&mut self) -> impl Iterator<Item = &mut Ring> {
        let (a, b) = match self {
            Bound::Ring(r) => (Some(r), None),
            Bound::Band { lo, hi } => (lo.ring_mut(), hi.ring_mut()),
            Bound::Circle { .. } => (None, None),
        };
        a.into_iter().chain(b)
    }

    /// The polygon ring, asserted — for consumers whose population cannot carry curved bounds
    /// (and tests). Panics with the caller's location.
    #[cfg_attr(not(test), allow(dead_code))]
    #[track_caller]
    pub(crate) fn expect_ring(&self) -> &Ring {
        match self {
            Bound::Ring(r) => r,
            Bound::Circle { cyl } => panic!("a polygon-only path got a circle bound (cyl {cyl})"),
            Bound::Band { lo, hi } => {
                panic!("a polygon-only path got a band bound ({lo:?}..{hi:?})")
            }
        }
    }
}

/// The rim a curved bound is assembled on: `(group, cylinder class, plane class) → (seam vertex,
/// rim edge)`. The group is in the key for the reason it is in the seam vertices' — two result
/// solids never share a handle — and the rest is what makes a cap face and the band that meets it
/// pick up the *same* edge. A **cut** circle's entry carries its seam vertex and `None`: the
/// closed `[v, v]` edge spelling is false for it, and its boundary is assembled from arc pieces
/// instead.
pub(super) type RimTable = HashMap<(usize, usize, usize), (Handle<Vertex>, Option<Handle<Edge>>)>;

/// A reconstructed result face: which combined plane it is on, its boundaries, and whether to
/// flip it (cut's inside-A B-pieces).
#[derive(Clone, Debug)]
pub(crate) struct LocalFace {
    /// Which class table this face's surface lives in — a plane class, or a cylinder
    /// class whose band this face is. Plane-only consumers project with [`ClassIx::plane`], whose
    /// panic is the upstream-filter-bug detector the type was introduced with.
    pub(crate) surf: crate::planes::ClassIx,
    pub(crate) outer: Bound,
    /// Hole boundaries, each already wound so the kept material stays on its left
    /// about the face's outward normal.
    pub(crate) inner: Vec<Bound>,
    pub(crate) flip: bool,
}

impl LocalFace {
    /// Every **polygon** ring of the face, outer first. Circle bounds carry no nodes and are
    /// deliberately absent — node-level consumers (the seam table, edge scans, self-touch)
    /// have nothing to read off them.
    pub(crate) fn poly_rings(&self) -> impl Iterator<Item = &Ring> {
        std::iter::once(&self.outer)
            .chain(self.inner.iter())
            .flat_map(Bound::rings)
    }

    pub(crate) fn poly_rings_mut(&mut self) -> impl Iterator<Item = &mut Ring> {
        std::iter::once(&mut self.outer)
            .chain(self.inner.iter_mut())
            .flat_map(Bound::rings_mut)
    }
}

/// **A ring step's passage of θ = 0, signed** — the one seam predicate, from which both the
/// loop builder's split rule ([`wrapping_rim`]) and the cleaning pass's winding count derive.
///
/// `Some((cyl, plane, sign))` when the step is an arc on a **cut** rim whose CCW parametrization
/// passes θ = 0 — `sign` is `+1` walked CCW, `−1` walked CW — and `None` for a ruling, a plane
/// wall, or an arc that stays clear of the seam. Two facts, neither re-derived by a caller:
///
/// - **which circle's plane the arc rides** — a *plane* face's own class (a cap), or, for a face
///   on the cylinder class (a panel), the plane its two end names share. A panel arc's ends may
///   share **both** planes (one wall cuts a circle twice, so each end is named
///   {wall, circle-plane}), so the unique-share rule a ruling edge uses cannot apply: the circle's
///   plane is the shared candidate that carries a cut-rim record. Two such candidates is a naming
///   this ladder does not arrange (`RulingBoundNotYet`).
/// - **the directed passage test** — with two pierce nodes the two complementary arcs share one
///   unordered endpoint pair, so the step is oriented by the `ccw` bit the wall carries and
///   compared against the split's own θ order: the CCW arc `nodes.last() → nodes[0]` is the one
///   holding θ = 0. ★ **Half-open.** When the seam *is* a node (`CutRim::seam_is_node`), that arc
///   is the one **ending** at the seam node, and it passes θ = 0; the arc leaving the seam node
///   does not. So a loop that runs *along* the seam ruling counts its arrival and its departure
///   once between them, never twice — which is what makes Σ sign over a cycle its winding number
///   about the axis (`+1` a lower boundary walked forward, `−1` an upper one walked backward, `0`
///   a hole), read off the split's order table and no coordinate.
pub(super) fn seam_step(
    surf: ClassIx,
    a: NodeId,
    b: NodeId,
    wall: Wall,
    cut_rims: &crate::draft::CutRims,
) -> Result<Option<(usize, usize, i8)>, BoolError> {
    let Wall::Arc { cyl, ccw } = wall else {
        return Ok(None);
    };
    let c = match surf {
        ClassIx::Plane(c) => c,
        ClassIx::Cyl(_) => {
            let shared = |x| {
                let (pa, _, _) = combinatorics::pierce_name(a)?;
                let (pb, _, _) = combinatorics::pierce_name(b)?;
                (pa.contains(&x) && pb.contains(&x)).then_some(x)
            };
            let mut hits = combinatorics::pierce_name(a)
                .map(|(pa, _, _)| pa)
                .into_iter()
                .flatten()
                .filter_map(shared)
                .filter(|&c| cut_rims.contains_key(&(cyl, c)));
            let c = hits
                .next()
                .ok_or_else(|| reject(RejectReason::RulingBoundNotYet))?;
            if hits.next().is_some() {
                return Err(reject(RejectReason::RulingBoundNotYet));
            }
            c
        }
    };
    let Some(cr) = cut_rims.get(&(cyl, c)) else {
        return Ok(None);
    };
    let ccw_pair = if ccw { (a, b) } else { (b, a) };
    let last = *cr.nodes.last().expect("a cut circle has pierce nodes");
    Ok((ccw_pair == (last, cr.nodes[0])).then_some((cyl, c, if ccw { 1 } else { -1 })))
}

/// **The cut rim a ring edge wraps past θ = 0 on**, `None` if it does not — [`seam_step`] for
/// the loop builder, which splits such an edge at the rim's seam vertex (and the rim table has to
/// have minted that vertex first). A rim whose seam **is** a node splits nothing, so it answers
/// `None` here even though the step passes the seam.
pub(super) fn wrapping_rim(
    surf: ClassIx,
    a: NodeId,
    b: NodeId,
    wall: Wall,
    cut_rims: &crate::draft::CutRims,
) -> Result<Option<(usize, usize)>, BoolError> {
    Ok(seam_step(surf, a, b, wall, cut_rims)?
        .filter(|&(cyl, c, _)| !cut_rims[&(cyl, c)].seam_is_node)
        .map(|(cyl, c, _)| (cyl, c)))
}

/// Which boolean to compute.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolKind {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

/// The boolean keep predicate on one chamber's `(inA, inB)`.
pub(crate) fn keep(kind: BoolKind, in_a: bool, in_b: bool) -> bool {
    match kind {
        BoolKind::Fuse => in_a || in_b,
        BoolKind::Cut => in_a && !in_b,
        BoolKind::Common => in_a && in_b,
    }
}

/// **A cut circle's seam datum — carried from the split, never re-derived.** `split_circles`
/// already orders a cut circle's pierce nodes by θ about the seam and classifies a seam-incident
/// node by name (`circular_order_about_seam`), so the one fact the assembly cannot re-derive
/// cheaply — *where θ = 0 sits among the arcs* — travels from the place that computed it.
#[derive(Clone, Debug)]
pub(crate) struct CutRim {
    /// The circle's pierce nodes in θ order (CCW about the axis). When `seam_is_node`, the
    /// seam-incident node is first; otherwise θ = 0 lies inside the wrap arc
    /// `nodes.last() → nodes[0]`.
    pub(crate) nodes: Vec<combinatorics::NodeId>,
    pub(crate) seam_is_node: bool,
}

/// Per `(cylinder class, plane class)`, the cut circles — presence in this map **is** the one
/// source of "this rim is cut" (the assembly's rim skip and curved arms all read it).
pub(crate) type CutRims = HashMap<(usize, usize), CutRim>;

/// What each input face's plane became: its **plane class's representative surface**, which is
/// what every result face on that plane carries — and the `Orientation` the face would carry stated
/// on that surface (`None` where the judge could not say which way the two run, read off the
/// truth by `planes::face_facing`).
pub(crate) type ClassOf =
    std::collections::HashMap<Handle<Face>, (Handle<Surface>, Option<nacre_topo::Orientation>)>;
