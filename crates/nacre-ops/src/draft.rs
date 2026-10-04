//! **What the engine says a result face is, before anything is built.**
//!
//! The arrangement decides a face's boundary long before a `Handle` exists for any of it, and the
//! assembly turns that description into b-rep. Both halves speak this vocabulary, so it belongs to
//! neither: [`LocalFace`] with its [`Bound`]s, a [`Ring`] of [`NodeId`]s and the [`Wall`] each edge
//! rides, a [`Rim`], a [`SeamVertex`] and the [`Def`] a node's vertex is minted from, and the
//! record a cut circle carries out of the split ([`CutRim`]) with the part of it the result's
//! faces keep ([`HeldRims`]). Alongside them the two terms that decide whether a draft face survives at all --
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
use nacre_store::Handle;
use nacre_topo::{Edge, Vertex};
use std::collections::HashMap;
/// **A result vertex's definition, in class space** — what becomes a `Vertex` once the classes
/// are read as surface handles ([`Def::vertex`]).
///
/// ★ Two readers ask for it, and they must ask the same thing: the seam table realizes a node's
/// coordinate from the definition the minting will push, so the two cannot describe one vertex
/// two ways. `Three` is a triple the naming derives from the faces through the node (or, in the
/// seam table, the node's own triple — the same point); `Pierce` is a **declaration, not a
/// derivation**: `NodeId::Pierce` already names two result plane classes and the cylinder, so its
/// def is the name's own payload ([`Def::of_pierce_name`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Def {
    Three(combinatorics::Canon3),
    Pierce {
        planes: [usize; 2],
        cyl: usize,
        root: nacre_topo::QuadRoot,
    },
}

impl Def {
    /// A pierce node's definition — its name's payload, or `None` for a three-plane node.
    pub(crate) fn of_pierce_name(node: NodeId) -> Option<Def> {
        let (planes, cyl, root) = combinatorics::pierce_name(node)?;
        Some(Def::Pierce { planes, cyl, root })
    }

    /// The definition a node's own name states — its triple, or its pierce payload. The seam table
    /// realizes this one; the naming may define the vertex by another triple through the same
    /// point where four planes concur.
    pub(crate) fn of_name(node: NodeId) -> Option<Def> {
        Def::of_pierce_name(node).or_else(|| {
            combinatorics::three_plane_name(node)
                .map(|t| Def::Three(combinatorics::Canon3::three(t)))
        })
    }

    /// The model's vertex for this definition — the class → handle mapping and nothing else.
    ///
    /// ★★★ **That mapping is where `QuadRoot::canonical` answers a second time**: `NodeId::Pierce`
    /// is canonical in *class* order, `Vertex::Pierce` in *handle* order, and the class → handle
    /// map is not monotone in general — a re-sort must carry the root through (`transform`'s
    /// remap already locks the same rule on the way back out).
    /// ★ The flip is unexercised **at the minting**: dropping it leaves every fence green, because
    /// both fixtures' class order happens to match their handle order. ★★★★★ **The rule is not
    /// unexercised, though — its inverse is red.** The population that exercises it is a
    /// boolean's *result used as the next operand*, where a second boolean builds its classes
    /// afresh and their order is not the handles': dropping the same restatement in
    /// `combinatorics::pierce_name_from_def` turns
    /// `an_operand_bounded_by_a_cylinder_is_named_in_class_space` red. What is still owed is only a
    /// fixture that reaches the minting, not the rule.
    pub(crate) fn vertex(
        &self,
        planes: &[WorkingPlane],
        cyls: &[crate::planes::WorkingCyl],
    ) -> Vertex {
        match *self {
            Def::Three(t) => {
                let t = t.planes();
                Vertex::ThreePlane([planes[t[0]].surf, planes[t[1]].surf, planes[t[2]].surf])
            }
            Def::Pierce {
                planes: p2,
                cyl,
                root,
            } => {
                let (pair, root) =
                    nacre_topo::QuadRoot::canonical([planes[p2[0]].surf, planes[p2[1]].surf], root);
                Vertex::Pierce {
                    planes: pair,
                    cylinder: cyls[cyl].surf,
                    root,
                }
            }
        }
    }
}

/// A seam vertex — a point on both `∂A` and `∂B` (a three-plane corner or a pierce point), shared
/// (one `Handle`) by every incident result piece.
///
/// ★ **It carries the cache its vertex receives, and nothing measured beside it.** `Bounded`
/// says the truth lies within `coord ± bound` — the one thing the self-touch sieve needs to keep
/// its boxes conservative — and any other variant says no bound is proven, which the sieve reads
/// as "cannot rule this out". A residual (the point's distance to its carriers' `f64` planes) is
/// not such a bound: measured against the nearest `f64` of the truth it fell short on 3,415 of the
/// suite's seam points, 155 of them by more than a thousandfold.
pub(crate) struct SeamVertex {
    /// The cache the vertex of this node's own definition receives ([`Def::of_name`] through
    /// `realize::point_cache`).
    pub(crate) cache: nacre_topo::PointCache,
    pub(crate) triple: NodeId,
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
/// Keys live in *handle* space, the space `edge_of` always keyed. A **straight** piece is its
/// unordered endpoint pair — a plane-pair line and a ruling alike: two points bound one straight
/// segment, so a line one face states in the plane vocabulary and the face across states as a
/// ruling (a tangent wall ending on a cylinder's ruling) is one edge, and its carriers are read
/// off the faces that use it, not off the key. Between one pair of pierce vertices a circle
/// offers **two** complementary pieces, so the endpoints alone cannot name an arc — the key
/// carries them **in CCW order about the axis** (`from → to`), and the two complementary arcs
/// get the two orders. The minted edge stores its vertices in that same order, which is what the
/// `[A, B]`-CCW convention means downstream (`Model::derive_edge_curve`'s circle arm).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum EdgeKey {
    Line((usize, usize)),
    Arc { cyl: usize, from: usize, to: usize },
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
    /// table's def - the same clone convention every producer follows), a pierce end on a plane
    /// edge by `pin_for` (its cylinder when its pair is the edge's line) and elsewhere by its
    /// cylinder, and a three-plane end by `pin_on_line` - so a mixed ring yields edges its
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
                // ★ A pierce end is pinned by its cylinder only on the line its pair spells; a
                // point on a line two classes share on the cylinder can arrive under the other
                // pair's name, and on a plane edge its pin is then derived — `pin_for`, the one
                // rule («a pin is a function of the name and the line»).
                return match (wall, p) {
                    (Wall::Plane(c), Some(p)) => combinatorics::pin_for(jd, p, *c, n)
                        .ok_or_else(|| reject(RejectReason::RingNaming)),
                    _ => Ok(combinatorics::EndPin::Cylinder),
                };
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

/// One boundary of a result face: a polygon of seam nodes, or a **whole circle** — on a plane
/// face, the circle a cylinder class cuts in it (`Circle`); on a lateral, the rim a plane class
/// cuts on it (`Rim`). A whole circle has no nodes or walls and is assembled through the rim
/// machinery (`push_edge([lateral, plane], [v, v])` + an `OnSeam` vertex), never through the
/// seam-vertex table. A lateral that wraps the axis is bounded by two of its own rims — each a
/// `Rim`, or a ring of arcs and rulings where the rim is cut (a chain) — and no seam edge.
#[derive(Clone, Debug)]
pub(crate) enum Bound {
    Ring(Ring),
    Circle {
        cyl: usize,
    },
    /// A lateral's whole rim on plane class `plane`; `ccw` when the face lies above it on the
    /// chart — walked counter-clockwise about the axis (the lower rim), against it for the upper.
    /// The face it bounds is the cylinder itself, so the class is on [`LocalFace::surf`].
    Rim {
        plane: usize,
        ccw: bool,
    },
}

impl Bound {
    /// The polygon ring, `None` for a whole circle — the consumers that want a *polygon
    /// boundary* filter on this, and the node-level ones (the seam table, vertex minting, the
    /// node join, the naming tables) read every bound's through it.
    pub(crate) fn ring(&self) -> Option<&Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } | Bound::Rim { .. } => None,
        }
    }

    pub(crate) fn ring_mut(&mut self) -> Option<&mut Ring> {
        match self {
            Bound::Ring(r) => Some(r),
            Bound::Circle { .. } | Bound::Rim { .. } => None,
        }
    }

    /// The polygon ring, asserted — for consumers whose population cannot carry curved bounds
    /// (and tests). Panics with the caller's location.
    #[cfg_attr(not(test), allow(dead_code))]
    #[track_caller]
    pub(crate) fn expect_ring(&self) -> &Ring {
        match self {
            Bound::Ring(r) => r,
            Bound::Circle { cyl } => panic!("a polygon-only path got a circle bound (cyl {cyl})"),
            Bound::Rim { plane, .. } => {
                panic!("a polygon-only path got a lateral's rim (plane {plane})")
            }
        }
    }
}

/// The whole circle a curved bound is assembled on: `(group, cylinder class, plane class) →
/// (seam vertex, closed rim edge)`. The group is in the key for the reason it is in the seam
/// vertices' — two result solids never share a handle — and the rest is what makes a cap face and
/// the lateral that meets it pick up the *same* edge. Only an uncut circle has one: a cut circle's
/// boundary is its arcs between nodes, minted by the rings that walk them.
pub(super) type RimTable = HashMap<(usize, usize, usize), (Handle<Vertex>, Handle<Edge>)>;

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
            .filter_map(Bound::ring)
    }

    pub(crate) fn poly_rings_mut(&mut self) -> impl Iterator<Item = &mut Ring> {
        std::iter::once(&mut self.outer)
            .chain(self.inner.iter_mut())
            .filter_map(Bound::ring_mut)
    }
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

/// **A cut circle's nodes in their order round it — carried from the split, never re-derived.**
/// `split_circles` orders a cut circle's pierce nodes by θ (`circular_order_about_seam`), so the
/// arcs between consecutive nodes — the circle's pieces — travel from the place that computed
/// them.
#[derive(Clone, Debug)]
pub(crate) struct CutRim {
    /// The circle's pierce nodes in θ order (CCW about the axis).
    pub(crate) nodes: Vec<combinatorics::NodeId>,
}

/// Per `(cylinder class, plane class)`, how the arrangement **split** each circle it cut — the
/// table the lateral chart reads its arc labels against (`Curved::split_rims`). It is not the
/// result's rim: the cleaning pass merges cap cells back and dissolves the nodes between them, so
/// what the result's faces and edges meet on is [`HeldRims`].
pub(crate) type CutRims = HashMap<(usize, usize), CutRim>;

/// **The result's cut rims — the split's nodes that the cleaned faces still hold.** Presence here
/// is the one source of "this rim is cut" for everything that builds the result (the lateral's
/// rim pieces, the assembly's rim skip, arc join and seam predicate). A type of its own because
/// the defect it closes was one of two same-typed tables read in the wrong place: a Cut whose
/// tool rests on a cap across its rim merged the cap back into one arc while the lateral, reading
/// the split, still cut its rim at the tool's planes — one arc on the cap, three on the lateral,
/// an edge used once.
#[derive(Clone, Debug, Default)]
pub(crate) struct HeldRims(CutRims);

impl HeldRims {
    pub(crate) fn get(&self, key: &(usize, usize)) -> Option<&CutRim> {
        self.0.get(key)
    }

    pub(crate) fn contains_key(&self, key: &(usize, usize)) -> bool {
        self.0.contains_key(key)
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (&(usize, usize), &CutRim)> {
        self.0.iter()
    }

    /// This table re-read against a later face list — for the readers after the per-solid
    /// straight-angle pass, which can drop nodes of its own. Only removes: a node the later list
    /// adds is not a rim node this table knew. Blind to groups (the key is `(cyl, plane)`), so a
    /// node one body dropped and another kept stays for both — no per-solid pass drops a rim node
    /// in the suite or the census (the two pierce nodes it drops there are ruling points).
    pub(crate) fn prune(&self, faces: &[LocalFace]) -> HeldRims {
        held_rims(faces, &self.0)
    }
}

/// [`HeldRims`] from the faces as they stand: for each split circle `(cyl, c)`, the nodes some
/// face on plane class `c` still has in a ring, in the split's θ order. A circle none of whose
/// nodes survive is **whole** again and
/// leaves the table — the lateral emits its rim as the closed edge, the cap holds it as a circle.
///
/// Any node of a class-`c` ring, not only an arc's ends: the inward wedge with its apex on the
/// seam has the cap turn at the apex between two lines, and the arc-ends rule dropped that node
/// — the lateral's corner there — and moved that `A − B`'s refusal to `ArcBoundNotYet` (measured).
///
/// ★ A circle left with **one** node keeps the split's nodes — the record as it was before this
/// derivation, so an input it cannot improve answers as it did rather than by a new refusal (one
/// node cannot state a rim as arcs, and the rim may be one no lateral reaches). No suite or census
/// boolean reaches it: the rims this derivation shortens go from three or four nodes to two.
pub(crate) fn held_rims(faces: &[LocalFace], split: &CutRims) -> HeldRims {
    if split.is_empty() {
        return HeldRims::default();
    }
    let mut held: std::collections::HashSet<(usize, combinatorics::NodeId)> =
        std::collections::HashSet::new();
    for lf in faces {
        let ClassIx::Plane(c) = lf.surf else {
            continue;
        };
        for ring in lf.poly_rings() {
            held.extend(ring.nodes.iter().map(|&n| (c, n)));
        }
    }
    let mut out = CutRims::new();
    for (&(cyl, c), cr) in split {
        let nodes: Vec<_> = cr
            .nodes
            .iter()
            .copied()
            .filter(|&n| held.contains(&(c, n)))
            .collect();
        let rim = match nodes.len() {
            0 => continue,
            1 => cr.clone(),
            _ => CutRim { nodes },
        };
        out.insert((cyl, c), rim);
    }
    HeldRims(out)
}
