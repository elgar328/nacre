//! **What the engine says a result face is, before anything is built.**
//!
//! The arrangement decides a face's boundary long before a `Handle` exists for any of it, and the
//! assembly turns that description into b-rep. Both halves speak this vocabulary, so it belongs to
//! neither: [`LocalFace`] with its [`Bound`]s, a [`Ring`] of [`NodeId`]s and the [`Wall`] each edge
//! rides, a [`RimTable`], a [`SeamVertex`] and the [`Def`] a node's vertex is minted from, and the
//! record a cut circle carries out of the split ([`CutRim`]) with the part of it the result's
//! faces keep ([`HeldRims`]). Alongside them the two terms that decide whether a draft face survives at all --
//! [`BoolKind`] and [`keep`], one predicate on one chamber's `(inA, inB)`.
//!
//! **This module is under the engine, not beside it.** It names `combinatorics` (a ring is a ring
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
/// Two readers ask for it, and they must ask the same thing: the seam table realizes a node's
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
        root: crate::combinatorics::QuadRoot,
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

    /// The model's vertex for this definition — the class → handle mapping, and for a pierce
    /// point the class name read as the result's.
    ///
    /// **A pierce name says which crossing along its classes' meet line; the vertex says which
    /// place about the cylinder** (`Vertex::Pierce`'s `rank`). So the meet is solved on the class
    /// statements, its crossings ranked about the cylinder ([`nacre_exact::quad::seam_ranked`]),
    /// and the name's crossing found among them. The planes' handles are then stored ascending,
    /// which the rank does not read.
    ///
    /// Refused by name where the class statements give no exact meet: a class with no world
    /// statement (`CylinderGateUndecided` — one the gate let through silent, whose pierce points
    /// are interval boxes) or one `Rat` cannot hold (`WitnessNotRational`, as for an overflow); a
    /// name whose crossing the meet does not hold is an earlier stage's record broken
    /// (`CylinderStagesDisagree`).
    pub(crate) fn vertex(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        planes: &[WorkingPlane],
        cyls: &[crate::planes::WorkingCyl],
    ) -> Result<Vertex, BoolError> {
        Ok(match *self {
            Def::Three(t) => {
                let t = t.planes();
                Vertex::ThreePlane([planes[t[0]].surf, planes[t[1]].surf, planes[t[2]].surf])
            }
            Def::Pierce {
                planes: p2,
                cyl,
                root,
            } => {
                let width = || reject(RejectReason::WitnessNotRational);
                let coeffs = |k: usize| {
                    jd.world_rat(k).ok_or_else(|| {
                        reject(if planes[k].world.is_some() {
                            RejectReason::WitnessNotRational
                        } else {
                            RejectReason::CylinderGateUndecided
                        })
                    })
                };
                let def = &cyls[cyl].def;
                let (o, m) = (def.origin(), def.dir());
                let meet = nacre_exact::quad::plane_plane_cylinder(
                    &coeffs(p2[0])?,
                    &coeffs(p2[1])?,
                    &o,
                    &m,
                    def.r2(),
                )
                .ok_or_else(width)?;
                // The name was minted from this same meet on these same statements, so a meet
                // without its crossing is an earlier stage's record broken.
                let (_, named) = root
                    .crossing(&meet)
                    .ok_or_else(|| reject(RejectReason::CylinderStagesDisagree))?;
                let ranked = nacre_exact::quad::seam_ranked(meet, &o, &m, &def.ref_dir())
                    .ok_or_else(width)?;
                let mut rank = None;
                for (k, (_, s)) in ranked.iter().enumerate() {
                    if s.checked_sub(&named).ok_or_else(width)?.sign() == nacre_exact::Orient::Zero
                    {
                        rank = Some(k);
                    }
                }
                let rank = rank.ok_or_else(|| reject(RejectReason::CylinderStagesDisagree))?;
                let [a, b] = [planes[p2[0]].surf, planes[p2[1]].surf];
                Vertex::Pierce {
                    planes: if b < a { [b, a] } else { [a, b] },
                    cylinder: cyls[cyl].surf,
                    rank: u8::try_from(rank).expect("a line crosses a cylinder twice at most"),
                }
            }
        })
    }
}

/// A seam vertex — a point on both `∂A` and `∂B` (a three-plane corner or a pierce point), shared
/// (one `Handle`) by every incident result piece.
///
/// **It carries the cache its vertex receives, and nothing measured beside it.** `Realized`
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
/// **The walls are carried, not derived.** Reading an edge's supporting plane back out of its two
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
/// Keys live in *handle* space, the space `edge_of` keys. A **straight** piece is its
/// unordered endpoint pair — a plane-pair line and a ruling alike: two points bound one straight
/// segment, so a line one face states in the plane vocabulary and the face across states as a
/// ruling (a tangent wall ending on a cylinder's ruling) is one edge, and its carriers are read
/// off the faces that use it, not off the key. An arc is the curve its cylinder and plane meet in,
/// between two ends, so the key names all four: the **plane** — a cap's circle and an oblique
/// class's section through the same two pierce vertices are two curves, and so are two sections —
/// and the ends **in CCW order about the axis** (`from → to`), because between one pair of ends a
/// curve offers two complementary pieces, which get the two orders. The minted edge stores its
/// vertices in that same order, which is what the `[A, B]`-CCW convention means downstream
/// (`Model::derive_edge_curve`'s circle arm). Two faces that state one arc on different planes
/// would mint two edges, each used once, and the closed-shell guard would name the open shell.
/// The plane is locked by a tilted cylinder against a box (`every_kind_with_an_oblique_cylinder_builds`):
/// the box's floor cuts the cylinder's cap, and the cap's arc and the floor's section run between
/// one pair of points.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum EdgeKey {
    Line((usize, usize)),
    Arc {
        cyl: usize,
        plane: usize,
        from: usize,
        to: usize,
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

    /// This ring's edges, ready for the exact predicates — [`combinatorics::ring_edges`] of its
    /// nodes and walls; `p` is the plane class the ring lies in, `None` on a lateral.
    pub(crate) fn edges(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        p: Option<usize>,
    ) -> Result<Vec<combinatorics::RingEdge>, BoolError> {
        combinatorics::ring_edges(jd, cyls, &self.nodes, &self.walls, p)
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
    /// Which class table this face's surface lives in — a plane class, or a cylinder class whose
    /// band this face is. Plane-only consumers project with [`ClassIx::plane`], whose panic is an
    /// upstream-filter-bug detector.
    pub(crate) surf: crate::planes::ClassIx,
    pub(crate) outer: Bound,
    /// Hole boundaries, each already wound so the kept material stays on its left
    /// about the face's outward normal.
    pub(crate) inner: Vec<Bound>,
    pub(crate) flip: bool,
}

impl LocalFace {
    /// Every **polygon** ring of the face, outer first. Circle bounds carry no nodes and are
    /// deliberately absent — node-level consumers (the seam table, edge scans) have nothing to
    /// read off them. A reader of the face's **region** takes every bound (`comp_face`).
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

/// **The material about one line, as the sectors between the half-faces that leave it** — the
/// one answer to «does a result face run on across this line, or end there» and «does the kept
/// material meet itself only along it».
///
/// Round a line where faces meet, each face leaves along one half-plane per side of the line it
/// lies on — a *spoke* — and the material between two neighbouring spokes is one *sector*, in
/// each operand or not. This is the radial order of the faces about an edge (CGAL `Nef_3`'s sphere
/// map, Mäntylä's sectors on a shared line); OCCT's `GetFaceOff` picks the next face by the same
/// angle. Only [`keep`] decides a sector, so the operation enters here and nowhere else.
///
/// - **A face runs on across the line exactly when every sector between its two spokes, on the
///   side it keeps, is kept** ([`Self::continues`]) — the material it bounds passes round the line
///   from one half to the other. A face with one spoke — the line is its edge — ends there.
/// - **The kept material meets itself only along the line when its sectors fall into two runs or
///   more** ([`Self::kept_runs`]).
///
/// Spokes are counted counter-clockwise from any one of them; `kept[i]` is the sector
/// counter-clockwise after spoke `i`.
#[derive(Clone, Debug)]
pub(crate) struct LineStar {
    /// The face each spoke is a half of — two spokes of one face carry one key.
    faces: Vec<usize>,
    kept: Vec<bool>,
}

/// One half-face leaving a line, for [`LineStar::from_spokes`]: its face, its solid, and whether
/// that solid's material lies counter-clockwise after it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Spoke {
    pub(crate) face: usize,
    pub(crate) solid: crate::planes::SolidSide,
    pub(crate) material_ccw: bool,
}

impl LineStar {
    /// A star whose sectors' keep the caller has already read (a class's labels).
    pub(crate) fn new(faces: Vec<usize>, kept: Vec<bool>) -> Self {
        debug_assert_eq!(faces.len(), kept.len(), "one sector after each spoke");
        Self { faces, kept }
    }

    /// **The sectors' `(in_a, in_b)`, walked from the spokes** — a sector lies in solid `X`
    /// exactly as the nearest spoke of `X` clockwise before it says. `None` unless each solid
    /// brings exactly two spokes whose material sides face each other: about a line, a valid
    /// solid is one face running through it or two faces ending on it, and its material is the
    /// sectors on one side of that pair.
    pub(crate) fn sectors_of(spokes: &[Spoke]) -> Option<Vec<(bool, bool)>> {
        use crate::planes::SolidSide;
        let n = spokes.len();
        let mut read = [None; 2];
        for (slot, solid) in read.iter_mut().zip([SolidSide::A, SolidSide::B]) {
            let own: Vec<usize> = (0..n).filter(|&i| spokes[i].solid == solid).collect();
            let [p, q] = own[..] else { return None };
            if spokes[p].material_ccw == spokes[q].material_ccw {
                return None;
            }
            // Sector `i` lies after `p` and before `q` (counter-clockwise) or after `q`.
            *slot = Some(move |i: usize| {
                if (p..q).contains(&i) {
                    spokes[p].material_ccw
                } else {
                    spokes[q].material_ccw
                }
            });
        }
        let [Some(a), Some(b)] = read else {
            return None;
        };
        Some((0..n).map(|i| (a(i), b(i))).collect())
    }

    /// The star of `spokes` under `kind` — `None` where [`Self::sectors_of`] is.
    pub(crate) fn from_spokes(kind: BoolKind, spokes: &[Spoke]) -> Option<Self> {
        let sectors = Self::sectors_of(spokes)?;
        Some(Self::new(
            spokes.iter().map(|s| s.face).collect(),
            sectors.iter().map(|&(a, b)| keep(kind, a, b)).collect(),
        ))
    }

    /// **Does `face` run on across the line in the result?** `None` where no spoke of it is in the
    /// result (the sectors on its two sides are kept alike) or it has no spoke here; `Some(false)`
    /// where it ends — the line is its edge (one spoke), only one of its two spokes is in the
    /// result, or a sector between them on the side it keeps is not kept; `Some(true)` where the
    /// material it bounds passes round the line from one spoke to the other.
    pub(crate) fn continues(&self, face: usize) -> Option<bool> {
        let n = self.kept.len();
        let own: Vec<usize> = (0..n).filter(|&i| self.faces[i] == face).collect();
        let in_result = |i: usize| self.kept[(i + n - 1) % n] != self.kept[i];
        match own[..] {
            [p] => in_result(p).then_some(false),
            [p, q] => match (in_result(p), in_result(q)) {
                (false, false) => None,
                (true, true) => {
                    // The side `p` keeps runs counter-clockwise to `q`, or the other way round.
                    let arc: Vec<usize> = if self.kept[p] {
                        (p..q).collect()
                    } else {
                        (q..p + n).map(|i| i % n).collect()
                    };
                    Some(arc.iter().all(|&i| self.kept[i]))
                }
                _ => Some(false),
            },
            _ => {
                debug_assert!(
                    own.is_empty(),
                    "a face leaves a line along two half-planes at most"
                );
                None
            }
        }
    }

    /// How many runs the kept sectors fall into round the line — `0` with none kept, `1` with all.
    pub(crate) fn kept_runs(&self) -> usize {
        let n = self.kept.len();
        if self.kept.iter().all(|&k| k) {
            return usize::from(n > 0);
        }
        (0..n)
            .filter(|&i| self.kept[i] && !self.kept[(i + n - 1) % n])
            .count()
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
/// table the lateral chart reads its arc labels against and cuts its rims at
/// (`Curved::split_rims`). It is not the result's rim: the cleaning pass merges cap cells back and
/// dissolves the nodes no ring turns at, so what the assembly meets the result's edges on is
/// [`HeldRims`].
pub(crate) type CutRims = HashMap<(usize, usize), CutRim>;

/// **The result's cut rims — the split's nodes that the cleaned faces still hold.** Presence here
/// is the one source of "this rim is cut" for the assembly (per solid, its rim skip —
/// [`HeldRims::cut_per_group`]). A type of its own, apart from [`CutRims`]: with two same-typed
/// tables, the one read in the wrong place would let a Cut whose tool rests on a cap across its
/// rim merge the cap back into a whole circle while the assembly, reading the split, would still
/// take the rim for cut at the tool's planes.
#[derive(Clone, Debug, Default)]
pub(crate) struct HeldRims(CutRims);

impl HeldRims {
    /// The held nodes of one rim — the tests' reading; the assembly asks per group
    /// ([`Self::cut_per_group`]).
    #[cfg(test)]
    pub(crate) fn get(&self, key: &(usize, usize)) -> Option<&CutRim> {
        self.0.get(key)
    }

    /// **Which result solid still cuts which rim** — `(group, cyl, plane)` for every rim of this
    /// table some face of that group on the rim's plane still holds a node of, read off a later
    /// face list. For the assembly, after the per-solid straight-angle pass: that pass decides
    /// per solid, so one body can dissolve a rim to its whole circle while another, touching it,
    /// keeps a node of it as a corner of its own. Each body's handles are its own, so the rim is
    /// whole for the first and cut for the second, and asking with the group is what lets the
    /// two answers stand side by side. One held node is enough, as in [`held_rims`].
    pub(crate) fn cut_per_group(
        &self,
        faces: &[LocalFace],
        group_of: &[usize],
    ) -> std::collections::HashSet<(usize, usize, usize)> {
        let mut held: std::collections::HashSet<(usize, usize, NodeId)> =
            std::collections::HashSet::new();
        for (lf, &g) in faces.iter().zip(group_of) {
            let ClassIx::Plane(c) = lf.surf else {
                continue;
            };
            for ring in lf.poly_rings() {
                held.extend(ring.nodes.iter().map(|&n| (g, c, n)));
            }
        }
        let groups: std::collections::BTreeSet<usize> = group_of.iter().copied().collect();
        let mut out = std::collections::HashSet::new();
        for (&(cyl, c), cr) in &self.0 {
            for &g in &groups {
                if cr.nodes.iter().any(|&n| held.contains(&(g, c, n))) {
                    out.insert((g, cyl, c));
                }
            }
        }
        out
    }
}

/// [`HeldRims`] from the faces as they stand: for each split circle `(cyl, c)`, the nodes some
/// face on plane class `c` still has in a ring, in the split's θ order. A circle none of whose
/// nodes survive is **whole** again and
/// leaves the table — the cap holds it as a circle and the lateral, cleaned by the same pass, as
/// its closed rim.
///
/// Any node of a class-`c` ring, not only an arc's ends: the inward wedge with its apex on the seam
/// has the cap turn at the apex between two lines, and an arc-ends rule would drop that node — the
/// lateral's corner there — moving that `A − B`'s refusal to `ArcBoundNotYet` (measured).
///
/// A circle left **one** node — the split cut it at one point, or the cleaning dissolved all its
/// others — is cut at that one node: the cap holds it as the arc from it back to itself
/// (`dissolve_straight_angles`), and so does the lateral's ring, cleaned by the same pass.
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
        if nodes.is_empty() {
            continue;
        }
        out.insert((cyl, c), CutRim { nodes });
    }
    HeldRims(out)
}

#[cfg(test)]
#[path = "tests/draft.rs"]
mod tests;
