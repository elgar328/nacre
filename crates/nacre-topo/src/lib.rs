//! b-rep topology and the truth-only `Model` aggregate (design.md §2, §4).
//!
//! The topology (`Vertex`/`Edge`/`Face`/`Loop`/`HalfEdge`/`Shell`/`Solid`)
//! references exact geometry only by `Handle` — geometry never knows about
//! topology, topology never inspects coordinates (overview 절대원칙 2).
//!
//! [`Model`] is the **truth**: exact geometry stores + topology stores + the
//! derived [`Adjacency`] cache. It holds no tessellation and no operation log —
//! a mesh cache and the op log are companions owned at higher layers (§0: "the
//! model is the replay result"; putting them here would make topo depend on
//! tess/ops and break the truth/cache split).

mod adjacency;
mod topology;

pub use adjacency::{Adjacency, nonmanifold_vertices};
pub use topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

use nacre_geom::{Circle, Curve, Cylinder, Line, Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Angle, Axis, Rat};
use nacre_store::{Handle, Store};
use std::collections::{HashMap, HashSet};

/// How a discovered vertex is *defined* — the primitives whose intersection it
/// is (design §4). This definition is the **truth**; the vertex's `f64` point is
/// a cache derived from it. Sign decisions (in/out, orientation) feed this
/// definition to the indirect predicates rather than the cached coordinate
/// (design §3, §8 M5), so a discovered vertex must carry it.
///
/// M5 polyhedral vertices are three-plane intersections; later milestones add
/// variants (a line∩plane point, quadric intersections). `Handle<Surface>` is
/// `Copy` regardless of `Surface`, so this stays `Copy`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VertexDef {
    /// The intersection of three planes (their surface handles).
    ThreePlane([Handle<Surface>; 3]),
}

/// One motion in a history — what a [`MotionNode`] carries.
///
/// Motions do **not** commute, so a history is one ordered chain, never a store per kind: "turn
/// then place" and "place then turn" are different motions and must stay tellable apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Motion {
    /// One axis-aligned rotation about the line through the rational pivot `point`.
    Rotate {
        axis: Axis,
        point: [Rat; 3],
        angle: Angle,
    },
    /// One exact rational translation.
    Translate { offset: [Rat; 3] },
    /// One reflection in the coordinate plane `axis = offset` (`x ↦ 2·offset − x` on that axis).
    ///
    /// **Improper** — the only motion here with `det = −1`. It preserves lengths, distances and
    /// incidence like the others, but it *negates* every determinant of its images rather than
    /// leaving them alone. Judgments that answer a determinant question in a shared pre-motion
    /// frame must account for that; see `nacre-ops`' `BaseFrame` and `nacre-cip`'s `shared_base`.
    Mirror { axis: Axis, offset: Rat },
    /// A change of basis **into a plane's own frame** — what makes a sketch on a tilted face
    /// exact. Coordinates written against this node are read as `(u, v, w)` in that plane's
    /// frame; the motion carries them out into the frame the plane itself lives in.
    ///
    /// ★★★★ **The normal is named, not spelled.** A wall raised on a tilted face has an
    /// irrational world normal — there is no `[Rat; 3]` for it — so the node points at the
    /// *plane*, whose coefficients are rational **in its own frame**. The recursion terminates
    /// because it walks down to a plane with no frame, which is where the world is.
    ///
    /// ★★★ **The frame is a function of the plane, so `plane` is the whole node.** Both of the
    /// other things a frame needs are derived and neither is stored:
    ///
    /// ```text
    /// origin  = the world origin projected onto the plane
    /// ref_dir = ẑ × n   (ŷ × n when n is vertical — the arbitrary-axis branch)
    /// ```
    ///
    /// That is what makes canonicalization unnecessary: one plane has one frame node, so a handle
    /// is already the key. Widening to a user-named frame later means adding the fields back and
    /// treating today's as the specialization.
    ///
    /// ★★★ **`flip` is what makes the node a whole frame and not half of one.** Canonical plane
    /// coefficients carry **no direction** — the first nonzero component is forced positive,
    /// because their question is *"are these the same plane"*. A frame's `ŵ` is a direction, and
    /// two faces of one plane can face opposite ways; `flip` says to negate the coefficients, so
    /// the node names the sense as well as the plane. Without it a sketch on a reversed face comes
    /// out mirrored in `u` with its sweep running inward (measured: a tilted second boss came back
    /// `PadMissesFace`).
    ///
    /// It costs nothing that mattered: the frame is still a pure function of what the node holds,
    /// so it is still reconstructible from the handle, and two sketches on **one face** still
    /// intern to one node — which is the sharing the design needs.
    ///
    /// ★★★ **`origin` and `ref_dir` are named, not derived.** A *face* has no name, so its frame
    /// is derived by convention (`nacre_scalar::plane_frame_default` — the world origin projected
    /// onto the plane, and `ẑ × n`). A plane a caller **named** can insist on something no
    /// derivation produces: the script layer's `ZX` has `+u = +ẑ`, while `ẑ × n` there is `−x̂`.
    /// So the node carries both, and the derived case is the one that fills them in from the
    /// convention. `ref_dir` lies in the plane, need not be unit length, and is reduced to its
    /// primitive form so two spellings of one direction name one node.
    ///
    /// ★ **Proper** (`det = +1`) — `(u, v, w)` is right-handed by construction (`v = w × u`), so
    /// unlike [`Mirror`](Self::Mirror) it contributes nothing to a chain's parity. Negating the
    /// coefficients flips `ŵ` and `û` together and leaves `v̂`, which is a half-turn: still proper.
    Frame {
        plane: Handle<Surface>,
        origin: [Rat; 3],
        ref_dir: [Rat; 3],
        flip: bool,
    },
}

/// A node in the motion-history forest (design §CIP ⑦): one [`Motion`] applied to a solid, with a
/// parent link so several points can share a history's tail.
///
/// Stored in [`Model::motions`]; a moved vertex's [`Origin::Moved`] and a moved surface's
/// [`SurfaceDef::Moved`] name their leaf node. The tol a motion contributes is
/// application-point-dependent, so it is **not** stored here — judgment computes it by traversing
/// to the root.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MotionNode {
    pub motion: Motion,
    pub parent: Option<Handle<MotionNode>>,
}

/// Provenance of a vertex or edge (design §4, overview 절대원칙 4).
///
/// `Constructed` elements know their identity by construction and carry no
/// tolerance; `Discovered` elements come from an intersection and hold both the
/// [`VertexDef`] that defines them (the truth) and the *measured* accuracy the
/// relaxation/closed-form achieved (`tol` — the point is a within-`tol` cache of
/// the definition). In M1–M4 every element is `Constructed`; M5's
/// `PolyhedralBoolean` is the first `Discovered` producer. `Rotated` (overhaul stage
/// 1b) names a vertex that is another vertex (`base`) moved by a motion node — its
/// point is a cache; the tol is judgment-time (§CIP ⑦), so no tol slot here.
///
/// Holds an `f64`, so `PartialEq` only — no `Eq`/`Hash` (identity is by
/// `Handle`, never by value).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Origin {
    Constructed,
    Discovered {
        tol: f64,
        definition: VertexDef,
    },
    Moved {
        base: Handle<Vertex>,
        motion: Handle<MotionNode>,
    },
}

/// What makes two surfaces the same plane, for [`Model::surface_ids`]: the canonical rational
/// coefficients and **the motion they are stated in**. Identical arrays under different motions
/// are different planes, because `Constructed` coefficients speak about the world and `Moved`
/// ones about the pre-motion frame.
pub type SurfaceKey = ([nacre_scalar::Rat; 4], Option<Handle<MotionNode>>);

/// Provenance of a **surface** — the same question [`Origin`] answers for a vertex.
///
/// The kernel's rule is that exact geometry is the truth and f64 is a cache. A surface's
/// coefficients are that truth only while nothing irrational has been applied to them: a
/// non-90° rotation turns a plane's normal into an irrational direction, and the stored
/// coefficients become a *rounded copy*. Without this, the kernel has no way to tell the two
/// apart, so it treats the copy as exact — and two rounded copies of one wall then fail to be
/// the same plane, which is how a chained boolean loses a coplanar contact.
///
/// Recorded through [`Model::push_surface`] and enforced by `nacre-validate`: every face of a
/// live solid must lie on a surface that has a definition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SurfaceDef {
    /// The coefficients **are** the truth — a constructed surface, or one moved by a motion
    /// that preserves exactness (a translation, a 90°-family rotation).
    Constructed,
    /// The image of an exact surface under a motion history. The truth is `(witness, motion)`;
    /// the coefficients are a cache.
    ///
    /// `motion` is the leaf of the history in [`Model::motions`] — the chain lives in the forest,
    /// exactly as [`Origin::Moved`] uses it. **What is moved** is the plane's own three points in
    /// [`Model::surface_points`], stated in the pre-motion frame.
    ///
    /// ★★★★★ **This used to carry a `witness: [Point3; 3]` and it was the last place a rounded
    /// coordinate served as a definition.** The doc claimed those three points were "exactly
    /// representable in f64"; measured, a third of them were not — `11/10` is not an f64, so
    /// lifting the realization back recovered a different rational and the judge described the
    /// plane through three rounded points. A moved plane whose source has no exact points is now
    /// [`Inexact`](Self::Inexact), which is what it always was.
    Moved { motion: Handle<MotionNode> },
    /// **Not exactly describable** — a history the forest cannot express. The kernel does not
    /// pretend the coefficients are exact; consumers that need the truth reject honestly.
    Inexact,
}

/// Whether a face uses its surface normal as-is (`Forward`) or flipped
/// (`Reversed`). A pure tag — full derives.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reversed,
}

impl Orientation {
    /// The opposite orientation.
    #[inline]
    pub fn flipped(self) -> Orientation {
        match self {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
        }
    }
}

/// The truth-only aggregate: exact geometry + topology stores + the derived
/// adjacency cache. No `tess`, no `ops` (see the crate docs).
#[derive(Debug, Default)]
pub struct Model {
    // exact geometry (truth)
    pub surfaces: Store<Surface>,
    pub curves: Store<Curve>,
    /// The motion-history forest (§CIP ⑦): motion definitions named by `Origin::Moved` vertices
    /// and `SurfaceDef::Moved` surfaces. Not geometry — a definition store.
    ///
    /// **Interned** — see [`Model::push_motion`]; write through it, never through `Store::push`.
    pub motions: Store<MotionNode>,
    /// Interning table for [`Model::push_motion`]: the handle already issued for a given
    /// `(motion, parent)`. Not iterated (a `HashMap`'s order must never reach a result).
    pub motion_ids: HashMap<MotionNode, Handle<MotionNode>>,
    /// Each surface's provenance, keyed by handle. Not geometry — the twin of the vertex
    /// `Origin`, kept beside the store because [`Surface`] is a `nacre-geom` type and cannot
    /// name a `Handle<MotionNode>`. Written only by [`Model::push_surface`].
    ///
    /// Iterate this **through the faces**, never over the map: a `HashMap`'s order is not
    /// deterministic and replay determinism (DNA 3) forbids letting it reach a result.
    pub surface_defs: HashMap<Handle<Surface>, SurfaceDef>,
    /// Each surface's plane as **exact rational coefficients**, canonicalized — present only for
    /// surfaces whose producer had a rational description to record.
    ///
    /// ★ **The point is that two faces of one plane get the same array.** `nacre_geom::Plane`
    /// keeps an un-normalized normal whose length follows the *face's size*, so the same plane
    /// reaches `coefficients()` as `[2.2, 0, 0, −6.6000000000000005]` from one face and
    /// `[13.2, 0, 0, −39.599999999999994]` from another — not exactly proportional, because `d`
    /// is a rounded product. Canonical rationals have no scale to disagree about.
    ///
    /// ★★ **These are built from the dimensions the user wrote, never lifted from the f64
    /// coefficients above.** Lifting is lossless and useless here: it preserves the rounding, so
    /// the two vectors stay different (`nacre_scalar::canonical_plane_coeffs`).
    ///
    /// ★★★ **The coefficients are stated in the frame this surface's [`SurfaceDef`] names** —
    /// the world for [`SurfaceDef::Constructed`], and the **pre-motion** frame for
    /// [`SurfaceDef::Moved`], whose world coefficients are irrational and so cannot be written
    /// down at all. A moved surface therefore inherits its source's array unchanged: the motion
    /// is recorded beside it, not folded into it. (That pairing — exact coefficients plus a
    /// motion — is the shape `docs/truth-and-cache.md` builds toward.)
    ///
    /// Absent is ordinary — an f64 construction path, or `i128` overflow.
    /// Iterate through the faces, never over the map.
    pub surface_coeffs: HashMap<Handle<Surface>, [nacre_scalar::Rat; 4]>,
    /// ★★★★★ **Three exact points the plane passes through** — stated in the same frame as
    /// [`Model::surface_coeffs`], and carried out by the same [`SurfaceDef`] motion.
    ///
    /// **The coefficients cannot be the whole truth, and this is why.** A plane's exact
    /// coefficients are a *product* of two point differences, so they need about twice the bits
    /// the points do — measured, a wall plane's coefficients run to a median of 159 bits and a
    /// maximum of 213, and only 25.8% fit `i128`, while the ring coordinates they come from top
    /// out at **107** with 20 bits to spare. Points fit where their products do not, so the truth
    /// stores points and the judging layer — already arbitrary-precision — does the products.
    ///
    /// ★★★ **The recursion has a base case**: a boolean never creates a surface (asserted by
    /// `nacre-ops`), so a plane is never defined by discovered points, only by construction ones.
    /// Points are defined by planes and planes by construction points — two levels, not a cycle.
    ///
    /// ★★ **Interning keeps the first pusher's triple.** Any non-collinear triple on the plane
    /// names the same plane, so the survivor is as good a witness as the newcomer; and the
    /// separate invariant that a *face's* vertices lie on its own plane is carried by
    /// `outer_tri`, not by this.
    ///
    /// Absent is ordinary while the producers are still being wired, and permanently for a prism
    /// whose ring is not rational at all. Iterate through the faces, never over the map.
    pub surface_points: HashMap<Handle<Surface>, [[nacre_scalar::Rat; 3]; 3]>,
    /// Interning table for [`Model::push_surface_with_coeffs`]: the handle already issued for a
    /// plane, keyed by its canonical coefficients **and the motion they are stated in**. The twin
    /// of [`Model::motion_ids`]; not iterated (a `HashMap`'s order must never reach a result).
    ///
    /// ★ **The motion belongs in the key.** `Constructed` coefficients speak about the world and
    /// `Moved` ones about the pre-motion frame, so two identical arrays under different motions
    /// are different planes. `Inexact` has no coefficients and never interns.
    ///
    /// ★★ **The witness does not.** Two faces of one plane sharing one `SurfaceDef::Moved`
    /// witness is already how this works — `collect_planes` says the witness "was captured from
    /// whichever face first reached this surface" and winds it to *each* face's own outward
    /// normal — so it is not part of what makes two planes the same.
    pub surface_ids: HashMap<SurfaceKey, Handle<Surface>>,
    // topology (references geometry by Handle only)
    pub vertices: Store<Vertex>,
    pub edges: Store<Edge>,
    pub faces: Store<Face>,
    pub shells: Store<Shell>,
    pub solids: Store<Solid>,
    /// The live solids — the "current model" (design §2 supersede semantics).
    /// Editing ops supersede topology by pushing new cells and updating this
    /// list; the old cells stay in the append-only arena but, unreferenced by
    /// any live solid, drop out of the reachable closure. Producers register
    /// through [`Model::push_solid`]; `validate`/`Adjacency`/`nacre-step`
    /// traverse [`Model::reachable`], not the whole store.
    pub live_solids: Vec<Handle<Solid>>,
    // derived cache (rebuilt on demand)
    pub adj: Adjacency,
}

/// The handles reachable from a model's live solids — the live model (design §2).
///
/// Only the sets consumers need today: `validate`'s Euler counts vertices/edges/
/// faces/shells, and `Adjacency`/loop/incidence walk faces/edges. Surfaces and
/// curves are never orphaned by the M4 face ops, so they are not tracked.
#[derive(Debug, Default)]
pub struct Reachable {
    pub vertices: HashSet<Handle<Vertex>>,
    pub edges: HashSet<Handle<Edge>>,
    pub faces: HashSet<Handle<Face>>,
    pub shells: HashSet<Handle<Shell>>,
}

/// Whether a handle indexes inside its store (guards the reachable traversal
/// against dangling/out-of-range handles).
#[inline]
fn in_bounds<T>(h: Handle<T>, store: &Store<T>) -> bool {
    (h.index() as usize) < store.len()
}

impl Loop {
    /// This loop with its winding reversed: half-edges in reverse order, each traversed the
    /// opposite way. Both parts matter — reversing the *order* flips the normal the winding
    /// implies, and flipping each `forward` keeps every edge used once in each direction, so two
    /// faces that share an edge stay opposed when both are reversed (still a valid 2-manifold).
    ///
    /// Two callers want it for opposite reasons. [`Model::reversed_shell`] pairs it with an
    /// [`Orientation`] toggle to turn a boundary inward (a cavity). A **reflection** pairs it with
    /// nothing: mirroring negates the normal a winding implies, so rewinding restores it and the
    /// orientation flag stays as it was.
    pub fn reversed(&self) -> Loop {
        Loop {
            half_edges: self
                .half_edges
                .iter()
                .rev()
                .map(|he| HalfEdge {
                    edge: he.edge,
                    forward: !he.forward,
                })
                .collect(),
        }
    }
}

impl Model {
    /// An empty model.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Recompute the [`Adjacency`] cache from the current topology stores.
    /// Call once after a batch of additions (the cache is otherwise stale).
    pub fn rebuild_adjacency(&mut self) {
        // Build against an immutable borrow, then move into place — avoids
        // borrowing `self.adj` mutably while iterating the other stores.
        let adj = Adjacency::rebuild(&*self);
        self.adj = adj;
    }

    /// Push a solid into the store **and mark it live** (design §2). This is the
    /// blessed way for a producer to add a solid; the reachable closure
    /// ([`Model::reachable`]) grows to include it. Editing ops instead mutate
    /// [`Model::live_solids`] directly (drop the superseded solid, add the new).
    pub fn push_solid(&mut self, solid: Solid) -> Handle<Solid> {
        let h = self.solids.push(solid);
        self.live_solids.push(h);
        h
    }

    /// Push a motion node, **interned**: the same `(motion, parent)` always yields the same
    /// handle.
    ///
    /// The handle is the canonical name of "which motion", and judgments use it to decide whether
    /// a set of points shares one rigid motion — which lets the whole judgement be answered
    /// exactly in the pre-motion frame. That identity has to be *both* collision-free and free of
    /// false misses:
    ///
    /// - it used to be a 64-bit hash of the chain's contents, and a collision would hand the exact
    ///   predicate two incompatible frames and answer a different question with full confidence;
    /// - a raw `Store::push` per transform is collision-free but *misses*: turning two solids by
    ///   the same 30° would make two nodes, and their shared motion would stop cancelling — so
    ///   rotating a model would turn its exact questions into assumed ones (measured: it breaks
    ///   `a_shared_rotation_still_assumes_nothing`).
    ///
    /// Interning gives both. It also keeps the forest small, since a chain shared by many solids
    /// is stored once.
    pub fn push_motion(
        &mut self,
        motion: Motion,
        parent: Option<Handle<MotionNode>>,
    ) -> Handle<MotionNode> {
        let node = MotionNode { motion, parent };
        if let Some(&h) = self.motion_ids.get(&node) {
            return h;
        }
        let h = self.motions.push(node);
        self.motion_ids.insert(node, h);
        h
    }

    /// Push a surface **and state its provenance** ([`SurfaceDef`]) — the blessed way to add
    /// one. Pushing into `surfaces` directly leaves the surface undefined, which
    /// `nacre-validate` reports for any face that lands on it.
    pub fn push_surface(&mut self, surface: Surface, def: SurfaceDef) -> Handle<Surface> {
        let h = self.surfaces.push(surface);
        self.surface_defs.insert(h, def);
        h
    }

    /// [`Model::push_surface`], also recording the surface's exact rational description when the
    /// producer has it: the canonical coefficients ([`Model::surface_coeffs`]) and the three
    /// points the plane passes through ([`Model::surface_points`]). `None` records nothing, which
    /// is what a producer without a rational description passes.
    ///
    /// ★ **Both are stated in the frame this `def` names**, and an interned plane keeps whichever
    /// the first pusher supplied.
    ///
    /// ★ **The `bool` says the returned surface's normal points the *other* way** from the one
    /// handed in, and a caller that meets it must record its face `Orientation::flipped()`.
    /// It is always `false` today, because every call mints a fresh handle. It exists because a
    /// plane's canonical form has no direction — `[0,0,1,−3]` and `[0,0,−1,3]` are one plane —
    /// so once identical planes start sharing a handle the direction has to be reconciled
    /// somewhere, and the honest place is where the caller still knows what it asked for.
    ///
    /// ★★ **The relative answer, not an absolute one.** A caller could instead read the returned
    /// surface's normal and derive its orientation from its own outward direction, which is what
    /// `extrude`'s pad-on-face path already does — but `transform_solid` *copies* faces without
    /// recomputing an outward normal, so it has nothing to compare against. "Is it flipped from
    /// what I asked for?" is the question every caller can answer.
    pub fn push_surface_with_coeffs(
        &mut self,
        surface: Surface,
        def: SurfaceDef,
        coeffs: Option<[nacre_scalar::Rat; 4]>,
        points: Option<[[nacre_scalar::Rat; 3]; 3]>,
    ) -> (Handle<Surface>, bool) {
        // ★★★★★ **A plane may not be stated two ways.** The points are the truth; the
        // coefficients are its *name*, and a name that does not fit the thing is not issued.
        //
        // ★ **Not a `debug_assert`.** A producer that disagreed would then be caught in tests and
        // waved through in release — which is exactly the silence this kernel keeps closing. The
        // established answer is the one used for `i128` overflow everywhere else: **demote, do not
        // crash**. A plane with no name still has its points, and every consumer of the name
        // already handles its absence.
        //
        // ★★ **Inconclusive counts as "no".** `c · p` can overflow `Rat` when the coefficients are
        // wide (measured 774 of 73,369 — the same "wide normal" population `plane_frame_named`
        // already declines), and a check that could not be evaluated is not a check that passed.
        //
        // ★★★ **Measured, and it caught one.** Across the suite: **0 violations**, and 1,546 of
        // 76,135 (2.0%) inconclusive — the wide-coefficient population, whose names it therefore
        // does not issue. That costs nothing observable: those planes already could not host a
        // frame (`plane_frame_named` declines on the same width) and none of them interned with
        // anything (census bit-identical).
        //
        // ★★★★★ The violation it did catch was real and was **introduced by the points work
        // itself**: a prism's base cap recorded the ring in *frame* coordinates while its
        // coefficients spoke about the world. Axis-aligned models never showed it, because there
        // the frame **is** the world. See `ops::named_plane_points`.
        let coeffs = coeffs.filter(|c| match points {
            None => true,
            Some(ps) => ps.iter().all(|p| {
                (|| {
                    let mut acc = c[3];
                    for k in 0..3 {
                        acc = acc.checked_add(c[k].checked_mul(p[k])?)?;
                    }
                    Some(acc == nacre_scalar::Rat::from_int(0))
                })() == Some(true)
            }),
        });
        let key = coeffs.map(|c| {
            let motion = match def {
                SurfaceDef::Moved { motion, .. } => Some(motion),
                _ => None,
            };
            (c, motion)
        });
        if let Some(k) = key {
            if let Some(&h) = self.surface_ids.get(&k) {
                // Same plane, already issued. The canonical form says nothing about direction, so
                // report whether the survivor points the other way and let the caller spell its
                // outward the other way round.
                let dir = |s: &Surface| match s {
                    Surface::Plane(p) => Some(p.normal()),
                    Surface::Cylinder(_) => None,
                };
                let flipped = match (dir(self.surfaces.get(h)), dir(&surface)) {
                    (Some(a), Some(b)) => a.dot(b) < 0.0,
                    _ => false,
                };
                return (h, flipped);
            }
        }
        let h = self.push_surface(surface, def);
        if let Some((c, _)) = key {
            self.surface_coeffs.insert(h, c);
        }
        if let Some(p) = points {
            self.surface_points.insert(h, p);
        }
        if let Some(k) = key {
            self.surface_ids.insert(k, h);
        }
        (h, false)
    }

    /// A new shell whose faces are copies of `src`'s with their outward normals
    /// flipped inward: every loop's winding is reversed and every face
    /// `orientation` is toggled. Pushes fresh [`Face`] cells and a fresh
    /// [`Shell`], but **reuses** `src`'s surfaces, edges, curves, and vertices —
    /// which stay valid handles after the source solid is superseded
    /// (append-only). This is the building block for a cavity (void) shell: the
    /// boundary of a solid whose interior becomes empty space (design §8 M5
    /// containment). The reversed winding keeps each shared edge used with
    /// opposed half-edges (a valid 2-manifold), and the toggled orientation makes
    /// the outward normal point into the void.
    pub fn reversed_shell(&mut self, src: Handle<Shell>) -> Handle<Shell> {
        let src_faces = self.shells.get(src).faces.clone();
        let faces: Vec<Handle<Face>> = src_faces
            .iter()
            .map(|&fh| {
                let face = self.faces.get(fh).clone();
                self.faces.push(Face {
                    surface: face.surface,
                    outer: face.outer.reversed(),
                    inner: face.inner.iter().map(Loop::reversed).collect(),
                    orientation: face.orientation.flipped(),
                })
            })
            .collect();
        self.shells.push(Shell { faces })
    }

    /// The vertex a half-edge starts at: its edge's `bounds[0]` when the use runs
    /// forward, `bounds[1]` when it runs back.
    ///
    /// A traversal accessor, not an analysis — the same kind of thing as
    /// [`Model::reachable`], and the reason it lives here: walking a loop's
    /// corners is the first thing every consumer above does, and it was written
    /// twice (with two different failure policies) before this existed.
    ///
    /// `None` for an edge with no endpoints — the standalone full circle of §4,
    /// which is a legitimate form but has no start. Callers pick their own
    /// policy: a solid's loop edge is always bounded, so `nacre-ops` unwraps with
    /// that invariant while `nacre-props` reports it as unsupported input.
    #[inline]
    pub fn he_start(&self, he: HalfEdge) -> Option<Handle<Vertex>> {
        let [a, b] = self.edges.get(he.edge).bounds?;
        Some(if he.forward { a } else { b })
    }

    /// The handles reachable from the live solids — the live model (design §2).
    ///
    /// Superseded cells left in the append-only arena are excluded (nothing live
    /// references them). Every step is bounds-guarded, so this is safe even on a
    /// corrupt or partially-built model (a dangling handle simply prunes that
    /// branch; `validate`'s reference-integrity check reports it separately).
    pub fn reachable(&self) -> Reachable {
        let mut r = Reachable::default();
        for &solid_h in &self.live_solids {
            if !in_bounds(solid_h, &self.solids) {
                continue;
            }
            let solid = self.solids.get(solid_h);
            for &shell_h in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
                if !in_bounds(shell_h, &self.shells) || !r.shells.insert(shell_h) {
                    continue;
                }
                for &face_h in &self.shells.get(shell_h).faces {
                    if !in_bounds(face_h, &self.faces) || !r.faces.insert(face_h) {
                        continue;
                    }
                    let face = self.faces.get(face_h);
                    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                        for he in &lp.half_edges {
                            if !in_bounds(he.edge, &self.edges) || !r.edges.insert(he.edge) {
                                continue;
                            }
                            if let Some(bounds) = self.edges.get(he.edge).bounds {
                                for v in bounds {
                                    if in_bounds(v, &self.vertices) {
                                        r.vertices.insert(v);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        r
    }

    /// Add an axis-aligned box `min`..`max` to this model and return its solid.
    ///
    /// Requires `max[i] > min[i]` on every axis (a degenerate box is a caller
    /// bug → panic). Every element is [`Origin::Constructed`] (design §8); each
    /// face is wound so its plane normal points outward, so all faces are
    /// [`Orientation::Forward`]. Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cuboid(&mut self, min: Point3, max: Point3) -> Handle<Solid> {
        let [x0, y0, z0] = min.as_array();
        let [x1, y1, z1] = max.as_array();
        debug_assert!(
            x1 > x0 && y1 > y0 && z1 > z0,
            "cuboid needs max > min on every axis"
        );

        // 8 corners, numbered by bits (i·x, j·y, k·z).
        let corners = [
            Point3::from_array([x0, y0, z0]), // V0
            Point3::from_array([x1, y0, z0]), // V1
            Point3::from_array([x1, y1, z0]), // V2
            Point3::from_array([x0, y1, z0]), // V3
            Point3::from_array([x0, y0, z1]), // V4
            Point3::from_array([x1, y0, z1]), // V5
            Point3::from_array([x1, y1, z1]), // V6
            Point3::from_array([x0, y1, z1]), // V7
        ];
        // 6 faces: (first 3 loop vertices for the outward plane, half-edges as
        // (edge index, forward)). Loops wound CCW seen from outside → outward
        // normal. See the design plan's winding table (hand-verified).
        type FaceDef = ([usize; 3], [(usize, bool); 4]);
        let faces_def: [FaceDef; 6] = [
            ([0, 3, 2], [(3, false), (2, false), (1, false), (0, false)]), // Bottom −Z
            ([4, 5, 6], [(4, true), (5, true), (6, true), (7, true)]),     // Top +Z
            ([0, 1, 5], [(0, true), (9, true), (4, false), (8, false)]),   // Front −Y
            ([2, 3, 7], [(2, true), (11, true), (6, false), (10, false)]), // Back +Y
            ([0, 4, 7], [(8, true), (7, false), (11, false), (3, true)]),  // Left −X
            ([1, 2, 6], [(1, true), (10, true), (5, false), (9, false)]),  // Right +X
        ];
        // ★ **Surfaces before vertices**, so each corner can name the three it lies on. Separate
        // arenas, so the interleaving does not move any handle; the surfaces' order among
        // themselves is what matters and it is unchanged.
        let surf: [(Handle<Surface>, bool); 6] = core::array::from_fn(|i| {
            let (tri, _) = &faces_def[i];
            // The same three corners in rationals. `from_decimal` because a corner is a value the
            // caller *wrote* — lifting the f64 bit pattern instead would carry its drift in and
            // defeat the whole point (see `Model::surface_coeffs`).
            let rat_corner = |k: usize| -> Option<[nacre_scalar::Rat; 3]> {
                let c = corners[k].as_array();
                Some([
                    nacre_scalar::Rat::from_decimal(c[0])?,
                    nacre_scalar::Rat::from_decimal(c[1])?,
                    nacre_scalar::Rat::from_decimal(c[2])?,
                ])
            };
            let coeffs = (|| {
                nacre_scalar::plane_through_points(
                    rat_corner(tri[0])?,
                    rat_corner(tri[1])?,
                    rat_corner(tri[2])?,
                )
            })();
            self.push_surface_with_coeffs(
                Surface::Plane(
                    Plane::through_points(corners[tri[0]], corners[tri[1]], corners[tri[2]])
                        .expect("non-degenerate box"),
                ),
                SurfaceDef::Constructed,
                coeffs,
                // The same three corners, exactly — a box's corners are integers, so this is
                // recorded even where `coeffs` overflows (it cannot here, but the rule is the
                // rule: points are the truth, coefficients are the name).
                (|| {
                    Some([
                        rat_corner(tri[0])?,
                        rat_corner(tri[1])?,
                        rat_corner(tri[2])?,
                    ])
                })(),
            )
        });

        let vh: [Handle<Vertex>; 8] = core::array::from_fn(|i| {
            // Which three of the six faces meet at corner `i`. The corner order is
            // `V0..V3` round the bottom then `V4..V7` round the top, so the high bit picks the cap
            // and the position round the ring picks the two walls.
            let m = i % 4;
            let cap = if i < 4 { 0 } else { 1 }; // Bottom −Z / Top +Z
            let along_x = if m == 1 || m == 2 { 5 } else { 4 }; // Right +X / Left −X
            let along_y = if m >= 2 { 3 } else { 2 }; // Back +Y / Front −Y
            self.vertices.push(Vertex {
                point: corners[i],
                origin: Origin::Constructed,
                definition: Some(VertexDef::ThreePlane([
                    surf[cap].0,
                    surf[along_y].0,
                    surf[along_x].0,
                ])),
            })
        });

        // 12 edges as (start, end) vertex indices: bottom ring, top ring, verticals.
        const EDGES: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let eh: [Handle<Edge>; 12] = core::array::from_fn(|i| {
            let (a, b) = EDGES[i];
            let curve = self.curves.push(Curve::Line(
                Line::through_points(corners[a], corners[b]).expect("non-degenerate box"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([vh[a], vh[b]]),
                origin: Origin::Constructed,
            })
        });

        let fh: [Handle<Face>; 6] = core::array::from_fn(|i| {
            let (_, hes) = &faces_def[i];
            let (surface, flipped) = surf[i];
            let outer = Loop {
                half_edges: hes
                    .iter()
                    .map(|&(e, forward)| HalfEdge {
                        edge: eh[e],
                        forward,
                    })
                    .collect(),
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                // The winding table below is written for a surface whose normal this face uses as-is;
                // a shared surface may point the other way, and then the same outward direction is
                // spelled `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            })
        });

        let shell = self.shells.push(Shell { faces: fh.to_vec() });
        self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        })
    }

    /// Add a closed cylinder solid and return it: the axis runs from `base` along
    /// `axis` for `height`, with the given `radius`.
    ///
    /// Built as the standard **seam b-rep** (V2/E3/F3): two seam vertices, two
    /// full-circle rim edges (`bounds: Some([seam, seam])`, start == end), one
    /// straight seam edge, a cylindrical lateral face whose loop uses the seam
    /// edge twice (opposite orientation), and two planar caps (each a
    /// single-half-edge rim loop). This forms a valid CW-complex (Euler χ = 2)
    /// that [`validate`](../nacre_validate/fn.validate.html) accepts — a closed
    /// periodic surface needs a seam vertex, so the rims are `Some([v, v])`, not
    /// `bounds: None` (that form is for a standalone full circle; §4).
    ///
    /// `radius`/`height` must be positive and `axis` nonzero (caller bug →
    /// panic). Every element is [`Origin::Constructed`]. Does **not** rebuild
    /// adjacency — call [`Model::rebuild_adjacency`] once after all additions.
    pub fn add_cylinder(
        &mut self,
        base: Point3,
        axis: Vector3,
        radius: f64,
        height: f64,
    ) -> Handle<Solid> {
        debug_assert!(
            radius > 0.0 && height > 0.0,
            "cylinder needs positive radius and height"
        );
        let d = axis.normalize().expect("cylinder axis must be nonzero");
        let u = d
            .any_perpendicular()
            .expect("a unit axis has a perpendicular");
        let c0 = base;
        let c1 = base + d * height;
        let p_bot = c0 + u * radius; // seam point on the bottom rim (angle 0)
        let p_top = c1 + u * radius; // seam point on the top rim

        let v_bot = self.vertices.push(Vertex {
            point: p_bot,
            origin: Origin::Constructed,
            // A cylinder's seam vertex sits on a curved surface, which `ThreePlane`
            // cannot name (M6). Its coordinate stays the only description it has.
            definition: None,
        });
        let v_top = self.vertices.push(Vertex {
            point: p_top,
            origin: Origin::Constructed,
            // A cylinder's seam vertex sits on a curved surface, which `ThreePlane`
            // cannot name (M6). Its coordinate stays the only description it has.
            definition: None,
        });

        // Rims are full circles seamed at their vertex (start == end); the seam is
        // a straight edge joining the two rim seam points.
        let bottom = {
            let curve = self.curves.push(Curve::Circle(
                Circle::from_center_normal(c0, d, u, radius).expect("non-degenerate rim"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_bot, v_bot]),
                origin: Origin::Constructed,
            })
        };
        let top = {
            let curve = self.curves.push(Curve::Circle(
                Circle::from_center_normal(c1, d, u, radius).expect("non-degenerate rim"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_top, v_top]),
                origin: Origin::Constructed,
            })
        };
        let seam = {
            let curve = self.curves.push(Curve::Line(
                Line::through_points(p_bot, p_top).expect("positive height"),
            ));
            self.edges.push(Edge {
                curve,
                bounds: Some([v_bot, v_top]),
                origin: Origin::Constructed,
            })
        };

        // Lateral cylindrical face: one loop wrapping the seam twice (opposite).
        let lateral = {
            let surface = self.push_surface(
                Surface::Cylinder(
                    Cylinder::from_axis(c0, d, u, radius).expect("non-degenerate cylinder"),
                ),
                SurfaceDef::Constructed,
            );
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: bottom,
                        forward: true,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: true,
                    },
                    HalfEdge {
                        edge: top,
                        forward: false,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: false,
                    },
                ],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Bottom cap: outward normal −d, the bottom rim reversed.
        let bottom_cap = {
            let surface = self.push_surface(
                Surface::Plane(Plane::from_point_normal(c0, -d).expect("nonzero axis")),
                SurfaceDef::Constructed,
            );
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: bottom,
                    forward: false,
                }],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Top cap: outward normal +d, the top rim forward.
        let top_cap = {
            let surface = self.push_surface(
                Surface::Plane(Plane::from_point_normal(c1, d).expect("nonzero axis")),
                SurfaceDef::Constructed,
            );
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: top,
                    forward: true,
                }],
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };

        let shell = self.shells.push(Shell {
            faces: vec![lateral, bottom_cap, top_cap],
        });
        self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;
    use proptest::prelude::*;

    /// ★★★ **The defect the rational coefficients exist to remove, pinned from both sides.**
    ///
    /// Two boxes meet on the plane `x = 3` with faces of different size. `Plane` stores an
    /// un-normalized normal whose length follows that size, so `d = −raw·origin` is a differently
    /// rounded product on each side and the two coefficient vectors are **not exactly
    /// proportional** — the f64 test says "different planes" about one plane. Measured across the
    /// census, 18 pairs are merged only because a second test looks at the faces' coordinates
    /// instead (`docs/dev-log.md`).
    ///
    /// The rational coefficients are built from the corners the caller wrote and canonicalized, so
    /// they have no scale to disagree about and come out **equal**.
    ///
    /// Both halves are load-bearing. If the first assertion ever fails the f64 defect was fixed
    /// somewhere else and this test should be re-read, not deleted; if the second fails the
    /// rational path stopped reaching these surfaces.
    #[test]
    fn two_faces_of_one_plane_disagree_in_f64_and_agree_in_the_rationals() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 2.2, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([3.0, 0.0, 0.0]),
            Point3::from_array([5.0, 13.2, 1.0]),
        );
        // Each solid's face on x = 3: a's outward +X, b's outward −X.
        let face_on_x3 = |s: Handle<Solid>, want_x: f64| -> Handle<Surface> {
            let shell = m.solids.get(s).outer;
            *m.shells
                .get(shell)
                .faces
                .iter()
                .map(|&fh| &m.faces.get(fh).surface)
                .find(|&&sh| match m.surfaces.get(sh) {
                    Surface::Plane(p) => {
                        let [a, b, c, d] = p.coefficients();
                        b == 0.0 && c == 0.0 && a != 0.0 && (-d / a - want_x).abs() < 1e-12
                    }
                    Surface::Cylinder(_) => false,
                })
                .expect("a face on x = 3")
        };
        let (sa, sb) = (face_on_x3(a, 3.0), face_on_x3(b, 3.0));

        // ★ **They are one handle now** — that is what the rational coefficients bought.
        assert_eq!(sa, sb, "one plane, one surface");
        assert!(m.surface_coeffs.contains_key(&sa), "and it is recorded");

        // The f64 defect that made this necessary, shown on the planes themselves rather than
        // through the model, since the model no longer holds two of them. `Plane` keeps an
        // un-normalized normal whose length follows the face's size, so `d = −raw·origin` is a
        // differently rounded product on each side and the two vectors are not exactly
        // proportional — the f64 test says "different planes" about one plane.
        let wall = |dy: f64| {
            Plane::through_points(
                Point3::from_array([3.0, 0.0, 0.0]),
                Point3::from_array([3.0, dy, 0.0]),
                Point3::from_array([3.0, 0.0, 1.0]),
            )
            .expect("non-degenerate")
        };
        let (pa, pb) = (wall(2.2), wall(13.2));
        assert!(
            !nacre_geom::intersect::planes_coplanar(&pa, &pb),
            "f64 coefficients of one plane at two face sizes: {:?} vs {:?}",
            pa.coefficients(),
            pb.coefficients()
        );
        let rat = |dy: f64| {
            let r = |x: f64| nacre_scalar::Rat::from_decimal(x).expect("decimal");
            nacre_scalar::plane_through_points(
                [r(3.0), r(0.0), r(0.0)],
                [r(3.0), r(dy), r(0.0)],
                [r(3.0), r(0.0), r(1.0)],
            )
        };
        assert_eq!(
            rat(2.2),
            rat(13.2),
            "the rationals have no scale to disagree about"
        );
        assert!(rat(2.2).is_some());
    }

    fn build(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
        m.rebuild_adjacency();
        m
    }

    fn he_start(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).bounds.unwrap();
        if he.forward { a } else { b }
    }
    fn he_end(m: &Model, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = m.edges.get(he.edge).bounds.unwrap();
        if he.forward { b } else { a }
    }
    fn face_plane_normal(m: &Model, f: &Face) -> Vector3 {
        match m.surfaces.get(f.surface) {
            Surface::Plane(p) => p.normal(),
            // Planar-only helper: callers filter to plane faces (caps), never cylinders.
            Surface::Cylinder(_) => unreachable!("face_plane_normal called on a curved face"),
        }
    }
    fn face_centroid(m: &Model, f: &Face) -> Point3 {
        let pts: Vec<Point3> = f
            .outer
            .half_edges
            .iter()
            .map(|he| m.vertices.get(he_start(m, *he)).point)
            .collect();
        Point3::centroid(&pts).unwrap()
    }

    // --- golden ---

    #[test]
    fn cuboid_counts() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.edges.len(), 12);
        assert_eq!(m.faces.len(), 6);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        assert_eq!(m.surfaces.len(), 6);
        assert_eq!(m.curves.len(), 12);
    }

    #[test]
    fn cuboid_corner_points() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        let expected = vec![
            [-2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 4.0, 0.0],
            [-2.0, 4.0, 0.0],
            [-2.0, 1.0, 10.0],
            [3.0, 1.0, 10.0],
            [3.0, 4.0, 10.0],
            [-2.0, 4.0, 10.0],
        ];
        let got: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn face_normals_point_outward() {
        let m = build([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let center = Point3::origin().lerp(Point3::from_array([2.0, 3.0, 4.0]), 0.5);
        for (_, f) in m.faces.iter() {
            let outward = face_plane_normal(&m, f).dot(face_centroid(&m, f) - center);
            assert!(outward > 0.0);
        }
    }

    #[test]
    fn every_edge_used_twice_opposite() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1);
        }
    }

    #[test]
    fn every_vertex_incident_to_three_edges() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(m.adj.vertex_edges.len(), 8);
        for edges in m.adj.vertex_edges.values() {
            assert_eq!(edges.len(), 3);
        }
    }

    #[test]
    fn outer_loops_are_closed() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        for (_, f) in m.faces.iter() {
            let hes = &f.outer.half_edges;
            assert_eq!(hes.len(), 4);
            for i in 0..hes.len() {
                assert_eq!(he_end(&m, hes[i]), he_start(&m, hes[(i + 1) % hes.len()]));
            }
        }
    }

    #[test]
    fn edge_endpoints_lie_on_their_curve() {
        let m = build([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        for (_, e) in m.edges.iter() {
            let curve = m.curves.get(e.curve);
            let [a, b] = e.bounds.unwrap();
            assert!(curve.contains(m.vertices.get(a).point, 1e-9));
            assert!(curve.contains(m.vertices.get(b).point, 1e-9));
        }
    }

    #[test]
    fn euler_poincare_holds() {
        let m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (v, e, f) = (
            m.vertices.len() as i64,
            m.edges.len() as i64,
            m.faces.len() as i64,
        );
        // V − E + F = 2(S − G) + L_i, with S=1, G=0, L_i=0. Formal validate:
        // nacre-validate (next unit).
        assert_eq!(v - e + f, 2);
    }

    // --- cylinder (seam b-rep) --- (validate-clean lives in nacre-validate)

    fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
        let mut m = Model::new();
        m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
        m.rebuild_adjacency();
        m
    }

    #[test]
    fn cylinder_counts_and_euler() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        assert_eq!(m.vertices.len(), 2);
        assert_eq!(m.edges.len(), 3);
        assert_eq!(m.faces.len(), 3);
        assert_eq!(m.shells.len(), 1);
        assert_eq!(m.solids.len(), 1);
        // Euler χ = V − E + F = 2 (one shell, genus 0, no inner loops).
        assert_eq!(
            m.vertices.len() as i64 - m.edges.len() as i64 + m.faces.len() as i64,
            2
        );
    }

    #[test]
    fn cylinder_geometry() {
        // +Z axis, r=2, h=5. Seam direction is X (least-aligned axis of +Z).
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let pts: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        // Seam direction for +Z is any_perpendicular([0,0,1]) = X×Z = [0,-1,0], so
        // the seam vertices sit at radius 2 along −Y, at z=0 and z=5.
        assert_eq!(pts, vec![[0.0, -2.0, 0.0], [0.0, -2.0, 5.0]]);
        // Two rim circles carry a Circle; the straight seam a Line.
        let mut circles = 0;
        for (_, e) in m.edges.iter() {
            if let Curve::Circle(c) = m.curves.get(e.curve) {
                assert_eq!(c.radius(), 2.0);
                circles += 1;
            }
        }
        assert_eq!(circles, 2);
    }

    #[test]
    fn cylinder_caps_point_outward() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        // Planar caps only (the lateral cylindrical face has no single normal).
        let mut caps = 0;
        for (_, f) in m.faces.iter() {
            if matches!(m.surfaces.get(f.surface), Surface::Plane(_)) {
                let n = face_plane_normal(&m, f).as_array();
                // Bottom cap → −Z, top cap → +Z (outward along the axis).
                assert!(n == [0.0, 0.0, -1.0] || n == [0.0, 0.0, 1.0]);
                caps += 1;
            }
        }
        assert_eq!(caps, 2);
    }

    #[test]
    fn cylinder_seam_edge_is_self_adjacent() {
        // The novel topology: the seam edge is used twice by the SAME lateral
        // face with opposite orientation (a valid non-manifold-looking but
        // manifold seam). Every edge is still used exactly twice, opposite.
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let mut seam_uses = None;
        for (eh, e) in m.edges.iter() {
            let uses = &m.adj.edge_uses[&eh];
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1); // opposite orientation
            if matches!(m.curves.get(e.curve), Curve::Line(_)) {
                seam_uses = Some(uses.clone());
            }
        }
        let uses = seam_uses.expect("a seam line edge exists");
        assert_eq!(uses[0].0, uses[1].0); // both uses are the same (lateral) face
    }

    // --- proptest ---

    fn box_strategy() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
        (
            prop::array::uniform3(-1e3f64..1e3),
            prop::array::uniform3(1e-2f64..1e3),
        )
            .prop_map(|(min, ext)| {
                let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
                (min, max)
            })
    }

    proptest! {
        #[test]
        fn prop_cuboid_structural((min, max) in box_strategy()) {
            let m = build(min, max);
            prop_assert_eq!(m.vertices.len(), 8);
            prop_assert_eq!(m.edges.len(), 12);
            prop_assert_eq!(m.faces.len(), 6);
            prop_assert_eq!(m.adj.edge_uses.len(), 12);
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
            prop_assert_eq!(m.adj.vertex_edges.len(), 8);
            for edges in m.adj.vertex_edges.values() {
                prop_assert_eq!(edges.len(), 3);
            }
        }

        #[test]
        fn prop_cuboid_outward_normals((min, max) in box_strategy()) {
            let m = build(min, max);
            let center = Point3::from_array(min).lerp(Point3::from_array(max), 0.5);
            for (_, f) in m.faces.iter() {
                prop_assert!(face_plane_normal(&m, f).dot(face_centroid(&m, f) - center) > 0.0);
            }
        }

        #[test]
        fn prop_cuboid_endpoints_on_curves((min, max) in box_strategy()) {
            let m = build(min, max);
            let scale = 1e-6 * (max.iter().map(|x| x.abs()).fold(0.0, f64::max) + 1.0);
            for (_, e) in m.edges.iter() {
                let curve = m.curves.get(e.curve);
                let [a, b] = e.bounds.unwrap();
                prop_assert!(curve.contains(m.vertices.get(a).point, scale));
                prop_assert!(curve.contains(m.vertices.get(b).point, scale));
            }
        }

        #[test]
        fn prop_cylinder_structural(
            base in prop::array::uniform3(-1e3f64..1e3),
            axis in prop::array::uniform3(-1.0f64..1.0),
            r in 0.5f64..10.0,
            h in 0.1f64..10.0,
        ) {
            let axis = Vector3::from_array(axis);
            prop_assume!(axis.norm() > 0.1); // skip near-zero axes
            let mut m = Model::new();
            m.add_cylinder(Point3::from_array(base), axis, r, h);
            m.rebuild_adjacency();
            prop_assert_eq!(m.vertices.len(), 2);
            prop_assert_eq!(m.edges.len(), 3);
            prop_assert_eq!(m.faces.len(), 3);
            // Every edge used exactly twice, opposite orientation (seam included).
            for uses in m.adj.edge_uses.values() {
                prop_assert_eq!(uses.len(), 2);
                prop_assert_ne!(uses[0].1, uses[1].1);
            }
        }
    }

    // --- reversed_shell (M5 containment building block) ---

    #[test]
    fn orientation_flip_is_involution() {
        assert_eq!(Orientation::Forward.flipped(), Orientation::Reversed);
        assert_eq!(Orientation::Reversed.flipped(), Orientation::Forward);
        assert_eq!(
            Orientation::Forward.flipped().flipped(),
            Orientation::Forward
        );
    }

    #[test]
    fn reversed_shell_toggles_orientation_and_reverses_loops() {
        let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let outer = m.solids.get(m.live_solids[0]).outer;
        let f0 = m.shells.get(outer).faces[0];
        let orig = m.faces.get(f0).clone();

        let rev_shell = m.reversed_shell(outer);
        // Fresh cells (not reused faces), fresh shell.
        assert_ne!(rev_shell, outer);
        let rf0 = m.shells.get(rev_shell).faces[0];
        assert_ne!(rf0, f0);
        let rev = m.faces.get(rf0);

        assert_eq!(rev.surface, orig.surface); // surface reused
        assert_eq!(rev.orientation, orig.orientation.flipped());
        let n = orig.outer.half_edges.len();
        assert_eq!(rev.outer.half_edges.len(), n);
        // Reversed winding: he[i] mirrors orig[n-1-i] with the edge reused and
        // the traversal direction flipped.
        for i in 0..n {
            let o = orig.outer.half_edges[n - 1 - i];
            let r = rev.outer.half_edges[i];
            assert_eq!(r.edge, o.edge);
            assert_ne!(r.forward, o.forward);
        }
    }

    #[test]
    fn reversed_shell_is_a_valid_manifold() {
        // Reversing every face's winding preserves the b-rep manifold: each edge
        // is still used by exactly two faces with opposed half-edges. The
        // reversed shell reuses the cube's edges, but the source solid is
        // superseded, so `Adjacency` (reachable-scoped) counts only the reversed
        // faces — no 4-use false positive.
        let mut m = build([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let outer = m.solids.get(m.live_solids[0]).outer;
        let rev = m.reversed_shell(outer);
        let s = m.push_solid(Solid {
            outer: rev,
            cavities: vec![],
        });
        m.live_solids.retain(|&h| h == s); // supersede the original cube
        m.rebuild_adjacency();
        assert_eq!(m.adj.edge_uses.len(), 12);
        for uses in m.adj.edge_uses.values() {
            assert_eq!(uses.len(), 2);
            assert_ne!(uses[0].1, uses[1].1); // opposite forward
        }
    }
}
