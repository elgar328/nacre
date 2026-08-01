//! Feature operations (design §7): the public sketch/extrude/pad/pocket API and the `apply`/
//! `replay` driver. The top layer — it composes the boolean engine ([`crate::boolean`]) and rigid
//! transform ([`crate::transform`]) over the plane substrate below.

use crate::boolean::boolean;
use crate::exact::Swept;
use crate::planes::outer_tri;
use crate::transform::transform;
use crate::{BoolError, he_start};
use nacre_geom::intersect::{
    RingSide, plane_side, point_in_ring_2d, ring_self_intersection, rings_cross,
};
use nacre_geom::{Curve, Line, Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::{Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Orientation, Origin, Shell, Solid, SurfaceDef, Vertex,
    VertexDef,
};

/// A sketch-plane frame: a 2-D point `(u, v)` maps to `origin + u·x + v·y`.
/// `x_axis`/`y_axis` are assumed unit and orthogonal (the constructors ensure
/// it); the plane normal is `x × y`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchPlane {
    pub origin: Point3,
    pub x_axis: Vector3,
    pub y_axis: Vector3,
}

/// A closed planar region: one outer ring and any number of hole rings (straight segments only,
/// at least 3 points each).
///
/// **Winding is not the caller's business, and not this type's either.** The rings are stored as
/// given; the prism builder is the single place that decides orientation, because it is the only
/// one that knows the sweep direction (a pocket sweeps *into* a face, which flips what
/// "counter-clockwise" means). Fixing the winding here as well would put that decision in two
/// places, which is exactly how the holes and the outer ring come to disagree.
///
/// The constructors cost nothing and check nothing; [`Profile2d::check`] states the contract and
/// **every operation that consumes a profile runs it first**, so the kernel never works from an
/// unverified one. `sketch::from_rings` is the checking constructor for loose rings.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile2d {
    outer: Vec<Point2>,
    inners: Vec<Vec<Point2>>,
}

/// Which ring of a [`Profile2d`] an error is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileRing {
    /// The outer ring.
    Outer,
    /// The hole at this index in [`Profile2d::inners`].
    Hole(usize),
}

impl Profile2d {
    /// A simple polygon — no holes.
    pub fn polygon(points: Vec<Point2>) -> Profile2d {
        Profile2d {
            outer: points,
            inners: Vec::new(),
        }
    }

    /// An outer ring with holes. Nothing is verified here; [`Profile2d::check`] is where the
    /// contract is enforced, and every operation that consumes a profile runs it.
    pub fn with_holes(outer: Vec<Point2>, inners: Vec<Vec<Point2>>) -> Profile2d {
        Profile2d { outer, inners }
    }

    /// Verify the contract: every ring is a **simple polygon** of at least three points, the rings
    /// are pairwise disjoint, and each hole lies inside the outer ring and inside no other hole.
    ///
    /// Every operation that consumes a profile calls this first, so an unverified profile never
    /// reaches the topology. It is public because an app can ask before it builds.
    ///
    /// **Simplicity is a contract on authored input, not an invariant of kernel data.** A drawn
    /// ring's inside is *defined* by even-odd, which needs simplicity to mean anything; the
    /// contours a boolean produces get their inside from the arrangement instead, and those may
    /// legitimately be non-simple (a figure-8 pinch — see `loop_winding`). Do not "unify" the two.
    ///
    /// Exact: every decision is an `orient2d` sign on the coordinates as written. The 2-D ring is
    /// then placed into 3-D in `f64`, which is a separate exactness question (the same one every
    /// constructed coordinate has), so this rejects what the *author* drew, not every ring that
    /// could conceivably self-intersect after placement.
    ///
    /// `O(n²)` in the ring size, and it runs on every extrude and every replay of one. Measured
    /// on a convex ring (the worst case — nothing short-circuits): 47 µs at 100 points, 4.6 ms at
    /// 1 000, 115 ms at 5 000. Hand-written sketches are nowhere near that; a generator or an
    /// import that emits thousands of points per ring would feel it, and that is when a
    /// sweep-line is worth its robustness cost — not before.
    pub fn check(&self) -> Result<(), OpError> {
        let rings = || {
            std::iter::once((ProfileRing::Outer, &self.outer)).chain(
                self.inners
                    .iter()
                    .enumerate()
                    .map(|(i, r)| (ProfileRing::Hole(i), r)),
            )
        };
        for (id, r) in rings() {
            if r.len() < 3 {
                return Err(OpError::DegenerateProfile);
            }
            // Simplicity first: `point_in_ring_2d` below is only meaningful on a simple ring.
            // The predicate reports a zero-length edge as a pair with itself.
            if let Some((a, b)) = ring_self_intersection(r) {
                return Err(if a == b {
                    OpError::ZeroLengthProfileEdge { ring: id, edge: a }
                } else {
                    OpError::SelfIntersectingProfile {
                        ring: id,
                        edges: (a, b),
                    }
                });
            }
        }
        let all: Vec<(ProfileRing, &Vec<Point2>)> = rings().collect();
        for (i, (ida, a)) in all.iter().enumerate() {
            for (idb, b) in &all[i + 1..] {
                if rings_cross(a, b) {
                    return Err(OpError::ProfileRingsMeet { a: *ida, b: *idb });
                }
            }
        }
        // The rings are disjoint, so any one vertex answers for a whole ring.
        for (h, hole) in self.inners.iter().enumerate() {
            if point_in_ring_2d(hole[0], &self.outer) != RingSide::Inside {
                return Err(OpError::HoleNotInsideOuter { hole: h });
            }
            for (k, other) in self.inners.iter().enumerate() {
                if k != h && point_in_ring_2d(hole[0], other) == RingSide::Inside {
                    // A ring inside a hole is an island — material again, so it belongs to a
                    // profile of its own. `sketch::from_rings` is what splits those out.
                    return Err(OpError::NestedHole { outer: k, inner: h });
                }
            }
        }
        Ok(())
    }

    /// The outer ring, as given.
    pub fn outer(&self) -> &[Point2] {
        &self.outer
    }

    /// The hole rings, as given.
    pub fn inners(&self) -> &[Vec<Point2>] {
        &self.inners
    }
}

/// A modelling operation.
#[derive(Clone, Debug, PartialEq)]
pub enum Operation {
    /// Extrude `profile` (on `plane`) by `dist` along the plane normal.
    Extrude {
        plane: SketchPlane,
        profile: Profile2d,
        dist: f64,
    },
    /// Pad a boss: extrude `profile` on a planar `face` into a tool prism (height
    /// `dist`) and `Fuse` it onto the solid — boolean sugar over [`Operation::Boolean`],
    /// not a direct face-split. No "profile inside the face" constraint: an overhanging
    /// footprint is handled by the boolean's coplanar-contact / overhang path. Adds
    /// material (design §6).
    PadOnFace {
        face: Handle<Face>,
        profile: Profile2d,
        dist: f64,
    },
    /// Carve a blind pocket: extrude `profile` on a planar `face` into a tool prism
    /// (depth `dist`) and `Cut` it from the solid — boolean sugar over
    /// [`Operation::Boolean`], not a direct face-split. No "profile inside the face"
    /// constraint (overhang footprints route through the boolean). A cut that would
    /// punch through is rejected as not-blind. Removes material (design §6).
    PocketOnFace {
        face: Handle<Face>,
        profile: Profile2d,
        dist: f64,
    },
    /// Boolean of two live solids (design §8 M5). M5-c3 implements only
    /// `Common` (intersection) of convex planar solids; other kinds/inputs are
    /// rejected with [`BoolError`].
    Boolean {
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    },
    /// Rigid-body transform: supersede `solid` by its image under `isometry`
    /// (overhaul stage 1). 1a realizes a rational translation; 1b adds rotation.
    /// The `Isometry` is the exact definition (op-log truth); the geometry is a
    /// realized cache.
    Transform {
        solid: Handle<Solid>,
        isometry: Isometry,
    },
    /// Reflect `solid` in the coordinate plane `axis = offset`, superseding it (pair with
    /// [`Operation::Copy`] to keep the original — the usual move, since mirroring exists to build
    /// the other half of a symmetric part). Lengths are preserved and handedness is reversed:
    /// this is a reflection, not a negative scale.
    Mirror {
        solid: Handle<Solid>,
        axis: Axis,
        offset: Rat,
    },
    /// Duplicate `solid` in place, keeping the original live — **the only operation that adds to
    /// `live_solids` without removing anything** (design §2). Every other edit supersedes its
    /// input, so this is what makes "cut with the same tool twice", "keep the original and a moved
    /// copy", and pattern/mirror sugar expressible at all.
    Copy { solid: Handle<Solid> },
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

/// A failure while applying an operation.
#[derive(Debug, PartialEq)]
pub enum OpError {
    /// A profile ring with fewer than 3 points.
    DegenerateProfile,
    /// A profile ring meets itself — it crosses, touches, doubles back over one of its own edges,
    /// or doubles back over one of its own edges. Such a ring has no unambiguous inside, so the
    /// region it is supposed to bound is undefined. `edges` are the two offending edge indices
    /// within that ring. The input is wrong; this is not a missing capability.
    SelfIntersectingProfile {
        /// Which ring of the profile.
        ring: ProfileRing,
        /// The offending pair of edge indices within that ring.
        edges: (usize, usize),
    },
    /// A profile ring repeats a point, so one of its edges has zero length. It would leave a
    /// degenerate edge in the topology. (Reported apart from
    /// [`OpError::SelfIntersectingProfile`] because "you typed the same point twice" and "your
    /// outline crosses itself" are different mistakes to the author, though the same predicate
    /// finds them.)
    ZeroLengthProfileEdge {
        /// Which ring of the profile.
        ring: ProfileRing,
        /// The zero-length edge's index within that ring.
        edge: usize,
    },
    /// Two rings of a profile touch or cross. A hole must lie strictly inside the outer ring and
    /// strictly outside its siblings; anything else has no unambiguous inside.
    ProfileRingsMeet {
        /// The first ring.
        a: ProfileRing,
        /// The second ring.
        b: ProfileRing,
    },
    /// A hole ring lies outside the outer ring, so it cuts nothing. Left unchecked this is a
    /// *silent* wrong: the prism builds, `validate` is clean, and the volume comes out reduced by
    /// a hole that is not there.
    HoleNotInsideOuter {
        /// The index in [`Profile2d::inners`].
        hole: usize,
    },
    /// A hole ring lies inside another hole. That region is material again — an island — and
    /// belongs to a profile of its own; `sketch::from_rings` splits those out by nesting depth.
    NestedHole {
        /// The containing hole's index.
        outer: usize,
        /// The contained hole's index.
        inner: usize,
    },
    /// A non-positive extrusion distance.
    NonPositiveDistance,
    /// A curve/surface construction collapsed (collinear/coincident points, a
    /// zero-length profile edge).
    DegenerateGeometry,
    /// A pad/pocket target face is not planar (only planar faces carry a sketch frame;
    /// curved-face features arrive with the quadric milestones).
    NonPlanarFace,
    /// A pad/pocket target face belongs to no live solid's outer shell (a stale
    /// or non-live handle).
    FaceNotInLiveSolid,
    /// A pocket's depth reaches through the solid: the carved prism is not blind, so `Cut`
    /// produced a through-hole with no floor face. `pocket` requires `dist` less than the
    /// thickness at the face (the boolean pocket path honestly rejects instead of the old
    /// direct path's silent invalid result).
    PocketNotBlind,
    /// A pad's footprint does not meet the face at all: the `Fuse` came back severed, which two
    /// one-shell solids can only do if they never touched. Like [`OpError::PocketNotBlind`] this is
    /// the *operation's* premise breaking, not a boolean failure — the boolean answered correctly
    /// (a base and a detached boss). Use `Operation::Boolean` directly if two disjoint solids are
    /// what you want. An overhanging footprint still touches and is not this error.
    PadMissesFace,
    /// A boolean operation failed (design §8 M5).
    Boolean(BoolError),
    /// A `Transform` input solid is not live (a stale or non-live handle).
    SolidNotLive,
    /// A `Discovered` vertex or edge of a `Transform`/`Copy` input names a surface that is not one
    /// of that solid's face surfaces, so its exact definition cannot be carried onto the duplicate.
    /// The solid's provenance is inconsistent — declined rather than aborting the kernel, and
    /// unreached in the suite (the definitions and the result faces both name the plane class's
    /// representative surface).
    OriginNotOnSolid,
    /// A `Mirror` input carries curved geometry (a cylindrical face, a circular edge). A
    /// reflection reverses a circle's parametrisation, and which convention a mirrored quadric
    /// should take is a curved-geometry decision, so it is declined rather than guessed.
    MirrorNotPlanar,
}

/// The handles an operation produced. Not `Copy`: `Extrude` carries a `Vec`.
#[derive(Clone, Debug, PartialEq)]
pub enum OpOutput {
    /// The created solid and its faces in push order: `faces[0]` base cap,
    /// `faces[1]` top cap, then one side face per profile edge.
    Extrude {
        solid: Handle<Solid>,
        faces: Vec<Handle<Face>>,
    },
    /// The superseding solid and the boss's top cap face.
    PadOnFace {
        solid: Handle<Solid>,
        top_face: Handle<Face>,
    },
    /// The superseding solid and the pocket's floor face.
    PocketOnFace {
        solid: Handle<Solid>,
        bottom_face: Handle<Face>,
    },
    /// The boolean result solids (supersede both inputs). Usually one; a boolean that severs the
    /// body yields several (cell 0.4), and `Cut(A, A)` (deferred) would yield none.
    Boolean { solids: Vec<Handle<Solid>> },
    /// The transformed solid (supersedes the input).
    Transform { solid: Handle<Solid> },
    /// The duplicate. Unlike every other output, the input stays live alongside it.
    Copy { solid: Handle<Solid> },
    /// The reflected solid (supersedes the input).
    Mirror { solid: Handle<Solid> },
}

/// Apply one operation to `model`, returning the handles it created. Does not
/// rebuild the adjacency cache (do that once after a batch — see [`replay`]).
pub fn apply(model: &mut Model, op: &Operation) -> Result<OpOutput, OpError> {
    match op {
        Operation::Extrude {
            plane,
            profile,
            dist,
        } => {
            let (solid, faces) = extrude(model, plane, profile, *dist)?;
            Ok(OpOutput::Extrude { solid, faces })
        }
        Operation::PadOnFace {
            face,
            profile,
            dist,
        } => {
            let (solid, top_face) = pad(model, *face, profile, *dist)?;
            Ok(OpOutput::PadOnFace { solid, top_face })
        }
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        } => {
            let (solid, bottom_face) = pocket(model, *face, profile, *dist)?;
            Ok(OpOutput::PocketOnFace { solid, bottom_face })
        }
        Operation::Boolean { kind, a, b } => {
            let solids = boolean(model, *kind, *a, *b).map_err(OpError::Boolean)?;
            Ok(OpOutput::Boolean { solids })
        }
        Operation::Transform { solid, isometry } => {
            let out = transform(model, *solid, isometry)?;
            Ok(OpOutput::Transform { solid: out })
        }
        Operation::Copy { solid } => {
            let out = crate::transform::copy(model, *solid)?;
            Ok(OpOutput::Copy { solid: out })
        }
        Operation::Mirror {
            solid,
            axis,
            offset,
        } => {
            let out = crate::transform::mirror(model, *solid, *axis, *offset)?;
            Ok(OpOutput::Mirror { solid: out })
        }
    }
}

/// Replay an operation log into a fresh model. Deterministic: the same log
/// reproduces the same model down to handle indices.
pub fn replay(ops: &[Operation]) -> Result<Model, OpError> {
    let mut model = Model::new();
    for op in ops {
        apply(&mut model, op)?;
    }
    model.rebuild_adjacency();
    Ok(model)
}

fn push_line_edge(
    model: &mut Model,
    a: Handle<Vertex>,
    ap: Point3,
    b: Handle<Vertex>,
    bp: Point3,
) -> Result<Handle<Edge>, OpError> {
    let curve = model.curves.push(Curve::Line(
        Line::through_points(ap, bp).ok_or(OpError::DegenerateGeometry)?,
    ));
    Ok(model.edges.push(Edge {
        curve,
        bounds: Some([a, b]),
        origin: Origin::Constructed,
    }))
}

/// Build a prism: the profile forms the base and (translated by `normal·dist`)
/// the top; each profile edge grows a side quad. Winding generalizes the M1
/// cuboid — base loop reversed (normal −N, outward), top forward (+N), side
/// `(B_i, B_{i+1}, T_{i+1}, T_i)`; every edge is used twice with opposite flags.
pub(crate) fn extrude(
    model: &mut Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    profile.check()?;
    let (outer, holes) = swept_profile(plane, profile, dist);
    build_prism(model, outer, holes, plane.normal(), None)
}

/// Place a profile on its plane and sweep it — **exactly where the plane admits it**.
///
/// The rational path is not an optimization: it is what makes `extrude(7.7)` and
/// `extrude(1.1)` then `extrude(6.6)` put their caps on the same plane rather than an
/// ulp apart. Where it does not apply — a rotated frame, a dimension outside the
/// decimal window — this falls back to the f64 arithmetic that was here before, which
/// is no worse than it was.
///
/// **Mapping only — no winding decision, and no containment check.** Forcing the outer
/// ring CCW here would be a second opinion on a question `build_prism` already answers
/// from the sweep, and two opinions is how an outer ring and its holes end up wound the
/// same way (the pocket case, where the sweep runs `−n` and flips the outer ring). A
/// profile may also reach past the face boundary; an overhanging footprint routes to
/// the overhang boolean sidecars, which reject honestly what they do not cover.
fn swept_profile(plane: &SketchPlane, profile: &Profile2d, dist: f64) -> (Swept, Vec<Swept>) {
    if let Some(rings) = crate::exact::prism_rings(plane, profile, dist) {
        return rings;
    }
    let sweep = plane.normal() * dist;
    let place =
        |ring: &[Point2]| Swept::along(ring.iter().map(|p| plane.point(*p)).collect(), sweep);
    (
        place(profile.outer()),
        profile.inners().iter().map(|h| place(h)).collect(),
    )
}

/// Sweep a profile's rings along `sweep` into a prism solid: caps, side walls, and — for each
/// hole ring — a wall of its own plus an inner loop on each cap.
///
/// **This is the only place winding is decided**, because it is the only place that knows the
/// sweep. The outer ring is normalized CCW **about `sweep`** (area vector dotted with the sweep
/// normal — frame-independent, unlike `proj2`+`signed_area`, which mis-signs when the sweep runs
/// along a negative dominant axis, e.g. a pocket into a `+z` face sweeping `−z`). Each hole is
/// then normalized to the **opposite** sense *of that normalized outer ring*, never of the input:
/// get that backwards and a pocket — where the sweep flips the outer ring — silently produces
/// holes wound the same way as the outer, which is not a hole at all.
///
/// Opposite winding is all a hole needs. The wall quads are built from the ring's own traversal,
/// so a reversed ring yields walls facing into the hole (out of the material), and the cap loops
/// come out opposed to the cap's outer loop, which is what makes them holes.
///
/// All vertices are `Origin::Constructed`. Returns the solid and its faces: `faces[0]` = base cap
/// (at the ring, normal `−ŝ`), `faces[1]` = far cap, then the outer walls, then each hole's walls.
/// Shared by [`extrude`] (a boss) and the pocket (`sweep = −n`).
pub(crate) fn build_prism(
    model: &mut Model,
    outer_ring: Swept,
    inner_rings: Vec<Swept>,
    normal: Vector3,
    base_cap_surface: Option<Handle<Surface>>,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if outer_ring.base.len() < 3 || inner_rings.iter().any(|h| h.base.len() < 3) {
        return Err(OpError::DegenerateProfile);
    }

    let outer_pts = oriented_ring(outer_ring, normal, true);
    // `false` = opposite to the normalized outer ring, whichever way that ended up.
    let hole_pts: Vec<Swept> = inner_rings
        .into_iter()
        .map(|h| oriented_ring(h, normal, false))
        .collect();

    // ★★★ **Surfaces before topology.** A vertex is defined by the three faces that meet at it,
    // and `Store` is append-only, so the handles have to exist before the vertex does. The two
    // arenas are separate, so interleaving them differently does not shift either one's numbering
    // — but the surfaces' order *among themselves* is what the numbering depends on, and it is
    // preserved exactly: base cap, then top cap, then walls, outer ring before the holes.

    // Base cap: outward normal −N, loops reversed.
    // When padding/pocketing on a face, reuse that face's `Surface` handle (explicit sharing) so
    // the flush contact is a shared-handle coplanar pair the boolean can recognize by `Handle`
    // identity; otherwise push a fresh plane. The materialized outward normal must stay −N, so the
    // face orientation is chosen from the shared surface's stored normal — `surface` and
    // `orientation` travel together, and the reconstruction copies both.
    let (base_surface, base_orient) = match base_cap_surface {
        Some(h) => {
            let n_h = match model.surfaces.get(h) {
                Surface::Plane(p) => p.normal(),
                Surface::Cylinder(_) => return Err(OpError::DegenerateGeometry),
            };
            let orient = if n_h.dot(-normal) > 0.0 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            (h, orient)
        }
        None => {
            // The two caps' rational coefficients come from the same frame normal the f64 pair
            // above uses, so they agree with the walls that meet them.
            let caps = outer_pts.exact.as_ref().and_then(|e| e.cap_planes());
            let (s, flipped) = model.push_surface_with_coeffs(
                Surface::Plane(
                    Plane::from_point_normal(outer_pts.base[0], -normal)
                        .ok_or(OpError::DegenerateGeometry)?,
                ),
                SurfaceDef::Constructed,
                caps.map(|(base, _)| base),
            );
            // The plane was built with `−N` as its normal, so `Forward` is what states an outward
            // `−N` — unless a shared surface points the other way, which `flipped` reports.
            let orient = if flipped {
                Orientation::Forward.flipped()
            } else {
                Orientation::Forward
            };
            (s, orient)
        }
    };
    // Top cap: outward normal +N.
    let (top_surface, top_flipped) = model.push_surface_with_coeffs(
        Surface::Plane(
            Plane::from_point_normal(outer_pts.top[0], normal)
                .ok_or(OpError::DegenerateGeometry)?,
        ),
        SurfaceDef::Constructed,
        outer_pts
            .exact
            .as_ref()
            .and_then(|e| e.cap_planes())
            .map(|(_, top)| top),
    );
    let top_orient = if top_flipped {
        Orientation::Forward.flipped()
    } else {
        Orientation::Forward
    };

    // Wall surfaces, in the same order the faces will be emitted: the outer ring's, then each
    // hole's (facing into the hole).
    let outer_walls = wall_surfaces(model, &outer_pts)?;
    let hole_walls: Vec<Vec<(Handle<Surface>, bool)>> = hole_pts
        .iter()
        .map(|h| wall_surfaces(model, h))
        .collect::<Result<_, _>>()?;

    // ── Topology. Every surface it needs already exists.
    let caps = (base_surface, top_surface);
    let outer = sweep_ring(model, &outer_pts, &outer_walls, caps)?;
    let holes: Vec<RingCells> = hole_pts
        .iter()
        .zip(hole_walls.iter())
        .map(|(h, w)| sweep_ring(model, h, w, caps))
        .collect::<Result<_, _>>()?;

    let mut faces =
        Vec::with_capacity(2 + outer.len() + holes.iter().map(|h| h.len()).sum::<usize>());
    faces.push(model.faces.push(Face {
        surface: base_surface,
        outer: outer.cap_loop(Cap::Base),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Base)).collect(),
        orientation: base_orient,
    }));
    faces.push(model.faces.push(Face {
        surface: top_surface,
        outer: outer.cap_loop(Cap::Top),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Top)).collect(),
        orientation: top_orient,
    }));
    for (ring, walls) in
        std::iter::once((&outer, &outer_walls)).chain(holes.iter().zip(hole_walls.iter()))
    {
        ring.push_walls(model, walls, &mut faces);
    }

    let shell = model.shells.push(Shell {
        faces: faces.clone(),
    });
    let solid = model.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    Ok((solid, faces))
}

/// Which cap a ring's loop is being built for.
#[derive(Clone, Copy)]
enum Cap {
    Base,
    Top,
}

/// One ring swept into cells: the two rings of vertices and the three edge families that join
/// them. Built the same way for the outer ring and for a hole — the difference is only which way
/// the ring runs, which the caller has already decided.
struct RingCells {
    /// Kept for [`RingCells::len`]; the geometry itself is read off the `Swept` these came from,
    /// which is also where the wall planes were built (`wall_surfaces`).
    base_pts: Vec<Point3>,
    be: Vec<Handle<Edge>>, // base  B_i -> B_{i+1}
    te: Vec<Handle<Edge>>, // top   T_i -> T_{i+1}
    ve: Vec<Handle<Edge>>, // riser B_i -> T_i
}

impl RingCells {
    fn len(&self) -> usize {
        self.base_pts.len()
    }

    /// The cap loop for this ring. The base cap faces `−ŝ`, so its loops run backwards.
    fn cap_loop(&self, cap: Cap) -> Loop {
        let n = self.len();
        match cap {
            Cap::Base => Loop {
                half_edges: (0..n)
                    .rev()
                    .map(|i| HalfEdge {
                        edge: self.be[i],
                        forward: false,
                    })
                    .collect(),
            },
            Cap::Top => Loop {
                half_edges: (0..n)
                    .map(|i| HalfEdge {
                        edge: self.te[i],
                        forward: true,
                    })
                    .collect(),
            },
        }
    }

    /// One quad per ring segment. The quad's winding follows the ring's, so a ring wound against
    /// the outer one yields walls whose normals point into the hole.
    fn push_walls(
        &self,
        model: &mut Model,
        walls: &[(Handle<Surface>, bool)],
        faces: &mut Vec<Handle<Face>>,
    ) {
        let n = self.len();
        for (i, &(surface, flipped)) in walls.iter().enumerate().take(n) {
            let j = (i + 1) % n;
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: self.be[i],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.ve[j],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.te[i],
                        forward: false,
                    },
                    HalfEdge {
                        edge: self.ve[i],
                        forward: false,
                    },
                ],
            };
            faces.push(model.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                // The quad's winding is the ring's, which is what the plane above was built from;
                // a shared surface pointing the other way spells the same outward as `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            }));
        }
    }
}

/// One plane per ring segment, pushed **before** any of the ring's topology exists — see
/// `build_prism`. The same three points `push_walls` used to build them from, and in the same
/// order, so the surface arena's numbering is untouched.
fn wall_surfaces(model: &mut Model, ring: &Swept) -> Result<Vec<(Handle<Surface>, bool)>, OpError> {
    let n = ring.base.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            Ok(model.push_surface_with_coeffs(
                Surface::Plane(
                    Plane::through_points(ring.base[i], ring.base[j], ring.top[i])
                        .ok_or(OpError::DegenerateGeometry)?,
                ),
                SurfaceDef::Constructed,
                // Same three points, in rationals — so this wall and any other face of the same
                // plane record one array. `None` here is the f64 path or an i128 overflow.
                ring.exact.as_ref().and_then(|e| e.wall_plane(i)),
            ))
        })
        .collect()
}

/// `pts` wound counter-clockwise about `normal` when `ccw`, clockwise when not. The test is the
/// polygon's area vector against `normal`, so it does not care which axis dominates.
/// **Requires a nonzero area.** The winding is read from the sign of the area vector, and a ring
/// that encloses nothing gives zero — the comparison below would then pick a side by accident.
/// [`Profile2d::check`] is what guarantees it: a simple polygon cannot have zero area, and a ring
/// that folds back on itself (a symmetric bowtie cancels to exactly zero) is not simple.
fn oriented_ring(ring: Swept, normal: Vector3, ccw: bool) -> Swept {
    let v = &ring.base;
    let k = v.len();
    let area_vec = (0..k)
        .map(|i| (v[i] - Point3::origin()).cross(v[(i + 1) % k] - Point3::origin()))
        .fold(Vector3::from_array([0.0; 3]), |a, b| a + b);
    if (area_vec.dot(normal) < 0.0) == ccw {
        ring.reversed()
    } else {
        ring
    }
}

/// Push one ring's vertices and edges (base ring, top ring, risers).
///
/// The top ring arrives already computed rather than being derived here as
/// `base + sweep`: where the frame allows it that arithmetic is done in exact rationals
/// (see [`crate::exact`]), and a dimension split into two then lands on the same points
/// as the undivided one instead of an ulp away.
fn sweep_ring(
    model: &mut Model,
    ring: &Swept,
    walls: &[(Handle<Surface>, bool)],
    caps: (Handle<Surface>, Handle<Surface>),
) -> Result<RingCells, OpError> {
    let n = ring.base.len();
    let base_pts: Vec<Point3> = ring.base.clone();
    let top_pts: Vec<Point3> = ring.top.clone();
    // Corner `i` is where the wall before it, the wall after it, and the cap meet.
    //
    // ★ **Two walls that are one plane name a line, not a point.** A profile with a collinear
    // vertex makes exactly that — and since surfaces are interned, "one plane" *is* "one handle",
    // so the check is a comparison. Such a vertex has no three-plane definition, and saying so is
    // more useful than inventing one.
    let define = |i: usize, cap: Handle<Surface>| -> Option<VertexDef> {
        let prev = walls[(i + n - 1) % n].0;
        let here = walls[i].0;
        (prev != here && prev != cap && here != cap)
            .then_some(VertexDef::ThreePlane([prev, here, cap]))
    };
    let push_verts =
        |model: &mut Model, ps: &[Point3], cap: Handle<Surface>| -> Vec<Handle<Vertex>> {
            ps.iter()
                .enumerate()
                .map(|(i, p)| {
                    model.vertices.push(Vertex {
                        point: *p,
                        origin: Origin::Constructed,
                        definition: define(i, cap),
                    })
                })
                .collect()
        };
    let bv = push_verts(model, &base_pts, caps.0);
    let tv = push_verts(model, &top_pts, caps.1);

    let (mut be, mut te, mut ve) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        let j = (i + 1) % n;
        be.push(push_line_edge(
            model,
            bv[i],
            base_pts[i],
            bv[j],
            base_pts[j],
        )?);
        te.push(push_line_edge(model, tv[i], top_pts[i], tv[j], top_pts[j])?);
        ve.push(push_line_edge(
            model,
            bv[i],
            base_pts[i],
            tv[i],
            top_pts[i],
        )?);
    }
    Ok(RingCells {
        base_pts,
        be,
        te,
        ve,
    })
}

/// A planar face's live solid, its in-plane right-handed frame (`x × y = n`, centred on the face
/// centroid so a profile's `(0,0)` lands there), and its loops — the shared setup for placing a
/// profile on a face (pad / pocket).
struct FaceFrame {
    solid_h: Handle<Solid>,
    surface_h: Handle<Surface>,
    n: Vector3, // outward normal
    x: Vector3,
    y: Vector3,
    origin: Point3, // face centroid
}

/// The sketch plane of a planar face — **the very frame [`Operation::PadOnFace`] and
/// [`Operation::PocketOnFace`] place their profile in**, so a caller can work out where its
/// `(0, 0)` will land before it builds anything.
///
/// That equality is the contract, not a coincidence: this is a projection of the frame those
/// operations use, never a second derivation. A test pins a hand-placed profile against a pad to
/// keep it that way.
///
/// `NonPlanarFace` for a curved surface (only a plane carries a frame); `FaceNotInLiveSolid` if no
/// live solid's outer shell holds the face.
pub fn face_plane(model: &Model, face: Handle<Face>) -> Result<SketchPlane, OpError> {
    let f = face_frame(model, face)?;
    Ok(SketchPlane {
        origin: f.origin,
        x_axis: f.x,
        y_axis: f.y,
    })
}

/// Locate `face`'s live solid and build its planar frame. `NonPlanarFace` for a curved surface,
/// `FaceNotInLiveSolid` if no live outer shell holds it.
fn face_frame(model: &Model, face: Handle<Face>) -> Result<FaceFrame, OpError> {
    let (solid_h, _) = model
        .live_solids
        .iter()
        .map(|&s| (s, model.solids.get(s).outer))
        .find(|&(_, sh)| model.shells.get(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.faces.get(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
    let plane = match model.surfaces.get(surface_h) {
        Surface::Plane(p) => *p,
        Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };
    let sign = match orientation {
        Orientation::Forward => 1.0,
        Orientation::Reversed => -1.0,
    };
    let n = plane.normal() * sign;
    let x = n.any_perpendicular().ok_or(OpError::DegenerateGeometry)?;
    let y = n.cross(x);
    // The origin is the face **region's** area centroid, holes included in the
    // subtraction. It was the mean of the outer loop's corners, which is not the same point on a
    // reflex face and — worse — *moves when a vertex is added along a straight edge*, so the same
    // shape could seat a boss in two different places. An area centroid is a property of the
    // region, so it does not care how the boundary is subdivided.
    let loop_pts = |lp: &Loop| -> Vec<Point3> {
        lp.half_edges
            .iter()
            .map(|he| model.vertices.get(he_start(model, *he)).point)
            .collect()
    };
    let outer_pts = loop_pts(&f.outer);
    let holes: Vec<Vec<Point3>> = f.inner.iter().map(&loop_pts).collect();
    let hole_refs: Vec<&[Point3]> = holes.iter().map(|h| h.as_slice()).collect();
    let (_, origin) = nacre_geom::planar_region_area_centroid(&outer_pts, &hole_refs)
        .ok_or(OpError::DegenerateGeometry)?;
    Ok(FaceFrame {
        solid_h,
        surface_h,
        n,
        x,
        y,
        origin,
    })
}

/// A face-local feature built as **tool body + boolean**: the profile
/// extrudes off `face` into a top-flush prism, then `kind` fuses/cuts it against the face's solid.
/// A **contained** footprint takes the contained-coplanar path (empty seam → all
/// `Origin::Constructed`); one that **reaches past the face** routes to the overhang boolean
/// sidecars (a boss cantilever / an edge slot; Discovered seam vertices). `Fuse` sweeps **outward**
/// (a boss); `Cut` sweeps **inward** (a blind pocket). Returns the result solid and the feature's
/// exposed cap — the boss top or the pocket floor, the outer-shell face on the prism's far-cap plane
/// with outward normal `+n`. `Option::None` there ⇒ the far cap did not survive (a through-cut with
/// no floor); callers map it to their own error. `NonPositiveDistance`/`DegenerateProfile` propagate
/// from the frame; overhang configurations the boolean does not cover surface as `Boolean(_)`.
fn extrude_and_boolean(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
    kind: BoolKind,
) -> Result<(Handle<Solid>, Option<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    profile.check()?;
    let frame = face_frame(model, face)?;
    // No containment check — an overhanging footprint routes to the overhang boolean sidecars.
    let n = frame.n;
    // Cut carves inward, Fuse raises outward; either way the prism's near cap is flush on the face.
    let signed = if matches!(kind, BoolKind::Cut) {
        -dist
    } else {
        dist
    };
    // The face's own frame, so a pad or pocket takes the same exact-rational path an
    // extrude does: its axes are `{0, ±1}` exactly whenever the face is axis-aligned.
    let plane = SketchPlane {
        origin: frame.origin,
        x_axis: frame.x,
        y_axis: frame.y,
    };
    let (outer, holes) = swept_profile(&plane, profile, signed);
    let (prism, prism_faces) = build_prism(
        model,
        outer,
        holes,
        n * signed.signum(),
        Some(frame.surface_h),
    )?;
    let solids = boolean(model, kind, frame.solid_h, prism).map_err(|e| {
        model.live_solids.retain(|&s| s != prism); // drop the transient prism (atomic on failure)
        OpError::Boolean(e)
    })?;
    // A pad's prism must actually meet the face. Two live solids are each one outer shell, so a
    // `Fuse` of them can only come back severed if they never touched — the footprint missed the
    // face entirely. Returning the piece that carries the cap would hand back a floating boss and
    // silently drop the base, so this is `PadMissesFace`: the boolean succeeded and answered
    // correctly (two solids); it is the *pad's* premise that broke. A footprint that merely
    // overhangs still touches, fuses into one solid, and takes the normal path.
    //
    // `assemble_fuse_cut` already retired the inputs, so restoring `live_solids` is what keeps the
    // "no reject-after-commit" contract true from the outside: the model the caller sees is the one
    // it had before. Only `live_solids` is touched — the store stays append-only.
    if matches!(kind, BoolKind::Fuse) && solids.len() > 1 {
        model.live_solids.retain(|s| !solids.contains(s));
        model.live_solids.push(frame.solid_h);
        return Err(OpError::PadMissesFace);
    }
    // Exposed cap = the result face on the prism's far-cap plane (face plane offset by n·signed),
    // its outward normal +n (the opening side for a pocket, the boss top for a boss). A pocket (Cut)
    // that severs leaves several solids — scan them all for the cap and return the piece that
    // carries it, leaving the others live; that is a valid multi-solid model, not a failure.
    // `build_prism` returns the far cap as `faces[1]`; the store is append-only, so it is still
    // readable after the boolean retired the prism, and it names the cap plane exactly.
    let far_cap = prism_faces[1];
    match solids
        .iter()
        .find_map(|&s| find_face_coplanar_with(model, s, far_cap, n).map(|c| (s, c)))
    {
        Some((solid, cap)) => Ok((solid, Some(cap))),
        // Nothing carries the cap. If anything survived at all, hand it back capless and let
        // `pad`/`pocket` decide; if the boolean came back empty the prism removed the whole solid,
        // which is `PocketNotBlind` taken to its limit — not merely floorless, but nothing left.
        // Only `Cut` can empty a result: `Fuse` of two non-empty solids is never empty.
        None => match solids.first() {
            Some(&primary) => Ok((primary, None)),
            None => {
                debug_assert!(
                    matches!(kind, BoolKind::Cut),
                    "a Fuse cannot produce an empty result"
                );
                Err(OpError::PocketNotBlind)
            }
        },
    }
}

/// Pad a boss on a planar `face`: extrude the profile **outward** by `dist` and `Fuse` it onto the
/// solid, adding `profile_area · dist` of material. Returns `(new solid, top cap face)`. A boss
/// always yields its top cap, so the `None` guard is an unreachable internal-invariant defense.
pub(crate) fn pad(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    let (solid, top) = extrude_and_boolean(model, face, profile, dist, BoolKind::Fuse)?;
    Ok((solid, top.ok_or(OpError::DegenerateGeometry)?))
}

/// Carve a blind pocket on a planar `face`: extrude the profile **inward** by `dist` and `Cut` it
/// from the solid, removing `profile_area · dist` of material. Returns `(new solid, floor face)`.
/// `PocketNotBlind` if `dist` reaches through the solid (the far cap is not blind → a through-cut
/// with no floor face).
pub(crate) fn pocket(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    let (solid, floor) = extrude_and_boolean(model, face, profile, dist, BoolKind::Cut)?;
    Ok((solid, floor.ok_or(OpError::PocketNotBlind)?))
}

/// The outer-shell face of `solid` that lies on `reference`'s plane with its outward normal on
/// `want`'s side — how a pad/pocket recovers its own exposed cap (the boss top, the pocket floor)
/// from the boolean result. `None` if there is none (a through-pocket has no floor).
///
/// **`reference` is a real face, not a `(point, normal)` pair, and that is the point.** Naming the
/// plane by coefficients meant comparing `d = −n·origin` computed at *different* points of the same
/// plane: exact only when the dot happens to reproduce bit for bit, which for an axis-aligned frame
/// it does (`n·p` is one coordinate) and for a slanted one it does not. With a face in hand the
/// question is answered the way the kernel answers identity everywhere else:
///
/// 1. **the same `Surface` handle** — integers, not coordinates (overview §2). `assemble_fuse_cut`
///    gives a result face the surface of the operand plane it came from, so the surviving cap
///    normally lands here.
/// 2. **the faces' own coordinates, exactly** — every `outer_tri` point of the candidate lies on
///    `reference`'s tri plane (`plane_side`, an exact `orient3d` on the points the user gave).
///    Needed because `plane_idx` names a *class representative*: if the cap plane merged with a
///    coplanar face of the other operand, the survivor can carry that operand's surface instead.
///    A `None` from `outer_tri` (no non-collinear triple) means no evidence — that face is skipped,
///    and a degenerate `reference` leaves only branch 1.
///
/// The direction filter reads the **candidate's** outward normal against `want`, never
/// `reference`'s: a pocket's tool cap faces along the sweep (`−n`) while the floor it becomes faces
/// back into the void (`+n`). Coplanarity is settled by then, so the two are parallel and the dot
/// is a full magnitude away from zero — an f64 read whose sign cannot round the wrong way.
///
/// If the cap survives as several faces they all satisfy this, and the first is returned; the
/// coefficient test had the same ambiguity.
///
/// **Measured (2026-07-22): the corpus does not separate the two branches** — disabling either one
/// leaves the whole suite at 207 passed / 23 failed. So branch 2 has no firing test today and is a
/// documented backstop (cf. `NON_MANIFOLD_EDGE`); branch 1 is kept because handle identity is the
/// strongest answer available and is the path a surviving cap normally takes.
pub(crate) fn find_face_coplanar_with(
    model: &Model,
    solid: Handle<Solid>,
    reference: Handle<Face>,
    want: Vector3,
) -> Option<Handle<Face>> {
    let ref_surf = model.faces.get(reference).surface;
    let ref_tri = outer_tri(model, model.faces.get(reference)).map(|(tri, _)| tri);
    let shell = model.solids.get(solid).outer;
    model.shells.get(shell).faces.iter().copied().find(|&fh| {
        let Some((tri, _)) = outer_tri(model, model.faces.get(fh)) else {
            return false;
        };
        let coplanar = model.faces.get(fh).surface == ref_surf
            || ref_tri.is_some_and(|r| tri.iter().all(|&q| plane_side(r, q) == 0));
        coplanar && (tri[1] - tri[0]).cross(tri[2] - tri[0]).dot(want) > 0.0
    })
}
