//! Rigid-body transform of a solid: copy a solid under an isometry, remapping every
//! surface/curve and carrying each vertex's exact definition forward so a rotated operand stays
//! exactly defined. [`copy`] is the same walk with no motion at all.

use crate::OpError;
use nacre_exact::{Axis, Isometry, Rat};
use nacre_geom::{AxisMirror, Cylinder, Plane};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::PointCache;
use nacre_topo::{
    Edge, Face, HalfEdge, Loop, Model, Motion, MotionNode, Shell, Solid, Surface, Vertex,
};
use std::collections::{HashMap, HashSet};

/// Supersede `solid` by its image under `isometry` (a rotation, then a translation). Clones
/// the solid's cells with moved geometry ([`transform_solid`])
/// and drops the input from `live_solids` — the op-log is the truth.
pub(crate) fn transform(
    model: &mut Model,
    solid: Handle<Solid>,
    isometry: &Isometry,
) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids().contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let out = transform_solid(model, solid, &Xform::Rigid(isometry))?;
    model.supersede_live(&[solid]);
    Ok(out)
}

/// An independent twin of `solid` at the same place — **the one operation that only adds to
/// `live_solids`** (supersede semantics). Every other edit supersedes: `transform` and
/// `boolean` drop their inputs, so without this the kernel can *move* a solid but never *copy* one,
/// and "cut with the same tool twice" or "keep the original and a moved copy" cannot be expressed.
///
/// It is [`transform_solid`] under a zero translation: the cells are duplicated (two live solids
/// must not share cells — a shared edge would read as four face uses and break the manifold check)
/// while the geometry is rebuilt bit-for-bit (a pure translation keeps a plane's exact `raw`).
///
/// **Nothing of the source's point cache is carried.** Like every walk, the copy pushes each
/// vertex as [`PointCache::Unrealized`] holding the moved `f64` figure and re-realizes it from
/// the remapped definition; the figure stands only where that realization does not answer (see
/// pass 3).
///
/// A non-live input is rejected rather than resurrected: reusing a superseded handle is a caller
/// bug, and letting it succeed would hide it.
pub(crate) fn copy(model: &mut Model, solid: Handle<Solid>) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids().contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let zero = Isometry::translation([Rat::from_int(0); 3]);
    // `push_solid` registers the twin as live; the input is *not* retained away — that missing
    // line is the whole difference from `transform`.
    transform_solid(model, solid, &Xform::Rigid(&zero))
}

/// Whether every vertex **definition** in `solid` names surfaces the walk will remap — that is,
/// face surfaces of this solid. The remap cannot proceed otherwise, and a kernel must decline
/// rather than abort ("honest-reject > silent-wrong"; `assemble_fuse_cut` takes the same line at
/// its own naming failure), so the three entry points check first and reject.
///
/// ★ **Every definition is checked** — the definition is
/// the vertex's identity, so every def must remap. The gate is expected to fire **zero** times
/// (every producer writes definitions from its own face surfaces — the coverage gate measures
/// 100%), and the suite staying green is that evidence: a firing here is a new rejection, which
/// is a failing test. The positive control (`a_foreign_definition_is_rejected`) shows the gate
/// actually bites.
///
/// ★★ **That expectation is enforced where definitions are made.** The assembly names a result
/// vertex by the canonical triple of its incident faces' planes — not the arrangement's, which
/// can include a plane the result keeps no face on — and refuses, in release builds too, a result
/// whose definitions still name an absent surface (`VertexNamesAbsentSurface`). A solid that could
/// not move is refused by the boolean that made it, not by this gate two operations later.
pub(crate) fn defs_are_remappable(model: &Model, solid: Handle<Solid>) -> bool {
    foreign_named_vertex(model, solid).is_none()
}

/// The offending vertex behind a [`defs_are_remappable`] refusal — the lowest-handle vertex
/// whose definition names a surface the solid keeps no face on (`min` so the witness is the
/// same whatever order the face walk visits it in), `None` when every definition re-solves.
/// One walk answers both spellings so the two cannot drift.
pub(crate) fn foreign_named_vertex(model: &Model, solid: Handle<Solid>) -> Option<Handle<Vertex>> {
    let src = model.solid(solid);
    let shells: Vec<Handle<Shell>> = std::iter::once(src.outer)
        .chain(src.cavities.iter().copied())
        .collect();
    let mut surfs: HashSet<Handle<Surface>> = HashSet::new();
    for &sh in &shells {
        for &fh in &model.shell(sh).faces {
            surfs.insert(model.face(fh).surface);
        }
    }
    let named = |def: &Vertex| def.carriers().all(|s| surfs.contains(&s));
    let mut worst: Option<Handle<Vertex>> = None;
    for &sh in &shells {
        for &fh in &model.shell(sh).faces {
            let face = model.face(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    let edge = model.edge(he.edge);
                    for vh in edge.vertices.iter() {
                        if !named(model.vertex(*vh)) && worst.is_none_or(|w| vh.index() < w.index())
                        {
                            worst = Some(*vh);
                        }
                    }
                }
            }
        }
    }
    worst
}

/// Supersede `solid` by its reflection in the coordinate plane `axis = offset`.
///
/// Lengths are preserved and **handedness is reversed** — this is a reflection, not a negative
/// scale (the kernel has no scale, and `Isometry` cannot hold an improper motion). The plane is
/// axis-aligned, the same restriction rotation already has: a general plane's unit normal is
/// irrational, so the image could not be reproduced from an exact definition.
///
/// Like `transform` this consumes its input; pair it with [`copy`] to keep the original — which
/// is the usual move, since mirroring exists to build the other half of a symmetric part.
pub(crate) fn mirror(
    model: &mut Model,
    solid: Handle<Solid>,
    axis: Axis,
    offset: Rat,
) -> Result<Handle<Solid>, OpError> {
    // A foreign handle answers "yes" here: `Handle`'s equality is its index, deliberately (it has
    // no `T: Eq` bound to give). The guard is one step further in — the first `get` dies with the
    // cross-store message. Re-anchoring in `replay` restores the premise for a log; a handle
    // passed straight to `apply` from another model is a caller bug and stays one.
    if !model.live_solids().contains(&solid) {
        return Err(OpError::SolidNotLive);
    }
    if !defs_are_remappable(model, solid) {
        return Err(OpError::OriginNotOnSolid);
    }
    let out = transform_solid(model, solid, &Xform::mirror(axis, offset))?;
    model.supersede_live(&[solid]);
    Ok(out)
}

/// The motion nodes this `motion` appends to `parent`, or `None` when it records nothing.
///
/// **One rule for all three motions.** An `Isometry` is "rotate about a pivot, then translate",
/// so it contributes up to two nodes in that order — **order is the definition**, since the two
/// do not commute. A reflection contributes one. Carrying a reflection a second way
/// (conjugating the input's chain, `M ∘ R = (M R M⁻¹) ∘ M`, and mirroring its root vertex) would
/// be two mechanisms for one question and a chain that cannot hold the reflection it has just
/// performed. It is a motion; it goes in the chain.
///
/// A node is *omitted* only when the motion changes nothing the definition needs to say: a zero
/// translation is the identity (`copy` is `transform_solid` under one), a 90°-family rotation of
/// an exact datum keeps it exact, and a translation that lands every coordinate back on an exact
/// `f64` does too (`realizes_exactly`, which asks whether the realized image is the exact one).
///
/// **Both exceptions lapse once the datum already has a history.** "This motion kept the
/// coordinates exact" is a statement about *this step*; it says nothing about the chain, and a
/// chain missing a link describes the datum as it was *before* that link — silently. So once
/// `leaf` is `Some`, every motion is recorded.
///
/// The translation half of that used not to be reachable — an inexact rotation leaves coordinates
/// using the full 53-bit mantissa, so no nonzero translation of such a solid is exact. A
/// reflection in a dyadic plane *is* exactness-preserving and can carry a history, which is what
/// makes the rule load-bearing rather than merely correct.
fn chain_motion(
    model: &mut Model,
    parent: Option<Handle<MotionNode>>,
    motion: &Xform<'_>,
    carry: Carry,
) -> Option<Handle<MotionNode>> {
    // ★ **Record what the statements did not absorb — and everything once a history exists**
    // (the exceptions lapse the moment a chain is there to extend). One decision, made per solid
    // by [`carry_of`], answers for both parts of a rigid motion; this function keeps no table of
    // its own (its old `is_exact`/`leaf` filters were the two halves of the half-recorded chain).
    let all = parent.is_some();
    let mut leaf = parent;
    match motion {
        Xform::Rigid(iso) => {
            if let Some(r) = iso.rotate
                && (all || carry == Carry::None)
            {
                leaf = Some(model.push_motion(
                    Motion::Rotate {
                        axis: r.axis,
                        pivot: r.pivot,
                        angle: r.angle,
                    },
                    leaf,
                ));
            }
            if iso.translate.iter().any(|r| r.numer() != 0) && (all || carry != Carry::Full) {
                leaf = Some(model.push_motion(
                    Motion::Translate {
                        offset: iso.translate,
                    },
                    leaf,
                ));
            }
        }
        Xform::Mirror { axis, offset, .. } => {
            if all || carry == Carry::None {
                leaf = Some(model.push_motion(
                    Motion::Mirror {
                        axis: *axis,
                        offset: *offset,
                    },
                    leaf,
                ));
            }
        }
    }
    leaf
}

/// **What of a motion a solid carries into its statements exactly** — the rest is recorded as
/// motion nodes. The law: `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` — one
/// operation behaves as its two would, so an exact turn is transported into the statements and a
/// translation that rounds is recorded behind it, and the chain never describes a datum as it
/// was *before a part the statements already absorbed*.
///
/// ★★★★★ **Why this exists — the half-recorded chain.** The rule this replaced asked two
/// questions in two places: `Isometry::is_exact` for the turn (recorded only when irrational or
/// when a history already existed) and a data probe for the translation (recorded when the
/// data rounded). An exact turn with a rounding translation recorded the translation **alone**,
/// while the statements stayed pre-motion — the truth then said «unturned, shifted» and the
/// cache «turned, shifted», and `world_cylinder_def` folded the chain into a cylinder standing
/// somewhere else (measured: the offset boss under `rz90 + t(5,−3,2)` — its seam vertex
/// `4.8 + 5` rounds — met its plate as *disjoint* in release, the postcondition catching it only
/// in debug). One decision, per solid, answers both questions now.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Carry {
    /// Every statement is transported; nothing is recorded (an exact motion of a fresh solid).
    Full,
    /// The turn is transported and the translation recorded — a rigid motion whose turn is exact
    /// on this data and whose translation rounds.
    Rotation,
    /// Nothing is transported: the whole motion is recorded and the statements stay pre-motion.
    None,
}

/// **Does the realized image of `p` equal its exact image?** — the one spelling of «exact»: the
/// producer's own `f64` arithmetic ([`Xform::point`]) against the rational transport
/// ([`transport_points`]' `each`), coordinate by coordinate. `false` when the exact image does
/// not exist (an irrational turn, or `i128` overflow) — those are recorded.
///
/// ★ The realized side is the producer's arithmetic itself, not the translation tested on the
/// **pre-turn** coordinates (`p_i + t_i`): that test calls a datum exact that is exact before the
/// turn and rounds after it — `(4.8, 0.5) + (−4, 8)` under `rz90`: `0.8`, `8.5` exact, `12.8`
/// rounds — and its cache comes out one ulp off its truth.
fn realizes_exactly(motion: &Xform<'_>, p: Point3) -> bool {
    let q = p.as_array();
    let Some(exact) = (|| {
        let r = [
            Rat::try_from_f64(q[0])?,
            Rat::try_from_f64(q[1])?,
            Rat::try_from_f64(q[2])?,
        ];
        match motion {
            Xform::Rigid(iso) => iso.point_rat(r),
            Xform::Mirror { axis, offset, .. } => nacre_exact::mirror_point_rat(r, *axis, *offset),
        }
    })() else {
        return false;
    };
    let img = motion.point(p).as_array();
    (0..3).all(|i| Rat::try_from_f64(img[i]) == Some(exact[i]))
}

/// **What of `motion` this solid carries into its statements exactly** — per solid, in one pass.
///
/// Each candidate prefix of the motion — the whole of it, its turn alone (same pivot, no
/// translation), nothing — is asked the same question of every datum a judgment reads (vertex
/// coordinates, plane origins, cylinder axis origins) and of every statement pass 1 will
/// transport (`transport_points`, `transport_cylinder` — the very functions, so «this will not
/// overflow» is structural): does its realized image equal its exact image
/// ([`realizes_exactly`])? The first candidate every datum passes is carried; the rest of the
/// motion is recorded ([`chain_motion`]).
///
/// **Per solid, not per face.** `transform_solid` requires a solid's boundary vertices to share
/// one node ("uniform motion"), and a dyadic parameter does split by magnitude — measured, `0.5`
/// is exact at `p = 1.0` and rounds at `p = 1e17`. Deciding per face would leave one solid with
/// some vertices carrying a node and some not, and break that invariant; so one datum failing a
/// candidate fails it for the whole solid — a `Through` plane (a statement of handles, which
/// no transport can move) puts the whole solid on the recorded path, turn included.
///
/// (The per-plane **surface-statement** exemption — pass 1's `invariant` flag, a plane the
/// motion fixes as a set — is a different question and does not touch this one: the vertices
/// of such a plane still move and still share the solid's node story; only the plane's own
/// statement needed no new spelling.)
fn carry_of(model: &Model, solid: Handle<Solid>, motion: &Xform<'_>) -> Carry {
    let turn_only: Option<Isometry> = match motion {
        Xform::Rigid(iso) => iso.rotate.map(Isometry::rotation),
        Xform::Mirror { .. } => None,
    };
    let turn_xform: Option<Xform<'_>> = turn_only.as_ref().map(Xform::Rigid);
    let candidates: Vec<(Carry, &Xform<'_>)> = std::iter::once((Carry::Full, motion))
        .chain(turn_xform.as_ref().map(|x| (Carry::Rotation, x)))
        .collect();
    let src = model.solid(solid);
    // ★ The surfaces' exact points must survive the no-node path too. An exact motion
    // carries a `Constructed` surface's rational triple through `point_rat`/`mirror_point_rat`,
    // and that arithmetic can overflow `i128` even when every f64 above is exact (the two
    // conditions are independent). Dropping the points would leave a plane with no exact
    // statement; recording a node instead keeps the original triple as the pre-motion truth.
    let points_move = |m: &Xform<'_>, s: Handle<Surface>| -> bool {
        // ★★★★★ **An exhaustive `match`, deliberately — this used to be the most dangerous
        // `let`-`else` in the file, and the compiler could not see it.** Adding
        // `PlanePoints::Through` produced exactly two non-exhaustive-match errors and **not**
        // one here — the fallback arm would have swallowed the new variant silently, answering
        // for geometry it had never seen. Spelled as a match, the next variant (M6's) is a
        // compile error at exactly this decision.
        match model.surface(s) {
            nacre_topo::Surface::Plane {
                points: nacre_topo::PlanePoints::Known(p),
                ..
            } => transport_points(m, *p).is_some(),
            // A `Through` plane's truth carries geometry **by reference** — no transport can
            // move it while the f64 cache moves, so it refuses every carrying candidate.
            nacre_topo::Surface::Plane {
                points: nacre_topo::PlanePoints::Through(_),
                ..
            } => false,
            nacre_topo::Surface::Cylinder { def, .. } => transport_cylinder(m, def).is_some(),
        }
    };
    let mut ok: Vec<bool> = vec![true; candidates.len()];
    let probe = |ok: &mut Vec<bool>, p: Point3| {
        for (k, (_, m)) in candidates.iter().enumerate() {
            if ok[k] && !realizes_exactly(m, p) {
                ok[k] = false;
            }
        }
    };
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &model.shell(sh).faces {
            let face = model.face(fh);
            // Exhaustive for the same reason as `points_move`: a new `Surface` variant (M6's
            // sphere/cone) must be a compile error here, not a silently skipped probe.
            match model.surface_cache(face.surface) {
                nacre_geom::Surface::Plane(pl) => probe(&mut ok, pl.origin()),
                // ★ The cylinder's axis origin is a datum a judgment reads since M6 — the gate's
                // clearance and `world_cylinder_def`'s postcondition compare it against the
                // cache — so it is probed like a plane's origin (the direction rides `dir_rat`,
                // exact under any turn the rationals can state).
                nacre_geom::Surface::Cylinder(cy) => probe(&mut ok, cy.axis().origin()),
            }
            for (k, (_, m)) in candidates.iter().enumerate() {
                if ok[k] && !points_move(m, face.surface) {
                    ok[k] = false;
                }
            }
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in model.edge(he.edge).vertices.iter() {
                        probe(&mut ok, model.vertex_point(vh));
                    }
                }
            }
        }
    }
    candidates
        .iter()
        .zip(&ok)
        .find(|(_, o)| **o)
        .map_or(Carry::None, |((c, _), _)| *c)
}

/// A plane's exact triple carried through the part of a motion the statements absorb — the one
/// transport both `carry_of`'s probe and pass 1 use, so the probe's feasibility answer and the
/// actual transport cannot disagree. `None` on `i128` overflow, which the probe turns
/// into "record a node instead".
fn transport_points(motion: &Xform<'_>, p: [[Rat; 3]; 3]) -> Option<[[Rat; 3]; 3]> {
    let each = |q: [Rat; 3]| -> Option<[Rat; 3]> {
        match motion {
            Xform::Rigid(iso) => iso.point_rat(q),
            Xform::Mirror { axis, offset, .. } => nacre_exact::mirror_point_rat(q, *axis, *offset),
        }
    };
    Some([each(p[0])?, each(p[1])?, each(p[2])?])
}

/// The cylinder twin of [`transport_points`] — pass 1 transports with the very function the
/// probe checked, so "this will not overflow" is structural, not a parallel re-derivation.
///
/// Rigid: the origin rides `point_rat`; the two directions ride the rotation alone
/// (`dir_rat` — a direction is a difference of points, so pivot and translation cancel); the
/// radius is invariant under a rigid motion. The image is rebuilt through the checked
/// constructor, which a rigid motion cannot fail except on overflow — then `None`, the
/// conservative answer the probe turns into "take the recorded path".
///
/// Mirror: `None`, deliberately — pass 1's `MirrorNotPlanar` rejection preempts (a mirrored
/// cylinder has no cache image, `Xform::surface`), so this arm only ever steers the probe and
/// must refuse quietly rather than panic.
fn transport_cylinder(
    motion: &Xform<'_>,
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_topo::CylinderDef> {
    match motion {
        Xform::Rigid(iso) => nacre_topo::CylinderDef::new(
            iso.point_rat(def.origin())?,
            iso.dir_rat(def.dir())?,
            iso.dir_rat(def.ref_dir())?,
            def.r2().clone(),
        ),
        Xform::Mirror { .. } => None,
    }
}

/// The motion history a moved surface's image carries — the surface twin of the vertex
/// definition rules, held in the surface's own `motion` field:
///
/// | source motion | `Carry::None` — the whole motion recorded | `Carry::Rotation` — the turn carried, the translation recorded | `Carry::Full` — nothing recorded |
/// |---|---|---|---|
/// | `None` (world) | `Some(leaf)` — the points stay pre-motion | `Some(leaf)` — the points are **turned**, the node says the translation | `None` (points transported exactly) |
/// | `Some(m)` | `Some(leaf)`, hanging off `m` — replaying from the root applies every motion once | as `Carry::None` (the exceptions lapse once a history exists) | as `Carry::None` |
///
/// ★ The middle column is the law `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` (cell
/// ④): a recorded node states the plane **before the recorded part**, not before the motion.
///
/// The `None (world)` row has a third outcome decided **before** this function is asked: a
/// rigid motion that *fixes the plane* (`Isometry::fixes_plane` — pass 1's `invariant` flag)
/// records nothing and carries the points **verbatim**, so the image interns back onto the
/// source handle. Per plane, not per solid — the caps of a turned block take it while the
/// walls land in this table.
///
/// **One table for all three motions.** A reflection appends a node like everything else, so the
/// reflection is *in* the definition rather than folded into it.
///
/// A cylinder rides the same rows through its own motion slot: its history is recorded.
fn moved_surface_motion(
    model: &mut Model,
    src: Handle<Surface>,
    motion: &Xform<'_>,
    carry: Carry,
    surf_rot: &mut HashMap<Option<Handle<MotionNode>>, Option<Handle<MotionNode>>>,
) -> Option<Handle<MotionNode>> {
    // **Each surface chains from its own leaf, not the solid's.** One solid does not have one
    // history: a boolean between differently-moved operands hands back walls that came from
    // different ones. Memoized per distinct parent so surfaces that did share a history still do.
    let parent = model.plane_motion(src);
    let leaf = match surf_rot.get(&parent) {
        Some(&h) => h,
        None => {
            // `entry` cannot hold a `&mut Model` across the closure, so look up then insert.
            let h = chain_motion(model, parent, motion, carry);
            surf_rot.insert(parent, h);
            h
        }
    };
    // Nothing recorded: the motion kept the data exact and the image keeps the source's own
    // history. (`leaf` is always `Some` when `parent` is — the exceptions lapse once a history
    // exists — so `or` never resurrects a stale parent past a recorded node.)
    leaf.or(parent)
}

/// How [`transform_solid`] maps a solid's cells. One walker serves both kinds so the seven
/// passes are not duplicated; the kinds differ in exactly three places — how a point/direction
/// maps, whether a curved surface can be carried at all, and whether loops must be rewound.
pub(crate) enum Xform<'a> {
    /// A proper motion: rotation then translation. Preserves handedness.
    Rigid(&'a Isometry),
    /// A reflection in a coordinate plane. Reverses handedness, so `det = −1`.
    ///
    /// The plane is carried **twice**: `m` is the `f64` map the coordinates actually go through,
    /// `(axis, offset)` the exact statement of the same plane that the motion history records.
    /// One variant carries both, so neither half can arrive without the other (as separate
    /// arguments, one of them an `Option`, they could). [`Xform::mirror`] is the only
    /// constructor, so the two halves cannot be made to disagree either.
    Mirror {
        m: AxisMirror,
        axis: Axis,
        offset: Rat,
    },
}

impl Xform<'_> {
    /// A direction carried through the motion's linear part, exactly — `None` where the motion
    /// has no exact statement (an irrational turn).
    fn dir_rat(&self, v: [Rat; 3]) -> Option<[Rat; 3]> {
        match self {
            Xform::Rigid(iso) => {
                let moved = iso.point_rat(v)?;
                let origin = iso.point_rat([Rat::from_int(0); 3])?;
                Some([
                    moved[0].checked_sub(origin[0])?,
                    moved[1].checked_sub(origin[1])?,
                    moved[2].checked_sub(origin[2])?,
                ])
            }
            Xform::Mirror { axis, .. } => {
                let mut out = v;
                out[axis.index()] = Rat::from_int(0).checked_sub(v[axis.index()])?;
                Some(out)
            }
        }
    }

    /// The reflection in `axis = offset`, with its `f64` map **derived** from the exact plane
    /// rather than handed in beside it.
    fn mirror(axis: Axis, offset: Rat) -> Xform<'static> {
        Xform::Mirror {
            m: AxisMirror::new(axis.index(), offset.to_f64()).expect("axis index is 0..3"),
            axis,
            offset,
        }
    }

    fn point(&self, p: Point3) -> Point3 {
        match self {
            Xform::Rigid(iso) => Point3::from_array(iso.apply_point(p.as_array())),
            Xform::Mirror { m, .. } => m.point(p),
        }
    }

    /// `None` when the variant has no image under this motion — a mirrored cylinder, whose
    /// parametrisation handedness is a curved-geometry decision (see `Surface::mirrored`).
    fn surface(&self, s: &nacre_geom::Surface, offset: Vector3) -> Option<nacre_geom::Surface> {
        match self {
            Xform::Rigid(iso) => Some(transform_surface(s, iso, offset)),
            Xform::Mirror { m, .. } => s.mirrored(*m),
        }
    }

    /// A reflection negates the normal a loop's winding implies (`R(a) × R(b) = −R(a × b)`), so
    /// every loop is rewound to put it back — and then the `Orientation` flag needs no change,
    /// because a reflection preserves dot products.
    fn reverses_orientation(&self) -> bool {
        matches!(self, Xform::Mirror { .. })
    }

    /// The isometry, for the rotation-forest bookkeeping that only proper motion does.
    fn rigid(&self) -> Option<&Isometry> {
        match self {
            Xform::Rigid(iso) => Some(iso),
            Xform::Mirror { .. } => None,
        }
    }
}

/// A surface moved by `isometry`. A pure translation uses `translated` (the normal —
/// and its exact `raw` — is unchanged). A rotation rebuilds from the moved
/// origin/normal via the constructor (rotation makes `raw` irrational, as expected —
/// the plane then carries tol, judged by CIP later).
fn transform_surface(
    s: &nacre_geom::Surface,
    iso: &Isometry,
    offset: Vector3,
) -> nacre_geom::Surface {
    if iso.rotate.is_none() {
        return s.translated(offset);
    }
    let p = |q: Point3| Point3::from_array(iso.apply_point(q.as_array()));
    let d = |v: Vector3| Vector3::from_array(iso.apply_dir(v.as_array()));
    match s {
        nacre_geom::Surface::Plane(pl) => nacre_geom::Surface::Plane(
            Plane::from_point_normal(p(pl.origin()), d(pl.normal()))
                .expect("rotation preserves a nonzero normal"),
        ),
        nacre_geom::Surface::Cylinder(cy) => {
            let ax = cy.axis();
            nacre_geom::Surface::Cylinder(
                Cylinder::from_axis(
                    p(ax.origin()),
                    d(ax.direction()),
                    d(cy.ref_dir()),
                    cy.radius(),
                )
                .expect("rotation preserves a valid cylinder"),
            )
        }
    }
}

/// Clone `solid` into a new solid with every cell's geometry mapped by `motion`,
/// preserving topology, shared cells (surfaces/curves/vertices/edges are deduped),
/// inner-loop holes, cavity shells, and each vertex's definition (its plane handles are
/// remapped to the moved surfaces). Cells are pushed in a
/// **deterministic traversal order** (shell → face → loop) with per-cell dedup maps, so the same
/// op-log reproduces identical handles (replay determinism).
///
/// Face orientation flags are carried unchanged for **both** kinds of motion. A rigid motion
/// turns the normal and the winding together; a reflection negates the winding's implied normal,
/// which [`Xform::reverses_orientation`] undoes by rewinding every loop — and since a reflection
/// preserves dot products, `sign(plane.normal · n_out)` is then unchanged too, which is exactly
/// what the `Orientation` flag records.
///
/// `Err` only when the motion has no image for some cell's geometry (a mirrored cylinder).
/// **Does the moved pierce pair's meet line run the other way?** `ℓ = n₁ × n₂` of the new pair
/// (in its stored, ascending order) against the old pair's `ℓ` carried through the motion as a
/// direction — exact, in the world names. `None` when a name is not narrow or the motion has no
/// exact statement; the caller then keeps the order-only restatement.
fn pierce_line_reversed(
    model: &Model,
    old: [Handle<Surface>; 2],
    new: [Handle<Surface>; 2],
    motion: &Xform<'_>,
) -> Option<bool> {
    let normal = |h: Handle<Surface>| -> Option<[Rat; 3]> {
        let c = *model.world_plane_name(h)?.narrow()?;
        Some([c[0], c[1], c[2]])
    };
    let cross = |a: [Rat; 3], b: [Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            a[1].checked_mul(b[2])?
                .checked_sub(a[2].checked_mul(b[1])?)?,
            a[2].checked_mul(b[0])?
                .checked_sub(a[0].checked_mul(b[2])?)?,
            a[0].checked_mul(b[1])?
                .checked_sub(a[1].checked_mul(b[0])?)?,
        ])
    };
    let l_old = cross(normal(old[0])?, normal(old[1])?)?;
    let l_new = cross(normal(new[0])?, normal(new[1])?)?;
    let carried = motion.dir_rat(l_old)?;
    let dot = l_new[0]
        .checked_mul(carried[0])?
        .checked_add(l_new[1].checked_mul(carried[1])?)?
        .checked_add(l_new[2].checked_mul(carried[2])?)?;
    // The same line either way, so a zero here is a contradiction, not an answer.
    if dot.numer() == 0 {
        return None;
    }
    Some(dot.numer() < 0)
}

fn transform_solid(
    model: &mut Model,
    solid: Handle<Solid>,
    motion: &Xform<'_>,
) -> Result<Handle<Solid>, OpError> {
    let offset = motion
        .rigid()
        .map(|iso| Vector3::from_array(iso.offset_f64()))
        .unwrap_or_else(Vector3::zero);
    let src = model.solid(solid).clone();

    // Decided once, for the whole solid — see [`carry_of`].
    let carry = carry_of(model, solid, motion);
    // The part of the motion pass 1 transports into the statements — `None` when the whole
    // motion is recorded and the statements stay pre-motion.
    let turn_only: Option<Isometry> = match (carry, motion) {
        (Carry::Rotation, Xform::Rigid(iso)) => iso.rotate.map(Isometry::rotation),
        _ => None,
    };
    let carry_xform: Option<Xform<'_>> = match carry {
        Carry::Full => Some(match motion {
            Xform::Rigid(iso) => Xform::Rigid(iso),
            Xform::Mirror { m, axis, offset } => Xform::Mirror {
                m: *m,
                axis: *axis,
                offset: *offset,
            },
        }),
        Carry::Rotation => turn_only.as_ref().map(Xform::Rigid),
        Carry::None => None,
    };
    // Whether the source already carries a motion history (any face surface's truth records
    // one). There is no vertex-side motion (the faces record it themselves, and a
    // vertex's motion is its faces'), so this is the whole of the
    // "exceptions lapse once the solid has a history" test: an exact move of a fresh solid
    // records nothing; anything else is recorded.
    // Deterministic order: outer shell then cavities; each shell's faces in order.
    let shell_order: Vec<Handle<Shell>> = std::iter::once(src.outer)
        .chain(src.cavities.iter().copied())
        .collect();
    let face_order: Vec<Handle<Face>> = shell_order
        .iter()
        .flat_map(|&sh| model.shell(sh).faces.clone())
        .collect();

    // Pass 1 — surfaces (dedup, moved): needed before the vertex definitions are remapped.
    //
    // A moved surface's coefficients are only the truth while the motion kept them exact, so each
    // one states its provenance here (its `motion`). The witness for a rotation is the face's own
    // **pre-rotation** triangle, read on the spot — the definition must not depend on anything
    // else in the model surviving, and those points are exactly the ones a rotated operand's
    // `tri_pt3` is built from.
    //
    // **A surface chains from its own leaf, not from the solid's.** One solid does not have one
    // rotation history: a boolean between differently-rotated operands hands back a result whose
    // walls came from different rotations, so there is no solid-wide rotation to answer for it.
    // `surf_rot` memoizes one new forest
    // node per distinct parent leaf, so surfaces that did share a history still share it.
    let mut surf_rot: HashMap<Option<Handle<MotionNode>>, Option<Handle<MotionNode>>> =
        HashMap::new();
    let mut surf_map: HashMap<Handle<Surface>, Handle<Surface>> = HashMap::new();
    // ★ Set when the surface the model handed back points the other way from the one built here.
    // A copied face keeps its `Orientation` because its surface moved with it — but a *shared*
    // surface did not, so the same outward direction has to be spelled the other way. This is a
    // different field from the mirror's rewind below, which turns loop *winding*; the two do not
    // interact (`Loop::reversed`'s doc spells out why a reflection touches only the winding).
    let mut surf_flip: HashMap<Handle<Surface>, bool> = HashMap::new();
    for &fh in &face_order {
        let s = model.face(fh).surface;
        if surf_map.contains_key(&s) {
            continue;
        }
        let moved = motion
            .surface(model.surface_cache(s), offset)
            .ok_or(OpError::MirrorNotPlanar)?;
        let src_truth = model.surface(s).clone();
        // ★ **A motion that fixes this plane restates nothing — the source statement already
        // states the image.** The per-plane sibling of `carry_of`'s whole-solid
        // decision: a rigid motion mapping this plane onto itself *as a set* (axis ∥
        // normal, in-plane translation — `Isometry::fixes_plane`, exact) leaves the same
        // points and no node to record, and re-pushing that statement interns back onto
        // the **source handle** — the road `Copy` already takes for the identity. The caps
        // of a solid turned about their own normal are the population; their frames stay
        // world-spoken (`face_sketch_frame` verifies instead of declining) and their
        // predicates keep the exact roads.
        //
        // `Known` + narrow name only: a `Through` plane's truth is vertex handles in a
        // separate intern table, and a source that already carries a motion keeps its
        // recorded path — which costs a fixed plane nothing, since it never gains a history
        // and a second turn about the same normal restates again. A mirror never restates here.
        let invariant = matches!(
            &src_truth,
            nacre_topo::Surface::Plane {
                points: nacre_topo::PlanePoints::Known(_),
                motion: None,
            }
        ) && motion.rigid().is_some_and(|iso| {
            model
                .surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| iso.fixes_plane(c))
        });
        let new_motion = if invariant {
            None
        } else {
            moved_surface_motion(model, s, motion, carry, &mut surf_rot)
        };
        // ★★★★★ **Only the points move.** The image's canonical name is derived from them by
        // `Model::push_plane`, so there is no second description to keep in step.
        //
        // A recorded node states its plane **before** the motion, and the image's base is the
        // source's base — `moved_surface_motion` chains from the source's own leaf for the same
        // reason — so the triple is inherited verbatim. Otherwise the carried part of the motion
        // moves the points in the world; the transport is the very function `carry_of` probed,
        // so it cannot fail here.
        // ★★ A `Through` plane's statement is *handles*, and handles do not move. Its image keeps
        // the same three vertices and gains the node — the definition composes as "the plane
        // through those, then this motion". Transporting them is not an option (they may not even
        // belong to this solid), and leaving them without a node would move the cache while the
        // truth stayed put. `points_move` refuses the no-node path for exactly this reason, so
        // `new_motion` is always `Some` here.
        let (new_s, flipped) = match (moved, &src_truth) {
            (
                nacre_geom::Surface::Plane(pl),
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Known(p),
                    ..
                },
            ) => {
                // ★ `invariant` must gate first: there `new_motion` and the source's motion are
                // both `None`, and
                // the else arm's expect names a probe (`carry_of`) that never ran on
                // this road — an irrational turn would panic in `transport_points`. The
                // statement is the image's verbatim, so there is nothing to transport.
                // ★ Verbatim when the statement already has a history (the node hangs off it
                // and says the plane before this motion) or when nothing of the motion is
                // carried; otherwise the carried part moves the triple — and only that part,
                // so a turn the statements absorbed is never described twice.
                let carried = if invariant || model.plane_motion(s).is_some() {
                    *p
                } else {
                    match &carry_xform {
                        Some(cx) => transport_points(cx, *p).expect("probed by carry_of"),
                        None => *p,
                    }
                };
                model.push_plane(pl, carried, new_motion)
            }
            (
                nacre_geom::Surface::Plane(pl),
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Through(vs),
                    ..
                },
            ) => {
                // ★ `new_motion` may be `None`, and only for a motion that moves nothing —
                // `Copy` is `transform_solid` with the identity. Then the plane is unchanged, the
                // same statement re-pushed interns back onto the source handle, and that is the
                // right answer. `points_move` refusing the transporting path is what rules out
                // the other case, where a node is missing because the walk thought it could carry
                // points that do not exist.
                //
                // ★★ **A second road to `None` exists and has not been reached
                // here — recorded rather than assumed away.** The `invariant` test above reads
                // the plane's *name*, and the fixed-carrier
                // licence gives a turned solid's `Through` datums one. So `invariant = true` is
                // possible for this arm. It needs the datum's plane to be
                // fixed by the motion — for a rotation that means its normal lies along the
                // axis — **and** to have stayed a `Through` truth.
                //
                // ★ **What is measured and what is argued, kept apart.** *Measured*: the assert
                // below did not fire anywhere in the suite, `replay`'s proptest included, and that
                // proptest does emit `Copy`/`Rotate` beside `DatumThroughVertices`. *Argued* (not
                // swept): for a cuboid every axis-normal plane through three of its vertices is a
                // face plane, so it interns onto that `Known` surface and never reaches this arm —
                // which would explain the silence, but no probe has confirmed it is the reason.
                // Shapes that ought to reach it: a stepped solid with three co-planar vertices off
                // any face, or a translation along a datum's own plane. Nobody has built one, so
                // "unreachable" is **not** what this says.
                let out = model.push_plane_through(pl, *vs, new_motion);
                debug_assert!(
                    new_motion.is_some() || out.0 == s,
                    "a Through plane gained no node yet changed handle — its truth would be \
                     describing the plane it used to be"
                );
                out
            }
            (nacre_geom::Surface::Cylinder(cy), nacre_topo::Surface::Cylinder { def, .. }) => {
                // The same fork as the `Known` plane above, minus the invariant road (an
                // invariant-cylinder restatement — a turn about its own axis — is deliberately
                // deferred; the condition is narrower than a plane's because `ref_dir` turns).
                // The same sentence as the plane's: verbatim behind a history or when nothing
                // is carried; otherwise the carried part moves the def (the very transport the
                // probe checked).
                let carried = if model.plane_motion(s).is_some() {
                    def.clone()
                } else {
                    match &carry_xform {
                        Some(cx) => transport_cylinder(cx, def).expect("probed by carry_of"),
                        None => def.clone(),
                    }
                };
                (model.push_cylinder(cy, carried, new_motion), false)
            }
            (nacre_geom::Surface::Plane(_), nacre_topo::Surface::Cylinder { .. })
            | (nacre_geom::Surface::Cylinder(_), nacre_topo::Surface::Plane { .. }) => {
                unreachable!("a surface's cache and truth cannot disagree about its kind")
            }
        };
        surf_map.insert(s, new_s);
        surf_flip.insert(s, flipped);
    }

    // Edge order (deterministic dedup) — used by passes 2/3/4.
    let mut edge_order: Vec<Handle<Edge>> = Vec::new();
    let mut edge_seen: HashSet<Handle<Edge>> = HashSet::new();
    for &fh in &face_order {
        let face = model.face(fh).clone();
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                if edge_seen.insert(he.edge) {
                    edge_order.push(he.edge);
                }
            }
        }
    }

    // (There is no pass 2 — curves are not stored: the moved edge's
    // curve derives from its moved carriers and endpoints in pass 4's `push_edge`.)

    // Pass 3 — vertices (dedup): the definition's handles re-pointed onto the moved surfaces,
    // the coordinate moved as a fallback figure, the point re-realized from the definition.
    let mut vert_order: Vec<Handle<Vertex>> = Vec::new();
    let mut vert_seen: HashSet<Handle<Vertex>> = HashSet::new();
    for &eh in &edge_order {
        for vh in model.edge(eh).vertices {
            if vert_seen.insert(vh) {
                vert_order.push(vh);
            }
        }
    }
    let mut vert_map: HashMap<Handle<Vertex>, Handle<Vertex>> = HashMap::new();
    for &vh in &vert_order {
        // `defs_are_remappable` cleared every definition handle before the walk began.
        let remap = |s: Handle<Surface>| {
            *surf_map
                .get(&s)
                .expect("a definition's surface must be a face surface of the solid")
        };
        let def = match *model.vertex(vh) {
            Vertex::ThreePlane(planes) => Vertex::ThreePlane(planes.map(remap)),
            Vertex::OnSeam(pair) => Vertex::OnSeam(pair.map(remap)),
            // ★ `.map(remap)` alone would be wrong here: pass 1 issues new surface handles in
            // face-traversal order, so the two planes' handle order can invert, and the root is
            // defined against the meet line of the *stored* order. `QuadRoot::canonical` owns
            // that restatement — including the tangency, whose single point a swap fixes — and
            // this site reads it rather than spelling it again.
            Vertex::Pierce {
                planes: [p0, p1],
                cylinder,
                root,
            } => {
                // ★★★ **The sign half.** `canonical` counts
                // the swap; but a restatement may spell a moved plane with the *opposite* normal
                // (its world name is canonical over four coefficients), and a reflection reverses
                // every cross product — each reverses `ℓ = n₁ × n₂`, and `Lo`/`Hi` trade names
                // once per reversal (`QuadRoot::canonical`'s doc: *"flip once per reversal"*).
                // Decided exactly from the world names: the new pair's `ℓ` against the old pair's
                // `ℓ` carried through the motion, which folds the swap in as well. Measured, by the
                // commutation oracle: counting the swap alone, a boss turned 90° names the *other*
                // crossing on 18 quadrantal cells, and a stored `f64` hides the wrong label.
                //
                // Unavailable (a chain the world cannot state exactly) means nothing was restated
                // — the spellings travelled whole — so the order-only answer stands.
                let new_pair = {
                    let [a, b] = [remap(p0), remap(p1)];
                    if b < a { [b, a] } else { [a, b] }
                };
                let root = match pierce_line_reversed(model, [p0, p1], new_pair, motion) {
                    Some(true) => root.flipped(),
                    Some(false) => root,
                    None => nacre_topo::QuadRoot::canonical([remap(p0), remap(p1)], root).1,
                };
                Vertex::Pierce {
                    planes: new_pair,
                    cylinder: remap(cylinder),
                    root,
                }
            }
        };
        let coord = motion.point(model.vertex_point(vh));
        // ★ **The moved figure is a fallback, and only a fallback.** A *measured residual* would
        // survive a rigid motion, but the cache holds a *bound*, which does not — so there is
        // nothing to carry across: this coordinate is `f64` arithmetic on the old one, never the
        // nearest `f64` of the moved definition, so "no realization stands behind it" is the whole
        // truth about it.
        // The realization of the moved definition runs next and replaces it wherever it answers.
        let cache = PointCache::Unrealized { coord };
        // ★ The one site that extends a chain: this vertex *is* the image of `vh`, so the prefix
        // `vh` already folded is exactly what this realization would otherwise walk again.
        vert_map.insert(
            vh,
            crate::realize::push_vertex_realized(
                model,
                def,
                cache,
                crate::realize::ChainLink::Extends,
            ),
        );
    }

    // Pass 4 — edges (carrier/vertex handles remapped; the curve cache derives from them).
    let mut edge_map: HashMap<Handle<Edge>, Handle<Edge>> = HashMap::new();
    for &eh in &edge_order {
        let e = *model.edge(eh);
        // The carriers move with the surfaces (pass 1 mapped every reachable one, so the
        // lookups cannot miss); `push_edge` re-canonicalizes the pair. A rigid image of a
        // non-degenerate edge cannot degenerate, so the `None` is unreachable in practice —
        // mapped to the honest reject rather than a panic all the same.
        let new_e = model
            .push_edge(
                [surf_map[&e.surfaces[0]], surf_map[&e.surfaces[1]]],
                e.vertices.map(|v| vert_map[&v]),
            )
            .ok_or(OpError::DegenerateGeometry)?;
        edge_map.insert(eh, new_e);
    }

    // Pass 5 — faces (loops rebuilt onto the new edges; orientation unchanged).
    // A reflection rewinds every loop (see the fn doc); a rigid motion keeps the winding.
    let rewind = motion.reverses_orientation();
    let map_loop = |lp: &Loop| {
        let mapped = Loop {
            half_edges: lp
                .half_edges
                .iter()
                .map(|he| HalfEdge {
                    edge: edge_map[&he.edge],
                    forward: he.forward,
                })
                .collect(),
        };
        if rewind { mapped.reversed() } else { mapped }
    };
    let mut face_map: HashMap<Handle<Face>, Handle<Face>> = HashMap::new();
    for &fh in &face_order {
        let face = model.face(fh).clone();
        let new_f = Face {
            surface: surf_map[&face.surface],
            outer: map_loop(&face.outer),
            inner: face.inner.iter().map(&map_loop).collect(),
            orientation: if surf_flip[&face.surface] {
                face.orientation.flipped()
            } else {
                face.orientation
            },
        };
        face_map.insert(fh, model.push_face(new_f));
    }

    // Pass 6 — shells; Pass 7 — solid.
    let mut shell_map: HashMap<Handle<Shell>, Handle<Shell>> = HashMap::new();
    for &sh in &shell_order {
        let faces: Vec<Handle<Face>> = model
            .shell(sh)
            .faces
            .iter()
            .map(|fh| face_map[fh])
            .collect();
        shell_map.insert(sh, model.push_shell(Shell { faces }));
    }
    let new_solid = Solid {
        outer: shell_map[&src.outer],
        cavities: src.cavities.iter().map(|sh| shell_map[sh]).collect(),
    };
    Ok(model.push_solid(new_solid))
}

#[cfg(test)]
#[path = "tests/transform.rs"]
mod tests;
