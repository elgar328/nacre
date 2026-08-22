//! The plane/face substrate: per-face (`FaceInfo`) and per-plane-class (`WorkingPlane`) tables and
//! their construction. Everything the boolean engine and its combinatorial queries build on.

use crate::combinatorics;
use crate::{BoolError, RejectReason, he_start, reject};
use nacre_cip::predicate::{Judge, Notes};
use nacre_cip::{Standard, WitnessPoint};
use nacre_geom::intersect::{plane_plane, planes_coplanar};
use nacre_geom::{Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_scalar::Mag;
use nacre_store::Handle;
use nacre_topo::{Edge, Face, HalfEdge, Model, Shell, Solid, Vertex};
use std::collections::HashMap;

/// A face's supporting plane plus the exact in/out data the seam path needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` is the face's *stated*
/// outward — `plane.normal()` × `orientation`, the reading props and STEP trust
/// — and the loop's winding is held to it by a `debug_assert` at construction
/// and by `validate`'s `FaceMisoriented` at every op. Every sign test here reads
/// `n_out` (or `tri`), and the two agree by that enforcement.
/// One row of the boolean's face table — the face vocabulary the engine reads (M6-2a).
///
/// The table used to be `Vec<FaceInfo>` with a hard `CylinderFace` reject at the door; the row
/// is now an enum so a cylinder face can *sit in the table* (keeping the face-index space that
/// `surf_ix`/`EdgeFaces`/`plane_ix` share) while the plane data keeps its own struct — plane
/// consumers read through [`FaceRow::plane`], and the population gate decides what flows.
/// ★ The size gap is the `SurfaceTruth` trade taken again: planes dominate every table (a
/// prism is all planes; a cylinder contributes one lateral row), so boxing the plane data
/// would put an allocation and a pointer chase on the common row to shrink the rare one.
#[allow(clippy::large_enum_variant)]
#[derive(Clone)]
pub(crate) enum FaceRow {
    Plane(FaceInfo),
    Cylinder(CylFaceInfo),
}

impl FaceRow {
    /// The plane data of this row. **Panics on a cylinder row** — every caller sits behind the
    /// population gate or a kind filter (`loop_triples`/`trace_one` skip cylinder rows), so a
    /// cylinder here is an upstream filter bug, and a loud panic beats a silently wrong plane.
    #[inline]
    #[track_caller]
    pub(crate) fn plane(&self) -> &FaceInfo {
        match self {
            FaceRow::Plane(p) => p,
            FaceRow::Cylinder(c) => panic!(
                "a plane-only path reached a cylinder row (face {:?}) — upstream filter bug",
                c.face
            ),
        }
    }

    /// The row's surface, whichever kind it is.
    // Test consumers today; the first production consumer is C2's population gate.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn surf(&self) -> Handle<Surface> {
        match self {
            FaceRow::Plane(p) => p.surf,
            FaceRow::Cylinder(c) => c.surf,
        }
    }

    /// The row's model face, whichever kind it is.
    #[inline]
    pub(crate) fn face(&self) -> Option<Handle<nacre_topo::Face>> {
        match self {
            FaceRow::Plane(p) => p.face,
            FaceRow::Cylinder(c) => c.face,
        }
    }
}

/// A cylinder face's table row (M6-2a): what `collect_planes` can state about a lateral face
/// without pretending it has a plane's four-piece description (`plane`/`tri`/`n_out`/`tri_pt3`
/// are constant-normal vocabulary — a dummy would be the type lying). The class index and the
/// exact def arrive with the cylinder class table (C2).
#[derive(Clone, Debug)]
pub(crate) struct CylFaceInfo {
    pub(crate) surf: Handle<Surface>,
    /// See [`FaceInfo::face`].
    pub(crate) face: Option<Handle<nacre_topo::Face>>,
    /// The `Orientation` flag as a sign — same reading as [`FaceInfo::orient_sign`]: `+1` when
    /// the stored surface normal (radially outward) is this face's outward.
    pub(crate) orient_sign: i8,
    /// The motion-history leaf of the cylinder's truth, `None` for a constructed one.
    pub(crate) motion: Option<Handle<nacre_topo::MotionNode>>,
    /// The cylinder's exact statement — cloned here so the tracer (which works off the face
    /// table, never the `Model`) can ask ⊥-ness and axis parameters.
    pub(crate) def: nacre_topo::CylinderDef,
    /// This lateral **face**'s span in the axis parameter `t` (of the raw `def.dir()`), read
    /// off its two rim carrier planes: `t = −(n·o + d)/(n·m)` per rim plane, ordered. `None`
    /// when a rim carrier has no narrow rational name — the transversal-circle producer then
    /// declines the face rather than guessing.
    pub(crate) span: Option<[nacre_scalar::Rat; 2]>,
}

#[derive(Clone)]
pub(crate) struct FaceInfo {
    pub(crate) surf: Handle<Surface>,
    /// The face this plane came from. Distinguishes two coplanar faces that share one
    /// `Surface` (a Cut splits one face into disjoint pieces reusing its surface —
    /// cell coplanar-narrow), which `surf` alone collapses. `surf_ix` keys on this.
    ///
    /// ★ **`None` has no producer today.** It was the cap a half-space clip put on an operand — a
    /// face that lives for one boolean and is never emitted — and the subdivision that minted those
    /// is gone (see `docs/dev-log.md`). Every constructor writes `Some`, so every read `expect`s.
    /// Kept as an `Option` because the table is the natural home for a face the model does not own,
    /// and the next engine that needs one should not have to re-thread the type.
    pub(crate) face: Option<Handle<Face>>,
    pub(crate) plane: Plane,
    /// Three non-collinear outer-loop points, **ordered so their RH normal is outward**.
    /// The order need not follow the loop: at a reflex corner it is reversed.
    pub(crate) tri: [Point3; 3],
    /// Outward normal — `plane.normal()` × the face's stated `orientation`; the single
    /// source of "outward" for both the in/out sign test and face ordering. Read off
    /// the b-rep's statement, never re-derived from loop geometry (the re-derivation's
    /// conditioning was the pad-eats-material defect).
    pub(crate) n_out: Vector3,
    /// `+1` when this face's stored plane normal already points out of its solid, `-1` when the
    /// face is `Reversed` and the two oppose — the `Orientation` flag as a sign.
    ///
    /// **This face's**, not its plane class's. The class-frame twin is [`WorkingPlane::frame_sign`],
    /// and the two used to be one function called with either kind of index — the single place the
    /// face/plane convention could not be asserted, because both readings were legitimate
    /// (dev-log, normalization cell). Separate names, separate questions.
    pub(crate) orient_sign: i8,
    /// The three `tri` points as **exact `WitnessPoint` definitions**, in the same order as `tri`.
    /// Built once here and borrowed by every predicate (`plane_def`) — it used to be rebuilt
    /// per judgment, which dominated the boolean's runtime.
    pub(crate) tri_pt3: [WitnessPoint; 3],
    /// The motion-history leaf this face's plane was moved by, or `None` for a constructed one.
    /// **The canonical identity of "which motion"** — see [`BaseFrame`].
    pub(crate) motion: Option<Handle<nacre_topo::MotionNode>>,
    /// The surface's plane as **exact rational coefficients in the frame its truth names**
    /// (`Model::surface_name`) — the world when unmoved, the pre-motion frame when moved.
    /// `None` when the producer had no rational description. Read by [`BaseFrame`], which would
    /// otherwise re-derive a moved plane from its pre-motion triangle and round `d`.
    pub(crate) base_rat: Option<[nacre_scalar::Rat; 4]>,
    /// The surface's full canonical name (`Model::surface_name`), **any width** — what
    /// [`WorkingPlane::name_ints`] is folded from. `base_rat` above is its narrow projection,
    /// kept beside it because the narrow consumers (`BaseFrame`, the composed-rotation route)
    /// read `[Rat; 4]` directly.
    pub(crate) name: Option<nacre_scalar::PlaneName>,
    /// Whether this face's plane is a *moved image* — the predicate-routing signal, read from
    /// the surface's own truth (`Model::surface_truth`).
    ///
    /// **Set together with `tri_pt3`, and only here.** It used to be decided per solid, by asking
    /// the vertices — which a boolean's result cannot answer, since its vertices are all
    /// `Discovered`. The surface answers for itself, and one solid can hold both kinds at once
    /// (fuse an axis-aligned hub with a turned fin).
    pub(crate) rotated: bool,
}

/// The supporting planes of a solid's outer shell. `Rejected` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<FaceRow>, BoolError> {
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let plane = match model.surface(face.surface) {
                Surface::Plane(p) => *p,
                // A cylinder face sits in the table (M6-2a) — its row keeps the shared facts
                // (surface, face, stated outward sign, motion leaf) and none of the plane
                // vocabulary. Whether it may *flow* is the population gate's question, asked
                // in `plane_index_setup`, not a door slam here.
                Surface::Cylinder(_) => {
                    let nacre_topo::SurfaceTruth::Cylinder { def, motion } =
                        model.surface_truth(face.surface)
                    else {
                        unreachable!("a cylinder cache carries a cylinder truth")
                    };
                    out.push(FaceRow::Cylinder(CylFaceInfo {
                        surf: face.surface,
                        face: Some(fh),
                        orient_sign: face.orientation.sign(),
                        motion: *motion,
                        span: lateral_axis_span(model, face, def),
                        def: def.clone(),
                    }));
                    continue;
                }
            };
            // **Outward is read off the face's statement, not re-derived from its loop.**
            // `orientation` relates the stored surface normal to "out of the solid" — the
            // same reading props and STEP already trust — so `n_out` is that product,
            // exact in direction by construction. The triangle's cross used to be the
            // source, and its conditioning was the pad-eats-material defect: on rotated
            // near-collinear corners the direction that came back was the rounding. The
            // winding is still consulted — as the cross-check below, not as the answer.
            let orient_sign = face.orientation.sign();
            let n_out = plane.normal() * f64::from(orient_sign);
            let tri = match outer_tri(model, face) {
                Some((tri, _)) => tri,
                // A **disk** face (a cylinder cap): its outer loop is one circle edge with a
                // single seam vertex, so no three loop points exist — the plane's own truth
                // points state the triangle instead (the cap's construction points,
                // realized), wound to this face's outward like every other `tri`.
                None if face.outer.half_edges.len() == 1
                    && matches!(
                        model.edge_curve(face.outer.half_edges[0].edge),
                        nacre_geom::Curve::Circle(_)
                    ) =>
                {
                    let nacre_topo::SurfaceTruth::Plane {
                        points: nacre_topo::PlanePoints::Known(pts),
                        ..
                    } = model.surface_truth(face.surface)
                    else {
                        return Err(reject(RejectReason::DegenerateFace));
                    };
                    let mut tri = pts.map(|p| Point3::from_array(p.map(|x| x.to_f64())));
                    let wound = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
                    if wound.dot(n_out) < 0.0 {
                        tri.swap(1, 2);
                    }
                    tri
                }
                None => return Err(reject(RejectReason::DegenerateFace)),
            };
            // **The plane's exact definition comes from the surface, not from the vertices.**
            //
            // Both used to be decided per *solid* ("is this solid rotated?"), which a boolean's
            // result cannot answer — it carries no rotation provenance, so every result face was
            // described by `WitnessPoint::exact` of its rounded triangle and one wall became two plane
            // classes on the next operation. The surface knows (its truth — points + motion,
            // S6b), and a result face reuses its operand's surface handle, so the answer
            // survives a chain of booleans.
            //
            // `tri_pt3` is an *oriented* plane witness, but the recorded triple belongs to the
            // *plane* — two faces sharing it can face opposite ways, and the implicit-point
            // `orient3d` reads the side `tri_pt3` spans. So both arms wind it to agree with
            // *this* face's `n_out`.
            let wind = |mut w: [WitnessPoint; 3]| -> [WitnessPoint; 3] {
                let e1 = Vector3::from_array(w[1].coord) - Vector3::from_array(w[0].coord);
                let e2 = Vector3::from_array(w[2].coord) - Vector3::from_array(w[0].coord);
                if e1.cross(e2).dot(n_out) < 0.0 {
                    w.swap(1, 2);
                }
                w
            };
            let (tri_pt3, rotated, motion) = match model.surface_truth(face.surface) {
                // ★★★★★ **The plane's own points state the plane — not the face's triangle.**
                //
                // The face's triangle is where this used to read from, and after a chain of
                // booleans those corners are `Discovered`: seam points the kernel itself annotated
                // with a tol, handed to `WitnessPoint::exact`, which states `tol = 0` **by construction**.
                // Measured, 1,938 of 68,350 triangles carried one, in 162 plane tables, 77 of
                // which also held a rotated plane — where the claim is actually consulted. The
                // plane's recorded points are construction points, so no such claim is made.
                //
                // ★ `WitnessPoint::exact` is kept for a point that **is** an f64: it states the same thing
                // and skips the nine BigFloat operations `WitnessPoint::at` spends measuring a zero. The
                // round-trip test is the one `BaseFrame` already uses.
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    motion: None,
                } => {
                    let w = pts.map(|b| {
                        let f = b.map(|r| r.to_f64());
                        match WitnessPoint::exact(f).filter(|_| {
                            b.iter()
                                .zip(f)
                                .all(|(&r, x)| nacre_scalar::Rat::try_from_f64(x) == Some(r))
                        }) {
                            Some(p) => p,
                            // ★★★ **The bound is free; measuring it is not.** `WitnessPoint::at` reads
                            // the rounding at 120 bits — nine BigFloat operations per point,
                            // and 38.3% of these points are not f64, so the suite paid 6% for
                            // it. `Rat::to_f64` is documented as *"the **nearest** f64 … ties
                            // to even"*, so `|r − to_f64(r)| ≤ ½ ulp` holds by its contract and
                            // `|x|·2⁻⁵³` is an upper bound on that for every normal `x`.
                            //
                            // ★ The cost is looseness: a filter interval slightly wider than
                            // the truth escalates slightly more often. Never unsound — a tol
                            // that overstates the error can only make a definite answer
                            // indefinite, never the other way round.
                            None => {
                                let bound = |x: f64| {
                                    (x.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE)
                                };
                                WitnessPoint::at_with_tol(
                                    b,
                                    [bound(f[0]), bound(f[1]), bound(f[2])],
                                )
                            }
                        }
                    });
                    (wind(w), false, None)
                }
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    motion: Some(motion),
                } => {
                    let motion = *motion;
                    let _t = Watch::new(); // charged at the arm's end
                    // The pre-motion description, carried through the recorded chain — the same
                    // computation, in the same order, that a moved vertex's `WitnessPoint` performs.
                    let chain = crate::rotated_vertex::motion_chain(model, motion)
                        .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                    let turn = |base: [nacre_scalar::Rat; 3]| -> Result<WitnessPoint, BoolError> {
                        crate::rotated_vertex::replay(WitnessPoint::at(base), &chain)
                            .ok_or_else(|| reject(RejectReason::FrameOutOfRange))
                    };
                    // ★★★★★ **The exact points the plane was built from — the only description.**
                    //
                    // This arm used to fall back to a `witness: [Point3; 3]` beside the motion: the
                    // same points *realized*, and a third of them do not survive the trip
                    // (measured 34.5% — `11/10` is not an f64, so lifting the realization back
                    // recovers a different rational). The judge then described the plane through
                    // three rounded points, which is how two caps that are one plane were told
                    // apart with full confidence.
                    let w = [turn(pts[0])?, turn(pts[1])?, turn(pts[2])?];
                    _t.charge(Sub::TriPt3);
                    (wind(w), true, Some(motion))
                }
                // ★★★ **A `Through` plane is solved into the same witness triangle here.**
                //
                // Its truth is handles, and the judging layer takes points — so the points are
                // *derived* at the boundary, once per plane per operation, exactly like the name
                // was derived once at push. That is not a rule-1 violation: this table lives for
                // one operation and is a mirror, not truth.
                //
                // The rational-closure branch is stage 1: three vertices that solve to `Rat`
                // give a triangle indistinguishable from a stated one, so every predicate below
                // runs unchanged.
                //
                // ★★★ **The nameless branch (open item 16, second wall — heterogeneous half).**
                // A pure-mixed datum's vertices are each exact *in their own frame*, so its
                // witness triangle exists — three `WitnessPoint`s whose chains simply differ.
                // The judging layer never required them to agree: `plane_iv`/`plane_hp` realize
                // each point independently, so the whole toleranced route (C4) runs unchanged,
                // and `standard_for` sizes the operation from these very points. What such a
                // plane has none of is exact f64 coefficients — so it is flagged `rotated`
                // (the routing signal means "no exact description", not "carries a motion"),
                // which makes every exact shortcut decline and `reconcile` carry nothing.
                //
                // A *straddling*-vertex datum has no witness triangle **from its vertices** —
                // its witness is the judged frame's probes, built in the branch below (16-3).
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Through(vs),
                    motion,
                } => {
                    let motion = *motion;
                    let _t = Watch::new();
                    let replay_all = |mut w: [WitnessPoint; 3],
                                      m: Handle<nacre_topo::MotionNode>|
                     -> Result<[WitnessPoint; 3], BoolError> {
                        let chain = crate::rotated_vertex::motion_chain(model, m)
                            .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                        for o in w.iter_mut() {
                            *o = crate::rotated_vertex::replay(o.clone(), &chain)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                        }
                        Ok(w)
                    };
                    let (w, rotated) = match model.through_points_rat(*vs) {
                        Some(base) => {
                            let w = base.map(WitnessPoint::at);
                            let w = match motion {
                                None => w,
                                Some(m) => replay_all(w, m)?,
                            };
                            (w, motion.is_some())
                        }
                        None => {
                            let j = crate::rotated_vertex::through_judged_points(model, *vs)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                            // ★ All-pure unwraps to 16-1's heterogeneous triangle, letter for
                            // letter — that road is locked by its own tests and stays.
                            //
                            // ★★★★ **Any implicit point: the witness is the plane's own judged
                            // frame** (16-3). The table's contract has been "three exact points
                            // *on the plane*, wound to n_out" since S6b — never "the face's
                            // corners" — and a judged plane has such points by definition: its
                            // canonical frame's probes `(0,0,0)·(1,0,0)·(0,1,0)`, the same ones
                            // `frame_world_basis` realizes. The origin is the foot of the
                            // perpendicular (on the plane exactly), û and v̂ are in-plane by
                            // construction, and `frame_chain` already appends the plane's own
                            // later motion — so the probes are exact definitions the escalation
                            // realizes at any precision. No `WorkingPlaneDef::Through` was ever
                            // needed; this is where that assumption died.
                            let all_pure = j
                                .iter()
                                .all(|p| matches!(p, nacre_cip::JudgedPoint::Pure(_)));
                            if all_pure {
                                let mut w: [Option<WitnessPoint>; 3] = [None, None, None];
                                for (o, jp) in w.iter_mut().zip(j) {
                                    if let nacre_cip::JudgedPoint::Pure(wp) = jp {
                                        *o = Some(wp);
                                    }
                                }
                                let w = w.map(|o| o.expect("all pure"));
                                let w = match motion {
                                    // ★ The plane's own later motion appends to each point's
                                    // own chain — motions add, never multiply (S5(ii)-1).
                                    None => w,
                                    Some(m) => replay_all(w, m)?,
                                };
                                (w, true)
                            } else {
                                let chain = crate::rotated_vertex::frame_chain(
                                    model,
                                    face.surface,
                                    &nacre_topo::FramePlacement::Canonical,
                                    false,
                                )
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                                let r = nacre_scalar::Rat::from_int;
                                let probe = |u: i128, v: i128| {
                                    crate::rotated_vertex::replay(
                                        WitnessPoint::at([r(u), r(v), r(0)]),
                                        &chain,
                                    )
                                    .ok_or_else(|| reject(RejectReason::FrameOutOfRange))
                                };
                                // `frame_chain` already carries the plane's own later motion in
                                // its tail, so no `replay_all` here — appending it twice would
                                // move the witness off the plane.
                                let w = [probe(0, 0)?, probe(1, 0)?, probe(0, 1)?];
                                (w, true)
                            }
                        }
                    };
                    _t.charge(Sub::TriPt3);
                    (wind(w), rotated, motion)
                }
                nacre_topo::SurfaceTruth::Cylinder { .. } => {
                    unreachable!("a plane cache cannot carry a cylinder truth")
                }
            };
            // The loop's winding, cross-checked against the flag it stated. `validate`
            // pins this same invariant as `FaceMisoriented`; here it runs at every
            // boolean in debug. A failure is a producer bug — a lying flag or a
            // mis-wound loop — never a conditioning artifact (`outer_tri` picks the
            // widest corner).
            //
            // ★ This assertion needs three points, so it says nothing about a **disk**
            // face (one closed rim) or about a face's **holes**. Those belong to the
            // same invariant and `validate` owns them — it reads a rim's circle and
            // every inner loop. Replicating that here would be a second spelling of one
            // rule, which is the shape this kernel keeps having to undo.
            debug_assert!(
                (tri[1] - tri[0])
                    .cross(tri[2] - tri[0])
                    .normalize()
                    .is_some_and(|w| w.dot(n_out) > 0.5),
                "a face's outer winding must agree with its stated orientation"
            );
            let name = model.surface_name.get(&face.surface).cloned();
            out.push(FaceRow::Plane(FaceInfo {
                base_rat: name.as_ref().and_then(|n| n.narrow()).copied(),
                name,
                surf: face.surface,
                face: Some(fh),
                plane,
                tri,
                n_out,
                orient_sign,
                tri_pt3,
                rotated,
                motion,
            }));
        }
    }
    // ★★★ **The mirror takes the strongest description — a restated plane rejoins its
    // solid's chain here.** The truth keeps a motion-fixed plane world-stated (the
    // invariant-plane restatement: one plane, one handle, no node), but this table is a
    // per-operation mirror, and a chain-borne description is strictly stronger when the
    // rest of the solid carries one: with every plane on one chain the shared-motion
    // roads answer exactly — a shared rotation must not turn exact questions into
    // assumed ones (`a_shared_rotation_still_assumes_nothing`). A fixed plane admits it
    // verbatim: its world points replayed through the chain still lie on the plane (the
    // chain maps the plane onto itself), and its world coefficients ARE its pre-motion
    // coefficients — bit-identical to the description the moved twin carried before the
    // restatement. Row-verbatim fixedness (`chain_preserves_plane_row`), not mere
    // set-fixedness: a sign-flipped row would poison determinant reads.
    //
    // The chain is chosen deterministically (lowest node index that qualifies) so replay
    // reproduces the table; a plane no chain fixes keeps its world description.
    let leaves = {
        let mut v: Vec<Handle<nacre_topo::MotionNode>> = out
            .iter()
            .filter_map(|r| match r {
                FaceRow::Plane(f) => f.motion,
                FaceRow::Cylinder(c) => c.motion,
            })
            .collect();
        v.sort_unstable_by_key(|h| h.index());
        v.dedup();
        v
    };
    if !leaves.is_empty() {
        for row in out.iter_mut() {
            // The restatement mirror is a plane story — a motion-fixed *cylinder*'s
            // restatement is separately deferred (M6-0's handover list).
            let FaceRow::Plane(f) = row else { continue };
            if f.motion.is_some() {
                continue;
            }
            let Some(name) = f.base_rat else { continue };
            let nacre_topo::SurfaceTruth::Plane {
                points: nacre_topo::PlanePoints::Known(pts),
                motion: None,
            } = model.surface_truth(f.surf)
            else {
                continue;
            };
            let Some(&c) = leaves
                .iter()
                .find(|&&c| model.chain_preserves_plane_row(c, &name))
            else {
                continue;
            };
            let Some(chain) = crate::rotated_vertex::motion_chain(model, c) else {
                continue;
            };
            let turned: Option<Vec<WitnessPoint>> = pts
                .iter()
                .map(|&b| crate::rotated_vertex::replay(WitnessPoint::at(b), &chain))
                .collect();
            let Some(w) = turned else { continue }; // conservative: keep the world description
            let mut w: [WitnessPoint; 3] = [w[0].clone(), w[1].clone(), w[2].clone()];
            let e1 = Vector3::from_array(w[1].coord) - Vector3::from_array(w[0].coord);
            let e2 = Vector3::from_array(w[2].coord) - Vector3::from_array(w[0].coord);
            if e1.cross(e2).dot(f.n_out) < 0.0 {
                w.swap(1, 2);
            }
            f.tri_pt3 = w;
            f.rotated = true;
            f.motion = Some(c);
        }
    }
    Ok(out)
}

/// A lateral face's axis-parameter span, read off its rim carrier planes (M6-2a).
///
/// Each rim edge's carrier pair is `[lateral, cap-plane]`; the cap plane meets the axis
/// `o + t·m` at `t = −(n·o + d)/(n·m)` — rational whenever the plane has a narrow name (its
/// ⊥-ness guarantees `n·m ≠ 0`). Two distinct rim planes give the span; anything else (a
/// nameless rim carrier, a non-⊥ rim, fewer or more than two distinct rims) answers `None`
/// and the consumer declines the face.
/// **Where a plane meets a cylinder's axis, as the axis parameter `t`** — the one spelling of
/// `t = −(n·o + d)/(n·m)` for `axis(t) = o + t·m`.
///
/// `None` when the plane is parallel to the axis (`n·m = 0`, no meeting point) or when the
/// checked `Rat` arithmetic overflows. ★ Three consumers ask this question — the lateral's rim
/// span, the transversal-circle test, and the band pass — and a rule that lives inlined in one
/// place while a second site spells a reduced version of it is this repo's dominant defect
/// shape, so it lives here once.
/// **Does an increasing axis parameter move toward this class's "above"?** — where "above" is the
/// side of the class's **stored** plane normal, which is the frame every cell label is written in
/// (`arrangement`'s `w_normal`).
///
/// ★ The `f64` dot is exact enough by construction: this is only ever asked of a class ⊥ to the
/// axis, so the dot is `±|n||m|` — a full magnitude from the sign boundary, not a near-zero
/// comparison.
///
/// ★★ **Ask this, do not re-derive it from the class's rational name.** `base_rat` is the
/// *canonical* name (first nonzero component positive), which points the other way from the stored
/// normal on half the classes; a rule spelled against it reads "above" backwards exactly there.
/// That mistake, made while adding the second consumer below, turned 36 tests red at once.
pub(crate) fn plus_t_is_above(wp: &WorkingPlane, def: &nacre_topo::CylinderDef) -> bool {
    let m = def.dir();
    let axis = Vector3::from_array([m[0].to_f64(), m[1].to_f64(), m[2].to_f64()]);
    wp.plane.normal().dot(axis) > 0.0
}

pub(crate) fn axis_param_of_plane(
    coeffs: &[nacre_scalar::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_scalar::Rat> {
    use nacre_scalar::Rat;
    let (o, m) = (def.origin(), def.dir());
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let dot3 = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    let nm = dot3(&n, &m)?;
    if nm == Rat::from_int(0) {
        return None;
    }
    let no_d = dot3(&n, &o)?.checked_add(coeffs[3])?;
    Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)
}

fn lateral_axis_span(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_scalar::Rat; 2]> {
    use nacre_scalar::Rat;
    let mut ts: Vec<Rat> = Vec::new();
    for he in &face.outer.half_edges {
        let e = model.edges.get(he.edge);
        let [a, b] = e.surfaces;
        let cap = if a == face.surface { b } else { a };
        if cap == face.surface {
            continue; // the seam edge is self-adjacent — not a rim
        }
        let coeffs = *model.surface_name.get(&cap)?.narrow()?;
        // `None` here is a rim carrier parallel to the axis (or overflow) — outside this
        // vocabulary either way.
        let t = axis_param_of_plane(&coeffs, def)?;
        if !ts.contains(&t) {
            ts.push(t);
        }
    }
    let [a, b] = ts[..] else { return None };
    Some(if a < b { [a, b] } else { [b, a] })
}

/// All shells of a solid — outer first, then cavities. The boolean seam
/// front-end walks these so a cavitied operand's void walls are seen (cell
/// (5c-in)); a non-hollow solid yields just its outer shell, unchanged.
pub(crate) fn solid_shell_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Shell>> {
    let s = model.solids.get(solid);
    std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect()
}

/// Three **well-spread** points of a face's outer loop — with the **vertex handle** each
/// point came from — ordered so their right-hand normal points **out** of the solid.
/// The handles let the toleranced predicates rebuild each point as a `WitnessPoint` (overhaul
/// stage 3); the coordinates alone drive the axis-aligned path.
///
/// "Well spread" rather than "non-collinear" is the whole contract: the triangle is what states
/// this face's outward direction *and* what `WorkingPlane::tri` carries into the predicates, so a
/// nearly-flat one is not a lesser answer but a wrong one. See the corner choice below.
pub(crate) fn outer_tri(model: &Model, face: &Face) -> Option<([Point3; 3], [Handle<Vertex>; 3])> {
    let verts: Vec<Handle<Vertex>> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let pts: Vec<Point3> = verts.iter().map(|&vh| model.vertex_point(vh)).collect();
    let n = pts.len();
    // The turn at one corner does not know which way the ring winds. Every b-rep loop is
    // CCW about its face's outward normal, but at a *reflex* corner the local turn
    // opposes the global winding, so three consecutive points can hand back an inward
    // normal. The Newell sum has no single corner to be fooled by.
    let newell = (0..n).fold(Vector3::zero(), |acc, i| {
        acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
    });
    // ★★★ **The widest corner, not the first non-flat one.** A loop carries vertices that do not
    // turn: two faces sharing an edge must list the same vertices along it, so a pad that splits a
    // neighbour's face leaves this loop with points strung along one straight line. Three of those
    // span **exactly** zero area — but only in exact arithmetic. On rotated coordinates the f64
    // cancellation leaves ~2⁻⁵³ instead, which passes a `> 0.0` gate and `normalize`, and the
    // direction that comes back is the rounding, not the plane: measured 90° off its own surface.
    //
    // **`0.0` is not a threshold in floating point** — "not zero" is not "well conditioned". The
    // sibling that answers this same question already says so: `exact::cap_points` takes "the
    // widest turn ... so a nearly-collinear pair is not chosen when a better one exists", and
    // records there why the f64 realization is the right instrument for a *selection* (the points
    // kept are exact; only "which three are spread out" is being asked, and answering it in
    // rationals would risk an `i128` overflow for nothing).
    //
    // `None` when every corner is degenerate, exactly as before — that is `DegenerateFace`.
    let (i, best) = (0..n).fold((0usize, 0.0f64), |acc, i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        let spread = (b - a).cross(c - a).norm();
        if spread > acc.1 { (i, spread) } else { acc }
    });
    if best <= 0.0 {
        return None;
    }
    let (i0, i1, i2) = (i, (i + 1) % n, (i + 2) % n);
    let (a, b, c) = (pts[i0], pts[i1], pts[i2]);
    // Same b/c swap for coords and handles, so `tri[k]` and `tri_verts[k]` stay aligned.
    Some(if (b - a).cross(c - a).dot(newell) < 0.0 {
        ([a, c, b], [verts[i0], verts[i2], verts[i1]])
    } else {
        ([a, b, c], [verts[i0], verts[i1], verts[i2]])
    })
}

/// Max distance of `p` to its 3 planes and 3 pairwise lines (the measured
/// vertex cache's measured tolerance).
pub(crate) fn vertex_tol(p: Point3, a: &Plane, b: &Plane, c: &Plane) -> f64 {
    let mut tol = a.distance(p).max(b.distance(p)).max(c.distance(p));
    for (x, y) in [(a, b), (a, c), (b, c)] {
        if let Some(line) = plane_plane(x, y) {
            tol = tol.max(line.distance(p));
        }
    }
    tol
}

/// [`vertex_tol`]'s branch sibling — the measured tolerance of a realized `plane ∩ plane ∩
/// cylinder` vertex: max distance of `p` to its two planes, the cylinder **surface**, and the
/// planes' meet line.
///
/// ★★ **The meet line is the one pairwise curve this covers, on purpose.** The seam's tolerance
/// means "surfaces *and* pairwise meets" — `boolean`'s vertex minting says the pairwise part "is
/// a real part of what this number means", and [`vertex_tol`] above covers all three of its
/// lines. Here only `a ∩ b` has a closed form that is always a line: a branch plane can be
/// **parallel to the axis**, where `plane ∩ cylinder` is a pair of ruling lines — a per-case
/// curve family this deliberately does not chase. The cylinder-*surface* distance already bounds
/// the radial part of that error; the meet line is also exactly the line the point is defined on
/// (`combinatorics::branch_point` realizes from `(line, s)`), so its residual is the first-class
/// question about the realization.
pub(crate) fn branch_vertex_tol(
    p: Point3,
    a: &Plane,
    b: &Plane,
    cyl: &nacre_geom::Cylinder,
) -> f64 {
    let mut tol = a.distance(p).max(b.distance(p)).max(cyl.distance(p));
    if let Some(line) = plane_plane(a, b) {
        tol = tol.max(line.distance(p));
    }
    tol
}

/// Which operand a face or segment came from — the boolean's per-cell label needs both solids'
/// material above and below, so provenance cannot be merged away.
///
/// It lives with the plane table rather than with the arrangement because it labels **an
/// operand**, not a segment: this is where the two solids first share one index space, so it is
/// where "whose is this?" first has an answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SolidSide {
    A,
    B,
}

/// Whose faces each plane class carries: `Some(side)` when every face in the class is that
/// operand's, `None` when both operands have one there.
///
/// **A `None` class is a coplanar contact** — the two solids meet on that plane — so nothing may
/// treat it as belonging to one side.
fn class_owners(plane_ix: &[ClassIx], n_a: usize, n_class: usize) -> Vec<Option<SolidSide>> {
    let mut out: Vec<Option<SolidSide>> = vec![None; n_class];
    let mut seen = vec![false; n_class];
    for (fi, &ci) in plane_ix.iter().enumerate() {
        // Cylinder classes get their own owner table with the cylinder class table (C2).
        let ClassIx::Plane(c) = ci else { continue };
        let side = if fi < n_a { SolidSide::A } else { SolidSide::B };
        if !seen[c] {
            seen[c] = true;
            out[c] = Some(side);
        } else if out[c] != Some(side) {
            out[c] = None;
        }
    }
    out
}

/// A face's class in the arrangement — **which index space the face's surface lives in** (M6-2a).
/// The old `plane_ix: Vec<usize>` presumed every class is a plane; the enum makes a cylinder
/// class unrepresentable as a plane index instead of smuggling it through a sentinel.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum ClassIx {
    /// An index into the dense plane-class table (`PlaneSetup::geom`).
    Plane(usize),
    /// An index into the cylinder-class table (arrives with C2; the index space exists first so
    /// the type is total).
    Cyl(usize),
}

impl ClassIx {
    /// The plane-class index. **Panics on a cylinder class** — same contract as
    /// [`FaceRow::plane`]: plane-only paths sit behind kind filters, and a cylinder here is an
    /// upstream filter bug.
    #[inline]
    #[track_caller]
    pub(crate) fn plane(self) -> usize {
        match self {
            ClassIx::Plane(i) => i,
            ClassIx::Cyl(k) => panic!("a plane-only path got cylinder class {k}"),
        }
    }

    /// The cylinder-class index, `None` on a plane — the curved paths' filter (they are the ones
    /// that *choose* per kind rather than assuming one).
    #[inline]
    pub(crate) fn cyl(self) -> Option<usize> {
        match self {
            ClassIx::Cyl(k) => Some(k),
            ClassIx::Plane(_) => None,
        }
    }
}

/// A cylinder class of one boolean (M6-2a): the exact statement the population gate reasons
/// about, beside its f64 cache. One entry per distinct lateral surface, in [`ClassIx::Cyl`]
/// numbering order.
pub(crate) struct WorkingCyl {
    #[allow(dead_code)] // the arrangement's circle elements read these from C3 on
    pub(crate) surf: Handle<Surface>,
    pub(crate) def: nacre_topo::CylinderDef,
    /// The f64 twin of `def` — what a *measurement* reads (`branch_vertex_tol` measures a branch
    /// realization against this surface), while every decision reads `def`.
    pub(crate) cache: nacre_geom::Cylinder,
}

/// **The M6-2a population gate** — decides, exactly, whether this operand pair stays inside
/// the axis-perpendicular population the cylinder arrangement serves, and names the refusal
/// otherwise. All arithmetic is checked `Rat` on world-stated descriptions; anything the gate
/// cannot decide exactly is [`RejectReason::CylinderGateUndecided`] — a conservative honest
/// refusal, never a guess.
///
/// Per (plane class, cylinder) pair, with `n` the class's rational normal and `m`/`o`/`r` the
/// cylinder's raw axis/origin/radius:
/// - `n × m = 0` — a perpendicular cut. Passes, a cap seated flush on the other body included.
/// - `n · m = 0` — a wall parallel to the axis. It must provably miss the **rectangle** the
///   cylinder occupies in that plane. The infinite plane clearing the axis by more than `r`
///   (`(n·o + d)² > r²·|n|²`) settles it outright and decides most inputs; otherwise each face on
///   the class answers for itself, across the strip or along a lateral face's span
///   ([`face_clears_footprint`]). A face not shown to miss is
///   [`RejectReason::WallMeetsLateral`] (M6-2b's rulings and arcs).
/// - anything else — [`RejectReason::ObliqueCylinderCut`] (an ellipse, M6-3).
///
/// Per cylinder pair: axes clear of each other (`dist > r₁+r₂`, whatever their orientation)
/// pass; a pair that touches or overlaps is [`RejectReason::CylinderPairContact`] (M6b, where
/// the quartic intersection curve lives).
pub(crate) fn cylinder_gate(
    model: &Model,
    cyl_surfs: &[Handle<Surface>],
    geom: &[WorkingPlane],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
) -> Result<Vec<WorkingCyl>, BoolError> {
    // ★ Every question below is a **sign**, and the scalar layer answers signs totally: the
    // local checked-`Rat` closures this used to carry declined on overflow, which put a width
    // limit inside `CylinderGateUndecided` and made that name say less than it claimed.
    use nacre_scalar::Orient;
    let undecided = || reject(RejectReason::CylinderGateUndecided);

    let mut cyls = Vec::with_capacity(cyl_surfs.len());
    for &surf in cyl_surfs {
        let nacre_topo::SurfaceTruth::Cylinder { def, motion } = model.surface_truth(surf) else {
            unreachable!("a cylinder row carries a cylinder truth")
        };
        // A moved cylinder's def is stated before its motion; the gate compares world
        // geometry, and realizing a def through its chain is machinery this population
        // does not have yet.
        if motion.is_some() {
            return Err(undecided());
        }
        let Surface::Cylinder(cache) = model.surface(surf) else {
            unreachable!("a cylinder truth carries a cylinder cache")
        };
        cyls.push(WorkingCyl {
            surf,
            def: def.clone(),
            cache: *cache,
        });
    }

    // ★ **One row per `cyl_surfs` entry, in order** — the loop above either pushes or returns, so
    // this table is indexed by the very `ClassIx::Cyl` number that named the surface. Consumers
    // index it directly (`merge_circles`), which is only sound while that holds.
    debug_assert_eq!(
        cyls.len(),
        cyl_surfs.len(),
        "the cylinder class table is index-aligned with the class numbering"
    );

    for cyl in cyls.iter() {
        let (o, m, r) = (cyl.def.origin(), cyl.def.dir(), cyl.def.radius());
        // The footprint's second axis, gathered **lazily and at most once** per cylinder: it does
        // not depend on the plane class the loop below walks, but almost no boolean ever asks for
        // it — the plane-level test decides first. ★ Measured before this was made lazy: one cut
        // over a plate with 16 bores built the table 16 times and read it 0, which is exactly the
        // shape `wall_faces_clear` warns about two doc comments below.
        let mut spans: Option<Vec<[nacre_scalar::Rat; 2]>> = None;
        for (c, wp) in geom.iter().enumerate() {
            if wp.rotated {
                return Err(undecided());
            }
            let Some(coeffs) = wp.base_rat else {
                return Err(undecided());
            };
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            // A perpendicular cut is the circle population, and it passes — **including a cap
            // seated flush on the other body's face**. That seating used to be refused, and the
            // refusal was wider than anything it could name: what makes a seated circle hard is
            // its boundary meeting the counterpart's boundary, and a boundary is either an edge
            // on a plane (that plane is parallel to the axis → the wall rule below, or oblique →
            // the oblique rule; both are already conservative because the wall rule judges the
            // *infinite* plane, not the face) or another cylinder's rim (two circles can only
            // overlap when the axes stand closer than r₁+r₂ → the pair rule below). So what the
            // seated rule turned away was exactly the population whose circle lies wholly inside
            // the counterpart's face — measured across the family (through hole, blind hole, boss
            // fuse, common, drilling a pocket floor): every one exact and validating clean.
            //
            // ★ **Nothing downstream has to catch a degenerate seating, because this gate still
            // does** — by the two rules below rather than by a rule about seating. A seated
            // circle can only reach the counterpart's boundary through a plane parallel to the
            // axis (which must prove clearance > r) or an oblique one (refused outright) or
            // another cylinder's rim (the pair rule), so a tangency or a crossing is named
            // before the arrangement ever sees it.
            if !nacre_scalar::parallel_rat(&n, &m) {
                if nacre_scalar::dot_sign_rat(&n, &m) != Orient::Zero {
                    return Err(reject(RejectReason::ObliqueCylinderCut));
                }
                // A parallel wall must provably miss the lateral surface. ★ **The question is
                // about the wall's *faces*, not its plane** — the uniform-slab theorem this feeds
                // says so in its own words ("the other operand's **boundary** does not meet the
                // open cylinder slab"). Judging the infinite plane is a cheaper *sufficient*
                // condition, so it is asked first and still decides most inputs; when it fails,
                // the faces on this class get to answer for themselves. Refusing on the plane
                // alone turned away a whole family the engine serves — a boss standing far away
                // whose wall plane, extended, happens to pass through a hole.
                //
                // ★★ What a face is asked is whether it misses the **rectangle** this cylinder
                // occupies in that plane: the strip across, the lateral face's span along. See
                // [`face_clears_footprint`].
                if nacre_scalar::point_plane_clearance_rat(&coeffs, &o, r) != Orient::Positive {
                    // ★ The class's coefficients come from its root's *name*, which is the world
                    // only while that plane carries no motion; the cylinder's `def` is world by
                    // the check above. Comparing across those two frames would be a silent wrong
                    // answer, so a moved class is refused rather than measured.
                    if model.plane_motion(wp.surf).is_some() {
                        return Err(undecided());
                    }
                    let spans = spans.get_or_insert_with(|| lateral_spans(faces, cyl.surf));
                    if !wall_faces_clear(model, faces, plane_ix, c, &coeffs, &o, &m, r, spans)? {
                        return Err(reject(RejectReason::WallMeetsLateral));
                    }
                }
            }
        }
    }

    for (i, a) in cyls.iter().enumerate() {
        for b in &cyls[i + 1..] {
            // Clear iff the distance between the two axes exceeds the radius **sum** — one
            // proposition for any pair, parallel or not. ★ This used to demand parallel axes
            // first and refuse everything else, which made a limit of the *arithmetic* read as
            // a limit of the kernel: a drill crossing a bore at a safe distance came back as
            // "two cylinders touch". The predicate now carries both spellings of the distance,
            // so the refusal below means what it says.
            if nacre_scalar::cylinders_clear(
                &a.def.origin(),
                &a.def.dir(),
                a.def.radius(),
                &b.def.origin(),
                &b.def.dir(),
                b.def.radius(),
            ) != Orient::Positive
            {
                return Err(reject(RejectReason::CylinderPairContact));
            }
        }
    }
    Ok(cyls)
}

/// **Does every face on plane class `c` provably miss this cylinder?** — the boundary question
/// [`cylinder_gate`]'s wall rule asks once the cheaper plane-level one has failed.
///
/// `Err` is the honest "could not decide exactly"; `Ok(false)` means some face was not shown to
/// clear, which is the wall refusal's whole content.
///
/// ★ **No table is built for this.** The scan runs only on the class that failed the plane test —
/// rare — so walking the face rows there costs nothing on the common path and allocates nothing.
/// A per-class face table computed for every boolean would be the shape M6-2's `caps` had, built
/// for everyone and read by almost no one. ★★ The `spans` this takes is held to the same rule: the
/// caller builds it on first use, not per cylinder — measured, an ordinary cut over a 16-bore
/// plate wants it **zero** times.
#[allow(clippy::too_many_arguments)]
fn wall_faces_clear(
    model: &Model,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    c: usize,
    coeffs: &[nacre_scalar::Rat; 4],
    o: &[nacre_scalar::Rat; 3],
    m: &[nacre_scalar::Rat; 3],
    r: nacre_scalar::Rat,
    spans: &[[nacre_scalar::Rat; 2]],
) -> Result<bool, BoolError> {
    let mut seen = 0usize;
    for (i, row) in faces.iter().enumerate() {
        let (ClassIx::Plane(k), FaceRow::Plane(fi)) = (plane_ix[i], row) else {
            continue;
        };
        if k != c {
            continue;
        }
        let Some(fh) = fi.face else {
            return Err(reject(RejectReason::CylinderGateUndecided));
        };
        seen += 1;
        if !face_clears_footprint(model, model.faces.get(fh), coeffs, o, m, r, spans)? {
            return Ok(false);
        }
    }
    // A class exists because faces made it, so finding none is a wiring failure rather than an
    // input — and "every one of no faces clears" is a pass this must not hand out by default.
    if seen == 0 {
        return Err(reject(RejectReason::CylinderGateUndecided));
    }
    Ok(true)
}

/// The axis-parameter spans this cylinder's lateral faces occupy — one per face, in the scale
/// `axis_param_of_plane` produces.
///
/// ★★ **Empty means "unusable", never "nothing in the way".** A span this returns is a rectangle
/// the wall rule must miss, and "every one of no rectangles is missed" is a pass no caller should
/// hand out by default — so a face whose span could not be stated, and a cylinder with no lateral
/// face in the table at all, both come back empty and leave the axis unused. The second could be
/// argued safe (no lateral face, no band, nothing to protect), but that argument rests on the face
/// table being complete here, which is a separate premise from the one this function is about.
fn lateral_spans(faces: &[FaceRow], surf: Handle<Surface>) -> Vec<[nacre_scalar::Rat; 2]> {
    let mut out = Vec::new();
    for row in faces {
        let FaceRow::Cylinder(cf) = row else { continue };
        // ★ Matched by **handle**, not by geometry: two operands may state the same cylinder
        // twice, and each statement gets its own rows, its own bands clipped to its own spans,
        // and its own turn through the gate's cylinder loop.
        if cf.surf != surf {
            continue;
        }
        let Some(span) = cf.span else {
            return Vec::new();
        };
        // ★ The reader below asks "every vertex at or below `span[0]`, or every one at or above
        // `span[1]`", which is only the interval's outside while the pair is ordered — reversed,
        // that same phrasing reads "outside the *union* of two half-lines" and would clear a face
        // sitting squarely in the band. `lateral_axis_span` orders it; this is where that is
        // relied on, so this is where it is said.
        debug_assert!(span[0] <= span[1], "a span is stated low end first");
        out.push(span);
    }
    out
}

/// Whether one planar face misses the **rectangle** this cylinder occupies in the face's plane.
///
/// ★★ **One question, two separating axes.** A plane parallel to the axis meets the solid cylinder
/// in a rectangle: the strip across ([`nacre_scalar::cylinder_strip_side`]) and a lateral face's
/// axis-parameter span along ([`nacre_scalar::point_axis_side`]). A rectangle is the intersection
/// of those two bands, so clearing *either* axis clears it — these are not two rules to be
/// weighed, they are the two axes of one. With several lateral faces there are several rectangles
/// (the strip is shared, the spans are not), and the face must miss them all:
///
/// ```text
/// clear  ⟺  clear across the strip  ∨  clear along **every** span
/// ```
///
/// ★★★ **"Every" over an empty list is a pass, and that one is not on offer.** With no usable
/// span the axis along is skipped outright rather than answered vacuously. Not a hypothetical:
/// measured with the emptiness guard removed, a wall that genuinely crosses the bore and one
/// tangent to it were both waved straight through.
///
/// ★ **The span reading is an open interval** — a face resting exactly on a cap plane is clear,
/// because the theorem being fed speaks of the *open* slab. What such a face touches is the rim's
/// own plane, and whether its edge crosses the rim circle there is a question the arrangement asks
/// where the circles and segments are, by name.
///
/// ★ **The two axes are only *complete* for a face that is an axis-aligned rectangle** — then both
/// rectangles' edge normals coincide and there are no other separating axes to try. An extruded
/// wall is exactly that shape (its edges run along the axis or across it), chamfers included: what
/// a chamfer tilts is the plane, not the edges within it. A face a previous boolean took a bite
/// out of is not, and neither is a slanted or L-shaped one; those are refused as "not shown to
/// clear", which is true. Completing the test means adding the face's own edge normals as further
/// axes — the same test with more axes, not a different machine — and waits for a shape that
/// needs it.
///
/// **Only the outer loop is walked**, and that is sound: a face is contained in the convex hull of
/// its outer loop's vertices, a half-space is convex, and inner loops only *remove* material. It
/// is also what keeps the test reachable — a wall drilled by a crosswise bore carries that bore's
/// rim as an inner loop, whose vertices have no rational meet at all.
///
/// ★★ **The straight-edge demand is soundness, not convenience.** This plane may be some *other*
/// cylinder's cap plane (one whose axis is perpendicular to it), and such a cap face's outer loop
/// is a single circle with one seam vertex — "all vertices on one side" would then be satisfied by
/// a single point and would pass a disk that crosses the strip. Today `vertex_meet` happens to
/// decline that vertex, but relying on the coincidence would leave the barrier to vanish silently
/// the day seam coordinates become solvable.
#[allow(clippy::too_many_arguments)]
fn face_clears_footprint(
    model: &Model,
    face: &Face,
    coeffs: &[nacre_scalar::Rat; 4],
    o: &[nacre_scalar::Rat; 3],
    m: &[nacre_scalar::Rat; 3],
    r: nacre_scalar::Rat,
    spans: &[[nacre_scalar::Rat; 2]],
) -> Result<bool, BoolError> {
    use nacre_scalar::{Orient, StripSide};
    let undecided = || reject(RejectReason::CylinderGateUndecided);
    let mut side: Option<StripSide> = None;
    let mut across = true;
    // Per span, whether every vertex so far has stayed at or below its start, and at or above its
    // end. Either one surviving the walk clears that rectangle along the axis.
    let mut along: Vec<(bool, bool)> = vec![(true, true); spans.len()];
    let mut vertices = 0usize;
    for he in &face.outer.half_edges {
        if !matches!(model.edge_curve(he.edge), nacre_geom::Curve::Line(_)) {
            return Err(undecided());
        }
        let Some((p, frame)) = model.vertex_meet(he_start(model, *he)) else {
            return Err(undecided());
        };
        // The point is stated in `frame`; the cylinder and the coefficients are world.
        if frame.is_some() {
            return Err(undecided());
        }
        // ★ **The class's coefficients must actually describe *this* face's plane.** Classes merge
        // on three exact witnesses, one of which compares *rounded* coefficients — so a face can
        // sit in a class whose exact name its own vertices do not satisfy (the two-descriptions
        // hazard `FaceInfo::exact_coeffs` documents). The strip decomposition takes the plane's
        // distance from the axis as the point's, so judging a point against a plane it is not on
        // would answer about geometry that is not there. Checked, not assumed: a `debug_assert`
        // would say nothing in the build that ships.
        if !nacre_scalar::point_on_plane_exact(coeffs, &p) {
            return Err(undecided());
        }
        vertices += 1;
        // The axis across the strip.
        match nacre_scalar::cylinder_strip_side(coeffs, &p, o, m, r) {
            StripSide::Inside => across = false,
            s => match side {
                None => side = Some(s),
                Some(prev) if prev != s => across = false, // the face straddles the strip
                Some(_) => {}
            },
        }
        // The axis along it — one rectangle per lateral face, and the span is open at both ends.
        for (i, span) in spans.iter().enumerate() {
            if nacre_scalar::point_axis_side(&p, o, m, span[0]) == Orient::Positive {
                along[i].0 = false;
            }
            if nacre_scalar::point_axis_side(&p, o, m, span[1]) == Orient::Negative {
                along[i].1 = false;
            }
        }
    }
    // ★ An empty outer loop names no half-space and no interval, so it proves nothing — and
    // "every vertex of none stayed below" would otherwise be a vacuous pass for any wall.
    if vertices == 0 {
        return Ok(false);
    }
    if across && side.is_some() {
        return Ok(true);
    }
    // `spans` empty means the axis is unusable, not that every rectangle was missed.
    Ok(!spans.is_empty() && along.iter().all(|(below, above)| *below || *above))
}

/// The minimal per-op plane table two solids share: the
/// concatenated plane list (`a`'s then `b`'s), the face→index map, and each solid's
/// [`combinatorics::EdgeFaces`]. Built once and shared: indices into the returned `planes`/`surf_ix`
/// are common to both solids, so a vertex of `a` and a face of `b` compose in one index space.
/// Destructure it with `..` (`let PlaneSetup { planes: faces_tab, geom: planes, plane_ix, .. } = …`):
/// the tables here grow as the arrangement learns to say "plane" and "face" in different index
/// spaces, and a positional tuple made every one of those steps touch all ~25 call sites.
///
/// The plane classes (`canon`) are computed here to build `geom`/`plane_ix` and then dropped — the
/// dense `plane_ix` is the only face→plane map anything downstream needs, so the sparse union-find
/// output does not escape.
pub(crate) struct PlaneSetup {
    pub(crate) planes: Vec<FaceRow>,
    pub(crate) surf_ix: HashMap<Handle<Face>, usize>,
    pub(crate) inc_a: combinatorics::EdgeFaces,
    pub(crate) inc_b: combinatorics::EdgeFaces,
    /// Where `a`'s faces end and `b`'s begin in `planes`. The concatenation always created this
    /// boundary; it was just never written down, so every later "whose face is this?" had to
    /// rebuild it.
    pub(crate) n_a: usize,
    /// The arrangement's planes, densely indexed — see [`dense_planes`].
    pub(crate) geom: Vec<WorkingPlane>,
    /// `plane_ix[face]` is that face's class — a plane index into `geom`, or a cylinder class
    /// ([`ClassIx`]).
    pub(crate) plane_ix: Vec<ClassIx>,
    /// Whose faces each plane class carries — see [`class_owners`].
    pub(crate) class_owner: Vec<Option<SolidSide>>,
    /// The cylinder classes, in [`ClassIx::Cyl`] numbering order — empty for an all-planar
    /// boolean. Filled by the population gate, which is also what refuses the interactions this
    /// milestone does not build.
    pub(crate) cyls: Vec<WorkingCyl>,
    /// How this operation judges, and where its evidence goes — the two facts that belong to the
    /// operation rather than to any one plane. The caller pairs them with a table to make a
    /// [`Judge`].
    pub(crate) standard: Standard,
    pub(crate) notes: Notes,
}

/// Which sub-phase of [`plane_index_setup`] a [`Watch`] charges — spike instrumentation, and only
/// in a test build (see `arrangement::phase`).
pub(crate) enum Sub {
    TriPt3,
    Std,
    Collect,
    Edges,
    Classes,
    Dense,
}

/// Times the scope it is charged from, or does nothing at all in a release build.
pub(crate) struct Watch(#[cfg(test)] std::time::Instant);

impl Watch {
    pub(crate) fn new() -> Self {
        Watch(
            #[cfg(test)]
            std::time::Instant::now(),
        )
    }
    #[allow(unused_variables)]
    pub(crate) fn charge(self, which: Sub) {
        #[cfg(test)]
        {
            use crate::arrangement::phase;
            let c = match which {
                Sub::TriPt3 => &phase::S_TRIPT3,
                Sub::Std => &phase::S_STD,
                Sub::Collect => &phase::S_COLLECT,
                Sub::Edges => &phase::S_EDGES,
                Sub::Classes => &phase::S_CLASSES,
                Sub::Dense => &phase::S_DENSE,
            };
            c.fetch_add(
                self.0.elapsed().as_nanos() as u64,
                std::sync::atomic::Ordering::Relaxed,
            );
        }
    }
}

pub(crate) fn plane_index_setup(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<PlaneSetup, BoolError> {
    let (mut setup, cyl_surfs) = plane_index_setup_inner(model, a, b)?;
    // ★ The cylinder door: the population gate decides **by name** what stands in the way (an
    // oblique cut, a wall touching the lateral surface, an undecidable pair), and what passes now
    // goes on to be arranged. The `CylinderBooleanNotYet` stopper that stood here from C2 to
    // C4b-2 is gone — the bands and the assembly that serve this population landed.
    if !cyl_surfs.is_empty() {
        let cyls = cylinder_gate(
            model,
            &cyl_surfs,
            &setup.geom,
            &setup.planes,
            &setup.plane_ix,
        )?;
        setup.cyls = cyls;
    }
    Ok(setup)
}

fn plane_index_setup_inner(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(PlaneSetup, Vec<Handle<Surface>>), BoolError> {
    let t = Watch::new();
    let mut planes = collect_planes(model, a)?;
    let n_a = planes.len();
    planes.extend(collect_planes(model, b)?);
    t.charge(Sub::Collect);
    let t = Watch::new();
    let standard = standard_for(&planes);
    t.charge(Sub::Std);
    let notes = Notes::new();
    if standard.prec > JUDGE_PREC_CAP {
        return Err(reject(RejectReason::PrecisionBudget {
            needed: standard.prec,
            cap: JUDGE_PREC_CAP,
        }));
    }
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        // Synthetic faces are appended later, after this table is built; every entry here is real.
        surf_ix.insert(pi.face().expect("collect_planes yields real faces"), i);
    }
    let t = Watch::new();
    let inc_a = combinatorics::edge_faces(model, a, &surf_ix)?;
    let inc_b = combinatorics::edge_faces(model, b, &surf_ix)?;
    t.charge(Sub::Edges);
    // One judging context for the whole operation: the witnesses, the standard they are held to,
    // and where the evidence goes. The face table judges first (it is what *defines* the plane
    // classes), then the dense plane table inherits the same three.
    let t = Watch::new();
    let canon = plane_classes(&Judge::new(&planes, standard, &notes));
    t.charge(Sub::Classes);
    let t = Watch::new();
    let (geom, plane_ix, cyl_surfs) = dense_planes(&planes, &canon);
    let class_owner = class_owners(&plane_ix, n_a, geom.len());
    t.charge(Sub::Dense);
    Ok((
        PlaneSetup {
            planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            geom,
            plane_ix,
            class_owner,
            cyls: Vec::new(),
            standard,
            notes,
        },
        cyl_surfs,
    ))
}

/// **How precisely this operation's rotated definitions must be realized.**
///
/// The judges' error radius is `C · 2⁻ᵖʳᵉᶜ`, and `C` belongs to the model — it grows about one
/// bit per turn of rotation history and with the coordinate magnitudes. A fixed precision
/// therefore decides, silently, how long a model's history may be: at 256 bits a solid turned 245
/// times stops building, with a reject that names a symptom rather than the cause. So the
/// precision is read off the model instead.
///
/// The target is the **coincidence precision**: two things closer than this are treated as
/// coincident, and the kernel will only say so once it has *proved* the separation is below it.
/// Its default is derived rather than chosen —
///
/// - `output_precision = scale · 2⁻⁵²`, the finest distinction the `f64` coordinates this kernel
///   emits can carry. Below it nothing survives export, so distinguishing is meaningless.
/// - `coincidence_precision = output_precision · 2⁻¹²⁸`, two whole words further down. Erring low
///   only costs bits, while erring high merges features that were genuinely apart, so the
///   asymmetry says push it down; and a word is the natural unit because astro-float allocates
///   whole words anyway.
///
/// `scale` is the largest coordinate magnitude in either operand, taken over the whole table so
/// the result does not depend on traversal order (replay must reproduce it exactly).
///
/// The precision that reaches the target is then [`nacre_cip::judge_precision`]'s to compute;
/// [`JUDGE_PREC_CAP`] is where the kernel stops and says so instead, and [`CLIMB_HEADROOM`] is
/// what a single hard judgement may spend on top of it.
///
/// **Only the coincidence limit is a candidate for a setting.** Everything else here — the
/// precision, the cap, the headroom — is derived from it and from the model, because a bit count
/// means a different physical thing in every model ("256 bits" is `1e-76` for a solid turned once
/// and `1e+15` for one turned three hundred times).
/// The plane data of a row, `None` for a cylinder — the kind filter the plane-only sweeps
/// share.
#[inline]
pub(crate) fn plane_of(r: &FaceRow) -> Option<&FaceInfo> {
    match r {
        FaceRow::Plane(p) => Some(p),
        FaceRow::Cylinder(_) => None,
    }
}

// A cylinder row contributes no witness points to the standard — structurally right, not an
// omission: an axis-aligned-grade rational cylinder is the `rotated == false` case (its exact
// def needs no high-precision realization); the rotated-cylinder story is M6-3's (CIP).
fn standard_for(rows: &[FaceRow]) -> Standard {
    // ★ **A face that was never moved contributes exactly nothing, so it is not asked.**
    //
    // Its `tri_pt3` are `WitnessPoint::exact` of the face's own f64 triangle: base = `mantissa · 2^exp`, so
    // the denominator is a power of two and the numerator fits `Rat`'s 127 bits, and the chain is
    // empty — which is precisely when `rat_to_hp` returns an *exact* interval. The realization has
    // no error to report (`an_exact_point_demands_no_precision` in `nacre-cip`, and the const
    // assert at `TRIAL_PREC` that keeps it true).
    //
    // So the loop below used to spend a full high-precision replay per point to compute a zero —
    // measured, an axis-aligned 60-fin fold did that 24,120 times for 15.7ms and a `worst` of
    // exactly `Mag::ZERO`. The same shape was removed one level down when `WitnessPoint::exact` replaced
    // `WitnessPoint::at` for these points ("nine BigFloat operations to compute a zero").
    //
    // `max` over the empty set is `Mag::ZERO`, which is the right answer for a model with no
    // rotation history — `precision_for` reads that as "nothing to size" and returns `TRIAL_PREC`.
    let worst = worst_trial(
        rows.iter()
            .filter_map(plane_of)
            .filter(|p| p.rotated)
            .flat_map(|p| p.tri_pt3.iter()),
    );
    // `scale`, by contrast, is every point's business: it is the model's size, and an unmoved face
    // is as far from the origin as any other.
    standard_from(
        rows.iter()
            .filter_map(plane_of)
            .flat_map(|p| p.tri_pt3.iter()),
        worst,
    )
}

/// **How deep a model may be before the operation is rejected instead.**
///
/// Not a resolution limit — the arithmetic is correct at any depth — but a **cost** limit, so it
/// is set from measured cost. A judgement's realization is quadratic-ish in the precision, and the
/// cap is placed where a single boolean's judging stays in the seconds rather than the minutes:
/// 4096 bits covers a rotation history of roughly four thousand turns (measured: `C` grows one bit
/// per turn), which is far past any real model, and a model that does exceed it is told *why*
/// rather than handed a wrong answer or an unbounded wait.
pub(crate) const JUDGE_PREC_CAP: usize = 4096;

/// **How thin a witness the kernel will still judge**, expressed as the bits a single judgement
/// may ask for *beyond* what the model itself needed.
///
/// This is a **separate budget from [`JUDGE_PREC_CAP`], and it has to be.** Sharing one absolute
/// ceiling would mean a deeply-turned model — already near the cap — leaves a hard judgement no
/// room at all, so the same sliver would be judged in a fresh model and abandoned in a turned one.
/// The model's depth and a judgement's difficulty are different quantities; only the second
/// belongs here.
///
/// It has a physical reading. A judgement's uncertainty is `(C / |cofactor|) · 2⁻ᵖʳᵉᶜ`, and the
/// model already chose `prec` so that `C · 2⁻ᵖʳᵉᶜ` clears the coincidence limit; what is left is
/// `log₂(1 / cofactor)` — the **thinness of the witness**, a needle triangle or three planes that
/// almost share a line. Two words says: a witness up to `2¹²⁸` (≈ 3·10³⁸) times more degenerate
/// than the model's own size is still judged to the end.
///
/// Two words, and not a measured number, for the same reason the coincidence limit is two words
/// below the output resolution: the error is asymmetric. Too small abandons a judgement that had
/// an answer; too large only spends bits. And measurement says there is nothing to tune — across
/// the rotation corpus and models turned 100 and 800 times, **no judgement asked for even one bit
/// beyond the model's own precision** (measured with the headroom forced to zero).
pub(crate) const CLIMB_HEADROOM: usize = 128;

/// A judging context over a hand-built table, for fixtures.
///
/// The standard is the derived default for a unit-scale model, and the collector is leaked so a
/// fixture is a one-liner — a handful of `Vec`s per test run, and nothing reads them. A fixture
/// that *does* want the evidence builds its own [`Notes`] and calls [`Judge::new`].
#[cfg(test)]
pub(crate) fn test_judge<W>(planes: &[W]) -> Judge<'_, W> {
    let notes: &'static Notes = Box::leak(Box::new(Notes::new()));
    Judge::new(
        planes,
        Standard {
            prec: 256,
            coincidence: Mag::pow2(-180),
            scale: Mag::of(1.0),
            cap: 256 + CLIMB_HEADROOM,
        },
        notes,
    )
}

/// [`standard_for`] over a bare set of definitions, with no plane table — **a test helper.**
///
/// It once served witness selection in `rotated_vertex`; that consumer is gone, and splitting
/// [`worst_trial`] out of [`standard_from`] is what surfaced it. Kept because a fixture that asks
/// "what precision does *this* point demand" wants exactly the two halves in order, and spelling
/// them out at every call site says less than the name does.
#[cfg(test)]
pub(crate) fn standard_for_points<'a>(
    pts: impl IntoIterator<Item = &'a WitnessPoint> + Clone,
) -> Standard {
    standard_from(pts.clone(), worst_trial(pts))
}

/// **The realization depth this set of definitions demands** — `max` over their trial bounds.
///
/// Split from [`standard_from`] because the two halves have nothing in common but the answer: this
/// one is **all of the cost** (a full high-precision replay per point), and the other is f64
/// arithmetic on already-known numbers. Keeping them apart is what lets a caller that already knows
/// this maximum skip straight to the second half.
///
/// **This is where a boolean spends most of what is left after the arrangement went parallel**
/// (measured: 76% of setup, and setup is 43% of the largest booleans once the trace is off the
/// critical path). Each point's trial realization is independent and they combine by **maximum**,
/// which is associative and exact — so evaluating them across cores cannot move the answer the way
/// a reassociated sum would.
fn worst_trial<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>) -> Mag {
    let pts: Vec<&WitnessPoint> = pts.into_iter().collect();
    let bounds = crate::par::map_range(pts.len(), |i| nacre_cip::trial_bound(pts[i]));
    bounds
        .into_iter()
        .fold(Mag::ZERO, |w, b| if w.lt(b) { b } else { w })
}

/// The standard for points whose worst trial bound is already known: `scale` off the f64
/// coordinates, the coincidence limit derived from it, and the precision that reaches it.
fn standard_from<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>, worst: Mag) -> Standard {
    let mut scale = 1.0f64;
    for p in pts {
        for c in p.coord {
            scale = scale.max(c.abs());
        }
    }
    let scale = Mag::of(scale);
    let output_precision = scale.times(Mag::pow2(-52));
    let coincidence = output_precision.times(Mag::pow2(-128));
    let prec = nacre_cip::precision_for(worst, coincidence);
    Standard {
        prec,
        coincidence,
        scale,
        // Relative to this model's own depth — see [`CLIMB_HEADROOM`]. The absolute ceiling that
        // leaves is `JUDGE_PREC_CAP + CLIMB_HEADROOM`, since a model deeper than the first is
        // rejected before any judging starts.
        cap: prec + CLIMB_HEADROOM,
    }
}

/// One plane of the arrangement, indexed by a **dense** class id.
///
/// The face table cannot answer "which plane" without a convention: a class holds faces from both
/// operands, and two of them can face opposite ways, so there is no such thing as *the* plane's
/// outward normal. What a plane has is a **frame** — the class root's stored normal — and the only
/// direction fact anyone needs from it is [`WorkingPlane::frame_sign`]. Everything else here is a
/// witness: three points known to lie on this plane, used to reconstruct it exactly.
/// The pre-rotation twin of a witness triangle: its `chain_id`, base points and base plane.
///
/// A rigid motion preserves the determinants the predicates take, so a judgement whose inputs all
/// carry **one** motion can be answered on these instead — exactly, off the toleranced path
/// entirely. `None` for the base data when a base coordinate is not `f64`-representable, since the
/// exact predicate takes `f64`; the judgement then stays toleranced (slower, never wrong).
#[derive(Clone, Copy, Debug)]
pub(crate) struct BaseFrame {
    /// `0` = no motion. Equal only for structurally identical chains.
    pub(crate) chain_id: u64,
    pub(crate) tri: Option<[Point3; 3]>,
    pub(crate) coeffs: Option<[f64; 4]>,
}

impl BaseFrame {
    /// No motion to cancel — for a plane with no history: a hand-built table in a test, or the
    /// synthetic split plane a subdivided boolean cuts with. Identical to what `of` returns for an
    /// unmoved face, so such a plane takes the same predicate routes an axis-aligned model does.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Self {
            chain_id: 0,
            tri: None,
            coeffs: None,
        }
    }

    fn of(
        tri_pt3: &[WitnessPoint; 3],
        motion: Option<Handle<nacre_topo::MotionNode>>,
        frame_sign: i8,
        base_rat: Option<[nacre_scalar::Rat; 4]>,
    ) -> Self {
        // **Identity by handle, not by hash.** This used to fold the chain into a 64-bit
        // `DefaultHasher` digest and compare digests — and a collision does not make a judgement
        // slow, it makes `shared_base` hand the *exact* predicate two incompatible pre-motion
        // frames and answer a different question with full confidence. The motion-history leaf is
        // the canonical name of "which motion": equal handles are the same chain by construction,
        // and two structurally-equal chains under different handles are a conservative miss.
        //
        // The three witnesses share one chain by construction — `collect_planes` builds all three
        // from the same surface truth — so there is nothing to cross-check here either.
        let Some(leaf) = motion else {
            return Self {
                chain_id: 0,
                tri: None,
                coeffs: None,
            }; // no motion to cancel
        };
        // 0 is reserved for "no motion", so never hand it out as an id.
        let chain_id = leaf.index() as u64 + 1;
        let exact = |r: nacre_scalar::Rat| nacre_scalar::Rat::try_from_f64(r.to_f64()) == Some(r);
        if !tri_pt3.iter().all(|p| p.base.iter().all(|&r| exact(r))) {
            return Self {
                chain_id,
                tri: None,
                coeffs: None,
            };
        }
        let pt = |p: &WitnessPoint| {
            Point3::from_array([p.base[0].to_f64(), p.base[1].to_f64(), p.base[2].to_f64()])
        };
        let tri = [pt(&tri_pt3[0]), pt(&tri_pt3[1]), pt(&tri_pt3[2])];
        // ★ **An improper chain is corrected here, on the points, before anything is derived
        // from them.** An odd number of reflections leaves the base frame with the opposite
        // handedness, and the shortcuts' whole licence is that the motion preserves the
        // determinants they take. One more reflection puts the handedness back, and then the base
        // is related to the moved frame by a *proper* motion again — which is what every consumer
        // of this struct assumes. A sign flip on x is exact for every finite `f64`, and it is the
        // same convention `nacre_cip::frame3::shared_base` applies to its own points.
        let improper = nacre_cip::chain_parity(&tri_pt3[0].chain) < 0;
        let tri = if improper {
            tri.map(|p| {
                let [x, y, z] = p.as_array();
                Point3::from_array([-x, y, z])
            })
        } else {
            tri
        };
        // ★ The base plane must carry the **stored** orientation, not the triangle's. A class's
        // stored normal and its witness triangle's `cross` can oppose — that is exactly what
        // `frame_sign` records — and `through_points` gives the triangle's. A proper motion
        // preserves the cross product (`det = 1`), so multiplying by `frame_sign` reproduces the
        // same relation in the base frame. Without it the exact path answers with a flipped sign,
        // which the suite caught immediately.
        //
        // ★ **And it is derived from the corrected triangle, not corrected afterwards.** A plane
        // is not a bag of points: reflecting a triangle and re-deriving its normal is *not* the
        // same as reflecting the normal, because the cross product is a pseudovector — the two
        // differ by a global sign, and `plane_pair_dir_sign` reads exactly that sign. Deriving
        // last removes the question: the points are corrected once, and everything downstream is
        // the ordinary derivation from them. (The earlier spelling corrected the plane separately
        // and was off by that one sign; `a_reflected_spelling_takes_the_same_direction_signs`
        // is what found it.)
        let derived = Plane::through_points(tri[0], tri[1], tri[2]).map(|pl| {
            let c = pl.coefficients();
            let k = f64::from(frame_sign);
            [c[0] * k, c[1] * k, c[2] * k, c[3] * k]
        });
        // ★★★ **Take `d` from the record and the direction from the triangle.**
        //
        // The derivation above is the two-descriptions problem in miniature: `d` comes out of an
        // f64 dot product, so the plane it names is not quite the one `tri` lies on — measured, for
        // 27% of the census's rotated classes and 40% of the fin sweep's. The surface's recorded
        // pre-motion coefficients (`Model::surface_name`) *are* that plane, exactly, with no
        // triangle in the derivation at all.
        //
        // ★ Only the **direction** still comes from the triangle, and that is deliberate. The
        // record is canonicalized, so its sign is a normal form, not this face's outward sense;
        // and the reflection correction above cannot simply be applied to a normal, because the
        // cross product is a pseudovector and reflecting-then-deriving differs from
        // deriving-then-reflecting by a global sign (the comment above, and the test
        // `a_reflected_spelling_takes_the_same_direction_signs` that found it). Orienting the
        // exact plane to agree with the derived one reproduces whatever convention the derivation
        // had, without re-deriving the convention — and *direction* is the half where the two
        // descriptions do not part.
        let exact_coeffs = base_rat.and_then(|c| {
            // ★ **The same correction, applied to the plane.** When the chain is improper the
            // triangle above was reflected in `x`, so everything derived from it lives in the
            // reflected base frame — and the record does not. Orienting the normals afterwards
            // cannot repair that: an unreflected plane and a reflected one are *different planes*,
            // not the same plane spelled with the opposite sign, so the two mirror fixtures fail
            // outright. Reflect the plane, then let the orientation step below settle the sign
            // (which is where reflecting-then-deriving and deriving-then-reflecting differ).
            let c = if improper {
                nacre_scalar::mirror_plane_coeffs(
                    c,
                    nacre_scalar::Axis::X,
                    nacre_scalar::Rat::from_int(0),
                )?
            } else {
                c
            };
            let f = c.map(|r| r.to_f64());
            // A canonicalized vector is integral; if it does not survive the round trip the
            // realization is a rounding and buys nothing over the derivation.
            c.iter()
                .zip(f)
                .all(|(&r, x)| nacre_scalar::Rat::try_from_f64(x) == Some(r))
                .then_some(f)
        });
        let coeffs = match (exact_coeffs, derived) {
            (Some(e), Some(d)) => {
                let dot = e[0] * d[0] + e[1] * d[1] + e[2] * d[2];
                Some(if dot < 0.0 { e.map(|x| -x) } else { e })
            }
            _ => derived,
        };
        Self {
            chain_id,
            tri: Some(tri),
            coeffs,
        }
    }
}

#[derive(Clone)]
pub(crate) struct WorkingPlane {
    pub(crate) plane: Plane,
    /// The class root's exact rational coefficients — see [`FaceInfo::base_rat`].
    pub(crate) base_rat: Option<[nacre_scalar::Rat; 4]>,
    /// The class's representative surface — what `assemble_fuse_cut` records in a
    /// `VertexDef::ThreePlane`.
    pub(crate) surf: Handle<Surface>,
    /// Witness points on this plane (the root face's `tri`), outward-ordered for that face.
    pub(crate) tri: [Point3; 3],
    /// The witness as exact `WitnessPoint` definitions (the root face's), borrowed by every predicate.
    pub(crate) tri_pt3: [WitnessPoint; 3],
    /// Whether the root face's solid is rotated — copied from it together with `tri_pt3` so the
    /// pair cannot disagree. See [`FaceInfo::rotated`].
    pub(crate) rotated: bool,
    /// `+1` when the plane's stored normal agrees with the root face's outward normal, `-1` when
    /// they oppose. This *is* the label frame: `[A_above, A_below, …]` is defined about the class
    /// root's stored normal, and this sign is what relates it to material. Precomputed here so the
    /// two `debug_assert`s that guard the convention run once, at construction.
    pub(crate) frame_sign: i8,
    /// The pre-rotation twin — see [`BaseFrame`].
    pub(crate) base: BaseFrame,
    /// ★★★ **The coefficients, but only where they describe the same plane as [`tri`](Self::tri).**
    ///
    /// A plane has two exact descriptions and they need not agree: `d` is the `f64` product
    /// `raw·origin`, so a face at `y = −0.2` gets a coefficient plane `2⁻⁵⁴` from the one its own
    /// witness spans (`Plane::coefficients` has the numbers). A predicate that describes one plane
    /// by its coefficients in one question and by its triangle in the next composes answers about
    /// **two different planes**, and what comes out is not even an order.
    ///
    /// So the disagreement is resolved here rather than guarded against at every call: when the
    /// two do not agree, this is `None` and there is nothing to describe the plane with except its
    /// triangle. **The same shape [`BaseFrame`] already uses** — a description that cannot be
    /// trusted is not carried, so no consumer has to remember to check it.
    ///
    /// `None` for a rotated plane too: there are no exact `f64` coefficients for one.
    pub(crate) exact_coeffs: Option<[f64; 4]>,
    /// The **normal** under the weaker agreement — parallel to what `tri` spans, direction not
    /// required (`frame_sign` records that separately).
    ///
    /// ★ `d` is where the two descriptions part, so a predicate that never reads it can keep its
    /// exact route on a plane [`exact_coeffs`](Self::exact_coeffs) has to refuse. Measured:
    /// demanding the full agreement for those cost 4.7x on the axis-aligned fold and bought
    /// nothing.
    pub(crate) exact_normal: Option<[f64; 3]>,
    /// The class root's name integers, folded to the stored orientation
    /// ([`nacre_cip::predicate::name_stored_ints`]) — what gives a **wide** name its exact
    /// judging shortcuts back (truth-and-cache open item 15). `None` when the root's surface
    /// has no name.
    pub(crate) name_ints: Option<nacre_cip::predicate::NameInts>,
}

impl WorkingPlane {
    /// **Reconcile a plane's two exact descriptions, once, at construction.**
    ///
    /// Returns what may be carried: the coefficients when they describe the same plane the witness
    /// spans, and the normal under the weaker agreement (parallel — the direction is `frame_sign`'s
    /// to record). `None` for a rotated plane, which has no exact `f64` coefficients at all.
    ///
    /// ★ **One function, so a test fixture cannot route differently from the arrangement.** Filling
    /// the two fields by hand at a second construction site is how the fixture and the engine come
    /// to disagree about which planes are describable — and this whole item exists because two
    /// descriptions of one plane disagreed.
    pub(crate) fn reconcile(
        plane: &Plane,
        tri: [Point3; 3],
        rotated: bool,
    ) -> (Option<[f64; 4]>, Option<[f64; 3]>) {
        if rotated {
            return (None, None);
        }
        let c = plane.coefficients();
        (
            plane.spans_exactly(tri).then_some(c),
            plane.normal_spans(tri).then(|| [c[0], c[1], c[2]]),
        )
    }

    /// The class's outward normal — the root face's, which is what `tri` is wound for and what
    /// `emit_faces` winds its rings about. Not normalized: only its direction is ever read.
    pub(crate) fn tri_n_out(&self) -> Vector3 {
        (self.tri[1] - self.tri[0]).cross(self.tri[2] - self.tri[0])
    }
}

/// Dense plane ids for a face table: `(geom, plane_ix)` where `plane_ix[face]` indexes `geom`.
///
/// **The numbering is monotone in `canon`.** Roots are ranked in increasing order, so
/// `canon[i] < canon[j]` iff `plane_ix[i] < plane_ix[j]` — every comparison, sort and lex-min over
/// plane indices is order-isomorphic to the sparse form. Nothing found in the engine turns out to
/// depend on that (the two candidates — `loop_winding`'s lex-min node and `crossings`' pre-dedup
/// sort — are by coordinate and by set, respectively), but the audit cannot be proved exhaustive
/// over ~175 sites, so the numbering removes the question instead of answering it.
pub(crate) fn dense_planes(
    planes: &[FaceRow],
    canon: &[usize],
) -> (Vec<WorkingPlane>, Vec<ClassIx>, Vec<Handle<Surface>>) {
    // Plane classes densify from the union-find roots; cylinder rows never entered the
    // union-find (their identity question is the cylinder class table's, C2), so here they
    // number by first-seen surface — the index space exists before the table does.
    let mut roots: Vec<usize> = canon
        .iter()
        .enumerate()
        .filter(|&(i, _)| matches!(planes[i], FaceRow::Plane(_)))
        .map(|(_, &c)| c)
        .collect();
    roots.sort_unstable();
    roots.dedup();
    let mut cyl_ix: HashMap<Handle<Surface>, usize> = HashMap::new();
    let plane_ix = canon
        .iter()
        .enumerate()
        .map(|(i, c)| match &planes[i] {
            FaceRow::Plane(_) => ClassIx::Plane(
                roots
                    .binary_search(c)
                    .expect("a class root is in the root set"),
            ),
            FaceRow::Cylinder(cf) => {
                let n = cyl_ix.len();
                ClassIx::Cyl(*cyl_ix.entry(cf.surf).or_insert(n))
            }
        })
        .collect();
    let geom = roots
        .iter()
        .map(|&r| {
            let pi = planes[r].plane();
            // ★ The two descriptions are reconciled **here**, once, and a description that loses
            // is simply not carried. A rotated plane has no exact `f64` coefficients at all, so
            // both are `None` there — which is also what makes the predicates stop asking
            // "is it rotated?" and ask "did I get coefficients?" instead.
            let (exact_coeffs, exact_normal) =
                WorkingPlane::reconcile(&pi.plane, pi.tri, pi.rotated);
            WorkingPlane {
                base_rat: pi.base_rat,
                base: BaseFrame::of(&pi.tri_pt3, pi.motion, pi.orient_sign, pi.base_rat),
                name_ints: nacre_cip::predicate::name_stored_ints(
                    pi.name.as_ref(),
                    &pi.tri_pt3,
                    pi.orient_sign,
                ),
                plane: pi.plane,
                surf: pi.surf,
                tri: pi.tri,
                tri_pt3: pi.tri_pt3.clone(),
                rotated: pi.rotated,
                frame_sign: pi.orient_sign,
                exact_coeffs,
                exact_normal,
            }
        })
        .collect();
    // The cylinder surfaces in `ClassIx::Cyl` numbering order (first seen) — rebuilt from the
    // same map so the two cannot drift.
    let mut cyls: Vec<(usize, Handle<Surface>)> = cyl_ix.into_iter().map(|(s, k)| (k, s)).collect();
    cyls.sort_unstable_by_key(|&(k, _)| k);
    (geom, plane_ix, cyls.into_iter().map(|(_, s)| s).collect())
}

/// Whether three `WitnessPoint` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected `RAY_DEGENERATE`, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, so this is the exact
/// zero-area (collinear) test on the vertices' rotation definitions.
#[cfg(test)]
pub(crate) fn pt3_base_collinear(a: &WitnessPoint, b: &WitnessPoint, c: &WitnessPoint) -> bool {
    use nacre_scalar::Rat;
    let (a, b, c) = (&a.base, &b.base, &c.base);
    let proj_zero = |i: usize, j: usize| -> Option<bool> {
        let det = b[i]
            .checked_sub(a[i])?
            .checked_mul(c[j].checked_sub(a[j])?)?
            .checked_sub(
                b[j].checked_sub(a[j])?
                    .checked_mul(c[i].checked_sub(a[i])?)?,
            )?;
        Some(det == Rat::from_int(0))
    };
    matches!(
        (proj_zero(1, 2), proj_zero(2, 0), proj_zero(0, 1)),
        (Some(true), Some(true), Some(true))
    )
}

/// Each outer-shell edge with its bound vertices and the two combined-plane
/// indices of its adjacent faces, in first-seen (deterministic) order.
///
/// Every loop of every face is walked, holes included: a hole-ring edge is used
/// once by the holed face's inner loop and once by the neighbouring wall's outer
/// loop, so it too has exactly two incident faces. Walking `outer` before `inner`
/// on each face leaves the order of a hole-free solid untouched.
///
/// The pair is returned as `[usize; 2]`, so no caller can index a third slot: an
/// edge with any other incidence count is a non-manifold shell and rejects here.
#[allow(clippy::type_complexity)]
pub(crate) fn edge_incidence(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<Vec<(Handle<Edge>, [Handle<Vertex>; 2], [usize; 2])>, BoolError> {
    let mut order: Vec<Handle<Edge>> = Vec::new();
    let mut map: HashMap<Handle<Edge>, ([Handle<Vertex>; 2], Vec<usize>)> = HashMap::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&fh];
            for he in face_half_edges(face) {
                let bounds = model.edges.get(he.edge).vertices;
                let entry = map.entry(he.edge).or_insert_with(|| {
                    order.push(he.edge);
                    (bounds, Vec::new())
                });
                entry.1.push(pidx);
            }
        }
    }
    order
        .into_iter()
        .map(|e| {
            let (b, p) = map.remove(&e).unwrap();
            match p[..] {
                [x, y] => Ok((e, b, [x, y])),
                // `validate` would call this `NonOpposedEdge`, but `boolean` never runs
                // `validate` on its inputs, so the guard stays. No firing test.
                _ => Err(reject(RejectReason::NonManifoldEdge)),
            }
        })
        .collect()
}

/// Every half-edge of a face: its outer loop first, then each hole ring in order.
pub(crate) fn face_half_edges(face: &Face) -> impl Iterator<Item = &HalfEdge> {
    face.outer
        .half_edges
        .iter()
        .chain(face.inner.iter().flat_map(|l| l.half_edges.iter()))
}

/// Two faces lie on the same plane — by a **shared `Surface` handle** (explicit
/// sharing: O(1) `Handle` identity, exact, rotation-independent) or, as a fallback,
/// by the geometric rank-1 `planes_coplanar` test. A referenced coplanar contact —
/// a pad/pocket cap that reuses its face's surface — is caught by the handle path
/// without any coordinate test. On the axis-aligned M5 corpus the handle path is
/// redundant with `planes_coplanar` (same handle ⇒ same plane), so the geometric
/// fallback is what keeps independently-built coplanar contacts working; the handle
/// path's real payoff is rotated frames, where the geometric test would need the
/// rotation-exact judgment.
pub(crate) fn shares_or_coplanar(jd: &Judge<'_, FaceRow>, i: usize, j: usize) -> bool {
    let (pa, pb) = (jd.planes[i].plane(), jd.planes[j].plane());
    // Three independent witnesses, OR-ed, so this can only ever merge *more* than before:
    //  1. the same `Surface` handle — coplanar by reference (what an ops-built tool's base cap and
    //     its target face share, and what a chained operand's split coplanar faces share);
    //  2. exactly proportional coefficients — the original test, kept;
    //  3. the faces' own coordinates, exactly (`Judge::planes_coplanar`) — the only one of the three
    //     that does not read a *derived* value, and the one that catches two independently built
    //     solids whose walls coincide (`add_cuboid` stacked on `add_cuboid`), where the rounded
    //     coefficients of differently-sized faces are not exactly proportional.
    pa.surf == pb.surf || planes_coplanar(&pa.plane, &pb.plane) || jd.planes_coplanar(i, j)
}

/// Union-find root of `x` in `parent` (with path compression). Roots are the smallest index
/// of their class, so the result is deterministic (replay, DNA §absolute-3).
/// Union-find root with path compression. Drives component grouping in
/// [`crate::boolean::unify_coplanar_faces`].
pub(crate) fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut c = x;
    while parent[c] != r {
        let next = parent[c];
        parent[c] = r;
        c = next;
    }
    r
}

/// Canonicalize the combined plane table by coplanarity: two planes that are the same plane
/// (shared `Surface` handle, or exact rank-1 [`planes_coplanar`]) are merged into one class, so
/// a wall of `a` coplanar with a wall of `b` names a **single line** in a shared plane π. This is
/// the one thing the seam engine cannot do (it rejects `order_along(R,R)==0` as `FOURPLANE`);
/// canonicalizing turns that self-comparison into a real order. Returns `canon` where `canon[i]`
/// is the class root (the smallest index in the class). Every decision is exact
/// (`shares_or_coplanar`) — no coordinate. O(n²) scan over the (small) face count.
pub(crate) fn plane_classes(jd: &Judge<'_, FaceRow>) -> Vec<usize> {
    let planes = jd.planes;
    let n = planes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    // Cylinder rows stay self-rooted: their identity question belongs to the cylinder class
    // table (C2), not the coplanarity union-find — they simply never enter `reps` below.

    // **Merge by `Surface` handle first, then compare only one face per distinct surface.**
    //
    // Coplanarity is an equivalence, and a shared handle *is* the same plane — exactly, by
    // identity, in O(1). So faces that share one are already one class, and asking the
    // geometric question of each of them separately asks the same question many times.
    //
    // The answer cannot move: union-find's result is the transitive closure of the pairs it was
    // given, and a class's root stays its minimum index, so dropping *redundant* pairs leaves
    // the partition and its roots alone. No tolerance is involved — a shared handle is identity.
    //
    // **Measured, and smaller than the pair count suggests.** On the 80-fin fold's largest
    // boolean, 406 faces carry 172 distinct surfaces: 82,215 pairs become 14,706, and over the
    // fold 7.52M become 2.30M — 3.3x fewer. But the scan only got ~1.4x faster (0.96s → 0.69s
    // over the fold, ~1.03x end to end), because **the pairs this drops are the cheapest ones**:
    // they matched on the handle and returned at the first `||`. What is left is the
    // geometrically distinct pairs, which are the ones that were expensive all along. Counting
    // removed operations overstates the saving whenever the removed ones are the cheap ones.
    let mut rep: HashMap<Handle<Surface>, usize> = HashMap::new();
    let mut reps: Vec<usize> = Vec::new();
    for (i, row) in planes.iter().enumerate() {
        let FaceRow::Plane(p) = row else { continue };
        match rep.get(&p.surf) {
            Some(&r) => {
                let (ri, rj) = (uf_find(&mut parent, r), uf_find(&mut parent, i));
                if ri != rj {
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
            None => {
                rep.insert(p.surf, i);
                reps.push(i);
            }
        }
    }

    for a in 0..reps.len() {
        for b in (a + 1)..reps.len() {
            let (i, j) = (reps[a], reps[b]);
            if shares_or_coplanar(jd, i, j) {
                let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                if ri != rj {
                    // Attach the larger root under the smaller so a class's root is its min index.
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }
    (0..n).map(|i| uf_find(&mut parent, i)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_scalar::{Angle, Axis, Rat};

    /// **The branch tolerance is a measurement, and each of its terms can move.**
    ///
    /// Hand geometry — no class indices, so nothing here is copied from an engine run: the turned
    /// boss's second crossing, where the meet line of `y = 0` and `x = 4` pierces a cylinder of
    /// radius `0.5` about the `+X` axis through `(y, z) = (0.25, 2)`. The true point is
    /// `(4, 0, 2 − √3/4)` (irrational on purpose — every term sees real rounding, not exact zeros).
    ///
    /// ★ The second probe moves **along the meet line**: both planes and the line stay at zero, so
    /// only the cylinder-surface term can see it — which is how this asserts that term is
    /// *exercised*, not merely present ([[instrument-decides-the-answer]]'s negative control).
    #[test]
    fn branch_vertex_tol_measures_and_each_term_moves() {
        let pa = Plane::from_point_normal(
            Point3::from_array([0.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        )
        .unwrap();
        let pb = Plane::from_point_normal(
            Point3::from_array([4.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        )
        .unwrap();
        let cyl = nacre_geom::Cylinder::from_axis(
            Point3::from_array([0.0, 0.25, 2.0]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            0.5,
        )
        .unwrap();
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        let p = Point3::from_array([4.0, 0.0, s]);
        assert!(
            branch_vertex_tol(p, &pa, &pb, &cyl) < 1e-12,
            "the true crossing measures at rounding scale: {}",
            branch_vertex_tol(p, &pa, &pb, &cyl)
        );
        // Along the meet line: planes and line stay zero, the cylinder term alone answers.
        let along = Point3::from_array([4.0, 0.0, s - 1e-6]);
        assert!(
            branch_vertex_tol(along, &pa, &pb, &cyl) > 1e-7,
            "the cylinder-surface term is exercised"
        );
        // Off a plane: the instrument moves by the full offset.
        let off = Point3::from_array([4.0, 1e-6, s]);
        assert!(
            branch_vertex_tol(off, &pa, &pb, &cyl) > 0.9e-6,
            "a plane term is exercised"
        );
    }

    /// **A loop with points that do not turn still states its own outward direction.**
    ///
    /// Two faces sharing an edge list the same vertices along it, so a pad that covers part of a
    /// face leaves the neighbouring walls with points strung along one line. Those points cannot
    /// be removed — that would leave a T-vertex — so the population is every partially-covering
    /// pad, and `outer_tri` has to survive it.
    ///
    /// ★ **Turned coordinates are the whole difficulty.** Axis-aligned collinear points cancel to
    /// exactly `0.0`, which the old "first corner that is not flat" rule skipped correctly. Turned,
    /// the cancellation leaves ~2⁻⁵³ — not zero, so it was taken, and the direction that came back
    /// was rounding: measured 90° off its own plane, which is `orient_sign`, the winding of the
    /// witness triangle every predicate borrows, and which side of the plane holds material.
    ///
    /// Asserted on `collect_planes`' own output rather than through `debug_assert`, so it is a
    /// measurement in release too. `tests/collinear_loop_points.rs` sweeps the same proposition
    /// from outside the crate over a band of angles.
    #[test]
    fn a_loop_whose_points_do_not_all_turn_still_faces_the_right_way() {
        use crate::{OpOutput, Operation, Profile2d, SketchFrame, apply};
        use nacre_math::Point2;
        use nacre_scalar::{Isometry, Rotation};

        let rect = |a: f64, b: f64, c: f64, d: f64| {
            Profile2d::polygon(vec![
                Point2::from_array([a, b]),
                Point2::from_array([c, b]),
                Point2::from_array([c, d]),
                Point2::from_array([a, d]),
            ])
            .expect("a rectangle is a fair profile")
        };

        let mut m = Model::new();
        let world = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: world,
                profile: rect(0.0, 0.0, 2.0, 2.0),
                dist: 1.0,
            },
        )
        .expect("the block") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
                }),
            },
        )
        .expect("the turn") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        // A pad covering part of the face: the walls around it inherit the split edge.
        let face = m.shells.get(m.solids.get(solid).outer).faces[0];
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face,
                profile: rect(0.25, 0.25, 1.25, 1.25),
                dist: 1.0,
            },
        )
        .expect("the pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();

        let faces = collect_planes(&m, solid).expect("the planes of a padded block");
        // The population has to be present, or this asserts over nothing.
        let straight = faces
            .iter()
            .filter(|f| {
                let n = m
                    .faces
                    .get(f.face().expect("a model face"))
                    .outer
                    .half_edges
                    .len();
                n > 4
            })
            .count();
        assert!(
            straight > 0,
            "no face carries a split edge — the fixture stopped reaching the population"
        );
        // The lock is on the *winding* side now: `n_out` is read off the stored
        // orientation, so `plane.normal()·n_out` is ±1 by construction and asserts
        // nothing. What the widest-corner rule still owns is the triangle — revert
        // `outer_tri` to "first non-flat corner" and the cross below is rounding
        // noise again, failing this in release.
        for row in &faces {
            let f = row.plane();
            let cos = (f.tri[1] - f.tri[0])
                .cross(f.tri[2] - f.tri[0])
                .normalize()
                .expect("a widest corner spans area")
                .dot(f.n_out);
            assert!(
                cos > 0.5,
                "a face's outer triangle is {cos:.3} of the way to perpendicular against its \
                 stated outward"
            );
        }
    }

    /// **A judgement's headroom is relative to its model, not carved out of a shared ceiling.**
    ///
    /// The two budgets answer different questions — how deep the model is, and how thin a witness
    /// it may still judge — and sharing one absolute number silently couples them: a solid turned
    /// a thousand times would leave a hard judgement no room, so the same sliver would be judged
    /// in a fresh model and abandoned in a turned one. Turning a model must not change what
    /// counts as judgeable, so this pins the headroom to a constant *above the model's own
    /// precision*, at both ends of the depth range.
    #[test]
    fn the_climbing_headroom_survives_a_deep_model() {
        let deg = Angle::from_deg(Rat::from_int(37)).expect("angle");
        let mut p = WitnessPoint::at([Rat::from_int(1), Rat::from_int(2), Rat::from_int(3)]);
        let mut seen = Vec::new();
        for turn in 0..=200 {
            if turn == 0 || turn == 20 || turn == 200 {
                let j = standard_for_points(std::slice::from_ref(&p));
                assert_eq!(
                    j.cap,
                    j.prec + CLIMB_HEADROOM,
                    "turn {turn}: headroom is not the model's own precision plus a constant"
                );
                seen.push(j.prec);
            }
            p = p.rotate(Axis::Z, deg);
        }
        // …and the precision really did grow with the history, or the test above is vacuous.
        assert!(
            seen[0] < seen[2],
            "precision did not grow with the rotation history: {seen:?}"
        );
    }
}
