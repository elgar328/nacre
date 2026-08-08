//! The plane/face substrate: per-face (`FaceInfo`) and per-plane-class (`PlaneGeom`) tables and
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
use nacre_topo::{Edge, Face, HalfEdge, Model, Orientation, Shell, Solid, Vertex};
use std::collections::HashMap;

/// A face's supporting plane plus the exact in/out data the seam path needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` happens to equal that RH normal
/// only because `tri` is taken outer-CCW; `plane.normal()` is the *surface's*
/// normal and may point inward on a `Reversed` face. Every sign test here reads
/// `n_out` (or `tri`), and none reads `plane.normal()`.
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
    /// Outward normal, `(tri[1]−tri[0])×(tri[2]−tri[0])` normalized — the single
    /// source of "outward" for both the in/out sign test and face ordering.
    pub(crate) n_out: Vector3,
    /// `+1` when this face's stored plane normal already points out of its solid, `-1` when the
    /// face is `Reversed` and the two oppose.
    ///
    /// **This face's**, not its plane class's. The class-frame twin is [`PlaneGeom::frame_sign`],
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
    /// Whether this face's plane is a *moved image* — the predicate-routing signal, read from
    /// the surface's own truth (`Model::surface_truth`).
    ///
    /// **Set together with `tri_pt3`, and only here.** It used to be decided per solid, by asking
    /// the vertices — which a boolean's result cannot answer, since its vertices are all
    /// `Discovered`. The surface answers for itself, and one solid can hold both kinds at once
    /// (fuse an axis-aligned hub with a turned fin).
    pub(crate) rotated: bool,
}

/// The supporting planes of a solid's outer shell. `Unsupported` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<FaceInfo>, BoolError> {
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let plane = match model.surface(face.surface) {
                Surface::Plane(p) => *p,
                Surface::Cylinder(_) => return Err(reject(RejectReason::CylinderFace)),
            };
            let (tri, _) =
                outer_tri(model, face).ok_or_else(|| reject(RejectReason::DegenerateFace))?;
            let n_out = (tri[1] - tri[0])
                .cross(tri[2] - tri[0])
                .normalize()
                .ok_or_else(|| reject(RejectReason::DegenerateNormal))?;
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
                // The rational-closure branch is the whole of stage 1: three vertices that solve
                // to `Rat` give a triangle indistinguishable from a stated one, so every predicate
                // below runs unchanged. The other branch — mixed frames, irrational motion — has
                // no witness triangle at all; the producer refuses to build such a plane, so
                // reaching here means the invariant broke rather than the user asked for something
                // unsupported.
                //
                // ★★★ **Corrected 2026-08-08: this is the *second* wall, not the first.** It used
                // to say the homogeneous-lifting stage opens that branch. It does not on its own —
                // such a plane has no name, so it never gets a `SketchFrame`, never becomes a base
                // cap, and never reaches this table (`OpError::VerticesInMixedFrames` has the
                // chain). The order is: give a nameless plane a frame, *then* this arm needs the
                // homogeneous route, and that is where `nacre-cip`'s `plane_iv_through` finally
                // has a caller.
                nacre_topo::SurfaceTruth::Plane {
                    points: nacre_topo::PlanePoints::Through(vs),
                    motion,
                } => {
                    let motion = *motion;
                    let _t = Watch::new();
                    let base = model
                        .through_points_rat(*vs)
                        .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                    let w = match motion {
                        None => base.map(WitnessPoint::at),
                        Some(m) => {
                            let chain = crate::rotated_vertex::motion_chain(model, m)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                            let mut out = [
                                WitnessPoint::at(base[0]),
                                WitnessPoint::at(base[1]),
                                WitnessPoint::at(base[2]),
                            ];
                            for (o, b) in out.iter_mut().zip(base) {
                                *o = crate::rotated_vertex::replay(WitnessPoint::at(b), &chain)
                                    .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                            }
                            out
                        }
                    };
                    _t.charge(Sub::TriPt3);
                    (wind(w), motion.is_some(), motion)
                }
                // The cache match above already rejected cylinders.
                nacre_topo::SurfaceTruth::Cylinder { .. } => {
                    return Err(reject(RejectReason::CylinderFace));
                }
            };
            // `orient_sign`, precomputed: the two invariants it used to re-check on every call
            // are properties of this face, so they are decided once, here.
            let dot = plane.normal().dot(n_out);
            debug_assert!(
                dot.abs() > 0.5,
                "a plane's normal must be parallel to n_out"
            );
            debug_assert_eq!(
                dot > 0.0,
                face.orientation == Orientation::Forward,
                "n_out's sign against the surface normal is the face's orientation"
            );
            out.push(FaceInfo {
                base_rat: model
                    .surface_name
                    .get(&face.surface)
                    .and_then(|n| n.narrow())
                    .copied(),
                surf: face.surface,
                face: Some(fh),
                plane,
                tri,
                n_out,
                orient_sign: if dot > 0.0 { 1 } else { -1 },
                tri_pt3,
                rotated,
                motion,
            });
        }
    }
    Ok(out)
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

/// Three non-collinear points of a face's outer loop — with the **vertex handle** each
/// point came from — ordered so their right-hand normal points **out** of the solid.
/// The handles let the toleranced predicates rebuild each point as a `WitnessPoint` (overhaul
/// stage 3); the coordinates alone drive the axis-aligned path.
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
    let i = (0..n).find(|&i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        (b - a).cross(c - a).norm() > 0.0
    })?;
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
fn class_owners(plane_ix: &[usize], n_a: usize, n_class: usize) -> Vec<Option<SolidSide>> {
    let mut out: Vec<Option<SolidSide>> = vec![None; n_class];
    let mut seen = vec![false; n_class];
    for (fi, &c) in plane_ix.iter().enumerate() {
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
    pub(crate) planes: Vec<FaceInfo>,
    pub(crate) surf_ix: HashMap<Handle<Face>, usize>,
    pub(crate) inc_a: combinatorics::EdgeFaces,
    pub(crate) inc_b: combinatorics::EdgeFaces,
    /// Where `a`'s faces end and `b`'s begin in `planes`. The concatenation always created this
    /// boundary; it was just never written down, so every later "whose face is this?" had to
    /// rebuild it.
    pub(crate) n_a: usize,
    /// The arrangement's planes, densely indexed — see [`dense_planes`].
    pub(crate) geom: Vec<PlaneGeom>,
    /// `plane_ix[face]` is that face's plane, as an index into `geom`.
    pub(crate) plane_ix: Vec<usize>,
    /// Whose faces each plane class carries — see [`class_owners`].
    pub(crate) class_owner: Vec<Option<SolidSide>>,
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
        surf_ix.insert(pi.face.expect("collect_planes yields real faces"), i);
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
    let (geom, plane_ix) = dense_planes(&planes, &canon);
    let class_owner = class_owners(&plane_ix, n_a, geom.len());
    t.charge(Sub::Dense);
    Ok(PlaneSetup {
        planes,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        geom,
        plane_ix,
        class_owner,
        standard,
        notes,
    })
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
fn standard_for(planes: &[FaceInfo]) -> Standard {
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
        planes
            .iter()
            .filter(|p| p.rotated)
            .flat_map(|p| p.tri_pt3.iter()),
    );
    // `scale`, by contrast, is every point's business: it is the model's size, and an unmoved face
    // is as far from the origin as any other.
    standard_from(planes.iter().flat_map(|p| p.tri_pt3.iter()), worst)
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
/// direction fact anyone needs from it is [`PlaneGeom::frame_sign`]. Everything else here is a
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

pub(crate) struct PlaneGeom {
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
}

impl PlaneGeom {
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
pub(crate) fn dense_planes(planes: &[FaceInfo], canon: &[usize]) -> (Vec<PlaneGeom>, Vec<usize>) {
    let mut roots: Vec<usize> = canon.to_vec();
    roots.sort_unstable();
    roots.dedup();
    let plane_ix = canon
        .iter()
        .map(|c| {
            roots
                .binary_search(c)
                .expect("a class root is in the root set")
        })
        .collect();
    let geom = roots
        .iter()
        .map(|&r| {
            let pi = &planes[r];
            // ★ The two descriptions are reconciled **here**, once, and a description that loses
            // is simply not carried. A rotated plane has no exact `f64` coefficients at all, so
            // both are `None` there — which is also what makes the predicates stop asking
            // "is it rotated?" and ask "did I get coefficients?" instead.
            let (exact_coeffs, exact_normal) = PlaneGeom::reconcile(&pi.plane, pi.tri, pi.rotated);
            PlaneGeom {
                base_rat: pi.base_rat,
                base: BaseFrame::of(&pi.tri_pt3, pi.motion, pi.orient_sign, pi.base_rat),
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
    (geom, plane_ix)
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
pub(crate) fn shares_or_coplanar(jd: &Judge<'_, FaceInfo>, i: usize, j: usize) -> bool {
    let (pa, pb) = (&jd.planes[i], &jd.planes[j]);
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
/// Union-find root with path compression. Drives component grouping in [`unify_coplanar_faces`].
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
// Wired into the unified coplanar handler's dispatch in a later cell; used by tests now.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn plane_classes(jd: &Judge<'_, FaceInfo>) -> Vec<usize> {
    let planes = jd.planes;
    let n = planes.len();
    let mut parent: Vec<usize> = (0..n).collect();

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
    for (i, p) in planes.iter().enumerate() {
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
