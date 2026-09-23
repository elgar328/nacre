//! The plane/face substrate: per-face (`FaceInfo`) and per-plane-class (`WorkingPlane`) tables and
//! their construction. Everything the boolean engine and its combinatorial queries build on.

use crate::{BoolError, RejectReason, he_start, reject};
use nacre_exact::Mag;
use nacre_geom::Plane;
use nacre_geom::intersect::{plane_plane, planes_coplanar};
use nacre_judge::predicate::Judge;
use nacre_judge::{Standard, WitnessPoint};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, HalfEdge, Model, Shell, Solid, Surface, Vertex};
use std::collections::HashMap;

mod corners;
mod cyl_gate;
mod cyl_geom;
mod frames;
mod incidence;
mod standard;
mod table;

pub(crate) use corners::*;
pub(crate) use cyl_gate::*;
pub(crate) use cyl_geom::*;
pub(crate) use frames::*;
pub(crate) use incidence::*;
pub(crate) use standard::*;
pub(crate) use table::*;

/// A face's supporting plane plus the exact in/out data the seam path needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` is the face's *stated*
/// outward — `plane.normal()` × `orientation`, the reading props and STEP trust
/// — and the loop's winding is held to it by a `debug_assert` at construction
/// and by `validate`'s `FaceMisoriented` at every op. Every sign test here reads
/// `n_out` (or `tri`), and the two agree by that enforcement.
/// One row of the boolean's face table — the face vocabulary the engine reads.
///
/// The row is an enum so a cylinder face can *sit in the table* (keeping the face-index space that
/// `surf_ix`/`EdgeFaces`/`plane_ix` share) while the plane data keeps its own struct — plane
/// consumers read through [`FaceRow::plane`], and the population gate decides what flows.
/// ★ The size gap is the `Surface` trade taken again: planes dominate every table (a
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

    /// The row's surface, whichever kind it is. Only tests ask this of a row.
    #[cfg(test)]
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

/// A cylinder face's table row: what `collect_planes` can state about a lateral face
/// without pretending it has a plane's four-piece description (`plane`/`tri`/`n_out`/`tri_pt3`
/// are constant-normal vocabulary — a dummy would be the type lying). The class index and the
/// exact def live in the cylinder class table (`ClassIx::Cyl`).
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
    /// The cylinder's exact statement **in the world** — cloned here so the tracer (which works
    /// off the face table, never the `Model`) can ask ⊥-ness and axis parameters. `None` when the
    /// truth is written in a frame no fold carries out exactly (a frame node, a turn off the
    /// quarters, or overflow); consumers decline by their own names, and the population gate
    /// refuses such a boolean before the arrangement runs. See
    /// [`nacre_topo::Model::world_cylinder_def`].
    pub(crate) def: Option<nacre_topo::CylinderDef>,
    /// This lateral **face**'s extent in the axis parameter `t` (of the raw `def.dir()`): the
    /// least and greatest `t` of the ⊥ carriers on its outer loop — `t = −(n·o + d)/(n·m)` per cap
    /// plane. `None` when a ⊥ carrier has no narrow rational name, or the loop has fewer than two
    /// distinct ⊥ stations.
    ///
    /// ★★★★★ **This is a range, not a promise about the loop's shape.** Whether the face is a
    /// band — two whole rims and holes — is `combinatorics::FaceLoops::cycles`' question
    /// (`arrangement::lateral_shape`), and the two lateral roads ask it there before they read anything; a face
    /// that is not a band (a panel, a chain rim) declines by name (`CylSpan`) until the chart can
    /// hold it. The range is what the population gate's rectangle reads
    /// (`lateral_spans`) — a wider rectangle only refuses more — and what the chart's rows carry
    /// once the tracer has cut every circle the face covers only partly.
    ///
    /// ★★ It is read from the ⊥ carriers' stations, not from "the two rims of the outer loop":
    /// a hole spliced into the outer walk, a panel or a chain rim has no such pair.
    ///
    /// ★ This range is the first axis of the face's [`Footprint`] on its own chart; the
    /// second, θ, joins with the angular extent.
    pub(crate) footprint: Footprint,
}

/// **A lateral face's footprint on its own chart `(t, θ)`** — the bounding rectangle of
/// the region the face occupies there (a lateral *is* a region of its chart). Conservative for
/// a chain rim or a notched face: wider, never narrower, so a clearance proved against it holds
/// for the face. Every clearance the gate asks of a lateral is a two-axis question against this
/// rectangle — the shape [`face_clears_footprint`] already has for a plane's faces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Footprint {
    /// The axis-parameter extent — [`CylFaceInfo::footprint`]'s doc says how it is read.
    pub(crate) span: Option<[nacre_exact::Rat; 2]>,
    /// The angular extent — the arc the face's rims trace, as two radial vectors of the
    /// cylinder's radius from the axis, `from → to` counter-clockwise about the axis direction
    /// ([`lateral_theta_extent`]). `None` is a whole circle, or an extent this road could not
    /// state (an irrational corner, rims that do not chain into one arc) — read as the whole
    /// circle, which is conservative.
    pub(crate) theta: Option<RimArc>,
}

/// An arc of a cylinder's cross-section, by two **radial vectors** of the cylinder's radius — a
/// rim point minus the axis point on its cap — `from → to` counter-clockwise about the axis
/// direction, `from == to` never (a whole circle is [`Footprint::theta`]'s `None`). Rational, as
/// every corner a prism's rim has is: a seam point when the reference direction's norm is
/// (`inv_sqrt_exact`), a pierce corner when its root is ([`nacre_exact::quad::QuadVal::as_rat`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct RimArc {
    pub(crate) from: [nacre_exact::Rat; 3],
    pub(crate) to: [nacre_exact::Rat; 3],
}

/// Whether the radial direction `x` lies on the arc `from → to` (counter-clockwise about `m`,
/// ends included). The sign of `(u × v)·m` says whether `v` is counter-clockwise of `u` by less
/// than a half turn; an arc of more than a half turn is the complement of the shorter way round —
/// the same three-way reading `nacre_geom::mixed`'s arc containment makes in 2D.
fn arc_contains(
    arc: &RimArc,
    x: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
) -> Option<bool> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let turn = |u: &[Rat; 3], v: &[Rat; 3]| -> Option<Rat> {
        nacre_exact::dot3_rat(&nacre_exact::cross3_rat(u, v)?, m)
    };
    let ft = turn(&arc.from, &arc.to)?;
    let fx = turn(&arc.from, x)?;
    let xt = turn(x, &arc.to)?;
    Some(if ft > zero {
        fx >= zero && xt >= zero
    } else if ft < zero {
        !(turn(&arc.to, x)? > zero && turn(x, &arc.from)? > zero)
    } else if nacre_exact::dot3_rat(&arc.from, &arc.to)? > zero {
        true // `from` and `to` point the same way: the whole circle
    } else {
        fx >= zero // exactly a half turn
    })
}

#[derive(Clone)]
pub(crate) struct FaceInfo {
    pub(crate) surf: Handle<Surface>,
    /// The face this plane came from. Distinguishes two coplanar faces that share one
    /// `Surface` (a Cut splits one face into disjoint pieces reusing its surface —
    /// cell coplanar-narrow), which `surf` alone collapses. `surf_ix` keys on this.
    ///
    /// ★ **`None` has no producer today.** It stands for a cap a half-space clip puts on an
    /// operand — a
    /// face that lives for one boolean and is never emitted — and nothing mints those.
    /// Every constructor writes `Some`, so every read `expect`s.
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
    /// and one function called with either kind of index would be the single place the
    /// face/plane convention cannot be asserted, because both readings are legitimate.
    /// Separate names, separate questions.
    pub(crate) orient_sign: i8,
    /// The three `tri` points as **exact `WitnessPoint` definitions**, in the same order as `tri`.
    /// Built once here and borrowed by every predicate (`plane_def`) — rebuilt per judgment, it
    /// dominates the boolean's runtime.
    pub(crate) tri_pt3: [WitnessPoint; 3],
    /// The motion-history leaf this face's plane was moved by, or `None` for a constructed one.
    /// **The canonical identity of "which motion"** — see [`BaseFrame`].
    pub(crate) motion: Option<Handle<nacre_topo::MotionNode>>,
    /// The surface's plane as **exact rational coefficients in the frame its truth names**
    /// (`Model::surface_name`) — the world when unmoved, the pre-motion frame when moved.
    /// `None` when the producer had no rational description. Read by [`BaseFrame`], which would
    /// otherwise re-derive a moved plane from its pre-motion triangle and round `d`.
    pub(crate) base_rat: Option<[nacre_exact::Rat; 4]>,
    /// The same plane as **exact rational coefficients in the world**, whatever frame the truth
    /// is written in — see [`world_plane_coeffs`]. `None` when no exact world description exists
    /// (a rotation, a frame, a wide name, an overflow). This is what the cylinder roads compare
    /// against a world axis; `base_rat` above answers the *other* question (the description in
    /// the frame the provenance names, which is what `BaseFrame` cancels).
    pub(crate) world_rat: Option<[nacre_exact::Rat; 4]>,
    /// The surface's full canonical name (`Model::surface_name`), **any width** — what
    /// [`WorkingPlane::name_ints`] is folded from. `base_rat` above is its narrow projection,
    /// kept beside it because the narrow consumers (`BaseFrame`, the composed-rotation route)
    /// read `[Rat; 4]` directly.
    pub(crate) name: Option<nacre_exact::PlaneName>,
    /// Whether this face's plane is a *moved image* — the predicate-routing signal, read from
    /// the surface's own truth (`Model::surface`).
    ///
    /// **Set together with `tri_pt3`, and only here.** The surface answers for itself, not its
    /// solid's vertices, and one solid can hold both kinds at once (fuse an axis-aligned hub with
    /// a turned fin).
    pub(crate) rotated: bool,
}

pub(crate) fn solid_shell_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Shell>> {
    let s = model.solid(solid);
    std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect()
}

/// Three **well-spread** points of a face's outer loop — with the **vertex handle** each
/// point came from — ordered so their right-hand normal points **out** of the solid.
/// ★ No production consumer reads the handles any more (the toleranced predicates take their
/// witnesses from the surface truth, `tri_pt3`); every caller keeps the coordinates only, and
/// the `None` of a loop that spreads no three points is answered by the caller — the plane's
/// truth points for a circle-bounded face (`collect_planes`), "no evidence" for
/// `find_face_coplanar_with`.
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
    // sibling that answers this same question already says so: `construct::Swept::cap_points` takes "the
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
    // Same b/c swap for coords and handles, so the k-th point and the k-th handle stay aligned.
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

/// [`vertex_tol`]'s pierce sibling — the measured tolerance of a realized `plane ∩ plane ∩
/// cylinder` vertex: max distance of `p` to its two planes, the cylinder **surface**, and the
/// planes' meet line.
///
/// ★★ **The meet line is the one pairwise curve this covers, on purpose.** The seam's tolerance
/// means "surfaces *and* pairwise meets" — `boolean`'s vertex minting says the pairwise part "is
/// a real part of what this number means", and [`vertex_tol`] above covers all three of its
/// lines. Here only `a ∩ b` has a closed form that is always a line: a pierce plane can be
/// **parallel to the axis**, where `plane ∩ cylinder` is a pair of ruling lines — a per-case
/// curve family this deliberately does not chase. The cylinder-*surface* distance already bounds
/// the radial part of that error; the meet line is also exactly the line the point is defined on
/// (`combinatorics::pierce_point` realizes from `(line, s)`), so its residual is the first-class
/// question about the realization.
///
/// ★ **The meet-line term is unexercised in today's corpus, measured** — with it removed, every
/// assertion stays green. Both fixtures' pierce planes are perpendicular, and for a
/// perpendicular pair the line residual never exceeds `√2 ×` the larger plane residual, so the
/// `max` cannot turn on it. It earns its keep the day a pierce pair meets **obliquely** (a
/// turned wall), where the line residual outgrows both plane residuals near the line.
pub(crate) fn pierce_vertex_tol(
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
pub(crate) fn class_owners(
    plane_ix: &[ClassIx],
    n_a: usize,
    n_class: usize,
) -> Vec<Option<SolidSide>> {
    let mut out: Vec<Option<SolidSide>> = vec![None; n_class];
    let mut seen = vec![false; n_class];
    for (fi, &ci) in plane_ix.iter().enumerate() {
        // Cylinder classes get their own owner table with the cylinder class table.
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

/// A face's class in the arrangement — **which index space the face's surface lives in**.
/// The old `plane_ix: Vec<usize>` presumed every class is a plane; the enum makes a cylinder
/// class unrepresentable as a plane index instead of smuggling it through a sentinel.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum ClassIx {
    /// An index into the dense plane-class table (`arrangement::PlaneSetup::geom`).
    Plane(usize),
    /// An index into the cylinder-class table. The index space predates the table — it
    /// existed first so the type was total before the table had entries.
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

/// A cylinder class of one boolean: the exact statement the population gate reasons
/// about, beside its f64 cache. One entry per distinct lateral surface, in [`ClassIx::Cyl`]
/// numbering order.
pub(crate) struct WorkingCyl {
    pub(crate) surf: Handle<Surface>,
    pub(crate) def: nacre_topo::CylinderDef,
    /// The f64 twin of `def` — its realization, what a *measurement* reads (`pierce_vertex_tol`
    /// measures a pierce realization against this surface), while every decision reads `def`. Not
    /// a cache of the model: it lives for one operation, like `WitnessPoint::realized`.
    pub(crate) realized: nacre_geom::Cylinder,
    /// Which operand states this class. The pair loop asks only pairs of **different**
    /// owners: two classes of one valid solid keep their faces apart by construction, and the
    /// arrangement has no cylinder–cylinder road that would need the proof. One surface stated by
    /// both solids — the coincident pair by another spelling — never reaches this table: the class
    /// loop refuses it by name before pushing.
    pub(crate) owner: SolidSide,
}

/// **The interval a lateral face reaches along a direction `d`**, in `d·p` units.
///
/// A point of the face is `p = o + s·m + r·u` with `u ⊥ m` a unit vector and `s` over the face's
/// own span, so `d·p = d·o + s·(d·m) + r·(d·u)` with `(d·u)² ≤ |d⊥|² = d·d − (d·m)²/(m·m)`: the
/// face's projected span widened by a radial reach `ρ = r·|d⊥|`, carried as its square so no
/// root is ever formed. It is the one clearance question a lateral face answers for any other
/// surface: a cylinder pair asks two of these against each other, once per direction its shape
/// can state ([`separating_dirs`], [`reaches_apart`]); an oblique plane reads `d = n` against its
/// single station.
///
/// `span = None` is a face whose span could not be stated ([`lateral_spans`] empty). The reach
/// is still bounded when `d·m = 0` — the projection is a point whatever `s` is — and unbounded
/// otherwise, which is `None`: nothing proved. ★ That `d·m = 0` case is not a branch of its own;
/// it is the general formula with the `s` term vanishing — and it is what makes the **common
/// perpendicular** answerable at all (`d ⊥ m` for both cylinders, so neither needs a span).
/// The `d·m ≠ 0` arm is what the oblique plane arm reads (`d = n`) and what a cylinder
/// reads along **its own** axis, where the radial term vanishes instead and the reach is the span
/// itself. `None` is also `Rat` overflow.
/// The reach is `[lo − √rho2_lo, hi + √rho2_hi]`: each end is a rational base and a radical the
/// arc may or may not add. ★ With an angular extent the radial term `r·(d·û)` over the
/// face's arc peaks at the direction of `d⊥` when the arc holds it — `√ρ²`, as before — and at an
/// **end** of the arc otherwise, where it is `d·v` for the end's radial vector, a rational folded
/// into the base with a zero radical. One root at most on each end, and no new arithmetic.
struct Reach {
    lo: nacre_exact::Rat,
    hi: nacre_exact::Rat,
    rho2_lo: nacre_exact::Rat,
    rho2_hi: nacre_exact::Rat,
}

#[cfg(test)]
#[path = "../tests/planes.rs"]
mod tests;
