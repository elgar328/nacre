//! The plane/face substrate: per-face (`FaceInfo`) and per-plane-class (`WorkingPlane`) tables and
//! their construction. Everything the boolean engine and its combinatorial queries build on.

use crate::{BoolError, RejectReason, he_start, reject};
use nacre_exact::Mag;
use nacre_geom::Plane;
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
/// An `orient3d` against a face answers `+1` when the implicit point lies on its **witness
/// triangle's right-hand-normal side** — the convention is tied to the triangle (`tri_pt3`), never
/// to `plane`. The triangle is wound to the face's *stated* outward — the way its plane faces (the
/// truth's `sense`) × `orientation`, the reading props and STEP trust — by those truth bits
/// (`collect_planes`), and `validate`'s `FaceMisoriented` holds the loop's winding to the same
/// statement at every op.
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
/// rectangle — the shape [`face_clears_footprint`] already has for a plane's faces. Whether a line
/// along the axis lies on the face asks θ alone ([`Footprint::theta_holds_line`]).
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

impl Footprint {
    /// **Does the line along the axis through `p` lie within this face's angular extent?** — the
    /// face half of «this line lies on a lateral *face*, not only on its surface», which both of
    /// the gate's line records ask (`SharedRuling`, `Tangency`). Ends included: a line on the face's
    /// end ruling is the face's edge. `theta: None` holds every line (the whole circle, or an
    /// extent this road could not state). `None` is arithmetic that could not answer. The axis
    /// span is not asked here — only θ.
    ///
    /// ★ **A `false` here is true of the face.** `theta` is the union of the outer loop's rim arcs
    /// ([`lateral_theta_extent`]), and a lateral face on this road is bounded on its chart by arcs
    /// and rulings alone, so at every θ it covers, its upper and lower boundaries are arcs: the
    /// face's θ-projection is the union of its arcs', and the extent is never narrower than the
    /// face. What it does not see is a **hole** in the face — a line through a lateral's hole reads
    /// `true`, the conservative side.
    pub(super) fn theta_holds_line(
        &self,
        p: &[nacre_exact::Rat; 3],
        o: &[nacre_exact::Rat; 3],
        m: &[nacre_exact::Rat; 3],
    ) -> Option<bool> {
        let Some(arc) = &self.theta else {
            return Some(true);
        };
        // The radial direction: `p − o` less its component along the axis.
        let mut w = [nacre_exact::Rat::from_int(0); 3];
        for k in 0..3 {
            w[k] = p[k].checked_sub(o[k])?;
        }
        let along = nacre_exact::dot3_rat(&w, m)?;
        let mm = nacre_exact::dot3_rat(m, m)?;
        let q = along.checked_mul(nacre_exact::Rat::new(mm.denom(), mm.numer())?)?;
        let mut x = [nacre_exact::Rat::from_int(0); 3];
        for k in 0..3 {
            x[k] = w[k].checked_sub(q.checked_mul(m[k])?)?;
        }
        arc_contains(arc, &x, m)
    }
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

/// **A plane's statement in the world: its canonical name there, and which way that name faces.**
///
/// The name answers «where» and carries no direction (its first nonzero component is positive);
/// `sense` answers «which way» against it, from the truth
/// ([`nacre_topo::Model::world_plane_name_sense`]). One value, so a row cannot hold a world name
/// without its sense — a reader that needs a direction asks this, never the plane cache.
#[derive(Clone, Debug)]
pub(crate) struct WorldName {
    pub(crate) name: nacre_exact::PlaneName,
    /// `Forward` when `name`'s normal points the way the plane faces (its truth's `sense`).
    pub(crate) sense: nacre_topo::Orientation,
}

impl WorldName {
    /// The surface's world statement, or `None` where it has none.
    ///
    /// ★ **The name is never dropped for want of a sense.** Both come from the same truth — the
    /// sense exists wherever the world name does (the name is derived from the plane's points, and
    /// collinear points have no name) — so a missing sense beside a present name is a broken door,
    /// not an input. Dropping the name there would silently move the class merge off its name road.
    pub(crate) fn of(model: &Model, surf: Handle<Surface>) -> Option<WorldName> {
        let name = model.world_plane_name(surf)?;
        let sense = model
            .world_plane_name_sense(surf)
            .expect("a plane with a world name faces some way in the world");
        Some(WorldName { name, sense })
    }

    /// The name's exact coefficients **turned to face the way the plane does** — the name times
    /// its sense; `None` for a wide name or on overflow.
    pub(crate) fn oriented_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        let c = *self.name.narrow()?;
        match self.sense {
            nacre_topo::Orientation::Forward => Some(c),
            nacre_topo::Orientation::Reversed => {
                let zero = nacre_exact::Rat::from_int(0);
                Some([
                    zero.checked_sub(c[0])?,
                    zero.checked_sub(c[1])?,
                    zero.checked_sub(c[2])?,
                    zero.checked_sub(c[3])?,
                ])
            }
        }
    }
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
    /// `+1` when this face's stored plane normal already points out of its solid, `-1` when the
    /// face is `Reversed` and the two oppose — the `Orientation` flag as a sign.
    ///
    /// **This face's**, not its plane class's. The class-frame twin is [`WorkingPlane::frame_sign`],
    /// and one function called with either kind of index would be the single place the
    /// face/plane convention cannot be asserted, because both readings are legitimate.
    /// Separate names, separate questions.
    pub(crate) orient_sign: i8,
    /// The witness: three exact points on the plane (**`WitnessPoint` definitions** — the plane's
    /// own points, or its frame's probes), wound to this face's outward by the truth's sense and
    /// the face's orientation. Built once here and borrowed by every predicate (`plane_def`) —
    /// rebuilt per judgment, it dominates the boolean's runtime.
    pub(crate) tri_pt3: [WitnessPoint; 3],
    /// The motion-history leaf this face's plane was moved by, or `None` for a constructed one.
    /// **The canonical identity of "which motion"** — see [`BaseFrame`].
    pub(crate) motion: Option<Handle<nacre_topo::MotionNode>>,
    /// The surface's plane as **exact rational coefficients in the frame its truth names**
    /// (`Model::surface_name`) — the world when unmoved, the pre-motion frame when moved.
    /// `None` when the producer had no rational description. The judge's pre-motion routes read it
    /// (`Witness::base_coeffs_rat` — the composed-rotation route), and so does the restatement
    /// mirror in `collect_planes`.
    pub(crate) base_rat: Option<[nacre_exact::Rat; 4]>,
    /// The same plane **in the world** — its canonical name there and which way that name faces,
    /// whatever frame the truth is written in ([`WorldName`]). `None` when no exact world
    /// statement exists (a turn off the quarters, a frame, a moved wide name, an overflow). The
    /// class table reads its narrow projection ([`WorkingPlane::world_rat`]); `base_rat` above
    /// answers the *other* question (the description in the frame the provenance names, which is
    /// what `BaseFrame` cancels).
    pub(crate) world: Option<WorldName>,
    /// The surface's full canonical name (`Model::surface_name`), **any width** — what
    /// [`WorkingPlane::name_ints`] is folded from and [`BaseFrame`] orients its pre-motion
    /// coefficients by. `base_rat` above is its narrow projection, kept beside it because the
    /// narrow consumers (the composed-rotation route, the restatement mirror) read `[Rat; 4]`
    /// directly.
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
/// truth points for a circle-bounded face (`collect_planes`).
///
/// "Well spread" rather than "non-collinear" is the whole contract: the triangle is what states
/// the face's area — `collect_planes` refuses a loop that spreads none (`DegenerateFace`) and
/// checks the face's winding against it — so a nearly-flat one is not a lesser answer but a wrong
/// one. See the corner choice below.
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
    /// The f64 twin of `def` — its realization, what a construction figure reads (the rim table
    /// solves a seam vertex's figure on it), while every decision reads `def`. Not a cache of the
    /// model: it lives for one operation, like `WitnessPoint::realized`.
    pub(crate) realized: nacre_geom::Cylinder,
    /// Which operand states this class. The pair loop asks only pairs of **different**
    /// owners: two classes of one valid solid keep their faces apart by construction, and the
    /// arrangement has no cylinder–cylinder road that would need the proof. One surface stated by
    /// both solids — the coincident pair by another spelling — never reaches this table: the class
    /// loop refuses it by name before pushing.
    pub(crate) owner: SolidSide,
    /// **The lines of two plane classes that lie on this cylinder's lateral face** — the gate's
    /// record ([`SharedRuling`]), carried with the class like `owner` so every road that holds
    /// the table reads it rather than deciding it again.
    pub(crate) shared: Vec<SharedRuling>,
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
