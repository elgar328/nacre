//! The arrangement engine — the sole boolean path ([`boolean`], which `crate::boolean`
//! delegates to). It has **two arms on one job**: every plane class here, and every cylinder
//! class on its own chart ([`cyl_chart`]). Both end at the same `LocalFace` list.
//!
//! The plane arm traces each face of both solids onto every plane class as **line segments**
//! (not closed loops), splits them at their crossings, extracts the cells, nests holes, labels
//! in/out per operand, and emits the result faces. Geometry lying in a plane is *input*, not a
//! degeneracy.
//!
//! Two face cases feed the trace. **Seated** — a face whose own class is the cut class — lies wholly
//! in the plane, so its whole boundary is trace. **Transversal / mixed** — a face crossing the
//! plane, possibly with an edge lying in it — is 2D polygon-vs-line clipping on `L = W ∩ fp`, with
//! a three-valued (−/0/+) scan that collapses on-line edges instead of rejecting them. Faces with
//! **holes** are declined by the trace brick (a forced-covered run inside a hole would claim its
//! void as material, and the corpus has no case to verify it). Declining is a per-face record, never
//! a whole-solid abort — a non-empty `declined` means "this trace is incomplete, do not conclude
//! from it".

use super::*;
use crate::assembly::*;
use crate::combinatorics::{Canon3, NodeId, NodeKind, three_plane_name};
use crate::draft::*;
#[cfg(test)]
use crate::phase;
use crate::planes::*;
use crate::tolerant::Judge;
#[cfg(test)]
use crate::transform::transform;
use nacre_geom::intersect::three_planes;
use nacre_judge::Decision;
use nacre_judge::predicate::Notes;

/// `phase::timed` where the counters exist, and a plain call where they do not.
macro_rules! timed {
    ($c:ident, $e:expr) => {{
        #[cfg(test)]
        {
            crate::phase::timed(&crate::phase::$c, || $e)
        }
        #[cfg(not(test))]
        {
            $e
        }
    }};
}

/// A `phase::Watch` over the rest of the enclosing scope, and nothing at all in a release build.
/// (The module is `cfg(test)`, so there is nothing to link to in a doc build.)
macro_rules! watch {
    ($c:ident) => {
        #[cfg(test)]
        let _w = crate::phase::Watch::new(&crate::phase::$c);
    };
}

mod aliases;
#[cfg(test)]
#[path = "../tests/probes/arc_probe.rs"]
pub(crate) mod arc_probe;
#[cfg(test)]
#[path = "../tests/arrangement/audits.rs"]
mod audits;
mod cells;
#[cfg(test)]
#[path = "../tests/probes/crossing_probe.rs"]
pub(crate) mod crossing_probe;
#[cfg(test)]
#[path = "../tests/probes/cycle_probe.rs"]
pub(crate) mod cycle_probe;
/// The cylinder chart: the lateral faces' arrangement, emitted from cells —
/// this engine's second arm, run on the cylinder's own chart rather than a plane class.
pub(crate) mod cyl_chart;
mod cyl_trace;
#[cfg(test)]
#[path = "../tests/probes/decline_probe.rs"]
pub(crate) mod decline_probe;
#[cfg(test)]
#[path = "../tests/probes/disk_side_probe.rs"]
pub(crate) mod disk_side_probe;
mod emit;
#[cfg(test)]
#[path = "../tests/probes/extent_probe.rs"]
pub(crate) mod extent_probe;
#[cfg(test)]
#[path = "../tests/probes/order_probe.rs"]
pub(crate) mod order_probe;
mod per_class;
mod result;
#[cfg(test)]
#[path = "../tests/probes/ruling_probe.rs"]
pub(crate) mod ruling_probe;
mod rulings;
mod setup;
mod sides;
mod split;
mod split_circles;
mod split_rulings;
mod trace_class;
mod trace_plane;
mod traces;

pub(crate) use aliases::*;
#[cfg(test)]
pub(crate) use audits::*;
pub(crate) use cells::*;
pub(crate) use cyl_trace::*;
use emit::*;
// ★ Not a glob, because the two names crossing here cross for different reasons. The stage
// entry is this module's own (`pub(super)`), so a plain `use` carries it to the children.
// The mixed-class audit leaves `arrangement` for exactly one reader -- a test, and one that
// only runs where `debug_assert!` does -- so it is re-exported under both conditions. A
// `pub(crate)` glob claimed both at the wider visibility, which made it claim nothing at all
// in a release build and say so.
#[cfg(all(test, debug_assertions))]
pub(crate) use per_class::MIXED_CLASS_AUDIT;
use per_class::per_class;
pub(crate) use result::*;
pub(crate) use rulings::*;
pub(crate) use setup::*;
pub(crate) use sides::*;
pub(crate) use split::*;
pub(crate) use split_circles::*;
use split_rulings::*;
pub(crate) use trace_class::*;
use trace_plane::*;
pub(crate) use traces::*;

/// One feature node on the line `L = W ∩ fp` — **its own name, and what pins it there**.
///
/// ★★ Not a bare third plane (`r`), which is the shape a three-plane point has and no other: the
/// name is carried whole ([`NodeId`], so a corner a cylinder made can be one) and the pin says
/// **which kind** it is ([`combinatorics::EndPin`]) — the same pair a segment's ends travel as.
/// The crossing arm writes the other pair — `Pierce`/`Cylinder` — where the
/// ring crosses the class line on a ruling ([`crossing_on_ruling`]).
struct Node {
    /// The point's own name — what the arrangement calls it, and what a segment's end records.
    id: NodeId,
    /// What pins it on `L`: the third plane class, or the cylinder whose crossing it is.
    pin: combinatorics::EndPin,
    /// Whether crossing this node toggles inside/outside. Strict crossings and single-vertex runs
    /// set it in Phase A; a two-vertex run's *hi* node has it set after the sort.
    flip: bool,
    /// A two-vertex on-line edge's run id; the two nodes share it, and the gap between them is
    /// forced-covered (it is boundary of `F`, hence of `F ∩ W`).
    run: Option<usize>,
    /// For run nodes: whether the flanks (the off-line neighbours) sit on opposite sides.
    flanks_differ: bool,
    /// A single on-line vertex with equal flanks — a tangential touch, not a crossing.
    single_touch: bool,
    /// `Some(body_above)` for **every** run — the side of `W` the face occupies along it, in the
    /// label frame, from [`run_body_above`]. A run is an edge of the face lying on `L`, so the face
    /// is on one side of it no matter what the ring does afterwards; `flanks_differ` is a parity
    /// fact, not an occupancy one, and must not gate this. `None` only for a strict crossing or a
    /// single-vertex touch, which are points rather than intervals.
    graze_above: Option<bool>,
}

/// A lateral face as its cycles state it: its axial range `[lo, hi]`, its whole rims `(station,
/// plane class)` in station order (0, 1 or 2 of them), and every other boundary cycle with its
/// kind.
struct LateralShape {
    range: [Rat; 2],
    rims: Vec<(Rat, usize)>,
    cycles: Vec<(combinatorics::CycleKind, combinatorics::LoopRing)>,
}

/// One angular extent of a lateral face's mark on a ⊥ plane class — see [`circle_on_class`].
///
/// ★ `arc` is `None` for the whole circle, and `Some([a, b])` runs **counter-clockwise** from `a`
/// to `b` about the cylinder's axis direction — [`MergedArc::end`]'s convention, because these are
/// what that becomes.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CircleSpan {
    pub arc: Option<[NodeId; 2]>,
    pub on: CylOnClass,
}

/// What one hole ring takes away from the circle: an extent, and what is left there.
struct Carved {
    /// CCW `[from, to]`.
    arc: [NodeId; 2],
    /// `None` — the face is **absent** over this extent (the class runs through the hole's
    /// interior). `Some(g)` — the extent is one of the hole's own rims, so the face grazes.
    on: Option<CylOnClass>,
}

/// Split the segments on plane class `wc` into the arrangement's 1-skeleton by a **per-wall 1D
/// interval overlay**: on each wall `W`, the line `wc ∩ W` is cut at the sorted union of every
/// W-segment's endpoints and every different-wall crossing, and one `MergedSeg` is emitted per
/// non-empty sub-interval carrying the **union** of the contributions of the W-segments whose
/// closed extent covers it. A crossing of two segments riding `W`/`o.wall` is the plane triple
/// `{wc, W, o.wall}` — no new point species.
///
/// This resolves same-wall partial overlap: where two segments overlap, the shared sub-interval
/// carries both their contributions, which the label brick reads per solid. Because sub-intervals
/// run between **distinct** points, no zero-length piece is ever emitted (a crossing coinciding with
/// an endpoint is the same plane class, deduped away). A clean arrangement — no same-wall overlap —
/// has exactly one segment covering each sub-interval, so the output is identical to the naive
/// per-segment split. No re-merge is needed: each `(wall, sub-interval)` is emitted once with its
/// full union, so no two outputs can coincide.
///
/// `Err` (honest reject) on: a degenerate endpoint name (`ThreePlanes`), or two distinct plane
/// classes coincident on a wall's line — a four-plane concurrency `{wc, W, a, b}` the 3-plane DCEL
/// cannot name (`FourPlane`).
/// One wall of a class's arrangement: which plane class it is, the **direction family** its line on
/// that class falls in, and the segments riding it.
///
/// Two walls' lines meet in a point iff their families **differ**, which is why the crossing
/// collector compares `dir` instead of asking a predicate per pair (see [`split_at_crossings`]).
///
/// ★ Everything a wall knows travels with it. The alternative — a `Vec<usize>` of classes beside a
/// `HashMap<class, Vec<usize>>` of segments beside a `Vec<usize>` of families indexed by *position* —
/// is two or three index spaces read in one breath, which is the shape this crate has been bitten by
/// four times (`combinatorics`'s "One `usize`, two meanings"). Same rule as `draft::Ring` carrying
/// the wall each edge rides rather than deriving it.
struct Wall {
    class: usize,
    dir: usize,
    segs: Vec<usize>,
}

/// **What one plane class's arrangement produced, before anything refuses it.**
///
/// ★★★ **Taking `emit_faces` in made this *wider*, and "a stage moving in shrinks the interface"
/// was never the law.** When `label_cells` moved in, `face_of` stopped crossing the boundary and
/// that read like a sign the boundary was right. It was not a general rule — that field simply had
/// one consumer. Here production reads only `faces` and `disk_labels`, but `frame_audit` is a
/// legitimate *second* consumer of `cells`, `nesting` and `labels` (the counts and the labels it
/// pins are exactly what the previous rung locked), so they stay. The test is not "did it shrink"
/// but **"does each field have a consumer"** — and the `cfg(test)` says which consumer, since that
/// second one is an instrument (`ClassAudit`, `#[cfg(test)]` like everything it touches).
struct Staged {
    #[cfg(test)]
    cells: Vec<Cell>,
    #[cfg(test)]
    nesting: Nesting,
    #[cfg(test)]
    labels: Vec<Label>,
    faces: Vec<LocalFace>,
    disk_labels: Vec<(usize, Label)>,
    arc_labels: Vec<(usize, ArcLabel)>,
    /// This class's ruling pieces, keyed by cylinder class — see [`RulingExtent`].
    ruling_extents: Vec<(usize, RulingExtent)>,
}

/// **One plane class's edges, after the arc split — and the only thing that knows the half-edge
/// numbering.**
///
/// Half-edge encoding: edge `i` gives `he = 2i` (forward, `end[0]→end[1]`) and `2i+1` (reverse);
/// `twin(he) = he ^ 1`. The ranges are, in order:
///
/// - **segments** `[0, 2·ns)` — an ordinary two-ended edge on a plane's meet with `P`;
/// - **arcs** `[2·ns, 2·(ns+na))` — a piece of a circle a segment cut, equally two-ended;
/// - **uncut circles**, on pseudo-half-edges `2·(ns+na) + 2i` — a closed curve with no vertex is
///   not an orbit the walk can express, so [`walk_cells`] appends those cells *after* the walk.
///
/// ★★★★ **The kind is a `match`, not `he >= 2 * segs.len()`.** There are three ranges, not two —
/// a cut circle's arcs sit between — so the arithmetic comparison reads an arc as a circle, and
/// indexing `segs[he / 2]` with it is not a wrong answer but an **index out of bounds**
/// (measured). One type owns the numbering, so the question is a `match` the compiler checks
/// rather than an arithmetic comparison each reader restates.
///
/// ★★ **It borrows or owns** (`Cow`). The split makes new `Vec`s; if this only borrowed, the
/// caller would have to keep them alive and the whole preamble — split, build, walk — would be
/// written once per caller, which is the shape `decline_to_reject`'s doc calls *"two copies … would
/// let them drift"*. Owning lets [`ClassEdges::of`] be the single constructor both callers use.
struct ClassEdges<'a> {
    segs: std::borrow::Cow<'a, [MergedSeg]>,
    arcs: std::borrow::Cow<'a, [MergedArc]>,
    /// Ruling pieces — a fourth range, empty until a tracer contributes
    /// rulings on a ∥ class.
    rulings: std::borrow::Cow<'a, [MergedRuling]>,
    circles: std::borrow::Cow<'a, [MergedCircle]>,
    /// Each cut circle's seam datum, `(cylinder class, rim)` — carried out to the assembly
    /// (keyed by plane class where the per-class products are aggregated).
    cut_rims: Vec<(usize, CutRim)>,
}

/// Which of the four ranges a half-edge is in, and its index within that range.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum HalfEdgeKind {
    Seg(usize),
    Arc(usize),
    Ruling(usize),
    Circle(usize),
}

/// The nesting of an arrangement's cells: which winding `-1` contour is the unbounded root, and
/// which are holes of a `+1` cell. A `+1` cell and its holes form one **group** (a face with
/// holes); every other `+1` cell and the root are singleton groups.
///
/// A contour `c` (winding `-1`) is a **hole** of a `+1` cell `r` iff `c` lies inside `r`, tested
/// only against an `r` it shares **no** node with. A shared node means `c` bounds `r` from the
/// other side (they are adjacent, not nested) — this also excludes `c`'s own `+1` partner, which
/// carries the same ring. Two vertex-disjoint simple loops are nested-or-separate (a crossing
/// would be a shared split node), so one representative vertex settles containment via
/// [`combinatorics::point_in_ring`]. A `c` inside no `+1` cell bounds the unbounded region: the root.
///
/// `c` may lie inside **several** `+1` rings at once, and that is not a degeneracy: nesting two
/// levels deep (`A⁺ ⊃ D⁻ ⊃ B⁺ ⊃ c⁻`) puts `c` inside `B`'s ring *and* `A`'s, because
/// [`combinatorics::point_in_ring`] asks about a ring, not about the material it bounds. The owner is the
/// **innermost** candidate ([`innermost_host`]).
struct Nesting {
    /// `group_of[cell]` = the cell's group representative (the `+1` host for a hole group, else
    /// the cell itself).
    group_of: Vec<usize>,
    /// The groups of the unbounded contours — `label_cells`' seeds. There can be **several**: one
    /// outside region may be bounded by more than one cycle (two boxes side by side on the plane
    /// give `0:A⁺ 1:A⁻ 2:B⁺ 3:B⁻`, and `1` and `3` bound the same outside). Every one of them is
    /// void, so every one seeds.
    root_groups: Vec<usize>,
    /// `holes[host]` = the `-1` cells nested in that `+1` host, emitted as its inner rings.
    holes: HashMap<usize, Vec<usize>>,
}

/// Drive the arrangement pipeline over **every** plane class and assemble the result solid — the
/// engine the public `crate::boolean` delegates to.
///
/// Vertex welding is automatic: every result vertex is a sorted triple of three canon plane
/// classes, so a corner shared by three planes gets one identical [`combinatorics::NodeId`]
/// regardless of which plane was the cut class W — `assemble_fuse_cut` welds them to one vertex.
/// A four-plane concurrency would name it inconsistently, which is what the alias table settles
/// and what `draft::Ring`'s carried walls keep out of the naming in the first place.
///
/// A class that declines (holes, degenerate) aborts the whole boolean: skipping it would drop real
/// faces and silently produce a non-manifold or wrong-volume solid.
pub(crate) fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(Vec<Handle<Solid>>, Notes, ClassOf), BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        geom,
        plane_ix,
        class_owner,
        cyls,
        crossings,
        tangencies,
        standard,
        notes,
        ..
    } = timed!(SETUP, plane_index_setup(model, a, b))?;
    // The operation's judging, made once: the dense plane table, the standard it is held to, and
    // the collector. Everything below reaches predicates through this, so there is exactly one
    // place where "how this boolean judges" is decided.
    let jd = Judge::new(&geom, standard, &notes);
    // The plane classes are already decided at this point — `plane_index_setup` runs
    // `Judge::planes_coplanar` to build them — so a judgement that could not be made has already
    // shaped everything downstream. Say so before doing the work it would invalidate.
    undecided_reject(&notes)?;
    // ★★★ **What every input face's plane became.** The classes are decided by now, and
    // `assemble_fuse_cut` gives each result face `geom[..].surf` — its class's representative — so
    // this is the only honest answer to *"which surface did my face's plane end up as?"*. A caller
    // that asks it afterwards, by comparing handles or coordinates, is re-deciding a question this
    // engine already settled with evidence (see `ops::find_face_coplanar_with`).
    // Plane rows only: "which plane class did my face end up on" is not a question a lateral
    // face has an answer to, and the report's consumers ask it of planar faces
    // (`ops::find_face_coplanar_with`).
    let class_of: ClassOf = surf_ix
        .iter()
        .filter_map(|(&f, &i)| match plane_ix[i] {
            ClassIx::Plane(c) => Some((f, geom[c].surf)),
            ClassIx::Cyl(_) => None,
        })
        .collect();
    // ★ **A closure so a `?` inside cannot skip the evidence check below.** Everything from here on
    // may reject, and every one of those rejects has to pass through `undecided_reject` first —
    // otherwise a symptom is reported where a precision failure is the cause. Returning early from
    // the function body would be exactly that bug.
    let run = |model: &mut Model| -> Result<Vec<Handle<Solid>>, BoolError> {
        let trace_in = timed!(
            TRACE_IN,
            combinatorics::trace_input(
                model,
                [(a, &inc_a), (b, &inc_b)],
                &surf_ix,
                faces_tab.len(),
                &jd,
                &plane_ix,
                &cyls,
                crossings.clone(),
            )
        );
        // `disk_labels`: per (cylinder class, plane class), the four bits of that circle's disk
        // cell — what the band pass reads instead of casting a witness ray.
        // `deferred`: an arc class's stopper reject, made per class and raised below **after** the
        // seam stretch, so the seam's pierce arm runs before the population is refused.
        let (faces, curved, deferred) = trace_result_faces(
            model,
            kind,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            &cyls,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &trace_in,
        )?;
        // **Every debug run answers the same boolean twice and requires the same answer.**
        //
        // Not tracing a class is the one thing here that changes what the engine *knows*: a class
        // never traced never reports the concurrencies it would have found, and 17% of the alias
        // table is discovered only in classes the other operand cannot reach. That those names are
        // wanted by nobody is an argument, so it is checked rather than believed — and at the
        // whole-boolean level, because per class there is no arrangement to compare to.
        //
        // The reference gets its own `Notes`: evidence is a side effect `undecided_reject` reads,
        // and a second run must not double it.
        #[cfg(debug_assertions)]
        {
            let plain_notes = Notes::new();
            let plain_jd = Judge::new(&geom, standard, &plain_notes);
            let plain_in = combinatorics::trace_input(
                model,
                [(a, &inc_a), (b, &inc_b)],
                &surf_ix,
                faces_tab.len(),
                &plain_jd,
                &plane_ix,
                &cyls,
                crossings.clone(),
            );
            let plain = trace_result_faces(
                model,
                kind,
                a,
                b,
                &plain_jd,
                &faces_tab,
                &plane_ix,
                &cyls,
                n_a,
                &class_owner,
                crate::reuse::ClassReuse::Off,
                &plain_in,
            );
            match &plain {
                // ★ **Faces only.** The disk labels are the same arrangement's product, so
                // comparing them would widen this differential's proposition ("reuse does not
                // change the faces") into one it was not built to make.
                Ok((p, _, _)) => assert_eq!(
                    crate::reuse::canonical(&faces),
                    crate::reuse::canonical(p),
                    "reuse changed the faces this boolean emits"
                ),
                // A reject only the reference reaches means reuse skipped the class that could not
                // be judged. Arguably better, definitely a change — so it fails here rather than
                // passing quietly.
                Err(e) => panic!("reuse turned a {e:?} into a result"),
            }
        }
        // ★★★ **The seam stretch, held as one closed result — the deferred stopper's second
        // interception layer.** Everything from the cleaning pass to the seam table runs even for
        // an arc input (that is the point: the seam's pierce arm is exercised), and then the
        // deferred reject wins over whatever the stretch produced, `Ok` *or* `Err`. Without the
        // `Err` half, an arc input would carry out whichever of the stretch's five fallible
        // steps — `unify_coplanar_faces`, `cyl_rows`, `emit_lateral`, the seam fill, the
        // `SeamAlias` scan — happened to fail first, and the population's name would depend on
        // how far the pipeline got: the same property the per-class `deferred.unwrap_or(e)` in
        // `arrange` protects, one level down. (`SeamAlias` is the sharp case: its class says
        // "report a bug", and its own doc records having mis-named a population once before.)
        let stretch = || -> Result<(Vec<LocalFace>, Vec<SeamVertex>), BoolError> {
            // Clean the raw arrangement output: merge coplanar, same-normal faces that share a full edge
            // (e.g. the split side walls a fused coincident interface leaves) so the result is a minimal,
            // chainable solid — a second boolean on it then sees no redundant coplanar planes.
            let faces = timed!(
                UNIFY,
                crate::assembly::unify_coplanar_faces(faces, &jd, &cyls)
            )?;

            // ★ **The lateral bands, appended after the differential above**: reuse can
            // only change what the *plane* arrangement emits, so the two routes are compared on that
            // list; the bands are a separate pass over the same operands and belong to neither route.
            // From here on there is one face list — the grouping, the closed-shell guard and the
            // assembly all read it.
            let faces = if cyls.is_empty() {
                faces
            } else {
                let rows = crate::bands::cyl_rows(&faces_tab, &plane_ix, n_a)?;
                // ★ **The lateral faces come from the chart**: every cylinder
                // class's cells, read off the plane arrangement's own labels and emitted in the
                // band road's vocabulary.
                let lateral = cyl_chart::emit_lateral(kind, &jd, &cyls, &faces, &curved, &rows);
                // ★★ The census sits *before* the curved cleaning pass on purpose: the cleaning
                // merges pieces, which would blur the face-by-face question being asked.
                // ★ The census holds the chart against its own rules (and the emitter against
                // the census's independent count), in test builds; the emitter's refusal is
                // handed to it before `?` decides, so a refused class is still recorded.
                #[cfg(test)]
                cyl_chart::census(&jd, &cyls, kind, &faces, &curved, &rows, &lateral);
                let mut faces = faces;
                faces.extend(lateral?);
                // ★ No curved cleaning pass follows: the emitter's lateral faces are the
                // regions of each chart already, so there is no phantom seam to erase (measured:
                // a cleaning pass here merges nothing over the suite).
                faces
            };

            let seam = seam_table(&faces, &cyls, &jd)?;
            Ok((faces, seam))
        };
        // ★ The raise moved into `reconstruct` (after the vertex naming, before any minting);
        // what stays here is the interception — a stretch failure on an arc input still carries
        // the stopper's name out, not its own.
        let (faces, seam) = match stretch() {
            Ok(x) => x,
            Err(e) => return Err(deferred.unwrap_or(e)),
        };

        timed!(
            ASSEMBLE,
            assemble_fuse_cut(
                model,
                a,
                b,
                &jd,
                &seam,
                &faces,
                &cyls,
                &curved.cut_rims,
                deferred,
                crate::assembly::Tangencies {
                    kind,
                    rows: &tangencies
                },
            )
        )
    };
    let out = run(model);
    // **The cause outranks the symptom, on both paths.** An undecided judgement has already been
    // read as a `0` by everything downstream, so whatever the engine then complains about — a
    // loop that will not orient, a trace that will not close — is a consequence being reported as
    // if it were the problem — a precision shortage hiding behind a ring-orientation symptom — and
    // checking the evidence *before* returning the symptom is what keeps it shut.
    undecided_reject(&notes)?;
    out.map(|solids| (solids, notes, class_of))
}

/// **A judgement that could not be made is a reject, not a zero.**
///
/// `Decision::orient` collapses every inconclusive outcome to `Zero`, which the arrangement reads
/// as "these are the same thing" — so an undecided judgement that reaches the geometry has
/// already merged something the kernel never established. The two causes get their own names,
/// because they call for opposite responses: one is a budget (raise it, or simplify the model),
/// the other is the arrangement itself (nothing to raise).
fn undecided_reject(notes: &Notes) -> Result<(), BoolError> {
    for e in notes.sorted() {
        match e.outcome {
            Decision::Exhausted { .. } => return Err(reject(RejectReason::JudgeExhausted)),
            Decision::Degenerate => return Err(reject(RejectReason::DegenerateWitness)),
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "../tests/arrangement/mod.rs"]
mod tests;
