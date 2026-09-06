//! The plane-class arrangement engine — the sole boolean path ([`boolean`], which `crate::boolean`
//! delegates to). It traces each face of both solids onto every plane class as **line segments**
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
use crate::boolean::*;
use crate::combinatorics::{Canon3, NodeId, three_plane_name};
use crate::planes::*;
use crate::tolerant::Judge;
#[cfg(test)]
use crate::transform::transform;
use nacre_cip::Decision;
use nacre_cip::predicate::Notes;
use nacre_geom::intersect::three_planes;

/// **Phase timers, on the production path.**
///
/// ★ They live *inside* `boolean` rather than in a spike that replays it, because a replica measures
/// the proposition next to the one that matters — the last profile of these phases was taken with
/// `ClassReuse::Off` and read as if it were production's.
///
/// Read them with `--no-default-features`: parallel accumulation would sum CPU across threads, and a
/// share of a wall-clock whole computed from a CPU sum is not a share of anything. (That mistake
/// once made a part measure larger than its whole.) `cfg(test)` so release binaries carry nothing.
#[cfg(test)]
pub(crate) mod phase {
    use std::sync::atomic::{AtomicU64, Ordering};

    macro_rules! counters {
        ($($name:ident = $label:literal),+ $(,)?) => {
            $(pub(crate) static $name: AtomicU64 = AtomicU64::new(0);)+
            /// Every counter with its label, in report order.
            pub(crate) fn all() -> Vec<(&'static str, u64)> {
                vec![$(($label, $name.load(Ordering::Relaxed))),+]
            }
            pub(crate) fn reset() { $($name.store(0, Ordering::Relaxed);)+ }
        };
    }

    counters! {
        SETUP      = "plane_index_setup",
        S_TRIPT3   = "    collect: tri_pt3 (chain replay)",
        S_COLLECT  = "    collect: the per-face loop",
        S_STD      = "    standard_for",
        S_EDGES    = "    edge_faces x2",
        S_CLASSES  = "    plane_classes (pairwise)",
        S_DENSE    = "    dense_planes + owners",
        TRACE_IN   = "trace_input (face table)",
        TRACE_ON   = "  trace_on_class",
        MERGE      = "  merge_coincident",
        SPLIT      = "  split_at_crossings",
        S_PART     = "    build the direction partition",
        COLLECT    = "    (1) collect crossings",
        SORT       = "    (2) sort + flush groups",
        COVER      = "    (3) cover sub-intervals",
        CELLS      = "  split + cells + nest + label + emit",
        C_SPLIT    = "    split_circles (ClassEdges::of)",
        C_EXTRACT  = "    walk_cells",
        E_ORDER    = "      per-half-edge order_along",
        E_ANGULAR  = "      angular_order",
        E_WALK     = "      the rest (walk + cells)",
        C_NEST     = "    nest_cells",
        C_LABEL    = "    label_cells",
        C_EMIT     = "    emit_faces",
        REUSE      = "  reuse pass-through",
        UNIFY      = "unify_coplanar_faces",
        SEAM       = "seam table + alias scan",
        ASSEMBLE   = "assemble_fuse_cut",
    }

    /// **Scale, not time.** The two ratios that decide whether S3b/S3c are worth building: how many
    /// segments a wall carries (a hull test over wall *pairs* buys nothing at 1), and how big the
    /// covering loop is against the collecting one.
    pub(crate) mod scale {
        use std::sync::atomic::{AtomicU64, Ordering};
        pub(crate) static SEGS: AtomicU64 = AtomicU64::new(0);
        pub(crate) static WALLS: AtomicU64 = AtomicU64::new(0);
        pub(crate) static PTS: AtomicU64 = AtomicU64::new(0);
        /// `Σ_w |pts on w| × |segs on w|` — the covering loop's actual trip count.
        pub(crate) static COVER_TRIPS: AtomicU64 = AtomicU64::new(0);
        /// `Σ_w |segs|` — the collecting loop's actual trip count.
        pub(crate) static COLLECT_TRIPS: AtomicU64 = AtomicU64::new(0);
        pub(crate) fn add(c: &AtomicU64, n: usize) {
            c.fetch_add(n as u64, Ordering::Relaxed);
        }
        pub(crate) fn get(c: &AtomicU64) -> u64 {
            c.load(Ordering::Relaxed)
        }
        pub(crate) fn reset() {
            for c in [&SEGS, &WALLS, &PTS, &COVER_TRIPS, &COLLECT_TRIPS] {
                c.store(0, Ordering::Relaxed);
            }
        }
    }

    /// Run `f`, adding its wall time to `c`. Returns what `f` returned.
    pub(crate) fn timed<R>(c: &AtomicU64, f: impl FnOnce() -> R) -> R {
        let t = std::time::Instant::now();
        let r = f();
        c.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        r
    }

    /// Charges the enclosing scope to `c` **on drop** — for a block whose `?` or early `return`
    /// would jump past the end of a closure.
    pub(crate) struct Watch(&'static AtomicU64, std::time::Instant);

    impl Watch {
        pub(crate) fn new(c: &'static AtomicU64) -> Self {
            Self(c, std::time::Instant::now())
        }
    }

    impl Drop for Watch {
        fn drop(&mut self) {
            self.0
                .fetch_add(self.1.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
    }
}

/// `phase::timed` where the counters exist, and a plain call where they do not.
macro_rules! timed {
    ($c:ident, $e:expr) => {{
        #[cfg(test)]
        {
            phase::timed(&phase::$c, || $e)
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
        let _w = phase::Watch::new(&phase::$c);
    };
}

/// **Which arcs of the circle around an edge a wall fills**, over one sub-interval of that edge.
///
/// An arrangement edge on plane class `W` is the meeting of `W` with one other plane, so the little
/// circle around it is cut into just two arcs — above `W` and below. A contributing face either
/// fills both (it passes through `W` there) or exactly one (its own boundary lies along the edge
/// there, and it hangs off to one side). That is the whole content of this enum, and `edge_mask`
/// turns it into the label flip. It is the BRep spelling of a Nef local pyramid restricted to two
/// planes, and the flip itself is binary winding propagation (Zhou et al. 2016).
///
/// ★ **The occupancy is a fact about the sub-interval, not about the face.** One face can fill both
/// arcs over one stretch of `W ∩ fp` and one arc over the next — a notched or holed face does
/// exactly that where its boundary rides the line. So a segment must be *homogeneous* in this kind:
/// [`split_at_crossings`] subdivides segments later and every piece inherits the kind verbatim, so
/// a mixed segment silently mislabels the pieces that disagree with it. Phase C therefore closes and
/// reopens wherever the kind changes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SegKind {
    /// `L` runs through the face's **interior** here, so the face is on both sides of `W` and
    /// crossing this segment flips the above-label and the below-label **together**. `mat` is which
    /// side of the directed line `W ∩ fp` the solid's material lies on (`+1`/`-1`, relative to
    /// `d = n_W × n_fp`).
    Transversal { mat: i8 },
    /// A boundary edge of a face lying **in** `W`. The face is the boundary between material and
    /// void on one side of `W`, so crossing it flips **one** label — the side the body occupies
    /// (`body_above`, from the face's own outward normal against `W`'s stored normal).
    Seated { body_above: bool },
    /// `L` runs along the face's **boundary** here, so the face fills one arc only and crossing this
    /// segment flips **one** label — `body_above`, the side it fills, from [`run_body_above`].
    ///
    /// Not "the face merely touches `W`": a face may pass clean through `W` elsewhere and still fill
    /// one arc over *this* stretch, which is what the name `Graze` originally missed. Kept distinct
    /// from `Seated` so `edge_mask` can let it override a coincident seated rim at a reflex
    /// dihedral, where the two disagree on which side flips.
    Graze { body_above: bool },
    /// The **tangent ruling** of a face lying in `W` (cell ⑩): the face ends where its cylinder
    /// begins, tangent to `W`. Beyond the line what `W` separates is decided by which side the
    /// cylinder bends toward — its **axis** side: a convex fillet's axis is on the body side and
    /// that side turns void past the line (the cusp under the arc); a hole's axis is on the void
    /// side and that side turns material. So crossing this edge flips the axis side's label, one
    /// bit, whichever side the body was on.
    Tangent { axis_above: bool },
}

/// One segment of a solid's trace on a plane class, named entirely in plane-class triples.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Seg {
    /// Canon plane class the segment rides besides the cut plane (its line is `W ∩ wall`).
    pub wall: usize,
    /// Endpoint identities, canonized.
    pub end: [NodeId; 2],
    /// The same two endpoints as **handles on this segment's own line**: the third plane pinning
    /// each on `wc ∩ wall`. Carried, not recovered from `end`, because a canonical name need not
    /// mention either of this line's planes.
    pub end_h: [combinatorics::EndPin; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// A **full circle** of a solid's trace on a plane class (M6-2a): a cylinder's mark, closed —
/// no endpoints, no wall, no place in the segment machinery. Seated circles come from disk
/// faces and circular holes lying in the class; transversal circles from a lateral surface
/// crossing it, and circles join the arrangement only at the cell stage.
///
/// ★★ **That it is *full* is checked, not inherited.** The population gate's wall rule used to
/// keep this as a side effect — every ∥ wall face stands clear of the lateral, and a segment on
/// this class is that face's own trace, so it inherited the clearance — which made a promise about
/// *circles* rest on a rule about *walls*. [`split_circles`] asks it of the segments
/// themselves now, so a crossed circle is **split into arcs** rather than treated as the
/// closed cell it is not, and the wall rule is free to become precise about its own question.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CircleTrace {
    /// The cylinder class ([`ClassIx::Cyl`] payload) whose circle this is.
    pub cyl: usize,
    pub solid: SolidSide,
    pub kind: SegKind,
    /// **The angular extent this contribution covers**, counter-clockwise `[from, to]` about the
    /// cylinder's axis — `None` for the whole circle.
    ///
    /// ★★ A lateral face with a **hole** contributes over part of its circle and nothing over the
    /// rest, and a rim of that hole contributes a *different* kind from the band around it. The
    /// straight sibling [`RulingTrace`] has carried ends since it was written, for the same
    /// reason; a circle carried none because before holes there was no partial circle to state.
    /// [`split_circles`] cuts at these ends and hands each arc only the contributions that cover
    /// it.
    pub arc: Option<[NodeId; 2]>,
}

/// One circle of the class's arrangement after merging: both operands' contributions on one
/// cylinder class, plus the cylinder's exact statement (the containment predicates read it).
/// The circle twin of [`MergedSeg`].
#[derive(Clone, Debug)]
pub(crate) struct MergedCircle {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    /// Each contribution with the **angular extent** it covers — see [`CircleTrace::arc`]. `None`
    /// is the whole circle, which is what every contribution was before a lateral face could have
    /// a hole.
    pub merged: Vec<(SolidSide, SegKind, Option<[NodeId; 2]>)>,
}

impl MergedCircle {
    /// The contributions as `edge_mask` reads them.
    ///
    /// ★ A circle only reaches a mask **uncut**, and a partial contribution forces its own cut
    /// ([`split_circles`] feeds every extent end into the split), so every extent here is `None`.
    /// A `Some` would mean a whole-circle mask was about to be built from a partial trace — the
    /// silently wrong answer this vessel exists to prevent — so it is refused by name.
    fn whole_marks(&self) -> Result<Vec<(SolidSide, SegKind)>, BoolError> {
        self.merged
            .iter()
            .map(|&(s, k, arc)| match arc {
                None => Ok((s, k)),
                Some(_) => Err(reject(RejectReason::PartialCircleUncut)),
            })
            .collect()
    }
}

/// **One arc of a circle a segment cut** — the DCEL edge a crossed circle becomes.
///
/// ★★ **A circle leaves the parallel road here.** An uncut circle is a one-edge closed cell that
/// the orbit walk cannot express (orbits need ≥ 2 half-edges and a vertex to turn at), which is why
/// [`walk_cells`] appends it *after* the walk on pseudo-half-edges. Cut, it has endpoints and is
/// an ordinary edge — so the machinery it used to stand beside is the machinery it now uses.
///
/// ★ `end` is in the arc's **own travel order**, `end[0]` → `end[1]` counter-clockwise about the
/// circle's normal. That is not a convention this crate invents: STEP's `EDGE_CURVE` runs from
/// `edge_start` to `edge_end` along increasing parameter, and `nacre-step` already writes a
/// circular edge that way — so the exporter needs nothing, and the STEP round-trip measures it.
#[derive(Clone, Debug)]
pub(crate) struct MergedArc {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    pub end: [NodeId; 2],
    /// The contributions that **cover this arc** — no longer the circle's list entire.
    ///
    /// ★★ It used to read "inherited whole from the circle: an arc is a piece of the same trace,
    /// so it carries the same contributions". That was true while a lateral face could only mark a
    /// class over its whole circle, and a face with a **hole** makes it false: over the hole the
    /// face is not there at all, and along the hole's rim it grazes where the band around it
    /// crosses. [`split_circles`] selects per arc; `edge_mask` still reads the result unchanged.
    ///
    /// ★ Read by `label_cells`' `mask_of`: crossing an arc flips the same bits crossing its circle
    /// would, which is what "a piece of the same trace" means. (It was carried before that consumer
    /// existed, because the split is the only place that knows which circle an arc came from.)
    pub merged: Vec<(SolidSide, SegKind)>,
}

/// **One ruling piece of the class's arrangement** — the straight sibling of [`MergedArc`], for a
/// class **parallel** to a cylinder's axis (the M6-2 rulings ladder): the wall plane meets the
/// lateral surface in up to two axis-parallel lines, and a piece of one is an ordinary two-ended
/// edge whose carrier is the cylinder, not a plane pair.
///
/// ★ `end` is in the piece's **own travel order**: `end[0]` → `end[1]` runs along `+m` (the
/// cylinder's axis direction) — the straight reading of [`MergedArc`]'s CCW convention, read the
/// same way (`edge_at`: the even half-edge travels the stated way, its twin the other).
///
/// `side` names which of the two parallel rulings ([`combinatorics::RulingCarrier::side`] — the
/// sign of `(x − o) · (m × n̂)` against the class's canonical coefficients).
#[derive(Clone, Debug)]
pub(crate) struct MergedRuling {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    pub side: i8,
    pub end: [NodeId; 2],
    /// The `(solid, kind)` contributions, like a segment's — `edge_mask` reads it unchanged.
    pub merged: Vec<(SolidSide, SegKind)>,
    /// The lateral's `orient_sign` (`+1` material inside the cylinder, `−1` outside — a boss or a
    /// bore), for the vertical answer's content check; an instrument's field, like the label.
    #[cfg(test)]
    pub orient: i8,
}

/// Group a class's circle traces by cylinder class, in ascending class order (deterministic —
/// replay mints handles from this order). The def comes from the class table, which a
/// `ClassIx::Cyl` index indexes directly.
fn merge_circles(
    circles: &[CircleTrace],
    cyls: &[crate::planes::WorkingCyl],
    aliases: &Aliases,
) -> Result<Vec<MergedCircle>, BoolError> {
    // ★ Cell ⑫: the merge layer is where names become canonical — `merge_coincident` for a
    // segment's ends, here for an extent's, `merge_rulings` for a ruling's — so everything after
    // it (splits, the walk, the chart, assembly) sees one name per point.
    let canon = |e: [NodeId; 2]| [aliases.canon_point(e[0]), aliases.canon_point(e[1])];
    let mut out: Vec<MergedCircle> = Vec::new();
    let mut sorted: Vec<&CircleTrace> = circles.iter().collect();
    sorted.sort_by_key(|c| c.cyl);
    for c in sorted {
        if let Some(last) = out.last_mut() {
            if last.cyl == c.cyl {
                last.merged.push((c.solid, c.kind, c.arc.map(canon)));
                continue;
            }
        }
        // ★ The class table answers directly. This used to scan the face rows for a lateral face
        // of the same class and `expect` one — a runtime claim where the truth is structural: a
        // `ClassIx::Cyl` index is an index *into this table*, built from the same map that issued
        // it ("rebuilt from the same map so the two cannot drift", `dense_planes`).
        let def = cyls[c.cyl].def.clone();
        out.push(MergedCircle {
            cyl: c.cyl,
            def,
            merged: vec![(c.solid, c.kind, c.arc.map(canon))],
        });
    }
    Ok(out)
}

/// **One piece of a solid's ruling trace on a ∥ class** (M6-2 rulings ladder): a lateral face's
/// mark on a plane class that contains the cylinder's axis, along one of the two axis-parallel
/// lines. The ruling sibling of [`CircleTrace`], with ends because a ruling is not closed:
/// `end[0]` → `end[1]` ascends the axis, and both are Branch names (`{wc, ⊥ plane, cyl, root}` —
/// real classes, found by the tracer).
///
/// ★★ **It used to say "spanning the face's own rims", and a hole makes that false.** Where the
/// face is buried in the other body it does not *cross* the wall its hole's vertical edges lie on
/// — it ends at it — so one ruling comes in pieces of different [`SegKind`]s and each is its own
/// trace. See [`ruling_sweep`]. A face with no hole still yields exactly one piece per ruling,
/// rim to rim.
#[derive(Clone, Debug)]
pub(crate) struct RulingTrace {
    /// The cylinder class, with `side` the ruling's identity
    /// ([`combinatorics::RulingCarrier::side`]).
    pub cyl: usize,
    pub side: i8,
    pub end: [NodeId; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
    /// The lateral's `orient_sign`, carried to [`MergedRuling::orient`].
    #[cfg(test)]
    pub orient: i8,
}

/// Group the class's ruling traces into [`MergedRuling`]s — the ruling twin of
/// [`merge_circles`]: dedupe identical rulings (two operands stating one), collect their
/// contributions, def from the class table. Deterministic order: `(cyl, side descending, extent)` —
/// a fixed rule, ascending cyl with the `+1` ruling first.
///
/// ★★★★★ **The extent is part of the key, and leaving it out is a silent duplicate waiting for a
/// population.** The fold below is *adjacency*-based: it merges a trace into the previous one when
/// the two name the same edge. While each `(cyl, side)` carries exactly one trace per operand, the
/// sort puts those two side by side and the fold works. The moment one ruling comes in **pieces** —
/// a lateral face with a hole grazes part of its own ruling instead of crossing it — a stable sort
/// on `(cyl, side)` alone leaves operand A's pieces before operand B's, so the two statements of one
/// piece are **not adjacent** and each becomes its own `MergedRuling`: one edge appearing twice,
/// with half its contributions each. Sorting by the extent too puts identical pieces together
/// whatever order the tracer emitted them in.
///
/// ☑ Measured over the suite: **no `(cyl, side)` group carries more than one trace today**, so the
/// fold below never actually fires and this key change moves nothing. Both facts are written for
/// the population that has not arrived yet — a cylinder belongs to one operand, so two statements
/// of one ruling need a shape nothing builds so far.
fn merge_rulings(
    rulings: &[RulingTrace],
    cyls: &[crate::planes::WorkingCyl],
    aliases: &Aliases,
) -> Vec<MergedRuling> {
    let mut out: Vec<MergedRuling> = Vec::new();
    // Canonical ends (cell ⑫, see `merge_circles`) — the key below is by name.
    let mut sorted: Vec<(&RulingTrace, [NodeId; 2])> = rulings
        .iter()
        .map(|r| {
            (
                r,
                [aliases.canon_point(r.end[0]), aliases.canon_point(r.end[1])],
            )
        })
        .collect();
    sorted.sort_by_key(|(r, end)| (r.cyl, -r.side, *end));
    for (r, end) in sorted {
        if let Some(last) = out.last_mut() {
            if last.cyl == r.cyl && last.side == r.side && last.end == end {
                last.merged.push((r.solid, r.kind));
                continue;
            }
        }
        out.push(MergedRuling {
            cyl: r.cyl,
            def: cyls[r.cyl].def.clone(),
            side: r.side,
            end,
            merged: vec![(r.solid, r.kind)],
            #[cfg(test)]
            orient: r.orient,
        });
    }
    out
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// The closed circle elements beside the segments (M6-2a) — see [`CircleTrace`].
    pub circles: Vec<CircleTrace>,
    /// The rulings beside them (M6-2 rulings ladder) — see [`RulingTrace`].
    pub rulings: Vec<RulingTrace>,
    /// Names this trace found to denote one feature — see [`Aliases`].
    pub aliases: Aliases,
    /// Single-point tangential contacts — real arrangement vertices, but not segments (a
    /// zero-length chord would abort at `Line::through_points`).
    pub touches: Vec<NodeId>,
    /// `(face index, reason)` for every face this brick could not trace.
    pub declined: Vec<(usize, DeclineKind)>,
}

/// Sort a triple and reject it if two of its planes coincide.
///
/// The triples already carry dense plane ids — `loop_triples` maps face indices through `plane_ix`
/// at the source — so there is nothing to canonize here; the point is the collapse check. `None`
/// when two names are equal: such a triple defines no point (`three_planes` answers `None`, and the
/// exact predicates need `D ≠ 0`), so the face is declined rather than fed a degenerate meet. A
/// face→plane map still collapses — two faces of one solid meeting a vertex on one plane — so this
/// guard outlives the old face/plane ambiguity it was born with.
/// The reject a declined face reports to the caller.
///
/// Most kinds *are* the answer — the tracer says what it could not do, and `face` says where.
/// [`DeclineKind::FourPlane`] is the exception: the naming failure is a symptom, and reporting it
/// as one would hide the substrate limit that caused it, so it is raised as the cause. The face
/// handle is dropped there because [`RejectReason::FourPlane`] carries no payload; per-face detail
/// remains in the class audit, which is where `TraceDeclined`'s own docs put it.
///
/// **One function for both consumers.** The boolean's error and the audit's `failed_at` must agree
/// — two copies of this mapping would let them drift.
fn decline_to_reject(kind: DeclineKind, face: Option<Handle<Face>>) -> RejectReason {
    match kind {
        DeclineKind::FourPlane => RejectReason::FourPlane,
        kind => RejectReason::TraceDeclined { kind, face },
    }
}

/// **The alias table learns from the operands first** (cell ⑪): every concurrency an operand's
/// own topology knows — a vertex with four or more incident plane classes, carried by its rings as
/// `NamedRing::concurrencies` — is recorded before any class is traced. The representative the
/// ring already named the vertex by (`canonical_triple`) is then the representative the
/// arrangement's own discoveries (`{wc} ∪ t` at a run vertex, a wall family's line) fold onto.
fn seed_from_operands(
    aliases: &mut Aliases,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    trace_in: &combinatorics::TraceInput,
) {
    let mut corners: Vec<NodeId> = Vec::new();
    for side in &trace_in.faces {
        for (_, loops) in side {
            let rings = loops
                .outer
                .iter()
                .chain(loops.holes.iter().flatten())
                .chain(loops.cycles.iter().flatten().map(|(_, r)| r));
            for lr in rings {
                if let Some(nr) = lr.poly() {
                    for s in &nr.concurrencies {
                        aliases.record(jd, s);
                    }
                    corners.extend(
                        nr.triples
                            .iter()
                            .copied()
                            .filter(|&n| combinatorics::branch_name(n).is_some()),
                    );
                }
            }
        }
    }
    // ★ Cell ⑫ — **the operands' cylinder corners, against every class.** A branch corner
    // (`Branch{[p0, p1], cyl, root}` — a fillet's tangent corner, a boss's foot) is a point some
    // *other* class may pass through: the gusset's side plane through the fillet axis contains
    // the tangent ruling, and so the corner. Then the point has three names — its own, the
    // three-plane `[p0, p1, wc]`, and the class's crossing of that ruling — and
    // [`Aliases::record_on_cylinder`] joins them. Asked here, once, of the operands' own
    // topology and the class table, so every round and every class trace starts knowing it:
    // asked during a trace instead (where `third_on_l` meets the corner on `wc`), a class traced
    // in the same round could refuse the coincidence before the round that learnt it — a round
    // that declines returns, and there is no next. `side_of` is exact (the quad tower), and the
    // question is corners × classes, both small.
    corners.sort_unstable();
    corners.dedup();
    for &corner in &corners {
        let Some((planes, _, _)) = combinatorics::branch_name(corner) else {
            continue;
        };
        for c in 0..jd.planes.len() {
            if planes.contains(&c) {
                continue;
            }
            if combinatorics::side_of(jd, cyls, corner, c) == Some(0) {
                aliases.record_on_cylinder(jd, cyls, corner, c);
            }
        }
    }
}

/// The names that turned out to denote **one** feature, learned while tracing.
///
/// Both kinds come from the same discovery. When a producer learns the full set `S` of planes
/// through a point:
///
/// - **A point** with `|S| > 3` has one valid name per 3-subset, so the tables would enter it once
///   per name. [`Aliases::point`] folds them onto the lexicographically first subset that actually
///   names a point (three planes sharing a line name none, so such a subset cannot be the winner).
/// - **A line** shows up as a 3-subset of `S` whose planes share a line rather than meeting at the
///   point. On each of those three classes the other two are *the same wall* — the engine calls a
///   line by the plane it rides besides the cut plane, so that is two names for one line, and
///   [`Aliases::wall`] folds them onto the smallest.
///
/// Neither fold is optional if the other happens: `merge_coincident` keys an edge by
/// `(wall, endpoints)`, so a duplicate edge merges only when **both** its wall and its endpoint
/// names agree. That is why the two live in one table and are applied together.
///
/// ★ Cell ⑪: the table has **two sources**. The operands seed it before any class is traced
/// ([`seed_from_operands`] — every vertex with four or more incident plane classes, which the
/// operand's own topology knows), and the tracer adds what it discovers (`{wc} ∪ t` at a run
/// vertex, a wall family's line). The representative is [`combinatorics::canonical_triple`]'s
/// answer, which is also the name the operand's ring already gave the point; with four planes
/// through a point every record is the whole set, so the two sources land in one component.
// `Clone` so a round can hand every class the table as it stood when the round began, and
// merge their discoveries afterwards — see `trace_result_faces`. In the ordinary model the
// maps are empty, so the copy costs nothing.
#[derive(Default, Debug, Clone)]
pub(crate) struct Aliases {
    /// Union-find over vertex names.
    point: HashMap<NodeId, NodeId>,
    /// Union-find over `(class, wall)` — walls of one class that carry the same line.
    wall: HashMap<(usize, usize), usize>,
}

impl Aliases {
    /// Record that every plane in `s` (sorted, deduped) passes through one point.
    fn record(&mut self, jd: &Judge<'_, WorkingPlane>, s: &[usize]) {
        if s.len() < 4 {
            return; // three planes meeting at a point is the ordinary case, and names nothing new
        }
        // ★ Cell ⑪: the representative is [`combinatorics::canonical_triple`]'s answer — the one
        // rule every producer of a point's name calls — so a name an operand's ring already gave
        // the point is the representative it folds onto here.
        let rep = combinatorics::canonical_triple(jd, s).map(NodeId::three_planes);
        for i in 0..s.len() {
            for j in (i + 1)..s.len() {
                for k in (j + 1)..s.len() {
                    let t = [s[i], s[j], s[k]];
                    if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                        // Shares a line: names no point, and tells us two walls are one line.
                        self.union_wall(t[0], t[1], t[2]);
                        self.union_wall(t[1], t[0], t[2]);
                        self.union_wall(t[2], t[0], t[1]);
                    } else if let Some(rep) = rep {
                        self.union_point(rep, NodeId::three_planes(Canon3::three(t)));
                    }
                }
            }
        }
    }

    /// **Two names for one point, learned by a producer that knows both** (cell ⑫) — the door
    /// [`Aliases::record_on_cylinder`] and the tests use; [`Aliases::record`] is the plane-set
    /// spelling of the same fold. The representative is the union-find's minimum, and it is a
    /// **key**: consumers that need a name's geometry (a branch's meet, its θ) keep asking with
    /// the name they hold.
    /// **A branch corner found on a further plane class** (cell ⑫) — the cylinder twin of
    /// [`Aliases::record`]. The corner's planes `p0, p1` and the class `wc` all pass through the
    /// point, and so does the cylinder; every name that set can produce denotes it: the corner's
    /// own, the three-plane name of `{p0, p1, wc}` (when independent), and `wc`'s ruling crossing
    /// with the corner's cap at the root that lies on the corner's side of `wc`. All of them are
    /// folded here, in one place, so the ruling sweep and the plane roads — which mint the names
    /// on their own — never have to decide "same point" themselves.
    ///
    /// The crossing's side is the corner's own: [`ruling_side`] of the corner's meet against
    /// `wc`'s stored normal, the predicate the sweep uses to tell `wc`'s two rulings apart.
    pub(crate) fn record_on_cylinder(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        corner: NodeId,
        wc: usize,
    ) {
        let Some((planes, cyl, _)) = combinatorics::branch_name(corner) else {
            return;
        };
        let mut s = vec![planes[0], planes[1], wc];
        s.sort_unstable();
        s.dedup();
        if s.len() < 3 {
            return; // `wc` is one of the corner's own planes: the ordinary corner, nothing to fold
        }
        if let Some(t) = combinatorics::canonical_triple(jd, &s) {
            self.union_point(corner, NodeId::three_planes(t));
        }
        let Some(wcy) = cyls.get(cyl) else { return };
        let def = &wcy.def;
        let Some(w) = combinatorics::class_coeffs_rat(jd, wc) else {
            return;
        };
        let Some(meet) = combinatorics::branch_meet(jd, cyl, def, corner) else {
            return;
        };
        let Some(side) = ruling_side(&w, def, (&meet.0, &meet.1)) else {
            return; // the corner sits on `wc`'s axis plane's own line: no ruling to name
        };
        for fc in planes {
            if let Ok(id) = crossing_on_ruling(jd, def, fc, wc, cyl, side) {
                self.union_point(corner, id);
            }
        }
    }

    fn find_point(&self, t: NodeId) -> NodeId {
        let mut x = t;
        while let Some(&p) = self.point.get(&x) {
            if p == x {
                break;
            }
            x = p;
        }
        x
    }

    /// **The representative of a point's names** — the one every key uses (`canon_point`).
    ///
    /// Among three-plane names it is the least, which is [`combinatorics::canonical_triple`]'s
    /// own answer (cell ⑪: the operand's name and the arrangement's discoveries meet there).
    /// ★ Cell ⑫: **a point on a cylinder is represented on the cylinder.** A `Branch` name
    /// locates the point exactly *and* says which cylinder it lies on — the lateral chart's
    /// geometry (θ about the axis, `branch_meet`) reads that from the name — while a three-plane
    /// name of the same point (a class through a tangent corner: `ThreePlane([cap, t, wc])`) is
    /// a key only. So a branch name outranks a three-plane name; among branch names the least.
    /// Either way the representative is a function of the class alone.
    fn rep_rank(n: NodeId) -> (u8, NodeId) {
        (u8::from(matches!(n, NodeId::ThreePlane(_))), n)
    }

    fn union_point(&mut self, a: NodeId, b: NodeId) {
        let (ra, rb) = (self.find_point(a), self.find_point(b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if Self::rep_rank(ra) < Self::rep_rank(rb) {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.point.insert(hi, lo);
        self.point.entry(lo).or_insert(lo);
    }

    fn find_wall(&self, class: usize, w: usize) -> usize {
        let mut x = w;
        while let Some(&p) = self.wall.get(&(class, x)) {
            if p == x {
                break;
            }
            x = p;
        }
        x
    }

    fn union_wall(&mut self, class: usize, a: usize, b: usize) {
        let (ra, rb) = (self.find_wall(class, a), self.find_wall(class, b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.wall.insert((class, hi), lo);
        self.wall.entry((class, lo)).or_insert(lo);
    }

    /// The identity to key a vertex on — itself when nothing was merged with it.
    pub(crate) fn canon_point(&self, t: NodeId) -> NodeId {
        if self.point.is_empty() {
            return t; // the overwhelmingly common case pays nothing
        }
        self.find_point(t)
    }

    /// The name to call a line by on `class` — `w` itself when no other wall carries it.
    pub(crate) fn canon_wall(&self, class: usize, w: usize) -> usize {
        if self.wall.is_empty() {
            return w;
        }
        self.find_wall(class, w)
    }

    /// The walls that carry **one line** with `w` on `class` — `w` alone when none does.
    ///
    /// The wall fold says two classes name one line; this reads that back out, because a point on
    /// that line lies on **every** plane in the family, and the point fold needs to hear about it
    /// (see the report in [`split_at_crossings`]).
    fn wall_family(&self, class: usize, w: usize) -> Vec<usize> {
        if self.wall.is_empty() {
            return vec![w];
        }
        let rep = self.find_wall(class, w);
        let mut fam: Vec<usize> = self
            .wall
            .keys()
            .filter(|(c, _)| *c == class)
            .map(|&(_, x)| x)
            .filter(|&x| self.find_wall(class, x) == rep)
            .collect();
        fam.push(w);
        fam.sort_unstable();
        fam.dedup();
        fam
    }

    fn absorb(&mut self, other: &Aliases) {
        for (&k, &v) in &other.point {
            self.union_point(k, v);
        }
        for (&(c, w), &v) in &other.wall {
            self.union_wall(c, w, v);
        }
    }

    /// How much has been learned — the fixed-point loop's measure.
    fn len(&self) -> usize {
        self.point.len() + self.wall.len()
    }
}

/// Why a ring of names could not be carried.
///
/// ★★★★★ **It had three variants and now has one, and that is the shape of what changed.** The
/// other two — a corner with no three-plane name, a carrier riding a cylinder — were never
/// *decisions*: the ring came out as plane ids and a curved ring could not be **described** in it,
/// so the type ran out and the refusal was dressed up as a cause. Two rungs widened the vessel and
/// both refusals went with the projections they were standing in for.
///
/// What is left is a name that is degenerate **as a name**, which is a fact about the triple and
/// has nothing to do with cylinders — so it survives, alone, and says only itself.
#[derive(Debug)]
enum RingFail {
    /// Two of a vertex's three planes coincide, so the name denotes no point.
    Collapsed,
}

/// The decline kind for a [`RingFail`] at a caller that has no finer fact to add.
///
/// ★ Two of the four callers *do* have one — they know which loop of the face failed and report
/// `OuterRing`/`HoleRing` — so this is not "the mapping", only the default.
fn decline_of(f: RingFail) -> DeclineKind {
    match f {
        RingFail::Collapsed => DeclineKind::CollapsedTriple,
    }
}

/// ★ **The one place a ring of *names* becomes a ring of *plane data*** for the tracer, which asks
/// only "which side, which wall, which third plane" of each vertex. The sort that used to stand
/// here is gone — [`NodeId`]'s only constructor sorts — so what remains is the collapse check, and
/// that is the whole reason this function exists.
fn plane_ring(
    nr: &combinatorics::NamedRing,
) -> Result<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>), RingFail> {
    // ★★★★★ **Both refusals are gone, and the vessel is why they could go.** The corner used to
    // be projected to `[usize; 3]` and the carrier to a plane class, so a ring a cylinder touched
    // could not be *described* — the refusals here were the type running out, dressed as
    // decisions. Two rungs widened the vessel; this one stops asking. What a curved ring now meets
    // is the roads themselves, which answer or decline on their own terms.
    let ts: Vec<combinatorics::NodeId> = nr
        .triples
        .iter()
        .map(|&n| {
            // A collapsed name is still a collapsed name — but only a three-plane one can be.
            if let Some(c) = three_plane_name(n) {
                // Sorted by construction, so equal neighbours catch every duplicate.
                if c[0] == c[1] || c[1] == c[2] {
                    return Err(RingFail::Collapsed);
                }
            }
            Ok(n)
        })
        .collect::<Result<_, _>>()?;
    Ok((ts, nr.walls.clone()))
}

/// One feature node on the line `L = W ∩ fp` — **its own name, and what pins it there**.
///
/// ★★ It used to be a bare third plane (`r`), which is the shape a three-plane point has and no
/// other. The name is now carried whole ([`NodeId`], so a corner a cylinder made can be one) and
/// the pin says **which kind** it is ([`combinatorics::EndPin`]) — the same pair a segment's ends
/// already travel as. The crossing arm writes the other pair — `Branch`/`Cylinder` — where the
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

/// Trace one **transversal or mixed** face `f` (plane `fp ≠ W`) onto plane class `wc`, as the
/// segments where `f` meets `W`. This is 2D polygon-vs-line clipping on `L = W ∩ fp`, handling
/// on-line edges (a run of `side == 0` vertices) instead of rejecting them. A **holed** face is
/// handled by scanning its inner rings into the same node list: each hole's two crossings of `L` toggle the parity sweep back to void between them,
/// carving the hole out of the emitted chord.
///
/// Every emitted segment rides `wall = fp` (its line is `W ∩ fp`); the side walls appear only as
/// the third plane naming each endpoint. `mat = orient_sign(fp)` is a per-face constant.
#[allow(clippy::too_many_arguments)]
fn trace_transversal_face(
    fp: usize,
    loops: &combinatorics::FaceLoops,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    crossings: &std::collections::HashSet<(usize, usize)>,
    out: &mut Trace,
) {
    // `fp` names a *face* (`n_out`, `orient`, the declined log); `fc` names the *plane class* it
    // lies on (triples, comparisons, predicate arguments). Every ring below is in class form, so
    // the two must not be confused — see `canon_ring`.
    let fc = plane_ix[fp].plane();
    // Each ring pairs its class-form triples with the **carried walls** the producer read off
    // the model's edges (`NamedRing`) — the Crossing arm below names an edge's wall from the
    // ride, never from the two endpoint names.
    // The flip nodes a **circle** of this face leaves on `L` — its chord's two ends — for the
    // outer and the holes alike (E3-b for the holes; the outer joined in cell ⑩), joined to the
    // scan's node list below.
    let mut circle_nodes: Vec<Node> = Vec::new();
    // The outer ring the scan walks — `None` for a circular outer, whose two chord ends are
    // already in `circle_nodes` and which has no polygon to walk.
    let outer: Option<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>)> =
        match loops.outer.as_ref().and_then(|r| r.poly()) {
            Some(nr) => match plane_ring(nr) {
                Ok(pair) => Some(pair),
                Err(f) => {
                    out.declined.push((fp, decline_of(f)));
                    return;
                }
            },
            // ★ A **circular outer is a disk — and a disk can have holes.** A recorded wall
            // class cuts the circle in a chord (a diameter when it runs through the axis), and
            // the chord's two ends are flip nodes of the same parity sweep every polygon face
            // runs, so the face's **holes carve the chord** exactly as they carve a polygon's
            // section: a tube's annular cap sectioned by a wall is two pieces, not the whole
            // chord across the bore. ★★ This arm used to emit the whole chord on a vessel of its
            // own (`ChordTrace`) and *return before reading the holes* — "a circular outer is a
            // disk face", true only while the gate refused every solid with two coaxial cylinders
            // (annulus, tube); the day the pair rule read faces, the tube's bore came through as
            // a chord piece the cap does not cover and the class refused `LabelConflict`
            // (measured, the annulus × box rows). One road now; the vessel is gone.
            //
            // Any other class misses or is left silent: another ⊥ class is parallel to the disk's
            // own and meets it nowhere, and a ∥-axis wall the gate did not record either clears
            // the disk (the gate's plane test) or cuts an **irrational** chord this road cannot
            // state yet (the face-cleared family — its wrongly-silent trace lands in cells no
            // face reaches, today's recorded state). A disk the wall clears has holes the wall
            // clears too, so the silence covers the holes as well.
            //
            // `None` is also how an *unnamed* ring arrives, but not for a face the tracer
            // reaches: the only loops `loop_triples` leaves unnamed are the ones it rejects for,
            // and those are `Err` at the source. A circle is the one `None` that means "named,
            // and not a polygon".
            None if matches!(loops.outer, Some(combinatorics::LoopRing::Circle { .. })) => {
                let Some(combinatorics::LoopRing::Circle { cyl }) = loops.outer else {
                    unreachable!("the guard just matched a circle outer");
                };
                match chord_nodes(jd, faces, plane_ix, wc, fc, cyl, crossings) {
                    Ok(Some(pair)) => circle_nodes.extend(pair),
                    Ok(None) | Err(ChordFail::NoCoefficients) => return,
                    Err(ChordFail::Unstatable) => {
                        out.declined.push((fp, DeclineKind::Ruling));
                        return;
                    }
                }
                None
            }
            None => {
                #[cfg(test)]
                decline_probe::mark(decline_probe::Site::NoOuter);
                out.declined.push((fp, DeclineKind::OuterRing));
                return;
            }
        };
    let mut holes: Vec<(Vec<combinatorics::NodeId>, Vec<crate::boolean::Wall>)> = Vec::new();
    // A hole whose ring cannot be named is not "no hole" — swallowing the error would trace the
    // face as solid where it is pierced, which is a silent wrong answer rather than a reject.
    let Some(raw_holes) = &loops.holes else {
        out.declined.push((fp, DeclineKind::HoleRing));
        return;
    };
    for r in raw_holes {
        match r {
            // A lateral's own rim is a cycle of a lateral row (`FaceLoops::cycles`), never a
            // plane face's hole: reaching here is a producer inconsistency, declined by the
            // loop it arrived as.
            combinatorics::LoopRing::Rim { .. } => {
                out.declined.push((fp, DeclineKind::HoleRing));
                return;
            }
            // ★ A **circular hole meets the line in a chord, and then it flips parity twice**
            // (E3-b). The population gate admits a ∥-axis wall clear of the hole's cylinder
            // (skip: the circle cannot meet `L`) or **recorded** — within the radius, through
            // the axis or offset from it (cell ③) — and there the line cuts the hole in the
            // chord's two branch points, which used to be skipped "by proof": the proof covered
            // the clear half only, and the planted full-width chord surfaced as `LabelConflict`
            // on the through-family × through-axis tool (measured). The two roots are pushed as
            // flip nodes — the hole's own chord, named the same way a disk outer's is
            // ([`chord_nodes`], the one spelling for both).
            //
            // ★ **The gate's record, the same one the rulings road reads first.** A pair the
            // gate did not list was proven clear — no emitted face reaches the hole's footprint
            // on this class, and the rulings road is silent for it — so planting flip nodes
            // would break the line against rulings that are (rightly) absent: measured, the
            // d = 0 boss-and-bore fixture walked its spur out and back and refused
            // `StraightAngle`. Listed pairs get the nodes; unlisted pairs keep today's skip, and
            // the two roads stay one rule.
            combinatorics::LoopRing::Circle { cyl } => {
                match chord_nodes(jd, faces, plane_ix, wc, fc, *cyl, crossings) {
                    Ok(Some(pair)) => circle_nodes.extend(pair),
                    Ok(None) => {}
                    Err(_) => {
                        out.declined.push((fp, DeclineKind::HoleRing));
                        return;
                    }
                }
                continue;
            }
            combinatorics::LoopRing::Poly(r) => match plane_ring(r) {
                Ok(pair) => holes.push(pair),
                Err(f) => {
                    out.declined.push((fp, decline_of(f)));
                    return;
                }
            },
        }
    }

    // `L`'s third plane naming a point of the on-line edge: the vertex triple `{fc, W-class, r}`.
    //
    // Failing to name it has two *different* causes, and the caller must be able to tell them
    // apart. `t`'s three classes are distinct (`plane_ring` rejects a collapsed triple first) and
    // `fc != wc` is this function's precondition, so:
    //
    // - `wc ∉ t` — every caller passes a **run** vertex, whose point lies on `wc` (that is what
    //   `side == 0` says). A point on `wc` that `wc` does not name, plus `t`'s three distinct
    //   classes, is **four distinct planes through one point** — an identity the plane-triple
    //   substrate cannot express, whether or not `fc` is among them.
    // - `wc ∈ t` but `fc ∉ t` — the face's own plane does not name its own vertex. That is a
    //   naming anomaly, not a concurrency, and keeps the generic `RunName`.
    //
    // Counting to the end rather than returning early matters: a four-plane vertex has *two*
    // off-plane classes, so an early "two off-planes" bail would exit before `wc`'s absence is
    // ever noticed and report every such point as `RunName`.
    let third_on_l =
        |n: NodeId, out: &mut Trace| -> Result<(NodeId, combinatorics::EndPin), DeclineKind> {
            // ★★ **A corner a cylinder made needs no third plane, and has none.** It is already the
            // point's own name, and what pins it on `L` is the quadric — measured, every such corner
            // that lands on a cut line is pinned by that very cut plane (`far == wc`), so `wc` is in
            // its name and there is no fourth-plane alias to record either.
            let (n, t) = match three_plane_name(n) {
                Some(t) => (n, t),
                None => match combinatorics::branch_name(n) {
                    // Pinned by the quadric: its own pair is this line's.
                    Some((planes, _, _)) if planes.contains(&wc) => {
                        return Ok((n, combinatorics::EndPin::Cylinder));
                    }
                    // ★ **On `wc` by coincidence** (cell ⑩): a wall through a fillet's axis runs
                    // through the fillet's tangent corner, so the corner is on `L` while neither of
                    // its own planes is `wc`. Three planes pass through the point — its own two and
                    // `wc` — so it is a **three-plane point**, named by them like any other; the
                    // cylinder is a fourth carrier that names nothing here. (Restating it as a
                    // branch of a new pair kept a root that pair does not have — measured, the
                    // fold's gusset beside the plate's fillet.)
                    Some((planes, _, _)) => {
                        // ★ Cell ⑫: the **discovery event** for a point on a cylinder — this
                        // corner has a third plane through it, so every name that set yields is
                        // one point; the table learns it here, once, and the ruling sweep on
                        // this same class reads it (`trace_one` walks the plane faces first).
                        // The identity of this corner with `[p0, p1, wc]` and with the class's crossing of
                        // its ruling is the seed's (`seed_from_operands`, cell ⑫), known before any trace.
                        let mut t = [planes[0], planes[1], wc];
                        t.sort_unstable();
                        (NodeId::three_planes(Canon3::three(t)), t)
                    }
                    None => return Err(DeclineKind::RunName),
                },
            };
            // ★ Cell ⑪: the face's own plane need not appear in the name — a ring vertex is on
            // its face by topology, and a concurrency's canonical triple may name it by three
            // *other* planes. One question remains: does `wc` name it? If so the pin is whichever
            // plane of `t` cuts `L` (`pin_on_line`; with `t = {fc, wc, r}` that is `r`, the old
            // ordinary arm, measured identical). If not, `wc` is a further plane through the
            // point: record the concurrency, pin the same way, and call the point by the
            // canonical triple of what is now known through it — the rule the operand's ring
            // and the alias table use, so the three agree on the representative.
            if t.contains(&wc) {
                combinatorics::pin_on_line(jd, wc, fc, t)
                    .map(|r| (n, combinatorics::EndPin::Class(r)))
                    .ok_or(DeclineKind::NoPinOnLine)
            } else {
                let mut set = t.to_vec();
                set.push(wc);
                set.sort_unstable();
                set.dedup();
                out.aliases.record(jd, &set);
                // ★ A handle has a duty the identity does not: it must **cut** `L` — see
                // [`combinatorics::pin_on_line`], where that rule and its reason live.
                // ★ `FourPlane` and not `NoPinOnLine`: this arm has already established that
                // `wc` is a fourth plane through the point, so the substrate limit is the cause
                // and a missing pin is its symptom.
                combinatorics::pin_on_line(jd, wc, fc, t)
                    .ok_or(DeclineKind::FourPlane)
                    .map(|r| {
                        let name = combinatorics::canonical_triple(jd, &set)
                            .map(NodeId::three_planes)
                            .unwrap_or(n);
                        (name, combinatorics::EndPin::Class(r))
                    })
            }
        };

    // Phase A — scan the outer ring and every hole ring, collecting feature nodes into ONE list.
    // A hole ring keeps its stored CW winding, but the Phase-A predicates are winding-agnostic; the
    // parity sweep (Phase C) then carves the hole because its two crossings of `L` toggle parity
    // back to void between them. A ring that does not meet `W` (every vertex one side) yields no
    // nodes; `run_counter` is shared so run ids stay unique across rings.
    let mut nodes: Vec<Node> = Vec::new();
    nodes.append(&mut circle_nodes);
    let mut run_counter = 0usize;
    let mut declined: Option<DeclineKind> = None;
    'rings: for (ring, walls) in outer.iter().chain(holes.iter()) {
        let n = ring.len();
        // Where this ring meets `L`, and whether it crosses or only touches — the walk the ray
        // caster shares (`combinatorics::ring_against_plane`). What is done with a feature is this
        // function's own business: naming it, recording a four-plane alias, deciding occupancy.
        // ★ **What "on the meet" means here**: `W` cuts this *planar* face in a straight **line**,
        // and two points fix a line — so a straight edge between two on-line nodes is on it and a
        // curved one is not. (The lateral road's meet is a **circle**, where an arc *can* lie on it,
        // so it tests the carrier instead — see `cycle_on_class`.)
        // An arc between two on-line nodes leaves the line to the side its tangent points
        // (E3-c, [`combinatorics::arc_departure_side`]); a straight edge stays on it.
        let on_meet = |i: usize| match walls[i] {
            crate::boolean::Wall::Arc { cyl, ccw } => {
                combinatorics::arc_departure_side(jd, cyls, ring[i], wc, cyl, ccw)
                    .map(combinatorics::EdgeMeet::Departs)
            }
            crate::boolean::Wall::Plane(_) | crate::boolean::Wall::Ruling { .. } => {
                Some(combinatorics::EdgeMeet::On)
            }
        };
        let features = match combinatorics::ring_against_plane(jd, cyls, ring, wc, on_meet) {
            combinatorics::RingWalk::Met(f) => f,
            // Every vertex on `W`: a ring lying in the cut plane is degenerate here.
            combinatorics::RingWalk::AllOn => {
                declined = Some(DeclineKind::AllOnPlane);
                break;
            }
            // ★ Not "a corner a cylinder made" any more — the walk reads those. This is the
            // walk's own `None`: a side it could not form exactly, or an arc tangent to the
            // class at its end (E3-c). ☑ Measured 0 raises across the workspace suite; the name
            // is kept because the walk can still say it.
            combinatorics::RingWalk::Unnameable => {
                declined = Some(DeclineKind::BranchNode);
                break;
            }
        };
        for feature in features {
            match feature {
                combinatorics::Feature::Crossing { edge, .. } => {
                    // The crossed edge's wall, **carried** from the producer — it used to be
                    // re-derived from the two endpoint names (`ring_from_names`), which is
                    // sound only while every vertex lies on exactly three planes and could
                    // hand back a plane the edge does not ride at a concurrency.
                    //
                    //
                    // ★ A crossing on a **ruling** is a point on the cylinder — a branch node
                    // `wc ∩ fc ∩ cyl`, pinned by the quadric ([`crossing_on_ruling`]); the second
                    // operation on a wall boss makes one wherever a ⊥ cap crosses the plate's
                    // wall face along the boss's rulings (the crossing census's mid slab). A
                    // crossing on an **arc** still declines: the lateral's ruling road cannot
                    // yet cut a ruling at a hole it does not carry on the class (E2-2's sweep),
                    // so naming the arc's crossing here would meet a phantom ruling there.
                    let (id, pin) = match walls[edge] {
                        crate::boolean::Wall::Plane(w) => (
                            NodeId::three_planes(Canon3::three([wc, fc, w])),
                            combinatorics::EndPin::Class(w),
                        ),
                        crate::boolean::Wall::Ruling { cyl, side, .. } => {
                            let Some(wcy) = cyls.get(cyl) else {
                                declined = Some(DeclineKind::CurvedRingWall);
                                break 'rings;
                            };
                            match crossing_on_ruling(jd, &wcy.def, wc, fc, cyl, side) {
                                Ok(id) => {
                                    crossing_probe::record(
                                        jd, &wcy.def, cyl, wc, fc, side, false, id,
                                    );
                                    (id, combinatorics::EndPin::Cylinder)
                                }
                                Err(d) => {
                                    declined = Some(d);
                                    break 'rings;
                                }
                            }
                        }
                        crate::boolean::Wall::Arc { cyl, ccw } => {
                            let Some(wcy) = cyls.get(cyl) else {
                                declined = Some(DeclineKind::CurvedRingWall);
                                break 'rings;
                            };
                            let b = ring[(edge + 1) % ring.len()];
                            match crossing_on_arc(jd, &wcy.def, wc, fc, cyl, ccw, ring[edge], b) {
                                Ok(id) => {
                                    crossing_probe::record(jd, &wcy.def, cyl, wc, fc, 0, true, id);
                                    (id, combinatorics::EndPin::Cylinder)
                                }
                                Err(d) => {
                                    declined = Some(d);
                                    break 'rings;
                                }
                            }
                        }
                    };
                    // A crossing is strict on both ends and its carrier is a line or a ruling —
                    // straight either way — so it is met exactly once: parity flips.
                    nodes.push(Node {
                        id,
                        pin,
                        flip: true,
                        run: None,
                        flanks_differ: false,
                        single_touch: false,
                        graze_above: None,
                    });
                }
                combinatorics::Feature::Run {
                    first,
                    len: m,
                    flanks_differ,
                    ..
                } => {
                    // ★ A run's end may be a corner a cylinder made (a half-disk cap's chord, E3-c);
                    // `third_on_l` names that one by itself, pinned by the quadric.
                    let name = |k: usize, out: &mut Trace| third_on_l(ring[(first + k) % n], out);
                    if m == 1 {
                        match name(0, out) {
                            Ok((id, pin)) => nodes.push(Node {
                                id,
                                pin,
                                flip: flanks_differ,
                                run: None,
                                flanks_differ,
                                single_touch: !flanks_differ,
                                // A single-vertex tangential touch is a point (`touches`), never a
                                // graze segment; a flank-differing single vertex is a strict
                                // crossing.
                                graze_above: None,
                            }),
                            Err(kind) => declined = Some(kind),
                        }
                    } else {
                        // m >= 2: one on-line interval. Only its two ends and the flanks decide
                        // anything; the interior points are names the arrangement may split at.
                        let id = run_counter;
                        run_counter += 1;
                        let names: Result<Vec<(NodeId, combinatorics::EndPin)>, DeclineKind> =
                            (0..m).map(|k| name(k, out)).collect();
                        match names {
                            Ok(rs) => {
                                // Every run is one-sided: it is an *edge* of `f` lying on `L`, so
                                // `f` is on one side of it whatever the ring does afterwards.
                                // `flanks_differ` says only whether the sweep's parity toggles here
                                // (Phase B gives it to `flip`) — it is not an occupancy fact, and
                                // gating the side on it left a crossing run classified as a
                                // straddling transversal.
                                let Some(ba) = run_body_above(jd, cyls, faces, wc, fc, fp, &rs)
                                else {
                                    declined = Some(DeclineKind::BranchNode);
                                    break 'rings;
                                };
                                let graze_above = Some(ba);
                                for (nid, pin) in rs {
                                    nodes.push(Node {
                                        id: nid,
                                        pin,
                                        flip: false,
                                        run: Some(id),
                                        flanks_differ,
                                        single_touch: false,
                                        graze_above,
                                    });
                                }
                            }
                            Err(kind) => declined = Some(kind),
                        }
                    }
                }
            }
            if declined.is_some() {
                break 'rings;
            }
        }
    }
    if let Some(reason) = declined {
        out.declined.push((fp, reason));
        return;
    }

    // Phase B — order the nodes along L and fix run structure.
    // ★ A comparator cannot decline, so it raises a flag the caller reads. What can make it fire
    // has narrowed: the second road landed, so a cylinder-pinned node is *ordered* here rather
    // than refused, and the only `None` left is a description that could not be formed exactly.
    // ☑ Still unexercised in the corpus — measured through `reject_census`, which does not move.
    let mut unordered = false;
    let order = |a: &Node, b: &Node, unordered: &mut bool| -> std::cmp::Ordering {
        let Some(o) = combinatorics::order_pinned(jd, cyls, wc, fc, (a.id, a.pin), (b.id, b.pin))
        else {
            *unordered = true;
            return std::cmp::Ordering::Equal;
        };
        match o {
            -1 => std::cmp::Ordering::Less,
            1 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        }
    };
    nodes.sort_by(|a, b| order(a, b, &mut unordered));
    if unordered {
        out.declined.push((fp, DeclineKind::BranchNode));
        return;
    }
    for w in nodes.windows(2) {
        if order(&w[0], &w[1], &mut unordered) == std::cmp::Ordering::Equal {
            out.declined.push((fp, DeclineKind::CoincidentFeatures));
            return;
        }
    }
    // Each run's two nodes must be adjacent after the sort; the hi node carries the flip.
    for k in 0..nodes.len() {
        if let Some(id) = nodes[k].run {
            let prev = k > 0 && nodes[k - 1].run == Some(id);
            let next = k + 1 < nodes.len() && nodes[k + 1].run == Some(id);
            if !prev && !next {
                out.declined.push((fp, DeclineKind::RunSplit));
                return;
            }
            // The whole run is one interval, so its flip happens once — at the far end.
            if prev && !next {
                nodes[k].flip = nodes[k].flanks_differ;
            }
        }
    }

    // Phase C — sweep left to right; covered gaps merge into segments.
    let mut parity = 0i8;
    // `seg_start` also carries the graze side of the segment being built: `Some(ba)` when the
    // segment is a pure graze gap (opened by a forced run while outside material), else `None`.
    let mut seg_start: Option<((NodeId, combinatorics::EndPin), Option<bool>)> = None;
    let emit = |a: (NodeId, combinatorics::EndPin),
                b: (NodeId, combinatorics::EndPin),
                graze: Option<bool>,
                out: &mut Trace| {
        out.segs.push(Seg {
            wall: plane_ix[fp].plane(),
            end: [a.0, b.0],
            end_h: [a.1, b.1],
            solid: which,
            kind: match graze {
                Some(body_above) => SegKind::Graze { body_above },
                None => SegKind::Transversal {
                    mat: faces[fp].plane().orient_sign,
                },
            },
        });
    };
    for k in 0..nodes.len() {
        if nodes[k].single_touch && parity == 0 {
            out.touches.push(nodes[k].id);
        }
        if nodes[k].flip {
            parity ^= 1;
        }
        let forced =
            k + 1 < nodes.len() && nodes[k].run.is_some() && nodes[k].run == nodes[k + 1].run;
        let covered = parity == 1 || forced;
        // The kind of the gap *after* node `k`, as a fact about that gap alone: riding a run means
        // `L` is on `f`'s boundary there, so `f` occupies one side; any other covered gap has `L`
        // in `f`'s interior, so `f` straddles. `parity` says whether the gap is covered, not which
        // arcs it fills — reading occupancy off it classified a hole's edge as a straddle.
        let gap_graze = forced.then(|| {
            nodes[k]
                .graze_above
                .expect("every run node carries its occupied side")
        });
        match (seg_start, covered) {
            (None, true) => seg_start = Some(((nodes[k].id, nodes[k].pin), gap_graze)),
            // The kind changes here, so close and reopen: a segment must be homogeneous, because
            // `split_at_crossings` subdivides it later and every piece inherits its kind.
            (Some((a, graze)), true) if graze != gap_graze => {
                emit(a, (nodes[k].id, nodes[k].pin), graze, out);
                seg_start = Some(((nodes[k].id, nodes[k].pin), gap_graze));
            }
            (Some((a, graze)), false) => {
                emit(a, (nodes[k].id, nodes[k].pin), graze, out);
                seg_start = None;
            }
            _ => {}
        }
    }
    if seg_start.is_some() || parity != 0 {
        // A line enters and leaves a bounded region equally; an unbalanced sweep is degenerate.
        out.declined.push((fp, DeclineKind::OddParity));
    }
}

/// Which side of `W` face `fp` occupies along an on-line run — **the left of the ring's travel**,
/// stated in the label frame (`true` = above).
///
/// A run is an *edge* of `f`'s boundary lying on `L = W ∩ fp`, so along it `f` is on exactly one
/// side, and which side is the universal boundary convention: **material is on the left of the
/// ring's direction of travel**. That holds for an outer ring, a hole ring, a notch, and a reflex
/// corner alike — it does not care how the ring continues past the run.
///
/// The earlier version read the side off the run's *flank* (the neighbouring off-line vertex)
/// instead. That is only a proxy for "which way the ring bulges", and it is **inverted wherever the
/// ring turns away from material**: on a hole ring the neighbours point into the hole, and on an
/// outer-ring notch they point across the notch — both void. Measured: over the whole corpus the
/// two agree on 251 of 252 runs, and the one disagreement (a fused boss's `z=1` annulus, where the
/// run is the hole's edge) is the one the flank gets wrong.
///
/// **Derivation.** `order_along` runs along `d = n_s(wc) × n_s(fc)` (stored normals — see
/// [`combinatorics::dir_sign`] and [`Judge::plane_pair_dir_sign`](nacre_cip::predicate::Judge)), and the label frame's "above" is
/// `n_s(wc)` (see [`combinatorics::side_of`]). With `t = order_along(wc, fc, first, last)` the travel
/// is `−t·d`, so material points along `n_out(fp) × (−t·d)`. Since `fp` is coplanar with its class
/// root, `n_out(fp) = σ·n_s(fc)`, and the triple product collapses to
/// `−t·σ·(1 − (n_s(fc)·n_s(wc))²)` — whose bracket is positive because the two planes are not
/// parallel. Hence `body_above = (t·σ < 0)`. Note `orient_sign(wc)` **cancels**: the ordering
/// direction and the label frame are defined by the same stored normal.
///
/// `σ` is an f64 dot of two **parallel** unit vectors (`fp` and `fc` are the same plane class), so
/// `|σ| ≈ 1` — a full unit from the sign boundary, the same robustness [`FaceInfo::orient_sign`] and
/// `trace_seated_face` already rely on. Everything else here is exact.
#[allow(clippy::too_many_arguments)]
fn run_body_above(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    fc: usize,
    fp: usize,
    rs: &[(NodeId, combinatorics::EndPin)],
) -> Option<bool> {
    let planes = jd.planes;
    // ★ The run's own direction along `L`, through the one rule both roads share — a run whose
    // ends a cylinder pinned orders by the `a + b√c` tower, and by the integer predicates
    // otherwise. `None` is a missing description, which the caller turns into its own decline.
    let t = combinatorics::order_pinned(jd, cyls, wc, fc, rs[0], rs[rs.len() - 1])?;
    let sigma = faces[fp].plane().n_out.dot(planes[fc].plane.normal());
    Some((t < 0) == (sigma > 0.0))
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
///
/// `faces_in` is that solid's `(slot, rings)` pairs — see [`combinatorics::TraceInput`]. This walk
/// used to read the shells out of the `Model` and derive the rings here, once per class; nothing
/// about either depends on `wc`.
#[allow(clippy::too_many_arguments)]
fn trace_one(
    faces_in: &[(usize, combinatorics::FaceLoops)],
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    crossings: &std::collections::HashSet<(usize, usize)>,
    aliases: &Aliases,
    out: &mut Trace,
) {
    let planes = jd.planes;
    let w_normal = planes[wc].plane.normal();
    // ★ Every curved carrier a seated ring rode and this walk declined to re-emit, with the face
    // that rode it. Checked once at the end — a lateral face may be visited after the seated one,
    // so the question is only answerable when the solid's whole contribution is in.
    let mut curves_owed: Vec<(usize, crate::boolean::Wall)> = Vec::new();
    for (fp, fl) in faces_in {
        let fp = *fp;
        // A **lateral face** (a cylinder row): what it leaves on a ⊥ class is a **transversal**
        // circle where the class cuts it in two and a **graze** where the class carries one of its
        // rims — see [`circle_on_class`]. Beyond the span it leaves nothing. Non-⊥ classes cannot
        // carry a circle, and the population gate has already named every such interaction.
        if let ClassIx::Cyl(k) = plane_ix[fp] {
            let cf = match &faces[fp] {
                FaceRow::Cylinder(cf) => cf,
                FaceRow::Plane(_) => unreachable!("ClassIx::Cyl marks a cylinder row"),
            };
            match circle_on_class(jd, cyls, cf, fl, wc, k) {
                Ok(spans) if !spans.is_empty() => {
                    out.circles.extend(spans.into_iter().map(|s| CircleTrace {
                        cyl: k,
                        solid: which,
                        kind: match s.on {
                            CylOnClass::Crosses => SegKind::Transversal {
                                mat: cf.orient_sign,
                            },
                            CylOnClass::Grazes { body_above } => SegKind::Graze { body_above },
                        },
                        arc: s.arc,
                    }));
                }
                // No circle: a wall class **within the radius** leaves two rulings instead (the
                // M6-2 rulings road, the offset wall since cell ③); any other non-⊥ class still
                // leaves nothing, silently — the population gate names those interactions.
                Ok(_) => match rulings_on_class(jd, cf, fl, wc, k, which, crossings, aliases) {
                    Ok((v, grazes)) => {
                        out.rulings.extend(v);
                        // The lateral's side of a plane-pair line it touches (cell ⑫).
                        out.segs.extend(grazes);
                    }
                    Err(kind) => out.declined.push((fp, kind)),
                },
                Err(kind) => out.declined.push((fp, kind)),
            }
            continue;
        }
        if plane_ix[fp].plane() != wc {
            trace_transversal_face(fp, fl, which, wc, jd, cyls, faces, plane_ix, crossings, out);
            continue;
        }
        // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
        // side of W — `n_out` points away from the body, so the body is above W exactly when
        // `n_out · n_W < 0`. The f64 sign is robust even rotated: seated means `canon[fp]==wc`,
        // so `n_out ∥ n_W` (both unit) and the dot is ≈ ±1, a full unit from the sign boundary
        // (the rotated-tunnel tests exercise this seated path through the cube's own caps).
        let body_above = faces[fp].plane().n_out.dot(w_normal) < 0.0;
        let kind = SegKind::Seated { body_above };
        // Collect every ring in class form first: a collapsed name declines the whole face, and
        // deciding that before the emitting closure exists keeps the two borrows apart.
        // A **circular** boundary — a disk face's outer, or a circular drill hole — is a
        // closed trace element ([`CircleTrace`]) rather than a segment ring.
        // ★ The **walls travel with the triples**. `plane_ring` has handed both back since the
        // vessel widened, and this caller dropped one of them on the floor.
        let mut rings: Vec<(Vec<NodeId>, Vec<crate::boolean::Wall>)> = Vec::new();
        match &fl.outer {
            Some(combinatorics::LoopRing::Poly(nr)) => match plane_ring(nr) {
                Ok(r) => rings.push(r),
                // ★ `OuterRing` here, not `CollapsedTriple`: this caller reports *which loop*
                // failed, which is the finer fact when a face has several. The two causes that
                // used to keep their own names here were the vessel's, and went with it.
                Err(RingFail::Collapsed) => {
                    #[cfg(test)]
                    decline_probe::mark(decline_probe::Site::PolyCollapsed);
                    out.declined.push((fp, DeclineKind::OuterRing));
                    continue;
                }
            },
            Some(combinatorics::LoopRing::Circle { cyl }) => {
                out.circles.push(CircleTrace {
                    cyl: *cyl,
                    solid: which,
                    kind,
                    // A seated circle is a whole loop of the face lying in the class: the disk
                    // cap's own rim, or a circular hole through it. It has no partial extent.
                    arc: None,
                });
            }
            // A lateral's rim never names a plane face's outer loop (see `FaceLoops::cycles`).
            Some(combinatorics::LoopRing::Rim { .. }) | None => {
                #[cfg(test)]
                decline_probe::mark(decline_probe::Site::RimOrNone);
                out.declined.push((fp, DeclineKind::OuterRing));
                continue;
            }
        }
        let mut failed: Option<DeclineKind> = None;
        // As above: an unnameable hole is a reject, not "no hole".
        let Some(raw) = &fl.holes else {
            out.declined.push((fp, DeclineKind::HoleRing));
            continue;
        };
        for r in raw {
            match r {
                combinatorics::LoopRing::Circle { cyl } => out.circles.push(CircleTrace {
                    cyl: *cyl,
                    solid: which,
                    kind,
                    arc: None,
                }),
                combinatorics::LoopRing::Rim { .. } => failed = Some(DeclineKind::HoleRing),
                combinatorics::LoopRing::Poly(nr) => match plane_ring(nr) {
                    Ok(r) => rings.push(r),
                    Err(f) => failed = Some(decline_of(f)),
                },
            }
        }
        if let Some(kind) = failed {
            out.declined.push((fp, kind));
            continue;
        }
        let fc = plane_ix[fp].plane();
        let mut emit_ring = |ns: &[combinatorics::NodeId], ws: &[crate::boolean::Wall]| {
            // ★★★★★ **The carrier is taken, not derived.** It used to be re-read out of the two
            // endpoint *names* — "the class they share besides `fc`" — which is sound only while
            // every vertex lies on exactly three planes, and at a four-plane concurrency can hand
            // back a plane the edge does not ride. The tracer's scan road stopped doing that on
            // 2026-08-17 and this road, its sibling, was left on the old derivation.
            // `NamedRing.walls` is what the producer read off the **edge's own two faces**, total
            // even where the vertex names fall back. Measured across the suite before the swap:
            // the two agree 810,925 times and differ 0, and the carried one additionally answers
            // 40 edges the derivation cannot name at all.
            //
            // ★ The corner keeps its own name too: `NodeId::three_planes(t)` was a round trip
            // through a projection, and `three_plane_name` is pure extraction of an already-sorted
            // triple, so the two are the same value.
            let n = ns.len();
            for i in 0..n {
                // ★★★★★ **A curved edge emits nothing, and that is a refusal to say a thing
                // twice rather than a hole.** A [`Seg`] is a straight edge on `wc`; an arc or a
                // ruling is not that species. The element it needs is already on this class,
                // contributed by the **cylinder's own lateral face** of this same solid — the
                // circle where a ⊥ class cuts it, the two rulings where a class carries its axis.
                // Measured on every class where this fires: `circles = 1` on the cap classes and
                // `rulings = 2` on the wall class, all present before the seated face is reached.
                // ★ The claim is checked at the end of this walk rather than trusted — see
                // `curves_owed` — because an element silently missing does not decline, it
                // corrupts every label on the class.
                // ★ Except a **tangent** ruling (`side == 0`, cell ⑩): the lateral touches this
                // class along it and crosses nowhere, so no ruling trace will ever back it — the
                // seated face is the only face on this class with that line as an edge, and it
                // states the piece itself, `Seated` like its straight edges. `up` says which end
                // is the lower one along the axis (`MergedRuling::end` ascends).
                if let crate::boolean::Wall::Ruling { cyl, side: 0, up } = ws[i] {
                    let (a, b) = (ns[i], ns[(i + 1) % n]);
                    // The axis side, in the label frame (`W`'s stored normal): the one bit a
                    // tangent edge flips ([`SegKind::Tangent`]).
                    let axis_above = (|| -> Option<bool> {
                        let w = combinatorics::stored_coeffs_rat(jd, wc)?;
                        let o = cyls.get(cyl)?.def.origin();
                        let mut acc = w[3];
                        for k in 0..3 {
                            acc = acc.checked_add(w[k].checked_mul(o[k])?)?;
                        }
                        Some(acc > Rat::from_int(0))
                    })();
                    let Some(axis_above) = axis_above else {
                        out.declined.push((fp, DeclineKind::Ruling));
                        continue;
                    };
                    out.rulings.push(RulingTrace {
                        cyl,
                        side: 0,
                        end: if up { [a, b] } else { [b, a] },
                        solid: which,
                        kind: SegKind::Tangent { axis_above },
                        #[cfg(test)]
                        orient: faces[fp].plane().orient_sign,
                    });
                    curves_owed.push((fp, ws[i]));
                    continue;
                }
                let crate::boolean::Wall::Plane(wall) = ws[i] else {
                    curves_owed.push((fp, ws[i]));
                    continue;
                };
                // ★★★★★ **This test comes first, and the order is what keeps two causes apart.**
                // A carrier that *is* this face's own class makes `fc ∩ wall` not a line, and
                // `plane_pair_dir_sign(fc, fc, ·)` is a determinant with two equal rows — zero
                // against every candidate. Ask the pin first and the pencil case arrives wearing
                // `NoPinOnLine`, which is the wrong sentence for it.
                // ☑ Measured 0 in the corpus. The derivation this replaced could not produce it
                // (it filtered `c != fc`), so the guard is what carries that property across.
                // ☑ This arm and the pin's below were both made `unreachable!()` with the workspace
                // suite and the ignored sweep green — that is what "unexercised" means here, rather
                // than an argument that they cannot fire. (The curved arm above is *not* in that
                // set: it fires, 20 times, and `curves_owed` is what checks it.)
                if wall == fc {
                    out.declined.push((fp, DeclineKind::SeatedEdgeNaming));
                    continue;
                }
                // ★★ **What pins each end, in the vocabulary that names it.** A three-plane
                // corner is pinned by a plane — [`combinatorics::pin_on_line`], the same rule the
                // scan road and the ray caster ask. A corner a cylinder made has no third plane
                // and needs none: the quadric pins it, and [`combinatorics::EndPin::Cylinder`]
                // says so while the name beside it says which point.
                let pin = |k: usize| -> Option<combinatorics::EndPin> {
                    match three_plane_name(ns[k]) {
                        Some(t) => combinatorics::pin_on_line(jd, fc, wall, t)
                            .map(combinatorics::EndPin::Class),
                        None => Some(combinatorics::EndPin::Cylinder),
                    }
                };
                let (Some(h0), Some(h1)) = (pin(i), pin((i + 1) % n)) else {
                    out.declined.push((fp, DeclineKind::NoPinOnLine));
                    continue;
                };
                out.segs.push(Seg {
                    wall,
                    end: [ns[i], ns[(i + 1) % n]],
                    end_h: [h0, h1],
                    solid: which,
                    kind,
                });
            }
        };
        for (ring, walls) in &rings {
            emit_ring(ring, walls);
        }
    }
    // ★★★★ **The seated walk's silence over a curved edge, checked rather than argued.** It emits
    // no segment there because the cylinder's own lateral face already put the element on this
    // class. That follows from the population gate — a class meeting a cylinder is ⊥ to its axis
    // (a circle), carries its axis (two rulings), or was refused long before here — but an
    // argument is not a measurement, and an element that goes missing does not raise: it silently
    // relabels the class.
    //
    // ★ Asked as narrowly as the carrier speaks. A `(cylinder, class)` pair has **two** rulings and
    // a ring rides one of them, so the side is part of the question; a plane meets a cylinder in
    // one circle, so there the cylinder alone names it.
    for (fp, w) in curves_owed {
        let backed = match w {
            crate::boolean::Wall::Ruling { cyl, side, .. } => out
                .rulings
                .iter()
                .any(|r| r.cyl == cyl && r.side == side && r.solid == which),
            crate::boolean::Wall::Arc { cyl, .. } => {
                out.circles.iter().any(|c| c.cyl == cyl && c.solid == which)
            }
            crate::boolean::Wall::Plane(_) => unreachable!("only curved carriers are collected"),
        };
        if !backed {
            out.declined.push((fp, DeclineKind::SeatedCurveUnbacked));
        }
    }
}

/// What a lateral face leaves on a ⊥ plane class — see [`circle_on_class`], which answers it
/// **per angular extent** ([`CircleSpan`]) rather than once for the whole circle.
///
/// The empty answer is "this class leaves nothing here": a non-⊥ class carries no circle at all,
/// and a ⊥ class beyond both rims never meets the face. `Err` is "the face cannot answer" (no rim
/// span, or a class with no exact description the gate would already have refused) — declining
/// beats a silently missing circle, which would corrupt every label on the class.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CylOnClass {
    /// The class cuts the face in two: the solid straddles it here.
    Crosses,
    /// The class carries one of the face's **rims**: the wall touches it along that circle with
    /// its body to one side. `body_above` is that side, about the class's **stored** normal —
    /// the frame every cell label is written in.
    Grazes { body_above: bool },
}

fn circle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    cyl: usize,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    let wp = &jd.planes[wc];
    // The world description — the cylinder's statement is world, and a comparison across two
    // frames is a silently wrong answer, not a slow one.
    let Some(coeffs) = wp.world_rat else {
        return Err(DeclineKind::CylSpan);
    };
    // No world statement for this lateral (a rotated or frame-borne truth): the same decline
    // the class side makes — both sides of this comparison have to speak about the world.
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::CylSpan);
    };
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let m = def.dir();
    // ⊥ to the axis, decided **totally** (`parallel_rat` clears denominators and answers in
    // integers, so this cannot decline for want of bits). A non-⊥ class carries no circle at
    // all — a miss, like a parallel plane, not a decline.
    if !nacre_scalar::parallel_rat(&n, &m) {
        return Ok(Vec::new());
    }
    // ★★★★★ **The face's boundary as the tracer names it — its cycles — and nothing else.** A
    // band's two rims and its holes come from the outer loop cut at its slits
    // (`combinatorics::lateral_cycles`), so a hole the assembly spliced into the outer walk is
    // a hole here like any other; a panel or a chain rim (E2-2) is a cycle like a hole, carved
    // out of an outer answer that may be **absent** rather than out of a whole circle.
    let LateralShape {
        range,
        rims,
        cycles,
    } = lateral_shape(jd, cf, fl, def)?;
    let Some(t) = crate::planes::axis_param_of_plane(&coeffs, def) else {
        return Err(DeclineKind::CylSpan);
    };
    // The **outer** answer, which the holes below carve out of. Unchanged from when it was the
    // whole answer:
    //
    // ★★ **A rim is a graze, not a miss.** This used to answer "the lateral does not reach this
    // plane", on the premise that coplanarity would have folded a rim into a seated class — but a
    // *cylinder* face never becomes seated on a plane class, so the touch was simply dropped.
    // Every planar wall meeting a class along an edge contributes a `Graze`; the lateral
    // contributed nothing, and `edge_mask`'s precedence (`Graze > Seated`) exists exactly for the
    // corner where the two disagree. On a convex cap they agree and the omission is invisible; at
    // the **reflex dihedral of a blind bore's ceiling** they disagree, the seated rule then flipped
    // the wrong label bit, and the *next* boolean on that solid came back
    // `CylinderGateUndecided` — while a through bore, which has no ceiling, was fine.
    //
    // Which side the body is on: the face runs from `span[0]` toward `span[1]`, so at the low rim
    // it lies toward `+t` and at the high rim toward `−t`.
    let up = crate::planes::plus_t_is_above(wp, def);
    if t < range[0] || range[1] < t {
        // Beyond the face: this class does not meet it at all, and every cycle lives inside the
        // range, so there is nothing for the walk below to find either.
        return Ok(Vec::new());
    }
    // ★★★★★ **The outer answer is optional (E2-2).** A whole rim at this station grazes as a
    // whole circle; strictly inside the range the face crosses the class *somewhere* — a
    // connected face covers every interior station at some θ — and the cycles below carve out
    // where it does not; at an end of the range that is **no** rim (a panel's arc, a chain's top)
    // nothing stands: only the runs the cycles lay on the class speak, and the rest is absent.
    // ★ `Crosses` is exact only because a run whose flanks differ is declined below
    // (`cycle_on_class`): such a run is a parity boundary in θ, and the face beyond it at this
    // station would be absent with no crossing to carve it — a staircase chain, not today's
    // population. Pairing run ends with crossings is that road's generalization, not this one.
    let outer: Option<CylOnClass> = if rims.iter().any(|(s, _)| *s == t) {
        Some(CylOnClass::Grazes {
            body_above: if t == range[0] { up } else { !up },
        })
    } else if range[0] < t && t < range[1] {
        Some(CylOnClass::Crosses)
    } else {
        None
    };
    // ★★★★★ **The hole is read by the walk every other face's boundary is read by** — the face's
    // own loop, against this class, through [`combinatorics::ring_against_plane`]. It used to be an
    // interval derived from the loop's ⊥ carriers, which is exact for a chart rectangle and a
    // *premise* for anything else; the walk asks the ring instead and needs no premise about the
    // hole's shape.
    let mut carved: Vec<Carved> = Vec::new();
    for (kind, ring) in &cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair, not
        // this road's. It has no ring to walk, so it declines rather than passing unread.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        cycle_on_class(jd, cyls, cf, def, wc, cyl, up, *kind, nr, &mut carved)?;
    }
    let spans = if carved.is_empty() {
        match outer {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        }
    } else {
        assemble_spans(jd, cyl, &cyls[cyl].def, &carved, outer)?
    };
    cycle_probe::record(t, outer, rims.len(), &cycles, carved.len(), spans.len());
    Ok(spans)
}

/// The ⊥ road's answers, recorded: one entry per lateral face × ⊥ class it meets — the station,
/// the outer answer, what kinds of cycles the face has, and how many extents were carved and how
/// many spans came out. What says the road ran on a panel or a chain, and what it said there.
/// **Where an `OuterRing` decline was produced** (cell ⑩, S0). Three sites push that one kind and
/// the census key stops at the kind, so a fixture could say *that* it was declined but not by
/// which sentence. Test-only tally by producing site, keyed by thread name like `tie_probe` — a
/// test reads its own rows and never a parallel test's.
#[cfg(test)]
pub(crate) mod decline_probe {
    use std::sync::Mutex;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub(crate) enum Site {
        /// The transversal walk named no outer loop for the face.
        NoOuter,
        /// A seated face's polygon ring collapsed in `plane_ring` — a carrier it could not spell.
        PolyCollapsed,
        /// A seated face whose outer loop is a lateral rim, or unnamed.
        RimOrNone,
    }

    pub(crate) static ROWS: Mutex<Vec<(String, Site)>> = Mutex::new(Vec::new());

    pub(crate) fn mark(site: Site) {
        let name = std::thread::current().name().unwrap_or("?").to_string();
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push((name, site));
    }

    /// The sites produced on threads whose name contains `tag`, in order. The trace runs on
    /// rayon workers under `parallel`, so a test attributes rows by running its boolean in a pool
    /// named after itself; a sequential build's test thread carries the test's name already.
    pub(crate) fn named_like(tag: &str) -> Vec<Site> {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .iter()
            .filter(|(n, _)| n.contains(tag))
            .map(|(_, s)| *s)
            .collect()
    }
}

pub(crate) mod cycle_probe {
    use super::{CylOnClass, combinatorics};
    use std::sync::Mutex;

    #[derive(Clone, Copy, Debug)]
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) struct Hit {
        /// The class's axis station, realized.
        pub t: f64,
        pub outer: Option<CylOnClass>,
        /// `[rims, chains, panels, holes]` among the face's cycles (rims are not carved).
        pub kinds: [usize; 4],
        pub carved: usize,
        pub spans: usize,
    }

    pub(crate) static HITS: Mutex<Vec<Hit>> = Mutex::new(Vec::new());

    pub(crate) fn record(
        t: nacre_scalar::Rat,
        outer: Option<CylOnClass>,
        rims: usize,
        cycles: &[(combinatorics::CycleKind, combinatorics::LoopRing)],
        carved: usize,
        spans: usize,
    ) {
        let mut kinds = [rims, 0, 0, 0];
        for (k, _) in cycles {
            kinds[match k {
                combinatorics::CycleKind::Rim => 0,
                combinatorics::CycleKind::Chain => 1,
                combinatorics::CycleKind::Panel => 2,
                combinatorics::CycleKind::Hole => 3,
            }] += 1;
        }
        HITS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(Hit {
                t: t.to_f64(),
                outer,
                kinds,
                carved,
                spans,
            });
    }
}

/// **A lateral face's shape as its cycles state it** (E2-2) — what both lateral roads read of a
/// lateral's loops: its axial `range`, its whole rims (station and plane class, in station order),
/// and every other cycle — a panel, a chain rim, a hole — with its kind. The stations are read
/// from the planes' world coefficients, so `range` is a second derivation of
/// `CylFaceInfo::footprint.span` (the model side's: the outer loop's ⊥ carriers), and the two are held
/// equal where they meet.
///
/// Declines by the one name: `None` cycles (the outer loop could not be cut or named), a rim kind
/// that is not a whole circle, a class with no exact station, more than two rims, two rims at one
/// station, or a rim that is not at an end of the range — a whole circle bounds the face, so it is
/// an extreme of it.
fn lateral_shape(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    def: &nacre_topo::CylinderDef,
) -> Result<LateralShape, DeclineKind> {
    use combinatorics::{CycleKind, LoopRing};
    let Some(cycles) = &fl.cycles else {
        return Err(DeclineKind::CylSpan);
    };
    // One spelling of "a class's axis station" (`bands::param_opt`), not a third.
    let station = |c: usize| crate::bands::param_opt(jd, c, def).ok_or(DeclineKind::CylSpan);
    let mut rims: Vec<(Rat, usize)> = Vec::new();
    let mut rest: Vec<(CycleKind, LoopRing)> = Vec::new();
    let mut range: Option<[Rat; 2]> = None;
    let mut widen = |t: Rat| {
        range = Some(match range {
            None => [t, t],
            Some([lo, hi]) => [if t < lo { t } else { lo }, if hi < t { t } else { hi }],
        });
    };
    for (kind, ring) in cycles {
        match (kind, ring) {
            (CycleKind::Rim, LoopRing::Rim { plane }) => {
                let t = station(*plane)?;
                widen(t);
                rims.push((t, *plane));
            }
            (CycleKind::Rim, _) => return Err(DeclineKind::CylSpan),
            (_, ring) => {
                // The cycle's arcs (an edge with a stated sense) each lie on a ⊥ class; its
                // rulings have no station of their own.
                if let Some(nr) = ring.poly() {
                    for (i, w) in nr.walls.iter().enumerate() {
                        if nr.arc_ccw[i].is_none() {
                            continue;
                        }
                        let crate::boolean::Wall::Plane(c) = *w else {
                            return Err(DeclineKind::CylSpan);
                        };
                        widen(station(c)?);
                    }
                }
                rest.push((*kind, ring.clone()));
            }
        }
    }
    let Some(range) = range else {
        return Err(DeclineKind::CylSpan);
    };
    if range[0] == range[1] || rims.len() > 2 {
        return Err(DeclineKind::CylSpan);
    }
    rims.sort_by_key(|r| r.0);
    if rims.windows(2).any(|w| w[0].0 == w[1].0)
        || rims.iter().any(|(t, _)| *t != range[0] && *t != range[1])
    {
        return Err(DeclineKind::CylSpan);
    }
    debug_assert_eq!(
        Some(range),
        cf.footprint.span,
        "the cycles' stations by name and the face's range by model disagree"
    );
    Ok(LateralShape {
        range,
        rims,
        cycles: rest,
    })
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

/// **Which way a hole's arcs actually run, realized** — the independent oracle for
/// [`cycle_on_class`]'s winding rule.
///
/// ★★★★★ **The rule is derived from signs and nothing downstream reads it yet.** A reversed arc
/// would put the band's answer inside the hole and the hole's on the band — the exactly-opposite
/// answer — and today's fixtures stop at `loop_winding` before `label_cells` could notice. So every
/// extent a hole carves is **realized here** and a lock judges it against the fixture's own
/// geometry. The kernel reads no coordinate to choose an arc; this reads one afterwards, to check.
#[cfg(test)]
pub(crate) mod arc_probe {
    use super::{Judge, NodeId, WorkingPlane, combinatorics};
    use std::sync::Mutex;

    /// The direction, from the circle's centre, of the **midpoint of the stated counter-clockwise
    /// arc** — one entry per extent a cycle has carved out of a circle in this binary, whether the
    /// face grazes there (the cycle's own arc) or is absent (a hole's interior, a panel's outside).
    /// Beside it the cycle's kind and the cylinder's origin, so a reader can pick **its own**
    /// fixture's holes out of a ledger every test in the binary writes to.
    pub(crate) static MIDS: Mutex<Vec<Mid>> = Mutex::new(Vec::new());

    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Mid {
        pub dir: [f64; 3],
        pub kind: combinatorics::CycleKind,
        pub origin: [f64; 3],
    }

    pub(crate) fn record(
        jd: &Judge<'_, WorkingPlane>,
        cyl: usize,
        def: &nacre_topo::CylinderDef,
        kind: combinatorics::CycleKind,
        arc: [NodeId; 2],
    ) {
        let (Some(pa), Some(pb)) = (
            combinatorics::branch_point(jd, cyl, def, arc[0]),
            combinatorics::branch_point(jd, cyl, def, arc[1]),
        ) else {
            return;
        };
        let o = def.origin().map(|x| x.to_f64());
        let raw = def.dir().map(|x| x.to_f64());
        let ml = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
        let m: [f64; 3] = core::array::from_fn(|i| raw[i] / ml);
        // The circle's centre is not `origin` — that is a point on the axis — so the axial part
        // comes off first.
        let perp = |p: [f64; 3]| -> [f64; 3] {
            let v: [f64; 3] = core::array::from_fn(|i| p[i] - o[i]);
            let h = v[0] * m[0] + v[1] * m[1] + v[2] * m[2];
            core::array::from_fn(|i| v[i] - h * m[i])
        };
        let va = perp(pa);
        let vb = perp(pb);
        let ra = (va[0] * va[0] + va[1] * va[1] + va[2] * va[2]).sqrt();
        let u: [f64; 3] = core::array::from_fn(|i| va[i] / ra);
        // `w` completes a right-handed frame with `u` about `m`, so +90° about `m` takes `u` to
        // `w` — which is the direction θ increases in.
        let w = [
            m[1] * u[2] - m[2] * u[1],
            m[2] * u[0] - m[0] * u[2],
            m[0] * u[1] - m[1] * u[0],
        ];
        // The counter-clockwise sweep from `a` to `b`, then half of it. Written as a rotation
        // rather than as `va + vb`, which vanishes when the arc is exactly a half — and the hole
        // this measures **is** exactly a half.
        let mut phi = (vb[0] * w[0] + vb[1] * w[1] + vb[2] * w[2])
            .atan2(vb[0] * u[0] + vb[1] * u[1] + vb[2] * u[2]);
        if phi <= 0.0 {
            phi += core::f64::consts::TAU;
        }
        let (c, s) = ((phi / 2.0).cos(), (phi / 2.0).sin());
        MIDS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(Mid {
                dir: core::array::from_fn(|i| c * u[i] + s * w[i]),
                kind,
                origin: o,
            });
    }
}

/// **Lay the carved extents over the outer answer and hand back *exclusive* spans.**
///
/// ★★★★★ **Exclusive, not "the whole circle plus overrides".** The cheaper spelling would emit the
/// outer answer over the whole circle and let the hole's arcs sit on top, leaning on `edge_mask`'s
/// `Graze > Transversal` precedence to pick the winner. That precedence is for **different faces
/// meeting on one edge**, not for one face's own extent — borrowing it would be a category error,
/// and it cannot express "absent" at all. So the circle is cut at every boundary and each piece
/// carries exactly one answer.
///
/// Adjacent pieces with the same answer are **merged**, so a boundary where nothing changes leaves
/// no cut behind: an extra split point would be a degree-2 vertex the winding walk has to have a
/// turn for.
fn assemble_spans(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    carved: &[Carved],
    outer: Option<CylOnClass>,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    let mut nodes: Vec<NodeId> = carved.iter().flat_map(|c| c.arc).collect();
    nodes.sort_unstable();
    nodes.dedup();
    let (order, _) = circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
        CircleOrderFail::Undecided => DeclineKind::CylSpan,
        CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
    })?;
    let n = order.len();
    // θ position of each node, keyed by its index in the sorted `nodes`.
    let mut at = vec![0usize; n];
    for (k, &i) in order.iter().enumerate() {
        at[i] = k;
    }
    let place = |x: NodeId| {
        at[nodes
            .binary_search(&x)
            .expect("every boundary node is in the set")]
    };
    // The answer over elementary arc `k` (from `order[k]` to `order[k+1]`), or `None` where no
    // cycle reaches: there the outer answer stands — which may itself be "absent" (E2-2).
    let mut ans: Vec<Option<Option<CylOnClass>>> = vec![None; n];
    for c in carved {
        let (a, b) = (place(c.arc[0]), place(c.arc[1]));
        if a == b {
            return Err(DeclineKind::CylHoleFeature); // an extent of no length, or of the whole circle
        }
        let mut k = a;
        while k != b {
            if ans[k].is_some() {
                return Err(DeclineKind::CylHoleFeature); // two cycles claiming one extent
            }
            ans[k] = Some(c.on);
            k = (k + 1) % n;
        }
    }
    let answer = |k: usize| ans[k].unwrap_or(outer);
    // One answer everywhere: the circle was never divided, whatever the holes touched.
    let Some(start) = (0..n).find(|&k| answer(k) != answer((k + n - 1) % n)) else {
        return Ok(match answer(0) {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        });
    };
    let mut out = Vec::new();
    let mut k = 0;
    while k < n {
        let i = (start + k) % n;
        let a = nodes[order[i]];
        let this = answer(i);
        let mut len = 1;
        while k + len < n && answer((start + k + len) % n) == this {
            len += 1;
        }
        if let Some(on) = this {
            out.push(CircleSpan {
                arc: Some([a, nodes[order[(start + k + len) % n]]]),
                on,
            });
        }
        k += len;
    }
    Ok(out)
}

/// **Which way round a lateral face's material lies, for a boundary edge travelling along the
/// axis** — the one place that sign is made, and the third application of one convention.
///
/// *Material is on the left of the ring's direction of travel* — the sentence [`run_body_above`]
/// states for a planar face's on-line run. On a lateral face, with `r̂` the outward radial
/// direction, `m̂` the axis and `θ̂` increasing θ (counter-clockwise about `m̂`), cylindrical
/// coordinates are right-handed so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The face's outward is `σ·r̂`
/// with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], travel is `τ·m̂`, and left of travel is
///
/// ```text
///   n_out × travel = (σ·r̂) × (τ·m̂) = σ·τ·(r̂ × m̂) = −σ·τ·θ̂
/// ```
///
/// so the answer is `−σ·τ`: the θ direction the face occupies. Two readers need it — the ⊥ road,
/// where a ruling edge crossing the circle puts the **hole** in the opposite θ direction, and the
/// ∥ road, where a ruling edge lying on the class puts the **face** to one side of the wall — and
/// they must not spell it twice.
fn material_theta_sign(orient_sign: i8, travel_up: i8) -> i8 {
    -(orient_sign * travel_up)
}

/// **How the class's rational name is oriented against its stored normal** — `+1` when
/// [`combinatorics::class_coeffs_rat`] points the same way as `jd.planes[c].plane`, `-1` when it
/// opposes.
///
/// ★★★ **`world_rat` is a *name*, not an oriented normal** — it may be any nonzero multiple of the
/// stored one, negative included. A predicate built on it answers about *identity* (which of two
/// rulings, which side of a pair) frame-freely, because the same spelling is used on both sides of
/// the comparison; a **label** is different, because "above" is defined by the stored normal. This
/// is the correction [`combinatorics::side_of`]'s branch arm makes inline, lifted so the ∥ ruling
/// road can make the same one without spelling it a second time.
/// **Which way `+θ̂` points across a ∥ wall, in the frame a label is written in** — `+1` when
/// leaving a ruling counter-clockwise about the axis enters the wall's **stored-normal** side.
///
/// ★★★ **One atom, three readers.** [`combinatorics::RulingCarrier::side`] is
/// `sign((x − o) · (m̂ × n̂_r))` against the class's *rational* name, and the scalar triple product
/// gives `(x − o) · (m̂ × n̂) = n̂ · ((x − o) × m̂) = −r·(n̂ · θ̂)`, so `sign(n̂_r · θ̂) = −side`;
/// [`world_rat_sense`] (`κ`) carries that to the stored normal. The three consumers are the
/// ruling sweep's graze side, [`ruling_interior_is_even`], and the chart's vertical read — and the
/// whole point of naming it is that none of them spells the product a second time. ★ Two of them
/// read a **label** and stay in the stored frame; the one that picks a **cell** —
/// [`ruling_interior_is_even`] — crosses into the chart's frame with `frame_sign`, and it is the
/// only reader that does.
pub(crate) fn plus_theta_is_above(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    side: i8,
) -> Option<bool> {
    // A tangent ruling (`side == 0`) is never a recorded crossing, so nothing labelled by this
    // sign reaches it; the day one does, the answer is a measurement, not a sign.
    debug_assert_ne!(side, 0, "a tangent ruling has no +θ side to be above");
    Some(-(world_rat_sense(jd, wc)? * side) == 1)
}

fn world_rat_sense(jd: &Judge<'_, WorkingPlane>, c: usize) -> Option<i8> {
    let co = combinatorics::class_coeffs_rat(jd, c)?;
    let raw = jd.planes[c].plane.coefficients();
    let zero = nacre_scalar::Rat::from_int(0);
    // Both must be nonzero, not just the rational one: they are proportional so their zero sets
    // agree exactly, but `raw` is `f64` and a component it rounds to zero would hand back a sign
    // with nothing behind it.
    let i = (0..4).find(|&i| co[i] != zero && raw[i] != 0.0)?;
    Some(if (co[i] > zero) == (raw[i] > 0.0) {
        1
    } else {
        -1
    })
}

/// **One boundary cycle — a hole, a panel, a chain rim — read against one ⊥ class**: the walk's
/// features turned into extents. A hole is carved out of a whole-circle outer answer; a panel or a
/// chain is carved out of one that may be absent (E2-2) — the same arms, because neither reads
/// "inside" or "outside": only the ring's winding, below.
///
/// ★★★★★ **Which way round is decided by the ring's own winding, and by nothing else.** A plane
/// cuts a circle in two points, and "which of the two arcs is the hole" is the whole difficulty:
/// a hole whose two ruling edges lie on **one** plane cannot be told apart by that plane's sides.
/// The universal boundary convention answers it locally instead — *material is on the left of the
/// ring's direction of travel* — the same sentence [`run_body_above`] states for a planar face's
/// on-line run, and it holds for an outer ring, a hole ring, a notch and a reflex corner alike.
///
/// **Derivation.** Write `r̂` for the cylinder's outward radial direction, `m̂` for its axis and `θ̂`
/// for increasing θ (counter-clockwise about `m̂`, which is what
/// `nacre_scalar::quad::circular_order_about_seam` ranks and what [`MergedArc::end`] runs along).
/// Cylindrical coordinates are right-handed, so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The face's outward
/// is `σ·r̂` with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], and "left of travel" is
/// `n_out × travel`.
///
/// - **A ruling edge crossing the class** travels `τ·m̂`. Material lies along
///   `(σ·r̂) × (τ·m̂) = −σ·τ·θ̂`, so the **hole** lies along `+σ·τ·θ̂`: the hole runs counter-clockwise
///   from the crossing exactly when `σ·τ = +1`.
/// - **An arc edge lying on the class** (one of the hole's own rims) travels `ν·θ̂`. Material lies
///   along `(σ·r̂) × (ν·θ̂) = σ·ν·m̂`, so `ν = σ·μ` where `μ` is the axis direction the face occupies
///   — which is the opposite of the side the run's flanks sit on. That fixes the rim arc's CCW
///   orientation without reading a coordinate.
///
/// `up` is `plus_t_is_above`: whether the class's **stored** normal points along `+m̂`. It is the
/// bridge between the two frames here, because [`combinatorics::side_of`] answers in the class
/// root's **outward** frame and a label's "above" is the stored one — the correction is
/// [`crate::planes::WorkingPlane::frame_sign`], and forgetting it is a silently mirrored answer.
#[allow(clippy::too_many_arguments)]
fn cycle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    cyl: usize,
    up: bool,
    kind: combinatorics::CycleKind,
    nr: &combinatorics::NamedRing,
    out: &mut Vec<Carved>,
) -> Result<(), DeclineKind> {
    let _ = kind; // the ledger's, under test
    // ★★ **`true` here is not a shrug — the meet on this road is a *circle*, and an arc can lie on
    // one.** A rim arc of the hole at this very axis parameter *is* the class's meet, so "not an
    // arc" would cut runs that are genuinely continuous. What decides it is the carrier, and the
    // `Run` arm below asks that directly (declining, rather than splitting, is this road's scope).
    let features = match combinatorics::ring_against_plane(jd, cyls, &nr.triples, wc, |_| {
        Some(combinatorics::EdgeMeet::On)
    }) {
        combinatorics::RingWalk::Met(f) => f,
        // A hole ring lying wholly in the class has no thickness to bound anything with, and a
        // node the walk cannot name is the same refusal every other consumer makes of it.
        combinatorics::RingWalk::AllOn | combinatorics::RingWalk::Unnameable => {
            return Err(DeclineKind::CylFaceHole);
        }
    };
    let ring = &nr.triples;
    let n = ring.len();
    let sigma = i32::from(cf.orient_sign);
    let fs = i32::from(jd.planes[wc].frame_sign);
    let axis_of = |stored_side: i32| if up { stored_side } else { -stored_side };
    // Crossings, each with the direction the hole runs in from it. Paired below, after the θ sort.
    let mut cuts: Vec<(NodeId, bool)> = Vec::new();
    for f in features {
        match f {
            combinatorics::Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } => {
                // A run whose flanks differ is a rim the hole continues *through* — it is both an
                // extent and a parity toggle, and no fixture makes one. Declined and counted
                // rather than guessed at.
                if flanks_differ {
                    return Err(DeclineKind::CylHoleFeature);
                }
                // A single on-line vertex with equal flanks is the hole's corner touching the
                // circle at a point. A point takes no extent away, so there is nothing to carve.
                if len < 2 {
                    continue;
                }
                // ★★★★★ **Which way the arc runs, and so which side the face is on, is the
                // producer's to say — not the flank's.** This used to read the flank (the side
                // the ring's off-class neighbours sit on) and place the face opposite it, which
                // is the ring's *interior* side — right for a convex hole and wrong for a
                // wrapping rim, whose highest arc has both neighbours below it and the face
                // below too. The arc's own direction about the axis (`NamedRing::arc_ccw`, the
                // producer's convention) settles both: walked with sense `ν`, material lies
                // along `σ·ν·m̂`, so the face is above the class exactly when that agrees with
                // the class's stored normal (`up`), and the carved extent is the arc as walked.
                // ☑ E2-0 measured the two readings equal on every on-class arc of today's holes
                // (50 of 50) before this rule replaced the flank's.
                for k in 0..len - 1 {
                    // ★★★★★ **"Both ends on the class" is not "the edge is on the circle."** The
                    // walk answers about *nodes*; on a plane the edge between two on-line nodes
                    // follows, and on a **cylinder** it does not — a tilted carrier meets the
                    // lateral in an ellipse, which can cross this class at both ends without
                    // lying on it. Taking that for a rim arc would state an extent along a curve
                    // that is not there. The carrier says it directly: an edge on the circle
                    // `wc ∩ cylinder` lies in `wc`, so its far face is of that class.
                    // ★ This is the check `lateral_inner_loops` used to make from the row
                    // ("neither ⊥ nor ∥: an ellipse, outside this vocabulary") and that went with
                    // it; it belongs here, per edge, where the answer is actually used.
                    let edge = (first + k) % n;
                    if nr.walls[edge] != crate::boolean::Wall::Plane(wc) {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    // An arc on the class with no stated sense is not a lateral's arc at all.
                    let Some(nu) = nr.arc_ccw[edge] else {
                        return Err(DeclineKind::CylHoleFeature);
                    };
                    let body_above = (sigma * if nu { 1 } else { -1 } > 0) == up;
                    let a = ring[edge];
                    let b = ring[(edge + 1) % n];
                    let arc = if nu { [a, b] } else { [b, a] };
                    #[cfg(test)]
                    arc_probe::record(jd, cyl, def, kind, arc);
                    out.push(Carved {
                        arc,
                        on: Some(CylOnClass::Grazes { body_above }),
                    });
                }
            }
            combinatorics::Feature::Crossing { edge, from } => {
                // The crossed edge's carrier, **carried** from the producer rather than re-derived
                // from the two endpoint names.
                let crate::boolean::Wall::Plane(j) = nr.walls[edge] else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                let start = ring[edge];
                let Some((planes, ncyl, _)) = combinatorics::branch_name(start) else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                if !planes.contains(&j) {
                    return Err(DeclineKind::CylHoleFeature);
                }
                // ★★★★★ **The crossing is the *same ruling*, restated for this class — through
                // the one door the planar scan uses** ([`crossing_on_ruling`]): the ruling's side,
                // measured against the wall `j` at the hole's own corner ([`node_ruling_side`],
                // the predicate `curved_wall` carries on a ruling edge), picks the root of
                // `{wc, j}` on the cylinder. Its previous spelling restated the corner's *root* by
                // the axis senses of the two ⊥ classes (`ε·sign(k)`, the order of the roots along
                // `ℓ`) — the same point derived the other way round, and two spellings of one
                // rule were one too many. `j` must hold a ruling (a wall plane within the radius)
                // and `wc` be ⊥, which the door checks; anything else is not this feature.
                let cj = combinatorics::class_coeffs_rat(jd, j).ok_or(DeclineKind::CylSpan)?;
                let wdef = &cyls.get(ncyl).ok_or(DeclineKind::CylHoleFeature)?.def;
                let side =
                    node_ruling_side(jd, wdef, &cj, start).ok_or(DeclineKind::CylHoleFeature)?;
                let cut = crossing_on_ruling(jd, wdef, wc, j, ncyl, side)
                    .map_err(|_| DeclineKind::CylHoleFeature)?;
                // Which way the hole runs from here: the travel's axis sense, then the winding.
                // The side the edge leaves comes from the walk, for the same reason the run's
                // flank does.
                let s0 = i32::from(from);
                let tau = axis_of(if s0 * fs < 0 { 1 } else { -1 });
                // The face occupies `material_theta_sign`; the hole is the other way.
                let hole_ccw = material_theta_sign(cf.orient_sign, tau as i8) < 0;
                cuts.push((cut, hole_ccw));
            }
        }
    }
    // A closed ring enters and leaves the circle equally often.
    if cuts.len() % 2 != 0 {
        return Err(DeclineKind::CylHoleFeature);
    }
    // The crossings alternate in θ, so sorting them and pairing each "hole ahead" with the next
    // boundary is the whole assignment — the circle twin of the parity sweep
    // `trace_transversal_face` runs along a line.
    if !cuts.is_empty() {
        let nodes: Vec<NodeId> = cuts.iter().map(|c| c.0).collect();
        let (order, _) = circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
            CircleOrderFail::Undecided => DeclineKind::CylSpan,
            CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
        })?;
        for k in 0..order.len() {
            let (a, ahead) = cuts[order[k]];
            let (b, next_ahead) = cuts[order[(k + 1) % order.len()]];
            if ahead == next_ahead {
                return Err(DeclineKind::CylHoleFeature); // the boundaries did not alternate
            }
            if ahead {
                let arc = [a, b];
                #[cfg(test)]
                arc_probe::record(jd, cyl, def, kind, arc);
                out.push(Carved { arc, on: None });
            }
        }
    }
    Ok(())
}

/// Why a circle's chord on a class could not be stated — [`chord_nodes`]' two refusals, which its
/// two callers read differently: a **circular outer** falls silent on missing coefficients (the
/// disk's old silence for a class it cannot read) and declines the rest, a **circular hole**
/// declines both.
enum ChordFail {
    /// The class `wc` has no rational world coefficients (a rotated class).
    NoCoefficients,
    /// No cylinder statement, no coefficients for the face's own class, or the meet is not a pair
    /// of roots — a piece this road cannot state.
    Unstatable,
}

/// **The two ends of the chord a recorded wall class `wc` cuts on a circle of a planar face** (the
/// face lies in class `fc`, the circle rides cylinder `cyl`), as flip nodes of
/// [`trace_transversal_face`]'s parity sweep — a diameter when the wall runs through the axis,
/// any chord within the radius otherwise (cell ③). `Ok(None)` for a pair the gate did not record:
/// the wall clears the circle (the gate's plane test), or cuts an **irrational** chord this road
/// cannot state yet — today's silence, and the same trigger the rulings road reads.
///
/// ★ **Written once because a disk's outer and a bored face's hole are the same circle asked the
/// same question** — each used to have its own spelling (a `ChordTrace` vessel for the outer, an
/// inline pair of nodes for the hole), and the outer's returned before the holes were read.
///
/// The cylinder's statement comes from any face row on its class — the same table `merge_circles`
/// indexes — so a caller that traces with no cylinder table (the test shims) is served too.
fn chord_nodes(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    wc: usize,
    fc: usize,
    cyl: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
) -> Result<Option<[Node; 2]>, ChordFail> {
    use nacre_scalar::quad::CylinderMeet;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(ChordFail::NoCoefficients)?;
    if !crossings.contains(&(wc, cyl)) {
        return Ok(None);
    }
    let def = faces
        .iter()
        .zip(plane_ix)
        .find_map(|(row, ix)| match (row, ix) {
            (FaceRow::Cylinder(cf), ClassIx::Cyl(k)) if *k == cyl => cf.def.as_ref(),
            _ => None,
        })
        .ok_or(ChordFail::Unstatable)?;
    // ★ The record is the crossing statement (cell ③): a listed pair's plane runs within the
    // radius, so `L` enters the circle at one root and leaves at the other. The premise is the
    // gate's, asserted rather than re-derived. (A *tangent* plane has one root and is never
    // listed; see `planes::Tangency`.)
    debug_assert_eq!(
        nacre_scalar::point_plane_clearance_rat(&w, &def.origin(), def.radius()),
        nacre_scalar::Orient::Negative,
        "a recorded pair's plane runs within the radius"
    );
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(ChordFail::Unstatable)?;
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let Some(CylinderMeet::Pair { .. }) =
        nacre_scalar::quad::plane_plane_cylinder(&w, &v, &o, &m, r)
    else {
        return Err(ChordFail::Unstatable);
    };
    let node = |root| Node {
        id: NodeId::branch(wc, fc, cyl, root),
        pin: combinatorics::EndPin::Cylinder,
        flip: true,
        run: None,
        flanks_differ: false,
        single_touch: false,
        graze_above: None,
    };
    Ok(Some([
        node(nacre_topo::QuadRoot::Lo),
        node(nacre_topo::QuadRoot::Hi),
    ]))
}

/// Which of the two rulings a point of the lateral lies on: the sign of `(x − o) · (m × n̂)`,
/// with `n̂` the class's **canonical** coefficients ([`combinatorics::RulingCarrier::side`] — the
/// one spelling). The point arrives as `(line, s)` from [`nacre_scalar::quad::plane_plane_cylinder`],
/// so the sign is one [`nacre_scalar::quad::plane_side`] against the plane through `o` with
/// normal `m × n̂`. `None`: overflow, or the point is on the axis plane itself (no side — the
/// tangent shape). ★ **A tangent ruling's side is `0`, and it is not read here** (cell ⑩, S3): the
/// two places that build a ruling carrier ([`combinatorics::curved_wall`], [`node_ruling_side`])
/// derive it from the corner's root (`Double`), because this function's `None` is a *sign*
/// fact the ray caster (`departs_across`) and the arc departure read as "abstain" — answering
/// `Some(0)` here would drop a ray into the `Negative` arm in silence. A fillet's or a slot's
/// own walls are the callers that feed the tangent shape.
pub(crate) fn ruling_side(
    w: &[Rat; 4],
    def: &nacre_topo::CylinderDef,
    at: (&nacre_scalar::quad::MeetLine, &nacre_scalar::quad::QuadVal),
) -> Option<i8> {
    ruling_side_signed(w, def, at).filter(|&s| s != 0)
}

/// **Which ruling of `w` a point of the lateral is on, `0` being `w`'s tangent ruling** — the
/// same sign as [`ruling_side`], with the axis plane answered rather than abstained on.
///
/// ★★★★★ Cell ⑫ — **a side is a fact about the point and the wall, not about the point's own
/// name.** Two sites used to read it off the name's root instead (`QuadRoot::Double ⇒ 0`:
/// [`node_ruling_side`] and [`combinatorics::curved_wall`]), which is the same statement *only*
/// when the wall asked about is the very wall the name pairs — the tangent wall the corner was
/// minted on. It is not the same statement once the alias table can hand back a representative
/// with another pair (a tangent corner a class through the axis also passes through, this cell's
/// whole subject) or once a Double-rooted corner is asked about a *different* wall: then the
/// root says `0` for a point that sits squarely on one of that wall's two rulings — silently, and
/// with the sign the caller needs. Asked here, of the point, both cases answer correctly and the
/// tangent wall still answers `0`.
///
/// `None` is checked-`Rat` overflow only.
pub(crate) fn ruling_side_signed(
    w: &[Rat; 4],
    def: &nacre_topo::CylinderDef,
    at: (&nacre_scalar::quad::MeetLine, &nacre_scalar::quad::QuadVal),
) -> Option<i8> {
    let (o, m) = (def.origin(), def.dir());
    let n = [w[0], w[1], w[2]];
    let c = combinatorics::cross3_rat(&m, &n)?;
    let mut d = Rat::from_int(0);
    for k in 0..3 {
        d = d.checked_sub(c[k].checked_mul(o[k])?)?;
    }
    Some(
        match nacre_scalar::quad::plane_side(&[c[0], c[1], c[2], d], at.0, at.1) {
            nacre_scalar::Orient::Positive => 1,
            nacre_scalar::Orient::Negative => -1,
            nacre_scalar::Orient::Zero => 0,
        },
    )
}

/// **The planar scan's crossing on a ruling** — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a cylinder's ruling, named as the branch node
/// `wc ∩ fc ∩ cyl` at the root that *is* this ruling.
///
/// A ruling edge lies in the face's own plane `fc` (a plane holding a ruling runs through the
/// axis — or is tangent, which the gate passes but does **not** record, so no ruling of it ever
/// reaches here), so the pair `{wc, fc}` cuts the cylinder in two
/// points, one on each of `fc`'s two rulings, and `(cyl, side)` — the identity
/// [`crate::boolean::Wall::Ruling`] carries, measured by [`ruling_side`] against `fc` when the
/// ring was named — says which. The same predicate asked of each root picks it. `wc` must be ⊥
/// to the axis for the class to cross a ruling in a point at all (∥ contains it; a tilt is
/// refused at the gate).
///
/// ★ Solved in the caller's order `(wc, fc)`, as [`rulings_on_class`] and [`chord_on_class`]
/// spell it — [`NodeId::branch`] canonicalizes the pair and the root together, so this is the one
/// name every road gives the point. The lateral face names the same point when its hole ring
/// crosses the class ([`cycle_on_class`]'s crossing arm restates the hole corner's own root to
/// `wc` by the axis senses of the two ⊥ classes); the two spellings agree because the roots of
/// `{⊥, wall}` are ordered along `ε·k·(m̂ × n_wall)`, so «which root» and «which side of `wall`»
/// are the same question — a derivation the exact-volume rows of the crossing census check.
///
/// `Err(CurvedRingWall)` is every shape this does not state: no exact description, a wall that
/// is not parallel to the axis, no pair of roots, or both roots on one side (the two rulings of
/// a wall within the radius — through the axis or offset from it, cell ③ — are symmetric about
/// the plane through the axis with normal `m × n̂`, so `ruling_side` tells them apart).
pub(crate) fn crossing_on_ruling(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    side: i8,
) -> Result<NodeId, DeclineKind> {
    use nacre_scalar::quad::CylinderMeet;
    let no = DeclineKind::CurvedRingWall;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(no)?;
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(no)?;
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    if !nacre_scalar::parallel_rat(&[w[0], w[1], w[2]], &m) {
        return Err(no);
    }
    // ★ A **tangent** wall (`side == 0`, cell ⑩) has one ruling and one root — `Double`.
    let meet = nacre_scalar::quad::plane_plane_cylinder(&w, &v, &o, &m, r);
    if side == 0 {
        return match meet {
            Some(CylinderMeet::Tangent { .. }) => {
                Ok(NodeId::branch(wc, fc, cyl, nacre_topo::QuadRoot::Double))
            }
            _ => Err(no),
        };
    }
    let Some(CylinderMeet::Pair { line, s }) = meet else {
        return Err(no);
    };
    let mut found = None;
    for (root, sv) in [
        (nacre_topo::QuadRoot::Lo, &s[0]),
        (nacre_topo::QuadRoot::Hi, &s[1]),
    ] {
        if ruling_side(&v, def, (&line, sv)) == Some(side) {
            if found.is_some() {
                return Err(no); // both roots on one side: not a pair of rulings
            }
            found = Some(root);
        }
    }
    Ok(NodeId::branch(wc, fc, cyl, found.ok_or(no)?))
}

/// **The planar scan's crossing on an arc** (E3-b) — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a circle of `cyl`, named as the branch node
/// `wc ∩ fc ∩ cyl` at the root that lies **inside** the arc.
///
/// The arc lies in the face's own plane `fc` (a cap's ⊥ plane), so the pair `{wc, fc}` cuts the
/// cylinder in two points on the arc's circle, and the crossing is the one the travelled arc
/// `a → b` strictly contains ([`theta_between`], the ruling sweep's containment — the arc's CCW
/// pair is `ccw ? [a, b] : [b, a]`). Exactly one must: none is a crossing the walk mis-read, and
/// both is an arc meeting the line twice — a shape the gate keeps out (a > π arc against its own
/// diameter's class), refused rather than guessed.
///
/// ★ The same canonical name the lateral's ruling sweep gives this point as a station
/// (`crossing_on_ruling(fc, wc, …)` — [`NodeId::branch`] folds the pair and root together), so
/// [`merge_coincident`] reads one point, not two.
#[allow(clippy::too_many_arguments)]
pub(crate) fn crossing_on_arc(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    ccw: bool,
    a: NodeId,
    b: NodeId,
) -> Result<NodeId, DeclineKind> {
    use nacre_scalar::quad::CylinderMeet;
    let no = DeclineKind::CurvedRingWall;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(no)?;
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(no)?;
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let Some(CylinderMeet::Pair { s, .. }) =
        nacre_scalar::quad::plane_plane_cylinder(&w, &v, &o, &m, r)
    else {
        return Err(no);
    };
    let _ = s;
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let mut found = None;
    for root in [nacre_topo::QuadRoot::Lo, nacre_topo::QuadRoot::Hi] {
        let node = NodeId::branch(wc, fc, cyl, root);
        if theta_between(jd, cyl, def, lo, hi, node).map_err(|_| no)? {
            if found.is_some() {
                return Err(no); // both roots inside: the arc meets the line twice
            }
            found = Some(node);
        }
    }
    found.ok_or(no)
}

/// The scan's crossings on rulings, realized — the second road to the sign
/// [`crossing_on_ruling`] chooses by. One entry per crossing named in this binary: the point,
/// how far it sits from the class plane, the face plane and the cylinder (all should be 0), and
/// the ruling side read off the realization beside the side the ring carried.
/// **The disk-side rule's premise, watched where the rule is applied** (`emit_faces`' arc labels).
///
/// A cut circle's per-arc label must be the cell on the arc's **disk** side, and that side is read
/// off the cells themselves: a cell cannot straddle the circle — the circle is an arrangement edge
/// — so any rational corner of a cell with a definite radial side names the side the whole cell is
/// on. This watches the premise rather than the conclusion: **no cell has corners on both sides**,
/// and **an arc's two cells are never on the same side**. Both would make the answer a coin toss,
/// and neither is checkable from the label afterwards (a holed lateral's two sectors can carry a
/// literally identical label — see [`ArcLabel`]).
///
/// ★ It replaced a rule derived from the class's **stored** frame (`axis_up`), which two classes
/// with identical stored *and* canonical normals were measured to disagree about — the same plane,
/// the same circle, the same two arcs, opposite half-edges. That rule had been set by measuring
/// two fixtures and generalising; this one asks the geometry every time.
/// **Which side of a circle a cell lies on, from its own corners** — the geometry that watches
/// the disk-side rule. A cell cannot straddle the circle (the circle is an arrangement edge), so
/// any rational corner with a definite radial side names the side the whole cell is on; a branch
/// corner sits *on* the circle and says nothing, and so does a three-plane corner that happens to
/// land there.
#[cfg(test)]
fn corner_sides<'a>(
    jd: &'a Judge<'_, WorkingPlane>,
    edges: &'a ClassEdges<'_>,
    cell: &'a Cell,
    def: &'a nacre_topo::CylinderDef,
) -> impl Iterator<Item = bool> + 'a {
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    cell.half_edges.iter().filter_map(move |&h| {
        let p = combinatorics::node_coords_rat(jd, edges.origin(h))?;
        match nacre_scalar::cylinder_radial_side(&p, &o, &m, r) {
            nacre_scalar::Orient::Negative => Some(true),
            nacre_scalar::Orient::Positive => Some(false),
            nacre_scalar::Orient::Zero => None,
        }
    })
}

/// **Which side of a circle a cell lies on, by its rational corners — when they agree.** The
/// instrument used to take the *first* corner, on the premise that no cell has corners on both
/// sides of a circle it borders. ★ A fillet refutes the premise (cell ⑩): the cap face's cell is
/// bounded by a **quarter** of the circle and reaches far beyond it, so its corners lie outside
/// while the arc bounds it from the disk side. Such a cell has no single side and the witness
/// abstains; a cell all of whose corners agree answers as before.
#[cfg(test)]
fn cell_side(
    jd: &Judge<'_, WorkingPlane>,
    edges: &ClassEdges<'_>,
    cell: &Cell,
    def: &nacre_topo::CylinderDef,
) -> Option<bool> {
    let mut sides = corner_sides(jd, edges, cell, def);
    let first = sides.next()?;
    sides.all(|s| s == first).then_some(first)
}

#[cfg(test)]
pub(crate) mod disk_side_probe {
    use std::sync::Mutex;

    /// The sentence the disk-side check panics with — one spelling, shared with the commuting
    /// oracle's `KNOWN` list (cell ④).
    pub(crate) const NOT_OWN_SOLID: &str = "a cut circle's disk side does not carry its own solid";

    /// One entry per arc.
    #[derive(Clone, Copy, Debug, Default)]
    pub(crate) struct Row {
        /// A corner named the side, so the rule was checked against the geometry here.
        pub(crate) checked: bool,
        /// Both sides spoke — the strongest form of the check.
        pub(crate) both_spoke: bool,
        /// The class's `frame_sign` was `-1`, so the factor decided this arc.
        pub(crate) frame_negative: bool,
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    /// Assert the premises **and** the rule, and record the row. Each side arrives as its
    /// **first** corner that names a side (`None` = no corner spoke, or the cell was not found).
    ///
    /// Two propositions, where the fact is made: an arc's two cells are never on the same side of
    /// its circle, and where the geometry speaks it **agrees with the derived rule**. The second
    /// is what keeps the derivation honest — the rule was a guess once, and the guess was wrong
    /// by one factor.
    ///
    /// ★ A third premise — *no cell has corners on both sides of a circle it cannot straddle* —
    /// is what licenses taking the **first** corner as the cell's answer. Scanning every corner
    /// to assert it costs a rational three-plane solve per corner per arc: the lib suite measured
    /// **526 s** that way against **58 s** taking the first (one session, warm builds). It was
    /// measured **0 violations** over the whole suite twice — when the geometric road was built
    /// and again here — so it is recorded rather than re-proved on every run.
    pub(crate) fn record(
        even: Option<bool>,
        odd: Option<bool>,
        rule_says_even: bool,
        frame_sign: i8,
    ) {
        let (a, b) = (even, odd);
        // ★ **Two cells on one side is the witness failing, not the arc** (cell ⑩): a **convex**
        // arc — a fillet's quarter, a slot's end — bounds a cell that lies inside the circle at
        // the arc and reaches far outside it, so every rational corner of that cell is outside
        // while the arc bounds it from the disk side. The corner witness cannot see a side
        // there; where the two cells' corners agree, both abstain and the rule goes unchecked
        // for that arc (the volume oracles of the tangent fixtures are what measure it). The
        // bite population — arcs concave into a plate — keeps its witness exactly as before.
        let (a, b) = match (a, b) {
            (Some(x), Some(y)) if x == y => (None, None),
            other => other,
        };
        // The geometry's verdict on which half-edge borders the disk, where it has one.
        let witness = a.or_else(|| b.map(|inside| !inside));
        if let Some(even_is_disk) = witness {
            assert_eq!(
                even_is_disk, rule_says_even,
                "the disk-side rule and the cell's own corners disagree (frame_sign {frame_sign}, even {a:?}, odd {b:?})"
            );
        }
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(Row {
                checked: witness.is_some(),
                both_spoke: a.is_some() && b.is_some(),
                frame_negative: frame_sign < 0,
            });
    }
}

pub(crate) mod crossing_probe {
    use super::{Judge, NodeId, WorkingPlane, combinatorics};
    use std::sync::Mutex;

    #[derive(Clone, Copy, Debug)]
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) struct Hit {
        /// Read through `Debug` in the probe's messages.
        #[allow(dead_code)]
        pub point: [f64; 3],
        /// Distances to the class plane, the face plane, and the cylinder's surface.
        pub off: [f64; 3],
        /// `sign((x − o) · (m̂ × n_fc))` of the realized point — the f64 twin of
        /// [`super::ruling_side`].
        pub side_f64: i8,
        /// The side the ring's edge carried.
        pub side: i8,
        /// `true` for a crossing on an arc (E3-b), whose carried side is 0 — the ruling-side
        /// twin check is the rulings' alone.
        pub arc: bool,
    }

    pub(crate) static HITS: Mutex<Vec<Hit>> = Mutex::new(Vec::new());

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn record(
        jd: &Judge<'_, WorkingPlane>,
        def: &nacre_topo::CylinderDef,
        cyl: usize,
        wc: usize,
        fc: usize,
        side: i8,
        arc: bool,
        id: NodeId,
    ) {
        let Some(p) = combinatorics::branch_point(jd, cyl, def, id) else {
            return;
        };
        let coeffs = |c: usize| -> Option<[f64; 4]> {
            combinatorics::class_coeffs_rat(jd, c).map(|w| w.map(|x| x.to_f64()))
        };
        let (Some(w), Some(v)) = (coeffs(wc), coeffs(fc)) else {
            return;
        };
        let plane_off = |w: [f64; 4]| -> f64 {
            let n = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
            (w[0] * p[0] + w[1] * p[1] + w[2] * p[2] + w[3]).abs() / n
        };
        let o = def.origin().map(|x| x.to_f64());
        let raw = def.dir().map(|x| x.to_f64());
        let ml = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
        let m: [f64; 3] = core::array::from_fn(|i| raw[i] / ml);
        let d: [f64; 3] = core::array::from_fn(|i| p[i] - o[i]);
        let h = d[0] * m[0] + d[1] * m[1] + d[2] * m[2];
        let perp: [f64; 3] = core::array::from_fn(|i| d[i] - h * m[i]);
        let rho = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
        let cyl_off = (rho - def.radius().to_f64()).abs();
        let c = [
            m[1] * v[2] - m[2] * v[1],
            m[2] * v[0] - m[0] * v[2],
            m[0] * v[1] - m[1] * v[0],
        ];
        let dot = c[0] * d[0] + c[1] * d[1] + c[2] * d[2];
        let side_f64 = if dot > 0.0 {
            1
        } else if dot < 0.0 {
            -1
        } else {
            0
        };
        HITS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(Hit {
                point: p,
                off: [plane_off(w), plane_off(v), cyl_off],
                side_f64,
                side,
                arc,
            });
    }
}

/// **What a lateral face leaves on a recorded wall class** — the ∥ sibling of
/// [`circle_on_class`], and like it, an answer **per extent** rather than one for the whole
/// ruling.
///
/// Asked only after the circle road answered "no circle"; **empty** (silently, today's state)
/// unless the gate **recorded** the pair (the wall within the radius — through the axis or
/// offset from it, cell ③), and declining rather than contributing
/// *partially* when it does but a piece cannot be stated — a half-contributed rectangle would
/// leave the class's 1-skeleton dangling, which is a worse lie than an honest incomplete-trace
/// mark.
///
/// The class cuts **two** rulings, and each is [`ruling_sweep`]'s: the face's boundary cycles
/// against the line `θ = θ_side` of the `(θ, z)` chart, swept along the axis — `Transversal`
/// where the solid straddles the wall, `Graze` where a cycle's own edge lies on the ruling, and
/// nothing where the face is not there (E2-2: a panel, a chain rim, a band with holes, all by the
/// one rule).
#[allow(clippy::too_many_arguments)]
fn rulings_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    k: usize,
    which: SolidSide,
    crossings: &std::collections::HashSet<(usize, usize)>,
    aliases: &Aliases,
) -> Result<(Vec<RulingTrace>, Vec<Seg>), DeclineKind> {
    // ★ **The gate's answer, first** — see [`combinatorics::TraceInput::crossings`]. A pair not
    // listed there was proven clear (or never in question), and contributing its rectangle
    // anyway is what broke the straddling family: the boss's own wall class is a d=0 pair, and
    // its coplanar bottom cap's chord lands on a line the plate already traces. Silence for
    // unlisted pairs is today's exact behavior; the gate-opening cell lists exactly the pairs
    // that need the rectangle.
    if !crossings.contains(&(wc, k)) {
        return Ok((Vec::new(), Vec::new()));
    }
    let Some(w) = combinatorics::class_coeffs_rat(jd, wc) else {
        return Ok((Vec::new(), Vec::new())); // no exact description: the circle road's decline covers ⊥ classes
    };
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::Ruling);
    };
    // ★ The record is the crossing statement (cell ③): a listed pair's plane runs within the
    // radius, so the wall meets the lateral in two rulings; the sweep reads the face's cycles
    // and says where. The premise is the gate's, asserted rather than re-derived — and it survived
    // the tangent wall opening (cell ⑥) precisely because a tangency is written to `tangencies`
    // and **not** to `crossings`: "within" still means within.
    debug_assert_eq!(
        nacre_scalar::point_plane_clearance_rat(&w, &def.origin(), def.radius()),
        nacre_scalar::Orient::Negative,
        "a recorded pair's plane runs within the radius"
    );
    // The face's shape from its cycles — the same reading the ⊥ road makes.
    let LateralShape { rims, cycles, .. } = lateral_shape(jd, cf, fl, def)?;
    let mat = SegKind::Transversal {
        mat: cf.orient_sign,
    };
    let mut out = Vec::with_capacity(2);
    let mut grazes = Vec::new();
    for side in [1i8, -1i8] {
        let (pieces, g) = ruling_sweep(
            jd, cf, def, &w, wc, k, side, &rims, &cycles, mat, which, aliases,
        )?;
        grazes.extend(g);
        #[cfg(test)]
        if pieces
            .iter()
            .any(|(_, kind)| matches!(kind, SegKind::Graze { .. }))
        {
            let o = def.origin().map(|x| x.to_f64());
            let at = |n: NodeId| node_axis_param(jd, def, n).map_or(f64::NAN, |t| t.to_f64());
            ruling_probe::CARVED
                .lock()
                .expect("the probe's lock is never held across a panic")
                .push(ruling_probe::Carved {
                    origin: o,
                    span: [at(pieces[0].0[0]), at(pieces[pieces.len() - 1].0[1])],
                    kinds: pieces.iter().map(|&(_, k)| k).collect(),
                });
        }
        for (end, kind) in pieces {
            out.push(RulingTrace {
                cyl: k,
                side,
                end,
                solid: which,
                kind,
                #[cfg(test)]
                orient: cf.orient_sign,
            });
        }
    }
    Ok((out, grazes))
}

/// The θ side a cycle's arc lies on, seen from the node where it meets a ruling: an arc arriving
/// counter-clockwise came from smaller θ (`−1`), one departing counter-clockwise goes to larger θ
/// (`+1`) — the flank of an on-line run, read from the arc's own sense (`arc_ccw`) rather than
/// from a coordinate, because a global "side" does not exist on a cylinder.
fn arc_flank(ccw: bool, arriving: bool) -> i8 {
    if ccw == arriving { -1 } else { 1 }
}

/// **Is `node` strictly inside the counter-clockwise arc from `lo` to `hi`?** — cyclic order about
/// the axis, by name ([`circular_order`]). `Err` when the order cannot be formed or two of the
/// three are one point wearing two names — the caller's population, not a tie to shrug at.
fn theta_between(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    lo: NodeId,
    hi: NodeId,
    node: NodeId,
) -> Result<bool, DeclineKind> {
    let nodes = [lo, hi, node];
    let (order, _) = circular_order(jd, k, def, &nodes).map_err(|_| DeclineKind::Ruling)?;
    let pos = |x: usize| {
        order
            .iter()
            .position(|&i| i == x)
            .expect("a permutation of the three")
    };
    let (p_lo, p_hi, p_n) = (pos(0), pos(1), pos(2));
    Ok((p_n + 3 - p_lo) % 3 < (p_hi + 3 - p_lo) % 3)
}

/// What one side's sweep states: the ruling's pieces (its own vocabulary) and, for a plane-pair
/// line the face touches on this ruling, the face's side of that line (the plane vocabulary).
type SweepOut = (Vec<([NodeId; 2], SegKind)>, Vec<Seg>);

/// **One ruling of a through-axis class against the face's boundary cycles, swept along the
/// axis** (E2-2) — the ∥ twin of the circle road's carving, and the polygon-against-a-line rule
/// of `ring_against_plane` stated on the `(θ, z)` chart for the line `θ = θ_side`.
///
/// Every cycle contributes **stations** on the ruling, each with an event:
/// - a whole **rim** toggles the face (it starts or ends there);
/// - an **arc** that strictly contains the ruling's θ ([`theta_between`]) toggles it — the boundary
///   crosses the line there, at the branch node `{arc's plane, wc}` ([`crossing_on_ruling`], the
///   one door every road names such a point through);
/// - a **run** of the cycle's own edges lying on the ruling is a `Graze` over its extent, the side
///   the face occupies along it `σ·τ·κ·side` (the rulings ladder's derivation, unchanged), and it
///   toggles the face exactly when the arcs at its two ends lie on different θ flanks
///   ([`arc_flank`]) — a collinear piece the boundary passes *through* rather than touches, the
///   `flanks_differ` of an on-line run. A panel's arcs both lie on its own side (no toggle: absent
///   above and below); a chain's arrive from the outer half and leave into the inner (toggle:
///   transversal below, absent above); a band's notch hole has both on the hole's side (no toggle).
///
/// Then the sweep: absent below the first station, a piece per gap between stations — `Graze`
/// inside a run, `Transversal` while the face is present, nothing while it is not — and adjacent
/// equal pieces merged. Declined by name: a run whose neighbours are not arcs, an arc ending on the
/// ruling where no run of its cycle sits (the boundary would cross at a vertex — a shape
/// `loop_triples` does not produce), two names at one station, a run inside a run, a toggle inside
/// a run (the successor of "a graze past a rim"), and a face still present past the last station.
#[allow(clippy::too_many_arguments)]
fn ruling_sweep(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    w: &[Rat; 4],
    wc: usize,
    k: usize,
    side: i8,
    rims: &[(Rat, usize)],
    cycles: &[(combinatorics::CycleKind, combinatorics::LoopRing)],
    mat: SegKind,
    which: SolidSide,
    aliases: &Aliases,
) -> Result<SweepOut, DeclineKind> {
    use crate::boolean::Wall;
    // ★ Cell ⑫: **the lateral's half of a plane-pair line.** An edge of this cycle can lie on
    // this ruling with the face *across* it being another plane `t` — the fillet's tangent edge,
    // when `wc` passes through the axis and so contains the tangent ruling. That line is a
    // plane-pair line (`wc ∩ t`), and the carrier rule says the plane vocabulary states it: the
    // tangent wall's own run already does (`Feature::Run` → `Graze`), from *its* side. But the
    // solid's material near that line lies on **both** sides of `wc` — the tangent wall's below,
    // the fillet face's above — and an edge's mask is the sum of its faces' statements
    // (`edge_mask`: two grazes of one solid, one per side, flip both bits). So the lateral states
    // its side too — as a `Seg` on that line, in the plane vocabulary, so `merge_coincident`
    // joins it with the wall's — and states nothing in the ruling vocabulary (`quiet_nodes`).
    // The side is read exactly as a run's on this ruling is (`body_side`): the same cycle, the
    // same ruling, the same traversal.
    let mut grazes: Vec<Seg> = Vec::new();
    // ★ Cell ⑫: **identity is the table's** — a station made here (`wc`'s crossing with a cap)
    // and a corner the cycle carries can be one point under two names, and the plane roads have
    // already told the table so (`Aliases::record_on_cylinder`). Every "same node?" below asks
    // through it; the sweep never decides coincidence on its own.
    let same = |x: NodeId, y: NodeId| aliases.canon_point(x) == aliases.canon_point(y);
    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum Event {
        // Ordered as they apply at one station: a run ends, the face toggles, a run starts.
        GrazeEnd { toggle: bool },
        Toggle,
        GrazeStart { body_above: bool },
    }
    let on_ruling =
        |c: usize| crossing_on_ruling(jd, def, c, wc, k, side).map_err(|_| DeclineKind::Ruling);
    let mut stations: Vec<(Rat, Event, NodeId)> = Vec::new();
    for &(t, c) in rims {
        stations.push((t, Event::Toggle, on_ruling(c)?));
    }
    // The side of `wc` the face's body lies on along an edge of this cycle that runs on the
    // ruling from `start` to `end` — see [`plus_theta_is_above`], the one spelling of the product.
    let body_side = |start: NodeId, end: NodeId| -> Result<(Rat, Rat, i8, bool), DeclineKind> {
        let (Some(ta), Some(tb)) = (
            node_axis_param(jd, def, start),
            node_axis_param(jd, def, end),
        ) else {
            return Err(DeclineKind::Ruling);
        };
        if ta == tb {
            return Err(DeclineKind::CylHoleFeature); // an extent of no length
        }
        let tau: i8 = if tb > ta { 1 } else { -1 };
        let theta_dot_stored: i8 =
            if plus_theta_is_above(jd, wc, side).ok_or(DeclineKind::Ruling)? {
                1
            } else {
                -1
            };
        let body_above =
            i32::from(material_theta_sign(cf.orient_sign, tau)) * i32::from(theta_dot_stored) > 0;
        Ok((ta, tb, tau, body_above))
    };
    for (kind, ring) in cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair.
        let Some(nr) = ring.poly() else {
            return Err(DeclineKind::CylFaceHole);
        };
        let n = nr.triples.len();
        #[cfg(not(test))]
        let _ = kind;
        // Which edges lie on this ruling. Two shapes, one rule — **the line's carrier says who
        // states it** (cell ⑫):
        // - a straight edge carried by `wc` itself (its far face is *seated* on the class): the
        //   line is the cylinder's alone, and this sweep states it as a run (`on`);
        // - a straight edge carried by another plane `t` whose two ends are this ruling's cap
        //   crossings (by name, through the table): the line is the plane pair `wc ∩ t`, and
        //   `t`'s own trace states it as a graze run on this class — this sweep stays **silent**
        //   there (`quiet`), and only remembers the ends so an arc ending on them is "an arc that
        //   ends on this ruling", not a crossing.
        let mut on = vec![false; n];
        let mut quiet_nodes: Vec<NodeId> = Vec::new();
        for (i, slot) in on.iter_mut().enumerate() {
            if nr.arc_ccw[i].is_some() {
                continue;
            }
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            match nr.walls[i] {
                Wall::Plane(t) if t == wc => {
                    let (Some(sa), Some(sb)) = (
                        node_ruling_side(jd, def, w, a),
                        node_ruling_side(jd, def, w, b),
                    ) else {
                        return Err(DeclineKind::Ruling);
                    };
                    // An edge whose two ends answer differently is not a ruling at all — it
                    // would have to cross the plane through the axis. Named rather than silently
                    // mis-placed.
                    if sa != sb {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    *slot = sa == side;
                }
                Wall::Plane(t) => {
                    // Both ends this ruling's crossings with the caps the ends name? Then the edge
                    // lies on the ruling and `wc ∩ t` carries it.
                    let on_ruling_end = |x: NodeId| -> bool {
                        combinatorics::branch_name(x).is_some_and(|(pl, _, _)| {
                            pl.iter().any(|&c| {
                                c != t
                                    && crossing_on_ruling(jd, def, c, wc, k, side)
                                        .is_ok_and(|id| same(id, x))
                            })
                        })
                    };
                    if on_ruling_end(a) && on_ruling_end(b) {
                        quiet_nodes.push(a);
                        quiet_nodes.push(b);
                        let (_, _, _, body_above) = body_side(a, b)?;
                        let pin = |x: NodeId| {
                            combinatorics::pin_for(jd, wc, t, x).ok_or(DeclineKind::Ruling)
                        };
                        grazes.push(Seg {
                            wall: t,
                            end: [a, b],
                            end_h: [pin(a)?, pin(b)?],
                            solid: which,
                            kind: SegKind::Graze { body_above },
                        });
                    }
                }
                _ => {}
            }
        }
        if n > 0 && on.iter().all(|&x| x) {
            return Err(DeclineKind::Ruling); // a cycle lying wholly on one ruling is no face
        }
        let mut run_nodes: Vec<NodeId> = Vec::new();
        for i0 in 0..n {
            if !on[i0] || on[(i0 + n - 1) % n] {
                continue; // not the first edge of a run
            }
            let mut i1 = i0;
            while on[(i1 + 1) % n] {
                i1 = (i1 + 1) % n;
            }
            let (start, end) = (nr.triples[i0], nr.triples[(i1 + 1) % n]);
            let (prev, next) = ((i0 + n - 1) % n, (i1 + 1) % n);
            let (Some(nu_in), Some(nu_out)) = (nr.arc_ccw[prev], nr.arc_ccw[next]) else {
                return Err(DeclineKind::Ruling); // a run's neighbours are arcs
            };
            let (ta, tb, tau, body_above) = body_side(start, end)?;
            let toggle = arc_flank(nu_in, true) != arc_flank(nu_out, false);
            let ((t_lo, n_lo), (t_hi, n_hi)) = if tau > 0 {
                ((ta, start), (tb, end))
            } else {
                ((tb, end), (ta, start))
            };
            #[cfg(test)]
            ruling_probe::GRAZE_SIDE
                .lock()
                .expect("the probe's lock is never held across a panic")
                .push(ruling_probe::GrazeSide {
                    kind: *kind,
                    body_above,
                    ny: jd.planes[wc].plane.normal().as_array()[1],
                    origin: def.origin().map(|x| x.to_f64()),
                    span: [t_lo.to_f64(), t_hi.to_f64()],
                });
            stations.push((t_lo, Event::GrazeStart { body_above }, n_lo));
            stations.push((t_hi, Event::GrazeEnd { toggle }, n_hi));
            run_nodes.push(start);
            run_nodes.push(end);
        }
        // Arcs strictly containing the ruling's θ: crossings of the line.
        for i in 0..n {
            let Some(nu) = nr.arc_ccw[i] else { continue };
            let Wall::Plane(c) = nr.walls[i] else {
                return Err(DeclineKind::Ruling);
            };
            let node = on_ruling(c)?;
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            if same(node, a) || same(node, b) {
                // The arc ends on this ruling: an edge of this cycle must lie there — a run this
                // sweep states, or a plane-pair line another face states (`quiet_nodes`) — or the
                // boundary crosses the ruling at a vertex: not a shape the ring producer makes,
                // and not one to answer "no" to silently.
                if !run_nodes
                    .iter()
                    .chain(quiet_nodes.iter())
                    .any(|&r| same(r, node))
                {
                    return Err(DeclineKind::Ruling);
                }
                continue;
            }
            let (lo, hi) = if nu { (a, b) } else { (b, a) };
            if theta_between(jd, k, def, lo, hi, node)? {
                let t = crate::bands::param_opt(jd, c, def).ok_or(DeclineKind::Ruling)?;
                stations.push((t, Event::Toggle, node));
            }
        }
    }
    stations.sort_by(|x, y| x.0.cmp(&y.0).then(x.1.cmp(&y.1)));
    for pair in stations.windows(2) {
        if pair[0].0 == pair[1].0 && pair[0].2 != pair[1].2 {
            return Err(DeclineKind::CylHoleFeature); // two names for one station
        }
    }
    let mut pieces: Vec<([NodeId; 2], SegKind)> = Vec::new();
    let (mut present, mut graze): (bool, Option<bool>) = (false, None);
    let mut i = 0;
    while i < stations.len() {
        let (t, _, node) = stations[i];
        let mut j = i;
        while j < stations.len() && stations[j].0 == t {
            match stations[j].1 {
                Event::GrazeEnd { toggle } => {
                    if graze.take().is_none() {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    if toggle {
                        present = !present;
                    }
                }
                Event::Toggle => {
                    if graze.is_some() {
                        return Err(DeclineKind::CylHoleFeature); // a toggle inside a run
                    }
                    present = !present;
                }
                Event::GrazeStart { body_above } => {
                    if graze.replace(body_above).is_some() {
                        return Err(DeclineKind::CylHoleFeature); // a run inside a run
                    }
                }
            }
            j += 1;
        }
        if j < stations.len() {
            let next = stations[j].2;
            let kind = match graze {
                Some(body_above) => Some(SegKind::Graze { body_above }),
                None if present => Some(mat),
                None => None,
            };
            if let Some(kind) = kind {
                match pieces.last_mut() {
                    Some((end, k)) if *k == kind && end[1] == node => end[1] = next,
                    _ => pieces.push(([node, next], kind)),
                }
            }
        }
        i = j;
    }
    if present || graze.is_some() {
        return Err(DeclineKind::CylHoleFeature); // an open boundary
    }
    Ok((pieces, grazes))
}

/// **What a holed lateral's ruling actually came out as** — the lock for [`ruling_grazes`].
///
/// ★ This cell's result lives *inside* an operation that still refuses further down, so there is no
/// solid to open and count faces on. The trace's own answer is the thing to hold, the way
/// `arrangement`'s `sides == [-1, 1]` lock already holds the hole-free one.
#[cfg(test)]
pub(crate) mod ruling_probe {
    use super::{SegKind, combinatorics};
    use std::sync::Mutex;

    /// The sentence the ruling label's postcondition panics with — one spelling, so the
    /// commuting oracle's `KNOWN` list (cell ④) names the site by the same constant the
    /// `assert` prints.
    pub(crate) const WRONG_SIDE: &str = "the ruling label took the wrong side of the wall";

    /// One entry per ruling a cycle grazed, in emission order: the kinds of its pieces, beside
    /// the cylinder's origin and the ruling's first and last station — so a reader can pick its
    /// own fixture out of a ledger every test in the binary writes to (a panel's ruling is one
    /// graze, a chain's a transversal then a graze).
    pub(crate) static CARVED: Mutex<Vec<Carved>> = Mutex::new(Vec::new());

    #[derive(Clone, Debug)]
    pub(crate) struct Carved {
        pub origin: [f64; 3],
        pub span: [f64; 2],
        pub kinds: Vec<SegKind>,
    }

    /// **The ruling label's postcondition, one entry per ruling piece** (capability D, D2a).
    ///
    /// ★★★★★ **A second, independent description of the side the derivation picked.**
    /// [`super::ruling_interior_is_even`] derives it from `side · κ · frame_sign`; the check asks the *content*
    /// instead — on the side of the surface the lateral's material lies (inside for a boss,
    /// outside for a bore or a notch: `MergedRuling::orient`) its own solid has material, and on
    /// the other side it does not.
    /// The two share no step, so a disagreement is real. (`ArcLabels`' doc set its own side by
    /// this same content rule, measured.)
    ///
    /// `Some(true)` the two agree · `Some(false)` they contradict · `None` the content does not
    /// distinguish, so the check is **blind** there and only the derivation speaks. A boss whose
    /// own plate surrounds it is exactly such a case, which is why this is recorded rather than
    /// asserted.
    pub(crate) static SIDE_CHECK: Mutex<Vec<Option<bool>>> = Mutex::new(Vec::new());

    /// One entry per ruling piece: whether it got a label at all (`world_rat_sense` may decline).
    pub(crate) static LABELLED: Mutex<Vec<bool>> = Mutex::new(Vec::new());

    /// One entry per graze: the stated `body_above`, beside the **realized** stored normal of the
    /// wall class it is stated against.
    ///
    /// ★★★★★ **The sign has no oracle downstream yet** — the operation this exercises refuses at
    /// the labelling for an unrelated incompleteness, so flipping any factor of `body_above`
    /// changes nothing a test can see (☑ measured: all four controls green). So the claim is
    /// checked against the fixture's own geometry instead: the buried half of the boss lies at
    /// `y < 0` of the wall, so the face is on the stored-normal side exactly when that normal
    /// points at `−y`.
    /// Beside them the cylinder's origin and the run's stations, the fixture's identity.
    pub(crate) static GRAZE_SIDE: Mutex<Vec<GrazeSide>> = Mutex::new(Vec::new());

    #[derive(Clone, Copy, Debug)]
    pub(crate) struct GrazeSide {
        /// The cycle whose run this is — a hole's face lies on one side of the wall, a panel's on
        /// the other, so a reader must say which it is asking about.
        pub kind: combinatorics::CycleKind,
        pub body_above: bool,
        /// The realized stored normal's `y` component of the wall class.
        pub ny: f64,
        pub origin: [f64; 3],
        pub span: [f64; 2],
    }
}

/// Which of the two rulings of `w` this branch node sits on, `0` being a **tangent** wall's
/// single ruling — [`ruling_side_signed`] asked of the point the name denotes. The one spelling:
/// the chart's `ruling_name` reads it too (D5, 1a — it used to carry a twin).
///
/// ★ Cell ⑩ read the `0` off the name's root (`Double`) and cell ⑫ moved it here, to the point:
/// see [`ruling_side_signed`] for why the two stopped agreeing.
pub(crate) fn node_ruling_side(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    w: &[nacre_scalar::Rat; 4],
    n: NodeId,
) -> Option<i8> {
    let (_, cyl, _) = combinatorics::branch_name(n)?;
    let (line, s) = combinatorics::branch_meet(jd, cyl, def, n)?;
    ruling_side_signed(w, def, (&line, &s))
}

/// A branch point's axis parameter: one of the two planes in its name is ⊥ the axis (a cap, a
/// rim), and **that class's** parameter is the point's — a rational, so the ordering along a
/// ruling needs no quadratic comparison at all.
///
/// ★ Cell ⑫: it used to require `wc` in the name — a station the sweep minted. A corner the
/// operand named (`Branch{[cap, t], Double}`, a fillet's tangent corner) lies on `wc`'s ruling
/// too when `wc` contains it, and its parameter is read the same way: from its cap.
fn node_axis_param(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<nacre_scalar::Rat> {
    let (planes, _, _) = combinatorics::branch_name(n)?;
    planes.iter().find_map(|&c| {
        crate::planes::axis_param_of_plane(&combinatorics::class_coeffs_rat(jd, c)?, def)
    })
}

/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
fn trace_on_class(
    input: &combinatorics::TraceInput,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    aliases: &Aliases,
) -> Trace {
    let mut out = Trace::default();
    for (side, which) in [SolidSide::A, SolidSide::B].into_iter().enumerate() {
        trace_one(
            &input.faces[side],
            which,
            wc,
            jd,
            cyls,
            faces,
            plane_ix,
            &input.crossings,
            aliases,
            &mut out,
        );
    }
    out
}

/// Test shims: trace straight from the two solids, deriving [`combinatorics::TraceInput`] on the
/// spot. Production derives it once per boolean (`trace_result_faces`) because it is the same for
/// every class; a test that traces a single class should not have to say so.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn trace_on_class_of(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &combinatorics::EdgeFaces,
    inc_b: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
    crossings: std::collections::HashSet<(usize, usize)>,
) -> Trace {
    let input = combinatorics::trace_input(
        model,
        [(a, inc_a), (b, inc_b)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
        &[],
        crossings,
    );
    trace_on_class(&input, wc, jd, cyls, faces, plane_ix, &Aliases::default())
}

/// One solid only — the second operand slot is filled with the same solid, whose loops are
/// identical, and only `faces[0]` is read.
///
/// ★ It carries the cylinder table since cell ⑩: a disk cap's chord is two branch-pinned nodes of
/// the parity sweep, and ordering them along the line asks the cylinder's statement — the old
/// chord vessel bypassed the sweep, which is why this shim could trace with none.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn trace_one_of(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
    crossings: std::collections::HashSet<(usize, usize)>,
    out: &mut Trace,
) {
    let input = combinatorics::trace_input(
        model,
        [(solid, inc), (solid, inc)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
        cyls,
        crossings,
    );
    trace_one(
        &input.faces[0],
        which,
        wc,
        jd,
        cyls,
        faces,
        plane_ix,
        &input.crossings,
        &Aliases::default(),
        out,
    );
}

/// A coincident-merged arrangement edge on a plane class. When the cut plane is a solid's **cap
/// plane**, every side wall traces the same edge twice — once as the cap's boundary (`Seated`) and
/// once as the wall's top edge (`Transversal`) — because the cap's rim *is* the wall's top edge.
/// These are one geometric edge, and the arrangement (and its DCEL) needs them as one: two
/// coincident `(fp, s)` edges at a vertex would collapse in `angular_order`'s zero bucket.
///
/// The merge keeps the geometry (`wall`, `end`) as one and **preserves every contribution** rather
/// than deciding a single `kind`: which label rule a coincident edge follows (a cap-rim edge flips
/// only the below-bit, seated-style) is verified in the label brick, not guessed here.
#[derive(Clone, Debug)]
pub(crate) struct MergedSeg {
    pub wall: usize,
    /// Read by the next brick (crossings + split); kept here so the merged edge carries its
    /// geometry, not just its contributions.
    pub end: [NodeId; 2],
    /// What pins each endpoint on this edge's line — see [`Seg::end_h`] for the plane case and
    /// [`combinatorics::EndPin`] for why the arc split needed a second arm.
    pub end_h: [combinatorics::EndPin; 2],
    /// Every `(solid, kind)` that produced this one geometric edge. Length 1 when nothing was
    /// coincident.
    pub merged: Vec<(SolidSide, SegKind)>,
    /// The travel sense from `end[0]` to `end[1]`, stated by whoever cut the piece — see
    /// [`combinatorics::Carrier::Plane`]. `None` everywhere except a sub-segment the arc or ruling
    /// split cut, whose ends it already knows the order of.
    pub sense: Option<i8>,
}

/// Merge segments that are the **same geometric edge** — same `wall` and same endpoint-triple set
/// (direction-independent) — into one `MergedSeg`, collecting their contributions. Partial overlap
/// (same `wall`, *different* extent — the E5 case) is left for the per-wall interval overlay in
/// `split_at_crossings` to resolve into non-overlapping sub-segments with unioned contributions.
fn merge_coincident(
    jd: &Judge<'_, WorkingPlane>,
    segs: &[Seg],
    wc: usize,
    aliases: &Aliases,
) -> Vec<MergedSeg> {
    // Key an edge by (wall, sorted endpoint pair) — both folded onto their canonical names first,
    // because two producers can describe one edge with a different wall *and* different endpoint
    // names when planes are concurrent, and it takes both folds for the two keys to coincide.
    let key = |s: &Seg| -> (usize, [NodeId; 2]) {
        let (mut a, mut b) = (aliases.canon_point(s.end[0]), aliases.canon_point(s.end[1]));
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        (aliases.canon_wall(wc, s.wall), [a, b])
    };
    let mut order: Vec<(usize, [NodeId; 2])> = Vec::new();
    let mut groups: HashMap<(usize, [NodeId; 2]), MergedSeg> = HashMap::new();
    for s in segs {
        let k = key(s);
        groups
            .entry(k)
            .or_insert_with(|| {
                order.push(k);
                {
                    // The canonical name of the line, so every producer on it agrees. A handle
                    // stays as recorded while the endpoint's name does: an aliased wall carries
                    // the *same* line, so a handle that pinned an endpoint there still pins it
                    // here. ★ Cell ⑫: when the **name** folds (a tangent corner's `Branch` onto
                    // the `ThreePlane` of its planes with this class) the pin is derived again
                    // from the representative and the line — [`combinatorics::pin_for`], the one
                    // rule — because a pin is a fact about the name beside it, not luggage.
                    let wall = aliases.canon_wall(wc, s.wall);
                    let mut end = [NodeId::three_planes(Canon3::three([0, 1, 2])); 2];
                    let mut end_h = s.end_h;
                    for k in 0..2 {
                        end[k] = aliases.canon_point(s.end[k]);
                        if end[k] != s.end[k] {
                            if let Some(pin) = combinatorics::pin_for(jd, wc, wall, end[k]) {
                                end_h[k] = pin;
                            }
                        }
                    }
                    MergedSeg {
                        wall,
                        end,
                        end_h,
                        merged: Vec::new(),
                        sense: None,
                    }
                }
            })
            .merged
            .push((s.solid, s.kind));
    }
    // Deterministic order: first appearance.
    order
        .into_iter()
        .map(|k| groups.remove(&k).unwrap())
        .collect()
}

/// Split the segments on plane class `wc` into the arrangement's 1-skeleton by a **per-wall 1D
/// interval overlay**: on each wall `W`, the line `wc ∩ W` is cut at the sorted union of every
/// W-segment's endpoints and every different-wall crossing, and one `MergedSeg` is emitted per
/// non-empty sub-interval carrying the **union** of the contributions of the W-segments whose
/// closed extent covers it. A crossing of two segments riding `W`/`o.wall` is the plane triple
/// `{wc, W, o.wall}` — no new point species.
///
/// This resolves E5 (same-wall partial overlap): where two segments overlap, the shared sub-interval
/// carries both their contributions, which the label brick reads per solid. Because sub-intervals
/// run between **distinct** points, no zero-length piece is ever emitted (a crossing coinciding with
/// an endpoint is the same plane class, deduped away). A clean arrangement — no same-wall overlap —
/// has exactly one segment covering each sub-interval, so the output is identical to the naive
/// per-segment split. No re-merge is needed: each `(wall, sub-interval)` is emitted once with its
/// full union, so no two outputs can coincide.
///
/// `Err` (honest reject) on: a degenerate endpoint name (`THREE_PLANES`), or two distinct plane
/// classes coincident on a wall's line — a four-plane concurrency `{wc, W, a, b}` the 3-plane DCEL
/// cannot name (`FOURPLANE`).
/// One wall of a class's arrangement: which plane class it is, the **direction family** its line on
/// that class falls in, and the segments riding it.
///
/// Two walls' lines meet in a point iff their families **differ**, which is why the crossing
/// collector compares `dir` instead of asking a predicate per pair (see [`split_at_crossings`]).
///
/// ★ Everything a wall knows travels with it. The alternative — a `Vec<usize>` of classes beside a
/// `HashMap<class, Vec<usize>>` of segments beside a `Vec<usize>` of families indexed by *position* —
/// is two or three index spaces read in one breath, which is the shape this crate has been bitten by
/// four times (`combinatorics`'s "One `usize`, two meanings"). Same rule as `boolean::Ring` carrying
/// the wall each edge rides rather than deriving it.
struct Wall {
    class: usize,
    dir: usize,
    segs: Vec<usize>,
}

/// **A split point on one wall's line, named by what is common to every point there.**
///
/// The line is `wc ∩ w`, so those two planes are the same for all of them and only the third thing
/// differs: a plane class that cuts the line, or which root of which cylinder. The full
/// [`NodeId`] is derived by [`Split::name`] where it is needed — the same trade
/// [`combinatorics::OnLine::Class`] makes, and for the same measured reason: these vectors are
/// rebuilt once per wall and a rotated fold spends most of itself allocating.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Split {
    /// ★ First, so `Ord` still runs ascending in the class — the tie groups below hand
    /// `group.first()` to the output as the point's name, and that is the "smallest name wins"
    /// rule the arrangement replays on.
    Class(usize),
    Branch {
        cyl: usize,
        root: nacre_topo::QuadRoot,
    },
}

impl Split {
    /// ★ **The branch arm keeps the root the producer wrote, not a re-canonicalized one.** A
    /// segment on `w` in class `wc` has branch ends whose plane pair *is* `{wc, w}`, so the stored
    /// name is already canonical for it and [`Split::name`] can put it back verbatim. Restating it
    /// through `NodeId::branch` would re-run the pair ordering and could flip the root — the trap
    /// cell ⑩ recorded.
    /// ★★★★ **A cylinder pin arrives with a branch name, and that is a producer's invariant, so
    /// this panics rather than naming a refusal.** Every site that writes an `EndPin::Cylinder`
    /// writes the `NodeId::Branch` beside it — the seated walk, the arc split, the chord pass — so
    /// the two disagreeing is a defect in *this* kernel, not a property of the model. The
    /// alternative was `BranchVertexUnnamed`, and its sentence is the mirror image of this case:
    /// *"the point is exactly named; what is missing is that these paths have no other name to
    /// carry it by."* Here the point is **not** exactly named — its two halves contradict. A false
    /// sentence in a reject is worse than a loud stop, which is what `ClassIx::plane` says one
    /// door over: *"a loud panic beats a silently wrong plane."*
    /// ☑ Measured unexercised over the whole suite and the ignored sweep before it was made loud.
    fn of(name: NodeId, pin: combinatorics::EndPin) -> Split {
        match (name, pin) {
            (_, combinatorics::EndPin::Class(r)) => Split::Class(r),
            (NodeId::Branch { cyl, root, .. }, combinatorics::EndPin::Cylinder) => {
                Split::Branch { cyl, root }
            }
            (n, combinatorics::EndPin::Cylinder) => unreachable!(
                "a cylinder pin was written beside a three-plane name: {n:?} — the producer that \
                 made this segment set its two halves apart"
            ),
        }
    }

    /// The point in the form the ordering rule takes, on the line `p ∩ q` it was collected for.
    /// ★ A plane split point costs nothing here — its name *is* its pin, and the rule reads the pin.
    fn on(self, p: usize, q: usize) -> combinatorics::PointOn {
        match self {
            Split::Class(r) => combinatorics::PointOn::Class(r),
            Split::Branch { .. } => combinatorics::PointOn::Branch(self.name(p, q)),
        }
    }

    /// The name this point ships under.
    fn name(self, p: usize, q: usize) -> NodeId {
        match self {
            Split::Class(r) => NodeId::three_planes(Canon3::three([p, q, r])),
            Split::Branch { cyl, root } => {
                let mut planes = [p, q];
                planes.sort_unstable();
                NodeId::Branch { planes, cyl, root }
            }
        }
    }

    /// What pins it — the arrangement ships this beside the name.
    fn pin(self) -> combinatorics::EndPin {
        match self {
            Split::Class(r) => combinatorics::EndPin::Class(r),
            Split::Branch { .. } => combinatorics::EndPin::Cylinder,
        }
    }
}

fn split_at_crossings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    aliases: &mut Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    // ★ **The direction sign of each endpoint, once per segment.** `order_along(wc, wall, i, j)`
    // factors into `orient3d(wc, wall, i, j) × dir_sign(wc, wall, j)`, and the second factor does
    // **not** mention `i` — for a segment's endpoints it is a property of the segment alone. The
    // containment test below sweeps `i` over every wall, so leaving it inside asked the same
    // question `|walls|` times over. (Measured: 1.5M `dir_sign` calls where 113k are distinct.)
    // ★★★★★ **This pass is not plane-only, and the wall that said so is gone.** A split point is
    // a [`Split`] — a plane class **or** which root of which cylinder — and the order comes from
    // [`combinatorics::order_located`], which takes both. The endpoints are not read as class ids
    // up front either: they are located as [`combinatorics::OnLine`]s, so an end a cylinder pinned
    // has a form to be compared by. The refusal that used to stand here
    // (`RejectReason::BranchVertexUnnamed`) went with it — [`Split::of`] and
    // `combinatorics::PointOn::of` **panic** on a cylinder pin beside a three-plane name, because
    // those two halves disagreeing is this kernel's defect and not a shape a model can have.
    // ★★ **Both ends of every segment, located on that segment's own line, once.** This is where
    // `end_c` (the endpoint pins) and `end_ds` (their `dir_sign`s) used to live separately; a
    // [`combinatorics::Located`] carries both, built for the pair `(wc, s.wall)` it will be asked
    // about. Building it here rather than per question is the hoist the containment loop below
    // rests on — measured, that loop runs 1.4–3.2M times per 60-fin fold.
    let end_l: Vec<[combinatorics::OnLine; 2]> = segs
        .iter()
        .map(|s| {
            let mut it = (0..2).map(|k| {
                let on = Split::of(s.end[k], s.end_h[k]).on(wc, s.wall);
                combinatorics::on_line(jd, cyls, wc, s.wall, on)
                    .ok_or_else(|| reject(RejectReason::WitnessNotRational))
            });
            Ok([it.next().unwrap()?, it.next().unwrap()?])
        })
        .collect::<Result<_, BoolError>>()?;

    // The extent question is [`combinatorics::closed_contains`]; what stays here is *which* pair to
    // ask it for.
    // ★★★★★ **The line is `wc ∩ segs[si].wall`, not `wc ∩` the wall being processed.** The
    // collector below asks *other* walls' segments about this wall's crossing, and each of those
    // segments rides its own line — which is also the line `end_l[si]` was located on, and which
    // the rule requires the pair to be.
    // ★★★★★ `q` is a **call-site constant**, not a per-segment lookup: the collector asks one
    // wall's crossing of every segment on **one other wall**, and the cover loop asks about
    // segments on **this** wall. Reading `segs[si].wall` at the question instead cost a touch of a
    // wide `MergedSeg` on a loop that runs millions of times — measured, and this is the fix.
    let closed_contains = |si: usize, at: &combinatorics::Located<'_>, q: usize| -> Option<bool> {
        combinatorics::closed_contains(jd, wc, q, at, &end_l[si])
    };

    // The walls, in first-appearance order for deterministic output — each with its **direction
    // family** and the segments riding it. See [`Wall`].
    //
    // ★ **The family replaces a predicate per wall pair, and that is `|W|²` of them.** Two walls'
    // lines on `wc` meet in a point iff `plane_pair_dir_sign(wc, w, r) != 0`, which is the sign of
    // `det(n_wc, n_w, n_r)` — zero exactly when `n_wc × n_w` lies in `r`, i.e. when `wc ∩ w` is
    // **parallel** to `wc ∩ r`. Parallelism of two lines in one plane is an **equivalence relation**,
    // so the walls partition, and each wall needs only to find its family: ask the representatives
    // until one matches. `|W| × families` questions instead of `|W|²`, and the collector then asks
    // none at all.
    //
    // ★ The relation needs `wc ∩ w` to *be* a line. It is: a wall comes from a `MergedSeg` that
    // rides it, so the meet carries a segment by construction. (Were `n_w ∥ n_wc` the determinant
    // would vanish against every `r` and the partition would collapse to one family.)
    //
    // ★ **Not union-find**, though `Aliases` next door is one. That structure merges *given* pairs;
    // the whole point here is to **not ask** the pairs — transitivity is what lets one question per
    // family stand in for all of them.
    let walls: Vec<Wall> = {
        watch!(S_PART);
        let mut walls: Vec<Wall> = Vec::new();
        // Class → its slot in `walls`, for construction only. Nothing below reads it: a wall's
        // segments and family travel with the wall, so the loops have one index space.
        let mut slot: HashMap<usize, usize> = HashMap::new();
        // One representative class per direction family, in first-appearance order.
        let mut reps: Vec<usize> = Vec::new();
        for (i, s) in segs.iter().enumerate() {
            let at = *slot.entry(s.wall).or_insert_with(|| {
                let dir = reps
                    .iter()
                    .position(|&rep| combinatorics::parallel_carriers(jd, wc, rep, s.wall))
                    .unwrap_or_else(|| {
                        reps.push(s.wall);
                        reps.len() - 1
                    });
                walls.push(Wall {
                    class: s.wall,
                    dir,
                    segs: Vec::new(),
                });
                walls.len() - 1
            });
            walls[at].segs.push(i);
        }
        walls
    };

    #[cfg(test)]
    {
        phase::scale::add(&phase::scale::SEGS, segs.len());
        phase::scale::add(&phase::scale::WALLS, walls.len());
    }

    let mut out = Vec::new();
    for wall in &walls {
        let w = wall.class;
        // Split points on W's line: every W-segment endpoint, plus every real different-wall
        // crossing (a segment on `o.wall` whose closed extent reaches W's line).
        // ★★★★★ **A point, not a plane class.** These used to be `usize` — the third plane naming
        // the point — which cannot say "a cylinder pins this one", and a `Vec<usize>` silently
        // dropped such an endpoint: the two pieces either side of a boss then covered no interval
        // and vanished, leaving the class's boundary open. Measured that way before it was fixed.
        //
        // ★★★★ **`Split`, not a `NodeId` — the name is derived.** On *this* line a point is either
        // a third plane or a root of a cylinder, so `{wc, w}` is common to every one of them and
        // storing it 114k–134k times is waste. Carrying the full name instead cost the collecting
        // loop 36% on a rotated fold (75.8→103.5ms, same trip counts): these vectors are rebuilt
        // per wall, and allocation is where that fixture spends most of its time. The same lesson
        // as [`combinatorics::OnLine::Class`] dropping its `name`, one file over.
        let mut pts: Vec<Split> = Vec::new();
        {
            watch!(COLLECT);
            #[cfg(test)]
            phase::scale::add(&phase::scale::COLLECT_TRIPS, segs.len());
            // ★ A segment's own endpoints are split points on its line — whatever pins them.
            for &i in &wall.segs {
                for k in 0..2 {
                    pts.push(Split::of(segs[i].end[k], segs[i].end_h[k]));
                }
            }
            // ★ **Wall-major, because the question is about walls.** What lands in `pts` is a *wall*
            // — the class naming the crossing — so the loop that used to sweep every segment asked
            // "is `w` parallel to this segment's wall?" once per segment when the answer depends
            // only on the pair. And once one segment on `r` reaches `w`'s line, the rest cannot add
            // anything: `r` is already a split point.
            for other in &walls {
                // ★ Same family ⇒ the two lines are parallel and meet in no point. This also
                // subsumes `other.class == w`: a wall is in its own family, so one test does what
                // an identity check and a predicate call used to.
                if other.dir == wall.dir {
                    continue;
                }
                let r = other.class;
                // ★ One handle per wall pair — `r`'s segments all ask about this same crossing,
                // and they all ride `r`'s line, which is the pair this point is located for.
                let at = combinatorics::locate(jd, cyls, wc, r, combinatorics::PointOn::Class(w))
                    .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
                let mut reaches = false;
                for &i in &other.segs {
                    match closed_contains(i, &at, r) {
                        Some(true) => {
                            reaches = true;
                            break;
                        }
                        Some(false) => {}
                        None => return Err(reject(RejectReason::WitnessNotRational)),
                    }
                }
                if reaches {
                    // ★ A crossing of two plane lines is a three-plane point; only a segment's own
                    // endpoint can be one a cylinder pins.
                    pts.push(Split::Class(r));
                }
            }
        }
        let pts = {
            watch!(SORT);
            // Distinct names (same name = same point), then ordered along the line.
            // ★ Sorting by `NodeId` first keeps the old key's order: `{wc, w, r}`'s sorted triple
            // is monotone in `r`, so a group of tied points still hands its **smallest class** to
            // `group.first()` below, which is the representative and so the emitted name.
            pts.sort_unstable();
            pts.dedup();
            let mut bad = false;
            pts.sort_by(|&x, &y| {
                match combinatorics::order_on(jd, cyls, wc, w, x.on(wc, w), y.on(wc, w)) {
                    Some(-1) => std::cmp::Ordering::Less,
                    Some(1) => std::cmp::Ordering::Greater,
                    Some(_) => std::cmp::Ordering::Equal,
                    // ★ A comparator cannot decline, so it raises a flag the caller reads — the
                    // idiom `split_circles` uses for the same reason.
                    None => {
                        bad = true;
                        std::cmp::Ordering::Equal
                    }
                }
            });
            if bad {
                return Err(reject(RejectReason::WitnessNotRational));
            }
            // ★ Two DISTINCT points ordering equal are one point wearing two names — a four-plane
            // concurrency `{wc, w, ·, ·}`. Record it, and keep one representative as a split point:
            // splitting at both would emit a zero-length piece between them.
            // ★★★★ **A tie that involves a cylinder-pinned point is one the alias table must
            // already know.** The plane fold below is a statement about *plane classes* — it
            // hands `Aliases` a set of them — and a branch name has none to contribute; nor could
            // a fold made here be trusted, since the DCEL keys vertices by name and the two
            // handles would still ship two names. ★ Cell ⑫: the table *does* know such a point
            // when it is an operand's corner a class passes through (the seed,
            // `Aliases::record_on_cylinder`): every handle's name canonicalizes to one
            // representative, the pieces below are emitted under it (`sorted`), and one handle
            // is kept as the split point. A tie the table does not know is refused — the honest
            // floor, as `split_segments_at` and `split_circles` refuse the same shape.
            // ☑ Measured: the refusal was unexercised until cell ⑫'s fixtures reached it — and
            // only with the slab as the first operand, where the corner's own name is not the
            // representative and the tie is between the three-plane name and the class's root.
            let mut reps: Vec<Split> = Vec::with_capacity(pts.len());
            let mut group: Vec<Split> = Vec::new();
            let mut tied_branch = false;
            let flush = |group: &mut Vec<Split>,
                         reps: &mut Vec<Split>,
                         al: &mut Aliases,
                         tied_branch: &mut bool| {
                if group.len() > 1 && group.iter().any(|s| matches!(s, Split::Branch { .. })) {
                    let rep = al.canon_point(group[0].name(wc, w));
                    if group.iter().all(|g| al.canon_point(g.name(wc, w)) == rep) {
                        reps.push(group[0]);
                    } else {
                        *tied_branch = true;
                    }
                    group.clear();
                    return;
                }
                if let Some(&rep) = group.first() {
                    reps.push(rep);
                }
                if group.len() > 1 {
                    let mut set: Vec<usize> = vec![wc, w];
                    set.extend(group.iter().filter_map(|s| match s {
                        Split::Class(r) => Some(*r),
                        Split::Branch { .. } => None,
                    }));
                    set.sort_unstable();
                    set.dedup();
                    al.record(jd, &set);
                }
                group.clear();
            };
            for &r in &pts {
                let same = group.first().is_some_and(|&g| {
                    combinatorics::order_on(jd, cyls, wc, w, g.on(wc, w), r.on(wc, w)) == Some(0)
                });
                if !same {
                    flush(&mut group, &mut reps, aliases, &mut tied_branch);
                }
                group.push(r);
            }
            flush(&mut group, &mut reps, aliases, &mut tied_branch);
            if tied_branch {
                return Err(reject(RejectReason::CoincidentNodes));
            }
            reps
        };
        // ★★ **The fourth plane that *carries* this line, rather than crossing it.**
        //
        // The report above sees a concurrency when two classes **cross** `W`'s line at one point:
        // they order equal, so one point wears two handles. A plane that *contains* the line cannot
        // be seen that way — it is in `W`'s own direction family, which the collector skips
        // (`other.dir == wall.dir`), so it never becomes a handle and never orders against anything.
        //
        // But `Aliases::wall` already knows it: folding two walls onto one line is exactly the
        // statement that a second plane carries it. Every split point on this line therefore lies on
        // every plane of the family, and that is a four-plane concurrency the point fold must hear
        // about — otherwise one class names the point with the canonical wall and another names it
        // with a plane that cuts the line, and the two names reach the seam table unlinked
        // (`SeamAlias`, measured: `[3,11,12]` vs `[3,6,11]` at one point).
        //
        // ★ The report happens on the class whose wall was **folded**, and the class that produces
        // the other name reads the alias afterwards — the fixpoint loop above runs another round
        // whenever this grows the table, which is what carries it across.
        {
            let family = aliases.wall_family(wc, w);
            if family.len() > 1 {
                // ★ Only a point three planes name joins a plane-class fold — a branch name has
                // no third class to contribute to the set.
                for r in pts.iter().filter_map(|s| match s {
                    Split::Class(r) => Some(*r),
                    Split::Branch { .. } => None,
                }) {
                    let mut set: Vec<usize> = vec![wc, r];
                    set.extend(family.iter().copied());
                    set.sort_unstable();
                    set.dedup();
                    aliases.record(jd, &set);
                }
            }
        }
        // A split point carries its own name — folded, because a point rebuilt from a handle would
        // otherwise re-introduce the very alias `merge_coincident` just removed. Canonicalizing
        // here and in the merge means every name **downstream** is already the canonical one, and
        // no later stage has to know the table exists.
        // ★ The **raw** name is what `dedup` and the tie groups above key on, and the canonical one
        // is what ships: fold first and a four-plane concurrency would collapse before it is
        // recorded, and the symptom surfaces much later as `SeamAlias`.
        let sorted = |s: Split, al: &Aliases| al.canon_point(s.name(wc, w));
        // Each sub-interval [p, q] carries the union of the W-segments that cover it. Every segment
        // endpoint is itself a split point, so "covers both ends" means "spans the whole interval".
        watch!(COVER);
        #[cfg(test)]
        {
            phase::scale::add(&phase::scale::PTS, pts.len());
            phase::scale::add(
                &phase::scale::COVER_TRIPS,
                pts.len().saturating_sub(1) * wall.segs.len(),
            );
        }
        for pair in pts.windows(2) {
            let (p, q) = (pair[0], pair[1]);
            let mut merged: Vec<(SolidSide, SegKind)> = Vec::new();
            // The same hoist as the collector: every `w`-segment asks about these two points, and
            // they all ride `w`'s line.
            let loc = |s: Split| {
                combinatorics::locate(jd, cyls, wc, w, s.on(wc, w))
                    .ok_or_else(|| reject(RejectReason::WitnessNotRational))
            };
            let (at_p, at_q) = (loc(p)?, loc(q)?);
            for &i in &wall.segs {
                // ★★★★★ **Short-circuiting, deliberately.** `A && B` never asks the second
                // question when the first says no, and this loop runs 1.6M times per fold —
                // evaluating both eagerly doubled it (48.6→89.9ms, same trip counts). It also
                // keeps the *declines* identical: a second question that cannot be formed is never
                // reached when the first already ruled the segment out.
                let undecided = || reject(RejectReason::WitnessNotRational);
                if !closed_contains(i, &at_p, w).ok_or_else(undecided)? {
                    continue;
                }
                if closed_contains(i, &at_q, w).ok_or_else(undecided)? {
                    merged.extend(segs[i].merged.iter().copied());
                }
            }
            if !merged.is_empty() {
                // The pin is a fact about the name beside it (`pin_for`, cell ⑫): a handle whose
                // point the table represents under another name — a class's root that is an
                // operand's corner — is pinned as that name is on this line.
                let (np, nq) = (sorted(p, aliases), sorted(q, aliases));
                let pin = |s: Split, n: NodeId| {
                    combinatorics::pin_for(jd, wc, w, n).unwrap_or_else(|| s.pin())
                };
                out.push(MergedSeg {
                    wall: w,
                    end: [np, nq],
                    end_h: [pin(p, np), pin(q, nq)],
                    merged,
                    sense: None,
                });
            }
        }
    }
    Ok(out)
}

/// CCW cyclic order of the edges around one arrangement vertex on plane class `w`, **read from no
/// coordinate**. Each edge rides a line `w ∩ fp` and runs in direction `s·(n_w × n_fp)`; it is
/// given as `(fp, s)`. The signed turn between edges i and j is `turn_at`'s atom
/// `s_i·s_j·plane_pair_dir_sign(w, fp_i, fp_j)·orient_sign(w)` (combinatorics.rs).
///
/// The turn sign is transitive only *within an open half-plane* (span < π), where it reproduces
/// the angle order exactly — the textbook Graham-scan fact. So: bucket every edge by its turn
/// against a reference `r = edges[0]` into the two open half-planes, plus the `0`/`π` pole where
/// the turn is 0 (collinear with `r`). The pole splits by same-`fp` opposite direction: same fp,
/// opposite `s` is angle π; otherwise angle 0. Each open bucket is then a real sort. Returned
/// indices are the CCW order starting at `r`.
///
/// **What it relies on, and who supplies it.** Two edges leaving one vertex with a `0` turn are on
/// one line — parallel plus a shared point — and then they have no cyclic order to read. Upstream
/// makes that impossible in two places: `merge_coincident` folds an edge traced twice (its own doc
/// names this bucket as the reason), and `Aliases::union_wall` folds two walls that carry one line
/// so they arrive under a single `fp`. Collinear same-direction overlap (E5) is resolved earlier
/// still, by `split_at_crossings`.
///
/// ★ **That reliance is checked, not assumed.** `Aliases::record` only folds where four or more
/// planes meet at a point, and that every same-line wall pair meets such a point is not
/// established — so a `0` turn that survives to here is an honest `UnorderedEdges` reject rather
/// than an order picked arbitrarily. Measured: it fires nowhere in the suite, and disabling
/// `union_wall` makes the four-plane concurrency model reject, which is what says the reliance is
/// real and the check reaches it.
fn angular_order(
    jd: &Judge<'_, WorkingPlane>,
    w: usize,
    edges: &[combinatorics::EdgeDir],
) -> Result<Vec<usize>, BoolError> {
    // ★ The same atom the winding reads (`combinatorics::turn`), which is what keeps the two from
    // drifting; the `0` it returns is *this* function's to interpret — here it is the `0`/π pole,
    // not the straight angle a ring's turn would reject.
    let cross = |i: usize, j: usize| combinatorics::turn(jd, w, &edges[i], &edges[j]);
    let (mut zero, mut pos, mut pole, mut neg) = (vec![], vec![], vec![], vec![]);
    for i in 0..edges.len() {
        match cross(0, i)? {
            c if c > 0 => pos.push(i),
            c if c < 0 => neg.push(i),
            // Collinear with the reference: angle 0 (same ray) or π (opposite ray). ★ A line and
            // an arc tangent at this vertex with the arc leaving the other way (a fillet's smooth
            // corner, cell ⑩) are π apart too — the structural test cannot see it, the geometric
            // one can ([`combinatorics::tangent_pole`]).
            _ if combinatorics::antiparallel(&edges[i], &edges[0])
                || combinatorics::tangent_pole(jd, w, &edges[i], &edges[0])? =>
            {
                pole.push(i)
            }
            _ => zero.push(i),
        }
    }
    // ★★ **The π pole is a question about the direction's *representation*, not about the turn**
    // — "same wall class, opposite travel sign", and for arcs "same circle, opposite tangent".
    // Both live in `combinatorics::antiparallel` now. ★ Measured: the arc/arc arm never fires here,
    // because segments are numbered before arcs and every crossing has one, so `edges[0]` is always
    // a line; the two arcs at a node sit π apart and land in the two *open* buckets instead.
    //
    // ★★★ **The order below exists only if no two edges share an angle, and that is an
    // assumption about *upstream*, not about this function.** Two edges leaving one vertex whose
    // turn is zero are on one line — parallel plus a shared point — so their cyclic order has no
    // answer, and `next` is built from exactly that order: a guess there is a wrong face, with
    // nothing to say so.
    //
    // Upstream is supposed to make it impossible, in two places, and both are load-bearing:
    // `merge_coincident` folds an edge that was traced twice (its own doc names this bucket as
    // the reason), and `Aliases::union_wall` folds two walls that carry one line. Measured, the
    // second fires 122 times over the coverage suite and this check never does.
    //
    // ★ **But "measured never" is not "cannot".** `Aliases::record` only folds where four or more
    // planes meet at a point (`s.len() < 4` returns early), and that every same-line wall pair
    // meets such a point is *not* established. So the assumption is checked rather than trusted,
    // and a failure is an honest reject — this is the only thing standing where the argument
    // stops.
    //
    // The pairs are checked before sorting rather than inside the comparator: `sort_by` does not
    // compare every pair, so a tie it happens not to evaluate would slip through — and that is
    // precisely the tie that reorders the result. At the degrees this sees (2-3 measured) a
    // bucket holds at most two entries, so the check costs the one comparison the sort would
    // have made anyway.
    for b in [&pos, &neg] {
        for x in 0..b.len() {
            for y in (x + 1)..b.len() {
                if cross(b[x], b[y])? == 0 {
                    return Err(reject(RejectReason::UnorderedEdges));
                }
            }
        }
    }
    // The reference's own bucket: anything else collinear with it that the pole test did not
    // claim is the same ray, which is an overlap `split_at_crossings` should have resolved.
    if zero.len() > 1 {
        return Err(reject(RejectReason::UnorderedEdges));
    }
    // Within an open half-plane, `a` precedes `b` (smaller angle) iff `d_a × d_b > 0`.
    // `turn` can now fail rather than answer — an arc's side is `a + b√c` arithmetic — and a
    // comparator cannot carry that out, so the failure is caught in a flag and raised after the
    // sort (the idiom this file already uses where an exact comparison rides a `sort_by`).
    let by_turn = |v: &mut Vec<usize>| -> Result<(), BoolError> {
        let mut bad = None;
        v.sort_by(|&a, &b| match cross(a, b) {
            Ok(c) if c > 0 => std::cmp::Ordering::Less,
            Ok(c) if c < 0 => std::cmp::Ordering::Greater,
            // Unreachable after the check above, and spelled anyway: a comparator that answers
            // `Greater` both ways is a contract violation, and `sort_by` is then free to produce
            // any permutation.
            Ok(_) => std::cmp::Ordering::Equal,
            Err(e) => {
                bad = Some(e);
                std::cmp::Ordering::Equal
            }
        });
        match bad {
            Some(e) => Err(e),
            None => Ok(()),
        }
    };
    by_turn(&mut pos)?;
    by_turn(&mut neg)?;
    let mut out = zero;
    out.extend(pos);
    out.extend(pole);
    out.extend(neg);
    Ok(out)
}

/// Keep only the edges that carry news: an edge whose [`edge_mask`] is all-false changes no
/// label when crossed, so the cells on its two sides are the same region — of the arrangement
/// *and* of the result, since equal labels mean an equal keep decision. Left in the skeleton it
/// only damages structure: a knife-edge tangency far from everything else is a dangling edge the
/// face walk cannot close (a two-half-edge orbit), which surfaced as `RingOrientation` two layers
/// from its cause.
///
/// This is the 1D twin of `touches`: a single-point tangency is recorded as a vertex but never
/// becomes a segment, and a tangency along a line is now an edge that never enters the DCEL.
/// Where a knife edge coincides with the *other* solid's real boundary, that solid's
/// contributions make the mask non-false and the edge stays — a true boundary cannot be dropped
/// by this filter.
///
/// Mask evaluation moves ahead of the walk here, so a conflict (`EdgeOccupancyConflict`) that
/// used to hide behind a walk-stage reject can now surface first — nearer its cause.
fn drop_newsless(segs: Vec<MergedSeg>) -> Result<Vec<MergedSeg>, BoolError> {
    let mut kept = Vec::with_capacity(segs.len());
    for s in segs {
        if edge_mask(&s.merged)? != [false; 4] {
            kept.push(s);
            continue;
        }
        // The only way to all-false is an even count of same-side grazes per solid (seated flips
        // one bit, a transversal flips two, and an edge has at least one contribution) — anything
        // else dropped here would be a hole in that argument, so it is checked, not assumed.
        debug_assert!(
            s.merged
                .iter()
                .all(|(_, k)| matches!(k, SegKind::Graze { .. })),
            "a newsless edge that is not a graze pair: {:?}",
            s.merged
        );
    }
    Ok(kept)
}

/// One face of the arrangement: the cyclic list of half-edges bounding it, and its winding
/// (`+1` a bounded island, `-1` the unbounded outer contour).
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    pub half_edges: Vec<usize>,
    pub winding: i8,
}

/// Number of connected components of the 1-skeleton (union-find over vertex names joined by each
/// edge — segment, arc, or ruling). The `-1`-winding face count must equal this.
fn component_count(segs: &[MergedSeg], arcs: &[MergedArc], rulings: &[MergedRuling]) -> usize {
    let mut idx: HashMap<NodeId, usize> = HashMap::new();
    let mut id = |t: NodeId, parent: &mut Vec<usize>| -> usize {
        let n = idx.len();
        *idx.entry(t).or_insert_with(|| {
            parent.push(n);
            n
        })
    };
    let mut parent: Vec<usize> = Vec::new();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    // ★ **The arcs join too.** The `-1` count the walk compares against is a count of the
    // 1-skeleton's components, and an arc is an edge of that skeleton like any other — leaving them
    // out would make a crossed circle look like a component of its own and the walk's check would
    // then be off by one wherever a circle was cut.
    let ends = segs
        .iter()
        .map(|s| s.end)
        .chain(arcs.iter().map(|a| a.end))
        .chain(rulings.iter().map(|r| r.end));
    for e in ends {
        let (a, b) = (id(e[0], &mut parent), id(e[1], &mut parent));
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..parent.len())
        .filter(|&i| find(&mut parent, i) == i)
        .count()
}

/// `split_circles`' product: `None` when nothing crossed (the caller keeps its slices), else the
/// split segments, the still-whole circles, the arcs, and each cut circle's [`CutRim`].
type SplitCircles = Option<(
    Vec<MergedSeg>,
    Vec<MergedCircle>,
    Vec<MergedArc>,
    Vec<(usize, CutRim)>,
)>;

/// **Cut every circle a segment crosses into arcs, and the segments with it.**
///
/// ★★★ This is the arc split. What it does *not* do is find the crossings — [`circle_crossings`]
/// already names them, and has since the guard that used to refuse this population was written.
/// The work here is ordering: around the circle (θ, [`nacre_scalar::quad::circular_order_about_seam`])
/// to make arcs, and along each segment (the line parameter, [`cmp_along`]) to make sub-segments.
///
/// ★ A circle nothing crosses is returned **whole**, on the road it has always taken. The parallel
/// road shrinks to the population it is actually about.
fn split_circles(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
    aliases: &Aliases,
) -> Result<SplitCircles, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    // ★ Nothing to cross. This used to be answered further down, after the loop had asked every
    // segment for nothing; the ends below are built for the whole slice at once, so the empty case
    // has to be answered before them or an empty class would start declining. `split_rulings`
    // opens with the same line. ☑ Measured: **61993** times over the suite — this is the common
    // case, not an edge one.
    if circles.is_empty() {
        return Ok(None);
    }
    // ★★ **Both ends of every segment, located on that segment's own line, once** — the hoist
    // `split_at_crossings` makes, for the same reason and with the same obligation: each pair is
    // located for `(wc, sg.wall)`, which is the pair [`combinatorics::closed_contains`] must be
    // called with below.
    let seg_ends: Vec<[combinatorics::OnLine; 2]> = segs
        .iter()
        .map(|s| {
            let mut it = (0..2).map(|k| {
                combinatorics::on_line(
                    jd,
                    cyls,
                    wc,
                    s.wall,
                    Split::of(s.end[k], s.end_h[k]).on(wc, s.wall),
                )
                .ok_or_else(undecided)
            });
            Ok([it.next().unwrap()?, it.next().unwrap()?])
        })
        .collect::<Result<_, BoolError>>()?;
    // Which nodes land on each circle, and which on each segment. Collected together because one
    // crossing is a point of both — splitting only one of them would leave the other's edge running
    // through a vertex it does not have.
    let canon = |x: NodeId| aliases.canon_point(x);
    let mut on_circle: Vec<Vec<NodeId>> = vec![Vec::new(); circles.len()];
    let mut on_seg: Vec<Vec<NodeId>> = vec![Vec::new(); segs.len()];
    // ★★★★★ **A contribution's own ends are cut points, and that is a rule rather than a hope.**
    // Where a lateral face's hole starts is where the material a class sees changes, and the
    // arrangement's own principle is that such a place is a vertex. Feeding both ends of every
    // extent in here is what makes it one — after which no arc straddles a boundary and the
    // distribution below is a containment test with no partial case to arbitrate.
    for (ci, circ) in circles.iter().enumerate() {
        for &(_, _, arc) in &circ.merged {
            if let Some(ends) = arc {
                on_circle[ci].extend(ends);
            }
        }
    }
    // ★★★★★ **Where a segment's ends are, asked once per segment — and a missing one is no longer
    // fatal.** This used to sit inside the circle loop, re-solving the same two points once per
    // circle, and it *stopped* the whole class when either end was a branch point, which has no
    // rational coordinate at any width. Nothing below needs a coordinate to decide anything now;
    // the one thing left that reads it is a **filter**.
    let seg_coords: Vec<[Option<[Rat; 3]>; 2]> = segs
        .iter()
        .map(|s| s.end.map(|t| combinatorics::node_coords_rat(jd, t)))
        .collect();
    for (ci, circ) in circles.iter().enumerate() {
        let (o, m, r) = (circ.def.origin(), circ.def.dir(), circ.def.radius());
        for (si, sg) in segs.iter().enumerate() {
            // A segment whose two ends are one point cuts nothing. Asked by **name**, which is the
            // identity: two spellings of one point are refused upstream, not tolerated here.
            // ☑ Measured unexercised, as its coordinate-comparing predecessor presumably was — the
            // census is bit-identical across the change, which is what says the two agree here.
            if sg.end[0] == sg.end[1] {
                continue;
            }
            // ★★★★★ **The cheap rejection runs where it can be formed, and skipping it asks
            // *more*, not less.** `segment_meets_cylinder` needs both endpoints' coordinates; where
            // one is a branch point the question goes straight to `circle_crossings`, which answers
            // `Miss` exactly when the line misses. Losing the filter costs a solve, never an answer.
            if let [Some(p0), Some(p1)] = &seg_coords[si]
                && (p0 == p1 || !nacre_scalar::segment_meets_cylinder(p0, p1, &o, &m, r))
            {
                continue;
            }
            let Some(xs) = circle_crossings(jd, wc, circ, sg) else {
                return Err(undecided());
            };
            for n in xs {
                // A tangency touches without separating — `segment_meets_cylinder` above already
                // let that shape through, and cutting there would make a zero-length arc.
                //
                // ★★★★★ **That question is answered: the skip is right, and a tangency mints no
                // vertex** (cell ⑥, which opened the tangent wall). The note 60 lines down reads
                // *"a circle cut at exactly one point is **slit**, not divided"* and the two used
                // to disagree; they do not now. Three reasons, none of them new: cell ⑤ measured
                // that the link at such a touch is a *single* circle and that OCCT returns the same
                // body, so minting a vertex would invent structure that is not there; the
                // `UnorderedEdges` this skip's lifting once reached is the **rulings** vocabulary
                // (`side` ±1), which the tangent wall never enters because the gate records no
                // crossing for it; and this arm is not gated on `crossings`, so ⊥ classes behave
                // exactly as they always did once the gate opens. ☑ 15 of 21 measured cells
                // assemble with it in place (the other 6 hold a third plane on the tangent line
                // and are refused earlier, by `CoincidentNodes`).
                if matches!(
                    combinatorics::branch_name(n),
                    Some((_, _, nacre_topo::QuadRoot::Double))
                ) {
                    continue;
                }
                // ★★★★★ **Whether the crossing is on this segment is asked here, in the one
                // vocabulary that can answer it for both kinds of end.** It used to be two plane
                // fences inside `circle_crossings`, which needed a *plane* through each endpoint
                // and so could only be built where every end was three-plane named.
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Branch(n))
                        .ok_or_else(undecided)?;
                #[cfg(test)]
                if let [Some(p0), Some(p1)] = &seg_coords[si] {
                    extent_probe::against_the_fences(
                        jd,
                        sg,
                        [p0, p1],
                        n,
                        &circ.def,
                        &at,
                        &seg_ends[si],
                        wc,
                    );
                }
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, &seg_ends[si])
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                // The crossing's name, as the table knows it (cell ⑫): a class through a
                // corner names the corner again here, and the tables key by name.
                let n = canon(n);
                on_circle[ci].push(n);
                on_seg[si].push(n);
            }
        }
    }

    // ★ Nothing crossed: the caller keeps its own slices and no edge is copied. The common case
    // pays for the question and not for an answer it does not need.
    if on_circle.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    let mut out_circles = Vec::new();
    let mut arcs = Vec::new();
    let mut cut_rims: Vec<(usize, CutRim)> = Vec::new();

    // ★★★★★ **The circle's θ order runs first, because it decides which crossings are cut points
    // at all** (cell ⑫). `circle_crossings` names where the segment's *line* meets the whole
    // circle; but the class's circle is only where its arcs are (`MergedCircle::merged` — a
    // fillet's quarter, a slot's half), and a crossing on the circle's **continuation** lies on
    // no edge of this class. It used to be cut anyway — a vertex on the segment where nothing
    // crosses it (☑ measured: 4 per op on the `tangentline slab 1.6` controls, 16 on
    // `rrect-box`, 0 on the 353 older census rows) — and where that point *is* a vertex of the
    // segment already (a face through the top of a fillet's circle: the user's gusset foot,
    // `y = 0 = −2 + r`) the segment split refused it as two names for one point.
    // **The rule: a crossing is a cut point iff it lies on an arc the class carries.** An arc's
    // own ends are cut points by construction (fed in above); a crossing is on an arc iff the
    // piece of circle before or after it is covered — being on an arc's boundary is being its
    // end, and ends are kept by name.
    //
    // ★ The old order — segments first, "the sharper question before the vaguer one" — was
    // unmeasured and said so; both halves refuse a point wearing two names with the same reject,
    // and identity on both halves is now the alias table's (`canon`), so what reaches either
    // refusal is a coincidence no discovery event recorded — the honest floor.
    struct Ordered {
        nodes: Vec<NodeId>,
        order: Vec<usize>,
        seam_is_node: bool,
    }
    let mut ordered: Vec<Option<Ordered>> = Vec::with_capacity(circles.len());
    for (ci, circ) in circles.iter().enumerate() {
        let raw = std::mem::take(&mut on_circle[ci]);
        if raw.is_empty() {
            ordered.push(None);
            continue;
        }
        let is_end = |n: NodeId| {
            circ.merged
                .iter()
                .any(|&(_, _, arc)| arc.is_some_and(|e| e.contains(&n)))
        };
        // One node per point, asked of the table; an arc's own end keeps its name, because the
        // arc is looked up by it below.
        let mut nodes: Vec<NodeId> = Vec::with_capacity(raw.len());
        for n in raw {
            match nodes.iter().position(|&m| canon(m) == canon(n)) {
                Some(i) => {
                    if is_end(n) && !is_end(nodes[i]) {
                        nodes[i] = n;
                    }
                }
                None => nodes.push(n),
            }
        }
        nodes.sort_unstable();
        let (order, seam_is_node) = match circular_order(jd, circ.cyl, &circ.def, &nodes) {
            Ok(o) => o,
            Err(CircleOrderFail::Undecided) => return Err(undecided()),
            Err(CircleOrderFail::Coincident) => {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        };
        let len = order.len();
        let mut at = vec![0usize; len];
        for (k, &i) in order.iter().enumerate() {
            at[i] = k;
        }
        let place = |x: NodeId| {
            nodes
                .iter()
                .position(|&m| canon(m) == canon(x))
                .map(|i| at[i])
                .expect("an extent's ends were fed into this very split")
        };
        let covers = |k: usize| {
            circ.merged.iter().any(|&(_, _, arc)| match arc {
                None => true,
                Some([p, q]) => {
                    let (a, b) = (place(p), place(q));
                    if a <= b {
                        a <= k && k < b
                    } else {
                        a <= k || k < b
                    }
                }
            })
        };
        let keep: Vec<bool> = (0..len)
            .map(|k| is_end(nodes[order[k]]) || covers(k) || covers((k + len - 1) % len))
            .collect();
        if keep.iter().all(|&b| b) {
            ordered.push(Some(Ordered {
                nodes,
                order,
                seam_is_node,
            }));
            continue;
        }
        // The crossings on no arc leave — the circle and every segment they were put on.
        let dropped: Vec<NodeId> = (0..len)
            .filter(|&k| !keep[k])
            .map(|k| nodes[order[k]])
            .collect();
        for v in on_seg.iter_mut() {
            v.retain(|&n| !dropped.iter().any(|&d| canon(d) == canon(n)));
        }
        let kept: Vec<NodeId> = nodes
            .iter()
            .copied()
            .filter(|n| !dropped.contains(n))
            .collect();
        // The reduced order is the old one filtered: the relative θ order is unchanged, and it
        // still runs from the seam — which is a node only if it was one and stayed.
        let seam_is_node = seam_is_node && keep[0];
        let order: Vec<usize> = order
            .iter()
            .filter(|&&i| keep[at[i]])
            .map(|&i| kept.iter().position(|&m| m == nodes[i]).expect("kept"))
            .collect();
        ordered.push(Some(Ordered {
            nodes: kept,
            order,
            seam_is_node,
        }));
    }
    // ---- segments → sub-segments, in line order ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;
    // ---- circles → arcs, in θ order about the seam ----
    for (ci, circ) in circles.iter().cloned().enumerate() {
        let Some(Ordered {
            nodes,
            order,
            seam_is_node,
        }) = ordered[ci].take()
        else {
            out_circles.push(circ);
            continue;
        };
        // θ position of each node, so an extent can be tested against an arc by rank alone.
        let mut at = vec![0usize; order.len()];
        for (k, &i) in order.iter().enumerate() {
            at[i] = k;
        }
        let place = |x: NodeId| {
            nodes
                .iter()
                .position(|&m| canon(m) == canon(x))
                .map(|i| at[i])
                .expect("an extent's ends were fed into this very split")
        };
        for k in 0..order.len() {
            // ★★★★★ **Each arc takes only the contributions that cover it.** A lateral face with a
            // hole marks its class over an *arc*, and the arc inside the hole is a piece of circle
            // that face does not bound at all — copying the circle's whole contribution list onto
            // it would flip the bits of a face that is not there. Cyclic containment by rank:
            // `[p, q)` runs counter-clockwise from `p`, so arc `k` is inside it exactly when `k`
            // sits between their ranks the same way round.
            let merged = circ
                .merged
                .iter()
                .filter(|&&(_, _, arc)| match arc {
                    None => true,
                    Some([p, q]) => {
                        let (a, b) = (place(p), place(q));
                        if a <= b {
                            a <= k && k < b
                        } else {
                            a <= k || k < b
                        }
                    }
                })
                .map(|&(s, kind, _)| (s, kind))
                .collect::<Vec<_>>();
            // ★ **A piece no contribution covers is not an edge** (cell ⑩, S3). A partial rim — a
            // fillet's quarter, a slot's half — states its arc alone, and the rest of the circle
            // bounds no face on this class; emitted, it stood as a phantom edge that tied with
            // the tangent wall's line at the fillet's corner (`UnorderedEdges`, measured). The
            // straight road drops its newsless pieces (`drop_newsless`); this is the arc's twin.
            // The rim's cut nodes stay in `cut_rims` either way.
            if merged.is_empty() {
                continue;
            }
            arcs.push(MergedArc {
                cyl: circ.cyl,
                def: circ.def.clone(),
                end: [nodes[order[k]], nodes[order[(k + 1) % order.len()]]],
                merged,
            });
        }
        cut_rims.push((
            circ.cyl,
            CutRim {
                nodes: order.iter().map(|&k| nodes[k]).collect(),
                seam_is_node,
            },
        ));
    }

    Ok(Some((out_segs, out_circles, arcs, cut_rims)))
}

/// Why an order around a circle could not be formed — two causes, because two callers turn them
/// into different words (the split into a reject, the tracer into a decline), the same split
/// [`RingFail`]/`decline_of` makes on the segment side.
pub(crate) enum CircleOrderFail {
    /// A θ comparison could not be formed exactly (checked-`Rat` overflow, a missing description).
    Undecided,
    /// Two of the nodes are one point wearing two names.
    Coincident,
}

/// **Order nodes around one circle by θ, the seam first** — the one place that order is decided.
///
/// ★★★ **A crossing on the seam is ordered, not refused — it is the cut point.**
/// `circular_order_about_seam` ranks θ ∈ (0, 2π) and answers `SeamIncident` **by name** for a point
/// at θ = 0, because that point is outside the chart's *total* order. But what arcs need is the
/// **cyclic** order, and a cyclic order tolerates one cut anywhere: the seam point is simply first.
/// (Measured: the very first fixture puts a crossing there — a boss on a plate's edge cuts its own
/// rim exactly on the seam generator, so this is the common case, not an exotic one.)
///
/// ★ At most one node can be seam-incident: two would be the same point, which the adjacency check
/// below refuses as the two-names-for-one-point it is.
///
/// ★★ **Two names for one point.** Two walls crossing the circle at one point are two *different*
/// branch names — a `dedup` cannot see it, and the θ sort puts them adjacent — so an arc of zero
/// length would follow. The segment side asks this question already; asking it here too is what
/// keeps the two sides from disagreeing about what "one point" means.
///
/// Returns the permutation of `nodes` in θ order and whether the first of them is the seam point.
/// `nodes` must already be deduped by name.
pub(crate) fn circular_order(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    nodes: &[NodeId],
) -> Result<(Vec<usize>, bool), CircleOrderFail> {
    use nacre_scalar::quad::{SeamOrder, circular_order_about_seam};
    let meets = nodes
        .iter()
        .map(|&n| combinatorics::branch_meet(jd, cyl, def, n).ok_or(CircleOrderFail::Undecided))
        .collect::<Result<Vec<_>, _>>()?;
    let cmp = |i: usize, j: usize| {
        circular_order_about_seam(
            &def.origin(),
            &def.dir(),
            &def.ref_dir(),
            (&meets[i].0, &meets[i].1),
            (&meets[j].0, &meets[j].1),
        )
    };
    let mut seam: Vec<usize> = Vec::new();
    let mut chart: Vec<usize> = Vec::new();
    for i in 0..meets.len() {
        let on_seam = match cmp(i, i) {
            Some(SeamOrder::SeamIncident { first, .. }) => first,
            Some(SeamOrder::Ordered(_)) => false,
            None => return Err(CircleOrderFail::Undecided),
        };
        if on_seam { seam.push(i) } else { chart.push(i) }
    }
    if seam.len() > 1 {
        return Err(CircleOrderFail::Coincident);
    }
    let mut bad_theta = false;
    chart.sort_by(|&i, &j| match cmp(i, j) {
        Some(SeamOrder::Ordered(o)) => o,
        _ => {
            bad_theta = true;
            core::cmp::Ordering::Equal
        }
    });
    if bad_theta {
        return Err(CircleOrderFail::Undecided);
    }
    let seam_is_node = !seam.is_empty();
    let order: Vec<usize> = seam.into_iter().chain(chart).collect();
    for w in order.windows(2) {
        if matches!(
            cmp(w[0], w[1]),
            Some(SeamOrder::Ordered(core::cmp::Ordering::Equal))
        ) {
            return Err(CircleOrderFail::Coincident);
        }
    }
    Ok((order, seam_is_node))
}

/// **Cut every segment at the points collected on it** — the one emitter the arc split and the
/// ruling split share.
///
/// ★★★★★ **It was written twice, almost to the character.** Both roads collect crossings onto their
/// segments and then have to put them in order, refuse the ones that coincide, and hand each piece
/// the sense the whole segment ran in; the two copies drifted only in the names of their locals. One
/// consequence of joining them is that the arc split's rung is the ruling split's rung too — an end
/// a cylinder pins passes both or neither, rather than one road being taught and the other
/// forgotten, which is the shape of defect this file has already had.
///
/// `on_seg[si]` is drained: a segment nobody crossed is passed through untouched, and the common
/// case copies nothing.
fn split_segments_at(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    on_seg: &mut [Vec<NodeId>],
    aliases: &Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let mut out_segs = Vec::with_capacity(segs.len());
    for (si, sg) in segs.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_seg[si]);
        if nodes.is_empty() {
            out_segs.push(sg);
            continue;
        }
        nodes.sort_unstable();
        // ★ Cell ⑫: identity is the table's — two crossings the table knows as one point are one.
        nodes.dedup_by(|a, b| aliases.canon_point(*a) == aliases.canon_point(*b));
        // ★★★★★ **The ruler is gone; the points are ordered by the rule, not by a parameter.**
        // This used to solve every crossing's `(line, s)`, put both endpoints on that same line with
        // `along`, and sort the parameters — which is why an end a cylinder pinned stopped it: it
        // has no rational coordinate for `along` to take. [`combinatorics::order_located`] answers
        // the same question from the pair `{wc, sg.wall}` and the points' own names, so there is no
        // ruler to lay and nothing to convert.
        //
        // ★★ **The pair is the sorted one, and that is not cosmetic.** The retired ruler was
        // `branch_meet`'s canonical line, whose direction is `n_min × n_max`; asking in call order
        // would flip every comparison on a class where `wc > wall`, taking the emitted pieces and
        // the sense with it. ☑ Differenced against the ruler below.
        //
        // ★★ **The guard that checked all crossings share one line becomes a check on their
        // *names*.** It was there because the sort silently mixed two rulers if it stopped holding;
        // there is no ruler now, but the proposition it stood for — every crossing on this segment
        // is a point of `wc ∩ sg.wall` — still has to hold, and the name says so more cheaply than
        // two `MeetLine` comparisons did. ☑ Measured **unexercised** over the suite, like the ruler
        // check before it: it is here because the sort has no way to notice, not because a fixture
        // is red.
        let pair = {
            let mut p = [wc, sg.wall];
            p.sort_unstable();
            p
        };
        let mut keyed: Vec<(combinatorics::PointOn, NodeId, combinatorics::EndPin)> = Vec::new();
        for &n in &nodes {
            if combinatorics::branch_name(n).is_none() {
                return Err(reject(RejectReason::RingNaming));
            }
            // ★ Cell ⑫: the crossing arrives under the table's name, which may carry another
            // pair — a corner's own — so its pin on *this* line is derived, not assumed
            // (`pin_for`: the cylinder when the pair is the line's, else the plane of its pair
            // that cuts the line). `PointOn::Branch` locates the point by its name either way.
            let Some(pin) = combinatorics::pin_for(jd, pair[0], pair[1], n) else {
                return Err(reject(RejectReason::RingNaming));
            };
            keyed.push((combinatorics::PointOn::Branch(n), n, pin));
        }
        // ★★★★★ **A crossing that *is* an endpoint is one point with one name, so it is deduped
        // rather than refused.** The old sentence here — "one point wearing two names, a three-plane
        // one and a branch one" — is still true and still refused, but only for the case it
        // describes: a crossing at a *three-plane* end really does carry a second name, and the
        // equality check below catches it. Where the end was pinned by a cylinder, its name **is**
        // the crossing's `Branch{planes, cyl, root}` — the same vertex, arrived at twice — and
        // shipping it twice would put a zero-length piece between a point and itself.
        //
        // ★ Cell ⑫: **"the same name" is asked of the alias table.** A tangent corner is a point
        // the cylinder crossing names `Branch{[cap, wc], root}` and the plane road names
        // `ThreePlane([cap, t, wc])`; `Aliases::record_on_cylinder` has joined them, so the
        // crossing folds onto the end here as it would had the names been equal. The end's own
        // name and pin stay on the slot — they are the canonical ones (`merge_coincident`) and
        // the pieces emitted below are keyed by them. What the table does *not* know still
        // meets the equality check below, and that refusal is the honest floor.
        //
        // ★ This is why the ends' indices are remembered instead of read off the tail of the
        // vector: after a dedup the last two entries are no longer the two ends.
        // ☑ Measured over the suite: the dedup fires **400** times and the equality check below
        // still refuses **1** — the two cases are separated by `NodeId`, and both happen.
        let mut ends = (0usize, 0usize);
        for k in 0..2 {
            let same = |e: &(combinatorics::PointOn, NodeId, combinatorics::EndPin)| {
                aliases.canon_point(e.1) == aliases.canon_point(sg.end[k])
            };
            let slot = match keyed.iter().position(same) {
                Some(ix) => {
                    if keyed[ix].1 != sg.end[k] {
                        keyed[ix] = (
                            Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
                            sg.end[k],
                            sg.end_h[k],
                        );
                    }
                    ix
                }
                None => {
                    keyed.push((
                        Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
                        sg.end[k],
                        sg.end_h[k],
                    ));
                    keyed.len() - 1
                }
            };
            if k == 0 { ends.0 = slot } else { ends.1 = slot }
        }
        // ★★ **Nothing between them: every crossing was an end, so the segment is not cut.** It goes
        // back whole rather than being re-emitted as a single piece — which would recompute a sense
        // it already carries, and could decline where the untouched segment does not.
        if keyed.len() == 2 {
            out_segs.push(sg);
            continue;
        }
        // The comparison's **second**-argument form, once per point — the hoist `Located`'s own doc
        // measured: a first argument is built a few times, a second is stored and asked repeatedly.
        let on: Vec<combinatorics::OnLine> = keyed
            .iter()
            .map(|k| combinatorics::on_line(jd, cyls, pair[0], pair[1], k.0).ok_or_else(undecided))
            .collect::<Result<_, BoolError>>()?;
        let mut bad = false;
        let cmp = |i: usize, j: usize| -> Option<core::cmp::Ordering> {
            let a = combinatorics::locate(jd, cyls, pair[0], pair[1], keyed[i].0)?;
            Some(
                match combinatorics::order_located(jd, pair[0], pair[1], &a, &on[j])? {
                    -1 => core::cmp::Ordering::Less,
                    1 => core::cmp::Ordering::Greater,
                    _ => core::cmp::Ordering::Equal,
                },
            )
        };
        let mut order: Vec<usize> = (0..keyed.len()).collect();
        order.sort_by(|&i, &j| match cmp(i, j) {
            Some(o) => o,
            None => {
                bad = true;
                core::cmp::Ordering::Equal
            }
        });
        if bad {
            return Err(undecided());
        }
        #[cfg(test)]
        order_probe::against_the_ruler(jd, cyls, wc, &sg, &keyed, pair);
        // ★ A crossing that lands **on** an endpoint is one point wearing two names — a three-plane
        // one and a branch one — and the DCEL keys vertices by name, so shipping both would make
        // two vertices where there is one. Names the alias table knows as one point were folded
        // above (cell ⑫); what reaches here is a coincidence no discovery event recorded, and
        // refusing it is honest — there is no ground to invent the identity on.
        for w in order.windows(2) {
            if cmp(w[0], w[1]) == Some(core::cmp::Ordering::Equal) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        // ★★★ **The sub-segments' travel sense, taken from the whole segment's own two named
        // ends — exactly, with no frame arithmetic at all.** Every sub-segment runs the way the
        // segment ran, and the pieces below are emitted in **ascending** order; so the only
        // question is whether ascending *is* the segment's direction, and the two ends answer it.
        // (The alternative — turning a meet line's `dir` into an `EdgeDir` sense — needs the
        // canonical-vs-stored turn on two classes and an `f64` dot to decide it. This needs
        // neither: `edge_dir` is still the one place a sense is made.)
        let forward = match cmp(ends.0, ends.1) {
            Some(core::cmp::Ordering::Less) => true,
            Some(core::cmp::Ordering::Greater) => false,
            // Equal is the coincident-endpoint case the windows check above already refused, and
            // `None` is a width decline.
            _ => return Err(undecided()),
        };
        let whole = combinatorics::edge_dir(
            jd,
            cyls,
            wc,
            sg.wall,
            (sg.end[0], sg.end_h[0]),
            (sg.end[1], sg.end_h[1]),
        )?
        .sense()
        .ok_or_else(undecided)?;
        let sense = if forward { whole } else { -whole };
        for w in order.windows(2) {
            let (a, b) = (&keyed[w[0]], &keyed[w[1]]);
            out_segs.push(MergedSeg {
                wall: sg.wall,
                end: [a.1, b.1],
                end_h: [a.2, b.2],
                merged: sg.merged.clone(),
                sense: Some(sense),
            });
        }
    }
    Ok(out_segs)
}

/// **Where a rational point sits along a line** — the retired ruler's half for a three-plane end
/// (a branch end brought its own parameter from [`combinatorics::branch_meet`]).
///
/// ★★ **Production does not lay a ruler any more**, so this survives only inside
/// [`order_probe`], which differences the order rule against what it used to compute. That it needs
/// a *rational* point is the whole reason it had to go: an end a cylinder pins has none.
#[cfg(test)]
fn along(line: &nacre_scalar::quad::MeetLine, p: &[Rat; 3]) -> Option<nacre_scalar::quad::QuadVal> {
    use nacre_scalar::quad::QuadVal;
    let (base, dir) = (line.base(), line.dir());
    let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    let mut rel = [Rat::from_int(0); 3];
    for k in 0..3 {
        rel[k] = p[k].checked_sub(base[k])?;
    }
    let dd = dot(&dir, &dir)?;
    let t = dot(&rel, &dir)?.checked_mul(Rat::new(dd.denom(), dd.numer())?)?;
    Some(QuadVal::from_rat(t))
}

/// The order of two line parameters. Same radicand is a subtraction and a sign; different ones
/// (two cylinders cutting one segment) land on `biquad_sign`, which answers there without
/// declining. `None` is checked-`Rat` overflow — the road's own name, not a shape answer.
fn cmp_along(
    a: &nacre_scalar::quad::QuadVal,
    b: &nacre_scalar::quad::QuadVal,
) -> Option<core::cmp::Ordering> {
    use core::cmp::Ordering;
    use nacre_scalar::Orient;
    let orient = if a.c() == b.c() {
        a.checked_sub(b)?.sign()
    } else {
        nacre_scalar::quad::biquad_sign(
            a.a().checked_sub(b.a())?,
            a.b(),
            Rat::from_int(0).checked_sub(b.b())?,
            Rat::from_int(0),
            a.c(),
            b.c(),
        )?
    };
    Some(match orient {
        Orient::Negative => Ordering::Less,
        Orient::Positive => Ordering::Greater,
        Orient::Zero => Ordering::Equal,
    })
}

/// **Where a segment crosses a circle**, exactly — the points a `VertexDef::Branch` names.
///
/// The segment rides `wc ∩ sg.wall` and the circle is `cylinder ∩ wc`, so a crossing is
/// `plane ∩ plane ∩ cylinder` — the very shape [`nacre_scalar::quad::plane_plane_cylinder`]
/// answers and
/// [`nacre_topo::VertexDef::Branch`] names. Solving along the segment instead would be shorter and
/// would yield a point with **no name**, which the next rung (splitting the circle into arcs)
/// would have to re-derive.
///
/// ★ **Only three of `CylinderMeet`'s variants can arrive here, and that is provable.** A plane
/// cuts a cylinder in a *circle* only when it is ⊥ to the axis, so `wc ∩ sg.wall` is always ⊥ to
/// the axis too:
/// - `CoincidentPlanes`/`ParallelPlanes` — the segment lies on that line, so the line exists;
/// - `OnRuling` — a line on the cylinder and on `wc` would have to lie inside `cylinder ∩ wc`,
///   which is a circle, and a circle contains no line;
/// - `AxisParallelMiss` — a line ⊥ to the axis is not ∥ to it.
///
/// The quadratic cannot degenerate either: its leading coefficient is `|d|² > 0` because `d ⊥ m`.
///
/// `None` is overflow (the caller drops the witness, never the rejection).
///
/// ★★ **It returns the crossings' *names*, and the coordinate is derived from the name** (see
/// [`combinatorics::branch_point`]). The pair is solved in this function's own call order —
/// `(wc, sg.wall)` — and `NodeId::branch` puts it in canonical order, restating the root with it.
/// Solving in ascending order instead would make the correspondence true by construction and leave
/// the canonicalization unexercised, which is where a wrong rule hides; the next mint site (the arc
/// split, walking segments in DCEL order) will not have that luxury either.
/// ★★★★★ **It names where the *line* crosses, and the caller says which of those are on the
/// segment.** Both roots come back; two plane fences through the endpoints used to drop the ones
/// outside, and a fence is a *plane*, which an end a cylinder pinned does not have. Extent is an
/// ordering question, so it belongs with the rule that answers ordering for both kinds of end
/// ([`combinatorics::closed_contains`]) — and the caller already post-filters here anyway, for the
/// tangency.
fn circle_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circ: &MergedCircle,
    sg: &MergedSeg,
) -> Option<Vec<combinatorics::NodeId>> {
    use nacre_scalar::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let w = combinatorics::class_coeffs_rat(jd, wc)?;
    let v = combinatorics::class_coeffs_rat(jd, sg.wall)?;
    let (o, m, r) = (circ.def.origin(), circ.def.dir(), circ.def.radius());
    let roots = match nacre_scalar::quad::plane_plane_cylinder(&w, &v, &o, &m, r)? {
        CylinderMeet::Pair { s, .. } => vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])],
        // ★ `Double`, not `Lo`: the two roots coincide, so a re-sort must leave the name alone.
        CylinderMeet::Tangent { s, .. } => vec![(QuadRoot::Double, QuadVal::from_rat(s))],
        // ★★ **Reachable, and it always was the honest answer.** It used to be unreachable behind
        // the caller's `segment_meets_cylinder` — a segment that meets the solid cylinder has a line
        // that meets its surface. That filter needs both endpoints' coordinates, so it is skipped
        // where one end is a branch point, and the line genuinely can miss. "No crossings" is what
        // a miss means; nothing about the arm changes but the sentence above it.
        // ☑ Measured: **6** times over the suite, where the doc used to say never.
        CylinderMeet::Miss(_) => return Some(Vec::new()),
        other => unreachable!(
            "a circle's class is ⊥ to the axis, so its meet with any wall is ⊥ to the axis and \
             can only miss, touch or cross the cylinder — got {other:?}"
        ),
    };
    Some(
        roots
            .into_iter()
            .map(|(root, _)| combinatorics::NodeId::branch(wc, sg.wall, circ.cyl, root))
            .collect(),
    )
}

/// The retired ruler, kept so the order rule that replaced it can be **differenced against it**.
///
/// ★★★★★ **The two do not agree pointwise, and that is not a defect — it is measured, and it is
/// the reason the difference has to be stated as a *shape* rather than as equality.** The ruler's
/// direction is `n₁ × n₂` built from the classes' **rational coefficients**;
/// [`combinatorics::order_located`] takes its axis sign from the judge's **stored** planes. Those
/// two spellings name the same plane but not the same *side*, so on a class where one of them is
/// stored negated the whole comparison flips. ☑ Measured over the suite: 1222 of 3571 comparisons
/// read the other way.
///
/// So what must hold is not "same answer" but **"same or exactly opposite, per segment"** — a
/// wholesale reversal cancels, because the pieces are emitted in the comparator's own ascending
/// order and `forward` is taken with that same comparator, while a *partial* disagreement would be
/// a genuine reshuffle and would put the sub-segments in the wrong places. That is what is counted.
///
/// ★ The pair handed to the rule is the **sorted** one, matching `NodeId::Branch`'s own convention
/// (`planes` ascending, and `root` defined against that order). ☑ Measured: **36 of 256** segments
/// have `wc > wall`, so the corpus does reach the case where the two orders are opposite calls —
/// and by the paragraph above it does not matter which is taken, because swapping the pair can only
/// turn "same" into "reversed" for a whole segment at once.
#[cfg(test)]
pub(crate) mod order_probe {
    use super::{Judge, MergedSeg, NodeId, WorkingPlane, along, cmp_along, combinatorics};
    use core::sync::atomic::{AtomicUsize, Ordering};

    /// Segments where the ruler could be laid and every point compared.
    pub(crate) static SEGMENTS: AtomicUsize = AtomicUsize::new(0);
    /// …of which the rule read the ruler's order **exactly backwards** — benign, and the common case.
    pub(crate) static REVERSED: AtomicUsize = AtomicUsize::new(0);
    /// …of which the rule agreed with the ruler on some pairs and not others. **This is the defect.**
    pub(crate) static SCRAMBLED: AtomicUsize = AtomicUsize::new(0);
    /// Comparisons both roads called a tie — no direction in them, so they are counted apart.
    pub(crate) static EQ_BOTH: AtomicUsize = AtomicUsize::new(0);
    /// Comparisons where one road called two points **the same place** and the other did not — a
    /// claim about coincidence, not about sequence, so it is kept out of the two counts above.
    pub(crate) static EQUALITY_DISAGREED: AtomicUsize = AtomicUsize::new(0);
    /// Segments reached with `wc < wall`, and with `wc > wall` — the relation that decides whether
    /// the corpus can distinguish the sorted pair from the call-order one at all.
    pub(crate) static WC_BELOW_WALL: AtomicUsize = AtomicUsize::new(0);
    pub(crate) static WC_ABOVE_WALL: AtomicUsize = AtomicUsize::new(0);

    /// The old key: a crossing's parameter on the canonical meet line, an endpoint's `along` on the
    /// same line. `None` where the ruler could not be laid — which is the very shape this cell is
    /// removing, so it is skipped rather than counted.
    fn ruler_keys(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        keyed: &[(combinatorics::PointOn, NodeId, combinatorics::EndPin)],
    ) -> Option<Vec<nacre_scalar::quad::QuadVal>> {
        let mut line = None;
        let mut out = vec![None; keyed.len()];
        for (i, k) in keyed.iter().enumerate() {
            // The endpoints are three-plane named and take the second pass; only a crossing lays
            // the ruler. (Written with `?` at first, which made every call return `None` — the
            // aliveness check below is what said so.)
            let Some((_, cyl, _)) = combinatorics::branch_name(k.1) else {
                continue;
            };
            let def = &cyls.get(cyl)?.def;
            let (l, s) = combinatorics::branch_meet(jd, cyl, def, k.1)?;
            line = Some(l);
            out[i] = Some(s);
        }
        let line = line?;
        for (i, k) in keyed.iter().enumerate() {
            if out[i].is_none() {
                out[i] = Some(along(&line, &combinatorics::node_coords_rat(jd, k.1)?)?);
            }
        }
        out.into_iter().collect()
    }

    pub(crate) fn against_the_ruler(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        wc: usize,
        sg: &MergedSeg,
        keyed: &[(combinatorics::PointOn, NodeId, combinatorics::EndPin)],
        pair: [usize; 2],
    ) {
        if wc < sg.wall {
            WC_BELOW_WALL.fetch_add(1, Ordering::Relaxed);
        } else {
            WC_ABOVE_WALL.fetch_add(1, Ordering::Relaxed);
        }
        let Some(keys) = ruler_keys(jd, cyls, keyed) else {
            return;
        };
        let rule = |p: [usize; 2], i: usize, j: usize| -> Option<core::cmp::Ordering> {
            let a = combinatorics::locate(jd, cyls, p[0], p[1], keyed[i].0)?;
            let b = combinatorics::on_line(jd, cyls, p[0], p[1], keyed[j].0)?;
            Some(
                match combinatorics::order_located(jd, p[0], p[1], &a, &b)? {
                    -1 => core::cmp::Ordering::Less,
                    1 => core::cmp::Ordering::Greater,
                    _ => core::cmp::Ordering::Equal,
                },
            )
        };
        use core::cmp::Ordering::Equal;
        let (mut same, mut opposite) = (0usize, 0usize);
        for i in 0..keyed.len() {
            for j in 0..keyed.len() {
                if i == j {
                    continue;
                }
                let (Some(want), Some(got)) = (cmp_along(&keys[i], &keys[j]), rule(pair, i, j))
                else {
                    return;
                };
                // ★ **Coincidence is a different disagreement from order, and is counted apart.**
                // One side calling two points the same place while the other separates them says
                // nothing about the sequence; it says the `CoincidentNodes` refusal is reachable
                // from one road and not the other, which is its own fact.
                if want == Equal && got == Equal {
                    // ★ Both roads put the two points in the same place. That is agreement about
                    // *coincidence*, and it carries no direction, so counting it as "same order"
                    // makes a wholly reversed segment look scrambled. (It did: one segment in the
                    // suite, and this is what it was.)
                    EQ_BOTH.fetch_add(1, Ordering::Relaxed);
                } else if (want == Equal) != (got == Equal) {
                    EQUALITY_DISAGREED.fetch_add(1, Ordering::Relaxed);
                } else if got == want {
                    same += 1;
                } else {
                    opposite += 1;
                }
            }
        }
        SEGMENTS.fetch_add(1, Ordering::Relaxed);
        match (same, opposite) {
            (0, n) if n > 0 => {
                REVERSED.fetch_add(1, Ordering::Relaxed);
            }
            (_, 0) => {}
            _ => {
                SCRAMBLED.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}

/// The retired plane-fence extent test, kept for one commit so the rule that replaced it can be
/// **differenced against it** rather than argued equal to it.
///
/// The two are the same predicate wherever the fence is well posed: the fence plane is the class
/// that pins one end, so it crosses the line exactly *at* that end, and "the side the far endpoint
/// is on" is the half-line from there. Where it is **not** well posed — the far endpoint sitting on
/// the fence plane, so `want == 0` and every nonzero side is read as outside — the two disagree, and
/// that is what the counters below are for.
///
/// ☑ It expires with the gate: it needs both endpoints' rational coordinates, which is the very
/// demand this cell is removing.
#[cfg(test)]
pub(crate) mod extent_probe {
    use super::{Judge, MergedSeg, WorkingPlane, combinatorics};
    use core::sync::atomic::{AtomicUsize, Ordering};

    pub(crate) static ASKED: AtomicUsize = AtomicUsize::new(0);
    pub(crate) static DISAGREED: AtomicUsize = AtomicUsize::new(0);

    fn by_fences(
        jd: &Judge<'_, WorkingPlane>,
        sg: &MergedSeg,
        ends: [&[nacre_scalar::Rat; 3]; 2],
        meet: &(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal),
    ) -> Option<bool> {
        for k in 0..2 {
            let e = combinatorics::class_coeffs_rat(jd, sg.end_h[k].class()?)?;
            let far = ends[1 - k];
            let mut at_far = e[3];
            for i in 0..3 {
                at_far = at_far.checked_add(e[i].checked_mul(far[i])?)?;
            }
            let want = at_far.numer().signum();
            let got = match nacre_scalar::quad::plane_side(&e, &meet.0, &meet.1) {
                nacre_scalar::Orient::Positive => 1,
                nacre_scalar::Orient::Negative => -1,
                nacre_scalar::Orient::Zero => 0,
            };
            if got != 0 && got != want {
                return Some(false);
            }
        }
        Some(true)
    }

    /// Both verdicts for one crossing, counted. `None` from either side is not a disagreement —
    /// it is a description that could not be formed, and the two roads decline for different causes.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn against_the_fences(
        jd: &Judge<'_, WorkingPlane>,
        sg: &MergedSeg,
        ends: [&[nacre_scalar::Rat; 3]; 2],
        n: combinatorics::NodeId,
        def: &nacre_topo::CylinderDef,
        at: &combinatorics::Located<'_>,
        seg_ends: &[combinatorics::OnLine; 2],
        wc: usize,
    ) {
        let Some((_, cyl, _)) = combinatorics::branch_name(n) else {
            return;
        };
        let Some(meet) = combinatorics::branch_meet(jd, cyl, def, n) else {
            return;
        };
        let (a, b) = (
            by_fences(jd, sg, ends, &meet),
            combinatorics::closed_contains(jd, wc, sg.wall, at, seg_ends),
        );
        if let (Some(a), Some(b)) = (a, b) {
            ASKED.fetch_add(1, Ordering::Relaxed);
            if a != b {
                DISAGREED.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}

/// **Where a segment crosses a cylinder's rulings on a ∥ class** — the straight sibling of
/// [`circle_crossings`], same machinery ([`nacre_scalar::quad::plane_plane_cylinder`], the same
/// fences, the same names), different reachable arms: on a class **parallel** to the axis the
/// meet can be tangent to the lateral — **skipped**, because a touch divides nothing and cutting
/// there would make a zero-length piece — or lie on it, which the ladder refuses by name rather
/// than reasoning it unreachable.
///
/// Returns each crossing with the exact point it was solved at (`(line, s)` — the side and
/// extent tests below read it without re-solving). `None` is overflow.
fn lateral_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    w: &[Rat; 4],
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    sg: &MergedSeg,
) -> Result<Option<Vec<LateralCrossing>>, BoolError> {
    use nacre_scalar::quad::CylinderMeet;
    use nacre_topo::QuadRoot;
    let Some(v) = combinatorics::class_coeffs_rat(jd, sg.wall) else {
        return Ok(None);
    };
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let (line, roots) = match nacre_scalar::quad::plane_plane_cylinder(w, &v, &o, &m, r) {
        Some(CylinderMeet::Pair { line, s }) => {
            (line, vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])])
        }
        // A tangent line touches without separating; cutting there would make a zero-length
        // piece — the same skip the circle side's `Double` arm makes.
        // ★ A **tangent** class touches the cylinder along one ruling (cell ⑩, S3): the
        // segment's line meets it in one point, root `Double` — the crossing the seated face's own
        // tangent piece (`side == 0`) is cut at, named as the scan names it.
        Some(CylinderMeet::Tangent { line, s }) => (
            line,
            vec![(QuadRoot::Double, nacre_scalar::quad::QuadVal::from_rat(s))],
        ),
        Some(CylinderMeet::Miss(_)) | Some(CylinderMeet::AxisParallelMiss(_)) => {
            return Ok(Some(Vec::new()));
        }
        // The segment's line lies on the lateral (a wall meeting the cylinder exactly along a
        // ruling), or the planes degenerate: shapes this ladder does not arrange yet.
        Some(CylinderMeet::OnRuling(_))
        | Some(CylinderMeet::CoincidentPlanes)
        | Some(CylinderMeet::ParallelPlanes) => {
            return Err(reject(RejectReason::RulingBoundNotYet));
        }
        None => return Ok(None),
    };
    // ★ Whether a crossing is **on this segment** is the caller's question now, asked in the one
    // vocabulary that answers it for both kinds of end — the same move the circle side makes.
    Ok(Some(
        roots
            .into_iter()
            .map(|(root, s)| {
                (
                    combinatorics::NodeId::branch(wc, sg.wall, cyl, root),
                    line.clone(),
                    s,
                )
            })
            .collect(),
    ))
}

/// One crossing [`lateral_crossings`] found: its Branch name, and the exact `(line, s)` it was
/// solved at (the side and extent tests read it without re-solving).
type LateralCrossing = (
    combinatorics::NodeId,
    nacre_scalar::quad::MeetLine,
    nacre_scalar::quad::QuadVal,
);

/// `split_rulings`' product — `None` when nothing crossed.
type SplitRulings = Option<(Vec<MergedSeg>, Vec<MergedRuling>)>;

/// **Cut every ruling a segment crosses into pieces, and the segments with it** — the straight
/// twin of [`split_circles`], run after it (the two populations share no class today: a class is
/// ∥ or ⊥ to one cylinder's axis, and the mixed cross-axis case is refused before this).
///
/// The ordering vocabulary is the established one: along a segment, [`cmp_along`] on the
/// canonical meet line; along a ruling, the same comparator on the **axis coordinate**
/// (`m · x`, evaluated exactly from each node's re-solved name).
fn split_rulings(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    segs: &[MergedSeg],
    rulings: &[MergedRuling],
    aliases: &Aliases,
) -> Result<SplitRulings, BoolError> {
    use nacre_scalar::quad::QuadVal;
    let undecided = || reject(RejectReason::WitnessNotRational);
    if rulings.is_empty() {
        return Ok(None);
    }
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
    let dot = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    // A point's axis coordinate, from the exact `(line, s)` it was solved at.
    let coord_at = |def: &nacre_topo::CylinderDef,
                    line: &nacre_scalar::quad::MeetLine,
                    s: &QuadVal|
     -> Option<QuadVal> {
        let m = def.dir();
        let (b, d) = (line.base(), line.dir());
        QuadVal::from_rat(dot(&b, &m)?).checked_add(&s.checked_mul_rat(dot(&d, &m)?)?)
    };
    // Each ruling's end coordinates, re-solved from the names once.
    let ruling_ends: Vec<[QuadVal; 2]> = rulings
        .iter()
        .map(|r| {
            let coord = |n| {
                let (line, s) = combinatorics::branch_meet(jd, r.cyl, &r.def, n)?;
                coord_at(&r.def, &line, &s)
            };
            match (coord(r.end[0]), coord(r.end[1])) {
                (Some(a), Some(b)) => Ok([a, b]),
                _ => Err(undecided()),
            }
        })
        .collect::<Result<_, BoolError>>()?;
    // Distinct cylinders carrying rulings here, each solved once per segment.
    let mut cyl_list: Vec<(usize, &nacre_topo::CylinderDef)> = Vec::new();
    for r in rulings {
        if !cyl_list.iter().any(|&(c, _)| c == r.cyl) {
            cyl_list.push((r.cyl, &r.def));
        }
    }
    let mut on_ruling: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); rulings.len()];
    let mut on_seg: Vec<Vec<combinatorics::NodeId>> = vec![Vec::new(); segs.len()];
    for (si, sg) in segs.iter().enumerate() {
        // ★★★★★ **The skip is gone with the reason it gave.** It said a segment whose end a
        // cylinder pinned "carries no fence planes to test extent with" — true of the fences, which
        // are retired; the extent question is `closed_contains` now and it answers for both kinds of
        // end. The rest of that sentence ("such a segment shares only endpoints with a ruling") was
        // an argument about one narrow population, never a measurement, and it is exactly the sort
        // of claim that goes quietly false when the population widens. So it is asked instead of
        // assumed.
        // ☑ Measured: **200** segments now reach here that the skip dropped, and **392** crossings
        // on them are examined. The census is bit-identical, so the skip's *conclusion* held — none
        // of those crossings splits a ruling. It is a measurement now rather than an argument.
        //
        // A segment whose two ends are one point cuts nothing — asked by **name**, the identity,
        // rather than by a coordinate a branch end does not have. ☑ Measured unexercised.
        if sg.end[0] == sg.end[1] {
            continue;
        }
        // Both ends on this segment's own line, once — the pair `closed_contains` is asked with.
        let mut it = (0..2).map(|k| {
            combinatorics::on_line(
                jd,
                cyls,
                wc,
                sg.wall,
                Split::of(sg.end[k], sg.end_h[k]).on(wc, sg.wall),
            )
            .ok_or_else(undecided)
        });
        let seg_ends = [it.next().unwrap()?, it.next().unwrap()?];
        for &(cyl, def) in &cyl_list {
            let Some(xs) = lateral_crossings(jd, wc, &w, cyl, def, sg)? else {
                return Err(undecided());
            };
            for (n, line, s) in xs {
                // ★ On this segment? The caller's question since the fences went — same rule, same
                // vocabulary as the circle side.
                let at =
                    combinatorics::locate(jd, cyls, wc, sg.wall, combinatorics::PointOn::Branch(n))
                        .ok_or_else(undecided)?;
                if !combinatorics::closed_contains(jd, wc, sg.wall, &at, &seg_ends)
                    .ok_or_else(undecided)?
                {
                    continue;
                }
                // A tangent crossing's side is `0`, from its root — `ruling_side`'s `None` is
                // the sign reading and stays so (cell ⑩).
                let side = match combinatorics::branch_name(n) {
                    Some((_, _, nacre_topo::QuadRoot::Double)) => 0,
                    _ => ruling_side(&w, def, (&line, &s)).ok_or_else(undecided)?,
                };
                let Some(c) = coord_at(def, &line, &s) else {
                    return Err(undecided());
                };
                // Which piece of that ruling (several faces of one solid can put several
                // pieces on one side): closed-extent test on the axis coordinate.
                for (ri, r) in rulings.iter().enumerate() {
                    if r.cyl != cyl || r.side != side {
                        continue;
                    }
                    let [lo, hi] = &ruling_ends[ri];
                    let (a, b) = (
                        cmp_along(&c, lo).ok_or_else(undecided)?,
                        cmp_along(&c, hi).ok_or_else(undecided)?,
                    );
                    use core::cmp::Ordering::*;
                    match (a, b) {
                        (Greater, Less) => {
                            // Strictly inside: the ruling splits here, and the segment with it.
                            on_ruling[ri].push(aliases.canon_point(n));
                            on_seg[si].push(aliases.canon_point(n));
                        }
                        (Equal, _) | (_, Equal) => {
                            // On an end: the T-junction. One point, and it must wear **one**
                            // name — the crossing's `{wc, wall}` pair must be the very pair
                            // that named the ruling's end, or (cell ⑫) the table must know the
                            // two names as one point — else two names share the point.
                            let end_n = if a == Equal { r.end[0] } else { r.end[1] };
                            if aliases.canon_point(end_n) != aliases.canon_point(n) {
                                return Err(reject(RejectReason::CoincidentNodes));
                            }
                            on_seg[si].push(aliases.canon_point(n));
                        }
                        _ => {} // beyond this piece: the surface continues, the face does not
                    }
                }
            }
        }
    }
    if on_ruling.iter().all(|v| v.is_empty()) && on_seg.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    // ---- segments → sub-segments, the arc split's own idiom ----
    let out_segs = split_segments_at(jd, cyls, wc, segs, &mut on_seg, aliases)?;
    // ---- rulings → pieces, in axis order ----
    let mut out_rulings = Vec::with_capacity(rulings.len());
    for (ri, r) in rulings.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_ruling[ri]);
        if nodes.is_empty() {
            out_rulings.push(r);
            continue;
        }
        nodes.sort_unstable();
        nodes.dedup();
        let mut keyed: Vec<(QuadVal, combinatorics::NodeId)> = Vec::new();
        for &n in &nodes {
            let (line, s) =
                combinatorics::branch_meet(jd, r.cyl, &r.def, n).ok_or_else(undecided)?;
            keyed.push((coord_at(&r.def, &line, &s).ok_or_else(undecided)?, n));
        }
        let [lo, hi] = ruling_ends[ri];
        keyed.push((lo, r.end[0]));
        keyed.push((hi, r.end[1]));
        let mut order: Vec<usize> = (0..keyed.len()).collect();
        let mut bad = false;
        order.sort_by(|&i, &j| match cmp_along(&keyed[i].0, &keyed[j].0) {
            Some(o) => o,
            None => {
                bad = true;
                core::cmp::Ordering::Equal
            }
        });
        if bad {
            return Err(undecided());
        }
        for pair in order.windows(2) {
            if cmp_along(&keyed[pair[0]].0, &keyed[pair[1]].0) == Some(core::cmp::Ordering::Equal) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        // Ascending axis coordinate **is** the stored convention (`end[0] → end[1]` ascends
        // `+m`), so each window is already a well-formed piece.
        for pair in order.windows(2) {
            out_rulings.push(MergedRuling {
                cyl: r.cyl,
                def: r.def.clone(),
                side: r.side,
                end: [keyed[pair[0]].1, keyed[pair[1]].1],
                merged: r.merged.clone(),
                #[cfg(test)]
                orient: r.orient,
            });
        }
    }
    Ok(Some((out_segs, out_rulings)))
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

/// **The per-class stages, in one place** — the walk and the nesting.
///
/// ★★★★ **It exists because the sequence was written twice.** `trace_result_faces`' `arrange` and
/// `frame_audit`'s replay ran the same stages against the same edges, and the second is what makes
/// a reject debuggable — so a drift between them is, in `decline_to_reject`'s words, *"the worst
/// possible time to be lying"*. One copy, and the audit is auditing what the boolean ran.
///
/// ★★★ **And it puts every stage behind one `Result`, which is the failure this shape is really
/// for.** While a stopper occupies the caller's socket, it must run before this `Result` is
/// unwrapped, so an out-of-coverage class carries the same name out however far the pipeline
/// got. With the stages held separately that order had to be repeated per stage, and getting it
/// wrong is silent — the suite stays green and only the reject's *name* changes, which is
/// exactly what happened once (`f8eb935`). One `Result` leaves one place to put the `?`, and it
/// is after the socket.
///
/// ★ **`emit_faces` is the last stage, and it cannot fail.** Every other stage returns a
/// `Result`; this one returns its product outright — so a probe of the socket's interception
/// cannot be planted from outside by failing this stage.
fn per_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    kind: BoolKind,
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<Staged, BoolError> {
    let (cells, face_of) = timed!(C_EXTRACT, walk_cells(jd, cyls, wc, edges))?;
    let nesting = timed!(C_NEST, nest_cells(jd, cyls, wc, &cells, edges))?;
    // ★ **The seed is `[false; 4]`, and the argument is why it stays an argument.** The
    // arrangement covers all of space, so its unbounded cells reach infinity, where neither solid
    // is. That is a fact about arranging the *whole* model — restrict the input to a region of
    // space and the unbounded cells become an artifact of the restriction, which is what the
    // parameter records.
    let labels = timed!(
        C_LABEL,
        label_cells(&cells, &face_of, edges, &nesting, [false; 4])
    )?;
    let (faces, disk_labels, arc_labels) = timed!(
        C_EMIT,
        emit_faces(kind, &labels, &cells, edges, jd, wc, &nesting.holes)
    );
    // The cylinder chart's **vertical** lines, read off the pieces this class already made. Both
    // ends carry an axis parameter or the piece is not stated: a ruling whose end has no ⊥ partner
    // has no place on a chart's `z` axis, and inventing one would be worse than leaving it out —
    // so it is dropped here and **counted** by the chart's own census rather than guessed at.
    // ★ Walked by index rather than `enumerate`: the half-edge number below is `2·(ns+na) + 2i`,
    // and indexing keeps `i` used in every build (a discarded `enumerate` index is a lint, and an
    // `allow` for it would be a second thing to keep true).
    let ruling_extents: Vec<(usize, RulingExtent)> = (0..edges.rulings.len())
        .filter_map(|i| {
            let r = &edges.rulings[i];
            let z = [
                node_axis_param(jd, &r.def, r.end[0])?,
                node_axis_param(jd, &r.def, r.end[1])?,
            ];
            // ★★ The chart's vertical answer — see [`RulingExtent::label`]. It is an `Option`
            // rather than a `?` deliberately: dropping the piece here would make the recorded
            // population differ between a test build and a release one, and then the census would
            // be measuring a corpus production never sees.
            let label = {
                let base = 2 * (edges.segs.len() + edges.arcs.len());
                let cell_at = |even: bool| -> Option<Label> {
                    Some(labels[*face_of.get(&(base + 2 * i + usize::from(!even)))?])
                };
                // ★ A **tangent** ruling (`side == 0`, cell ⑩) has no interior side on this class:
                // the plane touches the cylinder, and both cells beside the line lie outside its
                // disk. Its extent carries no label — the chart reads it as a station where the
                // face ends, not as a wall with a chamber behind it.
                let inside = if r.side == 0 {
                    None
                } else {
                    ruling_interior_is_even(jd, wc, r.side)
                };
                let out = inside.and_then(cell_at);
                #[cfg(test)]
                {
                    // ★ The postcondition, checked rather than assumed — see `ruling_probe::SIDE_CHECK`.
                    // ★★ It asks whether a **cell** carries the lateral's solid, not whether the two
                    // cells *differ*: crossing a ruling on this wall crosses the **lateral**, so they
                    // differ always (☑ measured while designing this rung: every piece, without
                    // exception) and that says nothing about the side.
                    // ★ Recorded only where a label was formed, so `None` means one thing — the content
                    // did not distinguish — and never "there was nothing to check".
                    if let (Some(inside), Some(s)) = (inside, r.merged.first().map(|m| m.0)) {
                        let has = |l: Label| {
                            let b = match s {
                                crate::planes::SolidSide::A => 0,
                                crate::planes::SolidSide::B => 2,
                            };
                            l[b] || l[b + 1]
                        };
                        if let (Some(i_l), Some(e_l)) = (cell_at(inside), cell_at(!inside)) {
                            // ★ The lateral's material is inside the cylinder for a boss and
                            // **outside** for a bore or a notch (`orient` — E2-2's re-operated Cut
                            // results were the first such population on a through-axis class), so
                            // the content agrees with the derivation when the cell on the material's
                            // side has the solid and the other does not.
                            let material_inside = r.orient > 0;
                            let verdict = match (has(i_l), has(e_l)) {
                                (a, b) if a != b => Some(a == material_inside),
                                _ => None,
                            };
                            // ★★★★★ **Asserted where the fact is made**, so the coverage is total and
                            // the panic names the offending test — a test reading the ledger afterwards
                            // sees only what ran before it (this ladder measured that hole twice).
                            assert_ne!(
                                verdict,
                                Some(false),
                                "{}: class {wc}, side {}",
                                ruling_probe::WRONG_SIDE,
                                r.side
                            );
                            ruling_probe::SIDE_CHECK
                                .lock()
                                .expect("the probe's lock is never held across a panic")
                                .push(verdict);
                        }
                    }
                    // A tangent piece (`side == 0`) carries no label by design (cell ⑩); the
                    // ledger's proposition is about the rulings a chamber lies behind.
                    if r.side != 0 {
                        ruling_probe::LABELLED
                            .lock()
                            .expect("the probe's lock is never held across a panic")
                            .push(out.is_some());
                    }
                }
                out
            };
            Some((
                r.cyl,
                RulingExtent {
                    wall: wc,
                    side: r.side,
                    end: r.end,
                    z,
                    label,
                    #[cfg(test)]
                    marks: r.merged.clone(),
                },
            ))
        })
        .collect();
    Ok(Staged {
        #[cfg(test)]
        cells,
        #[cfg(test)]
        nesting,
        #[cfg(test)]
        labels,
        faces,
        disk_labels,
        arc_labels,
        ruling_extents,
    })
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
/// ★★★★ **The kind used to be asked with `he >= 2 * segs.len()`, in five places, and that sentence
/// is now false.** Two ranges became three the day a circle could be cut, and every one of those
/// five would read an arc as a circle — `nest_cells`' `ring_of` does it by indexing `segs[he / 2]`,
/// which is not a wrong answer but an **index out of bounds** (measured). One type owns the
/// numbering now, so the question is a `match` the compiler checks rather than an arithmetic
/// comparison five readers each restate.
///
/// ★★ **It borrows or owns** (`Cow`). The split makes new `Vec`s; if this only borrowed, the
/// caller would have to keep them alive and the whole preamble — split, build, walk — would be
/// written once per caller, which is the shape `decline_to_reject`'s doc calls *"two copies … would
/// let them drift"*. Owning lets [`ClassEdges::of`] be the single constructor both callers use.
struct ClassEdges<'a> {
    segs: std::borrow::Cow<'a, [MergedSeg]>,
    arcs: std::borrow::Cow<'a, [MergedArc]>,
    /// Ruling pieces (M6-2 rulings ladder) — a fourth range, empty until a tracer contributes
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

impl<'a> ClassEdges<'a> {
    /// Run the arc split — and, when the class carries them, the chord injection and the ruling
    /// split — and hold the result: **the one place a class's edges are assembled.**
    fn of(
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        wc: usize,
        segs: &'a [MergedSeg],
        circles: &'a [MergedCircle],
        rulings: &'a [MergedRuling],
        aliases: &Aliases,
    ) -> Result<Self, BoolError> {
        use std::borrow::Cow;
        // A class carrying both a circle (⊥ one cylinder's axis) and rulings (∥ another's) is
        // the cross-axis pair population — arrangeable by neither split alone, refused by the
        // ladder's name until its cell.
        if !rulings.is_empty() && !circles.is_empty() {
            return Err(reject(RejectReason::RulingBoundNotYet));
        }
        type AfterCircles<'a> = (
            Cow<'a, [MergedSeg]>,
            Cow<'a, [MergedArc]>,
            Cow<'a, [MergedCircle]>,
            Vec<(usize, CutRim)>,
        );
        let (segs2, arcs, circles2, cut_rims): AfterCircles<'a> =
            match split_circles(jd, cyls, wc, segs, circles, aliases)? {
                Some((s, c, a, r)) => {
                    // ★ The names below rest on this: `split_circles` answers `Some` only when some
                    // circle collected a crossing, and a crossing yields at least one arc (a tangency
                    // is skipped before it is collected). So "was split" and "has arcs" are one fact.
                    debug_assert!(!a.is_empty(), "a split that cut no circle into arcs");
                    (Cow::Owned(s), Cow::Owned(a), Cow::Owned(c), r)
                }
                None => (
                    Cow::Borrowed(segs),
                    Cow::Borrowed(&[][..]),
                    Cow::Borrowed(circles),
                    Vec::new(),
                ),
            };

        let (segs3, rulings2) = match split_rulings(jd, cyls, wc, &segs2, rulings, aliases)? {
            Some((s, r)) => (Cow::Owned(s), Cow::Owned(r)),
            None => (segs2, Cow::Borrowed(rulings)),
        };
        Ok(ClassEdges {
            segs: segs3,
            arcs,
            rulings: rulings2,
            circles: circles2,
            cut_rims,
        })
    }

    /// Where the walk's half-edges end and the circles' pseudo-half-edges begin.
    fn he_count(&self) -> usize {
        2 * (self.segs.len() + self.arcs.len() + self.rulings.len())
    }

    /// **The one classifier.** Everything that used to compare against `2 * segs.len()` asks this.
    fn kind(&self, he: usize) -> HalfEdgeKind {
        let ns = self.segs.len();
        let na = self.arcs.len();
        if he < 2 * ns {
            HalfEdgeKind::Seg(he / 2)
        } else if he < 2 * (ns + na) {
            HalfEdgeKind::Arc((he - 2 * ns) / 2)
        } else if he < self.he_count() {
            HalfEdgeKind::Ruling((he - 2 * (ns + na)) / 2)
        } else {
            HalfEdgeKind::Circle((he - self.he_count()) / 2)
        }
    }

    /// The vertex a half-edge leaves. `he % 2 == 0` takes `end[0]` — one rule, all ranges.
    fn origin(&self, he: usize) -> NodeId {
        match self.kind(he) {
            HalfEdgeKind::Seg(i) => self.segs[i].end[he % 2],
            HalfEdgeKind::Arc(i) => self.arcs[i].end[he % 2],
            HalfEdgeKind::Ruling(i) => self.rulings[i].end[he % 2],
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge has no vertex to leave")
            }
        }
    }

    /// **The one place a `RingEdge` is made from a half-edge.**
    ///
    /// ★★★★ It used to be written twice, and the two spellings **differed**: the walk passed the
    /// sense a split carried onto its sub-segments, and `nest_cells`' ring builder passed
    /// `Carrier::plane(wall)` — that is, `sense: None`. The answers agreed only because the second
    /// never saw a split segment; the day it does, it would drop the one fact the endpoints can no
    /// longer supply. One spelling, so there is nothing to drift.
    fn edge_at(&self, he: usize) -> combinatorics::RingEdge {
        let carrier = match self.kind(he) {
            // ★ `MergedArc::end` runs counter-clockwise about the axis, so the even half-edge
            // travels that way and its twin the other. The convention is stated once, at the
            // split; this is the only place it is read.
            HalfEdgeKind::Arc(i) => {
                let a = &self.arcs[i];
                combinatorics::Carrier::Arc(Box::new(combinatorics::ArcCarrier {
                    cyl: a.cyl,
                    def: a.def.clone(),
                    ccw: he % 2 == 0,
                }))
            }
            HalfEdgeKind::Seg(i) => combinatorics::Carrier::Plane {
                wall: self.segs[i].wall,
                sense: self.segs[i].sense.map(|s| if he % 2 == 0 { s } else { -s }),
            },
            // ★ `MergedRuling::end` runs along `+m`, so the even half-edge travels up and its
            // twin down — the straight reading of the arc convention above.
            HalfEdgeKind::Ruling(i) => {
                let r = &self.rulings[i];
                combinatorics::Carrier::Ruling(Box::new(combinatorics::RulingCarrier {
                    cyl: r.cyl,
                    def: r.def.clone(),
                    side: r.side,
                    up: he % 2 == 0,
                }))
            }
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge is not a ring edge")
            }
        };
        // The endpoints as handles on this edge's line, carried by the producer. Not recovered from
        // the names: a canonical name need not mention `wc` or the wall (see
        // `combinatorics::RingEdge`). ★ An arc's ends are branch points by construction — and a
        // ruling's too — which is what `EndPin::Cylinder` says.
        let (from_h, to_h) = match self.kind(he) {
            HalfEdgeKind::Seg(i) => (self.segs[i].end_h[he % 2], self.segs[i].end_h[1 - he % 2]),
            HalfEdgeKind::Arc(_) | HalfEdgeKind::Ruling(_) | HalfEdgeKind::Circle(_) => (
                combinatorics::EndPin::Cylinder,
                combinatorics::EndPin::Cylinder,
            ),
        };
        combinatorics::RingEdge {
            node: self.origin(he),
            to: self.origin(he ^ 1),
            carrier,
            from_h,
            to_h,
        }
    }
}

/// The DCEL walk itself: half-edges into cells, over the **three** ranges an arrangement can hold.
///
/// - segments `[0, 2·ns)` — an ordinary two-ended edge on a plane's meet with `P`;
/// - arcs `[2·ns, 2·(ns+na))` — a piece of a circle a segment cut, equally two-ended;
/// - uncut circles, appended **after** the walk on pseudo-half-edges `2·(ns+na) + 2i`, because a
///   closed curve with no vertex is not an orbit the walk can express.
fn walk_cells(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<(Vec<Cell>, HashMap<usize, usize>), BoolError> {
    let (segs, arcs, circles) = (&edges.segs, &edges.arcs, &edges.circles);
    let n = segs.len();
    let he_count = edges.he_count();
    let is_arc = |he: usize| matches!(edges.kind(he), HalfEdgeKind::Arc(_));
    let origin = |he: usize| edges.origin(he);
    let edge_of = |he: usize| edges.edge_at(he);

    // Outgoing half-edges per vertex.
    let mut outgoing: HashMap<NodeId, Vec<usize>> = HashMap::new();
    for he in 0..he_count {
        outgoing.entry(origin(he)).or_default().push(he);
    }

    // For each vertex, the cyclic order of its outgoing half-edges (indices into its `outs` list).
    let mut cyclic: HashMap<NodeId, (Vec<usize>, Vec<usize>)> = HashMap::new();
    for (&v, outs) in &outgoing {
        let mut edges = Vec::with_capacity(outs.len());
        {
            watch!(E_ORDER);
            for &he in outs {
                // Direction sign away from `v` toward the far end. ★ This used to call
                // `order_along(target, origin)` — the swapped-argument spelling of what
                // `edge_dir` makes by inverting the result. The two agreed only by the
                // antisymmetry of a difference cancelling an inversion; now there is one.
                // ★ The maker pairs the carrier with the sense, so this site no longer can.
                // ★ The half-edge **leaves** `v`, so its direction is read at `v` — and `dir_at`
                // is where that is checked. Built once per half-edge, which is the hoist the
                // sort below rests on.
                edges.push(combinatorics::dir_at(jd, cyls, wc, &edge_of(he), v)?);
            }
        }
        let ord = timed!(E_ANGULAR, angular_order(jd, wc, &edges))?;
        cyclic.insert(v, (outs.clone(), ord));
    }

    watch!(E_WALK);
    let components = component_count(segs, arcs, &edges.rulings);

    // Try both step directions (predecessor / successor); keep the one whose bounded/outer split
    // is right. Handedness is fixed by `orient_sign(w)` and unknown up front.
    for predecessor in [true, false] {
        let mut next = vec![usize::MAX; he_count];
        for (outs, ord) in cyclic.values() {
            let len = ord.len();
            let step = if predecessor { len - 1 } else { 1 };
            for (local, &he_out) in outs.iter().enumerate() {
                let incoming = he_out ^ 1; // twin arrives at this vertex
                let pos = ord.iter().position(|&x| x == local).unwrap();
                next[incoming] = outs[ord[(pos + step) % len]];
            }
        }

        // Walk next-orbits into faces.
        let mut face_of: HashMap<usize, usize> = HashMap::new();
        let mut cells: Vec<Cell> = Vec::new();
        let mut ok = true;
        for start in 0..he_count {
            if face_of.contains_key(&start) {
                continue;
            }
            let mut cyc = Vec::new();
            let mut he = start;
            loop {
                if face_of.contains_key(&he) {
                    ok = false;
                    break;
                }
                face_of.insert(he, cells.len());
                cyc.push(he);
                he = next[he];
                if he == usize::MAX || cyc.len() > he_count {
                    ok = false;
                    break;
                }
                if he == start {
                    break;
                }
            }
            // ★★ **Three half-edges is the *straight* floor, and it is a fact about lines.** Two
            // straight edges between two points are one edge traced twice; two *arcs* between two
            // points are a lens, and an arc and a chord are a circular segment — both are honest
            // cells. So the floor reads the carriers, and a one-half-edge orbit (a circle a
            // segment only grazed, slit but not divided) is refused in either.
            let floor = if cyc.iter().any(|&h| is_arc(h)) { 2 } else { 3 };
            if !ok || cyc.len() < floor {
                ok = false;
                break;
            }
            // The cell's edges come from the walk, which knows each one's carrier and both
            // handles. The endpoint names cannot supply them: a canonical name need not mention
            // this line's planes at all, which is what retired the derivation check that used to
            // stand here — and an arc's carrier is not a plane class it could name anyway.
            let ring: Vec<combinatorics::RingEdge> = cyc.iter().map(|&h| edge_of(h)).collect();
            let w = combinatorics::loop_winding(jd, cyls, wc, &ring)?;
            cells.push(Cell {
                half_edges: cyc,
                winding: w,
            });
        }
        // ★★★★ **Euler, and it is not decoration — it is the half of "this is a subdivision" the
        // contour count cannot see.** The outer-contour rule counts one number and a walk can
        // satisfy it while attaching the wrong edges: measured, dropping the canonical-to-stored
        // turn in `combinatorics::arc_side` makes this class trace **one 8-half-edge orbit where
        // there are four cells**, and the contour count passes it — a silently wrong subdivision,
        // which is the thing this file exists to refuse. `V − E + F = 2C` sees it (`0 ≠ 2`).
        //
        // `F` is the walk's cell count, which is not the topological face count: each extra
        // component contributes a second contour bounding the same outer region, so `F = f + C − 1`
        // and `V − E + f = 1 + C` becomes this. Measured across the whole suite before it was
        // asserted: **111,750 accepted walks, slack 0 in every one** — so it decides nothing that
        // was passing, and the census is bit-identical.
        let subdivides = |cells: &[Cell]| -> bool {
            let (v, e, f) = (
                outgoing.len() as i64,
                (n + arcs.len() + edges.rulings.len()) as i64,
                cells.len() as i64,
            );
            v - e + f == 2 * components as i64
        };
        if ok
            && cells.iter().filter(|c| c.winding == -1).count() == components
            && subdivides(&cells)
        {
            // ★ **Circle cells, appended after the segment orbits** (M6-2a). Each circle is a
            // 1-edge component the DCEL walk cannot express (orbits need ≥3): a `+1` disk cell
            // and a `−1` contour, with pseudo-half-edges numbered past the segment range —
            // `2n + 2i` (disk side) and its `^1` twin (outside), so the label propagation's
            // twin arithmetic works unmodified (`2n` is even). No crossing machinery is owed: a
            // circle and a segment cannot **cross** on this class — the gate cleared the pair or
            // recorded it — and the one shape that touches, a tangency, is skipped by
            // `split_circles`' `Double` arm rather than cut (a touch divides nothing).
            for (i, _) in circles.iter().enumerate() {
                let he_in = he_count + 2 * i;
                face_of.insert(he_in, cells.len());
                cells.push(Cell {
                    half_edges: vec![he_in],
                    winding: 1,
                });
                face_of.insert(he_in + 1, cells.len());
                cells.push(Cell {
                    half_edges: vec![he_in + 1],
                    winding: -1,
                });
            }
            return Ok((cells, face_of));
        }
    }
    Err(reject(RejectReason::RingOrientation))
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

/// Classify every cell as root / hole / plain `+1` (see [`Nesting`]). A face may carry **any number
/// of holes** (`holes[host]` is a list — e.g. a slab pierced by a U's two prongs), nested **any
/// number of levels** deep (a pocket sealed by a slab, a boss cut after fusing), and the plane may
/// carry **any number of separate bodies** (each contributes its own unbounded contour).
///
/// A root is a `-1` contour inside no `+1` ring at all, so the region it bounds lies outside every
/// cell — that is the one unbounded region, however many cycles bound it. Having none of them is
/// the only impossibility (a closed figure always has an outside), and that is `HOLE_ROOTS`. A hole
/// whose owner is not uniquely determined is `HOLE_DEPTH` (see [`innermost_host`]).
fn nest_cells(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    cells: &[Cell],
    edges: &ClassEdges<'_>,
) -> Result<Nesting, BoolError> {
    let circles = &edges.circles;
    let n = cells.len();
    // Which circle a cell is, `None` for a cell the walk built. ★ The range test is
    // `ClassEdges::kind`'s now: an arc's half-edge is past the segments too, and reading one as a
    // circle would ask the *disk's* exact statement about a shape that is not a disk.
    let circle_of = |c: &Cell| -> Option<usize> {
        match edges.kind(*c.half_edges.first()?) {
            HalfEdgeKind::Circle(i) => Some(i),
            HalfEdgeKind::Seg(_) | HalfEdgeKind::Arc(_) | HalfEdgeKind::Ruling(_) => None,
        }
    };
    // The walk knows each edge's carrier and both handles, so the ring carries them instead of
    // leaving them to be re-derived from the endpoint names. A circle cell has no ring — its
    // containment questions run on the cylinder's exact statement instead.
    let ring_of = |c: &Cell| -> Vec<combinatorics::RingEdge> {
        if circle_of(c).is_some() {
            return Vec::new();
        }
        c.half_edges.iter().map(|&he| edges.edge_at(he)).collect()
    };
    let rings: Vec<Vec<combinatorics::RingEdge>> = cells.iter().map(ring_of).collect();
    let pos: Vec<usize> = (0..n).filter(|&i| cells[i].winding == 1).collect();

    // Union-find over cells (the pattern of `component_count`, but joining cells, not vertices).
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    let circle_ix: Vec<Option<usize>> = cells.iter().map(circle_of).collect();
    let mut holes: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut roots: Vec<usize> = Vec::new();
    for c in (0..n).filter(|&i| cells[i].winding == -1) {
        let mut hosts: Vec<usize> = Vec::new();
        for &r in &pos {
            if cell_in_cell(jd, cyls, wc, &rings, circles, &circle_ix, c, r)? == Some(true) {
                hosts.push(r);
            }
        }
        if hosts.is_empty() {
            roots.push(c);
        } else {
            // ★ **A disk may host** (2026-08-21). This used to `debug_assert` that it could not
            // and then `retain` the disks away — and the assertion's stated reason was already
            // wrong: what kept a polygon out of a disk was not the gate's clearance proof but
            // `circles_meet_no_segment`, which refused the whole shape. Narrowing that guard to
            // *separating* breaks brings the population here, and dropping the host would send the
            // contour to `roots`, where `label_cells` seeds it **void** — the boss over a bore
            // would lose its base face and the shell would open.
            let host = innermost_host(jd, cyls, wc, &rings, circles, &circle_ix, &hosts)?;
            let (rc, rr) = (find(&mut parent, c), find(&mut parent, host));
            parent[rc] = rr;
            holes.entry(host).or_default().push(c);
        }
    }
    if roots.is_empty() {
        return Err(reject(RejectReason::HoleRoots)); // no unbounded contour: not a closed arrangement
    }
    let group_of: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    let mut root_groups: Vec<usize> = roots.iter().map(|&r| group_of[r]).collect();
    root_groups.sort_unstable();
    root_groups.dedup();
    Ok(Nesting {
        group_of,
        root_groups,
        holes,
    })
}

/// Whether a circle's **center** lies inside a polygon ring of the class — the containment
/// witness `nest_cells` uses for a circle contour (the loops cannot **cross** — the gate proved
/// clearance or recorded the crossing, and a tangency touches at most at a point — so one point
/// decides). Exact: the center is `axis ∩ W` (rational), the
/// ring corners are rational meets, and the parity runs in a rational 2D basis of `W`
/// (`point_in_ring_2d_rat` — parity is invariant under the affine projection).
///
/// ★ It takes the cylinder's **statement**, not an arrangement element: the coplanar merge asks
/// the same question of a `Bound::Circle` it is carrying into a merged region, and one spelling
/// serves both.
/// **Where a cylinder's axis meets this class** — the circle's centre on the plane, exact.
///
/// ★ It is rational **whatever way the axis points**: the class has rational coefficients or this
/// says nothing, and the meet is one division. That is why a circle can always name a witness of
/// its own where a *ring* cannot — a ring's corners are branch points and carry radicals.
fn circle_centre_rat(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_scalar::Rat; 3]> {
    use nacre_scalar::Rat;
    let coeffs = combinatorics::class_coeffs_rat(jd, wc)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let (o, m) = (def.origin(), def.dir());
    let dot3 = crate::planes::dot3;
    let nm = dot3(&n, &m)?;
    let no_d = dot3(&n, &o)?.checked_add(coeffs[3])?;
    let t = Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)?;
    let mut p = o;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(p)
}

/// **The circle a ring *is***, when every one of its edges is an arc of one cylinder and they
/// chain the whole way round — `None` otherwise.
///
/// ☑ **Its clauses beyond "the first edge is an arc" are guards, and none of them fired** over the
/// suite (88 acceptances, 0 rejections past that first test). They are kept because each states a
/// proposition proved somewhere else — two circles on one class cannot meet (the gate), a ring is
/// a chain (measured over 153,798 rings), a mixed sense would retrace one arc — and a guard that
/// stops holding is how a producer change is meant to surface here rather than two layers down.
///
/// ★★★★★ **A cut circle's cell is still a disk, and must be asked a disk's question.** The
/// dispatch below reaches [`circle_center_in_ring`] only for a cell whose half-edge is literally a
/// `Circle`; a circle a wall has **split into arcs** takes the polygon road instead, where the
/// probes are plane-triple names and an arc's corners are branch points — so the list comes out
/// **empty** and the road refuses with a name about rays it never cast. The shape is the same
/// circle either way, and this is what says so.
///
/// **Why the chain is the whole circle.** The gate proves every pair of cylinders clear — their
/// surfaces by more than the radius sum, or their faces along an axis
/// ([`crate::planes::cylinder_gate`], else `CylinderPairContact`) — so two
/// circles on one class **cannot meet** — arcs that chain head-to-tail therefore all ride the same
/// circle, and running one way round (`ccw` all equal — a mixed pair would retrace one arc) closes
/// it exactly once. ☑ Measured before this was written: **every** ring the road refused for an
/// empty probe list is either such a chain — and each of those the circle's own centre answers —
/// or a wall panel with no circle at all, which still refuses and now by its own name
/// ([`crate::RejectReason::RingHasNoWitness`]). In cells: the crossing census's `NoClearRay` 44
/// became **8 built** and **20 renamed**, the rest being the other road's (`boolean.rs`).
/// ★ Counted as *cells*, not as raises — a traced boolean runs twice in `debug`, and
/// `reject_census`'s own note forbids reading raise counts as populations.
fn ring_own_circle<'a>(ring: &'a [combinatorics::RingEdge]) -> Option<&'a nacre_topo::CylinderDef> {
    let n = ring.len();
    if n < 2 {
        return None;
    }
    let arc = |e: &'a combinatorics::RingEdge| match &e.carrier {
        combinatorics::Carrier::Arc(ac) => Some(&**ac),
        _ => None,
    };
    let first = arc(&ring[0])?;
    for (i, e) in ring.iter().enumerate() {
        let a = arc(e)?;
        if a.cyl != first.cyl || a.ccw != first.ccw || e.to != ring[(i + 1) % n].node {
            return None;
        }
    }
    Some(&first.def)
}

pub(crate) fn circle_center_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    def: &nacre_topo::CylinderDef,
    ring: &[combinatorics::RingEdge],
) -> Result<bool, BoolError> {
    // ★ Every refusal below is a **value** that could not be formed exactly — a class with no
    // narrow description, a coordinate past `Rat` — which is the road's name, not the gate's.
    // (The gate's own questions are signs and were made total; borrowing its name here pointed
    // at a layer that had already answered.)
    let undecided = || reject(RejectReason::WitnessNotRational);
    let center = circle_centre_rat(jd, wc, def).ok_or_else(undecided)?;
    rational_point_in_ring(jd, cyls, wc, &center, ring)?.ok_or_else(undecided)
}

/// **A rational point strictly inside a straight edge whose two ends are branch corners on one
/// line** (cell ⑩). The chord witness needs the two ends to be one solve's two roots; a cap's
/// section between the rulings of two *coaxial* cylinders — a bore inside a fillet, cut by a wall
/// within both radii — has its ends on two solves, one radical each, and every corner of that cell
/// irrational. Both ends still lie on one rational line (the pair `{wc, wall}`'s meet, the same
/// parametrization from either solve), so a **rational parameter between the two** names a point
/// of the edge's interior exactly: chosen by the realized midpoint, then **verified** against each
/// end in its own radical ([`rational_between`]). `None` for any other edge shape, or when the
/// two solves do not parametrize one line.
fn edge_interior_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &combinatorics::RingEdge,
) -> Option<[nacre_scalar::Rat; 3]> {
    if !matches!(e.carrier, combinatorics::Carrier::Plane { .. }) {
        return None;
    }
    let (pa, ca, _) = combinatorics::branch_name(e.node)?;
    let (pb, cb, _) = combinatorics::branch_name(e.to)?;
    if pa != pb {
        return None;
    }
    let (la, sa) = combinatorics::branch_meet(jd, ca, &cyls.get(ca)?.def, e.node)?;
    let (lb, sb) = combinatorics::branch_meet(jd, cb, &cyls.get(cb)?.def, e.to)?;
    if la.base() != lb.base() || la.dir() != lb.dir() {
        return None;
    }
    let t = rational_between(&sa, &sb)?;
    let (b, d) = (la.base(), la.dir());
    let mut p = b;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(d[k])?)?;
    }
    Some(p)
}

/// A rational strictly between two quadratic values that need not share a radical: the realized
/// midpoint, taken exactly as the `f64` it is (`Rat::try_from_f64`), then **verified** against
/// each end in that end's own radical — a comparison with a rational is always formable. If the
/// midpoint lands outside (ends closer than the realization resolves), a few bisections toward
/// the realized interval's middle are tried; `None` when none is inside.
fn rational_between(
    a: &nacre_scalar::quad::QuadVal,
    b: &nacre_scalar::quad::QuadVal,
) -> Option<nacre_scalar::Rat> {
    use nacre_scalar::{Orient, Rat, quad::QuadVal};
    let (mut lo, mut hi) = (a.to_f64(), b.to_f64());
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    let inside = |t: Rat| -> Option<bool> {
        let q = QuadVal::from_rat(t);
        let da = q.checked_sub(a)?.sign();
        let db = q.checked_sub(b)?.sign();
        // strictly between: on opposite sides of the two ends
        Some(matches!(
            (da, db),
            (Orient::Positive, Orient::Negative) | (Orient::Negative, Orient::Positive)
        ))
    };
    let mut mid = (lo + hi) / 2.0;
    for _ in 0..8 {
        let t = Rat::try_from_f64(mid)?;
        if inside(t)? {
            return Some(t);
        }
        // the realization put it outside: pull toward the interval's middle
        mid = (mid + (lo + hi) / 2.0) / 2.0;
    }
    None
}

/// **The midpoint of a ring edge whose two ends are one solve's two roots** — rational, exactly,
/// and strictly between them.
///
/// ★★★★★ **A chord names its own middle.** `plane_plane_cylinder` builds the pair as
/// `lo = (mid, −half, disc)` and `hi = (mid, +half, disc)` — **one `mid`, shared** — so a segment
/// whose ends are that pair has `base + s.a()·dir` for its midpoint whichever end is asked, with
/// no second solve and no approximation. `disc > 0` for a `Pair`, so it is strictly inside.
///
/// **Conjugacy is a question about names, not values**: the two ends must carry the same canonical
/// plane pair, the same cylinder, and the two roots. That is also what keeps a *piece* of a chord
/// out — an edge cut short by another feature has a different node at one end, and
/// `split_at_crossings` states that "whether a crossing is on this segment is the caller's
/// question", which this answers by refusing to guess.
fn chord_midpoint_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    e: &combinatorics::RingEdge,
) -> Option<[nacre_scalar::Rat; 3]> {
    use nacre_topo::QuadRoot::{Hi, Lo};
    let (pa, ca, ra) = combinatorics::branch_name(e.node)?;
    let (pb, cb, rb) = combinatorics::branch_name(e.to)?;
    if pa != pb || ca != cb || !matches!((ra, rb), (Lo, Hi) | (Hi, Lo)) {
        return None;
    }
    // ★★★★★ **The edge must *be* the segment between its ends, and an arc is not.** Two ends can
    // be one solve's two roots and still be joined by a **curve**: a plane cutting a circle names
    // both crossings, and *either* arc between them carries that same pair of names. The chord's
    // midpoint is then a point strictly inside the circle and **not on this ring at all** — and a
    // point off the ring is not a witness for it, since containment is read from a point *of* `a`
    // and a point in `a`'s interior answers a different question wherever `b` nests inside it.
    // ☑ Measured over the whole lib suite: 120 acceptances, **not one** curved carrier — an
    // all-arc ring is answered by `ring_own_circle` one arm up, and every ring that reaches here
    // offered two straight chords. The guard states the precondition; it does not describe a
    // population.
    if matches!(e.carrier, combinatorics::Carrier::Arc(_)) {
        return None;
    }
    let def = &cyls.get(ca)?.def;
    let (line, s) = combinatorics::branch_meet(jd, ca, def, e.node)?;
    let (b, d) = (line.base(), line.dir());
    let t = s.a();
    let mut p = [b[0], b[1], b[2]];
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(d[k])?)?;
    }
    // ★★★★★ **Asserted where the fact is made, not where it is consumed.** Both claims this
    // function rests on are checkable here and nowhere cheaper: that `MeetLine`'s `base`/`dir`
    // really do parameterize the two planes' meet (so `base + t·dir` is on both), and that the
    // shared `a()` lands **strictly between** the two roots (so it is strictly inside the
    // cylinder, which is what `disc > 0` buys). A producer change that broke either would
    // otherwise surface as a wrong containment answer two layers up.
    debug_assert!(
        [pa[0], pa[1]].iter().all(|&k| {
            combinatorics::class_coeffs_rat(jd, k).is_none_or(|c| {
                let n = [c[0], c[1], c[2]];
                combinatorics::dot3_rat(&n, &p)
                    .and_then(|v| v.checked_add(c[3]))
                    .is_none_or(|v| v == nacre_scalar::Rat::from_int(0))
            })
        }),
        "a chord midpoint is on both of its planes"
    );
    debug_assert_eq!(
        nacre_scalar::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.radius()),
        nacre_scalar::Orient::Negative,
        "a chord midpoint is strictly inside the cylinder"
    );
    Some(p)
}

/// **Is a rational point inside a ring?** — `Ok(None)` when *this point* cannot answer, `Err` when
/// a **value** could not be formed exactly.
///
/// ★★★★★ **The two are different facts and the split is the point of this function.** The body
/// below used to sit inside [`circle_center_in_ring`], where one witness is all there is, so a
/// point that landed *on* the ring could be folded into the same refusal as a class with no
/// rational description. A caller with **several** witnesses must not read them the same way: a
/// point on the ring has the **next witness as its remedy**, a value that cannot be formed does
/// not — the lesson `point_in_component`'s doc states for the road one dimension up.
///
/// ★ **And it dispatches to the strongest road it can.** A ring the rational chart can name takes
/// [`nacre_geom::intersect::point_in_ring_2d_rat`], which is **half-open in y** and so decides
/// even where a ring corner sits on the ray; only a ring with branch corners or arc steps falls
/// to `point_in_mixed_ring`, which abstains there. Anything that hands this a witness gets the
/// better answer for free.
fn rational_point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    p: &[nacre_scalar::Rat; 3],
    ring: &[combinatorics::RingEdge],
) -> Result<Option<bool>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let coeffs = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let center = *p;
    // The class's rational chart — the one copy of that rule ([`combinatorics::Chart2dRat`]);
    // parity is affine-invariant, so the basis need not be orthonormal.
    // ★ **A ring the chart road cannot name takes the mixed road** (M6-2b chaining ladder,
    // wall 3): branch corners have no rational coordinates and arc steps no straight chart
    // image, so the parity walks the ring step by step in ℚ(√c) instead. Rings the old road
    // could always name still take it — the mixed arm activates on exactly the population the
    // old road refused, which is what keeps every green census row bit-identical.
    if combinatorics::ring_is_mixed(ring) {
        return Ok(combinatorics::point_in_mixed_ring(
            jd, cyls, &coeffs, &center, ring,
        ));
    }
    let chart = combinatorics::Chart2dRat::of_normal(&n).ok_or_else(undecided)?;
    let p2 = chart.project(&center).ok_or_else(undecided)?;
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let ring2 = chart.ring(jd, &nodes).ok_or_else(undecided)?;
    Ok(
        match nacre_geom::intersect::point_in_ring_2d_rat(p2, &ring2) {
            nacre_geom::intersect::RingSide::Inside => Some(true),
            nacre_geom::intersect::RingSide::Outside => Some(false),
            // On the boundary this witness says nothing — the caller's next one may.
            nacre_geom::intersect::RingSide::OnBoundary => None,
        },
    )
}

/// Whether a polygon contour lies inside a disk — **all but impossible** (its edges ride wall
/// faces, and a wall face that meets the lateral is either recorded as a crossing, which puts its
/// edges on the ruling road, or *tangent*, which since cell ⑥ passes: a corner of such a face can
/// sit exactly on the tangent line and land a ring node exactly on the circle). Computed honestly
/// from one node's radial side rather than assumed, and the `Zero` arm below is what that leftover
/// reaches.
pub(crate) fn node_in_circle(
    jd: &Judge<'_, WorkingPlane>,
    ring: &[combinatorics::RingEdge],
    def: &nacre_topo::CylinderDef,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    // Any node decides (disjoint loops put every node on one side), so take the first
    // *rational* one — a bitten ring's branch corners have none, but its wall-meet corners do.
    //
    // An all-branch ring (no rational corner at all) still refuses.
    //
    // ★ **A fallback for that was written here and then measured away.** The plan for the road
    // above predicted it would revive contours that arrive here with no rational corner, and the
    // remedy looked free — a ring that *is* a circle can hand over its centre
    // ([`ring_own_circle`] + [`circle_centre_rat`]). ☑ It fired **0** times over the suite: the
    // contours that road revives are asked against *polygons*, which have rational corners. So the
    // note above stands as it was, and the machinery is not here waiting for a population that
    // does not exist.
    let p = ring
        .iter()
        .find_map(|e| combinatorics::node_coords_rat(jd, e.node))
        .ok_or_else(undecided)?;
    match nacre_scalar::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.radius()) {
        nacre_scalar::Orient::Negative => Ok(true),
        nacre_scalar::Orient::Positive => Ok(false),
        // On the surface: a contact the gate admits but this road cannot rank — a *geometric*
        // degeneracy, so it keeps the gate's name ("the gate could not decide exactly") while the
        // width causes above take the road's. ★ It used to say "gate-impossible"; the tangent arm
        // opened in cell ⑥ and a face corner on the tangent line reaches here.
        nacre_scalar::Orient::Zero => Err(reject(RejectReason::CylinderGateUndecided)),
    }
}

/// **Is polygon ring `ra` inside polygon ring `rb`?** — the two-polygon arm of [`cell_in_cell`],
/// written once because the coplanar merge asks the same question of a hole and an outer
/// (cell ⑩: it used to mirror this arm by hand, minus the witness supplies, and named a hole with
/// no three-plane corner `NoClearRay` for a ray never cast). `Ok(None)` is adjacency (a shared
/// node): not nested, not comparable.
pub(crate) fn ring_in_ring_by_witness(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    ra: &[combinatorics::RingEdge],
    rb: &[combinatorics::RingEdge],
) -> Result<Option<bool>, BoolError> {
    if ra.iter().any(|e| rb.iter().any(|f| f.node == e.node)) {
        return Ok(None);
    }
    if combinatorics::ring_is_mixed(rb) {
        let undecided = || reject(RejectReason::WitnessNotRational);
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
        // From each rational witness of `a` until one ray is clear (`None` is the probe
        // on the ring, a tangent ray, or a seam-incident root — a corner on the ray is
        // decided, cell ②); an exhausted ring is the same degeneracy `ring_in_ring`
        // names. ★ The witnesses, in order: the three-plane corners, the **rational
        // branch corners** ([`combinatorics::branch_coords_rat`] — a half-cylinder
        // prism's cap has no other kind, cell ⑩), and the chords' midpoints
        // ([`chord_midpoint_rat`]) — the supplies the named road below has, so the two
        // roads offer the same points.
        let mut asked = false;
        let witnesses = ra.iter().flat_map(|e| {
            combinatorics::node_coords_rat(jd, e.node)
                .or_else(|| combinatorics::branch_coords_rat(jd, cyls, e.node))
                .into_iter()
                .chain(chord_midpoint_rat(jd, cyls, e))
                .chain(edge_interior_rat(jd, cyls, e))
        });
        for p in witnesses {
            asked = true;
            if let Some(hit) = combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, rb) {
                return Ok(Some(hit));
            }
        }
        // ★ Same split as the named road below: a list that was **empty** is not a list
        // that ran out. ☑ Measured (cell ⑩): the half-cylinder pair reached this arm with
        // an empty list before the branch corners joined it.
        return Err(reject(if asked {
            RejectReason::NoClearRay
        } else {
            RejectReason::RingHasNoWitness
        }));
    }
    // ★ `ring_in_ring` casts from each of `a`'s nodes until one gives a clear ray; an
    // exhausted ring is the genuine degeneracy it rejects for.
    //
    // ★★★★★ **This used to end «— which is also why dropping a branch node from the probe
    // list here is honest», and that was the defect wearing the word.** Dropping them left
    // *nothing*: measured, every refusal this road raised came from a list that was empty
    // before the first cast. What is honest is to say so, which the two arms above now do.
    let probes = combinatorics::three_plane_probes(ra.iter().map(|e| e.node));
    // ★★★★★ **A ring that *is* a circle is asked the circle's question.** The probe list
    // above is plane-triple **names**, and a circle a wall has split into arcs has none —
    // every corner is a branch point — so it comes out empty and `ring_in_ring` refuses
    // with a name about rays it never cast. The cell is a disk either way, and the arm one
    // match-arm up already answers disks exactly ([`circle_center_in_ring`], whose witness
    // is the centre and so is rational whatever way the axis points).
    //
    // ★ **Only where the names run out.** Where they do not, today's road and today's
    // order are untouched — a coordinate is not a name, and a rotated class has names but
    // no rational coefficients (`WorkingPlane::world_rat` is `None` there), so keying this
    // on "no rational witness" instead would divert the rotation sweep's own population.
    // Asking the circle question *always* is the tidier end state and should be measured
    // as an agreement first; this is the strict extension.
    if probes.is_empty() {
        if let Some(def) = ring_own_circle(ra) {
            return Ok(Some(
                circle_center_in_ring(jd, cyls, wc, def, rb)? && !node_in_circle(jd, rb, def)?,
            ));
        }
        // ★★★★★ **And when the corners cannot name a witness, an edge can.** A ring edge
        // whose two ends are one solve's two roots has a **rational midpoint**
        // ([`chord_midpoint_rat`]) — the pair is built from a shared `mid`, so it costs one
        // `branch_meet` and no approximation. That is the wall panel's case: four branch
        // corners, no circle to take a centre from, and two perpendicular traces that are
        // each a whole chord.
        //
        // ★ **`Ok(None)` is this witness's abstention, not the ring's** — the next edge's
        // midpoint may still answer, which is why [`rational_point_in_ring`] hands the two
        // apart. ☑ Measured before this was written: every ring here offers **two**
        // midpoints and the two always agree, and none of them abstains.
        //
        // ★ It sits after [`ring_own_circle`] for the reader's sake only: the two supplies
        // are **disjoint by construction** — that one needs every edge to be an arc, this
        // one needs an edge that is not.
        for e in ra {
            let Some(mid) = chord_midpoint_rat(jd, cyls, e) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &mid, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And a **rational branch corner** is a witness too (cell ⑩): a slot's stadium
        // has four, all tangent points, and neither a circle nor a whole chord — the
        // supply the mixed road above offers ([`combinatorics::branch_coords_rat`]).
        for e in ra {
            let Some(p) = combinatorics::branch_coords_rat(jd, cyls, e.node) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &p, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And a rational point **inside an edge** whose ends are two solves' roots
        // ([`edge_interior_rat`]) — the cell between a bore's and a fillet's rulings.
        for e in ra {
            let Some(p) = edge_interior_rat(jd, cyls, e) else {
                continue;
            };
            if let Some(hit) = rational_point_in_ring(jd, cyls, wc, &p, rb)? {
                return Ok(Some(hit));
            }
        }
        // ★ And when neither a circle nor a chord names one, say **that** —
        // `ring_in_ring` below would report an exhausted probe list, which is a different
        // fact and one that never happened here.
        return Err(reject(RejectReason::RingHasNoWitness));
    }
    combinatorics::ring_in_ring(jd, wc, &probes, rb).map(Some)
}

/// **Is cell `a`'s loop inside cell `b`'s?** — the one dispatch, four arms by carrier kind.
///
/// `Ok(None)` is *not comparable*, and it covers two shapes that both mean "these are neighbours,
/// not nested": two circles (see below), and two polygons that **share a node** — a shared node is
/// a split point, so the loops touch rather than one wrapping the other (that also excludes a
/// contour's own `+1` partner, which carries the same ring).
///
/// ★★ **Written once because it is asked from two directions.** [`nest_cells`] asks it of
/// (contour, `+1` cell) to find hosts, and [`innermost_host`] asks it of (host, host) to order
/// them. The four arms are the same question either way; two spellings of it would be free to
/// drift, which is this repository's most-repeated defect.
///
/// The arms:
/// - **two circles** — [`disk_in_disk`] on the radii and the centre distance (cell ⑩ S1: two disks
///   of one class *do* nest — a pin fused onto a boss — the gate keeps only different classes apart);
/// - **circle in polygon** — the centre (`axis ∩ wc`, rational) inside the ring, **and the ring not
///   inside the circle**. ★★ That second clause is not belt-and-braces: with disjoint loops the
///   centre test alone says "inside" for *both* nestings when the polygon happens to straddle the
///   centre — a boss standing over a bore, footprint `[7,9]×[9,11]` around the axis `(8,10)`, put
///   the **circle** inside the **square**. `circle_center_in_ring`'s own doc says one point decides
///   *"the loops cannot cross"*, and cross they do not; what it does not settle is **which way
///   round**. The other direction does, so the pair is the predicate;
/// - **polygon in circle** — one node's radial side ([`node_in_circle`]), decisive on its own:
///   disjoint loops put every node on one side;
/// - **polygon in polygon** — a ray from each of `a`'s nodes until one is clear
///   ([`combinatorics::ring_in_ring`]); when `b` is a **mixed** ring (branch corners, arc steps —
///   a plate's section bitten by a boss, once the scan names crossings on rulings), the ray is
///   [`combinatorics::point_in_mixed_ring`]'s from each of `a`'s rational nodes, the predicate the circle arm
///   already asks of a centre. ★ The old road *always* exhausts on a mixed `b` (every ray answers
///   `Unnameable` at the first branch corner), so this arm activates on exactly the population it
///   refused — every other row stays bit-identical.
#[allow(clippy::too_many_arguments)]
fn cell_in_cell(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circles: &[MergedCircle],
    circle_ix: &[Option<usize>],
    a: usize,
    b: usize,
) -> Result<Option<bool>, BoolError> {
    match (circle_ix[a], circle_ix[b]) {
        // ★ **Two disks** (cell ⑩). A disk lies inside another iff its rim does, and the gate
        // keeps the rims of two classes apart, so the question is one rational inequality on
        // the plane: `r_a < r_b` and `(r_b − r_a)² > |c_b − c_a|²` (a boss's cap under a pin's, a
        // tube's two rims around a pin's trace — concentric or not). ★ This arm used to answer
        // `None`, which left two disk cells unnested and *silently* kept both operands' caps: a
        // pin stacked on a boss fused into **two** untouched bodies (measured the day the pair
        // rule stopped refusing nested parallel cylinders). The old refusal was masking a gap.
        (Some(ca), Some(cb)) => disk_in_disk(jd, wc, &circles[ca].def, &circles[cb].def).map(Some),
        (Some(ci), None) => Ok(Some(
            circle_center_in_ring(jd, cyls, wc, &circles[ci].def, &rings[b])?
                && !node_in_circle(jd, &rings[b], &circles[ci].def)?,
        )),
        (None, Some(ri)) => node_in_circle(jd, &rings[a], &circles[ri].def).map(Some),
        (None, None) => ring_in_ring_by_witness(jd, cyls, wc, &rings[a], &rings[b]),
    }
}

/// **Whether the disk cylinder `a` cuts on class `wc` lies inside the disk `b` cuts there** — one
/// rational inequality on the plane: `r_a < r_b` and `(r_b − r_a)² > |c_b − c_a|²`. A disk lies
/// inside another iff its rim does, and the rims of two classes never meet (the gate proves the
/// faces clear, or the solid is valid), so the two centres and radii decide — concentric or not.
///
/// ★ **Written once because it is asked from two directions** (cell ⑩): the arrangement's nesting
/// ([`cell_in_cell`], a pin's trace under a boss's cap) and the coplanar merge (a region whose
/// outer bound is a circle asking which circle holes it owns).
pub(crate) fn disk_in_disk(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    a: &nacre_topo::CylinderDef,
    b: &nacre_topo::CylinderDef,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let pa = circle_centre_rat(jd, wc, a).ok_or_else(undecided)?;
    let pb = circle_centre_rat(jd, wc, b).ok_or_else(undecided)?;
    let (ra, rb) = (a.radius(), b.radius());
    if ra >= rb {
        return Ok(false);
    }
    let dr = rb.checked_sub(ra).ok_or_else(undecided)?;
    let mut dist2 = nacre_scalar::Rat::from_int(0);
    for k in 0..3 {
        let d = pb[k].checked_sub(pa[k]).ok_or_else(undecided)?;
        dist2 = dist2
            .checked_add(d.checked_mul(d).ok_or_else(undecided)?)
            .ok_or_else(undecided)?;
    }
    let dr2 = dr.checked_mul(dr).ok_or_else(undecided)?;
    Ok(dr2 > dist2)
}

/// Which of `hosts` owns the hole: the **innermost** one — the candidate contained in all the
/// others.
///
/// Being inside several `+1` rings at once is ordinary two-level nesting, not a degeneracy: for
/// `A⁺ ⊃ D⁻ ⊃ B⁺ ⊃ c⁻`, `c` lies inside `B`'s ring *and* `A`'s, since `A`'s ring encloses its own
/// hole. Only `B` actually wraps `c` in material, and `emit_faces` must hang `c` off `B` — hanging
/// it off `A` would punch a hole through a face the hole is not even on.
///
/// **The innermost candidate always exists.** Cell rings are simple closed curves that do not cross
/// (a crossing would have been split into a node), so containment among them is a *total* order;
/// the candidates are a chain and its minimum is unique. Measured 2026-07-22 across the OCCT corpus:
/// every one of the 8 multi-host cases was a chain of exactly two. So the reject below is a net for
/// a broken invariant — if it ever fires, rings are crossing and the fault is upstream in
/// `split_at_crossings`, not here.
///
/// Containment is read the same way [`nest_cells`] reads it: one representative vertex through
/// [`combinatorics::point_in_ring`], and candidates sharing a node are adjacent rather than nested, so
/// they cannot be ordered and the honest answer is to reject.
fn innermost_host(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circles: &[MergedCircle],
    circle_ix: &[Option<usize>],
    hosts: &[usize],
) -> Result<usize, BoolError> {
    // ★ Two answers that used to be one. `None` meant *either* "adjacent, so not comparable" or
    // "every ray was spoiled", and both left as `HoleDepth` — whose own doc reads "if this fires,
    // rings are crossing and the fault is upstream", a diagnosis that is simply wrong for a
    // spoiled ray. Adjacency stays `None` here; an exhausted ring is now `NoClearRay`, its cause.
    //
    // ★★ **A disk can be a host, so this asks the same four-way question `nest_cells` does** —
    // it used to hold a polygon-only copy of one arm, which was sound only while a disk was
    // filtered out before it got here. It is [`cell_in_cell`] now.
    let inside = |a: usize, b: usize| cell_in_cell(jd, cyls, wc, rings, circles, circle_ix, a, b);
    let mut found = None;
    for &h in hosts {
        let mut wraps_all = true;
        for &o in hosts {
            if o != h && inside(h, o)? != Some(true) {
                wraps_all = false;
                break;
            }
        }
        if wraps_all {
            if found.is_some() {
                return Err(reject(RejectReason::HoleDepth)); // two minima: not a chain
            }
            found = Some(h);
        }
    }
    found.ok_or_else(|| reject(RejectReason::HoleDepth)) // no minimum: not a chain
}

/// A per-solid, per-side material label of one cell: `[A_above, A_below, B_above, B_below]`.
///
/// ★★ **"Above" is the side of the class's *stored* plane normal** — `trace_one`'s `w_normal`,
/// which is what every `body_above` in this file is measured against. Not the class's canonical
/// rational name (that points the other way on half the classes) and not the root face's outward
/// (`frame_sign` relates the two, and `emit`'s `flip` is what folds it back in). The absence of
/// this sentence is what let a lateral face's rim contribution be written in the wrong frame.
pub(crate) type Label = [bool; 4];

/// The flip mask an edge applies when crossed, grouping its `merged` contributions **per solid**.
/// Crossing the edge XORs this into the cell label.
///
/// **The circle around the edge.** An arrangement edge is the intersection of exactly two planes —
/// `W` and the wall — so the little circle around it has only two arcs to fill, above `W` and below.
/// A solid's material fills an arc or it does not, and crossing the edge inside `W` flips a label
/// exactly for the arcs the solid's boundary separates there:
///
/// ```text
///        above (n_s(W))          [T,F] one arc  → flip the above bit
///     ────────┼────────  W       [F,T] one arc  → flip the below bit
///        below                   [T,T] both     → flip both
/// ```
///
/// This is the BRep form of a Nef local pyramid (Hachenberger & Kettner, CGAL `Nef_3`) collapsed to
/// two planes, and the XOR itself is binary winding-number propagation (Zhou et al. 2016, "Mesh
/// Arrangements for Solid Geometry"; libigl `propagate_winding_numbers`).
///
/// Precedence per solid is **Graze > Seated > Transversal**:
/// - A `Graze` is a wall touching W along this edge with its body on one side — the solid's true
///   material boundary here. It overrides a coincident `Seated` cap-rim: the two agree on a convex
///   cap (same side) but disagree at a **reflex dihedral** in W, where the graze side is correct and
///   seated-wins would flip the wrong bit.
/// - A `Seated` face is ground truth where no graze coincides: it directly says which side the body
///   is on, so it **wins over a coincident `Transversal`** (the transversal's "flip both" assumes the
///   solid straddles W, which is false exactly where a cap crosses it).
/// - A lone `Transversal` genuinely straddles → flips both sides.
///
/// `> 1 Transversal`, disagreeing graze/seated sides, or a graze coincident with a same-solid true
/// crossing is a coincident-wall degeneracy outside the corpus → honest reject.
fn edge_mask(merged: &[(SolidSide, SegKind)]) -> Result<Label, BoolError> {
    let mut mask = [false; 4];
    for (solid, base) in [(SolidSide::A, 0usize), (SolidSide::B, 2)] {
        let kinds: Vec<SegKind> = merged
            .iter()
            .filter(|(s, _)| *s == solid)
            .map(|(_, k)| *k)
            .collect();
        let side_of_kind = |want_graze: bool| -> Vec<bool> {
            kinds
                .iter()
                .filter_map(|k| match k {
                    SegKind::Graze { body_above } if want_graze => Some(*body_above),
                    SegKind::Seated { body_above } if !want_graze => Some(*body_above),
                    // A tangent ruling flips its axis side — the seated rule with that side.
                    SegKind::Tangent { axis_above } if !want_graze => Some(*axis_above),
                    _ => None,
                })
                .collect()
        };
        let grazes = side_of_kind(true);
        let seated = side_of_kind(false);
        let transversals = kinds
            .iter()
            .filter(|k| matches!(k, SegKind::Transversal { .. }))
            .count();
        if !grazes.is_empty() {
            // Graze wins: it is the real boundary. A same-solid true crossing must not coincide.
            if transversals > 0 {
                // A graze coincident with a same-solid true crossing: the two disagree about
                // whether the solid straddles here, and nothing says which to believe.
                return Err(reject(RejectReason::EdgeOccupancyConflict));
            }
            // A valid solid brings exactly two faces to an edge, so more than two grazes on one
            // merged edge is degenerate input, not a case to arbitrate.
            if grazes.len() > 2 {
                return Err(reject(RejectReason::EdgeOccupancyConflict));
            }
            // ★ Grazes combine by **parity per side**, not by union. One graze is a step (the
            // face fills that arc — flip it). An opposite pair is a solid whose edge lies in `W`
            // with material above on one in-plane side and below on the other — both flip, the
            // `[T,T]` of the diagram, which is what a four-plane concurrency is made of. And a
            // **same-side pair flips nothing**: whether it is a knife edge (the wedge between the
            // two faces is the material, which pinches to measure zero at `W`) or its reflex
            // complement (the wedge is the void), what is immediately above `W` is the same on
            // both sides of this line. The union rule read that pair as one fill and flipped a
            // bit that changes nothing — a mask the propagation check then caught two layers
            // later as LabelConflict.
            for side in [true, false] {
                if grazes.iter().filter(|&&a| a == side).count() % 2 == 1 {
                    mask[base + usize::from(!side)] = true;
                }
            }
        } else if !seated.is_empty() {
            if seated.iter().any(|&b| b != seated[0]) {
                return Err(reject(RejectReason::EdgeOccupancyConflict)); // disagreeing seated sides
            }
            // seated wins over a coincident transversal: flip above if body_above, else below.
            mask[base + usize::from(!seated[0])] ^= true;
        } else if transversals == 1 {
            // pure transversal: the solid straddles W, flip both.
            mask[base] ^= true;
            mask[base + 1] ^= true;
        } else if transversals > 1 {
            return Err(reject(RejectReason::EdgeOccupancyConflict)); // >1 transversal, same solid
        }
        // no contributions ⇒ solid absent from this edge ⇒ no flip.
    }
    Ok(mask)
}

/// **The disk side of a cylinder's trace carries that cylinder's own solid** — `false` when it
/// does not, which is a labelling that cannot be right whatever the frame conventions are.
///
/// ★★★★★ **The absolute claim the arrangement had none of.** `label_cells` seeds the unbounded
/// contours void and flips per edge, then verifies **every** edge — and that verification is blind
/// to a global flip of one solid's two bits, because `l ^ mask == l'` survives flipping both sides
/// together. So a wrong seed inverts a whole component's material for one solid and every check
/// passes. Nothing in the crate stated an *absolute* fact about a label until this.
///
/// The fact: where the mask sets **both** of a solid's bits from a lone transversal
/// (`edge_mask`'s *"the solid straddles W, flip both"*), that solid's material is on the inside of
/// the cylinder exactly when the lateral face's outward is radially outward —
/// [`crate::planes::CylFaceInfo::orient_sign`] `= +1`, which is the `mat` the trace carries. So
/// the disk-side cell has the solid on **both** sides for a boss and on **neither** for a bore.
///
/// ★ It reads `edge_mask` rather than re-deriving the branch, so graze and seated keep their
/// priority; and it says nothing where both bits come from an **opposite graze pair** instead,
/// which carries no `mat` to compare against.
#[cfg(test)]
fn disk_side_agrees(merged: &[(SolidSide, SegKind)], label: &Label) -> Result<bool, BoolError> {
    let mask = edge_mask(merged)?;
    for (solid, base) in [(SolidSide::A, 0usize), (SolidSide::B, 2)] {
        if !(mask[base] && mask[base + 1]) {
            continue;
        }
        let Some(mat) = merged.iter().find_map(|(s, k)| match k {
            SegKind::Transversal { mat } if *s == solid => Some(*mat),
            _ => None,
        }) else {
            continue;
        };
        let want = mat > 0;
        if label[base] != want || label[base + 1] != want {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Label every cell by propagating from the unbounded contours across edges, flipping per
/// `edge_mask`. The cell interior cannot be point-queried (no plane-triple name), so propagation is
/// the only route. After propagating, **every** edge's flip
/// relation is verified (`label[c] XOR mask == label[neighbour]`); a violation means the trace was
/// incomplete and is an honest reject.
///
/// ★ **`seed` is the one place global information enters the arrangement.** Everything else here
/// is a flip relation between neighbours — purely local — so a plane's labelling is determined by
/// its own segments *plus* one known classification to start from. Whole-model arrangements pass
/// `[false; 4]`, which is the statement "the region outside every contour is void" and is true
/// exactly because the arrangement covers the whole plane: `Nesting`'s unbounded region reaches
/// infinity, where neither solid is. An arrangement restricted to a region of space cannot say
/// that — its unbounded region is an artifact of the restriction — and must be told instead. That
/// is the parameter's whole reason to exist; today's only caller still passes the old constant.
fn label_cells(
    cells: &[Cell],
    face_of: &HashMap<usize, usize>,
    edges: &ClassEdges<'_>,
    nesting: &Nesting,
    seed: Label,
) -> Result<Vec<Label>, BoolError> {
    // Crossing a circle flips like crossing any edge — its occupancy list is the mask's input;
    // the pseudo-half-edge numbering (`≥ 2·segs.len()`) picks the table.
    let mask_of = |he: usize| -> Result<Label, BoolError> {
        match edges.kind(he) {
            HalfEdgeKind::Seg(i) => edge_mask(&edges.segs[i].merged),
            // An arc is a piece of its circle's trace, so it carries the same contributions —
            // and a ruling piece is a piece of the lateral's trace, likewise.
            HalfEdgeKind::Arc(i) => edge_mask(&edges.arcs[i].merged),
            HalfEdgeKind::Ruling(i) => edge_mask(&edges.rulings[i].merged),
            HalfEdgeKind::Circle(i) => edge_mask(&edges.circles[i].whole_marks()?),
        }
    };
    // A face-with-holes is one region: label its group as a unit. Group representative → members.
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &g) in nesting.group_of.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut label = vec![None; cells.len()];
    let mut queue = std::collections::VecDeque::new();
    // Seed every unbounded contour's group, enqueuing every member. The outside region may be
    // bounded by more than one cycle, and they all describe the *same* region — so one seed serves
    // them all, whatever it says.
    for g in &nesting.root_groups {
        for &i in &members[g] {
            label[i] = Some(seed);
            queue.push_back(i);
        }
    }
    while let Some(c) = queue.pop_front() {
        let lc = label[c].unwrap();
        for &he in &cells[c].half_edges {
            let nb = face_of[&(he ^ 1)];
            if label[nb].is_none() {
                let mask = mask_of(he)?;
                let lab: Label = std::array::from_fn(|i| lc[i] ^ mask[i]);
                // A hole and its host bound the same region: label the whole group at once and
                // enqueue every member, so the hole cell's edges bridge to the cell inside it —
                // the two components share no edge, so nothing else reaches the interior one.
                for &mem in &members[&nesting.group_of[nb]] {
                    if label[mem].is_none() {
                        label[mem] = Some(lab);
                        queue.push_back(mem);
                    }
                }
            }
        }
    }
    let out: Vec<Label> = label
        .into_iter()
        .collect::<Option<_>>()
        .ok_or_else(|| reject(RejectReason::UnreachedCell))?; // a cell never reached
    // Verify every edge (tree and non-tree): the flip relation must hold everywhere.
    for (he, &c) in face_of {
        let nb = face_of[&(he ^ 1)];
        let mask = mask_of(*he)?;
        if std::array::from_fn::<bool, 4, _>(|i| out[c][i] ^ mask[i]) != out[nb] {
            return Err(reject(RejectReason::LabelConflict)); // inconsistent propagation
        }
    }
    Ok(out)
}

/// The boolean keep predicate on one chamber's `(inA, inB)`.
pub(crate) fn keep(kind: BoolKind, in_a: bool, in_b: bool) -> bool {
    match kind {
        BoolKind::Fuse => in_a || in_b,
        BoolKind::Cut => in_a && !in_b,
        BoolKind::Common => in_a && in_b,
    }
}

/// Emit the result faces on plane class `wc` for a boolean `kind`. A `+1` cell is a face of the
/// result iff its two chambers disagree under `keep` (material on one side of W, void on the
/// other); its `-1` holes (`nesting.holes`) ride along as inner rings. The DCEL cell ring is
/// already CCW about `n_out(wc)` (the walk stored `winding == +1`) and a hole cell is CW
/// (`winding == -1`) — exactly the `LocalFace.inner` contract ("kept material on the loop's left"),
/// so both are emitted verbatim; `flip` alone carries the chamber and `assemble_fuse_cut` reverses
/// outer and inner together, making the result normal point out of the kept solid:
/// `flip = keep_above == (orient_sign(wc) > 0)`.
///
/// Only `+1` cells are hosts — a `-1` cell is either a hole (emitted as some host's inner ring) or
/// the void root — so `-1` cells are skipped, never emitted as their own face.
///
/// **Output contract:** every ring vertex is a [`combinatorics::NodeId`] triple, including triples that coincide
/// with an original A/B vertex (a cap corner); the weld table canonicalizes such a W-triple onto
/// the same result vertex, or `assemble_fuse_cut`'s manifold guard rejects. This brick emits
/// all-`Seam`; the reconciliation and assembly are later.
#[allow(clippy::too_many_arguments)]
fn emit_faces(
    kind: BoolKind,
    labels: &[Label],
    cells: &[Cell],
    edges: &ClassEdges<'_>,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    holes: &HashMap<usize, Vec<usize>>,
) -> EmitOut {
    let planes = jd.planes;
    // [`combinatorics::NodeId`] is the vertex-identity enum; the local `Node` (this module's
    // three-valued-scan struct, a *different* type that happens to share the word) shadows the
    // glob-imported name here, which is why the identity is always spelled out in full.
    // ★ The wall travels with the ring. A half-edge was *told* which plane its edge rides, and
    // that is the one fact a name cannot always give back (see `boolean::Ring`). A circle cell
    // (pseudo-half-edge past the segment range) has no nodes — its boundary is the cylinder
    // class itself.
    let bound_of = |cell: &Cell| -> crate::boolean::Bound {
        if let Some(&he) = cell.half_edges.first()
            && let HalfEdgeKind::Circle(i) = edges.kind(he)
        {
            return crate::boolean::Bound::Circle {
                cyl: edges.circles[i].cyl,
            };
        }
        // ★★ `Ring.walls` is a carrier now (`boolean::Wall`), so an arc half-edge has something
        // true to put here at last: its cylinder class and which way this side travels. What is
        // still missing is downstream — `edge_for` refuses an arc carrier by the population's
        // name until the arc-casting cell teaches it the ordered circle key.
        crate::boolean::Bound::Ring(crate::boolean::Ring::new(
            cell.half_edges.iter().map(|&he| edges.origin(he)).collect(),
            cell.half_edges
                .iter()
                .map(|&he| match edges.kind(he) {
                    HalfEdgeKind::Seg(i) => crate::boolean::Wall::Plane(edges.segs[i].wall),
                    // `edge_at`'s own convention, carried not re-derived: `MergedArc::end` runs
                    // counter-clockwise about the axis, so the even half-edge travels that way
                    // and its twin the other.
                    HalfEdgeKind::Arc(i) => crate::boolean::Wall::Arc {
                        cyl: edges.arcs[i].cyl,
                        ccw: he % 2 == 0,
                    },
                    // Same shape for a ruling: `MergedRuling::end` ascends the axis, the even
                    // half-edge travels up. `edge_for` refuses this wall by the population's
                    // name (`RulingBoundNotYet`) until the panel cell teaches it the key.
                    HalfEdgeKind::Ruling(i) => {
                        let r = &edges.rulings[i];
                        crate::boolean::Wall::Ruling {
                            cyl: r.cyl,
                            side: r.side,
                            up: he % 2 == 0,
                        }
                    }
                    // Mirrors `origin`'s statement: the nodes map above already refused it.
                    HalfEdgeKind::Circle(_) => {
                        unreachable!("a circle's pseudo-half-edge has no vertex to leave")
                    }
                })
                .collect(),
        ))
    };
    let mut out = Vec::new();
    for (c, cell) in cells.iter().enumerate() {
        if cell.winding != 1 {
            continue; // a -1 cell is a hole or the void root, never a face of its own
        }
        let l = labels[c];
        let keep_above = keep(kind, l[0], l[2]);
        let keep_below = keep(kind, l[1], l[3]);
        if keep_above == keep_below {
            continue; // material the same on both sides ⇒ not a result face here
        }
        let flip = keep_above == (planes[wc].frame_sign > 0);
        let inner: Vec<crate::boolean::Bound> = holes
            .get(&c)
            .map(|hs| hs.iter().map(|&h| bound_of(&cells[h])).collect())
            .unwrap_or_default();
        out.push(LocalFace {
            surf: crate::planes::ClassIx::Plane(wc),
            outer: bound_of(cell),
            inner,
            flip,
        });
    }
    // ★★ **Every circle's disk label, whether or not a face was kept** (M6-2a K1). The four bits
    // of a disk cell say which solid's material lies immediately above and below this plane
    // *inside the circle* — which is exactly the chamber of the cylinder slab that starts here,
    // so the band pass reads them instead of casting a witness ray.
    //
    // ★ Collected from the **cells**, deliberately outside the keep filter above. A through
    // hole's outermost bands end on the cylinder's own cap classes, and `Cut` drops those cap
    // faces — following the faces would leave those bands with no label to read, while the cell
    // (and its label) is right here.
    let disk_labels = edges
        .circles
        .iter()
        .enumerate()
        .filter_map(|(i, mc)| {
            let he = edges.he_count() + 2 * i;
            let c = cells.iter().position(|cell| cell.half_edges == [he])?;
            #[cfg(test)]
            assert!(
                mc.whole_marks()
                    .ok()
                    .and_then(|w| disk_side_agrees(&w, &labels[c]).ok())
                    .unwrap_or(true),
                "a whole circle's disk does not carry its own solid: {:?}",
                labels[c]
            );
            Some((mc.cyl, labels[c]))
        })
        .collect();
    // ★ **A cut circle's labels, per arc** ([`ArcLabel`]) — same discipline as the disk labels
    // above: collected from the cells, outside the keep filter.
    //
    // ★★★★★ **Which half-edge borders the disk side, derived — and the factor that was missing.**
    // `MergedArc::end` runs counter-clockwise about the *axis*, so the even half-edge travels
    // `+θ̂`. The walk keeps a cell on the **left of its travel in the root face's frame**, whose
    // outward is `n_out = frame_sign · n_P` ([`crate::planes::WorkingPlane::frame_sign`] — the
    // sentence's one home; the ∥ road's `ruling_interior_is_even` cites it too), so
    //
    // ```text
    //   left = n_out × θ̂ = frame_sign·sign(n_P·m̂)·(m̂ × θ̂) = −frame_sign·axis_up·r̂
    //   ⇒ the even half-edge borders the disk  ⟺  axis_up · frame_sign = +1
    // ```
    //
    // The same product is already written one module over — `combinatorics`'s
    // `smooth_extremum_winding` computes `winding = ccw · axis_up · frame_sign`, and a circle
    // cell's boundary winds `+1` on the disk side, which is this statement rearranged.
    //
    // ★★★★★ **The rule this replaced was the same one missing `frame_sign`, and that is exactly
    // what refuted it.** It read `axis_up` alone. Measured 2026-08-31: a boss cutting a plate and
    // the **notch it leaves** put the same plane (z = 1), the same stored normal and the same two
    // arcs on **opposite** half-edges — the boss's top cap is `Forward` there and the notch's
    // ceiling is `Reversed`, so `frame_sign` is the one thing that differs. The label that came
    // back was the annulus cell's, all four bits true, which no disk-side cell there can be.
    // Over the lib suite the factor decides **112** of some 4,300 arcs, and it was wrong on every
    // one of them, silently: the only consumer is the chart's `read_cell`, which mostly refuses
    // before reading a cut end.
    //
    // ★★ **The geometry still watches it, on every arc** ([`disk_side_probe`]): a cell cannot
    // straddle the circle, so a rational corner's radial side names the side that whole cell is
    // on — and where a corner speaks it must agree with this rule. Measured 3,260 of 3,260 on the
    // census corpus. That check is what makes this a *derivation* rather than a third guess.
    let ns_arcs = 2 * edges.segs.len();
    let arc_labels = edges
        .arcs
        .iter()
        .enumerate()
        .filter_map(|(i, ma)| {
            // ★ `plus_t_is_above` and nothing spelled beside it: the inline `normal().dot(axis)`
            // that used to stand here was that function's second spelling, letter for letter.
            let axis_up = crate::planes::plus_t_is_above(&jd.planes[wc], &ma.def);
            let even = ns_arcs + 2 * i;
            let he = even + usize::from(axis_up != (jd.planes[wc].frame_sign > 0));
            #[cfg(test)]
            disk_side_probe::record(
                cells
                    .iter()
                    .position(|q| q.half_edges.contains(&even))
                    .and_then(|cx| cell_side(jd, edges, &cells[cx], &ma.def)),
                cells
                    .iter()
                    .position(|q| q.half_edges.contains(&(even + 1)))
                    .and_then(|cx| cell_side(jd, edges, &cells[cx], &ma.def)),
                he == even,
                jd.planes[wc].frame_sign,
            );
            let c = cells
                .iter()
                .position(|cell| cell.half_edges.contains(&he))?;
            #[cfg(test)]
            assert!(
                disk_side_agrees(&ma.merged, &labels[c]).unwrap_or(true),
                "{}: {:?} {:?}",
                disk_side_probe::NOT_OWN_SOLID,
                ma.merged,
                labels[c]
            );
            // ★ The label and the trace that made it, **from one visit to one arc** — see
            // [`ArcLabel`] for why they may not become two maps.
            Some((
                ma.cyl,
                ArcLabel {
                    ends: ma.end,
                    label: labels[c],
                    marks: ma.merged.clone(),
                },
            ))
        })
        .collect();
    (out, disk_labels, arc_labels)
}

/// [`emit_faces`]' product: the kept faces, the disk labels, and the cut circles' per-arc
/// disk-side labels.
type EmitOut = (Vec<LocalFace>, Vec<(usize, Label)>, Vec<(usize, ArcLabel)>);

/// A face's box, from its vertices' cached coordinates.
///
/// **`f64` is the right precision here and that is not a compromise.** Nothing this decides is
/// part of an answer — the culling spike below only *counts* what a box would skip.
#[cfg(test)]
fn face_box(model: &Model, fh: Handle<Face>) -> [[f64; 2]; 3] {
    let mut b = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    for p in face_points(model, fh) {
        for (k, s) in b.iter_mut().enumerate() {
            s[0] = s[0].min(p[k]);
            s[1] = s[1].max(p[k]);
        }
    }
    b
}

/// A face's vertex coordinates, from the caches.
///
#[cfg(test)]
fn face_points(model: &Model, fh: Handle<Face>) -> Vec<[f64; 3]> {
    let f = model.faces.get(fh);
    let mut out = Vec::new();
    for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
        for he in &lp.half_edges {
            for &vh in model.edges.get(he.edge).vertices.iter() {
                out.push(model.vertex_point(vh).as_array());
            }
        }
    }
    out
}

/// [`trace_result_faces`] for the band pass's tests, with the reuse mode fixed at `Proved` so
/// the cylinder guard inside is what turns it off.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_result_faces_full_for_test(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    n_a: usize,
    class_owner: &[Option<SolidSide>],
    trace_in: &combinatorics::TraceInput,
) -> Result<(Vec<LocalFace>, Curved, Option<BoolError>), BoolError> {
    trace_result_faces(
        model,
        kind,
        a,
        b,
        jd,
        faces,
        plane_ix,
        cyls,
        n_a,
        class_owner,
        crate::reuse::ClassReuse::Proved,
        trace_in,
    )
}

/// Per `(cylinder class, plane class)`, the four bits of that circle's disk cell — see the band
/// pass. Keyed rather than positional because a class carries a circle only when a cylinder cuts
/// it, and a cylinder cuts only some classes.
pub(crate) type DiskLabels = HashMap<(usize, usize), Label>;

/// **A cut circle's per-arc disk-side labels** — the sector sibling of [`DiskLabels`] (the
/// rulings ladder): a cut circle bounds no whole disk, so "the material immediately above and
/// below this plane inside the circle" is answered **per arc**, by the label of the cell on the
/// arc's disk side. Keyed like `DiskLabels`; each entry lists one [`ArcLabel`] per arc in the
/// split's own arc order. The frame is the disk labels' (the class's stored normal).
pub(crate) type ArcLabels = HashMap<(usize, usize), Vec<ArcLabel>>;

/// One arc's row in [`ArcLabels`] — what the arc bounds, and **who traced it**.
///
/// ★★★★★ **Two different questions live here and they are not one question.** `label` answers
/// *membership* — "is this side material" — and that is all a label can say. Whether this lateral
/// face is even **here** to bound it is a different proposition, and the trace answers that one:
/// `marks` is the contribution list that covered this arc ([`MergedArc::merged`]), where a face
/// running through says `Transversal` and one whose boundary stops at the arc says `Graze`.
/// A **holed** lateral re-entering a boolean makes the difference visible — its two sectors can
/// carry a literally identical label, and only the kind says which one is a face at all.
///
/// ★ The two are set together, in one pass over one arc, for the reason `tri_pt3` and `rotated`
/// are: split into two maps they could disagree, and then nothing could say which was the truth.
#[derive(Clone, Debug)]
pub(crate) struct ArcLabel {
    /// The arc's end pair, in the split's own arc order.
    pub(crate) ends: [NodeId; 2],
    /// The four bits of the cell on the arc's disk side, in the class's stored frame.
    pub(crate) label: Label,
    /// The `(solid, kind)` contributions covering this arc — [`MergedArc::merged`] verbatim.
    ///
    /// ★★★ **A lateral face's mark here is never `Seated`, and the type already says so.** The
    /// circle arm of [`trace_one`] builds its kind from [`CylOnClass`], whose only two answers are
    /// `Crosses` and `Grazes` — there is no `Seated` arm to take, because a cylinder's lateral
    /// surface cannot lie *in* a plane. `Seated` circle contributions come from the seated arm one
    /// branch further down, which a face enters only when its own class **is** this class: a
    /// planar face's rim. So "skip `Seated`" is not a guess about producers — it reads the one
    /// thing on this list that a lateral could not have written.
    pub(crate) marks: Vec<(SolidSide, SegKind)>,
}

/// **A cut circle's seam datum — carried from the split, never re-derived.** `split_circles`
/// already orders a cut circle's branch nodes by θ about the seam and classifies a seam-incident
/// node by name (`circular_order_about_seam`), so the one fact the assembly cannot re-derive
/// cheaply — *where θ = 0 sits among the arcs* — travels from the place that computed it.
#[derive(Clone, Debug)]
pub(crate) struct CutRim {
    /// The circle's branch nodes in θ order (CCW about the axis). When `seam_is_node`, the
    /// seam-incident node is first; otherwise θ = 0 lies inside the wrap arc
    /// `nodes.last() → nodes[0]`.
    pub(crate) nodes: Vec<combinatorics::NodeId>,
    pub(crate) seam_is_node: bool,
}

/// Per `(cylinder class, plane class)`, the cut circles — presence in this map **is** the one
/// source of "this rim is cut" (the assembly's rim skip and curved arms all read it).
pub(crate) type CutRims = HashMap<(usize, usize), CutRim>;

/// **One ruling piece, as the cylinder's own chart needs it** — the wall it rides and the
/// **axis interval** it spans.
///
/// ★★★★★ **Carried out rather than recomputed.** `rulings_on_class` already decides this extent
/// (the rulings ladder), and a chart that worked it out again would be a **second source of one
/// fact** — the defect shape this repository names first. So the value leaves by the road
/// [`CutRims`] already takes: made per class, folded once, keyed in the global class space.
///
/// ★ `end` is the piece's own two branch nodes and `z` their axis parameters
/// ([`node_axis_param`]). Both are carried because the chart needs the **order** (from the node's
/// `(MeetLine, QuadVal)` name, via `circular_order_about_seam`) and the **position** (`z`), and
/// deriving one from the other twice is how a name loses a sign.
/// ☑ Measured over the whole suite: **826 ruling pieces, every one with both ends named** — so a
/// piece whose ends have no ⊥ partner (which would have no `z` at all) is not a population today.
///
/// ★★★★ **There is no `side` here, and that is the point.** The first spelling carried one, and
/// when the cells arrived (D1b) nothing read it: a ruling's identity is `(wall class, root)`, and
/// which of the two parallel rulings that is derives from it through the one production spelling
/// ([`ruling_side`], which is how `cyl_chart::emit_lateral` gets it at emission time). A carried copy of a
/// derived value is the second spelling this ladder keeps being bitten by.
///
/// ★★★★★ **`cfg_attr(not(test), ...)`, not a blanket allow.** The only reader is `cyl_chart`,
/// which is `#[cfg(test)]`, so outside a test build these fields have none — but *inside* one they
/// must genuinely be read, and that is what caught `side`. A plain `#[allow(dead_code)]` would
/// have kept carrying it silently, which is how the field survived a whole rung.
#[derive(Clone, Debug)]
pub(crate) struct RulingExtent {
    pub(crate) wall: usize,
    /// Which of the wall's two rulings — or `0`, a **tangent** wall's single one (cell ⑩), the
    /// one station that carries no label: nothing changes across it, the face ends there.
    pub(crate) side: i8,
    pub(crate) end: [NodeId; 2],
    pub(crate) z: [nacre_scalar::Rat; 2],
    /// **The chart's vertical answer** (capability D, third rung): the label of the cell this
    /// ruling borders **inside** the cylinder, on the wall class this ruling lies on.
    ///
    /// ★★★★★ **This is the half `disk_labels`/`arc_labels` do not carry.** A [`Label`] holds the
    /// material on **both** sides of its plane, so a ⊥ class's disk label already answers the
    /// chart's *horizontal* crossings entire — which is why `cyl_chart::Chart::read_cell` reads two
    /// of them and asks them to agree. Crossing a *vertical* line is crossing the **wall**, and the
    /// cell inside the strip is where that plane's two sides are stated at the lateral.
    ///
    /// `None` when the wall's rational name and its stored normal cannot be related
    /// ([`world_rat_sense`]) — the ruling is then left without an answer and **counted**, never
    /// guessed at.
    ///
    /// ★ It was `#[cfg(test)]` through D3 — an instrument, because *"the emitter reads the
    /// horizontal lines only (a face has a circle at both ends)"*. **That day came**: a cut end
    /// the reader cannot pair with its rim leaves a cell with no horizontal answer at all, and
    /// this is what answers it ([`crate::cyl_chart::Chart::read_cell`]).
    pub(crate) label: Option<Label>,
    /// The `(solid, kind)` contributions that covered this piece — [`MergedRuling::merged`], the
    /// same list [`ArcLabel::marks`] carries for an arc.
    ///
    /// ★★★★★ **A label answers *membership*; this answers *existence*.** Cell ㉒ named that split
    /// on the ⊥ side: the chart collects **every** perpendicular class, and a ruling's extent is
    /// set by its *wall*, not by this lateral face — so a vertical line can be in the chart while
    /// the face is not there at all. `bands::face_spans` already states the rule for exactly this
    /// list, and it is carried beside the label for the reason `ArcLabel` carries both: split into
    /// two maps they could disagree, and then nothing could say which was the truth.
    #[cfg(test)]
    pub(crate) marks: Vec<(crate::planes::SolidSide, SegKind)>,
}

/// **Which of a ruling's two half-edges borders the cell inside the cylinder** — `true` for the
/// even one (the piece's own `end[0] → end[1]`, ascending the axis).
///
/// ★★★★★ **Derived from the definitions, not re-spelled.** Three sentences already in this crate
/// compose to the answer, and none of them is written a second time here:
///
/// * [`combinatorics::RulingCarrier::side`] is `sign((x − o) · (m̂ × n̂))` against the class's
///   **canonical rational** name, so the strip's interior lies along `−side · (m̂ × n̂_r)`; lifted
///   to the stored normal by [`world_rat_sense`] (`κ`) that is `−side·κ·(m̂ × n̂_P)` — the
///   product [`plus_theta_is_above`] already spells.
/// * The walk keeps a cell **on the left of its travel in the root face's outward frame**,
///   `n_out = frame_sign · n̂_P` ([`crate::planes::WorkingPlane::frame_sign`] — the sentence lives
///   there), and the even half-edge travels `+m̂` ([`MergedRuling::end`]) — so its cell lies
///   along `n_out × m̂ = −frame_sign · (m̂ × n̂_P)`.
///
/// ```text
///   even half-edge's cell is interior  ⟺  −frame_sign·(m̂ × n̂_P) ∥₊ −side·κ·(m̂ × n̂_P)
///                                      ⟺  side · κ · frame_sign = +1
///                                      ⟺  plus_theta_is_above ≠ (frame_sign > 0)
/// ```
///
/// — letter for letter the ⊥ road's rule for a cut circle's disk side (`axis_up ≠ (frame_sign >
/// 0)`, `emit_faces`), which is the other half-edge rule the chart frame enters.
///
/// ★★★★★ **The factor that was missing, and what found it.** This read `side · κ` alone until
/// cell ④, and the note here argued the product was frame-free — true of *which ruling* it names,
/// false of *which half-edge's cell*: the chart's left is the root face's, not the stored
/// normal's. The unmoved corpus never told them apart because a wall class is `frame_sign = −1`
/// only when its face lies on a **seed plane** with its outward along +axis (`Model::new` plants
/// x = 0, y = 0, z = 0 with cache direction −axis), which no fixture did until the commuting
/// oracle put the plate's max faces there by a translation and turned walls onto them by a
/// rotation — 105 cells, every one caught by `ruling_probe::SIDE_CHECK` at the fact.
fn ruling_interior_is_even(jd: &Judge<'_, WorkingPlane>, wc: usize, side: i8) -> Option<bool> {
    // The strip's interior is the side `+θ̂` does **not** enter, read in the chart's frame.
    Some(plus_theta_is_above(jd, wc, side)? != (jd.planes[wc].frame_sign > 0))
}

/// What the arrangement learned about the curved boundary, bundled: the band pass reads
/// `disk_labels`, the assembly reads `cut_rims`. One struct so the trace's return does not grow
/// element by element (it was widened once already, for `deferred`).
pub(crate) struct Curved {
    pub(crate) disk_labels: DiskLabels,
    pub(crate) arc_labels: ArcLabels,
    pub(crate) cut_rims: CutRims,
    /// ★ Cell ⑫: the alias table the class world settled on — every name in the labels, rims
    /// and rulings above is its representative, and the lateral chart, which mints station
    /// names of its own, asks it (`canon_point`) so a station on a corner *is* the corner.
    pub(crate) aliases: Aliases,
    /// Per cylinder class, the ruling pieces on it — the **vertical** lines of that cylinder's
    /// chart, where `disk_labels`/`arc_labels` carry the horizontal ones' answers.
    ///
    /// ★ **Its consumer is `cyl_chart`** — production since the D2b cutover: `emit_lateral` names
    /// every ruling from these (`Chart::node_on`, `Chart::ruling_name`), and the census reads their
    /// `label`/`marks`. Built as an instrument on capability D's first rung, measured against the
    /// two hand-written roads, then promoted.
    pub(crate) rulings: HashMap<usize, Vec<RulingExtent>>,
}

/// **The seam table — every node the result faces reference, realized to a coordinate and a
/// measured tolerance.** The weld table `assemble_fuse_cut` reads; built directly from the
/// emitted rings (no `build_seam`: that is raw-index and pierce-only), rejecting rather than
/// panicking on a degenerate meet.
///
/// ★ A named function rather than a block for the same reason `per_class` is one: the arc fence
/// calls it on the very faces production feeds it. The deferred stopper intercepts the whole
/// stretch this runs in, so a failure *here* never reaches an arc population's caller — which
/// means no reject name can testify that the branch arm works, and only a direct second consumer
/// can (measured: with the arm disabled wholesale, every boolean-level fence stays green).
pub(crate) fn seam_table(
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<SeamVertex>, BoolError> {
    watch!(SEAM);
    let geom = jd.planes;
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seen: HashMap<NodeId, ()> = HashMap::new();
    for f in faces {
        for loop_ in f.poly_rings() {
            for &node in loop_.iter() {
                if seen.insert(node, ()).is_some() {
                    continue;
                }
                // ★★ **A branch node is realized from its name, like everything else
                // here: the truth is the definition, the coordinate its cache.** The
                // coordinate is `a + b√c` — `branch_point` re-solves it from the name's
                // `(line, s)`; a rational road cannot hold it (`node_coords_rat`'s doc
                // calls that a type fact, not a width decline). The tolerance is the same
                // rule as the three-plane arm below: how far the realized point sits from
                // each surface that defines it, plus the closed-form pairwise meet
                // (`branch_vertex_tol` carries the argument for which pairwise curves are
                // in and out).
                //
                // ★ The **vertex minting** past this table is `boolean`'s
                // `def_triple`/`node_handle`, which names a branch node's vertex as
                // `VertexDef::Branch` — a cut rim's node and the scan's crossing on a
                // ruling both travel that road.
                if let Some(([p0, p1], cyl, _)) = combinatorics::branch_name(node) {
                    let wcy = &cyls[cyl];
                    let arr = combinatorics::branch_point(jd, cyl, &wcy.def, node)
                        .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                    let point = nacre_math::Point3::from_array(arr);
                    seam.push(SeamVertex {
                        point,
                        triple: node,
                        tol: crate::planes::branch_vertex_tol(
                            point,
                            &geom[p0].plane,
                            &geom[p1].plane,
                            &wcy.cache,
                        ),
                    });
                    continue;
                }
                // Declining, never `continue`: a skipped seam entry surfaces downstream as
                // `MissingSeam`, whose class is `SuspectedDefect` and whose sentence is
                // "a reconstruction dropped a crossing" — a wrong diagnosis for an input
                // the kernel simply does not build yet.
                let t = three_plane_name(node)
                    .ok_or_else(|| reject(RejectReason::BranchVertexUnnamed))?;
                let point = three_planes(&geom[t[0]].plane, &geom[t[1]].plane, &geom[t[2]].plane)
                    .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                seam.push(SeamVertex {
                    point,
                    triple: node,
                    tol: vertex_tol(
                        point,
                        &geom[t[0]].plane,
                        &geom[t[1]].plane,
                        &geom[t[2]].plane,
                    ),
                });
            }
        }
    }
    // **Two names, one point.** Every arrangement vertex is a distinct plane triple, and the
    // materialized coordinate is only its cache — so two *different* triples landing on the same
    // coordinate means the exact substrate and the f64 cache disagree about how many vertices
    // there are. Downstream that becomes a zero-length edge, so catch it here, where both triples
    // are still in hand, instead of letting `assemble_fuse_cut` discover it as a degenerate line.
    //
    // The usual cause is a **split plane table**: one geometric plane carried by two classes, whose
    // triples then name one point twice (measured 2026-07-22 — two `add_cuboid` walls at the same
    // x that `planes_coplanar` could not prove coplanar because their un-normalized coefficients
    // are not exactly proportional). A genuine 4-plane concurrency does the same.
    for (i, u) in seam.iter().enumerate() {
        for v in &seam[i + 1..] {
            if u.point == v.point {
                return Err(reject(RejectReason::SeamAlias));
            }
        }
    }
    Ok(seam)
}

/// Every result face across all plane classes, before assembly (the driver's risky half, testable
/// by face count without mutating the model). A declining class aborts the whole boolean.
///
/// Classes are visited in order, so the faces are produced in a sequence that is a function of the
/// input — which is what keeps `assemble_fuse_cut`'s handle minting replayable.
#[allow(clippy::too_many_arguments)]
fn trace_result_faces(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    n_a: usize,
    class_owner: &[Option<SolidSide>],
    reuse: crate::reuse::ClassReuse,
    trace_in: &combinatorics::TraceInput,
) -> Result<(Vec<LocalFace>, Curved, Option<BoolError>), BoolError> {
    let planes = jd.planes;
    let mut local_faces: Vec<LocalFace> = Vec::new();

    // **What each class contributes, for the classes the other operand cannot reach.** Everything
    // else stays `Arrange`, which is the whole engine as it was.
    // ★ **A cylinder in the operands turns the reuse shortcut off entirely** (M6-2a C4b).
    // `reuse::pass_through` moves faces **per plane class** (`plane_ix[fi].plane() != wc →
    // continue`), and a lateral face belongs to no plane class — so a class decided
    // `PassThrough` would carry the planes across and leave the cylinder's own faces behind,
    // silently. The guard also keeps `VertexClasses::of` (which walks face slots and asks each
    // for its plane row) away from cylinder rows, so it does double duty.
    let has_cyl = plane_ix.iter().any(|c| matches!(c, ClassIx::Cyl(_)));
    let reuse = if has_cyl {
        crate::reuse::ClassReuse::Off
    } else {
        reuse
    };
    let plans = crate::reuse::class_plans(model, reuse, kind, a, b, planes, class_owner);

    // ★ Two passes, because an identity must not depend on the order classes happen to be visited.
    // Pass A traces and splits every class, learning aliases as it goes; pass B builds the cells.
    // Doing both in one loop would key an early class's tables before a later class had reported
    // the alias that renames one of its vertices — the same point would then be assembled under two
    // names and rejected as a `SeamAlias`, for a reason that is really "we asked too early".
    //
    // `split_at_crossings`' output does not depend on the table (it only writes to it), so pass A
    // is a pure prefix of the old work rather than extra work — except that a *new* alias found
    // during a split leaves the merge that ran before it stale, so pass A repeats until the table
    // stops growing. Discoveries only accumulate and are bounded, so this terminates; a model with
    // no concurrency at all makes exactly one round.
    // ★ **Every class that has anything to arrange, in a fixed order.** Only the classes an operand
    // face lies on: a result face on `W` is part of `∂A` or `∂B`, so `W` carries an operand face.
    let work: Vec<usize> = {
        // ★ **Plane classes only.** A cylinder row rides in `trace_in` too (its lateral face is
        // what leaves circles on the ⊥ classes), and a lateral surface *is* not a plane class —
        // asking which one it is has no answer, which is exactly what `ClassIx::plane`'s panic
        // says. Filtering here is the upstream filter that panic is a detector for; without it
        // the first cylinder boolean past the C2 stopper aborts the kernel.
        let mut c: Vec<usize> = trace_in
            .faces
            .iter()
            .flatten()
            .filter_map(|(fp, _)| match plane_ix[*fp] {
                ClassIx::Plane(i) => Some(i),
                ClassIx::Cyl(_) => None,
            })
            .collect();
        c.sort_unstable();
        c.dedup();
        c
    };
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, jd, cyls, trace_in);
    #[allow(clippy::type_complexity)]
    let mut splits: Vec<(Vec<MergedSeg>, Vec<MergedCircle>, Vec<MergedRuling>)> = Vec::new();
    loop {
        let before = aliases.len();
        // **Every class in the round sees the table as it stood when the round began**, and
        // its own discoveries on top — where the sequential loop also showed it whatever the
        // lower-numbered classes had found meanwhile. The fixed point is the same, and for
        // two reasons that are both properties of `Aliases` rather than of the schedule:
        // a merged class's representative is its **minimum** element, so the final partition
        // does not depend on the order the unions happened in; and discoveries only
        // accumulate, so a round that learns something later than it used to just costs one
        // more round. The last round — the one whose splits are kept — runs on a table that
        // has stopped growing either way.
        let snapshot = aliases.clone();
        // Class-major, so the reject a decline raises is still the lowest-numbered class's.
        let round = crate::par::try_map_range(work.len(), |k| {
            let wc = work[k];
            // **Pass A runs for every class, including the ones pass B will not arrange.**
            //
            // It used to skip them, and that was unsound: pass B's reuse can *decline* — a vertex
            // where four planes meet has four possible names and only the alias table settles
            // which, so the class falls back to arranging — and the fallback reads `splits[wc]`,
            // which skipping pass A leaves empty. The class's faces would then vanish.
            //
            // Nothing in the code stopped that; it simply needed a model with a concurrency in a
            // region the other operand cannot reach, and the corpus has none. The cheap repair is
            // to keep the fallback a real one, which is what this does.
            let mut tr = timed!(
                TRACE_ON,
                trace_on_class(trace_in, wc, jd, cyls, faces, plane_ix, &snapshot)
            );
            let mut local = snapshot.clone();
            local.absorb(&std::mem::take(&mut tr.aliases));
            // An incomplete trace ⇒ honest reject, naming what the tracer could not do and on
            // which operand face. A class can decline several faces; the first is the one
            // reported, and `try_map_range` picks the lowest-numbered class, which is the one
            // the sequential loop returned at.
            if let Some(&(fp, kind)) = tr.declined.first() {
                // The witness is the face handle, read **kind-agnostically**: a lateral face can
                // decline too (`DeclineKind::CylSpan` is exactly that), and asking it for its
                // plane row would abort where an honest reject belongs.
                return Err(reject(decline_to_reject(kind, faces[fp].face())));
            }
            // ★ A disk cap's chord is one of these segments since cell ⑩ (both ends branch-pinned,
            // from the same parity sweep as every polygon's section) — it used to travel on a
            // vessel of its own and join here.
            let merged = timed!(MERGE, merge_coincident(jd, &tr.segs, wc, &local));
            let split = timed!(SPLIT, split_at_crossings(jd, cyls, wc, &merged, &mut local))?;
            let split = drop_newsless(split)?;
            let circles = merge_circles(&tr.circles, cyls, &local)?;
            let rulings = merge_rulings(&tr.rulings, cyls, &local);
            Ok(((split, circles, rulings), local))
        })?;
        splits.clear();
        for (split, local) in round {
            splits.push(split);
            // Absorbing the snapshot back is a no-op; only the round's discoveries are new.
            aliases.absorb(&local);
        }
        if aliases.len() == before {
            break;
        }
    }

    // Pass B is independent per class — it reads `splits[wc]` and the judging context, and
    // returns owned faces — so it is evaluated across cores. `try_map_range` is what keeps
    // that from being observable: the faces are consumed in class order, so the handles
    // `assemble_fuse_cut` mints are the ones a single thread would have minted, and a class
    // that declines surfaces the same rejection the sequential loop returned (the lowest
    // index, not whichever worker got there first).
    // The vertex→class map of each operand, built once and only if some class needs it.
    let vc_a = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::A))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, 0..n_a));
    let vc_b = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::B))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, n_a..faces.len()));

    let per_class = crate::par::try_map_range(splits.len(), |k| {
        let wc = work[k];
        let (split, circles, rulings) = &splits[k];
        // The per-class product: the faces, the disk labels the band pass reads (empty for a
        // class with no circles — and for a reused class, which is why cylinders switch reuse
        // off), the cut rims, and the **stopper socket's** deferred reject (empty since M6-2b
        // went green — see the socket note below).
        type ClassOut = (
            Vec<LocalFace>,
            Vec<(usize, Label)>,
            Vec<(usize, ArcLabel)>,
            Vec<(usize, CutRim)>,
            Vec<(usize, RulingExtent)>,
            Option<BoolError>,
        );
        let arrange = |wc: usize| -> Result<ClassOut, BoolError> {
            watch!(CELLS);
            // ★★ **One `ClassEdges` for the whole pipeline.** The split renumbers half-edges, so
            // everything below must read the *same* edges the walk did — building it here is what
            // makes that structural instead of a promise. (`frame_audit` runs its own copy of this
            // pipeline and must build it the same way; the arc fence locks that they agree.)
            let edges = timed!(
                C_SPLIT,
                ClassEdges::of(jd, cyls, wc, split, circles, rulings, &aliases)
            )?;
            // ★★★ **A stopper is *made* here — and only made.** M6-2b's arc stopper lived in
            // this socket: built where the class-level fact lives, raised inside `reconstruct`
            // at the assembly's very end, intercepting every stage between via the
            // `deferred.unwrap_or(e)` below — so a refused population carried one name out
            // however far the pipeline got. The population went green and the socket holds
            // `None`; the ladder stays for the next out-of-coverage class (M6-3's ellipses).
            //
            // ★★ `the_audit_and_the_boolean_agree_about_an_arc_class` compares the audit
            // *with* the boolean, so it stays green when both slide to the same wrong name; asking
            // it about interception measures the proposition next door.
            let staged = per_class(jd, cyls, kind, wc, &edges);
            // ★ **The stopper socket.** M6-2b's arc stopper was made here — per class, raised at
            // the assembly's very end, intercepting every stage between (seven layers of
            // `deferred.unwrap_or` down the whole pipeline). The arc population is supported
            // now, so the socket holds `None`; the next out-of-coverage class (M6-3's ellipses,
            // say) plugs its own deferred reject in here and inherits the entire interception
            // ladder instead of re-plumbing it.
            let deferred: Option<BoolError> = None;
            let s = match staged {
                Ok(s) => s,
                // The socket is empty, so clippy sees a literal `None` being unwrapped — the
                // yield's *shape* is the point (a plugged stopper wins here), kept as is.
                #[allow(clippy::unnecessary_literal_unwrap)]
                Err(e) => return Err(deferred.unwrap_or(e)),
            };
            Ok((
                s.faces,
                s.disk_labels,
                s.arc_labels,
                edges.cut_rims.clone(),
                s.ruling_extents,
                deferred,
            ))
        };
        // **The plan decides, and only ever downwards.** A `PassThrough` that cannot name one of
        // its vertices falls back to arranging, so this can lose the shortcut but never the answer.
        let reused = match plans[wc] {
            crate::reuse::ClassPlan::Arrange => None,
            crate::reuse::ClassPlan::Empty => Some(Vec::new()),
            crate::reuse::ClassPlan::PassThrough(side) => {
                watch!(REUSE);
                let (vc, range) = match side {
                    SolidSide::A => (vc_a.as_ref(), 0..n_a),
                    SolidSide::B => (vc_b.as_ref(), n_a..faces.len()),
                };
                vc.and_then(|vc| {
                    crate::reuse::pass_through(
                        model,
                        wc,
                        &planes[wc],
                        faces,
                        plane_ix,
                        range,
                        vc,
                        |t| aliases.canon_point(t),
                    )
                })
            }
        };
        match reused {
            Some(f) => Ok((f, Vec::new(), Vec::new(), Vec::new(), Vec::new(), None)),
            None => arrange(wc),
        }
    })?;
    let mut curved = Curved {
        disk_labels: HashMap::new(),
        arc_labels: HashMap::new(),
        cut_rims: HashMap::new(),
        rulings: HashMap::new(),
        aliases,
    };
    // The first arc class's deferred reject, in `work` order — the map above may run its classes
    // in parallel, but this fold reads the vec in order, so the choice is deterministic.
    let mut deferred: Option<BoolError> = None;
    for (k, (faces, labels, arcs, rims, ruls, d)) in per_class.into_iter().enumerate() {
        local_faces.extend(faces);
        // ★ `work[k]` is the translation from this arrangement's k-th class to the global plane
        // class index — the curved maps are keyed in the global space the assembly speaks.
        for (cyl, label) in labels {
            curved.disk_labels.insert((cyl, work[k]), label);
        }
        for (cyl, al) in arcs {
            curved
                .arc_labels
                .entry((cyl, work[k]))
                .or_default()
                .push(al);
        }
        for (cyl, rim) in rims {
            curved.cut_rims.insert((cyl, work[k]), rim);
        }
        // ★ `wall` was this arrangement's k-th class; restate it in the global space the rest of
        // the map already speaks, exactly as the three keys above do.
        for (cyl, mut r) in ruls {
            r.wall = work[k];
            curved.rulings.entry(cyl).or_default().push(r);
        }
        deferred = deferred.or(d);
    }
    Ok((local_faces, curved, deferred))
}

/// One arrangement vertex a boolean named, with **every** plane through it.
///
/// `planes` is ground truth: it is found by asking every plane in the table whether it passes
/// through the point, not by the rule the engine uses to notice concurrencies. That is the whole
/// value of it — the discovery rule can then be measured against something other than itself.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct Concurrency {
    /// The plane class whose arrangement used the name.
    pub wc: usize,
    /// The name that class used for the point.
    pub triple: NodeId,
    /// Every plane through that point, sorted. Longer than 3 exactly when the point is concurrent.
    pub planes: Vec<usize>,
    /// The sub-triples of `planes` that share a **line** rather than meeting at the point — the
    /// second face of the same degeneracy, and what a line-identity rule would have to fold. Read
    /// off `planes`, which is why one discovery answers both questions.
    pub lines: Vec<[usize; 3]>,
}

/// **Trace every class of the two solids' arrangement, past any decline** — for the ledgers.
///
/// The production driver returns at the first declined class, and in serial mode (no default
/// features) it never traces the classes after it; a lock that reads a road's ledger for a face
/// whose boolean is still refused one road later would then see different records per traversal
/// mode. This runs the tracer on every class, like [`concurrency_audit`], and keeps nothing but
/// what the roads recorded on the way.
#[cfg(test)]
pub(crate) fn trace_every_class(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(), BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings,
    );
    // The table production starts from (cell ⑫): the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    for wc in 0..geom.len() {
        let _ = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
    }
    Ok(())
}

/// Every **concurrent** vertex (four or more planes) the two solids' arrangement would name.
///
/// Runs the per-class front half — trace, merge, split — over **all** classes and keeps going past
/// a decline, which the production driver cannot do: it returns at the first declined class, so on
/// exactly the models that have concurrencies it would see one and stop. Both halves are called,
/// not reimplemented, so what is observed is what the engine does.
#[cfg(test)]
pub(crate) fn concurrency_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Concurrency>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut out = Vec::new();
    // The table production starts from (cell ⑫): the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        // Names this class used: segment endpoints, single-point touches, and — since a crossing
        // the arrangement mints is a vertex too — the split's endpoints where it got that far.
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let mut names: Vec<NodeId> = tr.touches.clone();
        for s in merged.iter().chain(
            split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default())
                .as_deref()
                .unwrap_or(&[])
                .iter(),
        ) {
            names.extend(s.end);
        }
        // ★ And the input the *trace's* rule actually sees: the operand faces' ring vertices.
        // Emitted segment endpoints are named `[wc, wall, r]`, so they always mention `wc` and can
        // never exercise the `wc ∉ t` branch — collecting only those measured nothing, which is
        // what the "was the rule exercised?" counter in the driver test caught.
        for (solid, inc) in [(a, &inc_a), (b, &inc_b)] {
            let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
                .into_iter()
                .flat_map(|sh| model.shells.get(sh).faces.clone())
                .collect();
            for fh in face_handles {
                let Some(&fp) = surf_ix.get(&fh) else {
                    continue;
                };
                let mut rings =
                    combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                        .map(|lr| lr.poly().map(|nr| nr.triples.clone()).unwrap_or_default())
                        .unwrap_or_default();
                rings.extend(
                    combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                        .unwrap_or_default()
                        .into_iter()
                        .filter_map(|lr| lr.poly().map(|nr| nr.triples.clone()))
                        .flatten(),
                );
                // Only vertices the class would actually name: those lying on it (`side == 0`),
                // which is exactly the run condition the trace's rule fires under.
                //
                // ★ The re-canonicalization that used to stand on both sides of this filter is
                // gone: a ring's nodes are `NodeId`s, and the only constructor sorts.
                // ★★ **The projection here was an `expect`, and it is gone with it.** A corner a
                // cylinder made is *excluded* — the honest answer for a four-plane concurrency
                // hunt — instead of panicking.
                // ★★★★★ **And the exclusion is said, not left to an argument value.** It used to
                // happen only because `side_of` was handed an empty cylinder table and so declined
                // a branch node; give that call a real table one day "for consistency" and the
                // `expect` below turns into a panic. A concurrency is a fact about plane triples,
                // so scope is the reason and it belongs in the filter — after which `&[]` is
                // provably never read, which is what the other two sites say by wrapping their own
                // triple.
                // ☑ Measured: it drops **0** nodes across the workspace suite — this audit's own
                // corpus is plane-only, so the filter is a precondition made explicit, not a
                // behaviour change.
                // ★ The call is one line because the source-text meta-test
                // `no_production_code_walks_a_ring_past_the_shared_walk` allow-lists this site by
                // its argument text, and it reads line by line.
                names.extend(
                    rings
                        .into_iter()
                        .filter(|&n| three_plane_name(n).is_some())
                        .filter(|&n| combinatorics::side_of(&jd, &[], n, wc) == Some(0)),
                );
            }
        }
        names.sort_unstable();
        names.dedup();
        for n in names {
            let t = three_plane_name(n).expect("a three-plane node");
            if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                continue; // names no point, so "the planes through it" is not a question
            }
            let mut planes: Vec<usize> = t.to_vec();
            planes.extend(
                (0..geom.len())
                    .filter(|q| !t.contains(q) && jd.orient3d(t[0], t[1], t[2], *q) == 0),
            );
            if planes.len() > 3 {
                planes.sort_unstable();
                let mut lines = Vec::new();
                for i in 0..planes.len() {
                    for j in (i + 1)..planes.len() {
                        for k in (j + 1)..planes.len() {
                            let tri = [planes[i], planes[j], planes[k]];
                            if jd.plane_pair_dir_sign(tri[0], tri[1], tri[2]) == 0 {
                                lines.push(tri);
                            }
                        }
                    }
                }
                out.push(Concurrency {
                    wc,
                    triple: n,
                    planes,
                    lines,
                });
            }
        }
    }
    Ok(out)
}

/// **One operand vertex, as every producer names it** (cell ⑪'s instrument).
///
/// A point's identity is meant to be a function of the set of planes through it. This report
/// holds, for one vertex of one operand, what the three sources say that set is and what names
/// come out: the **topology** (the plane classes of the faces incident to the vertex — what the
/// operand itself knows), the **geometry** (every class of the arrangement whose plane passes
/// through the point, asked plane by plane), and the **names** the ring road hands the tracer
/// for that vertex — one per face loop that visits it, with whether each is an independent
/// triple (three planes sharing a line name no point) and what the alias table folds it to after
/// every class has been traced.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct OperandVertexReport {
    /// Which operand (`0` = A, `1` = B) and the vertex.
    pub side: usize,
    pub vertex: Handle<Vertex>,
    pub point: [f64; 3],
    /// Plane classes of the incident faces, sorted; cylinder classes are not counted.
    pub topo: Vec<usize>,
    /// Every class through the point, from the first independent name; empty when no name is
    /// independent (then no point exists to ask about).
    pub geom: Vec<usize>,
    /// `(face class, name)` per loop visit; `dependent[i]` says the name's planes share a line.
    pub names: Vec<(usize, NodeId)>,
    pub dependent: Vec<bool>,
    /// `names[i]` after the alias fold of a full trace over every class.
    pub folded: Vec<NodeId>,
    /// Whether `folded[i]` names a point on every plane of `topo` — a fold that lands on another
    /// point is the silent shape of a naming defect (a dependent fold counts as off).
    pub folded_on_vertex: Vec<bool>,
}

/// **A branch corner of an operand seen from a class plane that passes through it** (cell ⑫'s
/// instrument). The point has more names than the corner's own: the three-plane name of its two
/// planes with the class, and the class's ruling crossing with the corner's cap plane at whichever
/// root is this point. Whether the alias table knows they are one point is what this reports.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct BranchCornerReport {
    pub side: usize,
    pub vertex: Handle<Vertex>,
    pub point: [f64; 3],
    /// The corner's own name (its two planes, its cylinder, its root).
    pub corner: NodeId,
    /// The class plane through the point that is not one of the corner's own.
    pub class: usize,
    /// The candidate names for the same point: the corner, the three-plane name, and the two
    /// ruling-crossing names (`side ±1`) of `class` with the corner's planes — one of the two is
    /// this point, the other the class's other ruling on the same cap.
    pub candidates: Vec<NodeId>,
    /// `candidates` after the alias fold of a full trace over every class.
    pub folded: Vec<NodeId>,
}

/// One declined (class, face, kind) of a full trace — [`trace_declines`]' row.
#[cfg(test)]
pub(crate) type ClassDecline = (usize, Option<Handle<Face>>, DeclineKind);

/// Every declined (class, face) of a full trace over every class — the audits' companion, since
/// the production driver stops at the first.
#[cfg(test)]
pub(crate) fn trace_declines(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<ClassDecline>, BoolError> {
    use crate::planes::{PlaneSetup, plane_index_setup};
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut out = Vec::new();
    // The table production starts from (cell ⑫): the operands' own concurrencies and corners.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        for &(fp, kind) in &tr.declined {
            out.push((wc, faces_tab[fp].face(), kind));
        }
    }
    Ok(out)
}

/// Every branch corner of `a` and `b` that lies on a class plane not its own, as
/// [`BranchCornerReport`]s — the cylinder twin of [`operand_vertex_audit`].
#[cfg(test)]
pub(crate) fn branch_corner_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<BranchCornerReport>, BoolError> {
    use crate::planes::{PlaneSetup, plane_index_setup};
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        aliases.absorb(&tr.aliases);
        let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
        let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
    }
    // The corners: every branch name a plane face's ring carries, once per vertex.
    let mut out = Vec::new();
    for (side, (solid, inc)) in [(a, &inc_a), (b, &inc_b)].into_iter().enumerate() {
        let mut seen: Vec<Handle<Vertex>> = Vec::new();
        let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
            .into_iter()
            .flat_map(|sh| model.shells.get(sh).faces.clone())
            .collect();
        for fh in face_handles {
            let Some(&fp) = surf_ix.get(&fh) else {
                continue;
            };
            let Ok(lr) =
                combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
            else {
                continue;
            };
            let Some(nr) = lr.poly() else {
                continue;
            };
            let face = model.faces.get(fh);
            if nr.triples.len() != face.outer.half_edges.len() {
                continue;
            }
            for (he, &n) in face.outer.half_edges.iter().zip(nr.triples.iter()) {
                let Some((planes, cyl, _root)) = combinatorics::branch_name(n) else {
                    continue;
                };
                let vh = crate::he_start(model, *he);
                if seen.contains(&vh) {
                    continue;
                }
                seen.push(vh);
                let def = &cyls[cyl].def;
                for class in 0..geom.len() {
                    if planes.contains(&class) {
                        continue;
                    }
                    if combinatorics::side_of(&jd, &cyls, n, class) != Some(0) {
                        continue;
                    }
                    let mut cands = vec![n];
                    let mut s = vec![planes[0], planes[1], class];
                    s.sort_unstable();
                    if let Some(t) = combinatorics::canonical_triple(&jd, &s) {
                        cands.push(NodeId::three_planes(t));
                    }
                    for fc in planes {
                        for side_r in [1i8, -1i8] {
                            // `crossing_on_ruling` takes the cap (normal ∥ axis) first and the
                            // class holding the axis second; the corner's wall plane errs out.
                            if let Ok(id) = crossing_on_ruling(&jd, def, fc, class, cyl, side_r) {
                                if !cands.contains(&id) {
                                    cands.push(id);
                                }
                            }
                        }
                    }
                    let folded = cands.iter().map(|&c| aliases.canon_point(c)).collect();
                    let p = model.vertex_point(vh);
                    out.push(BranchCornerReport {
                        side,
                        vertex: vh,
                        point: [p[0], p[1], p[2]],
                        corner: n,
                        class,
                        candidates: cands,
                        folded,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// Every plane-only operand vertex of `a` and `b`, as [`OperandVertexReport`]s — cell ⑪'s
/// audit of the ring road's naming rule against the topology and the geometry.
///
/// Like [`concurrency_audit`] it drives the real front half (trace, merge, split) over **every**
/// class and keeps going past a decline, accumulating the alias table so the fold it reports is
/// the one production would reach. Vertices a cylinder touches are skipped: their names are
/// branch points, a different vocabulary.
#[cfg(test)]
pub(crate) fn operand_vertex_audit(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<OperandVertexReport>, BoolError> {
    use crate::planes::{ClassIx, PlaneSetup, plane_index_setup};
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        crossings,
        cyls,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    // The alias fold production would reach: every class traced, merged and split, discoveries
    // accumulated across classes.
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, &jd, &cyls, &trace_in);
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        aliases.absorb(&tr.aliases);
        let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
        let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
    }
    let plane_class = |k: usize| match plane_ix[k] {
        ClassIx::Plane(c) => Some(c),
        ClassIx::Cyl(_) => None,
    };
    let mut out = Vec::new();
    for (side, (solid, inc)) in [(a, &inc_a), (b, &inc_b)].into_iter().enumerate() {
        // Topology: vertex → plane classes of its incident faces; a vertex on any cylinder face
        // is left out (its name is a branch point).
        let mut topo: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
        let mut curved: std::collections::HashSet<Handle<Vertex>> = Default::default();
        for (bounds, pair) in inc.edges() {
            for &vh in bounds {
                for &k in pair {
                    match plane_class(k) {
                        Some(c) => topo.entry(vh).or_default().push(c),
                        None => {
                            curved.insert(vh);
                        }
                    }
                }
            }
        }
        // Names: per face loop, the ring road's triples, position `i` being the start of
        // half-edge `i` (a seam joint pushes no name, so a loop whose lengths differ is skipped).
        let mut names: HashMap<Handle<Vertex>, Vec<(usize, NodeId)>> = HashMap::new();
        let face_handles: Vec<Handle<Face>> = solid_shell_handles(model, solid)
            .into_iter()
            .flat_map(|sh| model.shells.get(sh).faces.clone())
            .collect();
        for fh in face_handles {
            let Some(&fp) = surf_ix.get(&fh) else {
                continue;
            };
            let Some(fc) = plane_class(fp) else {
                continue;
            };
            let face = model.faces.get(fh);
            let mut loops: Vec<(Vec<nacre_topo::HalfEdge>, Option<Vec<NodeId>>)> = Vec::new();
            loops.push((
                face.outer.half_edges.clone(),
                combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix, &cyls)
                    .ok()
                    .and_then(|lr| lr.poly().map(|nr| nr.triples.clone())),
            ));
            if let Ok(holes) = combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix, &cyls)
            {
                for (lp, lr) in face.inner.iter().zip(holes) {
                    loops.push((
                        lp.half_edges.clone(),
                        lr.poly().map(|nr| nr.triples.clone()),
                    ));
                }
            }
            for (hes, triples) in loops {
                let Some(triples) = triples else {
                    continue;
                };
                if triples.len() != hes.len() {
                    continue;
                }
                for (he, &n) in hes.iter().zip(triples.iter()) {
                    let vh = crate::he_start(model, *he);
                    names.entry(vh).or_default().push((fc, n));
                }
            }
        }
        let mut vertices: Vec<Handle<Vertex>> = topo.keys().copied().collect();
        vertices.sort_by_key(|v| v.index());
        for vh in vertices {
            if curved.contains(&vh) {
                continue;
            }
            let mut t = topo[&vh].clone();
            t.sort_unstable();
            t.dedup();
            let vnames = names.get(&vh).cloned().unwrap_or_default();
            let dependent: Vec<bool> = vnames
                .iter()
                .map(|(_, n)| match three_plane_name(*n) {
                    Some(tr) => jd.plane_pair_dir_sign(tr[0], tr[1], tr[2]) == 0,
                    None => false,
                })
                .collect();
            let geom_set =
                vnames
                    .iter()
                    .zip(&dependent)
                    .find(|(_, d)| !**d)
                    .and_then(|((_, n), _)| three_plane_name(*n))
                    .map(|tr| {
                        let mut g: Vec<usize> = tr.to_vec();
                        g.extend((0..geom.len()).filter(|q| {
                            !tr.contains(q) && jd.orient3d(tr[0], tr[1], tr[2], *q) == 0
                        }));
                        g.sort_unstable();
                        g
                    })
                    .unwrap_or_default();
            let folded: Vec<NodeId> = vnames
                .iter()
                .map(|(_, n)| aliases.canon_point(*n))
                .collect();
            let folded_on_vertex = folded
                .iter()
                .map(|n| match three_plane_name(*n) {
                    Some(tr) => {
                        jd.plane_pair_dir_sign(tr[0], tr[1], tr[2]) != 0
                            && t.iter().all(|&q| {
                                tr.contains(&q) || jd.orient3d(tr[0], tr[1], tr[2], q) == 0
                            })
                    }
                    None => false,
                })
                .collect();
            let p = model.vertex_point(vh);
            out.push(OperandVertexReport {
                side,
                vertex: vh,
                point: [p[0], p[1], p[2]],
                topo: t,
                geom: geom_set,
                names: vnames,
                dependent,
                folded,
                folded_on_vertex,
            });
        }
    }
    Ok(out)
}

/// One plane class's **label-frame audit** (family #2 diagnostic): what each producer says about
/// "above" on this class, plus how far the per-class pipeline gets. `#[cfg(test)]`, `pub(crate)`
/// so the driver test can live in `crate::tests` where the two-solid fixtures are (the same reason
/// [`boolean`] is `pub(crate)`).
#[cfg(test)]
#[derive(Debug)]
pub(crate) struct ClassAudit {
    pub wc: usize,
    /// The root plane resolved to real coordinates — a point on it and its stored normal. Plane
    /// indices alone have twice carried a wrong geometric story into a write-up.
    pub root_point: [f64; 3],
    pub root_normal: [f64; 3],
    /// `sign(stored normal · n_out)` — `-1` exactly when the class root face is `Reversed`, which
    /// is when the stored-normal label frame and `side_of`'s outward frame disagree.
    pub orient_sign: i8,
    /// `body_above` of every seated segment, and `body_above` of every graze segment.
    pub seated: Vec<bool>,
    pub grazes: Vec<bool>,
    pub transversals: usize,
    pub declined: Vec<(usize, DeclineKind)>,
    /// The reason the per-class pipeline raised, if any (`None` = the class went through).
    pub failed_at: Option<RejectReason>,
    /// **What the arrangement produced, before anything refused it** — cells walked, arcs the
    /// split cut, and the nesting's two counts.
    ///
    /// ★★★ Born while the arc stopper swallowed every arrangement stage's answer (a probe was
    /// the only way to see the population, and a probe is deleted before the commit); the
    /// population is green now, and this stays as the audit's direct read of what a class
    /// produced. `None` when the class stopped before those stages ran.
    pub produced: Option<Produced>,
    /// **Every emitted face's outer ring, as coordinates** — one `Vec` per emitted face whose outer
    /// bound is a polygon ring (a face bounded by an uncut circle has no nodes and contributes
    /// nothing), sorted by their rotation-normalized first coordinate. `None` when the class
    /// stopped before `emit_faces`.
    ///
    /// ★★★★★ **Coordinates, not node identities, and the reason is that the lock must not be
    /// circular.** A node name is `ThreePlane([0, 2, 5])` — plane **class** indices, which no
    /// fixture derives; writing the expected value means running the engine and copying what it
    /// said, and a test whose oracle is its subject measures nothing. Coordinates come from the
    /// fixture's own numbers.
    ///
    /// ★★★★★ **Rotation is normalized, reversal deliberately is not.** The walk tries both
    /// handednesses and pushes cells in its own order, so *which* vertex a ring starts at is not a
    /// fact about the geometry — the minimum coordinate goes first. Reversal *is* a fact: it is
    /// what a wrong `sense` on a split segment produces, and it is invisible in every
    /// order-independent summary the previous rungs measured (cell counts, nesting, sorted
    /// labels). Folding it here would delete the one thing this field exists to see.
    ///
    /// ★★★ **And rotation-normalization is itself reversal-blind at `n <= 2`**: reversing a
    /// two-cycle *is* a rotation of it. Such a ring cannot carry this lock at all — which is why
    /// the fence's red probe is read on the turned boss (5- and 3-rings), not the straddling one
    /// (6- and **2**-rings).
    ///
    /// ★ It is a sibling of [`Produced`] rather than a field in it because `[f64; 3]` is not `Eq`
    /// (the derive would break) and one of these coordinates is irrational, so it needs `near()`
    /// rather than `==` — the two could not ride the same `assert_eq!` regardless.
    ///
    /// ★★★ **Two things are left out on purpose, and saying so is the point** — a silent omission
    /// reads as an oversight to whoever comes next:
    /// - **`flip`.** It is `keep_above == (frame_sign > 0)`, so it turns on the *class's stored
    ///   normal direction* — an implementation fact, not one of the fixture's numbers. Carrying it
    ///   would put one bit in here that a fence could only fill by copying a run, and reversal —
    ///   the thing this field exists for — is already caught by the sequence itself.
    /// - **Inner rings.** A face's holes are not carried. Today's arc classes have none
    ///   (`nesting.holes` is empty, measured), so there is nothing to lose *yet*; an arc class with
    ///   a hole would go unchecked here, and that is the day this grows.
    pub outer_rings: Option<Vec<Vec<[f64; 3]>>>,
}

/// The per-class arrangement's output, for [`ClassAudit`].
#[cfg(test)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Produced {
    pub cells: usize,
    pub arcs: usize,
    pub roots: usize,
    pub holes: usize,
    /// Every `+1` cell's label, **sorted** — the labelling stage's answer, in a form the walk's
    /// handedness cannot reorder.
    ///
    /// ★★ `Copy` goes with this, and that was the choice rather than the discovery: a sibling field
    /// on [`ClassAudit`] would keep it, but then one class's two facts live in two places and the
    /// fence has to join them by index — a seam where they can drift. `PartialEq` stays, so the
    /// fence remains a single `assert_eq!` against a value derived from the geometry.
    pub pos_labels: Vec<Label>,
}

/// Lexicographic order on realized coordinates. `total_cmp` rather than `partial_cmp`: a node that
/// failed to realize comes through as `NaN` below, and this must still be a total order.
#[cfg(test)]
fn cmp_pt(a: &[f64; 3], b: &[f64; 3]) -> std::cmp::Ordering {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.total_cmp(y))
        .find(|o| o.is_ne())
        .unwrap_or(std::cmp::Ordering::Equal)
}

/// Every emitted face's outer ring as coordinates, for [`ClassAudit::outer_rings`] — that field
/// carries the argument for coordinates, for normalizing rotation, and for leaving reversal alone.
///
/// ★★ **A ring here can need both roads to a coordinate** — three of its nodes are plane triples
/// and two are branch points — and the dispatch between them lives in
/// [`combinatorics::node_point_f64`], not here: spelling a [`combinatorics::NodeId`] variant
/// outside that file is what `three_plane_name`'s gate forbids, and the first draft of this
/// function broke it with the suite green.
///
/// ★ A node that cannot be realized becomes `NaN`, not a dropped element: the ring keeps its
/// length, so the fence reads "this vertex had no coordinate" instead of "the ring is short".
#[cfg(test)]
fn face_ring_coords(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
) -> Vec<Vec<[f64; 3]>> {
    let point = |n: NodeId| combinatorics::node_point_f64(jd, cyls, n).unwrap_or([f64::NAN; 3]);
    let mut out: Vec<Vec<[f64; 3]>> = faces
        .iter()
        .filter_map(|f| f.outer.ring())
        .map(|r| {
            let pts: Vec<[f64; 3]> = r.nodes.iter().map(|&n| point(n)).collect();
            // Rotation only. Where the walk started is not a fact about the geometry; which way it
            // went is, and it is the one this instrument was added to see.
            match pts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| cmp_pt(a, b))
                .map(|(i, _)| i)
            {
                Some(i) => pts[i..].iter().chain(&pts[..i]).copied().collect(),
                None => pts,
            }
        })
        .collect();
    // Between faces the order *is* the walk's, so it is normalized — lexicographically over the
    // whole (already rotated) sequence.
    //
    // ★ **Total on purpose, not just "by the first vertex".** Two faces of one class can share
    // their minimum vertex — they meet there — and ordering on that alone would leave such a pair
    // in walk order, which is a handedness coin flip inside an instrument whose whole job is to be
    // stable. Today's two fixtures have distinct minima, so this costs nothing and removes the
    // class of flake rather than relying on the fixture to avoid it.
    out.sort_by(|a, b| {
        a.iter()
            .zip(b)
            .map(|(x, y)| cmp_pt(x, y))
            .find(|o| o.is_ne())
            .unwrap_or_else(|| a.len().cmp(&b.len()))
    });
    out
}

/// Audit every plane class of one boolean input pair. Runs the same per-class pipeline as
/// [`trace_result_faces`] but never aborts, so one failing class does not hide the rest.
#[cfg(test)]
pub(crate) fn frame_audit(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<ClassAudit>, BoolError> {
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
        cyls,
        crossings,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    let jd = Judge::new(&geom, standard, &notes);
    let trace_in = combinatorics::trace_input(
        model,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        crossings.clone(),
    );
    // The alias fixpoint, exactly as the boolean runs it (sequentially): classes discover names
    // for one another's features, so auditing each class against an empty table is auditing a
    // *different* pipeline — it diverged from the boolean's answer the day the rounds arrived,
    // and the divergence surfaced when this file's test moved to a fixture that needs them.
    let mut aliases = Aliases::default();
    loop {
        let before = aliases.len();
        for wc in 0..geom.len() {
            let mut tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
            aliases.absorb(&std::mem::take(&mut tr.aliases));
            if tr.declined.is_empty() {
                let merged = merge_coincident(&jd, &tr.segs, wc, &aliases);
                let _ = split_at_crossings(&jd, &cyls, wc, &merged, &mut aliases);
            }
        }
        if aliases.len() == before {
            break;
        }
    }
    let mut out = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &cyls, &faces_tab, &plane_ix, &aliases);
        let side = |f: fn(&SegKind) -> Option<bool>| -> Vec<bool> {
            tr.segs.iter().filter_map(|s| f(&s.kind)).collect()
        };
        let mut audit = ClassAudit {
            wc,
            root_point: geom[wc].tri[0].as_array(),
            root_normal: geom[wc].plane.normal().as_array(),
            orient_sign: geom[wc].frame_sign,
            seated: side(|k| match k {
                SegKind::Seated { body_above } => Some(*body_above),
                _ => None,
            }),
            grazes: side(|k| match k {
                SegKind::Graze { body_above } => Some(*body_above),
                _ => None,
            }),
            transversals: tr
                .segs
                .iter()
                .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
                .count(),
            declined: tr.declined.clone(),
            failed_at: None,
            produced: None,
            outer_rings: None,
        };
        // Run the rest of the per-class pipeline, recording where it stops.
        audit.failed_at = if let Some(&(fp, kind)) = audit.declined.first() {
            Some(decline_to_reject(kind, faces_tab[fp].face()))
        } else {
            let mut produced = None;
            let mut outer_rings = None;
            let mut run = || -> Result<(), BoolError> {
                // A copy, so one class's split discoveries cannot leak into the next class's
                // audit — the fixpoint above already holds everything the boolean would know.
                let mut local = aliases.clone();
                let merged = merge_coincident(&jd, &tr.segs, wc, &local);
                let split = split_at_crossings(&jd, &cyls, wc, &merged, &mut local)?;
                let split = drop_newsless(split)?;
                let circles = merge_circles(&tr.circles, &cyls, &local)?;
                let rulings = merge_rulings(&tr.rulings, &cyls, &local);
                // ★★★★ **The same edges the boolean uses.** This copy of the pipeline is what
                // makes the audit an instrument; feeding it un-split edges would let it run to
                // the end and report `failed_at: None` for an input the boolean refuses —
                // *"the worst possible time to be lying"* (`decline_to_reject`). The arc fence
                // in `bands.rs` locks the two together — since M6-2b went green, as «both
                // succeed» (a stopper plugged into the socket would raise here per class while
                // the boolean defers it to the assembly's end; same classes, same name).
                let edges = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &rulings, &local)?;
                let staged = per_class(&jd, &cyls, kind, wc, &edges);
                if let Ok(s) = &staged {
                    let mut pos_labels: Vec<Label> = s
                        .cells
                        .iter()
                        .zip(&s.labels)
                        .filter(|(c, _)| c.winding == 1)
                        .map(|(_, &l)| l)
                        .collect();
                    pos_labels.sort_unstable();
                    produced = Some(Produced {
                        cells: s.cells.len(),
                        arcs: edges.arcs.len(),
                        roots: s.nesting.root_groups.len(),
                        holes: s.nesting.holes.len(),
                        pos_labels,
                    });
                    outer_rings = Some(face_ring_coords(&jd, &cyls, &s.faces));
                }
                let _ = staged?;
                Ok(())
            };
            let stopped = run().err().and_then(|e| match e {
                BoolError::Rejected { reason, .. } => Some(reason),
                // No live-set check runs inside the pipeline, so this arm is unreachable.
                BoolError::InputNotLive => None,
            });
            audit.produced = produced;
            audit.outer_rings = outer_rings;
            stopped
        };
        out.push(audit);
    }
    Ok(out)
}

/// Drive the arrangement pipeline over **every** plane class and assemble the result solid — the
/// engine the public `crate::boolean` delegates to.
///
/// Vertex welding is automatic: every result vertex is a sorted triple of three canon plane
/// classes (`canon3`), so a corner shared by three planes gets one identical [`combinatorics::NodeId`]
/// regardless of which plane was the cut class W — `assemble_fuse_cut` welds them to one vertex.
/// A four-plane concurrency would name it inconsistently, which is what the alias table settles
/// and what `boolean::Ring`'s carried walls keep out of the naming in the first place.
///
/// A class that declines (holes, degenerate) aborts the whole boolean: skipping it would drop real
/// faces and silently produce a non-manifold or wrong-volume solid.
/// What each input face's plane became: its **plane class's representative surface**, which is
/// what every result face on that plane carries.
pub(crate) type ClassOf = std::collections::HashMap<Handle<Face>, Handle<nacre_geom::Surface>>;

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
        // cell — what the band pass reads instead of casting a witness ray (M6-2a K1).
        // `deferred`: an arc class's stopper reject, made per class and raised below **after** the
        // seam stretch, so the seam's branch arm runs before the population is refused.
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
        // whole-boolean level, because per class there is no longer an arrangement to compare to.
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
        // an arc input (that is the point: the seam's branch arm is exercised), and then the
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
                crate::boolean::unify_coplanar_faces(faces, &jd, &cyls)
            )?;

            // ★ **The lateral bands, appended after the differential above** (M6-2a C4b): reuse can
            // only change what the *plane* arrangement emits, so the two routes are compared on that
            // list; the bands are a separate pass over the same operands and belong to neither route.
            // From here on there is one face list — the grouping, the closed-shell guard and the
            // assembly all read it.
            let faces = if cyls.is_empty() {
                faces
            } else {
                let rows = crate::bands::cyl_rows(&faces_tab, &plane_ix, n_a)?;
                // ★ **The lateral faces come from the chart** (D2b cutover): every cylinder
                // class's cells, read off the plane arrangement's own labels and emitted in the
                // band road's vocabulary (the road itself was deleted in D3).
                let lateral =
                    crate::cyl_chart::emit_lateral(kind, &jd, &cyls, &faces, &curved, &rows);
                // ★★ The census sits *before* the curved cleaning pass on purpose: the cleaning
                // merges pieces, which would blur the face-by-face question being asked.
                // ★ The census holds the chart against its own rules (and the emitter against
                // the census's independent count), in test builds; the emitter's refusal is
                // handed to it before `?` decides, so a refused class is still recorded.
                #[cfg(test)]
                crate::cyl_chart::census(&jd, &cyls, kind, &faces, &curved, &rows, &lateral);
                let mut faces = faces;
                faces.extend(lateral?);
                // ★ No curved cleaning pass follows (D5): the emitter's lateral faces are the
                // regions of each chart already, so there is no phantom seam left to erase — the
                // pass that used to run here measured `merged 0` over the suite and was deleted.
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
                crate::boolean::Tangencies {
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
    // if it were the problem. That is the `LoopOrientMismatch`-hiding-precision-exhaustion trap,
    // and checking the evidence *before* returning the symptom is what keeps it shut.
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
mod tests {
    use super::*;
    use crate::SketchFrame;
    use nacre_scalar::Axis;

    /// ★ **A fixture with no cylinders, said as a fact rather than left as a hole.** The table is
    /// how a `NodeId::Branch` reaches its definition, so an all-plane fixture has nothing to put
    /// in it — and a bare `&[]` at a call site reads like something forgotten.
    const NO_CYLS: &[crate::planes::WorkingCyl] = &[];

    /// ★★★ **Where a named curved ring stops, and under which name.**
    ///
    /// ★★★★★ **The vessel carries a curved ring through, and only a *collapsed* name is refused.**
    ///
    /// This used to assert the opposite: a branch corner was `RingFail::Branch` and a curved
    /// carrier `RingFail::CurvedWall`, because the ring came out as plane ids and could not hold
    /// either. Both refusals were the type running out rather than a decision, and they are gone
    /// with the projections. What still stops here is a name that is degenerate *as a name* — two
    /// of its three planes equal — which is a fact about the triple and not about cylinders.
    ///
    /// ★ Rewritten rather than patched: the old proposition ("the plane road cannot use such a
    /// ring") became false, and a lock whose sentence is false is worse than no lock.
    #[test]
    fn a_named_curved_ring_rides_through_and_only_a_collapsed_name_stops() {
        let three =
            |a, b, c| combinatorics::NodeId::three_planes(combinatorics::Canon3::three([a, b, c]));
        let branch = combinatorics::NodeId::branch(0, 1, 0, nacre_topo::QuadRoot::Lo);
        let plane = |c| crate::boolean::Wall::Plane(c);
        let ruling = crate::boolean::Wall::Ruling {
            cyl: 0,
            side: 1,
            up: true,
        };
        // A curved carrier rides through, carried as itself.
        let curved_carrier = combinatorics::NamedRing {
            triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
            walls: vec![plane(1), ruling, plane(3)],
            arc_ccw: vec![None; 3],
            concurrencies: vec![],
        };
        let (_, walls) = plane_ring(&curved_carrier).expect("a curved carrier is describable");
        assert_eq!(
            walls[1], ruling,
            "the carrier is the producer's, unflattened"
        );
        // So does a branch corner, as its own name.
        let curved_corner = combinatorics::NamedRing {
            triples: vec![three(0, 1, 2), branch, three(0, 1, 3)],
            walls: vec![plane(1), ruling, plane(3)],
            arc_ccw: vec![None; 3],
            concurrencies: vec![],
        };
        let (ts, _) = plane_ring(&curved_corner).expect("a branch corner is describable");
        assert_eq!(ts[1], branch, "the corner is the producer's, unflattened");
        // ★ What is still refused, and for a reason that has nothing to do with cylinders.
        let collapsed = combinatorics::NamedRing {
            triples: vec![three(0, 1, 2), three(0, 0, 3), three(0, 1, 3)],
            walls: vec![plane(1), plane(2), plane(3)],
            arc_ccw: vec![None; 3],
            concurrencies: vec![],
        };
        assert!(matches!(plane_ring(&collapsed), Err(RingFail::Collapsed)));
        // And the plane-only ring still comes back with both halves.
        let plain = combinatorics::NamedRing {
            triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
            walls: vec![plane(1), plane(2), plane(3)],
            arc_ccw: vec![None; 3],
            concurrencies: vec![],
        };
        let (ts, walls) = plane_ring(&plain).expect("a plane ring");
        assert_eq!(ts.len(), 3);
        assert_eq!(walls, vec![plane(1), plane(2), plane(3)]);
    }

    /// The graze rows of `edge_mask`'s algebra, one by one. The missing row was the
    /// same-side pair: a knife edge's two faces both leave `W` upward, the material
    /// between them pinches to measure zero at the plane, and flipping anything there
    /// hands the label propagation a contradiction it reports two layers later.
    #[test]
    fn edge_mask_graze_algebra() {
        let g = |above: bool| (SolidSide::A, SegKind::Graze { body_above: above });
        let flips = |m: &[(SolidSide, SegKind)]| edge_mask(m).unwrap();

        // One graze: a step — the face fills that arc, its bit flips.
        assert_eq!(flips(&[g(true)]), [true, false, false, false]);
        assert_eq!(flips(&[g(false)]), [false, true, false, false]);
        // An opposite pair: the solid's edge lies in W with material above on one
        // in-plane side and below on the other — both flip (the [T,T] picture).
        assert_eq!(flips(&[g(true), g(false)]), [true, true, false, false]);
        // ★ A same-side pair: knife edge (or its reflex complement) — nothing flips,
        // whichever side the pair is on.
        assert_eq!(flips(&[g(true), g(true)]), [false; 4]);
        assert_eq!(flips(&[g(false), g(false)]), [false; 4]);
        // The two solids are independent lanes.
        let gb = (SolidSide::B, SegKind::Graze { body_above: true });
        assert_eq!(flips(&[g(true), g(true), gb]), [false, false, true, false]);

        // More than two grazes from one solid is degenerate input, not a case.
        assert!(edge_mask(&[g(true), g(true), g(true)]).is_err());
        // Unchanged refusals: a graze coincident with a same-solid true crossing.
        let t = (SolidSide::A, SegKind::Transversal { mat: 1 });
        assert!(edge_mask(&[g(true), t]).is_err());
    }

    /// The engine entry with the evidence dropped — these tests assert geometry, and the report
    /// has its own tests. Shadows [`super::boolean`] so the call sites read as they always did.
    fn boolean(
        model: &mut Model,
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> Result<Vec<Handle<Solid>>, BoolError> {
        super::boolean(model, kind, a, b).map(|(solids, ..)| solids)
    }

    /// Point of a named vertex, for asserting geometry by hand.
    fn pt(n: NodeId, jd: &Judge<'_, WorkingPlane>) -> [f64; 3] {
        let t = three_plane_name(n).expect("a three-plane node");
        let planes = jd.planes;
        three_planes(
            &planes[t[0]].plane,
            &planes[t[1]].plane,
            &planes[t[2]].plane,
        )
        .unwrap()
        .as_array()
    }

    /// A seated face's whole boundary is emitted, with `body_above` matching the geometry.
    ///
    /// Two stacked cubes sharing the plane z=1: `a=[0,1]³`, `b=[0,1]²×[1,2]`. On the shared
    /// class, `a`'s top cap is seated with its body **below** (body_above=false) and `b`'s bottom
    /// cap is seated with its body **above** (body_above=true). Each emits its 4 boundary edges.
    #[test]
    fn a_seated_cap_emits_its_boundary_with_the_right_body_side() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);

        // The shared class z=1: the class both solids seat a cap on.
        let wc = (0..planes.len())
            .find(|&c| {
                let seats = |s: Handle<Solid>| {
                    solid_shell_handles(&m, s).into_iter().any(|sh| {
                        m.shells.get(sh).faces.iter().any(|fh| {
                            plane_ix[surf_ix[fh]].plane() == c
                                && face_on_z1(*fh, &surf_ix, &faces_tab)
                        })
                    })
                };
                seats(a) && seats(b)
            })
            .expect("a shared cap class");

        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        // Two seated caps × 4 edges each. (The side walls only *graze* z=1 — each cube's body is
        // entirely on one side — so they trace as `Graze` chords, not `Transversal`; the two far
        // caps are parallel to z=1 and contribute nothing — misses, not declines.)
        let seated: Vec<&Seg> = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Seated { .. }))
            .collect();
        assert_eq!(seated.len(), 8, "two seated caps, four edges each: {tr:?}");
        // Four side walls per cube each graze z=1 in one chord, body on the cube's side.
        let grazes: Vec<&Seg> = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Graze { .. }))
            .collect();
        assert_eq!(grazes.len(), 8, "eight side-wall graze chords: {tr:?}");
        assert_eq!(
            tr.segs
                .iter()
                .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
                .count(),
            0,
            "no wall crosses z=1: {tr:?}"
        );
        for s in &grazes {
            match (s.solid, s.kind) {
                (SolidSide::A, SegKind::Graze { body_above }) => assert!(!body_above, "a below"),
                (SolidSide::B, SegKind::Graze { body_above }) => assert!(body_above, "b above"),
                _ => unreachable!(),
            }
        }
        assert!(
            tr.declined.is_empty(),
            "axis-aligned cubes decline nothing: {tr:?}"
        );
        assert!(tr.touches.is_empty());
        for s in &seated {
            // On z=1, riding a side wall (not W itself), body on the geometrically correct side.
            for e in &s.end {
                assert!((pt(*e, &jd)[2] - 1.0).abs() < 1e-12, "endpoint on z=1");
            }
            assert_ne!(s.wall, wc, "a seated edge rides a side wall, not W");
            match (s.solid, s.kind) {
                (SolidSide::A, SegKind::Seated { body_above }) => assert!(!body_above, "a below"),
                (SolidSide::B, SegKind::Seated { body_above }) => assert!(body_above, "b above"),
                _ => unreachable!(),
            }
        }
    }

    /// Round a point for set membership by hand.
    fn rp(p: [f64; 3]) -> [i64; 3] {
        [
            (p[0] * 1e6).round() as i64,
            (p[1] * 1e6).round() as i64,
            (p[2] * 1e6).round() as i64,
        ]
    }

    /// The chord across a non-convex reflex plane — the cap chord — which the seated brick provably
    /// cannot produce (it declines this face). l_prism (profile (0,0),(2,0),(2,1),(1,1),(1,2),(0,2),
    /// extruded z∈[0,1]) cut at y=1: the z=0 cap covers x∈[0,2], as a **transversal** stretch
    /// x∈[0,1] where `y=1` runs through the cap's interior plus a **graze** x∈[1,2] where the cap's
    /// own boundary edge rides the line and the cap lies on one side of it.
    #[test]
    fn the_cap_chord_stops_where_the_on_line_edge_begins() {
        let profile = Profile2d::polygon(vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([2.0, 0.0]),
            Point2::from_array([2.0, 1.0]),
            Point2::from_array([1.0, 1.0]),
            Point2::from_array([1.0, 2.0]),
            Point2::from_array([0.0, 2.0]),
        ])
        .unwrap();
        let mut m = replay(&[Operation::Extrude {
            frame: SketchFrame::world(&Model::new(), Axis::Z),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        // A far cube so plane_index_setup has two solids; it never touches the y=1 class.
        m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        let a2 = m.live_solids[0];
        let b2 = m.live_solids[1];
        let PlaneSetup {
            cyls,
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a2, b2).unwrap();
        let jd = Judge::new(&planes, standard, &notes);

        // The y=1 class: a plane through all-y=1 points that a's reflex face sits on.
        let wc = (0..planes.len())
            .find(|&c| {
                planes[c]
                    .tri
                    .iter()
                    .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
            })
            .expect("a y=1 class");

        let mut out = Trace::default();
        trace_one_of(
            &m,
            a2,
            SolidSide::A,
            wc,
            &jd,
            &cyls,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &plane_ix,
            Default::default(),
            &mut out,
        );

        // Collect transversal segments as rounded endpoint-pairs.
        let chords: Vec<([i64; 3], [i64; 3])> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .collect();

        // The cap chords: the z=0 (and z=1) cap ∩ (y=1) covers x∈[0,2]. A cap is a horizontal
        // segment (both endpoints share z), which distinguishes it from the reflex wall x=1's own
        // vertical chord (1,1,0)-(1,1,1), a legitimate different face's trace.
        let horiz_z0: Vec<_> = chords
            .iter()
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        // ★ Falsifiable core: the on-line edge does not merge into the chord that runs through the
        // cap's interior. x∈[0,1] is interior (the cap straddles y=1) and stays transversal;
        // x∈[1,2] is the cap's own boundary edge, so the cap is on one side there and it leaves as
        // a graze — see `run_body_above`. A single chord spanning [0,2] would be the old,
        // occupancy-blind bridging.
        assert_eq!(
            horiz_z0.len(),
            1,
            "one transversal chord, the interior stretch: {chords:?}"
        );
        assert_eq!(
            *horiz_z0[0],
            (rp([0.0, 1.0, 0.0]), rp([1.0, 1.0, 0.0])),
            "the interior chord is x∈[0,1]: {chords:?}"
        );
        // The seated brick provably cannot produce this: it declines the z=0 cap as `not-seated`
        // (its class is z=0, not y=1). So a transversal chord riding fp=z=0 is new capability.
        // mat is a per-face constant, never 0.
        for s in &out.segs {
            if let SegKind::Transversal { mat } = s.kind {
                assert!(mat == 1 || mat == -1, "mat is a clean side: {mat}");
            }
        }
    }

    /// The flanks-equal (tangential on-line edge) branch: u_prism's notch bottom on y=1. The run
    /// x∈[1,2] is flanked by material on **both** sides (the two prongs), so it does not toggle
    /// parity yet must still appear inside the trace.
    ///
    /// The cap ∩ y=1 covers x∈[0,3], but **not as one kind**: over the notch bottom the line is on
    /// the cap's *boundary* and the cap lies below it, while on either side the line runs through
    /// the cap's *interior* and the cap straddles. So the trace is `T[0,1] · G[1,2] · T[2,3]`, and
    /// the graze's side is **below** — note the run's flanks (2,2.3) and (1,2) are both *above*,
    /// which is why reading the side off a flank gets a notch backwards.
    #[test]
    fn a_tangential_on_line_edge_spans_as_transversal_graze_transversal() {
        let profile = Profile2d::polygon(vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([3.0, 0.0]),
            Point2::from_array([3.0, 2.3]),
            Point2::from_array([2.0, 2.3]),
            Point2::from_array([2.0, 1.0]),
            Point2::from_array([1.0, 1.0]),
            Point2::from_array([1.0, 2.0]),
            Point2::from_array([0.0, 2.0]),
        ])
        .unwrap();
        let mut m = replay(&[Operation::Extrude {
            frame: SketchFrame::world(&Model::new(), Axis::Z),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        let a2 = m.live_solids[0];
        let b2 = m.live_solids[1];
        let PlaneSetup {
            cyls,
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a2, b2).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = (0..planes.len())
            .find(|&c| {
                planes[c]
                    .tri
                    .iter()
                    .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
            })
            .expect("a y=1 class");
        let mut out = Trace::default();
        trace_one_of(
            &m,
            a2,
            SolidSide::A,
            wc,
            &jd,
            &cyls,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &plane_ix,
            Default::default(),
            &mut out,
        );

        let mut horiz_z0: Vec<_> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        horiz_z0.sort();
        assert_eq!(
            horiz_z0,
            vec![
                (rp([0.0, 1.0, 0.0]), rp([1.0, 1.0, 0.0])),
                (rp([2.0, 1.0, 0.0]), rp([3.0, 1.0, 0.0])),
            ],
            "the interior stretches straddle; the notch bottom is not one of them"
        );

        // ★ The run itself is covered, as a one-sided graze — coverage is not lost, the kind
        //   differs. Its body is **below** (y < 1), the side the notch's material is on.
        let graze_z0: Vec<_> = out
            .segs
            .iter()
            .filter_map(|s| match s.kind {
                SegKind::Graze { body_above } => {
                    let (mut u, mut v) = (rp(pt(s.end[0], &jd)), rp(pt(s.end[1], &jd)));
                    if u > v {
                        std::mem::swap(&mut u, &mut v);
                    }
                    (u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
                        .then_some((u, v, body_above))
                }
                _ => None,
            })
            .collect();
        assert_eq!(
            graze_z0,
            vec![(rp([1.0, 1.0, 0.0]), rp([2.0, 1.0, 0.0]), false)],
            "the notch bottom x∈[1,2], body below"
        );
    }

    /// ★ Spike: a coordinate-free exact cyclic order of edges around an arrangement vertex is
    /// buildable — the "DNA question" the winding-based engine design flagged as a possible death
    /// condition for the per-plane arrangement. A degree-5 vertex with directions +x, +y, −x, −y and one **oblique**
    /// (the oblique is non-optional: without it every open half-plane bucket holds one element and
    /// transitivity never fires). `angular_order` returns the CCW order reading no coordinate; we
    /// check it against the CCW order computed *with* coordinates (atan2), which is the oracle.
    #[test]
    fn angular_order_around_a_vertex_is_coordinate_free_and_ccw() {
        // A pentagon prism giving y-, x-, and diagonal-normal side faces plus z caps.
        let profile = Profile2d::polygon(vec![
            Point2::from_array([0.0, 0.0]), // (0,0)-(3,0): y=0
            Point2::from_array([3.0, 0.0]), // (3,0)-(3,3): x=3
            Point2::from_array([3.0, 3.0]), // (3,3)-(2,3): y=3
            Point2::from_array([2.0, 3.0]), // (2,3)-(0,1): diagonal y=x+1
            Point2::from_array([0.0, 1.0]), // (0,1)-(0,0): x=0
        ])
        .unwrap();
        let m = replay(&[Operation::Extrude {
            frame: SketchFrame::world(&Model::new(), Axis::Z),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        let a = m.live_solids[0];
        let faces_tab = collect_planes(&m, a).unwrap();
        // One prism: no two faces are coplanar, so `plane_ix` is the identity and a face index and
        // its plane id coincide. Built through the real path anyway, so the test cannot drift.
        let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces_tab));
        let (planes, _plane_ix, _cyls) = crate::planes::dense_planes(&faces_tab, &canon);
        assert_eq!(planes.len(), faces_tab.len(), "no coplanar pair in a prism");

        // Find a face by its outward normal direction (z cap, y-wall, x-wall, diagonal wall).
        let axis = |i: usize| {
            let n = planes[i].plane.normal();
            [n[0], n[1], n[2]]
        };
        let find = |f: &dyn Fn([f64; 3]) -> bool| (0..planes.len()).find(|&i| f(axis(i)));
        let w = find(&|n| n[0].abs() < 1e-9 && n[1].abs() < 1e-9).expect("z cap"); // z-normal
        let fpy = find(&|n| n[0].abs() < 1e-9 && n[2].abs() < 1e-9).expect("y wall"); // y-normal
        let fpx = find(&|n| n[1].abs() < 1e-9 && n[2].abs() < 1e-9).expect("x wall"); // x-normal
        let fpd = find(&|n| {
            n[2].abs() < 1e-9 && n[0].abs() > 1e-6 && (n[0].abs() - n[1].abs()).abs() < 1e-6
        })
        .expect("diagonal wall");

        // The direction an edge (fp, s) actually runs, WITH coordinates — the oracle only.
        let n_w = planes[w].plane.normal();
        let dir = |fp: usize, s: i8| {
            let d = n_w.cross(planes[fp].plane.normal());
            [d[0] * s as f64, d[1] * s as f64, d[2] * s as f64]
        };
        // Five edges: both directions on the y-wall and x-wall lines, one on the diagonal.
        //
        // ★ **The pairs stay, and the oracle below reads *them*, not the `EdgeDir`s.** A direction's
        // representation is private to `combinatorics` on purpose (see `EdgeDir`); a test that needs
        // coordinates is an oracle, and an oracle reads its own inputs — reading them back out of
        // the value under test would derive the oracle from the answer.
        let raw = [(fpy, 1i8), (fpy, -1), (fpx, 1), (fpx, -1), (fpd, 1)];
        let edges: Vec<combinatorics::EdgeDir> = raw
            .iter()
            .map(|&(c, s)| combinatorics::EdgeDir::new(c, s))
            .collect();

        let order = angular_order(&crate::planes::test_judge(&planes), w, &edges)
            .expect("five distinct directions leave one vertex — nothing to refuse here");
        assert_eq!(
            order.len(),
            edges.len(),
            "every edge placed exactly once: {order:?}"
        );
        assert_eq!(order[0], 0, "order starts at the reference edge");

        // Oracle: the angular order of the actual directions (atan2), starting from edge 0.
        let ang = |e: (usize, i8)| {
            let d = dir(e.0, e.1);
            d[1].atan2(d[0])
        };
        let a0 = ang(raw[0]);
        let mut want: Vec<usize> = (0..edges.len()).collect();
        want.sort_by(|&i, &j| {
            let (ci, cj) = (
                (ang(raw[i]) - a0).rem_euclid(std::f64::consts::TAU),
                (ang(raw[j]) - a0).rem_euclid(std::f64::consts::TAU),
            );
            ci.partial_cmp(&cj).unwrap()
        });
        // The order is a consistent cyclic order — CW or CCW depending on orient_sign(w)'s
        // convention; both are correct rotations. So it matches the oracle, or the oracle with its
        // tail reversed (the same cycle traversed the other way, reference fixed).
        let mut rev_tail = order.clone();
        rev_tail[1..].reverse();
        assert!(
            order == want || rev_tail == want,
            "coordinate-free order is the atan2 cycle (either direction): got {order:?}, want {want:?}"
        );

        // ★ Transitivity fired: the oblique lands strictly between the +x and +y axis directions in
        // the cyclic order — a bucket-sort no-op could not place it. Checked as a cyclic adjacency
        // so it holds regardless of traversal direction.
        let cyc = if order == want { &order } else { &rev_tail };
        let pos = |e: usize| cyc.iter().position(|&x| x == e).unwrap();
        // edges: 0=(y,+) 1=(y,-) 2=(x,+) 3=(x,-) 4=(diag,+). The oblique (4) is 45°, between the
        // two axis directions flanking it. Identify its neighbours are axis edges, not each other.
        let obl = pos(4);
        let lo = cyc[(obl + cyc.len() - 1) % cyc.len()];
        let hi = cyc[(obl + 1) % cyc.len()];
        assert!(
            [0, 1, 2, 3].contains(&lo) && [0, 1, 2, 3].contains(&hi),
            "oblique is flanked by axis edges (transitivity placed it): {cyc:?}"
        );
    }

    /// On a cap plane, the cap rim and each side wall's top edge are one geometric edge traced
    /// twice (seated + transversal). Merge collapses coincident edges by (wall, endpoint set) while
    /// preserving contributions. Two footprint cases pin the count by hand.
    #[test]
    fn coincident_cap_edges_merge_by_footprint() {
        // (a) Same footprint: stacked cubes share the z=1 plane AND the same [0,1]² rim. Each of
        // the 4 rim edges is produced 4× (a-seated, a-wall, b-seated, b-wall) → 4 merged edges.
        {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([1.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.0, 0.0, 1.0]),
                Point3::from_array([1.0, 1.0, 2.0]),
            );
            m.rebuild_adjacency();
            let PlaneSetup {
                planes: faces_tab,
                geom: planes,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                standard,
                notes,
                ..
            } = plane_index_setup(&m, a, b).unwrap();
            let jd = Judge::new(&planes, standard, &notes);
            let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
            assert_eq!(merged.len(), 4, "one merged edge per rim edge: {merged:#?}");
            for e in &merged {
                assert_eq!(
                    e.merged.len(),
                    4,
                    "a·b × seated·transversal: {:?}",
                    e.merged
                );
            }
        }
        // (b) Different footprint: a cross. a's cap [0,3]×[1,2] and b's cap [1,2]×[0,3] are
        // different rims, so no a↔b coincidence: a's 4 pairs → 4, b's 4 pairs → 4 = 8.
        {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0, 1.0, 0.0]),
                Point3::from_array([3.0, 2.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([1.0, 0.0, 0.0]),
                Point3::from_array([2.0, 3.0, 1.0]),
            );
            m.rebuild_adjacency();
            let PlaneSetup {
                planes: faces_tab,
                geom: planes,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                standard,
                notes,
                ..
            } = plane_index_setup(&m, a, b).unwrap();
            let jd = Judge::new(&planes, standard, &notes);
            let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
            assert_eq!(
                merged.len(),
                8,
                "4 per box, no a↔b coincidence: {merged:#?}"
            );
            for e in &merged {
                assert_eq!(
                    e.merged.len(),
                    2,
                    "one solid × seated·transversal: {:?}",
                    e.merged
                );
            }
        }
    }

    /// Split at crossings: a cross fixture makes four strict interior crossings; the chords cut
    /// into the right pieces, ending at the right triples, and the crossing vertex is a clean
    /// degree-4 the DCEL can order.
    #[test]
    fn split_cuts_the_cross_into_the_right_pieces() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        assert_eq!(merged.len(), 8, "8 merged edges before split");

        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
        // 4 chords crossed twice → 3 pieces each = 12; 4 outer walls uncrossed = 4; total 16.
        assert_eq!(split.len(), 16, "16 sub-segments: {}", split.len());

        // a's y=1 chord splits into exactly 3, ending at x=0,1,2,3 — checked by coordinate.
        let seg_pts = |s: &MergedSeg| {
            let mut u = pt(s.end[0], &jd);
            let mut v = pt(s.end[1], &jd);
            if u[0] > v[0] {
                std::mem::swap(&mut u, &mut v);
            }
            (
                [(u[0] * 1e6).round() as i64, (u[1] * 1e6).round() as i64],
                [(v[0] * 1e6).round() as i64, (v[1] * 1e6).round() as i64],
            )
        };
        // a's y=1 chord pieces: horizontal (both y=1), z=1, spanning consecutive x's.
        let y1: Vec<_> = split
            .iter()
            .filter(|s| {
                let (u, v) = seg_pts(s);
                u[1] == 1_000_000 && v[1] == 1_000_000 // both endpoints at y=1
                    && s.merged.iter().any(|(sd, _)| *sd == SolidSide::A)
            })
            .map(&seg_pts)
            .collect();
        assert_eq!(y1.len(), 3, "a's y=1 chord → 3 pieces: {y1:?}");
        let mut xs: Vec<[i64; 2]> = y1.iter().map(|(u, v)| [u[0], v[0]]).collect();
        xs.sort();
        assert_eq!(
            xs,
            vec![
                [0, 1_000_000],
                [1_000_000, 2_000_000],
                [2_000_000, 3_000_000]
            ],
            "pieces are [0,1][1,2][2,3], adjacent and endpoint-exact"
        );

        // The crossing vertex (1,1,1) is a clean degree-4 the DCEL can order: 4 incident pieces,
        // and their (wall, far-R) pairs — i.e. (fp, direction) — are all distinct (no coincidence
        // that angular_order's zero bucket would collapse). This is what the merge earned.
        let is_v = |n: NodeId| {
            let p = pt(n, &jd);
            (p[0] - 1.0).abs() < 1e-9 && (p[1] - 1.0).abs() < 1e-9 && (p[2] - 1.0).abs() < 1e-9
        };
        let incident: Vec<&MergedSeg> = split
            .iter()
            .filter(|s| is_v(s.end[0]) || is_v(s.end[1]))
            .collect();
        assert_eq!(incident.len(), 4, "degree-4 at (1,1,1): {}", incident.len());
        let mut keys = std::collections::HashSet::new();
        for s in &incident {
            // The far end as a handle on this edge's own line — carried by the segment now.
            let far_h = if is_v(s.end[0]) {
                s.end_h[1]
            } else {
                s.end_h[0]
            };
            keys.insert((s.wall, far_h));
        }
        assert_eq!(keys.len(), 4, "no coincident (fp,direction) — DCEL-ready");
    }

    /// Same-wall partial overlap (E5) is **resolved** by the per-wall overlay: a=[0,2], b=[1,3]
    /// share y=1, their chords overlap on x∈[1,2]. The overlay splits the y=1 wall into three
    /// non-overlapping pieces `[0,1] [1,2] [2,3]`, and the shared middle `[1,2]` carries a
    /// contribution from **both** solids (which the label brick then reads per solid), while the
    /// flanks are single-solid.
    #[test]
    fn partial_overlap_is_resolved() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();

        // The shared y=1 wall (a face at y=1).
        let y1 = planes
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12))
            .expect("a y=1 face");
        // x-extent of a y=1 sub-segment, plus which solids contribute.
        let piece = |s: &MergedSeg| -> ([i64; 2], bool, bool) {
            let mut u = pt(s.end[0], &jd)[0];
            let mut v = pt(s.end[1], &jd)[0];
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            let has = |sd: SolidSide| s.merged.iter().any(|(x, _)| *x == sd);
            (
                [(u * 1e6).round() as i64, (v * 1e6).round() as i64],
                has(SolidSide::A),
                has(SolidSide::B),
            )
        };
        let mut pieces: Vec<([i64; 2], bool, bool)> =
            split.iter().filter(|s| s.wall == y1).map(piece).collect();
        pieces.sort();
        assert_eq!(
            pieces,
            vec![
                ([0, 1_000_000], true, false),         // [0,1] A only
                ([1_000_000, 2_000_000], true, true),  // [1,2] BOTH solids
                ([2_000_000, 3_000_000], false, true), // [2,3] B only
            ],
            "y=1 overlap resolved into three pieces, middle carries both solids: {pieces:?}"
        );
    }

    /// DCEL face-walk on the cross: 6 cells (1 outer + 5 bounded), exactly one winding -1, and the
    /// center square is the unique cell with 2 A-edges + 2 B-edges. All checkable without labels.
    #[test]
    fn the_cross_arrangement_has_six_cells() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();

        let edges =
            ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();

        // ★★★ **An independent reading of every winding: the shoelace sign.**
        //
        // `loop_winding` reads the turn at one extreme node, which is correct only if the
        // comparison it picks that node with is an order. When it was not, exactly one ring of
        // five came back with the wrong sign and stayed *confident* about it — the engine only
        // noticed two layers later, as "no outer contour", and named the symptom. Realizing the
        // nodes and summing cross products says the same thing directly, so a disagreement points
        // at the ring rather than at the count.
        {
            let pl = &jd.planes[wc].plane;
            let n = pl.normal();
            let ax = if n[0].abs() < 0.9 {
                Vector3::from_array([1.0, 0.0, 0.0])
            } else {
                Vector3::from_array([0.0, 1.0, 0.0])
            };
            let u = n.cross(ax).normalize().expect("in-plane axis");
            let v = n.cross(u);
            let o = Point3::from_array([0.0; 3]);
            for c in &cells {
                let pts: Vec<[f64; 2]> = c
                    .half_edges
                    .iter()
                    .filter_map(|&h| {
                        let t =
                            three_plane_name(split[h / 2].end[h % 2]).expect("a three-plane node");
                        three_planes(
                            &jd.planes[t[0]].plane,
                            &jd.planes[t[1]].plane,
                            &jd.planes[t[2]].plane,
                        )
                        .map(|q| [u.dot(q - o), v.dot(q - o)])
                    })
                    .collect();
                assert_eq!(pts.len(), c.half_edges.len(), "every node realizes");
                let area: f64 = (0..pts.len())
                    .map(|i| {
                        let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                        a[0] * b[1] - b[0] * a[1]
                    })
                    .sum();
                assert!(
                    area.abs() > 1e-9,
                    "a degenerate ring gives the oracle nothing to say"
                );
                assert_eq!(
                    if area > 0.0 { 1i8 } else { -1 },
                    c.winding,
                    "shoelace {area:+.6e} disagrees with the reported winding {}",
                    c.winding
                );
            }
        }
        assert_eq!(cells.len(), 6, "1 outer + 5 bounded: {}", cells.len());
        assert_eq!(
            cells.iter().filter(|c| c.winding == -1).count(),
            1,
            "exactly one outer (winding -1)"
        );
        assert_eq!(
            cells.iter().filter(|c| c.winding == 1).count(),
            5,
            "five bounded (winding +1)"
        );
        // Every half-edge belongs to exactly one cell.
        assert_eq!(face_of.len(), 2 * split.len(), "2E half-edges all placed");

        // The center square [1,2]²: the unique bounded cell with 2 A-edges and 2 B-edges. (a and b
        // ride disjoint walls, so there are no both-solid edges — provenance, not coincidence.)
        let solids_of = |he: usize| -> Vec<SolidSide> {
            split[he / 2].merged.iter().map(|(sd, _)| *sd).collect()
        };
        let mut center = 0;
        for c in cells.iter().filter(|c| c.winding == 1) {
            let (mut na, mut nb) = (0, 0);
            for &he in &c.half_edges {
                let s = solids_of(he);
                if s.contains(&SolidSide::A) {
                    na += 1;
                }
                if s.contains(&SolidSide::B) {
                    nb += 1;
                }
            }
            if na == 2 && nb == 2 {
                center += 1;
            }
        }
        assert_eq!(
            center, 1,
            "exactly one center square (2 A-edges + 2 B-edges)"
        );
    }

    /// Cell labels propagate from the void, and match footprint containment (an independent
    /// oracle). Cross: W=z=1, a=[0,3]×[1,2], b=[1,2]×[0,3]. Both bodies below z=1 → aboves all F;
    /// each cell's X_below = "cell centroid inside X's footprint".
    #[test]
    fn cell_labels_propagate_and_match_footprints() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
        let edges =
            ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
        let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

        for (i, c) in cells.iter().enumerate() {
            if c.winding == -1 {
                assert_eq!(labels[i], [false; 4], "unbounded is void");
                continue;
            }
            // Centroid of the cell (convex here) via the average of its vertex points.
            let verts: Vec<[f64; 3]> = {
                let mut vs: Vec<NodeId> = c
                    .half_edges
                    .iter()
                    .map(|&h| split[h / 2].end[h % 2])
                    .collect();
                vs.dedup();
                vs.iter().map(|&t| pt(t, &jd)).collect()
            };
            let cx = verts.iter().map(|p| p[0]).sum::<f64>() / verts.len() as f64;
            let cy = verts.iter().map(|p| p[1]).sum::<f64>() / verts.len() as f64;
            let in_a = (0.0..=3.0).contains(&cx) && (1.0..=2.0).contains(&cy);
            let in_b = (1.0..=2.0).contains(&cx) && (0.0..=3.0).contains(&cy);
            // aboves F (bodies below the cap); belows = footprint containment.
            assert_eq!(
                labels[i],
                [false, in_a, false, in_b],
                "cell centroid ({cx},{cy}) label vs footprint"
            );
        }
    }

    /// A straddle: a=[0,1]²×[0,2] passes THROUGH z=1 (transversal, no face there), b=[0,1]²×[1,2]
    /// caps at z=1 from above. Same [0,1]² footprint → one square, 2 cells. The inner cell exercises
    /// transversal-flips-both (a's above reaches T) and seated-wins on b in the same edge.
    #[test]
    fn a_straddle_drives_the_above_bit_and_seated_wins() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        // a passes through z=1 (no seated face); find the class b caps at z=1.
        let wc = (0..planes.len())
            .find(|&c| {
                solid_shell_handles(&m, b).into_iter().any(|sh| {
                    m.shells.get(sh).faces.iter().any(|fh| {
                        plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, &surf_ix, &faces_tab)
                    })
                })
            })
            .expect("b's z=1 cap class");
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
        let edges =
            ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
        let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

        assert_eq!(cells.len(), 2, "one square: inner + unbounded");
        let inner = cells.iter().position(|c| c.winding == 1).unwrap();
        // "above/below" is relative to +n_w (the class plane's stored normal), whose sign is a
        // convention the boolean's keep-rule treats symmetrically. Here n_w points -z, so +n_w
        // ("above") is physically below. A straddles z=1 → material on BOTH sides (T,T); b caps
        // from z>1 → material only on -n_w ("below") → (B_above=F, B_below=T). The point of the
        // fixture stands: the transversal drives A's above-bit to T (impossible on a pure cap),
        // and B is seated-wins (a single bit) on the very same rim edge.
        assert_eq!(
            labels[inner],
            [true, true, false, true],
            "A straddles both sides; B is one-sided (seated-wins); transversal reached the above-bit"
        );
    }

    /// Boolean keep-decision + result-face emission on the cross. Fuse keeps 5 cells (the plus
    /// cap), Cut a−b keeps the 2 a-arms, Common keeps the 1 center — hand-verified — and the
    /// emitted loops share interior edges in opposite directions with a consistent winding and a
    /// flip that points the normal out of the kept solid.
    #[test]
    fn boolean_keep_and_result_faces_on_the_cross() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 1.0, 0.0]),
            Point3::from_array([3.0, 2.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, NO_CYLS, wc, &merged, &mut Aliases::default()).unwrap();
        let edges =
            ClassEdges::of(&jd, NO_CYLS, wc, &split, &[], &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, NO_CYLS, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, NO_CYLS, wc, &cells, &edges).unwrap();
        let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();

        // Centroid of a face's ring (convex cells here).
        let centroid = |f: &LocalFace| -> [f64; 2] {
            let ps: Vec<[f64; 3]> = f.outer.expect_ring().iter().map(|&n| pt(n, &jd)).collect();
            [
                ps.iter().map(|p| p[0]).sum::<f64>() / ps.len() as f64,
                ps.iter().map(|p| p[1]).sum::<f64>() / ps.len() as f64,
            ]
        };
        let near = |c: [f64; 2], x: f64, y: f64| (c[0] - x).abs() < 1e-9 && (c[1] - y).abs() < 1e-9;

        // Fuse: all 5 bounded cells (the plus cap).
        let (fuse, _, _) = emit_faces(
            BoolKind::Fuse,
            &labels,
            &cells,
            &edges,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(fuse.len(), 5, "Fuse keeps the whole plus cap");
        // Cut a−b: exactly the two a-arms.
        let (cut, _, _) = emit_faces(
            BoolKind::Cut,
            &labels,
            &cells,
            &edges,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(cut.len(), 2, "Cut keeps the a-arms");
        let cut_c: Vec<[f64; 2]> = cut.iter().map(centroid).collect();
        assert!(
            cut_c.iter().any(|c| near(*c, 0.5, 1.5)) && cut_c.iter().any(|c| near(*c, 2.5, 1.5)),
            "the two survivors are the a-arms at (0.5,1.5),(2.5,1.5): {cut_c:?}"
        );
        // Common: exactly the center square.
        let (common, _, _) = emit_faces(
            BoolKind::Common,
            &labels,
            &cells,
            &edges,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(common.len(), 1, "Common keeps the center");
        assert!(near(centroid(&common[0]), 1.5, 1.5), "center at (1.5,1.5)");

        // Each emitted loop winds +1 (CCW about n_out(wc)) and has ≥3 distinct nodes.
        for f in fuse.iter().chain(&cut).chain(&common) {
            let ring: Vec<[usize; 3]> = f
                .outer
                .expect_ring()
                .iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect();
            assert!(ring.len() >= 3);
            assert_eq!(
                combinatorics::loop_winding(
                    &jd,
                    NO_CYLS,
                    wc,
                    &combinatorics::ring_from_names(wc, &ring).unwrap()
                )
                .unwrap(),
                1,
                "CCW about n_out"
            );
        }

        // flip oracle (coordinate, test-only): the result normal points away from the kept
        // chamber. Fuse keeps below (bodies below the cap), so n_result·n_w > 0.
        let n_w = planes[wc].plane.normal();
        let os = planes[wc].frame_sign as f64;
        for f in &fuse {
            let n_result = os * if f.flip { -1.0 } else { 1.0 };
            let dot = n_result * n_w.dot(n_w); // n_result·n_w, |n_w|²>0
            assert!(
                dot > 0.0,
                "Fuse (keep below) normal points +n_w side: flip={}",
                f.flip
            );
        }

        // Edge-parity: count how many times each undirected edge is emitted, and the net
        // direction. An edge shared by two survivors (count 2) must net to 0 (opposite directions
        // — assemble's "each edge twice"); a plus-cap boundary edge (count 1) is excluded (it
        // closes against a wall face only at assembly). At least one interior edge must exist.
        let mut edges: HashMap<([usize; 3], [usize; 3]), (i32, i32)> = HashMap::new();
        for f in &fuse {
            let ns: Vec<[usize; 3]> = f
                .outer
                .expect_ring()
                .iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect();
            for w in ns
                .windows(2)
                .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
            {
                let (key, sign) = if w[0] < w[1] {
                    ((w[0], w[1]), 1)
                } else {
                    ((w[1], w[0]), -1)
                };
                let e = edges.entry(key).or_insert((0, 0));
                e.0 += 1;
                e.1 += sign;
            }
        }
        let mut interior = 0;
        for (e, (count, net)) in &edges {
            if *count == 2 {
                interior += 1;
                assert_eq!(
                    *net, 0,
                    "interior edge {e:?} emitted in opposite directions"
                );
            }
        }
        assert!(
            interior > 0,
            "the plus cap has interior edges between arms and center"
        );
    }

    /// Intermediate lock (before assembly, A/B isolation): the driver emits exactly 10 faces for
    /// stacked Fuse (z=0 cap + z=2 cap + 4 walls × 2 z-split), no z=1 face, and every undirected
    /// edge appears exactly twice across all faces — so a later NON_MANIFOLD reject is a weld gap,
    /// not an emit gap.
    #[test]
    fn stacked_fuse_emits_ten_closed_faces() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            plane_ix,
            class_owner,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let trace_in = combinatorics::trace_input(
            &m,
            [(a, &inc_a), (b, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &cyls,
            Default::default(),
        );
        let (faces, _, _) = trace_result_faces(
            &m,
            BoolKind::Fuse,
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
        )
        .unwrap();
        assert_eq!(faces.len(), 10, "z=0 + z=2 + 4 walls×2");

        // Every undirected edge (sorted triple pair) is used exactly twice — a closed shell.
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            let ns: Vec<[usize; 3]> = f
                .outer
                .expect_ring()
                .iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect();
            for w in ns
                .windows(2)
                .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
            {
                let key = if w[0] < w[1] {
                    (w[0], w[1])
                } else {
                    (w[1], w[0])
                };
                *count.entry(key).or_insert(0) += 1;
            }
        }
        assert!(
            count.values().all(|&c| c == 2),
            "every edge used exactly twice (closed shell): {:?}",
            count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
        );
    }

    /// First end-to-end: drive every plane class, assemble, check volume + manifold. Stacked cubes
    /// Fuse = a 1×1×2 box, volume 2.0; the z=1 interface face vanishes (kept both sides).
    #[test]
    fn end_to_end_stacked_fuse_is_a_tall_box() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        assert_eq!(solids.len(), 1, "one connected solid");
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "manifold: {vs:?}");
        let vol = nacre_props::mass_props(&m, solids[0]).unwrap().volume;
        assert!(
            (vol - 2.0).abs() < 1e-9,
            "stacked Fuse volume 2.0, got {vol}"
        );
    }

    /// Non-degenerate all-three: overlapping cubes a=[0,1]³, b=[0.5,1.5]³.
    /// Fuse = 2 − 0.5³ = 1.875, Common = 0.5³ = 0.125, Cut = 1 − 0.125 = 0.875.
    #[test]
    fn end_to_end_overlapping_cubes_all_three() {
        let vol_of = |kind: BoolKind| -> f64 {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([1.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.5, 0.5, 0.5]),
                Point3::from_array([1.5, 1.5, 1.5]),
            );
            m.rebuild_adjacency();
            let solids = boolean(&mut m, kind, a, b).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
            solids
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                .sum()
        };
        assert!((vol_of(BoolKind::Fuse) - 1.875).abs() < 1e-9, "Fuse 1.875");
        assert!((vol_of(BoolKind::Cut) - 0.875).abs() < 1e-9, "Cut 0.875");
        assert!(
            (vol_of(BoolKind::Common) - 0.125).abs() < 1e-9,
            "Common 0.125"
        );
    }

    /// cube[0,3]³ and a square prism [1,2]²×[−1,4] piercing it in z; run the boolean, validate,
    /// and sum (volume, area) over the result solids.
    fn tunnel_vol_area(kind: BoolKind) -> (f64, f64) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        m.rebuild_adjacency();
        let solids = boolean(&mut m, kind, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
        solids.iter().fold((0.0, 0.0), |acc, &s| {
            let p = nacre_props::mass_props(&m, s).unwrap();
            (acc.0 + p.volume, acc.1 + p.area)
        })
    }

    /// The first result with a hole. On the z=0/z=3 caps the arrangement is two disjoint loops
    /// (cube rim `O`, tunnel mouth `I`); the mouth is a hole of the annular cap, which cannot fall
    /// out as an absent cell (tiling the annulus needs a rim node the adjacent wall lacks) — it
    /// must emit as an inner ring. Cut volume 24 (27−3); AREA 64 (not 54) proves a bore, not a
    /// filled cap — volume alone cannot tell a through-hole from a blind dent.
    #[test]
    fn tunnel_cut_is_bored_not_dented() {
        let (v, area) = tunnel_vol_area(BoolKind::Cut);
        assert!((v - 24.0).abs() < 1e-9, "Cut vol 24 (27−1·1·3), got {v}");
        assert!(
            (area - 64.0).abs() < 1e-9,
            "Cut area 64 (36 walls + 8+8 punched caps + 12 tunnel walls; a filled mouth is 54): {area}"
        );
    }

    /// Fuse: the bar protrudes z∈[−1,0] and z∈[3,4], so each cap is still annular (the peg passes
    /// through the mouth, kept on both sides) — volume 29 (27+5−3), area 62.
    #[test]
    fn tunnel_fuse_has_annular_caps() {
        let (v, area) = tunnel_vol_area(BoolKind::Fuse);
        assert!((v - 29.0).abs() < 1e-9, "Fuse vol 29 (27+5−3), got {v}");
        assert!(
            (area - 62.0).abs() < 1e-9,
            "Fuse area 62 (cube 52 + two 1×1×1 pegs at 5 each): {area}"
        );
    }

    /// Common: `nest_cells` still finds the same annulus grouping, but `keep` leaves only the mouth
    /// (inside `I`); the annular host is skipped, its inner ring never leaks, and the mouth emits a
    /// plain `[1,2]²` cap. Net a plain 1×1×3 box — vol 3, area 14. A stronger control than a
    /// hole-free fixture: it exercises the nesting path and proves it does not leak a spurious hole.
    #[test]
    fn tunnel_common_is_the_bar_box_control() {
        let (v, area) = tunnel_vol_area(BoolKind::Common);
        assert!(
            (v - 3.0).abs() < 1e-9,
            "Common vol 3 ([1,2]²×[0,3]), got {v}"
        );
        assert!(
            (area - 14.0).abs() < 1e-9,
            "Common area 14 (2·1 caps + 4·3 walls): {area}"
        );
    }

    /// Intermediate lock (before assembly, A/B isolation): the Cut tunnel emits exactly 10 faces
    /// (2 annular caps + 4 cube walls + 4 tunnel walls), exactly 2 carry a non-empty inner ring
    /// (the caps), and every undirected edge across all rings (outer + inner) is used exactly twice
    /// — a closed shell, so a later NON_MANIFOLD is a weld gap, not an emit gap.
    #[test]
    fn tunnel_cut_emits_ten_faces_two_annular() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            plane_ix,
            class_owner,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let trace_in = combinatorics::trace_input(
            &m,
            [(a, &inc_a), (b, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &cyls,
            Default::default(),
        );
        let (faces, _, _) = trace_result_faces(
            &m,
            BoolKind::Cut,
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
        )
        .unwrap();
        assert_eq!(
            faces.len(),
            10,
            "2 annular caps + 4 cube walls + 4 tunnel walls"
        );
        assert_eq!(
            faces.iter().filter(|f| !f.inner.is_empty()).count(),
            2,
            "exactly the two annular caps carry a hole"
        );

        // Every undirected edge across outer + inner rings is used exactly twice (closed shell).
        let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in f.poly_rings() {
                let ns = triples(ring);
                for w in ns
                    .windows(2)
                    .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
                {
                    let key = if w[0] < w[1] {
                        (w[0], w[1])
                    } else {
                        (w[1], w[0])
                    };
                    *count.entry(key).or_insert(0) += 1;
                }
            }
        }
        assert!(
            count.values().all(|&c| c == 2),
            "every edge used exactly twice (closed shell): {:?}",
            count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
        );
    }

    /// Multi-hole lock: a slab fused with a U (its two prongs pierce the slab) emits exactly one
    /// face with **two** inner rings — the y=2 slab annulus, holed by both prongs — and every
    /// undirected edge across all rings is used exactly twice (closed shell). This exercises
    /// `nest_cells`' multi-hole support (the dropped `HOLE_MULTI` reject).
    #[test]
    fn u_slab_fuse_emits_a_two_hole_face() {
        let u_profile = Profile2d::polygon(vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([3.0, 0.0]),
            Point2::from_array([3.0, 2.3]),
            Point2::from_array([2.0, 2.3]),
            Point2::from_array([2.0, 1.0]),
            Point2::from_array([1.0, 1.0]),
            Point2::from_array([1.0, 2.0]),
            Point2::from_array([0.0, 2.0]),
        ])
        .unwrap();
        let mut m = replay(&[Operation::Extrude {
            frame: SketchFrame::world(&Model::new(), Axis::Z),
            profile: u_profile,
            dist: 1.0,
        }])
        .unwrap();
        let u = m.live_solids[0];
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.5, -0.5]),
            Point3::from_array([3.5, 2.5, 1.5]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            plane_ix,
            class_owner,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, u, slab).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let trace_in = combinatorics::trace_input(
            &m,
            [(u, &inc_a), (slab, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &cyls,
            Default::default(),
        );
        let (faces, _, _) = trace_result_faces(
            &m,
            BoolKind::Fuse,
            u,
            slab,
            &jd,
            &faces_tab,
            &plane_ix,
            &cyls,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &trace_in,
        )
        .unwrap();

        // Exactly one face carries two inner rings: the y=2 slab annulus, holed by both prongs.
        assert_eq!(
            faces.iter().filter(|f| f.inner.len() == 2).count(),
            1,
            "the slab face on y=2 has two holes (the U's two prongs)"
        );

        // Every undirected edge across outer + inner rings is used exactly twice (closed shell).
        let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in f.poly_rings() {
                let ns = triples(ring);
                for w in ns
                    .windows(2)
                    .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
                {
                    let key = if w[0] < w[1] {
                        (w[0], w[1])
                    } else {
                        (w[1], w[0])
                    };
                    *count.entry(key).or_insert(0) += 1;
                }
            }
        }
        assert!(
            count.values().all(|&c| c == 2),
            "every edge used exactly twice (closed shell): {:?}",
            count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
        );
    }

    // --- Rotation-generality: the trace engine's decisions are coordinate-free. Rigidly rotating
    // both operands by the same isometry must leave the result invariant. A non-90° angle flags the
    // solid, whose faces record the motion, routing every predicate to the exact frame3 backend. `rot30` (lib.rs)
    // is in a sibling test module and unreachable here, so the isometries are built inline.

    /// 30° about `axis` through (1,1,0) — non-90°, so the motion is recorded (exact frame3 path).
    fn rot_iso(axis: nacre_scalar::Axis) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation};
        Isometry::rotation(Rotation {
            axis,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        })
    }

    /// Rotate a solid by each axis in turn. A compound tilt needs `rebuild_adjacency` BETWEEN the
    /// transforms (matching the oracle's `rotated_boolean_matches_occt`), or the second reads a
    /// stale topology.
    fn tilt(m: &mut Model, mut s: Handle<Solid>, axes: &[nacre_scalar::Axis]) -> Handle<Solid> {
        for &ax in axes {
            s = transform(m, s, &rot_iso(ax)).unwrap();
            m.rebuild_adjacency();
        }
        s
    }

    /// The boolean's error and the class audit's `failed_at` are the **same** reject.
    ///
    /// They are two consumers of one `DeclineKind → RejectReason` mapping, and before
    /// `decline_to_reject` they were two copies of it. A copy that drifts makes the audit — the
    /// tool used to debug a reject — disagree with the reject being debugged, which is the worst
    /// possible time to be lying.
    ///
    /// The model is **chosen by measurement, not by taste**: a sweep over the fixture corpus found
    /// no boolean that declines inside the arrangement at all (every fixture builds), so the
    /// agreement had to be pinned on an input that still stops after the classes have run. If a
    /// later capability makes it build, the fix is to take whatever still declines — not to weaken
    /// the assertion.
    #[test]
    fn the_audit_does_not_invent_failures() {
        // The audit's scope is the per-class pipeline, and its duty is to run **the pipeline the
        // boolean runs** — with the alias fixpoint. Audited against an empty alias table it
        // reported `UnorderedEdges` for three classes of this input (names two classes discover
        // for each other were missing), failures the boolean never had: an instrument that
        // invents readings. The boolean's own reject here (`StraightAngle`) comes from the
        // assembly's vertex naming, after every class pipeline has run and outside the audit's
        // scope — so the audit's honest answer for this input is "no class failed".
        //
        // ★ **The fixture has moved twice, exactly as the note above prescribes**, and each move
        // is a capability the kernel gained. First it was the same pair at 30° with `Cut`,
        // rejecting `DegenerateWitness` — a reject from the component outwardness test, which the
        // nesting-parity label replaced. Then it was that pair at 60°, rejecting `StraightAngle`:
        // two unit cubes that only *touch*, which now come back as the two bodies they are
        // (measured: `Common` empty, fused volume 2.0, `validate` clean).
        //
        // So the fixture is now a pinch that **cannot** part: A and B meet only along the line
        // `x = 2, y = 2`, and a bridge overlapping both runs the material around the contact, so
        // cutting there leaves one piece. The closed-shell guard in `boolean::reconstruct` rejects
        // it — after every class pipeline has run, which is the property this test needs.
        //
        // (This test once asserted the audit reports the boolean's *class-level* reject, on a
        // fixture chosen as "some input that rejects" — a bar rotated through an L-shaped
        // target. The knife-edge fix turned that family, and every class-level-rejecting valid
        // input we could construct, into answers; the shared mapping the old test guarded,
        // `decline_to_reject`, is one function called by both consumers, so it cannot drift.)
        let build = || -> (Model, Handle<Solid>, Handle<Solid>) {
            let mut m = Model::new();
            let cub = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
                let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
                m.rebuild_adjacency();
                s
            };
            let a = cub(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
            let bridge = cub(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
            let b = cub(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
            let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
            m.rebuild_adjacency();
            (m, ab[0], b)
        };
        let (mut m, a, b) = build();
        let err = boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err();
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::NonManifoldResultEdge,
                    ..
                }
            ),
            "the fixture's premise: a reject from outside the class pipeline (got {err:?})"
        );
        let (m, a, b) = build();
        let audits = frame_audit(&m, BoolKind::Fuse, a, b).unwrap();
        let failed: Vec<RejectReason> = audits.iter().filter_map(|x| x.failed_at).collect();
        assert_eq!(
            failed,
            vec![],
            "every class runs clean under the boolean's own alias table — a failure invented \
             here is an artifact of running a different pipeline than the boolean runs"
        );
    }

    /// Rigidly rotating both operands (same single-Z tilt) leaves all three booleans' volumes
    /// invariant — the axis values (hand-anchored by `end_to_end_overlapping_cubes_all_three`)
    /// transfer to the rotated case. A coordinate-dependent decision that flipped under rotation
    /// would add/drop a cell and move the volume by O(0.1), far past 1e-9.
    #[test]
    fn rotated_overlapping_cubes_all_three() {
        use nacre_scalar::Axis;
        let vol_of = |kind: BoolKind| -> f64 {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([1.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.5, 0.5, 0.5]),
                Point3::from_array([1.5, 1.5, 1.5]),
            );
            let a = tilt(&mut m, a, &[Axis::Z]);
            let b = tilt(&mut m, b, &[Axis::Z]);
            let solids = boolean(&mut m, kind, a, b).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{kind:?} manifold: {vs:?}");
            solids
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                .sum()
        };
        assert!(
            (vol_of(BoolKind::Fuse) - 1.875).abs() < 1e-9,
            "rotated Fuse invariant 1.875"
        );
        assert!(
            (vol_of(BoolKind::Cut) - 0.875).abs() < 1e-9,
            "rotated Cut invariant 0.875"
        );
        assert!(
            (vol_of(BoolKind::Common) - 0.125).abs() < 1e-9,
            "rotated Common invariant 0.125"
        );
    }

    /// The through-tunnel Cut is invariant under every orientation: single Z, X, Y, and compound
    /// Z∘X. The axis-DEPENDENCE was the tell of the bug (Y worked; Z/X declined before the
    /// `Judge::orient3d` on-plane fix), so all four orientations returning 24 is the fix's direct
    /// regression lock. The cube's own z=0/z=3 caps exercise the seated path under rotation.
    #[test]
    fn rotated_tunnel_cut_all_orientations() {
        use nacre_scalar::Axis;
        let vol_of = |axes: &[Axis]| -> f64 {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0, 0.0, 0.0]),
                Point3::from_array([3.0, 3.0, 3.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([1.0, 1.0, -1.0]),
                Point3::from_array([2.0, 2.0, 4.0]),
            );
            let a = tilt(&mut m, a, axes);
            let b = tilt(&mut m, b, axes);
            let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "manifold: {vs:?}");
            solids
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                .sum()
        };
        for axes in [
            &[Axis::Z][..],
            &[Axis::X][..],
            &[Axis::Y][..],
            &[Axis::Z, Axis::X][..],
        ] {
            let v = vol_of(axes);
            assert!(
                (v - 24.0).abs() < 1e-9,
                "tunnel Cut vol 24 for {axes:?}, got {v}"
            );
        }
    }

    /// Compound-tilted (no face normal axis-aligned) tunnel Cut: AREA 64 is invariant (the bore
    /// discriminator volume cannot see), and the pre-assembly face set is combinatorially identical
    /// to the axis case (`tunnel_cut_emits_ten_faces_two_annular`) — 10 faces, 2 with an inner ring,
    /// every undirected edge twice — proving the arrangement itself survived rotation.
    #[test]
    fn rotated_tunnel_area_and_faces() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
        let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);

        // Combinatorial invariant (pre-assembly, A/B isolation).
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            plane_ix,
            class_owner,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let trace_in = combinatorics::trace_input(
            &m,
            [(a, &inc_a), (b, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &cyls,
            Default::default(),
        );
        let (faces, _, _) = trace_result_faces(
            &m,
            BoolKind::Cut,
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
        )
        .unwrap();
        assert_eq!(faces.len(), 10, "rotated arrangement keeps 10 faces");
        assert_eq!(
            faces.iter().filter(|f| !f.inner.is_empty()).count(),
            2,
            "the two annular caps survive rotation"
        );
        let triples = |ns: &[combinatorics::NodeId]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|&n| three_plane_name(n).expect("a three-plane node"))
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in f.poly_rings() {
                let ns = triples(ring);
                for w in ns
                    .windows(2)
                    .chain(std::iter::once(&[ns[ns.len() - 1], ns[0]][..]))
                {
                    let key = if w[0] < w[1] {
                        (w[0], w[1])
                    } else {
                        (w[1], w[0])
                    };
                    *count.entry(key).or_insert(0) += 1;
                }
            }
        }
        assert!(
            count.values().all(|&c| c == 2),
            "every edge used exactly twice: {:?}",
            count.iter().filter(|(_, c)| **c != 2).collect::<Vec<_>>()
        );

        // Area (assembled).
        let solids = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let area: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().area)
            .sum();
        assert!(
            (area - 64.0).abs() < 1e-9,
            "rotated Cut area 64 (bore survives rotation), got {area}"
        );
    }

    /// The sharp tripwire: a compound-tilted tunnel must decline nothing, exactly like the axis
    /// baseline `axis_aligned_cubes_decline_nothing`. This is the EXACT site the bug broke
    /// (`coincident-features` from a `wall == wc` degenerate), so it fails immediately if the
    /// `Judge::orient3d` on-plane fix is reverted.
    #[test]
    fn rotated_tunnel_declines_nothing() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let a = tilt(&mut m, a, &[Axis::Z, Axis::X]);
        let b = tilt(&mut m, b, &[Axis::Z, Axis::X]);
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        for wc in 0..planes.len() {
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            assert!(
                tr.declined.is_empty(),
                "rotated class {wc} declined: {tr:?}"
            );
        }
    }

    /// The z=1 class both solids seat a cap on.
    fn shared_cap_class(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        surf_ix: &HashMap<Handle<Face>, usize>,
        faces: &[FaceRow],
        plane_ix: &[ClassIx],
    ) -> usize {
        let n_planes = plane_ix
            .iter()
            .map(|c| c.plane())
            .max()
            .map_or(0, |m| m + 1);
        (0..n_planes)
            .find(|&c| {
                let seats = |s: Handle<Solid>| {
                    solid_shell_handles(m, s).into_iter().any(|sh| {
                        m.shells.get(sh).faces.iter().any(|fh| {
                            plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, surf_ix, faces)
                        })
                    })
                };
                seats(a) && seats(b)
            })
            .expect("a shared z=1 cap class")
    }

    /// ★★ **The production path, run past the C2 stopper on a cylinder-bearing model.**
    ///
    /// Three `.plane()` calls sat on the production road with cylinder rows flowing into them —
    /// the `work` class list, the decline witness, and the report's `class_of` — and every one
    /// of them would have aborted the kernel on the first drill the moment C4b-3 removes the
    /// stopper. None was reachable while the stopper stood, so neither the suite nor the census
    /// said a word; reading the ~30 call sites did not find them either. **Running the road is
    /// what finds them**, which is why this test exists before the bands do.
    ///
    /// What it asserts is deliberately weak on geometry and strong on survival: the tracer
    /// completes, every face it emits is still a plane face (bands arrive in C4b-2), and the
    /// drilled cap carries its circular hole.
    #[test]
    fn the_production_road_survives_a_cylinder_past_the_stopper() {
        let mut m = Model::new();
        let (a, b) = drilled(&mut m, -1.0, 4.0); // a through-hole: circles on both box caps
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let trace_in = combinatorics::trace_input(
            &m,
            [(a, &inc_a), (b, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &cyls,
            Default::default(),
        );
        let (faces, _, _) = trace_result_faces(
            &m,
            BoolKind::Cut,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            &cyls,
            n_a,
            &class_owner,
            // `Proved` on purpose: the reuse guard must be what turns the shortcut off, not the
            // caller. Without it `pass_through` would carry planes across and drop the lateral.
            crate::reuse::ClassReuse::Proved,
            &trace_in,
        )
        .expect("the drill population traces");
        assert!(
            faces.iter().all(|f| matches!(f.surf, ClassIx::Plane(_))),
            "the plane arrangement emits plane faces only; bands are C4b-2"
        );
        // The box: 4 walls + 2 caps, and each cap carries the drill's circular hole.
        assert_eq!(faces.len(), 6, "{faces:?}");
        let holed = faces
            .iter()
            .filter(|f| {
                f.inner
                    .iter()
                    .any(|b| matches!(b, crate::boolean::Bound::Circle { .. }))
            })
            .count();
        assert_eq!(holed, 2, "both caps are drilled: {faces:?}");
    }

    /// The rulings-road harness (M6-2 rulings ladder, first cell): the through-boss geometry
    /// the ladder will open — a plate and a cylinder overlapping its full height, axis exactly
    /// on the `x = 40` wall plane — assembled **past the standing gate** from production parts:
    /// `plane_index_setup_inner` (the gate-free half) plus the cylinder table built the way the
    /// gate builds it. Returns everything the two locks below read, and the record the
    /// gate-opening cell will produce: `(wall class, cylinder)` marked not-proven-clear.
    #[allow(clippy::type_complexity)]
    fn armed_through_boss() -> (
        Model,
        Handle<Solid>,
        Handle<Solid>,
        crate::planes::PlaneSetup,
        usize,
        std::collections::HashSet<(usize, usize)>,
    ) {
        armed_through_boss_z(-10.0, 50.0)
    }

    /// [`armed_through_boss`] with the boss's axial extent chosen: `z_lo` and height `h`. The
    /// default runs through the plate (`−10`, `50`); a boss whose lower cap sits **inside** the
    /// plate (`10`, `30`) has a lateral whose lower boundary is a **chain** — arcs at z = 10
    /// (outside the plate) and z = 20 (inside) joined by rulings — the staircase the corner
    /// rule is watched on.
    #[allow(clippy::type_complexity)]
    fn armed_through_boss_z(
        z_lo: f64,
        h: f64,
    ) -> (
        Model,
        Handle<Solid>,
        Handle<Solid>,
        crate::planes::PlaneSetup,
        usize,
        std::collections::HashSet<(usize, usize)>,
    ) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([40.0, 20.0, z_lo]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            h,
        );
        m.rebuild_adjacency();
        let (mut setup, cyl_surfs) =
            crate::planes::plane_index_setup_inner(&m, plate, boss).unwrap();
        for &surf in &cyl_surfs {
            let nacre_topo::SurfaceTruth::Cylinder { def, .. } = m.surface_truth(surf) else {
                unreachable!("a cylinder row carries a cylinder truth")
            };
            let nacre_geom::Surface::Cylinder(cache) = m.surface(surf) else {
                unreachable!("a cylinder truth carries a cylinder cache")
            };
            setup.cyls.push(crate::planes::WorkingCyl {
                surf,
                def: def.clone(),
                cache: *cache,
                owner: crate::planes::SolidSide::A,
            });
        }
        let wc = setup
            .geom
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[0] - 40.0).abs() < 1e-12))
            .expect("the x = 40 wall class");
        let crossings: std::collections::HashSet<(usize, usize)> = [(wc, 0)].into_iter().collect();
        (m, plate, boss, setup, wc, crossings)
    }

    /// ★ **The gate's record arms the rulings road — and only the record.** With the
    /// through-boss pair listed, the boss's trace on the wall class is the rectangle: two
    /// rulings (one per side, Branch-named ends, no plain segments) and two cap chords
    /// (`Lo`/`Hi` roots on the cap classes). With the record empty — every production call
    /// today — the same trace is empty: the negative control that pins "an empty record
    /// changes nothing", which is the very population an unconditionally-firing arm broke
    /// (14 arc-family tests red, measured).
    #[test]
    fn the_gates_record_arms_the_ruling_trace() {
        let (m, _plate, boss, setup, wc, crossings) = armed_through_boss();
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let mut tr = Trace::default();
        trace_one_of(
            &m,
            boss,
            SolidSide::B,
            wc,
            &jd,
            &setup.cyls,
            &setup.planes,
            &setup.surf_ix,
            &setup.inc_b,
            &setup.plane_ix,
            crossings.clone(),
            &mut tr,
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        assert_eq!(
            tr.segs.len(),
            2,
            "one chord per cap — a segment of the sweep since cell ⑩"
        );
        let mut sides: Vec<i8> = tr.rulings.iter().map(|r| r.side).collect();
        sides.sort_unstable();
        assert_eq!(sides, [-1, 1], "one ruling per side");
        for r in &tr.rulings {
            for nd in r.end {
                let (_, cyl, _) = combinatorics::branch_name(nd).expect("a Branch end");
                assert_eq!(cyl, 0);
            }
            assert!(matches!(r.kind, SegKind::Transversal { .. }));
        }
        for c in &tr.segs {
            assert_eq!(
                c.end_h,
                [combinatorics::EndPin::Cylinder; 2],
                "a chord's ends are the two roots"
            );
            // The chord rides the cap's own class — a ⊥ plane at z = −10 or z = 40.
            let z = setup.geom[c.wall].tri[0].as_array()[2];
            assert!(
                setup.geom[c.wall]
                    .tri
                    .iter()
                    .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
                    && ((z + 10.0).abs() < 1e-12 || (z - 40.0).abs() < 1e-12),
                "cap class at z = {z}"
            );
        }
        // The negative control: today's record.
        let mut tr0 = Trace::default();
        trace_one_of(
            &m,
            boss,
            SolidSide::B,
            wc,
            &jd,
            &setup.cyls,
            &setup.planes,
            &setup.surf_ix,
            &setup.inc_b,
            &setup.plane_ix,
            Default::default(),
            &mut tr0,
        );
        assert!(
            tr0.rulings.is_empty() && tr0.segs.is_empty(),
            "empty record, empty road"
        );
    }

    /// ★ **The armed arrangement digests the rectangle** — production bricks end to end on the
    /// wall class. Each ruling is cut at its T-junctions with the plate's `z = 0` and `z = 20`
    /// lines (three pieces each), those lines split at the same Branch nodes, the chords close
    /// the far ends, and the walk closes the subdivision: one unbounded contour and five
    /// bounded cells — plate-left, plate-right, the overlap band, and the rectangle's two
    /// overhangs.
    #[test]
    fn the_armed_wall_class_walks_to_closed_cells() {
        let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let tr = trace_on_class_of(
            &m,
            plate,
            boss,
            wc,
            &jd,
            &setup.cyls,
            &setup.planes,
            &setup.surf_ix,
            &setup.inc_a,
            &setup.inc_b,
            &setup.plane_ix,
            crossings,
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let split =
            split_at_crossings(&jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
        let split = drop_newsless(split).unwrap();
        let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
        assert!(circles.is_empty(), "a ∥ class carries no circles");
        let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
        assert_eq!(rulings.len(), 2);
        let edges = ClassEdges::of(
            &jd,
            &setup.cyls,
            wc,
            &split,
            &circles,
            &rulings,
            &Aliases::default(),
        )
        .unwrap();
        assert_eq!(
            edges.rulings.len(),
            6,
            "each ruling cut at z = 0 and z = 20"
        );
        let (cells, _face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
        assert_eq!(cells.len(), 6, "five bounded cells and the outer");
        assert_eq!(
            cells.iter().filter(|c| c.winding == -1).count(),
            1,
            "one connected skeleton, one outer contour"
        );
    }

    /// One armed class's arrangement, through the production bricks (the chain
    /// `the_armed_wall_class_walks_to_closed_cells` spells out) — for the locks that need
    /// several classes' products at once.
    fn armed_class_edges<'a>(
        m: &Model,
        plate: Handle<Solid>,
        boss: Handle<Solid>,
        setup: &'a crate::planes::PlaneSetup,
        jd: &Judge<'a, WorkingPlane>,
        wc: usize,
        crossings: &std::collections::HashSet<(usize, usize)>,
    ) -> ClassEdges<'static> {
        let tr = trace_on_class_of(
            m,
            plate,
            boss,
            wc,
            jd,
            &setup.cyls,
            &setup.planes,
            &setup.surf_ix,
            &setup.inc_a,
            &setup.inc_b,
            &setup.plane_ix,
            crossings.clone(),
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        let merged = merge_coincident(jd, &tr.segs, wc, &Aliases::default());
        let split =
            split_at_crossings(jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
        let split = drop_newsless(split).unwrap();
        let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
        let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
        let e = ClassEdges::of(
            jd,
            &setup.cyls,
            wc,
            &split,
            &circles,
            &rulings,
            &Aliases::default(),
        )
        .unwrap();
        // Owned copies so the borrows above may end — a test convenience, not a production shape.
        ClassEdges {
            segs: std::borrow::Cow::Owned(e.segs.into_owned()),
            arcs: std::borrow::Cow::Owned(e.arcs.into_owned()),
            rulings: std::borrow::Cow::Owned(e.rulings.into_owned()),
            circles: std::borrow::Cow::Owned(e.circles.into_owned()),
            cut_rims: e.cut_rims,
        }
    }

    /// **A lateral face answers the crossing question by the same parity, on its own chart**
    /// (cell ②-b) — the digon oracle's lateral twin, on the through-boss Fuse's lateral. The
    /// truth is known in world coordinates, so every crossing a rational ray makes with the
    /// cylinder is checked: rays `{y = y₀, z = z₀}` over a half-step lattice (two roots each,
    /// at x = 40 ± √(25 − (y₀ − 20)²) — irrational θ) and the **station column**
    /// `{x = 40, z = z₀}`, whose roots are exactly the rulings' points (40, 15) and (40, 25):
    /// on the ruling within the notch's z (a corner at its ends), on the face beyond it.
    ///
    /// Two builds. **Through** (caps at −10 and 40): a band between the caps' whole circles
    /// with the plate's notch as its one hole (the plate-side half-circle × z ∈ (0, 20), two
    /// arcs on the cut rims and two rulings on the wall x = 40). **Staircase** (caps at 10 and
    /// 40, the lower cap inside the plate): the lower boundary is a chain — the arc at z = 10
    /// outside the plate, the arc at z = 20 inside, the two rulings between — so at the
    /// station (40, 25) one arc ends `hi` and the other starts `lo`: the corner rule's one
    /// discriminating population. ☑ In a notch both arcs share their ends, and «both ends
    /// count» is invisible there — measured; the staircase is why the second build exists.
    /// The other station, (40, 15), is the seam (`add_cylinder`'s `ref_dir` is −y): a
    /// seam-incident root with an arc above it is the one tie the loops road keeps.
    ///
    /// ☑ The through build's rays were the **banded arm's shadow** before that arm was deleted
    /// (cell ②-b, 1a): on the two whole-circle rims alone, its «opposite sides of the two rim
    /// planes» and the loops road's «exactly one rim above» agreed on all 2,180 rays — every
    /// root, the graze on a rim included. The suite's own band hits were 0 (the gate refuses a
    /// wall inside the strip, so no probe ray crosses a band there); this lattice was the
    /// only evidence, and is why the arm could go.
    #[test]
    fn the_lateral_parity_agrees_with_the_notch_it_bounds() {
        for staircase in [false, true] {
            lateral_lattice(staircase);
        }
    }

    fn lateral_lattice(staircase: bool) {
        use nacre_scalar::quad::{CylinderMeet, plane_plane_cylinder, plane_side};
        use nacre_scalar::{Orient, Rat};
        let caps = if staircase {
            [10.0, 40.0]
        } else {
            [-10.0, 40.0]
        };
        let (m, plate, boss, setup, wc, crossings) =
            armed_through_boss_z(caps[0], caps[1] - caps[0]);
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let (curved, plane_faces, rows, _) =
            armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, caps);
        let fuse = crate::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &plane_faces,
            &curved,
            &rows,
        )
        .unwrap();
        assert_eq!(fuse.len(), 1);
        let cf = crate::boolean::comp_face(&jd, &setup.cyls, &fuse[0]).unwrap();
        let combinatorics::CompSurf::Cylinder(def) = &cf.surf else {
            panic!("a lateral")
        };
        let combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
            panic!("loops, got {:?}", cf.outer)
        };
        assert!(cf.inner.is_empty(), "a lateral's holes are among its loops");
        let rims: Vec<usize> = loops
            .iter()
            .filter_map(|l| match l {
                combinatorics::LateralLoop::Circle(c) => Some(*c),
                combinatorics::LateralLoop::Ring(_) => None,
            })
            .collect();
        if staircase {
            assert_eq!(rims.len(), 1, "the upper cap's whole circle");
            assert_eq!(loops.len(), 2, "and the chain below");
        } else {
            assert_eq!(rims.len(), 2, "the caps' two whole circles");
            assert_eq!(loops.len(), 3, "and the notch as one ring");
        }
        let (o, mm, r) = (def.origin(), def.dir(), def.radius());
        let rat = |k: i128| Rat::new(k, 2).unwrap();
        let ri = Rat::from_int;
        let zero = ri(0);
        let neg = |v: Rat| zero.checked_sub(v).unwrap();
        let x40 = [ri(1), zero, zero, neg(ri(40))];
        let y20 = [zero, ri(1), zero, neg(ri(20))];
        let (cap_lo, cap_hi) = (ri(caps[0] as i128), ri(caps[1] as i128));
        // The notch's z range: the plate's own (0, 20) through the boss, or — with the lower
        // cap inside the plate — the chain's step from the cap (10) to the plate's top (20).
        let (n_lo, n_hi) = (if staircase { cap_lo } else { zero }, ri(20));
        // The truth on the cylinder, by the root's side of the wall and the ray's z. `seam`:
        // the root is the station (40, 15), the seam generator; a seam-incident root with an
        // arc above it is the tie the loops road keeps (`arc_span`'s `SeamRoot`; z is asked
        // first, so with no arc above the rims still decide).
        let truth = |x_side: Orient, z0: Rat, seam: bool| -> Option<bool> {
            if seam && z0 < n_hi {
                return None;
            }
            if z0 >= cap_hi {
                return if z0 == cap_hi { None } else { Some(false) };
            }
            if z0 < cap_lo {
                return Some(false);
            }
            if z0 == cap_lo {
                // The lower cap: through the plate a whole rim (a tie everywhere); in the
                // staircase an arc outside the plate only — a tie there, a corner on the wall,
                // and nothing inside, where the face starts at the plate's top.
                return match (staircase, x_side) {
                    (true, Orient::Negative) => Some(false),
                    _ => None,
                };
            }
            match x_side {
                // Outside the plate: the band, whole.
                Orient::Positive => Some(true),
                // Inside the plate's footprint: the notch (a hole, or the chain's step), its
                // arcs the boundary — both arcs through the plate, the upper one alone in the
                // staircase (its lower boundary there is the cap, handled above).
                Orient::Negative => {
                    if (!staircase && z0 == n_lo) || z0 == n_hi {
                        None
                    } else {
                        Some(!(z0 > n_lo && z0 < n_hi))
                    }
                }
                // On the wall: a ruling for z within the notch (corners at its ends), the
                // face beyond.
                Orient::Zero => {
                    if z0 >= n_lo && z0 <= n_hi {
                        None
                    } else {
                        Some(true)
                    }
                }
            }
        };
        // on face, off (hole), boundary, tangent, station on-ruling, station on-face, seam ties
        let mut n = [0usize; 7];
        let mut ask = |pa: [Rat; 4], pb: [Rat; 4], z0: Rat, station: bool| {
            let roots = match plane_plane_cylinder(&pa, &pb, &o, &mm, r).unwrap() {
                CylinderMeet::Pair { line, s } => (line, s),
                CylinderMeet::Tangent { .. } => {
                    n[3] += 1;
                    return;
                }
                other => panic!("an ordinary ray: {other:?}"),
            };
            let (line, s) = roots;
            for root in &s {
                let x_side = plane_side(&x40, &line, root);
                let seam = station && plane_side(&y20, &line, root) == Orient::Negative;
                let want = truth(x_side, z0, seam);
                let got = combinatorics::loop_parity(&jd, def, loops, &line, root);
                assert_eq!(
                    got, want,
                    "staircase {staircase} z0 {z0:?} side {x_side:?} station {station} seam {seam}"
                );
                match (got, station, seam) {
                    (Some(true), false, _) => n[0] += 1,
                    (Some(false), false, _) => n[1] += 1,
                    (None, false, _) => n[2] += 1,
                    (None, true, true) if z0 < n_lo || z0 > n_hi => n[6] += 1,
                    (None, true, _) => n[4] += 1,
                    (Some(_), true, _) => n[5] += 1,
                }
            }
        };
        for k in 0..=108 {
            let z0 = ri(-12).checked_add(rat(k)).unwrap();
            let pz = [zero, zero, ri(1), neg(z0)];
            for j in 0..=20 {
                let y0 = ri(15).checked_add(rat(j)).unwrap();
                ask([zero, ri(1), zero, neg(y0)], pz, z0, false);
            }
            ask(x40, pz, z0, true);
        }
        eprintln!(
            "lateral lattice (staircase {staircase}): on {} hole {} boundary {} tangent {} \
             station on-ruling {} station on-face {} seam ties {}",
            n[0], n[1], n[2], n[3], n[4], n[5], n[6]
        );
        assert!(
            n.iter().all(|&c| c > 0),
            "every arm has a population: {n:?}"
        );
    }

    /// ★ **A cut circle is a band boundary** (rulings ladder, cell 2 — commit ①). On the armed
    /// through-boss, the per-class products of the four ⊥ classes are collected the production
    /// way (`per_class` on each), and the chart must break the lateral at the two **cut**
    /// circles (z = 0, z = 20) as well as the rims: three intervals, not one full-height band
    /// (`boundary_lines` and `chart_of`'s z-lines agree on the four). The middle interval is
    /// both-cut and is emitted as panel rings; blinding **both** of the chart's axes (the arc
    /// labels and the rulings') leaves the reader no answer there and `emit_lateral` refuses by
    /// name (`CylinderGateUndecided`) — blinding only the arcs does not, because the vertical
    /// lines answer in their place.
    /// The through-boss's curved carriers and plane faces, collected the production way
    /// (`per_class` on each of the four ⊥ classes and the wall class), with the cylinder's rows
    /// and the four ⊥ classes `[z0, z20, cap_lo, cap_hi]`. ★ This mirrors production's fold by
    /// hand (`trace_result_faces`' accumulation). A drift between them is not caught by
    /// anything: a test would simply start measuring a map the boolean never builds.
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    fn armed_curved(
        m: &Model,
        plate: Handle<Solid>,
        boss: Handle<Solid>,
        setup: &crate::planes::PlaneSetup,
        jd: &Judge<'_, WorkingPlane>,
        wc: usize,
        crossings: &std::collections::HashSet<(usize, usize)>,
        caps: [f64; 2],
    ) -> (
        Curved,
        Vec<LocalFace>,
        Vec<crate::bands::CylRow>,
        [usize; 4],
    ) {
        let z_class = |z: f64| -> usize {
            setup
                .geom
                .iter()
                .position(|p| p.tri.iter().all(|q| (q.as_array()[2] - z).abs() < 1e-12))
                .unwrap_or_else(|| panic!("a class at z = {z}"))
        };
        let (z0, z20, cap_lo, cap_hi) = (
            z_class(0.0),
            z_class(20.0),
            z_class(caps[0]),
            z_class(caps[1]),
        );
        // The production carriers: disk labels from every ⊥ class, cut rims from the cut ones.
        let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
        let mut cut_rims: CutRims = HashMap::new();
        let mut plane_faces: Vec<LocalFace> = Vec::new();
        let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
        // ★ The wall class too: the chart's rulings are the wall's pieces, and without them the
        // both-cut middle is one whole-circle cell the emitter cannot read per sector.
        let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
        for c in [z0, z20, cap_lo, cap_hi, wc] {
            let edges = armed_class_edges(m, plate, boss, setup, jd, c, crossings);
            let staged = per_class(jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
            for (cyl, label) in &staged.disk_labels {
                disk_labels.insert((*cyl, c), *label);
            }
            for (cyl, al) in &staged.arc_labels {
                arc_labels.entry((*cyl, c)).or_default().push(al.clone());
            }
            for (cyl, rim) in &edges.cut_rims {
                cut_rims.insert((*cyl, c), rim.clone());
            }
            for (cyl, r) in staged.ruling_extents {
                rulings.entry(cyl).or_default().push(r);
            }
            plane_faces.extend(staged.faces);
        }
        let curved = Curved {
            aliases: Aliases::default(),
            disk_labels,
            arc_labels,
            cut_rims,
            rulings,
        };
        let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
        (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi])
    }

    #[test]
    fn a_cut_circle_bounds_the_bands() {
        let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi]) =
            armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, [-10.0, 40.0]);
        let (disk_labels, arc_labels, cut_rims) =
            (&curved.disk_labels, &curved.arc_labels, &curved.cut_rims);
        assert!(cut_rims.contains_key(&(0, z0)) && cut_rims.contains_key(&(0, z20)));
        assert!(disk_labels.contains_key(&(0, cap_lo)) && disk_labels.contains_key(&(0, cap_hi)));
        assert_eq!(arc_labels[&(0, z0)].len(), 2, "two arcs, two sector labels");
        assert_eq!(arc_labels[&(0, z20)].len(), 2);
        assert_eq!(rows.len(), 1);
        // ★ A cut circle is a band boundary (the rulings ladder): the chart's boundary rule names
        // the two cut rims beside the caps, and the chart's own lines are exactly those four.
        let def = &setup.cyls[0].def;
        let t_of = |c: usize| crate::bands::axis_param(&jd, c, def).unwrap();
        let expect: Vec<Rat> = [cap_lo, z0, z20, cap_hi].map(t_of).to_vec();
        assert_eq!(
            crate::cyl_chart::boundary_lines(&jd, 0, def, &plane_faces, &curved, &rows).unwrap(),
            expect,
            "three intervals, cut circles included"
        );
        let chart = crate::cyl_chart::chart_of(&jd, &setup.cyls, 0, &plane_faces, &curved).unwrap();
        assert_eq!(
            chart.z_lines.iter().map(|l| l.t).collect::<Vec<_>>(),
            expect,
            "the chart's lines are the boundaries and nothing else here"
        );
        // ★ **The emitter answers with regions (D5).** Under `keep` the wall splits the middle
        // interval into two sectors, one kept and one not, and the kept one joins the whole
        // bands above and below it: **fuse** emits one face — a `Band` between the two cap
        // circles with the unkept sector as its one **hole** (a 4-node ring: two arcs on the
        // cut rims' own nodes, two rulings); **cut** keeps the other sector alone — one 4-node
        // `Ring` face. Which sector is the outer one is the assembly's and the volume oracle's to
        // measure (the gate-opening cell), not this harness's to re-derive.
        let faces_for = |kind: BoolKind| -> Vec<LocalFace> {
            crate::cyl_chart::emit_lateral(kind, &jd, &setup.cyls, &plane_faces, &curved, &rows)
                .unwrap()
        };
        let (fuse, cut) = (faces_for(BoolKind::Fuse), faces_for(BoolKind::Cut));
        assert_eq!(fuse.len(), 1, "fuse: one lateral face, a band with a hole");
        assert_eq!(cut.len(), 1, "cut: one lateral face, the kept sector");
        let rim = &cut_rims[&(0, z0)];
        // The ccw arc a ring carries on the lower cut rim, as an ordered node pair.
        let arc_on_z0 = |r: &crate::boolean::Ring| -> [NodeId; 2] {
            let n = r.nodes.len();
            assert_eq!(n, 4, "two arcs and two rulings");
            let arcs = r
                .walls
                .iter()
                .filter(|w| matches!(w, crate::boolean::Wall::Arc { .. }))
                .count();
            assert_eq!(arcs, 2);
            for i in 0..n {
                let (a, b) = (r.nodes[i], r.nodes[(i + 1) % n]);
                if let crate::boolean::Wall::Arc { ccw, .. } = r.walls[i]
                    && rim.nodes.contains(&a)
                    && rim.nodes.contains(&b)
                {
                    return if ccw { [a, b] } else { [b, a] };
                }
            }
            panic!("no arc on the lower cut rim");
        };
        let crate::boolean::Bound::Band { lo, hi } = &fuse[0].outer else {
            panic!("fuse: a band between the caps, got {:?}", fuse[0].outer);
        };
        assert!(
            matches!(
                (lo, hi),
                (
                    crate::boolean::Rim::Circle(_),
                    crate::boolean::Rim::Circle(_)
                )
            ),
            "the caps' whole circles are the band's rims"
        );
        assert_eq!(
            fuse[0].inner.len(),
            1,
            "the unkept sector is the band's one hole"
        );
        let crate::boolean::Bound::Ring(hole) = &fuse[0].inner[0] else {
            panic!("a hole is a ring");
        };
        let crate::boolean::Bound::Ring(panel) = &cut[0].outer else {
            panic!("cut: the kept sector is a ring, got {:?}", cut[0].outer);
        };
        assert!(cut[0].inner.is_empty());
        // ★ **Fuse's hole is cut's panel**: the sector fuse drops (inside the plate) is exactly
        // the sector cut keeps (the groove's wall), so the two rings carry the same ccw arc on
        // the lower cut rim — the old «complementary panels» claim, restated for regions.
        let (pf, pc) = (arc_on_z0(hole), arc_on_z0(panel));
        assert_eq!(pf, pc, "fuse's hole is cut's panel");
        // Negative control: without the sector labels the both-cut interval's cells have no
        // speaking end, and the emitter refuses the class rather than guess a chamber. Called
        // directly, not through the census — this hand-broken input violates the very premise
        // (`src0_present == 0`) the census asserts where the fact is made.
        let mut blind = Curved {
            aliases: Aliases::default(),
            disk_labels: curved.disk_labels.clone(),
            arc_labels: HashMap::new(),
            cut_rims: curved.cut_rims.clone(),
            rulings: curved.rulings.clone(),
        };
        blind.arc_labels.clear();
        // ★ **Blinding one axis is no longer blinding the reader** — the chart has two, and the
        // rulings answer where the rims are silent (the vertical read's own cell). So this now
        // says the weaker, truer thing: with the arcs gone the emitter still gets an answer, and
        // it is only when **both** axes are blinded that it refuses by name.
        assert!(
            crate::cyl_chart::emit_lateral(
                BoolKind::Fuse,
                &jd,
                &setup.cyls,
                &plane_faces,
                &blind,
                &rows
            )
            .is_ok(),
            "the vertical lines answer what the blinded rims cannot"
        );
        for v in blind.rulings.values_mut() {
            for r in v.iter_mut() {
                r.label = None;
            }
        }
        assert!(matches!(
            crate::cyl_chart::emit_lateral(
                BoolKind::Fuse,
                &jd,
                &setup.cyls,
                &plane_faces,
                &blind,
                &rows
            ),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            })
        ));
    }

    /// The armed through-boss fuse **assembled to the end of the road** — the pipeline the
    /// cell-3 lock walks, extracted so the consumer locks (props' θ-range integrals, tess's
    /// open-rim merge) measure the very same solid through one spelling: every class's
    /// arrangement, the coplanar unify, the band/panel pass, the seam table, `reconstruct`.
    /// Returns the model with the one welded solid live, plus the setup and the wall class.
    fn armed_assembled_through_boss() -> (Model, crate::planes::PlaneSetup, usize) {
        let (mut m, plate, boss, setup, wc, crossings) = armed_through_boss();
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let mut faces: Vec<LocalFace> = Vec::new();
        let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
        let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
        let mut cut_rims: CutRims = HashMap::new();
        let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
        for c in 0..setup.geom.len() {
            let edges = armed_class_edges(&m, plate, boss, &setup, &jd, c, &crossings);
            let staged = per_class(&jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
            for (cyl, label) in &staged.disk_labels {
                disk_labels.insert((*cyl, c), *label);
            }
            // ★ This mirrors production's fold by hand (`trace_result_faces`' accumulation). A
            // drift between them is not caught by anything: the test would simply start measuring
            // a map the boolean never builds.
            for (cyl, al) in &staged.arc_labels {
                arc_labels.entry((*cyl, c)).or_default().push(al.clone());
            }
            for (cyl, rim) in &edges.cut_rims {
                cut_rims.insert((*cyl, c), rim.clone());
            }
            for (cyl, r) in staged.ruling_extents {
                rulings.entry(cyl).or_default().push(r);
            }
            faces.extend(staged.faces);
        }
        let curved = Curved {
            aliases: Aliases::default(),
            disk_labels,
            arc_labels,
            cut_rims,
            rulings,
        };
        let faces = crate::boolean::unify_coplanar_faces(faces, &jd, &setup.cyls).unwrap();
        let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
        let mut faces = faces;
        faces.extend(
            crate::cyl_chart::emit_lateral(
                BoolKind::Fuse,
                &jd,
                &setup.cyls,
                &faces,
                &curved,
                &rows,
            )
            .unwrap(),
        );
        let seam = seam_table(&faces, &setup.cyls, &jd).unwrap();
        let out = crate::boolean::reconstruct(
            &mut m,
            &jd,
            &seam,
            &faces,
            &setup.cyls,
            &curved.cut_rims,
            None,
            crate::boolean::Tangencies::none(),
        )
        .unwrap();
        assert_eq!(out.len(), 1, "one welded solid");
        m.live_solids = out;
        m.rebuild_adjacency();
        (m, setup, wc)
    }

    /// ★ **The armed assembly welds the panels** (rulings ladder, cell 3 — the lock): the whole
    /// through-boss fuse, from production parts past the standing gate — every class's
    /// arrangement, the coplanar unify, the band/panel pass, the seam table, and
    /// `reconstruct` — comes back one solid with a clean `validate`. The ruling edges are
    /// pinned structurally: exactly two straight lateral edges (the kept outer panel's), each
    /// used exactly twice, carriers stated as **the edge's own fact** — the cylinder and the
    /// wall plane — and the outer panel's seam-holding arc is split at an `OnSeam` vertex
    /// (the panel road runs the wrap-arc split the Band road already had).
    #[test]
    fn the_armed_assembly_welds_the_panels() {
        let (m, setup, wc) = armed_assembled_through_boss();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        // The ruling edges, structurally: straight lateral edges = carriers {cylinder, a plane}
        // with a Line curve (the seam's [cyl, cyl] spelling is excluded by the mixed pair).
        let cyl_surf = setup.cyls[0].surf;
        let wall_surf = setup.geom[wc].surf;
        let reach = m.reachable();
        let mut rulings = 0;
        for (eh, e) in m.edges.iter() {
            if !reach.edges.contains(&eh) {
                continue;
            }
            let mixed = (e.surfaces[0] == cyl_surf) != (e.surfaces[1] == cyl_surf);
            if !mixed {
                continue;
            }
            let curve = m.derive_edge_curve(e.surfaces, e.vertices).unwrap();
            if matches!(curve, nacre_geom::Curve::Line(_)) {
                rulings += 1;
                assert!(
                    e.surfaces.contains(&wall_surf),
                    "a ruling's plane carrier is the wall: {:?}",
                    e.surfaces
                );
            }
        }
        assert_eq!(rulings, 2, "the kept outer panel's two rulings");
        // The seam split ran on the panel: an OnSeam vertex is reachable.
        let on_seam = reach
            .vertices
            .iter()
            .filter(|&&v| matches!(m.vertices.get(v).def, nacre_topo::VertexDef::OnSeam(_)))
            .count();
        assert!(on_seam >= 1, "the outer panel's wrap arc split at the seam");
    }

    /// ★ **The armed solid's mass properties are exact** (rulings ladder, cell 4 — the props
    /// instrument): `mass_props` on the assembled through-boss fuse, with the θ-range lateral
    /// integrals live, answers the derived closed forms — volume `32000 + 1000π` (plate plus
    /// the boss outside it), area `6200 + 425π` (walls 3000 + split x = 40 wall 600, caps
    /// bitten `3200 − 25π`, boss disks `50π`, full bands `300π`, the outer half-panel `100π`).
    /// The lifted probe's full-2π mis-answer (7849.34 = truth + the panel counted whole,
    /// `+100π`) cross-checks the area derivation. This is the volume oracle standing while the
    /// gate still refuses production input.
    #[test]
    fn the_armed_solids_mass_properties_are_exact() {
        let (m, _setup, _wc) = armed_assembled_through_boss();
        let props = nacre_props::mass_props(&m, m.live_solids[0]).unwrap();
        let volume = 32000.0 + 1000.0 * std::f64::consts::PI;
        let area = 6200.0 + 425.0 * std::f64::consts::PI;
        assert!(
            (props.volume - volume).abs() <= 1e-9 * volume,
            "volume {} != {volume}",
            props.volume
        );
        assert!(
            (props.area - area).abs() <= 1e-9 * area,
            "area {} != {area}",
            props.area
        );
    }

    /// ★ **The armed solid tessellates watertight** (rulings ladder, cell 4 — the tess
    /// instrument): the open-rim merge triangulates the θ-panels of the assembled through-boss
    /// fuse — the outer panel's rims are two arcs split at the seam vertex, so this one solid
    /// exercises the multi-polyline open chain *and* the seam-crossing unwrap — and every
    /// undirected triangle edge is used exactly twice.
    #[test]
    fn the_armed_solid_tessellates_watertight() {
        let (m, _setup, _wc) = armed_assembled_through_boss();
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).unwrap();
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "the mesh is watertight");
    }

    /// **The mixed parity reads a bitten ring — exactly, on the production pieces.**
    ///
    /// The straddling boss cuts the plate-top ring at `(40, 15)` and `(40, 25)`, so that ring
    /// carries two branch corners and one arc step; the overhang digon is chord + outer arc.
    /// Probes are derived, not read back: the chart for `n = +z` picks `e1 = [0, −1, 0]`, so
    /// the ray runs toward −y at fixed x. `[37, 30]` is inside and its ray crosses the **arc
    /// twice** (`(37−40)² + (y−20)² = 25` → y = 16, 24) before the bottom edge — the arc arm is
    /// what that probe measures, and a chord-minded parity would answer it wrong. `[37, 18]`
    /// sits inside the bite (outside the face), `[50, 20]` outside everything; `[43, 20]` is
    /// inside the overhang (one arc crossing at y = 16). Cells are identified structurally,
    /// not by size: the outside cell (winding −1) owns the outer arc piece's twin, so that
    /// piece's forward cell is the overhang digon; the other digon is the bite.
    #[test]
    fn the_mixed_parity_reads_a_bitten_ring() {
        with_bitten_rings(|jd, cyls, wc, big, overhang, bite| {
            let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
            let rat = |x: i128, y: i128| {
                [
                    nacre_scalar::Rat::from_int(x),
                    nacre_scalar::Rat::from_int(y),
                    nacre_scalar::Rat::from_int(20),
                ]
            };
            let ask = |ring: &[combinatorics::RingEdge], p: [nacre_scalar::Rat; 3]| {
                combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring)
            };
            assert_eq!(ask(big, rat(12, 12)), Some(true), "plain interior");
            assert_eq!(
                ask(big, rat(37, 30)),
                Some(true),
                "interior whose ray crosses the inner arc twice"
            );
            assert_eq!(ask(big, rat(37, 18)), Some(false), "inside the bite");
            assert_eq!(ask(big, rat(50, 20)), Some(false), "outside everything");
            assert_eq!(
                ask(overhang, rat(43, 20)),
                Some(true),
                "inside the overhang"
            );
            assert_eq!(
                ask(overhang, rat(37, 18)),
                Some(false),
                "the bite is not the overhang"
            );
            assert_eq!(ask(bite, rat(37, 18)), Some(true), "inside the bite digon");
            assert_eq!(
                ask(bite, rat(43, 20)),
                Some(false),
                "the overhang is not the bite"
            );
        });
    }

    /// ★★★★★ **The digon's truth is known independently, so the parity can be swept rather than
    /// spot-checked.** A digon of a chord and an arc is `disk ∩ half-space`, and both halves are
    /// exact rational predicates the arrangement already owns —
    /// [`nacre_scalar::quad::cylinder_radial_side`] and the sign of the wall's plane equation. A
    /// grid over the circle's neighbourhood therefore checks **every** answer, and it is the only
    /// control here that crosses the straight arm, the arc arm **and their junction**: the
    /// non-mixed road's oracle cannot see an arc, and a whole circle's arm collapses to
    /// `Ordering::Equal` on a single step.
    ///
    /// ☑ **What it does and does not reach, measured.** On this fixture the chord's ends are the
    /// seam and the tangent columns, and a ray along the chord meets **both** arc ends at once —
    /// so this sweep is green under a flipped arc-end departure sign (two ends on one ray flip
    /// together and keep the parity) and green under the straight arm's old corner abstention
    /// (cell ②'s negative controls, both measured). It is the oracle for the mixed road as a
    /// whole — the straight arm, the arc arm and their junction on ordinary roots — and the
    /// corner-bitten plate (`a_root_at_an_arc_end_is_a_corner_on_the_ray`) is the watch on the
    /// arc-end arm, where a ray meets one end alone.
    #[test]
    fn the_mixed_parity_agrees_with_the_digon_it_bounds() {
        use nacre_scalar::{Orient, Rat};
        with_bitten_rings(|jd, cyls, wc, _big, overhang, bite| {
            let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
            let def = cyls
                .iter()
                .map(|c| &c.def)
                .find(|d| d.radius() == Rat::from_int(5))
                .expect("the bitten circle");
            let (mut swept, mut on_boundary) = (0usize, 0usize);
            let (mut abstained, mut inside_seen) = (0usize, 0usize);
            let mut on_chord = 0usize;
            for i in 60..=100 {
                for j in 20..=60 {
                    let p = [
                        Rat::new(i.into(), 2).unwrap(),
                        Rat::new(j.into(), 2).unwrap(),
                        Rat::from_int(20),
                    ];
                    let radial = nacre_scalar::quad::cylinder_radial_side(
                        &p,
                        &def.origin(),
                        &def.dir(),
                        def.radius(),
                    );
                    // The chord is the plate's wall `x = 40`; the bite keeps `x < 40`.
                    let wall = p[0].checked_sub(Rat::from_int(40)).unwrap();
                    if radial == Orient::Zero || wall == Rat::from_int(0) {
                        on_boundary += 1;
                        // ★★★★★ **A point strictly inside the circle and *on* the chord is on
                        // both digons' boundary, and "inside" has no answer there.** These are the
                        // grid's sharpest points: the ray along the chart's first axis leaves such
                        // a point straight through **both arc endpoints** — the chord is a step
                        // lying along the ray, and the straight arm names the probe between its
                        // ends as the boundary before any arc is asked (cell ②).
                        if radial == Orient::Negative && wall == Rat::from_int(0) {
                            for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                                assert_eq!(
                                    combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                    None,
                                    "{who} must abstain on its own chord at {p:?}"
                                );
                            }
                            on_chord += 1;
                        }
                        continue;
                    }
                    let inside_bite = radial == Orient::Negative && wall < Rat::from_int(0);
                    let inside_over = radial == Orient::Negative && wall > Rat::from_int(0);
                    for (ring, want, who) in [
                        (bite, inside_bite, "bite"),
                        (overhang, inside_over, "overhang"),
                    ] {
                        // ★ **An abstention is allowed and a wrong answer is not.** The ray runs
                        // along the chart's first axis, so a grid point whose ray meets a ring
                        // corner has no parity to report — the caller's remedy is another point,
                        // which is exactly what `coord_probes` does with its candidate list. What
                        // the sweep locks is that every answer it *does* give is the truth.
                        match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                            Some(got) => {
                                assert_eq!(got, want, "{who} at {p:?}");
                                swept += 1;
                                inside_seen += usize::from(want);
                            }
                            None => abstained += 1,
                        }
                    }
                }
            }
            // Neither vacuous nor all-outside: the grid straddles the circle and the chord, and
            // the two digons between them own a real interior.
            assert!(swept > 1_000, "answers swept: {swept}");
            assert!(inside_seen > 100, "points inside a digon: {inside_seen}");
            assert!(
                on_boundary > 0,
                "the grid meets the boundary: {on_boundary}"
            );
            // Recorded rather than bounded: every abstention here is the tangent ray (the rows
            // `y = 15, 25` graze the circle — measured, cell ②: 160 of 160, no corner among
            // them), and its count is a property of this grid, not of the rule.
            assert!(abstained > 0, "abstentions: {abstained}");
            // What those abstentions were, by kind.
            {
                let rows = combinatorics::tie_probe::ROWS
                    .lock()
                    .expect("the probe's lock is never held across a panic");
                let me = std::thread::current().name().unwrap_or("?").to_string();
                let mut hist: Vec<(combinatorics::tie_probe::Tie, usize)> = Vec::new();
                for (_, t) in rows.iter().filter(|(n, _)| *n == me) {
                    match hist.iter_mut().find(|(k, _)| k == t) {
                        Some((_, c)) => *c += 1,
                        None => hist.push((*t, 1)),
                    }
                }
                eprintln!("P2 digon abstained {abstained} on_chord {on_chord} kinds {hist:?}");
            }
            assert!(
                on_chord > 0,
                "points on the chord inside the circle: {on_chord}"
            );
        });
    }

    /// The bitten-plate fixture, handed to `f` as `(judge, cyls, class, big, overhang, bite)`.
    /// The bitten fixture: one wall (`x = 40`) cuts the circle into two digons.
    fn with_bitten_rings(
        f: impl FnOnce(
            &Judge<'_, WorkingPlane>,
            &[crate::planes::WorkingCyl],
            usize,
            &[combinatorics::RingEdge],
            &[combinatorics::RingEdge],
            &[combinatorics::RingEdge],
        ),
    ) {
        with_cut_rings([40.0, 40.0, 20.0], false, f);
    }

    /// A plate `[0, plate]` with a z-cylinder (r = 5, h = 10) standing on its top at
    /// `(40, 20)`: the class carrying the cut circle, its judge, and the three bounded rings —
    /// the bitten top, the overhang (disk ∖ plate) and the bite (disk ∩ plate). The plate's
    /// extent picks the cut: `[40, 40, 20]` bites with one wall (two digons); `[40, 24, 20]`
    /// with two — the corner `(40, 24)` sits inside the circle and both rings are trigons
    /// whose arc ends at the rational `(37, 24)`, off the seam and off the tangent columns.
    /// `seam_off` states the cylinder exactly with `ref_dir = x`, so the seam sits at
    /// `(45, 20)` — on no ring corner — instead of `add_cylinder`'s `−y`, which puts it on the
    /// chord end `(40, 15)`.
    fn with_cut_rings(
        plate: [f64; 3],
        seam_off: bool,
        f: impl FnOnce(
            &Judge<'_, WorkingPlane>,
            &[crate::planes::WorkingCyl],
            usize,
            &[combinatorics::RingEdge],
            &[combinatorics::RingEdge],
            &[combinatorics::RingEdge],
        ),
    ) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(plate));
        let b = if seam_off {
            use nacre_scalar::Rat;
            let r = Rat::from_int;
            m.add_cylinder_exact(
                [r(40), r(20), r(20)],
                [r(0), r(0), r(1)],
                [r(1), r(0), r(0)],
                r(5),
                r(10),
                None,
            )
            .expect("an exact cylinder on an axis frame")
            .0
        } else {
            m.add_cylinder(
                Point3::from_array([40.0, 20.0, 20.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                5.0,
                10.0,
            )
        };
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        // The one class whose circle is cut: find it by running the split everywhere.
        let mut found = None;
        for wc in 0..planes.len() {
            if !matches!(plane_ix.get(wc), Some(ClassIx::Plane(_)) | None) && wc < plane_ix.len() {
                continue;
            }
            if combinatorics::class_coeffs_rat(&jd, wc).is_none() {
                continue;
            }
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
            let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default())
            else {
                continue;
            };
            let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
                continue;
            };
            let Ok(edges) =
                ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
            else {
                continue;
            };
            if edges.arcs.is_empty() {
                continue;
            }
            let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
            let ns = edges.segs.len();
            let na = edges.arcs.len();
            // Four cells share this class: the bitten plate-top face (+1, many steps), the
            // outside (−1 — it also borders the arc, via the outer bulge), and two digons
            // (chord + arc): the overhang (disk minus plate) and the bite (disk ∩ plate).
            // Tell them apart structurally: the outside owns the outer piece's twin, so the
            // outer piece's forward cell is the overhang; the remaining digon is the bite.
            let outside = cells
                .iter()
                .position(|c| c.winding == -1)
                .expect("the outside cell");
            let bounded: Vec<usize> = (0..cells.len())
                .filter(|&i| cells[i].winding == 1)
                .collect();
            assert_eq!(
                bounded.len(),
                3,
                "the bitten top, the overhang and the bite"
            );
            let big_ix = *bounded
                .iter()
                .max_by_key(|&&i| cells[i].half_edges.len())
                .expect("the bitten top ring");
            let outer_ai = (0..na)
                .find(|ai| face_of[&(2 * (ns + ai) + 1)] == outside)
                .expect("the outer piece borders the outside");
            let overhang_ix = face_of[&(2 * (ns + outer_ai))];
            let bite_ix = *bounded
                .iter()
                .find(|&&i| i != big_ix && i != overhang_ix)
                .expect("the bite");
            assert_eq!(
                cells[overhang_ix].half_edges.len(),
                cells[bite_ix].half_edges.len(),
                "the same walls cut the overhang and the bite"
            );
            let ring = |ix: usize| -> Vec<combinatorics::RingEdge> {
                cells[ix]
                    .half_edges
                    .iter()
                    .map(|&he| edges.edge_at(he))
                    .collect()
            };
            found = Some((wc, ring(big_ix), ring(overhang_ix), ring(bite_ix)));
            break;
        }
        let (wc, big, overhang, bite) = found.expect("one class carries the cut circle");
        f(&jd, &cyls, wc, &big, &overhang, &bite);
    }

    /// **A root at an arc's own end is a corner on the ray, and the arc's tangent there says
    /// whether the arc counts it** (cell ②) — the half-open rule in the arc arm, on the one
    /// fixture that reaches it.
    ///
    /// The bitten fixture cannot: its chord's ends are the seam and the tangent columns, and
    /// each abstains first (measured over the whole suite before this cell: `ArcRootAtEnd` 0
    /// while the straight arm's corner tie spoke for the same corners). The plate's second wall
    /// `y = 24` puts the corner `(40, 24)` inside the circle, so the bite's arc ends at
    /// `(37, 24)` — a 3-4-5 point, rational, off the seam `(40, 15)` and off the tangent
    /// columns `x = 35, 45` — and the lattice column `x = 37` sends its rays through that end
    /// **alone**, beside an ordinary second root at `(37, 16)`. That is the single-end crossing
    /// a global sign error cannot hide: two ends on one ray flip together and keep the parity
    /// (the bitten fixture's chord), one end on a ray does not. Every answer the sweep gives
    /// is the truth, the boundary abstains, and the arm decided the end exactly once per point
    /// on the shooting side of that column, for both rings — under both seam placements, so
    /// the seam-incident arms and the `(false, false)` arm each decide an end.
    #[test]
    fn a_root_at_an_arc_end_is_a_corner_on_the_ray() {
        use nacre_scalar::{Orient, Rat};
        for seam_off in [false, true] {
            with_cut_rings(
                [40.0, 24.0, 20.0],
                seam_off,
                |jd, cyls, wc, _big, overhang, bite| {
                    let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
                    let def = cyls
                        .iter()
                        .map(|c| &c.def)
                        .find(|d| d.radius() == Rat::from_int(5))
                        .expect("the cut circle");
                    let decided0 = combinatorics::tie_probe::arc_end_decisions_here();
                    let (mut swept, mut abstained, mut inside_seen, mut boundary) =
                        (0usize, 0, 0, 0);
                    let mut column = [0usize; 2];
                    let rat = |k: i128| Rat::new(k, 2).unwrap();
                    for i in 60..=100 {
                        for j in 20..=60 {
                            let p = [rat(i), rat(j), Rat::from_int(20)];
                            let radial = nacre_scalar::quad::cylinder_radial_side(
                                &p,
                                &def.origin(),
                                &def.dir(),
                                def.radius(),
                            );
                            let (x, y) = (p[0], p[1]);
                            let (wx, wy) = (x == Rat::from_int(40), y == Rat::from_int(24));
                            // Both rings' boundary: the circle and the two chords, `x = 40` for
                            // `15 ≤ y ≤ 24` and `y = 24` for `37 ≤ x ≤ 40`, ends included.
                            let on_chord = (wx && y >= Rat::from_int(15) && y <= Rat::from_int(24))
                                || (wy && x >= Rat::from_int(37) && x <= Rat::from_int(40));
                            if radial == Orient::Zero || on_chord {
                                boundary += 1;
                                for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                                    assert_eq!(
                                        combinatorics::point_in_mixed_ring(
                                            jd, cyls, &coeffs, &p, ring
                                        ),
                                        None,
                                        "{who} must abstain on its boundary at {p:?}"
                                    );
                                }
                                continue;
                            }
                            let in_disk = radial == Orient::Negative;
                            let on_plate = x < Rat::from_int(40) && y < Rat::from_int(24);
                            for (k, (ring, want, who)) in [
                                (bite, in_disk && on_plate, "bite"),
                                (overhang, in_disk && !on_plate, "overhang"),
                            ]
                            .into_iter()
                            .enumerate()
                            {
                                match combinatorics::point_in_mixed_ring(
                                    jd, cyls, &coeffs, &p, ring,
                                ) {
                                    Some(got) => {
                                        assert_eq!(got, want, "{who} at {p:?}");
                                        swept += 1;
                                        inside_seen += usize::from(want);
                                        if x == Rat::from_int(37) && radial == Orient::Positive {
                                            column[k] += 1;
                                        }
                                    }
                                    None => abstained += 1,
                                }
                            }
                        }
                    }
                    assert!(swept > 1_000, "answers swept: {swept}");
                    assert!(inside_seen > 100, "points inside a ring: {inside_seen}");
                    assert!(boundary > 0, "the grid meets the boundary: {boundary}");
                    // The `x = 37` column outside the circle: 24 lattice points, every one answered by
                    // both rings — the shooting side through the end, the other side by a miss — and
                    // the end decided exactly once per point on the shooting side: 12 × 2 rings.
                    assert_eq!(column, [24, 24], "the single-end column is answered");
                    let decided = combinatorics::tie_probe::arc_end_decisions_here() - decided0;
                    eprintln!(
                        "arc-end sweep (seam_off {seam_off}): swept {swept} abstained {abstained} \
                 boundary {boundary} arc-end decisions {decided}"
                    );
                    // With the seam on the chord end `(40, 15)` (the `add_cylinder` build) the
                    // `(true, false)` / `(false, true)` arms decide the `x = 37` column's shooting side:
                    // 12 points × 2 rings. With the seam off every corner (`ref_dir = x`) both ends are
                    // ordinary and the `(false, false)` arm decides — and the `x = 40` column's rays now
                    // meet `(40, 15)` alone as well (its other root `(40, 25)` is the overhang's
                    // interior), 12 more per ring, one of them a call that then abstains at the corner
                    // `(40, 24)`; the column's seam-root abstentions are gone (182 → 160).
                    assert_eq!(
                        decided,
                        if seam_off { 48 } else { 24 },
                        "the arc-end arm decided the single-end rays (seam_off {seam_off})"
                    );
                },
            );
        }
    }

    /// **On a ring whose every corner is rational, the mixed road and the chart road are one
    /// rule** (cell ②): `point_in_ring_2d_rat` says Inside/Outside/OnBoundary and
    /// `point_in_mixed_ring` `Some(true)`/`Some(false)`/`None`, and the pairs match at every
    /// lattice point of every bounded ring of every class of two overlapping plates — squares
    /// and L-shapes with reflex corners. A lattice at half steps shares a row with every corner
    /// and every horizontal step, so every corner-on-the-ray configuration the half-open rule
    /// distinguishes (both steps up, both down, one each, a step along the ray, the probe at
    /// the corner) is on the sweep. Before this cell the mixed road abstained at every such
    /// corner (5 abstentions where the chart road answered, over the whole suite's rational
    /// rings; 0 disagreements — P3); a whole-suite shadow of the two roads cost 61 % of the
    /// serial sweep and is not kept — this lattice is the lock.
    #[test]
    fn the_two_roads_agree_on_every_rational_ring() {
        use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
        use nacre_scalar::Rat;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 2.0, 0.0]),
            Point3::from_array([6.0, 6.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let (mut rings_swept, mut points, mut inside, mut boundary, mut reflex) =
            (0usize, 0usize, 0usize, 0usize, 0usize);
        for wc in 0..planes.len() {
            let Some(coeffs) = combinatorics::class_coeffs_rat(&jd, wc) else {
                continue;
            };
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
            let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default())
            else {
                continue;
            };
            let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
                continue;
            };
            let Ok(edges) =
                ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
            else {
                continue;
            };
            let Ok((cells, _)) = walk_cells(&jd, &cyls, wc, &edges) else {
                continue;
            };
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            let chart = combinatorics::Chart2dRat::of_normal(&n).unwrap();
            let (e1, e2) = chart.axes();
            // Axis-aligned classes: unit chart axes and a unit normal, so a chart point
            // `(X, Y)` lifts to `X·e1 + Y·e2 − d·n`.
            for e in [e1, e2, &n] {
                assert_eq!(combinatorics::dot3_rat(e, e), Some(Rat::from_int(1)));
            }
            let lift = |xy: [Rat; 2]| -> [Rat; 3] {
                let mut p = [Rat::from_int(0); 3];
                for k in 0..3 {
                    p[k] = xy[0]
                        .checked_mul(e1[k])
                        .unwrap()
                        .checked_add(xy[1].checked_mul(e2[k]).unwrap())
                        .unwrap()
                        .checked_sub(coeffs[3].checked_mul(n[k]).unwrap())
                        .unwrap();
                }
                p
            };
            for c in cells.iter().filter(|c| c.winding == 1) {
                let ring: Vec<combinatorics::RingEdge> =
                    c.half_edges.iter().map(|&he| edges.edge_at(he)).collect();
                assert!(!combinatorics::ring_is_mixed(&ring), "a planar fixture");
                let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
                let ring2 = chart.ring(&jd, &nodes).expect("rational corners");
                if ring2.len() > 4 {
                    reflex += 1;
                }
                let bound = |k: usize, lo: bool| -> i128 {
                    let it = ring2.iter().map(|q| {
                        assert_eq!(q[k].denom(), 1, "integer corners");
                        q[k].numer()
                    });
                    if lo {
                        it.min().unwrap() - 1
                    } else {
                        it.max().unwrap() + 1
                    }
                };
                for xi in 2 * bound(0, true)..=2 * bound(0, false) {
                    for yi in 2 * bound(1, true)..=2 * bound(1, false) {
                        let q = [Rat::new(xi, 2).unwrap(), Rat::new(yi, 2).unwrap()];
                        let p = lift(q);
                        assert_eq!(
                            chart.project(&p),
                            Some(q),
                            "the lift is the chart's inverse"
                        );
                        let plane = point_in_ring_2d_rat(q, &ring2);
                        let mixed =
                            combinatorics::point_in_mixed_ring(&jd, &cyls, &coeffs, &p, &ring);
                        match (plane, mixed) {
                            (RingSide::Inside, Some(true)) => inside += 1,
                            (RingSide::Outside, Some(false)) => {}
                            (RingSide::OnBoundary, None) => boundary += 1,
                            other => panic!("class {wc} ring {ring2:?} at {q:?}: {other:?}"),
                        }
                        points += 1;
                    }
                }
                rings_swept += 1;
            }
        }
        assert!(rings_swept >= 3, "rings swept: {rings_swept}");
        assert!(reflex > 0, "an L-shaped ring: {reflex}");
        assert!(
            inside > 0 && boundary > 0,
            "points {points} inside {inside} boundary {boundary}"
        );
        eprintln!(
            "two roads: rings {rings_swept} points {points} inside {inside} boundary {boundary}"
        );
    }

    /// The gated drill population's fixture: a `[0,2]³` box and an axis-aligned cylinder at
    /// `(1,1)`, r=0.5 — every wall is a full unit from the axis, so the population gate passes
    /// and [`plane_index_setup`] hands the arrangement bricks a cylinder-bearing table. (It stood
    /// behind a stopper until C4b removed it; the fixture outlived the stopper.)
    fn drilled(m: &mut Model, z0: f64, h: f64) -> (Handle<Solid>, Handle<Solid>) {
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cylinder(
            Point3::from_array([1.0, 1.0, z0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
        m.rebuild_adjacency();
        (a, b)
    }

    /// The plane class whose defining triangle lies wholly at `z` — a z-cap.
    fn class_at_z(planes: &[WorkingPlane], z: f64) -> usize {
        planes
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[2] - z).abs() < 1e-12))
            .expect("a z-cap class")
    }

    /// A cylinder cap alone on its class: the disk's circular outer traces as a **seated
    /// circle** — and the lateral, whose rim lies on this very plane, adds its **graze** beside it.
    /// The cell bricks turn the pair into a disk `+1` / contour `−1` pair whose labels say "body on
    /// the cap's inside", with no segments anywhere. The emitted face's outer is the circle itself
    /// ([`crate::boolean::Bound::Circle`]), the vocabulary C4b's assembly will consume.
    #[test]
    fn a_cylinder_cap_is_a_seated_circle_and_a_disk_face() {
        let mut m = Model::new();
        let (a, b) = drilled(&mut m, -1.0, 4.0); // caps at z=-1 and z=3, clear of the box
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_b,
            plane_ix,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = class_at_z(&planes, 3.0);

        let mut tr = Trace::default();
        trace_one_of(
            &m,
            b,
            SolidSide::B,
            wc,
            &jd,
            &cyls,
            &faces_tab,
            &surf_ix,
            &inc_b,
            &plane_ix,
            Default::default(),
            &mut tr,
        );
        assert!(
            tr.segs.is_empty(),
            "a circle owes the segment machinery nothing"
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        // The class's normal is the cap's own outward +z (the cap is its only member), so the
        // body — below z=3 — is `body_above: false`. ★ **Two contributions, not one**: the cap is
        // seated here and the lateral's rim grazes here, and on this convex cap they agree.
        assert!(planes[wc].plane.normal().as_array()[2] > 0.0);
        let mut kinds: Vec<SegKind> = tr
            .circles
            .iter()
            .inspect(|c| assert!(c.cyl == 0 && c.solid == SolidSide::B, "{:?}", tr.circles))
            .map(|c| c.kind)
            .collect();
        kinds.sort_by_key(|k| matches!(k, SegKind::Graze { .. }));
        assert!(
            matches!(
                kinds[..],
                [
                    SegKind::Seated { body_above: false },
                    SegKind::Graze { body_above: false },
                ]
            ),
            "{:?}",
            tr.circles
        );

        let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
        let edges =
            ClassEdges::of(&jd, &cyls, wc, &[], &circles, &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
        // Pseudo-half-edges 0 and 1 (no segments): the disk (+1) and its contour (−1).
        assert_eq!(cells.len(), 2, "{cells:?}");
        assert_eq!(
            (cells[0].half_edges.as_slice(), cells[0].winding),
            (&[0][..], 1)
        );
        assert_eq!(
            (cells[1].half_edges.as_slice(), cells[1].winding),
            (&[1][..], -1)
        );
        let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
        assert_eq!(
            nesting.root_groups,
            vec![1],
            "the contour bounds the unbounded region"
        );
        assert!(nesting.holes.is_empty());
        let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
        assert_eq!(labels[1], [false; 4]);
        assert_eq!(
            labels[0],
            [false, false, false, true],
            "inside the circle, B below the plane only"
        );
        // Fuse keeps below and not above across the disk → the disk is a result face, and its
        // outer boundary is the circle — no ring, no nodes.
        let (out, _, _) = emit_faces(
            BoolKind::Fuse,
            &labels,
            &cells,
            &edges,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(matches!(
            out[0].outer,
            crate::boolean::Bound::Circle { cyl: 0 }
        ));
        assert!(out[0].inner.is_empty());
    }

    /// The through-drill's cap arrangement, brick by brick: the box cap is a seated 4-ring, the
    /// lateral crosses the cap plane in a **transversal circle**, and the bricks nest the circle
    /// as the cap cell's hole, label the disk with the cylinder straddling `W`, and emit — for
    /// `Cut` — exactly one face: the cap with a circular hole. The per-kind `edge_mask`
    /// difference is what the two labels measure (seated flips one side, transversal flips
    /// both).
    #[test]
    fn a_drill_circle_is_a_hole_of_the_cap_ring() {
        let mut m = Model::new();
        let (a, b) = drilled(&mut m, -1.0, 4.0); // z∈[-1,3]: through both caps of the box
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            cyls,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = class_at_z(&planes, 0.0);

        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        assert!(
            matches!(
                tr.circles[..],
                [CircleTrace {
                    cyl: 0,
                    solid: SolidSide::B,
                    kind: SegKind::Transversal { .. },
                    // No hole in this band, so the mark is the whole circle.
                    arc: None,
                }]
            ),
            "{:?}",
            tr.circles
        );

        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
        let split = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()).unwrap();
        assert_eq!(
            split.len(),
            4,
            "the cap ring alone — the circle owes the splitter nothing"
        );

        let edges =
            ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default()).unwrap();
        let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
        assert_eq!(cells.len(), 4, "cap ±1 and circle ±1: {cells:?}");
        let at = |he: usize| cells.iter().position(|c| c.half_edges == [he]).unwrap();
        let (disk, contour) = (at(2 * split.len()), at(2 * split.len() + 1));
        let cap = cells
            .iter()
            .position(|c| c.winding == 1 && c.half_edges.len() == 4)
            .unwrap();

        let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
        assert_eq!(
            nesting.holes.get(&cap).map(Vec::as_slice),
            Some(&[contour][..]),
            "the circle contour is the cap cell's hole"
        );

        let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
        // The class is the box's **bottom** cap, so which of `W`'s two sides carries the box is
        // read off the class normal rather than assumed: the box occupies z>0.
        let up = planes[wc].plane.normal().as_array()[2] > 0.0;
        let box_side = usize::from(!up); // 0 = above `n_W`, 1 = below
        let mut cap_label = [false; 4];
        cap_label[box_side] = true;
        assert_eq!(
            labels[cap], cap_label,
            "outside the circle: box material on the z>0 side only"
        );
        assert_eq!(
            labels[contour], labels[cap],
            "a hole bounds its host's region"
        );
        let mut disk_label = cap_label;
        disk_label[2] = true;
        disk_label[3] = true;
        assert_eq!(
            labels[disk], disk_label,
            "inside the circle the cylinder straddles W, the box side is unchanged"
        );

        let (out, _, _) = emit_faces(
            BoolKind::Cut,
            &labels,
            &cells,
            &edges,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(out.len(), 1, "one face: the drilled cap — {out:?}");
        assert_eq!(out[0].outer.expect_ring().len(), 4);
        assert!(
            matches!(out[0].inner[..], [crate::boolean::Bound::Circle { cyl: 0 }]),
            "{:?}",
            out[0].inner
        );
    }

    /// The hole arm of the seated tracer: a face whose inner loop is a circle (a two-hole
    /// plate's input shape — no producer builds one until C4b/C5, so the loops are doctored by
    /// hand) emits its polygon segments **and** a seated circle that inherits the face's own
    /// body side.
    #[test]
    fn a_circular_hole_ring_traces_as_a_seated_circle() {
        let mut m = Model::new();
        let (a, b) = drilled(&mut m, -1.0, 4.0);
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = class_at_z(&planes, 0.0);

        let input = combinatorics::trace_input(
            &m,
            [(a, &inc_a), (b, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
            &[],
            Default::default(),
        );
        let (fp, fl) = input.faces[0]
            .iter()
            .find(|(fp, _)| matches!(plane_ix[*fp], ClassIx::Plane(c) if c == wc))
            .expect("the box's bottom cap sits on wc");
        let doctored = combinatorics::FaceLoops {
            outer: fl.outer.clone(),
            holes: Some(vec![combinatorics::LoopRing::Circle { cyl: 0 }]),
            cycles: None,
        };
        let mut tr = Trace::default();
        trace_one(
            &[(*fp, doctored)],
            SolidSide::A,
            wc,
            &jd,
            &[],
            &faces_tab,
            &plane_ix,
            &Default::default(),
            &Aliases::default(),
            &mut tr,
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        assert_eq!(tr.segs.len(), 4, "the polygon outer still emits its edges");
        let SegKind::Seated { body_above } = tr.segs[0].kind else {
            panic!("a face on wc traces seated: {:?}", tr.segs[0]);
        };
        assert!(
            matches!(
                tr.circles[..],
                [CircleTrace {
                    cyl: 0,
                    solid: SolidSide::A,
                    kind: SegKind::Seated { body_above: ba },
                    arc: None
                }]
                if ba == body_above
            ),
            "the hole circle inherits the face's body side: {:?}",
            tr.circles
        );
    }

    /// The existence condition, negatively: a cap plane **outside** the lateral's rim span gets
    /// no circle (and no decline — a miss, like a parallel plane), while the rim-interior cap
    /// still does, and a rim-**coincident** plane carries the cap's seated circle rather than a
    /// transversal one. A ghost circle here would corrupt every label on the class.
    #[test]
    fn no_ghost_circle_outside_the_rim_span() {
        let mut m = Model::new();
        // z∈[-1,1.5]: through the box's bottom cap, short of its top cap at z=2.
        let (a, b) = drilled(&mut m, -1.0, 2.5);
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);

        // The box's top cap at z=2 is past the rim span [−1, 1.5]: no circle, no decline.
        let top = trace_on_class_of(
            &m,
            a,
            b,
            class_at_z(&planes, 2.0),
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert!(top.circles.is_empty(), "{:?}", top.circles);
        assert!(top.declined.is_empty(), "{:?}", top.declined);

        // The bottom cap at z=0 is strictly inside the span: the transversal circle is there.
        let bottom = trace_on_class_of(
            &m,
            a,
            b,
            class_at_z(&planes, 0.0),
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert_eq!(
            bottom
                .circles
                .iter()
                .filter(|c| matches!(c.kind, SegKind::Transversal { .. }))
                .count(),
            1,
            "{:?}",
            bottom.circles
        );

        // The cylinder's own top cap at z=1.5 (inside the box): the rim-coincident plane carries
        // **two** contributions — the cap's `Seated` circle and the lateral's `Graze`, because a
        // rim is a touch, not a miss. It adds no *transversal* twin (t = span end, not strictly
        // inside). ★ On a **convex** cap the two agree about which side the body is on; the whole
        // point of carrying both is the reflex corner (a blind bore's ceiling) where they do not,
        // and `edge_mask`'s `Graze > Seated` then picks the right one.
        let cap = trace_on_class_of(
            &m,
            a,
            b,
            class_at_z(&planes, 1.5),
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert!(cap.declined.is_empty(), "{:?}", cap.declined);
        let sides: Vec<(bool, bool)> = cap
            .circles
            .iter()
            .filter_map(|c| match c.kind {
                SegKind::Seated { body_above } => Some((false, body_above)),
                SegKind::Graze { body_above } => Some((true, body_above)),
                SegKind::Transversal { .. } | SegKind::Tangent { .. } => None,
            })
            .collect();
        assert_eq!(sides.len(), 2, "seated and graze: {:?}", cap.circles);
        assert!(
            sides.iter().any(|(g, _)| *g) && sides.iter().any(|(g, _)| !*g),
            "one of each kind: {:?}",
            cap.circles
        );
        assert_eq!(
            sides[0].1, sides[1].1,
            "a convex cap's two contributions agree on the body's side: {:?}",
            cap.circles
        );
    }

    /// ★★ **The corner the graze exists for.** On a convex cap the lateral's rim-graze and the
    /// cap's seated circle say the same thing, so carrying both changes nothing — the test above
    /// locks that. Here they **disagree**: at a blind bore's ceiling the plate's material is above
    /// the cap while the bore's wall hangs below it, a reflex dihedral in the plane. `edge_mask`'s
    /// `Graze > Seated` precedence then picks the wall's side, which is the whole mechanism of the
    /// repair; before the graze existed, the seated rule flipped the wrong label bit and the next
    /// boolean on that solid came back `CylinderGateUndecided`.
    #[test]
    fn at_a_blind_bores_ceiling_the_graze_and_the_seated_circle_disagree() {
        let mut m = Model::new();
        let plate = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let hole = m.add_cylinder(
            Point3::from_array([5.0, 5.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            3.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // A second operand so the arrangement runs on the bored solid as an operand.
        let boss = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 10.0]),
            Point3::from_array([2.0, 2.0, 12.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, bored, boss).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let tr = trace_on_class_of(
            &m,
            bored,
            boss,
            class_at_z(&planes, 3.0), // the bore's ceiling
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let seated: Vec<bool> = tr
            .circles
            .iter()
            .filter_map(|c| match c.kind {
                SegKind::Seated { body_above } => Some(body_above),
                _ => None,
            })
            .collect();
        let grazes: Vec<bool> = tr
            .circles
            .iter()
            .filter_map(|c| match c.kind {
                SegKind::Graze { body_above } => Some(body_above),
                _ => None,
            })
            .collect();
        assert_eq!(
            seated.len(),
            1,
            "the ceiling is seated here: {:?}",
            tr.circles
        );
        assert_eq!(grazes.len(), 1, "the wall grazes here: {:?}", tr.circles);
        assert_ne!(
            seated[0], grazes[0],
            "a reflex dihedral: the cap and the wall put the body on opposite sides — this is the \
             only place the precedence matters, and the only reason the graze is emitted at all"
        );
        // The wall hangs below the ceiling, and that is the side `edge_mask` believes.
        assert!(!grazes[0], "the bore's wall is below its ceiling");
    }

    /// Partial overlap (E5) is NOT merged — different endpoints mean different edges. a and b share
    /// the y=1 plane; a's y=1 chord is x∈[0,2], b's is x∈[1,3] — overlapping on x∈[1,2] but not
    /// coincident. They must stay separate.
    #[test]
    fn partial_overlap_is_not_merged() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
        // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
        // seated≡transversal coincidence).
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        // The shared y=1 wall class (a face at y=1).
        let y1 = planes
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12))
            .expect("a y=1 face");
        // a's chord x∈[0,2] and b's chord x∈[1,3] ride y=1 but have different endpoints, so they
        // stay as two distinct MergedSegs. A merge that ignored extent would collapse them to one.
        let on_y1 = merged.iter().filter(|e| e.wall == y1).count();
        assert!(
            on_y1 >= 2,
            "a's and b's y=1 chords stay distinct (partial overlap not merged): {on_y1}"
        );
    }

    fn face_on_z1(
        fh: Handle<Face>,
        surf_ix: &HashMap<Handle<Face>, usize>,
        faces: &[FaceRow],
    ) -> bool {
        let p = &faces[surf_ix[&fh]];
        // A cap in the z=1 plane: all three defining points at z=1.
        p.plane()
            .tri
            .iter()
            .all(|q| (q.as_array()[2] - 1.0).abs() < 1e-12)
    }

    /// Axis-aligned cubes never decline: every face is seated, a clean transversal chord, or a
    /// parallel miss. This is falsifiable — a `declined` entry would mean the producer hit a
    /// degeneracy it should not on this input.
    #[test]
    fn axis_aligned_cubes_decline_nothing() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            standard,
            notes,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        for wc in 0..planes.len() {
            let tr = trace_on_class_of(
                &m,
                a,
                b,
                wc,
                &jd,
                &[],
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plane_ix,
                Default::default(),
            );
            assert!(tr.declined.is_empty(), "class {wc} declined: {tr:?}");
        }
    }

    fn tilt_by(m: &mut Model, s: Handle<Solid>, deg: nacre_scalar::Rat) -> Handle<Solid> {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let out = transform(
            m,
            s,
            &Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(deg).unwrap(),
            }),
        )
        .unwrap();
        m.rebuild_adjacency();
        out
    }

    /// **The crossing collector's direction families really are equivalence classes.**
    ///
    /// ★ `split_at_crossings` replaced a predicate per wall pair with "same family?", which is sound
    /// only because `wc ∩ w ∥ wc ∩ r` is transitive. The argument is in that function; this is the
    /// **check**, over every wall pair of every class of a rotated and an axis-aligned fold: the
    /// predicate's own answer must agree with the partition, everywhere.
    ///
    /// It matters that this is a test and not a spike. A violation would not reject or panic — it
    /// would invent a crossing point where the lines never meet, or lose one where they do, and the
    /// census would drift by one vertex somewhere. **A silent wrong answer is exactly what a corpus
    /// this size can hide**, so the invariant is asserted rather than inspected.
    ///
    /// Separate from the timing spike on purpose: this asks the predicate about every pair, which
    /// fills `WitnessPoint`'s realization cells and would make the phase timers read the collector 1.76x
    /// cheaper than it is.
    #[test]
    #[ignore = "slow: every wall pair of every class, two folds"]
    fn direction_families_partition_the_walls() {
        for rotated in [true, false] {
            let n = 24i128;
            let mut m = Model::new();
            let mut acc = if rotated {
                m.add_cuboid(
                    Point3::from_array([-3.0, -3.0, 0.0]),
                    Point3::from_array([3.0, 3.0, 2.0]),
                )
            } else {
                m.add_cuboid(
                    Point3::from_array([-1.0, -1.0, 0.0]),
                    Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
                )
            };
            m.rebuild_adjacency();
            let mut pairs = 0usize;
            for i in 0..n {
                let fin = if rotated {
                    let f = m.add_cuboid(
                        Point3::from_array([2.0, -0.4, 0.0]),
                        Point3::from_array([8.0, 0.4, 1.0]),
                    );
                    m.rebuild_adjacency();
                    tilt_by(&mut m, f, nacre_scalar::Rat::new(360 * i, n).unwrap())
                } else {
                    let x = i as f64 * 0.5;
                    let f = m.add_cuboid(
                        Point3::from_array([x, 0.5, 0.0]),
                        Point3::from_array([x + 0.2, 4.0, 2.0]),
                    );
                    m.rebuild_adjacency();
                    f
                };
                // The audit needs the *same* judging context the boolean uses, and the walls of each
                // class as the tracer finds them — so it re-runs the setup and the trace, then checks
                // the partition the collector would build.
                let setup = plane_index_setup(&m, acc, fin).expect("setup");
                let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
                let trace_in = combinatorics::trace_input(
                    &m,
                    [(acc, &setup.inc_a), (fin, &setup.inc_b)],
                    &setup.surf_ix,
                    setup.planes.len(),
                    &jd,
                    &setup.plane_ix,
                    &setup.cyls,
                    Default::default(),
                );
                for wc in 0..setup.geom.len() {
                    let tr = trace_on_class(
                        &trace_in,
                        wc,
                        &jd,
                        &setup.cyls,
                        &setup.planes,
                        &setup.plane_ix,
                        &Aliases::default(),
                    );
                    let mut walls: Vec<usize> = Vec::new();
                    for s in &tr.segs {
                        if !walls.contains(&s.wall) {
                            walls.push(s.wall);
                        }
                    }
                    let par = |a: usize, b: usize| jd.plane_pair_dir_sign(wc, a, b) == 0;
                    // The partition, exactly as `split_at_crossings` builds it.
                    let mut reps: Vec<usize> = Vec::new();
                    let dir: Vec<usize> = walls
                        .iter()
                        .map(|&w| {
                            reps.iter().position(|&rep| par(rep, w)).unwrap_or_else(|| {
                                reps.push(w);
                                reps.len() - 1
                            })
                        })
                        .collect();
                    for (i, &a) in walls.iter().enumerate() {
                        for (j, &b) in walls.iter().enumerate().skip(i + 1) {
                            pairs += 1;
                            assert_eq!(
                                par(a, b),
                                dir[i] == dir[j],
                                "class {wc}: walls {a} and {b} — the predicate and the partition \
                                 disagree, so parallelism is not transitive here"
                            );
                        }
                    }
                }
                acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
                    .expect("fuse")
                    .0[0];
                m.rebuild_adjacency();
            }
            println!(
                "  {} fold: {pairs} wall pairs audited, all agree",
                if rotated { "rotated" } else { "axis-aligned" }
            );
            assert!(
                pairs > 10_000,
                "{} fold audited only {pairs} pairs — too few to mean anything",
                if rotated { "rotated" } else { "axis-aligned" }
            );
        }
    }

    /// **S2a: where does a whole boolean's time go?** Every earlier profile answered a share *of a
    /// phase* — and the one that mattered was never taken.
    ///
    /// ★ Two ways to be wrong that this is built against:
    ///
    /// 1. **Mixed clocks.** The last breakdown summed per-thread timers (CPU) and read the result
    ///    against the fold's wall clock. The part came out larger than the whole — 3,719ms of a
    ///    2,540ms fold — and the ratio it implied was meaningless. **Run this with
    ///    `--no-default-features`**, where every phase and the total are the same clock.
    /// 2. **Phases that do not add up.** The sum is printed against the measured whole, so anything
    ///    unaccounted for shows as a gap rather than hiding inside a phase's share.
    /// 3. ★★★ **Another test adding to the same counters.** `phase::` are process-global atomics
    ///    and `cargo test` runs this binary's tests on parallel threads, so `reset()` here does not
    ///    fence anything: every concurrently running test that calls `boolean` — and in an
    ///    `--ignored` pass that is the *other* spikes and the two rotation stress tests — lands in
    ///    the same buckets, while `whole` below is this thread's wall clock alone. The share
    ///    percentages then exceed 100 for reasons that have nothing to do with the code.
    ///    **Run it alone**: `cargo test -p nacre-ops --no-default-features spike_where -- --ignored
    ///    --nocapture --test-threads=1`. (Found 2026-08-22 while auditing an unrelated change; the
    ///    numbers already printed were taken that way, so they stand — but nothing said so.)
    ///
    /// Timers live on the production path (`phase::` in this module), not in a replica of it.
    #[test]
    #[ignore = "spike"]
    fn spike_where_the_boolean_spends_it() {
        for rotated in [true, false] {
            spend(60, rotated);
        }
    }

    /// One fold's phase breakdown. `rotated` picks which predicate routes the judgements take:
    /// an axis-aligned fold answers on the exact path, a rotated one mostly on the certified one,
    /// and the difference between the two breakdowns is what the certification actually costs.
    fn spend(n: i128, rotated: bool) {
        let mut m = Model::new();
        // ★ **The two folds must be the same *shape* of work, not the same model minus a rotation.**
        // Dropping the tilt stacks all 60 fins on one another — 7.6x fewer segments, a degenerate
        // model, and a comparison that says nothing. The axis-aligned arm places each fin at its
        // own x instead, so both arms fuse `n` distinct blades onto a growing solid and the only
        // difference is which predicate route the judgements take.
        let mut acc = if rotated {
            m.add_cuboid(
                Point3::from_array([-3.0, -3.0, 0.0]),
                Point3::from_array([3.0, 3.0, 2.0]),
            )
        } else {
            m.add_cuboid(
                Point3::from_array([-1.0, -1.0, 0.0]),
                Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
            )
        };
        m.rebuild_adjacency();
        let blade = |m: &mut Model, i: i128| {
            let f = if rotated {
                m.add_cuboid(
                    Point3::from_array([2.0, -0.4, 0.0]),
                    Point3::from_array([8.0, 0.4, 1.0]),
                )
            } else {
                let x = i as f64 * 0.5;
                m.add_cuboid(
                    Point3::from_array([x, 0.5, 0.0]),
                    Point3::from_array([x + 0.2, 4.0, 2.0]),
                )
            };
            m.rebuild_adjacency();
            if rotated {
                tilt_by(m, f, nacre_scalar::Rat::new(360 * i, n).unwrap())
            } else {
                f
            }
        };
        // Warm the code paths, then zero the counters: the first boolean pays for lazily-built
        // caches that the other n do not.
        {
            let f = blade(&mut m, if rotated { 1 } else { n });
            acc = super::boolean(&mut m, BoolKind::Fuse, acc, f)
                .expect("warm")
                .0[0];
            m.rebuild_adjacency();
        }
        phase::reset();
        phase::scale::reset();

        let mut whole = std::time::Duration::ZERO;
        for i in 0..n {
            let fin = blade(&mut m, i);
            let t = std::time::Instant::now();
            acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
                .expect("fuse")
                .0[0];
            whole += t.elapsed();
            m.rebuild_adjacency();
        }

        let rows = phase::all();
        // ★ Indentation is nesting: depth 0 are `boolean`'s own phases, depth 1 the inside of
        // `trace_result_faces`, depth 2 the inside of `split_at_crossings`. Depth 1 partitions the
        // work depth 0 does not name, so 0 and 1 add to the whole — and adding depth 2 on top would
        // count `split_at_crossings` twice.
        let depth = |l: &str| (l.len() - l.trim_start().len()) / 2;
        let sum: u64 = rows
            .iter()
            .filter(|(l, _)| depth(l) <= 1)
            .map(|(_, ns)| ns)
            .sum();
        let total = whole.as_nanos() as u64;
        println!(
            "\n{n}-fin fold ({}), whole boolean, serial build:",
            if rotated { "rotated" } else { "axis-aligned" }
        );
        for (label, ns) in &rows {
            println!(
                "  {label:<32} {:>8.1?}  {:>5.1}%",
                std::time::Duration::from_nanos(*ns),
                100.0 * *ns as f64 / total as f64
            );
        }
        println!(
            "  {:<32} {:>8.1?}  {:>5.1}%   ← accounted",
            "sum of the above",
            std::time::Duration::from_nanos(sum),
            100.0 * sum as f64 / total as f64
        );
        println!("  {:<32} {whole:>8.1?}  100.0%   ← measured", "the fold");
        println!(
            "\n  ★ split_at_crossings is {:.1}% of the whole boolean. The plan continues at 25%.",
            100.0
                * rows
                    .iter()
                    .find(|(l, _)| l.trim() == "split_at_crossings")
                    .map(|(_, ns)| *ns)
                    .unwrap_or(0) as f64
                / total as f64
        );

        // ★ **Scale, which decides whether the structural fixes are worth building.** A hull test
        // over wall *pairs* replaces `2 × |segs on r|` predicate calls with `2` — worth nothing when
        // a wall carries one segment, and worth the ratio when it carries many.
        use phase::scale as sc;
        let (segs, walls) = (sc::get(&sc::SEGS), sc::get(&sc::WALLS));
        println!("\n  scale, summed over every class of every boolean:");
        println!(
            "    segments / walls        {segs:>10} / {walls:<10} = {:.2}   ← S3b lives or dies here",
            segs as f64 / walls.max(1) as f64
        );
        println!("    split points (reps)     {:>10}", sc::get(&sc::PTS));
        println!(
            "    (1) collect trips       {:>10}",
            sc::get(&sc::COLLECT_TRIPS)
        );
        println!(
            "    (3) cover trips         {:>10}   = {:.2}x the collecting loop",
            sc::get(&sc::COVER_TRIPS),
            sc::get(&sc::COVER_TRIPS) as f64 / sc::get(&sc::COLLECT_TRIPS).max(1) as f64
        );
    }

    /// S3a gate: **what did hoisting the loop invariant do to the evidence?** `plane_pair_dir_sign`
    /// records on the rotated path, and `BoolReport::coincidences` is a *count*, so asking the same
    /// question fewer times moves it. The answer must not move; the count may.
    #[test]
    #[ignore = "spike"]
    fn spike_report_after_hoisting() {
        let n = 24i128;
        let mut m = Model::new();
        let mut acc = m.add_cuboid(
            Point3::from_array([-3.0, -3.0, 0.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        let (mut total, mut loosest) = (0usize, String::new());
        for i in 0..n {
            let fin = m.add_cuboid(
                Point3::from_array([2.0, -0.4, 0.0]),
                Point3::from_array([8.0, 0.4, 1.0]),
            );
            m.rebuild_adjacency();
            let fin = tilt_by(&mut m, fin, nacre_scalar::Rat::new(360 * i, n).unwrap());
            let (solids, report) =
                crate::boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
            total += report.coincidences;
            if let Some(e) = &report.loosest {
                loosest = format!("{e:?}");
            }
            acc = solids[0];
            m.rebuild_adjacency();
        }
        let v = nacre_props::mass_props(&m, acc).unwrap().volume;
        println!("\n  coincidences over {n} booleans : {total}");
        println!("  loosest (last)                : {loosest}");
        println!("  volume                        : {v:.9}");
    }

    /// B0: what is actually left to save, **in production's configuration**?
    ///
    /// The earlier phase timing used `ClassReuse::Off`, so it counted work production never does —
    /// `reuse.rs` skips most classes outright. This counts what survives that, how much of it a
    /// per-class bounding box would cull, and how often the two things that would make the cull
    /// unsound actually occur. Counting only; nothing is built and nothing is timed (the machine is
    /// not quiet).
    #[test]
    #[ignore = "spike"]
    fn spike_cull_potential() {
        let n = 60i128;
        let mut m = Model::new();
        let mut acc = m.add_cuboid(
            Point3::from_array([-3.0, -3.0, 0.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        let mut tot = Counts::default();
        for i in 0..n {
            let fin = m.add_cuboid(
                Point3::from_array([2.0, -0.4, 0.0]),
                Point3::from_array([8.0, 0.4, 1.0]),
            );
            m.rebuild_adjacency();
            let fin = tilt_by(&mut m, fin, nacre_scalar::Rat::new(360 * i, n).unwrap());
            tot.add(&count_cull(&m, acc, fin));
            acc = super::boolean(&mut m, BoolKind::Fuse, acc, fin)
                .expect("fuse")
                .0[0];
            m.rebuild_adjacency();
        }
        println!("over the fold, production config (ClassReuse::Proved):");
        println!("  classes                 {:>9}", tot.classes);
        println!(
            "  arranged after reuse    {:>9}  ({:.1}%)",
            tot.arranged,
            100.0 * tot.arranged as f64 / tot.classes as f64
        );
        println!(
            "  (face, class) pairs     {:>9}   in arranged classes",
            tot.pairs
        );
        println!(
            "  ★ cullable by box       {:>9}  ({:.1}%)",
            tot.cullable,
            100.0 * tot.cullable as f64 / tot.pairs.max(1) as f64
        );
        println!(
            "  seated footprint is one piece: {} of {} arranged classes",
            tot.one_piece, tot.arranged
        );
        // ★ Recorded, and **not** the operative number. The plan restricted the cull to classes
        // whose seated footprint is one connected piece, fearing that several pieces would need
        // several seeds. They do not: what the correction needs is that the culled chords' parity
        // be *constant over the footprint*, and the cull criterion (`box(F)` disjoint from the
        // footprint box) already guarantees no culled chord enters that box. A box is connected, so
        // the parity is one constant however many pieces the seated faces form.
        println!(
            "  (one-piece classes only:  {:>6}  of {} — not the limit, see the note)",
            tot.cullable_1p, tot.pairs_1p
        );
        println!(
            "  (class, solid) needing a seed: {}  of {}   — of those, {} are free (solid box misses)",
            tot.needs_seed,
            2 * tot.arranged,
            tot.seed_free
        );
    }

    #[derive(Default)]
    struct Counts {
        classes: usize,
        arranged: usize,
        pairs: usize,
        cullable: usize,
        one_piece: usize,
        needs_seed: usize,
        /// Pairs in classes whose footprint is one piece — the realistic cull, since the
        /// multi-piece ones cannot take a single seed.
        pairs_1p: usize,
        cullable_1p: usize,
        /// Of `needs_seed`, the ones where the solid's whole box misses the footprint, so it
        /// **cannot** enclose it and the seed is false without any work.
        seed_free: usize,
    }

    impl Counts {
        fn add(&mut self, o: &Counts) {
            self.classes += o.classes;
            self.arranged += o.arranged;
            self.pairs += o.pairs;
            self.cullable += o.cullable;
            self.one_piece += o.one_piece;
            self.needs_seed += o.needs_seed;
            self.pairs_1p += o.pairs_1p;
            self.cullable_1p += o.cullable_1p;
            self.seed_free += o.seed_free;
        }
    }

    fn count_cull(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> Counts {
        let mut c = Counts::default();
        let Ok(setup) = plane_index_setup(m, a, b) else {
            return c;
        };
        let PlaneSetup {
            planes: faces_tab,
            n_a,
            geom,
            plane_ix,
            class_owner,
            standard,
            notes,
            ..
        } = setup;
        let jd = Judge::new(&geom, standard, &notes);
        let _ = &jd;
        let plans = crate::reuse::class_plans(
            m,
            crate::reuse::ClassReuse::Proved,
            BoolKind::Fuse,
            a,
            b,
            &geom,
            &class_owner,
        );
        let boxes: Vec<[[f64; 2]; 3]> = faces_tab
            .iter()
            .map(|f| face_box(m, f.face().expect("real")))
            .collect();
        let overlap = |x: &[[f64; 2]; 3], y: &[[f64; 2]; 3]| {
            (0..3).all(|k| x[k][0] <= y[k][1] && y[k][0] <= x[k][1])
        };
        c.classes = geom.len();
        for (wc, plan) in plans.iter().enumerate() {
            if *plan != crate::reuse::ClassPlan::Arrange {
                continue;
            }
            c.arranged += 1;
            let seated: Vec<usize> = (0..faces_tab.len())
                .filter(|&f| plane_ix[f].plane() == wc)
                .collect();
            if seated.is_empty() {
                continue;
            }
            // The class's region of interest, and whether it is one connected piece (box graph).
            let mut foot = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
            for &s in &seated {
                for k in 0..3 {
                    foot[k][0] = foot[k][0].min(boxes[s][k][0]);
                    foot[k][1] = foot[k][1].max(boxes[s][k][1]);
                }
            }
            let mut parent: Vec<usize> = (0..seated.len()).collect();
            for i in 0..seated.len() {
                for j in (i + 1)..seated.len() {
                    if overlap(&boxes[seated[i]], &boxes[seated[j]]) {
                        let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                        parent[ri] = rj;
                    }
                }
            }
            let pieces = (0..seated.len())
                .map(|i| uf_find(&mut parent, i))
                .collect::<std::collections::HashSet<_>>()
                .len();
            if pieces == 1 {
                c.one_piece += 1;
            }
            // How many of this class's traced faces the box would cull, and whether a solid that
            // only crosses `W` has any culled face (that is the one needing a seed).
            let mut culled_side = [false; 2];
            let mut seated_side = [false; 2];
            for &s in &seated {
                seated_side[usize::from(s >= n_a)] = true;
            }
            for f in 0..faces_tab.len() {
                c.pairs += 1;
                if pieces == 1 {
                    c.pairs_1p += 1;
                }
                if plane_ix[f].plane() == wc {
                    continue; // seated: never culled
                }
                if !overlap(&boxes[f], &foot) {
                    c.cullable += 1;
                    if pieces == 1 {
                        c.cullable_1p += 1;
                    }
                    culled_side[usize::from(f >= n_a)] = true;
                }
            }
            // Each operand's whole box: if it misses the footprint entirely it cannot enclose it,
            // so its parity is false for free.
            let mut solid_box = [[[f64::INFINITY, f64::NEG_INFINITY]; 3]; 2];
            for (f, fb) in boxes.iter().enumerate() {
                let side = usize::from(f >= n_a);
                for k in 0..3 {
                    solid_box[side][k][0] = solid_box[side][k][0].min(fb[k][0]);
                    solid_box[side][k][1] = solid_box[side][k][1].max(fb[k][1]);
                }
            }
            for side in 0..2 {
                if culled_side[side] && !seated_side[side] {
                    c.needs_seed += 1;
                    if !overlap(&solid_box[side], &foot) {
                        c.seed_free += 1;
                    }
                }
            }
        }
        c
    }
}
