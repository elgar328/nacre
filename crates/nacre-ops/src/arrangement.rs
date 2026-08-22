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
use crate::combinatorics::{NodeId, three_plane_name};
use crate::planes::*;
use crate::tolerant::{ImplicitPoint, Judge};
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
    pub end_h: [usize; 2],
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
/// *circles* rest on a rule about *walls*. [`arc_split_witness`] asks it of the segments
/// themselves now, so a crossed circle is **split into arcs** rather than treated as the
/// closed cell it is not, and the wall rule is free to become precise about its own question.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CircleTrace {
    /// The cylinder class ([`ClassIx::Cyl`] payload) whose circle this is.
    pub cyl: usize,
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// One circle of the class's arrangement after merging: both operands' contributions on one
/// cylinder class, plus the cylinder's exact statement (the containment predicates read it).
/// The circle twin of [`MergedSeg`].
#[derive(Clone, Debug)]
pub(crate) struct MergedCircle {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    pub merged: Vec<(SolidSide, SegKind)>,
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
    /// Inherited whole from the circle: an arc is a piece of the same trace, so it carries the same
    /// `(solid, kind)` contributions and `edge_mask` reads it unchanged.
    ///
    /// ★ Read by `label_cells`' `mask_of`: crossing an arc flips the same bits crossing its circle
    /// would, which is what "a piece of the same trace" means. (It was carried before that consumer
    /// existed, because the split is the only place that knows which circle an arc came from.)
    pub merged: Vec<(SolidSide, SegKind)>,
}

/// Group a class's circle traces by cylinder class, in ascending class order (deterministic —
/// replay mints handles from this order). The def comes from the class table, which a
/// `ClassIx::Cyl` index indexes directly.
fn merge_circles(
    circles: &[CircleTrace],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Vec<MergedCircle>, BoolError> {
    let mut out: Vec<MergedCircle> = Vec::new();
    let mut sorted: Vec<&CircleTrace> = circles.iter().collect();
    sorted.sort_by_key(|c| c.cyl);
    for c in sorted {
        if let Some(last) = out.last_mut() {
            if last.cyl == c.cyl {
                last.merged.push((c.solid, c.kind));
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
            merged: vec![(c.solid, c.kind)],
        });
    }
    Ok(out)
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// The closed circle elements beside the segments (M6-2a) — see [`CircleTrace`].
    pub circles: Vec<CircleTrace>,
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
        let mut names = Vec::new();
        for i in 0..s.len() {
            for j in (i + 1)..s.len() {
                for k in (j + 1)..s.len() {
                    let t = [s[i], s[j], s[k]];
                    if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                        // Shares a line: names no point, and tells us two walls are one line.
                        self.union_wall(t[0], t[1], t[2]);
                        self.union_wall(t[1], t[0], t[2]);
                        self.union_wall(t[2], t[0], t[1]);
                    } else {
                        names.push(NodeId::three_planes(t));
                    }
                }
            }
        }
        if let Some((&rep, rest)) = names.split_first() {
            for &t in rest {
                self.union_point(rep, t);
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

    fn union_point(&mut self, a: NodeId, b: NodeId) {
        let (ra, rb) = (self.find_point(a), self.find_point(b));
        if ra == rb {
            return;
        }
        // Smallest name wins, so the representative is a function of the class alone.
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
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

/// Why a ring of names could not become a ring of plane data.
///
/// ★ **It is not one label.** The two causes are different propositions — a triple that came out
/// degenerate, versus a vertex that has no triple at all — and the four callers of [`plane_ring`]
/// already map the old single failure onto *three* different [`DeclineKind`]s, so the label is the
/// caller's to choose and this only says which cause it was.
enum RingFail {
    /// Two of a vertex's three planes coincide, so the name denotes no point.
    Collapsed,
    /// A branch point: named exactly, but not by three planes.
    Branch,
}

/// The decline kind for a [`RingFail`] at a caller that has no finer fact to add.
///
/// ★ Two of the four callers *do* have one — they know which loop of the face failed and report
/// `OuterRing`/`HoleRing` — so this is not "the mapping", only the default.
fn decline_of(f: RingFail) -> DeclineKind {
    match f {
        RingFail::Collapsed => DeclineKind::CollapsedTriple,
        RingFail::Branch => DeclineKind::BranchNode,
    }
}

/// ★ **The one place a ring of *names* becomes a ring of *plane data*** for the tracer, which asks
/// only "which side, which wall, which third plane" of each vertex. The sort that used to stand
/// here is gone — [`NodeId`]'s only constructor sorts — so what remains is the collapse check, and
/// that is the whole reason this function exists.
fn plane_ring(ts: &[NodeId]) -> Result<Vec<[usize; 3]>, RingFail> {
    ts.iter()
        .map(|&n| {
            let c = three_plane_name(n).ok_or(RingFail::Branch)?;
            // Sorted by construction, so equal neighbours catch every duplicate.
            if c[0] == c[1] || c[1] == c[2] {
                return Err(RingFail::Collapsed);
            }
            Ok(c)
        })
        .collect()
}

/// One feature node on the line `L = W ∩ fp`: a point named by its third plane `r` (raw index).
struct Node {
    /// Third plane naming this point on `L` (the point is `{W, fp, r}`).
    r: usize,
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
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    out: &mut Trace,
) {
    // `fp` names a *face* (`n_out`, `orient`, the declined log); `fc` names the *plane class* it
    // lies on (triples, comparisons, predicate arguments). Every ring below is in class form, so
    // the two must not be confused — see `canon_ring`.
    let fc = plane_ix[fp].plane();
    // Each ring pairs its class-form triples with the **carried walls** the producer read off
    // the model's edges (`NamedRing`) — the Crossing arm below names an edge's wall from the
    // ride, never from the two endpoint names.
    let outer = match loops.outer.as_ref().and_then(|r| r.poly()) {
        Some(nr) => match plane_ring(&nr.triples) {
            Ok(ts) => (ts, nr.walls.as_slice()),
            Err(f) => {
                out.declined.push((fp, decline_of(f)));
                return;
            }
        },
        // ★ A **circular outer is a disk face, and it misses this class by proof** — the same
        // clearance the circular-hole arm below rests on. Two classes can carry `W`: another ⊥
        // class, which is parallel to the disk's own and meets it nowhere; or a ∥-axis wall,
        // which the population gate proved farther from the axis than r, so the disk (radius r
        // about that axis) does not reach the line `W ∩ fc`. Either way the face contributes no
        // chord — a miss, like a parallel plane, not a decline.
        //
        // `None` is also how an *unnamed* ring arrives, but not for a face the tracer reaches:
        // the only loops `loop_triples` leaves unnamed are the ones it rejects for, and those
        // are `Err` at the source. A circle is the one `None` that means "named, and not a
        // polygon".
        None if matches!(loops.outer, Some(combinatorics::LoopRing::Circle { .. })) => return,
        None => {
            out.declined.push((fp, DeclineKind::OuterRing));
            return;
        }
    };
    let mut holes: Vec<(Vec<[usize; 3]>, &[usize])> = Vec::new();
    // A hole whose ring cannot be named is not "no hole" — swallowing the error would trace the
    // face as solid where it is pierced, which is a silent wrong answer rather than a reject.
    let Some(raw_holes) = &loops.holes else {
        out.declined.push((fp, DeclineKind::HoleRing));
        return;
    };
    for r in raw_holes {
        match r {
            // ★ A **circular hole is skipped, by proof, not by hope**: the population gate
            // established every ∥-axis wall clear of the circle's cylinder by more than r, and
            // this trace's line `L = W ∩ fc` lies in such a wall — so the circle cannot meet
            // `L`, and a disjoint hole contributes no crossings to the 3-valued scan.
            combinatorics::LoopRing::Circle { .. } => continue,
            combinatorics::LoopRing::Poly(r) => match plane_ring(&r.triples) {
                Ok(ts) => holes.push((ts, r.walls.as_slice())),
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
    let third_on_l = |t: [usize; 3], out: &mut Trace| -> Result<usize, DeclineKind> {
        let (mut offs, mut has_fp, mut has_w) = (Vec::<usize>::new(), false, false);
        for &x in &t {
            if x == fc {
                has_fp = true;
            } else if x == wc {
                has_w = true;
            } else {
                offs.push(x);
            }
        }
        match (has_fp, has_w, offs.len()) {
            // The ordinary point: `t` carries both of this line's planes, and the third both pins
            // and names it.
            (true, true, 1) => Ok(offs[0]),
            // On `wc` (this is a run vertex) yet `wc` does not name it. `t`'s classes are distinct
            // (`plane_ring`) and `fc != wc` here, so `wc` is a **fourth** plane through the point.
            (_, false, _) => {
                let mut set = t.to_vec();
                set.push(wc);
                set.sort_unstable();
                set.dedup();
                out.aliases.record(jd, &set);
                // ★ A handle has a duty the identity does not: it must **cut** `L`. One parallel
                // to it names no point there, and `order_along` — being `orient3d × dir_sign` —
                // would read 0 against everything, fabricating a coincidence rather than missing
                // one.
                offs.sort_unstable();
                offs.iter()
                    .copied()
                    .find(|&r| jd.plane_pair_dir_sign(wc, fc, r) != 0)
                    .ok_or(DeclineKind::FourPlane)
            }
            _ => Err(DeclineKind::RunName),
        }
    };

    // Phase A — scan the outer ring and every hole ring, collecting feature nodes into ONE list.
    // A hole ring keeps its stored CW winding, but the Phase-A predicates are winding-agnostic; the
    // parity sweep (Phase C) then carves the hole because its two crossings of `L` toggle parity
    // back to void between them. A ring that does not meet `W` (every vertex one side) yields no
    // nodes; `run_counter` is shared so run ids stay unique across rings.
    let mut nodes: Vec<Node> = Vec::new();
    let mut run_counter = 0usize;
    let mut declined: Option<DeclineKind> = None;
    'rings: for (ring, walls) in std::iter::once(&outer).chain(holes.iter()) {
        let n = ring.len();
        // Where this ring meets `L`, and whether it crosses or only touches — the walk the ray
        // caster shares (`combinatorics::ring_against_plane`). What is done with a feature is this
        // function's own business: naming it, recording a four-plane alias, deciding occupancy.
        let Some(features) = combinatorics::ring_against_plane(jd, ring, wc) else {
            // Every vertex on `W`: a ring lying in the cut plane is degenerate here.
            declined = Some(DeclineKind::AllOnPlane);
            break;
        };
        for feature in features {
            match feature {
                combinatorics::Feature::Crossing { edge } => {
                    // The crossed edge's wall, **carried** from the producer — it used to be
                    // re-derived from the two endpoint names (`ring_from_names`), which is
                    // sound only while every vertex lies on exactly three planes and could
                    // hand back a plane the edge does not ride at a concurrency.
                    nodes.push(Node {
                        r: walls[edge],
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
                } => {
                    let name = |k: usize, out: &mut Trace| third_on_l(ring[(first + k) % n], out);
                    if m == 1 {
                        match name(0, out) {
                            Ok(r) => nodes.push(Node {
                                r,
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
                        let names: Result<Vec<usize>, DeclineKind> =
                            (0..m).map(|k| name(k, out)).collect();
                        match names {
                            Ok(rs) => {
                                // Every run is one-sided: it is an *edge* of `f` lying on `L`, so
                                // `f` is on one side of it whatever the ring does afterwards.
                                // `flanks_differ` says only whether the sweep's parity toggles here
                                // (Phase B gives it to `flip`) — it is not an occupancy fact, and
                                // gating the side on it left a crossing run classified as a
                                // straddling transversal.
                                let graze_above = Some(run_body_above(jd, faces, wc, fc, fp, &rs));
                                for r in rs {
                                    nodes.push(Node {
                                        r,
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
    nodes.sort_by(
        |a, b| match combinatorics::order_along(jd, wc, fc, a.r, b.r) {
            -1 => std::cmp::Ordering::Less,
            1 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        },
    );
    for w in nodes.windows(2) {
        if combinatorics::order_along(jd, wc, fc, w[0].r, w[1].r) == 0 {
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
    let mut seg_start: Option<(usize, Option<bool>)> = None;
    let emit = |a: usize, b: usize, graze: Option<bool>, out: &mut Trace| {
        out.segs.push(Seg {
            wall: plane_ix[fp].plane(),
            end: [
                NodeId::three_planes([wc, fc, a]),
                NodeId::three_planes([wc, fc, b]),
            ],
            end_h: [a, b],
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
            out.touches.push(NodeId::three_planes([wc, fc, nodes[k].r]));
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
            (None, true) => seg_start = Some((nodes[k].r, gap_graze)),
            // The kind changes here, so close and reopen: a segment must be homogeneous, because
            // `split_at_crossings` subdivides it later and every piece inherits its kind.
            (Some((a, graze)), true) if graze != gap_graze => {
                emit(a, nodes[k].r, graze, out);
                seg_start = Some((nodes[k].r, gap_graze));
            }
            (Some((a, graze)), false) => {
                emit(a, nodes[k].r, graze, out);
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
fn run_body_above(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    wc: usize,
    fc: usize,
    fp: usize,
    rs: &[usize],
) -> bool {
    let planes = jd.planes;
    let t = combinatorics::order_along(jd, wc, fc, rs[0], rs[rs.len() - 1]);
    let sigma = faces[fp].plane().n_out.dot(planes[fc].plane.normal());
    (t < 0) == (sigma > 0.0)
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
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    out: &mut Trace,
) {
    let planes = jd.planes;
    let w_normal = planes[wc].plane.normal();
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
            match circle_on_class(cf, &planes[wc]) {
                Ok(Some(on)) => out.circles.push(CircleTrace {
                    cyl: k,
                    solid: which,
                    kind: match on {
                        CylOnClass::Crosses => SegKind::Transversal {
                            mat: cf.orient_sign,
                        },
                        CylOnClass::Grazes { body_above } => SegKind::Graze { body_above },
                    },
                }),
                Ok(None) => {}
                Err(kind) => out.declined.push((fp, kind)),
            }
            continue;
        }
        if plane_ix[fp].plane() != wc {
            trace_transversal_face(fp, fl, which, wc, jd, faces, plane_ix, out);
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
        let mut rings = Vec::new();
        match &fl.outer {
            Some(combinatorics::LoopRing::Poly(nr)) => match plane_ring(&nr.triples) {
                Ok(ts) => rings.push(ts),
                // ★ `OuterRing` here, not `CollapsedTriple`: this caller reports *which loop*
                // failed, which is the finer fact when a face has several. A branch node is a
                // different cause and keeps its own name either way.
                Err(RingFail::Collapsed) => {
                    out.declined.push((fp, DeclineKind::OuterRing));
                    continue;
                }
                Err(RingFail::Branch) => {
                    out.declined.push((fp, DeclineKind::BranchNode));
                    continue;
                }
            },
            Some(combinatorics::LoopRing::Circle { cyl }) => {
                out.circles.push(CircleTrace {
                    cyl: *cyl,
                    solid: which,
                    kind,
                });
            }
            None => {
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
                }),
                combinatorics::LoopRing::Poly(nr) => match plane_ring(&nr.triples) {
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
        let mut emit_ring = |tris: &[[usize; 3]]| {
            let n = tris.len();
            for i in 0..n {
                // Edge i runs vertex i → vertex i+1; the wall it rides is the plane the two
                // endpoint triples share besides `fc`. The triples are already dense plane ids
                // (`face_vertex_triples` mapped them through `plane_ix`), so this is a plain set
                // intersection — no second remap, which under a non-idempotent `plane_ix` would
                // index the table with a value that is already an index.
                let (t0, t1) = (tris[i], tris[(i + 1) % n]);
                let mut shared: Vec<usize> = t0
                    .iter()
                    .copied()
                    .filter(|&c| c != fc && t1.contains(&c))
                    .collect();
                shared.sort_unstable();
                shared.dedup();
                let [wall] = shared[..] else {
                    out.declined.push((fp, DeclineKind::SeatedEdgeNaming));
                    continue;
                };
                // The handle on `wc ∩ wall` is what the triple carries besides those two. A
                // seated face lies in `wc`, so `fc == wc` here and `shared` being a singleton
                // is the same fact as there being exactly one such element.
                let handle = |t: [usize; 3]| {
                    t.into_iter()
                        .find(|&c| c != fc && c != wall)
                        .expect("a seated edge's endpoint has a third plane")
                };
                out.segs.push(Seg {
                    wall,
                    end: [NodeId::three_planes(t0), NodeId::three_planes(t1)],
                    end_h: [handle(t0), handle(t1)],
                    solid: which,
                    kind,
                });
            }
        };
        for ring in &rings {
            emit_ring(ring);
        }
    }
}

/// Whether a lateral face crosses the plane class `wp` in a full circle (M6-2a).
///
/// `Ok(Some(()))` — the class is ⊥ this cylinder's axis and its axis parameter lies strictly
/// inside the face's rim span. `Ok(None)` — no circle (a non-⊥ class carries none, and a ⊥
/// class outside the span never meets this face; a rim-coincident plane would have been one
/// class with the rim's cap and traced seated). `Err` — the face cannot answer (no rim span,
/// or the class has no exact description the gate would already have refused); declining
/// beats a silently missing circle, which would corrupt every label on the class.
/// What a lateral face leaves on a ⊥ plane class — see [`circle_on_class`].
pub(crate) enum CylOnClass {
    /// The class cuts the face in two: the solid straddles it here.
    Crosses,
    /// The class carries one of the face's **rims**: the wall touches it along that circle with
    /// its body to one side. `body_above` is that side, about the class's **stored** normal —
    /// the frame every cell label is written in.
    Grazes { body_above: bool },
}

fn circle_on_class(
    cf: &crate::planes::CylFaceInfo,
    wp: &WorkingPlane,
) -> Result<Option<CylOnClass>, DeclineKind> {
    let Some(coeffs) = wp.base_rat else {
        return Err(DeclineKind::CylSpan);
    };
    if wp.rotated {
        return Err(DeclineKind::CylSpan);
    }
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let m = cf.def.dir();
    // ⊥ to the axis, decided **totally** (`parallel_rat` clears denominators and answers in
    // integers, so this cannot decline for want of bits). A non-⊥ class carries no circle at
    // all — a miss, like a parallel plane, not a decline.
    if !nacre_scalar::parallel_rat(&n, &m) {
        return Ok(None);
    }
    let Some(span) = cf.span else {
        return Err(DeclineKind::CylSpan);
    };
    let Some(t) = crate::planes::axis_param_of_plane(&coeffs, &cf.def) else {
        return Err(DeclineKind::CylSpan);
    };
    if span[0] < t && t < span[1] {
        return Ok(Some(CylOnClass::Crosses));
    }
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
    let up = crate::planes::plus_t_is_above(wp, &cf.def);
    if t == span[0] {
        return Ok(Some(CylOnClass::Grazes { body_above: up }));
    }
    if t == span[1] {
        return Ok(Some(CylOnClass::Grazes { body_above: !up }));
    }
    Ok(None)
}

/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
fn trace_on_class(
    input: &combinatorics::TraceInput,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
) -> Trace {
    let mut out = Trace::default();
    for (side, which) in [SolidSide::A, SolidSide::B].into_iter().enumerate() {
        trace_one(&input.faces[side], which, wc, jd, faces, plane_ix, &mut out);
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
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &combinatorics::EdgeFaces,
    inc_b: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
) -> Trace {
    let input = combinatorics::trace_input(
        model,
        [(a, inc_a), (b, inc_b)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
    );
    trace_on_class(&input, wc, jd, faces, plane_ix)
}

/// One solid only — the second operand slot is filled with the same solid, whose loops are
/// identical, and only `faces[0]` is read.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn trace_one_of(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
    out: &mut Trace,
) {
    let input = combinatorics::trace_input(
        model,
        [(solid, inc), (solid, inc)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
    );
    trace_one(&input.faces[0], which, wc, jd, faces, plane_ix, out);
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
    /// The travel sense from `end[0]` to `end[1]`, when the endpoints can no longer supply it —
    /// see [`combinatorics::Carrier::Plane`]. `None` everywhere except a sub-segment the arc split
    /// cut, whose cut end is a branch point with no third plane to order by.
    pub sense: Option<i8>,
}

/// Merge segments that are the **same geometric edge** — same `wall` and same endpoint-triple set
/// (direction-independent) — into one `MergedSeg`, collecting their contributions. Partial overlap
/// (same `wall`, *different* extent — the E5 case) is left for the per-wall interval overlay in
/// `split_at_crossings` to resolve into non-overlapping sub-segments with unioned contributions.
fn merge_coincident(segs: &[Seg], wc: usize, aliases: &Aliases) -> Vec<MergedSeg> {
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
                MergedSeg {
                    // The canonical name of the line, so every producer on it agrees. The handles
                    // stay as recorded: an aliased wall carries the *same* line, so a handle that
                    // pinned an endpoint there still pins it here.
                    wall: aliases.canon_wall(wc, s.wall),
                    end: [aliases.canon_point(s.end[0]), aliases.canon_point(s.end[1])],
                    end_h: s.end_h.map(combinatorics::EndPin::Class),
                    merged: Vec::new(),
                    sense: None,
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

fn split_at_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    segs: &[MergedSeg],
    aliases: &mut Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    // ★ **The direction sign of each endpoint, once per segment.** `order_along(wc, wall, i, j)`
    // factors into `orient3d(wc, wall, i, j) × dir_sign(wc, wall, j)`, and the second factor does
    // **not** mention `i` — for a segment's endpoints it is a property of the segment alone. The
    // containment test below sweeps `i` over every wall, so leaving it inside asked the same
    // question `|walls|` times over. (Measured: 1.5M `dir_sign` calls where 113k are distinct.)
    // ★ **This pass is plane-only, and says so once.** It runs *before* the arc split, so every
    // endpoint is still pinned by a third plane; a cylinder-pinned one here would be a wiring
    // failure, and `RingNaming` is already the sentence for it ("a node has no third plane to be
    // named by"). Reading the classes up front keeps the body below in the one vocabulary it can
    // answer in, rather than threading an `Option` it has no arm for.
    let end_c: Vec<[usize; 2]> = segs
        .iter()
        .map(|s| {
            let [a, b] = s.end_h;
            match (a.class(), b.class()) {
                (Some(a), Some(b)) => Ok([a, b]),
                _ => Err(reject(RejectReason::RingNaming)),
            }
        })
        .collect::<Result<_, BoolError>>()?;
    let end_ds: Vec<[i8; 2]> = segs
        .iter()
        .zip(&end_c)
        .map(|(s, c)| {
            [
                combinatorics::dir_sign(jd, wc, s.wall, c[0]),
                combinatorics::dir_sign(jd, wc, s.wall, c[1]),
            ]
        })
        .collect();

    // Is the point named by plane class `r` on segment `si`'s line within its CLOSED extent
    // (endpoints included)? On an endpoint (`r` is one of the two endpoint classes) it is
    // contained — checked by integer identity, because `order_along(x, x)` is not defined to return
    // 0 (the old `strictly_inside` never compared a class with itself). Otherwise it is contained
    // iff it is strictly between the two endpoints (opposite `order_along` signs).
    // ★ **The point is a parameter, and that is what lets the caller hoist it.** `r` is the class
    // naming the crossing, and the point asked about is `{wc, segs[si].wall, r}` — constant for every
    // segment on one wall, so the caller makes the handle once per **wall pair** and every segment
    // there shares its Cramer parts. Taking `r` instead would cap the sharing at one segment's two
    // endpoints, which is what a `_pair` predicate did.
    let closed_contains = |si: usize, at: &ImplicitPoint<'_, WorkingPlane>, r: usize| -> bool {
        let [r0, r1] = end_c[si];
        // ★ **Asked before the predicates, not after.** Integer identity is a *sufficient* condition
        // for containment, so answering it first skips both orientations. It is not the whole test:
        // where four planes meet, one point wears two handles, and `r` may be the group's
        // representative while the segment still remembers the other — which is what the `== 0`
        // arms below catch. Subsumed, not dropped; the order between them is free.
        if r == r0 || r == r1 {
            return true;
        }
        let ds = end_ds[si];
        let (a, b) = (at.orient3d(r0) * ds[0], at.orient3d(r1) * ds[1]);
        a == 0 || b == 0 || a != b
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
        // Split-point plane classes on W's line: every W-segment endpoint, plus every real
        // different-wall crossing (a segment on `o.wall` whose closed extent reaches W's line).
        let mut pts: Vec<usize> = Vec::new();
        {
            watch!(COLLECT);
            #[cfg(test)]
            phase::scale::add(&phase::scale::COLLECT_TRIPS, segs.len());
            for &i in &wall.segs {
                pts.extend(end_c[i]);
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
                // ★ One handle per wall pair — `r`'s segments all ask about this same crossing.
                let at = jd.point(wc, r, w);
                if other.segs.iter().any(|&i| closed_contains(i, &at, w)) {
                    pts.push(r);
                }
            }
        }
        let pts = {
            watch!(SORT);
            // Distinct classes (same class = same point), then ordered along the line.
            pts.sort_unstable();
            pts.dedup();
            pts.sort_by(|&x, &y| match combinatorics::order_along(jd, wc, w, x, y) {
                -1 => std::cmp::Ordering::Less,
                1 => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            });
            // ★ Two DISTINCT classes ordering equal are one point wearing two handles — a four-plane
            // concurrency `{wc, w, ·, ·}`. Record it, and keep one representative as a split point:
            // splitting at both would emit a zero-length piece between them.
            let mut reps: Vec<usize> = Vec::with_capacity(pts.len());
            let mut group: Vec<usize> = Vec::new();
            let flush = |group: &mut Vec<usize>, reps: &mut Vec<usize>, al: &mut Aliases| {
                if let Some(&rep) = group.first() {
                    reps.push(rep);
                }
                if group.len() > 1 {
                    let mut set: Vec<usize> = vec![wc, w];
                    set.extend(group.iter().copied());
                    set.sort_unstable();
                    set.dedup();
                    al.record(jd, &set);
                }
                group.clear();
            };
            for &r in &pts {
                let same = group
                    .first()
                    .is_some_and(|&g| combinatorics::order_along(jd, wc, w, g, r) == 0);
                if !same {
                    flush(&mut group, &mut reps, aliases);
                }
                group.push(r);
            }
            flush(&mut group, &mut reps, aliases);
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
                for &r in &pts {
                    let mut set: Vec<usize> = vec![wc, r];
                    set.extend(family.iter().copied());
                    set.sort_unstable();
                    set.dedup();
                    aliases.record(jd, &set);
                }
            }
        }
        // A split point is named `{wc, w, third}` — then folded, because this rebuilds the name
        // from a handle and so would otherwise re-introduce the very alias `merge_coincident` just
        // removed. Canonicalizing here and in the merge means every name **downstream** is already
        // the canonical one, and no later stage has to know the table exists.
        let sorted = |r: usize, al: &Aliases| al.canon_point(NodeId::three_planes([wc, w, r]));
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
            // The same hoist as the collector: every `w`-segment asks about these two points.
            let (at_p, at_q) = (jd.point(wc, w, p), jd.point(wc, w, q));
            for &i in &wall.segs {
                if closed_contains(i, &at_p, p) && closed_contains(i, &at_q, q) {
                    merged.extend(segs[i].merged.iter().copied());
                }
            }
            if !merged.is_empty() {
                out.push(MergedSeg {
                    wall: w,
                    end: [sorted(p, aliases), sorted(q, aliases)],
                    end_h: [
                        combinatorics::EndPin::Class(p),
                        combinatorics::EndPin::Class(q),
                    ],
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
            // Collinear with the reference: angle 0 (same ray) or π (opposite ray).
            _ if combinatorics::antiparallel(&edges[i], &edges[0]) => pole.push(i),
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
/// edge — segment or arc). The `-1`-winding face count must equal this.
fn component_count(segs: &[MergedSeg], arcs: &[MergedArc]) -> usize {
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
    let ends = segs.iter().map(|s| s.end).chain(arcs.iter().map(|a| a.end));
    for e in ends {
        let (a, b) = (id(e[0], &mut parent), id(e[1], &mut parent));
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..parent.len())
        .filter(|&i| find(&mut parent, i) == i)
        .count()
}

/// **The premise the circle shortcut below rests on, checked.**
///
/// [`walk_cells`] appends each circle as a **one-edge closed cell**, on the pseudo-half-edges past
/// which is only a cell at all while the circle meets no segment: a crossed circle is arcs with
/// endpoints, and those belong in the segment machinery. That premise used to be kept as a side
/// effect of the population gate's wall rule — a rule about something else — so this asks it
/// directly, of the very segments the cells are about to be built from.
///
/// ★ **The segments, not the faces they came from.** A face-level question would be a proxy: what
/// must hold is a fact about the objects the cell/label machinery consumes, and asking anything
/// wider refuses shapes whose segments never come near the circle (measured: ten edges in today's
/// corpus have a line within `r` whose nearest approach lies off the segment).
///
/// ★ Both copies of the per-class pipeline — the boolean's and the audit's replay — reach the
/// cells through [`walk_cells`], so living here is what keeps the two from disagreeing about
/// where a class stopped.
///
/// ★★ **It says where, and the where is the durable half.** [`circle_crossings`] locates what this
/// rejects on: the points where a segment crosses the circle, which are `plane ∩ plane ∩ cylinder`
/// — the shape [`nacre_topo::VertexDef::Branch`] names, and the points the *next* rung will split
/// the circle into arcs at. When that lands, this rejection goes away and the locator stays. The
/// witness is today's reachable consumer of it, not its purpose.
///
/// ★ **A trap left named for that next rung.** `Branch`'s `root` (`Lo`/`Hi`) is defined against
/// the meet line of the two planes taken in **ascending handle order**, and its own doc warns that
/// anything re-sorting the pair must toggle `root`. The arrangement works in *class indices*, so
/// the roots recorded here are ordered by **this function's call order** (`wc`, then the
/// segment's wall). Minting a `Branch` from one means establishing that correspondence, not
/// assuming it.
fn arc_split_witness(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
) -> Option<RejectWhere> {
    // ★ **Collected, not returned at the first hit.** The verdict does not care which pair came
    // first, but the *witness* would: reporting whichever the loops reached first makes the
    // location depend on iteration order. Gathering them all lets the witness be chosen by name
    // (`Break::key`), the rule `NonManifoldResultEdge` established. Only the failing path pays,
    // and that path is a rejection.
    let mut breaks: Vec<Break> = Vec::new();
    for circ in circles {
        let (o, m, r) = (circ.def.origin(), circ.def.dir(), circ.def.radius());
        for sg in segs {
            let ends = sg.end.map(|t| combinatorics::node_coords_rat(jd, t));
            // A node whose coordinates do not fit the road's vessel has no witness to give; the
            // split above already declined for the same reason, so this only shapes the message.
            let [Some(p0), Some(p1)] = ends else {
                return None;
            };
            // ★ A collapsed edge is a *point*, and «is this point inside the disk» is still a
            // real question — skipping it would let the one case it can express slip through.
            // (the split's own `CoincidentNodes` owns the malformed-edge story; this only avoids
            // asking a segment predicate about something that is not a segment.)
            let collapsed = p0 == p1;
            let hit = if collapsed {
                nacre_scalar::cylinder_radial_side(&p0, &o, &m, r) != nacre_scalar::Orient::Positive
            } else {
                nacre_scalar::segment_meets_cylinder(&p0, &p1, &o, &m, r)
            };
            if !hit {
                continue;
            }
            // ★★ **The verdict above is not re-asked below.** "Does a segment meet the solid
            // cylinder" and "where does it cross the circle" are different questions — a segment
            // lying wholly *inside* the disk answers yes to the first and has no crossing at all,
            // and it still breaks the closed-cell premise. So the locator only decides the
            // witness's shape; `hit` alone decides the rejection.
            if collapsed {
                // The line through a collapsed edge is not where the point is; its own coordinate
                // is. Running the locator here would name a place the geometry never visits — the
                // trap is that a collapsed edge still *has* a wall and a line, so the locator
                // would answer, plausibly and wrongly.
                //
                // ★ **Unmeasured, and said so.** A probe across the whole ops suite reaches this
                // branch **zero** times; it is defensive (as it already was before it had a
                // witness to choose). So this rule is *chosen*, not measured — which is why it is
                // written out rather than left to whatever the locator would have returned.
                breaks.push(Break {
                    key: (circ.cyl, sg.wall, sg.end, None),
                    crossing: false,
                    at: RejectWhere::Point(realize(&p0)),
                });
                continue;
            }
            match circle_crossings(jd, wc, circ, sg, [&p0, &p1]) {
                // Overflow in the locator leaves the rejection witnessless rather than turning it
                // into a different rejection — this pass only ever *adds* a place.
                None => {}
                Some(xs) if xs.is_empty() => breaks.push(Break {
                    key: (circ.cyl, sg.wall, sg.end, None),
                    crossing: false,
                    at: RejectWhere::Segment([realize(&p0), realize(&p1)]),
                }),
                // ★ **The place comes from the name**, not from the solve that found it — so a
                // name whose root did not follow its pair through canonical order points here at
                // the other crossing, visibly. A realization that overflows drops that one break
                // (the same rule the `None` arm above follows) rather than moving the rejection.
                Some(xs) => breaks.extend(xs.into_iter().filter_map(|n| {
                    let at = combinatorics::branch_point(jd, circ.cyl, &circ.def, n)?;
                    Some(Break {
                        key: (circ.cyl, sg.wall, sg.end, Some(n)),
                        crossing: true,
                        at: RejectWhere::Point(Point3::from_array(at)),
                    })
                })),
            }
        }
    }
    // ★★★ **Only what *separates* the circle is a break.** The premise above is that a circle
    // stays a *closed cell*, and exactly one kind of contact costs it that — a **transversal
    // crossing**, which cuts the loop into arcs. The other two leave it closed:
    //
    // - a **tangency** meets it at one point (`QuadRoot::Double`, the variant that exists because
    //   "the two roots coincide" had to be sayable);
    // - a **containment** — a segment wholly inside the disk — never meets the circle at all
    //   (no root, so no name), it only puts a polygon in the disk cell. `nest_cells` hosts that
    //   now, which is what let this arm widen.
    //
    // So the test is the crossing's own name: a `Lo`/`Hi` root separates, anything else does not.
    let separates = |b: &Break| {
        matches!(
            b.key.3,
            Some(NodeId::Branch { root, .. }) if root != nacre_topo::QuadRoot::Double
        )
    };
    if !breaks.iter().any(separates) {
        return None;
    }
    witness(breaks)
}

/// **Cut every circle a segment crosses into arcs, and the segments with it.**
///
/// ★★★ This is the arc split. What it does *not* do is find the crossings — [`circle_crossings`]
/// already names them, and has since the guard that used to refuse this population was written.
/// The work here is ordering: around the circle (θ, [`nacre_scalar::quad::circular_order_about_seam`])
/// to make arcs, and along each segment (the line parameter, [`cmp_along`]) to make sub-segments.
///
/// ★ A circle nothing crosses is returned **whole**, on the road it has always taken. The parallel
/// road shrinks to the population it is actually about.
#[allow(clippy::type_complexity)]
fn split_circles(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
) -> Result<Option<(Vec<MergedSeg>, Vec<MergedCircle>, Vec<MergedArc>)>, BoolError> {
    use nacre_scalar::quad::QuadVal;
    let undecided = || reject(RejectReason::WitnessNotRational);
    // Which nodes land on each circle, and which on each segment. Collected together because one
    // crossing is a point of both — splitting only one of them would leave the other's edge running
    // through a vertex it does not have.
    let mut on_circle: Vec<Vec<NodeId>> = vec![Vec::new(); circles.len()];
    let mut on_seg: Vec<Vec<NodeId>> = vec![Vec::new(); segs.len()];
    for (ci, circ) in circles.iter().enumerate() {
        let (o, m, r) = (circ.def.origin(), circ.def.dir(), circ.def.radius());
        for (si, sg) in segs.iter().enumerate() {
            let ends = sg.end.map(|t| combinatorics::node_coords_rat(jd, t));
            let [Some(p0), Some(p1)] = ends else {
                return Err(undecided());
            };
            if p0 == p1 || !nacre_scalar::segment_meets_cylinder(&p0, &p1, &o, &m, r) {
                continue;
            }
            let Some(xs) = circle_crossings(jd, wc, circ, sg, [&p0, &p1]) else {
                return Err(undecided());
            };
            for n in xs {
                // A tangency touches without separating — `separates` above already let that
                // shape through, and cutting there would make a zero-length arc.
                if matches!(
                    n,
                    NodeId::Branch {
                        root: nacre_topo::QuadRoot::Double,
                        ..
                    }
                ) {
                    continue;
                }
                on_circle[ci].push(n);
                on_seg[si].push(n);
            }
        }
    }

    // The cylinder each crossing names, from the very list that made it — a segment's crossings
    // all come from circles on this class, so nothing has to be looked up in the class table.
    let def_of: HashMap<usize, nacre_topo::CylinderDef> =
        circles.iter().map(|c| (c.cyl, c.def.clone())).collect();
    // ★ Nothing crossed: the caller keeps its own slices and no edge is copied. The common case
    // pays for the question and not for an answer it does not need.
    if on_circle.iter().all(|v| v.is_empty()) {
        return Ok(None);
    }
    let mut out_segs = Vec::with_capacity(segs.len());
    let mut out_circles = Vec::new();
    let mut arcs = Vec::new();

    // ★★ **The segment split runs first: the sharper question before the vaguer one.** Both halves
    // refuse a point wearing two names, and the segment side sees the case the circle side cannot —
    // a crossing that lands on an **endpoint**, a three-plane name and a branch name for one point.
    //
    // ★ **Unmeasured, and said so.** Swapping the two leaves every fixture green: the only one that
    // reaches the two-names case cuts its circle at a *single* point, so the circle side's adjacent
    // pairs are empty and it never asks. The order is right on principle and today it decides
    // nothing — what it does decide is that a degenerate `[n, n]` arc is not built before a refusal
    // that was going to happen anyway.
    //
    // ★★ **And it is why a one-node circle needs no name of its own here.** A circle cut at
    // exactly one point is *slit*, not divided, and the arc below comes out `[n, n]` — the closed
    // form a rim has. Measured: the only fixture that reaches it is
    // `a_crossing_on_a_segments_endpoint`, where the single crossing **is** the two-names case and
    // the segment half now refuses it first. A genuine slit (a segment ending strictly inside the
    // disk) would fall through to the walk, whose orbit-length rule refuses a one-edge cycle by
    // name — loudly, and where the sentence is true.
    // ---- segments → sub-segments, in line order ----
    for (si, sg) in segs.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_seg[si]);
        if nodes.is_empty() {
            out_segs.push(sg);
            continue;
        }
        nodes.sort_unstable();
        nodes.dedup();
        // ★★ **Every crossing on this segment shares one line, and that is checked rather than
        // asserted — an unfired guard, measured to be.** (4 of 4 segment splits in the suite carry
        // one line; the check is here because the *sort* below silently mixes two rulers if it ever
        // stops holding, not because a fixture is red.) `branch_meet` re-solves from the *name*, whose plane pair is `{wc, wall}`
        // sorted — the same two planes for every crossing here — and `plane_plane_cylinder` builds
        // the meet line by Cramer from those two planes **before** it looks at the cylinder. So the
        // parameters of crossings made by *different* cylinders are still comparable. That is the
        // whole reason the ends and the crossings can be sorted into one sequence, so a comment
        // is the wrong place for it: if it ever stops holding, the sort silently mixes two rulers.
        let mut line: Option<nacre_scalar::quad::MeetLine> = None;
        let mut keyed: Vec<(QuadVal, NodeId, combinatorics::EndPin)> = Vec::new();
        for &n in &nodes {
            let NodeId::Branch { cyl, .. } = n else {
                return Err(reject(RejectReason::RingNaming));
            };
            let def = def_of.get(&cyl).ok_or_else(undecided)?;
            let (l, s) = combinatorics::branch_meet(jd, cyl, def, n).ok_or_else(undecided)?;
            if let Some(first) = &line {
                if first.base() != l.base() || first.dir() != l.dir() {
                    return Err(reject(RejectReason::RingNaming));
                }
            }
            keyed.push((s, n, combinatorics::EndPin::Cylinder));
            line = Some(l);
        }
        let line = line.ok_or_else(undecided)?;
        for k in 0..2 {
            let p = combinatorics::node_coords_rat(jd, sg.end[k]).ok_or_else(undecided)?;
            keyed.push((
                along(&line, &p).ok_or_else(undecided)?,
                sg.end[k],
                sg.end_h[k],
            ));
        }
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
        // ★ A crossing that lands **on** an endpoint is one point wearing two names — a three-plane
        // one and a branch one — and the DCEL keys vertices by name, so shipping both would make
        // two vertices where there is one. Refusing is honest; folding them is its own step.
        for w in order.windows(2) {
            if cmp_along(&keyed[w[0]].0, &keyed[w[1]].0) == Some(core::cmp::Ordering::Equal) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        // ★★★ **The sub-segments' travel sense, taken from the whole segment's own two named
        // ends — exactly, with no frame arithmetic at all.** Every sub-segment runs the way the
        // segment ran, and the pieces below are emitted in **ascending** `along`; so the only
        // question is whether ascending `along` *is* the segment's direction, and the two ends'
        // own parameters answer it. (The alternative — turning `line.dir` into an `EdgeDir` sense
        // — needs the canonical-vs-stored turn on two classes and an `f64` dot to decide it. This
        // needs neither: `edge_dir` is still the one place a sense is made.)
        let ends = (keyed.len() - 2, keyed.len() - 1); // the two pushed just above, in order
        let forward = match cmp_along(&keyed[ends.0].0, &keyed[ends.1].0) {
            Some(core::cmp::Ordering::Less) => true,
            Some(core::cmp::Ordering::Greater) => false,
            // Equal is the coincident-endpoint case the windows check above already refused, and
            // `None` is a width decline.
            _ => return Err(undecided()),
        };
        let whole = combinatorics::edge_dir(jd, wc, sg.wall, sg.end_h[0], sg.end_h[1])?
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
    // ---- circles → arcs, in θ order about the seam ----
    for (ci, circ) in circles.iter().cloned().enumerate() {
        let mut nodes = std::mem::take(&mut on_circle[ci]);
        if nodes.is_empty() {
            out_circles.push(circ);
            continue;
        }
        nodes.sort_unstable();
        nodes.dedup();
        let meets: Vec<(nacre_scalar::quad::MeetLine, QuadVal)> = nodes
            .iter()
            .map(|&n| combinatorics::branch_meet(jd, circ.cyl, &circ.def, n).ok_or_else(undecided))
            .collect::<Result<_, BoolError>>()?;
        // ★★★ **A crossing on the seam is ordered, not refused — it is the cut point.**
        // `circular_order_about_seam` ranks θ ∈ (0, 2π) and answers `SeamIncident` **by name** for
        // a point at θ = 0, because that point is outside the chart's *total* order. But what arcs
        // need is the **cyclic** order, and a cyclic order tolerates one cut anywhere: the seam
        // point is simply first. (Measured: the very first fixture puts a crossing there — a boss
        // on a plate's edge cuts its own rim exactly on the seam generator, so this is the common
        // case, not an exotic one.)
        //
        // ★ At most one node can be seam-incident: two would be the same point, and a crossing
        // that coincides with another is refused below as the two-names-for-one-point it is.
        let mut seam: Vec<usize> = Vec::new();
        let mut chart: Vec<usize> = Vec::new();
        for (i, meet) in meets.iter().enumerate() {
            let on_seam = match nacre_scalar::quad::circular_order_about_seam(
                &circ.def.origin(),
                &circ.def.dir(),
                &circ.def.ref_dir(),
                (&meet.0, &meet.1),
                (&meet.0, &meet.1),
            ) {
                Some(nacre_scalar::quad::SeamOrder::SeamIncident { first, .. }) => first,
                Some(nacre_scalar::quad::SeamOrder::Ordered(_)) => false,
                None => return Err(undecided()),
            };
            if on_seam { seam.push(i) } else { chart.push(i) }
        }
        if seam.len() > 1 {
            return Err(reject(RejectReason::CoincidentNodes));
        }
        let mut bad_theta = false;
        chart.sort_by(|&i, &j| {
            match nacre_scalar::quad::circular_order_about_seam(
                &circ.def.origin(),
                &circ.def.dir(),
                &circ.def.ref_dir(),
                (&meets[i].0, &meets[i].1),
                (&meets[j].0, &meets[j].1),
            ) {
                Some(nacre_scalar::quad::SeamOrder::Ordered(o)) => o,
                _ => {
                    bad_theta = true;
                    core::cmp::Ordering::Equal
                }
            }
        });
        if bad_theta {
            return Err(undecided());
        }
        let order: Vec<usize> = seam.into_iter().chain(chart).collect();
        // ★ **Two names for one point, on the circle side.** Two walls crossing the circle at one
        // point are two *different* branch names — `dedup` cannot see it, and the θ sort puts them
        // adjacent — so an arc of zero length would follow. The segment side asks this question
        // already; asking it here too is what keeps the two sides from disagreeing about what
        // "one point" means.
        for w in order.windows(2) {
            if matches!(
                nacre_scalar::quad::circular_order_about_seam(
                    &circ.def.origin(),
                    &circ.def.dir(),
                    &circ.def.ref_dir(),
                    (&meets[w[0]].0, &meets[w[0]].1),
                    (&meets[w[1]].0, &meets[w[1]].1),
                ),
                Some(nacre_scalar::quad::SeamOrder::Ordered(
                    core::cmp::Ordering::Equal
                ))
            ) {
                return Err(reject(RejectReason::CoincidentNodes));
            }
        }
        for k in 0..order.len() {
            arcs.push(MergedArc {
                cyl: circ.cyl,
                def: circ.def.clone(),
                end: [nodes[order[k]], nodes[order[(k + 1) % order.len()]]],
                merged: circ.merged.clone(),
            });
        }
    }

    Ok(Some((out_segs, out_circles, arcs)))
}

/// **Where a point sits along a line**, as the one parameter both kinds of point can state: a
/// three-plane end is rational, a branch end is `a + b√c`.
///
/// ★ The line is the **canonical** one — [`combinatorics::branch_meet`] re-solves from the name, so
/// every crossing on one `(class, wall)` pair shares a `base`/`dir` whatever cylinder made it. That
/// is what lets the two kinds be sorted into one sequence at all.
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

/// One place a circle's closed-cell premise is broken, carrying where.
struct Break {
    /// **The identity the witness is chosen by** — the circle's cylinder class, the segment's
    /// wall and canonical endpoint names, and the crossing's own name (`None` where there is no
    /// crossing: a containment has no root, and spelling that `0` was a small lie the type now
    /// refuses). Never the coordinate: choosing the
    /// smallest `f64` would let rounding pick what the user is shown, and never the `Vec`
    /// position, which is the array spelling of the "map's first hit" `NonManifoldResultEdge`
    /// warns about.
    ///
    /// ★ **The wall is in the key so the key is total.** Segments merge by wall *and* endpoint
    /// set, so two of them can carry the same endpoints on different walls; without the wall the
    /// two would tie, and a stable sort would hand the choice straight back to the `Vec` order
    /// this key exists to escape.
    /// ★ The **name**, not its root alone: extracting the root would need a second door out of
    /// the identity, which is exactly what this design has one of. The derived `Ord` does the work
    /// and the order is unchanged — within a `(cyl, wall, ends)` tie the name's `planes` and `cyl`
    /// are fixed, so it discriminates on `root`, and `Lo < Hi` reproduces the old `0 < 1`.
    key: (usize, usize, [NodeId; 2], Option<NodeId>),
    /// A crossing outranks a containment: it names a point of the geometry, where a containment
    /// can only point at the edge that sits inside.
    crossing: bool,
    at: RejectWhere,
}

/// The witness of the whole class — `None` when nothing broke.
fn witness(mut breaks: Vec<Break>) -> Option<RejectWhere> {
    breaks.sort_by_key(|b| (!b.crossing, b.key));
    breaks.into_iter().next().map(|b| b.at)
}

fn realize(p: &[nacre_scalar::Rat; 3]) -> Point3 {
    Point3::from_array([p[0].to_f64(), p[1].to_f64(), p[2].to_f64()])
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
fn circle_crossings(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circ: &MergedCircle,
    sg: &MergedSeg,
    ends: [&[nacre_scalar::Rat; 3]; 2],
) -> Option<Vec<combinatorics::NodeId>> {
    use nacre_scalar::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let w = combinatorics::class_coeffs_rat(jd, wc)?;
    let v = combinatorics::class_coeffs_rat(jd, sg.wall)?;
    let (o, m, r) = (circ.def.origin(), circ.def.dir(), circ.def.radius());
    let (line, roots) = match nacre_scalar::quad::plane_plane_cylinder(&w, &v, &o, &m, r)? {
        CylinderMeet::Pair { line, s } => (line, vec![(QuadRoot::Lo, s[0]), (QuadRoot::Hi, s[1])]),
        // ★ `Double`, not `Lo`: the two roots coincide, so a re-sort must leave the name alone.
        CylinderMeet::Tangent { line, s } => (line, vec![(QuadRoot::Double, QuadVal::from_rat(s))]),
        // ★ Unreachable while the caller's `hit` holds — a segment that meets the solid
        // cylinder has a line that meets its surface, and this line is ⊥ to the axis so it
        // cannot pass inside without crossing. Left returning "no crossings" rather than made
        // loud: the cost of being wrong here is a witness of the wrong shape, and the cost of
        // being wrong the other way is a panic on a valid model.
        CylinderMeet::Miss(_) => return Some(Vec::new()),
        other => unreachable!(
            "a circle's class is ⊥ to the axis, so its meet with any wall is ⊥ to the axis and \
             can only miss, touch or cross the cylinder — got {other:?}"
        ),
    };
    // The two planes pinning the segment's ends on this line. ★ Which side of each is "inside"
    // is read off the **opposite endpoint**, never assumed: the sign conventions of a class's
    // stored coefficients are not this function's to guess.
    let mut fences = Vec::with_capacity(2);
    for k in 0..2 {
        // ★ The fence is a **plane** through the segment's end. An end the arc split pinned with a
        // cylinder has none, and this locator answers `None` — the caller's own policy for a
        // witness it cannot form (it drops the witness, never the verdict).
        let e = combinatorics::class_coeffs_rat(jd, sg.end_h[k].class()?)?;
        let far = ends[1 - k];
        let mut at_far = e[3];
        for i in 0..3 {
            at_far = at_far.checked_add(e[i].checked_mul(far[i])?)?;
        }
        fences.push((e, at_far.numer().signum()));
    }
    let mut out = Vec::new();
    for (root, s) in roots {
        let mut inside = true;
        for (e, want) in &fences {
            // An endpoint-coincident crossing (`Zero`) counts as inside: it is a real point of
            // both the circle and the segment, and the arc split will need it.
            let side = nacre_scalar::quad::plane_side(e, &line, &s);
            let got = match side {
                nacre_scalar::Orient::Positive => 1,
                nacre_scalar::Orient::Negative => -1,
                nacre_scalar::Orient::Zero => 0,
            };
            if got != 0 && got != *want {
                inside = false;
                break;
            }
        }
        if !inside {
            continue;
        }
        out.push(combinatorics::NodeId::branch(wc, sg.wall, circ.cyl, root));
    }
    Some(out)
}

/// Extract the arrangement's cells (faces) from the split 1-skeleton by a DCEL face-walk. Returns
/// the cells plus `face_of[he] = cell index`, which the label brick uses to reach a neighbour cell
/// across an edge as `face_of[twin(he)]`.
///
/// Half-edge encoding: edge `i` gives `he = 2i` (forward, `end[0]→end[1]`) and `2i+1` (reverse);
/// `twin(he) = he ^ 1`. Segments are numbered first and the arcs a split made after them, so `he`
/// alone says which kind it is — see [`walk_cells`]. `next` is set per vertex from `angular_order`: a half-edge arriving at `v`
/// leaves as `twin`, and `next` is `twin`'s **one-step** neighbour in the cyclic order. The step
/// direction (predecessor vs successor) is `angular_order`'s handedness — unknown up front, so both
/// are tried and the one giving exactly `component_count` faces of winding `-1` is kept.
/// **The stopper**: an arc-bounded cell is walked, and then refused.
///
/// ★★★ Everything above it is live — the split, the walk over three half-edge ranges, the winding
/// read at a branch point — and what is still missing is the *assembly*: `Ring` has no carrier for
/// an arc, the rim table keys one closed edge per `(cylinder, plane)`, and `edge_for` keys edges by
/// an unordered vertex pair, which would fold a chord and its two complementary arcs into one.
///
/// ★★ **It reads `has_arcs`, and the condition used to be "was split".** That stricter reading was
/// right while the cells were extracted by a function that split *inside* itself: what it returned
/// was numbered against a
/// slice the caller did not have, so any `Ok` from a split arrangement was unusable whether or not
/// arcs came out of it. The caller holds the `ClassEdges` now and hands the same one to everything
/// below, so that mismatch cannot arise and the condition is the stopper's own sentence again.
/// (`ClassEdges::of` still asserts the two coincide.)
///
/// ★ The witness is read from the **uncut** edges — `arc_split_witness` finds the crossing that
/// made the split, and after the split there is nothing left crossing.
fn arc_stopper(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    edges: &ClassEdges<'_>,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
) -> Result<(), BoolError> {
    if !edges.has_arcs() {
        return Ok(());
    }
    Err(match arc_split_witness(jd, wc, segs, circles) {
        Some(w) => crate::reject_at(RejectReason::ArcBoundNotYet, w),
        None => reject(RejectReason::ArcBoundNotYet),
    })
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
}

/// **The per-class stages, in one place** — the walk and the nesting.
///
/// ★★★★ **It exists because the sequence was written twice.** `trace_result_faces`' `arrange` and
/// `frame_audit`'s replay ran the same stages against the same edges, and the second is what makes
/// a reject debuggable — so a drift between them is, in `decline_to_reject`'s words, *"the worst
/// possible time to be lying"*. One copy, and the audit is auditing what the boolean ran.
///
/// ★★★ **And it puts every stage behind one `Result`, which is the failure this shape is really
/// for.** The caller must run [`arc_stopper`] before unwrapping, so an arc class carries the same
/// name out however far the pipeline got. With the stages held separately that order had to be
/// repeated per stage, and getting it wrong is silent — the suite stays green and only the reject's
/// *name* changes, which is exactly what happened once (`f8eb935`). One `Result` leaves one place
/// to put the `?`, and it is after the stopper.
///
/// ★ **`emit_faces` is the last stage, and it cannot fail.** Every other stage returns a `Result`,
/// so the stopper's interception is what decides the reject's *name*; this one returns its product
/// outright. That is worth reading precisely — it means "the stopper still intercepts" cannot be
/// probed from outside by making this stage fail, and the probe has to be planted here instead.
fn per_class(
    jd: &Judge<'_, WorkingPlane>,
    kind: BoolKind,
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<Staged, BoolError> {
    let (cells, face_of) = timed!(C_EXTRACT, walk_cells(jd, wc, edges))?;
    let nesting = timed!(C_NEST, nest_cells(jd, wc, &cells, edges))?;
    // ★ **The seed is `[false; 4]`, and the argument is why it stays an argument.** The
    // arrangement covers all of space, so its unbounded cells reach infinity, where neither solid
    // is. That is a fact about arranging the *whole* model — restrict the input to a region of
    // space and the unbounded cells become an artifact of the restriction, which is what the
    // parameter records.
    let labels = timed!(
        C_LABEL,
        label_cells(&cells, &face_of, edges, &nesting, [false; 4])
    )?;
    let (faces, disk_labels) = timed!(
        C_EMIT,
        emit_faces(kind, &labels, &cells, edges, jd, wc, &nesting.holes)
    );
    Ok(Staged {
        #[cfg(test)]
        cells,
        #[cfg(test)]
        nesting,
        #[cfg(test)]
        labels,
        faces,
        disk_labels,
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
    circles: std::borrow::Cow<'a, [MergedCircle]>,
    /// Whether [`split_circles`] cut anything. Equal to `!arcs.is_empty()` — see [`Self::of`].
    split: bool,
}

/// Which of the three ranges a half-edge is in, and its index within that range.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum HalfEdgeKind {
    Seg(usize),
    Arc(usize),
    Circle(usize),
}

impl<'a> ClassEdges<'a> {
    /// Run the arc split and hold its result — **the one place a class's edges are assembled.**
    fn of(
        jd: &Judge<'_, WorkingPlane>,
        wc: usize,
        segs: &'a [MergedSeg],
        circles: &'a [MergedCircle],
    ) -> Result<Self, BoolError> {
        use std::borrow::Cow;
        Ok(match split_circles(jd, wc, segs, circles)? {
            Some((s, c, a)) => {
                // ★ The names below rest on this: `split_circles` answers `Some` only when some
                // circle collected a crossing, and a crossing yields at least one arc (a tangency
                // is skipped before it is collected). So "was split" and "has arcs" are one fact.
                debug_assert!(!a.is_empty(), "a split that cut no circle into arcs");
                ClassEdges {
                    segs: Cow::Owned(s),
                    arcs: Cow::Owned(a),
                    circles: Cow::Owned(c),
                    split: true,
                }
            }
            None => ClassEdges {
                segs: Cow::Borrowed(segs),
                arcs: Cow::Borrowed(&[]),
                circles: Cow::Borrowed(circles),
                split: false,
            },
        })
    }

    /// Where the walk's half-edges end and the circles' pseudo-half-edges begin.
    fn he_count(&self) -> usize {
        2 * (self.segs.len() + self.arcs.len())
    }

    fn has_arcs(&self) -> bool {
        debug_assert_eq!(self.split, !self.arcs.is_empty(), "see `of`");
        !self.arcs.is_empty()
    }

    /// **The one classifier.** Everything that used to compare against `2 * segs.len()` asks this.
    fn kind(&self, he: usize) -> HalfEdgeKind {
        let ns = self.segs.len();
        if he < 2 * ns {
            HalfEdgeKind::Seg(he / 2)
        } else if he < self.he_count() {
            HalfEdgeKind::Arc((he - 2 * ns) / 2)
        } else {
            HalfEdgeKind::Circle((he - self.he_count()) / 2)
        }
    }

    /// The vertex a half-edge leaves. `he % 2 == 0` takes `end[0]` — one rule, both ranges.
    fn origin(&self, he: usize) -> NodeId {
        match self.kind(he) {
            HalfEdgeKind::Seg(i) => self.segs[i].end[he % 2],
            HalfEdgeKind::Arc(i) => self.arcs[i].end[he % 2],
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
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge is not a ring edge")
            }
        };
        // The endpoints as handles on this edge's line, carried by the producer. Not recovered from
        // the names: a canonical name need not mention `wc` or the wall (see
        // `combinatorics::RingEdge`). ★ An arc's ends are branch points by construction, which is
        // what `EndPin::Cylinder` says.
        let (from_h, to_h) = match self.kind(he) {
            HalfEdgeKind::Seg(i) => (self.segs[i].end_h[he % 2], self.segs[i].end_h[1 - he % 2]),
            _ => (
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
                edges.push(combinatorics::dir_at(jd, wc, &edge_of(he), v)?);
            }
        }
        let ord = timed!(E_ANGULAR, angular_order(jd, wc, &edges))?;
        cyclic.insert(v, (outs.clone(), ord));
    }

    watch!(E_WALK);
    let components = component_count(segs, arcs);

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
            let w = combinatorics::loop_winding(jd, wc, &ring)?;
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
                (n + arcs.len()) as i64,
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
            // twin arithmetic works unmodified (`2n` is even). The gate proves a circle meets
            // no segment, so no crossing machinery is owed.
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
            HalfEdgeKind::Seg(_) | HalfEdgeKind::Arc(_) => None,
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
            if cell_in_cell(jd, wc, &rings, circles, &circle_ix, c, r)? == Some(true) {
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
            let host = innermost_host(jd, wc, &rings, circles, &circle_ix, &hosts)?;
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
/// witness `nest_cells` uses for a circle contour (the loops are disjoint by the gate's
/// clearance proof, so one point decides). Exact: the center is `axis ∩ W` (rational), the
/// ring corners are rational meets, and the parity runs in a rational 2D basis of `W`
/// (`point_in_ring_2d_rat` — parity is invariant under the affine projection).
fn circle_center_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    circle: &MergedCircle,
    ring: &[combinatorics::RingEdge],
) -> Result<bool, BoolError> {
    use nacre_scalar::Rat;
    // ★ Every refusal below is a **value** that could not be formed exactly — a class with no
    // narrow description, a coordinate past `Rat` — which is the road's name, not the gate's.
    // (The gate's own questions are signs and were made total; borrowing its name here pointed
    // at a layer that had already answered.)
    let undecided = || reject(RejectReason::WitnessNotRational);
    let zero = Rat::from_int(0);
    let coeffs = combinatorics::class_coeffs_rat(jd, wc).ok_or_else(undecided)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let (o, m) = (circle.def.origin(), circle.def.dir());
    let dot3 = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<Rat> {
        x[0].checked_mul(y[0])?
            .checked_add(x[1].checked_mul(y[1])?)?
            .checked_add(x[2].checked_mul(y[2])?)
    };
    let center = (|| -> Option<[Rat; 3]> {
        let nm = dot3(&n, &m)?;
        let no_d = dot3(&n, &o)?.checked_add(coeffs[3])?;
        let t = zero
            .checked_sub(no_d)?
            .checked_mul(Rat::new(nm.denom(), nm.numer())?)?;
        let mut p = o;
        for k in 0..3 {
            p[k] = p[k].checked_add(t.checked_mul(m[k])?)?;
        }
        Some(p)
    })()
    .ok_or_else(undecided)?;
    // The class's rational chart — the one copy of that rule ([`combinatorics::Chart2dRat`]);
    // parity is affine-invariant, so the basis need not be orthonormal.
    let chart = combinatorics::Chart2dRat::of_normal(&n).ok_or_else(undecided)?;
    let p2 = chart.project(&center).ok_or_else(undecided)?;
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let ring2 = chart.ring(jd, &nodes).ok_or_else(undecided)?;
    match nacre_geom::intersect::point_in_ring_2d_rat(p2, &ring2) {
        nacre_geom::intersect::RingSide::Inside => Ok(true),
        nacre_geom::intersect::RingSide::Outside => Ok(false),
        // On the boundary is the gate-impossible contact; refusing is the honest answer.
        nacre_geom::intersect::RingSide::OnBoundary => Err(undecided()),
    }
}

/// Whether a polygon contour lies inside a disk — population-impossible (its edges ride wall
/// faces, all of them clear of the lateral), but computed honestly from one node's
/// radial side rather than assumed.
fn node_in_circle(
    jd: &Judge<'_, WorkingPlane>,
    ring: &[combinatorics::RingEdge],
    circle: &MergedCircle,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let node = ring.first().ok_or_else(undecided)?.node;
    let p = combinatorics::node_coords_rat(jd, node).ok_or_else(undecided)?;
    match nacre_scalar::cylinder_radial_side(
        &p,
        &circle.def.origin(),
        &circle.def.dir(),
        circle.def.radius(),
    ) {
        nacre_scalar::Orient::Negative => Ok(true),
        nacre_scalar::Orient::Positive => Ok(false),
        // On the surface: the gate-impossible contact — a *geometric* degeneracy, so it keeps
        // the gate's name while the width causes above take the road's.
        nacre_scalar::Orient::Zero => Err(reject(RejectReason::CylinderGateUndecided)),
    }
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
/// - **two circles** — never nested here: one circle's disk and its own contour are adjacent, and
///   distinct cylinders are pairwise clear by the gate (`dist > r₁+r₂` forbids containment);
/// - **circle in polygon** — the centre (`axis ∩ wc`, rational) inside the ring, **and the ring not
///   inside the circle**. ★★ That second clause is not belt-and-braces: with disjoint loops the
///   centre test alone says "inside" for *both* nestings when the polygon happens to straddle the
///   centre — a boss standing over a bore, footprint `[7,9]×[9,11]` around the axis `(8,10)`, put
///   the **circle** inside the **square**. `circle_center_in_ring`'s own doc says one point decides
///   *"the loops are disjoint by the gate's clearance proof"*, and disjoint they are; what it does
///   not settle is **which way round**. The other direction does, so the pair is the predicate;
/// - **polygon in circle** — one node's radial side ([`node_in_circle`]), decisive on its own:
///   disjoint loops put every node on one side;
/// - **polygon in polygon** — a ray from each of `a`'s nodes until one is clear
///   ([`combinatorics::ring_in_ring`]).
fn cell_in_cell(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circles: &[MergedCircle],
    circle_ix: &[Option<usize>],
    a: usize,
    b: usize,
) -> Result<Option<bool>, BoolError> {
    match (circle_ix[a], circle_ix[b]) {
        (Some(_), Some(_)) => Ok(None),
        (Some(ci), None) => Ok(Some(
            circle_center_in_ring(jd, wc, &circles[ci], &rings[b])?
                && !node_in_circle(jd, &rings[b], &circles[ci])?,
        )),
        (None, Some(ri)) => node_in_circle(jd, &rings[a], &circles[ri]).map(Some),
        (None, None) => {
            if rings[a]
                .iter()
                .any(|e| rings[b].iter().any(|f| f.node == e.node))
            {
                return Ok(None);
            }
            // ★ `ring_in_ring` casts from each of `a`'s nodes until one gives a clear ray; an
            // exhausted ring is the genuine degeneracy it rejects for — which is also why
            // dropping a branch node from the probe list here is honest.
            let probes = combinatorics::three_plane_probes(rings[a].iter().map(|e| e.node));
            combinatorics::ring_in_ring(jd, wc, &probes, &rings[b]).map(Some)
        }
    }
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
    let inside = |a: usize, b: usize| cell_in_cell(jd, wc, rings, circles, circle_ix, a, b);
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
type Label = [bool; 4];

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
            // An arc is a piece of its circle's trace, so it carries the same contributions.
            HalfEdgeKind::Arc(i) => edge_mask(&edges.arcs[i].merged),
            HalfEdgeKind::Circle(i) => edge_mask(&edges.circles[i].merged),
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
) -> (Vec<LocalFace>, Vec<(usize, Label)>) {
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
        // ★★ **`Ring` has no carrier for an arc, and this is where that runs out.** A ring's walls
        // are plane classes and `assemble_fuse_cut` indexes the class table with them, so an arc
        // half-edge has nothing true to put here; widening `Ring` to a carrier is the next cell's
        // own item. Until then the sentinel goes in and is checked **where it is read**
        // (`boolean::edge_for`), not here: this is the producer, and the true proposition is not
        // "no arc reaches me" — an arc class reaches this function on every run now that the
        // stopper stands behind it — but "nobody indexes the class table with `usize::MAX`".
        crate::boolean::Bound::Ring(crate::boolean::Ring::new(
            cell.half_edges.iter().map(|&he| edges.origin(he)).collect(),
            cell.half_edges
                .iter()
                .map(|&he| match edges.kind(he) {
                    HalfEdgeKind::Seg(i) => edges.segs[i].wall,
                    _ => usize::MAX,
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
            Some((mc.cyl, labels[c]))
        })
        .collect();
    (out, disk_labels)
}

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
) -> Result<(Vec<LocalFace>, DiskLabels), BoolError> {
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
) -> Result<(Vec<LocalFace>, DiskLabels), BoolError> {
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
    let mut splits: Vec<(Vec<MergedSeg>, Vec<MergedCircle>)> = Vec::new();
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
            let mut tr = timed!(TRACE_ON, trace_on_class(trace_in, wc, jd, faces, plane_ix));
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
            let merged = timed!(MERGE, merge_coincident(&tr.segs, wc, &local));
            let split = timed!(SPLIT, split_at_crossings(jd, wc, &merged, &mut local))?;
            let split = drop_newsless(split)?;
            let circles = merge_circles(&tr.circles, cyls)?;
            Ok(((split, circles), local))
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
        let (split, circles) = &splits[k];
        // The per-class product: the faces, and the disk labels the band pass reads (empty for a
        // class with no circles — and for a reused class, which is why cylinders switch reuse off).
        type ClassOut = (Vec<LocalFace>, Vec<(usize, Label)>);
        let arrange = |wc: usize| -> Result<ClassOut, BoolError> {
            watch!(CELLS);
            // ★★ **One `ClassEdges` for the whole pipeline.** The split renumbers half-edges, so
            // everything below must read the *same* edges the walk did — building it here is what
            // makes that structural instead of a promise. (`frame_audit` runs its own copy of this
            // pipeline and must build it the same way; the arc fence locks that they agree.)
            let edges = timed!(C_SPLIT, ClassEdges::of(jd, wc, split, circles))?;
            // ★★★ **The stages run, then the stopper, then the `?`.** The stopper *intercepts*:
            // an arc class must carry the same name out **however far the pipeline got**, or the
            // fences' `ArcBoundNotYet` + witness would become whatever a stage said and the reject
            // census would gain a raise site. Writing `per_class(..)?` puts the stages' failure
            // ahead of the stopper and loses exactly that — measured: with a stage stubbed to fail
            // on an arc class, the fence sees that stage's reason instead.
            //
            // ★★ **The fences that pin the *name* are `bands.rs`'
            // `a_boss_overhanging_the_plates_edge_is_still_refused` and
            // `a_turned_boss_over_the_plates_corner_names_the_crossing_on_the_segment`** — read the
            // probe there. `the_audit_and_the_boolean_agree_about_an_arc_class` compares the audit
            // *with* the boolean, so it stays green when both slide to the same wrong name; asking
            // it about interception measures the proposition next door.
            let staged = per_class(jd, kind, wc, &edges);
            arc_stopper(jd, wc, &edges, split, circles)?;
            let Staged {
                faces, disk_labels, ..
            } = staged?;
            Ok((faces, disk_labels))
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
            Some(f) => Ok((f, Vec::new())),
            None => arrange(wc),
        }
    })?;
    let mut disk_labels: DiskLabels = HashMap::new();
    for (k, (faces, labels)) in per_class.into_iter().enumerate() {
        local_faces.extend(faces);
        for (cyl, label) in labels {
            disk_labels.insert((cyl, work[k]), label);
        }
    }
    Ok((local_faces, disk_labels))
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
    );
    let mut out = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &faces_tab, &plane_ix);
        // Names this class used: segment endpoints, single-point touches, and — since a crossing
        // the arrangement mints is a vertex too — the split's endpoints where it got that far.
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let mut names: Vec<NodeId> = tr.touches.clone();
        for s in merged.iter().chain(
            split_at_crossings(&jd, wc, &merged, &mut Aliases::default())
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
                    combinatorics::face_vertex_triples(model, fh, fp, inc, &jd, &plane_ix)
                        .map(|lr| lr.poly().map(|nr| nr.triples.clone()).unwrap_or_default())
                        .unwrap_or_default();
                rings.extend(
                    combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix)
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
                names.extend(rings.into_iter().filter(|&n| {
                    // The name is hoisted so the call stays on one line: the source-text meta-test
                    // `no_production_code_walks_a_ring_past_the_shared_walk` allow-lists this site
                    // by its argument text, and it reads line by line.
                    let name = three_plane_name(n).expect("a three-plane node");
                    combinatorics::side_of(&jd, name, wc) == 0
                }));
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
    /// ★★★ Without this the arc population is measured and thrown away: the stopper stands after
    /// every arrangement stage and swallows their answer, so a probe is the only way to see it —
    /// and a probe is deleted before the commit. Then the next cell rediscovers a break here while
    /// debugging something else, which is exactly the blame this stage was split to isolate.
    /// `None` when the class stopped before those stages ran.
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
    );
    // The alias fixpoint, exactly as the boolean runs it (sequentially): classes discover names
    // for one another's features, so auditing each class against an empty table is auditing a
    // *different* pipeline — it diverged from the boolean's answer the day the rounds arrived,
    // and the divergence surfaced when this file's test moved to a fixture that needs them.
    let mut aliases = Aliases::default();
    loop {
        let before = aliases.len();
        for wc in 0..geom.len() {
            let mut tr = trace_on_class(&trace_in, wc, &jd, &faces_tab, &plane_ix);
            aliases.absorb(&std::mem::take(&mut tr.aliases));
            if tr.declined.is_empty() {
                let merged = merge_coincident(&tr.segs, wc, &aliases);
                let _ = split_at_crossings(&jd, wc, &merged, &mut aliases);
            }
        }
        if aliases.len() == before {
            break;
        }
    }
    let mut out = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for wc in 0..geom.len() {
        let tr = trace_on_class(&trace_in, wc, &jd, &faces_tab, &plane_ix);
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
                let merged = merge_coincident(&tr.segs, wc, &local);
                let split = split_at_crossings(&jd, wc, &merged, &mut local)?;
                let split = drop_newsless(split)?;
                let circles = merge_circles(&tr.circles, &cyls)?;
                // ★★★★ **The same edges and the same stopper the boolean uses.** This copy of the
                // pipeline is what makes the audit an instrument; feeding it un-split edges would
                // let it run to the end and report `failed_at: None` for an input the boolean
                // refuses — *"the worst possible time to be lying"* (`decline_to_reject`). The
                // arc fence in `bands.rs` locks the two together.
                let edges = ClassEdges::of(&jd, wc, &split, &circles)?;
                // ★ The same stages, the same order, the same stopper the boolean runs — the arc
                // fence in `bands.rs` locks that the two agree.
                let staged = per_class(&jd, kind, wc, &edges);
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
                arc_stopper(&jd, wc, &edges, &split, &circles)?;
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
            )
        );
        // `disk_labels`: per (cylinder class, plane class), the four bits of that circle's disk
        // cell — what the band pass reads instead of casting a witness ray (M6-2a K1).
        let (faces, disk_labels) = trace_result_faces(
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
                Ok((p, _)) => assert_eq!(
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
        // Clean the raw arrangement output: merge coplanar, same-normal faces that share a full edge
        // (e.g. the split side walls a fused coincident interface leaves) so the result is a minimal,
        // chainable solid — a second boolean on it then sees no redundant coplanar planes.
        let faces = timed!(UNIFY, crate::boolean::unify_coplanar_faces(faces, &jd))?;

        // Build the SeamVertex weld table directly from the emitted triples (no `build_seam`: that is
        // raw-index and pierce-only). Reject rather than panic on a degenerate meet.
        // ★ **The lateral bands, appended after the differential above** (M6-2a C4b): reuse can
        // only change what the *plane* arrangement emits, so the two routes are compared on that
        // list; the bands are a separate pass over the same operands and belong to neither route.
        // From here on there is one face list — the grouping, the closed-shell guard and the
        // assembly all read it.
        let faces = if cyls.is_empty() {
            faces
        } else {
            let rows = crate::bands::cyl_rows(&faces_tab, &plane_ix, n_a)?;
            let mut faces = faces;
            faces.extend(crate::bands::band_faces(
                kind,
                &faces,
                &rows,
                &jd,
                &disk_labels,
            )?);
            faces
        };

        let seam = {
            watch!(SEAM);
            let mut seam: Vec<SeamVertex> = Vec::new();
            let mut seen: HashMap<NodeId, ()> = HashMap::new();
            for f in &faces {
                for loop_ in f.poly_rings() {
                    for &node in loop_.iter() {
                        if seen.insert(node, ()).is_some() {
                            continue;
                        }
                        // ★ **The one site behind the door that materializes a coordinate and
                        // mints a `VertexDef`** — so this is where the next rung starts minting
                        // `VertexDef::Branch`, and where the class-order/handle-order
                        // correspondence has to be established a *second* time
                        // (`NodeId::Branch` is canonical in plane **classes**, `VertexDef::Branch`
                        // in `Handle<Surface>` index, and the class→surf map is not monotone).
                        //
                        // Declining, never `continue`: a skipped seam entry surfaces downstream as
                        // `MissingSeam`, whose class is `SuspectedDefect` and whose sentence is
                        // "a reconstruction dropped a crossing" — a wrong diagnosis for an input
                        // the kernel simply does not build yet.
                        let t = three_plane_name(node)
                            .ok_or_else(|| reject(RejectReason::BranchVertexUnnamed))?;
                        let point =
                            three_planes(&geom[t[0]].plane, &geom[t[1]].plane, &geom[t[2]].plane)
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
            seam
        };

        timed!(
            ASSEMBLE,
            assemble_fuse_cut(model, a, b, &jd, &seam, &faces, &cyls)
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &plane_ix,
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &plane_ix,
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
                &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
            );
            let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
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
                &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
            );
            let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        assert_eq!(merged.len(), 8, "8 merged edges before split");

        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();

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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();

        let edges = ClassEdges::of(&jd, wc, &split, &[]).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();

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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
        let edges = ClassEdges::of(&jd, wc, &split, &[]).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &edges).unwrap();
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
        let edges = ClassEdges::of(&jd, wc, &split, &[]).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &edges).unwrap();
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
        let edges = ClassEdges::of(&jd, wc, &split, &[]).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &edges).unwrap();
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
        let (fuse, _) = emit_faces(
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
        let (cut, _) = emit_faces(
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
        let (common, _) = emit_faces(
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
        );
        let (faces, _) = trace_result_faces(
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
        );
        let (faces, _) = trace_result_faces(
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
        );
        let (faces, _) = trace_result_faces(
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
        );
        let (faces, _) = trace_result_faces(
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
                &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
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
        );
        let (faces, _) = trace_result_faces(
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
            &faces_tab,
            &surf_ix,
            &inc_b,
            &plane_ix,
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

        let circles = merge_circles(&tr.circles, &cyls).unwrap();
        let edges = ClassEdges::of(&jd, wc, &[], &circles).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();
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
        let nesting = nest_cells(&jd, wc, &cells, &edges).unwrap();
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
        let (out, _) = emit_faces(
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        assert!(tr.declined.is_empty(), "{:?}", tr.declined);
        assert!(
            matches!(
                tr.circles[..],
                [CircleTrace {
                    cyl: 0,
                    solid: SolidSide::B,
                    kind: SegKind::Transversal { .. },
                }]
            ),
            "{:?}",
            tr.circles
        );

        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let circles = merge_circles(&tr.circles, &cyls).unwrap();
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
        assert_eq!(
            split.len(),
            4,
            "the cap ring alone — the circle owes the splitter nothing"
        );

        let edges = ClassEdges::of(&jd, wc, &split, &circles).unwrap();
        let (cells, face_of) = walk_cells(&jd, wc, &edges).unwrap();
        assert_eq!(cells.len(), 4, "cap ±1 and circle ±1: {cells:?}");
        let at = |he: usize| cells.iter().position(|c| c.half_edges == [he]).unwrap();
        let (disk, contour) = (at(2 * split.len()), at(2 * split.len() + 1));
        let cap = cells
            .iter()
            .position(|c| c.winding == 1 && c.half_edges.len() == 4)
            .unwrap();

        let nesting = nest_cells(&jd, wc, &cells, &edges).unwrap();
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

        let (out, _) = emit_faces(
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
        );
        let (fp, fl) = input.faces[0]
            .iter()
            .find(|(fp, _)| matches!(plane_ix[*fp], ClassIx::Plane(c) if c == wc))
            .expect("the box's bottom cap sits on wc");
        let doctored = combinatorics::FaceLoops {
            outer: fl.outer.clone(),
            holes: Some(vec![combinatorics::LoopRing::Circle { cyl: 0 }]),
        };
        let mut tr = Trace::default();
        trace_one(
            &[(*fp, doctored)],
            SolidSide::A,
            wc,
            &jd,
            &faces_tab,
            &plane_ix,
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
                [CircleTrace { cyl: 0, solid: SolidSide::A, kind: SegKind::Seated { body_above: ba } }]
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
        );
        assert!(cap.declined.is_empty(), "{:?}", cap.declined);
        let sides: Vec<(bool, bool)> = cap
            .circles
            .iter()
            .filter_map(|c| match c.kind {
                SegKind::Seated { body_above } => Some((false, body_above)),
                SegKind::Graze { body_above } => Some((true, body_above)),
                SegKind::Transversal { .. } => None,
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
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
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
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
        // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
        // seated≡transversal coincidence).
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
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
                &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
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
                );
                for wc in 0..setup.geom.len() {
                    let tr = trace_on_class(&trace_in, wc, &jd, &setup.planes, &setup.plane_ix);
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
