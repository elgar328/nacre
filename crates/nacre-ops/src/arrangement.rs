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
use crate::planes::*;
use crate::tolerant::Judge;
#[cfg(test)]
use crate::transform::transform;
use nacre_cip::Decision;
#[cfg(test)]
use nacre_cip::Pt3;
use nacre_cip::predicate::Notes;
#[cfg(test)]
use nacre_geom::Plane;
use nacre_geom::intersect::three_planes;
#[cfg(test)]
use std::collections::HashSet;

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
    pub end: [[usize; 3]; 2],
    /// The same two endpoints as **handles on this segment's own line**: the third plane pinning
    /// each on `wc ∩ wall`. Carried, not recovered from `end`, because a canonical name need not
    /// mention either of this line's planes.
    pub end_h: [usize; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// Names this trace found to denote one feature — see [`Aliases`].
    pub aliases: Aliases,
    /// Single-point tangential contacts — real arrangement vertices, but not segments (a
    /// zero-length chord would abort at `Line::through_points`).
    pub touches: Vec<[usize; 3]>,
    /// `(face index, reason)` for every face this brick could not trace.
    pub declined: Vec<(usize, DeclineKind)>,
}

/// Sort a triple whose elements are already dense plane ids.
fn sorted3(mut t: [usize; 3]) -> [usize; 3] {
    t.sort_unstable();
    t
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
fn decline_to_reject(kind: DeclineKind, face: Handle<Face>) -> RejectReason {
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
    point: HashMap<[usize; 3], [usize; 3]>,
    /// Union-find over `(class, wall)` — walls of one class that carry the same line.
    wall: HashMap<(usize, usize), usize>,
}

impl Aliases {
    /// Record that every plane in `s` (sorted, deduped) passes through one point.
    fn record(&mut self, jd: &Judge<'_, PlaneGeom>, s: &[usize]) {
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
                        names.push(t);
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

    fn find_point(&self, t: [usize; 3]) -> [usize; 3] {
        let mut x = t;
        while let Some(&p) = self.point.get(&x) {
            if p == x {
                break;
            }
            x = p;
        }
        x
    }

    fn union_point(&mut self, a: [usize; 3], b: [usize; 3]) {
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
    pub(crate) fn canon_point(&self, t: [usize; 3]) -> [usize; 3] {
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

fn plane_ring(ts: &[[usize; 3]]) -> Option<Vec<[usize; 3]>> {
    ts.iter()
        .map(|&t| {
            let mut c = t;
            c.sort_unstable(); // sorted, so equal neighbours catch every duplicate
            (c[0] != c[1] && c[1] != c[2]).then_some(c)
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
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    plane_ix: &[usize],
    out: &mut Trace,
) {
    // `fp` names a *face* (`n_out`, `orient`, the declined log); `fc` names the *plane class* it
    // lies on (triples, comparisons, predicate arguments). Every ring below is in class form, so
    // the two must not be confused — see `canon_ring`.
    let fc = plane_ix[fp];
    let outer = match &loops.outer {
        Some(r) => match plane_ring(r) {
            Some(r) => r,
            None => {
                out.declined.push((fp, DeclineKind::CollapsedTriple));
                return;
            }
        },
        None => {
            out.declined.push((fp, DeclineKind::OuterRing));
            return;
        }
    };
    let mut holes: Vec<Vec<[usize; 3]>> = Vec::new();
    // A hole whose ring cannot be named is not "no hole" — swallowing the error would trace the
    // face as solid where it is pierced, which is a silent wrong answer rather than a reject.
    let Some(raw_holes) = &loops.holes else {
        out.declined.push((fp, DeclineKind::HoleRing));
        return;
    };
    for r in raw_holes {
        match plane_ring(r) {
            Some(r) => holes.push(r),
            None => {
                out.declined.push((fp, DeclineKind::CollapsedTriple));
                return;
            }
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
    'rings: for ring in std::iter::once(&outer).chain(holes.iter()) {
        let n = ring.len();
        let side: Vec<i8> = (0..n)
            .map(|i| combinatorics::side_of(jd, ring[i], wc))
            .collect();
        let Some(start) = side.iter().position(|&s| s != 0) else {
            // Every vertex on `W`: a ring lying in the cut plane is degenerate here.
            declined = Some(DeclineKind::AllOnPlane);
            break;
        };
        let mut j = 0;
        while j < n {
            let i = (start + j) % n;
            if side[i] != 0 {
                let ni = (i + 1) % n;
                if side[ni] != 0 && side[ni] != side[i] {
                    // Strict crossing on edge i; its wall is the plane the edge rides besides fp.
                    match combinatorics::ring_from_names(fc, ring).map(|es| es[i].wall) {
                        Ok(wall) => nodes.push(Node {
                            r: wall,
                            flip: true,
                            run: None,
                            flanks_differ: false,
                            single_touch: false,
                            graze_above: None,
                        }),
                        Err(_) => declined = Some(DeclineKind::CrossingName),
                    }
                }
                j += 1;
            } else {
                // A maximal run of side==0 vertices. Two is the common case, but a vertex whose
                // name had to be taken from its touching planes (`loop_triples`) stays in the ring
                // even when the loop runs straight through it, so a run can be longer: an on-line
                // *interval* with any number of named points inside it.
                let run_start = i;
                let mut m = 0;
                while j < n && side[(start + j) % n] == 0 {
                    m += 1;
                    j += 1;
                }
                let before = side[(run_start + n - 1) % n];
                let after = side[(start + j) % n];
                let flanks_differ = before != after;
                let name = |k: usize, out: &mut Trace| third_on_l(ring[(run_start + k) % n], out);
                if m == 1 {
                    match name(0, out) {
                        Ok(r) => nodes.push(Node {
                            r,
                            flip: flanks_differ,
                            run: None,
                            flanks_differ,
                            single_touch: !flanks_differ,
                            // A single-vertex tangential touch is a point (`touches`), never a graze
                            // segment; a flank-differing single vertex is a strict crossing.
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
                            // Every run is one-sided: it is an *edge* of `f` lying on `L`, so `f`
                            // is on one side of it whatever the ring does afterwards.
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
            wall: plane_ix[fp],
            end: [sorted3([wc, fc, a]), sorted3([wc, fc, b])],
            end_h: [a, b],
            solid: which,
            kind: match graze {
                Some(body_above) => SegKind::Graze { body_above },
                None => SegKind::Transversal {
                    mat: faces[fp].orient_sign,
                },
            },
        });
    };
    for k in 0..nodes.len() {
        if nodes[k].single_touch && parity == 0 {
            out.touches.push(sorted3([wc, fc, nodes[k].r]));
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
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    wc: usize,
    fc: usize,
    fp: usize,
    rs: &[usize],
) -> bool {
    let planes = jd.planes;
    let t = combinatorics::order_along(jd, wc, fc, rs[0], rs[rs.len() - 1]);
    let sigma = faces[fp].n_out.dot(planes[fc].plane.normal());
    (t < 0) == (sigma > 0.0)
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
///
/// `side_faces` is that solid's face slots and `loops` every face's rings — see
/// [`combinatorics::TraceInput`]. This walk used to read the shells out of the `Model` and derive
/// the rings here, once per class; nothing about either depends on `wc`.
#[allow(clippy::too_many_arguments)]
fn trace_one(
    side_faces: &[usize],
    loops: &[combinatorics::FaceLoops],
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    plane_ix: &[usize],
    out: &mut Trace,
) {
    let planes = jd.planes;
    let w_normal = planes[wc].plane.normal();
    for &fp in side_faces {
        let fl = &loops[fp];
        if plane_ix[fp] != wc {
            trace_transversal_face(fp, fl, which, wc, jd, faces, plane_ix, out);
            continue;
        }
        // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
        // side of W — `n_out` points away from the body, so the body is above W exactly when
        // `n_out · n_W < 0`. The f64 sign is robust even rotated: seated means `canon[fp]==wc`,
        // so `n_out ∥ n_W` (both unit) and the dot is ≈ ±1, a full unit from the sign boundary
        // (the rotated-tunnel tests exercise this seated path through the cube's own caps).
        let body_above = faces[fp].n_out.dot(w_normal) < 0.0;
        let kind = SegKind::Seated { body_above };
        // Collect every ring in class form first: a collapsed name declines the whole face, and
        // deciding that before the emitting closure exists keeps the two borrows apart.
        let Some(outer) = fl.outer.as_deref().and_then(plane_ring) else {
            out.declined.push((fp, DeclineKind::OuterRing));
            continue;
        };
        let mut rings = vec![outer];
        let mut collapsed = false;
        // As above: an unnameable hole is a reject, not "no hole".
        let Some(raw) = &fl.holes else {
            out.declined.push((fp, DeclineKind::HoleRing));
            continue;
        };
        for r in raw {
            match plane_ring(r) {
                Some(r) => rings.push(r),
                None => collapsed = true,
            }
        }
        if collapsed {
            out.declined.push((fp, DeclineKind::CollapsedTriple));
            continue;
        }
        let fc = plane_ix[fp];
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
                    end: [sorted3(t0), sorted3(t1)],
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

/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
fn trace_on_class(
    input: &combinatorics::TraceInput,
    wc: usize,
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    plane_ix: &[usize],
) -> Trace {
    let mut out = Trace::default();
    for (side, which) in [SolidSide::A, SolidSide::B].into_iter().enumerate() {
        trace_one(
            &input.side_faces[side],
            &input.loops,
            which,
            wc,
            jd,
            faces,
            plane_ix,
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
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &combinatorics::EdgeFaces,
    inc_b: &combinatorics::EdgeFaces,
    plane_ix: &[usize],
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
/// identical, and only `side_faces[0]` is read.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn trace_one_of(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &combinatorics::EdgeFaces,
    plane_ix: &[usize],
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
    trace_one(
        &input.side_faces[0],
        &input.loops,
        which,
        wc,
        jd,
        faces,
        plane_ix,
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
    #[cfg_attr(test, allow(dead_code))]
    pub end: [[usize; 3]; 2],
    /// The endpoints as handles on this edge's line — see [`Seg::end_h`].
    pub end_h: [usize; 2],
    /// Every `(solid, kind)` that produced this one geometric edge. Length 1 when nothing was
    /// coincident.
    pub merged: Vec<(SolidSide, SegKind)>,
}

/// Merge segments that are the **same geometric edge** — same `wall` and same endpoint-triple set
/// (direction-independent) — into one `MergedSeg`, collecting their contributions. Partial overlap
/// (same `wall`, *different* extent — the E5 case) is left for the per-wall interval overlay in
/// `split_at_crossings` to resolve into non-overlapping sub-segments with unioned contributions.
fn merge_coincident(segs: &[Seg], wc: usize, aliases: &Aliases) -> Vec<MergedSeg> {
    // Key an edge by (wall, sorted endpoint pair) — both folded onto their canonical names first,
    // because two producers can describe one edge with a different wall *and* different endpoint
    // names when planes are concurrent, and it takes both folds for the two keys to coincide.
    let key = |s: &Seg| -> (usize, [[usize; 3]; 2]) {
        let (mut a, mut b) = (aliases.canon_point(s.end[0]), aliases.canon_point(s.end[1]));
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        (aliases.canon_wall(wc, s.wall), [a, b])
    };
    let mut order: Vec<(usize, [[usize; 3]; 2])> = Vec::new();
    let mut groups: HashMap<(usize, [[usize; 3]; 2]), MergedSeg> = HashMap::new();
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
                    end_h: s.end_h,
                    merged: Vec::new(),
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
fn split_at_crossings(
    jd: &Judge<'_, PlaneGeom>,
    wc: usize,
    segs: &[MergedSeg],
    aliases: &mut Aliases,
) -> Result<Vec<MergedSeg>, BoolError> {
    // Is the point named by plane class `r` on `s`'s line within `s`'s CLOSED extent (endpoints
    // included)? On an endpoint (`r` is one of the two endpoint classes) it is contained — checked
    // by integer identity, because `order_along(x, x)` is not defined to return 0 (the old
    // `strictly_inside` never compared a class with itself). Otherwise it is contained iff it is
    // strictly between the two endpoints (opposite `order_along` signs).
    let closed_contains = |s: &MergedSeg, r: usize| -> Option<bool> {
        let (r0, r1) = (s.end_h[0], s.end_h[1]);
        let (a, b) = (
            combinatorics::order_along(jd, wc, s.wall, r, r0),
            combinatorics::order_along(jd, wc, s.wall, r, r1),
        );
        // An endpoint is contained. Integer identity is not the whole test: where four planes meet,
        // one point wears two handles, and `r` may be the group's representative while the segment
        // still remembers the other. Ordering equal is the test — and `order_along(x, x)` is 0 by
        // the same predicate, so the integer check is subsumed, not dropped.
        if r == r0 || r == r1 || a == 0 || b == 0 {
            return Some(true);
        }
        Some(a != b)
    };

    // Group segment indices by wall (walls in first-appearance order for deterministic output).
    let mut walls: Vec<usize> = Vec::new();
    let mut by_wall: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        by_wall
            .entry(s.wall)
            .or_insert_with(|| {
                walls.push(s.wall);
                Vec::new()
            })
            .push(i);
    }

    let mut out = Vec::new();
    for &w in &walls {
        // Split-point plane classes on W's line: every W-segment endpoint, plus every real
        // different-wall crossing (a segment on `o.wall` whose closed extent reaches W's line).
        let mut pts: Vec<usize> = Vec::new();
        for &i in &by_wall[&w] {
            pts.extend(segs[i].end_h);
        }
        for o in segs {
            if o.wall == w {
                continue;
            }
            if jd.plane_pair_dir_sign(wc, w, o.wall) == 0 {
                continue; // walls meet wc in no point (parallel)
            }
            if closed_contains(o, w) == Some(true) {
                pts.push(o.wall);
            }
        }
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
        let pts = reps;
        // A split point is named `{wc, w, third}` — then folded, because this rebuilds the name
        // from a handle and so would otherwise re-introduce the very alias `merge_coincident` just
        // removed. Canonicalizing here and in the merge means every name **downstream** is already
        // the canonical one, and no later stage has to know the table exists.
        let sorted = |r: usize, al: &Aliases| {
            let mut t = [wc, w, r];
            t.sort_unstable();
            al.canon_point(t)
        };
        // Each sub-interval [p, q] carries the union of the W-segments that cover it. Every segment
        // endpoint is itself a split point, so "covers both ends" means "spans the whole interval".
        for pair in pts.windows(2) {
            let (p, q) = (pair[0], pair[1]);
            let mut merged: Vec<(SolidSide, SegKind)> = Vec::new();
            for &i in &by_wall[&w] {
                let s = &segs[i];
                if closed_contains(s, p) == Some(true) && closed_contains(s, q) == Some(true) {
                    merged.extend(s.merged.iter().copied());
                }
            }
            if !merged.is_empty() {
                out.push(MergedSeg {
                    wall: w,
                    end: [sorted(p, aliases), sorted(q, aliases)],
                    end_h: [p, q],
                    merged,
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
/// **Coverage boundary — it does not handle the general case:** two *different* fp's whose lines
/// are parallel (a `0`-turn that is not same-fp) fall into the angle-0 bucket unresolved, and
/// collinear same-direction overlap (E5) is out of scope here. The corpus's arrangement vertices
/// are degree ≥ 3 with distinct fp's per real direction, which is what this handles.
fn angular_order(jd: &Judge<'_, PlaneGeom>, w: usize, edges: &[(usize, i8)]) -> Vec<usize> {
    let planes = jd.planes;
    let os = planes[w].frame_sign;
    let cross = |i: usize, j: usize| -> i8 {
        edges[i].1 * edges[j].1 * jd.plane_pair_dir_sign(w, edges[i].0, edges[j].0) * os
    };
    let (mut zero, mut pos, mut pole, mut neg) = (vec![], vec![], vec![], vec![]);
    for i in 0..edges.len() {
        match cross(0, i) {
            c if c > 0 => pos.push(i),
            c if c < 0 => neg.push(i),
            // Collinear with the reference: angle 0 (same ray) or π (opposite ray).
            _ if edges[i].0 == edges[0].0 && edges[i].1 != edges[0].1 => pole.push(i),
            _ => zero.push(i),
        }
    }
    // Within an open half-plane, `a` precedes `b` (smaller angle) iff `d_a × d_b > 0`.
    let by_turn = |v: &mut Vec<usize>| {
        v.sort_by(|&a, &b| {
            if cross(a, b) > 0 {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        })
    };
    by_turn(&mut pos);
    by_turn(&mut neg);
    let mut out = zero;
    out.extend(pos);
    out.extend(pole);
    out.extend(neg);
    out
}

/// One face of the arrangement: the cyclic list of half-edges bounding it, and its winding
/// (`+1` a bounded island, `-1` the unbounded outer contour).
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    pub half_edges: Vec<usize>,
    pub winding: i8,
}

/// Number of connected components of the 1-skeleton (union-find over vertex triples joined by each
/// segment). The `-1`-winding face count must equal this.
fn component_count(segs: &[MergedSeg]) -> usize {
    let mut idx: HashMap<[usize; 3], usize> = HashMap::new();
    let mut id = |t: [usize; 3], parent: &mut Vec<usize>| -> usize {
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
    for s in segs {
        let (a, b) = (id(s.end[0], &mut parent), id(s.end[1], &mut parent));
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..parent.len())
        .filter(|&i| find(&mut parent, i) == i)
        .count()
}

/// Extract the arrangement's cells (faces) from the split 1-skeleton by a DCEL face-walk. Returns
/// the cells plus `face_of[he] = cell index`, which the label brick uses to reach a neighbour cell
/// across an edge as `face_of[twin(he)]`.
///
/// Half-edge encoding: segment `i` gives `he = 2i` (forward, `end[0]→end[1]`) and `2i+1` (reverse);
/// `twin(he) = he ^ 1`. `next` is set per vertex from `angular_order`: a half-edge arriving at `v`
/// leaves as `twin`, and `next` is `twin`'s **one-step** neighbour in the cyclic order. The step
/// direction (predecessor vs successor) is `angular_order`'s handedness — unknown up front, so both
/// are tried and the one giving exactly `component_count` faces of winding `-1` is kept.
fn extract_cells(
    jd: &Judge<'_, PlaneGeom>,
    wc: usize,
    segs: &[MergedSeg],
) -> Result<(Vec<Cell>, HashMap<usize, usize>), BoolError> {
    let n = segs.len();
    let he_count = 2 * n;
    let origin = |he: usize| segs[he / 2].end[he % 2]; // he%2==0: end[0]; ==1: end[1]
    // The endpoints as handles on this edge's line, carried by the segment. Not recovered from the
    // names: a canonical name need not mention `wc` or the wall (see `combinatorics::RingEdge`).
    let origin_h = |he: usize| segs[he / 2].end_h[he % 2];
    let target_h = |he: usize| segs[he / 2].end_h[1 - he % 2];
    let wall = |he: usize| segs[he / 2].wall;

    // Outgoing half-edges per vertex.
    let mut outgoing: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for he in 0..he_count {
        outgoing.entry(origin(he)).or_default().push(he);
    }

    // For each vertex, the cyclic order of its outgoing half-edges (indices into its `outs` list).
    let mut cyclic: HashMap<[usize; 3], (Vec<usize>, Vec<usize>)> = HashMap::new();
    for (&v, outs) in &outgoing {
        let mut edges = Vec::with_capacity(outs.len());
        for &he in outs {
            let (rv, rf) = (origin_h(he), target_h(he));
            // Direction sign away from v toward the far end (edge_sign convention).
            let s = combinatorics::order_along(jd, wc, wall(he), rf, rv);
            if s == 0 {
                return Err(reject(RejectReason::CoincidentNodes));
            }
            edges.push((wall(he), s));
        }
        let ord = angular_order(jd, wc, &edges);
        cyclic.insert(v, (outs.clone(), ord));
    }

    let components = component_count(segs);

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
            if !ok || cyc.len() < 3 {
                ok = false;
                break;
            }
            // The cell's edges come from the walk, which knows each one's wall and both handles.
            // The endpoint names cannot supply them: a canonical name need not mention this line's
            // planes at all, which is what retired the derivation check that used to stand here.
            let ring: Vec<combinatorics::RingEdge> = cyc
                .iter()
                .map(|&h| combinatorics::RingEdge {
                    node: origin(h),
                    wall: wall(h),
                    from_h: origin_h(h),
                    to_h: target_h(h),
                })
                .collect();
            let w = combinatorics::loop_winding(jd, wc, &ring)?;
            cells.push(Cell {
                half_edges: cyc,
                winding: w,
            });
        }
        if ok && cells.iter().filter(|c| c.winding == -1).count() == components {
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
    jd: &Judge<'_, PlaneGeom>,
    wc: usize,
    cells: &[Cell],
    segs: &[MergedSeg],
) -> Result<Nesting, BoolError> {
    let n = cells.len();
    // As in `extract_cells`: the walk knows each edge's wall and handles, so the ring carries them
    // instead of leaving them to be re-derived from the endpoint names.
    let ring_of = |c: &Cell| -> Vec<combinatorics::RingEdge> {
        c.half_edges
            .iter()
            .map(|&he| combinatorics::RingEdge {
                node: segs[he / 2].end[he % 2],
                wall: segs[he / 2].wall,
                from_h: segs[he / 2].end_h[he % 2],
                to_h: segs[he / 2].end_h[1 - he % 2],
            })
            .collect()
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

    let mut holes: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut roots: Vec<usize> = Vec::new();
    for c in (0..n).filter(|&i| cells[i].winding == -1) {
        let mut hosts: Vec<usize> = Vec::new();
        for &r in &pos {
            // Shares a node ⇒ adjacent (or `c`'s own partner) ⇒ not a hole of `r`.
            if rings[c]
                .iter()
                .any(|e| rings[r].iter().any(|f| f.node == e.node))
            {
                continue;
            }
            // Vertex-disjoint: one clear ray settles it. Retry past a spoiled (ring-node) ray;
            // all of `c`'s vertices spoiled against `r` is a genuine degeneracy → honest reject.
            let mut inside = None;
            for v in rings[c].iter().map(|e| e.node) {
                if let Ok(hit) = combinatorics::point_in_ring(jd, wc, v, &rings[r]) {
                    inside = Some(hit);
                    break;
                }
            }
            match inside {
                Some(true) => hosts.push(r),
                Some(false) => {}
                None => return Err(reject(RejectReason::NoClearRay)),
            }
        }
        if hosts.is_empty() {
            roots.push(c);
        } else {
            let host = innermost_host(jd, wc, &rings, &hosts)?;
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
    jd: &Judge<'_, PlaneGeom>,
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    hosts: &[usize],
) -> Result<usize, BoolError> {
    let inside = |a: usize, b: usize| -> Option<bool> {
        if rings[a]
            .iter()
            .any(|e| rings[b].iter().any(|f| f.node == e.node))
        {
            return None; // adjacent, not nested — not comparable
        }
        rings[a]
            .iter()
            .find_map(|e| combinatorics::point_in_ring(jd, wc, e.node, &rings[b]).ok())
    };
    let mut found = None;
    for &h in hosts {
        if hosts.iter().all(|&o| o == h || inside(h, o) == Some(true)) {
            if found.is_some() {
                return Err(reject(RejectReason::HoleDepth)); // two minima: not a chain
            }
            found = Some(h);
        }
    }
    found.ok_or_else(|| reject(RejectReason::HoleDepth)) // no minimum: not a chain
}

/// A per-solid, per-side material label of one cell: `[A_above, A_below, B_above, B_below]`.
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
            // ★ Grazes fill arcs, so take their **union** rather than requiring them to agree.
            // Two faces of one solid meeting *at* this edge each fill one side, and together they
            // fill both — which is the same occupancy the diagram above calls `[T,T]`, and the same
            // thing a lone transversal means. Requiring agreement read that as a contradiction and
            // rejected it; it is the ordinary picture wherever a solid's edge lies in `W`, which is
            // exactly what a four-plane concurrency is made of.
            for &above in &grazes {
                mask[base + usize::from(!above)] = true;
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
    segs: &[MergedSeg],
    nesting: &Nesting,
    seed: Label,
) -> Result<Vec<Label>, BoolError> {
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
                let mask = edge_mask(&segs[he / 2].merged)?;
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
        let mask = edge_mask(&segs[he / 2].merged)?;
        if std::array::from_fn::<bool, 4, _>(|i| out[c][i] ^ mask[i]) != out[nb] {
            return Err(reject(RejectReason::LabelConflict)); // inconsistent propagation
        }
    }
    Ok(out)
}

/// The boolean keep predicate on one chamber's `(inA, inB)`.
fn keep(kind: BoolKind, in_a: bool, in_b: bool) -> bool {
    match kind {
        BoolKind::Fuse => in_a || in_b,
        BoolKind::Cut => in_a && !in_b,
        BoolKind::Common => in_a && in_b,
    }
}

/// Emit the result faces on plane class `wc` for a boolean `kind`. A `+1` cell is a face of the
/// result iff its two chambers disagree under `keep` (material on one side of W, void on the
/// other); its `-1` holes (`nesting.holes`) ride along as inner rings. The DCEL cell ring is
/// already CCW about `n_out(wc)` (extract_cells stored `winding == +1`) and a hole cell is CW
/// (`winding == -1`) — exactly the `LocalFace.inner` contract ("kept material on the loop's left"),
/// so both are emitted verbatim; `flip` alone carries the chamber and `assemble_fuse_cut` reverses
/// outer and inner together, making the result normal point out of the kept solid:
/// `flip = keep_above == (orient_sign(wc) > 0)`.
///
/// Only `+1` cells are hosts — a `-1` cell is either a hole (emitted as some host's inner ring) or
/// the void root — so `-1` cells are skipped, never emitted as their own face.
///
/// **Output contract:** every ring vertex is a `Node::Seam` triple, including triples that coincide
/// with an original A/B vertex (a cap corner); the weld table canonicalizes such a W-triple onto
/// the same result vertex, or `assemble_fuse_cut`'s manifold guard rejects. This brick emits
/// all-`Seam`; the reconciliation and assembly are later.
fn emit_faces(
    kind: BoolKind,
    labels: &[Label],
    cells: &[Cell],
    segs: &[MergedSeg],
    jd: &Judge<'_, PlaneGeom>,
    wc: usize,
    holes: &HashMap<usize, Vec<usize>>,
) -> Vec<LocalFace> {
    let planes = jd.planes;
    // `crate::boolean::Node` is lib.rs's arrangement node enum; the local `Node` (this module's
    // three-valued-scan struct) shadows it here.
    let ring_of = |cell: &Cell| -> Vec<crate::boolean::Node> {
        cell.half_edges
            .iter()
            .map(|&he| crate::boolean::Node::Seam(segs[he / 2].end[he % 2]))
            .collect()
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
        let inner: Vec<Vec<crate::boolean::Node>> = holes
            .get(&c)
            .map(|hs| hs.iter().map(|&h| ring_of(&cells[h])).collect())
            .unwrap_or_default();
        out.push(LocalFace {
            plane_idx: wc,
            loop_nodes: ring_of(cell),
            inner,
            flip,
        });
    }
    out
}

/// One subproblem the arrangement runs over: the tracer's input, and what that input's unbounded
/// cells are known to classify as.
///
/// A whole-model boolean is **one** region — every face, seed `[false; 4]`, which is true because
/// the arrangement then covers all of space and its unbounded cells reach infinity, where neither
/// solid is. Restricting the input to a region of space makes several; each is arranged
/// independently and their result faces concatenated, and each has to be *told* its seed because
/// its unbounded cells are an artifact of the restriction rather than a fact about the model.
pub(crate) struct Region {
    pub input: combinatorics::TraceInput,
    pub seed: Label,
}

/// How space is divided before arranging.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Subdivision {
    /// One region over everything — the arrangement as it has always been run.
    Whole,
    /// **V2a spike: partition the face set, without clipping.** Each face is dealt to one
    /// half-space of `coord[axis] = at` by which side it lies wholly on, and the two halves are
    /// arranged independently.
    ///
    /// It cuts nothing, so no new vertex is named and each half's outside is still genuinely void
    /// — the seed stays `[false; 4]` and no classification machinery is needed. What it does prove
    /// is the plumbing the subdivided engine rests on: that arranging a *subset* of the faces and
    /// concatenating the result faces gives the same solid. A face straddling the plane cannot be
    /// dealt, so a fixture must choose a plane that cuts none; it panics rather than guess.
    ///
    /// A non-cutting split can only ever separate *disconnected* pieces, so this is as far as the
    /// partition goes. [`Subdivision::ClipAt`] is where the question lives.
    #[cfg(test)]
    SplitFaces { axis: usize, at: f64 },
    /// **V2b spike: clip the faces to each half-space.** The split plane joins the class table as
    /// a synthetic class that names the cut vertices and is never itself arranged — there are no
    /// caps and no box walls, which is what keeps the halves from needing to be sewn back
    /// together.
    ///
    /// `at` must miss every model plane on `axis`. A split coinciding with one would merge into
    /// that class, and the artificial vertices would then be indistinguishable from real ones —
    /// `dissolve_straight_angles` could delete a genuine corner, silently.
    #[cfg(test)]
    ClipAt { axis: usize, at: f64 },
}

/// The split plane as a class of its own: axis-aligned, at `at`, its stored normal along `+axis`.
///
/// Wound so the witness triangle's right-hand normal **is** the stored normal, hence `frame_sign`
/// `+1`. The base frame is the one an unmoved model plane already carries (`chain_id 0`), which is
/// what a synthetic plane with no motion history is.
#[cfg(test)]
fn split_plane_class(model: &mut Model, axis: usize, at: f64) -> PlaneGeom {
    use nacre_math::Vector3;
    let e = |i: usize| {
        let mut v = [0.0; 3];
        v[i] = 1.0;
        v
    };
    // `u x v == e(axis)` for each axis, so the triangle winds to the stored normal.
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let mut o = [0.0; 3];
    o[axis] = at;
    let corner = |d: Option<usize>| {
        let mut p = o;
        if let Some(k) = d {
            p[k] = 1.0;
        }
        p
    };
    let tri_c = [corner(None), corner(Some(u)), corner(Some(v))];
    let tri = tri_c.map(Point3::from_array);
    let plane = Plane::from_point_normal(tri[0], Vector3::from_array(e(axis)))
        .expect("an axis-aligned plane has a unit normal");
    let surf = model.push_surface(
        nacre_geom::Surface::Plane(plane),
        nacre_topo::SurfaceDef::Constructed,
    );
    PlaneGeom {
        base: BaseFrame::none(),
        plane,
        surf,
        tri,
        tri_pt3: tri_c
            .map(|c| Pt3::exact(c).expect("a dyadic split plane is exactly representable")),
        rotated: false,
        frame_sign: 1,
    }
}

/// One surviving piece of a clipped face: its outer ring and its holes, in class-triple form.
#[cfg(test)]
pub(crate) type ClipPiece = (Vec<[usize; 3]>, Vec<Vec<[usize; 3]>>);

/// One face's rings clipped to a half-space.
#[cfg(test)]
#[derive(Debug, Default)]
pub(crate) struct Clipped {
    /// The surviving pieces — **there can be several**. Each is `(outer, holes)`.
    pub pieces: Vec<ClipPiece>,
    /// The directed edges the clip laid **on the split plane**, oriented as this face's boundary
    /// runs. Collected over a solid's faces they are the cap's edges — see `cap_rings`.
    pub cap_edges: Vec<([usize; 3], [usize; 3])>,
}

/// Clip a face's rings (outer + holes, plane class `fc`) to the `keep` side of class `split`.
///
/// **Every decision is an existing exact predicate, and no point is built.** A vertex's side is
/// [`combinatorics::side_of`]; a crossing's *name* is read straight off the two triples the edge
/// already carries — an edge between `t_i` and `t_{i+1}` rides the class they share besides `fc`,
/// so the new point is `{fc, wall, split}`; and the ordering along the clip line is
/// [`combinatorics::order_along`].
///
/// ★ **The crossings are paired in *line* order, not in ring order.** Sutherland-Hodgman pairs them
/// in ring order, which for a non-convex ring bridges two disjoint pieces along the clip line and
/// swallows the gap between them — measured on the U-prism's cap, where it turned two rectangles
/// into one eight-vertex ring covering the notch. That is a wrong answer, not a reject. Sorting the
/// chain ends along the line and pairing consecutive ones is the standard fix, and it is exact here
/// because `order_along` is.
///
/// A vertex *on* `split` is declined (`FourPlane`): a named vertex lying on the split plane means
/// four planes through one point, which this substrate cannot name — and the split position rule
/// exists to prevent it, so a firing means the rule was not applied.
#[cfg(test)]
fn clip_face(
    jd: &Judge<'_, PlaneGeom>,
    fc: usize,
    split: usize,
    keep: i8,
    rings: &[&[[usize; 3]]],
) -> Result<Clipped, BoolError> {
    // Which class a ring edge rides: the one its two endpoint triples share besides `fc`.
    let wall_of = |a: [usize; 3], b: [usize; 3]| -> Result<usize, BoolError> {
        let mut shared: Vec<usize> = a
            .iter()
            .copied()
            .filter(|&c| c != fc && b.contains(&c))
            .collect();
        shared.sort_unstable();
        shared.dedup();
        match shared[..] {
            [w] => Ok(w),
            _ => Err(reject(RejectReason::RingNaming)),
        }
    };

    let mut edges: Vec<([usize; 3], [usize; 3])> = Vec::new();
    // Chain ends on the clip line: `(node, is_start)`. A chain *starts* where the ring re-enters
    // the kept half and *ends* where it leaves.
    let mut ends: Vec<([usize; 3], bool)> = Vec::new();

    for ring in rings {
        let n = ring.len();
        if n < 3 {
            return Err(reject(RejectReason::DegenerateRing));
        }
        let mut side = Vec::with_capacity(n);
        for &t in ring.iter() {
            match combinatorics::side_of(jd, t, split) * keep {
                0 => return Err(reject(RejectReason::FourPlane)),
                s => side.push(s),
            }
        }
        if side.iter().all(|&s| s < 0) {
            continue; // wholly on the far side
        }
        for i in 0..n {
            let j = (i + 1) % n;
            let (a, b) = (ring[i], ring[j]);
            match (side[i] > 0, side[j] > 0) {
                (true, true) => edges.push((a, b)),
                (true, false) => {
                    let x = sorted3([fc, wall_of(a, b)?, split]);
                    edges.push((a, x));
                    ends.push((x, false));
                }
                (false, true) => {
                    let x = sorted3([fc, wall_of(a, b)?, split]);
                    edges.push((x, b));
                    ends.push((x, true));
                }
                (false, false) => {}
            }
        }
    }

    // ★ Pair the chain ends **along the line**. Every end is `{fc, w, split}`, so its position on
    // `L = fc ∩ split` is named by `w` alone and `order_along` compares two of them exactly.
    let mut cap_edges: Vec<([usize; 3], [usize; 3])> = Vec::new();
    if !ends.is_empty() {
        if ends.len() % 2 != 0 {
            return Err(reject(RejectReason::RingNaming)); // a chain with no partner
        }
        let handle = |t: [usize; 3]| -> Result<usize, BoolError> {
            t.into_iter()
                .find(|&c| c != fc && c != split)
                .ok_or_else(|| reject(RejectReason::RingNaming))
        };
        let mut order: Vec<(usize, ([usize; 3], bool))> = Vec::with_capacity(ends.len());
        for &e in &ends {
            order.push((handle(e.0)?, e));
        }
        // Two ends sharing a handle are one point reached twice — the ring pinches there and the
        // order is not defined. Honest reject rather than an arbitrary tie-break.
        let mut fail = None;
        order.sort_by(|x, y| {
            if x.0 == y.0 {
                fail = Some(());
                return std::cmp::Ordering::Equal;
            }
            match combinatorics::order_along(jd, fc, split, x.0, y.0) {
                -1 => std::cmp::Ordering::Less,
                1 => std::cmp::Ordering::Greater,
                _ => {
                    fail = Some(());
                    std::cmp::Ordering::Equal
                }
            }
        });
        if fail.is_some() {
            return Err(reject(RejectReason::CoincidentNodes));
        }
        for pair in order.chunks(2) {
            let [(_, (p, p_start)), (_, (q, q_start))] = pair else {
                unreachable!("even count, chunks of two")
            };
            // One of the two closes a chain and the other opens the next; the segment between them
            // is inside the face, and it runs from the closing end to the opening one.
            match (p_start, q_start) {
                (false, true) => cap_edges.push((*p, *q)),
                (true, false) => cap_edges.push((*q, *p)),
                _ => return Err(reject(RejectReason::RingNaming)),
            }
        }
        edges.extend(cap_edges.iter().copied());
    }

    // Thread the directed edges into cycles. Two edges leaving one node means the pieces meet at a
    // point and the cycles are not determined.
    let mut next: HashMap<[usize; 3], [usize; 3]> = HashMap::new();
    for &(a, b) in &edges {
        if next.insert(a, b).is_some() {
            return Err(reject(RejectReason::CoincidentNodes));
        }
    }
    let mut starts: Vec<[usize; 3]> = next.keys().copied().collect();
    starts.sort_unstable(); // deterministic cycle order
    let mut seen: HashSet<[usize; 3]> = HashSet::new();
    let mut cycles: Vec<Vec<[usize; 3]>> = Vec::new();
    for start in starts {
        if seen.contains(&start) {
            continue;
        }
        let mut cyc = vec![start];
        seen.insert(start);
        let mut cur = start;
        while let Some(&nx) = next.get(&cur) {
            if nx == start {
                break;
            }
            if !seen.insert(nx) {
                return Err(reject(RejectReason::RingNaming)); // threads into another cycle
            }
            cyc.push(nx);
            cur = nx;
        }
        if cyc.len() >= 3 {
            cycles.push(cyc);
        }
    }
    if cycles.is_empty() {
        return Ok(Clipped::default());
    }

    // Outer or hole, by winding **against the input's own outer ring** — so this never has to know
    // how the face is oriented against its plane class.
    let winding = |ring: &[[usize; 3]]| -> Result<i8, BoolError> {
        combinatorics::loop_winding(jd, fc, &combinatorics::ring_from_names(fc, ring)?)
    };
    let outward = winding(rings[0])?;
    let mut outers: Vec<Vec<[usize; 3]>> = Vec::new();
    let mut holes: Vec<Vec<[usize; 3]>> = Vec::new();
    for cyc in cycles {
        if winding(&cyc)? == outward {
            outers.push(cyc);
        } else {
            holes.push(cyc);
        }
    }
    let mut pieces: Vec<ClipPiece> = outers.into_iter().map(|o| (o, Vec::new())).collect();
    for hole in holes {
        let probe = hole[0];
        let mut owner = None;
        for (i, (outer, _)) in pieces.iter().enumerate() {
            let ring = combinatorics::ring_from_names(fc, outer)?;
            if combinatorics::point_in_ring(jd, fc, probe, &ring)? {
                if owner.is_some() {
                    return Err(reject(RejectReason::HoleDepth)); // nested deeper than one
                }
                owner = Some(i);
            }
        }
        pieces[owner.ok_or_else(|| reject(RejectReason::HoleRoots))?]
            .1
            .push(hole);
    }
    Ok(Clipped { pieces, cap_edges })
}

/// The regions one boolean is arranged over.
#[allow(clippy::too_many_arguments)]
pub(crate) fn regions_for(
    how: Subdivision,
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    faces: &[FaceInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &combinatorics::EdgeFaces,
    inc_b: &combinatorics::EdgeFaces,
    jd: &Judge<'_, PlaneGeom>,
    plane_ix: &[usize],
    split_class: Option<usize>,
) -> Vec<Region> {
    let _ = split_class;
    let input = combinatorics::trace_input(
        model,
        [(a, inc_a), (b, inc_b)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
    );
    match how {
        Subdivision::Whole => vec![Region {
            input,
            seed: [false; 4],
        }],
        #[cfg(test)]
        Subdivision::SplitFaces { axis, at } => {
            let mut halves = [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]];
            for (side, slots) in input.side_faces.iter().enumerate() {
                for &fp in slots {
                    let (lo, hi) = face_span(model, faces[fp].face, axis);
                    assert!(
                        hi <= at || lo >= at,
                        "the split plane {axis}={at} cuts face slot {fp} (span {lo}..{hi}); \
                         V2a partitions, it does not clip"
                    );
                    halves[usize::from(lo >= at)][side].push(fp);
                }
            }
            let [low, high] = halves;
            vec![
                Region {
                    input: combinatorics::TraceInput {
                        loops: input.loops.clone(),
                        side_faces: low,
                    },
                    seed: [false; 4],
                },
                Region {
                    input: combinatorics::TraceInput {
                        loops: input.loops,
                        side_faces: high,
                    },
                    seed: [false; 4],
                },
            ]
        }
        #[cfg(test)]
        Subdivision::ClipAt { .. } => {
            let split = split_class.expect("a clipping subdivision names its split class");
            // ★ **The seed stays `[false; 4]`, and that is a claim the test is here to falsify.**
            // A clipped ring closes along the split plane, so the restricted arrangement's
            // unbounded cells are outside the *clipped* operand too — where it genuinely has no
            // material. If that reasoning is wrong, it is wrong as a volume, not as a rounding
            // error.
            [1i8, -1]
                .into_iter()
                .map(|keep| Region {
                    input: clip_input(jd, &input, plane_ix, split, keep),
                    seed: [false; 4],
                })
                .collect()
        }
    }
}

/// Every face's loops, clipped to one side of `split`; faces left with nothing are dropped.
#[cfg(test)]
fn clip_input(
    jd: &Judge<'_, PlaneGeom>,
    input: &combinatorics::TraceInput,
    plane_ix: &[usize],
    split: usize,
    keep: i8,
) -> combinatorics::TraceInput {
    let mut loops = vec![combinatorics::FaceLoops::default(); input.loops.len()];
    let mut side_faces = [Vec::new(), Vec::new()];
    for (side, slots) in input.side_faces.iter().enumerate() {
        for &fp in slots {
            let fc = plane_ix[fp];
            let src = &input.loops[fp];
            let Some(outer) = src.outer.as_ref() else {
                // The whole-model run could not name this loop either; carry the failure through
                // so the region declines for the same reason rather than for a missing face.
                side_faces[side].push(fp);
                loops[fp] = src.clone();
                continue;
            };
            let mut rings: Vec<&[[usize; 3]]> = vec![outer];
            if let Some(hs) = src.holes.as_ref() {
                rings.extend(hs.iter().map(Vec::as_slice));
            }
            let clipped = clip_face(jd, fc, split, keep, &rings)
                .unwrap_or_else(|e| panic!("clip of face slot {fp} declined: {e:?}"));
            let mut it = clipped.pieces.into_iter();
            let Some((outer, holes)) = it.next() else {
                continue; // wholly on the other side
            };
            assert!(
                it.next().is_none(),
                "face slot {fp} clipped into several pieces; each needs its own `faces_tab` slot, \
                 which C3 adds — C1 only ships the clipper"
            );
            side_faces[side].push(fp);
            loops[fp] = combinatorics::FaceLoops {
                outer: Some(outer),
                holes: Some(holes),
            };
        }
    }
    combinatorics::TraceInput { loops, side_faces }
}

/// A face's extent along one axis, from its vertices' cached coordinates.
///
/// f64 is enough for what it decides — *which* half-space to deal a face to — because the split
/// plane is chosen away from every model coordinate. It decides nothing about the answer: a face
/// dealt to the wrong half would show up as a wrong solid, not as a rounding error.
#[cfg(test)]
fn face_span(model: &Model, fh: Handle<Face>, axis: usize) -> (f64, f64) {
    let f = model.faces.get(fh);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
        for he in &lp.half_edges {
            for &vh in model.edges.get(he.edge).bounds.iter().flatten() {
                let c = model.vertices.get(vh).point.as_array()[axis];
                lo = lo.min(c);
                hi = hi.max(c);
            }
        }
    }
    (lo, hi)
}

/// Every result face across all plane classes **of every region**, before assembly (the driver's
/// risky half, testable by face count without mutating the model). A declining class aborts the
/// whole boolean.
///
/// Regions are visited in order, classes within a region in order, so the faces are produced in a
/// sequence that is a function of the input — which is what keeps `assemble_fuse_cut`'s handle
/// minting replayable.
#[allow(clippy::too_many_arguments)]
fn trace_result_faces(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, PlaneGeom>,
    faces: &[FaceInfo],
    plane_ix: &[usize],
    n_a: usize,
    class_owner: &[Option<SolidSide>],
    reuse: crate::reuse::ClassReuse,
    regions: &[Region],
    n_classes: usize,
) -> Result<Vec<LocalFace>, BoolError> {
    let planes = jd.planes;
    let mut local_faces: Vec<LocalFace> = Vec::new();

    // **What each class contributes, for the classes the other operand cannot reach.** Everything
    // else stays `Arrange`, which is the whole engine as it was.
    //
    // The shortcut is a statement about a class over the *whole model* — "the other operand cannot
    // reach it, so its faces pass through unchanged" — which a region cannot make about its own
    // slice. Subdivided runs arrange everything; V5 retires the shortcut outright.
    let reuse = if regions.len() == 1 {
        reuse
    } else {
        crate::reuse::ClassReuse::Off
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
    let mut aliases = Aliases::default();
    let mut splits: Vec<Vec<MergedSeg>> = Vec::new();
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
        // Flattened `(region, class)`, region-major — so a single region is exactly the old index
        // space and the reject a decline raises is still the lowest-numbered class's.
        let round = crate::par::try_map_range(regions.len() * n_classes, |k| {
            let (trace_in, wc) = (&regions[k / n_classes].input, k % n_classes);
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
            let mut tr = trace_on_class(trace_in, wc, jd, faces, plane_ix);
            let mut local = snapshot.clone();
            local.absorb(&std::mem::take(&mut tr.aliases));
            // An incomplete trace ⇒ honest reject, naming what the tracer could not do and on
            // which operand face. A class can decline several faces; the first is the one
            // reported, and `try_map_range` picks the lowest-numbered class, which is the one
            // the sequential loop returned at.
            if let Some(&(fp, kind)) = tr.declined.first() {
                return Err(reject(decline_to_reject(kind, faces[fp].face)));
            }
            let merged = merge_coincident(&tr.segs, wc, &local);
            let split = split_at_crossings(jd, wc, &merged, &mut local)?;
            Ok((split, local))
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
        let (region, wc) = (&regions[k / n_classes], k % n_classes);
        let split = &splits[k];
        let arrange = |wc: usize| -> Result<Vec<LocalFace>, BoolError> {
            // ★ A region need not touch every class. One whose faces all lie elsewhere traces
            // nothing here, and an empty arrangement has no cells, hence no unbounded contour and
            // nothing to emit. `nest_cells` reads that as `HoleRoots` — "not a closed arrangement"
            // — which is the right verdict for the whole model, where every class carries at least
            // one face by construction and an empty trace really would be a defect, and the wrong
            // one for a region. Whole-model runs never reach this line.
            if split.is_empty() {
                return Ok(Vec::new());
            }
            let (cells, face_of) = extract_cells(jd, wc, split)?;
            let nesting = nest_cells(jd, wc, &cells, split)?;
            let labels = label_cells(&cells, &face_of, split, &nesting, region.seed)?;
            Ok(emit_faces(
                kind,
                &labels,
                &cells,
                split,
                jd,
                wc,
                &nesting.holes,
            ))
        };
        // **The plan decides, and only ever downwards.** A `PassThrough` that cannot name one of
        // its vertices falls back to arranging, so this can lose the shortcut but never the answer.
        let reused = match plans[wc] {
            crate::reuse::ClassPlan::Arrange => None,
            crate::reuse::ClassPlan::Empty => Some(Vec::new()),
            crate::reuse::ClassPlan::PassThrough(side) => {
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
            Some(f) => Ok(f),
            None => arrange(wc),
        }
    })?;
    for faces in per_class {
        local_faces.extend(faces);
    }
    Ok(local_faces)
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
    pub triple: [usize; 3],
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
        let mut names: Vec<[usize; 3]> = tr.touches.clone();
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
                        .unwrap_or_default();
                rings.extend(
                    combinatorics::hole_rings(model, fh, fp, inc, &jd, &plane_ix)
                        .unwrap_or_default()
                        .into_iter()
                        .flatten(),
                );
                // Only vertices the class would actually name: those lying on it (`side == 0`),
                // which is exactly the run condition the trace's rule fires under.
                names.extend(
                    rings
                        .into_iter()
                        .filter(|&t| combinatorics::side_of(&jd, sorted3(t), wc) == 0)
                        .map(sorted3),
                );
            }
        }
        names.sort_unstable();
        names.dedup();
        for t in names {
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
                    triple: t,
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
        };
        // Run the rest of the per-class pipeline, recording where it stops.
        audit.failed_at = if let Some(&(fp, kind)) = audit.declined.first() {
            Some(decline_to_reject(kind, faces_tab[fp].face))
        } else {
            let run = || -> Result<(), BoolError> {
                let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
                let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default())?;
                let (cells, face_of) = extract_cells(&jd, wc, &split)?;
                let nesting = nest_cells(&jd, wc, &cells, &split)?;
                let labels = label_cells(&cells, &face_of, &split, &nesting, [false; 4])?;
                let _ = emit_faces(kind, &labels, &cells, &split, &jd, wc, &nesting.holes);
                Ok(())
            };
            run().err().and_then(|e| match e {
                BoolError::Unsupported { reason } => Some(reason),
                // No live-set check runs inside the pipeline, so this arm is unreachable.
                BoolError::InputNotLive => None,
            })
        };
        out.push(audit);
    }
    Ok(out)
}

/// Drive the arrangement pipeline over **every** plane class and assemble the result solid — the
/// first end-to-end boolean from the trace engine. Isolated: nothing in production calls this.
///
/// Vertex welding is automatic: every result vertex is a sorted triple of three canon plane
/// classes (`canon3`), so a corner shared by three planes gets one identical `Node::Seam`
/// regardless of which plane was the cut class W — `assemble_fuse_cut` welds them to one vertex.
/// (This holds while every arrangement vertex is a 3-plane point; a 4-plane concurrency would name
/// it inconsistently and is out of scope here — axis-aligned boxes never produce one.)
///
/// A class that declines (holes, degenerate) aborts the whole boolean: skipping it would drop real
/// faces and silently produce a non-manifold or wrong-volume solid.
///
/// `pub(crate)` because the public `crate::boolean` delegates to it.
pub(crate) fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(Vec<Handle<Solid>>, Notes), BoolError> {
    boolean_over(model, kind, a, b, Subdivision::Whole)
}

/// [`boolean`], with the space division named — the seam the subdivided engine is built along.
pub(crate) fn boolean_over(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    how: Subdivision,
) -> Result<(Vec<Handle<Solid>>, Notes), BoolError> {
    // `mut` only for the synthetic split class, which is a test-only subdivision today.
    #[allow(unused_mut)]
    let PlaneSetup {
        planes: faces_tab,
        surf_ix,
        inc_a,
        inc_b,
        n_a,
        mut geom,
        plane_ix,
        mut class_owner,
        standard,
        notes,
        ..
    } = plane_index_setup(model, a, b)?;
    // **The classes the arrangement runs over.** A synthetic split plane is appended *after* them
    // and never arranged: it exists only to name the vertices clipping mints, and arranging it
    // would put a cap on every region boundary — the very thing the seed exists to avoid.
    let n_model_classes = geom.len();
    let split_class = match how {
        Subdivision::Whole => None,
        #[cfg(test)]
        Subdivision::SplitFaces { .. } => None,
        #[cfg(test)]
        Subdivision::ClipAt { axis, at } => {
            geom.push(split_plane_class(model, axis, at));
            class_owner.push(None);
            Some(geom.len() - 1)
        }
    };
    // The operation's judging, made once: the dense plane table, the standard it is held to, and
    // the collector. Everything below reaches predicates through this, so there is exactly one
    // place where "how this boolean judges" is decided.
    let jd = Judge::new(&geom, standard, &notes);
    // The plane classes are already decided at this point — `plane_index_setup` runs
    // `Judge::planes_coplanar` to build them — so a judgement that could not be made has already
    // shaped everything downstream. Say so before doing the work it would invalidate.
    undecided_reject(&notes)?;
    let run = |model: &mut Model| -> Result<Vec<Handle<Solid>>, BoolError> {
        let whole = regions_for(
            how,
            model,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            split_class,
        );
        let faces = trace_result_faces(
            model,
            kind,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &whole,
            n_model_classes,
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
            let whole = regions_for(
                how,
                model,
                a,
                b,
                &faces_tab,
                &surf_ix,
                &inc_a,
                &inc_b,
                &plain_jd,
                &plane_ix,
                split_class,
            );
            let plain = trace_result_faces(
                model,
                kind,
                a,
                b,
                &plain_jd,
                &faces_tab,
                &plane_ix,
                n_a,
                &class_owner,
                crate::reuse::ClassReuse::Off,
                &whole,
                n_model_classes,
            );
            match &plain {
                Ok(p) => assert_eq!(
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
        let faces = crate::boolean::unify_coplanar_faces(faces, &jd)?;

        // Build the SeamVertex weld table directly from the emitted triples (no `build_seam`: that is
        // raw-index and pierce-only). Reject rather than panic on a degenerate meet.
        let mut seam: Vec<SeamVertex> = Vec::new();
        let mut seen: HashMap<[usize; 3], ()> = HashMap::new();
        for f in &faces {
            for loop_ in std::iter::once(&f.loop_nodes).chain(f.inner.iter()) {
                for node in loop_ {
                    let crate::boolean::Node::Seam(t) = node;
                    if seen.insert(*t, ()).is_some() {
                        continue;
                    }
                    let point =
                        three_planes(&geom[t[0]].plane, &geom[t[1]].plane, &geom[t[2]].plane)
                            .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                    seam.push(SeamVertex {
                        point,
                        triple: *t,
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

        assemble_fuse_cut(model, a, b, &jd, &seam, &faces)
    };
    let out = run(model);
    // **The cause outranks the symptom, on both paths.** An undecided judgement has already been
    // read as a `0` by everything downstream, so whatever the engine then complains about — a
    // loop that will not orient, a trace that will not close — is a consequence being reported as
    // if it were the problem. That is the `LoopOrientMismatch`-hiding-precision-exhaustion trap,
    // and checking the evidence *before* returning the symptom is what keeps it shut.
    undecided_reject(&notes)?;
    out.map(|solids| (solids, notes))
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

    /// The engine entry with the evidence dropped — these tests assert geometry, and the report
    /// has its own tests. Shadows [`super::boolean`] so the call sites read as they always did.
    fn boolean(
        model: &mut Model,
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> Result<Vec<Handle<Solid>>, BoolError> {
        super::boolean(model, kind, a, b).map(|(solids, _)| solids)
    }

    /// Point of a canon triple, for asserting geometry by hand.
    fn pt(t: [usize; 3], jd: &Judge<'_, PlaneGeom>) -> [f64; 3] {
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
                            plane_ix[surf_ix[fh]] == c && face_on_z1(*fh, &surf_ix, &faces_tab)
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
        ]);
        let mut m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
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
        ]);
        let mut m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
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
        ]);
        let m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        let a = m.live_solids[0];
        let faces_tab = collect_planes(&m, a).unwrap();
        // One prism: no two faces are coplanar, so `plane_ix` is the identity and a face index and
        // its plane id coincide. Built through the real path anyway, so the test cannot drift.
        let canon = crate::planes::plane_classes(&crate::planes::test_judge(&faces_tab));
        let (planes, _plane_ix) = crate::planes::dense_planes(&faces_tab, &canon);
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
        let edges = [(fpy, 1i8), (fpy, -1), (fpx, 1), (fpx, -1), (fpd, 1)];

        let order = angular_order(&crate::planes::test_judge(&planes), w, &edges);
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
        let a0 = ang(edges[0]);
        let mut want: Vec<usize> = (0..edges.len()).collect();
        want.sort_by(|&i, &j| {
            let (ci, cj) = (
                (ang(edges[i]) - a0).rem_euclid(std::f64::consts::TAU),
                (ang(edges[j]) - a0).rem_euclid(std::f64::consts::TAU),
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
        let is_v = |t: [usize; 3]| {
            let p = pt(t, &jd);
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

        let (cells, face_of) = extract_cells(&jd, wc, &split).unwrap();
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
        let (cells, face_of) = extract_cells(&jd, wc, &split).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting, [false; 4]).unwrap();

        for (i, c) in cells.iter().enumerate() {
            if c.winding == -1 {
                assert_eq!(labels[i], [false; 4], "unbounded is void");
                continue;
            }
            // Centroid of the cell (convex here) via the average of its vertex points.
            let verts: Vec<[f64; 3]> = {
                let mut vs: Vec<[usize; 3]> = c
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
                        plane_ix[surf_ix[fh]] == c && face_on_z1(*fh, &surf_ix, &faces_tab)
                    })
                })
            })
            .expect("b's z=1 cap class");
        let tr = trace_on_class_of(
            &m, a, b, wc, &jd, &faces_tab, &surf_ix, &inc_a, &inc_b, &plane_ix,
        );
        let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
        let split = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()).unwrap();
        let (cells, face_of) = extract_cells(&jd, wc, &split).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting, [false; 4]).unwrap();

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
        let (cells, face_of) = extract_cells(&jd, wc, &split).unwrap();
        let nesting = nest_cells(&jd, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting, [false; 4]).unwrap();

        // Centroid of a face's ring (convex cells here).
        let centroid = |f: &LocalFace| -> [f64; 2] {
            let ps: Vec<[f64; 3]> = f
                .loop_nodes
                .iter()
                .map(|n| match n {
                    crate::boolean::Node::Seam(t) => pt(*t, &jd),
                })
                .collect();
            [
                ps.iter().map(|p| p[0]).sum::<f64>() / ps.len() as f64,
                ps.iter().map(|p| p[1]).sum::<f64>() / ps.len() as f64,
            ]
        };
        let near = |c: [f64; 2], x: f64, y: f64| (c[0] - x).abs() < 1e-9 && (c[1] - y).abs() < 1e-9;

        // Fuse: all 5 bounded cells (the plus cap).
        let fuse = emit_faces(
            BoolKind::Fuse,
            &labels,
            &cells,
            &split,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(fuse.len(), 5, "Fuse keeps the whole plus cap");
        // Cut a−b: exactly the two a-arms.
        let cut = emit_faces(
            BoolKind::Cut,
            &labels,
            &cells,
            &split,
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
        let common = emit_faces(
            BoolKind::Common,
            &labels,
            &cells,
            &split,
            &jd,
            wc,
            &nesting.holes,
        );
        assert_eq!(common.len(), 1, "Common keeps the center");
        assert!(near(centroid(&common[0]), 1.5, 1.5), "center at (1.5,1.5)");

        // Each emitted loop winds +1 (CCW about n_out(wc)) and has ≥3 distinct nodes.
        for f in fuse.iter().chain(&cut).chain(&common) {
            let ring: Vec<[usize; 3]> = f
                .loop_nodes
                .iter()
                .map(|crate::boolean::Node::Seam(t)| *t)
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
                .loop_nodes
                .iter()
                .map(|crate::boolean::Node::Seam(t)| *t)
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
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let whole = regions_for(
            Subdivision::Whole,
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            None,
        );
        let faces = trace_result_faces(
            &m,
            BoolKind::Fuse,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &whole,
            planes.len(),
        )
        .unwrap();
        assert_eq!(faces.len(), 10, "z=0 + z=2 + 4 walls×2");

        // Every undirected edge (sorted triple pair) is used exactly twice — a closed shell.
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            let ns: Vec<[usize; 3]> = f
                .loop_nodes
                .iter()
                .map(|n| match n {
                    crate::boolean::Node::Seam(t) => *t,
                })
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
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let whole = regions_for(
            Subdivision::Whole,
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            None,
        );
        let faces = trace_result_faces(
            &m,
            BoolKind::Cut,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &whole,
            planes.len(),
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
        let triples = |ns: &[crate::boolean::Node]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|n| match n {
                    crate::boolean::Node::Seam(t) => *t,
                })
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in std::iter::once(&f.loop_nodes).chain(f.inner.iter()) {
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
        ]);
        let mut m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
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
            ..
        } = plane_index_setup(&m, u, slab).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let whole = regions_for(
            Subdivision::Whole,
            &m,
            u,
            slab,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            None,
        );
        let faces = trace_result_faces(
            &m,
            BoolKind::Fuse,
            u,
            slab,
            &jd,
            &faces_tab,
            &plane_ix,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &whole,
            planes.len(),
        )
        .unwrap();

        // Exactly one face carries two inner rings: the y=2 slab annulus, holed by both prongs.
        assert_eq!(
            faces.iter().filter(|f| f.inner.len() == 2).count(),
            1,
            "the slab face on y=2 has two holes (the U's two prongs)"
        );

        // Every undirected edge across outer + inner rings is used exactly twice (closed shell).
        let triples = |ns: &[crate::boolean::Node]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|n| match n {
                    crate::boolean::Node::Seam(t) => *t,
                })
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in std::iter::once(&f.loop_nodes).chain(f.inner.iter()) {
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
    // solid `Origin::Moved`, routing every predicate to the exact frame3 backend. `rot30` (lib.rs)
    // is in a sibling test module and unreachable here, so the isometries are built inline.

    /// 30° about `axis` through (1,1,0) — non-90°, so `Origin::Moved` (exact frame3 path).
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
    /// agreement had to be pinned on a four-plane variant that still stops there. `half_z = 0.5`
    /// at `120°` is one; if a later capability makes it build, the fix is to re-run that sweep and
    /// take whatever still declines — not to weaken the assertion.
    #[test]
    fn the_audit_reports_the_same_reject_as_the_boolean() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let build = || -> (Model, Handle<Solid>, Handle<Solid>) {
            let mut m = Model::new();
            let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let block = m.add_cuboid(
                Point3::from_array([0.5, 0.0, 1.0]),
                Point3::from_array([1.0, 1.0, 2.0]),
            );
            m.rebuild_adjacency();
            let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("block fuses on")[0];
            m.rebuild_adjacency();
            let bar = m.add_cuboid(
                Point3::from_array([0.3, -0.5, 0.5]),
                Point3::from_array([0.7, 1.5, 1.5]),
            );
            m.rebuild_adjacency();
            let bar = transform(
                &mut m,
                bar,
                &Isometry::rotation(Rotation {
                    axis: Axis::Y,
                    point: [Rat::new(1, 2).unwrap(), Rat::from_int(0), Rat::from_int(1)],
                    angle: Angle::from_deg(Rat::from_int(120)).unwrap(),
                }),
            )
            .unwrap();
            m.rebuild_adjacency();
            (m, target, bar)
        };

        let (mut m, target, bar) = build();
        let err = boolean(&mut m, BoolKind::Cut, target, bar).unwrap_err();
        let BoolError::Unsupported { reason } = err else {
            panic!("expected an Unsupported rejection, got {err:?}");
        };

        let (m, target, bar) = build();
        let audits = frame_audit(&m, BoolKind::Cut, target, bar).unwrap();
        let failed: Vec<RejectReason> = audits.iter().filter_map(|a| a.failed_at).collect();
        assert!(
            failed.contains(&reason),
            "the audit must report the boolean's reject ({reason:?}), got {failed:?}"
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
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&planes, standard, &notes);
        let whole = regions_for(
            Subdivision::Whole,
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            None,
        );
        let faces = trace_result_faces(
            &m,
            BoolKind::Cut,
            a,
            b,
            &jd,
            &faces_tab,
            &plane_ix,
            n_a,
            &class_owner,
            crate::reuse::ClassReuse::Proved,
            &whole,
            planes.len(),
        )
        .unwrap();
        assert_eq!(faces.len(), 10, "rotated arrangement keeps 10 faces");
        assert_eq!(
            faces.iter().filter(|f| !f.inner.is_empty()).count(),
            2,
            "the two annular caps survive rotation"
        );
        let triples = |ns: &[crate::boolean::Node]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|n| match n {
                    crate::boolean::Node::Seam(t) => *t,
                })
                .collect()
        };
        let mut count: HashMap<([usize; 3], [usize; 3]), usize> = HashMap::new();
        for f in &faces {
            for ring in std::iter::once(&f.loop_nodes).chain(f.inner.iter()) {
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
        faces: &[FaceInfo],
        plane_ix: &[usize],
    ) -> usize {
        let n_planes = plane_ix.iter().copied().max().map_or(0, |m| m + 1);
        (0..n_planes)
            .find(|&c| {
                let seats =
                    |s: Handle<Solid>| {
                        solid_shell_handles(m, s).into_iter().any(|sh| {
                            m.shells.get(sh).faces.iter().any(|fh| {
                                plane_ix[surf_ix[fh]] == c && face_on_z1(*fh, surf_ix, faces)
                            })
                        })
                    };
                seats(a) && seats(b)
            })
            .expect("a shared z=1 cap class")
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
        faces: &[FaceInfo],
    ) -> bool {
        let p = &faces[surf_ix[&fh]];
        // A cap in the z=1 plane: all three defining points at z=1.
        p.tri.iter().all(|q| (q.as_array()[2] - 1.0).abs() < 1e-12)
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

    // ── V2a: the arrangement, run over two regions instead of one ─────────────────────────────
    //
    // The engine's work unit is an unbounded plane, so a plane class's arrangement has always
    // needed the whole model. These prove the first half of undoing that: that the *same* engine,
    // shown a **subset** of the faces and asked for its result faces, composes with a second run
    // shown the rest.
    //
    // The split cuts nothing, which is what keeps this stage free of new machinery: no vertex is
    // named that was not named before, and each half's outside is genuinely void, so the seed
    // stays the constant `label_cells` always used. The seed itself is V2b's question and the
    // clipping is V3's; what is on trial here is only the composition.

    /// A solid's faces, and the coordinates it was built from, in a form two runs can be compared
    /// by. Sorted by bits, so it depends on the geometry and not on handle order.
    fn shape_digest(m: &Model, solids: &[Handle<Solid>]) -> Vec<(usize, f64, Vec<[u64; 3]>)> {
        solids
            .iter()
            .map(|&s| {
                let src = m.solids.get(s).clone();
                let mut bits: Vec<[u64; 3]> = Vec::new();
                let mut faces = 0usize;
                for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
                    for &fh in &m.shells.get(sh).faces {
                        faces += 1;
                        let face = m.faces.get(fh);
                        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                            for he in &lp.half_edges {
                                for &vh in m.edges.get(he.edge).bounds.iter().flatten() {
                                    let p = m.vertices.get(vh).point.as_array();
                                    bits.push(p.map(f64::to_bits));
                                }
                            }
                        }
                    }
                }
                bits.sort_unstable();
                bits.dedup();
                let v = nacre_props::mass_props(m, s).expect("props").volume;
                (faces, v, bits)
            })
            .collect()
    }

    /// Two boxes apart along x that **share four plane classes** (both span y,z ∈ [0,1]), so a
    /// split between them gives each region a *part* of those classes rather than a class of its
    /// own. That is the case the whole design turns on: a class's arrangement restricted to a
    /// region of space.
    fn two_bars(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        (a, b)
    }

    /// ★ The check is not vacuous: the split really does deal the faces into two non-empty halves,
    /// and every face lands in exactly one of them. Without this, "the answers agree" would also
    /// pass if `SplitFaces` quietly produced one whole region.
    #[test]
    fn splitting_the_face_set_deals_every_face_to_exactly_one_half() {
        let mut m = Model::new();
        let (a, b) = two_bars(&mut m);
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
        } = plane_index_setup(&m, a, b).unwrap();
        let jd = Judge::new(&geom, standard, &notes);
        let regions = regions_for(
            Subdivision::SplitFaces { axis: 0, at: 1.5 },
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            None,
        );
        assert_eq!(regions.len(), 2);
        let count = |r: &Region| r.input.side_faces.iter().map(Vec::len).sum::<usize>();
        assert_eq!(count(&regions[0]), 6, "a's six faces are below x=1.5");
        assert_eq!(count(&regions[1]), 6, "b's six faces are above x=1.5");
        // A shared class is *split* between the regions rather than owned by one of them — the
        // property that makes this fixture worth more than two unrelated solids.
        let classes = |r: &Region| -> Vec<usize> {
            let mut c: Vec<usize> = r
                .input
                .side_faces
                .iter()
                .flatten()
                .map(|&fp| plane_ix[fp])
                .collect();
            c.sort_unstable();
            c.dedup();
            c
        };
        let shared: Vec<usize> = classes(&regions[0])
            .into_iter()
            .filter(|c| classes(&regions[1]).contains(c))
            .collect();
        assert_eq!(shared.len(), 4, "z=0, z=1, y=0, y=1 are in both halves");
    }

    /// **The composition itself.** Every kind, answered twice — once by the engine over the whole
    /// model, once over two regions — and required to agree down to the coordinate bits.
    #[test]
    fn arranging_two_regions_gives_the_same_solid_as_arranging_one() {
        // `1.5` separates the two bars. `-5.0` and `9.0` put *everything* on one side, so the
        // other region is empty — the degenerate end of the same mechanism, and the one that
        // says an empty region contributes nothing rather than declining.
        for at in [1.5, -5.0, 9.0] {
            for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                one_split_agrees_with_the_whole(kind, at);
            }
        }
    }

    fn one_split_agrees_with_the_whole(kind: BoolKind, at: f64) {
        {
            let mut m1 = Model::new();
            let (a1, b1) = two_bars(&mut m1);
            let whole = super::boolean(&mut m1, kind, a1, b1).map(|(s, _)| s);

            let mut m2 = Model::new();
            let (a2, b2) = two_bars(&mut m2);
            let split = boolean_over(
                &mut m2,
                kind,
                a2,
                b2,
                Subdivision::SplitFaces { axis: 0, at },
            )
            .map(|(s, _)| s);

            match (&whole, &split) {
                (Ok(w), Ok(s)) => {
                    m1.rebuild_adjacency();
                    m2.rebuild_adjacency();
                    assert_eq!(
                        shape_digest(&m1, w),
                        shape_digest(&m2, s),
                        "{kind:?} at {at}: two regions disagreed with one"
                    );
                }
                (Err(w), Err(s)) => assert_eq!(
                    format!("{w:?}"),
                    format!("{s:?}"),
                    "{kind:?} at {at}: the two runs rejected differently"
                ),
                _ => panic!(
                    "{kind:?} at {at}: one run answered and the other did not: {whole:?} vs {split:?}"
                ),
            }
        }
    }

    // ── V2b: the same, with the faces actually cut ────────────────────────────────────────────
    //
    // The substantive stage. The split plane joins the class table as a synthetic class, every
    // face's rings are clipped to each half-space in **triple space** — a crossing vertex is
    // named `{face plane, wall, split}`, read straight off the two triples the edge already
    // carries — and the two halves are arranged and concatenated as in V2a.
    //
    // There are no caps and no box walls. The claim under test is that none are needed: a clipped
    // ring closes along the split plane, so a region's unbounded cells lie outside the *clipped*
    // operand, where it has no material, and the seed stays the constant `label_cells` always
    // used. If that is wrong it is wrong as a volume.

    /// Two cubes that overlap, so no plane separates them and every split cuts real faces.
    fn overlapping_cubes(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        m.rebuild_adjacency();
        (a, b)
    }

    /// ★ Not vacuous: the split really cuts. Every face of both cubes spans `x = 0.75`, so each
    /// half keeps a *clipped* copy of all twelve, with vertices that name the split class.
    #[test]
    fn clipping_at_a_plane_through_material_really_cuts() {
        let mut m = Model::new();
        let (a, b) = overlapping_cubes(&mut m);
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
        } = plane_index_setup(&m, a, b).unwrap();
        let mut geom = geom;
        let n_model_classes = geom.len();
        geom.push(split_plane_class(&mut m, 0, 0.75));
        let split = geom.len() - 1;
        let jd = Judge::new(&geom, standard, &notes);
        let regions = regions_for(
            Subdivision::ClipAt { axis: 0, at: 0.75 },
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            Some(split),
        );
        assert_eq!(regions.len(), 2);
        for (i, r) in regions.iter().enumerate() {
            let slots: Vec<usize> = r.input.side_faces.iter().flatten().copied().collect();
            // The two `x` walls of each cube lie wholly on one side; the other ten straddle.
            assert!(
                (10..=11).contains(&slots.len()),
                "half {i} kept {} faces",
                slots.len()
            );
            let cut = slots
                .iter()
                .filter(|&&fp| {
                    r.input.loops[fp]
                        .outer
                        .as_ref()
                        .is_some_and(|o| o.iter().any(|t| t.contains(&split)))
                })
                .count();
            assert_eq!(cut, 8, "half {i}: eight faces carry split-plane vertices");
        }
        assert!(n_model_classes < geom.len(), "the split class is appended");
    }

    /// ★ **Why clipping alone does not work, as a property rather than an anecdote.**
    ///
    /// A face **seated** on class `W` closes by itself: its own ring was clipped, and
    /// Sutherland-Hodgman put an edge on the split plane to close it. A solid that merely
    /// **straddles** `W` does not: its footprint there is bounded by the chords of its walls, and
    /// the wall at the split plane is the cap that "no caps" declined to make. So its trace is an
    /// open chain, the cells do not close, and `label_cells` reports the contradiction it should.
    ///
    /// Measured on every class of the clipped fixture that has a straddling operand: **exactly two
    /// odd-degree nodes, and both name the split class** — the two ends of the missing cap chord.
    /// That is the shape of the repair, not just of the failure: one segment on `wall = split`
    /// joining them, which is what tracing the clipped solid's cap face would have produced.
    #[test]
    fn clipping_leaves_a_straddling_operand_open() {
        let mut m = Model::new();
        let (a, b) = overlapping_cubes(&mut m);
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
        } = plane_index_setup(&m, a, b).unwrap();
        let mut geom = geom;
        let n_model = geom.len();
        geom.push(split_plane_class(&mut m, 0, 0.75));
        let split = geom.len() - 1;
        let jd = Judge::new(&geom, standard, &notes);
        let regions = regions_for(
            Subdivision::ClipAt { axis: 0, at: 0.75 },
            &m,
            a,
            b,
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &jd,
            &plane_ix,
            Some(split),
        );
        let mut open_classes = 0usize;
        for (ri, r) in regions.iter().enumerate() {
            for wc in 0..n_model {
                let tr = trace_on_class(&r.input, wc, &jd, &faces_tab, &plane_ix);
                if tr.segs.is_empty() {
                    continue;
                }
                let merged = merge_coincident(&tr.segs, wc, &Aliases::default());
                let Ok(sp) = split_at_crossings(&jd, wc, &merged, &mut Aliases::default()) else {
                    continue;
                };
                let mut deg: HashMap<[usize; 3], usize> = HashMap::new();
                for s in &sp {
                    for e in s.end {
                        *deg.entry(sorted3(e)).or_default() += 1;
                    }
                }
                let odd: Vec<[usize; 3]> = deg
                    .iter()
                    .filter(|&(_, &d)| d % 2 == 1)
                    .map(|(n, _)| *n)
                    .collect();
                if odd.is_empty() {
                    continue;
                }
                open_classes += 1;
                assert_eq!(
                    odd.len(),
                    2,
                    "region {ri} class {wc}: an open chain has two ends"
                );
                for n in &odd {
                    assert!(
                        n.contains(&split),
                        "region {ri} class {wc}: the open end {n:?} is not on the split plane, \
                         so the missing edge is not the cap chord"
                    );
                }
            }
        }
        // ★ Not vacuous: the failure has to actually occur, or the assertions above prove nothing.
        assert!(
            open_classes >= 4,
            "only {open_classes} classes were left open; the fixture stopped exercising the case"
        );
    }

    // ── C1: the half-space clipper ────────────────────────────────────────────────────────────

    /// A U-prism and a far-away box, so the plane table has two operands and the U's caps are
    /// eight-vertex **non-convex** rings.
    fn u_prism_and_far(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
        let u_profile = Profile2d::polygon(vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([3.0, 0.0]),
            Point2::from_array([3.0, 2.3]),
            Point2::from_array([2.0, 2.3]),
            Point2::from_array([2.0, 1.0]),
            Point2::from_array([1.0, 1.0]),
            Point2::from_array([1.0, 2.0]),
            Point2::from_array([0.0, 2.0]),
        ]);
        let out = apply(
            m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: u_profile,
                dist: 1.0,
            },
        )
        .expect("extrude");
        let OpOutput::Extrude { solid: u, .. } = out else {
            unreachable!()
        };
        let far = m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        (u, far)
    }

    /// The clipper's answer for one of the U's caps, as rounded 2D vertex sets per piece.
    fn clip_u_cap(at: f64, keep: i8) -> Vec<(Vec<[i64; 2]>, usize)> {
        let mut m = Model::new();
        let (u, far) = u_prism_and_far(&mut m);
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
        } = plane_index_setup(&m, u, far).unwrap();
        let mut geom = geom;
        geom.push(split_plane_class(&mut m, 1, at));
        let split = geom.len() - 1;
        let jd = Judge::new(&geom, standard, &notes);
        let input = combinatorics::trace_input(
            &m,
            [(u, &inc_a), (far, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
        );
        let mut out = Vec::new();
        for &fp in &input.side_faces[0] {
            let ring = input.loops[fp].outer.as_ref().unwrap();
            if ring.len() != 8 {
                continue; // only the two U caps are the eight-vertex non-convex rings
            }
            let fc = plane_ix[fp];
            let clipped = clip_face(&jd, fc, split, keep, &[ring]).expect("clip");
            for (outer, holes) in &clipped.pieces {
                assert!(holes.is_empty(), "a U cap clip has no holes");
                let mut xy: Vec<[i64; 2]> = outer
                    .iter()
                    .map(|&t| {
                        let p = pt(t, &jd);
                        [
                            (p[0] * 1000.0).round() as i64,
                            (p[1] * 1000.0).round() as i64,
                        ]
                    })
                    .collect();
                xy.sort_unstable();
                out.push((xy, clipped.cap_edges.len()));
            }
            // Both caps must answer identically; one sample is enough for the assertions below.
            break;
        }
        out
    }

    /// ★ **The refutation, turned into a gate.** A line through the U's notch crosses its cap four
    /// times, and the answer is two rectangles. Sutherland-Hodgman pairs the crossings in *ring*
    /// order and returns one eight-vertex ring that swallows the notch — a wrong answer, not a
    /// reject (measured 2026-07-29). Pairing along the **line** is what fixes it.
    #[test]
    fn clipping_a_non_convex_ring_gives_the_pieces_not_a_bridge() {
        let pieces = clip_u_cap(1.5, 1);
        assert_eq!(pieces.len(), 2, "the notch separates the two prongs");
        let sets: Vec<Vec<[i64; 2]>> = pieces.iter().map(|(v, _)| v.clone()).collect();
        let left = vec![[0, 1500], [0, 2000], [1000, 1500], [1000, 2000]];
        let right = vec![[2000, 1500], [2000, 2300], [3000, 1500], [3000, 2300]];
        assert!(sets.contains(&left), "left prong: {sets:?}");
        assert!(sets.contains(&right), "right prong: {sets:?}");
        // ★ Negative control: the answer Sutherland-Hodgman gave — one ring over all eight
        // vertices, covering the notch — must not be reachable.
        let bridged: Vec<[i64; 2]> = {
            let mut v = left.clone();
            v.extend(right.iter().copied());
            v.sort_unstable();
            v
        };
        assert!(
            !sets.contains(&bridged),
            "the bridged eight-vertex ring is the wrong answer this test exists to exclude"
        );
        // Two pieces are closed by two segments on the split plane — the cap's edges.
        assert_eq!(pieces[0].1, 2, "one cap edge per piece");
    }

    /// The other side of the same cut is **connected** (the U's bar joins the prongs below the
    /// notch), so the same machinery must give exactly one piece — while still laying **two**
    /// segments on the split plane, because between `x=1` and `x=2` the boundary runs along the
    /// notch floor at `y=1` and not along the cut. One piece and two cap edges is the pair of
    /// facts that says the pairing is by position on the line and not by piece.
    #[test]
    fn the_connected_side_of_the_same_cut_stays_one_piece() {
        let pieces = clip_u_cap(1.5, -1);
        assert_eq!(pieces.len(), 1, "below y=1.5 the U is connected");
        assert_eq!(
            pieces[0].0,
            vec![
                [0, 0],
                [0, 1500],
                [1000, 1000],
                [1000, 1500],
                [2000, 1000],
                [2000, 1500],
                [3000, 0],
                [3000, 1500],
            ]
        );
        assert_eq!(
            pieces[0].1, 2,
            "two segments on the split plane, one per prong foot"
        );
    }

    /// A cut that misses the ring entirely keeps it whole on one side and drops it on the other,
    /// and lays no edge on the split plane either way.
    #[test]
    fn a_cut_that_misses_keeps_or_drops_the_whole_ring() {
        let above = clip_u_cap(-1.0, 1);
        assert_eq!(above.len(), 1, "wholly above y=-1");
        assert_eq!(above[0].0.len(), 8, "untouched");
        assert_eq!(above[0].1, 0, "no cap edge");
        assert!(clip_u_cap(-1.0, -1).is_empty(), "nothing below y=-1");
    }

    /// A slab with a square hole bored through it, and a far box for the second operand. The
    /// through-hole leaves the slab's two caps with an inner ring, which is what makes this the
    /// holed-face fixture.
    fn holed_slab_and_far(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
        let slab = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([4.0, 4.0, 1.0]),
        );
        let drill = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        );
        m.rebuild_adjacency();
        let holed = boolean(m, BoolKind::Cut, slab, drill).expect("bore")[0];
        m.rebuild_adjacency();
        let far = m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([11.0, 11.0, 11.0]),
        );
        m.rebuild_adjacency();
        (holed, far)
    }

    /// Clip the holed slab's first holed cap at `x = at`, as `(outer vertex count, hole count)`
    /// per piece plus the number of edges laid on the split plane.
    fn clip_holed_cap(at: f64, keep: i8) -> (Vec<(usize, usize)>, usize) {
        let mut m = Model::new();
        let (holed, far) = holed_slab_and_far(&mut m);
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
        } = plane_index_setup(&m, holed, far).unwrap();
        let mut geom = geom;
        geom.push(split_plane_class(&mut m, 0, at));
        let split = geom.len() - 1;
        let jd = Judge::new(&geom, standard, &notes);
        let input = combinatorics::trace_input(
            &m,
            [(holed, &inc_a), (far, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
        );
        for &fp in &input.side_faces[0] {
            let fl = &input.loops[fp];
            let (Some(outer), Some(holes)) = (fl.outer.as_ref(), fl.holes.as_ref()) else {
                continue;
            };
            if holes.len() != 1 {
                continue; // the two caps are the only holed faces
            }
            let mut rings: Vec<&[[usize; 3]]> = vec![outer];
            rings.extend(holes.iter().map(Vec::as_slice));
            let clipped = clip_face(&jd, plane_ix[fp], split, keep, &rings).expect("clip");
            return (
                clipped
                    .pieces
                    .iter()
                    .map(|(o, h)| (o.len(), h.len()))
                    .collect(),
                clipped.cap_edges.len(),
            );
        }
        panic!("the bored slab has a holed cap");
    }

    /// A cut that misses the hole leaves it whole on one side and gives a plain rectangle on the
    /// other — the hole must not be dropped, duplicated, or promoted to an outer.
    #[test]
    fn clipping_beside_a_hole_keeps_it_on_the_side_it_belongs_to() {
        assert_eq!(clip_holed_cap(0.5, -1), (vec![(4, 0)], 1), "x<0.5: no hole");
        assert_eq!(clip_holed_cap(0.5, 1), (vec![(4, 1)], 1), "x>0.5: the hole");
    }

    /// ★ A cut **through** the hole: the inner ring stops being a hole and becomes part of the
    /// outer boundary. This is the case a per-ring clipper cannot express — it has to see the
    /// outer and the hole as one crossing set, which is why `clip_face` takes all the rings at
    /// once and pairs across them.
    #[test]
    fn a_cut_through_a_hole_merges_it_into_the_outer_ring() {
        let (pieces, cap_edges) = clip_holed_cap(2.0, -1);
        assert_eq!(pieces, vec![(8, 0)], "a C shape: one ring, no hole");
        assert_eq!(cap_edges, 2, "the cut meets material twice along x=2");
    }

    /// ★ **A conservation identity, stronger than a piece count.** Clipping a ring set by a plane
    /// splits every crossed edge in two and gives each half a copy of the crossing vertex, so
    ///
    /// ```text
    ///   verts(keep +1) + verts(keep -1) == verts(original) + 2 * crossings
    /// ```
    ///
    /// and `crossings == 2 * cap_edges` on either side. It holds whatever the piece structure is,
    /// which is what makes it usable on a rotated fixture where the answer's shape is not obvious
    /// by hand.
    fn clip_conserves_vertices(
        jd: &Judge<'_, PlaneGeom>,
        fc: usize,
        split: usize,
        rings: &[&[[usize; 3]]],
    ) -> bool {
        let original: usize = rings.iter().map(|r| r.len()).sum();
        let mut total = 0usize;
        let mut cap = [0usize; 2];
        for (i, keep) in [1i8, -1].into_iter().enumerate() {
            let c = clip_face(jd, fc, split, keep, rings).expect("clip");
            total += c
                .pieces
                .iter()
                .map(|(o, h)| o.len() + h.iter().map(Vec::len).sum::<usize>())
                .sum::<usize>();
            cap[i] = c.cap_edges.len();
        }
        assert_eq!(cap[0], cap[1], "both sides are closed by the same segments");
        total == original + 4 * cap[0]
    }

    /// The clipper on a **rotated** solid: the split plane is axis-aligned but the model's planes
    /// are not, so every `side_of` and `order_along` leaves the exact-shortcut path and is answered
    /// by the toleranced kernel. Nothing may decline, and the conservation identity must hold at
    /// every cut position that meets the solid.
    #[test]
    fn clipping_survives_a_rotated_model() {
        let mut m = Model::new();
        let (u, far) = u_prism_and_far(&mut m);
        let u = tilt(&mut m, u, &[nacre_scalar::Axis::Z]);
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
        } = plane_index_setup(&m, u, far).unwrap();
        // All four cut positions join the table at once — `PlaneGeom` is not `Clone`, and one
        // table is what production has anyway.
        let cuts = [0.5f64, 1.25, 1.75, 2.5];
        let mut geom = geom;
        let n_model = geom.len();
        for &at in &cuts {
            geom.push(split_plane_class(&mut m, 1, at));
        }
        let jd = Judge::new(&geom, standard, &notes);
        let input = combinatorics::trace_input(
            &m,
            [(u, &inc_a), (far, &inc_b)],
            &surf_ix,
            faces_tab.len(),
            &jd,
            &plane_ix,
        );
        let mut fired = 0usize;
        for (k, at) in cuts.into_iter().enumerate() {
            let split = n_model + k;
            for &fp in &input.side_faces[0] {
                let Some(outer) = input.loops[fp].outer.as_ref() else {
                    continue;
                };
                let rings: Vec<&[[usize; 3]]> = vec![outer];
                let sides: Vec<i8> = outer
                    .iter()
                    .map(|&t| combinatorics::side_of(&jd, t, split))
                    .collect();
                if sides.iter().all(|&s| s >= 0) || sides.iter().all(|&s| s <= 0) {
                    continue; // this face is not cut here
                }
                fired += 1;
                assert!(
                    clip_conserves_vertices(&jd, plane_ix[fp], split, &rings),
                    "rotated clip at y={at} lost or invented vertices on face slot {fp}"
                );
            }
        }
        // ★ Not vacuous: some face has to actually be cut, or the loop asserted nothing.
        assert!(fired >= 4, "only {fired} rotated faces were cut");
    }
}
