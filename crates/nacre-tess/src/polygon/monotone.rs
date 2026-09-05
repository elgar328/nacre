//! **y-monotone decomposition of a polygon with holes** — the plane sweep of
//! de Berg, Cheong, van Kreveld & Overmars, *Computational Geometry* §3.2.
//!
//! The point of running a sweep instead of clipping ears is that it **never merges
//! rings**. Ear clipping needs one ring, so a hole has to be bridged into the outer
//! boundary by a zero-width slit — and that repeats a vertex, which makes the ring
//! non-simple, which is exactly the hypothesis Meisters' two-ears theorem needs. A
//! bridged ring can have *no ear at all*, and the kernel measured that: 27.6% of
//! random rectilinear faces with two to four rectangular holes could not be meshed.
//!
//! The sweep has no such gap. Hole edges are ordinary edges; a hole's topmost vertex
//! is a *split* vertex and its bottommost a *merge* vertex, and the diagonals those
//! two cases add are what join the hole to the rest. Every diagonal runs between two
//! **distinct** vertices, so no piece it produces ever repeats one.
//!
//! Everything decided here is a sign: which of two vertices comes first ([`lex_less`],
//! exact because the projected coordinates are copied `f64`s rather than computed
//! ones) and which side of a line a point falls on (`orient2d`, Shewchuk adaptive).
//! There is no tolerance anywhere in this file.

use super::sos::{Twins, lex_less_idx, side_idx};
use super::{P2, lex_less, strictly_between};
use crate::TessError;
use nacre_predicates::orient2d;

/// Which side of the directed line `a → b` the point `c` falls on: `+1` right,
/// `-1` left, `0` exactly on it.
///
/// **Right and left are named for a sweep edge**, which always points downward, so
/// "right" is the larger-`u` side. That is the orientation the status list is sorted
/// in and the one every comparison below reads.
fn side(a: P2, b: P2, c: P2) -> i8 {
    match orient2d(a, b, c) {
        d if d > 0.0 => 1,
        d if d < 0.0 => -1,
        _ => 0,
    }
}

/// The five cases the sweep distinguishes, by where a vertex's neighbours lie and
/// whether the interior angle is convex.
///
/// Only **split** and **merge** need a diagonal: they are the vertices where the
/// region locally stops being monotone, one opening a new interior span downward and
/// the other closing two of them. The other three are bookkeeping.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Start,
    End,
    Split,
    Merge,
    Regular,
}

/// The decomposition: `rings[0]` is the outer boundary (CCW) and the rest are holes
/// (CW), so **every directed ring edge has the interior on its left** — one uniform
/// rule the whole sweep leans on.
///
/// Returns the monotone pieces, each a CCW ring of the same indices. No vertex is
/// created, and no vertex is repeated within a piece.
///
/// Indices must be `< uv.len()` and each must appear in **exactly one** ring exactly
/// once. That is the b-rep's own invariant (a face's loops are disjoint and its holes
/// are siblings, never nested), but this file does not get to assume it — a violation
/// comes back as [`TessError::DegenerateRing`] rather than as a quietly wrong mesh.
pub(super) fn decompose(
    uv: &[P2],
    rings: &[&[usize]],
    twins: Option<&Twins>,
) -> Result<Vec<Vec<usize>>, TessError> {
    let n = uv.len();
    let (prev, next) = link(n, rings)?;
    // "Wrong" outranks "cannot": a crossing means every triangulation of these rings is wrong,
    // a touch only that this decomposition has no answer for them. Neither says a word about
    // the solid — that is `validate`'s business, one crate over.
    let meets = self_touch(uv, &prev, &next);
    if meets.crossing.is_some() {
        return Err(TessError::DegenerateRing);
    }
    // The bridge's twins touch each other's outgoing segment at their shared point — exactly
    // two `AtEnd` records, `(o on h→…)` and `(h on o→…)` — and those two are the boundary this
    // decomposition *was* given an order for. Any other touch is what it always was.
    let sanctioned = |tc: &Touch| {
        twins.is_some_and(|tw| {
            tc.kind == TouchKind::AtEnd
                && ((tc.vertex == tw.o && tc.segment == tw.h)
                    || (tc.vertex == tw.h && tc.segment == tw.o))
        })
    };
    if meets.touches.iter().any(|tc| !sanctioned(tc)) {
        return Err(TessError::SelfTouchingBoundary);
    }

    let kinds = classify(uv, &prev, &next)?;
    let diagonals = sweep(uv, &prev, &next, &kinds, twins)?;
    trace(uv, &next, &diagonals)
}

/// Ring order as `prev`/`next` arrays, and the precondition check that makes them
/// well defined: **every vertex is used exactly once.** A repeated index would give
/// one vertex two successors, and the sweep would have no way to know which chain it
/// was on.
pub(super) fn link(n: usize, rings: &[&[usize]]) -> Result<(Vec<usize>, Vec<usize>), TessError> {
    const NONE: usize = usize::MAX;
    let (mut prev, mut next) = (vec![NONE; n], vec![NONE; n]);
    for r in rings {
        for k in 0..r.len() {
            let (a, b) = (r[k], r[(k + 1) % r.len()]);
            if a >= n || b >= n || next[a] != NONE || prev[b] != NONE {
                return Err(TessError::DegenerateRing);
            }
            next[a] = b;
            prev[b] = a;
        }
    }
    if next.contains(&NONE) {
        return Err(TessError::DegenerateRing);
    }
    Ok((prev, next))
}

/// **Every way the boundary meets itself, and what each meeting is worth.**
///
/// The geometric twin of [`link`]'s check one function up. That one refuses two rings that share
/// an **index** — *"a pinch, and the sweep would have no way to know which chain it was on"* — and
/// this one finds the same pinch spelled in **coordinates**, which is the spelling a chart
/// actually produces: `face_rings` gives every ring its own index range, so two rings that touch
/// do it with distinct indices at one point.
///
/// Two scans, and **both always run** — the second used to be skipped whenever the first found
/// anything, so one touch hid every crossing on the face:
///
/// 1. **A vertex on a segment it does not belong to** (`strictly_between`, or coincident with
///    one of the segment's ends). The two segments incident to the vertex are skipped — every
///    vertex touches those. Exact throughout: `orient2d` decides collinearity and the rest is
///    `f64` comparison on stored coordinates, like everything else in this file.
/// 2. **Two segments passing through each other**, ends strictly straddling each other's line.
///
/// ★★★★★ **A touch vouches for nothing by itself — so each one is asked to.** A vertex sitting on
/// another segment is either a *tangency* (the boundary reaches the segment and turns back) or a
/// *crossing that happens to land on a vertex* (the boundary passes through). One `orient2d`
/// pair tells them apart: the vertex's two ring neighbours lie on the **same** side of the
/// segment's line for a tangency and on **opposite** sides for a crossing. The test is complete
/// only where the other ring is straight at the touch — a vertex strictly *inside* a segment
/// ([`TouchKind::Interior`]). Where two rings share a vertex coordinate ([`TouchKind::AtEnd`])
/// the other ring bends there too and one line is not enough; that case, and a neighbour lying
/// exactly on the line, are left [`Witness::Abstained`] — recorded, never promoted. ☑ Measured:
/// a diamond hole crossing the outer ring at its own two vertices came back as a *touch* before
/// this witness existed, and a bridge laid on that evidence would have meshed a region outside
/// the outer ring.
///
/// ★ **Why a vertex is enough for the touching case.** A tangency is one point, and the sampled
/// boundary reproduces it only when a sample happens to land there — which both measured fixtures
/// do, their seams sitting on the touch. A tangency whose sampling *misses* the point produces a
/// boundary that is genuinely separated, meshes, and is wrong only by the sampling error it was
/// always allowed. So this looks for the vertex, not for the tangency.
///
/// ★★★★★ **A *crossing* is the other answer, and it is a different sentence.** Two segments that
/// pass through each other with no vertex at the crossing make the ring set not a polygon at all —
/// [`TessError::DegenerateRing`]: any triangulation of it would be wrong — where a touch is a
/// boundary this decomposition cannot draw. Neither name is a verdict on the solid.
///
/// ☑ **The crossing scan is strict, and must be.** With the scan no longer short-circuited by a
/// touch, an endpoint sitting exactly on the other segment's line (`orient2d == 0`) reaches it —
/// and the old form, which folded a zero in with the negative side, called the measured tangency
/// fixture a crossing. A zero at an endpoint means the segments meet the line only there, which
/// the touch scan already names; a straddle needs four non-zero signs.
///
/// ☑ **That branch is here because the gap was measured, not imagined.** A self-crossing single
/// ring (a bow-tie) already came back `DegenerateRing` from the sweep, but a **hole crossing its
/// outer ring** came back `Ok` — eight confident, wrong triangles. This file's charter is that a
/// wrong cache is worse than none, so the check that was going to *document* that gap closes it
/// instead.
pub(super) fn self_touch(uv: &[P2], prev: &[usize], next: &[usize]) -> Meets {
    let box_misses = |p: P2, u: P2, v: P2| {
        p[0] < u[0].min(v[0])
            || p[0] > u[0].max(v[0])
            || p[1] < u[1].min(v[1])
            || p[1] > u[1].max(v[1])
    };
    let mut meets = Meets {
        touches: Vec::new(),
        crossing: None,
    };
    // 1. A vertex sitting on a segment it does not belong to — the touch, asked for its witness.
    for (i, &p) in uv.iter().enumerate() {
        for a in 0..uv.len() {
            let b = next[a];
            if a == i || b == i {
                continue;
            }
            let (u, v) = (uv[a], uv[b]);
            // The segment's box first: four comparisons instead of an exact predicate.
            if box_misses(p, u, v) {
                continue;
            }
            if orient2d(u, v, p) != 0.0 {
                continue;
            }
            let kind = if strictly_between(u, v, p) {
                TouchKind::Interior
            } else if p == u {
                TouchKind::AtEnd
            } else if p == v {
                // The same shared coordinate is the start of `next[a]`, recorded there — one
                // record per touch, or a bridge laid on it would be laid twice.
                continue;
            } else {
                continue;
            };
            let witness = match kind {
                TouchKind::AtEnd => Witness::Abstained,
                TouchKind::Interior => {
                    let (s_prev, s_next) = (side(u, v, uv[prev[i]]), side(u, v, uv[next[i]]));
                    if s_prev == 0 || s_next == 0 {
                        Witness::Abstained
                    } else if s_prev == s_next {
                        Witness::Tangent
                    } else {
                        // The boundary passes through the segment here: a crossing that landed
                        // on a vertex. It is not a touch at all.
                        if meets.crossing.is_none() {
                            meets.crossing = Some((i, a));
                        }
                        continue;
                    }
                }
            };
            meets.touches.push(Touch {
                vertex: i,
                segment: a,
                segment_end: b,
                kind,
                witness,
            });
        }
    }
    // 2. Two segments passing through each other — always looked for, never hidden by a touch.
    for a in 0..uv.len() {
        let (b, c_lo) = (next[a], uv[a]);
        let c_hi = uv[b];
        for c in (a + 1)..uv.len() {
            let d = next[c];
            if c == a || c == b || d == a || d == b {
                continue;
            }
            let (e, f) = (uv[c], uv[d]);
            if c_lo[0].max(c_hi[0]) < e[0].min(f[0])
                || e[0].max(f[0]) < c_lo[0].min(c_hi[0])
                || c_lo[1].max(c_hi[1]) < e[1].min(f[1])
                || e[1].max(f[1]) < c_lo[1].min(c_hi[1])
            {
                continue;
            }
            // Strict: four non-zero signs. An endpoint on the other line is a touch, named above.
            let straddles = |p: P2, q: P2, r: P2, s: P2| {
                let (x, y) = (orient2d(p, q, r), orient2d(p, q, s));
                x != 0.0 && y != 0.0 && ((x > 0.0) != (y > 0.0))
            };
            if straddles(c_lo, c_hi, e, f) && straddles(e, f, c_lo, c_hi) {
                if meets.crossing.is_none() {
                    meets.crossing = Some((a, c));
                }
                return meets;
            }
        }
    }
    meets
}

/// Where on its segment a touching vertex sits — the bit a bridge builder branches on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum TouchKind {
    /// The vertex coincides with one of the segment's ends: two rings share a coordinate, and
    /// both bend there. A bridge between them costs nothing.
    AtEnd,
    /// The vertex sits strictly inside the segment: the other ring is straight there. A bridge
    /// needs that segment split at the vertex first — and the split must reach every face that
    /// shares the edge, or the mesh cracks.
    Interior,
}

/// What the tangency witness said about a touch.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Witness {
    /// Both ring neighbours on one side of the touched line: the boundary turns back here.
    Tangent,
    /// Not decidable from one line — an `AtEnd` touch, or a neighbour exactly on the line.
    /// Kept as a touch (this decomposition still has no answer), never promoted to a crossing.
    Abstained,
}

/// One vertex on one segment it does not belong to — see [`self_touch`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Touch {
    pub(crate) vertex: usize,
    /// The segment `segment → segment_end` (`segment_end == next[segment]`, carried so a caller
    /// without the ring arrays can find the segment in its own polylines).
    pub(crate) segment: usize,
    pub(crate) segment_end: usize,
    pub(crate) kind: TouchKind,
    pub(crate) witness: Witness,
}

/// Everything [`self_touch`] found. `touches` holds one record per touch — a vertex on a shared
/// endpoint is recorded for the segment that *starts* there only. A crossing, whether found by
/// the segment scan or by a touch whose witness said the boundary passes through, outranks every
/// touch.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct Meets {
    pub(crate) touches: Vec<Touch>,
    pub(crate) crossing: Option<(usize, usize)>,
}

fn classify(uv: &[P2], prev: &[usize], next: &[usize]) -> Result<Vec<Kind>, TessError> {
    (0..uv.len())
        .map(|i| {
            let (p, c, q) = (uv[prev[i]], uv[i], uv[next[i]]);
            let below = |x: P2| lex_less(c, x); // `c` comes first ⇒ `x` is below it
            let turn = side(p, c, q);
            Ok(match (below(p) && below(q), below(p) || below(q), turn) {
                // A vertex whose neighbours are both on one side and yet collinear
                // with it is a spike: the boundary doubles back, so the ring is not
                // simple and no triangulation of it exists.
                (true, _, 0) => return Err(TessError::DegenerateRing),
                (_, false, 0) => return Err(TessError::DegenerateRing),
                // `turn < 0` is reflex, because a ring edge keeps the interior on its
                // left: turning right means the interior angle opened past π.
                (true, _, 1) => Kind::Start,
                (true, _, _) => Kind::Split,
                (false, false, 1) => Kind::End,
                (false, false, _) => Kind::Merge,
                _ => Kind::Regular,
            })
        })
        .collect()
}

/// One entry of the sweep status: an edge crossing the sweep line, and the lowest
/// vertex seen so far that the region to that edge's *left* is responsible for.
///
/// The helper is the whole trick. A merge vertex leaves the region above it split in
/// two with no diagonal yet drawn; recording it as a helper means the next vertex that
/// can see it draws that diagonal, and the "can see it" question is answered by the
/// status order rather than by any visibility test.
struct Edge {
    /// The edge runs `upper → next[upper]`, downward in sweep order.
    upper: usize,
    helper: usize,
}

fn sweep(
    uv: &[P2],
    prev: &[usize],
    next: &[usize],
    kinds: &[Kind],
    twins: Option<&Twins>,
) -> Result<Vec<[usize; 2]>, TessError> {
    let n = uv.len();
    // The one tie the order allows is the sanctioned pair, decided once here (a comparator
    // cannot fail) and read by the sort below.
    let twin_first = match twins {
        Some(tw) => Some(lex_less_idx(tw.o, tw.h, uv, twins).ok_or(TessError::DegenerateRing)?),
        None => None,
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        if let (Some(tw), Some(first)) = (twins, twin_first) {
            if tw.is_pair(a, b) {
                return if (a == tw.o) == first {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                };
            }
        }
        if lex_less(uv[a], uv[b]) {
            std::cmp::Ordering::Less
        } else if lex_less(uv[b], uv[a]) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    // Two vertices at the same point make `lex_less` a non-strict order, and every
    // "which is above" question below would then have no answer — except the sanctioned
    // twins, whose order `sos` supplies.
    if order
        .windows(2)
        .any(|w| uv[w[0]] == uv[w[1]] && !twins.is_some_and(|tw| tw.is_pair(w[0], w[1])))
    {
        return Err(TessError::DegenerateRing);
    }

    // Sorted left to right by where each edge crosses the sweep line. Edges in the
    // status never cross, so an order established at insertion stays correct.
    let mut status: Vec<Edge> = Vec::new();
    let mut out: Vec<[usize; 2]> = Vec::new();

    // The edge of `status` immediately left of `v` — the one whose region `v` falls
    // into. `None` means `v` is outside every span, which a valid polygon cannot be.
    // A zero here is structural only with the twins (a twin queried against the other twin's
    // edge); `side_idx` decides it, and a zero it cannot decide is a refusal, not a guess.
    let left_of = |status: &Vec<Edge>, v: usize, uv: &[P2]| -> Result<Option<usize>, TessError> {
        for (i, e) in status.iter().enumerate().rev() {
            let s =
                side_idx(e.upper, next[e.upper], v, uv, twins).ok_or(TessError::DegenerateRing)?;
            if s > 0 {
                return Ok(Some(i));
            }
        }
        Ok(None)
    };
    let insert =
        |status: &mut Vec<Edge>, upper: usize, helper: usize, uv: &[P2]| -> Result<(), TessError> {
            let mut at = status.len();
            for (i, e) in status.iter().enumerate() {
                let s = side_idx(e.upper, next[e.upper], upper, uv, twins)
                    .ok_or(TessError::DegenerateRing)?;
                if s < 0 {
                    at = i;
                    break;
                }
            }
            status.insert(at, Edge { upper, helper });
            Ok(())
        };
    let remove = |status: &mut Vec<Edge>, upper: usize| -> Result<usize, TessError> {
        let at = status
            .iter()
            .position(|e| e.upper == upper)
            .ok_or(TessError::DegenerateRing)?;
        Ok(status.remove(at).helper)
    };

    for &v in &order {
        match kinds[v] {
            Kind::Start => insert(&mut status, v, v, uv)?,
            Kind::End => {
                let h = remove(&mut status, prev[v])?;
                if kinds[h] == Kind::Merge {
                    out.push([v, h]);
                }
            }
            Kind::Split => {
                let at = left_of(&status, v, uv)?.ok_or(TessError::DegenerateRing)?;
                out.push([v, status[at].helper]);
                status[at].helper = v;
                insert(&mut status, v, v, uv)?;
            }
            Kind::Merge => {
                let h = remove(&mut status, prev[v])?;
                if kinds[h] == Kind::Merge {
                    out.push([v, h]);
                }
                let at = left_of(&status, v, uv)?.ok_or(TessError::DegenerateRing)?;
                if kinds[status[at].helper] == Kind::Merge {
                    out.push([v, status[at].helper]);
                }
                status[at].helper = v;
            }
            // The two chains are not symmetric: on the one where the boundary runs
            // *downward* through `v`, the interior lies to `v`'s right and `v` owns
            // the edge above it; on the other it owns nothing and only updates the
            // helper of the edge to its left.
            Kind::Regular => {
                if lex_less(uv[prev[v]], uv[v]) {
                    let h = remove(&mut status, prev[v])?;
                    if kinds[h] == Kind::Merge {
                        out.push([v, h]);
                    }
                    insert(&mut status, v, v, uv)?;
                } else {
                    let at = left_of(&status, v, uv)?.ok_or(TessError::DegenerateRing)?;
                    if kinds[status[at].helper] == Kind::Merge {
                        out.push([v, status[at].helper]);
                    }
                    status[at].helper = v;
                }
            }
        }
    }
    if !status.is_empty() {
        return Err(TessError::DegenerateRing);
    }
    Ok(out)
}

/// Walk the subdivision made by the ring edges plus the diagonals, and hand back the
/// face on the left of every directed edge.
///
/// At each step the successor of `a → b` is the outgoing edge of `b` that comes
/// **first clockwise** from `b → a` — the standard rule, and the reason it needs no
/// geometry beyond an angular order: turning as sharply as possible keeps the walk
/// hugging the same face.
fn trace(
    uv: &[P2],
    next: &[usize],
    diagonals: &[[usize; 2]],
) -> Result<Vec<Vec<usize>>, TessError> {
    let n = uv.len();
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); n];
    for a in 0..n {
        outgoing[a].push(next[a]);
    }
    for d in diagonals {
        outgoing[d[0]].push(d[1]);
        outgoing[d[1]].push(d[0]);
    }
    // Counter-clockwise by direction, starting from +u. Exact: the half decides on a
    // coordinate comparison and the tie on one `orient2d`.
    let half = |b: usize, c: usize| -> u8 {
        let (p, q) = (uv[b], uv[c]);
        u8::from(!(q[1] > p[1] || (q[1] == p[1] && q[0] > p[0])))
    };
    for (b, outs) in outgoing.iter_mut().enumerate() {
        outs.sort_by(|&c1, &c2| {
            half(b, c1)
                .cmp(&half(b, c2))
                .then_with(|| side(uv[b], uv[c1], uv[c2]).cmp(&0).reverse())
        });
    }

    let mut used: Vec<Vec<bool>> = outgoing.iter().map(|o| vec![false; o.len()]).collect();
    let mut pieces = Vec::new();
    for a0 in 0..n {
        for k0 in 0..outgoing[a0].len() {
            if used[a0][k0] {
                continue;
            }
            let mut piece = Vec::new();
            let (mut a, mut k) = (a0, k0);
            loop {
                used[a][k] = true;
                piece.push(a);
                let b = outgoing[a][k];
                // The predecessor of `b → a` in the CCW order around `b` is the first
                // edge clockwise from it.
                let outs = &outgoing[b];
                let at = outs
                    .iter()
                    .position(|&c| c == a)
                    .map(|i| (i + outs.len() - 1) % outs.len())
                    .unwrap_or_else(|| {
                        // `a` is not a neighbour of `b` in this direction (a one-way
                        // ring edge), so locate where it would sit and step back.
                        let cmp = |&c: &usize| {
                            half(b, c)
                                .cmp(&half(b, a))
                                .then_with(|| side(uv[b], uv[c], uv[a]).cmp(&0).reverse())
                        };
                        let i = outs.partition_point(|c| cmp(c) == std::cmp::Ordering::Less);
                        (i + outs.len() - 1) % outs.len()
                    });
                if piece.len() > 4 * n + 4 {
                    return Err(TessError::DegenerateRing);
                }
                (a, k) = (b, at);
                if (a, k) == (a0, k0) {
                    break;
                }
            }
            pieces.push(piece);
        }
    }
    Ok(pieces)
}

/// Triangulate one **y-monotone** piece by the stack algorithm (de Berg §3.3).
///
/// Linear, and it owes that to monotonicity: the two chains can be merged into one
/// sweep-ordered sequence, and a stack of the vertices not yet triangulated is enough
/// — no search for a cuttable corner and no containment test over the rest of the
/// polygon, which is the work ear clipping repeats for every triangle it emits.
///
/// **Winding is normalized at emission rather than reasoned about per case.** Which of
/// the four cases produced a triangle does not change what it *is*, only the order the
/// three indices arrive in, and one `orient2d` settles that. A triangle that comes out
/// with no area at all is not a presentation problem — the algorithm was not supposed
/// to make one — so it is an error.
pub(super) fn triangulate_monotone(
    uv: &[P2],
    piece: &[usize],
    out: &mut Vec<[usize; 3]>,
    twins: Option<&Twins>,
) -> Result<(), TessError> {
    let m = piece.len();
    if m < 3 {
        return Err(TessError::DegenerateRing);
    }
    // The sweep's order, twins included — the same order `sweep` used, or the chains below
    // would not be the chains it cut.
    let below = |a: usize, b: usize| lex_less_idx(a, b, uv, twins).ok_or(TessError::DegenerateRing);
    let (mut top, mut bot) = (0, 0);
    for i in 1..m {
        if below(piece[i], piece[top])? {
            top = i;
        }
        if below(piece[bot], piece[i])? {
            bot = i;
        }
    }

    // Walking the CCW ring forward from the top descends the **left** chain, so the
    // interior lies to its right; walking backward descends the right chain.
    let walk = |mut i: usize, step: usize| -> Vec<usize> {
        let mut v = Vec::new();
        while i != bot {
            i = (i + step) % m;
            v.push(i);
        }
        v.pop(); // `bot` closes both chains and is merged in last
        v
    };
    let (left, right) = (walk(top, 1), walk(top, m - 1));
    // The precondition this function never checked in production: each chain descends. A
    // piece that does not is not the sweep's — and with twins in play, the one place a wrong
    // order would otherwise pass silently.
    debug_assert!(
        [&left, &right].iter().all(|chain| {
            std::iter::once(top)
                .chain(chain.iter().copied())
                .chain(std::iter::once(bot))
                .collect::<Vec<_>>()
                .windows(2)
                .all(|w| below(piece[w[0]], piece[w[1]]).unwrap_or(false))
        }),
        "a monotone piece's chains must descend in sweep order"
    );

    // One sweep-ordered sequence, each vertex tagged with the chain it came from.
    let mut seq: Vec<(usize, bool)> = Vec::with_capacity(m);
    seq.push((top, true));
    let (mut li, mut ri) = (0, 0);
    while li < left.len() || ri < right.len() {
        let take_left = if li >= left.len() {
            false
        } else if ri >= right.len() {
            true
        } else {
            below(piece[left[li]], piece[right[ri]])?
        };
        if take_left {
            seq.push((left[li], true));
            li += 1;
        } else {
            seq.push((right[ri], false));
            ri += 1;
        }
    }
    seq.push((bot, false));

    let mut emit = |a: usize, b: usize, c: usize| -> Result<(), TessError> {
        let (x, y, z) = (piece[a], piece[b], piece[c]);
        match side(uv[x], uv[y], uv[z]) {
            1 => out.push([x, y, z]),
            -1 => out.push([x, z, y]),
            _ => return Err(TessError::DegenerateRing),
        }
        Ok(())
    };

    let mut stack: Vec<(usize, bool)> = vec![seq[0], seq[1]];
    for j in 2..seq.len() - 1 {
        let (vj, cj) = seq[j];
        if cj != stack[stack.len() - 1].1 {
            // The opposite chain is reachable from every vertex still on the stack, so
            // the whole fan comes off at once.
            while stack.len() > 1 {
                let (a, _) = stack.pop().expect("len > 1");
                let (b, _) = stack[stack.len() - 1];
                emit(a, b, vj)?;
            }
            stack.clear();
            stack.push(seq[j - 1]);
            stack.push(seq[j]);
        } else {
            // Same chain: cut while the boundary bends *away* from the interior. On the
            // left chain the interior is to the right, so a left turn at the popped
            // vertex means the triangle clears the boundary; on the right chain it is
            // the mirror of that.
            let want = if cj { 1 } else { -1 };
            let mut last = stack.pop().expect("two seeded");
            while let Some(&(t, _)) = stack.last() {
                if side(uv[piece[t]], uv[piece[last.0]], uv[piece[vj]]) != want {
                    break;
                }
                emit(t, last.0, vj)?;
                last = stack.pop().expect("checked");
            }
            stack.push(last);
            stack.push(seq[j]);
        }
    }
    let (vn, _) = seq[seq.len() - 1];
    while stack.len() > 1 {
        let (a, _) = stack.pop().expect("len > 1");
        let (b, _) = stack[stack.len() - 1];
        emit(a, b, vn)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uv(v: &[[f64; 2]]) -> Vec<P2> {
        v.to_vec()
    }

    /// Every piece is `v`-monotone: walking it, the sweep order rises to exactly one
    /// peak and falls to exactly one trough. That is the whole proposition of this
    /// module — a piece that fails it cannot be triangulated by the stack algorithm,
    /// and a decomposition that produces one has not decomposed anything.
    fn assert_monotone(uv: &[P2], piece: &[usize]) {
        let n = piece.len();
        assert!(n >= 3, "piece {piece:?} is not a polygon");
        let turns = (0..n)
            .filter(|&i| {
                let (p, c, q) = (
                    uv[piece[(i + n - 1) % n]],
                    uv[piece[i]],
                    uv[piece[(i + 1) % n]],
                );
                lex_less(c, p) == lex_less(c, q)
            })
            .count();
        assert_eq!(turns, 2, "piece {piece:?} has {turns} extrema, want 2");
    }

    fn check(pts: &[[f64; 2]], rings: &[&[usize]]) -> Vec<Vec<usize>> {
        let uv = uv(pts);
        let pieces = decompose(&uv, rings, None).expect("decomposes");
        for p in &pieces {
            assert_monotone(&uv, p);
        }
        // The pieces partition the polygon, so their areas sum to its area.
        let area = |r: &[usize]| -> f64 {
            let mut s = 0.0;
            for k in 0..r.len() {
                let (a, b) = (uv[r[k]], uv[r[(k + 1) % r.len()]]);
                s += a[0] * b[1] - b[0] * a[1];
            }
            0.5 * s
        };
        let want: f64 = rings.iter().map(|r| area(r)).sum();
        let got: f64 = pieces.iter().map(|p| area(p)).sum();
        assert!(
            (got - want).abs() < 1e-9,
            "piece area {got} vs polygon {want} ({pieces:?})"
        );
        pieces
    }

    /// **The case this module exists for, and the one CAD actually produces**: every
    /// edge axis-aligned, so the sweep coordinate is shared by pairs of vertices all
    /// over. If the total order were not total, this is where it would show.
    #[test]
    fn a_rectangle_with_a_rectangular_hole() {
        let pieces = check(
            &[
                [0.0, 0.0],
                [4.0, 0.0],
                [4.0, 3.0],
                [0.0, 3.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [3.0, 2.0],
                [3.0, 1.0],
            ],
            &[&[0, 1, 2, 3], &[4, 5, 6, 7]],
        );
        assert_eq!(pieces.len(), 2, "one hole splits it in two: {pieces:?}");
    }

    /// A convex ring is already monotone, so the sweep must add nothing.
    #[test]
    fn a_square_is_left_alone() {
        let pieces = check(
            &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            &[&[0, 1, 2, 3]],
        );
        assert_eq!(pieces.len(), 1);
    }

    /// One split and one merge vertex, and no holes — the notch supplies both.
    #[test]
    fn a_u_shape_needs_its_diagonals() {
        check(
            &[
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 3.0],
                [2.0, 3.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 3.0],
                [0.0, 3.0],
            ],
            &[&(0..8).collect::<Vec<_>>()],
        );
    }

    /// **The measured failure**, transcribed: the two holes whose bridges collided on
    /// one outer corner and left ear clipping with no ear at all.
    #[test]
    fn the_wall_with_two_windows() {
        let pieces = check(
            &[
                [0.0, 1.0],
                [0.0, -1.0],
                [3.0, -1.0],
                [3.0, 1.0],
                [1.0, 0.25959136597258015],
                [1.0, 0.7035578716424774],
                [2.0, 0.7035578716424774],
                [2.0, 0.25959136597258015],
                [1.0, -0.703557871642477],
                [1.0, -0.25959136597258003],
                [2.0, -0.25959136597258003],
                [2.0, -0.703557871642477],
            ],
            &[&[0, 1, 2, 3], &[7, 6, 5, 4], &[11, 10, 9, 8]],
        );
        assert!(pieces.len() >= 3, "two holes need at least two cuts");
    }

    /// **The corpus that measured the old triangulator's 27.6% failure rate**, pointed
    /// at the decomposition instead: random axis-aligned rectangles with one to four
    /// non-overlapping rectangular holes on an integer grid.
    ///
    /// Hand-picked fixtures only cover the degeneracies someone thought of. This one
    /// produces shared sweep coordinates, horizontal edges and collinear vertices by
    /// the thousand, which is precisely where a sweep goes wrong — and every sample
    /// must come back as monotone pieces whose areas add up.
    #[test]
    fn random_rectilinear_faces_decompose() {
        let mut st = 0x5EED_1234_ABCD_0001u64;
        let mut lcg = move || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            st >> 33
        };
        let mut rng = move |lo: i64, hi: i64| lo + (lcg() % ((hi - lo + 1) as u64)) as i64;
        let (w, h) = (24i64, 16i64);
        let mut samples = 0;
        for _ in 0..3000 {
            let k = rng(1, 4);
            let mut boxes: Vec<[i64; 4]> = Vec::new();
            for _ in 0..k {
                for _try in 0..20 {
                    let x0 = rng(1, w - 3);
                    let y0 = rng(1, h - 3);
                    let x1 = x0 + rng(1, (w - 1 - x0).min(5));
                    let y1 = y0 + rng(1, (h - 1 - y0).min(5));
                    if boxes
                        .iter()
                        .any(|b| x0 <= b[2] && b[0] <= x1 && y0 <= b[3] && b[1] <= y1)
                    {
                        continue;
                    }
                    boxes.push([x0, y0, x1, y1]);
                    break;
                }
            }
            if boxes.is_empty() {
                continue;
            }
            let mut pts = vec![
                [0.0, 0.0],
                [w as f64, 0.0],
                [w as f64, h as f64],
                [0.0, h as f64],
            ];
            let mut holes: Vec<Vec<usize>> = Vec::new();
            for b in &boxes {
                let base = pts.len();
                // Clockwise, so the interior stays on every edge's left.
                for c in [[b[0], b[1]], [b[0], b[3]], [b[2], b[3]], [b[2], b[1]]] {
                    pts.push([c[0] as f64, c[1] as f64]);
                }
                holes.push((base..base + 4).collect());
            }
            let mut rings: Vec<&[usize]> = vec![&[0, 1, 2, 3]];
            rings.extend(holes.iter().map(|h| h.as_slice()));
            check(&pts, &rings);
            samples += 1;
        }
        assert!(samples > 2500, "only {samples} samples reached the gate");
    }

    /// A vertex used by two rings is not a polygon with sibling holes — it is a pinch,
    /// and the sweep would have no way to know which chain it was on.
    #[test]
    fn a_shared_vertex_is_rejected() {
        let pts = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 3.0],
            [0.0, 3.0],
            [1.0, 1.0],
            [1.0, 2.0],
        ];
        assert!(matches!(
            decompose(&uv(&pts), &[&[0, 1, 2, 3], &[4, 5, 0]], None),
            Err(TessError::DegenerateRing)
        ));
    }
}

#[cfg(test)]
mod touch_tests {
    //! The tangency witness itself, read off `self_touch`'s output — `polygon`'s tests only see
    //! the error name, and the name cannot tell `Tangent` from `Abstained`.
    use super::*;

    fn meets(uv: &[P2], rings: &[&[usize]]) -> Meets {
        let (prev, next) = link(uv.len(), rings).expect("rings are well formed");
        self_touch(uv, &prev, &next)
    }

    /// A hole's apex on the outer ring's edge: one record, strictly inside the segment, and both
    /// neighbours on the inner side — the witness says tangent.
    #[test]
    fn an_apex_on_an_edge_is_an_interior_tangency() {
        let uv = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [2.0, 0.0],
            [1.0, 2.0],
            [3.0, 2.0],
        ];
        let m = meets(&uv, &[&[0, 1, 2, 3], &[4, 5, 6]]);
        assert_eq!(m.crossing, None);
        assert_eq!(
            m.touches,
            vec![Touch {
                vertex: 4,
                segment: 0,
                segment_end: 1,
                kind: TouchKind::Interior,
                witness: Witness::Tangent,
            }]
        );
    }

    /// A diamond crossing the outer ring at its own two vertices: the witness sees the neighbours
    /// straddle the edge, so neither vertex is a touch — it is a crossing.
    #[test]
    fn a_vertex_the_boundary_passes_through_is_a_crossing_not_a_touch() {
        let uv = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [3.0, 2.0],
            [4.0, 3.0],
            [5.0, 2.0],
            [4.0, 1.0],
        ];
        let m = meets(&uv, &[&[0, 1, 2, 3], &[4, 5, 6, 7]]);
        assert!(m.crossing.is_some(), "{m:?}");
        assert!(m.touches.is_empty(), "{m:?}");
    }

    /// Two rings sharing a corner coordinate: recorded once per (vertex, other ring), never for
    /// both segments meeting at the shared point; and the witness abstains, since the other ring
    /// bends there too.
    #[test]
    fn a_shared_corner_is_recorded_once_and_abstains() {
        let uv = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [4.0, 0.0],
            [2.0, 1.0],
            [3.0, 2.0],
        ];
        let m = meets(&uv, &[&[0, 1, 2, 3], &[4, 5, 6]]);
        assert_eq!(m.crossing, None, "{m:?}");
        assert_eq!(m.touches.len(), 2, "{m:?}");
        for t in &m.touches {
            assert_eq!(t.kind, TouchKind::AtEnd, "{t:?}");
            assert_eq!(t.witness, Witness::Abstained, "{t:?}");
        }
        // The hole's vertex on the outer's segment starting at the corner, and the outer's corner
        // on the hole's segment starting there.
        assert!(
            m.touches.iter().any(|t| t.vertex == 4 && t.segment == 1),
            "{m:?}"
        );
        assert!(
            m.touches.iter().any(|t| t.vertex == 1 && t.segment == 4),
            "{m:?}"
        );
    }

    /// A hole edge lying along the outer edge: each of its ends is on the segment, but the other
    /// neighbour is on the line too — the witness cannot decide and says so.
    #[test]
    fn a_neighbour_on_the_line_makes_the_witness_abstain() {
        let uv = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 0.0],
            [2.0, 2.0],
            [3.0, 0.0],
        ];
        let m = meets(&uv, &[&[0, 1, 2, 3], &[4, 5, 6]]);
        assert_eq!(m.crossing, None, "{m:?}");
        assert_eq!(m.touches.len(), 2, "{m:?}");
        for t in &m.touches {
            assert_eq!(t.kind, TouchKind::Interior, "{t:?}");
            assert_eq!(t.witness, Witness::Abstained, "{t:?}");
        }
    }
}
