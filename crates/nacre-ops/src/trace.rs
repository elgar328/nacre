#![cfg_attr(not(test), allow(dead_code))]
//! Segment trace producer for the plane-class arrangement (the B engine). Isolated and unwired:
//! nothing outside this module and its tests calls it, following the `winding.rs` precedent.
//!
//! Unlike [`crate::section_of_solid`], which rejects the whole solid the moment any vertex lies
//! on the cut plane and then assembles **closed loops** (needing per-face even parity, degree-2,
//! and an alternating walk), this produces **line segments** — the pieces where each face of the
//! solid meets a plane class — and leaves assembly to the arrangement. Geometry lying in the
//! plane is *input*, not a degeneracy.
//!
//! Two face cases are handled. **Seated** — a face whose own class is the cut class — lies wholly
//! in the plane, so its whole boundary is trace. **Transversal / mixed** — a face crossing the
//! plane, possibly with an edge lying in it — is 2D polygon-vs-line clipping on `L = W ∩ fp`, with
//! a three-valued (−/0/+) scan that collapses on-line edges instead of rejecting them. Faces with
//! **holes** are declined this brick (a forced-covered run inside a hole would claim its void as
//! material, and the corpus has no case to verify it). Declining is a per-face record, never a
//! whole-solid abort — a non-empty `declined` means "this trace is incomplete, do not conclude
//! from it", which is the defect `section_of_solid` had.

use super::*;

/// Which operand a segment came from — the boolean's per-cell label needs both solids' material
/// above and below, so provenance cannot be merged away.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SolidSide {
    A,
    B,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SegKind {
    /// A face crossing the plane transversally. Crossing this segment in the arrangement flips
    /// **both** the above-label and the below-label together. `mat` is which side of the directed
    /// line `W ∩ fp` the solid's material lies on (`+1`/`-1`, relative to `d = n_W × n_fp`).
    Transversal { mat: i8 },
    /// A boundary edge of a face lying in the plane. Crossing it flips **one** label — the side
    /// the seated face's body occupies (`body_above`).
    Seated { body_above: bool },
}

/// One segment of a solid's trace on a plane class, named entirely in plane-class triples.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Seg {
    /// Canon plane class the segment rides besides the cut plane (its line is `W ∩ wall`).
    pub wall: usize,
    /// Endpoint triples, canonized.
    pub end: [[usize; 3]; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// Single-point tangential contacts — real arrangement vertices, but not segments (a
    /// zero-length chord would abort at `Line::through_points`).
    pub touches: Vec<[usize; 3]>,
    /// `(face plane index, reason)` for every face this brick could not trace.
    pub declined: Vec<(usize, &'static str)>,
}

/// Canonize a raw triple and sort it.
fn canon3(t: [usize; 3], canon: &[usize]) -> [usize; 3] {
    let mut c = [canon[t[0]], canon[t[1]], canon[t[2]]];
    c.sort_unstable();
    c
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
}

/// Trace one **transversal or mixed** face `f` (plane `fp ≠ W`) onto plane class `wc`, as the
/// segments where `f` meets `W`. This is 2D polygon-vs-line clipping on `L = W ∩ fp`, handling
/// on-line edges (a run of `side == 0` vertices) instead of rejecting them the way
/// `section_of_solid` does. Holes are declined this brick (a forced-covered run inside a hole
/// would wrongly claim the hole's void as material, and the corpus has no case to verify it).
///
/// Every emitted segment rides `wall = fp` (its line is `W ∩ fp`); the side walls appear only as
/// the third plane naming each endpoint. `mat = orient_sign(fp)` is a per-face constant.
#[allow(clippy::too_many_arguments)]
fn trace_transversal_face(
    model: &Model,
    fh: Handle<Face>,
    fp: usize,
    which: SolidSide,
    wc: usize,
    planes: &[PlaneInfo],
    inc: &arrange::EdgePlanes,
    canon: &[usize],
    out: &mut Trace,
) {
    if arrange::hole_rings(model, fh, fp, inc).is_ok_and(|h| !h.is_empty()) {
        out.declined.push((fp, "has-holes"));
        return;
    }
    let ring = match arrange::face_vertex_triples(model, fh, fp, inc) {
        Ok(r) => r,
        Err(_) => {
            out.declined.push((fp, "outer-ring"));
            return;
        }
    };
    let n = ring.len();
    let side: Vec<i8> = (0..n)
        .map(|i| arrange::side_of(planes, ring[i], wc))
        .collect();
    let Some(start) = side.iter().position(|&s| s != 0) else {
        // A face off `W` cannot have every vertex on `W`; if it does, do not guess.
        out.declined.push((fp, "all-on-plane"));
        return;
    };

    // `L`'s third plane naming a point of the on-line edge: the vertex triple `{fp, W-class, r}`.
    let third_on_l = |t: [usize; 3]| -> Option<usize> {
        let (mut r, mut has_fp, mut has_w) = (None, false, false);
        for &x in &t {
            if x == fp {
                has_fp = true;
            } else if canon[x] == wc {
                has_w = true;
            } else if r.replace(x).is_some() {
                return None; // two off-planes: not a clean point on L
            }
        }
        (has_fp && has_w).then_some(r).flatten()
    };

    // Phase A — one pass around the ring collecting feature nodes.
    let mut nodes: Vec<Node> = Vec::new();
    let mut run_counter = 0usize;
    let mut declined: Option<&'static str> = None;
    let mut j = 0;
    while j < n {
        let i = (start + j) % n;
        if side[i] != 0 {
            let ni = (i + 1) % n;
            if side[ni] != 0 && side[ni] != side[i] {
                // Strict crossing on edge i; its wall is the plane the edge rides besides fp.
                match arrange::ring_edge(fp, &ring, i) {
                    Ok((wall, _, _)) => nodes.push(Node {
                        r: wall,
                        flip: true,
                        run: None,
                        flanks_differ: false,
                        single_touch: false,
                    }),
                    Err(_) => declined = Some("crossing-name"),
                }
            }
            j += 1;
        } else {
            // A maximal run of side==0 vertices (at most two, since three would be a straight
            // angle already rejected by `loop_triples`).
            let run_start = i;
            let mut m = 0;
            while j < n && side[(start + j) % n] == 0 {
                m += 1;
                j += 1;
            }
            let before = side[(run_start + n - 1) % n];
            let after = side[(start + j) % n];
            let flanks_differ = before != after;
            let name = |k: usize| third_on_l(ring[(run_start + k) % n]);
            if m == 1 {
                match name(0) {
                    Some(r) => nodes.push(Node {
                        r,
                        flip: flanks_differ,
                        run: None,
                        flanks_differ,
                        single_touch: !flanks_differ,
                    }),
                    None => declined = Some("run-name"),
                }
            } else if m == 2 {
                let id = run_counter;
                run_counter += 1;
                match (name(0), name(1)) {
                    (Some(ra), Some(rb)) => {
                        for r in [ra, rb] {
                            nodes.push(Node {
                                r,
                                flip: false,
                                run: Some(id),
                                flanks_differ,
                                single_touch: false,
                            });
                        }
                    }
                    _ => declined = Some("run-name"),
                }
            } else {
                declined = Some("long-run");
            }
        }
        if declined.is_some() {
            break;
        }
    }
    if let Some(reason) = declined {
        out.declined.push((fp, reason));
        return;
    }

    // Phase B — order the nodes along L and fix run structure.
    nodes.sort_by(
        |a, b| match arrange::order_along(planes, wc, fp, a.r, b.r) {
            -1 => std::cmp::Ordering::Less,
            1 => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        },
    );
    for w in nodes.windows(2) {
        if arrange::order_along(planes, wc, fp, w[0].r, w[1].r) == 0 {
            out.declined.push((fp, "coincident-features"));
            return;
        }
    }
    // Each run's two nodes must be adjacent after the sort; the hi node carries the flip.
    for k in 0..nodes.len() {
        if let Some(id) = nodes[k].run {
            let prev = k > 0 && nodes[k - 1].run == Some(id);
            let next = k + 1 < nodes.len() && nodes[k + 1].run == Some(id);
            if !prev && !next {
                out.declined.push((fp, "run-split"));
                return;
            }
            if prev {
                nodes[k].flip = nodes[k].flanks_differ;
            }
        }
    }

    // Phase C — sweep left to right; covered gaps merge into segments.
    let mut parity = 0i8;
    let mut seg_start: Option<usize> = None;
    let emit = |a: usize, b: usize, out: &mut Trace| {
        out.segs.push(Seg {
            wall: canon[fp],
            end: [canon3([wc, fp, a], canon), canon3([wc, fp, b], canon)],
            solid: which,
            kind: SegKind::Transversal {
                mat: arrange::orient_sign(planes, fp),
            },
        });
    };
    for k in 0..nodes.len() {
        if nodes[k].single_touch && parity == 0 {
            out.touches.push(canon3([wc, fp, nodes[k].r], canon));
        }
        if nodes[k].flip {
            parity ^= 1;
        }
        let forced =
            k + 1 < nodes.len() && nodes[k].run.is_some() && nodes[k].run == nodes[k + 1].run;
        let covered = parity == 1 || forced;
        match (seg_start, covered) {
            (None, true) => seg_start = Some(nodes[k].r),
            (Some(a), false) => {
                emit(a, nodes[k].r, out);
                seg_start = None;
            }
            _ => {}
        }
    }
    if seg_start.is_some() || parity != 0 {
        // A line enters and leaves a bounded region equally; an unbalanced sweep is degenerate.
        out.declined.push((fp, "odd-parity"));
    }
}

/// Trace one solid on plane class `wc` (a canon root, i.e. an index into `planes`). This brick:
/// seated faces → their boundary as segments; every other face → `declined`.
#[allow(clippy::too_many_arguments)]
fn trace_one(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &arrange::EdgePlanes,
    canon: &[usize],
    out: &mut Trace,
) {
    let w_normal = planes[wc].plane.normal();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let fp = surf_ix[&fh];
            if canon[fp] != wc {
                trace_transversal_face(model, fh, fp, which, wc, planes, inc, canon, out);
                continue;
            }
            // Seated: the face lies in W, so its whole boundary is trace. The body lies on one
            // side of W — `n_out` points away from the body, so the body is above W exactly when
            // `n_out · n_W < 0`. The f64 sign is robust even rotated: seated means `canon[fp]==wc`,
            // so `n_out ∥ n_W` (both unit) and the dot is ≈ ±1, a full unit from the sign boundary
            // (the rotated-tunnel tests exercise this seated path through the cube's own caps).
            let body_above = planes[fp].n_out.dot(w_normal) < 0.0;
            let kind = SegKind::Seated { body_above };
            let mut emit_ring = |tris: &[[usize; 3]]| {
                let n = tris.len();
                for i in 0..n {
                    // Edge i runs vertex i → vertex i+1; the wall it rides is the plane the two
                    // endpoint triples share besides `fp`.
                    let (t0, t1) = (tris[i], tris[(i + 1) % n]);
                    let shared: Vec<usize> = t0
                        .iter()
                        .copied()
                        .filter(|&x| x != fp && t1.contains(&x))
                        .collect();
                    let [wall] = shared[..] else {
                        out.declined.push((fp, "seated-edge-naming"));
                        continue;
                    };
                    out.segs.push(Seg {
                        wall: canon[wall],
                        end: [canon3(t0, canon), canon3(t1, canon)],
                        solid: which,
                        kind,
                    });
                }
            };
            match arrange::face_vertex_triples(model, fh, fp, inc) {
                Ok(ts) => emit_ring(&ts),
                Err(_) => {
                    out.declined.push((fp, "outer-ring"));
                    continue;
                }
            }
            if let Ok(rings) = arrange::hole_rings(model, fh, fp, inc) {
                for r in rings {
                    emit_ring(&r);
                }
            }
        }
    }
}

/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
#[allow(clippy::too_many_arguments)]
fn trace_on_class(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    wc: usize,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &arrange::EdgePlanes,
    inc_b: &arrange::EdgePlanes,
    canon: &[usize],
) -> Trace {
    let mut out = Trace::default();
    trace_one(
        model,
        a,
        SolidSide::A,
        wc,
        planes,
        surf_ix,
        inc_a,
        canon,
        &mut out,
    );
    trace_one(
        model,
        b,
        SolidSide::B,
        wc,
        planes,
        surf_ix,
        inc_b,
        canon,
        &mut out,
    );
    out
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
    /// Every `(solid, kind)` that produced this one geometric edge. Length 1 when nothing was
    /// coincident.
    pub merged: Vec<(SolidSide, SegKind)>,
}

/// Merge segments that are the **same geometric edge** — same `wall` and same endpoint-triple set
/// (direction-independent) — into one `MergedSeg`, collecting their contributions. Partial overlap
/// (same `wall`, *different* extent — the E5 case) is left alone: those are different edges.
fn merge_coincident(segs: &[Seg]) -> Vec<MergedSeg> {
    // Key an edge by (wall, sorted endpoint pair). Endpoints are canon triples, so the sorted pair
    // is a direction-independent identity.
    let key = |s: &Seg| -> (usize, [[usize; 3]; 2]) {
        let (mut a, mut b) = (s.end[0], s.end[1]);
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        (s.wall, [a, b])
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
                    wall: s.wall,
                    end: s.end,
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

/// The third plane naming a `MergedSeg`'s endpoint on its line `wc ∩ wall`: the canon-triple
/// element that is neither `wc` nor the wall. `None` on a degenerate `[wc, wall, wall]` triple.
fn endpoint_third(end: [usize; 3], wc: usize, wall: usize) -> Option<usize> {
    let mut r = None;
    for x in end {
        if x != wc && x != wall && r.replace(x).is_some() {
            return None; // two off-planes: degenerate naming
        }
    }
    r
}

/// Split the segments on plane class `wc` at their mutual crossings, returning the arrangement's
/// 1-skeleton: segments cut so no crossing lies in an edge's interior. A crossing of two segments
/// riding `wall1`/`wall2` is the plane triple `{wc, wall1, wall2}` — no new point species.
///
/// `Err` (honest reject) on: a degenerate endpoint name, or two crossings coincident on one
/// segment (a four-plane concurrency the DCEL cannot yet represent — `FOURPLANE`, following
/// `ordered_on_edge`). The bool is **`overlap`**: whether any same-`wall` partial overlap (E5) was
/// detected — those segments are not split against each other, so the arrangement is *incomplete*
/// and a caller must not conclude from it.
fn split_at_crossings(
    planes: &[PlaneInfo],
    wc: usize,
    segs: &[MergedSeg],
) -> Result<(Vec<MergedSeg>, bool), BoolError> {
    // `r` is strictly inside segment `s` (between its two endpoint thirds along the line).
    let strictly_inside = |s: &MergedSeg, r: usize| -> Option<bool> {
        let (r0, r1) = (
            endpoint_third(s.end[0], wc, s.wall)?,
            endpoint_third(s.end[1], wc, s.wall)?,
        );
        let (a, b) = (
            arrange::order_along(planes, wc, s.wall, r, r0),
            arrange::order_along(planes, wc, s.wall, r, r1),
        );
        Some(a != 0 && b != 0 && a != b)
    };
    // E5: two segments on the same wall whose extents partially overlap (an endpoint of one lies
    // strictly inside the other). Detected, not resolved.
    let mut overlap = false;
    for (i, s1) in segs.iter().enumerate() {
        for s2 in &segs[i + 1..] {
            if s1.wall != s2.wall {
                continue;
            }
            for (a, b) in [(s1, s2), (s2, s1)] {
                for e in b.end {
                    if let Some(r) = endpoint_third(e, wc, a.wall) {
                        if strictly_inside(a, r) == Some(true) {
                            overlap = true;
                        }
                    }
                }
            }
        }
    }

    let mut out = Vec::new();
    for (i, s) in segs.iter().enumerate() {
        // Per-segment: collect the crossings strictly interior to `s` (re-test, do not rely on a
        // deduped global set — a crossing may split several segments of different extent).
        let mut interior: Vec<usize> = Vec::new();
        for (j, o) in segs.iter().enumerate() {
            if i == j || o.wall == s.wall {
                continue;
            }
            if tolerant::t_plane_pair_dir_sign(planes, wc, s.wall, o.wall) == 0 {
                continue; // walls meet wc in no point (parallel)
            }
            // The crossing of `s` and `o` is named by `o.wall` along `s`'s line.
            if strictly_inside(s, o.wall) == Some(true)
                && strictly_inside(o, s.wall) == Some(true)
                && !interior.contains(&o.wall)
            {
                interior.push(o.wall);
            }
        }
        let (r0, r1) = (
            endpoint_third(s.end[0], wc, s.wall).ok_or_else(|| reject(tag::THREE_PLANES))?,
            endpoint_third(s.end[1], wc, s.wall).ok_or_else(|| reject(tag::THREE_PLANES))?,
        );
        if interior.is_empty() {
            out.push(s.clone());
            continue;
        }
        // Sort interior crossings along the line; two coinciding is a four-plane concurrency.
        interior.sort_by(
            |&a, &b| match arrange::order_along(planes, wc, s.wall, a, b) {
                -1 => std::cmp::Ordering::Less,
                1 => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            },
        );
        for w in interior.windows(2) {
            if arrange::order_along(planes, wc, s.wall, w[0], w[1]) == 0 {
                return Err(reject(tag::FOURPLANE));
            }
        }
        // Emit sub-segments over the R-sequence [r0, interior…, r1], inheriting wall/merged.
        let seq: Vec<usize> = std::iter::once(r0)
            .chain(interior.iter().copied())
            .chain(std::iter::once(r1))
            .collect();
        // `wc`, `s.wall`, and each R are already canon classes; the endpoint triple is just sorted.
        let sorted = |r: usize| {
            let mut t = [wc, s.wall, r];
            t.sort_unstable();
            t
        };
        for w in seq.windows(2) {
            out.push(MergedSeg {
                wall: s.wall,
                end: [sorted(w[0]), sorted(w[1])],
                merged: s.merged.clone(),
            });
        }
    }
    Ok((out, overlap))
}

/// CCW cyclic order of the edges around one arrangement vertex on plane class `w`, **read from no
/// coordinate**. Each edge rides a line `w ∩ fp` and runs in direction `s·(n_w × n_fp)`; it is
/// given as `(fp, s)`. The signed turn between edges i and j is `turn_at`'s atom
/// `s_i·s_j·t_plane_pair_dir_sign(w, fp_i, fp_j)·orient_sign(w)` (arrange.rs:524-556).
///
/// The turn sign is transitive only *within an open half-plane* (span < π), where it reproduces
/// the angle order exactly — the textbook Graham-scan fact. So: bucket every edge by its turn
/// against a reference `r = edges[0]` into the two open half-planes, plus the `0`/`π` pole where
/// the turn is 0 (collinear with `r`). The pole splits by same-`fp` opposite direction: same fp,
/// opposite `s` is angle π; otherwise angle 0. Each open bucket is then a real sort. Returned
/// indices are the CCW order starting at `r`.
///
/// **This brick proves the predicate is buildable on the exact substrate — it does not yet handle
/// the general case:** two *different* fp's whose lines are parallel (a `0`-turn that is not
/// same-fp) fall into the angle-0 bucket unresolved, and collinear same-direction overlap (E5) is
/// out of scope. The corpus's arrangement vertices are degree ≥ 3 with distinct fp's per real
/// direction, which is what the spike exercises.
fn angular_order(planes: &[PlaneInfo], w: usize, edges: &[(usize, i8)]) -> Vec<usize> {
    let os = arrange::orient_sign(planes, w);
    let cross = |i: usize, j: usize| -> i8 {
        edges[i].1
            * edges[j].1
            * tolerant::t_plane_pair_dir_sign(planes, w, edges[i].0, edges[j].0)
            * os
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
    planes: &[PlaneInfo],
    wc: usize,
    segs: &[MergedSeg],
) -> Result<(Vec<Cell>, HashMap<usize, usize>), BoolError> {
    let n = segs.len();
    let he_count = 2 * n;
    let origin = |he: usize| segs[he / 2].end[he % 2]; // he%2==0: end[0]; ==1: end[1]
    let target = |he: usize| segs[he / 2].end[1 - he % 2];
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
            let (rv, rf) = (
                endpoint_third(origin(he), wc, wall(he))
                    .ok_or_else(|| reject(tag::THREE_PLANES))?,
                endpoint_third(target(he), wc, wall(he))
                    .ok_or_else(|| reject(tag::THREE_PLANES))?,
            );
            // Direction sign away from v toward the far end (edge_sign convention).
            let s = arrange::order_along(planes, wc, wall(he), rf, rv);
            if s == 0 {
                return Err(reject(tag::LOOP_ORIENT_MISMATCH));
            }
            edges.push((wall(he), s));
        }
        let ord = angular_order(planes, wc, &edges);
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
            let ring: Vec<[usize; 3]> = cyc.iter().map(|&h| origin(h)).collect();
            let w = arrange::loop_winding(planes, wc, &ring)?;
            cells.push(Cell {
                half_edges: cyc,
                winding: w,
            });
        }
        if ok && cells.iter().filter(|c| c.winding == -1).count() == components {
            return Ok((cells, face_of));
        }
    }
    Err(reject(tag::LOOP_ORIENT_MISMATCH))
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
/// [`arrange::point_in_ring`]. A `c` inside no `+1` cell bounds the unbounded region: the root.
struct Nesting {
    /// `group_of[cell]` = the cell's group representative (the `+1` host for a hole group, else
    /// the cell itself).
    group_of: Vec<usize>,
    /// The group of the single unbounded contour — `label_cells`' seed.
    root_group: usize,
    /// `holes[host]` = the `-1` cells nested in that `+1` host, emitted as its inner rings.
    holes: HashMap<usize, Vec<usize>>,
}

/// Classify every cell as root / hole / plain `+1` (see [`Nesting`]). Single-hole scope: more than
/// one hole in one face (`HOLE_MULTI`), a hole nested deeper than one level (`HOLE_DEPTH`), or more
/// than one unbounded contour — several disjoint bodies on the plane — (`HOLE_ROOTS`) are honest
/// rejects, each a distinct tag so a refactor cannot silently merge them.
fn nest_cells(
    planes: &[PlaneInfo],
    wc: usize,
    cells: &[Cell],
    segs: &[MergedSeg],
) -> Result<Nesting, BoolError> {
    let n = cells.len();
    let ring_of = |c: &Cell| -> Vec<[usize; 3]> {
        c.half_edges
            .iter()
            .map(|&he| segs[he / 2].end[he % 2])
            .collect()
    };
    let rings: Vec<Vec<[usize; 3]>> = cells.iter().map(ring_of).collect();
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
            if rings[c].iter().any(|t| rings[r].contains(t)) {
                continue;
            }
            // Vertex-disjoint: one clear ray settles it. Retry past a spoiled (ring-node) ray;
            // all of `c`'s vertices spoiled against `r` is a genuine degeneracy → honest reject.
            let mut inside = None;
            for &v in &rings[c] {
                if let Ok(hit) = arrange::point_in_ring(planes, wc, v, &rings[r]) {
                    inside = Some(hit);
                    break;
                }
            }
            match inside {
                Some(true) => hosts.push(r),
                Some(false) => {}
                None => return Err(reject(tag::NO_CLEAR_RAY)),
            }
        }
        match hosts.len() {
            0 => roots.push(c),
            1 => {
                let (rc, rr) = (find(&mut parent, c), find(&mut parent, hosts[0]));
                parent[rc] = rr;
                holes.entry(hosts[0]).or_default().push(c);
            }
            _ => return Err(reject(tag::HOLE_DEPTH)), // inside >1 +1 cell: nested deeper
        }
    }
    if roots.len() != 1 {
        return Err(reject(tag::HOLE_ROOTS)); // several disjoint bodies on one plane
    }
    if holes.values().any(|hs| hs.len() > 1) {
        return Err(reject(tag::HOLE_MULTI)); // more than one hole in one face
    }
    let group_of: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    let root_group = group_of[roots[0]];
    Ok(Nesting {
        group_of,
        root_group,
        holes,
    })
}

/// A per-solid, per-side material label of one cell: `[A_above, A_below, B_above, B_below]`.
type Label = [bool; 4];

/// The flip mask an edge applies when crossed, grouping its `merged` contributions **per solid**.
/// Crossing the edge XORs this into the cell label.
///
/// A `Seated` face is ground truth for its solid's material next to W (it directly says which side
/// the body is on), so **seated wins**: if a solid contributes any `Seated`, that fully determines
/// its flip and a same-solid `Transversal` is discarded (the transversal's "flip both" assumes the
/// solid straddles W, which is false exactly where it is capped). A solid with only `Transversal`
/// genuinely straddles → flips both sides. `> 1 Transversal` or disagreeing seated body sides is a
/// coincident-wall degeneracy outside the corpus → honest reject.
fn edge_mask(merged: &[(SolidSide, SegKind)]) -> Result<Label, BoolError> {
    let mut mask = [false; 4];
    for (solid, base) in [(SolidSide::A, 0usize), (SolidSide::B, 2)] {
        let kinds: Vec<SegKind> = merged
            .iter()
            .filter(|(s, _)| *s == solid)
            .map(|(_, k)| *k)
            .collect();
        let seated: Vec<bool> = kinds
            .iter()
            .filter_map(|k| match k {
                SegKind::Seated { body_above } => Some(*body_above),
                _ => None,
            })
            .collect();
        if !seated.is_empty() {
            if seated.iter().any(|&b| b != seated[0]) {
                return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // disagreeing seated sides
            }
            // seated wins: flip above if body_above, else below.
            mask[base + usize::from(!seated[0])] ^= true;
        } else if kinds.len() == 1 {
            // pure transversal: the solid straddles W, flip both.
            mask[base] ^= true;
            mask[base + 1] ^= true;
        } else if kinds.len() > 1 {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // >1 transversal, same solid
        }
        // kinds empty ⇒ solid absent from this edge ⇒ no flip.
    }
    Ok(mask)
}

/// Label every cell by propagating from the unbounded cell (all void) across edges, flipping per
/// `edge_mask`. The cell interior cannot be point-queried (no plane-triple name), so propagation is
/// the only route — the 2D-face analogue of `run_classes`. After propagating, **every** edge's flip
/// relation is verified (`label[c] XOR mask == label[neighbour]`); a violation means the trace was
/// incomplete and is an honest reject.
fn label_cells(
    cells: &[Cell],
    face_of: &HashMap<usize, usize>,
    segs: &[MergedSeg],
    nesting: &Nesting,
) -> Result<Vec<Label>, BoolError> {
    // A face-with-holes is one region: label its group as a unit. Group representative → members.
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &g) in nesting.group_of.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut label = vec![None; cells.len()];
    let mut queue = std::collections::VecDeque::new();
    // Seed the unbounded root group with all-void, enqueuing every member.
    for &i in &members[&nesting.root_group] {
        label[i] = Some([false; 4]);
        queue.push_back(i);
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
        .ok_or_else(|| reject(tag::LOOP_ORIENT_MISMATCH))?; // a cell never reached
    // Verify every edge (tree and non-tree): the flip relation must hold everywhere.
    for (he, &c) in face_of {
        let nb = face_of[&(he ^ 1)];
        let mask = edge_mask(&segs[he / 2].merged)?;
        if std::array::from_fn::<bool, 4, _>(|i| out[c][i] ^ mask[i]) != out[nb] {
            return Err(reject(tag::LOOP_ORIENT_MISMATCH)); // inconsistent propagation
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
/// with an original A/B vertex (a cap corner). The abutting wall faces name that point
/// `Node::Orig`, so the next brick's SeamVertex weld table must canonicalize such a W-triple onto
/// the same result vertex, or `assemble_fuse_cut`'s manifold guard rejects. This brick emits
/// all-`Seam`; the reconciliation and assembly are later.
fn emit_faces(
    kind: BoolKind,
    labels: &[Label],
    cells: &[Cell],
    segs: &[MergedSeg],
    planes: &[PlaneInfo],
    wc: usize,
    holes: &HashMap<usize, Vec<usize>>,
) -> Vec<LocalFace> {
    // `crate::Node` is lib.rs's arrangement node enum; the local `Node` (this module's
    // three-valued-scan struct) shadows it here.
    let ring_of = |cell: &Cell| -> Vec<crate::Node> {
        cell.half_edges
            .iter()
            .map(|&he| crate::Node::Seam(segs[he / 2].end[he % 2]))
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
        let flip = keep_above == (arrange::orient_sign(planes, wc) > 0);
        let inner: Vec<Vec<crate::Node>> = holes
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

/// Every result face across all plane classes, before assembly (the driver's risky half, testable
/// by face count without mutating the model). A declining class aborts the whole boolean.
#[allow(clippy::too_many_arguments)]
fn trace_result_faces(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    planes: &[PlaneInfo],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &arrange::EdgePlanes,
    inc_b: &arrange::EdgePlanes,
    canon: &[usize],
) -> Result<Vec<LocalFace>, BoolError> {
    let mut faces: Vec<LocalFace> = Vec::new();
    for wc in (0..planes.len()).filter(|&i| canon[i] == i) {
        let tr = trace_on_class(model, a, b, wc, planes, surf_ix, inc_a, inc_b, canon);
        if !tr.declined.is_empty() {
            return Err(reject(tag::COPLANAR_PAIR)); // incomplete trace ⇒ honest reject
        }
        let merged = merge_coincident(&tr.segs);
        let (split, overlap) = split_at_crossings(planes, wc, &merged)?;
        if overlap {
            return Err(reject(tag::COPLANAR_OVERLAP_MULTI)); // E5 not resolved yet
        }
        let (cells, face_of) = extract_cells(planes, wc, &split)?;
        let nesting = nest_cells(planes, wc, &cells, &split)?;
        let labels = label_cells(&cells, &face_of, &split, &nesting)?;
        faces.extend(emit_faces(
            kind,
            &labels,
            &cells,
            &split,
            planes,
            wc,
            &nesting.holes,
        ));
    }
    Ok(faces)
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
fn boolean_via_trace(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(model, a, b)?;
    let faces = trace_result_faces(model, kind, a, b, &planes, &surf_ix, &inc_a, &inc_b, &canon)?;

    // Build the SeamVertex weld table directly from the emitted triples (no `build_seam`: that is
    // raw-index and pierce-only). Reject rather than panic on a degenerate meet.
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seen: HashMap<[usize; 3], ()> = HashMap::new();
    for f in &faces {
        for loop_ in std::iter::once(&f.loop_nodes).chain(f.inner.iter()) {
            for node in loop_ {
                let crate::Node::Seam(t) = node else { continue };
                if seen.insert(*t, ()).is_some() {
                    continue;
                }
                let point = three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .ok_or_else(|| reject(tag::THREE_PLANES))?;
                seam.push(SeamVertex {
                    point,
                    triple: *t,
                    tol: vertex_tol(
                        point,
                        &planes[t[0]].plane,
                        &planes[t[1]].plane,
                        &planes[t[2]].plane,
                    ),
                });
            }
        }
    }

    assemble_fuse_cut(model, a, b, &planes, &seam, &faces)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Point of a canon triple, for asserting geometry by hand.
    fn pt(t: [usize; 3], planes: &[PlaneInfo]) -> [f64; 3] {
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();

        // The shared class z=1: the class both solids seat a cap on.
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                let seats =
                    |s: Handle<Solid>| {
                        solid_shell_handles(&m, s).into_iter().any(|sh| {
                            m.shells.get(sh).faces.iter().any(|fh| {
                                canon[surf_ix[fh]] == c && face_on_z1(*fh, &surf_ix, &planes)
                            })
                        })
                    };
                seats(a) && seats(b)
            })
            .expect("a shared cap class");

        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        // Two seated caps × 4 edges each. (The side walls now also trace as transversal chords;
        // the two far caps are parallel to z=1 and contribute nothing — misses, not declines.)
        let seated: Vec<&Seg> = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Seated { .. }))
            .collect();
        assert_eq!(seated.len(), 8, "two seated caps, four edges each: {tr:?}");
        // Four side walls per cube each cross z=1 in one chord.
        let trans = tr
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .count();
        assert_eq!(trans, 8, "eight side-wall chords: {tr:?}");
        assert!(
            tr.declined.is_empty(),
            "axis-aligned cubes decline nothing: {tr:?}"
        );
        assert!(tr.touches.is_empty());
        for s in &seated {
            // On z=1, riding a side wall (not W itself), body on the geometrically correct side.
            for e in &s.end {
                assert!((pt(*e, &planes)[2] - 1.0).abs() < 1e-12, "endpoint on z=1");
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

    /// The transversal chord across a non-convex reflex plane — the cap chord — is produced whole,
    /// collapsing the on-line edge, which the seated brick provably cannot do (it declines this
    /// face). l_prism (profile (0,0),(2,0),(2,1),(1,1),(1,2),(0,2), extruded z∈[0,1]) cut at y=1:
    /// the z=0 cap's trace is the single segment x∈[0,2], **not** [0,1]+[1,2].
    #[test]
    fn the_cap_chord_collapses_the_on_line_edge() {
        let profile = Profile2d {
            points: vec![
                Point2::from_array([0.0, 0.0]),
                Point2::from_array([2.0, 0.0]),
                Point2::from_array([2.0, 1.0]),
                Point2::from_array([1.0, 1.0]),
                Point2::from_array([1.0, 2.0]),
                Point2::from_array([0.0, 2.0]),
            ],
        };
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
        let (planes, surf_ix, inc_a, _inc_b, canon) = plane_index_setup(&m, a2, b2).unwrap();

        // The y=1 class: a plane through all-y=1 points that a's reflex face sits on.
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                planes.iter().enumerate().any(|(i, p)| {
                    canon[i] == c && p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
                })
            })
            .expect("a y=1 class");

        let mut out = Trace::default();
        trace_one(
            &m,
            a2,
            SolidSide::A,
            wc,
            &planes,
            &surf_ix,
            &inc_a,
            &canon,
            &mut out,
        );

        // Collect transversal segments as rounded endpoint-pairs.
        let chords: Vec<([i64; 3], [i64; 3])> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &planes)), rp(pt(s.end[1], &planes)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .collect();

        // The cap chords: the z=0 (and z=1) cap ∩ (y=1) is the single segment x∈[0,2]. A cap is a
        // horizontal segment (both endpoints share z), which distinguishes it from the reflex wall
        // x=1's own vertical chord (1,1,0)-(1,1,1), a legitimate different face's trace.
        let horiz_z0: Vec<_> = chords
            .iter()
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        // ★ Falsifiable core: a failed collapse yields TWO horizontal chords ([0,1]+[1,2], an x=1
        // endpoint appears); a correct collapse yields exactly ONE spanning [0,2].
        assert_eq!(
            horiz_z0.len(),
            1,
            "the on-line edge collapsed to one chord: {chords:?}"
        );
        assert_eq!(
            *horiz_z0[0],
            (rp([0.0, 1.0, 0.0]), rp([2.0, 1.0, 0.0])),
            "cap chord spans x∈[0,2]: {chords:?}"
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
    /// parity yet must still appear inside the trace. The cap ∩ y=1 is the single span x∈[0,3].
    #[test]
    fn a_tangential_on_line_edge_still_spans() {
        let profile = Profile2d {
            points: vec![
                Point2::from_array([0.0, 0.0]),
                Point2::from_array([3.0, 0.0]),
                Point2::from_array([3.0, 2.3]),
                Point2::from_array([2.0, 2.3]),
                Point2::from_array([2.0, 1.0]),
                Point2::from_array([1.0, 1.0]),
                Point2::from_array([1.0, 2.0]),
                Point2::from_array([0.0, 2.0]),
            ],
        };
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
        let (planes, surf_ix, inc_a, _inc_b, canon) = plane_index_setup(&m, a2, b2).unwrap();
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                planes.iter().enumerate().any(|(i, p)| {
                    canon[i] == c && p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
                })
            })
            .expect("a y=1 class");
        let mut out = Trace::default();
        trace_one(
            &m,
            a2,
            SolidSide::A,
            wc,
            &planes,
            &surf_ix,
            &inc_a,
            &canon,
            &mut out,
        );

        let horiz_z0: Vec<_> = out
            .segs
            .iter()
            .filter(|s| matches!(s.kind, SegKind::Transversal { .. }))
            .map(|s| {
                let (mut u, mut v) = (rp(pt(s.end[0], &planes)), rp(pt(s.end[1], &planes)));
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                (u, v)
            })
            .filter(|(u, v)| u[2] == v[2] && u[2] == 0 && u[1] == 1_000_000 && v[1] == 1_000_000)
            .collect();
        assert_eq!(
            horiz_z0.len(),
            1,
            "tangential run bridged into one span: {horiz_z0:?}"
        );
        assert_eq!(
            horiz_z0[0],
            (rp([0.0, 1.0, 0.0]), rp([3.0, 1.0, 0.0])),
            "spans x∈[0,3]"
        );
    }

    /// ★ Spike: a coordinate-free exact cyclic order of edges around an arrangement vertex is
    /// buildable — the "DNA question" winding-engine.md:66/115 flagged as a possible death
    /// condition for B. A degree-5 vertex with directions +x, +y, −x, −y and one **oblique**
    /// (the oblique is non-optional: without it every open half-plane bucket holds one element and
    /// transitivity never fires). `angular_order` returns the CCW order reading no coordinate; we
    /// check it against the CCW order computed *with* coordinates (atan2), which is the oracle.
    #[test]
    fn angular_order_around_a_vertex_is_coordinate_free_and_ccw() {
        // A pentagon prism giving y-, x-, and diagonal-normal side faces plus z caps.
        let profile = Profile2d {
            points: vec![
                Point2::from_array([0.0, 0.0]), // (0,0)-(3,0): y=0
                Point2::from_array([3.0, 0.0]), // (3,0)-(3,3): x=3
                Point2::from_array([3.0, 3.0]), // (3,3)-(2,3): y=3
                Point2::from_array([2.0, 3.0]), // (2,3)-(0,1): diagonal y=x+1
                Point2::from_array([0.0, 1.0]), // (0,1)-(0,0): x=0
            ],
        };
        let m = replay(&[Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist: 1.0,
        }])
        .unwrap();
        let a = m.live_solids[0];
        let planes = collect_planes(&m, a).unwrap();

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

        let order = angular_order(&planes, w, &edges);
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
            let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
            let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
            let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
            let merged = merge_coincident(&tr.segs);
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
            let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
            let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
            let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
            let merged = merge_coincident(&tr.segs);
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        assert_eq!(merged.len(), 8, "8 merged edges before split");

        let (split, overlap) = split_at_crossings(&planes, wc, &merged).unwrap();
        assert!(!overlap, "cross has no same-wall partial overlap");
        // 4 chords crossed twice → 3 pieces each = 12; 4 outer walls uncrossed = 4; total 16.
        assert_eq!(split.len(), 16, "16 sub-segments: {}", split.len());

        // a's y=1 chord splits into exactly 3, ending at x=0,1,2,3 — checked by coordinate.
        let seg_pts = |s: &MergedSeg| {
            let mut u = pt(s.end[0], &planes);
            let mut v = pt(s.end[1], &planes);
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
            let p = pt(t, &planes);
            (p[0] - 1.0).abs() < 1e-9 && (p[1] - 1.0).abs() < 1e-9 && (p[2] - 1.0).abs() < 1e-9
        };
        let incident: Vec<&MergedSeg> = split
            .iter()
            .filter(|s| is_v(s.end[0]) || is_v(s.end[1]))
            .collect();
        assert_eq!(incident.len(), 4, "degree-4 at (1,1,1): {}", incident.len());
        let mut keys = std::collections::HashSet::new();
        for s in &incident {
            let far = if is_v(s.end[0]) { s.end[1] } else { s.end[0] };
            keys.insert((s.wall, endpoint_third(far, wc, s.wall).unwrap()));
        }
        assert_eq!(keys.len(), 4, "no coincident (fp,direction) — DCEL-ready");
    }

    /// Same-wall partial overlap (E5) is detected (not resolved): a=[0,2], b=[1,3] share y=1, their
    /// chords overlap on x∈[1,2]. `overlap` must be true — a vacuous "false" would mean the
    /// detector is blind.
    #[test]
    fn partial_overlap_is_detected() {
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        let (_, overlap) = split_at_crossings(&planes, wc, &merged).unwrap();
        assert!(
            overlap,
            "a's x∈[0,2] and b's x∈[1,3] on y=1 overlap — must be flagged"
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        let (split, _) = split_at_crossings(&planes, wc, &merged).unwrap();

        let (cells, face_of) = extract_cells(&planes, wc, &split).unwrap();
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        let (split, _) = split_at_crossings(&planes, wc, &merged).unwrap();
        let (cells, face_of) = extract_cells(&planes, wc, &split).unwrap();
        let nesting = nest_cells(&planes, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting).unwrap();

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
                vs.iter().map(|&t| pt(t, &planes)).collect()
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        // a passes through z=1 (no seated face); find the class b caps at z=1.
        let wc = (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                solid_shell_handles(&m, b).into_iter().any(|sh| {
                    m.shells
                        .get(sh)
                        .faces
                        .iter()
                        .any(|fh| canon[surf_ix[fh]] == c && face_on_z1(*fh, &surf_ix, &planes))
                })
            })
            .expect("b's z=1 cap class");
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        let (split, _) = split_at_crossings(&planes, wc, &merged).unwrap();
        let (cells, face_of) = extract_cells(&planes, wc, &split).unwrap();
        let nesting = nest_cells(&planes, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting).unwrap();

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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        let merged = merge_coincident(&tr.segs);
        let (split, _) = split_at_crossings(&planes, wc, &merged).unwrap();
        let (cells, face_of) = extract_cells(&planes, wc, &split).unwrap();
        let nesting = nest_cells(&planes, wc, &cells, &split).unwrap();
        let labels = label_cells(&cells, &face_of, &split, &nesting).unwrap();

        // Centroid of a face's ring (convex cells here).
        let centroid = |f: &LocalFace| -> [f64; 2] {
            let ps: Vec<[f64; 3]> = f
                .loop_nodes
                .iter()
                .map(|n| match n {
                    crate::Node::Seam(t) => pt(*t, &planes),
                    _ => unreachable!("all-Seam"),
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
            &planes,
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
            &planes,
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
            &planes,
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
                .map(|n| match n {
                    crate::Node::Seam(t) => *t,
                    _ => unreachable!(),
                })
                .collect();
            assert!(ring.len() >= 3);
            assert_eq!(
                arrange::loop_winding(&planes, wc, &ring).unwrap(),
                1,
                "CCW about n_out"
            );
        }

        // flip oracle (coordinate, test-only): the result normal points away from the kept
        // chamber. Fuse keeps below (bodies below the cap), so n_result·n_w > 0.
        let n_w = planes[wc].plane.normal();
        let os = arrange::orient_sign(&planes, wc) as f64;
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
                .map(|n| match n {
                    crate::Node::Seam(t) => *t,
                    _ => unreachable!(),
                })
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let faces = trace_result_faces(
            &m,
            BoolKind::Fuse,
            a,
            b,
            &planes,
            &surf_ix,
            &inc_a,
            &inc_b,
            &canon,
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
                    crate::Node::Seam(t) => *t,
                    _ => unreachable!("all-Seam"),
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
        let solids = boolean_via_trace(&mut m, BoolKind::Fuse, a, b).unwrap();
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
            let solids = boolean_via_trace(&mut m, kind, a, b).unwrap();
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
        let solids = boolean_via_trace(&mut m, kind, a, b).unwrap();
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let faces = trace_result_faces(
            &m,
            BoolKind::Cut,
            a,
            b,
            &planes,
            &surf_ix,
            &inc_a,
            &inc_b,
            &canon,
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
        let triples = |ns: &[crate::Node]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|n| match n {
                    crate::Node::Seam(t) => *t,
                    _ => unreachable!("all-Seam"),
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
    // solid `Origin::Rotated`, routing every predicate to the exact frame3 backend. `rot30` (lib.rs)
    // is in a sibling test module and unreachable here, so the isometries are built inline.

    /// 30° about `axis` through (1,1,0) — non-90°, so `Origin::Rotated` (exact frame3 path).
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
            let solids = boolean_via_trace(&mut m, kind, a, b).unwrap();
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
    /// `t_orient3d` on-plane fix), so all four orientations returning 24 is the fix's direct
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
            let solids = boolean_via_trace(&mut m, BoolKind::Cut, a, b).unwrap();
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let faces = trace_result_faces(
            &m,
            BoolKind::Cut,
            a,
            b,
            &planes,
            &surf_ix,
            &inc_a,
            &inc_b,
            &canon,
        )
        .unwrap();
        assert_eq!(faces.len(), 10, "rotated arrangement keeps 10 faces");
        assert_eq!(
            faces.iter().filter(|f| !f.inner.is_empty()).count(),
            2,
            "the two annular caps survive rotation"
        );
        let triples = |ns: &[crate::Node]| -> Vec<[usize; 3]> {
            ns.iter()
                .map(|n| match n {
                    crate::Node::Seam(t) => *t,
                    _ => unreachable!("all-Seam"),
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
        let solids = boolean_via_trace(&mut m, BoolKind::Cut, a, b).unwrap();
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
    /// `t_orient3d` on-plane fix is reverted.
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        for wc in {
            let mut c: Vec<usize> = canon.clone();
            c.sort_unstable();
            c.dedup();
            c
        } {
            let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
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
        planes: &[PlaneInfo],
        canon: &[usize],
    ) -> usize {
        (0..planes.len())
            .map(|i| canon[i])
            .find(|&c| {
                let seats =
                    |s: Handle<Solid>| {
                        solid_shell_handles(m, s).into_iter().any(|sh| {
                            m.shells.get(sh).faces.iter().any(|fh| {
                                canon[surf_ix[fh]] == c && face_on_z1(*fh, surf_ix, planes)
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        let wc = shared_cap_class(&m, a, b, &surf_ix, &planes, &canon);
        let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
        // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
        // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
        // seated≡transversal coincidence).
        let merged = merge_coincident(&tr.segs);
        // The shared y=1 wall class (a face at y=1).
        let y1 = canon[planes
            .iter()
            .position(|p| p.tri.iter().all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12))
            .expect("a y=1 face")];
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
        planes: &[PlaneInfo],
    ) -> bool {
        let p = &planes[surf_ix[&fh]];
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
        let (planes, surf_ix, inc_a, inc_b, canon) = plane_index_setup(&m, a, b).unwrap();
        for wc in {
            let mut c: Vec<usize> = canon.clone();
            c.sort_unstable();
            c.dedup();
            c
        } {
            let tr = trace_on_class(&m, a, b, wc, &planes, &surf_ix, &inc_a, &inc_b, &canon);
            assert!(tr.declined.is_empty(), "class {wc} declined: {tr:?}");
        }
    }
}
