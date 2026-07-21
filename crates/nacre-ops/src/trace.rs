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
            // `n_out · n_W < 0`. (f64 sign is exact for axis-aligned normals; rotated inputs are
            // outside this corpus and this brick does not claim them.)
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
) -> Result<Vec<Label>, BoolError> {
    let seed = cells
        .iter()
        .position(|c| c.winding == -1)
        .ok_or_else(|| reject(tag::LOOP_ORIENT_MISMATCH))?;
    let mut label = vec![None; cells.len()];
    label[seed] = Some([false; 4]);
    let mut queue = std::collections::VecDeque::from([seed]);
    while let Some(c) = queue.pop_front() {
        let lc = label[c].unwrap();
        for &he in &cells[c].half_edges {
            let nb = face_of[&(he ^ 1)];
            if label[nb].is_none() {
                let mask = edge_mask(&segs[he / 2].merged)?;
                label[nb] = Some(std::array::from_fn(|i| lc[i] ^ mask[i]));
                queue.push_back(nb);
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
        let labels = label_cells(&cells, &face_of, &split).unwrap();

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
        let labels = label_cells(&cells, &face_of, &split).unwrap();

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
