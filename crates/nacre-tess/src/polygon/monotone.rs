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

use super::{P2, lex_less};
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
pub(super) fn decompose(uv: &[P2], rings: &[&[usize]]) -> Result<Vec<Vec<usize>>, TessError> {
    let n = uv.len();
    let (prev, next) = link(n, rings)?;

    let kinds = classify(uv, &prev, &next)?;
    let diagonals = sweep(uv, &prev, &next, &kinds)?;
    trace(uv, &next, &diagonals)
}

/// Ring order as `prev`/`next` arrays, and the precondition check that makes them
/// well defined: **every vertex is used exactly once.** A repeated index would give
/// one vertex two successors, and the sweep would have no way to know which chain it
/// was on.
fn link(n: usize, rings: &[&[usize]]) -> Result<(Vec<usize>, Vec<usize>), TessError> {
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
) -> Result<Vec<[usize; 2]>, TessError> {
    let n = uv.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        if lex_less(uv[a], uv[b]) {
            std::cmp::Ordering::Less
        } else if lex_less(uv[b], uv[a]) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    // Two vertices at the same point make `lex_less` a non-strict order, and every
    // "which is above" question below would then have no answer.
    if order.windows(2).any(|w| uv[w[0]] == uv[w[1]]) {
        return Err(TessError::DegenerateRing);
    }

    // Sorted left to right by where each edge crosses the sweep line. Edges in the
    // status never cross, so an order established at insertion stays correct.
    let mut status: Vec<Edge> = Vec::new();
    let mut out: Vec<[usize; 2]> = Vec::new();

    // The edge of `status` immediately left of `v` — the one whose region `v` falls
    // into. `None` means `v` is outside every span, which a valid polygon cannot be.
    let left_of = |status: &Vec<Edge>, v: usize, uv: &[P2]| -> Option<usize> {
        status
            .iter()
            .rposition(|e| side(uv[e.upper], uv[next[e.upper]], uv[v]) > 0)
    };
    let insert = |status: &mut Vec<Edge>, upper: usize, helper: usize, uv: &[P2]| {
        let at = status
            .iter()
            .position(|e| side(uv[e.upper], uv[next[e.upper]], uv[upper]) < 0)
            .unwrap_or(status.len());
        status.insert(at, Edge { upper, helper });
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
            Kind::Start => insert(&mut status, v, v, uv),
            Kind::End => {
                let h = remove(&mut status, prev[v])?;
                if kinds[h] == Kind::Merge {
                    out.push([v, h]);
                }
            }
            Kind::Split => {
                let at = left_of(&status, v, uv).ok_or(TessError::DegenerateRing)?;
                out.push([v, status[at].helper]);
                status[at].helper = v;
                insert(&mut status, v, v, uv);
            }
            Kind::Merge => {
                let h = remove(&mut status, prev[v])?;
                if kinds[h] == Kind::Merge {
                    out.push([v, h]);
                }
                let at = left_of(&status, v, uv).ok_or(TessError::DegenerateRing)?;
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
                    insert(&mut status, v, v, uv);
                } else {
                    let at = left_of(&status, v, uv).ok_or(TessError::DegenerateRing)?;
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
) -> Result<(), TessError> {
    let m = piece.len();
    if m < 3 {
        return Err(TessError::DegenerateRing);
    }
    let (mut top, mut bot) = (0, 0);
    for i in 1..m {
        if lex_less(uv[piece[i]], uv[piece[top]]) {
            top = i;
        }
        if lex_less(uv[piece[bot]], uv[piece[i]]) {
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
            lex_less(uv[piece[left[li]]], uv[piece[right[ri]]])
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
        let pieces = decompose(&uv, rings).expect("decomposes");
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
            decompose(&uv(&pts), &[&[0, 1, 2, 3], &[4, 5, 0]]),
            Err(TessError::DegenerateRing)
        ));
    }
}
