//! **A symbolic order for one sanctioned coincident pair — the twins a bridge or a pinch makes.**
//!
//! A face's boundary can pass one point twice. When a hole's ring is spliced into the ring it
//! touches, the touching vertex appears twice in the merged ring: once as the split ring's
//! vertex, once as the hole's, at the same coordinate. When one ring passes a point twice — a
//! pinch, a groove whose tip touches the face's own rim — it already does. Every exact predicate
//! the sweep asks about that pair answers zero, and the sweep has no order for a zero. The
//! classic answer (Edelsbrunner–Mücke) is to *simulate* a displacement: decide every tie as if
//! each twin had moved by an infinitesimal `ε` in a direction of its own, and read the sign off
//! the first-order term. No coordinate moves — the triangles come out on the exact points — only
//! the ties are decided, consistently, as one actual small displacement would decide them.
//!
//! ★★★★★ **Each twin slides along one of its own edges**: `d_o = uv[o_along] − T`, `d_h =
//! uv[h_along] − T`. On a bridge the two are the ends of the split straight segment — the
//! split ring's copy slides toward its predecessor, the hole's toward its successor — so both
//! twins stay *on* that line for every `ε` and keep every collinearity they had with the rest of
//! the ring. On a pinch the earlier visit slides along its incoming edge and the later along its
//! outgoing one, each into its own material wedge, which is what separates them
//! (`polygon::pinch` checks the wedges allow it).
//!
//! ★★ **The order needs nothing more — the sweep's own tie-break carries the rest.** The sweep's
//! order (`v` descending, then `u` ascending) is itself a perturbation: the sweep line turned by
//! an infinitesimal `δ`. The twins' displacement is taken infinitely smaller still, `ε ≪ δ`, so
//! a twin against a non-twin at the same `v` is decided by `u` exactly as any two points are, and
//! only the twins' own tie needs `ε`. `orient2d` does not turn with the sweep, so `δ` never enters
//! a side test. A displacement with a normal component is therefore no harm to the order; what
//! the tangential slide buys a bridge is the kept collinearity.
//!
//! ★ **Twin-only.** Two facts keep this from touching anything else. Every zero the sweep can
//! meet without the twins is already a refusal (`self_touch`: a vertex on another edge's open
//! segment), so once the duplicate-coordinate guard is kept for every pair but the sanctioned
//! one, the *only* new zeros are those with a twin in them. And the two twins are the only
//! points that share a coordinate, so `d` is zero for every other index and the first-order
//! term reduces to a handful of exact products.
//!
//! **The first-order term.** `orient2d` is affine in each argument, so
//! `orient2d(a+εd_a, b+εd_b, c+εd_c) = orient2d(a,b,c) + ε·[(b−a)×(d_c−d_a) + (d_b−d_a)×(c−a)] + O(ε²)`.
//! It is evaluated with Shewchuk expansions, exactly. The lexicographic tie is simpler still:
//! the twins share `T`, so `d_o − d_h = uv[o_along] − uv[h_along]`, a difference of two stored
//! points that is never zero — the twins' order is `lex_less` of the two neighbours they slide
//! toward.
//!
//! What is *not* symbolic: `classify`'s turn (a zero there is a real spike), the touch
//! detector and its witness (a zero there *is* the touch), the crossing scan, `emit`'s
//! zero-area refusal, and `delaunay::refine`'s convexity test. Each of those must go on
//! answering about the real numbers.

use super::{P2, lex_less};
use nacre_predicates::{Expansion, orient2d};

/// The sanctioned pair: two indices at one coordinate, each with the index of the ring
/// neighbour it slides toward (on a bridge the two ends of the split segment).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Twins {
    /// The split ring's copy of the touching vertex, or a pinch's earlier visit; slides toward
    /// `o_along` (its `prev`).
    pub(super) o: usize,
    /// The hole's copy, or a pinch's later visit; slides toward `h_along` (its `next`).
    pub(super) h: usize,
    pub(super) o_along: usize,
    pub(super) h_along: usize,
}

impl Twins {
    fn is_twin(&self, i: usize) -> bool {
        i == self.o || i == self.h
    }

    /// Whether `a` and `b` are the two twins (in either order).
    pub(super) fn is_pair(&self, a: usize, b: usize) -> bool {
        (a == self.o && b == self.h) || (a == self.h && b == self.o)
    }

    /// The index `i` slides toward, if it is a twin.
    fn along(&self, i: usize) -> Option<usize> {
        if i == self.o {
            Some(self.o_along)
        } else if i == self.h {
            Some(self.h_along)
        } else {
            None
        }
    }

    /// `d_i` as exact per-component expansions: `uv[along] − uv[i]` for a twin, zero otherwise.
    fn displacement(&self, i: usize, uv: &[P2]) -> [Expansion; 2] {
        match self.along(i) {
            Some(j) => [
                Expansion::two_diff(uv[j][0], uv[i][0]),
                Expansion::two_diff(uv[j][1], uv[i][1]),
            ],
            None => [Expansion::two_diff(0.0, 0.0), Expansion::two_diff(0.0, 0.0)],
        }
    }
}

/// `lex_less(uv[a], uv[b])`, with the twins ordered as their displacement orders them.
///
/// Distinct coordinates are compared exactly, as any two points are. The twins' tie is decided by
/// where they slide: `T + ε(u − T)` precedes `T + ε(v − T)` exactly when `u` precedes `v`, so the
/// answer is `lex_less` of the two neighbours they slide toward — a comparison of stored
/// coordinates, never a tie, because only the twins share a coordinate. `None` is a coincident
/// pair that is not the sanctioned one.
pub(super) fn lex_less_idx(a: usize, b: usize, uv: &[P2], twins: Option<&Twins>) -> Option<bool> {
    if a == b {
        return Some(false);
    }
    if uv[a] != uv[b] {
        return Some(lex_less(uv[a], uv[b]));
    }
    let tw = twins?;
    if !tw.is_pair(a, b) {
        return None;
    }
    let (ua, ub) = (tw.along(a)?, tw.along(b)?);
    Some(lex_less(uv[ua], uv[ub]))
}

/// `side(a, b, c)` — the sign of `orient2d` — with a zero that involves a twin decided by the
/// first-order term of the displacement.
///
/// A non-zero determinant is returned as is. A zero with **no** twin among the three points is
/// returned as `Some(0)`: it is one of today's zeros (a collinear triple) and every caller
/// already has a rule for it. A zero with a twin in it is the structural `orient2d(p, q, p)` the
/// twins create; its first-order term is exact, and `None` says even that vanished — the one
/// case this order does not decide, which a caller turns into a refusal rather than a guess.
pub(super) fn side_idx(
    a: usize,
    b: usize,
    c: usize,
    uv: &[P2],
    twins: Option<&Twins>,
) -> Option<i8> {
    let det = orient2d(uv[a], uv[b], uv[c]);
    if det > 0.0 {
        return Some(1);
    }
    if det < 0.0 {
        return Some(-1);
    }
    let Some(tw) = twins else {
        return Some(0);
    };
    if !(tw.is_twin(a) || tw.is_twin(b) || tw.is_twin(c)) {
        return Some(0);
    }
    // (b−a)×(d_c−d_a) + (d_b−d_a)×(c−a), exactly.
    let (da, db, dc) = (
        tw.displacement(a, uv),
        tw.displacement(b, uv),
        tw.displacement(c, uv),
    );
    let ba = [
        Expansion::two_diff(uv[b][0], uv[a][0]),
        Expansion::two_diff(uv[b][1], uv[a][1]),
    ];
    let ca = [
        Expansion::two_diff(uv[c][0], uv[a][0]),
        Expansion::two_diff(uv[c][1], uv[a][1]),
    ];
    let dca = [dc[0].sub(&da[0]), dc[1].sub(&da[1])];
    let dba = [db[0].sub(&da[0]), db[1].sub(&da[1])];
    let cross = |p: &[Expansion; 2], q: &[Expansion; 2]| p[0].mul(&q[1]).sub(&p[1].mul(&q[0]));
    let term = cross(&ba, &dca).add(&cross(&dba, &ca));
    match term.sign() {
        0 => None,
        s => Some(s),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tw(o: usize, h: usize, o_along: usize, h_along: usize) -> Twins {
        Twins {
            o,
            h,
            o_along,
            h_along,
        }
    }

    /// Distinct coordinates: exactly `lex_less`, twins or not.
    #[test]
    fn distinct_points_are_compared_exactly() {
        let uv = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]];
        let t = tw(2, 3, 0, 1);
        for a in 0..2 {
            for b in 0..4 {
                if uv[a] != uv[b] {
                    assert_eq!(
                        lex_less_idx(a, b, &uv, Some(&t)),
                        Some(lex_less(uv[a], uv[b]))
                    );
                    assert_eq!(lex_less_idx(a, b, &uv, None), Some(lex_less(uv[a], uv[b])));
                }
            }
        }
    }

    /// The twins' order is their on-line neighbours' order: on a vertical segment the copy that
    /// slides up comes first; on a horizontal one the copy that slides left comes first.
    #[test]
    fn the_twins_order_as_their_neighbours_do() {
        // vertical segment u=(0,4) → v=(0,0), T=(0,2); o slides toward u, h toward v
        let uv = [[0.0, 4.0], [0.0, 0.0], [0.0, 2.0], [0.0, 2.0]];
        let t = tw(2, 3, 0, 1);
        assert_eq!(lex_less_idx(2, 3, &uv, Some(&t)), Some(true));
        assert_eq!(lex_less_idx(3, 2, &uv, Some(&t)), Some(false));
        // horizontal segment u=(0,0) → v=(4,0), T=(2,0); smaller x first
        let uv = [[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 0.0]];
        let t = tw(2, 3, 0, 1);
        assert_eq!(lex_less_idx(2, 3, &uv, Some(&t)), Some(true));
        assert_eq!(lex_less_idx(3, 2, &uv, Some(&t)), Some(false));
        // a coincident pair that is not the sanctioned one has no order
        assert_eq!(lex_less_idx(2, 3, &uv, None), None);
        assert_eq!(lex_less_idx(2, 3, &uv, Some(&tw(0, 1, 2, 3))), None);
    }

    /// The structural zero `side(T_o, h1, T_h)` — the later twin's `insert` against the earlier
    /// twin's edge — resolves at first order: `(h1−T)×(v−u)` here, and flipping the twins'
    /// directions flips it.
    #[test]
    fn a_structural_zero_resolves_at_first_order() {
        // left edge u=(0,4) → v=(0,0), T=(0,2); h1=(0.5,1.5) on v's side
        let uv = [[0.0, 4.0], [0.0, 0.0], [0.0, 2.0], [0.0, 2.0], [0.5, 1.5]];
        let (o, h, h1) = (2, 3, 4);
        assert_eq!(orient2d(uv[o], uv[h1], uv[h]), 0.0, "exactly zero today");
        let t = tw(o, h, 0, 1);
        // (h1−T)×(v−u) = (0.5,−0.5)×(0,−4) = −2
        assert_eq!(side_idx(o, h1, h, &uv, Some(&t)), Some(-1));
        // swapping which end each twin slides toward flips the sign
        let t = tw(o, h, 1, 0);
        assert_eq!(side_idx(o, h1, h, &uv, Some(&t)), Some(1));
        // a zero with no twin in it stays today's zero
        let uv2 = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [5.0, 5.0], [5.0, 5.0]];
        assert_eq!(side_idx(0, 1, 2, &uv2, Some(&tw(3, 4, 0, 1))), Some(0));
        assert_eq!(side_idx(0, 1, 2, &uv2, None), Some(0));
        // a non-zero determinant is untouched
        assert_eq!(side_idx(0, 1, 3, &uv2, Some(&tw(3, 4, 0, 1))), Some(1));
    }

    /// A twin queried against an edge lying along its own line stays collinear after sliding —
    /// the first-order term is zero too, and that is reported, not guessed.
    #[test]
    fn a_zero_that_survives_the_displacement_is_reported() {
        // u=(0,4), v=(0,0), T twins at (0,2), w=(0,3) on the same line
        let uv = [[0.0, 4.0], [0.0, 0.0], [0.0, 2.0], [0.0, 2.0], [0.0, 3.0]];
        let t = tw(2, 3, 0, 1);
        assert_eq!(side_idx(0, 2, 4, &uv, Some(&t)), None);
    }
}
