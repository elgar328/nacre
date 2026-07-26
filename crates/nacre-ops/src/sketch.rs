//! Turning loose rings into profiles — the nesting policy on top of geom's exact predicates.
//!
//! A sketch is drawn as closed paths; which of them are material, which are holes, and which are
//! islands inside holes is **not** declared by the author but decided from containment, the way
//! every CAD sketcher does it. The classification itself (is this point inside that ring? do two
//! rings meet?) is a robustness-sensitive sign question and lives in `nacre-geom::intersect`
//! (design §1 isolates that); what is here is the policy — depth parity — and the assembly.

use crate::Profile2d;
use nacre_geom::intersect::{RingSide, point_in_ring_2d, ring_self_intersection, rings_cross};
use nacre_math::Point2;

/// Why a set of rings or edges is not a valid set of profiles.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum SketchError {
    /// A ring with fewer than three points encloses nothing.
    DegenerateRing { ring: usize },
    /// Two rings touch or cross. A hole must lie strictly inside its outer ring and strictly
    /// outside its siblings; anything else has no unambiguous inside.
    RingsMeet { a: usize, b: usize },
    /// An edge starts where it ends.
    ZeroLengthEdge { edge: usize },
    /// The same edge appears twice (same endpoints, either way round).
    DuplicateEdge { a: usize, b: usize },
    /// A chain ran out of edges before closing. `at` is the endpoint left dangling and `gap` the
    /// distance to the nearest other free endpoint — reported because "you meant to close this
    /// and missed by 1e-9" is the likely story, and the kernel will not close it for you: an
    /// endpoint either *is* the same point or is not (overview 절대원칙 4).
    OpenChain { at: [f64; 2], gap: Option<f64> },
    /// Three or more edges meet at one point, so the chain has no unambiguous continuation.
    BranchingVertex { at: [f64; 2] },
    /// A ring meets itself — it crosses, touches, or doubles back over its own edge. Such a ring
    /// has no unambiguous inside, so it cannot be sorted into a profile at all.
    ///
    /// Reported as the two offending edges' **midpoints**, not as indices: [`from_edges`] chains
    /// the edges into rings in walk order, so a ring's edge index says nothing about where the
    /// author's input went wrong. Points are what the rest of this enum reports too.
    RingSelfIntersects { ring: usize, at: [[f64; 2]; 2] },
}

/// Sort closed rings into profiles by containment depth.
///
/// Depth is how many rings a ring lies inside. **Even depth is material, odd depth is a hole** —
/// so a ring inside a hole is material again (an island), and each becomes a profile of its own.
/// A hole is attached to its *immediate* container, the deepest ring that contains it.
///
/// There is no fill-rule parameter. Even-odd is the rule; offering the choice would mean asking
/// callers to understand a distinction the syntax deliberately hides, and for the simple,
/// non-touching rings this accepts, even-odd and non-zero agree anyway.
pub fn from_rings(rings: Vec<Vec<Point2>>) -> Result<Vec<Profile2d>, SketchError> {
    for (i, r) in rings.iter().enumerate() {
        if r.len() < 3 {
            return Err(SketchError::DegenerateRing { ring: i });
        }
    }
    // Simplicity comes before nesting, and not merely for a better message: `point_in_ring_2d`
    // decides containment by even-odd parity, which only means "inside" on a simple ring.
    for (i, r) in rings.iter().enumerate() {
        if let Some((a, b)) = ring_self_intersection(r) {
            let mid = |e: usize| {
                let (p, q) = (r[e], r[(e + 1) % r.len()]);
                [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0]
            };
            return Err(SketchError::RingSelfIntersects {
                ring: i,
                at: [mid(a), mid(b)],
            });
        }
    }
    let n = rings.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if rings_cross(&rings[i], &rings[j]) {
                return Err(SketchError::RingsMeet { a: i, b: j });
            }
        }
    }

    // `inside[i]` = the rings containing ring `i`. The rings are disjoint by the check above, so
    // any vertex of `i` answers for the whole ring.
    let inside: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && point_in_ring_2d(rings[i][0], &rings[j]) == RingSide::Inside)
                .collect()
        })
        .collect();

    // A ring's parent is the deepest ring containing it; holes hang off the material ring they
    // are cut from, not off some distant ancestor.
    let depth = |i: usize| inside[i].len();
    let parent: Vec<Option<usize>> = (0..n)
        .map(|i| inside[i].iter().copied().max_by_key(|&j| depth(j)))
        .collect();

    let mut out = Vec::new();
    for i in 0..n {
        if depth(i) % 2 != 0 {
            continue; // odd depth: a hole, carried by its parent below
        }
        let holes: Vec<Vec<Point2>> = (0..n)
            .filter(|&h| depth(h) % 2 == 1 && parent[h] == Some(i))
            .map(|h| rings[h].clone())
            .collect();
        out.push(Profile2d::with_holes(rings[i].clone(), holes));
    }
    Ok(out)
}

/// A 2-D sketch curve. Only straight segments are wired; the enum exists now so that adding arcs
/// with the curved-geometry milestone extends the API instead of breaking it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Curve2d {
    Line,
}

/// One drawn segment: its curve and its two endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge2d {
    pub curve: Curve2d,
    pub start: Point2,
    pub end: Point2,
}

impl Edge2d {
    /// A straight segment.
    pub fn line(start: Point2, end: Point2) -> Edge2d {
        Edge2d {
            curve: Curve2d::Line,
            start,
            end,
        }
    }
}

/// Chain loose edges into closed rings, then sort those into profiles ([`from_rings`]).
///
/// The edges may arrive in any order and pointing either way — this is the shape a generator or
/// an imported file produces. Endpoints must **coincide exactly**; nothing here snaps a near-miss
/// shut. Tolerance is for intersections the kernel *discovers*, never for what a caller
/// constructs (overview 절대원칙 4), and a sketch that silently welds a 1e-9 gap is a sketch
/// whose author does not know what they built. A dangling endpoint is reported with the distance
/// to the nearest free endpoint, which is the number needed to fix it.
pub fn from_edges(edges: Vec<Edge2d>) -> Result<Vec<Profile2d>, SketchError> {
    for (i, e) in edges.iter().enumerate() {
        if e.start == e.end {
            return Err(SketchError::ZeroLengthEdge { edge: i });
        }
    }
    for i in 0..edges.len() {
        for j in (i + 1)..edges.len() {
            let (a, b) = (&edges[i], &edges[j]);
            if (a.start == b.start && a.end == b.end) || (a.start == b.end && a.end == b.start) {
                return Err(SketchError::DuplicateEdge { a: i, b: j });
            }
        }
    }

    // Endpoint → the edges touching it. Exactly two is a corner; more has no single continuation.
    let mut at: Vec<(Point2, Vec<usize>)> = Vec::new();
    let slot = |at: &mut Vec<(Point2, Vec<usize>)>, p: Point2| -> usize {
        match at.iter().position(|(q, _)| *q == p) {
            Some(k) => k,
            None => {
                at.push((p, Vec::new()));
                at.len() - 1
            }
        }
    };
    let ends: Vec<(usize, usize)> = edges
        .iter()
        .map(|e| (slot(&mut at, e.start), slot(&mut at, e.end)))
        .collect();
    for (i, (s, t)) in ends.iter().enumerate() {
        at[*s].1.push(i);
        at[*t].1.push(i);
    }
    // Most specific first: a branch says *which* junction is ambiguous, while a dangling end only
    // says the outline is open — and a branch always leaves an odd end somewhere, so checking in
    // the other order would report the vaguer of the two. (`check_result_topology` orders its
    // defects the same way.)
    if let Some((p, _)) = at.iter().find(|(_, touching)| touching.len() > 2) {
        return Err(SketchError::BranchingVertex { at: p.as_array() });
    }
    if let Some((p, _)) = at.iter().find(|(_, touching)| touching.len() == 1) {
        return Err(SketchError::OpenChain {
            at: p.as_array(),
            gap: nearest_free_gap(*p, &at),
        });
    }

    // Walk each cycle: from an unused edge, hop endpoint to endpoint until back at the start.
    let mut used = vec![false; edges.len()];
    let mut rings: Vec<Vec<Point2>> = Vec::new();
    for seed in 0..edges.len() {
        if used[seed] {
            continue;
        }
        let (first, mut here) = ends[seed];
        let (mut edge, mut ring) = (seed, vec![at[first].0]);
        loop {
            used[edge] = true;
            ring.push(at[here].0);
            let Some(&next) = at[here].1.iter().find(|&&e| !used[e]) else {
                break;
            };
            let (s, t) = ends[next];
            here = if s == here { t } else { s };
            edge = next;
        }
        // The walk returns to its start, so the last point repeats the first.
        if ring.last() == ring.first() {
            ring.pop();
        }
        rings.push(ring);
    }
    from_rings(rings)
}

/// The distance from `p` to the nearest *other* endpoint that is also dangling — the size of the
/// gap the author probably meant to close.
fn nearest_free_gap(p: Point2, at: &[(Point2, Vec<usize>)]) -> Option<f64> {
    at.iter()
        .filter(|(q, touching)| *q != p && touching.len() == 1)
        .map(|(q, _)| {
            let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
            (dx * dx + dy * dy).sqrt()
        })
        .min_by(|a, b| a.partial_cmp(b).unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sq(a: f64, b: f64) -> Vec<Point2> {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    }

    #[test]
    fn a_lone_ring_is_one_profile() {
        let p = from_rings(vec![sq(0.0, 4.0)]).unwrap();
        assert_eq!(p.len(), 1);
        assert!(p[0].inners().is_empty());
    }

    /// The ring order must not matter: the hole is recognised by containment, not by position.
    #[test]
    fn a_ring_inside_a_ring_is_a_hole_either_way() {
        for rings in [
            vec![sq(0.0, 4.0), sq(1.0, 3.0)],
            vec![sq(1.0, 3.0), sq(0.0, 4.0)],
        ] {
            let p = from_rings(rings).unwrap();
            assert_eq!(p.len(), 1, "one body");
            assert_eq!(p[0].inners().len(), 1, "one hole");
            assert_eq!(p[0].outer().len(), 4);
        }
    }

    /// Depth 2 is material again — the island is its own body, and it does *not* become a hole of
    /// the outer ring.
    #[test]
    fn a_ring_inside_a_hole_is_an_island() {
        let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0)]).unwrap();
        assert_eq!(p.len(), 2, "outer body + island");
        let outer = p.iter().find(|q| q.outer()[0][0] == 0.0).unwrap();
        let island = p.iter().find(|q| q.outer()[0][0] == 2.0).unwrap();
        assert_eq!(outer.inners().len(), 1, "the depth-1 ring is its hole");
        assert!(island.inners().is_empty(), "nothing inside the island");
    }

    /// A hole belongs to the ring it is actually cut from, not to a distant ancestor.
    #[test]
    fn a_hole_attaches_to_its_immediate_container() {
        // outer(0..9) ⊃ hole(1..8) ⊃ island(2..7) ⊃ island-hole(3..6)
        let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0), sq(3.0, 6.0)]).unwrap();
        assert_eq!(p.len(), 2);
        let island = p.iter().find(|q| q.outer()[0][0] == 2.0).unwrap();
        assert_eq!(island.inners().len(), 1, "the depth-3 ring is the island's");
        assert_eq!(island.inners()[0][0][0], 3.0);
    }

    /// Two separate bodies, neither inside the other.
    #[test]
    fn disjoint_rings_are_separate_bodies() {
        let p = from_rings(vec![sq(0.0, 1.0), sq(5.0, 6.0)]).unwrap();
        assert_eq!(p.len(), 2);
        assert!(p.iter().all(|q| q.inners().is_empty()));
    }

    #[test]
    fn overlapping_rings_are_rejected() {
        let err = from_rings(vec![sq(0.0, 4.0), sq(2.0, 6.0)]).unwrap_err();
        assert!(matches!(err, SketchError::RingsMeet { .. }), "{err:?}");
    }

    /// Touching counts as meeting: a hole flush against its outer ring has no strict inside.
    #[test]
    fn touching_rings_are_rejected() {
        let err = from_rings(vec![sq(0.0, 4.0), sq(4.0, 8.0)]).unwrap_err();
        assert!(matches!(err, SketchError::RingsMeet { .. }), "{err:?}");
    }

    #[test]
    fn a_two_point_ring_is_rejected() {
        let err = from_rings(vec![vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([1.0, 0.0]),
        ]])
        .unwrap_err();
        assert_eq!(err, SketchError::DegenerateRing { ring: 0 });
    }

    fn seg(a: [f64; 2], b: [f64; 2]) -> Edge2d {
        Edge2d::line(Point2::from_array(a), Point2::from_array(b))
    }

    /// Edges in scrambled order, some drawn backwards — the shape a generator emits.
    #[test]
    fn loose_edges_chain_into_a_ring() {
        let p = from_edges(vec![
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([0.0, 4.0], [0.0, 0.0]),
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([0.0, 4.0], [4.0, 4.0]), // drawn right-to-left
        ])
        .unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].outer().len(), 4);
        assert!(p[0].inners().is_empty());
    }

    /// Two cycles at once, one inside the other: chaining and nesting compose.
    #[test]
    fn two_chains_become_a_profile_with_a_hole() {
        let ring = |a: f64, b: f64| {
            vec![
                seg([a, a], [b, a]),
                seg([b, a], [b, b]),
                seg([b, b], [a, b]),
                seg([a, b], [a, a]),
            ]
        };
        let mut edges = ring(0.0, 4.0);
        edges.extend(ring(1.0, 3.0));
        let p = from_edges(edges).unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].inners().len(), 1);
    }

    /// A gap is not closed for the author — it is measured and handed back.
    #[test]
    fn an_unclosed_chain_reports_the_gap() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([4.0, 4.0], [0.0, 4.0]),
            seg([0.0, 4.0], [0.0, 0.25]), // 0.25 short of the start
        ])
        .unwrap_err();
        let SketchError::OpenChain { gap, .. } = err else {
            panic!("{err:?}");
        };
        assert!((gap.unwrap() - 0.25).abs() < 1e-12, "{gap:?}");
    }

    #[test]
    fn a_zero_length_edge_is_rejected() {
        let err = from_edges(vec![seg([1.0, 1.0], [1.0, 1.0])]).unwrap_err();
        assert_eq!(err, SketchError::ZeroLengthEdge { edge: 0 });
    }

    #[test]
    fn a_repeated_edge_is_rejected() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [0.0, 0.0]), // the same segment, reversed
            seg([4.0, 0.0], [4.0, 4.0]),
        ])
        .unwrap_err();
        assert!(matches!(err, SketchError::DuplicateEdge { .. }), "{err:?}");
    }

    /// A T-junction has no single continuation, so the walk refuses rather than picking one.
    #[test]
    fn a_branching_vertex_is_rejected() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([4.0, 0.0], [8.0, 0.0]),
        ])
        .unwrap_err();
        assert!(
            matches!(err, SketchError::BranchingVertex { .. }),
            "{err:?}"
        );
    }

    /// Containment must not depend on the rings' winding — the author draws in whatever direction
    /// is convenient, and `build_prism` fixes winding later anyway.
    #[test]
    fn winding_does_not_affect_nesting() {
        let mut hole = sq(1.0, 3.0);
        hole.reverse();
        let p = from_rings(vec![sq(0.0, 4.0), hole]).unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].inners().len(), 1);
    }
}
