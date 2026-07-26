//! Turning loose rings into profiles — the nesting policy on top of geom's exact predicates.
//!
//! A sketch is drawn as closed paths; which of them are material, which are holes, and which are
//! islands inside holes is **not** declared by the author but decided from containment, the way
//! every CAD sketcher does it. The classification itself (is this point inside that ring? do two
//! rings meet?) is a robustness-sensitive sign question and lives in `nacre-geom::intersect`
//! (design §1 isolates that); what is here is the policy — depth parity — and the assembly.

use crate::Profile2d;
use nacre_geom::intersect::{RingSide, point_in_ring_2d, rings_cross};
use nacre_math::Point2;

/// Why a set of rings is not a valid set of profiles.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SketchError {
    /// A ring with fewer than three points encloses nothing.
    DegenerateRing { ring: usize },
    /// Two rings touch or cross. A hole must lie strictly inside its outer ring and strictly
    /// outside its siblings; anything else has no unambiguous inside.
    RingsMeet { a: usize, b: usize },
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
