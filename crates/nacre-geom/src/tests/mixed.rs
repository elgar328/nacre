use super::*;
use crate::intersect::{point_in_ring_2d_rat, ring_self_intersection_rat, rings_cross_rat};
use proptest::prelude::*;

fn r(n: i128) -> Rat {
    Rat::from_int(n)
}
fn rq(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn pt(x: i128, y: i128) -> [Rat; 2] {
    [r(x), r(y)]
}
fn arc(cx: i128, cy: i128, radius: i128, ccw: bool) -> Edge2d {
    Edge2d::Arc {
        center: pt(cx, cy),
        r2: r(radius * radius),
        ccw,
    }
}

/// A slot: centres `(0,0)` and `(30,0)`, radius 5, counter-clockwise.
fn slot() -> (Vec<[Rat; 2]>, Vec<Edge2d>) {
    (
        vec![pt(0, -5), pt(30, -5), pt(30, 5), pt(0, 5)],
        vec![
            Edge2d::Line,
            arc(30, 0, 5, true),
            Edge2d::Line,
            arc(0, 0, 5, true),
        ],
    )
}

/// A whole circle of radius 5 about `(cx, cy)`, seam at `+x`.
fn circle(cx: i128, cy: i128, radius: i128, ccw: bool) -> (Vec<[Rat; 2]>, Vec<Edge2d>) {
    (vec![pt(cx + radius, cy)], vec![arc(cx, cy, radius, ccw)])
}

fn mk<'a>(v: &'a [[Rat; 2]], s: &'a [Edge2d]) -> MixedRing<'a> {
    MixedRing {
        vertices: v,
        segs: s,
    }
}

#[test]
fn a_whole_circle_classifies_points_by_distance() {
    let (v, s) = circle(0, 0, 5, true);
    let ring = mk(&v, &s);
    assert_eq!(point_in_mixed_ring(pt(0, 0), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(3, 3), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(6, 0), ring), Ok(RingSide::Outside));
    assert_eq!(point_in_mixed_ring(pt(-7, 1), ring), Ok(RingSide::Outside));
    assert_eq!(
        point_in_mixed_ring(pt(5, 0), ring),
        Ok(RingSide::OnBoundary)
    );
    assert_eq!(
        point_in_mixed_ring(pt(3, 4), ring),
        Ok(RingSide::OnBoundary)
    );
    // The probe level with the seam (the ray passes through the seam vertex) and level with
    // the top (a horizontal tangency) — both decided, neither a tie.
    assert_eq!(point_in_mixed_ring(pt(-1, 0), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(-9, 0), ring), Ok(RingSide::Outside));
    assert_eq!(point_in_mixed_ring(pt(-1, 5), ring), Ok(RingSide::Outside));
    assert_eq!(point_in_mixed_ring(pt(1, -5), ring), Ok(RingSide::Outside));
    // A clockwise circle is the same set of points.
    let (v2, s2) = circle(0, 0, 5, false);
    let cw = mk(&v2, &s2);
    assert_eq!(point_in_mixed_ring(pt(3, 3), cw), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(-1, 0), cw), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(6, 0), cw), Ok(RingSide::Outside));
    assert!(mixed_ring_self_intersection(cw).unwrap().is_none());
}

#[test]
fn a_slot_is_simple_and_classifies_its_points() {
    let (v, s) = slot();
    let ring = mk(&v, &s);
    assert_eq!(mixed_ring_self_intersection(ring), Ok(None));
    assert_eq!(point_in_mixed_ring(pt(15, 0), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(33, 0), ring), Ok(RingSide::Inside)); // in the round end
    assert_eq!(point_in_mixed_ring(pt(-3, 0), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(36, 0), ring), Ok(RingSide::Outside));
    assert_eq!(point_in_mixed_ring(pt(15, 6), ring), Ok(RingSide::Outside));
    assert_eq!(
        point_in_mixed_ring(pt(35, 0), ring),
        Ok(RingSide::OnBoundary)
    );
    assert_eq!(
        point_in_mixed_ring(pt(15, 5), ring),
        Ok(RingSide::OnBoundary)
    );
    assert_eq!(
        point_in_mixed_ring(pt(33, 4), ring),
        Ok(RingSide::OnBoundary)
    ); // 3-4-5
    // Probes level with the arc ends `y = ±5` (the ray runs along the straight sides'
    // line): the half-open rule decides them.
    assert_eq!(point_in_mixed_ring(pt(-10, 5), ring), Ok(RingSide::Outside));
    assert_eq!(
        point_in_mixed_ring(pt(-10, -5), ring),
        Ok(RingSide::Outside)
    );
    assert_eq!(point_in_mixed_ring(pt(-2, 0), ring), Ok(RingSide::Inside));
    // Level with the centres: the ray crosses the left arc once and the right arc once.
    assert_eq!(point_in_mixed_ring(pt(-6, 0), ring), Ok(RingSide::Outside));
}

#[test]
fn a_rounded_rectangle_is_simple_and_its_fillet_arcs_are_tangent_junctions() {
    // 40 × 20 with r = 5 at every corner: four lines, four quarter arcs, all ccw.
    let v = vec![
        pt(5, 0),
        pt(35, 0),
        pt(40, 5),
        pt(40, 15),
        pt(35, 20),
        pt(5, 20),
        pt(0, 15),
        pt(0, 5),
    ];
    let s = vec![
        Edge2d::Line,
        arc(35, 5, 5, true),
        Edge2d::Line,
        arc(35, 15, 5, true),
        Edge2d::Line,
        arc(5, 15, 5, true),
        Edge2d::Line,
        arc(5, 5, 5, true),
    ];
    let ring = mk(&v, &s);
    assert_eq!(mixed_ring_self_intersection(ring), Ok(None));
    assert_eq!(point_in_mixed_ring(pt(20, 10), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(1, 1), ring), Ok(RingSide::Outside)); // cut corner
    assert_eq!(point_in_mixed_ring(pt(2, 2), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(38, 2), ring), Ok(RingSide::Inside)); // (3,−3) from (35,5): 18 < 25
    assert_eq!(point_in_mixed_ring(pt(39, 1), ring), Ok(RingSide::Outside)); // 16 + 16 > 25
}

#[test]
fn rings_meet_when_arcs_cross_touch_or_coincide_and_not_otherwise() {
    let (c1v, c1s) = circle(0, 0, 5, true);
    let (c2v, c2s) = circle(8, 0, 5, true); // crossing
    let (c3v, c3s) = circle(10, 0, 5, true); // externally tangent at (5, 0)
    let (c4v, c4s) = circle(20, 0, 5, true); // apart
    let (c5v, c5s) = circle(0, 0, 3, true); // concentric inside
    let (c6v, c6s) = circle(0, 0, 5, false); // the same circle, drawn the other way
    let c1 = mk(&c1v, &c1s);
    assert_eq!(mixed_rings_cross(c1, mk(&c2v, &c2s)), Ok(true));
    assert_eq!(mixed_rings_cross(c1, mk(&c3v, &c3s)), Ok(true));
    assert_eq!(mixed_rings_cross(c1, mk(&c4v, &c4s)), Ok(false));
    assert_eq!(mixed_rings_cross(c1, mk(&c5v, &c5s)), Ok(false));
    assert_eq!(mixed_rings_cross(c1, mk(&c6v, &c6s)), Ok(true));
    // A square through a circle, a square around it, a square inside it.
    let sq = |a: i128, b: i128| -> (Vec<[Rat; 2]>, Vec<Edge2d>) {
        (
            vec![pt(a, a), pt(b, a), pt(b, b), pt(a, b)],
            vec![Edge2d::Line; 4],
        )
    };
    let (tv, ts) = sq(3, 9);
    let (av, as_) = sq(-9, 9);
    let (iv, is) = sq(-2, 2);
    assert_eq!(mixed_rings_cross(c1, mk(&tv, &ts)), Ok(true));
    assert_eq!(mixed_rings_cross(c1, mk(&av, &as_)), Ok(false));
    assert_eq!(mixed_rings_cross(c1, mk(&iv, &is)), Ok(false));
    // A square whose corner touches the circle at (3, 4): touching counts.
    let (kv, ks) = (
        vec![pt(3, 4), pt(9, 4), pt(9, 9), pt(3, 9)],
        vec![Edge2d::Line; 4],
    );
    assert_eq!(mixed_rings_cross(c1, mk(&kv, &ks)), Ok(true));
    // A line tangent to the circle at its top: touching counts too.
    let (lv, ls) = (
        vec![pt(-9, 5), pt(9, 5), pt(9, 9), pt(-9, 9)],
        vec![Edge2d::Line; 4],
    );
    assert_eq!(mixed_rings_cross(c1, mk(&lv, &ls)), Ok(true));
}

#[test]
fn a_ring_that_doubles_back_or_crosses_its_own_arc_is_named() {
    // Two adjacent arcs of one circle, the second running back over the first: the
    // clockwise arc (0,5) → (0,−5) about the origin passes through (5,0), the first arc's
    // other end — an overlap beyond the shared vertex.
    let v = vec![pt(5, 0), pt(0, 5), pt(0, -5)];
    let s = vec![arc(0, 0, 5, true), arc(0, 0, 5, false), Edge2d::Line];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 1))));
    // The slot with its two straight sides crossed (a bow tie with round ends).
    let v = vec![pt(0, -5), pt(30, 5), pt(30, -5), pt(0, 5)];
    let s = vec![
        Edge2d::Line,
        arc(30, 0, 5, false),
        Edge2d::Line,
        arc(0, 0, 5, true),
    ];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 2))));
    // A square whose right side is replaced by an arc bulging *inward* (clockwise about
    // (10,5), reaching x = 5): still simple.
    let v = vec![pt(0, 0), pt(10, 0), pt(10, 10), pt(0, 10)];
    let s = vec![
        Edge2d::Line,
        arc(10, 5, 5, false),
        Edge2d::Line,
        Edge2d::Line,
    ];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(None));
    // Pull the left side in to x = 7: it now cuts that bulge.
    let v = vec![pt(10, 0), pt(10, 10), pt(7, 10), pt(7, 0)];
    let s = vec![
        arc(10, 5, 5, false),
        Edge2d::Line,
        Edge2d::Line,
        Edge2d::Line,
    ];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 2))));
    // Two steps, the second running back along the first: one arc, out and back. Both ends are
    // the vertices the pair shares, so no end can witness the overlap — the arcs are one arc.
    let v = vec![pt(5, 0), pt(0, 5)];
    let s = vec![arc(0, 0, 5, true), arc(0, 0, 5, false)];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 1))));
    // The control: two half circles the same way round share both ends too, and are the
    // circle — they meet at their ends and nowhere else.
    let v = vec![pt(5, 0), pt(-5, 0)];
    let s = vec![arc(0, 0, 5, true), arc(0, 0, 5, true)];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(None));
}

#[test]
fn a_whole_circle_step_inside_a_longer_ring_and_a_zero_length_line_are_degenerate() {
    let v = vec![pt(5, 0), pt(5, 0)];
    let s = vec![arc(0, 0, 5, true), Edge2d::Line];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 0))));
    let v = vec![pt(0, 0), pt(0, 0), pt(1, 1)];
    let s = vec![Edge2d::Line; 3];
    assert_eq!(mixed_ring_self_intersection(mk(&v, &s)), Ok(Some((0, 0))));
}

#[test]
fn an_irrational_meeting_point_is_still_decided_exactly() {
    // The line y = 1 against the unit-ish circle r = 5 about the origin meets it at
    // x = ±√24 — irrational — and the probe (0, 1) is inside, (6, 1) outside.
    let (v, s) = circle(0, 0, 5, true);
    let ring = mk(&v, &s);
    assert_eq!(point_in_mixed_ring(pt(0, 1), ring), Ok(RingSide::Inside));
    assert_eq!(point_in_mixed_ring(pt(6, 1), ring), Ok(RingSide::Outside));
    assert_eq!(
        point_in_mixed_ring([rq(-489, 100), r(1)], ring),
        Ok(RingSide::Inside)
    ); // √24 ≈ 4.899
    assert_eq!(
        point_in_mixed_ring([rq(-491, 100), r(1)], ring),
        Ok(RingSide::Outside)
    );
    // Two circles crossing at irrational points: r = 5 about (0,0) and about (7,0) meet at
    // x = 3.5, y = ±√(25 − 12.25) — irrational.
    let (c2v, c2s) = circle(7, 0, 5, true);
    assert_eq!(
        mixed_rings_cross(
            ring,
            MixedRing {
                vertices: &c2v,
                segs: &c2s
            }
        ),
        Ok(true)
    );
    // Quarter arc (5,0)→(0,5) ccw and a line from (1,1) to (9,9): the crossing is at
    // (5/√2, 5/√2) — irrational — and it is on the arc.
    let av = vec![pt(5, 0), pt(0, 5), pt(0, 0)];
    let as_ = vec![arc(0, 0, 5, true), Edge2d::Line, Edge2d::Line];
    let lv = vec![pt(1, 1), pt(9, 9), pt(9, 1)];
    let ls = vec![Edge2d::Line; 3];
    assert_eq!(
        mixed_rings_cross(
            MixedRing {
                vertices: &av,
                segs: &as_
            },
            MixedRing {
                vertices: &lv,
                segs: &ls
            }
        ),
        Ok(true)
    );
    // The same line against the quarter arc (5,0)→(0,−5) clockwise — the other quadrant —
    // misses the arc although it crosses the circle.
    let bv = vec![pt(5, 0), pt(0, -5), pt(0, 0)];
    let bs = vec![arc(0, 0, 5, false), Edge2d::Line, Edge2d::Line];
    assert_eq!(
        mixed_rings_cross(
            MixedRing {
                vertices: &bv,
                segs: &bs
            },
            MixedRing {
                vertices: &lv,
                segs: &ls
            }
        ),
        Ok(false)
    );
}

fn poly_strategy() -> impl Strategy<Value = Vec<[Rat; 2]>> {
    prop::collection::vec((-6i128..=6, -6i128..=6), 3..=7)
        .prop_map(|v| v.into_iter().map(|(x, y)| pt(x, y)).collect())
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// On rings without arcs the mixed predicates are the polygon predicates, answer for
    /// answer — the arc road did not move the straight one.
    #[test]
    fn polygon_rings_answer_exactly_as_before(a in poly_strategy(), b in poly_strategy(), px in -8i128..=8, py in -8i128..=8) {
        let sa = vec![Edge2d::Line; a.len()];
        let sb = vec![Edge2d::Line; b.len()];
        let ra = MixedRing { vertices: &a, segs: &sa };
        let rb = MixedRing { vertices: &b, segs: &sb };
        prop_assert_eq!(mixed_ring_self_intersection(ra).unwrap(), ring_self_intersection_rat(&a));
        prop_assert_eq!(mixed_rings_cross(ra, rb).unwrap(), rings_cross_rat(&a, &b));
        if ring_self_intersection_rat(&a).is_none() {
            prop_assert_eq!(point_in_mixed_ring(pt(px, py), ra).unwrap(), point_in_ring_2d_rat(pt(px, py), &a));
        }
    }
}
