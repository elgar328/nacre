use super::*;

fn pt(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

fn sq(a: f64, b: f64) -> Vec<Point2> {
    vec![pt(a, a), pt(b, a), pt(b, b), pt(a, b)]
}

#[test]
fn a_lone_ring_is_one_profile() {
    let p = from_rings(vec![sq(0.0, 4.0)]).unwrap();
    assert_eq!(p.len(), 1);
    assert!(p[0].holes().is_empty());
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
        assert_eq!(p[0].holes().len(), 1, "one hole");
        assert_eq!(p[0].outer().vertices().len(), 4);
    }
}

/// Depth 2 is material again — the island is its own body, and it does *not* become a hole of
/// the outer ring.
#[test]
fn a_ring_inside_a_hole_is_an_island() {
    let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0)]).unwrap();
    assert_eq!(p.len(), 2, "outer body + island");
    let outer = p
        .iter()
        .find(|q| q.outer().vertices()[0][0] == Rat::from_int(0))
        .unwrap();
    let island = p
        .iter()
        .find(|q| q.outer().vertices()[0][0] == Rat::from_int(2))
        .unwrap();
    assert_eq!(outer.holes().len(), 1, "the depth-1 ring is its hole");
    assert!(island.holes().is_empty(), "nothing inside the island");
}

/// A hole belongs to the ring it is actually cut from, not to a distant ancestor.
#[test]
fn a_hole_attaches_to_its_immediate_container() {
    // outer(0..9) ⊃ hole(1..8) ⊃ island(2..7) ⊃ island-hole(3..6)
    let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0), sq(3.0, 6.0)]).unwrap();
    assert_eq!(p.len(), 2);
    let island = p
        .iter()
        .find(|q| q.outer().vertices()[0][0] == Rat::from_int(2))
        .unwrap();
    assert_eq!(island.holes().len(), 1, "the depth-3 ring is the island's");
    assert_eq!(island.holes()[0].vertices()[0][0], Rat::from_int(3));
}

/// Two separate bodies, neither inside the other.
#[test]
fn disjoint_rings_are_separate_bodies() {
    let p = from_rings(vec![sq(0.0, 1.0), sq(5.0, 6.0)]).unwrap();
    assert_eq!(p.len(), 2);
    assert!(p.iter().all(|q| q.holes().is_empty()));
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
    let err = from_rings(vec![vec![pt(0.0, 0.0), pt(1.0, 0.0)]]).unwrap_err();
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
    assert_eq!(p[0].holes().len(), 1);
}

// ---- the ring door ----

fn r(n: i128) -> Rat {
    Rat::from_int(n)
}

/// A ring stated as vertices and steps — the same square `from_rings` builds, bit for bit.
#[test]
fn a_stated_ring_is_the_polygon_door_s_ring() {
    let stated = Ring2d::new(
        vec![[r(0), r(0)], [r(4), r(0)], [r(4), r(4)], [r(0), r(4)]],
        vec![Edge2d::Line; 4],
    )
    .unwrap();
    assert_eq!(stated, Ring2d::polygon_decimal(sq(0.0, 4.0)).unwrap());
    let p = from_paths(vec![stated]).unwrap();
    assert_eq!(p, from_rings(vec![sq(0.0, 4.0)]).unwrap());
}

/// What a stated ring can get wrong is refused at the door, by name.
#[test]
fn a_stated_ring_is_checked_at_the_door() {
    let v = vec![[r(0), r(0)], [r(4), r(0)], [r(4), r(4)]];
    assert_eq!(
        Ring2d::new(v.clone(), vec![Edge2d::Line; 2]).unwrap_err(),
        SketchError::UnevenRing {
            vertices: 3,
            edges: 2
        }
    );
    assert_eq!(
        Ring2d::new(
            vec![[r(1), r(1)], [r(1), r(1)], [r(4), r(0)]],
            vec![Edge2d::Line; 3]
        )
        .unwrap_err(),
        SketchError::ZeroLengthEdge { edge: 0 }
    );
    // An arc step whose stated r² is not |start − centre|².
    let bad_radius = Edge2d::Arc {
        center: [r(0), r(0)],
        r2: r(3),
        ccw: true,
    };
    assert!(matches!(
        Ring2d::new(
            vec![[r(4), r(0)], [r(0), r(4)], [r(0), r(0)]],
            vec![bad_radius, Edge2d::Line, Edge2d::Line]
        ),
        Err(SketchError::ArcRadiusMismatch { .. })
    ));
    // An arc step whose next vertex is off its circle.
    let (quarter, _) = arc_turns(pt(0.0, 0.0), pt(4.0, 0.0), 1).unwrap();
    assert!(matches!(
        Ring2d::new(
            vec![[r(4), r(0)], [r(0), r(3)], [r(0), r(0)]],
            vec![quarter, Edge2d::Line, Edge2d::Line]
        ),
        Err(SketchError::ArcEndOffCircle { .. })
    ));
}

#[test]
fn a_circle_is_a_ring_of_one_vertex_at_its_seam() {
    let p = from_paths(vec![Ring2d::circle(pt(1.0, 2.0), 3.0).unwrap()]).unwrap();
    assert_eq!(p.len(), 1);
    let outer = p[0].outer();
    assert_eq!(outer.len(), 1);
    assert_eq!(
        outer.vertices()[0],
        [Rat::from_int(4), Rat::from_int(2)],
        "seam at +x"
    );
    assert!(matches!(outer.edges()[0], Edge2d::Arc { ccw: true, .. }));
    assert!(p[0].has_arcs());
    p[0].check().unwrap();
}

#[test]
fn concentric_circles_nest_like_squares() {
    let c = |r: f64| Ring2d::circle(pt(0.0, 0.0), r).unwrap();
    let p = from_paths(vec![c(4.0), c(20.0), c(12.0)]).unwrap();
    assert_eq!(p.len(), 2, "a ring and an island");
    let ring = p
        .iter()
        .find(|q| q.outer().vertices()[0][0] == Rat::from_int(20))
        .unwrap();
    assert_eq!(ring.holes().len(), 1);
    assert_eq!(ring.holes()[0].vertices()[0][0], Rat::from_int(12));
    let island = p
        .iter()
        .find(|q| q.outer().vertices()[0][0] == Rat::from_int(4))
        .unwrap();
    assert!(island.holes().is_empty());
}

/// A slot: two lines and two half circles, stated in order, the arcs' ends handed back by the
/// step door and used as the next vertices.
#[test]
fn a_slot_of_lines_and_arcs_is_one_ring() {
    // Centres (0,0) and (30,0), r 5. Start at (0,−5), east along the bottom.
    let (right, top_right) = arc_turns(pt(30.0, 0.0), pt(30.0, -5.0), 2).unwrap(); // → (30,5)
    let (left, bottom_left) = arc_turns(pt(0.0, 0.0), pt(0.0, 5.0), 2).unwrap(); // → (0,−5)
    assert_eq!(top_right, [r(30), r(5)]);
    assert_eq!(bottom_left, [r(0), r(-5)]);
    let ring = Ring2d::new(
        vec![bottom_left, [r(30), r(-5)], top_right, [r(0), r(5)]],
        vec![Edge2d::Line, right, Edge2d::Line, left],
    )
    .unwrap();
    let p = from_paths(vec![ring]).unwrap();
    assert_eq!(p.len(), 1);
    let o = p[0].outer();
    assert_eq!(o.len(), 4);
    assert_eq!(
        o.edges()
            .iter()
            .filter(|s| matches!(s, Edge2d::Arc { .. }))
            .count(),
        2
    );
    p[0].check().unwrap();
}

#[test]
fn two_half_circles_are_one_circle_in_normal_form() {
    let (a, mid) = arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 2).unwrap(); // (5,0) → (−5,0)
    let (b, back) = arc_turns(pt(0.0, 0.0), pt(-5.0, 0.0), 2).unwrap(); // (−5,0) → (5,0)
    assert_eq!(back, [r(5), r(0)]);
    let ring = Ring2d::new(vec![[r(5), r(0)], mid], vec![a, b]).unwrap();
    let p = from_paths(vec![ring]).unwrap();
    assert_eq!(p[0].outer().len(), 1, "merged into a whole circle");
    assert_eq!(
        p[0].outer().vertices()[0],
        [Rat::from_int(5), Rat::from_int(0)]
    );
    let c = from_paths(vec![Ring2d::circle(pt(0.0, 0.0), 5.0).unwrap()]).unwrap();
    assert_eq!(p, c, "the same profile `circle` states");
}

#[test]
fn equal_fillets_on_a_short_side_merge_into_a_half_circle() {
    // A 30 × 10 rectangle rounded r = 5 at every corner: the short sides vanish into half
    // circles — the slot, drawn as four lines and four quarter arcs.
    let q = |c: [f64; 2], s: [f64; 2]| arc_turns(pt(c[0], c[1]), pt(s[0], s[1]), 1).unwrap();
    let (q1, e1) = q([25.0, 5.0], [25.0, 0.0]); // (25,0) → (30,5)
    let (q2, e2) = q([25.0, 5.0], [30.0, 5.0]); // (30,5) → (25,10)
    let (q3, e3) = q([5.0, 5.0], [5.0, 10.0]); // (5,10) → (0,5)
    let (q4, e4) = q([5.0, 5.0], [0.0, 5.0]); // (0,5) → (5,0)
    assert_eq!((e1, e2), ([r(30), r(5)], [r(25), r(10)]));
    assert_eq!((e3, e4), ([r(0), r(5)], [r(5), r(0)]));
    let ring = Ring2d::new(
        vec![[r(5), r(0)], [r(25), r(0)], e1, e2, [r(5), r(10)], e3],
        vec![Edge2d::Line, q1, q2, Edge2d::Line, q3, q4],
    )
    .unwrap();
    let p = from_paths(vec![ring]).unwrap();
    assert_eq!(
        p[0].outer().len(),
        4,
        "line, half circle, line, half circle"
    );
    assert_eq!(
        p[0].outer()
            .edges()
            .iter()
            .filter(|s| matches!(s, Edge2d::Arc { .. }))
            .count(),
        2
    );
}

#[test]
fn arcs_that_cannot_be_stated_are_refused_by_name() {
    // r² = 2: no rational radius, and no refusal either — the truth is the square. The
    // quarter turn about (1, 1) from the origin lands on (2, 0).
    let (arc, end) = arc_turns(pt(1.0, 1.0), pt(0.0, 0.0), 1).expect("r² = 2 is a circle");
    assert_eq!(end, [r(2), r(0)]);
    assert!(matches!(arc, Edge2d::Arc { r2, .. } if r2 == r(2)));
    assert!(matches!(
        arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 0),
        Err(SketchError::ArcTurnsOutOfRange { turns: 0 })
    ));
    assert!(matches!(
        arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 4),
        Err(SketchError::ArcTurnsOutOfRange { turns: 4 })
    ));
    assert!(matches!(
        Ring2d::circle(pt(0.0, 0.0), 0.0),
        Err(SketchError::NonPositiveRadius { .. })
    ));
    assert!(matches!(
        arc_to_rat([r(0), r(0)], [r(5), r(0)], [r(0), r(4)], true),
        Err(SketchError::ArcEndOffCircle { .. })
    ));
    assert!(matches!(
        arc_to_rat([r(0), r(0)], [r(5), r(0)], [r(5), r(0)], true),
        Err(SketchError::ZeroLengthArc { .. })
    ));
    // 3-4-5: a rational radius off the axes is fine.
    assert!(arc_to_rat([r(0), r(0)], [r(3), r(4)], [r(-4), r(3)], true).is_ok());
}

#[test]
fn a_lens_of_two_arcs_is_a_valid_region() {
    // Two arcs between (0,0) and (6,0), centres (3,4) and (3,−4), r 5 — a lens.
    let up = arc_to_rat([r(3), r(4)], [r(0), r(0)], [r(6), r(0)], true).unwrap();
    let down = arc_to_rat([r(3), r(-4)], [r(6), r(0)], [r(0), r(0)], true).unwrap();
    let ring = Ring2d::new(vec![[r(0), r(0)], [r(6), r(0)]], vec![up, down]).unwrap();
    let p = from_paths(vec![ring]).unwrap();
    assert_eq!(p[0].outer().len(), 2);
    p[0].check().unwrap();
}

#[test]
fn hole_containment_reads_arcs_on_either_side() {
    let square = |a: f64, b: f64| Ring2d::polygon_decimal(sq(a, b)).unwrap();
    // A round hole in a square plate.
    let p = from_paths(vec![
        square(0.0, 10.0),
        Ring2d::circle(pt(5.0, 5.0), 2.0).unwrap(),
    ])
    .unwrap();
    assert_eq!((p.len(), p[0].holes().len()), (1, 1));
    assert!(!p[0].holes()[0].is_polygon());
    // A square hole in a round disk.
    let p = from_paths(vec![
        square(-2.0, 2.0),
        Ring2d::circle(pt(0.0, 0.0), 10.0).unwrap(),
    ])
    .unwrap();
    assert_eq!((p.len(), p[0].holes().len()), (1, 1));
    assert!(!p[0].outer().is_polygon());
    assert!(p[0].holes()[0].is_polygon());
    // A circle crossing the square: the rings meet.
    let rings = vec![
        square(0.0, 10.0),
        Ring2d::circle(pt(10.0, 5.0), 2.0).unwrap(),
    ];
    assert!(matches!(
        from_paths(rings),
        Err(SketchError::RingsMeet { .. })
    ));
}

/// A ring of two arcs that runs out along a circle and back is a spike of no area, and the sketch
/// doors name it — `from_paths` as a self-intersection, and an extrude handed it anyway as
/// `SelfIntersectingProfile`, before any solid is built.
#[test]
fn an_arc_and_its_retrace_are_a_self_intersecting_ring() {
    use crate::{OpError, Operation, SketchFrame, apply};
    let arc = |ccw| Edge2d::Arc {
        center: [r(0), r(0)],
        r2: r(25),
        ccw,
    };
    let ring = || {
        Ring2d::new(
            vec![[r(5), r(0)], [r(0), r(5)]],
            vec![arc(true), arc(false)],
        )
        .expect("each step is a stated arc")
    };
    assert!(matches!(
        from_paths(vec![ring()]),
        Err(SketchError::RingSelfIntersects { ring: 0, .. })
    ));
    let mut m = nacre_topo::Model::new();
    let frame = SketchFrame::world(&m, nacre_exact::Axis::Z);
    let profile = Profile2d::from_normalized_rings(ring(), vec![]);
    assert!(matches!(
        apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 1.0
            }
        ),
        Err(OpError::SelfIntersectingProfile { .. })
    ));
}

/// The normal form dissolves flat corners to a fixpoint and keeps everything `check` must still
/// see: four points on one line leave its two ends; a repeated point and a spike's tip survive;
/// a ring that is all one line comes back as its two ends, for `check` to call degenerate.
#[test]
fn the_normal_form_dissolves_flat_corners_and_keeps_every_reportable_defect() {
    let ring =
        |pts: &[[i128; 2]]| -> Vec<[Rat; 2]> { pts.iter().map(|p| [r(p[0]), r(p[1])]).collect() };
    let normal = |pts: &[[i128; 2]]| Ring2d::polygon(ring(pts)).vertices().to_vec();
    let square = ring(&[[0, 0], [3, 0], [3, 3], [0, 3]]);
    assert_eq!(
        normal(&[[0, 0], [1, 0], [2, 0], [3, 0], [3, 3], [0, 3]]),
        square
    );
    assert_eq!(
        normal(&[[0, 0], [3, 0], [3, 3], [0, 3]]),
        square,
        "a clean ring is untouched"
    );
    let dup = [[0, 0], [4, 0], [4, 0], [0, 4]];
    assert_eq!(normal(&dup), ring(&dup), "a repeated point survives");
    let spike = [[0, 0], [4, 0], [2, 0], [2, 4]];
    assert_eq!(normal(&spike), ring(&spike), "a spike's tip survives");
    assert_eq!(normal(&[[0, 0], [1, 0], [2, 0]]).len(), 2, "all one line");
}
