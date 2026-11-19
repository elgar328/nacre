//! **A plane face of a result meets its own boundary only where one loop passes a vertex more than
//! once** ([`crate::assembly::face_bounds_meet`]) — on hand-built faces of the plane `z = 0`, every
//! corner the meeting of that plane with one plane `x = i` and one `y = j`, so each meet the
//! arrangement could have failed to cut is laid out by hand. A refusal is `FaceBoundsMeet`; the
//! faces that pass are the ones the arrangement does build. The arcs' cases run end to end:
//! `tangent_lumps.rs`'s slab between two slanted walls, whose touch points are vertices.

use super::coplanar::mk_axis_plane;
use super::*;

/// The class `z = 0` (index `0`), the planes `x = 0..=4` (`1..=5`) and `y = 0..=4` (`6..=10`), two
/// slanted ones, `x + y = 2` (`11`) and `x − y = 2` (`12`), and the three lines of the triangle
/// `(1, 1), (7/3, 5/3), (5/3, 7/3)`: through `(1, 1)` and `(3, 2)` (`13`), through `(1, 1)` and
/// `(2, 3)` (`14`), and `x + y = 4` (`15`).
fn grid() -> Vec<WorkingPlane> {
    let mut m = Model::new();
    let mut p = vec![mk_axis_plane(&mut m, 2, 0.0, true)];
    for axis in [0, 1] {
        for d in 0..=4 {
            p.push(mk_axis_plane(&mut m, axis, d as f64, true));
        }
    }
    let pt = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    p.push(plane_through(
        &mut m,
        [pt(2.0, 0.0, 0.0), pt(0.0, 2.0, 0.0), pt(2.0, 0.0, 1.0)],
    ));
    p.push(plane_through(
        &mut m,
        [pt(2.0, 0.0, 0.0), pt(3.0, 1.0, 0.0), pt(2.0, 0.0, 1.0)],
    ));
    for [a, b] in [
        [(1.0, 1.0), (3.0, 2.0)],
        [(1.0, 1.0), (2.0, 3.0)],
        [(4.0, 0.0), (0.0, 4.0)],
    ] {
        p.push(plane_through(
            &mut m,
            [pt(a.0, a.1, 0.0), pt(b.0, b.1, 0.0), pt(a.0, a.1, 1.0)],
        ));
    }
    p
}

/// The corner `(x, y)` of the grid.
fn at(x: usize, y: usize) -> NodeId {
    NodeId::three_planes(Canon3::three([0, 1 + x, 6 + y]))
}

fn corners(xy: &[(usize, usize)]) -> Vec<NodeId> {
    xy.iter().map(|&(x, y)| at(x, y)).collect()
}

/// What the arm answers of one face: `None` for a face it passes, the reason otherwise.
fn asked(outer: &[(usize, usize)], holes: &[&[(usize, usize)]]) -> Option<RejectReason> {
    ask(corners(outer), holes.iter().map(|h| corners(h)).collect())
}

fn ask(outer: Vec<NodeId>, holes: Vec<Vec<NodeId>>) -> Option<RejectReason> {
    let p = grid();
    let lf = face(0, outer, holes);
    match crate::assembly::face_bounds_meet(
        &crate::planes::test_judge(&p),
        &[],
        &Default::default(),
        &[lf],
    ) {
        Ok(()) => None,
        Err(BoolError::Rejected { reason, .. }) => Some(reason),
        Err(e) => panic!("not a rejection: {e:?}"),
    }
}

const SQUARE: &[(usize, usize)] = &[(0, 0), (4, 0), (4, 4), (0, 4)];

/// Each way a face's boundary can meet itself away from a vertex both pieces carry.
#[test]
fn a_face_whose_boundary_meets_itself_elsewhere_is_refused() {
    // Each case answers on its own, so switching one check off shows every case it alone holds.
    let mut rows: Vec<(&str, Option<RejectReason>)> = Vec::new();
    // A hole crossing the outer ring's edge `x = 4` between two of its vertices.
    let crossing: &[(usize, usize)] = &[(1, 1), (1, 2), (4, 2), (4, 1)];
    rows.push((
        "a crossing",
        asked(&[(0, 0), (3, 0), (3, 3), (0, 3)], &[crossing]),
    ));
    // A hole whose corner sits inside the outer ring's edge — and runs along it.
    let along: &[(usize, usize)] = &[(3, 1), (3, 2), (4, 2), (4, 1)];
    rows.push(("a vertex inside an edge", asked(SQUARE, &[along])));
    // A triangular hole whose corner `(2, 0)` rests on the outer ring's edge `y = 0`, its two edges
    // leaving it at an angle — a vertex inside an edge with nothing along it.
    let (slant, back, y1) = (11, 12, 7);
    let three = |a: usize, b: usize| NodeId::three_planes(Canon3::three([0, a, b]));
    let triangle = vec![three(slant, back), three(slant, y1), three(back, y1)];
    rows.push((
        "a vertex on an edge at an angle",
        ask(corners(SQUARE), vec![triangle]),
    ));
    // Two holes, each its own loop, through one vertex.
    let h1: &[(usize, usize)] = &[(1, 1), (1, 2), (2, 2), (2, 1)];
    let h2: &[(usize, usize)] = &[(2, 2), (2, 3), (3, 3), (3, 2)];
    rows.push(("two loops through one vertex", asked(SQUARE, &[h1, h2])));
    // One outer loop through `(1, 1)` twice, both parts wound outward: two lobes meeting there.
    let lobes = [
        (0, 0),
        (1, 0),
        (1, 1),
        (2, 1),
        (2, 2),
        (1, 2),
        (1, 1),
        (0, 1),
    ];
    rows.push(("two lobes", asked(&lobes, &[])));
    // One hole loop through `(2, 2)` twice, one part wound as a hole and one as material.
    let island = [
        (1, 1),
        (1, 2),
        (2, 2),
        (3, 2),
        (3, 3),
        (2, 3),
        (2, 2),
        (2, 1),
    ];
    rows.push(("a hole holding a lobe", asked(SQUARE, &[&island])));
    // An outer loop running out along `y = 1` and back over the same edge.
    let spike = [(0, 0), (2, 0), (2, 1), (3, 1), (2, 1), (2, 2), (0, 2)];
    rows.push(("an edge passed twice", asked(&spike, &[])));
    // An outer loop through `(1, 1)` twice, once up `x = 1` and once across `y = 1`: it crosses
    // itself there. Its parts wind `−1` and `+1`, as a notch's do.
    let crossed = [
        (0, 0),
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 2),
        (2, 1),
        (1, 1),
        (0, 1),
    ];
    rows.push(("a loop crossing itself at a vertex", asked(&crossed, &[])));
    // `(2, 0)` twice, under two names: a corner of the outer ring and the hole triangle's corner
    // `x + y = 2 ∩ x − y = 2` — one point, a vertex of neither's edge under the other's name.
    let named_twice = corners(&[(0, 0), (2, 0), (4, 0), (4, 4), (0, 4)]);
    let triangle = vec![three(slant, back), three(slant, y1), three(back, y1)];
    rows.push((
        "an end standing on another end under another name",
        ask(named_twice, vec![triangle]),
    ));
    let wrong: Vec<_> = rows
        .iter()
        .filter(|(_, got)| *got != Some(RejectReason::FaceBoundsMeet))
        .collect();
    assert!(wrong.is_empty(), "cases answered otherwise: {wrong:?}");
}

/// **A hole loop whose two parts nest and cross at the vertex they share** — the square `(1, 1),
/// (1, 3), (3, 3), (3, 1)` and the triangle `(1, 1), (5/3, 7/3), (7/3, 5/3)` strictly inside it,
/// both wound as holes, one loop through `(1, 1)` twice: it arrives from `(3, 1)` and leaves for
/// `(5/3, 7/3)`, and arrives from `(7/3, 5/3)` and leaves for `(1, 3)`, so the two passes cross
/// there. Its parts wind
/// `−1` and `−1`, as two holes touching do.
#[test]
fn a_hole_loop_crossing_itself_at_a_vertex_is_refused() {
    use crate::combinatorics::Wall;
    let three = |a: usize, b: usize| NodeId::three_planes(Canon3::three([0, a, b]));
    let outer = vec![
        (at(0, 0), Wall::Plane(6)),
        (at(4, 0), Wall::Plane(5)),
        (at(4, 4), Wall::Plane(10)),
        (at(0, 4), Wall::Plane(1)),
    ];
    let hole = vec![
        (at(1, 1), Wall::Plane(2)),
        (at(1, 3), Wall::Plane(9)),
        (at(3, 3), Wall::Plane(4)),
        (at(3, 1), Wall::Plane(7)),
        (at(1, 1), Wall::Plane(14)),
        (three(14, 15), Wall::Plane(15)),
        (three(13, 15), Wall::Plane(13)),
    ];
    assert_eq!(
        ask_rings(outer, vec![hole]),
        Some(RejectReason::FaceBoundsMeet)
    );
}

/// The faces the arrangement builds pass: a plain ring with a hole, and two holes touching at a
/// vertex, written as one loop through it twice — both parts wound as holes.
#[test]
fn a_face_pinched_into_one_region_passes() {
    let hole: &[(usize, usize)] = &[(1, 1), (1, 3), (3, 3), (3, 1)];
    assert_eq!(asked(SQUARE, &[hole]), None, "a ring and its hole");
    let touching = [
        (1, 1),
        (1, 2),
        (2, 2),
        (2, 3),
        (3, 3),
        (3, 2),
        (2, 2),
        (2, 1),
    ];
    assert_eq!(
        asked(SQUARE, &[&touching]),
        None,
        "two holes touching at a vertex"
    );
}

/// A face of hand-laid rings: `outer` and `holes` as `(node, wall)` runs, each wall the one the
/// edge from that node rides.
fn ask_rings(
    outer: Vec<(NodeId, crate::combinatorics::Wall)>,
    holes: Vec<Vec<(NodeId, crate::combinatorics::Wall)>>,
) -> Option<RejectReason> {
    use crate::draft::{Bound, Ring};
    let ring = |r: Vec<(NodeId, crate::combinatorics::Wall)>| {
        let (nodes, walls) = r.into_iter().unzip();
        Bound::Ring(Ring::new(nodes, walls))
    };
    let p = grid();
    let lf = LocalFace {
        surf: crate::planes::ClassIx::Plane(0),
        outer: ring(outer),
        inner: holes.into_iter().map(ring).collect(),
        flip: false,
    };
    match crate::assembly::face_bounds_meet(
        &crate::planes::test_judge(&p),
        &[],
        &Default::default(),
        &[lf],
    ) {
        Ok(()) => None,
        Err(BoolError::Rejected { reason, .. }) => Some(reason),
        Err(e) => panic!("not a rejection: {e:?}"),
    }
}

/// **A pinch whose ring cannot be read is no proof of a meet.** The outer ring is pinched at
/// `(2, 2)` into the face's outline and a notch wound as a hole — one region — and the edge arriving
/// at the pinch from `(4, 2)` is an arc of a cylinder the table does not hold, so the ring's edges,
/// and with them which way each visit arrives and leaves, are unread. The verdict is the held
/// failure, not a meet.
#[test]
fn an_unread_pinch_is_held_not_a_meet() {
    use crate::combinatorics::Wall;
    let (y0, y1, y2) = (6, 7, 8);
    let (x0, x1, x2, x4) = (1, 2, 3, 5);
    let y4 = 10;
    let outer = vec![
        (at(0, 0), Wall::Plane(y0)),
        (at(4, 0), Wall::Plane(x4)),
        (
            at(4, 2),
            Wall::Arc {
                cyl: 0,
                ccw: true,
                plane: 0,
            },
        ),
        (at(2, 2), Wall::Plane(x2)),
        (at(2, 1), Wall::Plane(y1)),
        (at(1, 1), Wall::Plane(x1)),
        (at(1, 2), Wall::Plane(y2)),
        (at(2, 2), Wall::Plane(x2)),
        (at(2, 4), Wall::Plane(y4)),
        (at(0, 4), Wall::Plane(x0)),
    ];
    assert_eq!(
        ask_rings(outer, Vec::new()),
        Some(RejectReason::WitnessNotRational)
    );
}

/// **A pinch whose way out cannot be read is no meet.** The notch of the test above, all straight:
/// the pinch is named as where the line `x = 2` pierces a cylinder the table does not hold. The
/// ring's edges are built from the names and read; the two arrivals along `y = 2` are pinned by
/// the plane `x = 2`, but the two departures along `x = 2` itself are pinned by the cylinder, so
/// their directions — which way each visit leaves — are unread. The verdict is the held failure.
/// Counted as a pass instead, the pinch would still be held here: the node's other readings need
/// the same cylinder.
#[test]
fn a_pinch_whose_direction_is_unread_is_not_a_meet() {
    use crate::combinatorics::{QuadRoot, Wall};
    let (y0, y1, y2) = (6, 7, 8);
    let (x0, x1, x2, x4) = (1, 2, 3, 5);
    let y4 = 10;
    let pinch = NodeId::pierce(x2, 0, 0, QuadRoot::Lo);
    let outer = vec![
        (at(0, 0), Wall::Plane(y0)),
        (at(4, 0), Wall::Plane(x4)),
        (at(4, 2), Wall::Plane(y2)),
        (pinch, Wall::Plane(x2)),
        (at(2, 1), Wall::Plane(y1)),
        (at(1, 1), Wall::Plane(x1)),
        (at(1, 2), Wall::Plane(y2)),
        (pinch, Wall::Plane(x2)),
        (at(2, 4), Wall::Plane(y4)),
        (at(0, 4), Wall::Plane(x0)),
    ];
    // Premise: the edges read, and a departure from the pinch does not — the reading part 1 holds.
    let p = grid();
    let jd = crate::planes::test_judge(&p);
    let (nodes, walls): (Vec<NodeId>, Vec<Wall>) = outer.iter().cloned().unzip();
    let edges = crate::draft::Ring::new(nodes.clone(), walls)
        .edges(&jd, &[], Some(0))
        .expect("premise: the ring's edges are read from its names");
    let departures: Vec<bool> = (0..nodes.len())
        .filter(|&i| nodes[i] == pinch)
        .map(|i| crate::combinatorics::dir_at(&jd, &[], 0, &edges[i], pinch).is_ok())
        .collect();
    assert_eq!(
        departures,
        [false, false],
        "premise: no departure from the pinch is read"
    );
    assert_eq!(
        ask_rings(outer, Vec::new()),
        Some(RejectReason::WitnessNotRational)
    );
}

/// **A straight edge and a ruling with the same two ends are one point set** — both pieces are
/// straight, so the loop goes out along a line and back over it. Asked by the ends alone: the
/// ruling's cylinder is not in the table, so nothing else about it can be read.
#[test]
fn a_segment_and_a_ruling_with_the_same_ends_meet() {
    use crate::combinatorics::Wall;
    let outer = vec![
        (at(0, 0), Wall::Plane(6)),
        (
            at(4, 0),
            Wall::Ruling {
                cyl: 0,
                side: 1,
                up: true,
                plane: 0,
            },
        ),
    ];
    assert_eq!(
        ask_rings(outer, Vec::new()),
        Some(RejectReason::FaceBoundsMeet)
    );
}
