use super::*;
use crate::predicate::tests::{RW, jd};

/// `r = 2` about the `z` axis, its seam at `+x`.
fn cylinder() -> CylinderDef {
    let r = |n| Rat::from_int(n);
    CylinderDef::new(
        [r(0), r(0), r(0)],
        [r(0), r(0), r(1)],
        [r(1), r(0), r(0)],
        nacre_exact::BigRat::from_rat(r(4)),
    )
    .expect("a cylinder")
}

/// The wall `x = 1`, the cap `z = 3`, the plane `y = 0`, the cap `z = 5`, the wall `x = 5` (beyond
/// the radius), the wall `x = 2` (tangent) and the wall `x = 1` again through other points, each
/// stated in the world or not.
fn planes(stated: bool) -> Vec<RW> {
    let rows = [
        [[1, 0, 0], [1, 1, 0], [1, 0, 1]],
        [[0, 0, 3], [1, 0, 3], [0, 1, 3]],
        [[0, 0, 0], [0, 0, 1], [1, 0, 0]],
        [[0, 0, 5], [1, 0, 5], [0, 1, 5]],
        [[5, 0, 0], [5, 1, 0], [5, 0, 1]],
        [[2, 0, 0], [2, 1, 0], [2, 0, 1]],
        [[1, 0, 0], [1, 2, 0], [1, 0, 2]],
    ]
    .map(|pts| RW::moved(pts, |p| p));
    rows.into_iter()
        .map(|w| if stated { w.stated() } else { w })
        .collect()
}

/// **The two roads put a pierce name's two points where each other does.** The line `x = 1, z = 3`
/// meets the cylinder at `y = ±√3`. Stated, the points are exact `a + b√c` pairs; with no world
/// statement they are interval boxes solved from the witness triangles. Each box holds one exact
/// point and the two boxes hold different ones — the roots' names may differ between the tables
/// (a stated plane orders them along its name, an unstated one along its stored facing), the points
/// may not.
#[test]
fn the_interval_road_holds_the_points_the_exact_road_names() {
    let def = cylinder();
    let (ps, qs) = (planes(true), planes(false));
    let (exact, judged) = (jd(&ps), jd(&qs));
    let prec = judged.standard.prec;
    let at = |j: &Judge<'_, RW>, root| {
        let p = j.pierce([0, 1], 0, &def, root).expect("the line crosses");
        let hp = p.point_hp(j).expect("a box");
        (p, hp)
    };
    let (e_lo, e_lo_hp) = at(&exact, QuadRoot::Lo);
    let (e_hi, e_hi_hp) = at(&exact, QuadRoot::Hi);
    assert!(e_lo.exact().is_some() && e_hi.exact().is_some());
    let (b_lo, b_lo_hp) = at(&judged, QuadRoot::Lo);
    let (_, b_hi_hp) = at(&judged, QuadRoot::Hi);
    assert!(
        b_lo.exact().is_none(),
        "no world statement — an interval box"
    );
    let overlaps = |a: &[HpBounded; 3], b: &[HpBounded; 3]| {
        (0..3).all(|k| a[k].sub(&b[k], prec).sign().is_none())
    };
    for b in [&b_lo_hp, &b_hi_hp] {
        assert_eq!(
            [&e_lo_hp, &e_hi_hp]
                .iter()
                .filter(|e| overlaps(e, b))
                .count(),
            1,
            "each box holds exactly one of the exact points"
        );
    }
    assert!(
        !overlaps(&b_lo_hp, &b_hi_hp),
        "the two boxes are two points"
    );
    // The plane `y = 0` parts the two points, on either road; a point is on its own planes.
    let sides = |j: &Judge<'_, RW>, p: &PiercePoint| (p.side(j, 2), p.side(j, 0));
    assert_eq!(e_lo.side(&exact, 2).map(i8::abs), Some(1));
    assert_eq!(e_lo.side(&exact, 2), e_hi.side(&exact, 2).map(|s| -s));
    assert_eq!(sides(&judged, &b_lo).1, Some(0));
    assert_eq!(b_lo.side(&judged, 2).map(i8::abs), Some(1));
}

/// **Which ruling of the wall a point is on** answers the same on both roads: the wall `x = 1`
/// cuts the cylinder in two rulings, and the two crossings of the cap line lie one on each.
#[test]
fn the_ruling_a_point_is_on_is_one_per_crossing_on_either_road() {
    let def = cylinder();
    for stated in [true, false] {
        let ps = planes(stated);
        let j = jd(&ps);
        let side = |root| {
            j.pierce([0, 1], 0, &def, root)
                .expect("the line crosses")
                .ruling_side(&j, 0, &def)
        };
        let (lo, hi) = (side(QuadRoot::Lo), side(QuadRoot::Hi));
        assert_eq!(lo.map(i8::abs), Some(1), "stated: {stated}");
        assert_eq!(lo, hi.map(|s| -s), "stated: {stated}");
    }
}

/// A stated coordinate where there is one, else the point.
fn value<'p>(p: &'p PiercePoint, q: &'p Option<QuadVal>) -> AxisValue<'p> {
    match q {
        Some(q) => AxisValue::Stated(q),
        None => AxisValue::Point(CylPoint::Pierce(p)),
    }
}

/// **The three orders agree with where the points are**, on both roads. On the wall `x = 1` the
/// cap `z = 3` crosses at `y = ±√3`: their world `y` order is the two roots' and opposite under a
/// swap; the three-plane corner `(1, 0, 3)` lies between them; along the line `x = 1, z = 3` the
/// order is the `y` order times the line's direction; and the cap `z = 5`'s crossing is higher up
/// the axis than the cap `z = 3`'s.
#[test]
fn the_orders_read_the_points_where_they_are() {
    let def = cylinder();
    for stated in [true, false] {
        let ps = planes(stated);
        let j = jd(&ps);
        let at = |pair, root| j.pierce(pair, 0, &def, root).expect("the line crosses");
        let (lo, hi) = (at([0, 1], QuadRoot::Lo), at([0, 1], QuadRoot::Hi));
        let y = j.world_order(CylPoint::Pierce(&lo), CylPoint::Pierce(&hi), 1);
        assert_eq!(y.map(i8::abs), Some(1), "stated: {stated}");
        assert_eq!(
            j.world_order(CylPoint::Pierce(&hi), CylPoint::Pierce(&lo), 1),
            y.map(|s| -s)
        );
        let corner = CylPoint::Three([0, 1, 2]);
        let c_lo = j.world_order(corner, CylPoint::Pierce(&lo), 1);
        let c_hi = j.world_order(corner, CylPoint::Pierce(&hi), 1);
        assert_eq!(c_lo, y.map(|s| -s), "the corner is between: {stated}");
        assert_eq!(c_hi, y, "the corner is between: {stated}");
        let along = j.line_order(0, 1, &lo, &hi);
        assert_eq!(along.map(i8::abs), Some(1), "stated: {stated}");
        assert_eq!(j.line_order(0, 1, &hi, &lo), along.map(|s| -s));
        let up = at([0, 3], QuadRoot::Lo);
        let (vl, vu) = (lo.stated_axis(&def), up.stated_axis(&def));
        assert_eq!(
            j.axis_order(&def, value(&lo, &vl), value(&up, &vu)),
            Some(core::cmp::Ordering::Less),
            "stated: {stated}"
        );
    }
}

/// **The meet's shape and the gate's two readings** on both roads: the cap line on the wall `x = 1`
/// crosses, on the wall `x = 5` it misses; a line along the axis inside or outside the cylinder is
/// off its surface, one on the tangent wall is not shown off; the wall `x = 1` cuts the lateral, the
/// wall `x = 5` does not, the tangent wall cannot be said.
#[test]
fn the_meet_and_the_gate_read_where_the_line_is() {
    let def = cylinder();
    for stated in [true, false] {
        let ps = planes(stated);
        let j = jd(&ps);
        assert_eq!(j.line_crosses(0, 1, &def), Some(true), "stated: {stated}");
        assert_eq!(j.line_crosses(4, 1, &def), Some(false), "stated: {stated}");
        assert!(j.parallel_line_off(0, 2, &def), "inside: {stated}");
        assert!(j.parallel_line_off(4, 2, &def), "outside: {stated}");
        assert!(!j.parallel_line_off(5, 2, &def), "on the surface: {stated}");
        assert_eq!(j.plane_cuts(0, &def), Some(true), "stated: {stated}");
        assert_eq!(j.plane_cuts(4, &def), Some(false), "stated: {stated}");
        assert_eq!(j.plane_cuts(5, &def), None, "tangent: {stated}");
    }
}

/// **Three points of one circle turn one way or the other**: the cap `z = 3` meets the wall
/// `x = 1` at `y = ±√3` and the plane `y = 0` at `x = ±2`; any three of the four are in one
/// counter-clockwise order or its reverse, and two names of distinct points are proved apart.
#[test]
fn three_points_of_a_circle_turn_one_way_on_either_road() {
    let def = cylinder();
    for stated in [true, false] {
        let ps = planes(stated);
        let j = jd(&ps);
        let at = |pair, root| j.pierce(pair, 0, &def, root).expect("the line crosses");
        let (a, b) = (at([0, 1], QuadRoot::Lo), at([0, 1], QuadRoot::Hi));
        let c = at([1, 2], QuadRoot::Lo);
        let (pa, pb, pc) = (
            CylPoint::Pierce(&a),
            CylPoint::Pierce(&b),
            CylPoint::Pierce(&c),
        );
        let turn = j.circle_turn(&def, pa, pc, pb);
        assert!(turn.is_some(), "stated: {stated}");
        assert_eq!(j.circle_turn(&def, pb, pc, pa), turn.map(|t| !t));
        assert!(j.points_apart(pa, pb));
        assert!(!j.points_apart(pa, pa));
    }
}

/// **A point with no exact meet is its definition, so it realizes deeper.** The unstated cap line's
/// crossing re-solved at twice the precision is a box inside the first one's reach — the two
/// overlap — and narrower.
#[test]
fn a_bounded_point_realizes_deeper_from_its_definition() {
    let def = cylinder();
    let ps = planes(false);
    let j = jd(&ps);
    let p = j
        .pierce([0, 1], 0, &def, QuadRoot::Lo)
        .expect("the line crosses");
    let prec = j.standard.prec;
    let (here, deep) = (
        p.point_hp_at(&j, prec).expect("a box"),
        p.point_hp_at(&j, 2 * prec).expect("a deeper box"),
    );
    for k in 0..3 {
        assert!(
            here[k].sub(&deep[k], 2 * prec).sign().is_none(),
            "the two boxes hold one point"
        );
        assert!(
            !here[k].error.lt(deep[k].error) || deep[k].error.is_zero(),
            "the deeper box is no wider"
        );
    }
}

/// **The line of two planes meets the cylinder in one vocabulary on both roads.** Stated, every
/// shape is named: two crossings, a miss, parallel caps, one plane twice, a line along the axis
/// inside the cylinder (off the surface) and one on it (a ruling), a tangency. Unstated, the interval
/// road names what it can — two crossings (realized in the order the name solves them), a miss, a
/// line along the axis the caller states, off the surface — and says it cannot say the rest.
#[test]
fn the_line_meets_the_cylinder_in_one_vocabulary() {
    let def = cylinder();
    let ps = planes(true);
    let j = jd(&ps);
    let shape = |p, q, along| j.line_meet(p, q, &def, along);
    assert!(matches!(
        shape(0, 1, false),
        LineMeet::Pair(Roots::Exact(..))
    ));
    assert!(matches!(shape(4, 1, false), LineMeet::Miss));
    assert!(matches!(shape(1, 3, false), LineMeet::ParallelPlanes));
    assert!(matches!(shape(0, 6, false), LineMeet::CoincidentPlanes));
    assert!(matches!(shape(0, 2, false), LineMeet::AxisParallelMiss));
    assert!(matches!(shape(5, 2, false), LineMeet::OnRuling(_)));
    assert!(matches!(shape(5, 1, false), LineMeet::Tangent(..)));
    let qs = planes(false);
    let k = jd(&qs);
    let shape = |p, q, along| k.line_meet(p, q, &def, along);
    assert!(matches!(
        shape(0, 1, false),
        LineMeet::Pair(Roots::Boxes { .. })
    ));
    assert!(matches!(shape(4, 1, false), LineMeet::Miss));
    assert!(matches!(shape(0, 2, true), LineMeet::AxisParallelMiss));
    assert!(matches!(shape(5, 2, true), LineMeet::Undecided));
    assert!(matches!(shape(5, 1, false), LineMeet::Undecided));
}

/// The wall `x = 1`, the cap `z = 3` and the plane `y = 0`, stated or not, each **facing either
/// way**: bit 0 turns the wall's stored coefficients against its triangle (`frame_sign = −1`, and a
/// stated name faces against them), bit 1 the cap's, bit 2 winds the plane `y = 0`'s triangle the
/// other way round (its outward flips with it).
fn faced(stated: bool, mask: u8) -> Vec<RW> {
    let wall = RW::moved([[1, 0, 0], [1, 1, 0], [1, 0, 1]], |p| p);
    let cap = RW::moved([[0, 0, 3], [1, 0, 3], [0, 1, 3]], |p| p);
    let plane = if mask & 4 == 0 {
        RW::moved([[0, 0, 0], [0, 0, 1], [1, 0, 0]], |p| p)
    } else {
        RW::moved([[0, 0, 0], [1, 0, 0], [0, 0, 1]], |p| p)
    };
    let wall = if mask & 1 == 0 { wall } else { wall.flipped() };
    let cap = if mask & 2 == 0 { cap } else { cap.flipped() };
    [wall, cap, plane]
        .into_iter()
        .map(|w| if stated { w.stated() } else { w })
        .collect()
}

fn sgn(x: f64) -> i8 {
    if x > 0.0 { 1 } else { -1 }
}

/// **Every sign answers the true geometry, whichever way the planes face.** On the wall `x = 1` the
/// cap `z = 3` crosses at `y = ±√3`; the crossing with `y > 0` (read off its box) lies on the
/// outward side of `y = 0` its triangle names, on the ruling of the wall the canonical normal
/// names, ahead of the other along the stored line direction, and `Lo` is the lower of the two
/// along the canonical one — on stated and unstated tables, with each plane's stored facing and
/// winding turned every way. Expectations come from the fixture's own facing (`coeffs`,
/// `tri_normal`, the names' canonical normals), never from the code under test.
#[test]
fn the_signs_read_the_true_geometry_whatever_the_facing() {
    let def = cylinder();
    for stated in [true, false] {
        for mask in 0..8u8 {
            let ps = faced(stated, mask);
            let j = jd(&ps);
            let at = |root| j.pierce([0, 1], 0, &def, root).expect("the line crosses");
            let (lo, hi) = (at(QuadRoot::Lo), at(QuadRoot::Hi));
            let y = |p: &PiercePoint| p.point_hp(&j).expect("a box")[1].sign().expect("off y = 0");
            let (up, down) = if y(&lo) { (&lo, &hi) } else { (&hi, &lo) };
            assert!(
                y(up) && !y(down),
                "one crossing on each side ({stated}, {mask})"
            );
            // The canonical normals: a stated plane's name (first nonzero component positive), an
            // unstated one's stored facing.
            let canon = |k: usize| {
                if stated {
                    [1.0, 1.0, 1.0][k]
                } else {
                    sgn(ps[k].coeffs()[[0, 2, 1][k]]) as f64
                }
            };
            let case = format!("stated {stated}, mask {mask}");
            assert_eq!(
                up.side(&j, 2),
                Some(sgn(ps[2].tri_normal()[1])),
                "outward side of y = 0 ({case})"
            );
            assert_eq!(up.side(&j, 0), Some(0), "on its own plane ({case})");
            assert_eq!(
                up.ruling_side(&j, 0, &def),
                Some(sgn(canon(0))),
                "the wall's ruling the canonical normal names ({case})"
            );
            let (c0, c1) = (ps[0].coeffs(), ps[1].coeffs());
            assert_eq!(
                j.line_order(0, 1, up, down),
                Some(sgn(-c0[0] * c1[2])),
                "ahead along the stored line direction ({case})"
            );
            assert_eq!(
                y(&hi),
                -canon(0) * canon(1) > 0.0,
                "`Hi` is the higher along the canonical line direction ({case})"
            );
            assert_eq!(
                j.world_order(CylPoint::Pierce(up), CylPoint::Pierce(down), 1),
                Some(1),
                "({case})"
            );
            // The meet's two roots name the same points in either plane order, `Lo` the lower
            // along that order's canonical direction.
            for pair in [[0, 1], [1, 0]] {
                let LineMeet::Pair(roots) = j.line_meet(pair[0], pair[1], &def, false) else {
                    panic!("two crossings ({case})");
                };
                let hi_up =
                    -canon(pair[0]) * canon(pair[1]) * if pair == [0, 1] { 1.0 } else { -1.0 };
                let p_hi = roots.point(&j, 0, &def, QuadRoot::Hi).expect("a point");
                assert_eq!(
                    y(&p_hi),
                    hi_up > 0.0,
                    "the meet's `Hi`, pair {pair:?} ({case})"
                );
            }
        }
    }
}

/// **Three points of one circle turn the way they stand**: about `+z`, from the crossing at
/// `y = −√3` through the plane `y = 0`'s crossing at `x = +2` to the one at `y = +√3` is
/// counter-clockwise, whichever way the planes face.
#[test]
fn a_circles_turn_reads_the_true_order_whatever_the_facing() {
    let def = cylinder();
    for stated in [true, false] {
        for mask in 0..8u8 {
            let ps = faced(stated, mask);
            let j = jd(&ps);
            let pick = |pair: [usize; 2], k: usize| {
                [QuadRoot::Lo, QuadRoot::Hi]
                    .into_iter()
                    .map(|r| j.pierce(pair, 0, &def, r).expect("the line crosses"))
                    .find(|p| p.point_hp(&j).expect("a box")[k].sign() == Some(true))
                    .expect("one crossing on the positive side")
            };
            let up = pick([0, 1], 1);
            let east = pick([1, 2], 0);
            let down = [QuadRoot::Lo, QuadRoot::Hi]
                .into_iter()
                .map(|r| j.pierce([0, 1], 0, &def, r).expect("the line crosses"))
                .find(|p| p.point_hp(&j).expect("a box")[1].sign() == Some(false))
                .expect("one crossing below");
            let (d, e, u) = (
                CylPoint::Pierce(&down),
                CylPoint::Pierce(&east),
                CylPoint::Pierce(&up),
            );
            assert_eq!(
                j.circle_turn(&def, d, e, u),
                Some(true),
                "stated {stated}, mask {mask}"
            );
            assert_eq!(
                j.circle_turn(&def, u, e, d),
                Some(false),
                "stated {stated}, mask {mask}"
            );
        }
    }
}
