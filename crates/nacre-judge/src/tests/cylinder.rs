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

/// The wall `x = 1`, the cap `z = 3` and the plane `y = 0`, each stated in the world or not.
fn planes(stated: bool) -> Vec<RW> {
    let rows = [
        [[1, 0, 0], [1, 1, 0], [1, 0, 1]],
        [[0, 0, 3], [1, 0, 3], [0, 1, 3]],
        [[0, 0, 0], [0, 0, 1], [1, 0, 0]],
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
