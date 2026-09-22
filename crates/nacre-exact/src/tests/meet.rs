//! Unit tests of the three-plane solve's narrow road (`three_planes_rat_narrow`).

use super::*;

/// ★★★★★ **An overflowing intermediate does not cost the answer** —
/// locked on the fixture that documents the overflow.
///
/// These three planes meet at `(1, 1, 1)`, which fits `Rat` with room to spare. But their
/// coefficients carry coprime denominators — a power of two and a power of five, what decimal
/// arithmetic produces once it reduces — and the determinant is their product: `2⁹³·5³⁴`, a
/// ~172-bit denominator. The narrow Cramer multiplies before it can reduce, so it overflows;
/// this fixture asserts the fallback answers through the integer core — the hand-known point, so
/// this is an independent oracle
/// and not the two routes agreeing with each other.
///
/// ★ The private narrow route still declines here (asserted), so the fixture keeps proving
/// the fallback is *reached*, not merely present.
#[test]
fn an_overflowing_intermediate_no_longer_costs_the_answer() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let rows = [
        [r(1, 1 << 53), r(0, 1), r(0, 1), r(-1, 1 << 53)],
        [r(0, 1), r(1, 5i128.pow(23)), r(0, 1), r(-1, 5i128.pow(23))],
        [
            r(0, 1),
            r(0, 1),
            r(1, (1 << 40) * 5i128.pow(11)),
            r(-1, (1 << 40) * 5i128.pow(11)),
        ],
    ];
    assert_eq!(
        three_planes_rat_narrow(rows),
        None,
        "the narrow route was expected to overflow on coprime denominators — \
             without that this fixture no longer exercises the fallback"
    );
    assert_eq!(
        three_planes_rat(rows),
        Some([Rat::from_int(1); 3]),
        "the point fits `Rat`, and the solve now says so"
    );
    let names = rows.map(PlaneName::Narrow);
    let found = three_planes_big([&names[0], &names[1], &names[2]])
        .expect("the wide solve finds the point");
    assert_eq!(
        found,
        MeetPoint::Narrow([Rat::from_int(1); 3]),
        "the point fits `Rat` — the invariant demands it come back Narrow"
    );
    assert_eq!(found.width_bits(), 1, "1/1 is one bit wide");
}

/// ★ `three_planes_rat` — the exact corner of three rational planes. Hand-checkable
/// fixture, a degenerate (line-sharing) triple, and an i128-overflow decline.
#[test]
fn three_rational_planes_meet_where_they_should() {
    let r = Rat::from_int;
    // x = 2, y = 3, z = 5  (as a·x + d = 0 rows: [1,0,0,-2] etc.)
    let pt = three_planes_rat([
        [r(1), r(0), r(0), r(-2)],
        [r(0), r(1), r(0), r(-3)],
        [r(0), r(0), r(1), r(-5)],
    ])
    .expect("axis planes meet in one point");
    assert_eq!(pt, [r(2), r(3), r(5)]);
    // A tilted but exact triple: x+y=1, x−y=0, z=7 ⇒ (1/2, 1/2, 7).
    let pt = three_planes_rat([
        [r(1), r(1), r(0), r(-1)],
        [r(1), r(-1), r(0), r(0)],
        [r(0), r(0), r(1), r(-7)],
    ])
    .expect("a fair triple");
    assert_eq!(pt, [Rat::new(1, 2).unwrap(), Rat::new(1, 2).unwrap(), r(7)]);
    // Three planes through one line (z = 0, y = 0, y + z = 0): det = 0 — no unique point.
    assert_eq!(
        three_planes_rat([
            [r(0), r(0), r(1), r(0)],
            [r(0), r(1), r(0), r(0)],
            [r(0), r(1), r(1), r(0)],
        ]),
        None
    );
    // Coefficients near the i128 edge overflow the narrow Cramer — and the answer still
    // fits `Rat` (denominator 2¹²⁶ + 1), so the fallback answers it: declining here would be a defect.
    let big = Rat::from_int(1 << 126);
    let edge = [
        [big, r(1), r(0), r(-1)],
        [r(1), big, r(0), r(-1)],
        [r(0), r(0), big, big],
    ];
    assert_eq!(
        three_planes_rat_narrow(edge),
        None,
        "the narrow route overflows here"
    );
    let inv = Rat::new(1, (1 << 126) + 1).unwrap();
    assert_eq!(three_planes_rat(edge), Some([inv, inv, r(-1)]));
}
