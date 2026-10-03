//! Unit tests of the lateral face's chart geometry — the half-open rim-arc test the cover reads.

use super::*;
use nacre_exact::Rat;

fn q(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}

fn rational(v: [Rat; 3]) -> QuadVec {
    QuadVec {
        r0: v,
        r1: [Rat::from_int(0); 3],
        c: Rat::from_int(0),
    }
}

/// **The half-open test is the closed one less its `to` end** — over every ordered pair of a
/// ring of rational directions, so arcs shorter than, longer than and exactly a half turn are all
/// read, with `x` on both ends and between. The oracle is [`arc_contains`], the rational reading
/// the angular extent already uses, which owns the three-way split this one restates over
/// [`QuadVec`] ends.
#[test]
fn the_half_open_arc_is_the_closed_arc_less_its_end() {
    let z = Rat::from_int(0);
    let m = [z, z, Rat::from_int(1)];
    // Unit directions from Pythagorean triples, all the way round, with antipodal pairs.
    let ring: Vec<[Rat; 3]> = [
        (1, 1, 0, 1),
        (3, 5, 4, 5),
        (0, 1, 1, 1),
        (-4, 5, 3, 5),
        (-1, 1, 0, 1),
        (-3, 5, -4, 5),
        (0, 1, -1, 1),
        (4, 5, -3, 5),
        (5, 13, 12, 13),
        (-5, 13, -12, 13),
    ]
    .iter()
    .map(|&(a, b, c, d)| [q(a, b), q(c, d), z])
    .collect();
    let same_way = |u: &[Rat; 3], v: &[Rat; 3]| {
        nacre_exact::cross3_rat(u, v).unwrap() == [z; 3] && nacre_exact::dot3_rat(u, v).unwrap() > z
    };
    let mut half_turns = 0;
    let mut long_arcs = 0;
    for from in &ring {
        for to in &ring {
            if same_way(from, to) {
                continue; // `from == to` is a whole rim, which the cover reads before asking
            }
            let turn =
                nacre_exact::dot3_rat(&nacre_exact::cross3_rat(from, to).unwrap(), &m).unwrap();
            if turn == z {
                half_turns += 1;
            } else if turn < z {
                long_arcs += 1;
            }
            for x in &ring {
                let want = arc_contains(
                    &RimArc {
                        from: *from,
                        to: *to,
                    },
                    x,
                    &m,
                )
                .unwrap()
                    && !same_way(x, to);
                let got = arc_holds_half_open(&rational(*from), &rational(*to), x, &m).unwrap();
                assert_eq!(got, want, "from {from:?} to {to:?} x {x:?}");
            }
        }
    }
    // The population the branches need is there.
    assert!(half_turns > 0 && long_arcs > 0, "{half_turns} {long_arcs}");
}

/// **Two different radicals in one arc.** An end like `(√3/2, 1/2)` puts `√3` in the sign and the
/// other end's `√2` beside it — the second storey, [`nacre_exact::biquad_sign`]. Read against
/// `f64` angles at directions clear of the ends, each arc both ways round:
///
/// - `30° → 135°`, the plain case;
/// - `60° → 234.7…°` (`(1, √3)` to `(−1, −√2)`), a hair under a half turn: its turn is
///   `√3 − √2`, so a sign that took the radicals the wrong way round (`√2 − √3`) would read the
///   long way.
#[test]
fn an_arc_with_two_radicals_reads_like_its_angles() {
    let z = Rat::from_int(0);
    let m = [z, z, Rat::from_int(1)];
    let qv = |r0: [Rat; 3], r1: [Rat; 3], c: i128| QuadVec {
        r0,
        r1,
        c: Rat::from_int(c),
    };
    let arcs = [
        (
            qv([z, q(1, 2), z], [q(1, 2), z, z], 3),
            qv([z; 3], [q(-1, 2), q(1, 2), z], 2),
            30.0f64,
            135.0f64,
        ),
        (
            qv([q(1, 1), z, z], [z, q(1, 1), z], 3),
            qv([q(-1, 1), z, z], [z, q(-1, 1), z], 2),
            60.0,
            180.0 + 2f64.sqrt().atan().to_degrees(),
        ),
    ];
    let within = |a: f64, lo: f64, hi: f64| {
        let span = (hi - lo).rem_euclid(360.0);
        (a - lo).rem_euclid(360.0) < span
    };
    for (from, to, a_from, a_to) in &arcs {
        // `(5, 8)` and `(−5, −8)` (58°, 238°) sit in the two slivers where the half-turn reading
        // decides — between `54.7…°` and `60°`, and between `234.7…°` and `240°`.
        for (x, y) in [
            (1, 0),
            (3, 4),
            (5, 8),
            (0, 1),
            (-4, 3),
            (-1, 0),
            (-5, -8),
            (-3, -4),
            (0, -1),
            (4, -3),
        ] {
            let d = ((x * x + y * y) as f64).sqrt();
            let x_dir = [q(x, 1), q(y, 1), z];
            let angle = (y as f64 / d)
                .atan2(x as f64 / d)
                .to_degrees()
                .rem_euclid(360.0);
            let got = arc_holds_half_open(from, to, &x_dir, &m).unwrap();
            assert_eq!(
                got,
                within(angle, *a_from, *a_to),
                "{a_from}°→{a_to}° at {angle}°"
            );
            let back = arc_holds_half_open(to, from, &x_dir, &m).unwrap();
            assert_eq!(
                back,
                within(angle, *a_to, *a_from),
                "{a_to}°→{a_from}° at {angle}°"
            );
        }
    }
}
