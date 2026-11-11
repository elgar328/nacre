//! Unit tests of the lateral face's chart geometry — the half-open rim-arc test the cover reads,
//! and the band a closed trace sweeps.

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
                let got = arc_holds_beside(&rational(*from), &rational(*to), x, &m, true).unwrap();
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
            let got = arc_holds_beside(from, to, &x_dir, &m, true).unwrap();
            assert_eq!(
                got,
                within(angle, *a_from, *a_to),
                "{a_from}°→{a_to}° at {angle}°"
            );
            let back = arc_holds_beside(to, from, &x_dir, &m, true).unwrap();
            assert_eq!(
                back,
                within(angle, *a_to, *a_from),
                "{a_to}°→{a_from}° at {angle}°"
            );
        }
    }
}

/// **The band a closed trace sweeps along its own axis.** The cylinder stands on the `z` axis,
/// radius 2. On `z = 5` the trace is the circle at station 5. On `z = y/4 + 5/2` it is an ellipse
/// over `z ∈ [2, 3]`, `t ± Δ` with `Δ = 2·(1/4)/1`. With the raw axis `(0, 0, 2)` the stations
/// halve and so does `Δ` — `|n·m|` doubles. With radius `√2` and `z = y + 5` the band is
/// `5 ± √2`, rounded outward to rationals. A plane along the axis sweeps no band.
#[test]
fn a_traces_band_is_its_axis_meet_widened_by_the_ellipses_reach() {
    let z = |v: i128| Rat::from_int(v);
    let cylinder = |dir: [Rat; 3], r2: Rat| {
        nacre_topo::CylinderDef::new(
            [z(0); 3],
            dir,
            [z(1), z(0), z(0)],
            nacre_exact::BigRat::from(r2),
        )
        .unwrap()
    };
    let up = cylinder([z(0), z(0), z(1)], z(4));
    assert_eq!(
        trace_band(&[z(0), z(0), z(1), z(-5)], &up),
        Some([z(5), z(5)])
    );
    let slant = [z(0), q(-1, 4), z(1), q(-5, 2)];
    assert_eq!(trace_band(&slant, &up), Some([z(2), z(3)]));
    let long = cylinder([z(0), z(0), z(2)], z(4));
    assert_eq!(trace_band(&slant, &long), Some([z(1), q(3, 2)]));
    let thin = cylinder([z(0), z(0), z(1)], z(2));
    let [lo, hi] = trace_band(&[z(0), z(-1), z(1), z(-5)], &thin).unwrap();
    let (below, above) = (z(5).checked_sub(lo).unwrap(), hi.checked_sub(z(5)).unwrap());
    for d in [below, above] {
        assert!(d.checked_mul(d).unwrap() >= z(2), "never narrower than √2");
        assert!(d < q(1415, 1000), "and near it");
    }
    assert_eq!(
        trace_band(&[z(1), z(0), z(0), z(-1)], &up),
        None,
        "along the axis"
    );
}

/// `body`, cut by the slant `z = 5 + x/2` (the obliquely cut rod's prism), and its one lateral face
/// with the cylinder's world statement.
fn slanted_lateral(
    body: fn(&mut Model) -> Handle<nacre_topo::Solid>,
) -> (Model, Handle<nacre_topo::Face>, nacre_topo::CylinderDef) {
    let mut m = Model::new();
    let b = body(&mut m);
    let p = crate::tests::slant_prism(&mut m, 0.5);
    let out = crate::boolean(&mut m, crate::BoolKind::Cut, b, p).expect("the slant cuts");
    let laterals: Vec<_> = m
        .reachable_from(&out)
        .faces
        .into_iter()
        .filter(|&f| {
            matches!(
                m.surface(m.face(f).surface),
                nacre_topo::Surface::Cylinder { .. }
            )
        })
        .collect();
    assert_eq!(laterals.len(), 1, "one lateral face");
    let def = m
        .world_cylinder_def(m.face(laterals[0]).surface)
        .expect("a world statement");
    (m, laterals[0], def)
}

/// **A lateral bounded by a slanted rim reaches the rim's band, and has no heights** — the rod
/// (radius 2 about `(0, 3)`, `z ∈ [0, 10]`) cut along `z = 5 + x/2`: its rim sweeps `[4, 6]`, so
/// the face reaches `[0, 6]` (`lateral_axial_reach`), while its range by rim heights is `None` —
/// an ellipse has no one height (`lateral_t_range`). A window cut in the side first (`x ≥ 1`,
/// `z ∈ [2, 3]`) adds two circles at `z = 2, 3`; read past the slanted rim they would make a range
/// `[0, 3]` that stops short of the face. Both rims are whole, so the angle is the whole turn.
#[test]
fn a_lateral_with_a_slanted_rim_reaches_its_band_and_has_no_heights() {
    let notched = |m: &mut Model| {
        let r = crate::tests::slant_rod(m);
        let b = crate::fixtures::cuboid(m, [1.0, 2.0, 2.0], [3.0, 4.0, 3.0]);
        let out = crate::boolean(m, crate::BoolKind::Cut, r, b).expect("the window");
        out[0]
    };
    for body in [crate::tests::slant_rod as fn(&mut Model) -> _, notched] {
        let (m, f, def) = slanted_lateral(body);
        let face = m.face(f);
        assert_eq!(
            lateral_axial_reach(&m, face, &def),
            Some([Rat::from_int(0), Rat::from_int(6)])
        );
        assert_eq!(lateral_t_range(&m, face, &def), None);
        assert_eq!(lateral_theta_extent(&m, face, &def), Theta::Whole);
    }
}

/// **A half rod cut by the slant spans its half turn** — the rod less `y < 3`: its lateral's rims
/// are a half circle below and a half ellipse above, both from the generator through `(2, 3)` to
/// the one through `(−2, 3)`, so the angular extent is the arc from `+x̂` to `−x̂` counter-clockwise
/// about `+z`. The slanted rim's corners stand at `z = 6` and `z = 4`, off its axis meet `z = 5`,
/// so each is read by its direction off the axis (`rim_radial`), not by its offset from that meet.
#[test]
fn a_half_rod_cut_by_the_slant_spans_its_half_turn() {
    let half = |m: &mut Model| {
        let r = crate::tests::slant_rod(m);
        let b = crate::fixtures::cuboid(m, [-3.0, -1.0, -1.0], [3.0, 3.0, 11.0]);
        let out = crate::boolean(m, crate::BoolKind::Cut, r, b).expect("the half");
        out[0]
    };
    let (m, f, def) = slanted_lateral(half);
    let z = Rat::from_int(0);
    assert_eq!(
        lateral_theta_extent(&m, m.face(f), &def),
        Theta::Arc(RimArc {
            from: [Rat::from_int(2), z, z],
            to: [Rat::from_int(-2), z, z],
        })
    );
}
