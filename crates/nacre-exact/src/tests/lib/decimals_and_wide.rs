//! Decimal realization and the wide (arbitrary-precision) twins of the narrow solves.

use super::*;

/// **The property `from_decimal` rests on**: the shortest decimal of an f64,
/// read as an exact rational and realized again, gives back the same bits.
///
/// It holds *because* `to_f64` is correctly rounded, and only because of that. The
/// decimal is by construction a value whose nearest f64 is `x`; nearest rounding
/// therefore has no choice. Rounding numerator and denominator separately first —
/// what this function used to do — fails **17.4%** of the values below (measured),
/// since a 17-digit decimal has a numerator past 2⁵³. That is the whole reason this
/// cell touches `to_f64` at all.
#[test]
fn a_shortest_decimal_realizes_back_to_its_own_f64() {
    let mut state = 0x243f_6a88_85a3_08d3u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let (mut tried, mut declined) = (0, 0);
    for _ in 0..200_000 {
        // A random sign and mantissa over a CAD-plausible exponent range.
        let bits = next();
        let exp = 1023 - 60 + bits % 121;
        let x = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (exp << 52));
        match Rat::from_decimal(x) {
            Some(r) => {
                tried += 1;
                assert_eq!(r.to_f64(), x, "{x:?} → {r:?} → {:?}", r.to_f64());
            }
            None => declined += 1,
        }
    }
    // Nothing declines over this exponent range; the guard is here so a corpus that
    // drifted out of i128 could not turn this test vacuous.
    assert!(tried > 100_000, "{tried} tried, {declined} declined");
}

/// The same property, but over **every** finite f64 rather than the CAD-plausible
/// band: subnormals, `MIN_POSITIVE`, `MAX`, and the whole exponent range in between.
/// Where the power of ten leaves i128 the answer is `None` and the caller keeps its
/// f64 — what must never happen is a `Some` that realizes back to a *different*
/// number, because that would move a coordinate this cell has no business moving.
#[test]
fn from_decimal_never_lies_anywhere_in_the_finite_range() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut corpus = vec![
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::MIN,
        f64::from_bits(1), // the smallest subnormal
        1.1,
        7.7,
        1e-30,
        1e300,
    ];
    for _ in 0..200_000 {
        let bits = next();
        // Every exponent a finite f64 can have, subnormals included.
        let exp = (bits % 2047) << 52;
        corpus.push(f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | exp));
    }
    let mut declined = 0;
    for x in corpus {
        match Rat::from_decimal(x) {
            Some(r) => assert_eq!(r.to_f64(), x, "{x:e} → {r:?}"),
            None => declined += 1,
        }
    }
    assert!(
        declined > 0,
        "the i128 limit should bite somewhere out here"
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(Rat::from_decimal(bad), None, "{bad}");
    }
}

/// Where the i128 limit actually bites, pinned so the doc comment stays true and a
/// caller can tell whether its dimensions are anywhere near it. They are not: a CAD
/// model in millimetres sits around `1e0`, sixty orders of magnitude inside.
///
/// The two edges differ because a 17-digit mantissa spends 16 of its own powers of
/// ten on the fraction, so it runs out on the small side first. Note that exponent
/// *notation* is no obstacle in itself — `1e-30` is accepted; only the magnitude is.
#[test]
fn the_decimal_window_is_wide_and_its_edges_are_where_they_should_be() {
    let at = |s: &str| Rat::from_decimal(s.parse::<f64>().unwrap()).is_some();
    for e in -22..=38 {
        assert!(at(&format!("1.2345678901234567e{e}")), "17 digits at 1e{e}");
    }
    assert!(!at("1.2345678901234567e-23"), "17 digits at 1e-23");
    assert!(!at("1.2345678901234567e39"), "17 digits at 1e39");
    assert!(
        at("1e-38") && at("1e-30") && !at("1e-39"),
        "a short decimal reaches further down"
    );
}

/// The ends of the `i128` range, where the two operands' bounds differ.
///
/// `i128::MIN.unsigned_abs()` is `2¹²⁷` **exactly**, one past what the denominator
/// may be — an asymmetry easy to assert away and, when it was, a debug-only panic
/// on a value the algorithm handles correctly. The answers here are exact powers of
/// two and their neighbours, so they can be written down rather than approximated.
#[test]
fn the_ends_of_the_range_realize_exactly() {
    assert_eq!(Rat::from_int(i128::MIN).to_f64(), -(2f64.powi(127)));
    assert_eq!(Rat::from_int(i128::MAX).to_f64(), 2f64.powi(127)); // rounds up to 2¹²⁷
    assert_eq!(Rat::new(i128::MIN, 2).unwrap().to_f64(), -(2f64.powi(126)));
    assert_eq!(Rat::new(i128::MIN, i128::MAX).unwrap().to_f64(), -1.0);
    assert_eq!(Rat::new(1, i128::MAX).unwrap().to_f64(), 2f64.powi(-127));
    // Every one is finite and signed the way its numerator is.
    for (n, d) in [
        (i128::MIN, 3),
        (i128::MAX, 7),
        (-1, i128::MAX),
        (i128::MIN + 1, 1),
    ] {
        let q = Rat::new(n, d).unwrap().to_f64();
        assert!(q.is_finite() && (q < 0.0) == (n < 0), "{n}/{d} -> {q:e}");
    }
}

/// **The identity this whole cell exists for.** A dimension split into two and
/// stacked must land where the undivided one does. Lifting the f64 *values* cannot
/// give that — the drift is already inside them — and reading the decimals can.
#[test]
fn a_split_dimension_stacks_back_to_the_whole_one() {
    let split = Rat::from_decimal(1.1)
        .unwrap()
        .checked_add(Rat::from_decimal(6.6).unwrap())
        .unwrap();
    assert_eq!(split, Rat::new(77, 10).unwrap());
    assert_eq!(split, Rat::from_decimal(7.7).unwrap());
    assert_eq!(split.to_f64(), 7.7);
    // The f64 arithmetic this replaces, and the exact-but-binary lift that does not
    // help either — both land an ulp away.
    assert_ne!(1.1 + 6.6, 7.7);
    assert_ne!(
        Rat::try_from_f64(1.1)
            .unwrap()
            .checked_add(Rat::try_from_f64(6.6).unwrap())
            .unwrap(),
        Rat::try_from_f64(7.7).unwrap()
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Nearest-ness, checked against the definition rather than against a second
    /// implementation of it: `q` is the nearest f64 to `n/d` exactly when no
    /// neighbour of `q` is closer, and the comparison `|n/d − a/b| ≤ |n/d − c/e|`
    /// is decidable in exact rationals. Held to the range where those stay in i128.
    #[test]
    fn to_f64_lands_on_the_nearest_f64(n in -(1i128 << 60)..(1i128 << 60), d in 1i128..(1i128 << 60)) {
        prop_assume!(n != 0);
        let r = Rat::new(n, d).unwrap();
        let q = r.to_f64();
        let zero = Rat::from_int(0);
        let err = |y: f64| {
            let e = Rat::try_from_f64(y).and_then(|yr| r.checked_sub(yr))?;
            if e < zero { zero.checked_sub(e) } else { Some(e) }
        };
        // ★ The *check* has a range, and it is narrower than the generator: comparing exactly
        // needs `r`'s denominator times `q`'s power of two, which for a small enough value
        // leaves `i128` (measured: `n = 1, d = 500930446045` — 2^39 · 2^91). The two
        // neighbours below have always skipped on the same limit; `here` said `expect` and
        // so turned a limit of the instrument into a failure of the thing measured.
        prop_assume!(err(q).is_some());
        let here = err(q).expect("just assumed representable");
        for nb in [f64::from_bits(q.to_bits() + 1), f64::from_bits(q.to_bits() - 1)] {
            if let Some(there) = err(nb) {
                prop_assert!(here <= there, "{r:?}: {q:?} is not nearest ({nb:?} is closer)");
            }
        }
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// ★★★★★ **The total predicate must answer what the checked one answered.**
    ///
    /// `parallel_rat` replaced a checked-`Rat` cross that declined on overflow, and the
    /// change is only sound if the two agree wherever the old one could speak at all. The
    /// old spelling is kept here as the oracle — an independent derivation, not a call back
    /// into the code under test — and compared on every input it can answer.
    ///
    /// ★ **The width is drawn, not fixed**, because both halves of the range have to be
    /// visited: at `bits ≈ 30` the checked cross answers everything (so the two must agree),
    /// past `bits ≈ 63` its products leave `i128` and it falls silent (so the total one is
    /// alone). Measured with a counter before this spelling settled — a fixed narrow range
    /// left the oracle answering 100% of the time (differential, no coverage of the new
    /// behaviour) and a fixed wide one left it silent 100% of the time (coverage, no
    /// differential). Either way the test would have been half a test.
    #[test]
    fn the_total_cross_answers_what_the_checked_one_did(
        m in 0u32..40,
        xs in prop::array::uniform6(-(1i128 << 20)..(1i128 << 20)),
        ds in prop::array::uniform6(1i128..(1i128 << 20)),
    ) {
        // ★ The denominator's width is what decides whether the checked route survives, and
        // it is **drawn** so both halves of the range are visited. The scale is a power of
        // **three**, not two: scaling numerator and denominator by the same power of two
        // reduces straight back out (measured — the first spelling did exactly that and left
        // every case narrow), so the fraction has to be widened where the gcd cannot undo it.
        let scale = 3i128.saturating_pow(m);
        let r = |i: usize| {
            Rat::new(xs[i], (2 * ds[i] + 1).saturating_mul(scale)).unwrap()
        };
        let (a, b) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)]);
        // The retired spelling, verbatim: checked `Rat`, `None` on overflow.
        let checked = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<bool> {
            let term = |i: usize, j: usize| -> Option<Rat> {
                x[i].checked_mul(y[j])?.checked_sub(x[j].checked_mul(y[i])?)
            };
            let zero = Rat::from_int(0);
            Some(term(1, 2)? == zero && term(2, 0)? == zero && term(0, 1)? == zero)
        };
        if let Some(want) = checked(&a, &b) {
            prop_assert_eq!(parallel_rat(&a, &b), want, "{:?} x {:?}", a, b);
        }
        // Parallelism is symmetric and scale-invariant — properties the checked route could
        // not always demonstrate, and the total one must. (The scaling is asserted only when
        // every component survived it: a `Rat` whose numerator is already near the ceiling
        // has no ×3, and asserting through an `unwrap_or` fallback would compare a *different
        // vector* — the shape of a test that measures nothing.)
        prop_assert_eq!(parallel_rat(&a, &b), parallel_rat(&b, &a));
        let three = Rat::from_int(3);
        let scaled: Option<Vec<Rat>> = a.iter().map(|c| c.checked_mul(three)).collect();
        if let Some(s) = scaled {
            let s = [s[0], s[1], s[2]];
            prop_assert_eq!(parallel_rat(&s, &b), parallel_rat(&a, &b), "scale-invariance");
        }
    }

    /// **The clearance predicate must not drop the point's scale** — the negative control
    /// for the one decision in this section that is not a matter of taste.
    ///
    /// `(n·p + d)² − r²|n|²` is *not* homogeneous in `p`: clearing `p`'s denominators
    /// multiplies `n·p` while leaving `d` where it was, which states a different
    /// proposition. The naive spelling is written out here and compared with the truth in
    /// exact rationals; on mixed denominators the two part company (measured: 1.3% of
    /// random cases, so this test bites without needing a hand-picked adversary).
    #[test]
    fn clearance_keeps_the_points_scale(
        ns in prop::array::uniform3(-10i128..10),
        ps in prop::array::uniform3(-1000i128..1000),
        pd in prop::array::uniform3(1i128..(1i128 << 24)),
        rn in 1i128..100,
        rd in 1i128..100_000,
    ) {
        // ★ The generator is aimed at the **straddle**: an integer normal through the
        // origin, a point a hair off the plane (denominators up to 2²⁴), and a radius of
        // comparable size — so the answer is genuinely `Negative` about as often as
        // `Positive`. A generator whose points sit far outside every radius would compare
        // two implementations that always say `Positive`, which is how the first spelling
        // of this test passed while the scale-dropping defect was installed (measured: the
        // probe went green here and was caught only by a fixture two crates away).
        let coeffs: [Rat; 4] = [
            Rat::from_int(ns[0]),
            Rat::from_int(ns[1]),
            Rat::from_int(ns[2]),
            Rat::from_int(0),
        ];
        prop_assume!(coeffs[..3].iter().any(|c| *c != Rat::from_int(0)));
        let p: [Rat; 3] = core::array::from_fn(|i| Rat::new(ps[i], pd[i]).unwrap());
        let r = Rat::new(rn, rd).unwrap();

        // The truth, in exact rationals — no lifting, no scales, just the definition.
        let truth = (|| -> Option<Orient> {
            let dot = (0..3).try_fold(Rat::from_int(0), |a, i| {
                a.checked_add(coeffs[i].checked_mul(p[i])?)
            })?.checked_add(coeffs[3])?;
            let nn = (0..3).try_fold(Rat::from_int(0), |a, i| {
                a.checked_add(coeffs[i].checked_mul(coeffs[i])?)
            })?;
            let val = dot.checked_mul(dot)?
                .checked_sub(r.checked_mul(r)?.checked_mul(nn)?)?;
            Some(match val.cmp(&Rat::from_int(0)) {
                core::cmp::Ordering::Less => Orient::Negative,
                core::cmp::Ordering::Equal => Orient::Zero,
                core::cmp::Ordering::Greater => Orient::Positive,
            })
        })();
        if let Some(want) = truth {
            prop_assert_eq!(point_plane_clearance_rat(&coeffs, &p, &BigRat::square_of(r)), want,
                "coeffs {:?} p {:?} r {:?}", coeffs, p, r);
        }
    }

    /// The dot sign is the sign of the rational dot product, wherever that can be formed at
    /// all — and is total where it cannot.
    #[test]
    fn dot_sign_is_the_sign_of_the_dot(
        xs in prop::array::uniform6(-(1i128 << 30)..(1i128 << 30)),
        ds in prop::array::uniform6(1i128..(1i128 << 30)),
    ) {
        let r = |i: usize| Rat::new(xs[i], ds[i]).unwrap();
        let (a, b) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)]);
        let checked = (0..3).try_fold(Rat::from_int(0), |acc, i| {
            acc.checked_add(a[i].checked_mul(b[i])?)
        });
        if let Some(v) = checked {
            let want = match v.cmp(&Rat::from_int(0)) {
                core::cmp::Ordering::Less => Orient::Negative,
                core::cmp::Ordering::Equal => Orient::Zero,
                core::cmp::Ordering::Greater => Orient::Positive,
            };
            prop_assert_eq!(dot_sign_rat(&a, &b), want);
        }
        // Perpendicularity is symmetric, whatever the widths.
        prop_assert_eq!(dot_sign_rat(&a, &b), dot_sign_rat(&b, &a));
    }

    /// A vector is parallel to itself, to its multiples, and to zero — and **not** to a
    /// vector off its line. Without this the test above passes for a predicate that always
    /// answers `false` on the inputs the oracle cannot check.
    #[test]
    fn parallel_knows_a_line_from_a_plane(
        xs in prop::array::uniform3(-(1i128 << 40)..(1i128 << 40)),
        ds in prop::array::uniform3(1i128..(1i128 << 40)),
        k in 1i128..(1i128 << 20),
    ) {
        let v = [
            Rat::new(xs[0], ds[0]).unwrap(),
            Rat::new(xs[1], ds[1]).unwrap(),
            Rat::new(xs[2], ds[2]).unwrap(),
        ];
        let zero = [Rat::from_int(0); 3];
        prop_assert!(parallel_rat(&v, &v));
        prop_assert!(parallel_rat(&v, &zero), "zero is parallel to everything");
        let scaled: [Rat; 3] = core::array::from_fn(|i| {
            Rat::new(xs[i], ds[i]).unwrap().checked_mul(Rat::new(k, 1).unwrap()).unwrap_or(v[i])
        });
        prop_assert!(parallel_rat(&v, &scaled), "a multiple stays on the line");
        // A vector with one coordinate moved off the line is not parallel — unless `v` was
        // itself degenerate in that coordinate's plane, which the cross decides exactly.
        let mut off = v;
        off[0] = v[0].checked_add(Rat::from_int(1)).unwrap_or(v[0]);
        off[1] = v[1].checked_sub(Rat::from_int(1)).unwrap_or(v[1]);
        let cross_zero = {
            let t = |i: usize, j: usize| {
                v[i].checked_mul(off[j])
                    .and_then(|a| a.checked_sub(v[j].checked_mul(off[i])?))
            };
            match (t(1, 2), t(2, 0), t(0, 1)) {
                (Some(a), Some(b), Some(c)) => {
                    let z = Rat::from_int(0);
                    Some(a == z && b == z && c == z)
                }
                _ => None,
            }
        };
        if let Some(want) = cross_zero {
            prop_assert_eq!(parallel_rat(&v, &off), want);
        }
    }

    /// ★★★★★ **The two derivations must be the same function.**
    ///
    /// `plane_name_exact` runs the `Rat` route first and only falls back, so wherever the
    /// narrow one answers, the wide one is never consulted — and an error in it would sit
    /// there unseen until the day it *is* consulted, on inputs no test covers. Calling both on
    /// the same inputs is the only way to say they agree.
    ///
    /// ★ `plane_name_big` is `pub(crate)` for exactly this reason: a fallback hidden behind
    /// its filter cannot be tested against it.
    #[test]
    fn the_wide_derivation_answers_what_the_narrow_one_does(
        xs in prop::array::uniform9(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform9(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = (
            [r(0), r(1), r(2)],
            [r(3), r(4), r(5)],
            [r(6), r(7), r(8)],
        );
        let narrow = plane_through_points(a, b, c);
        let wide = plane_name_big(a, b, c);
        if let Some(n) = narrow {
            // ★ The invariant rides along: an answer the narrow route reached fits `i128`
            // by construction, so the wide route must store it `Narrow` — same value, same
            // representation, structural equality.
            prop_assert_eq!(wide.clone(), Some(PlaneName::Narrow(n)),
                "narrow answered but wide disagrees");
        }
        // ★ And the wide one, whenever it answers narrowly at all, answers about a plane
        // these points are actually on — checked in the rationals, no tolerance.
        //
        // ★★★★★ **The residual can overflow even when the name and the points both fit**, and
        // this test found that by asserting it could not. `c · p` multiplies a canonical
        // coefficient by a point coordinate, so it needs the *sum* of their widths — which is
        // exactly the population `Model::push_surface_with_coeffs` cannot verify either. An
        // unevaluable check is not a failed one, here as there: skip it, never fail on it.
        if let Some(w) = wide.as_ref().and_then(|n| n.narrow()) {
            for p in [a, b, c] {
                let residual = (|| {
                    let mut acc = w[3];
                    for k in 0..3 {
                        acc = acc.checked_add(w[k].checked_mul(p[k])?)?;
                    }
                    Some(acc)
                })();
                if let Some(acc) = residual {
                    prop_assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
                }
            }
        }
    }
}

/// ★★★★ **The width the narrow route cannot reach, and the wide one can.**
///
/// A triple whose reduced denominators are coprime — one a power of two, one a power of five,
/// which is what decimal arithmetic produces once it reduces — needs their product to state
/// the plane, and `(b − a) × (c − a)` needs it squared. `plane_through_points` gives up there;
/// the answer is small, and this is the case that says so.
#[test]
fn a_plane_the_narrow_route_gives_up_on_is_still_named() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let a = [r(1, 1 << 53), r(0, 1), r(0, 1)];
    let b = [r(0, 1), r(1, 5i128.pow(23)), r(0, 1)];
    let c = [r(0, 1), r(0, 1), r(1, (1 << 40) * 5i128.pow(11))];
    assert_eq!(
        plane_through_points(a, b, c),
        None,
        "the narrow route was expected to overflow on coprime denominators"
    );
    let name = plane_name_exact(a, b, c).expect("the wide route names it");
    // The canonical answer here is small (the doc above says so) — the invariant demands it
    // come back `Narrow`.
    let wide = name.narrow().expect("a small answer must be stored Narrow");
    for p in [a, b, c] {
        let mut acc = wide[3];
        for k in 0..3 {
            acc = acc
                .checked_add(wide[k].checked_mul(p[k]).expect("no overflow"))
                .expect("no overflow");
        }
        assert_eq!(acc, Rat::from_int(0), "a point is off the derived plane");
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// ★★★★★ **The two solves must be the same function** — the point-side twin of
    /// `the_wide_derivation_answers_what_the_narrow_one_does`, and for the same reason: a
    /// wide route that is only ever consulted where the narrow one declined would carry an
    /// error unseen until the day it is consulted. Calling both on the same inputs is the
    /// only way to say they agree.
    #[test]
    fn the_wide_solve_answers_what_the_narrow_one_does(
        xs in prop::array::uniform12(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform12(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let rows: [[Rat; 4]; 3] =
            core::array::from_fn(|k| core::array::from_fn(|j| r(4 * k + j)));
        let names = rows.map(PlaneName::Narrow);
        let wide = three_planes_big([&names[0], &names[1], &names[2]]);

        if let Some(n) = three_planes_rat(rows) {
            // ★ The invariant rides along: a point the narrow route reached fits `Rat` by
            // construction, so the wide route must store it `Narrow` — same value, same
            // representation, structural equality.
            prop_assert_eq!(wide.clone(), Some(MeetPoint::Narrow(n)),
                "narrow answered but wide disagrees");
        }

        // ★★ And wherever the wide one answers **at all**, the point it names is on all three
        // planes — checked in the rationals, no tolerance. This half covers the `Wide` arm,
        // which the differential above cannot reach: the narrow route is silent there by
        // definition, so agreement says nothing and only the residual does.
        if let Some(w) = wide.as_ref() {
            use num_bigint::BigInt;
            use num_rational::Ratio;
            use num_traits::Zero;
            let big = |n: i128, d: i128| Ratio::new(BigInt::from(n), BigInt::from(d));
            let coord = |i: usize| -> Ratio<BigInt> {
                match w {
                    MeetPoint::Narrow(p) => big(p[i].numer(), p[i].denom()),
                    MeetPoint::Wide(p) => Ratio::new(p[i].0.clone(), p[i].1.clone()),
                }
            };
            let point = [coord(0), coord(1), coord(2)];
            for row in &rows {
                let mut acc = big(row[3].numer(), row[3].denom());
                for (j, x) in point.iter().enumerate() {
                    acc += big(row[j].numer(), row[j].denom()) * x;
                }
                prop_assert!(acc.is_zero(), "the solved point is off one of the planes");
            }
        }
    }
}

/// The two remaining meanings of [`three_planes_rat`]'s `None`, each still honest:
/// no unique point (parallel planes), and a point that truly does not fit `Rat` — which
/// [`three_planes_big`] tells apart by answering `Wide`.
#[test]
fn the_solve_still_declines_what_it_should() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // Parallel pair: x = 0 and x = 1 — no unique point, both routes say so.
    let parallel = [
        [r(1, 1), r(0, 1), r(0, 1), r(0, 1)],
        [r(1, 1), r(0, 1), r(0, 1), r(-1, 1)],
        [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
    ];
    assert_eq!(
        three_planes_rat(parallel),
        None,
        "parallel planes meet nowhere"
    );
    // Narrow rows whose meeting point is wider than `Rat`: x + y = 2⁻¹⁰⁰, x − y = 5⁻⁵⁰
    // put a denominator of 2¹⁰¹·5⁵⁰ (~217 bits) on x. `three_planes_rat` declines;
    // the wide twin answers, and answers `Wide` — the causes stay told apart.
    let wide_point = [
        [r(1, 1), r(1, 1), r(0, 1), r(-1, 1 << 100)],
        [r(1, 1), r(-1, 1), r(0, 1), r(-1, 5i128.pow(50))],
        [r(0, 1), r(0, 1), r(1, 1), r(0, 1)],
    ];
    assert_eq!(
        three_planes_rat(wide_point),
        None,
        "a point no `Rat` can hold is a decline, not an answer"
    );
    let names = wide_point.map(PlaneName::Narrow);
    assert!(
        matches!(
            three_planes_big([&names[0], &names[1], &names[2]]),
            Some(MeetPoint::Wide(_))
        ),
        "the twin names the cause: the point exists and is wide"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Totality: wherever the integer core answers `Narrow`, [`three_planes_rat`] answers
    /// the same. ★ A **wiring** lock, stated as such: after the refactor both routes
    /// converge on `three_planes_int`, so this pins the `.or_else` plumbing (lift,
    /// `narrow()` return) rather than serving as an independent oracle — that role belongs
    /// to `an_overflowing_intermediate_no_longer_costs_the_answer`'s hand-known point and
    /// to the agreement direction the ops-side `point_width` invariant keeps.
    #[test]
    fn the_solve_answers_wherever_the_answer_is_narrow(
        rows in proptest::array::uniform3(proptest::array::uniform4((-9i128..=9, 1u32..=60))),
    ) {
        let rat_rows: [[Rat; 4]; 3] =
            rows.map(|row| row.map(|(n, e)| Rat::new(n, 1i128 << e).unwrap()));
        let names = rat_rows.map(PlaneName::Narrow);
        let big = three_planes_big([&names[0], &names[1], &names[2]]);
        let expect = match &big {
            Some(MeetPoint::Narrow(p)) => Some(*p),
            _ => None, // no unique point, or truly wide — the honest declines
        };
        prop_assert_eq!(three_planes_rat(rat_rows), expect);
    }
}

/// ★★★★★ **The negative control for [`MeetPoint::width_bits`]** — without it, a corpus that
/// reports "nothing over 127 bits" is indistinguishable from a dead probe.
///
/// A `Wide` carrier states a plane at `x = 2²⁰⁰`, which no `Rat` can hold; the meeting point
/// inherits that width and the fork has to take its other branch. So the instrument can say
/// the other thing, and `over127 = 0` in a measurement means the population, not the meter.
#[test]
fn the_width_meter_reports_a_point_no_rat_can_hold() {
    use num_bigint::BigInt;
    let far = BigInt::from(1) << 200;
    let wide = PlaneName::Wide([BigInt::from(1), BigInt::from(0), BigInt::from(0), -&far]);
    let r = |v: [i128; 4]| PlaneName::Narrow(v.map(Rat::from_int));
    let (py, pz) = (r([0, 1, 0, 0]), r([0, 0, 1, 0]));
    let found = three_planes_big([&wide, &py, &pz]).expect("three planes meet at (2²⁰⁰, 0, 0)");
    assert!(
        matches!(found, MeetPoint::Wide(_)),
        "a 201-bit coordinate cannot be Narrow"
    );
    assert_eq!(
        found.width_bits(),
        201,
        "the meter reads the coordinate's width"
    );
    assert_eq!(found.narrow(), None);
}

/// ★★★★★ **A canonical answer wider than `i128` is a name now, not a `None`** — and two
/// statements of that plane are one value. The
/// interning consequence is locked on the model side
/// (`a_wide_plane_interns_but_opens_no_shortcut` in nacre-topo).
#[test]
fn a_plane_too_wide_for_i128_is_named_wide() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // Cross-product terms multiply two coordinates' numerators, so two ~2^90 coprime
    // numerators push the canonical coefficients past i128 with no content to divide out.
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [r(big1, 3), r(big2, 7), r(0, 1)];
    let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
    let c = [r(1, 13), r(1, 17), r(1, 19)];
    assert_eq!(
        plane_through_points(a, b, c),
        None,
        "expected the narrow route to overflow"
    );
    let name = plane_name_exact(a, b, c).expect("collinear it is not — it must be named");
    // ★ Wide carries identity only: the arithmetic shortcuts' door stays shut
    // (a wide name must NOT open a frame).
    assert!(
        name.narrow().is_none(),
        "an answer past i128 must be stored Wide"
    );
    // ★ Same plane, two spellings (a permuted triple) — one value, structurally.
    let permuted = plane_name_exact(b, c, a).expect("the same plane, permuted");
    assert_eq!(
        name, permuted,
        "two statements of one wide plane must be one value"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// **The two orientation arms must be the same function** — same argument as
    /// `the_wide_derivation_answers_what_the_narrow_one_does`: `orient2d_rat` runs the `Rat`
    /// route first, so wherever it answers the `BigInt` arm is never consulted, and an error
    /// there would wait for the first overflowing input. `orient2d_big` is `pub(crate)` for
    /// exactly this call.
    #[test]
    fn the_big_orientation_answers_what_the_narrow_one_does(
        xs in prop::array::uniform6(-(1i64 << 20)..(1i64 << 20)),
        ds in prop::array::uniform6(1i64..(1i64 << 20)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = ([r(0), r(1)], [r(2), r(3)], [r(4), r(5)]);
        // Small operands: the narrow route always answers, so this compares the arms.
        prop_assert_eq!(orient2d_rat(a, b, c), orient2d_big(a, b, c));
    }
}

/// ★★ **A `Wide` name still judges its point** — the population `plane_residual_sign` exists
/// for. Same fixture as `a_plane_too_wide_for_i128_is_named_wide`: canonical coefficients
/// past `i128`, no `Rat` route to fall back on.
#[test]
fn a_residual_against_a_wide_name_still_gets_its_sign() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [r(big1, 3), r(big2, 7), r(0, 1)];
    let b = [r(-big2, 5), r(big1, 11), r(0, 1)];
    let c = [r(1, 13), r(1, 17), r(1, 19)];
    let name = plane_name_exact(a, b, c).expect("a genuine plane");
    assert!(name.narrow().is_none(), "the fixture must actually be Wide");
    for p in [a, b, c] {
        assert_eq!(plane_residual_sign(&name, p), 0, "a naming point is on it");
    }
    // A point off the plane along ±z (the naming triangle is not vertical: a and b span
    // z = 0 and c leaves it, so the normal has a z-component): opposite pushes must get
    // opposite, nonzero signs.
    let one = Rat::from_int(1);
    let up = [c[0], c[1], c[2].checked_add(one).unwrap()];
    let down = [c[0], c[1], c[2].checked_sub(one).unwrap()];
    let (su, sd) = (
        plane_residual_sign(&name, up),
        plane_residual_sign(&name, down),
    );
    assert_ne!(su, 0);
    assert_eq!(su, -sd, "opposite sides, opposite signs");
}

/// The core rational-representation property in miniature: exact rational accumulation does not
/// drift, where the f64 control does. `(1/10)` summed ten times is exactly
/// `1`, but `0.1_f64` summed ten times is not `1.0`.
#[test]
fn rational_accumulation_is_exact_where_f64_drifts() {
    let tenth = Rat::new(1, 10).unwrap();
    let mut acc = Rat::from_int(0);
    for _ in 0..10 {
        acc = acc.checked_add(tenth).unwrap();
    }
    assert_eq!(acc, Rat::from_int(1));

    let mut f = 0.0_f64;
    for _ in 0..10 {
        f += 0.1;
    }
    assert_ne!(f, 1.0); // 0.9999999999999999 — the drift Rat avoids
}

/// The "thin film" example: `1.1 × 7` must be exactly `7.7`. In rationals
/// `11/10 × 7 = 77/10`; in f64 `1.1 * 7.0` is not `7.7`.
#[test]
fn one_point_one_times_seven_is_exact() {
    let a = Rat::new(11, 10).unwrap();
    let seven = Rat::from_int(7);
    assert_eq!(a.checked_mul(seven).unwrap(), Rat::new(77, 10).unwrap());

    assert_ne!(1.1_f64 * 7.0, 7.7); // f64 cannot represent 7.7 exactly
}

/// Measurement: with fixed-width i128, chained coprime-denominator
/// accumulation *does* overflow (the finite-precision cliff handled by
/// downgrading). Summing `1/p` over successive primes forces the denominator
/// toward the primorial, which exceeds i128. This test pins two facts:
///   (a) a modest sum (first 8 primes) stays exact — normal use is fine;
///   (b) accumulation eventually overflows within the prime list — proving
///       the downgrade trigger fires, and reporting *where* (onset index).
#[test]
fn coprime_accumulation_overflows_and_reports_onset() {
    const PRIMES: [i128; 30] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113,
    ];

    // (a) first 8 primes stay exact.
    let mut acc = Rat::from_int(0);
    for &p in &PRIMES[..8] {
        acc = acc
            .checked_add(Rat::new(1, p).unwrap())
            .expect("first 8 primes must not overflow");
    }

    // (b) full accumulation eventually overflows; record the onset.
    let mut acc = Rat::from_int(0);
    let mut onset = None;
    let mut last_bits = 0;
    for (i, &p) in PRIMES.iter().enumerate() {
        match acc.checked_add(Rat::new(1, p).unwrap()) {
            Some(next) => {
                acc = next;
                last_bits = acc.bit_width();
            }
            None => {
                onset = Some((i, last_bits));
                break;
            }
        }
    }
    let (idx, bits) = onset.expect("i128 rational accumulation must overflow within 30 primes");
    // Measurement record (run with `-- --nocapture`).
    eprintln!(
        "[rational] coprime 1/p accumulation overflows at prime index {idx} (p={}); \
             denominator bit-width just before onset = {bits}",
        PRIMES[idx]
    );
    // Onset is comfortably past normal use and near the i128 ceiling (~127 bits).
    assert!(idx >= 8, "overflow onset {idx} should be past modest use");
    assert!(
        bits > 100,
        "denominator should be near the i128 ceiling at onset, got {bits} bits"
    );
}
