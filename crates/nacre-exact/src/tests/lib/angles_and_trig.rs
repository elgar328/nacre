//! Rational angles, trigonometric and inverse-sqrt realization with their error bounds, and `Rat`
//! algebra.

use super::*;

/// The depth these tests realize at. **Test-local on purpose**: production has no default
/// precision — it is computed from the model — so a constant here must not be reachable from
/// outside the fixtures that chose it.
const GT_PREC: usize = 160;

/// Closure: a rational angle accumulates exactly, so a full turn lands back
/// on exactly `0` — while the f64 control drifts. `360/7` degrees added seven
/// times is exactly `360 → 0`; `360.0/7.0` summed seven times in f64 is not
/// `360.0`.
#[test]
fn rational_angle_closes_exactly_where_f64_drifts() {
    let step = Rat::new(360, 7).unwrap();
    let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
    for _ in 0..7 {
        a = a.checked_add(step).unwrap();
    }
    assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());

    let mut f = 0.0_f64;
    for _ in 0..7 {
        f += 360.0 / 7.0;
    }
    assert_ne!(f, 360.0); // f64 drifts off the full turn
}

/// Many small rational steps also close: `1/3` degree added 1080 times is
/// exactly one full turn → `0`.
#[test]
fn many_small_steps_close_to_zero() {
    let step = Rat::new(1, 3).unwrap();
    let mut a = Angle::from_deg(Rat::from_int(0)).unwrap();
    for _ in 0..1080 {
        a = a.checked_add(step).unwrap();
    }
    assert_eq!(a, Angle::from_deg(Rat::from_int(0)).unwrap());
}

/// The angle value stays exact, and its cos/sin realization is f64 — but a *correctly rounded*
/// one, so the irrational-realization boundary now costs at most half an ulp rather than
/// whatever the platform's libm happened to do.
///
/// **`cos 45°` is the case to look at**: `√2/2` cannot be an f64, so the realization is
/// genuinely lossy, and it must land on the nearest f64 to the truth. Checked by asking a far
/// deeper realization whether anything is closer.
#[test]
fn realization_is_f64_while_angle_stays_exact() {
    let a = Angle::from_deg(Rat::from_int(45)).unwrap();
    assert_eq!(a.deg(), Rat::from_int(45)); // angle exact
    let (c, _) = a.cos_sin_f64();
    assert!(a.try_exact_cos_sin().is_none()); // √2/2 is not rational

    // Nothing is nearer: both neighbours are further from the deep truth than `c` is.
    let (deep, _) = a.cos_sin_at(512);
    let dist = |f: f64| BigFloat::from_f64(f, 512).sub(&deep, 512, HP_RM).abs();
    let here = dist(c);
    for nb in [
        f64::from_bits(c.to_bits() - 1),
        f64::from_bits(c.to_bits() + 1),
    ] {
        assert!(
            here.cmp(&dist(nb)).is_some_and(|s| s < 0),
            "a neighbour of {c:e} is nearer the truth"
        );
    }
}

/// A 90°-family angle yields exact rational cos/sin, so rotating a rational
/// point stays exact (tol 0): `(x, y)` rotated 90° is `(-y, x)`. A 45° angle
/// has no exact rational cos/sin (√2/2), so it returns `None` and would fall
/// to the f64/dd realization. Note the f64 path is *not* exact here:
/// `cos 90°` realizes to ~6e-17, not `0`.
#[test]
fn exact_cos_sin_only_for_quadrantal_angles() {
    let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
    let (c, s) = a90.try_exact_cos_sin().unwrap();
    assert_eq!((c, s), (Rat::from_int(0), Rat::from_int(1)));

    // Exact rotation of (3, 5) by 90° → (-5, 3), all rational (tol 0):
    // x' = x·cos − y·sin,  y' = x·sin + y·cos.
    let (x, y) = (Rat::from_int(3), Rat::from_int(5));
    let xr = x
        .checked_mul(c)
        .unwrap()
        .checked_sub(y.checked_mul(s).unwrap())
        .unwrap();
    let yr = x
        .checked_mul(s)
        .unwrap()
        .checked_add(y.checked_mul(c).unwrap())
        .unwrap();
    assert_eq!((xr, yr), (Rat::from_int(-5), Rat::from_int(3)));

    // 45° has no exact rational realization → None (falls to the rounded high-precision path).
    let a45 = Angle::from_deg(Rat::from_int(45)).unwrap();
    assert!(a45.try_exact_cos_sin().is_none());
    // The high-precision realization of 90° is *not* zero either — it is ~2⁻¹²⁸ — which is
    // exactly why the branch above exists rather than being an optimization.
    assert!(!a90.cos_sin_at(128).0.is_zero());
}

/// `cos_sin_f64` snaps the 90°-family to exact `0.0`/`±1.0`, and the exactly-representable
/// values Niven allows off that family come out exact too.
///
/// ★★ **`cos 60° == 0.5` is the visible proof that libm left.** `1/2` is one of the three
/// rational values a rational-degree cosine can take, and it *is* an f64 — but reaching it
/// through `(60.0 * PI / 180.0).cos()` does not land on it. Rounding the high-precision value
/// does.
#[test]
fn cos_sin_f64_is_exact_where_the_true_value_is_representable() {
    let deg = |d| Angle::from_deg(Rat::from_int(d)).unwrap();
    assert_eq!(deg(0).cos_sin_f64(), (1.0, 0.0));
    assert_eq!(deg(90).cos_sin_f64(), (0.0, 1.0));
    assert_eq!(deg(180).cos_sin_f64(), (-1.0, 0.0));
    assert_eq!(deg(270).cos_sin_f64(), (0.0, -1.0));
    for (d, want) in [(60, 0.5), (120, -0.5), (240, -0.5), (300, 0.5)] {
        assert_eq!(deg(d).cos_sin_f64().0, want, "cos {d}°");
    }
    for (d, want) in [(30, 0.5), (150, 0.5), (210, -0.5), (330, -0.5)] {
        assert_eq!(deg(d).cos_sin_f64().1, want, "sin {d}°");
    }
}

/// **The 90°-family realizes with no error at all, and that zero is load-bearing.**
///
/// `WitnessPoint::rotate_about` reads [`Angle::realization_error_of`] to decide whether a rotation
/// contributes any tolerance. For `cos`/`sin` in `{0, ±1}` the f64 values *are* the true ones,
/// and every product and difference downstream is exact too — which is why a quadrantal origin
/// rotation stays at tol 0 and an axis-aligned model never leaves the exact predicate path.
/// A nonzero here would not fail loudly; it would quietly move those models.
///
/// Everything else must report *something* — but only where the true value is *not*
/// representable. Since the realization became correctly rounded, `cos 60°` really is `0.5`
/// exactly, so the corpus below has to avoid the four angles where that happens or it would be
/// asserting a nonzero error that does not exist.
#[test]
fn only_the_quadrantal_family_realizes_exactly() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    let err = |a: Angle| {
        let (c, s) = a.cos_sin_f64();
        a.realization_error_of(c, s)
    };
    for d in [0, 90, 180, 270] {
        assert_eq!(err(deg(d, 1)), (0.0, 0.0), "{d} deg");
    }
    // ★ And the zero is of the *pair*, not of the angle: hand a 90°-family angle some other
    // realization and it must be measured like anything else. Otherwise a caller that got its
    // cos/sin from somewhere other than `cos_sin_f64` would be handed an exactness claim that
    // does not hold of what it is holding.
    let a90 = deg(90, 1);
    assert!(a90.realization_error_of(1e-17, 1.0).0 > 0.0);
    for (n, d) in [(1, 1), (37, 1), (45, 1), (337, 1), (1, 3), (359999, 1000)] {
        let (dc, ds) = err(deg(n, d));
        assert!(
            dc > 0.0 && ds > 0.0,
            "{n}/{d} deg claimed an exact realization"
        );
        // ε-scale: a bound this large would mean the measurement, not the platform, is wrong.
        assert!(
            dc < 64.0 * f64::EPSILON && ds < 64.0 * f64::EPSILON,
            "{n}/{d} deg: {dc:e}"
        );
    }
}

/// **The f64 read-out is a re-encoding, not a computation — checked by round trip.**
///
/// `set_precision(53, ToEven)` is where a value loses bits; everything after it is supposed to
/// be pure bookkeeping over the mantissa words. So any `f64` put in must come back out
/// unchanged. **The 32-bit-`Word` path is the one this is really for** — the mantissa is
/// assembled across two words there, wasm is a 32-bit target, and no amount of reading the
/// crate source substitutes for running it on the target that ships.
#[test]
fn the_f64_readout_round_trips() {
    let mut cases = vec![
        1.0,
        0.5,
        -0.5,
        0.9999999999999999,
        1e-300,
        -3.7e17,
        f64::MIN_POSITIVE,
    ];
    let mut st = 0x2545_F491_4F6C_DD1Du64;
    for _ in 0..2000 {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        // Any finite normal double; the exponent is squeezed into the normal range.
        let bits = (st & !(0x7ffu64 << 52)) | ((1 + (st >> 53) % 2045) << 52);
        let x = f64::from_bits(bits);
        if x.is_normal() {
            cases.push(x);
        }
    }
    for x in cases {
        let back = to_f64_exact(&BigFloat::from_f64(x, 128));
        assert_eq!(back, Some(x), "{x:e} did not survive the round trip");
    }
    assert_eq!(to_f64_exact(&BigFloat::from_f64(0.0, 128)), Some(0.0));
}

/// **The rounding check has to be able to say no.**
///
/// [`round_to_f64`] answers `None` when the interval straddles a rounding boundary, and that
/// branch is the whole reason the function is not just "round the midpoint". A check that can
/// only say yes is indistinguishable from no check — and there is a specific way to build one
/// here, by forming `mid ± rad` at the realization's own precision so the radius rounds away.
/// So: a radius wide enough to be undecidable must be refused, and a tight one accepted.
#[test]
fn the_rounding_check_refuses_an_undecidable_interval() {
    let a = Angle::from_deg(Rat::new(37, 1).unwrap()).unwrap();
    let (c, _) = a.cos_sin_bounded(128);
    let (c, rc) = (c.value, c.error);
    assert!(
        round_to_f64(&c, rc, 128).is_some(),
        "a 2^-128 radius is decidable"
    );
    // An ulp-wide radius cannot be: it reaches both neighbours.
    assert!(round_to_f64(&c, Mag::pow2(-52), 128).is_none());
    // And zero is the case no precision resolves — cos 90 is exactly 0, so its interval
    // straddles zero forever. This is why `cos_sin_f64` resolves the family first.
    let a90 = Angle::from_deg(Rat::from_int(90)).unwrap();
    for prec in [128usize, 256, 512] {
        let (c90, _) = a90.cos_sin_bounded(prec);
        let (c90, r90) = (c90.value, c90.error);
        assert!(
            round_to_f64(&c90, r90, prec).is_none(),
            "cos 90 became decidable at {prec}, which would make the quadrantal branch optional"
        );
    }
}

/// The escalation and fallback rungs, counted — over a corpus that reaches for them.
///
/// Near an axis `|cos|` is tiny while its error bound is absolute, so the relative radius grows
/// and 128 bits stops being enough; that is the only place the second rung is reachable. **The
/// count is reported rather than asserted nonzero**: with `Rat` bounded by `i128` an angle
/// cannot get closer to 90° than ~6e-39 degrees, so it is entirely possible that nothing in a
/// finite corpus needs it — but a silent zero and an unreachable branch look identical, and
/// this at least says which corpus produced the zero.
#[test]
fn the_escalation_rung_is_reachable() {
    let before = ROUND_ESCALATED.with_borrow(|c| *c);
    let mut asked = 0usize;
    for k in 1..400i128 {
        // Just off 90 degrees, by ever smaller amounts.
        for d in [
            10i128.pow(9),
            10i128.pow(18),
            10i128.pow(30),
            i128::MAX / 91,
        ] {
            if let Some(a) = Rat::new(90 * d + k, d).and_then(Angle::from_deg) {
                let (c, s) = a.cos_sin_f64();
                assert!(
                    c.is_finite() && s.is_finite(),
                    "{a:?} realized to a non-number"
                );
                // Whatever rung answered, the answer must still bound its own error.
                let (dc, ds) = a.realization_error_of(c, s);
                assert!(dc >= 0.0 && ds >= 0.0);
                asked += 1;
            }
        }
    }
    let after = ROUND_ESCALATED.with_borrow(|c| *c);
    eprintln!(
        "[round_to_f64] {} escalated to 256, {} fell back, over {asked} near-axis angles",
        after.0 - before.0,
        after.1 - before.1
    );
    assert!(asked > 500, "corpus shrank to {asked}");
    assert_eq!(
        after.1, before.1,
        "the undecidable fallback fired, which is derived not to"
    );
}

/// The memo answers the second ask, keys on the angle's *value* rather than its spelling — same
/// failure mode as [`TRIG`]'s, showing up as work done twice rather than a wrong answer — **and
/// keys on the realized pair**, which is the part that is about correctness rather than cost.
#[test]
fn the_realization_error_memo_keys_on_the_angle_and_the_pair() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    let entries = || F64_ERR.with_borrow(|m| m.len());
    // An angle no other test in this thread asks for.
    let a = deg(1234567, 9973);
    let (c, s) = a.cos_sin_f64();
    let before = entries();
    let first = a.realization_error_of(c, s);
    assert_eq!(entries(), before + 1, "the first ask must insert");
    assert_eq!(a.realization_error_of(c, s), first);
    assert_eq!(
        entries(),
        before + 1,
        "the second ask must be answered, not recomputed"
    );
    // The same angle, unreduced, with the same pair, is the same entry.
    assert_eq!(deg(2469134, 19946).realization_error_of(c, s), first);
    assert_eq!(
        entries(),
        before + 1,
        "the key is the angle, not its spelling"
    );
    // ★ A neighbouring realization of that same angle is a *different question* and must get
    // its own answer — an entry keyed on the angle alone would hand back `first`, an error
    // measured against a pair this caller is not holding.
    let nudged = f64::from_bits(c.to_bits() + 1);
    assert_ne!(a.realization_error_of(nudged, s), first);
    assert_eq!(entries(), before + 2, "the pair is part of the key");
}

/// `apply_point`/`apply_dir` realize a 90°-family rotation bit-exactly: no ~6e-17
/// spurious cross-term. (3,5,z) about Z by 90° → exactly (-5,3,z); a rational-pivot
/// rotation is exact too; a non-quadrantal angle is unchanged from the f64 path.
#[test]
fn apply_point_exact_for_quadrantal() {
    let iso = |d| {
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(d)).unwrap(),
        })
    };
    assert_eq!(iso(90).apply_point([3.0, 5.0, 7.0]), [-5.0, 3.0, 7.0]);
    assert_eq!(iso(180).apply_point([3.0, 5.0, 7.0]), [-3.0, -5.0, 7.0]);
    assert_eq!(iso(270).apply_point([3.0, 5.0, 7.0]), [5.0, -3.0, 7.0]);
    assert_eq!(iso(90).apply_dir([0.0, 1.0, 0.0]), [-1.0, 0.0, 0.0]);

    // Non-origin pivot (2,2): 90° maps (3,5)→pivot+R(1,3)=(2-3, 2+1)=(-1,3).
    let piv = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(2), Rat::from_int(2), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    });
    assert_eq!(piv.apply_point([3.0, 5.0, 0.0]), [-1.0, 3.0, 0.0]);

    // Non-quadrantal: the same arithmetic on the same realized pair, so this pins the *route*
    // (`px + u·c − v·s`, in that order) rather than the values.
    let a = Angle::from_deg(Rat::from_int(37)).unwrap();
    let (c, s) = a.cos_sin_f64();
    assert_eq!(
        iso(37).apply_point([3.0, 5.0, 0.0]),
        [3.0 * c - 5.0 * s, 3.0 * s + 5.0 * c, 0.0]
    );
}

/// H1.5: the arbitrary-precision cos/sin realization must be far more accurate
/// than any f64/double-double — the accuracy gate astro-float passes and the
/// double-double `twofloat` failed (its trig degraded to ~1e-16 near zero-
/// crossings). Error at rational angles must beat `2^-100` (~1e-30).
#[test]
fn high_precision_trig_meets_gate() {
    const GATE_EXP: i32 = -100; // 2^-100 ≈ 7.9e-31
    // (deg, exact_cos, exact_sin | None where irrational, e.g. sin 60°)
    let cases = [
        (0i128, 1.0, Some(0.0)),
        (60, 0.5, None),
        (90, 0.0, Some(1.0)),
        (120, -0.5, None),
        (180, -1.0, Some(0.0)),
        (270, 0.0, Some(-1.0)),
    ];
    let mut worst = i32::MIN;
    for (deg, exact_cos, exact_sin) in cases {
        let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
        let ec = hp_err_exp(&a.cos_sin_at(GT_PREC).0, exact_cos);
        assert!(ec < GATE_EXP, "cos {deg}° error 2^{ec} exceeds gate");
        worst = worst.max(ec);
        if let Some(s) = exact_sin {
            let es = hp_err_exp(&a.cos_sin_at(GT_PREC).1, s);
            assert!(es < GATE_EXP, "sin {deg}° error 2^{es} exceeds gate");
            worst = worst.max(es);
        }
    }
    eprintln!("[H1.5] astro-float {GT_PREC}-bit worst cos/sin error at rational angles: 2^{worst}");
}

/// A large angle normalizes in one step, not one step per turn.
///
/// Normalization used to subtract 360° in a loop, so `2⁶⁰` degrees needed ~3·10¹⁵ iterations —
/// the kernel did not reject that input, it stopped responding to it. Found by execution: a
/// trig corpus reached for a big numerator and the test never returned.
#[test]
fn a_huge_angle_normalizes_without_counting_turns() {
    let huge = Rat::new(1i128 << 60, 7).unwrap();
    let a = Angle::from_deg(huge).expect("a large rational angle is representable");
    assert!(a.deg() >= Rat::from_int(0) && a.deg() < Rat::from_int(360));
    // Same residue class as the input, so the reduction is `− 360k`, not a different angle.
    let back = a.deg().checked_sub(huge).unwrap();
    let turns = back.to_f64() / -360.0;
    assert_eq!(
        turns.fract(),
        0.0,
        "the reduction was not a whole number of turns"
    );
    // Negatives land in range too, and exactly on 0 at a full turn.
    assert_eq!(
        Angle::from_deg(Rat::from_int(-720)).unwrap().deg(),
        Rat::from_int(0)
    );
    assert_eq!(
        Angle::from_deg(Rat::new(-1, 2).unwrap()).unwrap().deg(),
        Rat::new(719, 2).unwrap()
    );
}

/// **The memo memoizes, and its key is the *value* of the angle rather than its spelling.**
///
/// Two failure modes, and neither shows up as a wrong answer — [`Angle::cos_sin_bounded`] is a
/// pure function either way, so a broken memo is only slow. That is exactly why it needs a
/// test: a rotated boolean spent 22% of its time re-running Ziv's loop for angles it had
/// already realized, and nothing but a measurement would say so.
///
/// - **It does not memoize** (a rewrite drops the lookup): asking twice would insert twice.
/// - ★ **The key fragments**: `90/1` and `180/2` are the same angle. If they hashed apart the
///   memo would still be *correct* and still show a high hit rate, while paying twice for every
///   angle a caller happened to spell in two ways. `Angle::from_deg` reduces through
///   `Rat::new`, so they are one entry — this is the guard on that.
///
/// Deltas, not absolute counts: the memo is a `thread_local`, and the harness gives each test
/// its own thread, but nothing here should depend on which tests ran first.
#[test]
fn the_trig_memo_keys_on_the_angle_not_its_spelling() {
    let deg = |n, d| Angle::from_deg(Rat::new(n, d).unwrap()).unwrap();
    // A precision no other test asks for, so this thread's map cannot be pre-warmed for it.
    let prec = 704;
    let before = trig_entries();
    let first = deg(37, 1).cos_sin_bounded(prec);
    assert_eq!(trig_entries(), before + 1, "the first ask must insert");

    let again = deg(37, 1).cos_sin_bounded(prec);
    assert_eq!(
        trig_entries(),
        before + 1,
        "the second ask must be answered, not recomputed"
    );
    assert_eq!(
        (first.0.value.clone(), first.0.error),
        (again.0.value.clone(), again.0.error),
        "and answered with the same value"
    );

    // ★ The same angle, spelled as an unreduced ratio, is the same entry.
    let spelled = deg(74, 2).cos_sin_bounded(prec);
    assert_eq!(
        trig_entries(),
        before + 1,
        "74/2 is 37/1: a second entry means the key is the spelling, not the angle"
    );
    assert_eq!(
        (first.0.value.clone(), first.0.error),
        (spelled.0.value.clone(), spelled.0.error)
    );

    // …and precision *is* part of the key: a different depth is a different answer, so reusing
    // an entry across depths would hand back coordinates realized at the wrong one.
    let deeper = deg(37, 1).cos_sin_bounded(prec + 64);
    assert_eq!(trig_entries(), before + 2, "precision must key the memo");
    assert_ne!(
        deeper.0.error, again.0.error,
        "a deeper realization has a smaller bound"
    );
}

/// **The seed of every error radius, checked against a realization far deeper than itself.**
///
/// [`Angle::cos_sin_bounded`] derives its bound from two things the crate does not promise in
/// writing: that `Consts::pi` is correctly rounded, and that `cos`/`sin` are too (they run
/// Ziv's loop, which is how one builds a correctly-rounded transcendental — but an
/// implementation detail, not a documented contract). If either weakens, every interval above
/// it is unsound, so the claim is measured: at each rung the ladder uses, the value must sit
/// within its own bound of the same value realized with 512 extra bits.
///
/// The reference is not independent code — it is the same routine at higher precision — so
/// this cannot catch an error that grows with precision in the same shape. What it does catch
/// is the failure that matters here: a bound that is simply too small.
#[test]
fn the_trig_bound_holds_against_a_far_deeper_realization() {
    // Angles spanning the quadrants, plus rationals with awkward denominators and one whose
    // numerator is large enough to exercise the `i128 → f64` term.
    let angles = [
        (0i128, 1i128),
        (30, 1),
        (45, 1),
        (60, 1),
        (90, 1),
        (135, 1),
        (180, 1),
        (271, 1),
        (359, 1),
        (1, 7),
        (22, 7),
        (1000, 3),
        (1, 1_000_000),
        // A denominator past 2⁵³ that is still exact in `f64` (a power of two), so the bound
        // must *not* charge it the conversion term.
        (1, 1i128 << 60),
        // …and one that genuinely does not round-trip, where the angle itself is only known
        // to a relative `2⁻⁵³` and no working precision can recover it.
        (1, (1i128 << 60) + 1),
    ];
    // Slack is tracked per rung so the two stories stay separable: a word-aligned precision
    // is delivered as asked, while `200` is silently rounded up to 256 and the extra bits show
    // up as slack that is astro-float's, not this bound's.
    let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
    for prec in [128usize, 200, 256, 512, 1024] {
        for (num, den) in angles {
            let a = Angle::from_deg(Rat::new(num, den).unwrap()).unwrap();
            let (c, s) = a.cos_sin_bounded(prec);
            let (bc, bs) = (c.error, s.error);
            let (c, s) = (c.value, s.value);
            let deep = prec + 512;
            let (rc, rs) = a.cos_sin_at(deep);
            for (got, reference, bound, what) in [(&c, &rc, bc, "cos"), (&s, &rs, bs, "sin")] {
                let diff = got.sub(reference, deep, HP_RM);
                let Some(de) = (if diff.is_zero() {
                    None
                } else {
                    diff.exponent()
                }) else {
                    continue; // exactly equal — nothing to bound
                };
                // `|diff| < 2^de`; the bound must be at least that.
                let observed = Mag::pow2(de as i64);
                assert!(
                    !bound.lt(observed),
                    "{what} {num}/{den}° at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                    bound.exp2()
                );
                // Track how much slack the bound carries, so a bound that is merely
                // enormous does not pass as a bound that is right.
                // Slack is only meaningful where the angle converts exactly. Where it does
                // not, the `2⁻⁵³` term dominates by design and the gap to the observed error
                // is the honest cost of an unrepresentable angle, not looseness.
                let angle_exact = (num as f64) as i128 == num && (den as f64) as i128 == den;
                // And only where the two realizations actually disagree above the
                // *reference's* own resolution. `cos 60° = 1/2` is exact at both precisions,
                // so their difference measures the reference, not this bound.
                let above_reference_noise = (de as i64) > -(deep as i64) + 8;
                if let (Some(be), true) = (bound.exp2(), angle_exact && above_reference_noise) {
                    let slack = be - de as i64;
                    let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                    if slack > e.0 {
                        *e = (
                            slack,
                            format!("{what} {num}/{den}°, bound 2^{be} vs error 2^{de}"),
                        );
                    }
                }
            }
        }
    }
    for (prec, (slack, at)) in &worst {
        eprintln!("[cip] {prec}-bit trig bound slack: 2^{slack}  ({at})");
    }
    for (&prec, (slack, at)) in &worst {
        if prec % 64 == 0 {
            assert!(
                *slack <= 16,
                "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                     that is a fudge factor wearing a derivation's clothes, not a tight bound"
            );
        }
    }
    // …and the odd precision out proves why the realization ladder's rungs (`LADDER` in
    // `nacre-ops`) are multiples of 64: asking for 200 bits buys 256, so ~56 bits of the result
    // are paid for and then claimed away.
    let (slack_200, _) = &worst[&200];
    assert!(
        (40..=80).contains(slack_200),
        "expected ~56 bits of unclaimed precision at 200 bits (the word-size round-up), got \
             2^{slack_200}"
    );
}

/// Squared lengths a frame normal actually produces, plus awkward ones.
///
/// `(0,0,1)` and `(1,1,0)` give `1` and `2`; a profile edge `(1,2)` gives a wall normal
/// `(2,−1,0)` and so `5`; the fractions are what a normal reduced by its own content leaves;
/// and the last two exercise the `i128 → f64` term at and past 2⁵³.
const INV_SQRT_CASES: [(i128, i128); 12] = [
    (1, 1),
    (2, 1),
    (3, 1),
    (5, 1),
    (4, 1),
    (9, 25),
    (1, 2),
    (13, 7),
    (1_000_003, 3),
    (1, 1_000_000),
    (1, 1i128 << 60),
    (1, (1i128 << 60) + 1),
];

#[test]
fn the_inverse_sqrt_bound_holds_against_a_far_deeper_realization() {
    let mut worst = std::collections::BTreeMap::<usize, (i64, String)>::new();
    for prec in [128usize, 256, 512] {
        for (num, den) in INV_SQRT_CASES {
            let v = Rat::new(num, den).unwrap();
            let HpBounded {
                value: z,
                error: bound,
            } = inv_sqrt_bounded(v, prec).unwrap();
            let deep = prec + 512;
            let rz = inv_sqrt_bounded(v, deep).unwrap().value;
            let diff = z.sub(&rz, deep, HP_RM);
            let Some(de) = (if diff.is_zero() {
                None
            } else {
                diff.exponent()
            }) else {
                continue; // exactly equal — nothing to bound
            };
            // `|diff| < 2^de`; the bound must be at least that.
            assert!(
                !bound.lt(Mag::pow2(de as i64)),
                "1/sqrt({num}/{den}) at {prec} bits: error 2^{de} exceeds its bound 2^{:?}",
                bound.exp2()
            );
            // ★★ Slack is only meaningful where the observation actually *measures* this
            // bound. A `prec`-bit realization of a value of magnitude `2^zexp` must carry an
            // error near `2^(zexp − prec)`; when the observed gap is far below that, the two
            // realizations agreed better than either one's own accuracy warrants — measured,
            // `1/√(1/(2⁶⁰+1))` at 128 bits agrees to 154 bits where 98 is all that is earned,
            // and `1/√(1/10⁶)` is exactly `1000` at every precision. Those samples understate
            // the true error, so counting them would read a *lucky observation* as a loose
            // bound. Same shape as the trig test's reference-noise filter.
            let earned = z.exponent().unwrap_or(0) as i64 - prec as i64 - 8;
            if let (Some(be), true) = (bound.exp2(), (de as i64) >= earned) {
                let slack = be - de as i64;
                let e = worst.entry(prec).or_insert((i64::MIN, String::new()));
                if slack > e.0 {
                    *e = (slack, format!("{num}/{den}, bound 2^{be} vs error 2^{de}"));
                }
            }
        }
    }
    for (prec, (slack, at)) in &worst {
        eprintln!("[cip] {prec}-bit inverse-sqrt bound slack: 2^{slack}  ({at})");
    }
    for (prec, (slack, at)) in &worst {
        assert!(
            *slack <= 8,
            "at {prec} bits the bound is 2^{slack} above the worst observed error ({at}) — \
                 that is a fudge factor wearing a derivation's clothes, not a tight bound"
        );
    }
}

/// ★ The claim the axis-aligned corpus rests on: a frame whose normal is a coordinate
/// direction realizes its `1/|n|` **exactly**, so it never reaches arbitrary precision and
/// carries no realization error at all.
#[test]
fn an_axis_aligned_frame_needs_no_arbitrary_precision() {
    // `n·n` for (0,0,1) and (2,0,0) — an axis-aligned face and a wall — then the Pythagorean
    // ones, where the answer is exactly rational but (`5/3`) need not be an f64.
    for (v, want) in [
        ((1, 1), 1.0),
        ((4, 1), 0.5),
        ((25, 1), 0.2),
        ((9, 25), 5.0 / 3.0),
    ] {
        let r = Rat::new(v.0, v.1).unwrap();
        assert_eq!(
            inv_sqrt_exact(r).map(Rat::to_f64),
            Some(want),
            "1/sqrt{v:?}"
        );
        assert_eq!(inv_sqrt_f64(r), Some(want));
    }
    // …and a tilted one is honestly irrational, so it falls to the ladder.
    for v in [(2, 1), (3, 1), (5, 1), (1, 2)] {
        assert_eq!(inv_sqrt_exact(Rat::new(v.0, v.1).unwrap()), None, "{v:?}");
    }
    // There is no frame with a non-positive squared length; that is a broken premise.
    for v in [(0, 1), (-1, 1)] {
        let r = Rat::new(v.0, v.1).unwrap();
        assert_eq!(inv_sqrt_exact(r), None);
        assert_eq!(inv_sqrt_f64(r), None);
        assert!(inv_sqrt_bounded(r, 128).is_none());
    }
}

/// **Correctly rounded, which is what keeps debug and release identical** — the failure this
/// crate already paid for once in the trig path. Checked against a realization deep enough
/// that its own rounding cannot reach the 53rd bit.
#[test]
fn inv_sqrt_f64_lands_on_the_nearest_f64() {
    // Deltas, not absolute counts: the memo is a `thread_local` and the harness may have run
    // other tests on this thread first.
    let before = INV_SQRT_ESCALATED.with_borrow(|c| *c);
    for (num, den) in INV_SQRT_CASES {
        let v = Rat::new(num, den).unwrap();
        let deep = inv_sqrt_bounded(v, 1024).unwrap().value;
        let want = to_f64_exact(&deep).unwrap();
        assert_eq!(
            inv_sqrt_f64(v),
            Some(want),
            "1/sqrt({num}/{den}) is not the nearest f64"
        );
    }
    let after = INV_SQRT_ESCALATED.with_borrow(|c| *c);
    let (deeper, undecided) = (after.0 - before.0, after.1 - before.1);
    eprintln!("[cip] inverse-sqrt escalations: {deeper} deeper, {undecided} undecided");
    // ★ 128 bits sufficed for every case here, and the second rung is a proven cap rather
    // than a rung anything reaches: `1/√v` is never zero for a positive `v`, so its interval
    // cannot straddle zero the way `cos 90°`'s does and the rounding always decides.
    assert_eq!(
        (deeper, undecided),
        (0, 0),
        "a case escalated — the ladder's first rung no longer covers the population"
    );
}

/// **The reported realization error really does cover the error that is there** — checked
/// against a realization far deeper than the one the measurement itself uses (128 bits), so
/// the ruler and the thing being measured are not the same instrument.
#[test]
fn the_inverse_sqrt_realization_error_covers_the_error_that_is_there() {
    let mut worst = (i32::MIN, String::new());
    for (num, den) in INV_SQRT_CASES {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        let reported = inv_sqrt_error_of(v, f).unwrap();
        let HpBounded {
            value: deep,
            error: deep_rad,
        } = inv_sqrt_bounded(v, 1024).unwrap();
        let true_err = hp_err_exp(&deep, f);
        // A `BigFloat` is `m · 2^e` with `m ∈ [0.5, 1)`, so the exponent gives
        // `2^(e−1) ≤ |f − deep| < 2^e` — the *lower* end is what a bound has to clear. Using
        // `2^e` would demand the reported value exceed an over-estimate.
        //
        // ★ **And the reference is not exact either**, so its own radius comes off:
        // `|f − true| ≥ |f − deep| − rad`. Without that, `1/√(1/10⁶)` — which is exactly
        // `1000.0`, correctly reported as a zero error — fails against the 1024-bit
        // realization's `2⁻¹⁰¹³` residue, and the test would be measuring its own ruler.
        let ref_rad = deep_rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32));
        let floor = if true_err == i32::MIN {
            0.0
        } else {
            (2f64.powi(true_err - 1) - ref_rad).max(0.0)
        };
        assert!(
            reported >= floor,
            "1/sqrt({num}/{den}): reported {reported:e} is below the true error 2^{true_err}"
        );
        if true_err > worst.0 {
            worst = (true_err, format!("{num}/{den}"));
        }
    }
    // ★★★★ **Exact is not the same as dyadic, and only dyadic earns the zero.** An
    // axis-aligned frame lands on `1`/`½` and must report a literal `0.0` — that is what keeps
    // it on the exact predicate path. A Pythagorean one lands on `1/5` or `5/3`, which are
    // exactly rational and still not f64, so they must report the rounding that really
    // happened. Reporting `0` there (the shape the trig path can safely use) is the bug this
    // pair of assertions exists to keep out.
    // `1/9` is here rather than below because `1/√(1/9)` is `3` — an integer is dyadic too.
    for (num, den) in [(1i128, 1i128), (4, 1), (1, 4), (1, 64), (1, 9)] {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        assert_eq!(inv_sqrt_error_of(v, f), Some(0.0), "1/sqrt({num}/{den})");
    }
    // `(3,4,0)` gives `n·n = 25` and `1/|n| = 1/5`; `(3,4,0)` scaled gives `9/25` and `5/3`.
    for (num, den) in [(25i128, 1i128), (9, 25), (49, 1)] {
        let v = Rat::new(num, den).unwrap();
        let f = inv_sqrt_f64(v).unwrap();
        assert!(
            inv_sqrt_exact(v).is_some() && inv_sqrt_error_of(v, f).unwrap() > 0.0,
            "1/sqrt({num}/{den}) is exactly rational but not an f64 — it must be charged"
        );
    }
    eprintln!(
        "[cip] worst inverse-sqrt f64 realization error: 2^{} ({})",
        worst.0, worst.1
    );
}

/// ★★ **The naive f64 route is not good enough, measured** — `1.0 / v.to_f64().sqrt()` is
/// three roundings and lands on the wrong `f64` often enough to see. Without a case that
/// actually differs, the exact realization above would be a cost with no effect, and this
/// test would be indistinguishable from one that is not running.
#[test]
fn the_naive_f64_route_gets_the_last_bit_wrong() {
    let (mut n, mut differ) = (0usize, 0usize);
    for num in 1i128..400 {
        for den in 1i128..7 {
            let v = Rat::new(num, den).unwrap();
            let naive = 1.0 / v.to_f64().sqrt();
            n += 1;
            if inv_sqrt_f64(v) != Some(naive) {
                differ += 1;
            }
        }
    }
    eprintln!("[cip] naive 1/sqrt disagrees with the correctly rounded one: {differ}/{n}");
    assert!(
        differ > 0,
        "the naive route agreed everywhere on {n} cases — either the correctly rounded path \
             is not running, or this population cannot tell them apart"
    );
}

/// Error of a high-precision realization `a` against an exact f64 `b`, as a
/// power-of-two exponent (`i32::MIN` when exactly equal).
fn hp_err_exp(a: &BigFloat, b: f64) -> i32 {
    let err = a.sub(&BigFloat::from_f64(b, GT_PREC), GT_PREC, HP_RM);
    if err.is_zero() {
        i32::MIN
    } else {
        err.exponent().unwrap_or(i32::MIN)
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Rational `+` and `×` are commutative (exact, no drift).
    #[test]
    fn add_and_mul_commute(a in small_rat(), b in small_rat()) {
        prop_assert_eq!(a.checked_add(b), b.checked_add(a));
        prop_assert_eq!(a.checked_mul(b), b.checked_mul(a));
    }

    /// `a + b + c` associates regardless of grouping.
    #[test]
    fn add_associates(a in small_rat(), b in small_rat(), c in small_rat()) {
        let left = a.checked_add(b).and_then(|ab| ab.checked_add(c));
        let right = b.checked_add(c).and_then(|bc| a.checked_add(bc));
        prop_assert_eq!(left, right);
    }

    /// `from_deg` always normalizes into `[0, 360)`.
    #[test]
    fn from_deg_normalizes_into_range(deg in -3600i128..=3600) {
        let a = Angle::from_deg(Rat::from_int(deg)).unwrap();
        prop_assert!(a.deg() >= Rat::from_int(0));
        prop_assert!(a.deg() < Rat::from_int(360));
    }
}
