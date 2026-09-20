use super::*;
use crate::Angle;

/// An upper `f64` reading of a [`Mag`], for comparisons in tests only.
fn mag_f64(m: Mag) -> f64 {
    m.exp2().map_or(0.0, |e| 2f64.powi(e as i32))
}

/// |x − y| at the two values' own precision, as f64 — a test-side residual reader.
fn resid(x: &BigFloat, y: &BigFloat, prec: usize) -> f64 {
    let d = x.sub(y, prec, HP_RM).abs();
    if d.is_zero() {
        return 0.0;
    }
    d.exponent().map_or(f64::INFINITY, |e| 2f64.powi(e))
}

/// ★ The differential the two new interval operations are locked by: on **exact** inputs they
/// must agree with the exact routes that already exist, and on **interval** inputs the radius
/// must contain the truth — probed at the interval's endpoints, which is where a monotone
/// function's image is extreme (so the corners are the whole question, not a sample).
#[test]
fn the_interval_division_agrees_with_the_exact_route_and_contains_the_corners() {
    let prec = 192;
    let gt = 512;
    for (a, ra, b, rb) in [
        (3.75, 0.0, 1.5, 0.0),
        (-7.0, 1e-20, 0.3, 1e-21),
        (1e12, 1e-6, -2.5e-3, 1e-12),
        (0.1, 1e-30, 12345.678, 1e-18),
    ] {
        let av = HpBounded::new(big(a, prec), Mag::of(ra));
        let bv = HpBounded::new(big(b, prec), Mag::of(rb));
        let q = av.div(&bv, prec).expect("divisor clears zero");
        if ra == 0.0 && rb == 0.0 {
            // The exact route answers the same value, bit for bit — same BigFloat division.
            let qe = av.div_exact(&bv, prec).expect("exact divisor");
            assert!(
                q.value.sub(&qe.value, prec, HP_RM).is_zero(),
                "interval and exact division disagree on exact inputs"
            );
        }
        // Soundness at the corners: the true quotient for any (a', b') in the box lies
        // within value ± error. A quotient is monotone in each argument on a sign-constant
        // box, so the four corners bound the image.
        for (da, db) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
            let at = big(a + da * ra, gt);
            let bt = big(b + db * rb, gt);
            let truth = at.div(&bt, gt, HP_RM);
            let err = resid(&q.value, &truth, gt);
            assert!(
                err <= mag_f64(q.error) || err == 0.0,
                "corner ({da},{db}) of {a}±{ra} / {b}±{rb}: off by {err:e}, radius {:e}",
                mag_f64(q.error)
            );
        }
    }
}

#[test]
fn the_interval_inv_sqrt_agrees_with_the_exact_route_and_contains_the_endpoints() {
    let prec = 192;
    let gt = 512;
    for (x, r) in [(2.0, 0.0), (0.09, 1e-22), (1e20, 1.0), (5.0e-7, 1e-27)] {
        let xv = HpBounded::new(big(x, prec), Mag::of(r));
        let s = xv.inv_sqrt(prec).expect("bounded away from zero");
        if r == 0.0 {
            // Against the exact-rational route on a value both can state.
            let rat = Rat::from_decimal(x).expect("in the window");
            let exact = crate::inv_sqrt_bounded(rat, prec).expect("positive");
            let err = resid(&s.value, &exact.value, prec);
            assert!(
                err <= mag_f64(s.error) + mag_f64(exact.error),
                "exact input {x}: the two routes disagree by {err:e}"
            );
        }
        // 1/√x is monotone, so the endpoints are the extreme truths.
        for d in [1.0, -1.0] {
            let xt = big(x + d * r, gt);
            let truth = big(1.0, gt).div(&xt.sqrt(gt, HP_RM), gt, HP_RM);
            let err = resid(&s.value, &truth, gt);
            assert!(
                err <= mag_f64(s.error) || err == 0.0,
                "endpoint {d} of {x}±{r}: off by {err:e}, radius {:e}",
                mag_f64(s.error)
            );
        }
    }
}

/// The refusals: a divisor or radicand whose interval may reach zero is `None`, never a
/// guessed bound — and a negative radicand is refused outright.
#[test]
fn a_quantity_that_may_reach_zero_is_refused_not_bounded() {
    let prec = 128;
    let wide = HpBounded::new(big(1e-10, prec), Mag::of(1.0)); // straddles zero
    let a = HpBounded::exact(big(1.0, prec));
    assert!(a.div(&wide, prec).is_none(), "divisor may reach zero");
    assert!(wide.inv_sqrt(prec).is_none(), "radicand may reach zero");
    let neg = HpBounded::exact(big(-4.0, prec));
    assert!(neg.inv_sqrt(prec).is_none(), "a negative has no real root");
}

fn big(x: f64, prec: usize) -> BigFloat {
    BigFloat::from_f64(x, prec)
}

/// **Does the error `Angle::realization_error_of` reports actually bound the real one?**
///
/// `WitnessPoint::rotate_about` charges that number for the one input here without a rounding contract:
/// neither Rust nor any libm promises an accuracy for `f64::cos`. It used to charge a constant
/// measured once and written into a doc — sound only on machines like the one it was taken on,
/// which a kernel that ships to browsers cannot assume. Now it measures, so **a worse platform
/// reports a bigger number and the tolerance grows to match** rather than silently under-stating.
///
/// That moves the question. There is no budget left to overrun; what has to hold is that the
/// measurement is an *upper* bound. So this compares production's answer — taken at 128 bits —
/// against a realization 384 bits deeper.
///
/// **Not `#[ignore]`d, and that is the point.** The slow astro-float suites are skipped because
/// they are slow; this one exists so that *every build checks its own platform*, and a few
/// hundred angles is milliseconds.
///
/// ★★★ **It asks about the pair, because that is what the answer is about.** The values go in
/// as arguments — the same ones `rotate_about` writes into `coord` — rather than being re-derived
/// inside. Written the other way, this test asked whether the error of *some* realization bounded
/// the error of *another*, and it failed at 27°: a debug and a release build realize `sin 27°`
/// one ulp apart, and so do two call sites within one release build. That divergence is a real
/// and separate defect in the kernel's determinism story; what this test now checks is the
/// property `tol` actually needs, which is about the pair in hand.
///
/// ★ **The 90°-family is included, not skipped, and that is the load-bearing part.** There
/// `realization_error_of` claims *exactly zero*, and `rotate_about` reads that zero to contribute
/// no tolerance at all — which is what keeps a quadrantal origin rotation at tol 0 and keeps
/// axis-aligned models on the exact predicate path. A claimed zero that was not really zero
/// would not fail loudly; it would quietly move those models. Here it has to be real.
///
/// ★★ **It also still reports the worst error, though nothing depends on it.** Charging a
/// measured value means a degraded platform no longer fails anything — it just escalates more
/// and runs slower, invisibly. This line is what keeps that visible.
///
/// What it does not do: it checks one input. `nacre-judge`'s `tol_bounds_error_over_random_chains`
/// checks the conclusion — that `tol` bounds the error, second-order terms and all.
#[test]
fn the_measured_trig_error_bounds_the_real_one() {
    const GT: usize = 512; // the reference it is checked against
    let mut st = 0x9E37_79B9_7F4A_7C15u64;
    let mut rng = |lo: i128, hi: i128| {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + ((st >> 33) as i128) % (hi - lo).max(1)
    };
    let mut angles: Vec<(i128, i128)> = (0..360).map(|d| (d, 1)).collect();
    for _ in 0..200 {
        angles.push((rng(0, 360_000), rng(1, 9973)));
        // ★ Near an axis, where the realized cos or sin is close to zero — the region that
        // refuted the tangential tol model, and where a relative error reads worst.
        angles.push((90 * rng(0, 4) * 1000 + rng(-3, 4), 1000));
        // ★ Past `2⁵³`, where `Rat::to_f64` leaves its fast path. The radian argument is what
        // dominates the error, so an angle whose conversion takes the other route has to be here.
        angles.push((rng(1, 1 << 62) * 360, rng(1 << 53, 1 << 62)));
    }
    // ★ **Reported to the nearest ε, which needs a comparison rather than an octave reader.**
    // `2^exponent` is an upper *octave*, up to 2× above the true magnitude — so a figure read
    // from it blurs exactly the range that matters. The smallest integer `k` with `x < k·ε` is
    // one comparison per candidate and says the thing plainly.
    //
    // ★★★ **The range starts at 1, and it must**: `0·ε` is `0.0`, and astro-float's `cmp`
    // answers `Some(-1)` in *both* directions against a zero — so `x.cmp(&zero)` reports "less
    // than" for a positive `x`, every inexact realization reads `k = 0`, and `worst_k` is then
    // only ever raised by the *exact* ones, which read the floor. That prints a headline of
    // "worst < 1ε" for a platform whose real figure is several times that. A `k` of 0 would be
    // meaningless anyway — nothing is below `0·ε` — which is what makes the widening look free.
    let eps_multiple = |x: &BigFloat| -> usize {
        (1..=64)
            .find(|k| {
                // astro-float's `cmp` yields a sign as `Option<i128>`, not an `Ordering`.
                x.abs()
                    .cmp(&big(*k as f64 * f64::EPSILON, GT))
                    .is_some_and(|sign| sign < 0)
            })
            .unwrap_or(usize::MAX)
    };
    // ★★★ **The reading is pinned before it is trusted.** Three separate figures in this
    // kernel's error work turned out to be artefacts of how they were printed rather than of
    // what was measured, and each time the tell was a number that would not move. A reading
    // that cannot be shown to respond to a known input is not evidence.
    assert_eq!(eps_multiple(&big(0.0, GT)), 1, "zero reads as the floor");
    assert_eq!(eps_multiple(&big(3.5 * f64::EPSILON, GT)), 4);
    assert_eq!(eps_multiple(&big(-9.0 * f64::EPSILON, GT)), 10);
    assert_eq!(eps_multiple(&big(64.0 * f64::EPSILON, GT)), usize::MAX);
    let (mut worst_k, mut worst_at) = (0usize, (0i128, 1i128));
    let (mut checked, mut exact_seen) = (0usize, 0usize);
    for (n, d) in angles {
        let Some(a) = Rat::new(n, d).and_then(Angle::from_deg) else {
            continue; // outside Rat's range — `from_deg` declines, and so does the kernel
        };
        let (c, s) = a.cos_sin_f64();
        let (dc, ds) = a.realization_error_of(c, s);
        if dc == 0.0 && ds == 0.0 {
            exact_seen += 1;
        }
        let (hc, hs) = a.cos_sin_bounded(GT);
        for (which, f, h, reported) in [("cos", c, &hc, dc), ("sin", s, &hs, ds)] {
            // |f64 − true| ≤ |f64 − deep midpoint| + that realization's own radius.
            let truth = big(f, GT).sub(&h.value, GT, HP_RM).abs().add(
                &big(h.error.exp2().map_or(0.0, |e| 2f64.powi(e as i32)), GT),
                GT,
                HP_RM,
            );
            assert!(
                truth.cmp(&big(reported, GT)).is_some_and(|sign| sign <= 0),
                "realization_error_of under-states at {n}/{d} deg [{which}]: \
                     f = {f:e}, reported {reported:e}, true error 2^{:?}",
                truth.exponent(),
            );
            // ★ The worst is tracked over the *inexact* realizations only. An exact one has
            // error zero, which reads as the floor `k = 1` — and since the corpus opens at 0°,
            // including them let the family claim the headline and named an angle with no
            // error at all as this platform's worst case.
            let k = eps_multiple(&truth);
            if reported != 0.0 && k > worst_k {
                (worst_k, worst_at) = (k, (n, d));
            }
            checked += 1;
        }
    }
    assert!(checked > 700, "corpus shrank to {checked} realizations");
    assert!(
        exact_seen > 0,
        "no quadrantal angle reached the zero-error claim — the strongest case went untested"
    );
    eprintln!(
        "[trig realization] worst error < {worst_k}ε (at {}/{} deg), {exact_seen} exact, \
             over {checked} realizations",
        worst_at.0, worst_at.1,
    );
}

/// The point of the type: a value that is only rounding residue must not report a sign, no
/// matter how confidently its `value` is nonzero.
#[test]
fn a_residue_left_by_cancellation_reports_no_sign() {
    let prec = 200;
    // Two large values differing by far less than the radius they carry. The separation has
    // to stay inside `prec` bits of the operands, or the subtraction is exactly zero and the
    // test proves nothing.
    let error = Mag::of(1.0e-10);
    let a = HpBounded::new(big(1.0e30, prec), error);
    let b = HpBounded::new(
        big(1.0e30, prec).sub(&big(1.0e-20, prec), prec, HP_RM),
        error,
    );
    let d = a.sub(&b, prec);
    assert!(
        !d.value.is_zero(),
        "the midpoint is nonzero — that is the trap"
    );
    assert_eq!(
        d.sign(),
        None,
        "a difference under its own radius claimed a sign"
    );
}

/// …and a value genuinely clear of its radius still decides, so the type is not merely
/// refusing to answer.
#[test]
fn a_value_clear_of_its_radius_still_decides() {
    let prec = 200;
    let a = HpBounded::new(big(3.0, prec), Mag::pow2(-100));
    let b = HpBounded::new(big(2.0, prec), Mag::pow2(-100));
    assert_eq!(a.sub(&b, prec).sign(), Some(true));
    assert_eq!(b.sub(&a, prec).sign(), Some(false));
    assert_eq!(a.mul(&b, prec).sign(), Some(true));
}

/// The radius survives a rung deep enough to flush an `f64` radius to zero — the reason
/// [`Mag`] exists. With an `f64` radius this product would report a confident sign.
#[test]
fn a_deep_rung_does_not_lose_the_radius() {
    let prec = 2048;
    let tiny = HpBounded::new(big(1.0, prec), Mag::pow2(-(prec as i64)));
    let p = tiny.mul(&tiny, prec);
    assert!(!p.error.is_zero(), "the radius vanished at {prec} bits");
    // A difference of exactly that size is therefore undecided, not positive.
    let q = HpBounded::new(
        big(1.0, prec).add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
        Mag::pow2(-(prec as i64) + 4),
    );
    let r = HpBounded::new(big(2.0, prec), Mag::pow2(-(prec as i64) + 4));
    assert_eq!(q.sub(&r, prec).sign(), None);
}

/// An exact input contributes no radius of its own, so the only radius an exact computation
/// carries is the half-ulp charged per operation.
///
/// That charge is unconditional: astro-float does not report whether an operation was exact,
/// so a subtraction of two dyadics — exact in fact — is still charged. Conservative in the
/// safe direction, and negligible in size, but it does mean **this path never proves zero**.
/// Judgements that need a proved zero take the exact route above this kernel instead.
#[test]
fn exact_inputs_keep_a_radius_only_from_rounding() {
    let prec = 200;
    let a = HpBounded::exact(big(0.5, prec));
    let b = HpBounded::exact(big(0.25, prec));
    assert!(
        a.error.is_zero() && b.error.is_zero(),
        "an exact input carried a radius"
    );
    let d = a.sub(&b, prec);
    assert!(
        d.error.lt(Mag::pow2(-190)),
        "an exact subtraction picked up more than a half-ulp: 2^{:?}",
        d.error.exp2()
    );
    assert_eq!(d.sign(), Some(true));
    // Exactly equal inputs are the case that cannot be settled here: the difference is zero
    // and zero has no sign to read.
    assert_eq!(a.sub(&a, prec).sign(), None);
}

/// The one place the zero guard shows: rounding an exact zero costs nothing, so a product or
/// sum that lands exactly on zero from exact inputs reports the radius it earned — none.
/// [`Mag::above`] answers `ZERO` there where a guard-less exponent reader would not.
#[test]
fn an_exact_zero_is_charged_no_rounding() {
    let prec = 128;
    let zero = HpBounded::exact(big(0.0, prec));
    let one = HpBounded::exact(big(1.0, prec));
    assert!(
        zero.mul(&one, prec).error.is_zero(),
        "0·1 charged a rounding"
    );
    assert!(
        one.sub(&one, prec).error.is_zero(),
        "1−1 charged a rounding"
    );
    assert_eq!(HpBounded::round_off(&big(0.0, prec), prec), Mag::ZERO);
}
