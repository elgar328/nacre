use super::*;
impl Angle {
    /// Normalize `deg` into `[0, 360)` (exact), in one division: the numerator's
    /// remainder modulo one turn (`360·denom`), so any number of turns reduces in
    /// one step. The result keeps no turn count — this is a direction. `None` when
    /// `360·denom` overflows `i128`.
    pub fn from_deg(deg: Rat) -> Option<Self> {
        // Exact reduction mod 360, in one division. Subtracting a turn at a time is the same
        // arithmetic but costs one iteration per turn, so an input like `2⁶⁰` degrees — a value a
        // script can produce without meaning anything unusual by it — does not return.
        let (n, d) = (deg.numer(), deg.denom()); // `d > 0` after reduction
        let full = 360i128.checked_mul(d)?; // one turn, in units of `1/d`
        let rem = n.rem_euclid(full); // `0 ≤ rem < full`, so the result is in `[0, 360)`
        Rat::new(rem, d).map(Angle)
    }

    /// Turn by `delta_deg` (exact) and renormalize into `[0, 360)`.
    pub fn checked_add(self, delta_deg: Rat) -> Option<Self> {
        Self::from_deg(self.0.checked_add(delta_deg)?)
    }

    /// The normalized degree value (exact, in `[0, 360)`).
    pub fn deg(self) -> Rat {
        self.0
    }

    /// `(cos, sin)` realized in arbitrary precision at `prec` bits — the judgment path.
    /// astro-float, not twofloat: twofloat's trig is only f64-level near zero-crossings. The depth
    /// is the caller's: a judgement's precision is a property of the model it judges, not of this
    /// crate. (numer/denom pass through f64, exact
    /// for the small values used here; a general large-rational path would build from a string.)
    pub fn cos_sin_at(self, prec: usize) -> (BigFloat, BigFloat) {
        let (c, s) = self.cos_sin_bounded(prec);
        (c.value, s.value)
    }

    /// `(cos, sin)` at `prec` bits **with an upper bound on how far each may be from the true
    /// value** — the seed every error radius in the judgment path grows from.
    ///
    /// The bound is *derived*, term by term, not chosen:
    ///
    /// - `numer`/`denom` enter as `i128`, exactly. Routing them through `f64` would cost a
    ///   relative `2⁻⁵³` past 2⁵³ that no working precision recovers, so it is not done.
    /// - `n/d`, `·π`, `/180` are three round-to-nearest operations at `prec` bits, each a relative
    ///   `2⁻ᵖʳᵉᶜ`, and `π` itself carries one more.
    /// - The argument's absolute error `δθ = |θ|·ρ` passes through the **derivative**:
    ///   `d cos = −sin·dθ` and `d sin = cos·dθ`. Slope 1 would also be sound, but near a zero
    ///   crossing the true slope is what keeps a tiny result from being swamped by its own bound.
    /// - astro-float's `cos`/`sin` run **Ziv's loop** — `cos_series` at a working precision, then
    ///   `try_set_precision(p, rm, p_wrk)`, retrying with more bits when the value sits too close
    ///   to a rounding boundary to decide. That is the standard construction for a *correctly
    ///   rounded* transcendental, so the realization adds at most a half-ulp — `|result| · 2⁻ᵖʳᵉᶜ`,
    ///   *relative* to the value, which is why the two functions get separate bounds. The crate
    ///   does not document this, so `the_trig_bound_holds_against_a_far_deeper_realization` checks
    ///   it rather than trusting it.
    ///
    /// **Memoized by `(self, prec)`** in [`TRIG`] — the value below is a pure function of those
    /// two, and a rotation's realization asks for the same angle once per point. See [`TRIG`] for
    /// why a process-wide memo is sound here and a handle-keyed one would not be.
    ///
    /// Returns `(cos, sin)`, each with the radius it carries.
    pub fn cos_sin_bounded(self, prec: usize) -> (HpBounded, HpBounded) {
        if let Some(hit) = TRIG.with_borrow(|t| t.get(&(self, prec)).cloned()) {
            return hit;
        }
        let out = self.realize_cos_sin(prec);
        TRIG.with_borrow_mut(|t| t.insert((self, prec), out.clone()));
        out
    }

    /// [`cos_sin_bounded`](Self::cos_sin_bounded) without the memo — the evaluation itself.
    ///
    /// **Separate so no `TRIG` borrow is held across it.** `HP_CONSTS` is borrowed for the whole
    /// realization and the trig calls are the slow part; nesting the memo's borrow around that is
    /// how a re-entrant call would panic rather than merely be slow.
    fn realize_cos_sin(self, prec: usize) -> (HpBounded, HpBounded) {
        HP_CONSTS.with_borrow_mut(|cc| {
            let pi = cc.pi(prec, HP_RM);
            let d180 = BigFloat::from_f64(180.0, prec);
            // As `i128`, not through `f64`: past 2⁵³ the conversion would cost a relative
            // `2⁻⁵³` that no working precision recovers, because the loss happens before
            // astro-float sees the value.
            let ip = prec.max(128);
            let n = BigFloat::from_i128(self.0.numer(), ip);
            let d = BigFloat::from_i128(self.0.denom(), ip);
            let rad = n
                .div(&d, prec, HP_RM)
                .mul(&pi, prec, HP_RM)
                .div(&d180, prec, HP_RM);
            let u = Mag::pow2(-(prec as i64));
            // Relative error of the argument: the two `i128 → f64` conversions, then four
            // rounded high-precision operations (the division, the product, the division, and π).
            // Four rounded operations build the argument (the division, the product, the
            // division, and π itself). The integers contribute nothing — they go in exactly.
            let rel = u.times(Mag::of(4.0));
            // `|θ|` in radians, over-estimated from its exponent (`|x| < 2^exponent`).
            let d_theta = Mag::above(&rad).times(rel);
            let (c, s) = (rad.cos(prec, HP_RM, cc), rad.sin(prec, HP_RM, cc));
            // `|x| < 2^exponent` — the slope of the *other* function, and the scale of the
            // half-ulp of this one.
            let (uc, us) = (Mag::above(&c), Mag::above(&s));
            let err_cos = us.times(d_theta).plus(uc.times(u));
            let err_sin = uc.times(d_theta).plus(us.times(u));
            (HpBounded::new(c, err_cos), HpBounded::new(s, err_sin))
        })
    }

    /// Exact `(cos, sin)` as rationals — `Some` only for the quadrantal angles
    /// (0/90/180/270°), the sole angles where *both* are rational (`{0, ±1}`, by
    /// Niven); `None` otherwise, so the caller falls to the f64/dd realization.
    ///
    /// The values `0`/`±1` are exact in f64 too, so this returns `Rat` not for
    /// representability but to keep a rotated *point* exact: `x·cos − y·sin` must
    /// stay in rational arithmetic, and multiplying a rational coordinate by an
    /// f64 (even an exact `0.0`/`1.0`) would drop the point into f64 and lose the
    /// very exactness this path exists for. So a 90°-family rotation of a rational
    /// point stays tol 0.
    pub fn try_exact_cos_sin(self) -> Option<(Rat, Rat)> {
        let zero = Rat::from_int(0);
        let one = Rat::from_int(1);
        let neg_one = Rat::from_int(-1);
        if self.0 == Rat::from_int(0) {
            Some((one, zero))
        } else if self.0 == Rat::from_int(90) {
            Some((zero, one))
        } else if self.0 == Rat::from_int(180) {
            Some((neg_one, zero))
        } else if self.0 == Rat::from_int(270) {
            Some((zero, neg_one))
        } else {
            None
        }
    }

    /// `(cos, sin)` realized in f64 — **exact** (`0.0`/`±1.0`) for the 90°-family, and everywhere
    /// else the arbitrary-precision value **correctly rounded**. The single source of truth for
    /// realizing a rotation angle into f64: every path that turns a point or direction by an angle
    /// goes through here.
    ///
    /// ★★★ **No libm.** This used to be `(deg.to_f64() * PI / 180.0).cos()`, and the error of that
    /// is not something anyone contracts: neither Rust nor any platform promises an accuracy for
    /// `f64::cos`, measured here at `< 5ε`. Worse, it was not a function of the angle — `sin 27°`
    /// came out one ulp apart between a debug and a release build, and between two call sites
    /// *within one release build*, because LLVM evaluates a visible constant angle at compile time
    /// and its answer differs from the runtime library's. Rounding the high-precision realization
    /// instead makes the result **unique**: the same bits on every platform, profile and call site.
    ///
    /// ★★ **The 90°-family branch is not an optimization, it is what makes this terminate.**
    /// `cos 90°` is exactly `0`, so its interval straddles zero and the two ends never round to the
    /// same f64 no matter how deep the realization goes — [`round_to_f64`] would answer `None`
    /// forever. Niven's theorem says the only rational values are `{0, ±1/2, ±1}`, and only the
    /// zeros have this problem; they are exactly the quadrantal ones caught here. (`±1/2` at 60°
    /// and friends is exactly representable and comes out of the general path just fine.)
    ///
    /// **The realization is memoized per angle** ([`F64_TRIG`]), because one costs ~32µs against
    /// libm's 22ns and `Isometry::apply_point` runs per vertex. The memo cannot change an answer —
    /// a correctly rounded value is unique — so it is a cost question only.
    ///
    /// A caller accounting for the realization's *error* hands what it got back to
    /// [`realization_error_of`](Self::realization_error_of), which returns zero for exactly this
    /// family — so nothing here has to report which branch ran.
    pub fn cos_sin_f64(self) -> (f64, f64) {
        if let Some((cr, sr)) = self.try_exact_cos_sin() {
            return (cr.to_f64(), sr.to_f64());
        }
        if let Some(hit) = F64_TRIG.with_borrow(|m| m.get(&self).copied()) {
            return hit;
        }
        // Outside the borrow: the realization below takes `TRIG`'s and `HP_CONSTS`' in turn.
        let out = self.realize_rounded_f64();
        F64_TRIG.with_borrow_mut(|m| m.insert(self, out));
        out
    }

    /// [`cos_sin_f64`](Self::cos_sin_f64)'s general branch without the memo.
    ///
    /// **The ladder is `TRIAL_PREC` then twice that, and both rungs are derived rather than tried.**
    /// The realization's error is dominated by the degrees→radians conversion, not by the cosine:
    /// `≈ 30 · 2⁻ᵖʳᵉᶜ`. Against a result of magnitude `2^e` that has to clear a half-ulp of `2^(e-54)`,
    /// so `prec > 54 - e + 5`.
    ///
    /// - **128** covers every `|cos| > 2⁻⁶⁹`, which is every angle a model has ever held. It is also
    ///   `nacre_judge`'s trial precision, so a model that goes on to be judged **shares this exact
    ///   realization** rather than paying for a second one at a different depth. That sharing is
    ///   why the ladder does not start lower: 64 bits would satisfy the inequality and measured no
    ///   cheaper (31.9µs against 32.3µs — the cost is setup, not bit count), but it would be a
    ///   different `TRIG` key and so pure duplication for anything judged.
    /// - **256** is the proven cap. An `Angle` holds `Ratio<i128>`, so a normalized angle cannot
    ///   come closer to 90° than `1/denominator ≥ 5.9e-39` degrees; `|cos|` is therefore never
    ///   below `~2⁻¹³³`, which needs `prec > 192`.
    ///
    /// ★ **If even 256 leaves it undecided the answer is still returned, not a panic.** The value is
    /// then *faithfully* rounded (within an ulp) instead of correctly rounded, which stays sound
    /// because [`realization_error_of`](Self::realization_error_of) measures the error that is
    /// actually there and the tolerance grows to match — and it stays deterministic, because a
    /// 256-bit midpoint is. It is counted so that "can't happen" does not quietly become "happens".
    fn realize_rounded_f64(self) -> (f64, f64) {
        for (i, prec) in [128usize, 256].into_iter().enumerate() {
            let (c, s) = self.cos_sin_bounded(prec);
            if let (Some(cf), Some(sf)) = (
                round_to_f64(&c.value, c.error, prec),
                round_to_f64(&s.value, s.error, prec),
            ) {
                if i > 0 {
                    ROUND_ESCALATED.with_borrow_mut(|(e, _)| *e += 1);
                }
                return (cf, sf);
            }
        }
        ROUND_ESCALATED.with_borrow_mut(|(_, f)| *f += 1);
        let (c, s) = self.cos_sin_bounded(256);
        (
            to_f64_exact(&c.value).unwrap_or(f64::NAN),
            to_f64_exact(&s.value).unwrap_or(f64::NAN),
        )
    }

    /// **How far the `(cos, sin)` the caller was handed sits from the true ones** — measured
    /// against an arbitrary-precision realization, not assumed from a constant.
    ///
    /// `|f64 − true| ≤ |f64 − hp midpoint| + hp's own radius`, which is the ruler this kernel
    /// already uses for a rational's realization in `WitnessPoint`'s `translate`, `mirror` and pivot terms.
    /// Memoized in [`F64_ERR`], keyed by the angle *and the pair* — see there for why the pair.
    ///
    /// ★★★ **The caller passes the values in rather than letting this re-realize them**, and that
    /// is the whole soundness argument. The consumer is `WitnessPoint::rotate_about`, whose `tol` must bound
    /// the error in the `coord` it just wrote from *its* `cos_sin_f64()` result. An error measured
    /// against a second, independent realization would bound a number nobody stored — and those two
    /// realizations are measured to differ (see [`F64_ERR`]). Taking `c` and `s` as arguments makes
    /// "the error describes the value that was used" hold by construction instead of by hope.
    ///
    /// ★★ **This is what lets the error accounting stop guessing.** `f64::cos` has no accuracy
    /// contract — neither Rust nor any libm promises one — so a measured-once constant with margin
    /// would be sound only on platforms like the one it was taken on, which a kernel that ships to
    /// browsers cannot know. Measuring instead means a worse libm
    /// simply reports a larger error and the tolerance grows to match: **the kernel adapts rather
    /// than assumes.**
    ///
    /// ★ **Exactly zero for the 90°-family**, because `cos_sin_f64` returns `0.0`/`±1.0` there and
    /// those *are* the true values. Callers rely on that zero: a rotation whose realization carries
    /// no error also performs no rounding downstream (`u·(±1)` and `u·0` are exact), which is what
    /// keeps a quadrantal origin rotation at tol 0.
    ///
    /// ★ **An `f64`, not a [`Mag`].** `Mag` exists because a deep ladder's `2⁻ᵖʳᵉᶜ` underflows
    /// `f64` to zero and a zero radius claims exactness; this quantity is always ε-scale, so that
    /// hazard is absent — and the consumer is `WitnessPoint::tol`, which is `f64`.
    ///
    /// The reading is at octave granularity (`2^exponent` of the residual), so it can sit up to 2×
    /// above the true error. Conservative in the sound direction, and still a measurement.
    pub fn realization_error_of(self, c: f64, s: f64) -> (f64, f64) {
        let key = (self, c.to_bits(), s.to_bits());
        if let Some(hit) = F64_ERR.with_borrow(|m| m.get(&key).copied()) {
            return hit;
        }
        // Outside the borrow: the computation below takes `TRIG`'s and `HP_CONSTS`' in turn.
        let out = self.measure_realization_error(c, s);
        F64_ERR.with_borrow_mut(|m| m.insert(key, out));
        out
    }

    /// [`realization_error_of`](Self::realization_error_of) without the memo — the measurement.
    fn measure_realization_error(self, c: f64, s: f64) -> (f64, f64) {
        const P: usize = 128; // the hp radius is then ~2⁻¹²⁸ against an ε-scale quantity
        // ★ The zero is claimed of *these* values, not of the angle: the exactness that callers
        // depend on is "the pair I am holding is the true cos/sin", and only comparing the pair
        // says that. An angle in the family whose caller somehow realized it the general way is
        // then measured like any other rather than being handed a zero it has not earned.
        if self
            .try_exact_cos_sin()
            .is_some_and(|(cr, sr)| cr.to_f64() == c && sr.to_f64() == s)
        {
            return (0.0, 0.0);
        }
        let (hc, hs) = self.cos_sin_bounded(P);
        let (rc, rs) = (hc.error, hs.error);
        let (hc, hs) = (hc.value, hs.value);
        let gap = |f: f64, h: &BigFloat, rad: Mag| {
            let diff = BigFloat::from_f64(f, P).sub(h, P, HP_RM);
            let mag = if diff.is_zero() {
                0.0
            } else {
                2f64.powi(diff.exponent().unwrap_or(0))
            };
            // ★ **The sum is rounded *up*, and that is not pedantry.** `mag` is an octave bound, so
            // it usually sits well above the truth — but when `|diff|` is itself a power of two the
            // slack is exactly zero, and then `mag + rad` rounds back down to `mag` in f64 and the
            // "bound" is short by the radius. `sin 30°` is that case: it misses 0.5 by exactly
            // `2⁻⁵⁴`, and the ground-truth test caught the missing ulp the day this was written.
            (mag + rad.exp2().map_or(0.0, |e| 2f64.powi(e as i32))) * (1.0 + 2.0 * f64::EPSILON)
        };
        (gap(c, &hc, rc), gap(s, &hs, rs))
    }
}
