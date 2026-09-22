use super::*;
/// **The standard of proof a judgement is held to: how deep to realize, and how close counts as
/// one thing.**
///
/// The threshold here is a **length in the model's own units**, never a bit count. "256 bits"
/// means a resolution of `1e-76` for a solid turned once and `1e+15` for one turned three hundred
/// times — the same setting meaning entirely different things per model, which is not something
/// anyone can reason about. Bits are the implementation detail the kernel computes per model
/// ([`judge_precision`]); the length is the physics.
///
/// `coincidence` is not a tolerance in the usual sense, and deliberately not called one: a global
/// tol says *"anything closer than this, snap together"* (merging without knowing), while this
/// says *"a coincidence must be **proved** to be closer than this"*. Nothing is merged on
/// ignorance; a judgement that cannot prove it climbs, and then says so.
#[derive(Clone, Copy, Debug)]
pub struct Standard {
    /// The precision the escalation realizes definitions at — chosen for this model, uniform
    /// across the operation so [`WitnessPoint`]'s realization cache stays warm.
    pub prec: usize,
    /// Two things **proved** to lie within this distance of each other are one thing.
    pub coincidence: Mag,
    /// The model's size — what turns the length limit into an angle for the one judgement whose
    /// question is about directions ([`dir_sign_judge`]).
    pub scale: Mag,
    /// The most bits an escalation may ask for. Past it the judgement is [`Decision::Exhausted`]:
    /// answerable in principle, too expensive in practice, and said out loud rather than guessed.
    pub cap: usize,
}

/// **What a judgement established** — the answer *and* what backs it.
///
/// `Orient::Zero` alone cannot say whether a zero was proved or assumed, which is how a kernel
/// ends up quietly guessing. These four cases are the honest partition:
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Decision {
    /// A **proved** sign. `Zero` appears here only from a path that can prove one — an exact
    /// predicate, or a rotation that cancelled — never from a determinant that merely failed to
    /// separate from zero.
    Sign(Orient),
    /// The two things were shown to lie within `within` of each other, and that is at or below
    /// the coincidence limit. They are treated as one, and `within` is the evidence for it.
    Coincident { within: Mag },
    /// Still straddling zero at `at` bits, the cap — and the separation it might stand for is
    /// *larger* than the coincidence limit, so calling it a coincidence would be a guess. More
    /// precision would decide it; this judgement has outgrown the budget.
    ///
    /// **`within` is what it did establish**, and it is the whole point of reporting this rather
    /// than a bare "undecided": the truth lies in `±within`, above the limit but bounded. A
    /// judgement stopped at `1e-30` and one stopped at `1e-3` are the same variant and *not* the
    /// same news — the first is far below anything the `f64` output can carry, the second is a
    /// gap somebody would see. `None` when even that bound could not be formed: the cofactor was
    /// still unresolved at the cap, so there is no distance to quote at all.
    Exhausted { at: usize, within: Option<Mag> },
    /// The quantity that turns this determinant into a distance cannot be bounded away from zero
    /// — a degenerate witness triangle, or three planes with no well-defined meeting point. **A
    /// different cause from [`Self::Exhausted`], and the distinction matters**: more precision
    /// never helps here, because the question has no metric answer to sharpen.
    Degenerate,
}

impl Decision {
    /// The sign the geometry consumes, with every inconclusive outcome collapsing to `Zero`.
    ///
    /// A proved coincidence and an exhausted judgement are the *same instruction* to the
    /// arrangement — "treat these as equal" — and differ only in what can be said about it
    /// afterwards. Keeping the difference out of the control flow is what lets the reporting
    /// channel be added without touching a single geometric decision.
    pub fn orient(self) -> Orient {
        match self {
            Decision::Sign(o) => o,
            _ => Orient::Zero,
        }
    }
}

/// **The separation an undecided determinant stands for — or why it could not be formed.**
///
/// Turning a determinant into a distance means dividing by a cofactor, and that division needs the
/// cofactor kept away from zero. It can fail two ways, and they are **not the same answer**:
#[derive(Clone, Copy, Debug)]
pub(super) enum Gap {
    /// An upper bound on the separation, in the unit the judgement's limit is in.
    Of(Mag),
    /// A cofactor is nonzero but has not been separated from its own error radius yet. `short`
    /// bits more would separate it — **so this climbs**, exactly like a gap that is merely too
    /// wide. Collapsing it into the case below is what would make the kernel quietly merge two
    /// things it never looked at closely enough.
    Unresolved { short: usize },
    /// A cofactor came out **exactly zero**, so there is no distance to sharpen: a witness
    /// triangle with no plane, or three planes with no meeting point. Depth does not change it —
    /// measured directly on the case this distinction was found in, where the determinant was
    /// still bit-exactly zero 8192 bits deeper.
    Vanished,
}

/// A denominator's lower bound, or which of the two failures it is.
///
/// The mantissa is deliberately not read: `lb` is an exponent bound, so `short` is generous by up
/// to a bit — in the direction that costs a word, not correctness.
pub(super) fn denom_lo(v: &HpBounded) -> Result<Mag, Gap> {
    let Some(lo) = Mag::below(&v.value) else {
        return Err(Gap::Vanished); // the midpoint is exactly zero
    };
    match lo.minus(v.error) {
        Some(b) => Ok(b),
        // Nonzero, but the radius swallows it: the shortfall is how far the radius has to fall.
        None => Err(Gap::Unresolved {
            short: match (v.error.exp2(), lo.exp2()) {
                (Some(r), Some(m)) => (r - m + 1).max(1) as usize,
                _ => 1,
            },
        }),
    }
}

pub(super) fn escalate(
    j: Standard,
    limit: Mag,
    mut attempt: impl FnMut(usize) -> Result<Orient, Gap>,
) -> Decision {
    climb_census::CLIMBS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    // A zero precision means the context never got stamped: astro-float would be asked for a
    // realization with no bits, and every judgement would come back exhausted. That is a wiring
    // mistake, not a geometry one, so it fails loudly here rather than quietly answering `Zero`.
    debug_assert!(
        j.prec > 0,
        "escalate at zero precision — unstamped Standard"
    );
    // ★ Charged on **every** way out, not just the resolved one: most climbs end in `Coincident`
    // or a proved zero, and counting only the `Ok` path reported a mean of exactly 0 bits — a
    // number that looked like "no cost" and was really "no instrument".
    let done = |p: usize, d: Decision| -> Decision {
        climb_census::BITS.fetch_add(p as u64, std::sync::atomic::Ordering::Relaxed);
        d
    };
    let mut prec = j.prec;
    loop {
        let gap = match attempt(prec) {
            Ok(o) => return done(prec, Decision::Sign(o)),
            Err(g) => g,
        };
        // How many bits short this judgement is, and the bound it managed — the second is what
        // gets reported if the budget runs out. From a gap the shortfall is `log₂(gap / limit)`;
        // from an unresolved cofactor the cofactor itself named it, and there is no bound to
        // quote. Either way the number is *computed*, never guessed at by doubling.
        let (short, within) = match gap {
            Gap::Vanished => return done(prec, Decision::Degenerate),
            Gap::Unresolved { short } => (short, None),
            // **A zero gap is a proof, not a near miss.** The radius bounds how far the computed
            // value is from the true one, so a zero radius says the midpoint *is* the value — and
            // the sign came back undecided only because that midpoint is zero. Determinants reach
            // this honestly: a row that is exactly zero (a plane perpendicular to the rotation
            // axis keeps its coordinate exactly) multiplies every error term to nothing.
            Gap::Of(g) if g.is_zero() => return done(prec, Decision::Sign(Orient::Zero)),
            Gap::Of(g) if !limit.lt(g) => return done(prec, Decision::Coincident { within: g }),
            Gap::Of(g) => match (g.exp2(), limit.exp2()) {
                (Some(gx), Some(lx)) => ((gx - lx).max(1) as usize, Some(g)),
                // A limit with no exponent (a zero bound) is a target no depth reaches.
                _ => {
                    climb_census::EXHAUSTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    return done(
                        prec,
                        Decision::Exhausted {
                            at: prec,
                            within: Some(g),
                        },
                    );
                }
            },
        };
        if prec >= j.cap {
            climb_census::EXHAUSTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return done(prec, Decision::Exhausted { at: prec, within });
        }
        // Up to a whole word, and never a standstill: a jump that rounded back to `prec` would
        // spin here forever.
        prec = (prec + short)
            .div_ceil(WORD)
            .max(prec / WORD + 1)
            .saturating_mul(WORD)
            .min(j.cap);
    }
}

/// The precision at which a realization's error is *measured*, before the real one is chosen.
///
/// Nothing is judged here — this only has to be deep enough that the radius it produces is a
/// meaningful reading of `C` (see [`judge_precision`]), and cheap.
const TRIAL_PREC: usize = 128;

/// ★ **And deep enough that an exactly-representable point reads exactly zero.**
///
/// A second requirement, once a caller is allowed to *skip* [`trial_bound`] for a definition it
/// knows is exact. An f64-representable base (`Rat::try_from_f64` = `mantissa · 2^exp`, the case
/// `at_nearest` states with tol 0) has a power-of-two denominator and a numerator of at most
/// `Rat`'s own 127 bits — and
/// `HpBounded::of_rat` returns an exact interval exactly when both hold *at this precision*. Below
/// 127 a
/// large exact coordinate would start carrying a bound again, and a caller that skipped the call on
/// the strength of the zero would read a precision the model had not earned.
/// See `an_exact_point_demands_no_precision`.
const _: () = assert!(TRIAL_PREC >= 127);

/// Bits per word: astro-float allocates whole words, so asking for less than a multiple of 64
/// pays for the round-up and then throws the difference away.
const WORD: usize = 64;

/// **The precision this model needs, in bits.**
///
/// A realization's error radius is `C · 2⁻ᵖʳᵉᶜ`, where `C` depends on the *model* — its rotation
/// history and its coordinate magnitudes — and **not on the precision** (measured: identical `C`
/// at 256, 320, 384, 512 and 1024 bits). So one reading of `C` at any precision fixes the
/// precision needed to bring the radius under `limit`:
///
/// ```text
///     C · 2⁻ⁿᵉᵉᵈ ≤ limit    ⇒    need = log₂(C / limit)
/// ```
///
/// rounded up to a whole word. That round-up is free and its leftover is real confidence: asking
/// for 130 bits costs the same as 192, so take the 192.
///
/// `C` grows about **one bit per turn** in the rotation history, which is why a fixed precision
/// cannot work — it silently decides how long a model's history may be. At 256 bits a solid
/// turned 245 times stops building.
///
/// This is an *estimate*, and correctness does not rest on it: every judgement checks its own
/// interval, so an under-estimate costs a re-run and never an answer. It is deliberately
/// generous by one word to cover the determinant arithmetic stacked on top of the coordinates.
pub fn judge_precision<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>, limit: Mag) -> usize {
    let mut worst = Mag::ZERO;
    for p in pts {
        let b = trial_bound(p);
        if worst.lt(b) {
            worst = b;
        }
    }
    precision_for(worst, limit)
}

/// One point's share of [`judge_precision`] — its realization error at the trial precision.
///
/// Split out because a model has hundreds of these and they do not depend on each other, so a
/// caller can evaluate them however it likes. **Combining them is a maximum**, which is
/// associative and exact, so no order of combination — and no schedule — can change the
/// answer. That is what makes this safe to hand out, where exposing a partial *sum* would not
/// be: a reassociated floating-point sum is a different number.
pub fn trial_bound(p: &WitnessPoint) -> Mag {
    let mut worst = Mag::ZERO;
    // **Uncached on purpose.** `hp_coord` fills a point's realization cell with whatever
    // precision asks first, and this measurement runs before the real precision is known —
    // so going through it would fill every cell at `TRIAL_PREC` and make every later
    // judgement miss and re-realize. That is the exact cost the cell exists to remove.
    for c in p.compute_hp(TRIAL_PREC) {
        if worst.lt(c.error) {
            worst = c.error;
        }
    }
    worst
}

/// The working precision that brings a model whose worst trial-precision error is `worst`
/// within `limit` — the arithmetic half of [`judge_precision`], once the maximum is known.
pub fn precision_for(worst: Mag, limit: Mag) -> usize {
    // `worst = C · 2⁻ᵗʳⁱᵃˡ`, so `C = worst · 2ᵗʳⁱᵃˡ` and `need = log₂C − log₂limit`.
    let (Some(w), Some(l)) = (worst.exp2(), limit.exp2()) else {
        return TRIAL_PREC; // an exact model, or no limit to reach — nothing to size
    };
    let need = (w + TRIAL_PREC as i64 - l).max(0) as usize;
    let words = need.div_ceil(WORD) + 1; // + one word for the determinants above the coordinates
    (words * WORD).max(TRIAL_PREC)
}

/// **A determinant is not a length, and the coincidence limit is.**
///
/// `orient3d(a, b, c, d) = det[a−d, b−d, c−d]` is a signed volume: six times the tetrahedron's.
/// Divide it by the area term `|(b−d) × (c−d)|` and what is left is the **height of `a` above the
/// plane through `b, c, d`** — a distance, in the model's own units, which is the only thing a
/// limit like "closer than 1e-54" can be compared against. Comparing the raw determinant instead
/// would make the threshold scale with the triangle's size, which is how a tolerance becomes a
/// number nobody can reason about.
///
/// Returns an upper bound on that distance, given the determinant's own interval: the value is
/// somewhere inside `±error`, so the distance is at most `error / |cross|` — and `|cross|` is itself
/// uncertain, so its **lower** bound is what divides. When the area term cannot be bounded away
/// from zero the answer is the [`Gap`] saying which failure it is: a triangle that has collapsed
/// has no plane to be a distance from, while one that is merely unresolved is a matter of depth.
pub(super) fn distance_bound(det_rad: Mag, cross: &[HpBounded; 3], prec: usize) -> Gap {
    // `|cross|² = Σ cross[k]²`, and a lower bound on the norm needs a lower bound on the sum.
    let mut lo = HpBounded::exact(BigFloat::from_f64(0.0, prec));
    for c in cross {
        lo = lo.add(&c.mul(c, prec), prec);
    }
    // `|cross| ≥ √(value − error)`, and the square root only halves the exponent, so working in
    // exponents avoids needing a high-precision sqrt at all.
    let sq_lo = match denom_lo(&lo) {
        Ok(b) => b,
        // The shortfall is on the *square*, so half of it separates the norm — and the halving
        // rounds up, since a bit too many costs a word and a bit too few costs another round.
        Err(Gap::Unresolved { short }) => {
            return Gap::Unresolved {
                short: short.div_ceil(2),
            };
        }
        Err(g) => return g,
    };
    let Some(e) = sq_lo.exp2() else {
        return Gap::Vanished;
    };
    // `√(m · 2^e) ≥ 2^(⌊e/2⌋ − 1)` for `m ∈ [0.5, 1)`, which is the bound we need below.
    let norm_lo = Mag::pow2(e.div_euclid(2) - 1);
    match det_rad.over(norm_lo) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}
