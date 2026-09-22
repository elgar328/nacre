//! Exact rational scalars — the overhaul's exact rational value engine (`Rat`/`Angle`) and
//! the axis/rotation/isometry value types the kernel builds on. (The toleranced-sign frame
//! judgment this value engine enables now lives in `nacre-judge`.)
//!
//! The truth layer for user-input dimensions and angles: a value the user typed
//! is preserved *exactly*, so `1.1` stays `11/10` and `1.1 × 7` is exactly `7.7`
//! (the "thin film" problem's root fix). This is the exact-*value* counterpart to
//! `nacre-predicates`, which decides exact *signs* of geometric determinants —
//! complementary, not redundant.
//!
//! - [`Rat`] — a rational scalar (tol 0). Fixed-width `Ratio<i128>` with
//!   **checked** arithmetic: overflow is a *signal* (the downgrade trigger),
//!   not a panic or silent wrap. The caller downgrades that value's cache to
//!   f64 and records its tol; the definition is never lost, so a
//!   judgement can still realize it at whatever precision it needs.
//! - [`Angle`] — rational degrees, normalized mod-360, with exact accumulation so
//!   a full turn lands back on exactly `0` (no f64 drift — a sketch closes). The
//!   `cos`/`sin` realization crosses into f64 (the irrational boundary); the
//!   90°-family (`0/90/180/270°`) realizes to exact rationals `{0, ±1}` (Niven),
//!   so those rotations of a rational point stay tol 0; and `cos_sin_at`
//!   realize in arbitrary precision (astro-float) for the judgment path.
//!
//! Scope: the exact value engine (the toleranced-sign frame judgment lives in `nacre-judge`).
//! [`Rat::from_decimal`] carries this into *construction*: a prism's
//! placement and sweep are done in the rationals the dimensions were **written** as, so
//! `1.1` then `6.6` reaches the same plane as `7.7` (`nacre_ops::exact`). That holds where
//! the sketch frame is exactly orthonormal — a rotated frame's axes are irrational, and
//! there this crate has nothing to offer; `nacre-judge` is what keeps *judgments* sound
//! there. Still deferred: the unified `Scalar { value, tol }` wrapper.
//! An unprovable sign leaves as a proved coincidence carrying its evidence, or as a reject named
//! for its cause (`nacre-judge`) — not as a declared zero awaiting user confirmation, which
//! measurement refuted.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
pub mod bounded;
pub mod mag;
pub mod quad;
pub use bounded::{Bounded, HpBounded, rat_to_big};
pub use mag::Mag;
pub use quad::{
    QuadVal, biquad_sign, cylinder_radial_side, cylinders_clear, cylinders_nested,
    exceeds_root_sum, segment_meets_cylinder,
};

use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedDiv, CheckedMul, CheckedSub};

use astro_float::{BigFloat, Consts, RoundingMode};
use std::cell::RefCell;
use std::collections::HashMap;

mod angle;
mod bigkernel;
mod cylinder;
mod frame;
mod isometry;
mod meet;
mod orient;
mod plane_name;
mod rat;
mod realize;
mod sqrt;
mod winding;

pub(crate) use bigkernel::*;
pub use cylinder::*;
pub use frame::*;
pub use isometry::*;
pub use meet::*;
pub use orient::*;
pub use plane_name::*;
pub use rat::*;
pub use realize::*;
pub use sqrt::*;
pub use winding::*;

/// Rounding for the high-precision realization layer (astro-float).
///
/// **There is no default precision here, on purpose.** A fixed one would be a hand-picked depth
/// waiting to be wired into a judgement whose precision belongs to the *model*
/// (`nacre_judge::judge_precision`). Callers pass `prec`.
///
/// Public because the judge (`nacre-judge`) rounds with it too — one mode in two places is one
/// more thing that can drift.
pub const HP_RM: RoundingMode = RoundingMode::ToEven;

thread_local! {
    /// Transcendental-constant cache (π, …) for the high-precision layer.
    ///
    /// **Constants, not results.** `astro_float::Consts` holds `pi/e/ln2/ln10/tenpowers` and
    /// nothing else, and `cos`/`sin` take it only to reach π (degrees → radians, and argument
    /// reduction). So it saves computing π once per call — a real saving, and orthogonal to
    /// [`TRIG`] below, which is what saves the *evaluation*.
    static HP_CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float consts"));

    /// **Realized `cos`/`sin` per `(angle, precision)`** — see [`Angle::cos_sin_bounded`].
    ///
    /// A rotation's realization asks for the *same* angle once per point: a 60-fin fold evaluated
    /// **29,037 of them for 112 distinct `(angle, prec)` pairs**, and a solid turned 4,200 times
    /// re-evaluated its one angle 4,200 times per point. astro-float's `cos`/`sin` run Ziv's loop
    /// (a series, then a retry at more bits when the value sits too close to a rounding boundary),
    /// which is what that repetition was paying for — measured, **22% of a rotated boolean**.
    ///
    /// **The key is a value, not a handle**, which is what makes a process-wide memo sound here:
    /// `Angle` is an exact rational number of degrees, so two models asking for the same angle are
    /// asking the same question. (A cache keyed by a `Handle` could not be shared this way — a
    /// handle is an index into *one* model's store.) And the result is a pure function of the key,
    /// so the memo cannot move an answer; only how often Ziv's loop runs. `prec` is *in* the key,
    /// so a model whose judging precision grows simply lands on a different entry rather than
    /// reading one realized at the wrong depth.
    static TRIG: RefCell<HashMap<(Angle, usize), TrigAt>> = RefCell::new(HashMap::new());

    /// **How far this platform's `f64` cos/sin land from the truth** — see
    /// [`Angle::realization_error_of`].
    ///
    /// ★ **No precision in the key, and that is the point.** [`TRIG`] is keyed by
    /// `(angle, precision)` because it holds a value *realized at* a depth. This holds how wrong
    /// the **f64** realization is, and f64 has one precision. The `P` used to measure it says how
    /// finely the answer is read, not what the answer is.
    ///
    /// ★★★ **The realized pair is in the key, not just the angle.** The quantity is how far *this
    /// `(cos, sin)`* sits from the truth, so the pair that was measured is part of the question.
    /// Keying on the angle alone would hand a caller the error of a realization it did not use —
    /// which is not hypothetical: `(deg.to_f64() * PI / 180.0).sin()` is measured to differ by
    /// 1 ulp between a debug and a release build, and between two call sites within one release
    /// build (LLVM folds the literal-angle site at compile time and leaves the other to libm).
    ///
    /// ★★ **So a second entry under one angle is a signal, not waste**: it says two realizations of
    /// that angle are live in this process. The tests read the count for exactly that reason.
    static F64_ERR: RefCell<HashMap<RealizedAt, (f64, f64)>> = RefCell::new(HashMap::new());

    /// **The f64 realization of an angle, keyed by the angle alone** — see [`Angle::cos_sin_f64`].
    ///
    /// ★ Unlike [`F64_ERR`], the pair is *not* in the key, and it must not be: the value here **is**
    /// the pair, and it is the correctly rounded one, which is unique. This memo therefore cannot
    /// change an answer — it only stops the arbitrary-precision realization from running per vertex.
    ///
    /// **Measured, which is why it is here from the start**: one realization costs ~32µs against
    /// libm's 22ns, and `Isometry::apply_point` runs per vertex.
    static F64_TRIG: RefCell<HashMap<Angle, (f64, f64)>> = RefCell::new(HashMap::new());

    /// How often the f64 realization needed a second, deeper pass, and how often even that left
    /// the rounding undecided — see [`Angle::cos_sin_f64`]. **Read by tests**: a checker that
    /// never fires is indistinguishable from one that is not running.
    static ROUND_ESCALATED: RefCell<(usize, usize)> = const { RefCell::new((0, 0)) };

    /// **Realized `1/√v` per `(rational, precision)`** — see [`inv_sqrt_bounded`].
    ///
    /// The same argument that makes [`TRIG`] sound applies unchanged: the key is a *value*, not a
    /// handle, so two models asking for the same squared length are asking the same question; the
    /// result is a pure function of the key, so the memo cannot move an answer; and `prec` is in
    /// the key, so a model judged more deeply lands on a different entry rather than reading one
    /// realized too shallowly.
    static INV_SQRT: RefCell<HashMap<(Rat, usize), HpBounded>> =
        RefCell::new(HashMap::new());

    /// **The f64 realization of `1/√v`, keyed by the rational alone** — see [`inv_sqrt_f64`].
    ///
    /// As with [`F64_TRIG`], the value is the correctly rounded one, which is unique, so this memo
    /// is a cost question only.
    static F64_INV_SQRT: RefCell<HashMap<Rat, f64>> = RefCell::new(HashMap::new());

    /// How often `1/√v`'s f64 realization needed the deeper rung, and how often even that left the
    /// rounding undecided — the counterpart of [`ROUND_ESCALATED`], read by tests for the same
    /// reason.
    static INV_SQRT_ESCALATED: RefCell<(usize, usize)> = const { RefCell::new((0, 0)) };
}

/// One `(angle, precision)` realization: `(cos, sin, |Δcos|, |Δsin|)`.
type TrigAt = (HpBounded, HpBounded);

/// An angle **together with one f64 realization of it** — `(angle, cos.to_bits(), sin.to_bits())`.
/// The bits, not the floats, because the key has to be `Hash` and `Eq`.
type RealizedAt = (Angle, u64, u64);

/// A rational scalar (exact, tol 0). Arithmetic returns `None` on i128 overflow
/// so the caller sees the downgrade trigger explicitly; on overflow the kernel
/// switches that value's cache to f64 and tags its `Origin` with the resulting
/// tol, while the definition (the op-log of input rationals) is preserved. Overflow is far from normal use — adversarial coprime-denominator
/// accumulation reaches it near the i128 ceiling (~122 bits).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rat(Ratio<i128>);

/// **A rational of any width** — the truth's spelling for a quantity that is the *square* of a
/// stated number: a cylinder's `r²`. Every `Rat` has a square, and this is where it fits when
/// `i128` does not (a 16-digit decimal below `1e-4` has a denominator whose square leaves `i128`),
/// so a statement is never refused for its square being wide — the reason `MeetPoint::Wide`
/// stands beside `Narrow`. Lowest terms, positive denominator (`Ratio`'s invariant).
///
/// A *storage and door* type: the exact predicates lift to `BigInt` anyway and take it directly;
/// the few `Rat` arithmetic sites on a squared radius ask [`BigRat::narrow`] and decline exactly
/// where `r·r` would overflow.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BigRat(Ratio<num_bigint::BigInt>);

/// A *direction* angle in degrees, kept normalized to `[0, 360)` exactly
/// (rational). Accumulation is exact: turning by a rational angle repeatedly and
/// completing a full turn lands back on exactly `0` — no f64 drift. `cos`/`sin`
/// realization crosses into f64 (deg→rad via π): the irrational-realization
/// boundary, past which a *judgement* realizes at arbitrary precision instead
/// ([`Angle::cos_sin_at`]). The angle stays exact; only its realized coordinate
/// carries tol.
///
/// This is the direction type. A multi-turn *amount* (helix pitch × turns, revolve
/// sweep > 360°) must preserve the turn count, so it belongs to a separate
/// *unnormalized* `Sweep` type — flat 2D sketches never need it, so it is left as
/// a documented companion, not built here.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Angle(Rat); // invariant: 0 <= inner < 360

#[cfg(test)]
#[path = "tests/lib/mod.rs"]
mod tests;

#[cfg(test)]
mod symprobe {
    use crate::{Angle, Rat};
    #[test]
    #[ignore = "probe"]
    fn measure_mirror_pairs() {
        let a = |d: i128| Angle::from_deg(Rat::from_int(d)).unwrap();
        for (x, y) in [(72i128, 288i128), (9, 351), (117, 243), (45, 315)] {
            let (cx, sx) = a(x).cos_sin_f64();
            let (cy, sy) = a(y).cos_sin_f64();
            println!(
                "[sym] {x} vs {y}: cos equal {}  sin negated {}   ({cx:.20e} / {cy:.20e})",
                cx == cy,
                sx == -sy
            );
        }
    }
}
