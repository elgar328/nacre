//! The toleranced point + its `orient3d` judgment.
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (`new tol = |R|·old + mix`). [`orient3d_judge`]
//! consumes that tol: an f64 determinant filter with a sound error bound
//! ([`det3_bound`]) decides the easy cases, the ambiguous ones **escalate**
//! to astro-float from the point definitions, and one whose interval still straddles zero is
//! **normalized into the distance it stands for** and held against the operation's coincidence
//! limit ([`Standard`]): below it the coincidence is proved, above it the judgement climbs to the
//! precision the shortfall names, and past the cap it is reported ([`Decision`]) instead of
//! assumed. This is the judgment **layer only**
//! — wired into the boolean since the CIP stages that followed (`nacre-ops`'
//! `tolerant` module is the seam), and it is the *tol > 0*
//! path: a tol-0 (`Constructed`) config is faster/exact via `nacre-predicates`
//! (Shewchuk), routed by a higher layer, not here.
//!
//! Validated before the port by an isolated 3D experiment (verdict: GO): H-a (`det3_bound`
//! soundness), H-d/H-f (chain + arbitrary-pivot tol) — the bound never under-estimates the true
//! error (astro-float ground truth) over random heterogeneous-rotation configs.

use super::HP_RM;
use astro_float::BigFloat;
use nacre_exact::{Angle, Axis, Bounded, HpBounded, Mag, Orient, Rat, rat_to_big};
#[cfg(feature = "parallel")]
use std::sync::{Arc as HpRc, OnceLock as HpOnce};
#[cfg(not(feature = "parallel"))]
use std::{cell::OnceCell as HpOnce, rc::Rc as HpRc};

mod chain;
mod coord;
mod det3;
mod frames;
mod indirect;
mod orient3d;
mod policy;
mod witness;

pub use chain::*;
pub use coord::*;
pub use det3::*;
pub use frames::*;
pub use indirect::*;
pub use orient3d::*;
pub use policy::*;
pub use witness::*;

/// Shared, lazily-initialized cell for the memoized high-precision realization.
///
/// Under `parallel` it is `Arc<OnceLock>` — `Send + Sync`. **What needs that is the shared
/// borrow**: the boolean hands every worker the same `&[WorkingPlane]`, so `WitnessPoint` must be `Sync`
/// or the plane table cannot cross the closure at all. The cache being shared rather than
/// per-thread is the second benefit: a hot definition point is realized once for all workers.
/// Two workers racing to fill one cell compute the same value and one wins, so the answer
/// does not depend on who did.
///
/// Otherwise it is `Rc<OnceCell>` — single-threaded, no atomic overhead. `get_or_init` has
/// the identical signature on both, so the consumer ([`WitnessPoint::hp_coord`]) is unchanged.
type HpCell = HpRc<HpOnce<(usize, [HpBounded; 3])>>;

/// One motion in a point's definition. `Rotate` turns about `axis` (the line through the
/// rational pivot `point`) by the rational `angle` — `point = [0,0,0]` is the origin-pivot case.
/// Motions do not commute, so the chain's **order is the definition**.
/// ★ **`Frame` is the large variant and it is boxed nowhere.** It holds four rational vectors and
/// three rational lengths against `Rotate`'s one vector — but a `MoveNode` lives in a shared
/// `Rc<[MoveNode]>` that the judgment path clones by refcount, and a chain is a handful of nodes
/// per point at most. Boxing would trade a flat read for a pointer chase on the hot path to save
/// bytes nothing is short of.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoveNode {
    Rotate {
        axis: Axis,
        angle: Angle,
        pivot: [Rat; 3],
    },
    /// An exact rational translation. Realized by adding the offset — see [`WitnessPoint::compute_hp`].
    Translate { offset: [Rat; 3] },
    /// An exact reflection in `axis = offset` (`x ↦ 2·offset − x` on that axis).
    ///
    /// **Improper** (`det = −1`): unlike the other two it negates a determinant of its images.
    /// A judgement that cancels a shared motion out of a determinant must bring the pre-motion
    /// data into the same handedness first — see [`shared_base`].
    Mirror { axis: Axis, offset: Rat },
    /// A change of basis **into a plane's own frame**: the point's coordinates are read as
    /// `(u, v, w)` there and carried out into the frame the plane lives in.
    ///
    /// ```text
    /// ŵ = n / |n|        û = u_raw / |u_raw|        v̂ = ŵ × û
    /// out = origin + u·û + v·v̂ + w·ŵ
    /// ```
    ///
    /// ★★★ **Every input is exact and rational; the only irrational step is dividing by a
    /// length.** `u_raw` is `ẑ × n` (or `ŷ × n`), which is *already in the plane* — the cross
    /// product is perpendicular to both its arguments, so the projection a general `ref_dir`
    /// would need is not here, and `u_raw` need not be a unit vector for the frame to be exact.
    /// So the whole realization is **two `1/√(rational)` scalars**
    /// ([`nacre_exact::inv_sqrt_f64`]), which is why this can be judged at all.
    ///
    /// ★★ **`v̂` is exact too, and that matters.** Realizing it as `ŵ × û` in f64 costs two
    /// roundings that do not cancel — a wall whose `v` is exactly `ẑ` came out three ulps short,
    /// which is a frame that is not quite orthonormal. `n ⊥ u_raw` by construction, so
    /// `|v_raw|² = |n|²·|u_raw|²` exactly, and one inverse square root of a rational lands it on
    /// the nose. [`nacre_exact::plane_frame`] checks that product fits `i128` before handing the
    /// frame out, so the overflow is refused at the source rather than handled here.
    ///
    /// ★ **Proper** (`det = +1`), so it contributes nothing to a chain's mirror parity.
    Frame { frame: nacre_exact::PlaneFrame },
    /// [`MoveNode::Frame`] for a plane whose exact data does not fit `Rat` — a `Wide`
    /// name, or a narrow one whose squared lengths overflow `i128`. Same realization shape,
    /// arbitrary-precision integers instead: **nothing here can overflow**, so unlike
    /// `PlaneFrame` there is no partial (`v: None`) form.
    ///
    /// ★ This variant is why `MoveNode` is no longer `Copy` — `BigInt` owns heap. The chain is
    /// shared by `Rc`, so nothing hot copies nodes.
    ///
    /// ★ **Proper** (`det = +1`), like [`MoveNode::Frame`].
    FrameWide(WideFrame),
    /// [`MoveNode::Frame`] for a plane that has **no name at all** — a mixed-frame datum, whose
    /// exact world coefficients are irrational so neither exact vessel can hold them. The node
    /// carries the plane's three defining points instead (exact each in its own frame), and the
    /// basis is *judged*: derived as intervals at whatever precision the realization runs at.
    /// See [`FrameThrough`].
    ///
    /// ★ **Proper** (`det = +1`) — `v̂ = ŵ × û` by construction, like the other frame nodes, so
    /// [`chain_parity`]'s mirrors-only count stays correct.
    ///
    /// Boxed: the three points dwarf every other variant.
    FrameThrough(Box<FrameThrough>),
}

/// A high-precision interval narrowed to an `(f64 value, f64 error)` pair — how a wide frame's
/// realization reaches the f64 cache. Correctly rounded when the interval decides it (then the
/// error is the radius plus the value's own half-ulp); an undecided or out-of-normal-range value
/// flushes to zero with its whole magnitude charged — sound, and reachable only for axis
/// components below f64's normal floor.
fn narrow_hp(x: &HpBounded) -> (f64, f64) {
    match nacre_exact::round_to_f64(&x.value, x.error, 128) {
        Some(v) => (v, rad_f64(x.error) + v.abs() * f64::EPSILON),
        None => (0.0, rad_f64(x.error) + bf_mag(&x.value)),
    }
}

/// An upper `f64` for a [`Mag`] radius: `2^e` from the exponent. Below f64's floor it lands on
/// the smallest positive value rather than a zero that would claim exactness; above the range it
/// is honestly infinite (an infinite tol declines, never lies).
fn rad_f64(r: Mag) -> f64 {
    match r.exp2() {
        None => 0.0, // a genuinely zero radius
        Some(e) if e < -1074 => f64::MIN_POSITIVE,
        Some(e) if e > 1023 => f64::INFINITY,
        Some(e) => 2f64.powi(e as i32),
    }
}

/// Run one judgement at `j.prec` and, while it neither decides a sign nor proves a coincidence,
/// again with more bits — up to the cap.
///
/// `attempt` answers at a given precision with the sign, if it has one, and otherwise (`Err`) the
/// [`Gap`] its undecided determinant stands for, normalized into the unit `limit` is in. That separation
/// is the whole point: a determinant is not a length, and a threshold applied to one directly
/// would move with the size of the witness triangle.
///
/// **The next precision is computed, not doubled.** The gap is `C · 2⁻ᵖʳᵉᶜ` over a cofactor and
/// `C` does not depend on the precision (measured), so `log₂(gap / limit)` *is* the number of bits
/// missing — the same derivation [`judge_precision`] uses to size the model in the first place.
/// One jump lands there, rounded up to a whole word because astro-float allocates whole words
/// anyway. Doubling would either overshoot (paying for bits nobody asked for) or, on a model that
/// starts deep, undershoot and realize everything twice for nothing.
/// **How often a judgement had to climb, and how far.**
///
/// Every production escalation goes through [`escalate`], so one place here answers "what
/// fraction of judgements the f64 filter could not settle" and "what they cost" without
/// instrumentation scattered across the predicates.
///
/// ★ Unconditional, like `nacre-topo`'s `WIDE_PLANES` — and affordable for the same kind of
/// reason: this path has already decided to realize in arbitrary precision, so a relaxed atomic
/// add is noise beside the BigFloat work it is about to do. (`#[cfg(test)]` would not serve: the
/// measurements that read these live in another crate.)
pub mod climb_census {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

    /// Judgements that reached [`super::escalate`] at all.
    pub static CLIMBS: AtomicU64 = AtomicU64::new(0);
    /// Sum of the precisions those judgements finished at — divide by `CLIMBS` for the mean.
    pub static BITS: AtomicU64 = AtomicU64::new(0);
    /// Judgements the budget could not settle.
    pub static EXHAUSTED: AtomicU64 = AtomicU64::new(0);

    /// `(climbs, bits, exhausted)` — read, and reset so the next window is its own.
    pub fn take() -> (u64, u64, u64) {
        (
            CLIMBS.swap(0, Relaxed),
            BITS.swap(0, Relaxed),
            EXHAUSTED.swap(0, Relaxed),
        )
    }
}

#[cfg(test)]
#[path = "../tests/frame3/mod.rs"]
mod tests;
