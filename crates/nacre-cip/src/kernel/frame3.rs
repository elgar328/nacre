//! The toleranced point + its `orient3d` judgment (design.md §9 CIP, stage 2).
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** (§CIP ⑤) that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (§CIP ①: `new tol = |R|·old + mix`). [`orient3d_judge`]
//! consumes that tol: an f64 determinant filter with a sound error bound
//! ([`det3_bound`], §CIP ②) decides the easy cases, the ambiguous ones **escalate**
//! to astro-float from the point definitions, and one whose interval still straddles zero is
//! **normalized into the distance it stands for** and held against the operation's coincidence
//! limit ([`Standard`]): below it the coincidence is proved, above it the judgement climbs to the
//! precision the shortfall names, and past the cap it is reported ([`Decision`]) instead of
//! assumed. This is the judgment **layer only**
//! ("층만") — not yet wired into boolean (that is stage 3), and it is the *tol > 0*
//! path: a tol-0 (`Constructed`) config is faster/exact via `nacre-predicates`
//! (Shewchuk), routed by a higher layer, not here.
//!
//! Validated before the port by an isolated 3D experiment (verdict: GO): H-a (`det3_bound`
//! soundness), H-d/H-f (chain + arbitrary-pivot tol) — the bound never under-estimates the true
//! error (astro-float ground truth) over random heterogeneous-rotation configs.

use super::HP_RM;
use super::interval::{DA_F64, bf_mag, rat_to_big, rat_to_hp};
use super::interval::{HpIv, Iv};
use astro_float::BigFloat;
use nacre_scalar::{Angle, Axis, Bound, Orient, Rat};
#[cfg(feature = "parallel")]
use std::sync::{Arc as HpRc, OnceLock as HpOnce};
#[cfg(not(feature = "parallel"))]
use std::{cell::OnceCell as HpOnce, rc::Rc as HpRc};

/// Shared, lazily-initialized cell for the memoized high-precision realization.
///
/// Under `parallel` it is `Arc<OnceLock>` — `Send + Sync`. **What needs that is the shared
/// borrow**: the boolean hands every worker the same `&[PlaneGeom]`, so `Pt3` must be `Sync`
/// or the plane table cannot cross the closure at all. The cache being shared rather than
/// per-thread is the second benefit: a hot definition point is realized once for all workers.
/// Two workers racing to fill one cell compute the same value and one wins, so the answer
/// does not depend on who did.
///
/// Otherwise it is `Rc<OnceCell>` — single-threaded, no atomic overhead. `get_or_init` has
/// the identical signature on both, so the consumer ([`Pt3::hp_coord`]) is unchanged.
type HpCell = HpRc<HpOnce<(usize, [HpIv; 3])>>;

/// One motion in a point's definition. `Rotate` turns about `axis` (the line through the
/// rational pivot `point`) by the rational `angle` — `point = [0,0,0]` is the origin-pivot case.
/// Motions do not commute, so the chain's **order is the definition**.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoveNode {
    Rotate {
        axis: Axis,
        angle: Angle,
        point: [Rat; 3],
    },
    /// An exact rational translation. Realized by adding the offset — see [`Pt3::compute_hp`].
    Translate { offset: [Rat; 3] },
    /// An exact reflection in `axis = offset` (`x ↦ 2·offset − x` on that axis).
    ///
    /// **Improper** (`det = −1`): unlike the other two it negates a determinant of its images.
    /// A judgement that cancels a shared motion out of a determinant must bring the pre-motion
    /// data into the same handedness first — see [`shared_base`].
    Mirror { axis: Axis, offset: Rat },
}

/// A rational base point carried through a chain of axis rotations (§CIP ⑦ rotation
/// history). `base` + `chain` are the exact **definition** (never lost); `coord` is
/// the f64 realization (a cache), and `tol` bounds its error as a **direction-wise
/// xyz vector** (§CIP ⑤). [`hp_coord`](Self::hp_coord) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct Pt3 {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,
    pub coord: [f64; 3],
    pub tol: [f64; 3],
    /// Memoized `hp_coord` at the boolean's chosen precision — the astro-float realization
    /// once per definition and shared across clones (`Rc`). A judge escalates the *same*
    /// definition-point dozens of times per boolean (`plane_def` clones `tri_pt3` per call);
    /// without this each escalation replays the rotation's cos/sin at 200 bits, which dominates
    /// the rotated-boolean cost. `base`/`chain` never change after construction except through
    /// [`rotate_about`], which resets this cell, so the cached value always matches the
    /// definition (a pure, path-independent function). See [`HpCell`] for the
    /// `Arc<OnceLock>` (parallel) vs `Rc<OnceCell>` (serial) choice.
    hp: HpCell,
}

impl Pt3 {
    /// A point at `base`, tol seeded with the base→f64 rounding (a division; exactly
    /// 0 for an f64-representable base, positive otherwise). An axis a later chain
    /// never rotates keeps exactly this, which a per-axis tol check needs.
    pub fn at(base: [Rat; 3]) -> Self {
        let coord = [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()];
        // Actual rounding = |coord − base| at high precision; 2× the power-of-two
        // magnitude is a sound upper bound (0 when the base is exact).
        let round_tol = |b: Rat, c: f64| {
            let e = BigFloat::from_f64(c, 120).sub(&rat_to_big(b, 120), 120, HP_RM);
            2.0 * bf_mag(&e.abs())
        };
        Self::at_with_tol(
            base,
            [
                round_tol(base[0], coord[0]),
                round_tol(base[1], coord[1]),
                round_tol(base[2], coord[2]),
            ],
        )
    }

    /// A point whose coordinates are **exactly representable** as f64 — the axis-aligned case.
    ///
    /// `base` is that f64 (`try_from_f64` is exact: `mantissa · 2^exp`), so realizing it back
    /// yields the same f64 and the rounding tol is **`0` by construction**. `None` if a
    /// coordinate falls outside `Rat`'s exponent range.
    ///
    /// **Use this, not [`at`](Self::at), when the coordinates came from f64.** `at` *measures*
    /// the rounding at 120 bits; for this case that is nine BigFloat operations to compute a
    /// zero, and the boolean's hot path used to pay it hundreds of thousands of times.
    pub fn exact(coord: [f64; 3]) -> Option<Self> {
        let b = |x: f64| Rat::try_from_f64(x);
        Some(Self::at_with_tol(
            [b(coord[0])?, b(coord[1])?, b(coord[2])?],
            [0.0; 3],
        ))
    }

    /// A point at `base` with an explicit initial `tol` — for a root that already
    /// carries tol (a `Discovered` boolean seam), whose tol the chain then transports
    /// (`|R|·old`). [`at`](Self::at) is the Constructed case (initial tol = base
    /// rounding).
    pub fn at_with_tol(base: [Rat; 3], tol: [f64; 3]) -> Self {
        Pt3 {
            coord: [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
            base,
            chain: HpRc::from([] as [MoveNode; 0]),
            tol,
            hp: HpCell::default(),
        }
    }

    /// Extend the definition by a rotation about `axis` through the origin. See
    /// [`rotate_about`](Self::rotate_about).
    pub fn rotate(self, axis: Axis, angle: Angle) -> Self {
        self.rotate_about(axis, angle, [Rat::from_int(0); 3])
    }

    /// Extend the definition by a rotation about `axis` through the rational pivot
    /// `point`, updating the f64 cache and its direction-wise tol. The two in-plane
    /// coords are taken relative to the pivot, rotated, and shifted back. Same-axis
    /// 90°-family angles about the origin rotate exactly (tol 0, Niven). A non-origin
    /// pivot adds its own f64 rounding (`coord − p`, `p + …`, the pivot's Rat→f64) —
    /// present even for an exact angle — covered by the `piv` term (validated H-f).
    pub fn rotate_about(mut self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let (i, j) = axis.plane();
        let (px, py) = (point[i].to_f64(), point[j].to_f64());
        let (ci, cj) = (self.coord[i], self.coord[j]); // pre-rotation magnitudes for tol
        let (u, v) = (ci - px, cj - py);
        // cos/sin — exact (rational) for the 90°-family, else f64 (with realization tol).
        let (c, s, exact) = match angle.try_exact_cos_sin() {
            Some((cr, sr)) => (cr.to_f64(), sr.to_f64(), true),
            None => (angle.cos(), angle.sin(), false),
        };
        self.coord[i] = px + u * c - v * s;
        self.coord[j] = py + u * s + v * c;
        // Rotation-realization error (coordinate-mixing), 0 for an exact angle; plus
        // the pivot arithmetic (subtract p, re-add p, the pivot's own rounding),
        // exactly 0 for an origin pivot and present even for an exact angle otherwise.
        let rot = if exact {
            0.0
        } else {
            (u.abs() + v.abs()) * DA_F64
        };
        let piv = if px != 0.0 || py != 0.0 {
            (ci.abs() + cj.abs() + px.abs() + py.abs()) * DA_F64
        } else {
            0.0
        };
        let mix = rot + piv;
        let (ti, tj) = (self.tol[i], self.tol[j]);
        self.tol[i] = c.abs() * ti + s.abs() * tj + mix;
        self.tol[j] = s.abs() * ti + c.abs() * tj + mix;
        // Rebuild the shared slice with the new node appended. This runs when a solid is
        // *transformed*, never on the judgment path, so the copy is not hot — and in exchange
        // `clone` becomes a refcount bump instead of an allocation, which the judgment path
        // does hundreds of thousands of times.
        let mut nodes = self.chain.to_vec();
        nodes.push(MoveNode::Rotate { axis, angle, point });
        self.chain = HpRc::from(nodes);
        // The definition changed — invalidate the memoized hp of the old definition. A fresh
        // (unshared) cell, so clones made before this rotation keep their own cached value.
        self.hp = HpCell::default();
        self
    }

    /// The point reflected in `axis = offset` — one more link in the definition.
    ///
    /// **`coord` follows the producer's own route** (`AxisMirror::point`), so a replay of the
    /// definition reproduces the stored coordinate bit for bit. The exact offset lives in the
    /// chain, where [`compute_hp`](Self::compute_hp) realizes it.
    pub fn mirror(mut self, axis: Axis, offset: Rat) -> Self {
        let k = axis.index();
        let c = offset.to_f64();
        // `2·c` is exact (a power-of-two multiply), so the new error is the offset's own
        // realization plus the subtraction's half-ulp. Negation itself is exact, so the incoming
        // tol passes through unscaled.
        //
        // **The offset enters doubled** (`2·c`), so its realization error does too — the factor
        // here is `4 = 2 (the doubling) × 2 (the same safety margin every other term carries)`.
        // Copying `translate`'s `2.0` looks right and is not: there the offset enters once, so its
        // `2.0` *was* the margin. Measured — a chain of two reflections overran the bound by 4%.
        self.tol[k] += 4.0
            * bf_mag(&rat_to_big(offset, 120).sub(&BigFloat::from_f64(c, 120), 120, HP_RM)).abs()
            + f64::EPSILON * (2.0 * c.abs() + self.coord[k].abs());
        // The producer's own route (`AxisMirror::point`), operation for operation — a replay of
        // the definition has to reproduce the stored coordinate bit for bit.
        self.coord[k] = 2.0 * c - self.coord[k];
        let mut nodes = self.chain.to_vec();
        nodes.push(MoveNode::Mirror { axis, offset });
        self.chain = HpRc::from(nodes);
        self.hp = HpCell::default();
        self
    }

    /// The point translated by an exact rational `offset` — one more link in the definition.
    ///
    /// **`coord` is updated exactly as the producer does it** (`Isometry::apply_point`: add the
    /// offset's f64 image), so a replay of the definition reproduces the stored coordinate bit for
    /// bit. The exact offset lives in the chain, where [`compute_hp`](Self::compute_hp) realizes
    /// it; the two roundings this f64 step takes (the offset's own, and the add's) are what `tol`
    /// grows by.
    ///
    /// That split is the whole point: two placements that reach the same real wall by different
    /// routes keep f64 coordinates an ulp apart, but their *definitions* realize to the same
    /// value, and it is the definition the judgment reads.
    pub fn translate(mut self, offset: [Rat; 3]) -> Self {
        for (k, &off) in offset.iter().enumerate() {
            let t = off.to_f64();
            // The offset's realization error, plus the add's own half-ulp on the result.
            self.tol[k] += 2.0
                * bf_mag(&rat_to_big(off, 120).sub(&BigFloat::from_f64(t, 120), 120, HP_RM)).abs()
                + f64::EPSILON * (self.coord[k].abs() + t.abs());
            self.coord[k] += t;
        }
        let mut nodes = self.chain.to_vec();
        nodes.push(MoveNode::Translate { offset });
        self.chain = HpRc::from(nodes);
        self.hp = HpCell::default();
        self
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain, each node about its pivot) — path-independent ground truth /
    /// escalation realization. The result is memoized in
    /// [`Pt3::hp`] and shared across clones of the same definition, so a definition-point pays
    /// the astro-float cos/sin once per boolean rather than once per predicate.
    pub(crate) fn hp_coord(&self, prec: usize) -> [HpIv; 3] {
        // **Keyed by precision, and that key is load-bearing.** The precision is chosen per
        // boolean, so within one operation every call arrives with the same value and the cell is
        // filled once — which is the whole point, since re-realizing a definition per predicate
        // was 77% of a rotated boolean's runtime. A call at a different precision (the rare
        // per-judgement fallback) recomputes without disturbing the cached value, rather than
        // silently returning coordinates realized at the wrong precision.
        let cached = self.hp.get_or_init(|| (prec, self.compute_hp(prec)));
        if cached.0 == prec {
            return cached.1.clone();
        }
        self.compute_hp(prec)
    }

    /// The uncached realization (the body of [`hp_coord`]), **with the error it carries**.
    ///
    /// This is the same walk as [`rotate_about`](Self::rotate_about)'s tol propagation, one level
    /// up: the definition is exact, the realization is not, and the radius is what the realization
    /// cost. Nothing here is a chosen constant — the base contributes its own rounding (zero when
    /// the rational lands on a `prec`-bit dyadic), each `cos`/`sin` contributes the bound
    /// [`Angle::cos_sin_bounded`] derives, and every arithmetic operation adds its half-ulp.
    fn compute_hp(&self, prec: usize) -> [HpIv; 3] {
        let mut p = [
            rat_to_hp(self.base[0], prec),
            rat_to_hp(self.base[1], prec),
            rat_to_hp(self.base[2], prec),
        ];
        for node in self.chain.iter() {
            match node {
                MoveNode::Rotate { axis, angle, point } => {
                    let (i, j) = axis.plane();
                    let (c, s, bc, bs) = angle.cos_sin_bounded(prec);
                    let (c, s) = (HpIv::new(c, bc), HpIv::new(s, bs));
                    let (px, py) = (rat_to_hp(point[i], prec), rat_to_hp(point[j], prec));
                    // pivot-relative: u = p − pivot, rotate, shift back.
                    let u = p[i].sub(&px, prec);
                    let v = p[j].sub(&py, prec);
                    p[i] = px.add(&u.mul(&c, prec).sub(&v.mul(&s, prec), prec), prec);
                    p[j] = py.add(&u.mul(&s, prec).add(&v.mul(&c, prec), prec), prec);
                }
                // A translation is exact input: the only error is `rat_to_hp`'s own division
                // rounding and the add's half-ulp, both of which the interval carries.
                MoveNode::Translate { offset } => {
                    for k in 0..3 {
                        p[k] = p[k].add(&rat_to_hp(offset[k], prec), prec);
                    }
                }
                // `2·offset − x` on one coordinate: exact input, so the interval carries only its
                // own rounding. The other two coordinates are untouched.
                MoveNode::Mirror { axis, offset } => {
                    let k = axis.index();
                    let c = rat_to_hp(*offset, prec);
                    p[k] = c.add(&c, prec).sub(&p[k], prec);
                }
            }
        }
        p
    }
}

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
    /// across the operation so [`Pt3`]'s realization cache stays warm.
    pub prec: usize,
    /// Two things **proved** to lie within this distance of each other are one thing.
    pub coincidence: Bound,
    /// The model's size — what turns the length limit into an angle for the one judgement whose
    /// question is about directions ([`dir_sign_judge`]).
    pub scale: Bound,
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
    Coincident { within: Bound },
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
    Exhausted { at: usize, within: Option<Bound> },
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
enum Gap {
    /// An upper bound on the separation, in the unit the judgement's limit is in.
    Of(Bound),
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
fn denom_lo(v: &HpIv) -> Result<Bound, Gap> {
    let Some(lo) = super::interval::lb(&v.mid) else {
        return Err(Gap::Vanished); // the midpoint is exactly zero
    };
    match lo.minus(v.rad) {
        Some(b) => Ok(b),
        // Nonzero, but the radius swallows it: the shortfall is how far the radius has to fall.
        None => Err(Gap::Unresolved {
            short: match (v.rad.exp2(), lo.exp2()) {
                (Some(r), Some(m)) => (r - m + 1).max(1) as usize,
                _ => 1,
            },
        }),
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
fn escalate(
    j: Standard,
    limit: Bound,
    mut attempt: impl FnMut(usize) -> Result<Orient, Gap>,
) -> Decision {
    // A zero precision means the context never got stamped: astro-float would be asked for a
    // realization with no bits, and every judgement would come back exhausted. That is a wiring
    // mistake, not a geometry one, so it fails loudly here rather than quietly answering `Zero`.
    debug_assert!(
        j.prec > 0,
        "escalate at zero precision — unstamped Standard"
    );
    let mut prec = j.prec;
    loop {
        let gap = match attempt(prec) {
            Ok(o) => return Decision::Sign(o),
            Err(g) => g,
        };
        // How many bits short this judgement is, and the bound it managed — the second is what
        // gets reported if the budget runs out. From a gap the shortfall is `log₂(gap / limit)`;
        // from an unresolved cofactor the cofactor itself named it, and there is no bound to
        // quote. Either way the number is *computed*, never guessed at by doubling.
        let (short, within) = match gap {
            Gap::Vanished => return Decision::Degenerate,
            Gap::Unresolved { short } => (short, None),
            // **A zero gap is a proof, not a near miss.** The radius bounds how far the computed
            // value is from the true one, so a zero radius says the midpoint *is* the value — and
            // the sign came back undecided only because that midpoint is zero. Determinants reach
            // this honestly: a row that is exactly zero (a plane perpendicular to the rotation
            // axis keeps its coordinate exactly) multiplies every error term to nothing.
            Gap::Of(g) if g.is_zero() => return Decision::Sign(Orient::Zero),
            Gap::Of(g) if !limit.lt(g) => return Decision::Coincident { within: g },
            Gap::Of(g) => match (g.exp2(), limit.exp2()) {
                (Some(gx), Some(lx)) => ((gx - lx).max(1) as usize, Some(g)),
                // A limit with no exponent (a zero bound) is a target no depth reaches.
                _ => {
                    return Decision::Exhausted {
                        at: prec,
                        within: Some(g),
                    };
                }
            },
        };
        if prec >= j.cap {
            return Decision::Exhausted { at: prec, within };
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
pub fn judge_precision<'a>(pts: impl IntoIterator<Item = &'a Pt3>, limit: Bound) -> usize {
    let mut worst = Bound::ZERO;
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
pub fn trial_bound(p: &Pt3) -> Bound {
    let mut worst = Bound::ZERO;
    // **Uncached on purpose.** `hp_coord` fills a point's realization cell with whatever
    // precision asks first, and this measurement runs before the real precision is known —
    // so going through it would fill every cell at `TRIAL_PREC` and make every later
    // judgement miss and re-realize. That is the exact cost the cell exists to remove.
    for c in p.compute_hp(TRIAL_PREC) {
        if worst.lt(c.rad) {
            worst = c.rad;
        }
    }
    worst
}

/// The working precision that brings a model whose worst trial-precision error is `worst`
/// within `limit` — the arithmetic half of [`judge_precision`], once the maximum is known.
pub fn precision_for(worst: Bound, limit: Bound) -> usize {
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
/// somewhere inside `±rad`, so the distance is at most `rad / |cross|` — and `|cross|` is itself
/// uncertain, so its **lower** bound is what divides. When the area term cannot be bounded away
/// from zero the answer is the [`Gap`] saying which failure it is: a triangle that has collapsed
/// has no plane to be a distance from, while one that is merely unresolved is a matter of depth.
fn distance_bound(det_rad: Bound, cross: &[HpIv; 3], prec: usize) -> Gap {
    // `|cross|² = Σ cross[k]²`, and a lower bound on the norm needs a lower bound on the sum.
    let mut lo = HpIv::exact(BigFloat::from_f64(0.0, prec));
    for c in cross {
        lo = lo.add(&c.mul(c, prec), prec);
    }
    // `|cross| ≥ √(mid − rad)`, and the square root only halves the exponent, so working in
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
    let norm_lo = Bound::pow2(e.div_euclid(2) - 1);
    match det_rad.over(norm_lo) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

/// The three edge rows of `orient3d(a,b,c,d) = det[a−d, b−d, c−d]`.
fn rows(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] - d[0], a[1] - d[1], a[2] - d[2]],
        [b[0] - d[0], b[1] - d[1], b[2] - d[2]],
        [c[0] - d[0], c[1] - d[1], c[2] - d[2]],
    ]
}

/// f64 `orient3d` determinant `det[a−d, b−d, c−d]`.
fn det3_f64(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let r = rows(a, b, c, d);
    r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
}

/// Error radius of a product `a·b·c` given each factor's value and tol (≥0):
/// `Π(|v|+τ) − Π|v|` — the interval product radius (worst-case, sound).
fn prod_err(va: f64, ta: f64, vb: f64, tb: f64, vc: f64, tc: f64) -> f64 {
    (va.abs() + ta) * (vb.abs() + tb) * (vc.abs() + tc) - va.abs() * vb.abs() * vc.abs()
}

/// Sound error bound on the f64 `orient3d` determinant from each point's directional
/// tol (§CIP ②). The determinant is six signed triple-products of the edge entries;
/// each entry `(a−d)[k]` carries tol `tol_a[k] + tol_d[k]`. The bound sums the six
/// product radii (input-tol propagation, triangle-inequality worst case) plus a term
/// for the f64 rounding of the determinant's own arithmetic. The 3D analogue of
/// [`crate::frame2`]'s 2D `det_bound`. Validated in exact3d (H-a).
fn det3_bound(p: [[f64; 3]; 4], t: [[f64; 3]; 4]) -> f64 {
    let r = rows(p[0], p[1], p[2], p[3]);
    let td = t[3]; // apex tol adds to every edge on subtraction
    let t = [
        [t[0][0] + td[0], t[0][1] + td[1], t[0][2] + td[2]],
        [t[1][0] + td[0], t[1][1] + td[1], t[1][2] + td[2]],
        [t[2][0] + td[0], t[2][1] + td[1], t[2][2] + td[2]],
    ];
    // The six signed triple products of the cofactor expansion (indices into r/t).
    let terms = [
        [(0, 0), (1, 1), (2, 2)],
        [(0, 0), (1, 2), (2, 1)],
        [(0, 1), (1, 0), (2, 2)],
        [(0, 1), (1, 2), (2, 0)],
        [(0, 2), (1, 0), (2, 1)],
        [(0, 2), (1, 1), (2, 0)],
    ];
    let mut input_tol = 0.0;
    let mut mag = 0.0;
    for tri in terms {
        let [(r0, c0), (r1, c1), (r2, c2)] = tri;
        input_tol += prod_err(
            r[r0][c0], t[r0][c0], r[r1][c1], t[r1][c1], r[r2][c2], t[r2][c2],
        );
        mag += r[r0][c0].abs() * r[r1][c1].abs() * r[r2][c2].abs();
    }
    // f64 rounding of the determinant's own ~17 operations. Round-to-nearest costs `ε/2` each, so
    // the accumulation is about `8.5·ε·mag`; `16·ε` is that with a factor of two over it.
    input_tol + 16.0 * f64::EPSILON * mag
}

/// 3×3 determinant of high-precision interval rows at `prec` bits — the same expression as
/// [`det3_iv`], one precision up, so the filter and the escalation cannot drift apart.
///
/// The entries come in **by reference**, and that is not a micro-optimization: Cramer's four
/// matrices are the same nine or twelve values in different arrangements, so a signature that
/// owns its rows makes the caller copy each value up to three times over. An `HpIv` copy is a
/// mantissa heap allocation, and this determinant is the hot path — [`cramer_hp`] alone made
/// forty-five of them per call for values it only ever read.
fn det3_big(r: [[&HpIv; 3]; 3], prec: usize) -> HpIv {
    let mul = |x: &HpIv, y: &HpIv| x.mul(y, prec);
    let m0 = mul(r[1][1], r[2][2]).sub(&mul(r[1][2], r[2][1]), prec);
    let m1 = mul(r[1][0], r[2][2]).sub(&mul(r[1][2], r[2][0]), prec);
    let m2 = mul(r[1][0], r[2][1]).sub(&mul(r[1][1], r[2][0]), prec);
    mul(r[0][0], &m0)
        .sub(&mul(r[0][1], &m1), prec)
        .add(&mul(r[0][2], &m2), prec)
}

/// [`det3_big`] over rows the caller already owns — the borrow, spelled once.
fn det3_big_rows(r: &[[HpIv; 3]; 3], prec: usize) -> HpIv {
    det3_big(
        [
            [&r[0][0], &r[0][1], &r[0][2]],
            [&r[1][0], &r[1][1], &r[1][2]],
            [&r[2][0], &r[2][1], &r[2][2]],
        ],
        prec,
    )
}

/// `orient3d` determinant realized at `prec` bits (astro-float) from the point
/// definitions — path-independent ground truth / escalation realization.
fn det3_hp(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3, prec: usize) -> HpIv {
    let (a, b, c, d) = (
        pa.hp_coord(prec),
        pb.hp_coord(prec),
        pc.hp_coord(prec),
        pd.hp_coord(prec),
    );
    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
    let rows = [
        [sub(&a[0], &d[0]), sub(&a[1], &d[1]), sub(&a[2], &d[2])],
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    ];
    det3_big_rows(&rows, prec)
}

/// The points' pre-motion coordinates, **brought into the moved frame's handedness**, when they
/// all carry one and the same motion and every base is `f64`-representable.
///
/// A motion preserves the determinants these judges take *up to its own determinant*, so under
/// that condition the answer on the bases is the answer — exactly. Chains are compared
/// structurally, so two spellings of one motion read as different: a missed cancellation is
/// slower, never wrong. Pivots must match too, since a rotation about `c` is `Rx + (c − Rc)` and
/// two pivots leave two translations that do not cancel in a difference.
///
/// **A reflection is improper (`det = −1`), and the correction belongs here rather than at each
/// caller.** Every point carries the same chain, so an odd number of reflections flips the
/// determinant *uniformly* — the shortcut would return a confidently wrong sign, not a
/// conservative miss. Rather than hand each judge a parity to multiply in (four sites, each
/// taking a different quantity, one of them already carrying a second sign convention), the bases
/// are handed back **already reflected once** when the parity is odd: the base frame then has the
/// moved frame's handedness and every determinant question transfers unchanged.
///
/// The canonical reflection is a sign flip on x, which is exact for every finite `f64`.
fn shared_base<const N: usize>(pts: &[&Pt3; N]) -> Option<[[f64; 3]; N]> {
    let first = pts[0];
    if first.chain.is_empty() {
        return None; // nothing to cancel; the caller's own fast path already handled this
    }
    for p in pts {
        if p.chain.len() != first.chain.len() {
            return None;
        }
        // Structural equality over the **whole** node, variant included: two motions that compare
        // equal here are declared one motion and the exact predicate then answers in their shared
        // pre-motion frame. A comparison that ignored the variant would answer a different
        // question with full confidence — silently wrong, not slow.
        if p.chain.iter().zip(first.chain.iter()).any(|(a, b)| a != b) {
            return None;
        }
    }
    let exact = |r: Rat| Rat::try_from_f64(r.to_f64()) == Some(r);
    if !pts.iter().all(|p| p.base.iter().all(|&r| exact(r))) {
        return None; // the exact predicate takes f64; a base that does not round-trip cannot go
    }
    // Odd parity ⇒ hand back the mirror image, so the base frame is handed like the moved one.
    let flip = if chain_parity(&first.chain) < 0 {
        -1.0
    } else {
        1.0
    };
    Some(std::array::from_fn(|i| {
        [
            flip * pts[i].base[0].to_f64(),
            pts[i].base[1].to_f64(),
            pts[i].base[2].to_f64(),
        ]
    }))
}

/// `−1` when a chain contains an odd number of reflections, `+1` otherwise — the factor an
/// improper motion puts on any determinant of its images.
pub fn chain_parity(chain: &[MoveNode]) -> i8 {
    let mirrors = chain
        .iter()
        .filter(|n| matches!(n, MoveNode::Mirror { .. }))
        .count();
    if mirrors % 2 == 0 { 1 } else { -1 }
}

/// **Everything [`orient3d_judge`] can answer without spending unbounded precision** — the f64
/// filter, then the shared-motion exact path. `None` means "not settled here".
///
/// **It never escalates, and that is the point.** A caller that only wants to *skip work when it
/// can prove the work is pointless* must not pay an astro-float climb to find out; and because
/// `None` sends it back to doing the work, a missed proof costs time and nothing else. So this
/// can only ever lose an optimization — never change an answer.
///
/// Split out of [`orient3d_judge`] rather than written beside it: two implementations of one
/// filter would be two error budgets to keep in step, and this kernel has already been bitten
/// twice by a question with two answering sites. The judge calls this, so they cannot drift, and
/// soundness here is not a property to test but a consequence of being the same code.
pub fn orient3d_filter(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3) -> Option<Orient> {
    let (a, b, c, d) = (pa.coord, pb.coord, pc.coord, pd.coord);
    let det = det3_f64(a, b, c, d);
    let bound = det3_bound([a, b, c, d], [pa.tol, pb.tol, pc.tol, pd.tol]);
    if det > bound {
        return Some(Orient::Positive);
    }
    if det < -bound {
        return Some(Orient::Negative);
    }
    // One shared rigid motion cancels out of `det[a−d, b−d, c−d]` (`det(R) = 1`), so the same
    // question is answered exactly on the pre-rotation coordinates — no tolerance, no escalation.
    // Its zero is a **proved** zero, which is why it is a sign and not a coincidence.
    let base = shared_base(&[pa, pb, pc, pd])?;
    Some(
        match nacre_predicates::orient3d(base[0], base[1], base[2], base[3]) {
            x if x > 0.0 => Orient::Positive,
            x if x < 0.0 => Orient::Negative,
            _ => Orient::Zero,
        },
    )
}

/// CIP `orient3d`: f64 filter (`|det| > bound` → trust the sign), else escalate to
/// astro-float at the operation's precision; a determinant that still straddles zero there is
/// turned into a distance and answered by [`escalate`]. Path-independent (a
/// function of the four point definitions). This is the *tol > 0* path — a tol-0
/// config is exact/faster via `nacre-predicates` (Shewchuk), routed above this crate.
pub fn orient3d_judge(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3, j: Standard) -> Decision {
    if let Some(o) = orient3d_filter(pa, pb, pc, pd) {
        return Decision::Sign(o);
    }
    escalate(j, j.coincidence, |prec| {
        let det = det3_hp(pa, pb, pc, pd, prec);
        match det.sign() {
            Some(pos) => Ok(orient_of(pos)),
            // Only now is the area term worth forming: the sign decides on the first try in all
            // but the coincident cases, and the cross product is six high-precision products.
            None => Err(distance_bound(det.rad, &cross_of(pb, pc, pd, prec), prec)),
        }
    })
}

/// `(b−d) × (c−d)` at `prec` bits — the area term that turns [`orient3d_judge`]'s determinant
/// into a height above the plane through `b, c, d`.
fn cross_of(pb: &Pt3, pc: &Pt3, pd: &Pt3, prec: usize) -> [HpIv; 3] {
    let (b, c, d) = (pb.hp_coord(prec), pc.hp_coord(prec), pd.hp_coord(prec));
    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
    let (u, v) = (
        [sub(&b[0], &d[0]), sub(&b[1], &d[1]), sub(&b[2], &d[2])],
        [sub(&c[0], &d[0]), sub(&c[1], &d[1]), sub(&c[2], &d[2])],
    );
    [
        u[1].mul(&v[2], prec).sub(&u[2].mul(&v[1], prec), prec),
        u[2].mul(&v[0], prec).sub(&u[0].mul(&v[2], prec), prec),
        u[0].mul(&v[1], prec).sub(&u[1].mul(&v[0], prec), prec),
    ]
}

/// How far `pa` may be from the plane through `pb, pc, pd`, in the model's own length units.
///
/// This is [`orient3d_judge`]'s determinant turned into a distance by [`distance_bound`]. When
/// the judge cannot decide a sign, this is the honest statement of what it *did* establish: not
/// "these are the same", but "`pa` is within **this much** of that plane". `None` when the three
/// plane points are too near collinear for a distance to mean anything.
pub fn orient3d_distance(pa: &Pt3, pb: &Pt3, pc: &Pt3, pd: &Pt3, prec: usize) -> Option<Bound> {
    let det = det3_hp(pa, pb, pc, pd, prec);
    // The *value* the judge could not separate from zero is somewhere in `±rad`, so the distance
    // it bounds is `rad / |cross|`.
    match distance_bound(det.rad, &cross_of(pb, pc, pd, prec), prec) {
        Gap::Of(b) => Some(b),
        _ => None,
    }
}

/// CIP `dir_orient3d`: the sign of `det[d, x−base, y−base] = d·((x−base)×(y−base))` — the
/// orientation of the ray direction `d` against the edge fan `(base→x, base→y)`. The
/// **direction analogue** of [`orient3d_judge`]: `point_in_solid`'s ray-triangle test asks
/// `orient3d(p, p+d, ·, ·)`, but `p+d` (a rotated point plus a rational offset) has no exact
/// `base+chain` `Pt3` (`R⁻¹d` is irrational). Every such determinant reduces to this form,
/// where `d` enters as one **exact** (rad-0) column and only `x, y, base` carry rotation tol.
/// Interval filter → astro-float escalation, exactly like the indirect judges; a below-floor
/// determinant that stays undecided is [`Orient::Zero`], absorbed by the caller's ray retry.
///
/// **The one judge with no metric normalization**, and the reason is that its caller does not
/// need one: a zero here means the ray runs along the face's plane, and `boolean` answers that by
/// casting a different ray rather than by asking how close it came.
///
/// `d` is realized through an unrotated [`Pt3`] purely to reuse `pt_iv`/`hp_coord` for its
/// coord/tol/hp — the row is `d` itself, never `d − base`.
pub fn dir_orient3d_judge(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3, prec: usize) -> Orient {
    let dp = Pt3::at(d);
    let (bi, xi, yi, di) = (pt_iv(base), pt_iv(x), pt_iv(y), pt_iv(&dp));
    let sub_iv = |u: [Iv; 3], v: [Iv; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
    if let Some(pos) = det3_iv([di, sub_iv(xi, bi), sub_iv(yi, bi)]).sign() {
        return orient_of(pos);
    }
    // Escalate: the same determinant at prec from the exact definitions.
    let (bh, xh, yh, dh) = (
        base.hp_coord(prec),
        x.hp_coord(prec),
        y.hp_coord(prec),
        dp.hp_coord(prec),
    );
    let sub_hp = |u: &HpIv, v: &HpIv| u.sub(v, prec);
    let rows = [
        dh,
        [
            sub_hp(&xh[0], &bh[0]),
            sub_hp(&xh[1], &bh[1]),
            sub_hp(&xh[2], &bh[2]),
        ],
        [
            sub_hp(&yh[0], &bh[0]),
            sub_hp(&yh[1], &bh[1]),
            sub_hp(&yh[2], &bh[2]),
        ],
    ];
    match det3_big_rows(&rows, prec).sign() {
        Some(pos) => orient_of(pos),
        None => Orient::Zero,
    }
}

/// CIP `orient3d(base, base+dir, x, y)` — an `orient3d` whose 2nd point is the ideal point
/// in direction `dir` (a ray from `base`). This is exactly the shape `ray_triangle_cross`'s
/// edge tests take (`orient3d(p, p+d, ·, ·)`), so a caller mirrors the f64 predicate
/// argument-for-argument instead of hand-reducing the direction column.
///
/// The determinant `det[base−y, (base+dir)−y, x−y]` column-reduces (`R1−R0 = dir`) and, after
/// the two swaps that move `dir` to the front, is `dir·((x−y)×(base−y))` — i.e. a plain
/// [`dir_orient3d_judge`] with the points permuted, no sign fix needed.
pub fn orient3d_ray(base: &Pt3, dir: [Rat; 3], x: &Pt3, y: &Pt3, prec: usize) -> Orient {
    dir_orient3d_judge(dir, y, x, base, prec)
}

// ---- indirect orient3d: three rotated planes meet at an implicit point ----
//
// A `Discovered` seam vertex is `∩` of three planes, each a plane through three
// rotated points, so its coordinates are never exact — an **implicit point** `V`.
// The kernel decides `orient3d(V, q, r, s)` without materializing `V`:
// `sign = sign(D)·sign(M)` where `D = det(normals)` and
// `M = (Dvec − D·s)·((q−s)×(r−s))` (Cramer, no division). The f64 filter is over
// **intervals** (value ± tol) — the design's "dynamic filter" (§CIP ⑨): interval
// arithmetic is a sound worst-case bound by construction, so no per-predicate bound
// formula is hand-derived. An interval straddling 0 escalates to astro-float from the
// point definitions, exactly as [`orient3d_judge`]. This is the *tol > 0* path, the
// judgment **layer only** ("층만") — the boolean wiring (seam → plane `Pt3`s) is
// stage 3. Validated before the port by the isolated 3D experiment (H-b coefficient
// tol, H-c indirect soundness); the constant `mag`-floor policy is indirect-only (distinct from the
// explicit `16·scale³` floor of [`orient3d_judge`]).

/// 3×3 determinant of interval rows.
fn det3_iv(r: [[Iv; 3]; 3]) -> Iv {
    let m0 = r[1][1].mul(r[2][2]).sub(r[1][2].mul(r[2][1]));
    let m1 = r[1][0].mul(r[2][2]).sub(r[1][2].mul(r[2][0]));
    let m2 = r[1][0].mul(r[2][1]).sub(r[1][1].mul(r[2][0]));
    r[0][0].mul(m0).sub(r[0][1].mul(m1)).add(r[0][2].mul(m2))
}

/// A point's coord+tol as an interval per component (the pivot is already folded into
/// `coord`/`tol` by [`Pt3::rotate_about`]).
fn pt_iv(p: &Pt3) -> [Iv; 3] {
    [
        Iv::new(p.coord[0], p.tol[0]),
        Iv::new(p.coord[1], p.tol[1]),
        Iv::new(p.coord[2], p.tol[2]),
    ]
}

/// Plane `[a,b,c,d]` (`n·X + d = 0`) through three points, as intervals: `n =
/// (p1−p0)×(p2−p0)`, `d = −n·p0`. Coefficient tol propagates from the point tols
/// through the subtraction/cross/dot — "coefficient tol is a corollary of point tol"
/// (§CIP ②, design.md §539). Validated H-b.
pub(crate) fn plane_iv(p0: &Pt3, p1: &Pt3, p2: &Pt3) -> [Iv; 4] {
    let (a, b, c) = (pt_iv(p0), pt_iv(p1), pt_iv(p2));
    let e1 = [b[0].sub(a[0]), b[1].sub(a[1]), b[2].sub(a[2])];
    let e2 = [c[0].sub(a[0]), c[1].sub(a[1]), c[2].sub(a[2])];
    let n = [
        e1[1].mul(e2[2]).sub(e1[2].mul(e2[1])),
        e1[2].mul(e2[0]).sub(e1[0].mul(e2[2])),
        e1[0].mul(e2[1]).sub(e1[1].mul(e2[0])),
    ];
    let d = Iv::new(0.0, 0.0)
        .sub(n[0].mul(a[0]))
        .sub(n[1].mul(a[1]))
        .sub(n[2].mul(a[2]));
    [n[0], n[1], n[2], d]
}

/// The plane's four coefficients realized at `prec` bits from the definitions
/// (ground truth for the coefficient tol; no trig of its own — consumes point
/// `hp_coord`).
fn plane_hp(p0: &Pt3, p1: &Pt3, p2: &Pt3, prec: usize) -> [HpIv; 4] {
    let (a, b, c) = (p0.hp_coord(prec), p1.hp_coord(prec), p2.hp_coord(prec));
    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
    let mul = |x: &HpIv, y: &HpIv| x.mul(y, prec);
    let e1 = [sub(&b[0], &a[0]), sub(&b[1], &a[1]), sub(&b[2], &a[2])];
    let e2 = [sub(&c[0], &a[0]), sub(&c[1], &a[1]), sub(&c[2], &a[2])];
    let n = [
        mul(&e1[1], &e2[2]).sub(&mul(&e1[2], &e2[1]), prec),
        mul(&e1[2], &e2[0]).sub(&mul(&e1[0], &e2[2]), prec),
        mul(&e1[0], &e2[1]).sub(&mul(&e1[1], &e2[0]), prec),
    ];
    let d = HpIv::exact(BigFloat::from_f64(0.0, prec))
        .sub(&mul(&n[0], &a[0]), prec)
        .sub(&mul(&n[1], &a[1]), prec)
        .sub(&mul(&n[2], &a[2]), prec);
    let [n0, n1, n2] = n; // `d` is already built from them, so the normal moves out rather than copying
    [n0, n1, n2, d]
}

/// `sign(D)·sign(M)` combined into an orientation (both must be definite).
fn combine(dsign: Option<bool>, msign: Option<bool>) -> Option<Orient> {
    match (dsign, msign) {
        (Some(dp), Some(mp)) => Some(if dp == mp {
            Orient::Positive
        } else {
            Orient::Negative
        }),
        _ => None,
    }
}

/// The interval Cramer parts of an implicit point `V = ∩(planes)`: `D = det(normals)`
/// and the numerator vector `Dvec` (column `j` replaced by `h = −d`), so `V[j] =
/// Dvec[j]/D` (no division taken here). Shared by [`indirect_filter`] (orient3d) and
/// [`cmp_filter`] (cmp_coord).
fn cramer_iv(planes: [[Iv; 4]; 3]) -> (Iv, [Iv; 3]) {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    let h = |k: usize| Iv::new(0.0, 0.0).sub(planes[k][3]); // n·X = h, h = −d
    let (n0, n1, n2) = (n(0), n(1), n(2));
    let d = det3_iv([n0, n1, n2]);
    // Cramer numerators Dvec: replace column j with h.
    let col_h = [h(0), h(1), h(2)];
    let dvec = [
        det3_iv([
            [col_h[0], n0[1], n0[2]],
            [col_h[1], n1[1], n1[2]],
            [col_h[2], n2[1], n2[2]],
        ]),
        det3_iv([
            [n0[0], col_h[0], n0[2]],
            [n1[0], col_h[1], n1[2]],
            [n2[0], col_h[2], n2[2]],
        ]),
        det3_iv([
            [n0[0], n0[1], col_h[0]],
            [n1[0], n1[1], col_h[1]],
            [n2[0], n2[1], col_h[2]],
        ]),
    ];
    (d, dvec)
}

/// The Cramer `D` **alone** — the determinant of the three planes' normals.
///
/// ★ **A caller that wants only `D` must not go through [`cramer_iv`]**, which builds four
/// determinants and the `h` column for the three it does not need. [`dir_sign_judge`] asks exactly
/// this question ("how does the line `a ∩ b` run relative to `c`"), and it used to throw away
/// three-quarters of the work on every call.
fn normals_det_iv(planes: [[Iv; 4]; 3]) -> Iv {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    det3_iv([n(0), n(1), n(2)])
}

/// [`normals_det_iv`] at `prec` bits — the `d` half of [`cramer_hp`] without its `h` column or its
/// three `Dvec` determinants. The escalation is where a wasted determinant costs the most.
fn normals_det_hp(planes: &[[HpIv; 4]; 3], prec: usize) -> HpIv {
    let n = |k: usize, j: usize| &planes[k][j];
    det3_big(
        [
            [n(0, 0), n(0, 1), n(0, 2)],
            [n(1, 0), n(1, 1), n(1, 2)],
            [n(2, 0), n(2, 1), n(2, 2)],
        ],
        prec,
    )
}

/// The interval f64 filter for `orient3d(V, q, r, s)`, `V = ∩(planes)`. `None` if
/// either `D` or `M` straddles 0 (escalate). Coefficient-direct (no division).
#[cfg(test)]
fn indirect_filter(planes: [[Iv; 4]; 3], q: [Iv; 3], r: [Iv; 3], s: [Iv; 3]) -> Option<Orient> {
    filter_rest(cramer_iv(planes), q, r, s)
}

/// [`indirect_filter`] with the implicit point's Cramer parts **already in hand**.
///
/// ★ **The split exists because `(D, Dvec)` *is* the point** — it does not mention `q, r, s`, so a
/// caller asking about one point against several query triangles pays for it once. The arrangement's
/// crossing collector asks twice in a row (a segment's two endpoints), which is what
/// [`indirect_orient3d2_judge_pre`] exploits.
fn filter_rest((d, dvec): (Iv, [Iv; 3]), q: [Iv; 3], r: [Iv; 3], s: [Iv; 3]) -> Option<Orient> {
    // row1 = Dvec − D·s ; cross = (q−s)×(r−s) ; M = row1·cross.
    let row1 = [
        dvec[0].sub(d.mul(s[0])),
        dvec[1].sub(d.mul(s[1])),
        dvec[2].sub(d.mul(s[2])),
    ];
    let dq = [q[0].sub(s[0]), q[1].sub(s[1]), q[2].sub(s[2])];
    let dr = [r[0].sub(s[0]), r[1].sub(s[1]), r[2].sub(s[2])];
    let cross = [
        dq[1].mul(dr[2]).sub(dq[2].mul(dr[1])),
        dq[2].mul(dr[0]).sub(dq[0].mul(dr[2])),
        dq[0].mul(dr[1]).sub(dq[1].mul(dr[0])),
    ];
    let m = row1[0]
        .mul(cross[0])
        .add(row1[1].mul(cross[1]))
        .add(row1[2].mul(cross[2]));
    combine(d.sign(), m.sign())
}

/// The Cramer parts of `V = ∩(planes)` at `prec` bits: `(D, Dvec)` — the determinant and its
/// numerator vector, each carrying the error radius accumulated along the way. Shared by
/// [`indirect_hp`] (orient3d) and [`cmp_hp`] (cmp_coord).
///
/// The magnitude bounds these used to return alongside are gone: the radius rides *with* the
/// value now, so there is nothing left for a caller to forget to use — which is exactly how the
/// indirect judge came to bound `M` by a scale that had already cancelled.
fn cramer_hp(planes: &[[HpIv; 4]; 3], prec: usize) -> (HpIv, [HpIv; 3]) {
    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
    let zero = HpIv::exact(BigFloat::from_f64(0.0, prec));
    // The four matrices are the same twelve values rearranged, so they are built as *views* of
    // the coefficients rather than copies of them — see [`det3_big`].
    let n = |k: usize, j: usize| &planes[k][j];
    let hc = [
        sub(&zero, &planes[0][3]),
        sub(&zero, &planes[1][3]),
        sub(&zero, &planes[2][3]),
    ];
    let d = det3_big(
        [
            [n(0, 0), n(0, 1), n(0, 2)],
            [n(1, 0), n(1, 1), n(1, 2)],
            [n(2, 0), n(2, 1), n(2, 2)],
        ],
        prec,
    );
    let dvec = [
        det3_big(
            [
                [&hc[0], n(0, 1), n(0, 2)],
                [&hc[1], n(1, 1), n(1, 2)],
                [&hc[2], n(2, 1), n(2, 2)],
            ],
            prec,
        ),
        det3_big(
            [
                [n(0, 0), &hc[0], n(0, 2)],
                [n(1, 0), &hc[1], n(1, 2)],
                [n(2, 0), &hc[2], n(2, 2)],
            ],
            prec,
        ),
        det3_big(
            [
                [n(0, 0), n(0, 1), &hc[0]],
                [n(1, 0), n(1, 1), &hc[1]],
                [n(2, 0), n(2, 1), &hc[2]],
            ],
            prec,
        ),
    ];
    (d, dvec)
}

/// The same `sign(D)·sign(M)` realized at `prec` bits (astro-float) — ground truth /
/// escalation. Returns the two determinants as intervals; their radii are what decides whether
/// either sign may be used.
///
/// **This is where the cancellation bug lived.** `row1 = Dvec − D·s` collapses to nothing when
/// the query point coincides with the implicit point, and the old code sized the declare-0 floor
/// off that collapsed value, so a floor of `1e-131` let a `8.6e-78` rounding residue through as a
/// confident sign while two other ways of asking the same question answered zero. An interval
/// cannot make that mistake: the radius of `row1` is the sum of what went into it, and a
/// subtraction that cancels leaves the radius behind.
fn indirect_hp(
    planes: [[HpIv; 4]; 3],
    q: [HpIv; 3],
    r: [HpIv; 3],
    s: [HpIv; 3],
    prec: usize,
) -> (HpIv, HpIv, [HpIv; 3]) {
    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
    let mul = |x: &HpIv, y: &HpIv| x.mul(y, prec);
    let add = |x: &HpIv, y: &HpIv| x.add(y, prec);
    let (d, dvec) = cramer_hp(&planes, prec);
    let row1 = [
        sub(&dvec[0], &mul(&d, &s[0])),
        sub(&dvec[1], &mul(&d, &s[1])),
        sub(&dvec[2], &mul(&d, &s[2])),
    ];
    let dq = [sub(&q[0], &s[0]), sub(&q[1], &s[1]), sub(&q[2], &s[2])];
    let dr = [sub(&r[0], &s[0]), sub(&r[1], &s[1]), sub(&r[2], &s[2])];
    let cross = [
        sub(&mul(&dq[1], &dr[2]), &mul(&dq[2], &dr[1])),
        sub(&mul(&dq[2], &dr[0]), &mul(&dq[0], &dr[2])),
        sub(&mul(&dq[0], &dr[1]), &mul(&dq[1], &dr[0])),
    ];
    let m = add(
        &add(&mul(&row1[0], &cross[0]), &mul(&row1[1], &cross[1])),
        &mul(&row1[2], &cross[2]),
    );
    // `cross` rides out with the two determinants rather than the finished gap: forming the gap
    // costs a division and an exponent walk, and the path that decides a sign never needs it.
    (d, m, cross)
}

/// **How far the implicit point may be from the triangle's plane**, in the model's own units.
///
/// `row1 = Dvec − D·s = D·(V − s)` and `M = row1 · cross`, so `M = D · (V−s)·cross`. The signed
/// distance from `V` to the plane through `q, r, s` is `(V−s)·n̂ = (V−s)·cross / |cross|`, hence
///
/// ```text
///     distance = M / (|D| · |cross|)
/// ```
///
/// Both denominators are needed: `|D|` because `V` is `Dvec/D` and the point's coordinates carry
/// that division, `|cross|` because the dot product carries the triangle's area. Bounded above
/// from `M`'s radius and below by the denominators' own lower bounds. `None` when either cannot
/// be kept away from zero — three near-parallel planes have no well-defined meeting point, and a
/// degenerate triangle no plane.
fn plane_gap(m: &HpIv, d: &HpIv, cross: &[HpIv; 3], prec: usize) -> Gap {
    let d_lo = match denom_lo(d) {
        Ok(b) => b,
        Err(g) => return g,
    };
    match distance_bound(m.rad, cross, prec) {
        Gap::Of(b) => match b.over(d_lo) {
            Some(b) => Gap::Of(b),
            None => Gap::Vanished,
        },
        g => g,
    }
}

/// CIP indirect `orient3d(V, q, r, s)`, `V = ∩(3 planes)` — each plane through three
/// rotated points, the triangle three rotated points. Interval filter → astro-float
/// escalation; an undecided `D` or `M` becomes the distance from the implicit point to the
/// triangle's plane and is answered by [`escalate`]. The `Discovered`-seam analogue of
/// [`orient3d_judge`]; boolean wiring is stage 3.
#[allow(clippy::too_many_arguments)]
pub fn indirect_orient3d_judge(
    plane_a: (&Pt3, &Pt3, &Pt3),
    plane_b: (&Pt3, &Pt3, &Pt3),
    plane_c: (&Pt3, &Pt3, &Pt3),
    q: &Pt3,
    r: &Pt3,
    s: &Pt3,
    j: Standard,
) -> Decision {
    let iv = [
        plane_iv(plane_a.0, plane_a.1, plane_a.2),
        plane_iv(plane_b.0, plane_b.1, plane_b.2),
        plane_iv(plane_c.0, plane_c.1, plane_c.2),
    ];
    indirect_orient3d_judge_pre(iv, plane_a, plane_b, plane_c, q, r, s, j)
}

/// [`indirect_orient3d_judge`] with the filter's three interval planes **already built**.
///
/// ★ The filter's inputs are a function of the three plane *definitions* alone, so a caller that
/// asks about one implicit point against many query triangles rebuilds them every time. Measured on
/// one rotated fixture, that rebuild is ~20-24% of a certified judgement, and in the arrangement's
/// crossing collector the same three definitions are handed in hundreds of times over. The `Pt3`
/// definitions stay in the signature because the escalation still needs them: `plane_hp` realizes
/// the coefficients afresh at each precision, and an interval cannot be sharpened after the fact.
#[allow(clippy::too_many_arguments)]
pub(crate) fn indirect_orient3d_judge_pre(
    planes: [[Iv; 4]; 3],
    plane_a: (&Pt3, &Pt3, &Pt3),
    plane_b: (&Pt3, &Pt3, &Pt3),
    plane_c: (&Pt3, &Pt3, &Pt3),
    q: &Pt3,
    r: &Pt3,
    s: &Pt3,
    j: Standard,
) -> Decision {
    orient3d_from_cramer(cramer_iv(planes), plane_a, plane_b, plane_c, q, r, s, j)
}

/// **Two query triangles against one implicit point, sharing its Cramer parts.**
///
/// `cramer_iv` is four `det3`s and does not mention the query, so asking about the same point twice
/// paid for it twice. The collector's containment test is exactly that shape: a segment's two
/// endpoints, both against the point where the query wall crosses the segment's line.
///
/// The caller's claim is only that the *three planes* are the same; that they name one point follows.
#[allow(clippy::too_many_arguments)]
pub(crate) fn indirect_orient3d2_judge_pre(
    planes: [[Iv; 4]; 3],
    plane_a: (&Pt3, &Pt3, &Pt3),
    plane_b: (&Pt3, &Pt3, &Pt3),
    plane_c: (&Pt3, &Pt3, &Pt3),
    t0: (&Pt3, &Pt3, &Pt3),
    t1: (&Pt3, &Pt3, &Pt3),
    j: Standard,
) -> (Decision, Decision) {
    let cr = cramer_iv(planes);
    let one = |t: (&Pt3, &Pt3, &Pt3)| {
        orient3d_from_cramer(cr, plane_a, plane_b, plane_c, t.0, t.1, t.2, j)
    };
    (one(t0), one(t1))
}

#[allow(clippy::too_many_arguments)]
fn orient3d_from_cramer(
    cr: (Iv, [Iv; 3]),
    plane_a: (&Pt3, &Pt3, &Pt3),
    plane_b: (&Pt3, &Pt3, &Pt3),
    plane_c: (&Pt3, &Pt3, &Pt3),
    q: &Pt3,
    r: &Pt3,
    s: &Pt3,
    j: Standard,
) -> Decision {
    if let Some(o) = filter_rest(cr, pt_iv(q), pt_iv(r), pt_iv(s)) {
        return Decision::Sign(o);
    }
    // Escalate: the same two determinants at prec, each carrying the radius accumulated
    // along its own computation. A determinant whose interval straddles zero is undecided,
    // and what it *did* establish is how far the implicit point may be from the triangle's plane.
    escalate(j, j.coincidence, |prec| {
        let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
        let (d, m, cross) = indirect_hp(
            [ph(plane_a), ph(plane_b), ph(plane_c)],
            q.hp_coord(prec),
            r.hp_coord(prec),
            s.hp_coord(prec),
            prec,
        );
        match combine(d.sign(), m.sign()) {
            Some(o) => Ok(o),
            None => Err(plane_gap(&m, &d, &cross, prec)),
        }
    })
}

// ---- indirect cmp_coord: order two implicit points along one axis ----
//
// Each implicit point's axis coordinate is `Dvec[axis]/D` (Cramer), so
// `sign(a[axis] − b[axis]) = sign(Dvec_a[axis]·D_b − Dvec_b[axis]·D_a)·sign(D_a)·
// sign(D_b)` — three exact signs, no division, matching `nacre_predicates::
// indirect_cmp_coord`. Interval filter → astro-float escalation, like orient3d. The
// result is invariant to each triple's plane-normal orientation: flipping a triple's
// normals negates both its `D` and `Dvec[axis]` (so the point is unchanged) and negates
// `M`, leaving `sign(M)·sign(D)` fixed — so a caller may define each plane by any three
// non-collinear points on it (the bridge relies on this).

/// Combine the three definite signs of `M·D_a·D_b` (`true` = positive) into an
/// ordering: `Positive` (`a[axis] > b[axis]`) iff an even number are negative; `None`
/// if any sign is indefinite.
fn cmp_combine(sm: Option<bool>, sda: Option<bool>, sdb: Option<bool>) -> Option<Orient> {
    match (sm, sda, sdb) {
        (Some(m), Some(a), Some(b)) => {
            let negatives = [m, a, b].iter().filter(|&&s| !s).count();
            Some(if negatives % 2 == 0 {
                Orient::Positive
            } else {
                Orient::Negative
            })
        }
        _ => None,
    }
}

/// The interval f64 filter for `cmp_coord(a, b, axis)` — the sign of `a[axis] −
/// b[axis]` between two implicit points. `None` if any of `M`, `D_a`, `D_b` straddles 0.
fn cmp_filter(a: [[Iv; 4]; 3], b: [[Iv; 4]; 3], axis: usize) -> Option<Orient> {
    let (da, dva) = cramer_iv(a);
    let (db, dvb) = cramer_iv(b);
    let m = dva[axis].mul(db).sub(dvb[axis].mul(da));
    cmp_combine(m.sign(), da.sign(), db.sign())
}

/// The astro-float escalation for `cmp_coord`: the three signs at `prec` bits, and — when they do
/// not combine into an ordering — **how far apart the two coordinates may be** in the model's own
/// units.
///
/// Each implicit point's coordinate is `Dvec[axis]/D` (Cramer), so their difference is
/// `M / (D_a · D_b)` — a length, once divided. `M` alone is not: it carries both denominators, so
/// a threshold applied to it would move with the planes' scaling. The gap is `None` when either
/// denominator cannot be bounded away from zero, which is the near-parallel-planes case where
/// the implicit point itself is not well defined.
fn cmp_hp_with_gap(
    a: [[HpIv; 4]; 3],
    b: [[HpIv; 4]; 3],
    axis: usize,
    prec: usize,
) -> Result<Orient, Gap> {
    let (da, dva) = cramer_hp(&a, prec);
    let (db, dvb) = cramer_hp(&b, prec);
    let m = dva[axis]
        .mul(&db, prec)
        .sub(&dvb[axis].mul(&da, prec), prec);
    match cmp_combine(m.sign(), da.sign(), db.sign()) {
        Some(o) => Ok(o),
        None => Err(coord_gap(&m, &da, &db)),
    }
}

/// `|M| / |D_a · D_b|`, the separation `M` stands for — bounded above from `M`'s own radius, and
/// below by the denominators' lower bounds (dividing by an upper bound would understate the gap).
fn coord_gap(m: &HpIv, da: &HpIv, db: &HpIv) -> Gap {
    let (lo_a, lo_b) = match (denom_lo(da), denom_lo(db)) {
        (Ok(a), Ok(b)) => (a, b),
        // Both may be short; the deeper shortfall is the one that has to be covered.
        (Err(Gap::Unresolved { short: x }), Err(Gap::Unresolved { short: y })) => {
            return Gap::Unresolved { short: x.max(y) };
        }
        (Err(g), _) | (_, Err(g)) => return g,
    };
    match m.rad.over(lo_a.times(lo_b)) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

/// CIP indirect `cmp_coord`: the sign of `a[axis] − b[axis]` where `a`, `b` are the
/// implicit points at which each three-plane triple meets — each plane through three
/// rotated points. Interval filter → astro-float escalation. `Positive` = `a[axis] >
/// b[axis]`, `Negative` = `<`, and a zero that is either **proved** or **proved within the
/// coincidence limit** (unlike the exact `nacre_predicates::indirect_cmp_coord`, whose `0` means
/// exactly equal — see [`Decision`] for which of the two this was). The
/// two-implicit companion of [`indirect_orient3d_judge`]; boolean wiring is stage 3.
pub fn indirect_cmp_coord_judge(
    a: [(&Pt3, &Pt3, &Pt3); 3],
    b: [(&Pt3, &Pt3, &Pt3); 3],
    axis: usize,
    j: Standard,
) -> Decision {
    let iv = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
        [
            plane_iv(t[0].0, t[0].1, t[0].2),
            plane_iv(t[1].0, t[1].1, t[1].2),
            plane_iv(t[2].0, t[2].1, t[2].2),
        ]
    };
    if let Some(o) = cmp_filter(iv(a), iv(b), axis) {
        return Decision::Sign(o);
    }
    escalate(j, j.coincidence, |prec| {
        let hp = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
            [
                plane_hp(t[0].0, t[0].1, t[0].2, prec),
                plane_hp(t[1].0, t[1].1, t[1].2, prec),
                plane_hp(t[2].0, t[2].1, t[2].2, prec),
            ]
        };
        cmp_hp_with_gap(hp(a), hp(b), axis, prec)
    })
}

/// A definite `bool` sign (`true` = positive) as an [`Orient`].
fn orient_of(pos: bool) -> Orient {
    if pos {
        Orient::Positive
    } else {
        Orient::Negative
    }
}

/// The sign of `det[n_a; n_b; n_c]`, the determinant of the three planes' normals — the
/// toleranced twin of `nacre_geom`'s `plane_pair_dir_sign` (`sign((n_a × n_b) · n_c)`),
/// which decides how the line `a ∩ b` runs relative to plane `c`. This is exactly the
/// Cramer `D` of the three planes ([`cramer_iv`] / [`cramer_hp`]); the interval filter
/// resolves the easy cases, an ambiguous one escalates, and a below-floor `D` is
/// [`Orient::Zero`] (the planes' normals are coincident — degenerate).
///
/// **Winding-dependent** (unlike [`orient3d_judge`] / [`indirect_cmp_coord_judge`], which
/// are normal-orientation invariant): each plane's normal is `(p1−p0)×(p2−p0)`, so the
/// result follows each triple's point order. The caller must pass the three points in a
/// consistent order (the ops wrapper passes each face's outward-oriented `tri`).
pub fn dir_sign_judge(
    a: (&Pt3, &Pt3, &Pt3),
    b: (&Pt3, &Pt3, &Pt3),
    c: (&Pt3, &Pt3, &Pt3),
    j: Standard,
) -> Decision {
    // ★ **Not fed from `Judge`'s interval-plane cache, and that is measured.** Wiring it here
    // bought 1.007x on the crossing collector: this judgement's cost is almost entirely the
    // escalation below, because a *true* zero — the walls meeting `wc` in no point, 8.9% of the
    // pairs the collector tests — can never be settled by an interval filter. See docs/dev-log.md.
    let d = normals_det_iv([
        plane_iv(a.0, a.1, a.2),
        plane_iv(b.0, b.1, b.2),
        plane_iv(c.0, c.1, c.2),
    ]);
    if let Some(pos) = d.sign() {
        return Decision::Sign(orient_of(pos));
    }
    // **The one judgement whose limit is an angle.** Its question is about directions, so the
    // coincidence length has to be converted: a deviation of `coincidence` across a model of size
    // `scale` subtends `coincidence / scale`. A model too small to divide by leaves no angle to
    // compare against, so the judgement is degenerate rather than guessed.
    let Some(limit) = j.coincidence.over(j.scale) else {
        return Decision::Degenerate;
    };
    escalate(j, limit, |prec| {
        let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
        let planes = [ph(a), ph(b), ph(c)];
        let dh = normals_det_hp(&planes, prec);
        match dh.sign() {
            Some(pos) => Ok(orient_of(pos)),
            None => Err(dir_gap(&dh, &planes, prec)),
        }
    })
}

/// **How far three plane normals are from being coplanar — as an angle, not a length.**
///
/// The other judges reduce to a distance; this one cannot, because it asks a question about
/// *directions*: whether three planes share a common line direction. Its determinant is the
/// triple product of the normals, which carries each normal's own magnitude — and those are
/// arbitrary, since a plane's coefficients may be scaled freely. Dividing by `|n_a||n_b||n_c|`
/// leaves the triple product of the **unit** normals, which is dimensionless and, near zero, is
/// the sine of the angle by which the third normal misses the other two's plane.
///
/// A caller compares it against `coincidence_precision / scale`: the angle a deviation of the
/// coincidence limit subtends across the model. `None` when a normal cannot be bounded away from
/// zero — a degenerate plane has no direction to be off by.
fn dir_gap(d: &HpIv, planes: &[[HpIv; 4]; 3], prec: usize) -> Gap {
    let mut denom = Bound::of(1.0);
    for p in planes {
        // `|n| ≥ max|n_k|`, which is enough and needs no square root. One component clearing zero
        // is all a normal needs, so a shortfall only counts when **every** component fell short —
        // and then the smallest of them is the cheapest way out.
        let mut lo = Bound::ZERO;
        let mut short: Option<usize> = None;
        let mut vanished = 0;
        for c in p.iter().take(3) {
            match denom_lo(c) {
                Ok(b) if lo.lt(b) => lo = b,
                Ok(_) => {}
                Err(Gap::Unresolved { short: s }) => {
                    short = Some(short.map_or(s, |t: usize| t.min(s)));
                }
                Err(_) => vanished += 1,
            }
        }
        if lo.is_zero() {
            return match short {
                Some(short) => Gap::Unresolved { short },
                // Every component of this normal is exactly zero: the plane has no direction.
                None => {
                    debug_assert_eq!(vanished, 3);
                    Gap::Vanished
                }
            };
        }
        denom = denom.times(lo);
    }
    let _ = prec;
    match d.rad.over(denom) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The precision these fixtures judge at. Production chooses it per model
    /// ([`judge_precision`]); a fixture pins one so its expectations stay fixed.
    const FIXTURE_PREC: usize = 256;

    /// The judgement context those fixtures use: that precision, plus the coincidence limit the
    /// default derivation gives a unit-scale model (output resolution `2⁻⁵²`, two words further
    /// down) and the production cap.
    fn fixture() -> Standard {
        Standard {
            prec: FIXTURE_PREC,
            coincidence: Bound::pow2(-180),
            scale: Bound::of(1.0),
            cap: 4096,
        }
    }

    fn ri(n: i128, d: i128) -> Rat {
        Rat::new(n, d).unwrap()
    }

    fn deg(n: i128, d: i128) -> Angle {
        Angle::from_deg(ri(n, d)).unwrap()
    }

    /// **The climb lands where the shortfall says, in one jump.**
    ///
    /// The gap is `C · 2⁻ᵖʳᵉᶜ` over a cofactor and `C` does not move with the precision, so
    /// `log₂(gap / limit)` is exactly the number of missing bits — the same reading
    /// [`judge_precision`] takes to size the model. This pins that the loop uses it: the fixture's
    /// gap shrinks bit-for-bit with the precision, it starts 100 bits short, and the next attempt
    /// must arrive at 384 — `256 + 100` rounded up to a word — not at 512 (doubling), and not at
    /// 320 (a fixed word step, which would need two more rounds).
    #[test]
    fn the_climb_lands_where_the_shortfall_names_in_one_jump() {
        let mut asked = Vec::new();
        let j = Standard {
            prec: 256,
            coincidence: Bound::pow2(-300),
            scale: Bound::of(1.0),
            cap: 4096,
        };
        let out = escalate(j, j.coincidence, |prec| {
            asked.push(prec);
            // `C = 2⁵⁶`: at 256 bits the gap is `2⁻²⁰⁰`, a hundred bits above the limit.
            Err(Gap::Of(Bound::pow2(56 - prec as i64)))
        });
        assert_eq!(asked, vec![256, 384], "the climb overshot or crept");
        assert!(
            matches!(out, Decision::Coincident { within } if within.lt(Bound::pow2(-299))),
            "{out:?} — the second attempt was inside the limit and had to be reported as proof"
        );
    }

    /// **At the cap the judgement says so; it does not quietly become a zero.**
    ///
    /// This is the failure the whole cell exists to remove — an undecided determinant that reads
    /// as "these coincide" and merges two things that are apart. A gap that never reaches the
    /// limit must come back [`Decision::Exhausted`], and `orient()` may collapse it to `Zero` for
    /// the geometry only because the outcome itself is still there to be reported.
    #[test]
    fn a_gap_that_never_reaches_the_limit_is_exhausted_not_coincident() {
        let j = Standard {
            prec: 256,
            coincidence: Bound::pow2(-300),
            scale: Bound::of(1.0),
            cap: 512,
        };
        let mut rounds = 0;
        // A gap that ignores the precision entirely — a cofactor collapsing as fast as the
        // radius shrinks. No depth settles it, which is what the cap is for.
        let out = escalate(j, j.coincidence, |_| {
            rounds += 1;
            Err(Gap::Of(Bound::pow2(-10)))
        });
        // **The bound it did establish rides out with it.** Without it the caller cannot tell a
        // judgement that stopped `2⁻¹⁰` short from one that stopped `2⁻³⁰⁰` short, and those are
        // not the same news about the model.
        assert_eq!(
            out,
            Decision::Exhausted {
                at: 512,
                within: Some(Bound::pow2(-10)),
            }
        );
        assert_eq!(out.orient(), Orient::Zero, "the geometry still gets a sign");
        assert!(rounds > 1 && rounds < 10, "climbed {rounds} times");

        // An unresolved cofactor that never resolves has no distance to quote at all — reporting
        // one would be inventing it.
        let out = escalate(j, j.coincidence, |_| Err(Gap::Unresolved { short: 64 }));
        assert_eq!(
            out,
            Decision::Exhausted {
                at: 512,
                within: None,
            }
        );
    }

    /// **A degenerate witness is a different answer from an exhausted one, and must not climb.**
    ///
    /// No gap at all means the quantity that turns the determinant into a distance could not be
    /// bounded away from zero — a collapsed witness triangle, or three planes with no meeting
    /// point. More bits sharpen a distance; they do not conjure one. Telling the two apart is the
    /// whole reason the outcome is not a single "could not decide": one cause is a budget, the
    /// other is the geometry, and a caller that cannot tell them apart cannot act on either.
    #[test]
    fn a_degenerate_witness_is_told_apart_from_an_exhausted_one() {
        let j = Standard {
            prec: 256,
            coincidence: Bound::pow2(-300),
            scale: Bound::of(1.0),
            cap: 4096,
        };
        let mut rounds = 0;
        let out = escalate(j, j.coincidence, |_| {
            rounds += 1;
            Err(Gap::Vanished)
        });
        assert_eq!(out, Decision::Degenerate);
        assert_eq!(rounds, 1, "a vanished cofactor must not be re-realized");

        // …but a cofactor that is merely *unresolved* is the other cause, and it must climb —
        // reporting it as degenerate is how a judgement nobody looked at closely enough turns
        // into a silent merge. It names its own shortfall, so the climb is one jump here too.
        let mut asked = Vec::new();
        let out = escalate(j, j.coincidence, |prec| {
            asked.push(prec);
            if prec < 512 {
                Err(Gap::Unresolved { short: 200 })
            } else {
                Ok(Orient::Positive)
            }
        });
        assert_eq!(out, Decision::Sign(Orient::Positive));
        assert_eq!(
            asked,
            vec![256, 512],
            "an unresolved cofactor did not climb"
        );
    }

    /// **A structural zero comes back proved, not assumed.**
    ///
    /// Four points that are coplanar *before* a rotation stay coplanar after it, but the judge
    /// cannot see that: one base is `1/3`, so the pre-rotation shortcut cannot hand the question
    /// to the exact predicate and the determinant is a transcendental zero no finite precision
    /// separates. Before the gap was wired in, that came back `Orient::Zero` with nothing behind
    /// it. Now it comes back with the distance it established — and that distance has to be
    /// **below the coincidence limit**, which is what makes "these are the same plane" a proof
    /// rather than a shrug.
    ///
    /// The plane is **tilted** (`z = x`), which the first version of this fixture was not: a plane
    /// perpendicular to the rotation axis keeps its `z` coordinate exactly, so every error term
    /// meets an exactly-zero factor and the determinant comes back a *proved* zero instead — a
    /// real outcome (and one this loop reports), but not the one under test here.
    #[test]
    fn a_structural_zero_comes_back_proved_not_assumed() {
        // The plane `z = x` through four points, one of them irrational, all turned together.
        let pt = |x: Rat, y: Rat| {
            Pt3::at([x, y, x]).rotate_about(Axis::Z, deg(30, 1), [ri(1, 3), ri(1, 7), ri(0, 1)])
        };
        let (a, b, c, d) = (
            pt(ri(1, 3), ri(1, 1)),
            pt(ri(0, 1), ri(0, 1)),
            pt(ri(1, 1), ri(0, 1)),
            pt(ri(0, 1), ri(1, 1)),
        );
        let j = fixture();
        let out = orient3d_judge(&a, &b, &c, &d, j);
        let Decision::Coincident { within } = out else {
            panic!("{out:?} — a structural zero must be proved, not assumed");
        };
        assert!(
            !j.coincidence.lt(within),
            "reported {:?} against a limit of {:?}",
            within.exp2(),
            j.coincidence.exp2()
        );
        assert_eq!(out.orient(), Orient::Zero);

        // …and the limit is genuinely consulted: ask for a coincidence a thousand bits finer than
        // the model can carry and the same judgement must refuse to call it one.
        let strict = Standard {
            coincidence: Bound::pow2(-2000),
            cap: 512,
            ..j
        };
        let strict_out = orient3d_judge(&a, &b, &c, &d, strict);
        let Decision::Exhausted {
            at: 512,
            within: Some(w),
        } = strict_out
        else {
            panic!("{strict_out:?} — the coincidence limit was ignored, any zero would pass");
        };
        // What it *could* show is still the honest number, and it is nowhere near the limit asked
        // for: the report says "within this much", not "these are the same".
        assert!(strict.coincidence.lt(w), "reported {:?}", w.exp2());
    }

    /// **The design principle, as a test: raising the precision must actually narrow the bound.**
    ///
    /// The old declare-0 floor was a constant, so a deeper escalation bought a smaller threshold
    /// only because `2⁻ᵖʳᵉᶜ` appeared in it — the *estimate* it multiplied never improved. A
    /// computed radius has to do better than that, and if it did not, every "a deeper rung settles
    /// nothing" measurement would be reporting a broken radius rather than a fact about the
    /// geometry. So this pins the response directly: doubling the precision must take roughly a
    /// factor of `2⁻ᵖʳᵉᶜ` off the radius, both for a realized coordinate and for a determinant
    /// built out of several of them.
    #[test]
    fn a_deeper_precision_actually_narrows_the_radius() {
        // A rotation with an irrational cos/sin, about a non-origin pivot, so the radius has
        // every term in it: the trig bound, the pivot arithmetic, and the per-operation rounding.
        let pt = |x: i128, y: i128, z: i128| {
            Pt3::at([ri(x, 10), ri(y, 10), ri(z, 10)]).rotate_about(
                Axis::Z,
                deg(37, 1),
                [ri(1, 3), ri(1, 7), ri(0, 1)],
            )
        };
        let (a, b, c, d) = (pt(3, 5, 0), pt(11, 2, 4), pt(-7, 9, 13), pt(1, -6, 2));

        let mut prev_coord: Option<i64> = None;
        let mut prev_det: Option<i64> = None;
        for prec in [256usize, 512, 1024, 2048] {
            let coord = a.hp_coord(prec)[0]
                .rad
                .exp2()
                .expect("a rotated coordinate is not exact, so its radius is not zero");
            let det = det3_hp(&a, &b, &c, &d, prec)
                .rad
                .exp2()
                .expect("a determinant over rotated points carries a radius");
            if let (Some(pc), Some(pd)) = (prev_coord, prev_det) {
                // Each step doubles `prec`, so the radius should drop by about that many binary
                // orders. Demand most of it — the seed terms and the operation count add a
                // constant offset that does not shrink, but it must not dominate.
                let want = (prec / 2) as i64 - 32;
                assert!(
                    pc - coord >= want,
                    "coordinate radius went 2^{pc} → 2^{coord} at {prec} bits: only {} orders, \
                     wanted {want}. A radius that ignores precision makes the ladder meaningless.",
                    pc - coord
                );
                assert!(
                    pd - det >= want,
                    "determinant radius went 2^{pd} → 2^{det} at {prec} bits: only {} orders, \
                     wanted {want}",
                    pd - det
                );
            }
            prev_coord = Some(coord);
            prev_det = Some(det);
        }
    }

    /// **Is the normalized determinant actually a distance?**
    ///
    /// The coincidence limit is a length, so the quantity compared against it has to be one too.
    /// `orient3d` is a signed volume; dividing by the area term is supposed to leave a height.
    /// Nothing else in the suite would notice if that division were wrong by a factor of the
    /// triangle's size — the sign would still be right, and only the *threshold* would silently
    /// become a number that scales with the model.
    ///
    /// So this builds a point a **known** height above a plane and checks the normalized value is
    /// that height: across heights spanning 12 orders of magnitude, model scales spanning 9, and
    /// a triangle deliberately made 1000× larger, which is exactly what an unnormalized
    /// determinant would be fooled by.
    #[test]
    fn the_normalized_determinant_is_a_height_in_model_units() {
        let prec = 256;
        for scale in [1i128, 1_000, 1_000_000_000] {
            for (hn, hd) in [(1i128, 1i128), (1, 1_000), (1, 1_000_000_000_000)] {
                for tri_span in [1i128, 1_000] {
                    let p = |x: i128, y: i128, z: (i128, i128)| {
                        Pt3::at([ri(x * scale, 1), ri(y * scale, 1), ri(z.0, z.1)])
                    };
                    // Plane z = 0 through three points, and the query a height h above it.
                    let (b, c, d) = (
                        p(tri_span, 0, (0, 1)),
                        p(0, tri_span, (0, 1)),
                        p(0, 0, (0, 1)),
                    );
                    let a = p(0, 0, (hn, hd));
                    let det = det3_hp(&a, &b, &c, &d, prec);
                    let sub = |x: &HpIv, y: &HpIv| x.sub(y, prec);
                    let (bh, ch, dh) = (b.hp_coord(prec), c.hp_coord(prec), d.hp_coord(prec));
                    let (u, v) = (
                        [
                            sub(&bh[0], &dh[0]),
                            sub(&bh[1], &dh[1]),
                            sub(&bh[2], &dh[2]),
                        ],
                        [
                            sub(&ch[0], &dh[0]),
                            sub(&ch[1], &dh[1]),
                            sub(&ch[2], &dh[2]),
                        ],
                    );
                    let cross = [
                        u[1].mul(&v[2], prec).sub(&u[2].mul(&v[1], prec), prec),
                        u[2].mul(&v[0], prec).sub(&u[0].mul(&v[2], prec), prec),
                        u[0].mul(&v[1], prec).sub(&u[1].mul(&v[0], prec), prec),
                    ];
                    // Squared, so no square root is needed: `det² / |cross|²` must be `h²`.
                    let sq = |x: &BigFloat| x.mul(x, prec, HP_RM);
                    let norm2 = cross.iter().fold(BigFloat::from_f64(0.0, prec), |acc, k| {
                        acc.add(&sq(&k.mid), prec, HP_RM)
                    });
                    let got = sq(&det.mid).div(&norm2, prec, HP_RM);
                    let h = BigFloat::from_i128(hn, prec).div(
                        &BigFloat::from_i128(hd, prec),
                        prec,
                        HP_RM,
                    );
                    let err = rel_err(&got, &sq(&h), prec);
                    assert!(
                        err < 1e-6,
                        "scale {scale}, h {hn}/{hd}, triangle span {tri_span}: the normalized \
                         value is not the height (relative error {err:e})"
                    );
                }
            }
        }
    }

    /// **Is `cmp_coord`'s normalized value actually a coordinate difference?**
    ///
    /// The same trap as the orient3d height: `M` alone carries both Cramer denominators, so a
    /// threshold applied to it moves when the planes are scaled — and the sign, which is all the
    /// suite checks, does not move at all. Two implicit points a **known** distance apart, with
    /// the plane coefficients deliberately scaled by 1000 (which multiplies `M` by 10⁹ and must
    /// leave the gap untouched).
    #[test]
    fn the_normalized_cmp_is_a_coordinate_difference() {
        let prec = 256;
        for (gn, gd) in [(1i128, 1i128), (7, 100), (1, 1_000_000_000)] {
            for k in [1i128, 1_000] {
                // x = 0, y = 0, z = 0  and  x = 0, y = 0, z = g: two points on the z axis, `g`
                // apart. Coefficients scaled by `k` — the same planes, differently written.
                let pl = |c: [i128; 4]| {
                    [
                        HpIv::exact(BigFloat::from_i128(c[0] * k, prec)),
                        HpIv::exact(BigFloat::from_i128(c[1] * k, prec)),
                        HpIv::exact(BigFloat::from_i128(c[2] * k, prec)),
                        // `d` carries the offset, which must not be scaled away: `g = gn/gd`.
                        HpIv::new(
                            BigFloat::from_i128(c[3] * k, prec).div(
                                &BigFloat::from_i128(gd, prec),
                                prec,
                                HP_RM,
                            ),
                            Bound::ZERO,
                        ),
                    ]
                };
                let x0 = pl([1, 0, 0, 0]);
                let y0 = pl([0, 1, 0, 0]);
                let z0 = pl([0, 0, 1, 0]);
                let zg = pl([0, 0, 1, -gn]); // z = gn/gd
                let (da, dva) = cramer_hp(&[x0.clone(), y0.clone(), z0], prec);
                let (db, dvb) = cramer_hp(&[x0, y0, zg], prec);
                let m = dva[2].mul(&db, prec).sub(&dvb[2].mul(&da, prec), prec);
                // The midpoint version of `coord_gap`: |M| / |D_a·D_b| must be the separation.
                let got = m
                    .mid
                    .div(&da.mid.mul(&db.mid, prec, HP_RM), prec, HP_RM)
                    .abs();
                let want =
                    BigFloat::from_i128(gn, prec).div(&BigFloat::from_i128(gd, prec), prec, HP_RM);
                let err = rel_err(&got, &want, prec);
                assert!(
                    err < 1e-6,
                    "gap {gn}/{gd}, coefficients scaled by {k}: the normalized value is not the \
                     separation (relative error {err:e})"
                );
            }
        }
    }

    /// **Is the indirect orient3d's normalized value a point-to-plane distance?**
    ///
    /// This one carries *two* denominators — `|D|` because the implicit point is `Dvec/D`, and
    /// `|cross|` because the dot product carries the triangle's area — so there are two separate
    /// ways for it to stop being a length. The fixture scales each independently: the plane
    /// coefficients by `k` (which moves `D`) and the query triangle by `t` (which moves `cross`),
    /// while the true distance stays put.
    #[test]
    fn the_normalized_indirect_orient3d_is_a_distance() {
        let prec = 256;
        let big = |v: i128| HpIv::exact(BigFloat::from_i128(v, prec));
        for (hn, hd) in [(1i128, 1i128), (3, 100), (1, 1_000_000)] {
            for k in [1i128, 1_000] {
                for t in [1i128, 1_000] {
                    // V = ∩(x=0, y=0, z=h) sits `h` above the plane z = 0 through the triangle
                    // (0,0,0), (t,0,0), (0,t,0).
                    let pl = |c: [i128; 3], d: (i128, i128)| {
                        [
                            big(c[0] * k),
                            big(c[1] * k),
                            big(c[2] * k),
                            HpIv::new(
                                BigFloat::from_i128(d.0 * k, prec).div(
                                    &BigFloat::from_i128(d.1, prec),
                                    prec,
                                    HP_RM,
                                ),
                                Bound::ZERO,
                            ),
                        ]
                    };
                    let planes = [
                        pl([1, 0, 0], (0, 1)),
                        pl([0, 1, 0], (0, 1)),
                        pl([0, 0, 1], (-hn, hd)),
                    ];
                    let pt = |x: i128, y: i128| [big(x), big(y), big(0)];
                    let (d, m, _) = indirect_hp(planes, pt(t, 0), pt(0, t), pt(0, 0), prec);
                    // |M| / (|D| · |cross|); `cross` here is (t,0,0)×(0,t,0) = (0,0,t²).
                    // `cross` here is (t,0,0)×(0,t,0) = (0,0,t²), so `|cross| = t²`.
                    let norm = BigFloat::from_i128(t * t, prec);
                    let got = m.mid.div(&d.mid.mul(&norm, prec, HP_RM), prec, HP_RM).abs();
                    let want = BigFloat::from_i128(hn, prec).div(
                        &BigFloat::from_i128(hd, prec),
                        prec,
                        HP_RM,
                    );
                    let err = rel_err(&got, &want, prec);
                    assert!(
                        err < 1e-6,
                        "h {hn}/{hd}, coefficients x{k}, triangle x{t}: the normalized value is \
                         not the distance (relative error {err:e})"
                    );
                }
            }
        }
    }

    /// **Is `dir_sign`'s normalized value an angle?**
    ///
    /// This is the one judge whose answer is not a length: it asks whether three plane normals
    /// are coplanar, and its determinant carries each normal's magnitude — which is arbitrary,
    /// since a plane's coefficients scale freely. Normals `(1,0,0)`, `(0,1,0)`, `(0,1,ε)` miss
    /// coplanarity by `ε/√(1+ε²)` in sine; scaling any coefficient set must leave that alone,
    /// and without the normalization it does not.
    #[test]
    fn the_normalized_dir_sign_is_an_angle() {
        let prec = 256;
        let big = |v: i128| HpIv::exact(BigFloat::from_i128(v, prec));
        for (en, ed) in [(1i128, 1i128), (1, 1_000), (1, 1_000_000_000)] {
            for k in [1i128, 1_000] {
                let eps = HpIv::new(
                    BigFloat::from_i128(en * k, prec).div(
                        &BigFloat::from_i128(ed, prec),
                        prec,
                        HP_RM,
                    ),
                    Bound::ZERO,
                );
                let planes = [
                    [big(1), big(0), big(0), big(0)],
                    [big(0), big(1), big(0), big(0)],
                    [big(0), big(k), eps, big(0)],
                ];
                let (d, _) = cramer_hp(&planes, prec);
                // |D| / (|n_a||n_b||n_c|), with each |n| taken as its largest component.
                // Squared again, so the three norms need no square root: `D² / Π|n|²` is `sin²`.
                let sq = |x: &BigFloat| x.mul(x, prec, HP_RM);
                let norm2 = planes.iter().fold(BigFloat::from_f64(1.0, prec), |acc, p| {
                    let n2 = (0..3).fold(BigFloat::from_f64(0.0, prec), |a, i| {
                        a.add(&sq(&p[i].mid), prec, HP_RM)
                    });
                    acc.mul(&n2, prec, HP_RM)
                });
                let got = sq(&d.mid).div(&norm2, prec, HP_RM);
                let sine =
                    BigFloat::from_i128(en, prec).div(&BigFloat::from_i128(ed, prec), prec, HP_RM);
                // The exact sine of the miss is `ε/√(1+ε²)` — `ε` only for small `ε`, and the
                // fixture spans up to `ε = 1` where the two differ by √2. Squared: `ε²/(1+ε²)`.
                let s2 = sq(&sine);
                let want2 = s2.div(
                    &s2.add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
                    prec,
                    HP_RM,
                );
                let err = rel_err(&got, &want2, prec);
                assert!(
                    err < 1e-6,
                    "sine {en}/{ed}, coefficients x{k}: the normalized value is not the angle \
                     (relative error {err:e})"
                );
            }
        }
    }

    /// Relative error of `got` against `want`, as an `f64` magnitude — computed **in astro-float**
    /// so the comparison is not limited by `bf_mag`'s power-of-two rounding. Used by the
    /// normalization tests, where a factor-of-two slop would hide a real units error.
    fn rel_err(got: &BigFloat, want: &BigFloat, prec: usize) -> f64 {
        let d = got.sub(want, prec, HP_RM);
        bf_mag(&d.abs()) / bf_mag(&want.abs())
    }

    /// Deterministic PRNG (splitmix64) for reproducible stress corpora.
    fn splitmix64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn rng(state: &mut u64, lo: i128, hi: i128) -> i128 {
        lo + (u128::from(splitmix64(state)) % (hi - lo + 1) as u128) as i128
    }
    fn axis_of(k: i128) -> Axis {
        match k.rem_euclid(3) {
            0 => Axis::X,
            1 => Axis::Y,
            _ => Axis::Z,
        }
    }
    fn rand_base(st: &mut u64) -> [Rat; 3] {
        [
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
            ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ]
    }
    /// The f64 value's distance from the high-precision realization's **midpoint** — what the
    /// f64 tol has to bound. (The realization's own radius is a separate, far smaller quantity;
    /// `GT` is deep enough that it does not enter these comparisons.)
    fn abs_err(f: f64, truth: &HpIv, gt: usize) -> f64 {
        bf_mag(&BigFloat::from_f64(f, gt).sub(&truth.mid, gt, HP_RM).abs())
    }

    /// **The shared-motion shortcut must answer the same question a reflection is in the chain.**
    ///
    /// `shared_base` cancels one shared motion out of `det[a−d, b−d, c−d]` and answers on the
    /// pre-motion coordinates, *exactly* — on the strength of "a rigid motion preserves this
    /// determinant". A **reflection does not**: it negates it. And the failure is not conservative,
    /// because all four points carry the same chain, so the flip is uniform and the shortcut
    /// returns a confidently wrong `Sign`.
    ///
    /// So the claim is measured rather than argued: random chains that mix rotations, translations
    /// and reflections, against a 512-bit realization of the same four definitions.
    ///
    /// `#[ignore]`: astro-float ground truth is slow.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn the_shared_motion_shortcut_agrees_with_ground_truth() {
        const GT: usize = 512;
        let mut st = 0x51D2_7E44_0C13_A001u64;
        let (mut took_shortcut, mut mirror_seen) = (0usize, false);
        for _ in 0..2000 {
            // Four points sharing one chain — the shortcut's precondition. **Dyadic** bases, or
            // `shared_base` declines before the sign is ever in question (its own guard is that a
            // base must round-trip through f64) and the test would pass vacuously.
            let dyadic = |st: &mut u64| -> [Rat; 3] {
                std::array::from_fn(|_| ri(rng(st, -100_000, 100_000), 1 << rng(st, 0, 6)))
            };
            let bases: [[Rat; 3]; 4] = std::array::from_fn(|_| dyadic(&mut st));
            let mut pts: Vec<Pt3> = bases.iter().map(|&b| Pt3::at(b)).collect();
            for _ in 0..rng(&mut st, 1, 4) {
                match rng(&mut st, 0, 2) {
                    0 => {
                        let off = dyadic(&mut st);
                        pts = pts.into_iter().map(|p| p.translate(off)).collect();
                    }
                    1 => {
                        mirror_seen = true;
                        let ax = axis_of(rng(&mut st, 0, 2));
                        let off = if rng(&mut st, 0, 2) == 0 {
                            ri(0, 1)
                        } else {
                            dyadic(&mut st)[0]
                        };
                        pts = pts.into_iter().map(|p| p.mirror(ax, off)).collect();
                    }
                    _ => {
                        let ax = axis_of(rng(&mut st, 0, 2));
                        let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 997));
                        let piv = dyadic(&mut st);
                        pts = pts
                            .into_iter()
                            .map(|p| p.rotate_about(ax, ang, piv))
                            .collect();
                    }
                }
            }
            let [pa, pb, pc, pd] = [&pts[0], &pts[1], &pts[2], &pts[3]];
            let Some(b) = shared_base(&[pa, pb, pc, pd]) else {
                continue; // a base that does not round-trip; the toleranced path answers
            };
            took_shortcut += 1;
            let shortcut = match nacre_predicates::orient3d(b[0], b[1], b[2], b[3]) {
                x if x > 0.0 => 1i8,
                x if x < 0.0 => -1,
                _ => 0,
            };
            // Ground truth: the same determinant, realized from the definitions at 512 bits. Deep
            // enough that only a genuinely near-zero case is ambiguous, and those are skipped.
            let det = det3_hp(pa, pb, pc, pd, GT);
            let Some(pos) = det.sign() else {
                continue;
            };
            let truth = if pos { 1i8 } else { -1 };
            assert_eq!(
                shortcut, truth,
                "the shortcut answered a different question than the definitions do"
            );
        }
        assert!(mirror_seen, "the sample must include a reflection");
        assert!(
            took_shortcut > 100,
            "the shortcut must actually fire ({took_shortcut} times)"
        );
    }

    /// Soundness: over random **motion** chains (mixed axes, arbitrary pivots, exact and inexact
    /// angles, **and rational translations and reflections interleaved**), the direction-wise tol
    /// must bound the true f64 error on every axis (astro-float 512-bit ground truth) — the
    /// production mirror of exact3d H-d/H-f. (An `err` below 1e-100 is 512-bit GT noise.)
    ///
    /// **Every node kind must appear here.** Each arrives with its own tol term, and nothing else
    /// in the suite checks that term is an upper bound — the whole judgment layer is sound only if
    /// it is. (The translate term was once added without this, and had to be back-filled.)
    /// `#[ignore]`: astro-float ground truth is slow; run with `--ignored` (+ CI). The
    /// per-predicate statistical validation is the `h_*` tests beside this one.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn tol_bounds_error_over_random_chains() {
        const GT: usize = 512;
        let mut st = 0x2A5C_1234_ABCD_9999u64;
        let (mut exact_seen, mut pivot_seen) = (false, false);
        let (mut translate_seen, mut mirror_seen) = (false, false);
        let mut worst_ratio = 0.0_f64;
        let mut worst_rot = 0.0_f64;
        for _ in 0..2000 {
            let base = rand_base(&mut st);
            let seed_free = Pt3::at(base).tol == [0.0; 3];
            let mut p = Pt3::at(base);
            // **From one node, not two.** A single origin rotation of an exact base is the case
            // the deleted 2D frame validated on its own; sampling it here is what makes this test
            // strictly cover that one, rather than merely resemble it.
            for _ in 0..rng(&mut st, 1, 5) {
                // Links of every kind, in both orders — which is where a tol term that only holds
                // "on its own" would show.
                match rng(&mut st, 0, 3) {
                    0 => {
                        translate_seen = true;
                        p = p.translate(rand_base(&mut st));
                        continue;
                    }
                    1 => {
                        mirror_seen = true;
                        let off = if rng(&mut st, 0, 2) == 0 {
                            Rat::from_int(0) // the exact case: a pure sign flip
                        } else {
                            rand_base(&mut st)[0]
                        };
                        p = p.mirror(axis_of(rng(&mut st, 0, 2)), off);
                        continue;
                    }
                    _ => {}
                }
                let ax = axis_of(rng(&mut st, 0, 2));
                let pivot = match rng(&mut st, 0, 2) {
                    0 => [Rat::from_int(0); 3],
                    1 => {
                        pivot_seen = true;
                        rand_base(&mut st)
                    }
                    _ => {
                        pivot_seen = true;
                        base // point on the pivot
                    }
                };
                let ang = if rng(&mut st, 0, 3) == 0 {
                    exact_seen = true;
                    deg(90 * rng(&mut st, 0, 3), 1)
                } else {
                    deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973))
                };
                p = p.rotate_about(ax, ang, pivot);
            }
            let hp = p.hp_coord(GT);
            for (axis, hp_a) in hp.iter().enumerate() {
                let err = abs_err(p.coord[axis], hp_a, GT);
                assert!(
                    err <= p.tol[axis] || err < 1e-100,
                    "tol must bound the error: axis {axis}, err {err:e} > tol {:e}",
                    p.tol[axis]
                );
                if p.tol[axis] > 0.0 {
                    worst_ratio = worst_ratio.max(err / p.tol[axis]);
                    if seed_free {
                        worst_rot = worst_rot.max(err / p.tol[axis]);
                    }
                }
            }
        }
        assert!(
            exact_seen && pivot_seen && translate_seen && mirror_seen,
            "corpus must mix exact angles, pivots, translations and reflections"
        );
        // **How much of the bound the real error actually uses.** `DA_F64` is the one constant in
        // this kernel that cannot be derived — `f64::cos`'s accuracy is not contracted by Rust or
        // by any libm this runs on — so it rests on measurement, and this is the measurement.
        eprintln!(
            "[tol tightness] worst err/tol: {worst_ratio:.4} all, {worst_rot:.4} rotation-only"
        );
        assert!(
            worst_ratio < 1.0,
            "the bound was reached exactly, which leaves nothing for a platform whose trig is \
             one ulp worse than this one's"
        );
    }

    /// A 90°-family chain about the origin stays exactly tol 0 (Niven exact realization).
    /// **Why the tol is a coordinate *mixing* term and not a tangential one.**
    ///
    /// The design's first directional tol was tangential — `tol_x = da·|y|`, `tol_y = da·|x|` —
    /// on the reasoning that rotation moves a point along its circle. It is **not sound**: `cos`
    /// and `sin` round independently, which puts a *radial* component into the error, and near an
    /// axis the tangential prediction goes to zero while the real error does not. The adopted
    /// `(|bx|+|by|)·da` per component covers both.
    ///
    /// This is the standing guard on that shape. It moved here when the 2D frame it was written
    /// in was deleted (nothing consumed that frame), because the *invariant* is about the tol
    /// formula, which 3D uses verbatim — a rotation about Z is the 2D case with `z` carried along.
    #[test]
    fn the_tangential_tol_is_unsound_and_the_mixing_one_is_not() {
        const GT: usize = 256;
        let (mut tangential_sound, mut mixing_sound) = (true, true);
        let mut first_break = None;
        for (bxn, bxd, byn, byd) in [(1, 1, 0, 1), (3, 1, 4, 1), (5, 2, 7, 3), (1, 1, 1, 1)] {
            // Near-axis angles are where the tangential model is thinnest, so they are in here.
            for (an, ad) in [
                (30, 1),
                (45, 1),
                (60, 1),
                (37, 1),
                (899, 10),
                (1, 10),
                (3, 1),
            ] {
                let base = [ri(bxn, bxd), ri(byn, byd), ri(0, 1)];
                let p = Pt3::at(base).rotate(Axis::Z, deg(an, ad));
                let hp = p.compute_hp(GT);
                for (k, hp_k) in hp.iter().enumerate().take(2) {
                    let err = abs_err(p.coord[k], hp_k, GT);
                    // (1) the refuted tangential prediction: the *other* coordinate's magnitude.
                    if err > DA_F64 * p.coord[1 - k].abs() {
                        tangential_sound = false;
                        first_break.get_or_insert((bxn, bxd, byn, byd, an, ad, k));
                    }
                    // (2) the adopted mixing tol, which is what `rotate_about` computes.
                    if err > p.tol[k] {
                        mixing_sound = false;
                    }
                }
            }
        }
        assert!(
            !tangential_sound,
            "the tangential formula bounded every sample — either the fixtures stopped reaching \
             near an axis, or the realization changed and this guard is now vacuous"
        );
        assert!(
            mixing_sound,
            "the adopted tol failed to bound a single origin rotation: {first_break:?}"
        );
    }

    #[test]
    fn quadrantal_origin_chain_is_tol_zero() {
        let p = Pt3::at([ri(3, 1), ri(5, 1), ri(7, 1)])
            .rotate(Axis::Z, deg(90, 1))
            .rotate(Axis::X, deg(180, 1))
            .rotate(Axis::Y, deg(270, 1));
        assert_eq!(p.tol, [0.0; 3]);
        // and the realized coords are the exact permutation/negation (no spurious term).
        assert_eq!(p.coord, [7.0, -3.0, -5.0]);
    }

    /// **`Pt3::exact` states the tol that `Pt3::at` would measure — the same value.**
    ///
    /// This equivalence is what licenses the substitution on the boolean's hot path, where `at`
    /// spent nine 120-bit BigFloat operations per call to arrive at zero. It rests on three
    /// links, and the third is the one worth a test: `try_from_f64` represents an f64 exactly,
    /// realizing that base back yields the same f64, and `bf_mag` of an exact zero is `0.0` (not
    /// a floor). If any link broke, `at` would report a nonzero tol here and the fast
    /// constructor would be silently changing geometry rather than skipping arithmetic.
    #[test]
    fn exact_states_the_tol_that_at_would_measure() {
        for c in [
            [0.0, 1.0, -1.0],         // integers, both signs
            [0.5, 0.25, -0.125],      // dyadic fractions
            [4.0, 0.2, 3.0],          // 0.2 is not decimal-exact but *is* an exact f64
            [1e-8, -2.5e-9, 7.5e-7],  // small
            [1e18, -4e17, 3.5e19],    // large, still inside Rat's exponent range
            [1e-20, -2.5e-21, 5e-18], // small, still inside it (the floor is 2^-74 ≈ 5.3e-23)
        ] {
            let fast = Pt3::exact(c).expect("representable");
            let measured = Pt3::at(c.map(|x| Rat::try_from_f64(x).expect("representable")));
            assert_eq!(fast.coord, c, "the coordinates round-trip: {c:?}");
            assert_eq!(fast.coord, measured.coord, "same coord for {c:?}");
            assert_eq!(
                measured.tol, [0.0; 3],
                "`at` must measure exactly zero for an f64-derived base: {c:?}"
            );
            assert_eq!(fast.tol, measured.tol, "same tol for {c:?}");
        }
    }

    /// Outside `Rat`'s exponent range there is no exact base, and `exact` says so instead of
    /// panicking — the caller decides (an operation turns it into a named reject).
    ///
    /// **The range is much narrower than f64's, at both ends** — easy to get wrong, and I did
    /// on the first attempt. `Rat` is `Ratio<i128>`, so `mantissa · 2^exp` must fit `i128`
    /// (`|x| ≲ 1.7e38`) *and* the denominator `2^k` must (`k = 1075 − exp_field ≤ 126`, i.e.
    /// `|x| ≳ 2^-74 ≈ 5.3e-23`). Exact zero is special-cased and always representable.
    /// A CAD model at either extreme is not real; the limit is.
    #[test]
    fn exact_declines_a_coordinate_it_cannot_represent() {
        assert!(Pt3::exact([1e300, 0.0, 0.0]).is_none(), "too large");
        assert!(
            Pt3::exact([1e-30, 0.0, 0.0]).is_none(),
            "below the 2^-74 floor"
        );
        assert!(
            Pt3::exact([f64::MIN_POSITIVE, 0.0, 0.0]).is_none(),
            "subnormal"
        );
        // Zero is not a boundary case — it is special-cased and exact.
        assert_eq!(Pt3::exact([0.0; 3]).expect("zero is exact").tol, [0.0; 3]);
    }

    /// `at` seeds the base→f64 rounding: 0 for an integer base, positive for a base
    /// that is not f64-representable (e.g. 1/3).
    #[test]
    fn at_seeds_base_rounding_tol() {
        assert_eq!(Pt3::at([ri(2, 1), ri(3, 1), ri(4, 1)]).tol, [0.0; 3]);
        let third = Pt3::at([ri(1, 3), ri(0, 1), ri(0, 1)]);
        assert!(third.tol[0] > 0.0 && third.tol[1] == 0.0);
    }

    /// `at_with_tol` seeds a nonzero root tol (a Discovered seam) and the chain
    /// transports it (`|R|·old`): a 90° rotation about Z swaps the x/y tol components.
    #[test]
    fn seeded_tol_transports_through_rotation() {
        let p = Pt3::at_with_tol([ri(1, 1), ri(0, 1), ri(0, 1)], [1e-9, 2e-9, 3e-9])
            .rotate(Axis::Z, deg(90, 1));
        // 90° about Z: |R| swaps x,y → tol[0]=old tol[1], tol[1]=old tol[0]; z unchanged.
        assert_eq!(p.tol, [2e-9, 1e-9, 3e-9]);
    }

    /// Same-axis bundling is tighter than incremental (H-e): K steps of θ amplify the
    /// tol (each transported by |c|+|s| ≥ 1) vs one step of Kθ.
    #[test]
    fn bundling_is_tighter_than_incremental() {
        let base = [ri(11, 1), ri(-7, 1), ri(4, 1)];
        let (k, theta) = (30i128, ri(1, 7));
        let mut incr = Pt3::at(base);
        for _ in 0..k {
            incr = incr.rotate(Axis::Z, Angle::from_deg(theta).unwrap());
        }
        let k_theta = theta.checked_mul(Rat::from_int(k)).unwrap();
        let bundled = Pt3::at(base).rotate(Axis::Z, Angle::from_deg(k_theta).unwrap());
        assert!(
            bundled.tol[0] < incr.tol[0] && bundled.tol[1] < incr.tol[1],
            "bundled {:?} must be tighter than incremental {:?}",
            bundled.tol,
            incr.tol
        );
    }

    // ---- orient3d_judge (2a-ii) ----

    /// A point rotated about one random axis (through a random pivot) by one **inexact**
    /// rational angle — generic non-degenerate, heterogeneous provenance (each point its
    /// own rotation), the H-a corpus shape. Inexact angles only (an exact/near-coplanar
    /// mix would expose the tol-0 GT noise floor — that is the declare-0 test's job).
    fn rand_point(st: &mut u64) -> Pt3 {
        let base = rand_base(st);
        let axis = axis_of(rng(st, 0, 2));
        let ang = deg(rng(st, 0, 360_000), rng(st, 1, 9973)); // inexact
        let pivot = if rng(st, 0, 1) == 0 {
            [Rat::from_int(0); 3]
        } else {
            rand_base(st) // non-origin pivot: exercises the pivot tol → det3_bound path
        };
        Pt3::at(base).rotate_about(axis, ang, pivot)
    }

    /// H-a — `det3_bound` soundness: over many random heterogeneous-rotation 4-point
    /// configs (inexact angles, arbitrary pivots), the bound must upper-bound the real
    /// error of the f64 determinant vs the astro-float truth — never exceeded. The
    /// production mirror of exact3d H-a. `#[ignore]`: slow (4× astro-float per sample).
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_a_det3_bound_soundness() {
        const GT: usize = 512;
        let mut st = 0x3D00_1234_ABCD_EF01u64;
        let (mut bad, mut worst) = (0usize, 0.0_f64);
        for _ in 0..5000 {
            let (pa, pb, pc, pd) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let det = det3_f64(pa.coord, pb.coord, pc.coord, pd.coord);
            let truth = det3_hp(&pa, &pb, &pc, &pd, GT);
            let err = abs_err(det, &truth, GT);
            let bound = det3_bound(
                [pa.coord, pb.coord, pc.coord, pd.coord],
                [pa.tol, pb.tol, pc.tol, pd.tol],
            );
            if err > bound {
                bad += 1;
            }
            if bound > 0.0 {
                worst = worst.max(err / bound);
            }
        }
        eprintln!("[H-a N=5000] det3_bound violations: {bad}; worst tightness: {worst:.3}");
        assert_eq!(
            bad, 0,
            "det3_bound must bound the determinant error on every sample"
        );
    }

    /// Fast port check (default suite): a handful of random configs must not violate
    /// `det3_bound` (catches a transcription bug; the full statistical soundness is the
    /// `#[ignore]`d H-a + exact3d).
    #[test]
    fn det3_bound_port_check() {
        const GT: usize = 512;
        let mut st = 0x0BAD_F00D_1234_5678u64;
        for _ in 0..40 {
            let (pa, pb, pc, pd) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let det = det3_f64(pa.coord, pb.coord, pc.coord, pd.coord);
            let err = abs_err(det, &det3_hp(&pa, &pb, &pc, &pd, GT), GT);
            let bound = det3_bound(
                [pa.coord, pb.coord, pc.coord, pd.coord],
                [pa.tol, pb.tol, pc.tol, pd.tol],
            );
            assert!(err <= bound, "det3_bound {bound:e} < err {err:e}");
        }
    }

    /// A known tetrahedron judges `Positive`, its mirror `Negative`, and the sign is
    /// preserved under a shared rotation (rotating all four points does not flip
    /// orientation). Exact (tol 0) points take the filter's fast path.
    #[test]
    fn orient3d_sanity_and_rotation_invariance() {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        // d at origin; a,b,c along +x,+y,+z → right-handed → Positive.
        let (a, b, c, d) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1), pt(0, 0, 0));
        assert_eq!(
            orient3d_judge(&a, &b, &c, &d, fixture()).orient(),
            Orient::Positive
        );
        assert_eq!(
            orient3d_judge(&b, &a, &c, &d, fixture()).orient(),
            Orient::Negative,
            "swap → mirror"
        );

        // Shared rotation of all four (37° about Z through a rational pivot) keeps the sign.
        let rot = |p: &Pt3| {
            p.clone()
                .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-3, 1), ri(0, 1)])
        };
        assert_eq!(
            orient3d_judge(&rot(&a), &rot(&b), &rot(&c), &rot(&d), fixture()).orient(),
            Orient::Positive,
            "orientation is rotation-invariant"
        );
    }

    /// Declare-0: four **exactly coplanar** rational points (tol 0) → the determinant is
    /// exactly 0 → `Orient::Zero` (the ask-the-user case).
    #[test]
    fn orient3d_coplanar_is_zero() {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        // all four in the plane z = 0.
        let (a, b, c, d) = (pt(0, 0, 0), pt(3, 0, 0), pt(0, 5, 0), pt(2, 7, 0));
        assert_eq!(
            orient3d_judge(&a, &b, &c, &d, fixture()).orient(),
            Orient::Zero
        );
    }

    // ---- indirect orient3d (2c-i) ----

    /// `Iv` arithmetic keeps a sound radius (worst-case interval), and `sign` is
    /// definite exactly when the interval clears 0.
    #[test]
    fn iv_arithmetic_is_sound_and_sign_decides() {
        let (a, b) = (Iv::new(3.0, 0.1), Iv::new(-2.0, 0.2));
        // sub radius ≥ sum of input radii.
        assert!(a.sub(b).rad >= 0.1 + 0.2);
        // mul radius ≥ |mid_a|·rad_b + |mid_b|·rad_a (+ rad·rad + rounding).
        assert!(a.mul(b).rad >= 3.0 * 0.2 + 2.0 * 0.1);
        assert_eq!(Iv::new(1.0, 0.5).sign(), Some(true));
        assert_eq!(Iv::new(-1.0, 0.5).sign(), Some(false));
        assert_eq!(Iv::new(0.3, 0.5).sign(), None); // straddles 0 → escalate
    }

    /// Sanity: three coordinate planes meet at the origin `V=(0,0,0)`; `orient3d(V,q,r,s)`
    /// judges a known sign, flips on a q/r swap, is `Zero` when `s` is coplanar with
    /// `V,q,r`, and is invariant under a shared rotation (rotations preserve orientation).
    #[test]
    fn indirect_sanity_and_rotation_invariance() {
        fn tri(t: &(Pt3, Pt3, Pt3)) -> (&Pt3, &Pt3, &Pt3) {
            (&t.0, &t.1, &t.2)
        }
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        let pa = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 1, 0)); // z = 0
        let pb = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 0, 1)); // y = 0
        let pc = (pt(0, 0, 0), pt(0, 1, 0), pt(0, 0, 1)); // x = 0  → V = (0,0,0)
        let (q, r, s) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1));
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s, fixture()).orient(),
            Orient::Negative
        );
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &r, &q, &s, fixture()).orient(),
            Orient::Positive,
            "q/r swap flips the sign"
        );
        // s coplanar with V,q,r (all z = 0) → orient exactly 0 → declare-0.
        let s0 = pt(1, 1, 0);
        assert_eq!(
            indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s0, fixture()).orient(),
            Orient::Zero
        );
        // Shared rotation of all twelve points keeps the definite sign.
        let rot = |p: &Pt3| {
            p.clone()
                .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-1, 1), ri(0, 1)])
        };
        let rp = |t: (&Pt3, &Pt3, &Pt3)| (rot(t.0), rot(t.1), rot(t.2));
        let (ra, rb, rc) = (rp(tri(&pa)), rp(tri(&pb)), rp(tri(&pc)));
        let (rq, rr, rs) = (rot(&q), rot(&r), rot(&s));
        assert_eq!(
            indirect_orient3d_judge(tri(&ra), tri(&rb), tri(&rc), &rq, &rr, &rs, fixture())
                .orient(),
            Orient::Negative,
            "indirect orient is rotation-invariant"
        );
    }

    /// The high-precision indirect truth with a stability flag: `None` if even the
    /// ground truth cannot resolve `D` or `M` above its floor (a genuine degeneracy).
    #[allow(clippy::too_many_arguments)]
    /// Ground truth for the indirect judge — and **not** by running the judge harder.
    ///
    /// It used to call `sign_with_floor` with the judge's own `mag`, differing only in precision.
    /// A floor that is wrong is then wrong identically in both, so the two agree and the test
    /// passes: the oracle shared the defect it existed to find, which is how the `mag`-collapse
    /// bug survived a soundness suite. Here the verdict comes from **comparing two precisions**
    /// instead: a sign is trusted only when `prec` and `2·prec` produce the same nonzero sign, and
    /// anything else is `None` (not asserted against). That borrows no formula from the judge.
    fn indirect_truth(
        pa: (&Pt3, &Pt3, &Pt3),
        pb: (&Pt3, &Pt3, &Pt3),
        pc: (&Pt3, &Pt3, &Pt3),
        q: &Pt3,
        r: &Pt3,
        s: &Pt3,
        prec: usize,
    ) -> Option<Orient> {
        let at = |prec: usize| -> Orient {
            let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
            let (d, m, _gap) = indirect_hp(
                [ph(pa), ph(pb), ph(pc)],
                q.hp_coord(prec),
                r.hp_coord(prec),
                s.hp_coord(prec),
                prec,
            );
            // Raw signs, no floor: the agreement between two precisions is what filters noise.
            let raw = |x: &BigFloat| {
                if x.is_zero() {
                    None
                } else {
                    Some(x.is_positive())
                }
            };
            combine(raw(&d.mid), raw(&m.mid)).unwrap_or(Orient::Zero)
        };
        let (lo, hi) = (at(prec), at(2 * prec));
        (lo == hi && lo != Orient::Zero).then_some(lo)
    }

    /// H-b — rotated plane coefficient tol soundness. A plane's four coefficients
    /// derive from three rotated points (subtraction/cross/dot); the interval `rad` on
    /// each must upper-bound the real f64 error vs the astro-float truth. Corpus mixes
    /// origin and arbitrary pivots (`rand_point`). `#[ignore]`: slow astro-float GT.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_b_plane_coefficient_tol_soundness() {
        const GT: usize = 512;
        const N: usize = 10_000;
        let mut st = 0xB0B0_5555_1111_2222u64;
        let (mut bad, mut worst) = (0usize, 0.0_f64);
        for _ in 0..N {
            let (p0, p1, p2) = (
                rand_point(&mut st),
                rand_point(&mut st),
                rand_point(&mut st),
            );
            let iv = plane_iv(&p0, &p1, &p2);
            let hp = plane_hp(&p0, &p1, &p2, GT);
            for k in 0..4 {
                let err = abs_err(iv[k].mid, &hp[k], GT);
                if err > iv[k].rad {
                    bad += 1;
                }
                if iv[k].rad > 0.0 {
                    worst = worst.max(err / iv[k].rad);
                }
            }
        }
        eprintln!("[H-b N={N}] coefficient tol violations: {bad}; worst tightness: {worst:.3}");
        assert_eq!(bad, 0, "plane coefficient tol must bound the error");
    }

    /// **Asking one question four ways must give one answer.**
    ///
    /// Four planes through a common point are concurrent or they are not, and that fact does not
    /// depend on which three of them you call "the point" and which one you call "the query". The
    /// judge is free to *abstain* — but not to say `Positive` for one splitting and `Zero` for
    /// another, because then two callers reading the same geometry disagree, and the arrangement
    /// enters one vertex twice.
    ///
    /// **The corpus is built, not sampled.** Random points never place a query vertex *on* the
    /// implicit point, and that coincidence is the whole difficulty: it makes `row1 = D·(V − s)`
    /// cancel to nothing. Here every plane is spanned by the shared point `p` and two others, so
    /// `p` is both the meet of any three and a defining point of the fourth — the exact shape the
    /// engine hits when a tool edge lands in a target plane.
    #[test]
    fn one_concurrency_read_four_ways_gives_one_answer() {
        let mut st = 0x5EED_1234_ABCD_9999u64;
        let mut checked = 0usize;
        for _ in 0..200 {
            // A rotated common point, and four planes each spanned by it and two more.
            let axis = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 1, 359_000), rng(&mut st, 1, 997)); // inexact ⇒ toleranced
            let pivot = rand_base(&mut st);
            let turn = |b: [Rat; 3], st: &mut u64| {
                let _ = st;
                Pt3::at(b).rotate_about(axis, ang, pivot)
            };
            let p = turn(rand_base(&mut st), &mut st);
            let spans: Vec<[Pt3; 2]> = (0..4)
                .map(|_| {
                    [
                        turn(rand_base(&mut st), &mut st),
                        turn(rand_base(&mut st), &mut st),
                    ]
                })
                .collect();
            let plane = |i: usize| (&p, &spans[i][0], &spans[i][1]);

            // Every way of choosing which three planes make the point and which one is queried.
            let mut verdicts = Vec::new();
            for q in 0..4 {
                let tri: Vec<usize> = (0..4).filter(|&i| i != q).collect();
                // Skip a splitting whose three planes do not meet at a single point at all.
                if dir_sign_judge(plane(tri[0]), plane(tri[1]), plane(tri[2]), fixture()).orient()
                    == Orient::Zero
                {
                    continue;
                }
                // `p` goes in the **third** slot: `indirect_hp` forms `row1 = Dvec − D·s` from
                // that one, so putting the shared point there is what makes the subtraction cancel
                // — the configuration the engine actually hit. With `p` first the term stays
                // healthy and the corpus measures nothing.
                let t = plane(q);
                verdicts.push(
                    indirect_orient3d_judge(
                        plane(tri[0]),
                        plane(tri[1]),
                        plane(tri[2]),
                        t.1,
                        t.2,
                        t.0,
                        fixture(),
                    )
                    .orient(),
                );
            }
            if verdicts.len() < 2 {
                continue; // nothing to compare
            }
            checked += 1;
            assert!(
                verdicts.iter().all(|v| *v == verdicts[0]),
                "the same concurrency read four ways: {verdicts:?}"
            );
        }
        // A sweep that quietly stops constructing anything reads as agreement.
        assert!(
            checked > 20,
            "only {checked} configurations were comparable"
        );
        eprintln!("[form-invariance] {checked} concurrent configurations, all splittings agreed");
    }

    /// H-c — indirect orient3d soundness over heterogeneous provenance (corpus A) and a
    /// near-coplanar escalation-forcing corpus (corpus B). The interval filter →
    /// astro-float judge must never claim a sign opposite to a GT-stable 512-bit truth.
    /// Asserts both paths are live: `escalated > 0` (filter defers) **and**
    /// `filter_resolved > 0` (fast path resolves — a filter stuck at `None` would pass
    /// wrong-sign 0 vacuously). `#[ignore]`: slow astro-float GT.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_c_indirect_orient3d_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut check = |a: (&Pt3, &Pt3, &Pt3),
                         b: (&Pt3, &Pt3, &Pt3),
                         c: (&Pt3, &Pt3, &Pt3),
                         q: &Pt3,
                         r: &Pt3,
                         s: &Pt3| {
            let judged = indirect_orient3d_judge(a, b, c, q, r, s, fixture()).orient();
            let planes = [
                plane_iv(a.0, a.1, a.2),
                plane_iv(b.0, b.1, b.2),
                plane_iv(c.0, c.1, c.2),
            ];
            if indirect_filter(planes, pt_iv(q), pt_iv(r), pt_iv(s)).is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match indirect_truth(a, b, c, q, r, s, GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        };

        // Corpus A — heterogeneous provenance (each of the twelve points its own
        // rotation, origin or arbitrary pivot). Generic → mostly filter-resolved.
        let mut st = 0xC0C0_9999_ABAB_CDCDu64;
        for _ in 0..2000 {
            let p: Vec<Pt3> = (0..12).map(|_| rand_point(&mut st)).collect();
            check(
                (&p[0], &p[1], &p[2]),
                (&p[3], &p[4], &p[5]),
                (&p[6], &p[7], &p[8]),
                &p[9],
                &p[10],
                &p[11],
            );
        }

        // Corpus B — one shared rotation about a shared pivot (affine → coplanarity
        // preserved), near-coplanar so the escalation's sign-resolution is exercised.
        let add3 = |a: [Rat; 3], b: [Rat; 3]| {
            [
                a[0].checked_add(b[0]).unwrap(),
                a[1].checked_add(b[1]).unwrap(),
                a[2].checked_add(b[2]).unwrap(),
            ]
        };
        let smul = |k: Rat, a: [Rat; 3]| {
            [
                a[0].checked_mul(k).unwrap(),
                a[1].checked_mul(k).unwrap(),
                a[2].checked_mul(k).unwrap(),
            ]
        };
        let cross = |a: [Rat; 3], b: [Rat; 3]| {
            [
                a[1].checked_mul(b[2])
                    .unwrap()
                    .checked_sub(a[2].checked_mul(b[1]).unwrap())
                    .unwrap(),
                a[2].checked_mul(b[0])
                    .unwrap()
                    .checked_sub(a[0].checked_mul(b[2]).unwrap())
                    .unwrap(),
                a[0].checked_mul(b[1])
                    .unwrap()
                    .checked_sub(a[1].checked_mul(b[0]).unwrap())
                    .unwrap(),
            ]
        };
        for _ in 0..2000 {
            let off = |st: &mut u64| {
                [
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                ]
            };
            let v = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
            ];
            let (u, w) = (off(&mut st), off(&mut st));
            let n = cross(u, w); // normal to the V-plane
            // ε off-plane push: 0 (exactly degenerate), or a wide window down to ~1e-15
            // where the interval filter is ambiguous but the true sign is definite.
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            let plane_pts = |st: &mut u64| [v, add3(v, off(st)), add3(v, off(st))];
            let a = plane_pts(&mut st);
            let b = plane_pts(&mut st);
            let c = plane_pts(&mut st);
            // Triangle q=v+u, r=v+w, s=v+u+w+ε·n (coplanar with V=v, s pushed ε off).
            let qb = add3(v, u);
            let rb = add3(v, w);
            let sb = add3(add3(add3(v, u), w), smul(eps, n));
            // One shared rotation (shared pivot) for all twelve points.
            let axis = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let rp = |p: [Rat; 3]| Pt3::at(p).rotate_about(axis, ang, piv);
            let (a0, a1, a2) = (rp(a[0]), rp(a[1]), rp(a[2]));
            let (b0, b1, b2) = (rp(b[0]), rp(b[1]), rp(b[2]));
            let (c0, c1, c2) = (rp(c[0]), rp(c[1]), rp(c[2]));
            let (q, r, s) = (rp(qb), rp(rb), rp(sb));
            check(
                (&a0, &a1, &a2),
                (&b0, &b1, &b2),
                (&c0, &c1, &c2),
                &q,
                &r,
                &s,
            );
        }

        // Corpus C — a **constructed** four-plane concurrency, because a random one never happens
        // and this is the configuration that matters: every plane is spanned by one shared point
        // `p` and two others, so `p` is both the meet of any three and a defining point of the
        // fourth. With `p` in the third slot, `indirect_hp`'s `row1 = Dvec − D·s` cancels to
        // nothing — the shape a tool edge lying in a target plane produces, and the one that let a
        // wrong sign through a suite of random corpora.
        for _ in 0..400 {
            let axis = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 1, 359_000), rng(&mut st, 1, 997)); // inexact
            let pivot = rand_base(&mut st);
            let turn = |b: [Rat; 3]| Pt3::at(b).rotate_about(axis, ang, pivot);
            let p = turn(rand_base(&mut st));
            let sp: Vec<[Pt3; 2]> = (0..4)
                .map(|_| [turn(rand_base(&mut st)), turn(rand_base(&mut st))])
                .collect();
            for qi in 0..4 {
                let t: Vec<usize> = (0..4).filter(|&i| i != qi).collect();
                check(
                    (&p, &sp[t[0]][0], &sp[t[0]][1]),
                    (&p, &sp[t[1]][0], &sp[t[1]][1]),
                    (&p, &sp[t[2]][0], &sp[t[2]][1]),
                    &sp[qi][0],
                    &sp[qi][1],
                    &p,
                );
            }
        }

        eprintln!(
            "[H-c] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "indirect judge must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path (else a stuck filter passes vacuously)"
        );
    }

    /// Fast port check (default suite): a handful of generic configs — the indirect
    /// judge must match a moderate-precision truth (catches a transcription bug; full
    /// statistical soundness is the `#[ignore]`d H-b/H-c + exact3d).
    #[test]
    fn indirect_port_check() {
        const GT: usize = 384;
        let mut st = 0x1DEA_2C11_9F00_5A5Au64;
        for _ in 0..16 {
            let p: Vec<Pt3> = (0..12).map(|_| rand_point(&mut st)).collect();
            let a = (&p[0], &p[1], &p[2]);
            let b = (&p[3], &p[4], &p[5]);
            let c = (&p[6], &p[7], &p[8]);
            let judged =
                indirect_orient3d_judge(a, b, c, &p[9], &p[10], &p[11], fixture()).orient();
            if let Some(truth) = indirect_truth(a, b, c, &p[9], &p[10], &p[11], GT) {
                assert!(
                    judged == truth || judged == Orient::Zero,
                    "indirect judge {judged:?} disagrees with truth {truth:?}"
                );
            }
        }
    }

    // ---- indirect cmp_coord (cmp-i) ----

    fn add3(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
        [
            a[0].checked_add(b[0]).unwrap(),
            a[1].checked_add(b[1]).unwrap(),
            a[2].checked_add(b[2]).unwrap(),
        ]
    }

    /// Borrow an owned plane-triple as the `&Pt3` tuples the judge takes.
    fn tr(t: &[[Pt3; 3]; 3]) -> [(&Pt3, &Pt3, &Pt3); 3] {
        [
            (&t[0][0], &t[0][1], &t[0][2]),
            (&t[1][0], &t[1][1], &t[1][2]),
            (&t[2][0], &t[2][1], &t[2][2]),
        ]
    }

    /// Three planes meeting at `v` (each through `v` + two small offsets), all rotated by
    /// `(ax, ang, piv)` — an implicit point at `rotate(v)` with heterogeneous provenance.
    fn triple_pts(v: [Rat; 3], st: &mut u64, ax: Axis, ang: Angle, piv: [Rat; 3]) -> [[Pt3; 3]; 3] {
        let plane = |st: &mut u64| {
            let off = |st: &mut u64| {
                [
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                    ri(rng(st, -20, 20), rng(st, 1, 5)),
                ]
            };
            let (o1, o2) = (off(st), off(st));
            [
                Pt3::at(v).rotate_about(ax, ang, piv),
                Pt3::at(add3(v, o1)).rotate_about(ax, ang, piv),
                Pt3::at(add3(v, o2)).rotate_about(ax, ang, piv),
            ]
        };
        [plane(st), plane(st), plane(st)]
    }

    /// A random rotated three-plane triple (own random center, axis, inexact angle,
    /// pivot) — heterogeneous provenance. Sequences the RNG draws so each `&mut st`
    /// borrow ends before the next.
    fn rand_triple(st: &mut u64) -> [[Pt3; 3]; 3] {
        let v = rand_base(st);
        let ax = axis_of(rng(st, 0, 2));
        let angle = deg(rng(st, 0, 360_000), rng(st, 1, 9973));
        let piv = rand_base(st);
        triple_pts(v, st, ax, angle, piv)
    }

    /// The three axis-perpendicular planes through integer point `p` (meet exactly at `p`,
    /// tol 0) — an exact axis-aligned implicit point for the sanity oracle.
    fn axis_planes(p: [i128; 3]) -> [[Pt3; 3]; 3] {
        let pt = |x, y, z| Pt3::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
        let [x, y, z] = p;
        [
            [pt(x, y, z), pt(x, y + 1, z), pt(x, y, z + 1)], // ⊥ x
            [pt(x, y, z), pt(x + 1, y, z), pt(x, y, z + 1)], // ⊥ y
            [pt(x, y, z), pt(x + 1, y, z), pt(x, y + 1, z)], // ⊥ z
        ]
    }

    /// The high-precision cmp truth (`None` when the coordinates are equal or below the
    /// GT floor — a genuine tie).
    fn cmp_truth(
        a: [(&Pt3, &Pt3, &Pt3); 3],
        b: [(&Pt3, &Pt3, &Pt3); 3],
        axis: usize,
        prec: usize,
    ) -> Option<Orient> {
        let hp = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
            [
                plane_hp(t[0].0, t[0].1, t[0].2, prec),
                plane_hp(t[1].0, t[1].1, t[1].2, prec),
                plane_hp(t[2].0, t[2].1, t[2].2, prec),
            ]
        };
        cmp_hp_with_gap(hp(a), hp(b), axis, prec).ok()
    }

    /// `cmp_combine` maps the parity of negative signs to the ordering.
    #[test]
    fn cmp_combine_counts_negatives() {
        assert_eq!(
            cmp_combine(Some(true), Some(true), Some(true)),
            Some(Orient::Positive)
        );
        assert_eq!(
            cmp_combine(Some(false), Some(true), Some(true)),
            Some(Orient::Negative)
        );
        assert_eq!(
            cmp_combine(Some(false), Some(false), Some(true)),
            Some(Orient::Positive)
        );
        assert_eq!(cmp_combine(None, Some(true), Some(true)), None);
    }

    /// Sanity: two exact axis-aligned implicit points order by the compared axis; a swap
    /// flips it; an equal coordinate is `Zero`.
    #[test]
    fn cmp_sanity_axis_aligned() {
        let a = axis_planes([1, 2, 3]);
        let b = axis_planes([1, 5, 3]);
        // y: 2 < 5 → a below b → Negative; swap → Positive.
        assert_eq!(
            indirect_cmp_coord_judge(tr(&a), tr(&b), 1, fixture()).orient(),
            Orient::Negative
        );
        assert_eq!(
            indirect_cmp_coord_judge(tr(&b), tr(&a), 1, fixture()).orient(),
            Orient::Positive
        );
        // x and z equal → Zero.
        assert_eq!(
            indirect_cmp_coord_judge(tr(&a), tr(&b), 0, fixture()).orient(),
            Orient::Zero
        );
        assert_eq!(
            indirect_cmp_coord_judge(tr(&a), tr(&b), 2, fixture()).orient(),
            Orient::Zero
        );
    }

    /// Fast port check: the cmp judge matches a moderate-precision truth on generic
    /// rotated configs (catches a transcription bug; full soundness is `#[ignore]`d H-g).
    #[test]
    fn cmp_port_check() {
        const GT: usize = 384;
        let mut st = 0xC301_7A5E_2266_9911u64;
        for _ in 0..16 {
            let a = rand_triple(&mut st);
            let b = rand_triple(&mut st);
            let axis = rng(&mut st, 0, 2) as usize;
            let judged = indirect_cmp_coord_judge(tr(&a), tr(&b), axis, fixture()).orient();
            if let Some(truth) = cmp_truth(tr(&a), tr(&b), axis, GT) {
                assert!(
                    judged == truth || judged == Orient::Zero,
                    "cmp judge {judged:?} disagrees with truth {truth:?}"
                );
            }
        }
    }

    /// H-g — indirect cmp_coord soundness over heterogeneous provenance (corpus A) and a
    /// near-tie corpus (corpus B: two points sharing an axis coordinate up to ε, rotated
    /// about that same axis so the near-tie survives). The judge must never disagree with
    /// a GT-stable 512-bit truth; both the fast filter and the escalation are exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_g_indirect_cmp_coord_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut check = |a: &[[Pt3; 3]; 3], b: &[[Pt3; 3]; 3], axis: usize| {
            let (ta, tb) = (tr(a), tr(b));
            let judged = indirect_cmp_coord_judge(ta, tb, axis, fixture()).orient();
            let iv = |t: [(&Pt3, &Pt3, &Pt3); 3]| {
                [
                    plane_iv(t[0].0, t[0].1, t[0].2),
                    plane_iv(t[1].0, t[1].1, t[1].2),
                    plane_iv(t[2].0, t[2].1, t[2].2),
                ]
            };
            if cmp_filter(iv(ta), iv(tb), axis).is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match cmp_truth(ta, tb, axis, GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        };

        // Corpus A — heterogeneous provenance (each triple its own rotation and pivot).
        let mut st = 0x6A11_C0DE_5151_2323u64;
        for _ in 0..2000 {
            let a = rand_triple(&mut st);
            let b = rand_triple(&mut st);
            check(&a, &b, rng(&mut st, 0, 2) as usize);
        }

        // Corpus B — near-tie in z, shared Z-rotation (a Z-rotation leaves z unchanged, so
        // the pre-rotation z near-equality survives → M near 0 → escalation forced).
        for _ in 0..2000 {
            let z0 = rng(&mut st, -200, 200);
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            let va = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(z0, 1),
            ];
            let vb = [
                ri(rng(&mut st, -200, 200), 1),
                ri(rng(&mut st, -200, 200), 1),
                ri(z0, 1).checked_add(eps).unwrap(),
            ];
            let angle = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let a = triple_pts(va, &mut st, Axis::Z, angle, piv);
            let b = triple_pts(vb, &mut st, Axis::Z, angle, piv);
            check(&a, &b, 2);
        }

        eprintln!(
            "[H-g] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "cmp judge must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }

    // ---- dir_sign_judge (3a-iii) ----

    fn rat_cross(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
        let m = |x: Rat, y: Rat| x.checked_mul(y).unwrap();
        let s = |x: Rat, y: Rat| x.checked_sub(y).unwrap();
        [
            s(m(a[1], b[2]), m(a[2], b[1])),
            s(m(a[2], b[0]), m(a[0], b[2])),
            s(m(a[0], b[1]), m(a[1], b[0])),
        ]
    }

    fn smul(k: Rat, a: [Rat; 3]) -> [Rat; 3] {
        [
            a[0].checked_mul(k).unwrap(),
            a[1].checked_mul(k).unwrap(),
            a[2].checked_mul(k).unwrap(),
        ]
    }

    /// Three base points defining a plane whose normal is parallel to `n` — `p0` and two
    /// in-plane edges `n × e1`, `n × e2` (single crosses, small magnitude).
    fn plane_norm(n: [Rat; 3], p0: [Rat; 3]) -> [[Rat; 3]; 3] {
        let e1 = [ri(1, 1), ri(2, 1), ri(3, 1)];
        let e2 = [ri(2, 1), ri(3, 1), ri(1, 1)];
        [p0, add3(p0, rat_cross(n, e1)), add3(p0, rat_cross(n, e2))]
    }

    fn rot_plane(b: [[Rat; 3]; 3], ax: Axis, ang: Angle, piv: [Rat; 3]) -> [Pt3; 3] {
        b.map(|p| Pt3::at(p).rotate_about(ax, ang, piv))
    }

    fn t3(p: &[Pt3; 3]) -> (&Pt3, &Pt3, &Pt3) {
        (&p[0], &p[1], &p[2])
    }

    /// The three points are far from collinear (`sin²` of the corner angle above a
    /// threshold) — a well-formed plane. A degenerate plane (tiny normal) is not the
    /// near-coplanar-*normals* regime under test, so the corpus skips it.
    fn well_conditioned(p: &[Pt3; 3]) -> bool {
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (e1, e2) = (sub(p[1].coord, p[0].coord), sub(p[2].coord, p[0].coord));
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        dot(n, n) > 1e-6 * dot(e1, e1) * dot(e2, e2)
    }

    /// The **raw** sign of `D` (det of the three normals) at `prec` bits — no radius, no limit.
    fn dir_d_sign_at(
        a: (&Pt3, &Pt3, &Pt3),
        b: (&Pt3, &Pt3, &Pt3),
        c: (&Pt3, &Pt3, &Pt3),
        prec: usize,
    ) -> Option<bool> {
        let ph = |t: (&Pt3, &Pt3, &Pt3)| plane_hp(t.0, t.1, t.2, prec);
        // The **raw** sign, with no radius and no floor — see `dir_orient_at` for why the oracle
        // must not borrow the judge's bound. `dir_sign_truth` gets its confidence from two
        // precisions agreeing instead.
        let (dh, _) = cramer_hp(&[ph(a), ph(b), ph(c)], prec);
        (!dh.mid.is_zero()).then(|| dh.mid.is_positive())
    }

    /// The **GT-stable** `dir_sign` truth: `Some` only when `prec` and `prec + 128` agree
    /// on a definite sign — otherwise the config is degenerate beyond what the ground
    /// truth itself resolves, so it is `None` (skip), not a spurious "wrong". (Mirrors
    /// exact3d's `indirect_truth` stability check.)
    fn dir_sign_truth(
        a: (&Pt3, &Pt3, &Pt3),
        b: (&Pt3, &Pt3, &Pt3),
        c: (&Pt3, &Pt3, &Pt3),
        prec: usize,
    ) -> Option<Orient> {
        match (
            dir_d_sign_at(a, b, c, prec),
            dir_d_sign_at(a, b, c, prec + 128),
        ) {
            (Some(x), Some(y)) if x == y => Some(orient_of(x)),
            _ => None,
        }
    }

    /// Sanity: the three axis-perpendicular planes have normals `+x, −y, +z`, so their
    /// determinant is `−1`; a swap flips it, and three coplanar normals give `Zero`.
    #[test]
    fn dir_sign_judge_sanity() {
        let ap = axis_planes([0, 0, 0]);
        assert_eq!(
            dir_sign_judge(t3(&ap[0]), t3(&ap[1]), t3(&ap[2]), fixture()).orient(),
            Orient::Negative,
            "det[+x, -y, +z] = -1"
        );
        assert_eq!(
            dir_sign_judge(t3(&ap[0]), t3(&ap[2]), t3(&ap[1]), fixture()).orient(),
            Orient::Positive,
            "one swap flips the sign"
        );
        // Three normals in the plane z = 0 → coplanar → D = 0 → Zero.
        let mk = |n: [i128; 3]| {
            plane_norm([ri(n[0], 1), ri(n[1], 1), ri(n[2], 1)], [ri(0, 1); 3]).map(Pt3::at)
        };
        let (a, b, c) = (mk([1, 0, 0]), mk([0, 1, 0]), mk([1, 1, 0]));
        assert_eq!(
            dir_sign_judge(t3(&a), t3(&b), t3(&c), fixture()).orient(),
            Orient::Zero,
            "coplanar normals → D = 0"
        );
    }

    /// `dir_sign` is a determinant of normals, invariant under a shared rotation
    /// (`det(R·n) = det(R)·det(n) = det(n)`).
    #[test]
    fn dir_sign_rotation_invariant() {
        let mut st = 0x0D12_5157_ABCD_0007u64;
        for _ in 0..40 {
            let mk = |st: &mut u64| plane_norm(rand_base(st), rand_base(st));
            let (ba, bb, bc) = (mk(&mut st), mk(&mut st), mk(&mut st));
            let un = |b: [[Rat; 3]; 3]| b.map(Pt3::at);
            let (ua, ub, uc) = (un(ba), un(bb), un(bc));
            let s = dir_sign_judge(t3(&ua), t3(&ub), t3(&uc), fixture()).orient();
            if s == Orient::Zero {
                continue;
            }
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let (ra, rb, rc) = (
                rot_plane(ba, ax, ang, piv),
                rot_plane(bb, ax, ang, piv),
                rot_plane(bc, ax, ang, piv),
            );
            assert_eq!(
                dir_sign_judge(t3(&ra), t3(&rb), t3(&rc), fixture()).orient(),
                s,
                "rotation-invariant"
            );
        }
    }

    /// H-i — dir_sign soundness over a **near-coplanar-normals** corpus (which H-c/H-g do
    /// not stress: they force `M ≈ 0`, not `D ≈ 0`). Three plane normals `n0, n1,
    /// n2 = α·n0 + β·n1 + ε·(n0×n1)` (ε tiny → `D ≈ ε` → escalation), shared rotation. The
    /// judge must never disagree with a GT-stable 512-bit `D` sign; both paths exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_i_dir_sign_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut st = 0x0D18_5160_C0C0_2323u64;
        let rand_n = |st: &mut u64| {
            [
                ri(rng(st, -20, 20), 1),
                ri(rng(st, -20, 20), 1),
                ri(rng(st, -20, 20), 1),
            ]
        };
        for _ in 0..2000 {
            let n0 = rand_n(&mut st);
            let n1 = rand_n(&mut st);
            let (alpha, beta) = (ri(rng(&mut st, -5, 5), 1), ri(rng(&mut st, -5, 5), 1));
            let eps = match rng(&mut st, 0, 2) {
                0 => ri(0, 1),
                1 => ri(1, rng(&mut st, 5_000, 200_000)),
                _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
            };
            // n2 = α·n0 + β·n1 + ε·(n0×n1) — near-coplanar with n0, n1.
            let n2 = add3(
                add3(smul(alpha, n0), smul(beta, n1)),
                smul(eps, rat_cross(n0, n1)),
            );
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let mk =
                |n: [Rat; 3], st: &mut u64| rot_plane(plane_norm(n, rand_base(st)), ax, ang, piv);
            let pa = mk(n0, &mut st);
            let pb = mk(n1, &mut st);
            let pc = mk(n2, &mut st);
            if !(well_conditioned(&pa) && well_conditioned(&pb) && well_conditioned(&pc)) {
                continue; // a degenerate plane is not the regime under test
            }
            let judged = dir_sign_judge(t3(&pa), t3(&pb), t3(&pc), fixture()).orient();
            let (d, _) = cramer_iv([
                plane_iv(&pa[0], &pa[1], &pa[2]),
                plane_iv(&pb[0], &pb[1], &pb[2]),
                plane_iv(&pc[0], &pc[1], &pc[2]),
            ]);
            if d.sign().is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match dir_sign_truth(t3(&pa), t3(&pb), t3(&pc), GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        }
        eprintln!(
            "[H-i] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(wrong, 0, "dir_sign must never disagree with GT (soundness)");
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }

    /// The `dir_orient3d` determinant `det[d, x−base, y−base]` at `prec` bits — the GT /
    /// escalation realization (mirrors the judge's hp path).
    fn dir_orient_at(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3, prec: usize) -> Option<bool> {
        let dp = Pt3::at(d);
        let sub = |u: &HpIv, v: &HpIv| u.sub(v, prec);
        let (bh, xh, yh, dh) = (
            base.hp_coord(prec),
            x.hp_coord(prec),
            y.hp_coord(prec),
            dp.hp_coord(prec),
        );
        let rows = [
            dh,
            [
                sub(&xh[0], &bh[0]),
                sub(&xh[1], &bh[1]),
                sub(&xh[2], &bh[2]),
            ],
            [
                sub(&yh[0], &bh[0]),
                sub(&yh[1], &bh[1]),
                sub(&yh[2], &bh[2]),
            ],
        ];
        let det = det3_big_rows(&rows, prec);
        // The **raw** sign at `prec` bits, with no radius and no floor. Independence from the
        // judge is the whole point of an oracle: `dir_orient_truth` gets its confidence from two
        // precisions agreeing, not from any bound this file also ships to production.
        (!det.mid.is_zero()).then(|| det.mid.is_positive())
    }

    /// GT-stable truth: `Some` only when `prec` and `prec + 128` agree (else too degenerate
    /// for the ground truth itself — skip, not a spurious "wrong"). Mirrors `dir_sign_truth`.
    fn dir_orient_truth(d: [Rat; 3], base: &Pt3, x: &Pt3, y: &Pt3, prec: usize) -> Option<Orient> {
        match (
            dir_orient_at(d, base, x, y, prec),
            dir_orient_at(d, base, x, y, prec + 128),
        ) {
            (Some(a), Some(b)) if a == b => Some(orient_of(a)),
            _ => None,
        }
    }

    /// Sanity: `det[d, x−base, y−base] = d·((x−base)×(y−base))`. With `base=0, x=e_x, y=e_y`
    /// the edge normal is `+e_z`, so `d=+e_z → Positive`, `−e_z → Negative`, in-plane →
    /// `Zero`; swapping the edge pair flips the sign.
    #[test]
    fn dir_orient3d_judge_sanity() {
        let o = Pt3::at([ri(0, 1); 3]);
        let x = Pt3::at([ri(1, 1), ri(0, 1), ri(0, 1)]);
        let y = Pt3::at([ri(0, 1), ri(1, 1), ri(0, 1)]);
        let e = |a: i128, b: i128, c: i128| [ri(a, 1), ri(b, 1), ri(c, 1)];
        assert_eq!(
            dir_orient3d_judge(e(0, 0, 1), &o, &x, &y, FIXTURE_PREC),
            Orient::Positive
        );
        assert_eq!(
            dir_orient3d_judge(e(0, 0, -1), &o, &x, &y, FIXTURE_PREC),
            Orient::Negative
        );
        assert_eq!(
            dir_orient3d_judge(e(1, 0, 0), &o, &x, &y, FIXTURE_PREC),
            Orient::Zero,
            "d in the edge plane → det 0"
        );
        assert_eq!(
            dir_orient3d_judge(e(0, 0, 1), &o, &y, &x, FIXTURE_PREC),
            Orient::Negative,
            "swapping x,y flips the sign"
        );
    }

    /// `orient3d_ray(base, dir, x, y)` is exactly `orient3d(base, base+dir, x, y)`: on
    /// unrotated points (where `base+dir` *is* an exact `Pt3`) it must equal `orient3d_judge`
    /// with the ideal point materialized — the argument-for-argument reduction the ops
    /// ray-triangle caller relies on.
    #[test]
    fn orient3d_ray_matches_materialized() {
        let mut st = 0x0D1B_7A44_0F0F_9001u64;
        let at = |b: [Rat; 3]| Pt3::at(b);
        for _ in 0..200 {
            let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
            let dir = [
                ri(rng(&mut st, -9, 9), 1),
                ri(rng(&mut st, -9, 9), 1),
                ri(rng(&mut st, -9, 9), 1),
            ];
            let q = [
                bb[0].checked_add(dir[0]).unwrap(),
                bb[1].checked_add(dir[1]).unwrap(),
                bb[2].checked_add(dir[2]).unwrap(),
            ];
            let (base, x, y) = (at(bb), at(xb), at(yb));
            assert_eq!(
                orient3d_ray(&base, dir, &x, &y, FIXTURE_PREC),
                orient3d_judge(&base, &at(q), &x, &y, fixture()).orient(),
                "orient3d_ray == orient3d(base, base+dir, x, y)"
            );
        }
    }

    /// Soundness — `dir_orient3d` over a rotated corpus, oracle = a GT-stable 512-bit
    /// determinant (NOT rotation-invariance: `d` is fixed while the points rotate, so the
    /// sign is not preserved). Two regimes: random directions, and **near-grazing** (`d`
    /// almost in the edge plane, `det ≈ 0`) to force escalation. The judge must never claim
    /// a sign opposite to the GT; both the fast filter and astro-float paths are exercised.
    #[test]
    #[ignore = "slow astro-float ground truth (run with --ignored)"]
    fn h_dir_orient3d_soundness() {
        const GT: usize = 512;
        let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
            (0usize, 0, 0, 0, 0, 0);
        let mut st = 0x0D1A_5170_BEEF_4242u64;
        for _ in 0..1500 {
            let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
            let grazing = rng(&mut st, 0, 1) == 1;
            let ax = axis_of(rng(&mut st, 0, 2));
            let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
            let piv = rand_base(&mut st);
            let rp = |b: [Rat; 3]| Pt3::at(b).rotate_about(ax, ang, piv);
            let (base, x, y) = (rp(bb), rp(xb), rp(yb));
            let tri = [base.clone(), x.clone(), y.clone()];
            if !well_conditioned(&tri) {
                continue; // a degenerate edge fan is not the regime under test
            }
            // Two direction regimes. Random: any integer `d`. near-grazing: a large integer
            // multiple of the **rotated** in-plane edge `x−base` (built from the f64 coords),
            // so `d` is nearly ⊥ the rotated normal → `det ≈ 0` → the fast filter fails and
            // the astro-float path decides. Built after rotation, since the rotated normal is
            // what `d` must graze.
            let d = if !grazing {
                [
                    ri(rng(&mut st, -19, 19), 1),
                    ri(rng(&mut st, -19, 19), 1),
                    ri(rng(&mut st, -19, 19), 1),
                ]
            } else {
                let big = rng(&mut st, 100_000_000, 9_000_000_000) as f64;
                let comp = |k: usize| ri((big * (x.coord[k] - base.coord[k])).round() as i128, 1);
                [comp(0), comp(1), comp(2)]
            };
            let judged = dir_orient3d_judge(d, &base, &x, &y, FIXTURE_PREC);
            // Recompute the Iv filter to tally which path resolved (mirrors the judge).
            let dp = Pt3::at(d);
            let (bi, xi, yi, di) = (pt_iv(&base), pt_iv(&x), pt_iv(&y), pt_iv(&dp));
            let subi = |u: [Iv; 3], v: [Iv; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
            if det3_iv([di, subi(xi, bi), subi(yi, bi)]).sign().is_none() {
                escalated += 1;
            } else {
                filter_resolved += 1;
            }
            match dir_orient_truth(d, &base, &x, &y, GT) {
                None => skipped += 1,
                Some(truth) => {
                    tested += 1;
                    if judged == Orient::Zero {
                        declined += 1;
                    } else if judged != truth {
                        wrong += 1;
                    }
                }
            }
        }
        eprintln!(
            "[dir_orient3d] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
        );
        assert_eq!(
            wrong, 0,
            "dir_orient3d must never disagree with GT (soundness)"
        );
        assert!(escalated > 0, "corpus must exercise the escalation path");
        assert!(
            filter_resolved > 0,
            "corpus must exercise the fast filter path"
        );
    }
}
