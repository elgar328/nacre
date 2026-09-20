//! The toleranced point + its `orient3d` judgment.
//!
//! A rotated point cannot be held exactly (cos/sin are irrational), but its f64
//! realization carries a **direction-wise xyz tol** that soundly bounds
//! the error, accumulated as the definition is turned through a chain of
//! axis-aligned rotations (§CIP ①: `new tol = |R|·old + mix`). [`orient3d_judge`]
//! consumes that tol: an f64 determinant filter with a sound error bound
//! ([`det3_bound`], §CIP ②) decides the easy cases, the ambiguous ones **escalate**
//! to astro-float from the point definitions, and one whose interval still straddles zero is
//! **normalized into the distance it stands for** and held against the operation's coincidence
//! limit ([`Standard`]): below it the coincidence is proved, above it the judgement climbs to the
//! precision the shortfall names, and past the cap it is reported ([`Decision`]) instead of
//! assumed. This is the judgment **layer only**
//! ("층만") — wired into the boolean since the CIP stages that followed (`nacre-ops`'
//! `tolerant` module is the seam), and it is the *tol > 0*
//! path: a tol-0 (`Constructed`) config is faster/exact via `nacre-predicates`
//! (Shewchuk), routed by a higher layer, not here.
//!
//! Validated before the port by an isolated 3D experiment (verdict: GO): H-a (`det3_bound`
//! soundness), H-d/H-f (chain + arbitrary-pivot tol) — the bound never under-estimates the true
//! error (astro-float ground truth) over random heterogeneous-rotation configs.

use super::HP_RM;
use astro_float::BigFloat;
use nacre_scalar::{Angle, Axis, Bounded, HpBounded, Mag, Orient, Rat, rat_to_big};
#[cfg(feature = "parallel")]
use std::sync::{Arc as HpRc, OnceLock as HpOnce};
#[cfg(not(feature = "parallel"))]
use std::{cell::OnceCell as HpOnce, rc::Rc as HpRc};

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
    /// ([`nacre_scalar::inv_sqrt_f64`]), which is why this can be judged at all.
    ///
    /// ★★ **`v̂` is exact too, and that matters.** Realizing it as `ŵ × û` in f64 costs two
    /// roundings that do not cancel — a wall whose `v` is exactly `ẑ` came out three ulps short,
    /// which is a frame that is not quite orthonormal. `n ⊥ u_raw` by construction, so
    /// `|v_raw|² = |n|²·|u_raw|²` exactly, and one inverse square root of a rational lands it on
    /// the nose. [`nacre_scalar::plane_frame`] checks that product fits `i128` before handing the
    /// frame out, so the overflow is refused at the source rather than handled here.
    ///
    /// ★ **Proper** (`det = +1`), so it contributes nothing to a chain's mirror parity.
    Frame { frame: nacre_scalar::PlaneFrame },
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

/// The exact data of a wide sketch frame — [`nacre_scalar::PlaneFrame`]'s arbitrary-precision
/// twin. Built once when a motion chain is flattened; realized on demand.
///
/// `origin = origin_num / origin_den` (one exact rational, common denominator `n·n`); the
/// `*_raw` vectors and squared lengths mirror `PlaneFrame` field for field. All integers, all
/// exact — the constructor cannot fail, which is the point.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WideFrame {
    pub origin_num: [num_bigint::BigInt; 3],
    pub origin_den: num_bigint::BigInt,
    pub u_raw: [num_bigint::BigInt; 3],
    pub n: [num_bigint::BigInt; 3],
    pub v_raw: [num_bigint::BigInt; 3],
    pub uu: num_bigint::BigInt,
    pub nn: num_bigint::BigInt,
    pub vv: num_bigint::BigInt,
}

impl WideFrame {
    /// The **canonical placement's** frame for the plane `n·x + d = 0`: origin at the world
    /// origin's projection (`(−d/n·n)·n`, kept as `num/den`), axes by the arbitrary-axis
    /// convention (`u_raw = ẑ×n`, or `ŷ×n` when the normal is exactly vertical — the same
    /// exact branch [`nacre_scalar::plane_frame_default`] takes). The narrow twin of this
    /// derivation is `plane_frame_default` + `plane_frame_named`; here nothing can overflow,
    /// so unlike them this cannot decline on width.
    ///
    /// The realization is scale-free (each axis is divided by its own length; the origin uses
    /// this `n` with this `n·n`), so `n` need not be primitive.
    ///
    /// `None` only for a zero normal — not a plane.
    pub fn canonical(n: [num_bigint::BigInt; 3], d: &num_bigint::BigInt) -> Option<WideFrame> {
        use num_bigint::BigInt;
        let zero = BigInt::from(0);
        if n.iter().all(|c| *c == zero) {
            return None;
        }
        let u_raw = if n[0] == zero && n[1] == zero {
            [n[2].clone(), zero.clone(), zero.clone()] // ŷ × n — the vertical-normal branch
        } else {
            [-&n[1], n[0].clone(), zero.clone()] // ẑ × n
        };
        let cross = |a: &[BigInt; 3], b: &[BigInt; 3]| -> [BigInt; 3] {
            [
                &a[1] * &b[2] - &a[2] * &b[1],
                &a[2] * &b[0] - &a[0] * &b[2],
                &a[0] * &b[1] - &a[1] * &b[0],
            ]
        };
        let dot = |a: &[BigInt; 3]| -> BigInt { a.iter().map(|c| c * c).sum() };
        let v_raw = cross(&n, &u_raw);
        let (uu, nn) = (dot(&u_raw), dot(&n));
        // `n ⊥ u_raw` by construction, so `|v_raw|² = |n|²·|u_raw|²` exactly — the same
        // identity the narrow frame relies on.
        let vv = &nn * &uu;
        let neg_d = -d;
        let origin_num = [&neg_d * &n[0], &neg_d * &n[1], &neg_d * &n[2]];
        Some(WideFrame {
            origin_num,
            origin_den: nn.clone(),
            u_raw,
            n,
            v_raw,
            uu,
            nn,
            vv,
        })
    }

    /// A **named placement's** frame in arbitrary precision — the fallback for a caller-stated
    /// `(origin, ref_dir)` whose narrow `PlaneFrame` overflows `i128` (its squared lengths do;
    /// the caller's values themselves are `Rat` by construction).
    ///
    /// Unlike the arbitrary-axis `u_raw`, a general `ref_dir` may have a normal component, so
    /// it is projected out exactly: `u_raw = (n·n)·r − (r·n)·n` — after which `n ⊥ u_raw` and
    /// the `|v|² = |n|²·|u|²` identity holds as in [`WideFrame::canonical`].
    ///
    /// `None` for a zero normal or a `ref_dir` parallel to it — the same declines the narrow
    /// route makes, never a width one.
    pub fn named(
        n: [num_bigint::BigInt; 3],
        d: &num_bigint::BigInt,
        origin: &[Rat; 3],
        ref_dir: &[Rat; 3],
    ) -> Option<WideFrame> {
        use num_bigint::BigInt;
        use num_integer::Integer;
        let zero = BigInt::from(0);
        if n.iter().all(|c| *c == zero) {
            return None;
        }
        // Clear a rational vector's denominators — scale-free, so only the direction matters.
        let lift = |v: &[Rat; 3]| -> ([BigInt; 3], BigInt) {
            let den = [
                BigInt::from(v[0].denom()),
                BigInt::from(v[1].denom()),
                BigInt::from(v[2].denom()),
            ];
            let l = den.iter().fold(BigInt::from(1), |acc, x| acc.lcm(x));
            (
                core::array::from_fn(|k| BigInt::from(v[k].numer()) * (&l / &den[k])),
                l,
            )
        };
        let (r, _) = lift(ref_dir);
        let dot2 = |a: &[BigInt; 3], b: &[BigInt; 3]| -> BigInt {
            a.iter().zip(b).map(|(x, y)| x * y).sum()
        };
        let nn = dot2(&n, &n);
        let rn = dot2(&r, &n);
        let u_raw: [BigInt; 3] = core::array::from_fn(|k| &nn * &r[k] - &rn * &n[k]);
        if u_raw.iter().all(|c| *c == zero) {
            return None; // ref_dir parallel to the normal
        }
        let cross = |a: &[BigInt; 3], b: &[BigInt; 3]| -> [BigInt; 3] {
            [
                &a[1] * &b[2] - &a[2] * &b[1],
                &a[2] * &b[0] - &a[0] * &b[2],
                &a[0] * &b[1] - &a[1] * &b[0],
            ]
        };
        let v_raw = cross(&n, &u_raw);
        let uu = dot2(&u_raw, &u_raw);
        let vv = &nn * &uu;
        // The caller's origin, over one common denominator. `d` is unused beyond the plane's
        // identity — the origin is stated, not derived, and its on-plane invariant is checked
        // where the claim is made: `SketchFrame::named` rejects an off-plane origin by exact
        // residual (`plane_residual_sign`) before any frame is built, and `PlaneDef` holds it
        // structurally (`origin = points[0]`).
        let _ = d;
        let (origin_num, origin_den) = lift(origin);
        Some(WideFrame {
            origin_num,
            origin_den,
            u_raw,
            n,
            v_raw,
            uu,
            nn,
            vv,
        })
    }

    /// [`WideFrame::canonical`] from a plane's stored name, with the frame's `flip` spent here
    /// (negating the coefficients — the same place the narrow route spends it). `Narrow` names
    /// lift; `Wide` ones are already the right width.
    pub fn canonical_of(name: &nacre_scalar::PlaneName, flip: bool) -> Option<WideFrame> {
        let [c0, c1, c2, c3] = name_bigints(name, flip);
        WideFrame::canonical([c0, c1, c2], &c3)
    }

    /// [`WideFrame::named`] from a plane's stored name — see [`WideFrame::canonical_of`].
    pub fn named_of(
        name: &nacre_scalar::PlaneName,
        origin: &[Rat; 3],
        ref_dir: &[Rat; 3],
        flip: bool,
    ) -> Option<WideFrame> {
        let [c0, c1, c2, c3] = name_bigints(name, flip);
        WideFrame::named([c0, c1, c2], &c3, origin, ref_dir)
    }
}

/// A [`nacre_scalar::PlaneName`]'s coefficients as `BigInt`s, negated when `flip` — the sign a
/// frame node carries is spent on the coefficients, exactly as the narrow route spends it
/// before `plane_frame_named`.
fn name_bigints(name: &nacre_scalar::PlaneName, flip: bool) -> [num_bigint::BigInt; 4] {
    let mut cs = name.coeff_ints();
    if flip {
        for c in &mut cs {
            *c = -&*c;
        }
    }
    cs
}

/// The judged frame of a plane that has **no name** — the canonical placement derived from the
/// plane through three exact points whose motion chains need not agree.
///
/// [`MoveNode::Frame`] and [`MoveNode::FrameWide`] both demand exact coefficients (`Rat`,
/// `BigInt`); a mixed-frame datum plane has neither, because its exact world coefficients are
/// irrational. What it does have is three defining points that are each exact **in their own
/// frame** — so the coefficients exist as *intervals* at any precision ([`plane_hp`] realizes
/// each point independently and never required the chains to agree), and the canonical
/// placement (foot of perpendicular + arbitrary axis — the frozen spec) is derived from those
/// intervals with its error carried, never assumed away.
///
/// ★ `flip` is the same measured direction every frame node carries, spent the same way: the
/// four coefficients are negated before the basis is derived. There is no canonical-sign field
/// beside it — the canonical origin is sign-and-scale invariant (the doc's own words), the axes
/// only read direction, so one global negation is the entire freedom and `flip` is its name.
///
/// ★★ `vertical` is the arbitrary-axis branch (`ŷ×n` instead of `ẑ×n`), **decided once at the
/// fixed rung and stored**. Deciding per realization precision could change the basis between
/// the f64 cache and an escalation, and a frame that changes basis with precision is not a
/// frame; the fixed-rung realization is deterministic, so the same statement always decides the
/// same way. For a *named* plane the branch reads `n₀ = n₁ = 0` exactly; here exact zero is not
/// provable, so the branch is "provably usable" instead — `|ẑ×n|²` bounded away from zero — and
/// a plane that can prove neither branch is refused by the producer, by name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameThrough {
    pub points: [JudgedPoint; 3],
    pub vertical: bool,
    pub flip: bool,
}

/// One defining point of a judged frame — **how the point is stated**, mirroring the two mixed
/// causes.
///
/// Equality is definitional throughout ([`WitnessPoint`]'s own `PartialEq`), which is what the
/// statement-level determinism tests and `shared_base`'s whole-node comparison consume.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum JudgedPoint {
    /// Exact in its own frame — a rational base carried by a chain.
    Pure(WitnessPoint),
    /// ★ **The meet of its three carrier planes — a point with no coordinate anywhere**.
    /// Each carrier is described by its own witness triangle; the point is realized as the
    /// homogeneous `[Dvec : D]` of the carriers' Cramer system and never divided until a
    /// realization explicitly asks for an affine value. Boxed: nine points dwarf one.
    Meet(Box<[[WitnessPoint; 3]; 3]>),
}

/// A defining point as a **homogeneous** `(D, Dvec)` pair at `prec` bits — the uniform currency
/// the projective join takes. A pure point is the special case `[p : 1]` (`D` exactly one); a
/// meet realizes each carrier's interval coefficients ([`plane_hp`] — three points, chains free
/// to differ) and runs Cramer ([`cramer_hp`]). No division anywhere.
fn judged_homog(p: &JudgedPoint, prec: usize) -> (HpBounded, [HpBounded; 3]) {
    match p {
        JudgedPoint::Pure(wp) => (
            HpBounded::exact(BigFloat::from_f64(1.0, prec)),
            wp.hp_coord(prec),
        ),
        JudgedPoint::Meet(carriers) => {
            let plane =
                |t: &[WitnessPoint; 3]| -> [HpBounded; 4] { plane_hp(&t[0], &t[1], &t[2], prec) };
            let planes = [
                plane(&carriers[0]),
                plane(&carriers[1]),
                plane(&carriers[2]),
            ];
            cramer_hp(&planes, prec)
        }
    }
}

impl FrameThrough {
    /// The fixed rung the branch decision and the f64 cache realize at — the ladder's first
    /// rung, the same 128 bits [`WitnessPoint::frame_wide`] narrows from.
    const RUNG: usize = 128;

    /// Build the node: decide the branch, and prove the whole basis stands at the fixed rung.
    ///
    /// `None` when no branch's in-plane axis can be bounded away from zero, or when the basis
    /// derivation fails there (a normal that may vanish — the points may be collinear — or an
    /// origin denominator that may reach zero, or a meet whose `D` sign the rung cannot decide).
    /// ★ The caller must reject **by name and must not claim degeneracy**: none of these are
    /// proofs that the plane is degenerate, only failures to prove it healthy, and
    /// `CollinearVertices` would be a lie here.
    pub fn of(points: [JudgedPoint; 3], flip: bool) -> Option<FrameThrough> {
        let p = Self::RUNG;
        let c = judged_coeffs(&points, false, p)?;
        let sq_sum = |a: &HpBounded, b: &HpBounded| a.mul(a, p).add(&b.mul(b, p), p);
        // |ẑ×n|² = n₀²+n₁², |ŷ×n|² = n₂²+n₀². `flip` negates every coefficient, which no
        // squared length can see, so the decision is flip-invariant by construction.
        let proven = |x: &HpBounded| x.sign() == Some(true);
        let vertical = if proven(&sq_sum(&c[0], &c[1])) {
            false
        } else if proven(&sq_sum(&c[2], &c[0])) {
            true
        } else {
            return None;
        };
        let f = FrameThrough {
            points,
            vertical,
            flip,
        };
        // The whole basis must stand at the deciding rung — normal length and origin included —
        // so that the realization methods' `None` arms are genuinely unreachable.
        judged_basis(&f, p)?;
        Some(f)
    }

    /// The rung realization of defining point 0, as the affine `[f64; 3]` a cache anchors at —
    /// the definition's own replay (`exact.rs`' rule), never a second route through a vertex
    /// cache. For a meet this is the one place the homogeneous point is divided, and the divisor
    /// is safe by construction: [`FrameThrough::of`] proved the rung derivation, which includes
    /// this point's `D` sign.
    pub fn anchor_coord(&self) -> Option<[f64; 3]> {
        match &self.points[0] {
            JudgedPoint::Pure(wp) => Some(wp.coord()),
            JudgedPoint::Meet(_) => {
                let p = Self::RUNG;
                let (d, dvec) = judged_homog(&self.points[0], p);
                let mut out = [0.0; 3];
                for (o, num) in out.iter_mut().zip(&dvec) {
                    (*o, _) = narrow_hp(&num.div(&d, p)?);
                }
                Some(out)
            }
        }
    }
}

/// The judged plane's four coefficients at `prec`, `flip` applied — the front half of
/// [`judged_basis`], split out because [`FrameThrough::of`] needs the coefficients *before* a
/// branch exists to build the rest of the basis with.
///
/// All-pure triples keep the affine route ([`plane_hp`] — the road its two-road differential
/// locks); any meet switches the whole triple to the **homogeneous route**: every point as
/// `(D, Dvec)` ([`judged_homog`]) and the plane through them by the projective join
/// ([`plane_hp_through`]). `None` when any meet's `D` sign is
/// undecided at this precision — the join's normalization has nothing to stand on then.
fn judged_coeffs(points: &[JudgedPoint; 3], flip: bool, prec: usize) -> Option<[HpBounded; 4]> {
    let mut c = if let [
        JudgedPoint::Pure(p0),
        JudgedPoint::Pure(p1),
        JudgedPoint::Pure(p2),
    ] = points
    {
        plane_hp(p0, p1, p2, prec)
    } else {
        let h = [
            judged_homog(&points[0], prec),
            judged_homog(&points[1], prec),
            judged_homog(&points[2], prec),
        ];
        plane_hp_through(
            [(&h[0].0, &h[0].1), (&h[1].0, &h[1].1), (&h[2].0, &h[2].1)],
            prec,
        )?
    };
    if flip {
        let zero = HpBounded::exact(BigFloat::from_f64(0.0, prec));
        for k in &mut c {
            *k = zero.sub(k, prec);
        }
    }
    Some(c)
}

/// The judged frame's basis at `prec` bits: `[origin, û, v̂, ŵ]`, every component carrying the
/// error its own derivation incurred. The one place the derivation is spelled, so the f64 cache
/// ([`WitnessPoint::frame_through`]) and the escalation ([`WitnessPoint::hp_coord`]) cannot
/// disagree about what the frame *is*.
fn judged_basis(f: &FrameThrough, prec: usize) -> Option<[[HpBounded; 3]; 4]> {
    let zero = || HpBounded::exact(BigFloat::from_f64(0.0, prec));
    let neg = |x: &HpBounded| zero().sub(x, prec);
    let c = judged_coeffs(&f.points, f.flip, prec)?;
    let [n0, n1, n2, d] = c;
    let u_raw = if f.vertical {
        [n2.clone(), zero(), neg(&n0)] // ŷ × n — the vertical-normal branch
    } else {
        [neg(&n1), n0.clone(), zero()] // ẑ × n
    };
    let n = [n0, n1, n2];
    let dot = |a: &[HpBounded; 3], b: &[HpBounded; 3]| {
        a[0].mul(&b[0], prec)
            .add(&a[1].mul(&b[1], prec), prec)
            .add(&a[2].mul(&b[2], prec), prec)
    };
    let uu = dot(&u_raw, &u_raw);
    let nn = dot(&n, &n);
    let iu = uu.inv_sqrt(prec)?;
    let iw = nn.inv_sqrt(prec)?;
    let scale = |v: &[HpBounded; 3], s: &HpBounded| [0, 1, 2].map(|k| v[k].mul(s, prec));
    let u = scale(&u_raw, &iu);
    let w = scale(&n, &iw);
    // v̂ = ŵ × û — right-handed by construction, so the node is proper (`det = +1`) and
    // contributes nothing to a chain's mirror parity, like the other two frame nodes.
    let v = [0, 1, 2].map(|k| {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        w[i].mul(&u[j], prec).sub(&w[j].mul(&u[i], prec), prec)
    });
    // Foot of the perpendicular from the world origin, `(−d·n) / (n·n)` — the frozen canonical
    // convention, through the interval division that exists for exactly this quotient.
    let nd = neg(&d);
    let origin = [
        nd.mul(&n[0], prec).div(&nn, prec)?,
        nd.mul(&n[1], prec).div(&nn, prec)?,
        nd.mul(&n[2], prec).div(&nn, prec)?,
    ];
    Some([origin, u, v, w])
}

/// A rational base point carried through a chain of axis rotations (§CIP ⑦ rotation
/// history). `base` + `chain` are the exact **definition** (never lost); `realized` is
/// the f64 realization (a cache) with, per axis, the bound on its error — a **direction-wise
/// xyz vector** (§CIP ⑤) held beside the value it bounds, so the two cannot drift apart.
/// [`hp_coord`](Self::hp_coord) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,
    // (equality is definitional — see the `PartialEq` impl below the struct)
    pub realized: [Bounded; 3],
    /// Memoized `hp_coord` at the boolean's chosen precision — the astro-float realization
    /// once per definition and shared across clones (`Rc`). A judge escalates the *same*
    /// definition-point dozens of times per boolean (`plane_def` clones `tri_pt3` per call);
    /// without this each escalation replays the rotation's cos/sin at 200 bits, which dominates
    /// the rotated-boolean cost. `base`/`chain` never change after construction except through
    /// [`Self::rotate_about`], which resets this cell, so the cached value always matches the
    /// definition (a pure, path-independent function). See [`HpCell`] for the
    /// `Arc<OnceLock>` (parallel) vs `Rc<OnceCell>` (serial) choice.
    hp: HpCell,
}

/// **Definitional equality — `base` and `chain`, nothing else.** `realized`/`hp` are caches,
/// and realization is a pure function of the definition (path independence is this file's root
/// soundness argument), so two points with equal definitions cannot honestly disagree in their
/// caches. The consumer this exists for is [`shared_base`]'s whole-node chain comparison — it
/// declares two motions to be *one* motion, a question about definitions, never about caches —
/// which [`FrameThrough`] joins by carrying points inside a node.
impl PartialEq for WitnessPoint {
    fn eq(&self, o: &Self) -> bool {
        self.base == o.base && *self.chain == *o.chain
    }
}
impl Eq for WitnessPoint {}

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero) — the `f64` reading of
/// the bound [`Mag::above`] carries, for the tol arithmetic in this file, which is `f64`. Kept as
/// its own reader rather than `2^(Mag::exp2)`: that reads one octave high for a power of two.
fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// **How far a rational's f64 image sits from the rational** — measured at high precision, `0`
/// when the value is exactly representable.
///
/// One spelling, because three call sites want the identical quantity ([`WitnessPoint::at`]'s seed,
/// [`WitnessPoint::translate`]'s offset, [`WitnessPoint::frame`]'s inputs) and "how a rational's rounding is
/// charged" is exactly the kind of thing that drifts when it is written out twice. `bf_mag` reads
/// at octave granularity, so the `2×` keeps the result above the truth.
fn rat_round_tol(r: Rat, f: f64) -> f64 {
    let e = BigFloat::from_f64(f, 120).sub(&rat_to_big(r, 120), 120, HP_RM);
    2.0 * bf_mag(&e.abs())
}

impl WitnessPoint {
    /// The realized coordinate — the values of [`Self::realized`], for readers that want the
    /// point whole.
    #[inline]
    pub fn coord(&self) -> [f64; 3] {
        self.realized.map(|b| b.value)
    }

    /// The per-axis error bounds — the radii of [`Self::realized`], for readers that want them
    /// as a vector.
    #[inline]
    pub fn tol(&self) -> [f64; 3] {
        self.realized.map(|b| b.error)
    }

    /// A point at `base`, tol seeded with the base→f64 rounding (a division; exactly
    /// 0 for an f64-representable base, positive otherwise). An axis a later chain
    /// never rotates keeps exactly this, which a per-axis tol check needs.
    pub fn at(base: [Rat; 3]) -> Self {
        let coord = [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()];
        Self::at_with_tol(
            base,
            [
                rat_round_tol(base[0], coord[0]),
                rat_round_tol(base[1], coord[1]),
                rat_round_tol(base[2], coord[2]),
            ],
        )
    }

    /// A point at a rational `base`, its coordinate the nearest `f64` and its tol **from that
    /// contract rather than measured**: exactly `0` where the coordinate is representable, else
    /// `|x|·2⁻⁵³` (floored at `f64::MIN_POSITIVE`) — `Rat::to_f64` is documented as *"the
    /// **nearest** f64 … ties to even"*, so `|r − to_f64(r)| ≤ ½ ulp` holds for every normal
    /// `x` without anyone computing it.
    ///
    /// [`at`](Self::at) learns the same number by measuring at 120 bits — nine BigFloat
    /// operations per point. This is the spelling `nacre-ops` uses where the base is a
    /// **definition** (a plane's own points, a solid's vertices solved from their carriers): the
    /// bound is free, and its looseness can only turn a definite filter answer into "escalate",
    /// never the other way round.
    ///
    /// ★ **This is how a point that came from a rational reaches the predicates** — not by
    /// lifting the rounded `f64` back into a `Rat` (that names a different point, with a
    /// tolerance of zero on top).
    pub fn at_nearest(base: [Rat; 3]) -> Self {
        let f = base.map(|r| r.to_f64());
        let representable = base
            .iter()
            .zip(f)
            .all(|(&r, x)| Rat::try_from_f64(x) == Some(r));
        if representable {
            return Self::at_with_tol(base, [0.0; 3]);
        }
        let bound = |x: f64| (x.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE);
        Self::at_with_tol(base, [bound(f[0]), bound(f[1]), bound(f[2])])
    }

    /// A point at `base` with an explicit initial `tol` — for a root that already
    /// carries tol (a `Discovered` boolean seam), whose tol the chain then transports
    /// (`|R|·old`). [`at`](Self::at) is the Constructed case (initial tol = base
    /// rounding).
    pub fn at_with_tol(base: [Rat; 3], tol: [f64; 3]) -> Self {
        WitnessPoint {
            realized: [0, 1, 2].map(|k| Bounded::new(base[k].to_f64(), tol[k])),
            base,
            chain: HpRc::from([] as [MoveNode; 0]),
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
    ///
    /// **Two error terms, and neither is a guess.** `rot` carries the one input here without a
    /// rounding contract (`f64::cos`) and charges what `Angle::realization_error_of` *measured* of
    /// the very pair used below; `piv` is nothing but round-to-nearest steps and a realization, so it
    /// is derived and measured the way [`translate`](Self::translate) and [`mirror`](Self::mirror)
    /// derive and measure theirs.
    ///
    /// ★★ **`rot` is per axis, and it has to be.** The two coordinates are not the same expression:
    /// `i` mixes `u·c − v·s` while `j` mixes `u·s + v·c`, so `u` pairs with `cos`'s error on one
    /// and with `sin`'s on the other. A single shared term is only sound when the two errors are
    /// charged the same amount — which a *constant* did, and a measurement does not. Measured on a
    /// 356.65° turn, `sin` was eight times further out than `cos`, and the axis whose `u` met `sin`
    /// needed 2.4× what the other did. The old shared charge hid that; the ground-truth test found
    /// it the moment the constant left.
    pub fn rotate_about(self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let mut p = self.rotated(axis, angle, point);
        p.remember(MoveNode::Rotate {
            axis,
            angle,
            pivot: point,
        });
        p
    }

    /// [`Self::rotate_about`]'s numbers, without remembering the node — the step
    /// [`Self::apply_chain`] folds.
    fn rotated(mut self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let (i, j) = axis.plane();
        let (px, py) = (point[i].to_f64(), point[j].to_f64());
        let (ci, cj) = (self.realized[i].value, self.realized[j].value); // pre-rotation magnitudes for tol
        let (u, v) = (ci - px, cj - py);
        // cos/sin — exact (rational) for the 90°-family, else f64 (with realization tol).
        //
        // **Through `Angle`'s single entry point, not re-spelled here.** This function's contract
        // is to redo the producer's f64 route operation for operation (`Isometry::apply_point`,
        // which calls the same thing), and a second spelling of "how an angle becomes f64" is
        // exactly how the two would drift — silently, at the 90°-family, where one route snaps to
        // `0.0`/`±1.0` and the other lands `cos(90°) ≈ 6e-17`.
        let (c, s) = angle.cos_sin_f64();
        self.realized[i].value = px + u * c - v * s;
        self.realized[j].value = py + u * s + v * c;
        // **The rotation's own error — measured for this angle, not charged from a constant.**
        //
        // `dc`/`ds` are how far this platform's `cos`/`sin` land from the truth. That used to be a
        // constant measured once and written into a doc: sound only on machines like the one it was
        // taken on, which a kernel that ships to browsers cannot assume. Measured, a worse platform
        // reports a bigger number and the tolerance grows to match.
        //
        // ★ **Per axis, because the two coordinates mix the pair differently** — `u` meets `c` on
        // one and `s` on the other. See this function's doc for what that cost when it was shared.
        //
        // **The arithmetic of the products and their combination lives here too, and only here.**
        // For an origin pivot `piv` below is skipped entirely, so nothing else would cover
        // `fl(u·c)`, `fl(v·s)` and their combination. Each is `≤ ε/2` of a magnitude bounded by
        // `|u| + |v|`, so `1·ε` of that sum covers all three.
        //
        // ★ **Both are gated on the realization being inexact, and that gate is the invariant.**
        // When `cos`/`sin` are `0`/`±1` the products and their combination are *exact*, not merely
        // small — so a quadrantal origin rotation contributes nothing at all, which is what keeps
        // it at tol 0 (`quadrantal_origin_chain_is_tol_zero`) and keeps axis-aligned models on the
        // exact predicate path. An unconditional arithmetic term would break that quietly.
        //
        // ★★★ **It measures `c` and `s` — the values three lines above — not the angle.** They are
        // handed in rather than re-derived because a second realization of one angle is measured to
        // be able to differ from the first (`Angle`'s `F64_ERR` says where and why). An error taken
        // against a pair that never reached `coord` would bound a number this kernel never stored.
        let (dc, ds) = angle.realization_error_of(c, s);
        let (rot_i, rot_j) = if dc == 0.0 && ds == 0.0 {
            (0.0, 0.0)
        } else {
            let arith = f64::EPSILON * (u.abs() + v.abs());
            (
                u.abs() * dc + v.abs() * ds + arith,
                u.abs() * ds + v.abs() * dc + arith,
            )
        };
        // **The pivot arithmetic, charged the way `translate` and `mirror` charge theirs.**
        //
        // It used to take the trig constant — which existed *because `f64::cos` has no accuracy
        // contract* — for three operations that all do have one. Its two siblings in this `impl`
        // already do the right thing, and this now matches them: **measure** the rational's
        // realization, **count** the round-to-nearest steps.
        //
        // - The pivot's `Rat → f64` is measured, not counted, for the same reason they measure it:
        //   `Rat::to_f64` is one rounding for a numerator and denominator under `2⁵³` and takes
        //   another path above, so counting would mean knowing which. It enters **twice per axis
        //   with opposite signs** (`px + (ci − px)·c − …`) and partly cancels; `|1 − c| ≤ 2` and
        //   `|s| ≤ 1` bound the pair plainly, and the outer `2.0` is the margin every sibling
        //   term carries. **Exactly zero for a dyadic pivot** — the common case, which the old
        //   lumped charge still billed.
        // - **The roundings the pivot itself adds are three**: the two differences `ci − px` and
        //   `cj − py`, and the final sum `px + …`. (The products and their combination belong to
        //   `rot` above — they happen whether or not there is a pivot.) Each is `≤ ε/2` of a
        //   magnitude bounded by `|ci| + |cj| + |px| + |py|`, so `1.5·ε` of that sum covers them;
        //   `3.0` is that with the same doubling.
        //
        // Origin pivot stays exactly 0: `ci − 0.0` and `0.0 + x` are exact, so there is nothing
        // to charge — and the whole term is skipped rather than measured. This term *is* symmetric
        // in the two axes, unlike `rot`: it is built from magnitudes, not from which of `cos`/`sin`
        // each coordinate met.
        let piv = if px != 0.0 || py != 0.0 {
            let realized = |r: Rat, f: f64| {
                bf_mag(&rat_to_big(r, 120).sub(&BigFloat::from_f64(f, 120), 120, HP_RM)).abs()
            };
            2.0 * (2.0 * realized(point[i], px) + realized(point[j], py))
                + 3.0 * f64::EPSILON * (ci.abs() + cj.abs() + px.abs() + py.abs())
        } else {
            0.0
        };
        let (ti, tj) = (self.realized[i].error, self.realized[j].error);
        self.realized[i].error = c.abs() * ti + s.abs() * tj + rot_i + piv;
        self.realized[j].error = s.abs() * ti + c.abs() * tj + rot_j + piv;
        self
    }

    /// **Append one node to the remembered definition.**
    ///
    /// The slice is shared (`HpRc`) so a *clone* is a refcount bump rather than an allocation —
    /// which the judgment path does hundreds of thousands of times. The price is here: appending
    /// copies the slice, so applying `n` nodes one at a time costs `O(n²)`.
    ///
    /// ⚠★★★★ **That price was measured against the wrong caller, and the comment said so**: *"this
    /// runs when a solid is transformed, never on the judgment path, so the copy is not hot."* Cell
    /// 52 broke the premise — a vertex is now realized as it is pushed, which **replays a chain per
    /// push** — and the copy became the dominant cost of a long history: a 4,200-turn fixture went
    /// from 5.6 s to 563 s, and a 1,600-turn loop from 17 ms to 114 s. A caller holding a whole
    /// chain must therefore use [`Self::apply_chain`], which folds the numbers and appends **once**.
    fn remember(&mut self, node: MoveNode) {
        let mut nodes = self.chain.to_vec();
        nodes.push(node);
        self.chain = HpRc::from(nodes);
        // The definition changed — invalidate the memo of the old one. A fresh (unshared) cell, so
        // clones made before this keep their own cached value.
        self.hp = HpCell::default();
    }

    /// **The whole chain applied in one pass** — the numbers folded node by node, the definition
    /// remembered once. Equivalent to applying each node through its own method, and measured so:
    /// same coordinate, same tol, same high-precision realization, bit for bit.
    ///
    /// This is what a *replay* wants. Rebuilding the remembered chain per node (what the public
    /// per-node methods must do, since each is a transform in its own right) makes a replay
    /// quadratic in the history's length — see [`Self::remember`].
    pub fn apply_chain(mut self, chain: &[MoveNode]) -> Option<Self> {
        for n in chain {
            self = match n {
                MoveNode::Rotate { axis, angle, pivot } => self.rotated(*axis, *angle, *pivot),
                MoveNode::Translate { offset } => self.translated(*offset),
                MoveNode::Mirror { axis, offset } => self.mirrored(*axis, *offset),
                MoveNode::Frame { frame } => self.framed(*frame)?,
                MoveNode::FrameWide(f) => self.framed_wide(f)?,
                MoveNode::FrameThrough(f) => self.framed_through(f)?,
            };
        }
        let mut nodes = self.chain.to_vec();
        nodes.extend(chain.iter().cloned());
        self.chain = HpRc::from(nodes);
        self.hp = HpCell::default();
        Some(self)
    }

    /// The point reflected in `axis = offset` — one more link in the definition.
    ///
    /// **`coord` follows the producer's own route** (`AxisMirror::point`), so a replay of the
    /// definition reproduces the stored coordinate bit for bit. The exact offset lives in the
    /// chain, where [`compute_hp`](Self::compute_hp) realizes it.
    pub fn mirror(self, axis: Axis, offset: Rat) -> Self {
        let mut p = self.mirrored(axis, offset);
        p.remember(MoveNode::Mirror { axis, offset });
        p
    }

    /// [`Self::mirror`]'s numbers, without remembering the node.
    fn mirrored(mut self, axis: Axis, offset: Rat) -> Self {
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
        self.realized[k].error += 4.0
            * bf_mag(&rat_to_big(offset, 120).sub(&BigFloat::from_f64(c, 120), 120, HP_RM)).abs()
            + f64::EPSILON * (2.0 * c.abs() + self.realized[k].value.abs());
        // The producer's own route (`AxisMirror::point`), operation for operation — a replay of
        // the definition has to reproduce the stored coordinate bit for bit.
        self.realized[k].value = 2.0 * c - self.realized[k].value;
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
    pub fn translate(self, offset: [Rat; 3]) -> Self {
        let mut p = self.translated(offset);
        p.remember(MoveNode::Translate { offset });
        p
    }

    /// [`Self::translate`]'s numbers, without remembering the node.
    fn translated(mut self, offset: [Rat; 3]) -> Self {
        for (k, &off) in offset.iter().enumerate() {
            let t = off.to_f64();
            // The offset's realization error, plus the add's own half-ulp on the result.
            self.realized[k].error += 2.0
                * bf_mag(&rat_to_big(off, 120).sub(&BigFloat::from_f64(t, 120), 120, HP_RM)).abs()
                + f64::EPSILON * (self.realized[k].value.abs() + t.abs());
            self.realized[k].value += t;
        }
        self
    }

    /// The point read as `(u, v, w)` in the frame of the plane `origin`/`u_raw`/`n` describe, and
    /// written out in that plane's own frame — one more link in the definition.
    ///
    /// **`coord` follows the producer's route operation for operation**, so a replay of the
    /// definition reproduces the stored coordinate bit for bit. The exact inputs live in the
    /// chain, where [`compute_hp`](Self::compute_hp) realizes them.
    ///
    /// `None` when a squared length is not positive — that is a degenerate frame, not an unusual
    /// one, and a caller that cannot build a basis must not get a point that pretends otherwise.
    ///
    /// ★★★★ **Everything charged here is measured, and the exact cases really do reach zero.**
    /// The three sources are the rational inputs' own `Rat → f64` rounding, the two `1/√` values'
    /// realization error ([`nacre_scalar::inv_sqrt_error_of`] — *not* a half-ulp taken on faith,
    /// see there), and the f64 arithmetic that combines them. An axis-aligned frame has a signed
    /// permutation for a basis and every one of those terms vanishes: `±1` and `0` are exact, so
    /// the products are exact and the sums pick out one coordinate each.
    pub fn frame(self, f: nacre_scalar::PlaneFrame) -> Option<Self> {
        let mut p = self.framed(f)?;
        p.remember(MoveNode::Frame { frame: f });
        Some(p)
    }

    /// [`Self::frame`]'s numbers, without remembering the node.
    fn framed(mut self, f: nacre_scalar::PlaneFrame) -> Option<Self> {
        // One inverse square root per axis, each of an exact rational — see `plane_frame` for why
        // `v̂` gets its own instead of being a cross product of the other two.
        let inv = |v: Rat| -> Option<(f64, f64)> {
            let x = nacre_scalar::inv_sqrt_f64(v)?;
            Some((x, nacre_scalar::inv_sqrt_error_of(v, x)?))
        };
        let ((iu, du), (iw, dw)) = (inv(f.uu)?, inv(f.nn)?);
        // A basis vector's component, and a bound on how far it lands from the true one: the
        // rational's own rounding scaled by the length, the length's error scaled by the
        // rational, and the product's half-ulp. All three vanish for a `0`/`±1` component.
        let axis_comp = |r: Rat, s: f64, ds: f64| {
            let rf = r.to_f64();
            let c = rf * s;
            let err = rat_round_tol(r, rf) * s.abs() + rf.abs() * ds;
            (
                c,
                if err == 0.0 {
                    0.0
                } else {
                    err + f64::EPSILON * c.abs()
                },
            )
        };
        let (mut uh, mut vh, mut wh) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut eu, mut ev, mut ew) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        for k in 0..3 {
            (uh[k], eu[k]) = axis_comp(f.u_raw[k], iu, du);
            (wh[k], ew[k]) = axis_comp(f.n[k], iw, dw);
        }
        // ★★ `v̂` exactly where the frame carries `v_raw`, and `ŵ × û` where it does not.
        //
        // The exact route is one rounding per component; the cross product is two that do not
        // cancel — a wall whose `v` is exactly `ẑ` comes out three ulps short through it. Both are
        // sound; the fallback is what this crate had before `v_raw` existed, and it is taken only
        // when `|n|²·|u_raw|²` does not fit `i128` (see `nacre_scalar::plane_frame`).
        match f.v {
            Some((v_raw, vv)) => {
                let (iv, dv) = inv(vv)?;
                for k in 0..3 {
                    (vh[k], ev[k]) = axis_comp(v_raw[k], iv, dv);
                }
            }
            None => {
                // `|ŵ|, |û| ≤ 1`, so each of the four products carries the other factor's error at
                // unit scale; `2·(max eu + max ew)` covers all four with room. Gated on the basis
                // being inexact for the reason `rotate_about` gates its own: with `0`/`±1` factors
                // the products *and* their differences are exact, not merely small.
                let mx = |e: [f64; 3]| e.iter().fold(0.0f64, |a, &b| a.max(b));
                let (eu_max, ew_max) = (mx(eu), mx(ew));
                let e = if eu_max == 0.0 && ew_max == 0.0 {
                    0.0
                } else {
                    2.0 * (eu_max + ew_max) + 3.0 * f64::EPSILON
                };
                for k in 0..3 {
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    vh[k] = wh[i] * uh[j] - wh[j] * uh[i];
                    ev[k] = e;
                }
            }
        }
        // The combination itself: three products and three sums per coordinate, each a
        // round-to-nearest of a magnitude bounded by the terms' sum — except where the basis is a
        // signed permutation *and* the frame sits on the origin, where every product is exact and
        // every sum picks out a single coordinate. That case is the axis-aligned sketch, and it
        // is the one that has to stay at tol 0.
        let exact_basis = eu.iter().chain(&ev).chain(&ew).all(|&e| e == 0.0);
        let perm = exact_basis
            && uh
                .iter()
                .chain(&vh)
                .chain(&wh)
                .all(|c| *c == 0.0 || c.abs() == 1.0);
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let ok = f.origin[k].to_f64();
            let terms = p[0] * uh[k] + p[1] * vh[k] + p[2] * wh[k];
            // The incoming tol turned by the same basis, plus what this step adds.
            let carried = uh[k].abs() * t[0] + vh[k].abs() * t[1] + wh[k].abs() * t[2];
            let realized = p[0].abs() * eu[k] + p[1].abs() * ev[k] + p[2].abs() * ew[k];
            let arith = if perm && ok == 0.0 {
                0.0
            } else {
                rat_round_tol(f.origin[k], ok)
                    + 3.0
                        * f64::EPSILON
                        * (ok.abs()
                            + (p[0] * uh[k]).abs()
                            + (p[1] * vh[k]).abs()
                            + (p[2] * wh[k]).abs())
            };
            self.realized[k].value = ok + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`WideFrame`] — the same propagation, with the axis and origin
    /// `(value, error)` pairs taken from a fixed-precision arbitrary-precision realization
    /// instead of `inv_sqrt_f64`/`axis_comp`. No new f64 error derivation exists here: every
    /// rounding on the way is inside an `HpBounded`, and the final narrowing to f64 charges itself.
    ///
    /// ★ The axes are realized per applied point. The wide population is a fraction of a percent
    /// of pushes, so that cost is accepted rather than memoized — a measure-later item.
    ///
    /// ★ No exact-permutation shortcut: a wide frame is never an axis permutation (its squared
    /// lengths exceed `i128`), so the tol-0 branch `WitnessPoint::frame` has cannot apply.
    pub fn frame_wide(self, f: &WideFrame) -> Option<Self> {
        let mut p = self.framed_wide(f)?;
        p.remember(MoveNode::FrameWide(f.clone()));
        Some(p)
    }

    /// [`Self::frame_wide`]'s numbers, without remembering the node.
    fn framed_wide(mut self, f: &WideFrame) -> Option<Self> {
        // The fixed rung for the f64 cache of a wide frame — the ladder's first rung, the same
        // one `inv_sqrt_f64` starts at. The judgment path re-realizes at its own precision.
        const P: usize = 128;
        let inv = |v: &num_bigint::BigInt| nacre_scalar::inv_sqrt_bigint_bounded(v, P);
        let (iu, iv2, iw) = (inv(&f.uu)?, inv(&f.vv)?, inv(&f.nn)?);
        let comp = |raw: &num_bigint::BigInt, s: &HpBounded| {
            narrow_hp(&HpBounded::of_bigint(raw, P).mul(s, P))
        };
        let (mut uh, mut vh, mut wh) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut eu, mut ev, mut ew) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut oh, mut eo) = ([0.0; 3], [0.0; 3]);
        let od = HpBounded::of_bigint(&f.origin_den, P);
        for k in 0..3 {
            (uh[k], eu[k]) = comp(&f.u_raw[k], &iu);
            (vh[k], ev[k]) = comp(&f.v_raw[k], &iv2);
            (wh[k], ew[k]) = comp(&f.n[k], &iw);
            (oh[k], eo[k]) =
                narrow_hp(&HpBounded::of_bigint(&f.origin_num[k], P).div_exact(&od, P)?);
        }
        // The same propagation as `WitnessPoint::frame`, with the realized origin's own error in place
        // of `rat_round_tol`.
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let terms = p[0] * uh[k] + p[1] * vh[k] + p[2] * wh[k];
            let carried = uh[k].abs() * t[0] + vh[k].abs() * t[1] + wh[k].abs() * t[2];
            let realized = p[0].abs() * eu[k] + p[1].abs() * ev[k] + p[2].abs() * ew[k];
            let arith = eo[k]
                + 3.0
                    * f64::EPSILON
                    * (oh[k].abs()
                        + (p[0] * uh[k]).abs()
                        + (p[1] * vh[k]).abs()
                        + (p[2] * wh[k]).abs());
            self.realized[k].value = oh[k] + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`FrameThrough`] — the judged basis realized
    /// at the fixed rung, narrowed to `(value, error)` pairs, then **the same propagation as
    /// [`WitnessPoint::frame_wide`]**: the incoming tol is turned by the basis, the basis's own
    /// realization error is scaled by the coordinates, and the combination charges its rounding.
    ///
    /// `None` is unreachable for a node built by [`FrameThrough::of`] — it proved this exact
    /// derivation at this exact rung — and stays an `Option` so that "unreachable" is a fact
    /// about the producer rather than an invariant this method asserts across a crate boundary.
    ///
    /// ★ Cost note: the node's points share their `hp` cells by `Rc`, so the cos/sin of their
    /// chains realize once per node, not once per applied point — what recurs per point is
    /// arithmetic, the same acceptance [`WitnessPoint::frame_wide`] records.
    pub fn frame_through(self, f: &FrameThrough) -> Option<Self> {
        let mut p = self.framed_through(f)?;
        p.remember(MoveNode::FrameThrough(Box::new(f.clone())));
        Some(p)
    }

    /// [`Self::frame_through`]'s numbers, without remembering the node.
    fn framed_through(mut self, f: &FrameThrough) -> Option<Self> {
        let basis = judged_basis(f, FrameThrough::RUNG)?;
        let [o, u, v, w] = basis.map(|row| row.map(|c| narrow_hp(&c)));
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let (uh, eu) = u[k];
            let (vh, ev) = v[k];
            let (wh, ew) = w[k];
            let (oh, eo) = o[k];
            let terms = p[0] * uh + p[1] * vh + p[2] * wh;
            let carried = uh.abs() * t[0] + vh.abs() * t[1] + wh.abs() * t[2];
            let realized = p[0].abs() * eu + p[1].abs() * ev + p[2].abs() * ew;
            let arith = eo
                + 3.0
                    * f64::EPSILON
                    * (oh.abs() + (p[0] * uh).abs() + (p[1] * vh).abs() + (p[2] * wh).abs());
            self.realized[k].value = oh + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain, each node about its pivot) — path-independent ground truth /
    /// escalation realization. The result is memoized in
    /// [`WitnessPoint::hp`] and shared across clones of the same definition, so a definition-point pays
    /// the astro-float cos/sin once per boolean rather than once per predicate.
    pub(crate) fn hp_coord(&self, prec: usize) -> [HpBounded; 3] {
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

    /// **The coordinate realized at `prec` bits from the definition, with the error it carries** —
    /// the door [`Self::hp_coord`] is behind, in the type this workspace hands out publicly
    /// (`nacre_scalar::HpBounded`, what `inv_sqrt_bounded` returns too).
    ///
    /// ★ **The radius is half the answer.** A value without the bound its realization cost cannot
    /// be rounded honestly — see [`HpBounded`]'s contract: a radius invented for convenience makes
    /// every sign above it unearned. Callers round with `nacre_scalar::round_to_f64` /
    /// `round_to_digits`, which report *undecided* rather than picking when `prec` is short.
    ///
    /// ★★ Calling this at a precision the boolean is not using is safe: [`Self::hp_coord`]'s memo
    /// is keyed by precision and recomputes rather than disturbing the cached value.
    pub fn realize(&self, prec: usize) -> [HpBounded; 3] {
        self.hp_coord(prec)
    }

    /// The uncached realization (the body of [`Self::hp_coord`]), **with the error it carries**.
    ///
    /// This is the same walk as [`rotate_about`](Self::rotate_about)'s tol propagation, one level
    /// up: the definition is exact, the realization is not, and the radius is what the realization
    /// cost. Nothing here is a chosen constant — the base contributes its own rounding (zero when
    /// the rational lands on a `prec`-bit dyadic), each `cos`/`sin` contributes the bound
    /// [`Angle::cos_sin_bounded`] derives, and every arithmetic operation adds its half-ulp.
    fn compute_hp(&self, prec: usize) -> [HpBounded; 3] {
        let p = [
            HpBounded::of_rat(self.base[0], prec),
            HpBounded::of_rat(self.base[1], prec),
            HpBounded::of_rat(self.base[2], prec),
        ];
        fold_suffix(p, &self.chain, prec)
    }
}

/// **Carry a realized point through `nodes`, continuing from wherever it already is.**
///
/// [`WitnessPoint::compute_hp`] is this function seeded from the definition's base, and that is the
/// whole relationship: the fold is a left fold whose entire state is the three coordinates, so the
/// value after `k` nodes *is* the prefix, and resuming from it reaches the same bits as folding the
/// prefix again. That is what lets a caller who kept the previous generation's value pay for one
/// node instead of the whole history.
///
/// ⚠ **One fold, not two.** The per-node arithmetic lives in [`fold_one`] and nothing else spells
/// it, because two copies of a rounding rule are two rules that drift.
pub fn fold_suffix(prefix: [HpBounded; 3], nodes: &[MoveNode], prec: usize) -> [HpBounded; 3] {
    nodes.iter().fold(prefix, |p, n| fold_one(p, n, prec))
}

/// One node applied to a realized point — the body of the fold, spelled once.
///
/// ⚠ The degenerate arms **return `p` untouched**. They read as `continue` when this was a loop
/// body, and "skip this node" is what they meant; extracted, the same meaning is "this node moves
/// nothing". All three are unreachable (their producers refuse a non-positive squared length), and
/// they stay total rather than inventing a bound.
fn fold_one(mut p: [HpBounded; 3], node: &MoveNode, prec: usize) -> [HpBounded; 3] {
    match node {
        MoveNode::Rotate { axis, angle, pivot } => {
            let (i, j) = axis.plane();
            let (c, s) = angle.cos_sin_bounded(prec);
            let (px, py) = (
                HpBounded::of_rat(pivot[i], prec),
                HpBounded::of_rat(pivot[j], prec),
            );
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
                p[k] = p[k].add(&HpBounded::of_rat(offset[k], prec), prec);
            }
        }
        // `2·offset − x` on one coordinate: exact input, so the interval carries only its
        // own rounding. The other two coordinates are untouched.
        MoveNode::Mirror { axis, offset } => {
            let k = axis.index();
            let c = HpBounded::of_rat(*offset, prec);
            p[k] = c.add(&c, prec).sub(&p[k], prec);
        }
        // The same derivation as [`WitnessPoint::frame`], in the domain that carries its own
        // error: two `1/√` intervals, two scaled basis vectors, their cross product, and
        // the combination. Nothing is charged by hand — every rounding is inside an
        // `HpBounded`, which is the point of realizing here rather than trusting the f64 tol.
        //
        // A degenerate frame cannot arise here: `WitnessPoint::frame` is the only producer of this
        // node and it refuses a non-positive squared length, so the chain never holds one.
        MoveNode::Frame { frame } => {
            let inv = |v: Rat| nacre_scalar::inv_sqrt_bounded(v, prec);
            let (Some(iu), Some(iw)) = (inv(frame.uu), inv(frame.nn)) else {
                return p; // unreachable — `plane_frame` refuses a non-positive length
            };
            let scaled = |v: [Rat; 3], s: &HpBounded| {
                [0, 1, 2].map(|k| HpBounded::of_rat(v[k], prec).mul(s, prec))
            };
            let (uh, wh) = (scaled(frame.u_raw, &iu), scaled(frame.n, &iw));
            // The same two routes `WitnessPoint::frame` takes, in the domain that carries its own
            // error interval.
            let vh = match frame.v.and_then(|(v_raw, vv)| Some((v_raw, inv(vv)?))) {
                Some((v_raw, iv)) => scaled(v_raw, &iv),
                None => [0, 1, 2].map(|k| {
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    wh[i].mul(&uh[j], prec).sub(&wh[j].mul(&uh[i], prec), prec)
                }),
            };
            p = [0, 1, 2].map(|k| {
                HpBounded::of_rat(frame.origin[k], prec)
                    .add(&p[0].mul(&uh[k], prec), prec)
                    .add(&p[1].mul(&vh[k], prec), prec)
                    .add(&p[2].mul(&wh[k], prec), prec)
            });
        }
        // [`MoveNode::Frame`]'s wide twin — same shape, arbitrary-precision inputs.
        // `v_raw`/`vv` always exist (nothing overflows a `BigInt`), so there is no
        // cross-product fallback branch here.
        MoveNode::FrameWide(f) => {
            let inv = |v: &num_bigint::BigInt| nacre_scalar::inv_sqrt_bigint_bounded(v, prec);
            let (Some(iu), Some(iv2), Some(iw)) = (inv(&f.uu), inv(&f.vv), inv(&f.nn)) else {
                return p; // unreachable — the builder derives positive squared lengths
            };
            let scaled = |v: &[num_bigint::BigInt; 3], s: &HpBounded| {
                [0, 1, 2].map(|k| HpBounded::of_bigint(&v[k], prec).mul(s, prec))
            };
            let (uh, vh, wh) = (
                scaled(&f.u_raw, &iu),
                scaled(&f.v_raw, &iv2),
                scaled(&f.n, &iw),
            );
            let od = HpBounded::of_bigint(&f.origin_den, prec);
            let origin = (|| {
                let o =
                    |k: usize| HpBounded::of_bigint(&f.origin_num[k], prec).div_exact(&od, prec);
                Some([o(0)?, o(1)?, o(2)?])
            })();
            let Some(o) = origin else {
                return p; // unreachable — the builder's denominator is n·n > 0
            };
            p = [0, 1, 2].map(|k| {
                o[k].add(&p[0].mul(&uh[k], prec), prec)
                    .add(&p[1].mul(&vh[k], prec), prec)
                    .add(&p[2].mul(&wh[k], prec), prec)
            });
        }
        // The judged frame — the same derivation the f64 rung took, at this precision.
        MoveNode::FrameThrough(f) => {
            // Construction proved the derivation at the fixed rung, and every input
            // radius shrinks as `prec` grows (they are the points' own realization
            // errors), so the primary call is expected to succeed; the rung fallback
            // keeps the arm total without inventing a bound if that expectation is ever
            // wrong — a 128-bit basis is a *sound*, merely wider, interval for the same
            // true frame.
            let Some(basis) = judged_basis(f, prec).or_else(|| judged_basis(f, FrameThrough::RUNG))
            else {
                return p; // unreachable — `FrameThrough::of` proved the rung derivation
            };
            let [o, u, v, w] = basis;
            p = [0, 1, 2].map(|k| {
                o[k].add(&p[0].mul(&u[k], prec), prec)
                    .add(&p[1].mul(&v[k], prec), prec)
                    .add(&p[2].mul(&w[k], prec), prec)
            });
        }
    }
    p
}

/// A high-precision interval narrowed to an `(f64 value, f64 error)` pair — how a wide frame's
/// realization reaches the f64 cache. Correctly rounded when the interval decides it (then the
/// error is the radius plus the value's own half-ulp); an undecided or out-of-normal-range value
/// flushes to zero with its whole magnitude charged — sound, and reachable only for axis
/// components below f64's normal floor.
fn narrow_hp(x: &HpBounded) -> (f64, f64) {
    match nacre_scalar::round_to_f64(&x.value, x.error, 128) {
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
enum Gap {
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
fn denom_lo(v: &HpBounded) -> Result<Mag, Gap> {
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

fn escalate(
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
/// `rat_to_hp` returns an exact interval exactly when both hold *at this precision*. Below 127 a
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
fn distance_bound(det_rad: Mag, cross: &[HpBounded; 3], prec: usize) -> Gap {
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
/// the 2D `det_bound` the retired `frame2` module carried (the name is kept because the
/// derivation is the same one; nothing links to it — that module is gone). Validated in
/// exact3d (H-a).
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
/// owns its rows makes the caller copy each value up to three times over. An `HpBounded` copy is a
/// mantissa heap allocation, and this determinant is the hot path — [`cramer_hp`] alone made
/// forty-five of them per call for values it only ever read.
fn det3_big(r: [[&HpBounded; 3]; 3], prec: usize) -> HpBounded {
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let m0 = mul(r[1][1], r[2][2]).sub(&mul(r[1][2], r[2][1]), prec);
    let m1 = mul(r[1][0], r[2][2]).sub(&mul(r[1][2], r[2][0]), prec);
    let m2 = mul(r[1][0], r[2][1]).sub(&mul(r[1][1], r[2][0]), prec);
    mul(r[0][0], &m0)
        .sub(&mul(r[0][1], &m1), prec)
        .add(&mul(r[0][2], &m2), prec)
}

/// [`det3_big`] over rows the caller already owns — the borrow, spelled once.
fn det3_big_rows(r: &[[HpBounded; 3]; 3], prec: usize) -> HpBounded {
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
fn det3_hp(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    prec: usize,
) -> HpBounded {
    let (a, b, c, d) = (
        pa.hp_coord(prec),
        pb.hp_coord(prec),
        pc.hp_coord(prec),
        pd.hp_coord(prec),
    );
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
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
fn shared_base<const N: usize>(pts: &[&WitnessPoint; N]) -> Option<[[f64; 3]; N]> {
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
pub fn orient3d_filter(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
) -> Option<Orient> {
    let (a, b, c, d) = (pa.coord(), pb.coord(), pc.coord(), pd.coord());
    let det = det3_f64(a, b, c, d);
    let bound = det3_bound([a, b, c, d], [pa.tol(), pb.tol(), pc.tol(), pd.tol()]);
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
pub fn orient3d_judge(
    pa: &WitnessPoint,
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    j: Standard,
) -> Decision {
    if let Some(o) = orient3d_filter(pa, pb, pc, pd) {
        return Decision::Sign(o);
    }
    escalate(j, j.coincidence, |prec| {
        let det = det3_hp(pa, pb, pc, pd, prec);
        match det.sign() {
            Some(pos) => Ok(orient_of(pos)),
            // Only now is the area term worth forming: the sign decides on the first try in all
            // but the coincident cases, and the cross product is six high-precision products.
            None => Err(distance_bound(det.error, &cross_of(pb, pc, pd, prec), prec)),
        }
    })
}

/// `(b−d) × (c−d)` at `prec` bits — the area term that turns [`orient3d_judge`]'s determinant
/// into a height above the plane through `b, c, d`.
fn cross_of(
    pb: &WitnessPoint,
    pc: &WitnessPoint,
    pd: &WitnessPoint,
    prec: usize,
) -> [HpBounded; 3] {
    let (b, c, d) = (pb.hp_coord(prec), pc.hp_coord(prec), pd.hp_coord(prec));
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
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

/// CIP `dir_orient3d`: the sign of `det[d, x−base, y−base] = d·((x−base)×(y−base))` — the
/// orientation of the ray direction `d` against the edge fan `(base→x, base→y)`. The
/// **direction analogue** of [`orient3d_judge`]: `point_in_solid`'s ray-triangle test asks
/// `orient3d(p, p+d, ·, ·)`, but `p+d` (a rotated point plus a rational offset) has no exact
/// `base+chain` `WitnessPoint` (`R⁻¹d` is irrational). Every such determinant reduces to this form,
/// where `d` enters as one **exact** (error-0) column and only `x, y, base` carry rotation tol.
/// Interval filter → astro-float escalation, exactly like the indirect judges; a below-floor
/// determinant that stays undecided is [`Orient::Zero`], absorbed by the caller's ray retry.
///
/// **The one judge with no metric normalization**, and the reason is that its caller does not
/// need one: a zero here means the ray runs along the face's plane, and `boolean` answers that by
/// casting a different ray rather than by asking how close it came.
///
/// `d` is realized through an unrotated [`WitnessPoint`] purely to reuse `pt_iv`/`hp_coord` for its
/// coord/tol/hp — the row is `d` itself, never `d − base`.
pub fn dir_orient3d_judge(
    d: [Rat; 3],
    base: &WitnessPoint,
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Orient {
    let dp = WitnessPoint::at(d);
    let (bi, xi, yi, di) = (base.realized, x.realized, y.realized, dp.realized);
    let sub_iv =
        |u: [Bounded; 3], v: [Bounded; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
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
    let sub_hp = |u: &HpBounded, v: &HpBounded| u.sub(v, prec);
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
pub fn orient3d_ray(
    base: &WitnessPoint,
    dir: [Rat; 3],
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Orient {
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
// judgment **layer only** ("층만") — the boolean wiring (seam → plane `WitnessPoint`s) is
// stage 3. Validated before the port by the isolated 3D experiment (H-b coefficient
// tol, H-c indirect soundness); the constant `mag`-floor policy is indirect-only (distinct from the
// explicit `16·scale³` floor of [`orient3d_judge`]).

/// 3×3 determinant of interval rows.
fn det3_iv(r: [[Bounded; 3]; 3]) -> Bounded {
    let m0 = r[1][1].mul(r[2][2]).sub(r[1][2].mul(r[2][1]));
    let m1 = r[1][0].mul(r[2][2]).sub(r[1][2].mul(r[2][0]));
    let m2 = r[1][0].mul(r[2][1]).sub(r[1][1].mul(r[2][0]));
    r[0][0].mul(m0).sub(r[0][1].mul(m1)).add(r[0][2].mul(m2))
}

/// Plane `[a,b,c,d]` (`n·X + d = 0`) through three points, as intervals: `n =
/// (p1−p0)×(p2−p0)`, `d = −n·p0`. Coefficient tol propagates from the point tols
/// through the subtraction/cross/dot — "coefficient tol is a corollary of point tol"
/// (§CIP ②). Validated H-b.
pub(crate) fn plane_iv(p0: &WitnessPoint, p1: &WitnessPoint, p2: &WitnessPoint) -> [Bounded; 4] {
    let (a, b, c) = (p0.realized, p1.realized, p2.realized);
    let e1 = [b[0].sub(a[0]), b[1].sub(a[1]), b[2].sub(a[2])];
    let e2 = [c[0].sub(a[0]), c[1].sub(a[1]), c[2].sub(a[2])];
    let n = [
        e1[1].mul(e2[2]).sub(e1[2].mul(e2[1])),
        e1[2].mul(e2[0]).sub(e1[0].mul(e2[2])),
        e1[0].mul(e2[1]).sub(e1[1].mul(e2[0])),
    ];
    let d = Bounded::new(0.0, 0.0)
        .sub(n[0].mul(a[0]))
        .sub(n[1].mul(a[1]))
        .sub(n[2].mul(a[2]));
    [n[0], n[1], n[2], d]
}

/// The plane's four coefficients realized at `prec` bits from the definitions
/// (ground truth for the coefficient tol; no trig of its own — consumes point
/// `hp_coord`).
fn plane_hp(
    p0: &WitnessPoint,
    p1: &WitnessPoint,
    p2: &WitnessPoint,
    prec: usize,
) -> [HpBounded; 4] {
    let (a, b, c) = (p0.hp_coord(prec), p1.hp_coord(prec), p2.hp_coord(prec));
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let e1 = [sub(&b[0], &a[0]), sub(&b[1], &a[1]), sub(&b[2], &a[2])];
    let e2 = [sub(&c[0], &a[0]), sub(&c[1], &a[1]), sub(&c[2], &a[2])];
    let n = [
        mul(&e1[1], &e2[2]).sub(&mul(&e1[2], &e2[1]), prec),
        mul(&e1[2], &e2[0]).sub(&mul(&e1[0], &e2[2]), prec),
        mul(&e1[0], &e2[1]).sub(&mul(&e1[1], &e2[0]), prec),
    ];
    let d = HpBounded::exact(BigFloat::from_f64(0.0, prec))
        .sub(&mul(&n[0], &a[0]), prec)
        .sub(&mul(&n[1], &a[1]), prec)
        .sub(&mul(&n[2], &a[2]), prec);
    let [n0, n1, n2] = n; // `d` is already built from them, so the normal moves out rather than copying
    [n0, n1, n2, d]
}

/// `true` when an odd number of the three scale factors is negative — the case where the join
/// comes out pointing the other way. `sign()` reports `true` for positive.
fn negatives_odd(signs: [bool; 3]) -> bool {
    signs.iter().filter(|positive| !**positive).count() % 2 == 1
}

/// **The plane through three points that have no coordinates** — the projective join, the dual
/// of [`cramer_hp`]: that one meets three planes into a homogeneous point `[Dvec : D]`, this
/// joins three such points back into a plane. Both are the same 3×3 determinant expansion read
/// the other way round, which is why the whole thing is four [`det3_big`] calls and no new
/// arithmetic. Consumed by [`judged_coeffs`]' homogeneous route (a judged frame whose defining
/// point is a [`JudgedPoint::Meet`]).
///
/// ★★★ **Nothing is divided, and that is the point.** A point named as an intersection has no
/// rational coordinates to hand out — dividing `Dvec/D` would manufacture an irrational and bake
/// its rounding into everything downstream, which is exactly what the homogeneous form exists to
/// avoid. The plane is built in the same currency: multiplies and subtracts.
///
/// ★★★★ **The scale's sign is removed here, not handed to the caller.** The join of
/// `P_i = D_i·[p_i : 1]` is multilinear in its rows, so the answer is `D0·D1·D2` times the plane
/// through the affine points. A positive multiple is harmless — every question asked of a plane
/// is a sign question — but a **negative** one silently reverses the plane's direction, and
/// direction is what `frame_sign`, outward normals and the whole label frame are built on; the
/// function negates its own result instead of reporting (the kernel's existing idiom:
/// `canonical_plane_coeffs` also settles sign *inside* the value).
///
/// `None` when any `D_i`'s sign is undecided — then it is not even known whether those three
/// planes meet in a point, so there is nothing to normalize against and the caller must climb.
/// Degree 9 in the nine input coefficients.
///
/// (It has no f64/`Bounded` filter twin: the judging table needs no interval-coefficient route
/// at all. The lock lives in the Pure-vs-Meet differential, which pins this
/// function's basis against the affine shortcut end to end.)
fn plane_hp_through(
    pts: [(&HpBounded, &[HpBounded; 3]); 3],
    prec: usize,
) -> Option<[HpBounded; 4]> {
    let zero = HpBounded::exact(BigFloat::from_f64(0.0, prec));
    let row = |i: usize| {
        let (d, v) = pts[i];
        [&v[0], &v[1], &v[2], d]
    };
    let (r0, r1, r2) = (row(0), row(1), row(2));
    let minor = |a: usize, b: usize, c: usize| {
        det3_big(
            [
                [r0[a], r0[b], r0[c]],
                [r1[a], r1[b], r1[c]],
                [r2[a], r2[b], r2[c]],
            ],
            prec,
        )
    };
    let mut pi = [
        minor(1, 2, 3),
        zero.sub(&minor(0, 2, 3), prec),
        minor(0, 1, 3),
        zero.sub(&minor(0, 1, 2), prec),
    ];
    if negatives_odd([pts[0].0.sign()?, pts[1].0.sign()?, pts[2].0.sign()?]) {
        for c in &mut pi {
            *c = zero.sub(c, prec);
        }
    }
    Some(pi)
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
/// Dvec[j]/D` (no division taken here). Shared by `indirect_filter` (orient3d) and
/// [`cmp_filter`] (cmp_coord).
pub(crate) fn cramer_iv(planes: [[Bounded; 4]; 3]) -> (Bounded, [Bounded; 3]) {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    let h = |k: usize| Bounded::new(0.0, 0.0).sub(planes[k][3]); // n·X = h, h = −d
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
fn normals_det_iv(planes: [[Bounded; 4]; 3]) -> Bounded {
    let n = |k: usize| [planes[k][0], planes[k][1], planes[k][2]];
    det3_iv([n(0), n(1), n(2)])
}

/// [`normals_det_iv`] at `prec` bits — the `d` half of [`cramer_hp`] without its `h` column or its
/// three `Dvec` determinants. The escalation is where a wasted determinant costs the most.
fn normals_det_hp(planes: &[[HpBounded; 4]; 3], prec: usize) -> HpBounded {
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
fn indirect_filter(
    planes: [[Bounded; 4]; 3],
    q: [Bounded; 3],
    r: [Bounded; 3],
    s: [Bounded; 3],
) -> Option<Orient> {
    filter_from_cramer(cramer_iv(planes), q, r, s)
}

/// `indirect_filter` with the implicit point's Cramer parts **already in hand**.
///
/// ★ **The split exists because `(D, Dvec)` *is* the point** — it does not mention `q, r, s`, so a
/// caller asking about one point against several query triangles pays for it once. The arrangement's
/// crossing collector asks twice in a row (a segment's two endpoints), which is what
/// [`crate::predicate::ImplicitPoint`] exploits.
fn filter_from_cramer(
    (d, dvec): (Bounded, [Bounded; 3]),
    q: [Bounded; 3],
    r: [Bounded; 3],
    s: [Bounded; 3],
) -> Option<Orient> {
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
/// `indirect_hp` (orient3d) and `cmp_hp_with_gap` (cmp_coord).
///
/// The magnitude bounds these used to return alongside are gone: the radius rides *with* the
/// value now, so there is nothing left for a caller to forget to use — which is exactly how the
/// indirect judge came to bound `M` by a scale that had already cancelled.
fn cramer_hp(planes: &[[HpBounded; 4]; 3], prec: usize) -> (HpBounded, [HpBounded; 3]) {
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let zero = HpBounded::exact(BigFloat::from_f64(0.0, prec));
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
    planes: [[HpBounded; 4]; 3],
    q: [HpBounded; 3],
    r: [HpBounded; 3],
    s: [HpBounded; 3],
    prec: usize,
) -> (HpBounded, HpBounded, [HpBounded; 3]) {
    let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
    let mul = |x: &HpBounded, y: &HpBounded| x.mul(y, prec);
    let add = |x: &HpBounded, y: &HpBounded| x.add(y, prec);
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
fn plane_gap(m: &HpBounded, d: &HpBounded, cross: &[HpBounded; 3], prec: usize) -> Gap {
    let d_lo = match denom_lo(d) {
        Ok(b) => b,
        Err(g) => return g,
    };
    match distance_bound(m.error, cross, prec) {
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
    plane_a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    j: Standard,
) -> Decision {
    let iv = [
        plane_iv(plane_a.0, plane_a.1, plane_a.2),
        plane_iv(plane_b.0, plane_b.1, plane_b.2),
        plane_iv(plane_c.0, plane_c.1, plane_c.2),
    ];
    orient3d_from_cramer(cramer_iv(iv), plane_a, plane_b, plane_c, q, r, s, j)
}

/// [`indirect_orient3d_judge`] with the implicit point's Cramer parts **already in hand** — the one
/// body of the certified `orient3d`, reached by both [`crate::predicate::Judge::orient3d`] and
/// [`crate::predicate::ImplicitPoint::orient3d`].
///
/// The `WitnessPoint` definitions stay in the signature because the escalation still needs them: `plane_hp`
/// realizes the coefficients afresh at each precision, and an interval cannot be sharpened after
/// the fact.
#[allow(clippy::too_many_arguments)]
pub(crate) fn orient3d_from_cramer(
    cr: (Bounded, [Bounded; 3]),
    plane_a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    plane_c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    j: Standard,
) -> Decision {
    if let Some(o) = filter_from_cramer(cr, q.realized, r.realized, s.realized) {
        return Decision::Sign(o);
    }
    // Escalate: the same two determinants at prec, each carrying the radius accumulated
    // along its own computation. A determinant whose interval straddles zero is undecided,
    // and what it *did* establish is how far the implicit point may be from the triangle's plane.
    escalate(j, j.coincidence, |prec| {
        let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
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
fn cmp_filter(a: [[Bounded; 4]; 3], b: [[Bounded; 4]; 3], axis: usize) -> Option<Orient> {
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
    a: [[HpBounded; 4]; 3],
    b: [[HpBounded; 4]; 3],
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
fn coord_gap(m: &HpBounded, da: &HpBounded, db: &HpBounded) -> Gap {
    let (lo_a, lo_b) = match (denom_lo(da), denom_lo(db)) {
        (Ok(a), Ok(b)) => (a, b),
        // Both may be short; the deeper shortfall is the one that has to be covered.
        (Err(Gap::Unresolved { short: x }), Err(Gap::Unresolved { short: y })) => {
            return Gap::Unresolved { short: x.max(y) };
        }
        (Err(g), _) | (_, Err(g)) => return g,
    };
    match m.error.over(lo_a.times(lo_b)) {
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
    a: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    b: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    axis: usize,
    j: Standard,
) -> Decision {
    let iv = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
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
        let hp = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
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
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    j: Standard,
) -> Decision {
    // ★ **Not fed from `Judge`'s interval-plane cache, and that is measured.** Wiring it here
    // bought 1.007x on the crossing collector: this judgement's cost is almost entirely the
    // escalation below, because a *true* zero — the walls meeting `wc` in no point, 8.9% of the
    // pairs the collector tests — can never be settled by an interval filter.
    //
    // ★★★ **And this is where the kernel's escalations come from.** Measured on `fin_fold(80)`:
    // this judgement is 5.4% of its own calls but **78% of every escalation in the boolean**, and
    // **92.6% of those escalations have a midpoint of exactly `0.0`** — genuine degeneracies the
    // arrangement builds by construction. That is why tightening a coordinate's tol does not buy
    // speed: setting the rotation's whole tol contribution to zero left the escalation count here
    // **bit-identical**. The lever, if one is ever wanted, is recognising a true zero *before*
    // climbing — not more accurate coordinates.
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
        let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
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
fn dir_gap(d: &HpBounded, planes: &[[HpBounded; 4]; 3], prec: usize) -> Gap {
    let mut denom = Mag::of(1.0);
    for p in planes {
        // `|n| ≥ max|n_k|`, which is enough and needs no square root. One component clearing zero
        // is all a normal needs, so a shortfall only counts when **every** component fell short —
        // and then the smallest of them is the cheapest way out.
        let mut lo = Mag::ZERO;
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
    match d.error.over(denom) {
        Some(b) => Gap::Of(b),
        None => Gap::Vanished,
    }
}

#[cfg(test)]
#[path = "tests/frame3/mod.rs"]
mod tests;
