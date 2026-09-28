use super::*;
/// The exact data of a wide sketch frame — [`nacre_exact::PlaneFrame`]'s arbitrary-precision
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
    /// exact branch [`nacre_exact::plane_frame_default`] takes). The narrow twin of this
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
    pub fn canonical_of(name: &nacre_exact::PlaneName, flip: bool) -> Option<WideFrame> {
        let [c0, c1, c2, c3] = name_bigints(name, flip);
        WideFrame::canonical([c0, c1, c2], &c3)
    }

    /// [`WideFrame::named`] from a plane's stored name — see [`WideFrame::canonical_of`].
    pub fn named_of(
        name: &nacre_exact::PlaneName,
        origin: &[Rat; 3],
        ref_dir: &[Rat; 3],
        flip: bool,
    ) -> Option<WideFrame> {
        let [c0, c1, c2, c3] = name_bigints(name, flip);
        WideFrame::named([c0, c1, c2], &c3, origin, ref_dir)
    }
}

/// A [`nacre_exact::PlaneName`]'s coefficients as `BigInt`s, negated when `flip` — the sign a
/// frame node carries is spent on the coefficients, exactly as the narrow route spends it
/// before `plane_frame_named`.
fn name_bigints(name: &nacre_exact::PlaneName, flip: bool) -> [num_bigint::BigInt; 4] {
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
/// `BigInt`); a mixed-frame datum plane has neither, because no one frame solves its three
/// points rationally. What it does have is three defining points that are each exact **in their own
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
    pub(super) const RUNG: usize = 128;

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
    /// the definition's own replay (`construct.rs`' rule), never a second route through a vertex
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
pub(super) fn judged_basis(f: &FrameThrough, prec: usize) -> Option<[[HpBounded; 3]; 4]> {
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
