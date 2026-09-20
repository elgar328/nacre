use super::*;
/// **The world origin projected onto a rational plane** — `p = (−d / n·n) · n` for
/// `a·x + b·y + c·z + d = 0`.
///
/// This is the point of the plane closest to `(0, 0, 0)`, and it is what a face's sketch frame
/// takes for its origin. The alternative — the face's area centroid — is computed in `f64` from the
/// face's own vertices, and lifting *that* back into a rational makes a rounded cache into the
/// truth, which is the one thing construction here must never do.
///
/// ★★★ **The answer does not depend on how the plane is spelled.** A plane has a family of
/// coefficient vectors and this formula is invariant across all of them:
///
/// ```text
/// sign:  ( −(−d) / n·n )·(−n)   =  ( −d / n·n )·n
/// scale: ( −λd / λ²(n·n) )·λn   =  ( −d / n·n )·n
/// ```
///
/// So it is a function of the *plane*, not of the vector describing it — which matters because
/// [`canonical_plane_coeffs`] deliberately carries no direction (`[0,0,1,−3]` and `[0,0,−1,3]`
/// canonicalize together), and because two faces of one plane may hold it either way round.
/// Canonicalizing first is therefore not needed for correctness, only for overflow headroom.
///
/// `None` when `n·n = 0` — the coefficients are not a plane — or on `i128` overflow, which is the
/// kernel's ordinary demotion signal: the caller keeps its f64 path.
pub fn plane_origin_projection(coeffs: [Rat; 4]) -> Option<[Rat; 3]> {
    let zero = Rat::from_int(0);
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let mut nn = zero;
    for c in n {
        nn = nn.checked_add(c.checked_mul(c)?)?;
    }
    if nn == zero {
        return None;
    }
    let neg_d = zero.checked_sub(coeffs[3])?;
    // `Rat` exposes no division — the sign and scale invariance above is what makes one here
    // legitimate, and `Ratio`'s checked division reports the overflow the rest of the crate does.
    let mut out = [zero; 3];
    for i in 0..3 {
        out[i] = Rat(neg_d.checked_mul(n[i])?.0.checked_div(&nn.0)?);
    }
    Some(out)
}

/// **The point where three rational planes meet**, exactly — the `Rat` Cramer route first, and
/// where an intermediate overflows `i128`, the same answer through [`three_planes_big`]'s
/// integer core, denominators cleared once per row (a plane row is scale-free, and scaling a
/// row of a linear system leaves its solution). Input rows are `[a, b, c, d]` for
/// `a·x + b·y + c·z + d = 0`.
///
/// `None` means exactly two things, **neither of them the arithmetic**: the determinant is zero
/// (no unique point — parallel or line-sharing planes), or the point itself does not fit `Rat`
/// ([`MeetPoint::Wide`]). It used to also mean "an intermediate overflowed" — the rational
/// cofactor expansion builds `a.num·b.den ± b.num·a.den` before it can reduce — and that
/// conflation silently cost a decimal-framed tool its whole class reuse: its constructed
/// corners solve to points that fit `Rat` (measured 8/8), but the road there overflowed.
/// The fallback runs only on the decline path, so the narrow
/// route's cost and answers are untouched.
///
/// This is what replaces a dissolved sketch-frame base vertex: the frame-shared triple of a
/// prism corner, solved in the frame the planes are stated in, is the corner's exact base —
/// measured bit-identical to the stored base-and-replay road (8/8).
pub fn three_planes_rat(p: [[Rat; 4]; 3]) -> Option<[Rat; 3]> {
    three_planes_rat_narrow(p).or_else(|| {
        use num_bigint::BigInt;
        use num_integer::Integer;
        let lift = |row: &[Rat; 4]| -> [BigInt; 4] {
            let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(row[i].denom()));
            let l = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
            core::array::from_fn(|i| BigInt::from(row[i].numer()) * (&l / &den[i]))
        };
        three_planes_int([lift(&p[0]), lift(&p[1]), lift(&p[2])])?
            .narrow()
            .copied()
    })
}

/// [`three_planes_rat`]'s narrow route — Cramer over `Rat` with checked arithmetic, whose
/// `None` still conflates "no unique point" with "an intermediate overflowed". That is fine
/// *here*: the caller above resolves the conflation by asking the integer core.
pub(super) fn three_planes_rat_narrow(p: [[Rat; 4]; 3]) -> Option<[Rat; 3]> {
    let zero = Rat::from_int(0);
    // 3×3 determinant by cofactor expansion, all checked.
    let det3 = |m: [[Rat; 3]; 3]| -> Option<Rat> {
        let minor = |r0: usize, r1: usize, c0: usize, c1: usize| -> Option<Rat> {
            m[r0][c0]
                .checked_mul(m[r1][c1])?
                .checked_sub(m[r0][c1].checked_mul(m[r1][c0])?)
        };
        m[0][0]
            .checked_mul(minor(1, 2, 1, 2)?)?
            .checked_sub(m[0][1].checked_mul(minor(1, 2, 0, 2)?)?)?
            .checked_add(m[0][2].checked_mul(minor(1, 2, 0, 1)?)?)
    };
    let rhs = [
        zero.checked_sub(p[0][3])?,
        zero.checked_sub(p[1][3])?,
        zero.checked_sub(p[2][3])?,
    ];
    let d = det3([
        [p[0][0], p[0][1], p[0][2]],
        [p[1][0], p[1][1], p[1][2]],
        [p[2][0], p[2][1], p[2][2]],
    ])?;
    if d == zero {
        return None;
    }
    let with_col = |c: usize| -> [[Rat; 3]; 3] {
        let mut m = [
            [p[0][0], p[0][1], p[0][2]],
            [p[1][0], p[1][1], p[1][2]],
            [p[2][0], p[2][1], p[2][2]],
        ];
        for r in 0..3 {
            m[r][c] = rhs[r];
        }
        m
    };
    let mut out = [zero; 3];
    for (i, o) in out.iter_mut().enumerate() {
        *o = Rat(det3(with_col(i))?.0.checked_div(&d.0)?);
    }
    Some(out)
}

/// Where three planes meet, at whatever width the answer needs — [`three_planes_big`]'s result.
///
/// The fork mirrors [`PlaneName`]'s: `Narrow` **whenever the answer fits**, so two routes to one
/// point are structurally equal.
///
/// ★ The reduction is **per coordinate**, because that is what `Rat` — a `Ratio<i128>` each — has
/// to hold. Cramer hands back `detᵢ/det` sharing one denominator, and that shared form is
/// systematically wider than the coordinates it names; measuring it would report the width of the
/// *arithmetic* rather than of the point.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MeetPoint {
    Narrow([Rat; 3]),
    /// `(numerator, denominator)` per coordinate, each in lowest terms with a positive
    /// denominator — at least one of which does not fit `i128`.
    Wide([(num_bigint::BigInt, num_bigint::BigInt); 3]),
}

impl MeetPoint {
    /// The `Rat` form when there is one — `None` (`Wide`) means the point cannot be a `Rat`
    /// triple at all, which is a different fact from [`three_planes_rat`] declining.
    #[inline]
    pub fn narrow(&self) -> Option<&[Rat; 3]> {
        match self {
            MeetPoint::Narrow(p) => Some(p),
            MeetPoint::Wide(_) => None,
        }
    }

    /// **(integer components, positive common denominator)** — `lift3`'s twin for a point whose
    /// width may exceed `Rat`.
    ///
    /// This is what lets a sign predicate take a meet *whatever* its width: the answer is a
    /// polynomial in these integers, so `Narrow` and `Wide` travel the same road and no
    /// consumer has to carry a width limit inside its refusal.
    pub fn lift(&self) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
        use num_bigint::BigInt;
        use num_integer::Integer;
        match self {
            MeetPoint::Narrow(p) => lift3(p),
            MeetPoint::Wide(p) => {
                let d = p.iter().fold(BigInt::from(1), |l, (_, den)| l.lcm(den));
                let num = core::array::from_fn(|i| {
                    let (n, den) = &p[i];
                    n * (&d / den)
                });
                (num, d)
            }
        }
    }

    /// **The width `Rat` would have to hold** — the widest of the six magnitudes (numerator and
    /// denominator of each coordinate) *after* the per-coordinate reduction. `≤ 127` is exactly
    /// the `Narrow` condition, so this is the quantity, not a proxy for it.
    ///
    /// ★ Reduced, deliberately: the unreduced `detᵢ/det` Cramer produces is systematically wider
    /// and would report how the answer was computed rather than how wide the answer is.
    pub fn width_bits(&self) -> u64 {
        match self {
            MeetPoint::Narrow(p) => p
                .iter()
                .flat_map(|r| [r.numer(), r.denom()])
                .map(|v| (128 - v.unsigned_abs().leading_zeros()) as u64)
                .max()
                .unwrap_or(0),
            MeetPoint::Wide(p) => p
                .iter()
                .flat_map(|(n, d)| [n.bits(), d.bits()])
                .max()
                .unwrap_or(0),
        }
    }
}

/// **The same meeting point [`three_planes_rat`] computes, without the `i128` ceiling** — the
/// wide twin that [`plane_name_exact`] has had on the plane-*name* side all along and the vertex
/// solve did not.
///
/// [`three_planes_rat`] is `checked_*` throughout, so its `None` conflates two different facts:
/// the point does not fit `Rat`, and an **intermediate** of the rational cofactor expansion
/// overflowed while the point itself would have fit. Only a route without the ceiling can tell
/// those apart — and telling them apart is what decides whether a limit belongs to the *type* or
/// to the *arithmetic*.
///
/// **Integers, not rationals** — [`plane_name_big`]'s argument, applied per row instead of per
/// point: clearing each row's denominators once up front leaves plain integer arithmetic, with
/// the reductions at the end where the content has to come out anyway. ★ **A plane row is
/// scale-free**, so scaling row `k` by its own denominators' lcm leaves the same plane, and
/// scaling a row of a linear system leaves the same solution.
///
/// Takes [`PlaneName`]s rather than `[Rat; 4]` rows so a `Wide` carrier — one whose canonical
/// name no longer fits `Rat` — is solvable too. That population is invisible to
/// [`three_planes_rat`], which reads [`PlaneName::narrow`] and so never sees it.
///
/// `None` for a zero determinant only: parallel or line-sharing planes, no unique point.
pub fn three_planes_big(p: [&PlaneName; 3]) -> Option<MeetPoint> {
    use num_bigint::BigInt;
    use num_integer::Integer;

    let row = |name: &PlaneName| -> [BigInt; 4] {
        match name {
            PlaneName::Wide(c) => c.clone(),
            PlaneName::Narrow(c) => {
                let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(c[i].denom()));
                let l = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
                core::array::from_fn(|i| BigInt::from(c[i].numer()) * (&l / &den[i]))
            }
        }
    };
    three_planes_int([row(p[0]), row(p[1]), row(p[2])])
}

/// [`three_planes_big`]'s integer core — Cramer over integer plane rows, one reduction per
/// coordinate, the [`MeetPoint`] fork at the end. Separate so [`three_planes_rat`]'s fallback
/// can reach it with rows that are *not* canonical names (wrapping those in
/// [`PlaneName::Narrow`] would make the type say something false — that variant means "a
/// canonical name", and [`PlaneName::coeff_ints`] asserts its denominators are 1).
fn three_planes_int(m: [[num_bigint::BigInt; 4]; 3]) -> Option<MeetPoint> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{Signed, ToPrimitive, Zero};

    let det3 = |a: &[[BigInt; 3]; 3]| -> BigInt {
        let minor = |r0: usize, r1: usize, c0: usize, c1: usize| -> BigInt {
            &a[r0][c0] * &a[r1][c1] - &a[r0][c1] * &a[r1][c0]
        };
        &a[0][0] * minor(1, 2, 1, 2) - &a[0][1] * minor(1, 2, 0, 2) + &a[0][2] * minor(1, 2, 0, 1)
    };

    let base: [[BigInt; 3]; 3] =
        core::array::from_fn(|r| core::array::from_fn(|c| m[r][c].clone()));
    let rhs: [BigInt; 3] = core::array::from_fn(|r| -&m[r][3]);

    let d = det3(&base);
    if d.is_zero() {
        return None; // no unique point — parallel or line-sharing planes
    }

    // Cramer, then **one reduction per coordinate**: `d` is nonzero, so every gcd here is at
    // least 1 and the division is total.
    let out: [(BigInt, BigInt); 3] = core::array::from_fn(|i| {
        let mut a = base.clone();
        for (r, x) in rhs.iter().enumerate() {
            a[r][i] = x.clone();
        }
        let (mut n, mut q) = (det3(&a), d.clone());
        let g = n.gcd(&q);
        n /= &g;
        q /= &g;
        if q.is_negative() {
            n = -n;
            q = -q;
        }
        (n, q)
    });

    // ★ The same normalization invariant `plane_name_big` states: narrow whenever the canonical
    // answer fits `i128`, `Wide` only when it does not — so equal points are structurally equal
    // whichever route derived them.
    let narrow = (|| {
        let mut r = [Rat::from_int(0); 3];
        for (o, (n, q)) in r.iter_mut().zip(&out) {
            *o = Rat::new(n.to_i128()?, q.to_i128()?)?;
        }
        Some(r)
    })();
    Some(match narrow {
        Some(r) => MeetPoint::Narrow(r),
        None => MeetPoint::Wide(out),
    })
}

/// **A plane pushed `t` along its own unit normal**, exactly — `d′ = d − t·|n|`.
///
/// What a prism's far cap *is*: the face it was raised from, moved out by the sweep. Deriving it
/// this way rather than from the realized geometry is what lets two prisms raised to one height
/// record **one plane** — `7.7` in a single step and `1.1` then `6.6` give the same `d′`, because
/// `11/10 + 66/10` is `77/10` in rationals — and what lets the far cap of a face on a turned solid
/// be exact at all: the coefficients live in that face's own pre-motion frame, where the rotation's
/// irrational numbers never appear.
///
/// ★★ **`None` unless `|n|` is rational**, which is the one thing that can stop this: `n·n` must be
/// a perfect square. It is `1` for a box face and `25` for a `3-4-5` normal; it is `3` for `[1,1,1]`,
/// and there the caller keeps whatever path it had. Also `None` on `i128` overflow and for `n = 0`,
/// which is not a plane.
///
/// The result is canonicalized, so it compares by `==` with any other exact description of the same
/// plane — including the one the sketch frame produces when it is exact, which is what makes the
/// two derivations checkable against each other.
pub fn plane_offset(coeffs: [Rat; 4], t: Rat) -> Option<[Rat; 4]> {
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let mut nn = Rat::from_int(0);
    for c in n {
        nn = nn.checked_add(c.checked_mul(c)?)?;
    }
    if nn == Rat::from_int(0) {
        return None; // not a plane
    }
    let len = rat_sqrt_exact(nn)?;
    canonical_plane_coeffs([
        n[0],
        n[1],
        n[2],
        coeffs[3].checked_sub(t.checked_mul(len)?)?,
    ])
}
