use super::*;
/// The canonical rational plane with normal `n` through `point` — `n·x − n·point = 0`.
///
/// See [`canonical_plane_coeffs`] for why the result is canonicalized and why the inputs must be
/// rational by construction rather than lifted from f64 coefficients. `None` on `i128` overflow.
pub fn plane_from_point_normal(n: [Rat; 3], point: [Rat; 3]) -> Option<[Rat; 4]> {
    let mut d = Rat::from_int(0);
    for i in 0..3 {
        d = d.checked_sub(n[i].checked_mul(point[i])?)?;
    }
    canonical_plane_coeffs([n[0], n[1], n[2], d])
}

/// The canonical rational plane through three points, normal `(b − a) × (c − a)` — the exact twin
/// of `nacre_geom::Plane::through_points`, so the two describe the same plane with the same
/// orientation. `None` if the points are collinear (no plane) or on `i128` overflow.
pub fn plane_through_points(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<[Rat; 4]> {
    let d = |p: [Rat; 3], q: [Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            p[0].checked_sub(q[0])?,
            p[1].checked_sub(q[1])?,
            p[2].checked_sub(q[2])?,
        ])
    };
    let (u, v) = (d(b, a)?, d(c, a)?);
    let term = |i: usize, j: usize| u[i].checked_mul(v[j])?.checked_sub(u[j].checked_mul(v[i])?);
    let n = [term(1, 2)?, term(2, 0)?, term(0, 1)?];
    if n.iter().all(|c| *c == Rat::from_int(0)) {
        return None; // collinear
    }
    plane_from_point_normal(n, a)
}

/// The canonical name of a rational plane — primitive integer coefficients, sign-fixed
/// (see [`canonical_plane_coeffs`] for the canonical form). One vessel, two widths.
///
/// ★ **Normalization invariant: a value that fits `i128` is ALWAYS stored `Narrow`** — the only
/// constructor ([`plane_name_exact`]) enforces it, so two statements of one plane are
/// structurally equal (`==`/`Hash`) across representations. That is what makes this the
/// interning key: identity never depends on which route derived the name.
///
/// ★★ `Wide` carries **identity only**. Arithmetic shortcuts (frames, Shewchuk transports,
/// `base_rat`) read [`PlaneName::narrow`] and decline on `None`, exactly as they declined on a
/// missing name before — a wide name does not open a frame.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum PlaneName {
    Narrow([Rat; 4]),
    Wide([num_bigint::BigInt; 4]),
}

impl PlaneName {
    /// The `i128` form when there is one — what the exact shortcuts (Shewchuk integer
    /// predicates, frame derivation, `Isometry` transport) consume. `None` (`Wide`) means those
    /// shortcuts decline and judgment takes the general route: slower, never wrong.
    #[inline]
    pub fn narrow(&self) -> Option<&[Rat; 4]> {
        match self {
            PlaneName::Narrow(c) => Some(c),
            PlaneName::Wide(_) => None,
        }
    }

    /// The canonical coefficients as integers, **whichever width holds them** — the one vessel
    /// the integer sign predicates ([`int_plane_side`], [`int_cmp_coord`], [`int_dir_sign`])
    /// consume. A canonical narrow name is a primitive integer vector (the invariant
    /// [`plane_name_exact`] normalizes to), so `numer()` is the value; a wide name already is
    /// the integers.
    ///
    /// ⚠ **Primitive over *four* coefficients — the normal alone is not**.
    /// The content divided out is `gcd(a,b,c,d)`, so `(a,b,c)` keeps a factor of
    /// `gcd(a,b,c)/gcd(a,b,c,d)`; for an axis-aligned plane at an offset needing a long decimal
    /// that factor is the offset's **denominator**, arbitrarily large on the plainest of planes.
    /// A consumer that reads only the normal **and multiplies** must divide it by its own gcd
    /// first, or it pays that factor — squared, if it takes a cross product.
    pub fn coeff_ints(&self) -> [num_bigint::BigInt; 4] {
        match self {
            PlaneName::Narrow(c) => {
                debug_assert!(c.iter().all(|r| r.denom() == 1));
                core::array::from_fn(|k| num_bigint::BigInt::from(c[k].numer()))
            }
            PlaneName::Wide(c) => c.clone(),
        }
    }
}

/// **The plane three points name, computed so that the arithmetic on the way cannot lose it.**
///
/// The same answer [`plane_through_points`] gives, in a vessel that always holds it: `Narrow`
/// when the canonical answer fits `i128`, `Wide` (arbitrary-precision integers) when it does
/// not. `None` means exactly one thing — the points are collinear and name no plane. It is
/// never a shrug about an intermediate, nor about the answer's width.
///
/// ★★★★★ **The distinction is not academic — it was most of the failures.** `plane_through_points`
/// works in `Rat`, so `(b − a) × (c − a)` multiplies the points' denominators together and
/// `−n · a` multiplies once more. Measured on the failing cases, that peak needs **271 bits** while
/// the canonical answer, after the content is divided out, comes back to about **50**. Recomputing
/// those at unbounded precision reproduced the stored name **1,545 times out of 1,545** — the two
/// vectors differ by a scalar factor and `canonical_plane_coeffs` removes exactly that freedom. So
/// what overflowed was the road, not the destination.
///
/// ★★★ **Why a second function rather than widening the first.** `plane_through_points` is what a
/// *caller* uses to state a plane they wrote down ([`crate::plane_frame`]'s callers, a sketch, a
/// box's corners); this is what the kernel uses to **derive** the name of a plane it already holds
/// three exact points for. Widening the shared one would change both at once, and they are
/// different propositions with different evidence. (Making the caller's path exact too is worth
/// doing and is its own measurement.)
///
/// ★ **Cost is paid only on the fallback.** The `Rat` route runs first and is the answer whenever
/// it fits; `BigInt` is reached on the rest. Nothing here is on a boolean's inner loop — a plane is
/// named once per `Model::push_plane`.
pub fn plane_name_exact(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<PlaneName> {
    plane_through_points(a, b, c)
        .map(PlaneName::Narrow)
        .or_else(|| plane_name_big(a, b, c))
}

/// [`plane_name_exact`]'s unbounded arm, always taken — the differential test needs to call it on
/// inputs the `Rat` route handles, which it cannot do through the filter.
///
/// **The same plane [`plane_through_points`] computes**, reached by clearing each point's
/// denominators first so the arithmetic is integer throughout. The two agreeing wherever the narrow
/// one answers is the correctness argument, and `the_wide_derivation_answers_what_the_narrow_one_does`
/// is what holds it. `None` is collinearity, or a canonical component that does not fit `Rat`.
pub(crate) fn plane_name_big(a: [Rat; 3], b: [Rat; 3], c: [Rat; 3]) -> Option<PlaneName> {
    // ★★★★ **Integers, not rationals.** The obvious spelling is `Ratio<BigInt>`, mirroring
    // `plane_through_points` term for term — and that is how this started. But `Ratio` reduces by a
    // gcd on *every* multiply and subtract, and there are a dozen of them, so the reduction work
    // dominates. Clearing each point's denominators once up front leaves plain integer arithmetic
    // and exactly one gcd at the end, where the content has to come out anyway.
    //
    // ★ **A plane is scale-free, which is what makes this legal.** Each point is scaled by its own
    // `D_k`, so `u'` and `v'` are the true edges times `D_a·D_b` and `D_a·D_c`; their cross product
    // is the true normal times a positive factor, and the canonical form divides all of it out.
    // The `d` term is `−N·a = −(N·P_a)/D_a`, so scaling the whole 4-vector by `D_a` clears it.
    plane_name_from_lifted([a, b, c].map(|p| lift_meet(&MeetPoint::Narrow(p))))
}

/// [`plane_name_big`]'s body after the lift — three points as `(integer coordinates, positive
/// denominator scale)` pairs. Separate so [`plane_name_from_meets`] can reach it with points
/// that never were `Rat` triples ([`MeetPoint::Wide`]).
fn plane_name_from_lifted(
    pts: [([num_bigint::BigInt; 3], num_bigint::BigInt); 3],
) -> Option<PlaneName> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{ToPrimitive, Zero};

    let n = lifted_normal(&pts)?; // `None`: collinear
    let [(pa, da), _, _] = pts;
    let dot: BigInt = (0..3).map(|i| &n[i] * &pa[i]).sum();
    let mut num = [&n[0] * &da, &n[1] * &da, &n[2] * &da, -dot];

    // Divide out the content and fix the sign of the first nonzero — `canonical_plane_coeffs`'
    // steps ② and ③; ① is already done, since these are integers.
    let g = num.iter().fold(BigInt::zero(), |g, x| g.gcd(x));
    if g.is_zero() {
        return None; // the zero vector is not a plane
    }
    for x in &mut num {
        *x /= &g;
    }
    if num
        .iter()
        .find(|x| !x.is_zero())
        .is_some_and(|x| *x < BigInt::zero())
    {
        for x in &mut num {
            *x = -&*x;
        }
    }

    // ★ The normalization invariant lives here: narrow whenever the canonical answer fits
    // `i128`, `Wide` only when it does not — so equal planes are structurally equal whichever
    // route derived them. (This used to be the one honest failure; now it is the fork.)
    let mut out = [Rat::from_int(0); 4];
    for (o, x) in out.iter_mut().zip(&num) {
        match x.to_i128() {
            Some(v) => *o = Rat::from_int(v),
            None => return Some(PlaneName::Wide(num)),
        }
    }
    Some(PlaneName::Narrow(out))
}

/// **The canonical name of the plane through three meeting points, at whatever width the points
/// needed** — what lets a datum through [`MeetPoint::Wide`] vertices keep a name.
/// The same lift-and-join [`plane_name_big`] performs, with each point's
/// denominators cleared from whichever vessel holds them; a plane is scale-free per point, so the
/// per-point scale cannot move the canonical answer.
///
/// `None` means the points are collinear and name no plane — width is never a cause here, on
/// either side: wide points are lifted the same way, and a canonical answer too wide for `Rat`
/// comes back [`PlaneName::Wide`].
pub fn plane_name_from_meets(points: [&MeetPoint; 3]) -> Option<PlaneName> {
    plane_name_from_lifted(points.map(lift_meet))
}
/// One point as integer coordinates over a positive common denominator — the lift every
/// width-free road in this file starts from (`Narrow` and `Wide` meets alike).
fn lift_meet(p: &MeetPoint) -> ([num_bigint::BigInt; 3], num_bigint::BigInt) {
    use num_bigint::BigInt;
    use num_integer::Integer;
    match p {
        MeetPoint::Narrow(p) => {
            let den = p.map(|r| BigInt::from(r.denom()));
            let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
            let num = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&d / &den[i]));
            (num, d)
        }
        MeetPoint::Wide(p) => {
            let d = p.iter().fold(BigInt::from(1), |l, (_, den)| l.lcm(den));
            let num = core::array::from_fn(|i| &p[i].0 * (&d / &p[i].1));
            (num, d)
        }
    }
}

/// `(b − a) × (c − a)` of lifted points, times the positive factor `D_a²·D_b·D_c`; `None` when
/// the three are collinear.
fn lifted_normal(
    pts: &[([num_bigint::BigInt; 3], num_bigint::BigInt); 3],
) -> Option<[num_bigint::BigInt; 3]> {
    use num_bigint::BigInt;
    use num_traits::Zero;
    let [(pa, da), (pb, db), (pc, dc)] = pts;
    let edge = |q: &[BigInt; 3], dq: &BigInt| -> [BigInt; 3] {
        core::array::from_fn(|i| &q[i] * da - &pa[i] * dq)
    };
    let (u, v) = (edge(pb, db), edge(pc, dc));
    let term = |i: usize, j: usize| &u[i] * &v[j] - &u[j] * &v[i];
    let n = [term(1, 2), term(2, 0), term(0, 1)];
    (!n.iter().all(Zero::is_zero)).then_some(n)
}

/// **The direction three points span, exactly** — `(b − a) × (c − a)` up to a positive factor,
/// at whatever width the points need. `None` when they are collinear.
///
/// The factor is why this is a direction and not a vector: it answers *which way*, which is all
/// a plane's sense asks, and it never overflows because nothing here is narrowed.
pub fn triple_normal(points: [&MeetPoint; 3]) -> Option<[num_bigint::BigInt; 3]> {
    lifted_normal(&points.map(lift_meet))
}

/// Which way `n` faces `toward` — the exact sign of `n · toward`.
pub fn normal_sense(n: &[num_bigint::BigInt; 3], toward: [Rat; 3]) -> Orient {
    let (t, _) = lift_meet(&MeetPoint::Narrow(toward));
    dot_sign(n, &t)
}

/// Whether two normals of **one plane** point the same way. Two such normals are parallel, so the
/// dot product is never zero and its sign is the whole answer.
pub fn same_sense(a: &[num_bigint::BigInt; 3], b: &[num_bigint::BigInt; 3]) -> bool {
    let s = dot_sign(a, b);
    debug_assert!(s != Orient::Zero, "two normals of one plane are parallel");
    s == Orient::Positive
}

/// **Whether a plane's name points the way its defining points turn** — the name's normal against
/// `(b − a) × (c − a)`, exactly and at any width. `None` when the points are collinear.
///
/// The canonical name carries no direction (its first nonzero component is positive), so this is
/// the comparison that recovers one against a statement. The points must be the plane's own —
/// on the named plane, in the frame the name speaks — so the two normals are parallel and the
/// answer is never a near call.
pub fn name_along_points(name: &PlaneName, points: [&MeetPoint; 3]) -> Option<bool> {
    // An interval filter, then the `Rat` route, then integers that cannot overflow — the judge's
    // shape and `plane_residual_sign`'s. The two normals are parallel, so the dot is a full
    // magnitude from zero and the filter answers unless the points are nearly collinear; every
    // radius is a sound bound (`to_f64` rounds to nearest — half an ulp, and the smallest normal
    // below that), so a filter answer is the exact answer.
    if let (
        PlaneName::Narrow(c),
        [
            MeetPoint::Narrow(a),
            MeetPoint::Narrow(b),
            MeetPoint::Narrow(p),
        ],
    ) = (name, points)
    {
        let iv = |r: Rat| {
            let x = r.to_f64();
            crate::Bounded::new(x, (x.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE))
        };
        let e = |q: &[Rat; 3]| [0, 1, 2].map(|i| iv(q[i]).sub(iv(a[i])));
        let (u, v) = (e(b), e(p));
        let t = |i: usize, j: usize| u[i].mul(v[j]).sub(u[j].mul(v[i]));
        let n = [t(1, 2), t(2, 0), t(0, 1)];
        let dot = (0..3).fold(crate::Bounded::new(0.0, 0.0), |acc, k| {
            acc.add(iv(c[k]).mul(n[k]))
        });
        if let Some(positive) = dot.sign() {
            return Some(positive);
        }
        let narrow = || -> Option<Option<bool>> {
            let e = |q: &[Rat; 3]| -> Option<[Rat; 3]> {
                Some([
                    q[0].checked_sub(a[0])?,
                    q[1].checked_sub(a[1])?,
                    q[2].checked_sub(a[2])?,
                ])
            };
            let (u, v) = (e(b)?, e(p)?);
            let t =
                |i: usize, j: usize| u[i].checked_mul(v[j])?.checked_sub(u[j].checked_mul(v[i])?);
            let n = [t(1, 2)?, t(2, 0)?, t(0, 1)?];
            let mut dot = Rat::from_int(0);
            for k in 0..3 {
                dot = dot.checked_add(c[k].checked_mul(n[k])?)?;
            }
            let zero = Rat::from_int(0);
            Some(if n.iter().all(|x| *x == zero) {
                None
            } else {
                debug_assert!(
                    dot != zero,
                    "a name's normal is parallel to its own points' turn"
                );
                Some(dot > zero)
            })
        };
        if let Some(answer) = narrow() {
            return answer;
        }
    }
    let n = triple_normal(points)?;
    let [a, b, c, _] = name.coeff_ints();
    Some(same_sense(&[a, b, c], &n))
}

fn dot_sign(a: &[num_bigint::BigInt; 3], b: &[num_bigint::BigInt; 3]) -> Orient {
    use num_traits::Zero;
    let d: num_bigint::BigInt = (0..3).map(|i| &a[i] * &b[i]).sum();
    if d.is_zero() {
        Orient::Zero
    } else if d > num_bigint::BigInt::zero() {
        Orient::Positive
    } else {
        Orient::Negative
    }
}

/// **The canonical representative of a rational plane** `a·x + b·y + c·z + d = 0`.
///
/// A plane is scale-invariant — `2x + 4y − 6 = 0` and `x + 2y − 3 = 0` are the same plane — so a
/// family of coefficient vectors describes it and equality has to pick one of them. Three steps:
/// clear the denominators, divide out the content, and fix the sign of the first nonzero component.
/// What comes back is a **primitive integer vector**, and two vectors describing the same plane
/// canonicalize to **bit-identical** arrays. That is what turns *"are these the same plane?"* from a
/// predicate into `==`.
///
/// ★★ **Scale-independence is the whole point, and it is what f64 coefficients cannot give.**
/// `nacre_geom::Plane` stores an un-normalized `raw` normal whose length follows the *face's size*,
/// so two faces of the plane `x = 3` come out as `[2.2, 0, 0, −6.6000000000000005]` and
/// `[13.2, 0, 0, −39.599999999999994]` — the same plane, not exactly proportional, because `d` was a
/// rounded product. Measured: 18 pairs in the census are merged only because a *second* test looks at
/// the faces' coordinates instead.
///
/// ★★★ **The input must be rational by construction, not lifted from those f64 coefficients.**
/// Lifting is lossless but it preserves the drift: `2.2` and `6.6000000000000005` are different
/// dyadics whose exact ratio is not 3, so canonicalizing them still gives two different vectors.
/// Rationals here come from the dimensions the user *wrote* ([`Rat::from_decimal`]) carried through
/// exact arithmetic — the same rule the construction path already follows.
///
/// `None` on `i128` overflow (the denominators' lcm, or a numerator scaled by it), which is the
/// kernel's ordinary demotion signal: the caller keeps the plain rational and the geometric tests
/// answer instead. Nothing is wrong, one shortcut is unavailable.
///
/// All-zero coefficients are not a plane; they canonicalize to themselves.
pub fn canonical_plane_coeffs(coeffs: [Rat; 4]) -> Option<[Rat; 4]> {
    // ① Clear the denominators: multiply through by their lcm.
    let mut lcm: i128 = 1;
    for c in coeffs {
        let d = c.denom();
        let g = gcd_u128(lcm.unsigned_abs(), d.unsigned_abs()) as i128;
        lcm = lcm.checked_div(g)?.checked_mul(d)?;
    }
    let mut num = [0i128; 4];
    for (i, c) in coeffs.iter().enumerate() {
        // Exact: `lcm` is a multiple of every denominator, so the division has no remainder.
        num[i] = c.numer().checked_mul(lcm.checked_div(c.denom())?)?;
    }

    // ② Divide out the content.
    let g = num.iter().fold(0u128, |g, n| gcd_u128(g, n.unsigned_abs()));
    if g == 0 {
        return Some(coeffs); // the zero vector is not a plane
    }
    let g = g as i128;
    for n in &mut num {
        *n /= g;
    }

    // ③ Fix the sign: the first nonzero component is positive.
    if num.iter().find(|n| **n != 0).is_some_and(|n| *n < 0) {
        for n in &mut num {
            *n = n.checked_neg()?;
        }
    }
    Some(num.map(Rat::from_int))
}
