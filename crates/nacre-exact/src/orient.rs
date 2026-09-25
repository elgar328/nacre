use super::*;
/// The sign of the 2-D orientation determinant `(b − a) × (c − a)` — positive when `abc` turns
/// counter-clockwise, zero when collinear. The exact-`Rat` twin of `nacre_predicates::orient2d`
/// (same convention), for coordinates that *are* rational rather than f64.
///
/// ★ **Total.** Unlike every other `Rat` derivation here, this cannot decline: a sign consumed by
/// a simplicity check that failed open on overflow would read "no intersection" where there is
/// one — silent-wrong, the exact shape [`plane_name_exact`] was built to kill. So the `Rat` route
/// runs first and an overflow falls through to integer arithmetic that cannot overflow.
///
/// ★ **The `BigInt` arm is load-bearing, not a corner case.** [`Rat::from_decimal`] denominators
/// reach `10^38`, and the determinant multiplies two coordinate *differences* — two 17-digit
/// dimensions already push the product denominator past `i128`. Same denominator-clearing shape
/// as [`plane_name_big`]: each point is scaled by its own positive denominator, which multiplies
/// the determinant by `Da²·Db·Dc > 0` and therefore cannot move the sign.
pub fn orient2d_rat(a: [Rat; 2], b: [Rat; 2], c: [Rat; 2]) -> i8 {
    let narrow = || -> Option<i8> {
        let u = [b[0].checked_sub(a[0])?, b[1].checked_sub(a[1])?];
        let v = [c[0].checked_sub(a[0])?, c[1].checked_sub(a[1])?];
        let det = u[0]
            .checked_mul(v[1])?
            .checked_sub(u[1].checked_mul(v[0])?)?;
        Some(rat_sign(det))
    };
    narrow().unwrap_or_else(|| orient2d_big(a, b, c))
}

/// `-1 / 0 / +1` of a rational.
fn rat_sign(x: Rat) -> i8 {
    match x.cmp(&Rat::from_int(0)) {
        core::cmp::Ordering::Less => -1,
        core::cmp::Ordering::Equal => 0,
        core::cmp::Ordering::Greater => 1,
    }
}

/// [`orient2d_rat`]'s unbounded arm, always taken — the differential test needs to call it on
/// inputs the `Rat` route handles, which it cannot do through the filter (the 2-D analogue of
/// [`plane_name_big`]'s arrangement).
pub(crate) fn orient2d_big(a: [Rat; 2], b: [Rat; 2], c: [Rat; 2]) -> i8 {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;

    // Clear each point's denominators once (positive by `Ratio`'s reduction invariant), then the
    // arithmetic is plain integers and the determinant is the true one times `Da²·Db·Dc`.
    let lift = |p: [Rat; 2]| -> ([BigInt; 2], BigInt) {
        let den = p.map(|r| BigInt::from(r.denom()));
        let d = den[0].lcm(&den[1]);
        let num = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&d / &den[i]));
        (num, d)
    };
    let ((pa, da), (pb, db), (pc, dc)) = (lift(a), lift(b), lift(c));
    let edge = |q: &[BigInt; 2], dq: &BigInt| -> [BigInt; 2] {
        core::array::from_fn(|i| &q[i] * &da - &pa[i] * dq)
    };
    let (u, v) = (edge(&pb, &db), edge(&pc, &dc));
    match (&u[0] * &v[1] - &u[1] * &v[0]).sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// **The sign of a plane's residual at a rational point** — `sign(a·x + b·y + c·z + d)` for the
/// plane a [`PlaneName`] names. `0` means the point lies exactly on the plane.
///
/// ★ **Total**, like [`orient2d_rat`] and for the same reason: this is what a *validating
/// constructor* consumes (is a caller's stated sketch origin on the plane they picked?), and a
/// check that failed open on overflow would accept an off-plane origin — silent-wrong. The `Rat`
/// route runs first; an overflow, or a `Wide` name, falls through to integer arithmetic that
/// cannot overflow.
///
/// The sign convention is the name's own: which side is positive depends on how the canonical
/// coefficients came out, so callers should compare against `0`, not against each other across
/// planes.
pub fn plane_residual_sign(name: &PlaneName, p: [Rat; 3]) -> i8 {
    match name {
        PlaneName::Narrow(c) => {
            let narrow = || -> Option<i8> {
                let mut acc = c[3];
                for k in 0..3 {
                    acc = acc.checked_add(c[k].checked_mul(p[k])?)?;
                }
                Some(rat_sign(acc))
            };
            narrow().unwrap_or_else(|| {
                // A canonical narrow name is a primitive integer vector (the invariant
                // `plane_name_exact` normalizes to), so `numer()` is the value.
                debug_assert!(c.iter().all(|r| r.denom() == 1));
                residual_sign_big(&c.map(|x| num_bigint::BigInt::from(x.numer())), p)
            })
        }
        PlaneName::Wide(c) => residual_sign_big(c, p),
    }
}

/// [`plane_residual_sign`]'s unbounded arm — integer coefficients (which both name forms reduce
/// to), a rational point. Clearing the point's denominators multiplies the residual by a
/// positive factor, which cannot move the sign.
fn residual_sign_big(ci: &[num_bigint::BigInt; 4], p: [Rat; 3]) -> i8 {
    use num_bigint::{BigInt, Sign};
    use num_integer::Integer;
    let den = p.map(|r| BigInt::from(r.denom()));
    let l = den.iter().fold(BigInt::from(1), |acc, x| acc.lcm(x));
    let num: [BigInt; 3] = core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&l / &den[i]));
    let res: BigInt = (0..3).map(|i| &ci[i] * &num[i]).sum::<BigInt>() + &ci[3] * &l;
    match res.sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}

/// **Whether two rational directions are parallel** — `a × b = 0`, exactly and always.
///
/// The cross is homogeneous in each argument, so each vector's own denominator is a positive
/// factor of the result and clearing it cannot move the zero test.
///
/// **A zero vector is parallel to everything**, this predicate included: `0 × b = 0`. Callers
/// that mean "these two span a plane" therefore get the answer they want — a zero direction
/// spans nothing — without a separate zero check.
pub fn parallel_rat(a: &[Rat; 3], b: &[Rat; 3]) -> bool {
    use num_traits::Zero;
    let (x, _) = lift3(a);
    let (y, _) = lift3(b);
    let term = |i: usize, j: usize| &x[i] * &y[j] - &x[j] * &y[i];
    term(1, 2).is_zero() && term(2, 0).is_zero() && term(0, 1).is_zero()
}

/// **The sign of `a · b`** for rational vectors, exactly and always — [`Orient::Zero`] exactly
/// when they are perpendicular. Bilinear, so each vector's positive denominator factors out of
/// the answer and the lift's scales are dropped.
pub fn dot_sign_rat(a: &[Rat; 3], b: &[Rat; 3]) -> Orient {
    let (x, _) = lift3(a);
    let (y, _) = lift3(b);
    orient_of(big_sign(
        &(0..3).map(|i| &x[i] * &y[i]).sum::<num_bigint::BigInt>(),
    ))
}

/// **How a plane stands to a cylinder's axis** — the one relation a plane × cylinder section
/// depends on, and so what kind of curve an edge on that pair is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AxisRelation {
    /// The plane crosses the axis squarely (normal ∥ axis): the section is a circle.
    Across,
    /// The plane runs along the axis (normal ⊥ axis): the section is rulings — straight lines.
    Along,
    /// Neither: the section is an ellipse.
    Oblique,
}

/// [`AxisRelation`] of a plane with normal `n` to an axis `d`, over integer vectors — exact at any
/// width, and scale-free: each vector's positive scale, and its sign, leave the answer alone.
///
/// ⚠ A zero vector is parallel to everything ([`parallel_rat`]'s convention), so it answers
/// `Across`; a statement's own constructor refuses a zero normal or axis before one gets here.
fn axis_relation_int(n: &[num_bigint::BigInt; 3], d: &[num_bigint::BigInt; 3]) -> AxisRelation {
    use num_traits::Zero;
    let term = |i: usize, j: usize| &n[i] * &d[j] - &n[j] * &d[i];
    if term(1, 2).is_zero() && term(2, 0).is_zero() && term(0, 1).is_zero() {
        return AxisRelation::Across;
    }
    let dot: num_bigint::BigInt = (0..3).map(|i| &n[i] * &d[i]).sum();
    if dot.is_zero() {
        AxisRelation::Along
    } else {
        AxisRelation::Oblique
    }
}

/// [`AxisRelation`] of a plane with normal `n` to an axis `d` — each lifted to integers by its own
/// positive denominator, which the relation cannot see.
pub fn axis_relation(n: &[Rat; 3], d: &[Rat; 3]) -> AxisRelation {
    axis_relation_int(&lift3(n).0, &lift3(d).0)
}

/// **A direction's line, in integers** — no sign and no scale, which is all a plane–axis
/// relation reads, carried through the motions a chain records. Exact at any width: a `Wide`
/// plane name's normal is one too.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LineDir([num_bigint::BigInt; 3]);

impl LineDir {
    /// The line a rational direction spans.
    pub fn of_rat(d: &[Rat; 3]) -> LineDir {
        LineDir(lift3(d).0)
    }

    /// The line of a plane name's normal.
    pub fn normal_of(name: &crate::PlaneName) -> LineDir {
        let [a, b, c, _] = name.coeff_ints();
        LineDir([a, b, c])
    }

    /// The line of `ẑ` — a sketch frame's normal, in the frame.
    pub fn z() -> LineDir {
        LineDir([0.into(), 0.into(), 1.into()])
    }

    /// Whether this is the line of `ẑ`.
    pub fn is_z(&self) -> bool {
        use num_traits::Zero;
        self.0[0].is_zero() && self.0[1].is_zero() && !self.0[2].is_zero()
    }

    /// The same line after (or before — a mirror is its own inverse) a reflection across `axis`.
    pub fn mirrored(&self, axis: Axis) -> LineDir {
        let mut out = self.0.clone();
        let k = axis_index(axis);
        out[k] = -out[k].clone();
        LineDir(out)
    }

    /// The line after the turn `r` — or before it, `back` — when that is a rational line: a
    /// quarter turn permutes it, and any turn leaves its own axis alone. `None` otherwise.
    pub fn turned(&self, r: Rotation, back: bool) -> Option<LineDir> {
        use num_traits::Zero;
        match crate::AxisAffine::rotation(r) {
            Some(f) if back => Some(LineDir(f.inverse().dir_int(&self.0))),
            Some(f) => Some(LineDir(f.dir_int(&self.0))),
            None => {
                let k = axis_index(r.axis);
                (0..3)
                    .all(|i| i == k || self.0[i].is_zero())
                    .then(|| self.clone())
            }
        }
    }

    /// How a plane with this normal stands to the axis `axis` ([`AxisRelation`]).
    pub fn relation_to(&self, axis: &LineDir) -> AxisRelation {
        axis_relation_int(&self.0, &axis.0)
    }
}

fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
}

/// `big_sign`'s answer as the geometry's [`Orient`].
pub(super) fn orient_of(s: i8) -> Orient {
    match s {
        1 => Orient::Positive,
        -1 => Orient::Negative,
        _ => Orient::Zero,
    }
}

/// **Which side of the plane `j` the implicit point `∩(p₀, p₁, p₂)` lies on**, over integer
/// coefficients — `sign(j·[Dvec : D]) · sign(D)`, the integer twin of
/// `nacre_predicates::indirect_plane_side` (same convention: `+1` on the side `j`'s normal
/// points to, and the sign is the plane's own, so a face whose stored normal opposes its
/// outward direction applies that relation itself).
///
/// Row-linear in each input, so any per-row **positive** scale — clearing a denominator,
/// un-reducing a canonical vector — cannot move the answer; negating a `p` row cannot either
/// (it negates `D` and the dot together), while negating `j` negates the answer. `D = 0`
/// (the three planes meet in no point) returns `0` with no branch, as in the twin.
pub fn int_plane_side(p: [&[num_bigint::BigInt; 4]; 3], j: &[num_bigint::BigInt; 4]) -> i8 {
    let (dvec, d) = cramer_big(p);
    let side: num_bigint::BigInt =
        (0..3).map(|i| &j[i] * &dvec[i]).sum::<num_bigint::BigInt>() + &j[3] * &d;
    big_sign(&side) * big_sign(&d)
}

/// The exact sign of `a[axis] − b[axis]` for two implicit points over integer coefficients —
/// `sign(Na·Db − Nb·Da) · sign(Da) · sign(Db)`, the integer twin of
/// `nacre_predicates::indirect_cmp_coord`. Orientation-invariant: negating any row negates a
/// numerator and its denominator together.
///
/// **Precondition:** both triples meet in a point (`Da, Db ≠ 0`), as in the twin.
pub fn int_cmp_coord(
    a: [&[num_bigint::BigInt; 4]; 3],
    b: [&[num_bigint::BigInt; 4]; 3],
    axis: usize,
) -> i8 {
    debug_assert!(axis < 3, "int_cmp_coord: axis must be 0, 1 or 2");
    let (na, da) = cramer_big(a);
    let (nb, db) = cramer_big(b);
    debug_assert!(
        big_sign(&da) != 0 && big_sign(&db) != 0,
        "int_cmp_coord: degenerate three-plane input (D = 0)"
    );
    big_sign(&(&na[axis] * &db - &nb[axis] * &da)) * big_sign(&da) * big_sign(&db)
}

/// `x·y` in checked [`Rat`] — `None` is overflow, which every reader takes as "not stated".
///
/// The value twin of [`dot_sign_rat`], which answers only the sign and never overflows because it
/// lifts to `BigInt`. A caller that needs the number takes this and handles the `None`.
pub fn dot3_rat(x: &[Rat; 3], y: &[Rat; 3]) -> Option<Rat> {
    x[0].checked_mul(y[0])?
        .checked_add(x[1].checked_mul(y[1])?)?
        .checked_add(x[2].checked_mul(y[2])?)
}

/// `x × y` in checked [`Rat`] — `None` is overflow.
///
/// Exact and in ℚ, so the result is a normal of the plane the two vectors span whenever they are
/// independent; a zero vector back means they are parallel, which [`parallel_rat`] answers
/// without the multiplications.
pub fn cross3_rat(x: &[Rat; 3], y: &[Rat; 3]) -> Option<[Rat; 3]> {
    Some([
        x[1].checked_mul(y[2])?
            .checked_sub(x[2].checked_mul(y[1])?)?,
        x[2].checked_mul(y[0])?
            .checked_sub(x[0].checked_mul(y[2])?)?,
        x[0].checked_mul(y[1])?
            .checked_sub(x[1].checked_mul(y[0])?)?,
    ])
}

/// `sign(det[n₀; n₁; n₂])` over the three planes' integer normals — how the line `p ∩ a` runs
/// relative to plane `b`, the integer twin of `nacre_predicates::det3_sign` on plane normals
/// (`d` never enters). Direction-sensitive in every row: negating one normal negates the
/// answer.
pub fn int_dir_sign(n: [&[num_bigint::BigInt; 4]; 3]) -> i8 {
    big_sign(&det3_big(n.map(|r| [&r[0], &r[1], &r[2]])))
}

#[cfg(test)]
#[path = "tests/orient.rs"]
mod tests;
