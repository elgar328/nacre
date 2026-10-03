use super::*;
/// **A cylinder's rational unit cross-section frame** — `û₁ = ref_dir/‖ref_dir‖`,
/// `û₂ = (dir × ref_dir)/(‖dir‖‖ref_dir‖)`, both exact rationals ⊥ the axis and to each other.
///
/// `None` — and the three causes are told apart on purpose:
/// * `ref_dir` is **parallel to the axis** (a zero vector included — [`parallel_rat`] says so):
///   the statement pins no cross-section at all;
/// * `ref_dir` is **not perpendicular** to the axis ([`dot_sign_rat`], exact and total). The
///   general recipe would project it, and that is deliberately **not** done here: a caller who
///   needs the frame of a `CylinderDef` gets a frame whose `û·dir = 0` is a *fact*, not a
///   derivation, because everything downstream (`quad::cylinder_radial_side` reduces to
///   `‖dir‖²·r²·(û·û − 1)`) is sound only when that holds exactly.
/// * a norm is **irrational**, or the arithmetic left `i128` ([`inv_sqrt_exact`]).
///
/// ★ Every cylinder the modelling road states satisfies the first two by construction — a sketch
/// frame is checked orthonormal exactly, and `ref_dir` is either its `x̂` or a rim chord divided
/// by its own radius — so `None` there means the third cause alone.
pub fn cyl_unit_frame(dir: &[Rat; 3], ref_dir: &[Rat; 3]) -> Option<([Rat; 3], [Rat; 3])> {
    if parallel_rat(ref_dir, dir) || dot_sign_rat(ref_dir, dir) != Orient::Zero {
        return None;
    }
    let dot = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<Rat> {
        let mut acc = Rat::from_int(0);
        for k in 0..3 {
            acc = acc.checked_add(a[k].checked_mul(b[k])?)?;
        }
        Some(acc)
    };
    let scale = |v: &[Rat; 3], k: Rat| -> Option<[Rat; 3]> {
        Some([
            v[0].checked_mul(k)?,
            v[1].checked_mul(k)?,
            v[2].checked_mul(k)?,
        ])
    };
    let cross = [
        dir[1]
            .checked_mul(ref_dir[2])?
            .checked_sub(dir[2].checked_mul(ref_dir[1])?)?,
        dir[2]
            .checked_mul(ref_dir[0])?
            .checked_sub(dir[0].checked_mul(ref_dir[2])?)?,
        dir[0]
            .checked_mul(ref_dir[1])?
            .checked_sub(dir[1].checked_mul(ref_dir[0])?)?,
    ];
    let inv_e = inv_sqrt_exact(dot(ref_dir, ref_dir)?)?;
    let inv_m = inv_sqrt_exact(dot(dir, dir)?)?;
    Some((
        scale(ref_dir, inv_e)?,
        scale(&cross, inv_e.checked_mul(inv_m)?)?,
    ))
}

/// **A cylinder's cache frame, correctly rounded** — the unit axis `dir/‖dir‖` and the unit
/// reference direction `e₁/‖e₁‖`, `e₁ = (dir·dir)·ref_dir − (ref_dir·dir)·dir` the part of `ref_dir`
/// perpendicular to the axis, each component the `f64` nearest the true value
/// ([`unit_vector_f64`]). The cache's twin of [`cyl_unit_frame`], which answers only where the
/// norms are rational; this one always does, because an `f64` names every unit vector.
///
/// Projected, not refused: the cache's `ref_dir` is the seam's direction across the axis, and a
/// statement whose `ref_dir` leans along the axis still names that direction. The two vectors are
/// lifted to integers first (a positive common denominator does not move a direction), so a wide
/// statement does not fall out of `Rat` on the way. `None` for a zero axis or a `ref_dir` parallel
/// to it.
pub fn cyl_unit_frame_f64(dir: &[Rat; 3], ref_dir: &[Rat; 3]) -> Option<([f64; 3], [f64; 3])> {
    use num_bigint::BigInt;
    let ints = |v: &[Rat; 3]| -> [BigInt; 3] {
        let den: BigInt = v.iter().map(|r| BigInt::from(r.denom())).product();
        (*v).map(|r| BigInt::from(r.numer()) * (&den / BigInt::from(r.denom())))
    };
    let (m, e) = (ints(dir), ints(ref_dir));
    let dot = |a: &[BigInt; 3], b: &[BigInt; 3]| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
    let (mm, em) = (dot(&m, &m), dot(&e, &m));
    let e1: [BigInt; 3] = core::array::from_fn(|k| &mm * &e[k] - &em * &m[k]);
    Some((unit_vector_f64(&m)?, unit_vector_f64(&e1)?))
}

/// **Where a plane meets the line `o + t·m`, as the parameter `t`** — `t = −(n·o + d)/(n·m)` for
/// the plane `n·x + d = 0`, in checked `Rat`. The one spelling of a cylinder axis meeting a cap:
/// every question of where a plane crosses an axis is this, and [`axis_plane_meet`] is the point
/// it names.
///
/// `None` when the plane is parallel to the line (`n·m = 0`, no meeting point) or the arithmetic
/// leaves `i128`.
pub fn axis_param_of_plane(coeffs: &[Rat; 4], o: &[Rat; 3], m: &[Rat; 3]) -> Option<Rat> {
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let nm = dot3_rat(&n, m)?;
    if nm == Rat::from_int(0) {
        return None;
    }
    let no_d = dot3_rat(&n, o)?.checked_add(coeffs[3])?;
    Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)
}

/// **The point where a plane meets the line `o + t·m`** — `o + t·m` at [`axis_param_of_plane`]'s
/// `t`, exact: a cap's circle centre on a cylinder's axis. `None` where that `t` is, or where a
/// coordinate leaves `i128`.
pub fn axis_plane_meet(coeffs: &[Rat; 4], o: &[Rat; 3], m: &[Rat; 3]) -> Option<[Rat; 3]> {
    let t = axis_param_of_plane(coeffs, o, m)?;
    let mut p = *o;
    for k in 0..3 {
        p[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(p)
}

/// **How a point's distance from a plane compares with `r`**, the radius stated as its square
/// `r2` — [`Orient::Negative`] inside the slab of half-width `r` about the plane, [`Orient::Zero`]
/// exactly at distance `r`, [`Orient::Positive`] clear of it. Exact and total.
///
/// `sign((n·p + d)² − r²|n|²)`, which is `sign(dist² − r²)` scaled by the positive `|n|²` — no
/// normalization and no square root. Stated in **plane** vocabulary on purpose: the cylinder
/// gate asks it about an axis point (a wall parallel to the axis is the same distance from every
/// point of it), but nothing here is about cylinders.
///
/// ★ **Not homogeneous in `p`** — the `d` term is why — so `p`'s denominator and the radius'
/// ride into the formula instead of dropping out: for `coeffs = C/Dc`, `p = P/Dp`, `r² = R/S`
/// the answer is `sign(S(C₀₋₂·P + C₃·Dp)² − R|C₀₋₂|²Dp²)` (`Dc²` *is* a positive common factor
/// and does cancel). Dropping `Dp` instead — the obvious spelling — states a different
/// proposition, and disagrees with the truth on 1.3% of mixed-denominator inputs (measured).
///
/// **Precondition:** `r2 ≥ 0`; a negative squared radius has no distance to compare with.
pub fn point_plane_clearance_rat(coeffs: &[Rat; 4], p: &[Rat; 3], r2: &BigRat) -> Orient {
    use num_bigint::BigInt;
    debug_assert!(
        !r2.is_negative(),
        "clearance compares against a non-negative squared radius"
    );
    let (c, _dc) = lift4(coeffs);
    let (pp, dp) = lift3(p);
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let dot: BigInt = (0..3).map(|i| &c[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp;
    let nn: BigInt = (0..3).map(|i| &c[i] * &c[i]).sum();
    orient_of(big_sign(&(&rd * (&dot * &dot) - &rn * nn * (&dp * &dp))))
}

/// **Does `p` satisfy these coefficients exactly?** — `n·p + c = 0`, at any width.
///
/// The consumer-side net for [`cylinder_strip_side`]'s on-plane precondition. A plane has two
/// exact descriptions — its coefficients and the triangle its points span — and they need not
/// agree; a face merged into a class by *rounded* coefficients can therefore have vertices that
/// do not satisfy the class root's exact name. Asking here turns that into a fact the caller can
/// refuse on, instead of a `debug_assert` that says nothing in a release build.
pub fn point_on_plane_exact(coeffs: &[Rat; 4], p: &MeetPoint) -> bool {
    use num_bigint::BigInt;
    let (c, _dc) = lift4(coeffs);
    let (pp, dp) = p.lift();
    // (C₀₋₂·P + C₃·Dp) / (Dc·Dp) — the denominators are positive, so only the numerator decides.
    (0..3).map(|i| &c[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp == BigInt::from(0)
}

/// Which side of the **strip** a cylinder cuts out of a plane parallel to its axis a point lies
/// on — see [`cylinder_strip_side`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StripSide {
    /// Clear of the strip, on the `+(n × m)` side.
    Plus,
    /// Clear of the strip, on the `−(n × m)` side.
    Minus,
    /// Inside the strip, boundary included — the point is within `r` of the axis.
    ///
    /// For a piece with **width** this is "reaches the strip but does not span a boundary":
    /// wholly between the two rulings, or touching one at a single point.
    Inside,
    /// **The piece spans a boundary** — its interior lies strictly on both sides of one ruling.
    ///
    /// ★ A point can never be this (it has no width), so every answer
    /// [`cylinder_strip_side`] and [`cylinder_strip_side_branch`] give is one of the three above,
    /// unchanged. Only [`cylinder_strip_side_margin`] with `rho > 0` produces it.
    Crosses,
}

/// **Where a point on a plane parallel to a cylinder's axis stands relative to the strip the
/// cylinder cuts out of that plane.** Exact and total, at any width.
///
/// A plane with `n · m = 0` meets the solid cylinder in a strip of half-width `h = √(r² − d²)`
/// (empty when the plane clears, `d ≥ r`), running along the axis direction. Decomposing
/// `p − o` into three **mutually orthogonal** parts — along `n` (magnitude `d`), along `m`, and
/// along `e = n × m` — gives `dist(p, axis)² = d² + t²` where `t` is the `e` component. So with
/// `U = (p − o) · e = t·|e|` and `|e| = |n||m|`:
///
/// ```text
/// clear of the strip  ⟺  U² > (r²|n|² − (n·o + c)²) · |m|²
/// ```
///
/// No normalization and **no square root** — `h` never has to be formed.
///
/// ★ **This is one of two axes, not a rule of its own.** What the cylinder occupies in this plane
/// is a *rectangle*: this strip across, and the lateral face's axis-parameter span along
/// ([`point_axis_side`]). Clearing either axis clears the rectangle, so a caller asking whether a
/// face misses the cylinder asks both and stops at the first that separates.
///
/// ★ **`U = 0` answers [`StripSide::Inside`] before the magnitude test.** With a non-empty strip
/// that is the truth (the point sits on the axis' own in-plane line). With an empty one it is
/// merely conservative — and a caller that cares has already passed the plane through
/// [`point_plane_clearance_rat`], which is the cheaper question and the one that decides
/// emptiness. That ordering makes this function total with no precondition to forget.
///
/// **Preconditions:** `r2 ≥ 0` (the radius stated as its square); `n · m = 0` (a plane that is not
/// parallel to the axis cuts a conic, not a strip); and **`p` lies on that plane** — the decomposition takes the perpendicular
/// distance from the *axis*, so an off-plane point would be judged against a distance that is not
/// its own. All three are `debug_assert`ed.
pub fn cylinder_strip_side(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    // ★ **A door, so the two can never drift.** A point is a disk of radius zero, and at that
    // radius the margin form's answers collapse onto this one exactly: `Crosses` needs a strictly
    // positive width to be possible at all, and the remaining comparison is term-for-term this
    // function's own.
    cylinder_strip_side_margin(coeffs, p, &BigRat::zero(), o, m, r2)
}

/// The three squared quantities the strip questions compare, in **one common positive scale** —
/// the whole derivation, shared by every door below.
///
/// With `e = n × m` (so `|e|² = |n|²|m|²`, the two being perpendicular):
///
/// ```text
/// U  = (p − o)·e        W² = (r²|n|² − (n·o + c)²)|m|²        ρ'² = ρ²|e|²
/// ```
///
/// A disk of radius `ρ` about `p` sweeps `U ± ρ'`, and the cylinder's two rulings stand at
/// `U = ±W`; every question below is a comparison among those three. The radii arrive **squared**
/// (`r2 = r²`, `rho2 = ρ²`, the form every radius takes in this family), so cleared of
/// denominators the scale is `dp²·d_o²·rd·sd` times the `1/(dc²dm²)` the plane's and direction's
/// own denominators contribute — positive throughout, so only signs survive.
struct StripScale {
    /// `sign(U)` — which side of the axis plane the centre is on.
    u_sign: Orient,
    /// `U²`, `W²` and `ρ'²` in the common scale. `ww` is negative exactly when the plane clears
    /// the cylinder altogether, and then there is no strip.
    uu: num_bigint::BigInt,
    ww: num_bigint::BigInt,
    rr: num_bigint::BigInt,
}

fn strip_scale(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripScale {
    use num_bigint::BigInt;
    debug_assert!(!r2.is_negative(), "a squared radius is not negative");
    debug_assert!(!rho2.is_negative(), "a squared margin is not negative");
    debug_assert!(
        dot_sign_rat(&[coeffs[0], coeffs[1], coeffs[2]], m) == Orient::Zero,
        "the strip only exists on a plane parallel to the axis"
    );
    // ★ `_dc` and `_dm` go unused, and that is the derivation showing: the plane's and the
    // direction's own denominators enter only as positive squares in the common denominator, so
    // they cancel. The point's and the origin's do **not** — the `c` term breaks homogeneity in
    // `p`, the same way it does in `point_plane_clearance_rat`.
    let (c, _dc) = lift4(coeffs);
    let (oo, d_o) = lift3(o);
    let (mm, _dm) = lift3(m);
    let (pp, dp) = p.lift();
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let (sn, sd): (BigInt, BigInt) = (rho2.numer().clone(), rho2.denom().clone());
    let n = [c[0].clone(), c[1].clone(), c[2].clone()];
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    debug_assert!(
        (0..3).map(|i| &n[i] * &pp[i]).sum::<BigInt>() + &c[3] * &dp == BigInt::from(0),
        "the point must lie on the plane — the decomposition takes its distance from the axis's, \
         so an off-plane point would be measured against a distance that is not its own"
    );
    // e = n × m, and w = (p − o) scaled by dp·d_o — both integer, both carrying their own
    // positive factor, which is why only the sign of the combination below matters.
    let e: [BigInt; 3] = core::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &n[j] * &mm[k] - &n[k] * &mm[j]
    });
    let w: [BigInt; 3] = core::array::from_fn(|i| &pp[i] * &d_o - &oo[i] * &dp);
    let u = dot(&w, &e);
    // `n·o + c` over the common denominator `dc·d_o`, and the two squared magnitudes.
    let g = dot(&n, &oo) + &c[3] * &d_o;
    let nn = dot(&n, &n);
    let m2 = dot(&mm, &mm);
    StripScale {
        u_sign: orient_of(big_sign(&u)),
        uu: &u * &u * &rd * &sd,
        ww: (&rn * &nn * (&d_o * &d_o) - &g * &g * &rd) * &m2 * (&dp * &dp) * &sd,
        rr: &sn * &nn * &m2 * (&dp * &dp) * (&d_o * &d_o) * &rd,
    }
}

/// **Where a disk of squared radius `rho2` about `p`, lying on the plane, stands relative to the
/// strip.** Exact and total, at any width and any margin.
///
/// The four answers are what a piece with **extent** can say where a point could only say three:
/// clear on the `+` side, clear on the `−` side, reaching the strip without spanning a boundary,
/// or **spanning** one. Consumers fold them differently — the footprint rule treats the last two
/// alike ("did not clear"), while the tangency road needs `Crosses` apart, because a piece that
/// spans the line by itself is exactly the straddle two separate corners would otherwise have to
/// witness between them.
///
/// ```text
/// clear    |U| > W + ρ'        crosses   | |U| − W | < ρ'        else inside
/// ```
///
/// Each is a comparison of a rational against `2√(xy)` for rational `x, y`, so **one sign case
/// and one further squaring** closes it — no radical tower, no new number type. `Crosses` is
/// strict on purpose: a disk touching a ruling at one point spans nothing, which is the same
/// answer a corner sitting on the line gives.
pub fn cylinder_strip_side_margin(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    cylinder_strip_side_extent(
        coeffs,
        &StripReach {
            lo: (p, rho2),
            hi: None,
        },
        o,
        m,
        r2,
    )
}

/// **A piece's reach across the strip, as its two ends** — each end a point on the plane and the
/// **square** of the margin the doors add to it (`ρ²`, the form every radius takes in this family).
///
/// `hi = None` is the symmetric piece a point or a disk is: one end answers both. An arc's ends
/// differ, and stating them is the only way that shape can ask these questions at all.
pub struct StripReach<'a> {
    pub lo: (&'a MeetPoint, &'a BigRat),
    pub hi: Option<(&'a MeetPoint, &'a BigRat)>,
}

/// **Where a piece whose extent across the strip is stated by its two *ends* stands** — the
/// general form of [`cylinder_strip_side_margin`], and the only one an **arc** can use.
///
/// An end is a point on the plane and a squared margin: the low end is `U(p_lo) − ρ_lo'`, the high end
/// `U(p_hi) + ρ_hi'`. A point is both ends with no margin, a disk is one point with the same
/// margin twice — those pass `hi = None` — and an arc's two ends differ, because its angular
/// extent reaches further one way than the other.
///
/// ★★ **`Plus`/`Minus` are exact for every shape; `Crosses` is claimed only where it is proved.**
/// Spanning a boundary is a *positive* fact the tangency road acts on, so under-reporting it is
/// the safe direction — and for an asymmetric extent this door reports [`StripSide::Inside`]
/// rather than prove it. The symmetric case (`hi = None`) keeps the complete answer it always
/// had. An arc that genuinely spans a ruling is therefore "reached the strip, spanned nothing",
/// which every consumer folds as "did not clear"; the day a caller needs the stronger claim from
/// an arc, the proof is `A < W < B` on the two ends and it goes here.
pub fn cylinder_strip_side_extent(
    coeffs: &[Rat; 4],
    reach: &StripReach<'_>,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    use num_bigint::BigInt;
    let (lo, hi) = (reach.lo, reach.hi);
    let side = |o: Orient| match o {
        Orient::Positive => StripSide::Plus,
        Orient::Negative => StripSide::Minus,
        Orient::Zero => StripSide::Inside,
    };
    let s = strip_scale(coeffs, lo.0, lo.1, o, m, r2);
    // The plane clears the cylinder: no strip exists, so nothing can reach it.
    if s.ww.sign() == num_bigint::Sign::Minus {
        return side(s.u_sign);
    }
    let Some(hi) = hi else {
        // Symmetric: one scale answers both boundaries, and `Crosses` with it.
        // `|U| > W + ρ'` — the whole clearance, in the one comparison the scalar family shares.
        if crate::quad::sqrt_exceeds_root_sum(&s.uu, &s.ww, &s.rr) {
            return side(s.u_sign);
        }
        let four = BigInt::from(4);
        let span = &s.uu + &s.ww - &s.rr;
        if span.sign() == num_bigint::Sign::Minus || &span * &span < &four * &s.uu * &s.ww {
            return StripSide::Crosses;
        }
        return StripSide::Inside;
    };
    // Asymmetric: each end answers in its own scale, which is sound because each verdict is a
    // self-contained comparison of that end against `±W`.
    //
    // `U_lo − ρ_lo' > W` is `|U_lo| > W + ρ_lo'` with `U_lo` on the `+` side — the same line,
    // read once per end.
    if s.u_sign == Orient::Positive && crate::quad::sqrt_exceeds_root_sum(&s.uu, &s.ww, &s.rr) {
        return StripSide::Plus;
    }
    let t = strip_scale(coeffs, hi.0, hi.1, o, m, r2);
    if t.ww.sign() == num_bigint::Sign::Minus {
        return side(t.u_sign);
    }
    if t.u_sign == Orient::Negative && crate::quad::sqrt_exceeds_root_sum(&t.uu, &t.ww, &t.rr) {
        return StripSide::Minus;
    }
    StripSide::Inside
}

/// **Does that disk reach *one named* ruling?** — `side` names the ruling in **this** family's
/// vocabulary: the sign of `(p − o)·(n × m)`, the same one [`StripSide::Plus`] is written about.
///
/// ⚠ **The arrangement's `side` is the opposite sign.** `arrangement::ruling_side_signed` measures
/// against `m × n`, so a caller carrying a `MergedRuling` must negate before asking here. Stated
/// rather than absorbed: a predicate that silently accepted either convention would answer about
/// the wrong ruling for whichever caller it was not written for.
///
/// The strip form above answers about **both** boundaries at once, which is what a face-clearance
/// question wants. A caller holding the class's actual edges may carry only one of the two
/// rulings, and asking about the strip would then count a circle that reaches the ruling **not
/// present**. Same derivation, one boundary: `|σW − U| ≤ ρ'`, closed (a touch counts).
pub fn cylinder_ruling_reached(
    coeffs: &[Rat; 4],
    p: &MeetPoint,
    rho2: &BigRat,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
    side: i8,
) -> bool {
    cylinder_ruling_reached_extent(
        coeffs,
        &StripReach {
            lo: (p, rho2),
            hi: None,
        },
        o,
        m,
        r2,
        side,
        true,
    )
}

/// **Does a piece whose reach across the strip is stated by its two *ends* touch the named
/// ruling?** — the general form of [`cylinder_ruling_reached`], and the only one an **arc** can use.
///
/// The piece occupies `[U_lo − ρ_lo', U_hi + ρ_hi']` across the strip; the ruling sits at `σW`. It
/// is touched unless the whole reach lies on one side, which is two comparisons of the same shape
/// — each a signed sum of three roots, and each falling to [`quad::sqrt_root_sum_cmp`] once the
/// signs of `U` and `σ` say which term goes where.
///
/// `hi = None` is the symmetric piece (a point, a disk), which is what this door has always been
/// handed; then one scale answers both ends and the verdict is the one it always gave.
///
/// ★★ **`touch_counts` names the boundary rather than assuming one**. Its predecessor is
/// closed — a piece touching the ruling at one point has *reached* it — and that is the right
/// reading for a clearance. The arrangement's net wants the other one: what it cannot mint is a
/// **crossing**, and an edge tangent to another divides nothing. Two propositions, one door, and
/// the caller says which.
pub fn cylinder_ruling_reached_extent(
    coeffs: &[Rat; 4],
    reach: &StripReach<'_>,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
    side: i8,
    touch_counts: bool,
) -> bool {
    use num_bigint::Sign;
    let (lo, hi) = (reach.lo, reach.hi);
    let s_lo = strip_scale(coeffs, lo.0, lo.1, o, m, r2);
    if s_lo.ww.sign() == Sign::Minus {
        return false; // the plane clears the cylinder — it has no ruling here at all
    }
    let sigma = side.signum();
    // `U_lo − √rr > σ√ww`: the whole reach starts above the ruling. When a touch does **not**
    // count, "starts above" includes starting exactly on it, so the comparison relaxes.
    if strip_end_beyond(&s_lo, s_lo.u_sign, sigma, touch_counts) {
        return false;
    }
    let s_hi = match hi {
        None => s_lo,
        Some(h) => {
            let t = strip_scale(coeffs, h.0, h.1, o, m, r2);
            if t.ww.sign() == Sign::Minus {
                return false;
            }
            t
        }
    };
    // `U_hi + √rr < σ√ww` is the same question with `U` and `σ` both negated.
    !strip_end_beyond(&s_hi, orient_neg(s_hi.u_sign), -sigma, touch_counts)
}

/// Is `su·√uu − √rr > σ·√ww` (or `≥`, when `strict` is false)? — one end of a reach against the
/// named ruling, with the signs of `U` and `σ` deciding which side of the identity each root
/// belongs on, and the boundary named rather than assumed.
fn strip_end_beyond(s: &StripScale, su: Orient, sigma: i8, strict: bool) -> bool {
    use crate::quad::sqrt_root_sum_cmp as cmp;
    use num_bigint::Sign;
    let nil = |x: &num_bigint::BigInt| x.sign() == Sign::NoSign;
    let ord = |a: &num_bigint::BigInt, b: &num_bigint::BigInt| if strict { a > b } else { a >= b };
    match (su, sigma) {
        // `√uu ? √rr + √ww`
        (Orient::Positive, 1) => cmp(&s.uu, &s.rr, &s.ww, strict),
        // `√uu ? √rr`
        (Orient::Positive, 0) => ord(&s.uu, &s.rr),
        // `√uu + √ww ? √rr`, i.e. not `√rr ?̄ √uu + √ww` with the boundary flipped
        (Orient::Positive, _) => !cmp(&s.rr, &s.uu, &s.ww, !strict),
        // `√ww ? √rr`
        (Orient::Zero, -1) => ord(&s.ww, &s.rr),
        // `0 ? √rr`, and `0 ? √rr + √ww` — only an equality can hold, and only when not strict.
        (Orient::Zero, 0) => !strict && nil(&s.rr),
        (Orient::Zero, _) => !strict && nil(&s.rr) && nil(&s.ww),
        // `√ww ? √uu + √rr`
        (Orient::Negative, -1) => cmp(&s.ww, &s.uu, &s.rr, strict),
        // `−√uu ? √rr (+ √ww)` — likewise, only as an equality of zeros.
        (Orient::Negative, 0) => !strict && nil(&s.uu) && nil(&s.rr),
        (Orient::Negative, _) => !strict && nil(&s.uu) && nil(&s.rr) && nil(&s.ww),
    }
}

/// `Orient`'s negation — the sign of `−x` from the sign of `x`.
fn orient_neg(o: Orient) -> Orient {
    match o {
        Orient::Positive => Orient::Negative,
        Orient::Negative => Orient::Positive,
        Orient::Zero => Orient::Zero,
    }
}

/// **Which side of the plane at axis parameter `t` a point stands on.** Exact and total, at any
/// width.
///
/// ★ **This is the second separating axis of one rectangle, not a second rule.** A wall parallel
/// to a cylinder's axis meets that cylinder in a rectangle of the wall's own plane: the strip
/// across ([`cylinder_strip_side`], the first axis) and the lateral face's axis-parameter span
/// along (this one). A face misses the rectangle as soon as it clears *either* axis, because a
/// rectangle is the intersection of the two bands — so the two are read together and neither
/// stands as a rule of its own.
///
/// The parameter is written in the **raw** direction's scale — `axis(t) = o + t·m` with `m`
/// unnormalized, the scale `axis_param_of_plane` produces and a lateral face's span is stored in.
/// Re-scaling to a unit axis here would silently mismatch those spans.
///
/// With `s = (p − o)·m / (m·m)` the point's own parameter and `m·m > 0`,
///
/// ```text
/// sign(s − t) = sign((p − o)·m − t·(m·m))
/// ```
///
/// so no division is formed and no square root ever appears. [`Orient::Zero`] is the point sitting
/// exactly on that plane — which the footprint reading treats as *clear*, the uniform-slab theorem
/// speaking of the **open** slab.
///
/// **Precondition:** `m ≠ 0` (`debug_assert`ed) — a zero direction names no axis. Unlike the strip
/// test there is no on-plane precondition: a point's axis parameter is defined wherever it sits.
pub fn point_axis_side(p: &MeetPoint, o: &[Rat; 3], m: &[Rat; 3], t: Rat) -> Orient {
    use num_bigint::BigInt;
    let (oo, d_o) = lift3(o);
    let (mm, d_m) = lift3(m);
    let (pp, dp) = p.lift();
    let (tn, td) = (BigInt::from(t.numer()), BigInt::from(t.denom()));
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    let m2 = dot(&mm, &mm);
    debug_assert!(
        m2 != BigInt::from(0),
        "a zero direction names no axis, so no parameter along it"
    );
    // `w = p − o` over the common denominator `dp·d_o`, integer throughout.
    let w: [BigInt; 3] = core::array::from_fn(|i| &pp[i] * &d_o - &oo[i] * &dp);
    // sign((p−o)·m − t·|m|²) with the positive common denominator `dp·d_o·d_m²·td` cleared:
    //   (W·M)·d_m·td  vs  tn·|M|²·dp·d_o
    // ★ Every cleared factor is positive — `lift3`/`lift` build denominators from an lcm of
    // `Rat` denominators, and `Rat` keeps its sign in the numerator — so the comparison is the
    // sign it claims to be.
    let lhs = dot(&w, &mm) * &d_m * &td;
    let rhs = &tn * &m2 * &dp * &d_o;
    match big_sign(&(lhs - rhs)) {
        1 => Orient::Positive,
        -1 => Orient::Negative,
        _ => Orient::Zero,
    }
}

/// **A branch point's coordinates as `(A + B√C) / D`** — integer throughout, `D > 0`, `C ≥ 0`.
///
/// The two footprint predicates below are the same questions [`cylinder_strip_side`] and
/// [`point_axis_side`] ask; only the point's *description* differs. A `Vertex::Pierce` has no
/// rational coordinates at all — it is `line.base() + s·line.dir()` with `s` quadratic-irrational
/// — so its coordinates live in `ℚ(√c)`, and every quantity those predicates form from a point is
/// a polynomial in them, hence of the shape `X + Y√C` whose sign the tower already answers.
///
/// ★ **The radicand is integerised here, once.** `√(p/q) = √(p·q)/q`, so `s` is rewritten over a
/// single positive denominator with an **integer** radicand `C = p·q`; every sum downstream then
/// carries one `C` and never has to reconcile two spellings of the same surd.
fn lift_branch(
    line: &quad::MeetLine,
    s: &quad::QuadVal,
) -> (
    [num_bigint::BigInt; 3],
    [num_bigint::BigInt; 3],
    num_bigint::BigInt,
    num_bigint::BigInt,
) {
    use num_bigint::BigInt;
    // `s = (Sa + Sb√C) / Ds`, integer and `Ds > 0`: fold the radical's denominator into the
    // rational part (`√(cp/cq) = √(cp·cq)/cq`), then clear both coefficients' denominators.
    let (cp, cq) = (BigInt::from(s.c().numer()), BigInt::from(s.c().denom()));
    let big_c = &cp * &cq;
    let (san, sad) = (BigInt::from(s.a().numer()), BigInt::from(s.a().denom()));
    let (sbn, sbd) = (BigInt::from(s.b().numer()), BigInt::from(s.b().denom()));
    let s_a = &san * &cq * &sbd;
    let s_b = &sbn * &sad;
    let d_s = &cq * &sad * &sbd;
    // `p = base + s·dir` over `db·Ds·dd`, which is positive because every factor is.
    let (bb, db) = lift3(&line.base());
    let (dd, d_d) = lift3(&line.dir());
    let den = &db * &d_s * &d_d;
    let a: [BigInt; 3] = core::array::from_fn(|i| &bb[i] * &d_s * &d_d + &s_a * &dd[i] * &db);
    let b: [BigInt; 3] = core::array::from_fn(|i| &s_b * &dd[i] * &db);
    (a, b, den, big_c)
}

/// **[`cylinder_strip_side`] asked of a branch point** — the same three-way answer, the same
/// expression, and **total** for the same reason: a sign has no width, so the road to it runs in
/// `BigInt` rather than checked `Rat` (the discipline [`quad::cylinder_radial_side`] states).
///
/// **Preconditions** are that function's, and the on-plane one is the caller's to establish —
/// `quad::plane_side(coeffs, line, s) == Orient::Zero` is the branch spelling of it.
pub fn cylinder_strip_side_branch(
    coeffs: &[Rat; 4],
    line: &quad::MeetLine,
    s: &quad::QuadVal,
    o: &[Rat; 3],
    m: &[Rat; 3],
    r2: &BigRat,
) -> StripSide {
    use num_bigint::BigInt;
    debug_assert!(!r2.is_negative(), "a squared radius is not negative");
    debug_assert!(
        dot_sign_rat(&[coeffs[0], coeffs[1], coeffs[2]], m) == Orient::Zero,
        "the strip only exists on a plane parallel to the axis"
    );
    let (c, _dc) = lift4(coeffs);
    let (oo, d_o) = lift3(o);
    let (mm, _dm) = lift3(m);
    let (pa, pb, dp, big_c) = lift_branch(line, s);
    let (rn, rd): (BigInt, BigInt) = (r2.numer().clone(), r2.denom().clone());
    let n = [c[0].clone(), c[1].clone(), c[2].clone()];
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    // Exactly the rational road's `e` and `w`, with `w` now a **pair**: the point's rational part
    // and its `√C` part travel together through every linear step.
    let e: [BigInt; 3] = core::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &n[j] * &mm[k] - &n[k] * &mm[j]
    });
    let wa: [BigInt; 3] = core::array::from_fn(|i| &pa[i] * &d_o - &oo[i] * &dp);
    let wb: [BigInt; 3] = core::array::from_fn(|i| &pb[i] * &d_o);
    let (ua, ub) = (dot(&wa, &e), dot(&wb, &e));
    let u_sign = quad::sign1_int(&ua, &ub, &big_c);
    if u_sign == Orient::Zero {
        return StripSide::Inside;
    }
    let g = dot(&n, &oo) + &c[3] * &d_o;
    let nn = dot(&n, &n);
    let m2 = dot(&mm, &mm);
    // `U² = (Ua² + Ub²C) + 2·Ua·Ub·√C`, and the right-hand side is rational — so the difference
    // is one `X + Y√C` and the tower reads its sign. `rd` is the squared radius' denominator.
    let rhs = (&rn * &nn * (&d_o * &d_o) - &g * &g * &rd) * &m2 * (&dp * &dp);
    let x = (&ua * &ua + &ub * &ub * &big_c) * &rd - rhs;
    let y = BigInt::from(2) * &ua * &ub * &rd;
    match quad::sign1_int(&x, &y, &big_c) {
        Orient::Positive if u_sign == Orient::Positive => StripSide::Plus,
        Orient::Positive => StripSide::Minus,
        _ => StripSide::Inside,
    }
}

/// **[`point_axis_side`] asked of a branch point** — same expression, same total contract.
pub fn point_axis_side_branch(
    line: &quad::MeetLine,
    s: &quad::QuadVal,
    o: &[Rat; 3],
    m: &[Rat; 3],
    t: Rat,
) -> Orient {
    use num_bigint::BigInt;
    let (oo, d_o) = lift3(o);
    let (mm, d_m) = lift3(m);
    let (pa, pb, dp, big_c) = lift_branch(line, s);
    let (tn, td) = (BigInt::from(t.numer()), BigInt::from(t.denom()));
    let dot = |x: &[BigInt; 3], y: &[BigInt; 3]| -> BigInt {
        (0..3).map(|i| &x[i] * &y[i]).sum::<BigInt>()
    };
    let m2 = dot(&mm, &mm);
    debug_assert!(
        m2 != BigInt::from(0),
        "a zero direction names no axis, so no parameter along it"
    );
    let wa: [BigInt; 3] = core::array::from_fn(|i| &pa[i] * &d_o - &oo[i] * &dp);
    let wb: [BigInt; 3] = core::array::from_fn(|i| &pb[i] * &d_o);
    let x = dot(&wa, &mm) * &d_m * &td - &tn * &m2 * &dp * &d_o;
    let y = dot(&wb, &mm) * &d_m * &td;
    quad::sign1_int(&x, &y, &big_c)
}

/// A rational 4-vector (plane coefficients) as **(integer components, positive denominator)**.
fn lift4(v: &[Rat; 4]) -> ([num_bigint::BigInt; 4], num_bigint::BigInt) {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let den: [BigInt; 4] = core::array::from_fn(|i| BigInt::from(v[i].denom()));
    let d = den.iter().fold(BigInt::from(1), |l, x| l.lcm(x));
    let num = core::array::from_fn(|i| BigInt::from(v[i].numer()) * (&d / &den[i]));
    (num, d)
}
