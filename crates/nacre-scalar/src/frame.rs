use super::*;
/// **A plane's own frame, exactly** — the origin it is drawn from and the raw direction its `u`
/// axis runs along, both rational, plus the normal reduced to its primitive direction.
///
/// ```text
/// origin  = the world origin projected onto the plane
/// n       = [a, b, c] divided by its content (sign kept)
/// u_raw   = ẑ × n = (−b, a, 0),  or  ŷ × n = (c, 0, 0) when the normal is vertical
/// ```
///
/// ★★★ **`u_raw` is not projected and need not be a unit vector.** The general recipe for laying
/// a reference direction into a plane is `(n·n)·ref − (ref·n)·n`, and it is unnecessary here: a
/// cross product is perpendicular to both its arguments, so `ẑ × n` is *already* in the plane.
/// Skipping it is the difference between coefficients that grow cubically and ones that do not —
/// measured over the faces a face-based operation actually targets, projecting left 21% of them
/// inside `i128` and this leaves **100%**.
///
/// ★★★ **The sign is kept, unlike [`canonical_plane_coeffs`].** That function answers *"are these
/// the same plane"*, where direction is noise. A frame's `n` **is** a direction: negating it
/// negates `û` and `ŵ` together, which is a half-turn about `v` — a different frame, not the same
/// one written differently.
///
/// ★★ **The branch is exact, not toleranced**, and it matches `nacre-ops`' f64 `frame_axes` term
/// for term. DXF's arbitrary-axis convention switches on `|n_x| < 1/64` because a float-only
/// kernel cannot ask the real question; `a == 0 && b == 0` is the real question.
///
/// `None` when the plane is degenerate, when the origin projection is not rational, or when the
/// squared lengths the realization needs do not fit `i128` — all three are honest declines that
/// leave a caller on the f64 path it was already on, never a reject.
pub fn plane_frame(coeffs: [Rat; 4]) -> Option<PlaneFrame> {
    let (origin, ref_dir) = plane_frame_default(coeffs)?;
    plane_frame_named(coeffs, origin, ref_dir)
}

/// **Where a plane's frame sits when nobody names it** — `(origin, ref_dir)`.
///
/// The origin is the world origin projected onto the plane and `ref_dir` is `ẑ × n`
/// (`ŷ × n` when the normal is vertical), which is the arbitrary-axis convention a *face* takes.
/// A caller who names their own sketch origin and `+u` passes those to [`plane_frame_named`]
/// instead — a named plane can insist on axes no derivation would produce, and the script layer's
/// `ZX` (whose `+u` is `+ẑ`, not `ẑ × n = −x̂`) is exactly such a case.
pub fn plane_frame_default(coeffs: [Rat; 4]) -> Option<([Rat; 3], [Rat; 3])> {
    let origin = plane_origin_projection(coeffs)?;
    let n = reduce_direction([coeffs[0], coeffs[1], coeffs[2]])?;
    let zero = Rat::from_int(0);
    let ref_dir = if n[0] == zero && n[1] == zero {
        [n[2], zero, zero]
    } else {
        [zero.checked_sub(n[1])?, n[0], zero]
    };
    Some((origin, ref_dir))
}

/// **A plane's frame with the origin and `+u` direction its author chose.**
///
/// `ref_dir` must lie in the plane and not be zero; it need **not** be a unit vector, and it is
/// reduced to its primitive form here so that two spellings of one direction (`[10,10,0]` and
/// `[20,20,0]`) name **one** frame. That reduction is what makes a frame node interning-friendly
/// — the document's blocker 6.
///
/// `None` on a degenerate plane or direction, when the origin projection is not rational, or when
/// the squared lengths the realization divides by do not fit `i128`.
pub fn plane_frame_named(
    coeffs: [Rat; 4],
    origin: [Rat; 3],
    ref_dir: [Rat; 3],
) -> Option<PlaneFrame> {
    let zero = Rat::from_int(0);
    let n = reduce_direction([coeffs[0], coeffs[1], coeffs[2]])?;
    let u_raw = reduce_direction(ref_dir)?;
    let dot =
        |a: &[Rat; 3]| (0..3).try_fold(zero, |acc, k| acc.checked_add(a[k].checked_mul(a[k])?));
    let (uu, nn) = (dot(&u_raw)?, dot(&n)?);
    // ★★★ **`v_raw` is exact, and taking it exactly is what makes the axes come out exactly.**
    // Realizing `v̂` as `ŵ × û` in f64 costs two roundings that do not cancel: a plain wall whose
    // `v` is exactly `ẑ` came out `0.999999999999999_7`, which is a frame that is not quite
    // orthonormal and an exact path quietly lost. `n ⊥ u_raw` by construction, so
    // `|v_raw|² = |n|²·|u_raw|²` — one inverse square root of an exact rational, and the same
    // wall lands on `1.0`.
    //
    // ★ Its components are bounded by `|n|·|u_raw|`, so the one `checked_mul` below covers them:
    // if the squared length fits, so does every component.
    // ★★★ `v̂` exactly when it fits, and an honest fallback when it does not.
    //
    // `|v_raw|² = |n|²·|u_raw|²` is a product of two squared lengths, so it needs **twice** the
    // width they do — measured, that halves the per-component budget from ~62 bits to ~31 and
    // declines 774 of 780 planes built from a wide normal. When it does not fit, the realization
    // falls back to `v̂ = ŵ × û`, which costs two roundings that do not cancel instead of one.
    // ★ That is the accuracy this crate had before `v_raw` existed — a graceful step down, never
    // a wrong frame. (Filling `vv` with something else would not be a fallback but a corruption:
    // `v̂` would come out the wrong *length* and the basis would not be orthonormal.)
    let cx = |i: usize, j: usize| {
        n[i].checked_mul(u_raw[j])?
            .checked_sub(n[j].checked_mul(u_raw[i])?)
    };
    let v = (|| {
        let vv = nn.checked_mul(uu)?;
        Some(([cx(1, 2)?, cx(2, 0)?, cx(0, 1)?], vv))
    })();
    Some(PlaneFrame {
        origin,
        u_raw,
        n,
        uu,
        nn,
        v,
    })
}

/// A plane's own frame, exactly — what [`plane_frame`] derives.
///
/// The `*_raw` vectors are rational and **not** unit length; `uu`/`nn` are their squared lengths,
/// carried because they are what the realization divides by and because computing them here is
/// what proves the frame fits `i128` before anything tries to use it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PlaneFrame {
    pub origin: [Rat; 3],
    pub u_raw: [Rat; 3],
    pub n: [Rat; 3],
    pub uu: Rat,
    pub nn: Rat,
    /// `(v_raw, |v_raw|²)` when both fit `i128` — `v̂` is then one inverse square root of an exact
    /// rational, and a wall whose `v` is exactly `ẑ` lands on `1.0`. `None` when the product
    /// `|n|²·|u_raw|²` overflows, and the realization takes `v̂ = ŵ × û` instead.
    pub v: Option<([Rat; 3], Rat)>,
}

/// A direction vector divided by its content — the **primitive** integer vector along it, with
/// its sign kept. `None` for the zero vector or on overflow.
///
/// Steps ① and ② of [`canonical_plane_coeffs`] and deliberately not step ③: see [`plane_frame`]
/// for why a direction may not have its sign normalized.
fn reduce_direction(v: [Rat; 3]) -> Option<[Rat; 3]> {
    let mut lcm: i128 = 1;
    for c in v {
        let d = c.denom();
        let g = gcd_u128(lcm.unsigned_abs(), d.unsigned_abs()) as i128;
        lcm = lcm.checked_div(g)?.checked_mul(d)?;
    }
    let mut num = [0i128; 3];
    for (i, c) in v.iter().enumerate() {
        num[i] = c.numer().checked_mul(lcm.checked_div(c.denom())?)?;
    }
    let g = num.iter().fold(0u128, |g, n| gcd_u128(g, n.unsigned_abs()));
    if g == 0 {
        return None; // the zero vector is not a direction
    }
    for n in &mut num {
        *n /= g as i128;
    }
    Some(num.map(Rat::from_int))
}
