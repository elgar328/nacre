//! **The exact motions a history folds into.**

use super::*;

/// **An affine map whose linear part permutes the axes with signs**: `x ↦ M·x + t`, every entry
/// of `M` in `{0, ±1}` with one nonzero per row and column, `t` rational.
///
/// It is what a chain of translations, quarter turns and axis reflections composes into, and
/// nothing wider: those three are closed under composition in this shape, and every one of them
/// maps rationals to rationals. [`Isometry`] cannot hold the composite — a reflection has
/// `det = −1`, and a turn and a reflection do not commute, so a turn beside a reflection bit
/// does not compose.
///
/// ★ **`M` and `t` are held apart.** `M` is `i8` and cannot overflow; `t` is `None` once a
/// composition overflows `i128`. A direction reads `M` alone, so it still answers there — only
/// points and planes need the offset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AxisAffine {
    m: [[i8; 3]; 3],
    t: Option<[Rat; 3]>,
}

const IDENTITY: [[i8; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

impl AxisAffine {
    /// `x ↦ x + t`.
    pub fn translation(t: [Rat; 3]) -> AxisAffine {
        AxisAffine {
            m: IDENTITY,
            t: Some(t),
        }
    }

    /// The reflection in the coordinate plane `axis = offset` — `x ↦ 2·offset − x` on that axis,
    /// the same map [`mirror_point_rat`] applies.
    pub fn mirror(axis: Axis, offset: Rat) -> AxisAffine {
        let a = axis.index();
        let mut m = IDENTITY;
        m[a][a] = -1;
        let t = offset.checked_mul(Rat::from_int(2)).map(|v| {
            let mut t = [Rat::from_int(0); 3];
            t[a] = v;
            t
        });
        AxisAffine { m, t }
    }

    /// A turn about the line through `pivot` along `axis` — the map [`Isometry::point_rat`]
    /// applies, `p + R·(x − p)` — when the turn is a quarter multiple; `None` otherwise (no other
    /// rational-degree angle has rational `cos`/`sin`, [`Angle::try_exact_cos_sin`]).
    pub fn rotation(r: Rotation) -> Option<AxisAffine> {
        let (cos, sin) = r.angle.try_exact_cos_sin()?;
        let unit = |v: Rat| -> Option<i8> { (v.denom() == 1).then(|| v.numer() as i8) };
        let (c, s) = (unit(cos)?, unit(sin)?);
        let (i, j) = r.axis.plane();
        let mut m = IDENTITY;
        m[i][i] = c;
        m[i][j] = -s;
        m[j][i] = s;
        m[j][j] = c;
        let linear = AxisAffine { m, t: None };
        let t = linear.dir_rat(r.pivot).and_then(|mp| {
            Some([
                r.pivot[0].checked_sub(mp[0])?,
                r.pivot[1].checked_sub(mp[1])?,
                r.pivot[2].checked_sub(mp[2])?,
            ])
        });
        Some(AxisAffine { m, t })
    }

    /// `outer ∘ self` — apply `self`, then `outer`.
    pub fn then(&self, outer: &AxisAffine) -> AxisAffine {
        let mut m = [[0i8; 3]; 3];
        for (r, row) in m.iter_mut().enumerate() {
            for (c, e) in row.iter_mut().enumerate() {
                *e = (0..3).map(|k| outer.m[r][k] * self.m[k][c]).sum();
            }
        }
        let t = (|| {
            let moved = outer.dir_rat(self.t?)?;
            let o = outer.t?;
            Some([
                moved[0].checked_add(o[0])?,
                moved[1].checked_add(o[1])?,
                moved[2].checked_add(o[2])?,
            ])
        })();
        AxisAffine { m, t }
    }

    /// The map applied to a direction: `M·d` (the offset plays no part).
    pub fn dir_rat(&self, d: [Rat; 3]) -> Option<[Rat; 3]> {
        let mut out = [Rat::from_int(0); 3];
        for (o, row) in out.iter_mut().zip(&self.m) {
            for (e, x) in row.iter().zip(d) {
                *o = match e {
                    0 => *o,
                    1 => o.checked_add(x)?,
                    _ => o.checked_sub(x)?,
                };
            }
        }
        Some(out)
    }

    /// The map applied to a point: `M·p + t`.
    pub fn point_rat(&self, p: [Rat; 3]) -> Option<[Rat; 3]> {
        let q = self.dir_rat(p)?;
        let t = self.t?;
        Some([
            q[0].checked_add(t[0])?,
            q[1].checked_add(t[1])?,
            q[2].checked_add(t[2])?,
        ])
    }

    /// The image of the plane `n·x + d = 0`, canonicalized ([`canonical_plane_coeffs`]):
    /// `n' = M·n`, `d' = d − n'·t` (`M` is orthogonal, so `n·x = n'·(M·x)`).
    pub fn plane_coeffs(&self, c: [Rat; 4]) -> Option<[Rat; 4]> {
        let n = self.dir_rat([c[0], c[1], c[2]])?;
        let t = self.t?;
        let mut d = c[3];
        for k in 0..3 {
            d = d.checked_sub(n[k].checked_mul(t[k])?)?;
        }
        canonical_plane_coeffs([n[0], n[1], n[2], d])
    }

    /// The linear part applied to an `f64` direction — exact, since it only permutes and negates.
    pub fn dir_f64(&self, v: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|r| (0..3).map(|c| f64::from(self.m[r][c]) * v[c]).sum())
    }

    /// The inverse map, `x ↦ Mᵀ·x − Mᵀ·t` — `M` permutes the axes with signs, so its inverse is its
    /// transpose, and the inverse is as exact as the map. A direction carried back reads `Mᵀ` alone,
    /// so it answers where the offset overflowed, as [`AxisAffine::dir_rat`] does.
    pub fn inverse(&self) -> AxisAffine {
        let mut m = [[0i8; 3]; 3];
        for (r, row) in m.iter_mut().enumerate() {
            for (c, e) in row.iter_mut().enumerate() {
                *e = self.m[c][r];
            }
        }
        let linear = AxisAffine { m, t: None };
        let t = self.t.and_then(|t| {
            let back = linear.dir_rat(t)?;
            let zero = Rat::from_int(0);
            Some([
                zero.checked_sub(back[0])?,
                zero.checked_sub(back[1])?,
                zero.checked_sub(back[2])?,
            ])
        });
        AxisAffine { m, t }
    }

    /// `det M`: `+1` for a proper motion, `−1` for one carrying an odd number of reflections —
    /// the factor on the cross product of the map's images of points.
    pub fn det(&self) -> i8 {
        let m = &self.m;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
}

#[cfg(test)]
#[path = "tests/axis_affine.rs"]
mod tests;
