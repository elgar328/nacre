//! The indirect orient3d's high-precision realization.

use super::*;
use crate::kernel::frame3::tests::rel_err;

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
    let big = |v: i128| HpBounded::exact(BigFloat::from_i128(v, prec));
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
                        HpBounded::new(
                            BigFloat::from_i128(d.0 * k, prec).div(
                                &BigFloat::from_i128(d.1, prec),
                                prec,
                                HP_RM,
                            ),
                            Mag::ZERO,
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
                let got = m
                    .value
                    .div(&d.value.mul(&norm, prec, HP_RM), prec, HP_RM)
                    .abs();
                let want =
                    BigFloat::from_i128(hn, prec).div(&BigFloat::from_i128(hd, prec), prec, HP_RM);
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
