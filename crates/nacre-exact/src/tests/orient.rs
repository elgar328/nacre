//! Unit tests of the orientation sign past `i128` — the wide arm (`residual_sign_big`) against
//! the narrow one.

use super::*;
use proptest::prelude::*;

/// ★★ **The population the narrow route cannot reach — and the reason the sign is total.**
///
/// Two decimal-window denominators of `10^22` put the determinant's product denominator at
/// `10^44`, past `i128` — precisely what `from_decimal`'d profile coordinates produce. The
/// answers are hand-computable: `(0,0) → (1/d,0) → (0,1/d)` turns counter-clockwise, its
/// mirror clockwise, and a doubled point on the same ray is collinear.
#[test]
fn an_orientation_past_i128_still_gets_its_sign() {
    let d = 10i128.pow(22);
    let r = |n: i128, den: i128| Rat::new(n, den).unwrap();
    let (o, x, y) = (
        [r(0, 1), r(0, 1)],
        [r(1, d), r(0, 1)],
        [r(0, 1), r(1, d + 1)],
    );
    // The narrow route must actually be dead here, or this pins nothing (self-qualification).
    assert!(
        r(1, d).checked_mul(r(1, d + 1)).is_none(),
        "the product denominator was expected to overflow i128"
    );
    assert_eq!(orient2d_rat(o, x, y), 1, "counter-clockwise");
    assert_eq!(orient2d_rat(o, y, x), -1, "clockwise");
    let far = [r(2, d), r(2, d)];
    let near = [r(1, d + 1), r(1, d + 1)];
    assert_eq!(orient2d_rat(o, near, far), 0, "one ray, three points");
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// **The residual's two arms must be the same function** — `plane_residual_sign` runs the
    /// `Rat` substitution first, so wherever it answers the `BigInt` arm is never consulted
    /// (the `orient2d` arrangement, one dimension up).
    #[test]
    fn the_big_residual_answers_what_the_narrow_one_does(
        // Small enough that the canonical name always fits i128 (the lifted integers stay
        // near 2^40 and the derivation's peak near 2^122), so `narrow()` below cannot shrug.
        xs in prop::array::uniform9(-(1i64 << 10)..(1i64 << 10)),
        ds in prop::array::uniform9(1i64..(1i64 << 10)),
    ) {
        let r = |i: usize| Rat::new(xs[i] as i128, ds[i] as i128).unwrap();
        let (a, b, c) = ([r(0), r(1), r(2)], [r(3), r(4), r(5)], [r(6), r(7), r(8)]);
        if let Some(name) = plane_name_exact(a, b, c) {
            let coeffs = name.narrow().expect("small operands stay narrow");
            let ci = coeffs.map(|x| num_bigint::BigInt::from(x.numer()));
            // The naming points themselves: both arms must call them on-plane...
            for p in [a, b, c] {
                prop_assert_eq!(plane_residual_sign(&name, p), 0);
                prop_assert_eq!(residual_sign_big(&ci, p), 0);
            }
            // ...and an off-plane probe (a naming point pushed along the normal) must get
            // the same nonzero sign from both.
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            if let Some(q) = (|| -> Option<[Rat; 3]> {
                Some([
                    a[0].checked_add(n[0])?,
                    a[1].checked_add(n[1])?,
                    a[2].checked_add(n[2])?,
                ])
            })() {
                let s = plane_residual_sign(&name, q);
                prop_assert_eq!(s, 1, "n·n > 0: the push is to the positive side");
                prop_assert_eq!(residual_sign_big(&ci, q), s);
            }
        }
    }
}
