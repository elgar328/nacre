use super::*;

/// **An interval that may hold zero proves no sign.** The walk's zero arms (a corner on the ray,
/// the probe on a step) stand only on exact values; read from a box, a value that straddles zero
/// must abstain — folded into either side, a corner truly on the ray would be counted by both of
/// its steps or by neither, and the parity would be wrong without a refusal.
#[test]
fn an_interval_proves_a_sign_only_apart_from_zero() {
    use nacre_exact::{HpBounded, Orient, Rat};
    let prec = 128;
    let third = HpBounded::of_rat(Rat::new(1, 3).expect("a third"), prec);
    let minus_third = HpBounded::of_rat(Rat::new(-1, 3).expect("a third"), prec);
    assert_eq!(interval_orient(&third), Some(Orient::Positive));
    assert_eq!(interval_orient(&minus_third), Some(Orient::Negative));
    // `1/3 − 1/3` at a rounding precision: a radius about zero, neither sign proved.
    let straddle = third.sub(
        &HpBounded::of_rat(Rat::new(1, 3).expect("a third"), prec),
        prec,
    );
    assert_eq!(interval_orient(&straddle), None);
}
