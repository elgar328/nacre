use super::*;

/// A witness at an `f64`-representable point, stated as the rational it is — the fixture
/// spelling of what `WitnessPoint::exact` used to be. That door is retired: no production
/// caller, and a precondition ("exactly representable") no caller could check — the one
/// production site that handed it a rounded cache named a different point (nacre-ops
/// reuse). `at_nearest` states the same tol `0` here and stays honest elsewhere.
fn exact(c: [f64; 3]) -> Option<WitnessPoint> {
    let b = |x: f64| Rat::try_from_f64(x);
    Some(WitnessPoint::at_nearest([b(c[0])?, b(c[1])?, b(c[2])?]))
}

/// The precision these fixtures judge at. Production chooses it per model
/// ([`judge_precision`]); a fixture pins one so its expectations stay fixed.
const FIXTURE_PREC: usize = 256;

/// The judgement context those fixtures use: that precision, plus the coincidence limit the
/// default derivation gives a unit-scale model (output resolution `2⁻⁵²`, two words further
/// down) and the production cap.
fn fixture() -> Standard {
    Standard {
        prec: FIXTURE_PREC,
        coincidence: Mag::pow2(-180),
        scale: Mag::of(1.0),
        cap: 4096,
    }
}

fn ri(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}

fn deg(n: i128, d: i128) -> Angle {
    Angle::from_deg(ri(n, d)).unwrap()
}

/// Deterministic PRNG (splitmix64) for reproducible stress corpora.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn rng(state: &mut u64, lo: i128, hi: i128) -> i128 {
    lo + (u128::from(splitmix64(state)) % (hi - lo + 1) as u128) as i128
}

fn axis_of(k: i128) -> Axis {
    match k.rem_euclid(3) {
        0 => Axis::X,
        1 => Axis::Y,
        _ => Axis::Z,
    }
}

fn rand_base(st: &mut u64) -> [Rat; 3] {
    [
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
    ]
}

/// The f64 value's distance from the high-precision realization's **midpoint** — what the
/// f64 tol has to bound. (The realization's own radius is a separate, far smaller quantity;
/// `GT` is deep enough that it does not enter these comparisons.)
fn abs_err(f: f64, truth: &HpBounded, gt: usize) -> f64 {
    bf_mag(&BigFloat::from_f64(f, gt).sub(&truth.value, gt, HP_RM).abs())
}

/// Borrow an owned plane-triple as the `&WitnessPoint` tuples the judge takes.
fn tr(t: &[[WitnessPoint; 3]; 3]) -> [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3] {
    [
        (&t[0][0], &t[0][1], &t[0][2]),
        (&t[1][0], &t[1][1], &t[1][2]),
        (&t[2][0], &t[2][1], &t[2][2]),
    ]
}

mod judges;
mod witness_points_and_frames;
