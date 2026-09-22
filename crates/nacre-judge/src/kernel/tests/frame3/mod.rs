use super::*;

/// A witness at an `f64`-representable point, stated as the rational it is. Production has no
/// such door: its precondition ("exactly representable") is one no caller can check, and a
/// rounded cache handed to it names a different point. `at_nearest` states the same tol `0`
/// here and stays honest elsewhere.
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
pub(in crate::kernel::frame3) fn fixture() -> Standard {
    Standard {
        prec: FIXTURE_PREC,
        coincidence: Mag::pow2(-180),
        scale: Mag::of(1.0),
        cap: 4096,
    }
}

pub(in crate::kernel::frame3) fn ri(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}

pub(in crate::kernel::frame3) fn deg(n: i128, d: i128) -> Angle {
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

pub(in crate::kernel::frame3) fn rng(state: &mut u64, lo: i128, hi: i128) -> i128 {
    lo + (u128::from(splitmix64(state)) % (hi - lo + 1) as u128) as i128
}

fn axis_of(k: i128) -> Axis {
    match k.rem_euclid(3) {
        0 => Axis::X,
        1 => Axis::Y,
        _ => Axis::Z,
    }
}

pub(in crate::kernel::frame3) fn rand_base(st: &mut u64) -> [Rat; 3] {
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
pub(in crate::kernel::frame3) fn tr(
    t: &[[WitnessPoint; 3]; 3],
) -> [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3] {
    [
        (&t[0][0], &t[0][1], &t[0][2]),
        (&t[1][0], &t[1][1], &t[1][2]),
        (&t[2][0], &t[2][1], &t[2][2]),
    ]
}

fn add3(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
    [
        a[0].checked_add(b[0]).unwrap(),
        a[1].checked_add(b[1]).unwrap(),
        a[2].checked_add(b[2]).unwrap(),
    ]
}

/// Three planes meeting at `v` (each through `v` + two small offsets), all rotated by
/// `(ax, ang, piv)` — an implicit point at `rotate(v)` with heterogeneous provenance.
pub(in crate::kernel::frame3) fn triple_pts(
    v: [Rat; 3],
    st: &mut u64,
    ax: Axis,
    ang: Angle,
    piv: [Rat; 3],
) -> [[WitnessPoint; 3]; 3] {
    let plane = |st: &mut u64| {
        let off = |st: &mut u64| {
            [
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
            ]
        };
        let (o1, o2) = (off(st), off(st));
        [
            WitnessPoint::at(v).rotate_about(ax, ang, piv),
            WitnessPoint::at(add3(v, o1)).rotate_about(ax, ang, piv),
            WitnessPoint::at(add3(v, o2)).rotate_about(ax, ang, piv),
        ]
    };
    [plane(st), plane(st), plane(st)]
}

/// A random rotated three-plane triple (own random center, axis, inexact angle,
/// pivot) — heterogeneous provenance. Sequences the RNG draws so each `&mut st`
/// borrow ends before the next.
pub(in crate::kernel::frame3) fn rand_triple(st: &mut u64) -> [[WitnessPoint; 3]; 3] {
    let v = rand_base(st);
    let ax = axis_of(rng(st, 0, 2));
    let angle = deg(rng(st, 0, 360_000), rng(st, 1, 9973));
    let piv = rand_base(st);
    triple_pts(v, st, ax, angle, piv)
}

/// Relative error of `got` against `want`, as an `f64` magnitude — computed **in astro-float**
/// so the comparison is not limited by `bf_mag`'s power-of-two rounding. Used by the
/// normalization tests, where a factor-of-two slop would hide a real units error.
pub(in crate::kernel::frame3) fn rel_err(got: &BigFloat, want: &BigFloat, prec: usize) -> f64 {
    let d = got.sub(want, prec, HP_RM);
    bf_mag(&d.abs()) / bf_mag(&want.abs())
}

mod judges;
pub(in crate::kernel::frame3) mod truth;
mod witness_points_and_frames;
