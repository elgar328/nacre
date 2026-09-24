use super::*;

#[test]
fn orient3d_sign_convention_is_det_a_minus_d() {
    // Reproduces geometry-predicates' own golden (det[a−d, b−d, c−d] = 10).
    assert_eq!(
        orient3d(
            [0.0, 1.0, 6.0],
            [2.0, 3.0, 4.0],
            [4.0, 5.0, 1.0],
            [6.0, 2.0, 5.3]
        ),
        10.0
    );
    // Unit tetra: a,b,c the CCW xy triangle, d one unit up (+z). det gives −1,
    // pinning the convention (positive ⇒ d on the −(b−a)×(c−a) side).
    assert_eq!(
        orient3d(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ),
        -1.0
    );
    // …and the opposite side flips the sign.
    assert_eq!(
        orient3d(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0]
        ),
        1.0
    );
}

#[test]
fn orient3d_coplanar_is_exactly_zero() {
    // Four points in the z = 0 plane.
    assert_eq!(
        orient3d(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0]
        ),
        0.0
    );
}

// Smoke test: the adaptive-arithmetic primitives the indirect predicates
// (M5-a) will be built from are exposed and usable. `two_product` / `two_sum`
// return `[lo, hi]` (tail first) with `hi` the rounded result and `hi + lo`
// the exact value — the roundoff tail is what makes exact expansions possible.
#[test]
fn adaptive_arithmetic_primitives_are_usable() {
    use geometry_predicates::predicates::{two_diff, two_product, two_sum};

    // Exact cases: small integers lose nothing, so the tail is zero.
    assert_eq!(two_product(3.0, 5.0), [0.0, 15.0]);
    assert_eq!(two_sum(1.0, 2.0), [0.0, 3.0]);
    assert_eq!(two_diff(5.0, 3.0), [0.0, 2.0]);

    // Roundoff cases: `hi` is the rounded result, and the tail recovers
    // exactly what an f64 result would drop (this is what enables exact
    // expansions). (1 + 2⁻³⁰)² = 1 + 2⁻²⁹ + 2⁻⁶⁰; the 2⁻⁶⁰ term is below the
    // ULP at 1, so f64 loses it but the tail keeps it.
    let a = 1.0 + 2f64.powi(-30);
    let [lo, hi] = two_product(a, a);
    assert_eq!(hi, a * a); // hi == fl(a·a)
    assert_ne!(lo, 0.0); // the lost bits are recovered, not dropped

    // 2⁵³ + 0.5 sits exactly halfway; round-to-even gives 2⁵³, tail 0.5.
    assert_eq!(two_sum(2f64.powi(53), 0.5), [0.5, 2f64.powi(53)]);
}

#[test]
fn det3_golden() {
    // Identity: det = 1 > 0.
    assert_eq!(
        det3_sign([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        1
    );
    // One row swap flips the sign.
    assert_eq!(
        det3_sign([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        -1
    );
}

#[test]
fn det3_singular_is_zero() {
    // Two equal rows.
    assert_eq!(
        det3_sign([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        0
    );
    // Row 3 = row 1 + row 2 (linearly dependent).
    assert_eq!(
        det3_sign([[1.5, -2.0, 4.0], [3.0, 7.0, -1.0], [4.5, 5.0, 3.0]]),
        0
    );
}

/// The exact integer determinant sign in i128 — an independent oracle for
/// integer inputs.
fn det3_i128(m: [[i128; 3]; 3]) -> i128 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

use proptest::prelude::*;

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The two machines agree on the sign of a determinant: `det3`'s exact
    /// expansion and `det3_sign`'s adaptive `orient3d`. Since `det3_sign` is now
    /// *defined* as `sgn(orient3d(..))`, comparing it to `orient3d` directly would
    /// be a tautology — so this compares it to the **expansion**, which `cramer`
    /// still depends on. Two paths, one sign.
    #[test]
    fn prop_det3_sign_matches_the_expansion(
        m in prop::array::uniform3(prop::array::uniform3(-1e6f64..1e6)),
    ) {
        prop_assert_eq!(det3_sign(m), det3(m).sign());
    }

    /// A rank-deficient matrix has determinant exactly zero, and both machines must
    /// read it. Rows are integers and the combining coefficients are integers, so
    /// `row2 = a·row0 + b·row1` is an exact integer vector — the matrix is singular in
    /// f64, not merely near it. This is the case the oracle above cannot reach on its
    /// own (random continuous rows never land exactly singular).
    #[test]
    fn prop_a_singular_matrix_reads_zero(
        r0 in prop::array::uniform3(-1000i32..1000),
        r1 in prop::array::uniform3(-1000i32..1000),
        a in prop::sample::select(vec![-2i32, -1, 1, 2]),
        b in prop::sample::select(vec![-2i32, -1, 1, 2]),
    ) {
        let row = |r: [i32; 3]| r.map(f64::from);
        let r2 = std::array::from_fn(|i| f64::from(a * r0[i] + b * r1[i]));
        let m = [row(r0), row(r1), r2];
        prop_assert_eq!(det3_sign(m), 0);
        prop_assert_eq!(sign_f64(orient3d(m[0], m[1], m[2], [0.0; 3])), 0);
    }

    /// Independent integer oracle: the exact i128 determinant sign.
    #[test]
    fn prop_det3_sign_matches_i128(
        e in prop::array::uniform3(prop::array::uniform3(-1_000_000i64..1_000_000)),
    ) {
        let mf = e.map(|row| row.map(|v| v as f64));
        let mi = e.map(|row| row.map(|v| v as i128));
        prop_assert_eq!(det3_sign(mf), det3_i128(mi).signum() as i8);
    }
}

// ---- indirect_plane_side ----

/// ★ **The coefficient form and the triangle form are the same question** — checked on a
/// fourth plane whose two descriptions genuinely agree, which is the only case where asking
/// both is meaningful. The triangle's right-hand normal is what `indirect_orient3d` measures
/// against, so the two agree when the triangle is wound to the coefficients' normal and
/// oppose when it is not. That relation is what a face's `frame_sign` records.
#[test]
fn the_coefficient_form_answers_what_the_triangle_form_does() {
    // x = 1, y = 2, z = 3 meet at (1, 2, 3).
    let tp = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -2.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    // Fourth plane z = 0, normal +z. Its triangle, wound so the RH normal is +z.
    let c = [0.0, 0.0, 1.0, 0.0];
    let (t0, t1, t2) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert_eq!(
        indirect_plane_side(&tp, c),
        indirect_orient3d(&tp, t0, t1, t2),
        "same winding ⇒ same sign"
    );
    assert_eq!(indirect_plane_side(&tp, c), 1, "(1,2,3) is above z = 0");
    assert_eq!(
        indirect_plane_side(&tp, c),
        -indirect_orient3d(&tp, t0, t2, t1),
        "reversed winding ⇒ opposite sign"
    );
    // Scaling the coefficients cannot move the point; negating them flips the side.
    assert_eq!(indirect_plane_side(&tp, [0.0, 0.0, 7.0, 0.0]), 1);
    assert_eq!(indirect_plane_side(&tp, [0.0, 0.0, -1.0, 0.0]), -1);
}

#[test]
fn a_point_on_the_plane_reads_zero() {
    let tp = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -2.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    // z = 3 passes through (1, 2, 3) exactly.
    assert_eq!(indirect_plane_side(&tp, [0.0, 0.0, 1.0, -3.0]), 0);
    // And so does a slanted plane through it: 2x − y + z − 3 = 0.
    assert_eq!(indirect_plane_side(&tp, [2.0, -1.0, 1.0, -3.0]), 0);
}

/// The implicit point is a *ratio*, so a triple whose determinant is negative must not flip
/// the answer — that is what the `sign(D)` factor is for.
#[test]
fn a_negatively_oriented_triple_names_the_same_point() {
    let a = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -2.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    // The same three planes, two swapped: D changes sign, the point does not.
    let b = ThreePlane([
        [0.0, 1.0, 0.0, -2.0],
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    for c in [
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0, 0.0],
        [2.0, -1.0, 1.0, -3.0],
        [1.0, 1.0, 1.0, -10.0],
    ] {
        assert_eq!(indirect_plane_side(&a, c), indirect_plane_side(&b, c));
    }
}

// ---- planes_coplanar ----

#[test]
fn planes_coplanar_names_the_same_plane() {
    // Same plane, and the same plane scaled by a negative (opposite normal).
    assert!(planes_coplanar(
        [0.0, 0.0, 1.0, -1.0],
        [0.0, 0.0, 1.0, -1.0]
    ));
    assert!(planes_coplanar(
        [0.0, 0.0, 1.0, -1.0],
        [0.0, 0.0, -2.0, 2.0]
    ));
    // Parallel but offset (z=1 vs z=2): not the same plane.
    assert!(!planes_coplanar(
        [0.0, 0.0, 1.0, -1.0],
        [0.0, 0.0, 1.0, -2.0]
    ));
    // Non-parallel normals.
    assert!(!planes_coplanar(
        [0.0, 0.0, 1.0, -1.0],
        [0.0, 1.0, 0.0, -1.0]
    ));
    // Through the origin (d = 0): coplanarity reduces to parallel normals, but a
    // parallel plane with d ≠ 0 is still distinct.
    assert!(planes_coplanar([1.0, 1.0, 0.0, 0.0], [2.0, 2.0, 0.0, 0.0]));
    assert!(!planes_coplanar(
        [1.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, -1.0]
    ));
}

/// Why `coplanar` is not an absolute `1e-9` tolerance: two
/// planes exactly `1e-9` apart (`z = 0` and `z = 1e-9`). An absolute
/// `distance ≤ 1e-9` test **false-merges** them; the exact rank-1 test splits them
/// (`minor(2,3) = 1e9·(−1) − 0·1e9 = −1e9 ≠ 0`). `two_product(1e9, 1)` is exact —
/// no overflow, no rounding.
#[test]
fn planes_coplanar_splits_a_1e_9_gap_the_absolute_tolerance_would_merge() {
    assert!(!planes_coplanar(
        [0.0, 0.0, 1e9, 0.0],
        [0.0, 0.0, 1e9, -1.0]
    ));
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Scale-invariance — the property an absolute-length coincidence tolerance
    /// lacks. Scaling either plane's coefficients by any nonzero λ names the same
    /// plane, so the decision is unchanged; λ is a power of two so the scaled
    /// coefficients are exact and the invariance is exact. A plane is always
    /// coplanar with its own scaling.
    #[test]
    fn prop_planes_coplanar_is_scale_invariant(
        a in prop::array::uniform4(-50.0f64..50.0),
        b in prop::array::uniform4(-50.0f64..50.0),
        lambda in prop::sample::select(vec![-4.0f64, -2.0, -0.5, 0.5, 2.0, 4.0]),
    ) {
        let scale = |p: [f64; 4], k: f64| p.map(|c| c * k);
        prop_assert_eq!(planes_coplanar(a, b), planes_coplanar(scale(a, lambda), b));
        prop_assert_eq!(planes_coplanar(a, b), planes_coplanar(a, scale(b, lambda)));
        prop_assert!(planes_coplanar(a, scale(a, lambda)));
    }
}

// ---- indirect orient3d (M5-a2) ----

/// The sign of an f64, `+1`/`-1`/`0`. Not `f64::signum`, which maps `0.0` to
/// `+1.0` — a coplanar `orient3d` (exactly `0.0`) must read as `0`.
fn sign_f64(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// The integer value of an expansion, exact when every component is an
/// integer within `i128` (each component is `< 2⁵³`, so `as i128` is lossless).
fn expansion_to_i128(e: &Expansion) -> i128 {
    e.0.iter().map(|&c| c as i128).sum()
}

#[test]
fn indirect_cmp_coord_orders_two_axis_points() {
    // (1,2,3) and (1,5,0), each cut out by three axis planes.
    let at = |p: [f64; 3]| {
        ThreePlane([
            [1.0, 0.0, 0.0, -p[0]],
            [0.0, 1.0, 0.0, -p[1]],
            [0.0, 0.0, 1.0, -p[2]],
        ])
    };
    let a = at([1.0, 2.0, 3.0]);
    let b = at([1.0, 5.0, 0.0]);
    assert_eq!(indirect_cmp_coord(&a, &b, 0), 0); // equal x
    assert_eq!(indirect_cmp_coord(&a, &b, 1), -1); // 2 < 5
    assert_eq!(indirect_cmp_coord(&a, &b, 2), 1); // 3 > 0
    assert_eq!(indirect_cmp_coord(&b, &a, 1), 1); // antisymmetric

    // A genuinely irrational-looking point: the planes are not axis-aligned, and the
    // ratio Dx/D has no exact f64 form. Only the sign is asked for.
    let tilted = ThreePlane([
        [3.0, 1.0, 0.0, -1.0],
        [0.0, 7.0, 1.0, -1.0],
        [1.0, 0.0, 5.0, -1.0],
    ]);
    assert_eq!(indirect_cmp_coord(&tilted, &tilted, 0), 0);
}

/// Two *different* points whose denominators have opposite signs. Cross-multiplying
/// `Na·Db − Nb·Da` flips with `D`, so the numerator alone reports the order backwards
/// half the time; `sign(Da)·sign(Db)` puts it right. Drop that factor and this test
/// says so — the same-point test cannot, since there the numerator is exactly zero.
#[test]
fn opposite_denominators_keep_the_order() {
    // (1,2,3): axis planes in order ⇒ D = +1.
    let a = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -2.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    // (4,2,3): the same planes with two rows swapped ⇒ D = −1.
    let b = ThreePlane([
        [0.0, 1.0, 0.0, -2.0],
        [1.0, 0.0, 0.0, -4.0],
        [0.0, 0.0, 1.0, -3.0],
    ]);
    assert_eq!(indirect_cmp_coord(&a, &b, 0), -1); // 1 < 4
    assert_eq!(indirect_cmp_coord(&b, &a, 0), 1);
    assert_eq!(indirect_cmp_coord(&a, &b, 1), 0); // 2 == 2
    assert_eq!(indirect_cmp_coord(&a, &b, 2), 0); // 3 == 3
}

/// The plane `[n₀, n₁, n₂, −n·p]` through integer point `p` with normal `n`.
fn plane_through(n: [i64; 3], p: [i64; 3]) -> [f64; 4] {
    let dot = n[0] * p[0] + n[1] * p[1] + n[2] * p[2];
    [n[0] as f64, n[1] as f64, n[2] as f64, -dot as f64]
}

/// The filter decides the everyday case, and the exact path is what a coplanar
/// input costs. Both halves matter: a filter that never fires buys nothing, and
/// one that fires on a zero would be wrong.
#[test]
fn the_filter_answers_a_clear_sign_and_declines_a_coplanar_one() {
    let planes = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -1.0],
        [0.0, 0.0, 1.0, -1.0],
    ]);
    // (1,1,1) against the plane x+y+z=3 it lies on: the bound cannot separate 0.
    assert_eq!(
        indirect_orient3d_filter(&planes, [3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]),
        None
    );
    // The same point against z=0: a clear sign, and no expansion is built.
    assert_eq!(
        indirect_orient3d_filter(&planes, [0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        Some(indirect_orient3d_exact(
            &planes,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ))
    );
}

#[test]
fn indirect_orient3d_coplanar_is_zero() {
    // Planes x=1, y=1, z=1 ⇒ implicit point (1,1,1).
    let planes = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -1.0],
        [0.0, 0.0, 1.0, -1.0],
    ]);
    // q,r,s span the plane x+y+z=3, which contains (1,1,1) ⇒ coplanar ⇒ 0.
    assert_eq!(
        indirect_orient3d(&planes, [3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]),
        0
    );
}

#[test]
fn indirect_orient3d_known_sign() {
    let planes = ThreePlane([
        [1.0, 0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0, -1.0],
        [0.0, 0.0, 1.0, -1.0],
    ]);
    // p=(1,1,1) above the CCW triangle in z=0 ⇒ det[p−s,q−s,r−s] = +1 (hand-computed).
    assert_eq!(
        indirect_orient3d(&planes, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        1
    );
    // Swapping two explicit points flips the sign.
    assert_eq!(
        indirect_orient3d(&planes, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        -1
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The filter's one obligation: when it answers, it answers correctly. It is a
    /// rounding-error bound, so a wrong sign here is not a slow path but a silently
    /// wrong b-rep. The scale factor spans nine decades so the bound is tested where
    /// it is tight, not only where it is slack.
    #[test]
    fn prop_the_filter_never_disagrees_with_the_exact_path(
        planes in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        q in prop::array::uniform3(-1e3f64..1e3),
        r in prop::array::uniform3(-1e3f64..1e3),
        s in prop::array::uniform3(-1e3f64..1e3),
        scale in -4i32..5,
    ) {
        let k = 10f64.powi(scale);
        let tp = ThreePlane(planes.map(|pl| pl.map(|c| c * k)));
        let (q, r, s) = (q.map(|c| c * k), r.map(|c| c * k), s.map(|c| c * k));
        if let Some(fast) = indirect_orient3d_filter(&tp, q, r, s) {
            prop_assert_eq!(fast, indirect_orient3d_exact(&tp, q, r, s));
            prop_assert_ne!(fast, 0); // the filter may never claim a zero
        }
    }

    /// The near-degenerate triple: the third plane is a *rounded* linear combination
    /// of the other two, so `D` is noise rather than an exact zero. The bound should
    /// swallow that noise and decline — but "should" is not "must", so this asserts
    /// only what is owed: whatever the filter says, the expansions say too.
    #[test]
    fn prop_the_filter_agrees_near_a_degenerate_triple(
        a in prop::array::uniform4(-100f64..100.0),
        b in prop::array::uniform4(-100f64..100.0),
        t in -3f64..3.0,
        q in prop::array::uniform3(-100f64..100.0),
        r in prop::array::uniform3(-100f64..100.0),
        s in prop::array::uniform3(-100f64..100.0),
    ) {
        let c: [f64; 4] = std::array::from_fn(|i| a[i] + t * b[i]);
        let tp = ThreePlane([a, b, c]);
        if let Some(fast) = indirect_orient3d_filter(&tp, q, r, s) {
            prop_assert_eq!(fast, indirect_orient3d_exact(&tp, q, r, s));
        }
    }

    /// Primary oracle: the exact i128 result. `M` is computed by a **generic**
    /// `det3_i128([Row1, q−s, r−s])`, a different path than the implementation's
    /// hand-factored `Row1·cross`, so a factoring/sign/cross bug shows up as a
    /// mismatch (no self-consistency trap). A `p`-verification pins the oracle's
    /// own Cramer to truth. Inputs in [−100,100] keep every i128 term ≪ 1.7×10³⁸.
    #[test]
    fn prop_indirect_orient3d_matches_i128(
        planes in prop::array::uniform3(prop::array::uniform4(-100i64..=100)),
        q in prop::array::uniform3(-100i64..=100),
        r in prop::array::uniform3(-100i64..=100),
        s in prop::array::uniform3(-100i64..=100),
    ) {
        let n = [
            [planes[0][0] as i128, planes[0][1] as i128, planes[0][2] as i128],
            [planes[1][0] as i128, planes[1][1] as i128, planes[1][2] as i128],
            [planes[2][0] as i128, planes[2][1] as i128, planes[2][2] as i128],
        ];
        let rhs = [
            -(planes[0][3] as i128),
            -(planes[1][3] as i128),
            -(planes[2][3] as i128),
        ];
        let d = det3_i128(n);
        prop_assume!(d != 0);
        let col = |k: usize| {
            let mut m = n;
            m[0][k] = rhs[0];
            m[1][k] = rhs[1];
            m[2][k] = rhs[2];
            m
        };
        let (dx, dy, dz) = (det3_i128(col(0)), det3_i128(col(1)), det3_i128(col(2)));
        // p-verification: (dx/d, dy/d, dz/d) lies on all three planes.
        for pl in &planes {
            let (a, b, c, dd) = (pl[0] as i128, pl[1] as i128, pl[2] as i128, pl[3] as i128);
            prop_assert_eq!(a * dx + b * dy + c * dz + dd * d, 0);
        }
        let si = [s[0] as i128, s[1] as i128, s[2] as i128];
        let row1 = [dx - d * si[0], dy - d * si[1], dz - d * si[2]];
        let dq = [(q[0] - s[0]) as i128, (q[1] - s[1]) as i128, (q[2] - s[2]) as i128];
        let dr = [(r[0] - s[0]) as i128, (r[1] - s[1]) as i128, (r[2] - s[2]) as i128];
        let m_int = det3_i128([row1, dq, dr]);
        let expected = (d.signum() * m_int.signum()) as i8;

        let planes_f = ThreePlane(planes.map(|pl| pl.map(|v| v as f64)));
        let got = indirect_orient3d(
            &planes_f,
            q.map(|v| v as f64),
            r.map(|v| v as f64),
            s.map(|v| v as f64),
        );
        prop_assert_eq!(got, expected);
    }

    /// Independent ground truth for the parts the i128 oracle *shares* with the
    /// implementation (Row1 assembly, the `sign(D)·sign(M)` decomposition, which
    /// point is subtracted): build an integer point `p` and integer planes
    /// through it, so `p_f64 = p` exactly, then compare to
    /// `orient3d(p_f64, q, r, s)` — Shewchuk-direct, no decomposition, no
    /// conditioning worry.
    #[test]
    fn prop_matches_materialized_via_integer_point(
        p in prop::array::uniform3(-20i64..=20),
        normals in prop::array::uniform3(prop::array::uniform3(-20i64..=20)),
        q in prop::array::uniform3(-50i64..=50),
        r in prop::array::uniform3(-50i64..=50),
        s in prop::array::uniform3(-50i64..=50),
    ) {
        let ni = normals.map(|nn| [nn[0] as i128, nn[1] as i128, nn[2] as i128]);
        prop_assume!(det3_i128(ni) != 0); // planes meet only at p
        let planes = ThreePlane([
            plane_through(normals[0], p),
            plane_through(normals[1], p),
            plane_through(normals[2], p),
        ]);
        let pf = [p[0] as f64, p[1] as f64, p[2] as f64];
        let qf = q.map(|v| v as f64);
        let rf = r.map(|v| v as f64);
        let sf = s.map(|v| v as f64);
        let expected = sign_f64(orient3d(pf, qf, rf, sf));
        prop_assert_eq!(indirect_orient3d(&planes, qf, rf, sf), expected);
    }

    /// The two-implicit comparator against the exact i128 rational. `Na/Da < Nb/Db`
    /// is compared by cross-multiplication *there too*, but through a different
    /// route: i128 integers rather than expansion arithmetic, and the sign of the
    /// product `Da·Db` rather than a product of two signs.
    #[test]
    fn prop_indirect_cmp_coord_matches_i128(
        pa in prop::array::uniform3(prop::array::uniform4(-100i64..=100)),
        pb in prop::array::uniform3(prop::array::uniform4(-100i64..=100)),
        axis in 0usize..3,
    ) {
        let cramer_i128 = |pl: [[i64; 4]; 3]| {
            let n = pl.map(|row| [row[0] as i128, row[1] as i128, row[2] as i128]);
            let rhs = [-(pl[0][3] as i128), -(pl[1][3] as i128), -(pl[2][3] as i128)];
            let col = |k: usize| {
                let mut m = n;
                m[0][k] = rhs[0];
                m[1][k] = rhs[1];
                m[2][k] = rhs[2];
                m
            };
            ([det3_i128(col(0)), det3_i128(col(1)), det3_i128(col(2))], det3_i128(n))
        };
        let (na, da) = cramer_i128(pa);
        let (nb, db) = cramer_i128(pb);
        prop_assume!(da != 0 && db != 0);
        // a[axis] − b[axis] = (Na·Db − Nb·Da) / (Da·Db).
        let expected = ((na[axis] * db - nb[axis] * da).signum() * (da * db).signum()) as i8;

        let ta = ThreePlane(pa.map(|row| row.map(|v| v as f64)));
        let tb = ThreePlane(pb.map(|row| row.map(|v| v as f64)));
        prop_assert_eq!(indirect_cmp_coord(&ta, &tb, axis), expected);
    }

    /// A point does not precede itself, however it is described. Swapping two of a
    /// triple's planes negates `D`, so this also says the comparator is invariant
    /// under the representation.
    ///
    /// It does **not** guard the `sign(Da)·sign(Db)` factor, which was measured: with
    /// two descriptions of one point the numerator is exactly zero, and zero times a
    /// wrong sign is still zero. `opposite_denominators_keep_the_order` below is what
    /// catches that, along with both cross-checks against the i128 and materialized
    /// oracles.
    #[test]
    fn prop_a_point_does_not_precede_itself(
        p in prop::array::uniform3(-20i64..=20),
        normals in prop::array::uniform3(prop::array::uniform3(-20i64..=20)),
    ) {
        let ni = normals.map(|nn| [nn[0] as i128, nn[1] as i128, nn[2] as i128]);
        prop_assume!(det3_i128(ni) != 0);
        let rows = [
            plane_through(normals[0], p),
            plane_through(normals[1], p),
            plane_through(normals[2], p),
        ];
        let straight = ThreePlane(rows);
        let swapped = ThreePlane([rows[1], rows[0], rows[2]]); // D → −D
        for axis in 0..3 {
            prop_assert_eq!(indirect_cmp_coord(&straight, &swapped, axis), 0);
            prop_assert_eq!(indirect_cmp_coord(&straight, &straight, axis), 0);
        }
    }

    /// Ground truth from the other side: build two integer points and integer planes
    /// through each, so both coordinates are exactly representable, then compare the
    /// f64 coordinates directly.
    #[test]
    fn prop_indirect_cmp_coord_matches_materialized(
        pa in prop::array::uniform3(-20i64..=20),
        pb in prop::array::uniform3(-20i64..=20),
        na in prop::array::uniform3(prop::array::uniform3(-20i64..=20)),
        nb in prop::array::uniform3(prop::array::uniform3(-20i64..=20)),
        axis in 0usize..3,
    ) {
        let det = |nn: [[i64; 3]; 3]| det3_i128(nn.map(|r| [r[0] as i128, r[1] as i128, r[2] as i128]));
        prop_assume!(det(na) != 0 && det(nb) != 0);
        let ta = ThreePlane([
            plane_through(na[0], pa),
            plane_through(na[1], pa),
            plane_through(na[2], pa),
        ]);
        let tb = ThreePlane([
            plane_through(nb[0], pb),
            plane_through(nb[1], pb),
            plane_through(nb[2], pb),
        ]);
        let expected = (pa[axis] - pb[axis]).signum() as i8;
        prop_assert_eq!(indirect_cmp_coord(&ta, &tb, axis), expected);
    }

    /// Scaling one plane's coefficients by λ (negative included) leaves the
    /// result unchanged (`D → λD`, `Row1 → λRow1`). λ is a power of two so the
    /// scaled coefficients are exact and the invariance is exact — also confirms
    /// geom need not normalize normals.
    #[test]
    fn prop_scaling_a_plane_is_invariant(
        planes in prop::array::uniform3(prop::array::uniform4(-50.0f64..50.0)),
        q in prop::array::uniform3(-50.0f64..50.0),
        r in prop::array::uniform3(-50.0f64..50.0),
        s in prop::array::uniform3(-50.0f64..50.0),
        which in 0usize..3,
        lambda in prop::sample::select(vec![-2.0f64, -1.0, -0.5, 0.5, 2.0, 4.0]),
    ) {
        let base = ThreePlane(planes);
        let mut scaled = planes;
        for coeff in &mut scaled[which] {
            *coeff *= lambda;
        }
        prop_assert_eq!(
            indirect_orient3d(&base, q, r, s),
            indirect_orient3d(&ThreePlane(scaled), q, r, s)
        );
    }

    /// Localize the expansion×expansion product: both factors are multi-component
    /// `det3` results, so the `Σⱼ scale + add` accumulation is exercised (a
    /// 2×2-component product would not catch accumulation bugs). Exact value is
    /// checked against the i128 product.
    #[test]
    fn prop_mul_matches_i128(
        a in prop::array::uniform3(prop::array::uniform3(-200i64..=200)),
        b in prop::array::uniform3(prop::array::uniform3(-200i64..=200)),
    ) {
        let ea = det3(a.map(|row| row.map(|v| v as f64)));
        let eb = det3(b.map(|row| row.map(|v| v as f64)));
        let ai = a.map(|row| row.map(|v| v as i128));
        let bi = b.map(|row| row.map(|v| v as i128));
        prop_assert_eq!(expansion_to_i128(&ea.mul(&eb)), det3_i128(ai) * det3_i128(bi));
    }
}

// ---- ray / segment vs triangle (M5-d1) ----
//
// Reference triangle in z = 0, wound CCW so its right-hand normal is +z:
//   v0 = (0,0,0), v1 = (1,0,0), v2 = (0,1,0).
const T0: [f64; 3] = [0.0, 0.0, 0.0];
const T1: [f64; 3] = [1.0, 0.0, 0.0];
const T2: [f64; 3] = [0.0, 1.0, 0.0];

#[test]
fn ray_forward_crossing_is_oriented_by_direction() {
    // From below, straight up through the interior: forward, +z aligned ⇒ Cross(+1).
    assert_eq!(
        ray_triangle_cross([0.25, 0.25, -1.0], [0.0, 0.0, 1.0], T0, T1, T2),
        RayCross::Cross(1)
    );
    // From above, straight down through the interior: forward, −z aligned ⇒ Cross(−1).
    assert_eq!(
        ray_triangle_cross([0.25, 0.25, 1.0], [0.0, 0.0, -1.0], T0, T1, T2),
        RayCross::Cross(-1)
    );
    // From above, going up (away): the triangle is behind ⇒ Miss.
    assert_eq!(
        ray_triangle_cross([0.25, 0.25, 1.0], [0.0, 0.0, 1.0], T0, T1, T2),
        RayCross::Miss
    );
}

#[test]
fn ray_missing_and_degenerate() {
    // Vertical line through (2,2) is outside the unit triangle ⇒ Miss.
    assert_eq!(
        ray_triangle_cross([2.0, 2.0, -1.0], [0.0, 0.0, 1.0], T0, T1, T2),
        RayCross::Miss
    );
    // Vertical line through (0.5,0) grazes edge v0-v1 (y = 0) ⇒ Degenerate.
    assert_eq!(
        ray_triangle_cross([0.5, 0.0, -1.0], [0.0, 0.0, 1.0], T0, T1, T2),
        RayCross::Degenerate
    );
    // Origin p on the triangle's plane (z = 0) ⇒ Degenerate (s0 == 0).
    assert_eq!(
        ray_triangle_cross([0.25, 0.25, 0.0], [1.0, 0.0, 0.0], T0, T1, T2),
        RayCross::Degenerate
    );
}

#[test]
fn segment_crossing_and_contacts() {
    // Below → above through the interior ⇒ Cross(+1).
    assert_eq!(
        segment_triangle_cross([0.25, 0.25, -1.0], [0.25, 0.25, 1.0], T0, T1, T2),
        SegCross::Cross(1)
    );
    // Both endpoints above ⇒ Miss.
    assert_eq!(
        segment_triangle_cross([0.25, 0.25, 1.0], [0.25, 0.25, 2.0], T0, T1, T2),
        SegCross::Miss
    );
    // An endpoint lands on the plane ⇒ Degenerate (touching contact).
    assert_eq!(
        segment_triangle_cross([0.25, 0.25, -1.0], [0.25, 0.25, 0.0], T0, T1, T2),
        SegCross::Degenerate
    );
    // Off to the side ⇒ Miss.
    assert_eq!(
        segment_triangle_cross([2.0, 2.0, -1.0], [2.0, 2.0, 1.0], T0, T1, T2),
        SegCross::Miss
    );
    // Parallel-coplanar to edge v0-v1 (both in plane y = 0) but above the
    // triangle's plane: no plane crossing ⇒ Miss, NOT a false Degenerate
    // (the ordering fix that matters for axis-aligned input).
    assert_eq!(
        segment_triangle_cross([0.2, 0.0, 0.5], [0.8, 0.0, 0.5], T0, T1, T2),
        SegCross::Miss
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// A forward ray and its reverse can never *both* be a forward crossing of
    /// the same triangle: the full line meets the triangle's plane once, so at
    /// most one half-line reaches it (`n_cross ∈ {0,1}`). Integer coords keep
    /// every `orient3d` exact.
    #[test]
    fn prop_opposite_rays_not_both_forward(
        p in prop::array::uniform3(-40i64..=40),
        d in prop::array::uniform3(-40i64..=40),
    ) {
        let pf = p.map(|v| v as f64);
        let df = d.map(|v| v as f64);
        let dn = [-df[0], -df[1], -df[2]];
        let (a, b, c) = ([0.0, 0.0, 0.0], [9.0, 0.0, 0.0], [0.0, 9.0, 0.0]);
        let fwd = ray_triangle_cross(pf, df, a, b, c);
        let bwd = ray_triangle_cross(pf, dn, a, b, c);
        prop_assume!(fwd != RayCross::Degenerate && bwd != RayCross::Degenerate);
        let n_cross = [fwd, bwd].iter().filter(|c| matches!(c, RayCross::Cross(_))).count();
        prop_assert!(n_cross <= 1, "fwd={:?} bwd={:?}", fwd, bwd);
    }

    /// A forward ray reaching the triangle and the segment from `p` to a point
    /// just past the plane agree on the oriented crossing sign — the ray and
    /// segment predicates share the same orientation convention. Aim from `p`
    /// through the triangle's interior so a crossing is guaranteed.
    #[test]
    fn prop_ray_and_segment_agree_when_aimed_through(
        p in prop::array::uniform3(-30i64..=30),
    ) {
        let pf = p.map(|v| v as f64);
        let (a, b, c) = ([0.0, 0.0, 0.0], [9.0, 0.0, 0.0], [0.0, 9.0, 0.0]);
        // Target the interior point (3,3,0); direction and a segment well past it.
        let target = [3.0, 3.0, 0.0];
        let d = [target[0] - pf[0], target[1] - pf[1], target[2] - pf[2]];
        let far = [pf[0] + 2.0 * d[0], pf[1] + 2.0 * d[1], pf[2] + 2.0 * d[2]];
        let ray = ray_triangle_cross(pf, d, a, b, c);
        let seg = segment_triangle_cross(pf, far, a, b, c);
        prop_assume!(ray != RayCross::Degenerate && seg != SegCross::Degenerate);
        // p is off the plane and aims through the interior ⇒ both cross, same sign.
        if let (RayCross::Cross(rs), SegCross::Cross(ss)) = (ray, seg) {
            prop_assert_eq!(rs, ss);
        } else {
            prop_assert!(false, "expected both to cross: ray={:?} seg={:?}", ray, seg);
        }
    }
}

// ---- the filters of indirect_plane_side and indirect_cmp_coord ----

/// Both halves of each new filter: a clear sign is answered here, and an exact zero — a point
/// on the plane, two points sharing a coordinate — is declined to the expansions.
#[test]
fn the_side_and_cmp_filters_answer_a_clear_sign_and_decline_a_zero() {
    let at = |p: [f64; 3]| {
        ThreePlane([
            [1.0, 0.0, 0.0, -p[0]],
            [0.0, 1.0, 0.0, -p[1]],
            [0.0, 0.0, 1.0, -p[2]],
        ])
    };
    let a = at([1.0, 2.0, 3.0]);
    assert_eq!(indirect_plane_side_filter(&a, [2.0, -1.0, 1.0, -3.0]), None);
    assert_eq!(
        indirect_plane_side_filter(&a, [0.0, 0.0, 1.0, 0.0]),
        Some(indirect_plane_side_exact(&a, [0.0, 0.0, 1.0, 0.0]))
    );
    let b = at([1.0, 5.0, 0.0]);
    assert_eq!(indirect_cmp_coord_filter(&a, &b, 0), None);
    assert_eq!(
        indirect_cmp_coord_filter(&a, &b, 1),
        Some(indirect_cmp_coord_exact(&a, &b, 1))
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// The side filter's one obligation, over nine decades of scale: when it answers, the
    /// expansions agree, and it never claims a zero.
    #[test]
    fn prop_the_side_filter_never_disagrees_with_the_exact_path(
        planes in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        c in prop::array::uniform4(-1e3f64..1e3),
        scale in -4i32..5,
    ) {
        let k = 10f64.powi(scale);
        let tp = ThreePlane(planes.map(|pl| pl.map(|x| x * k)));
        let c = c.map(|x| x * k);
        if let Some(fast) = indirect_plane_side_filter(&tp, c) {
            prop_assert_eq!(fast, indirect_plane_side_exact(&tp, c));
            prop_assert_ne!(fast, 0);
        }
    }

    /// Where the bound is tight: the fourth plane passes through the implicit point **as `f64`
    /// computes it**, so `S` is a few roundings from zero and its sign is the rounding's. A filter
    /// whose bound were too small would answer here with that sign; the expansions disagree.
    #[test]
    fn prop_the_side_filter_agrees_through_the_rounded_point(
        planes in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        m in prop::array::uniform3(-1e3f64..1e3),
    ) {
        let tp = ThreePlane(planes);
        let (num, d, _, _) = cramer_val(&tp);
        prop_assume!(d != 0.0);
        let x = num.map(|n| n / d);
        let c = [m[0], m[1], m[2], -(m[0] * x[0] + m[1] * x[1] + m[2] * x[2])];
        if let Some(fast) = indirect_plane_side_filter(&tp, c) {
            prop_assert_eq!(fast, indirect_plane_side_exact(&tp, c));
        }
    }

    /// The comparison filter's obligation, over nine decades of scale.
    #[test]
    fn prop_the_cmp_filter_never_disagrees_with_the_exact_path(
        a in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        b in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        axis in 0usize..3,
        scale in -4i32..5,
    ) {
        let k = 10f64.powi(scale);
        let (ta, tb) = (
            ThreePlane(a.map(|pl| pl.map(|x| x * k))),
            ThreePlane(b.map(|pl| pl.map(|x| x * k))),
        );
        if let Some(fast) = indirect_cmp_coord_filter(&ta, &tb, axis) {
            prop_assert_eq!(fast, indirect_cmp_coord_exact(&ta, &tb, axis));
            prop_assert_ne!(fast, 0);
        }
    }

    /// Where it is tight: the second point sits on the axis plane at the first point's
    /// coordinate **as `f64` computes it**, so the two coordinates differ by a rounding.
    #[test]
    fn prop_the_cmp_filter_agrees_at_the_rounded_coordinate(
        a in prop::array::uniform3(prop::array::uniform4(-1e3f64..1e3)),
        rest in prop::array::uniform2(prop::array::uniform4(-1e3f64..1e3)),
        axis in 0usize..3,
    ) {
        let ta = ThreePlane(a);
        let (num, d, _, _) = cramer_val(&ta);
        prop_assume!(d != 0.0);
        let mut axis_plane = [0.0; 4];
        axis_plane[axis] = 1.0;
        axis_plane[3] = -(num[axis] / d);
        let tb = ThreePlane([axis_plane, rest[0], rest[1]]);
        prop_assume!(cramer_val(&tb).1 != 0.0);
        if let Some(fast) = indirect_cmp_coord_filter(&ta, &tb, axis) {
            prop_assert_eq!(fast, indirect_cmp_coord_exact(&ta, &tb, axis));
        }
    }
}
