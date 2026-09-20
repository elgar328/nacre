use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn vclose<const D: usize>(a: Vector<D>, b: Vector<D>) -> bool {
    (0..D).all(|i| approx_eq(a[i], b[i], EPS, EPS))
}

// --- concrete-value ("golden") tests: pin actual numbers, not just laws ---

#[test]
fn from_array_roundtrips() {
    let v = Vector::from_array([1.0, 2.0, 3.0]);
    assert_eq!(v.as_array(), [1.0, 2.0, 3.0]);
    assert_eq!(v[0], 1.0);
    assert_eq!(v.get(2), Some(3.0));
    assert_eq!(v.get(3), None);
    assert_eq!(Vector::<3>::zero().as_array(), [0.0, 0.0, 0.0]);
}

#[test]
fn arithmetic_values() {
    let a = Vector::from_array([1.0, 2.0, 3.0]);
    let b = Vector::from_array([4.0, 5.0, 6.0]);
    assert_eq!((a + b).as_array(), [5.0, 7.0, 9.0]);
    assert_eq!((b - a).as_array(), [3.0, 3.0, 3.0]);
    assert_eq!((-a).as_array(), [-1.0, -2.0, -3.0]);
    assert_eq!((a * 2.0).as_array(), [2.0, 4.0, 6.0]);
    assert_eq!((2.0 * a).as_array(), [2.0, 4.0, 6.0]);
    assert_eq!((b / 2.0).as_array(), [2.0, 2.5, 3.0]);
}

#[test]
fn dot_value() {
    let a = Vector::from_array([1.0, 2.0, 3.0]);
    let b = Vector::from_array([4.0, 5.0, 6.0]);
    assert_eq!(a.dot(b), 32.0); // 4 + 10 + 18
}

#[test]
fn norm_and_normalize_values() {
    let v = Vector::from_array([3.0, 4.0]);
    assert_eq!(v.norm_squared(), 25.0);
    assert_eq!(v.norm(), 5.0);
    assert_eq!(v.normalize().unwrap().as_array(), [0.6, 0.8]);
    assert_eq!(Vector::<2>::zero().normalize(), None);
}

#[test]
fn cross_right_hand_rule() {
    let x = Vector::from_array([1.0, 0.0, 0.0]);
    let y = Vector::from_array([0.0, 1.0, 0.0]);
    let z = Vector::from_array([0.0, 0.0, 1.0]);
    assert_eq!(x.cross(y).as_array(), z.as_array());
    assert_eq!(y.cross(z).as_array(), x.as_array());
    assert_eq!(z.cross(x).as_array(), y.as_array());
    assert_eq!(x.cross(x).as_array(), [0.0, 0.0, 0.0]);
}

#[test]
fn any_perpendicular_values() {
    // +Z is least aligned with X → cross(X, Z) direction, normalized.
    let z = Vector::from_array([0.0, 0.0, 1.0]);
    let p = z.any_perpendicular().unwrap();
    assert!(approx_eq(p.norm(), 1.0, EPS, EPS));
    assert!(approx_eq(p.dot(z), 0.0, EPS, EPS));
    assert!(Vector::<3>::zero().any_perpendicular().is_none());
}

/// The rule is "cross with the axis of the **smallest** component", which is
/// the same one Onshape's `perpendicularVector` uses — and it is *discontinuous*
/// where two smallest components tie. No continuous choice exists (hairy ball),
/// so the tie-break is pinned instead: consumers place sketch frames with this,
/// and a change here would silently rotate every one of them by 90°.
#[test]
fn any_perpendicular_breaks_ties_towards_the_earlier_axis() {
    let unit = |v: Vector<3>| v.normalize().unwrap();
    // All three tie → x wins, so the perpendicular is x̂ × v̂.
    let v = Vector::from_array([1.0, 1.0, 1.0]);
    let want = unit(Vector::from_array([1.0, 0.0, 0.0]).cross(v));
    assert!((v.any_perpendicular().unwrap() - want).norm() < 1e-15);
    // y and z tie for smallest → y wins.
    let v = Vector::from_array([0.9, 0.3, 0.3]);
    let want = unit(Vector::from_array([0.0, 1.0, 0.0]).cross(v));
    assert!((v.any_perpendicular().unwrap() - want).norm() < 1e-15);
    // Nudging z below y flips the choice — the discontinuity, made visible.
    let v = Vector::from_array([0.9, 0.3, 0.299]);
    let want = unit(Vector::from_array([0.0, 0.0, 1.0]).cross(v));
    assert!((v.any_perpendicular().unwrap() - want).norm() < 1e-15);
    // Sign does not enter: the choice is on |components|.
    let v = Vector::from_array([0.9, -0.3, 0.3]);
    let want = unit(Vector::from_array([0.0, 1.0, 0.0]).cross(v));
    assert!((v.any_perpendicular().unwrap() - want).norm() < 1e-15);
}

// --- proptest algebraic laws ---

fn vec3() -> impl Strategy<Value = Vector<3>> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Vector::from_array)
}
fn vec2() -> impl Strategy<Value = Vector<2>> {
    prop::array::uniform2(-1e6f64..1e6f64).prop_map(Vector::from_array)
}
fn scalar() -> impl Strategy<Value = f64> {
    -1e3f64..1e3f64
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn add_commutes_exact(a in vec3(), b in vec3()) {
        prop_assert_eq!(a + b, b + a); // IEEE addition is exactly commutative
    }

    #[test]
    fn add_associates(a in vec3(), b in vec3(), c in vec3()) {
        prop_assert!(vclose((a + b) + c, a + (b + c)));
    }

    #[test]
    fn add_identity_and_inverse_exact(a in vec3()) {
        prop_assert_eq!(a + Vector::zero(), a);
        prop_assert_eq!(a + (-a), Vector::zero());
    }

    #[test]
    fn scalar_distributes_over_vec_add(a in vec3(), b in vec3(), s in scalar()) {
        prop_assert!(vclose((a + b) * s, a * s + b * s));
    }

    #[test]
    fn scalar_mul_order_exact(a in vec3(), s in scalar()) {
        prop_assert_eq!(a * s, s * a); // IEEE multiplication commutes exactly
    }

    #[test]
    fn dot_symmetric_exact(a in vec3(), b in vec3()) {
        // Same summation order + commutative products → bit-exact.
        prop_assert_eq!(a.dot(b), b.dot(a));
    }

    #[test]
    fn dot_bilinear(a in vec3(), b in vec3(), c in vec3()) {
        prop_assert!(approx_eq(a.dot(b + c), a.dot(b) + a.dot(c), 1e-6, 1e-6));
    }

    #[test]
    fn norm_squared_nonneg(a in vec3()) {
        prop_assert!(a.norm_squared() >= 0.0);
        prop_assert!(a.norm() >= 0.0);
    }

    #[test]
    fn normalize_yields_unit(a in vec3()) {
        prop_assume!(a.norm_squared() >= 1e-12);
        let u = a.normalize().unwrap();
        prop_assert!(approx_eq(u.norm(), 1.0, EPS, EPS));
        // direction preserved: u * ‖a‖ ≈ a
        prop_assert!(vclose(u * a.norm(), a));
    }

    #[test]
    fn cross_is_antisymmetric(a in vec3(), b in vec3()) {
        prop_assert!(vclose(a.cross(b), -(b.cross(a))));
    }

    #[test]
    fn cross_is_orthogonal(a in vec3(), b in vec3()) {
        let n = a.cross(b);
        // result ~0 regardless of magnitude → absolute tolerance, scaled by inputs
        let scale = a.norm() * b.norm() + 1.0;
        prop_assert!(a.dot(n).abs() <= 1e-6 * scale);
        prop_assert!(b.dot(n).abs() <= 1e-6 * scale);
    }

    #[test]
    fn any_perpendicular_is_unit_and_orthogonal(a in vec3()) {
        prop_assume!(a.norm_squared() >= 1e-12);
        let p = a.any_perpendicular().unwrap();
        prop_assert!(approx_eq(p.norm(), 1.0, EPS, EPS));
        let scale = a.norm() + 1.0;
        prop_assert!(a.dot(p).abs() <= 1e-6 * scale);
    }

    #[test]
    fn two_d_vectors_behave(a in vec2(), b in vec2()) {
        prop_assert_eq!(a + b, b + a);
        prop_assert!(a.norm_squared() >= 0.0);
    }
}
