use super::*;
use crate::approx_eq;
use proptest::prelude::*;

const EPS: f64 = 1e-9;

fn pclose<const D: usize>(a: Point<D>, b: Point<D>) -> bool {
    (0..D).all(|i| approx_eq(a[i], b[i], EPS, EPS))
}

// --- concrete-value ("golden") tests ---

#[test]
fn ctor_and_access() {
    let p = Point::from_array([1.0, 2.0, 3.0]);
    assert_eq!(p.as_array(), [1.0, 2.0, 3.0]);
    assert_eq!(p[1], 2.0);
    assert_eq!(p.get(2), Some(3.0));
    assert_eq!(p.get(9), None);
    assert_eq!(Point::<3>::origin().as_array(), [0.0, 0.0, 0.0]);
}

#[test]
fn affine_ops_values() {
    let p = Point::from_array([1.0, 2.0, 3.0]);
    let q = Point::from_array([4.0, 6.0, 3.0]);
    // Point − Point → Vector
    assert_eq!((q - p).as_array(), [3.0, 4.0, 0.0]);
    // Point + Vector → Point, Point − Vector → Point
    let v = Vector::from_array([1.0, 1.0, 1.0]);
    assert_eq!((p + v).as_array(), [2.0, 3.0, 4.0]);
    assert_eq!((p - v).as_array(), [0.0, 1.0, 2.0]);
}

#[test]
fn distance_values() {
    let a = Point::from_array([0.0, 0.0, 0.0]);
    let b = Point::from_array([1.0, 2.0, 2.0]);
    assert_eq!(a.distance_squared(b), 9.0);
    assert_eq!(a.distance(b), 3.0);
}

#[test]
fn lerp_values() {
    let a = Point::from_array([0.0, 0.0]);
    let b = Point::from_array([2.0, 4.0]);
    assert_eq!(a.lerp(b, 0.5).as_array(), [1.0, 2.0]);
    assert_eq!(a.lerp(b, 0.0), a); // exact at t=0
}

#[test]
fn centroid_values() {
    // Unit cube corners → center (0.5, 0.5, 0.5).
    let corners: Vec<Point<3>> = (0..8)
        .map(|k| Point::from_array([(k & 1) as f64, ((k >> 1) & 1) as f64, ((k >> 2) & 1) as f64]))
        .collect();
    let c = Point::centroid(&corners).unwrap();
    assert!(pclose(c, Point::from_array([0.5, 0.5, 0.5])));
    let single = Point::from_array([7.0, 8.0, 9.0]);
    assert_eq!(Point::centroid(&[single]), Some(single));
    assert_eq!(Point::<3>::centroid(&[]), None);
}

// --- proptest affine laws ---

fn pt3() -> impl Strategy<Value = Point<3>> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Point::from_array)
}
fn vec3() -> impl Strategy<Value = Vector<3>> {
    prop::array::uniform3(-1e6f64..1e6f64).prop_map(Vector::from_array)
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn point_plus_vector_minus_point(p in pt3(), v in vec3()) {
        prop_assert!((0..3).all(|i| approx_eq(((p + v) - p)[i], v[i], EPS, EPS)));
    }

    #[test]
    fn point_plus_displacement_recovers(p in pt3(), q in pt3()) {
        prop_assert!(pclose(p + (q - p), q));
    }

    #[test]
    fn distance_is_symmetric_exact(p in pt3(), q in pt3()) {
        // norm(p−q) == norm(−(p−q)): squares are identical, so bit-exact.
        prop_assert_eq!(p.distance(q), q.distance(p));
    }

    #[test]
    fn distance_nonneg_and_squared(p in pt3(), q in pt3()) {
        prop_assert!(p.distance(q) >= 0.0);
        let d = p.distance(q);
        prop_assert!(approx_eq(p.distance_squared(q), d * d, 1e-6, 1e-6));
    }

    #[test]
    fn lerp_endpoints(p in pt3(), q in pt3()) {
        prop_assert_eq!(p.lerp(q, 0.0), p);       // exact
        prop_assert!(pclose(p.lerp(q, 1.0), q));  // within tolerance
    }

    #[test]
    fn lerp_is_symmetric(p in pt3(), q in pt3(), t in -2.0f64..2.0) {
        prop_assert!(pclose(p.lerp(q, t), q.lerp(p, 1.0 - t)));
    }
}
