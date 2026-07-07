//! `Point<D>` — a position in `f64` coordinates, affinely distinct from `Vector`.

use crate::Vector;
use core::array;
use core::ops::{Add, AddAssign, Index, Sub, SubAssign};

/// A `D`-dimensional position with `f64` coordinates.
///
/// Affinely distinct from [`Vector`]: `Point − Point → Vector` and
/// `Point ± Vector → Point`, but **`Point + Point` does not compile** — adding
/// two positions is geometrically meaningless. Weighted combinations go through
/// the vector algebra (see [`Point::lerp`], [`Point::centroid`]).
///
/// `PartialEq` is **exact** `f64` comparison; use it only for tests and true
/// coincidence of literal coordinates. Geometric coincidence must be judged
/// with [`Point::distance`] against a tolerance, never `==` (overview 절대원칙
/// 2). `Eq`/`Hash` are deliberately not implemented.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<const D: usize> {
    coords: [f64; D],
}

impl<const D: usize> Point<D> {
    /// Construct from raw coordinates.
    #[inline]
    pub const fn from_array(coords: [f64; D]) -> Self {
        Self { coords }
    }

    /// The origin (all zeros).
    #[inline]
    pub const fn origin() -> Self {
        Self { coords: [0.0; D] }
    }

    /// Bounds-checked coordinate access.
    #[inline]
    pub fn get(&self, i: usize) -> Option<f64> {
        self.coords.get(i).copied()
    }

    /// The coordinates as an array (by value; `Point` is `Copy`).
    #[inline]
    pub const fn as_array(&self) -> [f64; D] {
        self.coords
    }

    /// Squared distance to `rhs` (no `sqrt`).
    #[inline]
    pub fn distance_squared(self, rhs: Point<D>) -> f64 {
        (self - rhs).norm_squared()
    }

    /// Euclidean distance to `rhs`.
    #[inline]
    pub fn distance(self, rhs: Point<D>) -> f64 {
        (self - rhs).norm()
    }

    /// Linear interpolation `self + (rhs − self)·t`.
    ///
    /// Exact at `t == 0.0` (returns `self`); `t` is unclamped, so values
    /// outside `[0, 1]` extrapolate.
    #[inline]
    pub fn lerp(self, rhs: Point<D>, t: f64) -> Point<D> {
        self + (rhs - self) * t
    }

    /// The centroid (arithmetic mean position) of `points`, or `None` if empty.
    ///
    /// Computed as `p₀ + (Σ(pᵢ − p₀))/n` so it never forms a `Point + Point`.
    #[inline]
    pub fn centroid(points: &[Point<D>]) -> Option<Point<D>> {
        let (first, rest) = points.split_first()?;
        let acc = rest
            .iter()
            .fold(Vector::<D>::zero(), |acc, &p| acc + (p - *first));
        Some(*first + acc / points.len() as f64)
    }
}

impl<const D: usize> Index<usize> for Point<D> {
    type Output = f64;

    /// Panics on out-of-bounds (slice-like); use [`Point::get`] for a checked
    /// access.
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.coords[i]
    }
}

impl<const D: usize> Sub for Point<D> {
    type Output = Vector<D>;

    /// `Point − Point → Vector` (the displacement from `rhs` to `self`).
    #[inline]
    fn sub(self, rhs: Point<D>) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| self.coords[i] - rhs.coords[i]))
    }
}

impl<const D: usize> Add<Vector<D>> for Point<D> {
    type Output = Point<D>;

    #[inline]
    fn add(self, v: Vector<D>) -> Point<D> {
        Point::from_array(array::from_fn(|i| self.coords[i] + v[i]))
    }
}

impl<const D: usize> Sub<Vector<D>> for Point<D> {
    type Output = Point<D>;

    #[inline]
    fn sub(self, v: Vector<D>) -> Point<D> {
        Point::from_array(array::from_fn(|i| self.coords[i] - v[i]))
    }
}

impl<const D: usize> AddAssign<Vector<D>> for Point<D> {
    #[inline]
    fn add_assign(&mut self, v: Vector<D>) {
        *self = *self + v;
    }
}

impl<const D: usize> SubAssign<Vector<D>> for Point<D> {
    #[inline]
    fn sub_assign(&mut self, v: Vector<D>) {
        *self = *self - v;
    }
}

#[cfg(test)]
mod tests {
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
            .map(|k| {
                Point::from_array([(k & 1) as f64, ((k >> 1) & 1) as f64, ((k >> 2) & 1) as f64])
            })
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
}
