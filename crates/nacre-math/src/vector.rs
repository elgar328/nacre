//! `Vector<D>` — a `D`-dimensional displacement in `f64` coordinates.

use core::array;
use core::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};

/// A `D`-dimensional vector (displacement) with `f64` components.
///
/// Distinct from [`Point`](crate::Point): a `Vector` is a free displacement,
/// so `Vector ± Vector → Vector` and scaling are defined, but it has no
/// position.
///
/// `PartialEq` is **exact** `f64` comparison. Coordinates are the cache side of
/// the truth/cache split, so exact `==` is meaningful only for tests and true-
/// zero checks; geometric coincidence must be judged with `distance`/tolerance,
/// never `==` (overview 절대원칙 2). Hence `Eq`/`Hash` are deliberately not
/// implemented (`f64` has neither, and identity lives on `Handle`, not coords).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<const D: usize> {
    coords: [f64; D],
}

impl<const D: usize> Vector<D> {
    /// Construct from raw components.
    #[inline]
    pub const fn from_array(coords: [f64; D]) -> Self {
        Self { coords }
    }

    /// The zero vector.
    #[inline]
    pub const fn zero() -> Self {
        Self { coords: [0.0; D] }
    }

    /// Bounds-checked component access.
    #[inline]
    pub fn get(&self, i: usize) -> Option<f64> {
        self.coords.get(i).copied()
    }

    /// The components as an array (by value; `Vector` is `Copy`).
    #[inline]
    pub const fn as_array(&self) -> [f64; D] {
        self.coords
    }

    /// Dot product `Σ selfᵢ·rhsᵢ`.
    #[inline]
    pub fn dot(self, rhs: Vector<D>) -> f64 {
        self.coords
            .iter()
            .zip(&rhs.coords)
            .map(|(a, b)| a * b)
            .sum()
    }

    /// Squared Euclidean length (no `sqrt`); prefer this when only comparisons
    /// are needed.
    #[inline]
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean length `‖self‖`.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Unit vector in the same direction, or `None` for the zero vector.
    ///
    /// Returns `None` iff `norm_squared` is not positive. Ordered `> 0.0` (not
    /// `== 0.0`) both avoids `clippy::float_cmp` and reads as "only normalize a
    /// positive-length vector". Note: components so small that `norm_squared`
    /// underflows to `0.0` also yield `None` — out of scope for M1's finite,
    /// normally-scaled inputs.
    #[inline]
    pub fn normalize(self) -> Option<Vector<D>> {
        let n2 = self.norm_squared();
        if n2 > 0.0 {
            Some(self / n2.sqrt())
        } else {
            None
        }
    }
}

impl Vector<3> {
    /// Cross product (3-D only; `Vector<2>::cross` does not exist).
    #[inline]
    pub fn cross(self, rhs: Vector<3>) -> Vector<3> {
        let a = self.coords;
        let b = rhs.coords;
        Vector::from_array([
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ])
    }

    /// An arbitrary **unit** vector perpendicular to `self`.
    ///
    /// Crosses `self` with the coordinate axis it is least aligned with (so the
    /// cross is never near-zero) and normalizes. The choice is arbitrary — any
    /// perpendicular will do — so this suits synthesizing a reference/seam
    /// direction. Returns `None` iff `self` is the zero vector.
    pub fn any_perpendicular(self) -> Option<Vector<3>> {
        let a = self.coords.map(f64::abs);
        let axis = if a[0] <= a[1] && a[0] <= a[2] {
            Vector::from_array([1.0, 0.0, 0.0])
        } else if a[1] <= a[2] {
            Vector::from_array([0.0, 1.0, 0.0])
        } else {
            Vector::from_array([0.0, 0.0, 1.0])
        };
        axis.cross(self).normalize()
    }
}

impl<const D: usize> Index<usize> for Vector<D> {
    type Output = f64;

    /// Panics on out-of-bounds (slice-like); use [`Vector::get`] for a checked
    /// access.
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.coords[i]
    }
}

impl<const D: usize> Add for Vector<D> {
    type Output = Vector<D>;

    #[inline]
    fn add(self, rhs: Vector<D>) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| self.coords[i] + rhs.coords[i]))
    }
}

impl<const D: usize> Sub for Vector<D> {
    type Output = Vector<D>;

    #[inline]
    fn sub(self, rhs: Vector<D>) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| self.coords[i] - rhs.coords[i]))
    }
}

impl<const D: usize> Neg for Vector<D> {
    type Output = Vector<D>;

    #[inline]
    fn neg(self) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| -self.coords[i]))
    }
}

impl<const D: usize> Mul<f64> for Vector<D> {
    type Output = Vector<D>;

    #[inline]
    fn mul(self, s: f64) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| self.coords[i] * s))
    }
}

/// `f64 * Vector` — the scalar-on-the-left form (orphan rule permits this
/// because `Vector` appears in the trait parameters).
impl<const D: usize> Mul<Vector<D>> for f64 {
    type Output = Vector<D>;

    #[inline]
    fn mul(self, v: Vector<D>) -> Vector<D> {
        v * self
    }
}

impl<const D: usize> Div<f64> for Vector<D> {
    type Output = Vector<D>;

    #[inline]
    fn div(self, s: f64) -> Vector<D> {
        Vector::from_array(array::from_fn(|i| self.coords[i] / s))
    }
}

impl<const D: usize> AddAssign for Vector<D> {
    #[inline]
    fn add_assign(&mut self, rhs: Vector<D>) {
        *self = *self + rhs;
    }
}

impl<const D: usize> SubAssign for Vector<D> {
    #[inline]
    fn sub_assign(&mut self, rhs: Vector<D>) {
        *self = *self - rhs;
    }
}

impl<const D: usize> MulAssign<f64> for Vector<D> {
    #[inline]
    fn mul_assign(&mut self, s: f64) {
        *self = *self * s;
    }
}

impl<const D: usize> DivAssign<f64> for Vector<D> {
    #[inline]
    fn div_assign(&mut self, s: f64) {
        *self = *self / s;
    }
}

#[cfg(test)]
mod tests {
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
}
