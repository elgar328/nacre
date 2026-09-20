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
/// never `==` (overview, principle 2). Hence `Eq`/`Hash` are deliberately not
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
#[path = "tests/vector.rs"]
mod tests;
