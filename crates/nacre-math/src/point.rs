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
#[path = "tests/point.rs"]
mod tests;
