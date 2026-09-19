//! Self-contained `f64` linear algebra for the nacre kernel.
//!
//! Ships only [`Point<D>`](Point) and [`Vector<D>`](Vector), both backed by
//! `[f64; D]`. This is the plain-`f64` layer: it does the
//! *construction* arithmetic on coordinates, which are the "cache" side of the
//! truth/cache split — never the source of geometric truth.
//!
//! Scope boundaries:
//! - **Scalar is concrete `f64`, not generic.** Extended precision does not live here:
//!   exact values are `nacre-scalar`'s `Rat`, and a coordinate that needs more than
//!   `f64` is realized from its definition at a chosen precision (`nacre-ops::realize`,
//!   climbing with `astro-float`). The old note named `nacre-geom`'s relaxation ladder
//!   and double-double; neither was built.
//! - **`Transform`/matrices are deferred** to when they are first needed
//!   (sketch-plane placement, M2); a future `transform.rs` will host them.
//! - **Inputs are assumed finite.** No NaN/inf guarding in M1.
//!
//! The affine distinction between points and vectors is enforced by the type
//! system: `Point − Point → Vector`, `Point + Vector → Point`, and
//! `Point + Point` does not compile.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod point;
mod vector;

pub use point::Point;
pub use vector::Vector;

/// A 2-D displacement vector.
pub type Vector2 = Vector<2>;
/// A 3-D displacement vector.
pub type Vector3 = Vector<3>;
/// A 2-D position.
pub type Point2 = Point<2>;
/// A 3-D position.
pub type Point3 = Point<3>;

/// Combined relative-or-absolute float comparison for tests.
///
/// The absolute term handles results that should be ~0 (e.g. `a·(a×b)`), the
/// relative term handles large-magnitude results. Test-only; production code
/// compares coordinates via `distance`/tolerance, never a bare `==`.
#[cfg(test)]
pub(crate) fn approx_eq(a: f64, b: f64, rel: f64, abs: f64) -> bool {
    let diff = (a - b).abs();
    diff <= abs || diff <= rel * a.abs().max(b.abs())
}
