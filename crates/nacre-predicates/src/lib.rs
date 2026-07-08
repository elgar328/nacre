//! Exact geometric predicates for the nacre CAD kernel (design §3, §8 M5).
//!
//! The M5 boolean ladder decides in/out and orientation by the **sign** of
//! determinants, and those signs must be exact and mutually consistent or the
//! combinatorial b-rep breaks (design §3 정밀도 분업). This crate is the sign
//! layer. It builds on [`geometry_predicates`] (a safe Rust port of Shewchuk's
//! adaptive-precision predicates, MIT/Apache), which exposes both the finished
//! predicates (`orient3d`) **and** the adaptive floating-point arithmetic
//! primitives (`two_product`, `two_sum`, `expansion_sum`, …). Those primitives
//! are what the coming **indirect** predicates (M5-a: a sign of a determinant
//! whose points are *implicit* — defined as plane intersections, never
//! materialized as coordinates; Attene 2020) are built from.
//!
//! **Pure numeric layer (design §9 predicate-cycle decision).** Everything here
//! takes plane coefficients and coordinates as plain `[f64; N]` arrays — never a
//! kernel `Handle`/`Surface`. This keeps the crate free of `nacre-geom`/`-topo`
//! (no `geom → predicates → geom` cycle) and extractable as a standalone crate
//! (the goal of being Rust's first open-source indirect-predicates crate, §1).
//! Callers (`nacre-geom`) convert their types to arrays at the boundary.

/// The exact sign of `orient3d` — the signed volume of the tetrahedron
/// `(a, b, c, d)`, computed as `det[a − d, b − d, c − d]`. Positive means `d`
/// lies on the negative side of the plane through `a, b, c` (i.e. `a, b, c` wind
/// counter-clockwise seen from `d`); zero means the four points are coplanar.
/// Robust to rounding (adaptive precision), so the sign is always exact.
///
/// A thin wrapper over [`geometry_predicates::orient3d`] — the seam where nacre
/// pins the base and the sign convention (which propagates to all M5 in/out
/// classification).
#[inline]
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    geometry_predicates::orient3d(a, b, c, d)
}

#[cfg(test)]
mod tests {
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
        use geometry_predicates::predicates::{two_product, two_sum};

        // Exact cases: small integers lose nothing, so the tail is zero.
        assert_eq!(two_product(3.0, 5.0), [0.0, 15.0]);
        assert_eq!(two_sum(1.0, 2.0), [0.0, 3.0]);

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
}
