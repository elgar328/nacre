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

/// A Shewchuk **nonoverlapping expansion**: a list of f64 components whose exact
/// sum is the represented value, most-significant last. Built from
/// [`geometry_predicates`]' adaptive-arithmetic primitives, it lets us evaluate
/// determinant polynomials *exactly* — the value the indirect predicates (M5-a2)
/// combine (a mere sign would not compose).
///
/// The value is always exact (no adaptive fast path yet — "make it work first";
/// a Shewchuk-style error-bounded fast path is a later optimization). The inner
/// list is never empty (a zero value is `[0.0]`), so the primitives — which read
/// the first component — are always safe.
#[derive(Clone, Debug, PartialEq)]
pub struct Expansion(Vec<f64>);

impl Expansion {
    /// The exact product `a · b` as a two-component expansion.
    #[inline]
    pub fn two_product(a: f64, b: f64) -> Expansion {
        Expansion(geometry_predicates::predicates::two_product(a, b).to_vec())
    }

    /// `self · b`, exact.
    pub fn scale(&self, b: f64) -> Expansion {
        let mut h = vec![0.0; 2 * self.0.len()];
        let n = geometry_predicates::predicates::scale_expansion_zeroelim(&self.0, b, &mut h);
        h.truncate(n);
        Expansion::nonempty(h)
    }

    /// `self + other`, exact.
    pub fn add(&self, other: &Expansion) -> Expansion {
        let mut h = vec![0.0; self.0.len() + other.0.len()];
        let n =
            geometry_predicates::predicates::fast_expansion_sum_zeroelim(&self.0, &other.0, &mut h);
        h.truncate(n);
        Expansion::nonempty(h)
    }

    /// `self − other`, exact.
    pub fn sub(&self, other: &Expansion) -> Expansion {
        let neg = Expansion(other.0.iter().map(|&x| -x).collect());
        self.add(&neg)
    }

    /// The exact sign of the represented value: `+1`, `-1`, or `0`.
    ///
    /// The components are nonoverlapping, so the most-significant nonzero one
    /// carries the sign of the whole sum. (Scanning from the top is defensive
    /// against a trailing zero the zero-elimination should already have removed.)
    pub fn sign(&self) -> i8 {
        for &x in self.0.iter().rev() {
            if x > 0.0 {
                return 1;
            }
            if x < 0.0 {
                return -1;
            }
        }
        0
    }

    /// Keep the never-empty invariant: a fully-cancelled result becomes `[0.0]`.
    fn nonempty(mut h: Vec<f64>) -> Expansion {
        if h.is_empty() {
            h.push(0.0);
        }
        Expansion(h)
    }
}

/// The exact 3×3 determinant of the matrix whose rows are `m[0], m[1], m[2]`, as
/// an [`Expansion`]. Cofactor expansion along the first row:
/// `m₀₀·(m₁₁m₂₂ − m₁₂m₂₁) − m₀₁·(m₁₀m₂₂ − m₁₂m₂₀) + m₀₂·(m₁₀m₂₁ − m₁₁m₂₀)`.
///
/// This is the building block of the M5 indirect predicates: a three-plane
/// implicit point is `(Dx/D, Dy/D, Dz/D)` (Cramer), and each of `D, Dx, Dy, Dz`
/// is a 3×3 determinant whose exact *value* the indirect `orient3d` combines.
pub fn det3(m: [[f64; 3]; 3]) -> Expansion {
    // 2×2 minor `p·q − r·s`.
    let minor = |p: f64, q: f64, r: f64, s: f64| {
        Expansion::two_product(p, q).sub(&Expansion::two_product(r, s))
    };
    let c0 = minor(m[1][1], m[2][2], m[1][2], m[2][1]); // m₁₁m₂₂ − m₁₂m₂₁
    let c1 = minor(m[1][0], m[2][2], m[1][2], m[2][0]); // m₁₀m₂₂ − m₁₂m₂₀
    let c2 = minor(m[1][0], m[2][1], m[1][1], m[2][0]); // m₁₀m₂₁ − m₁₁m₂₀
    c0.scale(m[0][0])
        .sub(&c1.scale(m[0][1]))
        .add(&c2.scale(m[0][2]))
}

/// The exact sign of the 3×3 determinant [`det3`]: `+1`, `-1`, or `0`.
#[inline]
pub fn det3_sign(m: [[f64; 3]; 3]) -> i8 {
    det3(m).sign()
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
        /// geometry-predicates' own `orient3d(r1, r2, r3, origin)` equals
        /// `det[r1, r2, r3]`, so it is an exact-sign oracle for any f64 rows.
        #[test]
        fn prop_det3_sign_matches_orient3d(
            m in prop::array::uniform3(prop::array::uniform3(-1e6f64..1e6)),
        ) {
            let expected = orient3d(m[0], m[1], m[2], [0.0, 0.0, 0.0]).signum() as i8;
            prop_assert_eq!(det3_sign(m), expected);
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
}
