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

    /// The exact difference `a − b` as a two-component expansion — the exact
    /// coordinate differences (`q − s`, …) the indirect predicates cross-multiply.
    #[inline]
    pub fn two_diff(a: f64, b: f64) -> Expansion {
        Expansion(geometry_predicates::predicates::two_diff(a, b).to_vec())
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

    /// `self · other`, exact. Distributes over `other`'s components —
    /// `Σⱼ self · other[j]` via [`scale`](Self::scale) + [`add`](Self::add),
    /// each step exact, so the whole product is exact. Unavoidable for the
    /// indirect predicates: the implicit-point determinant `M` multiplies two
    /// expansions (`Row1ᵢ · crossᵢ`) that no factoring can reduce to scalars.
    pub fn mul(&self, other: &Expansion) -> Expansion {
        // `other` is never empty (the type invariant), so `other.0[0]` exists.
        let mut acc = self.scale(other.0[0]);
        for &c in &other.0[1..] {
            acc = acc.add(&self.scale(c));
        }
        acc
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

/// Three planes, each `[a, b, c, d]` meaning `a·X + b·Y + c·Z + d = 0`. When they
/// meet in a single point that point is *implicit* — the indirect predicates
/// decide signs about it without ever materializing its (generally irrational)
/// coordinates. The result is invariant under scaling any plane's coefficients,
/// so the normals need not be unit length (design §8 M5; Attene 2020).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThreePlane(pub [[f64; 4]; 3]);

/// The exact sign of `orient3d(p, q, r, s)` where `p` is the implicit point at
/// which the three planes `p` meet and `q, r, s` are explicit points — same sign
/// convention as [`orient3d`] (`+1` ⇒ `p` on the negative side of the plane
/// through `q, r, s`; `0` ⇒ the four are coplanar).
///
/// By Cramer the implicit point is `(Dx/D, Dy/D, Dz/D)` where `D` is the
/// determinant of the plane-normal matrix and `Dx/Dy/Dz` replace its column
/// `0/1/2` with `−d`. Factoring `1/D` out of the first row of
/// `det[p−s, q−s, r−s]` gives `orient3d = (1/D)·M`, hence the sign is
/// `sign(D)·sign(M)` with `M = det[(Dx−D·sx, …), q−s, r−s]`. Every part is an
/// exact polynomial in the inputs, evaluated through [`Expansion`], so the sign
/// is exact.
///
/// **Precondition:** the three planes meet in a point (`D ≠ 0`). `D = 0`
/// (parallel or coincident planes) leaves the point undefined; the result is
/// then `0` (unspecified), which falls out of `sign(D)·sign(M)` with no branch —
/// `M` is division-free, so it stays finite even when `D = 0`. `PolyhedralBoolean`
/// (M5-c) only forms valid vertices, guaranteeing `D ≠ 0` in practice.
pub fn indirect_orient3d(p: &ThreePlane, q: [f64; 3], r: [f64; 3], s: [f64; 3]) -> i8 {
    let pl = p.0;
    // Plane-normal matrix N (rows = normals) and the right-hand side −d.
    let n = [
        [pl[0][0], pl[0][1], pl[0][2]],
        [pl[1][0], pl[1][1], pl[1][2]],
        [pl[2][0], pl[2][1], pl[2][2]],
    ];
    let rhs = [-pl[0][3], -pl[1][3], -pl[2][3]];
    // Cramer determinants: D and Dx/Dy/Dz (column k replaced by rhs).
    let col_replaced = |k: usize| {
        let mut m = n;
        m[0][k] = rhs[0];
        m[1][k] = rhs[1];
        m[2][k] = rhs[2];
        m
    };
    let d = det3(n);
    let dx = det3(col_replaced(0));
    let dy = det3(col_replaced(1));
    let dz = det3(col_replaced(2));

    // cross = (q − s) × (r − s), each component an exact expansion.
    let dq = [
        Expansion::two_diff(q[0], s[0]),
        Expansion::two_diff(q[1], s[1]),
        Expansion::two_diff(q[2], s[2]),
    ];
    let dr = [
        Expansion::two_diff(r[0], s[0]),
        Expansion::two_diff(r[1], s[1]),
        Expansion::two_diff(r[2], s[2]),
    ];
    let cross = [
        dq[1].mul(&dr[2]).sub(&dq[2].mul(&dr[1])),
        dq[2].mul(&dr[0]).sub(&dq[0].mul(&dr[2])),
        dq[0].mul(&dr[1]).sub(&dq[1].mul(&dr[0])),
    ];

    // Row1 = (Dx − D·sx, Dy − D·sy, Dz − D·sz); the D·sᵢ terms are cheap scales.
    let row1 = [
        dx.sub(&d.scale(s[0])),
        dy.sub(&d.scale(s[1])),
        dz.sub(&d.scale(s[2])),
    ];
    // M = Row1 · cross — three exact expansion×expansion products.
    let m = row1[0]
        .mul(&cross[0])
        .add(&row1[1].mul(&cross[1]))
        .add(&row1[2].mul(&cross[2]));

    debug_assert!(
        d.sign() != 0,
        "indirect_orient3d: degenerate three-plane input (D = 0)"
    );
    d.sign() * m.sign()
}

/// The exact sign of `x`: `+1`, `-1`, or `0` (unlike `f64::signum`, which maps
/// `0.0` to `+1.0`).
#[inline]
fn sgn(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// The outcome of casting a **forward** ray (half-line `p + t·d`, `t > 0`) at a
/// triangle. `Cross(s)` carries the *oriented* crossing sign `s = sign(d · n)`
/// (`n` = the triangle's right-hand normal `(v1−v0)×(v2−v0)`): summing these over
/// a triangulated closed surface is its winding number about `p` (design §8 M5
/// point-in-polyhedron). `Degenerate` means the ray grazes an edge/vertex or lies
/// in the triangle's plane — an incidence `orient3d` is exactly `0` — so the
/// caller must retry with another direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RayCross {
    Cross(i8),
    Miss,
    Degenerate,
}

/// Exact forward-ray/triangle crossing, decided entirely by `orient3d` **signs**
/// (no coordinate is materialized, no `orient3d` value is subtracted — so the
/// result is exact). `q = p + d`.
///
/// The ray *line* pierces the triangle interior iff `orient3d(p, q, ·, ·)` has the
/// same nonzero sign for all three edges. The crossing is *forward* (`t > 0`) iff
/// `p` is on the side the ray recedes from: `sign(orient3d(v0,v1,v2,p))` equals
/// the direction sign `sd = sign(d·n)`. `sd` is computed exactly as
/// `−sign(orient3d(v0,v1,v2, v0+d))` (since `orient3d(v0,v1,v2,v0+d) = −(d·n)` —
/// `orient3d` is affine in its 4th point). `v0 + d` is exact for small-integer
/// directions `d` and normal-range coordinates.
pub fn ray_triangle_cross(
    p: [f64; 3],
    d: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> RayCross {
    let q = add3(p, d);
    let e0 = sgn(orient3d(p, q, v1, v2));
    let e1 = sgn(orient3d(p, q, v2, v0));
    let e2 = sgn(orient3d(p, q, v0, v1));
    if e0 == 0 || e1 == 0 || e2 == 0 {
        return RayCross::Degenerate; // ray line through an edge/vertex
    }
    if e0 != e1 || e1 != e2 {
        return RayCross::Miss; // ray line misses the triangle
    }
    let s0 = sgn(orient3d(v0, v1, v2, p));
    if s0 == 0 {
        return RayCross::Degenerate; // p on the triangle's plane (excluded upstream by the gate)
    }
    let sd = -sgn(orient3d(v0, v1, v2, add3(v0, d)));
    if sd == 0 {
        return RayCross::Degenerate; // ray parallel to the plane (cannot co-occur with a pierce)
    }
    if s0 == sd {
        RayCross::Cross(sd) // forward crossing, oriented by sd
    } else {
        RayCross::Miss // the crossing is behind p
    }
}

/// The outcome of a **finite segment** `a→b` meeting a triangle. `Cross(s)` is the
/// oriented crossing sign `s = sign((b−a)·n)`; `Degenerate` means an endpoint lies
/// on the triangle's plane or the segment grazes an edge/vertex (an incidence
/// `orient3d` is `0`) — the caller rejects such contacts (design §8 M5 coplanar /
/// edge-edge casework).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegCross {
    Cross(i8),
    Miss,
    Degenerate,
}

/// Exact segment/triangle crossing, by `orient3d` signs only. The segment's line
/// pierces the triangle interior iff `orient3d(a, b, ·, ·)` shares one nonzero
/// sign over the three edges; the segment actually crosses (parameter in `(0,1)`)
/// iff its endpoints are on strictly opposite sides of the triangle's plane
/// (`sign(orient3d(v0,v1,v2,a)) = −sign(…,b)`, both nonzero). The oriented sign is
/// then `sign(orient3d(v0,v1,v2,a))` (= `sign((b−a)·n)`).
pub fn segment_triangle_cross(
    a: [f64; 3],
    b: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> SegCross {
    let e0 = sgn(orient3d(a, b, v1, v2));
    let e1 = sgn(orient3d(a, b, v2, v0));
    let e2 = sgn(orient3d(a, b, v0, v1));
    if e0 == 0 || e1 == 0 || e2 == 0 {
        return SegCross::Degenerate; // line through an edge/vertex
    }
    if e0 != e1 || e1 != e2 {
        return SegCross::Miss; // line misses the triangle
    }
    let sa = sgn(orient3d(v0, v1, v2, a));
    let sb = sgn(orient3d(v0, v1, v2, b));
    if sa == 0 || sb == 0 {
        return SegCross::Degenerate; // an endpoint on the plane (touching/coplanar contact)
    }
    if sa == sb {
        return SegCross::Miss; // both endpoints on the same side — no crossing
    }
    SegCross::Cross(sa) // opposite sides ⇒ crosses at t ∈ (0,1), oriented by sa
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

    /// The plane `[n₀, n₁, n₂, −n·p]` through integer point `p` with normal `n`.
    fn plane_through(n: [i64; 3], p: [i64; 3]) -> [f64; 4] {
        let dot = n[0] * p[0] + n[1] * p[1] + n[2] * p[2];
        [n[0] as f64, n[1] as f64, n[2] as f64, -dot as f64]
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
    }

    proptest! {
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
}
