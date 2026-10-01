//! Exact geometric predicates for the nacre CAD kernel.
//!
//! The boolean engine decides in/out and orientation by the **sign** of
//! determinants, and those signs must be exact and mutually consistent or the
//! combinatorial b-rep breaks. This crate is the sign
//! layer. It builds on [`geometry_predicates`] (a safe Rust port of Shewchuk's
//! adaptive-precision predicates, MIT/Apache), which exposes both the finished
//! predicates (`orient3d`) **and** the adaptive floating-point arithmetic
//! primitives (`two_product`, `two_sum`, `expansion_sum`, …). Those primitives
//! are what the **indirect** predicates are built from ([`indirect_cmp_coord`],
//! [`indirect_plane_side`]: the sign of a determinant whose points
//! are *implicit* — defined as plane intersections, never materialized as coordinates;
//! Attene 2020).
//!
//! **Pure numeric layer.** Everything here
//! takes plane coefficients and coordinates as plain `[f64; N]` arrays — never a
//! kernel `Handle`/`Surface`. This keeps the crate free of `nacre-geom`/`-topo`
//! (no `geom → predicates → geom` cycle) and extractable as a standalone crate
//! (the goal of being Rust's first open-source indirect-predicates crate).
//! Callers (`nacre-geom`) convert their types to arrays at the boundary.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
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

/// The signed area of triangle `abc` (positive when `abc` turns counter-clockwise);
/// zero means the three points are collinear. A thin wrapper over
/// [`geometry_predicates::orient2d`] — robust to rounding, so the sign is always exact.
/// Used for exact in-plane containment (drop the face normal's dominant axis first, so
/// the projection is exact — a coordinate is discarded, not recomputed).
#[inline]
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    geometry_predicates::orient2d(a, b, c)
}

/// Positive when `d` lies **strictly inside** the circle through `a`, `b`, `c` — which
/// must be given counter-clockwise — zero when the four are cocircular, negative when
/// `d` is outside. A thin wrapper over [`geometry_predicates::incircle`], so the sign is
/// exact.
///
/// This is the whole of the Delaunay condition: a triangulation is Delaunay exactly
/// when no triangle's circumcircle contains a fourth vertex, and *constrained* Delaunay
/// is the same statement over the edges left free. Flipping every free edge that fails
/// it maximises the smallest angle in the mesh (Lawson 1977; Chew 1989), which is what
/// `nacre-tess` uses it for — the sign has to be exact for the same reason `orient2d`
/// does, since a mis-signed test would flip an edge back and forth forever.
#[inline]
pub fn incircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    geometry_predicates::incircle(a, b, c, d)
}

/// A Shewchuk **nonoverlapping expansion**: a list of f64 components whose exact
/// sum is the represented value, most-significant last. Built from
/// [`geometry_predicates`]' adaptive-arithmetic primitives, it lets us evaluate
/// determinant polynomials *exactly* — the value the indirect predicates (M5-a2)
/// combine (a mere sign would not compose).
///
/// The value is always exact — the expansions carry every bit. The *sign* is what
/// gets a fast path: each indirect predicate filters in `f64` first and
/// only falls back to these expansions when the rounding bound cannot separate the
/// sign from zero. The expansion arithmetic itself is unfiltered; a filtered
/// coordinate representation would be a separate optimization. The inner list is
/// never empty (a zero value is `[0.0]`), so the primitives — which read the first
/// component — are always safe.
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

/// The exact sign of the 3×3 determinant of `m`: `+1`, `-1`, or `0`.
///
/// `det[r0, r1, r2] = orient3d(r0, r1, r2, 0)`, and `geometry_predicates::orient3d`
/// is already adaptive (filtered), so this borrows its fast path — no expansion is
/// built unless the sign is too close to call. The value (not the sign) still comes
/// from [`det3`], which [`cramer`] uses; `prop_det3_sign_matches_the_expansion` pins the
/// two together.
#[inline]
pub fn det3_sign(m: [[f64; 3]; 3]) -> i8 {
    sgn(orient3d(m[0], m[1], m[2], [0.0; 3]))
}

/// Three planes, each `[a, b, c, d]` meaning `a·X + b·Y + c·Z + d = 0`. When they
/// meet in a single point that point is *implicit* — the indirect predicates
/// decide signs about it without ever materializing its (generally irrational)
/// coordinates. The result is invariant under scaling any plane's coefficients,
/// so the normals need not be unit length (Attene 2020).
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
/// The implicit point's exact rational coordinates: numerators `(Dx, Dy, Dz)` over the
/// common denominator `D`, by Cramer. `D` is the determinant of the plane-normal matrix;
/// `Dx/Dy/Dz` replace its column `0/1/2` with `−d`.
///
/// Never public. The coordinates are only ever *compared*, never handed out — a caller
/// holding `Expansion`s would be holding a materialized point in all but name, and the
/// whole point of an implicit point is that it is never materialized.
fn cramer(p: &ThreePlane) -> ([Expansion; 3], Expansion) {
    let pl = p.0;
    // Plane-normal matrix N (rows = normals) and the right-hand side −d.
    let n = [
        [pl[0][0], pl[0][1], pl[0][2]],
        [pl[1][0], pl[1][1], pl[1][2]],
        [pl[2][0], pl[2][1], pl[2][2]],
    ];
    let rhs = [-pl[0][3], -pl[1][3], -pl[2][3]];
    let col_replaced = |k: usize| {
        let mut m = n;
        m[0][k] = rhs[0];
        m[1][k] = rhs[1];
        m[2][k] = rhs[2];
        m
    };
    (
        [
            det3(col_replaced(0)),
            det3(col_replaced(1)),
            det3(col_replaced(2)),
        ],
        det3(n),
    )
}

/// The exact sign of `a[axis] − b[axis]` for two implicit points: `-1`, `0`, `+1`.
///
/// Each coordinate is a ratio `N/D` of exact determinants ([`cramer`]), so
///
/// > `sign(Na/Da − Nb/Db) = sign(Na·Db − Nb·Da) · sign(Da) · sign(Db)`
///
/// — three exact signs, no division, nothing materialized. This is the *two*-implicit
/// predicate: the only place two implicit points meet inside one decision. A `0` means the
/// coordinates are exactly equal, which for distinct three-plane triples is real
/// information and not a tolerance question.
///
/// **Precondition:** both triples meet in a point (`D ≠ 0`).
pub fn indirect_cmp_coord(a: &ThreePlane, b: &ThreePlane, axis: usize) -> i8 {
    debug_assert!(axis < 3, "indirect_cmp_coord: axis must be 0, 1 or 2");
    if let Some(sign) = indirect_cmp_coord_filter(a, b, axis) {
        return sign;
    }
    indirect_cmp_coord_exact(a, b, axis)
}

/// [`indirect_cmp_coord`] with no filter: exact expansions, every time. The fallback, and the
/// oracle its filter is tested against.
fn indirect_cmp_coord_exact(a: &ThreePlane, b: &ThreePlane, axis: usize) -> i8 {
    let (na, da) = cramer(a);
    let (nb, db) = cramer(b);
    debug_assert!(
        da.sign() != 0 && db.sign() != 0,
        "indirect_cmp_coord: degenerate three-plane input (D = 0)"
    );
    na[axis].mul(&db).sub(&nb[axis].mul(&da)).sign() * da.sign() * db.sign()
}

/// The floating-point filter for [`indirect_cmp_coord`] — the same construction as
/// [`indirect_orient3d_filter`]. Both points' `D` and numerators are determinants (`≈ 5u`), and
/// `V = Na·Db − Nb·Da` multiplies two of them and subtracts:
///
/// | value | derived εₓ | constant used | margin |
/// |---|---|---|---|
/// | `Da`, `Db` | `≈ 5u` | `32·U` | 6× |
/// | `V = Na·Db − Nb·Da` | `≈ 5u + 5u + u + u = 12u` | `64·U` | 5× |
///
/// ★ **Why it exists** (measured): the exact route builds two Cramer expansions and a degree-six
/// product on every call — `≈ 2.7 µs`, against the filtered certified route's `≈ 1 µs` — and it
/// allocates, so under rayon the workers queue on the allocator. It was the largest single judging
/// cost of the axis-aligned fold (460 ms of about 1 s, one thread); with this filter that fold runs
/// in about 70% of its time.
#[inline]
fn indirect_cmp_coord_filter(a: &ThreePlane, b: &ThreePlane, axis: usize) -> Option<i8> {
    let (na, da, na_mag, da_mag) = cramer_val(a);
    let (nb, db, nb_mag, db_mag) = cramer_val(b);
    let v = na[axis] * db - nb[axis] * da;
    let v_mag = na_mag[axis] * db_mag + nb_mag[axis] * da_mag;
    // An overflow to infinity makes the bound meaningless; hand it to the expansions.
    if !v.is_finite() || !v_mag.is_finite() || !da_mag.is_finite() || !db_mag.is_finite() {
        return None;
    }
    (da.abs() > 32.0 * U * da_mag && db.abs() > 32.0 * U * db_mag && v.abs() > 64.0 * U * v_mag)
        .then(|| sgn(v) * sgn(da) * sgn(db))
}

/// Unit roundoff, `2^-53`: the relative error of one correctly rounded `f64` operation.
const U: f64 = f64::EPSILON / 2.0;

/// `det3` in plain `f64` — the same cofactor expansion [`det3`] evaluates exactly.
#[inline]
fn det3_val(m: [[f64; 3]; 3]) -> f64 {
    let minor = |p: f64, q: f64, r: f64, t: f64| p * q - r * t;
    m[0][0] * minor(m[1][1], m[2][2], m[1][2], m[2][1])
        - m[0][1] * minor(m[1][0], m[2][2], m[1][2], m[2][0])
        + m[0][2] * minor(m[1][0], m[2][1], m[1][1], m[2][0])
}

/// The same expression with every input replaced by its magnitude and every
/// subtraction by an addition. Two things at once: an upper bound on the exact
/// `|det3|`, and the scale against which the floating-point evaluation's rounding
/// error is measured (Higham, *Accuracy and Stability*, §3.1 — a straight-line
/// program of `n` rounded operations on exact inputs errs by at most `γₙ` times
/// this cancellation-free evaluation).
#[inline]
fn det3_mag(m: [[f64; 3]; 3]) -> f64 {
    let minor = |p: f64, q: f64, r: f64, t: f64| p.abs() * q.abs() + r.abs() * t.abs();
    m[0][0].abs() * minor(m[1][1], m[2][2], m[1][2], m[2][1])
        + m[0][1].abs() * minor(m[1][0], m[2][2], m[1][2], m[2][0])
        + m[0][2].abs() * minor(m[1][0], m[2][1], m[1][1], m[2][0])
}

/// [`cramer`] in `f64`, alongside the cancellation-free magnitudes of each part.
#[inline]
fn cramer_val(p: &ThreePlane) -> ([f64; 3], f64, [f64; 3], f64) {
    let pl = p.0;
    let n = [
        [pl[0][0], pl[0][1], pl[0][2]],
        [pl[1][0], pl[1][1], pl[1][2]],
        [pl[2][0], pl[2][1], pl[2][2]],
    ];
    let rhs = [-pl[0][3], -pl[1][3], -pl[2][3]];
    let col_replaced = |k: usize| {
        let mut m = n;
        m[0][k] = rhs[0];
        m[1][k] = rhs[1];
        m[2][k] = rhs[2];
        m
    };
    let num = [
        det3_val(col_replaced(0)),
        det3_val(col_replaced(1)),
        det3_val(col_replaced(2)),
    ];
    let num_mag = [
        det3_mag(col_replaced(0)),
        det3_mag(col_replaced(1)),
        det3_mag(col_replaced(2)),
    ];
    (num, det3_val(n), num_mag, det3_mag(n))
}

/// **Which side of the plane `c` the implicit point `p` lies on** — `+1` on the side its normal
/// `(c₀, c₁, c₂)` points to, `-1` on the other, `0` exactly on it.
///
/// Asked with the fourth plane's **coefficients**, not three of its points. That is not a
/// convenience: a plane described twice — by
/// coefficients and by a triangle — is described by two planes whenever `d` was a rounded product,
/// and handing one predicate both descriptions is how a boolean comes to disagree with itself. With
/// no triangle there is nothing left to disagree.
///
/// ★ **It is also cheaper.** `X = Dvec/D`, so
///
/// > `sign(c·X + c₃) = sign(c·Dvec + c₃·D) · sign(D)`
///
/// — a four-term dot product where the triangle form needs three differences, a cross product and
/// another dot. Division-free and exact, through [`Expansion`].
///
/// ★★ **The sign is the plane's own, not a face's.** A face whose stored normal opposes its
/// outward direction has to apply that relation itself (`nacre_ops`' `frame_sign`); this states
/// where the point is relative to the coefficients it was given, and nothing else.
///
/// **Precondition:** the three planes meet in a point (`D ≠ 0`).
pub fn indirect_plane_side(p: &ThreePlane, c: [f64; 4]) -> i8 {
    if let Some(sign) = indirect_plane_side_filter(p, c) {
        return sign;
    }
    indirect_plane_side_exact(p, c)
}

/// [`indirect_plane_side`] with no filter: exact expansions, every time. The fallback, and the
/// oracle its filter is tested against.
fn indirect_plane_side_exact(p: &ThreePlane, c: [f64; 4]) -> i8 {
    let ([dx, dy, dz], d) = cramer(p);
    let side = dx
        .scale(c[0])
        .add(&dy.scale(c[1]))
        .add(&dz.scale(c[2]))
        .add(&d.scale(c[3]));
    side.sign() * d.sign()
}

/// The floating-point filter for [`indirect_plane_side`] — the same shape as
/// [`indirect_orient3d_filter`], whose table this extends: `D` and the numerators are the same
/// determinants (`≈ 5u` each), and `S = c·N + c₃·D` adds one product and three additions on exact
/// `c`:
///
/// | value | derived εₓ | constant used | margin |
/// |---|---|---|---|
/// | `D` | `≈ 5u` | `32·U` | 6× |
/// | `S = c·N + c₃·D` | `≈ 5u + u + 3u = 9u` | `64·U` | 7× |
///
/// ★ **Why it exists** (measured): without it every call builds expansions, and once an unmoved
/// plane's `orient3d` asked its fourth plane by coefficients the axis-aligned fold ran **14×**
/// slower; with it 21.8 M calls of that fold decided here and none reached the expansions.
#[inline]
fn indirect_plane_side_filter(p: &ThreePlane, c: [f64; 4]) -> Option<i8> {
    let (num, d, num_mag, d_mag) = cramer_val(p);
    let side = c[0] * num[0] + c[1] * num[1] + c[2] * num[2] + c[3] * d;
    let side_mag = c[0].abs() * num_mag[0]
        + c[1].abs() * num_mag[1]
        + c[2].abs() * num_mag[2]
        + c[3].abs() * d_mag;
    // An overflow to infinity makes the bound meaningless; hand it to the expansions.
    if !side.is_finite() || !side_mag.is_finite() || !d_mag.is_finite() {
        return None;
    }
    (d.abs() > 32.0 * U * d_mag && side.abs() > 64.0 * U * side_mag).then(|| sgn(d) * sgn(side))
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
/// a triangulated closed surface is its winding number about `p` (the
/// point-in-polyhedron test). `Degenerate` means the ray grazes an edge/vertex or lies
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
/// `orient3d` is `0`) — the caller rejects such contacts (coplanar /
/// edge-edge casework).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegCross {
    Cross(i8),
    Miss,
    Degenerate,
}

/// Exact segment/triangle crossing, by `orient3d` signs only.
///
/// The plane test comes **first**: unless the endpoints straddle the triangle's
/// plane (`sign(orient3d(v0,v1,v2,a)) = −sign(…,b)`, both nonzero) the segment
/// cannot cross the triangle, so it is a `Miss` (or `Degenerate` if an endpoint
/// lies on the plane) — without consulting the edge tests. This ordering matters
/// for axis-aligned input: a segment parallel-coplanar to a triangle edge makes
/// `orient3d(a, b, ·, ·) = 0` yet is not a real contact, and must not be reported
/// as `Degenerate`. Only once the segment genuinely pierces the plane do the edge
/// tests decide whether the crossing point is inside the triangle; a zero there
/// is a true edge/vertex hit ⇒ `Degenerate`. The oriented sign is
/// `sign(orient3d(v0,v1,v2,a))` (= `sign((b−a)·n)`).
pub fn segment_triangle_cross(
    a: [f64; 3],
    b: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> SegCross {
    let sa = sgn(orient3d(v0, v1, v2, a));
    let sb = sgn(orient3d(v0, v1, v2, b));
    if sa == 0 || sb == 0 {
        return SegCross::Degenerate; // an endpoint on the plane (touching/coplanar contact)
    }
    if sa == sb {
        return SegCross::Miss; // both endpoints on one side — no plane crossing
    }
    // The segment pierces the plane at t ∈ (0,1); is that point inside the triangle?
    let e0 = sgn(orient3d(a, b, v1, v2));
    let e1 = sgn(orient3d(a, b, v2, v0));
    let e2 = sgn(orient3d(a, b, v0, v1));
    if e0 == 0 || e1 == 0 || e2 == 0 {
        return SegCross::Degenerate; // crossing point on a triangle edge/vertex
    }
    if e0 == e1 && e1 == e2 {
        SegCross::Cross(sa) // inside the triangle
    } else {
        SegCross::Miss // pierces the plane outside the triangle
    }
}

#[cfg(test)]
#[path = "tests/lib.rs"]
mod tests;
