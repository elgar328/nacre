//! [`QuadVal`] — a quadratic algebraic scalar `a + b·√c` (a, b, c rational, c ≥ 0), and the
//! exact sign tower over it (M6-1).
//!
//! **Why this type exists.** The intersection line of two rational planes meets a rational
//! cylinder where a rational quadratic `As² + Bs + C = 0` vanishes, so every coordinate of a
//! plane·plane·cylinder vertex is `a + b√c` with one shared radical. Point-versus-plane signs
//! close over [`QuadVal::sign`] (one radical); the circular order of events on a cylinder face
//! compares points born of *different* cutting planes — different discriminants — and closes
//! over [`biquad_sign`] (the 4-term element of ℚ(√u,√v)). Both reduce to rational sign
//! questions by the classic recursion; no general algebraic-number machinery is needed.
//!
//! **The sign tower is total.** Sign is the heart of judgment, so it must not refuse on
//! overflow: denominators are cleared and the recursion runs in `BigInt` ("integers, not
//! rationals" — the `three_planes_big` precedent; the lcm itself is a `BigInt`, so there is no
//! failure point). The *value* arithmetic (`checked_add`/`checked_mul`…) stays in checked
//! `Rat` and answers `None` on overflow — the caller sees the downgrade trigger, exactly like
//! `Rat` itself.
//!
//! **No `PartialEq`.** `1 + 2√2` and `1 + 1√8` are one value in two spellings, and collapsing
//! the spelling needs the square-free part of `c` — a factorization, whose cost is unbounded.
//! Value equality is the sign tower's job (`x.checked_sub(y)?.sign() == Orient::Zero`);
//! structural equality is deliberately not offered until a consumer needs it by name.

use crate::{Orient, Rat};
use num_bigint::BigInt;

/// A quadratic algebraic scalar `a + b·√c`, exact.
///
/// Invariants (kept by the constructors and every operation):
/// - `c ≥ 0` — a negative radicand means the value does not exist (a missed intersection has
///   no coordinates), so it is unrepresentable here; [`QuadVal::new`] refuses it.
/// - `b == 0 ⇔ c == 0` — a vanishing radical part is stored in one spelling (hygiene; this
///   does **not** make structural equality value equality, see the module docs).
///
/// Arithmetic is same-radical only: ℚ(√c) is closed under `+`/`−`/`×`, and mixing radicals is
/// a caller bug (debug-asserted, `None` in release). A rational (`c == 0`) is compatible with
/// every radical.
#[derive(Clone, Copy, Debug)]
pub struct QuadVal {
    a: Rat,
    b: Rat,
    c: Rat,
}

impl QuadVal {
    /// The checked constructor — `None` for `c < 0` (the value would not be real).
    pub fn new(a: Rat, b: Rat, c: Rat) -> Option<Self> {
        let zero = Rat::from_int(0);
        if c < zero {
            return None;
        }
        if b == zero || c == zero {
            return Some(QuadVal {
                a,
                b: zero,
                c: zero,
            });
        }
        Some(QuadVal { a, b, c })
    }

    /// A rational, embedded (`a + 0·√0`).
    pub fn from_rat(a: Rat) -> Self {
        let zero = Rat::from_int(0);
        QuadVal {
            a,
            b: zero,
            c: zero,
        }
    }

    /// The rational part `a`.
    pub fn a(&self) -> Rat {
        self.a
    }

    /// The radical coefficient `b` (zero iff the value is rational).
    pub fn b(&self) -> Rat {
        self.b
    }

    /// The radicand `c` (non-negative; zero iff the value is rational).
    pub fn c(&self) -> Rat {
        self.c
    }

    /// The shared radical two operands agree on — `None` when they genuinely differ (a caller
    /// bug: same-radical arithmetic is this type's contract). A rational operand (c = 0) is
    /// compatible with anything.
    fn common_radical(&self, rhs: &QuadVal) -> Option<Rat> {
        let zero = Rat::from_int(0);
        if self.c == zero {
            Some(rhs.c)
        } else if rhs.c == zero || self.c == rhs.c {
            Some(self.c)
        } else {
            debug_assert!(
                false,
                "same-radical arithmetic mixed √{:?} with √{:?}",
                self.c, rhs.c
            );
            None
        }
    }

    /// Exact addition; `None` on `Rat` overflow or a radical mismatch.
    pub fn checked_add(&self, rhs: &QuadVal) -> Option<QuadVal> {
        let c = self.common_radical(rhs)?;
        QuadVal::new(self.a.checked_add(rhs.a)?, self.b.checked_add(rhs.b)?, c)
    }

    /// Exact subtraction; `None` on `Rat` overflow or a radical mismatch.
    pub fn checked_sub(&self, rhs: &QuadVal) -> Option<QuadVal> {
        let c = self.common_radical(rhs)?;
        QuadVal::new(self.a.checked_sub(rhs.a)?, self.b.checked_sub(rhs.b)?, c)
    }

    /// Exact negation; `None` on overflow. ★ Checked, not infallible — the first spelling
    /// claimed "a reduced `Ratio` never holds `i128::MIN`", and that is false (measured:
    /// `Rat::from_int(i128::MIN)` is a legal reduced value whose negation overflows). The
    /// claim was a panic path dressed as a proof; the checked contract is the honest one.
    pub fn checked_neg(&self) -> Option<QuadVal> {
        let zero = Rat::from_int(0);
        QuadVal::new(zero.checked_sub(self.a)?, zero.checked_sub(self.b)?, self.c)
    }

    /// Exact multiplication within one radical:
    /// `(a₁+b₁√c)(a₂+b₂√c) = (a₁a₂ + b₁b₂c) + (a₁b₂ + a₂b₁)√c`.
    /// `None` on `Rat` overflow or a radical mismatch.
    pub fn checked_mul(&self, rhs: &QuadVal) -> Option<QuadVal> {
        let c = self.common_radical(rhs)?;
        let a = self
            .a
            .checked_mul(rhs.a)?
            .checked_add(self.b.checked_mul(rhs.b)?.checked_mul(c)?)?;
        let b = self
            .a
            .checked_mul(rhs.b)?
            .checked_add(rhs.a.checked_mul(self.b)?)?;
        QuadVal::new(a, b, c)
    }

    /// Exact scaling by a rational; `None` on overflow.
    pub fn checked_mul_rat(&self, k: Rat) -> Option<QuadVal> {
        QuadVal::new(self.a.checked_mul(k)?, self.b.checked_mul(k)?, self.c)
    }

    /// The f64 realization `a + b·√c` — a **cache** value, not exactly rounded (three
    /// roundings compose); consumers that need a certified value go through the sign tower.
    pub fn to_f64(&self) -> f64 {
        self.a.to_f64() + self.b.to_f64() * self.c.to_f64().sqrt()
    }

    /// The exact sign — total (the BigInt core cannot overflow and `c ≥ 0` by invariant).
    pub fn sign(&self) -> Orient {
        // Integerize the radical first: √(p/q) = √(p·q)/q, folding 1/q into b.
        let (big_a, big_b, big_c) = integerize(self.a, self.b, self.c);
        sign1_int(&big_a, &big_b, &big_c)
    }
}

/// `a + b·√c` with the radicand made an integer and all denominators cleared — the exact
/// integer triple whose sign equals the value's. Scaling by a positive rational preserves
/// sign, so: fold `1/q` (from `√(p/q) = √(pq)/q`) into `b`, then multiply through by the
/// (positive) product of the remaining denominators.
fn integerize(a: Rat, b: Rat, c: Rat) -> (BigInt, BigInt, BigInt) {
    debug_assert!(c >= Rat::from_int(0), "QuadVal invariant: c ≥ 0");
    // √(p/q) = √(p·q) / q  (q > 0 — `Ratio` keeps denominators positive).
    let (cp, cq) = (BigInt::from(c.numer()), BigInt::from(c.denom()));
    let big_c = &cp * &cq;
    // b/q √(pq): the radical's denominator joins b's.
    let (an, ad) = (BigInt::from(a.numer()), BigInt::from(a.denom()));
    let (bn, bd_raw) = (BigInt::from(b.numer()), BigInt::from(b.denom()));
    let bd = bd_raw * cq;
    // Multiply through by ad·bd (> 0): A = an·bd, B = bn·ad.
    (an * &bd, bn * ad, big_c)
}

fn sign_big(x: &BigInt) -> Orient {
    match x.sign() {
        num_bigint::Sign::Plus => Orient::Positive,
        num_bigint::Sign::Minus => Orient::Negative,
        num_bigint::Sign::NoSign => Orient::Zero,
    }
}

fn orient_mul(x: Orient, y: Orient) -> Orient {
    match (x, y) {
        (Orient::Zero, _) | (_, Orient::Zero) => Orient::Zero,
        (Orient::Positive, Orient::Positive) | (Orient::Negative, Orient::Negative) => {
            Orient::Positive
        }
        _ => Orient::Negative,
    }
}

/// `sign(A + B·√C)` for integers, `C ≥ 0` — the base of the tower.
///
/// Case ladder (each case *first* — the order is load-bearing):
/// - `C = 0` → `sign(A)` (`B·√0 = 0` whatever `B` is — answering `sign(B)` here was the
///   spec-level bug the plan review caught).
/// - `B = 0` → `sign(A)`; `A = 0` → `sign(B)` (√C > 0 here).
/// - same sign → that sign.
/// - opposite signs → `sign(A) · sign(A² − B²·C)`: the comparison alone does not carry the
///   answer's sign (A < 0, B > 0, A² < B²C is a *positive* value while A² − B²C is negative —
///   the product is what points the right way).
fn sign1_int(a: &BigInt, b: &BigInt, c: &BigInt) -> Orient {
    use num_bigint::Sign::NoSign;
    debug_assert!(c.sign() != num_bigint::Sign::Minus, "radicand must be ≥ 0");
    if c.sign() == NoSign || b.sign() == NoSign {
        return sign_big(a);
    }
    if a.sign() == NoSign {
        return sign_big(b);
    }
    let (sa, sb) = (sign_big(a), sign_big(b));
    if sa == sb {
        return sa;
    }
    orient_mul(sa, sign_big(&(a * a - b * b * c)))
}

/// `sign(A + B√u + C√v + D√(u·v))` for rationals with `u, v ≥ 0` — the second storey of the
/// tower, the sign a circular-order comparison on a cylinder face reduces to (two points born
/// of different cutting planes carry different discriminants).
///
/// `None` iff a radicand is negative (no such real value — the caller asked about a point
/// that does not exist). Otherwise total: the recursion runs in `BigInt`.
///
/// Recursion: write the element as `P + √v·Q` with `P = A + B√u`, `Q = C + D√u` — then
/// - `sign(Q) = 0` → `sign(P)`; `sign(P) = 0` → `sign(Q)` (√v > 0 when it matters);
/// - same sign → that sign;
/// - opposite → `sign(P) · sign(P² − v·Q²)`, where `P² − vQ² ∈ ℚ(√u)` exactly:
///   `(A²+B²u − v(C²+D²u)) + (2AB − 2vCD)·√u` — one more [`sign1_int`] call.
pub fn biquad_sign(a: Rat, b: Rat, c: Rat, d: Rat, u: Rat, v: Rat) -> Option<Orient> {
    let zero = Rat::from_int(0);
    if u < zero || v < zero {
        // A domain answer, not a caller bug (the `QuadVal::new` contract): a negative
        // radicand means the value asked about does not exist.
        return None;
    }
    // Integerize both radicals: √(pu/qu) = √(pu·qu)/qu folds 1/qu into b and d; likewise v.
    // √(u·v) integerizes consistently: √(U·V)/(qu·qv) with U = pu·qu, V = pv·qv.
    let (up, uq) = (BigInt::from(u.numer()), BigInt::from(u.denom()));
    let (vp, vq) = (BigInt::from(v.numer()), BigInt::from(v.denom()));
    let (big_u, big_v) = (&up * &uq, &vp * &vq);
    // Coefficient denominators after folding: a/1, b/qu, c/qv, d/(qu·qv) — then clear the
    // rational denominators by the (positive) product of all four.
    let (an, ad) = (BigInt::from(a.numer()), BigInt::from(a.denom()));
    let (bn, bd) = (BigInt::from(b.numer()), BigInt::from(b.denom()) * &uq);
    let (cn, cd) = (BigInt::from(c.numer()), BigInt::from(c.denom()) * &vq);
    let (dn, dd) = (BigInt::from(d.numer()), BigInt::from(d.denom()) * &uq * &vq);
    let big_a = &an * &bd * &cd * &dd;
    let big_b = &bn * &ad * &cd * &dd;
    let big_c = &cn * &ad * &bd * &dd;
    let big_d = &dn * &ad * &bd * &cd;
    Some(biquad_sign_int(
        &big_a, &big_b, &big_c, &big_d, &big_u, &big_v,
    ))
}

fn biquad_sign_int(
    a: &BigInt,
    b: &BigInt,
    c: &BigInt,
    d: &BigInt,
    u: &BigInt,
    v: &BigInt,
) -> Orient {
    use num_bigint::Sign::NoSign;
    // √v degenerate: the element is P = A + B√u.
    if v.sign() == NoSign {
        return sign1_int(a, b, u);
    }
    // P = A + B√u, Q = C + D√u.
    let sp = sign1_int(a, b, u);
    let sq = sign1_int(c, d, u);
    if sq == Orient::Zero {
        return sp;
    }
    if sp == Orient::Zero {
        return sq;
    }
    if sp == sq {
        return sp;
    }
    // Opposite signs: compare P² against v·Q² inside ℚ(√u).
    // P² − vQ² = (A² + B²u − v(C² + D²u)) + (2AB − 2vCD)√u.
    let r_rat = a * a + b * b * u - v * (c * c + d * d * u);
    let r_rad = BigInt::from(2) * (a * b - v * c * d);
    orient_mul(sp, sign1_int(&r_rat, &r_rad, u))
}

// ──────────────────── plane · plane · cylinder (M6-1, commit 2) ────────────────────
//
// Pure numeric, handle-ignorant (the S5-prep option-(a) precedent): planes arrive as
// `[Rat; 4]` (n·x + d = 0), the cylinder as its raw rational fields. The producers in
// topo/ops translate their handles down to these values.

type V3 = [Rat; 3];

fn dot3(x: &V3, y: &V3) -> Option<Rat> {
    x[0].checked_mul(y[0])?
        .checked_add(x[1].checked_mul(y[1])?)?
        .checked_add(x[2].checked_mul(y[2])?)
}

fn cross3(x: &V3, y: &V3) -> Option<V3> {
    Some([
        x[1].checked_mul(y[2])?
            .checked_sub(x[2].checked_mul(y[1])?)?,
        x[2].checked_mul(y[0])?
            .checked_sub(x[0].checked_mul(y[2])?)?,
        x[0].checked_mul(y[1])?
            .checked_sub(x[1].checked_mul(y[0])?)?,
    ])
}

fn sub3(x: &V3, y: &V3) -> Option<V3> {
    Some([
        x[0].checked_sub(y[0])?,
        x[1].checked_sub(y[1])?,
        x[2].checked_sub(y[2])?,
    ])
}

fn is_zero3(x: &V3) -> bool {
    let zero = Rat::from_int(0);
    x.iter().all(|c| *c == zero)
}

/// Exact rational division by a **positive** divisor; `None` for a zero divisor or overflow.
///
/// ★ The positivity is load-bearing, not cosmetic: `Rat::new(denom, numer)` with a negative
/// `numer` normalizes the sign by negating both parts *inside* `Ratio`, which is unchecked —
/// at `numer = i128::MIN` that negation panics. Every call site divides by `|ℓ|²` or `2A`
/// (both positive), and the assert keeps the next caller honest.
fn rat_div(x: Rat, y: Rat) -> Option<Rat> {
    debug_assert!(y > Rat::from_int(0), "rat_div is for positive divisors");
    x.checked_mul(Rat::new(y.denom(), y.numer())?)
}

/// The rational line two intersecting planes share — the structure a
/// plane·plane·cylinder point rides ("point = line + parameter"): holding the point this way
/// makes "all three coordinates share one radical" structural instead of an invariant to
/// police, and lets `plane_side` collapse to one `QuadVal` multiply-add.
///
/// Only [`plane_plane_cylinder`] constructs one. `base` is derived deterministically (Cramer
/// against the auxiliary plane `ℓ·x = 0`, i.e. the point of the line closest to the origin),
/// so recomputing the meet reproduces the same parameters `s` bit-for-bit; the *order* of the
/// two roots along `dir` is base-independent either way.
#[derive(Clone, Debug)]
pub struct MeetLine {
    base: V3,
    dir: V3,
}

impl MeetLine {
    /// A point of the line, exact.
    pub fn base(&self) -> V3 {
        self.base
    }

    /// The line's direction `n₁ × n₂`, raw (unnormalized), exact.
    pub fn dir(&self) -> V3 {
        self.dir
    }

    /// The f64 realization of the point at parameter `s` — cache only.
    pub fn point_f64(&self, s: &QuadVal) -> [f64; 3] {
        let sf = s.to_f64();
        core::array::from_fn(|k| self.base[k].to_f64() + sf * self.dir[k].to_f64())
    }
}

/// What two planes and a cylinder share — every outcome has its own name, no silent
/// fallback. `Pair`'s parameters are in ascending order along [`MeetLine::dir`] (`A > 0` by
/// Cauchy–Schwarz once the degenerate arms are peeled off, so the ordering is well-defined).
#[derive(Clone, Debug)]
pub enum CylinderMeet {
    /// The planes are one plane (proportional 4-vectors) — their meet is a plane, not a line.
    /// Reachable: a nameless/`Through` population can hold one geometric plane as two handles.
    CoincidentPlanes,
    /// Parallel and distinct — no line at all.
    ParallelPlanes,
    /// The line is parallel to the axis and lies **on** the cylinder — a ruling; the meet is
    /// the whole line, not points.
    OnRuling(MeetLine),
    /// Parallel to the axis, off the surface (inside or outside) — no intersection.
    AxisParallelMiss(MeetLine),
    /// The line crosses the cylinder's ambient quadric nowhere (negative discriminant).
    Miss(MeetLine),
    /// Tangent — a double root, rational.
    Tangent { line: MeetLine, s: Rat },
    /// Two points, `s[0] < s[1]` along `dir`.
    Pair { line: MeetLine, s: [QuadVal; 2] },
}

/// The exact meet of two planes and a cylinder's lateral surface.
///
/// `dir`/`radius` are the cylinder's **raw** rational axis direction and radius
/// (`CylinderDef`'s fields); preconditions `dir ≠ 0`, `radius > 0` are the caller's
/// (doc + debug_assert — this pure-numeric door is outside `CylinderDef::new`'s guard).
///
/// `None` means checked-`Rat` overflow — an honest decline, never a wrong variant.
pub fn plane_plane_cylinder(
    p1: &[Rat; 4],
    p2: &[Rat; 4],
    origin: &V3,
    dir: &V3,
    radius: Rat,
) -> Option<CylinderMeet> {
    let zero = Rat::from_int(0);
    debug_assert!(!is_zero3(dir), "cylinder axis must be nonzero");
    debug_assert!(radius > zero, "cylinder radius must be positive");
    let n1 = [p1[0], p1[1], p1[2]];
    let n2 = [p2[0], p2[1], p2[2]];
    let l = cross3(&n1, &n2)?;
    if is_zero3(&l) {
        // Parallel normals: coincident iff the full 4-vectors are proportional. Cross-multiply
        // against a nonzero normal component (both normals are nonzero — a plane's contract).
        let i = (0..3)
            .find(|&k| n1[k] != zero)
            .expect("a plane has a nonzero normal");
        let coincident = p2[3].checked_mul(n1[i])? == p1[3].checked_mul(n2[i])?;
        return Some(if coincident {
            CylinderMeet::CoincidentPlanes
        } else {
            CylinderMeet::ParallelPlanes
        });
    }
    // Deterministic base: Cramer for { n₁·x = −d₁, n₂·x = −d₂, ℓ·x = 0 }. The system's
    // determinant is ℓ·(n₁×n₂) = |ℓ|² > 0.
    let det = dot3(&l, &l)?;
    let rhs = [zero.checked_sub(p1[3])?, zero.checked_sub(p2[3])?, zero];
    let col = |k: usize| -> [Rat; 3] { [n1[k], n2[k], l[k]] };
    let det3 = |c0: &V3, c1: &V3, c2: &V3| -> Option<Rat> { dot3(c0, &cross3(c1, c2)?) };
    let mut base = [zero; 3];
    for (k, out) in base.iter_mut().enumerate() {
        // Replace column k by the rhs (columns here are the rows' k-th entries).
        let cols: [V3; 3] = core::array::from_fn(|j| if j == k { rhs } else { col(j) });
        *out = rat_div(det3(&cols[0], &cols[1], &cols[2])?, det)?;
    }
    let line = MeetLine { base, dir: l };

    // The quadratic A·s² + B·s + C = 0 of the line against |w|²|m|² − (w·m)² = r²|m|².
    let m = dir;
    let w0 = sub3(&base, origin)?;
    let mm = dot3(m, m)?;
    let ll = dot3(&line.dir, &line.dir)?;
    let lm = dot3(&line.dir, m)?;
    let w0m = dot3(&w0, m)?;
    let w0l = dot3(&w0, &line.dir)?;
    let w0w0 = dot3(&w0, &w0)?;
    let a = ll.checked_mul(mm)?.checked_sub(lm.checked_mul(lm)?)?;
    let b =
        Rat::from_int(2).checked_mul(w0l.checked_mul(mm)?.checked_sub(w0m.checked_mul(lm)?)?)?;
    let c = w0w0
        .checked_mul(mm)?
        .checked_sub(w0m.checked_mul(w0m)?)?
        .checked_sub(radius.checked_mul(radius)?.checked_mul(mm)?)?;

    if a == zero {
        // ℓ ∥ axis: substituting ℓ = k·m makes B vanish identically — checked, not assumed.
        debug_assert!(
            b == zero,
            "A = 0 forces B = 0 for a line parallel to the axis"
        );
        return Some(if c == zero {
            CylinderMeet::OnRuling(line)
        } else {
            CylinderMeet::AxisParallelMiss(line)
        });
    }
    debug_assert!(
        a > zero,
        "A = |ℓ|²|m|² − (ℓ·m)² is nonnegative (Cauchy–Schwarz)"
    );
    let disc = b
        .checked_mul(b)?
        .checked_sub(Rat::from_int(4).checked_mul(a)?.checked_mul(c)?)?;
    if disc < zero {
        return Some(CylinderMeet::Miss(line));
    }
    let two_a = Rat::from_int(2).checked_mul(a)?;
    let mid = rat_div(zero.checked_sub(b)?, two_a)?;
    if disc == zero {
        return Some(CylinderMeet::Tangent { line, s: mid });
    }
    let half = rat_div(Rat::from_int(1), two_a)?;
    let lo = QuadVal::new(mid, zero.checked_sub(half)?, disc)?;
    let hi = QuadVal::new(mid, half, disc)?;
    Some(CylinderMeet::Pair { line, s: [lo, hi] })
}

/// The side of `plane` a line-point lies on — `n·p + d = (n·base + d) + s·(n·dir)`, one
/// `QuadVal` multiply-add, then the sign tower. `None` on overflow.
pub fn plane_side(plane: &[Rat; 4], line: &MeetLine, s: &QuadVal) -> Option<Orient> {
    let n = [plane[0], plane[1], plane[2]];
    let base_term = dot3(&n, &line.base)?.checked_add(plane[3])?;
    let dir_term = dot3(&n, &line.dir)?;
    let val = QuadVal::from_rat(base_term).checked_add(&s.checked_mul_rat(dir_term)?)?;
    Some(val.sign())
}

/// Where a point sits on the seam-cut circle `(0, 2π)` — the four classes the circular order
/// speaks in. `θ` is measured right-handed about the axis from the seam (`+e₁`, the
/// axis-perpendicular part of `ref_dir`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SeamClass {
    /// On the seam generator itself (θ = 0) — **outside** the chart, surfaced by name.
    Seam,
    /// θ ∈ (0, π): the `w·e₂ > 0` half.
    Upper,
    /// θ = π exactly.
    Pi,
    /// θ ∈ (π, 2π).
    Lower,
}

/// The circular order of two cylinder-surface points about the seam — the comparator M6-2's
/// sweep will sort with.
///
/// ★ **The API shape is a candidate** (the math is not): the return type and calling
/// convention are the sweep's to finalize; what is fixed here is the class ladder and the
/// sign conventions, property-locked below.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeamOrder {
    /// At least one point lies on the seam generator (θ = 0) — outside the chart's total
    /// order, answered by name so a sweep cannot silently rank the excluded point.
    SeamIncident { first: bool, second: bool },
    /// Both points are chart-interior; their order along θ ∈ (0, 2π).
    Ordered(core::cmp::Ordering),
}

/// Compare two points of the cylinder's surface by their angle θ ∈ (0, 2π) about the seam.
///
/// Each point is a line-parameter pair from [`plane_plane_cylinder`] against the **same**
/// cylinder (`origin`/`dir`/`ref_dir` raw rational — the caller's `CylinderDef` fields;
/// points on the surface are a precondition, debug-asserted).
///
/// The two predicates (no chart-`t` arithmetic — its value drags norm radicals in):
/// - class: `sign(w·e₂)` splits upper/lower; on the boundary, `sign(w·e₁)` splits seam / π,
///   where `e₁ = (m·m)·ref − (ref·m)·m` (the axis-perpendicular part of `ref_dir`, scaled
///   positive-rationally) and `e₂ = m × e₁`;
/// - within one open half: `sign((w₁×w₂)·m)` — positive ⇔ θ₁ < θ₂ (right-handed, |Δθ| < π).
///
/// `None` on checked-`Rat` overflow (honest decline).
pub fn circular_order_about_seam(
    origin: &V3,
    dir: &V3,
    ref_dir: &V3,
    first: (&MeetLine, &QuadVal),
    second: (&MeetLine, &QuadVal),
) -> Option<SeamOrder> {
    let m = dir;
    // e₁ = (m·m)·ref − (ref·m)·m: rational, ⊥ m, positively proportional to the unit e₁.
    let mm = dot3(m, m)?;
    let rm = dot3(ref_dir, m)?;
    let mut e1 = [Rat::from_int(0); 3];
    for (k, out) in e1.iter_mut().enumerate() {
        *out = mm
            .checked_mul(ref_dir[k])?
            .checked_sub(rm.checked_mul(m[k])?)?;
    }
    // `ref_dir ∥ axis` pins no seam — `CylinderDef::new` forbids it, but this pure-numeric
    // door sits outside that guard, and without the check it would fall through to the
    // off-surface assert with the wrong diagnosis.
    if is_zero3(&e1) {
        debug_assert!(false, "ref_dir parallel to the axis pins no seam");
        return None;
    }
    let e2 = cross3(m, &e1)?;

    // w(s)·e as a QuadVal: (u·e) + s·(ℓ·e), u = base − origin.
    let dot_w = |(line, s): (&MeetLine, &QuadVal), e: &V3| -> Option<QuadVal> {
        let u = sub3(&line.base, origin)?;
        QuadVal::from_rat(dot3(&u, e)?).checked_add(&s.checked_mul_rat(dot3(&line.dir, e)?)?)
    };
    let class = |p: (&MeetLine, &QuadVal)| -> Option<SeamClass> {
        Some(match dot_w(p, &e2)?.sign() {
            Orient::Positive => SeamClass::Upper,
            Orient::Negative => SeamClass::Lower,
            Orient::Zero => match dot_w(p, &e1)?.sign() {
                Orient::Positive => SeamClass::Seam,
                Orient::Negative => SeamClass::Pi,
                Orient::Zero => {
                    // w ⊥-part vanishes: the point is on the axis, not the surface.
                    debug_assert!(false, "circular order asked about a point off the surface");
                    return None;
                }
            },
        })
    };
    let (c1, c2) = (class(first)?, class(second)?);
    if c1 == SeamClass::Seam || c2 == SeamClass::Seam {
        return Some(SeamOrder::SeamIncident {
            first: c1 == SeamClass::Seam,
            second: c2 == SeamClass::Seam,
        });
    }
    use core::cmp::Ordering;
    let rank = |c: SeamClass| match c {
        SeamClass::Upper => 0u8,
        SeamClass::Pi => 1,
        SeamClass::Lower => 2,
        SeamClass::Seam => unreachable!("surfaced above"),
    };
    if c1 != c2 {
        return Some(SeamOrder::Ordered(rank(c1).cmp(&rank(c2))));
    }
    if c1 == SeamClass::Pi {
        return Some(SeamOrder::Ordered(Ordering::Equal));
    }
    // Same open half: sign((w₁×w₂)·m) — the axial parts cancel algebraically, and within one
    // half |Δθ| < π, so this sign is the angle order. With wᵢ = uᵢ + sᵢℓᵢ the triple product
    // expands to k₀ + k₁s₁ + k₂s₂ + k₃s₁s₂ (rational kⱼ), and s₁s₂ opens into the 4-term
    // ℚ(√c₁,√c₂) element the biquadratic sign closes.
    let (l1, s1) = first;
    let (l2, s2) = second;
    let u1 = sub3(&l1.base, origin)?;
    let u2 = sub3(&l2.base, origin)?;
    let trip = |x: &V3, y: &V3| -> Option<Rat> { dot3(&cross3(x, y)?, m) };
    let k0 = trip(&u1, &u2)?;
    let k1 = trip(&l1.dir, &u2)?;
    let k2 = trip(&u1, &l2.dir)?;
    let k3 = trip(&l1.dir, &l2.dir)?;
    let (a1, b1, rad1) = (s1.a(), s1.b(), s1.c());
    let (a2, b2, rad2) = (s2.a(), s2.b(), s2.c());
    let big_a = k0
        .checked_add(k1.checked_mul(a1)?)?
        .checked_add(k2.checked_mul(a2)?)?
        .checked_add(k3.checked_mul(a1)?.checked_mul(a2)?)?;
    let big_b = b1.checked_mul(k1.checked_add(k3.checked_mul(a2)?)?)?;
    let big_c = b2.checked_mul(k2.checked_add(k3.checked_mul(a1)?)?)?;
    let big_d = k3.checked_mul(b1)?.checked_mul(b2)?;
    let cross_sign = biquad_sign(big_a, big_b, big_c, big_d, rad1, rad2)?;
    Some(SeamOrder::Ordered(match cross_sign {
        Orient::Positive => Ordering::Less,
        Orient::Negative => Ordering::Greater,
        Orient::Zero => Ordering::Equal,
    }))
}

#[cfg(test)]
#[path = "quad_tests.rs"]
mod tests;
