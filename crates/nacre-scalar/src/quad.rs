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

    /// Exact negation (cannot overflow: `Rat` keeps the sign in an `i128` numerator whose
    /// negation only fails at `i128::MIN`, which a reduced `Ratio` never holds).
    pub fn neg(&self) -> QuadVal {
        let zero = Rat::from_int(0);
        QuadVal {
            a: zero
                .checked_sub(self.a)
                .expect("negating a reduced rational"),
            b: zero
                .checked_sub(self.b)
                .expect("negating a reduced rational"),
            c: self.c,
        }
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

#[cfg(test)]
#[path = "quad_tests.rs"]
mod tests;
