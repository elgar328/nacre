use super::*;
/// A rational base point carried through a chain of axis rotations (a rotation history). `base` +
/// `chain` are the exact **definition** (never lost); `coord_iv` is the f64 coordinate (a cache)
/// with, per axis, the bound on its error — a **direction-wise xyz vector** held beside the value
/// it bounds, so the two cannot drift apart. [`coord_hp`](Self::coord_hp) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,
    // (equality is definitional — see the `PartialEq` impl below the struct)
    pub coord_iv: [Bounded; 3],
    /// Memoized `coord_hp` at the boolean's chosen precision — the astro-float realization
    /// once per definition and shared across clones (`Rc`). A judge escalates the *same*
    /// definition-point dozens of times per boolean (`plane_def` clones `tri_pt3` per call);
    /// without this each escalation replays the rotation's cos/sin at 200 bits, which dominates
    /// the rotated-boolean cost. `base`/`chain` never change after construction except through
    /// [`Self::rotate_about`], which resets this cell, so the cached value always matches the
    /// definition (a pure, path-independent function). See [`HpCell`] for the
    /// `Arc<OnceLock>` (parallel) vs `Rc<OnceCell>` (serial) choice.
    hp: HpCell,
}

/// **Definitional equality — `base` and `chain`, nothing else.** `coord_iv`/`hp` are caches,
/// and realization is a pure function of the definition (path independence is this file's root
/// soundness argument), so two points with equal definitions cannot honestly disagree in their
/// caches. The consumer this exists for is [`shared_base`]'s whole-node chain comparison — it
/// declares two motions to be *one* motion, a question about definitions, never about caches —
/// which [`FrameThrough`] joins by carrying points inside a node.
impl PartialEq for WitnessPoint {
    fn eq(&self, o: &Self) -> bool {
        self.base == o.base && *self.chain == *o.chain
    }
}
impl Eq for WitnessPoint {}

// ── The steps' arithmetic: `Bounded`, and the cases it cannot see are exact ──
//
// Every step below carries its coordinates through `Bounded`'s operations, whose radius is a proved
// bound per operation (value rounding, the radius's own rounding, the operands' radii — including
// what a realized `cos` or `1/√` contributes through the incoming tol). `Bounded` charges every
// result for its rounding, though, and some steps' operations do not round: adding an exact zero,
// or multiplying by an exact `0` or power of two (`±1`, a frame's `2⁻ᵏ`). Charging those would
// leave no model at tol `0` — a quadrantal turn about the origin, a reflection in a plane through
// it, an axis-aligned frame — and tol `0` is what keeps such a model on the exact predicate path.
// So these three ask first; the facts hold in `f64`'s normal range, the premise `Bounded` rests on.
//
// ★ **Here, not inside `Bounded`'s own operations.** The judge's interval filters meet the same
// cases often — 42–51% of their additions have an exact-zero operand and 31–40% of their products
// an exact `0` or power of two (census, perf) — but answering them there decides nothing: built
// both ways, the climbs and census rows are identical, and the check costs every filter operation
// (+8.2% instructions on the perf folds returning early, +16.9% choosing the charge). The fact
// matters only where tol `0` is a contract, which is here.

/// `x + y`, exact when either is an exact zero.
fn sum(x: Bounded, y: Bounded) -> Bounded {
    match (exact_zero(x), exact_zero(y)) {
        (_, true) => Bounded::new(x.value + y.value, x.error),
        (true, _) => Bounded::new(x.value + y.value, y.error),
        _ => x.add(y),
    }
}

/// `x − y`, exact when either is an exact zero.
fn diff(x: Bounded, y: Bounded) -> Bounded {
    match (exact_zero(x), exact_zero(y)) {
        (_, true) => Bounded::new(x.value - y.value, x.error),
        (true, _) => Bounded::new(x.value - y.value, y.error),
        _ => x.sub(y),
    }
}

/// `x · y`, exact when either is an exact `0` or power of two — then the other's radius scales by
/// the same power, exactly.
fn prod(x: Bounded, y: Bounded) -> Bounded {
    let scaling = |b: Bounded| {
        b.error == 0.0
            && (b.value == 0.0 || (b.value.is_normal() && b.value.to_bits() & ((1 << 52) - 1) == 0))
    };
    if scaling(y) {
        Bounded::new(x.value * y.value, x.error * y.value.abs())
    } else if scaling(x) {
        Bounded::new(x.value * y.value, y.error * x.value.abs())
    } else {
        x.mul(y)
    }
}

fn exact_zero(b: Bounded) -> bool {
    b.value == 0.0 && b.error == 0.0
}

/// A point read as `(p₀, p₁, p₂)` in a frame and written in the world: `o + ((p₀·û + p₁·v̂) + p₂·ŵ)`
/// per axis, in that order, with the basis rows `[û, v̂, ŵ]` and the origin `o` as intervals —
/// the three frame nodes' one combination. The incoming tol rides in `p` and is turned by the basis
/// through `prod`; the basis's own error is scaled by the coordinates the same way.
fn in_frame(p: [Bounded; 3], [u, v, w]: [[Bounded; 3]; 3], o: [Bounded; 3]) -> [Bounded; 3] {
    core::array::from_fn(|k| {
        let terms = sum(sum(prod(p[0], u[k]), prod(p[1], v[k])), prod(p[2], w[k]));
        sum(o[k], terms)
    })
}

impl WitnessPoint {
    /// The `f64` coordinate — the values of [`Self::coord_iv`], for readers that want the point
    /// whole. Approximate: a rotation is carried in `f64` arithmetic, so the value is not the
    /// nearest `f64` of the point; the proven distance from it is [`Self::tol`].
    #[inline]
    pub fn coord_approx(&self) -> [f64; 3] {
        self.coord_iv.map(|b| b.value)
    }

    /// The per-axis error bounds — the radii of [`Self::coord_iv`], for readers that want them
    /// as a vector.
    #[inline]
    pub fn tol(&self) -> [f64; 3] {
        self.coord_iv.map(|b| b.error)
    }

    /// A point at a rational `base`: each coordinate its nearest `f64`, with the error that
    /// rounding carries ([`Bounded::of_rat`]) — `0` exactly where the coordinate is an `f64`,
    /// else half an ulp, **per axis**: an axis a later chain never moves keeps exactly its own,
    /// which a per-axis tol check needs.
    ///
    /// The bound is `Rat::to_f64`'s contract (nearest, ties to even), not a measurement — a
    /// measurement at a fixed precision reads a rational nearer its `f64` than that precision
    /// resolves as exact, and costs BigFloat work per point. Its looseness can only turn a
    /// definite filter answer into "escalate", never the other way round.
    ///
    /// ★ **This is how a point that came from a rational reaches the predicates** — not by
    /// lifting the rounded `f64` back into a `Rat` (that names a different point, with a
    /// tolerance of zero on top).
    pub fn at(base: [Rat; 3]) -> Self {
        WitnessPoint {
            coord_iv: base.map(Bounded::of_rat),
            base,
            chain: HpRc::from([] as [MoveNode; 0]),
            hp: HpCell::default(),
        }
    }

    /// Extend the definition by a rotation about `axis` through the origin. See
    /// [`rotate_about`](Self::rotate_about).
    pub fn rotate(self, axis: Axis, angle: Angle) -> Self {
        self.rotate_about(axis, angle, [Rat::from_int(0); 3])
    }

    /// Extend the definition by a rotation about `axis` through the rational pivot
    /// `point`, updating the f64 cache and its direction-wise tol. The two in-plane
    /// coords are taken relative to the pivot, rotated, and shifted back. Same-axis
    /// 90°-family angles about the origin rotate exactly (tol 0, Niven).
    ///
    /// **The tol is `Bounded`'s, step by step** — `cos`/`sin` enter as the correctly rounded pair
    /// `Angle::cos_sin_f64` gives, each with the error `Angle::cos_sin_error_of` *measured* of it,
    /// and the pivot as `Bounded::of_rat`. Each axis is its own expression (`i` mixes `u·c − v·s`,
    /// `j` mixes `u·s + v·c`), so `u` meets `cos`'s error on one and `sin`'s on the other; the
    /// intervals keep that apart without being told (on a 356.65° turn `sin` is eight times
    /// further out than `cos`).
    pub fn rotate_about(self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let mut p = self.rotated(axis, angle, point);
        p.remember(MoveNode::Rotate {
            axis,
            angle,
            pivot: point,
        });
        p
    }

    /// [`Self::rotate_about`]'s numbers, without remembering the node — the step
    /// [`Self::apply_chain`] folds.
    fn rotated(mut self, axis: Axis, angle: Angle, point: [Rat; 3]) -> Self {
        let (i, j) = axis.plane();
        // **Through `Angle`'s single entry point, not re-spelled here.** At the 90°-family it
        // returns exactly `0.0`/`±1.0`, and the measured error is exactly `0` — the exactness
        // `prod` reads; a second spelling of "how an angle becomes f64" could land
        // `cos(90°) ≈ 6e-17` there instead.
        //
        // ★★★ **The error is measured of `c` and `s` — the values that reach `coord_approx` — not
        // of the angle**, so it bounds the number this kernel stores, whatever route produced it.
        // Measured here rather than taken as a constant: a platform whose trig lands further out
        // reports a bigger number.
        let (c, s) = angle.cos_sin_f64();
        let (dc, ds) = angle.cos_sin_error_of(c, s);
        let (c, s) = (Bounded::new(c, dc), Bounded::new(s, ds));
        let (pi, pj) = (Bounded::of_rat(point[i]), Bounded::of_rat(point[j]));
        // About the pivot and back: `(p + u·c) − v·s`, `(q + u·s) + v·c`. About the origin the
        // pivot is an exact zero and `sum`/`diff` pass the radius through — but the value still
        // takes `0.0 + …`, which is the value this point has always had.
        let (u, v) = (diff(self.coord_iv[i], pi), diff(self.coord_iv[j], pj));
        self.coord_iv[i] = diff(sum(pi, prod(u, c)), prod(v, s));
        self.coord_iv[j] = sum(sum(pj, prod(u, s)), prod(v, c));
        self
    }

    /// **Append one node to the remembered definition.**
    ///
    /// The slice is shared (`HpRc`) so a *clone* is a refcount bump rather than an allocation —
    /// which the judgment path does hundreds of thousands of times. The price is here: appending
    /// copies the slice, so applying `n` nodes one at a time costs `O(n²)`.
    ///
    /// ⚠★★★★ **That price was measured against the wrong caller, and the comment said so**: *"this
    /// runs when a solid is transformed, never on the judgment path, so the copy is not hot."* Cell
    /// 52 broke the premise — a vertex is now realized as it is pushed, which **replays a chain per
    /// push** — and the copy became the dominant cost of a long history: a 4,200-turn fixture went
    /// from 5.6 s to 563 s, and a 1,600-turn loop from 17 ms to 114 s. A caller holding a whole
    /// chain must therefore use [`Self::apply_chain`], which folds the numbers and appends **once**.
    fn remember(&mut self, node: MoveNode) {
        let mut nodes = self.chain.to_vec();
        nodes.push(node);
        self.chain = HpRc::from(nodes);
        // The definition changed — invalidate the memo of the old one. A fresh (unshared) cell, so
        // clones made before this keep their own cached value.
        self.hp = HpCell::default();
    }

    /// **The whole chain applied in one pass** — the numbers folded node by node, the definition
    /// remembered once. Equivalent to applying each node through its own method, and measured so:
    /// same coordinate, same tol, same high-precision realization, bit for bit.
    ///
    /// This is what a *replay* wants. Rebuilding the remembered chain per node (what the public
    /// per-node methods must do, since each is a transform in its own right) makes a replay
    /// quadratic in the history's length — see [`Self::remember`].
    pub fn apply_chain(mut self, chain: &[MoveNode]) -> Option<Self> {
        for n in chain {
            self = match n {
                MoveNode::Rotate { axis, angle, pivot } => self.rotated(*axis, *angle, *pivot),
                MoveNode::Translate { offset } => self.translated(*offset),
                MoveNode::Mirror { axis, offset } => self.mirrored(*axis, *offset),
                MoveNode::Frame { frame } => self.framed(*frame)?,
                MoveNode::FrameWide(f) => self.framed_wide(f)?,
                MoveNode::FrameThrough(f) => self.framed_through(f)?,
            };
        }
        let mut nodes = self.chain.to_vec();
        nodes.extend(chain.iter().cloned());
        self.chain = HpRc::from(nodes);
        self.hp = HpCell::default();
        Some(self)
    }

    /// The point reflected in `axis = offset` — one more link in the definition.
    ///
    /// **`coord_approx` is this step's `f64`** (`2c − x`) and `tol` is `Bounded`'s for it, so
    /// `coord_iv` holds the definition's value. The exact offset lives in the chain, where
    /// [`coord_hp_uncached`](Self::coord_hp_uncached) realizes it.
    pub fn mirror(self, axis: Axis, offset: Rat) -> Self {
        let mut p = self.mirrored(axis, offset);
        p.remember(MoveNode::Mirror { axis, offset });
        p
    }

    /// [`Self::mirror`]'s numbers, without remembering the node.
    fn mirrored(mut self, axis: Axis, offset: Rat) -> Self {
        let k = axis.index();
        // `2·c` is a power-of-two multiply — exact, its radius doubled with it. A reflection in a
        // plane through the origin is `0 − x`, a negation, exact as well.
        let twice = prod(Bounded::new(2.0, 0.0), Bounded::of_rat(offset));
        self.coord_iv[k] = diff(twice, self.coord_iv[k]);
        self
    }

    /// The point translated by an exact rational `offset` — one more link in the definition.
    ///
    /// **`coord_approx` adds the offset's `f64` image** and `tol` is `Bounded`'s for the sum (the
    /// offset's own rounding included), so `coord_iv` holds the definition's value. The exact offset lives in the chain, where
    /// [`coord_hp_uncached`](Self::coord_hp_uncached) realizes it.
    ///
    /// That split is the whole point: two placements that reach the same real wall by different
    /// routes keep f64 coordinates an ulp apart, but their *definitions* realize to the same
    /// value, and it is the definition the judgment reads.
    pub fn translate(self, offset: [Rat; 3]) -> Self {
        let mut p = self.translated(offset);
        p.remember(MoveNode::Translate { offset });
        p
    }

    /// [`Self::translate`]'s numbers, without remembering the node.
    fn translated(mut self, offset: [Rat; 3]) -> Self {
        for (x, &off) in self.coord_iv.iter_mut().zip(&offset) {
            *x = sum(*x, Bounded::of_rat(off));
        }
        self
    }

    /// The point read as `(u, v, w)` in the frame of the plane `origin`/`u_raw`/`n` describe, and
    /// written out in that plane's own frame — one more link in the definition.
    ///
    /// **`coord_approx` is this step's `f64`** and `tol` is `Bounded`'s for it, so `coord_iv` holds
    /// the definition's value. The exact inputs live in the chain, where
    /// [`coord_hp_uncached`](Self::coord_hp_uncached) realizes them.
    ///
    /// `None` when a squared length is not positive — that is a degenerate frame, not an unusual
    /// one, and a caller that cannot build a basis must not get a point that pretends otherwise.
    ///
    /// ★★★★ **Every charge has a source, and the exact cases really do reach zero.** The sources
    /// are the rational inputs' own `Rat → f64` rounding ([`Bounded::of_rat`] — `0` for an input
    /// that is an `f64`), the `1/√` values' realization error, measured
    /// ([`nacre_exact::inv_sqrt_error_of`] — *not* a half-ulp taken on faith, see there), and
    /// `Bounded`'s charge for the arithmetic that combines them. An axis-aligned frame through the
    /// origin has a signed permutation for a basis — integer components over an exact `1/√` of a
    /// power of four, so each component is an exact `0`/`±1` — and its products and single-term
    /// sums are exact, so the incoming tol passes through untouched.
    pub fn frame(self, f: nacre_exact::PlaneFrame) -> Option<Self> {
        let mut p = self.framed(f)?;
        p.remember(MoveNode::Frame { frame: f });
        Some(p)
    }

    /// [`Self::frame`]'s numbers, without remembering the node.
    fn framed(mut self, f: nacre_exact::PlaneFrame) -> Option<Self> {
        // One inverse square root per axis, each of an exact rational — see `plane_frame` for why
        // `v̂` gets its own instead of being a cross product of the other two. A basis component
        // is the rational over that length: `of_rat`'s rounding scaled by it, its measured error
        // scaled by the rational, and the product's rounding — `prod` sees the exact cases.
        let axis = |raw: [Rat; 3], squared: Rat| -> Option<[Bounded; 3]> {
            let x = nacre_exact::inv_sqrt_f64(squared)?;
            let len = Bounded::new(x, nacre_exact::inv_sqrt_error_of(squared, x)?);
            Some(raw.map(|r| prod(Bounded::of_rat(r), len)))
        };
        let (u, w) = (axis(f.u_raw, f.uu)?, axis(f.n, f.nn)?);
        // ★★ `v̂` exactly where the frame carries `v_raw`, and `ŵ × û` where it does not.
        //
        // The exact route is one rounding per component; the cross product is two that do not
        // cancel — a wall whose `v` is exactly `ẑ` comes out three ulps short through it. Both are
        // sound; the fallback is taken only when `|n|²·|u_raw|²` does not fit `i128` (see
        // `nacre_exact::plane_frame`), which no axis-aligned frame reaches.
        let v = match f.v {
            Some((v_raw, vv)) => axis(v_raw, vv)?,
            None => core::array::from_fn(|k| {
                let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                diff(prod(w[i], u[j]), prod(w[j], u[i]))
            }),
        };
        self.coord_iv = in_frame(self.coord_iv, [u, v, w], f.origin.map(Bounded::of_rat));
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`WideFrame`] — the same combination (`in_frame`), with the
    /// axis and origin intervals taken from a fixed-precision arbitrary-precision realization
    /// instead of `inv_sqrt_f64` over a rational. Every rounding on the way to them is inside an
    /// `HpBounded`, and the narrowing to `f64` charges itself.
    ///
    /// ★ The axes are realized per applied point. The wide population is a fraction of a percent
    /// of pushes, so that cost is accepted rather than memoized — a measure-later item.
    ///
    /// ★ A wide frame is never an axis permutation (its squared lengths exceed `i128`), so its
    /// components are never the exact `0`/`±1` that keep an axis-aligned `frame` at tol 0.
    pub fn frame_wide(self, f: &WideFrame) -> Option<Self> {
        let mut p = self.framed_wide(f)?;
        p.remember(MoveNode::FrameWide(f.clone()));
        Some(p)
    }

    /// [`Self::frame_wide`]'s numbers, without remembering the node.
    fn framed_wide(mut self, f: &WideFrame) -> Option<Self> {
        // The fixed rung for the f64 cache of a wide frame — the ladder's first rung, the same
        // one `inv_sqrt_f64` starts at. The judgment path re-realizes at its own precision.
        const P: usize = 128;
        let inv = |v: &num_bigint::BigInt| nacre_exact::inv_sqrt_hp_big(v, P);
        let axis = |raw: &[num_bigint::BigInt; 3], squared: &num_bigint::BigInt| {
            let len = inv(squared)?;
            Some(
                raw.each_ref()
                    .map(|r| narrow_hp(&HpBounded::of_bigint(r, P).mul(&len, P))),
            )
        };
        let (u, v, w) = (
            axis(&f.u_raw, &f.uu)?,
            axis(&f.v_raw, &f.vv)?,
            axis(&f.n, &f.nn)?,
        );
        let od = HpBounded::of_bigint(&f.origin_den, P);
        let mut o = [Bounded::new(0.0, 0.0); 3];
        for (o, num) in o.iter_mut().zip(&f.origin_num) {
            *o = narrow_hp(&HpBounded::of_bigint(num, P).div_exact(&od, P)?);
        }
        self.coord_iv = in_frame(self.coord_iv, [u, v, w], o);
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`FrameThrough`] — the judged basis realized
    /// at the fixed rung, narrowed to `f64` intervals, then **the same combination as
    /// [`WitnessPoint::frame_wide`]** (`in_frame`).
    ///
    /// `None` is unreachable for a node built by [`FrameThrough::of`] — it proved this exact
    /// derivation at this exact rung — and stays an `Option` so that "unreachable" is a fact
    /// about the producer rather than an invariant this method asserts across a crate boundary.
    ///
    /// ★ Cost note: the node's points share their `hp` cells by `Rc`, so the cos/sin of their
    /// chains realize once per node, not once per applied point — what recurs per point is
    /// arithmetic, the same acceptance [`WitnessPoint::frame_wide`] records.
    pub fn frame_through(self, f: &FrameThrough) -> Option<Self> {
        let mut p = self.framed_through(f)?;
        p.remember(MoveNode::FrameThrough(Box::new(f.clone())));
        Some(p)
    }

    /// [`Self::frame_through`]'s numbers, without remembering the node.
    fn framed_through(mut self, f: &FrameThrough) -> Option<Self> {
        let basis = judged_basis(f, FrameThrough::RUNG)?;
        let [o, u, v, w] = basis.map(|row| row.map(|c| narrow_hp(&c)));
        self.coord_iv = in_frame(self.coord_iv, [u, v, w], o);
        Some(self)
    }

    /// **The coordinate realized at `prec` bits from the definition** (base rotated through the
    /// chain, each node about its pivot), **with the error it carries** — path-independent
    /// ground truth / escalation realization, in the type this workspace hands out publicly
    /// (`nacre_exact::HpBounded`, what `inv_sqrt_hp` returns too). The result is memoized in
    /// [`WitnessPoint::hp`] and shared across clones of the same definition, so a definition-point
    /// pays the astro-float cos/sin once per boolean rather than once per predicate.
    ///
    /// ★ **The radius is half the answer.** A value without the bound its realization cost cannot
    /// be rounded honestly — see [`HpBounded`]'s contract: a radius invented for convenience makes
    /// every sign above it unearned. Callers round with `nacre_exact::round_to_f64` /
    /// `round_to_digits`, which report *undecided* rather than picking when `prec` is short.
    pub fn coord_hp(&self, prec: usize) -> [HpBounded; 3] {
        // **Keyed by precision, and that key is load-bearing.** The precision is chosen per
        // boolean, so within one operation every call arrives with the same value and the cell is
        // filled once — which is the whole point, since re-realizing a definition per predicate
        // was 77% of a rotated boolean's runtime. A call at a different precision (the rare
        // per-judgement fallback, or a caller outside the boolean) recomputes without disturbing
        // the cached value, rather than silently returning coordinates realized at the wrong
        // precision.
        let cached = self.hp.get_or_init(|| (prec, self.coord_hp_uncached(prec)));
        if cached.0 == prec {
            return cached.1.clone();
        }
        self.coord_hp_uncached(prec)
    }

    /// [`Self::coord_hp`] without the memo, **with the error it carries**.
    ///
    /// This is the walk the `f64` steps take, one level up — `HpBounded` where they carry
    /// `Bounded`, one machine at two precisions: the definition is exact, the realization is not,
    /// and the radius is what the realization cost. Nothing here is a chosen constant — the base contributes its own rounding (zero when
    /// the rational lands on a `prec`-bit dyadic), each `cos`/`sin` contributes the bound
    /// [`Angle::cos_sin_hp`] derives, and every arithmetic operation adds its half-ulp.
    pub(super) fn coord_hp_uncached(&self, prec: usize) -> [HpBounded; 3] {
        let p = [
            HpBounded::of_rat(self.base[0], prec),
            HpBounded::of_rat(self.base[1], prec),
            HpBounded::of_rat(self.base[2], prec),
        ];
        fold_suffix(p, &self.chain, prec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The steps' arithmetic passes a radius through only where the operation is exact, and
    /// charges where it rounds.** The chain locks see the first half (an axis-aligned model stays
    /// at tol `0`) but not the second: a product left uncharged hides inside a measured `1/√`
    /// error, which reads up to twice its truth — so the rule is held here, on exact operands
    /// whose results round, against the exact rational answer.
    #[test]
    fn the_exact_cases_pass_through_and_the_rest_are_charged() {
        let r = |x: f64| Rat::try_from_f64(x).expect("finite");
        let b = |x: f64| Bounded::new(x, 0.0);
        let covers = |got: Bounded, want: Rat, what: &str| {
            let gap = r(got.value).checked_sub(want).unwrap();
            let zero = Rat::from_int(0);
            let gap = if gap < zero {
                zero.checked_sub(gap).unwrap()
            } else {
                gap
            };
            assert!(gap > zero, "{what}: the fixture must round");
            assert!(
                r(got.error) >= gap,
                "{what}: radius {:e} under the rounding",
                got.error
            );
        };
        covers(
            prod(b(3.0), b(0.1)),
            r(3.0).checked_mul(r(0.1)).unwrap(),
            "3·0.1",
        );
        covers(
            sum(b(0.1), b(0.2)),
            r(0.1).checked_add(r(0.2)).unwrap(),
            "0.1 + 0.2",
        );
        covers(
            diff(b(1.0), b(0.1)),
            r(1.0).checked_sub(r(0.1)).unwrap(),
            "1 − 0.1",
        );

        let t = Bounded::new(0.1, 1e-20);
        let same = |got: Bounded, value: f64, error: f64, what: &str| {
            assert_eq!((got.value, got.error), (value, error), "{what}");
        };
        same(
            prod(b(0.1), b(0.5)),
            0.05,
            0.0,
            "a power of two scales exactly",
        );
        same(
            prod(t, b(-1.0)),
            -0.1,
            1e-20,
            "a sign flip carries the radius",
        );
        same(prod(t, b(0.0)), 0.0, 0.0, "an exact zero leaves nothing");
        same(sum(t, b(0.0)), 0.1, 1e-20, "adding an exact zero is exact");
        same(diff(b(0.0), t), -0.1, 1e-20, "a negation is exact");
    }
}
