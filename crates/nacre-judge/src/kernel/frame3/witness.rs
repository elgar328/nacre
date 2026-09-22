use super::*;
/// A rational base point carried through a chain of axis rotations (a rotation history). `base` + `chain` are the exact **definition** (never lost); `realized` is
/// the f64 realization (a cache) with, per axis, the bound on its error — a **direction-wise
/// xyz vector** held beside the value it bounds, so the two cannot drift apart.
/// [`hp_coord`](Self::hp_coord) realizes the chain at
/// arbitrary precision from the definition, so two points with the same definition
/// realize identically (path-independent — the soundness argument's root).
#[derive(Clone, Debug)]
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,
    // (equality is definitional — see the `PartialEq` impl below the struct)
    pub realized: [Bounded; 3],
    /// Memoized `hp_coord` at the boolean's chosen precision — the astro-float realization
    /// once per definition and shared across clones (`Rc`). A judge escalates the *same*
    /// definition-point dozens of times per boolean (`plane_def` clones `tri_pt3` per call);
    /// without this each escalation replays the rotation's cos/sin at 200 bits, which dominates
    /// the rotated-boolean cost. `base`/`chain` never change after construction except through
    /// [`Self::rotate_about`], which resets this cell, so the cached value always matches the
    /// definition (a pure, path-independent function). See [`HpCell`] for the
    /// `Arc<OnceLock>` (parallel) vs `Rc<OnceCell>` (serial) choice.
    hp: HpCell,
}

/// **Definitional equality — `base` and `chain`, nothing else.** `realized`/`hp` are caches,
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

/// Magnitude of a `BigFloat` as an f64 power of two (0 when exactly zero) — the `f64` reading of
/// the bound [`Mag::above`] carries, for the tol arithmetic in this file, which is `f64`. Kept as
/// its own reader rather than `2^(Mag::exp2)`: that reads one octave high for a power of two.
pub(super) fn bf_mag(bf: &BigFloat) -> f64 {
    if bf.is_zero() {
        0.0
    } else {
        2f64.powi(bf.exponent().unwrap_or(0))
    }
}

/// **How far a rational's f64 image sits from the rational** — measured at high precision, `0`
/// when the value is exactly representable.
///
/// One spelling, because three call sites want the identical quantity ([`WitnessPoint::at`]'s seed,
/// [`WitnessPoint::translate`]'s offset, [`WitnessPoint::frame`]'s inputs) and "how a rational's rounding is
/// charged" is exactly the kind of thing that drifts when it is written out twice. `bf_mag` reads
/// at octave granularity, so the `2×` keeps the result above the truth.
fn rat_round_tol(r: Rat, f: f64) -> f64 {
    let e = BigFloat::from_f64(f, 120).sub(&rat_to_big(r, 120), 120, HP_RM);
    2.0 * bf_mag(&e.abs())
}

impl WitnessPoint {
    /// The realized coordinate — the values of [`Self::realized`], for readers that want the
    /// point whole.
    #[inline]
    pub fn coord(&self) -> [f64; 3] {
        self.realized.map(|b| b.value)
    }

    /// The per-axis error bounds — the radii of [`Self::realized`], for readers that want them
    /// as a vector.
    #[inline]
    pub fn tol(&self) -> [f64; 3] {
        self.realized.map(|b| b.error)
    }

    /// A point at `base`, tol seeded with the base→f64 rounding (a division; exactly
    /// 0 for an f64-representable base, positive otherwise). An axis a later chain
    /// never rotates keeps exactly this, which a per-axis tol check needs.
    pub fn at(base: [Rat; 3]) -> Self {
        let coord = [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()];
        Self::at_with_tol(
            base,
            [
                rat_round_tol(base[0], coord[0]),
                rat_round_tol(base[1], coord[1]),
                rat_round_tol(base[2], coord[2]),
            ],
        )
    }

    /// A point at a rational `base`, its coordinate the nearest `f64` and its tol **from that
    /// contract rather than measured**: exactly `0` where the coordinate is representable, else
    /// `|x|·2⁻⁵³` (floored at `f64::MIN_POSITIVE`) — `Rat::to_f64` is documented as *"the
    /// **nearest** f64 … ties to even"*, so `|r − to_f64(r)| ≤ ½ ulp` holds for every normal
    /// `x` without anyone computing it.
    ///
    /// [`at`](Self::at) learns the same number by measuring at 120 bits — nine BigFloat
    /// operations per point. This is the spelling `nacre-ops` uses where the base is a
    /// **definition** (a plane's own points, a solid's vertices solved from their carriers): the
    /// bound is free, and its looseness can only turn a definite filter answer into "escalate",
    /// never the other way round.
    ///
    /// ★ **This is how a point that came from a rational reaches the predicates** — not by
    /// lifting the rounded `f64` back into a `Rat` (that names a different point, with a
    /// tolerance of zero on top).
    pub fn at_nearest(base: [Rat; 3]) -> Self {
        let f = base.map(|r| r.to_f64());
        let representable = base
            .iter()
            .zip(f)
            .all(|(&r, x)| Rat::try_from_f64(x) == Some(r));
        if representable {
            return Self::at_with_tol(base, [0.0; 3]);
        }
        let bound = |x: f64| (x.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE);
        Self::at_with_tol(base, [bound(f[0]), bound(f[1]), bound(f[2])])
    }

    /// A point at `base` with an explicit initial `tol` — for a root that already
    /// carries tol (a `Discovered` boolean seam), whose tol the chain then transports
    /// (`|R|·old`). [`at`](Self::at) is the Constructed case (initial tol = base
    /// rounding).
    pub fn at_with_tol(base: [Rat; 3], tol: [f64; 3]) -> Self {
        WitnessPoint {
            realized: [0, 1, 2].map(|k| Bounded::new(base[k].to_f64(), tol[k])),
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
    /// 90°-family angles about the origin rotate exactly (tol 0, Niven). A non-origin
    /// pivot adds its own f64 rounding (`coord − p`, `p + …`, the pivot's Rat→f64) —
    /// present even for an exact angle — covered by the `piv` term (checked by `tol_bounds_error_over_random_chains`).
    ///
    /// **Two error terms, and neither is a guess.** `rot` carries the one input here without a
    /// rounding contract (`f64::cos`) and charges what `Angle::realization_error_of` *measured* of
    /// the very pair used below; `piv` is nothing but round-to-nearest steps and a realization, so it
    /// is derived and measured the way [`translate`](Self::translate) and [`mirror`](Self::mirror)
    /// derive and measure theirs.
    ///
    /// ★★ **`rot` is per axis, and it has to be.** The two coordinates are not the same expression:
    /// `i` mixes `u·c − v·s` while `j` mixes `u·s + v·c`, so `u` pairs with `cos`'s error on one
    /// and with `sin`'s on the other. A single shared term is only sound when the two errors are
    /// charged the same amount — which a *constant* did, and a measurement does not. Measured on a
    /// 356.65° turn, `sin` was eight times further out than `cos`, and the axis whose `u` met `sin`
    /// needed 2.4× what the other did. The old shared charge hid that; the ground-truth test found
    /// it the moment the constant left.
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
        let (px, py) = (point[i].to_f64(), point[j].to_f64());
        let (ci, cj) = (self.realized[i].value, self.realized[j].value); // pre-rotation magnitudes for tol
        let (u, v) = (ci - px, cj - py);
        // cos/sin — exact (rational) for the 90°-family, else f64 (with realization tol).
        //
        // **Through `Angle`'s single entry point, not re-spelled here.** This function's contract
        // is to redo the producer's f64 route operation for operation (`Isometry::apply_point`,
        // which calls the same thing), and a second spelling of "how an angle becomes f64" is
        // exactly how the two would drift — silently, at the 90°-family, where one route snaps to
        // `0.0`/`±1.0` and the other lands `cos(90°) ≈ 6e-17`.
        let (c, s) = angle.cos_sin_f64();
        self.realized[i].value = px + u * c - v * s;
        self.realized[j].value = py + u * s + v * c;
        // **The rotation's own error — measured for this angle, not charged from a constant.**
        //
        // `dc`/`ds` are how far this platform's `cos`/`sin` land from the truth — measured, not a
        // constant measured once, which would be sound only on machines like the one it was taken
        // on, and a kernel that ships to browsers cannot assume that. A worse platform reports a
        // bigger number and the tolerance grows to match.
        //
        // ★ **Per axis, because the two coordinates mix the pair differently** — `u` meets `c` on
        // one and `s` on the other. See this function's doc for what that cost when it was shared.
        //
        // **The arithmetic of the products and their combination lives here too, and only here.**
        // For an origin pivot `piv` below is skipped entirely, so nothing else would cover
        // `fl(u·c)`, `fl(v·s)` and their combination. Each is `≤ ε/2` of a magnitude bounded by
        // `|u| + |v|`, so `1·ε` of that sum covers all three.
        //
        // ★ **Both are gated on the realization being inexact, and that gate is the invariant.**
        // When `cos`/`sin` are `0`/`±1` the products and their combination are *exact*, not merely
        // small — so a quadrantal origin rotation contributes nothing at all, which is what keeps
        // it at tol 0 (`quadrantal_origin_chain_is_tol_zero`) and keeps axis-aligned models on the
        // exact predicate path. An unconditional arithmetic term would break that quietly.
        //
        // ★★★ **It measures `c` and `s` — the values three lines above — not the angle.** They are
        // handed in rather than re-derived because a second realization of one angle is measured to
        // be able to differ from the first (`Angle`'s `F64_ERR` says where and why). An error taken
        // against a pair that never reached `coord` would bound a number this kernel never stored.
        let (dc, ds) = angle.realization_error_of(c, s);
        let (rot_i, rot_j) = if dc == 0.0 && ds == 0.0 {
            (0.0, 0.0)
        } else {
            let arith = f64::EPSILON * (u.abs() + v.abs());
            (
                u.abs() * dc + v.abs() * ds + arith,
                u.abs() * ds + v.abs() * dc + arith,
            )
        };
        // **The pivot arithmetic, charged the way `translate` and `mirror` charge theirs.**
        //
        // Not a trig constant — which would exist only *because `f64::cos` has no accuracy
        // contract* — for three operations that all do have one. Like its two siblings in this
        // `impl`, it **measures** the rational's realization and **counts** the round-to-nearest
        // steps.
        //
        // - The pivot's `Rat → f64` is measured, not counted, for the same reason they measure it:
        //   `Rat::to_f64` is one rounding for a numerator and denominator under `2⁵³` and takes
        //   another path above, so counting would mean knowing which. It enters **twice per axis
        //   with opposite signs** (`px + (ci − px)·c − …`) and partly cancels; `|1 − c| ≤ 2` and
        //   `|s| ≤ 1` bound the pair plainly, and the outer `2.0` is the margin every sibling
        //   term carries. **Exactly zero for a dyadic pivot** — the common case, which a lumped
        //   charge would still bill.
        // - **The roundings the pivot itself adds are three**: the two differences `ci − px` and
        //   `cj − py`, and the final sum `px + …`. (The products and their combination belong to
        //   `rot` above — they happen whether or not there is a pivot.) Each is `≤ ε/2` of a
        //   magnitude bounded by `|ci| + |cj| + |px| + |py|`, so `1.5·ε` of that sum covers them;
        //   `3.0` is that with the same doubling.
        //
        // Origin pivot stays exactly 0: `ci − 0.0` and `0.0 + x` are exact, so there is nothing
        // to charge — and the whole term is skipped rather than measured. This term *is* symmetric
        // in the two axes, unlike `rot`: it is built from magnitudes, not from which of `cos`/`sin`
        // each coordinate met.
        let piv = if px != 0.0 || py != 0.0 {
            let realized = |r: Rat, f: f64| {
                bf_mag(&rat_to_big(r, 120).sub(&BigFloat::from_f64(f, 120), 120, HP_RM)).abs()
            };
            2.0 * (2.0 * realized(point[i], px) + realized(point[j], py))
                + 3.0 * f64::EPSILON * (ci.abs() + cj.abs() + px.abs() + py.abs())
        } else {
            0.0
        };
        let (ti, tj) = (self.realized[i].error, self.realized[j].error);
        self.realized[i].error = c.abs() * ti + s.abs() * tj + rot_i + piv;
        self.realized[j].error = s.abs() * ti + c.abs() * tj + rot_j + piv;
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
    /// **`coord` follows the producer's own route** (`AxisMirror::point`), so a replay of the
    /// definition reproduces the stored coordinate bit for bit. The exact offset lives in the
    /// chain, where [`compute_hp`](Self::compute_hp) realizes it.
    pub fn mirror(self, axis: Axis, offset: Rat) -> Self {
        let mut p = self.mirrored(axis, offset);
        p.remember(MoveNode::Mirror { axis, offset });
        p
    }

    /// [`Self::mirror`]'s numbers, without remembering the node.
    fn mirrored(mut self, axis: Axis, offset: Rat) -> Self {
        let k = axis.index();
        let c = offset.to_f64();
        // `2·c` is exact (a power-of-two multiply), so the new error is the offset's own
        // realization plus the subtraction's half-ulp. Negation itself is exact, so the incoming
        // tol passes through unscaled.
        //
        // **The offset enters doubled** (`2·c`), so its realization error does too — the factor
        // here is `4 = 2 (the doubling) × 2 (the same safety margin every other term carries)`.
        // Copying `translate`'s `2.0` looks right and is not: there the offset enters once, so its
        // `2.0` *was* the margin. Measured — a chain of two reflections overran the bound by 4%.
        self.realized[k].error += 4.0
            * bf_mag(&rat_to_big(offset, 120).sub(&BigFloat::from_f64(c, 120), 120, HP_RM)).abs()
            + f64::EPSILON * (2.0 * c.abs() + self.realized[k].value.abs());
        // The producer's own route (`AxisMirror::point`), operation for operation — a replay of
        // the definition has to reproduce the stored coordinate bit for bit.
        self.realized[k].value = 2.0 * c - self.realized[k].value;
        self
    }

    /// The point translated by an exact rational `offset` — one more link in the definition.
    ///
    /// **`coord` is updated exactly as the producer does it** (`Isometry::apply_point`: add the
    /// offset's f64 image), so a replay of the definition reproduces the stored coordinate bit for
    /// bit. The exact offset lives in the chain, where [`compute_hp`](Self::compute_hp) realizes
    /// it; the two roundings this f64 step takes (the offset's own, and the add's) are what `tol`
    /// grows by.
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
        for (k, &off) in offset.iter().enumerate() {
            let t = off.to_f64();
            // The offset's realization error, plus the add's own half-ulp on the result.
            self.realized[k].error += 2.0
                * bf_mag(&rat_to_big(off, 120).sub(&BigFloat::from_f64(t, 120), 120, HP_RM)).abs()
                + f64::EPSILON * (self.realized[k].value.abs() + t.abs());
            self.realized[k].value += t;
        }
        self
    }

    /// The point read as `(u, v, w)` in the frame of the plane `origin`/`u_raw`/`n` describe, and
    /// written out in that plane's own frame — one more link in the definition.
    ///
    /// **`coord` follows the producer's route operation for operation**, so a replay of the
    /// definition reproduces the stored coordinate bit for bit. The exact inputs live in the
    /// chain, where [`compute_hp`](Self::compute_hp) realizes them.
    ///
    /// `None` when a squared length is not positive — that is a degenerate frame, not an unusual
    /// one, and a caller that cannot build a basis must not get a point that pretends otherwise.
    ///
    /// ★★★★ **Everything charged here is measured, and the exact cases really do reach zero.**
    /// The three sources are the rational inputs' own `Rat → f64` rounding, the two `1/√` values'
    /// realization error ([`nacre_exact::inv_sqrt_error_of`] — *not* a half-ulp taken on faith,
    /// see there), and the f64 arithmetic that combines them. An axis-aligned frame has a signed
    /// permutation for a basis and every one of those terms vanishes: `±1` and `0` are exact, so
    /// the products are exact and the sums pick out one coordinate each.
    pub fn frame(self, f: nacre_exact::PlaneFrame) -> Option<Self> {
        let mut p = self.framed(f)?;
        p.remember(MoveNode::Frame { frame: f });
        Some(p)
    }

    /// [`Self::frame`]'s numbers, without remembering the node.
    fn framed(mut self, f: nacre_exact::PlaneFrame) -> Option<Self> {
        // One inverse square root per axis, each of an exact rational — see `plane_frame` for why
        // `v̂` gets its own instead of being a cross product of the other two.
        let inv = |v: Rat| -> Option<(f64, f64)> {
            let x = nacre_exact::inv_sqrt_f64(v)?;
            Some((x, nacre_exact::inv_sqrt_error_of(v, x)?))
        };
        let ((iu, du), (iw, dw)) = (inv(f.uu)?, inv(f.nn)?);
        // A basis vector's component, and a bound on how far it lands from the true one: the
        // rational's own rounding scaled by the length, the length's error scaled by the
        // rational, and the product's half-ulp. All three vanish for a `0`/`±1` component.
        let axis_comp = |r: Rat, s: f64, ds: f64| {
            let rf = r.to_f64();
            let c = rf * s;
            let err = rat_round_tol(r, rf) * s.abs() + rf.abs() * ds;
            (
                c,
                if err == 0.0 {
                    0.0
                } else {
                    err + f64::EPSILON * c.abs()
                },
            )
        };
        let (mut uh, mut vh, mut wh) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut eu, mut ev, mut ew) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        for k in 0..3 {
            (uh[k], eu[k]) = axis_comp(f.u_raw[k], iu, du);
            (wh[k], ew[k]) = axis_comp(f.n[k], iw, dw);
        }
        // ★★ `v̂` exactly where the frame carries `v_raw`, and `ŵ × û` where it does not.
        //
        // The exact route is one rounding per component; the cross product is two that do not
        // cancel — a wall whose `v` is exactly `ẑ` comes out three ulps short through it. Both are
        // sound; the fallback is what this crate had before `v_raw` existed, and it is taken only
        // when `|n|²·|u_raw|²` does not fit `i128` (see `nacre_exact::plane_frame`).
        match f.v {
            Some((v_raw, vv)) => {
                let (iv, dv) = inv(vv)?;
                for k in 0..3 {
                    (vh[k], ev[k]) = axis_comp(v_raw[k], iv, dv);
                }
            }
            None => {
                // `|ŵ|, |û| ≤ 1`, so each of the four products carries the other factor's error at
                // unit scale; `2·(max eu + max ew)` covers all four with room. Gated on the basis
                // being inexact for the reason `rotate_about` gates its own: with `0`/`±1` factors
                // the products *and* their differences are exact, not merely small.
                let mx = |e: [f64; 3]| e.iter().fold(0.0f64, |a, &b| a.max(b));
                let (eu_max, ew_max) = (mx(eu), mx(ew));
                let e = if eu_max == 0.0 && ew_max == 0.0 {
                    0.0
                } else {
                    2.0 * (eu_max + ew_max) + 3.0 * f64::EPSILON
                };
                for k in 0..3 {
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    vh[k] = wh[i] * uh[j] - wh[j] * uh[i];
                    ev[k] = e;
                }
            }
        }
        // The combination itself: three products and three sums per coordinate, each a
        // round-to-nearest of a magnitude bounded by the terms' sum — except where the basis is a
        // signed permutation *and* the frame sits on the origin, where every product is exact and
        // every sum picks out a single coordinate. That case is the axis-aligned sketch, and it
        // is the one that has to stay at tol 0.
        let exact_basis = eu.iter().chain(&ev).chain(&ew).all(|&e| e == 0.0);
        let perm = exact_basis
            && uh
                .iter()
                .chain(&vh)
                .chain(&wh)
                .all(|c| *c == 0.0 || c.abs() == 1.0);
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let ok = f.origin[k].to_f64();
            let terms = p[0] * uh[k] + p[1] * vh[k] + p[2] * wh[k];
            // The incoming tol turned by the same basis, plus what this step adds.
            let carried = uh[k].abs() * t[0] + vh[k].abs() * t[1] + wh[k].abs() * t[2];
            let realized = p[0].abs() * eu[k] + p[1].abs() * ev[k] + p[2].abs() * ew[k];
            let arith = if perm && ok == 0.0 {
                0.0
            } else {
                rat_round_tol(f.origin[k], ok)
                    + 3.0
                        * f64::EPSILON
                        * (ok.abs()
                            + (p[0] * uh[k]).abs()
                            + (p[1] * vh[k]).abs()
                            + (p[2] * wh[k]).abs())
            };
            self.realized[k].value = ok + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`WideFrame`] — the same propagation, with the axis and origin
    /// `(value, error)` pairs taken from a fixed-precision arbitrary-precision realization
    /// instead of `inv_sqrt_f64`/`axis_comp`. No new f64 error derivation exists here: every
    /// rounding on the way is inside an `HpBounded`, and the final narrowing to f64 charges itself.
    ///
    /// ★ The axes are realized per applied point. The wide population is a fraction of a percent
    /// of pushes, so that cost is accepted rather than memoized — a measure-later item.
    ///
    /// ★ No exact-permutation shortcut: a wide frame is never an axis permutation (its squared
    /// lengths exceed `i128`), so the tol-0 branch `WitnessPoint::frame` has cannot apply.
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
        let inv = |v: &num_bigint::BigInt| nacre_exact::inv_sqrt_bigint_bounded(v, P);
        let (iu, iv2, iw) = (inv(&f.uu)?, inv(&f.vv)?, inv(&f.nn)?);
        let comp = |raw: &num_bigint::BigInt, s: &HpBounded| {
            narrow_hp(&HpBounded::of_bigint(raw, P).mul(s, P))
        };
        let (mut uh, mut vh, mut wh) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut eu, mut ev, mut ew) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let (mut oh, mut eo) = ([0.0; 3], [0.0; 3]);
        let od = HpBounded::of_bigint(&f.origin_den, P);
        for k in 0..3 {
            (uh[k], eu[k]) = comp(&f.u_raw[k], &iu);
            (vh[k], ev[k]) = comp(&f.v_raw[k], &iv2);
            (wh[k], ew[k]) = comp(&f.n[k], &iw);
            (oh[k], eo[k]) =
                narrow_hp(&HpBounded::of_bigint(&f.origin_num[k], P).div_exact(&od, P)?);
        }
        // The same propagation as `WitnessPoint::frame`, with the realized origin's own error in place
        // of `rat_round_tol`.
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let terms = p[0] * uh[k] + p[1] * vh[k] + p[2] * wh[k];
            let carried = uh[k].abs() * t[0] + vh[k].abs() * t[1] + wh[k].abs() * t[2];
            let realized = p[0].abs() * eu[k] + p[1].abs() * ev[k] + p[2].abs() * ew[k];
            let arith = eo[k]
                + 3.0
                    * f64::EPSILON
                    * (oh[k].abs()
                        + (p[0] * uh[k]).abs()
                        + (p[1] * vh[k]).abs()
                        + (p[2] * wh[k]).abs());
            self.realized[k].value = oh[k] + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// [`WitnessPoint::frame`] for a [`FrameThrough`] — the judged basis realized
    /// at the fixed rung, narrowed to `(value, error)` pairs, then **the same propagation as
    /// [`WitnessPoint::frame_wide`]**: the incoming tol is turned by the basis, the basis's own
    /// realization error is scaled by the coordinates, and the combination charges its rounding.
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
        let p = self.coord();
        let t = self.tol();
        for k in 0..3 {
            let (uh, eu) = u[k];
            let (vh, ev) = v[k];
            let (wh, ew) = w[k];
            let (oh, eo) = o[k];
            let terms = p[0] * uh + p[1] * vh + p[2] * wh;
            let carried = uh.abs() * t[0] + vh.abs() * t[1] + wh.abs() * t[2];
            let realized = p[0].abs() * eu + p[1].abs() * ev + p[2].abs() * ew;
            let arith = eo
                + 3.0
                    * f64::EPSILON
                    * (oh.abs() + (p[0] * uh).abs() + (p[1] * vh).abs() + (p[2] * wh).abs());
            self.realized[k].value = oh + terms;
            self.realized[k].error = carried + realized + arith;
        }
        Some(self)
    }

    /// The coordinate realized at `prec` bits from the **definition** (base rotated
    /// through the chain, each node about its pivot) — path-independent ground truth /
    /// escalation realization. The result is memoized in
    /// [`WitnessPoint::hp`] and shared across clones of the same definition, so a definition-point pays
    /// the astro-float cos/sin once per boolean rather than once per predicate.
    pub(crate) fn hp_coord(&self, prec: usize) -> [HpBounded; 3] {
        // **Keyed by precision, and that key is load-bearing.** The precision is chosen per
        // boolean, so within one operation every call arrives with the same value and the cell is
        // filled once — which is the whole point, since re-realizing a definition per predicate
        // was 77% of a rotated boolean's runtime. A call at a different precision (the rare
        // per-judgement fallback) recomputes without disturbing the cached value, rather than
        // silently returning coordinates realized at the wrong precision.
        let cached = self.hp.get_or_init(|| (prec, self.compute_hp(prec)));
        if cached.0 == prec {
            return cached.1.clone();
        }
        self.compute_hp(prec)
    }

    /// **The coordinate realized at `prec` bits from the definition, with the error it carries** —
    /// the door [`Self::hp_coord`] is behind, in the type this workspace hands out publicly
    /// (`nacre_exact::HpBounded`, what `inv_sqrt_bounded` returns too).
    ///
    /// ★ **The radius is half the answer.** A value without the bound its realization cost cannot
    /// be rounded honestly — see [`HpBounded`]'s contract: a radius invented for convenience makes
    /// every sign above it unearned. Callers round with `nacre_exact::round_to_f64` /
    /// `round_to_digits`, which report *undecided* rather than picking when `prec` is short.
    ///
    /// ★★ Calling this at a precision the boolean is not using is safe: [`Self::hp_coord`]'s memo
    /// is keyed by precision and recomputes rather than disturbing the cached value.
    pub fn realize(&self, prec: usize) -> [HpBounded; 3] {
        self.hp_coord(prec)
    }

    /// The uncached realization (the body of [`Self::hp_coord`]), **with the error it carries**.
    ///
    /// This is the same walk as [`rotate_about`](Self::rotate_about)'s tol propagation, one level
    /// up: the definition is exact, the realization is not, and the radius is what the realization
    /// cost. Nothing here is a chosen constant — the base contributes its own rounding (zero when
    /// the rational lands on a `prec`-bit dyadic), each `cos`/`sin` contributes the bound
    /// [`Angle::cos_sin_bounded`] derives, and every arithmetic operation adds its half-ulp.
    pub(super) fn compute_hp(&self, prec: usize) -> [HpBounded; 3] {
        let p = [
            HpBounded::of_rat(self.base[0], prec),
            HpBounded::of_rat(self.base[1], prec),
            HpBounded::of_rat(self.base[2], prec),
        ];
        fold_suffix(p, &self.chain, prec)
    }
}
