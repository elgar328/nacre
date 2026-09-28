//! **Asking a vertex for its coordinate at a precision you choose.**
//!
//! Everything the kernel shows — the viewport, the tessellation, a report row, the STEP file —
//! reads the f64 point cache. Since the cache became this module's memo (`push_vertex_realized`:
//! a vertex is realized from its definition the moment an operation pushes it) that cache *is* the
//! nearest f64 wherever the first rung decides it; printing more digits of it would still print
//! the rounding, not the point, which is what `realize_vertex_decimal` is for.
//!
//! This module goes the other way: it takes the vertex's **definition** and realizes a coordinate
//! from it, rounding exactly once at the end. Two roads meet here and they are chosen by the
//! *value*, not by the `Vertex` variant:
//!
//! - **rational** — a three-plane meet is an exact ratio ([`nacre_topo::Model::vertex_meet`]), and
//!   so is that meet carried through a chain that folds (`Model::chain_point_rat`), so its decimals
//!   come out of one long division and every digit printed is a digit of the coordinate itself.
//!   There is no rounding question to get wrong.
//! - **realized** — anything a turn off the quarters, a frame or a radical reaches is approached at
//!   `prec` bits with the error the realization cost, and a digit is printed only once the interval
//!   decides it.
//!
//! ★ **Where this sits.** `nacre-topo` does not depend on `nacre-judge`, so the composition
//! `vertex_meet → WitnessPoint::at → replay → realize` cannot live on `Model`; `nacre-ops` is the
//! first crate that sees both. The rounding itself lives one layer further down, in
//! `nacre-exact`, beside `round_to_f64` — cip realizes, scalar rounds, this module composes.
//!
//! ★★ **Migration 3b.** [`nacre_topo::Vertex::OnSeam`]'s doc records the missing piece as *"the
//! machinery that regenerates the cached coordinate"*, deferred with row 3b. This is the asking
//! half of it. It does **not** overwrite the cache — that is the row's other half.

use crate::rotated_vertex::{motion_chain, replay};
use nacre_exact::{HpBounded, Mag, MeetPoint};
use nacre_judge::WitnessPoint;
use nacre_math::Point3;
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, PrefixKey, Surface, Vertex};
use num_bigint::BigInt;

/// How precisely to realize — always stated, never defaulted.
///
/// ★ There is no `Digits` here on purpose: the digit count belongs to
/// [`Realized::to_decimal`] alone. Carrying it in both places lets the two disagree, and the way
/// they disagree is by printing digits the realization never determined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    /// Realize from the definition and read it out once ([`Realized::to_f64`] — the nearest
    /// `f64`, or `+0.0` within the coincidence limit), escalating until the interval names one.
    /// This is the value the cache *should* hold.
    NearestF64,
    /// An explicit working precision — for measurement, and for a caller driving its own ladder.
    Bits(usize),
}

/// A coordinate that has been realized, together with what it cost.
///
/// Opaque: the arms are an implementation detail, and a caller that could see them would be
/// coupled to astro-float's types through this crate as well as through `nacre-exact`.
#[derive(Clone, Debug)]
pub struct Realized(Arm);

#[derive(Clone, Debug)]
enum Arm {
    /// Exact: three numerators over one positive denominator ([`MeetPoint::lift`]'s shape, which
    /// is why a coordinate too wide for `Rat` is not a refusal here).
    Exact([BigInt; 3], BigInt),
    /// Approached: value and error radius per coordinate, **with the precision that produced
    /// them** — the rounding predicates need it, and re-deriving it at the rounding site is how a
    /// climbed ladder silently rounds at the bottom rung.
    Approached([HpBounded; 3], usize),
}

/// Why a vertex could not be realized. **Never falls back to the cache** — a measurement that
/// silently answers with the thing it is measuring reports nothing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RealizeError {
    /// [`nacre_topo::Model::vertex_meet`] declined: carriers carrying two motion histories, a
    /// carrier whose name is `Wide` (its licence reads narrow coefficients), a carrier with no
    /// recorded name, or three carriers meeting in no point.
    NoMeet,
    /// The meet is wider than `Rat` **and** the vertex carries a motion. The exact road is open
    /// for either alone; replaying a motion needs a rational base
    /// (`WitnessPoint::at` takes `[Rat; 3]`) and no wide constructor exists yet.
    WideUnderMotion,
    /// The vertex's motion chain could not be rebuilt exactly.
    NoMotionChain,
    /// A curved definition (`OnSeam`, `Pierce`) did not resolve into a point.
    ///
    /// ⚠ **This is a bag, and saying so is the point.** It covers: a carrier that cannot be
    /// stated in the world exactly, a cap plane that is not perpendicular to the axis (so it
    /// bounds no rim), a `root` that does not match the kind of crossing the carriers actually
    /// make (a `Double` asked of a `Pair`), and rational overflow on the way. Each deserves its
    /// own name; none of them is [`Self::NoMeet`], which is a *different* function declining —
    /// `Model::vertex_meet` is never called on this road.
    NoCurvedPoint,
    /// The realization ladder reached its ceiling without deciding what was asked.
    Undecided,
    /// The value is an exact rational, and **no `f64` names it**: it is outside the range one can
    /// hold. More bits cannot help — this is the one refusal a taller ladder does not answer.
    ///
    /// ⚠ Measured at both ends (`nacre_exact`'s own lock): the large end refuses at `2^1024`,
    /// where the scaled value stops being finite, and the small end at `2^-1149`. A carrier name
    /// cannot reach either — the widest this repository has measured is 168 bits — so the
    /// population is zero today. It has a name anyway, because the branch is real and
    /// [`Self::Undecided`] would be a lie about it ("more bits would do it" — they would not).
    Unrepresentable,
}

/// Why the **cache road** did not answer — two reasons, and only one of them is a failure.
///
/// ★ Deliberately not a `RealizeError` variant. That enum's own doc says it answers *"why a vertex
/// could not be realized"*, and a chain refused for cost was never attempted: "did not" is not
/// "cannot". Folding the two together would also hand every caller of [`realize_vertex`] a variant
/// they can never receive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CacheDecline {
    /// The motion history is deeper than [`CACHE_REPLAY_COST_CAP`] — the road did not walk it.
    /// A paid realization still can, which is what `refine_vertex_cache` is for.
    CostCap,
    /// The realization itself declined, by name.
    Cannot(RealizeError),
}

/// The ladder. Doubling, so a value needing `n` bits pays at most `2n`; capped, because an
/// unbounded climb on an undecidable ask is a hang rather than an answer.
///
/// ★ The last rung **is** the judgement's own cap ([`crate::planes::JUDGE_PREC_CAP`]) rather than
/// a second copy of the same digits: both name the precision this kernel is willing to pay for, so
/// they move together instead of agreeing by coincidence. ⚠ The price of tying them, stated: a
/// model past the judging cap can still be *exported*, so a realization deeper than this is not
/// meaningless in principle — it is simply not offered, and the caller is told `Undecided` rather
/// than made to wait.
const LADDER: [usize; 6] = [128, 256, 512, 1024, 2048, crate::planes::JUDGE_PREC_CAP];

impl Realized {
    /// The nearest `f64` per coordinate, with the error each carries — **or `+0.0` where the
    /// coordinate is proven within the coincidence limit of zero**. `None` where the realization
    /// names neither.
    ///
    /// ★★★ **The coincidence rule, as judging states it.** Two things proven closer than
    /// `2⁻¹⁸⁰` of the model's size are one (design 「숫자 규칙」 5: below that no split survives
    /// any output). A coordinate that is exactly `0` never decides by rounding — its interval
    /// straddles `0` at every rung, and `-tiny` and `+tiny` are different `f64`s — so an interval
    /// entirely inside `± max(1, largest |coordinate|) · 2⁻¹⁸⁰` is read as `+0.0`, with that
    /// interval's upper bound as its error, which keeps «the truth lies within `coord ± bound`»
    /// true. The scale is the point's own, which is never larger than the model's: the limit errs
    /// strict, and erring strict only costs bits. It is asked **before** rounding, so every rung
    /// that decides gives the same answer — asked after, a true `1e-70` would read `0` where
    /// rounding fails and `1e-70` a rung later. An exact arm is not asked: a rational rounds.
    ///
    /// ★★★ **A decided coordinate's error is half an ulp, on both arms** (`half_ulp`). The
    /// realization's own radius bounds its *interval*, not the `f64` read out of it: the truth
    /// sits within `2⁻¹²⁸` of the interval's middle and still up to half an ulp from the value
    /// the readout reports. What makes half an ulp *proven* rather than typical is the rounding
    /// predicate itself — [`nacre_exact::round_to_f64`] answers only when both ends of the
    /// interval round to the same `f64`, so the truth, which lies between them, rounds there too.
    /// The ladder radius handed over as the bound instead leaves `coord ± bound` short of the
    /// truth by up to half an ulp (`a_bounded_cache_contains_the_coordinate_its_definition_names`).
    ///
    /// ★★ **The error comes back as [`Mag`], not `f64`.** `nacre-exact`'s own test says why: *"a
    /// radius the ladder actually produces must not become zero. An `f64` cannot hold
    /// `2⁻²⁰⁴⁸`."* The coincidence arm carries exactly such a radius (an interval around `0` at
    /// 4096 bits), and an earlier spelling that converted it reported `0e0` — indistinguishable
    /// from the exact arm's honest zero.
    ///
    /// The *value* is an `f64` because that is what the caller asked for. The *bound* on it need
    /// not be one, and on the coincidence arm cannot be.
    pub fn to_f64(&self) -> Option<([f64; 3], [Mag; 3])> {
        match &self.0 {
            Arm::Exact(n, d) => {
                let mut v = [0.0; 3];
                let mut e = [Mag::ZERO; 3];
                for k in 0..3 {
                    // ⚠★★★ **An exact realization is not an exact `f64`.** The rational is the
                    // truth; reading it out at 53 bits rounds, and a 59-bit coordinate does not
                    // fit — measured on the tilted-frame family, where the value is
                    // `0.130864196953086372` and its `f64` is `0.13086419695308637578…`. An
                    // earlier spelling reported `[0.0; 3]` here, which is the cache's own lie in
                    // a new place, and the audit lock written for it *enforced* the lie.
                    let (val, no_loss) = nacre_exact::nearest_f64_big_exact(&n[k], d)?;
                    v[k] = val;
                    e[k] = match no_loss {
                        true => Mag::ZERO,
                        false => half_ulp(val),
                    };
                }
                Some((v, e))
            }
            Arm::Approached(p, prec) => {
                let mut v = [0.0; 3];
                let mut e = [Mag::ZERO; 3];
                let scale = p.iter().fold(Mag::of(1.0), |s, c| {
                    let m = Mag::above(&c.value);
                    if s.lt(m) { m } else { s }
                });
                let coincidence = scale.times(Mag::pow2(-180));
                for k in 0..3 {
                    let reach = Mag::above(&p[k].value).plus(p[k].error);
                    if reach.lt(coincidence) {
                        e[k] = reach;
                        continue;
                    }
                    v[k] = nacre_exact::round_to_f64(&p[k].value, p[k].error, *prec)?;
                    e[k] = half_ulp(v[k]);
                }
                Some((v, e))
            }
        }
    }

    /// `places` **decimal places** per coordinate, or `None` when this realization does not
    /// determine them — the signal to realize again at more bits.
    pub fn to_decimal(&self, places: usize) -> Option<[String; 3]> {
        match &self.0 {
            Arm::Exact(n, d) => Some(core::array::from_fn(|k| {
                nacre_exact::decimals_of_ratio(&n[k], d, places)
            })),
            Arm::Approached(p, _) => {
                let mut out = [const { String::new() }; 3];
                for k in 0..3 {
                    out[k] = nacre_exact::round_to_digits(&p[k].value, p[k].error, places)?;
                }
                Some(out)
            }
        }
    }

    /// True when the coordinate is an exact rational — every digit [`Self::to_decimal`] prints is
    /// then a digit of the point, not of an approximation to it.
    pub fn is_exact(&self) -> bool {
        matches!(self.0, Arm::Exact(..))
    }
}

/// **How far a correctly rounded `f64` can sit from the value it rounds** — `|v| · 2⁻⁵³`, which is
/// at least half an ulp of `v` on either side of a binade edge. As a [`Mag`], so the bound of a tiny
/// coordinate cannot vanish on the way out. One rule for both arms of [`Realized::to_f64`]: an exact
/// rational read at 53 bits and an interval that decided its `f64` are both a correct rounding.
fn half_ulp(v: f64) -> Mag {
    Mag::of(v).times(Mag::pow2(-53))
}

/// Realize one vertex's coordinate from its definition.
pub fn realize_vertex(
    model: &Model,
    v: Handle<Vertex>,
    p: Precision,
) -> Result<Realized, RealizeError> {
    realize_def(model, model.vertex(v), p)
}

/// [`realize_vertex`] on a definition that has not been pushed yet — the road every vertex an
/// operation makes takes *before* it exists, so its cache is the realization from the start.
pub(crate) fn realize_def(
    model: &Model,
    def: &Vertex,
    p: Precision,
) -> Result<Realized, RealizeError> {
    realize_def_tracked(model, def, p, &mut None)
}

/// [`realize_def`] that also reports what the accelerator could keep.
///
/// ★ Only the fixed-precision road fills `out`. A climb asks the same definition at rung after
/// rung, so whatever it learned at 128 bits is not what it returns, and filing the rung it happened
/// to pass through would put a value under a key the next reader asks at a different precision.
pub(crate) fn realize_def_tracked(
    model: &Model,
    def: &Vertex,
    p: Precision,
    out: &mut Option<PrefixWrite>,
) -> Result<Realized, RealizeError> {
    match p {
        Precision::Bits(bits) => build(model, def, bits, out),
        Precision::NearestF64 => climb(model, def, |r| r.to_f64().map(|_| r)),
    }
}

/// **How deep a motion history the cache road replays — a `cost` limit, not a precision rule.**
///
/// ⚠★★★★ **That distinction is the whole of this constant, and collapsing it cost a cell.** The
/// guard this replaces read as *"past here the arithmetic cannot decide anyway"*, and measured,
/// that is **false**: 150 rational translations and 80 turns of 7° decide perfectly at 128 bits,
/// while a single 37° turn stops deciding past the 63rd. What a node costs in bits depends on its
/// **angle** (`Angle::try_exact_cos_sin` answers only 0/90/180/270 — Niven, so those add no
/// radius), not on how many nodes came before it. Precision is decided by the first rung's own
/// verdict; this number decides only **how much work a push is willing to do**.
///
/// **Where the value comes from** (measured: one vertex realized at 128 bits, 200
/// repetitions). The replay is linear at ~1.35 µs per node:
///
/// | depth | 1 | 128 | 192 | 208 | 224 | 256 |
/// |---|---|---|---|---|---|---|
/// | ms per vertex | 0.014 | 0.180 | **0.264** | 0.291 | 0.314 | 0.355 |
///
/// The budget is **one vertex stays under 0.3 ms**, which lands between 208 and 224. 192 sits
/// inside it with room for the instrument's own 4 % spread, and is 18× what an ordinary vertex
/// (depth ≤ 2) pays. Change the budget and the number follows — that is why the formula is
/// written here and not just the digits.
///
/// ☑ **And the value carries no correctness load.** Real models measure depth ≤ 2 (census rows by
/// their deepest chain: 345 at 0, 53 at 1, 3 at 2), so every value above 2 behaves identically on
/// them; what this decides is how much is done *eagerly*, at push time. Past it, on a chain that
/// has to be replayed, the construction's own figure stands, and a caller who wants the point exactly still has
/// [`realize_vertex`], which climbs.
///
/// ⚠ **Not derived from [`crate::planes::JUDGE_PREC_CAP`].** That one caps the precision a
/// *judgement* will pay for and bites at roughly four thousand turns; this one caps what a *cache*
/// will pay for and bites two orders of magnitude earlier. Same kind of limit, different scale.
///
/// ☑ **The value survives the prefix accelerator — measured, not assumed.**
/// The expectation would be that this constant loses its ground: if a realization can
/// resume from a remembered prefix, one step should cost the same at any depth, and a budget in
/// milliseconds would hold everywhere. Measured (release, timing **one** transform at depth rather
/// than averaging a build, which halves every figure):
///
/// | depth | 190 | 1000 | 3000 | 4200 |
/// |---|---|---|---|---|
/// | ms per vertex, prefix hit | 0.044 | 0.187 | 0.583 | 0.826 |
/// | ms per vertex, prefix miss | 0.247 | 1.285 | 3.937 | 5.375 |
///
/// The miss column agrees with the 0.264 row above to within about 5 % (0.247 at depth 190, ~0.250
/// scaled to 192), which is what says the two instruments measure the same thing — two different
/// rigs landing on one number, not a reproduction to the digit. But the hit is **not** constant —
/// it grows linearly too, at about a sixth of the slope, so the 0.3 ms budget merely moves from
/// depth ~192 to ~1540. And the guard runs on the chain's depth **before** anything knows whether
/// a prefix will be there to hit, so it must bound the road it cannot rule out — the miss, which
/// the accelerator does not touch. ⇒ the derivation stands and the number does not move. What the
/// accelerator changes is the common case, not the worst one.
///
/// ★ **It bounds a replay, not a depth.** A vertex whose chain's fold answers
/// ([`read_without_replay`]) is read at any depth: the fold was composed as each node was born, so
/// reading it costs the same under two nodes and under two thousand. Measured on 2000 recorded
/// quarter turns (release, same session): the build goes 65 → 100 ms, and every vertex past the
/// cap is its exact point instead of the construction's figure.
const CACHE_REPLAY_COST_CAP: usize = 192;

/// **The cache's own road** — what [`Model::vertex_point`] holds for every vertex an operation
/// makes: the definition realized on the ladder's first two rungs (the second only where the
/// first does not decide) and read out by [`Realized::to_f64`], or `None` where neither names an
/// `f64` (a realization that declines by name, an interval too wide at 256 bits) — or where the
/// road did not even walk, because a replay is deeper than [`CACHE_REPLAY_COST_CAP`] is willing
/// to pay for.
///
/// Two rungs, deliberately: a coordinate a low rung decides is the same answer any higher rung
/// would give, so where this answers it agrees with [`realize_vertex`] at `NearestF64` bit for
/// bit; where it does not, the cost of climbing on every push is not paid, and the cache says so
/// by carrying the construction's figure instead. Public so an instrument
/// can ask the same question the push funnel asked and hold the cache to it.
///
/// ★ The `Err` says **which** of the two roads was not taken, and that is what lets the funnel
/// below sort a vertex into [`PointCache::Ceiling`] (ask again, pay more) or
/// [`PointCache::Unrealized`] (there is no road).
pub fn realize_cache(model: &Model, def: &Vertex) -> Result<([f64; 3], [Mag; 3]), CacheDecline> {
    realize_cache_tracked(model, def, &mut None)
}

/// **The realization [`build_three_plane`] reaches without replaying a chain** — an unmoved meet,
/// or one whose chain's fold answers; `None` where only a replay would. Asked of the fold itself,
/// not of whether the chain folds: a fold whose offset or point overflowed answers nothing, and a
/// replay would walk the whole history — the cost [`CACHE_REPLAY_COST_CAP`] exists to bound.
fn read_without_replay(model: &Model, def: &Vertex) -> Option<Realized> {
    let nacre_topo::Vertex::ThreePlane(_) = def else {
        return None;
    };
    let (meet, frame) = model.vertex_meet_of(def)?;
    let (n, d) = match frame {
        None => meet.lift(),
        Some(node) => MeetPoint::Narrow(model.chain_point_rat(node, *meet.narrow()?)?).lift(),
    };
    Some(Realized(Arm::Exact(n, d)))
}

/// [`realize_cache`] that also reports what the accelerator could keep.
///
/// ★ The public door stays narrow on purpose: an instrument asks it the same question the funnel
/// asks and holds the cache to the answer, and that question is about a coordinate, not about a
/// side table. Only the funnel — which holds `&mut Model` and therefore can act on it — takes the
/// wider one.
pub(crate) fn realize_cache_tracked(
    model: &Model,
    def: &Vertex,
    out: &mut Option<PrefixWrite>,
) -> Result<([f64; 3], [Mag; 3]), CacheDecline> {
    let deep = def
        .carriers()
        .filter_map(|h| model.plane_motion(h))
        .any(|leaf| model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP));
    let r = if deep {
        read_without_replay(model, def).ok_or(CacheDecline::CostCap)?
    } else {
        let first = realize_def_tracked(model, def, Precision::Bits(LADDER[0]), out)
            .map_err(CacheDecline::Cannot)?;
        // ★ **A second rung, and why 256.** A coordinate that is not `0` loses about a bit per
        // turn off the quarters, and the cost cap stops the replay at 192 nodes, so 256 bits keep
        // more than an `f64`'s 53 inside the cap. A coordinate that is exactly `0` needs its radius
        // under the coincidence limit (`2⁻¹⁸⁰`), which 256 bits reach through some sixty to eighty such
        // turns (about 76 at a bit a turn; a 37°-and-back fixture decides at 60 and not at 80);
        // deeper, it stays `Ceiling` until the paid door climbs. The second rung files nothing
        // (`&mut None`): the first rung's prefix is the one the next generation asks for.
        if first.is_exact() || first.to_f64().is_some() {
            first
        } else {
            realize_def_tracked(model, def, Precision::Bits(LADDER[1]), &mut None)
                .map_err(CacheDecline::Cannot)?
        }
    };
    r.to_f64().ok_or(CacheDecline::Cannot(match r.is_exact() {
        // An exact value no `f64` names: more bits are not the missing thing.
        true => RealizeError::Unrepresentable,
        // The interval was still too wide at the first rung: more bits are exactly the thing.
        false => RealizeError::Undecided,
    }))
}

/// Push a vertex whose cache is **realized from its definition** — the one road every vertex an
/// operation makes takes, so `Model::vertex_point` is the realization's memo from the moment the
/// vertex exists (design: *"the cache is what `realize` produced, never a second truth"*).
///
/// `fallback` is the coordinate the construction site computed, and it stands wherever
/// [`realize_cache`] declines — but the *variant* is this funnel's to choose, not the caller's
/// ([`point_cache`]).
pub(crate) fn push_vertex_realized(
    model: &mut Model,
    def: Vertex,
    fallback: PointCache,
    link: ChainLink,
) -> Handle<Vertex> {
    let mut prefix = None;
    let Ok(cache) = point_cache_tracked(
        model,
        &def,
        || Ok::<_, core::convert::Infallible>(fallback.coord()),
        &mut prefix,
    );
    // ★ **Take what was used and leave what was made** — the handover that keeps the table at one
    // live generation. A vertex minted here rather than carried forward files nothing: nobody will
    // ever ask for its prefix, and an entry no reader can hit is a leak with a slow fuse.
    if let (ChainLink::Extends, Some(w)) = (link, prefix) {
        model.hand_over_prefix_hp(w.used, w.key, w.value);
    }
    model.push_vertex(def, cache)
}

/// Push a vertex with the cache [`point_cache`] already chose **for this very definition** — the
/// funnel's answer, asked once by whoever needed the coordinate first (the seam table), so the
/// minting does not pay the realization a second time. Never a figure: a cache from anywhere else
/// goes through [`push_vertex_realized`].
pub(crate) fn push_vertex_asked(
    model: &mut Model,
    def: Vertex,
    cache: PointCache,
) -> Handle<Vertex> {
    model.push_vertex(def, cache)
}

/// **The cache a vertex with this definition receives** — asked before anything is pushed, by
/// whoever needs the coordinate first: the seam table realizes a node from the definition the
/// minting will push, so the arrangement's point and the model's cache are one answer.
///
/// `fallback` is asked only where the road declines, and its error is the caller's: a
/// construction figure the caller cannot produce either (the seam table's plane-cache solve on
/// planes whose caches are parallel) is a refusal of its own, not a cache.
///
/// ★ **An exhaustive `match`, so a new reason cannot be folded silently into an old name.** The
/// two reasons that mean "ask again and pay more" become [`PointCache::Ceiling`]; everything else
/// means the kernel has no road, and becomes [`PointCache::Unrealized`]. That split exists for the
/// reader: "this model is expensive" and "the kernel cannot do this" are different reports, and a
/// single label would hide the second behind the first.
pub(crate) fn point_cache<E>(
    model: &Model,
    def: &Vertex,
    fallback: impl FnOnce() -> Result<Point3, E>,
) -> Result<PointCache, E> {
    point_cache_tracked(model, def, fallback, &mut None)
}

/// [`point_cache`] that also reports what the accelerator could keep — the funnel's form.
fn point_cache_tracked<E>(
    model: &Model,
    def: &Vertex,
    fallback: impl FnOnce() -> Result<Point3, E>,
    out: &mut Option<PrefixWrite>,
) -> Result<PointCache, E> {
    Ok(match realize_cache_tracked(model, def, out) {
        Ok((coord, bound)) => PointCache::Bounded {
            coord: Point3::from_array(coord),
            bound,
        },
        // Bits or cost — either way a paid realization can still answer.
        Err(CacheDecline::CostCap | CacheDecline::Cannot(RealizeError::Undecided)) => {
            PointCache::Ceiling { coord: fallback()? }
        }
        Err(CacheDecline::Cannot(
            RealizeError::NoMeet
            | RealizeError::WideUnderMotion
            | RealizeError::NoMotionChain
            | RealizeError::NoCurvedPoint
            | RealizeError::Unrepresentable,
        )) => PointCache::Unrealized { coord: fallback()? },
    })
}

/// Push a plane whose cache is **realized from its truth** where the model cannot derive one — the
/// plane twin of [`push_vertex_realized`], and the road every constructed plane takes.
///
/// The push door derives a plane's cache itself wherever the plane has a world name (unmoved, or
/// carried by a chain that folds). What it cannot name is a plane under a turn off the quarters or
/// a frame, and there `figure` — the construction's own `f64` — would stand. This realizes that
/// plane's **normal** instead: the pre-motion name's normal carried by the chain (the difference
/// of the replayed origin and name vector — the translation cancels) times the name's sense
/// against the points and the chain's parity, normalized, and read out by [`Realized::to_f64`] at
/// 128 bits and then 256, so a component that is exactly `0` is `+0.0`. Where that does not
/// answer — a replay past [`CACHE_REPLAY_COST_CAP`], a pre-motion name wider than `Rat`, undecided
/// at 256 — `figure` stands, as a vertex's construction figure does.
///
/// ⚠ **The anchor stays `figure`'s.** Realizing it too (the first point replayed, read out the
/// same way) is the truer cache, and it cost a cell: the seam table solves the seam vertices the
/// realization has no road to from the three classes' plane caches (todo 「seam 표의 폴백과 tol 은 평면 캐시를 읽는다」), and the moved
/// `d` made two seam points alias — `collinear_loop_points`' `Y/305deg/inset0.5` refused `SeamAlias`
/// (107 of 108 cells built). The normal alone moved none.
///
/// ★ Not the rejected road (design 「가지 말 것」: realizing such a plane by replaying its three
/// points in `f64`, which was measured further from the truth than the producer's figure): this
/// realizes at 128 or 256 bits and rounds once.
///
/// ⚠ «Has a world name» is asked before the push, so without a handle: no motion, or a chain whose
/// fold answers. A `Wide` pre-motion name under a folding chain is named by neither side and keeps
/// `figure`.
pub(crate) fn push_plane_realized(
    model: &mut Model,
    figure: nacre_geom::Plane,
    points: [[nacre_exact::Rat; 3]; 3],
    motion: Option<Handle<nacre_topo::MotionNode>>,
    sense: nacre_topo::Orientation,
) -> (Handle<Surface>, bool) {
    let normal = motion.and_then(|leaf| realize_plane_normal(model, &points, leaf, sense));
    let cache = normal.map_or(figure, |n| {
        nacre_geom::Plane::from_point_unit_normal(
            figure.origin(),
            nacre_math::Vector3::from_array(n),
        )
    });
    model.push_plane(cache, points, motion, sense)
}

/// The normal [`push_plane_realized`] describes, `None` where it does not answer.
fn realize_plane_normal(
    model: &Model,
    points: &[[nacre_exact::Rat; 3]; 3],
    leaf: Handle<nacre_topo::MotionNode>,
    sense: nacre_topo::Orientation,
) -> Option<[f64; 3]> {
    let zero = [nacre_exact::Rat::from_int(0); 3];
    if model.chain_point_rat(leaf, zero).is_some()
        || model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP)
    {
        return None;
    }
    let name = nacre_exact::plane_name_exact(points[0], points[1], points[2])?;
    let n = *name.narrow()?;
    let along = nacre_exact::name_along_points(
        &name,
        points.each_ref().map(|p| MeetPoint::Narrow(*p)).each_ref(),
    )?;
    // The way the plane faces, in the frame its points are written in: the name's normal when it
    // runs with the points' turn times the sense, and a reflection in the chain turns it again.
    let facing = if along { sense.sign() } else { -sense.sign() };
    let turned = facing * crate::rotated_vertex::motion_parity(model, Some(leaf))? < 0;
    let chain = motion_chain(model, leaf)?;
    let tip = [n[0], n[1], n[2]];
    // `tail → head` is the carried normal; swapping them negates it exactly.
    let (tail, head) = if turned { (tip, zero) } else { (zero, tip) };
    let tail = replay(WitnessPoint::at(tail), &chain)?;
    let head = replay(WitnessPoint::at(head), &chain)?;
    for bits in [LADDER[0], LADDER[1]] {
        let (t, h) = (tail.realize(bits), head.realize(bits));
        let d: [HpBounded; 3] = core::array::from_fn(|k| h[k].sub(&t[k], bits));
        let nn = d[0]
            .mul(&d[0], bits)
            .add(&d[1].mul(&d[1], bits), bits)
            .add(&d[2].mul(&d[2], bits), bits);
        let inv = nn.inv_sqrt(bits)?;
        let unit: [HpBounded; 3] = core::array::from_fn(|k| d[k].mul(&inv, bits));
        if let Some((n, _)) = Realized(Arm::Approached(unit, bits)).to_f64() {
            return Some(n);
        }
    }
    None
}

/// What [`refine_vertex_cache`] did: how many coordinates it raised, and how many it could not.
///
/// ★ `left` is not an error count. A vertex stays behind when the full ladder still does not name
/// an `f64` for it — the only honest thing to report, and the number a caller watches if it wants
/// to know whether the model has coordinates no precision will settle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RefineReport {
    pub refined: usize,
    pub left: usize,
}

/// **Pay for the realizations the cache road would not** — raise every [`PointCache::Ceiling`] in
/// the live model to [`PointCache::Bounded`], climbing the whole ladder for each.
///
/// The cache road takes one rung and refuses a history past its cost cap, because it runs on
/// *every* push. This door runs when a caller asks, so it can afford what that one cannot: a
/// cost-capped chain is walked, and an undecided coordinate is climbed until it decides.
///
/// Only `Ceiling` is touched. `Bounded` is already the realization, and `Unrealized` means there is
/// no road to one — re-trying those would conflate "this model is expensive" with "the kernel
/// cannot do this", which is the distinction the two variants exist to keep.
///
/// ★★ **It changes caches, not truths — so it may run mid-log.** An operation after it derives its
/// own caches from the refined ones (a moved vertex's fallback figure, a datum's cache anchor), so
/// the caches differ from the unrefined road's; what becomes truth does not, because no operation
/// decides a truth from a vertex cache: `push_plane_through` reads the caches only to point the
/// plane's cache (its sense is the permutation's parity), a frame's `flip` comes from the truth
/// (`frame_toward`), and whether a motion is recorded is asked of the statements. Measured on a
/// log that moves, fuses, states a datum through vertices and pads after the door
/// (`refining_mid_log_leaves_every_later_truth_as_it_was`, in `tests/invariants/replay.rs`).
/// It is an export-time door because that is where the precision is wanted.
///
/// Idempotent: a second call finds nothing to raise. Edge curves are re-derived once at the end,
/// and only if something actually moved.
pub fn refine_vertex_cache(model: &mut Model) -> RefineReport {
    let todo: Vec<Handle<Vertex>> = model
        .reachable()
        .vertices
        .into_iter()
        .filter(|&vh| matches!(model.vertex_cache(vh), PointCache::Ceiling { .. }))
        .collect();
    let mut out = RefineReport {
        refined: 0,
        left: 0,
    };
    for vh in todo {
        // ⚠ `Ok` is not enough on its own: `climb` returns early for an exact arm without asking
        // whether an `f64` names it, so the answer can still be unrepresentable. That value cannot
        // be a `Ceiling` (it is classified `Unrealized`), so this branch should not arrive — but
        // that is an agreement between two functions, not something the types promise, and an
        // `expect` here would turn the disagreement into a panic instead of a count.
        match realize_vertex(model, vh, Precision::NearestF64).map(|r| r.to_f64()) {
            Ok(Some((coord, bound))) => {
                model.refine_vertex_cache(vh, Point3::from_array(coord), bound);
                out.refined += 1;
            }
            _ => out.left += 1,
        }
    }
    if out.refined > 0 {
        model.rebuild_edge_cache();
    }
    out
}

/// `places` decimal places, escalating until the realization determines them.
///
/// ★ This is the door the app and the kit call, and its contract is **"ask for `places`, get
/// `places`"**. The `None` on [`Realized::to_decimal`] is the rung-to-rung signal underneath, not
/// something a caller sees: the only ways this fails are a vertex that cannot be realized at all
/// and a ladder that ran out.
pub fn realize_vertex_decimal(
    model: &Model,
    v: Handle<Vertex>,
    places: usize,
) -> Result<[String; 3], RealizeError> {
    let out = climb(model, model.vertex(v), |r| r.to_decimal(places).map(|_| r))?;
    out.to_decimal(places).ok_or(RealizeError::Undecided)
}

/// Walk [`LADDER`] until `decided` accepts a realization. The exact arm decides on the first rung
/// whatever is asked, so a rational vertex never climbs.
fn climb(
    model: &Model,
    def: &Vertex,
    decided: impl Fn(Realized) -> Option<Realized>,
) -> Result<Realized, RealizeError> {
    for bits in LADDER {
        match build(model, def, bits, &mut None) {
            Ok(r) => {
                if r.is_exact() {
                    return Ok(r);
                }
                if let Some(r) = decided(r) {
                    return Ok(r);
                }
            }
            // A structural refusal does not improve with precision.
            Err(e) => return Err(e),
        }
    }
    // Every rung built, none decided — the only thing running out of ladder can mean.
    Err(RealizeError::Undecided)
}

/// One realization at `bits`, from the definition.
/// `out` is the accelerator's channel and only the three-plane road fills it — a curved definition
/// realizes from its own geometry, not by walking a chain, so it has no prefix to hand on.
fn build(
    model: &Model,
    def: &Vertex,
    bits: usize,
    out: &mut Option<PrefixWrite>,
) -> Result<Realized, RealizeError> {
    match *def {
        nacre_topo::Vertex::ThreePlane(_) => build_three_plane(model, def, bits, out),
        nacre_topo::Vertex::OnSeam([cyl, cap]) => curved(seam_point(model, cyl, cap, bits), bits),
        nacre_topo::Vertex::Pierce {
            planes,
            cylinder,
            root,
        } => curved(pierce_point(model, planes, cylinder, root, bits), bits),
    }
}

fn curved(p: Option<[HpBounded; 3]>, bits: usize) -> Result<Realized, RealizeError> {
    Ok(Realized(Arm::Approached(
        p.ok_or(RealizeError::NoCurvedPoint)?,
        bits,
    )))
}

/// **The seam vertex is the rim's `+ref_dir` point** — `centre + r·ê`, `ê` the reference
/// direction's unit part perpendicular to the axis. `Vertex::OnSeam`'s doc calls it *"a unique
/// point, exactly designated"*; this realizes that designation instead of reading the cache, which
/// is the half of migration row 3b that was recorded as missing.
fn seam_point(
    model: &Model,
    cyl: Handle<Surface>,
    cap: Handle<Surface>,
    bits: usize,
) -> Option<[HpBounded; 3]> {
    let def = model.world_cylinder_def(cyl)?;
    let coeffs = crate::planes::world_plane_coeffs(model, cap)?;
    let (o, m, r2, e) = (def.origin(), def.dir(), def.r2(), def.ref_dir());
    let t = crate::planes::axis_param_of_plane(&coeffs, &def)?;
    let mut centre = o;
    for k in 0..3 {
        centre[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
    }
    nacre_exact::realize_seam_point(centre, perp_component(&e, &m)?, r2, bits)
}

/// **A pierce vertex is the meet line's point at its root** — the two cutting planes give the
/// line, the cylinder gives the quadratic, and `root` names which crossing.
fn pierce_point(
    model: &Model,
    planes: [Handle<Surface>; 2],
    cylinder: Handle<Surface>,
    root: nacre_topo::QuadRoot,
    bits: usize,
) -> Option<[HpBounded; 3]> {
    let def = model.world_cylinder_def(cylinder)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let c0 = crate::planes::world_plane_coeffs(model, planes[0])?;
    let c1 = crate::planes::world_plane_coeffs(model, planes[1])?;
    let (line, sv) = pick_root(
        &nacre_exact::quad::plane_plane_cylinder(&c0, &c1, &o, &m, r2)?,
        root,
    )?;
    // `point = base + s·dir`, with `s` the quadratic root realized at `bits` — the only step that
    // is not exact rational arithmetic.
    let sb = nacre_exact::realize_quad(&sv, bits)?;
    let (b, d) = (line.base(), line.dir());
    let coord = |k: usize| nacre_exact::affine_bounded(b[k], d[k], &sb, bits);
    Some([coord(0)?, coord(1)?, coord(2)?])
}

fn build_three_plane(
    model: &Model,
    def: &Vertex,
    bits: usize,
    out: &mut Option<PrefixWrite>,
) -> Result<Realized, RealizeError> {
    let (meet, frame) = model.vertex_meet_of(def).ok_or(RealizeError::NoMeet)?;
    let Some(node) = frame else {
        // No motion: the meet *is* the coordinate, and `lift` states it as integers whatever its
        // width — so `Wide` is not a refusal on this road.
        let (n, d) = meet.lift();
        return Ok(Realized(Arm::Exact(n, d)));
    };
    let base = *narrow_or(&meet)?;
    // ★ **A chain that folds is read, not replayed.** Translations, quarter turns and axis
    // reflections compose to a signed axis permutation plus a rational offset, folded once per node
    // as it is born (`Model::chain_point_rat`), so the meet carried through it is an exact rational
    // — the same point the replay approaches, stated outright. A replay cannot decide a coordinate
    // that is exactly `0` (its interval straddles it at every rung), and read from the fold it is
    // the integer `0`. Where the fold does not answer (a frame, a turn off the quarters, an offset
    // past `Rat`) the replay does.
    if let Some(world) = model.chain_point_rat(node, base) {
        let (n, d) = MeetPoint::Narrow(world).lift();
        return Ok(Realized(Arm::Exact(n, d)));
    }
    let chain = motion_chain(model, node).ok_or(RealizeError::NoMotionChain)?;
    let (point, used) = match remembered_prefix(model, base, node, bits) {
        // The chain runs root-to-leaf, so what this point still owes is the tail past the prefix.
        Some((used, folded, prefix)) => (
            nacre_judge::fold_suffix(prefix, &chain[folded..], bits),
            Some(used),
        ),
        None => (
            replay(WitnessPoint::at(base), &chain)
                .ok_or(RealizeError::NoMotionChain)?
                .realize(bits),
            None,
        ),
    };
    *out = Some(PrefixWrite {
        used,
        key: (base, node, bits),
        value: (chain.len(), point.clone()),
    });
    Ok(Realized(Arm::Approached(point, bits)))
}

/// **The deepest remembered prefix of this vertex's chain, within two motion nodes of its leaf.**
///
/// Returns the key it hit (so the caller that owns `&mut Model` can consume it), how many chain
/// nodes that value already folded, and the value itself.
///
/// ⚠ **Two steps, and the bound is the producer's, not a round number.** `chain_motion` in
/// `transform.rs` records at most two nodes for one call — a rotation and a translation — so two
/// steps always reach the previous generation and a third could only reach a generation that was
/// already consumed. Without a bound a miss would walk the whole history doing a failed lookup per
/// node, and the paid door would pay that on every rung of the ladder.
fn remembered_prefix(
    model: &Model,
    base: [nacre_exact::Rat; 3],
    leaf: Handle<nacre_topo::MotionNode>,
    bits: usize,
) -> Option<(PrefixKey, usize, [HpBounded; 3])> {
    let mut anc = model.motion(leaf).parent;
    for _ in 0..2 {
        let a = anc?;
        if let Some((folded, p)) = model.prefix_hp(base, a, bits) {
            return Some(((base, a, bits), *folded, p.clone()));
        }
        anc = model.motion(a).parent;
    }
    None
}

/// What a realization learned that the accelerator could keep: the entry it consumed (if any) and
/// the one it produced. Filled on the way down, acted on by whoever holds `&mut Model`.
pub(crate) struct PrefixWrite {
    pub(crate) used: Option<PrefixKey>,
    pub(crate) key: PrefixKey,
    pub(crate) value: (usize, [HpBounded; 3]),
}

/// **Whether this vertex continues a chain the accelerator is already following.**
///
/// ⚠ Not a `bool`. At a call site `true` says nothing about what is being claimed, and the claim
/// here is specific: *the vertex being pushed is the image of one realized before, so its prefix is
/// worth keeping for the next motion*. Only `transform` can say that. An arrangement's result
/// vertex is nobody's prefix — remembering it would grow the table by an entry per boolean, for a
/// lookup that can never hit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChainLink {
    /// The image of a vertex that came before — file its prefix.
    Extends,
    /// Minted here. Read the table if it helps, but do not write to it.
    Fresh,
}

fn narrow_or(meet: &MeetPoint) -> Result<&[nacre_exact::Rat; 3], RealizeError> {
    meet.narrow().ok_or(RealizeError::WideUnderMotion)
}

/// **`e₁ = (m·m)e − (e·m)m`** — the part of `e` perpendicular to `m`, unnormalized (the caller
/// divides by its length, which is the one radical).
///
/// ⚠ **A no-op for every cylinder the kit builds today**, whose `ref_dir` is already perpendicular
/// to the axis — so no fixture can tell this from `mm·e`, and planting that very mutation left the
/// corpus green. It is kept, and unit-tested directly, because `CylinderDef::new` only requires
/// `ref_dir × dir ≠ 0`: a slanted reference direction is a *representable* statement, and this is
/// the step that stops it landing off the rim.
fn perp_component(
    e: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
) -> Option<[nacre_exact::Rat; 3]> {
    let (mm, em) = (nacre_exact::dot3_rat(m, m)?, nacre_exact::dot3_rat(e, m)?);
    let mut out = [nacre_exact::Rat::from_int(0); 3];
    for k in 0..3 {
        out[k] = mm.checked_mul(e[k])?.checked_sub(em.checked_mul(m[k])?)?;
    }
    Some(out)
}

/// **Which of the crossings `root` names.** `Lo`/`Hi` are ascending parameter along the meet
/// line's direction — `nacre_exact::quad`'s pair order, which is `Vertex::Pierce`'s stated
/// convention — and a tangency is one point spelled `Double`, never `Lo`.
///
/// ⚠ Split out because the integration oracle ("the point is on the cylinder") is satisfied by
/// **both** roots: swapping them left the whole corpus green when planted. This is the piece a
/// test can ask the distinguishing question of.
fn pick_root(
    meet: &nacre_exact::quad::CylinderMeet,
    root: nacre_topo::QuadRoot,
) -> Option<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)> {
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    match meet {
        CylinderMeet::Pair { line, s } => match root {
            QuadRoot::Lo => Some((line.clone(), s[0])),
            QuadRoot::Hi => Some((line.clone(), s[1])),
            QuadRoot::Double => None,
        },
        CylinderMeet::Tangent { line, s } => match root {
            QuadRoot::Double => Some((line.clone(), QuadVal::from_rat(*s))),
            _ => None,
        },
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_exact::Rat;

    fn r(n: i128) -> Rat {
        Rat::from_int(n)
    }

    /// A reference direction that leans along the axis must lose exactly that lean.
    #[test]
    fn the_perpendicular_component_removes_the_axial_lean() {
        let axis = [r(0), r(0), r(1)];
        // Already perpendicular: unchanged up to the `m·m` scale (1 here).
        assert_eq!(
            perp_component(&[r(1), r(0), r(0)], &axis).unwrap(),
            [r(1), r(0), r(0)]
        );
        // Leaning: `(1,0,5)` against `+z` keeps only its `x`.
        assert_eq!(
            perp_component(&[r(1), r(0), r(5)], &axis).unwrap(),
            [r(1), r(0), r(0)]
        );
        // A non-unit axis scales but still projects: `(0,0,2)` gives `m·m = 4`.
        assert_eq!(
            perp_component(&[r(3), r(0), r(7)], &[r(0), r(0), r(2)]).unwrap(),
            [r(12), r(0), r(0)]
        );
    }

    /// ★★ **`Lo` is the lower parameter along the meet line, and `Hi` the higher.**
    ///
    /// The oracle the integration test uses — "the point is on the cylinder" — is satisfied by
    /// *both* roots, so swapping them leaves it green (planted and measured). This asks the
    /// question that actually distinguishes them, in the vocabulary `Vertex::Pierce`'s doc
    /// defines: ascending parameter along `n₀ × n₁`.
    #[test]
    fn lo_and_hi_run_along_the_line() {
        use nacre_exact::quad::CylinderMeet;
        use nacre_topo::QuadRoot;
        let (o, m, rad) = ([r(2), r(2), r(0)], [r(0), r(0), r(1)], r(9)); // radius 3, as r²
        // x = 0 and z = 0: the meet line runs along +y and crosses the cylinder twice.
        let c0 = [r(1), r(0), r(0), r(0)];
        let c1 = [r(0), r(0), r(1), r(0)];
        let meet = nacre_exact::quad::plane_plane_cylinder(
            &c0,
            &c1,
            &o,
            &m,
            &nacre_exact::BigRat::from(rad),
        )
        .expect("a crossing");
        assert!(
            matches!(meet, CylinderMeet::Pair { .. }),
            "expected two roots"
        );
        let (line, lo_s) = pick_root(&meet, QuadRoot::Lo).expect("Lo");
        let (_, hi_s) = pick_root(&meet, QuadRoot::Hi).expect("Hi");
        assert!(
            pick_root(&meet, QuadRoot::Double).is_none(),
            "a pair is not a Double"
        );
        let at = |q: &nacre_exact::quad::QuadVal| {
            let sb = nacre_exact::realize_quad(q, 256).expect("realized");
            let (b, d) = (line.base(), line.dir());
            let c = |k: usize| {
                let s = nacre_exact::affine_bounded(b[k], d[k], &sb, 256).expect("affine");
                nacre_exact::round_to_digits(&s.value, s.error, 20)
                    .expect("decided")
                    .parse::<f64>()
                    .expect("a decimal")
            };
            [c(0), c(1), c(2)]
        };
        let (lo, hi) = (at(&lo_s), at(&hi_s));
        let dir = line.dir().map(|c| c.to_f64());
        let proj = |p: [f64; 3]| p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2];
        assert!(
            proj(lo) < proj(hi),
            "Lo must precede Hi along the line: {lo:?} then {hi:?}"
        );
        // And both must sit on the cylinder — the property the swap cannot break.
        for p in [lo, hi] {
            let rr = ((p[0] - 2.0).powi(2) + (p[1] - 2.0).powi(2)).sqrt();
            assert!((rr - 3.0).abs() < 1e-9, "{p:?} is {rr} from the axis");
        }
    }
}
