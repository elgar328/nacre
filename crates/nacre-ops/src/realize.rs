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
use nacre_topo::{Edge, EdgeGiven, Model, PointCache, PrefixKey, Surface, Vertex};
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
    /// A three-plane vertex no road reaches: its carriers share no frame the exact roads read
    /// ([`nacre_topo::Model::vertex_meet`] declines, or the meet is wider than `Rat` under a
    /// motion), **and** one of them has no witness triangle of its own for the mixed road
    /// (`rotated_vertex::surface_witness_triangle` — a nameless `Through` carrier, a chain outside
    /// the decimal window).
    NoMeet,
    /// The vertex's motion chain could not be rebuilt exactly.
    NoMotionChain,
    /// A curved definition (`OnSeam`, `Pierce`) did not resolve into a point.
    ///
    /// ⚠ **This is a bag, and saying so is the point.** It covers: a pierce corner whose carriers
    /// cannot be stated in the world exactly (a seam corner is met in the world instead), a chain
    /// a seam corner's carriers cannot replay, a cap plane parallel to the axis (so it bounds no
    /// rim), a `root` that does not match the kind of crossing the carriers actually
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
    realize_def_tracked(model, def, p, &mut Accel::default())
}

/// [`realize_def`] with the accelerators ([`Accel`]): what it may read to go faster, and what it
/// learned that they could keep.
///
/// ★ Only the fixed-precision road fills the prefix. A climb asks the same definition at rung after
/// rung, so whatever it learned at 128 bits is not what it returns, and filing the rung it happened
/// to pass through would put a value under a key the next reader asks at a different precision.
pub(crate) fn realize_def_tracked(
    model: &Model,
    def: &Vertex,
    p: Precision,
    acc: &mut Accel<'_>,
) -> Result<Realized, RealizeError> {
    match p {
        Precision::Bits(bits) => build(model, def, bits, acc),
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
    realize_cache_tracked(model, def, &mut Accel::default())
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
    acc: &mut Accel<'_>,
) -> Result<([f64; 3], [Mag; 3]), CacheDecline> {
    let deep = def
        .carriers()
        .filter_map(|h| model.plane_motion(h))
        .any(|leaf| model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP));
    // ★ **For a vertex, the cap bounds the nodes a replay walks, whichever road walks them.** The
    // shared-frame road replays one point through one chain; the mixed road ([`build_meet`])
    // replays three witness points through each carrier's chain, so under the same budget it
    // counts their sum. (The plane funnel guards the chain's depth instead — `push_plane_realized`
    // — and so does a seam met in the world, [`seam_point_met`]: three points through the
    // cylinder's chain and, for a cap without a world name, three through the cap's; the
    // per-carrier depth test above is its guard, `deep` sending it to `CostCap`.)
    if !deep && meet_road_over_cap(model, def) {
        return Err(CacheDecline::CostCap);
    }
    let r = if deep {
        read_without_replay(model, def).ok_or(CacheDecline::CostCap)?
    } else {
        // A first rung that could not decide (`Undecided` — a seam met after a long chain) is the
        // second rung's business, as a first rung whose interval names no `f64` is.
        let first = match realize_def_tracked(model, def, Precision::Bits(LADDER[0]), acc) {
            Err(RealizeError::Undecided) => None,
            other => Some(other.map_err(CacheDecline::Cannot)?),
        };
        // ★ **A second rung, and why 256.** A coordinate that is not `0` loses about a bit per
        // turn off the quarters, and the cost cap stops the replay at 192 nodes, so 256 bits keep
        // more than an `f64`'s 53 inside the cap. A coordinate that is exactly `0` needs its radius
        // under the coincidence limit (`2⁻¹⁸⁰`), which 256 bits reach through some sixty to eighty such
        // turns (about 76 at a bit a turn; a 37°-and-back fixture decides at 60 and not at 80);
        // deeper, it stays `Ceiling` until the paid door climbs. The second rung files nothing
        // (`&mut None`): the first rung's prefix is the one the next generation asks for.
        if let Some(first) = first.filter(|f| f.is_exact() || f.to_f64().is_some()) {
            first
        } else {
            let mut second = Accel {
                prefix: None,
                planes: acc.planes.as_deref_mut(),
            };
            realize_def_tracked(model, def, Precision::Bits(LADDER[1]), &mut second)
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
    let mut acc = Accel::default();
    let Ok(cache) = point_cache_tracked(
        model,
        &def,
        || Ok::<_, core::convert::Infallible>(fallback.coord()),
        &mut acc,
    );
    // ★ **Take what was used and leave what was made** — the handover that keeps the table at one
    // live generation. A vertex minted here rather than carried forward files nothing: nobody will
    // ever ask for its prefix, and an entry no reader can hit is a leak with a slow fuse.
    if let (ChainLink::Extends, Some(w)) = (link, acc.prefix) {
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
/// `planes` is the batch's plane memo ([`PlaneMemo`]) — a caller asking this of many vertices on
/// one model holds one, so a carrier plane shared by many corners is realized once.
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
    planes: &mut PlaneMemo,
    fallback: impl FnOnce() -> Result<Point3, E>,
) -> Result<PointCache, E> {
    let mut acc = Accel {
        prefix: None,
        planes: Some(planes),
    };
    point_cache_tracked(model, def, fallback, &mut acc)
}

/// [`point_cache`] with the accelerators ([`Accel`]) — the funnel's form.
fn point_cache_tracked<E>(
    model: &Model,
    def: &Vertex,
    fallback: impl FnOnce() -> Result<Point3, E>,
    acc: &mut Accel<'_>,
) -> Result<PointCache, E> {
    Ok(match realize_cache_tracked(model, def, acc) {
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
/// plane's cache instead, each half on its own:
/// * **Anchor** — the truth's first point replayed through the chain: the same point the door
///   anchors a named plane at, so the two kinds of plane agree about where a cache is pinned.
/// * **Normal** — the pre-motion name's normal carried by the chain (the difference of the
///   replayed origin and name vector — the translation cancels) times the name's sense against
///   the points and the chain's parity, normalized. A pre-motion name wider than `Rat` has no
///   vector to replay, so its normal is the one the three replayed points span
///   ([`nacre_judge::plane_hp`]) — the same plane, and where both roads decide the same bits
///   (measured on 21,543 narrow names); not the road for every plane, because at 256 bits it
///   leaves undecided 226 normals the name road decides.
///
/// Both are read out by [`Realized::to_f64`] at 128 bits and then 256, so a component that is
/// exactly `0` is `+0.0`. Where a half does not answer, `figure`'s half stands, as a vertex's
/// construction figure does: a replay past [`CACHE_REPLAY_COST_CAP`] (both halves), undecided at
/// 256 (that half). The door does the same for a
/// named plane whose truth cannot place its first point — a derived normal beside the producer's
/// anchor — so a cache with one realized half is not a new shape.
///
/// ★ **What the anchor buys** (measured over the suite, 21,769 planes whose normal this realizes):
/// the producer's anchor is more than an ulp off the true plane for 6,933 of them and up to 37 ulps
/// — a moved plane's figure is the previous cache moved in `f64`, so the rounding accumulates with
/// each generation — and the realized one is within an ulp for all of them.
///
/// ★ Not the rejected road (design 「가지 말 것」: realizing such a plane by replaying its three
/// points in `f64`, which was measured further from the truth than the producer's figure): this
/// realizes at 128 or 256 bits and rounds once.
///
/// ⚠ **The guard is the chain's depth, not the nodes replayed.** Three points walk the chain (the
/// anchor, the origin, the name vector — or the three stated points), so a plane under 192 nodes
/// costs up to three such
/// replays — where the vertex road's mixed arm counts the sum against the same cap
/// ([`meet_road_over_cap`]). Counting the sum here would send the 5,976 planes the suite realizes
/// at depths 65–192 back to `figure` (measured).
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
    let (anchor, normal) = match motion {
        Some(leaf) => realize_plane_cache(model, &points, leaf, sense),
        None => (None, None),
    };
    let cache = nacre_geom::Plane::from_point_unit_normal(
        anchor.map_or(figure.origin(), Point3::from_array),
        normal.map_or(figure.normal(), nacre_math::Vector3::from_array),
    );
    model.push_plane(cache, points, motion, sense)
}

/// Push a cylinder whose cache is **realized from its truth** where the model cannot derive one —
/// the cylinder twin of [`push_plane_realized`], and the road every constructed cylinder takes.
///
/// The push door derives a cylinder's cache itself wherever the cylinder is stated in the world
/// (unmoved, or carried by a chain that folds). What it cannot state is a cylinder under a turn
/// off the quarters or a frame, and there `figure` — the construction's own `f64` — would stand.
/// This realizes that cache instead from the chain ([`world_cylinder_hp`]): the origin, the unit
/// axis and the unit `ref_dir` replayed at 128 or 256 bits and each rounded once; the radius is
/// the statement's (the chain is an isometry). The same guard as the plane's: a chain deeper than
/// [`CACHE_REPLAY_COST_CAP`] keeps `figure`, as does one the two rungs leave undecided.
pub(crate) fn push_cylinder_realized(
    model: &mut Model,
    figure: nacre_geom::Cylinder,
    def: nacre_topo::CylinderDef,
    motion: Option<Handle<nacre_topo::MotionNode>>,
) -> Handle<Surface> {
    let zero = [nacre_exact::Rat::from_int(0); 3];
    let realized = motion
        .filter(|&leaf| {
            model.chain_point_rat(leaf, zero).is_none()
                && !model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP)
        })
        .and_then(|_| realize_cylinder_cache(model, &def, motion));
    model.push_cylinder(realized.unwrap_or(figure), def, motion)
}

/// The cache [`push_cylinder_realized`] describes, on the cache road's two rungs — `None` where
/// the chain does not answer or neither rung decides every component.
fn realize_cylinder_cache(
    model: &Model,
    def: &nacre_topo::CylinderDef,
    motion: Option<Handle<nacre_topo::MotionNode>>,
) -> Option<nacre_geom::Cylinder> {
    let read = |v: [HpBounded; 3], bits: usize| {
        Realized(Arm::Approached(v, bits)).to_f64().map(|(v, _)| v)
    };
    [LADDER[0], LADDER[1]].into_iter().find_map(|bits| {
        let hp = world_cylinder_hp(model, def, motion, bits)?;
        let origin = read(hp.origin, bits)?;
        let axis = read(unit(hp.axis, bits)?, bits)?;
        let ref_dir = read(unit(hp.seam, bits)?, bits)?;
        nacre_geom::Cylinder::from_unit_frame(
            Point3::from_array(origin),
            nacre_math::Vector3::from_array(axis),
            nacre_math::Vector3::from_array(ref_dir),
            def.radius_f64(),
        )
    })
}

/// **A cylinder in the world, carried by its motion chain at `bits`** — what both the cylinder
/// cache ([`push_cylinder_realized`]) and the seam vertex of a moved cylinder read.
#[derive(Clone)]
pub(crate) struct CylinderHp {
    /// The statement's origin, replayed.
    pub(crate) origin: [HpBounded; 3],
    /// The axis as replayed: `R(o + dir) − R(o)`, of the statement's length `|dir|`.
    pub(crate) axis: [HpBounded; 3],
    /// The seam direction scaled to the radius: `r·ê` for `ê` the unit part of `ref_dir` across
    /// the axis — so `origin + seam` is the seam's point on the rim through the origin.
    pub(crate) seam: [HpBounded; 3],
}

/// [`CylinderHp`] for the statement `def` under `motion` (`None`: the statement is the world).
///
/// ★ **A direction is the difference of two replayed points.** The chain is affine — turns,
/// translations, reflections, frames — so the image of a direction `v` is `R(o + v) − R(o)`: the
/// translation cancels, which is how the plane funnel carries a name's normal
/// ([`carried_name_normal`]). Replaying `v` as if it were a point would carry the translation
/// into the direction.
///
/// `None` where the chain cannot be read, the statement's arithmetic leaves `Rat` (the part of
/// `ref_dir` across the axis, `e₁ = (dir·dir)·ref_dir − (ref_dir·dir)·dir`), or `r²` is wider than
/// `Rat` — the caller keeps what it had.
pub(crate) fn world_cylinder_hp(
    model: &Model,
    def: &nacre_topo::CylinderDef,
    motion: Option<Handle<nacre_topo::MotionNode>>,
    bits: usize,
) -> Option<CylinderHp> {
    use nacre_exact::{Rat, dot3_rat};
    let chain = match motion {
        Some(leaf) => motion_chain(model, leaf)?,
        None => Vec::new(),
    };
    let (o, m) = (def.origin(), def.dir());
    let e1 = perp_component(&def.ref_dir(), &m)?;
    let plus = |v: &[Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            o[0].checked_add(v[0])?,
            o[1].checked_add(v[1])?,
            o[2].checked_add(v[2])?,
        ])
    };
    let at = |p: [Rat; 3]| replay(WitnessPoint::at(p), &chain).map(|w| w.realize(bits));
    let origin = at(o)?;
    let towards = |p: [HpBounded; 3]| -> [HpBounded; 3] {
        core::array::from_fn(|k| p[k].sub(&origin[k], bits))
    };
    let axis = towards(at(plus(&m)?)?);
    let across = towards(at(plus(&e1)?)?);
    // `r / |e₁|` read as `√(r² / e₁·e₁)` — one rational under one root.
    let ee = dot3_rat(&e1, &e1)?;
    let q = def
        .r2()
        .narrow()?
        .checked_mul(Rat::new(ee.denom(), ee.numer())?)?;
    let q = HpBounded::of_rat(q, bits);
    let k = q.mul(&q.inv_sqrt(bits)?, bits);
    let seam = core::array::from_fn(|i| across[i].mul(&k, bits));
    Some(CylinderHp { origin, axis, seam })
}

/// The anchor and normal [`push_plane_realized`] describes, each `None` where it does not answer.
fn realize_plane_cache(
    model: &Model,
    points: &[[nacre_exact::Rat; 3]; 3],
    leaf: Handle<nacre_topo::MotionNode>,
    sense: nacre_topo::Orientation,
) -> (Option<[f64; 3]>, Option<[f64; 3]>) {
    let zero = [nacre_exact::Rat::from_int(0); 3];
    if model.chain_point_rat(leaf, zero).is_some()
        || model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP)
    {
        return (None, None);
    }
    let Some(chain) = motion_chain(model, leaf) else {
        return (None, None);
    };
    let realize_point = |p: &WitnessPoint| on_the_cache_rungs(|bits| Some(p.realize(bits)));
    let Some(name) = nacre_exact::plane_name_exact(points[0], points[1], points[2]) else {
        return (None, None); // collinear: the points state no plane
    };
    let Some(&n) = name.narrow() else {
        let Some(tri) = crate::rotated_vertex::replayed_triangle(*points, &chain) else {
            return (None, None);
        };
        return (realize_point(&tri[0]), spanned_normal(&tri, sense));
    };
    let anchor = replay(WitnessPoint::at(points[0]), &chain).and_then(|p| realize_point(&p));
    (
        anchor,
        carried_name_normal(model, &name, n, points, leaf, sense, &chain),
    )
}

/// The plane's unit normal from its (narrow) pre-motion name `n`, carried by `chain` — `None`
/// where undecided at 256 bits.
fn carried_name_normal(
    model: &Model,
    name: &nacre_exact::PlaneName,
    n: [nacre_exact::Rat; 4],
    points: &[[nacre_exact::Rat; 3]; 3],
    leaf: Handle<nacre_topo::MotionNode>,
    sense: nacre_topo::Orientation,
    chain: &[nacre_judge::MoveNode],
) -> Option<[f64; 3]> {
    let zero = [nacre_exact::Rat::from_int(0); 3];
    let along = nacre_exact::name_along_points(
        name,
        points.each_ref().map(|p| MeetPoint::Narrow(*p)).each_ref(),
    )?;
    // The way the plane faces, in the frame its points are written in: the name's normal when it
    // runs with the points' turn times the sense, and a reflection in the chain turns it again.
    let facing = if along { sense.sign() } else { -sense.sign() };
    let turned = facing * crate::rotated_vertex::motion_parity(model, Some(leaf))? < 0;
    let tip = [n[0], n[1], n[2]];
    // `tail → head` is the carried normal; swapping them negates it exactly.
    let (tail, head) = if turned { (tip, zero) } else { (zero, tip) };
    let tail = replay(WitnessPoint::at(tail), chain)?;
    let head = replay(WitnessPoint::at(head), chain)?;
    on_the_cache_rungs(|bits| {
        let (t, h) = (tail.realize(bits), head.realize(bits));
        unit(core::array::from_fn(|k| h[k].sub(&t[k], bits)), bits)
    })
}

/// The unit normal three replayed points span, facing the way `sense` says: the truth's sense is
/// against the points' world turn `(p₁−p₀)×(p₂−p₀)` — a reflection in the chain is already in the
/// replayed points — and swapping the two edges negates it exactly (multiplying the `f64` by `−1`
/// would turn a `+0.0` into `−0.0`).
fn spanned_normal(tri: &[WitnessPoint; 3], sense: nacre_topo::Orientation) -> Option<[f64; 3]> {
    let (p1, p2) = if sense.sign() < 0 {
        (&tri[2], &tri[1])
    } else {
        (&tri[1], &tri[2])
    };
    on_the_cache_rungs(|bits| {
        let [a, b, c, _] = nacre_judge::plane_hp(&tri[0], p1, p2, bits);
        unit([a, b, c], bits)
    })
}

/// `d` scaled to unit length at `bits` — `None` where its squared length is not positive.
fn unit(d: [HpBounded; 3], bits: usize) -> Option<[HpBounded; 3]> {
    let nn = d[0]
        .mul(&d[0], bits)
        .add(&d[1].mul(&d[1], bits), bits)
        .add(&d[2].mul(&d[2], bits), bits);
    let inv = nn.inv_sqrt(bits)?;
    Some(core::array::from_fn(|k| d[k].mul(&inv, bits)))
}

/// Read a realization out on the cache road's two rungs — `LADDER[0]`, then `LADDER[1]` only where
/// the first does not name an `f64` ([`Realized::to_f64`]) — the plane cache's copy of the rule
/// [`realize_cache_tracked`] keeps for vertices.
fn on_the_cache_rungs(realize: impl Fn(usize) -> Option<[HpBounded; 3]>) -> Option<[f64; 3]> {
    [LADDER[0], LADDER[1]].into_iter().find_map(|bits| {
        Realized(Arm::Approached(realize(bits)?, bits))
            .to_f64()
            .map(|(v, _)| v)
    })
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
        rebuild_edges(model, Budget::Paid);
    }
    out
}

/// **How much an edge's realization may spend** — the two budgets a cache is realized under.
///
/// `Cache` is what a push pays, on every push: the cache road's two rungs, and no replay of a chain
/// deeper than [`CACHE_REPLAY_COST_CAP`] (the vertex funnel's own bargain, [`realize_cache`]).
/// `Paid` is what the refine door pays when a caller asks: the whole [`LADDER`], at any depth.
/// Correct rounding makes the answer unique, so where `Cache` answers, `Paid` answers the same bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Budget {
    Cache,
    Paid,
}

impl Budget {
    fn rungs(self) -> &'static [usize] {
        match self {
            Budget::Cache => &LADDER[..2],
            Budget::Paid => &LADDER,
        }
    }
}

/// **What an edge's pusher realizes of its truth** ([`EdgeGiven`]) — the pieces `nacre-topo`
/// cannot, because they sit behind a motion chain only this crate replays.
///
/// The line's direction: where one of its two planes has no world name (a turn off the quarters,
/// a frame), the direction the two planes meet in is the cross product of their normals, read off
/// the planes' coefficients realized at `bits` — the world name's where the plane has one, else the
/// plane its witness triangle spans after its chain ([`PlaneMemo`], the seam table's own road and
/// memo) — scaled to a unit and read out once ([`Realized::to_f64`]), unsigned: the derivation turns
/// it toward the end vertex. Where the names answer, nothing is given — `nacre-topo` reads them
/// itself ([`Model::line_direction_from_names`], asked here so «the names answer» has one spelling).
pub(crate) fn edge_given(
    model: &Model,
    surfaces: [Handle<Surface>; 2],
    budget: Budget,
    memo: &mut PlaneMemo,
) -> EdgeGiven {
    EdgeGiven {
        direction: line_direction(model, surfaces, budget, memo),
    }
}

/// [`edge_given`]'s line arm — `None` where the names answer, a carrier is not a plane, a carrier
/// has no road (a nameless `Through`), the budget does not walk the chain, or no rung decides.
fn line_direction(
    model: &Model,
    s: [Handle<Surface>; 2],
    budget: Budget,
    memo: &mut PlaneMemo,
) -> Option<[f64; 3]> {
    if !s
        .iter()
        .all(|&h| matches!(model.surface(h), Surface::Plane { .. }))
    {
        return None;
    }
    // ★ The guard runs before any chain is walked: the witness triangle replays its plane's whole
    // history, and a push must not pay for a history the cache road would not.
    if budget == Budget::Cache
        && s.iter().any(|&h| {
            crate::planes::world_plane_coeffs(model, h).is_none()
                && model
                    .plane_motion(h)
                    .is_some_and(|leaf| model.motion_deeper_than(leaf, CACHE_REPLAY_COST_CAP))
        })
    {
        return None;
    }
    for &bits in budget.rungs() {
        let a = plane_coeffs_hp(model, s[0], bits, memo).ok()?;
        let b = plane_coeffs_hp(model, s[1], bits, memo).ok()?;
        let x = |i: usize, j: usize| a[i].mul(&b[j], bits).sub(&a[j].mul(&b[i], bits), bits);
        // A cross whose length is not yet apart from zero at these bits climbs.
        let Some(u) = unit([x(1, 2), x(2, 0), x(0, 1)], bits) else {
            continue;
        };
        if let Some((v, _)) = Realized(Arm::Approached(u, bits)).to_f64() {
            return Some(v);
        }
    }
    None
}

/// A plane's coefficients at `bits`: its narrow world name's, exactly, where it has one; else the
/// plane its witness triangle spans after its chain, once per batch ([`PlaneMemo`]). The cap's
/// fork in [`seam_point_met`], spelled once for the edge road.
fn plane_coeffs_hp(
    model: &Model,
    h: Handle<Surface>,
    bits: usize,
    memo: &mut PlaneMemo,
) -> Result<[HpBounded; 4], RealizeError> {
    match crate::planes::world_plane_coeffs(model, h) {
        Some(c) => Ok(c.map(|r| HpBounded::of_rat(r, bits))),
        None => memo.plane(model, h, bits),
    }
}

/// **Push an edge whose cache is realized from its truth** — the edge twin of
/// [`push_vertex_realized`] and [`push_plane_realized`], and the road every edge an operation makes
/// takes: what `nacre-topo` cannot derive is realized here ([`edge_given`], [`Budget::Cache`]) and
/// handed to [`Model::push_edge`].
pub(crate) fn push_edge_realized(
    model: &mut Model,
    surfaces: [Handle<Surface>; 2],
    vertices: [Handle<Vertex>; 2],
    memo: &mut PlaneMemo,
) -> Result<Handle<Edge>, nacre_topo::EdgeDecline> {
    model.push_edge(surfaces, vertices, |m| {
        edge_given(m, surfaces, Budget::Cache, memo)
    })
}

/// Re-derive every live edge's curve with what `budget` realizes — the given pieces computed first
/// against the model as it stands, then handed to [`Model::rebuild_edge_cache`] as a table.
pub(crate) fn rebuild_edges(model: &mut Model, budget: Budget) {
    let mut memo = PlaneMemo::default();
    model.rebuild_edge_cache(|m, e| edge_given(m, m.edge(e).surfaces, budget, &mut memo));
}

/// **The edge cache as a push would derive it now** — every live edge re-derived on the push's own
/// budget. The «discard and regenerate» lock's door: on a model nothing has refined it changes no
/// bit.
#[cfg(any(test, feature = "test-util"))]
pub fn rebuild_edge_cache(model: &mut Model) {
    rebuild_edges(model, Budget::Cache);
}

/// **The edge cache as the refine door derives it** — every live edge re-derived on the door's
/// budget, which walks any chain. After the door this changes no bit; the push's budget
/// ([`rebuild_edge_cache`]) would, wherever the door realized a direction behind a chain deeper
/// than the push pays for.
#[cfg(any(test, feature = "test-util"))]
pub fn rebuild_edge_cache_paid(model: &mut Model) {
    rebuild_edges(model, Budget::Paid);
}

/// **A line edge's direction from its two endpoints' definitions** — the instrument's independent
/// road to the truth the edge road realizes: both endpoints realized at `bits` from their own
/// definitions ([`realize_vertex`]'s road), their difference scaled to a unit at `bits` and read out
/// by the one rule ([`Realized::to_f64`]), run from `vertices[0]` to `vertices[1]`. `None` where an
/// endpoint does not realize or the rung does not decide.
///
/// The second value says whether both endpoints were met in a shared frame (an exact meet, folded
/// or replayed) — the road that shares nothing with the edge's plane coefficients. An endpoint
/// realized by the mixed road ([`build_meet`]) reads the same witness planes the edge road does.
#[cfg(any(test, feature = "test-util"))]
pub fn line_direction_from_endpoints(
    model: &Model,
    e: Handle<Edge>,
    bits: usize,
) -> Option<([f64; 3], bool)> {
    let ends = model.edge(e).vertices;
    let hp = |v: Handle<Vertex>| -> Option<[HpBounded; 3]> {
        match build(model, model.vertex(v), bits, &mut Accel::default())
            .ok()?
            .0
        {
            Arm::Approached(p, _) => Some(p),
            Arm::Exact(n, d) => {
                let d = HpBounded::of_bigint(&d, bits);
                let mut out = Vec::with_capacity(3);
                for k in &n {
                    out.push(HpBounded::of_bigint(k, bits).div(&d, bits)?);
                }
                out.try_into().ok()
            }
        }
    };
    let (a, b) = (hp(ends[0])?, hp(ends[1])?);
    let u = unit(core::array::from_fn(|k| b[k].sub(&a[k], bits)), bits)?;
    let shared = ends.iter().all(|&v| {
        model
            .vertex_meet_of(model.vertex(v))
            .is_some_and(|(meet, _)| meet.narrow().is_some())
    });
    Some((Realized(Arm::Approached(u, bits)).to_f64()?.0, shared))
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
        match build(model, def, bits, &mut Accel::default()) {
            Ok(r) => {
                if r.is_exact() {
                    return Ok(r);
                }
                if let Some(r) = decided(r) {
                    return Ok(r);
                }
            }
            // Not decided at these bits: more bits are exactly the thing (a seam met in the world
            // whose run into its cap is not yet apart from zero after a long chain).
            Err(RealizeError::Undecided) => continue,
            // A structural refusal does not improve with precision.
            Err(e) => return Err(e),
        }
    }
    // Every rung built, none decided — the only thing running out of ladder can mean.
    Err(RealizeError::Undecided)
}

/// One realization at `bits`, from the definition.
/// `out` is the accelerator's channel and only the three-plane road fills it — a curved definition
/// realizes from its own geometry, or (a seam under a motion that does not fold) by replaying its
/// carriers' chains from their statements, never from a shared prefix, so it has none to hand on.
fn build(
    model: &Model,
    def: &Vertex,
    bits: usize,
    acc: &mut Accel<'_>,
) -> Result<Realized, RealizeError> {
    match *def {
        nacre_topo::Vertex::ThreePlane(_) => build_three_plane(model, def, bits, acc),
        nacre_topo::Vertex::OnSeam([cyl, cap]) => {
            seam_point(model, cyl, cap, bits).map(|p| Realized(Arm::Approached(p, bits)))
        }
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
/// point, exactly designated"*; this realizes that designation instead of reading the cache.
///
/// Two roads, one point: where the cylinder and its cap are stated in the world the rim's centre
/// is rational and only `r·ê` is not ([`seam_point_stated`]); where either is carried by a chain
/// that does not fold — a cylinder turned off the quarters, a circle extruded on a frame node —
/// the seam's line is replayed and met with the cap in the world ([`seam_point_met`]).
fn seam_point(
    model: &Model,
    cyl: Handle<Surface>,
    cap: Handle<Surface>,
    bits: usize,
) -> Result<[HpBounded; 3], RealizeError> {
    match seam_point_stated(model, cyl, cap, bits) {
        Some(p) => Ok(p),
        None => seam_point_met(model, cyl, cap, bits),
    }
}

/// [`seam_point`] where the cylinder and the cap are stated in the world.
fn seam_point_stated(
    model: &Model,
    cyl: Handle<Surface>,
    cap: Handle<Surface>,
    bits: usize,
) -> Option<[HpBounded; 3]> {
    let def = model.world_cylinder_def(cyl)?;
    let coeffs = crate::planes::world_plane_coeffs(model, cap)?;
    let (o, m, r2, e) = (def.origin(), def.dir(), def.r2(), def.ref_dir());
    let centre = nacre_exact::axis_plane_meet(&coeffs, &o, &m)?;
    nacre_exact::realize_seam_point(centre, perp_component(&e, &m)?, r2, bits)
}

/// [`seam_point`] met in the world: the seam's line — the cylinder's origin plus `r·ê`, along the
/// axis, both carried by the cylinder's own chain ([`world_cylinder_hp`]) — crossed with the cap's
/// plane in the world, at `bits`.
///
/// ★ **The cap and the cylinder are often on different chains**, so the meet cannot be solved
/// before the motion: a cylinder turned about its own axis takes a node while its caps, fixed by
/// the turn, stay as stated; a circle extruded on a frame node stands on the frame's own plane.
/// So the cap answers for itself — its world name where it has one, else the plane its witness
/// triangle spans after its own chain ([`nacre_judge::plane_hp`]).
///
/// [`RealizeError::NoCurvedPoint`] where a chain or the statement cannot be read;
/// [`RealizeError::Undecided`] where the line's run into the cap is not decided away from zero at
/// `bits` (`HpBounded::div`) — a long chain spends a bit a turn, and the ladder climbs.
fn seam_point_met(
    model: &Model,
    cyl: Handle<Surface>,
    cap: Handle<Surface>,
    bits: usize,
) -> Result<[HpBounded; 3], RealizeError> {
    let Surface::Cylinder { def, motion } = model.surface(cyl) else {
        return Err(RealizeError::NoCurvedPoint);
    };
    let hp = world_cylinder_hp(model, def, *motion, bits).ok_or(RealizeError::NoCurvedPoint)?;
    let on_rim: [HpBounded; 3] = core::array::from_fn(|k| hp.origin[k].add(&hp.seam[k], bits));
    let plane: [HpBounded; 4] = match crate::planes::world_plane_coeffs(model, cap) {
        Some(c) => c.map(|r| HpBounded::of_rat(r, bits)),
        None => {
            let tri = crate::rotated_vertex::surface_witness_triangle(model, cap)
                .ok_or(RealizeError::NoCurvedPoint)?;
            nacre_judge::plane_hp(&tri[0], &tri[1], &tri[2], bits)
        }
    };
    let dot = |v: &[HpBounded; 3]| {
        plane[0]
            .mul(&v[0], bits)
            .add(&plane[1].mul(&v[1], bits), bits)
            .add(&plane[2].mul(&v[2], bits), bits)
    };
    // `on_rim − s·axis` lies on `n·x + c = 0` for `s = (n·on_rim + c) / (n·axis)`.
    let s = dot(&on_rim)
        .add(&plane[3], bits)
        .div(&dot(&hp.axis), bits)
        .ok_or(RealizeError::Undecided)?;
    Ok(core::array::from_fn(|k| {
        on_rim[k].sub(&s.mul(&hp.axis[k], bits), bits)
    }))
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
    acc: &mut Accel<'_>,
) -> Result<Realized, RealizeError> {
    let Some((meet, frame)) = model.vertex_meet_of(def) else {
        return build_meet(model, def, bits, acc.planes.as_deref_mut());
    };
    let Some(node) = frame else {
        // No motion: the meet *is* the coordinate, and `lift` states it as integers whatever its
        // width — so `Wide` is not a refusal on this road.
        let (n, d) = meet.lift();
        return Ok(Realized(Arm::Exact(n, d)));
    };
    // Replaying a motion needs a rational base (`WitnessPoint::at` takes `[Rat; 3]`); a meet
    // wider than that takes the carriers' own witnesses instead.
    let Some(&base) = meet.narrow() else {
        return build_meet(model, def, bits, acc.planes.as_deref_mut());
    };
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
    acc.prefix = Some(PrefixWrite {
        used,
        key: (base, node, bits),
        value: (chain.len(), point.clone()),
    });
    Ok(Realized(Arm::Approached(point, bits)))
}

/// **A meet no shared frame states** — carriers whose motion histories differ (a turned block's
/// wall against an unturned one's), or a meet too wide to replay: each carrier's witness triangle
/// (`rotated_vertex::surface_witness_triangle`) is realized through its own chain into the plane's
/// coefficients ([`nacre_judge::plane_hp`]) and the three planes are met at `bits`
/// ([`nacre_judge::meet_hp`]), so the error is the realization's and the value is read out once.
/// A batch's [`PlaneMemo`] realizes each carrier plane once.
///
/// ★ Not the rejected road of replaying three points in `f64` to state a plane (design 「가지 말
/// 것」): that rounds at every step of the chain, and this rounds once, at the end.
///
/// ⚠ `Undecided` where `D` may be zero at this precision — which is also what three carriers
/// that truly meet in no point answer, rung after rung: a vertex is a point, so that population
/// is empty for a valid definition, and the ladder's ceiling is where it would stop.
fn build_meet(
    model: &Model,
    def: &Vertex,
    bits: usize,
    mut memo: Option<&mut PlaneMemo>,
) -> Result<Realized, RealizeError> {
    let nacre_topo::Vertex::ThreePlane(tri) = *def else {
        return Err(RealizeError::NoMeet);
    };
    let mut plane = |h: Handle<Surface>| match memo.as_deref_mut() {
        Some(m) => m.plane(model, h, bits),
        None => realize_plane(model, h, bits),
    };
    let planes = [plane(tri[0])?, plane(tri[1])?, plane(tri[2])?];
    let p = nacre_judge::meet_hp(&planes, bits).ok_or(RealizeError::Undecided)?;
    Ok(Realized(Arm::Approached(p, bits)))
}

/// One carrier plane's coefficients at `bits`, from its witness triangle. `NoMeet` for a carrier
/// that has none (a nameless `Through` carrier — depth — or a chain outside the decimal window).
fn realize_plane(
    model: &Model,
    h: Handle<Surface>,
    bits: usize,
) -> Result<[HpBounded; 4], RealizeError> {
    let [p0, p1, p2] =
        crate::rotated_vertex::surface_witness_triangle(model, h).ok_or(RealizeError::NoMeet)?;
    Ok(nacre_judge::plane_hp(&p0, &p1, &p2, bits))
}

/// **Carrier planes realized at a precision, for one batch of realizations on one model** — the
/// accelerator of the mixed road ([`build_meet`]) and of the edge road ([`edge_given`]), one memo
/// per operation so a boolean's seam table and its edges realize a plane once between them. A
/// turned solid's wall is a carrier of every corner it makes, and the seam table asks the corners
/// of a whole result at once: realizing the wall's coefficients again for each corner was the whole
/// cost (a fold of 80 turned fins: 1.1 s → 1.8 s; with this memo 1.2 s).
///
/// ★ It cannot change an answer: the value is a function of the plane's truth and the precision,
/// which is the key. Keyed by handle, so it lives no longer than one operation on one model.
#[derive(Default)]
pub(crate) struct PlaneMemo(std::collections::HashMap<(Handle<Surface>, usize), [HpBounded; 4]>);

impl PlaneMemo {
    fn plane(
        &mut self,
        model: &Model,
        h: Handle<Surface>,
        bits: usize,
    ) -> Result<[HpBounded; 4], RealizeError> {
        if let Some(p) = self.0.get(&(h, bits)) {
            return Ok(p.clone());
        }
        let p = realize_plane(model, h, bits)?;
        self.0.insert((h, bits), p.clone());
        Ok(p)
    }
}

/// Whether the mixed road ([`build_meet`]) would replay more chain nodes than
/// [`CACHE_REPLAY_COST_CAP`] pays for: three witness points through each carrier's chain. The road
/// is asked second, and only when the sum passes the cap — it is the same two branches
/// [`build_three_plane`] takes to `build_meet` (no shared frame, or a meet too wide to replay).
fn meet_road_over_cap(model: &Model, def: &Vertex) -> bool {
    let nodes: usize = def
        .carriers()
        .filter_map(|h| model.plane_motion(h))
        .map(|leaf| 3 * model.motion_depth_up_to(leaf, CACHE_REPLAY_COST_CAP))
        .sum();
    nodes > CACHE_REPLAY_COST_CAP
        && match model.vertex_meet_of(def) {
            None => true,
            Some((meet, frame)) => frame.is_some() && meet.narrow().is_none(),
        }
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

/// **The accelerators a realization may use** — what it may read to go faster ([`PlaneMemo`], for
/// a batch) and what it learned that the prefix table could keep ([`PrefixWrite`]). Neither
/// changes an answer; both are why this travels beside the definition instead of inside it.
#[derive(Default)]
pub(crate) struct Accel<'a> {
    pub(crate) prefix: Option<PrefixWrite>,
    pub(crate) planes: Option<&'a mut PlaneMemo>,
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

    /// **A turned cylinder's cache is its chain realized and rounded once — the same `f64` at twice
    /// the precision.** A cylinder turned 37° about `x` and one turned 30° about `y`, both about
    /// pivots off the origin (so a direction replayed as a point would carry the pivot's
    /// translation and miss), have no world statement; the cache their push left is read against
    /// [`world_cylinder_hp`] at 512 bits, rounded the same way.
    #[test]
    fn a_turned_cylinder_cache_is_its_chain_rounded_once() {
        use nacre_exact::{Angle, Axis, Isometry, Rotation};
        for (axis, deg, pivot) in [
            (Axis::X, 37, [r(1), r(0), r(0)]),
            (Axis::Y, 30, [r(0), r(1), r(2)]),
        ] {
            let mut m = Model::new();
            let c = crate::fixtures::cylinder(
                &mut m,
                Point3::from_array([1.0, 2.0, 0.5]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                1.5,
                2.0,
            );
            let Ok(crate::OpOutput::Transform { solid }) = crate::apply(
                &mut m,
                &crate::Operation::Transform {
                    solid: c.solid,
                    isometry: Isometry::rotation(Rotation {
                        axis,
                        pivot,
                        angle: Angle::from_deg(r(deg)).unwrap(),
                    }),
                },
            ) else {
                panic!("the turn")
            };
            let lateral = m
                .shell(m.solid(solid).outer)
                .faces
                .iter()
                .map(|&f| m.face(f).surface)
                .find(|&s| matches!(m.surface(s), Surface::Cylinder { .. }))
                .expect("a lateral face");
            let Surface::Cylinder { def, motion } = m.surface(lateral).clone() else {
                unreachable!()
            };
            assert!(
                motion.is_some() && m.world_cylinder_def(lateral).is_none(),
                "{axis:?}"
            );
            let hp = world_cylinder_hp(&m, &def, motion, 512).expect("the chain");
            let read = |v: [HpBounded; 3]| Realized(Arm::Approached(v, 512)).to_f64().unwrap().0;
            let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(lateral) else {
                unreachable!()
            };
            assert_eq!(
                cache.axis().origin().as_array(),
                read(hp.origin),
                "{axis:?}"
            );
            assert_eq!(
                cache.axis().direction().as_array(),
                read(unit(hp.axis, 512).unwrap()),
                "{axis:?}"
            );
            assert_eq!(
                cache.ref_dir().as_array(),
                read(unit(hp.seam, 512).unwrap()),
                "{axis:?}"
            );
            assert_eq!(cache.radius(), 1.5);
            // And independently: the original cylinder's cache turned in `f64`.
            let pv = pivot.map(|r| r.to_f64());
            let deg = f64::from(deg as i32);
            let at = |p: [f64; 3]| turned_f64(p, axis, deg, pv);
            let dir = |v: [f64; 3]| {
                let (a, o) = (at(v), at([0.0; 3]));
                [a[0] - o[0], a[1] - o[1], a[2] - o[2]]
            };
            assert!(
                near(cache.axis().origin().as_array(), at([1.0, 2.0, 0.5])),
                "{axis:?}"
            );
            assert!(
                near(cache.axis().direction().as_array(), dir([0.0, 0.0, 1.0])),
                "{axis:?}"
            );
            assert!(
                near(cache.ref_dir().as_array(), dir([1.0, 0.0, 0.0])),
                "{axis:?}"
            );
        }
    }

    /// `p` turned `deg` degrees about `axis` through `pivot`, in `f64` — the independent reading the
    /// turned-cylinder locks hold the realizations to (a derivation that shares nothing with the
    /// chain replay it checks).
    fn turned_f64(p: [f64; 3], axis: nacre_exact::Axis, deg: f64, pivot: [f64; 3]) -> [f64; 3] {
        let (s, c) = deg.to_radians().sin_cos();
        let v = [p[0] - pivot[0], p[1] - pivot[1], p[2] - pivot[2]];
        let w = match axis {
            nacre_exact::Axis::X => [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]],
            nacre_exact::Axis::Y => [c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]],
            nacre_exact::Axis::Z => [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]],
        };
        [w[0] + pivot[0], w[1] + pivot[1], w[2] + pivot[2]]
    }

    fn near(a: [f64; 3], b: [f64; 3]) -> bool {
        (0..3).all(|k| (a[k] - b[k]).abs() < 1e-12)
    }

    /// The seam vertices of `solid`'s faces, each with the `f64` of [`seam_point_met`] at 512 bits.
    fn seams_met(m: &Model, solid: Handle<nacre_topo::Solid>) -> Vec<(Handle<Vertex>, [f64; 3])> {
        let mut out = Vec::new();
        for &fh in &m.shell(m.solid(solid).outer).faces {
            let f = m.face(fh);
            for he in std::iter::once(&f.outer).flat_map(|lp| lp.half_edges.iter()) {
                for vh in m.edge(he.edge).vertices {
                    if let Vertex::OnSeam([cyl, cap]) = *m.vertex(vh)
                        && !out.iter().any(|&(v, _)| v == vh)
                    {
                        let p = seam_point_met(m, cyl, cap, 512).expect("the seam meets its cap");
                        out.push((vh, Realized(Arm::Approached(p, 512)).to_f64().unwrap().0));
                    }
                }
            }
        }
        out
    }

    /// **Two roads, one point.** Where the cylinder and its caps are stated in the world, the seam
    /// realized on the rational road (`centre + r·ê`) and met in the world agree to the bit — an
    /// axis-aligned cylinder off the origin and a 3-4-5 tilted one.
    #[test]
    fn the_seam_met_in_the_world_is_the_stated_seam() {
        for axis in [[0.0, 0.0, 1.0], [0.0, 0.6, 0.8]] {
            let mut m = Model::new();
            let c = crate::fixtures::cylinder(
                &mut m,
                Point3::from_array([1.0, 2.0, 0.5]),
                nacre_math::Vector3::from_array(axis),
                1.5,
                2.0,
            );
            let seams = seams_met(&m, c.solid);
            assert!(!seams.is_empty(), "{axis:?}");
            for (vh, met) in seams {
                let Vertex::OnSeam([cyl, cap]) = *m.vertex(vh) else {
                    unreachable!()
                };
                assert!(
                    m.world_cylinder_def(cyl).is_some(),
                    "{axis:?}: a stated cylinder"
                );
                let stated = seam_point_stated(&m, cyl, cap, 512).expect("the stated road");
                assert_eq!(
                    Realized(Arm::Approached(stated, 512)).to_f64().unwrap().0,
                    met,
                    "{axis:?}"
                );
            }
        }
    }

    /// **A turned cylinder's seam is met in the world.** Three placements where the cap and the
    /// cylinder ride different chains or none of them folds — a cylinder turned about its own axis
    /// (its caps, fixed by the turn, stay as stated), one on an axis whose frame is a node (its
    /// base cap is the frame's plane), and one turned twice: every seam vertex is `Bounded`, its
    /// cache the met point at 512 bits rounded the same way.
    #[test]
    fn a_turned_cylinder_seam_is_realized_where_it_meets_its_cap() {
        use nacre_exact::{Angle, Axis, Isometry, Rotation};
        let turn = |m: &mut Model, s, axis, deg: i128, pivot: [Rat; 3]| {
            let Ok(crate::OpOutput::Transform { solid }) = crate::apply(
                m,
                &crate::Operation::Transform {
                    solid: s,
                    isometry: Isometry::rotation(Rotation {
                        axis,
                        pivot,
                        angle: Angle::from_deg(r(deg)).unwrap(),
                    }),
                },
            ) else {
                panic!("the turn")
            };
            solid
        };
        let up = nacre_math::Vector3::from_array([0.0, 0.0, 1.0]);
        let base = Point3::from_array([1.0, 2.0, 0.5]);
        // The independent reading: the unturned seam points `(2.5, 2, 0.5)` and `(2.5, 2, 2.5)`
        // turned in `f64`; for the frame node, the cylinder's own definition.
        let seams_before = [[2.5, 2.0, 0.5], [2.5, 2.0, 2.5]];
        type Check = Box<dyn Fn([f64; 3]) -> bool>;
        let turned_from = |steps: Vec<(Axis, f64, [f64; 3])>| -> Check {
            Box::new(move |p| {
                seams_before.iter().any(|&q| {
                    let q = steps
                        .iter()
                        .fold(q, |q, &(a, d, pv)| turned_f64(q, a, d, pv));
                    near(p, q)
                })
            })
        };
        let on_the_tilted_rim: Check = Box::new(|p| {
            let k = 1.0 / 3f64.sqrt();
            let v = [p[0] - 1.0, p[1] - 2.0, p[2] - 0.5];
            let h = (v[0] + v[1] + v[2]) * k;
            let radial = (v.iter().map(|c| c * c).sum::<f64>() - h * h).sqrt();
            (radial - 1.5).abs() < 1e-12 && (h.abs() < 1e-12 || (h - 2.0).abs() < 1e-12)
        });
        let cases: Vec<(&str, Model, Handle<nacre_topo::Solid>, Check)> = vec![
            {
                let mut m = Model::new();
                let c = crate::fixtures::cylinder(&mut m, base, up, 1.5, 2.0);
                let s = turn(&mut m, c.solid, Axis::Z, 37, [r(1), r(0), r(0)]);
                let check = turned_from(vec![(Axis::Z, 37.0, [1.0, 0.0, 0.0])]);
                ("about its own axis", m, s, check)
            },
            {
                let mut m = Model::new();
                let tilted = nacre_math::Vector3::from_array([1.0, 1.0, 1.0]);
                let c = crate::fixtures::cylinder(&mut m, base, tilted, 1.5, 2.0);
                ("on a frame node", m, c.solid, on_the_tilted_rim)
            },
            {
                let mut m = Model::new();
                let c = crate::fixtures::cylinder(&mut m, base, up, 1.5, 2.0);
                let s = turn(&mut m, c.solid, Axis::X, 37, [r(0), r(1), r(0)]);
                let s = turn(&mut m, s, Axis::Y, 30, [r(2), r(0), r(1)]);
                let check = turned_from(vec![
                    (Axis::X, 37.0, [0.0, 1.0, 0.0]),
                    (Axis::Y, 30.0, [2.0, 0.0, 1.0]),
                ]);
                ("turned twice", m, s, check)
            },
        ];
        for (what, m, solid, check) in cases {
            let seams = seams_met(&m, solid);
            assert!(!seams.is_empty(), "{what}");
            for (vh, met) in seams {
                let Vertex::OnSeam([cyl, _]) = *m.vertex(vh) else {
                    unreachable!()
                };
                assert!(
                    m.world_cylinder_def(cyl).is_none(),
                    "{what}: no world statement"
                );
                let PointCache::Bounded { coord, .. } = m.vertex_cache(vh) else {
                    panic!("{what}: {:?}", m.vertex_cache(vh))
                };
                assert_eq!(coord.as_array(), met, "{what}");
                assert!(
                    check(met),
                    "{what}: {met:?} is not where the definition puts it"
                );
            }
        }
    }

    /// **Past the replay budget the seam waits for the paid door, and the door answers.** A cylinder
    /// turned 193 times — its chain one node past [`CACHE_REPLAY_COST_CAP`] — keeps its seam
    /// vertices `Ceiling` (the cache road declines by cost, not for want of a road), and
    /// [`refine_vertex_cache`] raises every one of them to `Bounded` through the same meet.
    #[test]
    fn a_seam_past_the_replay_budget_waits_for_the_paid_door() {
        use nacre_exact::{Angle, Axis, Isometry, Rotation};
        let mut m = Model::new();
        let mut solid = crate::fixtures::cylinder(
            &mut m,
            Point3::from_array([1.0, 2.0, 0.5]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            1.5,
            2.0,
        )
        .solid;
        for _ in 0..=CACHE_REPLAY_COST_CAP {
            let Ok(crate::OpOutput::Transform { solid: s }) = crate::apply(
                &mut m,
                &crate::Operation::Transform {
                    solid,
                    isometry: Isometry::rotation(Rotation {
                        axis: Axis::X,
                        pivot: [r(0); 3],
                        angle: Angle::from_deg(r(37)).unwrap(),
                    }),
                },
            ) else {
                panic!("the turn")
            };
            solid = s;
        }
        m.rebuild_adjacency();
        let seams: Vec<Handle<Vertex>> = m
            .reachable()
            .vertices
            .into_iter()
            .filter(|&v| matches!(*m.vertex(v), Vertex::OnSeam(_)))
            .filter(|&v| {
                m.shell(m.solid(solid).outer).faces.iter().any(|&f| {
                    m.face(f)
                        .outer
                        .half_edges
                        .iter()
                        .any(|he| m.edge(he.edge).vertices.contains(&v))
                })
            })
            .collect();
        assert!(!seams.is_empty());
        // On that chain the cylinder funnel leaves a new statement the figure it was handed, bit
        // for bit. ★ This pins the outcome, not the cost guard: past 192 turns the 256-bit rung
        // cannot decide either (the cap is sized so it can inside it), so the funnel without its
        // depth test answers the same — measured, that plant stays green. The guard's effect is
        // the replay not paid, which no result shows.
        let Some(Surface::Cylinder { def, motion }) = m
            .shell(m.solid(solid).outer)
            .faces
            .iter()
            .map(|&f| m.surface(m.face(f).surface).clone())
            .find(|s| matches!(s, Surface::Cylinder { .. }))
        else {
            panic!("a lateral face")
        };
        let wider = nacre_topo::CylinderDef::new(
            def.origin(),
            def.dir(),
            def.ref_dir(),
            nacre_exact::BigRat::from(r(4)),
        )
        .unwrap();
        let figure = nacre_geom::Cylinder::from_axis(
            Point3::from_array([0.1, 0.2, 0.3]),
            nacre_math::Vector3::from_array([0.0, 0.6, 0.8]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            2.0,
        )
        .unwrap();
        let h = push_cylinder_realized(&mut m, figure, wider, motion);
        assert_eq!(m.surface_cache(h), &nacre_geom::Surface::Cylinder(figure));
        for &v in &seams {
            assert!(
                matches!(m.vertex_cache(v), PointCache::Ceiling { .. }),
                "{:?}",
                m.vertex_cache(v)
            );
        }
        let report = refine_vertex_cache(&mut m);
        assert!(report.refined >= seams.len(), "{report:?}");
        for &v in &seams {
            assert!(
                matches!(m.vertex_cache(v), PointCache::Bounded { .. }),
                "{:?}",
                m.vertex_cache(v)
            );
        }
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
