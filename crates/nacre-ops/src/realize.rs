//! **Asking a vertex for its coordinate at a precision you choose.**
//!
//! Everything the kernel shows today — the viewport, the tessellation, a report row — reads the
//! f64 point cache, which carries its own error (`PointCache.tol`) and is *not* promised to be the
//! nearest f64 to the truth. Printing more digits of it would print the rounding, not the point.
//!
//! This module goes the other way: it takes the vertex's **definition** and realizes a coordinate
//! from it, rounding exactly once at the end. Two roads meet here and they are chosen by the
//! *value*, not by the `VertexDef` variant:
//!
//! - **rational** — a three-plane meet is an exact ratio ([`nacre_topo::Model::vertex_meet`]), so
//!   its decimals come out of one long division and every digit printed is a digit of the
//!   coordinate itself. There is no rounding question to get wrong.
//! - **realized** — anything a motion or a radical reaches is approached at `prec` bits with the
//!   error the realization cost, and a digit is printed only once the interval decides it.
//!
//! ★ **Where this sits.** `nacre-topo` does not depend on `nacre-cip`, so the composition
//! `vertex_meet → WitnessPoint::at → replay → realize` cannot live on `Model`; `nacre-ops` is the
//! first crate that sees both. The rounding itself lives one layer further down, in
//! `nacre-scalar`, beside `round_to_f64` — cip realizes, scalar rounds, this module composes.
//!
//! ★★ **Migration 3b.** [`nacre_topo::VertexDef::OnSeam`]'s doc records the missing piece as *"the
//! machinery that regenerates the cached coordinate"*, deferred with row 3b. This is the asking
//! half of it. It does **not** overwrite the cache — that is the row's other half.

use crate::rotated_vertex::{motion_chain, replay};
use nacre_cip::WitnessPoint;
use nacre_scalar::{Bounded, Mag, MeetPoint};
use nacre_store::Handle;
use nacre_topo::{Model, Surface, Vertex};
use num_bigint::BigInt;

/// How precisely to realize — always stated, never defaulted.
///
/// ★ There is no `Digits` here on purpose: the digit count belongs to
/// [`Realized::to_decimal`] alone. Carrying it in both places lets the two disagree, and the way
/// they disagree is by printing digits the realization never determined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    /// Realize from the definition and round once to the nearest `f64`, escalating until the
    /// interval names one. This is the value the cache *should* hold.
    NearestF64,
    /// An explicit working precision — for measurement, and for a caller driving its own ladder.
    Bits(usize),
}

/// A coordinate that has been realized, together with what it cost.
///
/// Opaque: the arms are an implementation detail, and a caller that could see them would be
/// coupled to astro-float's types through this crate as well as through `nacre-scalar`.
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
    Approached([Bounded; 3], usize),
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
}

/// The ladder. Doubling, so a value needing `n` bits pays at most `2n`; capped, because an
/// unbounded climb on an undecidable ask is a hang rather than an answer.
const LADDER: [usize; 6] = [128, 256, 512, 1024, 2048, 4096];

impl Realized {
    /// The nearest `f64` per coordinate, with the error each carries. `None` where the
    /// realization does not name one.
    ///
    /// ★★★ **The error comes back as [`Mag`], not `f64`, and that is what that type is for.**
    /// `nacre-scalar`'s own test says so: *"the reason this type exists: a radius the ladder
    /// actually produces must not become zero. An `f64` cannot hold `2⁻²⁰⁴⁸`."* A realization at
    /// 4096 bits carries exactly such a radius, so handing the bound over as an `f64` forces the
    /// conversion the type was built to prevent — an earlier spelling did, reported `0e0`, and
    /// became indistinguishable from the exact arm's honest zero. Clamping papered over that;
    /// staying in `Mag` removes it.
    ///
    /// The *value* is an `f64` because that is what the caller asked for. The *bound* on it need
    /// not be one, and at the top of the ladder cannot be.
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
                    let (val, no_loss) = nacre_scalar::nearest_f64_big_exact(&n[k], d)?;
                    v[k] = val;
                    // Half an ulp bounds a correct rounding — as a `Mag`, so a tiny coordinate's
                    // bound cannot vanish on the way out either.
                    e[k] = match no_loss {
                        true => Mag::ZERO,
                        false => Mag::of(val).times(Mag::pow2(-53)),
                    };
                }
                Some((v, e))
            }
            Arm::Approached(p, prec) => {
                let mut v = [0.0; 3];
                let mut e = [Mag::ZERO; 3];
                for k in 0..3 {
                    v[k] = nacre_scalar::round_to_f64(&p[k].0, p[k].1, *prec)?;
                    // The realization's own radius, handed over unchanged — no conversion, so
                    // nothing to underflow. That is what `Mag` is for: its own test records
                    // that an `f64` cannot hold `2⁻²⁰⁴⁸` and a deep rung's radius is exactly
                    // that small. An earlier spelling converted here and reported `0e0`.
                    e[k] = p[k].1;
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
                nacre_scalar::decimals_of_ratio(&n[k], d, places)
            })),
            Arm::Approached(p, _) => {
                let mut out = [const { String::new() }; 3];
                for k in 0..3 {
                    out[k] = nacre_scalar::round_to_digits(&p[k].0, p[k].1, places)?;
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

/// Realize one vertex's coordinate from its definition.
pub fn realize_vertex(
    model: &Model,
    v: Handle<Vertex>,
    p: Precision,
) -> Result<Realized, RealizeError> {
    match p {
        Precision::Bits(bits) => build(model, v, bits),
        Precision::NearestF64 => climb(model, v, |r| r.to_f64().map(|_| r)),
    }
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
    let out = climb(model, v, |r| r.to_decimal(places).map(|_| r))?;
    out.to_decimal(places).ok_or(RealizeError::Undecided)
}

/// Walk [`LADDER`] until `decided` accepts a realization. The exact arm decides on the first rung
/// whatever is asked, so a rational vertex never climbs.
fn climb(
    model: &Model,
    v: Handle<Vertex>,
    decided: impl Fn(Realized) -> Option<Realized>,
) -> Result<Realized, RealizeError> {
    for bits in LADDER {
        match build(model, v, bits) {
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
fn build(model: &Model, v: Handle<Vertex>, bits: usize) -> Result<Realized, RealizeError> {
    match model.vertices.get(v).def {
        nacre_topo::VertexDef::ThreePlane(_) => build_three_plane(model, v, bits),
        nacre_topo::VertexDef::OnSeam([cyl, cap]) => {
            curved(seam_point(model, cyl, cap, bits), bits)
        }
        nacre_topo::VertexDef::Pierce {
            planes,
            cylinder,
            root,
        } => curved(pierce_point(model, planes, cylinder, root, bits), bits),
    }
}

fn curved(p: Option<[Bounded; 3]>, bits: usize) -> Result<Realized, RealizeError> {
    Ok(Realized(Arm::Approached(
        p.ok_or(RealizeError::NoCurvedPoint)?,
        bits,
    )))
}

/// **The seam vertex is the rim's `+ref_dir` point** — `centre + r·ê`, `ê` the reference
/// direction's unit part perpendicular to the axis. `VertexDef::OnSeam`'s doc calls it *"a unique
/// point, exactly designated"*; this realizes that designation instead of reading the cache, which
/// is the half of migration row 3b that was recorded as missing.
fn seam_point(
    model: &Model,
    cyl: Handle<Surface>,
    cap: Handle<Surface>,
    bits: usize,
) -> Option<[Bounded; 3]> {
    let def = crate::planes::world_cylinder_def(model, cyl)?;
    let coeffs = crate::planes::world_plane_coeffs(model, cap)?;
    let (o, m, r, e) = (def.origin(), def.dir(), def.radius(), def.ref_dir());
    let t = crate::planes::axis_param_of_plane(&coeffs, &def)?;
    let mut centre = o;
    for k in 0..3 {
        centre[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
    }
    nacre_scalar::realize_seam_point(centre, perp_component(&e, &m)?, r, bits)
}

/// **A pierce vertex is the meet line's point at its root** — the two cutting planes give the
/// line, the cylinder gives the quadratic, and `root` names which crossing.
fn pierce_point(
    model: &Model,
    planes: [Handle<Surface>; 2],
    cylinder: Handle<Surface>,
    root: nacre_topo::QuadRoot,
    bits: usize,
) -> Option<[Bounded; 3]> {
    let def = crate::planes::world_cylinder_def(model, cylinder)?;
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let c0 = crate::planes::world_plane_coeffs(model, planes[0])?;
    let c1 = crate::planes::world_plane_coeffs(model, planes[1])?;
    let (line, sv) = pick_root(
        &nacre_scalar::quad::plane_plane_cylinder(&c0, &c1, &o, &m, r)?,
        root,
    )?;
    // `point = base + s·dir`, with `s` the quadratic root realized at `bits` — the only step that
    // is not exact rational arithmetic.
    let sb = nacre_scalar::realize_quad(&sv, bits)?;
    let (b, d) = (line.base(), line.dir());
    let coord = |k: usize| nacre_scalar::affine_bounded(b[k], d[k], &sb, bits);
    Some([coord(0)?, coord(1)?, coord(2)?])
}

fn build_three_plane(
    model: &Model,
    v: Handle<Vertex>,
    bits: usize,
) -> Result<Realized, RealizeError> {
    let (meet, frame) = model.vertex_meet(v).ok_or(RealizeError::NoMeet)?;
    let Some(node) = frame else {
        // No motion: the meet *is* the coordinate, and `lift` states it as integers whatever its
        // width — so `Wide` is not a refusal on this road.
        let (n, d) = meet.lift();
        return Ok(Realized(Arm::Exact(n, d)));
    };
    let base = *narrow_or(&meet)?;
    let chain = motion_chain(model, node).ok_or(RealizeError::NoMotionChain)?;
    let wp = replay(WitnessPoint::at(base), &chain).ok_or(RealizeError::NoMotionChain)?;
    Ok(Realized(Arm::Approached(wp.realize(bits), bits)))
}

fn narrow_or(meet: &MeetPoint) -> Result<&[nacre_scalar::Rat; 3], RealizeError> {
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
    e: &[nacre_scalar::Rat; 3],
    m: &[nacre_scalar::Rat; 3],
) -> Option<[nacre_scalar::Rat; 3]> {
    let (mm, em) = (crate::planes::dot3(m, m)?, crate::planes::dot3(e, m)?);
    let mut out = [nacre_scalar::Rat::from_int(0); 3];
    for k in 0..3 {
        out[k] = mm.checked_mul(e[k])?.checked_sub(em.checked_mul(m[k])?)?;
    }
    Some(out)
}

/// **Which of the crossings `root` names.** `Lo`/`Hi` are ascending parameter along the meet
/// line's direction — `nacre_scalar::quad`'s pair order, which is `VertexDef::Pierce`'s stated
/// convention — and a tangency is one point spelled `Double`, never `Lo`.
///
/// ⚠ Split out because the integration oracle ("the point is on the cylinder") is satisfied by
/// **both** roots: swapping them left the whole corpus green when planted. This is the piece a
/// test can ask the distinguishing question of.
fn pick_root(
    meet: &nacre_scalar::quad::CylinderMeet,
    root: nacre_topo::QuadRoot,
) -> Option<(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal)> {
    use nacre_scalar::quad::{CylinderMeet, QuadVal};
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
    use nacre_scalar::Rat;

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
    /// question that actually distinguishes them, in the vocabulary `VertexDef::Pierce`'s doc
    /// defines: ascending parameter along `n₀ × n₁`.
    #[test]
    fn lo_and_hi_run_along_the_line() {
        use nacre_scalar::quad::CylinderMeet;
        use nacre_topo::QuadRoot;
        let (o, m, rad) = ([r(2), r(2), r(0)], [r(0), r(0), r(1)], r(3));
        // x = 0 and z = 0: the meet line runs along +y and crosses the cylinder twice.
        let c0 = [r(1), r(0), r(0), r(0)];
        let c1 = [r(0), r(0), r(1), r(0)];
        let meet =
            nacre_scalar::quad::plane_plane_cylinder(&c0, &c1, &o, &m, rad).expect("a crossing");
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
        let at = |q: &nacre_scalar::quad::QuadVal| {
            let sb = nacre_scalar::realize_quad(q, 256).expect("realized");
            let (b, d) = (line.base(), line.dir());
            let c = |k: usize| {
                let (v, e) = nacre_scalar::affine_bounded(b[k], d[k], &sb, 256).expect("affine");
                nacre_scalar::round_to_digits(&v, e, 20)
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
