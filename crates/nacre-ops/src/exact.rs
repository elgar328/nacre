//! **Construction arithmetic in rationals**, so that dimensions which ought to
//! coincide actually do.
//!
//! A prism's corners come out of three multiplications and two additions —
//! `origin + x·u + y·v`, then `+ normal·dist` — and in f64 those do not associate.
//! Extruding `7.7` and extruding `1.1` then `6.6` from the same sketch put their top
//! caps one ulp apart, which is not enough to split the model and not little enough
//! to merge: the boolean fuses them into a solid carrying a face of area `9e-16`, and
//! every operation afterwards drags it along. Nothing rejects, nothing is wrong, and
//! nobody is told.
//!
//! Doing the same arithmetic in [`Rat`] fixes it — but only if the rationals are the
//! ones the dimensions were *written* as. Lifting the f64 with `try_from_f64` carries
//! the binary drift in exactly, and exact arithmetic then reproduces it perfectly;
//! [`Rat::from_decimal`] is what recovers `11/10 + 66/10 = 77/10`.
//!
//! **The frame has to be rational too, and often is not.** A world or axis-aligned
//! plane has axes in `{0, ±1}`, exactly representable. A plane at 45° does not: its
//! axes are irrational, and lifting their decimals would give vectors that are no
//! longer unit, so `normal·dist` would come out the wrong *length* — a worse error
//! than the one being fixed. So the entry point checks orthonormality exactly and
//! declines rather than approximating, and the caller keeps today's f64 path. Every
//! function here returns `Option` for that reason, and for i128 overflow, which is
//! the same answer for the same reason.

use crate::ops::{Profile2d, SketchPlane};
use nacre_geom::mixed::Edge2d;
use nacre_math::{Point3, Vector3};
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, MotionNode};

/// A sketch frame whose origin and axes are exact rationals, with the axes proved
/// orthonormal. Only [`SketchPlane::exact`] constructs one, so the proof cannot be
/// bypassed by building the struct directly.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RatFrame {
    origin: [Rat; 3],
    x: [Rat; 3],
    y: [Rat; 3],
}

fn dot(a: &[Rat; 3], b: &[Rat; 3]) -> Option<Rat> {
    let mut acc = Rat::from_int(0);
    for i in 0..3 {
        acc = acc.checked_add(a[i].checked_mul(b[i])?)?;
    }
    Some(acc)
}

fn cross(a: &[Rat; 3], b: &[Rat; 3]) -> Option<[Rat; 3]> {
    let term = |i: usize, j: usize| a[i].checked_mul(b[j])?.checked_sub(a[j].checked_mul(b[i])?);
    Some([term(1, 2)?, term(2, 0)?, term(0, 1)?])
}

fn scale(v: &[Rat; 3], k: Rat) -> Option<[Rat; 3]> {
    Some([
        v[0].checked_mul(k)?,
        v[1].checked_mul(k)?,
        v[2].checked_mul(k)?,
    ])
}

fn sub(a: &[Rat; 3], b: &[Rat; 3]) -> Option<[Rat; 3]> {
    Some([
        a[0].checked_sub(b[0])?,
        a[1].checked_sub(b[1])?,
        a[2].checked_sub(b[2])?,
    ])
}

fn add(a: &[Rat; 3], b: &[Rat; 3]) -> Option<[Rat; 3]> {
    Some([
        a[0].checked_add(b[0])?,
        a[1].checked_add(b[1])?,
        a[2].checked_add(b[2])?,
    ])
}

fn lift(p: [f64; 3]) -> Option<[Rat; 3]> {
    Some([
        Rat::from_decimal(p[0])?,
        Rat::from_decimal(p[1])?,
        Rat::from_decimal(p[2])?,
    ])
}

impl SketchPlane {
    /// This frame in exact rationals, or `None` if it does not have one.
    ///
    /// The test is orthonormality *of the lifted rationals*, checked exactly rather
    /// than within a tolerance. That is the whole condition: given `x·x = y·y = 1` and
    /// `x·y = 0` exactly, the cross product is a unit normal exactly, so `normal·dist`
    /// has exactly the length `dist` and the sweep lands where the dimension says.
    pub(crate) fn exact(&self) -> Option<RatFrame> {
        let f = RatFrame {
            origin: lift(self.origin().as_array())?,
            x: lift(self.x_axis().as_array())?,
            y: lift(self.y_axis().as_array())?,
        };
        let one = Rat::from_int(1);
        let zero = Rat::from_int(0);
        (dot(&f.x, &f.x)? == one && dot(&f.y, &f.y)? == one && dot(&f.x, &f.y)? == zero)
            .then_some(f)
    }
}

impl RatFrame {
    /// **A frame's basis asked in rationals, never realized** — the exact question behind
    /// [`SketchPlane::exact`], put to a frame instead of to a caller's `f64` axes.
    ///
    /// ★★★ **Realizing the axes first and lifting them back does not answer it.** Measured
    /// (`ops::frame_road`): `reduce_direction` turns a `(0.6, 0.8, 0)` axis into the primitive
    /// `(3, 4, 0)` with `uu = 25`, and a realization multiplies by a *numerically* computed `1/5`,
    /// landing on `0.6000000000000001`. `Rat::from_decimal` lifts that, orthonormality fails, and
    /// a plane that is perfectly rational takes the frame-node road — a different arena, not an
    /// ulp. Here `inv_sqrt_exact` answers in `Rat`: `Some` exactly when the axis lands on a
    /// rational, which is the whole question.
    ///
    /// It is the same rule [`nacre_scalar::plane_frame_named`] states for `v̂`, applied one level
    /// up: `v_raw` is carried exactly for this reason, so use it rather than crossing `ŵ × û`.
    ///
    /// `None` when either axis needs an irrational scale, when `v_raw` overflowed at construction
    /// (`PlaneFrame::v` is `None`), or on `i128` overflow — all of which mean the same thing here:
    /// this frame has no exact rational basis, so the sketch is written in the frame instead.
    pub(crate) fn of_plane_frame(pf: &nacre_scalar::PlaneFrame) -> Option<RatFrame> {
        let (v_raw, vv) = pf.v.as_ref()?;
        Some(RatFrame {
            origin: pf.origin,
            x: scale(&pf.u_raw, nacre_scalar::inv_sqrt_exact(pf.uu)?)?,
            y: scale(v_raw, nacre_scalar::inv_sqrt_exact(*vv)?)?,
        })
    }

    /// ★★★★ **The sketch frame of a plane, read in that plane's own frame — where it is the
    /// identity.**
    ///
    /// This is what makes a sketch on a *tilted* face exact. In world coordinates that face's
    /// axes are irrational and [`SketchPlane::exact`] declines; inside the frame the very same
    /// axes are `x̂` and `ŷ`, the origin is the origin, and orthonormality is not something to
    /// check but something to read off. The prism is then built by the arithmetic that was
    /// already here, on rationals that were already exact — the user's own profile decimals.
    ///
    /// What carries it back out to the world is the motion, not this frame: see
    /// [`nacre_topo::Motion::Frame`].
    pub(crate) fn identity() -> Self {
        let (one, zero) = (Rat::from_int(1), Rat::from_int(0));
        RatFrame {
            origin: [zero; 3],
            x: [one, zero, zero],
            y: [zero, one, zero],
        }
    }

    /// The unit normal `x × y`, exactly.
    pub(crate) fn normal(&self) -> Option<[Rat; 3]> {
        cross(&self.x, &self.y)
    }

    /// The three points of the parallel plane `w = d`, stated in **this frame's own coordinate
    /// system** — the exact form of "d away, along the normal".
    ///
    /// ★ This is what a datum offset uses when the frame lifts to exact rationals: the plane can
    /// then be said in the world, so it interns with every other statement of it (`push_plane`
    /// keys on `(name, motion)`, and a frame node would put it under a different key). `None` on
    /// `i128` overflow — which the caller must turn into a named reject rather than quietly taking
    /// the frame-node road, because that road is where the duplicate would appear.
    pub(crate) fn offset_plane_points(&self, d: Rat) -> Option<[[Rat; 3]; 3]> {
        let w = self.normal()?;
        let p0 = add(&self.origin, &scale(&w, d)?)?;
        Some([p0, add(&p0, &self.x)?, add(&p0, &self.y)?])
    }

    /// The sweep vector for a signed distance along the normal — the second of the two
    /// arithmetic steps a prism is built from.
    pub(crate) fn sweep(&self, dist: f64) -> Option<[Rat; 3]> {
        scale(&self.normal()?, Rat::from_decimal(dist)?)
    }

    /// A ring's sketch coordinates placed in space: `origin + x·u + y·v`, the first of
    /// the two steps.
    ///
    /// The coordinates arrive as the profile's stored rational truth — the decimal lift
    /// happened once, in `Profile2d`'s constructor. `None` here is i128 overflow in the
    /// placement arithmetic only.
    pub(crate) fn ring(&self, ring: &[[Rat; 2]]) -> Option<Vec<[Rat; 3]>> {
        ring.iter()
            .map(|p| {
                let u = scale(&self.x, p[0])?;
                let v = scale(&self.y, p[1])?;
                add(&add(&self.origin, &u)?, &v)
            })
            .collect()
    }
}

/// Translate a ring by the sweep — where the two dimensions have to meet.
pub(crate) fn swept(base: &[[Rat; 3]], sweep: &[Rat; 3]) -> Option<Vec<[Rat; 3]>> {
    base.iter().map(|b| add(b, sweep)).collect()
}

/// Down to the f64 cache, once, at the end. Equal rationals realize to equal bits —
/// that is what makes the two extrusion paths meet — and `to_f64` being correctly
/// rounded is what makes each of them land on the value its decimal named.
pub(crate) fn realize(pts: &[[Rat; 3]]) -> Vec<Point3> {
    pts.iter()
        .map(|p| Point3::from_array([p[0].to_f64(), p[1].to_f64(), p[2].to_f64()]))
        .collect()
}

/// A ring and the ring it sweeps to.
///
/// The two travel together because the prism builder may **reverse** a ring to fix its
/// winding, and the top has to follow — computing the top afterwards, as the builder
/// used to, is only safe while the top is a pure function of the base, which is exactly
/// what stops being true once the sweep is exact.
#[derive(Clone, Debug)]
pub(crate) struct Swept {
    pub base: Vec<Point3>,
    pub top: Vec<Point3>,
    /// The frame's unit normal realized in **world** f64 — the direction "counter-clockwise" is
    /// about, where the vertices above are (a motion frame turns it with them). What the builder
    /// compares the sweep direction against.
    pub normal: Vector3,
    /// The same two rings **before** realization — the truth the f64 above is a cache of.
    /// Carried so the prism's planes can state themselves in rationals ([`SweptRat`]).
    /// Not optional since S6b: the f64 fallback (`Swept::along`) is gone — a prism the exact
    /// arithmetic cannot state is a named reject at the operation, not a point-less build.
    pub exact: SweptRat,
}

/// A swept ring still in rationals: what the prism's faces are, before they are rounded.
///
/// The f64 realization is a cache. These are the numbers the planes come from, and computing a
/// plane here rather than from the realized points is what makes two faces of one plane carry
/// **the same coefficients** — see [`nacre_scalar::canonical_plane_coeffs`].
/// One step of a swept ring, `vertices[i] → vertices[i + 1]`: straight, or an arc around a point
/// of the base plane, counter-clockwise about the frame normal when `ccw`. The `cache` is the f64
/// cylinder the arc's wall will be stored under — realized in world coordinates where the ring's
/// vertices were, so that a motion frame's arc walls sit where its vertices sit.
// The arc carries its f64 cylinder cache beside the straight variant's nothing: a per-step record,
// a handful per ring, read once by the builder — boxing would buy nothing here.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug)]
pub(crate) enum Seg3 {
    Line,
    Arc {
        center: [Rat; 3],
        radius: Rat,
        ccw: bool,
        /// The seam direction the wall is stated with: the frame's `x̂`, or for a whole circle the
        /// direction from the centre to its one vertex (which is then the seam, at `+ref_dir`).
        ref_dir: [Rat; 3],
        cache: nacre_geom::Cylinder,
    },
}

impl Seg3 {
    /// The same step walked the other way.
    pub(crate) fn reversed(&self) -> Seg3 {
        match self {
            Seg3::Line => Seg3::Line,
            Seg3::Arc {
                center,
                radius,
                ccw,
                ref_dir,
                cache,
            } => Seg3::Arc {
                center: *center,
                radius: *radius,
                ccw: !ccw,
                ref_dir: *ref_dir,
                cache: *cache,
            },
        }
    }

    pub(crate) fn is_arc(&self) -> bool {
        matches!(self, Seg3::Arc { .. })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SweptRat {
    pub base: Vec<[Rat; 3]>,
    pub top: Vec<[Rat; 3]>,
    /// `segs[i]` is the step `base[i] → base[i + 1]` (and the same step on the top ring).
    pub segs: Vec<Seg3>,
    /// The frame's unit normal, exact and in the frame's own coordinates — the axis every arc
    /// wall is stated along.
    pub normal: [Rat; 3],
    /// Which way the base ring runs about that normal **as the ring stands in the world**,
    /// decided exactly on the profile's own 2-D coordinates ([`crate::Ring2d::winding_sign`]) —
    /// `Positive` is counter-clockwise — and carried through the frame's motion: a **reflection**
    /// in the chain reverses every loop's sense, so the 2-D reading is multiplied by the chain's
    /// parity ([`nacre_cip::chain_parity`]). Flipped by [`Swept::reversed`].
    ///
    /// ★ Written in the frame's coordinates alone it was wrong for every pad on a mirrored
    /// solid's tilted face: the ring realized clockwise about its world normal while this said
    /// counter-clockwise, and the prism went up inside-out (measured, a generated session —
    /// rotate, mirror, pad; the f64 cross-check in `oriented_ring` is what caught it).
    pub winding: Orient,
    /// ★★★ **Which frame the rationals above are written in.** `None` is the world, which is
    /// every case that existed before tilted faces became exact. `Some` says they are the
    /// coordinates of a plane's own frame, and that this motion is what carries them out — so
    /// every plane derived from them states itself *in that frame* and every vertex is
    /// written against it (its planes record the motion).
    pub motion: Option<Handle<MotionNode>>,
}

impl Swept {
    /// Reverse both rings, keeping them paired.
    pub(crate) fn reversed(mut self) -> Self {
        self.base.reverse();
        self.top.reverse();
        self.exact.base.reverse();
        self.exact.top.reverse();
        // Step `k` of the reversed ring is old step `n − 2 − k` walked backwards.
        let n = self.exact.segs.len();
        let segs: Vec<Seg3> = (0..n)
            .map(|k| self.exact.segs[(2 * n - 2 - k) % n].reversed())
            .collect();
        self.exact.segs = segs;
        self.exact.winding = match self.exact.winding {
            Orient::Positive => Orient::Negative,
            Orient::Negative => Orient::Positive,
            Orient::Zero => Orient::Zero,
        };
        self
    }
}

impl SweptRat {
    /// **Segment `i → i+1`'s wall, as the three points that define it** — `base[i]`, `base[j]`,
    /// `top[i]`, the same three `Plane::through_points` is given, so the exact record and the f64
    /// one describe the plane the same way round.
    ///
    /// ★★★ **There is no coefficient twin of this any more.** There used to be a `wall_plane`
    /// returning `(b − a) × (c − a)` canonicalized, which overflows `i128` far sooner than the
    /// points themselves do (measured — coefficients to 213 bits, ring coordinates to 107) and so
    /// left most walls unnamed. The name is now derived from these points where it is needed, at
    /// whatever precision the derivation takes.
    pub(crate) fn wall_points(&self, i: usize) -> [[Rat; 3]; 3] {
        let j = (i + 1) % self.base.len();
        [self.base[i], self.base[j], self.top[i]]
    }

    /// Three ring points spanning a cap — `base` when `top` is false, the top ring otherwise.
    ///
    /// ★ **The indices are chosen on the f64 realization, and that is not a compromise.** Whichever
    /// three indices come out, the points recorded are the exact rationals; the realization only
    /// answers *"which three are spread out"*, and a ring whose f64 image has no spread has no
    /// plane in the first place (`Plane::through_points` would already have refused it). Picking
    /// exactly, by rational cross products, would risk an `i128` overflow to answer a question
    /// that does not need exactness.
    pub(crate) fn cap_points(&self, top: bool) -> Option<[[Rat; 3]; 3]> {
        let ring = if top { &self.top } else { &self.base };
        if ring.len() < 3 {
            // Fewer than three vertices: an arc names the plane instead — its centre (lifted to
            // this cap), its vertex, and the vertex's radius turned a quarter about the normal.
            // For a circle seamed at `+x̂` that is `[c, c + r·x̂, c + r·ŷ]`, the cylinder
            // primitive's own cap triple.
            let (k, seg) = self.segs.iter().enumerate().find(|(_, s)| s.is_arc())?;
            let Seg3::Arc { center, .. } = seg else {
                return None;
            };
            let lifted = if top {
                add(center, &sub(&self.top[k], &self.base[k])?)?
            } else {
                *center
            };
            let v = ring[k];
            let radial = sub(&v, &lifted)?;
            return Some([lifted, v, add(&lifted, &cross(&self.normal, &radial)?)?]);
        }
        let f = |p: &[Rat; 3]| [p[0].to_f64(), p[1].to_f64(), p[2].to_f64()];
        let (a, b) = (f(&ring[0]), f(&ring[1]));
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        // The widest turn away from `ring[0] → ring[1]`, so a nearly-collinear pair is not chosen
        // when a better one exists.
        let (mut best, mut best_k) = (0.0f64, None);
        for (k, p) in ring.iter().enumerate().skip(2) {
            let c = f(p);
            let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let x = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let n = x[0].abs() + x[1].abs() + x[2].abs();
            if n > best {
                (best, best_k) = (n, Some(k));
            }
        }
        best_k.map(|k| [ring[0], ring[1], ring[k]])
    }
}

/// Every ring of a prism — placed on the plane and swept along it — computed in exact
/// rationals and realized once, at the end. Returns the outer ring and then the holes.
///
/// **Two frames can serve, and `frame` picks which.** `None` means the world: the sketch plane
/// has to lift to exact orthonormal rationals itself, which an axis-aligned face does and a
/// tilted one does not. `Some(motion)` means the plane's **own** frame, where the sketch frame is
/// the identity ([`RatFrame::identity`]) and the rationals below are its `(u, v, w)` — the same
/// arithmetic, on numbers that are exact by construction rather than by luck.
///
/// `None` when there is no exact form to compute in: a frame that is not exactly
/// orthonormal, or i128 overflow in the placement arithmetic. Both mean the same thing
/// to the caller, which is to keep its f64 path. (A dimension outside the decimal window
/// used to be a third reason; the profile constructor now names it before it gets here.)
pub(crate) fn prism_rings(
    model: &Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
    frame: Option<Handle<MotionNode>>,
) -> Option<(Swept, Vec<Swept>)> {
    let f = match frame {
        Some(_) => RatFrame::identity(),
        None => plane.exact()?,
    };
    prism_rings_in(model, f, profile, dist, frame)
}

/// [`prism_rings`] with the rational frame already in hand — for a caller that derived it from a
/// [`crate::SketchFrame`] rather than from a caller's `f64` axes. One implementation, two doors.
pub(crate) fn prism_rings_in(
    model: &Model,
    f: RatFrame,
    profile: &Profile2d,
    dist: f64,
    frame: Option<Handle<MotionNode>>,
) -> Option<(Swept, Vec<Swept>)> {
    let sweep = f.sweep(dist)?;
    // ★★★ **The realization must be the *definition's* own replay, not a second route to the
    // same real number.** A vertex written here is defined against `frame`, and a judge
    // reads that definition back through `replay`; if this rounded the coordinates some other
    // way the two would sit an ulp apart and the invariant that lets a coordinate be checked
    // against its definition would be false. In the world frame there is no motion and `to_f64`
    // *is* the replay, which is why that branch is the one that was always here.
    let chain = match frame {
        Some(leaf) => Some(crate::rotated_vertex::motion_chain(model, leaf)?),
        None => None,
    };
    let out = |pts: &[[Rat; 3]]| -> Option<Vec<Point3>> {
        match &chain {
            None => Some(realize(pts)),
            Some(c) => realize(pts)
                .iter()
                .map(|f| {
                    // ★★★★ **Realize to the frame's f64 first, and define against *that*.**
                    // The rationals here came from the user's written decimals, and `1/10` does
                    // not round-trip through f64 — so a definition holding `1/10` and a
                    // coordinate holding `0.1` would replay to different bits, and the invariant
                    // that lets a coordinate be checked against its definition would be false.
                    // The vertex's base is the f64 the frame coordinate realizes to, exactly as
                    // the world path already keeps only the realized f64 in the point cache.
                    //
                    // ★ Nothing is lost where it matters: the *plane coefficients* are still
                    // computed from the decimal rationals above, so `7.7` and `1.1 + 6.6` name
                    // one plane however their vertices round.
                    let b = crate::rotated_vertex::coord_rat(f.as_array()).ok()?;
                    let q = crate::rotated_vertex::replay(nacre_cip::WitnessPoint::at(b), c)?;
                    Some(Point3::from_array(q.coord()))
                })
                .collect(),
        }
    };
    let normal = f.normal()?;
    // The chain's handedness: a reflection reverses the sense of every loop it carries, and the
    // exact winding below is read on the profile's 2-D coordinates *before* the chain — so it is
    // turned into the world's sense here, once, where the chain is known.
    let parity = chain.as_deref().map_or(1, nacre_cip::chain_parity);
    let ring = |r: &crate::Ring2d| -> Option<Swept> {
        let base = f.ring(r.vertices())?;
        let top = swept(&base, &sweep)?;
        let base_f64 = out(&base)?;
        let normal_world = out(&[add(&base[0], &normal)?])?.pop()? - base_f64[0];
        let winding = match (r.winding_sign()?, parity) {
            (w, 1) => w,
            (Orient::Positive, _) => Orient::Negative,
            (Orient::Negative, _) => Orient::Positive,
            (Orient::Zero, _) => Orient::Zero,
        };
        let segs = r
            .edges()
            .iter()
            .enumerate()
            .map(|(i, seg)| match seg {
                Edge2d::Line => Some(Seg3::Line),
                Edge2d::Arc {
                    center,
                    radius,
                    ccw,
                } => {
                    let c = f.ring(&[*center])?.pop()?;
                    // The seam reference: a whole circle seams at its one vertex, an arc's wall
                    // is stated with the frame's `x̂` so every arc of one circle interns to one
                    // surface.
                    let ref_dir = if r.len() == 1 {
                        scale(
                            &sub(&base[i], &c)?,
                            Rat::new(radius.denom(), radius.numer())?,
                        )?
                    } else {
                        f.x
                    };
                    // The f64 cache, realized where the vertices were realized: directions as
                    // differences of realized points, so a motion frame turns them too.
                    let at = |p: &[Rat; 3]| -> Option<Point3> { out(&[*p])?.pop() };
                    let c_f64 = at(&c)?;
                    let axis = at(&add(&c, &normal)?)? - c_f64;
                    let refd = at(&add(&c, &ref_dir)?)? - c_f64;
                    let cache =
                        nacre_geom::Cylinder::from_axis(c_f64, axis, refd, radius.to_f64())?;
                    Some(Seg3::Arc {
                        center: c,
                        radius: *radius,
                        ccw: *ccw,
                        ref_dir,
                        cache,
                    })
                }
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Swept {
            base: base_f64,
            top: out(&top)?,
            normal: normal_world,
            exact: SweptRat {
                base,
                top,
                segs,
                normal,
                winding,
                motion: frame,
            },
        })
    };
    let outer = ring(profile.outer())?;
    let holes = profile
        .holes()
        .iter()
        .map(ring)
        .collect::<Option<Vec<_>>>()?;
    Some((outer, holes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;

    /// ★ **A pad on a mirrored solid's tilted face stands up, and a pocket goes down.** The face's
    /// sketch frame is carried by a motion chain with a reflection in it, and the exact winding
    /// read on the profile's 2-D coordinates is the *opposite* sense once the chain has carried
    /// the ring into the world — so [`SweptRat::winding`] folds the chain's parity in. Without
    /// that the prism went up inside-out (a generated session found it: rotate, mirror, pad; the
    /// f64 cross-check in `oriented_ring` fired). The control is the same solid without the
    /// mirror. Volumes are the oracle: a 2-cube plus a 1×1×1 pad, minus a 1×1×½ pocket.
    #[test]
    fn a_pad_on_a_mirrored_tilted_face_stands_up() {
        use crate::{OpOutput, Operation, SketchFrame, apply, from_rings};
        use nacre_math::Point2;
        use nacre_scalar::{Angle, Axis, Isometry, Rotation};
        use nacre_topo::Model;

        let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
        let rect = |x0: f64, y0: f64, x1: f64, y1: f64| {
            from_rings(vec![vec![p2(x0, y0), p2(x1, y0), p2(x1, y1), p2(x0, y1)]])
                .unwrap()
                .remove(0)
        };
        // A 2-cube centred on the origin, turned 15° about y (its planes leave the rational
        // world, so a face's sketch lives in the plane's own frame), then mirrored in x or not.
        let build = |mirror: bool| -> (Model, nacre_store::Handle<nacre_topo::Solid>) {
            let mut m = Model::new();
            let frame = SketchFrame::world(&m, Axis::Z);
            let OpOutput::Extrude { solid, .. } = apply(
                &mut m,
                &Operation::Extrude {
                    frame,
                    profile: rect(-1.0, -1.0, 1.0, 1.0),
                    dist: 2.0,
                },
            )
            .unwrap() else {
                unreachable!()
            };
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid,
                    isometry: Isometry::translation([
                        Rat::from_int(0),
                        Rat::from_int(0),
                        Rat::from_int(-1),
                    ]),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid,
                    isometry: Isometry::rotation(Rotation {
                        axis: Axis::Y,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(15)).unwrap(),
                    }),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            if !mirror {
                return (m, solid);
            }
            let OpOutput::Mirror { solid } = apply(
                &mut m,
                &Operation::Mirror {
                    solid,
                    axis: Axis::X,
                    offset: Rat::from_int(0),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            (m, solid)
        };
        // The face that was `x = +1`: its outward normal after the turn (and the mirror).
        let face_toward = |m: &Model, s: nacre_store::Handle<nacre_topo::Solid>, n: [f64; 3]| {
            let shell = m.solids.get(s).outer;
            *m.shells
                .get(shell)
                .faces
                .iter()
                .find(|&&fh| {
                    let f = m.faces.get(fh);
                    let nacre_geom::Surface::Plane(pl) = m.surface(f.surface) else {
                        return false;
                    };
                    let out = pl.normal() * f64::from(f.orientation.sign());
                    (0..3).all(|k| (out.as_array()[k] - n[k]).abs() < 1e-9)
                })
                .expect("the turned +x face")
        };
        let (s15, c15) = 15f64.to_radians().sin_cos();
        let volume = |m: &Model, s| nacre_props::mass_props(m, s).unwrap().volume;
        for (mirror, nx) in [(false, c15), (true, -c15)] {
            for (pocket, want) in [(false, 9.0), (true, 7.5)] {
                let (mut m, solid) = build(mirror);
                let face = face_toward(&m, solid, [nx, 0.0, -s15]);
                let profile = rect(-0.5, -0.5, 0.5, 0.5);
                let out = if pocket {
                    apply(
                        &mut m,
                        &Operation::PocketOnFace {
                            face,
                            profile,
                            dist: 0.5,
                        },
                    )
                } else {
                    apply(
                        &mut m,
                        &Operation::PadOnFace {
                            face,
                            profile,
                            dist: 1.0,
                        },
                    )
                };
                let solid = match out
                    .unwrap_or_else(|e| panic!("mirror {mirror} pocket {pocket}: {e:?}"))
                {
                    OpOutput::PadOnFace { solid, .. } | OpOutput::PocketOnFace { solid, .. } => {
                        solid
                    }
                    other => panic!("{other:?}"),
                };
                m.rebuild_adjacency();
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "mirror {mirror} pocket {pocket}: {:?}",
                    nacre_validate::validate(&m)
                );
                let v = volume(&m, solid);
                assert!(
                    (v - want).abs() < 1e-9,
                    "mirror {mirror} pocket {pocket}: {v} vs {want}"
                );
            }
        }
    }

    fn rotated(deg: f64) -> SketchPlane {
        let (s, c) = deg.to_radians().sin_cos();
        SketchPlane::from_axes(
            Point3::origin(),
            Vector3::from_array([c, s, 0.0]),
            Vector3::from_array([-s, c, 0.0]),
        )
    }

    /// The axes a sketch plane actually gets in the common cases are `{0, ±1}`, which
    /// are exact decimals and exactly orthonormal — so those frames lift.
    #[test]
    fn axis_aligned_frames_lift_and_rotated_ones_decline() {
        assert!(SketchPlane::world_xy().exact().is_some());
        for n in [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]] {
            let p = SketchPlane::from_origin_normal(
                Point3::from_array([2.5, -1.25, 0.0]),
                Vector3::from_array(n),
            )
            .unwrap();
            assert!(p.exact().is_some(), "axis-aligned normal {n:?}");
        }
        // Not a tolerance question: `cos 45°` is irrational, so its decimal is not a
        // unit vector and no amount of precision would make it one.
        for deg in [45.0, 30.0, 1.0, 0.1] {
            assert!(rotated(deg).exact().is_none(), "{deg}°");
        }
        // **And a quarter turn built through f64 trig declines too**, which is worth
        // knowing rather than assuming otherwise: `(90°).to_radians().sin_cos()` gives
        // `cos = 6.1e-17`, not `0`, so the axes are not orthonormal and there is nothing
        // here to recover — the exactness was lost before this module saw the frame.
        // The kernel's exact quarter turns come from `Angle` (Niven), not from `f64::cos`.
        for deg in [90.0, 180.0, 270.0] {
            assert!(rotated(deg).exact().is_none(), "{deg}° through f64 trig");
        }
        // Spelled exactly, the same rotation lifts.
        assert!(rotated(0.0).exact().is_some());
        assert!(
            SketchPlane::from_axes(
                Point3::origin(),
                Vector3::from_array([0.0, 1.0, 0.0]),
                Vector3::from_array([-1.0, 0.0, 0.0]),
            )
            .exact()
            .is_some(),
            "a quarter turn written down rather than computed"
        );
    }

    /// A lifted frame's normal is *exactly* unit, which is the property the sweep needs:
    /// `normal · dist` has length `dist` with no residue to accumulate.
    #[test]
    fn a_lifted_frames_normal_is_exactly_unit() {
        let f = SketchPlane::world_xy().exact().unwrap();
        let n = f.normal().unwrap();
        assert_eq!(n, [Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
        assert_eq!(dot(&n, &n).unwrap(), Rat::from_int(1));
        assert_eq!(f.sweep(7.7).unwrap()[2], Rat::new(77, 10).unwrap());
    }

    /// **The proposition this module exists for.** One sweep of `7.7`, and two sweeps
    /// of `1.1` and `6.6`, must put the top ring on the same points — as rationals and,
    /// after realization, in the same f64 bits.
    #[test]
    fn a_split_sweep_lands_on_the_same_points_as_the_whole_one() {
        let f = SketchPlane::world_xy().exact().unwrap();
        let d = |x: f64| Rat::from_decimal(x).unwrap();
        let ring = f
            .ring(&[
                [d(0.0), d(0.0)],
                [d(3.3), d(0.0)],
                [d(3.3), d(2.2)],
                [d(0.0), d(2.2)],
            ])
            .unwrap();

        let whole = swept(&ring, &f.sweep(7.7).unwrap()).unwrap();
        let lower = swept(&ring, &f.sweep(1.1).unwrap()).unwrap();
        let split = swept(&lower, &f.sweep(6.6).unwrap()).unwrap();

        assert_eq!(whole, split, "exact");
        assert_eq!(realize(&whole), realize(&split), "realized");
        assert_eq!(realize(&whole)[0][2], 7.7);

        // The f64 arithmetic this replaces does not agree with itself.
        assert_ne!(0.0 + 1.1 + 6.6, 0.0 + 7.7);
    }

    /// The escape hatch has to actually escape: a frame that cannot be lifted returns
    /// `None` at the entry point, not a wrong answer further in.
    #[test]
    fn a_frame_that_cannot_be_lifted_declines_before_any_arithmetic() {
        assert!(rotated(45.0).exact().is_none());
        // And so does a sweep distance outside the decimal window (design: i128) — the one
        // decimal lift still performed here. (A profile coordinate out of window no longer
        // reaches this module: `Profile2d`'s constructor names it.)
        let f = SketchPlane::world_xy().exact().unwrap();
        assert!(f.sweep(1e300).is_none());
        // Placement arithmetic that overflows i128 declines the same way: two in-window
        // factors (10^37 each) whose product (10^74) has no exact form to keep.
        let huge = Rat::from_decimal(1e37).unwrap();
        let stretched = RatFrame {
            origin: [Rat::from_int(0); 3],
            x: [huge, Rat::from_int(0), Rat::from_int(0)],
            y: [Rat::from_int(0), huge, Rat::from_int(0)],
        };
        assert!(stretched.ring(&[[huge, huge]]).is_none());
    }
}
