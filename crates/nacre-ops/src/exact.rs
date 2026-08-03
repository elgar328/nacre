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
use nacre_math::{Point2, Point3, Vector3};
use nacre_scalar::Rat;
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

    /// The sweep vector for a signed distance along the normal — the second of the two
    /// arithmetic steps a prism is built from.
    pub(crate) fn sweep(&self, dist: f64) -> Option<[Rat; 3]> {
        scale(&self.normal()?, Rat::from_decimal(dist)?)
    }

    /// A ring's sketch coordinates placed in space: `origin + x·u + y·v`, the first of
    /// the two steps.
    pub(crate) fn ring(&self, ring: &[Point2]) -> Option<Vec<[Rat; 3]>> {
        ring.iter()
            .map(|p| {
                let u = scale(&self.x, Rat::from_decimal(p[0])?)?;
                let v = scale(&self.y, Rat::from_decimal(p[1])?)?;
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
    /// The same two rings **before** realization, when the exact path produced them, plus the
    /// frame's normal. Carried so the prism's planes can state themselves in rationals — see
    /// [`SweptRat`]. `None` on the f64 fallback.
    pub exact: Option<SweptRat>,
}

/// A swept ring still in rationals: what the prism's faces are, before they are rounded.
///
/// The f64 realization is a cache. These are the numbers the planes come from, and computing a
/// plane here rather than from the realized points is what makes two faces of one plane carry
/// **the same coefficients** — see [`nacre_scalar::canonical_plane_coeffs`].
#[derive(Clone, Debug)]
pub(crate) struct SweptRat {
    pub base: Vec<[Rat; 3]>,
    pub top: Vec<[Rat; 3]>,
    /// The frame's exact unit normal. The caps face `∓` this **regardless of the sweep's sign**,
    /// which is why it is carried rather than recovered from `top − base`.
    pub normal: [Rat; 3],
    /// ★★★ **Which frame the rationals above are written in.** `None` is the world, which is
    /// every case that existed before tilted faces became exact. `Some` says they are the
    /// coordinates of a plane's own frame, and that this motion is what carries them out — so
    /// every plane derived from them states itself *in that frame* and every vertex is
    /// `Origin::Moved` against it.
    pub motion: Option<Handle<MotionNode>>,
}

impl Swept {
    /// A ring translated by a sweep vector, in f64 — the fallback for a frame or a
    /// dimension with no exact form, and what the builder did unconditionally before.
    pub(crate) fn along(base: Vec<Point3>, sweep: Vector3) -> Self {
        let top = base.iter().map(|b| *b + sweep).collect();
        Swept {
            base,
            top,
            exact: None,
        }
    }

    /// Reverse both rings, keeping them paired.
    pub(crate) fn reversed(mut self) -> Self {
        self.base.reverse();
        self.top.reverse();
        if let Some(e) = self.exact.as_mut() {
            e.base.reverse();
            e.top.reverse();
        }
        self
    }
}

impl SweptRat {
    /// ★★★ **How a plane built here states its provenance.**
    ///
    /// In the world frame it is `Constructed`: the coefficients are world coefficients and the
    /// face's own f64 triangle is the truth. Inside a plane's frame neither is so — the
    /// coefficients are the frame's, and saying `Constructed` would let the judgment read them as
    /// world coefficients, which is the *silent* half of getting this wrong. `Moved` names the
    /// motion that carries them out, and the witness is the same triangle **in frame
    /// coordinates**, which is where the replay starts.
    pub(crate) fn surface_def(&self, witness: [Point3; 3]) -> nacre_topo::SurfaceDef {
        match self.motion {
            None => nacre_topo::SurfaceDef::Constructed,
            Some(motion) => nacre_topo::SurfaceDef::Moved { witness, motion },
        }
    }

    /// The frame-coordinate f64 image of a base ring point — the witness's own frame, and the
    /// coordinate a vertex built here is defined against.
    pub(crate) fn base_f64(&self, i: usize) -> Point3 {
        realize(&self.base[i..=i])[0]
    }

    /// The same for a top ring point.
    pub(crate) fn top_f64(&self, i: usize) -> Point3 {
        realize(&self.top[i..=i])[0]
    }

    /// The base cap (`−normal`) and the top cap (`+normal`).
    pub(crate) fn cap_planes(&self) -> Option<([Rat; 4], [Rat; 4])> {
        let zero = Rat::from_int(0);
        let neg = [
            zero.checked_sub(self.normal[0])?,
            zero.checked_sub(self.normal[1])?,
            zero.checked_sub(self.normal[2])?,
        ];
        Some((
            nacre_scalar::plane_from_point_normal(neg, *self.base.first()?)?,
            nacre_scalar::plane_from_point_normal(self.normal, *self.top.first()?)?,
        ))
    }

    /// Segment `i → i+1`'s wall — the same three points `Plane::through_points` is given
    /// (`base[i]`, `base[j]`, `top[i]`), so the normal points the same way.
    pub(crate) fn wall_plane(&self, i: usize) -> Option<[Rat; 4]> {
        let j = (i + 1) % self.base.len();
        nacre_scalar::plane_through_points(self.base[i], self.base[j], self.top[i])
    }

    /// The same three points, as the wall plane's **exact witness**.
    ///
    /// ★ Deliberately not gated on [`wall_plane`] succeeding: the coefficients are a product of
    /// two point differences and overflow `i128` far sooner than the points themselves do
    /// (measured — coefficients to 213 bits, ring coordinates to 107). A wall with no rational
    /// coefficients still has an exact description, and this is it.
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
            return None;
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
/// orthonormal, a dimension outside the decimal window, or i128 overflow. All three
/// mean the same thing to the caller, which is to keep its f64 path.
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
    let sweep = f.sweep(dist)?;
    let normal = f.normal()?;
    // ★★★ **The realization must be the *definition's* own replay, not a second route to the
    // same real number.** A vertex written here is `Origin::Moved` against `frame`, and a judge
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
                    // the world path already keeps only the realized f64 in `Origin::Constructed`.
                    //
                    // ★ Nothing is lost where it matters: the *plane coefficients* are still
                    // computed from the decimal rationals above, so `7.7` and `1.1 + 6.6` name
                    // one plane however their vertices round.
                    let b = crate::rotated_vertex::coord_rat(f.as_array()).ok()?;
                    let q = crate::rotated_vertex::replay(nacre_cip::Pt3::at(b), c)?;
                    Some(Point3::from_array(q.coord))
                })
                .collect(),
        }
    };
    let ring = |r: &[Point2]| -> Option<Swept> {
        let base = f.ring(r)?;
        let top = swept(&base, &sweep)?;
        Some(Swept {
            base: out(&base)?,
            top: out(&top)?,
            exact: Some(SweptRat {
                base,
                top,
                normal,
                motion: frame,
            }),
        })
    };
    let outer = ring(profile.outer())?;
    let holes = profile
        .inners()
        .iter()
        .map(|h| ring(h))
        .collect::<Option<Vec<_>>>()?;
    Some((outer, holes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;

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
        let ring = f
            .ring(&[
                Point2::from_array([0.0, 0.0]),
                Point2::from_array([3.3, 0.0]),
                Point2::from_array([3.3, 2.2]),
                Point2::from_array([0.0, 2.2]),
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
        // And so does a dimension outside the decimal window (design: i128).
        let f = SketchPlane::world_xy().exact().unwrap();
        assert!(f.sweep(1e300).is_none());
        assert!(f.ring(&[Point2::from_array([1e300, 0.0])]).is_none());
    }
}
