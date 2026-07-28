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

// Built and tested here before `extrude`/`pad`/`pocket` call it, so that wiring it up
// is a commit whose only visible effect is the coordinates moving — which is the one
// thing that needs isolated blame. **The next commit wires it and removes this.**
#![allow(dead_code)]

use crate::ops::SketchPlane;
use nacre_math::{Point2, Point3};
use nacre_scalar::Rat;

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
            origin: lift(self.origin.as_array())?,
            x: lift(self.x_axis.as_array())?,
            y: lift(self.y_axis.as_array())?,
        };
        let one = Rat::from_int(1);
        let zero = Rat::from_int(0);
        (dot(&f.x, &f.x)? == one && dot(&f.y, &f.y)? == one && dot(&f.x, &f.y)? == zero)
            .then_some(f)
    }
}

impl RatFrame {
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

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Vector3;

    fn rotated(deg: f64) -> SketchPlane {
        let (s, c) = deg.to_radians().sin_cos();
        SketchPlane {
            origin: Point3::origin(),
            x_axis: Vector3::from_array([c, s, 0.0]),
            y_axis: Vector3::from_array([-s, c, 0.0]),
        }
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
            SketchPlane {
                origin: Point3::origin(),
                x_axis: Vector3::from_array([0.0, 1.0, 0.0]),
                y_axis: Vector3::from_array([-1.0, 0.0, 0.0]),
            }
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
