use super::*;
use nacre_math::Vector3;

/// A 2-cube centred on the origin, turned 15° about y (its planes leave the rational world, so a
/// face's sketch lives in the plane's own frame), then mirrored in x or not — and the face that was
/// `x = +1`, found by its outward normal after the turn (and the mirror).
fn turned_cube_face(
    mirror: bool,
) -> (
    nacre_topo::Model,
    nacre_store::Handle<nacre_topo::Solid>,
    nacre_store::Handle<nacre_topo::Face>,
) {
    use crate::{OpOutput, Operation, SketchFrame, apply};
    use nacre_exact::{Angle, Axis, Isometry, Rotation};
    use nacre_topo::Model;

    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(-1.0, 1.0),
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
    let OpOutput::Transform { mut solid } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Y,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(15)).unwrap(),
            }),
        },
    )
    .unwrap() else {
        unreachable!()
    };
    if mirror {
        let OpOutput::Mirror { solid: mirrored } = apply(
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
        solid = mirrored;
    }
    let (s15, c15) = 15f64.to_radians().sin_cos();
    let n = [if mirror { -c15 } else { c15 }, 0.0, -s15];
    let shell = m.solid(solid).outer;
    let face = *m
        .shell(shell)
        .faces
        .iter()
        .find(|&&fh| {
            let f = m.face(fh);
            let nacre_geom::Surface::Plane(pl) = m.surface_cache(f.surface) else {
                return false;
            };
            let out = pl.normal() * f64::from(f.orientation.sign());
            (0..3).all(|k| (out.as_array()[k] - n[k]).abs() < 1e-9)
        })
        .expect("the turned +x face");
    (m, solid, face)
}

/// The square `[lo, hi]²` as a profile.
fn square(lo: f64, hi: f64) -> crate::Profile2d {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    crate::from_rings(vec![vec![p2(lo, lo), p2(hi, lo), p2(hi, hi), p2(lo, hi)]])
        .unwrap()
        .remove(0)
}

/// ★ **A pad on a mirrored solid's tilted face stands up, and a pocket goes down.** The face's
/// sketch frame is carried by a motion chain with a reflection in it, and the exact winding
/// read on the profile's 2-D coordinates is the *opposite* sense once the chain has carried
/// the ring into the world — so [`SweptRat::winding`] folds the chain's parity in. Without
/// that the prism went up inside-out (a generated session found it: rotate, mirror, pad; the
/// f64 cross-check in `oriented_ring` fired). The control is the same solid without the
/// mirror. Volumes are the oracle: a 2-cube plus a 1×1×1 pad, minus a 1×1×½ pocket.
#[test]
fn a_pad_on_a_mirrored_tilted_face_stands_up() {
    use crate::{OpOutput, Operation, apply};
    let volume = |m: &nacre_topo::Model, s| nacre_props::mass_props(m, s).unwrap().volume;
    for mirror in [false, true] {
        for (pocket, want) in [(false, 9.0), (true, 7.5)] {
            let (mut m, _, face) = turned_cube_face(mirror);
            let profile = square(-0.5, 0.5);
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
                OpOutput::PadOnFace { solid, .. } | OpOutput::PocketOnFace { solid, .. } => solid,
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

/// ★ **A circle swept on a mirrored solid's tilted face stands up too.** The same frame as
/// [`a_pad_on_a_mirrored_tilted_face_stands_up`], a disk of radius ½ swept 1 on it. A circle's
/// ring has one vertex, so its direction is its arc's turn, not a vertex order — and the arc's
/// `ccw` must be read in the world as the winding is ([`Seg3::Arc`]). The control is the unmirrored
/// face; the oracle is `validate` and the volume π/4.
#[test]
fn a_circle_on_a_mirrored_tilted_face_stands_up() {
    use crate::{OpOutput, Operation, apply};
    for mirror in [false, true] {
        let (mut m, _, face) = turned_cube_face(mirror);
        let frame = crate::face_sketch_frame(&m, face).expect("the face's frame");
        let profile = crate::from_paths(vec![
            crate::Ring2d::circle(nacre_math::Point2::from_array([0.0, 0.0]), 0.5)
                .expect("a circle"),
        ])
        .expect("a disk")
        .remove(0);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 1.0,
            },
        )
        .unwrap_or_else(|e| panic!("mirror {mirror}: {e:?}")) else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "mirror {mirror}: {vs:?}");
        let v = nacre_props::mass_props(&m, solid).unwrap().volume;
        let want = std::f64::consts::PI / 4.0;
        assert!((v - want).abs() < 1e-9, "mirror {mirror}: {v} vs {want}");
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
    // decimal lift still performed here. (A profile coordinate out of window does not
    // reach this module: `Profile2d`'s constructor names it.)
    let f = SketchPlane::world_xy().exact().unwrap();
    assert!(f.sweep(1e300).is_none());
    // Placement arithmetic that overflows i128 declines the same way: two in-window
    // factors (10^37 each) whose product (10^74) has no exact form to keep.
    let huge = Rat::from_decimal(1e37).unwrap();
    let stretched = RatFrame {
        origin: [Rat::from_int(0); 3],
        x: [huge, Rat::from_int(0), Rat::from_int(0)],
        y: [Rat::from_int(0), huge, Rat::from_int(0)],
        parity: 1,
    };
    assert!(stretched.ring(&[[huge, huge]]).is_none());
}
