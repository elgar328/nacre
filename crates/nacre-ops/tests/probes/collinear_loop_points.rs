//! **A loop carries points that do not turn, and the face must still know which way it faces.**
//!
//! Two faces sharing an edge must list the same vertices along it, so a pad that covers part of a
//! face leaves the *neighbouring* walls with points strung along one straight line. Those points
//! cannot be dropped — removing one would leave a T-vertex — so every partially-covering pad or
//! pocket produces loops of this shape. It is a common population, not an exotic one.
//!
//! What used to break on it: `planes::outer_tri` took the first corner whose triangle was not
//! *exactly* flat. Three collinear points span exactly zero area only in exact arithmetic; on
//! rotated coordinates the f64 cancellation leaves ~2⁻⁵³, which passed that gate, and the
//! direction that came back was the rounding rather than the plane — measured 90° off. Downstream
//! that is the face's `orient_sign`, the winding of the exact witness triangle every predicate
//! borrows, and "which side of this plane is material".
//!
//! ★ **Two ingredients are needed to trigger it**, which is why it went unseen: rotation (an
//! axis-aligned loop's collinear cross is exactly `0.0`, and the old rule then skipped it
//! correctly) and the collinear run happening to come first in the loop's half-edge order.
//!
//! ★★ The proposition below is deliberately **not** "every cell builds": some angles may decline
//! for their own good reasons, now or later. It is that no cell comes back carrying a face whose
//! stored plane disagrees with its own outward direction — build correctly, or decline by name.

use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::{Point2, Point3};
use nacre_ops::{OpError, OpOutput, Operation, Profile2d, SketchFrame, apply, face_plane};
use nacre_store::Handle;
use nacre_topo::{Face, Model, Solid};

fn rect(a: f64, b: f64, c: f64, d: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, b]),
        Point2::from_array([c, b]),
        Point2::from_array([c, d]),
        Point2::from_array([a, d]),
    ])
    .expect("a rectangle is a fair profile")
}

fn faces_of(m: &Model, s: Handle<Solid>) -> Vec<Handle<Face>> {
    m.shell(m.solid(s).outer).faces.clone()
}

/// **The proposition, read from outside the kernel.**
///
/// For every face of `s`: the plane the model stores for it, and the direction its own outer loop
/// says is outward, must be parallel. `collect_planes` asserts this in debug and nothing else in
/// the kernel looks at it — so this is written as a plain assertion, which runs in release too.
///
/// The loop normal is taken by the Newell sum rather than from three points, on purpose: this is
/// the *check*, and it must not share its derivation with the code under test.
fn every_face_knows_which_way_it_faces(m: &Model, s: Handle<Solid>, cell: &str) {
    for (k, &fh) in faces_of(m, s).iter().enumerate() {
        let face = m.face(fh);
        let pts: Vec<Point3> = face
            .outer
            .half_edges
            .iter()
            .map(|&he| m.vertex_point(m.he_start(he)))
            .collect();
        let n = pts.len();
        let newell = (0..n).fold(nacre_math::Vector3::zero(), |acc, i| {
            acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
        });
        let Some(loop_dir) = newell.normalize() else {
            panic!("{cell}: face[{k}] has no area at all");
        };
        let plane = face_plane(m, fh).expect("planar");
        let cos = plane.normal().dot(loop_dir).abs();
        assert!(
            cos > 0.5,
            "{cell}: face[{k}]'s stored plane is {cos:.3} of the way to perpendicular against \
             what its own loop says is outward — the loop has {n} points"
        );
    }
}

/// A 1×1 footprint **centred on the face, then shifted by `off`** — the placement idiom the rest
/// of the suite uses (`tests/census.rs`'s `wf` family). A profile stated in raw sketch coordinates
/// lands wherever the face's frame origin happens to be, which for a turned solid is usually off
/// the face entirely: the first draft of this sweep declined 84 of 108 cells with `PadMissesFace`
/// and was measuring nothing.
fn footprint_on(m: &Model, face: Handle<Face>, off: f64) -> Profile2d {
    let sp = face_plane(m, face).expect("planar");
    let d = nacre_props::face_props(m, face).expect("props").centroid - sp.origin();
    let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
    rect(
        cu - 0.5 + off,
        cv - 0.5 + off,
        cu + 0.5 + off,
        cv + 0.5 + off,
    )
}

/// One cell: a 2×2×1 block turned by `deg` about `axis`, a pad that covers part of one face (so
/// the neighbouring walls gain non-turning points), then a second pad on a face of *that* result.
fn cell(axis: Axis, deg: i128, inset: f64) -> Result<(Model, Handle<Solid>, f64, f64), OpError> {
    let mut m = Model::new();
    let world = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: world,
            profile: rect(0.0, 0.0, 2.0, 2.0),
            dist: 1.0,
        },
    )?
    else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();

    let OpOutput::Transform { solid } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: Isometry::rotation(Rotation {
                axis,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("a whole-degree angle"),
            }),
        },
    )?
    else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();

    // First pad: partial cover, so the walls around it inherit points that do not turn.
    let f = faces_of(&m, solid)[0];
    let profile = footprint_on(&m, f, inset);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: f,
            profile,
            dist: 1.0,
        },
    )?
    else {
        unreachable!("pad yields PadOnFace output")
    };
    m.rebuild_adjacency();
    every_face_knows_which_way_it_faces(&m, solid, "after the first pad");
    let before = nacre_props::mass_props(&m, solid).expect("props").volume;

    // Second pad, on a face of the first pad's result — this is where the split edge is read.
    let f = faces_of(&m, solid)[0];
    let profile = footprint_on(&m, f, 0.0);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: f,
            profile,
            dist: 1.0,
        },
    )?
    else {
        unreachable!("pad yields PadOnFace output")
    };
    m.rebuild_adjacency();
    let after = nacre_props::mass_props(&m, solid).expect("props").volume;
    Ok((m, solid, before, after))
}

/// **The sweep.** Three axes × twelve angles × three insets.
///
/// ★ A single hand-picked case would only say "this one works". The population is every
/// partially-covering pad, so the band is what has to be swept — the same shape
/// `rotation_sweep.rs` uses for its own band.
#[test]
fn a_loop_with_points_that_do_not_turn_still_states_its_own_outward_direction() {
    let mut built = 0;
    let mut declined: Vec<(String, OpError)> = Vec::new();

    for axis in [Axis::X, Axis::Y, Axis::Z] {
        for deg in [5i128, 15, 30, 45, 60, 75, 100, 123, 150, 200, 250, 305] {
            for inset in [0.0, 0.25, 0.5] {
                let name = format!("{axis:?}/{deg}deg/inset{inset}");
                match cell(axis, deg, inset) {
                    Ok((m, solid, before, after)) => {
                        built += 1;
                        // The proposition: the result states its own faces correctly.
                        every_face_knows_which_way_it_faces(&m, solid, &name);
                        // And it is a solid, not merely a shape.
                        let bad = nacre_validate::validate(&m);
                        assert!(bad.is_empty(), "{name}: {bad:?}");
                        // ★ A bound from the definition, not an independent answer: a pad never
                        // removes material and cannot add more than the tool prism holds. The
                        // second boss overlaps the first and overhangs the block, so there is no
                        // hand-computable volume to compare against — this catches a grossly
                        // wrong answer and says so rather than implying more.
                        //
                        // ★★ The lower bound is `>=`, not `>`: a footprint that lands entirely
                        // inside material adds *nothing*, and that is a correct answer. Measured
                        // (`Z/100deg/inset0`), and it is why the first draft of this assertion
                        // failed on a cell that was right.
                        assert!(
                            after >= before - 1e-9 && after <= before + 1.0 + 1e-9,
                            "{name}: {before} -> {after} is outside what one 1x1x1 tool can add"
                        );
                    }
                    // A decline is a legitimate answer — the kernel says so by name. What must
                    // never happen is the middle case, and that is what the `Ok` arm asserts.
                    Err(e) => declined.push((name, e)),
                }
            }
        }
    }

    // ★★ Measured with the `outer_tri` fix removed: 32 cells built and **14 of them came back
    // with a volume a pad cannot produce** — `Z/5deg/inset0.5` went 5.0 → 3.21, a boss that ate
    // 1.79 of material, with `validate` clean. That is what this file is for, and it is why the
    // volume bound above is not decoration: the face-normal check cannot see it (its own normal
    // is the Newell sum, which was never the broken part), so the bound is the instrument that
    // moved.
    //
    // ★ **All 108 cells build.** They did not always: the X- and Y-turned cells (72 of them) used
    // to decline `PadMissesFace`, because `face_plane` reported a frame point-symmetric to the
    // one the pad realized in — the flip half-turn was written about `v̂` in the report and about
    // `û` in the realization — so a footprint centred on the face through `face_plane`'s own
    // coordinates was built on the opposite side, outside the face. The report reads the
    // realization itself now, and a decline reappearing here is news, not noise.
    assert_eq!(
        built, 108,
        "cells stopped building — a decline reappeared in this sweep: {declined:?}"
    );
}
