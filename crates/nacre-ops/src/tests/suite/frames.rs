//! Sketch frames and tilted faces: frame axes, the pole convention, pads and pockets on slanted
//! walls.

use super::*;

/// ★★★★ **Padding the same footprint twice must not leave zero-area faces.**
///
/// It did, and `validate` reported nothing. Padding `1.1` and then `6.6` on the top of a
/// cuboid left two faces of area `2.2e-16` on the plane `z = 2.1`, whose long edges sat one
/// ulp apart (`-0.6` against `-0.6000000000000001`).
///
/// The ulp came from the sketch frame's **origin**. `face_frame` takes it from the face's area
/// centroid, and `ring_area_centroid` was rounding twice more than it needed to, so the first
/// pad's top face reported its centre as `-1.11e-16` instead of `0`. The second pad then placed
/// the same profile one ulp away from where the first had placed it, and the kernel — correctly
/// — built the one-ulp-wide faces that answer describes.
///
/// ★ `SketchPlane::exact()` does not catch this: it checks only that the axes are orthonormal,
/// never the origin.
#[test]
fn padding_one_footprint_twice_leaves_no_zero_area_face() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let rect = || {
        Profile2d::polygon(vec![p(-1.0, -0.6), p(1.0, -0.6), p(1.0, 0.6), p(-1.0, 0.6)]).unwrap()
    };
    let mut m = Model::new();
    let mut solid = m.add_cuboid(
        Point3::from_array([-2.0, -2.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();
    // The topmost face pointing up. Everything here is axis-aligned, so this is unambiguous.
    let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        let shell = m.solid(s).outer;
        *m.shell(shell)
            .faces
            .iter()
            .max_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    crate::ops::face_plane(m, f)
                        .ok()
                        .filter(|sp| sp.normal().as_array()[2] > 0.5)
                        .map(|sp| sp.origin.as_array()[2])
                        .unwrap_or(f64::NEG_INFINITY)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face")
    };
    for dist in [1.1, 6.6] {
        let face = top(&m, solid);
        let OpOutput::PadOnFace { solid: out, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face,
                profile: rect(),
                dist,
            },
        )
        .expect("pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid = out;
    }
    let shell = m.solid(solid).outer;
    let degenerate: Vec<_> = m
        .shell(shell)
        .faces
        .iter()
        .filter_map(|&f| {
            let area = nacre_props::face_props(&m, f).ok()?.area;
            (area < 1e-9).then_some((f, area))
        })
        .collect();
    assert!(
        degenerate.is_empty(),
        "zero-area faces survived: {degenerate:?}"
    );
    // A plain box with one rib on top: 6 + 5 walls/cap, no leftovers from the seam.
    assert_eq!(m.shell(shell).faces.len(), 11);
}

/// ★★★★★ **The target: two ways of reaching one height land on one plane, far from the
/// origin.**
///
/// The frame's origin is where the drift used to enter — `face_frame` took it from the face's
/// area centroid, computed in f64 from the face's own vertices, and `exact.rs` lifted that as
/// truth. Padding the same footprint `1.1` then `6.6` put the second profile an ulp from the
/// first and left faces of area `2.2e-16`.
///
/// Placed **far from the world origin**, because that is where the projected origin is least
/// like the old centroid — if anything about the new rule were fragile with distance, a
/// hundred units of it would show here.
#[test]
fn two_routes_to_one_height_share_a_plane_far_from_the_origin() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    // On a lid, frame coordinates are world x and y — the origin is the world origin projected
    // onto the plane and the axes are `u = +x̂`, `v = +ŷ`. So this ring is world
    // x ∈ [99, 101], y ∈ [99.4, 100.6].
    let rect = || {
        Profile2d::polygon(vec![
            p(99.0, 100.6),
            p(99.0, 99.4),
            p(101.0, 99.4),
            p(101.0, 100.6),
        ])
        .unwrap()
    };
    let top = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        let shell = m.solid(s).outer;
        *m.shell(shell)
            .faces
            .iter()
            .max_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    crate::ops::face_plane(m, f)
                        .ok()
                        .filter(|sp| sp.normal().as_array()[2] > 0.5)
                        .map(|sp| sp.origin.as_array()[2])
                        .unwrap_or(f64::NEG_INFINITY)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face")
    };
    let build = |dists: &[f64]| -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let mut solid = m.add_cuboid(
            Point3::from_array([98.0, 98.0, 0.0]),
            Point3::from_array([102.0, 102.0, 1.0]),
        );
        m.rebuild_adjacency();
        for &dist in dists {
            let face = top(&m, solid);
            let OpOutput::PadOnFace { solid: out, .. } = apply(
                &mut m,
                &Operation::PadOnFace {
                    face,
                    profile: rect(),
                    dist,
                },
            )
            .expect("pad") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            solid = out;
        }
        (m, solid)
    };
    let (m1, one) = build(&[7.7]);
    let (m2, two) = build(&[1.1, 6.6]);
    // The two boss tops are the same plane, to the bit.
    let z = |m: &Model, s| {
        crate::ops::face_plane(m, top(m, s))
            .unwrap()
            .origin
            .as_array()[2]
    };
    assert_eq!(
        z(&m1, one),
        z(&m2, two),
        "7.7 and 1.1+6.6 disagree on the cap plane"
    );
    // And the two-step route left nothing degenerate behind.
    let shell = m2.solid(two).outer;
    let degenerate: Vec<_> = m2
        .shell(shell)
        .faces
        .iter()
        .filter_map(|&f| {
            let area = nacre_props::face_props(&m2, f).ok()?.area;
            (area < 1e-9).then_some((f, area))
        })
        .collect();
    assert!(
        degenerate.is_empty(),
        "zero-area faces survived: {degenerate:?}"
    );
    assert_eq!(
        m1.shell(m1.solid(one).outer).faces.len(),
        m2.shell(shell).faces.len(),
        "the two routes did not build the same solid"
    );
}

/// ★★★ **One plane, one origin** — even when two faces of it were made by different operations.
///
/// This is what the area centroid could not promise: it was a property of the *face*, so a
/// boolean that reshaped one face moved its sketch origin away from its coplanar neighbour's.
/// The projection is a property of the plane, and surfaces are interned, so the two cannot
/// disagree. Only the origin is asserted — the axes follow the face's `Orientation`, which is
/// today's behaviour and a separate question.
#[test]
fn two_faces_of_one_plane_share_a_sketch_origin() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    // A notch out of one end: the lid `z = 1` becomes two faces of the same plane.
    let cutter = m.add_cuboid(
        Point3::from_array([0.8, -1.0, 0.5]),
        Point3::from_array([1.2, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let r = boolean_one(&mut m, BoolKind::Cut, a, cutter).expect("cut");
    m.rebuild_adjacency();
    let lids: Vec<_> = m
        .shell(m.solid(r).outer)
        .faces
        .iter()
        .filter(|&&f| {
            crate::ops::face_plane(&m, f)
                .is_ok_and(|sp| sp.normal().as_array()[2] > 0.5 && sp.origin.as_array()[2] == 1.0)
        })
        .copied()
        .collect();
    assert_eq!(lids.len(), 2, "the notch should leave two lid faces");
    let o = |f| crate::ops::face_plane(&m, f).unwrap().origin.as_array();
    assert_eq!(
        o(lids[0]),
        o(lids[1]),
        "coplanar faces disagree on the origin"
    );
    assert_eq!(
        o(lids[0]),
        [0.0, 0.0, 1.0],
        "the lid's origin is the world origin projected"
    );
}

/// ★★★★★ **The convention, face by face.** This table *is* the rule — every property below is a
/// consequence of it, and pinning the consequences without pinning the table would let a
/// different rule that happens to satisfy them slip in.
#[test]
fn the_six_axis_directions_get_the_frames_the_convention_names() {
    let v = |a: [f64; 3]| Vector3::from_array(a);
    for (n, u, w) in [
        // The lid: this is the row that must equal `SketchPlane::world_xy`.
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        // Every wall: `v` is +ẑ.
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    ] {
        let (gu, gv) = crate::ops::frame_axes(v(n)).expect("a unit normal has frame axes");
        assert_eq!(gu.as_array(), u, "u for normal {n:?}");
        assert_eq!(gv.as_array(), w, "v for normal {n:?}");
        // ★ And the axes stay exactly representable, so the rational construction path still
        // fires — losing that would drop every axis-aligned model to f64 silently.
        let plane = SketchPlane::from_axes(Point3::origin(), gu, gv);
        assert!(plane.exact().is_some(), "exact path lost for normal {n:?}");
    }
}

/// The lid's frame and [`SketchPlane::world_xy`] name the same plane, so they must name it the
/// same way. They did not: `any_perpendicular` gave the lid `u = −ŷ, v = +x̂`, ninety degrees
/// round, and the kernel carried both spellings at once.
#[test]
fn a_lid_gets_the_same_frame_as_the_world_xy_plane() {
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let (u, v) = crate::ops::frame_axes(up).unwrap();
    let w = SketchPlane::world_xy();
    assert_eq!(u.as_array(), w.x_axis.as_array());
    assert_eq!(v.as_array(), w.y_axis.as_array());
}

/// **On anything but a horizontal face, `v` points up.** `u = ẑ × n` is horizontal, so
/// `v·ẑ = 1 − n_z² > 0` whenever `n` is not vertical — the reason a sketch on a wall has "up"
/// where a person expects it. Checked on tilts the axis-aligned table cannot reach.
#[test]
fn every_non_horizontal_face_has_its_v_pointing_up() {
    for n in [
        [1.0, 1.0, 0.0],
        [0.6, 0.0, 0.8],
        [-0.3, 0.5, -0.81],
        [0.0, 1.0, 0.001],
        [7.0, -13.0, 5.0],
    ] {
        let n = Vector3::from_array(n).normalize().unwrap();
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        assert_eq!(u.as_array()[2], 0.0, "u must be horizontal for {n:?}");
        assert!(v.as_array()[2] > 0.0, "v points down for {n:?}");
    }
}

/// `u ⊥ n`, `|u| = 1`, and `(u, v, n)` right-handed — for the tilted normals too, where the
/// table above says nothing.
#[test]
fn the_reference_axis_is_a_unit_normal_perpendicular() {
    for n in [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 1.0, 1.0],
        [-2.0, 0.5, 3.25],
        [1e-9, 0.0, 1.0],
    ] {
        let n = Vector3::from_array(n).normalize().unwrap();
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        assert!((u.norm() - 1.0).abs() < 1e-15, "|u| for {n:?}");
        assert!(u.dot(n).abs() < 1e-15, "u·n for {n:?}");
        assert!((u.cross(v) - n).norm() < 1e-15, "handedness for {n:?}");
    }
}

/// ★★ **The jump at the poles is intended, not a bug to be fixed later.**
///
/// No continuous tangent frame exists on the sphere, so some set of normals must jump; this
/// convention spends that budget on the two poles and nowhere else. A normal a billionth off
/// vertical takes the other branch and lands ninety degrees away — pinned here so the next
/// reader can see it was chosen.
#[test]
fn the_frame_jumps_at_the_poles_and_that_is_the_deal() {
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let tilted = Vector3::from_array([1e-9, 0.0, 1.0]).normalize().unwrap();
    assert_eq!(
        crate::ops::frame_axes(up).unwrap().0.as_array(),
        [1.0, 0.0, 0.0]
    );
    assert_eq!(
        crate::ops::frame_axes(tilted).unwrap().0.as_array(),
        [0.0, 1.0, 0.0]
    );
}

/// No `-0.0` reaches a caller. It compares equal to `0.0` and lifts to the same rational, so
/// this is presentation only — but a frame printed as `[-0.0, 1.0, 0.0]` reads like a defect.
#[test]
fn no_axis_component_is_negative_zero() {
    for n in [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ] {
        let n = Vector3::from_array(n);
        let (u, v) = crate::ops::frame_axes(n).unwrap();
        for c in u.as_array().iter().chain(v.as_array().iter()) {
            assert!(
                !(*c == 0.0 && c.is_sign_negative()),
                "negative zero in {n:?}'s frame"
            );
        }
    }
}

/// The zero vector has no frame — the only `None`.
#[test]
fn the_zero_vector_has_no_frame_axes() {
    assert!(crate::ops::frame_axes(Vector3::zero()).is_none());
}

/// ★★★★★ **A boss on a tilted face, and another beside it — `pad` used to lose the second one.**
///
/// The prism's far cap and the first boss's cap are one plane, and the boolean says so: its
/// classes are decided with evidence, and on a tilted face that evidence is a composed-rotation
/// proof or a coincidence within the limit, never a handle match — the cap's surface has no
/// rational coefficients to intern by, so it is minted fresh.
///
/// `find_face_coplanar_with` then had to guess which surface the class had collapsed to, from
/// handles and an exact `plane_side` on f64 points. Both miss: the survivor carries the *other*
/// operand's surface, and the two f64 planes sit `1.8e-15` apart. `pad` turned "cap not found"
/// into a hard error and **threw away a correct solid** — the volume was already right.
///
/// Now it asks the boolean instead.
#[test]
fn a_second_boss_on_a_tilted_face_keeps_its_cap() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let mut s = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    m.rebuild_adjacency();
    // Two turns, so the face's normal is off every world axis and its frame is not exact.
    for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::rotation(Rotation {
                    axis,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                }),
            },
        )
        .expect("turn") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
    }
    // Where the original +z went, so the face can be picked without a heuristic tie.
    let (sy, cy) = (53f64).to_radians().sin_cos();
    let (sz, cz) = (17f64).to_radians().sin_cos();
    let up = Vector3::from_array([cz * sy, sz * sy, cy]);
    // The **original** tilted face, not a boss raised on it: lowest along `up` among the faces
    // that point that way. Taking the highest would stack the second boss on the first.
    let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        *m.shell(m.solid(s).outer)
            .faces
            .iter()
            .filter(|&&f| crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99))
            .min_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .expect("a face along up")
    };
    // The frame is a function of the plane, so it survives the face being reshaped by the first
    // pad — both columns are placed from one reading.
    let (cu, cv) = {
        let f = facing(&m, s);
        let sp = crate::ops::face_plane(&m, f).expect("planar");
        let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
        (d.dot(sp.x_axis), d.dot(sp.y_axis))
    };
    let before = nacre_props::mass_props(&m, s).unwrap().volume;
    let mut caps = Vec::new();
    for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
        let f = facing(&m, s);
        let profile = Profile2d::polygon(vec![
            p(cu + lo, cv - 0.5),
            p(cu + hi, cv - 0.5),
            p(cu + hi, cv + 0.5),
            p(cu + lo, cv + 0.5),
        ])
        .unwrap();
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: f,
                profile,
                dist: 7.7,
            },
        )
        .expect("a boss on a tilted face keeps its cap") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
        caps.push(top_face);
    }
    // Each column is 0.6 × 1.0 × 7.7 = 4.62 of material.
    let after = nacre_props::mass_props(&m, s).unwrap().volume;
    assert!(
        (after - before - 2.0 * 4.62).abs() < 1e-9,
        "volume {before} -> {after}"
    );
    // ★ And the handle it returned is the boss top, not some other face that happened to pass:
    // area 0.6, outward along the face's normal.
    for cap in caps {
        let props = nacre_props::face_props(&m, cap).expect("a planar cap");
        assert!((props.area - 0.6).abs() < 1e-9, "cap area {}", props.area);
        assert!(
            props.normal.is_some_and(|n| n.dot(up) > 0.99),
            "cap faces the wrong way"
        );
    }
    assert!(nacre_validate::validate(&m).is_empty());
}

//
// ★★★ These pin **today's** behaviour, not a wish. The suite had almost none of them, so the
// stage that makes tilted frames exact would otherwise be built with nothing to measure against.
// Each says which part 2b is expected to change and which part must not move.

/// A profile edge running `(1, 2)` sweeps a wall whose normal is `(2, 1, 0)` — and **no
/// rational-degree rotation reaches it** (the angle is `atan(1/2)`). That is the case
/// `Motion::Frame` exists for.
///
/// What holds today and must keep holding:
///
/// * the wall's **world coefficients are rational** — `[2, 1, 0, −10]`, so a frame built on it
///   has exact data to derive from;
/// * its sketch frame follows the arbitrary-axis convention — origin at the world origin's
///   projection `(4, 2, 0) = (10/5)·(2,1,0)`, `u = ẑ × n` normalized, `v = +ẑ`;
/// * a pad on it works — **through `Motion::Frame`**, since stage 2 (an earlier line here
///   said "through the f64 path", which stopped being true then).
///
/// What the final assertion pins is narrower than the test's old name suggests: the
/// *reported* `SketchPlane`'s world axes are irrational, so `exact()` is `None` — that
/// is about the world lift, not about the frame road the
/// operation actually takes. ★ The axes must **not** move — this plane's own frame is the
/// world, so `ẑ × n` is the same vector before and after.
#[test]
fn a_sketch_on_a_prism_side_wall_takes_the_f64_path_today() {
    let (m, wall) = prism_with_a_slanted_wall();
    let sp = crate::ops::face_plane(&m, wall).expect("planar");
    let c = m
        .surface_name
        .get(&m.face(wall).surface)
        .expect("a world-frame wall has a name")
        .narrow()
        .expect("a world-frame wall's name is narrow");
    assert_eq!(c.map(|r| r.to_f64()), [2.0, 1.0, 0.0, -10.0]);
    assert_eq!(
        sp.origin.as_array(),
        [4.0, 2.0, 0.0],
        "the projected origin"
    );
    assert_eq!(
        sp.y_axis.as_array(),
        [0.0, 0.0, 1.0],
        "v points up on a wall"
    );
    assert!(
        (sp.x_axis - Vector3::from_array([-1.0, 2.0, 0.0]).normalize().unwrap()).norm() < 1e-15,
        "u is ẑ × n normalized, got {:?}",
        sp.x_axis.as_array()
    );
    // ★ The gap 2b closes: the frame is not exact, so nothing built here records coefficients.
    assert!(
        sp.exact().is_none(),
        "a tilted frame has no rational form today"
    );
}

/// The same wall, actually used: a boss on it comes out right through the f64 path. Pinned so
/// that making the frame exact cannot change the **answer**, only how it is recorded.
#[test]
fn a_boss_on_a_slanted_wall_is_correct_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PadOnFace { solid, top_face } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 1.0,
        },
    )
    .expect("pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // Base prism 15 × 3 = 45, plus a 1 × 1 × 1 boss.
    let props = nacre_props::mass_props(&m, solid).unwrap();
    assert!(
        (props.volume - 46.0).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(
        (nacre_props::face_props(&m, top_face).unwrap().area - 1.0).abs() < 1e-9,
        "the boss top is 1 × 1"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// ★★★ **A sketch on a wall raised from a sketch on a wall** — the nesting `Motion::Frame` was
/// redesigned for. The second wall's *world* normal is irrational, so its frame cannot be
/// written down by naming a normal; only by naming the plane.
///
/// It works today, through f64. Pinned because nesting is where a frame that names its plane
/// by handle must terminate its recursion.
#[test]
fn a_sketch_on_a_wall_raised_from_a_slanted_wall_works_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PadOnFace { solid, top_face } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 1.0,
        },
    )
    .expect("pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // A side face of that boss: not the cap, not on the original wall's plane.
    let cap_n = nacre_props::face_props(&m, top_face)
        .unwrap()
        .normal
        .unwrap();
    let side = *m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            f != top_face
                && nacre_props::face_props(&m, f).is_ok_and(|p| {
                    (p.area - 1.0).abs() < 1e-9
                        && p.normal.is_some_and(|n| n.dot(cap_n).abs() < 0.5)
                })
        })
        .expect("a boss side face");
    let profile = centred_on(&m, side, 0.2);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: side,
            profile,
            dist: 0.5,
        },
    )
    .expect("nested pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let props = nacre_props::mass_props(&m, solid).unwrap();
    // The nested boss is 0.4 × 0.4 × 0.5 = 0.08.
    assert!(
        (props.volume - 46.08).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A pocket on the slanted wall — the other face-based operation, so the sweep runs inward.
#[test]
fn a_pocket_in_a_slanted_wall_is_correct_today() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let profile = centred_on(&m, wall, 0.5);
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: wall,
            profile,
            dist: 0.5,
        },
    )
    .expect("pocket") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let props = nacre_props::mass_props(&m, solid).unwrap();
    assert!(
        (props.volume - 44.5).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A `2·half` square centred on `face`, in that face's own sketch frame. The frame is a
/// function of the plane, so reading it once is enough even if the face is reshaped later.
/// ★★★★★ **The payoff, and its limit — both measured.**
///
/// Two bosses of the same height on one tilted face used to be two plane records that agreed
/// only if their f64 coefficients happened to. Written in the plane's own frame they are both
/// `w = 7.7`, and `SurfaceKey` is `(coefficients, motion)` — so they are **one
/// `Handle<Surface>` at construction**, before anything is compared. And every face of the
/// result states itself exactly, where before a tilted sketch recorded nothing at all.
///
/// ★★★ **What this does *not* buy, stated plainly**: `7.7` against `1.1 + 6.6` — the target
/// the plan named. Stacking sketches the second boss on the **first boss's cap**, which is a
/// different plane and therefore a different frame, so the two caps come out `w = 7.7` and
/// `w = 6.6` — two exact descriptions of one plane that `SurfaceKey` cannot equate. That is
/// not a regression (before this they had no descriptions at all, and the merge still happens
/// through the judge), but the plan's headline claim only holds for sketches sharing a frame,
/// which is what this pins instead.
#[test]
fn two_bosses_on_one_tilted_face_share_a_cap_plane_by_name() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let __w0 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w0,
            profile: Profile2d::polygon(vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)])
                .unwrap(),
            dist: 3.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let mut s = solid;
    for (axis, deg) in [(Axis::Y, 53i128), (Axis::Z, 17)] {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::rotation(Rotation {
                    axis,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
                }),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
    }
    let (sy, cy) = (53f64).to_radians().sin_cos();
    let (sz, cz) = (17f64).to_radians().sin_cos();
    let up = Vector3::from_array([cz * sy, sz * sy, cy]);
    let facing = |m: &Model, s: Handle<Solid>| -> Handle<Face> {
        *m.shell(m.solid(s).outer)
            .faces
            .iter()
            .filter(|&&f| crate::ops::face_plane(m, f).is_ok_and(|sp| sp.normal().dot(up) > 0.99))
            .min_by(|&&a, &&b| {
                let h = |f: Handle<Face>| {
                    (nacre_props::face_props(m, f).unwrap().centroid - Point3::origin()).dot(up)
                };
                h(a).partial_cmp(&h(b)).unwrap()
            })
            .unwrap()
    };
    let mut caps = Vec::new();
    for (lo, hi) in [(-1.0f64, -0.4f64), (0.4, 1.0)] {
        let f = facing(&m, s);
        let sp = crate::ops::face_plane(&m, f).unwrap();
        let d = nacre_props::face_props(&m, f).unwrap().centroid - sp.origin;
        let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
        let profile = Profile2d::polygon(vec![
            p(cu + lo, cv - 0.5),
            p(cu + hi, cv - 0.5),
            p(cu + hi, cv + 0.5),
            p(cu + lo, cv + 0.5),
        ])
        .unwrap();
        let OpOutput::PadOnFace { solid, top_face } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: f,
                profile,
                dist: 7.7,
            },
        )
        .expect("boss") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        s = solid;
        let su = m.face(top_face).surface;
        assert_eq!(
            m.surface_name
                .get(&su)
                .and_then(|n| n.narrow())
                .map(|c| c.map(|r| r.to_f64())),
            Some([0.0, 0.0, 10.0, -77.0]),
            "★ a cap raised in a frame records `w = 7.7` there, exactly"
        );
        caps.push(su);
    }
    assert_eq!(caps[0], caps[1], "one plane, one handle, by name");
    // Rule audit: how many of the result's face surfaces state themselves exactly.
    let sh = m.solid(s).outer;
    let (mut with, mut without) = (0, 0);
    for &f in &m.shell(sh).faces {
        let su = m.face(f).surface;
        if m.surface_name.contains_key(&su) {
            with += 1
        } else {
            without += 1
        }
    }
    assert_eq!(
        (with, without),
        (16, 0),
        "★ every face of a twice-turned, twice-bossed result states itself exactly"
    );
}

/// ★★★★★ **A plane states itself exactly, and its f64 axes are the realization of that.**
///
/// This is the whole point of closing the struct: `(1, 1, 1)` is coefficients `[1, 1, 1, 0]`,
/// three integers, while the unit axes derived from it square to `0.9999999999999999…`. The
/// old API stored only the axes and threw the normal away, so nothing exact survived the door.
#[test]
fn a_named_plane_records_what_its_caller_stated() {
    // The canonical name is *derived* from the definition's points now; reading it back is
    // how the old coefficient assertions keep their meaning.
    let f = |d: &PlaneDef| {
        let p = d.points();
        nacre_scalar::plane_name_exact(p[0], p[1], p[2])
            .expect("a definition names a plane")
            .narrow()
            .expect("these fixtures are narrow")
            .map(|r| r.to_f64())
    };
    // The three world planes, with the axes the script layer documents.
    for (p, want, u) in [
        (
            SketchPlane::world_xy(),
            [0.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ),
        (
            SketchPlane::world_yz(),
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ),
        // ★ `ẑ × n` would give `−x̂` here; a named plane says `+u = ẑ` instead.
        (
            SketchPlane::world_zx(),
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ),
    ] {
        let d = p.def.expect("a world plane states itself");
        assert_eq!(f(&d), want);
        assert_eq!(d.ref_dir().map(|r| r.to_f64()), u);
        assert_eq!(p.x_axis().as_array(), u, "the axis follows the definition");
    }
    // A tilted normal the caller wrote: exact coefficients, though its axes never can be.
    let tilt =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    assert_eq!(f(&tilt.def.unwrap()), [1.0, 1.0, 1.0, 0.0]);
    assert!(
        tilt.exact().is_none(),
        "★ the axes still have no exact form — that is what the definition exists to replace"
    );
    // Three written points: the plane is exact and `+u` runs toward `x_point`.
    let tp = SketchPlane::through_points(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([1.0, 2.0, 0.0]),
        Point3::from_array([1.0, 0.0, 3.0]),
    )
    .unwrap();
    let d = tp.def.unwrap();
    assert_eq!(f(&d), [1.0, 0.0, 0.0, -1.0], "the plane x = 1");
    assert_eq!(
        d.ref_dir().map(|r| r.to_f64()),
        [0.0, 2.0, 0.0],
        "x_point − origin"
    );
    assert_eq!(d.origin().map(|r| r.to_f64()), [1.0, 0.0, 0.0]);
    // Moving the sketch origin keeps the plane and moves only `(0, 0)`.
    let moved = tp.with_origin(Point3::from_array([1.0, 5.0, 5.0]));
    let m = moved.def.unwrap();
    assert_eq!(f(&m), f(&d), "same plane");
    assert_eq!(m.origin().map(|r| r.to_f64()), [1.0, 5.0, 5.0]);
    // ★★★★★ **The invariant that used to tie the two halves together — the origin is *on*
    // the plane — is structural now: the origin IS `points[0]`, so there are no halves to
    // disagree (the failure this guards against cost a boolean 0.04 of volume once). The
    // loop keeps the check as a derivation audit: substituting the origin into the *derived*
    // name must still give zero, or the derivation itself is wrong.
    for p in [
        SketchPlane::world_xy(),
        SketchPlane::world_yz(),
        SketchPlane::world_zx(),
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
        SketchPlane::world_zx().with_origin(Point3::from_array([10.0, 20.0, 5.0])),
        SketchPlane::from_origin_normal(
            Point3::from_array([0.0, 1.3, 0.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
        )
        .unwrap(),
        tilt,
        tp,
        moved,
    ] {
        let d = p.def.expect("stated");
        let pts = d.points();
        let c = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2])
            .expect("a definition names a plane")
            .narrow()
            .copied()
            .expect("these fixtures are narrow");
        let mut s = c[3];
        for (ck, ok) in c.iter().zip(d.origin()) {
            s = s.checked_add(ck.checked_mul(ok).unwrap()).unwrap();
        }
        assert_eq!(
            s,
            nacre_scalar::Rat::from_int(0),
            "the origin must lie on the plane its points name: {:?} vs {:?}",
            c.map(|r| r.to_f64()),
            d.origin().map(|r| r.to_f64())
        );
    }
    // ★ And the axes-only route records the axes' decimal truth (S6a — this used to assert
    // `def.is_none()`, "honest about having no definition"; the honest statement now is the
    // definition itself, `[o, o + x, o + y]`).
    let axes = SketchPlane::from_axes(
        Point3::from_array([1.0, 2.0, 3.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    );
    let r = |v: [f64; 3]| v.map(|x| nacre_scalar::Rat::from_decimal(x).unwrap());
    assert_eq!(
        axes.def.expect("axes state their truth").points(),
        [r([1.0, 2.0, 3.0]), r([2.0, 2.0, 3.0]), r([1.0, 3.0, 3.0])],
        "o, o + x, o + y"
    );
    // A degenerate pair still names nothing.
    assert!(
        SketchPlane::from_axes(
            Point3::origin(),
            Vector3::from_array([1.0, 0.0, 0.0]),
            Vector3::from_array([2.0, 0.0, 0.0]),
        )
        .def
        .is_none(),
        "parallel axes name no plane"
    );
}

/// ★★★★★ **A prism raised on a named tilted plane states every one of its faces.**
///
/// In world coordinates that plane's axes are irrational, so the whole prism used to drop to
/// f64 and record nothing. Two things fixed it: the base cap **is** the plane the caller
/// named, so it states itself in the world; and the walls and far cap are built **inside that
/// plane's frame**, where the axes are `x̂`/`ŷ` and the profile's own decimals are the truth.
///
/// ★★★ **The base cap stays in the world on purpose.** Writing it as `[0,0,1,0]` in this
/// prism's frame would be a second exact description of one plane under a different
/// `SurfaceKey` — the duplication this work exists to remove. Stated in the world it is
/// `Constructed`, its judgment stays exact, and two extrudes share it whatever frames they chose.
#[test]
fn a_prism_on_a_named_tilted_plane_states_all_of_its_faces() {
    let plane =
        SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
            .unwrap();
    assert!(plane.exact().is_none(), "the axes have no exact form");
    let mut m = Model::new();
    let __f102 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f102,
            profile: square(),
            dist: 1.0,
        },
    )
    .expect("extrude on a tilted plane") else {
        unreachable!()
    };
    let coeffs = |f: Handle<Face>, m: &Model| {
        m.surface_name
            .get(&m.face(f).surface)
            .and_then(|n| n.narrow())
            .map(|c| c.map(|r| r.to_f64()))
    };
    assert_eq!(
        coeffs(faces[0], &m),
        Some([1.0, 1.0, 1.0, 0.0]),
        "★ the base cap is the caller's plane, in the world"
    );
    assert!(
        matches!(
            m.surface(m.face(faces[0]).surface),
            nacre_topo::Surface::Plane { motion: None, .. }
        ),
        "★ and it carries no motion, so its judgment stays exact"
    );
    assert_eq!(
        coeffs(faces[1], &m),
        Some([0.0, 0.0, 1.0, -1.0]),
        "★ the far cap is `w = dist` in the frame"
    );
    // ★ Every face now states itself — that is the whole measurement.
    for &f in &faces {
        assert!(coeffs(f, &m).is_some(), "a face with no exact plane");
    }

    // ★★★★★ **Two extrudes on one named plane put their far caps on one handle** — by name,
    // at construction, with no f64 comparison. That is what the frame buys over the f64 path,
    // where the two would agree only if their rounded coefficients happened to.
    let cap_of = |m: &mut Model, d: f64| -> Handle<Surface> {
        let __g201 = datum_frame(m, plane);
        let OpOutput::Extrude { faces, .. } = apply(
            m,
            &Operation::Extrude {
                frame: __g201,
                profile: square(),
                dist: d,
            },
        )
        .expect("extrude") else {
            unreachable!()
        };
        m.face(faces[1]).surface
    };
    let a = cap_of(&mut m, 2.5);
    let b = cap_of(&mut m, 2.5);
    assert_eq!(a, b, "one height on one plane is one plane");
    assert_ne!(a, cap_of(&mut m, 2.6), "and a different height is not");
}

/// ★★★★★ **The end-to-end lock: the f64-fallback chain is cut on the `n·n`-overflow
/// population** — the measured 1.6%, and the one today's producers actually reach.
///
/// A prism raised on a fully tilted, exactly-orthonormal decimal frame from a 16-digit
/// profile has walls whose names run to ~110 bits: **narrow names whose squared lengths
/// (~2^220) overflow `i128`**, so `plane_frame_default`/`plane_frame_named` hard-declined
/// them and a pad on such a wall fell to `Swept::along` — every new face point-less, the
/// chain that kept the `Inexact` state alive. Now the frame realizes through the wide
/// road and **every face of the result states its exact points**.
///
/// ★ Probed while building this fixture: a sketch→extrude wall's canonical name caps out
/// around ~115 bits (the profile's decimal window bounds the products), so **`Wide`-named
/// faces do not arise from today's construction route at all** — the `Wide` opening is
/// locked at the basis level (`a_wide_plane_hosts_a_canonical_frame`) and the end-to-end
/// chain here on the population that exists. The fixture qualifies itself rather than
/// assuming its population.
#[test]
fn a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    // A fully tilted, exactly-orthonormal decimal frame: u·u = v·v = 1 and u·v = 0 hold in
    // the lifted rationals, and the normal u×v = (0.64, −0.48, 0.6) is not axis-aligned —
    // so the walls' names mix all three coordinates.
    let plane = crate::ops::SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let __f101 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f101,
            profile: Profile2d::polygon(vec![
                p(0.1111111111111111, 0.1234567890123456),
                p(4.123456789012345, 0.2345678901234567),
                p(3.9876543210987654, 3.1234567890123459),
                p(0.2222222222222222, 2.765432109876543),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .expect("extrude on a Pythagorean frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // Fixture qualification: a wall whose name is narrow but whose squared lengths are not
    // — the exact population the narrow frame derivation hard-declines.
    let shell = m.solid(solid).outer;
    let wall = *m
        .shell(shell)
        .faces
        .iter()
        .find(|&&f| {
            let s = m.face(f).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
        })
        .expect("an nn-overflow wall — retune the fixture constants if this fails");
    let profile = centred_on(&m, wall, 0.3);
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 0.4,
        },
    )
    .expect("S4: a pad on a wide-named wall must build") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // The chain is cut: every face of the result states its exact points.
    let mut missing = 0;
    let mut total = 0;
    for &f in &m.shell(m.solid(solid).outer).faces {
        total += 1;
        if !matches!(
            m.surface(m.face(f).surface),
            nacre_topo::Surface::Plane { .. }
        ) {
            missing += 1;
        }
    }
    assert_eq!(
        missing, 0,
        "{missing} of {total} faces carry no exact points — the f64 fallback fired"
    );
}

/// ★★★★★ **An axes-only tilted frame extrudes on the exact road.**
///
/// The sibling of `a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road`, for the
/// population that one cannot reach: a `from_axes` plane whose axes never lift to exact
/// orthonormal rationals (a rotated frame — the caller who only has axes). With no definition
/// the whole prism would fall silently to f64 — every surface point-less, undemotable under
/// any later motion. The axes' decimal truth is the definition, the canonical name is
/// `Wide` (asserted — the full-width crosses exceed `i128`), `WideFrame::named_of` realizes
/// the frame, and **every face of the prism records its exact points**.
#[test]
fn a_prism_on_an_axes_only_tilted_frame_takes_the_exact_road() {
    let mut m = Model::new();
    let plane = crate::ops::SketchPlane::from_axes(
        Point3::from_array([0.2547863291057384, -0.5123456789012345, 1.5432109876543211]),
        Vector3::from_array([0.7123456789012345, 0.5876543210987654, 0.4098765432101234]),
        Vector3::from_array([-0.5876543210987654, 0.7123456789012345, 0.1234567890123456]),
    );
    // Fixture qualification: no world lift (the frame road is the only exact road), the
    // definition exists, and its name is genuinely Wide.
    assert!(
        plane.exact().is_none(),
        "the axes must not lift orthonormal"
    );
    let d = plane.def.expect("S6a: axes state their decimal truth");
    let pts = d.points();
    assert!(
        nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2])
            .expect("a plane")
            .narrow()
            .is_none(),
        "the fixture was chosen to have a Wide name — retune the axes if this fails"
    );
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let __f100 = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f100,
            profile: Profile2d::polygon(vec![
                p(0.1234567890123456, 0.2345678901234567),
                p(2.765432109876543, 0.3456789012345678),
                p(2.543210987654321, 1.9876543210987654),
                p(0.3456789012345678, 1.8765432109876543),
            ])
            .unwrap(),
            dist: 1.3,
        },
    )
    .expect("S6a: an axes-only tilted extrude must build exactly") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // The end-to-end claim: nothing fell to f64 — every face states its exact points.
    let mut missing = 0;
    let mut total = 0;
    for &f in &m.shell(m.solid(solid).outer).faces {
        total += 1;
        if !matches!(
            m.surface(m.face(f).surface),
            nacre_topo::Surface::Plane { .. }
        ) {
            missing += 1;
        }
    }
    assert_eq!(
        missing, 0,
        "{missing} of {total} faces carry no exact points — the f64 fallback fired"
    );
    // Geometry check: the realized frame is orthonormal, so the prism's volume is the
    // profile's own area times the sweep — independent of the tilt.
    let props = nacre_props::mass_props(&m, solid).unwrap();
    let ring = [
        [0.1234567890123456, 0.2345678901234567],
        [2.765432109876543, 0.3456789012345678],
        [2.543210987654321, 1.9876543210987654],
        [0.3456789012345678, 1.8765432109876543],
    ];
    let mut area2 = 0.0f64;
    for i in 0..4 {
        let (a, b) = (ring[i], ring[(i + 1) % 4]);
        area2 += a[0] * b[1] - b[0] * a[1];
    }
    let want = (area2 / 2.0).abs() * 1.3;
    assert!(
        (props.volume - want).abs() < 1e-9,
        "volume {} vs analytic {want}",
        props.volume
    );
}

/// ★★★ S6b: **what the f64 fallback used to build silently is a named reject now.** A plane
/// with no exact statement — axes outside the decimal window — and a sweep distance the
/// window cannot hold each get their own name at the operation's door. The prisms these
/// used to build recorded no exact points and could not survive a motion; the reject is the
/// honest form of the same fact.
#[test]
fn a_prism_the_exact_arithmetic_cannot_state_is_refused_by_name() {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let square =
        || Profile2d::polygon(vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]).unwrap();
    // ① Axes outside the decimal window: no def, no world lift — no exact form at all.
    let far = crate::ops::SketchPlane::from_axes(
        Point3::origin(),
        Vector3::from_array([1e300, 0.0, 0.0]),
        Vector3::from_array([0.0, 1e300, 0.0]),
    );
    assert!(
        far.exact().is_none() && far.def.is_none(),
        "fixture: no exact statement"
    );
    // A plane with no exact statement cannot even be stated as a datum, which is where the
    // rejection now lands — one step earlier than it used to, and by the same name.
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::DatumPlane {
                def: DatumDef::Stated(far)
            },
        ),
        Err(OpError::PlaneWithoutExactForm)
    );
    // ② A sweep distance outside the window — the profile-coordinate rule's sibling.
    assert_eq!(
        apply(
            &mut Model::new(),
            &Operation::Extrude {
                frame: SketchFrame::world(&Model::new(), Axis::Z),
                profile: square(),
                dist: 1e300,
            },
        ),
        Err(OpError::DistOutsideDecimalWindow)
    );
}

fn centred_on(m: &Model, face: Handle<Face>, half: f64) -> Profile2d {
    let sp = crate::ops::face_plane(m, face).expect("planar");
    let d = nacre_props::face_props(m, face).unwrap().centroid - sp.origin;
    let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    Profile2d::polygon(vec![
        p(cu - half, cv - half),
        p(cu + half, cv - half),
        p(cu + half, cv + half),
        p(cu - half, cv + half),
    ])
    .unwrap()
}
