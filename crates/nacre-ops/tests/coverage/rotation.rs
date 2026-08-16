//! Rotated operands — rotation-invariance of the boolean.

#![allow(unused_imports)]
use crate::common::*;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::SketchFrame;
use nacre_ops::{
    BoolError, BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean,
    replay,
};
use nacre_scalar::Axis;
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Orientation, Shell, Solid, Vertex};

#[test]
fn rotated_corner_bite_cut_is_rotation_invariant() {
    let (mut m, l, bx) = l_and_corner_box();
    let l = xf(&mut m, l, rot30());
    let bx = xf(&mut m, bx, rot30());
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_corner_bite_fuse_is_rotation_invariant() {
    let (mut m, l, bx) = l_and_corner_box();
    let l = xf(&mut m, l, rot30());
    let bx = xf(&mut m, bx, rot30());
    let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 + 0.924 - 0.224)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_sever_cut_severs_into_two() {
    let (mut m, l, rod) = l_and_rod();
    let l = xf(&mut m, l, rot30());
    let rod = xf(&mut m, rod, rot30());
    let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
    m.rebuild_adjacency();
    assert_eq!(solids.len(), 2, "rotated sever still yields two solids");
    assert!(nacre_validate::validate(&m).is_empty());
    let vol: f64 = solids
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
}

#[test]
fn rotated_containment_cut_makes_cavity() {
    let (mut m, a, b) = l_and_inner_box();
    let a = xf(&mut m, a, rot30());
    let b = xf(&mut m, b, rot30());
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(
        m.solids.get(r).cavities.len(),
        1,
        "the inner box becomes a cavity"
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.512)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn a_rotated_cut_result_can_be_cut_again() {
    // Rotating a boolean *result* (every vertex `Discovered`) and feeding it back into a boolean
    // used to reject (`ROTATED_UNSUPPORTED`) — `collect_planes` could not build the rotated seam
    // faces' plane witnesses. Now each such plane is witnessed through its provenance (its plane
    // is `R(π)` for an operand plane `π`). R = a cube minus a far-corner octant; rotate R and a
    // fresh slab, then Cut. A boolean commutes with a rigid motion, so the volume matches the
    // unrotated chain (7 − 2 = 5) and the result stays valid.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let r = xf(&mut m, r, rot30());
    // A slab that severs at x = 0.5 — no plane coincides with R's (avoids rotated coplanar).
    let c = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, -1.0]),
        Point3::from_array([0.5, 4.0, 4.0]),
    );
    let c = xf(&mut m, c, rot30());
    let r2 = boolean_one(&mut m, BoolKind::Cut, r, c).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r2).unwrap().volume;
    assert!((vol - 5.0).abs() < 1e-9, "volume {vol}");
}

#[test]
fn a_rotated_profile_is_the_same_solid_to_the_boolean() {
    // Rotating a profile's winding cannot change the solid. The cap here carries the
    // dimple's hole, so an inward `n_out` would reverse it — only `validate` and the
    // signed volume would ever say so.
    let (mut m, l) = rotated_l_prism();
    let stub = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
    m.rebuild_adjacency();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
}

#[test]
fn rotated_coplanar_contact_boss_and_pocket_solve() {
    // Rotated coplanar contact *away from a corner* is solved correctly (it was over-rejected by an
    // earlier rotated-coplanar guard, now removed — the Euler-parity self-check only rejects
    // genuinely malformed, odd-Euler results). A tilted boss fuses flush onto its base (z=1 shared);
    // a tilted pocket carves flush into it. Both must succeed with the exact volume and stay valid.
    let tilt = |m: &mut Model, s| xf(m, s, rot_iso(Axis::X, 30));
    // boss fuse — fused volume 1.25.
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let boss = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 1.0]),
        Point3::from_array([0.75, 0.75, 2.0]),
    );
    let (base, boss) = (tilt(&mut m, base), tilt(&mut m, boss));
    let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).expect("tilted boss fuse solves");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "solved boss is valid"
    );
    let v = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((v - 1.25).abs() < 1e-9, "solved boss fuse is correct: {v}");
    // pocket cut — carved volume 0.875.
    let mut m = Model::new();
    let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let prism = m.add_cuboid(
        Point3::from_array([0.25, 0.25, 0.5]),
        Point3::from_array([0.75, 0.75, 1.0]),
    );
    let (base, prism) = (tilt(&mut m, base), tilt(&mut m, prism));
    let r = boolean_one(&mut m, BoolKind::Cut, base, prism).expect("tilted pocket cut solves");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "solved pocket is valid"
    );
    let v = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!(
        (v - 0.875).abs() < 1e-9,
        "solved pocket cut is correct: {v}"
    );
}

#[test]
fn rotated_result_coplanar_reuse_under_a_general_rotation() {
    // The intersection of two separately-solved capabilities: reusing a boolean *result* (every
    // vertex `Discovered`) *and* coplanar contact, taken through a **compound** rotation so the
    // shared plane's normal is fully general (no zero component) — the axis-aligned or single-axis
    // tilt both leave one component zero, which under-tests the coplanar decision.
    //
    // Two unit boxes sit side by side sharing the vertical plane `x = 1`; the right one is a
    // `Cut` result (a corner nicked away from `x = 1`, so its seating face survives). Both are then
    // rotated `Z 35°` then `X 40°`: the shared normal `(1,0,0)` maps to
    // `(cos35, sin35·cos40, sin35·sin40) ≈ (0.819, 0.439, 0.369)`, all nonzero. A `Fuse` must weld
    // them flush — volume `1 + vol(right)` — in either operand order, and stay valid.
    let tilt = |m: &mut Model, s| {
        let s = xf(m, s, rot_iso(Axis::Z, 35));
        xf(m, s, rot_iso(Axis::X, 40))
    };
    for swap in [false, true] {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b_raw = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        let nick = m.add_cuboid(
            Point3::from_array([1.6, 0.6, 0.6]),
            Point3::from_array([2.5, 1.5, 1.5]),
        );
        let b = boolean_one(&mut m, BoolKind::Cut, b_raw, nick).unwrap(); // vol 1 − 0.064 = 0.936
        m.rebuild_adjacency();
        let bvol = nacre_props::mass_props(&m, b).unwrap().volume;
        let (a, b) = (tilt(&mut m, a), tilt(&mut m, b));
        let (x, y) = if swap { (b, a) } else { (a, b) };
        let r = boolean_one(&mut m, BoolKind::Fuse, x, y).unwrap_or_else(|e| {
            panic!("rotated result + coplanar contact (swap={swap}) rejected: {e:?}")
        });
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "swap={swap}: {vs:?}");
        let v = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!(
            (v - (1.0 + bvol)).abs() < 1e-9,
            "swap={swap}: fused volume {v} vs {}",
            1.0 + bvol
        );
    }
}

/// The near misses around the four-plane concurrency.
///
/// It needs **four** coincidences at once: 45°, a bar whose half-width equals its
/// pivot-to-bottom offset, a pivot in the plane those two measure from, and the fused block — which
/// is what supplies the plane `x = 0.5` that the bar's bottom corner edge lands in. Break any single
/// one and the concurrency is gone. Each row here was measured, and together they are what says the
/// case was a genuine concurrency rather than a rotation the engine could not handle: 44° and 46°
/// are every bit as irrational as 45°.
///
/// **The concurrent model itself builds now** (`nacre-oracle`'s `the_four_plane_cut_matches_occt`,
/// where OCCT agrees on volume and centroid), so these rows have changed job: they were the near
/// misses around a reject, and they are now the net that says normalising vertex identity did not
/// break the neighbourhood it was reached through.
#[test]
fn the_near_misses_around_the_four_plane_reject_all_build() {
    use nacre_scalar::Rat;
    let on_plane = Rat::from_int(1);
    for (what, half_z, deg, pivot_z, fuse) in [
        ("44°, one degree short", 0.2, 44, on_plane, true),
        ("46°, one degree past", 0.2, 46, on_plane, true),
        ("30°", 0.2, 30, on_plane, true),
        (
            "90°, the axis-aligned quarter turn",
            0.2,
            90,
            on_plane,
            true,
        ),
        ("a taller-than-wide bar", 0.3, 45, on_plane, true),
        ("a flatter-than-wide bar", 0.15, 45, on_plane, true),
        (
            "pivot lifted 0.05 off the plane",
            0.2,
            45,
            Rat::new(21, 20).unwrap(),
            true,
        ),
        (
            "pivot lifted 0.25 off the plane",
            0.2,
            45,
            Rat::new(5, 4).unwrap(),
            true,
        ),
        (
            "the plain cube, no block fused on",
            0.2,
            45,
            on_plane,
            false,
        ),
    ] {
        let (mut m, cube, bar) = cube_and_spun_bar(half_z, deg, pivot_z);
        let target = if fuse {
            let block = m.add_cuboid(
                Point3::from_array([0.5, 0.0, 1.0]),
                Point3::from_array([1.0, 1.0, 2.0]),
            );
            m.rebuild_adjacency();
            let t = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
            m.rebuild_adjacency();
            t
        } else {
            cube
        };
        let r = boolean_one(&mut m, BoolKind::Cut, target, bar)
            .unwrap_or_else(|e| panic!("{what}: {e:?}"));
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{what}: invalid");
        // The bar passes clean through, so the cut always removes material.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        let whole = if fuse { 1.5 } else { 1.0 };
        assert!(vol > 0.0 && vol < whole, "{what}: volume {vol}");
    }
}

/// ★ **The model that started the four-plane work now builds** — the switch the reject test used to
/// hold open.
///
/// A bar spun 45° whose bottom corner edge lands exactly in the plane the fused block supplies
/// (`x = 0.5`): three planes then share that edge's line, and every plane crossing it makes a vertex
/// with four planes through it. Naming such a vertex is what the substrate could not do.
///
/// **Two independent things are asserted, because this is the first time the engine answers here
/// rather than declining, and "it stopped rejecting" is not the same as "it is right".**
///
/// 1. The result is a **valid** 2-manifold and its volume is the hand figure `1.38` — the target is
///    `1 + 0.5` and the bar removes `0.12`.
/// 2. **It agrees with its non-degenerate neighbours.** Moving one operand coordinate by a single
///    ULP either way destroys the concurrency, and those two models already built before any of
///    this work. A degenerate answer that did not match them would be a discontinuity, and this is
///    the cheapest way to notice one.
#[test]
fn a_tool_edge_lying_in_a_target_plane_builds_and_agrees_with_its_neighbours() {
    use nacre_scalar::Rat;
    let cut_volume = |ulps: i64| -> (f64, bool) {
        let (mut m, cube, bar) = cube_and_spun_bar_ulp(0.2, 45, Rat::from_int(1), ulps);
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        m.rebuild_adjacency();
        let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
        m.rebuild_adjacency();
        let solids = boolean(&mut m, BoolKind::Cut, target, bar)
            .unwrap_or_else(|e| panic!("ulps {ulps}: {e:?}"));
        m.rebuild_adjacency();
        let vol = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        (vol, nacre_validate::validate(&m).is_empty())
    };

    let (vol, valid) = cut_volume(0);
    assert!(valid, "the four-plane result must be a valid 2-manifold");
    assert!((vol - 1.38).abs() < 1e-9, "volume {vol}, expected 1.38");

    for ulps in [1, -1] {
        let (near, valid) = cut_volume(ulps);
        assert!(valid, "the {ulps:+} ULP neighbour must stay valid");
        assert!(
            (near - vol).abs() < 1e-9,
            "the degenerate answer {vol} does not agree with its {ulps:+} ULP neighbour {near}"
        );
    }
}

/// **A rotation history longer than the judging budget is rejected by its cause, not a symptom.**
///
/// The error a rotated definition carries grows about one bit per turn, so the kernel sizes its
/// judging precision from the model — and past `JUDGE_PREC_CAP` it stops and says so. Measured
/// cost is what puts the cap where it is: this same `Cut` takes 1.2s at 400 turns, 16s at 1600 and
/// 90s at 3200 (all with the volume still exactly right), so the cap sits past any real model and
/// still inside the minutes.
///
/// The reject must arrive *before* the work: the precision is known once the plane tables exist,
/// which is why this test is fast while a boolean one turn under the cap is not.
#[test]
fn a_rotation_history_past_the_budget_is_rejected_by_name() {
    let (mut m, a, b) = two_boxes();
    let mut a = a;
    // ~1 bit per turn on top of the ~185 the coincidence limit already asks for. One operand is
    // enough: the precision is sized from the whole table, so the deepest history in it wins.
    for _ in 0..4200 {
        a = xf(&mut m, a, rot_iso(Axis::Z, 37));
    }
    match boolean(&mut m, BoolKind::Cut, a, b) {
        Err(BoolError::Rejected {
            reason: nacre_ops::RejectReason::PrecisionBudget { needed, cap },
            ..
        }) => {
            assert!(needed > cap, "needed {needed} bits, cap {cap}");
            assert!(needed > 4200, "one bit per turn: needed {needed}");
        }
        other => panic!("expected a named precision-budget reject, got {other:?}"),
    }
}

/// **A 40-fin ring, whose mirror-pair fins are exactly symmetric, builds.**
///
/// It did not, once the trig realization became correctly rounded: `cos 288° == cos 72°` exactly,
/// so mirror-pair fins landed on genuinely mirrored coordinates where `libm` had broken every such
/// pair by an ulp, and the arrangement met the symmetry for the first time. What it met there was
/// **not a symmetry problem** — it was a plane described two ways.
///
/// A `Constructed` plane carries its stored coefficients *and* a witness triangle, and the two need
/// not describe the same plane: `d = −(raw·origin)` is an `f64` product, and `3.5 × 0.2` lands on
/// `0.7000000000000001`. Predicates chose between the two descriptions by *the question* (does any
/// plane in it rotate?), so one plane sat in two places depending on what was asked, and comparisons
/// composed across the two were not even transitive. `loop_winding` then read its turn at a node
/// that was not extreme and returned a **confident wrong winding**; the `RingOrientation` reject
/// two layers later was the symptom.
///
/// The route now requires the descriptions to agree, so this builds. **It is a regression test for
/// the transitivity, not for the fins** — any model that puts two descriptions of one plane into a
/// single ordering would fail the same way.
#[test]
fn a_forty_fin_ring_builds_despite_exact_mirror_symmetry() {
    use nacre_scalar::{Angle, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let ring = |x0: f64, y0: f64, x1: f64, y1: f64| {
        vec![
            Point2::from_array([x0, y0]),
            Point2::from_array([x1, y0]),
            Point2::from_array([x1, y1]),
            Point2::from_array([x0, y1]),
        ]
    };
    let extrude = |m: &mut Model, prof: Vec<Point2>, dist: f64| {
        let out = apply(
            m,
            &Operation::Extrude {
                frame: SketchFrame::world(m, Axis::Z),
                profile: Profile2d::polygon(prof).unwrap(),
                dist,
            },
        )
        .expect("extrude");
        m.rebuild_adjacency();
        match out {
            OpOutput::Extrude { solid, .. } => solid,
            o => panic!("{o:?}"),
        }
    };
    let mut acc = extrude(&mut m, ring(-3.0, -3.0, 3.0, 3.0), 2.0);
    for i in 0..40i128 {
        let fin = extrude(&mut m, ring(2.0, -0.4, 8.0, 0.4), 1.0);
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: fin,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::new(9 * i, 1).unwrap()).expect("angle"),
                }),
            },
        )
        .expect("rotate");
        m.rebuild_adjacency();
        let OpOutput::Transform { solid: fin } = out else {
            unreachable!()
        };
        acc = boolean(&mut m, BoolKind::Fuse, acc, fin)
            .unwrap_or_else(|e| panic!("fin {i} ({}deg) declined: {e:?}", 9 * i))[0];
        m.rebuild_adjacency();
    }
}
