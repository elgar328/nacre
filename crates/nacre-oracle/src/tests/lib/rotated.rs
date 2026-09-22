//! Rotated booleans against OCCT.

use super::*;

// A boolean commutes with a rigid motion, so a rotated boolean is the rotated image of the
// unrotated one — and OCCT, given the rotated STEP, is an *independent* ground truth for its
// volume and area. 3d-i checked rotation-invariance (self-consistency); this cross-checks the
// rotated seam/reconstruction against a mature kernel, catching a same-volume-wrong-topology
// bug that invariance alone could miss.

/// An axis-aligned rotation about the line through `(1,1,0)` by `deg` degrees (a non-90°
/// degree makes cos/sin irrational, so the surfaces record the motion).
fn rot_about(axis: nacre_exact::Axis, deg: i128) -> nacre_exact::Isometry {
    use nacre_exact::{Angle, Isometry, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
    })
}

/// Rotate both operands by every isometry in `isos` (a chain), then assert nacre's boolean
/// matches OCCT on total volume and area. OCCT ground truth is taken *before* the nacre
/// boolean supersedes the rotated inputs. Returns the result solids for further assertions.
fn rotated_boolean_matches_occt(
    m: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    kind: BoolKind,
    isos: &[nacre_exact::Isometry],
) -> Vec<Handle<Solid>> {
    use nacre_ops::{OpOutput, Operation, apply};
    let rot = |m: &mut Model, mut s: Handle<Solid>| -> Handle<Solid> {
        for iso in isos {
            let out = apply(
                m,
                &Operation::Transform {
                    solid: s,
                    isometry: *iso,
                },
            )
            .unwrap();
            let OpOutput::Transform { solid } = out else {
                panic!("expected Transform output");
            };
            s = solid;
            m.rebuild_adjacency();
        }
        s
    };
    let a2 = rot(m, a);
    let b2 = rot(m, b);
    let occt_kind = match kind {
        BoolKind::Fuse => OcctBool::Fuse,
        BoolKind::Cut => OcctBool::Cut,
        BoolKind::Common => OcctBool::Common,
    };
    // OCCT ground truth on the rotated inputs, before nacre supersedes them.
    let occt = occt_boolean_of(m, occt_kind, a2, b2).unwrap();
    let solids = boolean(m, kind, a2, b2).unwrap();
    let vol: f64 = solids
        .iter()
        .map(|&s| mass_props(m, s).unwrap().volume)
        .sum();
    let area: f64 = solids.iter().map(|&s| mass_props(m, s).unwrap().area).sum();
    assert!(
        approx(vol, occt.volume),
        "volume {vol} vs occt {}",
        occt.volume
    );
    assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
    solids
}

/// A fully-tilted (Z then X rotation, every face normal irrational) corner-overlap `Cut`
/// matches OCCT — the strongest cross-check of the rotated arrangement/seam.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_overlap_cut_matches_occt() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    let isos = [rot_about(Axis::Z, 30), rot_about(Axis::X, 30)];
    rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Cut, &isos);
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_overlap_fuse_matches_occt() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Fuse, &[rot_about(Axis::Z, 30)]);
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_overlap_common_matches_occt() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
    rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Common, &[rot_about(Axis::Z, 30)]);
}

/// The sever — a bar cut clean through a cube into two solids (OCCT COMPOUND) — under a full
/// tilt. This is where a rotated outwardness bug would hide: both severed
/// pieces must read outward, so OCCT's total volume/area confirms two material solids, not one
/// with the other misjudged as a cavity.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_sever_cut_matches_occt() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let bar = m.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([2.0, 2.0, 4.0]),
    );
    let isos = [rot_about(Axis::Z, 30), rot_about(Axis::X, 30)];
    let solids = rotated_boolean_matches_occt(&mut m, bar, cube, BoolKind::Cut, &isos);
    assert_eq!(solids.len(), 2, "rotated sever yields two solids");
}

/// A rotated `Cut` of a strictly-contained box leaves a cavity (OCCT BREP_WITH_VOIDS): the
/// volume nets the void, and nacre records exactly one cavity.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_containment_cut_makes_cavity_matches_occt() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let solids =
        rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Cut, &[rot_about(Axis::Z, 30)]);
    assert_eq!(solids.len(), 1, "one solid");
    assert_eq!(
        m.solid(solids[0]).cavities.len(),
        1,
        "the contained box is a cavity"
    );
}

/// A **re-rotated** solid is still a valid b-rep: rotate a
/// cuboid 30° about Z, then 45° about X, so its vertices carry a two-node rotation
/// chain. A rigid re-rotation leaves volume/area invariant, so OCCT must agree with
/// nacre's `mass_props` — confirming the re-rotation multi-pass clone (chained
/// forest, base=root) produced a well-formed solid, not just an invariant-preserving
/// vertex shuffle. Single-solid export avoids summing the superseded intermediates.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rerotated_solid_props_match_occt() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    use nacre_ops::{Operation, apply};
    let rot = |axis, deg: i128| {
        Isometry::rotation(Rotation {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    };
    let mut m = Model::new();
    let c = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let step = |m: &mut Model, s, iso| {
        let out = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: iso,
            },
        )
        .unwrap();
        let nacre_ops::OpOutput::Transform { solid } = out else {
            panic!("expected Transform output");
        };
        solid
    };
    let c1 = step(&mut m, c, rot(Axis::Z, 30));
    let c2 = step(&mut m, c1, rot(Axis::X, 45));

    let occt = occt_props(&nacre_step::to_step_solid(&m, c2).unwrap()).unwrap();
    let nacre = mass_props(&m, c2).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "volume {} vs occt {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs occt {}",
        nacre.area,
        occt.area
    );
}

/// A **corner-flush** `Common` scored against OCCT: an L-prism and a box that both start at
/// the origin, so three of their face planes coincide (`z = 0`, `x = 0`, `y = 0`) and every
/// vertex of the shared corner sits exactly on the other solid's planes. The ops-side lock
/// `a_corner_flush_common_keeps_the_non_convex_overlap` pins the hand-derived volume 1.0 and
/// area 7.0, and this scores the same shape against an independent kernel.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn corner_flush_common_matches_occt() {
    use nacre_ops::{BoolKind, Operation, Profile2d, apply};
    let mut m = Model::new();
    let __w9 = SketchFrame::world(&m, Axis::Z);
    apply(
        &mut m,
        &Operation::Extrude {
            frame: __w9,
            profile: Profile2d::polygon(
                [
                    [0.0, 0.0],
                    [2.0, 0.0],
                    [2.0, 1.0],
                    [1.0, 1.0],
                    [1.0, 2.0],
                    [0.0, 2.0],
                ]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
            )
            .unwrap(),
            dist: 1.0,
        },
    )
    .unwrap();
    let l = *m.live_solids().first().unwrap();
    let b = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.5, 1.5, 0.5]),
    );
    let occt = occt_boolean_of(&m, OcctBool::Common, l, b).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, l, b).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "corner-flush common volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "corner-flush common area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// nacre's own `Common` result diffed against OCCT: build two overlapping
/// cubes, ask OCCT for the intersection volume, and compare it to
/// `mass_props` of the solid nacre's half-space enumeration produced.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn common_result_volume_matches_occt() {
    use nacre_ops::BoolKind;
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    // OCCT ground truth from the inputs (before nacre supersedes them).
    let occt = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
    let r = boolean_one(&mut m, BoolKind::Common, a, b).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// nacre's `Fuse` and `Cut` results diffed against OCCT for two overlapping
/// cubes (M5-c4 face clipping).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn fuse_cut_result_volume_matches_occt() {
    use nacre_ops::BoolKind;
    let boxes = || {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        (m, a, b)
    };
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        let occt_kind = match kind {
            BoolKind::Fuse => OcctBool::Fuse,
            BoolKind::Cut => OcctBool::Cut,
            BoolKind::Common => unreachable!(),
        };
        let (mut m, a, b) = boxes();
        let occt = occt_boolean_of(&m, occt_kind, a, b).unwrap();
        let r = boolean_one(&mut m, kind, a, b).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{kind:?}: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn containment_cut_cavity_matches_occt() {
    use nacre_ops::BoolKind;
    // A = [0,3]³ with B = [1,2]³ strictly inside ⇒ A − B is a hollow solid
    // (an internal void). OCCT diffs the two cavity-free inputs; nacre builds
    // the cavity independently, so the volume agreement is non-self-referential.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    // Capture OCCT's answer before the boolean supersedes the inputs.
    let occt = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "cavity cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn hollow_solid_step_volume_matches_occt() {
    // A = [0,4]³ (64) with a concentric B = [1,3]³ (8) void ⇒ material 56.
    // OCCT reads nacre's BREP_WITH_VOIDS export and computes the material
    // volume — the true gate on the exported void orientation. A flipped
    // void would read as 72 (= V_A + V_B), which the face-count round-trip
    // in nacre-step cannot catch.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
    let b_outer = m.solid(b).outer;
    let void = m.reversed_shell(b_outer);
    let a_outer = m.solid(a).outer;
    let hollow = m.push_solid(Solid {
        outer: a_outer,
        cavities: vec![void],
    });
    m.restore_live(vec![hollow]);

    let occt = occt_props_of(&m).unwrap();
    let nacre = mass_props(&m, hollow).unwrap();
    assert!(
        (nacre.volume - 56.0).abs() < 1e-9,
        "nacre volume {}",
        nacre.volume
    );
    assert!(
        approx(nacre.volume, occt.volume),
        "hollow volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn nonconvex_containment_cut_matches_occt() {
    use nacre_math::Point2;
    use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
    // A concave L-prism with a box strictly inside its bottom bar — a
    // non-convex containment cut ⇒ the L with an internal box void. OCCT reads
    // both (concave) inputs and cuts them independently of nacre's
    // point-in-polyhedron classification.
    let mut m = Model::new();
    let l_profile = Profile2d::polygon(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]
        .iter()
        .map(|&p| Point2::from_array(p))
        .collect(),
    )
    .unwrap();
    let __w8 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: l, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w8,
            profile: l_profile,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output");
    };
    let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "non-convex containment cut: {} vs {}",
        nacre.volume,
        occt.volume
    );
}
