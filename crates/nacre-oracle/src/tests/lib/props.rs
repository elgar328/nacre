//! Mass properties against OCCT: volume, area, centroid, moved and mirrored solids, cylinders, pads
//! and pockets.

use super::*;

/// Combined relative-or-absolute float comparison, sized to DRAWEXE's output
/// precision — it prints ~6 significant figures, so an exact value can land
/// ~5e-7 relative away (e.g. 20π → `62.8319`). A 1e-4 relative band clears
/// that rounding noise by 100× while still catching any real geometry error
/// (those miss by percents, not parts-per-thousand).
/// Faces of one live solid, outer shell and cavities — the figure `OcctProps::faces` is
/// compared against.
fn face_count(model: &Model, s: Handle<Solid>) -> usize {
    let sol = model.solid(s);
    std::iter::once(sol.outer)
        .chain(sol.cavities.iter().copied())
        .map(|sh| model.shell(sh).faces.len())
        .sum()
}

#[test]
fn parse_reads_all_fields_order_independent() {
    let stdout = "\
faces 6
volume 24
area 52
bbox_max 2 3 4
bbox_min 0 0 0
centroid 1 1.5 2
";
    let p = OcctProps::parse(stdout).unwrap();
    assert_eq!(p.volume, 24.0);
    assert_eq!(p.area, 52.0);
    assert_eq!(p.faces, 6);
    assert_eq!(p.bbox_min, [0.0, 0.0, 0.0]);
    assert_eq!(p.bbox_max, [2.0, 3.0, 4.0]);
    assert_eq!(p.centroid, [1.0, 1.5, 2.0]);
}

#[test]
fn parse_reports_missing_field() {
    let stdout = "volume 24\narea 52\nfaces 6\nbbox_min 0 0 0\n";
    match OcctProps::parse(stdout) {
        Err(OracleError::Parse(msg)) => assert!(msg.contains("bbox_max")),
        other => panic!("expected Parse error, got {other:?}"),
    }
}

#[test]
fn parse_reports_non_numeric() {
    let stdout = "volume oops\narea 52\nfaces 6\nbbox_min 0 0 0\nbbox_max 2 3 4\ncentroid 1 1 1\n";
    assert!(matches!(
        OcctProps::parse(stdout),
        Err(OracleError::Parse(_))
    ));
}

// End-to-end oracle tests: they shell out to DRAWEXE, so they are #[ignore]d
// (the pre-commit hook still compiles + clippy + fmt them). Run on a machine
// with `brew install opencascade`:  cargo test -p nacre-oracle -- --ignored
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cuboid_matches_occt() {
    let mut model = Model::new();
    model.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let p = occt_props_of(&model).unwrap();
    assert!(approx(p.volume, 24.0), "volume {}", p.volume); // 2·3·4
    assert!(approx(p.area, 52.0), "area {}", p.area); // 2(6+8+12)
    assert_eq!(p.faces, 6);
}

/// **The centroid's judge.** nacre derives the centroid from a cone decomposition;
/// OCCT's `vprops` computes it independently. Both shapes are deliberately
/// asymmetric — a symmetric solid lands its centroid in the middle whatever the
/// arithmetic does, so it would score nothing.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn centroid_matches_occt() {
    use nacre_math::Point2;
    use nacre_ops::{Operation, Profile2d, apply};

    // (a) An L-prism: reflex outline, centroid off both the bbox centre and the
    //     vertex average.
    let pts = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let mut model = Model::new();
    let op = Operation::Extrude {
        frame: SketchFrame::world(&model, Axis::Z),
        profile: Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap(),
        dist: 3.0,
    };
    apply(&mut model, &op).unwrap();
    model.rebuild_adjacency();
    let live = model.live_solids()[0];
    let p = occt_props_of(&model).unwrap();
    let c = nacre_props::centroid(&model, live).unwrap();
    for i in 0..3 {
        assert!(
            approx(c[i], p.centroid[i]),
            "L axis {i}: {} vs {}",
            c[i],
            p.centroid[i]
        );
    }

    // (b) A box with an off-centre void. The cavity shell carries the opposite
    //     sign; getting that wrong is invisible on a centred void.
    let mut model = Model::new();
    let outer = model.add_cuboid(Point3::origin(), Point3::from_array([10.0; 3]));
    let inner = model.add_cuboid(
        Point3::from_array([1.0; 3]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let r = boolean(&mut model, BoolKind::Cut, outer, inner).unwrap();
    model.rebuild_adjacency();
    assert_eq!(r.len(), 1);
    let p = occt_props_of(&model).unwrap();
    let c = nacre_props::centroid(&model, r[0]).unwrap();
    for i in 0..3 {
        assert!(
            approx(c[i], p.centroid[i]),
            "hollow axis {i}: {} vs {}",
            c[i],
            p.centroid[i]
        );
    }
}

/// **The bounding box's judge, on the shape that makes it hard.** A cylinder's
/// barrel bulges past its seam vertices, so a vertex hull would come out too
/// small — and OCCT knows the true extent.
///
/// The comparison is one-sided on purpose: DRAWEXE's `bounding` returns a
/// *conservative* box (measured: ~1e-7 of slack on a unit-scale part), so the
/// requirement is that nacre's box sits inside OCCT's and is not meaningfully
/// smaller — an equality assert would fail on OCCT's padding, not on a bug.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cylinder_bounds_match_occt() {
    let mut model = Model::new();
    let s = model.add_cylinder(
        Point3::origin(),
        Vector3::from_array([0.0, 3.0, 4.0]),
        2.0,
        5.0,
    );
    model.rebuild_adjacency();
    let p = occt_props_of(&model).unwrap();
    let (lo, hi) = nacre_props::bounds(&model, s).unwrap();
    let slack = 1e-5;
    for i in 0..3 {
        assert!(
            lo[i] >= p.bbox_min[i] - slack && lo[i] <= p.bbox_min[i] + slack,
            "min axis {i}: nacre {} vs occt {}",
            lo[i],
            p.bbox_min[i]
        );
        assert!(
            hi[i] <= p.bbox_max[i] + slack && hi[i] >= p.bbox_max[i] - slack,
            "max axis {i}: nacre {} vs occt {}",
            hi[i],
            p.bbox_max[i]
        );
    }
}

/// ★ **The first curved boolean, scored by an independent kernel**. A `[0,2]³`
/// box drilled through by a radius-0.5 bore: volume catches a bore that never opened, area
/// catches a missing or inverted wall, and the face count catches a cap whose inner loop was
/// dropped. Our own analytic figure (`8 − πr²h`) is not evidence about the *topology* — OCCT
/// reading the STEP back is.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn through_hole_cut_matches_occt() {
    let mut model = Model::new();
    let a = model.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = model.add_cylinder(
        Point3::from_array([1.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    model.rebuild_adjacency();
    let s = boolean_one(&mut model, BoolKind::Cut, a, b).unwrap();
    let ours = mass_props(&model, s).unwrap();
    let p = occt_props_of(&model).unwrap();
    // `approx`, not a hand-picked epsilon: the harness reports OCCT's figures to six
    // significant digits, so a tighter comparison measures the *printout* rather than the
    // two kernels (measured — 1e-6 absolute failed on 6.4292 vs 6.4292036732051026).
    assert!(
        approx(ours.volume, p.volume),
        "volume: nacre {} vs occt {}",
        ours.volume,
        p.volume
    );
    assert!(
        approx(ours.area, p.area),
        "area: nacre {} vs occt {}",
        ours.area,
        p.area
    );
    assert_eq!(p.faces, 7, "4 walls + 2 drilled caps + the bore wall");
    // And the analytic figure, so a *pair* of kernels agreeing on a wrong number would still
    // have to agree with arithmetic.
    assert!((ours.volume - (8.0 - PI * 0.25 * 2.0)).abs() < 1e-9);
}

/// **The tangency fixtures, scored by a second kernel**.
///
/// Two booleans in `nacre-ops` produce a solid whose *face* is pinched at one point — a boss
/// edge exactly tangent to a bore's rim, and a boss's base circle exactly tangent to the
/// plate's top edge. `validate` is clean and the volume is exact; the question
/// is whether the **solid** is valid.
///
/// The geometry answers that — the link of the boundary at the touch is a single
/// circle, so the surface is a 2-manifold there and only the *face* is pinched. **This is the
/// independent confirmation**: a second kernel is given the same two operands and asked for
/// the same union. If OCCT returns a body of the same volume, a commercial kernel agrees the
/// shape is a body.
///
/// ★ **Why the face count is recorded and not asserted.** OCCT splits periodic surfaces at
/// its own seams and merges or divides faces across a STEP round trip, so an absolute count
/// says nothing about the pinched face. The **ε-twin** (the same model with the tangency
/// broken by 0.01) is measured beside it so the difference has a control; both numbers are
/// printed rather than asserted.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_segment_tangent_to_a_rim_is_a_body_to_occt() {
    // The plate with a bore, and a boss whose `x = 11` face is exactly tangent to the rim.
    let build = |bx: f64| -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = boolean_one(&mut m, BoolKind::Cut, plate, hole).expect("the bore cuts");
        let boss = m.add_cuboid(
            Point3::from_array([bx, 8.0, 5.0]),
            Point3::from_array([15.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        (m, holed, boss)
    };
    // `11.0` is the tangency; `11.01` breaks it and is the control.
    for (what, bx) in [("tangent", 11.0), ("twin", 11.01)] {
        let (mut m, a, b) = build(bx);
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).expect("occt fuses");
        let s = boolean_one(&mut m, BoolKind::Fuse, a, b).expect("nacre fuses");
        let ours = mass_props(&m, s).unwrap();
        eprintln!(
            "TANOCCT seg-{what}: volume nacre {} occt {} | area nacre {} occt {} | faces nacre {} occt {}",
            ours.volume,
            occt.volume,
            ours.area,
            occt.area,
            face_count(&m, s),
            occt.faces
        );
        assert!(
            approx(ours.volume, occt.volume),
            "{what}: volume nacre {} vs occt {}",
            ours.volume,
            occt.volume
        );
    }
}

/// ★★★★★ **The tangent *wall* — the user's own script, scored against OCCT**. A unit
/// cube and a stud whose axis stands `0.3` from the origin with `r = 0.2`, so the wall
/// `x = 0.5` is exactly `r` away. The kernel used to refuse this outright; it now fuses, and
/// the question this asks is the same one asked of a point tangency: **does a second
/// kernel agree the shape is a body?**
///
/// ★ The twin (`0.35`, the stud pushed **through** the wall) is the control, and it is chosen
/// so the two volumes actually differ: a twin that merely clears the wall keeps the same
/// overlap and the same union, which would make "the volumes agree" a statement about
/// arithmetic that never moved. Face counts are recorded, never asserted — OCCT splits
/// periodic surfaces at its own seams.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn the_tangent_wall_fuses_to_a_body_occt_agrees_with() {
    let build = |ax: f64| -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let stud = m.add_cylinder(
            Point3::from_array([ax, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        (m, cube, stud)
    };
    // `0.3` is the tangency (`0.3 + 0.2 = 0.5`); `0.35` pushes the stud through the wall.
    for (what, ax) in [("tangent", 0.3), ("twin", 0.35)] {
        let (mut m, a, b) = build(ax);
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).expect("occt fuses");
        let s = boolean_one(&mut m, BoolKind::Fuse, a, b).expect("nacre fuses");
        let ours = mass_props(&m, s).unwrap();
        eprintln!(
            "TANOCCT wall-{what}: volume nacre {} occt {} | area nacre {} occt {} | faces nacre {} occt {}",
            ours.volume,
            occt.volume,
            ours.area,
            occt.area,
            face_count(&m, s),
            occt.faces
        );
        assert!(
            approx(ours.volume, occt.volume),
            "{what}: volume nacre {} vs occt {}",
            ours.volume,
            occt.volume
        );
    }
}

/// The second tangency fixture — a turned boss whose base circle touches the plate's top
/// edge at one point. Same question, same discipline as the test above.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_rim_tangent_to_a_plate_top_is_a_body_to_occt() {
    let build = |bz: f64| -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, bz]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        (m, plate, boss)
    };
    // `1.5` puts the circle's top exactly on `z = 2`; `1.49` clears it by 0.01.
    for (what, bz) in [("tangent", 1.5), ("twin", 1.49)] {
        let (mut m, a, b) = build(bz);
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).expect("occt fuses");
        let s = boolean_one(&mut m, BoolKind::Fuse, a, b).expect("nacre fuses");
        let ours = mass_props(&m, s).unwrap();
        eprintln!(
            "TANOCCT rim-{what}: volume nacre {} occt {} | area nacre {} occt {} | faces nacre {} occt {}",
            ours.volume,
            occt.volume,
            ours.area,
            occt.area,
            face_count(&m, s),
            occt.faces
        );
        assert!(
            approx(ours.volume, occt.volume),
            "{what}: volume nacre {} vs occt {}",
            ours.volume,
            occt.volume
        );
    }
}

/// ★ **Two bores in one plate, scored independently**. The second cut's counterpart
/// already carries a cylinder face, which is the case the band pass was rebuilt for — and the
/// case where assuming "a wall's own solid fills its cylinder" puts `keep` on the wrong
/// chamber. Volume catches that directly; OCCT reading the STEP back catches a topology that
/// merely *measures* right.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn two_bores_match_occt() {
    let mut model = Model::new();
    let a = model.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 2.0, 1.0]),
    );
    let bore = |m: &mut Model, x: f64| {
        m.add_cylinder(
            Point3::from_array([x, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        )
    };
    let d1 = bore(&mut model, 1.0);
    let d2 = bore(&mut model, 3.0);
    model.rebuild_adjacency();
    let one = boolean_one(&mut model, BoolKind::Cut, a, d1).unwrap();
    let two = boolean_one(&mut model, BoolKind::Cut, one, d2).unwrap();
    let ours = mass_props(&model, two).unwrap();
    let p = occt_props_of(&model).unwrap();
    assert!(
        approx(ours.volume, p.volume),
        "volume: nacre {} vs occt {}",
        ours.volume,
        p.volume
    );
    assert!(
        approx(ours.area, p.area),
        "area: nacre {} vs occt {}",
        ours.area,
        p.area
    );
    assert!((ours.volume - (8.0 - 2.0 * PI * 0.0625)).abs() < 1e-9);
}

/// A donut prism scored by an independent kernel: the hole has to be a hole to OCCT too.
/// Volume catches a hole that never opened, area catches missing or inverted hole walls, and
/// the face count catches a cap whose inner loop was dropped. It also exercises the STEP
/// writer's inner face bounds.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn swept_hole_matches_occt() {
    use nacre_math::Point2;
    use nacre_ops::{Operation, Profile2d, apply};

    let sq = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    let mut model = Model::new();
    let __w15 = SketchFrame::world(&model, Axis::Z);
    apply(
        &mut model,
        &Operation::Extrude {
            frame: __w15,
            profile: Profile2d::with_holes(sq(0.0, 4.0), vec![sq(1.0, 3.0)]).unwrap(),
            dist: 1.0,
        },
    )
    .unwrap();
    model.rebuild_adjacency();

    let p = occt_props_of(&model).unwrap();
    assert!(approx(p.volume, 12.0), "volume {}", p.volume); // 16 − 4
    // caps 2·(16 − 4) + outer walls 16·1 + hole walls 8·1
    assert!(approx(p.area, 48.0), "area {}", p.area);
    assert_eq!(p.faces, 10);
}

/// A mirrored solid, scored by an independent kernel — the only net for the one mistake a
/// reflection can make. Turning a solid inside out (reflecting coordinates without rewinding
/// the loops) leaves every per-cell invariant intact, so `validate` cannot see it; volume and
/// area are reflection-invariant, so OCCT reading our STEP must report the source values.
/// It checks the STEP export of mirrored geometry at the same time.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn mirrored_solid_matches_occt() {
    use nacre_math::Point2;
    use nacre_ops::{OpOutput, Operation, apply};
    use nacre_scalar::{Axis, Rat};

    // An L-prism, asymmetric so the mirror cannot be a no-op: area 4, height 1.
    let mut model = Model::new();
    let profile = nacre_ops::Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 0.0]),
        Point2::from_array([2.0, 1.0]),
        Point2::from_array([1.0, 1.0]),
        Point2::from_array([1.0, 3.0]),
        Point2::from_array([0.0, 3.0]),
    ])
    .unwrap();
    let __w14 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w14,
            profile,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    model.rebuild_adjacency();
    apply(
        &mut model,
        &Operation::Mirror {
            solid,
            axis: Axis::X,
            offset: Rat::from_int(0),
        },
    )
    .unwrap();
    model.rebuild_adjacency();

    let p = occt_props_of(&model).unwrap();
    assert!(approx(p.volume, 4.0), "volume {}", p.volume);
    assert!(approx(p.area, 18.0), "area {}", p.area); // caps 2·4 + perimeter 10 · height 1
    assert_eq!(p.faces, 8);
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cylinder_matches_occt() {
    let mut model = Model::new();
    model.add_cylinder(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        5.0,
    );
    let p = occt_props_of(&model).unwrap();
    // OCCT computes these analytically from the CYLINDRICAL_SURFACE, so a
    // match confirms our curved STEP really is a cylinder — and, since a
    // positive volume needs a closed outward-oriented solid, this automates
    // the orientation/validity check M3 deferred to manual FreeCAD.
    assert!(approx(p.volume, 20.0 * PI), "volume {}", p.volume); // π·r²·h
    assert!(approx(p.area, 28.0 * PI), "area {}", p.area); // 2πr² + 2πr·h
    assert_eq!(p.faces, 3); // lateral + 2 caps
}

// nacre-vs-OCCT diff: nacre computes volume/area analytically (nacre-props),
// OCCT computes them independently from the same STEP. Agreement cross-checks
// both — the two are wholly separate implementations. #[ignore]d like the
// other DRAWEXE tests.

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn rotated_result_reuse_diff_occt() {
    // Rotate a boolean *result* (all-`Discovered` geometry) and feed it back into a boolean —
    // the capability the provenance plane-witness enables. R = a cube minus a far-corner
    // octant; rotate R and a fresh severing slab by 30° about Z, then Cut. OCCT computes the
    // same Cut of the two rotated STEP solids; the volumes must agree (nacre's exact
    // arrangement vs OCCT's independent kernel).
    use nacre_ops::{OpOutput, Operation, apply};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    });
    let xf = |m: &mut Model, s: Handle<Solid>| -> Handle<Solid> {
        let out = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: iso,
            },
        )
        .unwrap();
        m.rebuild_adjacency();
        match out {
            OpOutput::Transform { solid } => solid,
            _ => unreachable!(),
        }
    };
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let r = xf(&mut m, r);
    let c = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, -1.0]),
        Point3::from_array([0.5, 4.0, 4.0]),
    );
    let c = xf(&mut m, c);
    // OCCT's boolean of the two rotated solids (taken while both are still live).
    let occt = occt_boolean_of(&m, OcctBool::Cut, r, c).unwrap();
    let out = boolean_one(&mut m, BoolKind::Cut, r, c).unwrap();
    m.rebuild_adjacency();
    let nacre = mass_props(&m, out).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cube_props_diff_occt() {
    let mut model = Model::new();
    let solid = model.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "{} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn cylinder_props_diff_occt() {
    let mut model = Model::new();
    let solid = model.add_cylinder(
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        5.0,
    );
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "{} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pad_diff_occt() {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};

    // Unit cube, then a 0.4-square boss of height 0.5 on the top face. The
    // padded solid's top face has a real hole (a FACE_BOUND in STEP); this
    // checks OCCT reads that holed boss as a closed solid and agrees on its
    // volume (1.08) and area (6.8) with nacre's analytic value.
    let sq = |s: f64| {
        Profile2d::polygon(
            [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut model = Model::new();
    let __w13 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w13,
            profile: sq(1.0),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let boss = Profile2d::polygon(
        [[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut model,
        &Operation::PadOnFace {
            face: faces[1],
            profile: boss,
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };

    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "{} vs {}",
        nacre.area,
        occt.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_diff_occt() {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};

    // Unit cube, then a 0.4-square pocket of depth 0.5 in the top face. The
    // inward walls remove material; this checks OCCT reads the holed,
    // concave solid and agrees (volume 0.92, area 6.8) with nacre.
    let sq = |s: f64| {
        Profile2d::polygon(
            [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut model = Model::new();
    let __w12 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w12,
            profile: sq(1.0),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let pocket = Profile2d::polygon(
        [[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut model,
        &Operation::PocketOnFace {
            face: faces[1],
            profile: pocket,
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };

    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "{} vs {}",
        nacre.area,
        occt.area
    );
}

/// A pad whose footprint overhangs one face edge (a boss cantilever). This is the first time
/// the `PadOnFace` pipe produces an overhang — it routes to the overhang Fuse sidecar. OCCT
/// scores the cantilevered boss.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn overhang_pad_matches_occt() {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    let mut model = Model::new();
    let __w11 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w11,
            profile: Profile2d::polygon(
                [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
                    .iter()
                    .map(|&p| nacre_math::Point2::from_array(p))
                    .collect(),
            )
            .unwrap(),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    // Footprint world x in [0.25,0.75], y in [-0.25,0.75] - overhangs the y=0 edge.
    let boss = Profile2d::polygon(
        [[-0.25, -0.25], [0.75, -0.25], [0.75, 0.25], [-0.25, 0.25]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut model,
        &Operation::PadOnFace {
            face: faces[1],
            profile: boss,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
}

/// A blind pocket whose footprint overhangs one face edge (an edge slot open to the side).
/// First overhang through the `PocketOnFace` pipe - routes to the overhang Cut sidecar. OCCT
/// scores the slotted solid.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn overhang_pocket_matches_occt() {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    let mut model = Model::new();
    let __w10 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w10,
            profile: Profile2d::polygon(
                [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
                    .iter()
                    .map(|&p| nacre_math::Point2::from_array(p))
                    .collect(),
            )
            .unwrap(),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let slot = Profile2d::polygon(
        [[-0.25, -0.25], [0.75, -0.25], [0.75, 0.25], [-0.25, 0.25]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut model,
        &Operation::PocketOnFace {
            face: faces[1],
            profile: slot,
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
}
