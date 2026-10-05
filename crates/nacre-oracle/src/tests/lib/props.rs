//! Mass properties against OCCT: volume, area, centroid, moved and mirrored solids, cylinders, pads
//! and pockets.

use super::*;

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
valid 0
";
    let p = OcctProps::parse(stdout).unwrap();
    assert_eq!(p.volume, 24.0);
    assert_eq!(p.area, 52.0);
    assert_eq!(p.faces, 6);
    assert_eq!(p.bbox_min, [0.0, 0.0, 0.0]);
    assert_eq!(p.bbox_max, [2.0, 3.0, 4.0]);
    assert_eq!(p.centroid, [1.0, 1.5, 2.0]);
    assert!(!p.valid);
    let malformed = stdout.replace("valid 0", "valid yes");
    assert!(matches!(
        OcctProps::parse(&malformed),
        Err(OracleError::Parse(_))
    ));
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
    nacre_ops::fixtures::cuboid(
        &mut model,
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
    let outer =
        nacre_ops::fixtures::cuboid(&mut model, Point3::origin(), Point3::from_array([10.0; 3]));
    let inner = nacre_ops::fixtures::cuboid(
        &mut model,
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

/// **A turned cylinder, scored by a second kernel.** Turned 37° about `x` through an off-origin
/// pivot, the cylinder has no world statement: its cache and its seam vertices are realized from
/// the motion chain. OCCT reads the STEP those realizations write; the volume and area must
/// match nacre's and the hand figures (`πr²h`, `2πr(r + h)`), and the box must agree.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_turned_cylinder_matches_occt() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let mut model = Model::new();
    let c = nacre_ops::fixtures::cylinder(
        &mut model,
        Point3::from_array([1.0, 2.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.5,
        2.0,
    );
    let Ok(nacre_ops::OpOutput::Transform { solid }) = nacre_ops::apply(
        &mut model,
        &nacre_ops::Operation::Transform {
            solid: c.solid,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::X,
                pivot: [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)],
                angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
            }),
        },
    ) else {
        panic!("the turn")
    };
    model.rebuild_adjacency();
    let ours = mass_props(&model, solid).unwrap();
    let occt = occt_props_of(&model).unwrap();
    let (vol, area) = (PI * 1.5 * 1.5 * 2.0, 2.0 * PI * 1.5 * (1.5 + 2.0));
    assert!(
        approx(ours.volume, occt.volume) && approx(ours.volume, vol),
        "volume: nacre {} occt {} hand {vol}",
        ours.volume,
        occt.volume
    );
    assert!(
        approx(ours.area, occt.area) && approx(ours.area, area),
        "area: nacre {} occt {} hand {area}",
        ours.area,
        occt.area
    );
    let (lo, hi) = nacre_props::bounds(&model, solid).unwrap();
    for i in 0..3 {
        assert!(
            (lo[i] - occt.bbox_min[i]).abs() < 1e-5 && (hi[i] - occt.bbox_max[i]).abs() < 1e-5,
            "axis {i}: nacre {lo:?}..{hi:?} vs occt {:?}..{:?}",
            occt.bbox_min,
            occt.bbox_max
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
    let s = nacre_ops::fixtures::cylinder(
        &mut model,
        Point3::origin(),
        Vector3::from_array([0.0, 3.0, 4.0]),
        2.0,
        5.0,
    )
    .solid;
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
    let a = nacre_ops::fixtures::cuboid(
        &mut model,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b = nacre_ops::fixtures::cylinder_with_seam(
        &mut model,
        Point3::from_array([1.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        4.0,
    )
    .solid;
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

/// **A tangent wall through a window in a boss, scored by a second kernel.** The slab's face
/// `x = 1` is tangent to the boss's cylinder along `x = 1, y = 0`, but only inside the window
/// [`nacre_ops::fixtures::windowed_boss`] cut through the lateral, so cutting the boss from the
/// slab leaves the slab's face whole — one body. OCCT is given the same two operands; and both
/// kernels must also agree with the hand figure, the slab less the cylinder's slice less what the
/// window took out of it (`w·√(1−w²) + asin(w) − w`), since the shape is the point here.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn a_tangent_wall_through_a_lateral_window_matches_occt() {
    for w in [0.3f64, 0.6] {
        for seam in [1.0, -1.0] {
            let mut m = Model::new();
            let boss = nacre_ops::fixtures::windowed_boss(&mut m, w, seam);
            let slab = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([-2.0, -2.0, 1.7]),
                Point3::from_array([1.0, 2.0, 2.3]),
            );
            m.rebuild_adjacency();
            let occt = occt_boolean_of(&m, OcctBool::Cut, slab, boss).unwrap();
            let r = boolean_one(&mut m, BoolKind::Cut, slab, boss).unwrap();
            let ours = mass_props(&m, r).unwrap();
            let window = w * (1.0 - w * w).sqrt() + w.asin() - w;
            let hand = 7.2 - (PI - window) * 0.6;
            assert!(
                approx(ours.volume, occt.volume) && approx(ours.volume, hand),
                "w {w} seam {seam}: volume nacre {} occt {} hand {hand}",
                ours.volume,
                occt.volume
            );
            assert!(
                approx(ours.area, occt.area),
                "w {w} seam {seam}: area nacre {} occt {}",
                ours.area,
                occt.area
            );
        }
    }
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
        let plate = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = nacre_ops::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            3.0,
            7.0,
        )
        .solid;
        m.rebuild_adjacency();
        let holed = boolean_one(&mut m, BoolKind::Cut, plate, hole).expect("the bore cuts");
        let boss = nacre_ops::fixtures::cuboid(
            &mut m,
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
/// `x = 0.5` is exactly `r` away. The kernel fuses it, and the question this asks is the same one
/// asked of a point tangency: **does a second
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
        let cube = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([-0.5; 3]),
            Point3::from_array([0.5; 3]),
        );
        let stud = nacre_ops::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([ax, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid;
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
        let plate = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = nacre_ops::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([4.0, 2.0, bz]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, -1.0]),
            0.5,
            1.0,
        )
        .solid;
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
    let a = nacre_ops::fixtures::cuboid(
        &mut model,
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 2.0, 1.0]),
    );
    let bore = |m: &mut Model, x: f64| {
        nacre_ops::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([x, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.25,
            3.0,
        )
        .solid
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
    use nacre_exact::{Axis, Rat};
    use nacre_math::Point2;
    use nacre_ops::{OpOutput, Operation, apply};

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
    nacre_ops::fixtures::cylinder_with_seam(
        &mut model,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
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
    // Rotate a boolean *result* (every vertex arrangement-named) and feed it back into a boolean —
    // the capability the provenance plane-witness enables. R = a cube minus a far-corner
    // octant; rotate R and a fresh severing slab by 30° about Z, then Cut. OCCT computes the
    // same Cut of the two rotated STEP solids; the volumes must agree (nacre's exact
    // arrangement vs OCCT's independent kernel).
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    use nacre_ops::{OpOutput, Operation, apply};
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
    let a = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0; 3]),
        Point3::from_array([3.0; 3]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    let r = xf(&mut m, r);
    let c = nacre_ops::fixtures::cuboid(
        &mut m,
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
    let solid = nacre_ops::fixtures::cuboid(
        &mut model,
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
    let solid = nacre_ops::fixtures::cylinder_with_seam(
        &mut model,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    )
    .solid;
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

    // Unit cube, then a 0.4-square boss of height 0.5 in the middle of the top face. The
    // padded solid's top face has a real hole (a FACE_BOUND in STEP); this checks OCCT reads
    // that holed boss as a closed solid and agrees on its volume (1.08) and area (6.8) with
    // nacre's analytic value.
    let mut model = Model::new();
    let __w13 = SketchFrame::world(&model, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut model,
        &Operation::Extrude {
            frame: __w13,
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
    let solid =
        nacre_ops::fixtures::pad(&mut model, faces[1], nacre_ops::fixtures::lid_square(), 0.5)
            .unwrap()
            .solid();

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
    assert!(
        approx(nacre.volume, 1.08) && approx(nacre.area, 6.8),
        "volume {}, area {} vs hand 1.08, 6.8",
        nacre.volume,
        nacre.area
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_diff_occt() {
    // `pocketed_cube`: a 0.4-square pocket of depth 0.5 in the middle of the top face. The
    // inward walls remove material; this checks OCCT reads the holed, concave solid and agrees
    // (volume 0.92, area 6.8) with nacre.
    let (model, solid) = nacre_ops::fixtures::pocketed_cube();

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
    assert!(
        approx(nacre.volume, 0.92) && approx(nacre.area, 6.8),
        "volume {}, area {} vs hand 0.92, 6.8",
        nacre.volume,
        nacre.area
    );
}

/// A pad whose footprint overhangs one face edge (a boss cantilever): the tool reaches past the
/// face, so the fuse meets an overhang. OCCT scores the cantilevered boss.
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
    // Footprint world x in [0.25,0.75], y in [-0.25,0.75] - overhangs the y=0 edge. On the lid
    // frame coordinates are world x and y (`nacre_ops::fixtures::lid_square`).
    let boss = Profile2d::polygon(
        [[0.25, 0.75], [0.25, -0.25], [0.75, -0.25], [0.75, 0.75]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let solid = nacre_ops::fixtures::pad(&mut model, faces[1], boss, 1.0)
        .unwrap()
        .solid();
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    // The whole boss stands above `z = 1`; area `5.625 + 0.5 + 0.125 + 3` — lid, boss top, the cantilever's underside, boss walls.
    assert!(
        approx(nacre.volume, 1.5) && approx(nacre.area, 9.25),
        "volume {}, area {} vs hand 1.5, 9.25",
        nacre.volume,
        nacre.area
    );
}

/// A blind pocket whose footprint overhangs one face edge (an edge slot open to the side): the
/// tool reaches past the face, so the cut meets an overhang. OCCT scores the slotted solid.
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
    // Footprint world x in [0.25,0.75], y in [-0.25,0.75] - overhangs the y=0 edge.
    let slot = Profile2d::polygon(
        [[0.25, 0.75], [0.25, -0.25], [0.75, -0.25], [0.75, 0.75]]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
    )
    .unwrap();
    let solid = nacre_ops::fixtures::pocket(&mut model, faces[1], slot, 0.5)
        .unwrap()
        .solid();
    let occt = occt_props_of(&model).unwrap();
    let nacre = mass_props(&model, solid).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{} vs {}",
        nacre.volume,
        occt.volume
    );
    // Only the on-face part `0.375` carves; area `4.75 + 0.625 + 0.375 + 1.0` — the faces below the lid less the slot's opening, lid, floor, walls.
    assert!(
        approx(nacre.volume, 0.8125) && approx(nacre.area, 6.75),
        "volume {}, area {} vs hand 0.8125, 6.75",
        nacre.volume,
        nacre.area
    );
}

/// ★★★★ **A seamless lateral is a valid shape to OCCT, in every shape the two-rim form takes.**
/// nacre writes a cylinder's lateral as its two rims and no seam edge; OCCT reads it, inserts the
/// seam it needs on import (`ShapeFix_Face::FixMissingSeam`) and must find the solid valid
/// (`checkshape`), with nacre's own face count and volume. The shapes: whole rims (a plain
/// cylinder, and a bore, which faces its axis); both rims cut, so the two rims' vertices stand on
/// different generators (`cut_at_both_rims`, both ways, the boxes off and on the plane through the
/// axis); a window on and off the rims' vertices (`windowed_boss`); a chain rim with a corner at
/// the whole rim's angle (a boss on a plate's edge); and a bore whose four slots cover every
/// angle. The helper's `valid` is calibrated both ways: this file's cylinder reads `1`, the same
/// file with the lateral's inner bound deleted reads `0` (`BRepCheck_NotClosed`).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn seamless_laterals_read_as_valid_solids() {
    use nacre_ops::{BoolKind, boolean, fixtures};
    let z = Vector3::from_array([0.0, 0.0, 1.0]);
    let mut shapes: Vec<(&str, Model)> = Vec::new();
    {
        let mut m = Model::new();
        fixtures::cylinder(&mut m, Point3::origin(), z, 2.0, 5.0);
        shapes.push(("a plain cylinder", m));
    }
    let bore = |m: &mut Model| {
        let plate = fixtures::cuboid(
            m,
            Point3::from_array([-2.0, -2.0, 0.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let pin = fixtures::cylinder(m, Point3::from_array([0.0, 0.0, -1.0]), z, 1.0, 6.0).solid;
        m.rebuild_adjacency();
        let s = boolean(m, BoolKind::Cut, plate, pin).expect("the bore")[0];
        m.rebuild_adjacency();
        s
    };
    {
        let mut m = Model::new();
        bore(&mut m);
        shapes.push(("a bore", m));
    }
    {
        let mut m = Model::new();
        let mut s = bore(&mut m);
        for (lo, hi) in [
            ([0.5, -0.8, 0.5], [3.0, 0.8, 1.0]),
            ([-0.8, 0.5, 1.3], [0.8, 3.0, 1.8]),
            ([-3.0, -0.8, 2.1], [-0.5, 0.8, 2.6]),
            ([-0.8, -3.0, 2.9], [0.8, -0.5, 3.4]),
        ] {
            let w = fixtures::cuboid(&mut m, Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            s = boolean(&mut m, BoolKind::Cut, s, w).expect("a slot")[0];
            m.rebuild_adjacency();
        }
        shapes.push(("a bore with four slots", m));
    }
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        for y_min in [-2.0, 0.0] {
            let mut m = Model::new();
            fixtures::cut_at_both_rims(&mut m, kind, y_min);
            m.rebuild_adjacency();
            shapes.push((
                match (kind, y_min < 0.0) {
                    (BoolKind::Fuse, true) => "both rims cut, fused",
                    (BoolKind::Fuse, false) => "both rims cut on the axis plane, fused",
                    (_, true) => "both rims cut, cut",
                    (_, false) => "both rims cut on the axis plane, cut",
                },
                m,
            ));
        }
    }
    for seam_x in [1.0, -1.0] {
        let mut m = Model::new();
        fixtures::windowed_boss(&mut m, 0.6, seam_x);
        m.rebuild_adjacency();
        shapes.push((
            if seam_x > 0.0 {
                "a window on the rims' vertices"
            } else {
                "a window off the rims' vertices"
            },
            m,
        ));
    }
    {
        let mut m = Model::new();
        let boss = fixtures::cylinder(&mut m, Point3::origin(), z, 1.0, 4.0).solid;
        let plate = fixtures::cuboid(
            &mut m,
            Point3::from_array([-2.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        m.rebuild_adjacency();
        boolean(&mut m, BoolKind::Fuse, boss, plate).expect("the boss on the edge");
        m.rebuild_adjacency();
        shapes.push(("a chain rim with a corner at the whole rim's angle", m));
    }
    assert_eq!(shapes.len(), 10, "the shapes");
    for (name, m) in &shapes {
        let [solid] = m.live_solids()[..] else {
            panic!("{name}: one solid");
        };
        let p = occt_props_of(m).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert!(p.valid, "{name}: OCCT's checkshape refuses it");
        let faces = m.shell(m.solid(solid).outer).faces.len();
        assert_eq!(p.faces, faces, "{name}: face count");
        let ours = nacre_props::mass_props(m, solid)
            .expect("mass props")
            .volume;
        assert!(
            approx(p.volume, ours),
            "{name}: OCCT {} vs nacre {ours}",
            p.volume
        );
    }
}
