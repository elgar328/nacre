//! Pocketed solids, holed profiles, non-convex overhangs, the fin array and the angle sweeps,
//! against OCCT.

use super::*;

/// The unit cube with a `0.4`-square pocket `0.5` deep in its top face — the fixture
/// whose lid carries an inner loop. OCCT reads the same solid from STEP; nothing here
/// asks it to reproduce `PocketOnFace`.
fn pocketed_cube() -> (Model, Handle<Solid>) {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    let prof = |pts: &[[f64; 2]]| {
        Profile2d::polygon(
            pts.iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut m = Model::new();
    let __w1 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w1,
            profile: prof(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: faces[1],
            profile: prof(&[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]),
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    (m, solid)
}

/// Score one boolean of a holed operand against OCCT, on volume and on area. Area
/// matters here: nacre's own gates all read the same rings, so only an independent
/// kernel makes the surviving hole's size a real claim.
fn diff_holed(name: &str, kind: OcctBool, boxes: [[f64; 3]; 2], swap: bool) {
    use nacre_ops::BoolKind;
    let (mut m, pc) = pocketed_cube();
    let bx = m.add_cuboid(Point3::from_array(boxes[0]), Point3::from_array(boxes[1]));
    let (x, y) = if swap { (bx, pc) } else { (pc, bx) };
    let occt = occt_boolean_of(&m, kind, x, y).unwrap();
    let bk = match kind {
        OcctBool::Cut => BoolKind::Cut,
        OcctBool::Fuse => BoolKind::Fuse,
        OcctBool::Common => BoolKind::Common,
    };
    let r = boolean_one(&mut m, bk, x, y).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "{name} volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "{name} area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// A boolean composing on a shape a boolean can make. The seam bites
/// the holed lid's corner, and the hole is placed inside the region left behind.
///
/// The box is the symmetric one: its vertical edge pierces the lid
/// at `(0.85, 0.85)`, on a fan diagonal from every apex. OCCT is asked the same question
/// on the same coordinates, and it does not fan.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_corner_cut_matches_occt() {
    diff_holed(
        "pocket corner cut",
        OcctBool::Cut,
        [[0.85, 0.85, 0.85], [1.15, 1.15, 1.15]],
        false,
    );
}

/// The same bite mirrored in `z`: the seam misses the lid, which rides out whole.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_bottom_corner_cut_matches_occt() {
    diff_holed(
        "pocket bottom corner cut",
        OcctBool::Cut,
        [[0.85, 0.85, -0.15], [1.15, 1.15, 0.15]],
        false,
    );
}

/// A corner box whose footprint overlaps the pocket, so its walls cross the lid's hole
/// rim. `∂(lid)` is two rings the seam threads into one notch, opening
/// the pocket to the outside. Only an independent kernel makes the absorbed hole's area a
/// real claim; OCCT is asked on the same coordinates. `Cut(pc, box) = 0.893`.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_rim_corner_cut_matches_occt() {
    diff_holed(
        "pocket rim corner cut",
        OcctBool::Cut,
        [[0.55, 0.55, 0.85], [1.15, 1.15, 1.15]],
        false,
    );
}

/// The same two solids named the other way — the box's face is arranged first and the
/// pocket's rim pierces it. `Cut(box, pc) = 0.081`.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_rim_corner_cut_either_way_matches_occt() {
    diff_holed(
        "pocket rim corner cut either way",
        OcctBool::Cut,
        [[0.55, 0.55, 0.85], [1.15, 1.15, 1.15]],
        true,
    );
}

/// A slab through the pocket between its floor and its lid. Its underside
/// carries the cube's cross-section as one loop with the pocket's nested inside it — a loop
/// within a loop. `Cut(slab, pc) = 1.488`, the pocket loop hung in the slab region as a hole
/// beside an island; the area is the independent claim that the nesting placed both rings.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_nests_pocket_cut_matches_occt() {
    diff_holed(
        "slab nests pocket cut",
        OcctBool::Cut,
        [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
        false,
    );
}

/// The other order: `Cut(pc, slab) = 0.668`, where the slab's dropped underside leaves no
/// region, so the cube section becomes an island carrying the pocket loop as its own hole.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_nests_pocket_cut_either_way_matches_occt() {
    diff_holed(
        "slab nests pocket cut either way",
        OcctBool::Cut,
        [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
        true,
    );
}

/// `Fuse` of the same two seals the pocket into an enclosed cavity: material
/// `2.408`, a void of `0.032`. OCCT builds the same hollow solid and its volume subtracts
/// the void while its area adds both surfaces — the independent claim that the seam path
/// assembled the void inward, not as a phantom outer piece.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_seals_pocket_fuse_matches_occt() {
    diff_holed(
        "slab seals pocket fuse",
        OcctBool::Fuse,
        [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
        false,
    );
}

/// A slab over the pocket, its underside below the pocket floor. `Cut` keeps the lid
/// as a reversed inside-B piece, hole and all — `flip` meeting `inner` for the first
/// time. `Fuse` drops it. The third scores the complement.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_cut_by_pocket_matches_occt() {
    diff_holed(
        "slab cut by pocket",
        OcctBool::Cut,
        [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
        true,
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn slab_and_pocket_fuse_matches_occt() {
    diff_holed(
        "slab and pocket fuse",
        OcctBool::Fuse,
        [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
        true,
    );
}

#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn pocket_cut_by_slab_matches_occt() {
    diff_holed(
        "pocket cut by slab",
        OcctBool::Cut,
        [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
        false,
    );
}

/// A rod drilled clean through the L's bar: the first genus-1 solid this kernel makes
/// Area is scored alongside volume — a tunnel's walls are area, and
/// nacre's own gates all read the same rings, so only an independent kernel makes the
/// hole's size a claim rather than a restatement.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn drilled_l_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, rod) = l_prism_and_box([0.3, 0.3, -0.5], [0.5, 0.6, 1.5]);
    let occt = occt_boolean_of(&m, OcctBool::Cut, l, rod).unwrap();
    let r = boolean_one(&mut m, BoolKind::Cut, l, rod).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "drilled L volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "drilled L area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The `Fuse`: the rod stands proud on both faces of the bar, and each of its walls
/// splits into the stub above and the stub below.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn l_and_rod_fuse_matches_occt() {
    use nacre_ops::BoolKind;
    let (mut m, l, rod) = l_prism_and_box([0.3, 0.3, -0.5], [0.5, 0.6, 1.5]);
    let occt = occt_boolean_of(&m, OcctBool::Fuse, l, rod).unwrap();
    let r = boolean_one(&mut m, BoolKind::Fuse, l, rod).unwrap();
    let nacre = mass_props(&m, r).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "l and rod fuse volume: {} vs {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "l and rod fuse area: {} vs {}",
        nacre.area,
        occt.area
    );
}

/// The non-convex overhang Cut deliverable (unified coplanar handler): a slot cut from a
/// top-pocketed cube, flush on the +x wall, breaking out the bottom, its top coplanar-disjoint
/// with the pocket floor. OCCT confirms the accepted volume (hand estimate 0.8575).
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn non_convex_overhang_cut_matches_occt() {
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    let sq = |pts: &[[f64; 2]]| {
        Profile2d::polygon(
            pts.iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut m = Model::new();
    let __w0 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w0,
            profile: sq(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let OpOutput::PocketOnFace { solid: pc, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: faces[1],
            // `[0.3, 0.7]²` of the lid. On a lid the sketch frame is the identity on world
            // x and y: the origin is the world origin projected onto `z = 1` and the axes are
            // `u = +x̂`, `v = +ŷ`.
            profile: sq(&[[0.3, 0.7], [0.3, 0.3], [0.7, 0.3], [0.7, 0.7]]),
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let slot = m.add_cuboid(
        Point3::from_array([0.75, 0.25, -0.25]),
        Point3::from_array([1.0, 0.75, 0.5]),
    );
    // OCCT ground truth from the inputs, before nacre supersedes them.
    let occt = occt_boolean_of(&m, OcctBool::Cut, pc, slot).unwrap();
    let solids = boolean(&mut m, BoolKind::Cut, pc, slot).unwrap();
    assert_eq!(solids.len(), 1);
    let nacre = mass_props(&m, solids[0]).unwrap();
    assert!(
        approx(nacre.volume, occt.volume),
        "volume {} vs occt {}",
        nacre.volume,
        occt.volume
    );
    assert!(
        approx(nacre.volume, 0.8575),
        "volume {} vs hand 0.8575",
        nacre.volume
    );
    assert!(
        approx(nacre.area, occt.area),
        "area {} vs occt {}",
        nacre.area,
        occt.area
    );
}

/// **The four-plane result, scored by a kernel that has never heard of our substrate.**
///
/// A bar spun 45° whose bottom corner edge lands exactly in the plane `x = 0.5` that a fused
/// block supplies: three planes share that edge's line, so every plane crossing it makes a
/// vertex with **four** planes through it. Naming such a vertex is what the arrangement could
/// not do until the concurrency work, and the coverage suite checks the answer against a hand
/// figure and against the models one ULP away.
///
/// Neither of those is independent. This is: OCCT builds the same cut from the same two
/// solids, through STEP, with its own arithmetic and its own idea of what a vertex is. A
/// degenerate case is exactly where two kernels are most likely to disagree, which is why the
/// one result born of a coincidence gets an outside opinion.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn the_four_plane_cut_matches_occt() {
    use nacre_math::Point2;
    use nacre_ops::{BoolKind, Operation, Profile2d, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let extrude = |m: &mut Model, poly: Vec<Point2>, z: f64, dist: f64| {
        let __g217 = datum_frame(
            m,
            SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
        );
        let out = apply(
            m,
            &Operation::Extrude {
                frame: __g217,
                profile: Profile2d::polygon(poly).unwrap(),
                dist,
            },
        )
        .expect("extrude");
        let nacre_ops::OpOutput::Extrude { solid, .. } = out else {
            unreachable!("extrude yields Extrude output")
        };
        solid
    };

    let mut m = Model::new();
    // The unit cube, plus a block fused on that carries the plane x = 0.5 up to z = 2.
    let cube = extrude(
        &mut m,
        vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)],
        0.0,
        1.0,
    );
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
    m.rebuild_adjacency();

    // A square-section bar (half-height 0.2) spun 45° about Y through (0.5, ·, 1): the corner
    // travels 0.2·√2 sideways and lands back on z = 1.
    let bar = extrude(
        &mut m,
        vec![p2(0.3, -0.5), p2(0.7, -0.5), p2(0.7, 1.5), p2(0.3, 1.5)],
        0.8,
        0.4,
    );
    m.rebuild_adjacency();
    let nacre_ops::OpOutput::Transform { solid: bar } = apply(
        &mut m,
        &Operation::Transform {
            solid: bar,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Y,
                pivot: [Rat::new(1, 2).unwrap(), Rat::from_int(0), Rat::from_int(1)],
                angle: Angle::from_deg(Rat::from_int(45)).unwrap(),
            }),
        },
    )
    .expect("transform") else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();

    // OCCT scores the same cut from the same operands — before nacre consumes them.
    let occt = occt_boolean_of(&m, OcctBool::Cut, target, bar).expect("occt cut");
    let r = boolean(&mut m, BoolKind::Cut, target, bar).expect("the four-plane cut builds");
    m.rebuild_adjacency();
    assert_eq!(r.len(), 1, "one solid");
    let ours = nacre_props::mass_props(&m, r[0]).expect("props");
    // **The fixture is reproduced by hand here, so pin it to the figure the coverage suite
    // knows** (`1 + 0.5` of target, `0.12` removed). Without this, building a *different*
    // model that both kernels agree on would read as a passing cross-check.
    assert!(
        (ours.volume - 1.38).abs() < 1e-9,
        "this is not the four-plane model: volume {}, expected 1.38",
        ours.volume
    );
    assert!(
        approx(ours.volume, occt.volume),
        "four-plane cut volume: nacre {} vs occt {}",
        ours.volume,
        occt.volume
    );
    let c = nacre_props::centroid(&m, r[0]).expect("centroid");
    for i in 0..3 {
        assert!(
            approx(c[i], occt.centroid[i]),
            "four-plane cut centroid axis {i}: nacre {} vs occt {}",
            c[i],
            occt.centroid[i]
        );
    }
}

/// **Placement, scored by a kernel that has never heard of a motion history.**
///
/// Two things this cell changed have no prior expectation for the suite to violate: a
/// placement that used to come back as *two* bodies now comes back as one, and "turn it, then
/// put it there" used to be refused outright. A suite cannot catch a reject becoming an
/// answer, so the answer gets an outside opinion — OCCT builds the same booleans from the same
/// operands, through STEP, with its own arithmetic.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn placement_matches_occt() {
    use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    let shift = |m: &mut Model, s, x: Rat| {
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::translation([x, Rat::from_int(0), Rat::from_int(0)]),
            },
        )
        .expect("translate") else {
            unreachable!("transform yields Transform output")
        };
        m.rebuild_adjacency();
        solid
    };

    // (a) The shared wall whose two f64 images differ by one ulp — including the offsets that
    // used to split the part in half, and dyadic ones that never did.
    for (n, d) in [
        (7i128, 11i128),
        (13, 23),
        (1, 3),
        (2, 5),
        (3, 10),
        (1, 2),
        (3, 1),
    ] {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        m.rebuild_adjacency();
        let a = shift(&mut m, a, Rat::new(n, d).expect("offset"));
        let b = shift(&mut m, b, Rat::new(n + d, d).expect("offset"));
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).expect("occt fuse");
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse");
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1, "{n}/{d}: one body");
        let ours = nacre_props::mass_props(&m, out[0]).expect("props");
        assert!(
            approx(ours.volume, occt.volume),
            "{n}/{d}: volume nacre {} vs occt {}",
            ours.volume,
            occt.volume
        );
        let c = nacre_props::centroid(&m, out[0]).expect("centroid");
        for i in 0..3 {
            assert!(
                approx(c[i], occt.centroid[i]),
                "{n}/{d}: centroid axis {i}: nacre {} vs occt {}",
                c[i],
                occt.centroid[i]
            );
        }
    }

    // (b) Turn it, then put it there — the whole sweep that used to be refused.
    for deg in [7i128, 17, 30, 45, 63] {
        for off in [[3i128, 0, 0], [5, -3, 2], [0, 4, 0]] {
            let mut m = Model::new();
            let base = m.add_cuboid(
                Point3::from_array([-2.0, -2.0, 0.0]),
                Point3::from_array([2.0, 2.0, 1.0]),
            );
            let tool = m.add_cuboid(
                Point3::from_array([-0.5, -0.5, -1.0]),
                Point3::from_array([0.5, 0.5, 2.0]),
            );
            m.rebuild_adjacency();
            let turn = |m: &mut Model, s, iso| {
                let OpOutput::Transform { solid } = apply(
                    m,
                    &Operation::Transform {
                        solid: s,
                        isometry: iso,
                    },
                )
                .expect("transform") else {
                    unreachable!("transform yields Transform output")
                };
                m.rebuild_adjacency();
                solid
            };
            let tool = turn(
                &mut m,
                tool,
                Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                }),
            );
            let place = Isometry::translation(off.map(Rat::from_int));
            let tool = turn(&mut m, tool, place);
            let base = turn(&mut m, base, place);
            let occt = occt_boolean_of(&m, OcctBool::Cut, base, tool).expect("occt cut");
            let out = boolean(&mut m, BoolKind::Cut, base, tool).expect("cut");
            m.rebuild_adjacency();
            let ours = nacre_props::mass_props(&m, out[0]).expect("props");
            assert!(
                approx(ours.volume, occt.volume),
                "{deg}° then {off:?}: volume nacre {} vs occt {}",
                ours.volume,
                occt.volume
            );
            let c = nacre_props::centroid(&m, out[0]).expect("centroid");
            for i in 0..3 {
                assert!(
                    approx(c[i], occt.centroid[i]),
                    "{deg}° then {off:?}: centroid axis {i}: nacre {} vs occt {}",
                    c[i],
                    occt.centroid[i]
                );
            }
        }
    }
}

/// **Reflection, scored by a kernel that has never heard of a motion history.**
///
/// A reflection is now a motion the chain records, which turned two things into answers the
/// suite has no prior expectation for: a wall reached by reflection and a wall reached by
/// translation used to be one ulp apart and split the part in two, and a mirrored *rotated*
/// solid used to be carried by conjugating its chain rather than extending it. Both changed
/// what the kernel decides, so both get an outside opinion.
///
/// The mirror planes are deliberately **non-dyadic** (`1/3`, `7/22`, `5/7`): `2c − x` is exact
/// for a dyadic `c`, so those are the cases where the reflection is recorded and the
/// definition — not the `f64` coordinate — is what answers. Dyadic planes ride along as
/// controls that must not have changed.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn reflection_matches_occt() {
    use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    let flip = |m: &mut Model, s, offset: Rat| {
        let OpOutput::Mirror { solid } = apply(
            m,
            &Operation::Mirror {
                solid: s,
                axis: Axis::X,
                offset,
            },
        )
        .expect("mirror") else {
            unreachable!("mirror yields Mirror output")
        };
        m.rebuild_adjacency();
        solid
    };
    let shift = |m: &mut Model, s, x: Rat| {
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::translation([x, Rat::from_int(0), Rat::from_int(0)]),
            },
        )
        .expect("translate") else {
            unreachable!("transform yields Transform output")
        };
        m.rebuild_adjacency();
        solid
    };
    let score = |m: &Model,
                 tag: &str,
                 ours: &[nacre_store::Handle<nacre_topo::Solid>],
                 occt: &OcctProps| {
        assert_eq!(ours.len(), 1, "{tag}: one body, got {}", ours.len());
        let mp = nacre_props::mass_props(m, ours[0]).expect("props");
        assert!(
            approx(mp.volume, occt.volume),
            "{tag}: volume nacre {} vs occt {}",
            mp.volume,
            occt.volume
        );
        let c = nacre_props::centroid(m, ours[0]).expect("centroid");
        for i in 0..3 {
            assert!(
                approx(c[i], occt.centroid[i]),
                "{tag}: centroid axis {i}: nacre {} vs occt {}",
                c[i],
                occt.centroid[i]
            );
        }
    };

    // (a) **The crossing.** One wall arrives by reflection, the other by translation, and the
    // two `f64` images differ in the last place. `[1, 2]` reflected in `x = p` puts its wall at
    // `2p − 1`; `[0, 1]` shifted by `2p − 1` puts its wall in the same place by another route.
    for (n, d) in [(1i128, 3i128), (7, 22), (5, 7), (1, 2), (3, 1)] {
        let p = Rat::new(n, d).expect("mirror plane");
        let t = p
            .checked_mul(Rat::from_int(2))
            .and_then(|q| q.checked_sub(Rat::from_int(1)))
            .expect("2p − 1");
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        m.rebuild_adjacency();
        let a = flip(&mut m, a, p);
        let b = shift(&mut m, b, t);
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).expect("occt fuse");
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse");
        m.rebuild_adjacency();
        score(&m, &format!("crossing {n}/{d} fuse"), &out, &occt);
    }

    // (b) **All three operations over a reflected operand that genuinely overlaps.** The
    // crossing above only ever asks about one shared wall; this asks about the whole
    // arrangement, with the mirrored solid on both sides of the operation.
    for (n, d) in [(1i128, 3i128), (7, 22), (5, 7)] {
        let p = Rat::new(n, d).expect("mirror plane");
        // `[1, 2]` reflected in `x = p` occupies `[2p − 2, 2p − 1]`; put `[0, 1]` half inside.
        let t = p
            .checked_mul(Rat::from_int(2))
            .and_then(|q| q.checked_sub(Rat::new(3, 2).expect("3/2")))
            .expect("2p − 3/2");
        for (kind, ok, name) in [
            (BoolKind::Fuse, OcctBool::Fuse, "fuse"),
            (BoolKind::Cut, OcctBool::Cut, "cut"),
            (BoolKind::Common, OcctBool::Common, "common"),
        ] {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([1.0, 0.0, 0.0]),
                Point3::from_array([2.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            m.rebuild_adjacency();
            let a = flip(&mut m, a, p);
            let b = shift(&mut m, b, t);
            let occt = occt_boolean_of(&m, ok, a, b).expect("occt boolean");
            let out = boolean(&mut m, kind, a, b).expect("boolean");
            m.rebuild_adjacency();
            score(&m, &format!("overlap {n}/{d} {name}"), &out, &occt);
        }
    }

    // (c) **Mirror of a rotated solid — the improper chain in production.** The tool carries
    // `[Rotate, Mirror]`, an odd parity, so every judgement among its own planes runs through
    // the canonicalised base frame. A sign error there is a confidently wrong answer, not a
    // slow one, and OCCT is what says it is not happening.
    for deg in [7i128, 30, 45, 63] {
        for (n, d) in [(1i128, 3i128), (5, 7)] {
            let mut m = Model::new();
            let base = m.add_cuboid(
                Point3::from_array([-2.0, -2.0, 0.0]),
                Point3::from_array([2.0, 2.0, 1.0]),
            );
            let tool = m.add_cuboid(
                Point3::from_array([-0.5, -0.5, -1.0]),
                Point3::from_array([0.5, 0.5, 2.0]),
            );
            m.rebuild_adjacency();
            let OpOutput::Transform { solid: tool } = apply(
                &mut m,
                &Operation::Transform {
                    solid: tool,
                    isometry: Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        pivot: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                    }),
                },
            )
            .expect("rotate") else {
                unreachable!("transform yields Transform output")
            };
            m.rebuild_adjacency();
            let tool = flip(&mut m, tool, Rat::new(n, d).expect("mirror plane"));
            let occt = occt_boolean_of(&m, OcctBool::Cut, base, tool).expect("occt cut");
            let out = boolean(&mut m, BoolKind::Cut, base, tool).expect("cut");
            m.rebuild_adjacency();
            score(
                &m,
                &format!("turned {deg}° mirrored {n}/{d} cut"),
                &out,
                &occt,
            );
        }
    }
}

/// **The fin array, scored by a kernel that has never heard of `SurfaceDef`.**
///
/// These eight arrangements did not build at all until surfaces carried their own provenance:
/// a boolean's result used to describe its faces by their rounded coordinates, so the third
/// fin's wall never merged with the wall the first two had already made. The coverage suite
/// pins the new volumes, but those numbers came out of the very engine the fix changed — and a
/// fix that turns a reject into an answer has no prior expectation to violate, which is exactly
/// the shape of change a suite cannot catch.
///
/// So OCCT builds the same three fuses from the same operands, and scores the result.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn the_fin_array_matches_occt() {
    use nacre_ops::{BoolKind, Operation, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    // The θ that each used to reject under: RunSplit, CoincidentNodes, UnreachedCell,
    // StraightAngle, LabelConflict ×3, TraceDeclined — the whole failing population.
    for theta in [16.0f64, 20.0, 23.0, 46.0, 50.0, 54.0, 59.0, 62.0] {
        let mut m = Model::new();
        let mut part = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([1.0, 1.0, 3.0]),
        );
        m.rebuild_adjacency();
        for deg in [theta, theta + 180.0, theta + 198.0] {
            let fin = m.add_cuboid(
                Point3::from_array([0.5, -0.2, 1.0]),
                Point3::from_array([4.0, 0.2, 3.0]),
            );
            m.rebuild_adjacency();
            let nacre_ops::OpOutput::Transform { solid: fin } = apply(
                &mut m,
                &Operation::Transform {
                    solid: fin,
                    isometry: Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        pivot: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(
                            Rat::new((deg * 100.0).round() as i128, 100).unwrap(),
                        )
                        .unwrap(),
                    }),
                },
            )
            .expect("transform") else {
                unreachable!("transform yields Transform output")
            };
            m.rebuild_adjacency();
            // OCCT scores this fuse from the same operands, before nacre consumes them.
            let occt = occt_boolean_of(&m, OcctBool::Fuse, part, fin).expect("occt fuse");
            part = boolean(&mut m, BoolKind::Fuse, part, fin).expect("the fin fuses on")[0];
            m.rebuild_adjacency();
            let ours = nacre_props::mass_props(&m, part).expect("props");
            assert!(
                approx(ours.volume, occt.volume),
                "θ={theta} fin {deg}°: volume nacre {} vs occt {}",
                ours.volume,
                occt.volume
            );
            let c = nacre_props::centroid(&m, part).expect("centroid");
            for i in 0..3 {
                assert!(
                    approx(c[i], occt.centroid[i]),
                    "θ={theta} fin {deg}°: centroid axis {i}: nacre {} vs occt {}",
                    c[i],
                    occt.centroid[i]
                );
            }
        }
    }
}

/// **The rotation sweep, scored fuse by fuse.**
///
/// The 30° step is where `sin 30° = ½` puts a rotated copy's corner exactly on one of the
/// original's face planes, and the union's vertex there is named by four concurrent planes. That
/// fuse *always* built — what it could not do was describe itself in its own surfaces, so the
/// next rotation refused it. Nothing ever checked whether the shape it built was **right**: the
/// defect was in the naming, the labelling sat one step away, and a wrong answer here would have
/// looked exactly like the correct one to every test that existed.
///
/// So OCCT scores every step of the sweep from the same operands, before nacre consumes them.
///
/// ★ 38°, 40° and 44° are here for the same reason at one remove: those are the sweeps that
/// only run at all because a ray is now allowed to graze a ring corner, and the containment
/// answers that ray gives are checked by **nothing else**. A hole hung on the wrong face is a
/// topology error, and volume and centroid see it.
#[test]
#[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
fn the_rotation_sweeps_match_occt() {
    for step in [30usize, 38, 40, 44] {
        sweep_matches_occt(step);
    }
}

fn sweep_matches_occt(step: usize) {
    use nacre_math::Point2;
    use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchFrame, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    let prism = |m: &mut Model, pts: &[[f64; 2]], axis: Axis, dist: f64| {
        let profile =
            Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap();
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("extrude") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid
    };
    let mut m = Model::new();
    let plate = prism(
        &mut m,
        &[
            [0.0, 0.0],
            [50.0, 0.0],
            [50.0, 25.0],
            [38.0, 25.0],
            [38.0, 50.0],
            [50.0, 50.0],
            [50.0, 75.0],
            [0.0, 75.0],
        ],
        Axis::Z,
        12.0,
    );
    let bar = prism(
        &mut m,
        &[[20.0, 12.0], [75.0, 12.0], [75.0, 37.0], [55.0, 37.0]],
        Axis::X,
        25.0,
    );
    let unit = boolean(&mut m, BoolKind::Fuse, plate, bar).expect("the part fuses")[0];
    m.rebuild_adjacency();

    // The template must stay live to be copied again, so the accumulator starts as a copy.
    let copy = |m: &mut Model, s| {
        let OpOutput::Copy { solid } = apply(m, &Operation::Copy { solid: s }).expect("copy")
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid
    };
    let mut part = copy(&mut m, unit);
    for deg in (step..360).step_by(step) {
        let c = copy(&mut m, unit);
        let OpOutput::Transform { solid: c } = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg as i128)).unwrap(),
                }),
            },
        )
        .expect("rotate a copy") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, part, c).expect("occt fuse");
        part = boolean(&mut m, BoolKind::Fuse, part, c).expect("the copy fuses on")[0];
        m.rebuild_adjacency();
        let ours = nacre_props::mass_props(&m, part).expect("props");
        assert!(
            approx(ours.volume, occt.volume),
            "{step}°/{deg}°: volume nacre {} vs occt {}",
            ours.volume,
            occt.volume
        );
        let c = nacre_props::centroid(&m, part).expect("centroid");
        for i in 0..3 {
            assert!(
                approx(c[i], occt.centroid[i]),
                "{step}°/{deg}°: centroid axis {i}: nacre {} vs occt {}",
                c[i],
                occt.centroid[i]
            );
        }
    }
}
