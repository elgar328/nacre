//! Cylinders in the arrangement: the ruling trace, lateral and mixed parity, cut circles, caps,
//! seated circles.

use super::*;

/// **The production path, run on a cylinder-bearing model.**
///
/// A `.plane()` call on the production road with cylinder rows flowing into it —
/// the `work` class list, the decline witness — aborts the
/// kernel on the first drill, and reading the ~30 call sites does not find one. **Running
/// the road is what finds them**, which is why this test exists.
///
/// What it asserts is deliberately weak on geometry and strong on survival: the tracer
/// completes, every face it emits is still a plane face, and the
/// drilled cap carries its circular hole.
#[test]
fn the_production_road_survives_a_cylinder() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // a through-hole: circles on both box caps
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        class_owner,
        n_a,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
    );
    let (faces, _) = trace_result_faces(
        &m,
        BoolKind::Cut,
        a,
        b,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        // `Proved` on purpose: the reuse guard must be what turns the shortcut off, not the
        // caller. Without it `pass_through` would carry planes across and drop the lateral.
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .expect("the drill population traces");
    assert!(
        faces.iter().all(|f| matches!(f.surf, ClassIx::Plane(_))),
        "the plane arrangement emits plane faces only; bands come from `bands`"
    );
    // The box: 4 walls + 2 caps, and each cap carries the drill's circular hole.
    assert_eq!(faces.len(), 6, "{faces:?}");
    let holed = faces
        .iter()
        .filter(|f| {
            f.inner
                .iter()
                .any(|b| matches!(b, crate::draft::Bound::Circle { .. }))
        })
        .count();
    assert_eq!(holed, 2, "both caps are drilled: {faces:?}");
}

/// The rulings-road harness: the through-boss geometry
/// — a plate and a cylinder overlapping its full height, axis exactly
/// on the `x = 40` wall plane — assembled **past the standing gate** from production parts:
/// `plane_index_setup_inner` (the gate-free half) plus the cylinder table built the way the
/// gate builds it. Returns everything the two locks below read, and the record the
/// gate produces: `(wall class, cylinder)` marked not-proven-clear.
#[allow(clippy::type_complexity)]
fn armed_through_boss() -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::arrangement::PlaneSetup,
    usize,
) {
    armed_through_boss_z(-10.0, 50.0)
}

/// [`armed_through_boss`] with the boss's axial extent chosen: `z_lo` and height `h`. The
/// default runs through the plate (`−10`, `50`); a boss whose lower cap sits **inside** the
/// plate (`10`, `30`) has a lateral whose lower boundary is a **chain** — arcs at z = 10
/// (outside the plate) and z = 20 (inside) joined by rulings — the staircase the corner
/// rule is watched on.
#[allow(clippy::type_complexity)]
fn armed_through_boss_z(
    z_lo: f64,
    h: f64,
) -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::arrangement::PlaneSetup,
    usize,
) {
    armed_boss_at(20.0, z_lo, h)
}

/// [`armed_through_boss_z`] with the axis at `(40, y)`: still on the `x = 40` wall, at any height
/// along it.
#[allow(clippy::type_complexity)]
fn armed_boss_at(
    y: f64,
    z_lo: f64,
    h: f64,
) -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::arrangement::PlaneSetup,
    usize,
) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(&mut m, [0.0; 3], [40.0, 40.0, 20.0]);
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([40.0, y, z_lo]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        5.0,
        h,
    )
    .solid;
    let (mut setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(&m, plate, boss).unwrap();
    for &surf in &cyl_surfs {
        let nacre_topo::Surface::Cylinder { def, .. } = m.surface(surf) else {
            unreachable!("a cylinder row carries a cylinder truth")
        };
        let Some(cache) = m.surface_cache(surf).cylinder() else {
            unreachable!("push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache")
        };
        setup.cyls.push(crate::planes::WorkingCyl {
            surf,
            def: def.clone(),
            realized: *cache,
            owner: crate::planes::SolidSide::A,
            shared: Vec::new(),
            grazed: Vec::new(),
            records: Vec::new(),
            boxes: Vec::new(),
        });
    }
    let wc = setup
        .geom
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[0] - 40.0).abs() < 1e-12)
        })
        .expect("the x = 40 wall class");
    setup.cyls[0]
        .records
        .push((wc, crate::planes::CylRecord::Cuts));
    (m, plate, boss, setup, wc)
}

/// **A closed arc holds every point of its circle but its node** (`theta_between` with
/// `lo == hi`). On the boss's circle in the plate's top class, the wall `x = 40` through the axis
/// meets it at two named points: the closed arc at either holds the other, and not itself — which
/// the three-point order alone would refuse as two names for two points.
#[test]
fn a_closed_arc_holds_every_point_of_its_circle_but_its_node() {
    let (_m, _plate, _boss, setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let def = &setup.cyls[0].def;
    let top = setup
        .geom
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - 20.0).abs() < 1e-12)
        })
        .expect("the plate's top class z = 20");
    let at = |root| combinatorics::NodeId::pierce(wc, top, 0, root);
    let (lo, hi) = (
        at(crate::combinatorics::QuadRoot::Lo),
        at(crate::combinatorics::QuadRoot::Hi),
    );
    for (node, other) in [(lo, hi), (hi, lo)] {
        let between = |n| combinatorics::theta_between(&jd, 0, def, node, node, n);
        assert_eq!(between(other), Some(true), "{node:?} holds {other:?}");
        assert_eq!(between(node), Some(false), "{node:?} does not hold itself");
    }
}

/// The plate's top class `z = 20` of [`armed_boss_at`]`(5.0, ..)`, its wall class `y = 0`, and three
/// names on the boss's circle there: the wall `y = 0`'s touch point, and the wall `x = 40`'s crossing at it (`near` —
/// one point with the touch, under another name) and away from it (`far`).
fn touch_and_crossings(
    jd: &Judge<'_, WorkingPlane>,
    setup: &crate::arrangement::PlaneSetup,
    x40: usize,
) -> (usize, usize, [combinatorics::NodeId; 3]) {
    let class = |axis: usize, at: f64, what: &str| {
        setup
            .geom
            .iter()
            .position(|p| {
                p.witness_coords()
                    .iter()
                    .all(|q| (q.as_array()[axis] - at).abs() < 1e-12)
            })
            .expect(what)
    };
    let top = class(2, 20.0, "the plate's top class z = 20");
    let y0 = class(1, 0.0, "the plate's wall class y = 0");
    let touch = combinatorics::NodeId::pierce(top, y0, 0, crate::combinatorics::QuadRoot::Double);
    let at = |root| combinatorics::NodeId::pierce(top, x40, 0, root);
    let (lo, hi) = (
        at(crate::combinatorics::QuadRoot::Lo),
        at(crate::combinatorics::QuadRoot::Hi),
    );
    let (near, far) = match combinatorics::circular_order(jd, 0, &setup.cyls[0].def, &[touch, lo]) {
        Err(_) => (lo, hi),
        Ok(_) => (hi, lo),
    };
    (top, y0, [touch, near, far])
}

/// **A crossing the three-point question cannot place stays for the order, and the round records
/// its second name.** The boss's axis is at `(40, 5)`, so its circle in the plate's top class passes
/// the plate's corner `(40, 0, 20)`: the wall `y = 0` touches it there (the tangency, skipped as a
/// crossing) and the wall `x = 40` crosses it there and at `(40, 10, 20)`. An arc from the touch
/// point, named on `y = 0` as a fillet's arc end is named on its tangent wall, to the far crossing
/// has the near crossing as another name of its end, which no alias states — so the three-point
/// order of the arc's ends and that crossing is `Coincident` (both lie on the seam, `(40, 0)`), and
/// `theta_between` is `None`. Dropped as if it were off the arc, the point would stay under two
/// names — the arc ending at the touch point, the segment `x = 40` at the plate's corner. Kept, the
/// circle's θ order ties the two names, and the split (a step of the alias rounds) records them as
/// one point and stands; the next round, on the grown table, reads one name and learns nothing —
/// red when the tie is refused. The record is made twice here: the segment `x = 40` ends at the
/// plate's corner, so its split ties the near crossing with that end and records the point with
/// every name its planes give it on the cylinder, the touch point's included — the circle's own
/// record is not what this test can see (measured: skipping it leaves the test green). Control: the
/// same names folded as the table's corner seeding folds them ([`Aliases::record_on_cylinder`])
/// leave the split nothing to learn.
#[test]
fn an_arc_ends_unaliased_second_name_is_recorded_by_the_round() {
    let (m, plate, boss, setup, x40) = armed_boss_at(5.0, -10.0, 50.0);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let def = setup.cyls[0].def.clone();
    let (top, _, [touch, near, far]) = touch_and_crossings(&jd, &setup, x40);
    let tr = trace_on_class_of(&m, plate, boss, top, &jd, &setup.cyls, &setup);
    let (split, _, _) =
        class_skeleton(&jd, &setup.cyls, top, &tr, &mut Aliases::default()).unwrap();
    assert!(
        combinatorics::circular_order(&jd, 0, &def, &[touch, near]).is_err()
            && combinatorics::circular_order(&jd, 0, &def, &[touch, far]).is_ok(),
        "one crossing of x = 40 is the touch point and the other is not"
    );
    let round = |table: &mut Aliases| {
        let circle = MergedCircle {
            cyl: 0,
            merged: vec![(
                SolidSide::B,
                SegKind::Seated { body_above: true },
                Some([touch, far].map(|n| table.canon_point(n))),
            )],
        };
        split_class_edges(
            &jd,
            &setup.cyls,
            &setup.planes,
            top,
            &split,
            std::slice::from_ref(&circle),
            &[],
            table,
        )
        .map(|_| ())
    };
    let mut table = Aliases::default();
    let first = round(&mut table);
    assert!(first.is_ok(), "the round records the tie: {first:?}");
    assert_eq!(
        table.canon_point(near),
        table.canon_point(touch),
        "the end's second name is the end"
    );
    let learned = table.learned();
    let next = round(&mut table);
    assert!(next.is_ok(), "one point under one name: {next:?}");
    assert_eq!(table.learned(), learned, "the next round learns nothing");
    let mut seeded = Aliases::default();
    seeded
        .record_on_cylinder(&jd, &setup.cyls, touch, x40)
        .expect("one cylinder");
    assert_eq!(
        seeded.canon_point(near),
        seeded.canon_point(touch),
        "the corner seeding folds the two"
    );
    let learned = seeded.learned();
    let folded = round(&mut seeded);
    assert!(folded.is_ok(), "one point under one name: {folded:?}");
    assert_eq!(
        seeded.learned(),
        learned,
        "seeded, the split has nothing to learn"
    );
}

/// **An extent whose two ends are one point is refused, whatever the table knows**
/// ([`CircleTrace::arc`]). The names are [`touch_and_crossings`]': the touch point and the near
/// crossing are one point under two names. Blind, the table does not know them, and the circle's θ
/// order ties them with no segment in the split — the refusal is reached only through the circle's
/// own record of that tie, by which the place of the end the group does not keep is found. Known,
/// the merge hands the extent on under the table's one name (`p == q`), which a rank test alone
/// reads as covering nothing. Control: the touch point and the far crossing are one arc, covered by
/// its contribution.
#[test]
fn an_extent_whose_ends_are_one_point_is_refused_whatever_the_table_knows() {
    let (_, _, _, setup, x40) = armed_boss_at(5.0, -10.0, 50.0);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (top, _, [touch, near, far]) = touch_and_crossings(&jd, &setup, x40);
    let split = |ends: [combinatorics::NodeId; 2], table: &mut Aliases| {
        let circle = MergedCircle {
            cyl: 0,
            merged: vec![(
                SolidSide::B,
                SegKind::Seated { body_above: true },
                Some(ends.map(|n| table.canon_point(n))),
            )],
        };
        match split_circles(&jd, &setup.cyls, top, &[], &[circle], table) {
            Ok(Some((_, _, arcs, _))) => Ok(arcs.iter().map(|a| a.end).collect::<Vec<_>>()),
            Ok(None) => panic!("the extent's ends are fed into the split"),
            Err(BoolError::Rejected { reason, .. }) => Err(reason),
            Err(e) => panic!("{e:?}"),
        }
    };
    let mut known = Aliases::default();
    known
        .record_same_point(&jd, &setup.cyls, touch, near)
        .expect("one point");
    let refused = Err(RejectReason::CoincidentNodes);
    assert_eq!(
        [
            split([touch, near], &mut Aliases::default()),
            split([touch, near], &mut known),
            split([touch, far], &mut Aliases::default()),
        ],
        [refused.clone(), refused, Ok(vec![[touch, far]])],
        "two names the table does not know, two it knows, and an arc"
    );
}

/// **A segment whose two ends are one point is refused by the overlay, whatever the table knows**
/// ([`split_at_crossings`]). The plate's top class of [`armed_boss_at`]`(5.0, ..)` with one more
/// segment on the wall `x = 40`: the plate's corner `{top, y = 0, x = 40}` is one point with the
/// near crossing of the boss's circle ([`touch_and_crossings`]). Its ends are that corner twice,
/// the corner and the near crossing under a table that knows them as one, and the same two names
/// under one that does not — the overlay's order ties them. Each covers no piece of the line, and
/// dropped it would vanish without a word. Control: the corner and the far crossing are a segment
/// that cuts the plate's edge on `x = 40`.
#[test]
fn a_segment_whose_ends_are_one_point_is_refused_whatever_the_table_knows() {
    let (m, plate, boss, setup, x40) = armed_boss_at(5.0, -10.0, 50.0);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (top, y0, [_, near, far]) = touch_and_crossings(&jd, &setup, x40);
    let corner =
        combinatorics::NodeId::three_planes(crate::combinatorics::Canon3::three([top, y0, x40]));
    let skeleton = |ends: [combinatorics::NodeId; 2], table: &mut Aliases| {
        let mut tr = trace_on_class_of(&m, plate, boss, top, &jd, &setup.cyls, &setup);
        let pin = |n| combinatorics::pin_for(&jd, top, x40, n).expect("an end on the line");
        tr.segs.push(Seg {
            wall: x40,
            end: ends,
            end_h: ends.map(pin),
            solid: SolidSide::B,
            kind: SegKind::Seated { body_above: true },
        });
        match class_skeleton(&jd, &setup.cyls, top, &tr, table) {
            Ok((segs, _, _)) => Ok(segs.len()),
            Err(BoolError::Rejected { reason, .. }) => Err(reason),
            Err(e) => panic!("{e:?}"),
        }
    };
    let mut known = Aliases::default();
    known
        .record_same_point(&jd, &setup.cyls, corner, near)
        .expect("one point");
    let refused = Err(RejectReason::CoincidentNodes);
    assert_eq!(
        [
            skeleton([corner, corner], &mut Aliases::default()),
            skeleton([corner, near], &mut known),
            skeleton([corner, near], &mut Aliases::default()),
            skeleton([corner, far], &mut Aliases::default()),
        ],
        [refused, refused, refused, Ok(5)],
        "one name, two the table knows, two it does not, and a segment"
    );
}

/// **A ray through a closed arc's node crosses it there once** (`point_in_mixed_ring`,
/// `ArcSpan::AtBothEnds`). The ring is the boss's circle in the plate's top class, cut at one point
/// where the wall `x = 40` meets it — the arc from that node back to it. The class's ray runs along
/// `y` at the probe's `x`, so a probe on `x = 40` sends it through both of the wall's points: the
/// node, where the arc both leaves and arrives, and the other point, inside the arc. From the
/// circle's centre the ray crosses once (inside); from beyond the circle on either side — one of
/// the two runs the ray through both points, whichever way the chart points it — twice (outside),
/// whichever of the two points is the node.
#[test]
fn a_ray_through_a_closed_arcs_node_crosses_it_once() {
    let (_m, _plate, _boss, setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let def = setup.cyls[0].def.clone();
    let top = setup
        .geom
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - 20.0).abs() < 1e-12)
        })
        .expect("the plate's top class z = 20");
    let coeffs = jd.world_rat(top).expect("a named class");
    let ri = nacre_exact::Rat::from_int;
    for root in [
        crate::combinatorics::QuadRoot::Lo,
        crate::combinatorics::QuadRoot::Hi,
    ] {
        let node = combinatorics::NodeId::pierce(wc, top, 0, root);
        for ccw in [true, false] {
            let ring = [combinatorics::RingEdge {
                node,
                to: node,
                carrier: combinatorics::Carrier::Arc(Box::new(combinatorics::ArcCarrier {
                    cyl: 0,
                    def: def.clone(),
                    plane: top,
                    ccw,
                })),
                from_h: combinatorics::EndPin::Cylinder,
                to_h: combinatorics::EndPin::Cylinder,
            }];
            for (probe, want) in [
                ([40, 20, 20], true),
                ([40, 30, 20], false),
                ([40, 10, 20], false),
            ] {
                assert_eq!(
                    combinatorics::point_in_mixed_ring(
                        &jd,
                        &setup.cyls,
                        &coeffs,
                        &probe.map(ri),
                        &ring
                    ),
                    Some(want),
                    "{root:?} ccw = {ccw}, probe {probe:?}"
                );
            }
        }
    }
}

/// **The gate's record arms the rulings road — and only the record.** With the through-boss pair
/// listed, the boss's trace on the wall class is the rectangle: two rulings (one per side,
/// Pierce-named ends, no plain segments) and two cap chords (`Lo`/`Hi` roots on the cap classes).
/// With the pair's record removed the same trace is empty: the negative control that pins "an
/// empty record changes nothing", which is the population an unconditionally-firing arm would
/// break (14 arc-family tests red, measured). With a sectioning record instead, the same rulings
/// and chords.
#[test]
fn the_gates_record_arms_the_ruling_trace() {
    let (m, _plate, boss, mut setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut tr = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(
        tr.segs.len(),
        2,
        "one chord per cap — a segment of the sweep"
    );
    let mut sides: Vec<i8> = tr.rulings.iter().map(|r| r.side).collect();
    sides.sort_unstable();
    assert_eq!(sides, [-1, 1], "one ruling per side");
    for r in &tr.rulings {
        for nd in r.end {
            let (_, cyl, _) = combinatorics::pierce_name(nd).expect("a Pierce end");
            assert_eq!(cyl, 0);
        }
        assert!(matches!(r.kind, SegKind::Transversal { .. }));
    }
    for c in &tr.segs {
        assert_eq!(
            c.end_h,
            [combinatorics::EndPin::Cylinder; 2],
            "a chord's ends are the two roots"
        );
        // The chord rides the cap's own class — a ⊥ plane at z = −10 or z = 40.
        let z = setup.geom[c.wall].witness_coords()[0].as_array()[2];
        assert!(
            setup.geom[c.wall]
                .witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
                && ((z + 10.0).abs() < 1e-12 || (z - 40.0).abs() < 1e-12),
            "cap class at z = {z}"
        );
    }
    // The negative control: today's record.
    setup.cyls[0].records.retain(|&(c, _)| c != wc);
    let mut tr0 = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        &mut tr0,
    );
    assert!(
        tr0.rulings.is_empty() && tr0.segs.is_empty(),
        "empty record, empty road"
    );
    // A sectioning record draws the rulings and the chords: both ask whether the class draws the
    // rulings (`WorkingCyl::draws_rulings`). The gate sections only a class with no world
    // statement, so the record is set by hand on the armed (named) pair, and this pins the exact
    // roots; a sectioned class's interval roots are the scoreboard's
    // (`a_plane_along_the_axis_sections_a_cylinder_with_a_whole_rim`).
    setup.cyls[0].records.push((
        wc,
        crate::planes::CylRecord::Silent(crate::planes::Silence::Sectioned),
    ));
    let mut tr1 = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        &mut tr1,
    );
    assert_eq!(
        tr1.rulings.len(),
        2,
        "a sectioning record draws the rulings"
    );
    assert_eq!(tr1.segs.len(), 2, "and a chord per cap");
    for c in &tr1.segs {
        assert_eq!(
            c.end_h,
            [combinatorics::EndPin::Cylinder; 2],
            "a sectioned chord's ends are the two roots"
        );
    }
}

/// **The armed arrangement digests the rectangle** — production bricks end to end on the
/// wall class. Each ruling is cut at its T-junctions with the plate's `z = 0` and `z = 20`
/// lines (three pieces each), those lines split at the same Pierce nodes, the chords close
/// the far ends, and the walk closes the subdivision: one unbounded contour and five
/// bounded cells — plate-left, plate-right, the overlap band, and the rectangle's two
/// overhangs.
#[test]
fn the_armed_wall_class_walks_to_closed_cells() {
    let (m, plate, boss, setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let tr = trace_on_class_of(&m, plate, boss, wc, &jd, &setup.cyls, &setup);
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    let (split, circles, rulings) =
        class_skeleton(&jd, &setup.cyls, wc, &tr, &mut Aliases::default()).unwrap();
    // This fixture's class is parallel to the only cylinder in it, so it traces rulings and
    // no circle. Not a general law: a class may carry both, when the circle
    // belongs to a *different*, perpendicular cylinder.
    assert!(circles.is_empty(), "this class meets no ⊥ cylinder");
    assert_eq!(rulings.len(), 2);
    let edges = split_edges(
        &jd,
        &setup.cyls,
        &setup.planes,
        wc,
        &split,
        &circles,
        &rulings,
        &Aliases::default(),
    )
    .unwrap();
    let edges = ClassEdges::of(&setup.cyls, wc, &edges);
    assert_eq!(
        edges
            .pieces
            .iter()
            .filter(|p| matches!(p.on, OnCurve::Ruling { .. }))
            .count(),
        6,
        "each ruling cut at z = 0 and z = 20"
    );
    // The numbering's order — segments, then arcs, then ruling pieces ([`ClassEdges`]): the walk and
    // the emitter visit edges in it, and the assembly mints handles after them.
    let rank = |p: &Piece| match p.on {
        OnCurve::Segment { .. } => 0,
        OnCurve::Arc { .. } => 1,
        OnCurve::Ruling { .. } => 2,
    };
    assert!(
        edges.pieces.windows(2).all(|w| rank(&w[0]) <= rank(&w[1])),
        "pieces out of kind order: {:?}",
        edges.pieces.iter().map(rank).collect::<Vec<_>>()
    );
    let (cells, _face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 6, "five bounded cells and the outer");
    assert_eq!(
        cells.iter().filter(|c| c.winding == -1).count(),
        1,
        "one connected skeleton, one outer contour"
    );
}

/// One armed class's split edges, through the production bricks (the chain
/// `the_armed_wall_class_walks_to_closed_cells` spells out) — for the locks that need
/// several classes' products at once.
fn armed_class_edges(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &crate::arrangement::PlaneSetup,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
) -> ClassSplit {
    let tr = trace_on_class_of(m, plate, boss, wc, jd, &setup.cyls, setup);
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    class_split(
        jd,
        &setup.cyls,
        &setup.planes,
        wc,
        &tr,
        &mut Aliases::default(),
    )
    .unwrap()
}

/// **A ruling is asked against the plane its side was measured on, and a tangent one too.**
/// The line `x = 1, y = 0` beside a cylinder `r = 1` about `z` is held by two planes: `T: x = 1`,
/// tangent there, and `P: y = 0`, through the axis. A lateral ring runs up that line and back
/// (`z = 1 → 3`), its corners named on `T` (the cap with the tangent wall, `QuadRoot::Double`) —
/// so the corner name's ∥ class is `T` whatever the wall says. Two walls: (i) the chart's spelling
/// of a station on `P`, `side` the sign of `(x − o)·(m × n̂_P)` read off coordinates; (ii) the
/// tangent wall's own, `side 0`. The face goes through `comp_face`, as a component face does, so
/// `Ring::edges` carries the plane on. The point of the line at `z = 2` is on the boundary
/// (`None`) under both; the point of `P`'s other ruling at `z = 2` is not, which keeps the reading
/// from widening into "any line of `P`".
#[test]
fn a_lateral_ruling_is_asked_against_its_own_plane() {
    use nacre_exact::quad::{CylinderMeet, plane_plane_cylinder, plane_side};
    use nacre_exact::{Orient, Rat};
    let mut m = Model::new();
    let block = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 1.0], [1.0, 2.0, 3.0]);
    let boss = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        4.0,
    );
    let surf = m.face(boss.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let Some(realized) = m.surface_cache(surf).cylinder() else {
        panic!("a cylinder")
    };
    let cyls = vec![crate::planes::WorkingCyl {
        surf,
        def: def.clone(),
        realized: *realized,
        owner: crate::planes::SolidSide::A,
        shared: Vec::new(),
        grazed: Vec::new(),
        records: Vec::new(),
        boxes: Vec::new(),
    }];
    let rows = crate::planes::collect_planes(&m, block).expect("the block's face table");
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    // The class whose plane is `n · x = d` — read off its own coefficients.
    let class = |n: [i128; 3], d: i128| -> usize {
        (0..planes.len())
            .find(|&c| {
                let w = jd.world_rat(c).expect("a named class");
                let k = (0..3).find(|&i| n[i] != 0).expect("a normal");
                (0..3).all(|i| w[i].checked_mul(ri(n[k])) == w[k].checked_mul(ri(n[i])))
                    && w[3].checked_mul(ri(n[k])) == w[k].checked_mul(ri(-d))
            })
            .expect("the block holds that plane")
    };
    let (t, p, z1, z3) = (
        class([1, 0, 0], 1),
        class([0, 1, 0], 0),
        class([0, 0, 1], 1),
        class([0, 0, 1], 3),
    );
    let n1 = NodeId::pierce(z1, t, 0, crate::combinatorics::QuadRoot::Double);
    let n3 = NodeId::pierce(z3, t, 0, crate::combinatorics::QuadRoot::Double);
    // `side` against `P`, from coordinates: the sign of `(x − o) · (m × n̂)` at `(1, 0, 2)`.
    let wp = jd.world_rat(p).unwrap();
    let n_hat = Vector3::from_array([wp[0].to_f64(), wp[1].to_f64(), wp[2].to_f64()]);
    let s =
        Vector3::from_array([1.0, 0.0, 0.0]).dot(Vector3::from_array([0.0, 0.0, 1.0]).cross(n_hat));
    let side_p: i8 = if s > 0.0 { 1 } else { -1 };
    let zero = ri(0);
    let (pz, py, px) = (
        [zero, zero, ri(1), ri(-2)],
        [zero, ri(1), zero, zero],
        [ri(1), zero, zero, zero],
    );
    let CylinderMeet::Pair { line, s: roots } =
        plane_plane_cylinder(&pz, &py, &def.origin(), &def.dir(), def.r2()).unwrap()
    else {
        panic!("the line `y = 0, z = 2` crosses the cylinder twice")
    };
    for (name, plane, side) in [("P", p, side_p), ("T", t, 0)] {
        let wall = |up| crate::combinatorics::Wall::Ruling {
            cyl: 0,
            side,
            up,
            plane,
        };
        let face = LocalFace {
            surf: crate::planes::ClassIx::Cyl(0),
            outer: crate::draft::Bound::Ring(crate::draft::Ring::new(
                vec![n1, n3],
                vec![wall(true), wall(false)],
            )),
            inner: Vec::new(),
            flip: false,
        };
        let cf = crate::assembly::comp_face(&jd, &cyls, &face).expect("a component face");
        let combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
            panic!("a lateral's loops, got {:?}", cf.outer)
        };
        for root in &roots {
            let got = combinatorics::loop_parity(&jd, &def, loops, &line, root);
            match plane_side(&px, &line, root) {
                Orient::Positive => assert_eq!(
                    got,
                    Some(nacre_geom::intersect::RingSide::OnBoundary),
                    "{name}: (1, 0, 2) is on the ruling"
                ),
                _ => assert!(
                    matches!(
                        got,
                        Some(
                            nacre_geom::intersect::RingSide::Inside
                                | nacre_geom::intersect::RingSide::Outside
                        )
                    ),
                    "{name}: (−1, 0, 2) is no ruling of this face: {got:?}"
                ),
            }
        }
    }
}

/// **A lateral face is asked as a plane face is: a point on it, an edge across it.** The band of
/// the cylinder `r = 1` about `z` between its rims `z = 1` and `z = 3`, read through `comp_face` as
/// the self-touch arms read it, against points and lines named by a staircase of three boxes
/// (planes `x = 0, 1`, `y = −1, 0, 1, 3/2, 2`, `z = 0 … 3`). The truth is the band's: a point of the
/// cylinder between the rims is on the face, a point on a rim its boundary, a point beside or above
/// it off the face, a point off the surface off it whatever its height. The line `x = 1, z = 2`
/// touches the cylinder at `(1, 0, 2)`, inside the band: the segment from `y = −1` to `y = 1`
/// holds the touch, the one from `y = 1` to `y = 2` does not; the line `y = 0, z = 2` meets the
/// surface at `x = ±1`, and the segment from the axis to `x = 1` meets it only at its end, which is
/// the vertex's question.
#[test]
fn a_lateral_face_is_asked_for_points_and_edges_as_a_plane_face_is() {
    use crate::combinatorics::{Canon3, QuadRoot, named_point_in_face, segment_meets_lateral};
    use nacre_exact::Rat;
    use nacre_geom::intersect::RingSide;
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(&mut m, [0.0, -1.0, 0.0], [1.0, 2.0, 1.0]);
    let b = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 1.0], [1.0, 1.5, 2.0]);
    let c = crate::fixtures::cuboid(&mut m, [0.0, -1.0, 2.0], [1.0, 1.0, 3.0]);
    let ab = crate::boolean(&mut m, crate::BoolKind::Fuse, a, b).unwrap()[0];
    let stair = crate::boolean(&mut m, crate::BoolKind::Fuse, ab, c).unwrap()[0];
    let boss = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        5.0,
    );
    let surf = m.face(boss.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let Some(realized) = m.surface_cache(surf).cylinder() else {
        panic!("a cylinder")
    };
    let cyls = vec![crate::planes::WorkingCyl {
        surf,
        def: def.clone(),
        realized: *realized,
        owner: crate::planes::SolidSide::A,
        shared: Vec::new(),
        grazed: Vec::new(),
        records: Vec::new(),
        boxes: Vec::new(),
    }];
    let rows = crate::planes::collect_planes(&m, stair).expect("the staircase's face table");
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    // The class whose plane is `n · x = d / k` — read off its own coefficients.
    let class = |n: [i128; 3], d: i128, k: i128| -> usize {
        (0..planes.len())
            .find(|&c| {
                let w = jd.world_rat(c).expect("a named class");
                let i = (0..3).find(|&i| n[i] != 0).expect("a normal");
                (0..3).all(|j| w[j].checked_mul(ri(n[i])) == w[i].checked_mul(ri(n[j])))
                    && w[3].checked_mul(ri(n[i]).checked_mul(ri(k)).unwrap())
                        == w[i].checked_mul(ri(-d))
            })
            .expect("the staircase holds that plane")
    };
    let x = |v: i128| class([1, 0, 0], v, 1);
    let y = |v: i128, k: i128| class([0, 1, 0], v, k);
    let z = |v: i128| class([0, 0, 1], v, 1);
    let face = LocalFace {
        surf: crate::planes::ClassIx::Cyl(0),
        outer: crate::draft::Bound::Rim {
            plane: z(1),
            ccw: true,
        },
        inner: vec![crate::draft::Bound::Rim {
            plane: z(3),
            ccw: false,
        }],
        flip: false,
    };
    let cf = crate::assembly::comp_face(&jd, &cyls, &face).expect("a component face");
    let at = |p: [usize; 3]| NodeId::three_planes(Canon3::three(p));
    for (what, p, want) in [
        ("on the band", [x(1), y(0, 1), z(2)], RingSide::Inside),
        (
            "on the lower rim",
            [x(1), y(0, 1), z(1)],
            RingSide::OnBoundary,
        ),
        (
            "on the upper rim",
            [x(1), y(0, 1), z(3)],
            RingSide::OnBoundary,
        ),
        (
            "on the cylinder below the band",
            [x(1), y(0, 1), z(0)],
            RingSide::Outside,
        ),
        ("on the axis", [x(0), y(0, 1), z(2)], RingSide::Outside),
        (
            "outside the cylinder",
            [x(1), y(1, 1), z(2)],
            RingSide::Outside,
        ),
    ] {
        assert_eq!(
            named_point_in_face(&jd, &cyls, at(p), &cf).expect("a lateral reads no error"),
            Some(want),
            "{what}"
        );
    }
    let crate::combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
        panic!("a lateral's loops, got {:?}", cf.outer)
    };
    use crate::combinatorics::EdgeOnLateral;
    // The touch itself is the line's double root, the point the arm reads.
    let touch = NodeId::pierce(x(1), z(2), 0, QuadRoot::Double);
    for (what, line, ends, want) in [
        (
            "a tangent line's segment over the touch",
            [x(1), z(2)],
            [at([x(1), y(-1, 1), z(2)]), at([x(1), y(1, 1), z(2)])],
            EdgeOnLateral::AtRoots(vec![(touch, RingSide::Inside)]),
        ),
        (
            "a tangent line's segment beside it",
            [x(1), z(2)],
            [at([x(1), y(1, 1), z(2)]), at([x(1), y(2, 1), z(2)])],
            EdgeOnLateral::Clear,
        ),
        (
            "a segment ending on the surface",
            [y(0, 1), z(2)],
            [at([x(0), y(0, 1), z(2)]), at([x(1), y(0, 1), z(2)])],
            EdgeOnLateral::Clear,
        ),
    ] {
        assert_eq!(
            segment_meets_lateral(&jd, &cyls, 0, loops, line, ends),
            Some(want),
            "{what}"
        );
    }
    assert_eq!(
        named_point_in_face(&jd, &cyls, touch, &cf).unwrap(),
        Some(RingSide::Inside)
    );
}

/// **An arc on a slanted plane is read on its own plane** — a lateral ring of one closed arc on the
/// plane `z − x = 1`, slanted to the axis `z` of the cylinder `r = 1`, its node where that plane
/// meets the cap plane `z = 1` on the cylinder. The arc is a section: it climbs from that node, so
/// no one height places it; which side of it a point of the cylinder is on is which side of its
/// plane, on the point's generator. On the generator `x = 1` the section stands at `z = 2`: the point
/// `(1, 0, 2)` is on it, `(1, 0, 3)` above it (the ray up the axis crosses nothing) and `(1, 0, 1.5)`
/// below (it crosses the arc once). Read off its node's name instead, the arc would stand at the
/// cap's `z = 1` and all three would get the wrong answer. The planes come from a wedge whose
/// faces lie on `z = 1`, `x = 1` and the slanted plane.
#[test]
fn an_arc_on_a_slanted_plane_is_read_on_its_own_plane() {
    use crate::combinatorics::{Carrier, EndPin, LateralLoop, QuadRoot, RingEdge, loop_parity};
    use nacre_exact::Rat;
    use nacre_exact::quad::{CylinderMeet, QuadVal, plane_plane_cylinder};
    let mut m = Model::new();
    let profile = crate::Profile2d::polygon(
        [(1.0, 0.0), (1.0, 1.0), (2.0, 1.0)]
            .map(|(u, v)| nacre_math::Point2::from_array([u, v]))
            .to_vec(),
    )
    .expect("a triangle");
    let frame = SketchFrame::world(&m, Axis::Y);
    let crate::OpOutput::Extrude { solid: wedge, .. } = crate::apply(
        &mut m,
        &crate::Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("the wedge extrudes") else {
        unreachable!()
    };
    let boss = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        4.0,
    );
    let surf = m.face(boss.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let rows = crate::planes::collect_planes(&m, wedge).expect("the wedge's face table");
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    let class = |n: [i128; 3], d: i128| -> usize {
        (0..planes.len())
            .find(|&c| {
                let w = jd.world_rat(c).expect("a named class");
                let i = (0..3).find(|&i| n[i] != 0).expect("a normal");
                (0..3).all(|j| w[j].checked_mul(ri(n[i])) == w[i].checked_mul(ri(n[j])))
                    && w[3].checked_mul(ri(n[i])) == w[i].checked_mul(ri(-d))
            })
            .expect("the wedge holds that plane")
    };
    let (cap, slant) = (class([0, 0, 1], 1), class([-1, 0, 1], 1));
    let node = crate::combinatorics::NodeId::pierce(slant, cap, 0, QuadRoot::Lo);
    let ring = vec![RingEdge {
        node,
        to: node,
        carrier: Carrier::Arc(Box::new(crate::combinatorics::ArcCarrier {
            cyl: 0,
            def: def.clone(),
            plane: slant,
            ccw: true,
        })),
        from_h: EndPin::Cylinder,
        to_h: EndPin::Cylinder,
    }];
    use nacre_geom::intersect::RingSide;
    let x1 = [ri(1), ri(0), ri(0), ri(-1)];
    for (z, want) in [
        (Rat::from_int(2), RingSide::OnBoundary),
        (Rat::from_int(3), RingSide::Outside),
        (Rat::new(3, 2).unwrap(), RingSide::Inside),
    ] {
        let zero = Rat::from_int(0);
        let at_z = [ri(0), ri(0), ri(1), zero.checked_sub(z).unwrap()];
        let CylinderMeet::Tangent { line, s } =
            plane_plane_cylinder(&x1, &at_z, &def.origin(), &def.dir(), def.r2()).unwrap()
        else {
            panic!("the line x = 1, z = {z:?} touches the cylinder")
        };
        assert_eq!(
            loop_parity(
                &jd,
                &def,
                &[LateralLoop::Ring(ring.clone())],
                &line,
                &QuadVal::from_rat(s)
            ),
            Some(want),
            "the point at z = {z:?} against the section's arc"
        );
    }
}

/// **A lateral face answers the crossing question by the same parity, on its own chart**
/// — the digon oracle's lateral twin, on the through-boss Fuse's lateral. The
/// truth is known in world coordinates, so every crossing a rational ray makes with the
/// cylinder is checked: rays `{y = y₀, z = z₀}` over a half-step lattice (two roots each,
/// at x = 40 ± √(25 − (y₀ − 20)²) — irrational θ) and the **station column**
/// `{x = 40, z = z₀}`, whose roots are exactly the rulings' points (40, 15) and (40, 25):
/// on the ruling within the notch's z (a corner at its ends), on the face beyond it.
///
/// Two builds. **Through** (caps at −10 and 40): a band between the caps' whole circles
/// with the plate's notch as its one hole (the plate-side half-circle × z ∈ (0, 20), two
/// arcs on the cut rims and two rulings on the wall x = 40). **Staircase** (caps at 10 and
/// 40, the lower cap inside the plate): the lower boundary is a chain — the arc at z = 10
/// outside the plate, the arc at z = 20 inside, the two rulings between — so at the
/// station (40, 25) one arc ends `hi` and the other starts `lo`: the corner rule's one
/// discriminating population. In a notch both arcs share their ends, and «both ends
/// count» is invisible there — measured; the staircase is why the second build exists.
/// The other station, (40, 15), is the seam (the fixture puts it on −y): a root there is the
/// cyclic order's first point, and below the notch the notch's arcs above decide it like any
/// other.
///
/// Measured on the through build's rays: on the two whole-circle rims alone, «opposite
/// sides of the two rim planes» and the loops road's «exactly one rim above» agree on all
/// 2,180 rays — every root, the graze on a rim included — which is why no banded arm
/// exists.
#[test]
fn the_lateral_parity_agrees_with_the_notch_it_bounds() {
    for staircase in [false, true] {
        lateral_lattice(staircase);
    }
}

fn lateral_lattice(staircase: bool) {
    use nacre_exact::quad::{CylinderMeet, plane_plane_cylinder, plane_side};
    use nacre_exact::{Orient, Rat};
    let caps = if staircase {
        [10.0, 40.0]
    } else {
        [-10.0, 40.0]
    };
    let (m, plate, boss, setup, wc) = armed_through_boss_z(caps[0], caps[1] - caps[0]);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, _, rows, _) = armed_curved(&m, plate, boss, &setup, &jd, wc, caps);
    let fuse = crate::arrangement::cyl_chart::emit_lateral(
        BoolKind::Fuse,
        &jd,
        &setup.cyls,
        &curved,
        &rows,
    )
    .unwrap();
    assert_eq!(fuse.len(), 1);
    let cf = crate::assembly::comp_face(&jd, &setup.cyls, &fuse[0]).unwrap();
    let combinatorics::CompSurf::Cylinder(def) = &cf.surf else {
        panic!("a lateral")
    };
    let combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
        panic!("loops, got {:?}", cf.outer)
    };
    assert!(cf.inner.is_empty(), "a lateral's holes are among its loops");
    let rims: Vec<usize> = loops
        .iter()
        .filter_map(|l| match l {
            combinatorics::LateralLoop::Rim(c) => Some(*c),
            combinatorics::LateralLoop::Ring(_) => None,
        })
        .collect();
    if staircase {
        assert_eq!(rims.len(), 1, "the upper cap's whole circle");
        assert_eq!(loops.len(), 2, "and the chain below");
    } else {
        assert_eq!(rims.len(), 2, "the caps' two whole circles");
        assert_eq!(loops.len(), 3, "and the notch as one ring");
    }
    let (o, mm, r2) = (def.origin(), def.dir(), def.r2());
    let rat = |k: i128| Rat::new(k, 2).unwrap();
    let ri = Rat::from_int;
    let zero = ri(0);
    let neg = |v: Rat| zero.checked_sub(v).unwrap();
    let x40 = [ri(1), zero, zero, neg(ri(40))];
    let y20 = [zero, ri(1), zero, neg(ri(20))];
    let (cap_lo, cap_hi) = (ri(caps[0] as i128), ri(caps[1] as i128));
    // The notch's z range: the plate's own (0, 20) through the boss, or — with the lower
    // cap inside the plate — the chain's step from the cap (10) to the plate's top (20).
    let (n_lo, n_hi) = (if staircase { cap_lo } else { zero }, ri(20));
    // The truth on the cylinder, by the root's side of the wall and the ray's z — the seam
    // station (40, 15) included: a root there is ordered as the cyclic order's first point, so
    // the notch's arcs above it say whether they hold it like any other.
    let truth = |x_side: Orient, z0: Rat| -> Option<bool> {
        if z0 >= cap_hi {
            return if z0 == cap_hi { None } else { Some(false) };
        }
        if z0 < cap_lo {
            return Some(false);
        }
        if z0 == cap_lo {
            // The lower cap: through the plate a whole rim (a tie everywhere); in the
            // staircase an arc outside the plate only — a tie there, a corner on the wall,
            // and nothing inside, where the face starts at the plate's top.
            return match (staircase, x_side) {
                (true, Orient::Negative) => Some(false),
                _ => None,
            };
        }
        match x_side {
            // Outside the plate: the band, whole.
            Orient::Positive => Some(true),
            // Inside the plate's footprint: the notch (a hole, or the chain's step), its
            // arcs the boundary — both arcs through the plate, the upper one alone in the
            // staircase (its lower boundary there is the cap, handled above).
            Orient::Negative => {
                if (!staircase && z0 == n_lo) || z0 == n_hi {
                    None
                } else {
                    Some(!(z0 > n_lo && z0 < n_hi))
                }
            }
            // On the wall: a ruling for z within the notch (corners at its ends), the
            // face beyond.
            Orient::Zero => {
                if z0 >= n_lo && z0 <= n_hi {
                    None
                } else {
                    Some(true)
                }
            }
        }
    };
    // on face, off (hole), boundary, tangent, station on-ruling, station on-face, and the seam
    // station below the notch, where the notch's arcs above decide a root on the seam
    let mut n = [0usize; 7];
    let mut ask = |pa: [Rat; 4], pb: [Rat; 4], z0: Rat, station: bool| {
        let roots = match plane_plane_cylinder(&pa, &pb, &o, &mm, r2).unwrap() {
            CylinderMeet::Pair { line, s } => (line, s),
            CylinderMeet::Tangent { .. } => {
                n[3] += 1;
                return;
            }
            other => panic!("an ordinary ray: {other:?}"),
        };
        let (line, s) = roots;
        for root in &s {
            let x_side = plane_side(&x40, &line, root);
            let seam = station && plane_side(&y20, &line, root) == Orient::Negative;
            use nacre_geom::intersect::RingSide;
            // The truth's `None` is the boundary: every root here is one the loops can read.
            let want = Some(truth(x_side, z0).map_or(RingSide::OnBoundary, |inside| {
                if inside {
                    RingSide::Inside
                } else {
                    RingSide::Outside
                }
            }));
            let got = combinatorics::loop_parity(&jd, def, loops, &line, root);
            assert_eq!(
                got, want,
                "staircase {staircase} z0 {z0:?} side {x_side:?} station {station} seam {seam}"
            );
            match (
                got.expect("asserted equal to a decided truth"),
                station,
                seam,
            ) {
                (RingSide::Inside, false, _) => n[0] += 1,
                (RingSide::Outside, false, _) => n[1] += 1,
                (RingSide::OnBoundary, false, _) => n[2] += 1,
                (RingSide::Inside | RingSide::Outside, true, true) if z0 < n_lo => n[6] += 1,
                (RingSide::OnBoundary, true, _) => n[4] += 1,
                (_, true, _) => n[5] += 1,
            }
        }
    };
    for k in 0..=108 {
        let z0 = ri(-12).checked_add(rat(k)).unwrap();
        let pz = [zero, zero, ri(1), neg(z0)];
        for j in 0..=20 {
            let y0 = ri(15).checked_add(rat(j)).unwrap();
            ask([zero, ri(1), zero, neg(y0)], pz, z0, false);
        }
        ask(x40, pz, z0, true);
    }
    eprintln!(
        "lateral lattice (staircase {staircase}): on {} hole {} boundary {} tangent {} \
             station on-ruling {} station on-face {} seam decided {}",
        n[0], n[1], n[2], n[3], n[4], n[5], n[6]
    );
    assert!(
        n.iter().all(|&c| c > 0),
        "every arm has a population: {n:?}"
    );
}

/// The through-boss's curved carriers and plane faces, collected the production way
/// (`per_class` on each of the four ⊥ classes and the wall class), with the cylinder's rows
/// and the four ⊥ classes `[z0, z20, cap_lo, cap_hi]`. This mirrors production's fold by
/// hand (`trace_result_faces`' accumulation). A drift between them is not caught by
/// anything: a test would simply start measuring a map the boolean never builds.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn armed_curved(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &crate::arrangement::PlaneSetup,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    caps: [f64; 2],
) -> (
    Curved,
    Vec<LocalFace>,
    Vec<crate::bands::CylRow>,
    [usize; 4],
) {
    let z_class = |z: f64| -> usize {
        setup
            .geom
            .iter()
            .position(|p| {
                p.witness_coords()
                    .iter()
                    .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
            })
            .unwrap_or_else(|| panic!("a class at z = {z}"))
    };
    let (z0, z20, cap_lo, cap_hi) = (
        z_class(0.0),
        z_class(20.0),
        z_class(caps[0]),
        z_class(caps[1]),
    );
    // The production carriers: disk labels from every ⊥ class, cut rims from the cut ones.
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut plane_faces: Vec<LocalFace> = Vec::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    // The wall class too: the chart's rulings are the wall's pieces, and without them the
    // both-cut middle is one whole-circle cell the emitter cannot read per sector.
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in [z0, z20, cap_lo, cap_hi, wc] {
        let edges = armed_class_edges(m, plate, boss, setup, jd, c);
        let edges = ClassEdges::of(&setup.cyls, c, &edges);
        let staged = per_class(jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        plane_faces.extend(staged.faces);
    }
    let curved = Curved {
        disk_labels,
        arc_labels,
        split_rims: cut_rims,
        rulings,
    };
    let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
    (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi])
}

/// **A cut circle is a band boundary**. On the armed through-boss, the per-class products of the
/// four ⊥ classes are collected the production way (`per_class` on each), and the chart must break
/// the lateral at the two **cut** circles (z = 0, z = 20) as well as the rims: three intervals, not
/// one full-height band (`boundary_lines` and the chart's lines agree on the four). The middle
/// interval is both-cut and is emitted as panel rings. Three controls on the reader: blinding only
/// the arcs still answers, because the vertical lines speak in their place; blinding **both** of
/// the chart's axes (the arc labels and the rulings') leaves the reader no answer there; and
/// turning every ruling label inside out while the arcs stand makes two sides of each cell
/// contradict each other — `emit_lateral` refuses both by name, neither side believed over the
/// other: a cell no side reads — its blinded ends read no circle, which a present cell always has
/// (`CylinderStagesDisagree`) — and labels that disagree (`LabelConflict`).
#[test]
fn a_cut_circle_bounds_the_bands() {
    let (m, plate, boss, setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi]) =
        armed_curved(&m, plate, boss, &setup, &jd, wc, [-10.0, 40.0]);
    let (disk_labels, arc_labels, cut_rims) =
        (&curved.disk_labels, &curved.arc_labels, &curved.split_rims);
    assert!(cut_rims.contains_key(&(0, z0)) && cut_rims.contains_key(&(0, z20)));
    assert!(disk_labels.contains_key(&(0, cap_lo)) && disk_labels.contains_key(&(0, cap_hi)));
    assert_eq!(arc_labels[&(0, z0)].len(), 2, "two arcs, two sector labels");
    assert_eq!(arc_labels[&(0, z20)].len(), 2);
    assert_eq!(rows.len(), 1);
    // A cut circle is a band boundary: the chart's boundary rule names the two cut rims beside the
    // caps, and the chart's own lines are exactly those four.
    let plane_of = |c: usize| jd.world_rat(c).unwrap();
    // Low to high, as a slab orders its lines.
    let expect: Vec<[Rat; 4]> = [cap_lo, z0, z20, cap_hi].map(plane_of).to_vec();
    let mut sorted = expect.clone();
    sorted.sort();
    assert_eq!(
        crate::arrangement::cyl_chart::boundary_lines(&jd, 0, &plane_faces, &curved, &rows)
            .unwrap(),
        sorted,
        "three intervals, cut circles included"
    );
    let chart = crate::arrangement::cyl_chart::chart_of(&jd, &setup.cyls, 0, &curved).unwrap();
    assert_eq!(
        chart.slabs[0]
            .iter()
            .map(|&l| chart.line_planes[l])
            .collect::<Vec<_>>(),
        expect,
        "the chart's lines are the boundaries and nothing else here"
    );
    // **The emitter answers with regions.** Under `keep` the wall splits the middle
    // interval into two sectors, one kept and one not, and the kept one joins the whole
    // bands above and below it: **fuse** emits one face — a band bounded by its two whole rims
    // (`Bound::Rim`) with the unkept sector as its one **hole** (a 4-node ring: two arcs on the
    // cut rims' own nodes, two rulings); **cut** keeps the other sector alone — one 4-node
    // `Ring` face. Which sector is the outer one is the assembly's and the volume oracle's to
    // measure (the gate-opening cell), not this harness's to re-derive.
    let faces_for = |kind: BoolKind| -> Vec<LocalFace> {
        crate::arrangement::cyl_chart::emit_lateral(kind, &jd, &setup.cyls, &curved, &rows).unwrap()
    };
    let (fuse, cut) = (faces_for(BoolKind::Fuse), faces_for(BoolKind::Cut));
    assert_eq!(fuse.len(), 1, "fuse: one lateral face, a band with a hole");
    assert_eq!(cut.len(), 1, "cut: one lateral face, the kept sector");
    let rim = &cut_rims[&(0, z0)];
    // The ccw arc a ring carries on the lower cut rim, as an ordered node pair.
    let arc_on_z0 = |r: &crate::draft::Ring| -> [NodeId; 2] {
        let n = r.nodes.len();
        assert_eq!(n, 4, "two arcs and two rulings");
        let arcs = r
            .walls
            .iter()
            .filter(|w| matches!(w, crate::combinatorics::Wall::Arc { .. }))
            .count();
        assert_eq!(arcs, 2);
        for i in 0..n {
            let (a, b) = (r.nodes[i], r.nodes[(i + 1) % n]);
            if let crate::combinatorics::Wall::Arc { ccw, .. } = r.walls[i]
                && rim.nodes.contains(&a)
                && rim.nodes.contains(&b)
            {
                return if ccw { [a, b] } else { [b, a] };
            }
        }
        panic!("no arc on the lower cut rim");
    };
    assert!(
        matches!(fuse[0].outer, crate::draft::Bound::Rim { ccw: true, .. }),
        "fuse: a band between the caps, its lower whole rim outside, got {:?}",
        fuse[0].outer
    );
    let (rims, holes): (Vec<_>, Vec<_>) = fuse[0]
        .inner
        .iter()
        .partition(|b| matches!(b, crate::draft::Bound::Rim { .. }));
    assert!(
        matches!(rims[..], [crate::draft::Bound::Rim { ccw: false, .. }]),
        "the upper whole rim inside: {rims:?}"
    );
    let [crate::draft::Bound::Ring(hole)] = holes[..] else {
        panic!("the unkept sector is the band's one hole, a ring: {holes:?}");
    };
    let crate::draft::Bound::Ring(panel) = &cut[0].outer else {
        panic!("cut: the kept sector is a ring, got {:?}", cut[0].outer);
    };
    assert!(cut[0].inner.is_empty());
    // **Fuse's hole is cut's panel**: the sector fuse drops (inside the plate) is exactly
    // the sector cut keeps (the groove's wall), so the two rings carry the same ccw arc on
    // the lower cut rim — «complementary panels», stated for regions.
    let (pf, pc) = (arc_on_z0(hole), arc_on_z0(panel));
    assert_eq!(pf, pc, "fuse's hole is cut's panel");
    // Negative control: without the sector labels the both-cut interval's cells have no
    // speaking end, and the emitter refuses the class rather than guess a chamber. Called
    // directly, not through the census — this hand-broken input violates the very premise
    // (`src0_present == 0`) the census asserts where the fact is made.
    let mut blind = Curved {
        disk_labels: curved.disk_labels.clone(),
        arc_labels: HashMap::new(),
        split_rims: curved.split_rims.clone(),
        rulings: curved.rulings.clone(),
    };
    blind.arc_labels.clear();
    // **Blinding one axis is not blinding the reader** — the chart has two, and the rulings
    // answer where the rims are silent (the vertical read's own cell). So this says the weaker,
    // truer thing: with the arcs gone the emitter still gets an answer, and
    // it is only when **both** axes are blinded that it refuses by name.
    assert!(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &blind,
            &rows,
        )
        .is_ok(),
        "the vertical lines answer what the blinded rims cannot"
    );
    for v in blind.rulings.values_mut() {
        for r in v.iter_mut() {
            r.label = None;
        }
    }
    assert!(matches!(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &blind,
            &rows,
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::CylinderStagesDisagree,
            ..
        })
    ));
    // **A wall that contradicts the rims is not outvoted.** With the arcs as they are and every
    // ruling's inside label inside out — the cell beside it inside the cylinder, the one a cell's
    // chamber reads — the walls read the opposite chamber on every cell they bound, and the reader
    // answers neither.
    let mut contrary = Curved {
        disk_labels: curved.disk_labels.clone(),
        arc_labels: curved.arc_labels.clone(),
        split_rims: curved.split_rims.clone(),
        rulings: curved.rulings.clone(),
    };
    for v in contrary.rulings.values_mut() {
        for r in v.iter_mut() {
            r.label = r.label.map(|l| l.map(|b| !b));
        }
    }
    assert!(matches!(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &contrary,
            &rows,
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::LabelConflict,
            ..
        })
    ));
}

/// The armed through-boss fuse **assembled to the end of the road** — the pipeline
/// [`the_armed_assembly_welds_the_panels`] walks, extracted so the consumer locks (props' θ-range
/// integrals, tess's open-rim merge) measure the very same solid through one spelling: every
/// class's arrangement, the coplanar unify, the band/panel pass, the seam table, `reconstruct`.
/// Returns the model with the one welded solid live, plus the setup and the wall class.
fn armed_assembled_through_boss() -> (Model, crate::arrangement::PlaneSetup, usize) {
    let (mut m, plate, boss, setup, wc) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut faces: Vec<LocalFace> = Vec::new();
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in 0..setup.geom.len() {
        let edges = armed_class_edges(&m, plate, boss, &setup, &jd, c);
        let edges = ClassEdges::of(&setup.cyls, c, &edges);
        let staged = per_class(&jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        // This mirrors production's fold by hand (`trace_result_faces`' accumulation). A
        // drift between them is not caught by anything: the test would simply start measuring
        // a map the boolean never builds.
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        faces.extend(staged.faces);
    }
    let curved = Curved {
        disk_labels,
        arc_labels,
        split_rims: cut_rims,
        rulings,
    };
    let (faces, rims) = crate::arrangement::whole_result_faces(
        BoolKind::Fuse,
        &jd,
        &setup.cyls,
        faces,
        &curved,
        &setup.planes,
        &setup.plane_ix,
        setup.n_a,
    )
    .unwrap();
    let mut memo = crate::realize::PlaneMemo::default();
    let seam = seam_table(&m, &faces, &setup.cyls, &jd, &mut memo).unwrap();
    let out = crate::assembly::reconstruct(
        &mut m,
        &jd,
        &seam,
        &faces,
        &setup.cyls,
        &rims,
        crate::assembly::Tangencies::none(),
        &mut memo,
    )
    .unwrap();
    assert_eq!(out.len(), 1, "one welded solid");
    m.restore_live(out);
    (m, setup, wc)
}

/// **The armed assembly welds the panels** (the lock): the whole
/// through-boss fuse, from production parts past the standing gate — every class's
/// arrangement, the coplanar unify, the band/panel pass, the seam table, and
/// `reconstruct` — comes back one solid with a clean `validate`. The ruling edges are
/// pinned structurally: exactly two straight lateral edges (the kept outer panel's), each
/// used exactly twice, carriers stated as **the edge's own fact** — the cylinder and the
/// wall plane — and only a whole rim holds an `OnSeam` vertex: the outer panel's cut rims are
/// arcs between their pierces, with none at θ = 0.
#[test]
fn the_armed_assembly_welds_the_panels() {
    let (m, setup, wc) = armed_assembled_through_boss();
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    // The ruling edges, structurally: straight lateral edges = carriers {cylinder, a plane}
    // with a Line curve.
    let cyl_surf = setup.cyls[0].surf;
    let wall_surf = setup.geom[wc].surf;
    let reach = m.reachable();
    let mut rulings = 0;
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        let mixed = (e.surfaces[0] == cyl_surf) != (e.surfaces[1] == cyl_surf);
        if !mixed {
            continue;
        }
        let curve = m
            .derive_edge_curve(e.surfaces, e.vertices, |_| nacre_topo::EdgeGiven::NONE)
            .unwrap();
        if matches!(curve, nacre_geom::Curve::Line(_)) {
            rulings += 1;
            assert!(
                e.surfaces.contains(&wall_surf),
                "a ruling's plane carrier is the wall: {:?}",
                e.surfaces
            );
        }
    }
    assert_eq!(rulings, 2, "the kept outer panel's two rulings");
    // The boss's two disks keep their whole rims, each a closed edge at its `OnSeam` vertex;
    // the panel's cut rims have no vertex at θ = 0.
    let closed: std::collections::HashSet<_> = reach
        .edges
        .iter()
        .map(|&e| m.edge(e).vertices)
        .filter(|[a, b]| a == b)
        .map(|[a, _]| a)
        .collect();
    let on_seam: Vec<_> = reach
        .vertices
        .iter()
        .filter(|&&v| matches!(*m.vertex(v), nacre_topo::Vertex::OnSeam(_)))
        .collect();
    assert_eq!(on_seam.len(), 2, "the two whole rims' vertices");
    assert!(
        on_seam.iter().all(|v| closed.contains(*v)),
        "an OnSeam vertex ends only a whole rim"
    );
}

/// **The armed solid's mass properties are exact** (the props
/// instrument): `mass_props` on the assembled through-boss fuse, with the θ-range lateral
/// integrals live, answers the derived closed forms — volume `32000 + 1000π` (plate plus
/// the boss outside it), area `6200 + 425π` (walls 3000 + split x = 40 wall 600, caps
/// bitten `3200 − 25π`, boss disks `50π`, full bands `300π`, the outer half-panel `100π`).
/// The lifted probe's full-2π mis-answer (7849.34 = truth + the panel counted whole,
/// `+100π`) cross-checks the area derivation. This is the volume oracle standing while the
/// gate still refuses production input.
#[test]
fn the_armed_solids_mass_properties_are_exact() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let props = nacre_props::mass_props(&m, m.live_solids()[0]).unwrap();
    let volume = 32000.0 + 1000.0 * std::f64::consts::PI;
    let area = 6200.0 + 425.0 * std::f64::consts::PI;
    assert!(
        (props.volume - volume).abs() <= 1e-9 * volume,
        "volume {} != {volume}",
        props.volume
    );
    assert!(
        (props.area - area).abs() <= 1e-9 * area,
        "area {} != {area}",
        props.area
    );
}

/// **The armed solid tessellates watertight** (the tess
/// instrument): the open-rim merge triangulates the θ-panels of the assembled through-boss
/// fuse — the outer panel's rims are arcs from θ = 0 to π, one ruling on θ = 0 (the boss's
/// reference direction is `−y`, along the wall), so this one solid exercises the multi-polyline
/// open chain *and* an unwrap that starts on θ = 0 — and every
/// undirected triangle edge is used exactly twice.
#[test]
fn the_armed_solid_tessellates_watertight() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    let open = uses.values().filter(|&&n| n != 2).count();
    assert_eq!(open, 0, "the mesh is watertight");
}

/// **The mixed parity reads a bitten ring — exactly, on the production pieces.**
///
/// The straddling boss cuts the plate-top ring at `(40, 15)` and `(40, 25)`, so that ring
/// carries two pierce corners and one arc step; the overhang digon is chord + outer arc.
/// Probes are derived, not read back: the chart for `n = +z` picks `e1 = [0, −1, 0]`, so
/// the ray runs toward −y at fixed x. `[37, 30]` is inside and its ray crosses the **arc
/// twice** (`(37−40)² + (y−20)² = 25` → y = 16, 24) before the bottom edge — the arc arm is
/// what that probe measures, and a chord-minded parity would answer it wrong. `[37, 18]`
/// sits inside the bite (outside the face), `[50, 20]` outside everything; `[43, 20]` is
/// inside the overhang (one arc crossing at y = 16). Cells are identified structurally,
/// not by size: the outside cell (winding −1) owns the outer arc piece's twin, so that
/// piece's forward cell is the overhang digon; the other digon is the bite.
#[test]
fn the_mixed_parity_reads_a_bitten_ring() {
    with_bitten_rings(|jd, cyls, wc, big, overhang, bite| {
        let coeffs = jd.world_rat(wc).unwrap();
        let rat = |x: i128, y: i128| {
            [
                nacre_exact::Rat::from_int(x),
                nacre_exact::Rat::from_int(y),
                nacre_exact::Rat::from_int(20),
            ]
        };
        let ask = |ring: &[combinatorics::RingEdge], p: [nacre_exact::Rat; 3]| {
            combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring)
        };
        assert_eq!(ask(big, rat(12, 12)), Some(true), "plain interior");
        assert_eq!(
            ask(big, rat(37, 30)),
            Some(true),
            "interior whose ray crosses the inner arc twice"
        );
        assert_eq!(ask(big, rat(37, 18)), Some(false), "inside the bite");
        assert_eq!(ask(big, rat(50, 20)), Some(false), "outside everything");
        assert_eq!(
            ask(overhang, rat(43, 20)),
            Some(true),
            "inside the overhang"
        );
        assert_eq!(
            ask(overhang, rat(37, 18)),
            Some(false),
            "the bite is not the overhang"
        );
        assert_eq!(ask(bite, rat(37, 18)), Some(true), "inside the bite digon");
        assert_eq!(
            ask(bite, rat(43, 20)),
            Some(false),
            "the overhang is not the bite"
        );
    });
}

/// **The digon's truth is known independently, so the parity can be swept rather than
/// spot-checked.** A digon of a chord and an arc is `disk ∩ half-space`, and both halves are
/// exact rational predicates the arrangement already owns —
/// [`nacre_exact::quad::cylinder_radial_side`] and the sign of the wall's plane equation. A
/// grid over the circle's neighbourhood therefore checks **every** answer, and it is the only
/// control here that crosses the straight arm, the arc arm **and their junction**: the
/// non-mixed road's oracle cannot see an arc, and a whole circle's arm collapses to
/// `Ordering::Equal` on a single step.
///
/// **What it does and does not reach, measured.** On this fixture the chord's ends are the
/// seam and the tangent columns, and a ray along the chord meets **both** arc ends at once —
/// so this sweep is green under a flipped arc-end departure sign (two ends on one ray flip
/// together and keep the parity) and green under the straight arm's old corner abstention
/// (negative controls, both measured). It is the oracle for the mixed road as a
/// whole — the straight arm, the arc arm and their junction on ordinary roots — and the
/// corner-bitten plate (`a_root_at_an_arc_end_is_a_corner_on_the_ray`) is the watch on the
/// arc-end arm, where a ray meets one end alone, and the quarter bite
/// (`a_tangent_ray_at_an_arc_end_counts_the_side_the_arc_bends_to`) the watch on an end a tangent
/// ray meets.
#[test]
fn the_mixed_parity_agrees_with_the_digon_it_bounds() {
    use nacre_exact::{Orient, Rat};
    with_bitten_rings(|jd, cyls, wc, _big, overhang, bite| {
        let coeffs = jd.world_rat(wc).unwrap();
        let def = cyls
            .iter()
            .map(|c| &c.def)
            .find(|d| *d.r2() == nacre_exact::BigRat::from(Rat::from_int(25)))
            .expect("the bitten circle");
        let (mut swept, mut on_boundary) = (0usize, 0usize);
        let (mut abstained, mut inside_seen) = (0usize, 0usize);
        let mut on_chord = 0usize;
        for i in 60..=100 {
            for j in 20..=60 {
                let p = [
                    Rat::new(i.into(), 2).unwrap(),
                    Rat::new(j.into(), 2).unwrap(),
                    Rat::from_int(20),
                ];
                let radial = nacre_exact::quad::cylinder_radial_side(
                    &p,
                    &def.origin(),
                    &def.dir(),
                    def.r2(),
                );
                // The chord is the plate's wall `x = 40`; the bite keeps `x < 40`.
                let wall = p[0].checked_sub(Rat::from_int(40)).unwrap();
                if radial == Orient::Zero || wall == Rat::from_int(0) {
                    on_boundary += 1;
                    // **A point strictly inside the circle and *on* the chord is on
                    // both digons' boundary, and "inside" has no answer there.** These are the
                    // grid's sharpest points: the ray along the chart's first axis leaves such
                    // a point straight through **both arc endpoints** — the chord is a step
                    // lying along the ray, and the straight arm names the probe between its
                    // ends as the boundary before any arc is asked.
                    if radial == Orient::Negative && wall == Rat::from_int(0) {
                        for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                            assert_eq!(
                                combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                None,
                                "{who} must abstain on its own chord at {p:?}"
                            );
                        }
                        on_chord += 1;
                    }
                    continue;
                }
                let inside_bite = radial == Orient::Negative && wall < Rat::from_int(0);
                let inside_over = radial == Orient::Negative && wall > Rat::from_int(0);
                for (ring, want, who) in [
                    (bite, inside_bite, "bite"),
                    (overhang, inside_over, "overhang"),
                ] {
                    // Every answer is the truth, and an abstention off the boundary is counted
                    // and must be none (below): a ray through a ring corner is decided by the
                    // half-open rule, and a ray touching the circle by its touch.
                    match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                        Some(got) => {
                            assert_eq!(got, want, "{who} at {p:?}");
                            swept += 1;
                            inside_seen += usize::from(want);
                        }
                        None => abstained += 1,
                    }
                }
            }
        }
        // Neither vacuous nor all-outside: the grid straddles the circle and the chord, and
        // the two digons between them own a real interior.
        assert!(swept > 1_000, "answers swept: {swept}");
        assert!(inside_seen > 100, "points inside a digon: {inside_seen}");
        assert!(
            on_boundary > 0,
            "the grid meets the boundary: {on_boundary}"
        );
        // No point off the boundary abstains: the rays of the columns `x = 35, 45` touch the
        // circle, and a touch inside an arc's span counts nothing.
        assert_eq!(abstained, 0, "abstentions off the boundary");
        assert!(
            on_chord > 0,
            "points on the chord inside the circle: {on_chord}"
        );
    });
}

/// The bitten-plate fixture, handed to `f` as `(judge, cyls, class, big, overhang, bite)`.
/// The bitten fixture: one wall (`x = 40`) cuts the circle into two digons.
fn with_bitten_rings(
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    with_cut_rings([40.0, 40.0, 20.0], false, f);
}

/// A plate `[0, plate]` with a z-cylinder (r = 5, h = 10) standing on its top at
/// `(40, 20)`: the class carrying the cut circle, its judge, and the three bounded rings —
/// the bitten top, the overhang (disk ∖ plate) and the bite (disk ∩ plate). The plate's
/// extent picks the cut: `[40, 40, 20]` bites with one wall (two digons); `[40, 24, 20]`
/// with two — the corner `(40, 24)` sits inside the circle and both rings are trigons
/// whose arc ends at the rational `(37, 24)`, off the seam and off the tangent columns.
/// `seam_off` puts the seam on `+x`, at `(45, 20)` — on no ring corner — instead of on `−y`,
/// which puts it on the chord end `(40, 15)`.
fn with_cut_rings(
    plate: [f64; 3],
    seam_off: bool,
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(&mut m, [0.0; 3], plate);
    // The seam off every ring corner (`+x`, at `(45, 20)`) or on the chord end (`−y`, `(40, 15)`).
    let seam = if seam_off {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, -1.0, 0.0]
    };
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([40.0, 20.0, 20.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array(seam),
        5.0,
        10.0,
    )
    .solid;
    let setup = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    // The one class whose circle is cut: find it by running the split everywhere.
    let mut found = None;
    for wc in 0..setup.geom.len() {
        if !matches!(setup.plane_ix.get(wc), Some(ClassIx::Plane(_)) | None)
            && wc < setup.plane_ix.len()
        {
            continue;
        }
        if jd.world_rat(wc).is_none() {
            continue;
        }
        let tr = trace_on_class_of(&m, a, b, wc, &jd, &[], &setup);
        let Ok((split, circles, _)) =
            class_skeleton(&jd, &setup.cyls, wc, &tr, &mut Aliases::default())
        else {
            continue;
        };
        let Ok(edges) = split_edges(
            &jd,
            &setup.cyls,
            &[],
            wc,
            &split,
            &circles,
            &[],
            &Aliases::default(),
        ) else {
            continue;
        };
        let edges = ClassEdges::of(&setup.cyls, wc, &edges);
        // The arcs' piece indices, in order.
        let arcs: Vec<usize> = (0..edges.pieces.len())
            .filter(|&i| matches!(edges.pieces[i].on, OnCurve::Arc { .. }))
            .collect();
        if arcs.is_empty() {
            continue;
        }
        let (cells, face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
        // Four cells share this class: the bitten plate-top face (+1, many steps), the
        // outside (−1 — it also borders the arc, via the outer bulge), and two digons
        // (chord + arc): the overhang (disk minus plate) and the bite (disk ∩ plate).
        // Tell them apart structurally: the outside owns the outer piece's twin, so the
        // outer piece's forward cell is the overhang; the remaining digon is the bite.
        let outside = cells
            .iter()
            .position(|c| c.winding == -1)
            .expect("the outside cell");
        let bounded: Vec<usize> = (0..cells.len())
            .filter(|&i| cells[i].winding == 1)
            .collect();
        assert_eq!(
            bounded.len(),
            3,
            "the bitten top, the overhang and the bite"
        );
        let big_ix = *bounded
            .iter()
            .max_by_key(|&&i| cells[i].half_edges.len())
            .expect("the bitten top ring");
        let outer_ai = arcs
            .iter()
            .copied()
            .find(|&ai| face_of[&(2 * ai + 1)] == outside)
            .expect("the outer piece borders the outside");
        let overhang_ix = face_of[&(2 * outer_ai)];
        let bite_ix = *bounded
            .iter()
            .find(|&&i| i != big_ix && i != overhang_ix)
            .expect("the bite");
        assert_eq!(
            cells[overhang_ix].half_edges.len(),
            cells[bite_ix].half_edges.len(),
            "the same walls cut the overhang and the bite"
        );
        let ring = |ix: usize| -> Vec<combinatorics::RingEdge> {
            cells[ix]
                .half_edges
                .iter()
                .map(|&he| edges.edge_at(he))
                .collect()
        };
        found = Some((wc, ring(big_ix), ring(overhang_ix), ring(bite_ix)));
        break;
    }
    let (wc, big, overhang, bite) = found.expect("one class carries the cut circle");
    f(&jd, &setup.cyls, wc, &big, &overhang, &bite);
}

/// **A root at an arc's own end is a corner on the ray, and the arc's tangent there says
/// whether the arc counts it** — the half-open rule in the arc arm, on the one
/// fixture that reaches it.
///
/// The bitten fixture cannot: its chord `x = 40` runs along its lattice rays, so a ray meets both
/// of the chord's ends at once and the two flip together. The plate's second wall
/// `y = 24` puts the corner `(40, 24)` inside the circle, so the bite's arc ends at
/// `(37, 24)` — a 3-4-5 point, rational, off the seam `(40, 15)` and off the tangent
/// columns `x = 35, 45` — and the lattice column `x = 37` sends its rays through that end
/// **alone**, beside an ordinary second root at `(37, 16)`. That is the single-end crossing
/// a global sign error cannot hide: two ends on one ray flip together and keep the parity
/// (the bitten fixture's chord), one end on a ray does not. Every answer the sweep gives
/// is the truth, each ring abstains on its own boundary (its arc and the chords) and answers
/// «outside» on the rest of the circle, and the arm decided the end exactly once per point on the
/// shooting side of that column, for both rings — under both seam placements, so an end on the
/// seam decides the way an ordinary end does.
#[test]
fn a_root_at_an_arc_end_is_a_corner_on_the_ray() {
    use nacre_exact::{Orient, Rat};
    for seam_off in [false, true] {
        with_cut_rings(
            [40.0, 24.0, 20.0],
            seam_off,
            |jd, cyls, wc, _big, overhang, bite| {
                let coeffs = jd.world_rat(wc).unwrap();
                let def = cyls
                    .iter()
                    .map(|c| &c.def)
                    .find(|d| *d.r2() == nacre_exact::BigRat::from(Rat::from_int(25)))
                    .expect("the cut circle");
                let decided0 = combinatorics::tie_probe::arc_end_decisions_mine();
                let (mut swept, mut abstained, mut inside_seen, mut boundary) = (0usize, 0, 0, 0);
                let mut column = [0usize; 2];
                let rat = |k: i128| Rat::new(k, 2).unwrap();
                for i in 60..=100 {
                    for j in 20..=60 {
                        let p = [rat(i), rat(j), Rat::from_int(20)];
                        let radial = nacre_exact::quad::cylinder_radial_side(
                            &p,
                            &def.origin(),
                            &def.dir(),
                            def.r2(),
                        );
                        let (x, y) = (p[0], p[1]);
                        let (wx, wy) = (x == Rat::from_int(40), y == Rat::from_int(24));
                        // Both rings' shared boundary: the two chords, `x = 40` for `15 ≤ y ≤ 24`
                        // and `y = 24` for `37 ≤ x ≤ 40`, ends included — the ends are both arcs'
                        // ends too, and (40, 15) is the seam when `seam_off` is false.
                        let on_chord = (wx && y >= Rat::from_int(15) && y <= Rat::from_int(24))
                            || (wy && x >= Rat::from_int(37) && x <= Rat::from_int(40));
                        if on_chord {
                            boundary += 1;
                            for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                                assert_eq!(
                                    combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                    None,
                                    "{who} must abstain on its chord at {p:?}"
                                );
                            }
                            continue;
                        }
                        // **A ring's boundary is its own arc, not the whole circle.** The bite's
                        // arc is the circle on the plate's side (`x ≤ 40 ∧ y ≤ 24`), the
                        // overhang's the rest: on its own arc a ring abstains, off it the point is
                        // outside that ring — on the tangent columns `x = 35, 45` too, where a ray
                        // along `y` only touches the circle and counts nothing.
                        if radial == Orient::Zero {
                            boundary += 1;
                            let on_bite_arc = x <= Rat::from_int(40) && y <= Rat::from_int(24);
                            for (ring, own, who) in [
                                (bite, on_bite_arc, "bite"),
                                (overhang, !on_bite_arc, "overhang"),
                            ] {
                                let got =
                                    combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring);
                                if own {
                                    assert_eq!(got, None, "{who} must abstain on its arc at {p:?}");
                                } else {
                                    assert_eq!(
                                        got,
                                        Some(false),
                                        "{who}: the circle off its arc is outside it at {p:?}"
                                    );
                                }
                            }
                            continue;
                        }
                        let in_disk = radial == Orient::Negative;
                        let on_plate = x < Rat::from_int(40) && y < Rat::from_int(24);
                        for (k, (ring, want, who)) in [
                            (bite, in_disk && on_plate, "bite"),
                            (overhang, in_disk && !on_plate, "overhang"),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                                Some(got) => {
                                    assert_eq!(got, want, "{who} at {p:?}");
                                    swept += 1;
                                    inside_seen += usize::from(want);
                                    if x == Rat::from_int(37) && radial == Orient::Positive {
                                        column[k] += 1;
                                    }
                                }
                                None => abstained += 1,
                            }
                        }
                    }
                }
                assert!(swept > 1_000, "answers swept: {swept}");
                assert!(inside_seen > 100, "points inside a ring: {inside_seen}");
                assert!(boundary > 0, "the grid meets the boundary: {boundary}");
                // The `x = 37` column outside the circle: 24 lattice points, every one answered by
                // both rings — the shooting side through the end, the other side by a miss — and
                // the end decided exactly once per point on the shooting side: 12 × 2 rings.
                assert_eq!(column, [24, 24], "the single-end column is answered");
                let decided = combinatorics::tie_probe::arc_end_decisions_mine() - decided0;
                eprintln!(
                    "arc-end sweep (seam_off {seam_off}): swept {swept} abstained {abstained} \
                 boundary {boundary} arc-end decisions {decided}"
                );
                // The `x = 37` column's shooting side, 12 points × 2 rings, and the `x = 40`
                // column's rays, which meet `(40, 15)` alone (its other root `(40, 25)` is the
                // overhang's interior), 12 more per ring, one of them a call that then abstains at
                // the corner `(40, 24)` — under either seam placement: with the seam on the chord
                // end `(40, 15)` that end is the cyclic order's first point and decides the same
                // way as an ordinary one.
                assert_eq!(
                    decided, 48,
                    "the arc-end arm decided the single-end rays (seam_off {seam_off})"
                );
            },
        );
    }
}

/// **A tangent ray that meets an arc at its own end counts that end by the side the arc bends
/// to** — the arc leaves the end along the ray, so its tangent says nothing, and the section lies
/// on its centre's side.
///
/// The plate's corner `(40, 20)` is the circle's centre, so the bite is a quarter disk: an arc from
/// `(35, 20)` to `(40, 15)` and two chords on the plate's walls. The rays run along `y`, and the
/// lattice column `x = 35` is tangent to the circle at the bite's arc end `(35, 20)`: a ray of that
/// column meets the corner there, where the chord `y = 20` leaves toward the centre's side and the
/// arc leaves along the ray and bends to the same side. Both count the corner, so the point is
/// outside both rings; read with the centre's side flipped, the arc would not count and the point
/// would come out inside. The column `x = 45` is tangent to the overhang's arc inside its span —
/// under `seam_off`, on its seam — where the ray touches and counts nothing. Every answer the sweep
/// gives is the truth, no point off the boundary abstains, and the column `x = 35` decides the
/// arc's end on every ray through it.
#[test]
fn a_tangent_ray_at_an_arc_end_counts_the_side_the_arc_bends_to() {
    use nacre_exact::{Orient, Rat};
    for seam_off in [false, true] {
        with_cut_rings(
            [40.0, 20.0, 20.0],
            seam_off,
            |jd, cyls, wc, _big, overhang, bite| {
                let coeffs = jd.world_rat(wc).unwrap();
                let def = cyls
                    .iter()
                    .map(|c| &c.def)
                    .find(|d| *d.r2() == nacre_exact::BigRat::from(Rat::from_int(25)))
                    .expect("the cut circle");
                let int = Rat::from_int;
                let (mut swept, mut inside_seen, mut boundary) = (0usize, 0usize, 0usize);
                let mut tangent_end = [0usize; 2];
                for i in 60..=100 {
                    for j in 20..=60 {
                        let p = [Rat::new(i, 2).unwrap(), Rat::new(j, 2).unwrap(), int(20)];
                        let (x, y) = (p[0], p[1]);
                        let radial = nacre_exact::quad::cylinder_radial_side(
                            &p,
                            &def.origin(),
                            &def.dir(),
                            def.r2(),
                        );
                        let on_chord = (x == int(40) && y >= int(15) && y <= int(20))
                            || (y == int(20) && x >= int(35) && x <= int(40));
                        let on_bite_arc = x <= int(40) && y <= int(20);
                        let on_circle = radial == Orient::Zero;
                        for (k, (ring, own_arc, who)) in [
                            (bite, on_bite_arc, "bite"),
                            (overhang, !on_bite_arc, "overhang"),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            let decided0 = combinatorics::tie_probe::arc_end_decisions_mine();
                            let got =
                                combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring);
                            if on_chord || (on_circle && own_arc) {
                                boundary += 1;
                                assert_eq!(
                                    got, None,
                                    "{who} must abstain on its boundary at {p:?}"
                                );
                                continue;
                            }
                            let in_disk = radial == Orient::Negative;
                            let on_plate = x < int(40) && y < int(20);
                            let want = if k == 0 {
                                in_disk && on_plate
                            } else {
                                in_disk && !on_plate
                            };
                            assert_eq!(got, Some(want), "{who} at {p:?}");
                            swept += 1;
                            inside_seen += usize::from(want);
                            if x == int(35)
                                && combinatorics::tie_probe::arc_end_decisions_mine() > decided0
                            {
                                tangent_end[k] += 1;
                            }
                        }
                    }
                }
                assert!(swept > 1_000, "answers swept: {swept}");
                assert!(inside_seen > 100, "points inside a ring: {inside_seen}");
                assert!(boundary > 0, "the grid meets the boundary: {boundary}");
                // The rays of the column `x = 35` that reach `(35, 20)`: the twenty lattice points
                // on its shooting side, each decided at the arc's end by both rings.
                assert_eq!(
                    tangent_end,
                    [20, 20],
                    "the tangent end decided (seam_off {seam_off})"
                );
            },
        );
    }
}

/// **On a ring whose every corner is rational, the mixed road and the chart road are one
/// rule**: `point_in_ring_2d_rat` says Inside/Outside/OnBoundary and
/// `point_in_mixed_ring` `Some(true)`/`Some(false)`/`None`, and the pairs match at every
/// lattice point of every bounded ring of every class of two overlapping plates — squares
/// and L-shapes with reflex corners. A lattice at half steps shares a row with every corner
/// and every horizontal step, so every corner-on-the-ray configuration the half-open rule
/// distinguishes (both steps up, both down, one each, a step along the ray, the probe at
/// the corner) is on the sweep. A mixed road that abstained at such a corner would do so where
/// the chart road answers (measured over the whole suite's rational rings: 5 such corners, 0
/// disagreements); a whole-suite shadow of the two roads costs 61 % of the serial sweep and is
/// not kept — this lattice is the lock.
#[test]
fn the_two_roads_agree_on_every_rational_ring() {
    use nacre_exact::Rat;
    use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(&mut m, [0.0; 3], [4.0, 4.0, 2.0]);
    let b = crate::fixtures::cuboid(&mut m, [2.0, 2.0, 0.0], [6.0, 6.0, 2.0]);
    let setup = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (mut rings_swept, mut points, mut inside, mut boundary, mut reflex) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    for wc in 0..setup.geom.len() {
        let Some(coeffs) = jd.world_rat(wc) else {
            continue;
        };
        let tr = trace_on_class_of(&m, a, b, wc, &jd, &[], &setup);
        let Ok((split, circles, _)) =
            class_skeleton(&jd, &setup.cyls, wc, &tr, &mut Aliases::default())
        else {
            continue;
        };
        let Ok(edges) = split_edges(
            &jd,
            &setup.cyls,
            &[],
            wc,
            &split,
            &circles,
            &[],
            &Aliases::default(),
        ) else {
            continue;
        };
        let edges = ClassEdges::of(&setup.cyls, wc, &edges);
        let Ok((cells, _)) = walk_cells(&jd, &setup.cyls, wc, &edges) else {
            continue;
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let chart = combinatorics::Chart2dRat::of_normal(&n).unwrap();
        let (e1, e2) = chart.axes();
        // Axis-aligned classes: unit chart axes and a unit normal, so a chart point
        // `(X, Y)` lifts to `X·e1 + Y·e2 − d·n`.
        for e in [e1, e2, &n] {
            assert_eq!(nacre_exact::dot3_rat(e, e), Some(Rat::from_int(1)));
        }
        let lift = |xy: [Rat; 2]| -> [Rat; 3] {
            let mut p = [Rat::from_int(0); 3];
            for k in 0..3 {
                p[k] = xy[0]
                    .checked_mul(e1[k])
                    .unwrap()
                    .checked_add(xy[1].checked_mul(e2[k]).unwrap())
                    .unwrap()
                    .checked_sub(coeffs[3].checked_mul(n[k]).unwrap())
                    .unwrap();
            }
            p
        };
        for c in cells.iter().filter(|c| c.winding == 1) {
            let ring: Vec<combinatorics::RingEdge> =
                c.half_edges.iter().map(|&he| edges.edge_at(he)).collect();
            assert!(!combinatorics::ring_is_mixed(&ring), "a planar fixture");
            let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
            let ring2 = chart.ring(&jd, &nodes).expect("rational corners");
            if ring2.len() > 4 {
                reflex += 1;
            }
            let bound = |k: usize, lo: bool| -> i128 {
                let it = ring2.iter().map(|q| {
                    assert_eq!(q[k].denom(), 1, "integer corners");
                    q[k].numer()
                });
                if lo {
                    it.min().unwrap() - 1
                } else {
                    it.max().unwrap() + 1
                }
            };
            for xi in 2 * bound(0, true)..=2 * bound(0, false) {
                for yi in 2 * bound(1, true)..=2 * bound(1, false) {
                    let q = [Rat::new(xi, 2).unwrap(), Rat::new(yi, 2).unwrap()];
                    let p = lift(q);
                    assert_eq!(
                        chart.project(&p),
                        Some(q),
                        "the lift is the chart's inverse"
                    );
                    let plane = point_in_ring_2d_rat(q, &ring2);
                    let mixed =
                        combinatorics::point_in_mixed_ring(&jd, &setup.cyls, &coeffs, &p, &ring);
                    match (plane, mixed) {
                        (RingSide::Inside, Some(true)) => inside += 1,
                        (RingSide::Outside, Some(false)) => {}
                        (RingSide::OnBoundary, None) => boundary += 1,
                        other => panic!("class {wc} ring {ring2:?} at {q:?}: {other:?}"),
                    }
                    points += 1;
                }
            }
            rings_swept += 1;
        }
    }
    assert!(rings_swept >= 3, "rings swept: {rings_swept}");
    assert!(reflex > 0, "an L-shaped ring: {reflex}");
    assert!(
        inside > 0 && boundary > 0,
        "points {points} inside {inside} boundary {boundary}"
    );
    eprintln!("two roads: rings {rings_swept} points {points} inside {inside} boundary {boundary}");
}

/// The gated drill population's fixture: a `[0,2]³` box and an axis-aligned cylinder at
/// `(1,1)`, r=0.5 — every wall is a full unit from the axis, so the population gate passes
/// and [`plane_index_setup`] hands the arrangement bricks a cylinder-bearing table.
fn drilled(m: &mut Model, z0: f64, h: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = crate::fixtures::cuboid(m, [0.0; 3], [2.0; 3]);
    let b = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([1.0, 1.0, z0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        h,
    )
    .solid;
    (a, b)
}

/// The plane class whose defining triangle lies wholly at `z` — a z-cap.
fn class_at_z(planes: &[WorkingPlane], z: f64) -> usize {
    planes
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
        })
        .expect("a z-cap class")
}

/// A cylinder cap alone on its class: the disk's circular outer traces as a **seated
/// circle** — and the lateral, whose rim lies on this very plane, adds its **graze** beside it.
/// The cell bricks turn the pair into a disk `+1` / contour `−1` pair whose labels say "body on
/// the cap's inside", with no segments anywhere. The emitted face's outer is the circle itself
/// ([`crate::draft::Bound::Circle`]), the vocabulary the assembly consumes.
#[test]
fn a_cylinder_cap_is_a_seated_circle_and_a_disk_face() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // caps at z=-1 and z=3, clear of the box
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 3.0);

    let mut tr = Trace::default();
    trace_one_of(
        &m,
        b,
        SolidSide::B,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_b,
        &plane_ix,
        &mut tr,
    );
    assert!(
        tr.segs.is_empty(),
        "a circle owes the segment machinery nothing"
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    // The class's normal is the cap's own outward +z (the cap is its only member), so the
    // body — below z=3 — is `body_above: false`. **Two contributions, not one**: the cap is
    // seated here and the lateral's rim grazes here, and on this convex cap they agree.
    assert!(planes[wc].plane.normal().as_array()[2] > 0.0);
    let mut kinds: Vec<SegKind> = tr
        .circles
        .iter()
        .inspect(|c| assert!(c.cyl == 0 && c.solid == SolidSide::B, "{:?}", tr.circles))
        .map(|c| c.kind)
        .collect();
    kinds.sort_by_key(|k| matches!(k, SegKind::Graze { .. }));
    assert!(
        matches!(
            kinds[..],
            [
                SegKind::Seated { body_above: false },
                SegKind::Graze { body_above: false },
            ]
        ),
        "{:?}",
        tr.circles
    );

    let circles = merge_circles(&tr.circles, &Aliases::default()).unwrap();
    let edges = split_edges(&jd, &cyls, &[], wc, &[], &circles, &[], &Aliases::default()).unwrap();
    let edges = ClassEdges::of(&cyls, wc, &edges);
    let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
    // Pseudo-half-edges 0 and 1 (no segments): the disk (+1) and its contour (−1).
    assert_eq!(cells.len(), 2, "{cells:?}");
    assert_eq!(
        (cells[0].half_edges.as_slice(), cells[0].winding),
        (&[0][..], 1)
    );
    assert_eq!(
        (cells[1].half_edges.as_slice(), cells[1].winding),
        (&[1][..], -1)
    );
    let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.root_groups,
        vec![1],
        "the contour bounds the unbounded region"
    );
    assert!(nesting.holes.is_empty());
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    assert_eq!(labels[1], [false; 4]);
    assert_eq!(
        labels[0],
        [false, false, false, true],
        "inside the circle, B below the plane only"
    );
    // Fuse keeps below and not above across the disk → the disk is a result face, and its
    // outer boundary is the circle — no ring, no nodes.
    let (out, _, _) = emit_faces(
        BoolKind::Fuse,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "{out:?}");
    assert!(matches!(
        out[0].outer,
        crate::draft::Bound::Circle { cyl: 0 }
    ));
    assert!(out[0].inner.is_empty());
}

/// The through-drill's cap arrangement, brick by brick: the box cap is a seated 4-ring, the
/// lateral crosses the cap plane in a **transversal circle**, and the bricks nest the circle
/// as the cap cell's hole, label the disk with the cylinder straddling `W`, and emit — for
/// `Cut` — exactly one face: the cap with a circular hole. The per-kind `edge_mask`
/// difference is what the two labels measure (seated flips one side, transversal flips
/// both).
#[test]
fn a_drill_circle_is_a_hole_of_the_cap_ring() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // z∈[-1,3]: through both caps of the box
    let setup = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let wc = class_at_z(&setup.geom, 0.0);

    let tr = trace_on_class_of(&m, a, b, wc, &jd, &[], &setup);
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::B,
                kind: SegKind::Transversal { .. },
                // No hole in this band, so the mark is the whole circle.
                arc: None,
            }]
        ),
        "{:?}",
        tr.circles
    );

    let (split, circles, _) =
        class_skeleton(&jd, &setup.cyls, wc, &tr, &mut Aliases::default()).unwrap();
    assert_eq!(
        split.len(),
        4,
        "the cap ring alone — the circle owes the splitter nothing"
    );

    let edges = split_edges(
        &jd,
        &setup.cyls,
        &[],
        wc,
        &split,
        &circles,
        &[],
        &Aliases::default(),
    )
    .unwrap();
    let edges = ClassEdges::of(&setup.cyls, wc, &edges);
    let (cells, face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 4, "cap ±1 and circle ±1: {cells:?}");
    let at = |he: usize| cells.iter().position(|c| c.half_edges == [he]).unwrap();
    let (disk, contour) = (at(2 * split.len()), at(2 * split.len() + 1));
    let cap = cells
        .iter()
        .position(|c| c.winding == 1 && c.half_edges.len() == 4)
        .unwrap();

    let nesting = nest_cells(&jd, &setup.cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.holes.get(&cap).map(Vec::as_slice),
        Some(&[contour][..]),
        "the circle contour is the cap cell's hole"
    );

    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    // The class is the box's **bottom** cap, so which of `W`'s two sides carries the box is
    // read off the class normal rather than assumed: the box occupies z>0.
    let up = setup.geom[wc].plane.normal().as_array()[2] > 0.0;
    let box_side = usize::from(!up); // 0 = above `n_W`, 1 = below
    let mut cap_label = [false; 4];
    cap_label[box_side] = true;
    assert_eq!(
        labels[cap], cap_label,
        "outside the circle: box material on the z>0 side only"
    );
    assert_eq!(
        labels[contour], labels[cap],
        "a hole bounds its host's region"
    );
    let mut disk_label = cap_label;
    disk_label[2] = true;
    disk_label[3] = true;
    assert_eq!(
        labels[disk], disk_label,
        "inside the circle the cylinder straddles W, the box side is unchanged"
    );

    let (out, _, _) = emit_faces(
        BoolKind::Cut,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "one face: the drilled cap — {out:?}");
    assert_eq!(out[0].outer.expect_ring().len(), 4);
    assert!(
        matches!(out[0].inner[..], [crate::draft::Bound::Circle { cyl: 0 }]),
        "{:?}",
        out[0].inner
    );
}

/// The hole arm of the seated tracer: a face whose inner loop is a circle (a two-hole plate's input
/// shape, which a bore `Cut` builds in production) emits its polygon segments **and** a seated
/// circle that inherits the face's own body side.
#[test]
fn a_circular_hole_ring_traces_as_a_seated_circle() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0);
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 0.0);

    let input = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &[],
    );
    let (fp, fl) = input.faces[0]
        .iter()
        .find(|(fp, _)| matches!(plane_ix[*fp], ClassIx::Plane(c) if c == wc))
        .expect("the box's bottom cap sits on wc");
    let doctored = combinatorics::FaceLoops {
        outer: fl.outer.clone(),
        holes: Some(vec![combinatorics::LoopRing::Circle { cyl: 0 }]),
        cycles: None,
    };
    let mut tr = Trace::default();
    trace_one(
        &[(*fp, doctored)],
        SolidSide::A,
        wc,
        &jd,
        &[],
        &faces_tab,
        &plane_ix,
        &Default::default(),
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(tr.segs.len(), 4, "the polygon outer still emits its edges");
    let SegKind::Seated { body_above } = tr.segs[0].kind else {
        panic!("a face on wc traces seated: {:?}", tr.segs[0]);
    };
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::A,
                kind: SegKind::Seated { body_above: ba },
                arc: None
            }]
            if ba == body_above
        ),
        "the hole circle inherits the face's body side: {:?}",
        tr.circles
    );
}

/// The existence condition, negatively: a cap plane **outside** the lateral's rim span gets
/// no circle (and no decline — a miss, like a parallel plane), while the rim-interior cap
/// still does, and a rim-**coincident** plane carries the cap's seated circle rather than a
/// transversal one. A ghost circle here would corrupt every label on the class.
#[test]
fn no_ghost_circle_outside_the_rim_span() {
    let mut m = Model::new();
    // z∈[-1,1.5]: through the box's bottom cap, short of its top cap at z=2.
    let (a, b) = drilled(&mut m, -1.0, 2.5);
    let setup = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);

    // The box's top cap at z=2 is past the rim span [−1, 1.5]: no circle, no decline.
    let top = trace_on_class_of(&m, a, b, class_at_z(&setup.geom, 2.0), &jd, &[], &setup);
    assert!(top.circles.is_empty(), "{:?}", top.circles);
    assert!(top.declined.is_empty(), "{:?}", top.declined);

    // The bottom cap at z=0 is strictly inside the span: the transversal circle is there.
    let bottom = trace_on_class_of(&m, a, b, class_at_z(&setup.geom, 0.0), &jd, &[], &setup);
    assert_eq!(
        bottom
            .circles
            .iter()
            .filter(|c| matches!(c.kind, SegKind::Transversal { .. }))
            .count(),
        1,
        "{:?}",
        bottom.circles
    );

    // The cylinder's own top cap at z=1.5 (inside the box): the rim-coincident plane carries
    // **two** contributions — the cap's `Seated` circle and the lateral's `Graze`, because a
    // rim is a touch, not a miss. It adds no *transversal* twin (t = span end, not strictly
    // inside). On a **convex** cap the two agree about which side the body is on; the whole
    // point of carrying both is the reflex corner (a blind bore's ceiling) where they do not,
    // and `edge_mask`'s `Graze > Seated` then picks the right one.
    let cap = trace_on_class_of(&m, a, b, class_at_z(&setup.geom, 1.5), &jd, &[], &setup);
    assert!(cap.declined.is_empty(), "{:?}", cap.declined);
    let sides: Vec<(bool, bool)> = cap
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some((false, body_above)),
            SegKind::Graze { body_above } => Some((true, body_above)),
            SegKind::Transversal { .. } | SegKind::Tangent { .. } => None,
        })
        .collect();
    assert_eq!(sides.len(), 2, "seated and graze: {:?}", cap.circles);
    assert!(
        sides.iter().any(|(g, _)| *g) && sides.iter().any(|(g, _)| !*g),
        "one of each kind: {:?}",
        cap.circles
    );
    assert_eq!(
        sides[0].1, sides[1].1,
        "a convex cap's two contributions agree on the body's side: {:?}",
        cap.circles
    );
}

/// **The corner the graze exists for.** On a convex cap the lateral's rim-graze and the cap's
/// seated circle say the same thing, so carrying both changes nothing — the test above locks that.
/// Here they **disagree**: at a blind bore's ceiling the plate's material is above the cap while
/// the bore's wall hangs below it, a reflex dihedral in the plane. `edge_mask`'s `Graze > Seated`
/// precedence then picks the wall's side, which is the whole mechanism; without the graze, the
/// seated rule would flip the wrong label bit and the next boolean on that solid would come back
/// `CylinderGateUndecided`.
#[test]
fn at_a_blind_bores_ceiling_the_graze_and_the_seated_circle_disagree() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(&mut m, [0.0; 3], [10.0; 3]);
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([5.0, 5.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        3.0,
    )
    .solid;
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // A second operand so the arrangement runs on the bored solid as an operand.
    let boss = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 10.0], [2.0, 2.0, 12.0]);
    let setup = plane_index_setup(&m, bored, boss).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let tr = trace_on_class_of(
        &m,
        bored,
        boss,
        class_at_z(&setup.geom, 3.0), // the bore's ceiling
        &jd,
        &[],
        &setup,
    );
    let seated: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    let grazes: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Graze { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    assert_eq!(
        seated.len(),
        1,
        "the ceiling is seated here: {:?}",
        tr.circles
    );
    assert_eq!(grazes.len(), 1, "the wall grazes here: {:?}", tr.circles);
    assert_ne!(
        seated[0], grazes[0],
        "a reflex dihedral: the cap and the wall put the body on opposite sides — this is the \
             only place the precedence matters, and the only reason the graze is emitted at all"
    );
    // The wall hangs below the ceiling, and that is the side `edge_mask` believes.
    assert!(!grazes[0], "the bore's wall is below its ceiling");
}

/// Partial overlap is NOT merged — different endpoints mean different edges. a and b share
/// the y=1 plane; a's y=1 chord is x∈[0,2], b's is x∈[1,3] — overlapping on x∈[1,2] but not
/// coincident. They must stay separate.
#[test]
fn partial_overlap_is_not_merged() {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
    let b = crate::fixtures::cuboid(&mut m, [1.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
    let setup = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let wc = shared_cap_class(&m, a, b, &setup.surf_ix, &setup.planes, &setup.plane_ix);
    let tr = trace_on_class_of(&m, a, b, wc, &jd, &[], &setup);
    // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
    // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
    // seated≡transversal coincidence).
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default()).expect("merge");
    // The shared y=1 wall class (a face at y=1).
    let y1 = setup
        .geom
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("a y=1 face");
    // a's chord x∈[0,2] and b's chord x∈[1,3] ride y=1 but have different endpoints, so they
    // stay as two distinct MergedSegs. A merge that ignored extent would collapse them to one.
    let on_y1 = merged.iter().filter(|e| e.wall == y1).count();
    assert!(
        on_y1 >= 2,
        "a's and b's y=1 chords stay distinct (partial overlap not merged): {on_y1}"
    );
}

/// **An arc that crosses the cut line once on its way into an end lying on it** — the half disk
/// over the arc `(1, 0) → (−1, 0)` (unit radius, `z ∈ [0, 2]`) traced on the slanted class
/// `x + y = 1`, which runs through the arc's end `(1, 0)` and crosses it again at `(0, 1)`. Each
/// cap's ring is read by the walk: the end on the line is flanked by the arc's last piece, which
/// lies across it, so the end is a crossing and the arc's interior crossing is another — one chord
/// per cap, from `(1, 0)` to `(0, 1)`. Read off the arc's two ends, the end would look like a touch
/// and the chord would not be there. (The booleans on this pair stop at `RulingBoundNotYet` before
/// they reach the arrangement, so the trace is where this is read.)
#[test]
fn an_arc_crossing_once_into_its_end_on_the_line_leaves_the_chord() {
    use nacre_exact::Rat;
    let r = Rat::from_int;
    let mut m = Model::new();
    let extrude = |m: &mut Model, profile: crate::Profile2d| {
        let frame = SketchFrame::world(m, Axis::Z);
        let Ok(OpOutput::Extrude { solid, .. }) = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 2.0,
            },
        ) else {
            panic!("the profile extrudes")
        };
        solid
    };
    let half = crate::Ring2d::new(
        vec![[r(1), r(0)], [r(-1), r(0)]],
        vec![
            crate::Edge2d::Arc {
                center: [r(0), r(0)],
                r2: r(1),
                ccw: true,
            },
            crate::Edge2d::Line,
        ],
    )
    .expect("a half disk");
    let a = extrude(&mut m, crate::from_paths(vec![half]).unwrap().remove(0));
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let tri = crate::Profile2d::polygon(vec![p2(-1.0, 2.0), p2(3.0, -2.0), p2(3.0, 2.0)]).unwrap();
    let b = extrude(&mut m, tri);
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = (0..planes.len())
        .find(|&c| {
            planes[c]
                .witness_coords()
                .iter()
                .all(|q| (q.as_array()[0] + q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("the x + y = 1 class");
    let mut tr = Trace::default();
    trace_one_of(
        &m,
        a,
        SolidSide::A,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_a,
        &plane_ix,
        &mut tr,
    );
    // The lateral's ruling road declines this class (the booleans' `RulingBoundNotYet`); the
    // caps are what is read here.
    assert!(
        tr.declined
            .iter()
            .all(|&(fp, _)| matches!(plane_ix[fp], ClassIx::Cyl(_))),
        "{:?}",
        tr.declined
    );
    let point = |n: NodeId| -> [f64; 3] {
        combinatorics::node_coords_rat(&jd, n)
            .map(|p| p.map(|x| x.to_f64()))
            .or_else(|| combinatorics::pierce_point(&jd, 0, &cyls[0].def, n))
            .expect("a point")
    };
    let mut chords: Vec<[[f64; 3]; 2]> = tr
        .segs
        .iter()
        .map(|s| s.end.map(point))
        .filter(|[p, q]| (p[2] - q[2]).abs() < 1e-12)
        .map(|mut e| {
            e.sort_by(|p, q| p[0].total_cmp(&q[0]));
            e
        })
        .collect();
    chords.sort_by(|p, q| p[0][2].total_cmp(&q[0][2]));
    assert_eq!(chords.len(), 2, "one chord per cap: {:?}", tr.segs);
    for (chord, z) in chords.iter().zip([0.0, 2.0]) {
        for (got, want) in chord.iter().zip([[0.0, 1.0, z], [1.0, 0.0, z]]) {
            assert!(
                (0..3).all(|k| (got[k] - want[k]).abs() < 1e-12),
                "{chord:?} at z = {z}"
            );
        }
    }
}

/// **A pencil names one crossing by any pair of its classes, and a split point keeps the name it
/// came with.** The classes `x = 0`, `y = 0` and `x + y = 0` hold one line, the `z` axis; a cylinder
/// about the line `y = 0, z = 2`, `r = 1`, crosses it at `z = 1` and `z = 3`. Each pair's roots run
/// along that pair's own direction, and one pair's direction is another's reversed, so one crossing
/// wears `Lo` under one pair and `Hi` under another (the third normal at 45° between the first two
/// makes that so whichever way the class numbers fall). Every name of a crossing places it at one place
/// on the line, its pin there is the cylinder, and a split point hands the name back verbatim —
/// restated on the line's own pair with its root kept, it would be the other crossing. Once the
/// alias table learns the line, every pair's name of a crossing is one point there, and a name
/// restated on another pair of the pencil is that pair's name of the same crossing. A round that
/// ties two pairs' names learns the line from that tie alone.
#[test]
fn a_pencil_names_one_crossing_by_any_pair_and_a_split_keeps_it() {
    use crate::combinatorics::{EndPin, PointOn, QuadRoot, order_on, pierce_name, pin_for};
    use nacre_exact::Rat;
    let mut m = Model::new();
    let prism = |m: &mut Model, tri: [[f64; 2]; 3]| {
        let p2 = |q: [f64; 2]| nacre_math::Point2::from_array(q);
        let profile = crate::Profile2d::polygon(tri.into_iter().map(p2).collect()).unwrap();
        let frame = SketchFrame::world(m, Axis::Z);
        let Ok(OpOutput::Extrude { solid, .. }) = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 4.0,
            },
        ) else {
            panic!("the prism extrudes")
        };
        solid
    };
    let a = prism(&mut m, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]); // x = 0, y = 0
    let b = prism(&mut m, [[0.0, 0.0], [-1.0, 1.0], [-1.0, 0.0]]); // x + y = 0
    let bar = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([-2.0, 0.0, 2.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
        4.0,
    );
    let surf = m.face(bar.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let Some(realized) = m.surface_cache(surf).cylinder() else {
        panic!("a cylinder")
    };
    let cyls = vec![crate::planes::WorkingCyl {
        surf,
        def,
        realized: *realized,
        owner: crate::planes::SolidSide::A,
        shared: Vec::new(),
        grazed: Vec::new(),
        records: Vec::new(),
        boxes: Vec::new(),
    }];
    let mut rows = crate::planes::collect_planes(&m, a).unwrap();
    rows.extend(crate::planes::collect_planes(&m, b).unwrap());
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    let class = |n: [i128; 3]| -> usize {
        (0..planes.len())
            .find(|&c| {
                let w = jd.world_rat(c).expect("a named class");
                let k = (0..3).find(|&i| n[i] != 0).expect("a normal");
                (0..3).all(|i| w[i].checked_mul(ri(n[k])) == w[k].checked_mul(ri(n[i])))
                    && w[3] == ri(0)
            })
            .expect("the prisms hold that plane")
    };
    let (cx, cy, cd) = (class([1, 0, 0]), class([0, 1, 0]), class([1, 1, 0]));
    let at =
        |p: NodeId, q: NodeId| order_on(&jd, &cyls, cx, cy, PointOn::Pierce(p), PointOn::Pierce(q));
    let root = |n: NodeId| pierce_name(n).expect("a pierce name").2;
    // The table learns the line where a fourth class passes through it — the caps' `z = 0` — and
    // from then on each crossing is one point under every pair's name.
    let mut table = Aliases::default();
    let mut four = vec![cx, cy, cd, class([0, 0, 1])];
    four.sort_unstable();
    table.record(&jd, &four).expect("no cylinder");
    let mut reversed = false;
    for r in [QuadRoot::Lo, QuadRoot::Hi] {
        let n = NodeId::pierce(cx, cy, 0, r);
        for (p, q) in [(cx, cd), (cy, cd)] {
            let same: Vec<NodeId> = [QuadRoot::Lo, QuadRoot::Hi]
                .into_iter()
                .map(|r2| NodeId::pierce(p, q, 0, r2))
                .filter(|&n2| at(n, n2) == Some(0))
                .collect();
            let [n2] = same[..] else {
                panic!("one root of ({p}, {q}) is {n:?}'s crossing: {same:?}")
            };
            reversed |= root(n2) != root(n);
            assert!(table.one_point(n, n2), "{n:?} and {n2:?} are one crossing");
            // A round that ties the two names records the pencil as a line too — the set is the
            // three classes of the pencil alone, with no point name among them.
            // The other crossing, named by this pair: once the table knows the line, its name and
            // `n` are two points whatever their labels say, and a tie between them is refused.
            let (pair, ..) = pierce_name(n2).expect("a pierce name");
            let flip = match root(n2) {
                QuadRoot::Lo => QuadRoot::Hi,
                _ => QuadRoot::Lo,
            };
            let other = NodeId::pierce(pair[0], pair[1], 0, flip);
            if root(n2) != QuadRoot::Double {
                assert_eq!(
                    table.clone().record_same_point(&jd, &cyls, n, other),
                    Err(RejectReason::CoincidentNodes),
                    "{n:?} and {other:?} are the line's two crossings"
                );
            }
            let mut tied = Aliases::default();
            assert_eq!(tied.record_same_point(&jd, &cyls, n, n2), Ok(()));
            assert_eq!(
                tied.canon_wall(cx, cy),
                tied.canon_wall(cx, cd),
                "one line on class {cx}"
            );
            let [lo, hi] = [p.min(q), p.max(q)];
            assert_eq!(
                crate::combinatorics::restate_pierce(&jd, &cyls, n, lo, hi),
                Some(n2),
                "{n:?} restated on ({lo}, {hi})"
            );
            assert_eq!(pin_for(&jd, cx, cy, n2), Some(EndPin::Cylinder), "{n2:?}");
            let split = Split::of(n2, EndPin::Cylinder);
            assert_eq!(
                order_on(&jd, &cyls, cx, cy, split.on(cx, cy), PointOn::Pierce(n)),
                Some(0),
                "{n2:?} split on ({cx}, {cy}) is still {n:?}'s crossing"
            );
        }
    }
    assert!(
        reversed,
        "some pair runs against the line's own: a root's label alone names no crossing"
    );
}

/// **A point on a cylinder has every name its planes give it there, and only those.** The unit
/// cube's corner `(1, 0, 0)` lies on the cylinder `x² + y² = 1` about `z`: its names are the
/// three-plane one, the tangent pair `x = 1, z = 0`'s double root, and the pair `y = 0, z = 0`'s
/// root at `x = 1` (not the one at `x = −1`). The pair `x = 1, y = 0` holds a ruling of the
/// cylinder — its line lies on it — and gives no pierce name. A round that proves two of those
/// names one point records the whole set ([`Aliases::record_same_point`]). Against an arc of the
/// circle `z = 0` the point is one θ: under any of its names it is `At` an end the arc names it by
/// ([`crate::combinatorics::theta_place`]).
#[test]
fn a_point_on_a_cylinder_has_the_names_its_planes_give_it() {
    use crate::combinatorics::{QuadRoot, names_on_cylinder};
    use nacre_exact::Rat;
    use nacre_exact::quad::{CylinderMeet, plane_plane_cylinder};
    let mut m = Model::new();
    let cube = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let bar = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        3.0,
    );
    let surf = m.face(bar.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let Some(realized) = m.surface_cache(surf).cylinder() else {
        panic!("a cylinder")
    };
    let cyls = vec![crate::planes::WorkingCyl {
        surf,
        def: def.clone(),
        realized: *realized,
        owner: crate::planes::SolidSide::A,
        shared: Vec::new(),
        grazed: Vec::new(),
        records: Vec::new(),
        boxes: Vec::new(),
    }];
    let rows = crate::planes::collect_planes(&m, cube).unwrap();
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    let class = |n: [i128; 3], d: i128| -> usize {
        (0..planes.len())
            .find(|&c| {
                let w = jd.world_rat(c).expect("a named class");
                let k = (0..3).find(|&i| n[i] != 0).expect("a normal");
                (0..3).all(|i| w[i].checked_mul(ri(n[k])) == w[k].checked_mul(ri(n[i])))
                    && w[3].checked_mul(ri(n[k])) == w[k].checked_mul(ri(-d))
            })
            .expect("the cube holds that plane")
    };
    let (x1, y0, z0) = (
        class([1, 0, 0], 1),
        class([0, 1, 0], 0),
        class([0, 0, 1], 0),
    );
    let mut s = vec![x1, y0, z0];
    s.sort_unstable();
    let corner = NodeId::three_planes(Canon3::three([x1, y0, z0]));
    let names = names_on_cylinder(&jd, &cyls, 0, corner, &s);
    // The root of `y = 0, z = 0` at `x = 1`: the crossing whose x coordinate is the larger.
    let (wy, wz) = (jd.world_rat(y0).unwrap(), jd.world_rat(z0).unwrap());
    let Some(CylinderMeet::Pair { line, .. }) =
        plane_plane_cylinder(&wy, &wz, &def.origin(), &def.dir(), def.r2())
    else {
        panic!("the x axis crosses the cylinder twice")
    };
    let [lo, hi] = [y0.min(z0), y0.max(z0)];
    // `ℓ = n_y × n_z` as solved; the name's pair runs along it or against it.
    let along_plus_x = (line.dir()[0] > ri(0)) == (lo == y0);
    let root = if along_plus_x {
        QuadRoot::Hi
    } else {
        QuadRoot::Lo
    };
    let mut want = vec![
        corner,
        NodeId::pierce(x1, z0, 0, QuadRoot::Double),
        NodeId::pierce(lo, hi, 0, root),
    ];
    let mut got = names.clone();
    want.sort_unstable();
    got.sort_unstable();
    assert_eq!(got, want, "the corner's names");
    // A round that proves the corner and the axis pair's root one point records it, and the point
    // then wears every name its planes give it there — the tangent pair's `Double` included.
    let mut table = Aliases::default();
    let axis_root = NodeId::pierce(lo, hi, 0, root);
    assert_eq!(
        table.record_same_point(&jd, &cyls, corner, axis_root),
        Ok(())
    );
    assert!(table.one_point(corner, axis_root));
    assert!(table.one_point(corner, NodeId::pierce(x1, z0, 0, QuadRoot::Double)));
    // Two roots of one pair are two points: a tie between them is the arithmetic's defect. Two
    // cylinders' crossings are where two cylinders meet, which no road names.
    assert_eq!(
        Aliases::default().record_same_point(
            &jd,
            &cyls,
            NodeId::pierce(lo, hi, 0, QuadRoot::Lo),
            NodeId::pierce(lo, hi, 0, QuadRoot::Hi)
        ),
        Err(RejectReason::CoincidentNodes)
    );
    assert_eq!(
        Aliases::default().record_same_point(
            &jd,
            &cyls,
            axis_root,
            NodeId::pierce(lo, hi, 1, root)
        ),
        Err(RejectReason::CylinderPairContact)
    );
    // The corner's component stands on cylinder 0: another cylinder's name for it is refused
    // whichever of the component's names is handed in.
    assert_eq!(
        table.record_same_point(&jd, &cyls, corner, NodeId::pierce(lo, hi, 1, root)),
        Err(RejectReason::CylinderPairContact)
    );
    // A trace's tie says the same under one name: its decline is raised as its cause.
    let kind = table
        .record_same_point_traced(&jd, &cyls, corner, NodeId::pierce(lo, hi, 1, root))
        .expect_err("two cylinders");
    assert_eq!(kind, DeclineKind::CylinderPair);
    assert_eq!(
        decline_to_reject(kind, None),
        RejectReason::CylinderPairContact
    );
    // One θ order places the point against an arc of `z = 0`: the tangent pair's `Double` is the
    // axis pair's root at `θ = 0` — `At` that end, on an open arc to `(0, 1, 0)` and on the closed
    // arc alike — and the root at `x = −1` is outside the arc, `θ = 0` inside the arc from it.
    use crate::combinatorics::{ThetaPlace, theta_between, theta_place};
    let y1 = class([0, 1, 0], 1);
    let top = NodeId::pierce(y1, z0, 0, QuadRoot::Double);
    let far = NodeId::pierce(
        lo,
        hi,
        0,
        if root == QuadRoot::Hi {
            QuadRoot::Lo
        } else {
            QuadRoot::Hi
        },
    );
    let double = NodeId::pierce(x1, z0, 0, QuadRoot::Double);
    let place = |lo, hi, node| theta_place(&jd, 0, &def, lo, hi, node);
    assert_eq!(place(axis_root, top, double), Some(ThetaPlace::At(0)));
    assert_eq!(place(top, axis_root, double), Some(ThetaPlace::At(1)));
    assert_eq!(place(axis_root, axis_root, double), Some(ThetaPlace::At(0)));
    assert_eq!(
        theta_between(&jd, 0, &def, axis_root, axis_root, double),
        Some(false)
    );
    assert_eq!(theta_between(&jd, 0, &def, axis_root, top, double), None);
    assert_eq!(place(axis_root, top, far), Some(ThetaPlace::Outside));
    assert_eq!(place(far, top, double), Some(ThetaPlace::Inside));
}

/// **A pencil's restatement asks the name's own pair, not its named planes alone.** The `z` axis
/// is held by `x = 0`, `y = 0` and `x + y = 0`, and crosses the cylinder about `y = 0.6, z = 2`
/// (`r = 1`) at `z = 6/5` and `z = 14/5` — off the circle's top and bottom, where any line through a
/// crossing enters the circle the same way up and a misread sense would come out right. Boxes'
/// floors and tops turned ±20° about the `x` direction through each crossing are planes with no
/// world statement through it, cutting the axis there (not ±37°: the 3-4-5 turn makes one of them
/// tangent to the circle at the crossing). The crossing's name by `x = 0` and that floor, restated
/// on another pair of the pencil, is that pair's name of the crossing: the named plane `x = 0`
/// holds the axis, so both of the pair's roots lie on it, and the floor — the plane that cuts — is
/// read by the intervals. Read as a pencil (the name's pair
/// restated by direction), the floor's pair is no pair of the axis and the sign means nothing.
#[test]
fn a_restatement_reads_a_cutting_plane_with_no_world_statement() {
    use crate::combinatorics::{PointOn, QuadRoot, order_on, pierce_of_pair};
    use nacre_exact::{Angle, Axis as EAxis, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let prism = |m: &mut Model, tri: [[f64; 2]; 3]| {
        let p2 = |q: [f64; 2]| nacre_math::Point2::from_array(q);
        let profile = crate::Profile2d::polygon(tri.into_iter().map(p2).collect()).unwrap();
        let frame = SketchFrame::world(m, Axis::Z);
        let Ok(OpOutput::Extrude { solid, .. }) = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 4.0,
            },
        ) else {
            panic!("the prism extrudes")
        };
        solid
    };
    let a = prism(&mut m, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let b = prism(&mut m, [[0.0, 0.0], [-1.0, 1.0], [-1.0, 0.0]]);
    // `below`: the box's top face, not its floor, passes through the crossing — the plane's stored
    // normal then faces the other way.
    // `z` in fifths: the crossings are `z = 6/5` and `z = 14/5`.
    let turned = |m: &mut Model, z: i128, deg: i128, below: bool| {
        let z0 = z - if below { 5 } else { 0 };
        let floor = crate::fixtures::cuboid(
            m,
            [-1.0, -1.0, z0 as f64 / 5.0],
            [1.0, 1.0, (z0 + 5) as f64 / 5.0],
        );
        crate::transform::transform(
            m,
            floor,
            &Isometry::rotation(Rotation {
                axis: EAxis::X,
                pivot: [
                    Rat::from_int(0),
                    Rat::from_int(0),
                    Rat::new(z, 5).expect("a fifth"),
                ],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
            }),
        )
        .expect("the box turns")
    };
    // Both crossings, both turns, both facings: the sense a pencil reading would take from the
    // floor's pair is the floor's, so some of these read it the wrong way.
    let floors: Vec<(i128, i128, Handle<Solid>)> = [6, 14]
        .into_iter()
        .flat_map(|z| {
            [
                (z, 20, false),
                (z, 20, true),
                (z, -20, false),
                (z, -20, true),
            ]
        })
        .map(|(z, deg, below)| (z, deg, turned(&mut m, z, deg, below)))
        .collect();
    let bar = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([-2.0, 0.6, 2.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
        4.0,
    );
    let surf = m.face(bar.lateral).surface;
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let Some(realized) = m.surface_cache(surf).cylinder() else {
        panic!("a cylinder")
    };
    let cyls = vec![crate::planes::WorkingCyl {
        surf,
        def: def.clone(),
        realized: *realized,
        owner: crate::planes::SolidSide::A,
        shared: Vec::new(),
        grazed: Vec::new(),
        records: Vec::new(),
        boxes: Vec::new(),
    }];
    let mut rows = crate::planes::collect_planes(&m, a).unwrap();
    rows.extend(crate::planes::collect_planes(&m, b).unwrap());
    for &(_, _, f) in &floors {
        rows.extend(crate::planes::collect_planes(&m, f).unwrap());
    }
    let canon = crate::planes::plane_classes(&crate::planes::test_judge(&rows));
    let (planes, _, _) = crate::planes::dense_planes(&rows, &canon);
    let jd = crate::planes::test_judge(&planes);
    let ri = Rat::from_int;
    let class = |n: [i128; 3]| -> usize {
        (0..planes.len())
            .find(|&c| {
                jd.world_rat(c).is_some_and(|w| {
                    let k = (0..3).find(|&i| n[i] != 0).expect("a normal");
                    (0..3).all(|i| w[i].checked_mul(ri(n[k])) == w[k].checked_mul(ri(n[i])))
                        && w[3] == ri(0)
                })
            })
            .expect("the prisms hold that plane")
    };
    let (cx, cy, cd) = (class([1, 0, 0]), class([0, 1, 0]), class([1, 1, 0]));
    let apart = |p: NodeId, q: NodeId| -> bool {
        let (pp, qq) = (
            crate::combinatorics::pierce_name(p).unwrap(),
            crate::combinatorics::pierce_name(q).unwrap(),
        );
        let (Some(x), Some(y)) = (
            jd.pierce(pp.0, 0, &def, pp.2),
            jd.pierce(qq.0, 0, &def, qq.2),
        ) else {
            panic!("both points are placed")
        };
        jd.points_apart(
            nacre_judge::cylinder::CylPoint::Pierce(&x),
            nacre_judge::cylinder::CylPoint::Pierce(&y),
        )
    };
    for &(z, deg, _) in &floors {
        // The turned floor: no world statement, through `(0, 0, z)` with normal
        // `(0, −sin θ, cos θ)`.
        let (s, c) = (
            (deg as f64).to_radians().sin(),
            (deg as f64).to_radians().cos(),
        );
        let tf = (0..planes.len())
            .find(|&k| {
                let n = planes[k].plane.normal().as_array();
                jd.world_rat(k).is_none()
                    && ((-n[1] * s + n[2] * c).abs() - 1.0).abs() < 1e-9
                    && planes[k]
                        .plane
                        .signed_distance(Point3::from_array([0.0, 0.0, z as f64 / 5.0]))
                        .abs()
                        < 1e-9
            })
            .expect("the turned floor is a class with no world statement");
        // The crossing on the axis by the floor and `x = 0`, and by the axis's own pair.
        let at = [QuadRoot::Lo, QuadRoot::Hi]
            .map(|r| NodeId::pierce(cx, tf, 0, r))
            .into_iter()
            .find(|&n| {
                [QuadRoot::Lo, QuadRoot::Hi]
                    .into_iter()
                    .any(|r| !apart(n, NodeId::pierce(cx, cy, 0, r)))
            })
            .expect("one of the floor pair's crossings is on the axis");
        let on_axis = [QuadRoot::Lo, QuadRoot::Hi]
            .map(|r| NodeId::pierce(cx, cy, 0, r))
            .into_iter()
            .find(|&n| !apart(at, n))
            .expect("the axis's crossing at the floor");
        for (p, q) in [(cy, cd), (cx, cd)] {
            let want = [QuadRoot::Lo, QuadRoot::Hi]
                .map(|r| NodeId::pierce(p, q, 0, r))
                .into_iter()
                .find(|&n2| {
                    order_on(
                        &jd,
                        &cyls,
                        cx,
                        cy,
                        PointOn::Pierce(on_axis),
                        PointOn::Pierce(n2),
                    ) == Some(0)
                })
                .expect("the pair's name of the crossing");
            assert_eq!(
                pierce_of_pair(&jd, &cyls, 0, at, p.min(q), p.max(q)),
                Some(want),
                "{at:?} restated on ({p}, {q})"
            );
        }
    }
}
