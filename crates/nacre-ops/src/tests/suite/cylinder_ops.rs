//! Cylinder operations end to end: axis and seam integers, bores, caps, windings, the mesh census.

use super::*;

/// ★★ The node-omission normalization's ground: **a seeded world plane's canonical frame
/// realizes bit-exactly on the unit axes** — origin `(0,0,0)`, every axis component `0.0` or
/// `±1.0`, no rounding anywhere. Such a frame's axes lift to exact orthonormal rationals, so
/// the operation's `exact()` gate elides the node (rule 207's normalization) and loses
/// nothing: the frame the node would state is the world statement already made.
///
/// The expected axes are the **arbitrary-axis convention's**, not `axis_plane`'s script
/// triples: for the ZX plane the convention's `+u` is `ẑ × ŷ = −x̂` where the script's is
/// `+ẑ` — the documented case a caller states through a `Named` placement instead. The seed
/// *points* carry the script triple; the *canonical frame* is the plane's own.
#[test]
fn a_seeded_planes_canonical_frame_is_the_world_basis_exactly() {
    use nacre_scalar::Axis;
    let m = Model::new();
    // (axis, expected û, v̂, ŵ) — ŵ is the +axis (the canonical name's own sign).
    let want = [
        (Axis::Z, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        (Axis::X, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        (Axis::Y, [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
    ];
    for (axis, u, v, w) in want {
        let h = m.world_plane(axis);
        let (o, ru, rv, rw) = crate::rotated_vertex::frame_world_basis(
            &m,
            h,
            &nacre_topo::FramePlacement::Canonical,
            false,
        )
        .expect("a seed always carries a name");
        assert_eq!(
            (o, ru, rv, rw),
            ([0.0; 3], u, v, w),
            "{axis:?}: not the world basis"
        );
    }
}

/// ★★ **The frame keeps the statement rational.** On world XY the axis is exactly `(0,0,1)` and
/// the seam reference exactly `(1,0,0)` — small integers, not decimals lifted back out of a
/// normalization. The centre and radius are the caller's own written decimals, lifted once.
///
/// This is the whole reason the op exists as it does: `add_cylinder` reaches the same shape by
/// normalizing an f64 axis, and the `1/3`-base gate in `nacre-topo` shows what that costs.
#[test]
fn a_cylinder_op_states_its_axis_and_seam_as_integers() {
    let mut m = Model::new();
    let op = cylinder_op(&m, [0.5, 0.25], 0.1, 2.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).expect("a world-XY cylinder builds")
    else {
        panic!("an extrude answers with an extrude");
    };
    let def = lone_cylinder_def(&m);
    let int = |n: i128| nacre_scalar::Rat::from_int(n);
    let rat = |n: i128, d: i128| nacre_scalar::Rat::new(n, d).unwrap();
    assert_eq!(
        def.dir(),
        [int(0), int(0), int(1)],
        "the frame's unit normal"
    );
    assert_eq!(def.ref_dir(), [int(1), int(0), int(0)], "the frame's +u");
    assert_eq!(
        def.origin(),
        [rat(1, 2), rat(1, 4), int(0)],
        "the centre, placed"
    );
    assert_eq!(
        *def.r2(),
        nacre_scalar::BigRat::from(rat(1, 100)),
        "r² of the stated 1/10"
    );
    assert_eq!(def.radius_exact(), Some(rat(1, 10)));
    assert_eq!(faces.len(), 3, "bottom cap, top cap, lateral");
}

/// A cylinder in a log replays to the same model — handles and coordinates both. The op is only
/// worth having if the log stays reproducible with it in there.
#[test]
fn a_cylinder_op_replays_deterministically() {
    let log = vec![
        extrude_log_op(square(), 1.0),
        cylinder_op(&Model::new(), [0.5, 0.5], 0.2, 1.0),
    ];
    let (m1, m2) = (replay(&log).unwrap(), replay(&log).unwrap());
    let pts = |m: &Model| {
        (0..m.vertex_count() as u32)
            .filter_map(|i| m.vertex_handle_at(i))
            .map(|h| (h, m.vertex(h)))
            .map(|(vh, _)| m.vertex_point(vh).as_array())
            .collect::<Vec<_>>()
    };
    assert_eq!(pts(&m1), pts(&m2));
    assert_eq!(m1.face_count(), m2.face_count());
    assert_eq!(m1.surface_count(), m2.surface_count());
}

/// ★★ **Every refusal is a name, and none of them leaves a cell behind.** The circle's own
/// numbers are refused where they are written — at the sketch, before any operation exists — and
/// the operation-level refusals are decided before anything is pushed: a refusal that had already
/// pushed would shift every later log index (`apply`'s own doc says so).
#[test]
fn a_refused_cylinder_op_is_named_and_leaves_nothing_behind() {
    let mut m = Model::new();
    // A cylinder surface to aim a frame at — the one thing a `SketchFrame` can name that is not
    // a plane.
    let fixture = cylinder_op(&m, [0.0, 0.0], 1.0, 1.0);
    let OpOutput::Extrude { faces, .. } =
        apply(&mut m, &fixture).expect("the fixture cylinder builds")
    else {
        panic!("an extrude answers with an extrude");
    };
    let on_lateral = SketchFrame::canonical(m.face(faces[2]).surface);

    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let wide = 1e300;
    for (center, radius) in [(p2(0.0, 0.0), 0.0), (p2(0.0, 0.0), -1.0)] {
        assert!(
            matches!(
                crate::Ring2d::circle(center, radius),
                Err(crate::SketchError::NonPositiveRadius { .. })
            ),
            "radius {radius}"
        );
    }
    for (center, radius) in [(p2(0.0, 0.0), wide), (p2(wide, 0.0), 1.0)] {
        assert!(
            matches!(
                crate::Ring2d::circle(center, radius),
                Err(crate::SketchError::OutsideDecimalWindow { .. })
            ),
            "a number outside the decimal window"
        );
    }

    let before = (
        m.surface_count(),
        m.vertex_count(),
        m.edge_count(),
        m.face_count(),
        m.live_solids().len(),
    );
    let cases: [(Operation, OpError); 2] = [
        (
            cylinder_op(&m, [0.0, 0.0], 1.0, 0.0),
            OpError::NonPositiveDistance,
        ),
        (
            Operation::Extrude {
                frame: on_lateral,
                profile: circle_profile([0.0, 0.0], 1.0),
                dist: 1.0,
            },
            OpError::NonPlanarFace,
        ),
    ];
    for (op, want) in cases {
        assert_eq!(apply(&mut m, &op).unwrap_err(), want, "op: {op:?}");
    }
    assert_eq!(
        (
            m.surface_count(),
            m.vertex_count(),
            m.edge_count(),
            m.face_count(),
            m.live_solids().len()
        ),
        before,
        "a refusal is decided before anything is pushed"
    );
}

/// ★★★ **What the app will do, done through the log alone**: a plate, a drill standing through
/// it, and a `Cut`.
///
/// The drill **overshoots** the plate, the way a "through all" hole is drawn. A drill that stops
/// exactly on the plate's own faces works too (`bands`' seated-cap block measures that family);
/// this fixture is the overshooting gesture a kit step emits.
#[test]
fn a_logged_cylinder_cuts_a_through_hole() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    // The log is assembled against a scratch model built the same way, so its handles are the
    // index vocabulary `replay` re-anchors (the `two_extrudes_make_two_solids` precedent).
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Extrude {
        frame,
        profile: circle_profile([2.0, 2.0], 0.5),
        dist: 3.0,
    };
    let OpOutput::Extrude { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let log = vec![
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ];
    let m = replay(&log).expect("a logged drill");
    assert_eq!(m.live_solids().len(), 1, "one drilled plate");
    let v = nacre_props::mass_props(&m, m.live_solids()[0])
        .expect("props")
        .volume;
    let want = 4.0 * 4.0 * 2.0 - std::f64::consts::PI * 0.25 * 2.0;
    assert!(
        (v - want).abs() < 1e-9,
        "volume {v} is not the plate minus the bore {want}"
    );
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
}

/// ★★ **The negative control the hole above needs**: the same drill, moved clear of the plate,
/// removes nothing. Without it, "the volume dropped by πr²t" could be read as "any cylinder in
/// the log makes a hole".
#[test]
fn a_cylinder_clear_of_the_plate_removes_nothing() {
    let mut m = Model::new();
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: plate_h, .. } = apply(&mut m, &extrude).expect("the plate")
    else {
        panic!("an extrude answers with an extrude");
    };
    let before = nacre_props::mass_props(&m, plate_h).expect("props").volume;
    let above = datum_frame(
        &mut m,
        crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 3.0])),
    );
    let clear = Operation::Extrude {
        frame: above,
        profile: circle_profile([2.0, 2.0], 0.5),
        dist: 1.0,
    };
    let OpOutput::Extrude { solid: drill, .. } = apply(&mut m, &clear).expect("a clear cylinder")
    else {
        panic!("an extrude answers with an extrude");
    };
    let solids = boolean(&mut m, BoolKind::Cut, plate_h, drill).expect("a disjoint cut answers");
    assert_eq!(solids.len(), 1, "cutting away nothing leaves one solid");
    let after = nacre_props::mass_props(&m, solids[0])
        .expect("props")
        .volume;
    assert!(
        (after - before).abs() < 1e-9,
        "a drill clear of the plate removed {} of material",
        before - after
    );
}

/// ★★ **A face's frame faces *out* of its solid, so a cylinder on it stands on top rather than
/// drilling in.** Asked as a placement question, which is where it is decided — the seam
/// vertices sit at the plate's top face and one unit above it, not below.
///
/// This is the honest statement of what the log expresses today: `flip` is *measured* by the
/// consuming operation (`measured_frame`), never stated by a caller, so drilling into a face is a
/// composite that measures `toward` inward — `PocketOnFace`'s sibling, and not yet written.
#[test]
fn a_cylinder_on_a_face_frame_stands_outward() {
    let mut m = Model::new();
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude).expect("the plate builds") else {
        panic!("an extrude answers with an extrude");
    };
    // faces[1] is the top cap; its frame is measured against that face's outward normal.
    let top = face_sketch_frame(&m, faces[1]).expect("a face has a frame");
    let before = m.vertex_count();
    let boss = Operation::Extrude {
        frame: top,
        profile: circle_profile([2.0, 2.0], 0.5),
        dist: 1.0,
    };
    apply(&mut m, &boss).expect("a cylinder on a face frame builds");
    let seam_z: Vec<f64> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .map(|h| (h, m.vertex(h)))
        .skip(before)
        .filter(|(_, v)| matches!(**v, Vertex::OnSeam(_)))
        .map(|(h, _)| m.vertex_point(h).as_array()[2])
        .collect();
    assert_eq!(seam_z.len(), 2, "one seam vertex per rim");
    assert!(
        seam_z.iter().all(|z| *z >= 2.0 - 1e-12),
        "the cylinder stands outward from the face, not into the plate: {seam_z:?}"
    );
    assert!(
        seam_z.iter().any(|z| (*z - 3.0).abs() < 1e-12),
        "and reaches its full height above it: {seam_z:?}"
    );
}

/// ★★ **A tilted frame builds a cylinder and cannot yet cut with it — and both halves are the
/// point.** The frame has no exact rational basis, so the statement is written in the plane's own
/// frame and a motion node carries it out (the road a tilted prism already takes). The gate then
/// declines a *moved* cylinder by name, because its truth is stated before its motion.
///
/// If the first half ever fails, the op stopped taking the frame-node road; if the second starts
/// succeeding, moved cylinders are supported and this test should be re-read, not deleted.
#[test]
fn a_cylinder_on_a_tilted_frame_is_built_and_honestly_declined() {
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: square(),
        dist: 1.0,
    };
    let OpOutput::Extrude { solid: cube_h, .. } = apply(&mut m, &extrude).expect("the cube builds")
    else {
        panic!("an extrude answers with an extrude");
    };
    let tilted = datum_frame(
        &mut m,
        crate::SketchPlane::from_origin_normal(
            Point3::from_array([0.5, 0.5, -1.0]),
            nacre_math::Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
        )
        .expect("a tilted plane"),
    );
    let op = Operation::Extrude {
        frame: tilted,
        profile: circle_profile([0.0, 0.0], 0.2),
        dist: 4.0,
    };
    let OpOutput::Extrude {
        solid: drill,
        faces,
    } = apply(&mut m, &op).expect("a tilted cylinder")
    else {
        panic!("an extrude answers with an extrude");
    };
    let lateral = m.face(faces[2]).surface;
    match m.surface(lateral) {
        nacre_topo::Surface::Cylinder { motion, .. } => assert!(
            motion.is_some(),
            "a tilted frame states its cylinder inside a motion node"
        ),
        _ => panic!("the third face is the lateral"),
    }
    // ★ The caches are realized from the same *local* rationals the truth is stated in — the
    // prism road's deal (`SweptRat::motion`). `validate` reads truth against cache, so a
    // world/local mix-up here would surface as `CylinderTruthCacheMismatch` rather than as a
    // wrong model much later.
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, cube_h, drill),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            })
        ),
        "a moved cylinder is declined by name, not silently mis-cut"
    );
}

/// ★★★ **A cylinder's two caps face opposite ways — even when their planes already exist.**
///
/// The bottom cap lies on the very plane the frame names, so it always interns; a plane's
/// canonical name has no direction, so the surface handed back can be the one that faces the
/// other way. Measured before the fix: a cylinder on a plate's top face came out with **both
/// caps facing +Z**, because the b-rep body dropped `push_plane`'s `flipped` bit (the bit
/// `add_cuboid` has always honoured). Not a solid at all, and no earlier fixture looked.
///
/// The three cases are the three ways a frame's plane can arrive: a seeded world plane, one a
/// `DatumPlane` stated, and one an earlier face already put there.
#[test]
fn a_cylinders_caps_face_opposite_ways_however_their_planes_arrived() {
    let outward = |m: &Model, fh: Handle<Face>| -> [f64; 3] {
        let f = m.face(fh);
        let nacre_geom::Surface::Plane(p) = m.surface_cache(f.surface) else {
            panic!("a cap is planar")
        };
        let s = f.orientation.sign() as f64;
        p.normal().as_array().map(|c| c * s)
    };
    let check = |m: &Model, faces: &[Handle<Face>], what: &str| {
        let (bot, top) = (outward(m, faces[0]), outward(m, faces[1]));
        assert!(
            bot[2] < 0.0,
            "{what}: the bottom cap must face down, got {bot:?}"
        );
        assert!(
            top[2] > 0.0,
            "{what}: the top cap must face up, got {top:?}"
        );
        assert!(
            nacre_validate::validate(m).is_empty(),
            "{what}: {:?}",
            nacre_validate::validate(m)
        );
    };

    // (a) a seeded world plane
    let mut m = Model::new();
    let op = cylinder_op(&m, [0.0, 0.0], 1.0, 2.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    check(&m, &faces, "world XY");

    // (b) a plane the log stated as a datum, facing +Z
    let mut m = Model::new();
    let below = datum_frame(
        &mut m,
        crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
    );
    let op = Operation::Extrude {
        frame: below,
        profile: circle_profile([0.0, 0.0], 1.0),
        dist: 2.0,
    };
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    check(&m, &faces, "stated datum");

    // (c) the plane of a face that already exists — the case that was wrong
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: square(),
        dist: 1.0,
    };
    let OpOutput::Extrude { faces: plate, .. } = apply(&mut m, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let top = face_sketch_frame(&m, plate[1]).expect("a face has a frame");
    let op = Operation::Extrude {
        frame: top,
        profile: circle_profile([0.5, 0.5], 0.2),
        dist: 1.0,
    };
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    check(&m, &faces, "an existing face's plane");
}

/// ★★ **A hole must wind the other way round, and now something checks it.**
///
/// `build_prism` is the single place that decides winding — outer CCW about the
/// sweep, holes the opposite — and until now nothing verified that decision:
/// `shell_signed_volume` takes an *unsigned* area and subtracts, so it assumes the
/// convention rather than reading it.
///
/// The defect is built by hand, and the twin's **orientation** is what gets flipped
/// rather than the loop's direction: reversing a rim half-edge would make
/// `NonOpposedEdge` fire first (the rim is shared with the cylinder's wall) and this
/// check's specificity would be lost. Flipping the flag leaves edge traversal alone,
/// so the two violations that come back are both this rule's — and **`cos`'s sign
/// says which loop spoke**: `−1` the outer polygon, `+1` the inner rim.
#[test]
fn a_holes_winding_is_checked_too() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Extrude {
        frame,
        profile: circle_profile([2.0, 2.0], 0.5),
        dist: 3.0,
    };
    let OpOutput::Extrude { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let mut m = replay(&[
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ])
    .expect("a logged drill");
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the drilled plate is clean to begin with: {:?}",
        nacre_validate::validate(&m)
    );

    // A drilled face: an outer polygon with one rim hole.
    let solid = m.live_solids()[0];
    let shell = m.solid(solid).outer;
    let faces = m.shell(shell).faces.clone();
    let victim = *faces
        .iter()
        .find(|&&fh| {
            let f = m.face(fh);
            f.inner.len() == 1 && f.inner[0].half_edges.len() == 1
        })
        .expect("a through hole leaves two drilled faces");
    let twin = {
        let f = m.face(victim).clone();
        m.push_face(nacre_topo::Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.push_shell(nacre_topo::Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(nacre_topo::Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![replaced]);
    m.rebuild_adjacency();

    let vs = nacre_validate::validate(&m);
    let cosines: Vec<f64> = vs
        .iter()
        .filter_map(|v| match v {
            nacre_validate::Violation::FaceMisoriented { face, cos } if *face == twin => Some(*cos),
            _ => None,
        })
        .collect();
    assert_eq!(vs.len(), 2, "both loops of that face should speak: {vs:?}");
    assert!(
        cosines.iter().any(|c| *c < -0.5),
        "the outer polygon must report a flipped flag: {cosines:?}"
    );
    assert!(
        cosines.iter().any(|c| *c > 0.5),
        "the hole must report that it now winds like an outer loop — this is the \
         assertion that says the inner branch ran: {cosines:?}"
    );
}

/// The polygonal twin of `a_holes_winding_is_checked_too`: a prism with a square
/// hole. Same rule, same convention, and — measured across the suite — a real
/// population (444 polygonal hole loops, all of them already honouring it).
///
/// Worth its own fixture because the two shapes reach the rule by different roads:
/// a rim answers from its circle, a polygon from its Newell sum, and only running
/// both says the shared `loop_winding` serves both.
#[test]
fn a_polygonal_holes_winding_is_checked_too() {
    let profile = Profile2d::with_holes(
        vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)],
        vec![vec![p2(1.0, 1.0), p2(3.0, 1.0), p2(3.0, 3.0), p2(1.0, 3.0)]],
    )
    .expect("a square ring");
    let mut m = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile,
        dist: 1.0,
    };
    let OpOutput::Extrude { solid, .. } = apply(&mut m, &extrude).expect("a ring prism") else {
        panic!("an extrude answers with an extrude");
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the ring prism is clean to begin with: {:?}",
        nacre_validate::validate(&m)
    );

    let shell = m.solid(solid).outer;
    let faces = m.shell(shell).faces.clone();
    let victim = *faces
        .iter()
        .find(|&&fh| {
            let f = m.face(fh);
            f.inner.len() == 1 && f.inner[0].half_edges.len() >= 3
        })
        .expect("a ring prism has two holed caps");
    let twin = {
        let f = m.face(victim).clone();
        m.push_face(nacre_topo::Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.push_shell(nacre_topo::Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(nacre_topo::Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![replaced]);
    m.rebuild_adjacency();

    let vs = nacre_validate::validate(&m);
    let cosines: Vec<f64> = vs
        .iter()
        .filter_map(|v| match v {
            nacre_validate::Violation::FaceMisoriented { face, cos } if *face == twin => Some(*cos),
            _ => None,
        })
        .collect();
    assert_eq!(vs.len(), 2, "both loops of that face should speak: {vs:?}");
    assert!(
        cosines.iter().any(|c| *c < -0.5) && cosines.iter().any(|c| *c > 0.5),
        "one report per loop, told apart by sign: {cosines:?}"
    );
}

/// ★★ **A boss's wall faces away from its axis; a bore's faces toward it.**
///
/// Same geometry, opposite sense — the material is inside the wall for one and outside it for
/// the other — and the boolean says so in its own words when it mints a curved face: *"a
/// cylinder's stored normal points away from its axis, and that is the face's outward normal
/// for a boss … `flip` (material outside the wall: a hole) is what reverses it."*
///
/// This asks that sentence back from a different layer: `props::face_normal_at` reads the
/// stored `orientation` and turns the radial direction by it. If the two ever disagreed, a
/// viewer would light a bore inside out and nothing else would notice.
#[test]
fn a_bores_wall_faces_its_axis_and_a_bosss_faces_away() {
    let plate =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 4.0), p2(0.0, 4.0)]).unwrap();
    let mut scratch = Model::new();
    let extrude = Operation::Extrude {
        frame: SketchFrame::world(&scratch, Axis::Z),
        profile: plate,
        dist: 2.0,
    };
    let OpOutput::Extrude { solid: a, .. } = apply(&mut scratch, &extrude).unwrap() else {
        panic!("an extrude answers with an extrude");
    };
    let below = Operation::DatumPlane {
        def: DatumDef::Stated(
            crate::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, -0.5])),
        ),
    };
    let OpOutput::DatumPlane { frame, .. } = apply(&mut scratch, &below).unwrap() else {
        panic!("a datum answers with a datum");
    };
    let drill = Operation::Extrude {
        frame,
        profile: circle_profile([2.0, 2.0], 0.5),
        dist: 3.0,
    };
    let OpOutput::Extrude { solid: b, .. } = apply(&mut scratch, &drill).unwrap() else {
        panic!("an extrude answers with an extrude");
    };

    // The axis of both the bore and the free-standing drill, as a line to measure against.
    let axis_at = |z: f64| Point3::from_array([2.0, 2.0, z]);
    let radial_sense = |m: &Model, solid: Handle<Solid>| -> f64 {
        let shell = m.solid(solid).outer;
        let wall = *m
            .shell(shell)
            .faces
            .iter()
            .find(|&&fh| {
                matches!(
                    m.surface_cache(m.face(fh).surface),
                    nacre_geom::Surface::Cylinder(_)
                )
            })
            .expect("a cylindrical wall");
        // A point on that wall: the seam vertex of one of its rims.
        let he = m.face(wall).outer.half_edges[0];
        let p = m.vertex_point(m.edge(he.edge).vertices[0]);
        let n = nacre_props::face_normal_at(m, wall, p).expect("off the axis");
        let out = p - axis_at(p.as_array()[2]);
        n.dot(out)
    };

    // The drill on its own is a boss: its wall faces outward.
    assert!(
        radial_sense(&scratch, b) > 0.0,
        "a free-standing cylinder's wall faces away from its axis"
    );

    // The same cylinder, once it is a hole, faces the other way.
    let m = replay(&[
        extrude,
        below,
        drill,
        Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    ])
    .expect("a logged drill");
    assert!(
        radial_sense(&m, m.live_solids()[0]) < 0.0,
        "a bore's wall faces its axis — the material is outside the wall"
    );
}

/// ★★ **A tunnel crossing a bore is not a tunnel touching a bore** — the model a user brought,
/// and the pair rule's real question.
///
/// A plate, a vertical bore through it, and a horizontal tunnel six away with radii summing to
/// four. Nothing touches, and the volume says so exactly. It used to be refused as "two
/// cylinders touch or overlap" because the gate could only measure the distance between
/// *parallel* axes and read its own blind spot as contact.
///
/// ★ The pair rule fires on the **second** boolean — the first result already carries a
/// cylindrical face — so a single cut would not exercise it at all.
///
/// The second half is what keeps the first honest: slide the tunnel onto the bore and the
/// refusal must come back, because two cylinders that really do meet cross in a quartic curve
/// nothing here computes. (Measured with the guard removed: the result kept the *same* volume
/// as the disjoint case — a silent wrong answer, validate clean.)
#[test]
fn a_tunnel_clear_of_a_bore_is_cut_and_one_through_it_is_refused() {
    let build = |x: f64| -> (Model, Handle<Solid>, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([-10.0, -10.0, -1.5]),
            Point3::from_array([10.0, 10.0, 1.5]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([0.0, 0.0, -10.0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            20.0,
        );
        let tunnel = m.add_cylinder(
            Point3::from_array([x, -25.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
            1.0,
            50.0,
        );
        m.rebuild_adjacency();
        (m, plate, bore, tunnel)
    };

    // Six apart: the two cuts go through, and the volume is the plate minus both voids.
    let (mut m, plate, bore, tunnel) = build(-6.0);
    let drilled = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    m.rebuild_adjacency();
    let both = boolean(&mut m, BoolKind::Cut, drilled, tunnel).expect("the tunnel cuts too")[0];
    let vol = nacre_props::mass_props(&m, both).expect("props").volume;
    let want = 20.0 * 20.0 * 3.0 - std::f64::consts::PI * 9.0 * 3.0 - std::f64::consts::PI * 20.0;
    assert!(
        (vol - want).abs() < 1e-9,
        "plate minus a bore and a tunnel: {vol} vs {want}"
    );
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );

    // Straight through the bore: still refused, and by the name that now states a fact.
    let (mut m, plate, bore, tunnel) = build(0.0);
    let drilled = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
    m.rebuild_adjacency();
    assert!(
        matches!(
            boolean(&mut m, BoolKind::Cut, drilled, tunnel),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderPairContact,
                ..
            })
        ),
        "a tunnel through the bore is a cylinder pair that meets"
    );
}

/// **The mesh census is live** — the invariant itself is asserted where the fact is made.
///
/// `boolean::tess_census::record` checks every result as it is produced ("every solid a boolean
/// returns meshes" — with no exemption), because a test reading the vector
/// afterwards sees only the booleans that ran before it. What is left for a test is that the
/// census is **running at all**: a hook that silently stopped recording would take the whole
/// guarantee with it and nothing would go red.
///
/// ★ The count is deliberately not asserted. It is whatever the suite happened to run before this
/// test, and it moves with every fixture added. The formerly refused population (the two tangency
/// fixtures in `bands::tests`) is now held to its faces' exact areas there, through
/// [`mesh_covers_faces`].
#[test]
fn the_mesh_census_is_running() {
    let seen = crate::boolean::tess_census::MESHED
        .lock()
        .expect("the census lock is never held across a panic")
        .clone();
    assert!(!seen.is_empty(), "the census never ran");
    assert!(
        seen.iter().any(|r| r.is_ok()),
        "the census recorded no mesh at all"
    );
}

/// **A ruling carries the label of the cell inside the cylinder** — capability D's vertical answer.
///
/// The chart's horizontal lines have had their answer since the band pass: a ⊥ class's
/// [`crate::arrangement::Label`] holds the material on **both** sides of its plane, which is why
/// the cell reader reads two disk labels and asks them to agree. The vertical lines had none, and
/// this is it — the label of the cell a ruling borders on the axis side of its wall.
///
/// ★★★★★ **The side is derived, and an independent description checks it.**
/// `ruling_interior_is_even` composes three sentences already in the crate (`RulingCarrier::side`,
/// the material-on-the-left convention in the **root face's frame** — `frame_sign` — and
/// `world_rat_sense`'s lift between the rational name and the stored normal) into
/// `side · κ · frame_sign`. The check shares no step with that: **inside the cylinder
/// the lateral's own solid has material and outside it does not**, which is the same content rule
/// `ArcLabels`' doc set its own side by. It is asserted at the record, in `per_class`.
///
/// ★★★★ **And the check is what watches this sign — the volume oracle cannot.** Everywhere else in
/// this ladder a side selector is guarded by the through-boss volume (a global flip passes the
/// relative locks), but that only works for a sign production *reads*. This label is
/// `#[cfg(test)]`, so no volume moves whatever it says. ☑ Flipping the product turns the check from
/// every agreement into a contradiction — it has eyes on 97% of the population.
///
/// ☑ Measured over the whole binary: **862 ruling pieces, 862 labelled** (`world_rat_sense`
/// declines none) · **838 agree, 0 contradict, 24 blind**. The blind ones are real and expected —
/// a boss surrounded by its own plate has that solid's material on *both* sides of the ruling, so
/// the content cannot tell, and only the derivation speaks there.
#[test]
fn a_ruling_labels_the_cell_inside_the_cylinder() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds");
    let lab = crate::arrangement::ruling_probe::LABELLED
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    let chk = crate::arrangement::ruling_probe::SIDE_CHECK
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    assert!(
        !lab.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    // ★ Universal over every recorded piece, which no interleaving can break — a filtered or
    // parallel run may bring more entries, never different ones.
    assert!(
        lab.iter().all(|&b| b),
        "a ruling piece got no label: {} of {} unlabelled",
        lab.iter().filter(|b| !**b).count(),
        lab.len()
    );
    assert!(
        !chk.is_empty(),
        "the side check never ran beside the labels"
    );
    assert!(
        chk.iter().all(|c| *c != Some(false)),
        "a ruling label contradicted the content check"
    );
    // ★ And the check is not vacuous: it decides for most of the population, not none of it.
    assert!(
        chk.contains(&Some(true)),
        "the content check decided nothing at all"
    );
}

/// **The census's cylinder corpus, in the lib suite.** `tests/census.rs` records these
/// families' digests, but it is an integration test and the lib is compiled without `cfg(test)`
/// there — so the chart's census, which asserts where the facts are made, never sees them there
/// (`wal corner-lo` slipped past the shadow comparison that way). What these
/// lock is only the **outcome** (built, or the refusal's name); the geometry stays the digest's.
/// Fixtures copied from `tests/census.rs` (`rul`, `wal`, `cap`, `ct2`, `trc` families).
///
/// ★★★★★ **`wal corner-lo` used to refuse here, and the reason written at this line was the
/// symptom.** It read: *"the plate classes' disk cells carry no B material at the (0,0) corner
/// while the caps do (an arrangement label defect) — the cells' ends disagree and the emitter
/// refuses rather than read either."* The missing material was real, but it was not a labelling
/// defect: that corner puts the boss's axis on **both** plate walls, so the chords through its
/// circle are radii and the sector outside the plate is **reflex** at the centre. `loop_winding`
/// read the turn at the lexicographically smallest **node** — the centre — which the arc bulges
/// past, so the sign came back inverted, the void seed landed on a bounded sector, and one
/// solid's material was flipped across the whole component. Reading the winding at the ring's own
/// extremum (`combinatorics::arc_extremum_winding`) fixes it, and all three kinds build with
/// their exact volumes.
///
/// ★ `rul flush`: a straddling boss whose caps sit flush with **both** plate planes. Cut
/// and Common assemble clean and build; the Fuse's two cap chords are each an interior boundary
/// between a half-disk 2-gon and its neighbour that the coplanar merge abstains on (the chord and
/// its arcs collide in the node-pair edge key), and the shipped pair spelled a stated plane-self
/// edge — validate's producer bug. The whole-result check refuses that shape by the cleaning's
/// own name (`CoplanarMerge`, the chord's corner as witness) until the merge learns
/// to thread mixed rings.
#[test]
fn census_corpus_cylinder_families_build_or_refuse_by_name() {
    fn tr(m: &mut Model, s: Handle<Solid>, v: [f64; 3]) -> Handle<Solid> {
        let s = transform(
            m,
            s,
            &nacre_scalar::Isometry::translation(v.map(|x| Rat::try_from_f64(x).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    }
    fn plate(m: &mut Model, hi: [f64; 3]) -> Handle<Solid> {
        m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(hi))
    }
    fn cyl(m: &mut Model, base: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
        m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            r,
            h,
        )
    }
    /// A `20 × 20 × 10` cell at `x0`, optionally pocketed, then bored (`ct2`).
    fn ct2cell(m: &mut Model, x0: f64, pocket: bool) -> Handle<Solid> {
        let plate = m.add_cuboid(
            Point3::from_array([x0, 0.0, 0.0]),
            Point3::from_array([x0 + 20.0, 20.0, 10.0]),
        );
        m.rebuild_adjacency();
        let mut out = plate;
        if pocket {
            let p = m.add_cuboid(
                Point3::from_array([x0 + 2.0, 2.0, 4.0]),
                Point3::from_array([x0 + 8.0, 8.0, 10.0]),
            );
            m.rebuild_adjacency();
            out = boolean(m, BoolKind::Cut, out, p).expect("the pocket cuts")[0];
            m.rebuild_adjacency();
        }
        let bore = cyl(m, [x0 + 14.0, 14.0, -1.0], 2.0, 12.0);
        m.rebuild_adjacency();
        out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    }
    type Fx = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);
    /// A family: its name, its two operands, and the outcome per kind (fuse, cut, common).
    type Family = (&'static str, Fx, [Result<usize, RejectReason>; 3]);
    let ok = |n: usize| Ok::<usize, RejectReason>(n);
    let fams: Vec<Family> = vec![
        (
            "rul offmid",
            |m| {
                let a = plate(m, [40.0, 40.0, 20.0]);
                let b = cyl(m, [40.0, 10.0, -10.0], 5.0, 50.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "rul flush",
            |m| {
                let a = plate(m, [40.0, 40.0, 20.0]);
                let b = cyl(m, [40.0, 20.0, 0.0], 5.0, 20.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "wal corner-lo",
            |m| {
                let a = plate(m, [4.0, 4.0, 2.0]);
                let b = cyl(m, [0.0, 0.0, -1.0], 0.5, 4.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "cap sunk",
            |m| {
                let a = plate(m, [4.0, 4.0, 2.0]);
                let b = cyl(m, [2.0, 2.0, 1.0], 0.5, 1.0);
                m.rebuild_adjacency();
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "trc boss",
            |m| {
                let a = plate(m, [40.0, 40.0, 10.0]);
                let b = cyl(m, [10.0, 20.0, 5.0], 3.0, 8.0);
                m.rebuild_adjacency();
                let b = tr(m, b, [15.3, 0.1, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(1)],
        ),
        (
            "trc onaxis",
            |m| {
                let p = plate(m, [40.0, 40.0, 10.0]);
                let bore = cyl(m, [18.0, 20.0, -5.0], 2.1, 30.0);
                m.rebuild_adjacency();
                let holed = boolean(m, BoolKind::Cut, p, bore).expect("the bore cuts")[0];
                m.rebuild_adjacency();
                let twin = cyl(m, [7.3, 20.0, -5.0], 2.1, 30.0);
                m.rebuild_adjacency();
                let twin = tr(m, twin, [10.7, 0.0, 0.0]);
                (holed, twin)
            },
            [Err(RejectReason::CylinderPairContact); 3],
        ),
        (
            "ct2 pocketed",
            |m| {
                let a = ct2cell(m, 0.0, true);
                let b = ct2cell(m, 20.0, true);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "ct2 chain",
            |m| {
                let a = ct2cell(m, 0.0, false);
                let b = ct2cell(m, 20.0, false);
                let ab = boolean(m, BoolKind::Fuse, a, b).expect("the first pair fuses")[0];
                m.rebuild_adjacency();
                let c = ct2cell(m, 40.0, false);
                (ab, c)
            },
            [ok(1), ok(1), ok(0)],
        ),
    ];
    for (name, fx, want) in &fams {
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
        {
            let mut m = Model::new();
            let (a, b) = fx(&mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Ok(n)) => {
                    assert_eq!(out.len(), *n, "{name} {kind:?}: solids");
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                }
                (Err(BoolError::Rejected { reason, .. }), Err(r)) => {
                    assert_eq!(&reason, r, "{name} {kind:?}: the refusal's name");
                }
                (got, _) => panic!("{name} {kind:?}: {got:?}, wanted {want:?}"),
            }
        }
    }
}

/// The `xy` families of the census corpus (rows of bored cells fused in two generations, and the
/// turned negative control), in the lib suite for the same reason as
/// [`census_corpus_cylinder_families_build_or_refuse_by_name`].
#[test]
fn census_corpus_xy_generations_build_or_refuse_by_name() {
    use nacre_scalar::{Angle, Isometry, Rotation};
    fn tr(m: &mut Model, s: Handle<Solid>, v: [f64; 3]) -> Handle<Solid> {
        let s = transform(
            m,
            s,
            &Isometry::translation(v.map(|x| Rat::try_from_f64(x).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    }
    fn cell(m: &mut Model) -> Handle<Solid> {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([20.0, 20.0, 10.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([6.3, 6.3, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.1,
            12.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    }
    fn row(m: &mut Model, n: usize, step: [f64; 3]) -> Handle<Solid> {
        let mut acc = cell(m);
        for i in 1..n {
            let c = cell(m);
            let c = tr(m, c, step.map(|v| v * i as f64));
            acc = boolean(m, BoolKind::Fuse, acc, c).expect("the row fuses")[0];
            m.rebuild_adjacency();
        }
        acc
    }
    type Fx = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);
    /// A family: its name, its two operands, and the outcome per kind (fuse, cut, common).
    type Family = (&'static str, Fx, [Result<usize, RejectReason>; 3]);
    let ok = |n: usize| Ok::<usize, RejectReason>(n);
    let fams: Vec<Family> = vec![
        (
            "xy grid",
            |m| {
                let a = row(m, 2, [20.0, 0.0, 0.0]);
                let b = row(m, 2, [20.0, 0.0, 0.0]);
                let b = tr(m, b, [0.0, 20.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy yx",
            |m| {
                let a = row(m, 2, [0.0, 20.0, 0.0]);
                let b = row(m, 2, [0.0, 20.0, 0.0]);
                let b = tr(m, b, [20.0, 0.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy three",
            |m| {
                let a = row(m, 3, [20.0, 0.0, 0.0]);
                let b = row(m, 3, [20.0, 0.0, 0.0]);
                let b = tr(m, b, [0.0, 20.0, 0.0]);
                (a, b)
            },
            [ok(1), ok(1), ok(0)],
        ),
        (
            "xy turned",
            |m| {
                let a = row(m, 2, [20.0, 0.0, 0.0]);
                let b = row(m, 2, [20.0, 0.0, 0.0]);
                let b = transform(
                    m,
                    b,
                    &Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        pivot: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                    }),
                )
                .unwrap();
                m.rebuild_adjacency();
                (a, b)
            },
            [Err(RejectReason::CylinderGateUndecided); 3],
        ),
    ];
    for (name, fx, want) in &fams {
        for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
        {
            let mut m = Model::new();
            let (a, b) = fx(&mut m);
            match (boolean(&mut m, kind, a, b), want) {
                (Ok(out), Ok(n)) => {
                    assert_eq!(out.len(), *n, "{name} {kind:?}: solids");
                    m.rebuild_adjacency();
                    let issues = nacre_validate::validate(&m);
                    assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                }
                (Err(BoolError::Rejected { reason, .. }), Err(r)) => {
                    assert_eq!(&reason, r, "{name} {kind:?}: the refusal's name");
                }
                (got, _) => panic!("{name} {kind:?}: {got:?}, wanted {want:?}"),
            }
        }
    }
}

/// ★★★★★ **A tool that reaches the recovered band — the lock the far cube cannot be.** The
/// re-operation census cuts each result with a cube far outside it, so a wrong band would pass
/// it unseen. Here the wall +x boss's result (a band whose notch met the seam and was spliced
/// into the outer walk; read back as a hole) is cut by a cuboid whose walls clear the
/// boss's strip and whose ⊥ caps meet the lateral. Volumes derived, not copied: the plate is
/// 4·4·2 = 32; the boss `r = 0.5` on the wall `x = 4` adds its outer half over its height 4
/// (`½·π/4·4 = π/2`) and its inner half outside the plate, below and above (`½·π/4·2 = π/4`) —
/// the inner half inside the plate is plate already — so the fuse is `32 + 3π/4`.
///
/// * A cuboid `[3, 6] × [−1, 5] × [2.5, 4]`: its lower cap crosses the band above the plate,
///   where the circle is whole (the band's rims and holes are read from its cycles, and the
///   notch below carves nothing here), and its upper cap is beyond the band. It removes the boss
///   over `z ∈ [2.5, 3]`, both halves: `π/8`. (A lower cap **on** the plate's top, `z = 2`,
///   refuses `NoClearRay` today — a coplanar-contact containment question, not this road's.)
/// * A cuboid `[3, 6] × [−1, 5] × [0.5, 1.5]`: its caps run **through** the notch, where the
///   band's circle is only the outer half (the Crossing arm) — and through the plate's wall face
///   `x = 4`, whose ring carries the boss's **rulings**; the scan names each crossing there as
///   a pierce node. It removes the plate's `x ∈ [3, 4]` slab and the boss's outer half over the
///   slab's height: `4 + π/8`.
#[test]
fn a_spliced_band_is_cut_across_its_notch() {
    let pi = std::f64::consts::PI;
    let build = || {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
        m.rebuild_adjacency();
        let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!(
            (v0 - (32.0 + 3.0 * pi / 4.0)).abs() < 1e-9,
            "the first op's volume: {v0}"
        );
        (m, out[0], v0)
    };
    // (a) caps across the band above the plate and beyond it.
    {
        let (mut m, r0, v0) = build();
        let tool = m.add_cuboid(
            Point3::from_array([3.0, -1.0, 2.5]),
            Point3::from_array([6.0, 5.0, 4.0]),
        );
        m.rebuild_adjacency();
        let r = boolean(&mut m, BoolKind::Cut, r0, tool)
            .unwrap_or_else(|e| panic!("the band is cut above the plate: {e:?}"));
        assert_eq!(r.len(), 1, "one solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
        assert!(
            (v - (v0 - pi / 8.0)).abs() < 1e-9,
            "{v} vs {}",
            v0 - pi / 8.0
        );
    }
    // (b) caps through the notch: the plate's wall face carries the boss's rulings, and the
    // caps' classes cross it there — the scan names each crossing as a pierce node. The
    // tool takes the plate's `x ∈ [3, 4]` slab (4) and the boss's outer half over the slab's
    // height (π/8); its wall `x = 3` clears the boss (1 > r). Beside the volume, the probe: every
    // crossing named in this binary lies on its class, its face plane and the cylinder, on the
    // ruling side the ring carried — the f64 road to the sign `crossing_on_ruling` chose by.
    {
        let (mut m, r0, v0) = build();
        let tool = m.add_cuboid(
            Point3::from_array([3.0, -1.0, 0.5]),
            Point3::from_array([6.0, 5.0, 1.5]),
        );
        m.rebuild_adjacency();
        let r = boolean(&mut m, BoolKind::Cut, r0, tool)
            .unwrap_or_else(|e| panic!("the notch is cut across its rulings: {e:?}"));
        assert_eq!(r.len(), 1, "one solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
        assert!(
            (v - (v0 - (4.0 + pi / 8.0))).abs() < 1e-9,
            "{v} vs {}",
            v0 - (4.0 + pi / 8.0)
        );
        let hits = crate::arrangement::crossing_probe::HITS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        // Two cap classes × the wall face's two halves, one ruling each.
        assert!(hits.len() >= 4, "the probe saw {} crossings", hits.len());
        for h in &hits {
            assert!(
                h.off.iter().all(|&d| d <= 1e-9),
                "a named crossing sits off its surfaces: {h:?}"
            );
            if !h.arc {
                assert_eq!(
                    h.side_f64, h.side,
                    "a named crossing lies on the other ruling: {h:?}"
                );
            }
        }
    }
}
