//! Cylinders in the arrangement: extent and order on a class, hole arcs, class-space names, disks.

use super::*;

/// **Two cylinder operations in a row.**
///
/// "A plate takes a second cylinder" has two halves, and this is the one with the first cylinder on
/// the plate's face, kept as its own guard.
///
/// A bore or a boss standing on the plate's face leaves the operand's rings plane-named (a hole is
/// a full circle, the one curved loop the tracer's road speaks), so the next cylinder is business
/// as usual. A boss standing on a **wall** does not;
/// [`a_chained_cylinder_bounded_by_the_first_builds`] holds that half.
///
/// The volumes are derived, not copied: the plate is 96, a bore through it removes `π/4·2`, a boss
/// on top adds `π/4·1`, and a wall boss adds `π/4·4` with half its section buried (`−½·π/4·2`).
#[test]
fn chained_cylinder_operations_on_a_face_build() {
    let quarter = std::f64::consts::PI * 0.25;
    let bore = ([2.0, 2.0, -1.0], 4.0, BoolKind::Cut);
    let top_boss = ([2.0, 2.0, 2.0], 1.0, BoolKind::Fuse);
    for (name, first, second, volume) in [
        (
            "bore then bore",
            bore,
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            96.0 - 2.0 * quarter * 2.0,
        ),
        (
            "bore then top boss",
            bore,
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 - quarter * 2.0 + quarter,
        ),
        (
            "top boss then top boss",
            top_boss,
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 + 2.0 * quarter,
        ),
        (
            "top boss then wall boss",
            top_boss,
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 + quarter + (quarter * 4.0 - 0.5 * quarter * 2.0),
        ),
        (
            "bore then wall boss",
            bore,
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 - quarter * 2.0 + (quarter * 4.0 - 0.5 * quarter * 2.0),
        ),
    ] {
        let (mut m, r) = chained(first, second);
        let out = r.unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(out.len(), 1, "{name}: one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{name}: validate");
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        assert!((v - volume).abs() < 1e-9, "{name}: {v} vs {volume}");
    }
    // A wall boss as a *middle* operation builds; its volumes live in
    // `a_chained_cylinder_bounded_by_the_first_builds`, beside the mechanism that opens it.
}

/// **The gate decides a boss-seated wall, and records the pair — and nothing else in the
/// crate would notice if it stopped.**
///
/// The clearance test's whole content is a `bool`, and on a chained operand its answer travels
/// exactly one place: a `false` sends the seated pair to the `d = 0` record-and-pass arm, which
/// puts it in `crossings`, which is what makes the tracer's ruling and chord arms fire at all. A
/// wrong `true` would refuse nothing, panic nowhere and change no output, so **the suite would stay
/// green with the pierce corner's answer thrown away entirely**. Measured, by throwing it away: 316
/// green. This is the lock that sees it.
///
/// **What it sees is that the corner is *readable*, not what it says** — recorded rather than
/// claimed. Forcing every pierce corner to one strip side leaves even this lock green, because the
/// verdict here is settled on the *other* separating axis: the wall's plate corners already sit
/// inside the boss's axial span, so `along` fails whatever the strip half answers. The corner's
/// side becomes decisive only for a face that clears along the axis and has to be judged across
/// it, and no fixture builds one yet. Removing the read itself is what turns this red
/// (`CurvedOperandBoundary`, measured).
#[test]
fn the_gate_records_a_wall_the_boss_is_seated_on() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let seam = Vector3::from_array([0.0, -1.0, 0.0]);
    let a = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 0.0, -1.0]),
        up,
        seam,
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([6.0, 2.0, -1.0]),
        up,
        seam,
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    // The gate *deciding* is half the claim — past the first corner the boss makes.
    let setup = crate::arrangement::plane_index_setup(&m, first, b).expect("the gate decides");
    // And recording is the other half. The boss sits **on** `y = 0`, so those faces are not
    // clear of it; that is what the record-and-pass arm is for, and a missing entry here means the
    // clearance test answered "clear" about a wall a cylinder is standing in.
    //
    // The pair is checked by **what it means**, not by its indices: the arm records exactly the
    // classes whose plane the cylinder's axis lies *on*, so that is the assertion — renumbering
    // the classes cannot make it pass or fail for the wrong reason.
    assert_eq!(setup.crossings.len(), 1, "one seated pair, not more");
    for &(c, ci) in &setup.crossings {
        let coeffs = setup.geom[c].world_rat().expect("a named wall class");
        assert_eq!(
            nacre_exact::point_plane_clearance_rat(
                &coeffs,
                &setup.cyls[ci].def.origin(),
                &nacre_exact::BigRat::zero()
            ),
            nacre_exact::Orient::Zero,
            "the recorded pair is a cylinder seated exactly on that class's plane"
        );
    }
}

/// **A holed lateral's ruling comes out in three pieces, and the middle one grazes.**
///
/// A wall through the boss's axis meets its lateral in two rulings, and where the boss is buried in
/// the plate the lateral does not *cross* that wall — it ends at it. So the mark is
/// `Transversal · Graze · Transversal` along the axis, not one full-height crossing.
///
/// **Held at the trace**: the operation builds, but a result says only that a ruling went
/// wrong, not which piece of it; the trace's own answer says that, the way `arrangement`'s
/// `sides == [-1, 1]` lock holds the hole-free case.
///
/// The claim is over every ruling **this fixture** carves on its boss, counted on the rows its
/// own operation recorded.
#[test]
fn a_holed_laterals_ruling_grazes_where_the_hole_is() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let seam = Vector3::from_array([0.0, -1.0, 0.0]);
    let a = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([2.0, 0.0, -1.0]),
        up,
        seam,
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([6.0, 2.0, -1.0]),
        up,
        seam,
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    crate::ledger::owned(|| boolean(&mut m, BoolKind::Cut, first, b).expect("the bore builds"));
    // This operation's own carvings — but it carves other classes of the same model too (a
    // panel's ruling is one graze, a chain's a transversal then a graze), so the boss is still
    // named: origin `(2, 0, −1)`, rulings spanning `t ∈ [0, 4]`.
    let carved: Vec<Vec<crate::arrangement::SegKind>> = crate::arrangement::ruling_probe::CARVED
        .mine()
        .iter()
        .filter(|c| c.origin == [2.0, 0.0, -1.0] && c.span == [0.0, 4.0])
        .map(|c| c.kinds.clone())
        .collect();
    assert!(
        !carved.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    for kinds in &carved {
        assert!(
            matches!(
                kinds[..],
                [
                    crate::arrangement::SegKind::Transversal { .. },
                    crate::arrangement::SegKind::Graze { .. },
                    crate::arrangement::SegKind::Transversal { .. }
                ]
            ),
            "a carved ruling came out as {kinds:?}"
        );
    }
    // Both rulings are carved, not just one: the hole has a vertical edge on each. A floor,
    // not a count: a `debug_assertions` build traces the arrangement twice (the reuse-off
    // reference).
    assert!(carved.len() >= 2, "only {} ruling carved", carved.len());
    // **And which side the face occupies, against the fixture's own geometry.** The buried
    // half of the boss is the `y > 0` one, so along the hole's vertical edges the lateral survives
    // at `y < 0` — the stored-normal side exactly when that normal points at `−y`. Flipping any
    // factor of the derivation turns this red; flipping `body_above` itself also fails the tests
    // of the booleans it feeds.
    let sides: Vec<(bool, f64)> = crate::arrangement::ruling_probe::GRAZE_SIDE
        .mine()
        .iter()
        // Only the **hole's** runs: this same operation puts a *panel* on the same wall and
        // span, and its face lies on the other side — measured (`body_above` false against the
        // hole's true), which is what a filter by origin alone let in.
        .filter(|g| g.inner && g.origin == [2.0, 0.0, -1.0] && g.span[1] - g.span[0] == 2.0)
        .map(|g| (g.body_above, g.ny))
        .collect();
    assert!(!sides.is_empty(), "the side probe never ran");
    for &(body_above, ny) in &sides {
        assert!(
            ny.abs() > 0.5,
            "the wall's stored normal is not axis-aligned: {ny}"
        );
        assert_eq!(
            body_above,
            ny < 0.0,
            "a hole's ruling graze states the wrong side of the wall: {body_above} against ny={ny}"
        );
    }
}

/// **Which of the two arcs is the hole — measured, not argued.**
///
/// A plane cuts a circle in two points, and the whole difficulty of a lateral face's hole is which
/// of the two arcs it is: the hole's two ruling edges lie on **one** plane here, so that plane's
/// sides cannot tell them apart. `cycle_on_class` answers from the ring's own winding (material is
/// on the left of the ring's travel), and gets there through three signs — the face's
/// `orient_sign`, the class's `frame_sign`, and whether the class's stored normal points along the
/// axis. Reversed, the answer is not approximate but **exactly opposite**: the graze would sit on
/// the band and the crossing inside the hole.
///
/// **Downstream notices only the symptom** — reversed, a chained operation is refused as
/// `LabelConflict` and emission's disk-side check fires, neither naming the arc — so the arc is
/// realized and judged here as well. The kernel reads no coordinate to choose it; this reads one
/// afterwards.
///
/// **Which factors this actually sees.** Turning each one off in turn: reversing either arm's
/// winding is **red**, dropping `plus_t_is_above` is **red**, and dropping the root restatement
/// (the ruling named for a different ⊥ plane) is **red**. Dropping `frame_sign` or the face's
/// `orient_sign` is **green** here — both are `+1` on a boss's holed lateral — and
/// [`a_bores_lateral_with_a_hole_takes_a_second_cylinder`] is the lock that turns them.
///
/// The claim is over this test's own arcs on the wall boss (origin `(2, 0, −1)`): the plate
/// stands on `y > 0` and the boss is centred on the wall `y = 0`, so the buried half — the hole —
/// is the `y > 0` one.
#[test]
fn a_holes_arcs_run_the_way_the_hole_lies() {
    // Two second operands, because they reach **different arms** of the walk. A cylinder rising
    // from `z = −1` puts its cap on `z = 0`, which is the hole's own low rim — an on-line *run*.
    // One rising from `z = 1` puts a cap **strictly inside** the hole's `z ∈ (0, 2)`, where the
    // ring crosses the class on its two ruling edges instead — the *crossing* arm, whose extra
    // step is restating the ruling's root for a different ⊥ plane.
    for base in [-1.0, 1.0] {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([12.0, 4.0, 2.0]),
        );
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let seam = Vector3::from_array([0.0, -1.0, 0.0]);
        let a = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([2.0, 0.0, -1.0]),
            up,
            seam,
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
        m.rebuild_adjacency();
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([6.0, 2.0, base]),
            up,
            seam,
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        crate::ledger::owned(|| {
            let _ = boolean(&mut m, BoolKind::Cut, first, b);
        });
    }
    // This test's own carvings, and among them only this boss (`(2, 0, −1)`) and only its
    // **hole** — the same models carve panels' and chains' arcs, which face any way, and the
    // sentence is about holes.
    let mids: Vec<[f64; 3]> = crate::arrangement::arc_probe::MIDS
        .mine()
        .iter()
        .filter(|m| m.inner && m.origin == [2.0, 0.0, -1.0])
        .map(|m| m.dir)
        .collect();
    assert!(
        !mids.is_empty(),
        "the probe never ran, so it measured nothing"
    );
    assert!(
        mids.iter().all(|d| d[1] > 0.0),
        "a hole's arc was stated running the wrong way round: midpoints {mids:?}"
    );
}

/// **A bore's lateral with a hole takes a second cylinder — the two signs the boss family
/// cannot turn.**
///
/// A tongue fused into a bored plate fills a patch of the bore's lateral, so that face comes back
/// with a hole whose material lies **outside** the cylinder (`orient_sign = −1`), where a boss's
/// lies inside. The second operation reads the holed lateral against its ⊥ classes in
/// `cycle_on_class`, through both arms: the tongue's own faces are the hole's rims (the run arm),
/// and a far cylinder's cap strictly inside the hole's height is crossed on its two rulings (the
/// crossing arm). The third row moves everything down by one so that cap lies on the seed plane
/// `z = 0` facing `+z`, which makes it a `Reversed` face — and the crossing class carries
/// `frame_sign = −1`, which the crossing arm reads to tell the travel's axis sense.
///
/// Measured by dropping each factor in `cycle_on_class`: without `orient_sign` every row is
/// refused (`LabelConflict`) and the hole's arcs come back on the `−y` side; without the class's
/// `frame_sign` the third row is refused (`CylinderFaceUndecided`), and it is the only row in the
/// suite, the ignored sweep and the census whose crossing arm meets a `frame_sign = −1` class.
///
/// The volumes are derived: the plate is 96, the bore removes `2π`, and the tongue adds the part
/// of the bore it fills — the strip `|x − 6| ≤ 0.6` of the half disk `y ≥ 2`, `0.6·0.8 + asin 0.6`,
/// one high (`0.6` is a 3-4-5 chord, so the rulings are rational). The second cylinder removes
/// `π/4` per unit of the plate it passes.
#[test]
fn a_bores_lateral_with_a_hole_takes_a_second_cylinder() {
    let q = std::f64::consts::PI * 0.25;
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let seam = Vector3::from_array([0.0, -1.0, 0.0]);
    let first_volume = 96.0 - 8.0 * q + 0.48 + 0.6f64.asin();
    // (name, the whole fixture's shift along z, the second cylinder's base, the plate it removes)
    for (name, dz, base, removed) in [
        ("the run arm", 0.0, -1.0, 2.0),
        ("the crossing arm", 0.0, 1.0, 1.0),
        ("the crossing arm on a seed plane", -1.0, -4.0, 1.0),
    ] {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 0.0, dz]),
            Point3::from_array([12.0, 4.0, 2.0 + dz]),
        );
        let origin = [6.0, 2.0, -1.0 + dz];
        let bore = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array(origin),
            up,
            seam,
            1.0,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        let bored = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore")[0];
        m.rebuild_adjacency();
        let tongue = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([5.4, 2.0, 0.5 + dz]),
            Point3::from_array([6.6, 3.5, 1.5 + dz]),
        );
        m.rebuild_adjacency();
        let first = boolean(&mut m, BoolKind::Fuse, bored, tongue)
            .unwrap_or_else(|e| panic!("{name}: the tongue: {e:?}"));
        assert_eq!(first.len(), 1, "{name}: the tongue makes one solid");
        m.rebuild_adjacency();
        let v = nacre_props::mass_props(&m, first[0]).expect("mass").volume;
        assert!(
            (v - first_volume).abs() < 1e-9,
            "{name}: {v} vs {first_volume}"
        );
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([2.0, 2.0, base]),
            up,
            seam,
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        // The premise of the third row: the cap that crosses the hole is `Reversed`, and on the
        // others it is not.
        let top = m
            .shell(m.solid(b).outer)
            .faces
            .iter()
            .map(|&f| m.face(f))
            .find(|f| match m.surface_cache(f.surface) {
                nacre_topo::SurfaceCache::Plane { plane: p, .. } => {
                    (p.normal() * f64::from(f.orientation.sign())).as_array() == [0.0, 0.0, 1.0]
                }
                _ => false,
            })
            .expect("the second cylinder's top cap");
        assert_eq!(
            top.orientation == nacre_topo::Orientation::Reversed,
            base + 4.0 == 0.0,
            "{name}: the cap on the seed plane, and only it, is Reversed"
        );
        let from = crate::arrangement::arc_probe::MIDS.len();
        let out = crate::ledger::owned(|| boolean(&mut m, BoolKind::Cut, first[0], b))
            .unwrap_or_else(|e| panic!("{name}: the second cylinder: {e:?}"));
        assert_eq!(out.len(), 1, "{name}: one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{name}: validate");
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        let want = first_volume - q * removed;
        assert!((v - want).abs() < 1e-9, "{name}: {v} vs {want}");
        // The tongue fills the bore on its `y > 2` side, so that is where the hole is. This row's
        // own arcs only, so an earlier row's cannot say the probe ran here.
        let mids: Vec<[f64; 3]> = crate::arrangement::arc_probe::MIDS
            .mine_since(from)
            .iter()
            .filter(|mid| mid.inner && mid.origin == origin)
            .map(|mid| mid.dir)
            .collect();
        assert!(!mids.is_empty(), "{name}: the probe never ran");
        assert!(
            mids.iter().all(|d| d[1] > 0.0),
            "{name}: a hole's arc was stated running the wrong way round: midpoints {mids:?}"
        );
    }
}

/// **A cycle that carves two extents on one class has both read.** The wall `[9,11]×[0,40]×[5,15]`
/// fused across the bored plate ([`walled_bore`]) opens two windows in the bore's lateral, one on
/// either side of the axis, from `z = 5` up to the plate's top, and the lateral's boundary around
/// them is one cycle. `cycle_on_class` reads it against each ⊥ class: on `z = 10` and on `z = 5` it
/// lays two arcs (the run arm), and a box standing on `z = 8` — far from the bore, there only for
/// its plane — crosses its four rulings, two stretches, one per window (the crossing arm).
///
/// Measured by dropping the last extent of a call that carved two or more: with only the arcs'
/// or only the stretches' dropped, every operation here is refused (`LabelConflict`). The
/// `nth ≥ 1` rows are the fixture's premise rather than a claim about the kernel — a change that
/// rightly read this cycle in one extent per class would turn them red without anything being
/// wrong, and then this case needs another fixture.
#[test]
fn a_lateral_cycle_with_two_extents_on_a_class_carves_both() {
    let walled = walled_bore_volume();
    for (kind, want) in [
        (BoolKind::Fuse, walled + 200.0),
        (BoolKind::Cut, walled - 50.0),
        (BoolKind::Common, 50.0),
    ] {
        let mut m = Model::new();
        let a = walled_bore(&mut m);
        let b = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([30.0, 30.0, 8.0]),
            Point3::from_array([35.0, 35.0, 18.0]),
        );
        m.rebuild_adjacency();
        let from = crate::arrangement::arc_probe::MIDS.len();
        let out = crate::ledger::owned(|| boolean(&mut m, kind, a, b))
            .unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        assert_eq!(out.len(), 1, "{kind:?}: one solid");
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{kind:?}: validate"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        assert!((v - want).abs() <= 1e-9 * want, "{kind:?}: {v} vs {want}");
        let second: Vec<bool> = crate::arrangement::arc_probe::MIDS
            .mine_since(from)
            .iter()
            .filter(|mid| mid.inner && mid.origin == [10.0, 10.0, -1.0] && mid.nth >= 1)
            .map(|mid| mid.grazes)
            .collect();
        for (arm, grazes) in [("an arc", true), ("a crossed stretch", false)] {
            assert!(
                second.contains(&grazes),
                "{kind:?}: no cycle carved {arm} as its second extent on a class — the fixture no \
                 longer reaches the case it locks"
            );
        }
    }
}

/// **A boolean's own result takes a second cylinder — the whole road, end to end.**
///
/// A boss standing on a **wall** buries half the boss's lateral in the plate, so that face comes
/// back as a band with a **hole**, and every layer downstream has to read that band right.
///
/// **A label does not answer two questions.** A label says where *material* is, not whether the
/// lateral face is **there**. On the first fusion the two sectors of the holed rim carry different
/// labels (the buried one is inside the *other* solid), so read as existence they would happen to
/// be right; on a second operation the plate is already own material and the two sectors' labels
/// are **literally identical** — read as existence, both would survive, filling the hole in and
/// leaving its four boundary edges claimed once each. The trace says which sector is a face, per
/// arc, by angular extent: `Graze` where the face ends, `Transversal` where it runs through.
/// `cyl_chart::read_cell`'s `face_spans` asks it.
///
/// Measured across these four fixtures: **four** sectors dropped for existence, exactly one per
/// operation — the buried half — and the volumes below are what says that count is right.
///
/// The volumes are derived, not copied: the plate is 96, the wall boss adds `π/4·4` with half its
/// section buried (`−½·π/4·2`), a through bore removes `π/4·2`, a boss on top adds `π/4·1`, and a
/// bore *on the wall* removes only the half-section the plate holds (`−½·π/4·2`).
#[test]
fn a_chained_cylinder_bounded_by_the_first_builds() {
    let q = std::f64::consts::PI * 0.25;
    let wall_boss = ([2.0, 0.0, -1.0], 4.0, BoolKind::Fuse);
    for (name, second, volume) in [
        (
            "bore",
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            96.0 + 3.0 * q - 2.0 * q,
        ),
        (
            "boss on top",
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            96.0 + 3.0 * q + q,
        ),
        (
            "another wall boss",
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            96.0 + 3.0 * q + 3.0 * q,
        ),
        (
            "a bore on the same wall",
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Cut),
            96.0 + 3.0 * q - q,
        ),
    ] {
        let (mut m, r) = chained(wall_boss, second);
        let out = r.unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(out.len(), 1, "{name}: one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "{name}: validate");
        let v = nacre_props::mass_props(&m, out[0]).expect("mass").volume;
        assert!((v - volume).abs() < 1e-9, "{name}: {v} vs {volume}");
    }
    // **And that the existence gate is what did it.** The chart's census holds it where the
    // cells are read:
    // `exist_marks_false` (a sector dropped because the face is not there — the four above, held
    // as growth ≥ 4 in `cyl_chart::tests::a_present_cell_has_a_horizontal_answer_and_the_ledger_is_total`)
    // and the census's one-mark assertions (a cut end read carries at most one lateral mark of
    // its own solid, and none only where the cell is absent — the doc that says why the `Seated`
    // skip is load-bearing lives on `face_spans`).
}

/// **A boolean's result carries names the tracer's plane-only road cannot read — and it
/// says so instead of aborting.**
///
/// A boss standing on a wall leaves the plate's caps bitten by an **arc** and the wall split in two
/// by the boss's **rulings**, so those faces' rings run along a cylinder. [`combinatorics`]'s ring
/// naming asks two plane-only questions of every edge — the carried wall class and the vertex's
/// three-plane name — and both go through `ClassIx::plane`, which **panics** on a cylinder. That
/// accessor is right to: for its forty-odd other callers a cylinder there is an upstream filter
/// bug. This is the caller whose input is an *operand*, so this is where the filter belongs.
///
/// **Without the filter this very call panics** ("a plane-only path got cylinder class 0",
/// measured). The population gate refuses such an operand long before the tracer sees it, so
/// nothing in production reaches this — which is exactly why the lock calls the function
/// **directly**, the same way the ring-naming goldens next door do.
///
/// The untouched walls are the negative control: they are still named, so the filter is a filter
/// and not a blanket refusal.
#[test]
fn an_operand_bounded_by_a_cylinder_is_named_in_class_space() {
    // **The corner boss is here because the wall boss is symmetric.** Its axis lies *on* the
    // wall, so the two rulings sit either side of that plane and a root read the wrong way round
    // lands on a mirror-image point — a lock that only ever saw this fixture could pass with the
    // restatement inverted. The corner's two walls break that symmetry.
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for at in [[2.0, 0.0, -1.0], [12.0, 0.0, -1.0]] {
        let (u, s) = named_in_class_space(at);
        ups.extend(u);
        sides.extend(s);
    }
    // **The direction assertions inside are only worth their ink if the fixtures move them** —
    // a boss whose rulings all climbed, or all sat one side, would let a constant pass for a
    // derivation. It takes **both** fixtures to move both bits, and that is the geometry rather
    // than a gap: the wall boss is symmetric about its wall, so its two rulings straddle that
    // plane and `side` takes both values; the corner boss keeps **one** ruling per wall (the
    // other is buried), and those two are measured against *different* planes, so nothing says
    // they should oppose. Measured — the corner run alone gives `[1, 1]`.
    assert!(
        ups.contains(&true) && ups.contains(&false),
        "the fixtures exercise both senses of `up`: {ups:?}"
    );
    assert!(
        sides.contains(&1) && sides.contains(&-1),
        "the fixtures exercise both sides of the axis plane: {sides:?}"
    );
}

fn named_in_class_space(at: [f64; 3]) -> (Vec<bool>, Vec<i8>) {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array(at),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
    m.rebuild_adjacency();
    let r = out[0];

    let faces_tab = collect_planes(&m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).expect("edge incidence");
    // The cylinder table the tracer would hold, built the way `cylinder_gate` builds it — the gate
    // itself cannot be asked here, because refusing this very operand is its job.
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = m.world_cylinder_def(surf).expect("a world cylinder");
            let Some(cache) = m.surface_cache(surf).cylinder() else {
                unreachable!("a cylinder class names a cylinder")
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                realized: *cache,
                owner: crate::planes::SolidSide::A,
                shared: Vec::new(),
                silent: Vec::new(),
                boxes: Vec::new(),
            }
        })
        .collect();

    let jd = crate::planes::test_judge(&planes);
    let (mut curved_walls, mut pierce_names) = (0usize, 0usize);
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        if matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            continue; // a lateral's loops are read below, every one of them
        }
        let ring = combinatorics::face_vertex_triples(&m, fh, fp, &inc, &jd, &plane_ix, &cyls)
            .unwrap_or_else(|e| panic!("face {:?}: {e:?}", fh.index()));
        let Some(nr) = ring.poly() else { continue };
        let hes = &m.face(fh).outer.half_edges;
        let n = hes.len();
        assert_eq!(
            nr.triples.len(),
            n,
            "face {:?}: one name per corner",
            fh.index()
        );
        let ends = |k: usize| m.edge(hes[k].edge).vertices;
        // The corner the walk stands on when it starts edge `k` — the vertex edge `k-1` and edge
        // `k` share. Read off the ring rather than off either edge's stored pair, because a
        // ruling's pair carries no order (`edge_for` keys it unordered) and this is the walk.
        let corner = |k: usize| {
            let (a, b) = (ends((k + n - 1) % n), ends(k));
            *a.iter()
                .find(|x| b.contains(x))
                .expect("consecutive edges share a corner")
        };
        for (i, &node) in nr.triples.iter().enumerate() {
            // The carrier and the corner are independent facts, and both are asserted: a curved
            // wall must be an arc or a ruling of the cylinder its far face is on, never a plane.
            let straight = matches!(m.edge_curve(hes[i].edge), nacre_geom::Curve::Line(_));
            match nr.walls[i] {
                crate::combinatorics::Wall::Plane(_) => {}
                crate::combinatorics::Wall::Arc { plane, .. } => {
                    assert!(
                        !straight,
                        "face {:?} edge {i}: an arc carrier on a straight curve",
                        fh.index()
                    );
                    assert_eq!(
                        plane,
                        plane_ix[fp].plane(),
                        "face {:?} edge {i}: a cap's arc rides the cap's own plane",
                        fh.index()
                    );
                    curved_walls += 1;
                }
                // **The direction bits, against the walk and against the geometry.** The
                // carrier alone says *which surface*; `up` and `side` say *which way* and *which
                // of the two rulings*, and a consumer that reads them wrong builds a face wound
                // backwards or seated on the far side of the cylinder. Both are measured off
                // realized coordinates on purpose — the code derives them without any (`side`
                // exactly through `quad::plane_side`, `up` from the cutting planes' axial
                // parameters), so the coordinate is a genuinely second road to the same bit.
                crate::combinatorics::Wall::Ruling {
                    cyl: k,
                    side,
                    up,
                    plane,
                } => {
                    assert!(
                        straight,
                        "face {:?} edge {i}: a ruling carrier on a curved curve",
                        fh.index()
                    );
                    let crate::planes::ClassIx::Plane(near) = plane_ix[fp] else {
                        unreachable!("a lateral face was skipped above")
                    };
                    assert_eq!(
                        plane,
                        near,
                        "face {:?} edge {i}: a plane face's ruling is measured against the face",
                        fh.index()
                    );
                    let axis = cyls[k].realized.axis();
                    let axial = |p: Point3| (p - axis.origin()).dot(axis.direction());
                    assert_eq!(
                        up,
                        axial(m.vertex_point(corner((i + 1) % n)))
                            > axial(m.vertex_point(corner(i))),
                        "face {:?} edge {i}: `up` disagrees with the walk",
                        fh.index()
                    );
                    // `side` is the sign of `(x − o) · (m × n̂)`. **`n̂` has to be the class's
                    // `world_rat` normal, and no other spelling of the class's plane will do** —
                    // measured: `plane.normal()` here would come out opposed (`[0,-1,0]` against
                    // `[0,1,0]`, and `frame_sign` is `+1`, so that is not the reconciliation
                    // either). A class has no outward normal to agree on: `world_rat` is the
                    // plane's *name*, whose sign is its own, while `plane` is a stored surface's.
                    // `side` is a **label** telling the two rulings apart, so it is well defined
                    // exactly as long as every road spells `n̂` the one way. Sharing that input
                    // leaves the two roads independent where it counts — the code decides in exact
                    // `quad::plane_side` on the pierce meet, this in `f64` on the realized vertex.
                    let wr = combinatorics::class_coeffs_rat(&jd, near).expect("a named class");
                    let n_hat =
                        Vector3::from_array([wr[0].to_f64(), wr[1].to_f64(), wr[2].to_f64()]);
                    let cross = axis.direction().cross(n_hat);
                    let s = (m.vertex_point(corner(i)) - axis.origin()).dot(cross);
                    assert_eq!(
                        side,
                        if s > 0.0 { 1 } else { -1 },
                        "face {:?} edge {i}: `side` disagrees with the geometry ({s})",
                        fh.index()
                    );
                    ups.push(up);
                    sides.push(side);
                    curved_walls += 1;
                }
            }
            let Some((_, cyl, _)) = combinatorics::pierce_name(node) else {
                continue;
            };
            pierce_names += 1;
            // **Two roads, one point.** The name is this boolean's classes; the coordinate is
            // what the *previous* boolean realized from *its* classes. Realizing the one and
            // measuring it against the other is what says the restatement — the pair's order and
            // each normal's sign — came out right: get either wrong and the name designates the
            // **other root**, a visibly different point on the far ruling.
            let got = combinatorics::pierce_point(&jd, cyl, &cyls[cyl].def, node)
                .expect("the name realizes");
            let want = m.vertex_point(corner(i)).as_array();
            let d: f64 = (0..3)
                .map(|k| (got[k] - want[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(
                d < 1e-9,
                "face {:?} corner {i}: the name realizes at {got:?}, the vertex is at {want:?}",
                fh.index()
            );
        }
    }
    // **A lateral's own curved walls, against the geometry** — every loop, for a lateral's
    // rulings sit in whichever loop meets the wall (here the hole the plate leaves). Nothing on the
    // operand side reads a lateral ruling's `up`, so this is the one lock on it; `side` is checked
    // against the plane the wall carries, and both ends of an arc or a ruling must be corners whose
    // names name that plane — exact, and what a wrong `plane` breaks however its sign falls.
    let mut lateral_rulings = 0usize;
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        let crate::planes::ClassIx::Cyl(k) = plane_ix[fp] else {
            continue;
        };
        let rings = combinatorics::lateral_cycles(&m, fh, fp, &inc, &jd, &plane_ix, &cyls)
            .unwrap_or_else(|e| panic!("lateral {:?}: {e:?}", fh.index()));
        let face = m.face(fh);
        // One ring per loop, or the `zip` below drops the loops past the shorter side unread.
        assert_eq!(
            rings.len(),
            1 + face.inner.len(),
            "lateral {:?}: a ring per loop",
            fh.index()
        );
        for (lp, ring) in std::iter::once(&face.outer).chain(&face.inner).zip(&rings) {
            let Some(nr) = ring.poly() else { continue };
            let hes = &lp.half_edges;
            let n = hes.len();
            let ends = |j: usize| m.edge(hes[j].edge).vertices;
            let corner = |j: usize| {
                let (a, b) = (ends((j + n - 1) % n), ends(j));
                *a.iter()
                    .find(|x| b.contains(x))
                    .expect("consecutive edges share a corner")
            };
            let names = |j: usize, c: usize| {
                combinatorics::pierce_name(nr.triples[j]).is_some_and(|(pl, _, _)| pl.contains(&c))
            };
            for (i, he) in hes.iter().enumerate() {
                let straight = matches!(m.edge_curve(he.edge), nacre_geom::Curve::Line(_));
                let (cyl, plane) = match nr.walls[i] {
                    crate::combinatorics::Wall::Plane(_) => {
                        panic!(
                            "lateral {:?} edge {i}: a lateral's ring has no plane wall",
                            fh.index()
                        )
                    }
                    crate::combinatorics::Wall::Arc { cyl, plane, .. } => {
                        assert!(
                            !straight,
                            "lateral {:?} edge {i}: an arc on a line",
                            fh.index()
                        );
                        (cyl, plane)
                    }
                    crate::combinatorics::Wall::Ruling {
                        cyl,
                        side,
                        up,
                        plane,
                    } => {
                        assert!(
                            straight,
                            "lateral {:?} edge {i}: a ruling on a curve",
                            fh.index()
                        );
                        let axis = cyls[cyl].realized.axis();
                        let axial = |p: Point3| (p - axis.origin()).dot(axis.direction());
                        assert_eq!(
                            up,
                            axial(m.vertex_point(corner((i + 1) % n)))
                                > axial(m.vertex_point(corner(i))),
                            "lateral {:?} edge {i}: `up` disagrees with the walk",
                            fh.index()
                        );
                        let wr =
                            combinatorics::class_coeffs_rat(&jd, plane).expect("a named class");
                        let n_hat =
                            Vector3::from_array([wr[0].to_f64(), wr[1].to_f64(), wr[2].to_f64()]);
                        let s = (m.vertex_point(corner(i)) - axis.origin())
                            .dot(axis.direction().cross(n_hat));
                        assert_eq!(
                            side,
                            if s > 0.0 { 1 } else { -1 },
                            "lateral {:?} edge {i}: `side` disagrees with the geometry ({s})",
                            fh.index()
                        );
                        ups.push(up);
                        sides.push(side);
                        lateral_rulings += 1;
                        (cyl, plane)
                    }
                };
                assert_eq!(
                    cyl,
                    k,
                    "lateral {:?} edge {i}: its own cylinder",
                    fh.index()
                );
                assert!(
                    names(i, plane) && names((i + 1) % n, plane),
                    "lateral {:?} edge {i}: both ends lie on the carried plane {plane}",
                    fh.index()
                );
            }
        }
    }
    // The plate's wall crosses the boss in two rulings, both in the lateral's one cut loop (the
    // other two are its whole rims).
    assert_eq!(
        lateral_rulings, 2,
        "the wall's two rulings on the boss's lateral"
    );
    // Four faces run along the boss — the two plate caps it bit an arc out of, and the two halves
    // its rulings split the wall into — and each contributes two curved edges and two pierce
    // corners. The boss's own caps are full circles (the one curved loop this road already spoke)
    // and the three untouched walls are plain.
    assert_eq!(curved_walls, 4, "edges riding the boss");
    assert_eq!(pierce_names, 8, "corners named as pierce points");
    (ups, sides)
}

/// **A cylinder-pinned end is ordered, not refused — and the sense is the geometry's.**
///
/// `edge_dir` is the one place a direction is made. A pierce point has no third plane, and the
/// integer predicate wants one, so its cylinder arm orders through the `a + b√c` tower; the 48
/// pins below fix what the answer must be (an arm that refuses by `RingNaming` refuses every one
/// of them, measured).
///
/// **The oracle is a second road, and only the half that matters is second.** The direction
/// `n_p × n_q` is read from the same raw coefficients the code reads — deliberately, the way the
/// ruling lock next door shares `world_rat`: a *label*'s reference frame has to be one spelling or
/// the two roads are not comparing the same thing. What is independent is the part under test —
/// the **order of the two points on that line**: the code decides it exactly (a rational meet
/// against a pierce root, or two roots against each other), the oracle realizes both points in
/// `f64` and subtracts. Measured `|t|` from 1.5 to 456, so nothing here is decided in the noise.
///
/// **The population is one-sided and that is a fact, not a blind instrument.** All 48 order
/// `-1`: a face's outer ring travels counter-clockwise about that face's own outward normal, which
/// is the direction `order_along` sorts by, so agreement is structural — two boss positions, one
/// of them with its axis reversed, and three `Cut` notches do not move it. What moves is **swapping the two ends**,
/// and the oracle swaps with it, so both signs are measured against something that could disagree.
#[test]
fn a_cylinder_pinned_end_orders_through_the_tower() {
    let mut pins = 0usize;
    for (at, dir, kind) in [
        ([2.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Fuse),
        ([12.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Fuse),
        ([2.0, 0.0, 3.0], [0.0, 0.0, -1.0], BoolKind::Fuse),
        ([2.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
        ([12.0, 0.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
        ([6.0, 4.0, -1.0], [0.0, 0.0, 1.0], BoolKind::Cut),
    ] {
        pins += pinned_ends_ordered(at, dir, kind);
    }
    // The count is the lock on the *population*: let the fixtures stop producing pierce-pinned
    // ends and every assertion below would pass vacuously. **48** — eight from each of the six
    // fixtures (measured), so any one of them going quiet moves the sum.
    assert_eq!(pins, 48, "cylinder-pinned ring ends exercised");
}

/// **A circle can be an interior boundary, and then the merge erases it.**
///
/// A tool that only *touches* removes nothing, and the planar engine says so plainly: a plate cut
/// by a box resting on its top face comes back as the plate — same volume, **six faces**, no
/// imprint left behind (measured, both for face contact and for a box straddling an edge). That
/// measurement is the oracle here; there is no sentence in the design documents that decides it.
///
/// The cylinder twin, locked against the top plane coming back as **two** faces: the plate's top
/// with a circular hole, and the disk filling it. A disk's boundary is a `Bound::Circle` with no
/// nodes at all, so a merge that reads only its edge table cannot see it; the corners left on that
/// circle then name the tool's cylinder after every face of it is gone, and the assembly refuses
/// the whole boolean (`VertexNamesAbsentSurface`).
///
/// What connects them is that the other face holds **the same cylinder class** as a hole — and one
/// plane class with one cylinder class names one circle, so that coincidence *is* adjacency.
#[test]
fn a_disk_merges_into_the_face_it_lies_in() {
    let plate_and_boss = |m: &mut Model, base: [f64; 3], h: f64| {
        let a = crate::fixtures::cuboid(
            m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            h,
        )
        .solid;
        m.rebuild_adjacency();
        (a, b)
    };
    let faces_of = |m: &Model, s: Handle<Solid>| -> usize {
        let sol = m.solid(s).clone();
        std::iter::once(&sol.outer)
            .chain(sol.cavities.iter())
            .map(|&sh| m.shell(sh).faces.len())
            .sum()
    };
    let run = |base: [f64; 3], h: f64, kind: BoolKind| -> (f64, usize) {
        let mut m = Model::new();
        let (a, b) = plate_and_boss(&mut m, base, h);
        let out = boolean(&mut m, kind, a, b).expect("the boolean builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // The merged face has to survive tessellation too: `validate` reads the topology store,
        // and a boundary this pass rewrote is exactly the kind a mesher can drop a piece of.
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
        (
            nacre_props::mass_props(&m, out[0]).unwrap().volume,
            faces_of(&m, out[0]),
        )
    };

    // ① A boss standing on the top face, cut: the plate, untouched — the planar twin's answer.
    assert_eq!(run([2.0, 2.0, 2.0], 1.0, BoolKind::Cut), (32.0, 6));
    // ② A boss sunk to the top face and fused: it adds nothing, so likewise.
    assert_eq!(run([2.0, 2.0, 1.0], 1.0, BoolKind::Fuse), (32.0, 6));
    // ③ A boss through the plate whose cap is flush with the bottom, fused. The bottom is **one**
    //    face (plate bottom + the boss's cap): 6 plate faces, of which the bottom absorbed the
    //    disk, plus the lateral and the boss's top cap = 8. This is the seam the user could see.
    let (v, f) = run([2.0, 2.0, 0.0], 3.0, BoolKind::Fuse);
    assert_eq!(f, 8, "the bottom is one face");
    // The boss spans z ∈ [0, 3] and the plate z ∈ [0, 2], so only its last millimetre of height
    // adds material.
    assert!(
        (v - (32.0 + std::f64::consts::PI * 0.25)).abs() < 1e-12,
        "{v}"
    );

    // Negative controls — a circle that separates *material from void* must survive.
    // ④ A through bore: its circle is shared with the cylinder, which is not in the plane's group,
    //    so nothing joins and the hole stays a hole.
    let (v, f) = run([2.0, 2.0, -1.0], 4.0, BoolKind::Cut);
    assert_eq!(f, 7, "plate faces with two mouths, plus the bore's lateral");
    assert!(
        (v - (32.0 - std::f64::consts::PI * 0.25 * 2.0)).abs() < 1e-12,
        "{v}"
    );
    // ⑤ A boss standing on the top face, **fused**: the contact disk is interior and is not a face,
    // so this pass has nothing to do and the answer must not move.
    let (v, f) = run([2.0, 2.0, 2.0], 1.0, BoolKind::Fuse);
    assert_eq!(f, 8, "annulus + 5 plate faces + lateral + cap");
    assert!(
        (v - (32.0 + std::f64::consts::PI * 0.25)).abs() < 1e-12,
        "{v}"
    );
}

/// **A wall whose plane passes near a hole, on a body that moved** — the face-level clearance
/// test. The infinite plane `y = 47.8` clears the bore at `(17.1, 48.6)` by 0.8 with `r = 2.12`,
/// so the cheap test fails and each face on that class has to answer for itself; the pocket wall
/// that carries it sits at `x ∈ [60.3, 73.2]`, nowhere near the bore, and says so. The user's
/// 12-up array is this shape. It runs on both roads a move takes. Plain, the move is exact on every
/// statement and carried: the moved cell is world-stated and its corners and axis share the world.
/// Notched (`fixtures::through_notched` — a `Through` face), the move is recorded: the wall's
/// corners are stated before the chain and the gate carries them out through its fold to meet the
/// axis, and the moved bore is a world cylinder only through the fold (`world_cylinder_def`).
#[test]
fn a_moved_face_answers_the_clearance_test() {
    for recorded in [false, true] {
        moved_face_answers_the_clearance_test(recorded);
    }
}

fn moved_face_answers_the_clearance_test(recorded: bool) {
    use nacre_exact::Rat;
    let cell = |m: &mut Model| {
        let plate = crate::fixtures::cuboid(
            m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([86.0, 86.0, 71.5]),
        );
        let pocket = crate::fixtures::cuboid(
            m,
            Point3::from_array([60.3, 47.8, 0.0]),
            Point3::from_array([73.2, 64.5, 68.5]),
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
        m.rebuild_adjacency();
        let bore = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([17.1, 48.6, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            2.12,
            80.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        if recorded {
            super::fixtures::through_notched(m, out, 10.0)
        } else {
            out
        }
    };
    let mut m = Model::new();
    let a = cell(&mut m);
    let b = cell(&mut m);
    let b = transform(
        &mut m,
        b,
        &nacre_exact::Isometry::translation([
            Rat::try_from_f64(100.3).unwrap(),
            Rat::from_int(0),
            Rat::from_int(0),
        ]),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert_eq!(carries_motion(&m, b), recorded, "recorded {recorded}");
    if recorded {
        // The population, or both roads go unmeasured: the bore stands on a chain that folds, and
        // the near-miss wall's corners answer in the chain's frame, not the world's.
        let faces = m.shell(m.solid(b).outer).faces.clone();
        let lateral = faces
            .iter()
            .map(|&f| m.face(f).surface)
            .find(|&s| matches!(m.surface(s), nacre_topo::Surface::Cylinder { .. }))
            .expect("the bore");
        let leaf = m
            .surface_motion(lateral)
            .expect("the bore records the move");
        assert!(
            m.chain_folds(leaf) && m.world_cylinder_def(lateral).is_some(),
            "the bore stands on a chain that folds to a world statement"
        );
        let corners = |f: Handle<Face>| -> Vec<Handle<nacre_topo::Vertex>> {
            m.face(f)
                .outer
                .half_edges
                .iter()
                .map(|he| m.he_start(*he))
                .collect()
        };
        let wall = faces
            .iter()
            .copied()
            .find(|&f| {
                corners(f)
                    .iter()
                    .all(|&v| m.vertex_point(v).as_array()[1] == 47.8)
            })
            .expect("the pocket wall on y = 47.8");
        assert!(
            corners(wall)
                .iter()
                .all(|&v| matches!(m.vertex_meet(v), Some((_, Some(_))))),
            "the near-miss wall's corners are stated before the chain"
        );
    }
    // The cells stand clear of each other (a *touching* pair is the contact family, another
    // cell), so this is an ordinary two-body result — what it pins is that the gate **decided**
    // the near-miss wall face by face.
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved cell is judged");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let pi = std::f64::consts::PI;
    let notch = if recorded { 10.0 * 298.0 / 75.0 } else { 0.0 };
    let one = 86.0 * 86.0 * 71.5 - 12.9 * 16.7 * 68.5 - pi * 2.12 * 2.12 * 71.5 - notch;
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!(
        (v - 2.0 * one).abs() <= 1e-9 * one,
        "recorded {recorded}: {v} vs {}",
        2.0 * one
    );
}

/// The negative control: a **rotated** cylinder has no exact world description, and the gate
/// says so by its own name rather than measuring across two frames. (A 90°-family turn keeps a
/// datum exact and records nothing, so the angle here is one that does record.)
#[test]
fn a_rotated_cylinder_is_still_undecided() {
    use nacre_exact::{Angle, Rat, Rotation};
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([18.0, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.1,
        30.0,
    )
    .solid;
    m.rebuild_adjacency();
    let tool = transform(
        &mut m,
        tool,
        &nacre_exact::Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(20), Rat::from_int(20), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(31)).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(carries_motion(&m, tool), "the turn records a chain");
    let live = m.live_solids().to_vec();
    let err = boolean(&mut m, BoolKind::Cut, plate, tool).expect_err("no world description");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            }
        ),
        "{err:?}"
    );
    assert_eq!(
        m.live_solids().to_vec(),
        live,
        "the live set survives the refusal"
    );
}
