//! Cylinders in the arrangement: extent and order on a class, hole arcs, class-space names, disks.

use super::*;

/// ★★★★ **Two cylinder operations in a row — measured for the first time.**
///
/// The corpus had **no chained-cylinder fixture at all**: every boss and bore stood on a plain
/// plate. So "a plate takes a second cylinder" was neither locked nor known, and it turned out to
/// be *half* true — this is the half that always built, kept as its own regression guard.
///
/// A bore or a boss standing on the plate's face leaves the operand's rings plane-named (a hole is
/// a full circle, the one curved loop the tracer's road speaks), so the next cylinder is business
/// as usual. A boss standing on a **wall** did not, for six wall-names running; it builds now, and
/// [`a_chained_cylinder_bounded_by_the_first_builds`] holds that half.
///
/// The volumes are derived, not copied: the plate is 96, a bore through it removes `π/4·2`, a boss
/// on top adds `π/4·1`, and a wall boss adds `π/4·4` with half its section buried (`−½·π/4·2`).
#[test]
fn chained_cylinder_operations_that_build_today_still_build() {
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
    // ★ A wall boss as a *middle* operation used to be the one that could not follow. It builds
    // now; its volumes live in `a_chained_cylinder_bounded_by_the_first_builds`, beside the
    // mechanism that opened it.
}

/// ★★★★★ **The gate decides a boss-seated wall, and records the pair — and nothing else in the
/// crate would notice if it stopped.**
///
/// The clearance test's whole content is a `bool`, and on a chained operand its answer travels
/// exactly one place: a `false` sends the seated pair to the `d = 0` record-and-pass arm, which
/// puts it in `crossings`, which is what makes the tracer's ruling and chord arms fire at all. A
/// wrong `true` would refuse nothing, panic nowhere and change no output — today's population
/// declines a step later either way (`PierceNode`, the lock below) — so **the suite would stay
/// green with the pierce corner's answer thrown away entirely**. Measured, by throwing it away:
/// 316 green. This is the lock that sees it.
///
/// ★★ **What it sees is that the corner is *readable*, not what it says** — recorded rather than
/// claimed. Forcing every pierce corner to one strip side leaves even this lock green, because the
/// verdict here is settled on the *other* separating axis: the wall's plate corners already sit
/// inside the boss's axial span, so `along` fails whatever the strip half answers. The corner's
/// side becomes decisive only for a face that clears along the axis and has to be judged across
/// it, and no fixture builds one yet. Removing the read itself is what turns this red
/// (`CurvedOperandBoundary`, measured).
#[test]
fn the_gate_records_a_wall_the_boss_is_seated_on() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array([6.0, 2.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    // ★ The gate *deciding* is half the claim — it used to stop at the first corner the boss made.
    let setup = crate::planes::plane_index_setup(&m, first, b).expect("the gate decides");
    // ★★ And recording is the other half. The boss sits **on** `y = 0`, so those faces are not
    // clear of it; that is what the record-and-pass arm is for, and a missing entry here means the
    // clearance test answered "clear" about a wall a cylinder is standing in.
    //
    // ★ The pair is checked by **what it means**, not by its indices: the arm records exactly the
    // classes whose plane the cylinder's axis lies *on*, so that is the assertion — renumbering
    // the classes cannot make it pass or fail for the wrong reason.
    assert_eq!(setup.crossings.len(), 1, "one seated pair, not more");
    for &(c, ci) in &setup.crossings {
        let coeffs = setup.geom[c].world_rat.expect("a named wall class");
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

/// ★★★★★ **The rule that replaced the two plane fences is *differenced* against them, not argued
/// equal to them.**
///
/// They are the same predicate wherever the fence is well posed — the fence plane is the class that
/// pins one end, so it crosses the line exactly at that end and "the side the far endpoint is on" is
/// the half-line from there. That argument has a hole: where the far endpoint sits *on* the fence
/// plane, `want` is `0` and every nonzero side reads as outside. So the arc split computes both
/// verdicts for every crossing it examines and counts the disagreements, and this reads the counter.
///
/// ★★★★★ **Three negative controls, and one of them says the green is thinner than it looks.**
/// ☑ Inverting `closed_contains`' verdict inside the probe: **red, 8 of 8** — both roads are really
/// read. ☑ Letting the fence side check only *one* of its two fences: **still green**. ☑ Letting it
/// check **neither**, so it answers "inside" unconditionally: **still green** — which says that on
/// this fixture *no crossing is out of extent*, so the two predicates agree by both saying yes. The
/// case that would separate them — a segment ending strictly inside the disk, whose line's far
/// crossing lies past its end — is what a **chained** operand produces, and that population is still
/// behind the gate. So: this locks that the roads are wired to the same question, and the census is
/// what locks the answers.
///
/// ★ The disagreement half is asserted **globally**, not as a delta: any fixture anywhere in this
/// binary that does reach an out-of-extent crossing has to agree too, whatever order it runs in.
///
/// ★ It expires with the gate: the fences need both endpoints' rational coordinates, which is the
/// demand the next rung removes.
#[test]
fn the_extent_rule_agrees_with_the_fences_it_replaced() {
    use crate::arrangement::extent_probe::{ASKED, DISAGREED};
    use core::sync::atomic::Ordering::Relaxed;
    let before = (ASKED.load(Relaxed), DISAGREED.load(Relaxed));
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
    let after = (ASKED.load(Relaxed), DISAGREED.load(Relaxed));
    assert!(
        after.0 > before.0,
        "the probe never ran, so it measured nothing: asked {} -> {}",
        before.0,
        after.0
    );
    assert_eq!(
        after.1, 0,
        "the line-order rule and the plane fences disagreed on {} of the {} crossings this binary \
         has examined",
        after.1, after.0
    );
}

/// ★★★★★ **The order rule may read the retired ruler backwards, but it may never reshuffle it.**
///
/// The arc split used to sort its points by their parameter on `pierce_meet`'s canonical meet line;
/// it asks [`combinatorics::order_located`] now. The two do **not** agree pointwise — the ruler's
/// direction comes from the classes' rational coefficients and the rule's axis sign from the judge's
/// stored planes, and those two spellings name the same plane without naming the same side. A
/// wholesale reversal is harmless: the pieces come out in the comparator's own ascending order and
/// `forward` is taken with that same comparator, so the two cancel. A **partial** disagreement would
/// not cancel — it would put sub-segments between the wrong pairs of points — and the bit census
/// cannot see it, because it sorts before it hashes.
///
/// So the proposition is per segment: **the same order, or exactly its reverse, and never in
/// between.** ☑ Flipping one comparison inside the probe turns this red.
#[test]
fn the_order_rule_never_reshuffles_the_ruler_it_replaced() {
    use crate::arrangement::order_probe::{
        EQUALITY_DISAGREED, REVERSED, SCRAMBLED, SEGMENTS, WC_ABOVE_WALL, WC_BELOW_WALL,
    };
    use core::sync::atomic::Ordering::Relaxed;
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
    assert!(
        SEGMENTS.load(Relaxed) > 0,
        "the probe never ran, so it measured nothing"
    );
    assert_eq!(
        SCRAMBLED.load(Relaxed),
        0,
        "{} of the {} segments this binary split came out neither in the ruler's order nor exactly \
         reversed ({} were reversed)",
        SCRAMBLED.load(Relaxed),
        SEGMENTS.load(Relaxed),
        REVERSED.load(Relaxed)
    );
    assert_eq!(
        EQUALITY_DISAGREED.load(Relaxed),
        0,
        "the two roads disagreed {} times about whether two split points are the same place — the \
         `CoincidentNodes` refusal is reachable from one and not the other",
        EQUALITY_DISAGREED.load(Relaxed)
    );
    // ★★★ **Both plane orders are exercised, and the invariant above is why neither is "the right
    // one".** ☑ Measured over the whole binary: **36 of 256** segments have `wc > wall`, where the
    // sorted pair and the call-order pair are opposite calls. Swapping the pair flips *every*
    // comparison on a segment, so it can only turn "same" into "reversed" — which the assertion
    // above already says is harmless. The sorted pair is chosen for agreeing with
    // `NodeId::Pierce`'s convention, not because a fixture prefers it.
    //
    // ★ That count is **recorded, not asserted**: these counters accumulate across the binary, and
    // this test cannot know what has run before it. Only the two "never" claims above are safe to
    // assert from here. The relation is read off `wc` and `wall`, which nothing in this cell moves.
    let _ = (WC_ABOVE_WALL.load(Relaxed), WC_BELOW_WALL.load(Relaxed));
}

/// ★★★★★ **A holed lateral's ruling comes out in three pieces, and the middle one grazes.**
///
/// A wall through the boss's axis meets its lateral in two rulings, and where the boss is buried in
/// the plate the lateral does not *cross* that wall — it ends at it. So the mark is
/// `Transversal · Graze · Transversal` along the axis, not one full-height crossing.
///
/// ★★ **Held at the trace, because there is no result to open**: the operation this exercises still
/// refuses further down (the assembly's own incompleteness), so the face count that would show the
/// ghost wall face gone cannot be read. The trace's own answer is what survives, the way
/// `arrangement`'s `sides == [-1, 1]` lock already holds the hole-free case.
///
/// ★ The claim is over **every** carved ruling this binary produces, in whatever order the tests
/// ran — all of them come from the wall-boss family.
#[test]
fn a_holed_laterals_ruling_grazes_where_the_hole_is() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array([6.0, 2.0, -1.0]), up, 0.5, 4.0);
    m.rebuild_adjacency();
    let _ = boolean(&mut m, BoolKind::Cut, first, b);
    // ★ The ledger is the binary's: a panel's ruling (one graze) and a chain's (a
    // transversal then a graze) are recorded too. Read this fixture's boss — origin `(2, 0, −1)`,
    // rulings spanning `t ∈ [0, 4]` — and nothing else.
    let carved: Vec<Vec<crate::arrangement::SegKind>> = crate::arrangement::ruling_probe::CARVED
        .lock()
        .expect("the probe's lock is never held across a panic")
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
    // ★ Both rulings are carved, not just one: the hole has a vertical edge on each.
    assert!(carved.len() >= 2, "only {} ruling carved", carved.len());
    // ★★★★★ **And which side the face occupies, against the fixture's own geometry.** The buried
    // half of the boss is the `y > 0` one, so along the hole's vertical edges the lateral survives
    // at `y < 0` — the stored-normal side exactly when that normal points at `−y`. ☑ Flipping any
    // factor of the derivation turns this red; nothing downstream does, yet.
    let sides: Vec<(bool, f64)> = crate::arrangement::ruling_probe::GRAZE_SIDE
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        // ★ Only the **hole's** runs: the same boss's Cut result (the census's) puts a *panel* on
        // the same wall and span, and its face lies on the other side — measured (`body_above`
        // false against the hole's true), which is exactly what a filter by origin alone let in.
        .filter(|g| {
            g.kind == combinatorics::CycleKind::Hole
                && g.origin == [2.0, 0.0, -1.0]
                && g.span[1] - g.span[0] == 2.0
        })
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

/// ★★★★★ **Which of the two arcs is the hole — measured, not argued.**
///
/// A plane cuts a circle in two points, and the whole difficulty of a lateral face's hole is which
/// of the two arcs it is: the hole's two ruling edges lie on **one** plane here, so that plane's
/// sides cannot tell them apart. `cycle_on_class` answers from the ring's own winding (material is
/// on the left of the ring's travel), and gets there through three signs — the face's
/// `orient_sign`, the class's `frame_sign`, and whether the class's stored normal points along the
/// axis. Reversed, the answer is not approximate but **exactly opposite**: the graze would sit on
/// the band and the crossing inside the hole.
///
/// ★★ **Nothing downstream notices yet** — these fixtures stop at `loop_winding` before
/// `label_cells` could judge a label — so the arc is realized and judged here instead. The kernel
/// reads no coordinate to choose it; this reads one afterwards.
///
/// ☑ **Which factors this actually sees.** Turning each one off in turn: reversing either arm's
/// winding is **red**, dropping `plus_t_is_above` is **red**, and dropping the root restatement
/// (the ruling named for a different ⊥ plane) is **red**. Dropping `frame_sign` or the face's
/// `orient_sign` is **green** — both are `+1` everywhere in today's holed population (an
/// `add_cuboid` face is never `Reversed`, and a boss's lateral faces outward), so this lock cannot
/// see them and their reasons stand on the derivation alone. A **bore** with a hole in its lateral
/// is what would exercise them.
///
/// ★ The claim is over **every** arc this binary carves, in whatever order the tests ran: all of
/// them come from the wall-boss family, where the plate stands on `y > 0` and the boss is centred
/// on the wall `y = 0`, so the buried half — the hole — is the `y > 0` one. A fixture that buries
/// a lateral somewhere else would need its own reading of "which side is the hole", and this
/// assertion is where that would show up.
#[test]
fn a_holes_arcs_run_the_way_the_hole_lies() {
    // ★ Two second operands, because they reach **different arms** of the walk. A cylinder rising
    // from `z = −1` puts its cap on `z = 0`, which is the hole's own low rim — an on-line *run*.
    // One rising from `z = 1` puts a cap **strictly inside** the hole's `z ∈ (0, 2)`, where the
    // ring crosses the class on its two ruling edges instead — the *crossing* arm, whose extra
    // step is restating the ruling's root for a different ⊥ plane.
    for base in [-1.0, 1.0] {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([12.0, 4.0, 2.0]),
        );
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let a = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
        m.rebuild_adjacency();
        let first = boolean(&mut m, BoolKind::Fuse, plate, a).expect("the wall boss builds")[0];
        m.rebuild_adjacency();
        let b = m.add_cylinder(Point3::from_array([6.0, 2.0, base]), up, 0.5, 4.0);
        m.rebuild_adjacency();
        let _ = boolean(&mut m, BoolKind::Cut, first, b);
    }
    // ★ The ledger is the binary's, not this test's: every fixture that carves a cycle writes to
    // it, and a panel's and a chain's arcs (which face any way) are carved too. Read
    // only this fixture's boss (`(2, 0, −1)`) and only its **hole** — the sentence is about holes.
    let mids: Vec<[f64; 3]> = crate::arrangement::arc_probe::MIDS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .iter()
        .filter(|m| m.kind == combinatorics::CycleKind::Hole && m.origin == [2.0, 0.0, -1.0])
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

/// ★★★★★ **A boolean's own result takes a second cylinder — the whole road, end to end.**
///
/// A boss standing on a **wall** was the one chained operand that did not build. It buries half the
/// boss's lateral in the plate, so that face comes back as a band with a **hole**, and every layer
/// downstream had a sentence that was false about it. The wall moved six times as those were
/// closed one at a time — `TraceDeclined { PierceNode }`, `WitnessNotRational`,
/// [`RejectReason::CurvedStraightRun`], [`RejectReason::OpenResultShell`],
/// [`RejectReason::LabelConflict`] when one of two cancelling falsehoods was fixed, then
/// `OpenResultShell` again with the arrangement whole — and this lock was the map of that walk.
/// It is now the map of the far side.
///
/// ★★★★★ **The last false sentence was that a label answers two questions.** A label says where
/// *material* is; the panel road also read it as saying whether the lateral face is **there**. On
/// the first fusion the two sectors of the holed rim carry different labels (the buried one is
/// inside the *other* solid) and the road happened to be right; on a second operation the plate is
/// already own material, the two sectors' labels are **literally identical**, and both survived —
/// filling the hole in and leaving its four boundary edges claimed once each. The trace had said
/// which sector was a face all along, per arc, since the hole taught it to mark by angular extent:
/// `Graze` where the face ends, `Transversal` where it runs through.
/// `cyl_chart::read_cell`'s `face_spans` asks it.
///
/// ☑ Measured across these four fixtures: **four** sectors dropped for existence, exactly one per
/// operation — the buried half — and the volumes below are what says that count is right.
///
/// ★ The volumes are derived, not copied: the plate is 96, the wall boss adds `π/4·4` with half its
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
    // ★★★ **And that the existence gate is what did it.** The chart's census holds it where the
    // cells are read:
    // `exist_marks_false` (a sector dropped because the face is not there — the four above, held
    // as growth ≥ 4 in `cyl_chart::tests::the_cells_read_their_chamber_from_the_horizontal_lines`)
    // and `arcs_no_mark`/`arcs_multi_mark` (every cut end read carries exactly one lateral mark of
    // its own solid — the doc that says why the `Seated` skip is load-bearing lives on `face_spans`).
}

/// ★★★★★ **A boolean's result carries names the tracer's plane-only road cannot read — and it
/// says so instead of aborting.**
///
/// A boss standing on a wall leaves the plate's caps bitten by an **arc** and the wall split in two
/// by the boss's **rulings**, so those faces' rings run along a cylinder. [`combinatorics`]'s ring
/// naming asks two plane-only questions of every edge — the carried wall class and the vertex's
/// three-plane name — and both go through `ClassIx::plane`, which **panics** on a cylinder. That
/// accessor is right to: for its forty-odd other callers a cylinder there is an upstream filter
/// bug. This is the caller whose input is an *operand*, so this is where the filter belongs.
///
/// ★ **Measured before the filter existed: this very call panicked** ("a plane-only path got
/// cylinder class 0"). The population gate refuses such an operand long before the tracer sees it,
/// so nothing in production reaches this today — which is exactly why the lock calls the function
/// **directly**, the same way the ring-naming goldens next door do. The gate opens in the cell that
/// builds the road; this is the net that has to be under it first.
///
/// The untouched walls are the negative control: they are still named, so the filter is a filter
/// and not a blanket refusal.
#[test]
fn an_operand_bounded_by_a_cylinder_is_named_in_class_space() {
    // ★★ **The corner boss is here because the wall boss is symmetric.** Its axis lies *on* the
    // wall, so the two rulings sit either side of that plane and a root read the wrong way round
    // lands on a mirror-image point — a lock that only ever saw this fixture could pass with the
    // restatement inverted. The corner's two walls break that symmetry.
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for at in [[2.0, 0.0, -1.0], [12.0, 0.0, -1.0]] {
        let (u, s) = named_in_class_space(at);
        ups.extend(u);
        sides.extend(s);
    }
    // ★★ **The direction assertions inside are only worth their ink if the fixtures move them** —
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
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array(at),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
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
            let def = crate::planes::world_cylinder_def(&m, surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(surf) else {
                unreachable!("a cylinder class names a cylinder")
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                realized: *cache,
                owner: crate::planes::SolidSide::A,
            }
        })
        .collect();

    let jd = crate::planes::test_judge(&planes);
    let (mut curved_walls, mut pierce_names) = (0usize, 0usize);
    let (mut ups, mut sides): (Vec<bool>, Vec<i8>) = (Vec::new(), Vec::new());
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        if matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            continue; // the tracer skips a lateral face; its loops are the rims
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
                crate::combinatorics::Wall::Arc { .. } => {
                    assert!(
                        !straight,
                        "face {:?} edge {i}: an arc carrier on a straight curve",
                        fh.index()
                    );
                    curved_walls += 1;
                }
                // ★★★★★ **The direction bits, against the walk and against the geometry.** The
                // carrier alone says *which surface*; `up` and `side` say *which way* and *which
                // of the two rulings*, and a consumer that reads them wrong builds a face wound
                // backwards or seated on the far side of the cylinder. Both are measured off
                // realized coordinates on purpose — the code derives them without any (`side`
                // exactly through `quad::plane_side`, `up` from the cutting planes' axial
                // parameters), so the coordinate is a genuinely second road to the same bit.
                crate::combinatorics::Wall::Ruling { cyl: k, side, up } => {
                    assert!(
                        straight,
                        "face {:?} edge {i}: a ruling carrier on a curved curve",
                        fh.index()
                    );
                    let crate::planes::ClassIx::Plane(near) = plane_ix[fp] else {
                        unreachable!("a lateral face was skipped above")
                    };
                    let axis = cyls[k].realized.axis();
                    let axial = |p: Point3| (p - axis.origin()).dot(axis.direction());
                    assert_eq!(
                        up,
                        axial(m.vertex_point(corner((i + 1) % n)))
                            > axial(m.vertex_point(corner(i))),
                        "face {:?} edge {i}: `up` disagrees with the walk",
                        fh.index()
                    );
                    // `side` is the sign of `(x − o) · (m × n̂)`. ★★ **`n̂` has to be the class's
                    // `world_rat` normal, and no other spelling of the class's plane will do** —
                    // measured, by writing `plane.normal()` here first and watching it come out
                    // opposed (`[0,-1,0]` against `[0,1,0]`, and `frame_sign` is `+1`, so that is
                    // not the reconciliation either). A class has no outward normal to agree on:
                    // `world_rat` is the plane's *name*, whose sign is its own, while `plane` is a
                    // stored surface's. `side` is a **label** telling the two rulings apart, so it
                    // is well defined exactly as long as every road spells `n̂` the one way.
                    // Sharing that input leaves the two roads independent where it counts — the
                    // code decides in exact `quad::plane_side` on the pierce meet, this in `f64`
                    // on the realized vertex.
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
            // ★★★★★ **Two roads, one point.** The name is this boolean's classes; the coordinate is
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
    // Four faces run along the boss — the two plate caps it bit an arc out of, and the two halves
    // its rulings split the wall into — and each contributes two curved edges and two pierce
    // corners. The boss's own caps are full circles (the one curved loop this road already spoke)
    // and the three untouched walls are plain.
    assert_eq!(curved_walls, 4, "edges riding the boss");
    assert_eq!(pierce_names, 8, "corners named as pierce points");
    (ups, sides)
}

/// ★★★★★ **A cylinder-pinned end is ordered, not refused — and the sense is the geometry's.**
///
/// `edge_dir` is the one place a direction is made, and its cylinder arm used to say `RingNaming`
/// by name: a pierce point has no third plane, and the integer predicate wants one. Measured
/// before this landed — every one of the 40 pins below came back refused. Now they order through
/// the `a + b√c` tower, and this fixes what the answer must be.
///
/// ★★★ **The oracle is a second road, and only the half that matters is second.** The direction
/// `n_p × n_q` is read from the same raw coefficients the code reads — deliberately, the way the
/// ruling lock next door shares `world_rat`: a *label*'s reference frame has to be one spelling or
/// the two roads are not comparing the same thing. What is independent is the part under test —
/// the **order of the two points on that line**: the code decides it exactly (a rational meet
/// against a pierce root, or two roots against each other), the oracle realizes both points in
/// `f64` and subtracts. Measured `|t|` from 1.5 to 456, so nothing here is decided in the noise.
///
/// ★★ **The population is one-sided and that is a fact, not a blind instrument.** All 40 order
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
    // ★ The count is the lock on the *population*: let the fixtures stop producing pierce-pinned
    // ends and every assertion below would pass vacuously. **48**, which counts the rings split at
    // a seam joint — the `[6, 4, −1]` boss sits on the +y wall, so
    // its bite on the plate's caps wraps the seam (θ = 0 is at −y), and those rings' ends join.
    assert_eq!(pins, 48, "cylinder-pinned ring ends exercised");
}

/// ★★★ **A circle can be an interior boundary, and then the merge erases it.**
///
/// A tool that only *touches* removes nothing, and the planar engine has always said so plainly:
/// a plate cut by a box resting on its top face comes back as the plate — same volume, **six
/// faces**, no imprint left behind (measured, both for face contact and for a box straddling an
/// edge). That measurement is the oracle here; there is no sentence in the design documents that
/// decides it.
///
/// The cylinder twin used to refuse. Its result was already right — exact volume, `validate`
/// clean — but the top plane came back as **two** faces: the plate's top with a circular hole, and
/// the disk filling it. Nothing joined them, because a disk's boundary is a `Bound::Circle` with
/// no nodes at all, so the merge's edge table could not see it; the corners left on that circle
/// went on naming the tool's cylinder after every face of it was gone, and the assembly refused
/// the whole boolean (`VertexNamesAbsentSurface`).
///
/// What connects them is that the other face holds **the same cylinder class** as a hole — and one
/// plane class with one cylinder class names one circle, so that coincidence *is* adjacency.
#[test]
fn a_disk_merges_into_the_face_it_lies_in() {
    let plate_and_boss = |m: &mut Model, base: [f64; 3], h: f64| {
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
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
        // ★ The merged face has to survive tessellation too: `validate` reads the topology store,
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

    // ★ Negative controls — a circle that separates *material from void* must survive.
    // ④ A through bore: its circle is shared with the cylinder, which is not in the plane's group,
    //    so nothing joins and the hole stays a hole.
    let (v, f) = run([2.0, 2.0, -1.0], 4.0, BoolKind::Cut);
    assert_eq!(f, 7, "plate faces with two mouths, plus the bore's lateral");
    assert!(
        (v - (32.0 - std::f64::consts::PI * 0.25 * 2.0)).abs() < 1e-12,
        "{v}"
    );
    // ⑤ A boss standing on the top face, **fused**: the contact disk is interior and never was a
    //    face, so this pass has nothing to do and the answer must not move.
    let (v, f) = run([2.0, 2.0, 2.0], 1.0, BoolKind::Fuse);
    assert_eq!(f, 8, "annulus + 5 plate faces + lateral + cap");
    assert!(
        (v - (32.0 + std::f64::consts::PI * 0.25)).abs() < 1e-12,
        "{v}"
    );
}

/// ★ **A wall whose plane passes near a hole, on a body that moved** — the face-level clearance
/// test's own frame question. The infinite plane `y = 47.8` clears the bore at `(17.1, 48.6)`
/// by 0.8 with `r = 2.12`, so the cheap test fails and each face on that class has to answer for
/// itself; the pocket wall that carries it sits at `x ∈ [60.3, 73.2]`, nowhere near the bore, and
/// says so — but only if its corners are read in the same frame as the axis. The user's 12-up
/// array is this shape, and it stopped here after the class descriptions were carried out.
#[test]
fn a_moved_face_answers_the_clearance_test() {
    use nacre_exact::Rat;
    let cell = |m: &mut Model| {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([86.0, 86.0, 71.5]),
        );
        let pocket = m.add_cuboid(
            Point3::from_array([60.3, 47.8, 0.0]),
            Point3::from_array([73.2, 64.5, 68.5]),
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
        m.rebuild_adjacency();
        let bore = m.add_cylinder(
            Point3::from_array([17.1, 48.6, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.12,
            80.0,
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
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
    assert!(carries_motion(&m, b), "the move records a chain");
    // The cells stand clear of each other (a *touching* pair is the contact family, another
    // cell), so this is an ordinary two-body result — what it pins is that the gate **decided**
    // at all, rather than declining because a corner was stated in another frame.
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved cell is judged");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let pi = std::f64::consts::PI;
    let one = 86.0 * 86.0 * 71.5 - 12.9 * 16.7 * 68.5 - pi * 2.12 * 2.12 * 71.5;
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
        .sum();
    assert!((v - 2.0 * one).abs() <= 1e-9 * one, "{v} vs {}", 2.0 * one);
}

/// ★ The negative control: a **rotated** cylinder has no exact world description, and the gate
/// says so by its own name rather than measuring across two frames. (A 90°-family turn keeps a
/// datum exact and records nothing, so the angle here is one that does record.)
#[test]
fn a_rotated_cylinder_is_still_undecided() {
    use nacre_exact::{Angle, Rat, Rotation};
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = m.add_cylinder(
        Point3::from_array([18.0, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.1,
        30.0,
    );
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
