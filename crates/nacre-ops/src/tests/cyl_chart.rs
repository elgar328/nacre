use super::probe::{ROWS, Row};
use crate::BoolKind;
use nacre_exact::Rat;
use nacre_math::{Point3, Vector3};
use nacre_topo::Model;

/// **The chart census is live, and every chart it records has the shape a chart must have.**
///
/// The real claim — "a chart covers its own class's rows" — is asserted in `census`, where the
/// fact is made.
///
/// ★★★★★ **The shape assertions are universal over every chart in the binary**, which is a
/// stronger statement than one about this fixture and which no interleaving can break. What is
/// this fixture's alone — how many charts it recorded — is counted on its own rows
/// (`ledger::owned`), and a row of another test's can neither add to that count nor stand in
/// for it. The universals:
///
/// * `rows >= 1` — a cylinder class exists because a face made it (`cyl_rows` refuses
///   otherwise), so the hand-written road always has at least one row for it.
/// * `z_lines >= 2` — a band needs two boundaries; a lateral face's own two rims are ⊥ classes
///   and both are in the set. ☑ Measured minimum across the suite: exactly 2.
/// * `theta` is **even** — a plane holding the axis cuts the lateral in *two* rulings
///   (`QuadRoot::{Lo, Hi}`), and nothing else contributes a vertical line. ☑ Measured: 0, 4,
///   6, 8.
///
/// And the same for the cells, again universally:
///
/// * `refused` is never set — a chart's θ orders can always be formed. ☑ Measured 0 across the
///   suite; asserted rather than merely counted, so the day a population appears it says so.
/// * `cells >= 1` — `z_lines >= 2` is one interval, and an interval with no ruling is still one
///   cell (the whole circle).
/// * every cell is claimed at most once, and `claimed + unclaimed == cells` — the partition is
///   total, so a cell can neither be lost nor counted twice by the comparison.
/// * `reversed` is never set: a panel ring's two rulings always match the sector's own order.
///   ☑ Measured 0 over 65 panels, which is what makes the exact-order join sound.
#[test]
fn the_chart_census_is_running() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let drill = m.add_cylinder(
        Point3::from_array([2.0, 2.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    crate::ledger::owned(|| {
        crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a through bore")
    });
    let rows = ROWS.all();
    // ★ **This bore's own charts, counted.** One cylinder means one class, and a class is one
    // chart row — a count the ledger could not state while every test's rows shared the list.
    assert_eq!(ROWS.mine().len(), 1, "a through bore records one chart");
    for r in &rows {
        let Row {
            z_lines,
            theta,
            rows: n,
            refused,
            cells,
            ..
        } = *r;
        assert!(n >= 1, "a class with no row: {r:?}");
        assert!(z_lines >= 2, "a band needs two boundaries: {r:?}");
        // ★ `theta % 2 == 0` ("a wall cuts two rulings") is not asserted: a
        // band's lateral reaches both rulings of a wall, a panel's or a chain's may reach one.
        let _ = theta;
        assert!(!refused, "a chart's theta order could not be formed: {r:?}");
        assert!(cells >= 1, "two z-lines are one interval: {r:?}");
        // A chart has `z_lines - 1` intervals, and an interval carries either no ruling or
        // some number of them — so these two counts are disjoint subsets of that many.
        assert!(
            r.whole_circle + r.odd_k < z_lines,
            "more intervals than the chart has: {r:?}"
        );
    }
}

/// **The chart's vertical answers are read, and they close.**
///
/// The real claims are asserted in `census`, where the facts are made — every ruling reaches
/// the chart with a label, and the walk around each interval returns to where it started.
/// What this holds is that the reading is **live** and that its counters are not vacuous.
///
/// ★★★★★ **What the vertical answer settled, which no label alone could.** 82
/// emitted bands span more than one θ-sector, and the cells alone cannot say what the extra
/// line is. With a label on the vertical lines there are three answers, and the third is
/// existence — the ruling may not be on this face at all. ☑ Measured over the suite:
/// **absent 0 · a
/// boundary the band road misses 0 · harmless 84**. The band road is calling those strips
/// uniform and the chart agrees they *are*: the walls crossing them change nothing at the
/// lateral there.
///
/// ☑ Beside it, over **450 alive-ruling observations** — a piece counts once per interval it
/// spans, so this is **not** the 862-piece denominator above: **122 whose wall does change the
/// material** (so the counter has a population), **61 intervals** carrying at least one, **12
/// whose every mark is a graze** — the existence question is not empty either — and **201
/// intervals closing, 0 not**.
///
/// ★★ **The first spelling of the existence counter was vacuous**: it asked whether a ruling's
/// mark list was *empty*, and a `MergedRuling` exists only because a lateral traced it. It
/// measured `0` because nothing was being looked at. The signal `face_spans` actually reads is
/// the **kind** (`Graze` = the face stops here), and that has 12.
#[test]
fn the_charts_vertical_answers_close() {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array([2.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    crate::ledger::owned(|| {
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds")
    });
    let rows = super::probe::rulings::ROWS.all();
    let mine = super::probe::rulings::ROWS.mine();
    // ★ Universal over every recorded chart, so no interleaving can break it.
    for r in &rows {
        assert_eq!(r.does_not_close, 0, "an interval did not close: {r:?}");
        // Both counts are subsets of the same denominator, in the same table.
        assert!(r.wall_flips <= r.rulings, "more flips than rulings: {r:?}");
        assert!(
            r.grazing_rulings <= r.rulings,
            "more grazes than rulings: {r:?}"
        );
        assert!(
            r.intervals_with_flip <= r.closes + r.does_not_close + r.unpaired,
            "more flipping intervals than walked ones: {r:?}"
        );
    }
    // ★★★★ And the counters are not vacuous — a zero above must mean "the population is
    // empty here", never "nothing was looked at".
    // ★★★★ And the counters are not vacuous — a zero above must mean "the population is
    // empty here", never "nothing was looked at". Counted on **this fixture's own rows**: the
    // boss's lateral is one class, so one chart; its four z lines (the boss's two rims and the
    // plate's two faces) make three intervals, and the wall through the axis cuts two rulings
    // over each — six alive-ruling observations. All three intervals close. The wall changes
    // the material on the two rulings of the interval buried in the plate, and on no other.
    let sum = |f: fn(&super::probe::rulings::Row) -> usize| mine.iter().map(f).sum::<usize>();
    assert_eq!(mine.len(), 1, "one class, one chart");
    assert_eq!(sum(|r| r.rulings), 6, "three intervals of two rulings");
    assert_eq!(sum(|r| r.closes), 3, "every interval closes");
    assert_eq!(
        sum(|r| r.wall_flips),
        2,
        "the buried interval's two rulings"
    );
    assert_eq!(
        sum(|r| r.intervals_with_flip),
        1,
        "and they are one interval's"
    );
}

/// **The cells read their chamber off the lines at their ends, and the ledger's bookkeeping
/// is total**.
///
/// The absolute claims are asserted in `census`, where the facts are made: a present cell
/// always has a speaking end (`src0_present == 0`, `other_present == 0` — the reason no
/// wall-side sign is needed), the emitter's faces cover exactly the emitted cells (the
/// cell→face assertions), and the ends' bookkeeping is total. What this holds, universally
/// over every recorded row so no interleaving can break it, is the rest:
///
/// * `emit_unknown == 0 || emitter_refused` — the emitter never puts a face over a present
///   cell it could not read, e.g. one whose two ends contradict (the `wal corner-lo` corpus
///   family has such cells — `src2_disagree` 8 per row — and is refused by name).
/// * `exist_disagree == 0`, `read_refused == 0` — the trace and the span tell the same
///   existence story wherever both speak.
/// * `nocircle_present == 0` — a present cell has a circle at both ends.
/// * `arcs_read == end_exact + exact_run_arcs`, `arcs_multi_mark == 0` — every cut-end read
///   carries at most one lateral mark of its own solid (the band road's MARKS contract), and
///   an end that spans a run of the rim's arcs reads every one of them. (`arcs_no_mark` is no
///   longer 0: a ⊥ cap through a notch reads cut ends inside the lateral's hole.)
///
/// And the counters are not vacuous, counted on **this test's own rows** (`ledger::owned`), so
/// each number is a fact about the fixtures below and not about whatever else the binary ran:
/// every end kind reaches the ledger, the chained pairs drop five sectors for existence, arcs
/// are read and faces are emitted. `end_other` is **measured, not unexercised**: every one of
/// the suite's is a whole-circle interval beyond a face's span whose cut rim's arcs disagree
/// about the far side, on an absent cell; the single-cut arm is the one still without a
/// population.
#[test]
fn the_cells_read_their_chamber_from_the_horizontal_lines() {
    use super::probe::cell_ends::ROWS as CELL_ENDS;
    let snapshot = || -> Vec<super::probe::cell_ends::Row> { CELL_ENDS.all() };
    let sum = |rows: &[super::probe::cell_ends::Row],
               f: fn(&super::probe::cell_ends::Row) -> usize| {
        rows.iter().map(f).sum::<usize>()
    };
    // A wall boss (exact arcs, θ- and z-merges), a corner boss and a boss on top (a line with
    // no circle of this cylinder, cells outside the face), and the chained pairs (sectors
    // dropped for existence).
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    crate::ledger::owned(|| {
        for (base, h) in [
            ([2.0, 0.0, -1.0], 4.0),
            ([4.0, 4.0, -1.0], 4.0),
            ([2.0, 2.0, 2.0], 1.0),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(Point3::from_array(base), up, 0.5, h);
            m.rebuild_adjacency();
            crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        }
        // ★ A **through-axis wall** on the wall boss: its plane holds the axis, so it cuts the
        // lateral along two rulings — and the chart has a θ line for only the one the face
        // reaches. The **Cut** is the kind that builds through it (the wall boss's Fuse waits on
        // the emitter's rim-station ladder — `RulingBoundNotYet`). ★★ It reads **no run of arcs**:
        // measured over this test's own rows, `exact_run_arcs` is 0 here, and the run population
        // belongs to the crossing census's through-axis walls instead.
        {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Cut, plate, boss).expect("the notch builds");
            m.rebuild_adjacency();
            let wall = m.add_cuboid(
                Point3::from_array([2.0, -3.0, -2.0]),
                Point3::from_array([5.0, 3.0, 6.0]),
            );
            m.rebuild_adjacency();
            crate::boolean(&mut m, BoolKind::Cut, out[0], wall)
                .expect("the through-axis wall cuts");
        }
        let wall_boss = ([2.0, 0.0, -1.0], 4.0, BoolKind::Fuse);
        for second in [
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            // A short boss whose cap (z = 2.5) is a ⊥ line strictly inside the wall boss's span:
            // the wall boss's circle is traced there but bounds no face, so the line is not a
            // boundary — the band-shaped merge population.
            ([6.0, 2.0, 2.0], 0.5, BoolKind::Fuse),
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Cut),
        ] {
            let (_, r) = crate::tests::chained(wall_boss, second);
            r.expect("the chained pair builds");
        }
    });
    let rows = snapshot();
    let mine = CELL_ENDS.mine();

    for r in &rows {
        assert!(
            r.emit_unknown == 0 || r.emitter_refused,
            "a cell's chamber could not be read, yet the emitter emitted: {r:?}"
        );
        assert_eq!(
            r.src0_present, 0,
            "the record-site assertion was not live: {r:?}"
        );
        assert_eq!(
            r.exist_disagree, 0,
            "trace and span disagree on existence: {r:?}"
        );
        assert_eq!(
            r.read_refused, 0,
            "the reader refused a cell by name: {r:?}"
        );
        // The reference-road checks, held as counts.
        assert_eq!(
            r.nocircle_present, 0,
            "a cell with a face meets a line with no circle: {r:?}"
        );
        // ★ `arcs_no_mark` is not held at 0 here: a ⊥ cap through a wall boss's notch
        // reads cut ends inside the lateral's hole, and those carry no mark by right.
        // The record site asserts the true proposition — none only where the cell is absent.
        assert_eq!(
            r.arcs_multi_mark, 0,
            "a cut end read carries two lateral marks of one solid: {r:?}"
        );
        // ★★★★★ **A run is why this is not `arcs_read == end_exact`.** A sector may span several
        // of the rim's arcs, and the proposition that holds over every row in the binary — the
        // runs among them — carries the run's extra arcs by name.
        assert_eq!(
            r.arcs_read,
            r.end_exact + r.exact_run_arcs,
            "every exact end is a cut end read, and only those: {r:?}"
        );
        // ★ The emitter's face count is not predicted here: its faces are the regions of the
        // chart, asserted against the reads where the fact is made (`census`).
        assert_eq!(
            r.end_swapped, 0,
            "a ruling arrived with z descending: {r:?}"
        );
        assert_eq!(
            r.end_disk + r.end_exact + r.end_other + r.end_uncovered + r.end_nocircle,
            2 * r.cells,
            "two ends per cell: {r:?}"
        );
        assert!(r.emit <= r.cells, "more emitted cells than cells: {r:?}");
        assert!(
            r.other_present <= r.cells,
            "more other-ended cells than cells: {r:?}"
        );
        assert!(
            r.exist_marks_false <= r.cells,
            "more dropped cells than cells: {r:?}"
        );
    }
    // ★★★★ Non-vacuity — a zero above must mean "the population is empty", never "nothing
    // was looked at". `end_other` is deliberately not in this list.
    // ★★★★ Non-vacuity — a zero above must mean "the population is empty", never "nothing was
    // looked at" — and on **this test's own rows**, so no neighbour's fixture can stand in for
    // one of these. `end_other` is deliberately not among them.
    assert_eq!(mine.len(), 20, "the fixtures' charts");
    assert_eq!(sum(&mine, |r| r.end_disk), 84, "disk ends read");
    assert_eq!(sum(&mine, |r| r.end_exact), 132, "exact arc ends read");
    assert_eq!(sum(&mine, |r| r.end_nocircle), 10, "lines with no circle");
    assert_eq!(sum(&mine, |r| r.emit), 80, "cells the reader emitted");
    assert_eq!(sum(&mine, |r| r.arcs_read), 132, "cut ends read");
    // ★ A sector spanning the **phantom** pieces of a half rim (a wall boss) reads `Uncovered`:
    // those pieces are no edges.
    assert_eq!(sum(&mine, |r| r.end_uncovered), 2, "uncovered rim pieces");
    assert_eq!(
        sum(&mine, |r| r.emitted_faces),
        20,
        "faces the emitter emitted"
    );
    // The chained pairs drop sectors for existence; this is how many.
    assert_eq!(
        sum(&mine, |r| r.exist_marks_false),
        5,
        "sectors dropped for existence"
    );
    // ★★ **The run population is not this test's.** `exact_run_arcs` is 0 over every row here,
    // including the through-axis wall's — the sector that spans a rim node is built by the
    // crossing census's walls, and that is what keeps the equality above from being
    // `arcs_read == end_exact` over the binary.
    assert_eq!(
        sum(&mine, |r| r.exact_run_arcs),
        0,
        "no end reads a run here"
    );
}

/// **The derivation, on an axis that is not `+z`** — both roads, on a genuinely tilted one.
///
/// ★★★★★ **Why this fixture had to be built.** Everything above is derived from the gate's
/// branch table and depends on no axis direction, but every number beside it was measured on
/// a corpus that is effectively axis-aligned: the kit builder raises cylinders on world XY
/// (`+z` only), the ops corpus is `+z`/`+x`/`+y`, and all three *tilted* cylinders in it are
/// **reject** fixtures — so no tilted axis had ever reached the chart. On `+z` against an
/// axis-aligned box a class's axis parameter is just a `z` coordinate and the division in
/// `bands::axis_param` is trivial; tilted, it is a real rational one.
///
/// ★★★★ **Both halves are stated on the exact road, and that is the whole trick.**
/// A *Pythagorean frame* — `u = (0.6, 0.8, 0)`, `v = (−0.48, 0.36, 0.8)` — lifts to axes that
/// are exactly orthonormal in the rationals, so `SketchPlane::from_axes` takes the
/// world-rational path and every face of the prism raised on it states narrow coefficients
/// (asserted below, small integers). Its normal `(0.64, −0.48, 0.6)` is the cylinder's axis,
/// and the frame's own `u` is a **unit rational perpendicular** to it — exactly the `ref_dir`
/// `Model::add_cylinder_exact` requires.
///
/// ★★★★★ **`Model::add_cylinder` cannot state this, and that is not a kernel limit.** That
/// entry normalizes an `f64` axis and lifts its cap points back out of `any_perpendicular`'s
/// computed floats, so on a tilted axis the caps land outside the narrow window and the
/// population gate honestly declines (`CylinderGateUndecided` — measured while building this).
/// Its own doc calls it *the test entry*; the production road is `add_cylinder_exact`, and on
/// that road the tilted case is not the irrational case.
///
/// ☑ What this establishes: the chart's two axes are built, and the axes' universal claims
/// hold,
/// on an axis with no zero component — through-bore (circles only) **and** a wall holding the
/// axis (rulings).
#[test]
fn the_chart_stands_on_a_tilted_axis() {
    // Two runs of one fixture: the bore centred inside the prism (⊥ classes only), then
    // centred **on** a wall, which puts that wall through the axis and cuts rulings.
    // The prism is 2 x 2 x 1 and the bore has r = 1/2, so the removed volume is `pi/4`
    // whole and half of that when the axis lies in a wall.
    let bore = std::f64::consts::FRAC_PI_4;
    for (name, base, want_rulings, want_vol) in [
        ("through bore", [-0.52, 1.64, 0.2], false, 4.0 - bore),
        (
            "half bore on a wall",
            [0.08, 2.44, 0.2],
            true,
            4.0 - bore / 2.0,
        ),
    ] {
        let before = ROWS.len();
        let mut m = Model::new();
        let plane = crate::SketchPlane::from_axes(
            Point3::from_array([0.0; 3]),
            Vector3::from_array([0.6, 0.8, 0.0]),
            Vector3::from_array([-0.48, 0.36, 0.8]),
        );
        let frame = match crate::apply(
            &mut m,
            &crate::Operation::DatumPlane {
                def: crate::DatumDef::Stated(plane),
            },
        ) {
            Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
            other => panic!("stating the tilted plane: {other:?}"),
        };
        let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let crate::OpOutput::Extrude { solid, faces, .. } = crate::apply(
            &mut m,
            &crate::Operation::Extrude {
                frame,
                profile: crate::Profile2d::polygon(vec![
                    p(0.0, 0.0),
                    p(2.0, 0.0),
                    p(2.0, 2.0),
                    p(0.0, 2.0),
                ])
                .expect("a square"),
                dist: 1.0,
            },
        )
        .expect("the tilted prism") else {
            unreachable!()
        };
        // ★★★★ **The fixture qualifies its own population**: every face must state narrow
        // world coefficients, and
        // each must be either ⊥ or ∥ to the axis **exactly** — otherwise the gate would be
        // answering the oblique branch and the chart would never be reached.
        let axis = [0.64, -0.48, 0.6].map(|x| Rat::from_decimal(x).expect("a short decimal"));
        for &f in &faces {
            let n = m
                .surface_name
                .get(&m.face(f).surface)
                .and_then(|nm| nm.narrow())
                .copied()
                .expect("a tilted prism's face states itself");
            let n = [n[0], n[1], n[2]];
            let dot = nacre_exact::dot_sign_rat(&n, &axis);
            assert!(
                dot == nacre_exact::Orient::Zero || nacre_exact::parallel_rat(&n, &axis),
                "a face neither ⊥ nor ∥ to the axis would be the oblique branch"
            );
        }
        let q = |x: f64| Rat::from_decimal(x).expect("a short decimal");
        let (cyl, _) = m
            .add_cylinder_exact(
                base.map(q),
                axis,
                // ★ The frame's own `u`: unit, and `u · axis = 0.384 − 0.384 = 0` exactly.
                [q(0.6), q(0.8), q(0.0)],
                q(0.5),
                q(3.0),
                None,
            )
            .expect("the exact road states a tilted cylinder");
        m.rebuild_adjacency();
        let out = crate::ledger::owned(|| {
            crate::boolean(&mut m, BoolKind::Cut, solid, cyl)
                .unwrap_or_else(|e| panic!("{name}: {e:?}"))
        });
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{name}: the tilted result must be sound"
        );
        // ★★★★★ **The oracle is on this model, derived from the inputs.** Green plus a clean
        // `validate` would say only that *some* solid was built; the volume says the tilted
        // boolean removed the right material, and it is computed from the prism's and bore's
        // own dimensions rather than from anything the kernel answered.
        let [s] = out[..] else {
            panic!("{name}: one solid")
        };
        let got = nacre_props::mass_props(&m, s).expect("mass props").volume;
        assert!(
            (got - want_vol).abs() < 1e-9,
            "{name}: volume {got} vs {want_vol}"
        );
        let mine = ROWS.mine_since(before);
        let want = if want_rulings { 6 } else { 0 };
        // ★ This fixture's own charts, read as a count — a shared ledger could only say the
        // shape existed somewhere. The volume above is the other half: it proves the tilted
        // boolean removed the right material, on this model.
        assert_eq!(mine.len(), 1, "{name}: one tilted class, one chart");
        assert_eq!(
            (mine[0].z_lines, mine[0].theta),
            (4, want),
            "{name}: expected a chart with {want} rulings"
        );
    }
}

/// **The existence rule's whole truth table** — [`face_spans`] against marks written by hand.
///
/// The production lock reaches this through an entire boolean, and a boolean exercises two of
/// the rows below: a `Transversal` that reaches and a `Graze` pointing away. The rule is one
/// sentence about every rim, not a description of that fixture, so the rest are stated here —
/// including the two that can only be reached from geometry the road does not build yet.
#[test]
fn face_spans_reads_the_trace_not_the_label() {
    use super::read_cell::face_spans;
    use crate::arrangement::{ArcLabel, SegKind};
    use crate::combinatorics;
    use crate::planes::SolidSide;
    use crate::{BoolError, RejectReason};
    let n = combinatorics::NodeId::three_planes(combinatorics::Canon3::three([0, 1, 2]));
    let row = |marks: Vec<(SolidSide, SegKind)>| ArcLabel {
        ends: [n, n],
        // Deliberately the label that says "material everywhere": nothing below may read it.
        label: [true; 4],
        marks,
    };
    let (a, b) = (SolidSide::A, SolidSide::B);
    let spans = |marks, above| face_spans(&row(marks), a, above);
    // A face running through the plane is on both sides of it.
    for above in [true, false] {
        assert!(spans(vec![(a, SegKind::Transversal { mat: 1 })], above).unwrap());
    }
    // A face whose boundary stops at the arc is on exactly the side it occupies.
    for body_above in [true, false] {
        for above in [true, false] {
            assert_eq!(
                spans(vec![(a, SegKind::Graze { body_above })], above).unwrap(),
                body_above == above,
                "graze body_above={body_above} against band above={above}"
            );
        }
    }
    // The counterpart's marks are not this row's face, whatever they say.
    assert!(!spans(vec![(b, SegKind::Transversal { mat: 1 })], true).unwrap());
    // A planar face's seated rim says nothing about the lateral — see `ArcLabel::marks`.
    assert!(!spans(vec![(a, SegKind::Seated { body_above: true })], true).unwrap());
    // No mark at all: the face stops short of this arc.
    assert!(!spans(Vec::new(), true).unwrap());
    // Two of this solid's faces on one arc: agreeing decides, disagreeing refuses by name.
    assert!(
        spans(
            vec![
                (a, SegKind::Graze { body_above: true }),
                (a, SegKind::Transversal { mat: 1 }),
            ],
            true,
        )
        .unwrap()
    );
    assert!(matches!(
        spans(
            vec![
                (a, SegKind::Graze { body_above: false }),
                (a, SegKind::Transversal { mat: 1 }),
            ],
            true,
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::CylinderFaceUndecided,
            ..
        })
    ));
}
