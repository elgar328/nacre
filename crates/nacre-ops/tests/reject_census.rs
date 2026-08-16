//! **The standing reject census: a frozen corpus of rejecting shapes, and what each one rings.**
//!
//! `nacre_ops::reject_census` records every guard that rings and every rejection that leaves the
//! kernel. This file is the population it records over — chosen shapes, run one at a time, so the
//! table says *which shape* produced *which reason at which file*.
//!
//! ★ **One `#[test]`, deliberately.** `reject_census::take()` resets a process-global, so a second
//! test in this binary would steal its counts. Binaries are separate processes, so the rest of the
//! suite is unaffected.
//!
//! ★★ **The fixtures are copied from the tests they came from, on purpose.** This population has
//! to be frozen: if `self_touch.rs` changes its wedge, the census must not move with it. It is the
//! same reason `tests/census.rs` hard-codes its coordinates instead of sharing helpers. What keeps
//! a copy from drifting silently is that **every fixture asserts its own reject reason** — a
//! drifted fixture then fails loudly instead of quietly re-baselining the table.
//!
//! ★★★ **What this gate does *not* see.** Its population is these seven shapes. A guard added for
//! M6 quadrics will not be exercised by planar fixtures and this stays green — the whole-suite
//! sweep (`--features reject-trace`, see the module doc) is what covers that, by hand. It also
//! cannot see a site moving *within* a file, and it does not contain `precision_budget` (a
//! 4200-turn chain) or `edge_occupancy_conflict` (measured: reachable only by calling `edge_mask`
//! directly, never through `boolean()`).

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::reject_census::{self, Census, ReasonId};
use nacre_ops::{
    BoolError, BoolKind, DatumDef, OpOutput, Operation, Profile2d, RejectReason, SketchFrame,
    SketchPlane, apply, boolean,
};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};
use std::collections::BTreeSet;

// ── shape helpers, copied from the tests each fixture came from ───────────────────────────────

fn cub(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
    let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
    m.rebuild_adjacency();
    s
}

/// A prism raised from `z` on the world XY plane (`contact_separates.rs`'s helper).
fn prism_z(m: &mut Model, pts: &[[f64; 2]], z: f64, dist: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("poly");
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(
                SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
            ),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    };
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
}

/// A prism raised along `axis` (`rotation_sweep.rs`'s helper).
fn prism_axis(m: &mut Model, pts: &[[f64; 2]], axis: Axis, dist: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, axis),
            profile,
            dist,
        },
    )
    .expect("extrude") else {
        panic!("extrude output")
    };
    m.rebuild_adjacency();
    solid
}

fn moved(m: &mut Model, s: Handle<Solid>, isometry: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } =
        apply(m, &Operation::Transform { solid: s, isometry }).expect("transform")
    else {
        panic!("transform output")
    };
    m.rebuild_adjacency();
    solid
}

fn copy_of(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let OpOutput::Copy { solid } = apply(m, &Operation::Copy { solid: s }).expect("copy") else {
        panic!("copy output")
    };
    m.rebuild_adjacency();
    solid
}

fn rot_z(m: &mut Model, s: Handle<Solid>, deg: i128) -> Handle<Solid> {
    moved(
        m,
        s,
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
        }),
    )
}

// ── the corpus ────────────────────────────────────────────────────────────────────────────────

/// A shape, what it must reject with, and the census it must produce.
struct Fixture {
    name: &'static str,
    /// `None` = it must succeed (the control).
    expect: Option<RejectReason>,
    run: fn(&mut Model) -> Result<Vec<Handle<Solid>>, BoolError>,
    /// **The pinned raise shape**: `(reason, detail, file)`. Read off a measurement, not predicted
    /// — and measured identical in all four of debug/release x parallel/serial, which is why it is
    /// asserted unconditionally.
    ///
    /// No line numbers: a gate keyed on them goes red on every unrelated edit above a guard, which
    /// is how a lock stops being read. Lines are printed instead.
    raised: &'static [(&'static str, Option<&'static str>, &'static str)],
    /// **The pinned surfaced shape** — the strong half. `par::try_map_range` returns the
    /// lowest-index error by construction, so which reason comes back does not depend on the
    /// schedule or the build.
    surfaced: &'static [(&'static str, Option<&'static str>)],
}

/// ① A fuse that pinches one solid at a single vertex — two bridges run around the contact, so the
/// material loops and there is no pair of solids to hand back (`coverage/rejects.rs`).
fn pinched_vertex(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let g1 = cub(m, [1.0, 0.3, 0.2], [4.5, 1.3, 0.8]);
    let g2 = cub(m, [3.5, 0.3, 0.2], [4.5, 3.8, 1.6]);
    let b = cub(m, [2.0, 2.0, 1.0], [4.0, 4.0, 2.0]);
    let t1 = boolean(m, BoolKind::Fuse, a, g1).expect("a and g1 overlap");
    m.rebuild_adjacency();
    let t2 = boolean(m, BoolKind::Fuse, t1[0], g2).expect("g1 and g2 overlap");
    m.rebuild_adjacency();
    boolean(m, BoolKind::Fuse, t2[0], b)
}

/// ② The edge twin: A and B meet only along `x = 2, y = 2`, and a bridge takes the material around
/// the contact, so four faces use that segment (`coverage/rejects.rs`).
fn pinched_edge(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let bridge = cub(m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
    let b = cub(m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
    let ab = boolean(m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
    m.rebuild_adjacency();
    boolean(m, BoolKind::Fuse, ab[0], b)
}

/// ③ A wedge sunk into a cube whose sharp tip lands exactly on the far wall: the cut's cavity
/// meets the outer shell along a line, so the result's surface touches itself (`self_touch.rs`).
fn wedge_tip_on_wall(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [1.0, 1.0, 1.0]);
    let b = prism_axis(m, &[[1.0, 0.5], [0.1, 0.1], [0.1, 0.9]], Axis::Z, 0.6);
    let b = moved(
        m,
        b,
        Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::new(2, 10).unwrap()]),
    );
    boolean(m, BoolKind::Cut, a, b)
}

/// ④ ★ **The corpus's carrier of swallowed raises.** A diamond void whose four corners all sit on
/// the block's walls: every candidate node for the nesting ray grazes, so the retry rings the
/// `no_clear_ray` guard once per node before the caller runs out of nodes and reports it
/// (`contact_separates.rs`). Without this shape the raised/surfaced split would measure nothing.
fn diamond_void(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let b = prism_z(
        m,
        &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
        0.0,
        3.0,
    );
    let diamond = prism_z(
        m,
        &[[0.0, 2.0], [2.0, 0.0], [4.0, 2.0], [2.0, 4.0]],
        1.0,
        1.0,
    );
    boolean(m, BoolKind::Cut, b, diamond)
}

/// ⑤ A non-planar operand — the coverage limit named at the plane table, not in the assembly
/// (`lib.rs::common_rejects_non_planar_input`).
fn cylinder_operand(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [2.0; 3]);
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Common, a, cyl)
}

/// ⑥ ★ **The target.** A part fused with a 45°-rotated copy of itself: one arm lands coplanar on
/// another, giving a solid of zero thickness. What comes back is `CoplanarPinch` — the merge
/// guard's own capability limit, under its own name since 2026-08-16 (this census measured that
/// it was the only firing site of the old shared `CoplanarMerge` label). The deeper truth (a
/// self-touch, `Impossible`) still needs the figure-8 merge: the measured detour of skipping the
/// merge died one check later as `StraightAngle` (`rotation_sweep.rs`, dev-log).
fn fold_45(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let plate = [
        [0.0, 0.0],
        [50.0, 0.0],
        [50.0, 25.0],
        [38.0, 25.0],
        [38.0, 50.0],
        [50.0, 50.0],
        [50.0, 75.0],
        [0.0, 75.0],
    ];
    let p1 = prism_axis(m, &plate, Axis::Z, 12.0);
    let bar = [[20.0, 12.0], [75.0, 12.0], [75.0, 37.0], [55.0, 37.0]];
    let p2 = prism_axis(m, &bar, Axis::X, 25.0);
    let unit = boolean(m, BoolKind::Fuse, p1, p2).expect("the part fuses")[0];
    m.rebuild_adjacency();

    let mut part = copy_of(m, unit);
    for deg in (45..360).step_by(45) {
        let c = copy_of(m, unit);
        let c = rot_z(m, c, deg as i128);
        match boolean(m, BoolKind::Fuse, part, c) {
            Ok(out) => {
                part = out[0];
                m.rebuild_adjacency();
            }
            Err(e) => return Err(e),
        }
    }
    Ok(vec![part])
}

/// ⑦ ★ **The control.** Two overlapping cubes fuse. A boolean that succeeds must leave the census
/// empty — otherwise every row above could be noise from the machinery rather than from the shape.
fn plain_fuse(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [1.0; 3]);
    let b = cub(m, [0.5; 3], [1.5; 3]);
    boolean(m, BoolKind::Fuse, a, b)
}

const BOOLEAN: &str = "crates/nacre-ops/src/boolean.rs";

const CORPUS: [Fixture; 7] = [
    Fixture {
        name: "pinched-vertex",
        expect: Some(RejectReason::NonManifoldVertex),
        run: pinched_vertex,
        raised: &[("non_manifold_vertex", None, BOOLEAN)],
        surfaced: &[("non_manifold_vertex", None)],
    },
    Fixture {
        name: "pinched-edge",
        expect: Some(RejectReason::NonManifoldResultEdge),
        run: pinched_edge,
        raised: &[("non_manifold_result_edge", None, BOOLEAN)],
        surfaced: &[("non_manifold_result_edge", None)],
    },
    Fixture {
        name: "wedge-tip-on-wall",
        expect: Some(RejectReason::SelfTouchingResult),
        run: wedge_tip_on_wall,
        raised: &[("self_touching_result", None, BOOLEAN)],
        surfaced: &[("self_touching_result", None)],
    },
    Fixture {
        name: "diamond-void",
        expect: Some(RejectReason::NoClearRay),
        run: diamond_void,
        // ★ One site now, and it is the one the caller is told about. This fixture used to be
        // the census's carrier of swallowed raises — its per-node probe rang `no_clear_ray` from
        // `combinatorics.rs` 24 times per boolean and surfaced once — until E (2026-08-17) made
        // the probe's "this node cannot decide" an abstention in its type (`Ok(None)`) instead
        // of an error for the retry to catch. What rings now is what surfaces.
        raised: &[("no_clear_ray", None, BOOLEAN)],
        surfaced: &[("no_clear_ray", None)],
    },
    Fixture {
        name: "cylinder-operand",
        expect: Some(RejectReason::CylinderFace),
        run: cylinder_operand,
        // The only corpus reject raised while reading the operands rather than in the assembly.
        raised: &[("cylinder_face", None, "crates/nacre-ops/src/planes.rs")],
        surfaced: &[("cylinder_face", None)],
    },
    Fixture {
        name: "fold-45",
        expect: Some(RejectReason::CoplanarPinch),
        run: fold_45,
        raised: &[("coplanar_pinch", None, BOOLEAN)],
        surfaced: &[("coplanar_pinch", None)],
    },
    Fixture {
        name: "plain-fuse",
        expect: None,
        run: plain_fuse,
        raised: &[],
        surfaced: &[],
    },
];

// ── the run ───────────────────────────────────────────────────────────────────────────────────

#[test]
fn the_reject_census() {
    // Anything the harness did before this point (nothing, in a single-test binary) is discarded,
    // so the first fixture's window is its own.
    let _ = reject_census::take();

    instrument_answers_for_itself();

    let mut table: Vec<(&'static str, Census)> = Vec::new();
    for f in &CORPUS {
        let mut m = Model::new();
        let got = (f.run)(&mut m);
        let census = reject_census::take();

        // ★ Assertion 1: the fixture states its own reason. A copied shape that drifts fails here
        // rather than quietly re-baselining the table below.
        match f.expect {
            Some(reason) => assert_eq!(
                got.err().map(|e| match e {
                    BoolError::Rejected { reason, .. } => reason,
                    other => panic!("{}: expected a named rejection, got {other:?}", f.name),
                }),
                Some(reason),
                "{}: the fixture stopped producing the reject it exists for",
                f.name
            ),
            None => assert!(
                got.is_ok(),
                "{}: the control stopped building — {got:?}",
                f.name
            ),
        }

        // ★ Assertion 3: which reason is raised in which file. This is what nothing else in the
        // suite pins — every reject test asserts the *reason*, none of them where it lives.
        let want: BTreeSet<(ReasonId, &str)> = f
            .raised
            .iter()
            .map(|&(reason, detail, file)| (ReasonId { reason, detail }, file))
            .collect();
        assert_eq!(
            census.raised_shape(),
            want,
            "{}: the guards that rang changed.\n{}\
             If that is intended, update this fixture's `raised`; if not, ask why a guard moved.",
            f.name,
            census.report()
        );

        // ★ Assertion 2: what the caller was told. The strong half — identical in every build.
        let want: BTreeSet<ReasonId> = f
            .surfaced
            .iter()
            .map(|&(reason, detail)| ReasonId { reason, detail })
            .collect();
        assert_eq!(
            census.surfaced_shape(),
            want,
            "{}: the reasons leaving the kernel changed.\n{}",
            f.name,
            census.report()
        );

        table.push((f.name, census));
    }

    #[allow(clippy::print_stdout)]
    {
        for (name, c) in &table {
            println!("── {name}\n{}", c.report());
        }
    }

    let by = |name: &str| -> &Census {
        &table
            .iter()
            .find(|(n, _)| *n == name)
            .expect("fixture present")
            .1
    };

    // ★★ Assertion 4: **the point-contact merge guard is one site, under its own name.** It was
    // one of `CoplanarMerge`'s ten sites — the only one that ever fired — and this census is what
    // measured that and licensed giving it its own name (`CoplanarPinch`, 2026-08-16). The
    // baseline going red on that change was the census working as designed: the vocabulary change
    // showed up in the diff. Distinct *sites*, not raise counts: counts move with the build
    // (`debug` re-runs every traced boolean, `parallel` evaluates the classes a failing input
    // would have skipped), and the claim is about how many guards speak, not how often.
    assert_eq!(
        by("fold-45").distinct_sites("coplanar_pinch"),
        1,
        "a second `coplanar_pinch` guard started ringing.\n{}",
        by("fold-45").report()
    );
    // ★ The positive control for that counter: it can answer something other than 1. The
    // two-site carrier this used to ride (diamond-void's swallowed probe raises) was removed by
    // E on purpose, so the control is the other direction now — a fixture where the reason
    // never rings at all.
    assert_eq!(
        by("plain-fuse").distinct_sites("coplanar_pinch"),
        0,
        "a successful boolean rang the pinch guard.\n{}",
        by("plain-fuse").report()
    );

    // ★ What used to be Assertion 5 — "at least one reason is raised more often than it
    // surfaces" — is gone **because the population it rode on was deliberately removed**: E
    // (2026-08-17) turned the swallowed-raise champion (`point_in_component`'s per-node probe,
    // 122 of the suite's 151 raises) into a typed abstention. The two columns still measure
    // something, but the reading inverted: their gap is now the count of *remaining* swallowed
    // raises, and near-zero is the goal state — a gap reopening here is news of a new
    // swallow, visible in the whole-suite sweep (`--features reject-trace`).

    // ★ Assertion 6: a successful boolean rings no guard at all.
    assert!(
        by("plain-fuse").is_empty(),
        "a boolean that succeeded still rang a guard: {}",
        by("plain-fuse").report()
    );
}

/// The census before the census: does the instrument point where it claims to?
fn instrument_answers_for_itself() {
    // `#[track_caller]` on `reject()` must report the *guard's* line, not `reject`'s own. Nothing
    // else in this file would notice if that attribute were dropped, and every row of the table
    // would silently collapse onto one line in `lib.rs`.
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    m.rebuild_adjacency();
    let _ = boolean(&mut m, BoolKind::Common, a, cyl);
    let c = reject_census::take();
    let site = c
        .raised
        .iter()
        .find(|(s, _)| s.id.reason == "cylinder_face")
        .map(|(s, _)| *s)
        .expect("the cylinder guard rang");
    assert!(
        site.file.ends_with("planes.rs"),
        "the raise site is reported as {}:{} — `#[track_caller]` is naming `reject()` itself, \
         not the guard",
        site.file,
        site.line
    );

    // `take()` resets: the window just read must not appear in the next one.
    assert!(
        reject_census::take().is_empty(),
        "take() did not reset — every fixture's window would include its predecessors'"
    );
}
