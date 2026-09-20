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
            pivot: [Rat::from_int(0); 3],
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

fn cylinder_notyet(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [2.0; 3]);
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Cut, a, cyl)
}

/// **The same drill, at a scale where the exact arithmetic used to run out.** Geometry decides
/// nothing new here — an axis-aligned box with an axis-aligned bore through it, walls clear of
/// the lateral surface — but every coordinate is a sub-micron value carried to a full f64's
/// digits, so its exact rational has a ~10²³ denominator and the gate's `(n·o+d)² vs r²|n|²`
/// squares it past `i128`.
///
/// ★ **Chosen from the algebra, not the story** (measured): a *tilted* wide axis cannot test
/// this — a real tilt makes the box's planes oblique and the gate rejects for that first — and
/// a unit-scale box with a wide radius does not overflow at all (`0.5000000000000001²` fits).
/// Smallness with full digits is the property that bites. Before the gate's predicates became
/// total this fixture answered `CylinderGateUndecided`: the arithmetic ran out, and the reject
/// named the symptom instead of the geometry.
fn cylinder_wide_axis(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(
        m,
        [0.0; 3],
        [
            2.0000000000000003e-7,
            2.0000000000000003e-7,
            2.0000000000000003e-7,
        ],
    );
    let cyl = m.add_cylinder(
        Point3::from_array([
            1.0000000000000002e-7,
            1.0000000000000002e-7,
            -1.0000000000000002e-7,
        ]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        5.000000000000001e-8,
        4.000000000000001e-7,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Cut, a, cyl)
}

fn cylinder_wall_contact(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [2.0; 3]);
    let cyl = m.add_cylinder(
        Point3::from_array([0.3, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Cut, a, cyl)
}

/// The tangent wall: the axis sits exactly `r = 0.5` from the `x = 0` wall, so the wall touches
/// the lateral along one line. The gate passes it; this `Cut` is refused by the
/// **verdict**, because what it leaves near the line is two wedges of one solid.
fn cylinder_wall_tangent(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = cub(m, [0.0; 3], [2.0; 3]);
    let cyl = m.add_cylinder(
        Point3::from_array([0.5, 1.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Cut, a, cyl)
}

fn cylinder_oblique(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    // A tilted rational axis against an axis-aligned box: the box's planes are neither ⊥ nor
    // ∥ to (0,1,1), and — unlike rotating the box — every description stays world-rational,
    // so the gate reaches the *oblique* verdict rather than declining as undecidable.
    let a = cub(m, [0.0; 3], [2.0; 3]);
    let cyl = m.add_cylinder(
        Point3::from_array([1.0, 1.0, -2.0]),
        Vector3::from_array([0.0, 1.0, 1.0]),
        0.5,
        6.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Cut, a, cyl)
}

fn cylinder_pair(m: &mut Model) -> Result<Vec<Handle<Solid>>, BoolError> {
    let a = m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        2.0,
    );
    let b = m.add_cylinder(
        Point3::from_array([0.6, 0.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        3.0,
    );
    m.rebuild_adjacency();
    boolean(m, BoolKind::Fuse, a, b)
}

/// ⑥ ★ **The target.** A part fused with a 45°-rotated copy of itself: one arm lands coplanar on
/// another, giving a solid of zero thickness — a self-touch. What comes back is its truth,
/// `SelfTouchingResult` (`Impossible`, the touching edge as the witness): the
/// merge abstains on a pinching group and the whole-result judgement runs before minting.
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

// Where the assembly's guards ring from, one file per stage.
const ENTRY: &str = "crates/nacre-ops/src/boolean/entry.rs";
const RECONSTRUCT: &str = "crates/nacre-ops/src/boolean/reconstruct.rs";
const SELF_TOUCH: &str = "crates/nacre-ops/src/boolean/self_touch.rs";

const CORPUS: [Fixture; 13] = [
    Fixture {
        name: "pinched-vertex",
        expect: Some(RejectReason::NonManifoldVertex),
        run: pinched_vertex,
        raised: &[("non_manifold_vertex", None, ENTRY)],
        surfaced: &[("non_manifold_vertex", None)],
    },
    Fixture {
        name: "pinched-edge",
        expect: Some(RejectReason::NonManifoldResultEdge),
        run: pinched_edge,
        raised: &[("non_manifold_result_edge", None, RECONSTRUCT)],
        surfaced: &[("non_manifold_result_edge", None)],
    },
    Fixture {
        name: "wedge-tip-on-wall",
        expect: Some(RejectReason::SelfTouchingResult),
        run: wedge_tip_on_wall,
        raised: &[("self_touching_result", None, SELF_TOUCH)],
        surfaced: &[("self_touching_result", None)],
    },
    Fixture {
        name: "diamond-void",
        // ★ One site, and it is the one the caller is told about. The per-node probe's "this
        // node cannot decide" is an abstention in its type (`Ok(None)`), not an error for the
        // retry to catch (as an error it rang `no_clear_ray` 24 times per boolean and surfaced
        // once). What rings is what surfaces.
        //
        // ★★ **And what surfaces is the shape's own name.** `no_clear_ray` was the
        // depth question running out of witnesses — a planar component was offered nothing but its
        // vertices, and here every one of them sits on the host's wall. Given the points its
        // **edges** name, the depth decides and the self-touch test says what this really is: a void
        // whose corner rides a wall makes the result's surface meet itself along that corner's
        // edge, which is exactly what the one-grazing-corner sibling in `contact_separates.rs` has
        // always reported. ⚠ This was `NoClearRay`'s only *surfacing* fixture; the name still
        // raises on the 2-D ring road and the `reject-trace` sweep is what keeps it visible.
        expect: Some(RejectReason::SelfTouchingResult),
        run: diamond_void,
        raised: &[("self_touching_result", None, SELF_TOUCH)],
        surfaced: &[("self_touching_result", None)],
    },
    // ── The cylinder population gate names its refusals — one fixture per cause. All are
    // raised while reading the operands (`plane_index_setup`), before any arrangement work.
    Fixture {
        name: "cylinder-operand",
        // ★ **The seated-cap row, and the second diff of that kind.** Flush caps: the cylinder's
        // z ∈ [0,2] caps intern onto the box's own cap planes, so every ⊥ class carries faces of
        // both operands. That used to be `SeatedCylinderCap` — a rule wider than anything it
        // could name, since what makes a seated circle hard is its *boundary*, and every way a
        // boundary can be met is already refused by the wall, oblique or pair rule. It **builds**
        // now, and the name it carried no longer exists in the code.
        expect: None,
        run: cylinder_operand,
        raised: &[],
        surfaced: &[],
    },
    Fixture {
        name: "cylinder-notyet",
        // ★ The drill population passes the gate and **builds**: a through hole, walls clear
        // (distance 1 > r = 0.5), caps unshared.
        expect: None,
        run: cylinder_notyet,
        raised: &[],
        surfaced: &[],
    },
    Fixture {
        name: "cylinder-wide-axis",
        // ★ The same population spelled with long decimals. The **gate** answers it (its
        // questions are total), and the wall this input met next was the **value** path: the circle
        // nesting projects the ring's corners into the class's rational chart, and a sub-micron
        // model at full f64 precision left `Rat` there.
        //
        // ★★ **It builds, and the decline it once met was not about this model at
        // all.** The chart's second axis is `n x e1`, which squares the normal — and that normal
        // was carrying a factor the plane-name canonicalisation had not removed, because that
        // canonicalisation divides the content out of **four** coefficients while the chart reads
        // three. On `z = s` with a long decimal `s`, the leftover factor is `s`'s denominator. It
        // is divided out now (`combinatorics::primitive_normal`), the axes come out along the same
        // directions as ever, and this cut is exact: one solid, `validate` clean, volume matching
        // `s^3 - pi r^2 s`. The fixture stays as coverage — a model three orders of magnitude
        // below anything else here.
        //
        // ⚠ It was `WitnessNotRational`'s only *surfacing* fixture. That name still raises inside
        // the engine and the `reject-trace` sweep is what keeps it visible; nothing was invented
        // to replace this row.
        expect: None,
        run: cylinder_wide_axis,
        raised: &[],
        surfaced: &[],
    },
    Fixture {
        name: "cylinder-wall-contact",
        // The axis sits 0.3 from the x = 0 wall with r = 0.5 — the wall pierces the lateral
        // surface in two rulings (`y = 1 ± 0.4`). ★ It **builds** (the offset wall):
        // the gate records the pair and the tracer cuts the lateral along the rulings and the
        // caps along the chord. The name it carried stays in the code for the tangent alone.
        expect: None,
        run: cylinder_wall_contact,
        raised: &[],
        surfaced: &[],
    },
    Fixture {
        name: "cylinder-wall-tangent",
        // The axis exactly r from the wall: a zero-thickness contact along one line, which
        // `validate` cannot see, which is why the answer is a verdict rather
        // than a fence. The reason comes from the
        // assembly, where the grouping can say the two wedges are one body.
        expect: Some(RejectReason::SelfTouchingResult),
        run: cylinder_wall_tangent,
        raised: &[("self_touching_result", None, SELF_TOUCH)],
        surfaced: &[("self_touching_result", None)],
    },
    Fixture {
        name: "cylinder-oblique",
        // A 30°-turned box: its planes are neither ⊥ nor ∥ to the axis — ellipses.
        expect: Some(RejectReason::ObliqueCylinderCut),
        run: cylinder_oblique,
        raised: &[(
            "oblique_cylinder_cut",
            None,
            "crates/nacre-ops/src/planes/cyl_gate.rs",
        )],
        surfaced: &[("oblique_cylinder_cut", None)],
    },
    Fixture {
        name: "cylinder-pair",
        // Two overlapping parallel cylinders: axis distance 0.6 < r₁+r₂ = 1 (M6b).
        expect: Some(RejectReason::CylinderPairContact),
        run: cylinder_pair,
        raised: &[(
            "cylinder_pair_contact",
            None,
            "crates/nacre-ops/src/planes/cyl_gate.rs",
        )],
        surfaced: &[("cylinder_pair_contact", None)],
    },
    Fixture {
        name: "fold-45",
        expect: Some(RejectReason::SelfTouchingResult),
        run: fold_45,
        raised: &[("self_touching_result", None, SELF_TOUCH)],
        surfaced: &[("self_touching_result", None)],
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

    // ★★ Assertion 4: **the fold's truth is one site, the whole-result judgement.** The
    // point-contact merge guard abstains — the merge emits a pinching group unmerged and
    // `self_touch_reject`, hoisted before minting, names the shape. Distinct
    // *sites*, not raise counts: counts move with the build (`debug` re-runs every traced
    // boolean, `parallel` evaluates the classes a failing input would have skipped), and the
    // claim is about how many guards speak, not how often.
    assert_eq!(
        by("fold-45").distinct_sites("self_touching_result"),
        1,
        "a second `self_touching_result` guard started ringing.\n{}",
        by("fold-45").report()
    );
    // ★ The positive control for that counter: it can answer something other than 1. The
    // two-site carrier this used to ride (diamond-void's swallowed probe raises) was removed by
    // E on purpose, so the control is the other direction now — a fixture where the reason
    // never rings at all.
    assert_eq!(
        by("plain-fuse").distinct_sites("self_touching_result"),
        0,
        "a successful boolean rang the self-touch guard.\n{}",
        by("plain-fuse").report()
    );

    // ★ There is no Assertion 5 ("at least one reason is raised more often than it
    // surfaces"): the swallowed-raise champion (`point_in_component`'s per-node probe,
    // 122 of the suite's 151 raises) is a typed abstention. The two columns still measure
    // something: their gap is the count of *remaining* swallowed
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
    // would silently collapse onto one line in `error.rs`.
    // ★ The oblique cut, not the tangent wall: the tangency's reason now comes from the
    // *assembly* (`boolean/self_touch.rs`), and what this checks is that a **gate** guard names its own line.
    let mut m = Model::new();
    let _ = cylinder_oblique(&mut m);
    let c = reject_census::take();
    let site = c
        .raised
        .iter()
        .find(|(s, _)| s.id.reason == "oblique_cylinder_cut")
        .map(|(s, _)| *s)
        .expect("the cylinder guard rang");
    assert!(
        site.file.ends_with("planes/cyl_gate.rs"),
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
