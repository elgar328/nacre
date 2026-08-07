//! **How wide a discovered coordinate actually is — measured.**
//!
//! `docs/truth-and-cache.md`'s rule 1 says the truth is a rational *or* a handle, and that values
//! which do not fit `Rat` — *"발견된 좌표 160–480 비트"* — are pointed at rather than stored. That
//! number decides whether the last migration row needs a new plane variant or not, and it had
//! never been measured in this kernel. This is that measurement.
//!
//! ★ **What the instrument had to be.** `three_planes_rat` cannot answer the question it is
//! about: it is `checked_*` throughout, so the widest thing it can report is ~127 bits and
//! everything it reports fits by construction. `three_planes_big` (scalar) removes the ceiling,
//! which also splits its predecessor's `None` into the two different facts it was conflating —
//! *the point does not fit* and *an intermediate overflowed on a point that would have*.
//!
//! ★★ **What is not measured here, and why.** Not the *world* coordinate. `Pt3` is
//! `{ base: [Rat; 3], chain, coord: [f64; 3], … }` — a moved point's world position is never
//! stored as a rational at all; the kernel keeps a narrow base and **points at** the motion,
//! realizing at whatever precision a judgment needs. So there is no world-frame width to have.
//! The width question is entirely about `base`, which is what a plane name pins.
//!
//! Populations come from the kernel's own discriminator rather than from a story: `vertex_tol` is
//! `Some` for a **discovered** vertex and `None` for a **constructed** one, which is exactly the
//! set rule 1 is about. Constructed vertices ride along as the free control.
//!
//! `OnSeam` vertices are excluded rather than counted as failures: that variant pins a curve, not
//! a point, and its doc already records the coordinate cache as load-bearing there.
//!
//! # What it measured (2026-08-07)
//!
//! | population | vertices | solved | **width max** | over 127 | narrow declined, point fits |
//! |---|---|---|---|---|---|
//! | `boolean_corner`   | 16 | 16 |  **7** | 0 | 0 |
//! | `boolean_twice`    | 32 | 32 |  **7** | 0 | 0 |
//! | `boolean_rotated`  | 20 |  8 |  **4** | 0 | 0 |
//! | `tilted_frame`     | 16 | 12 | **59** | 0 | **8 of 12** |
//! | `tilted_frame_x2`  | 24 | 16 | **59** | 0 | **8 of 16** |
//!
//! ★★★ **Two things came out, and neither was predicted.**
//!
//! 1. **Nothing is wide.** The prediction written before the run was 150–170 bits growing with
//!    depth, from the neighbouring 2026-08-04 measurement (coefficients ~50 bits, Cramer
//!    multiplying three of them). The corpus maximum is **59 bits**, and a second feature on the
//!    already-awkward tilted population moved it **not at all** (59 → 59). Depth does not
//!    compound here: each feature's planes are stated afresh in decimals, so the solve never
//!    stacks. ★ These are **corpus numbers, not bounds** — the meter's negative control lives in
//!    `nacre-scalar` (`the_width_meter_reports_a_point_no_rat_can_hold`, 201 bits), so `over127 =
//!    0` is a fact about this population and not about a clamped instrument.
//!
//! 2. ★★ **The declines are not about width at all.** On the tilted decimal family
//!    `three_planes_rat` gives up on two thirds of the vertices it is offered, and for **every
//!    one of them the point fits `Rat`** (`declined_and_wide = 0`). Its `None` there means its
//!    own rational cofactor expansion overflowed, not that the coordinate is unstorable.
//!
//! ★ **The corpus is small** — five models, 108 vertices — where the adjacent measurement it is
//! being compared against had 83,821 samples. It is enough to refute "everything is 160–480 bits"
//! (one counterexample population does that) and enough to establish the decline split. It is not
//! enough to say what the widest coordinate this kernel can produce is.

use nacre_geom::Surface;
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{BoolKind, DatumDef, OpOutput, Operation, Profile2d, SketchFrame, SketchPlane};
use nacre_ops::{apply, boolean};
use nacre_scalar::{Angle, Axis, Isometry, MeetPoint, PlaneName, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, SurfaceTruth, Vertex, VertexDef};

// ---------------------------------------------------------------------------------------------
// fixtures — the shapes `tests/points_coverage.rs` already uses
// ---------------------------------------------------------------------------------------------

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

fn cuboid(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
    m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
}

// ---------------------------------------------------------------------------------------------
// the tally
// ---------------------------------------------------------------------------------------------

#[derive(Default, Debug)]
struct Tally {
    /// Every vertex walked.
    vertices: usize,
    /// `vertex_tol` is `None` — a constructed vertex (the control).
    constructed: usize,
    /// `vertex_tol` is `Some(0.0)` — discovered, residual measured exactly zero.
    discovered_exact: usize,
    /// `vertex_tol` is `Some(t)`, `t > 0` — discovered, genuinely inexact.
    discovered_inexact: usize,
    /// `OnSeam` — excluded from the width question, not a failure.
    on_seam: usize,
    /// A carrier had no recorded name at all.
    unnamed_carrier: usize,
    /// ★ A carrier's name is `Wide` — invisible to `three_planes_rat`, solvable by the twin.
    wide_carrier: usize,
    /// The three carriers do not share one motion — a **structural** limit of the exact road,
    /// nothing to do with width. Counted apart so it cannot be read as "does not fit".
    mixed_motion: usize,
    /// `three_planes_big` declined: the determinant is zero, no unique point.
    singular: usize,
    /// Solved. `width_bits` of each, in the order walked.
    widths: Vec<u64>,
    /// Solved, and split by the tol class the vertex came from.
    widths_exact: Vec<u64>,
    widths_inexact: Vec<u64>,
    /// ★★ The two causes of `three_planes_rat`'s `None`, told apart by the twin's answer.
    narrow_answered: usize,
    declined_but_fits: usize,
    declined_and_wide: usize,
}

impl Tally {
    fn max(v: &[u64]) -> u64 {
        v.iter().copied().max().unwrap_or(0)
    }

    fn over(v: &[u64], bits: u64) -> usize {
        v.iter().filter(|w| **w > bits).count()
    }

    fn report(&self, what: &str) {
        println!(
            "stat {what:22} vertices={} constructed={} discovered(exact0={} inexact={}) seam={}",
            self.vertices,
            self.constructed,
            self.discovered_exact,
            self.discovered_inexact,
            self.on_seam
        );
        println!(
            "stat {what:22} solved={} singular={} mixed_motion={} unnamed={} wide_carrier={}",
            self.widths.len(),
            self.singular,
            self.mixed_motion,
            self.unnamed_carrier,
            self.wide_carrier
        );
        println!(
            "stat {what:22} width max={} over127={} | tol==0 max={} | tol>0 max={}",
            Self::max(&self.widths),
            Self::over(&self.widths, 127),
            Self::max(&self.widths_exact),
            Self::max(&self.widths_inexact),
        );
        println!(
            "stat {what:22} narrow_answered={} declined_but_fits={} declined_and_wide={}",
            self.narrow_answered, self.declined_but_fits, self.declined_and_wide
        );
    }
}

fn motion_of(m: &Model, h: Handle<Surface>) -> Option<Handle<nacre_topo::MotionNode>> {
    match m.surface_truth(h) {
        SurfaceTruth::Plane { motion, .. } => *motion,
        SurfaceTruth::Cylinder { motion } => *motion,
    }
}

fn live_vertices(m: &Model) -> Vec<Handle<Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &s in &m.live_solids {
        let sol = m.solids.get(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                let f = m.faces.get(fh);
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        let vh = m.he_start(he);
                        if seen.insert(vh) {
                            out.push(vh);
                        }
                    }
                }
            }
        }
    }
    out
}

/// Walk every live vertex and ask both solves what its base coordinate is.
fn measure(m: &Model) -> Tally {
    let mut t = Tally::default();
    for vh in live_vertices(m) {
        t.vertices += 1;
        let tol = m.vertex_tol(vh);
        match tol {
            None => t.constructed += 1,
            Some(0.0) => t.discovered_exact += 1,
            Some(_) => t.discovered_inexact += 1,
        }

        let VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
            t.on_seam += 1;
            continue;
        };

        let names: Vec<&PlaneName> = tri.iter().filter_map(|h| m.surface_name.get(h)).collect();
        if names.len() != 3 {
            t.unnamed_carrier += 1;
            continue;
        }
        if names.iter().any(|n| n.narrow().is_none()) {
            t.wide_carrier += 1;
        }

        let (a, b, c) = (
            motion_of(m, tri[0]),
            motion_of(m, tri[1]),
            motion_of(m, tri[2]),
        );
        if !(a == b && b == c) {
            // Structural, not width: the exact road needs one shared frame to replay from.
            t.mixed_motion += 1;
            continue;
        }

        let Some(point) = nacre_scalar::three_planes_big([names[0], names[1], names[2]]) else {
            t.singular += 1;
            continue;
        };
        let w = point.width_bits();
        t.widths.push(w);
        match tol {
            Some(0.0) => t.widths_exact.push(w),
            Some(_) => t.widths_inexact.push(w),
            None => {}
        }

        // ★★ The split the whole measurement turns on. The narrow route only runs when every
        // carrier is narrow; where it declines anyway, the twin says which fact that was.
        let narrow_in: Option<[[Rat; 4]; 3]> = (|| {
            let mut rows = [[Rat::from_int(0); 4]; 3];
            for (o, n) in rows.iter_mut().zip(&names) {
                *o = *n.narrow()?;
            }
            Some(rows)
        })();
        if let Some(rows) = narrow_in {
            match nacre_scalar::three_planes_rat(rows) {
                Some(got) => {
                    t.narrow_answered += 1;
                    assert_eq!(
                        MeetPoint::Narrow(got),
                        point,
                        "the two solves disagree on a vertex the narrow one answered"
                    );
                }
                None if w <= 127 => t.declined_but_fits += 1,
                None => t.declined_and_wide += 1,
            }
        }
    }
    t
}

// ---------------------------------------------------------------------------------------------
// populations
// ---------------------------------------------------------------------------------------------

/// Two axis-aligned cuboids cut — the plainest way to make a discovered corner.
fn boolean_corner() -> Model {
    let mut m = Model::new();
    let a = cuboid(&mut m, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let b = cuboid(&mut m, [3.3, 3.3, -1.0], [7.7, 7.7, 11.0]);
    boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();
    m
}

/// The same, then cut again — does depth widen the base?
fn boolean_twice() -> Model {
    let mut m = boolean_corner();
    let a = m.live_solids[0];
    let b = cuboid(&mut m, [-1.0, 4.4, 4.4], [11.0, 8.8, 8.8]);
    boolean(&mut m, BoolKind::Cut, a, b).expect("second cut");
    m.rebuild_adjacency();
    m
}

/// A turned operand — carriers share one motion, so the exact road applies and the base is
/// solved in the frame the planes are stated in.
fn boolean_rotated() -> Model {
    let mut m = Model::new();
    let a = cuboid(&mut m, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let b = cuboid(&mut m, [3.3, 3.3, -1.0], [7.7, 7.7, 11.0]);
    let iso = Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
    });
    let OpOutput::Transform { solid: b } = apply(
        &mut m,
        &Operation::Transform {
            solid: b,
            isometry: iso,
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();
    m
}

/// **The census `wf` family**, verbatim in shape: a prism on a tilted decimal orthonormal frame,
/// then a small pocket on one of its walls. That family qualifies itself on a wall whose name is
/// narrow but whose squared lengths overflow — the population the narrow frame derivation
/// declines — which is exactly where a plane's coefficients stop being small integers.
///
/// ★ The obvious fixture (a tilted prism cut through an axis-aligned box) is **not** usable here:
/// it rejects with `JudgeExhausted`. A feature on a face is the well-conditioned way to reach
/// this population, and it is the one the corpus already uses.
///
/// ★★ `passes` is the **depth** axis that matters. `boolean_twice` varies depth too, but every
/// input there is a short decimal, so a second cut has nothing to widen — it measures depth on a
/// population that cannot show the effect. Depth has to be varied where the coordinates are
/// already awkward, which is here.
fn tilted_frame(passes: usize) -> Model {
    let mut m = Model::new();
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(0.1111111111111111, 0.1234567890123456),
                p2(4.123456789012345, 0.2345678901234567),
                p2(3.9876543210987654, 3.1234567890123459),
                p2(0.2222222222222222, 2.765432109876543),
            ])
            .unwrap(),
            dist: 2.5,
        },
    )
    .expect("the wf base prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    let mut done: Vec<Handle<Surface>> = Vec::new();
    for pass in 0..passes {
        let live = m.live_solids[0];
        let wall = *m
            .shells
            .get(m.solids.get(live).outer)
            .faces
            .iter()
            .find(|&&f| {
                let s = m.faces.get(f).surface;
                !done.contains(&s)
                    && m.surface_name
                        .get(&s)
                        .and_then(|n| n.narrow())
                        .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
            })
            .unwrap_or_else(|| panic!("the wf population vanished on pass {pass}"));
        done.push(m.faces.get(wall).surface);
        let sp = nacre_ops::face_plane(&m, wall).expect("planar");
        let d = nacre_props::face_props(&m, wall).unwrap().centroid - sp.origin();
        let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
        apply(
            &mut m,
            &Operation::PocketOnFace {
                face: wall,
                profile: Profile2d::polygon(vec![
                    p2(cu - 0.3, cv - 0.3),
                    p2(cu + 0.3, cv - 0.3),
                    p2(cu + 0.3, cv + 0.3),
                    p2(cu - 0.3, cv + 0.3),
                ])
                .unwrap(),
                dist: 0.4,
            },
        )
        .unwrap_or_else(|e| panic!("the wf pocket, pass {pass}: {e:?}"));
        m.rebuild_adjacency();
    }
    m
}

// ---------------------------------------------------------------------------------------------
// the measurements
// ---------------------------------------------------------------------------------------------

/// ★★★★★ **The headline number: how wide is a discovered coordinate's base, really.**
///
/// Printed, not asserted. A corpus maximum is not a bound, and freezing one here would turn a
/// sample into a spec. What *is* asserted is that the instrument moved at all — a table of zeros
/// and a dead probe are indistinguishable otherwise.
#[test]
#[ignore = "measurement — run explicitly, prints the table"]
fn how_wide_a_discovered_coordinate_is() {
    let mut discovered = 0usize;
    let mut solved = 0usize;
    for (what, m) in [
        ("boolean_corner", boolean_corner()),
        ("boolean_twice", boolean_twice()),
        ("boolean_rotated", boolean_rotated()),
        ("tilted_frame", tilted_frame(1)),
        ("tilted_frame_x2", tilted_frame(2)),
    ] {
        let t = measure(&m);
        t.report(what);
        discovered += t.discovered_exact + t.discovered_inexact;
        solved += t.widths.len();
    }
    assert!(
        discovered > 0,
        "no discovered vertex in any population — the fixtures never reached the target"
    );
    assert!(solved > 0, "nothing solved — the instrument is dead");
}

/// ★★★★ **The invariant, cheap enough for every run**: wherever `three_planes_rat` answers, the
/// twin answers the same. `measure` asserts it per vertex; this runs it on the population where
/// both routes are live so a wiring error cannot hide behind an `#[ignore]`.
#[test]
fn the_two_solves_agree_on_every_vertex_the_narrow_one_answers() {
    let t = measure(&boolean_corner());
    assert!(
        t.narrow_answered > 0,
        "the narrow route answered nowhere — the agreement check measured nothing"
    );
}
