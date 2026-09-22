//! **How wide a discovered coordinate actually is — measured.**
//!
//! The truth is a rational *or* a handle, and values which do not fit `Rat` are pointed at rather
//! than stored. How wide a discovered coordinate is decides whether a new plane variant is needed
//! or not. This is its measurement.
//!
//! ★ **What the instrument has to be.** `three_planes_rat` cannot answer the question it is
//! about: it answers only where the point fits `Rat`, so everything it reports fits by
//! construction. `three_planes_big` has no ceiling, and it separates the two facts a checked
//! `None` would conflate — *the point does not fit* and *an intermediate overflowed on a point
//! that would have*.
//!
//! ★★ **What is not measured here, and why.** Not the *world* coordinate — a moved point does not
//! have one as a rational. A vertex's truth is `Vertex::ThreePlane`, three surface handles and
//! **no coordinate at all**; the motion rides on the surface, and a general
//! rotation's cos/sin are irrational. So the only rational a discovered point *could* be written
//! as is the one in the frame its carriers are stated in, which is what a plane name pins and what
//! this measures.
//!
//! Populations come from the kernel's own discriminator rather than from a story: a vertex is
//! `PointCache::Bounded` when it was realized from its definition, and anything else says the
//! cache has no realization behind it, as the cache itself says. Those ride along as the free
//! control.
//!
//! `OnSeam` vertices are excluded rather than counted as failures: that variant pins a curve, not
//! a point, and its doc already records the coordinate cache as load-bearing there.
//!
//! # What it measures
//!
//! | population | vertices | solved | **width max** | over 127 | narrow declined, point fits |
//! |---|---|---|---|---|---|
//! | `boolean_corner`   | 16 | 16 |  **7** | 0 | 0 |
//! | `boolean_twice`    | 32 | 32 |  **7** | 0 | 0 |
//! | `boolean_rotated`  | 20 | 12 |  **7** | 0 | 0 |
//! | `tilted_frame`     | 16 | 12 | **59** | 0 | 0 |
//! | `tilted_frame_x2`  | 24 | 16 | **59** | 0 | 0 |
//!
//! ★★★ **Nothing is wide.** The corpus maximum is **59 bits**; a second feature on the tilted
//! population does not move it (59 → 59), and stacking forty motions leaves the widest *name*
//! exactly where it started (`stacking_operations_does_not_widen_a_name`). ★ These are **corpus
//! numbers, not bounds** — the meter's negative control lives in `nacre-exact`
//! (`the_width_meter_reports_a_point_no_rat_can_hold`, 201 bits), so `over127 = 0` is a fact
//! about this population and not about a clamped instrument. The narrow route answers every
//! solved vertex (`declined_but_fits = 0`).
//!
//! ★ **The corpus is small** — five models, 108 vertices. It is enough to exhibit counterexamples.
//! It is **not** enough to say what the widest coordinate this kernel can produce is, and a corpus
//! maximum is not a bound in any case: the width a `Rat` must hold is a property of the *type*
//! (`PlaneName` is `Narrow | Wide`, so a solved coordinate can be arbitrarily wide), and these
//! numbers only say which branch today's models take.
//!
//! # And what the f64 road does with those coordinates
//!
//! | population | solved | triples agreeing | **different plane** |
//! |---|---|---|---|
//! | `boolean_corner` | 16 | 552 | **0** (negative control) |
//! | `tilted_frame`   | 12 | 0 | **220 — every one** |
//!
//! ★★★ Spelling "the plane through those three corners" in **coordinates** produces a plane with a
//! **different name** on tilted geometry — a different handle, and exact identity answers "no".
//! That is the capability gap a vertex-naming datum closes, and it is a fact about the population
//! (the control shows the two roads agreeing on every triple), not about the probe.
//!
//! # And how much of that vocabulary the frame question reaches
//!
//! A datum needs **three** vertices, so the vertex-level column above cannot answer it. Of the
//! triples the frame question decides (`accepted + accepted_nameless + accepted_straddle`;
//! collinear and undefined ones are refused whatever happens to frames):
//!
//! | population | vertices | pure | **pure frames** | named | pure-mixed | straddle | collinear |
//! |---|---|---|---|---|---|---|---|
//! | `boolean_corner`       | 16 | 16 | 1 | **100%** | 0 | 0 | 8 |
//! | `turned_after_the_cut` | 16 | 16 | 1 | **100%** | 0 | 0 | 8 |
//! | `boolean_rotated`      | 20 | 12 | **2** | **5.3%** | 14.0% | 80.7% | 0 |
//! | `tilted_frame`         | 16 | 12 | **2** | **10.7%** | 28.6% | 60.7% | 0 |
//! | `tilted_frame_x2`      | 24 | 16 | **3** | **3.2%** | 24.5% | 72.3% | 0 |
//!
//! ★ Every geometrically sound triple in every population classifies as acceptable — named, a
//! pure-mixed statement (each vertex exact in its own frame), or a straddling vertex taken as the
//! **meet of its carriers** (the op can still refuse an individual statement as
//! `ThroughFrameUndecided`). The boolean over such a datum's own face is a separate boundary.
//!
//! ★★ **Pure-mixed is its own cause.** Rotation against the world spans the still operand's frame
//! and the turned one's; a prism on a tilted frame with a feature on a *wall* spans **two** frames
//! of pure vertices, and two passes span three.
//!
//! ★★★★★ **`turned_after_the_cut` is letter-identical to `boolean_corner`**, which is the
//! invariance a rigid motion owes. A Z-turned block's caps stay world-stated while its walls carry
//! a node — a motion that *fixes* a plane restates nothing — so every corner is two walls and one
//! cap. Read as "all three carriers share one motion", all 16 corners would be straddling; the
//! rule that tells them apart is `Model::chain_fixes_plane` (a world-stated carrier is usable in
//! the chain's frame **iff** the chain fixes it), and this file asks `Model::vertex_meet` rather
//! than restating it. That is the shape to watch for here: a counter that restates a kernel rule
//! can shrink a population silently while every assertion stays green.

use nacre_exact::{Angle, Axis, Isometry, MeetPoint, PlaneName, Rat, Rotation};
use nacre_math::{Point3, Vector3};
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane};
use nacre_ops::{apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, Solid, Surface, Vertex};

use crate::fixtures::{datum_frame, p2};

// ---------------------------------------------------------------------------------------------
// fixtures — the shapes `tests/invariants/points_coverage.rs` already uses
// ---------------------------------------------------------------------------------------------

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
    /// No realization stands behind the coordinate — the construction's own figure (the control).
    unrealized: usize,
    /// The cache road stopped: too few bits on its one rung, or a history past its cost cap.
    ceiling: usize,
    /// Realized from the definition (`Bounded`): exact where every bound is zero.
    realized_exact: usize,
    /// Realized, with a nonzero bound on some axis (an irrational or non-representable point).
    realized_inexact: usize,
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
            "stat {what:22} vertices={} unrealized={} ceiling={} realized(exact0={} inexact={}) seam={}",
            self.vertices,
            self.unrealized,
            self.ceiling,
            self.realized_exact,
            self.realized_inexact,
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
    match m.surface(h) {
        Surface::Plane { motion, .. } => *motion,
        Surface::Cylinder { motion, .. } => *motion,
    }
}

fn live_vertices(m: &Model) -> Vec<Handle<Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &s in m.live_solids() {
        let sol = m.solid(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let f = m.face(fh);
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
        match m.vertex_cache(vh) {
            PointCache::Unrealized { .. } => t.unrealized += 1,
            PointCache::Ceiling { .. } => t.ceiling += 1,
            PointCache::Bounded { bound, .. } if bound.iter().all(|b| b.is_zero()) => {
                t.realized_exact += 1
            }
            PointCache::Bounded { .. } => t.realized_inexact += 1,
        }

        let Vertex::ThreePlane(tri) = *m.vertex(vh) else {
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

        // ★★★ **The kernel says which vertices have a frame; this table only labels the rest.**
        // Comparing the three carriers' motions here would be `vertex_meet`'s rule written a
        // second time — and with the invariant-plane restatement such a copy says "mixed" for
        // every corner of a turned solid and **drops them out of the width table entirely**,
        // green while measuring a smaller population than it claims. Asking the door keeps the
        // two in step.
        let Some((point, _frame)) = m.vertex_meet(vh) else {
            // Two causes are left (the def kind and a missing name were handled above): no
            // shared frame, or — same frame, no meet — carriers on one line. The split is for
            // the table's wording; the capability answer above is the kernel's.
            let f = m.plane_motion(tri[0]);
            if tri.iter().all(|h| m.plane_motion(*h) == f) {
                t.singular += 1;
            } else {
                t.mixed_motion += 1;
            }
            continue;
        };
        let w = point.width_bits();
        t.widths.push(w);
        // Split by what the cache knows: a realization with every bound zero is exact, one with a
        // nonzero bound is not; a coordinate with no realization behind it is counted in neither.
        match m.vertex_cache(vh) {
            PointCache::Bounded { bound, .. } if bound.iter().all(|b| b.is_zero()) => {
                t.widths_exact.push(w)
            }
            PointCache::Bounded { .. } => t.widths_inexact.push(w),
            PointCache::Ceiling { .. } | PointCache::Unrealized { .. } => {}
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
            match nacre_exact::three_planes_rat(rows) {
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
// how much of the datum vocabulary the frame wall costs — measured where the caller picks
// ---------------------------------------------------------------------------------------------

/// What `DatumDef::ThroughVertices` would find when it looks at **one** vertex — solved once and
/// reused across every triple it appears in.
///
/// ★ Solving per triple would be the same `three_planes_big` call `n²` times over; solving per
/// vertex makes the triple pass plain combination arithmetic, so the whole `C(n,3)` sweep runs
/// with no cap and nothing silently dropped.
#[allow(clippy::large_enum_variant)] // a per-vertex scratch value in a measurement walk
enum VertexReach {
    /// Carriers share one motion — the frame it is written in, and the meet **at whatever width
    /// it needs** (★ a meet wider than `Rat` is not its own kind: the named road accepts it,
    /// so width is not a reach distinction — exactly mirroring the producer).
    Pure(Option<Handle<nacre_topo::MotionNode>>, MeetPoint),
    /// The door cannot place this vertex in any one frame. **The kernel makes these**, not the
    /// caller: a cut between a turned operand and a still one leaves corners where an unmoved wall
    /// meets two turned ones.
    ///
    /// ★ Not simply "its three carriers carry different motions" — a *fixed* world-stated carrier
    /// sits beside moved ones and is placeable anyway (`Model::chain_fixes_plane`); counting those
    /// here would report a kernel that does not exist.
    Straddle,
    /// `OnSeam`, or a carrier with no name.
    Undefined,
}

/// Triple-level tallies — the unit a caller actually picks in.
#[derive(Default, Debug)]
struct Reach {
    vertices: usize,
    pure: usize,
    /// ★ How many **distinct frames** the pure vertices span. This is what decides whether the
    /// all-pure combination count is achieved: `through_points_rat`'s two conditions differ in
    /// kind — carriers sharing a motion is *per vertex*, three vertices sharing a frame is
    /// *pairwise* — so a triple of pure vertices in two frames is still refused.
    pure_frames: usize,
    triples: usize,
    accepted: usize,
    /// ★ A straddling vertex does not block its triple: the meet
    /// road accepts it (the op can still refuse an individual statement with
    /// `ThroughFrameUndecided`, so like the other buckets this is the classifier's answer — an
    /// upper bound on op-level acceptance).
    accepted_straddle: usize,
    /// ★ Pure vertices in differing frames: the judged frame accepts them, namelessly. The op can
    /// still refuse an individual one (`ThroughFrameUndecided`), so this is the *classifier's*
    /// answer, same as `accepted`.
    accepted_nameless: usize,
    /// Accepted triples in which at least one meet needed the wide vessel — visibility for the
    /// population of meets too wide for `Rat`, kept as a count so
    /// "still 0 in this corpus" stays a statement the table makes rather than an assumption.
    accepted_wide: usize,
    undefined: usize,
    /// Only decidable on the accepted road: a triple that solved and named no plane.
    collinear: usize,
}

impl Reach {
    fn report(&self, what: &str) {
        // The quantity: of the triples the *frame* is what blocks, how many would be freed.
        let blocked = 0usize;
        let denom = self.accepted + self.accepted_nameless + self.accepted_straddle;
        let pct = |n: usize| {
            if denom == 0 {
                0.0
            } else {
                100.0 * n as f64 / denom as f64
            }
        };
        println!(
            "stat datum_reach {what:18} vertices={} pure={} pure_frames={} triples={}",
            self.vertices, self.pure, self.pure_frames, self.triples
        );
        let _ = blocked;
        println!(
            "stat datum_reach {what:18} accepted={} ({:.1}%) accepted_nameless={} ({:.1}%) \
             accepted_straddle={} ({:.1}%)",
            self.accepted,
            pct(self.accepted),
            self.accepted_nameless,
            pct(self.accepted_nameless),
            self.accepted_straddle,
            pct(self.accepted_straddle),
        );
        println!(
            "stat datum_reach {what:18} collinear={} accepted_wide={} undefined={}",
            self.collinear, self.accepted_wide, self.undefined
        );
    }
}

/// Walk every live vertex once, then every distinct triple, and classify each triple by **what
/// the datum operation would answer**.
///
/// ★★★ **The denominator is stated, not assumed.** `collinear` and `undefined` are refused
/// whatever happens to frames, so they are outside the question; the reported percentage is over
/// `accepted + accepted_nameless + accepted_straddle` — *"of the triples the frame question
/// decides, how is each accepted"*. ★ A nameless triple is classified before anything asks
/// whether it is also collinear, so the nameless buckets can hold triples that would fail anyway.
/// `collinear`'s rate on the named road is the only estimate of that contamination, and it is
/// printed for exactly that reason.
///
/// ★ Causes are assigned by **priority**, not by which would fire first in the producer's
/// per-vertex loop — a triple can carry more than one, and a per-triple bucket has to pick. The
/// order is the producer's own: undefined, straddle, differing frames, width, collinear.
fn datum_reach(m: &Model) -> Reach {
    let mut t = Reach::default();
    let verts = live_vertices(m);
    let reach: Vec<VertexReach> = verts
        .iter()
        .map(|&vh| {
            let Vertex::ThreePlane(tri) = *m.vertex(vh) else {
                return VertexReach::Undefined;
            };
            let names: Vec<&PlaneName> = tri.iter().filter_map(|h| m.surface_name.get(h)).collect();
            if names.len() != 3 {
                return VertexReach::Undefined;
            }
            // ★★★ **Ask the door, do not re-derive its rule.** `VertexReach::Pure` carries
            // exactly `vertex_meet`'s return — this classifier *was* that function, written a
            // second time, and the copy is what kept reporting yesterday's kernel after the
            // invariant-plane restatement moved the line.
            match m.vertex_meet(vh) {
                Some((p, f)) => VertexReach::Pure(f, p),
                // The label only; the capability answer is the kernel's above.
                None => {
                    let f = m.plane_motion(tri[0]);
                    if tri.iter().all(|h| m.plane_motion(*h) == f) {
                        VertexReach::Undefined // carriers share a line — no unique point
                    } else {
                        VertexReach::Straddle
                    }
                }
            }
        })
        .collect();

    t.vertices = verts.len();
    let mut frames = std::collections::HashSet::new();
    for r in &reach {
        if let VertexReach::Pure(f, _) = r {
            t.pure += 1;
            frames.insert(*f);
        }
    }
    t.pure_frames = frames.len();

    let n = verts.len();
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                t.triples += 1;
                let three = [&reach[i], &reach[j], &reach[k]];
                if three.iter().any(|r| matches!(r, VertexReach::Undefined)) {
                    t.undefined += 1;
                } else if three.iter().any(|r| matches!(r, VertexReach::Straddle)) {
                    t.accepted_straddle += 1;
                } else {
                    let f = three.iter().map(|r| match r {
                        VertexReach::Pure(f, _) => *f,
                        _ => unreachable!("the other kinds were taken above"),
                    });
                    let pts: Vec<&MeetPoint> = three
                        .iter()
                        .map(|r| match r {
                            VertexReach::Pure(_, p) => p,
                            _ => unreachable!("the other kinds were taken above"),
                        })
                        .collect();
                    if f.clone().collect::<std::collections::HashSet<_>>().len() != 1 {
                        t.accepted_nameless += 1;
                    } else if nacre_exact::plane_name_from_meets([pts[0], pts[1], pts[2]]).is_none()
                    {
                        t.collinear += 1;
                    } else {
                        t.accepted += 1;
                        t.accepted_wide +=
                            usize::from(pts.iter().any(|p| matches!(p, MeetPoint::Wide(_))));
                    }
                }
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

/// ★★★★ **The negative control with teeth: motion everywhere, and all of it shared.**
///
/// `boolean_corner` reports zero mixed vertices, but it *cannot* report anything else — it has no
/// motion at all, so the instrument is not being asked. Turning the finished solid gives every
/// surface a motion (`transform_solid` memoizes one new node per distinct parent leaf, and every
/// parent here is `None`, so they all land on one) while leaving the geometry alone. A nonzero
/// mixed count here would mean either the measurement is reading something other than "do the
/// frames agree", or those motions did not intern to one node — and either is a finding.
///
/// ★ The source stays live after a `Transform` (it is not consumed the way a boolean's operands
/// are), and a sweep over *both* would readmit the caller-error shape — two solids in two frames —
/// which is precisely the variable this control is holding still. So the live set is narrowed to
/// the image.
fn turned_after_the_cut() -> Model {
    let mut m = boolean_corner();
    let src = m.live_solids()[0];
    let OpOutput::Transform { solid } = apply(
        &mut m,
        &Operation::Transform {
            solid: src,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
            }),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.restore_live(vec![solid]);
    m.rebuild_adjacency();
    m
}

/// The same, then cut again — does depth widen the base?
fn boolean_twice() -> Model {
    let mut m = boolean_corner();
    let a = m.live_solids()[0];
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
        pivot: [Rat::from_int(0); 3],
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
        let live = m.live_solids()[0];
        let wall = *m
            .shell(m.solid(live).outer)
            .faces
            .iter()
            .find(|&&f| {
                let s = m.face(f).surface;
                !done.contains(&s)
                    && m.surface_name
                        .get(&s)
                        .and_then(|n| n.narrow())
                        .is_some_and(|c| nacre_exact::plane_frame_default(*c).is_none())
            })
            .unwrap_or_else(|| panic!("the wf population vanished on pass {pass}"));
        done.push(m.face(wall).surface);
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
// (ii) — is the f64 road a different plane?
// ---------------------------------------------------------------------------------------------

/// Every vouched-for vertex — **realized from its definition** — that solves to a `Rat` triple,
/// paired with its handle.
fn solved_discovered(m: &Model) -> Vec<([Rat; 3], Handle<Vertex>)> {
    let mut out = Vec::new();
    for vh in live_vertices(m) {
        // ★ Stated positively on purpose. This read "not the bare figure" while there were two
        // ways to be vouched for; with `Ceiling` in the enum the negative form would quietly admit
        // a vertex nothing has proven anything about.
        if !matches!(m.vertex_cache(vh), PointCache::Bounded { .. }) {
            continue; // the construction's bare figure — the control lives in its own fixture below
        }
        let Vertex::ThreePlane(tri) = *m.vertex(vh) else {
            continue;
        };
        let names: Vec<&PlaneName> = tri.iter().filter_map(|h| m.surface_name.get(h)).collect();
        if names.len() != 3 {
            continue;
        }
        let ms = tri.map(|h| motion_of(m, h));
        if !(ms[0] == ms[1] && ms[1] == ms[2]) {
            continue; // mixed frames — no rational coordinate in any single frame
        }
        if let Some(p) = nacre_exact::three_planes_big([names[0], names[1], names[2]]) {
            if let Some(n) = p.narrow() {
                out.push((*n, vh));
            }
        }
    }
    out
}

/// `(same, different)` over every non-degenerate triple: the plane through the vertices' **exact**
/// coordinates vs the plane through their **f64 cache** coordinates lifted back by `from_decimal`.
fn triple_verdicts(m: &Model, solved: &[([Rat; 3], Handle<Vertex>)]) -> (usize, usize) {
    let lift = |vh: Handle<Vertex>| -> Option<[Rat; 3]> {
        let a = m.vertex_point(vh).as_array().map(Rat::from_decimal);
        Some([a[0]?, a[1]?, a[2]?])
    };
    let (mut same, mut diff) = (0, 0);
    for i in 0..solved.len() {
        for j in (i + 1)..solved.len() {
            for k in (j + 1)..solved.len() {
                let exact = nacre_exact::plane_name_exact(solved[i].0, solved[j].0, solved[k].0);
                let rounded = (|| {
                    nacre_exact::plane_name_exact(
                        lift(solved[i].1)?,
                        lift(solved[j].1)?,
                        lift(solved[k].1)?,
                    )
                })();
                match (exact, rounded) {
                    (Some(a), Some(b)) if a == b => same += 1,
                    (Some(_), Some(_)) => diff += 1,
                    _ => {} // collinear either way — carries no verdict
                }
            }
        }
    }
    (same, diff)
}

/// ★★★★★ **The capability gap, stated as a counterexample.**
///
/// "A plane through those three corners" is an ordinary CAD request, and today it can only be
/// spelled in coordinates. Where a vertex's exact coordinate is not an `f64` — and
/// the cache *is* its nearest `f64`, which is the closest a coordinate can come — the plane that
/// spelling produces is **a different plane**: not a nearby one, a different name, which interns
/// to a different handle and answers exact identity with "no".
///
/// ★ The assertion is **existence of a counterexample**, never "always different": a corpus
/// cannot carry a universal. The axis-aligned fixture is the negative control that keeps the
/// claim about the *population* rather than about the probe — there the two roads agree, because
/// there the cache **is** the exact coordinate.
#[test]
fn the_f64_road_names_a_different_plane() {
    let tilted = tilted_frame(1);
    let solved = solved_discovered(&tilted);
    let (same, diff) = triple_verdicts(&tilted, &solved);
    let exact_cache = solved
        .iter()
        .filter(|(r, vh)| {
            let c = tilted.vertex_point(*vh).as_array();
            (0..3).all(|t| r[t].to_f64() == c[t])
        })
        .count();
    println!(
        "stat ii tilted_frame     solved={} cache_is_exact={exact_cache} same={same} DIFFERENT={diff}",
        solved.len()
    );
    assert!(
        diff > 0,
        "no counterexample: every triple agreed, so this fixture cannot show the gap"
    );

    // ★ The negative control. Axis-aligned boolean corners land on exact f64s, so the same
    // comparison must come back *agreeing* — otherwise "different" is a property of the probe.
    let square = boolean_corner();
    let solved = solved_discovered(&square);
    let (same, diff) = triple_verdicts(&square, &solved);
    println!(
        "stat ii boolean_corner   solved={} same={same} DIFFERENT={diff}",
        solved.len()
    );
    assert!(same > 0, "the control measured nothing");
    assert_eq!(
        diff, 0,
        "the two roads disagreed where the cache is exact — the probe, not the population"
    );
}

// ---------------------------------------------------------------------------------------------
// accumulation — does repeated operation widen anything?
// ---------------------------------------------------------------------------------------------

/// The widest plane **name** in the model, and how the population is shaped.
fn name_census(m: &Model) -> (u64, usize, usize, usize) {
    let mut widest = 0;
    let mut moved = 0;
    let mut distinct = std::collections::HashSet::new();
    for (h, n) in m.surface_name.iter() {
        let w = match n {
            PlaneName::Narrow(c) => c
                .iter()
                .flat_map(|r| [r.numer(), r.denom()])
                .map(|v| (128 - v.unsigned_abs().leading_zeros()) as u64)
                .max()
                .unwrap_or(0),
            PlaneName::Wide(c) => c.iter().map(|x| x.bits()).max().unwrap_or(0),
        };
        widest = widest.max(w);
        distinct.insert(format!("{n:?}"));
        if let Surface::Plane {
            motion: Some(_), ..
        } = m.surface(*h)
        {
            moved += 1;
        }
    }
    (widest, m.surface_name.len(), distinct.len(), moved)
}

fn moved(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    )
    .expect("transform") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// ★★★★★ **«Operations accumulate, so surely it overflows eventually» — measured, and no.**
///
/// The width of a stored name does not depend on how many operations preceded it. Forty stacked
/// motions leave the maximum name width **exactly where it started**, and the reason is visible in
/// the same table: the surface count climbs by six per operation while the count of **distinct
/// names stays put**. Every moved face is the *same name* under a new motion handle — "the face
/// carries the motion", doing precisely what it says.
///
/// Three facts close the question together, and the other two are structural rather than
/// measured here:
/// - a boolean never calls `push_plane` at all (result faces reuse the operands' planes), so the
///   operation that composes most does not create coefficients;
/// - `plane_offset` fixes `n` and moves only `d`, and decimals added together share the factor
///   ten, so denominators meet at an lcm rather than multiplying (`1.1 + 6.6 = 7.7`).
///
/// ⇒ Width is set by **the geometry someone wrote down**, not by how much was done to it. It is a
/// single hop, and `Rat::from_decimal`'s window is what bounds that hop.
#[test]
#[ignore = "measurement — run explicitly, prints the table"]
fn stacking_operations_does_not_widen_a_name() {
    let dec = |x: f64| Rat::from_decimal(x).unwrap();
    let mut first = None;
    let mut last = None;

    for (what, turn) in [("translate", false), ("rot90_translate", true)] {
        let mut m = Model::new();
        let mut s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 2.0, 3.0]),
        );
        for i in 1..=40 {
            if turn {
                s = moved(
                    &mut m,
                    s,
                    Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        pivot: [dec(0.1), dec(0.7), Rat::from_int(0)],
                        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
                    }),
                );
            }
            s = moved(
                &mut m,
                s,
                Isometry::translation([dec(0.13), dec(0.3), dec(0.7)]),
            );
            if i % 20 == 0 {
                let (w, surfaces, names, with_node) = name_census(&m);
                println!(
                    "stat {what:18} step={i:3} width={w} surfaces={surfaces} distinct_names={names} with_motion_node={with_node}"
                );
                if what == "translate" {
                    if i == 20 {
                        first = Some(w);
                    } else {
                        last = Some(w);
                    }
                }
            }
        }
    }
    assert_eq!(
        first, last,
        "the name width moved between step 20 and step 40 — accumulation is not free after all"
    );
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
        discovered += t.realized_exact + t.realized_inexact;
        solved += t.widths.len();
    }
    assert!(
        discovered > 0,
        "no discovered vertex in any population — the fixtures never reached the target"
    );
    assert!(solved > 0, "nothing solved — the instrument is dead");
}

/// ★★★★★ **How much of the datum vocabulary the frame wall costs — in the unit a caller picks.**
///
/// The vertex-level number (12 of 20 mixed on `boolean_rotated`) cannot answer this: a datum needs
/// **three** vertices and `through_points_rat`'s two conditions differ in kind — carriers sharing
/// one motion is *per vertex*, the three vertices sharing one frame is *pairwise*. So the
/// all-pure combination count is only an upper bound, and `pure_frames` is what says whether it is
/// reached.
///
/// ★★ **Two negative controls, because one of them cannot fail.** `boolean_corner` has no motion,
/// so "zero mixed" there is true by construction and measures nothing.
/// [`turned_after_the_cut`] has motion on every surface that the turn moves, and every carrier
/// solvable in that one frame — the caps stay world-stated because the turn **fixes** them, and a
/// fixed plane's world equation is its pre-motion equation. That is the arm where a zero is
/// evidence, and it is the arm that went red for 169 commits when the kernel read the caps'
/// `motion: None` as a straddle (see the module doc's regression section).
/// ★★★ **Not `#[ignore]`d, unlike its two neighbours** — this one carries an *assertion*, and
/// hiding an assertion behind a flag the commit hook never passes is how the regression above
/// lived for 169 commits. It costs ~0.3s and its table is captured unless it fails, so the price
/// of the hook seeing it is nothing. The other two measurements print and assert nothing, so they
/// stay where they are.
#[test]
fn how_much_of_the_datum_vocabulary_the_frame_wall_costs() {
    let mut blocked_somewhere = false;
    for (what, m) in [
        ("boolean_corner", boolean_corner()),
        ("turned_after_the_cut", turned_after_the_cut()),
        ("boolean_rotated", boolean_rotated()),
        ("tilted_frame", tilted_frame(1)),
        ("tilted_frame_x2", tilted_frame(2)),
    ] {
        let t = datum_reach(&m);
        t.report(what);
        if matches!(what, "boolean_corner" | "turned_after_the_cut") {
            assert_eq!(
                t.accepted_straddle + t.accepted_nameless,
                0,
                "{what}: every carrier is solvable in one frame here, so nothing may be mixed \
                 at all — the turn fixes the caps, and a fixed plane's world equation is its \
                 pre-motion equation (`Model::chain_fixes_plane`)"
            );
            assert!(t.accepted > 0, "{what}: the control accepted nothing");
        } else {
            blocked_somewhere |= t.accepted_nameless + t.accepted_straddle > 0;
        }
    }
    assert!(
        blocked_somewhere,
        "no population exercised the judged roads — the fixtures stopped measuring them"
    );
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
