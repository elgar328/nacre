//! **What a wide-named datum costs the judge — open item 1, half of it.**
//!
//! `docs/truth-and-cache.md` §판정 says a rational-closure datum judges exactly and for free:
//! *"`Wide` 이름이 곧 정확 계수라 분모 털어 `Expansion` 으로 정확 판정 … 공짜"*.
//!
//! ★★★ **When this file was written that route did not exist** — the predicates asked for
//! coefficients through `exact_coeffs()`/`base_coeffs()`, both of which read
//! `PlaneName::narrow()`, so a `Wide` name answered `None` to every one and the judgement took
//! the toleranced route that climbs.
//!
//! ☑ **It exists now** (re-read 2026-09-09): `nacre_cip::predicate::NameInts { ints, wide }` and
//! `Judge::name_rescue` carry a wide name's integer coefficients into the three sign predicates
//! that read coefficients at all — `orient3d_cheap`, `cmp_coord`, `plane_pair_dir_sign` — and the
//! gate is written so that an all-narrow question keeps its existing route to the bit.
//! `planes_coplanar` is **not** a gap in that rescue: it runs *before* a plane table exists and
//! asks only for a [`Witness`], deliberately deciding on the faces' own coordinates rather than on
//! derived coefficients. ⇒ **the decision this file's header asks for has been made and shipped**,
//! and what the numbers below measure is not "the missing exact route" but the cost of the
//! toleranced route an **irrational motion** still takes. Re-measured today: `narrow_name` 246
//! climbs against `wide_name` 475 — **1.9×**, not the 3.2× recorded below.
//!
//! So the number below is not "how many `Expansion` pieces": it is **what the missing exact route
//! costs**, and it is the input to a decision the doc does not currently list — whether the
//! predicates should learn to read `Wide` coefficients at all.
//!
//! ★★★ **Run it with `--test-threads=1`, or the labels lie:**
//!
//! ```text
//! cargo test -p nacre-ops --test wide_datum_cost -- --ignored --nocapture --test-threads=1
//! ```
//!
//! `climb_census::take()` is a **process-global** take-and-reset, and both measurements here
//! bracket their own boolean with it. Run in parallel — the default — one test's `take()` walks
//! off with the other's accumulated climbs, and every printed number is still plausible while
//! being attributed to the wrong arm. Measured 2026-08-20: the same tree printed
//! `narrow_name=415 / wide_name=0` on one parallel run and `narrow_name=414 / wide_name=917` on
//! the next, which is how a comparison against a saved baseline invents a change that never
//! happened. Serialized, the same tree prints the same table twice.
//!
//! ★ What is **not** measured here: the irrational-motion branch (homogeneous lifting, degree
//! ~12). Its machinery is what S5(ii)-2 builds, so its cost cannot be taken before it exists —
//! the doc's "measure the cost, then build" ordering does not apply to that half.
//!
//! # What it measured (2026-08-08)
//!
//! **The population exists, and the datum is its first producer** — as open item 2b predicted.
//! Of 60 triples over the `wf` family's discovered vertices, **24 name a plane the wide vessel has
//! to hold**, the widest at **168 bits**. No construction path had ever produced one (`stat
//! wide_planes` is 0 across the corpus).
//!
//! ★ And the framed prism is the *wrong* place to look — 4 triples, **0 wide**, 5 bits. Inside its
//! own frame a tilted prism is axis-aligned, so its far cap is `w = dist`. Frames exist to make
//! that true; the width lives in what a *boolean* discovers, not in what a frame states.
//!
//! | the datum's name | climbs | mean bits | exhausted |
//! |---|---|---|---|
//! | `Narrow` | 278 | 256 | 0 |
//! | **`Wide`** | **887** | 256 | 0 |
//!
//! ⇒ **~3.2× the escalations** when the exact shortcuts decline. Both cuts still succeed and
//! nothing exhausts the budget, so this is a cost, not a cliff.
//!
//! ★★ **The residual confound, stated:** the two arms pick *different* vertex triples (that is
//! what makes one name wide and the other narrow), so their tools stand in different places. The
//! shape, size and motion history of the target are shared and the tool profile is identical, but
//! this is not a one-variable experiment — read 3.2× as "the same order of magnitude, several
//! times", not as a coefficient. A one-variable version would need one plane whose name is wide
//! and narrow at once, which is a contradiction.
//!
//! ★ Two earlier versions of this measurement reported **nothing** and looked fine doing it:
//! "datum vs no datum" gave identical climb counts, because a datum nobody cuts with never reaches
//! a plane table; and counting bits only on the resolved return path gave a mean of exactly 0,
//! because most climbs end in `Coincident` or a proved zero.

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{BoolKind, DatumDef, OpOutput, Operation, Profile2d, SketchFrame, SketchPlane};
use nacre_ops::{apply, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, Vertex};

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

fn live_verts(m: &Model) -> Vec<Handle<Vertex>> {
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

/// A tilted prism whose walls are written against a recorded frame (`n·n = 3` is not a perfect
/// square, so `exact_frame` declines), plus a box to cut against.
fn tilted_prism(m: &mut Model, dist: f64) -> Handle<Solid> {
    let plane = SketchPlane::from_origin_normal(
        Point3::from_array([0.1, 0.2, 0.3]),
        Vector3::from_array([1.0, 1.0, 1.0]),
    )
    .expect("a tilted plane");
    let frame = datum_frame(m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile: Profile2d::polygon(vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 2.0),
                p2(0.0, 2.0),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("a prism on the tilted frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// ★★★★★ **First: is there a population at all?**
///
/// A vertex-named datum's coefficients come from three solved points, so they are roughly three
/// times as wide as a coordinate — the prediction is that this is the kernel's first producer of a
/// `Wide` *name*, which open item 2b anticipated (*"the first producer is the datum (S5)"*) and no
/// construction path has ever reached (`stat wide_planes` is 0 across the corpus).
///
/// If it turns out narrow, this file stops here and says so rather than inventing a cost for a
/// population that does not exist.
#[test]
fn does_a_vertex_named_datum_produce_a_wide_name() {
    let mut m = Model::new();
    let _ = tilted_prism(&mut m, 1.7);
    let far_cap_verts: Vec<_> = live_verts(&m)
        .into_iter()
        .filter(|v| {
            let nacre_topo::VertexDef::ThreePlane(tri) = m.vertices.get(*v).def else {
                return false;
            };
            tri.iter().all(|h| m.plane_motion(*h).is_some())
        })
        .collect();
    assert!(
        far_cap_verts.len() >= 3,
        "the fixture must have vertices whose carriers share a recorded frame"
    );

    let mut widths = Vec::new();
    let mut wide = 0usize;
    for i in 0..far_cap_verts.len() {
        for j in (i + 1)..far_cap_verts.len() {
            for k in (j + 1)..far_cap_verts.len() {
                let mut t = [far_cap_verts[i], far_cap_verts[j], far_cap_verts[k]];
                t.sort_by_key(|v| v.index());
                let Some(name) = m.plane_name_through(t) else {
                    continue;
                };
                let w = match &name {
                    nacre_scalar::PlaneName::Narrow(c) => {
                        c.iter()
                            .flat_map(|r| [r.numer(), r.denom()])
                            .map(|v| 128 - v.unsigned_abs().leading_zeros())
                            .max()
                            .unwrap_or(0) as u64
                    }
                    nacre_scalar::PlaneName::Wide(c) => {
                        wide += 1;
                        c.iter().map(|x| x.bits()).max().unwrap_or(0)
                    }
                };
                widths.push(w);
            }
        }
    }
    widths.sort_unstable();
    println!(
        "stat datum_name framed_prism triples={} wide={wide} width_min={:?} width_max={:?}",
        widths.len(),
        widths.first(),
        widths.last()
    );
    assert!(!widths.is_empty(), "no nameable triple — nothing measured");

    // ★★★ The framed prism is the *wrong* population and the number says why: inside its own
    // frame the tilted geometry is axis-aligned, so the far cap is `w = dist` and its name is a
    // handful of bits. Frames exist to make exactly that true.
    //
    // The widest coordinates this kernel produces came from the decimal `wf` family's **discovered**
    // vertices (59 bits, `tests/point_width.rs`). Three of those should give a name near three
    // times as wide — that is the population to ask.
    let m = wf_family_with_pocket();
    let mut widths = Vec::new();
    let mut wide = 0usize;
    let solvable: Vec<_> = live_verts(&m)
        .into_iter()
        .filter(|v| m.vertex_tol(*v).is_some())
        .collect();
    for i in 0..solvable.len() {
        for j in (i + 1)..solvable.len() {
            for k in (j + 1)..solvable.len() {
                let mut t = [solvable[i], solvable[j], solvable[k]];
                t.sort_by_key(|v| v.index());
                let Some(name) = m.plane_name_through(t) else {
                    continue;
                };
                let w = match &name {
                    nacre_scalar::PlaneName::Narrow(c) => {
                        c.iter()
                            .flat_map(|r| [r.numer(), r.denom()])
                            .map(|v| 128 - v.unsigned_abs().leading_zeros())
                            .max()
                            .unwrap_or(0) as u64
                    }
                    nacre_scalar::PlaneName::Wide(c) => {
                        wide += 1;
                        c.iter().map(|x| x.bits()).max().unwrap_or(0)
                    }
                };
                widths.push(w);
            }
        }
    }
    widths.sort_unstable();
    println!(
        "stat datum_name wf_discovered triples={} wide={wide} width_min={:?} width_max={:?}",
        widths.len(),
        widths.first(),
        widths.last()
    );
}

/// The census `wf` family — a prism on a tilted decimal orthonormal frame with a pocket on one
/// wall, whose discovered vertices carry the widest coordinates this kernel makes (59 bits).
fn wf_family_with_pocket() -> Model {
    let mut m = Model::new();
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
        Vector3::from_array([0.6, 0.8, 0.0]),
        Vector3::from_array([-0.48, 0.36, 0.8]),
    );
    let frame = datum_frame(&mut m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
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
    let wall = *m
        .shells
        .get(m.solids.get(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            let s = m.faces.get(f).surface;
            m.surface_name
                .get(&s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
        })
        .expect("the wf population");
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
    .expect("the wf pocket");
    m.rebuild_adjacency();
    m
}

/// ★★★ **What the climb costs, split inside one model.**
///
/// Two models would have to be matched shape for shape; splitting *within* one leaves the shape,
/// the size and the motion history shared, so the only difference is whether the datum is there.
/// The counters live in `escalate`, the single funnel every production escalation passes.
///
/// ★ Reported, not asserted against a threshold — a corpus number is not a budget. What is
/// asserted is that the instrument moved at all.
#[test]
#[ignore = "measurement — run explicitly, prints the table"]
fn what_a_datum_bearing_boolean_costs() {
    use nacre_cip::kernel::frame3::climb_census;

    // ★★★ The control is **wide name vs narrow name**, not "datum vs no datum".
    //
    // The same fixture offers both (60 triples: 24 wide, 36 narrow), so the two runs share the
    // shape, the size and the motion history and differ only in whether the datum's canonical
    // name fits `Rat` — which is exactly the fork that decides whether the exact shortcuts run.
    // "Datum vs no datum" was the first attempt and it measured nothing: the two arms reported
    // identical climb counts, because a datum nobody cuts with never reaches a plane table.
    let run = |mode: &str| -> Option<(u64, u64, u64)> {
        let mut m = wf_family_with_pocket();
        let target = m.live_solids[0];
        // ★ The named arms measure over **discovered** vertices (their original question was the
        // width of discovered coordinates). The nameless arm lifts that filter: this family's
        // differ-population lives on *constructed* corners (prism corners in one frame, pocket
        // corners in another) — the discovered ones here all straddle, which stage 16-1 still
        // refuses, and an arm that only ever hit the reject would be an UNBUILDABLE lie.
        let solvable: Vec<_> = live_verts(&m)
            .into_iter()
            .filter(|v| mode == "nameless" || m.vertex_tol(*v).is_some())
            .collect();

        let mut tool = None;
        'pick: for i in 0..solvable.len() {
            for j in (i + 1)..solvable.len() {
                for k in (j + 1)..solvable.len() {
                    let mut t = [solvable[i], solvable[j], solvable[k]];
                    t.sort_by_key(|v| v.index());
                    match (mode, m.plane_name_through(t)) {
                        ("narrow_name", Some(n)) if n.narrow().is_some() => {}
                        ("wide_name", Some(n)) if n.narrow().is_none() => {}
                        // No name at all — the judged road; `apply` below still decides
                        // acceptance (a straddling triple rejects and the loop moves on).
                        ("nameless", None) => {}
                        _ => continue,
                    }
                    let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
                        &mut m,
                        &Operation::DatumPlane {
                            def: DatumDef::ThroughVertices([solvable[i], solvable[j], solvable[k]]),
                        },
                    ) else {
                        continue;
                    };
                    // The tool is raised **on** the datum, so the datum's plane is one of the
                    // boolean's operand faces — the only way it reaches a judgement at all.
                    if let Ok(OpOutput::Extrude { solid, .. }) = apply(
                        &mut m,
                        &Operation::Extrude {
                            frame,
                            profile: Profile2d::polygon(vec![
                                p2(-0.4, -0.4),
                                p2(0.4, -0.4),
                                p2(0.4, 0.4),
                                p2(-0.4, 0.4),
                            ])
                            .unwrap(),
                            dist: 1.2,
                        },
                    ) {
                        m.rebuild_adjacency();
                        tool = Some(solid);
                        break 'pick;
                    }
                }
            }
        }
        let tool = tool?;
        let _ = climb_census::take();
        let out = boolean(&mut m, BoolKind::Cut, target, tool);
        let (c, b, e) = climb_census::take();
        println!(
            "     (cut {} — {:?})",
            if out.is_ok() { "ok" } else { "declined" },
            out.err()
        );
        Some((c, b, e))
    };

    let mut seen = 0;
    for label in ["narrow_name", "wide_name", "nameless"] {
        match run(label) {
            Some((climbs, bits, exhausted)) => {
                seen += 1;
                println!(
                    "stat cost {label:12} climbs={climbs} mean_bits={} exhausted={exhausted}",
                    bits.checked_div(climbs).unwrap_or(0)
                );
            }
            // ★ Said, not silently skipped: an arm that could not be built is a fact about the
            // population, and pretending the comparison happened would be worse than no number.
            None => println!("stat cost {label:12} UNBUILDABLE — no such datum carried a tool"),
        }
    }
    assert!(seen > 0, "no arm built — the measurement is empty");
}

/// ★★★ **What the name-integer rescue actually buys, measured where its gates open.**
///
/// The first datum-bearing boolean (the table above) is a **mixed-frame** table: the wide name
/// speaks for the world (its carriers are discovered vertices — no chain), the target's walls
/// carry the wf frame, and the tool's walls carry the datum-derived `FrameWide` — no question
/// pairing them lives in one frame, so the integer gates stay shut there *by design* (and the
/// `FrameWide` walls' frame is irrational, so no integer description of those questions
/// exists to be found). The populations the rescue serves are the **second generation**:
///
/// - a solid whose faces are all world-named with the wide slice among them, cut again by a
///   world tool — the world gate;
/// - the same solid moved whole (rotate/mirror), cut by a tool moved with it — out of this
///   measurement's scope (the unit differential covers its correctness).
///
/// Two arms, one shape: the cuboid sliced by a **wide**-named datum plane vs a **narrow**-named
/// one, each then cut by the same world column. Before the rescue the wide arm had no exact
/// route at all; after it, the wide arm's cap-vs-column questions answer from the name's
/// integers. (The narrow arm's cap has a name too, but its *witness* is deep-rational, so its
/// f64 shortcuts may decline on rounding — the coefficient-mismatch population the plan lists
/// as its own follow-up. The gate deliberately requires a wide participant, so the narrow arm
/// measures the status quo.)
#[test]
#[ignore = "measurement — run explicitly, prints the table"]
fn what_a_second_generation_boolean_costs() {
    use nacre_cip::kernel::frame3::climb_census;

    let run = |want_wide: bool| -> Option<(u64, u64, u64)> {
        let mut m = wf_family_with_pocket();
        // A datum through the first triple of the wanted width, carrying a slab tool so large
        // that only the datum plane itself reaches the cuboid below.
        let solvable: Vec<_> = live_verts(&m)
            .into_iter()
            .filter(|v| m.vertex_tol(*v).is_some())
            .collect();
        let mut slab = None;
        'pick: for i in 0..solvable.len() {
            for j in (i + 1)..solvable.len() {
                for k in (j + 1)..solvable.len() {
                    let mut t = [solvable[i], solvable[j], solvable[k]];
                    t.sort_by_key(|v| v.index());
                    match m.plane_name_through(t) {
                        Some(n) if n.narrow().is_none() == want_wide => {}
                        _ => continue,
                    }
                    let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
                        &mut m,
                        &Operation::DatumPlane {
                            def: DatumDef::ThroughVertices([solvable[i], solvable[j], solvable[k]]),
                        },
                    ) else {
                        continue;
                    };
                    if let Ok(OpOutput::Extrude { solid, .. }) = apply(
                        &mut m,
                        &Operation::Extrude {
                            frame,
                            profile: Profile2d::polygon(vec![
                                p2(-50.0, -50.0),
                                p2(50.0, -50.0),
                                p2(50.0, 50.0),
                                p2(-50.0, 50.0),
                            ])
                            .unwrap(),
                            dist: 50.0,
                        },
                    ) {
                        m.rebuild_adjacency();
                        slab = Some(solid);
                        break 'pick;
                    }
                }
            }
        }
        let slab = slab?;
        let cub = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([5.0, 5.0, 4.0]),
        );
        let sliced = *nacre_ops::boolean(&mut m, BoolKind::Cut, cub, slab)
            .ok()?
            .first()?;
        m.rebuild_adjacency();
        // The second generation: every face of `sliced` is world-named, one of them by the
        // datum's name. Only this boolean is measured.
        let column = m.add_cuboid(
            Point3::from_array([0.5, 0.5, -2.0]),
            Point3::from_array([2.5, 2.5, 5.0]),
        );
        let _ = climb_census::take();
        let out = nacre_ops::boolean(&mut m, BoolKind::Cut, sliced, column);
        let (c, b, e) = climb_census::take();
        println!(
            "     (second cut {} — {:?})",
            if out.is_ok() { "ok" } else { "declined" },
            out.err()
        );
        Some((c, b, e))
    };

    let mut seen = 0;
    for (want_wide, label) in [(false, "narrow_slice"), (true, "wide_slice")] {
        match run(want_wide) {
            Some((climbs, bits, exhausted)) => {
                seen += 1;
                println!(
                    "stat cost2 {label:12} climbs={climbs} mean_bits={} exhausted={exhausted}",
                    bits.checked_div(climbs).unwrap_or(0)
                );
            }
            None => println!("stat cost2 {label:12} UNBUILDABLE — no such slice was cut"),
        }
    }
    assert!(seen > 0, "no arm built — the measurement is empty");
}
