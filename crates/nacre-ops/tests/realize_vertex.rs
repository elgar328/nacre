//! **A vertex answers for its own coordinate, at a precision the caller names.**
//!
//! The door under test is `nacre_ops::realize_vertex{,_decimal}`. What makes it worth locking is
//! not that it produces digits but that it produces *earned* ones: it realizes from the vertex's
//! definition and rounds once, where everything else in the kernel reads a cache that carries its
//! own error.
//!
//! # What the cache actually costs (measured here, `the_cache_is_not_always_nearest`)
//!
//! | population | vertex width | compared | **cache is nearest** |
//! |---|---|---|---|
//! | axis-aligned box        | 7 bits  | 16 | 16/16 |
//! | box, boolean'd twice    | 7 bits  | 32 | 32/32 |
//! | prism on a tilted frame | 59 bits | 12 | **0/12** — every one off, by up to 4 ulp |
//!
//! ⚠★★★ **The split is not representability, and saying it was is a mistake worth keeping
//! written down.** The first reading of this table said "7 bits land in an `f64` exactly, 59 do
//! not". Measured per coordinate: of `boolean_corner`'s 48, only **32** are dyadic — 16 are exact
//! rationals no `f64` can hold (the fixture cuts at 3.3 and 7.7, and `33/10` is not a binary
//! fraction) — and the cache agrees on **all 48**. So being unrepresentable does not make a
//! coordinate disagree.
//!
//! What actually separates the families: on the narrow ones both roads round *the same rational*
//! the same way, so they cannot differ. On the tilted ones the cached coordinate comes out of a
//! longer `f64` derivation that is not a correct rounding of anything — the 59-bit width is a
//! proxy for how much arithmetic happened, not the mechanism. This is the first measurement of
//! `truth-and-cache.md`'s *"최근접 f64 가 아닐 수 있다"* on vertices either way.

use nacre_geom::Surface;
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolKind, DatumDef, Edge2d, OpOutput, Operation, Precision, Profile2d, RealizeError,
    SketchFrame, SketchPlane, apply, boolean, from_edges, realize_vertex, realize_vertex_decimal,
};
use nacre_scalar::Axis;
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

fn cuboid(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
    m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
}

fn boolean_corner() -> Model {
    let mut m = Model::new();
    let a = cuboid(&mut m, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let b = cuboid(&mut m, [3.3, 3.3, -1.0], [7.7, 7.7, 11.0]);
    boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    m.rebuild_adjacency();
    m
}

fn boolean_twice() -> Model {
    let mut m = boolean_corner();
    let a = m.live_solids[0];
    let b = cuboid(&mut m, [-1.0, 4.4, 4.4], [11.0, 8.8, 8.8]);
    boolean(&mut m, BoolKind::Cut, a, b).expect("second cut");
    m.rebuild_adjacency();
    m
}

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
// (ii) — is the f64 road a different plane?
// ---------------------------------------------------------------------------------------------

/// A plain cylinder — a circle extruded. Its rim corners are `VertexDef::OnSeam`.
fn cylinder() -> Model {
    let mut m = Model::new();
    let profile = from_edges(vec![Edge2d::circle(p2(2.0, 2.0), 3.0).unwrap()])
        .expect("a circle is a profile")
        .remove(0);
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 4.0,
        },
    )
    .expect("the circle extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    m
}

/// A box with a bore through it — the wall/bore crossings are `VertexDef::Branch`.
fn bored_plate() -> Model {
    let mut m = Model::new();
    let profile = from_edges(vec![
        Edge2d::line(p2(-1.0, -1.0), p2(5.0, -1.0)).unwrap(),
        Edge2d::line(p2(5.0, -1.0), p2(5.0, 5.0)).unwrap(),
        Edge2d::line(p2(5.0, 5.0), p2(-1.0, 5.0)).unwrap(),
        Edge2d::line(p2(-1.0, 5.0), p2(-1.0, -1.0)).unwrap(),
        Edge2d::circle(p2(2.0, 2.0), 1.0).unwrap(),
    ])
    .expect("a plate with a bore")
    .remove(0);
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 3.0,
        },
    )
    .expect("the bored plate extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    m
}

/// **A hand-built `Branch` vertex** — the kernel does not mint these yet.
///
/// `VertexDef::Branch`'s own doc records why: *"The producer arrives with M6-2's boolean; until
/// then hand-built fixtures and validate are the consumers"* — the assembler declines to mint a
/// branch node because class order and handle order are canonical in different index spaces. So a
/// fixture that waited for a boolean to produce one would measure nothing, forever.
///
/// The shape: a radius-3 cylinder on the world Z axis at `(2, 2)`, and a cuboid whose faces supply
/// plane carriers. Every plane pair is offered to the door; the ones whose meet line actually
/// crosses the lateral surface answer, the rest refuse. Both outcomes are asserted, so this cannot
/// pass by refusing everything.
fn branch_vertices(m: &mut Model) -> Vec<Handle<Vertex>> {
    use nacre_geom::Surface;
    use nacre_topo::{QuadRoot, VertexDef};
    let mut cyl = None;
    let mut planes = Vec::new();
    for i in 0..m.surface_count() as u32 {
        let Some(h) = m.surface_handle_at(i) else {
            continue;
        };
        match m.surface(h) {
            Surface::Cylinder(_) => cyl = Some(h),
            Surface::Plane(_) => planes.push(h),
        }
    }
    let cylinder = cyl.expect("the fixture has a cylinder");
    let mut out = Vec::new();
    for a in 0..planes.len() {
        for b in (a + 1)..planes.len() {
            // Ascending handle order is `Branch`'s stored convention.
            let (p0, p1) = (planes[a], planes[b]);
            for root in [QuadRoot::Lo, QuadRoot::Hi, QuadRoot::Double] {
                out.push(m.push_vertex(
                    VertexDef::Branch {
                        planes: [p0, p1],
                        cylinder,
                        root,
                    },
                    Point3::from_array([0.0; 3]),
                    Some(0.0),
                ));
            }
        }
    }
    out
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

/// ★★★ **Lock 0 — the negative control, and the one that cannot be argued with.**
///
/// An axis-aligned box's corners are 7-bit rationals, so the truth *is* an `f64` and any correct
/// realization must reproduce the cache bit for bit. A door that quietly re-derived a slightly
/// different point would show up here and nowhere else.
#[test]
fn a_box_realizes_to_exactly_what_the_cache_holds() {
    let m = boolean_corner();
    let mut n = 0;
    for vh in live_vertices(&m) {
        let r = realize_vertex(&m, vh, Precision::NearestF64).expect("a box corner realizes");
        assert!(
            r.is_exact(),
            "a three-plane meet with no motion is a rational"
        );
        let (v, e) = r.to_f64().expect("an exact value names an f64");
        let cached = m.vertex_point(vh).as_array();
        assert_eq!(v, cached, "vertex {vh:?}");
        // ⚠ **Width is not representability.** These corners are 7 bits wide by
        // `tests/point_width.rs`'s metric, and `boolean_corner` cuts at 3.3 and 7.7 — `33/10` is
        // an exact `Rat` and no `f64` at all. So the readout rounds, and the error says so; a
        // coordinate that *is* dyadic (0, 10, …) reports zero. Asserting `[0.0; 3]` here passed
        // only because the arm used to claim it unconditionally.
        for k in 0..3 {
            let dyadic = format!("{:.60}", v[k]) == r.to_decimal(60).expect("exact")[k];
            assert_eq!(e[k] == 0.0, dyadic, "error {} vs dyadic {dyadic}", e[k]);
        }
        n += 1;
    }
    assert!(
        n >= 8,
        "swept {n} vertices — the fixture never reached the target"
    );
}

/// ★★★ **Lock 1 — idempotent and unique.** Realizing at more bits cannot change a decided answer,
/// because correct rounding makes the realization unique.
#[test]
fn more_bits_do_not_change_a_decided_answer() {
    let mut compared = 0;
    for m in [boolean_corner(), tilted_frame(1), boolean_twice()] {
        for vh in live_vertices(&m) {
            let (Ok(a), Ok(b)) = (
                realize_vertex(&m, vh, Precision::Bits(128)),
                realize_vertex(&m, vh, Precision::Bits(1024)),
            ) else {
                continue;
            };
            if let (Some((va, _)), Some((vb, _))) = (a.to_f64(), b.to_f64()) {
                assert_eq!(va, vb, "f64 moved with precision at {vh:?}");
                compared += 1;
            }
            for places in [5usize, 20, 40] {
                if let (Some(da), Some(db)) = (a.to_decimal(places), b.to_decimal(places)) {
                    assert_eq!(da, db, "digits moved with precision at {vh:?}");
                }
            }
        }
    }
    assert!(
        compared > 20,
        "only {compared} comparisons — instrument too small"
    );
}

/// ★★ **Lock 1's negative control.** A realization that cannot decide a digit must say so rather
/// than invent it. Without this, "always `Some`" passes every other test in this file.
#[test]
fn too_few_bits_report_undecided_rather_than_guessing() {
    let m = tilted_frame(1);
    let mut undecided = 0;
    let mut decided = 0;
    for vh in live_vertices(&m) {
        let Ok(r) = realize_vertex(&m, vh, Precision::Bits(128)) else {
            continue;
        };
        if r.is_exact() {
            continue; // an exact ratio decides every digit; nothing to withhold
        }
        match r.to_decimal(200) {
            None => undecided += 1,
            Some(_) => decided += 1,
        }
    }
    assert!(
        undecided > 0,
        "128 bits decided 200 places on every realized vertex ({decided}) — the guard is not \
         being exercised, so this file cannot see it go missing"
    );
}

/// ★★ **Lock 1′ — the door's own contract: ask for `places`, get `places`.**
///
/// The `None` on `to_decimal` is a rung-to-rung signal. A caller must never see it, and an
/// implementation that always reported undecided would satisfy every other lock but this one.
///
/// ★★★ **Counted per fixture, and refusals are counted too.** An earlier spelling wrote
/// `let Ok(d) = … else { continue }` and asserted a total — which drops exactly the vertices a
/// broken ladder would fail on, and stays green while the population shrinks under it. That is the
/// failure `tests/point_width.rs` records in its own module doc; this shape cannot have it.
#[test]
fn the_door_returns_the_places_it_was_asked_for() {
    // 120 places is past what the first rung decides (128 bits is ~38 decimal digits), so the
    // tilted family can only answer it by climbing — and a realization rounded at the wrong
    // precision shows up here as a refusal.
    for (what, m, floor) in [
        ("boolean_corner", boolean_corner(), 8usize),
        ("tilted_frame", tilted_frame(1), 8),
    ] {
        for places in [50usize, 120] {
            let (mut answered, mut refused) = (0usize, 0usize);
            for vh in live_vertices(&m) {
                match realize_vertex_decimal(&m, vh, places) {
                    Ok(d) => {
                        for s in &d {
                            let frac = s.split_once('.').expect("a decimal point").1;
                            assert_eq!(
                                frac.len(),
                                places,
                                "{what}: asked {places} places, got {}: {s}",
                                frac.len()
                            );
                        }
                        answered += 1;
                    }
                    Err(_) => refused += 1,
                }
            }
            assert!(
                answered >= floor,
                "{what} at {places} places: only {answered} answered ({refused} refused) — \
                 the ladder is not reaching what it reached before"
            );
        }
    }
}

/// ★★ **Lock 2 — every variant answers, and it is raised for real rather than tabulated.**
///
/// A table of variants proves nothing about a dispatch; these fixtures actually mint an `OnSeam`
/// rim corner and a `Branch` wall crossing, and the test fails if a population goes empty — which
/// is how deleting one arm shows up as something other than a smaller number nobody reads.
#[test]
fn every_vertex_variant_answers() {
    use nacre_topo::VertexDef;
    let (mut seam, mut branch, mut three) = (0usize, 0usize, 0usize);
    for m in [cylinder(), bored_plate(), boolean_corner()] {
        for vh in live_vertices(&m) {
            let kind = m.vertices.get(vh).def;
            let Ok(d) = realize_vertex_decimal(&m, vh, 25) else {
                continue;
            };
            for s in &d {
                assert_eq!(s.split_once('.').expect("a point").1.len(), 25);
            }
            match kind {
                VertexDef::OnSeam(_) => seam += 1,
                VertexDef::Branch { .. } => branch += 1,
                VertexDef::ThreePlane(_) => three += 1,
            }
        }
    }
    // `Branch` has no producer in the kernel yet, so its population is built here.
    let mut m = cylinder();
    let _ = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let mut refused = 0usize;
    for vh in branch_vertices(&mut m) {
        match realize_vertex_decimal(&m, vh, 25) {
            Ok(d) => {
                for s in &d {
                    assert_eq!(s.split_once('.').expect("a point").1.len(), 25);
                }
                // ★ An oracle that owes nothing to the derivation: a branch point is a point
                // **on the cylinder**, so its distance from the axis is the radius. The fixture's
                // bore is radius 3 about `(2, 2)` on the world Z axis.
                let (x, y) = (
                    d[0].parse::<f64>().expect("a decimal"),
                    d[1].parse::<f64>().expect("a decimal"),
                );
                let rr = ((x - 2.0).powi(2) + (y - 2.0).powi(2)).sqrt();
                assert!(
                    (rr - 3.0).abs() < 1e-9,
                    "a branch point is {rr} from the axis, not the radius 3"
                );
                branch += 1;
            }
            // ★ **Which refusal, not just that there was one.** A plane pair whose meet line
            // misses the cylinder is a *curved* definition failing to resolve; `vertex_meet` is
            // never called on this road, so answering `NoMeet` would name a decline that did not
            // happen. Planted: swapping the two left every other test here green.
            Err(e) => {
                assert_eq!(
                    e,
                    RealizeError::NoCurvedPoint,
                    "a branch pair that does not cross should say so in its own words"
                );
                refused += 1;
            }
        }
    }
    println!("answered: three_plane={three} on_seam={seam} branch={branch} (refused {refused})");
    assert!(three > 0, "no three-plane vertex answered");
    assert!(
        seam > 0,
        "no OnSeam vertex answered — that arm is not being exercised"
    );
    assert!(
        branch > 0,
        "no Branch vertex answered — that arm is not being exercised"
    );
    assert!(
        refused > 0,
        "every offered plane pair crossed the cylinder — the refusal side is untested"
    );
}

/// ★★★ **The curved arms are checked against a route that shares nothing with them.**
///
/// A seam or branch coordinate is assembled here out of `inv_sqrt_bounded` and exact rational
/// arithmetic; the cached one was computed at boolean time by different code entirely. Agreement
/// to the cache's own accuracy is therefore evidence, where comparing two spellings of the same
/// derivation would not be (`tests/` has that failure on record).
#[test]
fn a_curved_vertex_agrees_with_the_cache_it_did_not_use() {
    use nacre_topo::VertexDef;
    let mut checked = 0usize;
    let mut worst = 0.0f64;
    for m in [cylinder(), bored_plate()] {
        for vh in live_vertices(&m) {
            if matches!(m.vertices.get(vh).def, VertexDef::ThreePlane(_)) {
                continue;
            }
            let Ok(r) = realize_vertex(&m, vh, Precision::Bits(256)) else {
                continue;
            };
            let Some((v, _)) = r.to_f64() else { continue };
            let cached = m.vertex_point(vh).as_array();
            for k in 0..3 {
                worst = worst.max((v[k] - cached[k]).abs());
            }
            checked += 1;
        }
    }
    println!("curved vertices checked={checked} worst deviation={worst:e}");
    assert!(checked >= 2, "only {checked} curved vertices checked");
    assert!(
        worst < 1e-9,
        "a realized curved coordinate is {worst:e} from the cache — one of the two is wrong"
    );
}

/// ★★★ **What the cache costs — the measurement this door was built to make.**
///
/// Not `#[ignore]`d: it carries an assertion about the *narrow* populations, which is the half
/// that must never regress. The tilted row is printed rather than pinned, because it is a
/// statement about today's cache and the point of the door is to make it false eventually.
#[test]
fn the_cache_is_not_always_nearest() {
    for (what, m) in [
        ("axis-aligned box", boolean_corner()),
        ("boolean'd twice", boolean_twice()),
        ("tilted frame", tilted_frame(1)),
    ] {
        let (mut compared, mut nearest, mut worst_ulp) = (0usize, 0usize, 0i64);
        for vh in live_vertices(&m) {
            let Ok(r) = realize_vertex(&m, vh, Precision::NearestF64) else {
                continue;
            };
            let Some((v, _)) = r.to_f64() else { continue };
            compared += 1;
            let cached = m.vertex_point(vh).as_array();
            let mut ok = true;
            for k in 0..3 {
                if v[k] != cached[k] {
                    ok = false;
                    let d = (v[k].to_bits() as i64 - cached[k].to_bits() as i64).abs();
                    worst_ulp = worst_ulp.max(d);
                }
            }
            nearest += usize::from(ok);
        }
        println!("{what:<18} compared={compared:<3} nearest={nearest:<3} worst_ulp={worst_ulp}");
        if what != "tilted frame" {
            assert_eq!(
                nearest, compared,
                "{what}: a 7-bit coordinate is an f64 exactly, so the cache cannot differ"
            );
        }
    }
}

/// ★★★ **A realized coordinate never reports the error an exact one does.**
///
/// `to_f64` hands back `(value, error)`, and `0.0` there means *exact* — the rational arm's
/// promise. An approached coordinate must not be able to say it: at 4096 bits the radius is near
/// `2⁻⁴⁰⁹¹`, and reading it out through `2f64.powi` flushed every coordinate to `0e0`. Measured,
/// shipped, and caught only by this file's audit — so the lock lives here now.
#[test]
fn an_approached_coordinate_never_claims_to_be_exact() {
    let m = tilted_frame(1);
    let (mut approached, mut exact) = (0usize, 0usize);
    for vh in live_vertices(&m) {
        for bits in [128usize, 512, 1024, 2048, 4096] {
            let Ok(r) = realize_vertex(&m, vh, Precision::Bits(bits)) else {
                continue;
            };
            let Some((v, e)) = r.to_f64() else { continue };
            if r.is_exact() {
                // ⚠ **Not `[0.0; 3]`.** An exact *realization* still rounds when it is read out
                // at 53 bits, and the tilted family's coordinates are 59 bits wide. Zero is
                // allowed only where the rational really is an `f64`; where it is not, the error
                // must be positive and no wider than half an ulp. An earlier spelling of this
                // lock asserted zero and so pinned the lie it was written to catch.
                let d = r
                    .to_decimal(60)
                    .expect("an exact realization prints every place");
                for k in 0..3 {
                    let readout_is_the_value = format!("{:.60}", v[k]) == d[k];
                    match readout_is_the_value {
                        true => assert_eq!(e[k], 0.0, "a representable rational has no error"),
                        false => {
                            assert!(e[k] > 0.0, "a rounded readout reported no error: {}", d[k]);
                            assert!(
                                e[k] <= (v[k].abs() * f64::EPSILON).max(f64::MIN_POSITIVE),
                                "error {} is wider than an ulp of {}",
                                e[k],
                                v[k]
                            );
                        }
                    }
                }
                exact += 1;
                continue;
            }
            for c in e {
                assert!(
                    c > 0.0,
                    "a realization at {bits} bits reported error {c:e} — that is the exact \
                     arm's answer, and this coordinate is not exact"
                );
            }
            approached += 1;
        }
    }
    // Both arms have to be present, or this measures one of them and calls it both.
    assert!(approached > 8, "only {approached} approached readings");
    assert!(
        exact > 0,
        "no exact reading — the negative control is missing"
    );
}

/// ★★ **Refusals are named, and never fall back to the cache.**
///
/// A "there is none" lock needs its positive twin in the same test, or a door that refused
/// everything would pass it.
#[test]
fn a_refusal_is_named_and_a_success_stands_beside_it() {
    let m = tilted_frame(1);
    let (mut ok, mut refused) = (0usize, 0usize);
    for vh in live_vertices(&m) {
        match realize_vertex(&m, vh, Precision::NearestF64) {
            Ok(_) => ok += 1,
            Err(e) => {
                assert!(
                    matches!(
                        e,
                        RealizeError::NoMeet
                            | RealizeError::WideUnderMotion
                            | RealizeError::NoMotionChain
                            | RealizeError::NoCurvedPoint
                            | RealizeError::Undecided
                    ),
                    "unnamed refusal {e:?}"
                );
                refused += 1;
            }
        }
    }
    assert!(
        ok > 0,
        "nothing realized — the fixture is not exercising success"
    );
    assert!(
        refused > 0,
        "nothing refused on the tilted family — this lock is measuring nothing"
    );
}

/// ★ A rational vertex prints digits that are **the coordinate's own**, not an approximation's.
#[test]
fn an_exact_vertex_prints_exact_digits() {
    let m = boolean_corner();
    let vh = live_vertices(&m)[0];
    let d = realize_vertex_decimal(&m, vh, 30).expect("a box corner");
    for s in &d {
        let frac = s.split_once('.').expect("a decimal point").1;
        assert_eq!(frac.len(), 30);
        // A 7-bit rational is exactly representable in decimal, so past the significant digits
        // every place is a zero — an approximation would show noise instead.
        assert!(
            frac.trim_end_matches('0').len() <= 10,
            "an exact coordinate should not have 30 places of content: {s}"
        );
    }
}

/// ★★★★ **The deepest claim this door makes: a decided digit is a true digit.**
///
/// Everything else rests on the error radius really bounding the error — under-report it and
/// `round_to_digits` accepts a place the definition never determined, which is a *wrong digit
/// printed confidently*, the worst thing this feature could do. The radius comes from arithmetic
/// written for this cell (`mul_bounded`/`add_bounded`/`rat_bounded`), so it is not something to
/// take on faith.
///
/// The oracle is the same door at 8192 bits. Compared in **decimal**: an `f64` comparison
/// collapses every difference a high-precision radius is about, and an earlier spelling of this
/// probe reported "0 violations" while measuring nothing.
///
/// ★★ **Calibrated, because "0 wrong" is only worth what the instrument can see.** Shrinking the
/// radius by 2⁻²⁰ already produces 12 wrong digits here; 2⁻⁶⁰ gives 48, 2⁻¹²⁰ gives 192. So the
/// bound is honest to within about 20 bits of the true error — tight, not vacuously wide. (Three
/// smaller mutations — dropping the product's rounding, the operand cross-terms, or the rational
/// conversion's own half-ulp — leave it green, which says those terms are slack inside that
/// margin, not that the test is asleep.)
#[test]
fn a_digit_this_door_decides_is_the_coordinates_own() {
    let (mut checked, mut wrong) = (0usize, 0usize);
    for m in [cylinder(), tilted_frame(1)] {
        for vh in live_vertices(&m) {
            let Ok(truth) = realize_vertex(&m, vh, Precision::Bits(8192)) else {
                continue;
            };
            if truth.is_exact() {
                continue; // an exact ratio has no radius to be wrong about
            }
            for bits in [64usize, 80, 96, 128, 192, 256, 512] {
                let Ok(r) = realize_vertex(&m, vh, Precision::Bits(bits)) else {
                    continue;
                };
                for places in [5usize, 15, 25, 40, 70, 120] {
                    let (Some(got), Some(want)) = (r.to_decimal(places), truth.to_decimal(places))
                    else {
                        continue;
                    };
                    for k in 0..3 {
                        checked += 1;
                        wrong += usize::from(got[k] != want[k]);
                    }
                }
            }
        }
    }
    assert!(
        checked > 300,
        "only {checked} comparisons — instrument too small"
    );
    assert_eq!(
        wrong, 0,
        "{wrong} of {checked} decided digits were not the coordinate's"
    );
}
