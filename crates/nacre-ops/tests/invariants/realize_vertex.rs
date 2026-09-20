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
//! proxy for how much arithmetic happened, not the mechanism. This measures, on vertices, that
//! the cache may not be the nearest f64.

use crate::stated::*;
use nacre_math::{Point3, Vector3};
use nacre_ops::{
    BoolKind, CacheDecline, OpOutput, Operation, Precision, Profile2d, RealizeError, SketchFrame,
    SketchPlane, apply, boolean, realize_cache, realize_vertex, realize_vertex_decimal,
};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, Solid, Surface, Vertex};

use crate::fixtures::{datum_frame, p2};

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
    let a = m.live_solids()[0];
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
                        .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
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

/// A plain cylinder — a circle extruded. Its rim corners are `Vertex::OnSeam`.
fn cylinder() -> Model {
    let mut m = Model::new();
    let profile = stated(vec![circle(p2(2.0, 2.0), 3.0)])
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

/// A box with a bore through it — the wall/bore crossings are `Vertex::Pierce`.
fn bored_plate() -> Model {
    let mut m = Model::new();
    let profile = stated(vec![
        line(p2(-1.0, -1.0), p2(5.0, -1.0)),
        line(p2(5.0, -1.0), p2(5.0, 5.0)),
        line(p2(5.0, 5.0), p2(-1.0, 5.0)),
        line(p2(-1.0, 5.0), p2(-1.0, -1.0)),
        circle(p2(2.0, 2.0), 1.0),
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

/// **A hand-built `Pierce` vertex** — the kernel does not mint these yet.
///
/// Hand-built fixtures and validate are the consumers: the assembler declines to mint a
/// pierce node because class order and handle order are canonical in different index spaces. So a
/// fixture that waited for a boolean to produce one would measure nothing, forever.
///
/// The shape: a radius-3 cylinder on the world Z axis at `(2, 2)`, and a cuboid whose faces supply
/// plane carriers. Every plane pair is offered to the door; the ones whose meet line actually
/// crosses the lateral surface answer, the rest refuse. Both outcomes are asserted, so this cannot
/// pass by refusing everything.
fn pierce_vertices(m: &mut Model) -> Vec<Handle<Vertex>> {
    use nacre_topo::{QuadRoot, Vertex};
    let mut cyl = None;
    let mut planes = Vec::new();
    for i in 0..m.surface_count() as u32 {
        let Some(h) = m.surface_handle_at(i) else {
            continue;
        };
        match m.surface_cache(h) {
            nacre_geom::Surface::Cylinder(_) => cyl = Some(h),
            nacre_geom::Surface::Plane(_) => planes.push(h),
        }
    }
    let cylinder = cyl.expect("the fixture has a cylinder");
    let mut out = Vec::new();
    for a in 0..planes.len() {
        for b in (a + 1)..planes.len() {
            // Ascending handle order is `Pierce`'s stored convention.
            let (p0, p1) = (planes[a], planes[b]);
            for root in [QuadRoot::Lo, QuadRoot::Hi, QuadRoot::Double] {
                out.push(m.push_vertex(
                    Vertex::Pierce {
                        planes: [p0, p1],
                        cylinder,
                        root,
                    },
                    PointCache::Unrealized {
                        coord: Point3::from_array([0.0; 3]),
                    },
                ));
            }
        }
    }
    out
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
            assert_eq!(
                e[k].is_zero(),
                dyadic,
                "error {:?} vs dyadic {dyadic}",
                e[k]
            );
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
/// rim corner and a `Pierce` wall crossing, and the test fails if a population goes empty — which
/// is how deleting one arm shows up as something other than a smaller number nobody reads.
#[test]
fn every_vertex_variant_answers() {
    use nacre_topo::Vertex;
    let (mut seam, mut pierce, mut three) = (0usize, 0usize, 0usize);
    for m in [cylinder(), bored_plate(), boolean_corner()] {
        for vh in live_vertices(&m) {
            let kind = *m.vertex(vh);
            let Ok(d) = realize_vertex_decimal(&m, vh, 25) else {
                continue;
            };
            for s in &d {
                assert_eq!(s.split_once('.').expect("a point").1.len(), 25);
            }
            match kind {
                Vertex::OnSeam(_) => seam += 1,
                Vertex::Pierce { .. } => pierce += 1,
                Vertex::ThreePlane(_) => three += 1,
            }
        }
    }
    // `Pierce` has no producer in the kernel yet, so its population is built here.
    let mut m = cylinder();
    let _ = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let mut refused = 0usize;
    for vh in pierce_vertices(&mut m) {
        match realize_vertex_decimal(&m, vh, 25) {
            Ok(d) => {
                for s in &d {
                    assert_eq!(s.split_once('.').expect("a point").1.len(), 25);
                }
                // ★ An oracle that owes nothing to the derivation: a pierce point is a point
                // **on the cylinder**, so its distance from the axis is the radius. The fixture's
                // bore is radius 3 about `(2, 2)` on the world Z axis.
                let (x, y) = (
                    d[0].parse::<f64>().expect("a decimal"),
                    d[1].parse::<f64>().expect("a decimal"),
                );
                let rr = ((x - 2.0).powi(2) + (y - 2.0).powi(2)).sqrt();
                assert!(
                    (rr - 3.0).abs() < 1e-9,
                    "a pierce point is {rr} from the axis, not the radius 3"
                );
                pierce += 1;
            }
            // ★ **Which refusal, not just that there was one.** A plane pair whose meet line
            // misses the cylinder is a *curved* definition failing to resolve; `vertex_meet` is
            // never called on this road, so answering `NoMeet` would name a decline that did not
            // happen. Planted: swapping the two left every other test here green.
            Err(e) => {
                assert_eq!(
                    e,
                    RealizeError::NoCurvedPoint,
                    "a pierce pair that does not cross should say so in its own words"
                );
                refused += 1;
            }
        }
    }
    println!("answered: three_plane={three} on_seam={seam} pierce={pierce} (refused {refused})");
    assert!(three > 0, "no three-plane vertex answered");
    assert!(
        seam > 0,
        "no OnSeam vertex answered — that arm is not being exercised"
    );
    assert!(
        pierce > 0,
        "no Pierce vertex answered — that arm is not being exercised"
    );
    assert!(
        refused > 0,
        "every offered plane pair crossed the cylinder — the refusal side is untested"
    );
}

/// ★★★ **The curved arms are checked against a route that shares nothing with them.**
///
/// A seam or pierce coordinate is assembled here out of `inv_sqrt_bounded` and exact rational
/// arithmetic; the cached one was computed at boolean time by different code entirely. Agreement
/// to the cache's own accuracy is therefore evidence, where comparing two spellings of the same
/// derivation would not be (`tests/` has that failure on record).
#[test]
fn a_curved_vertex_agrees_with_the_cache_it_did_not_use() {
    use nacre_topo::Vertex;
    let mut checked = 0usize;
    let mut worst = 0.0f64;
    for m in [cylinder(), bored_plate()] {
        for vh in live_vertices(&m) {
            if matches!(*m.vertex(vh), Vertex::ThreePlane(_)) {
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

/// ★★★ **An operation's cache is the realization.** Every vertex the push funnel could realize
/// on the ladder's first rung carries `Bounded` with that realization bit for bit; every one it
/// could not is not `Bounded`.
#[test]
fn an_operations_cache_is_the_realization() {
    for (what, m) in [
        ("axis-aligned box", boolean_corner()),
        ("boolean'd twice", boolean_twice()),
        ("tilted frame", tilted_frame(1)),
    ] {
        let (mut realized, mut kept) = (0usize, 0usize);
        for vh in live_vertices(&m) {
            match realize_cache(&m, m.vertex(vh)) {
                Ok((v, bound)) => {
                    realized += 1;
                    let PointCache::Bounded { coord, bound: b } = *m.vertex_cache(vh) else {
                        panic!(
                            "{what}: a realizable vertex is not Bounded: {:?}",
                            m.vertex_cache(vh)
                        );
                    };
                    assert_eq!(
                        coord.as_array(),
                        v,
                        "{what}: the cache is the realization, bit for bit"
                    );
                    assert_eq!(b, bound, "{what}: and carries its bound");
                }
                Err(_) => {
                    kept += 1;
                    assert!(
                        !matches!(m.vertex_cache(vh), PointCache::Bounded { .. }),
                        "{what}: a refused vertex cannot be Bounded"
                    );
                }
            }
        }
        println!("{what:<18} realized={realized} kept={kept}");
        assert!(realized > 0, "{what}: nothing realized — the road is dead");
    }
}

/// ★ **The edge cache is already the derivation of the realized endpoints.** An edge is pushed
/// after its vertices, so `push_edge` derives its line from realized coordinates; rebuilding every
/// edge curve afterwards changes nothing — the derivation, standing on realized endpoints.
///
/// ⚠ This is a statement about a model **nothing has refined**. `refine_vertex_cache` moves
/// coordinates, and after it a rebuild is emphatically not a no-op — which is why that door
/// re-derives the curves itself, and why `the_refine_door_carries_the_edges_with_it` asserts the
/// rebuild is a no-op only *after* the door has already run.
#[test]
fn the_edge_cache_is_the_derivation_of_realized_endpoints() {
    let mut m = tilted_frame(1);
    let before: Vec<_> = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .map(|h| (h, m.edge(h)))
        .map(|(h, _)| m.edge_curve(h).clone())
        .collect();
    m.rebuild_edge_cache();
    let after: Vec<_> = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .map(|h| (h, m.edge(h)))
        .map(|(h, _)| m.edge_curve(h).clone())
        .collect();
    assert_eq!(before, after);
}

/// ★ **A moved solid is realized from its moved definition** — not from the moved `f64`. A
/// rational translation of the tilted prism: every moved vertex the road can realize is `Bounded`
/// and bit for bit the realization; the rest keep the construction's figure, counted.
#[test]
fn a_moved_solid_is_realized_from_its_moved_definition() {
    let mut m = tilted_frame(1);
    let solid = m.live_solids()[0];
    let OpOutput::Transform { solid: moved } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: nacre_scalar::Isometry::translation([
                nacre_scalar::Rat::new(7, 11).unwrap(),
                nacre_scalar::Rat::from_int(-2),
                nacre_scalar::Rat::new(1, 4).unwrap(),
            ]),
        },
    )
    .expect("a rational translation moves the prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert_eq!(m.live_solids().to_vec(), vec![moved]);
    let (mut bounded, mut kept) = (0usize, 0usize);
    for vh in live_vertices(&m) {
        match *m.vertex_cache(vh) {
            PointCache::Bounded { coord, .. } => {
                bounded += 1;
                let (v, _) = realize_vertex(&m, vh, Precision::NearestF64)
                    .expect("a Bounded vertex realizes")
                    .to_f64()
                    .expect("decided");
                assert_eq!(coord.as_array(), v);
            }
            // ★ Named, not a catch-all: `Ceiling` and `Unrealized` are different reports, and a
            // `_` here would have counted a cost-capped vertex as a refused one.
            PointCache::Ceiling { .. } | PointCache::Unrealized { .. } => kept += 1,
        }
    }
    println!("moved tilted prism: bounded={bounded} kept={kept}");
    assert!(
        bounded > 0,
        "nothing realized after the move — the road is dead"
    );
}

/// ★★★ **A realized coordinate never reports the error an exact one does.**
///
/// `to_f64` hands back `(value, error)`, and `0.0` there means *exact* — the rational arm's
/// promise. An approached coordinate must not be able to say it: at 4096 bits the radius is near
/// `2⁻⁴⁰⁹¹`, and reading it out through `2f64.powi` flushed every coordinate to `0e0`. Measured,
/// shipped, and caught only by this file's audit — so the lock lives here now.
///
/// One exception, earned rather than claimed: a coordinate that is exactly zero from
/// exact inputs carries a zero radius honestly — `Mag::above`'s zero guard leaves nothing to
/// charge — and zero is the only value this arm can earn that way, so a zero radius must sit on
/// a zero value. The cylinder's axis seam is the fixture that has one; the tilted frame has none.
#[test]
fn an_approached_coordinate_never_claims_to_be_exact() {
    let (mut approached, mut exact, mut zero_earned) = (0usize, 0usize, 0usize);
    for m in [cylinder(), tilted_frame(1)] {
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
                            true => {
                                assert!(e[k].is_zero(), "a representable rational has no error")
                            }
                            false => {
                                assert!(
                                    !e[k].is_zero(),
                                    "a rounded readout reported no error: {}",
                                    d[k]
                                );
                                assert!(
                                    e[k].lt(nacre_scalar::Mag::of(v[k])
                                        .times(nacre_scalar::Mag::pow2(-51))),
                                    "error {:?} is wider than an ulp of {}",
                                    e[k],
                                    v[k]
                                );
                            }
                        }
                    }
                    exact += 1;
                    continue;
                }
                // ★ A zero radius is allowed on this arm only where it is earned, and the machine
                // can earn exactly one value: zero. Every `add`/`mul` charges `round_off(result)`,
                // which `Mag::above`'s zero guard makes zero only for a zero result — so a seam
                // point's axis coordinate (exact centre 0, exact direction 0) is exact and says
                // so, while a flushed radius on any nonzero coordinate still fails here.
                for k in 0..3 {
                    if e[k].is_zero() {
                        assert_eq!(
                            v[k], 0.0,
                            "a realization at {bits} bits reported no error on a nonzero \
                             coordinate — the exact arm's answer, unearned"
                        );
                        zero_earned += 1;
                    }
                }
                approached += 1;
            }
        }
    }
    // Both arms have to be present, or this measures one of them and calls it both — and the
    // earned zero has to occur, or the rule above is never exercised.
    assert!(approached > 8, "only {approached} approached readings");
    assert!(
        zero_earned > 0,
        "no earned zero radius — the cylinder's axis seam is the positive control"
    );
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
                // ★★★ **Exhaustive `match`, not `matches!`** — and the difference is not style.
                // This site exists to notice when a refusal arrives without a name, so it has to
                // be the thing that breaks when a variant is added. `matches!(e, A | B | C)`
                // compiles happily forever: a new variant (`Unrepresentable`, say) leaves the list
                // stale, asserting "unnamed refusal" about a refusal that has a name. An
                // exhaustive match makes the compiler point here instead.
                match e {
                    RealizeError::NoMeet
                    | RealizeError::WideUnderMotion
                    | RealizeError::NoMotionChain
                    | RealizeError::NoCurvedPoint
                    | RealizeError::Undecided
                    | RealizeError::Unrepresentable => {}
                }
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
    eprintln!("[decided digits] {checked} checked, {wrong} wrong");
    assert!(
        checked > 300,
        "only {checked} comparisons — instrument too small"
    );
    assert_eq!(
        wrong, 0,
        "{wrong} of {checked} decided digits were not the coordinate's"
    );
}

fn moved(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    match apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    ) {
        Ok(OpOutput::Transform { solid }) => solid,
        other => panic!("transform: {other:?}"),
    }
}

fn turn(axis: Axis, deg: i128) -> Isometry {
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).expect("a whole-degree angle"),
    })
}

/// ★★★★ **A motion chain is held back for the bits it costs, never for its length** — the two
/// shapes that prove the difference.
///
/// The cache road was guarded with a depth constant (`CACHE_REPLAY_DEPTH = 64`) because a
/// 4,200-turn history realized on every push turned a 5.6 s test into minutes. That constant was
/// doing **two jobs** and only one of them was true. As a *cost* limit it was load-bearing —
/// removing it took the same fixture to 30 s. As a *precision* rule it was measurably wrong: what
/// costs precision is an **irrational turn**, and a rational translation costs nothing at all
/// (`Angle::try_exact_cos_sin` answers only 0/90/180/270 — Niven, so those add no radius). So the
/// guard held back chains that decide perfectly: rational translations by the hundred, and every
/// 7° turn between the 65th and the 90th.
///
/// The limit is now named for the job it does (`CACHE_REPLAY_COST_CAP`, sized from measured cost
/// against a stated per-vertex budget), and inside it the ladder's first rung decides. These two
/// come back `Bounded`.
///
/// ★ Each row also asserts that some vertex **really replays a chain** (`!is_exact`): a fixture
/// whose motion was folded away would pass this while measuring nothing, and that guard *bit twice
/// while this was written*.
///
/// ⚠ **A pure quadrantal chain is not a third row, and finding out why corrected the claim.** A
/// right-angle turn about the origin is exact, so `transform` carries it whole and restates the
/// planes rather than recording a node: a hundred of them leave *no chain at all*, and the guard
/// caught the empty measurement. Seeding one 37° turn to force a history then failed the other way
/// — at 128 bits that coordinate is undecidable, full walk and all. The lesson is that "a
/// quadrantal turn costs nothing" is true of the **angle** (`cos`/`sin` are exact, contributing no
/// radius) and false of the **node**: interval arithmetic still charges its roundings, so a radius
/// already in flight keeps compounding. Zero times anything is still zero, which is why a chain
/// with an exact base decides forever — but such a chain is never recorded, so the depth guard
/// never saw it either.
#[test]
fn a_chain_is_held_back_for_the_bits_it_costs_not_for_its_length() {
    for (what, n, rotate) in [
        ("150 rational translations", 150usize, false),
        ("80 turns of 7°", 80, true),
    ] {
        let mut m = Model::new();
        let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
        m.rebuild_adjacency();
        for _ in 0..n {
            let iso = if rotate {
                turn(Axis::Z, 7)
            } else {
                Isometry::translation([
                    Rat::new(1, 7).expect("1/7"),
                    Rat::from_int(0),
                    Rat::from_int(0),
                ])
            };
            s = moved(&mut m, s, iso);
        }
        m.rebuild_adjacency();

        let vs = live_vertices(&m);
        assert!(!vs.is_empty(), "{what}: no live vertices");
        for &vh in &vs {
            assert!(
                matches!(m.vertex_cache(vh), PointCache::Bounded { .. }),
                "{what}: a chain the first rung decides must be realized, got {:?}",
                m.vertex_cache(vh)
            );
        }
        // The fixture has to *have* a chain, or "the first rung decides it" is about nothing.
        let replayed = vs
            .iter()
            .filter(|&&vh| {
                realize_vertex(&m, vh, Precision::Bits(128)).is_ok_and(|r| !r.is_exact())
            })
            .count();
        assert!(
            replayed > 0,
            "{what}: every vertex realized exactly — no chain was replayed, so this measures nothing"
        );
    }
}

/// ★★★★ **`Ceiling` has two walls, and this builds one of each** — otherwise the variant has no
/// population at all to be wrong about (the census corpus measures **zero** of it: every real model
/// there is two motions deep at most, and every coordinate decides on the first rung).
///
/// - **Cost.** 300 rational translations put the history past `CACHE_REPLAY_COST_CAP`, so the cache
///   road never walks it — `CacheDecline::CostCap`.
/// - **Bits.** 70 turns of 37° stay well inside that cap, so the road *does* walk, and the ladder's
///   first rung cannot name an `f64` for the corners off the axis — `RealizeError::Undecided`.
///
/// ★★ **And the promise the variant makes is asserted, not assumed**: a road willing to pay still
/// answers. That is the whole difference between `Ceiling` and `Unrealized` — one says "ask again
/// with a bigger budget", the other says "there is no road" — and without this line the two would
/// be distinguishable only by which branch produced them.
///
/// ⚠ `validate` runs on the cost-capped model too. Its vertices carry the construction's figure
/// rather than a realization, which is exactly the population whose tolerance loosened when the
/// cache stopped storing a measured residual; the checker must still find it clean.
#[test]
fn a_ceiling_is_reached_by_cost_and_by_bits_and_the_paid_door_still_answers() {
    // ── the cost wall ──────────────────────────────────────────────────────────────────────
    let mut m = Model::new();
    let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    for _ in 0..300 {
        let iso = Isometry::translation([
            Rat::new(1, 7).expect("1/7"),
            Rat::from_int(0),
            Rat::from_int(0),
        ]);
        s = moved(&mut m, s, iso);
    }
    m.rebuild_adjacency();

    let vs = live_vertices(&m);
    assert!(!vs.is_empty(), "no live vertices");
    for &vh in &vs {
        assert!(
            matches!(m.vertex_cache(vh), PointCache::Ceiling { .. }),
            "past the cost cap the cache road stops: {:?}",
            m.vertex_cache(vh)
        );
        assert!(
            matches!(realize_cache(&m, m.vertex(vh)), Err(CacheDecline::CostCap)),
            "and says so by name"
        );
        assert!(
            realize_vertex(&m, vh, Precision::NearestF64)
                .expect("the paid door walks the chain")
                .to_f64()
                .is_some(),
            "a cost-capped vertex is not undecidable, only unpaid-for"
        );
    }
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "a model of Ceiling vertices is still a valid b-rep"
    );

    // ── the bits wall ──────────────────────────────────────────────────────────────────────
    let mut m = Model::new();
    let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    for _ in 0..70 {
        s = moved(&mut m, s, turn(Axis::Z, 37));
    }
    m.rebuild_adjacency();

    let undecided = live_vertices(&m)
        .into_iter()
        .filter(|&vh| {
            matches!(m.vertex_cache(vh), PointCache::Ceiling { .. })
                && matches!(
                    realize_cache(&m, m.vertex(vh)),
                    Err(CacheDecline::Cannot(RealizeError::Undecided))
                )
        })
        .count();
    assert!(
        undecided > 0,
        "70 turns of 37° must leave the first rung undecided somewhere"
    );
}

/// A box carrying `n` recorded translation nodes — `1/7` is not an `f64`, so the statements cannot
/// absorb it and every step records one.
fn translated_chain(n: usize) -> Model {
    let mut m = Model::new();
    let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    for _ in 0..n {
        let iso = Isometry::translation([
            Rat::new(1, 7).expect("1/7"),
            Rat::from_int(0),
            Rat::from_int(0),
        ]);
        s = moved(&mut m, s, iso);
    }
    m.rebuild_adjacency();
    m
}

/// ★★★★ **The paid door raises every `Ceiling`, to exactly what the expensive road says, and
/// touches nothing else.**
///
/// The fixture is 300 rational translations: past the cache road's cost cap, so every live vertex
/// is `Ceiling` — and decidable at the first rung once something is willing to walk the chain, so
/// the door must raise all of them rather than report them as left behind.
#[test]
fn the_refine_door_raises_every_ceiling_to_the_realization() {
    let mut m = translated_chain(300);
    let vs = live_vertices(&m);
    let ceilings = vs
        .iter()
        .filter(|&&vh| matches!(m.vertex_cache(vh), PointCache::Ceiling { .. }))
        .count();
    assert!(ceilings > 0, "the fixture must produce Ceilings to raise");

    // What the expensive road says, asked *before* the door runs so this is an independent oracle
    // rather than a restatement of what the door wrote.
    let want: Vec<[f64; 3]> = vs
        .iter()
        .map(|&vh| {
            realize_vertex(&m, vh, Precision::NearestF64)
                .expect("the paid road realizes a translated corner")
                .to_f64()
                .expect("and names an f64")
                .0
        })
        .collect();

    let report = nacre_ops::refine_vertex_cache(&mut m);
    assert_eq!(report.refined, ceilings, "every Ceiling is raised");
    assert_eq!(report.left, 0, "and none is left behind");
    for (&vh, w) in vs.iter().zip(&want) {
        let PointCache::Bounded { coord, .. } = *m.vertex_cache(vh) else {
            panic!("a raised vertex is Bounded: {:?}", m.vertex_cache(vh));
        };
        assert_eq!(
            coord.as_array(),
            *w,
            "and holds the realization, bit for bit"
        );
    }

    // Idempotent: nothing is a `Ceiling` any more, so there is nothing to pay for.
    assert_eq!(nacre_ops::refine_vertex_cache(&mut m).refined, 0);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "and the refined model is still a valid b-rep"
    );
}

/// ★★★★ **The door carries the edges with it** — the order problem, back again.
///
/// An edge's curve is derived from its endpoints' coordinates when the edge is pushed. That was
/// safe while coordinates never moved after the fact; this door moves them. So the door re-derives,
/// and the proof is that a rebuild *afterwards* finds nothing left to do.
#[test]
fn the_refine_door_carries_the_edges_with_it() {
    let mut m = translated_chain(300);
    let before: Vec<_> = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .map(|h| (h, m.edge(h)))
        .map(|(h, _)| m.edge_curve(h).clone())
        .collect();
    assert!(nacre_ops::refine_vertex_cache(&mut m).refined > 0);

    let after: Vec<_> = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .map(|h| (h, m.edge(h)))
        .map(|(h, _)| m.edge_curve(h).clone())
        .collect();
    assert_ne!(
        before, after,
        "if the curves did not move, this fixture proves nothing about carrying them"
    );

    m.rebuild_edge_cache();
    let again: Vec<_> = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .map(|h| (h, m.edge(h)))
        .map(|(h, _)| m.edge_curve(h).clone())
        .collect();
    assert_eq!(after, again, "the door already re-derived every live curve");
}

/// The write door raises a `Ceiling` and refuses anything else — the type keeping "only ever more
/// accurate" true rather than the caller remembering it.
#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "only a Ceiling is raised")]
fn the_write_door_refuses_a_vertex_that_is_not_a_ceiling() {
    let mut m = boolean_corner();
    let vh = live_vertices(&m)
        .into_iter()
        .find(|&vh| matches!(m.vertex_cache(vh), PointCache::Bounded { .. }))
        .expect("a boolean corner realizes its vertices");
    let coord = m.vertex_point(vh);
    m.refine_vertex_cache(vh, coord, [nacre_scalar::Mag::ZERO; 3]);
}

/// The same chain, built with the prefix table emptied before every step — so every realization
/// folds from the base, the way it did before the accelerator existed.
fn translated_chain_cold(n: usize) -> Model {
    let mut m = Model::new();
    let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    for _ in 0..n {
        m.clear_prefix_hp();
        let iso = Isometry::translation([
            Rat::new(1, 7).expect("1/7"),
            Rat::from_int(0),
            Rat::from_int(0),
        ]);
        s = moved(&mut m, s, iso);
    }
    m.rebuild_adjacency();
    m
}

/// ★★★★ **The accelerator changes the clock and nothing else.**
///
/// Two builds of one chain — one carrying prefixes forward, one emptying the table before every
/// step so each realization folds from the base — and every live vertex's cache must agree bit for
/// bit. That is the whole contract: the fold is a left fold, so resuming from the value after `k`
/// nodes reaches what folding all of them reaches.
///
/// ⚠ This lock is worth only as much as the accelerator's reach, so it also asserts the table was
/// **used**: a run where nothing was ever remembered would pass while measuring nothing.
#[test]
fn the_prefix_table_changes_no_coordinate() {
    let warm = translated_chain(100);
    let cold = translated_chain_cold(100);

    assert!(
        warm.prefix_hp_len() > 0,
        "nothing was remembered — this lock would pass without the accelerator running at all"
    );
    let (vw, vc) = (live_vertices(&warm), live_vertices(&cold));
    assert_eq!(vw.len(), vc.len(), "the two builds disagree on topology");
    assert!(!vw.is_empty());
    for (&a, &b) in vw.iter().zip(&vc) {
        assert_eq!(
            warm.vertex_cache(a),
            cold.vertex_cache(b),
            "vertex {}: carrying the prefix forward moved the coordinate",
            a.index()
        );
    }
}

/// **Emptying the table is always safe** — the realization that follows reaches the same value.
#[test]
fn clearing_the_prefix_table_costs_only_time() {
    let mut m = translated_chain(60);
    let before: Vec<_> = live_vertices(&m)
        .into_iter()
        .map(|vh| (vh, *m.vertex_cache(vh)))
        .collect();
    m.clear_prefix_hp();
    assert_eq!(m.prefix_hp_len(), 0);
    for (vh, cached) in before {
        let (v, _) = realize_vertex(&m, vh, Precision::NearestF64)
            .expect("a translated corner realizes")
            .to_f64()
            .expect("and names an f64");
        assert_eq!(cached.coord().as_array(), v, "vertex {}", vh.index());
    }
}

/// ★★★ **The table is bounded by the live generation, on both axes.**
///
/// Depth does not grow it: a prefix is read by exactly one successor, so the entry that was used is
/// taken as the new one is left. Nor do booleans: an arrangement's vertex is nobody's prefix, and
/// filing it would add an entry per boolean that no lookup could ever hit.
#[test]
fn the_prefix_table_stays_one_generation() {
    for depth in [20usize, 60, 100] {
        let m = translated_chain(depth);
        assert_eq!(
            m.prefix_hp_len(),
            live_vertices(&m).len(),
            "depth {depth}: the table should hold one entry per live vertex, not per node"
        );
    }

    let mut m = Model::new();
    for i in 0..8 {
        let a = cuboid(&mut m, [0.0; 3], [10.0, 10.0, 10.0]);
        let b = cuboid(
            &mut m,
            [3.3 + f64::from(i) * 0.01, 3.3, -1.0],
            [7.7, 7.7, 11.0],
        );
        boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
    }
    m.rebuild_adjacency();
    assert_eq!(
        m.prefix_hp_len(),
        0,
        "a boolean's result vertices are nobody's prefix and must file nothing"
    );
}

/// **Where a motion does not extend a prefix, the accelerator misses and the answer is unchanged.**
///
/// Two of the three shapes are buildable here: a translation the statements absorb (an exact `f64`
/// offset, so no node is recorded and the base itself moves) and a quarter turn about a box's own
/// normal (which fixes each plane, so `transform` restates rather than records). Both must give the
/// same coordinates as the same build with the table emptied throughout.
///
/// ⚠ The third — a triple that straddles two histories and solves through the world road, whose
/// answer carries no leaf — is **not built here**. It needs two operands with different recorded
/// chains meeting on one corner, and I did not get a fixture to stand; saying so is better than a
/// lock that quietly covers two cases while claiming three.
#[test]
fn a_motion_that_does_not_extend_a_prefix_is_unaffected() {
    for (what, exact_shift) in [("absorbed translation", true), ("quarter turn", false)] {
        let build = |clear: bool| {
            let mut m = Model::new();
            let mut s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
            m.rebuild_adjacency();
            for _ in 0..6 {
                if clear {
                    m.clear_prefix_hp();
                }
                let iso = if exact_shift {
                    Isometry::translation([
                        Rat::new(1, 2).expect("1/2"),
                        Rat::from_int(0),
                        Rat::from_int(0),
                    ])
                } else {
                    turn(Axis::Z, 90)
                };
                s = moved(&mut m, s, iso);
            }
            m.rebuild_adjacency();
            m
        };
        let (warm, cold) = (build(false), build(true));
        // ★ The guard that keeps this from being a third copy of the lock above. Neither motion
        // records a node — one is absorbed into the statements, the other fixes every plane it
        // touches — so there is no chain, nothing is ever filed, and the miss is structural rather
        // than incidental. If either premise were wrong the table would be non-empty here and this
        // test would be measuring ordinary hits while claiming to measure exceptions.
        assert_eq!(
            warm.prefix_hp_len(),
            0,
            "{what}: a motion that records no node cannot file a prefix"
        );
        let (vw, vc) = (live_vertices(&warm), live_vertices(&cold));
        assert!(!vw.is_empty(), "{what}: no live vertices");
        for (&a, &b) in vw.iter().zip(&vc) {
            assert_eq!(
                warm.vertex_cache(a),
                cold.vertex_cache(b),
                "{what}: vertex {} moved",
                a.index()
            );
        }
    }
}
