use super::*;
use nacre_math::Vector3;
use proptest::prelude::*;
use std::collections::HashSet;

fn cube(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(&mut m, Point3::from_array(min), Point3::from_array(max));
    m
}

fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
    let mut m = Model::new();
    nacre_ops::fixtures::cylinder(
        &mut m,
        Point3::from_array(base),
        Vector3::from_array(axis),
        r,
        h,
    );
    m
}

/// A closed manifold mesh has every undirected triangle edge shared by
/// exactly two triangles. Returns the count of edges that are not.
fn non_watertight_edges(t: &Tessellation) -> usize {
    let mut counts: std::collections::HashMap<(u32, u32), usize> = std::collections::HashMap::new();
    for (_, tri) in t.triangles.iter() {
        let [a, b, c] = tri.vertices.map(|h| h.index());
        for (x, y) in [(a, b), (b, c), (c, a)] {
            let key = (x.min(y), x.max(y));
            *counts.entry(key).or_default() += 1;
        }
    }
    counts.values().filter(|&&c| c != 2).count()
}

fn v_lines(obj: &str) -> Vec<&str> {
    obj.lines().filter(|l| l.starts_with("v ")).collect()
}
fn f_lines(obj: &str) -> Vec<&str> {
    obj.lines().filter(|l| l.starts_with("f ")).collect()
}
fn parse3(rest: &str) -> Vec<f64> {
    rest.split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

/// **What the OBJ writer owes, now that it is the only one.**
///
/// ★ This used to run the bootstrap `to_obj(&Model)`, whose contract was *"every model vertex
/// is written in store order, so an OBJ index stays a vertex handle's index"*. That writer is
/// gone and the contract with it — an index is now a **mesh** vertex, and a mesh vertex says
/// far more than a handle did ([`TessOrigin`] names the vertex, edge or face it came from).
/// What is left to check here is the text: one `v` per mesh vertex, one `f` per triangle,
/// 1-based, and nothing referenced that was not written.
#[test]
fn unit_cube_obj_shape() {
    let t = tessellate(
        &cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        &TessConfig::default(),
    )
    .unwrap();
    let obj = t.to_obj();
    assert_eq!(v_lines(&obj).len(), t.vertices.len());
    assert_eq!(f_lines(&obj).len(), t.triangles.len());
    assert!(obj.lines().next().unwrap().starts_with('#'));
    assert_eq!(v_lines(&obj).len(), 8);

    let faces = f_lines(&obj);
    assert_eq!(faces.len(), 12); // 6 quads × 2 triangles
    let mut used = HashSet::new();
    for f in faces {
        let idx = parse3(&f[2..]);
        assert_eq!(idx.len(), 3);
        for &i in &idx {
            assert!((1.0..=8.0).contains(&i));
            used.insert(i as u32);
        }
    }
    assert_eq!(used.len(), 8); // every vertex referenced
}

/// The eight corners come back as text — **as a set**.
///
/// ★ Ordered comparison would pass today (☑ measured: the mesh's vertices come out in the
/// model's order for a cuboid, because `sample_edge` walks the edge store and dedups on first
/// sight). It is not a contract, though — nothing promises that traversal — so asserting it
/// would pin an artifact.
#[test]
fn vertices_round_trip() {
    let obj = tessellate(
        &cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]),
        &TessConfig::default(),
    )
    .unwrap()
    .to_obj();
    let got: Vec<[f64; 3]> = v_lines(&obj)
        .into_iter()
        .map(|l| {
            let c = parse3(&l[2..]);
            [c[0], c[1], c[2]]
        })
        .collect();
    let mut got = got;
    got.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut expected = vec![
        [-2.0, 1.0, 0.0],
        [3.0, 1.0, 0.0],
        [3.0, 4.0, 0.0],
        [-2.0, 4.0, 0.0],
        [-2.0, 1.0, 10.0],
        [3.0, 1.0, 10.0],
        [3.0, 4.0, 10.0],
        [-2.0, 4.0, 10.0],
    ];
    expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(got, expected);
}

// --- provenance tessellation ---

#[test]
fn cube_tessellates_watertight() {
    let t = tessellate(
        &cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        &TessConfig::default(),
    )
    .unwrap();
    assert_eq!(t.vertices.len(), 8);
    assert_eq!(t.triangles.len(), 12); // 6 quads × 2
    assert_eq!(non_watertight_edges(&t), 0);
    for (_, v) in t.vertices.iter() {
        assert!(matches!(v.origin, TessOrigin::OnVertex(_)));
    }
}

#[test]
fn cylinder_tessellates_watertight() {
    let cfg = TessConfig::default();
    let n = circle_segments(&cfg, 2.0);
    let t = tessellate(&cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0), &cfg).unwrap();
    assert_eq!(t.vertices.len(), 2 * n); // two rim rings, seam vertices shared
    assert_eq!(t.triangles.len(), 4 * n - 4); // 2 caps (n−2) + band (2n)
    assert_eq!(non_watertight_edges(&t), 0);
}

#[test]
fn cylinder_provenance_matches_positions() {
    let m = cylinder([1.0, -2.0, 0.5], [0.0, 0.0, 1.0], 2.0, 5.0);
    let t = tessellate(&m, &TessConfig::default()).unwrap();
    for (_, tv) in t.vertices.iter() {
        let expected = match tv.origin {
            TessOrigin::OnVertex(v) => m.vertex_point(v),
            TessOrigin::OnEdge { edge, t: param } => match m.edge_curve(edge) {
                Curve::Circle(c) => c.point_at(param),
                Curve::Line(_) => unreachable!("cylinder rims are circles"),
            },
            TessOrigin::OnFace { .. } => unreachable!("cylinder uses no interior face samples"),
        };
        assert!((tv.pos - expected).norm() <= 1e-9);
    }
}

/// ★ **The radius is 20 so that the sagitta budget is the one being measured.**
/// At `r = 2`, with the angular budget in place, both tolerances answer the same
/// 180 segments — the test would compare a number with itself and pass on any
/// `tol` at all. The angular term leads until
/// `r ≈ 66·(tol/0.01)`, so a radius past that is where "finer tolerance" still
/// means something.
#[test]
fn finer_tolerance_adds_triangles() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 20.0, 5.0);
    let at = |tol| {
        tessellate(
            &m,
            &TessConfig {
                tol,
                ..Default::default()
            },
        )
        .unwrap()
        .triangles
        .len()
    };
    let (coarse, fine) = (at(0.5), at(0.001));
    assert!(fine > coarse, "fine {fine} should exceed coarse {coarse}");
}

/// ★★ **The promise, measured on the mesh rather than recomputed from the rule.**
///
/// Every turn between consecutive chords of a rim must be within the angular
/// budget — that is what "the circle does not look polygonal" means, and it is a
/// property of the *output*, so it survives a change of formula.
#[test]
fn no_chord_turns_more_than_the_angular_budget() {
    let cfg = TessConfig::default();
    for r in [0.1, 0.5, 3.0, 20.0] {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], r, 5.0);
        let t = tessellate(&m, &cfg).unwrap();
        for (&eh, ring) in t.by_edge.iter() {
            if ring.len() < 3 || matches!(m.edge_curve(eh), Curve::Line(_)) {
                continue; // a straight edge has no turn to measure, however many samples
            }
            // Wrapping is right because every multi-point *curved* polyline here is a
            // closed rim; an arc would need the two end turns left out. A straight
            // edge can carry a third sample too (the bridge pre-pass), and is skipped above.
            let p: Vec<Point3> = ring.iter().map(|&h| t.vertices.get(h).pos).collect();
            for i in 0..p.len() {
                let (a, b, c) = (p[i], p[(i + 1) % p.len()], p[(i + 2) % p.len()]);
                let (u, v) = ((b - a).normalize(), (c - b).normalize());
                let (Some(u), Some(v)) = (u, v) else { continue };
                let turn = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
                assert!(
                    turn <= cfg.max_angle_deg + 1e-9,
                    "r={r}: a chord turned {turn}°, past the {}° budget",
                    cfg.max_angle_deg
                );
            }
        }
    }
}

/// ★★★★ **No triangle may leave the surface it approximates.**
///
/// Watertightness cannot see this: a mesh whose triangles cut *through* the cylinder still
/// uses every edge twice and still counts right. What bounds it is the **sag** — how far a
/// chord spanning `Δθ` falls inside the surface, `r(1 − cos(Δθ/2))` — and the budget is the
/// one the edge sampler already works to, so a face may not undo with a diagonal what the
/// boundary paid for.
///
/// ★ This is the assertion the *chart* road owes. The sweep is free to draw a diagonal
/// between any two boundary vertices, and if the chart were laid out with the sweep running
/// **along the axis** instead of around it (the handedness repair swapping `u` and `v` would
/// do exactly that), the diagonals would span wide arcs and every other check here would
/// still pass.
#[test]
fn no_triangle_leaves_the_cylinder() {
    let cfg = TessConfig::default();
    for r in [0.1, 0.5, 3.0, 20.0] {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], r, 5.0);
        let t = tessellate(&m, &cfg).unwrap();
        let step = std::f64::consts::TAU / circle_segments(&cfg, r) as f64;
        let budget = r * (1.0 - (step / 2.0).cos());
        for (_, tri) in t.triangles.iter() {
            let f = m.face(tri.face);
            let Surface::Cylinder(cy) = m.surface_cache(f.surface) else {
                continue;
            };
            let (o, w) = (cy.axis().origin(), cy.axis().direction());
            let x = cy.ref_dir();
            let y = w.cross(x);
            let ang = |h: Handle<TessVertex>| {
                let d = t.vertices.get(h).pos - o;
                d.dot(y).atan2(d.dot(x))
            };
            for k in 0..3 {
                let (a, b) = (ang(tri.vertices[k]), ang(tri.vertices[(k + 1) % 3]));
                let mut d = (a - b).abs();
                if d > std::f64::consts::PI {
                    d = std::f64::consts::TAU - d;
                }
                let sag = r * (1.0 - (d / 2.0).cos());
                assert!(
                    sag <= budget + 1e-12,
                    "r={r}: a chord spanning {}° sags {sag}, past the {budget} the boundary \
                         is sampled to",
                    d.to_degrees()
                );
            }
        }
    }
}

/// Build the chart `triangulate_face` would build, for one face.
fn chart_of(t: &Tessellation, m: &Model, cfg: &TessConfig, fh: Handle<Face>) -> Chart {
    let face = m.face(fh);
    match m.surface_cache(face.surface) {
        Surface::Plane(p) => planar_chart(t, face, p).unwrap(),
        Surface::Cylinder(c) => {
            cylinder_chart(t, m, cfg, face, c, whole_rims_cut(t, m, face).map(Some)).unwrap()
        }
    }
}

/// ★★★★ **The chart follows the ring's own normal, on every plane — and the triangles follow
/// the ring.**
///
/// Two claims about the projection — *"the projector must follow the ring's own normal, not a
/// coordinate convention"* and *"reversing the ring flips the Newell normal, and the same
/// triangles come out with the opposite winding"* — stated where the projection rule lives
/// once ([`planar_chart`]), on a cube that states them **six times at once**, one per
/// axis-aligned plane, with `drop_axis` returning each of its three answers.
///
/// ★ Not just "a ring and its reversal disagree": this says which way is right — every face's
/// triangles wind about the face's **outward** normal, which is the property the mesh's
/// consumers actually rely on.
#[test]
fn every_planar_face_meshes_about_its_own_normal() {
    let cfg = TessConfig::default();
    let m = cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
    let t = tessellate(&m, &cfg).unwrap();
    let reach = m.reachable();
    let mut faces = 0;
    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
        if !reach.faces.contains(&fh) {
            continue;
        }
        let Surface::Plane(plane) = m.surface_cache(face.surface) else {
            continue;
        };
        let tris = &t.by_face[&fh];
        assert_eq!(tris.len(), 2, "a rectangle is two triangles");
        let sign = match face.orientation {
            nacre_topo::Orientation::Forward => 1.0,
            nacre_topo::Orientation::Reversed => -1.0,
        };
        let mut area = 0.0;
        for &th in tris {
            let [a, b, c] = t.triangles.get(th).vertices.map(|h| t.vertices.get(h).pos);
            let cr = (b - a).cross(c - a);
            assert!(
                cr.dot(plane.normal()) * sign > 0.0,
                "a triangle winds against its face's outward normal"
            );
            area += 0.5 * cr.norm();
        }
        // 5 × 3 × 10: two faces of each of the three rectangle shapes.
        assert!(
            [15.0, 50.0, 30.0].iter().any(|w| (area - w).abs() < 1e-12),
            "face area {area}"
        );
        faces += 1;
    }
    assert_eq!(faces, 6, "all six planes, so all three `drop_axis` answers");
}

/// ★★★★★ **A chart is a bijection, and the way back has to be the way back.**
///
/// The chart may now invent an interior point, and an invented point's position comes from
/// nowhere else — so [`ChartMap`] is load-bearing. Both arms undo a **handedness repair**, and
/// that is exactly the sort of thing that is silently wrong: a plane swaps its two axes, a
/// cylinder negates its height, and either undone in the wrong direction still produces a
/// plausible mesh in the wrong place. Measured where the answer is independently known — on the
/// boundary, whose points came from the model rather than from this function.
#[test]
fn a_chart_maps_back_to_the_point_it_flattened() {
    let cfg = TessConfig::default();
    for (name, m) in [
        ("cube", cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])),
        (
            "cylinder",
            cylinder([1.0, -2.0, 0.5], [0.0, 0.0, 1.0], 2.0, 5.0),
        ),
        (
            "tilted cylinder",
            cylinder([1.0, -2.0, 0.5], [1.0, 2.0, 3.0], 2.0, 5.0),
        ),
    ] {
        let t = tessellate(&m, &cfg).unwrap();
        let reach = m.reachable();
        let mut i = 0u32;
        while let Some(fh) = m.face_handle_at(i) {
            i += 1;
            if !reach.faces.contains(&fh) {
                continue;
            }
            let chart = chart_of(&t, &m, &cfg, fh);
            for (i, &h) in chart.handles.iter().enumerate() {
                let (params, back) = chart.map.invert(chart.uv[i]);
                let want = t.vertices.get(h).pos;
                assert!(
                    (back - want).norm() <= 1e-9,
                    "{name}: chart {:?} came back {back:?}, not {want:?}",
                    chart.uv[i]
                );
                // And the parameters name the same point through the surface's own evaluator.
                if let Surface::Cylinder(c) = m.surface_cache(m.face(fh).surface) {
                    assert!((c.point_at(params[0], params[1]) - want).norm() <= 1e-9);
                }
            }
        }
    }
}

/// ★★★★ **A band's rims already sample its curvature, so the chart offers nothing.**
///
/// Which is a claim about *stations*, and stations are where this can go wrong quietly: they
/// are read from each boundary circle's **centre** precisely because one rim's sampled points
/// disagree on that coordinate by ulps, and reading those instead would turn one rim into 181
/// stations and the lattice into 180 × 181 points. The lowest and highest are the face's own
/// rims and are dropped. An ordinary cylinder has exactly those two — so: zero.
#[test]
fn an_ordinary_band_is_offered_no_interior_points() {
    let cfg = TessConfig::default();
    for axis in [[0.0, 0.0, 1.0], [1.0, 2.0, 3.0]] {
        let m = cylinder([1.0, -2.0, 0.5], axis, 2.0, 5.0);
        let t = tessellate(&m, &cfg).unwrap();
        let mut i = 0u32;
        while let Some(fh) = m.face_handle_at(i) {
            i += 1;
            let face = m.face(fh);
            if !matches!(m.surface_cache(face.surface), Surface::Cylinder(_)) {
                continue;
            }
            let chart = chart_of(&t, &m, &cfg, fh);
            assert!(
                chart.interior.is_empty(),
                "axis {axis:?}: a plain band was offered {} points",
                chart.interior.len()
            );
        }
        // …and none was minted, which is the same claim read off the output.
        assert!(
            !t.vertices
                .iter()
                .any(|(_, v)| matches!(v.origin, TessOrigin::OnFace { .. }))
        );
    }
}

/// ★ **Scale invariance**: the angular budget is a *relative* one, so a small
/// circle and a large one are cut into the same number of pieces. Under an
/// absolute-only rule these are 10 and 100 — the asymmetry that makes small holes
/// look like polygons.
#[test]
fn a_small_circle_is_cut_as_finely_as_a_large_one() {
    let cfg = TessConfig::default();
    let n = |r| circle_segments(&cfg, r);
    assert_eq!(n(0.2), n(20.0));
    assert_eq!(n(0.2), 180, "360°/2°");
    assert!(n(0.1) > 8, "a small circle is not held to an octagon floor");
    // Past the crossover the sagitta budget leads and asks for more.
    assert!(n(1000.0) > n(20.0));
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn arbitrary_box_obj_shape(
        min in prop::array::uniform3(-1e3f64..1e3),
        ext in prop::array::uniform3(1e-2f64..1e3),
    ) {
        // Floor and height, as drawn: two random corners differ by no decimal an f64 carries.
        let mut m = Model::new();
        nacre_ops::fixtures::cuboid_on(&mut m, [min[0], min[1]], [min[0] + ext[0], min[1] + ext[1]], min[2], ext[2]);
        let obj = tessellate(&m, &TessConfig::default()).unwrap().to_obj();
        prop_assert_eq!(v_lines(&obj).len(), 8);
        let faces = f_lines(&obj);
        prop_assert_eq!(faces.len(), 12);
        for f in faces {
            let idx = parse3(&f[2..]);
            prop_assert_eq!(idx.len(), 3);
            for &i in &idx {
                prop_assert!((1.0..=8.0).contains(&i));
            }
        }
    }

    #[test]
    fn arbitrary_cylinder_watertight(
        axis in prop::array::uniform3(-1.0f64..1.0),
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        prop_assume!(Vector3::from_array(axis).norm() > 0.1);
        // ★ A **coarse** angle on purpose: this case is about the mesh's
        // *structure* — the counts and watertightness — which do not depend on how
        // finely the rim is cut. The 2° default puts 180 segments on every radius,
        // which made this property 8.5× slower (0.92s → 7.8s, measured) for no
        // coverage at all. Coarse also lets the sagitta budget lead at the larger
        // radii, so `n` still varies across the cases.
        let cfg = TessConfig {
            max_angle_deg: 20.0,
            ..Default::default()
        };
        let n = circle_segments(&cfg, r);
        let t = tessellate(&cylinder([0.0, 0.0, 0.0], axis, r, h), &cfg).unwrap();
        prop_assert_eq!(t.vertices.len(), 2 * n);
        prop_assert_eq!(t.triangles.len(), 4 * n - 4);
        prop_assert_eq!(non_watertight_edges(&t), 0);
    }
}

/// ★★★★★ **A band whose rims are both cut meshes, and its mesh is the surface.**
///
/// [`nacre_ops::fixtures::cut_at_both_rims`] builds a lateral that wraps the axis between two cut
/// rims with no seam edge — one rim the outer loop, the other an inner one. Its two rings unroll to
/// two open lines, not a polygon, so the chart cuts the band along one generator
/// ([`cut_seamless_bands`]). Three placements of that generator are asked for: an ordinary one
/// (`y_min = −2`), boxes standing on the plane through the axis (`y_min = 0` — each rim's `θ = 0`
/// point is a box corner's pierce, not a seam vertex), and a band with a hole cut where its seam
/// would be (`θ = 0`), which the generator must step around. Each must mesh watertight with the
/// lateral's mesh area equal to the analytic area (design 「메시로 재기」: `validate`, watertightness
/// and the volume all stay green when a curved face's mesh is wrong — the area is the oracle).
#[test]
fn a_band_cut_at_both_rims_meshes_as_its_surface() {
    use nacre_ops::{BoolKind, boolean, fixtures};
    let pi = std::f64::consts::PI;
    for built in [BoolKind::Fuse, BoolKind::Cut] {
        for (y_min, holed) in [(-2.0, false), (0.0, false), (-2.0, true)] {
            let mut m = Model::new();
            let mut band = fixtures::cut_at_both_rims(&mut m, built, y_min);
            // Each bite takes the rim's arc with `x > 0.5` and `y ≥ y_min`, half a unit tall.
            let bite = if y_min < 0.0 {
                2.0 * pi / 3.0
            } else {
                pi / 3.0
            };
            let mut want = 2.0 * pi * 4.0 - 2.0 * bite * 0.5;
            if holed {
                let window = fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.8, -0.2, 1.5]),
                    Point3::from_array([1.2, 0.2, 2.5]),
                );
                m.rebuild_adjacency();
                band = boolean(&mut m, BoolKind::Cut, band, window).expect("the window")[0];
                m.rebuild_adjacency();
                want -= 2.0 * 0.2f64.asin();
            }
            let lateral = m
                .shell(m.solid(band).outer)
                .faces
                .iter()
                .copied()
                .find(|&f| matches!(m.surface_cache(m.face(f).surface), Surface::Cylinder(_)))
                .expect("one lateral face");
            let face = m.face(lateral);
            // The premise: the shape this lock is about — two wrapping rims (and the hole) as
            // loops, no seam edge anywhere on the face.
            assert_eq!(
                face.inner.len(),
                if holed { 2 } else { 1 },
                "{built:?} {y_min} {holed}"
            );
            assert!(
                std::iter::once(&face.outer)
                    .chain(&face.inner)
                    .flat_map(|l| &l.half_edges)
                    .all(|he| {
                        let [a, b] = m.edge(he.edge).surfaces;
                        a != b
                    }),
                "no seam edge"
            );
            let t = tessellate(&m, &TessConfig::default())
                .unwrap_or_else(|e| panic!("{built:?} y_min {y_min} holed {holed}: {e:?}"));
            assert_eq!(non_watertight_edges(&t), 0, "{built:?} {y_min} {holed}");
            let area: f64 = t.by_face[&lateral]
                .iter()
                .map(|&th| {
                    let [a, b, c] = t.triangles.get(th).vertices.map(|v| t.vertices.get(v).pos);
                    (b - a).cross(c - a).norm() / 2.0
                })
                .sum();
            assert!(
                (area - want).abs() < 1e-3 * want,
                "{built:?} y_min {y_min} holed {holed}: mesh {area} vs {want}"
            );
        }
    }
}

/// ★★★★ **The band's generator stays out of a hole, and passes each rim once.** Synthetic turns,
/// where the choice is decided by the rule and not by where samples fall: both wrapping rings
/// sampled every 10°, so every candidate (the middle of a step, 5° off the grid) is equally far
/// from every sample; a hole spanning `[−100°, 100°]` with vertices only at its ends. The first
/// candidate (5°) lies inside the hole, and the rule must pass it over — a cut through a hole
/// crosses its boundary. Then a second ring that doubles back over `[40°, 60°]` passes those
/// generators three times, and the rule must skip them too.
#[test]
fn a_band_generator_keeps_out_of_holes_and_passes_each_rim_once() {
    let deg = |d: f64| d.to_radians();
    let tau = std::f64::consts::TAU;
    let up: Vec<f64> = (0..36).map(|k| deg(10.0 * f64::from(k))).collect();
    let down: Vec<f64> = (0..36).map(|k| deg(360.0 - 10.0 * f64::from(k))).collect();
    let hole = vec![deg(-100.0), deg(100.0)];
    let turns = vec![(up.clone(), tau), (down.clone(), -tau), (hole, 0.0)];
    let c = band_generator(&turns, [0, 1]).expect("a generator outside the hole");
    let c = c.to_degrees().rem_euclid(360.0);
    assert!(
        c > 100.0 && c < 260.0,
        "the generator {c}° lies in the hole"
    );

    // The second ring doubles back: 10° up to 30°, back to 10°, on to 360° — every θ in (10°, 30°)
    // is passed three times. The first ring starts at 10°, so its first candidates (15°, 25°) lie
    // in the fold and the rule must skip them.
    let from_ten: Vec<f64> = (1..=36).map(|k| deg(10.0 * f64::from(k))).collect();
    let mut doubled: Vec<f64> = [10.0, 20.0, 30.0, 20.0, 10.0].map(deg).to_vec();
    doubled.extend((4..=36).map(|k| deg(10.0 * f64::from(k)))); // 40° .. 360°
    let turns = vec![(from_ten, tau), (doubled, tau)];
    let c = band_generator(&turns, [0, 1]).expect("a generator off the fold");
    let c = c.to_degrees().rem_euclid(360.0);
    assert!(
        !(10.0..=30.0).contains(&c),
        "the generator {c}° is passed three times"
    );
}

/// A solid's cylindrical faces, in shell order.
fn laterals(m: &Model, solid: Handle<nacre_topo::Solid>) -> Vec<Handle<Face>> {
    let s = m.solid(solid);
    std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .flat_map(|sh| m.shell(sh).faces.clone())
        .filter(|&f| matches!(m.surface_cache(m.face(f).surface), Surface::Cylinder(_)))
        .collect()
}

/// The mesh area of one face.
fn face_area(t: &Tessellation, fh: Handle<Face>) -> f64 {
    t.by_face[&fh]
        .iter()
        .map(|&th| {
            let [a, b, c] = t.triangles.get(th).vertices.map(|v| t.vertices.get(v).pos);
            (b - a).cross(c - a).norm() / 2.0
        })
        .sum()
}

/// ★★★★★ **A lateral bounded by two whole rims and nothing else is cut at their own vertices.**
///
/// A bore through a plate: both rims are one closed edge, so the band is cut at their vertices
/// ([`cut_seamless_bands`]'s first rule) and nothing is inserted — the mesh has exactly the rims'
/// samples, as a plain cylinder's (`cylinder_tessellates_watertight`). Before that rule such a band
/// went to the generator, whose cut could land on a rim's closing step (a closed edge's polyline
/// does not repeat its first sample) and fail the whole model's mesh.
#[test]
fn a_bore_is_cut_at_its_rims_vertices() {
    use nacre_ops::{BoolKind, boolean, fixtures};
    let pi = std::f64::consts::PI;
    let mut m = Model::new();
    let plate = fixtures::cuboid(
        &mut m,
        Point3::from_array([-2.0, -2.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    let pin = fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let solid = boolean(&mut m, BoolKind::Cut, plate, pin).expect("the bore")[0];
    m.rebuild_adjacency();
    let [bore] = laterals(&m, solid)[..] else {
        panic!("one lateral");
    };
    let face = m.face(bore);
    assert!(
        face.inner.len() == 1
            && std::iter::once(&face.outer)
                .chain(&face.inner)
                .all(|l| l.half_edges.len() == 1),
        "the premise: two loops, each one closed rim"
    );
    let cfg = TessConfig::default();
    let t = tessellate(&m, &cfg).unwrap();
    assert_eq!(non_watertight_edges(&t), 0);
    let n = circle_segments(&cfg, 1.0);
    let on_lateral: HashSet<_> = t.by_face[&bore]
        .iter()
        .flat_map(|&th| t.triangles.get(th).vertices)
        .collect();
    assert_eq!(
        on_lateral.len(),
        2 * n,
        "the rims' samples and nothing else"
    );
    let want = 2.0 * pi;
    let area = face_area(&t, bore);
    assert!((area - want).abs() < 1e-3 * want, "{area} vs {want}");
}

/// ★★★★ **With a hole beside them, two whole rims are cut by the generator, not at their vertices.**
///
/// The windowed boss with its window on the seam and off it: three loops, so [`cut_seamless_bands`]'s first rule does not apply — where the window covers `θ = 0`,
/// a cut at the rims' own vertices would run through it. The generator steps around the window.
#[test]
fn a_seamless_lateral_with_a_window_is_cut_beside_it() {
    use nacre_ops::fixtures;
    let pi = std::f64::consts::PI;
    for seam_x in [1.0, -1.0] {
        let mut m = Model::new();
        let boss = fixtures::windowed_boss(&mut m, 0.6, seam_x);
        let [new] = laterals(&m, boss)[..] else {
            panic!("one lateral");
        };
        assert_eq!(m.face(new).inner.len(), 2, "seam {seam_x}: rim and window");
        let t = tessellate(&m, &TessConfig::default())
            .unwrap_or_else(|e| panic!("seam {seam_x}: {e:?}"));
        assert_eq!(non_watertight_edges(&t), 0, "seam {seam_x}");
        let want = 2.0 * pi * 4.0 - 2.0 * 0.6f64.asin();
        let area = face_area(&t, new);
        assert!(
            (area - want).abs() < 1e-3 * want,
            "seam {seam_x}: {area} vs {want}"
        );
    }
}

/// ★★★★ **A whole rim beside a cut one, with a pierce point on `θ = 0`.** A boss standing on a
/// plate's edge — the cylinder `r = 1` over `z ∈ [0, 4]` fused with the box `y ∈ [0, 2]`,
/// `z ∈ [0, 2]`, whose wall `y = 0` runs through the axis: the lateral below `z = 2` is the `y < 0`
/// half, so its lower rim is a chain with corners at `θ = 0` and `θ = π`, its upper rim whole. A cut
/// that put a sample into the chain at the whole rim's own `θ` would meet that corner.
#[test]
fn a_seamless_lateral_with_a_chain_rim_on_the_seam_is_cut_off_its_corner() {
    use nacre_ops::{BoolKind, boolean, fixtures};
    let pi = std::f64::consts::PI;
    let mut m = Model::new();
    let boss = fixtures::cylinder(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        4.0,
    )
    .solid;
    let plate = fixtures::cuboid(
        &mut m,
        Point3::from_array([-2.0, 0.0, 0.0]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    m.rebuild_adjacency();
    let fused = boolean(&mut m, BoolKind::Fuse, boss, plate).expect("the boss on the edge")[0];
    m.rebuild_adjacency();
    let [new] = laterals(&m, fused)[..] else {
        panic!("one lateral");
    };
    let t = tessellate(&m, &TessConfig::default()).unwrap();
    assert_eq!(non_watertight_edges(&t), 0);
    let want = 6.0 * pi;
    let area = face_area(&t, new);
    assert!((area - want).abs() < 1e-3 * want, "{area} vs {want}");
}

/// The tube with four windows — `±x`, `±y`, each `|·| < 0.8` across and half a unit tall, at four
/// heights — whose angular spans (`2·asin 0.8` ≈ 106° each, a quarter turn apart) together cover
/// every `θ`, so no generator misses them all.
fn windows_round_a_boss(m: &mut Model) -> Handle<nacre_topo::Solid> {
    use nacre_ops::{BoolKind, boolean, fixtures};
    let mut s = fixtures::cylinder(
        m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        4.0,
    )
    .solid;
    let windows: [([f64; 3], [f64; 3]); 4] = [
        ([0.5, -0.8, 0.5], [2.0, 0.8, 1.0]),
        ([-0.8, 0.5, 1.3], [0.8, 2.0, 1.8]),
        ([-2.0, -0.8, 2.1], [-0.5, 0.8, 2.6]),
        ([-0.8, -2.0, 2.9], [0.8, -0.5, 3.4]),
    ];
    for (lo, hi) in windows {
        let w = fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, s, w).expect("a window");
        assert_eq!(out.len(), 1, "a window leaves one body");
        s = out[0];
        m.rebuild_adjacency();
    }
    s
}

/// ★★★★★ **Windows that cover every angle between them are cut through.** [`windows_round_a_boss`]:
/// every generator crosses a window, so the cut crosses one — twice, its two runs spliced into the
/// outer ring ([`joined_band`]) — rather than refusing the face, which would fail the whole model's
/// mesh.
#[test]
fn a_seamless_lateral_whose_windows_cover_every_angle_is_cut_through_one() {
    let pi = std::f64::consts::PI;
    let mut m = Model::new();
    let boss = windows_round_a_boss(&mut m);
    let cfg = TessConfig::default();
    let [new] = laterals(&m, boss)[..] else {
        panic!("one lateral");
    };
    assert_eq!(m.face(new).inner.len(), 5, "a rim and four windows");
    let t = tessellate(&m, &cfg).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(non_watertight_edges(&t), 0);
    let want = 2.0 * pi * 4.0 - 4.0 * 2.0 * 0.8f64.asin() * 0.5;
    let area = face_area(&t, new);
    assert!((area - want).abs() < 1e-3 * want, "{area} vs {want}");
}

/// ★★★★ **A cut on a whole rim's closing step goes at the end of its polyline.** A closed edge's
/// polyline does not repeat its first sample, so the ring of a loop that is one closed edge has no
/// position the closing step arrives through. The plain cylinder, cut at the middle of its first rim's closing step — where the generator can land when every candidate is
/// as clear as the next.
#[test]
fn a_band_cut_on_a_closed_rims_closing_step_appends_to_its_polyline() {
    let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 4.0);
    let [fh] = laterals(&m, m.live_solids()[0])[..] else {
        panic!("one lateral");
    };
    let cfg = TessConfig::default();
    let mut t = Tessellation::default();
    let reach = m.reachable();
    sample_live_edges(&mut t, &m, &cfg, &reach);
    let face = m.face(fh).clone();
    let Surface::Cylinder(cyl) = m.surface_cache(face.surface) else {
        panic!("a lateral");
    };
    let traced: Vec<_> = std::iter::once(&face.outer)
        .chain(&face.inner)
        .map(|lp| traced_ring(&t, lp))
        .collect();
    assert!(traced[0].1.is_none(), "the premise: no closing position");
    let turns: Vec<_> = traced
        .iter()
        .map(|(r, _)| unwrapped_thetas(&t, cyl, &r.iter().map(|p| p.handle).collect::<Vec<_>>()))
        .collect();
    let (p, q) = *ring_steps(&turns[0]).last().unwrap();
    let rim = face.outer.half_edges[0].edge;
    let before = t.by_edge[&rim].len();
    let cut = band_cut(&mut t, &m, cyl, &traced, &turns, [0, 1], (p + q) / 2.0)
        .expect("a cut on the closing step");
    assert_eq!(t.by_edge[&rim].len(), before + 1);
    assert_eq!(*t.by_edge[&rim].last().unwrap(), cut.rims[0], "appended");
}

/// ★★★★ **The generator crosses the fewest holes it can, and none it cannot splice.** Synthetic
/// turns, as [`a_band_generator_keeps_out_of_holes_and_passes_each_rim_once`]: two holes whose
/// spans overlap and together cover the circle — every generator crosses one, and the cut must
/// cross exactly one, twice; and a hole that folds back over the whole circle, which every
/// generator crosses four times — no cut, because severing it so would leave more than one outer
/// polygon.
#[test]
fn a_band_generator_crosses_the_fewest_holes_and_never_one_four_times() {
    let deg = |d: f64| d.to_radians();
    let tau = std::f64::consts::TAU;
    let up: Vec<f64> = (0..36).map(|k| deg(10.0 * f64::from(k))).collect();
    let down: Vec<f64> = (0..36).map(|k| deg(360.0 - 10.0 * f64::from(k))).collect();
    let covering = vec![
        (up.clone(), tau),
        (down.clone(), -tau),
        (vec![deg(-100.0), deg(100.0)], 0.0),
        (vec![deg(80.0), deg(280.0)], 0.0),
    ];
    let c = band_generator(&covering, [0, 1]).expect("a generator through one hole");
    let c = c.to_degrees().rem_euclid(360.0);
    assert!(
        !(80.0..=100.0).contains(&c),
        "the generator {c}° crosses both holes"
    );
    let folded = vec![
        (up, tau),
        (down, -tau),
        (vec![deg(1.0), deg(359.0), deg(3.0), deg(357.0)], 0.0),
    ];
    assert_eq!(band_generator(&folded, [0, 1]), None);
}
