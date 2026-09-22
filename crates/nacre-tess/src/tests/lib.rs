use super::*;
use nacre_math::Vector3;
use proptest::prelude::*;
use std::collections::HashSet;

fn cube(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
    m
}

fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
    let mut m = Model::new();
    m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
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
        Surface::Cylinder(c) => cylinder_chart(t, m, cfg, face, c).unwrap(),
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
    assert!(n(0.1) > 8, "a small circle is no longer the octagon floor");
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
        let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
        let obj = tessellate(&cube(min, max), &TessConfig::default()).unwrap().to_obj();
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
