//! Cross-crate composition: ops (replay) → validate → step/tess export.
//!
//! The per-crate unit tests each cover one layer; this guards that the layers
//! compose. STEP format details (schema token, entity list, step_io round-trip)
//! stay owned by `nacre-step`'s own tests — here we only prove the pipeline hands
//! a valid model through to real STEP/OBJ output.

use nacre_math::{Point2, Point3};
use nacre_ops::{Operation, Profile2d, SketchPlane, replay};

/// A regular hexagon of the given radius, centred on the sketch origin.
fn hexagon(r: f64) -> Profile2d {
    let points = (0..6)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / 6.0;
            Point2::from_array([r * a.cos(), r * a.sin()])
        })
        .collect();
    Profile2d { points }
}

fn hex_extrude(plane: SketchPlane) -> Operation {
    Operation::Extrude {
        plane,
        profile: hexagon(10.0),
        dist: 5.0,
    }
}

#[test]
fn hexagon_extrude_exports_to_step_and_obj() {
    let model = replay(&[hex_extrude(SketchPlane::world_xy())]).unwrap();

    // The pipeline yields a valid closed b-rep.
    assert!(nacre_validate::validate(&model).is_empty());

    // ops model → a real STEP solid.
    let step = nacre_step::to_step(&model).expect("export to STEP");
    assert!(step.contains("MANIFOLD_SOLID_BREP"));

    // ops model → an OBJ mesh: 12 vertices (2n for a hexagon prism) plus faces.
    let obj = nacre_tess::to_obj(&model);
    let v_lines = obj.lines().filter(|l| l.starts_with("v ")).count();
    let f_lines = obj.lines().filter(|l| l.starts_with("f ")).count();
    assert_eq!(v_lines, 12);
    assert!(f_lines > 0);
}

#[test]
fn multi_op_log_composes_through_export() {
    // Two extrudes on parallel planes produce two independent solids that both
    // survive validation and export.
    let far = SketchPlane {
        origin: Point3::from_array([50.0, 0.0, 0.0]),
        ..SketchPlane::world_xy()
    };
    let model = replay(&[hex_extrude(SketchPlane::world_xy()), hex_extrude(far)]).unwrap();

    assert_eq!(model.solids.len(), 2);
    assert!(nacre_validate::validate(&model).is_empty());
    assert!(nacre_step::to_step(&model).is_ok());
}

// ---------------------------------------------------------------------------
// Mesh gate: watertight structure + area against an independent computation.
//
// These live here because `nacre-tess` cannot reach `extrude` (ops sits above it)
// and cannot reach `nacre-props` at all. The three checks catch different faults:
//
//   * `non_watertight` — an undirected triangle edge used other than twice. Blind
//     to a fan that escapes its polygon; that fan's edges pair up fine.
//   * `mesh_area` vs `props.area` — the *unsigned* sum. A triangle outside the
//     polygon is wound backwards, so a signed sum cancels and shoelace agrees with
//     itself. Unsigned, it overshoots.
//   * triangle count — pins how many the triangulator should emit.
//
// Several assertions below pin **today's bugs**, not today's contract. Each is
// inverted by the commit that fixes it, so the fix shows up in the diff rather
// than in a commit message (cf. `a_doubly_crossed_edge_breaks_the_transition_oracle`).
// ---------------------------------------------------------------------------

use nacre_ops::{BoolKind, OpOutput, apply, boolean};
use nacre_store::Handle;
use nacre_tess::{TessConfig, Tessellation, tessellate};
use nacre_topo::{Model, Solid};

/// Σ |triangle area| — unsigned on purpose.
fn mesh_area(t: &Tessellation) -> f64 {
    t.triangles
        .iter()
        .map(|(_, tri)| {
            let p = tri.vertices.map(|h| t.vertices.get(h).pos);
            0.5 * (p[1] - p[0]).cross(p[2] - p[0]).norm()
        })
        .sum()
}

/// Undirected triangle edges not shared by exactly two triangles.
fn non_watertight(t: &Tessellation) -> usize {
    let mut counts: std::collections::HashMap<(u32, u32), usize> = Default::default();
    for (_, tri) in t.triangles.iter() {
        let [a, b, c] = tri.vertices.map(|h| h.index());
        for (x, y) in [(a, b), (b, c), (c, a)] {
            *counts.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    counts.values().filter(|&&n| n != 2).count()
}

/// `(triangles, mesh area − props area, non-watertight edges)`.
fn mesh_vs_props(model: &Model, solid: Handle<Solid>) -> (usize, f64, usize) {
    let t = tessellate(model, &TessConfig::default());
    let props = nacre_props::mass_props(model, solid).expect("planar mass props");
    (
        t.triangles.len(),
        mesh_area(&t) - props.area,
        non_watertight(&t),
    )
}

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d {
        points: vec![p2(a, a), p2(b, a), p2(b, b), p2(a, b)],
    }
}

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

/// The `u_prism` of the boolean suite: prong tops at different heights so no two
/// faces are coplanar. Its cap is **not star-shaped from `(0,0)`**.
fn u_prism() -> (Model, Handle<Solid>) {
    let u = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(3.0, 0.0),
            p2(3.0, 2.3),
            p2(2.0, 2.3),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
    let m = replay(&[Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: u,
        dist: 1.0,
    }])
    .unwrap();
    let s = m.live_solids[0];
    (m, s)
}

/// The unit cube with a `0.4`-square pocket `0.5` deep in its top face.
fn pocketed_cube() -> (Model, Handle<Solid>) {
    let mut m = replay(&[Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: square(0.0, 1.0),
        dist: 1.0,
    }])
    .unwrap();
    let s = m.live_solids[0];
    let shell = m.solids.get(s).outer;
    let top = *m
        .shells
        .get(shell)
        .faces
        .iter()
        .find(|&&fh| {
            m.faces.get(fh).outer.half_edges.iter().all(|he| {
                let b = m.edges.get(he.edge).bounds.unwrap();
                b.iter()
                    .all(|&v| (m.vertices.get(v).point[2] - 1.0).abs() < 1e-9)
            })
        })
        .expect("the top face");
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: top,
            profile: square(0.3, 0.7),
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!("pocket yields PocketOnFace output")
    };
    (m, solid)
}

/// Two cubes fused across their shared face. `boolean` supersedes both operands —
/// it does not delete them, and `Store` is append-only by design.
fn stacked_fuse() -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

#[test]
fn a_convex_hole_free_solid_meshes_exactly() {
    // The shape the fan was written for. Nothing to catch here — this pins the gate
    // itself, so a later failure means the gate moved, not the mesh.
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let (tris, delta, leaks) = mesh_vs_props(&m, c);
    assert_eq!((tris, leaks), (12, 0));
    assert!(delta.abs() < 1e-12, "area delta {delta}");
}

#[test]
fn a_non_star_shaped_cap_overshoots_its_area_today() {
    // `triangulate_planar` fans from `ring[0]`. On the U's cap that drags a triangle
    // across the notch: it comes out clockwise, and the unsigned sum counts its 2.6
    // twice over instead of cancelling. Once per cap ⇒ +5.2.
    //
    // Watertight sees **nothing**: the fan's edges still pair up. Only area speaks.
    let (m, s) = u_prism();
    let (tris, delta, leaks) = mesh_vs_props(&m, s);
    assert_eq!((tris, leaks), (28, 0));
    assert!((delta - 5.2).abs() < 1e-9, "area delta {delta}");
}

#[test]
fn a_pocket_lid_is_filled_today() {
    // `triangulate_planar` ignores `face.inner`, so the lid's hole (0.16) is meshed
    // solid. Its rim edges are then used once instead of twice — the one bug that
    // watertight does catch. (The superseded top face used to ride along too, worth
    // another 1.0; the reachable walk retired that.)
    let (m, s) = pocketed_cube();
    let (tris, delta, leaks) = mesh_vs_props(&m, s);
    assert_eq!(tris, 22);
    // The four rim edges of the hole: used once by a pocket wall, never by the lid.
    assert_eq!(leaks, 4);
    assert!((delta - 0.16).abs() < 1e-9, "area delta {delta}");
}

#[test]
fn a_boolean_result_meshes_only_the_live_solid() {
    // `Store` is append-only; `boolean` supersedes rather than deletes. Walking the
    // face store meshed both original cubes (6 + 6) alongside the union.
    let (m, s) = stacked_fuse();
    let (tris, delta, leaks) = mesh_vs_props(&m, s);
    assert_eq!((tris, leaks), (20, 0));
    assert!(delta.abs() < 1e-12, "area delta {delta}");
}
