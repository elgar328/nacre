//! Cross-crate composition: ops (replay) → validate → step/tess export.
//!
//! The per-crate unit tests each cover one layer; this guards that the layers
//! compose. STEP format details (schema token, entity list, step_io round-trip)
//! stay owned by `nacre-step`'s own tests — here we only prove the pipeline hands
//! a valid model through to real STEP/OBJ output.

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{Operation, Profile2d, SketchPlane, replay};

/// A regular hexagon of the given radius, centred on the sketch origin.
fn hexagon(r: f64) -> Profile2d {
    let points = (0..6)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / 6.0;
            Point2::from_array([r * a.cos(), r * a.sin()])
        })
        .collect();
    Profile2d::polygon(points).unwrap()
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
    let obj = nacre_tess::to_obj(&model).expect("planar model meshes");
    let v_lines = obj.lines().filter(|l| l.starts_with("v ")).count();
    let f_lines = obj.lines().filter(|l| l.starts_with("f ")).count();
    assert_eq!(v_lines, 12);
    assert!(f_lines > 0);
}

#[test]
fn multi_op_log_composes_through_export() {
    // Two extrudes on parallel planes produce two independent solids that both
    // survive validation and export.
    let far = SketchPlane::world_xy().with_origin(Point3::from_array([50.0, 0.0, 0.0]));
    let model = replay(&[hex_extrude(SketchPlane::world_xy()), hex_extrude(far)]).unwrap();

    assert_eq!(model.solids.len(), 2);
    assert!(nacre_validate::validate(&model).is_empty());
    assert!(nacre_step::to_step(&model).is_ok());
}

// ---------------------------------------------------------------------------
// Mesh gate: watertight structure + area and volume against an independent
// computation.
//
// These live here because `nacre-tess` cannot reach `extrude` (ops sits above it)
// and cannot reach `nacre-props` at all. Four checks, four faults, none of them
// interchangeable:
//
//   * triangle count — pins how many the triangulator should emit.
//   * `non_watertight` — an undirected triangle edge used other than twice. Catches
//     a lost face or an unshared rim. Blind to a fan that escapes its polygon
//     (that fan's edges still pair up) and to a face wound backwards.
//   * `mesh_area` vs `props.area` — the *unsigned* sum. A triangle outside the
//     polygon is wound backwards, so a signed sum cancels and shoelace agrees with
//     itself. Unsigned, it overshoots. Blind to a face wound backwards.
//   * `mesh_volume` vs `props.volume` — the *signed* sum. This one is worth stating
//     precisely, because it looks circular and is not. `props` takes each face's
//     normal from `Face.orientation` and its area from `|A_vec|`; it never reads a
//     loop's winding. `tess` reads only the ring's Newell normal and never reads
//     `orientation`. So the two volumes agreeing says exactly:
//
//         on every face, `Face.orientation` agrees with its ring's winding.
//
//     Two different sources, compared. `validate` asks the same question
//     topologically (a shared edge used twice the same way); this asks it
//     geometrically, and neither stands in for the other.
//
// Several assertions below pin **today's bugs**, not today's contract. Each is
// inverted by the commit that fixes it, so the fix shows up in the diff rather
// than in a commit message (cf. `a_doubly_crossed_edge_breaks_the_transition_oracle`).
// ---------------------------------------------------------------------------

use nacre_ops::{BoolError, BoolKind, OpOutput, apply, boolean};
use nacre_store::Handle;
use nacre_tess::{TessConfig, Tessellation, tessellate};
use nacre_topo::{Face, Model, Shell, Solid};

/// Test shim: a boolean whose result is exactly one solid (cell 0.4 multi-solid).
fn boolean_one(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let solids = boolean(model, kind, a, b)?;
    assert_eq!(
        solids.len(),
        1,
        "boolean_one: expected one solid, got {}",
        solids.len()
    );
    Ok(solids[0])
}

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

/// Σ ⅙ p₀·(p₁ × p₂) — the divergence theorem on a closed triangle soup. Signed on
/// purpose: a face whose triangles wind the other way subtracts where it should add.
fn mesh_volume(t: &Tessellation) -> f64 {
    t.triangles
        .iter()
        .map(|(_, tri)| {
            let p = tri
                .vertices
                .map(|h| t.vertices.get(h).pos - Point3::origin());
            p[0].dot(p[1].cross(p[2])) / 6.0
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

/// What the gate measured. Named fields: four numbers read positionally is one
/// transposition away from a test that passes for the wrong reason.
struct Gate {
    tris: usize,
    /// mesh − props, unsigned area.
    area_delta: f64,
    /// mesh − props, signed volume.
    volume_delta: f64,
    leaks: usize,
}

fn mesh_vs_props(model: &Model, solid: Handle<Solid>) -> Gate {
    let t = tessellate(model, &TessConfig::default()).expect("planar model meshes");
    let props = nacre_props::mass_props(model, solid).expect("planar mass props");
    Gate {
        tris: t.triangles.len(),
        area_delta: mesh_area(&t) - props.area,
        volume_delta: mesh_volume(&t) - props.volume,
        leaks: non_watertight(&t),
    }
}

/// Both deltas are zero to a hair. Every fixture in this file must satisfy this;
/// the individual tests add their own triangle counts and hand values.
fn assert_agrees(g: &Gate, what: &str) {
    assert_eq!(g.leaks, 0, "{what}: non-watertight edges");
    assert!(
        g.area_delta.abs() < 1e-12,
        "{what}: area delta {}",
        g.area_delta
    );
    assert!(
        g.volume_delta.abs() < 1e-12,
        "{what}: volume delta {}",
        g.volume_delta
    );
}

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![p2(a, a), p2(b, a), p2(b, b), p2(a, b)]).unwrap()
}

/// `PocketOnFace`'s profile lives in a frame **derived from the face** (design §6), not in world
/// coordinates — but on the unit cube's lid the two coincide.
///
/// The sketch origin is the world origin projected onto the face's plane, and the arbitrary-axis
/// convention gives `n = ẑ` the axes `u = +x̂, v = +ŷ`, so a frame point `(a, b)` is world
/// `(a, b, 1)`. `half = 0.2` here is the `[0.3,0.7]²` void.
fn centred_on_the_cube_lid(half: f64) -> Profile2d {
    let (cx, cy) = (0.5, 0.5);
    Profile2d::polygon(vec![
        p2(cx - half, cy + half),
        p2(cx - half, cy - half),
        p2(cx + half, cy - half),
        p2(cx + half, cy + half),
    ])
    .unwrap()
}

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

/// The `u_prism` of the boolean suite: prong tops at different heights so no two
/// faces are coplanar. Its cap is **not star-shaped from `(0,0)`**.
fn u_prism() -> (Model, Handle<Solid>) {
    let u = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(3.0, 0.0),
        p2(3.0, 2.3),
        p2(2.0, 2.3),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
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
            profile: centred_on_the_cube_lid(0.2),
            dist: 0.5,
        },
    )
    .unwrap() else {
        unreachable!("pocket yields PocketOnFace output")
    };
    (m, solid)
}

/// The L-prism with a stub standing in its top face, footprint strictly inside it —
/// the `l_and_dimple` of the boolean suite. `Cut` leaves a blind pocket, `Fuse` a
/// boss; either way the L's top face gains an inner loop, and that loop's rim edges
/// must pair with the walls that drop or rise from them.
fn l_and_dimple(kind: BoolKind) -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: l,
        dist: 1.0,
    }])
    .unwrap();
    let a = m.live_solids[0];
    let b = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, kind, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// `Cut(stub, L)` — the operands of `l_and_dimple` reversed. The L's top face keeps only
/// the stub's footprint, so the answer is a `0.4 × 0.4 × 0.5` box whose floor is that
/// island face: an outer loop made of nothing but seam vertices (cell 3f-2).
fn island_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: l,
        dist: 1.0,
    }])
    .unwrap();
    let a = m.live_solids[0];
    let b = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, b, a).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// The L-prism and an L-shaped bar in its notch, biting two corners of the cap. Two
/// chords on one face; under `Cut` the bar's floor splits into the two bites' floors.
fn notch_bar_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let bar = Profile2d::polygon(vec![
        p2(1.8, 0.8),
        p2(2.1, 0.8),
        p2(2.1, 2.1),
        p2(0.8, 2.1),
        p2(0.8, 1.8),
        p2(1.8, 1.8),
    ])
    .unwrap();
    let mut m = replay(&[
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: l,
            dist: 1.0,
        },
        Operation::Extrude {
            plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
            profile: bar,
            dist: 1.0,
        },
    ])
    .unwrap();
    let (a, b) = (m.live_solids[0], m.live_solids[1]);
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// The L-prism with an L-shaped stub standing wholly inside its cap: a blind pocket whose
/// lid carries a **non-convex** inner loop, the first the bridging triangulator has seen.
fn ell_dimple_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let ell = Profile2d::polygon(vec![
        p2(0.2, 0.25),
        p2(0.85, 0.25),
        p2(0.85, 0.4),
        p2(0.35, 0.4),
        p2(0.35, 0.9),
        p2(0.2, 0.9),
    ])
    .unwrap();
    let mut m = replay(&[
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: l,
            dist: 1.0,
        },
        Operation::Extrude {
            plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
            profile: ell,
            dist: 1.0,
        },
    ])
    .unwrap();
    let (a, b) = (m.live_solids[0], m.live_solids[1]);
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// `Cut(u, slab)` — the slab shears off both prong tops. Its `y = 1.5` face has its whole
/// boundary dropped, so its two loops are **islands**: one input face, two output faces,
/// both flipped, and they become the result's new end caps (cell 3f-3).
fn u_cut_by_slab() -> (Model, Handle<Solid>) {
    let (mut m, u) = u_prism();
    let slab = m.add_cuboid(
        Point3::from_array([-0.5, 1.5, -0.5]),
        Point3::from_array([3.5, 2.5, 1.5]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// `Cut(staple, L)` — a П drawn in the XZ plane, extruded along `−y`, minus the L-prism.
/// The L's cap contributes two faces to the result: the kept region around the reflex corner,
/// and an **island** where the near leg's footprint sits in a dropped region (cell 3f-4).
fn staple_cut_by_l() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let staple = Profile2d::polygon(vec![
        p2(0.1, 0.5),
        p2(0.6, 0.5),
        p2(0.6, 1.3),
        p2(0.8, 1.3),
        p2(0.8, 0.45),
        p2(1.4, 0.45),
        p2(1.4, 1.5),
        p2(0.1, 1.5),
    ])
    .unwrap();
    let mut m = replay(&[
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: l,
            dist: 1.0,
        },
        Operation::Extrude {
            plane: SketchPlane::from_origin_normal(
                Point3::from_array([0.0, 1.3, 0.0]),
                Vector3::from_array([0.0, -1.0, 0.0]),
            )
            .expect("a unit normal"),
            profile: staple,
            dist: 0.65,
        },
    ])
    .unwrap();
    let (a, b) = (m.live_solids[0], m.live_solids[1]);
    let r = boolean_one(&mut m, BoolKind::Cut, b, a).unwrap();
    m.rebuild_adjacency();
    (m, r)
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
    let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

#[test]
fn a_convex_hole_free_solid_meshes_exactly() {
    // The shape the fan was written for. Nothing to catch here — this pins the gate
    // itself, so a later failure means the gate moved, not the mesh.
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let g = mesh_vs_props(&m, c);
    assert_eq!(g.tris, 12);
    assert_agrees(&g, "cube");
}

#[test]
fn a_non_star_shaped_cap_meshes_exactly() {
    // The fan used to drag a triangle across the U's notch: clockwise, so the
    // unsigned sum counted its 2.6 forwards instead of backwards, once per cap ⇒
    // +5.2. Watertight saw nothing — the fan's edges still paired up. Only area
    // spoke, which is why both checks are here.
    //
    // Ear clipping owes nothing to `ring[0]`. Each cap is `8 − 2 = 6` triangles.
    let (m, s) = u_prism();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 28);
    assert_agrees(&g, "u prism");
}

#[test]
fn a_pocket_lid_carries_its_hole() {
    // The lid used to be meshed solid (Δ +0.16), and its four rim edges were then
    // used once by a pocket wall and never by the lid — the one bug watertight
    // caught. Bridged and ear-clipped, the lid is `4 + 4 + 2·1 − 2 = 8` triangles;
    // the solid comes to 28.
    let (m, s) = pocketed_cube();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 28);
    assert_agrees(&g, "pocketed cube");
}

#[test]
fn the_bootstrap_obj_carries_holes_too() {
    // The OBJ a user actually looks at comes from `to_obj(&Model)`, not from
    // `tessellate`. It was a second fan with the same two bugs. Both now share one
    // triangulator, so the lid's hole survives to Quick Look: `f` lines match the
    // tessellation's triangles, and no triangle references a superseded vertex.
    let (m, s) = pocketed_cube();
    let obj = nacre_tess::to_obj(&m).expect("planar model meshes");
    let f_lines = obj.lines().filter(|l| l.starts_with("f ")).count();
    assert_eq!(f_lines, mesh_vs_props(&m, s).tris);
}

#[test]
fn a_boolean_result_carries_its_hole() {
    // The gate's first run on a boolean result. Cell 3f-1 gives the L's top face an
    // inner loop, and the pocket's walls (or the boss's) hang off that loop's rim
    // edges — so watertight is what proves the rim is shared rather than duplicated.
    //
    // 36 triangles: the holed lid is `6 + 4 + 2·1 − 2 = 10`, the bottom cap `6 − 2 = 4`,
    // and 2 apiece for the L's six sides, the four walls and the far cap. `Cut` and
    // `Fuse` only differ in which way the walls point, so both counts and both areas
    // agree — the same 14.8 = `(8·1 + 2·3) − 0.16 + 4·0.4·0.5 + 0.16`, though the
    // solids differ (2.92 against 3.08, pinned in the boolean suite).
    //
    // Meshing the lid solid would overshoot the area by the hole's 0.16 while staying
    // watertight; that is the fault this catches, and the fan had it until the tess cell.
    for kind in [BoolKind::Cut, BoolKind::Fuse] {
        let (m, s) = l_and_dimple(kind);
        let g = mesh_vs_props(&m, s);
        assert_eq!(g.tris, 36, "{kind:?}");
        assert_agrees(&g, &format!("{kind:?} dimple"));
        let area = nacre_props::mass_props(&m, s).unwrap().area;
        assert!((area - 14.8).abs() < 1e-9, "{kind:?} area {area}");
    }
}

#[test]
fn a_sealed_cavity_meshes_watertight() {
    // Cell 5c: fusing a slab across the pocket mouth seals it into an enclosed void, so the
    // result has two shells. The gate proves both close (watertight, no leaks between the
    // outer surface and the void's) and that the void's inward faces subtract in the mesh
    // volume exactly as they do in `props` — the cavity's `2.408` against the outer `2.44`.
    let (mut m, pc) = pocketed_cube();
    let slab = m.add_cuboid(
        Point3::from_array([-0.2, -0.25, 0.7]),
        Point3::from_array([1.3, 1.2, 1.5]),
    );
    let s = boolean_one(&mut m, BoolKind::Fuse, slab, pc).unwrap();
    assert_eq!(m.solids.get(s).cavities.len(), 1);
    let g = mesh_vs_props(&m, s);
    assert_agrees(&g, "sealed cavity");
    assert!((nacre_props::mass_props(&m, s).unwrap().volume - 2.408).abs() < 1e-9);
}

#[test]
fn a_split_face_meshes_like_two() {
    // Cell 3e-2's first shape whose bar floor becomes *two* b-rep faces. The signed volume
    // is the only check here that looks at their orientation: each floor got its own
    // `Orientation` from `flip` and its own ring from the stitcher, and nothing else in
    // this file compares those two sources.
    //
    // The bites are corner cuts, so `props.area` is the L's own 14.0 — the three faces
    // removed and the three added cancel. A gate that only summed area would see nothing.
    let (m, s) = notch_bar_cut();
    // 44: the bitten cap is a 10-gon (8), the bottom cap a hexagon (4), two untouched
    // walls (2 each), four bitten walls that are hexagons (4 each), then the bar's four
    // side pieces (2 each) and its two floors (2 each).
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 44);
    assert_agrees(&g, "notch bar cut");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!((props.area - 14.0).abs() < 1e-9, "area {}", props.area);
    assert!(
        (props.volume - 2.96).abs() < 1e-9,
        "volume {}",
        props.volume
    );
}

#[test]
fn a_non_convex_hole_bridges_and_meshes() {
    // Every inner loop the triangulator has met was a rectangle. This one has a reflex
    // node, so `bridge_holes` must find a mutually-visible pair across a ring that is not
    // star-shaped from anywhere obvious, and ear clipping must survive the slit.
    //
    // 44 triangles: the bitten lid is `6 + 6 + 2·1 − 2 = 12`, the bottom cap `6 − 2 = 4`,
    // 2 apiece for the L's six walls and the pocket's six, and `6 − 2 = 4` for the floor.
    let (m, s) = ell_dimple_cut();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 44);
    assert_agrees(&g, "ell dimple cut");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!((props.area - 15.3).abs() < 1e-9, "area {}", props.area);
    assert!(
        (props.volume - (3.0 - 0.1725 * 0.5)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
}

#[test]
fn two_islands_from_one_face_mesh_like_two() {
    // Cell 3f-3's payoff on the gate. The slab's face yields two islands, and each took its
    // `Orientation` from `flip` and its ring from `orient_seam_loop`. The signed volume is
    // the only check here that compares those two sources; watertight and the unsigned area
    // would wave a reversed island through.
    //
    // 28 triangles: the clipped U is an octagonal prism, so `6 + 6` for its caps and 2 for
    // each of eight walls — **two of those walls are the islands**, not extras. Area 18,
    // volume 4.0.
    let (m, s) = u_cut_by_slab();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 28);
    assert_agrees(&g, "u cut by slab");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!((props.area - 18.0).abs() < 1e-9, "area {}", props.area);
    assert!((props.volume - 4.0).abs() < 1e-9, "volume {}", props.volume);
}

#[test]
fn an_island_beside_an_arc_meshes() {
    // Cell 3f-4 on the gate. One input face — the L's cap — gives the result a kept region
    // *and* an island, and the two were told apart by containment alone. Both are flipped
    // (`Cut`'s B-piece), and the signed mesh volume is the only check in this file that
    // compares each face's `Orientation` against its own ring.
    let (m, s) = staple_cut_by_l();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 40);
    assert_agrees(&g, "staple cut by l");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!(
        (props.volume - (0.7605 - 0.311)).abs() < 1e-9,
        "volume {}",
        props.volume
    );
    assert!((props.area - 4.68).abs() < 1e-9, "area {}", props.area);
}

#[test]
fn only_the_signed_volume_sees_a_reversed_face() {
    // The check earns its place by being the only one in the gate that speaks.
    //
    // Reverse one face's outer loop and nothing else. The mesh keeps its 12 triangles;
    // every undirected edge is still used twice; every triangle keeps its area, so the
    // unsigned sum is unmoved. `props` is unmoved too — it reads `Face.orientation` for
    // the normal and `|A_vec|` for the area, and neither changed. Only the signed
    // volume, which reads the ring, dissents.
    //
    // `validate` also speaks (`NonOpposedEdge`), topologically, and it is not what is on
    // trial here. The point is that the *mesh* gate would have waved this through.
    //
    // The cuboid is `[1,2]³`, not `[0,1]³`: three faces of the latter pass through the
    // origin and contribute nothing to `Σ ⅙ p₀·(p₁ × p₂)`, so reversing one of those
    // would move the signed volume by exactly zero and this test would pass vacuously.
    let mut m = Model::new();
    let c = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    let good = mesh_vs_props(&m, c);
    assert_eq!(good.tris, 12);
    assert_agrees(&good, "offset cube");

    // `Store` is append-only: push a replacement face, swap it into a fresh shell and
    // solid, and move the live handle. The original falls out of `reachable()`.
    let faces = m.shells.get(m.solids.get(c).outer).faces.clone();
    let victim = faces[0];
    let f = m.faces.get(victim).clone();
    let mut outer = f.outer.clone();
    outer.half_edges.reverse();
    for he in &mut outer.half_edges {
        he.forward = !he.forward;
    }
    let bad = m.faces.push(Face { outer, ..f });
    let swapped = faces
        .iter()
        .map(|&x| if x == victim { bad } else { x })
        .collect();
    let shell = m.shells.push(Shell { faces: swapped });
    let solid = m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    m.live_solids = vec![solid];
    m.rebuild_adjacency();

    let g = mesh_vs_props(&m, solid);
    assert_eq!((g.tris, g.leaks), (12, 0), "watertight stayed quiet");
    assert!(
        g.area_delta.abs() < 1e-12,
        "area stayed quiet: {}",
        g.area_delta
    );
    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((vol - 1.0).abs() < 1e-12, "props stayed quiet: {vol}");
    // Exactly ⅔: `faces[0]` is the `x = 1` face, whose flux `∮ r·n̂ dA` is `−1`, so it
    // contributed `−⅓` and now contributes `+⅓`. A different face would give a
    // different number — if `add_cuboid` ever reorders, this says so rather than
    // shrugging at a loose bound.
    assert!(
        (g.volume_delta - 2.0 / 3.0).abs() < 1e-12,
        "the signed volume should have moved by twice the face's flux, got {}",
        g.volume_delta
    );
    assert!(
        nacre_validate::validate(&m)
            .iter()
            .any(|v| matches!(v, nacre_validate::Violation::NonOpposedEdge { .. })),
        "and validate speaks too, of its own accord"
    );
}

#[test]
fn an_island_face_meshes_like_any_other() {
    // The gate's first run on a face whose outer loop is *all* seam. It is a plain box —
    // 12 triangles, `2·0.16 + 4·0.4·0.5 = 1.12` of surface, 0.08 of volume — and that is
    // the claim: cell 3f-2's island is not a special kind of face downstream.
    //
    // The signed volume is the one check with something new to say. `props` reads
    // `Face.orientation`, which `assemble_fuse_cut` set from the `flip` flag; `tess`
    // reads the ring, which `orient_seam_loop` wound and `flip` then reversed. Two
    // routes to the same face's outward direction, and they have to meet.
    let (m, s) = island_cut();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 12);
    assert_agrees(&g, "island cut");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!((props.area - 1.12).abs() < 1e-9, "area {}", props.area);
    assert!(
        (props.volume - 0.08).abs() < 1e-9,
        "volume {}",
        props.volume
    );
}

#[test]
fn a_flipped_island_loop_is_caught() {
    // Cell 3f-1 showed that a flipped *hole* is caught by `validate` (topologically) and
    // by `tessellate` (geometrically, `HoleWinding`), while volume and OCCT see nothing.
    // The island inverts that, and the inversion was measured rather than assumed.
    //
    // `tessellate` never looks at an outer ring's winding — it takes its normal *from*
    // the ring — so it returns `Ok` and its triangles simply face inward. `props` never
    // looks at a ring at all, so its volume is still exactly 0.08. Watertight and the
    // unsigned area are blind by construction. `validate` speaks, and until this cell it
    // was the only thing that did.
    //
    // The signed volume moves by `2 · 0.16 / 3`: the island is the `z = 1` floor, area
    // 0.16, outward normal `−z`, flux `−0.16`, so its contribution goes from `−0.0533`
    // to `+0.0533`.
    let (mut m, s) = island_cut();

    // The island is the face whose whole boundary *is* the cut's seam ring, so it is the
    // one face lying on the other operand's plane — the L's top, `z = 1`.
    //
    // It used to be selected as "the only face all of whose vertices are `Discovered`".
    // That died with the arrangement: every result vertex is now named by a plane triple
    // (`an_unrotated_boolean_names_every_vertex_by_its_plane_triple`), so all six faces
    // match and `find` silently took the `z = 1.5` top instead — five of the six
    // assertions below pass on any face, so only the signed volume noticed. Hence the
    // geometric predicate, and hence the count: the uniqueness this relies on is asserted,
    // not narrated.
    let faces = m.shells.get(m.solids.get(s).outer).faces.clone();
    let isles: Vec<_> = faces
        .iter()
        .copied()
        .filter(|&f| {
            m.faces.get(f).outer.half_edges.iter().all(|he| {
                let b = m.edges.get(he.edge).bounds.unwrap();
                b.iter()
                    .all(|&v| (m.vertices.get(v).point[2] - 1.0).abs() < 1e-9)
            })
        })
        .collect();
    assert_eq!(isles.len(), 1, "the island face is not unique: {isles:?}");
    let isle = isles[0];

    // `Store` is append-only: push the reversed face, swap it into a fresh shell and
    // solid, move the live handle.
    let f = m.faces.get(isle).clone();
    let mut outer = f.outer.clone();
    outer.half_edges.reverse();
    for he in &mut outer.half_edges {
        he.forward = !he.forward;
    }
    let bad = m.faces.push(Face { outer, ..f });
    let swapped = faces
        .iter()
        .map(|&x| if x == isle { bad } else { x })
        .collect();
    let shell = m.shells.push(Shell { faces: swapped });
    let solid = m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    m.live_solids = vec![solid];
    m.rebuild_adjacency();

    assert!(
        nacre_validate::validate(&m)
            .iter()
            .any(|v| matches!(v, nacre_validate::Violation::NonOpposedEdge { .. })),
        "validate stayed quiet"
    );
    assert!(
        tessellate(&m, &TessConfig::default()).is_ok(),
        "tessellate has no opinion about an outer ring's winding"
    );
    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((vol - 0.08).abs() < 1e-12, "props stayed quiet: {vol}");

    let g = mesh_vs_props(&m, solid);
    assert_eq!((g.tris, g.leaks), (12, 0), "watertight stayed quiet");
    assert!(g.area_delta.abs() < 1e-12, "area stayed quiet");
    assert!(
        (g.volume_delta - 2.0 * 0.16 / 3.0).abs() < 1e-12,
        "signed volume delta {}",
        g.volume_delta
    );
}

#[test]
fn a_boolean_result_meshes_only_the_live_solid() {
    // `Store` is append-only; `boolean` supersedes rather than deletes. Walking the
    // face store meshed both original cubes (6 + 6) alongside the union. The union is a
    // clean 6-face 1×1×2 box (cell fuse-coplanar-merge dissolves the interface corners),
    // so 12 triangles — not the superseded operands' faces.
    let (m, s) = stacked_fuse();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 12);
    assert_agrees(&g, "stacked fuse");
}

/// The gate on a holed *operand*. Cell 3f-5 lets a pocketed cube be cut, and the lid
/// keeps its hole either by riding through `whole()` untouched or by being placed inside
/// the region a seam arc leaves. Both roads end here.
///
/// `Σ|area| == props.area` catches a hole filled in or a fan spilling out; `mesh_volume ==
/// props.volume` is the only one of the four that can see a hole ring reversed, since
/// `props` reads the face's `orientation` and never the ring's winding.
///
/// The bite is the symmetric `[0.85,1.15]³` cell (5a) gave back: its vertical edge pierces
/// the lid at `(0.85, 0.85)`, a point on the lid's fan diagonal from every apex, so the
/// boolean's exact crossing test is what carries it here. The lid that comes out is holed
/// *and* notched, and `nacre-tess` clips ears rather than fanning, so the gate is a second
/// opinion from machinery that never shared the fan's blind spot.
#[test]
fn a_holed_operand_keeps_its_hole_through_a_cut() {
    for (name, box_lo, box_hi) in [
        ("corner", [0.85, 0.85, 0.85], [1.15, 1.15, 1.15]),
        ("bottom corner", [0.85, 0.85, -0.15], [1.15, 1.15, 0.15]),
    ] {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(Point3::from_array(box_lo), Point3::from_array(box_hi));
        let s = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let g = mesh_vs_props(&m, s);
        assert_agrees(&g, name);
        let vol = nacre_props::mass_props(&m, s).unwrap().volume;
        assert!((vol - 0.916625).abs() < 1e-9, "{name} volume {vol}");
    }
}

/// A `Cut`'s inside-B piece with a hole in it: the pocket's lid, kept, reversed, hole and
/// all. `flip` and `inner` had never met. Only the signed volume tells a hole reversed
/// with its face from a hole reversed alone, and it is the one check `props` cannot make
/// for itself.
#[test]
fn a_flipped_face_keeps_its_hole() {
    let (mut m, pc) = pocketed_cube();
    let slab = m.add_cuboid(
        Point3::from_array([-0.2, -0.25, 0.3]),
        Point3::from_array([1.3, 1.2, 1.5]),
    );
    let s = boolean_one(&mut m, BoolKind::Cut, slab, pc).unwrap();
    m.rebuild_adjacency();
    let g = mesh_vs_props(&m, s);
    assert_agrees(&g, "slab cut by pocket");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!(
        (props.volume - 1.99).abs() < 1e-9,
        "volume {}",
        props.volume
    );
}

/// The gate on a solid of genus 1 — the first this kernel makes. A rod drilled clean
/// through the L's bar leaves each cap holed and the four tunnel walls joining the two
/// rims, so `watertight` is what proves the rims are shared rather than duplicated, and
/// `Σ|area| == props.area` is what proves neither hole was filled in.
///
/// `mesh_volume == props.volume` has the sharpest teeth here. A tunnel is a hole ring on
/// one cap and a hole ring on the other, wound oppositely about their own faces' normals;
/// reverse either and props, reading `orientation` and `|A_vec|`, would not notice.
#[test]
fn a_drilled_solid_meshes_watertight() {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let mut m = replay(&[Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile: l,
        dist: 1.0,
    }])
    .unwrap();
    let a = m.live_solids[0];
    let rod = m.add_cuboid(
        Point3::from_array([0.3, 0.3, -0.5]),
        Point3::from_array([0.5, 0.6, 1.5]),
    );
    let s = boolean_one(&mut m, BoolKind::Cut, a, rod).unwrap();
    m.rebuild_adjacency();

    let g = mesh_vs_props(&m, s);
    assert_agrees(&g, "drilled L");
    let props = nacre_props::mass_props(&m, s).unwrap();
    assert!(
        (props.volume - 2.94).abs() < 1e-9,
        "volume {}",
        props.volume
    );
}

/// A notch bitten from one edge of a cube: cell (5b) opened it (the convex path rejected
/// it as `poke_through`). The kept region's boundary is a seam chord plus part of one
/// edge — a run with no vertex — and that shape is new to the mesh gate, which had only
/// ever fanned squares-with-square-holes and reflex caps. `10³ − 4·1.4·1.2`.
#[test]
fn a_notched_cube_meshes() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let notch = m.add_cuboid(
        Point3::from_array([3.0, -1.0, -1.0]),
        Point3::from_array([7.0, 1.4, 1.2]),
    );
    let s = boolean_one(&mut m, BoolKind::Cut, a, notch).unwrap();
    m.rebuild_adjacency();
    let g = mesh_vs_props(&m, s);
    assert_agrees(&g, "notched cube");
    let vol = nacre_props::mass_props(&m, s).unwrap().volume;
    assert!((vol - 993.28).abs() < 1e-9, "volume {vol}");
}

/// **A wall with two windows — the case the old triangulator could not mesh.**
///
/// A hub box and two bars crossing its `x = −1` face, so that face carries two inner
/// loops. The b-rep was always fine here (`validate` clean, analytic mass properties);
/// only the mesh failed, and it failed by stalling.
///
/// The cause was a ring the triangulator built for itself: each hole was merged into
/// the outer ring by a zero-width bridge to the first mutually visible vertex, and both
/// holes picked the *same* corner, so the ring visited it three times. A repeated vertex
/// makes the ring non-simple, and Meisters' two-ears theorem — the only reason ear
/// clipping terminates — does not apply. Measured: all six convex vertices had a
/// diagonal that genuinely crossed the other hole.
///
/// **Nothing here is rotated.** The model that started the cell had seven fins at
/// 360/7°, which made it look like a precision problem; three boxes and two fuses are
/// enough. The monotone sweep never merges the rings, so the pinch is never built, and
/// the face now meshes into the same watertight surface everything else does.
#[test]
fn a_wall_with_two_windows_meshes() {
    let mut m = Model::new();
    let hub = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    m.rebuild_adjacency();
    let upper = m.add_cuboid(
        Point3::from_array([-4.0, 0.26, 1.0]),
        Point3::from_array([0.0, 0.70, 2.0]),
    );
    m.rebuild_adjacency();
    let part = boolean_one(&mut m, BoolKind::Fuse, hub, upper).unwrap();
    m.rebuild_adjacency();
    let lower = m.add_cuboid(
        Point3::from_array([-4.0, -0.70, 1.0]),
        Point3::from_array([0.0, -0.26, 2.0]),
    );
    m.rebuild_adjacency();
    let part = boolean_one(&mut m, BoolKind::Fuse, part, lower).unwrap();
    m.rebuild_adjacency();

    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the b-rep is correct"
    );
    let g = mesh_vs_props(&m, part);
    assert_agrees(&g, "wall with two windows");
}

/// **The model that started the cell**, meshed end to end.
///
/// Seven fins arrayed at 360/7° around a hub, then a square bore crossed with a 45°
/// copy — a code-CAD playground script that died on `tessellate` while `validate` and
/// `mass_props` were both happy. Two of the fins cross the same hub wall, which is how
/// that face ended up with the two holes whose bridges collided.
///
/// It is here for the shape, not the arithmetic: the two-window fixture above is the
/// minimal reproduction and the one that explains the bug. This one proves the fix
/// survives 22 booleans, rotations by an angle with no exact `f64`, and a face set built
/// entirely from `Discovered` vertices.
#[test]
fn the_fin_array_with_a_star_bore_meshes() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
    let turn = |m: &mut Model, s, deg: Rat| {
        let out = nacre_ops::apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(deg).expect("angle"),
                }),
            },
        )
        .expect("rotate");
        m.rebuild_adjacency();
        match out {
            nacre_ops::OpOutput::Transform { solid } => solid,
            o => panic!("{o:?}"),
        }
    };

    let fins = 7;
    let mut m = Model::new();
    let mut part = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    m.rebuild_adjacency();
    for i in 0..fins {
        // The script's `i * 360.0 / fins`, lifted to the exact rational it landed on.
        let deg = Rat::try_from_f64(i as f64 * 360.0 / fins as f64).expect("degrees");
        let fin = m.add_cuboid(
            Point3::from_array([0.5, -0.2, 1.0]),
            Point3::from_array([4.0, 0.2, 2.0]),
        );
        m.rebuild_adjacency();
        let fin = turn(&mut m, fin, deg);
        part = boolean_one(&mut m, BoolKind::Fuse, part, fin).unwrap();
        m.rebuild_adjacency();
    }
    for deg in [0, 45] {
        let bore = m.add_cuboid(
            Point3::from_array([-0.6, -0.6, 0.0]),
            Point3::from_array([0.6, 0.6, 2.0]),
        );
        m.rebuild_adjacency();
        let bore = turn(&mut m, bore, Rat::from_int(deg));
        part = boolean_one(&mut m, BoolKind::Cut, part, bore).unwrap();
        m.rebuild_adjacency();
    }

    assert!(nacre_validate::validate(&m).is_empty());
    let g = mesh_vs_props(&m, part);
    assert_agrees(&g, "fin array with a star bore");
}

/// **★ A dimension split into two lands exactly where the undivided one does.**
///
/// Every number below is a literal a user typed, and the only arithmetic is the kernel's own.
/// It used to be done in f64, where a prism raised from `z = 1.1` by `6.6` puts its top at
/// `7.699999999999999` — one ULP below the `7.7` that the block beside it reached in a single
/// step. Two planes where the model has one.
///
/// **That failure did not split anything and was not wrong**, which is what made it worth
/// fixing: the volume came out exactly right and the fuse yielded one body — a body carrying a
/// face of area `8.9e-16` and two faces more than the shape has. Nothing rejected, nothing
/// reported, and every later boolean, STEP export and mesh carried the sliver along.
///
/// Now the placement and the sweep are done in rationals read from the decimals as written
/// (`11/10 + 66/10 = 77/10`), so both paths reach the same plane and the fuse has one plane to
/// merge. Six faces, no sliver.
#[test]
fn a_split_dimension_meets_the_undivided_one() {
    fn rect(x0: f64, x1: f64) -> Profile2d {
        Profile2d::polygon(vec![p2(x0, 0.0), p2(x1, 0.0), p2(x1, 1.0), p2(x0, 1.0)]).unwrap()
    }
    fn raise(m: &mut Model, z: f64, x0: f64, x1: f64, dist: f64) -> Handle<Solid> {
        let out = nacre_ops::apply(
            m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
                profile: rect(x0, x1),
                dist,
            },
        )
        .expect("extrude");
        m.rebuild_adjacency();
        match out {
            nacre_ops::OpOutput::Extrude { solid, .. } => solid,
            o => panic!("{o:?}"),
        }
    }

    let mut m = Model::new();
    // The stack: 1.1, then 6.6 starting at 1.1 — the path that used to miss 7.7 by an ULP.
    let lower = raise(&mut m, 0.0, 3.0, 4.0, 1.1);
    let upper = raise(&mut m, 1.1, 3.0, 4.0, 6.6);
    let stack = boolean_one(&mut m, BoolKind::Fuse, lower, upper).unwrap();
    m.rebuild_adjacency();
    // The neighbour: the same height in one step, sharing the wall x = 4.
    let side = raise(&mut m, 0.0, 4.0, 5.0, 7.7);
    let part = boolean_one(&mut m, BoolKind::Fuse, stack, side).unwrap();
    m.rebuild_adjacency();

    assert!(
        nacre_validate::validate(&m).is_empty(),
        "the b-rep is valid"
    );
    let props = nacre_props::mass_props(&m, part).unwrap();
    assert!(
        (props.volume - 15.4).abs() < 1e-9,
        "volume {}",
        props.volume
    );

    let shell = m.solids.get(part).outer;
    let faces = &m.shells.get(shell).faces;
    let smallest = faces
        .iter()
        .filter_map(|&fh| nacre_props::face_props(&m, fh).ok())
        .map(|p| p.area)
        .fold(f64::INFINITY, f64::min);
    assert_eq!(faces.len(), 6, "the shape has six faces");
    assert!(
        smallest > 0.9,
        "no sliver: the smallest face is a whole wall, got {smallest}"
    );
}

/// **The exact-dimension guarantee reaches `pad`/`pocket` too**, and that is a claim about
/// `face_frame` — whose axes come from `any_perpendicular`, not from the caller — so it has
/// to be measured on that path rather than inferred from the extrude one.
///
/// **★ The obvious fixture proves nothing.** A base of `2.0` padded `1.1` then `6.6`, against
/// one padded `7.7`, agrees *already* in f64: `(2.0 + 1.1) + 6.6 == 2.0 + 7.7`. That test
/// passed on the unfixed kernel. The numbers below are chosen because they do not.
#[test]
fn a_pad_split_in_two_reaches_the_plane_the_whole_one_does() {
    fn base(m: &mut Model) -> Handle<Solid> {
        let OpOutput::Extrude { solid, .. } = nacre_ops::apply(
            m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: square(10.0),
                dist: 1.0,
            },
        )
        .expect("extrude") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid
    }
    fn square(n: f64) -> Profile2d {
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(n, 0.0), p2(n, n), p2(0.0, n)]).unwrap()
    }
    fn top_face(m: &Model, s: Handle<Solid>) -> Handle<Face> {
        let shell = m.solids.get(s).outer;
        let up = |f: Handle<Face>| nacre_props::face_props(m, f).unwrap();
        *m.shells
            .get(shell)
            .faces
            .iter()
            .filter(|&&f| up(f).normal.map(|n| n[2] > 0.5).unwrap_or(false))
            .max_by(|&&a, &&b| up(a).centroid[2].total_cmp(&up(b).centroid[2]))
            .expect("an upward face")
    }
    /// A 4×4 boss on the 10×10 lid — world `[3, 7]²`, which on a lid is also its frame
    /// coordinates (origin at the world origin's projection, `u = +x̂`, `v = +ŷ`).
    fn boss() -> Profile2d {
        Profile2d::polygon(vec![p2(3.0, 7.0), p2(3.0, 3.0), p2(7.0, 3.0), p2(7.0, 7.0)]).unwrap()
    }
    fn pad(m: &mut Model, s: Handle<Solid>, dist: f64) -> Handle<Solid> {
        let face = top_face(m, s);
        let OpOutput::PadOnFace { solid, .. } = nacre_ops::apply(
            m,
            &Operation::PadOnFace {
                face,
                profile: boss(),
                dist,
            },
        )
        .expect("pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        solid
    }

    assert_ne!(
        (1.0f64 + 0.1) + 0.1,
        1.0f64 + 0.2,
        "the fixture must discriminate, or it passes on the unfixed kernel too"
    );

    let (mut m1, mut m2) = (Model::new(), Model::new());
    let whole = base(&mut m1);
    let whole = pad(&mut m1, whole, 0.2);
    let split = base(&mut m2);
    let split = pad(&mut m2, split, 0.1);
    let split = pad(&mut m2, split, 0.1);

    let za = nacre_ops::face_plane(&m1, top_face(&m1, whole))
        .unwrap()
        .origin()[2];
    let zb = nacre_ops::face_plane(&m2, top_face(&m2, split))
        .unwrap()
        .origin()[2];
    assert_eq!(za, 1.2, "1.0 + 0.2");
    assert_eq!(za, zb, "the two pad paths reach different planes");
}

/// **A degenerate frame can no longer be handed in at all** — the refusal moved from the prism
/// builder to the door.
///
/// It used to be built as three raw axes and rejected inside `Plane::from_point_normal`, one
/// layer into `build_prism`. `SketchPlane`'s fields are private now and every public constructor
/// states a *plane*, so the degeneracies a caller can express are the ones checked here: a normal
/// with no direction, and three points that do not span one. The inner guard still stands for the
/// crate's own `from_axes`, but nothing outside can reach it.
#[test]
fn a_degenerate_frame_cannot_be_built() {
    let o = Point3::origin();
    assert!(
        SketchPlane::from_origin_normal(o, Vector3::from_array([0.0, 0.0, 0.0])).is_none(),
        "a zero normal is not a plane"
    );
    // `x_point` on the origin: no `+u` direction.
    assert!(SketchPlane::through_points(o, o, Point3::from_array([0.0, 1.0, 0.0])).is_none());
    // `y_hint` collinear with `+u`: the three points span a line, not a plane.
    assert!(
        SketchPlane::through_points(
            o,
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 0.0, 0.0]),
        )
        .is_none()
    );
}
