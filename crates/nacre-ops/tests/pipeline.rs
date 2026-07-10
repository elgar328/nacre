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

use nacre_ops::{BoolKind, OpOutput, apply, boolean};
use nacre_store::Handle;
use nacre_tess::{TessConfig, Tessellation, tessellate};
use nacre_topo::{Face, Model, Shell, Solid};

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
    Profile2d {
        points: vec![p2(a, a), p2(b, a), p2(b, b), p2(a, b)],
    }
}

/// `PocketOnFace`'s profile lives in a frame **derived from the face** (design §6),
/// centred on it — not in world coordinates. `±0.2` here is the `[0.3,0.7]²` void.
fn centred_square(half: f64) -> Profile2d {
    Profile2d {
        points: vec![
            p2(-half, -half),
            p2(half, -half),
            p2(half, half),
            p2(-half, half),
        ],
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
            profile: centred_square(0.2),
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
    let l = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
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
    let r = boolean(&mut m, kind, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// `Cut(stub, L)` — the operands of `l_and_dimple` reversed. The L's top face keeps only
/// the stub's footprint, so the answer is a `0.4 × 0.4 × 0.5` box whose floor is that
/// island face: an outer loop made of nothing but seam vertices (cell 3f-2).
fn island_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
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
    let r = boolean(&mut m, BoolKind::Cut, b, a).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// The L-prism and an L-shaped bar in its notch, biting two corners of the cap. Two
/// chords on one face; under `Cut` the bar's floor splits into the two bites' floors.
fn notch_bar_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
    let bar = Profile2d {
        points: vec![
            p2(1.8, 0.8),
            p2(2.1, 0.8),
            p2(2.1, 2.1),
            p2(0.8, 2.1),
            p2(0.8, 1.8),
            p2(1.8, 1.8),
        ],
    };
    let mut m = replay(&[
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: l,
            dist: 1.0,
        },
        Operation::Extrude {
            plane: SketchPlane {
                origin: Point3::from_array([0.0, 0.0, 0.5]),
                ..SketchPlane::world_xy()
            },
            profile: bar,
            dist: 1.0,
        },
    ])
    .unwrap();
    let (a, b) = (m.live_solids[0], m.live_solids[1]);
    let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
    m.rebuild_adjacency();
    (m, r)
}

/// The L-prism with an L-shaped stub standing wholly inside its cap: a blind pocket whose
/// lid carries a **non-convex** inner loop, the first the bridging triangulator has seen.
fn ell_dimple_cut() -> (Model, Handle<Solid>) {
    let l = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
    let ell = Profile2d {
        points: vec![
            p2(0.2, 0.25),
            p2(0.85, 0.25),
            p2(0.85, 0.4),
            p2(0.35, 0.4),
            p2(0.35, 0.9),
            p2(0.2, 0.9),
        ],
    };
    let mut m = replay(&[
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile: l,
            dist: 1.0,
        },
        Operation::Extrude {
            plane: SketchPlane {
                origin: Point3::from_array([0.0, 0.0, 0.5]),
                ..SketchPlane::world_xy()
            },
            profile: ell,
            dist: 1.0,
        },
    ])
    .unwrap();
    let (a, b) = (m.live_solids[0], m.live_solids[1]);
    let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
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
    let r = boolean(&mut m, BoolKind::Cut, u, slab).unwrap();
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

    // The island is the only face all of whose vertices are `Discovered` — the stub's
    // side pieces each keep two of the original box's corners.
    let faces = m.shells.get(m.solids.get(s).outer).faces.clone();
    let isle = *faces
        .iter()
        .find(|&&f| {
            m.faces.get(f).outer.half_edges.iter().all(|he| {
                let b = m.edges.get(he.edge).bounds.unwrap();
                b.iter().all(|&v| {
                    matches!(
                        m.vertices.get(v).origin,
                        nacre_topo::Origin::Discovered { .. }
                    )
                })
            })
        })
        .expect("the island face");

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
    // face store meshed both original cubes (6 + 6) alongside the union.
    let (m, s) = stacked_fuse();
    let g = mesh_vs_props(&m, s);
    assert_eq!(g.tris, 20);
    assert_agrees(&g, "stacked fuse");
}
