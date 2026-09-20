use super::*;
use nacre_geom::Plane;
use nacre_math::Vector3;
use nacre_store::Store;
use nacre_topo::PointCache;
use nacre_topo::{HalfEdge, Orientation, Shell, Solid};
use proptest::prelude::*;

fn cuboid(min: [f64; 3], max: [f64; 3]) -> Model {
    let mut m = Model::new();
    m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
    m
}

fn cylinder(base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Model {
    let mut m = Model::new();
    m.add_cylinder(Point3::from_array(base), Vector3::from_array(axis), r, h);
    m
}

/// A `Handle<Shell>` for `index` (minted from a throwaway store; reference
/// integrity reads only `.index()`, so the debug store-id guard is untouched).
fn shell_handle_at(index: u32) -> Handle<Shell> {
    let mut s: Store<Shell> = Store::new();
    let mut h = s.push(Shell { faces: vec![] });
    for _ in 0..index {
        h = s.push(Shell { faces: vec![] });
    }
    h
}

/// A `Handle<Surface>` for `index` (same throwaway-store trick). The store holds the
/// **truth**, as the model's does — only `.index()` is read, so any statement will do.
fn surface_handle_at(index: u32) -> Handle<Surface> {
    let plane = || Surface::Plane {
        points: nacre_topo::PlanePoints::Known(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]].map(|p| p.map(nacre_scalar::Rat::from_int)),
        ),
        motion: None,
    };
    let mut s: Store<Surface> = Store::new();
    let mut h = s.push(plane());
    for _ in 0..index {
        h = s.push(plane());
    }
    h
}

#[test]
fn cube_is_clean() {
    assert!(validate(&cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).is_empty());
}

#[test]
fn asymmetric_cuboid_is_clean() {
    assert!(validate(&cuboid([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).is_empty());
}

#[test]
fn cylinder_is_clean() {
    // The seam b-rep (V2/E3/F3 with a self-adjacent seam edge, single-half-edge
    // cap loops, and start==end circle edges) validates clean unmodified.
    let v = validate(&cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0));
    assert!(v.is_empty(), "{v:?}");
}

#[test]
fn slanted_offset_cylinder_is_clean() {
    let v = validate(&cylinder([3.0, -1.0, 2.0], [1.0, 2.0, 3.0], 0.7, 4.0));
    assert!(v.is_empty(), "{v:?}");
}

#[test]
fn dangling_reference_solid_shell() {
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    // shells.len() == 1; add a solid (index 1) pointing at shell index 5.
    m.push_solid_unlisted(Solid {
        outer: shell_handle_at(5),
        cavities: vec![],
    });
    assert_eq!(
        validate(&m),
        vec![Violation::DanglingReference {
            kind: RefKind::SolidShell,
            owner_index: 1,
            target_index: 5,
            target_len: 1,
        }]
    );
}

#[test]
fn dangling_reference_vertex_definition() {
    // A discovered vertex whose definition points past the surface store.
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]); // 6 surfaces, 8 vertices
    let vh = m.push_vertex(
        Vertex::ThreePlane([
            surface_handle_at(0),
            surface_handle_at(1),
            surface_handle_at(9), // out of bounds — only 6 surfaces
        ]),
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    assert_eq!(
        validate(&m),
        vec![Violation::DanglingReference {
            kind: RefKind::VertexSurface,
            owner_index: vh.index(),
            target_index: 9,
            target_len: 6,
        }]
    );
}

#[test]
fn vertex_def_is_copy() {
    let d = Vertex::ThreePlane([surface_handle_at(0); 3]);
    let copy = d; // move-or-copy
    let _again = d; // still usable ⇒ Copy, not moved
    assert_eq!(d, copy);
}

#[test]
fn stray_vertex_ignored() {
    // A vertex referenced by nothing is a dead arena item, not a defect: it
    // is unreachable from the live cube, so validate ignores it.
    // (A whole-store count would raise EulerParity{v:9}.)
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    m.push_vertex(
        Vertex::ThreePlane([
            m.world_plane(nacre_scalar::Axis::Z),
            m.world_plane(nacre_scalar::Axis::X),
            m.world_plane(nacre_scalar::Axis::Y),
        ]),
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    assert!(validate(&m).is_empty());
}

#[test]
fn supersede_solid_reuses_cells() {
    // What an editing op does: build a new solid that *reuses*
    // the old solid's cells, then drop the old solid from `live_solids`. The
    // old cells linger in the arena but, unreferenced by any live solid, fall
    // out of the reachable closure — so each shared edge is counted twice
    // (once per live face), not four times, and validate stays clean.
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let old = m.live_solids()[0];
    let old_shell = m.solid(old).outer;
    let faces = m.shell(old_shell).faces.clone();
    let new_shell = m.push_shell(Shell { faces });
    let new_solid = m.push_solid(Solid {
        outer: new_shell,
        cavities: vec![],
    });
    m.supersede_live(&[old]); // supersede: only the new one is live
    assert_eq!(m.live_solids().to_vec(), vec![new_solid]);
    let v = validate(&m);
    assert!(v.is_empty(), "{v:?}");
}

#[test]
fn orphaned_face_is_ignored() {
    // A face pushed into the store but not part of any live solid is
    // unreachable, so it adds nothing to Euler and does not pollute manifold
    // use-counts. It duplicates a cube face (reusing that face's surface and
    // edges), so reference integrity stays clean; were it counted, those
    // edges would be over-used (non-manifold) and F would rise by one.
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let dup = {
        let f0 = m.face(m.face_handle_at(0).expect("a face"));
        Face {
            surface: f0.surface,
            outer: f0.outer.clone(),
            inner: vec![],
            orientation: f0.orientation,
        }
    };
    m.push_face(dup);
    let v = validate(&m);
    assert!(v.is_empty(), "{v:?}");
}

// --- hand-built tetrahedron (smallest closed 2-manifold), verified winding ---

const TETRA_BASE: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETRA_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
/// (plane-defining vertex triple, half-edges as (edge index, forward)).
type TetraFace = ([usize; 3], [(usize, bool); 3]);
const TETRA_FACES: [TetraFace; 4] = [
    ([0, 2, 1], [(1, true), (3, false), (0, false)]),
    ([0, 1, 3], [(0, true), (4, true), (2, false)]),
    ([0, 3, 2], [(2, true), (5, false), (1, false)]),
    ([1, 2, 3], [(3, true), (5, true), (4, false)]),
];

#[derive(Default)]
struct TetraOpts {
    /// (vertex index, coordinate delta) — moves a vertex off its (un-moved) surfaces.
    ///
    /// ⚠ It used to carry a third element, a tolerance to install on that vertex, because the
    /// cache stored a measured residual and a checker read it. Nothing stores one now: every
    /// vertex is held to [`EPS_CONSTRUCTED`], so a per-vertex knob would drive nothing.
    nudge: Option<(usize, [f64; 3])>,
    drop_face: Option<usize>,
    flip_he: Option<(usize, usize)>, // (face, half-edge position)
    swap_he: Option<(usize, usize, usize)>, // (face, position a, position b)
}

/// Push a standard tetra translated by `t` (with `opts` defects) into `m`;
/// return its face handles. Curves/planes go through the *un-nudged* corners.
fn push_tetra(m: &mut Model, t: [f64; 3], opts: &TetraOpts) -> Vec<Handle<Face>> {
    let corner = |i: usize| {
        Point3::from_array([
            TETRA_BASE[i][0] + t[0],
            TETRA_BASE[i][1] + t[1],
            TETRA_BASE[i][2] + t[2],
        ])
    };
    // Surfaces first (they need only corners), so a discovered vertex's
    // definition can reference their handles. Built for all four faces even
    // if one is dropped — an orphan surface is outside the reachable set, so
    // it does not perturb Euler counts or reference checks.
    let sh: Vec<Handle<Surface>> = TETRA_FACES
        .iter()
        .map(|(tri, _)| {
            let lift = |p: Point3| {
                p.as_array()
                    .map(|x| nacre_scalar::Rat::from_decimal(x).expect("tetra corners"))
            };
            let (h, _) = m.push_plane(
                Plane::through_points(corner(tri[0]), corner(tri[1]), corner(tri[2])).unwrap(),
                [
                    lift(corner(tri[0])),
                    lift(corner(tri[1])),
                    lift(corner(tri[2])),
                ],
                None,
            );
            h
        })
        .collect();
    let vh: Vec<Handle<Vertex>> = (0..4)
        .map(|i| {
            let mut p = corner(i).as_array();
            if let Some((vi, d)) = opts.nudge {
                if vi == i {
                    p = [p[0] + d[0], p[1] + d[1], p[2] + d[2]];
                }
            }
            // Every vertex names the three tetra faces incident to it (the
            // definition is the vertex).
            let incident: Vec<Handle<Surface>> = TETRA_FACES
                .iter()
                .enumerate()
                .filter(|(_, (tri, _))| tri.contains(&i))
                .map(|(fi, _)| sh[fi])
                .collect();
            let coord = Point3::from_array(p);
            m.push_vertex(
                Vertex::ThreePlane(incident.try_into().expect("a tetra vertex is on 3 faces")),
                PointCache::Unrealized { coord },
            )
        })
        .collect();
    let eh: Vec<Handle<Edge>> = (0..6)
        .map(|i| {
            let (a, b) = TETRA_EDGES[i];
            // The two faces whose loops use edge `i` — its carriers, read off the same
            // table the loops are built from. `push_edge` derives the curve through the
            // (possibly nudged) vertex points — an endpoint cannot sit off its own
            // line (see `vertex_off_its_rim_circle`).
            let carriers: Vec<Handle<Surface>> = TETRA_FACES
                .iter()
                .enumerate()
                .filter(|(_, (_, hes))| hes.iter().any(|&(e, _)| e == i))
                .map(|(fi, _)| sh[fi])
                .collect();
            let [ca, cb] = carriers[..] else {
                panic!("a tetra edge is on exactly 2 faces")
            };
            m.push_edge([ca, cb], [vh[a], vh[b]])
                .expect("distinct tetra corners")
        })
        .collect();
    let mut faces = Vec::new();
    for (fi, (_tri, hes)) in TETRA_FACES.iter().enumerate() {
        if opts.drop_face == Some(fi) {
            continue;
        }
        let surface = sh[fi];
        let mut half_edges: Vec<HalfEdge> = hes
            .iter()
            .map(|&(e, forward)| HalfEdge {
                edge: eh[e],
                forward,
            })
            .collect();
        if let Some((face_i, pos)) = opts.flip_he {
            if face_i == fi {
                half_edges[pos].forward = !half_edges[pos].forward;
            }
        }
        if let Some((face_i, a, b)) = opts.swap_he {
            if face_i == fi {
                half_edges.swap(a, b);
            }
        }
        faces.push(m.push_face_unchecked(Face {
            surface,
            outer: Loop { half_edges },
            inner: vec![],
            orientation: Orientation::Forward,
        }));
    }
    faces
}

fn tetra_with(opts: TetraOpts) -> Model {
    let mut m = Model::new();
    let faces = push_tetra(&mut m, [0.0; 3], &opts);
    let sh = m.push_shell_unchecked(Shell { faces });
    m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    m
}

fn tetra() -> Model {
    tetra_with(TetraOpts::default())
}

#[test]
fn tetra_is_clean() {
    assert!(validate(&tetra()).is_empty());
}

#[test]
fn open_loop_from_swapped_half_edges() {
    let vs = validate(&tetra_with(TetraOpts {
        swap_he: Some((0, 0, 2)),
        ..Default::default()
    }));
    assert!(!vs.is_empty());
    assert!(vs.iter().all(|v| matches!(
        v,
        Violation::OpenLoop {
            loop_kind: LoopKind::Outer,
            ..
        }
    )));
}

#[test]
fn non_manifold_edge_from_dropped_face() {
    let vs = validate(&tetra_with(TetraOpts {
        drop_face: Some(0),
        ..Default::default()
    }));
    assert!(
        vs.iter()
            .any(|v| matches!(v, Violation::NonManifoldEdge { use_count: 1, .. }))
    );
}

#[test]
fn non_opposed_edge_from_flipped_half_edge() {
    let vs = validate(&tetra_with(TetraOpts {
        flip_he: Some((0, 0)),
        ..Default::default()
    }));
    assert!(
        vs.iter()
            .any(|v| matches!(v, Violation::NonOpposedEdge { .. }))
    );
}

#[test]
fn negative_genus_two_tetra_one_shell() {
    let mut m = Model::new();
    let mut faces = push_tetra(&mut m, [0.0; 3], &TetraOpts::default());
    faces.extend(push_tetra(&mut m, [5.0, 0.0, 0.0], &TetraOpts::default()));
    let sh = m.push_shell(Shell { faces });
    m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    assert_eq!(
        validate(&m),
        vec![Violation::NegativeGenus {
            v: 8,
            e: 12,
            f: 8,
            s: 1,
            inner_loops: 0,
            genus: -1,
        }]
    );
}

/// ★ Positive control for the **all-vertex** definition check: a vertex with *no* measured
/// tolerance — the constructed population — is caught when its
/// coordinate drifts off the surfaces its definition names. Without this, "the check
/// fired zero times across the suite" could equally mean "the check sees nothing there".
#[test]
fn a_built_vertex_off_its_definition_is_caught() {
    let m = tetra_with(TetraOpts {
        // Constructed (tol `None` ⇒ EPS_CONSTRUCTED), nudged well past that epsilon.
        nudge: Some((0, [0.0, 0.0, 1e-6])),
        ..Default::default()
    });
    let vs = validate(&m);
    assert!(
        vs.iter().any(|v| matches!(
            v,
            Violation::VertexOffDefinition { tol, .. } if *tol == EPS_CONSTRUCTED
        )),
        "a built vertex off its own definition must be caught: {vs:?}"
    );
}

#[test]
fn vertex_off_surface_when_nudged() {
    // Vertex 0 moved by 2·EPS_CONSTRUCTED; the planes stay on the un-moved corners, so
    // the vertex is off them. (`VertexOffCurve` cannot fire here: a line edge's
    // curve derives *through its endpoints*, so an endpoint is on its own line by
    // construction — the check's remaining teeth are the circles, see
    // `vertex_off_its_rim_circle`.)
    let vs = validate(&tetra_with(TetraOpts {
        nudge: Some((0, [0.0, 0.0, 2.0 * EPS_CONSTRUCTED])),
        ..Default::default()
    }));
    assert!(
        vs.iter()
            .any(|v| matches!(v, Violation::VertexOffSurface { .. }))
    );
    assert!(
        !vs.iter()
            .any(|v| matches!(v, Violation::VertexOffCurve { .. })),
        "a line endpoint is on its own derived line by construction"
    );
}

/// ★ Negative control: an edge whose stated carriers disagree with the two faces
/// actually using it is flagged — and so is a self-adjacent *plane* pair (a spelling
/// reserved for cylinder seams). Hand-built open surface: other violations fire too; the
/// assertion is only that the carrier one is among them.
#[test]
fn edge_carrier_mismatch_is_flagged() {
    let mut m = nacre_topo::Model::new();
    let r = nacre_scalar::Rat::from_int;
    let plane = |m: &mut nacre_topo::Model, n: [f64; 3], pts: [[i128; 3]; 3]| {
        m.push_plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array(n)).unwrap(),
            pts.map(|p| p.map(r)),
            None,
        )
        .0
    };
    let sa = plane(&mut m, [0.0, 0.0, 1.0], [[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
    let sb = plane(&mut m, [0.0, 1.0, 1.0], [[0, 0, 0], [1, 0, 0], [0, 1, -1]]);
    let sc = plane(&mut m, [1.0, 0.0, 1.0], [[0, 0, 0], [0, 1, 0], [1, 0, -1]]);
    let v = |m: &mut nacre_topo::Model, p: [f64; 3]| {
        m.push_vertex(
            Vertex::ThreePlane([sa, sb, sc]),
            PointCache::Unrealized {
                coord: Point3::from_array(p),
            },
        )
    };
    let v0 = v(&mut m, [0.0, 0.0, 0.0]);
    let v1 = v(&mut m, [1.0, 0.0, 0.0]);
    let v2 = v(&mut m, [0.0, 1.0, 0.0]);
    let v3 = v(&mut m, [0.0, -1.0, 0.0]);
    // The shared edge states carriers [sa, sc] — but the faces using it are on sa and sb.
    let e_shared = m.push_edge([sa, sc], [v0, v1]).unwrap();
    let ea1 = m.push_edge([sa, sb], [v1, v2]).unwrap();
    let ea2 = m.push_edge([sa, sb], [v2, v0]).unwrap();
    let eb1 = m.push_edge([sa, sb], [v1, v3]).unwrap();
    let eb2 = m.push_edge([sa, sb], [v3, v0]).unwrap();
    let he = |edge, forward| HalfEdge { edge, forward };
    let face = |m: &mut nacre_topo::Model, surface, hes: Vec<HalfEdge>| {
        m.push_face(Face {
            surface,
            outer: Loop { half_edges: hes },
            inner: vec![],
            orientation: Orientation::Forward,
        })
    };
    let fa = face(
        &mut m,
        sa,
        vec![he(e_shared, true), he(ea1, true), he(ea2, true)],
    );
    let fb = face(
        &mut m,
        sb,
        vec![he(e_shared, false), he(eb2, false), he(eb1, false)],
    );
    let shell = m.push_shell(Shell {
        faces: vec![fa, fb],
    });
    m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    let vs = validate(&m);
    assert!(
        vs.iter()
            .any(|x| matches!(x, Violation::EdgeCarrierMismatch { edge, .. } if *edge == e_shared)),
        "stated [sa, sc] vs observed [sa, sb] must be flagged: {vs:?}"
    );
    // The other shared-plane edges (1 use each) are non-manifold, not carrier mismatches.
    assert!(
        vs.iter()
            .any(|x| matches!(x, Violation::NonManifoldEdge { .. }))
    );
}

/// ★ The same-surface arm (the cylinder panel population): two faces on **one**
/// surface sharing an edge cannot spell the carrier pair between them, so the rule there is
/// membership — the shared surface must be one of the stated two (the *other* stated
/// carrier is what the curve derivation crosses it with). Both directions: membership
/// passes, non-membership is still flagged.
#[test]
fn same_surface_users_check_membership_not_equality() {
    // The arena is append-only, so each direction builds its own model.
    let build = |stated_second_is_sc: bool| {
        let mut m = nacre_topo::Model::new();
        let r = nacre_scalar::Rat::from_int;
        let plane = |m: &mut nacre_topo::Model, n: [f64; 3], pts: [[i128; 3]; 3]| {
            m.push_plane(
                Plane::from_point_normal(Point3::origin(), Vector3::from_array(n)).unwrap(),
                pts.map(|p| p.map(r)),
                None,
            )
            .0
        };
        let sa = plane(&mut m, [0.0, 0.0, 1.0], [[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
        let sb = plane(&mut m, [0.0, 1.0, 1.0], [[0, 0, 0], [1, 0, 0], [0, 1, -1]]);
        let sc = plane(&mut m, [1.0, 0.0, 1.0], [[0, 0, 0], [0, 1, 0], [1, 0, -1]]);
        let v = |m: &mut nacre_topo::Model, p: [f64; 3]| {
            m.push_vertex(
                Vertex::ThreePlane([sa, sb, sc]),
                PointCache::Unrealized {
                    coord: Point3::from_array(p),
                },
            )
        };
        let v0 = v(&mut m, [0.0, 0.0, 0.0]);
        let v1 = v(&mut m, [1.0, 0.0, 0.0]);
        let v2 = v(&mut m, [0.0, 1.0, 0.0]);
        let v3 = v(&mut m, [0.0, -1.0, 0.0]);
        // Both users of the shared edge sit on `sa`; the stated pair either contains it
        // ([sa, sc] — legal) or does not ([sb, sc] — a mismatch).
        let stated = if stated_second_is_sc {
            [sa, sc]
        } else {
            [sb, sc]
        };
        let e_shared = m.push_edge(stated, [v0, v1]).unwrap();
        let ea1 = m.push_edge([sa, sb], [v1, v2]).unwrap();
        let ea2 = m.push_edge([sa, sb], [v2, v0]).unwrap();
        let eb1 = m.push_edge([sa, sb], [v1, v3]).unwrap();
        let eb2 = m.push_edge([sa, sb], [v3, v0]).unwrap();
        let he = |edge, forward| HalfEdge { edge, forward };
        let face = |m: &mut nacre_topo::Model, surface, hes: Vec<HalfEdge>| {
            m.push_face(Face {
                surface,
                outer: Loop { half_edges: hes },
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        let fa = face(
            &mut m,
            sa,
            vec![he(e_shared, true), he(ea1, true), he(ea2, true)],
        );
        let fb = face(
            &mut m,
            sa,
            vec![he(e_shared, false), he(eb2, false), he(eb1, false)],
        );
        let shell = m.push_shell(Shell {
            faces: vec![fa, fb],
        });
        m.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        (m, e_shared)
    };
    let (m, e_shared) = build(true);
    let vs = validate(&m);
    assert!(
        !vs.iter()
            .any(|x| matches!(x, Violation::EdgeCarrierMismatch { edge, .. } if *edge == e_shared)),
        "same-surface users on a stated carrier are legal: {vs:?}"
    );
    let (m, e_shared) = build(false);
    let vs = validate(&m);
    assert!(
        vs.iter()
            .any(|x| matches!(x, Violation::EdgeCarrierMismatch { edge, .. } if *edge == e_shared)),
        "same-surface users off both stated carriers must be flagged: {vs:?}"
    );
}

/// ★ Negative controls: a definition variant whose carrier kinds contradict it is
/// flagged — `ThreePlane` naming a cylinder, and `OnSeam` naming two planes.
#[test]
fn a_contradictory_vertex_def_is_flagged() {
    let mut m = nacre_topo::Model::new();
    m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        5.0,
    );
    let cyl = (0..m.face_count() as u32)
        .filter_map(|i| m.face_handle_at(i))
        .map(|h| (h, m.face(h)))
        .map(|(_, f)| f.surface)
        .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("the lateral cylinder");
    let (z0, x0) = (
        m.world_plane(nacre_scalar::Axis::Z),
        m.world_plane(nacre_scalar::Axis::X),
    );
    let bad_three = m.push_vertex(
        Vertex::ThreePlane([z0, x0, cyl]),
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    let bad_seam = m.push_vertex(
        Vertex::OnSeam([z0, x0]),
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    let vs = validate(&m);
    for bad in [bad_three, bad_seam] {
        assert!(
            vs.iter().any(
                |v| matches!(v, Violation::VertexCarrierMismatch { vertex } if *vertex == bad)
            ),
            "a contradictory def must be flagged: {vs:?}"
        );
    }
}

/// ★★ **The positive control for a cap** — the loop with no polygon in it.
///
/// A cylinder cap's boundary is one closed rim, so the Newell path above has
/// nothing to wind. Skipping it would let a cylinder ship with
/// **both caps facing the same way** (a cap plane that interns with an
/// existing, oppositely-stated plane), a solid no other check here can see —
/// edge opposition, Euler and the signed volume are all blind to it (the
/// broken cylinder's signed volume stays positive, just wrong).
///
/// The witness the rim does have is its circle, and its direction comes from
/// the *cylinder's* axis rather than from this plane — which is what makes it
/// evidence rather than an echo.
#[test]
fn a_cap_whose_flag_lies_is_caught_by_its_rim() {
    let mut m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 2.0);
    let solid = m.live_solids()[0];
    let shell = m.solid(solid).outer;
    let faces = m.shell(shell).faces.clone();
    // The bottom cap: a planar face whose outer loop is a single half-edge.
    let victim = *faces
        .iter()
        .find(|&&fh| {
            let f = m.face(fh);
            matches!(m.surface_cache(f.surface), nacre_geom::Surface::Plane(_))
                && f.outer.half_edges.len() == 1
        })
        .expect("a cylinder has two disk caps");
    let twin = {
        let f = m.face(victim).clone();
        m.push_face(Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.push_shell(Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![replaced]);
    m.rebuild_adjacency();
    let vs = validate(&m);
    assert_eq!(vs.len(), 1, "only the flag lie should fire: {vs:?}");
    assert!(
        matches!(vs[0], Violation::FaceMisoriented { face, cos } if face == twin && cos < -0.5),
        "{vs:?}"
    );
}

/// The other half of the pair: the same cylinder, untouched, is clean — so the
/// check above is measuring the lie and not the shape.
#[test]
fn an_honest_cylinder_passes_the_rim_check() {
    let mut m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 2.0);
    m.rebuild_adjacency();
    assert!(validate(&m).is_empty(), "{:?}", validate(&m));
}

/// ★ The positive control for `FaceMisoriented`. The stores are sealed against
/// mutation, so the defect is *built*: a twin of one cuboid wall pushed with its
/// orientation flag flipped, seated in a hand-made shell that supersedes the
/// original solid. Every other cell is shared and the twin traverses the same
/// half-edges, so edge opposition, Euler and incidence all still hold — this
/// check is the only one that can see the lie.
#[test]
fn a_flipped_orientation_flag_is_caught() {
    let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let solid = m.live_solids()[0];
    let shell = m.solid(solid).outer;
    let faces = m.shell(shell).faces.clone();
    let victim = faces[0];
    let twin = {
        let f = m.face(victim).clone();
        m.push_face(Face {
            orientation: f.orientation.flipped(),
            ..f
        })
    };
    let sh = m.push_shell(Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![replaced]);
    let vs = validate(&m);
    assert_eq!(vs.len(), 1, "only the flag lie should fire: {vs:?}");
    assert!(
        matches!(vs[0], Violation::FaceMisoriented { face, cos } if face == twin && cos < -0.5),
        "{vs:?}"
    );
}

/// ★ Positive control: the def–cache net (`CylinderTruthCacheMismatch`) actually
/// bites. The stores are sealed, so the defect is *built* the `FaceMisoriented` way: a twin
/// of the lateral face whose surface is the same cache pushed under a **lying def** (radius
/// 2 where the cache says 1). Every geometric check still passes — the loop's vertices sit
/// on the same cache — so a def that lies about the radius is visible to this net alone.
/// (The twin's edges still name the original lateral as carrier, so `EdgeCarrierMismatch`
/// fires too; the assertion is that the truth net is among the findings, per-field.)
#[test]
fn a_lying_cylinder_def_is_caught() {
    let mut m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 2.0);
    let solid = m.live_solids()[0];
    let shell = m.solid(solid).outer;
    let faces = m.shell(shell).faces.clone();
    let victim = *faces
        .iter()
        .find(|&&fh| {
            matches!(
                m.surface_cache(m.face(fh).surface),
                nacre_geom::Surface::Cylinder(_)
            )
        })
        .expect("the lateral face");
    let cache = match m.surface_cache(m.face(victim).surface) {
        nacre_geom::Surface::Cylinder(c) => *c,
        _ => unreachable!(),
    };
    let r = |x: f64| nacre_scalar::Rat::from_decimal(x).expect("decimal");
    let lying = nacre_topo::CylinderDef::new(
        [r(0.0), r(0.0), r(0.0)],
        [r(0.0), r(0.0), r(1.0)],
        [r(0.0), r(-1.0), r(0.0)],
        nacre_scalar::BigRat::from(r(4.0)), // the lie: r² = 4, where the cache's radius is 1
    )
    .expect("well-formed statement");
    let liar = m.push_cylinder(cache, lying, None);
    let twin = {
        let f = m.face(victim).clone();
        m.push_face(Face { surface: liar, ..f })
    };
    let sh = m.push_shell(Shell {
        faces: faces
            .iter()
            .map(|&h| if h == victim { twin } else { h })
            .collect(),
    });
    let replaced = m.push_solid(Solid {
        outer: sh,
        cavities: vec![],
    });
    m.restore_live(vec![replaced]);
    let vs = validate(&m);
    assert!(
        vs.iter().any(|v| matches!(
            v,
            Violation::CylinderTruthCacheMismatch { field: "radius", def_value, cache_value, .. }
                if *def_value == 2.0 && *cache_value == 1.0
        )),
        "the lying radius must be caught: {vs:?}"
    );
}

/// ★ A well-formed `Pierce` vertex passes every check, and its coordinate is held
/// to **all three** carriers — the cylinder included, which is the carrier the variant
/// adds. The pierce points of {z = 0} ∧ {x = 0} against the r = 2 z-cylinder are
/// (0, ∓2, 0); the good one is clean, the one 1e−3 off the cylinder is flagged by
/// `VertexOffDefinition` through the cylinder carrier while both plane residuals stay 0
/// (the fixture is deliberately open — other violations may fire; the assertions are
/// per-proposition).
#[test]
fn a_pierce_vertex_is_held_to_its_cylinder() {
    use nacre_topo::QuadRoot;
    let mut m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let lateral = (0..m.face_count() as u32)
        .filter_map(|i| m.face_handle_at(i))
        .map(|h| (h, m.face(h)))
        .map(|(_, f)| f.surface)
        .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("lateral");
    let bottom = m.world_plane(nacre_scalar::Axis::Z);
    let x0 = m.world_plane(nacre_scalar::Axis::X);
    // The good statement, free-floating: reference integrity and the carrier-kind check
    // run over every vertex, and neither may fire.
    let _good = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: lateral,
            root: QuadRoot::Lo,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0, -2.0, 0.0]),
        },
    );
    assert_eq!(validate(&m), vec![], "a sound pierce statement is clean");
    // The lying coordinate, wired into a reachable face so the off-definition check
    // sees it.
    let bad = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: lateral,
            root: QuadRoot::Hi,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0, 2.001, 0.0]),
        },
    );
    let anchor = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: lateral,
            root: QuadRoot::Lo,
        },
        PointCache::Unrealized {
            coord: Point3::from_array([0.0, -2.0, 0.0]),
        },
    );
    let edge = m.push_edge([bottom, x0], [anchor, bad]).expect("a line");
    let face = m.push_face_unchecked(Face {
        surface: x0,
        outer: Loop {
            half_edges: vec![HalfEdge {
                edge,
                forward: true,
            }],
        },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    let shell = m.push_shell_unchecked(Shell { faces: vec![face] });
    m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    let vs = validate(&m);
    assert!(
        vs.iter().any(|v| matches!(
            v,
            Violation::VertexOffDefinition { vertex, surface_index, .. }
                if *vertex == bad && *surface_index == lateral.index()
        )),
        "the cylinder carrier must catch the lying coordinate: {vs:?}"
    );
    assert!(
        !vs.iter().any(|v| matches!(
            v,
            Violation::VertexOffDefinition { vertex, .. } if *vertex == anchor
        )),
        "the honest coordinate stays clean: {vs:?}"
    );
}

/// ★ Negative controls: a `Pierce` whose carrier kinds contradict the structure is
/// flagged — a cylinder in a plane slot, and a plane in the cylinder slot.
#[test]
fn a_contradictory_pierce_def_is_flagged() {
    use nacre_topo::QuadRoot;
    let mut m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
    let lateral = (0..m.face_count() as u32)
        .filter_map(|i| m.face_handle_at(i))
        .map(|h| (h, m.face(h)))
        .map(|(_, f)| f.surface)
        .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("lateral");
    let bottom = m.world_plane(nacre_scalar::Axis::Z);
    let x0 = m.world_plane(nacre_scalar::Axis::X);
    let cyl_in_plane_slot = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, lateral],
            cylinder: x0,
            root: QuadRoot::Lo,
        },
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    let plane_in_cyl_slot = m.push_vertex(
        Vertex::Pierce {
            planes: [bottom, x0],
            cylinder: m.world_plane(nacre_scalar::Axis::Y),
            root: QuadRoot::Lo,
        },
        PointCache::Unrealized {
            coord: Point3::origin(),
        },
    );
    let vs = validate(&m);
    for bad in [cyl_in_plane_slot, plane_in_cyl_slot] {
        assert!(
            vs.iter().any(
                |v| matches!(v, Violation::VertexCarrierMismatch { vertex } if *vertex == bad)
            ),
            "a contradictory pierce def must be flagged: {vs:?}"
        );
    }
}

/// ★ The population that refutes an ulp-based metric for the net, pinned. A
/// full-width random axis makes the cache's Gram–Schmidt smear ~1e−17 into the component
/// the def's shuffle keeps at exactly `0.0` — ulps are meaningless across zero, so the
/// bound is mixed absolute/relative ([`CYL_TRUTH_EPS`]). This exact axis is a proptest
/// minimal failing input under the ulp metric; it must stay clean.
#[test]
fn a_full_width_axis_survives_the_truth_net() {
    let m = cylinder(
        [0.0, 0.0, 0.0],
        [
            0.1181674690543037,
            -0.11258144547841069,
            -0.033712610662548576,
        ],
        0.5,
        0.1,
    );
    assert_eq!(validate(&m), vec![]);
}

/// ★ The check `VertexOffCurve` has teeth where the curve does NOT derive from
/// the vertex — a rim circle comes from the carriers (cylinder axis × cap), so a seam
/// vertex off the rim is caught. The stores are sealed against mutation, so the defect is
/// *built*: a disk face whose rim edge carries a seam vertex at the wrong radius. (The
/// fixture is deliberately not a closed solid — other violations fire too; the assertion
/// is only that this one is among them.)
#[test]
fn vertex_off_its_rim_circle() {
    let mut m = nacre_topo::Model::new();
    m.add_cylinder(
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        2.0,
        5.0,
    );
    let lateral = (0..m.face_count() as u32)
        .filter_map(|i| m.face_handle_at(i))
        .map(|h| (h, m.face(h)))
        .map(|(_, f)| f.surface)
        .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("the cylinder's lateral surface");
    let cap = (0..m.face_count() as u32)
        .filter_map(|i| m.face_handle_at(i))
        .map(|h| (h, m.face(h)))
        .map(|(_, f)| f.surface)
        .find(|&h| h != lateral && !matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("a cap plane");
    // A seam vertex at radius 2 + 1e-3 — off the derived radius-2 rim circle. Its
    // definition is `OnSeam` (the lateral cylinder and its cap), which is what a rim's
    // endpoint always is; the coordinate is the part that lies.
    let bad = m.push_vertex(
        Vertex::OnSeam([lateral, cap]),
        PointCache::Unrealized {
            coord: Point3::from_array([2.001, 0.0, 0.0]),
        },
    );
    let rim = m
        .push_edge([lateral, cap], [bad, bad])
        .expect("a rim derives from its carriers, not its vertices");
    let face = m.push_face(Face {
        surface: cap,
        outer: Loop {
            half_edges: vec![HalfEdge {
                edge: rim,
                forward: true,
            }],
        },
        inner: vec![],
        orientation: Orientation::Forward,
    });
    let shell = m.push_shell(Shell { faces: vec![face] });
    m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    let vs = validate(&m);
    assert!(
        vs.iter()
            .any(|v| matches!(v, Violation::VertexOffCurve { vertex, .. } if *vertex == bad)),
        "a seam vertex off its rim circle must be caught: {vs:?}"
    );
}

// ⚠★★★ **A third test stood here and its premise is gone, not its coverage.**
// `discovered_vertex_within_tolerance_is_clean` nudged a vertex by half of the tolerance it
// installed on that same vertex and asserted the model still validated. Nothing installs a
// per-vertex tolerance any more — every vertex is held to `EPS_CONSTRUCTED` — so the
// proposition "a vertex may sit within *its own recorded* tolerance" has no subject. The
// clean-model side is not lost: every fixture in this file that is *not* nudged asserts it.

#[test]
fn a_vertex_nudged_off_its_surfaces_flags() {
    let vs = validate(&tetra_with(TetraOpts {
        nudge: Some((0, [0.0, 0.0, 2e-6])),
        ..Default::default()
    }));
    assert!(vs.iter().any(|v| matches!(
        v,
        Violation::VertexOffSurface { .. } | Violation::VertexOffCurve { .. }
    )));
}

#[test]
fn a_vertex_nudged_off_its_definition_flags() {
    // Nudged past the epsilon, the vertex is off its three definition planes (its incident
    // faces). Assert the *definition* check fires specifically — not merely "some violation",
    // which `VertexOffSurface` alone would satisfy.
    let vs = validate(&tetra_with(TetraOpts {
        nudge: Some((0, [0.0, 0.0, 2e-6])),
        ..Default::default()
    }));
    assert!(
        vs.iter()
            .any(|v| matches!(v, Violation::VertexOffDefinition { .. })),
        "{vs:?}"
    );
}

// --- cavity orientation (M5 containment) ---

/// A hollow `outer`-cube with a concentric `inner`-cube void. `reverse`
/// selects a correct inward cavity (`reversed_shell`) or a mis-oriented one
/// (the inner shell used as-is, normals still pointing outward). Supersedes
/// the two source cubes so only the hollow solid is live.
fn hollow_cube(outer: f64, inner: f64, reverse: bool) -> Model {
    let mut m = Model::new();
    let ext = |s: f64| Vector3::from_array([s, s, s]);
    let a = m.add_cuboid(Point3::origin(), Point3::origin() + ext(outer));
    let inner_min = Point3::origin() + ext(0.5 * (outer - inner));
    let b = m.add_cuboid(inner_min, inner_min + ext(inner));
    let b_outer = m.solid(b).outer;
    let void = if reverse {
        m.reversed_shell(b_outer)
    } else {
        b_outer
    };
    let a_outer = m.solid(a).outer;
    let hollow = m.push_solid(Solid {
        outer: a_outer,
        cavities: vec![void],
    });
    m.restore_live(vec![hollow]);
    m
}

#[test]
fn correctly_oriented_cavity_is_clean() {
    let v = validate(&hollow_cube(4.0, 2.0, true));
    assert!(v.is_empty(), "{v:?}");
}

#[test]
fn misoriented_cavity_is_flagged() {
    // The un-reversed inner shell as a cavity: manifold and Euler still pass
    // (V16/E24/F12/S2), so the *only* violation is the orientation one.
    let v = validate(&hollow_cube(4.0, 2.0, false));
    assert_eq!(v.len(), 1, "{v:?}");
    assert!(
        matches!(v[0], Violation::CavityMisoriented { signed_volume, .. } if signed_volume > 0.0),
        "{v:?}"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn prop_random_box_is_clean(
        min in prop::array::uniform3(-1e3f64..1e3),
        ext in prop::array::uniform3(1e-2f64..1e3),
    ) {
        let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
        prop_assert!(validate(&cuboid(min, max)).is_empty());
    }

    #[test]
    fn prop_random_cylinder_is_clean(
        base in prop::array::uniform3(-1e3f64..1e3),
        axis in prop::array::uniform3(-1.0f64..1.0),
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        let axis = Vector3::from_array(axis);
        prop_assume!(axis.norm() > 0.1); // skip near-zero axes
        let mut m = Model::new();
        m.add_cylinder(Point3::from_array(base), axis, r, h);
        prop_assert!(validate(&m).is_empty());
    }
}
