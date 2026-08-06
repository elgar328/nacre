//! Mesh export for the nacre kernel.
//!
//! Two paths coexist during M3:
//! - [`to_obj`] — the M1 **bootstrap**: fan-triangulates a model's planar faces
//!   from edge endpoints (convex, planar only; panics on closed edges). Kept for
//!   the existing planar callers until they migrate.
//! - [`tessellate`] → [`Tessellation`] — the provenance-tagged layer (design §5):
//!   samples edges into shared polylines and triangulates faces (planar fans +
//!   ruled cylinder bands) crack-free, tagging every vertex with its origin.

mod polygon;

use nacre_geom::{Curve, Cylinder, Surface};
use nacre_math::Point3;
use nacre_store::{Handle, Store};
use nacre_topo::{Edge, Face, Loop, Model, Vertex};
use std::collections::HashMap;
use std::fmt::Write;

/// A face this tessellator cannot mesh. Returned rather than approximated: the
/// mesh is a cache, but a *wrong* cache is worse than none — silently fanning a
/// ring that is not star-shaped puts triangles outside the solid.
///
/// There is deliberately no `Fallback` variant. The old fan was one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TessError {
    /// The rings are not a polygon with sibling holes, so no triangulation of them
    /// exists: fewer than three vertices, a zero-area ring, a vertex used by two rings
    /// or repeated within one, two vertices at the same point, a spike, or a boundary
    /// that crosses itself.
    ///
    /// **The sweep detects these, where ear clipping used to notice them by accident**
    /// (it stalled, and that stall was reported as `NoEar`). It is checked rather than
    /// assumed because the b-rep guarantees it and this layer cannot: a decomposition
    /// handed a self-crossing ring would otherwise return a confident, wrong mesh.
    DegenerateRing,
    /// A hole wound the same way as its outer ring: the b-rep does not keep
    /// material on every loop's left. A broken solid, not a repairable mesh.
    HoleWinding,
    /// The bootstrap OBJ writer met a curved face. It used to fan the face's edge
    /// endpoints and emit nonsense; `tessellate` is the path that handles those.
    NonPlanarFace,
}

/// Export a model to Wavefront OBJ text (vertices shared; each face fan-
/// triangulated). Bootstrap — see the crate docs for the scope.
///
/// Only faces **reachable from `live_solids`** are emitted. `Store` is append-only
/// and `boolean`/`pocket`/`pad` supersede rather than delete, so iterating the face
/// store would mesh the operands alongside the result.
///
/// Every model vertex is still written, in store order, so an OBJ index stays a
/// vertex handle's index; a superseded vertex simply goes unreferenced.
///
/// Assumes every edge is bounded (M1); a closed edge (`bounds: None`, M3) would
/// panic. Faces are emitted in their loop winding, which for M1's outward-wound
/// `Orientation::Forward` faces yields outward-facing triangles.
pub fn to_obj(model: &Model) -> Result<String, TessError> {
    let mut out = String::new();
    // Writing to a String is infallible; unwrap keeps `unused_must_use` quiet.
    writeln!(
        out,
        "# nacre OBJ export (bootstrap: direct planar triangulation)"
    )
    .unwrap();

    // Vertices in store order → OBJ indices 1..=n (matches each handle's index).
    let pts: Vec<Point3> = model.vertices.iter().map(|(_, v)| v.point).collect();
    for p in &pts {
        let [x, y, z] = p.as_array();
        writeln!(out, "v {} {} {}", x, y, z).unwrap();
    }

    /// A loop's start vertices, as indices into the vertex store.
    fn ring(model: &Model, lp: &Loop) -> Vec<usize> {
        lp.half_edges
            .iter()
            .map(|he| {
                let bounds = model
                    .edges
                    .get(he.edge)
                    .bounds
                    .expect("M1: bounded edges only (closed edges arrive in M3)");
                let start = if he.forward { bounds[0] } else { bounds[1] };
                start.index() as usize
            })
            .collect()
    }

    let reach = model.reachable();
    // Store order, not `HashSet` order: the OBJ must be reproducible.
    for (fh, face) in model.faces.iter() {
        if !reach.faces.contains(&fh) {
            continue;
        }
        if !matches!(model.surface(face.surface), Surface::Plane(_)) {
            return Err(TessError::NonPlanarFace);
        }
        let outer = ring(model, &face.outer);
        let holes: Vec<Vec<usize>> = face.inner.iter().map(|lp| ring(model, lp)).collect();
        let holes: Vec<&[usize]> = holes.iter().map(|h| h.as_slice()).collect();
        for t in polygon::triangulate_polygon(&pts, &outer, &holes)? {
            writeln!(out, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1).unwrap();
        }
    }

    Ok(out)
}

// ---------------------------------------------------------------------------
// Provenance-tagged tessellation (design §5)
// ---------------------------------------------------------------------------

/// Where a tessellation vertex came from — the truth it can snap back to
/// (design §5). `t`/`uv` are the exact curve/surface parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TessOrigin {
    OnVertex(Handle<Vertex>),
    OnEdge { edge: Handle<Edge>, t: f64 },
    OnFace { face: Handle<Face>, uv: [f64; 2] },
}

/// A mesh vertex: a cached position plus its provenance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TessVertex {
    pub pos: Point3,
    pub origin: TessOrigin,
}

/// A mesh triangle, tagged with the face it approximates (hybrid-boolean fuel).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TessTriangle {
    pub vertices: [Handle<TessVertex>; 3],
    pub face: Handle<Face>,
}

/// A provenance-tagged triangle mesh derived from a [`Model`] (design §5).
///
/// `by_edge` holds each edge's shared polyline (the crack-free contract: faces
/// consume these, never re-sample), `by_face` the triangles per face. (The
/// incremental `stale` set of §5 is deferred — this is from-scratch, like
/// `Adjacency`.)
#[derive(Debug, Default)]
pub struct Tessellation {
    pub vertices: Store<TessVertex>,
    pub triangles: Store<TessTriangle>,
    pub by_edge: HashMap<Handle<Edge>, Vec<Handle<TessVertex>>>,
    pub by_face: HashMap<Handle<Face>, Vec<Handle<TessTriangle>>>,
}

/// Tessellation tolerance: the maximum chord deviation (sagitta) a sampled
/// polyline/mesh may have from the exact geometry.
#[derive(Clone, Copy, Debug)]
pub struct TessConfig {
    pub tol: f64,
}

impl Default for TessConfig {
    fn default() -> Self {
        Self { tol: 1e-2 }
    }
}

impl Tessellation {
    /// Wavefront OBJ text (vertices in store order, 1-based triangle indices).
    pub fn to_obj(&self) -> String {
        let mut out = String::new();
        writeln!(out, "# nacre OBJ export (provenance tessellation)").unwrap();
        for (_, v) in self.vertices.iter() {
            let [x, y, z] = v.pos.as_array();
            writeln!(out, "v {} {} {}", x, y, z).unwrap();
        }
        for (_, tri) in self.triangles.iter() {
            let [a, b, c] = tri.vertices;
            writeln!(
                out,
                "f {} {} {}",
                a.index() + 1,
                b.index() + 1,
                c.index() + 1
            )
            .unwrap();
        }
        out
    }
}

/// Tessellate a model into a provenance-tagged, crack-free triangle mesh.
///
/// **Specialized, not general**: planar faces are fan-triangulated (convex
/// only), cylindrical faces are sampled as a ruled band between their two
/// circular rims. Edge polylines are sampled once and shared, so adjacent faces
/// meet watertight (design §5). Assumes `Orientation::Forward` faces (every
/// current producer) — the loop winding already points outward.
///
/// Only cells **reachable from `live_solids`** are meshed. `Store` is append-only
/// and `boolean`/`pocket`/`pad` supersede rather than delete, so iterating the
/// stores would mesh the operands alongside the result.
pub fn tessellate(model: &Model, cfg: &TessConfig) -> Result<Tessellation, TessError> {
    let mut t = Tessellation::default();
    let mut vmap: HashMap<Handle<Vertex>, Handle<TessVertex>> = HashMap::new();
    // `Reachable` is a `HashSet`; walk the stores in their own order and merely ask
    // membership, so the mesh stays reproducible (design §2: replay).
    let reach = model.reachable();

    // 1. Sample every live edge into a shared polyline (the crack-free contract).
    for (eh, edge) in model.edges.iter() {
        if !reach.edges.contains(&eh) {
            continue;
        }
        let polyline = sample_edge(&mut t, &mut vmap, model, cfg, eh, edge);
        t.by_edge.insert(eh, polyline);
    }

    // 2. Triangulate each live face, reusing the shared edge polylines.
    for (fh, face) in model.faces.iter() {
        if !reach.faces.contains(&fh) {
            continue;
        }
        match model.surface(face.surface) {
            Surface::Plane(_) => triangulate_planar(&mut t, fh, face)?,
            Surface::Cylinder(cyl) => triangulate_cylinder(&mut t, model, fh, face, cyl),
        }
    }

    Ok(t)
}

/// A deduped mesh vertex for a topology vertex (tagged `OnVertex`).
fn vertex_of(
    t: &mut Tessellation,
    vmap: &mut HashMap<Handle<Vertex>, Handle<TessVertex>>,
    model: &Model,
    v: Handle<Vertex>,
) -> Handle<TessVertex> {
    if let Some(&h) = vmap.get(&v) {
        return h;
    }
    let pos = model.vertices.get(v).point;
    let h = t.vertices.push(TessVertex {
        pos,
        origin: TessOrigin::OnVertex(v),
    });
    vmap.insert(v, h);
    h
}

/// Sample one edge into a polyline of shared mesh vertices.
fn sample_edge(
    t: &mut Tessellation,
    vmap: &mut HashMap<Handle<Vertex>, Handle<TessVertex>>,
    model: &Model,
    cfg: &TessConfig,
    eh: Handle<Edge>,
    edge: &Edge,
) -> Vec<Handle<TessVertex>> {
    let [v0, v1] = edge
        .bounds
        .expect("M3 tess: bounded edges (a closed rim uses a seam vertex)");
    match model.edge_curve(eh) {
        Curve::Line(_) => vec![vertex_of(t, vmap, model, v0), vertex_of(t, vmap, model, v1)],
        Curve::Circle(c) => {
            // Full-circle rim (v0 == v1 = seam): point 0 is the seam vertex at
            // angle 0 (= centre + r·ref_dir), the rest are interior edge points.
            let n = circle_segments(cfg.tol, c.radius());
            let mut ring = Vec::with_capacity(n);
            ring.push(vertex_of(t, vmap, model, v0));
            for i in 1..n {
                let theta = std::f64::consts::TAU * (i as f64) / (n as f64);
                let h = t.vertices.push(TessVertex {
                    pos: c.point_at(theta),
                    origin: TessOrigin::OnEdge { edge: eh, t: theta },
                });
                ring.push(h);
            }
            ring
        }
    }
}

/// Segments to approximate a circle of `radius` within `tol` sagitta:
/// `r(1 − cos(π/n)) ≤ tol` ⟹ `n ≥ π / acos(1 − tol/r)`. Clamped for roundness
/// and against pathological inputs.
fn circle_segments(tol: f64, radius: f64) -> usize {
    const MIN: usize = 8;
    const MAX: usize = 4096;
    let ratio = tol / radius;
    if ratio <= 0.0 {
        return MAX; // tol ≤ 0: as fine as we allow
    }
    if ratio >= 2.0 {
        return MIN; // whole circle already within tol
    }
    let n = (std::f64::consts::PI / (1.0 - ratio).acos()).ceil();
    (n as usize).clamp(MIN, MAX)
}

/// Gather a loop's boundary as an ordered ring of shared mesh vertices,
/// concatenating each half-edge's oriented polyline and dropping consecutive
/// duplicates (shared endpoints) and the wrap-around duplicate.
fn boundary_ring(t: &Tessellation, lp: &Loop) -> Vec<Handle<TessVertex>> {
    let mut ring: Vec<Handle<TessVertex>> = Vec::new();
    for he in &lp.half_edges {
        let poly = &t.by_edge[&he.edge];
        // Append the oriented polyline, skipping a repeat of the previous vertex.
        let oriented: Box<dyn Iterator<Item = &Handle<TessVertex>>> = if he.forward {
            Box::new(poly.iter())
        } else {
            Box::new(poly.iter().rev())
        };
        for &h in oriented {
            if ring.last() != Some(&h) {
                ring.push(h);
            }
        }
    }
    if ring.len() > 1 && ring.first() == ring.last() {
        ring.pop();
    }
    ring
}

fn push_tri(t: &mut Tessellation, fh: Handle<Face>, vertices: [Handle<TessVertex>; 3]) {
    let th = t.triangles.push(TessTriangle { vertices, face: fh });
    t.by_face.entry(fh).or_default().push(th);
}

/// Triangulate a planar face: its outer ring, minus its holes.
///
/// The rings are shared edge polylines, and [`polygon::triangulate_polygon`] adds no
/// vertices — the sweep cuts along diagonals between existing ones, and the flip pass
/// only moves those — so adjacent faces still meet exactly (design §5).
fn triangulate_planar(
    t: &mut Tessellation,
    fh: Handle<Face>,
    face: &Face,
) -> Result<(), TessError> {
    let outer_h = boundary_ring(t, &face.outer);
    let holes_h: Vec<Vec<Handle<TessVertex>>> =
        face.inner.iter().map(|lp| boundary_ring(t, lp)).collect();

    // `triangulate_polygon` indexes a flat point slice; map the mesh handles onto one.
    let mut handles: Vec<Handle<TessVertex>> = outer_h.clone();
    handles.extend(holes_h.iter().flatten().copied());
    let pts: Vec<Point3> = handles.iter().map(|&h| t.vertices.get(h).pos).collect();
    let outer: Vec<usize> = (0..outer_h.len()).collect();
    let mut cut = outer_h.len();
    let holes: Vec<Vec<usize>> = holes_h
        .iter()
        .map(|h| {
            let r = (cut..cut + h.len()).collect();
            cut += h.len();
            r
        })
        .collect();
    let holes: Vec<&[usize]> = holes.iter().map(|h| h.as_slice()).collect();

    for tri in polygon::triangulate_polygon(&pts, &outer, &holes)? {
        push_tri(t, fh, tri.map(|i| handles[i]));
    }
    Ok(())
}

/// The centroid's projection onto `axis` (for ordering the two rims).
fn axial_centroid(t: &Tessellation, ring: &[Handle<TessVertex>], axis: nacre_math::Vector3) -> f64 {
    let sum: f64 = ring
        .iter()
        .map(|&h| (t.vertices.get(h).pos - Point3::origin()).dot(axis))
        .sum();
    sum / ring.len() as f64
}

/// Tessellate a cylindrical face as a ruled band between its two circular rims
/// (no interior samples — the surface is straight along the axis).
fn triangulate_cylinder(
    t: &mut Tessellation,
    model: &Model,
    fh: Handle<Face>,
    face: &Face,
    cyl: &Cylinder,
) {
    // The two distinct circular rim edges in the loop (seam lines excluded).
    let mut rim_edges: Vec<Handle<Edge>> = Vec::new();
    for he in &face.outer.half_edges {
        let is_circle = matches!(model.edge_curve(he.edge), Curve::Circle(_));
        if is_circle && !rim_edges.contains(&he.edge) {
            rim_edges.push(he.edge);
        }
    }
    debug_assert_eq!(rim_edges.len(), 2, "cylinder lateral face needs two rims");

    // Order by axial position so `a` is the lower rim → outward (radial) winding.
    let ring0 = t.by_edge[&rim_edges[0]].clone();
    let ring1 = t.by_edge[&rim_edges[1]].clone();
    let axis = cyl.axis().direction();
    let (a, b) = if axial_centroid(t, &ring0, axis) <= axial_centroid(t, &ring1, axis) {
        (ring0, ring1)
    } else {
        (ring1, ring0)
    };

    let n = a.len();
    debug_assert_eq!(b.len(), n, "rims sampled with matching segment counts");
    for i in 0..n {
        let j = (i + 1) % n;
        // Quad (a_i, a_j, b_j, b_i) split into two outward triangles.
        push_tri(t, fh, [a[i], a[j], b[j]]);
        push_tri(t, fh, [a[i], b[j], b[i]]);
    }
}

#[cfg(test)]
mod tests {
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
        let mut counts: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
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

    #[test]
    fn unit_cube_obj_shape() {
        let obj = to_obj(&cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).unwrap();
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

    #[test]
    fn vertices_round_trip() {
        let obj = to_obj(&cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0])).unwrap();
        let got: Vec<[f64; 3]> = v_lines(&obj)
            .into_iter()
            .map(|l| {
                let c = parse3(&l[2..]);
                [c[0], c[1], c[2]]
            })
            .collect();
        let expected = vec![
            [-2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 4.0, 0.0],
            [-2.0, 4.0, 0.0],
            [-2.0, 1.0, 10.0],
            [3.0, 1.0, 10.0],
            [3.0, 4.0, 10.0],
            [-2.0, 4.0, 10.0],
        ];
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
        let n = circle_segments(cfg.tol, 2.0);
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
                TessOrigin::OnVertex(v) => m.vertices.get(v).point,
                TessOrigin::OnEdge { edge, t: param } => match m.edge_curve(edge) {
                    Curve::Circle(c) => c.point_at(param),
                    Curve::Line(_) => unreachable!("cylinder rims are circles"),
                },
                TessOrigin::OnFace { .. } => unreachable!("cylinder uses no interior face samples"),
            };
            assert!((tv.pos - expected).norm() <= 1e-9);
        }
    }

    #[test]
    fn finer_tolerance_adds_triangles() {
        let m = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0, 5.0);
        let coarse = tessellate(&m, &TessConfig { tol: 0.5 })
            .unwrap()
            .triangles
            .len();
        let fine = tessellate(&m, &TessConfig { tol: 0.001 })
            .unwrap()
            .triangles
            .len();
        assert!(fine > coarse, "fine {fine} should exceed coarse {coarse}");
    }

    proptest! {
        #[test]
        fn arbitrary_box_obj_shape(
            min in prop::array::uniform3(-1e3f64..1e3),
            ext in prop::array::uniform3(1e-2f64..1e3),
        ) {
            let max = [min[0] + ext[0], min[1] + ext[1], min[2] + ext[2]];
            let obj = to_obj(&cube(min, max)).unwrap();
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
            let cfg = TessConfig::default();
            let n = circle_segments(cfg.tol, r);
            let t = tessellate(&cylinder([0.0, 0.0, 0.0], axis, r, h), &cfg).unwrap();
            prop_assert_eq!(t.vertices.len(), 2 * n);
            prop_assert_eq!(t.triangles.len(), 4 * n - 4);
            prop_assert_eq!(non_watertight_edges(&t), 0);
        }
    }
}
