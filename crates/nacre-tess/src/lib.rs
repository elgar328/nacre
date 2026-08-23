//! Mesh export for the nacre kernel.
//!
//! Two paths coexist during M3:
//! - [`to_obj`] — the M1 **bootstrap**: fan-triangulates a model's planar faces
//!   from edge endpoints (convex, planar only; panics on closed edges). Kept for
//!   the existing planar callers until they migrate.
//! - [`tessellate`] → [`Tessellation`] — the provenance-tagged layer (design §5):
//!   samples edges into shared polylines and triangulates faces (planar fans +
//!   ruled cylinder bands) crack-free, tagging every vertex with its origin.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
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
/// Every edge is bounded by type (S8); a boundless standalone circle would
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
    let pts: Vec<Point3> = model
        .vertices
        .iter()
        .map(|(vh, _)| model.vertex_point(vh))
        .collect();
    for p in &pts {
        let [x, y, z] = p.as_array();
        writeln!(out, "v {} {} {}", x, y, z).unwrap();
    }

    /// A loop's start vertices, as indices into the vertex store.
    fn ring(model: &Model, lp: &Loop) -> Vec<usize> {
        lp.half_edges
            .iter()
            .map(|he| {
                let bounds = model.edges.get(he.edge).vertices;
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

/// How finely a curve may be sampled — **two budgets, because they fail
/// differently.**
///
/// * `tol` is an **absolute** chord deviation (sagitta): the sampled polyline may
///   not stray further than this from the exact curve, in model units.
/// * `max_angle_deg` is a **relative** one: consecutive chords may not turn by more
///   than this. On a circle the two are one question asked twice — the sagitta is
///   `r(1 − cos(θ/2))` for a turn of `θ` — so the angular budget is the same rule
///   with the radius divided out, which is what makes it **scale-invariant**.
///
/// A single absolute budget is not "by curvature": it yields `n ≈ π√(r/2·tol)`, so a
/// small circle gets *fewer* segments (its absolute error was small to begin with) and
/// comes out visibly polygonal — a radius of `0.2` used to be a decagon. The angular
/// budget is what a small circle needs; the absolute one is what a circle much larger
/// than the tolerance needs. Neither substitutes for the other, so the segment count
/// takes whichever asks for more.
#[derive(Clone, Copy, Debug)]
pub struct TessConfig {
    pub tol: f64,
    /// Maximum turn between consecutive chords, in degrees.
    ///
    /// On a regular polygon this one number is three at once: the central angle of a
    /// segment, the kink at each vertex, and the angle between adjacent facet normals
    /// — so it bounds both how polygonal a silhouette looks and how banded its shading
    /// is. `2°` puts a circle that fills an 800px viewport within 0.06px of round.
    pub max_angle_deg: f64,
}

impl Default for TessConfig {
    fn default() -> Self {
        Self {
            tol: 1e-2,
            max_angle_deg: 2.0,
        }
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
/// meet watertight (design §5). Reads the loop winding, not the `Orientation`
/// flag — every producer winds loops outward, `Reversed` faces included
/// (booleans emit both), and `validate` holds the two in agreement.
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
    let pos = model.vertex_point(v);
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
    let [v0, v1] = edge.vertices;
    match model.edge_curve(eh) {
        Curve::Line(_) => vec![vertex_of(t, vmap, model, v0), vertex_of(t, vmap, model, v1)],
        Curve::Circle(c) => {
            if v0 == v1 {
                // Full-circle rim (v0 == v1 = seam): point 0 is the seam vertex at
                // angle 0 (= centre + r·ref_dir), the rest are interior edge points.
                let n = circle_segments(cfg, c.radius());
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
            } else {
                // ★ **An arc** (M6-2b): the stored `[v0, v1]` order is CCW about the axis — the
                // convention `derive_edge_curve`'s circle arm states — so the polyline walks
                // θ(v0) → θ(v0) + Δθ with `Circle::angle_of` as the one spelling of θ. The
                // segment count is the full circle's budget scaled by the arc's fraction (at
                // least one), and both endpoints are shared mesh vertices (crack-free with the
                // neighbouring faces' rings).
                let t0 = c.angle_of(model.vertex_point(v0));
                let dt =
                    (c.angle_of(model.vertex_point(v1)) - t0).rem_euclid(std::f64::consts::TAU);
                let n_full = circle_segments(cfg, c.radius()) as f64;
                let n = ((n_full * dt / std::f64::consts::TAU).ceil() as usize).max(1);
                let mut poly = Vec::with_capacity(n + 1);
                poly.push(vertex_of(t, vmap, model, v0));
                for i in 1..n {
                    let theta = t0 + dt * (i as f64) / (n as f64);
                    let h = t.vertices.push(TessVertex {
                        pos: c.point_at(theta),
                        origin: TessOrigin::OnEdge { edge: eh, t: theta },
                    });
                    poly.push(h);
                }
                poly.push(vertex_of(t, vmap, model, v1));
                poly
            }
        }
    }
}

/// Segments for a circle of `radius`: **whichever of the two budgets asks for more**
/// ([`TessConfig`] says why there are two).
///
/// * sagitta: `r(1 − cos(π/n)) ≤ tol` ⟹ `n ≥ π / acos(1 − tol/r)` — grows as `√r`, so
///   it is the term that speaks for a circle large in model units;
/// * angle: `360°/n ≤ Δθ` ⟹ `n ≥ 360/Δθ` — **the same for every radius**, so it is the
///   term that keeps a small circle from being an octagon.
///
/// With the default `2°` the angular term leads until `r ≈ 66` (that is `tol` divided by
/// `1 − cos(1°)`), which covers ordinary part sizes; past it the sagitta term takes over.
///
/// Pathological inputs are clamped rather than refused — a non-positive `tol` asks for
/// the finest we allow, a non-positive or absurd angle falls back on the other term and
/// the `[MIN, MAX]` clamp. The mesh is a cache; it answers rather than argues.
fn circle_segments(cfg: &TessConfig, radius: f64) -> usize {
    const MIN: usize = 8;
    const MAX: usize = 4096;
    let by_sagitta = {
        let ratio = cfg.tol / radius;
        if ratio <= 0.0 {
            MAX // tol ≤ 0: as fine as we allow
        } else if ratio >= 2.0 {
            MIN // whole circle already within tol
        } else {
            let n = (std::f64::consts::PI / (1.0 - ratio).acos()).ceil();
            (n as usize).clamp(MIN, MAX)
        }
    };
    let by_angle = {
        let n = (360.0 / cfg.max_angle_deg).ceil();
        if n.is_finite() && n >= 0.0 {
            (n as usize).clamp(MIN, MAX)
        } else {
            MIN
        }
    };
    by_sagitta.max(by_angle)
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

/// Tessellate a cylindrical face as a ruled band between its two rims — each either a single
/// closed circle edge or (M6-2b) a **chain of arc edges** — by a θ-merge walk over the rims'
/// shared polylines (no interior samples — the surface is straight along the axis).
///
/// ★ The walk generalizes the old equal-count quad pairing: with two closed rims the two rings
/// carry the same θs, every step is a tie, and the tie rule (advance the upper rim first)
/// reproduces the old quads' triangle set and windings — only the emission order inside each
/// quad swaps. With a chained rim the counts differ and the merge simply spends whichever ring's
/// next θ comes sooner; `n + m` triangles either way, every ring vertex consumed (crack-free
/// with the caps, which read the same polylines).
fn triangulate_cylinder(
    t: &mut Tessellation,
    model: &Model,
    fh: Handle<Face>,
    face: &Face,
    cyl: &Cylinder,
) {
    let axis = cyl.axis().direction();
    // Group the loop's circle-curve edges into the two rims by their circles' axial station.
    let mut rims: Vec<(f64, nacre_geom::Circle, Vec<Handle<Edge>>)> = Vec::new();
    for he in &face.outer.half_edges {
        let Curve::Circle(c) = model.edge_curve(he.edge) else {
            continue; // the seam (and, chained, nothing else) is straight
        };
        let key = (c.center() - Point3::origin()).dot(axis);
        match rims.iter_mut().find(|(k, ..)| (*k - key).abs() < 1e-9) {
            Some((.., edges)) => {
                if !edges.contains(&he.edge) {
                    edges.push(he.edge);
                }
            }
            None => rims.push((key, *c, vec![he.edge])),
        }
    }
    debug_assert_eq!(rims.len(), 2, "cylinder lateral face needs two rims");
    rims.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite axial stations"));

    // One rim as a CCW ring with an ascending θ per vertex, anchored at its θ-minimal vertex.
    // A near-seam angle that rounds to just under τ is folded to just under 0 first, so the
    // seam vertex anchors the ring whichever side of θ = 0 it realized on.
    let ring_of = |t: &Tessellation, c: &nacre_geom::Circle, edges: &[Handle<Edge>]| {
        let mut ring: Vec<Handle<TessVertex>> = Vec::new();
        if edges.len() == 1
            && model.edges.get(edges[0]).vertices[0] == model.edges.get(edges[0]).vertices[1]
        {
            ring = t.by_edge[&edges[0]].clone(); // a closed rim: the polyline is the ring
        } else {
            // Chain the arc polylines end-to-start by shared mesh vertices (the welding
            // guarantees each junction is one handle).
            let mut by_start: HashMap<Handle<TessVertex>, &Vec<Handle<TessVertex>>> =
                HashMap::new();
            for e in edges {
                let poly = &t.by_edge[e];
                by_start.insert(poly[0], poly);
            }
            let mut cur = *by_start.keys().next().expect("a rim has arcs");
            // Deterministic start: the handle-minimal polyline start.
            for &s in by_start.keys() {
                if s.index() < cur.index() {
                    cur = s;
                }
            }
            let chain_start = cur;
            for _ in 0..edges.len() {
                let poly = by_start[&cur];
                ring.extend(&poly[..poly.len() - 1]);
                cur = *poly.last().expect("arc polylines have two ends");
            }
            debug_assert_eq!(
                cur, chain_start,
                "the rim's arcs chain into one closed ring"
            );
        }
        let mut with_theta: Vec<(f64, Handle<TessVertex>)> = ring
            .iter()
            .map(|&h| {
                let raw = c.angle_of(t.vertices.get(h).pos);
                let th = if raw > std::f64::consts::TAU - 1e-9 {
                    raw - std::f64::consts::TAU
                } else {
                    raw
                };
                (th, h)
            })
            .collect();
        // Rotate to the θ-minimal vertex, then unwrap so θ ascends along the ring.
        let start = with_theta
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.0.partial_cmp(&b.1.0).expect("finite angles"))
            .map(|(i, _)| i)
            .expect("a rim has vertices");
        with_theta.rotate_left(start);
        for i in 1..with_theta.len() {
            if with_theta[i].0 < with_theta[i - 1].0 {
                with_theta[i].0 += std::f64::consts::TAU;
            }
        }
        with_theta
    };
    let (_, ref c_lo, ref e_lo) = rims[0];
    let (_, ref c_hi, ref e_hi) = rims[1];
    let a = ring_of(t, c_lo, e_lo);
    let b = ring_of(t, c_hi, e_hi);

    // The merge: from corner (aᵢ, bⱼ), spend whichever ring's next vertex comes first in θ
    // (the upper rim on a tie — the old quad split's diagonal), wrapping each ring once.
    let (n, m) = (a.len(), b.len());
    let theta_at = |ring: &[(f64, Handle<TessVertex>)], k: usize| {
        let (th, _) = ring[k % ring.len()];
        th + if k >= ring.len() {
            std::f64::consts::TAU
        } else {
            0.0
        }
    };
    let (mut i, mut j) = (0usize, 0usize);
    while i < n || j < m {
        let ai = a[i % n].1;
        let bj = b[j % m].1;
        let advance_b = if i >= n {
            true
        } else if j >= m {
            false
        } else {
            theta_at(&b, j + 1) <= theta_at(&a, i + 1)
        };
        if advance_b {
            push_tri(t, fh, [ai, b[(j + 1) % m].1, bj]);
            j += 1;
        } else {
            push_tri(t, fh, [ai, a[(i + 1) % n].1, bj]);
            i += 1;
        }
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
    /// It used to be 2, and with the angular budget in place both tolerances would
    /// answer the same 180 segments there — the test would have compared a number
    /// with itself and passed on any `tol` at all. The angular term leads until
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
            for (_, ring) in t.by_edge.iter() {
                if ring.len() < 3 {
                    continue; // a straight seam edge has no turn to measure
                }
                // Wrapping is right because every multi-point polyline here is a closed
                // rim; an arc would need the two end turns left out (M6-3).
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

    /// ★ **Scale invariance**: the angular budget is a *relative* one, so a small
    /// circle and a large one are cut into the same number of pieces. Under the old
    /// absolute-only rule these were 10 and 100 — the very asymmetry that made small
    /// holes look like polygons.
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
}
