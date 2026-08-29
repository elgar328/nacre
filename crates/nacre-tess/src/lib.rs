//! Mesh export for the nacre kernel.
//!
//! **One path**: [`tessellate`] → [`Tessellation`] (design §5). Edges are sampled once into shared
//! polylines and every face is triangulated in a chart of its own surface — a plane drops an axis,
//! a cylinder unrolls to `(z, r·θ)` — so adjacent faces meet crack-free and every mesh vertex
//! carries its origin ([`TessOrigin`]). OBJ text comes from [`Tessellation::to_obj`].
//!
//! ★★★★★ **There used to be two, and the second one lied.** An M1 bootstrap `to_obj(&Model)`
//! fan-triangulated planar faces **from half-edge start vertices**, which silently turns an arc
//! into a chord; this header promised it was "kept until the existing planar callers migrate"
//! (during M3). The migration landed in M6 instead, three milestones late, after that writer had
//! produced a wrong answer a third time — twice recorded in the dev-log and stepped in again. The
//! projection rule it carried (Newell → drop axis → repair handedness) lives once now, in
//! `planar_chart`; the sweep below it was always shared.

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
    /// **by index** or repeated within one, a spike, or two segments crossing with no
    /// vertex at the crossing.
    ///
    /// ★ **The cases where two vertices meet in *coordinates* moved out**, to
    /// [`Self::SelfTouchingBoundary`] — this doc used to list them here, and they are a different
    /// proposition: those rings are exactly what the b-rep asked for, and it is this decomposition
    /// that has no answer for them.
    ///
    /// **The sweep detects these, where ear clipping used to notice them by accident**
    /// (it stalled, and that stall was reported as `NoEar`). It is checked rather than
    /// assumed because the b-rep guarantees it and this layer cannot: a decomposition
    /// handed a self-crossing ring would otherwise return a confident, wrong mesh.
    DegenerateRing,
    /// **The face's boundary touches itself**: some vertex lies on the boundary somewhere other
    /// than at its own two edges — on another ring, or on a non-adjacent part of its own. The
    /// region is pinched there, so it is not a disk with sibling holes and this decomposition has
    /// no triangulation of it.
    ///
    /// ★★★★★ **This does not say the solid is wrong.** The measured population is an exact
    /// **tangency** — a hole touching another ring at one point, where a sampled circle vertex
    /// lands on the touch — and both fixtures that produce it are `validate`-clean with volumes
    /// exact to `1e-9`. Whether a face whose *interior* pinches is a non-manifold point on the
    /// surface is a real question, and `validate` cannot see it (there is no topology vertex
    /// there); it belongs to the capability that teaches the non-manifold test about tangential
    /// contact, not to this layer. So this name states only what this layer knows.
    ///
    /// ★ The **combinatorial** twin — two rings sharing an *index* — is refused one step earlier,
    /// by `monotone`'s `link`, and comes back as [`Self::DegenerateRing`]. That check has been
    /// there all along and calls the same thing a *pinch*; this is the spelling a chart actually
    /// produces, since `face_rings` gives every ring its own index range.
    SelfTouchingBoundary,
    /// A hole wound the same way as its outer ring: the b-rep does not keep
    /// material on every loop's left. A broken solid, not a repairable mesh.
    HoleWinding,
    /// A face's **interior** triangulation left an edge outside the declared budget — the mesh
    /// would misrepresent the surface there. See [`within_budget`]: the boundary is sampled to the
    /// budget by construction, so this only ever reports what the sweep chose, and it reports it
    /// rather than drawing a face that is quietly the wrong shape.
    OverBudget,
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
/// **One road for every face**: each is laid flat in a chart of its own surface (`Chart`) and
/// triangulated by the same sweep. Edge polylines are sampled once and shared, so adjacent faces
/// meet watertight (design §5). Reads the loop winding, not the `Orientation`
/// flag — every producer winds loops outward, `Reversed` faces included
/// (booleans emit both), and `validate` holds the two in agreement.
///
/// ★ This sentence used to read *"planar faces are fan-triangulated (convex only), cylindrical
/// faces are sampled as a ruled band between their two rims"* — both halves died with the chart
/// cell, and the doc outlived them.
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
        triangulate_face(&mut t, model, cfg, fh, face)?;
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

/// A face's boundary in the **chart its surface is meshed in**: `uv[i]` is where `handles[i]`
/// sits, and `rings` indexes into both — the outer ring first, then the holes.
///
/// ★★★★ **One road, one chart per surface.** Every face — planar or curved — is triangulated by
/// the same sweep ([`polygon::triangulate_uv`], design §5); what differs per surface is only how
/// its boundary is laid flat. A plane drops an axis; a cylinder unrolls to `(z, r·θ)`. When a
/// cone or a sphere arrives it adds a chart here and nothing else.
struct Chart {
    uv: Vec<polygon::P2>,
    handles: Vec<Handle<TessVertex>>,
    rings: Vec<Vec<usize>>,
    /// The way back — see [`ChartMap`].
    map: ChartMap,
    /// Points this chart offers the sweep's **interior**, in a fixed order. Empty for a plane,
    /// which needs none: a plane's chords lie on it.
    interior: Vec<polygon::P2>,
}

/// **A chart is a bijection, and this is the other direction** — from a point in the chart to the
/// surface parameters that name it and the position they evaluate to.
///
/// It exists because the chart is now allowed to *invent* a point ([`Chart::interior`]), and an
/// invented point has no model vertex to read a position off. Both arms therefore have to undo the
/// same handedness repair their forward direction applied, which is the one thing here that can be
/// silently wrong: the round-trip is measured on the boundary, where the forward answer is known.
enum ChartMap {
    /// Two world coordinates kept, one dropped — the dropped one is whatever the plane's own
    /// equation says it must be. `swapped` records the handedness repair.
    Plane {
        plane: nacre_geom::Plane,
        iu: usize,
        iv: usize,
        swapped: bool,
    },
    /// `(u, v) = (height, r·θ)`, with `negated` recording that the repair flipped the height.
    Cylinder { cyl: Cylinder, negated: bool },
}

impl ChartMap {
    /// The surface parameters at `uv`, and the point they name.
    ///
    /// For a cylinder the parameters are [`Cylinder::point_at`]'s own `(angle, height)`. For a
    /// plane they are the chart's two kept world coordinates — the plane has no parameterisation of
    /// its own, and the chart's is a perfectly good one.
    fn invert(&self, uv: polygon::P2) -> ([f64; 2], Point3) {
        match self {
            ChartMap::Plane {
                plane,
                iu,
                iv,
                swapped,
            } => {
                let (cu, cv) = if *swapped {
                    (uv[1], uv[0])
                } else {
                    (uv[0], uv[1])
                };
                let k = 3 - iu - iv;
                let n = plane.normal().as_array();
                let o = plane.origin().as_array();
                let mut c = [0.0; 3];
                c[*iu] = cu;
                c[*iv] = cv;
                // `n[k]` is the component the projection dropped, i.e. the largest — never zero.
                c[k] = o[k] - (n[*iu] * (cu - o[*iu]) + n[*iv] * (cv - o[*iv])) / n[k];
                ([cu, cv], Point3::from_array(c))
            }
            ChartMap::Cylinder { cyl, negated } => {
                let h = if *negated { -uv[0] } else { uv[0] };
                let theta = uv[1] / cyl.radius();
                ([theta, h], cyl.point_at(theta, h))
            }
        }
    }
}

/// The face's loops as shared mesh vertices, with `outer ++ holes` numbered in that order.
fn face_rings(t: &Tessellation, face: &Face) -> (Vec<Handle<TessVertex>>, Vec<Vec<usize>>) {
    let mut handles = boundary_ring(t, &face.outer);
    let mut rings = vec![(0..handles.len()).collect::<Vec<usize>>()];
    for lp in &face.inner {
        let h = boundary_ring(t, lp);
        rings.push((handles.len()..handles.len() + h.len()).collect());
        handles.extend(h);
    }
    (handles, rings)
}

/// A planar face's chart: drop the axis the face's own Newell normal is largest along.
///
/// The handedness repair is a **swap** of `u` and `v`: a plane's two axes are interchangeable, so
/// mirroring the frame costs nothing (a curved chart cannot do this — see [`cylinder_chart`]).
fn planar_chart(
    t: &Tessellation,
    face: &Face,
    plane: &nacre_geom::Plane,
) -> Result<Chart, TessError> {
    let (handles, rings) = face_rings(t, face);
    let pts: Vec<Point3> = handles.iter().map(|&h| t.vertices.get(h).pos).collect();
    let n = polygon::newell(&pts, &rings[0]);
    if n.norm() <= 0.0 {
        return Err(TessError::DegenerateRing);
    }
    let (iu, iv) = polygon::drop_axis(n);
    let mut uv: Vec<polygon::P2> = pts
        .iter()
        .map(|p| {
            let c = p.as_array();
            [c[iu], c[iv]]
        })
        .collect();
    let mut swapped = false;
    match polygon::ring_orientation(&rings[0], &uv) {
        1 => {}
        -1 => {
            for p in &mut uv {
                p.swap(0, 1);
            }
            swapped = true;
        }
        _ => return Err(TessError::DegenerateRing),
    }
    Ok(Chart {
        uv,
        handles,
        rings,
        map: ChartMap::Plane {
            plane: *plane,
            iu,
            iv,
            swapped,
        },
        interior: Vec::new(),
    })
}

/// A cylindrical face's chart: unroll to `(u, v) = (z, r·θ)` about the axis.
///
/// ★★ **`r·θ`, not `θ`.** A cylinder is developable, so scaling the angle by the radius makes
/// this an **isometry** — lengths and angles in the chart are lengths and angles on the surface.
/// That is what lets the Lawson pass ([`polygon`]'s flip) improve the mesh *on the surface* and
/// not merely in a distorted picture of it.
///
/// ★★ **`v` is the sweep axis, so the handedness repair negates `u` instead of swapping.**
/// Swapping would send the sweep along the axis instead of around it, and its diagonals would
/// then span wide arcs — chords that leave the surface. Negating `u` mirrors the frame just as
/// well and leaves the sweep going around.
///
/// ★ **θ is unwrapped along each loop, never taken absolutely.** A band's boundary walks its seam
/// **twice** — the same mesh vertices at `θ = 0` and at `θ = 2π` — and that is exactly what makes
/// the unrolled band a rectangle rather than a degenerate line. Holes are then shifted by whole
/// turns into the outer ring's range so they lie inside it.
fn cylinder_chart(
    t: &Tessellation,
    model: &Model,
    cfg: &TessConfig,
    face: &Face,
    cyl: &Cylinder,
) -> Result<Chart, TessError> {
    let (handles, rings) = face_rings(t, face);
    let axis = cyl.axis();
    let (o, w_dir) = (axis.origin(), axis.direction());
    let x_dir = cyl.ref_dir();
    let y_dir = w_dir.cross(x_dir);
    let r = cyl.radius();
    let mut uv: Vec<polygon::P2> = vec![[0.0, 0.0]; handles.len()];
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for ring in &rings {
        let mut prev_raw = 0.0;
        let mut theta = 0.0;
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for (k, &i) in ring.iter().enumerate() {
            let w = t.vertices.get(handles[i]).pos - o;
            let raw = w.dot(y_dir).atan2(w.dot(x_dir));
            if k == 0 {
                theta = raw;
            } else {
                let mut step = raw - prev_raw;
                while step > std::f64::consts::PI {
                    step -= std::f64::consts::TAU;
                }
                while step <= -std::f64::consts::PI {
                    step += std::f64::consts::TAU;
                }
                theta += step;
            }
            prev_raw = raw;
            uv[i] = [w.dot(w_dir), r * theta];
            lo = lo.min(uv[i][1]);
            hi = hi.max(uv[i][1]);
        }
        spans.push((lo, hi));
    }
    // Holes ride whole turns into the outer ring's window — an inner loop unwrapped on its own
    // may land a turn away from the material it is a hole in.
    let turn = std::f64::consts::TAU * r;
    for (ri, ring) in rings.iter().enumerate().skip(1) {
        let shift = ((spans[0].0 - spans[ri].0) / turn).ceil() * turn;
        if shift != 0.0 {
            for &i in ring {
                uv[i][1] += shift;
            }
        }
    }
    let mut negated = false;
    match polygon::ring_orientation(&rings[0], &uv) {
        1 => {}
        -1 => {
            for p in &mut uv {
                p[0] = -p[0];
            }
            negated = true;
        }
        _ => return Err(TessError::DegenerateRing),
    }
    let interior = interior_nodes(model, cfg, face, cyl, &uv, negated);
    Ok(Chart {
        uv,
        handles,
        rings,
        map: ChartMap::Cylinder { cyl: *cyl, negated },
        interior,
    })
}

/// ★★★★★ **The points a cylindrical face's boundary does not supply, and the mesh needs.**
///
/// A band's rims sample θ finely, and every diagonal the sweep draws between them lands on a
/// closely-sampled arc, so the interior was fine for free. It stops being free the moment a face's
/// boundary stops covering its θ range with arcs — a merged lateral face, whose erased phantom
/// seams had been the only arcs over one stretch, was meshed with chords spanning half a turn.
/// Nothing about the sweep was wrong there; it had no points to work with.
///
/// So the chart — the layer that knows the surface — supplies them, on a **lattice with no free
/// parameter in it**:
///
/// * **around** (`v`): the sampler's own step, `2πr / circle_segments`. Taken from that function
///   rather than re-derived from the budget, so the lattice lines up with every arc the boundary
///   already laid down instead of landing an ulp beside them. `i` covers the chart's whole `v`
///   range, because a chart's θ is *unwrapped* — a hole rides whole turns from the outer ring, and
///   a single turn's worth of lattice would miss it entirely.
/// * **along** (`u`): the axial coordinate of each boundary **circle's centre**. Read from the
///   curve, never from the sampled points: one rim's points agree on that coordinate mathematically
///   and differ by ulps in `f64`, so deduplicating *those* would turn one rim into 181 stations.
///   ★ Two arcs of one circle share a centre **bitwise**, and that is a guarantee rather than luck:
///   [`Model::push_edge`] derives the curve from the **canonicalized carrier pair alone** — the
///   circle arm never reads the endpoints — so two edges cut from the same circle by the same two
///   surfaces are handed identical inputs. The station count is therefore the number of planes the
///   face actually crosses, which is why `dedup` may compare `f64` with `==` here.
///
/// ★★ **The lowest and highest stations are dropped** — they are the face's own two rims (a
/// ruling's ends lie on arcs, so the arcs bound the `u` range), and a point there is a 1-ulp
/// duplicate of a boundary point rather than a new one. Measured: without this, an ordinary band
/// gained vertices and its triangle count moved.
fn interior_nodes(
    model: &Model,
    cfg: &TessConfig,
    face: &Face,
    cyl: &Cylinder,
    uv: &[polygon::P2],
    negated: bool,
) -> Vec<polygon::P2> {
    let (o, w_dir) = (cyl.axis().origin(), cyl.axis().direction());
    let mut stations: Vec<f64> = Vec::new();
    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
        for he in &lp.half_edges {
            if let Curve::Circle(c) = model.edge_curve(he.edge) {
                stations.push((c.center() - o).dot(w_dir));
            }
        }
    }
    stations.sort_by(f64::total_cmp);
    stations.dedup();
    if stations.len() <= 2 {
        return Vec::new();
    }
    let step = std::f64::consts::TAU * cyl.radius() / circle_segments(cfg, cyl.radius()) as f64;
    let (lo, hi) = uv.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
        (lo.min(p[1]), hi.max(p[1]))
    });
    let (i0, i1) = ((lo / step).ceil(), (hi / step).floor());
    let mut out = Vec::new();
    for &s in &stations[1..stations.len() - 1] {
        let u = if negated { -s } else { s };
        let mut i = i0;
        while i <= i1 {
            out.push([u, i * step]);
            i += 1.0;
        }
    }
    out
}

/// Triangulate one face: lay its boundary flat in the surface's chart, then run the one sweep.
fn triangulate_face(
    t: &mut Tessellation,
    model: &Model,
    cfg: &TessConfig,
    fh: Handle<Face>,
    face: &Face,
) -> Result<(), TessError> {
    let surface = model.surface(face.surface);
    let Chart {
        mut uv,
        mut handles,
        rings,
        map,
        interior,
    } = match surface {
        Surface::Plane(plane) => planar_chart(t, face, plane)?,
        Surface::Cylinder(cyl) => cylinder_chart(t, model, cfg, face, cyl)?,
    };
    let refs: Vec<&[usize]> = rings.iter().map(|r| r.as_slice()).collect();
    let boundary = handles.len();
    let tris = polygon::triangulate_uv(&mut uv, &refs, &interior)?;
    // The tail is exactly the candidates the sweep took, in the order it took them — the chart's
    // first minted vertices. Everything before it was already a shared boundary vertex.
    for &p in &uv[boundary..] {
        let (params, pos) = map.invert(p);
        handles.push(t.vertices.push(TessVertex {
            pos,
            origin: TessOrigin::OnFace {
                face: fh,
                uv: params,
            },
        }));
    }
    within_budget(t, cfg, surface, &handles, &rings, &tris)?;
    for tri in tris {
        push_tri(t, fh, tri.map(|i| handles[i]));
    }
    Ok(())
}

/// ★★★★★ **The mesh's one rule, asked of the face's *interior* — where it was never asked.**
///
/// `TessConfig` declares two budgets and [`circle_segments`] enforces **both** on every boundary
/// polyline. Nothing enforced them inside a face, which simply inherited whatever sampling its
/// boundary happened to supply: fine while the boundary samples the curvature (a band's rims do),
/// and silently wrong when it does not. A merged lateral face — whose erased phantom seams had
/// been the only thing sampling θ over one stretch — came out with four triangles spanning half a
/// turn as flat chords, losing a sixth of its area while `validate`, watertightness, the exact
/// volume and the face counts were all green.
///
/// So the same two questions are asked of every interior edge:
/// * **linear** — how far the surface strays from the chord, measured at its midpoint;
/// * **angular** — how far the surface *turns* between the ends, measured on its own normals.
///
/// Both come from [`nacre_geom::Surface`] and are exhaustive over the surface kinds, so a plane
/// answers zero to both without anyone declaring that a plane is flat, and a new surface kind is a
/// compile error rather than a missing case.
///
/// ★★ **A boundary edge cannot fail, so it is not asked.** An arc is cut into `circle_segments`
/// pieces, which is the *larger* of the two budgets' demands; a straight edge (a ruling, the seam,
/// a plane's side) strays zero and its ends share one normal. The two halves of the sentence meet
/// there — which is why this checks only what the sweep chose, never what the sampler laid down.
///
/// A violation is an error, not a repair: this layer's charter (see [`TessError`]) is that a wrong
/// cache is worse than none, and the defect above is exactly what a quiet one looks like.
fn within_budget(
    t: &Tessellation,
    cfg: &TessConfig,
    surface: &Surface,
    handles: &[Handle<TessVertex>],
    rings: &[Vec<usize>],
    tris: &[[usize; 3]],
) -> Result<(), TessError> {
    let boundary: std::collections::HashSet<(usize, usize)> = rings
        .iter()
        .flat_map(|r| {
            (0..r.len()).map(move |k| {
                let (a, b) = (r[k], r[(k + 1) % r.len()]);
                (a.min(b), a.max(b))
            })
        })
        .collect();
    // Today's populations sit **exactly** on the angular budget — an interior diagonal spanning one
    // sample turns by precisely the sampler's own step — so a bare `>` would be decided by the last
    // bit of two differently-rounded routes to the same number.
    let slack = 1.0 + 1e-9;
    let max_turn = cfg.max_angle_deg.to_radians() * slack;
    let max_sag = cfg.tol * slack;
    for tri in tris {
        for k in 0..3 {
            let (i, j) = (tri[k], tri[(k + 1) % 3]);
            if boundary.contains(&(i.min(j), i.max(j))) {
                continue;
            }
            let (a, b) = (
                t.vertices.get(handles[i]).pos,
                t.vertices.get(handles[j]).pos,
            );
            if surface.distance(a + (b - a) * 0.5) > max_sag {
                return Err(TessError::OverBudget);
            }
            let (Some(na), Some(nb)) = (surface.normal_at(a), surface.normal_at(b)) else {
                continue; // a normal with no name (a point on the axis) is not a budget failure
            };
            if na.dot(nb).clamp(-1.0, 1.0).acos() > max_turn {
                return Err(TessError::OverBudget);
            }
        }
    }
    Ok(())
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
    /// would pin an artifact. The old test could compare ordered because the bootstrap writer
    /// emitted the vertex *store*; that writer is gone.
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
                let f = m.faces.get(tri.face);
                let Surface::Cylinder(cy) = m.surface(f.surface) else {
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
        let face = m.faces.get(fh);
        match m.surface(face.surface) {
            Surface::Plane(p) => planar_chart(t, face, p).unwrap(),
            Surface::Cylinder(c) => cylinder_chart(t, m, cfg, face, c).unwrap(),
        }
    }

    /// ★★★★ **The chart follows the ring's own normal, on every plane — and the triangles follow
    /// the ring.**
    ///
    /// Two goldens used to say this one layer down, against `triangulate_polygon`: a square on
    /// `x = 5` wound CCW about `−x` (*"the projector must follow the ring's own normal, not a
    /// coordinate convention"*), and the same square handed both ways round (*"reversing the ring
    /// flips the Newell normal, and the same triangles come out with the opposite winding"*).
    /// That function is gone — the projection rule lives once now, in [`planar_chart`] — so the
    /// claims move here, where a cube states them **six times at once**, one per axis-aligned
    /// plane, with `drop_axis` returning each of its three answers.
    ///
    /// ★ Stronger than what it replaces: the old test compared a ring against its own reversal
    /// and could only say the two disagreed. This says which way is right — every face's
    /// triangles wind about the face's **outward** normal, which is the property the mesh's
    /// consumers actually rely on.
    #[test]
    fn every_planar_face_meshes_about_its_own_normal() {
        let cfg = TessConfig::default();
        let m = cube([-2.0, 1.0, 0.0], [3.0, 4.0, 10.0]);
        let t = tessellate(&m, &cfg).unwrap();
        let reach = m.reachable();
        let mut faces = 0;
        for (fh, face) in m.faces.iter() {
            if !reach.faces.contains(&fh) {
                continue;
            }
            let Surface::Plane(plane) = m.surface(face.surface) else {
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
            for (fh, _) in m.faces.iter() {
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
                    if let Surface::Cylinder(c) = m.surface(m.faces.get(fh).surface) {
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
            for (fh, face) in m.faces.iter() {
                if !matches!(m.surface(face.surface), Surface::Cylinder(_)) {
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
}
