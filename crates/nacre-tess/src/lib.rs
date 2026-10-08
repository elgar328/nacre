//! Mesh export for the nacre kernel.
//!
//! **One path**: [`tessellate`] (every live solid) or [`tessellate_solids`] (the solids given) →
//! [`Tessellation`]. Edges are sampled once into shared
//! polylines and every face is triangulated in a chart of its own surface — a plane drops an axis,
//! a cylinder unrolls to `(z, r·θ)` — so adjacent faces meet crack-free and every mesh vertex
//! carries its origin ([`TessOrigin`]). OBJ text comes from [`Tessellation::to_obj`].
//!
//! ★★★★★ **One, and only one.** A second writer that fan-triangulates planar faces **from
//! half-edge start vertices** silently turns an arc into a chord. The projection rule
//! (Newell → drop axis → repair handedness) lives once, in `planar_chart`; the sweep below it
//! is shared.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod polygon;

use nacre_geom::{Circle, Curve, Cylinder, Surface};
use nacre_math::Point3;
use nacre_store::{Handle, Store};
use nacre_topo::{Edge, Face, Loop, Model, Solid, Vertex};
use std::collections::HashMap;
use std::fmt::Write;

/// A face this tessellator cannot mesh. Returned rather than approximated: the
/// mesh is a cache, but a *wrong* cache is worse than none — silently fanning a
/// ring that is not star-shaped puts triangles outside the solid.
///
/// There is deliberately no `Fallback` variant. A fan would be one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TessError {
    /// The rings are not a polygon with sibling holes, so no triangulation of them
    /// exists: fewer than three vertices, a zero-area ring, a vertex used by two rings
    /// **by index** or repeated within one, a spike, **two segments passing through each
    /// other**, or a ring passing **through** a point it visits twice (`polygon::pinch`). On a
    /// cylinder also a face whose loops wrap the axis other than as two rims, or a
    /// band no generator cuts so that it stays one polygon — every generator crosses some hole
    /// other than twice (`cut_seamless_bands`).
    ///
    /// ★ **The cases where two vertices meet in *coordinates* are not here** but
    /// [`Self::SelfTouchingBoundary`] — a different proposition: those rings are exactly what the
    /// b-rep asked for, and it is this decomposition that has no answer for them.
    ///
    /// ☑ **The crossing clause covers a hole crossing its outer ring**, not only a ring crossing
    /// itself: without `monotone`'s `self_touch` that input comes back `Ok` with eight confident,
    /// wrong triangles (measured).
    ///
    /// **The sweep detects these** — ear clipping would notice them only by accident, as a
    /// stall. It is checked rather than
    /// assumed because the b-rep guarantees it and this layer cannot: a decomposition
    /// handed a self-crossing ring would otherwise return a confident, wrong mesh.
    DegenerateRing,
    /// **The face's boundary touches itself, and no bridge was laid there**: some vertex lies on
    /// the boundary somewhere other than at its own two edges — on another ring, or on a
    /// non-adjacent part of its own. The region is pinched there, so it is not a disk with sibling
    /// holes and this decomposition has no triangulation of it *as two rings*.
    ///
    /// ★★★★★ **The common case is drawn, not refused.** A hole whose sample sits
    /// strictly inside a straight edge shared by two planar faces (`self_touch`'s
    /// `Interior` + `Tangent`) has that sample put into the edge's polyline by [`tessellate`]'s
    /// bridge pre-pass, the two rings are spliced into one at that point, and the sweep orders the
    /// resulting coincident pair symbolically (`polygon::sos`). **One ring passing one mesh vertex
    /// twice** — a groove whose tip touches the face's own rim — is drawn the same way, its two
    /// visits the twins, when both visits' wedges are convex (`polygon::pinch`). What still comes
    /// back under this name is what those roads do not cover, the bridge's each counted by the
    /// pre-pass: a vertex on another ring's *vertex* (`AtEnd`), a neighbour exactly on the touched
    /// line, a touch on a curved edge or beside a curved face, two touches on one segment or at
    /// one point, a pinch with a reflex wedge, and a face with more than one bridge or pinch.
    ///
    /// ★★★★★ **This does not say the solid is wrong — and that is now measured, not hoped.** The
    /// population is an exact **tangency**: a hole touching another ring at one point, where a
    /// sampled circle vertex lands on the touch. Both fixtures are `validate`-clean with volumes
    /// exact to `1e-9`, and
    /// **the surface is a 2-manifold at the touch**: the link of the boundary on a small sphere
    /// there is a *single* circle, because the pinched face's two lobes are joined around through
    /// the neighbouring curved face; OCCT, given the same operands, returns a body of the same
    /// volume, area and face count. What is pinched is the **face**, not the surface, so the
    /// capability that owes an answer is **not** the non-manifold test (it was right) but this
    /// layer's own: a consistent symbolic order for coincident vertices in the sweep — which
    /// `polygon::sos` supplies. ★ For what is still refused, the
    /// loss is larger than one face: a mesh stops at the first refusal in its range ([`tessellate`]
    /// — every live solid — or [`tessellate_solids`]), and the playground's drawing stops at the
    /// first value whose mesh is refused, so one pinched face erases every body of that run.
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
// Provenance-tagged tessellation
// ---------------------------------------------------------------------------

/// Where a tessellation vertex came from — the truth it can snap back to.
/// `t`/`uv` are the exact curve/surface parameters.
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

/// A provenance-tagged triangle mesh derived from a [`Model`].
///
/// `by_edge` holds each edge's shared polyline (the crack-free contract: faces
/// consume these, never re-sample), `by_face` the triangles per face. (An
/// incremental `stale` set is deferred — this is from-scratch, like
/// `Adjacency`.)
///
/// ★★★★★ **A closed edge's polyline does not repeat its first point — the closure is
/// implicit.** [`sample_edge`]'s full-rim arm samples a circle at `0 .. (n−1)τ/n` and stops,
/// and `boundary_ring` drops the wrap-around duplicate for the same reason: these are mesh
/// vertices, and repeating a handle would give the face a degenerate triangle. So a consumer
/// that walks a polyline **pairwise** gets `n − 1` steps for a ring of `n`, and must add the
/// closing step itself; `edge.vertices[0] == edge.vertices[1]` is the fact to branch on — the
/// very test `sample_edge` uses. ⚠ Written here because a consumer that does not know it draws
/// every uncut rim with a gap in it (measured in the playground's viewport).
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
/// comes out visibly polygonal — a radius of `0.2` is a decagon under that budget alone. The
/// angular
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
    /// Wavefront OBJ text of the tessellation: vertices in store order, triangles in store order.
    /// What it holds is decided when it is built — every live solid ([`tessellate`]) or the ones
    /// given ([`tessellate_solids`]).
    ///
    /// ★ **Each corner carries its face's outward normal** (`vn`, `f a//na b//nb c//nc`) —
    /// [`nacre_props::face_normal_at`], the one the playground's viewport lights with. Without it a
    /// viewer either shades each triangle flat, and a cylinder reads as bands, or averages the
    /// normals of a shared vertex — and a rim's points *are* shared, cap and side meeting on the
    /// one polyline the edge was sampled into, so the sharp rim would be rounded off. With it a
    /// rim point is written once and paired with two normals: the cap's axis and the side's radius.
    ///
    /// `model` must be the model this tessellation was built from — the mesh holds only face
    /// handles, so another model would answer for other faces (or none, and panic).
    pub fn to_obj(&self, model: &Model) -> String {
        let vertices: Vec<Handle<TessVertex>> = self.vertices.iter().map(|(h, _)| h).collect();
        let triangles: Vec<Handle<TessTriangle>> = self.triangles.iter().map(|(h, _)| h).collect();
        self.write_obj(model, &vertices, &triangles)
    }

    /// The one OBJ writer: `v` lines for `vertices` (numbered 1.. in that order), then for each
    /// triangle one `vn` per corner not written before (equal bits, one line) and its `f` line.
    /// Every triangle's corners must be among `vertices`.
    fn write_obj(
        &self,
        model: &Model,
        vertices: &[Handle<TessVertex>],
        triangles: &[Handle<TessTriangle>],
    ) -> String {
        let mut out = String::new();
        writeln!(out, "# nacre OBJ export (provenance tessellation)").unwrap();
        let mut number: HashMap<Handle<TessVertex>, usize> = HashMap::new();
        for (i, &vh) in vertices.iter().enumerate() {
            let [x, y, z] = self.vertices.get(vh).pos.as_array();
            writeln!(out, "v {} {} {}", x, y, z).unwrap();
            number.insert(vh, i + 1);
        }
        let mut normals: HashMap<[u64; 3], usize> = HashMap::new();
        for &th in triangles {
            let tri = self.triangles.get(th);
            let p = tri.vertices.map(|vh| self.vertices.get(vh).pos);
            let mut corners = [(0, 0); 3];
            for (k, &vh) in tri.vertices.iter().enumerate() {
                let n = nacre_props::face_normal_at(model, tri.face, p[k])
                    .map(|v| v.as_array())
                    .unwrap_or_else(|| facet_normal(p));
                let next = normals.len() + 1;
                let ni = *normals.entry(n.map(f64::to_bits)).or_insert_with(|| {
                    writeln!(out, "vn {} {} {}", n[0], n[1], n[2]).unwrap();
                    next
                });
                corners[k] = (number[&vh], ni);
            }
            let [(a, na), (b, nb), (c, nc)] = corners;
            writeln!(out, "f {a}//{na} {b}//{nb} {c}//{nc}").unwrap();
        }
        out
    }
}

/// The plane of a triangle, unit — the normal for a corner whose face names no direction there (a
/// point on a cylinder's own axis, which the mesh does not produce); `+z` for a degenerate one.
fn facet_normal(p: [Point3; 3]) -> [f64; 3] {
    let n = (p[1] - p[0]).cross(p[2] - p[0]);
    n.normalize().map_or([0.0, 0.0, 1.0], |u| u.as_array())
}

/// Tessellate a model into a provenance-tagged, crack-free triangle mesh.
///
/// **One road for every face**: each is laid flat in a chart of its own surface (`Chart`) and
/// triangulated by the same sweep. Edge polylines are sampled once and shared, so adjacent faces
/// meet watertight. Reads the loop winding, not the `Orientation`
/// flag — every producer winds loops outward, `Reversed` faces included
/// (booleans emit both), and `validate` holds the two in agreement.
///
/// Only cells **reachable from `live_solids`** are meshed. `Store` is append-only
/// and `boolean`/`pocket`/`pad` supersede rather than delete, so iterating the
/// stores would mesh the operands alongside the result.
pub fn tessellate(model: &Model, cfg: &TessConfig) -> Result<Tessellation, TessError> {
    tessellate_solids(model, model.live_solids(), cfg)
}

/// [`tessellate`] over the cells `solids` (live solids) reach — **what is shown, and nothing else.**
///
/// An application's model holds more live solids than it draws: a script's every intermediate
/// value stays live when the layer above copies before each consumption. Measured, a cylinder moved
/// 300 times keeps 301 live solids and shows one — meshing all of them took 529 ms (native release)
/// for one cylinder's triangles.
///
/// The same road, walked over less: edges are sampled and faces triangulated one by one, two
/// solids share no edge, and the bridge pre-pass pairs only faces on one edge — so a face's
/// triangles are the same whichever other solids are meshed beside it.
pub fn tessellate_solids(
    model: &Model,
    solids: &[Handle<Solid>],
    cfg: &TessConfig,
) -> Result<Tessellation, TessError> {
    let mut t = Tessellation::default();
    // `Reachable` is a `HashSet`; walk the stores in their own order and merely ask
    // membership, so the mesh stays reproducible (replay).
    let reach = model.reachable_from(solids);

    // 1. Sample every live edge into a shared polyline (the crack-free contract).
    let live = sample_live_edges(&mut t, model, cfg, &reach);

    // 1¼. A cylindrical face whose boundary wraps the axis — a band, bounded by its two rims —
    // gets one generator to be cut along, its points put into the crossed loops' polylines.
    // Before the bridge pre-pass, which reads the polylines as final. See [`cut_seamless_bands`].
    let cuts = cut_seamless_bands(&mut t, model, &live);

    // 1½. Where a hole's sample lands exactly on a straight shared edge (an exact tangency, the
    // one boundary the sweep cannot decompose), put that sample into the edge's polyline so
    // every face on the edge sees it. Nothing is minted and no coordinate moves: an existing
    // vertex becomes a sample of a second edge it lies on. See [`bridge_shared_edges`].
    let _report = bridge_shared_edges(&mut t, model, &live);
    // 2. Triangulate each live face, reusing the shared edge polylines.
    for &(fh, face) in &live {
        triangulate_face(&mut t, model, cfg, fh, face, cuts.get(&fh).cloned())?;
    }

    Ok(t)
}

/// Phase 1 of [`tessellate`]: sample every live edge into `by_edge`, and hand back the live
/// faces in store order so the later phases walk exactly that set.
fn sample_live_edges<'m>(
    t: &mut Tessellation,
    model: &'m Model,
    cfg: &TessConfig,
    reach: &nacre_topo::Reachable,
) -> Vec<(Handle<Face>, &'m Face)> {
    let mut vmap: HashMap<Handle<Vertex>, Handle<TessVertex>> = HashMap::new();
    let mut i = 0u32;
    while let Some(eh) = model.edge_handle_at(i) {
        i += 1;
        let edge = model.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        let polyline = sample_edge(t, &mut vmap, model, cfg, eh, edge);
        t.by_edge.insert(eh, polyline);
    }
    (0..model.face_count() as u32)
        .filter_map(|i| model.face_handle_at(i))
        .map(|h| (h, model.face(h)))
        .filter(|(fh, _)| reach.faces.contains(fh))
        .collect()
}

/// One insertion the bridge pre-pass made: `vertex` now sits at `at` in `by_edge[edge]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Split {
    pub edge: Handle<Edge>,
    pub at: usize,
    pub vertex: Handle<TessVertex>,
}

/// Why the bridge pre-pass left a bridgeable-looking touch alone. Counted, never acted on:
/// every case is a shape this pass does not bridge, and the sweep goes on refusing it by name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Declined {
    /// The touched segment is a chord of a curved edge: a point on the chord is not on the curve.
    CurvedEdge,
    /// A face sharing the touched edge is not planar: a vertex planted there would seed diagonals
    /// on a curved chart, a population this pass has not measured.
    CurvedNeighbour,
    /// Two touches on one segment: three coincident vertices have no consistent order.
    SharedSegment,
    /// One vertex touching two different segments.
    MultiSegment,
    /// The touched segment was not found consecutive in any polyline of its face.
    NoEdge,
}

/// What the bridge pre-pass did and declined, for the caller that wants to read it.
#[derive(Clone, Debug, Default)]
pub struct BridgeReport {
    pub splits: Vec<Split>,
    pub declined: Vec<(Handle<Face>, Declined)>,
}

/// ★★★★★ **The bridge pre-pass — a tangency's sample becomes a sample of the edge it touches.**
///
/// A planar face whose hole touches another ring at exactly one point is a valid solid the
/// sweep cannot decompose: the touching vertex is a sample of the *curved* edge, and the
/// straight edge it lands on carries only its two ends (`sample_edge`'s `Line` arm), so the
/// square ring has no vertex there to bridge to. Inserting that existing vertex into the
/// straight edge's polyline gives the face the twin it needs — and, because `by_edge` is
/// shared and read live by every face, gives the neighbouring face the same vertex, which is
/// what keeps the mesh crack-free (a split one side does not know about is a T-vertex).
///
/// **Order matters twice.** Touches are collected over *every* live face before any polyline
/// is changed, so a face walked later does not see an already-split edge and misreport the
/// touch as `AtEnd`; and insertions on one edge go later-position-first, so earlier positions
/// stay valid. The pass is idempotent: a polyline that already holds the vertex is left alone.
///
/// **Edge → faces is built here, from the live set**, not read from the model's adjacency —
/// inside a boolean the adjacency is stale (its callers rebuild it afterwards), and this pass
/// runs on the model as handed to it.
///
/// What it declines, it counts (see [`Declined`]); the population it was measured on is four
/// planar caps, every touch `Interior` and `Tangent` on a straight edge shared by two planes.
fn bridge_shared_edges(
    t: &mut Tessellation,
    model: &Model,
    live: &[(Handle<Face>, &Face)],
) -> BridgeReport {
    use polygon::{TouchKind, Witness};
    struct Cand {
        face: Handle<Face>,
        vertex: Handle<TessVertex>,
        seg: [Handle<TessVertex>; 2],
    }
    struct Pending {
        edge: Handle<Edge>,
        at: usize,
        vertex: Handle<TessVertex>,
        face: Handle<Face>,
    }
    let mut report = BridgeReport::default();

    // 1. Every bridgeable touch, read before anything is changed.
    let mut cands: Vec<Cand> = Vec::new();
    for &(fh, face) in live {
        let plane = match model.surface_cache(face.surface) {
            Surface::Plane(plane) => plane,
            // Bridges are defined on planar charts only. An arm rather than an `else`, so a third
            // surface kind is a compile error here and its author decides whether it bridges.
            Surface::Cylinder(_) => continue,
        };
        let Ok(chart) = planar_chart(t, face, plane) else {
            continue; // the face's own triangulation will name this error
        };
        let refs: Vec<&[usize]> = chart.rings.iter().map(|r| r.as_slice()).collect();
        let Ok(m) = polygon::meets(&chart.uv, &refs) else {
            continue;
        };
        let bridgeable =
            |tc: &polygon::Touch| tc.kind == TouchKind::Interior && tc.witness == Witness::Tangent;
        let mut per_vertex: HashMap<usize, usize> = HashMap::new();
        for tc in m.touches.iter().filter(|tc| bridgeable(tc)) {
            *per_vertex.entry(tc.vertex).or_default() += 1;
        }
        for tc in m.touches.iter().filter(|tc| bridgeable(tc)) {
            if per_vertex[&tc.vertex] > 1 {
                report.declined.push((fh, Declined::MultiSegment));
                continue;
            }
            cands.push(Cand {
                face: fh,
                vertex: chart.handles[tc.vertex],
                seg: [chart.handles[tc.segment], chart.handles[tc.segment_end]],
            });
        }
    }
    if cands.is_empty() {
        return report;
    }

    // 2. Which faces share each edge, and which of them are planar.
    let planar: HashMap<Handle<Face>, bool> = live
        .iter()
        .map(|&(fh, f)| {
            (
                fh,
                // ★ The **truth**, path-qualified: in this file the bare `Surface` is geom's
                // (the cache), and aliasing the topo one in would give the truth a second
                // vocabulary here. A face's kind is a fact about what it is, so the truth
                // answers it — and a new surface kind lands there first.
                matches!(model.surface(f.surface), nacre_topo::Surface::Plane { .. }),
            )
        })
        .collect();
    let mut faces_on: HashMap<Handle<Edge>, Vec<Handle<Face>>> = HashMap::new();
    for &(fh, face) in live {
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                faces_on.entry(he.edge).or_default().push(fh);
            }
        }
    }

    // 3. Find each touched segment in its face's polylines and apply the guards.
    let mut pending: Vec<Pending> = Vec::new();
    for c in &cands {
        let face = live
            .iter()
            .find(|(fh, _)| *fh == c.face)
            .map(|(_, f)| *f)
            .expect("a candidate names a live face");
        let mut found = None;
        'search: for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                let poly = &t.by_edge[&he.edge];
                for k in 0..poly.len().saturating_sub(1) {
                    let (p, q) = (poly[k], poly[k + 1]);
                    if (p == c.seg[0] && q == c.seg[1]) || (p == c.seg[1] && q == c.seg[0]) {
                        found = Some((he.edge, k));
                        break 'search;
                    }
                }
            }
        }
        let Some((eh, k)) = found else {
            report.declined.push((c.face, Declined::NoEdge));
            continue;
        };
        if !matches!(model.edge_curve(eh), Curve::Line(_)) {
            report.declined.push((c.face, Declined::CurvedEdge));
            continue;
        }
        if faces_on[&eh].iter().any(|f| !planar[f]) {
            report.declined.push((c.face, Declined::CurvedNeighbour));
            continue;
        }
        pending.push(Pending {
            edge: eh,
            at: k,
            vertex: c.vertex,
            face: c.face,
        });
    }
    let mut on_segment: HashMap<(Handle<Edge>, usize), usize> = HashMap::new();
    for p in &pending {
        *on_segment.entry((p.edge, p.at)).or_default() += 1;
    }
    let (shared, mut ok): (Vec<_>, Vec<_>) = pending
        .into_iter()
        .partition(|p| on_segment[&(p.edge, p.at)] > 1);
    for p in shared {
        report.declined.push((p.face, Declined::SharedSegment));
    }

    // 4. Insert — per edge, later positions first, never twice.
    ok.sort_by_key(|p| (p.edge.index(), std::cmp::Reverse(p.at)));
    for p in ok {
        let poly = t.by_edge.get_mut(&p.edge).expect("the edge was sampled");
        if poly.contains(&p.vertex) {
            continue;
        }
        poly.insert(p.at + 1, p.vertex);
        report.splits.push(Split {
            edge: p.edge,
            at: p.at + 1,
            vertex: p.vertex,
        });
    }
    report
}

/// The bridge pre-pass alone, for a test that wants to see what it did before any face is
/// triangulated (the returned `Tessellation` holds the sampled, split edge polylines and no
/// triangles). Test-only (`test-util`).
#[cfg(feature = "test-util")]
pub fn bridge_report(model: &Model, cfg: &TessConfig) -> (BridgeReport, Tessellation) {
    let mut t = Tessellation::default();
    let reach = model.reachable();
    let live = sample_live_edges(&mut t, model, cfg, &reach);
    let report = bridge_shared_edges(&mut t, model, &live);
    (report, t)
}

/// A triangle on a cylinder's lateral that leans off the surface, and by how much: the angle in
/// degrees between its facet normal and the surface normal at its centroid — `None` when it has no
/// area, so no facet normal. Test-only (`test-util`).
#[cfg(any(test, feature = "test-util"))]
#[derive(Clone, Copy, Debug)]
pub struct Leaning {
    pub triangle: Handle<TessTriangle>,
    pub face: Handle<Face>,
    pub degrees: Option<f64>,
}

/// Every lateral triangle of `t` whose facet departs from its surface by more than `limit_deg`, or
/// that has no area. Test-only (`test-util`): the mesh locks call it with the budget's
/// `max_angle_deg`.
///
/// ★ **Why the budget bounds it.** Each edge of a lateral triangle turns about the axis by at most
/// the budget — an arc's chord by [`circle_segments`], an interior chord by `within_budget`, a
/// ruling not at all — so its three corners lie within one budget's turn of each other, and so does
/// its facet's normal of the surface's. A sliver the chart made out of a lattice point standing a few
/// ulps off the boundary leaned 2.6° to 90° (at 90° it lies in a cap plane).
#[cfg(any(test, feature = "test-util"))]
pub fn leaning_laterals(model: &Model, t: &Tessellation, limit_deg: f64) -> Vec<Leaning> {
    let mut out = Vec::new();
    for (triangle, tri) in t.triangles.iter() {
        let surface = model.surface_cache(model.face(tri.face).surface);
        if !matches!(surface, Surface::Cylinder(_)) {
            continue;
        }
        let [a, b, c] = tri.vertices.map(|h| t.vertices.get(h).pos);
        let facet = (b - a).cross(c - a);
        let len = facet.norm();
        let centroid = a + ((b - a) + (c - a)) * (1.0 / 3.0);
        let degrees = (len > 0.0)
            .then(|| surface.normal_at(centroid))
            .flatten()
            .map(|n| {
                (facet * (1.0 / len))
                    .dot(n)
                    .abs()
                    .clamp(0.0, 1.0)
                    .acos()
                    .to_degrees()
            });
        if degrees.is_none_or(|d| d > limit_deg) {
            out.push(Leaning {
                triangle,
                face: tri.face,
                degrees,
            });
        }
    }
    out
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
                // A full circle: point 0 is its vertex, the rest are interior edge points going
                // round from it. A rim's vertex is its seam point (`OnSeam`), at angle 0
                // (= centre + r·ref_dir) by definition; any other vertex is read off the curve.
                let n = circle_segments(cfg, c.radius());
                let t0 = match model.vertex(v0) {
                    Vertex::OnSeam(_) => 0.0,
                    _ => c.angle_of(model.vertex_point(v0)),
                };
                let mut ring = Vec::with_capacity(n);
                ring.push(vertex_of(t, vmap, model, v0));
                for i in 1..n {
                    let theta = t0 + std::f64::consts::TAU * (i as f64) / (n as f64);
                    let h = t.vertices.push(TessVertex {
                        pos: c.point_at(theta),
                        origin: TessOrigin::OnEdge { edge: eh, t: theta },
                    });
                    ring.push(h);
                }
                ring
            } else {
                // ★ **An arc**: the stored `[v0, v1]` order is CCW about the axis — the
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
/// the same sweep ([`polygon::triangulate_uv`]); what differs per surface is only how
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
/// ★ **θ is unwrapped along each loop, never taken absolutely.** A band is bounded by its two rims,
/// each turning once around the axis, and unrolled alone they are two open lines; `cut` joins
/// them at the generator [`cut_seamless_bands`] chose — first rim, its cut point again, up the
/// generator through any hole it crosses, the second rim, its cut point again, and back down
/// ([`joined_band`]) — so each cut point stands at `θ` and at `θ + 2π` and the unrolled band is a
/// rectangle rather than a degenerate line. The other holes are then shifted by whole turns into
/// the outer ring's range so they lie inside it.
///
/// ★ **`u` of a circle's point is the circle's place**, [`station`], not the point's own axial
/// coordinate — so every arc is level in the chart, at the height the lattice lays its row
/// ([`interior_nodes`]). Only a point on no circle — measured: none, over the lib suite and the
/// census; a ruling broken at a vertex by another face's edge would be one — reads its own.
fn cylinder_chart(
    t: &Tessellation,
    model: &Model,
    cfg: &TessConfig,
    face: &Face,
    cyl: &Cylinder,
    cut: Option<BandCut>,
) -> Result<Chart, TessError> {
    let (mut handles, mut rings) = face_rings(t, face);
    if let Some(cut) = cut {
        let cut = cut.ok_or(TessError::DegenerateRing)?;
        (handles, rings) = joined_band(&handles, &rings, &cut).ok_or(TessError::DegenerateRing)?;
    }
    let axis = cyl.axis();
    let (o, w_dir) = (axis.origin(), axis.direction());
    let r = cyl.radius();
    // Every point of a circle's polyline stands at that circle's place ([`station`]) — the number
    // the lattice lays its rows at, so a row point on an arc lies on the arc's chord exactly.
    let mut stations: Vec<f64> = Vec::new();
    let mut place: HashMap<Handle<TessVertex>, f64> = HashMap::new();
    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
        for he in &lp.half_edges {
            if let Curve::Circle(c) = model.edge_curve(he.edge) {
                let s = station(cyl, c);
                stations.push(s);
                for &h in &t.by_edge[&he.edge] {
                    place.insert(h, s);
                }
            }
        }
    }
    let mut uv: Vec<polygon::P2> = vec![[0.0, 0.0]; handles.len()];
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for ring in &rings {
        let ring_handles: Vec<Handle<TessVertex>> = ring.iter().map(|&i| handles[i]).collect();
        let (thetas, _) = unwrapped_thetas(t, cyl, &ring_handles);
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for (&i, &theta) in ring.iter().zip(&thetas) {
            let u = place
                .get(&handles[i])
                .copied()
                .unwrap_or_else(|| (t.vertices.get(handles[i]).pos - o).dot(w_dir));
            uv[i] = [u, r * theta];
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
    let interior = interior_nodes(cfg, cyl, &uv, &rings, stations, negated);
    Ok(Chart {
        uv,
        handles,
        rings,
        map: ChartMap::Cylinder { cyl: *cyl, negated },
        interior,
    })
}

/// One position of a loop's boundary ring and where in a polyline it came from — what a cut needs
/// to put a new sample into the edge the ring walks there.
#[derive(Clone, Copy)]
struct Traced {
    handle: Handle<TessVertex>,
    edge: Handle<Edge>,
    forward: bool,
    /// The sample's index in `by_edge[edge]`.
    at: usize,
}

/// [`boundary_ring`] with each position's polyline place, and the place the ring's closing step
/// arrives through (the repeat of the first sample that `boundary_ring` drops). The step from
/// position `k` to `k + 1` lies on the edge position `k + 1` was taken from.
fn traced_ring(t: &Tessellation, lp: &Loop) -> (Vec<Traced>, Option<Traced>) {
    let mut ring: Vec<Traced> = Vec::new();
    for he in &lp.half_edges {
        let poly = &t.by_edge[&he.edge];
        let n = poly.len();
        for k in 0..n {
            let at = if he.forward { k } else { n - 1 - k };
            if ring.last().map(|r| r.handle) != Some(poly[at]) {
                ring.push(Traced {
                    handle: poly[at],
                    edge: he.edge,
                    forward: he.forward,
                    at,
                });
            }
        }
    }
    let closing = (ring.len() > 1
        && ring.first().map(|r| r.handle) == ring.last().map(|r| r.handle))
    .then(|| ring.pop())
    .flatten();
    (ring, closing)
}

/// The cylinder chart's angle at each position of `ring`, unwrapped along it — every step taken the
/// short way round — and the ring's whole turn **including its closing step**: `±2π` for a ring
/// that wraps the axis, `0` for one that does not. The one spelling of `θ` on a cylinder's chart:
/// [`cylinder_chart`] lays its rings out with it and [`cut_seamless_bands`] counts their turns with
/// it, so the two cannot disagree about which rings wrap.
fn unwrapped_thetas(
    t: &Tessellation,
    cyl: &Cylinder,
    ring: &[Handle<TessVertex>],
) -> (Vec<f64>, f64) {
    use std::f64::consts::{PI, TAU};
    let axis = cyl.axis();
    let o = axis.origin();
    let (x, y) = (cyl.ref_dir(), axis.direction().cross(cyl.ref_dir()));
    let raw = |h: Handle<TessVertex>| {
        let d = t.vertices.get(h).pos - o;
        d.dot(y).atan2(d.dot(x))
    };
    let wrap = |mut step: f64| {
        while step > PI {
            step -= TAU;
        }
        while step <= -PI {
            step += TAU;
        }
        step
    };
    let mut out = Vec::with_capacity(ring.len());
    let (mut prev, mut theta) = (0.0, 0.0);
    for (k, &h) in ring.iter().enumerate() {
        let r = raw(h);
        theta = if k == 0 { r } else { theta + wrap(r - prev) };
        prev = r;
        out.push(theta);
    }
    let total = match (out.first(), ring.first()) {
        (Some(&first), Some(&h0)) => theta + wrap(raw(h0) - prev) - first,
        _ => 0.0,
    };
    (out, total)
}

/// ★★★★★ **The band-cut pre-pass — a lateral face that wraps the axis gets one generator to be cut
/// along.**
///
/// A lateral face that wraps the axis is bounded by two rims that each turn once around it, and
/// those two rings unroll to two open lines, not a polygon. So the face is cut along one generator
/// and the cut is joined back into the rectangle the chart can sweep ([`joined_band`]) — the
/// seamless modelers' faceting (Stallings, *Reduced Topology B-Reps*, §10), with the cut placed
/// by the face rather than fixed at the surface's `θ = 0`.
///
/// **One rule places it: the cut crosses as few loops as it can.** Each loop it crosses gets a new
/// sample in a polyline it shares with the face across that loop, and is severed and spliced into
/// the outer ring; a cut that needs neither is the cheapest and the safest. So, in order:
///
/// 1. **Two loops, each one closed rim edge** (a plain cylinder or bore): each rim is cut at its
///    own vertex — the first sample of its polyline, so nothing is inserted at all.
/// 2. Otherwise a generator from [`band_generator`]: it passes each wrapping ring exactly once and
///    each hole zero or two times, crosses as few holes as any such generator, and lies farthest
///    from every existing sample (a new point a few ulps beside an old one is a coincident pair the
///    sweep cannot order, and two existing samples at "the same" `θ` on two rings are one generator
///    only up to rounding). Its points go into the crossed loops' **shared** polylines (`by_edge`),
///    as samples of their arcs, so the face across each loop sees them too and the mesh stays
///    crack-free.
///
/// A face whose boundary wraps but not as exactly two rings, or that no generator cuts so — a
/// hole the generator would cross four times severs the face into more than one outer polygon,
/// and the sweep takes one — is recorded with `None`.
fn cut_seamless_bands(
    t: &mut Tessellation,
    model: &Model,
    live: &[(Handle<Face>, &Face)],
) -> HashMap<Handle<Face>, BandCut> {
    use std::f64::consts::PI;
    let mut cuts = HashMap::new();
    for &(fh, face) in live {
        let Surface::Cylinder(cyl) = model.surface_cache(face.surface) else {
            continue;
        };
        let loops: Vec<&Loop> = std::iter::once(&face.outer).chain(&face.inner).collect();
        let traced: Vec<(Vec<Traced>, Option<Traced>)> =
            loops.iter().map(|lp| traced_ring(t, lp)).collect();
        let turns: Vec<(Vec<f64>, f64)> = traced
            .iter()
            .map(|(ring, _)| {
                let handles: Vec<Handle<TessVertex>> = ring.iter().map(|p| p.handle).collect();
                unwrapped_thetas(t, cyl, &handles)
            })
            .collect();
        let wrapping: Vec<usize> = (0..turns.len())
            .filter(|&i| turns[i].1.abs() > PI)
            .collect();
        if wrapping.is_empty() {
            continue;
        }
        let cut = match wrapping[..] {
            [a, b] => whole_rims_cut(t, model, face, cyl).or_else(|| {
                // Step `k` arrives through the next position, the closing step through the
                // closing position, or the first where the ring has none.
                let (ring, closing) = &traced[a];
                let arc = |k: usize| {
                    let arriving = ring.get(k + 1).or(closing.as_ref()).unwrap_or(&ring[0]);
                    matches!(model.edge_curve(arriving.edge), Curve::Circle(_))
                };
                band_generator(&turns, [a, b], arc)
                    .and_then(|c| band_cut(t, model, cyl, &traced, &turns, [a, b], c))
            }),
            _ => None,
        };
        cuts.insert(fh, cut);
    }
    cuts
}

/// [`cut_seamless_bands`]'s first rule: a face whose two loops are each **one closed rim edge** is
/// cut at their vertices — the first sample of each polyline ([`sample_edge`] starts a closed
/// circle there). `None` for any other face, and for two vertices at different angles: the cut is
/// the segment between them, which lies on the surface only as a ruling, and the joined ring's
/// consecutive pairs are the one place the interior budget does not look. Every producer puts a
/// whole rim's vertex at `θ = 0`; the guard keeps a producer that does not from shipping a chord
/// through the solid — such a face goes to the generator instead.
fn whole_rims_cut(t: &Tessellation, model: &Model, face: &Face, cyl: &Cylinder) -> Option<Cut> {
    let [ref inner] = face.inner[..] else {
        return None;
    };
    let vertex = |lp: &Loop| {
        let [he] = lp.half_edges[..] else {
            return None;
        };
        let [v0, v1] = model.edge(he.edge).vertices;
        (v0 == v1).then(|| t.by_edge[&he.edge][0])
    };
    let rims = [vertex(&face.outer)?, vertex(inner)?];
    // One angle up to the noise of two realizations of `θ = 0` (a few ulps); anything wider is
    // another generator.
    let (o, x, y) = (
        cyl.axis().origin(),
        cyl.ref_dir(),
        cyl.axis().direction().cross(cyl.ref_dir()),
    );
    let [p, q] = rims.map(|h| t.vertices.get(h).pos - o);
    let turn = (p.dot(x) * q.dot(y) - p.dot(y) * q.dot(x))
        .atan2(p.dot(x) * q.dot(x) + p.dot(y) * q.dot(y));
    (turn.abs() < 1e-9).then_some(Cut {
        rims,
        holes: Vec::new(),
    })
}

/// A ring's steps as `(from, to)` in its unwrapped angle, the closing step last.
fn ring_steps((th, total): &(Vec<f64>, f64)) -> Vec<(f64, f64)> {
    let n = th.len();
    (0..n)
        .map(|k| (th[k], if k + 1 < n { th[k + 1] } else { th[0] + total }))
        .collect()
}

/// How many times the open step `(p, q)` passes the generator `c`, on any turn.
fn passes((p, q): (f64, f64), c: f64) -> i64 {
    use std::f64::consts::TAU;
    let (lo, hi) = (p.min(q), p.max(q));
    let first = ((lo - c) / TAU).floor() as i64 + 1;
    let last = ((hi - c) / TAU).ceil() as i64 - 1;
    (last - first + 1).max(0)
}

/// The generator a band is cut along ([`cut_seamless_bands`]'s second rule): among the middles of
/// the first wrapping ring's arc steps, in ring order, one that each wrapping ring passes exactly
/// once and each other ring (a hole) zero or two times, crossing **the fewest holes**, and of those
/// **farthest from every sample** (the first of equal clearance). `turns` is every ring's
/// [`unwrapped_thetas`]; `arc(k)` says whether step `k` of the first wrapping ring runs along an
/// arc — a ruling's step is never a candidate, whatever its two ends' angles round to (a ruling on
/// a tilted axis can have ends whose `f64` angles differ). `None` when no candidate passes.
///
/// A hole passed twice is severed into two runs the chart splices into the outer ring
/// ([`joined_band`]); one passed four or more times would sever the face into more than one outer
/// polygon. Fewest holes first, because a hole the cut avoids keeps its polylines and its ring as
/// they are — and so a face with a hole-free generator is cut where a face without holes would be.
fn band_generator(
    turns: &[(Vec<f64>, f64)],
    wrapping: [usize; 2],
    arc: impl Fn(usize) -> bool,
) -> Option<f64> {
    use std::f64::consts::TAU;
    let clearance = |c: f64| -> f64 {
        turns
            .iter()
            .flat_map(|(th, _)| th.iter())
            .map(|&a| {
                let d = (a - c).rem_euclid(TAU);
                d.min(TAU - d)
            })
            .fold(f64::INFINITY, f64::min)
    };
    let steps: Vec<Vec<(f64, f64)>> = turns.iter().map(ring_steps).collect();
    let mut best: Option<(usize, f64, f64)> = None; // (holes crossed, clearance, c)
    for (k, &(p, q)) in steps[wrapping[0]].iter().enumerate() {
        if !arc(k) {
            continue;
        }
        let c = (p + q) / 2.0;
        let mut holes = 0usize;
        let mut admissible = true;
        for (i, s) in steps.iter().enumerate() {
            let n: i64 = s.iter().map(|&st| passes(st, c)).sum();
            match (wrapping.contains(&i), n) {
                (true, 1) | (false, 0) => {}
                (false, 2) => holes += 1,
                _ => admissible = false,
            }
        }
        if !admissible {
            continue;
        }
        let gap = clearance(c);
        if best.is_none_or(|(h, g, _)| holes < h || (holes == h && gap > g)) {
            best = Some((holes, gap, c));
        }
    }
    best.map(|(_, _, c)| c)
}

/// What the band-cut pre-pass decided for a cylindrical face whose boundary wraps the axis, or
/// `None` when no generator cuts it (the chart then refuses the face as
/// [`TessError::DegenerateRing`]).
type BandCut = Option<Cut>;

/// Where the generator cuts a band: one point on each wrapping ring (loop order), and each hole it
/// crosses as its two points, ordered along the generator **from the first rim toward the second**.
#[derive(Clone, Debug)]
struct Cut {
    rims: [Handle<TessVertex>; 2],
    holes: Vec<[Handle<TessVertex>; 2]>,
}

/// The cut of one band at the generator `c` ([`cut_seamless_bands`]'s second rule): a sample put
/// into the polyline of every loop step that passes `c`. `None` when a ring passes `c` other than
/// as [`band_generator`] admits, a crossed step is not an arc, or two points on the generator do
/// not lie strictly in order along the axis.
///
/// ★ **Along the axis, by the circle's centre.** Every crossing lies on an arc (a ruling is
/// parallel to the generator, and `c` lies off every sample), and an arc's axial place is its
/// circle's — read from the curve, where every arc of one circle agrees bitwise, rather than from
/// the computed point, which carries rounding. Two points at one place would be two loops touching
/// on the cut.
#[allow(clippy::too_many_arguments)]
fn band_cut(
    t: &mut Tessellation,
    model: &Model,
    cyl: &Cylinder,
    traced: &[(Vec<Traced>, Option<Traced>)],
    turns: &[(Vec<f64>, f64)],
    wrapping: [usize; 2],
    c: f64,
) -> BandCut {
    let axis = cyl.axis();
    let (x, y) = (cyl.ref_dir(), axis.direction().cross(cyl.ref_dir()));
    // Every crossing is found before any insertion moves a polyline index (`Traced::at`).
    let mut found: Vec<(usize, Traced, bool)> = Vec::new(); // (ring, arriving, through the closing wrap)
    for (ring, (positions, closing)) in traced.iter().enumerate() {
        let steps = ring_steps(&turns[ring]);
        let n: i64 = steps.iter().map(|&s| passes(s, c)).sum();
        let at: Vec<usize> = (0..steps.len())
            .filter(|&k| passes(steps[k], c) == 1)
            .collect();
        let admitted = if wrapping.contains(&ring) {
            n == 1
        } else {
            n == 0 || n == 2
        };
        if !admitted || at.len() as i64 != n {
            return None;
        }
        for k in at {
            found.push(match (k + 1 < positions.len(), closing) {
                (true, _) => (ring, positions[k + 1], false),
                (false, Some(cl)) => (ring, *cl, false),
                // ★ **The closing step of a ring with no closing position** — one closed edge,
                // whose polyline does not repeat its first sample: the step from its last sample
                // back to the first. A new sample there goes at the polyline's end.
                (false, None) => (ring, positions[0], true),
            });
        }
    }
    let mut points: Vec<(usize, Handle<TessVertex>, f64)> = Vec::with_capacity(found.len());
    for (ring, arriving, wraps) in found {
        let Curve::Circle(circle) = model.edge_curve(arriving.edge) else {
            return None;
        };
        let pos = circle.center() + (x * c.cos() + y * c.sin()) * circle.radius();
        let h = t.vertices.push(TessVertex {
            pos,
            origin: TessOrigin::OnEdge {
                edge: arriving.edge,
                t: circle.angle_of(pos),
            },
        });
        let poly = t.by_edge.get_mut(&arriving.edge)?;
        // The step runs `poly[at − 1] → poly[at]` walked forward, `poly[at + 1] → poly[at]`
        // walked backward; the new sample goes between the two.
        let at = if wraps {
            poly.len()
        } else if arriving.forward {
            arriving.at
        } else {
            arriving.at + 1
        };
        poly.insert(at, h);
        points.push((ring, h, station(cyl, circle)));
    }
    let place_of = |ring: usize| points.iter().find(|p| p.0 == ring).map(|p| p.2);
    let (pa, pb) = (place_of(wrapping[0])?, place_of(wrapping[1])?);
    let up = pa < pb;
    let mut holes: Vec<([Handle<TessVertex>; 2], [f64; 2])> = Vec::new();
    for ring in (0..traced.len()).filter(|r| !wrapping.contains(r)) {
        let mut on: Vec<(Handle<TessVertex>, f64)> = points
            .iter()
            .filter(|p| p.0 == ring)
            .map(|p| (p.1, p.2))
            .collect();
        if on.is_empty() {
            continue;
        }
        on.sort_by(|p, q| p.1.total_cmp(&q.1));
        if !up {
            on.reverse();
        }
        holes.push(([on[0].0, on[1].0], [on[0].1, on[1].1]));
    }
    holes.sort_by(|p, q| p.1[0].total_cmp(&q.1[0]));
    if !up {
        holes.reverse();
    }
    let along: Vec<f64> = std::iter::once(pa)
        .chain(holes.iter().flat_map(|(_, k)| *k))
        .chain(std::iter::once(pb))
        .collect();
    let ordered = along
        .windows(2)
        .all(|w| if up { w[0] < w[1] } else { w[0] > w[1] });
    if !ordered {
        return None;
    }
    let a = points.iter().find(|p| p.0 == wrapping[0])?.1;
    let b = points.iter().find(|p| p.0 == wrapping[1])?.1;
    Some(Cut {
        rims: [a, b],
        holes: holes.into_iter().map(|(h, _)| h).collect(),
    })
}

/// A chart's boundary positions and its rings over them — [`face_rings`]'s shape.
type ChartRings = (Vec<Handle<TessVertex>>, Vec<Vec<usize>>);

/// A cut band's boundary as one ring: the first wrapping ring from its cut point round to that
/// point again; then up the generator through each crossed hole — from the hole's point nearer the
/// first ring, along the hole **in its own direction**, to its other point; the second ring from
/// its cut point round to that point again; then back down the generator through the holes in
/// reverse, each from its farther point along its own direction to the nearer. Each cut point thus
/// appears at both of its `θ`s ([`cylinder_chart`]).
/// Uncrossed rings (holes) follow unchanged.
///
/// ★ **Which run of a hole the climb takes is forced, not chosen.** The material lies on the left
/// of every loop, the hole's included, so leaving the generator at the hole's nearer point the walk
/// must follow the hole's own direction; the descent takes the other run the same way. A hole wound
/// against that would make the ring cross itself, which the sweep refuses.
fn joined_band(
    handles: &[Handle<TessVertex>],
    rings: &[Vec<usize>],
    cut: &Cut,
) -> Option<ChartRings> {
    let find = |c: Handle<TessVertex>| {
        rings
            .iter()
            .enumerate()
            .find_map(|(ri, r)| r.iter().position(|&i| handles[i] == c).map(|k| (ri, k)))
    };
    let [a, b] = cut.rims;
    let ((ra, ka), (rb, kb)) = (find(a)?, find(b)?);
    if ra == rb {
        return None;
    }
    let rotated = |ri: usize, k: usize| -> Vec<Handle<TessVertex>> {
        rings[ri][k..]
            .iter()
            .chain(&rings[ri][..k])
            .map(|&i| handles[i])
            .collect()
    };
    // The run of ring `ri` from position `from` to `to`, both included, in the ring's direction.
    let run = |ri: usize, from: usize, to: usize| -> Vec<Handle<TessVertex>> {
        let n = rings[ri].len();
        let len = (to + n - from) % n + 1;
        (0..len)
            .map(|s| handles[rings[ri][(from + s) % n]])
            .collect()
    };
    let mut crossed: Vec<(usize, usize, usize)> = Vec::with_capacity(cut.holes.len());
    for &[p, q] in &cut.holes {
        let ((rp, kp), (rq, kq)) = (find(p)?, find(q)?);
        if rp != rq || rp == ra || rp == rb {
            return None;
        }
        crossed.push((rp, kp, kq));
    }
    let mut out = rotated(ra, ka);
    out.push(a);
    for &(r, kp, kq) in &crossed {
        out.extend(run(r, kp, kq));
    }
    out.extend(rotated(rb, kb));
    out.push(b);
    for &(r, kp, kq) in crossed.iter().rev() {
        out.extend(run(r, kq, kp));
    }
    let mut joined = vec![(0..out.len()).collect::<Vec<usize>>()];
    for (ri, r) in rings.iter().enumerate() {
        if ri == ra || ri == rb || crossed.iter().any(|&(c, _, _)| c == ri) {
            continue;
        }
        let start = out.len();
        out.extend(r.iter().map(|&i| handles[i]));
        joined.push((start..out.len()).collect());
    }
    Some((out, joined))
}

/// The axial place of `circle` on `cyl`: its centre, measured along the axis.
///
/// ★★ **One spelling, three readers** — the band cut orders its points along the generator by it,
/// the chart gives every point of a circle's polyline it as `u`, and the lattice lays its rows at
/// it. Read from the curve, never from a point on it: every arc of one circle agrees on it bitwise
/// (see [`interior_nodes`]), while a point's own axial coordinate carries its rounding — on a
/// tilted axis an arc's samples sit a few ulps off their circle's place, and a lattice row computed
/// one way next to a boundary computed the other was a sliver leaning off the lateral — up to 90°,
/// lying in the cap plane.
fn station(cyl: &Cylinder, circle: &Circle) -> f64 {
    (circle.center() - cyl.axis().origin()).dot(cyl.axis().direction())
}

/// ★★★★★ **The points a cylindrical face's boundary does not supply, and the mesh needs.**
///
/// A band's rims sample θ finely, and every diagonal the sweep draws between them lands on a
/// closely-sampled arc, so the interior is fine for free. It stops being free where a face's
/// boundary leaves part of its θ range without arcs: the sweep then has only chords across that
/// stretch (measured on a lateral merged across erased seams: chords spanning half a turn).
///
/// So the chart — the layer that knows the surface — supplies them, on a **lattice with no free
/// parameter in it**:
///
/// * **around** (`v`): the sampler's own step, `2πr / circle_segments`, so the lattice is as fine
///   as the arcs the boundary laid down. It does **not** land on their samples: an arc is sampled
///   from its own start ([`sample_edge`]), the lattice from `θ = 0`. `i` covers the chart's whole
///   `v` range, because a chart's θ is *unwrapped* — a hole rides whole turns from the outer ring,
///   and a single turn's worth of lattice would miss it entirely.
/// * **along** (`u`): each boundary circle's place, [`station`] — the very number the chart gave
///   that circle's points. Read from the curve, never from the sampled points: one rim's points
///   agree on that coordinate mathematically and differ by ulps in `f64`, so deduplicating *those*
///   would turn one rim into 181 stations.
///   ★ Two arcs of one circle share a centre **bitwise**, and that is a guarantee rather than luck:
///   [`Model::push_edge`] derives the curve from the **canonicalized carrier pair alone** — the
///   circle arm never reads the endpoints — so two edges cut from the same circle by the same two
///   surfaces are handed identical inputs. The station count is therefore the number of planes the
///   face actually crosses, which is why `dedup` may compare `f64` with `==` here.
///
/// ★★★ **Where the boundary already stands, no lattice point enters the mesh** — a row lies at an
/// arc's height on purpose, so part of every row is boundary, and a point placed there a few ulps
/// into the material is a sliver leaning off the lateral (up to 90°). Two cases, and only one is a
/// tolerance:
///
/// * **On an arc**: the row and the arc's chord are at one `u`, bit for bit, because both are
///   [`station`] — so the point is *on* the chord and the sweep's exact predicates drop it. Nothing
///   to tune. (Read off the points instead, the two differed by ulps on a tilted axis.)
/// * **On a ruling**: a column `i·step` meets a ruling where the user's corner happens to stand at
///   a multiple of the lattice's angle (a wall at `x = 0.5` against `r = 1` is 30°, the fifteenth
///   2° column). That is not one number computed twice — it is two numbers that geometry makes
///   equal, and no one computation gives both — so a column within [`ON_RULING`] of a rising step
///   (a ruling, or a band's cut generator; the arcs are level) is withheld.
///
/// ★★ **The lowest and highest stations are not laid at all** — they are the face's own two rims,
/// where every row point lies on a rim's chord or outside the face (a ruling's ends lie on arcs, so
/// the arcs bound the `u` range). Laid anyway, an ordinary band is offered 358 points and the sweep
/// drops every one (measured): the mesh is the same, the work is not.
fn interior_nodes(
    cfg: &TessConfig,
    cyl: &Cylinder,
    uv: &[polygon::P2],
    rings: &[Vec<usize>],
    mut stations: Vec<f64>,
    negated: bool,
) -> Vec<polygon::P2> {
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
    // The boundary's steps that are not level — after the chart's one place per circle, only rulings
    // and a band's cut generator rise.
    let rising: Vec<(polygon::P2, polygon::P2)> = rings
        .iter()
        .flat_map(|ring| (0..ring.len()).map(move |k| (ring[k], ring[(k + 1) % ring.len()])))
        .map(|(a, b)| (uv[a], uv[b]))
        .filter(|(a, b)| a[0] != b[0])
        .collect();
    let near = step * ON_RULING;
    let mut out = Vec::new();
    for &s in &stations[1..stations.len() - 1] {
        let u = if negated { -s } else { s };
        let here: Vec<&(polygon::P2, polygon::P2)> = rising
            .iter()
            .filter(|(a, b)| a[0].min(b[0]) - near <= u && u <= a[0].max(b[0]) + near)
            .collect();
        let mut i = i0;
        while i <= i1 {
            let p = [u, i * step];
            if !here.iter().any(|&&(a, b)| distance_to_step(p, a, b) < near) {
                out.push(p);
            }
            i += 1.0;
        }
    }
    out
}

/// How close to a rising boundary step (a ruling) a lattice point may stand and still be offered,
/// in units of the lattice's own step — see [`interior_nodes`].
///
/// ★ **Not a knob**: over the lib suite's meshes and the census, every point withheld stands
/// within 2.8e-13 steps of a ruling and every point kept at least 0.095 steps from one, and the
/// census meshes are the same at `1e-12`, `1e-9`, `1e-6` and `1e-3`; at `1e-1` the kept 0.095
/// goes too and a mesh fails its budget. `1e-6` sits in the middle of that gap.
const ON_RULING: f64 = 1e-6;

/// The chart distance from `p` to the boundary step `a → b`.
fn distance_to_step(p: polygon::P2, a: polygon::P2, b: polygon::P2) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 {
        (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a[0] + t * dx - p[0]).hypot(a[1] + t * dy - p[1])
}

/// Triangulate one face: lay its boundary flat in the surface's chart, then run the one sweep.
fn triangulate_face(
    t: &mut Tessellation,
    model: &Model,
    cfg: &TessConfig,
    fh: Handle<Face>,
    face: &Face,
    cut: Option<BandCut>,
) -> Result<(), TessError> {
    // ★ The cache: a chart is sampled geometry for a mesh, which is what the realization is
    // for. Nothing here classifies — the kind was decided upstream, on the truth.
    let surface = model.surface_cache(face.surface);
    let Chart {
        mut uv,
        mut handles,
        rings,
        map,
        interior,
    } = match surface {
        Surface::Plane(plane) => planar_chart(t, face, plane)?,
        Surface::Cylinder(cyl) => cylinder_chart(t, model, cfg, face, cyl, cut)?,
    };
    let refs: Vec<&[usize]> = rings.iter().map(|r| r.as_slice()).collect();
    let boundary = handles.len();
    // A vertex shared by two rings is the bridge pre-pass's doing (a curved ring's sample put
    // into the straight edge it touches); one ring passing a vertex twice is a pinched face.
    let revisits = match surface {
        Surface::Plane(_) => shared_vertices(&handles, &rings),
        // Only on a plane is a shared handle one point in the chart — a band's cut point is the
        // same handle at two `θ`s. An arm rather than a boolean, so a third surface kind is a
        // compile error here rather than a silent "no bridges".
        Surface::Cylinder(_) => Vec::new(),
    };
    let (tris, rings_used) = polygon::triangulate_uv(&mut uv, &refs, &interior, &revisits)?;
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
    within_budget(t, cfg, surface, &handles, &rings_used, &tris)?;
    for tri in tris {
        push_tri(t, fh, tri.map(|i| handles[i]));
    }
    Ok(())
}

/// The mesh vertices the face's boundary passes twice: in two different rings, the bridges the
/// sweep should lay there; twice in one ring, a pinch. Handles only — which ring's straight edge
/// was split, and whether a pinch can be ordered, is `polygon`'s to read off the geometry.
fn shared_vertices(handles: &[Handle<TessVertex>], rings: &[Vec<usize>]) -> Vec<polygon::Revisit> {
    let mut seen: HashMap<Handle<TessVertex>, (usize, usize)> = HashMap::new();
    let mut out = Vec::new();
    for (ri, ring) in rings.iter().enumerate() {
        for (k, &i) in ring.iter().enumerate() {
            match seen.get(&handles[i]) {
                Some(&(rj, kj)) if rj != ri => {
                    out.push(polygon::Revisit::Bridge(polygon::Bridge {
                        ring_x: rj,
                        x: kj,
                        ring_y: ri,
                        y: k,
                    }))
                }
                Some(&(_, kj)) => out.push(polygon::Revisit::Pinch {
                    ring: ri,
                    x: kj,
                    y: k,
                }),
                None => {
                    seen.insert(handles[i], (ri, k));
                }
            }
        }
    }
    out
}

/// ★★★★★ **The mesh's one rule, asked of the face's *interior* — where it was never asked.**
///
/// `TessConfig` declares two budgets and [`circle_segments`] enforces **both** on every boundary
/// polyline. A face's interior otherwise inherits whatever sampling its boundary supplies: fine
/// while the boundary samples the curvature (a band's rims do), and silently wrong when it does
/// not — a lateral whose boundary leaves a θ stretch unsampled meshes as flat chords across it
/// (measured on one merged across erased seams: four triangles spanning half a turn, a sixth of
/// its area lost, while `validate`, watertightness, the exact volume and the face counts stay
/// green).
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
/// pieces, which is the *larger* of the two budgets' demands; a straight edge (a ruling, a
/// plane's side) strays zero and its ends share one normal. The two halves of the sentence meet
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
#[path = "tests/lib.rs"]
mod tests;
