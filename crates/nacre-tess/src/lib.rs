//! Mesh export for the nacre kernel.
//!
//! **One path**: [`tessellate`] (every live solid) or [`tessellate_solids`] (the solids given) →
//! [`Tessellation`]. Edges are sampled once into shared polylines and every face is triangulated
//! in a chart of its own surface — a plane drops an axis, a cylinder unrolls to `(z, r·θ)` — so
//! adjacent faces meet crack-free and every mesh vertex carries its origin ([`TessOrigin`]). OBJ
//! text comes from [`Tessellation::to_obj`].

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
mod polygon;

use nacre_geom::{Curve, Cylinder, Ellipse, Surface};
use nacre_math::Point3;
use nacre_store::{Handle, Store};
use nacre_topo::{Edge, Face, Loop, Model, Solid, SurfaceCache, Vertex};
use std::collections::HashMap;
use std::fmt::Write;

/// A face or a curve this tessellator cannot mesh. Returned rather than approximated: the
/// mesh is a cache, but a *wrong* cache is worse than none — silently fanning a
/// ring that is not star-shaped puts triangles outside the solid.
///
/// There is deliberately no `Fallback` variant. A fan would be one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TessError {
    /// The rings are not a polygon with sibling holes, so no triangulation of them
    /// exists: fewer than three vertices, a zero-area ring, a vertex used by two rings
    /// **by index** or repeated within one, a spike, **two segments passing through each
    /// other**, or a ring passing **through** a point it visits twice or holding one of its two
    /// sectors there inside the other (`polygon::pinch`). On a
    /// cylinder also a face whose loops wrap the axis other than as two rims, or a
    /// band no generator cuts so that it stays one polygon — every generator crosses some hole
    /// other than twice (`cut_seamless_bands`).
    ///
    /// **The cases where two vertices meet in *coordinates* are not here** but
    /// [`Self::SelfTouchingBoundary`] — a different proposition: those rings are exactly what the
    /// b-rep asked for, and it is this decomposition that has no answer for them.
    ///
    /// **The crossing clause covers a hole crossing its outer ring**, not only a ring crossing
    /// itself: without `monotone`'s `self_touch` that input comes back `Ok` with eight confident,
    /// wrong triangles (measured).
    ///
    /// **A crossing that is only the chords' doing is not here**: where an arc's chords sag past
    /// a curve near it, the arcs there are cut finer before the sweep sees the rings
    /// (`crossing_chords`). What comes back under this name is a crossing of straight segments, or
    /// one that arcs cut until their face's set holds `FINEST_CIRCLE` columns, or down to the `f64`
    /// floor, still have.
    ///
    /// **The sweep detects these** — ear clipping would notice them only by accident, as a
    /// stall. It is checked rather than
    /// assumed because the b-rep guarantees it and this layer cannot: a decomposition
    /// handed a self-crossing ring would otherwise return a confident, wrong mesh.
    DegenerateRing,
    /// **The face's boundary touches itself, and no bridge was laid there**: some vertex lies on
    /// the boundary somewhere other than at its own two edges — on another ring, or on a
    /// non-adjacent part of its own. The region is pinched there, so it is not a disk with sibling
    /// holes and this decomposition has no triangulation of it *as the rings it was given*.
    ///
    /// **The common case is drawn, not refused.** A hole whose sample sits
    /// strictly inside a straight edge shared by two planar faces (`self_touch`'s
    /// `Interior` + `Tangent`) has that sample put into the edge's polyline by [`tessellate`]'s
    /// bridge pre-pass, the two rings are spliced into one at that point, and the sweep orders the
    /// resulting coincident pair symbolically (`polygon::sos`). **One ring passing a mesh vertex
    /// twice** — a groove whose tip touches the face's own rim, a hole touching a boss's footprint
    /// at a corner — is drawn the same way, its two visits the twins, at as many points as the face
    /// has: where each visit bounds its own sector of material, the convex sectors slide and one
    /// that is not stays (`polygon::pinch`). A face whose material is two pieces meeting only at a
    /// point — a strip whose top a pin's foot takes away up to one point of its far edge — is drawn
    /// too, each piece its own ring there. What still comes back under this name is what those
    /// roads do not cover, the bridge's each counted by the pre-pass: a vertex on another ring's
    /// *vertex* (`AtEnd`), a neighbour exactly on the touched line, a touch on a curved edge or
    /// beside a curved face, two touches on one segment or at one point, a point passed three
    /// times, and a bridge beside any other bridge or pinch.
    ///
    /// **This does not say the solid is wrong — and that is measured, not hoped.** The
    /// population is an exact **tangency**: a hole touching another ring at one point, where a
    /// sampled circle vertex lands on the touch. Both fixtures are `validate`-clean with volumes
    /// exact to `1e-9`, and
    /// **the surface is a 2-manifold at the touch**: the link of the boundary on a small sphere
    /// there is a *single* circle, because the pinched face's two lobes are joined around through
    /// the neighbouring curved face; OCCT, given the same operands, returns a body of the same
    /// volume, area and face count. What is pinched is the **face**, not the surface, so the
    /// capability that owes an answer is **not** the non-manifold test (it was right) but this
    /// layer's own: a consistent symbolic order for coincident vertices in the sweep — which
    /// `polygon::sos` supplies. For what is still refused, the
    /// loss is larger than one face: a mesh stops at the first refusal in its range ([`tessellate`]
    /// — every live solid — or [`tessellate_solids`]), and the playground's drawing stops at the
    /// first value whose mesh is refused, so one pinched face erases every body of that run.
    ///
    /// The **combinatorial** twin — two rings sharing an *index* — is refused one step earlier, by
    /// `monotone`'s `link`, and comes back as [`Self::DegenerateRing`]. That check calls the same
    /// thing a *pinch*; this is the spelling a chart actually produces, since `face_rings` gives
    /// every ring its own index range.
    SelfTouchingBoundary,
    /// A hole wound the same way as its outer ring: the b-rep does not keep
    /// material on every loop's left. A broken solid, not a repairable mesh.
    HoleWinding,
    /// The mesh would leave the declared budget ([`TessConfig`]) somewhere — it would misrepresent
    /// the shape there, so it is reported rather than drawn quietly wrong. Two places:
    ///
    /// * **A lateral** whose curves the budget asks to cut into more than 4,096 chords a turn —
    ///   the finest this tessellator samples. Each curve is cut as finely as it bends where it
    ///   lies, so this is a circle past `r ≈ 34,000` at the default `tol`, or a steep ellipse on a
    ///   large cylinder: its chords sag like a circle of its semi-major axis near its major ends
    ///   (a rod of radius 20,000 cut at slope 5; at slope 2 it meshes). On a small cylinder a
    ///   steep cut adds few chords — its fast turns are short (a radius-2 rod cut at slope `10⁴`
    ///   meshes).
    /// * **A face's triangulation** left an interior edge or a facet outside the budget
    ///   ([`within_budget`]) — its boundary curves are within it by the first point.
    OverBudget,
}

impl std::fmt::Display for TessError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::DegenerateRing => "a face's rings are not a polygon with sibling holes",
            Self::SelfTouchingBoundary => "a face's boundary touches itself",
            Self::HoleWinding => "a hole winds the same way as its outer ring",
            Self::OverBudget => "the mesh cannot keep the budget: a curve would need more than the finest sampling, or a face's triangulation left an edge or a facet outside it",
        })
    }
}

impl std::error::Error for TessError {}

// ---------------------------------------------------------------------------
// Provenance-tagged tessellation
// ---------------------------------------------------------------------------

/// Where a tessellation vertex came from — the truth it can snap back to.
/// `t`/`uv` are the exact curve/surface parameters; a curve's `t` is its angle, in `[0, 2π)`.
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
/// **A closed edge's polyline does not repeat its first point — the closure is
/// implicit.** [`sample_curve`] samples a whole rim once round from its vertex and stops short of
/// it, and `boundary_ring` drops the wrap-around duplicate for the same reason: these are mesh
/// vertices, and repeating a handle would give the face a degenerate triangle. So a consumer
/// that walks a polyline **pairwise** gets `n − 1` steps for a ring of `n`, and must add the
/// closing step itself; `edge.vertices[0] == edge.vertices[1]` is the fact to branch on — the
/// very test [`sample_curve`] uses. Written here because a consumer that does not know it draws
/// every uncut rim with a gap in it (measured in the playground's viewport).
///
/// A violation planted here can go unseen twice over: a "no degree-1 vertex" check is
/// structurally blind to a spurious closing step — an extra segment raises a degree and never
/// leaves a 1 — and a fixture with no arc edge never exercises the closure at all.
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
///
/// **A budget asks only when it is a finite positive number**; one that is not asks nothing, and
/// a curve neither budget asks of is cut into eight. A budget the mesh cannot keep is refused
/// ([`TessError::OverBudget`]), never kept loosely.
#[derive(Clone, Copy, Debug)]
pub struct TessConfig {
    /// Furthest a chord may stray from its curve, and an interior edge's midpoint from its
    /// surface, in model units.
    pub tol: f64,
    /// Maximum turn between consecutive chords, in degrees.
    ///
    /// On a regular polygon this one number is three at once: the central angle of a
    /// segment, the kink at each vertex, and the angle between adjacent facet normals
    /// — so it bounds both how polygonal a silhouette looks and how banded its shading
    /// is. `2°` puts a circle that fills an 800px viewport within 0.06px of round. Every facet's
    /// normal, at its centroid, stays within it of its surface's too. An angle that is not finite
    /// and positive asks nothing: the curves fall back on `tol`, and nothing is asked of a mesh's
    /// angles.
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
    /// **Each corner carries its face's outward normal** (`vn`, `f a//na b//nb c//nc`) —
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
/// triangulated by the same sweep. Edge polylines are sampled — finished before any face reads
/// them — and shared, so adjacent faces meet watertight. Reads the loop winding, not the `Orientation`
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
/// The same road, walked over less: edges are sampled and faces triangulated one by one, a curve
/// on the θ set of the face it bounds ([`FaceThetas`]), two solids share no face or edge, a set
/// grows only by what its own solid's planar faces ask for ([`finish_polylines`] — a round another
/// solid needs lays the same samples again), and the bridge pre-pass pairs only faces on one edge
/// — so a face's triangles are the same whichever other solids are meshed beside it.
pub fn tessellate_solids(
    model: &Model,
    solids: &[Handle<Solid>],
    cfg: &TessConfig,
) -> Result<Tessellation, TessError> {
    // `Reachable` is a `HashSet`; walk the stores in their own order and merely ask
    // membership, so the mesh stays reproducible (replay).
    let reach = model.reachable_from(solids);

    // 1. Sample every live edge into a shared polyline (the crack-free contract), finished
    // before any face reads them.
    let (live, mut sets) = theta_sets(model, cfg, &reach)?;
    let (Finished { mut t, cuts, asked }, _report) =
        finish_polylines(model, &reach, &live, &mut sets);
    // 2. Triangulate each live face, reusing the shared edge polylines.
    for &(fh, face) in &live {
        triangulate_face(
            &mut t,
            model,
            cfg,
            fh,
            face,
            cuts.get(&fh).cloned(),
            sets.faces.get(&fh),
            asked.get(&fh),
        )?;
    }

    Ok(t)
}

/// The passes that finish the shared polylines before any face reads them, **in the one order they
/// must run** — every sample they leave is one every face on the edge sees.
///
/// Each round samples every curve afresh from its θ set ([`sample_edges`]) into a new
/// [`Tessellation`], then:
///
/// 1. A cylindrical face whose boundary wraps the axis — a band, bounded by its two rims — gets
///    one generator to be cut along, its points put into the crossed loops' polylines
///    ([`cut_seamless_bands`]).
/// 2. Each planar face's boundary is asked where it meets itself ([`crossing_chords`]): where two
///    chords pass through each other, each arc there asks its lateral's set for a column at the
///    chord's middle. **After the band cut**: a sample put into an arc pushes its chord out toward
///    the curve, so it can make a crossing that was not there.
///
/// A round in which no set grew is the last: its polylines are the ones shipped, and the bridge
/// pre-pass reads the touches its scans found ([`bridge_shared_edges`] — what it inserts lies
/// exactly on a straight edge, so it changes no polyline's shape). A round that grew a set is
/// thrown away whole: the store is append-only, so sampling a curve again into the same
/// tessellation would leave its old samples in it, and an OBJ writes every vertex the mesh holds.
/// Measured, that resampling and the band cut cost about 1.4 ms a round next to 100 bored disks.
///
/// **Every planar face's chart is built again each round, not only the crossing one's**: a column
/// a crossing on one cap asks for goes on every curve of its lateral, the rim another cap is
/// bounded by among them, so that cap's boundary moves too and can cross where it did not. **A
/// chart is asked once** ([`polygon::Asked`]): the scan reads the chart's points and rings and
/// nothing else, so a face whose chart is the previous round's, bit for bit, takes that round's
/// answer. What comes back is each planar face's last scan, which its triangulation reads the
/// same way — a face no bridge or pinch changed is not scanned again there. Measured next to 100
/// bored disks, a model taking nine rounds meshed in 711 ms asking every chart afresh and in
/// 457 ms asking each once (a one-round model, 437 ms). A model with no crossing has one round —
/// sampled, cut and scanned once.
///
/// **It ends.** Sets only grow, and a round goes on only when some column was laid. A column an
/// edge asks for stands inside the edge's own range, so its polyline either gains a chord or, where
/// the new column withholds a made column beside it ([`FaceThetas::columns`]), keeps its count —
/// once per made column. A curve asks only while its face's set holds fewer than
/// [`FINEST_CIRCLE`] columns ([`arc_midpoint`], [`FaceThetas::full`]), so each set grows at most to
/// that, and the rounds are bounded.
///
/// **A band's generator can move between rounds**: it is chosen clear of every sample
/// ([`band_generator`]), and a column laid beside it moves the clearest step elsewhere. It is
/// boundary, not a column of the set, so it never asks for one; what it moves is only where the
/// next round's crossings are looked for. Measured over the suite's crossing fixtures: no band's
/// generator moved.
fn finish_polylines(
    model: &Model,
    reach: &nacre_topo::Reachable,
    live: &[(Handle<Face>, &Face)],
    sets: &mut ThetaSets,
) -> (Finished, BridgeReport) {
    // Each planar face's last scan — a map only taken from and put into, never walked.
    let mut asked: HashMap<Handle<Face>, polygon::Asked> = HashMap::new();
    loop {
        let mut t = Tessellation::default();
        sample_edges(&mut t, model, reach, sets);
        let cuts = cut_seamless_bands(&mut t, model, live);
        let mut touches: Vec<Touch> = Vec::new();
        let mut report = BridgeReport::default();
        let mut grew = false;
        for &(fh, face) in live {
            let plane = match model.surface_cache(face.surface) {
                nacre_topo::SurfaceCache::Plane { plane, .. } => plane,
                // Crossings and bridges are asked on planar charts only — on a cylinder's chart a
                // circle is level, a ruling upright and an ellipse a graph of θ, so sampling cannot
                // cross them. An arm rather than an `else`, so a third surface kind is a compile
                // error here and its author decides.
                nacre_topo::SurfaceCache::Cylinder { .. } => continue,
            };
            // A face whose chart or rings are refused is left for its own triangulation to name.
            let Some(BoundaryScan {
                chart,
                scan,
                wanted,
            }) = crossing_chords(&t, model, face, plane, sets, asked.remove(&fh))
            else {
                continue;
            };
            for (edge, theta) in wanted {
                grew |= sets.of_mut(edge).add(theta);
            }
            bridgeable_touches(fh, &chart, scan.meets(), &mut touches, &mut report);
            asked.insert(fh, scan);
        }
        if !grew {
            let report = bridge_shared_edges(&mut t, model, live, touches, report);
            return (Finished { t, cuts, asked }, report);
        }
    }
}

/// What [`finish_polylines`] hands the faces' triangulation — beside the bridge pre-pass's report,
/// which only a test reads.
struct Finished {
    /// The finished edge polylines, and no triangles yet.
    t: Tessellation,
    /// Each band's cut ([`cut_seamless_bands`]).
    cuts: HashMap<Handle<Face>, BandCut>,
    /// Each planar face's last boundary scan, which its triangulation takes where its chart is the
    /// scan's ([`polygon::Asked`]).
    asked: HashMap<Handle<Face>, polygon::Asked>,
}

/// The live faces in store order, so the later phases walk exactly that set, with each curve's θ
/// set ([`ThetaSets`]) — the sets its curves are sampled on, the crossing pass grows, and a
/// lateral's lattice reads again.
///
/// A curve edge takes the set of the face it bounds on its cylinder. In a valid model there is
/// exactly one: a curve is derived only for a plane and a cylinder (`Model::derive_edge_curve`),
/// and an edge's two faces are its two carriers (`validate`'s `EdgeCarrierMismatch`). Where a
/// broken model gives two, the first in store order answers — the other face's bands may then
/// lean, and the facet judgement refuses them rather than draw them ([`facet_in_budget`]); where
/// it gives none, the curve stands on a set of its own, its count and its two ends.
///
/// [`TessError::OverBudget`] where a set would hold more than the finest sampling ([`walk`]) —
/// before any face is drawn, so a model refused there reports that first.
fn theta_sets<'m>(
    model: &'m Model,
    cfg: &TessConfig,
    reach: &nacre_topo::Reachable,
) -> Result<(LiveFaces<'m>, ThetaSets), TessError> {
    let live: LiveFaces<'m> = (0..model.face_count() as u32)
        .filter_map(|i| model.face_handle_at(i))
        .map(|h| (h, model.face(h)))
        .filter(|(fh, _)| reach.faces.contains(fh))
        .collect();
    let mut sets = ThetaSets::default();
    for &(fh, face) in &live {
        let Some(cyl) = model.surface_cache(face.surface).cylinder() else {
            continue;
        };
        let set = FaceThetas::of_face(cfg, model, face, cyl).ok_or(TessError::OverBudget)?;
        sets.faces.insert(fh, set);
        for he in loops_of(face).flat_map(|lp| &lp.half_edges) {
            if matches!(model.edge_curve(he.edge), Curve::Ellipse(_)) {
                sets.of_edge.entry(he.edge).or_insert(fh);
            }
        }
    }
    let mut i = 0u32;
    while let Some(eh) = model.edge_handle_at(i) {
        i += 1;
        if !reach.edges.contains(&eh) || sets.of_edge.contains_key(&eh) {
            continue;
        }
        if let Curve::Ellipse(c) = model.edge_curve(eh) {
            let set =
                FaceThetas::of_curve(cfg, model, model.edge(eh), c).ok_or(TessError::OverBudget)?;
            sets.own.insert(eh, set);
        }
    }
    Ok((live, sets))
}

/// **Every curve's θ set** — the face's it bounds on its cylinder, or its own where no live face
/// does ([`theta_sets`] says which). One set per face, never per cylinder ([`FaceThetas`]).
#[derive(Default)]
struct ThetaSets {
    faces: HashMap<Handle<Face>, FaceThetas>,
    /// The face whose set each curve edge stands on.
    of_edge: HashMap<Handle<Edge>, Handle<Face>>,
    /// A curve no live face of its cylinder bounds (a broken model's): its own set.
    own: HashMap<Handle<Edge>, FaceThetas>,
}

impl ThetaSets {
    /// The set curve edge `e` is sampled on.
    fn of(&self, e: Handle<Edge>) -> &FaceThetas {
        match self.of_edge.get(&e) {
            Some(f) => &self.faces[f],
            None => &self.own[&e],
        }
    }

    /// The set curve edge `e` asks a column of.
    fn of_mut(&mut self, e: Handle<Edge>) -> &mut FaceThetas {
        match self.of_edge.get(&e) {
            Some(f) => self
                .faces
                .get_mut(f)
                .expect("every face in `of_edge` has a set"),
            None => self.own.get_mut(&e).expect("every lone curve has a set"),
        }
    }
}

/// Sample every live edge into `by_edge` — a straight edge as its two ends, a curve on its θ set
/// ([`sample_curve`]) — sharing one mesh vertex per model vertex, so the faces on an edge meet it
/// crack-free.
fn sample_edges(
    t: &mut Tessellation,
    model: &Model,
    reach: &nacre_topo::Reachable,
    sets: &ThetaSets,
) {
    let mut vmap: HashMap<Handle<Vertex>, Handle<TessVertex>> = HashMap::new();
    let mut i = 0u32;
    while let Some(eh) = model.edge_handle_at(i) {
        i += 1;
        let edge = model.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        let polyline = match model.edge_curve(eh) {
            Curve::Line(_) => vec![
                vertex_of(t, &mut vmap, model, edge.vertices[0]),
                vertex_of(t, &mut vmap, model, edge.vertices[1]),
            ],
            Curve::Ellipse(c) => sample_curve(t, &mut vmap, model, eh, edge, c, sets.of(eh)),
        };
        t.by_edge.insert(eh, polyline);
    }
}

/// The first round's sampling — the sets as the model gives them, every live edge sampled on
/// them — for a test that reads the polylines before anything finishes them.
#[cfg(test)]
fn sample_live_edges<'m>(
    t: &mut Tessellation,
    model: &'m Model,
    cfg: &TessConfig,
    reach: &nacre_topo::Reachable,
) -> (LiveFaces<'m>, ThetaSets) {
    let (live, sets) = theta_sets(model, cfg, reach).expect("the budget can be kept");
    sample_edges(t, model, reach, &sets);
    (live, sets)
}

/// The live faces in store order, each with its handle.
type LiveFaces<'m> = Vec<(Handle<Face>, &'m Face)>;

/// A face's loops, outer first.
fn loops_of(face: &Face) -> impl Iterator<Item = &Loop> {
    std::iter::once(&face.outer).chain(face.inner.iter())
}

/// **The angle of `p` about `cyl`'s axis, from its `ref_dir`** — the angle the chart unrolls by
/// ([`unwrapped_thetas`]) and a face's θ set reads its vertices and its curves' spans by
/// ([`FaceThetas`], [`sample_curve`]). In `[−π, π]`.
///
/// It is also a curve's own parameter, up to rounding: a rim is derived from its cylinder's cached
/// `axis` and `ref_dir` as they are (`Ellipse::from_unit_frame` in `Model::derive_edge_curve`),
/// so its θ ([`Ellipse::angle_of`]) is this angle read about the rim's centre rather than the
/// axis origin — the same number, a few ulps apart. A sample is placed by the face's angle and
/// evaluated by the curve ([`Ellipse::point_at`]); nothing compares the two spellings.
fn angle_about(cyl: &Cylinder, p: Point3) -> f64 {
    let axis = cyl.axis();
    let d = p - axis.origin();
    d.dot(axis.direction().cross(cyl.ref_dir()))
        .atan2(d.dot(cyl.ref_dir()))
}

/// `x` brought into `[0, τ)`. `rem_euclid` alone can return `τ` itself — a value a hair below `0`
/// rounds up to it — and that is the column `0`.
fn in_turn(x: f64) -> f64 {
    let a = x.rem_euclid(std::f64::consts::TAU);
    if a == std::f64::consts::TAU { 0.0 } else { a }
}

/// The number halfway round, CCW, from column number `early` to `late` — itself a column's number,
/// in `[0, τ)`, so the check that judges the point there and the set that lays the column read one
/// number.
fn ccw_midpoint(early: f64, late: f64) -> f64 {
    in_turn(early + (late - early).rem_euclid(std::f64::consts::TAU) / 2.0)
}

/// **A cylindrical face's θ set: every curve of the face is sampled on it, and the face's
/// interior lattice lays its columns on it** — so the samples of any two curves stand on common
/// generators, and every facet across a band between them spans one.
///
/// Without it a band between two curves sampled on their own steps tilts by the offset samples'
/// sag over its height: an arc cut evenly between its ends against the full circle below it, or an
/// ellipse cut by its turn against the circle of the same cylinder. Thinner than about
/// `r(1 − cos(Δθ/2))/tan α` (`α` the angular budget) the face was past the budget and refused —
/// measured at `h = 0.1` on a disk of radius 50 with its rim cut down by a box, at a band `0.004`
/// high on a rod of radius 2 cut obliquely at `1/20`.
///
/// The set is the **made columns**, the **θ of every vertex of the face**, and **the columns
/// crossing chords ask for**:
///
/// * The made columns are where every curve of the face keeps the budget on each chord where the
///   chord lies ([`walk`]): the even lattice `τk/N` of the face's floor `N` (`count` — the
///   cylinder's own circle's) where nothing bends faster, finer only where an ellipse turns or
///   sags faster. The vertices' columns only shorten chords, so each curve keeps its budget.
/// * A vertex's θ is read once, about the axis ([`angle_about`]). Vertices at one θ are one
///   column, by [`ON_RULING`]: a ruling's two ends, which their own coordinates put a few ulps
///   apart, and vertices geometry makes stand on one generator — no stored number is their θ, so
///   there is no one place to read it from. A lattice column beside a vertex's is not laid either
///   ([`FaceThetas::columns`]).
/// * A chord of a planar face's boundary that crosses another's asks for a column at its arc's
///   middle ([`crossing_chords`], [`FaceThetas::add`]), and every curve of the face spanning it
///   and the lattice take it — so the partner of a cut arc is cut there too, and the band between
///   them still spans generators. Laid on the one arc alone, it leaves the partner a sample short:
///   with a column fused onto a disk past a hole that comes `5e-4` from the rim, its top `0.01`
///   above the disk, a band that thin leans past the budget and the mesh is refused (measured).
///
/// **A face, not a cylinder.** A boolean reuses its operands' cylinder handles, so a set per handle
/// would read the vertices of every solid on it and a face's samples would change with what else
/// is meshed — against [`tessellate_solids`]' promise.
///
/// **One sample a curve ends up with is not on it**: the band cut puts its generator's into every
/// ring it passes after sampling ([`finish_polylines`]), so the curves still share it — and the
/// lattice is withheld there, the generator being boundary.
struct FaceThetas {
    cyl: Cylinder,
    /// The face's floor — the count a turn the surface itself asks everywhere ([`circle_segments`]
    /// of the cylinder's radius); the made columns' unit is a power-of-two share of its step.
    count: usize,
    /// The made columns' unit, `count·2ˢ` a turn.
    grid: u64,
    /// The made columns ([`walk`]), ascending, in units of `grid`, in `[0, grid)`; the first is `0`.
    made: Vec<u64>,
    /// Ascending in `[0, τ)`, each with the vertices standing there — none for a column a crossing
    /// asked for.
    columns: Vec<(f64, Vec<Handle<Vertex>>)>,
}

impl FaceThetas {
    /// The set of `face`, a face on `cyl` — `None` where the budget asks more of its curves, or of
    /// the surface, than [`FINEST_CIRCLE`] columns ([`walk`]).
    fn of_face(cfg: &TessConfig, model: &Model, face: &Face, cyl: &Cylinder) -> Option<FaceThetas> {
        let half_edges = || loops_of(face).flat_map(|lp| &lp.half_edges);
        let curves: Vec<Ellipse> = half_edges()
            .filter_map(|he| match model.edge_curve(he.edge) {
                Curve::Ellipse(c) => Some(*c),
                Curve::Line(_) => None,
            })
            .collect();
        let vertices = half_edges().flat_map(|he| model.edge(he.edge).vertices);
        FaceThetas::new(cfg, model, *cyl, &curves, vertices)
    }

    /// The set a curve no live face of its cylinder bounds stands on alone — its own count and its
    /// two ends (a broken model's; see [`theta_sets`]) — `None` as in [`Self::of_face`].
    fn of_curve(cfg: &TessConfig, model: &Model, edge: &Edge, c: &Ellipse) -> Option<FaceThetas> {
        // The curve's own frame is its cylinder's (see [`angle_about`]); its radius is positive.
        let cyl = Cylinder::from_unit_frame(c.center(), c.axis(), c.ref_dir(), c.radius())
            .expect("a curve's radius is positive");
        FaceThetas::new(cfg, model, cyl, &[*c], edge.vertices.into_iter())
    }

    /// The set on `cyl` of `curves` and `vertices`: the made columns [`walk`] lays, and the
    /// vertices' columns — `None` where the budget asks more than [`FINEST_CIRCLE`] columns.
    fn new(
        cfg: &TessConfig,
        model: &Model,
        cyl: Cylinder,
        curves: &[Ellipse],
        vertices: impl Iterator<Item = Handle<Vertex>>,
    ) -> Option<FaceThetas> {
        use std::f64::consts::TAU;
        // The floor is the surface's own ask, `(r, 1)`: a chord's ends' normals part by its width,
        // and it stands off the cylinder as a chord of its circle does. What the curves ask beyond
        // it, where they bend, the walk reads (`Ellipse::pace_over`) — a curve's own circle is the
        // cylinder's, and a minor end turns slower than it.
        let count = circle_segments(cfg, cyl.radius())?;
        let grid = grid_of(count);
        let made = walk(cfg, curves, count, grid)?;
        let mut vs: Vec<Handle<Vertex>> = vertices.collect();
        vs.sort_by_key(|v| v.index());
        vs.dedup();
        let mut at: Vec<(f64, Vec<Handle<Vertex>>)> = vs
            .into_iter()
            .map(|v| (in_turn(angle_about(&cyl, model.vertex_point(v))), vec![v]))
            .collect();
        at.sort_by(|a, b| a.0.total_cmp(&b.0));
        // Vertices at one θ are one column (see the type's doc).
        let near = column_tolerance(count);
        let mut columns: Vec<(f64, Vec<Handle<Vertex>>)> = Vec::new();
        for (theta, vs) in at {
            match columns.last_mut() {
                Some(last) if theta - last.0 <= near => last.1.extend(vs),
                _ => columns.push((theta, vs)),
            }
        }
        if columns.len() > 1 && columns[0].0 + TAU - columns[columns.len() - 1].0 <= near {
            let (_, vs) = columns.pop().expect("more than one");
            columns[0].1.extend(vs);
        }
        Some(FaceThetas {
            cyl,
            count,
            grid,
            made,
            columns,
        })
    }

    /// **Lay a column at `theta`**, a column's number in `[0, τ)`, standing with no vertex —
    /// `false` where the set already holds that very number, and nothing changes.
    ///
    /// A crossing chord's middle ([`arc_midpoint`]) is placed apart on purpose, so it is not one of
    /// the two numbers [`ON_RULING`] merges (geometry made them equal, or rounding split them): a
    /// column is refused only where it would add nothing. Like a vertex's column, it withholds a
    /// made column within [`FaceThetas::made_near`] of it ([`FaceThetas::columns`]).
    fn add(&mut self, theta: f64) -> bool {
        debug_assert!(
            (0.0..std::f64::consts::TAU).contains(&theta),
            "a column's number is in [0, τ): {theta}"
        );
        let k = self.columns.partition_point(|c| c.0 < theta);
        if self.columns.get(k).is_some_and(|c| c.0 == theta) {
            return false;
        }
        self.columns.insert(k, (theta, Vec::new()));
        true
    }

    /// [`ON_RULING`] of the floor's step ([`FaceThetas::count`]) — the widest a made column's
    /// [`FaceThetas::made_near`] is, so the window held columns are read over; and how close two
    /// vertices stand to be one column.
    fn near(&self) -> f64 {
        column_tolerance(self.count)
    }

    /// **How close to a held column (a vertex's or a crossing's) made column `i` may stand and
    /// still be laid**, in radians — [`ON_RULING`] of its narrower step. Its own step, not the
    /// floor's: a column withheld moves its chord's end by up to this, so it stays the same share
    /// of the chord wherever the walk cut finer.
    fn made_near(&self, i: usize) -> f64 {
        let n = self.made.len();
        let (prev, next) = (
            if i == 0 {
                self.made[n - 1]
            } else {
                self.made[i - 1]
            },
            if i + 1 == n {
                self.grid
            } else {
                self.made[i + 1]
            },
        );
        let before = (self.made[i] + self.grid - prev) % self.grid;
        let before = if before == 0 { self.grid } else { before };
        let step = before.min(next - self.made[i]);
        std::f64::consts::TAU * step as f64 / self.grid as f64 * ON_RULING
    }

    /// **Whether the set holds as many columns as a face's set may** — [`FINEST_CIRCLE`], where
    /// [`walk`] refuses a face and [`arc_midpoint`] stops asking.
    fn full(&self) -> bool {
        self.made.len() + self.columns.len() >= FINEST_CIRCLE
    }

    /// `p`'s θ in the set's frame, in `[0, τ)`.
    ///
    /// The curves of the face share the frame: a rim is derived from its cylinder's cached `axis`
    /// and `ref_dir` bit for bit, so a sample placed by this angle is evaluated on the curve's own
    /// parameter ([`angle_about`]) — [`sample_curve`] asserts it.
    fn angle(&self, p: Point3) -> f64 {
        in_turn(angle_about(&self.cyl, p))
    }

    /// **Where curve `edge` of the face starts and how far it turns, CCW**, by the set's angle of
    /// its two vertices — the one range [`sample_curve`] lays its columns over and
    /// [`arc_midpoint`] bounds its cutting by. A closed edge turns once from its vertex; an arc's
    /// stored `[v0, v1]` order is CCW about the axis (the convention `Model::derive_edge_curve`'s
    /// circle arm states).
    fn span(&self, model: &Model, edge: &Edge) -> (f64, f64) {
        let [v0, v1] = edge.vertices;
        let t0 = self.angle(model.vertex_point(v0));
        let dt = if v0 == v1 {
            std::f64::consts::TAU
        } else {
            in_turn(self.angle(model.vertex_point(v1)) - t0)
        };
        (t0, dt)
    }

    /// **The number of the column vertex `v` stands at** — the stored θ, which is another
    /// vertex's where [`ON_RULING`] made them one column. `None` for a vertex the set does not
    /// hold; a curve's own two vertices it always holds (`of_face` reads every vertex of the face,
    /// `of_curve` the curve's two).
    fn column_of(&self, v: Handle<Vertex>) -> Option<f64> {
        self.columns.iter().find(|c| c.1.contains(&v)).map(|c| c.0)
    }

    /// **The set's columns in `[lo, hi]`**, ascending, each as `(place, number, vertices)` — the
    /// vertices standing there (none for a made column or a crossing's). `lo` and `hi` are
    /// unwrapped, so a range past `τ` or below `0` gets the columns whole turns away: `place` is
    /// where the column stands in that range, `number` the column's own, in `[0, τ)` — a held
    /// column's stored θ, a made column's `τu/G` with `u < G`.
    ///
    /// **A column has one number.** `place` only picks and orders the columns along a range (and
    /// sets a chart's coordinate); everything that names the column — a sample's position and its
    /// `t`, a crossing chord's middle — uses `number`. A held column's `place` is `number` plus
    /// whole turns in `f64`, which is not `number` again: two curves holding one column on
    /// different turns (a curve that starts a hair below `0` takes its range a turn up) would name
    /// it ulps apart. A made column's `place` is `τ(u + kG)/G`, exact in its integer.
    ///
    /// Every column the set holds (a vertex's or a crossing's) in the range is laid, and every
    /// made column that stands farther than [`FaceThetas::made_near`] from every held column —
    /// including one just outside the range, so a made column a hair inside a curve's end is not
    /// laid beside it. A made column is picked by its unit against `lo·G/τ` and `hi·G/τ`; `G` is
    /// the floor's count times a power of two, so where every made column is the floor's (a face
    /// with no ellipse) the pick is the floor lattice's `⌈lo·N/τ⌉ ..= ⌊hi·N/τ⌋`, bit for bit.
    fn columns(&self, lo: f64, hi: f64) -> Vec<(f64, f64, &[Handle<Vertex>])> {
        use std::f64::consts::TAU;
        let near = self.near();
        // The held columns over `[lo − near, hi + near]`: turn by turn, each turn's ascending.
        let mut anchors: Vec<(f64, f64, &[Handle<Vertex>])> = Vec::new();
        let mut turn = ((lo - near) / TAU).floor();
        while turn * TAU <= hi + near {
            for (theta, vs) in &self.columns {
                let c = theta + turn * TAU;
                if lo - near <= c && c <= hi + near {
                    anchors.push((c, *theta, vs));
                }
            }
            turn += 1.0;
        }
        let g = self.grid as f64;
        let (first, last) = ((lo * g / TAU).ceil(), (hi * g / TAU).floor());
        debug_assert!(
            first.abs() < 16.0 * g && last.abs() < 16.0 * g,
            "a range within the turns a made column's place is exact over: [{lo}, {hi}]"
        );
        let mut out: Vec<(f64, f64, &[Handle<Vertex>])> = anchors
            .iter()
            .filter(|a| (lo..=hi).contains(&a.0))
            .copied()
            .collect();
        let mut k = (first / g).floor();
        while k * g <= last {
            for (i, &u) in self.made.iter().enumerate() {
                let unit = u as f64 + k * g;
                if unit < first || last < unit {
                    continue;
                }
                let place = TAU * unit / g;
                let near = self.made_near(i);
                let at = anchors.partition_point(|a| a.0 < place);
                let clear = [at.checked_sub(1), Some(at)]
                    .into_iter()
                    .flatten()
                    .filter_map(|j| anchors.get(j))
                    .all(|a| (a.0 - place).abs() > near);
                if clear {
                    out.push((place, TAU * u as f64 / g, &[]));
                }
            }
            k += 1.0;
        }
        out.sort_by(|a, b| a.0.total_cmp(&b.0));
        out
    }
}

/// **Where two chords of a planar face pass through each other, the θs the arcs there ask their
/// laterals' sets for** — with the face's chart and how its boundary meets itself, as it stands.
///
/// A face's boundary curves do not meet — the b-rep says so — but each is drawn as chords, and an
/// arc's chords stand inside it by up to its sagitta. Two curves closer than that can have chords
/// that pass through each other: a hole whose rim comes within a sagitta of the outer rim, a groove
/// whose tip does, or a groove side leaving a rim point closer to the tangent than half a segment's
/// turn. The face is valid, and the sweep refuses its rings ([`TessError::DegenerateRing`]).
///
/// So the boundary is asked, with the sweep's own exact scan ([`polygon::meets`]), where it passes
/// through itself, and **every** arc chord in a crossing asks for a column at its arc's middle
/// ([`arc_midpoint`]). The column goes into the set of the lateral the arc bounds, so the next
/// round samples it on every curve of that lateral — the partner rim and its lattice too, not the
/// one arc alone ([`finish_polylines`] runs the rounds). A chord cut in two moves toward its curve,
/// and the curves are apart, so the crossings go; a crossing between two straight segments has
/// nothing to ask and is the sweep's to refuse. Nothing here changes a polyline.
///
/// Measured: at most 8 rounds that lay columns where a hole or a groove tip comes within `1e-8`
/// of a rim, or a groove side leaves a rim point `0.01°` off the tangent; a hole near-concentric
/// with its rim all the way round wants whole circles cut finely — twice as fine at a gap of
/// `1e-4`, about 2,300 samples a circle at `1e-6` — and at `1e-8` the cap is reached and the face
/// refused by name.
///
/// **One scan, three readers.** The last round's answer is the boundary as it is shipped, which is
/// what the bridge pre-pass reads its touches off ([`bridgeable_touches`]) — asked once, not again
/// by a second pass over the same polylines. Measured as a separate pass, the second ask was 7.7%
/// of the census's meshing time. The face's triangulation takes it too, where its chart is the
/// scan's ([`polygon::Asked`]); and `prior`, the face's scan from the round before, is this round's
/// answer where the chart has not moved.
///
/// **Faces that cannot cross are scanned and ask nothing**: one with no circle edge (straight
/// segments do not sag), and one bounded by one whole circle alone — its samples lie on the circle
/// in angular order, a convex polygon.
///
/// **A model with no crossing gets no column** — and a boundary with no crossing gets new samples
/// only where it shares a lateral with one that crossed. Measured over the workspace suite, the
/// census, the ignored sweep, nacre-kit and the playground's suite: no face crosses but those of
/// this crate's own crossing tests.
///
/// A crossing that lands on a vertex (`Crossing::Through`) asks for the chord it passes through and
/// both chords at the vertex. Every crossing measured was `Segments` — that arm is argued, its
/// population zero.
fn crossing_chords(
    t: &Tessellation,
    model: &Model,
    face: &Face,
    plane: &nacre_geom::Plane,
    sets: &ThetaSets,
    prior: Option<polygon::Asked>,
) -> Option<BoundaryScan> {
    use polygon::Crossing;
    let loops: Vec<&Loop> = loops_of(face).collect();
    let curved = loops
        .iter()
        .flat_map(|lp| &lp.half_edges)
        .any(|he| matches!(model.edge_curve(he.edge), Curve::Ellipse(_)));
    let one_circle = loops.len() == 1 && face.outer.half_edges.len() == 1;
    let chart = planar_chart(t, face, plane).ok()?;
    let scan = {
        let refs: Vec<&[usize]> = chart.rings.iter().map(|r| r.as_slice()).collect();
        polygon::meets(&chart.uv, &refs, prior).ok()?
    };
    let meets = scan.meets();
    if !curved || one_circle || meets.crossings.is_empty() {
        return Some(BoundaryScan {
            chart,
            scan,
            wanted: Vec::new(),
        });
    }
    // A chart position is `(ring, k)` of the same walk (`boundary_ring` is `traced_ring`'s
    // handles), the rings numbered in loop order over consecutive positions.
    let traced: Vec<(Vec<Traced>, Option<Traced>)> =
        loops.iter().map(|lp| traced_ring(t, lp)).collect();
    let place = |i: usize| -> (usize, usize) {
        let ring = chart
            .rings
            .iter()
            .position(|r| r[0] <= i && i < r[0] + r.len())
            .expect("a chart position lies in one ring");
        (ring, i - chart.rings[ring][0])
    };
    let before = |i: usize| {
        let (ring, k) = place(i);
        let r = &chart.rings[ring];
        r[(k + r.len() - 1) % r.len()]
    };
    // Each chord in a crossing, named by its first position.
    let mut chords: Vec<usize> = Vec::new();
    for c in &meets.crossings {
        match *c {
            Crossing::Segments(a, b) => chords.extend([a, b]),
            Crossing::Through { vertex, segment } => {
                chords.extend([segment, before(vertex), vertex])
            }
        }
    }
    chords.sort_unstable();
    chords.dedup();
    let mut wanted = Vec::new();
    for i in chords {
        let (ring, k) = place(i);
        let (positions, closing) = &traced[ring];
        let step = arriving(positions, *closing, k).0;
        if let Some(w) = arc_midpoint(t, model, sets, positions[k], step) {
            wanted.push(w);
        }
    }
    Some(BoundaryScan {
        chart,
        scan,
        wanted,
    })
}

/// What [`crossing_chords`] found on one planar face's boundary.
struct BoundaryScan {
    chart: Chart,
    /// How the boundary meets itself, kept with the chart it asked ([`polygon::Asked`]).
    scan: polygon::Asked,
    /// Each crossing arc's edge and the θ it asks its lateral's set for.
    wanted: Vec<(Handle<Edge>, f64)>,
}

/// The θ that cuts the chord from `from` to the step's arrival in two, on the edge that step runs
/// along — or `None` when the chord is not an arc's, its edge is already cut as finely as
/// [`FINEST_CIRCLE`] lets it be, or the midpoint is one of its ends in `f64`.
///
/// The polyline runs CCW about the circle ([`FaceThetas::span`]), so the chord's earlier end is
/// `from` walked forward and the arrival walked backward. **The θ stands halfway round between the
/// numbers of the two columns the chord joins** — a sample's `t`, a vertex's column
/// ([`FaceThetas::column_of`]) — so two curves of one lateral crossing between the same two
/// columns (a bore's two rims, each crossing on its cap) ask for one number, and the set lays one
/// column. That number is in `[0, τ)` ([`ccw_midpoint`]), so the `f64`-floor check below judges
/// the point the sample will stand at. Read off the ends' positions instead, about each rim's own
/// centre (`Ellipse::angle_of`), the two asks were ulps apart on a tilted axis: two columns
/// `1.1e-16` to `2.2e-16` apart left each rim a chord of next to no length, and 33 of 270 disks
/// bored or fused with a column near the rim, on four axes, were refused by the mesh, all on the
/// two 3-4-5 axes with a zero component
/// (measured); read off the samples' unwrapped places, the asks parted wherever one rim's range
/// started a hair below `0`, which far from the origin refused 57 of 270 (at `(1000, −2000, 3000)`).
///
/// The sets are read while the round's scans grow them ([`finish_polylines`]): a vertex's column
/// is never moved by [`FaceThetas::add`], which lays only columns with no vertex, and the samples'
/// `t` are the round's.
///
/// **At the `f64` floor the midpoint need not lie between its ends** — for chords under about
/// `1e-12` rad about 1% of the midpoints land outside them, whatever the formula (the angle
/// difference here, the turn read off both ends, or the chord's midpoint pushed out to the circle;
/// measured over random chords on three circles). Cutting reaches that only where two curves meet
/// within `f64`, and what such a sample leaves is a crossing or a spike, which the sweep's exact
/// checks refuse by name — not a wrong mesh.
fn arc_midpoint(
    t: &Tessellation,
    model: &Model,
    sets: &ThetaSets,
    from: Traced,
    step: Traced,
) -> Option<(Handle<Edge>, f64)> {
    let Curve::Ellipse(c) = model.edge_curve(step.edge) else {
        return None;
    };
    let set = sets.of(step.edge);
    if set.full() {
        return None;
    }
    // A chord's two ends are samples of this edge or its own vertices (`traced_ring` takes a
    // shared vertex once), so each has its number.
    let number = |h: Handle<TessVertex>| match t.vertices.get(h).origin {
        TessOrigin::OnEdge { t, .. } => t,
        TessOrigin::OnVertex(v) => set
            .column_of(v)
            .expect("a curve's set holds its two vertices"),
        TessOrigin::OnFace { .. } => unreachable!("a polyline holds no face sample"),
    };
    let (early, late) = if step.forward {
        (from.handle, step.handle)
    } else {
        (step.handle, from.handle)
    };
    let theta = ccw_midpoint(number(early), number(late));
    let (a, b) = (
        t.vertices.get(from.handle).pos,
        t.vertices.get(step.handle).pos,
    );
    let pos = c.point_at(theta);
    if pos == a || pos == b {
        return None;
    }
    Some((step.edge, theta))
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

/// A touch the bridge pre-pass may bridge: `vertex` of `face`'s boundary lies on the segment `seg`
/// of another of its rings.
struct Touch {
    face: Handle<Face>,
    vertex: Handle<TessVertex>,
    seg: [Handle<TessVertex>; 2],
}

/// The touches of one planar face's boundary that [`bridge_shared_edges`] may bridge, read off the
/// scan the crossing pass made of it ([`crossing_chords`]) — a vertex touching a straight run at an
/// interior point, tangentially; one touching two segments is declined (`MultiSegment`).
fn bridgeable_touches(
    face: Handle<Face>,
    chart: &Chart,
    meets: &polygon::Meets,
    out: &mut Vec<Touch>,
    report: &mut BridgeReport,
) {
    use polygon::{TouchKind, Witness};
    let bridgeable =
        |tc: &polygon::Touch| tc.kind == TouchKind::Interior && tc.witness == Witness::Tangent;
    let mut per_vertex: HashMap<usize, usize> = HashMap::new();
    for tc in meets.touches.iter().filter(|tc| bridgeable(tc)) {
        *per_vertex.entry(tc.vertex).or_default() += 1;
    }
    for tc in meets.touches.iter().filter(|tc| bridgeable(tc)) {
        if per_vertex[&tc.vertex] > 1 {
            report.declined.push((face, Declined::MultiSegment));
            continue;
        }
        out.push(Touch {
            face,
            vertex: chart.handles[tc.vertex],
            seg: [chart.handles[tc.segment], chart.handles[tc.segment_end]],
        });
    }
}

/// **The bridge pre-pass — a tangency's sample becomes a sample of the edge it touches.**
///
/// A planar face whose hole touches another ring at exactly one point is a valid solid the
/// sweep cannot decompose: the touching vertex is a sample of the *curved* edge, and the
/// straight edge it lands on carries only its two ends ([`sample_edges`]' `Line` arm), so the
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
/// **Edge → faces is built here, from the live set** — the model holds no reverse index, and this
/// pass runs on the model as handed to it.
///
/// What it declines, it counts (see [`Declined`]); the population it was measured on is four
/// planar caps, every touch `Interior` and `Tangent` on a straight edge shared by two planes.
fn bridge_shared_edges(
    t: &mut Tessellation,
    model: &Model,
    live: &[(Handle<Face>, &Face)],
    cands: Vec<Touch>,
    mut report: BridgeReport,
) -> BridgeReport {
    struct Pending {
        edge: Handle<Edge>,
        at: usize,
        vertex: Handle<TessVertex>,
        face: Handle<Face>,
    }
    // 1. Every bridgeable touch was read before anything is changed ([`bridgeable_touches`]).
    if cands.is_empty() {
        return report;
    }

    // 2. Which faces share each edge, and which of them are planar.
    let planar: HashMap<Handle<Face>, bool> = live
        .iter()
        .map(|&(fh, f)| {
            (
                fh,
                // The **truth**, path-qualified: in this file the bare `Surface` is geom's
                // (the cache), and aliasing the topo one in would give the truth a second
                // vocabulary here. A face's kind is a fact about what it is, so the truth
                // answers it — and a new surface kind lands there first.
                matches!(model.surface(f.surface), nacre_topo::Surface::Plane { .. }),
            )
        })
        .collect();
    let mut faces_on: HashMap<Handle<Edge>, Vec<Handle<Face>>> = HashMap::new();
    for &(fh, face) in live {
        for lp in loops_of(face) {
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
        'search: for lp in loops_of(face) {
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

/// What the bridge pre-pass did, for a test that wants to see it before any face is triangulated
/// (the returned `Tessellation` holds the finished edge polylines and no triangles) — the same
/// passes [`tessellate`] runs, in its order, and refused where it refuses a curve
/// ([`TessError::OverBudget`]). Test-only (`test-util`).
#[cfg(feature = "test-util")]
pub fn bridge_report(
    model: &Model,
    cfg: &TessConfig,
) -> Result<(BridgeReport, Tessellation), TessError> {
    let reach = model.reachable();
    let (live, mut sets) = theta_sets(model, cfg, &reach)?;
    let (Finished { t, .. }, report) = finish_polylines(model, &reach, &live, &mut sets);
    Ok((report, t))
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
/// **Why the budget bounds it.** The mesher asks every facet the same of its own accord
/// (`within_budget`, `facet_in_budget`), and refuses the face otherwise; this lock measures the
/// shipped mesh again without that judgement, so it bites if the judgement is bypassed or wrong.
/// Its spelling is not the mesher's: the angle to the surface's normal at the centroid, corners in
/// the triangle's order, a bare `>` with no slack (stricter at the budget's edge), no allowance for
/// a normal its points do not determine (looser for a sliver), laterals only. A sliver the chart
/// made out of a lattice point standing a few ulps off the boundary leaned 2.6° to 90° (at 90° it
/// lies in a cap plane).
#[cfg(any(test, feature = "test-util"))]
pub fn leaning_laterals(model: &Model, t: &Tessellation, limit_deg: f64) -> Vec<Leaning> {
    let mut out = Vec::new();
    for (triangle, tri) in t.triangles.iter() {
        let surface = model.surface_cache(model.face(tri.face).surface).to_geom();
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

/// A θ one curve of a cylindrical face is sampled at that another curve of the face, spanning it,
/// is not: the face, the curve that misses it, and the θ (radians about the cylinder's axis, from
/// its `ref_dir`). Test-only (`test-util`).
#[cfg(any(test, feature = "test-util"))]
#[derive(Clone, Copy, Debug)]
pub struct Unshared {
    pub face: Handle<Face>,
    pub edge: Handle<Edge>,
    pub theta: f64,
}

/// Every θ a curve of a cylindrical face of `t` stands a sample at — its ends included — that
/// another curve of the face spans and has no sample as near as one column. Test-only
/// (`test-util`): the mesh locks call it.
///
/// **What it holds the sampler to**: the curves of one lateral stand on common generators
/// ([`FaceThetas`]), so a band between two of them never tilts by their samples' offset. Its
/// spelling is not the sampler's: it reads the shipped polylines back — angle by angle, about the
/// axis — and knows nothing of the lattice, the vertices or the count. Two angles are one when
/// they are as close as the sampler ever takes two for one column ([`ON_RULING`] of the coarsest
/// lattice, [`COARSEST_CIRCLE`] a turn — `7.9e-7` rad), so a column the set merged is not reported.
/// That is far below the nearest two columns measured apart (`0.0169` of a step, `2.6e-5` rad at the
/// finest lattice; a crossing's column and its neighbour at least `1.4e-4` rad, `3.9e-3` of a step,
/// in the suite's crossing fixtures) — but not provably: two columns laid closer than it would be read as one, and a
/// curve missing one of them would go unreported.
#[cfg(any(test, feature = "test-util"))]
pub fn unshared_columns(model: &Model, t: &Tessellation) -> Vec<Unshared> {
    use std::f64::consts::TAU;
    const WITHIN: f64 = TAU / COARSEST_CIRCLE as f64 * ON_RULING;
    let gap = |a: f64, b: f64| {
        let d = (a - b).rem_euclid(TAU);
        d.min(TAU - d)
    };
    let mut faces: Vec<Handle<Face>> = t.by_face.keys().copied().collect();
    faces.sort_by_key(|f| f.index());
    let mut out = Vec::new();
    for fh in faces {
        let face = model.face(fh);
        let Some(cyl) = model.surface_cache(face.surface).cylinder() else {
            continue;
        };
        // Each curve's samples' angles, ascending, and the span it turns CCW from its first (none
        // for a whole rim).
        struct Samples {
            edge: Handle<Edge>,
            at: Vec<f64>,
            span: Option<(f64, f64)>,
        }
        let mut curves: Vec<Samples> = Vec::new();
        for he in loops_of(face).flat_map(|lp| &lp.half_edges) {
            if !matches!(model.edge_curve(he.edge), Curve::Ellipse(_)) {
                continue;
            }
            let poly = &t.by_edge[&he.edge];
            let at: Vec<f64> = poly
                .iter()
                .map(|&h| in_turn(angle_about(cyl, t.vertices.get(h).pos)))
                .collect();
            let [v0, v1] = model.edge(he.edge).vertices;
            let span = (v0 != v1).then(|| {
                let (a, b) = (at[0], at[at.len() - 1]);
                (a, (b - a).rem_euclid(TAU))
            });
            let mut sorted = at;
            sorted.sort_by(f64::total_cmp);
            curves.push(Samples {
                edge: he.edge,
                at: sorted,
                span,
            });
        }
        for (i, curve) in curves.iter().enumerate() {
            let (samples, span) = (&curve.at, curve.span);
            for (j, other) in curves.iter().enumerate() {
                if i == j {
                    continue;
                }
                for &theta in &other.at {
                    let inside = span.is_none_or(|(a, d)| {
                        let off = (theta - a).rem_euclid(TAU);
                        off > WITHIN && off < d - WITHIN
                    });
                    if !inside {
                        continue;
                    }
                    // The nearest sample about the circle is the one either side of `theta`, the
                    // ascending list wrapping round.
                    let k = samples.partition_point(|&s| s < theta);
                    let n = samples.len();
                    let near = [(k + n - 1) % n, k % n]
                        .into_iter()
                        .any(|i| gap(samples[i], theta) <= WITHIN);
                    if !near {
                        out.push(Unshared {
                            face: fh,
                            edge: curve.edge,
                            theta,
                        });
                    }
                }
            }
        }
    }
    out
}

/// Two consecutive chords of a curved edge's polyline that turn by more than the budget, and by how
/// much: the edge, the chord pair's middle sample (its place in the polyline), and the turn in
/// degrees — `None` when a chord has no length, so no direction. Test-only (`test-util`).
#[cfg(any(test, feature = "test-util"))]
#[derive(Clone, Copy, Debug)]
pub struct Kinked {
    pub edge: Handle<Edge>,
    pub at: usize,
    pub degrees: Option<f64>,
}

/// Every pair of consecutive chords in a curved edge's polyline (`by_edge`) that turns by more than
/// `limit_deg` — the closing chord of a whole rim with both its neighbours, an arc's chords only
/// among themselves (the turn across an edge's end vertex is the model's corner, not a chord's).
/// Test-only (`test-util`): the mesh locks call it with the budget's `max_angle_deg`.
///
/// **What it holds the sampler to.** A curve's polyline must keep the angular budget along the
/// curve itself, and nothing else asks: `within_budget` and [`leaning_laterals`] read surface
/// normals, which a boundary chord keeps by its θ step whatever its curve does. The slack is the
/// mesher's own ([`BUDGET_SLACK`]).
#[cfg(any(test, feature = "test-util"))]
pub fn kinked_edges(model: &Model, t: &Tessellation, limit_deg: f64) -> Vec<Kinked> {
    let limit = limit_deg.to_radians() * BUDGET_SLACK;
    let mut out = Vec::new();
    // `by_edge` is a `HashMap`; walk the edges in handle order so a report reads the same each run.
    let mut edges: Vec<_> = t.by_edge.iter().collect();
    edges.sort_by_key(|(e, _)| e.index());
    for (&edge, poly) in edges {
        if !matches!(model.edge_curve(edge), Curve::Ellipse(_)) {
            continue;
        }
        let pts: Vec<Point3> = poly.iter().map(|&h| t.vertices.get(h).pos).collect();
        let closed = {
            let [v0, v1] = model.edge(edge).vertices;
            v0 == v1
        };
        let n = pts.len();
        // A whole rim's polyline does not repeat its vertex: chord `k` runs from `k` to `k + 1`
        // modulo `n`, and every sample is some pair's middle. An arc's middles are its interior.
        let middles = if closed { 0..n } else { 1..n.saturating_sub(1) };
        for at in middles {
            let (a, b, c) = (pts[(at + n - 1) % n], pts[at], pts[(at + 1) % n]);
            let (u, v) = (b - a, c - b);
            let degrees =
                (u.norm() > 0.0 && v.norm() > 0.0).then(|| u.cross(v).norm().atan2(u.dot(v)));
            if degrees.is_none_or(|d| d > limit) {
                out.push(Kinked {
                    edge,
                    at,
                    degrees: degrees.map(f64::to_degrees),
                });
            }
        }
    }
    out
}

/// A chord of a curved edge's polyline that stands off its curve by more than the chord budget,
/// and by how much: the edge, the chord's first sample (its place in the polyline), and the
/// farthest the curve over it stands off the chord. Test-only (`test-util`).
#[cfg(any(test, feature = "test-util"))]
#[derive(Clone, Copy, Debug)]
pub struct Sagging {
    pub edge: Handle<Edge>,
    pub at: usize,
    pub standoff: f64,
}

/// Every chord of a curved edge's polyline (`by_edge`) the curve stands off by more than `cfg`'s
/// chord budget — a whole rim's closing chord included — or none where the budget asks nothing.
/// Test-only (`test-util`): the mesh locks call it beside [`kinked_edges`], which holds the other
/// budget to the same polylines.
///
/// **What it holds the sampler to.** The chord a polyline puts in for its curve must stay within
/// `tol` of it, and the mesh's own check does not ask: `within_budget` reads interior edges
/// against their surface, and a boundary chord's distance from the cylinder is not its distance
/// from the curve. A chord's ends are read as angles on the curve — a sample's own `t`, a vertex's
/// [`Ellipse::angle_of`] — and the curve between them is the shorter way round, a chord being
/// shorter than half a turn. The reading is [`Ellipse::chord_standoff`].
#[cfg(any(test, feature = "test-util"))]
pub fn sagging_edges(model: &Model, t: &Tessellation, cfg: &TessConfig) -> Vec<Sagging> {
    use std::f64::consts::{PI, TAU};
    let Some(limit) = chord_budget(cfg) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    // `by_edge` is a `HashMap`; walk the edges in handle order so a report reads the same each run.
    let mut edges: Vec<_> = t.by_edge.iter().collect();
    edges.sort_by_key(|(e, _)| e.index());
    for (&edge, poly) in edges {
        let Curve::Ellipse(c) = model.edge_curve(edge) else {
            continue;
        };
        let angles: Vec<f64> = poly
            .iter()
            .map(|&h| {
                let v = t.vertices.get(h);
                match v.origin {
                    TessOrigin::OnEdge { t, .. } => t,
                    _ => c.angle_of(v.pos),
                }
            })
            .collect();
        let closed = {
            let [v0, v1] = model.edge(edge).vertices;
            v0 == v1
        };
        let n = angles.len();
        let chords = if closed { n } else { n.saturating_sub(1) };
        for at in 0..chords {
            let (from, to) = (angles[at], angles[(at + 1) % n]);
            let mut span = (to - from).rem_euclid(TAU);
            if span > PI {
                span -= TAU;
            }
            let standoff = c.chord_standoff(from, from + span);
            if standoff > limit {
                out.push(Sagging { edge, at, standoff });
            }
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

/// **Sample one curve edge into a polyline of shared mesh vertices — on its face's θ set**
/// ([`FaceThetas`]): every column of the set strictly between its two ends, the ends being its
/// own two vertices (shared mesh vertices, so the neighbouring faces' rings meet it crack-free).
///
/// The span is read by the set's angle, the one that read the vertices' columns, so the curve's
/// own ends are known by handle rather than by how close they stand: a column holding either end
/// is the end. A whole rim turns once from its vertex. Point 0 is the first vertex and the rest go
/// round from it; an arc ends on its second vertex, a whole rim closes on its first without
/// repeating it. The made columns are absolute — a function of the face's curves, walked from the
/// cylinder's `ref_dir` ([`walk`]) — so on a face with no ellipse a plain rim's samples are where an
/// even cut from its seam vertex puts them.
///
/// Each sample stands at its column's own number and carries it as its `t`
/// ([`FaceThetas::columns`]), so two curves sampled at one column hold one number there, whatever
/// turn each one's range takes it on.
fn sample_curve(
    t: &mut Tessellation,
    vmap: &mut HashMap<Handle<Vertex>, Handle<TessVertex>>,
    model: &Model,
    eh: Handle<Edge>,
    edge: &Edge,
    c: &Ellipse,
    set: &FaceThetas,
) -> Vec<Handle<TessVertex>> {
    let [v0, v1] = edge.vertices;
    debug_assert!(
        c.axis() == set.cyl.axis().direction() && c.ref_dir() == set.cyl.ref_dir(),
        "a curve shares its face's frame: {c:?} on {:?}",
        set.cyl
    );
    let (t0, dt) = set.span(model, edge);
    let mut poly = vec![vertex_of(t, vmap, model, v0)];
    for (_, number, at) in set.columns(t0, t0 + dt) {
        if at.contains(&v0) || at.contains(&v1) {
            continue;
        }
        let h = t.vertices.push(TessVertex {
            pos: c.point_at(number),
            origin: TessOrigin::OnEdge {
                edge: eh,
                t: number,
            },
        });
        poly.push(h);
    }
    if v1 != v0 {
        poly.push(vertex_of(t, vmap, model, v1));
    }
    poly
}

/// The most columns a face's θ set holds — 4,096 a turn. [`walk`] refuses a face whose budget asks
/// for more, and [`arc_midpoint`] stops cutting a face's curves there ([`FaceThetas::full`]) — one
/// statement, read in two places. A circle too large for its `tol` (past `r ≈ 34,000` at `1e-2`)
/// asks more everywhere and is refused ([`TessError::OverBudget`]); an ellipse asks more only
/// where it turns or sags fast — short stretches on a small cylinder, where the angle leads, but
/// on a large one its sag near the major ends grows with the slope (`OverBudget` names both).
const FINEST_CIRCLE: usize = 4096;

/// The fewest segments a whole curve is cut into, whatever the budget — [`per_turn`]'s floor.
const COARSEST_CIRCLE: usize = 8;

/// The ceiling of a θ set's unit count a turn ([`FaceThetas::grid`]): `2⁴⁹`, so a made column's
/// place `u + kG` is an exact integer in `f64` for `|k| < 16`.
const GRID_CEILING: u64 = 1 << 49;

/// Segments for a whole circle of `radius` — [`per_turn`] for a curve that turns as a circle does,
/// or `None` past [`FINEST_CIRCLE`].
fn circle_segments(cfg: &TessConfig, radius: f64) -> Option<usize> {
    let n = per_turn(cfg, radius, 1.0);
    (n <= FINEST_CIRCLE as f64).then_some(n as usize)
}

/// **The pace a curve bending as `(radius, turn)` asks, as chords a turn**: whichever of the two
/// budgets asks for more ([`TessConfig`] says why there are two), for a curve that sags no more
/// than a circle of `radius` and whose tangent turns at most `turn` radians per radian of θ
/// ([`Ellipse::pace_over`] — `(r, 1)` for a circle). A chord no wider than a turn over this keeps
/// both. With no ceiling: `∞` where a budget cannot be met at all.
///
/// * sagitta: `r(1 − cos(π/n)) ≤ tol` ⟹ `n ≥ π / acos(1 − tol/r)` — grows as `√r`, so
///   it is the term that speaks for a curve large in model units;
/// * angle: `turn · 360°/n ≤ Δθ` ⟹ `n ≥ turn · 360/Δθ` — **the same for every radius**, so it is
///   the term that keeps a small circle from being an octagon, and an ellipse's major ends from
///   turning `a/b` times the budget.
///
/// A budget that is not a finite positive number asks nothing ([`asked`]); a curve neither asks
/// of gets the floor, [`COARSEST_CIRCLE`].
///
/// With the default `2°` the angular term leads a circle until `r ≈ 66` (that is `tol` divided by
/// `1 − cos(1°)`), which covers ordinary part sizes; past it the sagitta term takes over.
fn per_turn(cfg: &TessConfig, radius: f64, turn: f64) -> f64 {
    let by_sagitta = asked(cfg.tol).map_or(0.0, |tol| {
        let ratio = tol / radius;
        if ratio >= 2.0 {
            0.0 // the whole circle is already within tol
        } else {
            // A `tol` so small that `1 − ratio` rounds to `1` asks for infinitely many.
            (std::f64::consts::PI / (1.0 - ratio).acos()).ceil()
        }
    });
    let by_angle = asked(cfg.max_angle_deg).map_or(0.0, |deg| (360.0 / deg * turn).ceil());
    by_sagitta.max(by_angle).max(COARSEST_CIRCLE as f64)
}

/// **A face's made columns: every chord of every curve keeps the budget where it lies** — a
/// chord `h` wide in θ, with `(ρ, R)` the most its curves bend over it ([`Ellipse::pace_over`],
/// and the surface's own `(r, 1)`), has `h·R ≤ Δθ` and sags at most `ρ(1 − cos(h/2)) ≤ tol`
/// ([`per_turn`] over the chord's range). The tangent at a sample turns, between its two chords'
/// middles, by at most `R` of each half-step, so consecutive chords turn by at most `Δθ`; and the
/// chord's ends' surface normals part by `h ≤ Δθ`. A chord cut further keeps it — both its width
/// and its bend shrink — so the vertices' and crossings' columns laid among these leave it true.
/// The turn between consecutive chords would not serve as the rule: a column laid inside a chord
/// makes a new pair of chords, and their turn can be up to the whole chord's.
///
/// **A walk from `0`** (the cylinder's `ref_dir`), in units of `grid` (`count·2ˢ` a turn): from
/// each column the widest step `w` its own range allows ([`step_allowed`] over `[u, u + w]` at
/// least `w`). A wider range asks no less, so the steps that fit are every `w` up to the widest,
/// and a bisection between what the column's own point allows (which no step outruns) and a step
/// known to fit finds it; the lower end is always a step measured to fit, so a rounding that
/// bends that order by an ulp still leaves a step inside the budget. Where the next step reaches
/// the turn's end, the last column laid is taken back and the window from the one before to the
/// end is cut into the fewest equal pieces its ask allows — the window is at least a step, so no
/// piece is a sliver of one. On a face with no ellipse every ask is `count`, every step `2ˢ`, and
/// the columns are the even lattice `τk/count`, bit for bit. Where a curve turns or sags faster
/// the steps are shorter there and only there: an ellipse turning `20` times as fast as its
/// circle at its major ends gets 316 samples a rim, against the 3,605 an even cut at its sharpest
/// took (`an_ellipse_is_cut_where_it_turns`).
///
/// A face's set is the face's own, not its cylinder's ([`FaceThetas`]), and the walk reads only
/// the face's curves whole — not their spans — so it is a function of the face's geometry.
///
/// `None` where the set would hold more than [`FINEST_CIRCLE`] columns, or a step narrower than one
/// unit (a budget that cannot be met there).
fn walk(cfg: &TessConfig, curves: &[Ellipse], count: usize, grid: u64) -> Option<Vec<u64>> {
    let allowed = |lo: u64, hi: u64| step_allowed(cfg, curves, count, grid, lo, hi);
    let fits = |u: u64, w: u64| allowed(u, u + w) >= w;
    let mut made = vec![0u64];
    loop {
        let u = *made.last().expect("the walk starts at 0");
        let most = allowed(u, u);
        let w = if fits(u, most) {
            most
        } else {
            let first = allowed(u, u + most);
            let (mut lo, mut hi) = (if fits(u, first) { first } else { 0 }, most);
            while hi - lo > 1 {
                let mid = lo + (hi - lo) / 2;
                if fits(u, mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        if w == 0 {
            return None;
        }
        if u + w >= grid {
            if made.len() > 1 {
                made.pop();
            }
            let from = *made.last().expect("the walk starts at 0");
            let window = grid - from;
            let step = allowed(from, grid);
            if step == 0 {
                return None;
            }
            let pieces = window.div_ceil(step);
            made.extend((1..pieces).map(|i| from + window * i / pieces));
            return (made.len() <= FINEST_CIRCLE).then_some(made);
        }
        made.push(u + w);
        if made.len() > FINEST_CIRCLE {
            return None;
        }
    }
}

/// **The widest step, in units of `grid`, that the ask over `[lo, hi]` allows** — the ask being the
/// most any of `curves` bends there ([`Ellipse::pace_over`] read by [`per_turn`]) and at least the
/// face's floor `count`. [`walk`] lays a step `w` from `u` where this over `[u, u + w]` is at least
/// `w`; a wider range asks no less, so this never grows as `hi` does.
fn step_allowed(
    cfg: &TessConfig,
    curves: &[Ellipse],
    count: usize,
    grid: u64,
    lo: u64,
    hi: u64,
) -> u64 {
    use std::f64::consts::TAU;
    let theta = |u: u64| TAU * u as f64 / grid as f64;
    let ask = curves
        .iter()
        .map(|c| {
            let (radius, turn) = c.pace_over(theta(lo), theta(hi));
            per_turn(cfg, radius, turn)
        })
        .fold(count as f64, f64::max);
    (grid as f64 / ask).floor() as u64
}

/// The made columns' unit for a floor of `count` a turn: `count·2ˢ`, the largest within
/// [`GRID_CEILING`].
fn grid_of(count: usize) -> u64 {
    let mut grid = count as u64;
    while 2 * grid <= GRID_CEILING {
        grid *= 2;
    }
    grid
}

/// A loop's boundary as an ordered ring of shared mesh vertices — [`traced_ring`]'s walk, so a
/// chart position and the polyline place it came from are one numbering. The ring is the edges'
/// sampled polylines, not the half-edges' start vertices: a face triangulated from those start
/// vertices turns each arc into its chord.
fn boundary_ring(t: &Tessellation, lp: &Loop) -> Vec<Handle<TessVertex>> {
    traced_ring(t, lp).0.into_iter().map(|p| p.handle).collect()
}

fn push_tri(t: &mut Tessellation, fh: Handle<Face>, vertices: [Handle<TessVertex>; 3]) {
    let th = t.triangles.push(TessTriangle { vertices, face: fh });
    t.by_face.entry(fh).or_default().push(th);
}

/// A face's boundary in the **chart its surface is meshed in**: `uv[i]` is where `handles[i]`
/// sits, and `rings` indexes into both — the outer ring first, then the holes.
///
/// **One road, one chart per surface.** Every face — planar or curved — is triangulated by
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
/// It exists because the chart may *invent* a point ([`Chart::interior`]), and an
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

/// A planar face's chart: drop the axis the face's own Newell normal is largest along. The
/// projection rule (Newell → drop an axis → repair handedness) lives here alone, and every
/// planar face's sweep goes through it.
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
/// **`r·θ`, not `θ`.** A cylinder is developable, so scaling the angle by the radius makes
/// this an **isometry** — lengths and angles in the chart are lengths and angles on the surface.
/// That is what lets the Lawson pass ([`polygon`]'s flip) improve the mesh *on the surface* and
/// not merely in a distorted picture of it.
///
/// **`v` is the sweep axis, so the handedness repair negates `u` instead of swapping.**
/// Swapping would send the sweep along the axis instead of around it, and its diagonals would
/// then span wide arcs — chords that leave the surface. Negating `u` mirrors the frame just as
/// well and leaves the sweep going around.
///
/// **θ is unwrapped along each loop, never taken absolutely.** A band is bounded by its two rims,
/// each turning once around the axis, and unrolled alone they are two open lines; `cut` joins
/// them at the generator [`cut_seamless_bands`] chose — first rim, its cut point again, up the
/// generator through any hole it crosses, the second rim, its cut point again, and back down
/// ([`joined_band`]) — so each cut point stands at `θ` and at `θ + 2π` and the unrolled band is a
/// rectangle rather than a degenerate line. The other holes are then shifted by whole turns into
/// the outer ring's range so they lie inside it.
///
/// **`u` of a circle's point is the circle's place**, [`station`], not the point's own axial
/// coordinate — so every circle's arc is level in the chart, at the height the lattice lays its row
/// ([`interior_nodes`]). An ellipse's point stands at the curve's own height there
/// ([`curve_place`]) and lays no row — its arc rises, and the lattice treats it as it treats a
/// ruling; a point the ellipse shares with a circle keeps the circle's place, the row's. A point on
/// no curve stands part way along a ruling, where another face's edge breaks it; its `u` is
/// [`ruling_break_place`], and the lattice lays a row there too.
fn cylinder_chart(
    t: &Tessellation,
    model: &Model,
    face: &Face,
    cyl: &Cylinder,
    cut: Option<BandCut>,
    thetas: &FaceThetas,
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
    let mut circles: HashMap<Handle<nacre_topo::Surface>, f64> = HashMap::new();
    let curves = || {
        loops_of(face)
            .flat_map(|lp| &lp.half_edges)
            .filter_map(|he| match model.edge_curve(he.edge) {
                Curve::Ellipse(c) => Some((he.edge, c)),
                Curve::Line(_) => None,
            })
    };
    for (e, c) in curves().filter(|(_, c)| c.is_circle()) {
        let s = station(cyl, c);
        stations.push(s);
        for &h in &t.by_edge[&e] {
            place.insert(h, s);
        }
        for plane in model.edge(e).surfaces {
            circles.insert(plane, s);
        }
    }
    for (e, c) in curves().filter(|(_, c)| !c.is_circle()) {
        for &h in &t.by_edge[&e] {
            place
                .entry(h)
                .or_insert_with(|| curve_place(cyl, c, t.vertices.get(h).pos));
        }
    }
    let mut uv: Vec<polygon::P2> = vec![[0.0, 0.0]; handles.len()];
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for ring in &rings {
        let ring_handles: Vec<Handle<TessVertex>> = ring.iter().map(|&i| handles[i]).collect();
        let (thetas, _) = unwrapped_thetas(t, cyl, &ring_handles);
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for (&i, &theta) in ring.iter().zip(&thetas) {
            let u = place.get(&handles[i]).copied().unwrap_or_else(|| {
                let v = t.vertices.get(handles[i]);
                let u = ruling_break_place(model, face, cyl, &circles, v)
                    .unwrap_or_else(|| (v.pos - o).dot(w_dir));
                stations.push(u);
                u
            });
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
    let interior = interior_nodes(thetas, cyl, &uv, &rings, stations, negated);
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

/// A loop's boundary walked as a ring — each half-edge's oriented polyline in turn, a sample
/// shared by two consecutive edges taken once — with each position's polyline place, and the place
/// the ring's closing step arrives through (the repeat of the first sample, which the ring drops).
/// The step from position `k` to `k + 1` lies on the edge position `k + 1` was taken from;
/// [`arriving`] reads it.
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

/// Where the step from position `k` of a [`traced_ring`] arrives: the next position, the closing
/// one, or — `true` — the first, through the closing step of a ring with no closing position.
///
/// **That last is one closed edge**, whose polyline does not repeat its first sample: the step
/// from its last sample back to the first. A new sample there goes at the polyline's end
/// ([`between`]).
fn arriving(positions: &[Traced], closing: Option<Traced>, k: usize) -> (Traced, bool) {
    match (k + 1 < positions.len(), closing) {
        (true, _) => (positions[k + 1], false),
        (false, Some(cl)) => (cl, false),
        (false, None) => (positions[0], true),
    }
}

/// Where in `poly` a new sample goes to lie on the step arriving at `step` ([`arriving`]): the step
/// runs `poly[at − 1] → poly[at]` walked forward, `poly[at + 1] → poly[at]` walked backward, and
/// a closing step from the last sample back to the first.
fn between(poly: &[Handle<TessVertex>], step: Traced, wraps: bool) -> usize {
    if wraps {
        poly.len()
    } else if step.forward {
        step.at
    } else {
        step.at + 1
    }
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
    let raw = |h: Handle<TessVertex>| angle_about(cyl, t.vertices.get(h).pos);
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

/// **The band-cut pre-pass — a lateral face that wraps the axis gets one generator to be cut
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
/// and [`joined_band`] joins one — is recorded with `None`.
fn cut_seamless_bands(
    t: &mut Tessellation,
    model: &Model,
    live: &[(Handle<Face>, &Face)],
) -> HashMap<Handle<Face>, BandCut> {
    use std::f64::consts::PI;
    let mut cuts = HashMap::new();
    for &(fh, face) in live {
        let Some(cyl) = model.surface_cache(face.surface).cylinder() else {
            continue;
        };
        let loops: Vec<&Loop> = loops_of(face).collect();
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
                    matches!(model.edge_curve(arriving.edge), Curve::Ellipse(_))
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
/// cut at their vertices — the first sample of each polyline ([`sample_curve`] starts a closed
/// curve there). `None` for any other face, and for two vertices at different angles: the cut is
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
/// **Along the axis, by the curve.** Every crossing lies on an arc (a ruling is parallel to the
/// generator, and `c` lies off every sample), and an arc's axial place is its curve's at `c`
/// ([`curve_place`] — a circle's centre, where every arc of one circle agrees bitwise; an
/// ellipse's own height there) rather than the computed point's, which carries rounding. Two
/// points at one place would be two loops touching on the cut.
///
/// **Around it, one number.** Every ring's point stands at the generator's number in `[0, τ)` and
/// carries it as its `t`, as a sample of a column does ([`FaceThetas::columns`]) — so a crossing
/// chord ending there names the same number on every rim ([`arc_midpoint`]).
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
    // Every crossing is found before any insertion moves a polyline index (`Traced::at`).
    let mut found: Vec<(usize, Traced, bool)> = Vec::new(); // (ring, where the step arrives, closing wrap)
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
            let (step, wraps) = arriving(positions, *closing, k);
            found.push((ring, step, wraps));
        }
    }
    // The generator's one number, every ring's point standing at it.
    let number = in_turn(c);
    let mut points: Vec<(usize, Handle<TessVertex>, f64)> = Vec::with_capacity(found.len());
    for (ring, step, wraps) in found {
        let Curve::Ellipse(circle) = model.edge_curve(step.edge) else {
            return None;
        };
        let pos = circle.point_at(number);
        let h = t.vertices.push(TessVertex {
            pos,
            origin: TessOrigin::OnEdge {
                edge: step.edge,
                t: number,
            },
        });
        let poly = t.by_edge.get_mut(&step.edge)?;
        poly.insert(between(poly, step, wraps), h);
        points.push((ring, h, curve_place(cyl, circle, pos)));
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
/// **Which run of a hole the climb takes is forced, not chosen.** The material lies on the left
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
/// **One spelling, four readers** — the band cut orders its points along the generator by it,
/// the chart gives every point of a circle's polyline it as `u`, the lattice lays its rows at
/// it, and a ruling's break stands at its plane's ([`ruling_break_place`]). Read from the curve, never from a point on it: every arc of one circle agrees on it bitwise
/// (see [`interior_nodes`]), while a point's own axial coordinate carries its rounding — on a
/// tilted axis an arc's samples sit a few ulps off their circle's place, and a lattice row computed
/// one way next to a boundary computed the other was a sliver leaning off the lateral — up to 90°,
/// lying in the cap plane.
fn station(cyl: &Cylinder, circle: &Ellipse) -> f64 {
    (circle.center() - cyl.axis().origin()).dot(cyl.axis().direction())
}

/// **Where a point of a rim stands along the axis, read from the curve** — the rim's [`station`]
/// plus the curve's own height at the point's angle ([`Ellipse::z_at`], `0` for a circle, which
/// leaves the station itself). Read from the curve rather than the point, so every point one
/// carrier pair's edges share stands at one `f64`.
fn curve_place(cyl: &Cylinder, curve: &Ellipse, pos: Point3) -> f64 {
    let z = curve.z_at(curve.angle_of(pos));
    let s = station(cyl, curve);
    if z == 0.0 { s } else { s + z }
}

/// **Where a ruling broken part way stands along the axis** — the place of the plane across the
/// axis that breaks it, if one of the point's two planes is (two cannot be: they would be
/// parallel): the [`station`] of that plane's circle on the face where the face has one
/// (`circles`, keyed by the circle edges' carriers — the stored curve, so bit for bit its points'
/// `u`), else of the circle [`Model::derive_edge_curve`] derives for the cylinder and that plane.
/// Only where neither plane is across does a plane oblique to the axis answer, with its ellipse's
/// own height at the point ([`curve_place`]). Which plane is across the truth says — the derived
/// curve's kind is the exact axis relation's — not the order the planes are stored in: an oblique
/// plane read first would put the break at its ellipse's computed height, ulps off the station
/// the other plane's circle and breaks stand at (measured over the workspace suite, the ignored
/// sweep, the census and the kit: every break beside an across plane already read that plane
/// first, so this rule moves no mesh there). A plane along the axis derives a line, which two ends
/// at this one vertex refuse (`Coincident`), so it never answers.
///
/// One number twice is the hazard: two breaks on one plane are one height, and on a tilted axis
/// their own axial coordinates are one ulp apart — `1.0` and `0.9999999999999999` on the D-shaped
/// cylinder of `a_ruling_broken_part_way_lays_a_row_at_its_vertex` turned 37° and 23°, two rows
/// and 230 triangles leaning off the lateral. Read from the plane, both breaks and the plane's
/// circle on the face are one `f64`.
///
/// `None` where the point is not a [`Vertex::Pierce`] of this face's cylinder (measured: never,
/// over the workspace suite, the ignored sweep, the census and the kit) — the caller then reads the
/// point's own coordinate and lays a row there too. A [`Vertex::CylinderPair`] on a ruling two
/// cylinders share would be one; nothing mints it part way along that ruling, because `boolean`
/// refuses a solid with such an edge.
fn ruling_break_place(
    model: &Model,
    face: &Face,
    cyl: &Cylinder,
    circles: &HashMap<Handle<nacre_topo::Surface>, f64>,
    v: &TessVertex,
) -> Option<f64> {
    let TessOrigin::OnVertex(vh) = v.origin else {
        return None;
    };
    let Vertex::Pierce {
        planes, cylinder, ..
    } = *model.vertex(vh)
    else {
        return None;
    };
    if cylinder != face.surface {
        return None;
    }
    let derived = |plane| match model
        .derive_edge_curve([cylinder, plane], [vh, vh], |_| nacre_topo::EdgeGiven::NONE)
    {
        Ok(Curve::Ellipse(c)) => Some(c),
        _ => None,
    };
    planes
        .iter()
        .find_map(|&plane| {
            circles.get(&plane).copied().or_else(|| {
                derived(plane)
                    .filter(Ellipse::is_circle)
                    .map(|c| station(cyl, &c))
            })
        })
        .or_else(|| {
            planes
                .iter()
                .find_map(|&plane| derived(plane).map(|c| curve_place(cyl, &c, v.pos)))
        })
}

/// **The points a cylindrical face's boundary does not supply, and the mesh needs.**
///
/// A band's rims sample θ finely, and every diagonal the sweep draws between them lands on a
/// closely-sampled arc, so the interior is fine for free. It stops being free in two places. Where
/// a face's boundary leaves part of its θ range without arcs, the sweep has only chords across that
/// stretch (measured on a lateral merged across erased seams: chords spanning half a turn). And
/// where a ruling breaks part way between the rims, the break is a boundary point with no arc at
/// its height, and an interior edge from it leaves the angular budget though the rims cover every θ
/// (`OverBudget`, measured on a half bore with a box beside it). (A band between two curves needs
/// no point: the curves stand on common generators — [`FaceThetas`].)
///
/// So the chart — the layer that knows the surface — supplies them, on a **lattice with no free
/// parameter in it**:
///
/// * **around** (`v`): the face's θ set, `r·θ` of its columns ([`FaceThetas::columns`]) — the
///   columns its curves are sampled on, so a row stands on the generators the curves stand on and
///   a band between a row and a curve is spanned as a band between two curves is (a row on the
///   lattice's columns alone leaned its facets against a rim up to 1.16° where a vertex's column
///   stood between two of them, measured). The columns cover the chart's whole `v` range, because
///   a chart's θ is *unwrapped* — a hole rides whole turns from the outer ring, and a single turn's
///   worth of lattice would miss it entirely.
/// * **along** (`u`): each boundary circle's place, [`station`] — the very number the chart gave
///   that circle's points — and each ruling's break, at [`ruling_break_place`]. Read from the
///   curve, never from the sampled points: one rim's points
///   agree on that coordinate mathematically and differ by ulps in `f64`, so deduplicating *those*
///   would turn one rim into 181 stations.
///   Two arcs of one circle share a centre **bitwise**, and that is a guarantee rather than luck:
///   [`Model::push_edge`] derives the curve from the **canonicalized carrier pair alone** — the
///   circle arm never reads the endpoints — so two edges cut from the same circle by the same two
///   surfaces are handed identical inputs, and a break reads the same derivation for the plane that
///   breaks it. The station count is therefore the number of planes across the axis the face's
///   boundary meets, which is why `dedup` may compare `f64` with `==` here.
///
/// **Where the boundary already stands, no lattice point enters the mesh** — a row lies at an
/// arc's height on purpose, so part of every row is boundary, and a point placed there a few ulps
/// into the material is a sliver leaning off the lateral (up to 90°). Two cases, and only one is a
/// tolerance:
///
/// * **On an arc**: the row and the arc's chord are at one `u`, bit for bit, because both are
///   [`station`] — so the point is *on* the chord and the sweep's exact predicates drop it. Nothing
///   to tune. (Read off the points instead, the two differed by ulps on a tilted axis.)
/// * **On a ruling**: a column meets a ruling where it is a vertex's own (the ruling's), or where
///   the user's corner happens to stand at a multiple of the lattice's angle (a wall at `x = 0.5` against `r = 1` is 30°, the fifteenth
///   2° column). That is not one number computed twice — it is two numbers that geometry makes
///   equal, and no one computation gives both — so a column within [`ON_RULING`] of a rising step
///   (a ruling, an ellipse's arc, or a band's cut generator; a circle's arcs are level) is withheld.
///
/// **A station at the chart's own `u` extent is not laid at all** — there every row point lies on
/// a rim's chord or outside the face. Where only circles bound the face those are its lowest and
/// highest stations (a ruling runs from arc to arc, any break between, so the arcs bound the `u`
/// range); an ellipse reaches past its circle's place, so the extent is read off the chart, not off
/// the stations' order. Laid anyway, an ordinary band is offered 358 points and the sweep drops
/// every one (measured): the mesh is the same, the work is not.
fn interior_nodes(
    thetas: &FaceThetas,
    cyl: &Cylinder,
    uv: &[polygon::P2],
    rings: &[Vec<usize>],
    mut stations: Vec<f64>,
    negated: bool,
) -> Vec<polygon::P2> {
    stations.sort_by(f64::total_cmp);
    stations.dedup();
    // The chart's own `u` range — the rims' places where only circles bound it, the ellipses'
    // reach where they do.
    let (ulo, uhi) = uv.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
        (lo.min(p[0]), hi.max(p[0]))
    });
    let u_of = |s: f64| if negated { -s } else { s };
    stations.retain(|&s| ulo < u_of(s) && u_of(s) < uhi);
    if stations.is_empty() {
        return Vec::new();
    }
    let r = cyl.radius();
    let (lo, hi) = uv.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
        (lo.min(p[1]), hi.max(p[1]))
    });
    let columns = thetas.columns(lo / r, hi / r);
    // The boundary's steps that are not level — after the chart's one place per circle, rulings,
    // ellipses' arcs and a band's cut generator rise.
    let rising: Vec<(polygon::P2, polygon::P2)> = rings
        .iter()
        .flat_map(|ring| (0..ring.len()).map(move |k| (ring[k], ring[(k + 1) % ring.len()])))
        .map(|(a, b)| (uv[a], uv[b]))
        .filter(|(a, b)| a[0] != b[0])
        .collect();
    let near = r * thetas.near();
    let mut out = Vec::new();
    for &s in &stations {
        let u = u_of(s);
        let here: Vec<&(polygon::P2, polygon::P2)> = rising
            .iter()
            .filter(|(a, b)| a[0].min(b[0]) - near <= u && u <= a[0].max(b[0]) + near)
            .collect();
        for &(theta, _, _) in &columns {
            let p = [u, r * theta];
            if !here.iter().any(|&&(a, b)| distance_to_step(p, a, b) < near) {
                out.push(p);
            }
        }
    }
    out
}

/// How close two numbers that stand on one generator may stand and still be two, as a share of a
/// column's step — the one tolerance of the lateral's columns, for two numbers no one computation
/// gives both of: geometry made them equal, or rounding split them. A lattice point beside a rising
/// boundary step (a ruling) is not offered ([`interior_nodes`]) and two vertices at one θ — a
/// ruling's two ends among them — are one column ([`FaceThetas`]), both by the face's floor step
/// ([`column_tolerance`]); a made column beside a held column — a vertex's, or one a crossing chord
/// asked for — is not laid, by its own narrower step ([`FaceThetas::made_near`]). Two held columns are never merged by it unless both are vertices': a crossing's
/// column is placed apart on purpose ([`FaceThetas::add`]).
///
/// **Not a knob** — measured over the workspace suite, the census and the ignored sweep, against
/// the gap each reader leaves: beside a step, every point withheld stands within `6.1e-13` steps
/// and every point kept at least `0.0169`; beside a vertex's column, every lattice column withheld
/// within `5.1e-14` and every one kept at least `0.0169`; two vertices merged stand within
/// `4.3e-14` (a ruling's two ends within `2.6e-14`; in the census and the sweep both at exactly one
/// `f64`), and two vertex columns left apart at least `1.43` steps; beside a crossing's column (the
/// suite's crossing fixtures — the census and the sweep have none), no lattice column withheld and
/// every one kept at least `0.031`. The census meshes are the same
/// at `1e-12`, `1e-9`, `1e-6` and `1e-3`; at `1e-1` a mesh fails its budget. `1e-6` sits in the
/// middle of the gap.
///
/// **Its reach in model size.** Rounding parts a ruling's two ends by about `|p|·ε/r` radians
/// (`|p|` the distance from the origin — estimated, not measured), and the tolerance is a fixed
/// share of the step, so far enough out the ends would read as two columns. The facet judgement
/// refuses first: a unit cylinder with a ruling meshes `3e5` from the origin and is `OverBudget`
/// at `1e6` (measured), where the ends part by about `1e-10` rad against a tolerance of `3.5e-8`.
const ON_RULING: f64 = 1e-6;

/// [`ON_RULING`] of the step of a lattice of `count` columns a turn, in radians.
fn column_tolerance(count: usize) -> f64 {
    std::f64::consts::TAU / count as f64 * ON_RULING
}

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
/// `asked` is a planar face's last boundary scan ([`finish_polylines`]), which the sweep takes
/// where it is this chart's.
#[allow(clippy::too_many_arguments)]
fn triangulate_face(
    t: &mut Tessellation,
    model: &Model,
    cfg: &TessConfig,
    fh: Handle<Face>,
    face: &Face,
    cut: Option<BandCut>,
    thetas: Option<&FaceThetas>,
    asked: Option<&polygon::Asked>,
) -> Result<(), TessError> {
    // The cache: a chart is sampled geometry for a mesh, which is what the realization is
    // for. Nothing here classifies — the kind was decided upstream, on the truth.
    let surface = model.surface_cache(face.surface);
    let Chart {
        mut uv,
        mut handles,
        rings,
        map,
        interior,
    } = match surface {
        SurfaceCache::Plane { plane, .. } => planar_chart(t, face, plane)?,
        SurfaceCache::Cylinder { cylinder, .. } => {
            let thetas = thetas.expect("phase 1 gives every live cylindrical face its θ set");
            cylinder_chart(t, model, face, cylinder, cut, thetas)?
        }
    };
    let refs: Vec<&[usize]> = rings.iter().map(|r| r.as_slice()).collect();
    let boundary = handles.len();
    // A vertex shared by two rings is the bridge pre-pass's doing (a curved ring's sample put
    // into the straight edge it touches); one ring passing a vertex twice is a pinched face.
    let revisits = match surface {
        SurfaceCache::Plane { .. } => shared_vertices(&handles, &rings),
        // Only on a plane is a shared handle one point in the chart — a band's cut point is the
        // same handle at two `θ`s. An arm rather than a boolean, so a third surface kind is a
        // compile error here rather than a silent "no bridges".
        SurfaceCache::Cylinder { .. } => Vec::new(),
    };
    let geom = surface.to_geom();
    // The point the mesh will ship for each corner: a boundary vertex's shared position, or the
    // chart's inverse for a point this face mints.
    let at = |(i, p): polygon::Corner| {
        if i < boundary {
            t.vertices.get(handles[i]).pos
        } else {
            map.invert(p).1
        }
    };
    let within = |e: polygon::Element| match e {
        polygon::Element::Edge(c) => edge_in_budget(cfg, &geom, c.map(|k| (k.0, at(k)))),
        polygon::Element::Facet(c) => facet_in_budget(cfg, &geom, c.map(|k| (k.0, at(k)))),
    };
    let (tris, rings_used) =
        polygon::triangulate_uv(&mut uv, &refs, &interior, &revisits, &within, asked)?;
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
    within_budget(t, cfg, &geom, &handles, &rings_used, &tris)?;
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

/// **The mesh's one rule, asked of the face's *interior* and of every facet.**
///
/// `TessConfig` declares two budgets and [`walk`] enforces **both** on every chord of every curved
/// boundary polyline — or refuses the mesh where it cannot. A face's interior otherwise inherits whatever sampling its boundary supplies: fine
/// while the boundary samples the curvature (a band's rims do, on common generators —
/// [`FaceThetas`]), and silently wrong when it does not — a lateral whose boundary leaves a θ
/// stretch unsampled meshes as flat chords across it (measured on one merged across erased seams:
/// four triangles spanning half a turn, a sixth of its area lost, while `validate`, watertightness,
/// the exact volume and the face counts stay green).
///
/// So the same two questions are asked of every interior edge:
/// * **linear** — how far the surface strays from the chord, measured at its midpoint;
/// * **angular** — how far the surface *turns* between the ends, measured on its own normals;
///
/// and one of every facet, interior or not — how far its normal leans off the surface's
/// ([`facet_in_budget`]). Edges within the budget do not hold a facet to its surface: three
/// samples in a row along a curve oblique to a cylinder's axis make a facet in that curve's plane.
///
/// All of it comes from [`nacre_geom::Surface`] and is exhaustive over the surface kinds, so a
/// plane answers zero without anyone declaring that a plane is flat — only a sliver whose normal
/// its points do not determine is outside there (none in the corpus or the suites) — and a new
/// surface kind is a compile error rather than a missing case.
///
/// **A boundary edge cannot fail, so it is not asked.** A curve is cut so that each chord keeps the
/// *larger* of the two budgets' demands where it lies ([`walk`]), and a mesh whose curves would
/// need more than the finest sampling is refused before any face is drawn; a straight edge (a ruling, a
/// plane's side) strays zero and its ends share one normal. So of the edges this checks only what
/// the sweep chose, never what the sampler laid down — but every facet, since a facet with two
/// boundary edges is still the sweep's choice.
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
    let at = |i: usize| (i, t.vertices.get(handles[i]).pos);
    for tri in tris {
        if !facet_in_budget(cfg, surface, tri.map(at)) {
            return Err(TessError::OverBudget);
        }
        for k in 0..3 {
            let (i, j) = (tri[k], tri[(k + 1) % 3]);
            if boundary.contains(&(i.min(j), i.max(j))) {
                continue;
            }
            if !edge_in_budget(cfg, surface, [at(i), at(j)]) {
                return Err(TessError::OverBudget);
            }
        }
    }
    Ok(())
}

/// **The slack every budget comparison allows**, one number for the chord and the angle and for the
/// lock on a curve's turn (`kinked_edges`): today's populations sit **exactly** on the budget — an
/// interior diagonal spanning one sample turns by precisely the sampler's own step, a circle's two
/// chords by precisely the angle — so a bare `>` would be decided by the last bit of two
/// differently-rounded routes to the same number.
const BUDGET_SLACK: f64 = 1.0 + 1e-9;

/// **Whether a budget asks anything**: `budget` itself when it is a finite positive number, `None`
/// otherwise — one rule for both of [`TessConfig`]'s budgets, read by the sampler ([`per_turn`])
/// and by the mesh's check ([`angle_budget`], [`chord_budget`]).
fn asked(budget: f64) -> Option<f64> {
    (budget.is_finite() && budget > 0.0).then_some(budget)
}

/// The angular budget in radians, with [`BUDGET_SLACK`], or `None` when
/// [`TessConfig::max_angle_deg`] asks nothing ([`asked`]) — then nothing is asked of the angle.
fn angle_budget(cfg: &TessConfig) -> Option<f64> {
    asked(cfg.max_angle_deg).map(|deg| deg.to_radians() * BUDGET_SLACK)
}

/// The chord budget, with [`BUDGET_SLACK`], or `None` when [`TessConfig::tol`] asks nothing
/// ([`asked`]) — then nothing is asked of how far an edge strays.
fn chord_budget(cfg: &TessConfig) -> Option<f64> {
    asked(cfg.tol).map(|tol| tol * BUDGET_SLACK)
}

/// **Whether the edge `a — b` keeps both budgets on `surface`**: its midpoint within `tol` of it,
/// and the surface's normals at its ends within the angular budget of each other — each where its
/// budget asks ([`asked`]).
///
/// This and [`facet_in_budget`] are the one spelling of the question [`within_budget`] asks of the
/// finished mesh and the flip pass asks before it flips (`polygon::triangulate_uv`'s `within`).
/// Each corner comes with its index and the point the mesh will ship for it, and the corners are
/// put in index order first, so an element is judged in the same rounding whichever triangle or
/// flip asks about it — which is what lets the flip pass count.
fn edge_in_budget(cfg: &TessConfig, surface: &Surface, mut ends: [(usize, Point3); 2]) -> bool {
    ends.sort_unstable_by_key(|&(i, _)| i);
    let [a, b] = ends.map(|(_, p)| p);
    if chord_budget(cfg).is_some_and(|tol| surface.distance(a + (b - a) * 0.5) > tol) {
        return false;
    }
    let Some(limit) = angle_budget(cfg) else {
        return true;
    };
    let (Some(na), Some(nb)) = (surface.normal_at(a), surface.normal_at(b)) else {
        return true; // a normal with no name (a point on the axis) is not a budget failure
    };
    na.dot(nb).clamp(-1.0, 1.0).acos() <= limit
}

/// **Whether the facet `a b c` keeps the angular budget** (in [`edge_in_budget`]'s spelling): the
/// angle between its normal and the surface's normal at its centroid, plus how far rounding may
/// have turned its normal, within the budget. A facet whose normal its own points do not determine
/// — no area, or a sliver whose cross product is rounding — is outside.
///
/// **The rounding.** Each component of `d₁ × d₂`, `d₁ = b − a` and `d₂ = c − a` rounded, is a
/// difference of two products of rounded differences, so it is off by at most
/// `4u(|d₁ᵧ||d₂_z| + |d₁_z||d₂ᵧ|)` to first order (`u = ε/2`: two differences, a product, a
/// subtraction); over the three components that is at most `4√2·u|d₁||d₂| < 4ε|d₁||d₂|`. A vector
/// of length `|n|` moved by `|Δ| < |n|` turns by at most `asin(|Δ|/|n|)`. The points are the mesh's
/// own, so only this arithmetic is counted; the centroid, the surface's normal, the dot product
/// and the `acos` are rounded as [`edge_in_budget`]'s are, inside the budget's slack. A centroid on
/// the axis has no normal to compare with; it is outside — the facet's projection across the axis
/// would have to hold the axis, and then its edges are far outside the budget anyway.
///
/// **No chord budget here, and why none is missing.** A plane's facet lies in it. On a cylinder, a
/// facet whose corners span less than a half-turn projects across the axis to a triangle that does
/// not hold the axis, so its point nearest the axis — farthest from the surface — is on its
/// boundary, and on each edge, whose ends are both on the cylinder, at the edge's midpoint: which
/// [`edge_in_budget`] asks of an interior edge and the sampler lays down on a boundary one. A new
/// surface kind is where to ask again.
///
/// Both arguments — the centroid on the axis, and no chord budget — rest on the angular budget
/// holding every edge to a few degrees. Where nothing is asked of the angle ([`angle_budget`])
/// neither is asked, and a facet around the axis, each edge a third of a turn, passes a chord
/// budget of `r/2` with its centroid `r` off the surface — a pathological setting.
fn facet_in_budget(cfg: &TessConfig, surface: &Surface, mut corners: [(usize, Point3); 3]) -> bool {
    corners.sort_unstable_by_key(|&(i, _)| i);
    let [a, b, c] = corners.map(|(_, p)| p);
    let (d1, d2) = (b - a, c - a);
    let n = d1.cross(d2);
    let len = n.norm();
    let unsure = 4.0 * f64::EPSILON * d1.norm() * d2.norm() / len;
    if !unsure.is_finite() || unsure >= 1.0 {
        return false;
    }
    let Some(limit) = angle_budget(cfg) else {
        return true;
    };
    let Some(normal) = surface.normal_at(a + (d1 + d2) * (1.0 / 3.0)) else {
        return false;
    };
    (n * (1.0 / len)).dot(normal).abs().clamp(0.0, 1.0).acos() + unsure.asin() <= limit
}

#[cfg(test)]
#[path = "tests/lib.rs"]
mod tests;
