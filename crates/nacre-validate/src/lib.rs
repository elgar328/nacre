//! b-rep invariant checker for the nacre kernel (design.md §7).
//!
//! [`validate`] runs every M1 check over a [`Model`] and returns all
//! [`Violation`]s it finds (an empty `Vec` means the model is valid). Checks:
//! reference integrity, loop closure, half-edge manifold pairing, geometric
//! incidence (vertices on their curves/surfaces), and Euler-Poincaré. The
//! tessellation checks (§5, §7 — provenance coherence, crack-free) arrive in M3
//! when a `Tessellation` exists.

use nacre_geom::Surface;
use nacre_math::{Point3, Vector3};
use nacre_store::{Handle, Store};
use nacre_topo::{
    Adjacency, Edge, Face, Loop, Model, Orientation, Reachable, Shell, Solid, Vertex, VertexDef,
};

/// Residual bound for a vertex with **no measured tolerance** lying on its reference
/// curve/surface. Machine epsilon (~2.2e-16) is too tight — such a coordinate is the
/// output of a short floating-point construction chain (corner arithmetic, line/plane
/// fitting), so its residual against its own fitted geometry is ~magnitude·(a few ULP)
/// ≈ 1e-13 for coordinates up to ~1e3. `1e-9` sits comfortably above that floor yet
/// ~7 orders below any measured tolerance.
pub const EPS_CONSTRUCTED: f64 = 1e-9;

/// Which reference edge in the topology graph a dangling handle sits on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefKind {
    EdgeBoundVertex,
    FaceSurface,
    HalfEdgeEdge,
    ShellFace,
    SolidShell,
    /// A vertex definition's surface handle.
    VertexDefSurface,
}

/// Which loop of a face a defect was found in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LoopKind {
    Outer,
    Inner(usize),
}

/// A single broken invariant. [`validate`] returns every one it finds.
///
/// The f64 residual/tolerance and `Point3` fields mean this is `PartialEq` but
/// not `Eq`/`Hash` — identity of the offending element is carried by
/// `Handle`/index, not by value.
#[derive(Clone, Debug, PartialEq)]
pub enum Violation {
    /// A handle indexes at or past the end of its target store
    /// (`target_index >= target_len`). Guards hand-built / future deserialized
    /// models; a `push`-built model can never produce one. Type-erased to `u32`
    /// because the six [`RefKind`]s point at six different `Handle<T>` types.
    DanglingReference {
        kind: RefKind,
        owner_index: u32,
        target_index: u32,
        target_len: u32,
    },

    /// A loop's half-edge chain does not close: `end(he[at]) != start(he[at+1])`
    /// (indices mod loop length).
    OpenLoop {
        face: Handle<Face>,
        loop_kind: LoopKind,
        at: usize,
    },

    /// An edge is used by a number of half-edges other than two. A closed
    /// 2-manifold uses every edge exactly twice. `use_count == 0` = orphan;
    /// `1` = boundary/open surface; `>= 3` = non-manifold.
    NonManifoldEdge {
        edge: Handle<Edge>,
        use_count: usize,
    },

    /// An edge is used exactly twice but with the *same* `forward` flag — the
    /// two faces traverse it in the same direction (inconsistent orientation).
    NonOpposedEdge {
        edge: Handle<Edge>,
        faces: [Handle<Face>; 2],
    },

    /// A vertex whose surface link is not a single circle — a non-manifold "pinch" where two or
    /// more face-fans meet at one point (two solids touching only at a corner), even though every
    /// edge is manifold. Detected by [`nacre_topo::nonmanifold_vertices`].
    NonManifoldVertex { vertex: Handle<Vertex> },

    /// A vertex definition whose carrier kinds contradict its variant (S7): a `ThreePlane`
    /// naming a cylinder (three *planes* is the claim), or an `OnSeam` naming two planes
    /// (two planes meet in a line — a line's point IS a three-plane intersection, so the
    /// seam spelling would be hiding an expressible truth). The variants' invariants are
    /// per-variant (Q5), and this is the checker that keeps them so — the vertex sibling of
    /// [`Self::EdgeCarrierMismatch`].
    VertexDefCarrierMismatch { vertex: Handle<Vertex> },

    /// An edge's stated carriers disagree with adjacency (S8): the multiset of the two face
    /// surfaces using the edge is not the stored `Edge::surfaces` pair — or the pair is
    /// self-adjacent (`[s, s]`) on a *plane*, a spelling reserved for a cylinder seam. The
    /// carriers are stated, never derived, so a mismatch is a producer bug, and the curve
    /// cache derived from wrong carriers would be silently wrong geometry.
    EdgeCarrierMismatch {
        edge: Handle<Edge>,
        stated: [u32; 2],
        observed: [u32; 2],
    },

    /// A bound vertex does not lie on its edge's curve within tolerance.
    VertexOffCurve {
        edge: Handle<Edge>,
        vertex: Handle<Vertex>,
        point: Point3,
        residual: f64,
        tol: f64,
    },

    /// A loop vertex does not lie on its face's surface within tolerance.
    VertexOffSurface {
        face: Handle<Face>,
        vertex: Handle<Vertex>,
        point: Point3,
        residual: f64,
        tol: f64,
    },

    /// A **measured** vertex does not lie on one of its definition's surfaces within its
    /// measured tolerance — the definition is the truth, so the cached point must sit within
    /// `tol` of every surface it is defined as meeting (design §4). `surface_index` is
    /// type-erased (like [`Self::DanglingReference`]) so the checker never names geom's
    /// `Surface` (geom stays a dev-dependency).
    VertexOffDefinition {
        vertex: Handle<Vertex>,
        surface_index: u32,
        point: Point3,
        residual: f64,
        tol: f64,
    },

    /// `chi = V - E + F - L_i` is odd, so `2(S - G) = chi` has no integer
    /// solution: the boundary cannot be a valid closed 2-manifold.
    EulerParity {
        v: usize,
        e: usize,
        f: usize,
        s: usize,
        inner_loops: usize,
    },

    /// `chi` is even but the implied genus `G = S - chi/2` is negative — more
    /// handles than topologically possible (e.g. disjoint closed surfaces
    /// grouped under a single shell).
    NegativeGenus {
        v: usize,
        e: usize,
        f: usize,
        s: usize,
        inner_loops: usize,
        genus: i64,
    },

    /// A cavity (inner void) shell's faces do not point their outward normals
    /// into the void: its signed self-volume is `≥ 0` (a correct inward void is
    /// negative). Such a cavity would *add* to the solid's volume instead of
    /// subtracting it (design §8 M5 containment). The other checks miss it — a
    /// globally-reversed shell keeps edge opposition (so `check_manifold`
    /// passes) and the same V/E/F/S (so Euler passes). Planar cavities only; a
    /// non-planar void is not checked here.
    CavityMisoriented {
        solid: Handle<Solid>,
        cavity: Handle<Shell>,
        signed_volume: f64,
    },
}

/// Check every M1 invariant of `model`, returning all violations (empty = valid).
///
/// Reference integrity runs first and short-circuits: it uses only `.index()` /
/// `.len()` (never `Store::get`), so it is safe on a corrupt model; the later
/// checks dereference handles via `get`, which would panic on a dangling one.
/// The adjacency cache is rebuilt fresh here, so callers need not have called
/// [`Model::rebuild_adjacency`].
pub fn validate(model: &Model) -> Vec<Violation> {
    let mut out = Vec::new();

    check_reference_integrity(model, &mut out);
    if !out.is_empty() {
        return out;
    }

    // The live model, not the whole append-only arena (design §2): superseded
    // cells stay in the store but drop out here, so they neither break Euler nor
    // pollute manifold use-counts. Reference integrity ran first (and would have
    // short-circuited on a dangling handle), so this traversal is in-bounds.
    let reach = model.reachable();
    let adj = Adjacency::rebuild(model); // fresh; does not trust model.adj
    check_vertex_def_carriers(model, &mut out);
    check_loop_closure(model, &reach, &mut out);
    check_manifold(model, &adj, &reach, &mut out);
    check_cavity_orientation(model, &mut out);
    check_geometric_incidence(model, &reach, &mut out);
    check_euler_poincare(model, &reach, &mut out);
    out
}

/// The tolerance a vertex's provenance grants (S7: read from the point cache): a measured
/// tolerance where one exists (`Some` — a discovered vertex, `0.0` kept exact), else the
/// construction epsilon [`EPS_CONSTRUCTED`].
#[inline]
fn tol_of(m: &Model, vh: Handle<Vertex>) -> f64 {
    m.vertex_tol(vh).unwrap_or(EPS_CONSTRUCTED)
}

#[inline]
fn in_bounds<T>(h: Handle<T>, store: &Store<T>) -> bool {
    (h.index() as usize) < store.len()
}

fn check_reference_integrity(m: &Model, out: &mut Vec<Violation>) {
    for (eh, edge) in m.edges.iter() {
        // (The curve handle's own check died with `Store<Curve>` (S8): the curve is a cache
        // beside the store, not a reference an edge can dangle.)
        {
            for v in edge.vertices {
                if !in_bounds(v, &m.vertices) {
                    out.push(Violation::DanglingReference {
                        kind: RefKind::EdgeBoundVertex,
                        owner_index: eh.index(),
                        target_index: v.index(),
                        target_len: m.vertices.len() as u32,
                    });
                }
            }
        }
    }

    for (fh, face) in m.faces.iter() {
        // The surfaces store is private (S1) — bounds-check against its count.
        if face.surface.index() as usize >= m.surface_count() {
            out.push(Violation::DanglingReference {
                kind: RefKind::FaceSurface,
                owner_index: fh.index(),
                target_index: face.surface.index(),
                target_len: m.surface_count() as u32,
            });
        }
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                if !in_bounds(he.edge, &m.edges) {
                    out.push(Violation::DanglingReference {
                        kind: RefKind::HalfEdgeEdge,
                        owner_index: fh.index(),
                        target_index: he.edge.index(),
                        target_len: m.edges.len() as u32,
                    });
                }
            }
        }
    }

    for (sh, shell) in m.shells.iter() {
        for f in &shell.faces {
            if !in_bounds(*f, &m.faces) {
                out.push(Violation::DanglingReference {
                    kind: RefKind::ShellFace,
                    owner_index: sh.index(),
                    target_index: f.index(),
                    target_len: m.faces.len() as u32,
                });
            }
        }
    }

    for (soh, solid) in m.solids.iter() {
        for sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            if !in_bounds(*sh, &m.shells) {
                out.push(Violation::DanglingReference {
                    kind: RefKind::SolidShell,
                    owner_index: soh.index(),
                    target_index: sh.index(),
                    target_len: m.shells.len() as u32,
                });
            }
        }
    }

    // Every vertex's definition references surfaces by handle (S7: the definition is the
    // vertex, so this covers all of them, not just the discovered population).
    for (vh, vertex) in m.vertices.iter() {
        let surfaces: &[Handle<Surface>] = match &vertex.def {
            VertexDef::ThreePlane(s) => s,
            VertexDef::OnSeam(s) => s,
        };
        for &s in surfaces {
            if s.index() as usize >= m.surface_count() {
                out.push(Violation::DanglingReference {
                    kind: RefKind::VertexDefSurface,
                    owner_index: vh.index(),
                    target_index: s.index(),
                    target_len: m.surface_count() as u32,
                });
            }
        }
    }
}

fn check_euler_poincare(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Count the live model only (design §2), not the append-only store lengths.
    let v = reach.vertices.len();
    let e = reach.edges.len();
    let f = reach.faces.len();
    let s = reach.shells.len();
    let inner_loops: usize = reach
        .faces
        .iter()
        .map(|fh| m.faces.get(*fh).inner.len())
        .sum();

    // i64: E can exceed V + F. V - E + F - L_i = 2(S - G) for a closed 2-manifold.
    let chi = v as i64 - e as i64 + f as i64 - inner_loops as i64;
    if chi % 2 != 0 {
        out.push(Violation::EulerParity {
            v,
            e,
            f,
            s,
            inner_loops,
        });
    } else {
        let genus = s as i64 - chi / 2;
        if genus < 0 {
            out.push(Violation::NegativeGenus {
                v,
                e,
                f,
                s,
                inner_loops,
                genus,
            });
        }
    }
}

/// S7 structural check: each vertex definition's carrier kinds must match its variant —
/// `ThreePlane` names planes only, `OnSeam` includes a non-plane. Runs after reference
/// integrity (it dereferences surface handles).
fn check_vertex_def_carriers(m: &Model, out: &mut Vec<Violation>) {
    for (vh, vertex) in m.vertices.iter() {
        let bad = match &vertex.def {
            VertexDef::ThreePlane(planes) => planes
                .iter()
                .any(|&s| !matches!(m.surface(s), Surface::Plane(_))),
            VertexDef::OnSeam(pair) => pair
                .iter()
                .all(|&s| matches!(m.surface(s), Surface::Plane(_))),
        };
        if bad {
            out.push(Violation::VertexDefCarrierMismatch { vertex: vh });
        }
    }
}

fn check_loop_closure(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Store order (deterministic), filtered to the live faces.
    for (fh, face) in m.faces.iter() {
        if !reach.faces.contains(&fh) {
            continue;
        }
        check_loop(m, fh, LoopKind::Outer, &face.outer, out);
        for (i, lp) in face.inner.iter().enumerate() {
            check_loop(m, fh, LoopKind::Inner(i), lp, out);
        }
    }
}

fn check_loop(m: &Model, fh: Handle<Face>, kind: LoopKind, lp: &Loop, out: &mut Vec<Violation>) {
    let hes = &lp.half_edges;
    let n = hes.len();
    if n == 0 {
        out.push(Violation::OpenLoop {
            face: fh,
            loop_kind: kind,
            at: 0,
        });
        return;
    }
    // Resolve each half-edge to (start, end). (Every edge is bounded by type since S8 —
    // the `UnboundedEdgeInLoop` arm died with the `Option`.)
    let ends: Vec<Option<(Handle<Vertex>, Handle<Vertex>)>> = hes
        .iter()
        .map(|he| {
            let [a, b] = m.edges.get(he.edge).vertices;
            Some(if he.forward { (a, b) } else { (b, a) })
        })
        .collect();

    for i in 0..n {
        if let (Some((_, end)), Some((start, _))) = (ends[i], ends[(i + 1) % n]) {
            if end != start {
                out.push(Violation::OpenLoop {
                    face: fh,
                    loop_kind: kind,
                    at: i,
                });
            }
        }
    }
}

fn check_manifold(m: &Model, adj: &Adjacency, reach: &Reachable, out: &mut Vec<Violation>) {
    // Live edges only (design §2): a superseded edge left in the store has 0 uses
    // in the live adjacency but is not a defect — skip it. Every reachable edge
    // is referenced by a live face, so it must be used exactly twice.
    for (eh, _edge) in m.edges.iter() {
        if !reach.edges.contains(&eh) {
            continue;
        }
        let uses = adj.edge_uses.get(&eh).map(Vec::as_slice).unwrap_or(&[]);
        if uses.len() != 2 {
            out.push(Violation::NonManifoldEdge {
                edge: eh,
                use_count: uses.len(),
            });
        } else {
            if uses[0].1 == uses[1].1 {
                out.push(Violation::NonOpposedEdge {
                    edge: eh,
                    faces: [uses[0].0, uses[1].0],
                });
            }
            // S8: the stated carrier pair must be the two using faces' surfaces (as a
            // multiset), and `[s, s]` on a plane is a spelling reserved for cylinder seams.
            let stated = _edge.surfaces;
            let mut observed = [
                m.faces.get(uses[0].0).surface,
                m.faces.get(uses[1].0).surface,
            ];
            if observed[1].index() < observed[0].index() {
                observed.swap(0, 1);
            }
            let plane_self_pair = stated[0] == stated[1]
                && matches!(m.surface(stated[0]), nacre_geom::Surface::Plane(_));
            if stated != observed || plane_self_pair {
                out.push(Violation::EdgeCarrierMismatch {
                    edge: eh,
                    stated: [stated[0].index(), stated[1].index()],
                    observed: [observed[0].index(), observed[1].index()],
                });
            }
        }
    }
    // Non-manifold *vertices* (pinch points) — every edge can be manifold yet two face-fans meet
    // at one vertex. `adj` is already reachable-scoped (rebuilt fresh), so every returned vertex is
    // live.
    for vertex in nacre_topo::nonmanifold_vertices(&adj.vertex_edges, &adj.edge_uses) {
        out.push(Violation::NonManifoldVertex { vertex });
    }
}

/// Each live solid's cavity (void) shells must be inward-oriented: the shell's
/// signed self-volume (via each face's `plane.normal() × orientation`, the same
/// normal `nacre-props` integrates) must be negative. A `≥ 0` value means the
/// void points outward and would add to the solid's volume. Runs after the
/// reference-integrity short-circuit, so every dereferenced handle is in bounds.
fn check_cavity_orientation(m: &Model, out: &mut Vec<Violation>) {
    for &sh in &m.live_solids {
        if !in_bounds(sh, &m.solids) {
            continue;
        }
        for &cavity in &m.solids.get(sh).cavities {
            if let Some(v) = shell_signed_volume(m, cavity) {
                if v >= 0.0 {
                    out.push(Violation::CavityMisoriented {
                        solid: sh,
                        cavity,
                        signed_volume: v,
                    });
                }
            }
        }
    }
}

/// Signed volume of a planar shell via the divergence theorem, about a
/// shell-local reference `R` (the closed-surface identity `∮ n̂ dA = 0` makes it
/// `R`-independent; a local `R` avoids the cancellation an origin-far shell
/// would suffer). Outward-oriented shell ⇒ `+V`; a correct inward cavity ⇒ `−V`.
/// `None` if any face is non-planar or a loop is unbounded (M5 cavities are
/// planar; a curved void is simply not checked). Mirrors `nacre-props`'
/// `face_contribution`, duplicated to keep validate off the props layer.
fn shell_signed_volume(m: &Model, shell: Handle<Shell>) -> Option<f64> {
    let faces = &m.shells.get(shell).faces;
    let first = m.faces.get(*faces.first()?);
    let reference = m.vertex_point(loop_start(m, &first.outer)?);

    let mut flux = 0.0;
    for &fh in faces {
        let face = m.faces.get(fh);
        let Surface::Plane(plane) = m.surface(face.surface) else {
            return None;
        };
        let sign = match face.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        let normal = plane.normal() * sign;
        let (area, centroid) = loop_area_centroid(m, &face.outer)?;
        flux += normal.dot(centroid - reference) * area;
        for hole in &face.inner {
            let (a_in, c_in) = loop_area_centroid(m, hole)?;
            flux -= normal.dot(c_in - reference) * a_in;
        }
    }
    Some(flux / 3.0)
}

/// The start vertex of a loop's first half-edge (`None` for an empty loop).
fn loop_start(m: &Model, lp: &Loop) -> Option<Handle<Vertex>> {
    let he = lp.half_edges.first()?;
    let [a, b] = m.edges.get(he.edge).vertices;
    Some(if he.forward { a } else { b })
}

/// `(unsigned area, area-weighted centroid)` of a planar polygon loop, via a
/// signed triangle fan from the first vertex (exact for concave loops). `None`
/// if the loop has an unbounded edge or fewer than three vertices.
fn loop_area_centroid(m: &Model, lp: &Loop) -> Option<(f64, Point3)> {
    let pts: Vec<Point3> = lp
        .half_edges
        .iter()
        .map(|he| {
            let [a, b] = m.edges.get(he.edge).vertices;
            m.vertex_point(if he.forward { a } else { b })
        })
        .collect();
    if pts.len() < 3 {
        return None;
    }
    let base = pts[0];
    let mut area_vec = Vector3::zero();
    for w in pts[1..].windows(2) {
        area_vec += (w[0] - base).cross(w[1] - base);
    }
    let unit = area_vec.normalize()?;
    let mut weighted = Vector3::zero();
    let mut weight = 0.0;
    for w in pts[1..].windows(2) {
        let signed = (w[0] - base).cross(w[1] - base).dot(unit);
        let centroid_rel = ((w[0] - base) + (w[1] - base)) * (1.0 / 3.0);
        weighted += centroid_rel * signed;
        weight += signed;
    }
    Some((0.5 * area_vec.norm(), base + weighted * (1.0 / weight)))
}

fn check_geometric_incidence(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Each live edge's bound vertices must lie on the edge's curve.
    for (eh, edge) in m.edges.iter() {
        if !reach.edges.contains(&eh) {
            continue;
        }
        {
            let [a, b] = edge.vertices;
            let curve = m.edge_curve(eh);
            for vh in [a, b] {
                let residual = curve.distance(m.vertex_point(vh));
                // `.max(EPS_CONSTRUCTED)` is what `.max(tol_of(edge.origin))` always evaluated
                // to (every producer wrote `Constructed`), spelled as the constant it was after
                // `Edge.origin` died (S8). It is NOT redundant with the vertex term: a
                // `Discovered { tol: 0.0 }` vertex (an exact-zero measured residual — a real
                // population) relies on this floor to absorb the distance computation's own
                // machine-scale noise.
                let tol = tol_of(m, vh).max(EPS_CONSTRUCTED);
                if residual > tol {
                    out.push(Violation::VertexOffCurve {
                        edge: eh,
                        vertex: vh,
                        point: m.vertex_point(vh),
                        residual,
                        tol,
                    });
                }
            }
        }
    }
    // Each loop vertex must lie on the face's surface (Face/Surface carry no
    // Origin, so only the vertex's provenance relaxes the bound). Unbounded
    // edges have no start vertex here and are skipped (already flagged by
    // loop-closure).
    for (fh, face) in m.faces.iter() {
        if !reach.faces.contains(&fh) {
            continue;
        }
        let surface = m.surface(face.surface);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                {
                    let [a, b] = m.edges.get(he.edge).vertices;
                    let vh = if he.forward { a } else { b };
                    let residual = surface.distance(m.vertex_point(vh));
                    // ★ Plus what the residual's *own* arithmetic can produce. `tol_of` describes
                    // where the vertex may sit; it says nothing about `Surface::distance`, so a
                    // vertex exactly on the surface can still report a machine-scale residual and
                    // be flagged for it. Measured: a four-plane concurrency whose vertex is
                    // genuinely on the plane sat at residual *exactly equal* to its tolerance, and
                    // passed only because the comparison is strict — the slack that had been
                    // covering this term was the loose `Discovered` tolerances of the day.
                    let tol = tol_of(m, vh) + surface.distance_eps(m.vertex_point(vh));
                    if residual > tol {
                        out.push(Violation::VertexOffSurface {
                            face: fh,
                            vertex: vh,
                            point: m.vertex_point(vh),
                            residual,
                            tol,
                        });
                    }
                }
            }
        }
    }
    // Each live discovered vertex must lie within its tol of every plane its
    // definition claims it is the intersection of (design §4 — the definition is
    // the truth, the point is a within-tol cache). Reference integrity ran first,
    // so the definition's surface handles are in bounds.
    // (Population letter-preserved in the S7 swap: measured — tol-Some — vertices only, the
    // old `Discovered` set. The all-vertices extension is its own commit.)
    for (vh, vertex) in m.vertices.iter() {
        if !reach.vertices.contains(&vh) {
            continue;
        }
        let Some(tol) = m.vertex_tol(vh) else {
            continue;
        };
        let surfaces: &[Handle<Surface>] = match &vertex.def {
            VertexDef::ThreePlane(s) => s,
            VertexDef::OnSeam(s) => s,
        };
        for &sh in surfaces {
            let residual = m.surface(sh).distance(m.vertex_point(vh));
            if residual > tol {
                out.push(Violation::VertexOffDefinition {
                    vertex: vh,
                    surface_index: sh.index(),
                    point: m.vertex_point(vh),
                    residual,
                    tol,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_geom::{Plane, Surface};
    use nacre_math::Vector3;
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

    /// A `Handle<Surface>` for `index` (same throwaway-store trick).
    fn surface_handle_at(index: u32) -> Handle<Surface> {
        let plane = || {
            Surface::Plane(
                Plane::through_points(
                    Point3::from_array([0.0, 0.0, 0.0]),
                    Point3::from_array([1.0, 0.0, 0.0]),
                    Point3::from_array([0.0, 1.0, 0.0]),
                )
                .unwrap(),
            )
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
        m.solids.push(Solid {
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
            VertexDef::ThreePlane([
                surface_handle_at(0),
                surface_handle_at(1),
                surface_handle_at(9), // out of bounds — only 6 surfaces
            ]),
            Point3::origin(),
            Some(1e-9),
        );
        assert_eq!(
            validate(&m),
            vec![Violation::DanglingReference {
                kind: RefKind::VertexDefSurface,
                owner_index: vh.index(),
                target_index: 9,
                target_len: 6,
            }]
        );
    }

    #[test]
    fn vertex_def_is_copy() {
        let d = VertexDef::ThreePlane([surface_handle_at(0); 3]);
        let copy = d; // move-or-copy
        let _again = d; // still usable ⇒ Copy, not moved
        assert_eq!(d, copy);
    }

    #[test]
    fn stray_vertex_ignored() {
        // A vertex referenced by nothing is a dead arena item, not a defect: it
        // is unreachable from the live cube, so validate ignores it (design §2).
        // (Under the old whole-store count this raised EulerParity{v:9}.)
        let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        m.push_vertex(
            VertexDef::ThreePlane([
                m.world_plane(nacre_scalar::Axis::Z),
                m.world_plane(nacre_scalar::Axis::X),
                m.world_plane(nacre_scalar::Axis::Y),
            ]),
            Point3::origin(),
            None,
        );
        assert!(validate(&m).is_empty());
    }

    #[test]
    fn supersede_solid_reuses_cells() {
        // What an editing op does (design §2): build a new solid that *reuses*
        // the old solid's cells, then drop the old solid from `live_solids`. The
        // old cells linger in the arena but, unreferenced by any live solid, fall
        // out of the reachable closure — so each shared edge is counted twice
        // (once per live face), not four times, and validate stays clean.
        let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let old = m.live_solids[0];
        let old_shell = m.solids.get(old).outer;
        let faces = m.shells.get(old_shell).faces.clone();
        let new_shell = m.shells.push(Shell { faces });
        let new_solid = m.push_solid(Solid {
            outer: new_shell,
            cavities: vec![],
        });
        m.live_solids.retain(|&s| s != old); // supersede: only the new one is live
        assert_eq!(m.live_solids, vec![new_solid]);
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
            let (_, f0) = m.faces.iter().next().unwrap();
            Face {
                surface: f0.surface,
                outer: f0.outer.clone(),
                inner: vec![],
                orientation: f0.orientation,
            }
        };
        m.faces.push(dup);
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
        /// (vertex index, coordinate delta, measured tol) — moves a vertex off its
        /// (un-moved) surfaces; `Some(tol)` marks it measured (the old `Discovered`),
        /// `None` constructed.
        nudge: Option<(usize, [f64; 3], Option<f64>)>,
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
                let mut tol = None;
                if let Some((vi, d, t)) = opts.nudge {
                    if vi == i {
                        p = [p[0] + d[0], p[1] + d[1], p[2] + d[2]];
                        tol = t;
                    }
                }
                // Every vertex names the three tetra faces incident to it (S7: the
                // definition is the vertex).
                let incident: Vec<Handle<Surface>> = TETRA_FACES
                    .iter()
                    .enumerate()
                    .filter(|(_, (tri, _))| tri.contains(&i))
                    .map(|(fi, _)| sh[fi])
                    .collect();
                m.push_vertex(
                    VertexDef::ThreePlane(
                        incident.try_into().expect("a tetra vertex is on 3 faces"),
                    ),
                    Point3::from_array(p),
                    tol,
                )
            })
            .collect();
        let eh: Vec<Handle<Edge>> = (0..6)
            .map(|i| {
                let (a, b) = TETRA_EDGES[i];
                // The two faces whose loops use edge `i` — its carriers, read off the same
                // table the loops are built from. `push_edge` derives the curve through the
                // (possibly nudged) vertex points — an endpoint can no longer sit off its own
                // line, which is the S8 point (see `vertex_off_its_rim_circle`).
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
            faces.push(m.faces.push(Face {
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
        let sh = m.shells.push(Shell { faces });
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
        let sh = m.shells.push(Shell { faces });
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

    #[test]
    fn vertex_off_surface_when_nudged() {
        // Vertex 0 moved by 2·EPS_CONSTRUCTED; the planes stay on the un-moved corners, so
        // the vertex is off them. (`VertexOffCurve` cannot fire here since S8: a line edge's
        // curve derives *through its endpoints*, so an endpoint is on its own line by
        // construction — the check's remaining teeth are the circles, see
        // `vertex_off_its_rim_circle`.)
        let vs = validate(&tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 2.0 * EPS_CONSTRUCTED], None)),
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

    /// ★ S8 negative control: an edge whose stated carriers disagree with the two faces
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
                VertexDef::ThreePlane([sa, sb, sc]),
                Point3::from_array(p),
                None,
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
            m.faces.push(Face {
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
        let shell = m.shells.push(Shell {
            faces: vec![fa, fb],
        });
        m.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        let vs = validate(&m);
        assert!(
            vs.iter().any(
                |x| matches!(x, Violation::EdgeCarrierMismatch { edge, .. } if *edge == e_shared)
            ),
            "stated [sa, sc] vs observed [sa, sb] must be flagged: {vs:?}"
        );
        // The other shared-plane edges (1 use each) are non-manifold, not carrier mismatches.
        assert!(
            vs.iter()
                .any(|x| matches!(x, Violation::NonManifoldEdge { .. }))
        );
    }

    /// ★ S7 negative controls: a definition variant whose carrier kinds contradict it is
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
        let cyl = m
            .faces
            .iter()
            .map(|(_, f)| f.surface)
            .find(|&h| matches!(m.surface(h), Surface::Cylinder(_)))
            .expect("the lateral cylinder");
        let (z0, x0) = (
            m.world_plane(nacre_scalar::Axis::Z),
            m.world_plane(nacre_scalar::Axis::X),
        );
        let bad_three = m.push_vertex(VertexDef::ThreePlane([z0, x0, cyl]), Point3::origin(), None);
        let bad_seam = m.push_vertex(VertexDef::OnSeam([z0, x0]), Point3::origin(), None);
        let vs = validate(&m);
        for bad in [bad_three, bad_seam] {
            assert!(
                vs.iter().any(
                    |v| matches!(v, Violation::VertexDefCarrierMismatch { vertex } if *vertex == bad)
                ),
                "a contradictory def must be flagged: {vs:?}"
            );
        }
    }

    /// ★ S8: the check `VertexOffCurve` still has teeth where the curve does NOT derive from
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
        let lateral = m
            .faces
            .iter()
            .map(|(_, f)| f.surface)
            .find(|&h| matches!(m.surface(h), Surface::Cylinder(_)))
            .expect("the cylinder's lateral surface");
        let cap = m
            .faces
            .iter()
            .map(|(_, f)| f.surface)
            .find(|&h| h != lateral && !matches!(m.surface(h), Surface::Cylinder(_)))
            .expect("a cap plane");
        // A seam vertex at radius 2 + 1e-3 — off the derived radius-2 rim circle. Its
        // definition is `OnSeam` (the lateral cylinder and its cap), which is what a rim's
        // endpoint always is; the coordinate is the part that lies.
        let bad = m.push_vertex(
            VertexDef::OnSeam([lateral, cap]),
            Point3::from_array([2.001, 0.0, 0.0]),
            None,
        );
        let rim = m
            .push_edge([lateral, cap], [bad, bad])
            .expect("a rim derives from its carriers, not its vertices");
        let face = m.faces.push(Face {
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
        let shell = m.shells.push(Shell { faces: vec![face] });
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

    #[test]
    fn discovered_vertex_within_tolerance_is_clean() {
        let tol = 1e-6;
        let m = tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 0.5 * tol], Some(tol))),
            ..Default::default()
        });
        assert!(validate(&m).is_empty());
    }

    #[test]
    fn discovered_vertex_outside_tolerance_flags() {
        let tol = 1e-6;
        let vs = validate(&tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 2.0 * tol], Some(tol))),
            ..Default::default()
        }));
        assert!(vs.iter().any(|v| matches!(
            v,
            Violation::VertexOffSurface { .. } | Violation::VertexOffCurve { .. }
        )));
    }

    #[test]
    fn discovered_vertex_off_its_definition_flags() {
        // A discovered vertex nudged beyond tol is off its three definition
        // planes (its incident faces). Assert the definition check specifically
        // fires — not merely "some violation" (which VertexOffSurface satisfies).
        let tol = 1e-6;
        let vs = validate(&tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 2.0 * tol], Some(tol))),
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
        let b_outer = m.solids.get(b).outer;
        let void = if reverse {
            m.reversed_shell(b_outer)
        } else {
            b_outer
        };
        let a_outer = m.solids.get(a).outer;
        let hollow = m.push_solid(Solid {
            outer: a_outer,
            cavities: vec![void],
        });
        m.live_solids.retain(|&s| s == hollow);
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
}
