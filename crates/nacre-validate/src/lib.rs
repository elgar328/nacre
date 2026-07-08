//! b-rep invariant checker for the nacre kernel (design.md §7).
//!
//! [`validate`] runs every M1 check over a [`Model`] and returns all
//! [`Violation`]s it finds (an empty `Vec` means the model is valid). Checks:
//! reference integrity, loop closure, half-edge manifold pairing, geometric
//! incidence (vertices on their curves/surfaces), and Euler-Poincaré. The
//! tessellation checks (§5, §7 — provenance coherence, crack-free) arrive in M3
//! when a `Tessellation` exists.

use nacre_math::Point3;
use nacre_store::{Handle, Store};
use nacre_topo::{Adjacency, Edge, Face, Loop, Model, Origin, Reachable, Vertex};

/// Residual bound for a `Constructed` vertex lying on its reference
/// curve/surface. Machine epsilon (~2.2e-16) is too tight — a `Constructed`
/// coordinate is the output of a short floating-point construction chain
/// (corner arithmetic, line/plane fitting), so its residual against its own
/// fitted geometry is ~magnitude·(a few ULP) ≈ 1e-13 for coordinates up to
/// ~1e3. `1e-9` sits comfortably above that floor yet ~7 orders below any
/// `Discovered` tolerance.
pub const EPS_CONSTRUCTED: f64 = 1e-9;

/// Which reference edge in the topology graph a dangling handle sits on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefKind {
    EdgeCurve,
    EdgeBoundVertex,
    FaceSurface,
    HalfEdgeEdge,
    ShellFace,
    SolidShell,
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

    /// A half-edge in a loop resolves to an edge with `bounds == None`. M1
    /// expects every edge bounded; a closed edge (M3) here can't be checked for
    /// closure and is reported (its continuity check is skipped, not assumed).
    UnboundedEdgeInLoop {
        face: Handle<Face>,
        loop_kind: LoopKind,
        edge: Handle<Edge>,
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
    check_loop_closure(model, &reach, &mut out);
    check_manifold(model, &adj, &reach, &mut out);
    check_geometric_incidence(model, &reach, &mut out);
    check_euler_poincare(model, &reach, &mut out);
    out
}

/// The tolerance an element's provenance grants: `Constructed` is exact to
/// [`EPS_CONSTRUCTED`]; `Discovered` carries its own measured tolerance.
#[inline]
fn tol_of(o: Origin) -> f64 {
    match o {
        Origin::Constructed => EPS_CONSTRUCTED,
        Origin::Discovered { tol } => tol,
    }
}

#[inline]
fn in_bounds<T>(h: Handle<T>, store: &Store<T>) -> bool {
    (h.index() as usize) < store.len()
}

fn check_reference_integrity(m: &Model, out: &mut Vec<Violation>) {
    for (eh, edge) in m.edges.iter() {
        if !in_bounds(edge.curve, &m.curves) {
            out.push(Violation::DanglingReference {
                kind: RefKind::EdgeCurve,
                owner_index: eh.index(),
                target_index: edge.curve.index(),
                target_len: m.curves.len() as u32,
            });
        }
        if let Some(bounds) = edge.bounds {
            for v in bounds {
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
        if !in_bounds(face.surface, &m.surfaces) {
            out.push(Violation::DanglingReference {
                kind: RefKind::FaceSurface,
                owner_index: fh.index(),
                target_index: face.surface.index(),
                target_len: m.surfaces.len() as u32,
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
    // Resolve each half-edge to (start, end); unbounded edges are reported and
    // skipped in the continuity walk (not assumed closed).
    let ends: Vec<Option<(Handle<Vertex>, Handle<Vertex>)>> = hes
        .iter()
        .map(|he| match m.edges.get(he.edge).bounds {
            None => {
                out.push(Violation::UnboundedEdgeInLoop {
                    face: fh,
                    loop_kind: kind,
                    edge: he.edge,
                });
                None
            }
            Some([a, b]) => Some(if he.forward { (a, b) } else { (b, a) }),
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
        } else if uses[0].1 == uses[1].1 {
            out.push(Violation::NonOpposedEdge {
                edge: eh,
                faces: [uses[0].0, uses[1].0],
            });
        }
    }
}

fn check_geometric_incidence(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Each live edge's bound vertices must lie on the edge's curve.
    for (eh, edge) in m.edges.iter() {
        if !reach.edges.contains(&eh) {
            continue;
        }
        if let Some([a, b]) = edge.bounds {
            let curve = m.curves.get(edge.curve);
            for vh in [a, b] {
                let vertex = m.vertices.get(vh);
                let residual = curve.distance(vertex.point);
                let tol = tol_of(vertex.origin).max(tol_of(edge.origin));
                if residual > tol {
                    out.push(Violation::VertexOffCurve {
                        edge: eh,
                        vertex: vh,
                        point: vertex.point,
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
        let surface = m.surfaces.get(face.surface);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                if let Some([a, b]) = m.edges.get(he.edge).bounds {
                    let vh = if he.forward { a } else { b };
                    let vertex = m.vertices.get(vh);
                    let residual = surface.distance(vertex.point);
                    let tol = tol_of(vertex.origin);
                    if residual > tol {
                        out.push(Violation::VertexOffSurface {
                            face: fh,
                            vertex: vh,
                            point: vertex.point,
                            residual,
                            tol,
                        });
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_geom::{Curve, Line, Plane, Surface};
    use nacre_math::Vector3;
    use nacre_topo::{HalfEdge, Orientation, Origin, Shell, Solid};
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
    fn stray_vertex_ignored() {
        // A vertex referenced by nothing is a dead arena item, not a defect: it
        // is unreachable from the live cube, so validate ignores it (design §2).
        // (Under the old whole-store count this raised EulerParity{v:9}.)
        let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        m.vertices.push(Vertex {
            point: Point3::origin(),
            origin: Origin::Constructed,
        });
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
        /// (vertex index, coordinate delta, that vertex's origin) — moves a
        /// vertex off its (un-moved) curves/surfaces.
        nudge: Option<(usize, [f64; 3], Origin)>,
        drop_face: Option<usize>,
        unbind_edge: Option<usize>,
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
        let vh: Vec<Handle<Vertex>> = (0..4)
            .map(|i| {
                let mut p = corner(i).as_array();
                let mut origin = Origin::Constructed;
                if let Some((vi, d, o)) = opts.nudge {
                    if vi == i {
                        p = [p[0] + d[0], p[1] + d[1], p[2] + d[2]];
                        origin = o;
                    }
                }
                m.vertices.push(Vertex {
                    point: Point3::from_array(p),
                    origin,
                })
            })
            .collect();
        let eh: Vec<Handle<Edge>> = (0..6)
            .map(|i| {
                let (a, b) = TETRA_EDGES[i];
                let curve = m.curves.push(Curve::Line(
                    Line::through_points(corner(a), corner(b)).unwrap(),
                ));
                let bounds = if opts.unbind_edge == Some(i) {
                    None
                } else {
                    Some([vh[a], vh[b]])
                };
                m.edges.push(Edge {
                    curve,
                    bounds,
                    origin: Origin::Constructed,
                })
            })
            .collect();
        let mut faces = Vec::new();
        for (fi, (tri, hes)) in TETRA_FACES.iter().enumerate() {
            if opts.drop_face == Some(fi) {
                continue;
            }
            let surface = m.surfaces.push(Surface::Plane(
                Plane::through_points(corner(tri[0]), corner(tri[1]), corner(tri[2])).unwrap(),
            ));
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
    fn unbounded_edge_in_loop() {
        let vs = validate(&tetra_with(TetraOpts {
            unbind_edge: Some(0),
            ..Default::default()
        }));
        assert!(
            vs.iter()
                .any(|v| matches!(v, Violation::UnboundedEdgeInLoop { .. }))
        );
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
    fn vertex_off_curve_and_surface_when_nudged() {
        // Vertex 0 moved by 2·EPS_CONSTRUCTED; curves/planes stay on the
        // un-moved corners, so the vertex is off some of them.
        let vs = validate(&tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 2.0 * EPS_CONSTRUCTED], Origin::Constructed)),
            ..Default::default()
        }));
        assert!(
            vs.iter()
                .any(|v| matches!(v, Violation::VertexOffSurface { .. }))
        );
        assert!(
            vs.iter()
                .any(|v| matches!(v, Violation::VertexOffCurve { .. }))
        );
    }

    #[test]
    fn discovered_vertex_within_tolerance_is_clean() {
        let tol = 1e-6;
        let m = tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 0.5 * tol], Origin::Discovered { tol })),
            ..Default::default()
        });
        assert!(validate(&m).is_empty());
    }

    #[test]
    fn discovered_vertex_outside_tolerance_flags() {
        let tol = 1e-6;
        let vs = validate(&tetra_with(TetraOpts {
            nudge: Some((0, [0.0, 0.0, 2.0 * tol], Origin::Discovered { tol })),
            ..Default::default()
        }));
        assert!(vs.iter().any(|v| matches!(
            v,
            Violation::VertexOffSurface { .. } | Violation::VertexOffCurve { .. }
        )));
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
