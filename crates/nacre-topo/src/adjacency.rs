//! `Adjacency` — the reverse index over the topology stores.

use crate::Model;
use crate::topology::{Edge, Face, Vertex};
use nacre_store::Handle;
use std::collections::HashMap;

/// Reverse index over the topology stores — a **cache, not truth** (design §4,
/// overview 절대원칙 5). Rebuilt from scratch by scanning the stores; discard
/// and regenerate any time. Validation ("every edge used by exactly two faces
/// in opposite directions") and adjacency queries read this.
///
/// Uses plain `Vec` rather than the design's `SmallVec<[_; 2]>` / `[_; 4]` for
/// now — one fewer dependency; the inline-storage win is a future perf tweak.
#[derive(Clone, Debug, Default)]
pub struct Adjacency {
    /// Edge → every `(face, half-edge.forward)` that uses it. A manifold edge
    /// has exactly two entries with opposite `bool`s.
    pub edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>,
    /// Vertex → the edges bounded by it (either endpoint).
    pub vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
}

impl Adjacency {
    /// Build a fresh index from a model's topology stores. Scans every face's
    /// outer + inner loops → half-edges → edges (recording `forward`), and every
    /// bounded edge → its two endpoint vertices. From scratch — cheap, always correct.
    ///
    /// ★ **Rebuild is the whole story.** An earlier note promised incremental maintenance
    /// "with the operation log (M2)"; the log arrived and the incremental path did not, and
    /// nothing has needed it: [`Model::rebuild_adjacency`] replaces the index wholesale and
    /// callers run it once after a batch. Build it when a profile asks for it, not before.
    pub fn rebuild(model: &Model) -> Adjacency {
        // Only the live model is indexed (design §2): superseded cells left in
        // the arena must not pollute edge use-counts. Iterate the stores in
        // order (deterministic) filtered by the reachable closure.
        let reach = model.reachable();
        let mut edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>> = HashMap::new();
        let mut i = 0u32;
        while let Some(fh) = model.face_handle_at(i) {
            i += 1;
            let face = model.face(fh);
            if !reach.faces.contains(&fh) {
                continue;
            }
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    edge_uses.entry(he.edge).or_default().push((fh, he.forward));
                }
            }
        }
        let mut vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>> = HashMap::new();
        let mut i = 0u32;
        while let Some(eh) = model.edge_handle_at(i) {
            i += 1;
            let edge = model.edge(eh);
            if !reach.edges.contains(&eh) {
                continue;
            }
            let [a, b] = edge.vertices;
            vertex_edges.entry(a).or_default().push(eh);
            vertex_edges.entry(b).or_default().push(eh);
        }
        Adjacency {
            edge_uses,
            vertex_edges,
        }
    }
}

/// Union-find root of `x` (path compression); roots are the smaller index (deterministic).
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut c = x;
    while parent[c] != r {
        let next = parent[c];
        parent[c] = r;
        c = next;
    }
    r
}

fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra.max(rb)] = ra.min(rb);
    }
}

/// The vertices whose surface **link is not a single circle** — non-manifold "pinch" points where
/// two or more face-fans meet at one vertex (two cubes touching only at a corner; a cut whose tool
/// corner lands on a target's concave corner). Pure combinatorics on the reverse index, no
/// coordinates:
///
/// At a vertex `V`, each incident face contributes its two `V`-edges as an arc of the link (a face
/// using edge `e` at `V` has `V` as an endpoint of `e`, and `e` is one of the two loop edges the
/// face turns through at `V`). Union those arcs over `V`'s edges and count connected components. A
/// manifold vertex's link is one cycle → **one component**; a pinch is **≥2 components** (or a face
/// using ≠2 edges at `V` — a pinched face passing through `V` twice).
///
/// **Soundness** rests on the caller's edges being manifold (every edge used by exactly two faces —
/// the boolean's closed-shell guard, and what a closed solid always satisfies): the link is then
/// 2-regular, so its components are exactly its cycles (fans).
///
/// **Scope — planar (polyhedral) topology.** A polyhedral corner uses exactly two of its vertex's
/// edges per face; the function **abstains** (treats the vertex as manifold) at any vertex where a
/// face uses a different number — a cylinder's circular cap (one closed edge), its seam edge (used
/// twice by the lateral face), or a self-loop rim edge (`[v, v]`). These are parametric artifacts
/// this link analysis cannot classify (M6 quadrics), never a *planar* pinch, so real planar pinches
/// are unaffected. (A genuinely pinched *face* — a loop through `V` twice — is thus not flagged
/// here; the M5 boolean never emits one, and the edge check covers non-manifold edges.)
pub fn nonmanifold_vertices(
    vertex_edges: &HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
    edge_uses: &HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>,
) -> Vec<Handle<Vertex>> {
    let mut out = Vec::new();
    for (&v, raw_edges) in vertex_edges {
        // Distinct edges at `v` (a self-loop `[v, v]` is listed twice by the caller).
        let mut seen = std::collections::HashSet::new();
        let edges: Vec<Handle<Edge>> = raw_edges
            .iter()
            .copied()
            .filter(|e| seen.insert(*e))
            .collect();
        let idx: HashMap<Handle<Edge>, usize> =
            edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
        // face → the local indices of the `V`-edges it uses (its corner at `V`).
        let mut face_edges: HashMap<Handle<Face>, Vec<usize>> = HashMap::new();
        for &e in &edges {
            for &(f, _) in edge_uses.get(&e).map(Vec::as_slice).unwrap_or(&[]) {
                face_edges.entry(f).or_default().push(idx[&e]);
            }
        }
        // Abstain on non-polyhedral corners (see scope note).
        if face_edges.values().any(|es| es.len() != 2) {
            continue;
        }
        let mut parent: Vec<usize> = (0..edges.len()).collect();
        for es in face_edges.values() {
            uf_union(&mut parent, es[0], es[1]);
        }
        let components: std::collections::HashSet<usize> =
            (0..edges.len()).map(|i| uf_find(&mut parent, i)).collect();
        if components.len() > 1 {
            out.push(v);
        }
    }
    out.sort_unstable_by_key(|h| h.index());
    out
}

#[cfg(test)]
mod tests {
    use super::nonmanifold_vertices;
    use crate::topology::{Edge, Face, HalfEdge, Loop};
    use crate::{Handle, Model, Orientation, PointCache, Shell, Solid, Vertex};
    use nacre_geom::Plane;
    use nacre_math::Point3;
    use std::collections::HashMap;

    /// The pinch detector is per-vertex: a single face-fan at a vertex is manifold, two+ fans is a
    /// pinch — and **every** pinch is flagged independently, so an *even* number of them is caught
    /// (which an Euler-parity count cannot: two pinches flip χ back to even). Builds the reverse-index
    /// maps directly with minted handles (the detector reads only the maps, never the stores).
    #[test]
    fn nonmanifold_vertices_flags_each_pinch_even_count() {
        let mut m = Model::new();
        let r = nacre_scalar::Rat::from_int;
        let surface = m.push_plane_unregistered(
            Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        );
        // Hand-built cells: the def names the three world seeds — real, distinct planes the
        // detector never dereferences (it reads only the maps).
        let mk_v = |m: &mut Model| {
            let def = crate::VertexDef::ThreePlane([
                m.world_plane(nacre_scalar::Axis::Z),
                m.world_plane(nacre_scalar::Axis::X),
                m.world_plane(nacre_scalar::Axis::Y),
            ]);
            m.push_vertex(
                def,
                PointCache::Unrealized {
                    coord: Point3::origin(),
                },
            )
        };
        let mut vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>> = HashMap::new();
        let mut edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>> = HashMap::new();
        // Add one triangular face-fan (3 edges from `apex`, 3 faces cyclically joining them).
        let add_fan =
            |m: &mut Model,
             apex: Handle<Vertex>,
             vertex_edges: &mut HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
             edge_uses: &mut HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>| {
                let e: Vec<Handle<Edge>> = (0..3)
                    .map(|_| {
                        let other = mk_v(m);
                        // The detector reads only the maps — a self-pair of the one dummy
                        // surface is enough (never dereferenced as adjacency here; a raw
                        // `edges.push` skips the curve cache, which nothing here reads).
                        m.edges.push(Edge {
                            surfaces: [surface, surface],
                            vertices: [apex, other],
                        })
                    })
                    .collect();
                vertex_edges
                    .entry(apex)
                    .or_default()
                    .extend(e.iter().copied());
                for i in 0..3 {
                    let f = m.faces.push(Face {
                        surface,
                        outer: Loop { half_edges: vec![] },
                        inner: vec![],
                        orientation: Orientation::Forward,
                    });
                    edge_uses.entry(e[i]).or_default().push((f, true));
                    edge_uses.entry(e[(i + 1) % 3]).or_default().push((f, true));
                }
            };
        // Manifold apex: one fan. Two pinched apexes: two fans each.
        let ok = mk_v(&mut m);
        add_fan(&mut m, ok, &mut vertex_edges, &mut edge_uses);
        let pinch1 = mk_v(&mut m);
        add_fan(&mut m, pinch1, &mut vertex_edges, &mut edge_uses);
        add_fan(&mut m, pinch1, &mut vertex_edges, &mut edge_uses);
        let pinch2 = mk_v(&mut m);
        add_fan(&mut m, pinch2, &mut vertex_edges, &mut edge_uses);
        add_fan(&mut m, pinch2, &mut vertex_edges, &mut edge_uses);

        let mut got = nonmanifold_vertices(&vertex_edges, &edge_uses);
        got.sort_unstable_by_key(|h| h.index());
        let mut want = vec![pinch1, pinch2];
        want.sort_unstable_by_key(|h| h.index());
        assert_eq!(got, want, "both pinches flagged, the manifold fan is not");
    }

    /// Two triangles sharing edge v0–v1, wound so the shared edge is used once
    /// forward and once reversed. Exercises `rebuild` independent of the cube.
    #[test]
    fn rebuild_indexes_a_shared_edge() {
        let mut m = Model::new();
        let mk_v = |m: &mut Model, p: [f64; 3]| {
            let def = crate::VertexDef::ThreePlane([
                m.world_plane(nacre_scalar::Axis::Z),
                m.world_plane(nacre_scalar::Axis::X),
                m.world_plane(nacre_scalar::Axis::Y),
            ]);
            m.push_vertex(
                def,
                PointCache::Unrealized {
                    coord: Point3::from_array(p),
                },
            )
        };
        let v0 = mk_v(&mut m, [0.0, 0.0, 0.0]);
        let v1 = mk_v(&mut m, [1.0, 0.0, 0.0]);
        let v2 = mk_v(&mut m, [0.0, 1.0, 0.0]);
        let v3 = mk_v(&mut m, [0.0, -1.0, 0.0]);

        // rebuild never dereferences surface handles, so one dummy is enough for
        // every edge/face (raw `edges.push` — nothing here reads the curve cache).
        let r = nacre_scalar::Rat::from_int;
        let surface = m.push_plane_unregistered(
            Plane::from_point_normal(
                Point3::origin(),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            [[r(0); 3], [r(1), r(0), r(0)], [r(0), r(1), r(0)]],
        );
        let mk_e = |m: &mut Model, a, b| {
            m.edges.push(Edge {
                // rebuild never dereferences carriers either — the dummy pair mirrors the
                // dummy surface above.
                surfaces: [surface, surface],
                vertices: [a, b],
            })
        };
        let e_shared = mk_e(&mut m, v0, v1);
        let ea1 = mk_e(&mut m, v1, v2);
        let ea2 = mk_e(&mut m, v2, v0);
        let eb1 = mk_e(&mut m, v1, v3);
        let eb2 = mk_e(&mut m, v3, v0);

        let mk_face = |m: &mut Model, hes: Vec<HalfEdge>| {
            m.faces.push(Face {
                surface,
                outer: Loop { half_edges: hes },
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Triangle A: v0→v1→v2→v0, shared edge forward.
        let fa = mk_face(
            &mut m,
            vec![
                HalfEdge {
                    edge: e_shared,
                    forward: true,
                },
                HalfEdge {
                    edge: ea1,
                    forward: true,
                },
                HalfEdge {
                    edge: ea2,
                    forward: true,
                },
            ],
        );
        // Triangle B: v1→v0→v3→v1, shared edge reversed.
        let fb = mk_face(
            &mut m,
            vec![
                HalfEdge {
                    edge: e_shared,
                    forward: false,
                },
                HalfEdge {
                    edge: eb2,
                    forward: false,
                },
                HalfEdge {
                    edge: eb1,
                    forward: false,
                },
            ],
        );

        // Wrap the two faces in a live solid — adjacency now indexes only the
        // reachable closure (design §2), so loose faces would be invisible. This
        // pair is an open surface (validate would flag it), but the test only
        // inspects the rebuilt adjacency, which needs the faces reachable.
        let shell = m.shells.push(Shell {
            faces: vec![fa, fb],
        });
        m.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        m.rebuild_adjacency();

        // shared edge: two uses, opposite flags.
        let shared = &m.adj.edge_uses[&e_shared];
        assert_eq!(shared.len(), 2);
        assert_ne!(shared[0].1, shared[1].1);
        // a private edge: one use.
        assert_eq!(m.adj.edge_uses[&ea1].len(), 1);
        // v0 is incident to e_shared, ea2, eb2.
        assert_eq!(m.adj.vertex_edges[&v0].len(), 3);
        assert_eq!(m.adj.vertex_edges[&v2].len(), 2);
    }
}
