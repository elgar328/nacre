//! `Adjacency` — the reverse index over the topology stores.

use crate::Model;
use crate::Vertex;
use crate::topology::{Edge, Face};
use nacre_store::Handle;
use std::collections::HashMap;

/// Reverse index over the topology stores — derived, not truth. Built from scratch by scanning the
/// stores, by the consumer that needs it: validation, for "every edge used by exactly two faces in
/// opposite directions". The boolean's result check and tess keep their own edge tables, and the
/// model holds none.
///
/// Uses plain `Vec` rather than the design's `SmallVec<[_; 2]>` / `[_; 4]`: one fewer dependency,
/// and the inline-storage win is left for a profile to ask for.
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
    /// **Rebuild is the whole story.** There is no incremental path, and nothing needs one: each
    /// consumer builds the index wholesale when it asks.
    pub fn rebuild(model: &Model) -> Adjacency {
        // Only the live model is indexed: superseded cells left in
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

// How it reads: the link's nodes are edge ends, its arcs are the loops' corners. An end is an edge
// and one of its two endpoint slots (`vertices[0]` or `vertices[1]`); a corner is two consecutive
// half-edges of a loop, and joins the end the first arrives by to the end the second leaves by. At
// each vertex the ends meeting there fall into connected components: one is a manifold vertex, two
// or more a pinch. A closed edge `[v, v]` — a whole rim — is two ends at `v`, so a disk cap and the
// lateral beside it make one component through the rim's two ends, and a fan of other faces meeting
// the rim's vertex is a second. Two premises make the components the fans: every edge has two uses,
// so the link is 2-regular; and two consecutive half-edges of a loop bound one angular sector of
// their face — a face meeting a vertex with two separate loops loses which ends pair into sectors,
// and no reading of loops alone gets them back. A corner whose two half-edges do not meet at one
// vertex (an open loop) joins nothing.
/// The vertices of `faces` where two or more fans of faces meet — non-manifold "pinch" points (two
/// cubes touching only at a corner; a cut whose tool corner lands on a target's concave corner; a
/// box whose corner rests on a cap's rim beside it), in handle order. Read off the faces' loops; no
/// coordinates.
///
/// The faces must use every edge exactly twice, and a face must pass a vertex twice in **one**
/// loop, a ring through the point twice. Where a face meets one vertex with two separate loops the
/// answer there is not to be trusted either way: a pinch can be missed and a 2-manifold point
/// reported.
pub fn nonmanifold_vertices(
    model: &Model,
    faces: impl IntoIterator<Item = Handle<Face>>,
) -> Vec<Handle<Vertex>> {
    let mut index: HashMap<(Handle<Edge>, usize), usize> = HashMap::new();
    let mut at: Vec<Handle<Vertex>> = Vec::new();
    let mut parent: Vec<usize> = Vec::new();
    let mut node = |end: (Handle<Edge>, usize), v: Handle<Vertex>, parent: &mut Vec<usize>| {
        *index.entry(end).or_insert_with(|| {
            at.push(v);
            parent.push(parent.len());
            parent.len() - 1
        })
    };
    for fh in faces {
        let face = model.face(fh);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            let k = lp.half_edges.len();
            for i in 0..k {
                let (a, b) = (lp.half_edges[i], lp.half_edges[(i + 1) % k]);
                // A forward half-edge runs `vertices[0] → vertices[1]`: it arrives by slot 1 and
                // leaves by slot 0.
                let (sa, sb) = (usize::from(a.forward), usize::from(!b.forward));
                let (va, vb) = (
                    model.edge(a.edge).vertices[sa],
                    model.edge(b.edge).vertices[sb],
                );
                if va != vb {
                    continue;
                }
                let (na, nb) = (
                    node((a.edge, sa), va, &mut parent),
                    node((b.edge, sb), vb, &mut parent),
                );
                uf_union(&mut parent, na, nb);
            }
        }
    }
    let mut roots: HashMap<Handle<Vertex>, std::collections::HashSet<usize>> = HashMap::new();
    for (i, &v) in at.iter().enumerate() {
        let r = uf_find(&mut parent, i);
        roots.entry(v).or_default().insert(r);
    }
    let mut out: Vec<Handle<Vertex>> = roots
        .into_iter()
        .filter(|(_, r)| r.len() > 1)
        .map(|(v, _)| v)
        .collect();
    out.sort_unstable_by_key(|h| h.index());
    out
}

#[cfg(test)]
#[path = "tests/adjacency.rs"]
mod tests;
