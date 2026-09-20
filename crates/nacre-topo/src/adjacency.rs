//! `Adjacency` — the reverse index over the topology stores.

use crate::Model;
use crate::Vertex;
use crate::topology::{Edge, Face};
use nacre_store::Handle;
use std::collections::HashMap;

/// Reverse index over the topology stores — a **cache, not truth**
/// (overview 절대원칙 5). Rebuilt from scratch by scanning the stores; discard
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
#[path = "tests/adjacency.rs"]
mod tests;
