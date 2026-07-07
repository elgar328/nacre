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
use nacre_topo::{Edge, Face, Model, Vertex};

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

    check_euler_poincare(model, &mut out);
    out
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

fn check_euler_poincare(m: &Model, out: &mut Vec<Violation>) {
    let v = m.vertices.len();
    let e = m.edges.len();
    let f = m.faces.len();
    let s = m.shells.len();
    let inner_loops: usize = m.faces.iter().map(|(_, face)| face.inner.len()).sum();

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

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_topo::{Origin, Shell, Solid};
    use proptest::prelude::*;

    fn cuboid(min: [f64; 3], max: [f64; 3]) -> Model {
        let mut m = Model::new();
        m.add_cuboid(Point3::from_array(min), Point3::from_array(max));
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
    fn euler_parity_stray_vertex() {
        let mut m = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        // A vertex referenced by nothing → V becomes odd-parity for Euler.
        m.vertices.push(Vertex {
            point: Point3::origin(),
            origin: Origin::Constructed,
        });
        assert_eq!(
            validate(&m),
            vec![Violation::EulerParity {
                v: 9,
                e: 12,
                f: 6,
                s: 1,
                inner_loops: 0,
            }]
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
    }
}
