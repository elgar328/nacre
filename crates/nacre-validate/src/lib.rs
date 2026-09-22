//! b-rep invariant checker for the nacre kernel.
//!
//! [`validate`] runs every M1 check over a [`Model`] and returns all
//! [`Violation`]s it finds (an empty `Vec` means the model is valid). Checks:
//! reference integrity, loop closure, half-edge manifold pairing, face
//! orientation against loop winding, geometric incidence (vertices on their
//! curves/surfaces), and Euler-Poincaré. The
//! tessellation checks (provenance coherence, crack-free) arrive in M3
//! when a `Tessellation` exists.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Adjacency, Edge, Face, Loop, Model, Reachable, Shell, Solid, Surface, Vertex};

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
    VertexSurface,
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

    /// A vertex definition whose carrier kinds contradict its variant: a `ThreePlane`
    /// naming a cylinder (three *planes* is the claim), or an `OnSeam` naming two planes
    /// (two planes meet in a line — a line's point IS a three-plane intersection, so the
    /// seam spelling would be hiding an expressible truth). The variants' invariants are
    /// per-variant, and this is the checker that keeps them so — the vertex sibling of
    /// [`Self::EdgeCarrierMismatch`].
    VertexCarrierMismatch { vertex: Handle<Vertex> },

    /// An edge's stated carriers disagree with adjacency: the multiset of the two face
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

    /// A vertex does not lie on one of its definition's surfaces within tolerance — the
    /// definition is the truth, so the cached point must sit within `tol` of every surface it
    /// is defined as meeting; `tol` is the measured one where there is one, else
    /// [`EPS_CONSTRUCTED`]. Checked for **every** vertex. `surface_index` is
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
    /// subtracting it (M5 containment). The other checks miss it — a
    /// globally-reversed shell keeps edge opposition (so `check_manifold`
    /// passes) and the same V/E/F/S (so Euler passes). Planar cavities only; a
    /// non-planar void is not checked here.
    CavityMisoriented {
        solid: Handle<Solid>,
        cavity: Handle<Shell>,
        signed_volume: f64,
    },

    /// A planar face's **loop winding disagrees with its stated orientation**: the
    /// loop's own witness (a polygon's Newell area vector, or a closed rim's
    /// circle) and the outward normal the face states (`plane.normal()` ×
    /// `orientation`, the same normal [`shell_signed_volume`] integrates) do not
    /// stand as they must. `cos` is their unit dot, and it says which loop spoke:
    ///
    /// * an **outer** loop must *agree* — healthy is `cos ≈ +1`, so `≈ −1` is a
    ///   flipped orientation flag and `≈ 0` a loop that does not span its plane;
    /// * an **inner** loop must *oppose* (a hole winds the other way round, the
    ///   rule `build_prism` states once for every producer) — healthy is
    ///   `cos ≈ −1`, so a reported `≈ +1` is a hole wound like an outer loop.
    ///
    /// All of them are one defect ("the face lies about which way it faces"), so
    /// one variant carries the measurement.
    ///
    /// This is the release-side net for the invariant `collect_planes` guards
    /// with `debug_assert`s at every boolean: a wrong flag or winding survives
    /// every other check here (edge opposition and Euler are blind to it) and
    /// walks straight into "which side is material".
    FaceMisoriented { face: Handle<Face>, cos: f64 },

    /// A cylinder's exact truth (`CylinderDef`) and its f64 cache describe **different
    /// cylinders** — the two-descriptions net: coefficients-beside-witness shows that
    /// two exact-looking descriptions of one surface can drift apart, and the cure is a
    /// consumer-side postcondition, not trust in the producer. `field` names the disagreeing
    /// quantity (`"radius"`, `"origin"`, `"dir"`, `"ref_dir"`), the two values are the def's
    /// realization and the cache's, one component at a time.
    ///
    /// The failure modes this exists to net are enormous, not subtle — a wrong seam tie-break
    /// turns `ref_dir` ~90° (a component moves by ~1), a raw-vs-unit mix-up scales `dir` by
    /// its length — while every healthy realization path sits within a few machine epsilons;
    /// the threshold sits between at [`CYL_TRUTH_EPS`].
    CylinderTruthCacheMismatch {
        surface_index: u32,
        field: &'static str,
        def_value: f64,
        cache_value: f64,
    },
}

/// How far a cylinder def's realization may sit from the cache, per component, before
/// [`Violation::CylinderTruthCacheMismatch`] fires — `|def − cache| ≤ CYL_TRUTH_EPS ·
/// max(1, |def|, |cache|)`, a mixed absolute/relative bound. Not an ulp count: the cache's
/// Gram–Schmidt legitimately smears ~1e−17 into a slot the def's shuffle keeps at exactly
/// `0.0` (measured on a full-width random axis), and ulps are meaningless across zero.
/// Healthy paths measure ≤1e−17 absolute (`nacre-topo/tests/cylinder_truth.rs` sees 0–1 ulp
/// on the nonzero components); the defects the net exists for move a component by ~10¹⁵·ε.
pub const CYL_TRUTH_EPS: f64 = 8.0 * f64::EPSILON;

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

    // The live model, not the whole append-only arena: superseded
    // cells stay in the store but drop out here, so they neither break Euler nor
    // pollute manifold use-counts. Reference integrity ran first (and would have
    // short-circuited on a dangling handle), so this traversal is in-bounds.
    let reach = model.reachable();
    let adj = Adjacency::rebuild(model); // fresh; does not trust model.adj
    check_vertex_def_carriers(model, &mut out);
    check_loop_closure(model, &reach, &mut out);
    check_manifold(model, &adj, &reach, &mut out);
    check_cavity_orientation(model, &mut out);
    check_face_orientation(model, &reach, &mut out);
    check_cylinder_truth(model, &reach, &mut out);
    check_geometric_incidence(model, &reach, &mut out);
    check_euler_poincare(model, &reach, &mut out);
    out
}

// ★★★ **There is no per-vertex tolerance to read — every vertex takes the construction
// epsilon.** The cache stores no measured residual: what it knows is *whether the coordinate was
// realized* and, if so, a per-axis bound on it. The bound is deliberately not used here. It says
// how far the **coordinate**
// is from the truth (half an ulp, or the ladder's radius); this file asks how far the cached point
// sits from the cached **carriers**, and that distance also carries the carriers' own realization
// error, which nothing records yet (`SurfaceCache` has no `tol`).
//
// ⚠ **What that costs, stated.** Arrangement-born vertices, whose measured residual is **at most
// 1.07e-14** over the census corpus, are held to [`EPS_CONSTRUCTED`], five orders looser. What a
// residual would catch — a re-named definition wearing another triple's tolerance — cannot arise:
// the cache is the realization of the definition, and the census asserts that vertex by vertex.

/// ★★★★ **This check walks the whole arena on purpose — superseded cells included.**
///
/// Every other check in this file filters by [`Reachable`], and this one deliberately does not:
/// a torn page is torn whether or not anything still points at it. An editing op supersedes old
/// cells rather than removing them, so the dead share grows with history — measured
/// at 0% for a fresh solid, ~50% after a boolean, and **99%** at depth 120 (960 of 968 vertices).
/// Filtering here would make the check quieter, not cleaner.
///
/// ⚠ **The premise is not privacy, it is `Store` being append-only** — its public surface is
/// `new·push·get·handle_at·len·is_empty·iter`, with no removal at all. So an index that was once
/// valid stays valid forever and a bounds check cannot fire on a dead cell; the violations this
/// finds are real dangling references, never arena residue.
///
/// ⚠⚠ **Do not "optimize" this by adding a `reach` filter.** Four tests plant corruption in cells
/// that are *unreachable by construction* and go red the moment this stops looking:
/// `dangling_reference_solid_shell`, `dangling_reference_vertex_definition`,
/// `a_contradictory_pierce_def_is_flagged`, `a_contradictory_vertex_def_is_flagged`.
fn check_reference_integrity(m: &Model, out: &mut Vec<Violation>) {
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let edge = m.edge(eh);
        // (There is no curve handle to check: the curve is a cache beside the store, not a
        // reference an edge can dangle.)
        {
            for v in edge.vertices {
                if v.index() as usize >= m.vertex_count() {
                    out.push(Violation::DanglingReference {
                        kind: RefKind::EdgeBoundVertex,
                        owner_index: eh.index(),
                        target_index: v.index(),
                        target_len: m.vertex_count() as u32,
                    });
                }
            }
        }
    }

    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
        // The surfaces store is private — bounds-check against its count.
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
                if he.edge.index() as usize >= m.edge_count() {
                    out.push(Violation::DanglingReference {
                        kind: RefKind::HalfEdgeEdge,
                        owner_index: fh.index(),
                        target_index: he.edge.index(),
                        target_len: m.edge_count() as u32,
                    });
                }
            }
        }
    }

    let mut i = 0u32;
    while let Some(sh) = m.shell_handle_at(i) {
        i += 1;
        let shell = m.shell(sh);
        for f in &shell.faces {
            if f.index() as usize >= m.face_count() {
                out.push(Violation::DanglingReference {
                    kind: RefKind::ShellFace,
                    owner_index: sh.index(),
                    target_index: f.index(),
                    target_len: m.face_count() as u32,
                });
            }
        }
    }

    let mut i = 0u32;
    while let Some(soh) = m.solid_handle_at(i) {
        i += 1;
        let solid = m.solid(soh);
        for sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            if sh.index() as usize >= m.shell_count() {
                out.push(Violation::DanglingReference {
                    kind: RefKind::SolidShell,
                    owner_index: soh.index(),
                    target_index: sh.index(),
                    target_len: m.shell_count() as u32,
                });
            }
        }
    }

    // Every vertex's definition references surfaces by handle (the definition is the
    // vertex, so this covers all of them, not just the discovered population).
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let vertex = m.vertex(vh);
        for s in vertex.carriers() {
            if s.index() as usize >= m.surface_count() {
                out.push(Violation::DanglingReference {
                    kind: RefKind::VertexSurface,
                    owner_index: vh.index(),
                    target_index: s.index(),
                    target_len: m.surface_count() as u32,
                });
            }
        }
    }
}

fn check_euler_poincare(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Count the live model only, not the append-only store lengths.
    let v = reach.vertices.len();
    let e = reach.edges.len();
    let f = reach.faces.len();
    let s = reach.shells.len();
    let inner_loops: usize = reach.faces.iter().map(|fh| m.face(*fh).inner.len()).sum();

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

/// Structural check: each vertex definition's carrier kinds must match its variant —
/// `ThreePlane` names planes only, `OnSeam` includes a non-plane. Runs after reference
/// integrity (it dereferences surface handles).
///
/// ★ **Full-arena, like [`check_reference_integrity`] and for the same reason** — a definition
/// that contradicts its own carriers is wrong whether or not a live solid still names that
/// vertex. Two of the four planted-corruption tests land here.
fn check_vertex_def_carriers(m: &Model, out: &mut Vec<Violation>) {
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let vertex = m.vertex(vh);
        let bad = match vertex {
            Vertex::ThreePlane(planes) => planes
                .iter()
                .any(|&s| !matches!(m.surface(s), Surface::Plane { .. })),
            Vertex::OnSeam(pair) => pair
                .iter()
                .all(|&s| matches!(m.surface(s), Surface::Plane { .. })),
            // The structure says the kinds: two planes and one cylinder, positionally.
            Vertex::Pierce {
                planes, cylinder, ..
            } => {
                planes
                    .iter()
                    .any(|&s| !matches!(m.surface(s), Surface::Plane { .. }))
                    || !matches!(m.surface(*cylinder), Surface::Cylinder { .. })
            }
        };
        if bad {
            out.push(Violation::VertexCarrierMismatch { vertex: vh });
        }
    }
}

fn check_loop_closure(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Store order (deterministic), filtered to the live faces.
    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
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
    // Resolve each half-edge to (start, end). (Every edge is bounded by type, so
    // there is no unbounded-edge arm.)
    let ends: Vec<Option<(Handle<Vertex>, Handle<Vertex>)>> = hes
        .iter()
        .map(|he| {
            let [a, b] = m.edge(he.edge).vertices;
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
    // Live edges only: a superseded edge left in the store has 0 uses
    // in the live adjacency but is not a defect — skip it. Every reachable edge
    // is referenced by a live face, so it must be used exactly twice.
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let _edge = m.edge(eh);
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
            // The stated carriers must agree with adjacency. Two different using
            // surfaces must *be* the stated pair (as a multiset). Two
            // using faces on **one** surface (a cylinder panel and its neighbouring band
            // sharing an arc) cannot spell the
            // carrier pair between them: the shared surface must be one of the stated two,
            // and the *other* stated carrier is what the curve derivation crosses it with —
            // which is the whole purpose of the check (wrong carriers ⇒ a silently wrong
            // curve cache). `[s, s]` on a plane stays reserved for cylinder seams.
            let stated = _edge.surfaces;
            let mut observed = [m.face(uses[0].0).surface, m.face(uses[1].0).surface];
            if observed[1].index() < observed[0].index() {
                observed.swap(0, 1);
            }
            let plane_self_pair =
                stated[0] == stated[1] && matches!(m.surface(stated[0]), Surface::Plane { .. });
            let agrees = if observed[0] == observed[1] {
                stated.contains(&observed[0])
            } else {
                stated == observed
            };
            if !agrees || plane_self_pair {
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
    for &sh in m.live_solids() {
        if sh.index() as usize >= m.solid_count() {
            continue;
        }
        for &cavity in &m.solid(sh).cavities {
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
    let faces = &m.shell(shell).faces;
    let first = m.face(*faces.first()?);
    let reference = m.vertex_point(loop_start(m, &first.outer)?);

    let mut flux = 0.0;
    for &fh in faces {
        let face = m.face(fh);
        let plane = match m.surface_cache(face.surface) {
            nacre_geom::Surface::Plane(plane) => plane,
            // M5 cavities are planar; a curved void is simply not checked. An arm rather than an
            // `else`, so a third surface kind is a compile error here and its author decides
            // whether its flux can be summed.
            nacre_geom::Surface::Cylinder(_) => return None,
        };
        let sign = f64::from(face.orientation.sign());
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

/// **A loop's own witness of which way its face points**, unit — the normal the
/// loop winds CCW about.
///
/// Two shapes answer, and neither reads the face's own plane, which is what
/// makes the witness independent of the thing it checks:
///
/// * a **polygon**: the Newell area vector, which no single collinear corner can
///   fool;
/// * a **closed rim**: one half-edge whose curve is a full circle — a cylinder
///   cap's loop, where there is no polygon at all. Its direction comes from the
///   *cylinder's* axis ([`nacre_topo::Model::derive_edge_curve`] builds the
///   circle from `axis.direction()`; the plane only says where the centre sits),
///   so comparing it against the cap plane's normal compares two surfaces, not a
///   value with itself. `forward` reads as it must for a rim whose start is its
///   end: along the circle's own parameterization.
///
/// `None` where the loop witnesses nothing: an unclosed chain (`OpenLoop`
/// already names that face's defect, and measuring a broken chain would only add
/// noise beside it — the more-specific-defect-first rule), fewer than three
/// points with no circle to fall back on, or a Newell sum of exactly zero.
///
/// A loop mixing **arcs** with straight edges adds each arc's witness the same way a
/// full rim is one: the chord Newell sum alone reads an arc-dominated loop **backwards** (the
/// turned boss's crescent — the region between chord and long arc lies on the chord's other
/// side), so every circle half-edge contributes its **circular segment**'s area vector —
/// `Circle::segment_area` about the circle's own normal, signed by traversal, with Δθ read from
/// the edge's stored `[from, to]` order (CCW about the axis, the arc convention;
/// `Circle::angle_of` is the one spelling). A digon (chord + arc, two points) has an empty
/// Newell sum and one segment — which is why the too-few-points refusal only applies to
/// all-line loops.
fn loop_winding(m: &Model, lp: &Loop) -> Option<Vector3> {
    if let [he] = lp.half_edges[..] {
        if let nacre_geom::Curve::Circle(c) = m.edge_curve(he.edge) {
            return Some(c.normal() * if he.forward { 1.0 } else { -1.0 });
        }
    }
    let pts = loop_points(m, lp);
    let n = pts.len();
    let arcs: Vec<&nacre_topo::HalfEdge> = lp
        .half_edges
        .iter()
        .filter(|he| matches!(m.edge_curve(he.edge), nacre_geom::Curve::Circle(_)))
        .collect();
    if n < 3 && arcs.is_empty() {
        return None;
    }
    let ends = |he: &nacre_topo::HalfEdge| {
        let [a, b] = m.edge(he.edge).vertices;
        if he.forward { (a, b) } else { (b, a) }
    };
    let hes = &lp.half_edges;
    if (0..n).any(|i| ends(&hes[i]).1 != ends(&hes[(i + 1) % n]).0) {
        return None;
    }
    let chords = (0..n).fold(Vector3::zero(), |acc, i| {
        acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
    });
    let sum = arcs.iter().fold(chords, |acc, he| {
        let nacre_geom::Curve::Circle(c) = m.edge_curve(he.edge) else {
            unreachable!("filtered above");
        };
        let [va, vb] = m.edge(he.edge).vertices;
        let dt = (c.angle_of(m.vertex_point(vb)) - c.angle_of(m.vertex_point(va)))
            .rem_euclid(std::f64::consts::TAU);
        let sign = if he.forward { 1.0 } else { -1.0 };
        // The Newell fold above is twice the area vector; scale the segment to match.
        acc + c.normal() * (2.0 * sign * c.segment_area(dt))
    });
    sum.normalize()
}

/// Each live planar face's loops must stand the way the face says it faces: the
/// **outer** loop winds CCW about the stated outward normal, and every **inner**
/// loop the other way round (the rule `build_prism` states once, for every
/// producer). The loops are the one independent witness of "which way is out",
/// so [`loop_winding`] is compared against `plane.normal()` × `orientation`.
///
/// Skipped where a face cannot be asked: **non-planar** ones — a cylinder's
/// lateral normal changes from point to point, so "the face's normal" is not a
/// question — and loops that witness nothing ([`loop_winding`] returns `None`).
fn check_face_orientation(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
        if !reach.faces.contains(&fh) {
            continue;
        }
        let plane = match m.surface_cache(face.surface) {
            nacre_geom::Surface::Plane(plane) => plane,
            // A cylinder's lateral normal changes from point to point, so "the face's normal" is
            // not a question. An arm rather than an `else`, so a third surface kind is a compile
            // error here and its author decides whether it can be asked.
            nacre_geom::Surface::Cylinder(_) => continue,
        };
        let stated = plane.normal() * f64::from(face.orientation.sign());
        // Two aligned unit vectors sit at ±1; 0.5 is the same "a full unit from
        // the sign boundary" margin the boolean's own asserts use. A middling
        // `cos` is reported rather than skipped — for a rim it means the circle
        // does not lie in this face's plane, which is as much a defect as a
        // flipped flag, and one more silent skip is what let a cylinder ship
        // with both caps facing the same way.
        if let Some(dir) = loop_winding(m, &face.outer) {
            let cos = dir.dot(stated);
            if cos <= 0.5 {
                out.push(Violation::FaceMisoriented { face: fh, cos });
            }
        }
        for hole in &face.inner {
            if let Some(dir) = loop_winding(m, hole) {
                let cos = dir.dot(stated);
                if cos >= -0.5 {
                    out.push(Violation::FaceMisoriented { face: fh, cos });
                }
            }
        }
    }
}

/// The start vertex of a loop's first half-edge (`None` for an empty loop).
fn loop_start(m: &Model, lp: &Loop) -> Option<Handle<Vertex>> {
    let he = lp.half_edges.first()?;
    let [a, b] = m.edge(he.edge).vertices;
    Some(if he.forward { a } else { b })
}

/// A loop's traversal-order start points, one per half-edge.
fn loop_points(m: &Model, lp: &Loop) -> Vec<Point3> {
    lp.half_edges
        .iter()
        .map(|he| {
            let [a, b] = m.edge(he.edge).vertices;
            m.vertex_point(if he.forward { a } else { b })
        })
        .collect()
}

/// `(unsigned area, area-weighted centroid)` of a planar polygon loop, via a
/// signed triangle fan from the first vertex (exact for concave loops). `None`
/// if the loop has an unbounded edge or fewer than three vertices.
fn loop_area_centroid(m: &Model, lp: &Loop) -> Option<(f64, Point3)> {
    let pts = loop_points(m, lp);
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

/// **The def and the cache must describe one cylinder** — see
/// [`Violation::CylinderTruthCacheMismatch`]. The radius is compared on every cylinder (it is
/// invariant under every rigid motion, and a mirrored cylinder cannot exist); the frame —
/// origin, axis direction, seam direction — only where the statement is world-spoken
/// (`motion: None`): a recorded chain states the def *before* the motion, and realizing a
/// statement through its chain is machinery this net deliberately does not build.
fn check_cylinder_truth(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    let mut seen: std::collections::HashSet<Handle<Surface>> = std::collections::HashSet::new();
    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
        if !reach.faces.contains(&fh) || !seen.insert(face.surface) {
            continue;
        }
        let cy = match m.surface_cache(face.surface) {
            nacre_geom::Surface::Cylinder(cy) => cy,
            // A plane has no truth-versus-cache residual to check: its truth is exact
            // coefficients, and its orientation is `check_face_orientation`'s question. An arm
            // rather than an `else`, so a third surface kind is a compile error here and its
            // author decides what its truth check is.
            nacre_geom::Surface::Plane(_) => continue,
        };
        let (def, motion) = match m.surface(face.surface) {
            nacre_topo::Surface::Cylinder { def, motion } => (def, motion),
            // The stores are index-parallel with one entry door, so a kind mismatch cannot
            // arise; it is transform's `unreachable!`, not this check's proposition.
            nacre_topo::Surface::Plane { .. } => continue,
        };
        let surface_index = face.surface.index();
        let mut flag = |field: &'static str, dv: f64, cv: f64| {
            let scale = 1f64.max(dv.abs()).max(cv.abs());
            if (dv - cv).abs() > CYL_TRUTH_EPS * scale {
                out.push(Violation::CylinderTruthCacheMismatch {
                    surface_index,
                    field,
                    def_value: dv,
                    cache_value: cv,
                });
            }
        };
        flag("radius", def.radius_f64(), cy.radius());
        if motion.is_some() {
            continue;
        }
        let co = cy.axis().origin().as_array();
        for (k, o) in def.origin().iter().enumerate() {
            flag("origin", o.to_f64(), co[k]);
        }
        // The raw exact directions are positively parallel to the cache's unit ones, so their
        // own realizations, normalized, must land beside them. (No `Rat` in a signature here —
        // scalar types flow through topo's API, and this crate deliberately never names them.)
        for (raw, cache_v, name) in [
            (def.dir(), cy.axis().direction(), "dir"),
            (def.ref_dir(), cy.ref_dir(), "ref_dir"),
        ] {
            let realized = Vector3::from_array([raw[0].to_f64(), raw[1].to_f64(), raw[2].to_f64()]);
            if let Some(u) = realized.normalize() {
                for k in 0..3 {
                    flag(name, u[k], cache_v[k]);
                }
            }
        }
    }
}

fn check_geometric_incidence(m: &Model, reach: &Reachable, out: &mut Vec<Violation>) {
    // Each live edge's bound vertices must lie on the edge's curve.
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let edge = m.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        {
            let [a, b] = edge.vertices;
            let curve = m.edge_curve(eh);
            for vh in [a, b] {
                let residual = curve.distance(m.vertex_point(vh));
                // The floor and the vertex term are the same constant: no vertex carries a
                // measured residual.
                let tol = EPS_CONSTRUCTED;
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
    let mut i = 0u32;
    while let Some(fh) = m.face_handle_at(i) {
        i += 1;
        let face = m.face(fh);
        if !reach.faces.contains(&fh) {
            continue;
        }
        // ★ The **cache** is the proposition here, not a shortcut: this asks how far a cached
        // vertex coordinate sits from its cached carrier — one rounded description against
        // another. The truth-side question (does the definition match its carriers' kinds) is
        // `check_vertex_def_carriers`, and it asks the truth.
        let surface = m.surface_cache(face.surface);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                {
                    let [a, b] = m.edge(he.edge).vertices;
                    let vh = if he.forward { a } else { b };
                    let residual = surface.distance(m.vertex_point(vh));
                    // ★ Plus what the residual's *own* arithmetic can produce. The epsilon says
                    // where the vertex may sit; it says nothing about `Surface::distance`, so a
                    // vertex exactly on the surface can still report a machine-scale residual and
                    // be flagged for it. Measured: a four-plane concurrency whose vertex is
                    // genuinely on the plane sat at residual *exactly equal* to its tolerance, and
                    // passed only because the comparison is strict — the slack that had been
                    // covering this term was the loose measured tolerances of the day.
                    let tol = EPS_CONSTRUCTED + surface.distance_eps(m.vertex_point(vh));
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
    // **Every** live vertex must lie within tolerance of every surface its definition claims
    // it meets (the definition is the truth, the point is a within-tol cache).
    // Reference integrity ran first, so the definition's surface handles are in bounds.
    //
    // ★ This covers all of them, not just the measured population: the definition *is* the
    // vertex, so a constructed corner's claim is checked too, and the checker asks the same
    // question of every one — with the tolerance the vertex has (measured) or the construction
    // epsilon (built). The doctrine this serves: the surfaces a vertex is defined by
    // are self-evidently near it, and the real question is whether the *coordinate* still
    // matches the definition after everything that has happened to it.
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let vertex = m.vertex(vh);
        if !reach.vertices.contains(&vh) {
            continue;
        }
        let tol = EPS_CONSTRUCTED;
        for sh in vertex.carriers() {
            // ★ Cache against cache, deliberately — `vertex_point` is the cached coordinate,
            // so its residual has to be measured against the cached carrier for the two to be
            // the same description. Against the truth this would report the realization error
            // of both, which is a different quantity and has no threshold here.
            let residual = m.surface_cache(sh).distance(m.vertex_point(vh));
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
#[path = "tests/lib.rs"]
mod tests;
