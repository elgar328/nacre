//! The plane/face substrate: per-face (`FaceInfo`) and per-plane-class (`PlaneGeom`) tables and
//! their construction. Everything the boolean engine and its combinatorial queries build on.

use crate::combinatorics;
use crate::{BoolError, RejectReason, he_start, reject, tolerant};
use nacre_cip::Pt3;
use nacre_geom::intersect::{plane_plane, planes_coplanar};
use nacre_geom::{Plane, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, HalfEdge, Model, Orientation, Origin, Shell, Solid, Vertex};
use std::collections::HashMap;

/// A face's supporting plane plus the exact in/out data the seam path needs.
///
/// `three_plane_orient3d(.., tri[0], tri[1], tri[2])` returns `+1` when the
/// implicit point lies on **`tri`'s right-hand-normal side** — the convention is
/// tied to the triangle, never to `plane`. `n_out` happens to equal that RH normal
/// only because `tri` is taken outer-CCW; `plane.normal()` is the *surface's*
/// normal and may point inward on a `Reversed` face. Every sign test here reads
/// `n_out` (or `tri`), and none reads `plane.normal()`.
pub(crate) struct FaceInfo {
    pub(crate) surf: Handle<Surface>,
    /// The face this plane came from. Distinguishes two coplanar faces that share one
    /// `Surface` (a Cut splits one face into disjoint pieces reusing its surface —
    /// cell coplanar-narrow), which `surf` alone collapses. `surf_ix` keys on this.
    pub(crate) face: Handle<Face>,
    pub(crate) plane: Plane,
    /// Three non-collinear outer-loop points, **ordered so their RH normal is outward**.
    /// The order need not follow the loop: at a reflex corner it is reversed.
    pub(crate) tri: [Point3; 3],
    /// Outward normal, `(tri[1]−tri[0])×(tri[2]−tri[0])` normalized — the single
    /// source of "outward" for both the in/out sign test and face ordering.
    pub(crate) n_out: Vector3,
    /// `+1` when this face's stored plane normal already points out of its solid, `-1` when the
    /// face is `Reversed` and the two oppose.
    ///
    /// **This face's**, not its plane class's. The class-frame twin is [`PlaneGeom::frame_sign`],
    /// and the two used to be one function called with either kind of index — the single place the
    /// face/plane convention could not be asserted, because both readings were legitimate
    /// (dev-log, normalization cell). Separate names, separate questions.
    pub(crate) orient_sign: i8,
    /// The three `tri` points as **exact `Pt3` definitions**, in the same order as `tri`.
    /// Built once here and borrowed by every predicate (`plane_def`) — it used to be rebuilt
    /// per judgment, which dominated the boolean's runtime.
    pub(crate) tri_pt3: [Pt3; 3],
    /// Whether this face's solid is rotated — the predicate-routing signal, decided by
    /// [`solid_is_rotated`].
    ///
    /// **Set together with `tri_pt3`, and only here.** They used to be one field (`Option`),
    /// whose emptiness meant "not rotated"; that conflation is what stopped the definition from
    /// being cached. Do not derive this from the definition: a rotated solid's face can witness
    /// its plane through chain-less points, and a 90°-family rotation has tol exactly 0.
    pub(crate) rotated: bool,
}

/// The supporting planes of a solid's outer shell. `Unsupported` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<FaceInfo>, BoolError> {
    // A rotated operand's face coordinates are rounded, so each plane also carries its
    // exact `Pt3` definition (overhaul stage 3). Decided once per solid — the axis-aligned
    // path keeps `tri_pt3 = None` and pays nothing.
    let rotated = solid_is_rotated(model, solid);
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let plane = match model.surfaces.get(face.surface) {
                Surface::Plane(p) => *p,
                Surface::Cylinder(_) => return Err(reject(RejectReason::CylinderFace)),
            };
            let (tri, tri_verts) =
                outer_tri(model, face).ok_or_else(|| reject(RejectReason::DegenerateFace))?;
            let n_out = (tri[1] - tri[0])
                .cross(tri[2] - tri[0])
                .normalize()
                .ok_or_else(|| reject(RejectReason::DegenerateNormal))?;
            // `tri_pt3` and `rotated` are set here, together, and nowhere else. `plane_table`
            // copies the pair from a class root; nothing else constructs either.
            let tri_pt3 = if rotated {
                // Own vertices first: a surviving operand corner assembles directly, so a plain
                // rotated operand keeps `tri_pt3 == tri`. A seam-dominated face (its `tri` verts
                // are rotated seams) instead witnesses its plane through provenance — its plane is
                // `R(π)` for an operand plane `π`, recovered from the operand face still on `π`.
                let own = (|| {
                    Some([
                        crate::rotated_vertex::vertex_pt3(model, tri_verts[0]).ok()?,
                        crate::rotated_vertex::vertex_pt3(model, tri_verts[1]).ok()?,
                        crate::rotated_vertex::vertex_pt3(model, tri_verts[2]).ok()?,
                    ])
                })();
                match own {
                    Some(t) => t,
                    None => {
                        let mut w = crate::rotated_vertex::face_plane_witness(model, face)
                            .map_err(|_| reject(RejectReason::RotatedUnderdetermined))?;
                        // `tri_pt3` is an *oriented* plane witness: the own-vertex path inherits
                        // outward order from `outer_tri`, so a provenance witness must be wound to
                        // agree with this face's outward normal `n_out` too, or the implicit-point
                        // `orient3d` reads the plane's opposite side and flips every sign on it.
                        let e1 = Vector3::from_array(w[1].coord) - Vector3::from_array(w[0].coord);
                        let e2 = Vector3::from_array(w[2].coord) - Vector3::from_array(w[0].coord);
                        if e1.cross(e2).dot(n_out) < 0.0 {
                            w.swap(1, 2);
                        }
                        w
                    }
                }
            } else {
                // An axis-aligned face's `tri` coordinates are already exact f64, so the
                // definition is stated rather than measured (`Pt3::exact`). Building it here —
                // rather than per judgment — is the whole point of this field; it costs one
                // construction per face and no high-precision arithmetic.
                let e = |p: Point3| {
                    Pt3::exact(p.as_array())
                        .ok_or_else(|| reject(RejectReason::CoordinateOutOfRange))
                };
                [e(tri[0])?, e(tri[1])?, e(tri[2])?]
            };
            // `orient_sign`, precomputed: the two invariants it used to re-check on every call
            // are properties of this face, so they are decided once, here.
            let dot = plane.normal().dot(n_out);
            debug_assert!(
                dot.abs() > 0.5,
                "a plane's normal must be parallel to n_out"
            );
            debug_assert_eq!(
                dot > 0.0,
                face.orientation == Orientation::Forward,
                "n_out's sign against the surface normal is the face's orientation"
            );
            out.push(FaceInfo {
                surf: face.surface,
                face: fh,
                plane,
                tri,
                n_out,
                orient_sign: if dot > 0.0 { 1 } else { -1 },
                tri_pt3,
                rotated,
            });
        }
    }
    Ok(out)
}

/// All shells of a solid — outer first, then cavities. The boolean seam
/// front-end walks these so a cavitied operand's void walls are seen (cell
/// (5c-in)); a non-hollow solid yields just its outer shell, unchanged.
pub(crate) fn solid_shell_handles(model: &Model, solid: Handle<Solid>) -> Vec<Handle<Shell>> {
    let s = model.solids.get(solid);
    std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect()
}

/// Three non-collinear points of a face's outer loop — with the **vertex handle** each
/// point came from — ordered so their right-hand normal points **out** of the solid.
/// The handles let the toleranced predicates rebuild each point as a `Pt3` (overhaul
/// stage 3); the coordinates alone drive the axis-aligned path.
pub(crate) fn outer_tri(model: &Model, face: &Face) -> Option<([Point3; 3], [Handle<Vertex>; 3])> {
    let verts: Vec<Handle<Vertex>> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| he_start(model, he))
        .collect();
    let pts: Vec<Point3> = verts
        .iter()
        .map(|&vh| model.vertices.get(vh).point)
        .collect();
    let n = pts.len();
    // The turn at one corner does not know which way the ring winds. Every b-rep loop is
    // CCW about its face's outward normal, but at a *reflex* corner the local turn
    // opposes the global winding, so three consecutive points can hand back an inward
    // normal. The Newell sum has no single corner to be fooled by.
    let newell = (0..n).fold(Vector3::zero(), |acc, i| {
        acc + (pts[i] - pts[0]).cross(pts[(i + 1) % n] - pts[0])
    });
    let i = (0..n).find(|&i| {
        let (a, b, c) = (pts[i], pts[(i + 1) % n], pts[(i + 2) % n]);
        (b - a).cross(c - a).norm() > 0.0
    })?;
    let (i0, i1, i2) = (i, (i + 1) % n, (i + 2) % n);
    let (a, b, c) = (pts[i0], pts[i1], pts[i2]);
    // Same b/c swap for coords and handles, so `tri[k]` and `tri_verts[k]` stay aligned.
    Some(if (b - a).cross(c - a).dot(newell) < 0.0 {
        ([a, c, b], [verts[i0], verts[i2], verts[i1]])
    } else {
        ([a, b, c], [verts[i0], verts[i1], verts[i2]])
    })
}

/// Max distance of `p` to its 3 planes and 3 pairwise lines (the measured
/// `Origin::Discovered` tolerance).
pub(crate) fn vertex_tol(p: Point3, a: &Plane, b: &Plane, c: &Plane) -> f64 {
    let mut tol = a.distance(p).max(b.distance(p)).max(c.distance(p));
    for (x, y) in [(a, b), (a, c), (b, c)] {
        if let Some(line) = plane_plane(x, y) {
            tol = tol.max(line.distance(p));
        }
    }
    tol
}

/// The minimal per-op plane table two solids share: the
/// concatenated plane list (`a`'s then `b`'s), the face→index map, and each solid's
/// [`combinatorics::EdgeFaces`]. Built once and shared: indices into the returned `planes`/`surf_ix`
/// are common to both solids, so a vertex of `a` and a face of `b` compose in one index space.
/// Destructure it with `..` (`let PlaneSetup { planes: faces_tab, geom: planes, plane_ix, .. } = …`):
/// the tables here grow as the arrangement learns to say "plane" and "face" in different index
/// spaces, and a positional tuple made every one of those steps touch all ~25 call sites.
///
/// The plane classes (`canon`) are computed here to build `geom`/`plane_ix` and then dropped — the
/// dense `plane_ix` is the only face→plane map anything downstream needs, so the sparse union-find
/// output does not escape.
pub(crate) struct PlaneSetup {
    pub(crate) planes: Vec<FaceInfo>,
    pub(crate) surf_ix: HashMap<Handle<Face>, usize>,
    pub(crate) inc_a: combinatorics::EdgeFaces,
    pub(crate) inc_b: combinatorics::EdgeFaces,
    /// The arrangement's planes, densely indexed — see [`dense_planes`].
    pub(crate) geom: Vec<PlaneGeom>,
    /// `plane_ix[face]` is that face's plane, as an index into `geom`.
    pub(crate) plane_ix: Vec<usize>,
}

pub(crate) fn plane_index_setup(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<PlaneSetup, BoolError> {
    let mut planes = collect_planes(model, a)?;
    planes.extend(collect_planes(model, b)?);
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        surf_ix.insert(pi.face, i);
    }
    let inc_a = combinatorics::edge_faces(model, a, &surf_ix)?;
    let inc_b = combinatorics::edge_faces(model, b, &surf_ix)?;
    let canon = plane_classes(&planes);
    let (geom, plane_ix) = dense_planes(&planes, &canon);
    Ok(PlaneSetup {
        planes,
        surf_ix,
        inc_a,
        inc_b,
        geom,
        plane_ix,
    })
}

/// One plane of the arrangement, indexed by a **dense** class id.
///
/// The face table cannot answer "which plane" without a convention: a class holds faces from both
/// operands, and two of them can face opposite ways, so there is no such thing as *the* plane's
/// outward normal. What a plane has is a **frame** — the class root's stored normal — and the only
/// direction fact anyone needs from it is [`PlaneGeom::frame_sign`]. Everything else here is a
/// witness: three points known to lie on this plane, used to reconstruct it exactly.
/// The pre-rotation twin of a witness triangle: its `chain_id`, base points and base plane.
///
/// A rigid motion preserves the determinants the predicates take, so a judgement whose inputs all
/// carry **one** motion can be answered on these instead — exactly, off the toleranced path
/// entirely. `None` for the base data when a base coordinate is not `f64`-representable, since the
/// exact predicate takes `f64`; the judgement then stays toleranced (slower, never wrong).
#[derive(Clone, Copy, Debug)]
pub(crate) struct BaseFrame {
    /// `0` = no motion. Equal only for structurally identical chains.
    pub(crate) chain_id: u64,
    pub(crate) tri: Option<[Point3; 3]>,
    pub(crate) coeffs: Option<[f64; 4]>,
}

impl BaseFrame {
    /// No motion to cancel — for hand-built tables in tests.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Self {
            chain_id: 0,
            tri: None,
            coeffs: None,
        }
    }

    fn of(tri_pt3: &[Pt3; 3], frame_sign: i8) -> Self {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        // All three witnesses must carry the *same* chain for the plane to have one motion.
        for (k, p) in tri_pt3.iter().enumerate() {
            if k > 0 && p.chain.len() != tri_pt3[0].chain.len() {
                return Self {
                    chain_id: 0,
                    tri: None,
                    coeffs: None,
                };
            }
            if k > 0 {
                for (a, b) in p.chain.iter().zip(tri_pt3[0].chain.iter()) {
                    if a.axis != b.axis || a.angle != b.angle || a.point != b.point {
                        return Self {
                            chain_id: 0,
                            tri: None,
                            coeffs: None,
                        };
                    }
                }
            }
        }
        if tri_pt3[0].chain.is_empty() {
            return Self {
                chain_id: 0,
                tri: None,
                coeffs: None,
            }; // no motion to cancel
        }
        for n in tri_pt3[0].chain.iter() {
            (n.axis as u8).hash(&mut h);
            format!("{:?}", n.angle).hash(&mut h);
            for c in n.point {
                format!("{c:?}").hash(&mut h);
            }
        }
        // 0 is reserved for "no motion", so never hand it out as an id.
        let chain_id = h.finish() | 1;
        let exact = |r: nacre_scalar::Rat| nacre_scalar::Rat::try_from_f64(r.to_f64()) == Some(r);
        if !tri_pt3.iter().all(|p| p.base.iter().all(|&r| exact(r))) {
            return Self {
                chain_id,
                tri: None,
                coeffs: None,
            };
        }
        let pt = |p: &Pt3| {
            Point3::from_array([p.base[0].to_f64(), p.base[1].to_f64(), p.base[2].to_f64()])
        };
        let tri = [pt(&tri_pt3[0]), pt(&tri_pt3[1]), pt(&tri_pt3[2])];
        // ★ The base plane must carry the **stored** orientation, not the triangle's. A class's
        // stored normal and its witness triangle's `cross` can oppose — that is exactly what
        // `frame_sign` records — and `through_points` gives the triangle's. Rotation preserves the
        // cross product (`det(R) = 1`), so multiplying by `frame_sign` reproduces the same relation
        // in the base frame. Without it the exact path answers with a flipped sign, which the suite
        // caught immediately.
        let coeffs = Plane::through_points(tri[0], tri[1], tri[2]).map(|pl| {
            let c = pl.coefficients();
            let k = f64::from(frame_sign);
            [c[0] * k, c[1] * k, c[2] * k, c[3] * k]
        });
        Self {
            chain_id,
            tri: Some(tri),
            coeffs,
        }
    }
}

pub(crate) struct PlaneGeom {
    pub(crate) plane: Plane,
    /// The class's representative surface — what `assemble_fuse_cut` records in a
    /// `VertexDef::ThreePlane`.
    pub(crate) surf: Handle<Surface>,
    /// Witness points on this plane (the root face's `tri`), outward-ordered for that face.
    pub(crate) tri: [Point3; 3],
    /// The witness as exact `Pt3` definitions (the root face's), borrowed by every predicate.
    pub(crate) tri_pt3: [Pt3; 3],
    /// Whether the root face's solid is rotated — copied from it together with `tri_pt3` so the
    /// pair cannot disagree. See [`FaceInfo::rotated`].
    pub(crate) rotated: bool,
    /// `+1` when the plane's stored normal agrees with the root face's outward normal, `-1` when
    /// they oppose. This *is* the label frame: `[A_above, A_below, …]` is defined about the class
    /// root's stored normal, and this sign is what relates it to material. Precomputed here so the
    /// two `debug_assert`s that guard the convention run once, at construction.
    pub(crate) frame_sign: i8,
    /// The pre-rotation twin — see [`BaseFrame`].
    pub(crate) base: BaseFrame,
}

/// Dense plane ids for a face table: `(geom, plane_ix)` where `plane_ix[face]` indexes `geom`.
///
/// **The numbering is monotone in `canon`.** Roots are ranked in increasing order, so
/// `canon[i] < canon[j]` iff `plane_ix[i] < plane_ix[j]` — every comparison, sort and lex-min over
/// plane indices is order-isomorphic to the sparse form. Nothing found in the engine turns out to
/// depend on that (the two candidates — `loop_winding`'s lex-min node and `crossings`' pre-dedup
/// sort — are by coordinate and by set, respectively), but the audit cannot be proved exhaustive
/// over ~175 sites, so the numbering removes the question instead of answering it.
pub(crate) fn dense_planes(planes: &[FaceInfo], canon: &[usize]) -> (Vec<PlaneGeom>, Vec<usize>) {
    let mut roots: Vec<usize> = canon.to_vec();
    roots.sort_unstable();
    roots.dedup();
    let plane_ix = canon
        .iter()
        .map(|c| {
            roots
                .binary_search(c)
                .expect("a class root is in the root set")
        })
        .collect();
    let geom = roots
        .iter()
        .map(|&r| {
            let pi = &planes[r];
            PlaneGeom {
                base: BaseFrame::of(&pi.tri_pt3, pi.orient_sign),
                plane: pi.plane,
                surf: pi.surf,
                tri: pi.tri,
                tri_pt3: pi.tri_pt3.clone(),
                rotated: pi.rotated,
                frame_sign: pi.orient_sign,
            }
        })
        .collect();
    (geom, plane_ix)
}

/// Whether three `Pt3` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected `RAY_DEGENERATE`, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, so this is the exact
/// zero-area (collinear) test on the vertices' rotation definitions.
#[cfg(test)]
pub(crate) fn pt3_base_collinear(a: &Pt3, b: &Pt3, c: &Pt3) -> bool {
    use nacre_scalar::Rat;
    let (a, b, c) = (&a.base, &b.base, &c.base);
    let proj_zero = |i: usize, j: usize| -> Option<bool> {
        let det = b[i]
            .checked_sub(a[i])?
            .checked_mul(c[j].checked_sub(a[j])?)?
            .checked_sub(
                b[j].checked_sub(a[j])?
                    .checked_mul(c[i].checked_sub(a[i])?)?,
            )?;
        Some(det == Rat::from_int(0))
    };
    matches!(
        (proj_zero(1, 2), proj_zero(2, 0), proj_zero(0, 1)),
        (Some(true), Some(true), Some(true))
    )
}

/// Each outer-shell edge with its bound vertices and the two combined-plane
/// indices of its adjacent faces, in first-seen (deterministic) order.
///
/// Every loop of every face is walked, holes included: a hole-ring edge is used
/// once by the holed face's inner loop and once by the neighbouring wall's outer
/// loop, so it too has exactly two incident faces. Walking `outer` before `inner`
/// on each face leaves the order of a hole-free solid untouched.
///
/// The pair is returned as `[usize; 2]`, so no caller can index a third slot: an
/// edge with any other incidence count is a non-manifold shell and rejects here.
#[allow(clippy::type_complexity)]
pub(crate) fn edge_incidence(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<Vec<(Handle<Edge>, [Handle<Vertex>; 2], [usize; 2])>, BoolError> {
    let mut order: Vec<Handle<Edge>> = Vec::new();
    let mut map: HashMap<Handle<Edge>, ([Handle<Vertex>; 2], Vec<usize>)> = HashMap::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shells.get(sh).faces {
            let face = model.faces.get(fh);
            let pidx = surf_ix[&fh];
            for he in face_half_edges(face) {
                let bounds = model.edges.get(he.edge).bounds.expect("bounded");
                let entry = map.entry(he.edge).or_insert_with(|| {
                    order.push(he.edge);
                    (bounds, Vec::new())
                });
                entry.1.push(pidx);
            }
        }
    }
    order
        .into_iter()
        .map(|e| {
            let (b, p) = map.remove(&e).unwrap();
            match p[..] {
                [x, y] => Ok((e, b, [x, y])),
                // `validate` would call this `NonOpposedEdge`, but `boolean` never runs
                // `validate` on its inputs, so the guard stays. No firing test.
                _ => Err(reject(RejectReason::NonManifoldEdge)),
            }
        })
        .collect()
}

/// Every half-edge of a face: its outer loop first, then each hole ring in order.
pub(crate) fn face_half_edges(face: &Face) -> impl Iterator<Item = &HalfEdge> {
    face.outer
        .half_edges
        .iter()
        .chain(face.inner.iter().flat_map(|l| l.half_edges.iter()))
}

/// Two faces lie on the same plane — by a **shared `Surface` handle** (explicit
/// sharing: O(1) `Handle` identity, exact, rotation-independent) or, as a fallback,
/// by the geometric rank-1 `planes_coplanar` test. A referenced coplanar contact —
/// a pad/pocket cap that reuses its face's surface — is caught by the handle path
/// without any coordinate test. On the axis-aligned M5 corpus the handle path is
/// redundant with `planes_coplanar` (same handle ⇒ same plane), so the geometric
/// fallback is what keeps independently-built coplanar contacts working; the handle
/// path's real payoff is rotated frames, where the geometric test would need the
/// rotation-exact judgment.
pub(crate) fn shares_or_coplanar(planes: &[FaceInfo], i: usize, j: usize) -> bool {
    let (pa, pb) = (&planes[i], &planes[j]);
    // Three independent witnesses, OR-ed, so this can only ever merge *more* than before:
    //  1. the same `Surface` handle — coplanar by reference (what an ops-built tool's base cap and
    //     its target face share, and what a chained operand's split coplanar faces share);
    //  2. exactly proportional coefficients — the original test, kept;
    //  3. the faces' own coordinates, exactly (`t_planes_coplanar`) — the only one of the three
    //     that does not read a *derived* value, and the one that catches two independently built
    //     solids whose walls coincide (`add_cuboid` stacked on `add_cuboid`), where the rounded
    //     coefficients of differently-sized faces are not exactly proportional.
    pa.surf == pb.surf
        || planes_coplanar(&pa.plane, &pb.plane)
        || tolerant::t_planes_coplanar(planes, i, j)
}

/// Union-find root of `x` in `parent` (with path compression). Roots are the smallest index
/// of their class, so the result is deterministic (replay, DNA §absolute-3).
/// Union-find root with path compression. Drives component grouping in [`unify_coplanar_faces`].
pub(crate) fn uf_find(parent: &mut [usize], x: usize) -> usize {
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

/// Canonicalize the combined plane table by coplanarity: two planes that are the same plane
/// (shared `Surface` handle, or exact rank-1 [`planes_coplanar`]) are merged into one class, so
/// a wall of `a` coplanar with a wall of `b` names a **single line** in a shared plane π. This is
/// the one thing the seam engine cannot do (it rejects `order_along(R,R)==0` as `FOURPLANE`);
/// canonicalizing turns that self-comparison into a real order. Returns `canon` where `canon[i]`
/// is the class root (the smallest index in the class). Every decision is exact
/// (`shares_or_coplanar`) — no coordinate. O(n²) scan over the (small) face count.
// Wired into the unified coplanar handler's dispatch in a later cell; used by tests now.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn plane_classes(planes: &[FaceInfo]) -> Vec<usize> {
    let n = planes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            if shares_or_coplanar(planes, i, j) {
                let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                if ri != rj {
                    // Attach the larger root under the smaller so a class's root is its min index.
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }
    (0..n).map(|i| uf_find(&mut parent, i)).collect()
}

/// Whether `solid` was produced by a non-exact rotation — its vertices carry
/// `Origin::Rotated`. A `Transform` rotates a whole solid uniformly and `boolean`
/// rejects rotated inputs, so a solid is all-or-nothing rotated: one vertex decides
/// (O(1)). (90°-family rotations stay exact/`Constructed`, so this is false for them.)
pub(crate) fn solid_is_rotated(model: &Model, solid: Handle<Solid>) -> bool {
    let sh = model.solids.get(solid).outer;
    for &fh in &model.shells.get(sh).faces {
        for he in &model.faces.get(fh).outer.half_edges {
            if let Some(bounds) = model.edges.get(he.edge).bounds {
                return matches!(model.vertices.get(bounds[0]).origin, Origin::Rotated { .. });
            }
        }
    }
    false
}
