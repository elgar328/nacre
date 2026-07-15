//! Toleranced indirect predicates over the b-rep `Model` (overhaul §TIP, stage 2b).
//!
//! A read-only analysis above `nacre-topo` (like `nacre-validate` / `nacre-props`): it
//! assembles a rotated vertex's exact **definition** + **directional tol** from the
//! rotation forest ([`nacre_scalar::frame3::Pt3`]) and judges `orient3d` over kernel
//! vertices. tol-0 (`Constructed`, unrotated or 90°-family) configs take the exact
//! Shewchuk path (`nacre-predicates`); any tol > 0 (rotated) takes the filter →
//! astro-float escalation ([`nacre_scalar::frame3::orient3d_judge`]).
//!
//! **Scope (stage 2b + 2c-ii): direct predicate over Constructed-rooted,
//! pure-rotation points, plus the indirect bridge for a single `Discovered` seam.** A
//! `Discovered` seam vertex is a three-plane intersection (an implicit point), judged
//! by the *indirect* predicate from its planes' definitions, not its f64 cache. When
//! exactly one of the four `orient3d` vertices is `Discovered`, [`orient3d`]
//! reconstructs each of its three planes from three rotated points of a face on that
//! surface and routes to [`nacre_scalar::frame3::indirect_orient3d_judge`] (stage
//! 2c-ii). Two or more implicit points, a vertex translated after being rotated, and a
//! plane with too few direct vertices all surface as a typed [`TipError`] (honest
//! defer) rather than a silent wrong sign. Not yet wired into boolean (stage 3): the
//! boolean still rejects non-90° rotated inputs, so a rotated seam is exercised
//! synthetically.

use nacre_geom::Surface;
use nacre_scalar::frame3::{Pt3, RotNode, indirect_orient3d_judge, orient3d_judge};
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Rotation, Vertex, VertexDef};

/// Why a vertex could not be judged by the direct predicate here (honest defer).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TipError {
    /// A `Discovered` seam (or a point rooted on one): its exact position is a
    /// three-plane intersection, judged by the indirect predicate (stage 2c), not by a
    /// direct orient3d over its f64 cache.
    IndirectRequired,
    /// The vertex was translated after being rotated. The rotation forest records
    /// rotations only (§⑦), so base + chain no longer reproduce the kernel coordinate
    /// and the exact position is not recoverable here.
    TranslateInterleaved,
    /// An f64 coordinate did not fit an exact i128 rational (§4 downgrade).
    Downgrade,
    /// A `Discovered` seam's defining surface has fewer than three direct, non-collinear
    /// vertices on its faces, so its plane cannot be reconstructed for the indirect
    /// predicate (stage 2c-ii).
    PlaneUnderdetermined,
}

/// The toleranced point ([`Pt3`]) of a kernel vertex, assembled from its root
/// coordinate and rotation forest. `Err` when it is not a direct (Constructed-rooted,
/// pure-rotation) point — see [`TipError`].
pub fn vertex_pt3(model: &Model, vh: Handle<Vertex>) -> Result<Pt3, TipError> {
    let pt3 = build(model, vh)?;
    // Pure-rotation guard: the forest omits translations, so a rotate-then-translate
    // history leaves base + chain ≠ the kernel coord (and hp_coord would drop the
    // translation). A pure-rotation replay uses the same `apply_point` ops, so it is
    // bit-identical — any mismatch means a translation intervened.
    if pt3.coord != model.vertices.get(vh).point.as_array() {
        return Err(TipError::TranslateInterleaved);
    }
    Ok(pt3)
}

/// Assemble the `Pt3` from the definition (no guard — [`vertex_pt3`] guards the result).
fn build(model: &Model, vh: Handle<Vertex>) -> Result<Pt3, TipError> {
    let v = model.vertices.get(vh);
    match v.origin {
        Origin::Constructed => Ok(Pt3::at(coord_rat(v.point.as_array())?)),
        Origin::Discovered { .. } => Err(TipError::IndirectRequired),
        Origin::Rotated { base, rotation } => {
            // The root is a non-`Rotated` ancestor (1c invariant); build it, then replay
            // the rotation chain (root → leaf, via `parent` links).
            let mut p = build(model, base)?;
            for node in rotation_chain(model, rotation) {
                p = p.rotate_about(node.axis, node.angle, node.point);
            }
            Ok(p)
        }
    }
}

/// The rotation nodes from the root down to `leaf` (parent chain, reversed).
fn rotation_chain(model: &Model, leaf: Handle<Rotation>) -> Vec<RotNode> {
    let mut chain = Vec::new();
    let mut cur = Some(leaf);
    while let Some(h) = cur {
        let r: &Rotation = model.rotations.get(h);
        chain.push(RotNode {
            axis: r.axis,
            angle: r.angle,
            point: r.point,
        });
        cur = r.parent;
    }
    chain.reverse();
    chain
}

fn coord_rat(c: [f64; 3]) -> Result<[Rat; 3], TipError> {
    Ok([
        Rat::try_from_f64(c[0]).ok_or(TipError::Downgrade)?,
        Rat::try_from_f64(c[1]).ok_or(TipError::Downgrade)?,
        Rat::try_from_f64(c[2]).ok_or(TipError::Downgrade)?,
    ])
}

/// Sound `orient3d(a,b,c,d)` over four kernel vertices: `det[a−d, b−d, c−d]` — the
/// sign of a signed tetrahedron volume.
///
/// - **No `Discovered` vertex**: every point is direct. If all are exact (tol 0:
///   `Constructed`, unrotated or 90°-family), the exact Shewchuk predicate decides it;
///   if any carries rotation tol, the f64 filter → astro-float escalation judge does.
/// - **Exactly one `Discovered` vertex** (a three-plane seam, an implicit point): it is
///   moved to the `V` slot and its three planes are reconstructed from three rotated
///   points each ([`plane_pts`]), then [`indirect_orient3d_judge`] decides the sign
///   without materializing the seam's coordinate (stage 2c-ii).
/// - **Two or more `Discovered`** ([`TipError::IndirectRequired`]) — the ported
///   predicate carries a single implicit point; the two-implicit case is stage 3.
///
/// `Err` if any vertex is not recoverable here ([`TipError`]).
pub fn orient3d(
    model: &Model,
    a: Handle<Vertex>,
    b: Handle<Vertex>,
    c: Handle<Vertex>,
    d: Handle<Vertex>,
) -> Result<Orient, TipError> {
    let verts = [a, b, c, d];
    let discovered: Vec<usize> = verts
        .iter()
        .enumerate()
        .filter(|&(_, &vh)| matches!(model.vertices.get(vh).origin, Origin::Discovered { .. }))
        .map(|(i, _)| i)
        .collect();
    match discovered.as_slice() {
        [] => direct_orient3d(model, verts),
        [slot] => indirect_dispatch(model, verts, *slot),
        _ => Err(TipError::IndirectRequired), // two-plus implicit points (stage 3)
    }
}

/// The all-direct case of [`orient3d`]: build four `Pt3` and route by tol (tol-0 →
/// exact Shewchuk, tol > 0 → filter/escalation judge).
fn direct_orient3d(model: &Model, verts: [Handle<Vertex>; 4]) -> Result<Orient, TipError> {
    let p = [
        vertex_pt3(model, verts[0])?,
        vertex_pt3(model, verts[1])?,
        vertex_pt3(model, verts[2])?,
        vertex_pt3(model, verts[3])?,
    ];
    if p.iter().all(|q| q.tol == [0.0; 3]) {
        let det = nacre_predicates::orient3d(p[0].coord, p[1].coord, p[2].coord, p[3].coord);
        Ok(if det > 0.0 {
            Orient::Positive
        } else if det < 0.0 {
            Orient::Negative
        } else {
            Orient::Zero
        })
    } else {
        Ok(orient3d_judge(&p[0], &p[1], &p[2], &p[3]))
    }
}

/// `orient3d` with exactly one `Discovered` seam at `slot`. Move it to the `V` slot
/// (swap with slot 0; a single transposition, so flip the sign iff `slot != 0`),
/// reconstruct its three planes, build the other three points, and judge indirectly.
fn indirect_dispatch(
    model: &Model,
    verts: [Handle<Vertex>; 4],
    slot: usize,
) -> Result<Orient, TipError> {
    let mut v = verts;
    v.swap(0, slot);
    let flip = slot != 0;
    // v[0] is the seam (implicit V); v[1..4] are its explicit q, r, s.
    let surfs = seam_surfaces(model, v[0]).ok_or(TipError::IndirectRequired)?;
    let planes = [
        plane_pts(model, surfs[0])?,
        plane_pts(model, surfs[1])?,
        plane_pts(model, surfs[2])?,
    ];
    let q = vertex_pt3(model, v[1])?;
    let r = vertex_pt3(model, v[2])?;
    let s = vertex_pt3(model, v[3])?;
    let pl = |k: usize| (&planes[k][0], &planes[k][1], &planes[k][2]);
    let o = indirect_orient3d_judge(pl(0), pl(1), pl(2), &q, &r, &s);
    Ok(if flip { flip_orient(o) } else { o })
}

/// The three defining surfaces of a `Discovered` (`ThreePlane`) vertex, else `None`.
fn seam_surfaces(model: &Model, vh: Handle<Vertex>) -> Option<[Handle<Surface>; 3]> {
    match model.vertices.get(vh).origin {
        Origin::Discovered {
            definition: VertexDef::ThreePlane(surfs),
            ..
        } => Some(surfs),
        _ => None,
    }
}

/// Three rotated points defining the plane of surface `surf`, taken from the vertices of
/// any face on it. Only the exact provenance is needed, so the surface's own (rotated,
/// inexact) coefficients are never read — three of its face vertices carry the exact
/// definition. The whole `faces` store is scanned (not just live solids): an
/// append-only surface handle names one plane, and any face on it — even a superseded
/// one — has vertices that still lie on that plane. `Err` if fewer than three direct,
/// non-collinear vertices are found. The chosen triple's winding is irrelevant: the
/// indirect predicate is invariant to each plane's normal orientation.
fn plane_pts(model: &Model, surf: Handle<Surface>) -> Result<[Pt3; 3], TipError> {
    let mut seen = std::collections::HashSet::new();
    let mut cand: Vec<Pt3> = Vec::new();
    for (_, face) in model.faces.iter() {
        if face.surface != surf {
            continue;
        }
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                if let Some(bounds) = model.edges.get(he.edge).bounds {
                    for vh in bounds {
                        if seen.insert(vh) {
                            if let Ok(p) = vertex_pt3(model, vh) {
                                cand.push(p);
                            }
                        }
                    }
                }
            }
        }
    }
    select_three_noncollinear(&cand).ok_or(TipError::PlaneUnderdetermined)
}

/// Greedily pick three non-collinear points from `cand` (an f64 relative-area test — a
/// quality heuristic, not a soundness gate: an accidental collinear pick yields an exact
/// zero-normal plane, hence a declare-0, never a wrong sign).
fn select_three_noncollinear(cand: &[Pt3]) -> Option<[Pt3; 3]> {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let p0 = cand.first()?;
    // p1: distinct from p0.
    let p1 = cand.iter().find(|p| {
        let e = sub(p.coord, p0.coord);
        dot(e, e) > 0.0
    })?;
    // p2: not collinear with p0, p1 (relative area threshold).
    let e1 = sub(p1.coord, p0.coord);
    let p2 = cand.iter().find(|p| {
        let e2 = sub(p.coord, p0.coord);
        let n = cross(e1, e2);
        dot(n, n) > 1e-20 * dot(e1, e1) * dot(e2, e2)
    })?;
    Some([p0.clone(), p1.clone(), p2.clone()])
}

/// Flip an orientation for one argument transposition (`Zero` is fixed).
fn flip_orient(o: Orient) -> Orient {
    match o {
        Orient::Positive => Orient::Negative,
        Orient::Negative => Orient::Positive,
        Orient::Zero => Orient::Zero,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::Point3;
    use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
    use nacre_scalar::{Angle, Axis, Isometry, Rat as R, Rotation as SRot};

    fn rot30z() -> Isometry {
        Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [R::from_int(1), R::from_int(1), R::from_int(0)],
            angle: Angle::from_deg(R::from_int(30)).unwrap(),
        })
    }

    fn cuboid() -> (Model, Handle<nacre_topo::Solid>) {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        (m, s)
    }

    fn transformed(
        m: &mut Model,
        s: Handle<nacre_topo::Solid>,
        iso: &Isometry,
    ) -> Handle<nacre_topo::Solid> {
        let out = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: *iso,
            },
        )
        .unwrap();
        m.rebuild_adjacency();
        let OpOutput::Transform { solid } = out else {
            panic!("expected Transform");
        };
        solid
    }

    fn boundary_verts(m: &Model, s: Handle<nacre_topo::Solid>) -> Vec<Handle<Vertex>> {
        let mut seen = std::collections::HashSet::new();
        let mut vs = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            vs.push(vh);
                        }
                    }
                }
            }
        }
        vs
    }

    /// Four boundary vertices forming a non-degenerate tetrahedron (|orient3d| large).
    fn four_non_coplanar(m: &Model, vs: &[Handle<Vertex>]) -> [Handle<Vertex>; 4] {
        let c = |h: Handle<Vertex>| m.vertices.get(h).point.as_array();
        for i in 0..vs.len() {
            for j in i + 1..vs.len() {
                for k in j + 1..vs.len() {
                    for l in k + 1..vs.len() {
                        let d = nacre_predicates::orient3d(c(vs[i]), c(vs[j]), c(vs[k]), c(vs[l]));
                        if d.abs() > 1.0 {
                            return [vs[i], vs[j], vs[k], vs[l]];
                        }
                    }
                }
            }
        }
        panic!("no non-coplanar tetrahedron among {} vertices", vs.len());
    }

    /// A rotated Constructed vertex's `Pt3` reproduces the kernel coordinate bit-for-bit
    /// and carries a positive rotation tol.
    #[test]
    fn rotated_vertex_pt3_matches_and_has_tol() {
        let (mut m, s) = cuboid();
        let r = transformed(&mut m, s, &rot30z());
        for vh in boundary_verts(&m, r) {
            let p = vertex_pt3(&m, vh).unwrap();
            assert_eq!(
                p.coord,
                m.vertices.get(vh).point.as_array(),
                "coord bit-match"
            );
            assert!(p.tol.iter().any(|&t| t > 0.0), "rotation tol > 0");
        }
    }

    /// An unrotated Constructed vertex is exact (tol 0).
    #[test]
    fn unrotated_vertex_is_tol_zero() {
        let (m, s) = cuboid();
        let vh = boundary_verts(&m, s)[0];
        assert_eq!(vertex_pt3(&m, vh).unwrap().tol, [0.0; 3]);
    }

    /// Unrotated 4-vertex orient3d takes the exact Shewchuk path (matches
    /// `nacre-predicates` directly).
    #[test]
    fn orient3d_unrotated_matches_shewchuk() {
        let (m, s) = cuboid();
        let vs = boundary_verts(&m, s);
        let [a, b, c, d] = four_non_coplanar(&m, &vs);
        let cc = |h: Handle<Vertex>| m.vertices.get(h).point.as_array();
        let det = nacre_predicates::orient3d(cc(a), cc(b), cc(c), cc(d));
        let want = if det > 0.0 {
            Orient::Positive
        } else {
            Orient::Negative
        };
        assert_eq!(orient3d(&m, a, b, c, d).unwrap(), want);
    }

    /// Rotated 4-vertex orient3d takes the judge path and, on a definite (clearly
    /// non-coplanar) config, agrees with the exact Shewchuk sign of the same coords.
    #[test]
    fn orient3d_rotated_agrees_on_definite() {
        let (mut m, s) = cuboid();
        let r = transformed(&mut m, s, &rot30z());
        let vs = boundary_verts(&m, r);
        let [a, b, c, d] = four_non_coplanar(&m, &vs);
        let cc = |h: Handle<Vertex>| m.vertices.get(h).point.as_array();
        let det = nacre_predicates::orient3d(cc(a), cc(b), cc(c), cc(d));
        let want = if det > 0.0 {
            Orient::Positive
        } else {
            Orient::Negative
        };
        // routes to the judge (tol > 0) and resolves the definite config identically.
        assert_eq!(orient3d(&m, a, b, c, d).unwrap(), want);
    }

    /// A `Discovered` boolean-seam vertex is deferred to the indirect predicate.
    #[test]
    fn discovered_is_indirect_required() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        let cut = boolean(&mut m, BoolKind::Cut, a, b).unwrap()[0];
        m.rebuild_adjacency();
        // some seam vertex is Discovered → IndirectRequired.
        let deferred = boundary_verts(&m, cut)
            .into_iter()
            .any(|vh| matches!(vertex_pt3(&m, vh), Err(TipError::IndirectRequired)));
        assert!(
            deferred,
            "a boolean seam vertex must defer to the indirect predicate"
        );
    }

    /// A vertex translated after being rotated is deferred (the forest omits the
    /// translation, so base + chain no longer reproduce the coordinate).
    #[test]
    fn translate_after_rotate_is_deferred() {
        let (mut m, s) = cuboid();
        let r = transformed(&mut m, s, &rot30z());
        let t = transformed(
            &mut m,
            r,
            &Isometry::translation([R::from_int(5), R::from_int(-3), R::from_int(2)]),
        );
        let vh = boundary_verts(&m, t)[0];
        assert!(matches!(
            vertex_pt3(&m, vh),
            Err(TipError::TranslateInterleaved)
        ));
    }

    // ---- indirect bridge (2c-ii) ----

    /// The surface handles of the faces of solid `s` whose outer loop contains `v0`
    /// (three, for a cuboid corner).
    fn surfaces_at_vertex(
        m: &Model,
        s: Handle<nacre_topo::Solid>,
        v0: Handle<Vertex>,
    ) -> Vec<Handle<Surface>> {
        let sh = m.solids.get(s).outer;
        let mut surfs = Vec::new();
        for &fh in &m.shells.get(sh).faces {
            let face = m.faces.get(fh);
            let on = face.outer.half_edges.iter().any(|he| {
                m.edges
                    .get(he.edge)
                    .bounds
                    .is_some_and(|bd| bd.contains(&v0))
            });
            if on {
                surfs.push(face.surface);
            }
        }
        surfs
    }

    /// A synthetic `Discovered` seam at `at`'s coordinate, defined by three surfaces —
    /// the rotated-seam a real boolean cannot yet build (it rejects non-90° inputs).
    fn synth_discovered(
        m: &mut Model,
        at: Handle<Vertex>,
        surfs: [Handle<Surface>; 3],
    ) -> Handle<Vertex> {
        let point = m.vertices.get(at).point;
        m.vertices.push(Vertex {
            point,
            origin: Origin::Discovered {
                tol: 1e-9,
                definition: VertexDef::ThreePlane(surfs),
            },
        })
    }

    /// A cuboid corner where three faces meet, and its three surfaces.
    fn corner_with_surfaces(
        m: &Model,
        s: Handle<nacre_topo::Solid>,
    ) -> (Handle<Vertex>, [Handle<Surface>; 3]) {
        for v in boundary_verts(m, s) {
            let surfs = surfaces_at_vertex(m, s, v);
            if surfs.len() == 3 {
                return (v, [surfs[0], surfs[1], surfs[2]]);
            }
        }
        panic!("no corner with three faces");
    }

    /// Core check: a rotated seam judged indirectly (from its three planes) agrees with
    /// the direct judge on the same point, at every definite tetrahedron — with the
    /// seam in slot 0 (no sign flip) **and** slot 1 (one transposition → flip).
    #[test]
    fn indirect_matches_direct_on_rotated_corner() {
        let (mut m, s) = cuboid();
        let rs = transformed(&mut m, s, &rot30z());
        let (v0, surfs) = corner_with_surfaces(&m, rs);
        let vd = synth_discovered(&mut m, v0, surfs);
        let others: Vec<_> = boundary_verts(&m, rs)
            .into_iter()
            .filter(|&v| v != v0)
            .collect();
        let mut checked = 0usize;
        for i in 0..others.len() {
            for j in i + 1..others.len() {
                for k in j + 1..others.len() {
                    let (q, r, t) = (others[i], others[j], others[k]);
                    // slot 0: no flip. Oracle = direct judge over the explicit corner.
                    let truth0 = orient3d(&m, v0, q, r, t).unwrap();
                    if truth0 == Orient::Zero {
                        continue;
                    }
                    assert_eq!(orient3d(&m, vd, q, r, t).unwrap(), truth0, "slot 0");
                    // slot 1: one transposition → the dispatch must flip to match.
                    let truth1 = orient3d(&m, q, v0, r, t).unwrap();
                    assert_eq!(orient3d(&m, q, vd, r, t).unwrap(), truth1, "slot 1");
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no definite tetrahedron in the corpus");
    }

    /// An unrotated seam judged by this bridge (frame3 indirect) agrees with the
    /// independent exact-plane predicate (`nacre_geom::three_plane_orient3d`, keyed on
    /// exact coefficients) — two independent indirect implementations cross-checked.
    #[test]
    fn indirect_unrotated_matches_exact_plane() {
        use nacre_geom::intersect::three_plane_orient3d;
        let (mut m, s) = cuboid();
        m.rebuild_adjacency();
        let (v0, surfs) = corner_with_surfaces(&m, s);
        let plane = |h: Handle<Surface>| match m.surfaces.get(h) {
            Surface::Plane(p) => *p,
            _ => panic!("cuboid surfaces are planes"),
        };
        let (pa, pb, pc) = (plane(surfs[0]), plane(surfs[1]), plane(surfs[2]));
        let vd = synth_discovered(&mut m, v0, surfs);
        let others: Vec<_> = boundary_verts(&m, s)
            .into_iter()
            .filter(|&v| v != v0)
            .collect();
        let mut checked = 0usize;
        for i in 0..others.len() {
            for j in i + 1..others.len() {
                for k in j + 1..others.len() {
                    let c = |h: Handle<Vertex>| m.vertices.get(h).point;
                    let (q, r, t) = (c(others[i]), c(others[j]), c(others[k]));
                    let want = match three_plane_orient3d(&pa, &pb, &pc, q, r, t) {
                        1 => Orient::Positive,
                        -1 => Orient::Negative,
                        _ => continue,
                    };
                    assert_eq!(
                        orient3d(&m, vd, others[i], others[j], others[k]).unwrap(),
                        want
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no definite tetrahedron in the corpus");
    }

    /// A seam whose defining surface has no face (a fresh, unused surface handle) cannot
    /// reconstruct its plane → `PlaneUnderdetermined`.
    #[test]
    fn underdetermined_plane_is_deferred() {
        use nacre_geom::{Plane, Surface as GSurface};
        let (mut m, s) = cuboid();
        let vs = boundary_verts(&m, s);
        let dummy = m.surfaces.push(GSurface::Plane(
            Plane::through_points(
                Point3::from_array([0.0, 0.0, 0.0]),
                Point3::from_array([1.0, 0.0, 0.0]),
                Point3::from_array([0.0, 1.0, 0.0]),
            )
            .unwrap(),
        ));
        let vd = synth_discovered(&mut m, vs[0], [dummy, dummy, dummy]);
        assert!(matches!(
            orient3d(&m, vd, vs[1], vs[2], vs[3]),
            Err(TipError::PlaneUnderdetermined)
        ));
    }

    /// Two implicit points among the four defer (the ported predicate carries one).
    #[test]
    fn two_implicit_points_are_deferred() {
        let (mut m, s) = cuboid();
        let (v0, surfs) = corner_with_surfaces(&m, s);
        let vs: Vec<_> = boundary_verts(&m, s)
            .into_iter()
            .filter(|&v| v != v0)
            .collect();
        let vd1 = synth_discovered(&mut m, v0, surfs);
        let vd2 = synth_discovered(&mut m, vs[0], surfs);
        assert!(matches!(
            orient3d(&m, vd1, vd2, vs[1], vs[2]),
            Err(TipError::IndirectRequired)
        ));
    }
}
