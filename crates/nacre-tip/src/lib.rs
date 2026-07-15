//! Toleranced indirect predicates over the b-rep `Model` (overhaul §TIP, stage 2b).
//!
//! A read-only analysis above `nacre-topo` (like `nacre-validate` / `nacre-props`): it
//! assembles a rotated vertex's exact **definition** + **directional tol** from the
//! rotation forest ([`nacre_scalar::frame3::Pt3`]) and judges `orient3d` over kernel
//! vertices. tol-0 (`Constructed`, unrotated or 90°-family) configs take the exact
//! Shewchuk path (`nacre-predicates`); any tol > 0 (rotated) takes the filter →
//! astro-float escalation ([`nacre_scalar::frame3::orient3d_judge`]).
//!
//! **Scope (stage 2b): direct predicate over Constructed-rooted, pure-rotation
//! points.** `Discovered` seams need the *indirect* predicate (their exact position is
//! a three-plane intersection, not their f64 cache) — deferred to stage 2c. A vertex
//! translated after being rotated is deferred too (the forest records rotations only,
//! §⑦, so its position is not recoverable here). Both surface as a typed [`TipError`]
//! (honest defer) rather than a silent wrong sign. Not yet wired into boolean (stage 3).

use nacre_scalar::frame3::{Pt3, RotNode, orient3d_judge};
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Rotation, Vertex};

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
/// sign of a signed tetrahedron volume. If every point is exact (tol 0: `Constructed`,
/// unrotated or 90°-family), the exact Shewchuk predicate decides it; if any point
/// carries rotation tol, the f64 filter → astro-float escalation judge does. `Err` if
/// any vertex is not a direct point ([`TipError`]).
pub fn orient3d(
    model: &Model,
    a: Handle<Vertex>,
    b: Handle<Vertex>,
    c: Handle<Vertex>,
    d: Handle<Vertex>,
) -> Result<Orient, TipError> {
    let p = [
        vertex_pt3(model, a)?,
        vertex_pt3(model, b)?,
        vertex_pt3(model, c)?,
        vertex_pt3(model, d)?,
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
}
