//! Rotated-vertex `Pt3` assembly: the adapter that turns a b-rep `Model` vertex into the
//! [`Pt3`] the toleranced predicates ([`nacre_cip::predicate`]) consume.
//!
//! It reads the vertex's root coordinate and rotation forest and replays them into an exact
//! definition + directional tol. A `Constructed` root gives an exact (tol-0) point; each
//! rotation node transports the tol (`|R|·old`), so a non-90° rotation yields tol > 0. A
//! `Discovered` seam (an implicit three-plane point), or a vertex translated after being
//! rotated, cannot be reproduced from the rotation forest and is a typed [`Pt3Error`] (honest
//! defer) — the caller ([`crate::planes`]) treats either as "cannot judge here".

use nacre_cip::{Pt3, RotNode, orient3d_judge};
use nacre_geom::Surface;
use nacre_scalar::{Orient, Rat};
use nacre_store::Handle;
use nacre_topo::{Face, Model, Origin, Rotation, Vertex, VertexDef};

/// Why a vertex's `Pt3` could not be assembled here (honest defer).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pt3Error {
    /// A `Discovered` seam (or a point rooted on one): its exact position is a
    /// three-plane intersection, not recoverable from the rotation forest.
    IndirectRequired,
    /// The vertex was translated after being rotated. The rotation forest records
    /// rotations only, so base + chain no longer reproduce the kernel coordinate
    /// and the exact position is not recoverable here.
    TranslateInterleaved,
    /// An f64 coordinate did not fit an exact i128 rational (downgrade).
    Downgrade,
    /// A face's supporting plane could not be witnessed by three assemblable, non-collinear
    /// points — the plane is named only by seam (or rotated-seam) vertices, so no exact
    /// definition is recoverable here (the deep rotated-chain floor). Honest defer.
    PlaneUnderdetermined,
}

/// The [`Pt3`] (coordinate + sound directional tol) of a kernel vertex, assembled from its root
/// coordinate and rotation forest. `Err` when it is not a direct (Constructed-rooted,
/// pure-rotation) point — see [`Pt3Error`].
pub(crate) fn vertex_pt3(model: &Model, vh: Handle<Vertex>) -> Result<Pt3, Pt3Error> {
    let pt3 = build(model, vh)?;
    // Pure-rotation guard: the forest omits translations, so a rotate-then-translate
    // history leaves base + chain ≠ the kernel coord (and hp_coord would drop the
    // translation). A pure-rotation replay uses the same `apply_point` ops, so it is
    // bit-identical — any mismatch means a translation intervened.
    if pt3.coord != model.vertices.get(vh).point.as_array() {
        return Err(Pt3Error::TranslateInterleaved);
    }
    Ok(pt3)
}

/// Assemble the `Pt3` from the definition (no guard — [`vertex_pt3`] guards the result).
fn build(model: &Model, vh: Handle<Vertex>) -> Result<Pt3, Pt3Error> {
    let v = model.vertices.get(vh);
    match v.origin {
        Origin::Constructed => Ok(Pt3::at(coord_rat(v.point.as_array())?)),
        Origin::Discovered { .. } => Err(Pt3Error::IndirectRequired),
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

/// The coordinate a `Constructed`-rooted vertex would have with `base_point` as its root and
/// `leaf`'s chain applied — the *same* computation [`build`] performs, in the same order.
///
/// A producer of rotated vertices must store this, not an independently computed image of the
/// point: [`vertex_pt3`] guards on the replay matching the stored coordinate **bit for bit**, and
/// two float routes to the same real number do not agree in the last places. `Mirror` needs it
/// because it derives a vertex's definition by conjugating the chain, which is a different route
/// from reflecting the coordinate.
pub(crate) fn replay_chain_coord(
    model: &Model,
    base_point: [f64; 3],
    leaf: Handle<Rotation>,
) -> Result<[f64; 3], Pt3Error> {
    let mut p = Pt3::at(coord_rat(base_point)?);
    for node in rotation_chain(model, leaf) {
        p = p.rotate_about(node.axis, node.angle, node.point);
    }
    Ok(p.coord)
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

fn coord_rat(c: [f64; 3]) -> Result<[Rat; 3], Pt3Error> {
    Ok([
        Rat::try_from_f64(c[0]).ok_or(Pt3Error::Downgrade)?,
        Rat::try_from_f64(c[1]).ok_or(Pt3Error::Downgrade)?,
        Rat::try_from_f64(c[2]).ok_or(Pt3Error::Downgrade)?,
    ])
}

/// Three exact, non-collinear `Pt3` known to lie on `surf`, drawn from any face still on it in
/// the append-only arena. A boolean output face shares its operand's `Surface` handle, and an
/// operand face's own vertices are `Constructed` (assemblable), so a plane that any live-or-
/// superseded face names can be witnessed exactly even when the *querying* face's vertices are
/// all `Discovered` (a seam-only face). `Err(PlaneUnderdetermined)` if fewer than three
/// assemblable, non-collinear vertices are found — a plane witnessed only by seam / rotated-seam
/// points (the deep rotated-chain floor). The chosen triple's winding is irrelevant: the plane
/// witness is used winding-invariantly.
fn plane_pts(model: &Model, surf: Handle<Surface>) -> Result<[Pt3; 3], Pt3Error> {
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
    select_three_noncollinear(&cand).ok_or(Pt3Error::PlaneUnderdetermined)
}

/// Greedily pick three non-collinear points from `cand` — an f64 relative-area test, a quality
/// heuristic, not a soundness gate: an accidental collinear pick yields an exact zero-normal
/// plane, hence a declare-0 downstream, never a wrong sign.
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
    let p1 = cand.iter().find(|p| {
        let e = sub(p.coord, p0.coord);
        dot(e, e) > 0.0
    })?;
    let e1 = sub(p1.coord, p0.coord);
    let p2 = cand.iter().find(|p| {
        let e2 = sub(p.coord, p0.coord);
        let n = cross(e1, e2);
        dot(n, n) > 1e-20 * dot(e1, e1) * dot(e2, e2)
    })?;
    Some([p0.clone(), p1.clone(), p2.clone()])
}

/// Three exact `Pt3` on the supporting plane of a **rotated-result** face, for the plane witness
/// (`tri_pt3`) `collect_planes` needs. A boolean output face lies on an input plane, so a rotated
/// result face's plane is `R(π)` for an operand plane `π`; `π`'s exact witness comes from the
/// operand face still on it (append-only), rotated by the face's own chain `R`.
///
/// Two tiers, both exact:
///  1. the face's own assemblable (`Constructed`-rooted) vertices, if three are non-collinear —
///     a surviving operand wall, witnessed by `F`'s own points (keeps `tri_pt3 == tri` there);
///  2. else a seam-dominated face: `π` is the original surface common to every `Discovered`
///     corner's `ThreePlane` (each corner lies on `π`), witnessed via [`plane_pts`] on the
///     original operand face and rotated by `R`. A non-unique `π` is disambiguated exactly by
///     which candidate's rotated plane carries every tier-1 anchor point (no f64 tie-break).
pub(crate) fn face_plane_witness(model: &Model, face: &Face) -> Result<[Pt3; 3], Pt3Error> {
    let outer: Vec<Handle<Vertex>> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| crate::he_start(model, he))
        .collect();

    // Tier 1 — the face's own exact points (survivors of a rotated operand).
    let anchors: Vec<Pt3> = outer
        .iter()
        .filter_map(|&vh| vertex_pt3(model, vh).ok())
        .collect();
    if let Some(t) = select_three_noncollinear(&anchors) {
        return Ok(t);
    }

    // Tier 2 — provenance witness. The uniform rotation chain (any outer vertex names it).
    let chain = outer
        .iter()
        .find_map(|&vh| match model.vertices.get(vh).origin {
            Origin::Rotated { rotation, .. } => Some(rotation_chain(model, rotation)),
            _ => None,
        })
        .ok_or(Pt3Error::PlaneUnderdetermined)?;

    // Candidate original surfaces: the intersection of every `Discovered` base's `ThreePlane`.
    let mut candidates: Option<Vec<Handle<Surface>>> = None;
    for &vh in &outer {
        let Origin::Rotated { base, .. } = model.vertices.get(vh).origin else {
            continue;
        };
        let Origin::Discovered {
            definition: VertexDef::ThreePlane(surfs),
            ..
        } = model.vertices.get(base).origin
        else {
            continue;
        };
        candidates = Some(match candidates {
            None => surfs.to_vec(),
            Some(prev) => prev.into_iter().filter(|s| surfs.contains(s)).collect(),
        });
    }
    let candidates = candidates.ok_or(Pt3Error::PlaneUnderdetermined)?;

    // The exact witness of an original surface, rotated into the face's frame.
    let rotated_witness = |pi: Handle<Surface>| -> Result<[Pt3; 3], Pt3Error> {
        let pts = plane_pts(model, pi)?;
        Ok(pts.map(|p| {
            let mut q = p;
            for node in &chain {
                q = q.rotate_about(node.axis, node.angle, node.point);
            }
            q
        }))
    };

    match candidates.as_slice() {
        [] => Err(Pt3Error::PlaneUnderdetermined),
        [pi] => rotated_witness(*pi),
        many => {
            // π is the candidate whose rotated plane carries every anchor (exact `orient3d`).
            // With no anchors this cannot be decided here, so defer honestly.
            if anchors.is_empty() {
                return Err(Pt3Error::PlaneUnderdetermined);
            }
            // No plane table exists yet, so the precision comes from the points at hand by the
            // same derivation `plane_index_setup` uses — **once**, outside the loop. Deriving it
            // per candidate realizes every anchor again, and a long rotation history makes that
            // realization the expensive thing in the whole operation.
            let judge = crate::planes::judge_for_points(&anchors);
            for &pi in many {
                if let Ok(wit) = rotated_witness(pi) {
                    if anchors.iter().all(|a| {
                        // "On the plane" is any judgement that is not a definite side — a proved
                        // zero, or a coincidence within the limit. An exhausted or degenerate one
                        // lands here too, which is the same reading as before the outcomes were
                        // told apart; distinguishing them needs the reporting channel.
                        orient3d_judge(a, &wit[0], &wit[1], &wit[2], judge).orient() == Orient::Zero
                    }) {
                        return Ok(wit);
                    }
                }
            }
            Err(Pt3Error::PlaneUnderdetermined)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoolKind, OpOutput, Operation, apply, boolean};
    use nacre_math::Point3;
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

    /// A `Discovered` boolean-seam vertex is deferred (its position is not in the forest).
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
            .any(|vh| matches!(vertex_pt3(&m, vh), Err(Pt3Error::IndirectRequired)));
        assert!(
            deferred,
            "a boolean seam vertex must defer (its position is not in the rotation forest)"
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
            Err(Pt3Error::TranslateInterleaved)
        ));
    }
}
