//! Rotated-vertex `Pt3` assembly: the adapter that turns a b-rep `Model` vertex into the
//! [`Pt3`] the toleranced predicates ([`nacre_cip::predicate`]) consume.
//!
//! It reads the vertex's root coordinate and rotation forest and replays them into an exact
//! definition + directional tol. A `Constructed` root gives an exact (tol-0) point; each
//! rotation node transports the tol (`|R|·old`), so a non-90° rotation yields tol > 0. A
//! `Discovered` seam (an implicit three-plane point), or a vertex translated after being
//! rotated, cannot be reproduced from the rotation forest and is a typed [`Pt3Error`] (honest
//! defer) — the caller ([`crate::planes`]) treats either as "cannot judge here".

use nacre_cip::{Pt3, RotNode};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Rotation, Vertex};

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
