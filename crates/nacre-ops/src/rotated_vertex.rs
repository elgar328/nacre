//! Rotation-forest replay: turning a recorded rotation history back into coordinates.
//!
//! A rotation node names an axis, a pivot and an angle; a chain of them is the exact definition
//! of everything the rotation moved. This module walks that chain — [`rotation_chain`] reads it
//! root-to-leaf, [`replay_chain_coord`] turns a point through it in the *same* order and float
//! operations the producer used, which is what lets a consumer check a coordinate against its
//! definition bit for bit.
//!
//! It used to also *hunt* for a rotated face's exact plane, working back through vertices and
//! their `Discovered` three-plane definitions, because a `Surface` said nothing about where it
//! came from. Surfaces carry `SurfaceDef` now, so the plane is simply read (`crate::planes`) and
//! the hunt is gone.

use nacre_cip::{MoveNode, Pt3};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Model, Motion, MotionNode};

/// Why a coordinate could not be lifted to an exact rational.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pt3Error {
    /// An f64 coordinate did not fit an exact i128 rational (downgrade).
    Downgrade,
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
    leaf: Handle<MotionNode>,
) -> Result<[f64; 3], Pt3Error> {
    Ok(replay(Pt3::at(coord_rat(base_point)?), &motion_chain(model, leaf)).coord)
}

/// `p` carried through `chain`, in the producer's own order and float operations — the one
/// definition of "replay" this crate has, so a coordinate and its definition cannot drift apart.
pub(crate) fn replay(p: Pt3, chain: &[MoveNode]) -> Pt3 {
    chain.iter().fold(p, |q, n| match *n {
        MoveNode::Rotate { axis, angle, point } => q.rotate_about(axis, angle, point),
        MoveNode::Translate { offset } => q.translate(offset),
    })
}

/// The motion nodes from the root down to `leaf` (parent chain, reversed).
pub(crate) fn motion_chain(model: &Model, leaf: Handle<MotionNode>) -> Vec<MoveNode> {
    let mut chain = Vec::new();
    let mut cur = Some(leaf);
    while let Some(h) = cur {
        let n: &MotionNode = model.motions.get(h);
        chain.push(match n.motion {
            Motion::Rotate { axis, point, angle } => MoveNode::Rotate { axis, angle, point },
            Motion::Translate { offset } => MoveNode::Translate { offset },
        });
        cur = n.parent;
    }
    chain.reverse();
    chain
}

pub(crate) fn coord_rat(c: [f64; 3]) -> Result<[Rat; 3], Pt3Error> {
    Ok([
        Rat::try_from_f64(c[0]).ok_or(Pt3Error::Downgrade)?,
        Rat::try_from_f64(c[1]).ok_or(Pt3Error::Downgrade)?,
        Rat::try_from_f64(c[2]).ok_or(Pt3Error::Downgrade)?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{OpOutput, Operation, apply};
    use nacre_math::Point3;
    use nacre_scalar::{Angle, Axis, Isometry, Rat as R, Rotation as SRot};
    use nacre_topo::Origin;

    fn rot30z() -> Isometry {
        Isometry::rotation(SRot {
            axis: Axis::Z,
            point: [R::from_int(1), R::from_int(1), R::from_int(0)],
            angle: Angle::from_deg(R::from_int(30)).unwrap(),
        })
    }

    fn transformed(
        m: &mut Model,
        s: Handle<nacre_topo::Solid>,
        iso: &Isometry,
    ) -> Handle<nacre_topo::Solid> {
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: *iso,
            },
        )
        .unwrap() else {
            panic!("expected Transform");
        };
        m.rebuild_adjacency();
        solid
    }

    /// **A replay reproduces the stored coordinate bit for bit.**
    ///
    /// This is the contract that lets a consumer check a coordinate against its definition by
    /// equality rather than by tolerance: the replay must perform the *same* float operations, in
    /// the same order, that the producer did. `mirror` depends on it — it derives the definition
    /// by conjugating the chain, a different route from reflecting the point, so it has to store
    /// the replayed value.
    #[test]
    fn replay_reproduces_the_stored_coordinate() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        m.rebuild_adjacency();
        let r = transformed(&mut m, s, &rot30z());
        let sh = m.solids.get(r).outer;
        let mut checked = 0;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                for &vh in m.edges.get(he.edge).bounds.iter().flatten() {
                    let Origin::Moved {
                        base,
                        motion: rotation,
                    } = m.vertices.get(vh).origin
                    else {
                        panic!("a rotated solid's vertices carry their rotation");
                    };
                    let replayed =
                        replay_chain_coord(&m, m.vertices.get(base).point.as_array(), rotation)
                            .unwrap();
                    assert_eq!(replayed, m.vertices.get(vh).point.as_array());
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "walked no vertices");
    }

    /// A chain is read root-to-leaf, so replaying it applies the rotations in the order they
    /// happened — the reverse would be a different motion whenever the axes differ.
    #[test]
    fn a_chain_reads_root_to_leaf() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        m.rebuild_adjacency();
        let r = transformed(&mut m, s, &rot30z());
        let r = transformed(
            &mut m,
            r,
            &Isometry::rotation(SRot {
                axis: Axis::X,
                point: [R::from_int(0); 3],
                angle: Angle::from_deg(R::from_int(45)).unwrap(),
            }),
        );
        let sh = m.solids.get(r).outer;
        let fh = m.shells.get(sh).faces[0];
        let vh = m
            .edges
            .get(m.faces.get(fh).outer.half_edges[0].edge)
            .bounds
            .unwrap()[0];
        let Origin::Moved {
            motion: rotation, ..
        } = m.vertices.get(vh).origin
        else {
            panic!("rotated");
        };
        let axes: Vec<Axis> = motion_chain(&m, rotation)
            .iter()
            .filter_map(|n| match n {
                MoveNode::Rotate { axis, .. } => Some(*axis),
                MoveNode::Translate { .. } => None,
            })
            .collect();
        assert_eq!(axes, vec![Axis::Z, Axis::X]);
    }
}
