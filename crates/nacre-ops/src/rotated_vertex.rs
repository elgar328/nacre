//! Motion-forest replay: turning a recorded motion history back into coordinates.
//!
//! A node names one motion — a turn about an axis, a translation, or a reflection in a coordinate
//! plane — and a chain of them is the exact definition of everything that motion moved. This
//! module walks that chain: [`motion_chain`] reads it root-to-leaf and [`replay`] carries a point
//! through it in the *same* order and float operations the producer used, which is what lets a
//! consumer check a coordinate against its definition bit for bit.
//!
//! It used to also *hunt* for a rotated face's exact plane, working back through vertices and
//! their `Discovered` three-plane definitions, because a `Surface` said nothing about where it
//! came from. Surfaces carry `SurfaceDef` now, so the plane is simply read (`crate::planes`) and
//! the hunt is gone.

use nacre_cip::{MoveNode, Pt3};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Model, Motion, MotionNode, SurfaceDef};

/// Why a coordinate could not be lifted to an exact rational.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pt3Error {
    /// An f64 coordinate did not fit an exact i128 rational (downgrade).
    Downgrade,
}

/// The coordinate a `Constructed`-rooted vertex would have with `base_point` as its root and
/// `leaf`'s chain applied — the *same* computation [`build`] performs, in the same order.
///
/// **The invariant, not a producer.** Every producer reaches its coordinate by walking the motion
/// it just recorded, so no caller needs this to *build* anything — what needs it is the check that
/// the two agree **bit for bit**, which the predicates rely on and which two float routes to the
/// same real number would not satisfy. Reflection used to be the exception (it derived a
/// definition by conjugating the chain, a different route from reflecting the coordinate); it is
/// a chain node now, so the exception is gone and only the assertions remain.
#[cfg(test)]
pub(crate) fn replay_chain_coord(
    model: &Model,
    base_point: [f64; 3],
    leaf: Handle<MotionNode>,
) -> Result<[f64; 3], Pt3Error> {
    let chain = motion_chain(model, leaf).ok_or(Pt3Error::Downgrade)?;
    Ok(replay(Pt3::at(coord_rat(base_point)?), &chain)
        .ok_or(Pt3Error::Downgrade)?
        .coord)
}

/// `p` carried through `chain`, in the producer's own order and float operations — the one
/// definition of "replay" this crate has, so a coordinate and its definition cannot drift apart.
///
/// `None` only for a [`MoveNode::Frame`] whose squared lengths do not fit `i128`, which
/// [`motion_chain`] already refuses to emit — so in practice this is infallible, and the `Option`
/// is here so that "in practice" does not have to be an invariant spanning two crates.
pub(crate) fn replay(p: Pt3, chain: &[MoveNode]) -> Option<Pt3> {
    chain.iter().try_fold(p, |q, n| match *n {
        MoveNode::Rotate { axis, angle, point } => Some(q.rotate_about(axis, angle, point)),
        MoveNode::Translate { offset } => Some(q.translate(offset)),
        MoveNode::Mirror { axis, offset } => Some(q.mirror(axis, offset)),
        MoveNode::Frame { frame } => q.frame(frame),
    })
}

/// The motion nodes from the root down to `leaf` (parent chain, reversed).
///
/// ★★★ **A frame expands into more than one node, and that is the recursion.** `Motion::Frame`
/// names a *plane*, not a basis, because a wall raised on a tilted face has no rational world
/// normal to spell. Reading it means: take that plane's own rational coefficients, derive the
/// frame from them ([`nacre_scalar::plane_frame`] — one spelling, so the exact route and the f64
/// `frame_axes` cannot drift), and then keep going through **that plane's** motion, which is what
/// carries the result out of its frame and into the next one down. The walk terminates at a plane
/// with no frame, which is the world.
///
/// `None` when a frame cannot be built exactly — no coefficients recorded, a degenerate plane, or
/// squared lengths past `i128`. That is a decline, not a reject: the caller falls back to the f64
/// path it was on before frames existed.
pub(crate) fn motion_chain(model: &Model, leaf: Handle<MotionNode>) -> Option<Vec<MoveNode>> {
    let mut chain = Vec::new();
    let mut cur = Some(leaf);
    while let Some(h) = cur {
        let n: &MotionNode = model.motions.get(h);
        match n.motion {
            Motion::Rotate { axis, point, angle } => {
                chain.push(MoveNode::Rotate { axis, angle, point })
            }
            Motion::Translate { offset } => chain.push(MoveNode::Translate { offset }),
            Motion::Mirror { axis, offset } => chain.push(MoveNode::Mirror { axis, offset }),
            Motion::Frame {
                plane,
                origin,
                ref_dir,
                flip,
            } => {
                // Built root-to-leaf here and reversed at the end, so it goes on backwards.
                let mut c = frame_chain(model, plane, origin, ref_dir, flip)?;
                c.reverse();
                chain.append(&mut c);
            }
        }
        cur = n.parent;
    }
    chain.reverse();
    Some(chain)
}

/// **One plane's frame, as the motion nodes that carry it out to the world** — in reading order,
/// root first.
///
/// ★★★ **The frame node and the plane's own motion are one chain, not two.** A plane states
/// itself in *its* frame, so whatever moved that plane has to run after the frame — and the
/// recursion ends at a `Constructed` plane, which states itself in the world.
///
/// The single spelling of that, so [`motion_chain`] and the sketch frame `nacre-ops` reports to a
/// caller cannot describe different frames — which they must not, because `face_plane`'s contract
/// is that it names the frame `PadOnFace` actually places a profile in.
pub(crate) fn frame_chain(
    model: &Model,
    plane: Handle<nacre_geom::Surface>,
    origin: [Rat; 3],
    ref_dir: [Rat; 3],
    flip: bool,
) -> Option<Vec<MoveNode>> {
    // ★★★★ **Canonical coefficients carry no direction, and a frame needs one.**
    // `canonical_plane_coeffs` forces the first nonzero component positive, because its question
    // is *"are these the same plane"* — where direction is noise. A frame's `ŵ` **is** a
    // direction, so the node carries the sense in `flip` and this is where it is spent.
    let c = *model.surface_coeffs.get(&plane)?;
    let zero = Rat::from_int(0);
    let c = if flip {
        let mut neg = [zero; 4];
        for (k, x) in neg.iter_mut().enumerate() {
            *x = zero.checked_sub(c[k])?;
        }
        neg
    } else {
        c
    };
    let mut chain = vec![MoveNode::Frame {
        frame: nacre_scalar::plane_frame_named(c, origin, ref_dir)?,
    }];
    match model.surface_defs.get(&plane) {
        Some(SurfaceDef::Moved { motion, .. }) => chain.append(&mut motion_chain(model, *motion)?),
        Some(SurfaceDef::Constructed) => {}
        // `Inexact` has no exact definition to continue with, and an unrecorded surface has none
        // at all — both are declines, not rejects.
        _ => return None,
    }
    Some(chain)
}

/// A frame's origin and its three axes, realized in the world: `(origin, û, v̂, ŵ)`.
pub(crate) type WorldBasis = ([f64; 3], [f64; 3], [f64; 3], [f64; 3]);

/// **The world image of a frame's origin and axes** — what a caller sees as the sketch plane, and
/// what the operation places its profile in. `(origin, û, v̂, ŵ)`.
///
/// ★★★★ **Realized, because the sign and the direction cannot be reasoned out from the
/// coefficients.** They are canonical (no direction) *and* written in the plane's pre-motion
/// frame, so a dot product against a world normal compares two different frames — the mistake
/// that produced `PadMissesFace` on a twice-turned fixture, and one that reads as perfectly
/// plausible right up until the plane has a motion. Replaying the axes and looking at where they
/// land asks the question that is actually being asked.
///
/// The axes are differences of replayed points, so the origin and every translation cancel and
/// only the linear part is left.
pub(crate) fn frame_world_basis(
    model: &Model,
    plane: Handle<nacre_geom::Surface>,
    origin: [Rat; 3],
    ref_dir: [Rat; 3],
    flip: bool,
) -> Option<WorldBasis> {
    let chain = frame_chain(model, plane, origin, ref_dir, flip)?;
    let at = |p: [i128; 3]| -> Option<[f64; 3]> {
        Some(replay(Pt3::at(p.map(Rat::from_int)), &chain)?.coord)
    };
    let o = at([0, 0, 0])?;
    let axis = |p: [i128; 3]| -> Option<[f64; 3]> {
        let q = at(p)?;
        Some([q[0] - o[0], q[1] - o[1], q[2] - o[2]])
    };
    Some((o, axis([1, 0, 0])?, axis([0, 1, 0])?, axis([0, 0, 1])?))
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

    /// One step of a motion chain — the two operations that leave a `MotionNode` behind.
    enum Step {
        // Boxed: an `Isometry` dwarfs a `(Axis, Rat)`, and clippy is right that the unboxed
        // variant would make every `Flip` in the array pay for it.
        Move(Box<Isometry>),
        Flip(Axis, R),
    }

    fn reflected(
        m: &mut Model,
        s: Handle<nacre_topo::Solid>,
        axis: Axis,
        offset: R,
    ) -> Handle<nacre_topo::Solid> {
        let OpOutput::Mirror { solid } = apply(
            m,
            &Operation::Mirror {
                solid: s,
                axis,
                offset,
            },
        )
        .unwrap() else {
            panic!("expected Mirror");
        };
        m.rebuild_adjacency();
        solid
    }

    /// A turn about `axis` through an integer pivot, by `n/d` degrees.
    fn turn(axis: Axis, pivot: [i128; 3], (n, d): (i128, i128)) -> Isometry {
        Isometry::rotation(SRot {
            axis,
            point: pivot.map(R::from_int),
            angle: Angle::from_deg(R::new(n, d).unwrap()).unwrap(),
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
    /// the same order, that the producer did. Every motion node owes this, and `Mirror` is where
    /// it is easiest to lose — `Pt3::mirror` walks the same `2c − x` that `AxisMirror::point` did,
    /// rather than an algebraically equal rearrangement.
    ///
    /// ★★★ **It is also what stands between this kernel and a realization that is not a function
    /// of its angle.** `Angle`'s f64 route is `(deg.to_f64() * PI / 180.0).cos()`, and that is
    /// measured to give two answers one ulp apart — between a debug and a release build, and
    /// between two call sites within one release build, where LLVM evaluates a literal-angle site
    /// at compile time and leaves the other to libm. Both routes below realize the same angle, so a
    /// producer and a consumer that folded differently would land different coordinates here.
    ///
    /// ★★ **What holds it is that the angle crosses the model store.** It is written into an
    /// `Operation::Transform`, pushed, and read back out before either route realizes it — and no
    /// optimiser propagates a constant through a heap structure. That is *why* this passes with a
    /// literal angle, and it is worth knowing: a future route that realizes an `Angle` it never
    /// stored is not covered by this argument, only by this test happening to exercise it.
    ///
    /// **Measured, not argued:** the whole census is bit-identical between a debug and a release
    /// build (`tests/census.rs` documents the diff), which is the direct reading of the same claim
    /// over 130 cases rather than one.
    #[test]
    fn replay_reproduces_the_stored_coordinate() {
        // ★ A chain, not one turn, and none of it "nice": a pivot off the origin, an angle whose
        // realization is nowhere near a quadrantal one, then a translate and a reflection. A single
        // Z-turn about a rational pivot exercises one node and one of the two rotate axes.
        let mv = |iso| Step::Move(Box::new(iso));
        let chains: [Vec<Step>; 3] = [
            vec![mv(rot30z())],
            vec![
                mv(turn(Axis::Z, [1, 1, 0], (2749, 71))),
                mv(turn(Axis::X, [0, 3, 2], (617, 9))),
            ],
            vec![
                mv(turn(Axis::Y, [5, 0, 1], (89999, 1000))),
                mv(Isometry::translation([
                    R::new(7, 3).unwrap(),
                    R::from_int(-2),
                    R::from_int(0),
                ])),
                // ★ A reflection between two turns, because this is the node the contract is
                // easiest to lose on and the one whose parity the chain has to carry.
                Step::Flip(Axis::X, R::new(1, 2).unwrap()),
                mv(turn(Axis::Z, [0, 0, 0], (271, 4))),
            ],
        ];
        let mut checked = 0;
        for chain in &chains {
            let mut m = Model::new();
            let s = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            m.rebuild_adjacency();
            let mut r = s;
            for step in chain {
                r = match step {
                    Step::Move(iso) => transformed(&mut m, r, iso),
                    Step::Flip(axis, offset) => reflected(&mut m, r, *axis, *offset),
                };
            }
            let sh = m.solids.get(r).outer;
            for &fh in &m.shells.get(sh).faces {
                for he in &m.faces.get(fh).outer.half_edges {
                    for &vh in m.edges.get(he.edge).bounds.iter().flatten() {
                        let Origin::Moved {
                            base,
                            motion: rotation,
                        } = m.vertices.get(vh).origin
                        else {
                            panic!("a moved solid's vertices carry their motion");
                        };
                        let replayed =
                            replay_chain_coord(&m, m.vertices.get(base).point.as_array(), rotation)
                                .unwrap();
                        assert_eq!(replayed, m.vertices.get(vh).point.as_array());
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 100, "walked only {checked} vertices");
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
            .expect("an axis-aligned history holds no frame")
            .iter()
            .filter_map(|n| match n {
                MoveNode::Rotate { axis, .. } => Some(*axis),
                MoveNode::Translate { .. } | MoveNode::Mirror { .. } | MoveNode::Frame { .. } => {
                    None
                }
            })
            .collect();
        assert_eq!(axes, vec![Axis::Z, Axis::X]);
    }
}
