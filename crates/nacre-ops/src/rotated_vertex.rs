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
//! came from. Surfaces carry their truth now (points + motion), so the plane is simply read
//! (`crate::planes`) and
//! the hunt is gone.

use nacre_cip::{MoveNode, WitnessPoint};
use nacre_scalar::Rat;
use nacre_store::Handle;
use nacre_topo::{Model, Motion, MotionNode, Surface};

/// Why a coordinate could not be lifted to an exact rational.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum WitnessPointError {
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
) -> Result<[f64; 3], WitnessPointError> {
    let chain = motion_chain(model, leaf).ok_or(WitnessPointError::Downgrade)?;
    Ok(replay(WitnessPoint::at(coord_rat(base_point)?), &chain)
        .ok_or(WitnessPointError::Downgrade)?
        .coord)
}

/// `p` carried through `chain`, in the producer's own order and float operations — the one
/// definition of "replay" this crate has, so a coordinate and its definition cannot drift apart.
///
/// `None` only for a [`MoveNode::Frame`] whose squared lengths do not fit `i128`, which
/// [`motion_chain`] already refuses to emit — so in practice this is infallible, and the `Option`
/// is here so that "in practice" does not have to be an invariant spanning two crates.
pub(crate) fn replay(p: WitnessPoint, chain: &[MoveNode]) -> Option<WitnessPoint> {
    chain.iter().try_fold(p, |q, n| match n {
        MoveNode::Rotate { axis, angle, point } => Some(q.rotate_about(*axis, *angle, *point)),
        MoveNode::Translate { offset } => Some(q.translate(*offset)),
        MoveNode::Mirror { axis, offset } => Some(q.mirror(*axis, *offset)),
        MoveNode::Frame { frame } => q.frame(*frame),
        MoveNode::FrameWide(f) => q.frame_wide(f),
        MoveNode::FrameThrough(f) => q.frame_through(f),
    })
}

/// **A surface's plane as its exact witness triangle** — the per-carrier building block of an
/// implicit point (16-2). The same computation `collect_planes`' arms perform per face, spelled
/// once for the per-carrier consumer: `Known` points replay through the surface's own motion; a
/// **named** `Through` plane (rational closure) solves through `through_points_rat` and replays
/// the same way. ★ A named `Through` plane whose meets are wider than any witness base (open
/// item 17) gets its **own frame's canonical probes** instead — on-plane exact definitions, the
/// same witness `collect_planes`' probe branch builds; `frame_chain` already carries the
/// plane's later motion, so those return directly. `None` for a nameless `Through` carrier
/// (a datum on a nameless datum — depth, open item 16-3's question) and for a cylinder.
pub(crate) fn surface_witness_triangle(
    model: &Model,
    h: Handle<Surface>,
) -> Option<[WitnessPoint; 3]> {
    let (base, motion) = match model.surface_truth(h) {
        nacre_topo::Surface::Plane {
            points: nacre_topo::PlanePoints::Known(pts),
            motion,
        } => (*pts, *motion),
        nacre_topo::Surface::Plane {
            points: nacre_topo::PlanePoints::Through(vs),
            motion,
        } => match model.through_points_rat(*vs) {
            Some(base) => (base, *motion),
            None if model.surface_name.contains_key(&h) => {
                let chain = frame_chain(model, h, &nacre_topo::FramePlacement::Canonical, false)?;
                let r = nacre_scalar::Rat::from_int;
                let probe = |u: i128, v: i128| replay(WitnessPoint::at([r(u), r(v), r(0)]), &chain);
                return Some([probe(0, 0)?, probe(1, 0)?, probe(0, 1)?]);
            }
            None => return None, // nameless Through — depth
        },
        nacre_topo::Surface::Cylinder { .. } => return None,
    };
    let w = base.map(WitnessPoint::at);
    match motion {
        None => Some(w),
        Some(m) => {
            let chain = motion_chain(model, m)?;
            let mut out = w;
            for o in out.iter_mut() {
                *o = replay(o.clone(), &chain)?;
            }
            Some(out)
        }
    }
}

/// The three defining vertices of a **nameless** `Through` plane as judged points — the judged
/// road's twin of `Model::through_points_rat` (the named road's single solve). One spelling,
/// shared by the datum producer's pre-push validation and [`frame_chain`]'s node assembly, so
/// what was validated and what gets framed cannot drift.
///
/// Per vertex: **pure** (its three carriers share one motion, and the solve fits `Rat`) becomes
/// [`nacre_cip::JudgedPoint::Pure`] — 16-1's whole population; anything else that is still a
/// well-defined
/// three-plane meet becomes [`nacre_cip::JudgedPoint::Meet`] of its carriers' witness triangles — the
/// straddling population (16-2), **and** the pure-but-too-wide one, which a meet represents
/// without ever asking the coordinate to fit anything.
///
/// `None` for a seam vertex, or a carrier that has no witness triangle of its own (a nameless
/// `Through` carrier — depth — or a chain outside the decimal window). The producer tells those
/// apart by cause before any plane is pushed; here one answer suffices, because reaching this
/// through an *accepted* statement re-derives the same deterministic result.
pub(crate) fn through_judged_points(
    model: &Model,
    vs: [Handle<nacre_topo::Vertex>; 3],
) -> Option<[nacre_cip::JudgedPoint; 3]> {
    use nacre_cip::JudgedPoint;
    let mut out: [Option<JudgedPoint>; 3] = [None, None, None];
    for (o, vh) in out.iter_mut().zip(vs) {
        let tri = match model.vertices.get(vh).def {
            nacre_topo::VertexDef::ThreePlane(tri) => tri,
            // OnSeam pins a curve, not a point; a Branch point's coordinates are
            // quadratic-irrational, and this table's witnesses are rational by type —
            // both decline, per variant (M6-2's judging of branch points is new machinery,
            // not this road).
            nacre_topo::VertexDef::OnSeam(_) | nacre_topo::VertexDef::Branch { .. } => {
                return None;
            }
        };
        // ★★ **The frame question is the door's** ([`nacre_topo::Model::vertex_meet`]). This was
        // the last of the four places that compared `plane_motion(tri[0])` against the other two
        // itself, and that reading calls a turned solid's corner a straddle: the invariant-plane
        // restatement leaves its cap world-stated while the walls carry a node, and only the
        // chain-fixes licence tells that apart from a real straddle. Here the copy cost no
        // *capability* — a corner it misread still gets a correct `Meet` — but it made the judged
        // road hand out a weaker witness than it had to.
        //
        // ★ The old spelling also let one vertex's non-meeting `?` out of the whole call. A fact
        // about one vertex is not a fact about the triple; a vertex without a rational meet falls
        // to `Meet`, which is exactly what that variant is for.
        let pure = model.vertex_meet(vh).and_then(|(meet, frame)| {
            meet.narrow()
                .map(|p| WitnessPoint::at(*p))
                .and_then(|wp| match frame {
                    None => Some(wp),
                    Some(m) => replay(wp, &motion_chain(model, m)?),
                })
        });
        *o = Some(match pure {
            Some(wp) => JudgedPoint::Pure(wp),
            None => JudgedPoint::Meet(Box::new([
                surface_witness_triangle(model, tri[0])?,
                surface_witness_triangle(model, tri[1])?,
                surface_witness_triangle(model, tri[2])?,
            ])),
        });
    }
    Some([out[0].take()?, out[1].take()?, out[2].take()?])
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
/// `None` when a frame cannot be built exactly — no name recorded at all, or a degenerate plane.
/// ★ Width is **not** on that list since S4: a `Wide` name or overflowing squared lengths take the
/// arbitrary-precision road (`MoveNode::FrameWide`) instead of declining.
///
/// ★★ **This used to say the caller "falls back to the f64 path it was on before frames existed".
/// That is stale, and it misled a later plan into inventing a precondition.** No such fallback
/// remains: `planes.rs` turns this `None` into `RejectReason::FrameOutOfRange`, `exact.rs` into
/// `OpError::PlaneWithoutExactForm` (S6b deleted the f64 prism road), `reuse.rs` declines to the
/// arrangement — an equally correct road, not a degraded one — and `replay_chain_coord` names
/// `WitnessPointError::Downgrade`. Every consumer is honest; and since S2 the branch is unreachable anyway,
/// because every plane has a name.
pub(crate) fn motion_chain(model: &Model, leaf: Handle<MotionNode>) -> Option<Vec<MoveNode>> {
    let mut chain = Vec::new();
    let mut cur = Some(leaf);
    while let Some(h) = cur {
        let n: &MotionNode = model.motion(h);
        match n.motion {
            Motion::Rotate { axis, point, angle } => {
                chain.push(MoveNode::Rotate { axis, angle, point })
            }
            Motion::Translate { offset } => chain.push(MoveNode::Translate { offset }),
            Motion::Mirror { axis, offset } => chain.push(MoveNode::Mirror { axis, offset }),
            Motion::Frame {
                plane,
                placement,
                flip,
            } => {
                // Built root-to-leaf here and reversed at the end, so it goes on backwards.
                let mut c = frame_chain(model, plane, &placement, flip)?;
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
    plane: Handle<Surface>,
    placement: &nacre_topo::FramePlacement,
    flip: bool,
) -> Option<Vec<MoveNode>> {
    use nacre_topo::FramePlacement;
    let Some(name) = model.surface_name.get(&plane) else {
        // ★★★ **The judged road** (open item 16, first wall): a plane with no name at all — a
        // mixed-frame `Through` statement, whose exact world coefficients are irrational. Its
        // frame is derived from the defining points as intervals ([`nacre_cip::FrameThrough`]),
        // with the branch decided once at the fixed rung, so the same statement always frames
        // the same way. `Named` placement is refused here defensively (`SketchFrame::named`
        // already rejects it by name at construction — an on-plane claim needs a name to verify
        // against), and the plane's own later motion is appended by the shared tail below,
        // exactly as for the named roads.
        let nacre_topo::Surface::Plane {
            points: nacre_topo::PlanePoints::Through(vs),
            motion,
        } = model.surface_truth(plane)
        else {
            return None;
        };
        if !matches!(placement, FramePlacement::Canonical) {
            return None;
        }
        let pts = through_judged_points(model, *vs)?;
        let mut chain = vec![MoveNode::FrameThrough(Box::new(
            nacre_cip::FrameThrough::of(pts, flip)?,
        ))];
        if let Some(m) = motion {
            chain.append(&mut motion_chain(model, *m)?);
        }
        return Some(chain);
    };

    // ★★★★ **Canonical coefficients carry no direction, and a frame needs one.**
    // `canonical_plane_coeffs` forces the first nonzero component positive, because its question
    // is *"are these the same plane"* — where direction is noise. A frame's `ŵ` **is** a
    // direction, so the node carries the sense in `flip`, and both roads below spend it on the
    // coefficients just before building the frame.
    //
    // **The narrow road** — bit for bit the pre-S4 frame. For `Canonical` the placement pair is
    // derived from the canonical (unflipped) coefficients first and the sign is applied after —
    // the same split `face_frame`/`frame_chain` had before the derivation moved here, and a
    // correctness condition: `ref_dir = ẑ × n` is sign-sensitive.
    let narrow_road = || -> Option<MoveNode> {
        let c = *name.narrow()?;
        let (origin, ref_dir) = match placement {
            FramePlacement::Named { origin, ref_dir } => (*origin, *ref_dir),
            FramePlacement::Canonical => nacre_scalar::plane_frame_default(c)?,
        };
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
        Some(MoveNode::Frame {
            frame: nacre_scalar::plane_frame_named(c, origin, ref_dir)?,
        })
    };
    // **The wide road** (S4) — the same convention through arbitrary precision, where nothing
    // can overflow. Taken when the name is `Wide` or when any narrow derivation step hits
    // `i128` (the measured 1.6% `n·n` population). This is the frame arm of the f64-fallback
    // chain closing: a plane with a name can now always host a sketch.
    let wide_road = || -> Option<MoveNode> {
        let wf = match placement {
            FramePlacement::Canonical => nacre_cip::WideFrame::canonical_of(name, flip)?,
            FramePlacement::Named { origin, ref_dir } => {
                nacre_cip::WideFrame::named_of(name, origin, ref_dir, flip)?
            }
        };
        Some(MoveNode::FrameWide(wf))
    };
    let mut chain = vec![narrow_road().or_else(wide_road)?];
    match model.surface_truth(plane) {
        nacre_topo::Surface::Plane {
            motion: Some(m), ..
        }
        | nacre_topo::Surface::Cylinder {
            motion: Some(m), ..
        } => chain.append(&mut motion_chain(model, *m)?),
        nacre_topo::Surface::Plane { motion: None, .. }
        | nacre_topo::Surface::Cylinder { motion: None, .. } => {}
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
    plane: Handle<Surface>,
    placement: &nacre_topo::FramePlacement,
    flip: bool,
) -> Option<WorldBasis> {
    let chain = frame_chain(model, plane, placement, flip)?;
    let at = |p: [i128; 3]| -> Option<[f64; 3]> {
        Some(replay(WitnessPoint::at(p.map(Rat::from_int)), &chain)?.coord)
    };
    let o = at([0, 0, 0])?;
    let axis = |p: [i128; 3]| -> Option<[f64; 3]> {
        let q = at(p)?;
        Some([q[0] - o[0], q[1] - o[1], q[2] - o[2]])
    };
    Some((o, axis([1, 0, 0])?, axis([0, 1, 0])?, axis([0, 0, 1])?))
}

pub(crate) fn coord_rat(c: [f64; 3]) -> Result<[Rat; 3], WitnessPointError> {
    Ok([
        Rat::try_from_f64(c[0]).ok_or(WitnessPointError::Downgrade)?,
        Rat::try_from_f64(c[1]).ok_or(WitnessPointError::Downgrade)?,
        Rat::try_from_f64(c[2]).ok_or(WitnessPointError::Downgrade)?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{OpOutput, Operation, apply};
    use nacre_math::Point3;
    use nacre_scalar::{Angle, Axis, Isometry, Rat as R, Rotation as SRot};

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
    /// it is easiest to lose — `WitnessPoint::mirror` walks the same `2c − x` that `AxisMirror::point` did,
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
        let chains: [Vec<Step>; 4] = [
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
            // ★ A same-axis depth-2 chain: the caps are fixed by *both* turns (restated
            // twice, still world-stated), so their corners exercise the fixed-carrier
            // licence at depth 2 — walls [Z, Z], caps None, `chain_fixes_plane` proving
            // the caps ride the whole chain.
            vec![mv(rot30z()), mv(turn(Axis::Z, [2, -1, 3], (2749, 71)))],
        ];
        let mut checked = 0;
        let mut declined = 0;
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
                    for &vh in m.edges.get(he.edge).vertices.iter() {
                        // ★★ S7 promoted this from "replay the stored base vertex" to
                        // **"solve the definition"**: the corner's three planes share one
                        // motion, so solving their pre-motion names exactly (rational Cramer)
                        // and replaying that chain must reproduce the stored coordinate — the
                        // 8/8 measurement that let the base vertex die, now a permanent lock
                        // over every chain shape in this table.
                        let nacre_topo::VertexDef::ThreePlane(tri) = m.vertices.get(vh).def else {
                            panic!("a cuboid corner is a three-plane point");
                        };
                        let motion_of = |h| match m.surface_truth(h) {
                            nacre_topo::Surface::Plane { motion, .. } => *motion,
                            nacre_topo::Surface::Cylinder { motion, .. } => *motion,
                        };
                        // The solvable population mirrors `solid_points`' own criterion: the
                        // moved carriers share one leaf, and a world-stated carrier among
                        // them is provably fixed by that chain (the restatement licence).
                        // A corner whose carriers hold *different* leaves — a mixed-axis
                        // chain leaves a partially-fixed cap with a shorter history — is
                        // production's honest decline (Arrange), counted, not asserted.
                        let leaves: Vec<_> = tri.iter().filter_map(|&h| motion_of(h)).collect();
                        let leaf = *leaves.first().expect("a moved solid moves some carrier");
                        if !leaves.iter().all(|&l| l == leaf) {
                            declined += 1;
                            continue;
                        }
                        let mut coeffs = [[R::from_int(0); 4]; 3];
                        for (o, h) in coeffs.iter_mut().zip(tri) {
                            *o = *m.surface_name.get(&h).unwrap().narrow().unwrap();
                        }
                        for (c, h) in coeffs.iter().zip(tri) {
                            if motion_of(h).is_none() {
                                assert!(
                                    m.chain_fixes_plane(leaf, c),
                                    "a world-stated carrier must be provably fixed by the chain"
                                );
                            }
                        }
                        let base = nacre_scalar::three_planes_rat(coeffs)
                            .expect("three distinct planes of a cuboid corner meet");
                        let replayed = replay_chain_coord(
                            &m,
                            [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
                            leaf,
                        )
                        .unwrap();
                        assert_eq!(replayed, m.vertex_point(vh).as_array());
                        checked += 1;
                    }
                }
            }
        }
        // Both populations must be real, or half the lock is vacuous: the single-turn and
        // same-axis chains solve (fixed caps riding the licence), the mixed-axis chains
        // decline at their partially-fixed corners.
        assert!(checked > 80, "solved only {checked} vertices");
        assert!(
            declined > 0,
            "no mixed-leaf corner declined — the residual vanished?"
        );
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
        // The motion is the *face's* (S7: surfaces record it; vertices follow their planes).
        // Faces carry different depths since the restatement — the caps were fixed by the
        // Z turn (restated) and only joined at the X turn — so the order lock reads every
        // face: the walls' `[Z, X]` (never `[X, Z]` — the reversal this test exists to
        // forbid) and the caps' `[X]`.
        let sh = m.solids.get(r).outer;
        let mut counts: std::collections::HashMap<Vec<Axis>, usize> = Default::default();
        for &fh in &m.shells.get(sh).faces {
            let &nacre_topo::Surface::Plane {
                motion: Some(rotation),
                ..
            } = m.surface_truth(m.faces.get(fh).surface)
            else {
                panic!("every face of this solid records a motion");
            };
            let axes: Vec<Axis> = motion_chain(&m, rotation)
                .expect("an axis-aligned history holds no frame")
                .iter()
                .filter_map(|n| match n {
                    MoveNode::Rotate { axis, .. } => Some(*axis),
                    MoveNode::Translate { .. }
                    | MoveNode::Mirror { .. }
                    | MoveNode::Frame { .. }
                    | MoveNode::FrameWide(_)
                    | MoveNode::FrameThrough(_) => None,
                })
                .collect();
            *counts.entry(axes).or_insert(0) += 1;
        }
        assert_eq!(
            counts,
            [(vec![Axis::Z, Axis::X], 4), (vec![Axis::X], 2)]
                .into_iter()
                .collect(),
            "walls read root-to-leaf [Z, X]; the Z-fixed caps joined at X"
        );
    }

    /// Push a plane whose exact triple is `pts`, with an f64 `Plane` **consistent with it**
    /// (`Plane::through_points` of the realized corners) — the frame locks below ask about the
    /// realized basis's geometry, so unlike the interning lock the f64 form has to match.
    fn push_consistent(m: &mut Model, pts: [[R; 3]; 3]) -> Handle<Surface> {
        let f = |p: [R; 3]| Point3::from_array(p.map(|r| r.to_f64()));
        let pl = nacre_geom::Plane::through_points(f(pts[0]), f(pts[1]), f(pts[2]))
            .expect("a non-degenerate triple");
        let (h, _) = m.push_plane(pl, pts, None);
        h
    }

    /// The realized basis is a right-handed orthonormal frame whose origin sits on the plane
    /// and whose `ŵ` is parallel to the plane's normal — the sanity every frame lock needs.
    fn assert_frame_shape(
        m: &Model,
        h: Handle<Surface>,
        placement: &nacre_topo::FramePlacement,
        what: &str,
    ) {
        let (o, u, v, w) = frame_world_basis(m, h, placement, false)
            .unwrap_or_else(|| panic!("{what}: this frame must open"));
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        for (name, a) in [("u", u), ("v", v), ("w", w)] {
            assert!(
                (dot(a, a) - 1.0).abs() < 1e-12,
                "{what}: {name} is not unit ({:e})",
                dot(a, a) - 1.0
            );
        }
        assert!(dot(u, v).abs() < 1e-12, "{what}: u ⊥ v fails");
        assert!(dot(u, w).abs() < 1e-12, "{what}: u ⊥ w fails");
        assert!(dot(v, w).abs() < 1e-12, "{what}: v ⊥ w fails");
        let nacre_geom::Surface::Plane(pl) = m.surface(h) else {
            unreachable!()
        };
        assert!(
            pl.distance(Point3::from_array(o)).abs() < 1e-9,
            "{what}: the origin is off the plane by {:e}",
            pl.distance(Point3::from_array(o))
        );
        let n = pl.normal().as_array();
        let nn = dot(n, n).sqrt();
        let cross = [
            w[1] * n[2] - w[2] * n[1],
            w[2] * n[0] - w[0] * n[2],
            w[0] * n[1] - w[1] * n[0],
        ];
        assert!(
            dot(cross, cross).sqrt() / nn < 1e-9,
            "{what}: ŵ is not parallel to the plane normal"
        );
    }

    /// ★★★★★ S4: **a `Wide` name opens a frame.** The population `frame_chain` declined at
    /// `narrow()` — a plane whose canonical answer exceeds `i128` — now realizes its canonical
    /// placement through the arbitrary-precision road, and the basis is a real frame on the
    /// real plane. (What stays closed for `Wide` is the *narrow shortcuts* — `base_rat`,
    /// Shewchuk, `Isometry` transport — locked on the topo side.)
    #[test]
    fn a_wide_plane_hosts_a_canonical_frame() {
        let q = |n: i128, d: i128| R::new(n, d).unwrap();
        let big1 = (1i128 << 90) + 1;
        let big2 = (1i128 << 90) + 3;
        let pts = [
            [q(big1, 3), q(big2, 7), q(0, 1)],
            [q(-big2, 5), q(big1, 11), q(0, 1)],
            [q(1, 13), q(1, 17), q(1, 19)],
        ];
        // Fixture qualification (the S2 census lesson): genuinely wide.
        let name = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
        assert!(name.narrow().is_none(), "the fixture must be wide");
        let mut m = Model::new();
        let h = push_consistent(&mut m, pts);
        assert_frame_shape(&m, h, &nacre_topo::FramePlacement::Canonical, "wide plane");
    }

    /// ★★★★ S4: **a narrow name whose squared lengths overflow `i128` opens too** — the
    /// measured 1.6% population (`n·n` is a square, so it overflows long before the name).
    /// Before S4 this was `plane_frame_named`'s hard `None`; the v-fallback never applied.
    #[test]
    fn a_narrow_name_with_wide_squares_hosts_a_frame() {
        let q = |n: i128, d: i128| R::new(n, d).unwrap();
        // Intercept form: the plane through (1/p, 0, 0), (0, 1/q, 0), (0, 0, 1/r) has the
        // canonical name [p, q, r, −1] — narrow when p, q, r fit i128, while n·n = p²+q²+r²
        // does not (~2^140).
        let (p, q2, r) = ((1i128 << 70) + 1, (1i128 << 70) + 3, (1i128 << 70) + 7);
        let pts = [
            [q(1, p), q(0, 1), q(0, 1)],
            [q(0, 1), q(1, q2), q(0, 1)],
            [q(0, 1), q(0, 1), q(1, r)],
        ];
        // Fixture qualification: the name is narrow AND the narrow frame derivation dies on it.
        let name = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
        let c = *name.narrow().expect("the name itself fits i128");
        assert!(
            nacre_scalar::plane_frame_default(c).is_none(),
            "the fixture must be in the nn-overflow population"
        );
        let mut m = Model::new();
        let h = push_consistent(&mut m, pts);
        assert_frame_shape(
            &m,
            h,
            &nacre_topo::FramePlacement::Canonical,
            "nn-overflow plane",
        );
    }

    /// ★★★ S6a: **a `Named` placement opens on a `Wide` name too.** S4's locks covered the
    /// canonical road; the axes-only population (`from_axes`, a tilted full-width frame) is
    /// exactly the other pairing — a caller-stated origin/`ref_dir` on a plane whose canonical
    /// name exceeds `i128` — and it goes through `WideFrame::named_of`.
    ///
    /// ★ The fixture is that population, not the S2 `2^90` triple: full-width *decimal* axes at
    /// CAD scale. Their cross runs the denominators to `10^48`, so the canonical name is
    /// genuinely `Wide` (asserted), while the geometry stays near `1` — which matters, because
    /// a frame's unit axes are invisible in f64 next to a `2^90` origin (an ulp there is
    /// `~6e10`), and the first version of this test proved it by accident.
    #[test]
    fn a_wide_plane_hosts_a_named_frame() {
        let d = |x: f64| R::from_decimal(x).unwrap();
        let o = [
            d(0.2547863291057384),
            d(-0.5123456789012345),
            d(1.5432109876543211),
        ];
        let x = [
            d(0.7123456789012345),
            d(0.5876543210987654),
            d(0.4098765432101234),
        ];
        let y = [
            d(-0.5876543210987654),
            d(0.7123456789012345),
            d(0.1234567890123456),
        ];
        let add = |a: [R; 3], b: [R; 3]| -> [R; 3] {
            core::array::from_fn(|i| a[i].checked_add(b[i]).expect("decimal widths"))
        };
        let pts = [o, add(o, x), add(o, y)];
        // Fixture qualification: genuinely wide, at unit scale.
        let name = nacre_scalar::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
        assert!(
            name.narrow().is_none(),
            "the fixture must be wide — full-width crosses were expected to exceed i128"
        );
        let mut m = Model::new();
        let h = push_consistent(&mut m, pts);
        // The caller's statement, `PlaneDef`-style: origin = first point, +u toward the second.
        let placement = nacre_topo::FramePlacement::Named {
            origin: o,
            ref_dir: x,
        };
        assert_frame_shape(&m, h, &placement, "wide plane, named");
        // And `û` runs along the stated `ref_dir`, not some canonical direction.
        let (_, u, _, _) = frame_world_basis(&m, h, &placement, false).unwrap();
        let rd = x.map(|r| r.to_f64());
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let cos = dot(u, rd) / dot(rd, rd).sqrt();
        assert!(
            (cos - 1.0).abs() < 1e-9,
            "û does not follow the caller's ref_dir (cos = {cos})"
        );
    }
}
