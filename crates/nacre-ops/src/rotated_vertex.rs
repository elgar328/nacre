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

use nacre_judge::{MoveNode, WitnessPoint};
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
        .coord())
}

/// `p` carried through `chain`, in the producer's own order and float operations — the one
/// definition of "replay" this crate has, so a coordinate and its definition cannot drift apart.
///
/// `None` only for a [`MoveNode::Frame`] whose squared lengths do not fit `i128`, which
/// [`motion_chain`] already refuses to emit — so in practice this is infallible, and the `Option`
/// is here so that "in practice" does not have to be an invariant spanning two crates.
pub(crate) fn replay(p: WitnessPoint, chain: &[MoveNode]) -> Option<WitnessPoint> {
    p.apply_chain(chain)
}

/// **A surface's plane as its exact witness triangle** — the per-carrier building block of an
/// implicit point. The same computation `collect_planes`' arms perform per face, spelled
/// once for the per-carrier consumer: `Known` points replay through the surface's own motion; a
/// **named** `Through` plane (rational closure) solves through `through_points_rat` and replays
/// the same way. ★ A named `Through` plane whose meets are wider than any witness base
/// gets its **own frame's canonical probes** instead — on-plane exact definitions, the
/// same witness `collect_planes`' probe branch builds; `frame_chain` already carries the
/// plane's later motion, so those return directly. `None` for a nameless `Through` carrier
/// (a datum on a nameless datum — depth) and for a cylinder.
pub(crate) fn surface_witness_triangle(
    model: &Model,
    h: Handle<Surface>,
) -> Option<[WitnessPoint; 3]> {
    let (base, motion) = match model.surface(h) {
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
/// [`nacre_judge::JudgedPoint::Pure`]; anything else that is still a
/// well-defined
/// three-plane meet becomes [`nacre_judge::JudgedPoint::Meet`] of its carriers' witness triangles — the
/// straddling population, **and** the pure-but-too-wide one, which a meet represents
/// without ever asking the coordinate to fit anything.
///
/// `None` for a seam vertex, or a carrier that has no witness triangle of its own (a nameless
/// `Through` carrier — depth — or a chain outside the decimal window). The producer tells those
/// apart by cause before any plane is pushed; here one answer suffices, because reaching this
/// through an *accepted* statement re-derives the same deterministic result.
pub(crate) fn through_judged_points(
    model: &Model,
    vs: [Handle<nacre_topo::Vertex>; 3],
) -> Option<[nacre_judge::JudgedPoint; 3]> {
    use nacre_judge::JudgedPoint;
    let mut out: [Option<JudgedPoint>; 3] = [None, None, None];
    for (o, vh) in out.iter_mut().zip(vs) {
        let tri = match *model.vertex(vh) {
            nacre_topo::Vertex::ThreePlane(tri) => tri,
            // OnSeam pins a curve, not a point; a Pierce point's coordinates are
            // quadratic-irrational, and this table's witnesses are rational by type —
            // both decline, per variant (the judging of pierce points is its own machinery,
            // not this road).
            nacre_topo::Vertex::OnSeam(_) | nacre_topo::Vertex::Pierce { .. } => {
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
/// ★ Width is **not** on that list: a `Wide` name or overflowing squared lengths take the
/// arbitrary-precision road (`MoveNode::FrameWide`) instead of declining.
///
/// ★★ **No caller falls back to an f64 path on this `None`.**
/// `planes` turns this `None` into `RejectReason::FrameOutOfRange`, `exact.rs` into
/// `OpError::PlaneWithoutExactForm` (there is no f64 prism road), `reuse.rs` declines to the
/// arrangement — an equally correct road, not a degraded one — and `replay_chain_coord` names
/// `WitnessPointError::Downgrade`. Every consumer is honest.
pub(crate) fn motion_chain(model: &Model, leaf: Handle<MotionNode>) -> Option<Vec<MoveNode>> {
    let mut chain = Vec::new();
    let mut cur = Some(leaf);
    while let Some(h) = cur {
        let n: &MotionNode = model.motion(h);
        match n.motion {
            Motion::Rotate { axis, pivot, angle } => {
                chain.push(MoveNode::Rotate { axis, angle, pivot })
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
        // ★★★ **The judged road**: a plane with no name at all — a
        // mixed-frame `Through` statement, whose exact world coefficients are irrational. Its
        // frame is derived from the defining points as intervals ([`nacre_judge::FrameThrough`]),
        // with the branch decided once at the fixed rung, so the same statement always frames
        // the same way. `Named` placement is refused here defensively (`SketchFrame::named`
        // already rejects it by name at construction — an on-plane claim needs a name to verify
        // against), and the plane's own later motion is appended by the shared tail below,
        // exactly as for the named roads.
        let nacre_topo::Surface::Plane {
            points: nacre_topo::PlanePoints::Through(vs),
            motion,
        } = model.surface(plane)
        else {
            return None;
        };
        if !matches!(placement, FramePlacement::Canonical) {
            return None;
        }
        let pts = through_judged_points(model, *vs)?;
        let mut chain = vec![MoveNode::FrameThrough(Box::new(
            nacre_judge::FrameThrough::of(pts, flip)?,
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
    // **The narrow road.** For `Canonical` the placement pair is
    // derived from the canonical (unflipped) coefficients first and the sign is applied after —
    // a correctness condition: `ref_dir = ẑ × n` is sign-sensitive.
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
    // **The wide road** — the same convention through arbitrary precision, where nothing
    // can overflow. Taken when the name is `Wide` or when any narrow derivation step hits
    // `i128` (the measured 1.6% `n·n` population). So
    // a plane with a name can always host a sketch.
    let wide_road = || -> Option<MoveNode> {
        let wf = match placement {
            FramePlacement::Canonical => nacre_judge::WideFrame::canonical_of(name, flip)?,
            FramePlacement::Named { origin, ref_dir } => {
                nacre_judge::WideFrame::named_of(name, origin, ref_dir, flip)?
            }
        };
        Some(MoveNode::FrameWide(wf))
    };
    let mut chain = vec![narrow_road().or_else(wide_road)?];
    match model.surface(plane) {
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
        Some(replay(WitnessPoint::at(p.map(Rat::from_int)), &chain)?.coord())
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
#[path = "tests/rotated_vertex.rs"]
mod tests;
