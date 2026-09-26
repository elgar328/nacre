use super::*;
use crate::BoolKind;
/// A face-local feature built as **tool body + boolean**: the profile
/// extrudes off `face` into a top-flush prism, then `kind` fuses/cuts it against the face's solid.
/// `Fuse` sweeps **outward** (a boss); `Cut` sweeps **inward** (a blind pocket). Returns the result
/// solid and the feature's
/// exposed cap — the boss top or the pocket floor, the outer-shell face on the prism's far-cap plane
/// with outward normal `+n`. A cap that did not survive (a through-cut has no floor) is `no_cap`,
/// **the caller's error raised here** rather than an `Option` the caller turns into one: only this
/// scope holds the boolean's result solids, and restoring `live_solids` from them is what keeps a
/// reject from committing. `NonPositiveDistance`/`DegenerateProfile` propagate from the frame;
/// overhang configurations the boolean does not cover surface as `Boolean(_)`.
fn extrude_and_boolean(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
    kind: BoolKind,
    no_cap: OpError,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    if nacre_exact::Rat::from_decimal(dist).is_none() {
        return Err(OpError::DistOutsideDecimalWindow);
    }
    profile.check()?;
    refuse_non_quarter_arcs(profile)?;
    let frame = face_frame(model, face)?;
    // No containment check — an overhanging footprint routes to the overhang boolean sidecars.
    // The face's outward, realized — the direction the prism's caches are built along, a cache's
    // use; the prism itself is built from the frame.
    let n = frame_basis(model, &frame.frame)
        .map(|b| Vector3::from_array(b.3))
        .ok_or(OpError::PlaneWithoutExactForm)?;
    // Cut carves inward, Fuse raises outward; either way the prism's near cap is flush on the face.
    let signed = if matches!(kind, BoolKind::Cut) {
        -dist
    } else {
        dist
    };
    // ★★★ **The frame is not decided here — `face_frame` already decided it**, and that is the
    // point: `face_plane` promises a caller the frame this operation will use, so there must be
    // exactly one place that picks it. The rings take the one road from a frame to a prism the
    // extrude takes too.
    let (outer, holes) = frame_rings(model, &frame.frame, profile, signed)?;
    let (prism, prism_faces) = build_prism(
        model,
        outer,
        holes,
        n * signed.signum(),
        // The base cap is flush on the face and faces against the sweep: back into the solid on a
        // pad (the sweep runs along the face's outward), out of it on a pocket — so on the face's
        // own surface it is the face's orientation, reversed for a pad.
        Some((
            frame.surface_h,
            match kind {
                BoolKind::Cut => model.face(face).orientation,
                _ => model.face(face).orientation.flipped(),
            },
        )),
        // The pad reuses the face's own surface, so it pushes no base cap and states nothing.
        None,
    )?;
    let (solids, class_of) =
        crate::boolean::boolean_with_classes(model, kind, frame.solid_h, prism).map_err(|e| {
            model.supersede_live(&[prism]); // drop the transient prism (atomic on failure)
            OpError::Boolean(e)
        })?;
    // A pad's prism must actually meet the face. Two live solids are each one outer shell, so a
    // `Fuse` of them can only come back severed if they never touched — the footprint missed the
    // face entirely. Returning the piece that carries the cap would hand back a floating boss and
    // silently drop the base, so this is `PadMissesFace`: the boolean succeeded and answered
    // correctly (two solids); it is the *pad's* premise that broke. A footprint that merely
    // overhangs still touches, fuses into one solid, and takes the normal path.
    //
    // `assemble_fuse_cut` already retired the inputs, so restoring `live_solids` is what keeps the
    // "no reject-after-commit" contract true from the outside: the model the caller sees is the one
    // it had before. Only `live_solids` is touched — the store stays append-only.
    if matches!(kind, BoolKind::Fuse) && solids.len() > 1 {
        model.supersede_live(&solids);
        model.make_live(frame.solid_h);
        return Err(OpError::PadMissesFace);
    }
    // Exposed cap = the result face on the prism's far-cap plane (face plane offset by n·signed),
    // its outward normal +n (the opening side for a pocket, the boss top for a boss). A pocket (Cut)
    // that severs leaves several solids — scan them all for the cap and return the piece that
    // carries it, leaving the others live; that is a valid multi-solid model, not a failure.
    // `build_prism` returns the far cap as `faces[1]`; the store is append-only, so it is still
    // readable after the boolean retired the prism, and it names the cap plane exactly.
    let far_cap = prism_faces[1];
    // ★★★ **Which surface the cap's plane became** — the boolean's own answer, not a guess made
    // afterwards. Its plane classes are decided with evidence (an exact `orient3d`, a composed
    // rotation proof, or a coincidence within the limit), and a later comparison of handles or
    // coordinates can see none of that: on a tilted face the cap merges with a face of the other
    // operand and the survivor carries *that* surface.
    //
    // ★ **And which way it faces there, from the truth.** The cap the caller gets back faces `+n`:
    // the far cap itself on a pad (it was swept along `+n`), and its reverse on a pocket (the tool
    // cap faced along the sweep, `−n`; the floor it leaves faces back into the void). So the
    // expected orientation on the class's surface is the far cap's, flipped for a `Cut`.
    let (want_surf, far_orientation) = class_of.get(&far_cap).copied().unwrap_or_else(|| {
        let f = model.face(far_cap);
        (f.surface, Ok(f.orientation))
    });
    let want_orientation = match far_orientation {
        Ok(o) => match kind {
            BoolKind::Cut => o.flipped(),
            _ => o,
        },
        // ★ **The judge could not say which way the cap runs on its class's surface** — the cap is
        // there, its facing is not decided, so this is neither the pad's `no_cap` nor a found
        // cap: it is refused by the judgement's own name, the rule every undecided judgement the
        // result rests on follows (`arrangement`'s `undecided_reject`). The restore is the no-cap
        // arm's below, for the same reason.
        Err(outcome) => {
            model.supersede_live(&solids);
            model.make_live(frame.solid_h);
            return Err(OpError::Boolean(crate::reject(match outcome {
                nacre_judge::Decision::Exhausted { .. } => crate::RejectReason::JudgeExhausted,
                _ => crate::RejectReason::DegenerateWitness,
            })));
        }
    };
    match solids.iter().find_map(|&s| {
        find_face_coplanar_with(model, s, (want_surf, want_orientation)).map(|c| (s, c))
    }) {
        Some((solid, cap)) => Ok((solid, cap)),
        // Nothing carries the cap — either the prism reached through (a pocket with no floor) or
        // it removed the solid outright, which is the same verdict taken to its limit. Only `Cut`
        // can empty a result: a `Fuse` of two non-empty solids is never empty.
        //
        // ★ **The restore is the whole reason this arm lives here.** `assemble_fuse_cut` already
        // retired the operand and installed its own results, so returning an error now would hand
        // the caller a failure *and* a model it never asked for: its solid gone, a through-cut in
        // its place. Putting `live_solids` back is what makes the reject true from the outside —
        // the same move `PadMissesFace` makes above, for the same reason. (The arena keeps the
        // prism's cells; the store is append-only. That residue is why a session must rebuild
        // from its log before recording again — see `tests/invariants/replay.rs`.)
        None => {
            debug_assert!(
                !solids.is_empty() || matches!(kind, BoolKind::Cut),
                "a Fuse cannot produce an empty result"
            );
            model.supersede_live(&solids);
            model.make_live(frame.solid_h);
            Err(no_cap)
        }
    }
}

/// Pad a boss on a planar `face`: extrude the profile **outward** by `dist` and `Fuse` it onto the
/// solid, adding `profile_area · dist` of material. Returns `(new solid, top cap face)`. A boss
/// always yields its top cap, so the `None` guard is an unreachable internal-invariant defense.
pub(crate) fn pad(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    extrude_and_boolean(
        model,
        face,
        profile,
        dist,
        BoolKind::Fuse,
        OpError::DegenerateGeometry,
    )
}

/// Carve a blind pocket on a planar `face`: extrude the profile **inward** by `dist` and `Cut` it
/// from the solid, removing `profile_area · dist` of material. Returns `(new solid, floor face)`.
/// `PocketNotBlind` if `dist` reaches through the solid (the far cap is not blind → a through-cut
/// with no floor face).
pub(crate) fn pocket(
    model: &mut Model,
    face: Handle<Face>,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Handle<Face>), OpError> {
    extrude_and_boolean(
        model,
        face,
        profile,
        dist,
        BoolKind::Cut,
        OpError::PocketNotBlind,
    )
}

/// The outer-shell face of `solid` that lies on the cap's plane class and runs the cap's way — how
/// a pad/pocket recovers its own exposed cap (the boss top, the pocket floor) from the boolean
/// result. `None` if there is none (a through-pocket has no floor).
///
/// **Asked of the boolean's own answer, in integers.** `assemble_fuse_cut` gives a result face the
/// surface of its plane class's representative, which `class` names, and its direction is decided
/// exactly too: `class` carries the `Orientation` the cap must have on that surface, read off the
/// truth (`planes::face_facing`) — so this compares a handle and a flag and reads no coordinate
/// and no normal.
///
/// `class_of` carries every planar input face the arrangement touches, so no comparison of
/// coordinates is asked.
///
/// If the cap survives as several faces they all satisfy this, and the first is returned.
pub(crate) fn find_face_coplanar_with(
    model: &Model,
    solid: Handle<Solid>,
    class: (Handle<Surface>, Orientation),
) -> Option<Handle<Face>> {
    let shell = model.solid(solid).outer;
    model.shell(shell).faces.iter().copied().find(|&fh| {
        let face = model.face(fh);
        face.surface == class.0 && face.orientation == class.1
    })
}
