use super::*;
/// A [`SketchFrame`] with its `flip` measured — the placement's realized `ŵ` dotted against the
/// direction the sketch must face (a sweep's sense, a face's outward normal). This is the one
/// place `flip` is decided: every road calls it, so no two can measure differently. `None`
/// when the chain cannot realize a basis (a plane with no name).
///
/// ★★★ **Only the frame comes back — deliberately.** The basis realized here to measure `flip`
/// is the *unflipped* one, and it used to ride along "so deciding and looking cost one
/// realization, not two". That saving is what broke `face_plane`'s contract: the one caller who
/// wanted the axes combined the measured `flip` with the unflipped basis **by hand**, as a
/// half-turn about `v̂` — while the realization (`frame_chain`) half-turns about `û` — and every
/// flip=true face was reported a frame point-symmetric to the one the pad actually built in
/// (measured: a footprint centred on the face through `face_plane`'s own coordinates landed
/// outside it, `PadMissesFace` on 2 of 6 faces of a turned block). A caller that needs the axes
/// asks [`crate::rotated_vertex::frame_world_basis`] *with the measured flip*, so the geometry of
/// `flip` is written in exactly one place; the second 4-point replay is one plain f64 chain per
/// user operation, which is what the hand-combination was saving.
pub(super) fn measured_frame(
    model: &Model,
    plane: Handle<Surface>,
    placement: nacre_topo::FramePlacement,
    toward: Vector3,
) -> Option<SketchFrame> {
    let basis = crate::rotated_vertex::frame_world_basis(model, plane, &placement, false)?;
    let n = toward.as_array();
    let flip = (0..3).map(|k| basis.3[k] * n[k]).sum::<f64>() < 0.0;
    Some(SketchFrame {
        plane,
        placement,
        flip,
    })
}

/// Name a [`SketchFrame`] as the [`nacre_topo::Motion::Frame`] node the sweep writes coordinates
/// against — the one road from the frame value to a node, shared by the extrude and face paths.
///
/// ★★ **`push_motion` interns**, so two sketches in one frame name the *same* node — which is
/// what makes their surfaces intern too (`SurfaceKey` is `(name, motion)`): two routes to one
/// height become one `Handle<Surface>` at construction, with no f64 comparison anywhere. With
/// `Canonical` placement the node is `(plane, Canonical, flip)` — nothing per-sketch in the key.
pub(super) fn push_frame_node(
    model: &mut Model,
    frame: SketchFrame,
) -> Handle<nacre_topo::MotionNode> {
    model.push_motion(
        nacre_topo::Motion::Frame {
            plane: frame.plane,
            placement: frame.placement,
            flip: frame.flip,
        },
        None,
    )
}

/// **Which way is "right" and "up" on a face pointing `n`** — the `(u, v)` a sketch frame takes.
///
/// This is the **arbitrary-axis convention** (DXF/AutoCAD, and what most CAD puts on a face):
/// cross the world `ẑ` into the normal, unless the normal *is* vertical, in which case cross `ŷ`.
///
/// ★ **Not [`nacre_math::Vector::any_perpendicular`], and the difference is the point.** That one
/// answers a question of fact — *"give me a unit vector perpendicular to this"* — by crossing in
/// whichever world axis the normal is **least** aligned with, which keeps the cross far from zero
/// and is exactly right for a cylinder's seam or a STEP `ref_dir`. It is the wrong answer for a
/// frame a person draws in, because the axis it picks changes with the normal: a box lid comes out
/// `u = −ŷ, v = +x̂`, ninety degrees from world XY and from [`SketchPlane::world_xy`] itself.
///
/// What this convention buys, on top of the lid agreeing with `world_xy`:
///
/// * **On every face that is not horizontal, `v` points up.** `u = ẑ × n` is horizontal, so
///   `v·ẑ = (n × (ẑ × n))·ẑ = 1 − n_z²`, which is positive unless `n` is vertical. Sketching on a
///   wall, "up" is up.
/// * **Axis-aligned faces keep axes in `{0, ±1}`**, so `SketchPlane::exact` still fires and the
///   rational construction path is not lost.
///
/// ★★ **The branch is exact, not toleranced.** DXF switches on `|n_x| < 1/64`, a threshold only a
/// float-only kernel needs; `n_x == 0 && n_y == 0` is the real question and this kernel can ask it.
/// Nor is the non-vertical branch fragile near vertical: `ẑ × n = (−n_y, n_x, 0)` is a plain
/// rotation of `(n_x, n_y)` into the plane, with no cancellation to lose digits to.
///
/// ★ **A discontinuity is unavoidable and this convention chooses where to put it** — no continuous
/// tangent frame exists on the sphere. `any_perpendicular` breaks along whole arcs (wherever two
/// components tie for smallest, e.g. `(0.5, 0.5, 0.707)`, far from any pole); this breaks at the
/// two poles only. It is not smooth *at* the poles either — approaching `+ẑ` from different sides
/// gives different limits — but the set where that happens is two points instead of three arcs.
///
/// Returns both axes rather than just `u`, so the two call sites cannot disagree about which way
/// `v` runs. `None` only for the zero vector.
pub(crate) fn frame_axes(n: Vector3) -> Option<(Vector3, Vector3)> {
    let [nx, ny, nz] = n.as_array();
    let raw = if nx == 0.0 && ny == 0.0 {
        // ŷ × n for a vertical normal: (n_z, 0, 0).
        Vector3::from_array([nz, 0.0, 0.0])
    } else {
        // ẑ × n.
        Vector3::from_array([-ny, nx, 0.0])
    };
    // `-0.0` is worth nothing to keep: the negations above produce it whenever a component is
    // zero, and a frame reported to a caller as `[-0.0, 1.0, 0.0]` invites a double-take for no
    // reason. Adding zero is the identity on every other value.
    let tidy = |v: Vector3| Vector3::from_array(v.as_array().map(|c| c + 0.0));
    let u = tidy(raw.normalize()?);
    Some((u, tidy(n.cross(u))))
}

/// The sketch plane of a planar face — **the very frame [`Operation::PadOnFace`] and
/// [`Operation::PocketOnFace`] place their profile in**, so a caller can work out where its
/// `(0, 0)` will land before it builds anything.
///
/// That equality is the contract, not a coincidence: this is a projection of the frame those
/// operations use, never a second derivation. A test pins a hand-placed profile against a pad to
/// keep it that way.
///
/// `NonPlanarFace` for a curved surface (only a plane carries a frame); `FaceNotInLiveSolid` if no
/// live solid's outer shell holds the face.
pub fn face_plane(model: &Model, face: Handle<Face>) -> Result<SketchPlane, OpError> {
    let f = face_frame(model, face)?;
    Ok(realized_plane(f.origin, f.x, f.y))
}

/// **Where a frame is, in space** — its origin and its two axes, realized as f64.
///
/// [`face_plane`] answers this for a face; this answers it for any frame a caller holds, which
/// is what a viewer needs to draw a sketch where it was drawn: the sketch's coordinates are
/// `(u, v)` in this frame, and `origin + u·x + v·y` is the point.
///
/// ★ It is a **report**, not a truth. The frame's statement is the truth — this is that
/// statement realized, with all the rounding a realization carries, and nothing exact should be
/// decided from it.
///
/// `None` when the frame's chain cannot be realized at all (a plane with no name). A caller
/// that cannot place a thing should decline to draw it rather than draw it somewhere wrong.
pub fn frame_plane(model: &Model, frame: &SketchFrame) -> Option<SketchPlane> {
    let basis = crate::rotated_vertex::frame_world_basis(
        model,
        frame.plane(),
        frame.placement(),
        frame.flip(),
    )?;
    Some(realized_plane(
        Point3::from_array(basis.0),
        Vector3::from_array(basis.1),
        Vector3::from_array(basis.2),
    ))
}

/// A planar face's sketch frame **as a [`SketchFrame`]** — the plane handle, placement, and
/// measured flip that [`Operation::PadOnFace`] / [`Operation::PocketOnFace`] sketch in. Where
/// [`face_plane`] projects that frame to realized f64 axes for a caller to *look at*, this is
/// the exact vocabulary itself: the same value `face_frame` builds internally, not
/// thrown away at the boundary.
///
/// ★★★ **What comes back is verified against the pad's frame, by realization, to the bit.** On a
/// world-branch face (axes lifting exactly) the operation elides the frame node and sketches in
/// `face_frame`'s world axes — so this transcribes *those* axes into the vocabulary and returns a
/// candidate only if realizing it lands bit-identically on them. The old fallback assumed the
/// canonical frame realizes to the same axes ("the node-omission normalization"); that holds only
/// for flip=false faces of motion-free planes, and everywhere else the returned frame put a
/// sketch somewhere the pad does not (measured: point-symmetric on every flip=true axis-aligned
/// face). Verification is the contract now — no candidate can be returned wrong, whatever
/// population shows up next.
///
/// A face has no caller to name a placement, so the canonical frame is tried first (the stronger
/// normal form); where it realizes elsewhere, the pad's axes are transcribed as a `Named`
/// placement (origin + `ref_dir`) with `flip` measured as everywhere else.
///
/// Errors as [`face_plane`]: `NonPlanarFace`, `FaceNotInLiveSolid`; `PlaneWithoutExactForm` when
/// the plane carries no name to derive a frame from (a test-only unregistered surface); and
/// [`OpError::FrameNotRepresentable`] when the frame exists but no spelling realizes to it —
/// world-branch faces whose surface carries a motion; since the invariant-plane restatement
/// a plane its motion fixes carries none, so the residual population is the recorded one
/// (see the variant's doc).
pub fn face_sketch_frame(model: &Model, face: Handle<Face>) -> Result<SketchFrame, OpError> {
    let f = face_frame(model, face)?;
    if let Some(sf) = f.sketch_frame {
        return Ok(sf);
    }
    // The world-branch population: the operation will build no node and sketch in `f`'s axes.
    // Every candidate below must prove itself by realizing to exactly those axes — bits, not a
    // tolerance: both sides come from exact roads ({0,±1} axes, rational projections), so
    // agreement is exact when it holds and a threshold would only paper over a third derivation.
    let pad_frame =
        [f.origin.as_array(), f.x.as_array(), f.y.as_array()].map(|c| c.map(f64::to_bits));
    let verified = |sf: SketchFrame| -> Option<SketchFrame> {
        let (o, u, v, _) =
            crate::rotated_vertex::frame_world_basis(model, f.surface_h, sf.placement(), sf.flip)?;
        ([o, u, v].map(|c| c.map(f64::to_bits)) == pad_frame).then_some(sf)
    };
    let canonical = || {
        measured_frame(
            model,
            f.surface_h,
            nacre_topo::FramePlacement::Canonical,
            f.n,
        )
    };
    // The transcription: the pad's own origin and +u, said as a `Named` placement. `named`'s
    // exact checks (on-plane origin, non-degenerate ref_dir) ride along; any failure just drops
    // the candidate — the refusal below is the answer, never a silent wrong frame.
    let transcribed = || {
        let sf = SketchFrame::named(model, f.surface_h, f.origin, f.x).ok()?;
        measured_frame(model, f.surface_h, sf.placement, f.n)
    };
    if !model.surface_name.contains_key(&f.surface_h) {
        // No name at all: nothing can realize. The distinct, older proposition.
        return Err(OpError::PlaneWithoutExactForm);
    }
    canonical()
        .and_then(&verified)
        .or_else(|| transcribed().and_then(&verified))
        .ok_or(OpError::FrameNotRepresentable)
}

/// Locate `face`'s live solid and build its planar frame. `NonPlanarFace` for a curved surface,
/// `FaceNotInLiveSolid` if no live outer shell holds it.
pub(super) fn face_frame(model: &Model, face: Handle<Face>) -> Result<FaceFrame, OpError> {
    let (solid_h, _) = model
        .live_solids()
        .iter()
        .map(|&s| (s, model.solid(s).outer))
        // Index-only equality again: a face handle from another model can match here. The
        // shell lookup that follows is where the cross-store guard fires.
        .find(|&(_, sh)| model.shell(sh).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.face(face);
    let surface_h = f.surface;
    let orientation = f.orientation;
    let plane = match model.surface_cache(surface_h) {
        nacre_geom::Surface::Plane(p) => *p,
        nacre_geom::Surface::Cylinder(_) => return Err(OpError::NonPlanarFace),
    };
    let sign = f64::from(orientation.sign());
    let n = plane.normal() * sign;
    let (x, y) = frame_axes(n).ok_or(OpError::DegenerateGeometry)?;
    // ★ The origin is the **world origin projected onto the face's plane** — a property of the
    // plane, not of the face.
    //
    // It used to be the face region's area centroid, chosen over the mean of the outer loop's
    // corners because a centroid does not move when a vertex is added along a straight edge. Both
    // are computed in `f64` from the face's own vertices, though, and `construct.rs` lifts the frame
    // origin with `Rat::from_decimal` — so a rounded cache became the truth. Padding one footprint
    // twice then placed the second profile an ulp from the first and left faces of area `2.2e-16`
    // that `validate` did not report.
    //
    // The projection is `(−d / n·n)·n` from the plane's rational coefficients: one division, no
    // f64 in the derivation, and invariant under negating or scaling those coefficients — so two
    // faces of one plane cannot disagree about where `(0, 0)` is. Measured over the suite: for
    // every `Constructed` surface the realized point lies on the f64 plane at distance exactly `0`
    // (1121/1121), which the centroid did not always manage.
    //
    // Only a world-stated plane's coefficients are world truth (`motion: None`). A moved
    // surface records its **pre-motion** frame, so projecting those gives a pre-motion point —
    // measured `0.29` away from the world plane, not a rounding but a different place. Those
    // keep the f64 projection, which is the same rule computed from the description available.
    // `narrow()` gates the wide vessel out: a `Wide` name carries identity only, so it
    // keeps the f64 projection exactly as a missing name did.
    let world_stated = matches!(
        model.surface(surface_h),
        nacre_topo::Surface::Plane { motion: None, .. }
    );
    let origin = match (
        world_stated,
        model.surface_name.get(&surface_h).and_then(|n| n.narrow()),
    ) {
        (true, Some(&c)) => nacre_exact::plane_origin_projection(c)
            .map(|p| Point3::from_array([p[0].to_f64(), p[1].to_f64(), p[2].to_f64()]))
            .unwrap_or_else(|| plane.project(Point3::origin())),
        _ => plane.project(Point3::origin()),
    };
    // ★★★★★ **When the operation will sketch in the plane's own frame, report *that* frame.**
    //
    // `face_plane`'s contract is that it names the frame `PadOnFace` places a profile in, and a
    // tilted face is about to be sketched in its plane's frame rather than in world coordinates.
    // Reporting the world axes here and using the frame's there would put a caller's profile a
    // quarter turn from where it asked for it — measured, as a boss that missed its own face.
    //
    // ★★ **The gate is the same one the operation uses**, and it has to be the same expression,
    // not the same intent: take the frame only when the world axes do not lift to exact
    // orthonormal rationals. Axis-aligned faces therefore never go near it and are untouched.
    //
    // ★ **`flip` is measured, not derived** — by `measured_frame`, the one measuring place.
    // ★★★ **The axes are then realized *with* that flip, by the same function the operation's
    // prism replays through.** They used to be read off the unflipped basis with the sign applied
    // by hand here, as a half-turn about `v̂` — but the realization (`frame_chain`) half-turns
    // about `û` (its `ref_dir` is derived from the unflipped coefficients and survives the sign),
    // so every flip=true face was reported a frame point-symmetric to the one the pad built in.
    // Asking `frame_world_basis` with the measured flip leaves the geometry of `flip` written in
    // exactly one place; for flip=false the call is bit-identical to the measuring one.
    // ★★ A face has no caller to name a frame, so its placement is `Canonical` — derived
    // when the chain is flattened, stored nowhere. That is also what opens this branch for a
    // plane whose name is `Wide` or whose canonical values overflow `i128`: `frame_world_basis`
    // succeeds through the arbitrary-precision road where a narrow derivation would decline.
    let world = realized_plane(origin, x, y);
    let sketch = (world.exact().is_none())
        .then(|| measured_frame(model, surface_h, nacre_topo::FramePlacement::Canonical, n))
        .flatten()
        .and_then(|sf| {
            let (o, u, v, _) = crate::rotated_vertex::frame_world_basis(
                model,
                surface_h,
                &nacre_topo::FramePlacement::Canonical,
                sf.flip,
            )?;
            Some((
                sf,
                Point3::from_array(o),
                Vector3::from_array(u),
                // ★★ **`v̂` as realized, not as `ŵ × û` recomputed here.** It has its own exact
                // rational form (`plane_frame`), so realizing it costs one rounding where a cross
                // product costs two that do not cancel — measured, a wall whose `v` is exactly
                // `ẑ` came back three ulps short of `1.0` through the cross product.
                Vector3::from_array(v),
            ))
        });
    let (x, y, origin, sketch_frame) = match sketch {
        Some((sf, o, u, v)) => (u, v, o, Some(sf)),
        None => (x, y, origin, None),
    };
    Ok(FaceFrame {
        solid_h,
        surface_h,
        n,
        x,
        y,
        origin,
        sketch_frame,
    })
}

/// A [`SketchPlane`] from a **realized** (rounded) frame basis — kernel-internal, and
/// deliberately without a definition.
///
/// ★★ **Do not lift these axes.** The public [`SketchPlane::from_axes`] lifts what a *caller*
/// wrote, because written decimals are a statement. A realized basis is a cache of a plane that
/// already has an exact definition (points + motion); lifting it would mint a second,
/// ulp-different "truth" for the same wall — two exact descriptions of one plane, the defect
/// class the frame work exists to remove — and a sketch built on that lift would land on a
/// non-interned plane a hair off the face it means.
pub(super) fn realized_plane(origin: Point3, x_axis: Vector3, y_axis: Vector3) -> SketchPlane {
    SketchPlane {
        origin,
        x_axis,
        y_axis,
        def: None,
    }
}
