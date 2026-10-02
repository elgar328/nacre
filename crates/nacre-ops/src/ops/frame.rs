use super::*;
/// A [`SketchFrame`] whose `ŵ` faces `toward` — the plane's own facing when `Forward`, against it
/// when `Reversed` (a face's orientation, or which way a datum's statement faces against the
/// handle the door returned). Every road calls it, and its `flip` is [`flip_toward`]'s. `None`
/// when the plane has no frame (no name and no judged frame).
///
/// ★★★ **Only the frame comes back — deliberately.** A caller that needs the axes asks
/// [`crate::rotated_vertex::frame_world_basis`] *with this flip*, so the geometry of `flip` is
/// written in exactly one place. A flip turns a `Canonical` frame about `v̂` and a `Named` one
/// about `û` ([`crate::rotated_vertex::narrow_frame`]), so combining an unflipped basis with
/// `flip` **by hand** by any one rule is wrong on one of the two (measured when a report
/// half-turned about `v̂` and the realization about `û`: a footprint centred on the face through
/// the reported coordinates landed outside it — the pad missed its face — on 2 of 6 faces of a
/// turned block).
pub(super) fn frame_toward(
    model: &Model,
    plane: Handle<Surface>,
    placement: nacre_topo::FramePlacement,
    toward: Orientation,
) -> Option<SketchFrame> {
    crate::rotated_vertex::frame_chain(model, plane, &placement, false)?;
    let flip = flip_toward(model, plane, toward)?;
    Some(SketchFrame {
        plane,
        placement,
        flip,
    })
}

/// **Whether a frame on `plane` must be flipped for its `ŵ` to face `toward`** — the one place
/// `flip` is decided, so no two roads can decide differently. It does not depend on the
/// placement: which way `ŵ` runs is the plane's, and the placement only moves `û` within it.
///
/// `flip` is stated from the truth: `ŵ` runs with the plane's facing in its own frame by
/// [`crate::rotated_vertex::frame_normal_sense`], and the plane's motion carries `ŵ` as `L(ŵ)`
/// and the facing as `parity · L(·)`. So `ŵ` faces `toward` in the world exactly when
/// `sense · parity · toward` is `+1`. A realized `ŵ` dotted against an `f64` direction would ask
/// the same question of two roundings.
fn flip_toward(model: &Model, plane: Handle<Surface>, toward: Orientation) -> Option<bool> {
    let sense = crate::rotated_vertex::frame_normal_sense(model, plane)?;
    let parity = crate::rotated_vertex::motion_parity(model, model.plane_motion(plane))?;
    Some(sense * parity * toward.sign() < 0)
}

/// Name a [`SketchFrame`] as the [`nacre_topo::Motion::Frame`] node the sweep writes coordinates
/// against — the one road from the frame value to a node, shared by the extrude and face paths.
///
/// ★★ **`push_motion` interns**, so two sketches in one frame name the *same* node — which is
/// what makes their surfaces intern too (a plane's `SurfaceKey::Name` is `(name, motion)`): two
/// routes to one height become one `Handle<Surface>` at construction, with no f64 comparison
/// anywhere. With `Canonical` placement the node is `(plane, Canonical, flip)` — nothing
/// per-sketch in the key.
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
/// * **Axis-aligned faces keep axes in `{0, ±1}`**, exactly the rationals
///   [`nacre_exact::plane_frame_default`] states for a face — the same rule in exact arithmetic,
///   which is how a face's frame is chosen (`face_sketch_frame`).
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
/// Returns both axes rather than just `u`, so a caller cannot run `v` the other way. `None` only
/// for the zero vector.
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

/// The sketch plane of a planar face — **the very frame an [`Operation::Extrude`] on
/// [`face_sketch_frame`] places its profile in**, so a caller can work out where its `(0, 0)` will
/// land before it builds anything.
///
/// That equality is the contract, not a coincidence: this is [`face_sketch_frame`]'s frame
/// realized by [`frame_plane`], never a second derivation. A test pins a hand-placed profile
/// against a boss built in that frame to keep it that way.
///
/// `NonPlanarFace` for a curved surface (only a plane carries a frame); `FaceNotInLiveSolid` if no
/// live solid's outer shell holds the face; `PlaneWithoutExactForm` if its plane has no frame.
pub fn face_plane(model: &Model, face: Handle<Face>) -> Result<SketchPlane, OpError> {
    let frame = face_sketch_frame(model, face)?;
    frame_plane(model, &frame).ok_or(OpError::PlaneWithoutExactForm)
}

/// **Where a frame is, in space** — its origin and its two axes, realized as f64.
///
/// [`face_plane`] answers this for a face; this answers it for any frame a caller holds, which
/// is what a viewer needs to draw a sketch where it was drawn: the sketch's coordinates are
/// `(u, v)` in this frame, and `origin + u·x + v·y` is the point.
///
/// ★ It is a **report**, not a truth. The frame's statement is the truth — this is that
/// statement realized the way a prism built in it is realized (`frame_basis`), and nothing
/// exact should be decided from it.
///
/// `None` when the frame's chain cannot be realized at all (a plane with no name). A caller
/// that cannot place a thing should decline to draw it rather than draw it somewhere wrong.
pub fn frame_plane(model: &Model, frame: &SketchFrame) -> Option<SketchPlane> {
    let (o, u, v, _) = frame_basis(model, frame)?;
    Some(realized_plane(
        Point3::from_array(o),
        Vector3::from_array(u),
        Vector3::from_array(v),
    ))
}

/// **A frame realized the way a prism built in it is** — `(origin, û, v̂, ŵ)`.
///
/// The world road builds a prism's vertices in the frame's rational world basis and rounds each
/// once, so that basis rounded once is its realization; the frame-node road replays its vertices
/// through the chain, so the chain replayed is its realization
/// ([`crate::rotated_vertex::frame_world_basis`]). The two part in the last bit where an axis
/// needs normalizing — a replay multiplies by a numerically computed `1/5` and lands on
/// `0.6000000000000001` — so each road's report is its own road's realization.
pub(super) fn frame_basis(
    model: &Model,
    frame: &SketchFrame,
) -> Option<crate::rotated_vertex::WorldBasis> {
    match exact_frame(model, frame) {
        Some(rf) => rf.realized(),
        None => crate::rotated_vertex::frame_world_basis(
            model,
            frame.plane(),
            frame.placement(),
            frame.flip(),
        ),
    }
}

/// A planar face's sketch frame **as a [`SketchFrame`]** — the plane handle, placement, and flip a
/// feature on the face sketches in: an [`Operation::Extrude`] in it builds the boss (`dist > 0`)
/// or the pocket's tool (`dist < 0`) on the face's own plane handle; [`face_plane`] is its
/// realization.
///
/// ★★ **Which frame a face takes is a fact about where the face is, not about how its plane is
/// stored.** Whether a moved plane carries its motion as a node or had it carried into its points
/// is the transform's choice of what it can state exactly, so a rule that read the
/// representation would sketch two identical faces in two frames.
/// - A plane whose world equation the model states ([`Model::world_plane_name`] — unmoved, or
///   moved by a chain that folds) takes the **arbitrary-axis frame of its outward normal in the
///   world** (`frame_axes`' rule, stated exactly by [`nacre_exact::plane_frame_default`]): the
///   world origin's projection, `+u = ẑ × n` (`ŷ × n` when `n` is vertical). It is spelled as a
///   placement of the plane's own statement, carried back through the chain — which folds, so
///   exactly — and as `Canonical` wherever that is the same frame (the normal form). Through a
///   reflection the carried frame is left-handed: `û` is the world's, `v̂` its reverse, `ŵ`
///   still outward.
/// - Any other plane (a frame node or a turn off the quarters in its chain, a wide name, no name)
///   takes its own `Canonical` frame, carried out by its chain.
///
/// Errors: `NonPlanarFace`, `FaceNotInLiveSolid`; `PlaneWithoutExactForm` when the plane has no
/// frame at all (no name, and no judged frame).
pub fn face_sketch_frame(model: &Model, face: Handle<Face>) -> Result<SketchFrame, OpError> {
    model
        .live_solids()
        .iter()
        // Index-only equality again: a face handle from another model can match here. The
        // shell lookup that follows is where the cross-store guard fires.
        .find(|&&s| model.shell(model.solid(s).outer).faces.contains(&face))
        .ok_or(OpError::FaceNotInLiveSolid)?;
    let f = model.face(face);
    let surface_h = f.surface;
    let outward = f.orientation;
    match model.surface(surface_h) {
        nacre_topo::Surface::Plane { .. } => {}
        nacre_topo::Surface::Cylinder { .. } => return Err(OpError::NonPlanarFace),
    }
    let frame = world_placement(model, surface_h, outward)
        .and_then(|placement| frame_toward(model, surface_h, placement, outward))
        .or_else(|| {
            frame_toward(
                model,
                surface_h,
                nacre_topo::FramePlacement::Canonical,
                outward,
            )
        })
        .ok_or(OpError::PlaneWithoutExactForm)?;
    Ok(frame)
}

/// **The placement that puts a face's sketch in the arbitrary-axis frame of its outward normal
/// in the world**, written in the plane's own statement — `None` when the model does not state the
/// plane's world equation (a chain that does not fold, a wide name) or the arithmetic overflows.
///
/// The world equation, turned to face the face's outward, gives the frame's origin and `+u` in
/// the world ([`nacre_exact::plane_frame_default`] — the origin is the world origin's projection,
/// a property of the plane, so two faces of one plane agree about where `(0, 0)` is). A
/// placement speaks in the coordinates the plane's points are written in, so a moved plane takes
/// them carried back through its chain; the chain folds, so that is exact. Where the result is
/// the frame the plane's own name derives under the same `flip`, it is `Canonical` — asked of
/// [`crate::rotated_vertex::narrow_frame`], the derivation the frame will actually take, so the
/// normal form cannot name a `Canonical` that realizes elsewhere.
fn world_placement(
    model: &Model,
    plane: Handle<Surface>,
    outward: Orientation,
) -> Option<nacre_topo::FramePlacement> {
    let world = *model.world_plane_name(plane)?.narrow()?;
    let zero = Rat::from_int(0);
    let world = if model.world_plane_name_sense(plane)?.sign() * outward.sign() > 0 {
        world
    } else {
        [
            zero.checked_sub(world[0])?,
            zero.checked_sub(world[1])?,
            zero.checked_sub(world[2])?,
            zero.checked_sub(world[3])?,
        ]
    };
    let (origin, ref_dir) = nacre_exact::plane_frame_default(world)?;
    let (origin, ref_dir) = match model.plane_motion(plane) {
        None => (origin, ref_dir),
        Some(leaf) => (
            model.chain_point_rat_inverse(leaf, origin)?,
            model.chain_dir_rat_inverse(leaf, ref_dir)?,
        ),
    };
    let own = *model.surface_name.get(&plane)?.narrow()?;
    let flip = flip_toward(model, plane, outward)?;
    let named = nacre_topo::FramePlacement::Named { origin, ref_dir };
    let canonical = crate::rotated_vertex::narrow_frame(own, &named, flip)?
        == crate::rotated_vertex::narrow_frame(own, &nacre_topo::FramePlacement::Canonical, flip)?;
    Some(if canonical {
        nacre_topo::FramePlacement::Canonical
    } else {
        named
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
