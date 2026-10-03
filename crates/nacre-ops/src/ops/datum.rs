use super::*;
/// **The definition of a corner where a straight wall meets an arc wall on a cap**: the wall's
/// plane and the cap's plane meet in a line, and that line crosses the arc's cylinder at this
/// point — [`Vertex::Pierce`], the same definition the boolean mints for its pierce corners.
///
/// The root is read the way `QuadRoot` is defined: the two planes in **ascending handle order**,
/// each by its **stored canonical name** (`surface_name`, the sign convention that fixes the meet
/// line's direction), fed to [`nacre_exact::quad::plane_plane_cylinder`] with the cylinder's own
/// statement; the root whose parameter equals this point's is the name. A wall tangent to the
/// cylinder — a fillet's, a slot's straight side — is the double root, one point. Every statement
/// has to live in one frame for the meet to mean anything: a cap borrowed from another body in
/// another frame (a pad on a turned face) declines by name rather than reading a frame line
/// against a world cylinder.
pub(super) fn pierce_def(
    model: &Model,
    plane: Handle<Surface>,
    cap: Handle<Surface>,
    cylinder: Handle<Surface>,
    at: &[nacre_exact::Rat; 3],
) -> Result<Vertex, OpError> {
    use nacre_exact::quad::CylinderMeet;
    use nacre_topo::QuadRoot;
    let (a, b) = if plane.index() < cap.index() {
        (plane, cap)
    } else {
        (cap, plane)
    };
    let name = |h: Handle<Surface>| -> Result<[nacre_exact::Rat; 4], OpError> {
        model
            .surface_name
            .get(&h)
            .and_then(|n| n.narrow().copied())
            .ok_or(OpError::PlaneWithoutExactForm)
    };
    let (pa, pb) = (name(a)?, name(b)?);
    let nacre_topo::Surface::Cylinder { def, motion } = model.surface(cylinder) else {
        return Err(OpError::DegenerateGeometry);
    };
    if model.plane_motion(a) != *motion || model.plane_motion(b) != *motion {
        return Err(OpError::PlaneWithoutExactForm);
    }
    let meet =
        nacre_exact::quad::plane_plane_cylinder(&pa, &pb, &def.origin(), &def.dir(), def.r2())
            .ok_or(OpError::PlaneWithoutExactForm)?;
    // This point's parameter along the meet line: `(at − base)·dir / dir·dir`.
    let param = |line: &nacre_exact::quad::MeetLine| -> Option<nacre_exact::Rat> {
        let (bse, d) = (line.base(), line.dir());
        let mut num = nacre_exact::Rat::from_int(0);
        let mut den = nacre_exact::Rat::from_int(0);
        for k in 0..3 {
            num = num.checked_add(at[k].checked_sub(bse[k])?.checked_mul(d[k])?)?;
            den = den.checked_add(d[k].checked_mul(d[k])?)?;
        }
        num.checked_mul(nacre_exact::Rat::new(den.denom(), den.numer())?)
    };
    let root = match meet {
        CylinderMeet::Tangent { line, s } => {
            debug_assert_eq!(param(&line), Some(s), "the tangent point is this corner");
            QuadRoot::Double
        }
        CylinderMeet::Pair { line, s } => {
            let t = param(&line).ok_or(OpError::PlaneWithoutExactForm)?;
            // Compared as values: a rational root still arrives as `mid ± k·√disc` when the
            // discriminant is a perfect square, so `b == 0` is not the test — the difference's
            // sign is.
            let is = |q: &nacre_exact::quad::QuadVal| {
                q.checked_sub(&nacre_exact::quad::QuadVal::from_rat(t))
                    .is_some_and(|d| d.sign() == nacre_exact::Orient::Zero)
            };
            if is(&s[0]) {
                QuadRoot::Lo
            } else if is(&s[1]) {
                QuadRoot::Hi
            } else {
                return Err(OpError::PlaneWithoutExactForm);
            }
        }
        _ => return Err(OpError::DegenerateGeometry),
    };
    Ok(Vertex::Pierce {
        planes: [a, b],
        cylinder,
        root,
    })
}

pub(super) fn push_line_edge(
    model: &mut Model,
    a: Handle<Vertex>,
    b: Handle<Vertex>,
    carriers: [Handle<Surface>; 2],
    memo: &mut crate::realize::PlaneMemo,
) -> Result<Handle<Edge>, OpError> {
    crate::realize::push_edge_realized(model, carriers, [a, b], memo)
        .map_err(|_| OpError::DegenerateGeometry)
}

/// Put a stated plane in the arena and hand back its handle and the frame the statement implies.
///
/// ★★ **The cache faces `−normal`, and that is a convention with two independent reasons.**
/// (i) It is the sense a base cap gets: `extrude` pushes `−plane.normal()` and `Model::new` seeds
/// the world planes along `−axis`; measured, seeding `+axis` instead flips 781 stored
/// cap normals for nothing. (ii) `WorkingPlane::frame_sign` records whether a plane's *stored* normal
/// agrees with its root face's outward normal, and a base cap's outward is `−N` — so `−normal`
/// leaves that sign exactly where it is today. A datum that later becomes a base cap therefore
/// interns with `flipped == false` and nothing downstream has to compensate.
///
/// Which point anchors the cache is a choice the arena keeps (interning discards the newcomer's
/// cache), and `tests/instruments/plane_anchor.rs` measures what that choice costs: on the tilted
/// population
/// where anchors can disagree at all, the judgment path does not read the disagreeing part and the
/// result holds to 1.1e-15 at the worst anchor.
pub(super) fn datum_plane(
    model: &mut Model,
    def: &DatumDef,
) -> Result<(Handle<Surface>, SketchFrame), OpError> {
    /// The same solve `Model::through_points_rat` performs, with its single `None` split into the
    /// causes that have different owners.
    ///
    /// ★ The split is the point. `three_planes_rat`'s one `None` hid two different facts for as
    /// long as it existed, and the measurement that found it is the reason this variant is being
    /// added at all — repeating the shape here would make the next stage unable to read how much
    /// of the population it is opening.
    ///
    /// ★★ A carrier triple that does not meet is **not** on the list: those three planes met, or
    /// the vertex would not exist. That is asserted, not rejected, so a broken invariant cannot
    /// arrive disguised as a user error.
    ///
    /// ★★★★★ **The motion comes back with the points, because it *is* which frame they are
    /// written in.** Returning the points alone is what produced the defect this function was
    /// rewritten for: the caller then chose a motion, chose `None`, and a plane stated in a frame
    /// was filed as a world plane. `DatumDef::Offset` has always returned `(points, motion)` from
    /// one decision for exactly this reason.
    /// What the three vertices state, once the causes are told apart.
    ///
    /// (`Named` is 288 B beside a unit variant — the same shape and the same verdict as
    /// `PlanePoints`: it lives for one call on one stack frame, and boxing would buy nothing.)
    #[allow(clippy::large_enum_variant)]
    enum ThroughStatement {
        /// One shared frame — the named road: the vertices' meets in that frame
        /// (**any width** — a meet wider than `Rat` still names its plane
        /// through `plane_name_from_meets`), and the frame.
        Named(
            [nacre_exact::MeetPoint; 3],
            Option<Handle<nacre_topo::MotionNode>>,
        ),
        /// ★ No one frame holds all three: either a vertex the door cannot place at all (a
        /// straddle) or three placeable vertices whose frames differ. No
        /// rational triple, no name, and the plane takes the judged road.
        Nameless,
    }

    #[allow(clippy::type_complexity)]
    fn through_points_by_cause(
        model: &Model,
        vs: [Handle<Vertex>; 3],
    ) -> Result<ThroughStatement, OpError> {
        // ★★★★ **Which frame a vertex is solvable in is [`nacre_topo::Model::vertex_meet`]'s
        // answer, not a second copy of its rule.** Comparing `plane_motion(tri[0])` against the
        // other two here would call a turned solid's corner a straddle — the invariant-plane
        // restatement leaves its cap world-stated while the walls carry a node — and, worse,
        // could *disagree* with the door: `push_plane_through` derives the interning name through
        // `vertex_meet` and asserts the statement's motion carries that frame, so a producer that
        // says "no frame" while the door says "this one" stops there (it once filed a frame-local
        // name as a world plane — `a_datum_through_frame_local_vertices_is_not_a_world_plane`).
        // One decision, one place.
        let mut pts: [Option<nacre_exact::MeetPoint>; 3] = [None, None, None];
        let mut frames = [None; 3];
        for (i, vh) in vs.iter().enumerate() {
            let tri = match *model.vertex(*vh) {
                nacre_topo::Vertex::ThreePlane(tri) => tri,
                // A through-vertices datum needs three-plane meets; a seam vertex has no
                // point-meet at all and a pierce point has no rational one — the same honest
                // reject, spelled per variant.
                nacre_topo::Vertex::OnSeam(_) | nacre_topo::Vertex::Pierce { .. } => {
                    return Err(OpError::VertexNotThreePlane);
                }
            };
            if tri.iter().any(|h| !model.surface_name.contains_key(h)) {
                // ★★ The one population still refused here: a carrier that is itself a nameless
                // datum (a datum on a nameless datum — depth). Its triangle needs the judged
                // machinery recursively, which is the depth question. This is
                // the **only** thing `VerticesInMixedFrames` names; measure before splitting
                // the label further.
                return Err(OpError::VerticesInMixedFrames);
            }
            let Some((meet, frame)) = model.vertex_meet(*vh) else {
                // ★ A straddling vertex — no rational coordinate anywhere, but a complete
                // definition (the meet of its carriers). The judged road takes it, as
                // long as every carrier can hand out a witness triangle of its own.
                //
                // ★★★ **"The carriers meet" is not assertable here, and finding that out cost a
                // red run.** The old pass two `expect`ed it, and rightly: it ran only after the
                // frames agreed. Solving three names *stated in different frames* is not the
                // same question — they are not three planes in one coordinate system, so
                // `three_planes_big` declines for a straddling vertex as a matter of course. An
                // assertion here would measure the proposition next door. The invariant now
                // lives where it is meaningful: past this `else`, every `pts[i]` is `Some`, so
                // the named road below has the meets by construction rather than by `expect`.
                return Ok(ThroughStatement::Nameless);
            };
            pts[i] = Some(meet);
            frames[i] = Some(frame);
        }
        // ★★ Pure vertices in differing frames — the judged road. The collinearity question is
        // *not* asked there: with no shared frame there is no exact solve to ask it in, and the
        // judged constructor's failure is reported as undecided, never as proven collinear.
        if !(frames[1] == frames[0] && frames[2] == frames[0]) {
            return Ok(ThroughStatement::Nameless);
        }
        // ★ The meets are kept at whatever width they need — the name is a
        // function of the *plane*, and `plane_name_from_meets` derives it without ever asking a
        // coordinate to fit `Rat`.
        let meets = pts.map(|p| p.expect("filled above"));
        if nacre_exact::plane_name_from_meets([&meets[0], &meets[1], &meets[2]]).is_none() {
            return Err(OpError::CollinearVertices);
        }
        Ok(ThroughStatement::Named(meets, frames[0].flatten()))
    }

    /// Which way the caller's stated normal faces against the handle the door returned. Every arm
    /// pushes its statement facing `−stated` (the base cap's outward), so it is that handle's
    /// facing reversed — unless the door handed back a handle already facing the other way
    /// (`flipped`: two statements of one plane compared, the door's own answer).
    fn stated_side(flipped: bool) -> Orientation {
        if flipped {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }

    match def {
        DatumDef::Stated(sp) => {
            let d = sp.def.as_ref().ok_or(OpError::PlaneWithoutExactForm)?;
            let cache = Plane::from_point_normal(sp.origin(), -sp.normal())
                .ok_or(OpError::DegenerateGeometry)?;
            // Every `PlaneDef` orders its points so `u × v` is the stated normal; the cache faces
            // the other way (the base cap's outward).
            let (plane, flipped) = crate::realize::push_plane_realized(
                model,
                cache,
                d.points(),
                None,
                Orientation::Reversed,
            );
            // ★ `Named`, unconditionally — never derived. The canonical frame of the ZX plane has
            // `+u = −x̂` while the script convention (and `SketchPlane::world_zx`) says `+ẑ`, so a
            // placement inferred from the plane would silently turn some sketches. The values are
            // the `PlaneDef`'s own rationals, and `SketchFrame::named`'s checks hold structurally:
            // `origin = points[0]` is on the plane, and `ref_dir = points[1] − points[0]` lies in
            // it and is nonzero because the triple is not collinear.
            let placement = nacre_topo::FramePlacement::Named {
                origin: d.origin(),
                ref_dir: d.ref_dir(),
            };
            // ★★★ **`flip` makes the frame face the normal the caller stated** — the only place
            // the caller's *direction* can survive.
            //
            // A plane's canonical name has no direction, and planes intern: state `z = 0` facing
            // `+ẑ` and state it facing `−ẑ`, and both come back as **one handle whose `ŵ` is `+ẑ`**
            // (measured). So a frame built with `flip: false` would silently answer "up" to a
            // caller who said "down". Facing the direction their own point order fixes makes the
            // returned frame mean what they said, which is what lets an operation take a frame
            // rather than a plane and sweep the same way.
            let frame = frame_toward(model, plane, placement, stated_side(flipped))
                .ok_or(OpError::PlaneWithoutExactForm)?;
            Ok((plane, frame))
        }
        DatumDef::ThroughVertices(vs) => {
            // ★★★ **Every reject here happens before anything is pushed.** The causes are told
            // apart first, then the plane is built — so a refusal leaves the model exactly as it
            // found it, and the name of the refusal says which stage owns it.
            let mut sorted = *vs;
            sorted.sort_by_key(|v| v.index());
            if sorted[0] == sorted[1] || sorted[1] == sorted[2] {
                return Err(OpError::DuplicateVertex);
            }
            let statement = through_points_by_cause(model, sorted)?;

            // ★ **The caller's order is the stated normal**, by the right-hand rule — the one
            // place their choice of side can survive, since the stored triple is sorted and a
            // canonical name carries no direction. The pushed statement faces against it
            // (`through_sense` below, exactly), so the frame faces it as the `Stated` arm's does.
            //
            // `stated` itself is the vertices' world caches crossed — rounded, and read only to
            // point the plane's cache, which is a cache's use.
            let world = vs.map(|v| model.vertex_point(v));
            let stated = (world[1] - world[0])
                .cross(world[2] - world[0])
                .normalize()
                .ok_or(OpError::CollinearVertices)?;
            // `stated` is the caller's order; the plane stores the sorted order, which spans
            // `stated` times the permutation's sign — and the cache faces `−stated`. Exact: a
            // permutation's parity is a count, not a measurement.
            let inversions = (0..3)
                .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
                .filter(|&(i, j)| vs[i].index() > vs[j].index())
                .count();
            let through_sense = if inversions % 2 == 0 {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };

            let ThroughStatement::Named(meets, motion) = statement else {
                // ★★★ **The judged road**: every vertex pure, frames
                // differing — no name exists, and the frame is derived from the defining points
                // as intervals. **Validation comes before the push**: the judged constructor can
                // refuse (an undecidable basis), and rejecting after `push_plane_through` would
                // leave a frameless nameless plane in the arena — the reject-after-commit shape
                // this arm's own header forbids. The shared derivation
                // (`through_judged_points`) is the one `frame_chain` will re-run, so what was
                // validated is what gets framed.
                let jpts = crate::rotated_vertex::through_judged_points(model, sorted)
                    .ok_or(OpError::PlaneWithoutExactForm)?; // unreachable: causes told apart above
                let ft = nacre_judge::FrameThrough::of(jpts, false)
                    .ok_or(OpError::ThroughFrameUndecided)?;
                // The cache anchors at the validated realization of the first stored vertex —
                // the definition's own replay, same rule as the named road below.
                let anchor =
                    Point3::from_array(ft.anchor_coord().ok_or(OpError::ThroughFrameUndecided)?);
                let cache =
                    Plane::from_point_normal(anchor, -stated).ok_or(OpError::DegenerateGeometry)?;
                let (plane, flipped) = model.push_plane_through(cache, sorted, None, through_sense);
                let frame = frame_toward(
                    model,
                    plane,
                    nacre_topo::FramePlacement::Canonical,
                    stated_side(flipped),
                )
                .ok_or(OpError::PlaneWithoutExactForm)?;
                return Ok((plane, frame));
            };

            // ★★★ **The cache is a world description, and the meets are not world coordinates
            // unless the carriers share no motion.** Realizing them means walking the very chain
            // the definition names — `construct.rs` states the rule ("the realization must be the
            // definition's own replay, not a second route to the same real number"), and taking
            // any other road here is how a plane's cache and its truth end up describing
            // different planes.
            //
            // ★ **Split by width so the narrow bits stay put**: a `Narrow` meet takes the road
            // letter for letter. A `Wide` meet has no
            // `Rat` triple to replay — and no f64 realization of it survives an exact lift
            // either (a value like 5⁻⁴⁰ rounds to a dyadic whose denominator leaves `i128`) —
            // so its cache anchors at the **first stored vertex's own world cache**: the
            // definition's replay already performed by whoever made the vertex (the 8/8 lock is
            // what says re-solving and replaying lands on it), and a rounded cache like every
            // anchor (`plane_anchor.rs` measured what anchor wobble costs — nothing the judging
            // reads).
            let anchor = match (&meets[0], motion) {
                (nacre_exact::MeetPoint::Narrow(p), None) => {
                    Point3::from_array(p.map(|r| r.to_f64()))
                }
                (nacre_exact::MeetPoint::Narrow(p), Some(leaf)) => {
                    let chain = crate::rotated_vertex::motion_chain(model, leaf)
                        .ok_or(OpError::PlaneWithoutExactForm)?;
                    let w =
                        crate::rotated_vertex::replay(nacre_judge::WitnessPoint::at(*p), &chain)
                            .ok_or(OpError::PlaneWithoutExactForm)?;
                    Point3::from_array(w.coord())
                }
                (nacre_exact::MeetPoint::Wide(_), _) => model.vertex_point(sorted[0]),
            };
            let cache =
                Plane::from_point_normal(anchor, -stated).ok_or(OpError::DegenerateGeometry)?;
            let (plane, flipped) = model.push_plane_through(cache, sorted, motion, through_sense);
            let frame = frame_toward(
                model,
                plane,
                nacre_topo::FramePlacement::Canonical,
                stated_side(flipped),
            )
            .ok_or(OpError::PlaneWithoutExactForm)?;
            Ok((plane, frame))
        }
        DatumDef::Offset { frame, dist } => {
            // ★ Zero is the one offset that would duplicate a plane the model already holds; see
            // `OpError::ZeroOffset`. Checked before the lift so the reject names the real fault.
            if *dist == 0.0 {
                return Err(OpError::ZeroOffset);
            }
            let d =
                nacre_exact::Rat::from_decimal(*dist).ok_or(OpError::DistOutsideDecimalWindow)?;
            let base = frame.plane;

            // ★★ **Normalize to the plane's canonical frame, folding `flip` into the sign.** A
            // parallel plane is fixed by `(plane, signed distance)` — a placement's origin and
            // `+u` do not move it — so building in the caller's frame would give one geometric
            // plane as many handles as there are frames naming it. Every frame on one plane has
            // the same `ŵ` before its `flip` — the narrow, wide and judged roads all take it from
            // the plane, never from the placement — so the caller's side is their `flip`. (That
            // their frame realizes at all is the constructors' promise: `SketchFrame::named`
            // builds the frame before it hands one out.)
            let canonical = nacre_topo::FramePlacement::Canonical;
            let cb = crate::rotated_vertex::frame_world_basis(model, base, &canonical, false)
                .ok_or(OpError::PlaneWithoutExactForm)?;
            let d = if frame.flip() {
                nacre_exact::Rat::from_int(0)
                    .checked_sub(d)
                    .ok_or(OpError::DistOutsideDecimalWindow)?
            } else {
                d
            };

            // ★★ **Say it in the world when the world can hold it** — the same node-omission
            // normalization frames use. A plane stated under a frame node lives at the
            // key `(name, Some(node))`, so an offset of the world XY plane would *not* intern with
            // a box's cap on the same plane. Where the canonical frame has a rational world basis
            // (`exact_frame` — the one gate the extrude road asks, asked of the truth), the offset
            // plane is rational in the world and is stated there.
            let canonical_frame = SketchFrame {
                plane: base,
                placement: canonical,
                flip: false,
            };
            // The base plane's handedness: a frame carried through a reflection spans
            // `x × y = parity · ŵ`, on either road.
            let parity = nacre_judge::chain_parity(
                &crate::rotated_vertex::frame_chain(model, base, &canonical, false)
                    .ok_or(OpError::PlaneWithoutExactForm)?,
            );
            let (points, motion) = match exact_frame(model, &canonical_frame) {
                // ★ An overflowing pullback is a **named reject**, not a quiet switch to the frame
                // road — that switch is exactly where the duplicate handle would appear.
                // `rf`'s normal is the carried `ŵ`, a reflection included, so `d` along it lands on
                // the side the caller asked for.
                Some(rf) => (
                    rf.offset_plane_points(d)
                        .ok_or(OpError::PlaneWithoutExactForm)?,
                    None,
                ),
                None => {
                    let zero = nacre_exact::Rat::from_int(0);
                    let one = nacre_exact::Rat::from_int(1);
                    let node = push_frame_node(model, canonical_frame);
                    (
                        [[zero, zero, d], [one, zero, d], [zero, one, d]],
                        Some(node),
                    )
                }
            };

            // The cache rides the canonical realization: anchor `o + d·ŵ`, facing `−ŵ` — the same
            // `−normal` convention `Stated` and every base cap keep.
            //
            // ★ **`d` here is the *signed* distance, not the caller's `dist`.** Using the raw one
            // puts the cache on the far side of the plane its own truth names whenever the
            // caller's frame is flipped — an incoherent surface, and one that a `flip` positive
            // control caught rather than any amount of reading.
            let signed = d.to_f64();
            let anchor = Point3::from_array(core::array::from_fn(|k| cb.0[k] + signed * cb.3[k]));
            let w = Vector3::from_array(cb.3);
            let cache = Plane::from_point_normal(anchor, -w).ok_or(OpError::DegenerateGeometry)?;
            // Both roads' points span `x × y = parity · ŵ` in the world, and the cache faces `−ŵ`.
            let sense = if parity == 1 {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };
            let (plane, _flipped) =
                crate::realize::push_plane_realized(model, cache, points, motion, sense);
            Ok((plane, SketchFrame::canonical(plane)))
        }
    }
}
