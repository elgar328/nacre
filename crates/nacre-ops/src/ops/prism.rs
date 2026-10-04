use super::*;
/// **The frame's world basis in exact rationals, when it has one** — the one gate that decides
/// whether a sketch is written in world coordinates or inside a motion node.
///
/// The frame's own basis comes from its plane's name ([`RatFrame::of_plane_frame`] on the chain's
/// narrow `Frame` head — `None` when an axis needs an irrational scale); a plane with a motion of
/// its own then carries it through that motion's fold ([`RatFrame::carried`] — translations,
/// quarter turns and axis reflections keep it rational; a frame node or a turn off the quarters
/// does not fold, and declines). A wide or judged head declines too. The answer is asked of the
/// truth, never of a realization lifted back: a realized axis that needs normalizing comes back
/// `0.6000000000000001` and would read as irrational.
///
/// ★ Declining is not a failure — it is the frame-node road, the same one a tilted face takes.
///
/// [`RatFrame::of_plane_frame`]: crate::construct::RatFrame::of_plane_frame
/// [`RatFrame::carried`]: crate::construct::RatFrame::carried
pub(super) fn exact_frame(
    model: &Model,
    frame: &SketchFrame,
) -> Option<crate::construct::RatFrame> {
    let chain =
        crate::rotated_vertex::frame_chain(model, frame.plane(), frame.placement(), frame.flip())?;
    let [nacre_judge::MoveNode::Frame { frame: pf }, ..] = chain.as_slice() else {
        return None;
    };
    let own = crate::construct::RatFrame::of_plane_frame(pf)?;
    match model.plane_motion(frame.plane()) {
        None => Some(own),
        Some(leaf) => own.carried(model, leaf),
    }
}

/// The builder's exact winding needs quarter-turn arcs (`Ring2d::winding_sign`); an arc of any
/// other angle is refused by name here, before the builder reads the ring.
pub(super) fn refuse_non_quarter_arcs(profile: &Profile2d) -> Result<(), OpError> {
    let zero = Rat::from_int(0);
    let quarter = |r: &Ring2d| -> Result<(), OpError> {
        let n = r.len();
        for i in 0..n {
            let Edge2d::Arc { center, .. } = r.edges()[i] else {
                continue;
            };
            let (s0, e0) = (r.vertices()[i], r.vertices()[(i + 1) % n]);
            if s0 == e0 {
                continue; // a whole circle
            }
            let v = |p: [Rat; 2]| -> Option<[Rat; 2]> {
                Some([p[0].checked_sub(center[0])?, p[1].checked_sub(center[1])?])
            };
            let (a, b) = (
                v(s0).ok_or(OpError::ProfileUndecidable)?,
                v(e0).ok_or(OpError::ProfileUndecidable)?,
            );
            let dot = a[0]
                .checked_mul(b[0])
                .and_then(|x| x.checked_add(a[1].checked_mul(b[1])?));
            let cross = a[0]
                .checked_mul(b[1])
                .and_then(|x| x.checked_sub(a[1].checked_mul(b[0])?));
            let (Some(dot), Some(cross)) = (dot, cross) else {
                return Err(OpError::ProfileUndecidable);
            };
            // A quarter, a half, three quarters: perpendicular radii, or opposite ones.
            if !((dot == zero && cross != zero) || (dot < zero && cross == zero)) {
                return Err(OpError::ArcSweepNotQuarterTurn);
            }
        }
        Ok(())
    };
    quarter(profile.outer())?;
    profile.holes().iter().try_for_each(quarter)
}

/// **A profile's rings placed in `frame` and swept `dist` along its `ŵ`** (against it when `dist`
/// is negative) — the one road from a frame to a prism's rings: the world road when
/// the frame's world basis is rational ([`exact_frame`]), otherwise written inside the frame
/// behind its node ([`push_frame_node`]).
///
/// The rational path is not an optimization: it is what makes `extrude(7.7)` and `extrude(1.1)`
/// then `extrude(6.6)` put their caps on the same plane rather than an ulp apart. Where the
/// arithmetic cannot state the prism (`i128` overflow in the placement) the answer is a **named
/// reject**, never a silent f64 prism: a point-less solid cannot state itself and cannot survive
/// a motion.
///
/// **Mapping only — no winding decision, and no containment check.** Forcing the outer ring CCW
/// here would be a second opinion on a question `build_prism` already answers from the sweep, and
/// two opinions is how an outer ring and its holes end up wound the same way (a negative extrude,
/// where the sweep runs `−n` and flips the outer ring). A profile may also reach past the face
/// boundary; an overhanging footprint routes to the overhang boolean sidecars, which reject
/// honestly what they do not cover.
pub(super) fn frame_rings(
    model: &mut Model,
    frame: &SketchFrame,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Swept, Vec<Swept>), OpError> {
    let (rat, node) = match exact_frame(model, frame) {
        Some(f) => (f, None),
        None => (
            crate::construct::RatFrame::identity(),
            Some(push_frame_node(model, *frame)),
        ),
    };
    crate::construct::prism_rings_in(model, rat, profile, dist, node)
        .ok_or(OpError::PlaneWithoutExactForm)
}

/// **Extrude a profile on a frame the model already holds** — the handle vocabulary of
/// [`Operation::Extrude`].
///
/// The base cap **is** the frame's plane, so it is handed to [`build_prism`] as a surface rather
/// than as points: the flush contact then reads as one shared handle, which is what the boolean
/// recognizes. Nothing new is pushed for it, whichever way the prism runs.
///
/// `dist` is signed: `ŵ` is the frame's, turned by whoever built the frame (a datum toward the
/// caller's stated normal, a face toward its outward), and a negative `dist` sweeps against it in
/// the same axes. The rings take the sign as it is ([`frame_rings`]); the sweep normal and the base
/// cap's facing turn with it.
pub(crate) fn extrude_on_frame(
    model: &mut Model,
    frame: &SketchFrame,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist == 0.0 {
        return Err(OpError::ZeroDistance);
    }
    if nacre_exact::Rat::from_decimal(dist).is_none() {
        return Err(OpError::DistOutsideDecimalWindow);
    }
    profile.check()?;
    refuse_non_quarter_arcs(profile)?;
    // ★ A `SketchFrame` may name any surface — `SketchFrame::canonical` makes no claim and checks
    // nothing — so a cylinder can reach here, which a `SketchPlane` never could. Reject it by name
    // rather than letting the frame derivation fail later for a reason that reads as something
    // else ("no exact form" when the truth is "not a plane").
    match model.surface(frame.plane()) {
        nacre_topo::Surface::Plane { .. } => {}
        nacre_topo::Surface::Cylinder { .. } => return Err(OpError::NonPlanarFace),
    }
    let (_, _, _, w) = crate::rotated_vertex::frame_world_basis(
        model,
        frame.plane(),
        frame.placement(),
        frame.flip(),
    )
    .ok_or(OpError::PlaneWithoutExactForm)?;

    let along = dist > 0.0;
    let (outer, holes) = frame_rings(model, frame, profile, dist)?;
    let base = (frame.plane(), base_cap_orientation(model, frame, along)?);
    let w = Vector3::from_array(w);
    build_prism(
        model,
        outer,
        holes,
        if along { w } else { -w },
        Some(base),
        None,
    )
}

/// **How the base cap of a prism swept on `frame` faces its plane** — the cap's outward is against
/// the sweep: `−ŵ` when the prism runs `along` `ŵ`, `+ŵ` when it runs against it — read against the
/// way the frame's plane `h` faces, from the truth.
///
/// `ŵ` is a difference of carried points, so the plane's motion carries it as `L(ŵ)`; the plane's
/// facing is its points' turn times `sense`, which the motion carries as `parity · L(·)`. In `h`'s
/// own frame `ŵ` runs with or against that facing
/// ([`crate::rotated_vertex::frame_normal_sense`]); `flip` negates `ŵ`, and so does a sweep
/// against it. So the cap faces the plane's way exactly when
/// `−(that relation) · parity · flip · sweep` is `+1`.
fn base_cap_orientation(
    model: &Model,
    frame: &SketchFrame,
    along: bool,
) -> Result<Orientation, OpError> {
    let h = frame.plane();
    let nacre_topo::Surface::Plane { motion, .. } = model.surface(h) else {
        return Err(OpError::NonPlanarFace);
    };
    let relation = crate::rotated_vertex::frame_normal_sense(model, h)
        .ok_or(OpError::PlaneWithoutExactForm)?;
    let parity = crate::rotated_vertex::motion_parity(model, *motion)
        .ok_or(OpError::PlaneWithoutExactForm)?;
    let flip = if frame.flip() { -1 } else { 1 };
    let sweep = if along { 1 } else { -1 };
    Ok(if -relation * parity * flip * sweep > 0 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}

/// Sweep a profile's rings along `sweep` into a prism solid: caps, side walls, and — for each
/// hole ring — a wall of its own plus an inner loop on each cap.
///
/// **This is the only place winding is decided**, because it is the only place that knows the
/// sweep. The outer ring is normalized CCW **about `sweep`** (area vector dotted with the sweep
/// normal — frame-independent, unlike `proj2`+`signed_area`, which mis-signs when the sweep runs
/// along a negative dominant axis, e.g. a pocket into a `+z` face sweeping `−z`). Each hole is
/// then normalized to the **opposite** sense *of that normalized outer ring*, never of the input:
/// get that backwards and a pocket — where the sweep flips the outer ring — silently produces
/// holes wound the same way as the outer, which is not a hole at all.
///
/// Opposite winding is all a hole needs. The wall quads are built from the ring's own traversal,
/// so a reversed ring yields walls facing into the hole (out of the material), and the cap loops
/// come out opposed to the cap's outer loop, which is what makes them holes.
///
/// No vertex carries a measured tolerance. Returns the solid and its faces: `faces[0]` = base cap
/// (at the ring, normal `−ŝ`), `faces[1]` = far cap, then the outer walls, then each hole's walls.
/// [`extrude_on_frame`]'s builder, either way along the frame's normal.
pub(crate) fn build_prism(
    model: &mut Model,
    outer_ring: Swept,
    inner_rings: Vec<Swept>,
    normal: Vector3,
    base_cap_surface: Option<(Handle<Surface>, Orientation)>,
    // ★ The world points of the plane the caller named, when they named one. Its canonical name is
    // derived from these, so there is no second half that could travel separately.
    base_cap_points: Option<[[nacre_exact::Rat; 3]; 3]>,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    // A polygon needs three corners; a ring with an arc bounds area with one (a circle) or two
    // (a half disk, a slot's end drawn alone).
    let thin = |r: &Swept| r.base.len() < 3 && !r.exact.segs.iter().any(Seg3::is_arc);
    if thin(&outer_ring) || inner_rings.iter().any(thin) {
        return Err(OpError::DegenerateProfile);
    }

    // Whether the sweep runs along the frame normal (a boss) or against it (a pocket): the rings
    // are normalized about the *sweep*, the arcs' `ccw` is stated about the *normal*, and the
    // lateral orientation rule below compares the two. Both are the ring's exact statement, in
    // its own frame — a motion carries them alike, so their dot's sign is the world's.
    let sweep_up = {
        let e = &outer_ring.exact;
        let d = nacre_exact::dot3_rat(&e.normal, &e.sweep().ok_or(OpError::DegenerateGeometry)?)
            .ok_or(OpError::DegenerateGeometry)?;
        d > nacre_exact::Rat::from_int(0)
    };
    let outer_pts = oriented_ring(outer_ring, sweep_up, true)?;
    let hole_pts: Vec<Swept> = inner_rings
        .into_iter()
        .map(|h| oriented_ring(h, sweep_up, false))
        .collect::<Result<_, _>>()?;

    // ★★★ **Surfaces before topology.** A vertex is defined by the three faces that meet at it,
    // and `Store` is append-only, so the handles have to exist before the vertex does. The two
    // arenas are separate, so interleaving them differently does not shift either one's numbering
    // — but the surfaces' order *among themselves* is what the numbering depends on, and it is
    // preserved exactly: base cap, then top cap, then walls, outer ring before the holes.

    // Base cap: outward normal −N, loops reversed.
    // When the prism stands on a plane the model holds, reuse its `Surface` handle (explicit sharing) so
    // the flush contact is a shared-handle coplanar pair the boolean can recognize by `Handle`
    // identity; otherwise push a fresh plane. The materialized outward normal must stay −N, and
    // the caller states the orientation that makes it so — it knows the shared plane's truth
    // (a face's own orientation, a frame's `flip`), where this would have to read the surface's
    // cache. `surface` and `orientation` travel together, and the reconstruction copies both.
    let (base_surface, base_orient) = match base_cap_surface {
        Some(shared) => shared,
        None => {
            // ★★★★ **The caller's statement wins here, and the frame's is the fallback.**
            //
            // The base cap *is* the plane they named, and they named it in the world — so stating
            // it that way keeps its `motion` `None`, keeps its judgment exact, and lets two
            // extrudes
            // on one plane share it **whatever frames they chose**. Writing it as `[0,0,1,0]` in
            // this prism's frame instead would be a second exact description of one plane, under a
            // different `SurfaceKey`.
            //
            // The frame's answer is what a plane with no caller statement gets (an axis-aligned
            // sketch, where the two agree anyway).
            //
            // ★★★★★ **One frame or the other, and the `motion` goes with the points.**
            // A plane's points and its `motion` are two halves of one statement: `None` means
            // *"these speak about the world"* and `Some` means *"about the pre-motion frame"*.
            // There is no third half: the name is derived from whichever points are recorded, so
            // the two cannot come apart. (Coefficients chosen by their own `or_else` could — world
            // coefficients beside the prism's **frame** ring, whenever the points overflowed and
            // the coefficients did not.)
            //
            // ★ The frame's `motion` is the one the top cap takes, and where it is `None` the
            // frame's points speak about the world too — which is every case a missing caller
            // statement can produce.
            let (base_motion, cap_pts) = match base_cap_points {
                Some(p) => (None, p),
                None => {
                    let e = &outer_pts.exact;
                    (
                        e.motion,
                        e.cap_points(false).ok_or(OpError::DegenerateGeometry)?,
                    )
                }
            };
            // The base cap faces against the sweep. A caller's world triple is read against a
            // world sweep, which is every ring such a caller builds.
            debug_assert!(
                base_cap_points.is_none() || outer_pts.exact.motion.is_none(),
                "a world base-cap statement beside a ring stated in a frame"
            );
            let sense = crate::rotated_vertex::sense_toward(
                model,
                cap_pts,
                outer_pts
                    .exact
                    .sweep_back()
                    .ok_or(OpError::DegenerateGeometry)?,
                base_motion,
            )
            .ok_or(OpError::DegenerateGeometry)?;
            let (s, flipped) = crate::realize::push_plane_realized(
                model,
                Plane::from_point_normal(outer_pts.base[0], -normal)
                    .ok_or(OpError::DegenerateGeometry)?,
                cap_pts,
                base_motion,
                sense,
            );
            // The plane was built with `−N` as its normal, so `Forward` is what states an outward
            // `−N` — unless a shared surface points the other way, which `flipped` reports.
            let orient = if flipped {
                Orientation::Forward.flipped()
            } else {
                Orientation::Forward
            };
            (s, orient)
        }
    };
    // Top cap: outward normal +N.
    let top_motion = outer_pts.exact.motion;
    let top_points = outer_pts
        .exact
        .cap_points(true)
        .ok_or(OpError::DegenerateGeometry)?;
    // The top cap faces along the sweep; `cap_points` fixes no sense of its own.
    let top_sense = crate::rotated_vertex::sense_toward(
        model,
        top_points,
        outer_pts.exact.sweep().ok_or(OpError::DegenerateGeometry)?,
        top_motion,
    )
    .ok_or(OpError::DegenerateGeometry)?;
    let (top_surface, top_flipped) = crate::realize::push_plane_realized(
        model,
        Plane::from_point_normal(outer_pts.top[0], normal).ok_or(OpError::DegenerateGeometry)?,
        top_points,
        top_motion,
        top_sense,
    );
    let top_orient = if top_flipped {
        Orientation::Forward.flipped()
    } else {
        Orientation::Forward
    };

    // Wall surfaces, in the same order the faces will be emitted: the outer ring's, then each
    // hole's (facing into the hole).
    let outer_walls = wall_surfaces(model, &outer_pts)?;
    let hole_walls: Vec<Vec<(Handle<Surface>, bool)>> = hole_pts
        .iter()
        .map(|h| wall_surfaces(model, h))
        .collect::<Result<_, _>>()?;

    // ── Topology. Every surface it needs already exists.
    let caps = (base_surface, top_surface);
    let outer = sweep_ring(model, &outer_pts, &outer_walls, caps)?;
    let holes: Vec<RingCells> = hole_pts
        .iter()
        .zip(hole_walls.iter())
        .map(|(h, w)| sweep_ring(model, h, w, caps))
        .collect::<Result<_, _>>()?;

    let mut faces =
        Vec::with_capacity(2 + outer.len() + holes.iter().map(|h| h.len()).sum::<usize>());
    faces.push(model.push_face(Face {
        surface: base_surface,
        outer: outer.cap_loop(Cap::Base),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Base)).collect(),
        orientation: base_orient,
    }));
    faces.push(model.push_face(Face {
        surface: top_surface,
        outer: outer.cap_loop(Cap::Top),
        inner: holes.iter().map(|h| h.cap_loop(Cap::Top)).collect(),
        orientation: top_orient,
    }));
    for (ring, walls, swept_pts) in std::iter::once((&outer, &outer_walls, &outer_pts)).chain(
        holes
            .iter()
            .zip(hole_walls.iter())
            .zip(hole_pts.iter())
            .map(|((r, w), p)| (r, w, p)),
    ) {
        ring.push_walls(model, walls, &swept_pts.exact.segs, sweep_up, &mut faces);
    }

    let shell = model.push_shell(Shell {
        faces: faces.clone(),
    });
    let solid = model.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    Ok((solid, faces))
}

/// Which cap a ring's loop is being built for.
#[derive(Clone, Copy)]
enum Cap {
    Base,
    Top,
}

/// One ring swept into cells: the two rings of vertices and the three edge families that join
/// them. Built the same way for the outer ring and for a hole — the difference is only which way
/// the ring runs, which the caller has already decided.
struct RingCells {
    /// Kept for [`RingCells::len`]; the geometry itself is read off the `Swept` these came from,
    /// which is also where the wall planes were built (`wall_surfaces`).
    base_pts: Vec<Point3>,
    be: Vec<Handle<Edge>>, // base  B_i -> B_{i+1}
    te: Vec<Handle<Edge>>, // top   T_i -> T_{i+1}
    /// Riser `B_i → T_i`; none for a whole circle, whose lateral is bounded by its two rims
    /// alone (`push_walls`).
    ve: Vec<Handle<Edge>>,
    /// Whether step `i`'s edge, walked in its own direction, follows the ring. A straight edge is
    /// pushed in ring order; a circle edge's own direction is fixed by the kernel's convention
    /// (`[A, B]` counter-clockwise about the axis — for a whole circle `A == B`, so the vertex
    /// order cannot carry it), so a clockwise arc's edge runs *against* the ring and every loop
    /// that walks it flips `forward`.
    along: Vec<bool>,
}

impl RingCells {
    fn len(&self) -> usize {
        self.base_pts.len()
    }

    /// The cap loop for this ring. The base cap faces `−ŝ`, so its loops run backwards.
    fn cap_loop(&self, cap: Cap) -> Loop {
        let n = self.len();
        match cap {
            Cap::Base => Loop {
                half_edges: (0..n)
                    .rev()
                    .map(|i| HalfEdge {
                        edge: self.be[i],
                        forward: !self.along[i],
                    })
                    .collect(),
            },
            Cap::Top => Loop {
                half_edges: (0..n)
                    .map(|i| HalfEdge {
                        edge: self.te[i],
                        forward: self.along[i],
                    })
                    .collect(),
            },
        }
    }

    /// One quad per ring segment. The quad's winding follows the ring's, so a ring wound against
    /// the outer one yields walls whose normals point into the hole.
    fn push_walls(
        &self,
        model: &mut Model,
        walls: &[(Handle<Surface>, bool)],
        segs: &[Seg3],
        sweep_up: bool,
        faces: &mut Vec<Handle<Face>>,
    ) {
        let n = self.len();
        for (i, &(surface, flipped)) in walls.iter().enumerate().take(n) {
            // A plane wall's sense is the plane's (`flipped` from `push_plane`). A cylinder wall
            // is `Forward` iff the material lies **inside** the cylinder. Seen from `+ŵ`, the
            // material is on the ring's left iff the sweep runs along `ŵ` (the rings are
            // normalized about the sweep — outer and holes alike, since a hole runs the other way
            // and the material is outside it), and the arc's centre is on its left iff the arc is
            // counter-clockwise about `ŵ`. Material inside ⟺ the two sides agree: a boss's convex
            // circle is `Forward`, a bore's circle and a notch's concave arc `Reversed`.
            let flipped = match &segs[i] {
                Seg3::Line => flipped,
                Seg3::Arc { ccw, .. } => *ccw != sweep_up,
            };
            let j = (i + 1) % n;
            // ★ **A whole circle's lateral is bounded by its two rims and nothing else** — the
            // base rim the outer loop, the top rim an inner one, each walked the way the quad
            // below walks it. A riser here would join a vertex to its own twin across the face,
            // an edge whose two sides are the one lateral: a seam is where the surface's
            // parametrization closes, not where the face has a boundary.
            if n == 1 {
                faces.push(model.push_face(Face {
                    surface,
                    outer: Loop {
                        half_edges: vec![HalfEdge {
                            edge: self.be[0],
                            forward: self.along[0],
                        }],
                    },
                    inner: vec![Loop {
                        half_edges: vec![HalfEdge {
                            edge: self.te[0],
                            forward: !self.along[0],
                        }],
                    }],
                    orientation: if flipped {
                        Orientation::Forward.flipped()
                    } else {
                        Orientation::Forward
                    },
                }));
                continue;
            }
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: self.be[i],
                        forward: self.along[i],
                    },
                    HalfEdge {
                        edge: self.ve[j],
                        forward: true,
                    },
                    HalfEdge {
                        edge: self.te[i],
                        forward: !self.along[i],
                    },
                    HalfEdge {
                        edge: self.ve[i],
                        forward: false,
                    },
                ],
            };
            faces.push(model.push_face(Face {
                surface,
                outer,
                inner: vec![],
                // The quad's winding is the ring's, which is what the plane above was built from;
                // a shared surface pointing the other way spells the same outward as `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            }));
        }
    }
}

/// One plane per ring segment, pushed **before** any of the ring's topology exists — see
/// `build_prism`.
fn wall_surfaces(model: &mut Model, ring: &Swept) -> Result<Vec<(Handle<Surface>, bool)>, OpError> {
    let n = ring.base.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            match &ring.exact.segs[i] {
                // The cache is the carried ring's own three points in `wall_points`' order, so it
                // faces the way those points do — through any chain, reflections included.
                Seg3::Line => Ok(crate::realize::push_plane_realized(
                    model,
                    Plane::through_points(ring.base[i], ring.base[j], ring.top[i])
                        .ok_or(OpError::DegenerateGeometry)?,
                    ring.exact.wall_points(i),
                    ring.exact.motion,
                    Orientation::Forward,
                )),
                // The wall is the cylinder about the arc's centre along the frame normal, stated
                // exactly and interned by that statement: every arc of one circle in one sketch
                // lands on one surface. The `bool` is a plane's flip; a cylinder's sense is
                // decided by `push_walls` from the arc's turn.
                Seg3::Arc {
                    center,
                    r2,
                    ref_dir,
                    cache,
                    ..
                } => {
                    let def = CylinderDef::new(
                        *center,
                        ring.exact.normal,
                        *ref_dir,
                        nacre_exact::BigRat::from(*r2),
                    )
                    .ok_or(OpError::DegenerateGeometry)?;
                    Ok((
                        crate::realize::push_cylinder_realized(
                            model,
                            *cache,
                            def,
                            ring.exact.motion,
                        ),
                        false,
                    ))
                }
            }
        })
        .collect()
}

/// The ring wound counter-clockwise about the sweep when `ccw`, clockwise when not. The test is
/// the ring's **exact winding** (`ring.exact.winding`), so it does not care which axis dominates;
/// the polygon's f64 area vector appears only inside a `debug_assert` that cross-checks it on an
/// arc-free ring.
/// **A zero winding is rejected** (`DegenerateProfile`): a ring that encloses nothing has no
/// side to pick. [`Profile2d::check`] is what keeps it from arising: a simple polygon cannot have
/// zero area, and a ring
/// that folds back on itself (a symmetric bowtie cancels to exactly zero) is not simple.
/// Normalize a ring's direction about the **sweep**: `ccw` for the outer ring, its opposite for a
/// hole. The winding is read exactly ([`crate::construct::SweptRat::winding`], about the ring's world
/// normal — the frame's motion, a reflection included, already folded in) and turned to the
/// sweep's sense by `sweep_up`; a ring with arcs has no f64 polygon area to read, and a polygon's
/// reads the same as before (`debug_assert`ed — the check that caught the mirrored pad).
fn oriented_ring(ring: Swept, sweep_up: bool, ccw: bool) -> Result<Swept, OpError> {
    let about_normal = ring.exact.winding;
    if about_normal == nacre_exact::Orient::Zero {
        return Err(OpError::DegenerateProfile);
    }
    let ccw_about_normal = about_normal == nacre_exact::Orient::Positive;
    debug_assert!(
        ring.exact.segs.iter().any(Seg3::is_arc) || {
            let v = &ring.base;
            let k = v.len();
            let area_vec = (0..k)
                .map(|i| (v[i] - Point3::origin()).cross(v[(i + 1) % k] - Point3::origin()))
                .fold(Vector3::from_array([0.0; 3]), |a, b| a + b);
            (area_vec.dot(ring.normal) > 0.0) == ccw_about_normal
        },
        "the exact winding agrees with the f64 polygon area"
    );
    let ccw_about_sweep = ccw_about_normal == sweep_up;
    Ok(if ccw_about_sweep != ccw {
        ring.reversed()
    } else {
        ring
    })
}

/// Push one ring's vertices and edges (base ring, top ring, risers).
///
/// The top ring arrives already computed rather than being derived here as
/// `base + sweep`: where the frame allows it that arithmetic is done in exact rationals
/// (see [`crate::construct`]), and a dimension split into two then lands on the same points
/// as the undivided one instead of an ulp away.
fn sweep_ring(
    model: &mut Model,
    ring: &Swept,
    walls: &[(Handle<Surface>, bool)],
    caps: (Handle<Surface>, Handle<Surface>),
) -> Result<RingCells, OpError> {
    let n = ring.base.len();
    let base_pts: Vec<Point3> = ring.base.clone();
    let top_pts: Vec<Point3> = ring.top.clone();
    // Corner `i` is where the wall before it, the wall after it, and the cap meet.
    //
    // ★ **Two walls that are one plane would name a line, not a point** — and since surfaces are
    // interned, "one plane" *is* "one handle", so the check is a comparison. This is a
    // guard on an invariant, not a live case: the only producer of adjacent same-plane walls is
    // a profile with a collinear midpoint, and `Profile2d`'s constructor dissolves those, so
    // every corner gets its three-plane definition (locked by
    // `a_collinear_midpoint_profile_builds_its_clean_twin_bit_for_bit`). Non-adjacent walls may
    // still legitimately share a plane (a notch), which never lands `prev == here`.
    // ★ Total: the guard's `None` (adjacent same-plane walls) is unreachable — the
    // profile constructor dissolves collinear midpoints — and a wall can never equal a cap (a
    // wall contains the sweep direction, a cap has it as normal). The honest reject stands in
    // for the unreachable arm; the suite is what would refute "unreachable".
    //
    // (There is no frame base vertex: a frame-drawn corner is the intersection of
    // three planes sharing the frame's motion, and solving them in that frame and replaying
    // the chain reproduces the stored coordinate bit for bit — measured 8/8, re-measured
    // suite-wide by the reuse differential.)
    let define = |model: &Model,
                  i: usize,
                  cap: Handle<Surface>,
                  at: &[nacre_exact::Rat; 3]|
     -> Result<Vertex, OpError> {
        let prev = walls[(i + n - 1) % n].0;
        let here = walls[i].0;
        if here == cap || prev == cap {
            return Err(OpError::DegenerateGeometry);
        }
        let (prev_arc, here_arc) = (
            ring.exact.segs[(i + n - 1) % n].is_arc(),
            ring.exact.segs[i].is_arc(),
        );
        match (prev_arc, here_arc) {
            // Two straight walls and the cap: the corner every polygon prism has.
            (false, false) => {
                if prev == here {
                    return Err(OpError::DegenerateGeometry);
                }
                Ok(Vertex::ThreePlane([prev, here, cap]))
            }
            // A whole circle: one wall, one vertex — the rim's point at `+ref_dir`, the seam.
            (true, true) if n == 1 => Ok(Vertex::OnSeam([here, cap])),
            // Two arcs of one circle merged in the profile's normal form; two of different
            // circles are a corner the kernel does not define.
            (true, true) => Err(OpError::ArcsMeetAtVertex),
            // A straight wall meets an arc wall on the cap: the wall's plane and the cap's plane
            // meet in a line that crosses the arc's cylinder there.
            (true, false) => pierce_def(model, here, cap, prev, at),
            (false, true) => pierce_def(model, prev, cap, here, at),
        }
    };
    let push_verts = |model: &mut Model,
                      ps: &[Point3],
                      exact: &[[nacre_exact::Rat; 3]],
                      cap: Handle<Surface>|
     -> Result<Vec<Handle<Vertex>>, OpError> {
        ps.iter()
            .zip(exact.iter())
            .enumerate()
            .map(|(i, (p, at))| {
                let def = define(model, i, cap, at)?;
                Ok(crate::realize::push_vertex_realized(
                    model,
                    def,
                    PointCache::Unrealized { coord: *p },
                    crate::realize::ChainLink::Fresh,
                ))
            })
            .collect()
    };
    let bv = push_verts(model, &base_pts, &ring.exact.base, caps.0)?;
    let tv = push_verts(model, &top_pts, &ring.exact.top, caps.1)?;

    let along: Vec<bool> = ring
        .exact
        .segs
        .iter()
        .map(|s| !matches!(s, Seg3::Arc { ccw: false, .. }))
        .collect();
    let (mut be, mut te, mut ve) = (Vec::new(), Vec::new(), Vec::new());
    // The walls and caps are each a carrier of several edges, so they are realized once.
    let mut memo = crate::realize::PlaneMemo::default();
    for i in 0..n {
        let j = (i + 1) % n;
        // Carriers: the same expression `define` uses for the corner triples — a base/top edge
        // runs between wall `i` and its cap, a riser between the walls either side of corner `i`.
        // A clockwise arc's edge is stored the other way round so its own direction — the
        // kernel's `[A, B]` counter-clockwise — names the same points; the loops then walk it
        // backwards (`along`).
        let (p, q) = if along[i] { (i, j) } else { (j, i) };
        be.push(push_line_edge(
            model,
            bv[p],
            bv[q],
            [walls[i].0, caps.0],
            &mut memo,
        )?);
        te.push(push_line_edge(
            model,
            tv[p],
            tv[q],
            [walls[i].0, caps.1],
            &mut memo,
        )?);
        // A whole circle has no riser: corner `0` is between its one wall and itself.
        if n > 1 {
            ve.push(push_line_edge(
                model,
                bv[i],
                tv[i],
                [walls[(i + n - 1) % n].0, walls[i].0],
                &mut memo,
            )?);
        }
    }
    Ok(RingCells {
        base_pts,
        be,
        te,
        ve,
        along,
    })
}

/// **Does the gate read a frame the way its author wrote it?** — the measurement [`exact_frame`]
/// rests on.
///
/// A caller writes a frame as `f64` axes; their decimals lifted (`SketchPlane::exact`, test-only)
/// are what the caller meant. The operation asks the datum's frame in rationals instead
/// ([`exact_frame`]), and the two must answer alike — while the route that realizes the axes and
/// lifts them back must visibly lose the frames whose axes need normalizing.
///
/// ★★★ **A one-ulp disagreement there is not an ulp of error.** The plane would silently take the
/// frame-node road, and the arena would gain motion nodes and write its points in frame
/// coordinates: a *different but still valid* model. So the assertion order below matters —
/// **same road first**, values second. This crate has been bitten by exactly
/// this shape before: `nacre_exact::plane_frame_named` records `v̂` realized as `ŵ × û` coming out
/// `0.999999999999999_7`, "an exact path quietly lost".
///
/// The prediction is agreement, and it is structural rather than lucky: a datum's `ref_dir` is
/// `points[1] − points[0]`, which *is* the caller's `+u`; for a unit rational axis `|u_raw|² = 1`
/// so `inv_sqrt_exact` returns exactly one; and `v_raw` has its own exact form. But an argument is
/// not a gate.
#[cfg(test)]
#[path = "../tests/ops_frame_road.rs"]
mod frame_road;
