use super::*;
/// **The frame's basis in exact rationals, when it has one** — the gate that decides whether a
/// sketch is written in world coordinates or inside a motion node.
///
/// Only a frame on a plane with **no motion of its own** can be rational in the world: a moved
/// plane's axes carry the motion's irrational part, which is precisely why that road exists. So a
/// chain of one narrow frame node is the whole population, and anything longer (or wide) declines.
///
/// ★ Declining is not a failure — it is the frame-node road, the same one a tilted face takes.
pub(super) fn exact_frame(model: &Model, frame: &SketchFrame) -> Option<crate::exact::RatFrame> {
    let chain =
        crate::rotated_vertex::frame_chain(model, frame.plane(), frame.placement(), frame.flip())?;
    match chain.as_slice() {
        [nacre_cip::MoveNode::Frame { frame: pf }] => crate::exact::RatFrame::of_plane_frame(pf),
        _ => None,
    }
}

/// **Extrude a profile on a frame the model already holds** — the handle vocabulary of
/// [`Operation::Extrude`], and the same three steps `extrude_and_boolean` takes for a pad.
///
/// The base cap **is** the frame's plane, so it is handed to [`build_prism`] as a surface rather
/// than as points: the flush contact then reads as one shared handle, which is what the boolean
/// recognizes. Nothing new is pushed for it.
///
/// `dist > 0` is a **thickness**; which way it goes is the frame's `ŵ`, measured by whoever built
/// the frame (a datum against the caller's stated normal, a face against its outward). That is why
/// this can take a frame where the operation used to take a plane and sweep the same way.
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

pub(crate) fn extrude_on_frame(
    model: &mut Model,
    frame: &SketchFrame,
    profile: &Profile2d,
    dist: f64,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    if dist <= 0.0 {
        return Err(OpError::NonPositiveDistance);
    }
    if nacre_scalar::Rat::from_decimal(dist).is_none() {
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

    // World road when the frame's basis is rational, frame-node road otherwise — the same
    // question `SketchPlane::exact` asks a caller's axes, asked of the frame in `Rat`.
    let (rat, node) = match exact_frame(model, frame) {
        Some(f) => (f, None),
        None => (
            crate::exact::RatFrame::identity(),
            Some(push_frame_node(model, *frame)),
        ),
    };
    let (outer, holes) = crate::exact::prism_rings_in(model, rat, profile, dist, node)
        .ok_or(OpError::PlaneWithoutExactForm)?;
    build_prism(
        model,
        outer,
        holes,
        Vector3::from_array(w),
        Some(frame.plane()),
        None,
    )
}

/// Place a profile on its plane and sweep it — **exactly, or not at all** (S6b).
///
/// The rational path is not an optimization: it is what makes `extrude(7.7)` and
/// `extrude(1.1)` then `extrude(6.6)` put their caps on the same plane rather than an
/// ulp apart. Where it does not apply — a plane with no exact statement (axes outside the
/// decimal window, a degenerate pair), a frame the chain cannot realize, an i128 overflow in
/// the placement arithmetic — the answer is a **named reject**, not the silent f64 prism this
/// used to build: a point-less solid cannot state itself, cannot survive a motion, and is the
/// population `Inexact` grew from. (Profile coordinates and `dist` outside the decimal window
/// are named before this runs.)
///
/// **Mapping only — no winding decision, and no containment check.** Forcing the outer
/// ring CCW here would be a second opinion on a question `build_prism` already answers
/// from the sweep, and two opinions is how an outer ring and its holes end up wound the
/// same way (the pocket case, where the sweep runs `−n` and flips the outer ring). A
/// profile may also reach past the face boundary; an overhanging footprint routes to
/// the overhang boolean sidecars, which reject honestly what they do not cover.
pub(super) fn swept_profile(
    model: &Model,
    plane: &SketchPlane,
    profile: &Profile2d,
    dist: f64,
    frame: Option<Handle<MotionNode>>,
) -> Result<(Swept, Vec<Swept>), OpError> {
    crate::exact::prism_rings(model, plane, profile, dist, frame)
        .ok_or(OpError::PlaneWithoutExactForm)
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
/// Shared by [`extrude_on_frame`] (a boss) and the pocket (`sweep = −n`).
pub(crate) fn build_prism(
    model: &mut Model,
    outer_ring: Swept,
    inner_rings: Vec<Swept>,
    normal: Vector3,
    base_cap_surface: Option<Handle<Surface>>,
    // ★ The world points of the plane the caller named, when they named one. Its canonical name is
    // derived from these, so there is no second half that could travel separately.
    base_cap_points: Option<[[nacre_scalar::Rat; 3]; 3]>,
) -> Result<(Handle<Solid>, Vec<Handle<Face>>), OpError> {
    // A polygon needs three corners; a ring with an arc bounds area with one (a circle) or two
    // (a half disk, a slot's end drawn alone).
    let thin = |r: &Swept| r.base.len() < 3 && !r.exact.segs.iter().any(Seg3::is_arc);
    if thin(&outer_ring) || inner_rings.iter().any(thin) {
        return Err(OpError::DegenerateProfile);
    }

    // Whether the sweep runs along the frame normal (a boss) or against it (a pocket): the rings
    // are normalized about the *sweep*, the arcs' `ccw` is stated about the *normal*, and the
    // lateral orientation rule below compares the two.
    let sweep_up = outer_ring.normal.dot(normal) > 0.0;
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
    // When padding/pocketing on a face, reuse that face's `Surface` handle (explicit sharing) so
    // the flush contact is a shared-handle coplanar pair the boolean can recognize by `Handle`
    // identity; otherwise push a fresh plane. The materialized outward normal must stay −N, so the
    // face orientation is chosen from the shared surface's stored normal — `surface` and
    // `orientation` travel together, and the reconstruction copies both.
    let (base_surface, base_orient) = match base_cap_surface {
        Some(h) => {
            let n_h = match model.surface_cache(h) {
                nacre_geom::Surface::Plane(p) => p.normal(),
                nacre_geom::Surface::Cylinder(_) => return Err(OpError::DegenerateGeometry),
            };
            let orient = if n_h.dot(-normal) > 0.0 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            (h, orient)
        }
        None => {
            // ★★★★ **The caller's statement wins here, and the frame's is the fallback.**
            //
            // The base cap *is* the plane they named, and they named it in the world — so stating
            // it that way keeps it `Constructed`, keeps its judgment exact, and lets two extrudes
            // on one plane share it **whatever frames they chose**. Writing it as `[0,0,1,0]` in
            // this prism's frame instead would be a second exact description of one plane, under a
            // different `SurfaceKey` — the very duplication this work removes.
            //
            // The frame's answer is what a plane with no caller statement gets (an axis-aligned
            // sketch, where the two agree anyway).
            //
            // ★★★★★ **One frame or the other, and the `def` goes with the points.**
            // A plane's points and its `SurfaceDef` are two halves of one statement: `Constructed`
            // means *"these speak about the world"* and `Moved` means *"about the pre-motion
            // frame"*. There used to be a third half — coefficients, chosen by their own `or_else`,
            // so a caller whose points overflowed while their coefficients did not got world
            // coefficients beside the prism's **frame** ring. That half no longer exists: the name
            // is derived from whichever points are recorded, so the two can no longer come apart.
            //
            // ★ The frame's `surface_def()` is the same one the top cap takes, and it agrees with
            // `Constructed` wherever there is no motion — which is every case a missing caller
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
            let (s, flipped) = model.push_plane(
                Plane::from_point_normal(outer_pts.base[0], -normal)
                    .ok_or(OpError::DegenerateGeometry)?,
                cap_pts,
                base_motion,
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
    let (top_surface, top_flipped) = model.push_plane(
        Plane::from_point_normal(outer_pts.top[0], normal).ok_or(OpError::DegenerateGeometry)?,
        top_points,
        top_motion,
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
    ve: Vec<Handle<Edge>>, // riser B_i -> T_i
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
/// `build_prism`. The same three points `push_walls` used to build them from, and in the same
/// order, so the surface arena's numbering is untouched.
fn wall_surfaces(model: &mut Model, ring: &Swept) -> Result<Vec<(Handle<Surface>, bool)>, OpError> {
    let n = ring.base.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            match &ring.exact.segs[i] {
                Seg3::Line => Ok(model.push_plane(
                    Plane::through_points(ring.base[i], ring.base[j], ring.top[i])
                        .ok_or(OpError::DegenerateGeometry)?,
                    ring.exact.wall_points(i),
                    ring.exact.motion,
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
                        nacre_scalar::BigRat::from(*r2),
                    )
                    .ok_or(OpError::DegenerateGeometry)?;
                    Ok((model.push_cylinder(*cache, def, ring.exact.motion), false))
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
/// hole. The winding is read exactly ([`crate::exact::SweptRat::winding`], about the ring's world
/// normal — the frame's motion, a reflection included, already folded in) and turned to the
/// sweep's sense by `sweep_up`; a ring with arcs has no f64 polygon area to read, and a polygon's
/// reads the same as before (`debug_assert`ed — the check that caught the mirrored pad).
fn oriented_ring(ring: Swept, sweep_up: bool, ccw: bool) -> Result<Swept, OpError> {
    let about_normal = ring.exact.winding;
    if about_normal == nacre_scalar::Orient::Zero {
        return Err(OpError::DegenerateProfile);
    }
    let ccw_about_normal = about_normal == nacre_scalar::Orient::Positive;
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
/// (see [`crate::exact`]), and a dimension split into two then lands on the same points
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
                  at: &[nacre_scalar::Rat; 3]|
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
                      exact: &[[nacre_scalar::Rat; 3]],
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
    for i in 0..n {
        let j = (i + 1) % n;
        // Carriers: the same expression `define` uses for the corner triples — a base/top edge
        // runs between wall `i` and its cap, a riser between the walls either side of corner `i`.
        // A clockwise arc's edge is stored the other way round so its own direction — the
        // kernel's `[A, B]` counter-clockwise — names the same points; the loops then walk it
        // backwards (`along`).
        let (p, q) = if along[i] { (i, j) } else { (j, i) };
        be.push(push_line_edge(model, bv[p], bv[q], [walls[i].0, caps.0])?);
        te.push(push_line_edge(model, tv[p], tv[q], [walls[i].0, caps.1])?);
        ve.push(push_line_edge(
            model,
            bv[i],
            tv[i],
            [walls[(i + n - 1) % n].0, walls[i].0],
        )?);
    }
    Ok(RingCells {
        base_pts,
        be,
        te,
        ve,
        along,
    })
}
