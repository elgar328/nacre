use super::*;
/// The supporting planes of a solid's outer shell. `Rejected` if any face is
/// non-planar or lacks three non-collinear loop points.
pub(crate) fn collect_planes(
    model: &Model,
    solid: Handle<Solid>,
) -> Result<Vec<FaceRow>, BoolError> {
    let mut out = Vec::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shell(sh).faces {
            let face = model.face(fh);
            let plane = match model.surface_cache(face.surface) {
                nacre_geom::Surface::Plane(p) => *p,
                // A cylinder face sits in the table — its row keeps the shared facts
                // (surface, face, stated outward sign, motion leaf) and none of the plane
                // vocabulary. Whether it may *flow* is the population gate's question, asked
                // in `arrangement::plane_index_setup`, not a door slam here.
                nacre_geom::Surface::Cylinder(_) => {
                    let nacre_topo::Surface::Cylinder { motion, .. } = model.surface(face.surface)
                    else {
                        unreachable!(
                            "push_cylinder_raw pairs them, so a cylinder cache has a cylinder truth"
                        )
                    };
                    // ★ **The world statement, not the stated frame's** — a translated cylinder's
                    // truth is written before its motion, and every consumer of this row (the
                    // gate's clearance arithmetic, the arrangement's circles and rulings) compares
                    // it against *world* planes. A row whose statement cannot be carried out to
                    // the world keeps none: the gate then refuses by its own name rather than
                    // measuring across two frames.
                    let def = model.world_cylinder_def(face.surface);
                    out.push(FaceRow::Cylinder(CylFaceInfo {
                        surf: face.surface,
                        face: Some(fh),
                        orient_sign: face.orientation.sign(),
                        // ★ **The frame this row's description stands in**, which is the world
                        // once the statement has been carried out to it — not the provenance of
                        // the surface. The one reader is the restatement mirror below, whose
                        // question is exactly "which frames are in play here".
                        motion: motion.filter(|_| def.is_none()),
                        footprint: Footprint {
                            span: def.as_ref().and_then(|d| lateral_t_range(model, face, d)),
                            theta: def
                                .as_ref()
                                .and_then(|d| lateral_theta_extent(model, face, d)),
                        },
                        def,
                    }));
                    continue;
                }
            };
            // **Outward is read off the face's statement, not re-derived from its loop.**
            // `orientation` relates the stored surface normal to "out of the solid" — the
            // same reading props and STEP already trust — so `n_out` is that product,
            // exact in direction by construction. A triangle's cross is not the source: on
            // rotated near-collinear corners the direction it gives back is the rounding (the
            // pad-eats-material defect). The winding is consulted — as the cross-check below,
            // not as the answer.
            let orient_sign = face.orientation.sign();
            let n_out = plane.normal() * f64::from(orient_sign);
            let tri = match outer_tri(model, face) {
                Some((tri, _)) => tri,
                // ★ **A face bounded by a circle-curve edge has area whatever its vertices do.**
                // A disk cap's outer loop is one circle edge with a single seam vertex; a
                // half-disk's is an arc and its chord with two. Neither spreads three loop
                // points, and neither is degenerate — the plane's own truth points state the
                // triangle instead (the cap's construction points, realized), wound to this
                // face's stated outward like every other `tri`. A loop with no curved edge that
                // spreads nothing is still `DegenerateFace`: it bounds no area.
                None if face.outer.half_edges.iter().any(|he| {
                    matches!(model.edge_curve(he.edge), nacre_geom::Curve::Circle(_))
                }) =>
                {
                    let nacre_topo::Surface::Plane {
                        points: nacre_topo::PlanePoints::Known(pts),
                        motion: disk_motion,
                        ..
                    } = model.surface(face.surface)
                    else {
                        return Err(reject(RejectReason::DegenerateFace));
                    };
                    // ★ The points are stated in the frame the plane's motion names — a disk cap
                    // on a slanted wall's frame is in that frame's coordinates. Realized through
                    // the motion, as every other witness here is, before being wound to the
                    // face's *world* outward; naive f64 of the frame-stated points compared a
                    // frame triangle against a world normal and asserts on the first such face
                    // (a circle padded on a slanted wall).
                    let mut tri = match disk_motion {
                        None => pts.map(|p| Point3::from_array(p.map(|x| x.to_f64()))),
                        Some(m) => {
                            let chain = crate::rotated_vertex::motion_chain(model, *m)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                            let mut out = [Point3::origin(); 3];
                            for (o, p) in out.iter_mut().zip(pts.iter()) {
                                let q = crate::rotated_vertex::replay(WitnessPoint::at(*p), &chain)
                                    .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                                *o = Point3::from_array(q.coord());
                            }
                            out
                        }
                    };
                    let wound = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
                    if wound.dot(n_out) < 0.0 {
                        tri.swap(1, 2);
                    }
                    tri
                }
                None => return Err(reject(RejectReason::DegenerateFace)),
            };
            // **The plane's exact definition comes from the surface, not from the vertices.**
            //
            // Both are decided per surface, not per *solid* ("is this solid rotated?"): a
            // boolean's result carries no rotation provenance, so asking the solid would describe
            // every result face by its rounded triangle, and one wall would become two plane
            // classes on the next operation. The surface knows (its truth — points + motion),
            // and a result face reuses its operand's surface handle, so the answer
            // survives a chain of booleans.
            //
            // `tri_pt3` is an *oriented* plane witness, but the recorded triple belongs to the
            // *plane* — two faces sharing it can face opposite ways, and the implicit-point
            // `orient3d` reads the side `tri_pt3` spans. So both arms wind it to agree with
            // *this* face's `n_out`.
            let wind = |mut w: [WitnessPoint; 3]| -> [WitnessPoint; 3] {
                let e1 = Vector3::from_array(w[1].coord()) - Vector3::from_array(w[0].coord());
                let e2 = Vector3::from_array(w[2].coord()) - Vector3::from_array(w[0].coord());
                if e1.cross(e2).dot(n_out) < 0.0 {
                    w.swap(1, 2);
                }
                w
            };
            let (tri_pt3, rotated, motion) = match model.surface(face.surface) {
                // ★★★★★ **The plane's own points state the plane — not the face's triangle.**
                //
                // After a chain of booleans a face's corners are seam points the kernel itself
                // annotated with a bound, and stating them as tol-0 witnesses claims an exactness
                // they lack — measured, 1,938 of 68,350 triangles carried one, in 162 plane
                // tables, 77 of which also held a rotated plane, where the claim is consulted.
                // The plane's recorded points are construction points, so no such claim is made.
                //
                // ★ `WitnessPoint::at_nearest` states a representable point with tol 0 and a
                // non-representable one with the ½-ulp bound
                // `Rat::to_f64` promises — one spelling for a point that came from a rational.
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    motion: None,
                    ..
                } => {
                    // ★★★ **The bound is free; measuring it is not** — `at_nearest` states the
                    // rounding from `Rat::to_f64`'s contract (measured here once: 38.3% of these
                    // points are not f64, and measuring them at 120 bits cost the suite 6%).
                    let w = pts.map(WitnessPoint::at_nearest);
                    (wind(w), false, None)
                }
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    motion: Some(motion),
                    ..
                } => {
                    let motion = *motion;
                    let _t = Watch::new(); // charged at the arm's end
                    // The pre-motion description, carried through the recorded chain — the same
                    // computation, in the same order, that a moved vertex's `WitnessPoint` performs.
                    let chain = crate::rotated_vertex::motion_chain(model, motion)
                        .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                    let turn = |base: [nacre_exact::Rat; 3]| -> Result<WitnessPoint, BoolError> {
                        crate::rotated_vertex::replay(WitnessPoint::at(base), &chain)
                            .ok_or_else(|| reject(RejectReason::FrameOutOfRange))
                    };
                    // ★★★★★ **The exact points the plane was built from — the only description.**
                    //
                    // No realized `witness: [Point3; 3]` rides beside the motion: a third of
                    // realized points do not survive the trip (measured 34.5% — `11/10` is not an
                    // f64, so lifting the realization back recovers a different rational), and a
                    // plane described through three rounded points tells two caps that are one
                    // plane apart with full confidence.
                    let w = [turn(pts[0])?, turn(pts[1])?, turn(pts[2])?];
                    _t.charge(Sub::TriPt3);
                    (wind(w), true, Some(motion))
                }
                // ★★★ **A `Through` plane is solved into the same witness triangle here.**
                //
                // Its truth is handles, and the judging layer takes points — so the points are
                // *derived* at the boundary, once per plane per operation, exactly like the name
                // was derived once at push. That does not make the points a second truth: this
                // table lives for one operation and is a mirror.
                //
                // The rational-closure branch: three vertices that solve to `Rat` give a
                // triangle indistinguishable from a stated one, so every predicate below
                // runs unchanged.
                //
                // ★★★ **The nameless branch (the heterogeneous half).**
                // A pure-mixed datum's vertices are each exact *in their own frame*, so its
                // witness triangle exists — three `WitnessPoint`s whose chains simply differ.
                // The judging layer never required them to agree: `plane_iv`/`plane_hp` realize
                // each point independently, so the whole toleranced route runs unchanged,
                // and `standard_for` sizes the operation from these very points. What such a
                // plane has none of is exact f64 coefficients — so it is flagged `rotated`
                // (the routing signal means "no exact description", not "carries a motion"),
                // which makes every exact shortcut decline.
                //
                // A *straddling*-vertex datum has no witness triangle **from its vertices** —
                // its witness is the judged frame's probes, built in the branch below.
                nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Through(vs),
                    motion,
                    ..
                } => {
                    let motion = *motion;
                    let _t = Watch::new();
                    let replay_all = |mut w: [WitnessPoint; 3],
                                      m: Handle<nacre_topo::MotionNode>|
                     -> Result<[WitnessPoint; 3], BoolError> {
                        let chain = crate::rotated_vertex::motion_chain(model, m)
                            .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                        for o in w.iter_mut() {
                            *o = crate::rotated_vertex::replay(o.clone(), &chain)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                        }
                        Ok(w)
                    };
                    let (w, rotated) = match model.through_points_rat(*vs) {
                        Some(base) => {
                            let w = base.map(WitnessPoint::at);
                            let w = match motion {
                                None => w,
                                Some(m) => replay_all(w, m)?,
                            };
                            // ★ `rotated` is also the licence to read the name as a **world**
                            // description (`PlaneWitness::exact_coeffs`). A `Through` name is
                            // derived from the meets in the frame the three vertices share, which
                            // the plane's `motion` is meant to name; with no motion recorded and
                            // the vertices meeting in a pre-motion frame, the name speaks that
                            // frame. Asked here rather than trusted — the producers keep it, and
                            // nothing else checks it.
                            let frame_local = motion.is_none()
                                && vs.iter().any(|&v| {
                                    model.vertex_meet(v).is_some_and(|(_, f)| f.is_some())
                                });
                            (w, motion.is_some() || frame_local)
                        }
                        None => {
                            let j = crate::rotated_vertex::through_judged_points(model, *vs)
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                            // ★ All-pure unwraps to the heterogeneous triangle, letter for
                            // letter — that road is locked by its own tests and stays.
                            //
                            // ★★★★ **Any implicit point: the witness is the plane's own judged
                            // frame**. The table's contract is "three exact points *on the
                            // plane*, wound to n_out" — never "the face's corners" — and a judged plane has such points by definition: its
                            // canonical frame's probes `(0,0,0)·(1,0,0)·(0,1,0)`, the same ones
                            // `frame_world_basis` realizes. The origin is the foot of the
                            // perpendicular (on the plane exactly), û and v̂ are in-plane by
                            // construction, and `frame_chain` already appends the plane's own
                            // later motion — so the probes are exact definitions the escalation
                            // realizes at any precision.
                            let all_pure = j
                                .iter()
                                .all(|p| matches!(p, nacre_judge::JudgedPoint::Pure(_)));
                            if all_pure {
                                let mut w: [Option<WitnessPoint>; 3] = [None, None, None];
                                for (o, jp) in w.iter_mut().zip(j) {
                                    if let nacre_judge::JudgedPoint::Pure(wp) = jp {
                                        *o = Some(wp);
                                    }
                                }
                                let w = w.map(|o| o.expect("all pure"));
                                let w = match motion {
                                    // ★ The plane's own later motion appends to each point's
                                    // own chain — motions add, never multiply.
                                    None => w,
                                    Some(m) => replay_all(w, m)?,
                                };
                                (w, true)
                            } else {
                                let chain = crate::rotated_vertex::frame_chain(
                                    model,
                                    face.surface,
                                    &nacre_topo::FramePlacement::Canonical,
                                    false,
                                )
                                .ok_or_else(|| reject(RejectReason::FrameOutOfRange))?;
                                let r = nacre_exact::Rat::from_int;
                                let probe = |u: i128, v: i128| {
                                    crate::rotated_vertex::replay(
                                        WitnessPoint::at([r(u), r(v), r(0)]),
                                        &chain,
                                    )
                                    .ok_or_else(|| reject(RejectReason::FrameOutOfRange))
                                };
                                // `frame_chain` already carries the plane's own later motion in
                                // its tail, so no `replay_all` here — appending it twice would
                                // move the witness off the plane.
                                let w = [probe(0, 0)?, probe(1, 0)?, probe(0, 1)?];
                                (w, true)
                            }
                        }
                    };
                    _t.charge(Sub::TriPt3);
                    (wind(w), rotated, motion)
                }
                nacre_topo::Surface::Cylinder { .. } => {
                    unreachable!("a plane cache cannot carry a cylinder truth")
                }
            };
            // The loop's winding, cross-checked against the flag it stated. `validate`
            // pins this same invariant as `FaceMisoriented`; here it runs at every
            // boolean in debug. A failure is a producer bug — a lying flag or a
            // mis-wound loop — never a conditioning artifact (`outer_tri` picks the
            // widest corner).
            //
            // ★ This assertion needs three points, so it says nothing about a **disk**
            // face (one closed rim) or about a face's **holes**. Those belong to the
            // same invariant and `validate` owns them — it reads a rim's circle and
            // every inner loop. Replicating that here would be a second spelling of one
            // rule, which is the shape this kernel keeps having to undo.
            debug_assert!(
                (tri[1] - tri[0])
                    .cross(tri[2] - tri[0])
                    .normalize()
                    .is_some_and(|w| w.dot(n_out) > 0.5),
                "a face's outer winding must agree with its stated orientation: face {fh:?} on {:?}, {} outer half-edges, {} inner loops, orientation {:?}",
                face.surface,
                face.outer.half_edges.len(),
                face.inner.len(),
                face.orientation
            );
            let name = model.surface_name.get(&face.surface).cloned();
            // ★★ **The witness must span a plane** — three collinear points lie on every plane
            // through their line, and `Judge::planes_coplanar` would read them as "on" whatever it
            // asked. A name is the proof: it is derived from exactly these points (a motion keeps
            // them non-collinear), and `plane_name_exact` is `None` only for collinear ones. The
            // other arms have theirs — a judged frame's probes are `(0,0)·(1,0)·(0,1)`, and a
            // `Through` datum over collinear vertices is refused where it is stated. What is
            // left is a nameless `Known` plane: collinear in production, or a test's raw push
            // (`push_plane_unregistered` skips the name), so only there is the question asked.
            if name.is_none() {
                if let nacre_topo::Surface::Plane {
                    points: nacre_topo::PlanePoints::Known(pts),
                    ..
                } = model.surface(face.surface)
                {
                    if nacre_exact::plane_name_exact(pts[0], pts[1], pts[2]).is_none() {
                        return Err(reject(RejectReason::DegenerateFace));
                    }
                }
            }
            out.push(FaceRow::Plane(FaceInfo {
                base_rat: name.as_ref().and_then(|n| n.narrow()).copied(),
                world_name: model.world_plane_name(face.surface),
                name,
                surf: face.surface,
                face: Some(fh),
                plane,
                tri,
                n_out,
                orient_sign,
                tri_pt3,
                rotated,
                motion,
            }));
        }
    }
    // ★★★ **The mirror takes the strongest description — a restated plane rejoins its
    // solid's chain here.** The truth keeps a motion-fixed plane world-stated (the
    // invariant-plane restatement: one plane, one handle, no node), but this table is a
    // per-operation mirror, and a chain-borne description is strictly stronger when the
    // rest of the solid carries one: with every plane on one chain the shared-motion
    // roads answer exactly — a shared rotation must not turn exact questions into
    // assumed ones (`a_shared_rotation_still_assumes_nothing`). A fixed plane admits it
    // verbatim: its world points replayed through the chain still lie on the plane (the
    // chain maps the plane onto itself), and its world coefficients ARE its pre-motion
    // coefficients — bit-identical to the description the moved twin carried before the
    // restatement. Row-verbatim fixedness (`chain_preserves_plane_row`), not mere
    // set-fixedness: a sign-flipped row would poison determinant reads.
    //
    // The chain is chosen deterministically (lowest node index that qualifies) so replay
    // reproduces the table; a plane no chain fixes keeps its world description.
    let leaves = {
        let mut v: Vec<Handle<nacre_topo::MotionNode>> = out
            .iter()
            .filter_map(|r| match r {
                FaceRow::Plane(f) => f.motion,
                FaceRow::Cylinder(c) => c.motion,
            })
            .collect();
        v.sort_unstable_by_key(|h| h.index());
        v.dedup();
        v
    };
    if !leaves.is_empty() {
        for row in out.iter_mut() {
            // The restatement mirror is a plane story — a motion-fixed *cylinder*'s
            // restatement is separately deferred.
            let FaceRow::Plane(f) = row else { continue };
            if f.motion.is_some() {
                continue;
            }
            let Some(name) = f.base_rat else { continue };
            let nacre_topo::Surface::Plane {
                points: nacre_topo::PlanePoints::Known(pts),
                motion: None,
                ..
            } = model.surface(f.surf)
            else {
                continue;
            };
            let Some(&c) = leaves
                .iter()
                .find(|&&c| model.chain_preserves_plane_row(c, &name))
            else {
                continue;
            };
            let Some(chain) = crate::rotated_vertex::motion_chain(model, c) else {
                continue;
            };
            let turned: Option<Vec<WitnessPoint>> = pts
                .iter()
                .map(|&b| crate::rotated_vertex::replay(WitnessPoint::at(b), &chain))
                .collect();
            let Some(w) = turned else { continue }; // conservative: keep the world description
            let mut w: [WitnessPoint; 3] = [w[0].clone(), w[1].clone(), w[2].clone()];
            let e1 = Vector3::from_array(w[1].coord()) - Vector3::from_array(w[0].coord());
            let e2 = Vector3::from_array(w[2].coord()) - Vector3::from_array(w[0].coord());
            if e1.cross(e2).dot(f.n_out) < 0.0 {
                w.swap(1, 2);
            }
            f.tri_pt3 = w;
            f.rotated = true;
            f.motion = Some(c);
        }
    }
    Ok(out)
}
