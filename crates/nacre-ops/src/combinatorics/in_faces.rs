use super::*;
/// **Is the point `p` inside this component**, asked along the ray `p + t·dir` — the coordinate
/// twin of [`point_in_component`].
///
/// ★★★ **Two roads, one crossing rule.** The two differ only in how the *point* arrives: a
/// three-plane name (exact without coordinates, so it survives a rotated class) or rational
/// coordinates (the only thing a curved component can offer, since none of its faces carries a
/// vertex). Everything they ask about a **face** — a circle's radial side, a band's roots and
/// axial span, the abandon-on-boundary policy — is the same rule and is called, not restated.
/// Restating it is how the two would come to disagree about a graze.
///
/// ★★ **What is measured, and what is owed.** `an_enclosed_cylindrical_void_is_a_cavity` is the
/// fixture that makes this road say `true`, and it goes red if the road is stubbed to `false` —
/// two *disjoint* bodies would not, since their answer is "outside" whatever the road does.
///
/// ★ The **curved** arm is barely loaded from here: this road is entered a handful of times in
/// the suite and mostly looks at no cylinder face at all (the other component is all planes —
/// the void fixture's is a box); where it does (`two_cylinders_with_coplanar_caps_fuse_apart`)
/// the ray misses. Nor can a fixture with parallel axes do better: the gate refuses any boolean
/// whose two cylinders' faces it cannot prove apart (spans and arcs), so a ray from one
/// cap's centre crosses the other lateral **0 or 2 times** and the parity is the same. What loads the lateral
/// road ([`lateral_face_crossings`]) is the **named** road's probes — the crossing
/// census's corner Commons, whose vertex rays cross the other half's panel — and the lattice
/// oracle on the through-boss; a real `k = 1` from a coordinate probe wants a ∥ wall inside the
/// strip — a population the gate **serves** rather than refuses (a crossing or a
/// tangency), so what owes this arm a fixture is the rulings road, not a refusal.
///
/// `Ok(None)` = this ray grazed; the caller has other directions to try.
pub(crate) fn point_in_faces_rat(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: &[nacre_exact::Rat; 3],
    dir: &[nacre_exact::Rat; 3],
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    use nacre_exact::Rat;
    use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
    let not_rational = || reject(RejectReason::WitnessNotRational);
    let zero = Rat::from_int(0);
    if dir.iter().all(|c| *c == zero) {
        return Err(reject(RejectReason::DegenerateWitness));
    }
    // The half this road counts is **ahead** of the origin, so its cut plane's negative side is
    // `dir·(x − p) > 0` — the mirror of the named road, which counts behind. Either half has the
    // same parity; what matters is that one road picks one.
    let ahead = (|| {
        let n = [
            zero.checked_sub(dir[0])?,
            zero.checked_sub(dir[1])?,
            zero.checked_sub(dir[2])?,
        ];
        Some([
            n[0],
            n[1],
            n[2],
            zero.checked_sub(nacre_exact::dot3_rat(&n, p)?)?,
        ])
    })()
    .ok_or_else(not_rational)?;
    let line = planes_through_line(p, dir).ok_or_else(not_rational)?;

    let mut count = 0usize;
    for f in faces {
        let q = match &f.surf {
            CompSurf::Cylinder(def) => {
                #[cfg(test)]
                cylinder_asks::asked(f);
                // A circle or polygon outer on a cylinder face has no producer; refusing to
                // guess costs the caller another direction, never a wrong answer.
                let BoundEdges::Lateral(loops) = &f.outer else {
                    return Ok(None);
                };
                match lateral_face_crossings(jd, [&line[0], &line[1]], def, loops, &ahead) {
                    Some(CurvedHit::Counted(k)) => {
                        count += k;
                        continue;
                    }
                    Some(CurvedHit::Graze) | None => return Ok(None),
                }
            }
            CompSurf::Plane(q) => *q,
        };
        let coeffs = class_coeffs_rat(jd, q).ok_or_else(not_rational)?;
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let nd = nacre_exact::dot3_rat(&n, dir).ok_or_else(not_rational)?;
        let residual = nacre_exact::dot3_rat(&n, p)
            .and_then(|v| v.checked_add(coeffs[3]))
            .ok_or_else(not_rational)?;
        // ★ **Parallel is the origin's question, not the ray's.** `nd == 0` and a nonzero
        // residual is a plain miss. `nd == 0` with a zero residual is the ray lying **in** this
        // plane, where it never passes from one side of the face to the other — zero crossings —
        // *unless* the origin is on the face itself, and then "is the origin inside the component"
        // has no answer. That is exactly the `t == 0` policy below, so the two are one arm: a
        // coplanar ray is the origin's own crossing.
        let t = if nd == zero {
            if residual != zero {
                continue;
            }
            zero
        } else {
            (|| -> Option<Rat> {
                zero.checked_sub(residual)?
                    .checked_mul(Rat::new(nd.denom(), nd.numer())?)
            })()
            .ok_or_else(not_rational)?
        };
        if t < zero {
            continue; // behind the origin
        }
        let x = (|| -> Option<[Rat; 3]> {
            let mut x = *p;
            for k in 0..3 {
                x[k] = x[k].checked_add(t.checked_mul(dir[k])?)?;
            }
            Some(x)
        })()
        .ok_or_else(not_rational)?;
        let chart = Chart2dRat::of_normal(&n).ok_or_else(not_rational)?;
        let x2 = chart.project(&x).ok_or_else(not_rational)?;
        // Each bound answers by its own kind — the same split the named road makes.
        let inside = |b: &BoundEdges| -> Result<Option<bool>, BoolError> {
            match b {
                BoundEdges::Ring(r) => {
                    // A mixed ring forks to the rational walk here exactly as the named
                    // road forks: the crossing x is already rational, and the chart-ring
                    // derivation below has no spelling for a pierce corner or an arc step.
                    if ring_is_mixed(r) {
                        return Ok(point_in_mixed_ring(jd, cyls, &coeffs, &x, r));
                    }
                    // ★ Under three nodes `point_in_ring_2d_rat` answers `Outside` by contract,
                    // which would make a degenerate ring *invisible* to the parity instead of
                    // loud. A **polygon** face of a valid solid has no such ring, so saying so is
                    // free.
                    //
                    // ★★★★★ **It used to stand before the fork, and that read a curved face as
                    // degenerate.** A half-disc cap's ring is two edges — an arc and its chord —
                    // which is a perfectly good boundary and not a polygon at all; the sentence
                    // "a face of a valid solid has no such ring" was only ever true of the road
                    // *below*. Measured: the moment a wall boss's caps were given a witness, six
                    // census cells came here and were refused by name for being what they are.
                    if r.len() < 3 {
                        return Err(reject(RejectReason::DegenerateRing));
                    }
                    // The pierce check is spelled before the chart rather than left to
                    // `node_coords_rat`'s `None`, so an unnamed vertex is reported as itself and
                    // not as arithmetic that ran out of room — the same split `every_ray` makes.
                    let nodes: Vec<NodeId> = r
                        .iter()
                        .map(|e| three_plane_name(e.node).map(NodeId::ThreePlane))
                        .collect::<Option<_>>()
                        .ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?;
                    let ring2 = chart.ring(jd, &nodes).ok_or_else(not_rational)?;
                    Ok(match point_in_ring_2d_rat(x2, &ring2) {
                        RingSide::Inside => Some(true),
                        RingSide::Outside => Some(false),
                        RingSide::OnBoundary => None,
                    })
                }
                BoundEdges::Circle(def) => Ok(point_in_disk(&x, def)),
                BoundEdges::Lateral(_) => Ok(None),
            }
        };
        let Some(material) = material_of(f, inside)? else {
            return Ok(None);
        };
        if t == zero {
            // The crossing is the ray's own origin — either the ray meets this plane there, or it
            // lies in it. On the face's material the query sits on the component's boundary, where
            // "inside" has no answer; off it, the plane is touched (or run along) at points
            // outside the face and counts for nothing.
            if material {
                return Ok(None);
            }
            continue;
        }
        if material {
            count += 1;
        }
    }
    Ok(Some(count % 2 == 1))
}

/// **Is `p` strictly inside this circle?** — `None` when it lies *on* the circle (non-generic,
/// abandon) or when the arithmetic could not answer.
///
/// ★ It takes **coordinates**, so both roads ask it: the named probe realizes its three-plane
/// crossing first, the coordinate road already has one. The rule must not be written twice.
///
/// ★ The rule is `cylinder_radial_side`'s, the one the nesting engine reads for a disk target —
/// a circle bound is `cylinder ∩ plane`, so "inside the disk" is "inside the cylinder's radius".
fn point_in_disk(p: &[nacre_exact::Rat; 3], def: &nacre_topo::CylinderDef) -> Option<bool> {
    match nacre_exact::cylinder_radial_side(p, &def.origin(), &def.dir(), def.r2()) {
        nacre_exact::Orient::Negative => Some(true),
        nacre_exact::Orient::Positive => Some(false),
        nacre_exact::Orient::Zero => None,
    }
}

/// [`lateral_face_crossings`] for a probe ray named by three plane classes: it states the ray's
/// two planes and the "behind" cut plane from the class table and hands them over.
///
/// `Ok(None)` — as `Some(None)` here — is the abandon-this-ray answer; `Err` never happens because
/// every failure this can meet is arithmetic, and arithmetic that cannot answer is an abstention.
fn curved_count(
    jd: &Judge<'_, WorkingPlane>,
    a: usize,
    b: usize,
    c: usize,
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
) -> Result<Option<usize>, BoolError> {
    // ★ The cylinder gate refuses any class that is rotated or has no narrow rational name
    // (`planes`), so in a boolean that has a cylinder these are always `Some`. A `None` here
    // would mean that gate let something through — assert it, then abstain rather than guess.
    let Some(((ca, cb), cc)) = class_coeffs_rat(jd, a)
        .zip(class_coeffs_rat(jd, b))
        .zip(class_coeffs_rat(jd, c))
    else {
        debug_assert!(
            false,
            "the cylinder gate is supposed to make these rational"
        );
        return Ok(None);
    };
    // The cut plane through the ray's origin, normal `n_a × n_b` — the very direction
    // `plane_plane_cylinder` gives its meet line, so "negative side" is "behind the query".
    let Some(half) = (|| {
        let d = nacre_exact::cross3_rat(&[ca[0], ca[1], ca[2]], &[cb[0], cb[1], cb[2]])?;
        let at = nacre_exact::three_planes_rat([ca, cb, cc])?;
        let d0 = nacre_exact::Rat::from_int(0).checked_sub(nacre_exact::dot3_rat(&d, &at)?)?;
        Some([d[0], d[1], d[2], d0])
    })() else {
        return Ok(None);
    };
    Ok(
        match lateral_face_crossings(jd, [&ca, &cb], def, loops, &half) {
            Some(CurvedHit::Counted(k)) => Some(k),
            Some(CurvedHit::Graze) | None => None,
        },
    )
}

/// A rim's plane restated with the **axis** as its normal — `[m, −m·(o + t·m)]`, so
/// [`nacre_exact::quad::plane_side`] at a point of the cylinder is the sign of `z − t` along
/// the axis. The banded arm built it this way before the cutover; the loops road builds it for
/// every ⊥ class it meets.
fn rim_plane(def: &nacre_topo::CylinderDef, t: nacre_exact::Rat) -> Option<[nacre_exact::Rat; 4]> {
    let (o, m) = (def.origin(), def.dir());
    let mut d0 = nacre_exact::Rat::from_int(0);
    for k in 0..3 {
        let at = o[k].checked_add(t.checked_mul(m[k])?)?;
        d0 = d0.checked_sub(m[k].checked_mul(at)?)?;
    }
    Some([m[0], m[1], m[2], d0])
}

/// The axis parameter of the ⊥ class among a pierce corner's two naming planes — the `z` of
/// the arc that ends there, or of a ruling piece's end. Exactly one of the two is ⊥ in this
/// population (two ⊥ planes never meet, and [`crate::planes::axis_param_of_plane`] answers only
/// for `n · m ≠ 0`); `None` names a corner without one — a tilted cut or a three-plane
/// name where a pierce was expected.
fn corner_axis_param(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<nacre_exact::Rat> {
    let (planes, _, _) = pierce_name(n)?;
    planes
        .iter()
        .find_map(|&c| crate::planes::axis_param_of_plane(&class_coeffs_rat(jd, c)?, def))
}

/// The ∥ class among a pierce corner's naming planes — the wall a ruling piece ending there
/// rides: the one whose normal is ⊥ to the axis (`n · m = 0`), through the axis or offset from
/// it; the other name is the ⊥ class the arc rides.
fn corner_wall_class(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<usize> {
    let (planes, _, _) = pierce_name(n)?;
    let m = def.dir();
    planes.iter().copied().find(|&c| {
        class_coeffs_rat(jd, c).is_some_and(|w| {
            nacre_exact::dot3_rat(&[w[0], w[1], w[2]], &m) == Some(nacre_exact::Rat::from_int(0))
        })
    })
}

/// **Is the cylinder point `x` on the lateral face these loops bound?** — the parity of the
/// ray up the axis from `x` against the boundary loops, on the chart `(θ, z)`.
///
/// The rule is the planar rings' half-open rule read on this chart. A whole-circle rim is
/// crossed iff it is above `x`. A ring's **arc** (z = const, a CCW span in θ) is crossed iff
/// it is above `x` and `θ_x` is in its span — with a span end exactly at `θ_x` counted by the
/// end the arc *leaves upward* (`+θ`): `lo` counts, `hi` does not
/// ([`nacre_geom::intersect::ray_step_crossing`]`(Zero, ±, ±)`), and the ruling at that corner
/// runs along the ray and never counts. A **ruling** piece is crossed by nothing; `x` on it is
/// the boundary. Every comparison is one the arrangement already owns: `z` by
/// [`nacre_exact::quad::plane_side`] against [`rim_plane`], `θ` by [`arc_span`], a ruling by
/// its own name (`plane_side(wall) == 0 ∧ ruling_side == side`, the predicate
/// `crossing_on_ruling` names stations with).
///
/// `None` = `x` on the boundary (a rim, an arc, a ruling — `tie_probe` says which), a seam tie
/// inside `arc_span`, a loop the road cannot read (a tilted arc, a plane carrier), or checked
/// arithmetic running out. `z` is asked before `θ`: an arc below `x` needs no order, and that is
/// what lets a ray whose seam-incident root has no arc above it answer.
pub(crate) fn loop_parity(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
    meet: &nacre_exact::quad::MeetLine,
    s: &nacre_exact::quad::QuadVal,
) -> Option<bool> {
    use nacre_exact::Orient;
    use nacre_exact::quad::plane_side;
    let above =
        |t: nacre_exact::Rat| -> Option<Orient> { Some(plane_side(&rim_plane(def, t)?, meet, s)) };
    let mut crossings = 0usize;
    for lp in loops {
        match lp {
            LateralLoop::Circle(c) => {
                let t = crate::planes::axis_param_of_plane(&class_coeffs_rat(jd, *c)?, def)?;
                match above(t)? {
                    // `z_x < t`: the rim is above the point, and the ray up the axis crosses it.
                    Orient::Negative => crossings += 1,
                    Orient::Positive => {}
                    Orient::Zero => {
                        #[cfg(test)]
                        tie_probe::push(tie_probe::Tie::OnRim);
                        return None;
                    }
                }
            }
            LateralLoop::Ring(edges) => {
                for e in edges {
                    match &e.carrier {
                        Carrier::Arc(arc) => {
                            let Some(t) = corner_axis_param(jd, def, e.node) else {
                                #[cfg(test)]
                                tie_probe::push(tie_probe::Tie::TiltedArc);
                                return None;
                            };
                            let side = above(t)?;
                            if side == Orient::Positive {
                                continue; // the arc is below the point
                            }
                            let (lo_nd, hi_nd) = if arc.ccw {
                                (e.node, e.to)
                            } else {
                                (e.to, e.node)
                            };
                            let e_lo = pierce_meet(jd, arc.cyl, &arc.def, lo_nd)?;
                            let e_hi = pierce_meet(jd, arc.cyl, &arc.def, hi_nd)?;
                            let Some(span) = arc_span(&arc.def, &e_lo, &e_hi, &(meet.clone(), *s))
                            else {
                                // A seam tie (the root on the seam, two seam ends, a zero span)
                                // or arithmetic out — the mark says which.
                                #[cfg(test)]
                                tie_probe::flush_or(tie_probe::Tie::Other);
                                return None;
                            };
                            if side == Orient::Zero {
                                // At the arc's own z: on the arc iff within its closed span.
                                if span == ArcSpan::Outside {
                                    continue;
                                }
                                #[cfg(test)]
                                tie_probe::push(tie_probe::Tie::OnArc);
                                return None;
                            }
                            if matches!(span, ArcSpan::Inside | ArcSpan::AtLo) {
                                crossings += 1;
                            }
                        }
                        Carrier::Ruling(rl) => {
                            let (Some(t0), Some(t1)) = (
                                corner_axis_param(jd, def, e.node),
                                corner_axis_param(jd, def, e.to),
                            ) else {
                                #[cfg(test)]
                                tie_probe::push(tie_probe::Tie::Producer);
                                return None;
                            };
                            let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
                            // Outside the piece's closed axial range: not on it, whatever θ.
                            if above(lo)? == Orient::Negative || above(hi)? == Orient::Positive {
                                continue;
                            }
                            let Some(fc) = corner_wall_class(jd, def, e.node) else {
                                #[cfg(test)]
                                tie_probe::push(tie_probe::Tie::Producer);
                                return None;
                            };
                            let w = class_coeffs_rat(jd, fc)?;
                            if plane_side(&w, meet, s) == Orient::Zero
                                && ruling_side(&w, def, (meet, s)) == Some(rl.side)
                            {
                                #[cfg(test)]
                                tie_probe::push(tie_probe::Tie::OnRuling);
                                return None;
                            }
                        }
                        Carrier::Plane { .. } => {
                            #[cfg(test)]
                            tie_probe::push(tie_probe::Tie::Producer);
                            return None;
                        }
                    }
                }
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// **Where the line `line[0] ∩ line[1]` crosses one lateral face, and how many of those lie on
/// the half the caller is counting.**
///
/// ★★★ **One copy, two roads.** The named-point probe and the coordinate-point road both state
/// their ray as *two rational planes*, so both ask this. Writing the arms twice is how the two
/// would come to disagree about a graze. A `half` is a rational plane through the ray's origin,
/// and a crossing counts exactly when it is on that plane's **negative** side — so each road
/// states its own half and nothing here has a front or a back: the named road's normal is the
/// line's own direction `n_a × n_b` (which counts **behind** the query — `plane_plane_cylinder`
/// builds the meet line's `dir` as `cross3(n_a, n_b)`, the very `d` that `order_along`'s
/// `fwd == 1` measures against), the coordinate road's is `−dir` (**ahead** of the origin).
///
/// Each root asks the face by [`loop_parity`] **before** the half — a root off the face is not
/// a crossing whichever side it is on, and a root on the face at the ray's own origin is the
/// query on the other component's surface (`Graze`). ★ This used to be a banded arm reading a
/// band's two rims as an axial span (the crossing between them iff on opposite sides of the two
/// rim planes) and abstaining on every other lateral by name (`MissOnly`); the two whole-circle
/// loops say the same thing — measured identical on 2,180 lattice rays, every root and every
/// graze — and a panel, a chain rim or a hole is now read rather than passed over.
///
/// `None` is checked-`Rat` arithmetic that could not answer — an honest decline, never a guess.
pub(crate) fn lateral_face_crossings(
    jd: &Judge<'_, WorkingPlane>,
    line: [&[nacre_exact::Rat; 4]; 2],
    def: &nacre_topo::CylinderDef,
    loops: &[LateralLoop],
    half: &[nacre_exact::Rat; 4],
) -> Option<CurvedHit> {
    use nacre_exact::Orient;
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (meet, roots): (_, [QuadVal; 2]) =
        match nacre_exact::quad::plane_plane_cylinder(line[0], line[1], &o, &m, r2)? {
            CylinderMeet::Pair { line, s } => (line, s),
            CylinderMeet::Tangent { .. } | CylinderMeet::OnRuling(_) => {
                return Some(CurvedHit::Graze);
            }
            CylinderMeet::AxisParallelMiss(_) | CylinderMeet::Miss(_) => {
                return Some(CurvedHit::Counted(0));
            }
            CylinderMeet::CoincidentPlanes | CylinderMeet::ParallelPlanes => {
                return Some(CurvedHit::Graze);
            }
        };
    let mut count = 0usize;
    for s in &roots {
        match loop_parity(jd, def, loops, &meet, s) {
            // On the boundary, a seam tie, an unreadable loop, arithmetic out: this ray
            // cannot count this face — the banded arm's `Graze` on a rim, generalized.
            None => return Some(CurvedHit::Graze),
            Some(false) => continue,
            Some(true) => match nacre_exact::quad::plane_side(half, &meet, s) {
                Orient::Zero => return Some(CurvedHit::Graze),
                Orient::Negative => count += 1,
                Orient::Positive => {}
            },
        }
    }
    Some(CurvedHit::Counted(count))
}

/// A component as its faces, each boundary already carrying its edges' walls.
pub(crate) type ComponentFaces = Vec<CompFace>;

pub(crate) fn point_in_component(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    query: [usize; 3],
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    let mut vplanes = query.to_vec();
    vplanes.sort_unstable();
    vplanes.dedup();

    // One ray attempt along `L = a ∩ b`, located by the query point `{a,b,c}`. `Ok(None)` = the
    // line grazed and the caller should try the next plane pair.
    let attempt = |a: usize, b: usize, c: usize| -> Result<Option<bool>, BoolError> {
        let mut count = 0usize;
        for f in faces {
            // ★ A cylindrical face is counted by its own arm — the crossings are roots of a
            // quadratic, not three-plane points, and "inside the face" is an axial span when
            // the face states one, a bare miss-oracle when it does not.
            if let CompSurf::Cylinder(def) = &f.surf {
                #[cfg(test)]
                cylinder_asks::asked(f);
                // A circle or polygon outer on a cylinder face has no producer; refusing
                // to guess costs the caller another node, never a wrong answer.
                let BoundEdges::Lateral(loops) = &f.outer else {
                    return Ok(None);
                };
                match curved_count(jd, a, b, c, def, loops)? {
                    Some(k) => {
                        count += k;
                        continue;
                    }
                    None => return Ok(None),
                }
            }
            let CompSurf::Plane(q) = &f.surf else {
                return Ok(None);
            };
            let q = *q;
            // Inside `q`'s material at the point `x`: inside the outer bound, outside every hole.
            // ★ Each bound answers by its own kind — a polygon by a ring walk, a circle by its
            // radial side — and either can say "the point is *on* me", which is the same
            // non-generic abandon `point_on_ring` has always raised.
            let material = |x: [usize; 3]| -> Result<Option<bool>, BoolError> {
                let inside = |b: &BoundEdges| -> Result<Option<bool>, BoolError> {
                    match b {
                        BoundEdges::Ring(r) => {
                            // A mixed ring (an arc step, a pierce corner) takes the
                            // rational road: the crossing X is a rational three-plane
                            // point, the ring is walked by its carriers, and every tie
                            // abstains for the next probe. The named walk cannot read a
                            // pierce corner at all - letting it try would answer
                            // `RingNaming`, a false name for the cause.
                            if ring_is_mixed(r) {
                                let Some(coeffs) = class_coeffs_rat(jd, q) else {
                                    return Ok(None); // no exact class statement: abstain
                                };
                                let Some(px) =
                                    node_coords_rat(jd, NodeId::three_planes(Canon3::three(x)))
                                else {
                                    return Ok(None);
                                };
                                return Ok(point_in_mixed_ring(jd, cyls, &coeffs, &px, r));
                            }
                            if point_on_ring(jd, q, x, r)? {
                                return Ok(None);
                            }
                            Ok(every_ray(jd, q, x, r)?.first().copied())
                        }
                        BoundEdges::Circle(def) => {
                            Ok(node_coords_rat(jd, NodeId::three_planes(Canon3::three(x)))
                                .and_then(|p| point_in_disk(&p, def)))
                        }
                        // Loops bound a cylinder, never a plane — a producer error, not an input.
                        BoundEdges::Lateral(_) => Ok(None),
                    }
                };
                material_of(f, inside)
            };
            if jd.plane_pair_dir_sign(a, b, q) == 0 {
                // `L` is parallel to `q` — and possibly **in** it, which is not the same thing.
                // A coplanar ray never passes from one side of this face to the other, so zero
                // crossings is the right count and always was. What the old `continue` also
                // swallowed is the *query*: if it lies on this face, "is the query inside the
                // component" has no answer at all, and skipping the face answers it anyway. That
                // is the same proposition the `fwd == 0` arm below abandons for a transversal
                // plane — one rule that had only one of its two spellings.
                //
                // ★ **Unfired, and measured to be.** 26 rays in the suite lie in a face's plane
                // and **none** of them is on that face's material. (An earlier count said two;
                // it read the *outer bound* rather than the face, and both were in a hole.) It is
                // here because the sentence is true, not because a fixture is red.
                let mut vq = [query[0], query[1], query[2]];
                vq.sort_unstable();
                if side_of(jd, &[], NodeId::three_planes(Canon3::three(vq)), q) == Some(0)
                    && material(vq)? != Some(false)
                {
                    return Ok(None);
                }
                continue;
            }
            let mut x = [a, b, q];
            x.sort_unstable();
            let Some(in_g) = material(x)? else {
                return Ok(None);
            };
            let fwd = order_along(jd, a, b, c, q);
            if fwd == 0 {
                // The crossing is the ray origin itself (query on plane `q`). Inside `q`'s
                // material ⇒ query on the component's surface ⇒ undecidable → abandon.
                if in_g {
                    return Ok(None);
                }
                continue;
            }
            if fwd == 1 && in_g {
                count += 1; // one side of the line — parity is the same on either half
            }
        }
        Ok(Some(count % 2 == 1))
    };

    for i in 0..vplanes.len() {
        for j in (i + 1)..vplanes.len() {
            let (a, b) = (vplanes[i], vplanes[j]);
            // Locator `c`: a plane of the query off the line `a ∩ b`, so `{a,b,c}` is the query.
            let Some(&c) = vplanes
                .iter()
                .find(|&&x| x != a && x != b && jd.plane_pair_dir_sign(a, b, x) != 0)
            else {
                continue;
            };
            if let Some(inside) = attempt(a, b, c)? {
                return Ok(Some(inside));
            }
        }
    }
    // Every plane pair of this query was blocked: the node abstains. The inner `attempt`
    // already speaks this language per pair (`Ok(None)`); the boundary now keeps it instead of
    // dressing the abstention up as an error for the caller to catch and swallow.
    Ok(None)
}
