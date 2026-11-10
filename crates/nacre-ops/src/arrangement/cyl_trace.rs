use super::*;
/// What a lateral face leaves on a plane class that crosses its axis — see [`circle_on_class`], which
/// answers it **per angular extent** ([`CircleSpan`]) rather than once for the whole circle. On a ⊥
/// class the trace is a circle; on a class sectioning the cylinder obliquely
/// ([`crate::planes::CylRecord::Sections`]) it is the ellipse, read the same way — its extents are
/// angles about the axis, which order the ellipse's points as they order the circle's
/// ([`trace_on_face`]).
///
/// The empty answer is "this class leaves nothing here": a class along the axis carries no circle
/// at all, a ⊥ class beyond both rims never meets the face, an oblique section misses it, and a
/// class the gate proved silent about the cylinder stands clear of it. `Err` is "the face cannot answer" (no rim span) — declining beats a
/// silently missing circle, which would corrupt every label on the class.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CylOnClass {
    /// The class cuts the face in two: the solid straddles it here.
    Crosses,
    /// The class carries one of the face's **rims**: the wall touches it along that circle with
    /// its body to one side. `body_above` is that side, about the class's **stored** normal —
    /// the frame every cell label is written in.
    Grazes { body_above: bool },
}

pub(super) fn circle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    cyl: usize,
    found: &mut Aliases,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    // Two names the circle puts at one angle are one point — one circle, one angle; an exact proof
    // the round's table records ([`Aliases::record_same_point`]) whatever this face's answer is.
    let mut ties: Vec<(NodeId, NodeId)> = Vec::new();
    let spans = circle_spans(jd, cyls, cf, fl, wc, cyl, &mut ties);
    for (a, b) in ties {
        found.record_same_point_traced(jd, cyls, a, b)?;
    }
    spans
}

/// [`circle_on_class`] before the ties its traces found are recorded.
fn circle_spans(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    wc: usize,
    cyl: usize,
    ties: &mut Vec<(NodeId, NodeId)>,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    // A class the gate proved clear of this cylinder carries none of it — the gate's record, not a
    // question this road could ask of a class with no world statement.
    if cyls.get(cyl).is_some_and(|cy| cy.is_silent(wc)) {
        return Ok(Vec::new());
    }
    // The world description — the cylinder's statement is world, and a comparison across two
    // frames is a silently wrong answer, not a slow one.
    let Some(world) = jd.planes[wc].world.as_ref() else {
        return Err(DeclineKind::CylSpan);
    };
    let Some(coeffs) = world.name.narrow().copied() else {
        return Err(DeclineKind::CylSpan);
    };
    // No world statement for this lateral (a rotated or frame-borne truth): the same decline
    // the class side makes — both sides of this comparison have to speak about the world.
    let Some(def) = cf.def.as_ref() else {
        return Err(DeclineKind::CylSpan);
    };
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // A class along the axis carries no section — rulings, a miss here like a parallel plane, not a
    // decline. An oblique one carries the cylinder's ellipse where the gate recorded it as the
    // class's section ([`crate::planes::CylRecord::Sections`]); anywhere else the gate proved it
    // clear. Decided **totally** (`axis_relation` clears denominators and answers in integers).
    let oblique = match nacre_exact::axis_relation(&n, &def.dir()) {
        nacre_exact::AxisRelation::Along => return Ok(Vec::new()),
        nacre_exact::AxisRelation::Oblique => {
            if !cyls.get(cyl).is_some_and(|cy| cy.sections_oblique(wc)) {
                return Ok(Vec::new());
            }
            true
        }
        nacre_exact::AxisRelation::Across => false,
    };
    // **The face's boundary as the tracer names it — its cycles — and nothing else.** A band's
    // two rims and its holes are its loops (`combinatorics::lateral_cycles`); a panel or a chain
    // rim is a cycle like a hole, carved out of an outer answer that may be **absent** rather than
    // out of a whole trace.
    let shape = lateral_shape(jd, cf, fl, def)?;
    let mut rings: Vec<(bool, &combinatorics::NamedRing)> = Vec::with_capacity(shape.cycles.len());
    for (inner, ring) in &shape.cycles {
        // A one-edge loop whose far face is a cylinder is two laterals meeting: M6b's pair, not
        // this road's. It has no ring to walk, so it declines rather than passing unread.
        rings.push((*inner, ring.poly().ok_or(DeclineKind::CylFaceHole)?));
    }
    let read = trace_on_face(
        jd, cyls, cf, def, wc, cyl, world, &coeffs, &shape, &rings, ties,
    );
    // The face's shape is read; on an oblique section whatever the trace then cannot be drawn
    // from is the oblique coverage's end, one name for every way it ends.
    if oblique {
        read.map_err(|_| DeclineKind::ObliqueSection)
    } else {
        read
    }
}

/// **What a lateral face leaves on a plane class not along its axis** — the extents of the class's
/// section of the cylinder (a circle across the axis, an ellipse oblique to it) that the face
/// covers, past [`circle_spans`]' reading of the class and the face's shape.
///
/// **One carving for every section.** The trace goes once round the cylinder — every generator
/// meets the plane once — and the face's boundary is sections too (its rims, its cycles' arcs) and
/// rulings. Each rim is asked the one question of two sections
/// ([`nacre_exact::quad::section_meet`]): the class's own plane is the rim, and the face grazes it
/// along the whole trace, its body on the rim's face side; apart takes nothing; crossing cuts the
/// trace at the two points the line `wc ∩ rim` meets the cylinder in
/// ([`combinatorics::rim_crossings`], the names a cap's chord ends at too), the face present where
/// the trace, going counter-clockwise from a point, climbs into it ([`climb_over`] — the rim's
/// face side, read off which end of the face it is; two circles never cross); a touch declines (the
/// gate refuses it first). Every other cycle — a fillet's quarter, a hole, a panel — is walked by
/// [`cycle_on_class`], which names its crossings and lays the runs of its arcs on the trace.
///
/// **Out of reach first**: a trace beside every boundary section of the face, on one side, meets no
/// part of it — on each generator the face lies between its own boundary points — so nothing, and
/// no cycle is walked.
///
/// **The outer answer**, which the cuts carve: a graze of a rim is that graze; otherwise it is what
/// one point of the uncarved trace answers ([`outer_by_witness`]) — the uncarved part is all in the
/// face or all out of it.
#[allow(clippy::too_many_arguments)]
fn trace_on_face(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    cyl: usize,
    world: &crate::planes::WorldName,
    coeffs: &[Rat; 4],
    shape: &LateralShape,
    rings: &[(bool, &combinatorics::NamedRing)],
    ties: &mut Vec<(NodeId, NodeId)>,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    use core::cmp::Ordering::{Equal, Greater, Less};
    use nacre_exact::quad::{
        SectionMeet, SectionOrder, Toward, section_height_order, section_meet,
    };
    // Which side of the class the face's body is on is a label in the class's **stored** normal's
    // frame: the face runs from its low end toward `+t` and from its high end toward `−t`.
    let up = crate::planes::plus_t_is_above(world, def);
    let world_of = |c: usize| jd.world_rat(c).ok_or(DeclineKind::CylSpan);
    // Where the trace stands to a boundary section: on it (`Equal`), beside it (`Less`/`Greater`,
    // with no generator — two sections apart cross nowhere), or across it (`None`).
    let stand = |c: usize| -> Result<Option<core::cmp::Ordering>, DeclineKind> {
        let w = world_of(c)?;
        match section_meet(coeffs, &w, def) {
            Some(SectionMeet::Same) => Ok(Some(Equal)),
            Some(SectionMeet::Cross) => Ok(None),
            Some(SectionMeet::Apart) => {
                match section_height_order(coeffs, &w, def, None, Toward::At) {
                    Some(SectionOrder::Ordered(o)) => Ok(Some(o)),
                    _ => Err(DeclineKind::CylSpan),
                }
            }
            Some(SectionMeet::Touch) | None => Err(DeclineKind::ObliqueSection),
        }
    };
    let beside = |side| {
        shape
            .sections
            .iter()
            .all(|&c| stand(c).ok() == Some(Some(side)))
    };
    if beside(Less) || beside(Greater) {
        return Ok(Vec::new());
    }
    let mut carved: Vec<Carved> = Vec::new();
    let mut graze: Option<CylOnClass> = None;
    for &(c, above) in &shape.rims {
        match stand(c)? {
            Some(Equal) => {
                graze = Some(CylOnClass::Grazes {
                    body_above: if above { up } else { !up },
                });
            }
            Some(_) => {}
            None => {
                let Some(Some(roots)) = combinatorics::rim_crossings(jd, def, wc, c, cyl) else {
                    return Err(DeclineKind::ObliqueSection);
                };
                let toward = if above { 1 } else { -1 };
                let rim = world_of(c)?;
                let mut cuts: Vec<(NodeId, bool)> = Vec::with_capacity(2);
                for node in roots {
                    let climb = climb_over(jd, cyl, def, coeffs, &rim, node)?;
                    cuts.push((node, climb != toward));
                }
                carve_between(jd, cyl, def, &cuts, &mut carved)?;
            }
        }
    }
    let mut contacts: Vec<NodeId> = Vec::new();
    // **The hole is read by the walk every other face's boundary is read by** — the face's
    // own loop, against this class, through [`combinatorics::ring_against_plane`]. An interval
    // derived from the loop's carriers would be exact for a chart rectangle and a *premise* for
    // anything else; the walk asks the ring and needs no premise about the hole's shape.
    for &(inner, nr) in rings {
        cycle_on_class(
            jd,
            cyls,
            cf,
            def,
            wc,
            cyl,
            up,
            inner,
            nr,
            coeffs,
            &mut carved,
            &mut contacts,
        )?;
    }
    let outer = match graze {
        Some(g) => Some(g),
        None => {
            let loops = lateral_loops(jd, cyls, shape, rings)?;
            outer_by_witness(jd, def, cyl, coeffs, &loops, &carved, &contacts)?
        }
    };
    if carved.is_empty() {
        return Ok(outer
            .map(|on| vec![CircleSpan { arc: None, on }])
            .unwrap_or_default());
    }
    assemble_spans(jd, cyl, &cyls[cyl].def, &carved, outer, ties)
}

/// **The face's boundary as loops on its chart** — what [`combinatorics::loop_parity`] reads: each
/// rim (closed or split into arcs) as its plane class, every other cycle — the polygon rings
/// [`trace_on_face`] was handed, a lateral's cycles that are no rim — as its edges
/// ([`combinatorics::ring_edges`], with no plane — a lateral's corners are pierce names).
fn lateral_loops(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    shape: &LateralShape,
    rings: &[(bool, &combinatorics::NamedRing)],
) -> Result<Vec<combinatorics::LateralLoop>, DeclineKind> {
    use combinatorics::LateralLoop;
    let mut loops: Vec<LateralLoop> = shape
        .rims
        .iter()
        .map(|&(c, _)| LateralLoop::Rim(c))
        .collect();
    for (_, nr) in rings {
        loops.push(LateralLoop::Ring(
            combinatorics::ring_edges(jd, cyls, &nr.triples, &nr.walls, None)
                .map_err(|_| DeclineKind::CylSpan)?,
        ));
    }
    Ok(loops)
}

/// **What the trace's uncarved part answers** — one point of it, asked of the face's loops.
///
/// Every simple closed curve on a cylinder cuts it in two, so a lateral face is the meet of its
/// boundary loops' material sides. Carving took away, at every loop the trace crosses, the trace's
/// part on that loop's other side, and laid the runs along a loop as grazes; a touch takes nothing
/// and changes no side, and a run whose flanks differ is declined before this. So every point the
/// carving left is on the material side of each loop the trace crossed, and on one side — the same
/// along the whole trace — of each loop it did not: **the uncarved part is all in the face or all
/// out of it**, and one point of it says which, by the rule every lateral point is read by
/// ([`combinatorics::loop_parity`]).
///
/// The point is chosen where it cannot be on the boundary. The trace meets the face's boundary only
/// at `contacts` and the carved extents' ends (every node a cycle has on the class, every crossing)
/// — a touch inside a boundary arc is neither, and the gate refuses a plane that touches a boundary
/// section (`boundary_stand`); a witness on one would read `OnBoundary` and decline.
/// With none of them any point of the trace is off the boundary, and its generator is realized only
/// if a question needs one ([`combinatorics::ChartPoint`] — a face of whole rims apart from the trace
/// needs none). Otherwise the point stands strictly between two consecutive ones in θ, on an
/// elementary arc no extent claims: a root of the trace's plane and a plane through the axis whose
/// direction is the realized middle of that arc — a sample, not a truth — **verified** inside the
/// arc exactly ([`combinatorics::arc_span`]). `Ok(None)` where every elementary arc is claimed: the
/// answer is then read nowhere.
///
/// Declines `Width` where a point the question needs cannot be stated in `Rat` (a small model with
/// long decimals); `CylHoleFeature` where two of the nodes are one point
/// ([`combinatorics::CircleOrderFail::Coincident`]); `CylSpan` where their θ order is undecided, a
/// carved end is not among them, the realized middle does not land inside its arc, or the point
/// reads on the boundary or not at all.
fn outer_by_witness(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    cyl: usize,
    coeffs: &[Rat; 4],
    loops: &[combinatorics::LateralLoop],
    carved: &[Carved],
    contacts: &[NodeId],
) -> Result<Option<CylOnClass>, DeclineKind> {
    use combinatorics::{ArcEnds, ArcSpan, ChartPoint};
    use nacre_exact::quad::{CylinderMeet, MeetLine, QuadVal, plane_plane_cylinder};
    use nacre_geom::intersect::RingSide;
    let (o, m) = (def.origin(), def.dir());
    // The trace's two roots on the plane through the axis along `d`.
    let roots = |d: [Rat; 3]| -> Option<(MeetLine, [QuadVal; 2])> {
        let n = nacre_exact::cross3_rat(&m, &d)?;
        let d0 = Rat::from_int(0).checked_sub(nacre_exact::dot3_rat(&n, &o)?)?;
        match plane_plane_cylinder(coeffs, &[n[0], n[1], n[2], d0], &o, &m, def.r2())? {
            CylinderMeet::Pair { line, s } => Some((line, s)),
            _ => None,
        }
    };
    let answer = |x: &ChartPoint<'_>| match combinatorics::loop_parity(jd, def, loops, x) {
        Some(RingSide::Inside) => Ok(Some(CylOnClass::Crosses)),
        Some(RingSide::Outside) => Ok(None),
        None if x.unrealizable() => Err(DeclineKind::Width),
        Some(RingSide::OnBoundary) | None => Err(DeclineKind::CylSpan),
    };
    let mut names: Vec<NodeId> = carved
        .iter()
        .flat_map(|c| c.arc)
        .chain(contacts.iter().copied())
        .collect();
    names.sort_unstable();
    names.dedup();
    if names.is_empty() {
        let realize = || {
            let (line, [s, _]) = roots(seam_direction(def)?)?;
            Some((line, s))
        };
        return answer(&ChartPoint::lazy(*coeffs, &realize));
    }
    let groups = combinatorics::theta_groups(jd, cyl, def, &names).map_err(|e| match e {
        combinatorics::CircleOrderFail::Undecided => DeclineKind::CylSpan,
        combinatorics::CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
    })?;
    let n = groups.len();
    let group_of = |x: NodeId| groups.iter().position(|g| g.iter().any(|&i| names[i] == x));
    let mut claimed = vec![false; n];
    for c in carved {
        let (Some(a), Some(b)) = (group_of(c.arc[0]), group_of(c.arc[1])) else {
            return Err(DeclineKind::CylSpan);
        };
        let mut k = a;
        while k != b {
            claimed[k] = true;
            k = (k + 1) % n;
        }
    }
    let Some(k) = (0..n).find(|&k| !claimed[k]) else {
        return Ok(None);
    };
    let exact = |x: NodeId| {
        combinatorics::pierce_meet(jd, cyl, def, x)
            .and_then(combinatorics::PiercePoint::into_exact)
            .ok_or(DeclineKind::Width)
    };
    let (a, b, whole) = (
        exact(names[groups[k][0]])?,
        exact(names[groups[(k + 1) % n][0]])?,
        n == 1,
    );
    let d = middle_direction(def, &a, &b, whole).ok_or(DeclineKind::Width)?;
    let (line, s) = roots(d).ok_or(DeclineKind::Width)?;
    let ends = if whole {
        ArcEnds::Closed(&a)
    } else {
        ArcEnds::Open(&a, &b)
    };
    let mut inside = None;
    for r in s {
        match combinatorics::arc_span(def, ends, &(line.clone(), r)) {
            Some(ArcSpan::Inside) => inside = Some(r),
            Some(_) => {}
            None => return Err(DeclineKind::Width),
        }
    }
    let s = inside.ok_or(DeclineKind::CylSpan)?;
    answer(&ChartPoint::realized(*coeffs, line, s))
}

/// The seam's direction off the axis, `(m·m)·ref − (ref·m)·m` — rational and perpendicular to `m`,
/// as its smallest integer vector ([`primitive`]).
fn seam_direction(def: &nacre_topo::CylinderDef) -> Option<[Rat; 3]> {
    let (m, r) = (def.dir(), def.ref_dir());
    let (mm, rm) = (
        nacre_exact::dot3_rat(&m, &m)?,
        nacre_exact::dot3_rat(&r, &m)?,
    );
    let mut e = [Rat::from_int(0); 3];
    for k in 0..3 {
        e[k] = mm.checked_mul(r[k])?.checked_sub(rm.checked_mul(m[k])?)?;
    }
    primitive(e)
}

/// **A direction as its smallest integer vector** — denominators cleared, the content divided out.
/// A direction is unchanged by a positive rescale, and the planes built on it stay narrow.
fn primitive(v: [Rat; 3]) -> Option<[Rat; 3]> {
    fn gcd(a: i128, b: i128) -> i128 {
        let (mut a, mut b) = (a.checked_abs().unwrap_or(0), b.checked_abs().unwrap_or(0));
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }
    let mut l = 1i128;
    for c in &v {
        let d = c.denom();
        l = l.checked_mul(d / gcd(l, d).max(1))?;
    }
    let ints = [
        v[0].numer().checked_mul(l / v[0].denom())?,
        v[1].numer().checked_mul(l / v[1].denom())?,
        v[2].numer().checked_mul(l / v[2].denom())?,
    ];
    let g = gcd(gcd(ints[0], ints[1]), ints[2]).max(1);
    Some(ints.map(|x| Rat::from_int(x / g)))
}

/// **The direction off the axis halfway round from `a` to `b`, counter-clockwise** (the whole turn
/// past `a` when `whole`) — realized from the two points' coordinates, a sample for a witness the
/// caller verifies, never a truth.
///
/// Rounded to the grid the arc needs and no finer: a unit direction rounded at `2⁻ᵇ` per component
/// turns by at most `√3·2⁻ᵇ`, inside the half span once `2ᵇ > 2√3 / span`; a finer grid only widens
/// the planes built on it (taken exactly as the `f64` it is, the plane overflowed `i128` on 8,726
/// of the lib suite's traces).
fn middle_direction(
    def: &nacre_topo::CylinderDef,
    a: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    b: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    whole: bool,
) -> Option<[Rat; 3]> {
    let f = |v: [Rat; 3]| v.map(Rat::to_f64);
    let (o, m, e1) = (f(def.origin()), f(def.dir()), f(seam_direction(def)?));
    let unit = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        v.map(|x| x / l)
    };
    let (m, e1) = (unit(m), unit(e1));
    let e2 = [
        m[1] * e1[2] - m[2] * e1[1],
        m[2] * e1[0] - m[0] * e1[2],
        m[0] * e1[1] - m[1] * e1[0],
    ];
    let angle = |p: &(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)| {
        let (base, dir, s) = (f(p.0.base()), f(p.0.dir()), p.1.to_f64_approx());
        let q: [f64; 3] = core::array::from_fn(|k| base[k] + s * dir[k] - o[k]);
        let dot = |u: [f64; 3]| q[0] * u[0] + q[1] * u[1] + q[2] * u[2];
        dot(e2).atan2(dot(e1))
    };
    let tau = std::f64::consts::TAU;
    let (ta, tb) = (angle(a), angle(b));
    let span = if whole {
        tau
    } else {
        (tb - ta).rem_euclid(tau)
    };
    let mid = ta + span / 2.0;
    let d: [f64; 3] = core::array::from_fn(|k| mid.cos() * e1[k] + mid.sin() * e2[k]);
    let bits = ((-(span.log2())).ceil() + 2.0).clamp(2.0, 60.0);
    let scale = 2f64.powf(bits);
    let ints = d.map(|x| (x * scale).round());
    if ints
        .iter()
        .any(|x| !x.is_finite() || x.abs() >= 2f64.powi(62))
    {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    let v = ints.map(|x| Rat::from_int(x as i128));
    primitive(v)
}

/// **A lateral face's shape as its cycles state it** — what both lateral roads read of a
/// lateral's loops: its rims — each a closed rim edge or a rim split into arcs
/// ([`arcs_on_one_class`]) — by plane class with the side of it the face lies on, its two ends — the
/// classes of its lowest and highest boundary sections, ordered as sections are
/// ([`nacre_exact::quad::section_height_order`], with no generator: a face's boundary sections do
/// not cross its ends) — and every other cycle — a panel, a chain rim, a hole — with whether it is
/// an inner loop. The ends are a second derivation of `CylFaceInfo::ends` (the model side's, off the
/// face's own boundary planes), and the two are held equal where they meet.
///
/// Declines by the one name: `None` cycles (a loop could not be named), a class with no world
/// statement, a boundary section crossing an end, a single boundary section, more than two rims,
/// or a rim that is not an end — a whole rim bounds the face, so it is an extreme of it.
pub(super) fn lateral_shape(
    jd: &Judge<'_, WorkingPlane>,
    cf: &crate::planes::CylFaceInfo,
    fl: &combinatorics::FaceLoops,
    def: &nacre_topo::CylinderDef,
) -> Result<LateralShape, DeclineKind> {
    use combinatorics::LoopRing;
    let Some(cycles) = &fl.cycles else {
        return Err(DeclineKind::CylSpan);
    };
    let mut rims: Vec<usize> = Vec::new();
    let mut rest: Vec<(bool, LoopRing)> = Vec::new();
    let mut sections: Vec<usize> = Vec::new();
    let mut note = |c: usize| {
        if !sections.contains(&c) {
            sections.push(c);
        }
    };
    for (i, ring) in cycles.iter().enumerate() {
        match ring {
            LoopRing::Rim { plane } => {
                note(*plane);
                rims.push(*plane);
            }
            ring => {
                // **A rim split into arcs is a rim.** Corners on the rim's own plane (a pin over a
                // plate's edge: the edge's wall meets the rim there) cut the closed rim edge into
                // arcs, and the loop arrives as a ring — but a simple loop of arcs of one section
                // goes once round it, so it is the whole rim, as `LoopRing::Rim` is. Read as a hole
                // or a panel it lies wholly on its own class and has nothing to bound with. The
                // ring stays as the model has it (`combinatorics::FaceLoops::cycles`): «is this a
                // rim» is a question about the face's shape, asked here, where both lateral roads
                // read it — a different question from the gate's whole rim
                // (`planes::cyl_gate::only_whole_rims`).
                if let Some(nr) = ring.poly()
                    && let Some(plane) = arcs_on_one_class(&nr.walls)
                {
                    note(plane);
                    rims.push(plane);
                    continue;
                }
                // The cycle's arcs each lie on a boundary section; its rulings are none.
                if let Some(nr) = ring.poly() {
                    for w in &nr.walls {
                        if let crate::combinatorics::Wall::Arc { plane, .. } = *w {
                            note(plane);
                        }
                    }
                }
                rest.push((i > 0, ring.clone()));
            }
        }
    }
    let planes = sections
        .iter()
        .map(|&c| jd.world_rat(c).ok_or(DeclineKind::CylSpan))
        .collect::<Result<Vec<_>, _>>()?;
    let [lo, hi] = crate::planes::section_ends(&planes, def)
        .map(|[l, h]| [sections[l], sections[h]])
        .ok_or(DeclineKind::CylSpan)?;
    if lo == hi || rims.len() > 2 || rims.iter().any(|&c| c != lo && c != hi) {
        return Err(DeclineKind::CylSpan);
    }
    if let [a, b] = rims[..]
        && a == b
    {
        return Err(DeclineKind::CylSpan);
    }
    debug_assert_eq!(
        jd.world_rat(lo).zip(jd.world_rat(hi)).map(|(l, h)| [l, h]),
        cf.ends,
        "the cycles' ends by name and the face's ends by model disagree"
    );
    Ok(LateralShape {
        rims: rims.into_iter().map(|c| (c, c == lo)).collect(),
        sections,
        cycles: rest,
    })
}

/// The one plane class every wall of a ring is an arc on, `None` when a wall is not an arc or two
/// arcs ride different classes. A lateral's arcs are its own cylinder's, so they lie on one
/// circle, and a simple loop of arcs of one circle runs one way round it — asserted, as the
/// producer's premise.
fn arcs_on_one_class(walls: &[crate::combinatorics::Wall]) -> Option<usize> {
    use crate::combinatorics::Wall;
    let Some(&Wall::Arc { plane, .. }) = walls.first() else {
        return None;
    };
    if !walls
        .iter()
        .all(|w| matches!(*w, Wall::Arc { plane: c, .. } if c == plane))
    {
        return None;
    }
    debug_assert!(
        walls.windows(2).all(|w| match (w[0], w[1]) {
            (Wall::Arc { ccw: d0, .. }, Wall::Arc { ccw: d1, .. }) => d0 == d1,
            _ => false,
        }),
        "a rim split into arcs runs one way round: {walls:?}"
    );
    Some(plane)
}

/// **Lay the carved extents over the outer answer and hand back *exclusive* spans.**
///
/// **Exclusive, not "the whole circle plus overrides".** The cheaper spelling would emit the
/// outer answer over the whole circle and let the hole's arcs sit on top, leaning on `edge_mask`'s
/// `Graze > Transversal` precedence to pick the winner. That precedence is for **different faces
/// meeting on one edge**, not for one face's own extent — borrowing it would be a category error,
/// and it cannot express "absent" at all. So the circle is cut at every boundary and each piece
/// carries exactly one answer.
///
/// Adjacent pieces with the same answer are **merged**, so a boundary where nothing changes leaves
/// no cut behind: an extra split point would be a degree-2 vertex the winding walk has to have a
/// turn for.
fn assemble_spans(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    carved: &[Carved],
    outer: Option<CylOnClass>,
    ties: &mut Vec<(NodeId, NodeId)>,
) -> Result<Vec<CircleSpan>, DeclineKind> {
    let mut names: Vec<NodeId> = carved.iter().flat_map(|c| c.arc).collect();
    names.sort_unstable();
    names.dedup();
    // One θ order of the circle's nodes. Names at one angle of this circle are one point: each is
    // read as its group's first, and the identity goes to the round's table (`ties`).
    let groups = combinatorics::theta_groups(jd, cyl, def, &names).map_err(|e| match e {
        combinatorics::CircleOrderFail::Undecided => DeclineKind::CylSpan,
        combinatorics::CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
    })?;
    let mut first: HashMap<NodeId, NodeId> = HashMap::new();
    for g in &groups {
        for &i in &g[1..] {
            ties.push((names[g[0]], names[i]));
            first.insert(names[i], names[g[0]]);
        }
    }
    let one = |n: NodeId| first.get(&n).copied().unwrap_or(n);
    let carved: Vec<Carved> = carved
        .iter()
        .map(|c| Carved {
            arc: c.arc.map(one),
            on: c.on,
        })
        .collect();
    let carved = carved.as_slice();
    // The points in θ order, and the same points sorted by name with `order` between the two.
    let in_theta: Vec<NodeId> = groups.iter().map(|g| names[g[0]]).collect();
    let mut nodes = in_theta.clone();
    nodes.sort_unstable();
    let order: Vec<usize> = in_theta
        .iter()
        .map(|x| nodes.binary_search(x).expect("a group's first is a point"))
        .collect();
    let n = order.len();
    // θ position of each node, keyed by its index in the sorted `nodes`.
    let mut at = vec![0usize; n];
    for (k, &i) in order.iter().enumerate() {
        at[i] = k;
    }
    let place = |x: NodeId| {
        at[nodes
            .binary_search(&x)
            .expect("every boundary node is in the set")]
    };
    // The answer over elementary arc `k` (from `order[k]` to `order[k+1]`), or `None` where no
    // cycle reaches: there the outer answer stands — which may itself be "absent".
    let mut ans: Vec<Option<Option<CylOnClass>>> = vec![None; n];
    for c in carved {
        let (a, b) = (place(c.arc[0]), place(c.arc[1]));
        if a == b {
            return Err(DeclineKind::CylHoleFeature); // an extent of no length, or of the whole circle
        }
        let mut k = a;
        while k != b {
            if ans[k].is_some() {
                return Err(DeclineKind::CylHoleFeature); // two cycles claiming one extent
            }
            ans[k] = Some(c.on);
            k = (k + 1) % n;
        }
    }
    let answer = |k: usize| ans[k].unwrap_or(outer);
    // One answer everywhere: the circle is undivided, whatever the holes touch.
    let Some(start) = (0..n).find(|&k| answer(k) != answer((k + n - 1) % n)) else {
        return Ok(match answer(0) {
            Some(on) => vec![CircleSpan { arc: None, on }],
            None => Vec::new(),
        });
    };
    let mut out = Vec::new();
    let mut k = 0;
    while k < n {
        let i = (start + k) % n;
        let a = nodes[order[i]];
        let this = answer(i);
        let mut len = 1;
        while k + len < n && answer((start + k + len) % n) == this {
            len += 1;
        }
        if let Some(on) = this {
            out.push(CircleSpan {
                arc: Some([a, nodes[order[(start + k + len) % n]]]),
                on,
            });
        }
        k += len;
    }
    Ok(out)
}

/// **Which way round a lateral face's material lies, for a boundary edge travelling along the
/// axis** — the one place that sign is made, and the third application of one convention.
///
/// *Material is on the left of the ring's direction of travel* — the sentence
/// `trace_plane::run_body_above` states for a planar face's on-line run. On a lateral face, with
/// `r̂` the outward radial direction, `m̂` the axis and `θ̂` increasing θ (counter-clockwise about
/// `m̂`), cylindrical coordinates are right-handed so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The
/// face's outward is `σ·r̂` with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], travel is
/// `τ·m̂`, and left of travel is
///
/// ```text
///   n_out × travel = (σ·r̂) × (τ·m̂) = σ·τ·(r̂ × m̂) = −σ·τ·θ̂
/// ```
///
/// so the answer is `−σ·τ`: the θ direction the face occupies. Two readers need it — the roads
/// across the axis, where a ruling edge crossing the class or meeting it at a corner puts the
/// **hole** in the opposite θ direction ([`cycle_on_class`], [`corner_on_section`]), and the ∥
/// road, where a ruling edge lying on the class puts the **face** to one side of the wall — and
/// they must not spell it twice.
///
/// Flipping this sign turns the ⊥ road's lock and the ∥ road's lock red together.
pub(super) fn material_theta_sign(orient_sign: i8, travel_up: i8) -> i8 {
    -(orient_sign * travel_up)
}

/// The θ side a cycle's arc lies on, seen from the node where it meets a ruling: an arc arriving
/// counter-clockwise came from smaller θ (`−1`), one departing counter-clockwise goes to larger θ
/// (`+1`) — the flank of an on-line run, read from the arc's own sense (`Wall::Arc`'s `ccw`)
/// rather than from a coordinate, because a global "side" does not exist on a cylinder.
pub(super) fn arc_flank(ccw: bool, arriving: bool) -> i8 {
    if ccw == arriving { -1 } else { 1 }
}

/// **A corner of a lateral face's boundary — a ruling meeting an arc — lying on an oblique
/// section**: `(hole, convex)`, the θ direction the face is absent in beside the ruling (`hole`,
/// `±1`) and whether the face fills only the quarter between the two edges there (`convex`).
///
/// In the chart `(θ, z)` the ruling runs along `z` and the arc along `θ`, so they meet square: the
/// face at the corner is the quarter the two edges span or its complement. The ruling's own travel
/// (`Wall::Ruling`'s `up`, the producer's) puts the material on the θ side
/// [`material_theta_sign`] names, and the arc leaves the corner on the θ side [`arc_flank`] names;
/// the face is that quarter exactly when the two agree.
///
/// Where the section passes the corner transversally, one half of it runs into the quarter the
/// edges span or into the opposite one when the ring crosses it there (the walk's flanks differ),
/// and into the two side quarters when it touches. Crossing, exactly one half is in the face, and
/// on both kinds of corner it is the half away from `hole`; touching, both halves are out of a
/// convex face and in a concave one. Where the section is level at the corner (the plane holds the
/// rim circle's tangent there), its two halves bend to one side of the arc, so the arc's far end
/// lies on the other side of the plane; the flanks still say crossing or touch, and the same
/// reading holds.
pub(super) fn corner_on_section(
    orient_sign: i8,
    ruling_up: bool,
    arc_ccw: bool,
    arc_arriving: bool,
) -> (i8, bool) {
    let material = material_theta_sign(orient_sign, if ruling_up { 1 } else { -1 });
    (-material, material == arc_flank(arc_ccw, arc_arriving))
}

/// **Which way `+θ̂` points across a ∥ wall, in the frame a label is written in** — `+1` when
/// leaving a ruling counter-clockwise about the axis enters the side the wall's plane faces.
///
/// **One atom, three readers.** [`combinatorics::RulingCarrier::side`] is
/// `sign((x − o) · (m̂ × n̂_r))` against the class's *rational* name, and the scalar triple product
/// gives `(x − o) · (m̂ × n̂) = n̂ · ((x − o) × m̂) = −r·(n̂ · θ̂)`, so `sign(n̂_r · θ̂) = −side`;
/// the world name's sense (`κ`, [`crate::planes::WorldName`] — the truth's, since the name alone
/// carries no direction and the plane cache is a rounded image) carries that to the plane's own
/// facing. The three consumers are the
/// ruling sweep's graze side, [`ruling_interior_is_even`], and the chart's vertical read — and the
/// whole point of naming it is that none of them spells the product a second time. Two of them
/// read a **label** and stay in the stored frame; the one that picks a **cell** —
/// [`ruling_interior_is_even`] — crosses into the chart's frame with `frame_sign`, and it is the
/// only reader that does.
pub(crate) fn plus_theta_is_above(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    side: i8,
) -> Option<bool> {
    // A tangent ruling (`side == 0`) is never a recorded crossing, so nothing labelled by this
    // sign reaches it; the day one does, the answer is a measurement, not a sign.
    debug_assert_ne!(side, 0, "a tangent ruling has no +θ side to be above");
    let kappa = jd.canon_sense(wc);
    Some(-(kappa * side) == 1)
}

/// **One boundary cycle — a hole, a panel, a chain rim — read against one class crossing the
/// axis**: the walk's features turned into extents against the class's section (`trace` the class's
/// plane; [`trace_on_face`]) — a circle across the axis, an ellipse oblique to it. A hole is
/// carved out of a whole-trace outer answer; a panel or a
/// chain is carved out of one that may be absent — the same arms, because neither reads
/// "inside" or "outside": only the ring's winding, below.
///
/// **Which way round is decided by the ring's own winding, and by nothing else.** A plane cuts a
/// circle in two points, and "which of the two arcs is the hole" is the whole difficulty: a hole
/// whose two ruling edges lie on **one** plane cannot be told apart by that plane's sides. The
/// universal boundary convention answers it locally instead — *material is on the left of the
/// ring's direction of travel* — the same sentence `trace_plane::run_body_above` states for a
/// planar face's on-line run, and it holds for an outer ring, a hole ring, a notch and a reflex
/// corner alike.
///
/// Not "the side opposite the carrier's outward normal": the ring's winding needs no opposite
/// face, no `edge → face` table and no adjacency cache, and the same walk covers a stepped hole
/// and a U-shaped one.
///
/// **Derivation.** Write `r̂` for the cylinder's outward radial direction, `m̂` for its axis and `θ̂`
/// for increasing θ (counter-clockwise about `m̂`, which is what
/// `nacre_exact::quad::circular_order_about_seam` ranks and what an arc [`Piece`]'s `end` runs along).
/// Cylindrical coordinates are right-handed, so `r̂ × θ̂ = m̂` and `r̂ × m̂ = −θ̂`. The face's outward
/// is `σ·r̂` with `σ = ` [`crate::planes::CylFaceInfo::orient_sign`], and "left of travel" is
/// `n_out × travel`.
///
/// - **A ruling edge crossing the class** travels `τ·m̂` — `τ` is the ring's own travel along the
///   ruling, `Wall::Ruling`'s `up`. Material lies along `(σ·r̂) × (τ·m̂) = −σ·τ·θ̂`
///   ([`material_theta_sign`]), so the **hole** lies along `+σ·τ·θ̂`: the hole runs
///   counter-clockwise from the crossing exactly when `σ·τ = +1`.
/// - **An arc edge lying on the class** (one of the hole's own rims) travels `ν·θ̂`. Material lies
///   along `(σ·r̂) × (ν·θ̂) = σ·ν·m̂`, so `ν = σ·μ` where `μ` is the axis direction the face occupies
///   — which is the opposite of the side the run's flanks sit on. That fixes the rim arc's CCW
///   orientation without reading a coordinate.
/// - **An arc edge crossing the trace** travels `ν·θ̂` too, so material lies along
///   `σ·ν·m̂`; the section, going counter-clockwise from the crossing, climbs along `m̂` by the sign
///   of its tangent's axial part ([`climb_over`]) — into the material when the two agree, so the
///   hole then runs clockwise.
/// - **A corner on the trace** — a node of the ring on the class, where a ruling meets an
///   arc — is a crossing when the walk's flanks differ and a touch when they agree; the ruling's
///   material side and the arc's side read which ([`corner_on_section`]). A crossing is a cut at the
///   corner's own name, the hole on the ruling's empty side; a touch outside a convex corner takes
///   nothing away.
///
/// `up` is `plus_t_is_above`: whether the class's **stored** normal points along `+m̂` — the frame
/// a graze's `body_above` is a label in.
#[allow(clippy::too_many_arguments)]
fn cycle_on_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cf: &crate::planes::CylFaceInfo,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    cyl: usize,
    up: bool,
    inner: bool,
    nr: &combinatorics::NamedRing,
    trace: &[Rat; 4],
    out: &mut Vec<Carved>,
    contacts: &mut Vec<NodeId>,
) -> Result<(), DeclineKind> {
    use crate::combinatorics::Wall;
    let _ = inner; // the ledger's, under test
    #[cfg(test)]
    let start = out.len();
    let ring = &nr.triples;
    let n = ring.len();
    // **Every edge is read by its carrier against the class's section.** An arc lies in its
    // section's plane: the class's own plane is that section — the arc lies on the trace (`On`), a
    // run the arm below lays as a graze — and any other plane meets the class in a line, the arc
    // crossing the trace where it holds a root of that line ([`combinatorics::arc_crossings`], the
    // count the planar trace reads too; a plane parallel to the class crosses it nowhere); with both
    // ends on the class it leaves to the side its tangent points ([`combinatorics::arc_departure`]).
    // An open arc cannot leave along it — a line through two points of a section is no tangent — and
    // a loop of one section's arcs is a rim, read before this road; a closed arc inside a longer
    // ring that leaves along it touches the trace there, which the walk reads as a touch. A ruling
    // is a line along the axis and meets the plane once, at a sign change of its ends — a straight
    // edge, `0` inside. An end on the class is a corner on the trace, the root under the end's own
    // name, which the arc does not cross again there.
    //
    // The walk's count and the `Crossing` arm's node are one question — asked twice, the second
    // asking without the end would meet the end's root as one point under two names and refuse.
    let arc_ends = |i: usize| [ring[i], ring[(i + 1) % n]];
    let arc_roots = |i: usize, on_end: Option<NodeId>| match nr.walls[i] {
        Wall::Arc { cyl: c, ccw, plane } => {
            let def = &cyls.get(c)?.def;
            combinatorics::arc_crossings(jd, def, wc, plane, c, ccw, arc_ends(i), on_end)
        }
        Wall::Ruling { .. } | Wall::Plane(_) => None,
    };
    let features = match combinatorics::ring_against_plane(
        jd,
        cyls,
        ring,
        wc,
        |i| match nr.walls[i] {
            Wall::Arc { plane, .. } if plane == wc => Some(combinatorics::EdgeMeet::On),
            Wall::Arc { cyl: c, ccw, plane } => {
                combinatorics::arc_departure(jd, cyls, ring[i], plane, wc, c, ccw)
            }
            Wall::Ruling { .. } | Wall::Plane(_) => None,
        },
        |i, on_end| match nr.walls[i] {
            Wall::Arc { plane, .. } if plane == wc => Some(0),
            Wall::Arc { .. } => arc_roots(i, on_end.map(|k| ring[k]))
                .and_then(|roots| u8::try_from(roots.len()).ok()),
            Wall::Ruling { .. } | Wall::Plane(_) => Some(0),
        },
    ) {
        combinatorics::RingWalk::Met(f) => f,
        // A hole ring lying wholly in the class has no thickness to bound anything with (a rim
        // split into arcs lies wholly in its class too, and `lateral_shape` reads it as a rim
        // before this road), and a node the walk cannot name is the same refusal every other
        // consumer makes of it.
        combinatorics::RingWalk::AllOn | combinatorics::RingWalk::Unnameable => {
            return Err(DeclineKind::CylFaceHole);
        }
    };
    let sigma = i32::from(cf.orient_sign);
    // The nodes on the class: an arc crossing the trace with one of them as an end was counted with
    // that end (`arc_roots`).
    let on_class: Vec<usize> = features
        .iter()
        .filter_map(|f| match *f {
            combinatorics::Feature::Run { first, .. } => Some(first),
            combinatorics::Feature::Crossing { .. } => None,
        })
        .collect();
    contacts.extend(on_class.iter().map(|&i| ring[i]));
    // Crossings, each with the direction the hole runs in from it. Paired below, after the θ sort.
    let mut cuts: Vec<(NodeId, bool)> = Vec::new();
    for f in features {
        match f {
            combinatorics::Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } if len < 2 => {
                // **A corner of the boundary on the trace**, one node: the walk's flanks say
                // whether the ring crosses the trace there or touches it. A ruling meeting an arc is
                // read by [`corner_on_section`] — a crossing is a cut at the corner's own name, the
                // hole on the ruling's empty side; a touch outside a convex corner takes nothing
                // away, and one from inside the face is past this road. Any other corner has no
                // fixture: crossing there is declined rather than guessed at, and a touch takes no
                // extent away.
                let corner = match (nr.walls[(first + n - 1) % n], nr.walls[first]) {
                    (
                        Wall::Ruling {
                            cyl: rc, up: r_up, ..
                        },
                        Wall::Arc { cyl: ac, ccw, .. },
                    ) if rc == cyl && ac == cyl => {
                        Some(corner_on_section(cf.orient_sign, r_up, ccw, false))
                    }
                    (
                        Wall::Arc { cyl: ac, ccw, .. },
                        Wall::Ruling {
                            cyl: rc, up: r_up, ..
                        },
                    ) if rc == cyl && ac == cyl => {
                        Some(corner_on_section(cf.orient_sign, r_up, ccw, true))
                    }
                    _ => None,
                };
                match (flanks_differ, corner) {
                    (true, Some((hole, _))) => cuts.push((ring[first], hole > 0)),
                    (true, None) | (false, Some((_, false))) => {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    (false, Some((_, true)) | None) => {}
                }
            }
            combinatorics::Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } => {
                // A run whose flanks differ is a rim the hole continues *through* — it is both an
                // extent and a parity toggle, and no fixture makes one. Declined and counted
                // rather than guessed at.
                if flanks_differ {
                    return Err(DeclineKind::CylHoleFeature);
                }
                // **Which way the arc runs, and so which side the face is on, is the
                // producer's to say — not the flank's.** The flank (the side the ring's
                // off-class neighbours sit on) gives the ring's *interior* side — right for a
                // convex hole and wrong for a wrapping rim, whose highest arc has both
                // neighbours below it and the face below too. The arc's own direction about the
                // axis (`Wall::Arc`'s `ccw`, the
                // producer's convention) settles both: walked with sense `ν`, material lies
                // along `σ·ν·m̂`, so the face is above the class exactly when that agrees with
                // the class's stored normal (`up`), and the carved extent is the arc as walked.
                // Measured: this reading and the flank's are equal on every on-class arc of
                // today's holes (50 of 50).
                for k in 0..len - 1 {
                    // **"Both ends on the class" is not "the edge is on the trace."** The walk
                    // reads a run off its *nodes*; on a plane the edge between two on-line nodes
                    // follows, and on a **cylinder** it does not — an arc of another section can
                    // meet this class at both ends without lying on it. Taking that for an arc of
                    // the trace would state an extent along a curve that is not there. The carrier
                    // says it directly: an edge on the trace is an arc whose section's plane is that
                    // class. The check belongs here, per edge, where the answer is used.
                    let edge = (first + k) % n;
                    let crate::combinatorics::Wall::Arc { ccw: nu, plane, .. } = nr.walls[edge]
                    else {
                        return Err(DeclineKind::CylHoleFeature);
                    };
                    if plane != wc {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    let body_above = (sigma * if nu { 1 } else { -1 } > 0) == up;
                    let a = ring[edge];
                    let b = ring[(edge + 1) % n];
                    let arc = if nu { [a, b] } else { [b, a] };
                    out.push(Carved {
                        arc,
                        on: Some(CylOnClass::Grazes { body_above }),
                    });
                    #[cfg(test)]
                    arc_probe::record(jd, cyl, def, inner, &out[start..]);
                }
            }
            combinatorics::Feature::Crossing { edge, nth, .. } => {
                // The crossed edge's carrier, **carried** from the producer rather than re-derived
                // from the two endpoint names.
                //
                // **An arc crossing the trace** is the pierce node `wc ∩ plane ∩ cyl` the
                // arc holds — the `nth` in its travel ([`combinatorics::arc_crossings`]), the name
                // the planar trace gives the cap's own crossing there. The arc travels `ν·θ̂` and
                // material lies along `(σ·r̂) × (ν·θ̂) = σ·ν·m̂`; the section, going
                // counter-clockwise from the node, climbs along `m̂` by the sign of its tangent's
                // axial part ([`climb_over`]) — into the material when the two agree, so the
                // hole then runs clockwise.
                if let Wall::Arc {
                    cyl: c, ccw, plane, ..
                } = nr.walls[edge]
                {
                    if c != cyl {
                        return Err(DeclineKind::CylHoleFeature);
                    }
                    let end = [edge, (edge + 1) % n]
                        .into_iter()
                        .find(|k| on_class.contains(k))
                        .map(|k| ring[k]);
                    let cut = arc_roots(edge, end)
                        .and_then(|roots| roots.get(usize::from(nth)).copied())
                        .ok_or(DeclineKind::CylHoleFeature)?;
                    let mu = sigma * if ccw { 1 } else { -1 };
                    let curve = jd.world_rat(plane).ok_or(DeclineKind::ObliqueSection)?;
                    let climb = climb_over(jd, cyl, def, trace, &curve, cut)?;
                    cuts.push((cut, i32::from(climb) != mu));
                    continue;
                }
                // The ruling is read off the wall alone — its cylinder, its side and the plane the
                // side was measured against — not off the corner's name: a corner on three planes
                // is named by the pair its definition chose, which need not hold `j`.
                let crate::combinatorics::Wall::Ruling {
                    cyl: ncyl,
                    side,
                    up: r_up,
                    plane: j,
                } = nr.walls[edge]
                else {
                    return Err(DeclineKind::CylHoleFeature);
                };
                // **The crossing is the *same ruling*, restated for this class — through
                // the one door the planar scan uses** ([`combinatorics::crossing_on_ruling`]): the ruling's side,
                // measured against its wall `j` at the hole's own corner — the side the wall
                // carries (`curved_wall`) — picks the root of `{wc, j}` on the cylinder; restating
                // the corner's root by the axis senses of the two ⊥ classes would be a second
                // spelling of the same point. `j` must hold a ruling (a wall plane within the
                // radius) and `wc` cross the axis, which the door checks; anything else is not this
                // feature.
                let wdef = &cyls.get(ncyl).ok_or(DeclineKind::CylHoleFeature)?.def;
                let cut = combinatorics::crossing_on_ruling(jd, wdef, wc, j, ncyl, side)
                    .map_err(|_| DeclineKind::CylHoleFeature)?;
                // Which way the hole runs from here: the face occupies `material_theta_sign` of the
                // ruling's own travel (the producer's, as for an arc), and the hole is the other way.
                let hole_ccw = material_theta_sign(cf.orient_sign, if r_up { 1 } else { -1 }) < 0;
                cuts.push((cut, hole_ccw));
            }
        }
    }
    #[cfg(test)]
    let before = out.len();
    carve_between(jd, cyl, def, &cuts, out)?;
    #[cfg(test)]
    for k in before..out.len() {
        arc_probe::record(jd, cyl, def, inner, &out[start..=k]);
    }
    Ok(())
}

/// **Pair one closed boundary's crossings of the trace into the extents it takes away** — `cuts`
/// is each crossing with whether the absent side runs counter-clockwise from it. A closed boundary
/// enters and leaves the trace equally often and the crossings alternate in θ, so sorting them and
/// pairing each "absent ahead" with the next is the whole assignment — the circle twin of the sweep
/// along a line ([`sweep_line`]).
fn carve_between(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    cuts: &[(NodeId, bool)],
    out: &mut Vec<Carved>,
) -> Result<(), DeclineKind> {
    if cuts.len() % 2 != 0 {
        return Err(DeclineKind::CylHoleFeature);
    }
    if cuts.is_empty() {
        return Ok(());
    }
    let nodes: Vec<NodeId> = cuts.iter().map(|c| c.0).collect();
    let order = combinatorics::circular_order(jd, cyl, def, &nodes).map_err(|e| match e {
        combinatorics::CircleOrderFail::Undecided => DeclineKind::CylSpan,
        combinatorics::CircleOrderFail::Coincident => DeclineKind::CylHoleFeature,
    })?;
    for k in 0..order.len() {
        let (a, ahead) = cuts[order[k]];
        let (b, next_ahead) = cuts[order[(k + 1) % order.len()]];
        if ahead == next_ahead {
            return Err(DeclineKind::CylHoleFeature); // the boundaries did not alternate
        }
        if ahead {
            out.push(Carved {
                arc: [a, b],
                on: None,
            });
        }
    }
    Ok(())
}

/// **Which way the trace leaves a boundary curve it crosses at `node`** — `+1` when, travelling
/// counter-clockwise about the axis from the node, the trace (the class's section, plane `trace`)
/// stands above the curve (the section of plane `curve`) along `+m` (the raw `def.dir()`, the
/// direction the axis parameter grows in), `−1` below: the order of the two sections just after the
/// generator through the node ([`nacre_exact::quad::section_height_order`], the one spelling the
/// chart reads). Never level where the curve crosses the trace transversally — a touch there,
/// which `planes::boundary_stand` refuses at a boundary and [`combinatorics::arc_crossings`] reads as
/// no crossing elsewhere — so anything else declines rather than guessing.
fn climb_over(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    trace: &[Rat; 4],
    curve: &[Rat; 4],
    node: NodeId,
) -> Result<i8, DeclineKind> {
    use nacre_exact::quad::{SectionOrder, Toward, section_height_order};
    let (line, s) = combinatorics::pierce_meet(jd, cyl, def, node)
        .and_then(combinatorics::PiercePoint::into_exact)
        .ok_or(DeclineKind::ObliqueSection)?;
    match section_height_order(trace, curve, def, Some((&line, &s)), Toward::After) {
        Some(SectionOrder::Ordered(core::cmp::Ordering::Greater)) => Ok(1),
        Some(SectionOrder::Ordered(core::cmp::Ordering::Less)) => Ok(-1),
        _ => Err(DeclineKind::ObliqueSection),
    }
}
