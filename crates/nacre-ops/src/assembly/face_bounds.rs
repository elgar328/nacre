//! **The result contract's arm for a plane face's own bounds** — [`face_bounds_meet`]; asked after
//! the self-touch verdict ([`self_touch_reject`]), which leaves a face's own bounds out.

use super::*;
use nacre_geom::intersect::RingSide;

/// **One piece of a plane face's boundary** — a straight edge on the line `q ∩ wall`, an arc of
/// the face's own section of a cylinder (`ends` counter-clockwise about its axis), a ruling of a
/// cylinder lying in the face's plane (`side` which of the plane's two, `0` a tangent plane's one),
/// or a whole circle.
#[derive(Clone, Copy, Debug)]
enum BoundPiece {
    Segment {
        ends: [NodeId; 2],
        wall: usize,
    },
    Arc {
        ends: [NodeId; 2],
        cyl: usize,
    },
    Ruling {
        ends: [NodeId; 2],
        cyl: usize,
        side: i8,
    },
    Circle {
        cyl: usize,
    },
}

impl BoundPiece {
    fn ends(&self) -> &[NodeId] {
        match self {
            BoundPiece::Segment { ends, .. }
            | BoundPiece::Arc { ends, .. }
            | BoundPiece::Ruling { ends, .. } => ends,
            BoundPiece::Circle { .. } => &[],
        }
    }

    fn cyl(&self) -> Option<usize> {
        match *self {
            BoundPiece::Segment { .. } => None,
            BoundPiece::Arc { cyl, .. }
            | BoundPiece::Ruling { cyl, .. }
            | BoundPiece::Circle { cyl } => Some(cyl),
        }
    }

    /// A straight piece's ends — a segment's or a ruling's, both a line between them.
    fn straight_ends(&self) -> Option<[NodeId; 2]> {
        match *self {
            BoundPiece::Segment { ends, .. } | BoundPiece::Ruling { ends, .. } => Some(ends),
            BoundPiece::Arc { .. } | BoundPiece::Circle { .. } => None,
        }
    }

    /// The same point set as `other` — two straight pieces with the same two ends, a segment and a
    /// ruling alike, the same arc, the same circle. Two arcs of one circle sharing both ends are the
    /// circle's two halves unless they run alike.
    fn same_set(&self, other: &BoundPiece) -> bool {
        if let (Some(a), Some(b)) = (self.straight_ends(), other.straight_ends()) {
            return a == b || a == [b[1], b[0]];
        }
        match (*self, *other) {
            (BoundPiece::Arc { ends: a, cyl: c }, BoundPiece::Arc { ends: b, cyl: d }) => {
                c == d && a == b
            }
            (BoundPiece::Circle { cyl: c }, BoundPiece::Circle { cyl: d }) => c == d,
            _ => false,
        }
    }
}

/// **Where a face's boundary was found meeting itself** — at a node, or along a piece named by its
/// ends. Ordered so the witness is the same run after run.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Meet {
    Node(NodeId),
    Along([NodeId; 2]),
    /// A bound with no node to name it by (a whole circle, a rim).
    Unplaced,
}

/// **A plane face of the result is one connected region whose boundary meets itself only where one
/// loop passes a vertex more than once** — the result contract's arm for a face's own bounds, which
/// [`self_touch_reject`] leaves out (it asks what a face does not bound). Called after it: a body
/// whose surface touches itself is named by that verdict (`Impossible` — a 45° fold merges into a
/// face of two lobes, `assembly::coplanar`), and only a body it passes is asked here. On such a body
/// the arrangement cut every curve where another meets it, so a meet anywhere else is the kernel's
/// own failure (`FaceBoundsMeet`, `SuspectedDefect`).
///
/// Asked of each plane face, exactly, on the names:
/// 1. **A node on two bounds** of the face is a meet. Where one loop passes a node more than once,
///    every arrival there leaves by the way out the face's walk takes ([`combinatorics::way_out`] —
///    the first departure clockwise from the arrival turned back, the face on the walk's left):
///    then the loop is a boundary cycle of one face of the arrangement its own edges make, and the
///    face is one region pinched there. A notch whose tip touches the outer ring's vertex passes; two
///    lobes meeting at a point, an island in a hole, and a loop crossing itself at the vertex each
///    have an arrival whose way out is another visit's departure. This is the local form of the
///    rule exact arrangement libraries check a face's boundary by (CGAL's «every boundary half-edge
///    has the one face on its side»; JTS splits it into a crossing test and a material-side test),
///    and it counts nothing — splitting the loop into simple loops and counting positive windings
///    cannot tell a touch from a crossing. Its premise is the loops' roles: an outer loop wound
///    `+1` and a hole `−1` in the class frame, the face on each one's left, which is how the
///    arrangement keeps them — it sorts its rings by [`combinatorics::loop_winding`]. Whether a hole lies
///    inside the outer ring is the arrangement's nesting, not asked again.
/// 2. **A node inside a piece** it does not end — on a straight piece's closed span along its line,
///    on an arc's or circle's circle strictly inside the arc, on a ruling's closed span along the
///    axis, or at a curve's end under another name (two names, one point) — is a meet, and so are
///    two pieces that are one point set. This is every meet at an end, a T or a partial overlap,
///    so what is left is two interiors.
/// 3. **Two interiors** — pieces whose proven boxes overlap: two segments crossing (all four ends
///    strictly across); a segment and an arc or circle meeting inside both ([`Judge::line_meet`],
///    the predicate the arrangement cuts circles with); a segment and a ruling crossing inside both.
///    Two curves of one cylinder share a circle or are parallel lines, so only their ends meet
///    (part 2). **Two curves of two cylinders are not asked: on a result they never meet**, and
///    neither does a node of one cylinder with a curve of another. A plane face's curve of a
///    cylinder is an edge it shares with that cylinder's lateral, and a result's lateral is part of
///    an operand's. Two cylinders of two operands: the gate's pair rule proves their laterals apart
///    exactly or refuses the pair (`CylinderPairContact`). Two of one operand: two faces of a solid
///    meet only at its edges and vertices — here an edge on both cylinders, which the boolean's
///    door refuses (`CylinderPairEdge`), or a vertex on both, which a vertex names by both surfaces
///    (`Vertex::CylinderPair`) and the kernel builds only at such an edge. The day an operand with
///    one is taken in, this part must ask these pairs.
///
/// What cannot be decided is held, as the self-touch arms hold it: a proven meet anywhere is the
/// answer, and only a face with none raises the held failure (`WitnessNotRational` — the statement a
/// question needed is not there).
pub(crate) fn face_bounds_meet(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    pt: &NodeReach,
    faces: &[LocalFace],
) -> Result<(), BoolError> {
    let mut found: Option<Meet> = None;
    let mut held: Option<Meet> = None;
    let mut found_at = |m: Meet| found = Some(found.map_or(m, |f| f.min(m)));
    let mut hold = |m: Meet| held = Some(held.map_or(m, |h| h.min(m)));
    // One face at a time, into buffers the faces share.
    let mut pieces: Vec<BoundPiece> = Vec::new();
    let mut bound_of: HashMap<NodeId, usize> = HashMap::new();
    let mut seen: HashSet<NodeId> = HashSet::new();
    let mut repeats: Vec<NodeId> = Vec::new();
    let mut boxes: Vec<crate::planes::CacheReach> = Vec::new();
    for lf in faces {
        let ClassIx::Plane(q) = lf.surf else {
            continue;
        };
        // The pieces, and which bound passes each node.
        pieces.clear();
        bound_of.clear();
        for (bi, b) in std::iter::once(&lf.outer).chain(&lf.inner).enumerate() {
            match b {
                Bound::Ring(r) => {
                    let k = r.nodes.len();
                    seen.clear();
                    repeats.clear();
                    for i in 0..k {
                        let u = r.nodes[i];
                        if !seen.insert(u) && !repeats.contains(&u) {
                            repeats.push(u);
                        }
                        match bound_of.get(&u) {
                            Some(&other) if other != bi => found_at(Meet::Node(u)),
                            _ => {
                                bound_of.insert(u, bi);
                            }
                        }
                        let ends = [u, r.nodes[(i + 1) % k]];
                        pieces.push(match r.walls[i] {
                            Wall::Plane(wall) => {
                                debug_assert_ne!(wall, q, "a face's edge rides another plane");
                                BoundPiece::Segment { ends, wall }
                            }
                            Wall::Arc { cyl, ccw, plane } => {
                                // A plane face's arcs lie in its plane; another section is not
                                // one this reading has.
                                if plane != q {
                                    hold(Meet::Along(ends));
                                    continue;
                                }
                                let ends = if ccw { ends } else { [ends[1], ends[0]] };
                                BoundPiece::Arc { ends, cyl }
                            }
                            Wall::Ruling {
                                cyl, side, plane, ..
                            } => {
                                // A ruling's `side` is measured against `plane`, and this reading
                                // asks it against the face's own class.
                                if plane != q {
                                    hold(Meet::Along(ends));
                                    continue;
                                }
                                BoundPiece::Ruling { ends, cyl, side }
                            }
                        });
                    }
                    // 1. Where the ring passes a node more than once, each arrival leaves by the
                    //    way out the face's walk takes there.
                    if !repeats.is_empty() {
                        repeats.sort_unstable();
                        match r.edges(jd, cyls, Some(q)) {
                            Ok(edges) => {
                                for &u in &repeats {
                                    match leaves_as_the_face_walks(jd, cyls, q, &r.nodes, &edges, u)
                                    {
                                        Some(true) => {}
                                        Some(false) => found_at(Meet::Node(u)),
                                        None => hold(Meet::Node(u)),
                                    }
                                }
                            }
                            Err(_) => hold(Meet::Node(repeats[0])),
                        }
                    }
                }
                Bound::Circle { cyl } => pieces.push(BoundPiece::Circle { cyl: *cyl }),
                // A lateral's rim bounds the lateral, never a plane face.
                Bound::Rim { .. } => hold(Meet::Unplaced),
            }
        }
        boxes.clear();
        boxes.extend(pieces.iter().map(|p| piece_box(jd, cyls, q, p, pt)));
        let overlap = |(al, ah): ([f64; 3], [f64; 3]), (bl, bh): ([f64; 3], [f64; 3])| {
            (0..3).all(|k| al[k] <= bh[k] && bl[k] <= ah[k])
        };
        // 2 and 3, over the pairs whose boxes overlap: a node inside a piece lies in the box of every
        // piece it ends, so the pairs carry every node a piece can hold.
        for a in 0..pieces.len() {
            for b in a + 1..pieces.len() {
                if !overlap(boxes[a], boxes[b]) {
                    continue;
                }
                let (pa, pb) = (&pieces[a], &pieces[b]);
                if pa.same_set(pb) {
                    found_at(piece_witness(pa, pb));
                    continue;
                }
                let (near_a, near_b) = (
                    |n: NodeId| overlap(pt.reach(n), boxes[a]),
                    |n: NodeId| overlap(pt.reach(n), boxes[b]),
                );
                match pair_meets(jd, cyls, q, (pa, &near_a), (pb, &near_b)) {
                    Ok(()) => {}
                    Err(Some(m)) => found_at(m),
                    Err(None) => hold(piece_witness(pa, pb)),
                }
            }
        }
    }
    let witness = |m: Meet| match m {
        Meet::Node(n) => pt.witness_point(n),
        Meet::Along(e) => pt.witness_segment(e),
        Meet::Unplaced => None,
    };
    let raise = |reason: RejectReason, m: Meet| match witness(m) {
        Some(at) => crate::reject_at(reason, at),
        None => reject(reason),
    };
    if let Some(m) = found {
        return Err(raise(RejectReason::FaceBoundsMeet, m));
    }
    if let Some(m) = held {
        return Err(raise(RejectReason::WitnessNotRational, m));
    }
    Ok(())
}

/// **Do two pieces of one plane face meet anywhere but at a node they share?** — part 2 (an end of
/// one inside the other) and part 3 (their interiors) of [`face_bounds_meet`], `Err(Some(_))` with
/// where, `Err(None)` where it cannot be read. Two segments read each end's side of the other's line
/// once and answer both parts from it. Each piece comes with whether a node's proven reach meets its
/// box: a node outside it is not inside the piece, and is not asked.
fn pair_meets(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    (a, near_a): (&BoundPiece, &dyn Fn(NodeId) -> bool),
    (b, near_b): (&BoundPiece, &dyn Fn(NodeId) -> bool),
) -> Result<(), Option<Meet>> {
    if let (
        BoundPiece::Segment { ends: e1, wall: w1 },
        BoundPiece::Segment { ends: e2, wall: w2 },
    ) = (*a, *b)
    {
        let s = |n: NodeId, w: usize| side_on(jd, cyls, n, w);
        // One line, `q ∩ w`: an end of one on the other, other than as its end, is the overlap.
        let overlap = |w: usize| -> Result<(), Option<Meet>> {
            for (n, [u, v], near) in [
                (e2[0], e1, near_a),
                (e2[1], e1, near_a),
                (e1[0], e2, near_b),
                (e1[1], e2, near_b),
            ] {
                if near(n) && on_straight(jd, cyls, q, w, [u, v], n).ok_or(None)? {
                    return Err(Some(Meet::Node(n)));
                }
            }
            Ok(())
        };
        if w1 == w2 {
            return overlap(w1);
        }
        // Two lines through a node both pieces end at meet only there — unless they are one line,
        // which the far end of either says.
        if let Some(&shared) = e1.iter().find(|n| e2.contains(n)) {
            let far = |e: [NodeId; 2]| if e[0] == shared { e[1] } else { e[0] };
            return match (s(far(e2), w1), s(far(e1), w2)) {
                (Some(0), _) => overlap(w1),
                (_, Some(0)) => overlap(w2),
                (Some(_), _) | (_, Some(_)) => Ok(()),
                (None, None) => Err(None),
            };
        }
        let (b0, b1, a0, a1) = (s(e2[0], w1), s(e2[1], w1), s(e1[0], w2), s(e1[1], w2));
        // A piece wholly on one side of the other's line meets nothing of it.
        let off = |x: Option<i8>, y: Option<i8>| matches!((x, y), (Some(x), Some(y)) if x * y > 0);
        if off(b0, b1) || off(a0, a1) {
            return Ok(());
        }
        let (b0, b1, a0, a1) = (
            b0.ok_or(None)?,
            b1.ok_or(None)?,
            a0.ok_or(None)?,
            a1.ok_or(None)?,
        );
        // An end on the other's line, on its span other than as one of its ends.
        let inside = |n: NodeId,
                      side: i8,
                      [u, v]: [NodeId; 2],
                      w: usize,
                      near: &dyn Fn(NodeId) -> bool|
         -> Result<bool, Option<Meet>> {
            if side != 0 || !near(n) {
                return Ok(false);
            }
            on_straight(jd, cyls, q, w, [u, v], n).ok_or(None)
        };
        for (n, side, seg, w, near) in [
            (e2[0], b0, e1, w1, near_a),
            (e2[1], b1, e1, w1, near_a),
            (e1[0], a0, e2, w2, near_b),
            (e1[1], a1, e2, w2, near_b),
        ] {
            if inside(n, side, seg, w, near)? {
                return Err(Some(Meet::Node(n)));
            }
        }
        // Their interiors cross: every end strictly across the other's line.
        if b0 * b1 < 0 && a0 * a1 < 0 {
            return Err(Some(piece_witness(a, b)));
        }
        return Ok(());
    }
    // A segment and a ruling sharing an end meet only there where the segment's line does not run
    // along the axis: two lines of the face's plane that are not parallel meet once. Where it does
    // run along it, the segment lies on the ruling's own line, and the nodes answer the overlap
    // (measured over the workspace suite: 170 of 29,624 shared ends, none refused).
    if let (BoundPiece::Segment { ends: e, wall }, BoundPiece::Ruling { ends: r, cyl, .. })
    | (BoundPiece::Ruling { ends: r, cyl, .. }, BoundPiece::Segment { ends: e, wall }) = (*a, *b)
        && e.iter().any(|n| r.contains(n))
        && crosses_axis(jd, cyls, q, wall, cyl, &e) == Some(true)
    {
        return Ok(());
    }
    for (p, o, near) in [(a, b, near_a), (b, a, near_b)] {
        for &n in o.ends() {
            if p.ends().contains(&n) || !near(n) {
                continue;
            }
            match node_inside(jd, cyls, q, n, p) {
                Some(true) => return Err(Some(Meet::Node(n))),
                Some(false) => {}
                None => return Err(None),
            }
        }
    }
    match interiors_meet(jd, cyls, q, a, b) {
        Some(true) => Err(Some(piece_witness(a, b))),
        Some(false) => Ok(()),
        None => Err(None),
    }
}

/// **Which side of plane `w` node `n` lies on** — `0` where its name holds `w` (a point named by a
/// plane lies on it, whatever else can be read), [`node_side`] otherwise.
fn side_on(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
    w: usize,
) -> Option<i8> {
    if n.holds_plane(w) {
        return Some(0);
    }
    node_side(jd, cyls, n, w)
}

/// The witness for two pieces: the smaller of their names.
fn piece_witness(a: &BoundPiece, b: &BoundPiece) -> Meet {
    let of = |p: &BoundPiece| match p.ends() {
        [u, v] => Meet::Along([*u.min(v), *u.max(v)]),
        _ => Meet::Unplaced,
    };
    of(a).min(of(b))
}

/// **Does every arrival at `u` leave by the way out the face's walk takes there?** — part 1 of
/// [`face_bounds_meet`] at a node the ring `nodes` passes more than once (edge `i` runs `nodes[i] →
/// nodes[i + 1]`): of the ring's departures from `u`, [`combinatorics::way_out`] from each arrival
/// must be that visit's own. `None` where a direction or their order cannot be read — two
/// departures on one ray among them, an overlap part 2 proves.
fn leaves_as_the_face_walks(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    nodes: &[NodeId],
    edges: &[combinatorics::RingEdge],
    u: NodeId,
) -> Option<bool> {
    let k = nodes.len();
    let visits: Vec<usize> = (0..k).filter(|&i| nodes[i] == u).collect();
    let mut ways = Vec::with_capacity(visits.len());
    for &i in &visits {
        ways.push(combinatorics::dir_at(jd, cyls, q, &edges[i], u).ok()?);
    }
    for (w, &i) in visits.iter().enumerate() {
        let back = combinatorics::dir_at(jd, cyls, q, &edges[(i + k - 1) % k], u)
            .ok()?
            .reversed();
        if combinatorics::way_out(jd, q, &back, &ways).ok()? != w {
            return Some(false);
        }
    }
    Some(true)
}

/// **Is node `n` on the straight piece `[u, v]` of the line `q ∩ wall`, other than as one of its
/// ends?** — strictly between them, or at an end's point under another name (two names, one point:
/// a meet, as at an arc's end), by [`combinatorics::closed_contains`]. `None` where a point cannot
/// be placed on the line.
fn on_straight(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    wall: usize,
    [u, v]: [NodeId; 2],
    n: NodeId,
) -> Option<bool> {
    if n == u || n == v {
        return Some(false);
    }
    let pin = |m: NodeId| combinatorics::end_on_line(jd, q, wall, m);
    let at = combinatorics::locate(jd, cyls, q, wall, pin(n)?)?;
    let ends = [
        combinatorics::on_line(jd, cyls, q, wall, pin(u)?)?,
        combinatorics::on_line(jd, cyls, q, wall, pin(v)?)?,
    ];
    combinatorics::closed_contains(jd, q, wall, &at, &ends, &|m| m)
}

/// **Is node `n` inside piece `p`** — on a segment ([`on_straight`]) or a ruling ([`on_ruling_span`])
/// other than as one of its ends, on an arc's circle strictly inside the arc, anywhere on a whole
/// circle. `None` where it cannot be read. A node standing at a piece's end under another name is a
/// meet too: two names, one point.
fn node_inside(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    n: NodeId,
    p: &BoundPiece,
) -> Option<bool> {
    use combinatorics::ThetaPlace;
    let on_circle = |cyl: usize| -> Option<bool> {
        Some(matches!(
            combinatorics::named_point_in_disk(jd, cyls, n, &cyls.get(cyl)?.def)?,
            RingSide::OnBoundary
        ))
    };
    // A node of another cylinder's surface is a point of that cylinder's lateral, which never meets
    // this one's on a result ([`face_bounds_meet`]'s third part).
    if let (Some((_, c, _)), Some(k)) = (combinatorics::pierce_name(n), p.cyl())
        && c != k
    {
        return Some(false);
    }
    match *p {
        BoundPiece::Segment { ends: [u, v], wall } => match side_on(jd, cyls, n, wall)? {
            0 => on_straight(jd, cyls, q, wall, [u, v], n),
            _ => Some(false),
        },
        BoundPiece::Arc {
            ends: [lo, hi],
            cyl,
        } => {
            if !on_circle(cyl)? {
                return Some(false);
            }
            match combinatorics::theta_place(jd, cyl, &cyls[cyl].def, lo, hi, n)? {
                ThetaPlace::Inside | ThetaPlace::At(_) => Some(true),
                ThetaPlace::Outside => Some(false),
            }
        }
        BoundPiece::Circle { cyl } => on_circle(cyl),
        BoundPiece::Ruling {
            ends: [u, v],
            cyl,
            side,
        } => {
            if !on_circle(cyl)? {
                return Some(false);
            }
            let def = &cyls[cyl].def;
            if combinatorics::node_ruling_side(jd, def, q, n)? != side {
                return Some(false);
            }
            on_ruling_span(jd, cyls, def, n, [u, v])
        }
    }
}

/// **Is node `n` on the closed span of a ruling piece `[u, v]`** — its height along the axis between
/// the ends' or at one of them; `None` unread. Closed, because a node at an end's height on that
/// ruling under another name is the end's point — two names, one point, a meet as at an arc's end.
fn on_ruling_span(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    def: &nacre_topo::CylinderDef,
    n: NodeId,
    [u, v]: [NodeId; 2],
) -> Option<bool> {
    let h = |m: NodeId| combinatorics::axis_param_on_cylinder(jd, cyls, m, def);
    let (hn, hu, hv) = (h(n)?, h(u)?, h(v)?);
    let (a, b) = (hn.cmp_value(&hu)?, hn.cmp_value(&hv)?);
    Some(a == core::cmp::Ordering::Equal || a != b)
}

/// **Does the line `q ∩ wall` cross cylinder `cyl`'s axis direction** rather than run along it —
/// `None` unread. An end of the segment named as this cylinder's pierce of that very line proves it
/// (the line meets the cylinder at isolated points, so it is no ruling); otherwise the wall's world
/// normal against the axis decides.
fn crosses_axis(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    wall: usize,
    cyl: usize,
    ends: &[NodeId],
) -> Option<bool> {
    let named = ends.iter().any(|&n| {
        combinatorics::pierce_name(n)
            .is_some_and(|(planes, c, _)| c == cyl && planes.contains(&wall) && planes.contains(&q))
    });
    if named {
        return Some(true);
    }
    let w = jd.world_rat(wall)?;
    let axis = cyls.get(cyl)?.def.dir();
    Some(nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], &axis) != nacre_exact::Orient::Zero)
}

/// The proven box of a piece: its ends' reach for a straight piece, the arc's own box, the whole
/// circle's; unbounded where nothing is proven, which keeps every pair.
fn piece_box(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    p: &BoundPiece,
    pt: &NodeReach,
) -> crate::planes::CacheReach {
    let whole = |c: &crate::planes::WorkingCyl| {
        combinatorics::arc_box_rat(jd, q, &c.def, crate::planes::Theta::Whole)
    };
    match *p {
        BoundPiece::Segment { ends, .. } | BoundPiece::Ruling { ends, .. } => pt.edge_reach(ends),
        BoundPiece::Arc { ends, cyl } => crate::planes::rat_box_reach(
            cyls.get(cyl)
                .and_then(|c| combinatorics::arc_edge_box_rat(jd, cyls, q, &c.def, ends)),
        ),
        BoundPiece::Circle { cyl } => crate::planes::rat_box_reach(cyls.get(cyl).and_then(whole)),
    }
}

/// **Do the interiors of two pieces of one plane face meet?** — after every end, T and overlap has
/// been asked of the nodes ([`node_inside`]); two segments [`pair_meets`] answers itself. `None`
/// where it cannot be read.
fn interiors_meet(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    q: usize,
    a: &BoundPiece,
    b: &BoundPiece,
) -> Option<bool> {
    use nacre_judge::cylinder::LineMeet;
    match (*a, *b) {
        (BoundPiece::Segment { ends, wall }, BoundPiece::Arc { ends: arc, cyl })
        | (BoundPiece::Arc { ends: arc, cyl }, BoundPiece::Segment { ends, wall }) => {
            let def = &cyls.get(cyl)?.def;
            let [lo, hi] = arc;
            let on_line = |n: NodeId| side_on(jd, cyls, n, wall).map(|x| x == 0);
            let (l0, l1) = (on_line(lo)?, on_line(hi)?);
            // A line meets a circle twice at most: both ends on it, it meets the arc nowhere else.
            if lo != hi && l0 && l1 {
                return Some(false);
            }
            let on_end = if l0 {
                Some(lo)
            } else if l1 {
                Some(hi)
            } else {
                None
            };
            let cross = combinatorics::arc_crossings(jd, def, wall, q, cyl, true, arc, on_end)?;
            let touch = combinatorics::arc_touch(jd, cyls, wall, q, cyl, true, arc)?;
            let mut meet = false;
            for m in cross.into_iter().chain(touch) {
                meet |= on_straight(jd, cyls, q, wall, ends, m)?;
            }
            Some(meet)
        }
        (BoundPiece::Segment { ends, wall }, BoundPiece::Circle { cyl })
        | (BoundPiece::Circle { cyl }, BoundPiece::Segment { ends, wall }) => {
            let c = cyls.get(cyl)?;
            let def = &c.def;
            // The arrangement's cheap rejection, on its premise (a segment ⊥ the axis).
            if !c.sections_oblique(q)
                && let (Some(p0), Some(p1)) = (
                    combinatorics::node_coords_rat(jd, ends[0]),
                    combinatorics::node_coords_rat(jd, ends[1]),
                )
                && (p0 == p1
                    || !nacre_exact::segment_meets_cylinder(
                        &p0,
                        &p1,
                        &def.origin(),
                        &def.dir(),
                        def.r2(),
                    ))
            {
                return Some(false);
            }
            let roots: &[combinatorics::QuadRoot] = match jd.line_meet(wall, q, def, false) {
                LineMeet::Pair(_) => &[combinatorics::QuadRoot::Lo, combinatorics::QuadRoot::Hi],
                LineMeet::Tangent(..) => &[combinatorics::QuadRoot::Double],
                LineMeet::Miss | LineMeet::ParallelPlanes => &[],
                _ => return None,
            };
            let mut meet = false;
            for &root in roots {
                meet |= on_straight(jd, cyls, q, wall, ends, NodeId::pierce(wall, q, cyl, root))?;
            }
            Some(meet)
        }
        (BoundPiece::Segment { ends, wall }, BoundPiece::Ruling { ends: r, cyl, side })
        | (BoundPiece::Ruling { ends: r, cyl, side }, BoundPiece::Segment { ends, wall }) => {
            let def = &cyls.get(cyl)?.def;
            // A wall along the axis is no crossing, rightly: the line is parallel to the ruling —
            // another line, or the same one, whose overlap the nodes answer.
            let m = match combinatorics::crossing_on_ruling(jd, def, wall, q, cyl, side) {
                Ok(m) => m,
                Err(combinatorics::RulingNameFail::NoCrossing) => return Some(false),
                Err(_) => return None,
            };
            if !on_straight(jd, cyls, q, wall, ends, m)? {
                return Some(false);
            }
            on_ruling_span(jd, cyls, def, m, r)
        }
        // Two curves of one cylinder are one circle or parallel rulings, whose meets are their
        // ends; two curves of two cylinders never meet on a result ([`face_bounds_meet`]'s third
        // part says why).
        _ => Some(false),
    }
}
