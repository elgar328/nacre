use super::*;
/// **A ring node's coordinate, in whichever world names it** — the key
/// [`loop_winding`]'s lexicographic scan orders by. An arc's own minimum is ordered in the same
/// currency: it is the pierce point of the cap plane, a plane through the axis, and the cylinder
/// ([`arc_extremum_winding`]), though no ring node stands there.
///
/// ★★★ **The two arms are not two widths of one thing, they are two *kinds*.** A three-plane node
/// is the rational meet of three planes; a pierce node's coordinate is `a + b√c` and no rational
/// vessel holds it. That is why this is an enum and not a `[Rat; 3]` with a decline: the second
/// arm is not a precision failure to be lifted, it is a different number.
enum CoordKey {
    /// The three classes, handed to `Judge::cmp_coord` — which keeps its toleranced ladder and its
    /// escalation, so the existing population's answers are bit-identical to before.
    Three([usize; 3]),
    /// The point a plane pair cuts out of a cylinder: `base + s·dir` with `s = a + b√c`.
    /// Boxed: this arm is an order of magnitude wider than a name, and a ring of names is the
    /// common case.
    Pierce(Box<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)>),
}

/// Build one ring node's key.
///
/// ★ **The cylinder comes from the class table.** The name says which cylinder
/// (`NodeId::Pierce` carries the class), and the table's def is the statement every carrier's
/// def is a clone of. Searching the ring's own arc/ruling carriers instead — "a pierce
/// node is an arc endpoint" — is what a **chord** refutes: a cell bounded by a cap's chord
/// alone has pierce corners and only plane carriers, and the search refuses an honestly-named
/// point (`PierceVertexUnnamed` on the straddling flush corpus, measured).
fn coord_key(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    ring: &[RingEdge],
    i: usize,
) -> Result<CoordKey, BoolError> {
    let node = ring[i].node;
    let NodeKind::Pierce { cyl, .. } = node.kind() else {
        // A `match` and not a fallback: a third variant must light this up rather than fall in
        // here (`let`-`else` is what hid a new variant once already).
        return match node.kind() {
            NodeKind::ThreePlane(t) => Ok(CoordKey::Three(t)),
            NodeKind::Pierce { .. } => unreachable!("the let-else above took every pierce node"),
        };
    };
    let def = &cyls
        .get(cyl)
        .ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?
        .def;
    let (line, s) =
        pierce_meet(jd, cyl, def, node).ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
    Ok(CoordKey::Pierce(Box::new((line, s))))
}

/// Order two ring nodes along one world axis — `+1` when `a`'s coordinate is the larger.
///
/// The mixed pair is `nacre_exact::quad::cmp_coord_meet_branch`, which is exact and **total**:
/// both coordinates lift to one first-storey sign. It wants the three-plane side as a `MeetPoint`,
/// and this crate builds one directly — `MeetPoint::Narrow` is a public variant, so no door has to
/// be opened in `nacre-exact` for it.
fn cmp_key(
    jd: &Judge<'_, WorkingPlane>,
    a: &CoordKey,
    b: &CoordKey,
    axis: usize,
) -> Result<i8, BoolError> {
    use nacre_exact::{Orient, quad};
    let sign = |o: Orient| match o {
        Orient::Positive => 1i8,
        Orient::Negative => -1,
        Orient::Zero => 0,
    };
    // ★ Through [`node_coords_rat`], not a second copy of the three-plane solve — the rational
    // meet is stated once. The round-trip through the name is free: the solve is symmetric in its
    // three planes, so canonical order changes nothing.
    let meet = |t: [usize; 3]| {
        node_coords_rat(jd, NodeId::three_planes(Canon3::three(t)))
            .map(nacre_exact::MeetPoint::Narrow)
            .ok_or_else(|| reject(RejectReason::WitnessNotRational))
    };
    Ok(match (a, b) {
        (CoordKey::Three(x), CoordKey::Three(y)) => jd.cmp_coord(*x, *y, axis),
        (CoordKey::Three(x), CoordKey::Pierce(b)) => {
            sign(quad::cmp_coord_meet_branch(&meet(*x)?, &b.0, &b.1, axis))
        }
        (CoordKey::Pierce(b), CoordKey::Three(y)) => {
            -sign(quad::cmp_coord_meet_branch(&meet(*y)?, &b.0, &b.1, axis))
        }
        (CoordKey::Pierce(a), CoordKey::Pierce(b)) => {
            sign(quad::cmp_coord_branch((&a.0, &a.1), (&b.0, &b.1), axis))
        }
    })
}

/// **The pierce nodes lying strictly between two ring nodes, in ring order** — the exact half of
/// the split-twin subdivision (`assembly::name_result_vertices`' opening pass). The caller has
/// already matched the candidates' plane pair to the edge's `{own, wall}`, so by name every
/// candidate lies on the edge's own carrier line and single-axis order *is* order along it.
///
/// ★ **The axis is "wherever the endpoints differ", not the line's direction.** Two distinct
/// points of one line differ on some axis, the line is strictly monotone on that axis, and
/// betweenness is direction-blind — so no direction vector is read at all, and the choice is
/// deterministic (first differing axis). This is [`cmp_key`]'s vocabulary end to end; nothing new
/// is exact here.
///
/// `None` when an order cannot be formed (a coordinate outside the rational vessel, an
/// escalation, a class with no coefficients). The caller leaves such an edge **unsplit** — the
/// far-plane road starves there and the walls-fallback net answers — so the conservative arm
/// degrades to the unsplit state, never to something new.
pub(crate) fn pierce_between(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    a: NodeId,
    b: NodeId,
    candidates: &[NodeId],
) -> Option<Vec<NodeId>> {
    let key = |n: NodeId| -> Option<CoordKey> {
        match n.kind() {
            NodeKind::ThreePlane(t) => Some(CoordKey::Three(t)),
            NodeKind::Pierce { cyl, .. } => {
                let (line, s) = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?;
                Some(CoordKey::Pierce(Box::new((line, s))))
            }
        }
    };
    let (ka, kb) = (key(a)?, key(b)?);
    let axis = (0..3).find(|&ax| matches!(cmp_key(jd, &ka, &kb, ax), Ok(s) if s != 0))?;
    // `+1` = "the first is larger" (`cmp_key`), so `dir` is the a→b slope's sign on this axis and
    // "strictly between" is both hops running the same way.
    let dir = cmp_key(jd, &ka, &kb, axis).ok()?;
    let mut mid: Vec<(NodeId, CoordKey)> = Vec::new();
    for &c in candidates {
        if c == a || c == b {
            continue;
        }
        let kc = key(c)?;
        if cmp_key(jd, &ka, &kc, axis).ok()? == dir && cmp_key(jd, &kc, &kb, axis).ok()? == dir {
            mid.push((c, kc));
        }
    }
    // Ring order: nearer to `a` first — along `dir`, the larger-toward-`a` side leads. A pair
    // this cannot strictly order (an escalation, or two distinct nodes at one coordinate — which
    // on a shared line would be two names for one point) makes the whole edge unsplittable.
    let mut sortable = true;
    mid.sort_by(|(_, x), (_, y)| match cmp_key(jd, x, y, axis) {
        Ok(s) if s == dir => std::cmp::Ordering::Less,
        Ok(s) if s == -dir => std::cmp::Ordering::Greater,
        _ => {
            sortable = false;
            std::cmp::Ordering::Equal
        }
    });
    if !sortable {
        return None;
    }
    Some(mid.into_iter().map(|(n, _)| n).collect())
}

/// **The ring's own lexicographic minimum when it lies inside an arc**, and the winding read
/// there — `None` when every arc's minimum is at a node (then [`loop_winding`]'s `lo` is the
/// ring's minimum and its turn is the winding, as its doc argues).
///
/// ★★★★★ **Why this exists: `loop_winding`'s premise is about the *node set*, and a ring is not
/// its nodes.** Its doc reads *"the lexicographically smallest node … is an extreme point of the
/// node set, which is planar, so it is a vertex of the ring's hull"* — true for a polygon, and
/// false the moment an edge is an **arc**, because the arc can bulge past every node. Then the
/// turn at `lo` is read at a point the region does not support, and the sign comes back
/// **confident**: a 270° sector's cap, whose smallest node is the reflex centre, winds its two
/// orbits the wrong way round, and the walk's `−1`-cell count cannot see it (both orbits flip).
///
/// ★ **The answer was named three milestones ago** and is not a wider walk:
/// [`RejectReason::CurvedStraightRun`](crate::RejectReason::CurvedStraightRun)'s doc says *"read
/// the winding at the extremum of the **region**, which may lie in an arc's interior"*. This is
/// that reading, and the winding there is [`smooth_extremum_winding`]'s product — the ring is
/// smooth at an arc's interior point, so no turn is needed.
///
/// **The circle's minimum is a pierce point.** Let `ê_a` be the first world axis the circle
/// **spans** (`ê₀`, unless the axis *is* `ê₀` — then every point shares `x` and the minimum is taken
/// in `y`). The point of least coordinate `a` lies on the plane `H` through the axis with normal
/// `n_h = m × ê_a` (it is `c − s·p` for `p = ê_a − (m_a/|m|²)·m`, which lies in the span of `m` and
/// `ê_a`), so it is one of the two points the cap plane, `H` and the cylinder share — the
/// [`plane_plane_cylinder`](nacre_exact::quad::plane_plane_cylinder) meet every pierce node is
/// solved by, and a [`CoordKey::Pierce`] that [`cmp_key`] orders against the ring's nodes and
/// against another arc's minimum. The point is irrational whenever the axis tilts toward `ê_a` or
/// the radius is a surd (`c_x − 3√41/25` on a Pythagorean frame; `−√2` for a sector of `r² = 2` on
/// the world frame); a road that answered only for a rational `c − r·ê_a` left both to the node's
/// turn, and on a sector whose centre is the smallest node that turn was the wrong one.
///
/// ★ **What it cannot decide it refuses by name.** A minimum this cannot place or order is a ring
/// whose winding has no proof; reading the node's turn there would be the confident wrong sign
/// this function exists to catch. The meet's width limit is `WitnessNotRational`, as
/// [`coord_key`] names the same meet failing for a node.
fn arc_extremum_winding(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
    keys: &[CoordKey],
    lo: usize,
) -> Result<Option<i8>, BoolError> {
    use nacre_exact::quad::{CylinderMeet, plane_plane_cylinder};
    use nacre_exact::{Orient, Rat};
    let width = || reject(RejectReason::WitnessNotRational);
    let key = |i: usize| &keys[i];
    // `x` against `y`, axis by axis — `-1` when `x` is lexicographically the smaller, `0` when
    // they are the same point.
    let lex = |x: &CoordKey, y: &CoordKey| -> Result<i8, BoolError> {
        for ax in 0..3 {
            let o = cmp_key(jd, x, y, ax)?;
            if o != 0 {
                return Ok(o);
            }
        }
        Ok(0)
    };
    let n = ring.len();
    let zero = Rat::from_int(0);
    let mut best: Option<(CoordKey, i8)> = None;
    for (i, e) in ring.iter().enumerate() {
        let Carrier::Arc(ac) = &e.carrier else {
            continue;
        };
        // ★ Read first: `arc_at` is where an arc's start is required to be a pierce name
        // (`RingNaming` otherwise), and the halves below must not read an end that is not.
        let EdgeDir::Arc(ad) = dir_at(jd, cyls, p, e, e.node)? else {
            unreachable!("an arc carrier's direction is an arc")
        };
        let (o, m, r2) = (ac.def.origin(), ac.def.dir(), ac.def.r2());
        // The first world axis the circle spans.
        let a = usize::from(m[1] == zero && m[2] == zero);
        let mut e_a = [zero; 3];
        e_a[a] = Rat::from_int(1);
        // `H`: through the axis, so through the centre — `n_h · c = n_h · o` because `c − o ∥ m`.
        let n_h = nacre_exact::cross3_rat(&m, &e_a).ok_or_else(width)?;
        let h_plane = nacre_exact::dot3_rat(&n_h, &o)
            .and_then(|d| zero.checked_sub(d))
            .map(|d| [n_h[0], n_h[1], n_h[2], d])
            .ok_or_else(width)?;
        let cap = class_coeffs_rat(jd, p).ok_or_else(width)?;
        let ext = match plane_plane_cylinder(&cap, &h_plane, &o, &m, r2).ok_or_else(width)? {
            // `s` ascends along `line.dir()`, whose `a` component is `−k(|m|² − m_a²) ≠ 0`, so the
            // two roots differ in `a` and the smaller one is read off that one sign.
            CylinderMeet::Pair { line, s: [s0, s1] } => {
                let s = if line.dir()[a] > zero { s0 } else { s1 };
                CoordKey::Pierce(Box::new((line, s)))
            }
            // `H` holds the axis, which the cap plane crosses at the centre — inside the cylinder,
            // so the line through it meets the surface twice; and the cap plane is ⊥ `m` while `H`
            // contains `m`, so the two are neither one plane nor parallel.
            CylinderMeet::CoincidentPlanes
            | CylinderMeet::ParallelPlanes
            | CylinderMeet::OnRuling(_)
            | CylinderMeet::AxisParallelMiss(_)
            | CylinderMeet::Miss(_)
            | CylinderMeet::Tangent { .. } => {
                unreachable!("a line through a circle's centre in its plane meets it twice")
            }
        };
        // Is the circle's minimum lexicographically below `lo`? ★ Three answers, not two: every
        // axis equal means it **is** `lo`, and the premise holds.
        if lex(&ext, key(lo))? != -1 {
            continue;
        }
        // ★ **And is it in the arc's INTERIOR?** It cannot be an endpoint: `lo` is the smallest
        // ring **node** and this point is smaller still, so it is no node of this ring at all.
        // (That is also why a boss seated on a wall needs nothing special here — its circle is cut
        // by a **diameter**, so this point *is* a node, and the comparison above answered "not
        // below".)
        let (ka, kb) = (key(i), key((i + 1) % n));
        // The walk below runs counter-clockwise **about the axis** — the arc's own sense, so the
        // ends are taken in that order whichever way the ring traverses it.
        let (ka, kb) = if ac.ccw { (ka, kb) } else { (kb, ka) };
        // ★★ **The halves are the circle's own.** `H` splits the circle: at θ = 0 (`c + s·p`, the
        // maximum of coordinate `a`) counter-clockwise travel runs along `m × p = m × ê_a = +n_h`,
        // so the `+n_h` half is θ ∈ (0°, 180°) — where coordinate `a` falls — and the minimum
        // θ = 180° is where the walk **arrives from** the `+n_h` half and **leaves into** the
        // `−n_h` one. Nothing here reads a world picture, so no «as seen in a plane» correction is
        // needed whichever way the axis points. ★ An end *on* `H` is θ = 0° — θ = 180° is the
        // minimum itself, taken out above — and it belongs to the half the walk is in beside it:
        // a start leaves θ = 0° into the upper half, an end arrives at it from the lower.
        let half = |k: &CoordKey, is_start: bool| -> Result<bool, BoolError> {
            let o = match k {
                CoordKey::Three(t) => {
                    let q = node_coords_rat(jd, NodeId::three_planes(Canon3::three(*t)))
                        .ok_or_else(width)?;
                    let v = nacre_exact::dot3_rat(&n_h, &q)
                        .and_then(|v| v.checked_add(h_plane[3]))
                        .ok_or_else(width)?;
                    match v.cmp(&zero) {
                        std::cmp::Ordering::Greater => Orient::Positive,
                        std::cmp::Ordering::Less => Orient::Negative,
                        std::cmp::Ordering::Equal => Orient::Zero,
                    }
                }
                CoordKey::Pierce(b) => nacre_exact::quad::plane_side(&h_plane, &b.0, &b.1),
            };
            Ok(match o {
                Orient::Positive => true,
                Orient::Negative => false,
                Orient::Zero => is_start,
            })
        };
        let (ha, hb) = (half(ka, true)?, half(kb, false)?);
        // Walking CCW from the start, θ = 180° is reached iff the walk leaves the upper half, or
        // wraps the whole way round inside one half — and θ's order inside a half is read off
        // coordinate `a`: falling in the upper half, rising in the lower.
        let x_cmp = cmp_key(jd, ka, kb, a)?;
        let hit = match (ha, hb) {
            (true, false) => true,
            (false, true) => false,
            (true, true) => x_cmp <= 0,
            (false, false) => x_cmp >= 0,
        };
        if !hit {
            continue;
        }
        // The winding read **there**: the ring is smooth at an arc's interior point.
        let w = smooth_extremum_winding(jd, p, &ad);
        // ★ **The minimum, not the first.** Two arcs of one ring can each dip below `lo` only if
        // they ride different circles; the ring is supported at the lower of the two, and reading
        // the other would ask about a point the region is not extreme at.
        let take = match &best {
            None => true,
            Some((b, _)) => lex(&ext, b)? == -1,
        };
        if take {
            best = Some((ext, w));
        }
    }
    Ok(best.map(|(_, w)| w))
}

/// An ordered ring's winding about the face's outward normal: `-1` clockwise — the material
/// is *outside* the ring, so it bounds a hole — and `+1` counter-clockwise, an island.
///
/// The turn at a convex-hull vertex is the winding, and the lexicographically smallest node
/// is one: it is an extreme point of the node set, which is planar, so it is a vertex of the
/// ring's hull. Finding it is the **only** thing here that needs two implicit points in one
/// decision, and [`Judge::cmp_coord`] is that predicate.
///
/// A shortcut dies here, and is recorded so it is not walked twice: a *supporting edge* —
/// one whose plane `Q_j` has every other node on one side — would give a hull vertex from
/// a one-implicit side test (`Judge::orient3d`) alone. But a simple polygon need not have an edge
/// on its hull (fold each side of a pentagon slightly inward), so no such edge is guaranteed.
/// A hull *vertex* always exists.
///
/// ★★ **Where arcs touch this**: the turn is read at **one** node, so a curved edge
/// needs no angle sum — only its tangent's direction at that node. The other two sites that read
/// the direction's *representation* each got their own answer: the walk-back below asks
/// [`continuation`], which has a curved arm ("are the tangents parallel" is "same circle, same
/// travel"), and the lexicographic minimum above runs on [`CoordKey`], which holds a pierce node's
/// `a + b√c` coordinate beside a name and compares across the two through the quad tower.
///
/// A ring may be *non-simple* — visiting one node twice — and still be a legitimate face: the
/// unbounded contour of two cells that meet at a single point pinches through that point, tracing
/// a figure-8. The winding is read from the turn at the lexicographically smallest node, and a
/// coincidence elsewhere in the ring does not affect that turn, so a repeated node is not by itself
/// an error. (Rejecting on the first coincidence with the running minimum would make the verdict
/// depend on the ring's arbitrary start index — one operand order rejecting a pinch the other
/// accepts.) Only a pinch *at* the extreme node itself leaves the turn ambiguous; that is
/// `CoincidentNodes`, decided by exact equality rather than by a tolerance.
pub(crate) fn loop_winding(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
) -> Result<i8, BoolError> {
    // ★ **The floor reads the carriers.** Two straight edges between two points are one edge traced
    // twice; two arcs are a lens and an arc with its chord is a circular segment. See the same rule
    // at the walk's orbit length (`arrangement::walk_cells`).
    if ring.len()
        < if ring.iter().any(|e| matches!(e.carrier, Carrier::Arc(_))) {
            2
        } else {
            3
        }
    {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // ★ **The comparator, in one place, reading the identity directly.** This is the second of the
    // two sites that deliberately do not go through [`three_plane_name`]: `Judge::cmp_coord` speaks
    // three plane indices (it lives in `nacre-judge`, below this crate, so the name cannot travel
    // there), and a pierce point's coordinate is `a + b√c` with its own total comparators
    // (`nacre_exact::quad::cmp_coord_meet_branch` and `cmp_coord_branch`). **This dispatch is
    // where that decision belongs** — the door's single answer is the wrong one here, and
    // `Judge::cmp_coord`'s four-rung ladder speaks neither `MeetLine` nor `QuadVal`.
    //
    // ★★ **The keys are materialized, and that is a change of shape, not just of type.** The old
    // spelling read `[usize; 3]` out of the name on every comparison — free. A pierce key is a
    // *solve* ([`pierce_meet`] re-derives the point from the name), and the scan below asks for
    // each node's key on the order of six times, so re-solving per read would multiply the exact
    // work by that. One pass, one `Vec`.
    let keys: Vec<CoordKey> = (0..ring.len())
        .map(|i| coord_key(jd, cyls, ring, i))
        .collect::<Result<_, BoolError>>()?;
    let key = |i: usize| &keys[i];
    // Lexicographically smallest node — a hull vertex, hence a valid turn site. A coincidence with
    // the running minimum just means "not strictly smaller", so keep it; do not reject.
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let mut order = 0i8;
        for axis in 0..3 {
            order = cmp_key(jd, key(i), key(lo), axis)?;
            if order != 0 {
                break;
            }
        }
        if order == -1 {
            lo = i;
        }
    }
    // ★★★ **The scan above is only a minimum if the relation is an order, and that is an
    // assumption about the predicates, not about this loop.** It composes per-axis comparisons
    // lexicographically; a per-axis answer that is not a fact about the geometry — one plane
    // described two ways, say — makes the composition intransitive, and then a forward scan can
    // stop at a node with something smaller behind it. That node is not extreme, the turn read
    // there is not the winding, and the wrong sign comes back **confident**: the engine noticed
    // only two layers later, as "no outer contour", and named the symptom.
    //
    // ★ This is the postcondition the algorithm actually needs — cheaper than asking whether the
    // relation is transitive (`O(n)` against `O(n³)`, and rings here reach 95 nodes) and closer to
    // the point. **It holds however the predicates behave; it is the net under them.**
    // ★ A declining comparison makes no claim, so it cannot witness a violation either: `0` is
    // "says nothing" here, not "equal".
    let lex = |i: usize, j: usize| -> i8 {
        (0..3)
            .map(|axis| cmp_key(jd, key(i), key(j), axis).unwrap_or(0))
            .find(|&c| c != 0)
            .unwrap_or(0)
    };
    debug_assert!(
        !(0..ring.len()).any(|i| i != lo && lex(i, lo) == -1),
        "the lexicographic scan did not find a minimum — the comparison is not an order here"
    );
    // The turn is read at `lo`; if that exact point recurs the corner is a pinch and its turn is
    // ambiguous — honest-reject rather than guess.
    // ★ Stops at the first pinch. Running on would ask more comparisons, and one of those could
    // *decline* — turning a `CoincidentNodes` that is already decided into a width reject.
    for i in 0..ring.len() {
        let mut same = i != lo;
        for axis in 0..3 {
            same = same && cmp_key(jd, key(i), key(lo), axis)? == 0;
        }
        if same {
            return Err(reject(RejectReason::CoincidentNodes));
        }
    }
    // ★★★★★ **The ring's minimum may not be a node at all.** The scan above found the smallest
    // **node**; an arc can bulge past it, and then `lo` is not a hull vertex and the turn read
    // there is not the winding. [`arc_extremum_winding`] answers where that happens, from the arc
    // itself — the reading `CurvedStraightRun`'s doc named.
    if let Some(w) = arc_extremum_winding(jd, cyls, p, ring, &keys, lo)? {
        return Ok(w);
    }
    // **A ring node need not be a corner.** The arrangement names a point wherever another feature
    // crosses an edge, and `loop_triples` keeps such a vertex even when the loop runs straight
    // through it — so one edge of the polygon can arrive as several collinear ring edges. Reading
    // the turn at `lo` against its immediate predecessor then asks about two halves of one
    // straight edge, which has no turn to give.
    //
    // The turn to read is the one between the directions the loop **actually** arrives and leaves
    // on: walk back past the edges the loop runs straight through. `lo` stays a hull vertex — the
    // stretch lies on one **line** through it, so the region is still on one side of that line.
    // ★ That argument is the straight one, and the guard below is where it stops: two arcs of one
    // circle are also "straight through", and a *circle* through `lo` does not put the region on
    // one side of anything.
    //
    // **Only while the stretch keeps going the same way.** A collinear edge traversed the *other*
    // way means the ring doubles back along the line it came in on — an antenna, whose tip has no
    // turn and whose neighbours' turn belongs to a different vertex. Skipping past that would
    // read a turn from somewhere else and call it this vertex's: a wrong winding, silently.
    //
    // ★★★ **The comparison is between neighbours, at the node they share** — not between the
    // candidate and `ring[lo]`, which is the same answer for straight edges (parallel and
    // same-sense are both transitive along a chain of shared points) and **meaningless** once an
    // edge is curved: a far arc's tangent is not `lo`'s tangent, so comparing them asks about two
    // different places. Transitivity is also why neighbour-only does not weaken "the whole
    // stretch runs one way".
    //
    // ★★★★ **And the value carried out is the one already read at a shared node** — never an
    // edge's direction at its own far start. For a straight edge the two agree; on a **diameter**
    // chord they are exactly opposite, and that is what a boss straddling a plate edge measured
    // (both half-disks are `+1`; the half whose arc *arrives* came back `−1`). Where the stretch is curved the equality that licenses stepping at
    // all fails, and the walk refuses by name rather than reading a winding from the wrong place.
    let n = ring.len();
    let leaving = dir_at(jd, cyls, p, &ring[lo], ring[lo].node)?;
    let mut back = (lo + n - 1) % n;
    // ★ The direction the loop arrives on, carried out of the walk — the edge that ends it is the
    // one the turn is read against, and its direction is already in hand.
    let arriving = loop {
        let ahead = (back + 1) % n;
        let shared = ring[ahead].node;
        let earlier = dir_at(jd, cyls, p, &ring[back], shared)?;
        let later = dir_at(jd, cyls, p, &ring[ahead], shared)?;
        match continuation(jd, p, &earlier, &later)? {
            // ★★★★ **`earlier`, and not `ring[back]`'s direction at its own start.** The two are
            // the same value for a straight edge — a line's tangent does not change along it — and
            // that equality is what let the older spelling stand. On an arc they differ by the
            // whole turn of the arc, and on a *diameter* chord they are exactly opposite: measured,
            // the two half-disks of a boss straddling a plate edge came back `+1` and `−1` where
            // both are `+1`, because the half whose arc **arrives** read its tangent at the far
            // end. This one is read at the node the loop actually passes through.
            Continuation::Turns => break earlier,
            // ★ **A cusp at the extremum**: a line and an arc tangent at `lo`, the ring arriving
            // on one and leaving back along the other (a keyhole's cap less its bore: the box's
            // wall and the rim meet at the corner with no angle between them). The region is the
            // sliver between the two, on the arc's **convex** side — outside its circle — so the
            // winding is the opposite of a smooth join's ([`smooth_extremum_winding`]). Deeper in
            // the walk a doubling back is an antenna and stays refused.
            Continuation::DoublesBack
                if ahead == lo
                    && let (EdgeDir::Line { .. }, EdgeDir::Arc(a))
                    | (EdgeDir::Arc(a), EdgeDir::Line { .. }) = (&earlier, &later) =>
            {
                return Ok(-smooth_extremum_winding(jd, p, a));
            }
            Continuation::DoublesBack => return Err(reject(RejectReason::StraightAngle)),
            // ★★★ **The step is licensed by the stretch being a *line*.** What the walk carries out
            // is `ring[back]`'s direction at its own start, and that equals `lo`'s arriving
            // direction only because a line's tangent is the same everywhere on it. Two arcs of one
            // circle are tangent-continuous, so `continuation` answers `Straight` for them too —
            // and stepping there would read the winding from a different point of the ring.
            // ★★★★★ **At `lo` itself a smooth boundary still states a winding — by curvature.**
            // The step past a straight run is licensed by the stretch being a *line*, and two arcs
            // of one circle are tangent-continuous without being one; stepping there would read the
            // winding from a different point of the ring. But at `lo` there is nothing to step
            // past: the loop **is** smooth at the extreme node, and a smooth extremum's winding is
            // the arc's own rotation. See [`smooth_extremum_winding`].
            // Not curved here: an ordinary straight run through `lo`, walked back as before.
            Continuation::Straight
                if ahead == lo
                    && let (EdgeDir::Arc(e), EdgeDir::Arc(l)) = (&earlier, &later) =>
            {
                if e.cyl != l.cyl || e.ccw != l.ccw || e.axis_up != l.axis_up {
                    return Err(reject(RejectReason::CurvedStraightRun));
                }
                return Ok(smooth_extremum_winding(jd, p, l));
            }
            // ★ A smooth **line–arc** join at the extremum (a fillet's corner is
            // the rounded rectangle's extreme node): the ring turns there only at second order,
            // and the arc's bending is that turn.
            Continuation::Straight
                if ahead == lo
                    && let (EdgeDir::Line { .. }, EdgeDir::Arc(a))
                    | (EdgeDir::Arc(a), EdgeDir::Line { .. }) = (&earlier, &later) =>
            {
                return Ok(smooth_extremum_winding(jd, p, a));
            }
            Continuation::Straight
                if matches!(earlier, EdgeDir::Arc(_)) || matches!(later, EdgeDir::Arc(_)) =>
            {
                return Err(reject(RejectReason::CurvedStraightRun));
            }
            Continuation::Straight => {}
        }
        back = (back + n - 1) % n;
        if back == lo {
            // Every edge of the ring lies on one line: it bounds nothing. ★ This is the *loop's*
            // termination, not one of `continuation`'s answers — it is about having walked the
            // whole ring, not about what any one edge does.
            return Err(reject(RejectReason::DegenerateRing));
        }
    };
    turn_between(jd, p, &arriving, &leaving)
}

/// **The winding of a ring that runs *smooth* through its extreme node** — curvature, not a turn.
///
/// ★★★★★ **This is not [`turn`]'s question, and that is why it is not [`turn`]'s arm.** `turn`
/// answers "how much does the direction rotate at this node", and for two arcs of one circle the
/// honest answer is `0`: they are tangent-continuous, nothing rotates *at* the node. What
/// [`loop_winding`] needs there is a different fact — which way the boundary **curves** — and it is
/// available only because the node is the ring's lexicographic minimum.
///
/// **Why the minimum makes it answerable.** `lo` is a hull vertex: the whole ring lies on one side
/// of a supporting line through it. The boundary there is an arc, so the arc curves off that line
/// into the side the ring is on — the region is locally convex at `lo`, and a locally convex point's
/// turn carries the ring's orientation. For an arc that "turn" is spread along the arc rather than
/// concentrated at a vertex, but its **sign** is the arc's own rotation, which is exactly the
/// winding.
///
/// **The sign, in three factors.** Travel rotates about `s·m`, `s = +1` when `ccw`. The turn's
/// reference is the face's **outward** normal, `n_out = frame_sign · n_P`, and `n_P · m > 0` is
/// `axis_up`. So
///
/// ```text
///   (s·m) · n_out = s · frame_sign · (n_P · m)
///     ⇒  winding = ccw · axis_up · frame_sign
/// ```
///
/// — no coordinate, no predicate, three signs the directions already carry.
///
/// The two arcs come from [`Continuation::Straight`], which for arcs means *one cylinder and the
/// same travel sense*, so both agree on every factor; they are re-checked here rather than assumed,
/// because this function's answer is a sign and a wrong one is silent.
///
/// ☑ **Which factors the population locks.** Negating the product, dropping `ccw`, and dropping
/// `axis_up` each turn 60–75 of the crate's tests red, refused as
/// [`RejectReason::RingOrientation`] and [`RejectReason::LabelConflict`] (the plane cells and the
/// chart's labels no longer agreeing), so the lock names them. Dropping `frame_sign` changes nothing: it is `+1` on every class that reaches
/// this rule today, which is the same shape [`turn`]'s own note records for its factor — the
/// difference being that `turn`'s corpus does reach `Reversed` faces and this rule's does not yet.
/// ☑ Under the motion group it is still not caught — the commuting
/// oracle's 396 always-on cells stay green with the factor dropped, so a ring whose lexicographic
/// minimum is a smooth arc node on a `frame_sign = −1` class is a population no fixture has yet
/// (the arc extremum rung reads the smooth minimum *inside* an arc, [`arc_extremum_winding`],
/// which the oracle does exercise).
fn smooth_extremum_winding(jd: &Judge<'_, WorkingPlane>, p: usize, arc: &ArcDir) -> i8 {
    let sign = |b: bool| if b { 1i8 } else { -1 };
    sign(arc.ccw) * sign(arc.axis_up) * jd.planes[p].frame_sign
}

/// `sign((n_P × n_Q) · N_R)`, where `N_R` is the right-hand normal of `R`'s witness triangle.
///
/// [`Judge::plane_pair_dir_sign`] gives the sign against `R`'s *stored* normal, judged on the
/// truth. That normal is parallel to `N_R` but opposes it on a `Reversed` face, which the face's
/// `frame_sign` corrects.
///
/// The correction *is* the face's stated flag — since the stored-orientation
/// cutover, `frame_sign` is `Forward`/`Reversed` as a sign, and
/// "`Reversed` ⇔ `n_out = −plane.normal()`" holds by construction rather than by
/// hope. What keeps it honest is the winding: `collect_planes` winds the witness triangle by the
/// truth's sense and the face's orientation, and `validate` pins the loop itself as
/// `FaceMisoriented`.
pub(crate) fn dir_sign(jd: &Judge<'_, WorkingPlane>, p: usize, q: usize, r: usize) -> i8 {
    let planes = jd.planes;
    jd.plane_pair_dir_sign(p, q, r) * planes[r].frame_sign
}
