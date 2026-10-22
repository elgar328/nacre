use super::*;
/// **A ring node's coordinate, in whichever world names it** — the key
/// [`loop_winding`]'s lexicographic scan orders by. An arc's own minimum is ordered in the same
/// currency: it is the pierce point of the cap plane, a plane through the axis, and the cylinder
/// ([`arc_extremum_winding`]), though no ring node stands there.
///
/// **The two arms are not two widths of one thing, they are two *kinds*.** A three-plane node
/// is the rational meet of three planes; a pierce node's coordinate is `a + b√c` and no rational
/// vessel holds it. That is why this is an enum and not a `[Rat; 3]` with a decline: the second
/// arm is not a precision failure to be lifted, it is a different number.
enum CoordKey {
    /// The three classes, handed to `Judge::cmp_coord` — which keeps its toleranced ladder and its
    /// escalation.
    Three([usize; 3]),
    /// The point a plane pair cuts out of a cylinder: `base + s·dir` with `s = a + b√c`.
    /// Boxed: this arm is an order of magnitude wider than a name, and a ring of names is the
    /// common case.
    Pierce(Box<PierceMeet>),
}

/// Build one ring node's key.
///
/// **The cylinder comes from the class table.** The name says which cylinder
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
    let NodeKind::Pierce { planes, cyl, .. } = node.kind() else {
        // A `match` and not a fallback: a third variant must light this up rather than fall into
        // this `else` unnamed.
        return match node.kind() {
            NodeKind::ThreePlane(t) => Ok(CoordKey::Three(t)),
            NodeKind::Pierce { .. } => unreachable!("the let-else above took every pierce node"),
        };
    };
    let def = &cyls
        .get(cyl)
        .ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?
        .def;
    // A meet that cannot be formed is the width limit for a named pair, and for one with a plane of
    // no world statement a box the precision cannot form or a tangency — that refusal is the
    // gate's ([`undecided_reason`]).
    let meet =
        pierce_meet(jd, cyl, def, node).ok_or_else(|| reject(undecided_reason(cyls, &planes)))?;
    Ok(CoordKey::Pierce(Box::new(meet)))
}

/// **What ordering two ring nodes on one world axis needs besides the two keys** — the judge, and
/// the cylinder table a refusal reads its name from ([`undecided_reason`]).
struct KeyCmp<'a> {
    jd: &'a Judge<'a, WorkingPlane>,
    cyls: &'a [crate::planes::WorkingCyl],
}

/// The world axis a class's plane holds constant — the one nonzero component of its world normal.
/// `None` for a plane no world axis is normal to, or one with no world statement.
fn flat_axis(jd: &Judge<'_, WorkingPlane>, p: usize) -> Option<usize> {
    let zeros = jd.planes[p].world.as_ref()?.name.normal_zeros();
    let mut nonzero = (0..3).filter(|&k| !zeros[k]);
    let k = nonzero.next()?;
    nonzero.next().is_none().then_some(k)
}

/// The exact pair the `a + b√c` comparators read — a meet without one has no road in [`KeyCmp`].
/// Never handed an interval box: [`KeyCmp::cmp`] sends any pair with one to [`KeyCmp::bounded`]
/// first, so the refusal here is the width of a named meet.
fn exact(
    m: &PierceMeet,
) -> Result<(&nacre_exact::quad::MeetLine, &nacre_exact::quad::QuadVal), BoolError> {
    m.exact()
        .ok_or_else(|| reject(RejectReason::WitnessNotRational))
}

impl KeyCmp<'_> {
    /// Order two ring nodes along one world axis — `+1` when `a`'s coordinate is the larger.
    ///
    /// The mixed pair is `nacre_exact::quad::cmp_coord_meet_branch`, which is exact and **total**:
    /// both coordinates lift to one first-storey sign. It wants the three-plane side as a
    /// `MeetPoint`, and this crate builds one directly — `MeetPoint::Narrow` is a public variant, so
    /// no door has to be opened in `nacre-exact` for it. A three-plane node with no rational meet
    /// (a class with no world statement among its three) is ordered by intervals instead
    /// ([`Self::unstated`]).
    fn cmp(&self, a: &CoordKey, b: &CoordKey, axis: usize) -> Result<i8, BoolError> {
        use nacre_exact::{Orient, quad};
        let jd = self.jd;
        let bounded = |k: &CoordKey| matches!(k, CoordKey::Pierce(m) if m.exact().is_none());
        if bounded(a) || bounded(b) {
            return self.bounded(a, b, axis);
        }
        let sign = |o: Orient| match o {
            Orient::Positive => 1i8,
            Orient::Negative => -1,
            Orient::Zero => 0,
        };
        // Through [`node_coords_rat`], not a second copy of the three-plane solve — the rational
        // meet is stated once. The round-trip through the name is free: the solve is symmetric in
        // its three planes, so canonical order changes nothing.
        let meet = |t: [usize; 3]| {
            node_coords_rat(jd, NodeId::three_planes(Canon3::three(t)))
                .map(nacre_exact::MeetPoint::Narrow)
        };
        let order = match (a, b) {
            (CoordKey::Three(x), CoordKey::Three(y)) => jd.cmp_coord(*x, *y, axis),
            (CoordKey::Three(x), CoordKey::Pierce(m)) => match meet(*x) {
                Some(p) => {
                    let (line, s) = exact(m)?;
                    sign(quad::cmp_coord_meet_branch(&p, line, s, axis))
                }
                None => self.unstated(*x, m, axis)?,
            },
            (CoordKey::Pierce(m), CoordKey::Three(y)) => match meet(*y) {
                Some(p) => {
                    let (line, s) = exact(m)?;
                    -sign(quad::cmp_coord_meet_branch(&p, line, s, axis))
                }
                None => -self.unstated(*y, m, axis)?,
            },
            (CoordKey::Pierce(a), CoordKey::Pierce(b)) => {
                sign(quad::cmp_coord_branch(exact(a)?, exact(b)?, axis))
            }
        };
        Ok(order)
    }

    /// **A pair one of which is a pierce point with no exact meet** (a plane of its name has no
    /// world statement), in three rules and this order:
    /// * **stated**: both points state the coordinate exactly — a plane of the name (three-plane or
    ///   pierce) the axis is normal to, or an exact meet line running flat along the axis — and
    ///   they are compared exactly;
    /// * **one ruling**: two such pierce points on one nameless plane that **sections** the cylinder
    ///   ([`crate::planes::WorkingCyl::sections`] — along the axis, within the radius), on one ruling
    ///   of it, share every coordinate the ruling holds constant (the axis's component is `0`) —
    ///   equal by name, which no interval could prove. The section is what makes the rule true: an
    ///   oblique plane meets the cylinder in an ellipse, which `ruling_side`'s sign only halves, and
    ///   two points of one half are two points;
    /// * otherwise their interval boxes, and a pair the boxes cannot separate is refused, never
    ///   read as a tie ([`undecided_reason`] names it).
    fn bounded(&self, a: &CoordKey, b: &CoordKey, axis: usize) -> Result<i8, BoolError> {
        let jd = self.jd;
        let zero = nacre_exact::Rat::from_int(0);
        let flat_value = |c: usize| -> Option<nacre_exact::Rat> {
            if flat_axis(jd, c) != Some(axis) {
                return None;
            }
            let w = class_coeffs_rat(jd, c)?;
            zero.checked_sub(w[3])?
                .checked_mul(nacre_exact::Rat::new(w[axis].denom(), w[axis].numer())?)
        };
        let stated = |k: &CoordKey| -> Option<nacre_exact::Rat> {
            match k {
                CoordKey::Three(t) => t.iter().find_map(|&c| flat_value(c)),
                CoordKey::Pierce(m) => match m.exact() {
                    Some((line, _)) => (line.dir()[axis] == zero).then(|| line.base()[axis]),
                    None => m.name_planes()?.iter().find_map(|&c| flat_value(c)),
                },
            }
        };
        if let (Some(x), Some(y)) = (stated(a), stated(b)) {
            return Ok(i8::from(x > y) - i8::from(x < y));
        }
        if let (CoordKey::Pierce(x), CoordKey::Pierce(y)) = (a, b)
            && let (Some((px, cx)), Some((py, cy))) =
                (x.name_planes_and_cyl(), y.name_planes_and_cyl())
            && cx == cy
            && let Some(cyl) = self.cyls.get(cx)
            && cyl.def.dir()[axis] == zero
            && let Some(&w) = px.iter().find(|&&c| py.contains(&c) && cyl.sections(c))
            && let (Some(sx), Some(sy)) = (
                x.ruling_side(jd, w, &cyl.def),
                y.ruling_side(jd, w, &cyl.def),
            )
            && sx == sy
        {
            return Ok(0);
        }
        let classes: Vec<usize> = [a, b]
            .iter()
            .flat_map(|k| match k {
                CoordKey::Three(t) => t.to_vec(),
                CoordKey::Pierce(m) => m.name_planes().map_or(Vec::new(), |p| p.to_vec()),
            })
            .collect();
        let refuse = || reject(undecided_reason(self.cyls, &classes));
        let point = |k: &CoordKey| match k {
            CoordKey::Three(t) => three_point_hp(jd, *t),
            CoordKey::Pierce(m) => m.point_hp(jd),
        };
        let (x, y) = (point(a).ok_or_else(refuse)?, point(b).ok_or_else(refuse)?);
        match x[axis].sub(&y[axis], jd.standard.prec).sign() {
            Some(true) => Ok(1),
            Some(false) => Ok(-1),
            None => Err(refuse()),
        }
    }

    /// A three-plane node no rational meet states, against a pierce point: exactly where both
    /// state the coordinate ([`Self::both_stated`]), else the two as intervals ([`three_point_hp`],
    /// [`PierceMeet::point_hp`]) — two distinct names, so a comparison the intervals cannot separate is
    /// refused, never read as a tie ([`undecided_reason`] names it).
    fn unstated(&self, t: [usize; 3], m: &PierceMeet, axis: usize) -> Result<i8, BoolError> {
        if let Some(order) = m
            .exact()
            .and_then(|(line, _)| self.both_stated(t, line, axis))
        {
            return Ok(order);
        }
        let refuse = || reject(undecided_reason(self.cyls, &t));
        let x = three_point_hp(self.jd, t).ok_or_else(refuse)?;
        let y = m.point_hp(self.jd).ok_or_else(refuse)?;
        match x[axis].sub(&y[axis], self.jd.standard.prec).sign() {
            Some(true) => Ok(1),
            Some(false) => Ok(-1),
            None => Err(refuse()),
        }
    }

    /// **A coordinate both points state exactly, compared exactly.** The three-plane node states
    /// its coordinate on `axis` when one of its classes is a plane that axis is normal to
    /// (`n·x + d = 0`, so `x = −d/n`); the pierce point states it when its line runs flat along
    /// that axis (`x = base`). This is where **equality comes from**: two such points can share
    /// the coordinate exactly — every node of a ring on the wall `x = 9`, or a box's corner on that
    /// wall against a bore's pierce point on it while the ring lies on a cap — and an interval can
    /// prove two coordinates apart, never equal. `None` where either side does not state it, or on
    /// overflow.
    fn both_stated(
        &self,
        t: [usize; 3],
        line: &nacre_exact::quad::MeetLine,
        axis: usize,
    ) -> Option<i8> {
        let zero = nacre_exact::Rat::from_int(0);
        if line.dir()[axis] != zero {
            return None;
        }
        let w = t.iter().find_map(|&c| {
            (flat_axis(self.jd, c) == Some(axis))
                .then(|| class_coeffs_rat(self.jd, c))
                .flatten()
        })?;
        // `−d/n − base` has the sign of `−(n·base + d)·n`.
        let v = w[axis].checked_mul(line.base()[axis])?.checked_add(w[3])?;
        let sign = |r: nacre_exact::Rat| i8::from(r > zero) - i8::from(r < zero);
        Some(-sign(v) * sign(w[axis]))
    }
}

/// **The pierce nodes lying strictly between two ring nodes, in ring order** — the exact half of
/// the split-twin subdivision (`assembly::name_result_vertices`' opening pass). The caller has
/// already matched the candidates' plane pair to the edge's `{own, wall}`, so by name every
/// candidate lies on the edge's own carrier line and single-axis order *is* order along it.
///
/// **The axis is "wherever the endpoints differ", not the line's direction.** Two distinct
/// points of one line differ on some axis, the line is strictly monotone on that axis, and
/// betweenness is direction-blind — so no direction vector is read at all, and the choice is
/// deterministic (first differing axis). This is [`KeyCmp::cmp`]'s vocabulary end to end; nothing new
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
            NodeKind::Pierce { cyl, .. } => Some(CoordKey::Pierce(Box::new(pierce_meet(
                jd,
                cyl,
                &cyls.get(cyl)?.def,
                n,
            )?))),
        }
    };
    let kc = KeyCmp { jd, cyls };
    let (ka, kb) = (key(a)?, key(b)?);
    let axis = (0..3).find(|&ax| matches!(kc.cmp(&ka, &kb, ax), Ok(s) if s != 0))?;
    // `+1` = "the first is larger" (`KeyCmp::cmp`), so `dir` is the a→b slope's sign on this axis and
    // "strictly between" is both hops running the same way.
    let dir = kc.cmp(&ka, &kb, axis).ok()?;
    let mut mid: Vec<(NodeId, CoordKey)> = Vec::new();
    for &c in candidates {
        if c == a || c == b {
            continue;
        }
        let kx = key(c)?;
        if kc.cmp(&ka, &kx, axis).ok()? == dir && kc.cmp(&kx, &kb, axis).ok()? == dir {
            mid.push((c, kx));
        }
    }
    // Ring order: nearer to `a` first — along `dir`, the larger-toward-`a` side leads. A pair
    // this cannot strictly order (an escalation, or two distinct nodes at one coordinate — which
    // on a shared line would be two names for one point) makes the whole edge unsplittable.
    let mut order: Vec<usize> = (0..mid.len()).collect();
    try_sort_by(&mut order, |&i, &j| {
        kc.cmp(&mid[i].1, &mid[j].1, axis).map(|s| {
            if s == dir {
                std::cmp::Ordering::Less
            } else if s == -dir {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        })
    })
    .ok()?;
    // Two nodes tied on the axis sort next to each other, and a tie is unsplittable.
    for w in order.windows(2) {
        if kc.cmp(&mid[w[0]].1, &mid[w[1]].1, axis).ok()? == 0 {
            return None;
        }
    }
    Some(order.into_iter().map(|i| mid[i].0).collect())
}

/// **The ring's own lexicographic minimum when it lies inside an arc**, and the winding read
/// there — `None` when every arc's minimum is at a node (then [`loop_winding`]'s `lo` is the
/// ring's minimum and its turn is the winding, as its doc argues).
///
/// **Why this exists: `loop_winding`'s premise is about the *node set*, and a ring is not
/// its nodes.** Its doc reads *"the lexicographically smallest node … is an extreme point of the
/// node set, which is planar, so it is a vertex of the ring's hull"* — true for a polygon, and
/// false the moment an edge is an **arc**, because the arc can bulge past every node. Then the
/// turn at `lo` is read at a point the region does not support, and the sign comes back
/// **confident**: a 270° sector's cap, whose smallest node is the reflex centre, winds its two
/// orbits the wrong way round, and the walk's `−1`-cell count cannot see it (both orbits flip).
///
/// **This is where the region's extremum is found when it is not a node**, and the reading there
/// is [`loop_winding`]'s own: the ring is smooth at an arc's interior point, so its two departures
/// are exactly opposite and the outside half turn is the side the arc curves away from —
/// [`smooth_extremum_winding`]'s product.
///
/// **The circle's minimum is a pierce point.** Let `ê_a` be the first world axis the circle
/// **spans** (`ê₀`, unless the axis *is* `ê₀` — then every point shares `x` and the minimum is taken
/// in `y`). The point of least coordinate `a` lies on the plane `H` through the axis with normal
/// `n_h = m × ê_a` (it is `c − s·p` for `p = ê_a − (m_a/|m|²)·m`, which lies in the span of `m` and
/// `ê_a`), so it is one of the two points the cap plane, `H` and the cylinder share — the
/// [`plane_plane_cylinder`](nacre_exact::quad::plane_plane_cylinder) meet every pierce node is
/// solved by, and a [`CoordKey::Pierce`] that [`KeyCmp::cmp`] orders against the ring's nodes and
/// against another arc's minimum. The point is irrational whenever the axis tilts toward `ê_a` or
/// the radius is a surd (`c_x − 3√41/25` on a Pythagorean frame; `−√2` for a sector of `r² = 2` on
/// the world frame); a road that answers only for a rational `c − r·ê_a` leaves both to the node's
/// turn, which on a sector whose centre is the smallest node is the wrong one.
///
/// **What it cannot decide it refuses by name.** A minimum this cannot place or order is a ring
/// whose winding has no proof; reading the node's turn there would be the confident wrong sign
/// this function exists to catch. The meet's width limit is `WitnessNotRational`, as
/// [`coord_key`] names the same meet failing for a node. A class with no world statement never
/// asks: its own arrangement carries no arc (a class that sections a cylinder cuts it in rulings,
/// and a cylinder beside or apart from it leaves none on its plane), so `p` and every arc end here
/// are stated.
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
    let kc = KeyCmp { jd, cyls };
    let key = |i: usize| &keys[i];
    // `x` against `y`, axis by axis — `-1` when `x` is lexicographically the smaller, `0` when
    // they are the same point.
    let lex = |x: &CoordKey, y: &CoordKey| -> Result<i8, BoolError> {
        for ax in 0..3 {
            let o = kc.cmp(x, y, ax)?;
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
        // Read first: `arc_at` is where an arc's start is required to be a pierce name
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
                CoordKey::Pierce(Box::new(PierceMeet::Exact(line, s)))
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
        // Is the circle's minimum lexicographically below `lo`? Three answers, not two: every
        // axis equal means it **is** `lo`, and the premise holds.
        if lex(&ext, key(lo))? != -1 {
            continue;
        }
        // **And is it in the arc's INTERIOR?** It cannot be an endpoint: `lo` is the smallest
        // ring **node** and this point is smaller still, so it is no node of this ring at all.
        // (That is also why a boss seated on a wall needs nothing special here — its circle is cut
        // by a **diameter**, so this point *is* a node, and the comparison above answered "not
        // below".)
        let (ka, kb) = (key(i), key((i + 1) % n));
        // The walk below runs counter-clockwise **about the axis** — the arc's own sense, so the
        // ends are taken in that order whichever way the ring traverses it.
        let (ka, kb) = if ac.ccw { (ka, kb) } else { (kb, ka) };
        // **The halves are the circle's own.** `H` splits the circle: at θ = 0 (`c + s·p`, the
        // maximum of coordinate `a`) counter-clockwise travel runs along `m × p = m × ê_a = +n_h`,
        // so the `+n_h` half is θ ∈ (0°, 180°) — where coordinate `a` falls — and the minimum
        // θ = 180° is where the walk **arrives from** the `+n_h` half and **leaves into** the
        // `−n_h` one. Nothing here reads a world picture, so no «as seen in a plane» correction is
        // needed whichever way the axis points. An end *on* `H` is θ = 0° — θ = 180° is the
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
                CoordKey::Pierce(b) => {
                    let (line, s) = b.exact().ok_or_else(width)?;
                    nacre_exact::quad::plane_side(&h_plane, line, s)
                }
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
        let x_cmp = kc.cmp(ka, kb, a)?;
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
        // **The minimum, not the first.** Two arcs of one ring can each dip below `lo` only if
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

/// **The fewest edges a ring can have and still bound a cell** — the one rule the walk's orbit
/// length (`arrangement::walk_cells`), [`loop_winding`] and the coplanar merge's cycles
/// (`assembly::coplanar`) read, so they cannot part.
///
/// Two straight edges between two points are one edge traced twice, so a straight ring needs
/// three. Two arcs between two points are a lens, and an arc with its chord a circular segment —
/// honest cells at two. One arc whose two ends are the **same node** is a circle cut at exactly
/// that point — a polygon's vertex lying on it — and is a cell boundary on its own: the closed
/// edge with one vertex that OCCT (`BOPDS_DS::InitPaveBlocks`), Parasolid and ACIS all keep as a
/// self-loop, and Euler counts it as one vertex and one edge. A tangency never arrives here as
/// that ring: `split_circles` mints no vertex for a `Double` root, so the circle it grazes stays
/// whole.
pub(crate) fn ring_floor(any_arc: bool, self_loop_arc: bool) -> usize {
    if self_loop_arc {
        1
    } else if any_arc {
        2
    } else {
        3
    }
}

/// An ordered ring's winding about the face's outward normal: `-1` clockwise — the material
/// is *outside* the ring, so it bounds a hole — and `+1` counter-clockwise, an island.
///
/// **One rule, read at the ring's extreme point.** Let `P` be the lexicographically smallest point
/// of the ring — the smallest node `lo`, or a point inside an arc that bulges past every node
/// ([`arc_extremum_winding`]). The whole ring lies on one side of a line through `P`, so every
/// direction the ring leaves `P` along lies in one half-plane. Gather, for each passage of the ring
/// through `P`, the direction it leaves on and the reverse of the one it arrives on; order them
/// counter-clockwise ([`angular_order`]); `c1` is the first after the gap that faces outside. The
/// ring winds `+1` if `c1` leaves `P` and `−1` if it arrives — the winding of the sector beside
/// `c1`, fixed by the sector outside the hull (winding `0`) whatever the number of passages and
/// however they pair.
///
/// - A corner is two departures. A pinch — the ring passing `P` twice, as the contour of two cells
///   touching at a point does, or an outer boundary and a hole meeting there — is four or more,
///   and `c1` belongs to the outermost passage.
/// - An outgoing and an incoming departure along one ray (an antenna, its tip or its base at `P`)
///   add two equal and opposite contributions and cancel; if every departure cancels, the ring has
///   no side there ([`RejectReason::StraightAngle`]).
/// - Ties are the comparator's ([`departure_turn`]): a line and an arc leaving tangent the same way
///   (a cusp) are ordered by the arc's bending. Two departures exactly opposite — the ring smooth
///   at `P` — leave two half-turn gaps, and the outside one is the side the boundary curves away
///   from ([`smooth_extremum_winding`]).
///
/// This is the reading exact-arithmetic arrangement libraries use for curved boundaries (an extreme
/// vertex and the curves incident to it, ties broken by curvature); a signed area would need the
/// arcs' transcendental segments. Finding `P` is the **only** thing here that needs two implicit
/// points in one decision, and [`KeyCmp`] is that predicate. It is one rule rather than an arm per
/// local shape (design 「가지 말 것」: casework per failing configuration).
///
/// A shortcut dies here, and is recorded so it is not walked twice: a *supporting edge* —
/// one whose plane `Q_j` has every other node on one side — would give a hull vertex from
/// a one-implicit side test (`Judge::orient3d`) alone. But a simple polygon need not have an edge
/// on its hull (fold each side of a pentagon slightly inward), so no such edge is guaranteed.
/// A hull *vertex* always exists.
///
/// **Where arcs touch this**: the departures are read at **one** point, so a curved edge needs
/// no angle sum — only its tangent there. The lexicographic minimum runs on [`CoordKey`], which
/// holds a pierce node's `a + b√c` coordinate beside a name and compares across the two through
/// the quad tower.
pub(crate) fn loop_winding(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
) -> Result<i8, BoolError> {
    let self_loop =
        ring.len() == 1 && matches!(ring[0].carrier, Carrier::Arc(_)) && ring[0].node == ring[0].to;
    let any_arc = ring.iter().any(|e| matches!(e.carrier, Carrier::Arc(_)));
    if ring.len() < ring_floor(any_arc, self_loop) {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // **The comparator, in one place, reading the identity directly.** This is the second of the
    // two sites that deliberately do not go through [`three_plane_name`]: `Judge::cmp_coord` speaks
    // three plane indices (it lives in `nacre-judge`, below this crate, so the name cannot travel
    // there), and a pierce point's coordinate is `a + b√c` with its own total comparators
    // (`nacre_exact::quad::cmp_coord_meet_branch` and `cmp_coord_branch`). **This dispatch is
    // where that decision belongs** — the door's single answer is the wrong one here, and
    // `Judge::cmp_coord`'s ladder speaks neither `MeetLine` nor `QuadVal`.
    //
    // **The keys are materialized, once.** A three-plane key is free to read out of the name, but a
    // pierce key is a *solve* ([`pierce_meet`] re-derives the point from the name), and the scan
    // below asks for each node's key on the order of six times, so re-solving per read would
    // multiply the exact work by that. One pass, one `Vec`.
    let keys: Vec<CoordKey> = (0..ring.len())
        .map(|i| coord_key(jd, cyls, ring, i))
        .collect::<Result<_, BoolError>>()?;
    let key = |i: usize| &keys[i];
    let kc = KeyCmp { jd, cyls };
    // Lexicographically smallest node — a hull vertex, hence a valid turn site. A coincidence with
    // the running minimum just means "not strictly smaller", so keep it; do not reject.
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let mut order = 0i8;
        for axis in 0..3 {
            order = kc.cmp(key(i), key(lo), axis)?;
            if order != 0 {
                break;
            }
        }
        if order == -1 {
            lo = i;
        }
    }
    // **The scan above is only a minimum if the relation is an order, and that is an
    // assumption about the predicates, not about this loop.** It composes per-axis comparisons
    // lexicographically; a per-axis answer that is not a fact about the geometry — one plane
    // described two ways, say — makes the composition intransitive, and then a forward scan can
    // stop at a node with something smaller behind it. That node is not extreme, the turn read
    // there is not the winding, and the wrong sign comes back **confident**: the engine noticed
    // only two layers later, as "no outer contour", and named the symptom.
    //
    // This is the postcondition the algorithm actually needs — cheaper than asking whether the
    // relation is transitive (`O(n)` against `O(n³)`, and rings here reach 95 nodes) and closer to
    // the point. **It holds however the predicates behave; it is the net under them.**
    // A declining comparison makes no claim, so it cannot witness a violation either: `0` is
    // "says nothing" here, not "equal".
    let lex = |i: usize, j: usize| -> i8 {
        (0..3)
            .map(|axis| kc.cmp(key(i), key(j), axis).unwrap_or(0))
            .find(|&c| c != 0)
            .unwrap_or(0)
    };
    debug_assert!(
        !(0..ring.len()).any(|i| i != lo && lex(i, lo) == -1),
        "the lexicographic scan did not find a minimum — the comparison is not an order here"
    );
    // The region's extremum may lie inside an arc rather than at `lo`; the ring is smooth there
    // and the arc's rotation is the winding.
    if let Some(w) = arc_extremum_winding(jd, cyls, p, ring, &keys, lo)? {
        return Ok(w);
    }
    // `P` is `lo`. Every passage of the ring through it gives two departures: the direction it
    // leaves on, and the reverse of the one it arrives on. Every node is compared — a pinch is
    // one more passage, not a refusal.
    let n = ring.len();
    let mut found: Vec<(EdgeDir, bool)> = Vec::new();
    for i in 0..n {
        let mut at_p = true;
        for axis in 0..3 {
            at_p = at_p && kc.cmp(key(i), key(lo), axis)? == 0;
        }
        if !at_p {
            continue;
        }
        let node = ring[i].node;
        let arriving = dir_at(jd, cyls, p, &ring[(i + n - 1) % n], node)?;
        found.push((arriving.reversed(), false));
        found.push((dir_at(jd, cyls, p, &ring[i], node)?, true));
    }
    // An outgoing and an incoming departure that are **one ray** — one carrier, one way — add two
    // equal and opposite contributions and cancel: an antenna, whether its tip is `P` (one passage
    // in and back out) or its base (out on one passage, back on another). A line and an arc leaving
    // tangent (a cusp) are two curves, not one, and stay: the comparator orders them by the arc's
    // bending. Read as travel, one ray out and the same ray in are one carrier run both ways.
    let mut gone = vec![false; found.len()];
    for i in 0..found.len() {
        if gone[i] || !found[i].1 {
            continue;
        }
        for j in 0..found.len() {
            if gone[j] || found[j].1 {
                continue;
            }
            let back = found[j].0.reversed();
            if one_carrier_pole(jd, p, &found[i].0, &back)? {
                gone[i] = true;
                gone[j] = true;
                break;
            }
        }
    }
    let (dirs, leaves): (Vec<EdgeDir>, Vec<bool>) = found
        .into_iter()
        .zip(gone)
        .filter(|(_, g)| !g)
        .map(|(d, _)| d)
        .unzip();
    // Every departure at `P` cancelled: the ring has no side there.
    if dirs.is_empty() {
        return Err(reject(RejectReason::StraightAngle));
    }
    let order = angular_order(jd, p, &dirs)?;
    let m = order.len();
    // The gap that faces outside: the one counter-clockwise step between neighbours that turns
    // back — more than a half turn — and failing that, the half turn between the two departures on
    // the half-plane's edge.
    let mut first = None;
    for k in 0..m {
        let (a, b) = (order[k], order[(k + 1) % m]);
        if departure_turn(jd, p, &dirs[a], &dirs[b])? < 0 {
            first = Some(b);
            break;
        }
    }
    let c1 = match first {
        Some(b) => b,
        // Two departures exactly opposite: the ring is smooth at `P`, both gaps are half turns,
        // and the outside one is the side the boundary curves away from — the arc's own rotation,
        // read on the arc as travelled.
        None if m == 2 => {
            let arc = (0..2).find_map(|k| {
                let travelled = if leaves[k] {
                    dirs[k].clone()
                } else {
                    dirs[k].reversed()
                };
                match travelled {
                    EdgeDir::Arc(a) => Some(a),
                    _ => None,
                }
            });
            // Two straight departures exactly opposite would be a straight run through `P`,
            // which an extreme node cannot have.
            return match arc {
                Some(a) => Ok(smooth_extremum_winding(jd, p, &a)),
                None => Err(reject(RejectReason::StraightAngle)),
            };
        }
        None => {
            let mut edge = None;
            for k in 0..m {
                let (a, b) = (order[k], order[(k + 1) % m]);
                if departure_pole(jd, p, &dirs[a], &dirs[b])? {
                    edge = Some(b);
                }
            }
            edge.ok_or_else(|| reject(RejectReason::UnorderedEdges))?
        }
    };
    Ok(if leaves[c1] { 1 } else { -1 })
}

/// **The winding of a ring that runs *smooth* through its extreme node** — curvature, not a turn.
///
/// **This is not [`turn`]'s question, and that is why it is not [`turn`]'s arm.** `turn`
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
/// Its arc is the one the ring travels at its smooth extreme point: inside an arc
/// ([`arc_extremum_winding`]), or at a node where the two departures are exactly opposite — one
/// circle travelled one way through the node ([`antiparallel`] of the departure and the reversed
/// arrival), or a line running on into an arc tangent ([`tangent_pole`]).
///
/// **Which factors the population locks.** Negating the product, dropping `ccw`, and dropping
/// `axis_up` are each refused across the crate as [`RejectReason::RingOrientation`] and
/// [`RejectReason::LabelConflict`] (the plane cells and the chart's labels no longer agreeing), so
/// the lock names them. Dropping `frame_sign` is caught by the commuting oracle, the disk-side
/// census and the self-touch probes; dropped only where the minimum is a node (the arc extremum
/// rung, [`arc_extremum_winding`], keeping it), by the user's folds and fully rounded plates
/// (`sketch_models`, `user_models`) and the self-touch chord scoreboard.
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
/// The correction *is* the face's stated flag — `frame_sign` is `Forward`/`Reversed` as a sign, and
/// "`Reversed` ⇔ `n_out = −plane.normal()`" holds by construction rather than by
/// hope. What keeps it honest is the winding: `collect_planes` winds the witness triangle by the
/// truth's sense and the face's orientation, and `validate` pins the loop itself as
/// `FaceMisoriented`.
pub(crate) fn dir_sign(jd: &Judge<'_, WorkingPlane>, p: usize, q: usize, r: usize) -> i8 {
    let planes = jd.planes;
    jd.plane_pair_dir_sign(p, q, r) * planes[r].frame_sign
}
