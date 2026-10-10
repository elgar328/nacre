//! Mixed rings — straight steps and circular arcs — and the exact predicates the sketch policy
//! asks of them: does a point lie inside, do two rings meet, does a ring meet itself.
//!
//! [`crate::intersect`] answers those questions for polygons over `Rat`. An arc brings one square
//! root: a line meets a circle where a quadratic with rational coefficients vanishes, and two
//! circles meet where their radical line (rational) meets one of them, so every meeting point
//! here is `a + b√c` with **one radical per question** — exactly what
//! [`nacre_exact::quad::QuadVal`] carries and signs exactly. Nothing is realized in f64, and
//! nothing is snapped: an endpoint either is a point or is not.
//!
//! **The half-open ray rule is the one [`crate::intersect::ray_straddle`] states**, read on a
//! curve. A point on the ray is "not above"; the parity toggles at every transition between
//! *strictly above* and *not above* along the ring. For a straight step that is decided by its two
//! ends; on an arc a transition can also happen at an interior crossing of the ray's line, a
//! horizontal tangency is no transition, and an end on the ray toggles iff the arc leaves it
//! upward (start) or arrives at it from above (end) — which is the same sentence the straight
//! rule speaks about its ends. The kernel's arrangement applies this rule to its own arc-bearing
//! rings over named nodes (`nacre-ops`, the mixed-ring parity); this module is that rule over
//! rational 2-D data, spelled here because the two cannot share code.
//!
//! **Contract on arcs.** A [`Edge2d::Arc`] states its centre and squared radius `r2`; both of the
//! step's vertices lie on that circle and `r2` is positive. That is the caller's contract
//! (`nacre-ops` checks it once, at construction), `debug_assert`ed here.

use crate::intersect::{
    RingSide, on_segment_2d_rat, orient2d_rat, ray_step_crossing, segments_meet_2d_rat, spike_rat,
};
use nacre_exact::quad::QuadVal;
use nacre_exact::{Orient, Rat};

/// The step leaving a ring vertex: straight to the next vertex, or an arc around `center`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Edge2d {
    Line,
    /// The arc from this step's vertex to the next, around `center`, counter-clockwise when
    /// `ccw`. `r2` is the **squared** radius, derived once at the door and stored — `nacre-ops`'s
    /// step doors and `Ring2d::new` compute `|start − center|²`, which is rational for every circle
    /// through rational points (no square root is taken, so no radius is refused for being
    /// irrational), and check the next vertex against it (`ArcEndOffCircle`). It is the square
    /// because that is what every exact predicate compares against and what the cylinder truth
    /// carries; a radius wanted as a number is `√r2`, a realization. `Ring2d`'s private `edges` is
    /// what keeps a hand-built `Arc` from bypassing the door. A step whose two vertices are one
    /// point is the whole circle (its vertex is the seam).
    Arc {
        center: [Rat; 2],
        r2: Rat,
        ccw: bool,
    },
}

/// A ring as its vertices and the step leaving each: `segs[i]` runs
/// `vertices[i] → vertices[(i + 1) % n]`. `vertices.len() == segs.len()`.
#[derive(Clone, Copy, Debug)]
pub struct MixedRing<'a> {
    pub vertices: &'a [[Rat; 2]],
    pub segs: &'a [Edge2d],
}

/// Checked `Rat`/`QuadVal` arithmetic overflowed: the question has an answer this arithmetic
/// cannot state. Callers refuse by name; nothing here guesses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Undecidable;

type QPt = [QuadVal; 2];

/// One step with its two vertices in hand.
#[derive(Clone, Copy, Debug)]
struct Step {
    from: [Rat; 2],
    to: [Rat; 2],
    seg: Edge2d,
}

impl MixedRing<'_> {
    fn len(&self) -> usize {
        debug_assert_eq!(
            self.vertices.len(),
            self.segs.len(),
            "a ring's vertices and steps pair up"
        );
        self.vertices.len()
    }

    fn step(&self, i: usize) -> Step {
        let n = self.len();
        Step {
            from: self.vertices[i],
            to: self.vertices[(i + 1) % n],
            seg: self.segs[i],
        }
    }
}

impl Step {
    /// A step whose two vertices coincide: a zero-length line, or an arc that is a whole circle.
    fn is_closed(&self) -> bool {
        self.from == self.to
    }
}

// ─── small exact arithmetic ───────────────────────────────────────────────────────────────────

fn sub2(a: [Rat; 2], b: [Rat; 2]) -> Option<[Rat; 2]> {
    Some([a[0].checked_sub(b[0])?, a[1].checked_sub(b[1])?])
}

fn dot2(a: [Rat; 2], b: [Rat; 2]) -> Option<Rat> {
    a[0].checked_mul(b[0])?.checked_add(a[1].checked_mul(b[1])?)
}

fn dist2(a: [Rat; 2], b: [Rat; 2]) -> Option<Rat> {
    let d = sub2(a, b)?;
    dot2(d, d)
}

fn recip(x: Rat) -> Option<Rat> {
    Rat::new(x.denom(), x.numer())
}

fn q(p: [Rat; 2]) -> QPt {
    [QuadVal::from_rat(p[0]), QuadVal::from_rat(p[1])]
}

/// `p + t·d` in ℚ(√c), the point at parameter `t` on the line through `p` along `d`.
fn point_at(p: [Rat; 2], d: [Rat; 2], t: &QuadVal) -> Option<QPt> {
    Some([
        QuadVal::from_rat(p[0]).checked_add(&t.checked_mul_rat(d[0])?)?,
        QuadVal::from_rat(p[1]).checked_add(&t.checked_mul_rat(d[1])?)?,
    ])
}

/// The sign of the orientation determinant `(b − a) × (c − a)` over ℚ(√c) — the twin of
/// [`orient2d_rat`] for points that carry the question's one radical.
fn orient2d_q(a: &QPt, b: &QPt, c: &QPt) -> Option<Orient> {
    let ux = b[0].checked_sub(&a[0])?;
    let uy = b[1].checked_sub(&a[1])?;
    let vx = c[0].checked_sub(&a[0])?;
    let vy = c[1].checked_sub(&a[1])?;
    let det = ux.checked_mul(&vy)?.checked_sub(&uy.checked_mul(&vx)?)?;
    Some(det.sign())
}

fn eq_q(a: &QPt, b: [Rat; 2]) -> Option<bool> {
    let bq = q(b);
    Some(
        a[0].checked_sub(&bq[0])?.sign() == Orient::Zero
            && a[1].checked_sub(&bq[1])?.sign() == Orient::Zero,
    )
}

fn sign_rat(x: Rat) -> Orient {
    match x.cmp(&Rat::from_int(0)) {
        core::cmp::Ordering::Greater => Orient::Positive,
        core::cmp::Ordering::Less => Orient::Negative,
        core::cmp::Ordering::Equal => Orient::Zero,
    }
}

// ─── arcs ─────────────────────────────────────────────────────────────────────────────────────

/// An arc as the predicates read it: centre, squared radius, and its two ends in
/// **counter-clockwise** order (a clockwise step is the same set of points read from its other
/// end).
#[derive(Clone, Copy, Debug)]
struct Arc {
    c: [Rat; 2],
    r2: Rat,
    /// Counter-clockwise start and end. `s == e` is the whole circle.
    s: [Rat; 2],
    e: [Rat; 2],
}

impl Arc {
    fn of(step: &Step) -> Option<Arc> {
        match step.seg {
            Edge2d::Arc { center, r2, ccw } => {
                debug_assert!(r2 > Rat::from_int(0), "an arc's squared radius is positive");
                debug_assert!(
                    [step.from, step.to]
                        .iter()
                        .all(|v| dist2(*v, center).is_none() || dist2(*v, center) == Some(r2)),
                    "an arc's vertices lie on its circle (the caller's contract)"
                );
                let (s, e) = if ccw {
                    (step.from, step.to)
                } else {
                    (step.to, step.from)
                };
                Some(Arc {
                    c: center,
                    r2,
                    s,
                    e,
                })
            }
            Edge2d::Line => None,
        }
    }

    fn is_full(&self) -> bool {
        self.s == self.e
    }

    /// Is `p` — a point **on this circle** — on the arc, ends included?
    ///
    /// The counter-clockwise arc `s → e` is read against the chord: when `e` turns left of the ray
    /// `c → s` the arc is the short way round, and a point is on it iff it is left of (or on)
    /// `c → s` and right of (or on) `c → e`; when `e` turns right the arc is the long way round —
    /// everything but the interior of the short complement; when `c`, `s`, `e` are collinear the
    /// arc is a half circle, the closed half-plane left of `c → s`.
    fn contains(&self, p: &QPt) -> Option<bool> {
        if self.is_full() {
            return Some(true);
        }
        let (c, s, e) = (q(self.c), q(self.s), q(self.e));
        let turn = orient2d_rat(self.c, self.s, self.e);
        let cs = orient2d_q(&c, &s, p)?;
        let ce = orient2d_q(&c, p, &e)?;
        Some(match turn.cmp(&0) {
            core::cmp::Ordering::Greater => cs != Orient::Negative && ce != Orient::Negative,
            core::cmp::Ordering::Less => {
                // Not strictly inside the short complement `e → s`.
                !(orient2d_q(&c, &e, p)? == Orient::Positive
                    && orient2d_q(&c, p, &s)? == Orient::Positive)
            }
            core::cmp::Ordering::Equal => cs != Orient::Negative,
        })
    }

    /// Does the arc leave its end `v` (one of its two vertices, on the probe ray) **strictly
    /// upward**? Counter-clockwise motion at `v` has velocity `(−(v.y − c.y), v.x − c.x)`; a
    /// horizontal departure (`v.x == c.x`) is the top or bottom of the circle, where the arc
    /// bends toward the centre — up from the bottom, down from the top.
    fn departs_up(&self, v: [Rat; 2], ccw_from_v: bool) -> bool {
        match v[0].cmp(&self.c[0]) {
            core::cmp::Ordering::Equal => v[1] < self.c[1],
            core::cmp::Ordering::Greater => ccw_from_v,
            core::cmp::Ordering::Less => !ccw_from_v,
        }
    }
}

/// Where the line `p + t·d` meets the circle — `t` along `d`.
enum LineCircle {
    Miss,
    /// Tangent: one rational parameter.
    Touch(Rat),
    /// Two crossings, `t[0] < t[1]`, each `a ± k√disc` with the same radical.
    Cross([QuadVal; 2]),
}

fn line_circle(p: [Rat; 2], d: [Rat; 2], c: [Rat; 2], r2: Rat) -> Option<LineCircle> {
    // |p + t·d − c|² = r²  ⇒  (d·d)t² + 2(d·w)t + (w·w − r²) = 0,  w = p − c.
    let w = sub2(p, c)?;
    let aa = dot2(d, d)?;
    if aa == Rat::from_int(0) {
        return None; // a zero direction names no line
    }
    let bb = Rat::from_int(2).checked_mul(dot2(d, w)?)?;
    let cc = dot2(w, w)?.checked_sub(r2)?;
    let disc = bb
        .checked_mul(bb)?
        .checked_sub(Rat::from_int(4).checked_mul(aa)?.checked_mul(cc)?)?;
    let inv2a = recip(Rat::from_int(2).checked_mul(aa)?)?;
    let mid = Rat::from_int(0).checked_sub(bb)?.checked_mul(inv2a)?;
    Some(match sign_rat(disc) {
        Orient::Negative => LineCircle::Miss,
        Orient::Zero => LineCircle::Touch(mid),
        Orient::Positive => {
            let neg = Rat::from_int(0).checked_sub(inv2a)?;
            LineCircle::Cross([
                QuadVal::new(mid, neg, disc)?,
                QuadVal::new(mid, inv2a, disc)?,
            ])
        }
    })
}

/// The parameters at which the line meets the circle, as `QuadVal`s (a touch is one rational).
fn line_circle_params(p: [Rat; 2], d: [Rat; 2], c: [Rat; 2], r2: Rat) -> Option<Vec<QuadVal>> {
    Some(match line_circle(p, d, c, r2)? {
        LineCircle::Miss => vec![],
        LineCircle::Touch(t) => vec![QuadVal::from_rat(t)],
        LineCircle::Cross([t0, t1]) => vec![t0, t1],
    })
}

/// Does the closed segment `p → q` meet the arc anywhere except at the points in `skip`
/// (the vertices the two steps share)? Touching counts, as it does for two segments.
fn segment_meets_arc(p: [Rat; 2], qq: [Rat; 2], arc: &Arc, skip: &[[Rat; 2]]) -> Option<bool> {
    let d = sub2(qq, p)?;
    let one = QuadVal::from_rat(Rat::from_int(1));
    for t in line_circle_params(p, d, arc.c, arc.r2)? {
        if t.sign() == Orient::Negative || one.checked_sub(&t)?.sign() == Orient::Negative {
            continue; // outside the segment
        }
        let x = point_at(p, d, &t)?;
        if skipped(&x, skip)? {
            continue;
        }
        if arc.contains(&x)? {
            return Some(true);
        }
    }
    Some(false)
}

fn skipped(x: &QPt, skip: &[[Rat; 2]]) -> Option<bool> {
    for v in skip {
        if eq_q(x, *v)? {
            return Some(true);
        }
    }
    Some(false)
}

/// One circular arc of the plane, stated for [`arcs_share_a_point`]: `start → end` runs
/// counter-clockwise on the circle `(centre, √r2)` — the radius stated as its square, the form
/// the truth holds — and `start == end` is the whole circle. The points are on the circle — the
/// caller's contract, as for every arc in this module.
#[derive(Clone, Copy, Debug)]
pub struct ArcSpec {
    pub centre: [Rat; 2],
    pub r2: Rat,
    pub start: [Rat; 2],
    pub end: [Rat; 2],
}

/// **Do two circular arcs of one plane share a point?** — the predicate the mixed ring already
/// asks of two arc steps ([`arcs_meet`]), opened for the cylinder gate: two lateral faces
/// on parallel axes share a point iff their rims' arcs do in the common cross-section. Touching
/// counts as sharing. `None` is overflow.
pub fn arcs_share_a_point(a: &ArcSpec, b: &ArcSpec) -> Option<bool> {
    let arc = |x: &ArcSpec| Arc {
        c: x.centre,
        r2: x.r2,
        s: x.start,
        e: x.end,
    };
    arcs_meet(&arc(a), &arc(b), &[])
}

/// Do two arcs meet anywhere except at the points in `skip`? Touching counts.
fn arcs_meet(a: &Arc, b: &Arc, skip: &[[Rat; 2]]) -> Option<bool> {
    if a.c == b.c {
        if a.r2 != b.r2 {
            return Some(false); // concentric, distinct radii: disjoint circles
        }
        // One circle: the arcs overlap iff an end of one lies on the other (a full circle
        // contains every point) — except when they are **one arc** whose two ends are both
        // skipped (a ring of two steps running out along a circle and back), where no end is
        // left to witness an overlap that covers the whole arc.
        if !a.is_full() && a.s == b.s && a.e == b.e {
            return Some(true);
        }
        for (x, y) in [(a, b), (b, a)] {
            for v in [x.s, x.e] {
                if skip.contains(&v) {
                    continue;
                }
                if y.contains(&q(v))? {
                    return Some(true);
                }
            }
        }
        // Two whole circles with no vertex to name a point: they coincide entirely.
        return Some(a.is_full() && b.is_full());
    }
    // Distinct centres: the radical line `2(c_b − c_a)·x = |c_b|² − |c_a|² + r_a² − r_b²` carries
    // every common point. Its base `c_a + (h/|d|²)·d`, `h = (|d|² + r_a² − r_b²)/2`, and direction
    // `d⊥` are rational; where it meets circle `a` is the one radical of this question.
    let d = sub2(b.c, a.c)?;
    let d2 = dot2(d, d)?;
    let h = d2
        .checked_add(a.r2)?
        .checked_sub(b.r2)?
        .checked_mul(Rat::new(1, 2)?)?;
    let k = h.checked_mul(recip(d2)?)?;
    let base = [
        a.c[0].checked_add(k.checked_mul(d[0])?)?,
        a.c[1].checked_add(k.checked_mul(d[1])?)?,
    ];
    let perp = [Rat::from_int(0).checked_sub(d[1])?, d[0]];
    for t in line_circle_params(base, perp, a.c, a.r2)? {
        let x = point_at(base, perp, &t)?;
        if skipped(&x, skip)? {
            continue;
        }
        if a.contains(&x)? && b.contains(&x)? {
            return Some(true);
        }
    }
    Some(false)
}

/// Do two steps meet anywhere but at the vertices in `skip`? For two straight steps that share a
/// vertex the only other meeting is a collinear overlap — the spike [`crate::intersect`] names.
fn steps_meet(a: &Step, b: &Step, skip: &[[Rat; 2]]) -> Option<bool> {
    match (Arc::of(a), Arc::of(b)) {
        (None, None) => {
            debug_assert!(
                skip.is_empty(),
                "two straight steps sharing a vertex are the ring predicate's spike, by index"
            );
            Some(segments_meet_2d_rat(a.from, a.to, b.from, b.to))
        }
        (None, Some(arc)) => segment_meets_arc(a.from, a.to, &arc, skip),
        (Some(arc), None) => segment_meets_arc(b.from, b.to, &arc, skip),
        (Some(x), Some(y)) => arcs_meet(&x, &y, skip),
    }
}

// ─── the three predicates ─────────────────────────────────────────────────────────────────────

/// Where a point sits relative to a mixed ring — [`crate::intersect::point_in_ring_2d_rat`] for
/// rings with arcs. Even-odd parity of the rightward ray, half-open as documented on
/// [`crate::intersect::ray_straddle`]; only meaningful on a simple ring.
pub fn point_in_mixed_ring(p: [Rat; 2], ring: MixedRing<'_>) -> Result<RingSide, Undecidable> {
    point_in_mixed_ring_opt(p, ring).ok_or(Undecidable)
}

fn point_in_mixed_ring_opt(p: [Rat; 2], ring: MixedRing<'_>) -> Option<RingSide> {
    let n = ring.len();
    if n == 0 {
        return Some(RingSide::Outside);
    }
    let mut inside = false;
    for i in 0..n {
        let st = ring.step(i);
        match Arc::of(&st) {
            None => {
                let (a, b) = (st.from, st.to);
                let side = orient2d_rat(a, b, p);
                if side == 0 && on_segment_2d_rat(a, b, p) {
                    return Some(RingSide::OnBoundary);
                }
                match ray_step_crossing(
                    sign_rat(a[1].checked_sub(p[1])?),
                    sign_rat(b[1].checked_sub(p[1])?),
                    sign_rat(Rat::from_int(i128::from(side))),
                ) {
                    Some(true) => inside = !inside,
                    Some(false) => {}
                    None => return Some(RingSide::OnBoundary),
                }
            }
            Some(arc) => {
                let r2 = arc.r2;
                // The probe on the arc is boundary, whatever the ray says.
                if dist2(p, arc.c)? == r2 && arc.contains(&q(p))? {
                    return Some(RingSide::OnBoundary);
                }
                // The ray's line `y = p.y` against the circle: `x = c.x ± √(r² − dy²)`.
                let dy = p[1].checked_sub(arc.c[1])?;
                let disc = r2.checked_sub(dy.checked_mul(dy)?)?;
                let xs: Vec<QuadVal> = match sign_rat(disc) {
                    Orient::Negative => vec![],
                    Orient::Zero => vec![QuadVal::from_rat(arc.c[0])],
                    Orient::Positive => {
                        let neg = Rat::from_int(-1);
                        vec![
                            QuadVal::new(arc.c[0], neg, disc)?,
                            QuadVal::new(arc.c[0], Rat::from_int(1), disc)?,
                        ]
                    }
                };
                let tangent = sign_rat(disc) == Orient::Zero;
                for x in xs {
                    let pt: QPt = [x, QuadVal::from_rat(p[1])];
                    if !arc.contains(&pt)? {
                        continue;
                    }
                    // Strictly right of the probe; the probe itself was boundary above.
                    if x.checked_sub(&QuadVal::from_rat(p[0]))?.sign() != Orient::Positive {
                        continue;
                    }
                    // Which vertex of the *step* is this point, if any? The step runs
                    // `st.from → st.to`; the arc's ccw reading may have swapped them.
                    let at_from = eq_q(&pt, st.from)?;
                    let at_to = eq_q(&pt, st.to)?;
                    let ccw = matches!(st.seg, Edge2d::Arc { ccw: true, .. });
                    let mut toggles = 0u8;
                    if at_from {
                        // Leaves `from` upward?
                        if arc.departs_up(st.from, ccw) {
                            toggles += 1;
                        }
                    }
                    if at_to {
                        // Arrives at `to` from above ⇔ the reversed step leaves `to` upward.
                        if arc.departs_up(st.to, !ccw) {
                            toggles += 1;
                        }
                    }
                    if !at_from && !at_to && !tangent {
                        toggles += 1; // an interior crossing of the ray's line
                    }
                    if toggles % 2 == 1 {
                        inside = !inside;
                    }
                }
            }
        }
    }
    Some(if inside {
        RingSide::Inside
    } else {
        RingSide::Outside
    })
}

/// Do two mixed rings touch or cross — [`crate::intersect::rings_cross_rat`] for rings with arcs.
pub fn mixed_rings_cross(a: MixedRing<'_>, b: MixedRing<'_>) -> Result<bool, Undecidable> {
    let (n, m) = (a.len(), b.len());
    for i in 0..n {
        for j in 0..m {
            if steps_meet(&a.step(i), &b.step(j), &[]).ok_or(Undecidable)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Does a mixed ring meet itself — [`crate::intersect::ring_self_intersection_rat`] for rings
/// with arcs. Two steps that share a vertex may meet there and nowhere else; a straight step of
/// zero length, or a whole-circle step in a ring of more than one step, is reported against
/// itself. Returns the two offending step indices.
pub fn mixed_ring_self_intersection(
    ring: MixedRing<'_>,
) -> Result<Option<(usize, usize)>, Undecidable> {
    let n = ring.len();
    if n == 0 {
        return Ok(None);
    }
    for i in 0..n {
        let st = ring.step(i);
        let whole_circle = st.is_closed() && matches!(st.seg, Edge2d::Arc { .. });
        if st.is_closed() && !(whole_circle && n == 1) {
            return Ok(Some((i, i)));
        }
    }
    let pt = |i: usize| ring.vertices[i % n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (ring.step(i), ring.step(j));
            let hit = match (Arc::of(&a), Arc::of(&b)) {
                // Straight against straight: the polygon predicate's spelling, pair for pair —
                // adjacent steps meet elsewhere only as a collinear spike, the wrap-around pair
                // is read from the shared vertex `0`, everything else is the segment test.
                (None, None) => {
                    if j == i + 1 {
                        spike_rat(pt(i), pt(i + 1), pt(i + 2))
                    } else if i == 0 && j == n - 1 {
                        spike_rat(pt(1), pt(0), pt(n - 1))
                    } else {
                        segments_meet_2d_rat(pt(i), pt(i + 1), pt(j), pt(j + 1))
                    }
                }
                // An arc is involved: the steps may meet at the vertices they share **by index**
                // and nowhere else.
                _ => {
                    let mut shared: Vec<[Rat; 2]> = Vec::new();
                    if j == i + 1 {
                        shared.push(pt(i + 1));
                    }
                    if i == 0 && j == n - 1 && !shared.contains(&pt(0)) {
                        shared.push(pt(0));
                    }
                    steps_meet(&a, &b, &shared).ok_or(Undecidable)?
                }
            };
            if hit {
                return Ok(Some((i, j)));
            }
        }
    }
    Ok(None)
}

#[cfg(test)]
#[path = "tests/mixed.rs"]
mod tests;
