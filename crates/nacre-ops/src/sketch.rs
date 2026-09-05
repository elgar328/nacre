//! Turning loose edges into profiles — the nesting policy on top of geom's exact predicates.
//!
//! A sketch is drawn as closed paths; which of them are material, which are holes, and which are
//! islands inside holes is **not** declared by the author but decided from containment, the way
//! every CAD sketcher does it. The classification itself (is this point inside that ring? do two
//! rings meet?) is a robustness-sensitive sign question and lives in `nacre-geom` (design §1
//! isolates that — [`nacre_geom::intersect`] for polygons, [`nacre_geom::mixed`] once arcs are
//! in); what is here is the policy — depth parity — and the assembly.
//!
//! **Edges are stated, not derived.** A straight edge is its two ends; an arc is its centre, its
//! radius and its two ends, all rational, and the constructors here are the only place a
//! coordinate is computed — the end of a quarter-turn arc is the start rotated in `Rat`, never an
//! f64 that was rounded on the way. The f64 doors take what the author *wrote* (their decimals
//! become the rationals they spell); the `Rat` doors take what a caller computed exactly.

use crate::{Profile2d, Ring2d};
use nacre_geom::intersect::RingSide;
use nacre_geom::mixed::{
    MixedRing, Seg2d, Undecidable, mixed_ring_self_intersection, mixed_rings_cross,
    point_in_mixed_ring,
};
use nacre_math::Point2;
use nacre_scalar::{Rat, rat_sqrt_exact};

/// Why a set of rings or edges is not a valid set of profiles.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum SketchError {
    /// A ring with fewer than three points encloses nothing.
    DegenerateRing { ring: usize },
    /// Two rings touch or cross. A hole must lie strictly inside its outer ring and strictly
    /// outside its siblings; anything else has no unambiguous inside.
    RingsMeet { a: usize, b: usize },
    /// A straight edge starts where it ends.
    ZeroLengthEdge { edge: usize },
    /// The same edge appears twice (the same set of points, either way round).
    DuplicateEdge { a: usize, b: usize },
    /// A chain ran out of edges before closing. `at` is the endpoint left dangling and `gap` the
    /// distance to the nearest other free endpoint — reported because "you meant to close this
    /// and missed by 1e-9" is the likely story, and the kernel will not close it for you: an
    /// endpoint either *is* the same point or is not (overview 절대원칙 4).
    OpenChain { at: [f64; 2], gap: Option<f64> },
    /// Three or more edges meet at one point, so the chain has no unambiguous continuation.
    BranchingVertex { at: [f64; 2] },
    /// A ring meets itself — it crosses, touches, or doubles back over its own edge. Such a ring
    /// has no unambiguous inside, so it cannot be sorted into a profile at all.
    ///
    /// Reported as the two offending edges' **chord midpoints**, not as indices: [`from_edges`]
    /// chains the edges into rings in walk order, so a ring's edge index says nothing about where
    /// the author's input went wrong. Points are what the rest of this enum reports too.
    RingSelfIntersects { ring: usize, at: [[f64; 2]; 2] },
    /// A coordinate outside the decimal window (`~1e38` above, `~1e-22` below for a full-width
    /// value) has no rational truth for the kernel to keep — the sketch-layer twin of
    /// `OpError::ProfileOutsideDecimalWindow`. `at` is the offending point.
    OutsideDecimalWindow { at: [f64; 2] },
    /// An arc's radius is not rational: `|start − center|²` is not the square of a rational, so
    /// the cylinder this arc would stand has no exact radius. A start axis-aligned from the centre
    /// always passes; so does any `3-4-5`-like point.
    ArcRadiusNotRational { center: [f64; 2], start: [f64; 2] },
    /// An arc's end is not on the circle its centre and start define.
    ArcEndOffCircle {
        center: [f64; 2],
        start: [f64; 2],
        end: [f64; 2],
    },
    /// An arc from a point back to itself is a whole circle, which [`Edge2d::circle`] states;
    /// this door asked for a proper arc.
    ZeroLengthArc { at: [f64; 2] },
    /// A quarter-turn count outside `±1..=±3` — `0` is no arc, `±4` is a whole circle
    /// ([`Edge2d::circle`]).
    ArcTurnsOutOfRange { turns: i32 },
    /// A circle or arc with a radius that is not positive.
    NonPositiveRadius { center: [f64; 2], radius: f64 },
    /// Checked `Rat` arithmetic overflowed while classifying: the question has an answer the
    /// kernel cannot state here. Refused by name; nothing guesses.
    Undecidable,
}

/// A 2-D sketch edge, stated exactly. The `Line`/`Arc` split is the segment vocabulary the
/// kernel's prism builder and the profile predicates read; `Nurbs` will join it in M7.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Edge2d {
    /// Straight, `from → to`.
    Line { from: [Rat; 2], to: [Rat; 2] },
    /// Circular, `start → end` around `center`, counter-clockwise when `ccw`. `radius` is the
    /// rational `|start − center|`, checked at construction along with `end` being on the circle.
    /// `start == end` is the whole circle, whose one vertex is the seam.
    Arc {
        center: [Rat; 2],
        radius: Rat,
        start: [Rat; 2],
        end: [Rat; 2],
        ccw: bool,
    },
}

fn lift(p: Point2) -> Result<[Rat; 2], SketchError> {
    match (Rat::from_decimal(p[0]), Rat::from_decimal(p[1])) {
        (Some(x), Some(y)) => Ok([x, y]),
        _ => Err(SketchError::OutsideDecimalWindow { at: p.as_array() }),
    }
}

fn f2(p: [Rat; 2]) -> [f64; 2] {
    [p[0].to_f64(), p[1].to_f64()]
}

fn sub2(a: [Rat; 2], b: [Rat; 2]) -> Option<[Rat; 2]> {
    Some([a[0].checked_sub(b[0])?, a[1].checked_sub(b[1])?])
}

fn dist2(a: [Rat; 2], b: [Rat; 2]) -> Option<Rat> {
    let d = sub2(a, b)?;
    d[0].checked_mul(d[0])?.checked_add(d[1].checked_mul(d[1])?)
}

/// `center + rot(v, quarter turns)` — a quarter turn counter-clockwise is `(x, y) → (−y, x)`,
/// exact in `Rat`.
fn turned(center: [Rat; 2], v: [Rat; 2], quarter_turns: i32) -> Option<[Rat; 2]> {
    let zero = Rat::from_int(0);
    let mut d = v;
    let ccw = quarter_turns > 0;
    for _ in 0..quarter_turns.unsigned_abs() {
        d = if ccw {
            [zero.checked_sub(d[1])?, d[0]]
        } else {
            [d[1], zero.checked_sub(d[0])?]
        };
    }
    Some([center[0].checked_add(d[0])?, center[1].checked_add(d[1])?])
}

impl Edge2d {
    /// A straight edge between two written points.
    pub fn line(from: Point2, to: Point2) -> Result<Edge2d, SketchError> {
        Ok(Edge2d::Line {
            from: lift(from)?,
            to: lift(to)?,
        })
    }

    /// A straight edge between two exact points — the door for computed coordinates.
    pub fn line_rat(from: [Rat; 2], to: [Rat; 2]) -> Edge2d {
        Edge2d::Line { from, to }
    }

    /// The arc from `start` around `center` through `quarter_turns` right angles — positive is
    /// counter-clockwise (`+x → +y` in the sketch's own frame), `±1..=±3`. The end is the start
    /// **rotated in `Rat`**: it lies on the circle by construction, and no computed value passes
    /// through f64. Today's whole arc vocabulary is this door and [`Edge2d::circle`]; an arbitrary
    /// angle is a point the kernel would have to *name* (a start turned by θ), which is a later
    /// vocabulary, not a rounding.
    pub fn arc_turns(
        center: Point2,
        start: Point2,
        quarter_turns: i32,
    ) -> Result<Edge2d, SketchError> {
        if !(1..=3).contains(&quarter_turns.unsigned_abs()) {
            return Err(SketchError::ArcTurnsOutOfRange {
                turns: quarter_turns,
            });
        }
        let (c, s) = (lift(center)?, lift(start)?);
        let radius = radius_of(c, s)?;
        let v = sub2(s, c).ok_or(SketchError::Undecidable)?;
        let end = turned(c, v, quarter_turns).ok_or(SketchError::Undecidable)?;
        Ok(Edge2d::Arc {
            center: c,
            radius,
            start: s,
            end,
            ccw: quarter_turns > 0,
        })
    }

    /// The whole circle of `radius` about `center`, counter-clockwise, its seam at `center + (r, 0)`
    /// — the same point the cylinder primitive seams at (`+ref_dir`), which is what lets the two
    /// roads state one solid.
    pub fn circle(center: Point2, radius: f64) -> Result<Edge2d, SketchError> {
        let c = lift(center)?;
        let r = Rat::from_decimal(radius)
            .ok_or(SketchError::OutsideDecimalWindow { at: [radius, 0.0] })?;
        if r <= Rat::from_int(0) {
            return Err(SketchError::NonPositiveRadius {
                center: center.as_array(),
                radius,
            });
        }
        let seam = [c[0].checked_add(r).ok_or(SketchError::Undecidable)?, c[1]];
        Ok(Edge2d::Arc {
            center: c,
            radius: r,
            start: seam,
            end: seam,
            ccw: true,
        })
    }

    /// A proper arc from exact data — the door for computed coordinates (a fillet's tangent
    /// points and centre, say). Checks what the f64 doors guarantee by construction: `end` on the
    /// circle, the radius rational, `start ≠ end`.
    pub fn arc_rat(
        center: [Rat; 2],
        start: [Rat; 2],
        end: [Rat; 2],
        ccw: bool,
    ) -> Result<Edge2d, SketchError> {
        if start == end {
            return Err(SketchError::ZeroLengthArc { at: f2(start) });
        }
        let radius = radius_of(center, start)?;
        let r2 = radius.checked_mul(radius).ok_or(SketchError::Undecidable)?;
        if dist2(end, center).ok_or(SketchError::Undecidable)? != r2 {
            return Err(SketchError::ArcEndOffCircle {
                center: f2(center),
                start: f2(start),
                end: f2(end),
            });
        }
        Ok(Edge2d::Arc {
            center,
            radius,
            start,
            end,
            ccw,
        })
    }

    /// The edge's first vertex.
    pub fn start(&self) -> [Rat; 2] {
        match *self {
            Edge2d::Line { from, .. } => from,
            Edge2d::Arc { start, .. } => start,
        }
    }

    /// The edge's last vertex.
    pub fn end(&self) -> [Rat; 2] {
        match *self {
            Edge2d::Line { to, .. } => to,
            Edge2d::Arc { end, .. } => end,
        }
    }

    /// The step this edge is when walked `start → end`.
    pub fn seg(&self) -> Seg2d {
        match *self {
            Edge2d::Line { .. } => Seg2d::Line,
            Edge2d::Arc {
                center,
                radius,
                ccw,
                ..
            } => Seg2d::Arc {
                center,
                radius,
                ccw,
            },
        }
    }

    /// The step this edge is when walked `end → start`.
    fn seg_reversed(&self) -> Seg2d {
        match self.seg() {
            Seg2d::Line => Seg2d::Line,
            Seg2d::Arc {
                center,
                radius,
                ccw,
            } => Seg2d::Arc {
                center,
                radius,
                ccw: !ccw,
            },
        }
    }

    /// An arc from a point to itself: one vertex, the seam, and the whole circle.
    pub fn is_whole_circle(&self) -> bool {
        matches!(self, Edge2d::Arc { start, end, .. } if start == end)
    }

    /// The same set of points, walked either way.
    fn same_points(&self, other: &Edge2d) -> bool {
        match (self, other) {
            (Edge2d::Line { from: a, to: b }, Edge2d::Line { from: c, to: d }) => {
                (a == c && b == d) || (a == d && b == c)
            }
            (
                Edge2d::Arc {
                    center: c1,
                    radius: r1,
                    start: s1,
                    end: e1,
                    ccw: w1,
                },
                Edge2d::Arc {
                    center: c2,
                    radius: r2,
                    start: s2,
                    end: e2,
                    ccw: w2,
                },
            ) => {
                c1 == c2
                    && r1 == r2
                    && ((s1 == s2 && e1 == e2 && w1 == w2) || (s1 == e2 && e1 == s2 && w1 != w2))
            }
            _ => false,
        }
    }
}

/// The rational radius `|start − center|`, or the named refusal.
fn radius_of(center: [Rat; 2], start: [Rat; 2]) -> Result<Rat, SketchError> {
    let r2 = dist2(start, center).ok_or(SketchError::Undecidable)?;
    match rat_sqrt_exact(r2) {
        Some(r) if r > Rat::from_int(0) => Ok(r),
        Some(_) => Err(SketchError::NonPositiveRadius {
            center: f2(center),
            radius: 0.0,
        }),
        None => Err(SketchError::ArcRadiusNotRational {
            center: f2(center),
            start: f2(start),
        }),
    }
}

/// Sort closed polygon rings into profiles by containment depth — the straight-only door.
///
/// Depth is how many rings a ring lies inside. **Even depth is material, odd depth is a hole** —
/// so a ring inside a hole is material again (an island), and each becomes a profile of its own.
/// A hole is attached to its *immediate* container, the deepest ring that contains it.
///
/// There is no fill-rule parameter. Even-odd is the rule; offering the choice would mean asking
/// callers to understand a distinction the syntax deliberately hides, and for the simple,
/// non-touching rings this accepts, even-odd and non-zero agree anyway.
pub fn from_rings(rings: Vec<Vec<Point2>>) -> Result<Vec<Profile2d>, SketchError> {
    for (i, r) in rings.iter().enumerate() {
        if r.len() < 3 {
            return Err(SketchError::DegenerateRing { ring: i });
        }
    }
    // Lift to the rational truth first: classification judges the decimals the author wrote,
    // the same truth `Profile2d` stores — the f64 carriers can hold a *different* sign near
    // degeneracy, and it is the decimal that is the author's meaning.
    let lifted: Vec<Ring2d> = rings
        .iter()
        .map(|r| {
            r.iter()
                .map(|p| lift(*p))
                .collect::<Result<Vec<_>, _>>()
                .map(Ring2d::polygon)
        })
        .collect::<Result<_, _>>()?;
    classify(lifted)
}

/// Sort normalized rings — straight or arc-bearing — into profiles by containment depth. The
/// policy of [`from_rings`]; the predicates are [`nacre_geom::mixed`]'s.
fn classify(rings: Vec<Ring2d>) -> Result<Vec<Profile2d>, SketchError> {
    for (i, r) in rings.iter().enumerate() {
        if r.is_empty() || (r.is_polygon() && r.len() < 3) {
            return Err(SketchError::DegenerateRing { ring: i });
        }
    }
    let undecidable = |_: Undecidable| SketchError::Undecidable;
    // Simplicity comes before nesting, and not merely for a better message: the parity that
    // decides containment only means "inside" on a simple ring.
    for (i, r) in rings.iter().enumerate() {
        if let Some((a, b)) = mixed_ring_self_intersection(r.mixed()).map_err(undecidable)? {
            // The chord midpoint is diagnostic only, so f64 is the right material.
            return Err(SketchError::RingSelfIntersects {
                ring: i,
                at: [r.chord_midpoint(a), r.chord_midpoint(b)],
            });
        }
    }
    let n = rings.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if mixed_rings_cross(rings[i].mixed(), rings[j].mixed()).map_err(undecidable)? {
                return Err(SketchError::RingsMeet { a: i, b: j });
            }
        }
    }

    // `inside[i]` = the rings containing ring `i`. The rings are disjoint by the check above, so
    // any vertex of `i` answers for the whole ring.
    let mut inside: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if j != i
                && point_in_mixed_ring(rings[i].vertices()[0], rings[j].mixed())
                    .map_err(undecidable)?
                    == RingSide::Inside
            {
                inside[i].push(j);
            }
        }
    }

    // A ring's parent is the deepest ring containing it; holes hang off the material ring they
    // are cut from, not off some distant ancestor.
    let depth = |i: usize| inside[i].len();
    let parent: Vec<Option<usize>> = (0..n)
        .map(|i| inside[i].iter().copied().max_by_key(|&j| depth(j)))
        .collect();

    let mut out = Vec::new();
    for i in 0..n {
        if depth(i) % 2 != 0 {
            continue; // odd depth: a hole, carried by its parent below
        }
        let holes: Vec<Ring2d> = (0..n)
            .filter(|&h| depth(h) % 2 == 1 && parent[h] == Some(i))
            .map(|h| rings[h].clone())
            .collect();
        out.push(Profile2d::from_normalized_rings(rings[i].clone(), holes));
    }
    Ok(out)
}

/// Chain loose edges into closed rings, then sort those into profiles ([`from_rings`]).
///
/// The edges may arrive in any order and pointing either way — this is the shape a generator or
/// an imported file produces. Endpoints must **coincide exactly**; nothing here snaps a near-miss
/// shut. Tolerance is for intersections the kernel *discovers*, never for what a caller
/// constructs (overview 절대원칙 4), and a sketch that silently welds a 1e-9 gap is a sketch
/// whose meaning the kernel invented.
///
/// An arc walked against its stated direction is the same points the other way round, so the
/// ring it lands in carries it with `ccw` flipped. A whole circle is a ring by itself: it has no
/// end to chain. The rings then take their normal form ([`Ring2d::normalized`]) before the
/// predicates read them.
pub fn from_edges(edges: Vec<Edge2d>) -> Result<Vec<Profile2d>, SketchError> {
    for (i, e) in edges.iter().enumerate() {
        if matches!(e, Edge2d::Line { .. }) && e.start() == e.end() {
            return Err(SketchError::ZeroLengthEdge { edge: i });
        }
    }
    for i in 0..edges.len() {
        for j in (i + 1)..edges.len() {
            if edges[i].same_points(&edges[j]) {
                return Err(SketchError::DuplicateEdge { a: i, b: j });
            }
        }
    }

    let mut rings: Vec<Ring2d> = Vec::new();
    // Whole circles first, as they come: each is its own ring.
    let mut open: Vec<usize> = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        if e.is_whole_circle() {
            rings.push(Ring2d::normalized(vec![e.start()], vec![e.seg()]));
        } else {
            open.push(i);
        }
    }

    // Endpoint → the edges touching it. Exactly two is a corner; more has no single continuation.
    let mut at: Vec<([Rat; 2], Vec<usize>)> = Vec::new();
    let slot = |at: &mut Vec<([Rat; 2], Vec<usize>)>, p: [Rat; 2]| -> usize {
        match at.iter().position(|(q, _)| *q == p) {
            Some(k) => k,
            None => {
                at.push((p, Vec::new()));
                at.len() - 1
            }
        }
    };
    let ends: Vec<(usize, usize)> = open
        .iter()
        .map(|&i| {
            (
                slot(&mut at, edges[i].start()),
                slot(&mut at, edges[i].end()),
            )
        })
        .collect();
    for (k, (s, t)) in ends.iter().enumerate() {
        at[*s].1.push(k);
        at[*t].1.push(k);
    }
    // Most specific first: a branch says *which* junction is ambiguous, while a dangling end only
    // says the outline is open — and a branch always leaves an odd end somewhere, so checking in
    // the other order would report the vaguer of the two. (`check_result_topology` orders its
    // defects the same way.)
    if let Some((p, _)) = at.iter().find(|(_, touching)| touching.len() > 2) {
        return Err(SketchError::BranchingVertex { at: f2(*p) });
    }
    if let Some((p, _)) = at.iter().find(|(_, touching)| touching.len() == 1) {
        return Err(SketchError::OpenChain {
            at: f2(*p),
            gap: nearest_free_gap(*p, &at),
        });
    }

    // Walk each cycle: from an unused edge, hop endpoint to endpoint until back at the start,
    // recording each step in the direction it was walked.
    let mut used = vec![false; open.len()];
    for seed in 0..open.len() {
        if used[seed] {
            continue;
        }
        let (first, _) = ends[seed];
        let mut here = first;
        let mut k = seed;
        let (mut vertices, mut segs) = (vec![at[first].0], Vec::new());
        loop {
            used[k] = true;
            let e = &edges[open[k]];
            let (s, t) = ends[k];
            let (seg, next) = if s == here {
                (e.seg(), t)
            } else {
                (e.seg_reversed(), s)
            };
            segs.push(seg);
            here = next;
            let Some(&n) = at[here].1.iter().find(|&&e| !used[e]) else {
                break;
            };
            vertices.push(at[here].0);
            k = n;
        }
        // The walk returns to its start, so the vertex it closed on is the first one.
        debug_assert_eq!(here, first, "a closed chain returns to its seed");
        rings.push(Ring2d::normalized(vertices, segs));
    }
    classify(rings)
}

/// The distance from `p` to the nearest *other* endpoint that is also dangling — the size of the
/// gap the author probably meant to close.
fn nearest_free_gap(p: [Rat; 2], at: &[([Rat; 2], Vec<usize>)]) -> Option<f64> {
    let pf = f2(p);
    at.iter()
        .filter(|(q, touching)| touching.len() == 1 && *q != p)
        .map(|(q, _)| {
            let qf = f2(*q);
            ((qf[0] - pf[0]).powi(2) + (qf[1] - pf[1]).powi(2)).sqrt()
        })
        .min_by(|a, b| a.total_cmp(b))
}

/// A ring's view for the predicates.
impl Ring2d {
    pub(crate) fn mixed(&self) -> MixedRing<'_> {
        MixedRing {
            vertices: self.vertices(),
            segs: self.segs(),
        }
    }

    /// The f64 midpoint of step `i`'s chord — diagnostic material for a message.
    fn chord_midpoint(&self, i: usize) -> [f64; 2] {
        let n = self.len();
        let (p, q) = (f2(self.vertices()[i]), f2(self.vertices()[(i + 1) % n]));
        [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sq(a: f64, b: f64) -> Vec<Point2> {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    }

    #[test]
    fn a_lone_ring_is_one_profile() {
        let p = from_rings(vec![sq(0.0, 4.0)]).unwrap();
        assert_eq!(p.len(), 1);
        assert!(p[0].holes().is_empty());
    }

    /// The ring order must not matter: the hole is recognised by containment, not by position.
    #[test]
    fn a_ring_inside_a_ring_is_a_hole_either_way() {
        for rings in [
            vec![sq(0.0, 4.0), sq(1.0, 3.0)],
            vec![sq(1.0, 3.0), sq(0.0, 4.0)],
        ] {
            let p = from_rings(rings).unwrap();
            assert_eq!(p.len(), 1, "one body");
            assert_eq!(p[0].holes().len(), 1, "one hole");
            assert_eq!(p[0].outer().vertices().len(), 4);
        }
    }

    /// Depth 2 is material again — the island is its own body, and it does *not* become a hole of
    /// the outer ring.
    #[test]
    fn a_ring_inside_a_hole_is_an_island() {
        let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0)]).unwrap();
        assert_eq!(p.len(), 2, "outer body + island");
        let outer = p
            .iter()
            .find(|q| q.outer().vertices()[0][0] == Rat::from_int(0))
            .unwrap();
        let island = p
            .iter()
            .find(|q| q.outer().vertices()[0][0] == Rat::from_int(2))
            .unwrap();
        assert_eq!(outer.holes().len(), 1, "the depth-1 ring is its hole");
        assert!(island.holes().is_empty(), "nothing inside the island");
    }

    /// A hole belongs to the ring it is actually cut from, not to a distant ancestor.
    #[test]
    fn a_hole_attaches_to_its_immediate_container() {
        // outer(0..9) ⊃ hole(1..8) ⊃ island(2..7) ⊃ island-hole(3..6)
        let p = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0), sq(3.0, 6.0)]).unwrap();
        assert_eq!(p.len(), 2);
        let island = p
            .iter()
            .find(|q| q.outer().vertices()[0][0] == Rat::from_int(2))
            .unwrap();
        assert_eq!(island.holes().len(), 1, "the depth-3 ring is the island's");
        assert_eq!(island.holes()[0].vertices()[0][0], Rat::from_int(3));
    }

    /// Two separate bodies, neither inside the other.
    #[test]
    fn disjoint_rings_are_separate_bodies() {
        let p = from_rings(vec![sq(0.0, 1.0), sq(5.0, 6.0)]).unwrap();
        assert_eq!(p.len(), 2);
        assert!(p.iter().all(|q| q.holes().is_empty()));
    }

    #[test]
    fn overlapping_rings_are_rejected() {
        let err = from_rings(vec![sq(0.0, 4.0), sq(2.0, 6.0)]).unwrap_err();
        assert!(matches!(err, SketchError::RingsMeet { .. }), "{err:?}");
    }

    /// Touching counts as meeting: a hole flush against its outer ring has no strict inside.
    #[test]
    fn touching_rings_are_rejected() {
        let err = from_rings(vec![sq(0.0, 4.0), sq(4.0, 8.0)]).unwrap_err();
        assert!(matches!(err, SketchError::RingsMeet { .. }), "{err:?}");
    }

    #[test]
    fn a_two_point_ring_is_rejected() {
        let err = from_rings(vec![vec![
            Point2::from_array([0.0, 0.0]),
            Point2::from_array([1.0, 0.0]),
        ]])
        .unwrap_err();
        assert_eq!(err, SketchError::DegenerateRing { ring: 0 });
    }

    fn seg(a: [f64; 2], b: [f64; 2]) -> Edge2d {
        Edge2d::line(Point2::from_array(a), Point2::from_array(b)).unwrap()
    }

    /// Edges in scrambled order, some drawn backwards — the shape a generator emits.
    #[test]
    fn loose_edges_chain_into_a_ring() {
        let p = from_edges(vec![
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([0.0, 4.0], [0.0, 0.0]),
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([0.0, 4.0], [4.0, 4.0]), // drawn right-to-left
        ])
        .unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].outer().vertices().len(), 4);
        assert!(p[0].holes().is_empty());
    }

    /// Two cycles at once, one inside the other: chaining and nesting compose.
    #[test]
    fn two_chains_become_a_profile_with_a_hole() {
        let ring = |a: f64, b: f64| {
            vec![
                seg([a, a], [b, a]),
                seg([b, a], [b, b]),
                seg([b, b], [a, b]),
                seg([a, b], [a, a]),
            ]
        };
        let mut edges = ring(0.0, 4.0);
        edges.extend(ring(1.0, 3.0));
        let p = from_edges(edges).unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].holes().len(), 1);
    }

    /// A gap is not closed for the author — it is measured and handed back.
    #[test]
    fn an_unclosed_chain_reports_the_gap() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([4.0, 4.0], [0.0, 4.0]),
            seg([0.0, 4.0], [0.0, 0.25]), // 0.25 short of the start
        ])
        .unwrap_err();
        let SketchError::OpenChain { gap, .. } = err else {
            panic!("{err:?}");
        };
        assert!((gap.unwrap() - 0.25).abs() < 1e-12, "{gap:?}");
    }

    #[test]
    fn a_zero_length_edge_is_rejected() {
        let err = from_edges(vec![seg([1.0, 1.0], [1.0, 1.0])]).unwrap_err();
        assert_eq!(err, SketchError::ZeroLengthEdge { edge: 0 });
    }

    #[test]
    fn a_repeated_edge_is_rejected() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [0.0, 0.0]), // the same segment, reversed
            seg([4.0, 0.0], [4.0, 4.0]),
        ])
        .unwrap_err();
        assert!(matches!(err, SketchError::DuplicateEdge { .. }), "{err:?}");
    }

    /// A T-junction has no single continuation, so the walk refuses rather than picking one.
    #[test]
    fn a_branching_vertex_is_rejected() {
        let err = from_edges(vec![
            seg([0.0, 0.0], [4.0, 0.0]),
            seg([4.0, 0.0], [4.0, 4.0]),
            seg([4.0, 0.0], [8.0, 0.0]),
        ])
        .unwrap_err();
        assert!(
            matches!(err, SketchError::BranchingVertex { .. }),
            "{err:?}"
        );
    }

    /// Containment must not depend on the rings' winding — the author draws in whatever direction
    /// is convenient, and `build_prism` fixes winding later anyway.
    #[test]
    fn winding_does_not_affect_nesting() {
        let mut hole = sq(1.0, 3.0);
        hole.reverse();
        let p = from_rings(vec![sq(0.0, 4.0), hole]).unwrap();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].holes().len(), 1);
    }

    fn pt(x: f64, y: f64) -> Point2 {
        Point2::from_array([x, y])
    }

    #[test]
    fn a_circle_is_a_ring_of_one_vertex_at_its_seam() {
        let p = from_edges(vec![Edge2d::circle(pt(1.0, 2.0), 3.0).unwrap()]).unwrap();
        assert_eq!(p.len(), 1);
        let outer = p[0].outer();
        assert_eq!(outer.len(), 1);
        assert_eq!(
            outer.vertices()[0],
            [Rat::from_int(4), Rat::from_int(2)],
            "seam at +x"
        );
        assert!(matches!(outer.segs()[0], Seg2d::Arc { ccw: true, .. }));
        assert!(p[0].has_arcs());
        p[0].check().unwrap();
    }

    #[test]
    fn concentric_circles_nest_like_squares() {
        let c = |r: f64| Edge2d::circle(pt(0.0, 0.0), r).unwrap();
        let p = from_edges(vec![c(4.0), c(20.0), c(12.0)]).unwrap();
        assert_eq!(p.len(), 2, "a ring and an island");
        let ring = p
            .iter()
            .find(|q| q.outer().vertices()[0][0] == Rat::from_int(20))
            .unwrap();
        assert_eq!(ring.holes().len(), 1);
        assert_eq!(ring.holes()[0].vertices()[0][0], Rat::from_int(12));
        let island = p
            .iter()
            .find(|q| q.outer().vertices()[0][0] == Rat::from_int(4))
            .unwrap();
        assert!(island.holes().is_empty());
    }

    #[test]
    fn a_slot_chains_lines_and_arcs_whichever_way_they_were_drawn() {
        // Centres (0,0) and (30,0), r 5: two lines and two half circles, scrambled, one arc
        // stated clockwise — the walk turns it round.
        let edges = vec![
            Edge2d::arc_turns(pt(30.0, 0.0), pt(30.0, -5.0), 2).unwrap(), // (30,−5) → (30,5)
            Edge2d::line(pt(0.0, 5.0), pt(30.0, 5.0)).unwrap(),
            Edge2d::arc_turns(pt(0.0, 0.0), pt(0.0, -5.0), -2).unwrap(), // (0,−5) → (0,5), clockwise
            Edge2d::line(pt(0.0, -5.0), pt(30.0, -5.0)).unwrap(),
        ];
        let p = from_edges(edges).unwrap();
        assert_eq!(p.len(), 1);
        let o = p[0].outer();
        assert_eq!(o.len(), 4);
        assert_eq!(
            o.segs()
                .iter()
                .filter(|s| matches!(s, Seg2d::Arc { .. }))
                .count(),
            2
        );
        p[0].check().unwrap();
    }

    #[test]
    fn two_half_circles_are_one_circle_in_normal_form() {
        let e = vec![
            Edge2d::arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 2).unwrap(), // (5,0) → (−5,0)
            Edge2d::arc_turns(pt(0.0, 0.0), pt(-5.0, 0.0), 2).unwrap(), // (−5,0) → (5,0)
        ];
        let p = from_edges(e).unwrap();
        assert_eq!(p[0].outer().len(), 1, "merged into a whole circle");
        assert_eq!(
            p[0].outer().vertices()[0],
            [Rat::from_int(5), Rat::from_int(0)]
        );
        let c = from_edges(vec![Edge2d::circle(pt(0.0, 0.0), 5.0).unwrap()]).unwrap();
        assert_eq!(p, c, "the same profile `circle` states");
    }

    #[test]
    fn equal_fillets_on_a_short_side_merge_into_a_half_circle() {
        // A 30 × 10 rectangle rounded r = 5 at every corner: the short sides vanish into half
        // circles — the slot, drawn as four lines and four quarter arcs.
        let l = |a: [f64; 2], b: [f64; 2]| Edge2d::line(pt(a[0], a[1]), pt(b[0], b[1])).unwrap();
        let q = |c: [f64; 2], s: [f64; 2]| {
            Edge2d::arc_turns(pt(c[0], c[1]), pt(s[0], s[1]), 1).unwrap()
        };
        let e = vec![
            l([5.0, 0.0], [25.0, 0.0]),
            q([25.0, 5.0], [25.0, 0.0]), // (25,0) → (30,5)
            q([25.0, 5.0], [30.0, 5.0]), // (30,5) → (25,10)
            l([25.0, 10.0], [5.0, 10.0]),
            q([5.0, 5.0], [5.0, 10.0]), // (5,10) → (0,5)
            q([5.0, 5.0], [0.0, 5.0]),  // (0,5) → (5,0)
        ];
        let p = from_edges(e).unwrap();
        assert_eq!(
            p[0].outer().len(),
            4,
            "line, half circle, line, half circle"
        );
        assert_eq!(
            p[0].outer()
                .segs()
                .iter()
                .filter(|s| matches!(s, Seg2d::Arc { .. }))
                .count(),
            2
        );
    }

    #[test]
    fn arcs_that_cannot_be_stated_are_refused_by_name() {
        assert!(matches!(
            Edge2d::arc_turns(pt(1.0, 1.0), pt(0.0, 0.0), 1),
            Err(SketchError::ArcRadiusNotRational { .. })
        )); // r² = 2
        assert!(matches!(
            Edge2d::arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 0),
            Err(SketchError::ArcTurnsOutOfRange { turns: 0 })
        ));
        assert!(matches!(
            Edge2d::arc_turns(pt(0.0, 0.0), pt(5.0, 0.0), 4),
            Err(SketchError::ArcTurnsOutOfRange { turns: 4 })
        ));
        assert!(matches!(
            Edge2d::circle(pt(0.0, 0.0), 0.0),
            Err(SketchError::NonPositiveRadius { .. })
        ));
        let r = |n: i128| Rat::from_int(n);
        assert!(matches!(
            Edge2d::arc_rat([r(0), r(0)], [r(5), r(0)], [r(0), r(4)], true),
            Err(SketchError::ArcEndOffCircle { .. })
        ));
        assert!(matches!(
            Edge2d::arc_rat([r(0), r(0)], [r(5), r(0)], [r(5), r(0)], true),
            Err(SketchError::ZeroLengthArc { .. })
        ));
        // 3-4-5: a rational radius off the axes is fine.
        assert!(Edge2d::arc_rat([r(0), r(0)], [r(3), r(4)], [r(-4), r(3)], true).is_ok());
    }

    #[test]
    fn a_lens_of_two_arcs_is_a_valid_region() {
        // Two arcs between (0,0) and (6,0), centres (3,4) and (3,−4), r 5 — a lens.
        let r = |n: i128| Rat::from_int(n);
        let e = vec![
            Edge2d::arc_rat([r(3), r(4)], [r(0), r(0)], [r(6), r(0)], true).unwrap(),
            Edge2d::arc_rat([r(3), r(-4)], [r(6), r(0)], [r(0), r(0)], true).unwrap(),
        ];
        let p = from_edges(e).unwrap();
        assert_eq!(p[0].outer().len(), 2);
        p[0].check().unwrap();
    }

    #[test]
    fn hole_containment_reads_arcs_on_either_side() {
        let sq = |a: f64, b: f64| -> Vec<Edge2d> {
            vec![
                Edge2d::line(pt(a, a), pt(b, a)).unwrap(),
                Edge2d::line(pt(b, a), pt(b, b)).unwrap(),
                Edge2d::line(pt(b, b), pt(a, b)).unwrap(),
                Edge2d::line(pt(a, b), pt(a, a)).unwrap(),
            ]
        };
        // A round hole in a square plate.
        let mut e = sq(0.0, 10.0);
        e.push(Edge2d::circle(pt(5.0, 5.0), 2.0).unwrap());
        let p = from_edges(e).unwrap();
        assert_eq!((p.len(), p[0].holes().len()), (1, 1));
        assert!(!p[0].holes()[0].is_polygon());
        // A square hole in a round disk.
        let mut e = sq(-2.0, 2.0);
        e.push(Edge2d::circle(pt(0.0, 0.0), 10.0).unwrap());
        let p = from_edges(e).unwrap();
        assert_eq!((p.len(), p[0].holes().len()), (1, 1));
        assert!(!p[0].outer().is_polygon());
        assert!(p[0].holes()[0].is_polygon());
        // A circle crossing the square: the rings meet.
        let mut e = sq(0.0, 10.0);
        e.push(Edge2d::circle(pt(10.0, 5.0), 2.0).unwrap());
        assert!(matches!(from_edges(e), Err(SketchError::RingsMeet { .. })));
    }
}
