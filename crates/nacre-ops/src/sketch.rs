//! Rings into profiles — the nesting policy on top of geom's exact predicates.
//!
//! A sketch is drawn as closed paths; which of them are material, which are holes, and which are
//! islands inside holes is **not** declared by the author but decided from containment, the way
//! every CAD sketcher does it. The classification itself (is this point inside that ring? do two
//! rings meet?) is a robustness-sensitive sign question and lives in `nacre-geom` (the design
//! isolates that — [`nacre_geom::intersect`] for polygons, [`nacre_geom::mixed`] once arcs are
//! in); what is here is the policy — depth parity — and the doors.
//!
//! **The door is the ring, not the pen.** The kernel takes what every way of drawing ends in — an
//! ordered ring of exact vertices with the step leaving each ([`Ring2d::new`]) — and knows nothing
//! about how it was drawn. A pen, a constraint solver, an imported file all arrive here the same
//! way. There is no chaining of loose edges: order is the caller's, stated once, and an open or
//! branching outline is not a thing this door can be handed.
//!
//! **Edges are stated, not derived.** An arc is its centre, its radius and its two vertices, all
//! rational, and the step doors here are the only place a coordinate is computed — the end of a
//! quarter-turn arc is the start rotated in `Rat`, never an f64 that was rounded on the way. The
//! f64 doors take what the author *wrote* (their decimals become the rationals they spell); the
//! `Rat` doors take what a caller computed exactly.

use crate::{Profile2d, Ring2d};
use nacre_geom::intersect::RingSide;
use nacre_geom::mixed::{
    Edge2d, MixedRing, Undecidable, mixed_ring_self_intersection, mixed_rings_cross,
    point_in_mixed_ring,
};
use nacre_math::Point2;
use nacre_scalar::Rat;

/// Why a ring, or a set of rings, is not a valid set of profiles.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum SketchError {
    /// A ring with fewer than three points encloses nothing.
    DegenerateRing { ring: usize },
    /// Two rings touch or cross. A hole must lie strictly inside its outer ring and strictly
    /// outside its siblings; anything else has no unambiguous inside.
    RingsMeet { a: usize, b: usize },
    /// A ring's vertices and edges do not pair up — one edge leaves each vertex.
    UnevenRing { vertices: usize, edges: usize },
    /// A straight step whose two vertices are one point. `edge` is its index in the ring.
    ZeroLengthEdge { edge: usize },
    /// A ring meets itself — it crosses, touches, or doubles back over its own edge. Such a ring
    /// has no unambiguous inside, so it cannot be sorted into a profile at all.
    ///
    /// Reported as the two offending edges' **chord midpoints** — points are what the rest of
    /// this enum reports too, and a point is what an editor can put a marker on.
    RingSelfIntersects { ring: usize, at: [[f64; 2]; 2] },
    /// A coordinate outside the decimal window (`~1e38` above, `~1e-22` below for a full-width
    /// value) has no rational truth for the kernel to keep — the sketch-layer twin of
    /// `OpError::ProfileOutsideDecimalWindow`. `at` is the offending point.
    OutsideDecimalWindow { at: [f64; 2] },
    /// An arc's end is not on the circle its centre and start define.
    ArcEndOffCircle {
        center: [f64; 2],
        start: [f64; 2],
        end: [f64; 2],
    },
    /// An arc step was handed to [`Ring2d::new`] with an `r2` that is not the squared distance
    /// from its centre to the vertex it leaves — the stated step and the ring disagree about the
    /// circle. `stated` is the step's radius, `√r2`, as a number.
    ArcRadiusMismatch {
        center: [f64; 2],
        start: [f64; 2],
        stated: f64,
    },
    /// An arc from a point back to itself is a whole circle, which [`Ring2d::circle`] states;
    /// this door asked for a proper arc.
    ZeroLengthArc { at: [f64; 2] },
    /// A quarter-turn count outside `±1..=±3` — `0` is no arc, `±4` is a whole circle
    /// ([`Ring2d::circle`]).
    ArcTurnsOutOfRange { turns: i32 },
    /// A circle or arc with a radius that is not positive.
    NonPositiveRadius { center: [f64; 2], radius: f64 },
    /// Checked `Rat` arithmetic overflowed while classifying: the question has an answer the
    /// kernel cannot state here. Refused by name; nothing guesses.
    Undecidable,
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

/// The squared radius `|start − center|²` — rational for any two rational points, so no circle
/// is refused for its radius being irrational — or the named refusal for a zero one.
fn r2_of(center: [Rat; 2], start: [Rat; 2]) -> Result<Rat, SketchError> {
    let r2 = dist2(start, center).ok_or(SketchError::Undecidable)?;
    if r2 <= Rat::from_int(0) {
        return Err(SketchError::NonPositiveRadius {
            center: f2(center),
            radius: 0.0,
        });
    }
    Ok(r2)
}

// ---- step doors: an arc as the step it is, with the vertex it ends on ----

/// The arc from `start` around `center` through `quarter_turns` right angles — positive is
/// counter-clockwise (`+x → +y` in the sketch's own frame), `±1..=±3` — as the step it is and
/// **the vertex it ends on**. The end is the start **rotated in `Rat`**: it lies on the circle by
/// construction, and no computed value passes through f64. Today's whole arc vocabulary is this
/// door, [`arc_to_rat`] and [`Ring2d::circle`]; an arbitrary angle is a point the kernel would
/// have to *name* (a start turned by θ), which is a later vocabulary, not a rounding.
pub fn arc_turns(
    center: Point2,
    start: Point2,
    quarter_turns: i32,
) -> Result<(Edge2d, [Rat; 2]), SketchError> {
    arc_turns_rat(lift(center)?, lift(start)?, quarter_turns)
}

/// [`arc_turns`] for computed points: the centre and the start as `Rat` — the door a pen takes
/// when it already stands on an exact point (after a fillet's retreat, after an arc).
pub fn arc_turns_rat(
    center: [Rat; 2],
    start: [Rat; 2],
    quarter_turns: i32,
) -> Result<(Edge2d, [Rat; 2]), SketchError> {
    if !(1..=3).contains(&quarter_turns.unsigned_abs()) {
        return Err(SketchError::ArcTurnsOutOfRange {
            turns: quarter_turns,
        });
    }
    let r2 = r2_of(center, start)?;
    let v = sub2(start, center).ok_or(SketchError::Undecidable)?;
    let end = turned(center, v, quarter_turns).ok_or(SketchError::Undecidable)?;
    Ok((
        Edge2d::Arc {
            center,
            r2,
            ccw: quarter_turns > 0,
        },
        end,
    ))
}

/// A proper arc from exact data — the door for computed coordinates (a fillet's tangent points
/// and centre, say). Checks what [`arc_turns`] guarantees by construction: `end` on the circle,
/// `start ≠ end`. The step it returns runs `start → end` in a ring that states those two
/// vertices.
pub fn arc_to_rat(
    center: [Rat; 2],
    start: [Rat; 2],
    end: [Rat; 2],
    ccw: bool,
) -> Result<Edge2d, SketchError> {
    if start == end {
        return Err(SketchError::ZeroLengthArc { at: f2(start) });
    }
    let r2 = r2_of(center, start)?;
    if dist2(end, center).ok_or(SketchError::Undecidable)? != r2 {
        return Err(SketchError::ArcEndOffCircle {
            center: f2(center),
            start: f2(start),
            end: f2(end),
        });
    }
    Ok(Edge2d::Arc { center, r2, ccw })
}

impl Ring2d {
    /// A ring from its vertices and the step leaving each — `edges[i]` runs
    /// `vertices[i] → vertices[(i + 1) % n]` — checked and put in normal form.
    ///
    /// What is checked is what a stated step can get wrong: the two lists pair up
    /// ([`SketchError::UnevenRing`]); a straight step does not start where it ends
    /// ([`SketchError::ZeroLengthEdge`]); an arc's stated `r2` is the squared distance from its
    /// centre to the vertex it leaves ([`SketchError::ArcRadiusMismatch`]) and the vertex it
    /// arrives on lies on that circle ([`SketchError::ArcEndOffCircle`]); an arc between one and
    /// the same point is the whole circle only when it is the ring's sole step. Whether the ring
    /// encloses anything, meets itself, or meets another is [`from_paths`]'s question.
    pub fn new(vertices: Vec<[Rat; 2]>, edges: Vec<Edge2d>) -> Result<Ring2d, SketchError> {
        if vertices.len() != edges.len() {
            return Err(SketchError::UnevenRing {
                vertices: vertices.len(),
                edges: edges.len(),
            });
        }
        let n = vertices.len();
        for (i, edge) in edges.iter().enumerate() {
            let (a, b) = (vertices[i], vertices[(i + 1) % n]);
            match *edge {
                Edge2d::Line => {
                    if a == b {
                        return Err(SketchError::ZeroLengthEdge { edge: i });
                    }
                }
                Edge2d::Arc { center, r2, .. } => {
                    if r2 <= Rat::from_int(0) {
                        return Err(SketchError::NonPositiveRadius {
                            center: f2(center),
                            radius: r2.to_f64().signum() * r2.to_f64().abs().sqrt(),
                        });
                    }
                    if r2_of(center, a)? != r2 {
                        return Err(SketchError::ArcRadiusMismatch {
                            center: f2(center),
                            start: f2(a),
                            stated: r2.to_f64().sqrt(),
                        });
                    }
                    if a == b && n > 1 {
                        return Err(SketchError::ZeroLengthArc { at: f2(a) });
                    }
                    if dist2(b, center).ok_or(SketchError::Undecidable)? != r2 {
                        return Err(SketchError::ArcEndOffCircle {
                            center: f2(center),
                            start: f2(a),
                            end: f2(b),
                        });
                    }
                }
            }
        }
        Ok(Ring2d::normalized(vertices, edges))
    }

    /// A polygon ring from the points the author wrote — the per-ring half of [`from_rings`]
    /// (their decimals become the rationals they spell; a repeated point is left for
    /// [`Profile2d::check`] to name, as `from_rings` leaves it).
    pub fn polygon_decimal(points: Vec<Point2>) -> Result<Ring2d, SketchError> {
        let lifted = points
            .iter()
            .map(|p| lift(*p))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Ring2d::polygon(lifted))
    }

    /// The whole circle of `radius` about `center`, counter-clockwise, as the ring of one vertex
    /// it is — its seam at `center + (r, 0)`, the same point the cylinder primitive seams at
    /// (`+ref_dir`), which is what lets the two roads state one solid.
    pub fn circle(center: Point2, radius: f64) -> Result<Ring2d, SketchError> {
        let c = lift(center)?;
        let r = Rat::from_decimal(radius)
            .ok_or(SketchError::OutsideDecimalWindow { at: [radius, 0.0] })?;
        Ring2d::circle_rat(c, r)
    }

    /// [`Ring2d::circle`] for a computed radius — a diameter halved in `Rat`, say — and a
    /// computed centre.
    pub fn circle_rat(center: [Rat; 2], radius: Rat) -> Result<Ring2d, SketchError> {
        if radius <= Rat::from_int(0) {
            return Err(SketchError::NonPositiveRadius {
                center: f2(center),
                radius: radius.to_f64(),
            });
        }
        let seam = [
            center[0]
                .checked_add(radius)
                .ok_or(SketchError::Undecidable)?,
            center[1],
        ];
        Ok(Ring2d::normalized(
            vec![seam],
            vec![Edge2d::Arc {
                center,
                r2: radius.checked_mul(radius).ok_or(SketchError::Undecidable)?,
                ccw: true,
            }],
        ))
    }

    /// A ring's view for the predicates.
    pub(crate) fn mixed(&self) -> MixedRing<'_> {
        MixedRing {
            vertices: self.vertices(),
            segs: self.edges(),
        }
    }

    /// The f64 midpoint of step `i`'s chord — diagnostic material for a message.
    fn chord_midpoint(&self, i: usize) -> [f64; 2] {
        let n = self.len();
        let (p, q) = (f2(self.vertices()[i]), f2(self.vertices()[(i + 1) % n]));
        [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0]
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

/// Sort rings — straight or arc-bearing, each already checked and in normal form
/// ([`Ring2d::new`], [`Ring2d::circle`], [`Ring2d::polygon_decimal`]) — into profiles by
/// containment depth: [`from_rings`]'s policy for every kind of ring. Ring indices in the errors
/// are positions in `rings`.
pub fn from_paths(rings: Vec<Ring2d>) -> Result<Vec<Profile2d>, SketchError> {
    classify(rings)
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

#[cfg(test)]
#[path = "tests/sketch.rs"]
mod tests;
