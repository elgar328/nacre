use super::*;
/// One closed ring of a [`Profile2d`], stored as its rational truth: its vertices and the step
/// leaving each — straight to the next vertex, or an arc around a stated centre ([`Edge2d`],
/// `nacre-geom`'s vocabulary, which is also what its predicates read). `edges[i]` runs
/// `vertices[i] → vertices[(i + 1) % n]`. Built by `Ring2d::new` (a stated ring, checked), by
/// `Ring2d::circle`, or by the polygon doors; the fields are private so nothing bypasses them.
///
/// The coordinates are what the author's decimals *spelled* (`Rat::from_decimal`), not the f64s
/// that carried them — the same truth/cache split every dimension in the kernel gets.
/// A computed coordinate (an arc's far end) is computed in `Rat` and
/// never rounds through f64.
///
/// **Normal form** ([`Ring2d::normalized`]). A vertex strictly mid-run on a straight edge is
/// dissolved — its two walls would be one plane, so the corner it names has no three-plane
/// definition, and deleting it changes no geometry. Two adjacent arcs of one circle in one
/// direction are one arc — their meeting vertex names no corner either: a pair of fillets that ate
/// their whole edge is a half circle, two half circles are the circle. One vertex with one arc is a
/// whole circle, and that vertex is its seam.
#[derive(Clone, Debug, PartialEq)]
pub struct Ring2d {
    vertices: Vec<[Rat; 2]>,
    edges: Vec<Edge2d>,
}

impl Ring2d {
    /// A polygon ring, in normal form.
    pub(crate) fn polygon(points: Vec<[Rat; 2]>) -> Ring2d {
        let n = points.len();
        Ring2d::normalized(points, vec![Edge2d::Line; n])
    }

    /// The normal form of `edges[i]: vertices[i] → vertices[i + 1]` — see the type doc. Runs to a
    /// fixpoint; each pass removes at most one vertex, and a ring of one or two vertices is left
    /// for [`Profile2d::check`] to judge (one vertex with an arc is a whole circle).
    pub(crate) fn normalized(mut vertices: Vec<[Rat; 2]>, mut edges: Vec<Edge2d>) -> Ring2d {
        debug_assert_eq!(vertices.len(), edges.len(), "one step leaves each vertex");
        loop {
            let n = vertices.len();
            if n < 2 {
                break;
            }
            let mut dropped = false;
            // Later vertices first, the wrap-around pair last, so the ring keeps its first vertex
            // whenever it can (two half circles merge into the circle seamed at the first).
            for k in (1..n).chain(std::iter::once(0)) {
                let prev = (k + n - 1) % n;
                let next = (k + 1) % n;
                let flat = match (edges[prev], edges[k]) {
                    (Edge2d::Line, Edge2d::Line) => {
                        n >= 3 && flat_corner(vertices[prev], vertices[k], vertices[next])
                    }
                    (
                        Edge2d::Arc {
                            center: c1,
                            r2: r1,
                            ccw: w1,
                        },
                        Edge2d::Arc {
                            center: c2,
                            r2,
                            ccw: w2,
                        },
                    ) => c1 == c2 && r1 == r2 && w1 == w2,
                    _ => false,
                };
                if flat {
                    // The step arriving at `k` runs on to `next`; the vertex and the step that
                    // left it go.
                    vertices.remove(k);
                    edges.remove(k);
                    dropped = true;
                    break;
                }
            }
            if !dropped {
                break;
            }
        }
        Ring2d { vertices, edges }
    }

    pub fn vertices(&self) -> &[[Rat; 2]] {
        &self.vertices
    }

    pub fn edges(&self) -> &[Edge2d] {
        &self.edges
    }

    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Straight steps only.
    pub fn is_polygon(&self) -> bool {
        self.edges.iter().all(|s| matches!(s, Edge2d::Line))
    }

    /// **Which way the ring runs**, exactly: `Positive` is counter-clockwise (`+x → +y`) — the
    /// signed-area question over straight steps and quarter-turn arcs, answered in integers by
    /// [`nacre_exact::winding_sign_quarter_arcs`] on the profile's own coordinates (the author's
    /// decimals, not a lifted frame's widths). `None` for an arc that is not a quarter-turn
    /// multiple. Reachable by type, unreached by any producer today: `Edge2d::arc_rat` accepts any
    /// angle, but its one production caller (the kit's fillet) only rounds axis-aligned corners,
    /// so every arc a profile brings here is a quarter-turn multiple. A `None` reaches
    /// `prism_rings_in` and is refused as `PlaneWithoutExactForm` — a wall with no one at it.
    pub(crate) fn winding_sign(&self) -> Option<nacre_exact::Orient> {
        let n = self.vertices.len();
        let (mut lines, mut arcs) = (Vec::new(), Vec::new());
        for i in 0..n {
            let (s0, e0) = (self.vertices[i], self.vertices[(i + 1) % n]);
            match self.edges[i] {
                Edge2d::Line => lines.push([s0, e0]),
                Edge2d::Arc { center, r2, ccw } => arcs.push(nacre_exact::QuarterArc {
                    center,
                    r2,
                    start: s0,
                    end: e0,
                    ccw,
                }),
            }
        }
        nacre_exact::winding_sign_quarter_arcs(&lines, &arcs)
    }
}

/// A closed planar region: one outer ring and any number of hole rings — straight steps and
/// circular arcs ([`Ring2d`]).
///
/// **Winding is not the caller's business, and not this type's either.** The rings are stored in
/// author order; the prism builder is the single place that decides orientation, because it is
/// the only one that knows the sweep direction (a pocket sweeps *into* a face, which flips what
/// "counter-clockwise" means). Fixing the winding here as well would put that decision in two
/// places, which is exactly how the holes and the outer ring come to disagree.
///
/// **The constructors take the boundary f64s and keep the truth**:
/// each coordinate becomes the rational its shortest decimal spells, a dimension outside the
/// decimal window (`~1e38` above, `~1e-22` below for a 17-digit value) is a named error
/// **at construction** rather than a silent f64 fallback downstream, and every **flat corner**
/// (a vertex strictly mid-run on a straight edge) is dissolved — its two walls would be one
/// plane, so the corner it names has no three-plane definition, and deleting it changes no
/// geometry. What survives construction is the profile's normal form.
///
/// The constructors do **not** run [`Profile2d::check`]; that contract (simplicity, disjoint
/// rings, hole containment) is `O(n²)` and **every operation that consumes a profile runs it
/// first**, so the kernel never works from an unverified one. `sketch::from_rings` and
/// `sketch::from_paths` are the checking constructors — the sketch doors.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile2d {
    outer: Ring2d,
    holes: Vec<Ring2d>,
}

/// Which ring of a [`Profile2d`] an error is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileRing {
    /// The outer ring.
    Outer,
    /// The hole at this index in [`Profile2d::holes`].
    Hole(usize),
}

impl Profile2d {
    /// A simple polygon — no holes. `Err` when a coordinate has no rational truth (outside the
    /// decimal window); see the type doc.
    pub fn polygon(points: Vec<Point2>) -> Result<Profile2d, OpError> {
        Profile2d::with_holes(points, Vec::new())
    }

    /// An outer ring with holes. Lifts every coordinate to its rational truth and dissolves flat
    /// corners; the *contract* ([`Profile2d::check`]) is still enforced by every operation that
    /// consumes the profile.
    pub fn with_holes(outer: Vec<Point2>, inners: Vec<Vec<Point2>>) -> Result<Profile2d, OpError> {
        let lift = |ring: &[Point2], id: ProfileRing| -> Result<Vec<[Rat; 2]>, OpError> {
            ring.iter()
                .enumerate()
                .map(|(i, p)| {
                    let point = |x: f64| Rat::from_decimal(x);
                    match (point(p[0]), point(p[1])) {
                        (Some(x), Some(y)) => Ok([x, y]),
                        _ => Err(OpError::ProfileOutsideDecimalWindow { ring: id, point: i }),
                    }
                })
                .collect()
        };
        let outer = lift(&outer, ProfileRing::Outer)?;
        let holes = inners
            .iter()
            .enumerate()
            .map(|(h, ring)| lift(ring, ProfileRing::Hole(h)))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Profile2d::from_rat_rings(outer, holes))
    }

    /// The constructor behind the f64 boundary — rings already in their rational truth
    /// (`sketch::from_rings` arrives here after lifting and classifying on that truth).
    /// Normalization (the flat-corner dissolve) happens here, once, for every entry path.
    pub(crate) fn from_rat_rings(outer: Vec<[Rat; 2]>, holes: Vec<Vec<[Rat; 2]>>) -> Profile2d {
        Profile2d {
            outer: Ring2d::polygon(outer),
            holes: holes.into_iter().map(Ring2d::polygon).collect(),
        }
    }

    /// Rings already in normal form — `sketch`'s door once it has chained and sorted them.
    pub(crate) fn from_normalized_rings(outer: Ring2d, holes: Vec<Ring2d>) -> Profile2d {
        Profile2d { outer, holes }
    }

    /// Whether any ring carries an arc.
    pub fn has_arcs(&self) -> bool {
        !self.outer.is_polygon() || self.holes.iter().any(|h| !h.is_polygon())
    }

    /// Verify the contract: every ring is a **simple polygon** of at least three points, the rings
    /// are pairwise disjoint, and each hole lies inside the outer ring and inside no other hole.
    ///
    /// Every operation that consumes a profile calls this first, so an unverified profile never
    /// reaches the topology. It is public because an app can ask before it builds.
    ///
    /// **Simplicity is a contract on authored input, not an invariant of kernel data.** A drawn
    /// ring's inside is *defined* by even-odd, which needs simplicity to mean anything; the
    /// contours a boolean produces get their inside from the arrangement instead, and those may
    /// legitimately be non-simple (a figure-8 pinch — see `loop_winding`). Do not "unify" the two.
    ///
    /// Exact — on the truth: every decision is an `orient2d` sign on the rationals the author's
    /// decimals spelled, not on their f64 carriers. The two can disagree (three points collinear
    /// in decimal sit a hair off the line in binary), and it is the decimal that is the author's
    /// meaning. Edge indices in the errors refer to the normalized rings [`outer`](Profile2d::outer)
    /// returns — construction dissolved flat corners first.
    ///
    /// `O(n²)` in the ring size, and it runs on every extrude and every replay of one. Measured
    /// on a convex ring of 17-digit coordinates (the worst case — nothing short-circuits;
    /// `measure_profile_check_wall_clock` in `tests/perf.rs`): **34 ms at 100 points, 3.2 s at 1
    /// 000** —
    /// each rational sign costs ~1.7 µs (the narrow route's gcd reductions; such coordinates stay
    /// inside `i128`) against an f64 predicate's nanoseconds. Hand-drawn sketches are tens
    /// of points and pay well under a millisecond; a generator emitting thousands of points per
    /// ring is firmly outside this function's comfort, and the named follow-ups — an exact
    /// bounding-box prefilter for the segment pairs, or an f64 filter with a sound error bound
    /// escalating to `Rat` — are deliberately
    /// unbuilt until that population exists.
    pub fn check(&self) -> Result<(), OpError> {
        let rings = || {
            std::iter::once((ProfileRing::Outer, &self.outer)).chain(
                self.holes
                    .iter()
                    .enumerate()
                    .map(|(i, r)| (ProfileRing::Hole(i), r)),
            )
        };
        let undecidable = |_| OpError::ProfileUndecidable;
        for (id, r) in rings() {
            if r.is_empty() || (r.is_polygon() && r.len() < 3) {
                return Err(OpError::DegenerateProfile);
            }
            // Simplicity first: the parity below is only meaningful on a simple ring. The
            // predicate reports a zero-length edge as a pair with itself.
            if let Some((a, b)) = mixed_ring_self_intersection(r.mixed()).map_err(undecidable)? {
                return Err(if a == b {
                    OpError::ZeroLengthProfileEdge { ring: id, edge: a }
                } else {
                    OpError::SelfIntersectingProfile {
                        ring: id,
                        edges: (a, b),
                    }
                });
            }
        }
        let all: Vec<(ProfileRing, &Ring2d)> = rings().collect();
        for (i, (ida, a)) in all.iter().enumerate() {
            for (idb, b) in &all[i + 1..] {
                if mixed_rings_cross(a.mixed(), b.mixed()).map_err(undecidable)? {
                    return Err(OpError::ProfileRingsMeet { a: *ida, b: *idb });
                }
            }
        }
        for (h, hole) in self.holes.iter().enumerate() {
            let probe = hole.vertices()[0];
            if point_in_mixed_ring(probe, self.outer.mixed()).map_err(undecidable)?
                != RingSide::Inside
            {
                return Err(OpError::HoleNotInsideOuter { hole: h });
            }
            for (k, other) in self.holes.iter().enumerate() {
                if k != h
                    && point_in_mixed_ring(probe, other.mixed()).map_err(undecidable)?
                        == RingSide::Inside
                {
                    return Err(OpError::NestedHole { outer: k, inner: h });
                }
            }
        }
        Ok(())
    }

    /// The outer ring — the normalized rational truth.
    pub fn outer(&self) -> &Ring2d {
        &self.outer
    }

    /// The hole rings — the normalized rational truth.
    pub fn holes(&self) -> &[Ring2d] {
        &self.holes
    }
}
