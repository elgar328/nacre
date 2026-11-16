use super::*;
/// **A face corner, in whichever spelling it has** — the footprint test's currency.
///
/// Three plane-carried corners have rational coordinates; a corner a cylinder made does not, and
/// is `line.base() + s·line.dir()` with `s` quadratic-irrational; a corner no exact spelling
/// reaches is the box its realization proves ([`Corner::Bounded`]). Every spelling answers the
/// same three questions, so they are asked through one type rather than branched at each call.
pub(super) enum Corner {
    Rational(nacre_exact::MeetPoint),
    Pierce(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    /// **A piece of a circle in the face's own plane** — the whole disk a single-circle loop *is*
    /// (`theta: Whole`), or one **arc** of a loop that mixes lines and arcs (`theta: Arc`).
    ///
    /// The name «corner» is narrower than the type: what it carries is a **boundary piece's
    /// reach**, and a vertex is the piece with no width. Written this way the loops below do not
    /// learn a new shape; they ask the same three questions and one of the answers has extent.
    ///
    /// `axis` is the **carrier cylinder's** axis, and the arc's `from → to` is counter-clockwise
    /// about it ([`lateral_theta_extent`] states that convention). The face's own normal may run
    /// the other way; using it would name the **complementary** arc, which is not a bound on this
    /// piece but a different set.
    Round {
        centre: [nacre_exact::Rat; 3],
        /// The carrier's squared radius.
        rho2: nacre_exact::BigRat,
        axis: [nacre_exact::Rat; 3],
        theta: Theta,
    },
    /// **A corner known only as the box its realization proves** ([`rational_reach`]) — what a
    /// vertex is when no exact spelling reaches it: a three-plane corner whose planes share no frame
    /// (a wall turned off the quarter turns among them), a meet `Rat` does not hold, a chain that
    /// does not carry it, a seam corner, a pierce corner whose meet cannot be stated.
    ///
    /// The box is held as a ball: `centre` is its middle and `rho` the sum of its half-widths, which
    /// is at least their Euclidean length, so the ball holds the box. The true point lies on the
    /// face's plane, so it is no further from `centre`'s projection onto that plane than from
    /// `centre`: on the plane the corner is a disk of radius `rho`.
    ///
    /// **A box proves apart and nothing else.** It can show the corner on one side of a strip or
    /// past a station, never on a line or in a plane — so [`Corner::strip_side`] has no side for it
    /// where the disk reaches the strip, and [`Corner::on_plane`] asks the corner's planes rather
    /// than the box: `carriers` are the world names of its definition's planar carriers (`None` for
    /// a cylinder, or a plane with no narrow world name), and the corner lies on a class's plane
    /// exactly when one of them is that plane's name.
    ///
    /// `cause` is why no exact spelling reached the corner — the name a reader that needs the
    /// corner's value, or a side the box cannot prove, refuses under ([`Corner::refusal`]).
    Bounded {
        centre: [nacre_exact::Rat; 3],
        rho: nacre_exact::Rat,
        carriers: [Option<[nacre_exact::Rat; 4]>; 3],
        cause: CornerFail,
    },
}

impl Corner {
    /// Does this corner lie on the plane the class names? — [`nacre_exact::cylinder_strip_side`]'s precondition.
    pub(super) fn on_plane(&self, coeffs: &[nacre_exact::Rat; 4]) -> bool {
        match self {
            Self::Rational(p) => nacre_exact::point_on_plane_exact(coeffs, p),
            Self::Pierce(line, s) => {
                nacre_exact::quad::plane_side(coeffs, line, s) == nacre_exact::Orient::Zero
            }
            Self::Round { centre, .. } => {
                nacre_exact::point_on_plane_exact(coeffs, &nacre_exact::MeetPoint::Narrow(*centre))
            }
            Self::Bounded { carriers, .. } => carriers.iter().flatten().any(|c| c == coeffs),
        }
    }

    /// **Where this piece stands relative to the strip, where it can prove it** —
    /// [`nacre_exact::cylinder_strip_side`]'s answer in the piece's own spelling, and `None` for a
    /// side the piece cannot prove. The two readers fold `None` each to its own safe reading: the
    /// footprint test as «not clear», the tangency road as «cannot read this loop» — an exact
    /// corner's [`nacre_exact::StripSide::Inside`] is a corner **on** the line there, which an
    /// unproven side is not.
    pub(super) fn strip_side(
        &self,
        coeffs: &[nacre_exact::Rat; 4],
        o: &[nacre_exact::Rat; 3],
        m: &[nacre_exact::Rat; 3],
        r2: &nacre_exact::BigRat,
    ) -> Option<nacre_exact::StripSide> {
        use nacre_exact::StripSide;
        match self {
            Self::Rational(p) => Some(nacre_exact::cylinder_strip_side(coeffs, p, o, m, r2)),
            Self::Pierce(line, s) => Some(nacre_exact::cylinder_strip_side_branch(
                coeffs, line, s, o, m, r2,
            )),
            // The strip runs along `e = n × m`, so **that** is the direction the piece's extent
            // is asked for. A whole circle reaches `±ρ|e|` alike and takes the symmetric door,
            // which is its complete answer; an arc's two ends differ and take the
            // general one. `None` from the extent is an overflow.
            Self::Round { .. } => round_strip_side(self, coeffs, o, m, r2),
            // The disk about the centre's projection — the foot of the perpendicular from it, the
            // plane's meet with the line along its normal: only a disk clear of the strip proves a
            // side. Reaching the strip and spanning a ruling are both things the true point may or
            // may not do, so neither is claimed.
            Self::Bounded { centre, rho, .. } => {
                let n = [coeffs[0], coeffs[1], coeffs[2]];
                let p = nacre_exact::axis_plane_meet(coeffs, centre, &n)?;
                match nacre_exact::cylinder_strip_side_margin(
                    coeffs,
                    &nacre_exact::MeetPoint::Narrow(p),
                    &nacre_exact::BigRat::square_of(*rho),
                    o,
                    m,
                    r2,
                ) {
                    s @ (StripSide::Plus | StripSide::Minus) => Some(s),
                    StripSide::Inside | StripSide::Crosses => None,
                }
            }
        }
    }

    /// **The name a reader refuses this piece under where it needs what the piece cannot give** — a
    /// value on the line, or a side it could not prove: a box's own cause
    /// ([`CornerFail::reason_on`]), and for an exact spelling the name the face's boundary gives
    /// ([`unreadable_reason`]).
    pub(super) fn refusal(&self, model: &Model, face: &Face) -> RejectReason {
        match self {
            Self::Bounded { cause, .. } => cause.reason_on(model, face),
            Self::Rational(_) | Self::Pierce(..) | Self::Round { .. } => {
                unreadable_reason(model, face)
            }
        }
    }

    /// **Does this piece reach strictly past the plane at axis parameter `t`, to the `want` side?**
    ///
    /// **A point answers this with one sign; a piece with width cannot.** The caller asks
    /// twice — «did anything go above the span's start» and «did anything go below its end» — and
    /// for a vertex those are the two readings of a single `Orient`. A disk straddling a station
    /// is above it *and* below it, which no single `Orient` can say: answering `Zero` would make
    /// **both** questions false and let a face that covers the station be read as clearing it.
    /// So the question is named instead of inferred.
    ///
    /// **The span is read open** (`face_clears_footprint`'s own note: a corner sitting exactly on a
    /// cap plane must not be dropped), which is why the disk's comparison is spelled here rather
    /// than borrowed from `cyl_gate::reach_clears` — that one is closed against closed, and folding
    /// the two together would move one convention silently.
    fn reaches(
        &self,
        o: &[nacre_exact::Rat; 3],
        m: &[nacre_exact::Rat; 3],
        t: nacre_exact::Rat,
        want: nacre_exact::Orient,
    ) -> bool {
        match self {
            Self::Rational(p) => nacre_exact::point_axis_side(p, o, m, t) == want,
            Self::Pierce(line, s) => nacre_exact::point_axis_side_branch(line, s, o, m, t) == want,
            // The ball sweeps `q ± ρ|m|` with `q = (c − o)·m − t(m·m)`, so «reaches past» is
            // `q` on the wanted side or `ρ²(m·m) > q²` — squared once, exact in `BigRat`.
            // Overflow answers **`true`**, as for a round piece below.
            Self::Bounded { centre, rho, .. } => {
                let zero = nacre_exact::Rat::from_int(0);
                let q = (|| {
                    let mut rel = [zero; 3];
                    for k in 0..3 {
                        rel[k] = centre[k].checked_sub(o[k])?;
                    }
                    let mm = nacre_exact::dot3_rat(m, m)?;
                    Some((
                        nacre_exact::dot3_rat(&rel, m)?.checked_sub(t.checked_mul(mm)?)?,
                        mm,
                    ))
                })();
                match (q, want) {
                    (_, nacre_exact::Orient::Zero) => false,
                    (None, _) => true,
                    (Some((q, mm)), w) => {
                        let toward = if w == nacre_exact::Orient::Positive {
                            q > zero
                        } else {
                            q < zero
                        };
                        toward
                            || nacre_exact::BigRat::square_of(*rho).mul_rat(mm)
                                > nacre_exact::BigRat::square_of(q)
                    }
                }
            }
            // The comparison quantity is `q = (p − o)·m − t(m·m)`; over the piece it sweeps
            // `q_c + [lo_off − √rho2_lo, hi_off + √rho2_hi]` ([`arc_extent`]), so «reaches past»
            // is that end on the wanted side, or its radical covering the gap —
            // `rho2 > end²`, squared once and rational throughout. A whole circle is the
            // symmetric case (`off = 0`, `rho2 = ρ²|m⊥|²`), which is what this said before an
            // arc could be a piece.
            // Overflow answers **`true`**: "not shown to clear" is this test's safe direction.
            Self::Round {
                centre,
                rho2,
                axis,
                theta,
            } => {
                let zero = nacre_exact::Rat::from_int(0);
                let end = (|| {
                    let mut rel = [zero; 3];
                    for k in 0..3 {
                        rel[k] = centre[k].checked_sub(o[k])?;
                    }
                    let mm = nacre_exact::dot3_rat(m, m)?;
                    let q = nacre_exact::dot3_rat(&rel, m)?.checked_sub(t.checked_mul(mm)?)?;
                    let (lo_off, rho2_lo, hi_off, rho2_hi) = arc_extent(theta, rho2, axis, m)?;
                    Some(match want {
                        nacre_exact::Orient::Positive => (q.checked_add(hi_off)?, rho2_hi),
                        _ => (q.checked_add(lo_off)?, rho2_lo),
                    })
                })();
                match (end, want) {
                    (_, nacre_exact::Orient::Zero) => false,
                    (Some((e, rho2)), nacre_exact::Orient::Positive) => {
                        e > zero || e.checked_mul(e).is_none_or(|sq| rho2 > sq)
                    }
                    (Some((e, rho2)), nacre_exact::Orient::Negative) => {
                        e < zero || e.checked_mul(e).is_none_or(|sq| rho2 > sq)
                    }
                    (None, _) => true,
                }
            }
        }
    }
}

/// **A pierce vertex's exact point, solved from its own definition.**
///
/// The rank is the point's place about the cylinder ([`nacre_exact::quad::seam_ranked`]), which
/// reads only the cylinder's statement — so the planes' order and the way their names point
/// do not enter, and the world statements answer it directly.
///
/// **A failure is never a shape this road cannot spell** — the caller has already established
/// that the vertex is a `Pierce`, and a seam vertex never reaches here at all. It is a missing
/// description, or a meet with no crossing at the name's rank — the name was minted from this
/// same meet, so the stored definition and the re-solved meet disagreeing is a naming defect
/// ([`CornerFail::Inconsistent`]).
fn pierce_corner(
    model: &Model,
    planes: [Handle<Surface>; 2],
    cylinder: Handle<Surface>,
    rank: u8,
) -> Result<Corner, CornerFail> {
    let (line, s) = pierce_on(model, planes, &world_cylinder(model, cylinder)?, rank)?;
    Ok(Corner::Pierce(line, s))
}

/// [`pierce_corner`] on its cylinder's world statement `def`, for a reader that holds it already.
pub(super) fn pierce_on(
    model: &Model,
    planes: [Handle<Surface>; 2],
    def: &nacre_topo::CylinderDef,
    rank: u8,
) -> Result<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal), CornerFail> {
    let p1 = world_narrow(model, planes[0])?;
    let p2 = world_narrow(model, planes[1])?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let meet =
        nacre_exact::quad::plane_plane_cylinder(&p1, &p2, &o, &m, r2).ok_or(CornerFail::Width)?;
    nacre_exact::quad::seam_ranked(meet, &o, &m, &def.ref_dir())
        .ok_or(CornerFail::Width)?
        .into_iter()
        .nth(usize::from(rank))
        .ok_or(CornerFail::Inconsistent)
}

/// Why a boundary piece, or a statement it reads, could not be read — four causes, four names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CornerFail {
    /// A shape this road cannot spell exactly — a seam vertex, a three-plane corner whose planes
    /// share no frame, a circle on a face it is not a disk of, an arc end at an irrational root
    /// (its radial vector is no rational vector). A corner of the first two kinds is still read on
    /// its box ([`start_corner`]); this is what it refuses by where no box stands, or where a reader
    /// needs what a box cannot give. The name depends on the boundary ([`unreadable_reason`]).
    Shape,
    /// A statement that cannot be put in the world in one frame — a plane with no name, a chain
    /// that does not fold. Geometry: [`RejectReason::CylinderGateUndecided`].
    Unstated,
    /// A statement the world holds and `Rat` does not — a `Wide` name (moved or not) or meet, an
    /// overflow carrying or combining. Width: [`RejectReason::WitnessNotRational`].
    Width,
    /// A stored pierce definition its own re-solved meet does not offer as a point — the name was
    /// minted from that meet. A defect: [`RejectReason::CylinderStagesDisagree`].
    Inconsistent,
}

/// The cause a world-statement door names, as this road's cause.
impl From<nacre_topo::WorldFail> for CornerFail {
    fn from(f: nacre_topo::WorldFail) -> Self {
        match f {
            nacre_topo::WorldFail::Unstated => Self::Unstated,
            nacre_topo::WorldFail::Width => Self::Width,
        }
    }
}

impl CornerFail {
    /// The name a stated failure refuses under; `None` for [`Self::Shape`], whose name the caller
    /// reads off the boundary.
    pub(super) fn stated_reason(&self) -> Option<RejectReason> {
        match self {
            Self::Shape => None,
            Self::Unstated => Some(RejectReason::CylinderGateUndecided),
            Self::Width => Some(RejectReason::WitnessNotRational),
            Self::Inconsistent => Some(RejectReason::CylinderStagesDisagree),
        }
    }

    /// The name this failure refuses under on `face`: its stated name, else the one the face's
    /// boundary gives a shape this road cannot spell ([`unreadable_reason`]).
    pub(super) fn reason_on(self, model: &Model, face: &Face) -> RejectReason {
        self.stated_reason()
            .unwrap_or_else(|| unreadable_reason(model, face))
    }
}

/// **A plane's narrow world coefficients, or why it has none** — the narrow projection of
/// [`nacre_topo::Model::world_plane_name`], whose cause it keeps: geometry (`Unstated`) when the
/// truth puts the plane off the rationals in the world, width (`Width`) when the world name is
/// rational and `Rat` cannot hold it — an unmoved `Wide` name here, which the door hands back
/// verbatim.
pub(crate) fn world_narrow(
    model: &Model,
    plane: Handle<Surface>,
) -> Result<[nacre_exact::Rat; 4], CornerFail> {
    model
        .world_plane_name(plane)?
        .narrow()
        .copied()
        .ok_or(CornerFail::Width)
}

/// A cylinder's world statement, or why it has none — the door's own cause.
pub(super) fn world_cylinder(
    model: &Model,
    cylinder: Handle<Surface>,
) -> Result<nacre_topo::CylinderDef, CornerFail> {
    Ok(model.world_cylinder_def(cylinder)?)
}

/// Per axis, the low and high ends of where a vertex's truth can be.
pub(crate) type CacheReach = ([f64; 3], [f64; 3]);

/// **Where a vertex's truth can be** — per axis, the interval its cache proves. The self-touch
/// sieve reads it for seam vertices; the gate reads it in rationals ([`rational_reach`]), for the
/// corners of a face it asks to stand apart and for a corner no exact spelling names.
///
/// A `Realized` cache holds the truth within `coord ± bound`; the radius is widened by a few ulps
/// of both operands (and the smallest normal), since `coord − bound` computed in `f64` may round
/// inward by half an ulp of the result — which the widening outgrows on either side of a binade.
/// Widening only over-keeps. Any other variant proves nothing, and the answer is the whole line:
/// the sieve then keeps every candidate the vertex takes part in and leaves it to the exact test,
/// which is the only safe reading of "no bound". Measured: every seam point in the census and the
/// perf folds is `Realized`; the suite's 72 others (a history past the cost cap, a carrier with no
/// witness) are where this arm runs.
pub(crate) fn cache_reach(cache: &nacre_topo::PointCache) -> CacheReach {
    match *cache {
        nacre_topo::PointCache::Realized { coord, bound } => {
            let c = coord.as_array();
            let w: [f64; 3] = core::array::from_fn(|k| {
                bound[k].upper_f64() * (1.0 + 4.0 * f64::EPSILON) + round_slack(c[k])
            });
            (
                core::array::from_fn(|k| c[k] - w[k]),
                core::array::from_fn(|k| c[k] + w[k]),
            )
        }
        nacre_topo::PointCache::Ceiling { .. } | nacre_topo::PointCache::Unrealized { .. } => {
            ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3])
        }
    }
}

/// **A few ulps of `x`, and the smallest normal** — how far an `f64` computed near `x` may round
/// from the value it stands for. A bound read in `f64` is widened by it to hold the truth:
/// [`cache_reach`]'s radius, and the self-touch sieve's rational circle boxes.
pub(crate) fn round_slack(x: f64) -> f64 {
    x.abs() * (4.0 * f64::EPSILON) + f64::MIN_POSITIVE
}

/// **[`cache_reach`] in rationals, its ends moved outward** — `None` where the reach is not finite
/// or an end's magnitude is past `Rat`.
///
/// An end is its exact value when `Rat` holds it. What `Rat` cannot hold of a finite `f64` is
/// either too large, which no outward move fixes, or below `2⁻⁷⁴` in magnitude — the exact value
/// fails only past a denominator of `2¹²⁶`, and a significand has at most 53 bits — so such an end
/// moves out to `±2⁻⁶⁴` or `0`, whichever is outward. The floor [`round_slack`] adds puts every
/// exactly-zero coordinate there.
pub(crate) fn rational_reach(
    cache: &nacre_topo::PointCache,
) -> Option<([nacre_exact::Rat; 3], [nacre_exact::Rat; 3])> {
    use nacre_exact::Rat;
    let (lo, hi) = cache_reach(cache);
    let tiny = Rat::new(1, 1i128 << 64)?;
    let zero = Rat::from_int(0);
    let outward = |x: f64, up: bool| -> Option<Rat> {
        if let Some(r) = Rat::try_from_f64(x) {
            return Some(r);
        }
        if !(x.is_finite() && x.abs() < 1.0) {
            return None;
        }
        Some(match (up, x > 0.0) {
            (true, true) => tiny,
            (false, false) => zero.checked_sub(tiny)?,
            _ => zero,
        })
    };
    let mut out = ([zero; 3], [zero; 3]);
    for k in 0..3 {
        out.0[k] = outward(lo[k], false)?;
        out.1[k] = outward(hi[k], true)?;
    }
    Some(out)
}

/// **One boundary piece of a face, in whichever exact spelling it has — the single reader.**
///
/// **One reader for both roads.** `face_clears_footprint` and the tangency road's
/// `line_runs_through_face` both read a face's corners here, where motion chains fold and `Pierce`
/// corners solve. A second reader that took rational vertices and *skipped* everything else
/// would read a face ringed by tangent corners from nothing and answer «does not reach the
/// line», which acquits.
///
/// **A whole disk is a piece too.** A loop that is a single circle carried by a cylinder
/// perpendicular to the face's plane *is* that disk — the face itself when it is the outer loop, a
/// bore through the face when it is a hole — and its reach along any direction is `centre ± ρ`.
/// A straight-edge demand is not needed: «all vertices on one side» would be satisfied by a single
/// point and would pass a disk that crosses the strip, but the piece answers for its own extent
/// instead of for one point of it. `lp` is the loop `he` belongs to, because «a single circle» is
/// a fact about that loop and not about the face.
pub(super) fn corner_of(
    model: &Model,
    face: &Face,
    lp: &nacre_topo::Loop,
    he: &nacre_topo::HalfEdge,
) -> Result<Corner, CornerFail> {
    match model.edge_curve(he.edge) {
        nacre_geom::Curve::Line(_) => {}
        // A single circular edge is the whole loop, so the loop is a disk — anything else
        // curved bulges past a hull this road cannot state.
        nacre_geom::Curve::Ellipse(_) if lp.half_edges.len() == 1 => {
            return disk_of(model, face, he);
        }
        // **An arc is a piece too**. A loop that mixes lines and arcs — a filleted
        // outline is the everyday one — is spellable: the same circle a whole loop would be, cut
        // to the extent its two ends name (a face the reader cannot spell refuses the boolean).
        // A failure here says its own cause, the way the whole disk's does.
        nacre_geom::Curve::Ellipse(_) => {
            return arc_of(model, face, he);
        }
    }
    start_corner(model, he)
}

/// **The point a half-edge starts at, as a piece** — whatever the edge's curve: the corner
/// [`corner_of`] reads for a straight edge, and the end an arc leaves from for a reader that walks
/// a loop's points as well as its pieces (`cyl_gate`'s line reader).
///
/// **One door, whatever stopped the exact spelling.** A corner the exact spellings cannot reach
/// ([`exact_start`]) is still a point whose realization is proven, so it is read on that proof
/// ([`Corner::Bounded`]); only a cache that proves nothing finite in `Rat` leaves the cause
/// standing. [`CornerFail::Inconsistent`] is not a missing spelling but a stored pierce name its
/// own meet does not offer — a defect — so it refuses as it is.
pub(super) fn start_corner(model: &Model, he: &nacre_topo::HalfEdge) -> Result<Corner, CornerFail> {
    let vh = he_start(model, *he);
    match exact_start(model, vh) {
        Err(cause) if cause != CornerFail::Inconsistent => {
            bounded_corner(model, vh, cause).ok_or(cause)
        }
        read => read,
    }
}

/// **A vertex's exact spelling** — rational coordinates carried to the world, or a pierce
/// corner's line and crossing — or why it has none.
fn exact_start(model: &Model, vh: Handle<Vertex>) -> Result<Corner, CornerFail> {
    match model.vertex_meet(vh) {
        // The point is stated in `frame`; the cylinder and the coefficients are world. A
        // chain that folds carries it out exactly — the same move the class descriptions and
        // the cylinder statements make, so this test sees a moved body the way every other
        // reader does. A `Wide` meet has no narrow vessel to move and declines, as does a chain
        // that does not fold: honest, never a comparison across two frames.
        Some((p, None)) => Ok(Corner::Rational(p)),
        Some((p, Some(leaf))) => {
            let q = p.narrow().ok_or(CornerFail::Width)?;
            let w = model.chain_point_rat(leaf, *q)?;
            Ok(Corner::Rational(nacre_exact::MeetPoint::Narrow(w)))
        }
        // **A corner a cylinder made has no rational meet — and does not need one.**
        // `vertex_meet` declines a `Pierce` because its coordinates are quadratic-irrational,
        // which is a statement about *rationals*, not about knowability: the point is exactly
        // `line.base() + s·line.dir()`, and the questions asked of it read that spelling
        // directly. Any other vertex has no exact spelling on this road — an `OnSeam` vertex
        // pins a curve, a three-plane one whose planes share no frame has no meet — and
        // [`start_corner`] reads it on its box.
        None => {
            let nacre_topo::Vertex::Pierce {
                planes,
                cylinder,
                rank,
            } = *model.vertex(vh)
            else {
                return Err(CornerFail::Shape);
            };
            pierce_corner(model, planes, cylinder, rank)
        }
    }
}

/// **A corner as the box its realization proves** ([`Corner::Bounded`]), with `cause` the reason no
/// exact spelling reached it; `None` where the cache proves no box `Rat` can state.
fn bounded_corner(model: &Model, vh: Handle<Vertex>, cause: CornerFail) -> Option<Corner> {
    use nacre_exact::Rat;
    let (lo, hi) = rational_reach(model.vertex_cache(vh))?;
    let half = Rat::new(1, 2)?;
    let mut centre = [Rat::from_int(0); 3];
    let mut rho = Rat::from_int(0);
    for k in 0..3 {
        centre[k] = lo[k].checked_add(hi[k])?.checked_mul(half)?;
        rho = rho.checked_add(hi[k].checked_sub(lo[k])?.checked_mul(half)?)?;
    }
    // A carrier's world name, for a plane; a cylinder carries no plane equation to be on.
    let planar = |h: Handle<Surface>| match model.surface(h) {
        nacre_topo::Surface::Plane { .. } => world_narrow(model, h).ok(),
        nacre_topo::Surface::Cylinder { .. } => None,
    };
    let carriers = match *model.vertex(vh) {
        Vertex::ThreePlane(t) => t.map(planar),
        Vertex::OnSeam(s) => [planar(s[0]), planar(s[1]), None],
        Vertex::Pierce { planes, .. } => [planar(planes[0]), planar(planes[1]), None],
        Vertex::CylinderPair { plane, .. } => [planar(plane), None, None],
    };
    Some(Corner::Bounded {
        centre,
        rho,
        carriers,
        cause,
    })
}

/// **One arc of a face's loop, as the round piece it is**.
///
/// Centre, radius and axis are [`disk_of`]'s — the same derivation a whole circle takes — and the
/// extent is the arc's two ends as **radial vectors**, which is the vocabulary [`RimArc`] and
/// [`arc_extent`] speak. The ends come from the **edge's stored order**, because that is what
/// `derive_edge_curve` orders counter-clockwise about the axis; the half-edge's traversal
/// direction says nothing about which arc this is, and two faces sharing the edge see the same one.
fn arc_of(model: &Model, face: &Face, he: &nacre_topo::HalfEdge) -> Result<Corner, CornerFail> {
    let Corner::Round {
        centre, rho2, axis, ..
    } = disk_of(model, face, he)?
    else {
        unreachable!("disk_of builds a round piece")
    };
    let [a, b] = model.edge(he.edge).vertices;
    let radial = |vh: Handle<Vertex>| -> Result<[nacre_exact::Rat; 3], CornerFail> {
        let p = vertex_point(model, vh)?;
        let mut v = p;
        for k in 0..3 {
            v[k] = v[k].checked_sub(centre[k]).ok_or(CornerFail::Width)?;
        }
        Ok(v)
    };
    Ok(Corner::Round {
        centre,
        rho2,
        axis,
        theta: Theta::Arc(RimArc {
            from: radial(a)?,
            to: radial(b)?,
        }),
    })
}

/// **A vertex's rational coordinates, when it has them** — the *narrow* reading, for a caller that
/// needs a **vector** from this point rather than a point to judge.
///
/// [`corner_of`] deliberately does not narrow: a `Wide` meet still judges exactly, and a `Pierce`
/// corner answers through its own line-and-crossing spelling. An arc's radial vector is arithmetic
/// on coordinates, so it needs them — a `Wide` meet declines for its width, and a pierce whose
/// crossing is irrational for its shape: no rational vector reaches it.
fn vertex_point(model: &Model, vh: Handle<Vertex>) -> Result<[nacre_exact::Rat; 3], CornerFail> {
    match model.vertex_meet(vh) {
        Some((p, None)) => p.narrow().copied().ok_or(CornerFail::Width),
        Some((p, Some(leaf))) => {
            let q = p.narrow().ok_or(CornerFail::Width)?;
            Ok(model.chain_point_rat(leaf, *q)?)
        }
        None => {
            let nacre_topo::Vertex::Pierce {
                planes,
                cylinder,
                rank,
            } = *model.vertex(vh)
            else {
                return Err(CornerFail::Shape);
            };
            let Corner::Pierce(line, s) = pierce_corner(model, planes, cylinder, rank)? else {
                unreachable!("pierce_corner builds a pierce piece")
            };
            let sr = s.as_rat().ok_or(CornerFail::Shape)?;
            line.point_at(sr).ok_or(CornerFail::Width)
        }
    }
}

/// **Where a round piece stands relative to the strip** — the extent form of
/// [`nacre_exact::cylinder_strip_side_margin`], read along the strip's own direction `e = n × m`.
///
/// The piece's reach along `e` is [`arc_extent`]'s, and because the arc lies in a plane whose
/// normal is its own axis, `e ⊥ axis` makes `e⊥ = e` — so a side that reaches its full radial peak
/// has margin exactly `ρ`, the number the door already takes, and a side that stops at an arc end
/// has margin `0` and a **rational offset** along `e`. An offset is exact as a *point*: shifting
/// the centre by `off/(e·e) · e` stays on the plane and moves `U` by `off`.
fn round_strip_side(
    piece: &Corner,
    coeffs: &[nacre_exact::Rat; 4],
    o: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    r2: &nacre_exact::BigRat,
) -> Option<nacre_exact::StripSide> {
    use nacre_exact::MeetPoint;
    let Corner::Round {
        centre,
        rho2,
        axis,
        theta,
    } = piece
    else {
        return None;
    };
    if *theta == Theta::Whole {
        // A whole circle reaches alike both ways: the symmetric door, complete answer and all.
        return Some(nacre_exact::cylinder_strip_side_margin(
            coeffs,
            &MeetPoint::Narrow(*centre),
            rho2,
            o,
            m,
            r2,
        ));
    }
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let e = nacre_exact::cross3_rat(&n, m)?;
    let (lo, hi) = arc_ends_along(centre, rho2, axis, theta, &e)?;
    Some(nacre_exact::cylinder_strip_side_extent(
        coeffs,
        &nacre_exact::StripReach {
            lo: (&MeetPoint::Narrow(lo.0), &lo.1),
            hi: Some((&MeetPoint::Narrow(hi.0), &hi.1)),
        },
        o,
        m,
        r2,
    ))
}

/// One end of a round piece's reach along a direction: a point on the piece's plane and the
/// **squared** margin the doors add to it — `(p, ρ²)`, the pair every scalar door here takes.
pub(crate) type ArcEnd = ([nacre_exact::Rat; 3], nacre_exact::BigRat);

/// **The two ends of a round piece's reach along `d`, each as a point and a margin** — the form
/// both scalar doors take.
///
/// [`arc_extent`] states the reach as offsets and squared radicals about the centre; the doors
/// want `(point, ρ²)` pairs. The bridge is exact and needs no new arithmetic:
/// * an offset becomes a **moved point** — `centre + off/(d·d)·d` stays on the piece's plane and
///   moves `d·p` by exactly `off`, whatever origin the door measures from;
/// * a side that reaches its full radial peak has margin `ρ²` (because the piece's plane has the
///   carrier's axis as its normal, `d⊥ = d` there and `rho2 = ρ²|d|²`), and a side an arc end
///   stopped has margin `0` with the offset carrying it.
///
/// [`Theta::Whole`] is the whole circle: both ends are the centre with margin `ρ`, which is what
/// the symmetric doors have always been handed.
pub(crate) fn arc_ends_along(
    centre: &[nacre_exact::Rat; 3],
    rho2: &nacre_exact::BigRat,
    axis: &[nacre_exact::Rat; 3],
    theta: &Theta,
    d: &[nacre_exact::Rat; 3],
) -> Option<(ArcEnd, ArcEnd)> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let dd = nacre_exact::dot3_rat(d, d)?;
    let (lo_off, rho2_lo, hi_off, rho2_hi) = arc_extent(theta, rho2, axis, d)?;
    let shifted = |off: Rat| -> Option<[Rat; 3]> {
        if off == zero {
            return Some(*centre);
        }
        let k = off.checked_mul(Rat::new(dd.denom(), dd.numer())?)?;
        let mut p = *centre;
        for i in 0..3 {
            p[i] = p[i].checked_add(k.checked_mul(d[i])?)?;
        }
        Some(p)
    };
    let margin = |reach2: Rat| {
        if reach2 == zero {
            nacre_exact::BigRat::zero()
        } else {
            rho2.clone()
        }
    };
    Some((
        (shifted(lo_off)?, margin(rho2_lo)),
        (shifted(hi_off)?, margin(rho2_hi)),
    ))
}

/// **The disk a single-circle face is** — centre and radius from the exact statements, never the
/// cache.
///
/// **The centre is `axis ∩ the face's *own* plane`, not the class's.** The callers check the
/// piece against the class's coefficients before measuring it, because a face can be merged into a
/// class by *rounded* coefficients its own points do not satisfy. Deriving the centre from the
/// class would make that check pass by construction and quietly retire it.
///
/// [`CornerFail::Shape`] is any shape this is not: a face that is not planar, an edge no cylinder
/// carries, or an axis not perpendicular to the plane — that last one traces an **ellipse**, and
/// this piece would be claiming to know a shape it does not. The plane's and the cylinder's world
/// statements, and the arithmetic of the centre, fail by their own causes.
fn disk_of(model: &Model, face: &Face, he: &nacre_topo::HalfEdge) -> Result<Corner, CornerFail> {
    // Asked of the **truth**, like the `world_plane_name` on the very next line: a face's
    // kind is a fact about what it *is*, and the cache is a rounded copy of that. Before, this
    // one function asked the cache what kind it was and then the truth what it said.
    if !matches!(
        model.surface(face.surface),
        nacre_topo::Surface::Plane { .. }
    ) {
        return Err(CornerFail::Shape);
    }
    let plane = world_narrow(model, face.surface)?;
    let n = [plane[0], plane[1], plane[2]];
    let cylinder = *model
        .edge(he.edge)
        .surfaces
        .iter()
        .find(|&&s| matches!(model.surface(s), nacre_topo::Surface::Cylinder { .. }))
        .ok_or(CornerFail::Shape)?;
    let def = world_cylinder(model, cylinder)?;
    let (o, m) = (def.origin(), def.dir());
    if !class_carries_circle(&n, &m) {
        return Err(CornerFail::Shape);
    }
    // The axis is along the normal, so `n·m ≠ 0` and the meet fails only on width.
    let centre = nacre_exact::axis_plane_meet(&plane, &o, &m).ok_or(CornerFail::Width)?;
    Ok(Corner::Round {
        centre,
        rho2: def.r2().clone(),
        axis: m,
        theta: Theta::Whole,
    })
}

/// **The name a face this road cannot spell refuses under** — the boundary decides it, the answer
/// is not touched.
///
/// Not every abstention is `CylinderGateUndecided`, whose sentence is "the gate could not decide
/// exactly". For a face whose ring runs along a cylinder that sentence is **false**: the gate can
/// decide, and what cannot read the shape is the road behind it — `combinatorics::loop_triples`
/// declines such a ring by name (`CurvedOperandBoundary`). This population is a boolean's result
/// fed back in as an operand, which is why a *gate* is the first thing to meet it.
///
/// **An edge's carriers are the adjacency answer** — [`nacre_topo::Edge::surfaces`] says so in its
/// own doc ("the two faces that use the edge"), so a cylinder among them is exactly "the face
/// across this edge is a lateral". That is the same fact `loop_triples` reads through classes,
/// spelled here from the model because this side has no class table for the operand.
///
/// **And it carries that road's exception with it — a *single-edge* loop is a full circle, which
/// the road does speak** (`LoopRing::Circle`, named by the cylinder's class). Measured: without the
/// `len() > 1`, the crosswise bore's cap face — one arc and one seam vertex — took
/// `CurvedOperandBoundary`, and for that face its sentence is simply false. The two conditions have
/// to agree, or this one is claiming a limit the other does not have.
///
/// This is for a shape the road cannot spell ([`CornerFail::Shape`] — an arc end at an irrational
/// root, a corner whose box does not stand or cannot answer). A statement that cannot be put in
/// the world in one frame, or one the world holds and `Rat` does not, has its own name whatever the
/// boundary looks like ([`CornerFail::stated_reason`]) — and since a **pierce corner is
/// readable**, calling those "curved" would put a false sentence on a true refusal.
pub(super) fn unreadable_reason(model: &Model, face: &Face) -> RejectReason {
    let curved = face.outer.half_edges.len() > 1
        && face.outer.half_edges.iter().any(|he| {
            model
                .edge(he.edge)
                .surfaces
                .iter()
                .any(|&s| matches!(model.surface(s), nacre_topo::Surface::Cylinder { .. }))
        });
    if curved {
        RejectReason::CurvedOperandBoundary
    } else {
        RejectReason::CylinderGateUndecided
    }
}

/// Whether one planar face misses the **rectangle** this cylinder occupies in the face's plane.
///
/// **One question, two separating axes.** A plane parallel to the axis meets the solid cylinder
/// in a rectangle: the strip across ([`nacre_exact::cylinder_strip_side`]) and a lateral face's
/// axis-parameter span along ([`nacre_exact::point_axis_side`]). A rectangle is the intersection
/// of those two bands, so clearing *either* axis clears it — these are not two rules to be
/// weighed, they are the two axes of one. With several lateral faces there are several rectangles
/// (the strip is shared, the spans are not), and the face must miss them all:
///
/// ```text
/// clear  ⟺  clear across the strip  ∨  clear along **every** span
/// ```
///
/// **"Every" over an empty list is a pass, and that one is not on offer.** With no usable
/// span the axis along is skipped outright rather than answered vacuously. Not a hypothetical:
/// measured with the emptiness guard removed, a wall that genuinely crosses the bore and one
/// tangent to it were both waved straight through.
///
/// **The span reading is an open interval** — a face resting exactly on a cap plane is clear,
/// because the theorem being fed speaks of the *open* slab. What such a face touches is the rim's
/// own plane, and whether its edge crosses the rim circle there is a question the arrangement asks
/// where the circles and segments are, by name.
///
/// **The two axes are only *complete* for a face that is an axis-aligned rectangle** — then both
/// rectangles' edge normals coincide and there are no other separating axes to try. An extruded
/// wall is exactly that shape (its edges run along the axis or across it), chamfers included: what
/// a chamfer tilts is the plane, not the edges within it. A face a previous boolean took a bite
/// out of is not, and neither is a slanted or L-shaped one; those are refused as "not shown to
/// clear", which is true. Completing the test means adding the face's own edge normals as further
/// axes — the same test with more axes, not a different machine — and waits for a shape that
/// needs it. The cylinder **pair**'s version of this sentence is closed
/// (`cyl_gate::separating_dirs`); this arm, the plane's, still waits.
///
/// **Only the outer loop is walked**, and that is sound: a face is contained in the convex hull of
/// its outer loop's **pieces** (a vertex is the piece with no width), a half-space is convex, and
/// inner loops only *remove* material. So a hole is never read, and a hole whose pieces the reader
/// cannot spell cannot refuse the face.
///
/// **A straight-edge demand would be soundness, and it is paid for rather than kept.** This
/// plane may be some *other* cylinder's cap plane (one whose axis is perpendicular to it), and
/// such a cap face's outer loop is a single circle with one seam vertex — "all vertices on one
/// side" would be satisfied by a single **point** and would pass a disk that crosses the strip.
/// A barrier that stops it also stops every disk that genuinely clears, and folded to
/// «did not clear» writes tangency rows for contacts that are not there. So the piece answers for
/// its own **extent** ([`Corner::Round`]): the hull statement below holds because a face is
/// contained in the hull of its boundary, and each boundary piece reports how far it reaches.
/// Anything else curved is still refused — an arc that is only part of a loop bulges past a hull
/// this road cannot state.
#[allow(clippy::too_many_arguments)]
pub(super) fn face_clears_footprint(
    model: &Model,
    face: &Face,
    coeffs: &[nacre_exact::Rat; 4],
    o: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    r2: &nacre_exact::BigRat,
    spans: &[[nacre_exact::Rat; 2]],
) -> Result<bool, BoolError> {
    use nacre_exact::{Orient, StripSide};
    // A piece this road cannot read refuses under its cause's name
    // ([`CornerFail::reason_on`]); a class whose coefficients miss its own face is the gate's own
    // (`misfit` below).
    let misfit = || reject(RejectReason::CylinderGateUndecided);
    let mut side: Option<StripSide> = None;
    let mut across = true;
    // Per span, whether every vertex so far has stayed at or below its start, and at or above its
    // end. Either one surviving the walk clears that rectangle along the axis.
    let mut along: Vec<(bool, bool)> = vec![(true, true); spans.len()];
    let mut vertices = 0usize;
    for he in &face.outer.half_edges {
        // **The piece's description is chosen once, in one reader** — the three questions
        // below then ask it the same things whichever spelling it wears, and the straddle road
        // reads it the very same way.
        let corner = corner_of(model, face, &face.outer, he)
            .map_err(|e| reject(e.reason_on(model, face)))?;
        // **The class's coefficients must actually describe *this* face's plane.** Classes merge
        // on three exact witnesses, one of which compares *rounded* coefficients — so a face can
        // sit in a class whose exact name its own vertices do not satisfy (the two-descriptions
        // hazard). The strip decomposition takes the plane's
        // distance from the axis as the point's, so judging a point against a plane it is not on
        // would answer about geometry that is not there. Checked, not assumed: a `debug_assert`
        // would say nothing in the build that ships. It is also `cylinder_strip_side`'s own
        // precondition, so it comes first.
        if !corner.on_plane(coeffs) {
            return Err(misfit());
        }
        vertices += 1;
        // The axis across the strip.
        //
        // **The two clear sides are spelled out rather than caught.** A catch-all here would
        // take any *future* answer for "clear, on some side" — and the answer this test is about
        // to grow is the opposite one (a piece with width can straddle the strip by itself), which
        // a catch-all would record as a side and let the face pass. Written this way the compiler
        // asks the question at every new variant instead.
        match corner.strip_side(coeffs, o, m, r2) {
            // Reaching the strip and spanning it are different facts (the tangency road needs them
            // apart), but for *clearing a rectangle* they are the same one: not clear — and so is
            // a side the piece could not prove.
            None | Some(StripSide::Inside | StripSide::Crosses) => across = false,
            Some(s @ (StripSide::Plus | StripSide::Minus)) => match side {
                None => side = Some(s),
                Some(prev) if prev != s => across = false, // the face straddles the strip
                Some(_) => {}
            },
        }
        // The axis along it — one rectangle per lateral face, and the span is open at both ends.
        for (i, span) in spans.iter().enumerate() {
            if corner.reaches(o, m, span[0], Orient::Positive) {
                along[i].0 = false;
            }
            if corner.reaches(o, m, span[1], Orient::Negative) {
                along[i].1 = false;
            }
        }
    }
    // An empty outer loop names no half-space and no interval, so it proves nothing — and
    // "every vertex of none stayed below" would otherwise be a vacuous pass for any wall.
    if vertices == 0 {
        return Ok(false);
    }
    if across && side.is_some() {
        return Ok(true);
    }
    // `spans` empty means the axis is unusable, not that every rectangle was missed.
    Ok(!spans.is_empty() && along.iter().all(|(below, above)| *below || *above))
}

#[cfg(test)]
#[path = "../tests/corners.rs"]
mod tests;
