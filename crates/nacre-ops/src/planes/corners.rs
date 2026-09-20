use super::*;
/// Whether one planar face misses the **rectangle** this cylinder occupies in the face's plane.
///
/// ★★ **One question, two separating axes.** A plane parallel to the axis meets the solid cylinder
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
/// ★★★ **"Every" over an empty list is a pass, and that one is not on offer.** With no usable
/// span the axis along is skipped outright rather than answered vacuously. Not a hypothetical:
/// measured with the emptiness guard removed, a wall that genuinely crosses the bore and one
/// tangent to it were both waved straight through.
///
/// ★ **The span reading is an open interval** — a face resting exactly on a cap plane is clear,
/// because the theorem being fed speaks of the *open* slab. What such a face touches is the rim's
/// own plane, and whether its edge crosses the rim circle there is a question the arrangement asks
/// where the circles and segments are, by name.
///
/// ★ **The two axes are only *complete* for a face that is an axis-aligned rectangle** — then both
/// rectangles' edge normals coincide and there are no other separating axes to try. An extruded
/// wall is exactly that shape (its edges run along the axis or across it), chamfers included: what
/// a chamfer tilts is the plane, not the edges within it. A face a previous boolean took a bite
/// out of is not, and neither is a slanted or L-shaped one; those are refused as "not shown to
/// clear", which is true. Completing the test means adding the face's own edge normals as further
/// axes — the same test with more axes, not a different machine — and waits for a shape that
/// needs it. ★ The cylinder **pair**'s version of this sentence is closed
/// ([`separating_dirs`]); this arm, the plane's, still waits.
///
/// **Only the outer loop is walked**, and that is sound: a face is contained in the convex hull of
/// its outer loop's **pieces** (a vertex is the piece with no width), a half-space is convex, and
/// inner loops only *remove* material. It
/// is also what keeps the test reachable — a wall drilled by a crosswise bore carries that bore's
/// rim as an inner loop, whose vertices have no rational meet at all.
///
/// ★★ **A straight-edge demand would be soundness, and it is paid for rather than kept.** This
/// plane may be some *other* cylinder's cap plane (one whose axis is perpendicular to it), and
/// such a cap face's outer loop is a single circle with one seam vertex — "all vertices on one
/// side" would be satisfied by a single **point** and would pass a disk that crosses the strip.
/// A barrier that stops it also stops every disk that genuinely clears, and folded to
/// «did not clear» writes tangency rows for contacts that are not there. So the piece answers for
/// its own **extent** ([`Corner::Disk`]): the hull statement below holds because a face is
/// contained in the hull of its boundary, and each boundary piece reports how far it reaches.
/// Anything else curved is still refused — an arc that is only part of a loop bulges past a hull
/// this road cannot state.
/// **A face corner, in whichever exact spelling it has** — the footprint test's currency.
///
/// Three plane-carried corners have rational coordinates; a corner a cylinder made does not, and
/// is `line.base() + s·line.dir()` with `s` quadratic-irrational. The two spellings answer the
/// same three questions, so they are asked through one type rather than branched at each call.
pub(super) enum Corner {
    Rational(nacre_exact::MeetPoint),
    Pierce(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal),
    /// **A piece of a circle in the face's own plane** — the whole disk a single-circle loop *is*
    /// (`arc: None`), or one **arc** of a loop that mixes lines and arcs (`arc: Some`).
    ///
    /// ★ The name «corner» is historical: what this type carries is a **boundary piece's reach**,
    /// and a vertex is the piece with no width. Written this way the loops below do not learn a
    /// new shape; they ask the same three questions and one of the answers now has extent.
    ///
    /// ★★ `axis` is the **carrier cylinder's** axis, and `arc`'s `from → to` is counter-clockwise
    /// about it ([`lateral_theta_extent`] states that convention). The face's own normal may run
    /// the other way; using it would name the **complementary** arc, which is not a bound on this
    /// piece but a different set.
    Round {
        centre: [nacre_exact::Rat; 3],
        /// The carrier's squared radius.
        rho2: nacre_exact::BigRat,
        axis: [nacre_exact::Rat; 3],
        arc: Option<RimArc>,
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
        }
    }

    pub(super) fn strip_side(
        &self,
        coeffs: &[nacre_exact::Rat; 4],
        o: &[nacre_exact::Rat; 3],
        m: &[nacre_exact::Rat; 3],
        r2: &nacre_exact::BigRat,
    ) -> nacre_exact::StripSide {
        match self {
            Self::Rational(p) => nacre_exact::cylinder_strip_side(coeffs, p, o, m, r2),
            Self::Pierce(line, s) => {
                nacre_exact::cylinder_strip_side_branch(coeffs, line, s, o, m, r2)
            }
            // ★ The strip runs along `e = n × m`, so **that** is the direction the piece's extent
            // is asked for. A whole circle reaches `±ρ|e|` alike and takes the symmetric door,
            // which is the complete answer it always had; an arc's two ends differ and take the
            // general one. `None` from the extent is `Inside` — reached the strip, spanned
            // nothing — which is the safe reading for both consumers.
            Self::Round { .. } => {
                round_strip_side(self, coeffs, o, m, r2).unwrap_or(nacre_exact::StripSide::Inside)
            }
        }
    }

    /// **Does this piece reach strictly past the plane at axis parameter `t`, to the `want` side?**
    ///
    /// ★★★★★ **A point answers this with one sign; a piece with width cannot.** The caller asks
    /// twice — «did anything go above the span's start» and «did anything go below its end» — and
    /// for a vertex those are the two readings of a single `Orient`. A disk straddling a station
    /// is above it *and* below it, which no single `Orient` can say: answering `Zero` would make
    /// **both** questions false and let a face that covers the station be read as clearing it.
    /// So the question is named instead of inferred.
    ///
    /// ★ **The span is read open** (`face_clears_footprint`'s own note: a corner sitting exactly
    /// on a cap plane must not be dropped), which is why the disk's comparison is spelled here
    /// rather than borrowed from [`reach_clears`] — that one is closed against closed, and folding
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
            // The comparison quantity is `q = (p − o)·m − t(m·m)`; over the piece it sweeps
            // `q_c + [lo_off − √rho2_lo, hi_off + √rho2_hi]` ([`arc_extent`]), so «reaches past»
            // is that end on the wanted side, or its radical covering the gap —
            // `rho2 > end²`, squared once and rational throughout. A whole circle is the
            // symmetric case (`off = 0`, `rho2 = ρ²|m⊥|²`), which is what this said before an
            // arc could be a piece.
            // ★ Overflow answers **`true`**: "not shown to clear" is this test's safe direction.
            Self::Round {
                centre,
                rho2,
                axis,
                arc,
            } => {
                let zero = nacre_exact::Rat::from_int(0);
                let end = (|| {
                    let mut rel = [zero; 3];
                    for k in 0..3 {
                        rel[k] = centre[k].checked_sub(o[k])?;
                    }
                    let mm = dot3(m, m)?;
                    let q = dot3(&rel, m)?.checked_sub(t.checked_mul(mm)?)?;
                    let (lo_off, rho2_lo, hi_off, rho2_hi) =
                        arc_extent(arc.as_ref(), rho2, axis, m)?;
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
/// ★★★ **No restatement is owed here.** `QuadRoot` is written about the canonical direction of the
/// meet line of *the two planes the definition names, in the order it names them* — and this reads
/// exactly those, in that order, so the root applies directly. (The `ℓ` correction
/// `combinatorics::pierce_name_from_def` carries is the price of crossing from handle space into a
/// class table's order; this side has no class table for the operand at all.)
///
/// ★★ **`None` is only ever a missing *description***, never a shape this road cannot spell — the
/// caller has already established that the vertex is a `Pierce`, so the two refusals stay apart:
/// a plane whose world name is not narrow, a cylinder with no world statement, or a root the meet
/// does not offer are the gate's own arithmetic running out (`CylinderGateUndecided`), while a
/// seam vertex never reaches here at all.
fn pierce_corner(
    model: &Model,
    planes: [Handle<Surface>; 2],
    cylinder: Handle<Surface>,
    root: nacre_topo::QuadRoot,
) -> Option<Corner> {
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let p1 = *model.world_plane_name(planes[0])?.narrow()?;
    let p2 = *model.world_plane_name(planes[1])?.narrow()?;
    let def = world_cylinder_def(model, cylinder)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (line, s) = match (
        nacre_exact::quad::plane_plane_cylinder(&p1, &p2, &o, &m, r2)?,
        root,
    ) {
        (CylinderMeet::Pair { line, s }, QuadRoot::Lo) => (line, s[0]),
        (CylinderMeet::Pair { line, s }, QuadRoot::Hi) => (line, s[1]),
        (CylinderMeet::Tangent { line, s }, QuadRoot::Double) => (line, QuadVal::from_rat(s)),
        _ => return None,
    };
    Some(Corner::Pierce(line, s))
}

/// Why a boundary piece could not be read — the two causes the footprint road keeps apart.
pub(super) enum CornerFail {
    /// A shape this road cannot spell at all — a seam vertex. ★ An arc that
    /// is not a whole disk is a piece, so what is
    /// left here is a vertex with no exact name of any kind.
    Shape,
    /// The description ran out: a chain that will not fold, a name that is not narrow.
    Arithmetic,
}

/// **One boundary piece of a face, in whichever exact spelling it has — the single reader.**
///
/// ★★★★★ **It was written twice, and the second copy was an abbreviation.**
/// `face_clears_footprint` folded motion chains and solved `Pierce` corners; `face_straddles_line`
/// took rational vertices and *skipped* everything else — so a face ringed by tangent corners had
/// its straddle read from nothing and answered «does not straddle», which acquits. One reader, and
/// the two roads can no longer disagree about what a face says.
///
/// ★★ **A whole disk is a piece too.** A face whose outer loop is a single circle carried by a
/// cylinder perpendicular to its plane *is* that disk, and its reach along any direction is
/// `centre ± ρ`. The doc of [`face_clears_footprint`] named this shape as the reason its
/// straight-edge demand was soundness rather than convenience — «all vertices on one side would be
/// satisfied by a single point and would pass a disk that crosses the strip». The demand can go
/// now because the piece answers for its own extent instead of for one point of it.
pub(super) fn corner_of(
    model: &Model,
    face: &Face,
    he: &nacre_topo::HalfEdge,
) -> Result<Corner, CornerFail> {
    match model.edge_curve(he.edge) {
        nacre_geom::Curve::Line(_) => {}
        // A single circular edge is the whole boundary, so the face is a disk — anything else
        // curved bulges past a hull this road cannot state.
        nacre_geom::Curve::Circle(_) if face.outer.half_edges.len() == 1 => {
            return disk_of(model, face, he).ok_or(CornerFail::Shape);
        }
        // ★★★★★ **An arc is a piece too**. A loop that mixes lines and arcs — a filleted
        // outline is the everyday one — used to be unreadable here whatever the arc was or where
        // it sat, and a face the reader cannot spell refuses the boolean. The shape is spellable:
        // the same circle a whole loop would be, cut to the extent its two ends name. So a failure
        // now is `Arithmetic` — a value that could not be stated — and never `Shape`.
        nacre_geom::Curve::Circle(_) => {
            return arc_of(model, face, he).ok_or(CornerFail::Arithmetic);
        }
    }
    let vh = he_start(model, *he);
    match model.vertex_meet(vh) {
        // ★ The point is stated in `frame`; the cylinder and the coefficients are world. A
        // chain that folds to a rational translation carries it out exactly — the same move
        // the class descriptions and the cylinder statements make, so this test sees a
        // translated body the way every other reader does. A `Wide` meet has no narrow vessel
        // to shift and declines, as does any other chain: honest, never a comparison across
        // two frames.
        Some((p, None)) => Ok(Corner::Rational(p)),
        Some((p, Some(leaf))) => {
            let (Some(t), Some(q)) = (model.chain_translation(leaf), p.narrow()) else {
                return Err(CornerFail::Arithmetic);
            };
            let mut w = *q;
            for (c, d) in w.iter_mut().zip(t) {
                *c = c.checked_add(d).ok_or(CornerFail::Arithmetic)?;
            }
            Ok(Corner::Rational(nacre_exact::MeetPoint::Narrow(w)))
        }
        // ★★★★★ **A corner a cylinder made has no rational meet — and does not need one.**
        // `vertex_meet` declines a `Pierce` because its coordinates are quadratic-irrational,
        // which is a statement about *rationals*, not about knowability: the point is exactly
        // `line.base() + s·line.dir()`, and the questions asked of it read that spelling
        // directly. An `OnSeam` vertex is a different matter — it pins a curve, not a point —
        // so it stays unreadable.
        None => {
            let nacre_topo::Vertex::Pierce {
                planes,
                cylinder,
                root,
            } = *model.vertex(vh)
            else {
                return Err(CornerFail::Shape);
            };
            pierce_corner(model, planes, cylinder, root).ok_or(CornerFail::Arithmetic)
        }
    }
}

/// **One arc of a face's outer loop, as the round piece it is**.
///
/// Centre, radius and axis are [`disk_of`]'s — the same derivation a whole circle takes — and the
/// extent is the arc's two ends as **radial vectors**, which is the vocabulary [`RimArc`] and
/// [`arc_extent`] speak. The ends come from the **edge's stored order**, because that is what
/// `derive_edge_curve` orders counter-clockwise about the axis; the half-edge's traversal
/// direction says nothing about which arc this is, and two faces sharing the edge see the same one.
fn arc_of(model: &Model, face: &Face, he: &nacre_topo::HalfEdge) -> Option<Corner> {
    let Corner::Round {
        centre, rho2, axis, ..
    } = disk_of(model, face, he)?
    else {
        return None;
    };
    let [a, b] = model.edge(he.edge).vertices;
    let radial = |vh: Handle<Vertex>| -> Option<[nacre_exact::Rat; 3]> {
        let p = vertex_point(model, vh)?;
        let mut v = p;
        for k in 0..3 {
            v[k] = v[k].checked_sub(centre[k])?;
        }
        Some(v)
    };
    Some(Corner::Round {
        centre,
        rho2,
        axis,
        arc: Some(RimArc {
            from: radial(a)?,
            to: radial(b)?,
        }),
    })
}

/// **A vertex's rational coordinates, when it has them** — the *narrow* reading, for a caller that
/// needs a **vector** from this point rather than a point to judge.
///
/// [`corner_of`] deliberately does not narrow: a `Wide` meet still judges exactly, and a `Pierce`
/// corner answers through its own line-and-root spelling. An arc's radial vector is arithmetic on
/// coordinates, so it needs them — a `Wide` meet or a pierce whose root is irrational declines,
/// and the caller says so by name.
fn vertex_point(model: &Model, vh: Handle<Vertex>) -> Option<[nacre_exact::Rat; 3]> {
    match model.vertex_meet(vh) {
        Some((p, None)) => p.narrow().copied(),
        Some((p, Some(leaf))) => {
            let (t, q) = (model.chain_translation(leaf)?, p.narrow()?);
            let mut w = *q;
            for (c, d) in w.iter_mut().zip(t) {
                *c = c.checked_add(d)?;
            }
            Some(w)
        }
        None => {
            let nacre_topo::Vertex::Pierce {
                planes,
                cylinder,
                root,
            } = *model.vertex(vh)
            else {
                return None;
            };
            let Corner::Pierce(line, s) = pierce_corner(model, planes, cylinder, root)? else {
                return None;
            };
            let sr = s.as_rat()?;
            let (b, d) = (line.base(), line.dir());
            let mut out = b;
            for k in 0..3 {
                out[k] = out[k].checked_add(sr.checked_mul(d[k])?)?;
            }
            Some(out)
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
        arc,
    } = piece
    else {
        return None;
    };
    let arc = arc.as_ref();
    if arc.is_none() {
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
    let e = combinatorics::cross3_rat(&n, m)?;
    let (lo, hi) = arc_ends_along(centre, rho2, axis, arc, &e)?;
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
/// `arc = None` is the whole circle: both ends are the centre with margin `ρ`, which is what the
/// symmetric doors have always been handed.
pub(crate) fn arc_ends_along(
    centre: &[nacre_exact::Rat; 3],
    rho2: &nacre_exact::BigRat,
    axis: &[nacre_exact::Rat; 3],
    arc: Option<&RimArc>,
    d: &[nacre_exact::Rat; 3],
) -> Option<(ArcEnd, ArcEnd)> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let dd = dot3(d, d)?;
    let (lo_off, rho2_lo, hi_off, rho2_hi) = arc_extent(arc, rho2, axis, d)?;
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
/// ★★ **The centre is `axis ∩ the face's *own* plane`, not the class's.** The callers check the
/// piece against the class's coefficients before measuring it, because a face can be merged into a
/// class by *rounded* coefficients its own points do not satisfy. Deriving the centre from the
/// class would make that check pass by construction and quietly retire it.
///
/// `None` is any shape this is not: a face that is not planar, an edge no cylinder carries, or an
/// axis not perpendicular to the plane — that last one traces an **ellipse**, and this piece would
/// be claiming to know a shape it does not.
fn disk_of(model: &Model, face: &Face, he: &nacre_topo::HalfEdge) -> Option<Corner> {
    // ★ Asked of the **truth**, like the `world_plane_name` on the very next line: a face's
    // kind is a fact about what it *is*, and the cache is a rounded copy of that. Before, this
    // one function asked the cache what kind it was and then the truth what it said.
    if !matches!(
        model.surface(face.surface),
        nacre_topo::Surface::Plane { .. }
    ) {
        return None;
    }
    let plane = *model.world_plane_name(face.surface)?.narrow()?;
    let n = [plane[0], plane[1], plane[2]];
    let def = model
        .edge(he.edge)
        .surfaces
        .iter()
        .find(|&&s| matches!(model.surface(s), nacre_topo::Surface::Cylinder { .. }))
        .and_then(|&s| world_cylinder_def(model, s))?;
    let (o, m) = (def.origin(), def.dir());
    if !nacre_exact::parallel_rat(&n, &m) {
        return None;
    }
    // `n·(o + t·m) + d = 0`, and `n·m ≠ 0` because the axis is along the normal.
    let nm = dot3(&n, &m)?;
    let no_d = dot3(&n, &o)?.checked_add(plane[3])?;
    let t = nacre_exact::Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(nacre_exact::Rat::new(nm.denom(), nm.numer())?)?;
    let mut centre = o;
    for k in 0..3 {
        centre[k] = centre[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(Corner::Round {
        centre,
        rho2: def.r2().clone(),
        axis: m,
        arc: None,
    })
}

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
    // ★★★★★ **Which refusal this is, decided once — the answer is not touched.**
    //
    // Every abstention below turns into a refusal at the caller, and until now they all wore
    // `CylinderGateUndecided`, whose sentence is "the gate could not decide exactly". For a face
    // whose ring runs along a cylinder that sentence is **false**: the gate can decide, and what
    // cannot read the shape is the road behind it — `combinatorics::loop_triples` declines such a
    // ring by name (`CurvedOperandBoundary`). This population is a boolean's result fed back in as
    // an operand, which is why a *gate* was the first thing to meet it.
    //
    // ★★ **An edge's carriers are the adjacency answer** — [`nacre_topo::Edge::surfaces`] says so
    // in its own doc ("the two faces that use the edge"), so a cylinder among them is exactly "the
    // face across this edge is a lateral". That is the same fact `loop_triples` reads through
    // classes, spelled here from the model because this side has no class table for the operand.
    //
    // ★★★ **And it carries that road's exception with it — a *single-edge* loop is a full circle,
    // which the road does speak** (`LoopRing::Circle`, named by the cylinder's class). Measured:
    // without the `len() > 1`, the crosswise bore's cap face — one arc and one seam vertex — took
    // the new name, and for that face the new name's sentence is simply false. The two conditions
    // have to agree, or this one is claiming a limit the other does not have.
    let curved = face.outer.half_edges.len() > 1
        && face.outer.half_edges.iter().any(|he| {
            model
                .edge(he.edge)
                .surfaces
                .iter()
                .any(|&s| matches!(model.surface(s), nacre_topo::Surface::Cylinder { .. }))
        });
    // ★★★★★ **Two refusals, split by cause rather than by the flag.** `unreadable` is for a shape
    // this road cannot spell at all — an arc edge, a seam vertex — and there
    // `CurvedOperandBoundary`'s sentence ("the road behind cannot read this ring") is true.
    // `arithmetic` is for the gate's own description running out: a chain that will not fold, a
    // name that is not narrow, a class whose coefficients miss its own face. Those are
    // [`RejectReason::CylinderGateUndecided`] whatever the boundary looks like, because the road
    // behind has nothing to do with them — and since a **pierce corner is now readable**, calling
    // them "curved" would put a false sentence on a true refusal.
    let unreadable = || {
        reject(if curved {
            RejectReason::CurvedOperandBoundary
        } else {
            RejectReason::CylinderGateUndecided
        })
    };
    let arithmetic = || reject(RejectReason::CylinderGateUndecided);
    let mut side: Option<StripSide> = None;
    let mut across = true;
    // Per span, whether every vertex so far has stayed at or below its start, and at or above its
    // end. Either one surviving the walk clears that rectangle along the axis.
    let mut along: Vec<(bool, bool)> = vec![(true, true); spans.len()];
    let mut vertices = 0usize;
    for he in &face.outer.half_edges {
        // ★★★★ **The piece's description is chosen once, in one reader** — the three questions
        // below then ask it the same things whichever spelling it wears, and the straddle road
        // reads it the very same way.
        let corner = corner_of(model, face, he).map_err(|e| match e {
            CornerFail::Shape => unreadable(),
            CornerFail::Arithmetic => arithmetic(),
        })?; // ★ **The class's coefficients must actually describe *this* face's plane.** Classes merge
        // on three exact witnesses, one of which compares *rounded* coefficients — so a face can
        // sit in a class whose exact name its own vertices do not satisfy (the two-descriptions
        // hazard `FaceInfo::exact_coeffs` documents). The strip decomposition takes the plane's
        // distance from the axis as the point's, so judging a point against a plane it is not on
        // would answer about geometry that is not there. Checked, not assumed: a `debug_assert`
        // would say nothing in the build that ships. ★ It is also `cylinder_strip_side`'s own
        // precondition, so it comes first.
        if !corner.on_plane(coeffs) {
            return Err(arithmetic());
        }
        vertices += 1;
        // The axis across the strip.
        //
        // ★ **The two clear sides are spelled out rather than caught.** A catch-all here would
        // take any *future* answer for "clear, on some side" — and the answer this test is about
        // to grow is the opposite one (a piece with width can straddle the strip by itself), which
        // a catch-all would record as a side and let the face pass. Written this way the compiler
        // asks the question at every new variant instead.
        match corner.strip_side(coeffs, o, m, r2) {
            // Reaching the strip and spanning it are different facts (the tangency road needs them
            // apart), but for *clearing a rectangle* they are the same one: not clear.
            StripSide::Inside | StripSide::Crosses => across = false,
            s @ (StripSide::Plus | StripSide::Minus) => match side {
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
    // ★ An empty outer loop names no half-space and no interval, so it proves nothing — and
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
