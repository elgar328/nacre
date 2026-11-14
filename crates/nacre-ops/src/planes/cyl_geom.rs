use super::*;

/// **Does an increasing axis parameter move toward this class's "above"?** — where "above" is the
/// side the class's plane faces, which is the frame every cell label is written in.
///
/// Exact: the world name against the axis (integers, no overflow), turned by the name's sense.
/// It takes the class's [`WorldName`], so a class without one cannot ask — every reader already
/// stands on a class whose world coefficients it read.
///
/// **Ask this, do not re-derive it from the class's rational name alone.** The name is
/// *canonical* (first nonzero component positive), which points the other way from the plane on
/// half the classes; a rule spelled against it reads "above" backwards exactly there. Spelled that
/// way, it turned 36 tests red at once (measured).
///
/// One sign for the whole section, also oblique: along every generator the height grows toward
/// the side `n·m` points to, so the sign is the normal's against the axis. This is only asked of
/// a class not along the axis, so the product never vanishes.
pub(crate) fn plus_t_is_above(world: &WorldName, def: &nacre_topo::CylinderDef) -> bool {
    let [a, b, c, _] = world.name.coeff_ints();
    let along = nacre_exact::normal_sense(&[a, b, c], def.dir());
    debug_assert_ne!(
        along,
        nacre_exact::Orient::Zero,
        "asked of a class along the axis"
    );
    (along == nacre_exact::Orient::Positive) == (world.sense == nacre_topo::Orientation::Forward)
}

/// **Where a plane meets a cylinder's axis, as the axis parameter `t`** —
/// [`nacre_exact::axis_param_of_plane`] read on the statement's own origin and direction, so the
/// formula has one spelling for every reader; the readers that want the point itself (a rim's
/// centre) ask [`nacre_exact::axis_plane_meet`].
///
/// `None` when the plane is parallel to the axis (`n·m = 0`, no meeting point) or when the
/// checked `Rat` arithmetic overflows. A rule that lives inlined in one place while a second
/// site spells a reduced version of it is this repo's dominant defect shape, so every reader asks
/// this.
pub(crate) fn axis_param_of_plane(
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_exact::Rat> {
    nacre_exact::axis_param_of_plane(coeffs, &def.origin(), &def.dir())
}

/// **The band of its own axis a cylinder's closed trace on the plane `coeffs` sweeps** —
/// `[t − Δ, t + Δ]` in the axis parameter of the raw `def.dir()`, about the plane's axis meet `t`.
/// A point of the cylinder is `o + s·m + r·û` (`û` a unit radial vector), and on the plane
/// `s = t − r(n·û)/(n·m)`; `n·û` peaks at `|n⊥|`, the normal off the axis, so `Δ = r|n⊥|/|n·m|`,
/// formed as `Δ²` and rounded outward ([`sqrt_at_least`]) — never narrower than the
/// ellipse. Across the axis `Δ = 0`: the circle's own `[t, t]`, with no radius asked. `None` is a
/// plane along the axis, `Rat` overflow, or a radius wider than `Rat` on an oblique plane.
pub(crate) fn trace_band(
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 2]> {
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let m = def.dir();
    let t = axis_param_of_plane(coeffs, def)?;
    if class_carries_circle(&n, &m) {
        return Some([t, t]);
    }
    let nperp = nacre_exact::off_axis_rat(&n, &m)?;
    let nm = nacre_exact::dot3_rat(&n, &m)?;
    let nm2 = nm.checked_mul(nm)?;
    let delta2 = def
        .r2()
        .narrow()?
        .checked_mul(nacre_exact::dot3_rat(&nperp, &nperp)?)?
        .checked_mul(nacre_exact::Rat::new(nm2.denom(), nm2.numer())?)?;
    let delta = sqrt_at_least(delta2)?;
    Some([t.checked_sub(delta)?, t.checked_add(delta)?])
}

/// **Does this class carry that cylinder's *circle*?** — its normal is parallel to the axis, so
/// the section is a circle and not an ellipse. Exact and **total** ([`nacre_exact::axis_relation`]
/// clears denominators into `BigInt`), so a caller's `false` means the geometry, never the width.
/// It is the `Across` arm of the relation the gate reads three ways ([`cylinder_gate`]'s
/// `class_meet`), so a reader asking whether a class carries a circle asks this rather than
/// respelling the test.
///
/// **The one place this rule is named.** An arc of a cylinder on a class is the cylinder's section
/// by the class's plane — a circle when this answers `true`, an ellipse otherwise — and the readers
/// that only need the section's shape read it whatever the plane: an arc's direction, side and
/// extreme off the axis ([`nacre_exact::quad::axial_side`],
/// [`nacre_exact::quad::ccw_tangent_side`]), its box over the band the section sweeps
/// ([`trace_band`]), its order about the axis. The sites that read a **circle** ask this:
/// `circle_crosses_ruling`'s strip (an ellipse it does not measure), the circle crossings'
/// `segment_meets_cylinder` filter (its precondition `(p1 − p0)·m = 0` holds only across the
/// axis), [`crate::nesting`]'s rim witnesses and `disk_in_disk`, a planar face's own circle for
/// interior witnesses (`face_circle`), and the alias folds of a line along the axis, which a class
/// across it meets at its station.
///
/// **What holds the oblique arcs to what those roads can draw** ([`cylinder_gate`]): the gate runs
/// whenever the input has a cylinder at all and sweeps **every (plane class × cylinder)** pair; an
/// oblique pair whose lateral faces cannot be *proved* to miss the plane is refused
/// ([`crate::RejectReason::ObliqueCylinderCut`]) unless it is recorded as the class's section of
/// the cylinder ([`crate::planes::CylRecord::Sections`]) — drawn on the class's own arrangement,
/// whole or in arcs, which the roads above ask this of before they read a circle; a result face
/// bounded along that section carries an elliptical edge, and is an operand like any other. A seated
/// disk (`LoopRing::Circle`, minted from a **single closed edge** lying in the class) is the class's
/// section whole — a circle on a class across the axis, an ellipse on an oblique one (an operand's
/// slanted cap) — so its readers ask this too before they read a circle.
///
/// **The argument is easy to misread**,
/// which is why a consumer that would answer *silently wrong* asks instead of assuming.
/// The relation calls a **zero** vector parallel to everything (`Across`), so a zero normal or
/// axis answers `true` here — "carries" is the positive reading and a guard spelled `!` fails **open**.
/// Unreachable ([`nacre_topo::CylinderDef::new`] refuses a zero `dir`; the gate refuses a class
/// with no rational description), and the gate asks the same `axis_relation`, so it shares the
/// convention — but the name reads the convention backwards, so it is written down here.
pub(crate) fn class_carries_circle(n: &[nacre_exact::Rat; 3], dir: &[nacre_exact::Rat; 3]) -> bool {
    nacre_exact::axis_relation(n, dir) == nacre_exact::AxisRelation::Across
}

/// **The plane of a point's name that crosses the cylinder's axis** — the first of `planes`
/// whose normal is along the axis ([`class_carries_circle`]). The point lies on that plane and the
/// plane is ⊥ the axis, so its [`axis_param_of_plane`] is the point's height; two such planes are
/// parallel, so one name holds at most one. A plane oblique to the axis meets the axis in a point
/// too, but the named point is not there — which is why a reader that wants the height asks this
/// rather than taking the first plane [`axis_param_of_plane`] answers for. `None` when no plane of
/// the name crosses the axis.
///
/// It takes coefficients, not classes, because its readers hold them from different places: a
/// class's world statement, or a model vertex's plane names.
pub(crate) fn across_axis_plane(
    planes: impl IntoIterator<Item = [nacre_exact::Rat; 4]>,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 4]> {
    let m = def.dir();
    planes
        .into_iter()
        .find(|w| class_carries_circle(&[w[0], w[1], w[2]], &m))
}

/// **The plane of a point's name whose section places the point along its generator** — one not
/// along the axis: a plane across it where the name holds one ([`across_axis_plane`] — a circle
/// stands at one height on every generator, so two such points compare with no generator), else an
/// oblique one, whose section passes the point on its generator too. `None` when every plane of the
/// name runs along the axis.
pub(crate) fn section_plane(
    planes: impl IntoIterator<Item = [nacre_exact::Rat; 4]>,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 4]> {
    let m = def.dir();
    let planes: Vec<[nacre_exact::Rat; 4]> = planes.into_iter().collect();
    across_axis_plane(planes.iter().copied(), def).or_else(|| {
        planes.into_iter().find(|w| {
            nacre_exact::axis_relation(&[w[0], w[1], w[2]], &m)
                == nacre_exact::AxisRelation::Oblique
        })
    })
}

/// **How far along its axis a lateral face reaches** — the hull of the bands its boundary sections
/// sweep ([`trace_band`]), every carrier not along the axis: a circle's band is its own height, an
/// ellipse's the band about its axis meet, rounded outward. `None` where a band cannot be read or the
/// hull is a single height, and the consumer declines the face.
///
/// **A reach, not the face's extent.** It is never narrower than the face — a section bounds the face
/// only within its band — so a clearance proved against it holds; where every boundary section is a
/// circle it is exactly the face's range ([`lateral_t_range`]).
///
/// **Every loop, not the outer one.** A band's upper rim is an inner loop: read off the outer
/// loop alone, the span of a lateral over `[0, 4]`
/// bitten at both ends is `[0, 0.5]`, and the gates prove a cylinder or a slanted wall crossing it
/// at `z = 2` clear of it. A hole's bands lie inside the hull, so they never move its ends.
pub(super) fn lateral_axial_reach(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 2]> {
    let mut hull: Option<[nacre_exact::Rat; 2]> = None;
    for w in lateral_boundary_sections(model, face, def)? {
        let [lo, hi] = trace_band(&w, def)?;
        hull = Some(match hull {
            None => [lo, hi],
            Some([a, b]) => [a.min(lo), b.max(hi)],
        });
    }
    let [lo, hi] = hull?;
    (lo < hi).then_some([lo, hi])
}

/// **The planes whose sections bound a lateral face** — the world coefficients of every boundary
/// carrier not along the axis (a rim's, an arc's: a circle's or an ellipse's), each once, over every
/// loop. `None` where a carrier has no narrow world description.
pub(crate) fn lateral_boundary_sections(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<Vec<[nacre_exact::Rat; 4]>> {
    let m = def.dir();
    let mut out: Vec<[nacre_exact::Rat; 4]> = Vec::new();
    for he in std::iter::once(&face.outer)
        .chain(&face.inner)
        .flat_map(|l| &l.half_edges)
    {
        let coeffs = world_narrow(model, boundary_plane(model, face, he)).ok()?;
        if nacre_exact::axis_relation(&[coeffs[0], coeffs[1], coeffs[2]], &m)
            != nacre_exact::AxisRelation::Along
            && !out.contains(&coeffs)
        {
            out.push(coeffs);
        }
    }
    Some(out)
}

/// The plane on the other side of a lateral face's boundary edge.
fn boundary_plane(
    model: &Model,
    face: &nacre_topo::Face,
    he: &nacre_topo::HalfEdge,
) -> Handle<Surface> {
    let [a, b] = model.edge(he.edge).surfaces;
    if a == face.surface { b } else { a }
}

/// **A lateral face's range along its axis, where every boundary section is a circle** — the
/// lowest and highest of [`lateral_circle_heights`], the heights the face's rims stand at.
/// `None` where those cannot be read, are fewer than two, or a boundary section is not a circle:
/// an ellipse has no one height, and the face's reach is then [`lateral_axial_reach`]'s band.
pub(super) fn lateral_t_range(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 2]> {
    let ts = lateral_circle_heights(model, face, def)?;
    let (lo, hi) = (ts.iter().min()?, ts.iter().max()?);
    (lo < hi).then_some([*lo, *hi])
}

/// **The heights of every circle a lateral face's boundary rides** — its rims and the arcs of its
/// holes and panels, each once, as axis parameters.
///
/// Each such edge's carrier pair is `[lateral, cap-plane]`; the cap plane meets the axis
/// `o + t·m` at `t = −(n·o + d)/(n·m)` — rational whenever the plane has a narrow name (its
/// ⊥-ness guarantees `n·m ≠ 0`). A carrier parallel to the axis is a ruling's wall and gives no
/// height. `None` where a carrier has no world description, the arithmetic runs out, or a carrier
/// is oblique to the axis — its section is an ellipse with no one height, so the face's boundary is
/// not a set of circles at all.
pub(super) fn lateral_circle_heights(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<Vec<nacre_exact::Rat>> {
    use nacre_exact::Rat;
    let m = def.dir();
    let mut ts: Vec<Rat> = Vec::new();
    for he in std::iter::once(&face.outer)
        .chain(&face.inner)
        .flat_map(|l| &l.half_edges)
    {
        let cap = boundary_plane(model, face, he);
        // **The cap's world coefficients** — a plane's name is stated in the frame its own
        // motion names, and this parameter is read against a *world* axis. One whose chain does
        // not fold has no world description and the face declines (`None`), which is the same
        // answer the whole row already gives for a rim with no narrow name.
        let coeffs = world_narrow(model, cap).ok()?;
        // A carrier parallel to the axis is a ruling's wall: it has no station and is not one.
        // An oblique carrier's section has no one height, so the question has no answer.
        match nacre_exact::axis_relation(&[coeffs[0], coeffs[1], coeffs[2]], &m) {
            nacre_exact::AxisRelation::Along => continue,
            nacre_exact::AxisRelation::Oblique => return None,
            nacre_exact::AxisRelation::Across => {}
        }
        let t = axis_param_of_plane(&coeffs, def)?;
        if !ts.contains(&t) {
            ts.push(t);
        }
    }
    Some(ts)
}

/// **The angular extent of a lateral face** — the union of the arcs its outer loop's rim
/// edges trace on the cross-section, as one counter-clockwise arc of radial vectors
/// ([`Theta::Arc`]). A rim edge is one whose other carrier is a cap (a plane ⊥ the axis); its two
/// vertices are the arc's ends in the producer's own order (`derive_edge_curve`: `[A, B]` is A to
/// B counter-clockwise about the axis). A whole rim (`[v, v]`), an irrational corner, a carrier
/// this road cannot translate into the world, or rims that do not chain into one arc give
/// [`Theta::Whole`]. The outer loop is enough: a band's outer
/// loop is a wrapping rim, and that rim's arcs alone close the circle.
///
/// **The extent is never narrower than the face.** A lateral face is bounded on its chart by
/// rulings and arcs of plane sections — a boolean refuses two cylinders whose faces meet
/// (`CylinderPairContact`), so no operand carries anything else — and every plane section not along
/// the axis is a graph over θ, so at every angle the face covers, its upper and lower boundaries
/// are such arcs: the face's angular projection is the union of its arcs'. What the extent does not see is a **hole** in the face, the conservative
/// side. Every reader of a face's angle leans on this — the reaches ([`arc_extent`]), the box, the
/// common cross-section and a line off the arc.
///
/// A corner's radial vector is its point's direction off the axis ([`rim_radial`]): a pierce corner
/// is the meet line's point at its rank ([`nacre_exact::quad::QuadVal::as_rat`] —
/// rational for every wall through or perpendicular to the axis, and for a tangent wall's single
/// root).
pub(super) fn lateral_theta_extent(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Theta {
    rim_arc_union(model, face, def).map_or(Theta::Whole, Theta::Arc)
}

/// [`lateral_theta_extent`]'s union, `None` where it reads the whole circle.
fn rim_arc_union(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<RimArc> {
    use nacre_exact::Rat;
    use nacre_topo::Vertex;
    let m = def.dir();
    let radial =
        |vh: Handle<Vertex>, cap: [Rat; 4]| rim_radial(model, face, def, vh, &cap).ok()?.as_rat();
    let mut acc: Option<RimArc> = None;
    for he in &face.outer.half_edges {
        let coeffs = world_narrow(model, boundary_plane(model, face, he)).ok()?;
        if nacre_exact::axis_relation(&[coeffs[0], coeffs[1], coeffs[2]], &m)
            == nacre_exact::AxisRelation::Along
        {
            continue; // a ruling on a wall, not a rim
        }
        let [va, vb] = model.edge(he.edge).vertices;
        if va == vb {
            return None; // a whole rim: the whole circle
        }
        let arc = RimArc {
            from: radial(va, coeffs)?,
            to: radial(vb, coeffs)?,
        };
        acc = Some(match acc {
            None => arc,
            Some(cur) => {
                let (f_in, t_in) = (
                    arc_contains(&cur, &rat_dir(&arc.from), &m)?,
                    arc_contains(&cur, &rat_dir(&arc.to), &m)?,
                );
                if f_in && t_in {
                    // Inside the current arc — or the two together close the circle.
                    if arc != cur
                        && arc_contains(&arc, &rat_dir(&cur.from), &m)?
                        && arc_contains(&arc, &rat_dir(&cur.to), &m)?
                    {
                        return None;
                    }
                    cur
                } else if f_in {
                    RimArc {
                        from: cur.from,
                        to: arc.to,
                    }
                } else if t_in {
                    RimArc {
                        from: arc.from,
                        to: cur.to,
                    }
                } else if arc_contains(&arc, &rat_dir(&cur.from), &m)?
                    && arc_contains(&arc, &rat_dir(&cur.to), &m)?
                {
                    arc // the current arc lies inside this one
                } else {
                    return None; // two arcs that do not chain: no single extent
                }
            }
        });
    }
    acc
}

/// A radial vector `r₀ + √c·r₁` — rational vectors and one radical. A rim corner's radial vector
/// is rational for a seam vertex and for a pierce corner whose root is; a wall off the axis cuts
/// the rim at an irrational root, and this is how such a corner is held without rounding it.
#[derive(Clone, Copy, Debug)]
pub(super) struct QuadVec {
    r0: [nacre_exact::Rat; 3],
    r1: [nacre_exact::Rat; 3],
    c: nacre_exact::Rat,
}

impl QuadVec {
    /// The vector itself when it is rational (no radical part, or a square radicand).
    pub(super) fn as_rat(&self) -> Option<[nacre_exact::Rat; 3]> {
        let mut out = [nacre_exact::Rat::from_int(0); 3];
        for (k, o) in out.iter_mut().enumerate() {
            *o = nacre_exact::QuadVal::new(self.r0[k], self.r1[k], self.c)?.as_rat()?;
        }
        Some(out)
    }

    /// `(self × v)·m` for a rational `v`.
    fn turn_to(
        &self,
        v: &[nacre_exact::Rat; 3],
        m: &[nacre_exact::Rat; 3],
    ) -> Option<nacre_exact::QuadVal> {
        let a = nacre_exact::dot3_rat(&nacre_exact::cross3_rat(&self.r0, v)?, m)?;
        let b = nacre_exact::dot3_rat(&nacre_exact::cross3_rat(&self.r1, v)?, m)?;
        nacre_exact::QuadVal::new(a, b, self.c)
    }

    /// `self · v` for a rational `v`.
    fn dot(&self, v: &[nacre_exact::Rat; 3]) -> Option<nacre_exact::QuadVal> {
        nacre_exact::QuadVal::new(
            nacre_exact::dot3_rat(&self.r0, v)?,
            nacre_exact::dot3_rat(&self.r1, v)?,
            self.c,
        )
    }

    /// The sign of `f(self, other)` for a bilinear `f` — `(u × v)·m` or `u·v` — over two radicals:
    /// `f(r₀,s₀) + √c·f(r₁,s₀) + √c′·f(r₀,s₁) + √(c·c′)·f(r₁,s₁)`.
    fn bilinear_sign(
        &self,
        other: &QuadVec,
        f: impl Fn(&[nacre_exact::Rat; 3], &[nacre_exact::Rat; 3]) -> Option<nacre_exact::Rat>,
    ) -> Option<nacre_exact::Orient> {
        nacre_exact::biquad_sign(
            f(&self.r0, &other.r0)?,
            f(&self.r1, &other.r0)?,
            f(&self.r0, &other.r1)?,
            f(&self.r1, &other.r1)?,
            self.c,
            other.c,
        )
    }
}

/// **A rim corner's radial vector** — the corner's direction off the axis: the corner minus the axis
/// point on its cap `cap`, less any part along the axis (none on a cap across it; an oblique cap
/// meets the axis away from the corner's own height). A pierce corner is the meet line's point at
/// its rank, held with its radical — the meet [`pierce_on`] solves, so a corner it cannot read
/// fails by its cause. One rule for both readers: [`lateral_theta_extent`] folds it to a rational
/// vector, [`lateral_cover_on_ruling`] reads it as it is.
pub(super) fn rim_radial(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    vh: Handle<nacre_topo::Vertex>,
    cap: &[nacre_exact::Rat; 4],
) -> Result<QuadVec, CornerFail> {
    use nacre_exact::Rat;
    use nacre_topo::Vertex;
    let (o, m) = (def.origin(), def.dir());
    let zero = Rat::from_int(0);
    let width = CornerFail::Width;
    let centre = nacre_exact::axis_plane_meet(cap, &o, &m).ok_or(width)?;
    match *model.vertex(vh) {
        // A seam vertex ends only a whole rim (`[v, v]`), which both readers answer as the whole
        // circle before asking for a corner.
        Vertex::OnSeam(_) => Err(CornerFail::Shape),
        Vertex::Pierce {
            planes,
            cylinder,
            rank,
        } => {
            if cylinder != face.surface {
                return Err(CornerFail::Shape);
            }
            // `def` is this face's own cylinder's world statement.
            let (line, s) = pierce_on(model, planes, def, rank)?;
            let (b, d) = (line.base(), line.dir());
            let (mut r0, mut r1) = ([zero; 3], [zero; 3]);
            for k in 0..3 {
                r0[k] = b[k]
                    .checked_add(s.a().checked_mul(d[k]).ok_or(width)?)
                    .and_then(|x| x.checked_sub(centre[k]))
                    .ok_or(width)?;
                r1[k] = s.b().checked_mul(d[k]).ok_or(width)?;
            }
            // Both parts off the axis: the projection is linear, and leaves a vector across the
            // axis as it is.
            Ok(QuadVec {
                r0: nacre_exact::off_axis_rat(&r0, &m).ok_or(width)?,
                r1: nacre_exact::off_axis_rat(&r1, &m).ok_or(width)?,
                c: s.c(),
            })
        }
        // A boolean refuses an operand with a two-cylinder edge before the gate reads its rims
        // (`CylinderPairEdge`); unread, the corner answers as one the rule cannot state.
        Vertex::ThreePlane(_) | Vertex::CylinderPair { .. } => Err(CornerFail::Shape),
    }
}

/// Whether the rim arc `from → to` (counter-clockwise about `m`) holds the directions just beside
/// the rational radial direction `x` — counter-clockwise of it (`ccw`), the arc **half-open** with
/// `from` in and `to` out, or clockwise of it, `from` out and `to` in. Half-open so that where two
/// arcs of one rim meet on the line — at the node between them — the rim is counted once, and an
/// arc that ends on the line is counted on the side it lies. [`arc_contains`]'s three-way reading,
/// over [`QuadVec`] ends.
fn arc_holds_beside(
    from: &QuadVec,
    to: &QuadVec,
    x: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    ccw: bool,
) -> Option<bool> {
    use nacre_exact::Orient;
    let turn = |u: &[nacre_exact::Rat; 3], v: &[nacre_exact::Rat; 3]| {
        nacre_exact::dot3_rat(&nacre_exact::cross3_rat(u, v)?, m)
    };
    let ft = from.bilinear_sign(to, turn)?;
    let fx = from.turn_to(x, m)?.sign();
    // `(x × to)·m = −(to × x)·m`.
    let xt = match to.turn_to(x, m)?.sign() {
        Orient::Positive => Orient::Negative,
        Orient::Negative => Orient::Positive,
        Orient::Zero => Orient::Zero,
    };
    let held = match ft {
        Orient::Positive => fx != Orient::Negative && xt != Orient::Negative,
        Orient::Negative => !(xt == Orient::Negative && fx == Orient::Negative),
        Orient::Zero => {
            from.bilinear_sign(to, nacre_exact::dot3_rat)? == Orient::Positive
                || fx != Orient::Negative
        }
    };
    let at_end = |end: &QuadVec, turn: Orient| -> Option<bool> {
        Some(turn == Orient::Zero && end.dot(x)?.sign() == Orient::Positive)
    };
    let out = if ccw {
        at_end(to, xt)?
    } else {
        at_end(from, fx)?
    };
    Some(held && !out)
}

/// **Where a lateral face lies on one ruling** — the axis-parameter intervals the ruling through
/// `foot` (a point on the cylinder) shares with the face, in the parameter
/// [`axis_param_of_plane`] gives a station.
///
/// On its chart a lateral face is bounded by arcs of plane sections (each a graph over `θ`) and
/// rulings (`θ` fixed), so the vertical line `θ = θ₀` meets its boundary only where an arc's
/// angular range holds `θ₀` — at the height its section stands on that generator: those stations,
/// sorted, alternate in and out, and paired they are the face's stretches of the line. Every loop counts — a window is an inner loop, as the rim
/// opposite the outer loop is — which is what the face's angular extent and axial span, a
/// bounding rectangle, cannot say.
///
/// `None` is a face this cannot read, and the caller keeps its conservative answer: a corner off
/// the rational seam, a ruling edge on the line itself (its
/// plane holds the line — the face ends there, which is the gate's own question, not this one's),
/// or stations that do not pair.
pub(super) fn lateral_cover_on_ruling(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    foot: &[nacre_exact::Rat; 3],
) -> Option<Vec<[nacre_exact::Rat; 2]>> {
    lateral_stations_on_ruling(model, face, def, foot, None).ok()
}

/// **Where a lateral face lies on each side of one ruling** — `[ccw, cw]`: the axis-parameter
/// intervals over which the face holds the generators just counter-clockwise of the ruling through
/// `foot` (about the axis), and just clockwise of it. Both sides together are where the face runs
/// through the line; one side alone is where the line is the face's edge. The stations are
/// [`lateral_cover_on_ruling`]'s, each arc counted on the side it lies ([`arc_holds_beside`]); a
/// ruling edge on the line is the face's boundary there and a seam on it is no boundary at all, so
/// neither is a station. `Err` is a face this cannot read, by its cause
/// ([`lateral_stations_on_ruling`]).
pub(super) fn lateral_sides_on_ruling(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    foot: &[nacre_exact::Rat; 3],
) -> Result<[Vec<[nacre_exact::Rat; 2]>; 2], RejectReason> {
    Ok([
        lateral_stations_on_ruling(model, face, def, foot, Some(true))?,
        lateral_stations_on_ruling(model, face, def, foot, Some(false))?,
    ])
}

/// The stations both readings above pair: `beside` `None` is the cover (counter-clockwise
/// half-open, and a ruling edge on the line unread), `Some(ccw)` one side of the line.
///
/// `Err` names why the face could not be read, as the gate names everything it forms: a value
/// `Rat` does not hold is width ([`RejectReason::WitnessNotRational`]); a neighbour with no world
/// statement, or a corner [`rim_radial`] cannot state, by the corner reader's cause
/// ([`CornerFail::stated_reason`]); and a shape this reading does not spell — an edge to a surface
/// that is not a plane (in the cover, the face's own seam), a corner it cannot spell, stations that
/// do not pair, the cover's ruling edge on the line — is the lateral the gate could not read,
/// [`RejectReason::CylinderGateUndecided`].
fn lateral_stations_on_ruling(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    foot: &[nacre_exact::Rat; 3],
    beside: Option<bool>,
) -> Result<Vec<[nacre_exact::Rat; 2]>, RejectReason> {
    use nacre_exact::Rat;
    let width = RejectReason::WitnessNotRational;
    let unread = RejectReason::CylinderGateUndecided;
    let (o, m) = (def.origin(), def.dir());
    let zero = Rat::from_int(0);
    // The ruling's radial direction: `foot − o` off the axis.
    let mut w = [zero; 3];
    for k in 0..3 {
        w[k] = foot[k].checked_sub(o[k]).ok_or(width)?;
    }
    let x = nacre_exact::off_axis_rat(&w, &m).ok_or(width)?;
    // The generator's point at the axis origin's height, `o + x`: a section's height on this
    // generator is where the plane meets it — for a plane across the axis, its height everywhere.
    let mut g = o;
    for k in 0..3 {
        g[k] = o[k].checked_add(x[k]).ok_or(width)?;
    }
    let mut stations: Vec<Rat> = Vec::new();
    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
        for he in &lp.half_edges {
            let e = model.edge(he.edge);
            let [a, b] = e.surfaces;
            let other = if a == face.surface { b } else { a };
            if other == face.surface && beside.is_some() {
                continue; // a seam: the face runs on across it
            }
            if !matches!(model.surface(other), nacre_topo::Surface::Plane { .. }) {
                return Err(unread);
            }
            let coeffs =
                world_narrow(model, other).map_err(|f| f.stated_reason().unwrap_or(unread))?;
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            if nacre_exact::axis_relation(&n, &m) == nacre_exact::AxisRelation::Along {
                // A ruling edge: off the line it never meets it; on it, the face ends there.
                let side = nacre_exact::dot3_rat(&n, foot)
                    .and_then(|v| v.checked_add(coeffs[3]))
                    .ok_or(width)?;
                if side == zero && beside.is_none() {
                    return Err(unread);
                }
            } else {
                let [va, vb] = e.vertices;
                let radial = |vh| {
                    rim_radial(model, face, def, vh, &coeffs)
                        .map_err(|f| f.stated_reason().unwrap_or(unread))
                };
                let held = va == vb
                    || arc_holds_beside(&radial(va)?, &radial(vb)?, &x, &m, beside.unwrap_or(true))
                        .ok_or(width)?;
                if held {
                    stations.push(nacre_exact::axis_param_of_plane(&coeffs, &g, &m).ok_or(width)?);
                }
            }
        }
    }
    if stations.len() % 2 != 0 {
        return Err(unread);
    }
    stations.sort_unstable();
    Ok(stations.chunks(2).map(|p| [p[0], p[1]]).collect())
}

/// **How far an arc reaches either side of its centre, along `d`** — the rule for a
/// cylinder's lateral face, said once so a **planar** face's arc piece can ask it too.
///
/// `(lo_off, rho2_lo, hi_off, rho2_hi)`: the arc occupies
/// `[d·centre + lo_off − √rho2_lo, d·centre + hi_off + √rho2_hi]`. The radial term `r·(d·û)` peaks
/// at the direction of `d⊥` when the arc holds it — then the end is `√ρ²` with no offset — and at
/// an **end of the arc** otherwise, where it is `d·v` for that end's radial vector, a value folded
/// into the offset with a zero radical. One root at most on each end, and no new arithmetic.
///
/// **`axis` is the arc's own carrier axis, never a face's normal.** `arc.from → arc.to` is
/// counter-clockwise **about that axis** (`derive_edge_curve`'s convention, which
/// [`lateral_theta_extent`] states), so handing in a normal that runs the other way would silently
/// name the **complementary** arc — not a bound on this one but a different set, which would prove
/// clearances that are not there.
///
/// [`Theta::Whole`] is the whole circle, [`Theta::At`] one point of it (`d·x` both ways, no
/// radical), and `d ∥ axis` leaves no radial term at all.
pub(super) fn arc_extent(
    theta: &Theta,
    r2: &nacre_exact::BigRat,
    axis: &[nacre_exact::Rat; 3],
    d: &[nacre_exact::Rat; 3],
) -> Option<(
    nacre_exact::Rat,
    nacre_exact::Rat,
    nacre_exact::Rat,
    nacre_exact::Rat,
)> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    // `d⊥`, `d` off the axis: the direction of the radial term's peak.
    let dperp = nacre_exact::off_axis_rat(d, axis)?;
    let dperp2 = nacre_exact::dot3_rat(&dperp, &dperp)?;
    // The radial peak `ρ² = r²|d⊥|²`, formed only by a side that reaches it: a wide square (`r2`
    // past `Rat`) declines there, and an end an arc stops at — or one direction — needs no radius.
    let rho2 = || r2.narrow()?.checked_mul(dperp2);
    match theta {
        _ if dperp2 == zero => Some((zero, zero, zero, zero)), // `d ∥ axis`
        Theta::Whole => {
            let rho2 = rho2()?;
            Some((zero, rho2, zero, rho2))
        }
        Theta::At(x) => {
            let dx = nacre_exact::dot3_rat(d, x)?;
            Some((dx, zero, dx, zero))
        }
        Theta::Arc(arc) => {
            let (f, t) = (
                nacre_exact::dot3_rat(d, &arc.from)?,
                nacre_exact::dot3_rat(d, &arc.to)?,
            );
            let (hi_off, hi_rad) = if arc_contains(arc, &rat_dir(&dperp), axis)? {
                (zero, rho2()?)
            } else {
                (f.max(t), zero)
            };
            let mut neg = dperp;
            for x in neg.iter_mut() {
                *x = zero.checked_sub(*x)?;
            }
            let (lo_off, lo_rad) = if arc_contains(arc, &rat_dir(&neg), axis)? {
                (zero, rho2()?)
            } else {
                (f.min(t), zero)
            };
            Some((lo_off, lo_rad, hi_off, hi_rad))
        }
    }
}

#[cfg(test)]
#[path = "../tests/cyl_geom.rs"]
mod tests;
