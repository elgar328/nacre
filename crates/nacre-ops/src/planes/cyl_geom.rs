use super::*;
/// **A plane's exact rational coefficients in the world** — the narrow projection of
/// [`nacre_topo::Model::world_plane_name`], which is the one place the rule lives.
///
/// ★★★ **Not a restatement of the row.** A moved plane's `rotated` flag is the licence the exact
/// shortcuts read before taking the plane's *name* as a world statement, and a moved plane's name
/// speaks its pre-motion frame — so the flag stays set whatever this answers, and the world
/// description rides *beside* it. Clearing the flag once sent a merge through an exact-`f64`
/// triangle test on the face corners, and two unit cubes shifted by `7/11` and `18/11`, sharing a
/// wall whose two `f64` images differ in the last place, came back as two bodies; the class merge
/// reads world names now, never those corners.
pub(crate) fn world_plane_coeffs(
    model: &Model,
    surf: Handle<Surface>,
) -> Option<[nacre_exact::Rat; 4]> {
    model.world_plane_name(surf)?.narrow().copied()
}

/// **Does an increasing axis parameter move toward this class's "above"?** — where "above" is the
/// side the class's plane faces, which is the frame every cell label is written in.
///
/// Exact: the world name against the axis (integers, no overflow), turned by the name's sense.
/// It takes the class's [`WorldName`], so a class without one cannot ask — every reader already
/// stands on a class whose world coefficients it read.
///
/// ★★ **Ask this, do not re-derive it from the class's rational name alone.** The name is
/// *canonical* (first nonzero component positive), which points the other way from the plane on
/// half the classes; a rule spelled against it reads "above" backwards exactly there. That mistake,
/// made while adding the second consumer below, turned 36 tests red at once.
///
/// This is only asked of a class ⊥ to the axis, so the product never vanishes.
pub(crate) fn plus_t_is_above(world: &WorldName, def: &nacre_topo::CylinderDef) -> bool {
    let [a, b, c, _] = world.name.coeff_ints();
    let along = nacre_exact::normal_sense(&[a, b, c], def.dir());
    debug_assert_ne!(
        along,
        nacre_exact::Orient::Zero,
        "asked of a class that does not cross the axis"
    );
    (along == nacre_exact::Orient::Positive) == (world.sense == nacre_topo::Orientation::Forward)
}

/// **Where a plane meets a cylinder's axis, as the axis parameter `t`** —
/// [`nacre_exact::axis_param_of_plane`] read on the statement's own origin and direction, so the
/// formula has one spelling for every reader; the readers that want the point itself (a rim's
/// centre) ask [`nacre_exact::axis_plane_meet`].
///
/// `None` when the plane is parallel to the axis (`n·m = 0`, no meeting point) or when the
/// checked `Rat` arithmetic overflows. ★ A rule that lives inlined in one place while a second
/// site spells a reduced version of it is this repo's dominant defect shape, so every reader asks
/// this.
pub(crate) fn axis_param_of_plane(
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_exact::Rat> {
    nacre_exact::axis_param_of_plane(coeffs, &def.origin(), &def.dir())
}

/// A lateral face's axis-parameter span, read off its rim carrier planes.
///
/// Each rim edge's carrier pair is `[lateral, cap-plane]`; the cap plane meets the axis
/// `o + t·m` at `t = −(n·o + d)/(n·m)` — rational whenever the plane has a narrow name (its
/// ⊥-ness guarantees `n·m ≠ 0`). Two distinct rim planes give the span; anything else (a
/// nameless rim carrier, a non-⊥ rim, fewer or more than two distinct rims) answers `None`
/// and the consumer declines the face.
///
/// ★ **Every loop, not the outer one.** A band's upper rim is an inner loop: read off the outer
/// loop alone, the span of a lateral over `[0, 4]`
/// bitten at both ends is `[0, 0.5]`, and the gates prove a cylinder or a slanted wall crossing it
/// at `z = 2` clear of it. A hole's stations lie inside the span, so they never move its ends.
pub(super) fn lateral_t_range(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 2]> {
    use nacre_exact::Rat;
    let m = def.dir();
    let mut ts: Vec<Rat> = Vec::new();
    for he in std::iter::once(&face.outer)
        .chain(&face.inner)
        .flat_map(|l| &l.half_edges)
    {
        let e = model.edge(he.edge);
        let [a, b] = e.surfaces;
        let cap = if a == face.surface { b } else { a };
        // ★ **The cap's world coefficients** — a plane's name is stated in the frame its own
        // motion names, and this parameter is read against a *world* axis. One whose chain does
        // not fold has no world description and the face declines (`None`), which is the same
        // answer the whole row already gives for a rim with no narrow name.
        let coeffs = *model.world_plane_name(cap)?.narrow()?;
        // A carrier parallel to the axis is a ruling's wall: it has no station and is not one.
        // Only a ⊥ carrier — an arc's cap plane — speaks here; `None` past that is overflow.
        if !nacre_exact::parallel_rat(&[coeffs[0], coeffs[1], coeffs[2]], &m) {
            continue;
        }
        let t = axis_param_of_plane(&coeffs, def)?;
        if !ts.contains(&t) {
            ts.push(t);
        }
    }
    let (lo, hi) = (ts.iter().min()?, ts.iter().max()?);
    (lo < hi).then_some([*lo, *hi])
}

/// All shells of a solid — outer first, then cavities. The boolean seam
/// front-end walks these so a cavitied operand's void walls are seen;
/// a non-hollow solid yields just its outer shell, unchanged.
/// **The angular extent of a lateral face** — the union of the arcs its outer loop's rim
/// edges trace on the cross-section, as one counter-clockwise arc of radial vectors
/// ([`RimArc`]). A rim edge is one whose other carrier is a cap (a plane ⊥ the axis); its two
/// vertices are the arc's ends in the producer's own order (`derive_edge_curve`: `[A, B]` is A to
/// B counter-clockwise about the axis). A whole rim (`[v, v]`), an irrational corner, a carrier
/// this road cannot translate into the world, or rims that do not chain into one arc give `None`
/// — the whole circle, the reading before this existed. The outer loop is enough: a band's outer
/// loop is a wrapping rim, and that rim's arcs alone close the circle.
///
/// A corner's radial vector is its point minus the axis point on the cap: a pierce corner is the
/// meet line's point at its root ([`nacre_exact::quad::QuadVal::as_rat`] —
/// rational for every wall through or perpendicular to the axis, and for a tangent wall's single
/// root).
pub(super) fn lateral_theta_extent(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<RimArc> {
    use nacre_exact::Rat;
    use nacre_topo::Vertex;
    let m = def.dir();
    let world_coeffs = |plane: Handle<Surface>| -> Option<[Rat; 4]> {
        model.world_plane_name(plane)?.narrow().copied()
    };
    let radial =
        |vh: Handle<Vertex>, cap: [Rat; 4]| rim_radial(model, face, def, vh, &cap)?.as_rat();
    let mut acc: Option<RimArc> = None;
    for he in &face.outer.half_edges {
        let e = model.edge(he.edge);
        let [a, b] = e.surfaces;
        let cap = if a == face.surface { b } else { a };
        let coeffs = world_coeffs(cap)?;
        if !nacre_exact::parallel_rat(&[coeffs[0], coeffs[1], coeffs[2]], &m) {
            continue; // a ruling on a wall, not a rim
        }
        let [va, vb] = e.vertices;
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
                    arc_contains(&cur, &arc.from, &m)?,
                    arc_contains(&cur, &arc.to, &m)?,
                );
                if f_in && t_in {
                    // Inside the current arc — or the two together close the circle.
                    if arc != cur
                        && arc_contains(&arc, &cur.from, &m)?
                        && arc_contains(&arc, &cur.to, &m)?
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
                } else if arc_contains(&arc, &cur.from, &m)? && arc_contains(&arc, &cur.to, &m)? {
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

/// **A rim corner's radial vector** — the corner minus the axis point on its cap `cap`: a pierce
/// corner is the meet line's point at its root, held with its radical. One rule for both readers:
/// [`lateral_theta_extent`] folds it to a rational vector, [`lateral_cover_on_ruling`] reads it as
/// it is.
pub(super) fn rim_radial(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    vh: Handle<nacre_topo::Vertex>,
    cap: &[nacre_exact::Rat; 4],
) -> Option<QuadVec> {
    use nacre_exact::Rat;
    use nacre_exact::quad::CylinderMeet;
    use nacre_topo::{QuadRoot, Vertex};
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let zero = Rat::from_int(0);
    let centre = nacre_exact::axis_plane_meet(cap, &o, &m)?;
    match *model.vertex(vh) {
        // A seam vertex ends only a whole rim (`[v, v]`), which both readers answer as the whole
        // circle before asking for a corner.
        Vertex::OnSeam(_) => None,
        Vertex::Pierce {
            planes,
            cylinder,
            root,
        } => {
            if cylinder != face.surface {
                return None;
            }
            // `planes` are stored in ascending-handle order and the roots run along `n₀ × n₁`
            // (`Vertex::Pierce`'s convention); the world statement translates only the constants,
            // so the normals — and the order — are the stored ones.
            let (c0, c1) = (
                world_plane_coeffs(model, planes[0])?,
                world_plane_coeffs(model, planes[1])?,
            );
            let (line, s) = match nacre_exact::quad::plane_plane_cylinder(&c0, &c1, &o, &m, r2)? {
                CylinderMeet::Pair { line, s } => match root {
                    QuadRoot::Lo => (line, s[0]),
                    QuadRoot::Hi => (line, s[1]),
                    QuadRoot::Double => return None,
                },
                CylinderMeet::Tangent { line, s } => match root {
                    QuadRoot::Double => (line, nacre_exact::QuadVal::from_rat(s)),
                    _ => return None,
                },
                _ => return None,
            };
            let (b, d) = (line.base(), line.dir());
            let (mut r0, mut r1) = ([zero; 3], [zero; 3]);
            for k in 0..3 {
                r0[k] = b[k]
                    .checked_add(s.a().checked_mul(d[k])?)?
                    .checked_sub(centre[k])?;
                r1[k] = s.b().checked_mul(d[k])?;
            }
            Some(QuadVec { r0, r1, c: s.c() })
        }
        Vertex::ThreePlane(_) => None,
    }
}

/// Whether the rim arc `from → to` (counter-clockwise about `m`) holds the rational radial
/// direction `x`, **half-open**: `from` in, `to` out. Half-open so that where the arcs of one rim
/// meet on the line — at a seam vertex, which splits a cut circle into two arcs there — the rim
/// is counted once. [`arc_contains`]'s three-way reading, over [`QuadVec`] ends.
fn arc_holds_half_open(
    from: &QuadVec,
    to: &QuadVec,
    x: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
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
    let at_to = xt == Orient::Zero && to.dot(x)?.sign() == Orient::Positive;
    Some(held && !at_to)
}

/// **Where a lateral face lies on one ruling** — the axis-parameter intervals the ruling through
/// `foot` (a point on the cylinder) shares with the face, in the parameter
/// [`axis_param_of_plane`] gives a station.
///
/// On its chart a lateral face is bounded by arcs (a rim on a cap ⊥ the axis, `t` fixed) and
/// rulings (`θ` fixed), so the vertical line `θ = θ₀` meets its boundary only where an arc's
/// angular range holds `θ₀`: those stations, sorted, alternate in and out, and paired they are
/// the face's stretches of the line. Every loop counts — a window is an inner loop, or a notch in
/// the outer one where the seam runs through it — which is what the face's angular extent and
/// axial span, a bounding rectangle, cannot say.
///
/// `None` is a face this cannot read, and the caller keeps its conservative answer: a corner off
/// the rational seam, a carrier neither ⊥ nor ∥ the axis, a ruling edge on the line itself (its
/// plane holds the line — the face ends there, which is the gate's own question, not this one's),
/// or stations that do not pair.
pub(super) fn lateral_cover_on_ruling(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
    foot: &[nacre_exact::Rat; 3],
) -> Option<Vec<[nacre_exact::Rat; 2]>> {
    use nacre_exact::Rat;
    let (o, m) = (def.origin(), def.dir());
    let zero = Rat::from_int(0);
    // The ruling's radial direction: `foot − o` less its component along the axis.
    let mut w = [zero; 3];
    for k in 0..3 {
        w[k] = foot[k].checked_sub(o[k])?;
    }
    let mm = nacre_exact::dot3_rat(&m, &m)?;
    let q = nacre_exact::dot3_rat(&w, &m)?.checked_mul(Rat::new(mm.denom(), mm.numer())?)?;
    let mut x = [zero; 3];
    for k in 0..3 {
        x[k] = w[k].checked_sub(q.checked_mul(m[k])?)?;
    }
    let mut stations: Vec<Rat> = Vec::new();
    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
        for he in &lp.half_edges {
            let e = model.edge(he.edge);
            let [a, b] = e.surfaces;
            let other = if a == face.surface { b } else { a };
            let coeffs = world_plane_coeffs(model, other)?;
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            if nacre_exact::parallel_rat(&n, &m) {
                let [va, vb] = e.vertices;
                let held = va == vb
                    || arc_holds_half_open(
                        &rim_radial(model, face, def, va, &coeffs)?,
                        &rim_radial(model, face, def, vb, &coeffs)?,
                        &x,
                        &m,
                    )?;
                if held {
                    stations.push(axis_param_of_plane(&coeffs, def)?);
                }
            } else if nacre_exact::dot3_rat(&n, &m)? == zero {
                // A ruling edge: off the line it never meets it; on it, the face ends there.
                let side = nacre_exact::dot3_rat(&n, foot)?.checked_add(coeffs[3])?;
                if side == zero {
                    return None;
                }
            } else {
                return None;
            }
        }
    }
    if stations.len() % 2 != 0 {
        return None;
    }
    stations.sort_unstable();
    Some(stations.chunks(2).map(|p| [p[0], p[1]]).collect())
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
/// ★★ **`axis` is the arc's own carrier axis, never a face's normal.** `arc.from → arc.to` is
/// counter-clockwise **about that axis** (`derive_edge_curve`'s convention, which
/// [`lateral_theta_extent`] states), so handing in a normal that runs the other way would silently
/// name the **complementary** arc — not a bound on this one but a different set, which would prove
/// clearances that are not there.
///
/// `arc = None` is the whole circle, and `d ∥ axis` leaves no radial term at all.
pub(super) fn arc_extent(
    arc: Option<&RimArc>,
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
    let (dm, mm) = (
        nacre_exact::dot3_rat(d, axis)?,
        nacre_exact::dot3_rat(axis, axis)?,
    );
    // `d⊥ = d − (d·m / m·m) m`, the direction of the radial term's peak; `|d⊥|² = d·d − (d·m)²/m·m`.
    let k = dm.checked_mul(Rat::new(mm.denom(), mm.numer())?)?;
    let mut dperp = *d;
    for i in 0..3 {
        dperp[i] = dperp[i].checked_sub(k.checked_mul(axis[i])?)?;
    }
    let dperp2 = nacre_exact::dot3_rat(&dperp, &dperp)?;
    // A wide square (`r2` past `Rat`) declines here.
    let rho2 = r2.narrow()?.checked_mul(dperp2)?;
    match arc {
        _ if dperp2 == zero => Some((zero, zero, zero, zero)), // `d ∥ axis`
        None => Some((zero, rho2, zero, rho2)),
        Some(arc) => {
            let (f, t) = (
                nacre_exact::dot3_rat(d, &arc.from)?,
                nacre_exact::dot3_rat(d, &arc.to)?,
            );
            let (hi_off, hi_rad) = if arc_contains(arc, &dperp, axis)? {
                (zero, rho2)
            } else {
                (f.max(t), zero)
            };
            let mut neg = dperp;
            for x in neg.iter_mut() {
                *x = zero.checked_sub(*x)?;
            }
            let (lo_off, lo_rad) = if arc_contains(arc, &neg, axis)? {
                (zero, rho2)
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
