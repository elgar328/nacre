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

/// A lateral face's axis-parameter span, read off its rim carrier planes.
///
/// Each rim edge's carrier pair is `[lateral, cap-plane]`; the cap plane meets the axis
/// `o + t·m` at `t = −(n·o + d)/(n·m)` — rational whenever the plane has a narrow name (its
/// ⊥-ness guarantees `n·m ≠ 0`). Two distinct rim planes give the span; anything else (a
/// nameless rim carrier, a non-⊥ rim, fewer or more than two distinct rims) answers `None`
/// and the consumer declines the face.
/// **Where a plane meets a cylinder's axis, as the axis parameter `t`** — the one spelling of
/// `t = −(n·o + d)/(n·m)` for `axis(t) = o + t·m`.
///
/// `None` when the plane is parallel to the axis (`n·m = 0`, no meeting point) or when the
/// checked `Rat` arithmetic overflows. ★ Three consumers ask this question — the lateral's rim
/// span, the transversal-circle test, and the band pass — and a rule that lives inlined in one
/// place while a second site spells a reduced version of it is this repo's dominant defect
/// shape, so it lives here once.
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

pub(crate) fn axis_param_of_plane(
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<nacre_exact::Rat> {
    use nacre_exact::Rat;
    let (o, m) = (def.origin(), def.dir());
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let nm = nacre_exact::dot3_rat(&n, &m)?;
    if nm == Rat::from_int(0) {
        return None;
    }
    let no_d = nacre_exact::dot3_rat(&n, &o)?.checked_add(coeffs[3])?;
    Rat::from_int(0)
        .checked_sub(no_d)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)
}

pub(super) fn lateral_t_range(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 2]> {
    use nacre_exact::Rat;
    let m = def.dir();
    let mut ts: Vec<Rat> = Vec::new();
    for he in &face.outer.half_edges {
        let e = model.edge(he.edge);
        let [a, b] = e.surfaces;
        let cap = if a == face.surface { b } else { a };
        if cap == face.surface {
            continue; // a slit edge is self-adjacent — not a carrier
        }
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
/// — the whole circle, the reading before this existed.
///
/// A corner's radial vector is its point minus the axis point on the cap: a seam vertex is
/// `r·ê` for `ê` the reference direction's unit part ⊥ the axis (rational when the norm is), a
/// pierce corner is the meet line's point at its root ([`nacre_exact::quad::QuadVal::as_rat`] —
/// rational for every wall through or perpendicular to the axis, and for a tangent wall's single
/// root).
pub(super) fn lateral_theta_extent(
    model: &Model,
    face: &nacre_topo::Face,
    def: &nacre_topo::CylinderDef,
) -> Option<RimArc> {
    use nacre_exact::Rat;
    use nacre_exact::quad::CylinderMeet;
    use nacre_topo::{QuadRoot, Vertex};
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let world_coeffs = |plane: Handle<Surface>| -> Option<[Rat; 4]> {
        model.world_plane_name(plane)?.narrow().copied()
    };
    let sub = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            a[0].checked_sub(b[0])?,
            a[1].checked_sub(b[1])?,
            a[2].checked_sub(b[2])?,
        ])
    };
    let scaled = |v: &[Rat; 3], k: Rat| -> Option<[Rat; 3]> {
        Some([
            v[0].checked_mul(k)?,
            v[1].checked_mul(k)?,
            v[2].checked_mul(k)?,
        ])
    };
    // The radial vector of a rim corner on the cap `cap`.
    let radial = |vh: Handle<Vertex>, cap: [Rat; 4]| -> Option<[Rat; 3]> {
        let t = axis_param_of_plane(&cap, def)?;
        let mut centre = o;
        for k in 0..3 {
            centre[k] = centre[k].checked_add(t.checked_mul(m[k])?)?;
        }
        match *model.vertex(vh) {
            Vertex::OnSeam(_) => {
                let e = def.ref_dir();
                let (mm, em) = (
                    nacre_exact::dot3_rat(&m, &m)?,
                    nacre_exact::dot3_rat(&e, &m)?,
                );
                let mut e1 = [Rat::from_int(0); 3];
                for k in 0..3 {
                    e1[k] = mm.checked_mul(e[k])?.checked_sub(em.checked_mul(m[k])?)?;
                }
                // `r/|e₁|` read as `√(r²/|e₁|²)`: the same rational when both roots are, and a
                // rational where neither is alone (`r = √2` on `|e₁| = √2`).
                let ee = nacre_exact::dot3_rat(&e1, &e1)?;
                scaled(
                    &e1,
                    nacre_exact::rat_sqrt_exact_big(
                        &r2.mul_rat(Rat::new(ee.denom(), ee.numer())?),
                    )?,
                )
            }
            Vertex::Pierce {
                planes,
                cylinder,
                root,
            } => {
                if cylinder != face.surface {
                    return None;
                }
                // `planes` are stored in ascending-handle order and the roots run along
                // `n₀ × n₁` (`Vertex::Pierce`'s convention); the world statement translates
                // only the constants, so the normals — and the order — are the stored ones.
                let (c0, c1) = (world_coeffs(planes[0])?, world_coeffs(planes[1])?);
                let (line, sv) =
                    match nacre_exact::quad::plane_plane_cylinder(&c0, &c1, &o, &m, r2)? {
                        CylinderMeet::Pair { line, s } => match root {
                            QuadRoot::Lo => (line, s[0].as_rat()?),
                            QuadRoot::Hi => (line, s[1].as_rat()?),
                            QuadRoot::Double => return None,
                        },
                        CylinderMeet::Tangent { line, s } => match root {
                            QuadRoot::Double => (line, s),
                            _ => return None,
                        },
                        _ => return None,
                    };
                let (b, d) = (line.base(), line.dir());
                let mut p = b;
                for k in 0..3 {
                    p[k] = p[k].checked_add(sv.checked_mul(d[k])?)?;
                }
                sub(&p, &centre)
            }
            Vertex::ThreePlane(_) => None,
        }
    };
    let mut acc: Option<RimArc> = None;
    for he in &face.outer.half_edges {
        let e = model.edge(he.edge);
        let [a, b] = e.surfaces;
        let cap = if a == face.surface { b } else { a };
        if cap == face.surface {
            continue; // the seam
        }
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
