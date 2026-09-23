use super::*;
/// **A plane's exact rational coefficients in the world** — the narrow projection of
/// [`nacre_topo::Model::world_plane_name`], which is the one place the rule lives.
///
/// ★★★ **Not a restatement of the row.** A moved plane's `rotated` flag says two things at once,
/// and only one of them is about descriptions: it also says the row's `tri` is a *realized*
/// triangle rather than the truth, which is what routes `Judge::planes_coplanar` to the
/// high-precision road. Measured, by breaking it: two unit cubes shifted by `7/11` and `18/11`
/// share a wall whose two f64 images differ in the last place, and flipping such a plane to
/// unrotated sent the merge through the exact-f64 triangle test — one body came back as two. So
/// the world description rides *beside* the flag, and the judging road is left alone.
pub(crate) fn world_plane_coeffs(
    model: &Model,
    surf: Handle<Surface>,
) -> Option<[nacre_exact::Rat; 4]> {
    model.world_plane_name(surf)?.narrow().copied()
}

/// **A cylinder's exact statement in the world** — the one door between a cylinder's truth
/// (written in the frame its motion names) and every consumer that compares it against world
/// planes: the population gate's clearance arithmetic, the arrangement's circles and rulings,
/// the band pass.
///
/// Unmoved: the statement itself. Moved by a chain that is a pure rational translation
/// ([`nacre_topo::Model::chain_translation`]): the same statement with its origin shifted —
/// exact, because a rational translation maps a rational statement to a rational one, and the
/// axis direction, the seam reference and the radius are all invariant under it. Anything else
/// (a rotation, a frame, overflow): `None`, and the caller refuses rather than measuring across
/// two frames.
///
/// ★ **The postcondition is checked, not assumed.** The surface's `f64` cache is already the
/// *realized* world cylinder, so it is an independent second description of the very thing this
/// function claims to produce — a fold with the wrong sign, or one that walked the chain the
/// wrong way, disagrees with it by twice the offset. Consumer-side net, in the shape this kernel
/// keeps arriving at: check the postcondition rather than trusting the derivation.
pub(crate) fn world_cylinder_def(
    model: &Model,
    surf: Handle<Surface>,
) -> Option<nacre_topo::CylinderDef> {
    let nacre_topo::Surface::Cylinder { def, motion } = model.surface(surf) else {
        unreachable!("a cylinder surface carries a cylinder truth")
    };
    let out = match motion {
        None => def.clone(),
        Some(leaf) => {
            let t = model.chain_translation(*leaf)?;
            let mut o = def.origin();
            for (c, d) in o.iter_mut().zip(t) {
                *c = c.checked_add(d)?;
            }
            // The three invariants of a translation, so `new`'s checks cannot newly fail here —
            // it is called rather than bypassed because the type's constructor is the only way in.
            nacre_topo::CylinderDef::new(o, def.dir(), def.ref_dir(), def.r2().clone())?
        }
    };
    debug_assert!(
        {
            let nacre_geom::Surface::Cylinder(cache) = model.surface_cache(surf) else {
                unreachable!(
                    "push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache"
                )
            };
            let o = Point3::from_array(out.origin().map(|r| r.to_f64()));
            let scale = 1.0 + o.as_array().iter().fold(0.0, |m: f64, c| m.max(c.abs()));
            cache.axis().distance(o) <= 1e-9 * scale
                && (cache.radius() - out.radius_f64()).abs() <= 1e-9 * scale
        },
        "{}",
        ONE_CYLINDER
    );
    Some(out)
}

/// The sentence [`world_cylinder_def`]'s postcondition panics with — one spelling, shared with
/// the commuting oracle's `KNOWN` list, which names a panic site by its sentence.
pub(crate) const ONE_CYLINDER: &str =
    "the world statement and the realized cache describe one cylinder";

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
/// side of the class's **stored** plane normal, which is the frame every cell label is written in
/// (`arrangement`'s `w_normal`).
///
/// ★ The `f64` dot is exact enough by construction: this is only ever asked of a class ⊥ to the
/// axis, so the dot is `±|n||m|` — a full magnitude from the sign boundary, not a near-zero
/// comparison.
///
/// ★★ **Ask this, do not re-derive it from the class's rational name.** `base_rat` is the
/// *canonical* name (first nonzero component positive), which points the other way from the stored
/// normal on half the classes; a rule spelled against it reads "above" backwards exactly there.
/// That mistake, made while adding the second consumer below, turned 36 tests red at once.
pub(crate) fn plus_t_is_above(wp: &WorkingPlane, def: &nacre_topo::CylinderDef) -> bool {
    let m = def.dir();
    let axis = Vector3::from_array([m[0].to_f64(), m[1].to_f64(), m[2].to_f64()]);
    wp.plane.normal().dot(axis) > 0.0
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
