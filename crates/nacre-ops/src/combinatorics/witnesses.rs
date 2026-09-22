use super::*;
/// **The circle a planar face's boundary rides**, when its ring has one.
///
/// ★★★★★ **The question `ring_own_circle` asks, minus the part that made it too narrow.** That one
/// answers "is this ring *the whole* circle" — every edge an arc, chained the whole way round — and
/// so says `None` for a **cut** cap, which is a disk just the same with a chord across it. What a
/// face needs in order to name a point of its own interior is only *which circle bounds it*, and
/// that is: the ring's **arcs all ride one cylinder**, and the face's plane is **perpendicular to
/// that axis** (so the section really is a circle rather than an ellipse the gate would have
/// refused anyway). A ring with no arc at all has no circle — which is how a **wall panel**
/// (`[plane, ruling, plane, ruling]`) is turned away here rather than guessed at.
///
/// `ring_own_circle` is the special case with no chords; the two are kept apart because they answer
/// different questions — "is the ring a circle" versus "which circle bounds the face".
///
/// ☑ **Measured over the lib suite: 40 acceptances, and the only clause that ever refuses is the
/// first** — 20 rings with no arc at all (the wall panels). Two cylinders' arcs on one planar face
/// and a plane that is *not* perpendicular to the axis both refuse **0 times**, and they stay: the
/// first would take a circle that bounds only part of the ring, the second would call an ellipse a
/// circle and put the "centre" off the face — and a producer that stops holding either proposition
/// should be caught here rather than two layers down.
///
/// ⚠ **"The gate refuses an oblique cylinder cut" is false as a blanket sentence**:
/// the gate lets an oblique plane through when every lateral face
/// of the cylinder provably misses it. The conclusion stands on a narrower fact: this ring holds
/// an **arc** of that cylinder, so the plane does not miss it, and an oblique pair that meets is
/// what the gate refuses. [`class_carries_circle`] carries the argument.
fn face_circle<'a>(
    jd: &Judge<'_, WorkingPlane>,
    plane: usize,
    ring: &'a [RingEdge],
) -> Option<&'a nacre_topo::CylinderDef> {
    let mut arcs = ring.iter().filter_map(|e| match &e.carrier {
        Carrier::Arc(a) => Some(&**a),
        _ => None,
    });
    let first = arcs.next()?;
    if arcs.any(|a| a.cyl != first.cyl) {
        return None;
    }
    // The face's plane must be perpendicular to the axis: its normal is parallel to `m`.
    let coeffs = class_coeffs_rat(jd, plane)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let c = nacre_exact::cross3_rat(&n, &first.def.dir())?;
    let zero = nacre_exact::Rat::from_int(0);
    if c.iter().any(|v| *v != zero) {
        return None;
    }
    Some(&first.def)
}

/// **Points of the face's plane to offer the ring**, centre first — every one of them strictly
/// inside the circle, by derivation rather than by search.
///
/// The steps run along the class chart's own rational axes ([`Chart2dRat::axes`], the one spelling
/// for a rational basis of a plane) with `λ = r / (|e|² + 1)`. Then `λ²|e|² = r²·x/(x+1)²` for
/// `x = |e|²`, and `x/(x+1)²` is at most `1/4` (at `x = 1`), so `|λe| < r` — one rational
/// inequality, no magic constant and no halving loop. Which of them is inside the **face** is the
/// ring's question, not this one's.
///
/// ★ **And two points per chord**: a cap the wall cuts *off* the diameter can be a
/// segment thinner than any step from the centre reaches (an offset boss's Common, 0.2 deep
/// against `r/2 = 0.25`), and its corners are pierce names the vertex probe drops — the first
/// face with no witness at all. On the line through the centre along a chord's normal `n`,
/// `o − t·n`, the chord is at `t = q = (n·o + d)/|n|²` and the circle at `t² = T = r²/|n|²`
/// (both rational); the far side's point `t = 2qT/(q² + T)` lies beyond the chord (`|t| > |q|`
/// iff `T > q²`, the chord inside the circle) and inside the circle (`t² < T` iff
/// `(q² − T)² > 0`), the near side's `t = q/2` between the chord and the centre. A chord through
/// the centre (`q = 0`) is the axis steps' case and adds nothing. Which side is the face's is,
/// again, the ring's question. Not complete: a face whose second chord runs along that normal
/// line (a quarter of a segment) has the far point *on* that chord — named `RingHasNoWitness`
/// when it comes, in a multi-body result.
fn ring_interior_candidates(
    jd: &Judge<'_, WorkingPlane>,
    plane: usize,
    def: &nacre_topo::CylinderDef,
    centre: &[nacre_exact::Rat; 3],
    ring: &[RingEdge],
) -> Option<Vec<[nacre_exact::Rat; 3]>> {
    use nacre_exact::Rat;
    let coeffs = class_coeffs_rat(jd, plane)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let chart = Chart2dRat::of_normal(&n)?;
    let (e1, e2) = chart.axes();
    let mut out = vec![*centre];
    let r2 = def.r2();
    let sum = |a: &[Rat; 3], b: &[Rat; 3], neg: bool| -> Option<[Rat; 3]> {
        let mut v = [Rat::from_int(0); 3];
        for i in 0..3 {
            v[i] = if neg {
                a[i].checked_sub(b[i])?
            } else {
                a[i].checked_add(b[i])?
            };
        }
        Some(v)
    };
    // ★★★★★ **The chart's own lattice directions, not a search.** Two axes are not enough: a
    // chord can lie *along* one of them (then that step stays on the boundary) while the other's
    // ray crosses the circle exactly at the **seam**, where `circular_order_about_seam` has no
    // order to give and the arc step abstains. ☑ Measured: with `{±e1, ±e2}` alone, sixteen cap
    // faces answered `None` for every candidate, split exactly that way. The diagonals are off
    // both, and they cost one more derivation of the same inequality rather than a new rule.
    let diag: Vec<[Rat; 3]> = [false, true]
        .into_iter()
        .filter_map(|neg| sum(e1, e2, neg))
        .collect();
    let dirs: Vec<&[Rat; 3]> = [e1, e2].into_iter().chain(diag.iter()).collect();
    for e in dirs {
        let len2 = nacre_exact::dot3_rat(e, e)?;
        // A step strictly inside the circle along `e`. With a rational radius it is `r/(|e|²+1)`
        // — the spelling the corpus was measured with, kept verbatim so a stated radius walks the
        // same points it always did. Without one it is `min(r², 1)/(2(|e|²+1))`, inside because
        // `λ|e| ≤ 1/4 < 1 ≤ r` when `r ≥ 1` and `λ|e| ≤ r²/4 < r` when `r < 1`; the assertion
        // below is the judge either way.
        let lam = match nacre_exact::rat_sqrt_exact_big(r2) {
            Some(r) => Rat::new(
                r.numer().checked_mul(len2.denom())?,
                r.denom()
                    .checked_mul(len2.numer().checked_add(len2.denom())?)?,
            )?,
            None => {
                let one = Rat::from_int(1);
                let top = if *r2 < nacre_exact::BigRat::from(one) {
                    r2.narrow()?
                } else {
                    one
                };
                let bottom = Rat::from_int(2).checked_mul(len2.checked_add(one)?)?;
                top.checked_mul(Rat::new(bottom.denom(), bottom.numer())?)?
            }
        };
        for sign in [
            Rat::from_int(1),
            Rat::from_int(0).checked_sub(Rat::from_int(1))?,
        ] {
            let k = lam.checked_mul(sign)?;
            let mut p = *centre;
            for i in 0..3 {
                p[i] = p[i].checked_add(k.checked_mul(e[i])?)?;
            }
            debug_assert_eq!(
                nacre_exact::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), r2),
                nacre_exact::Orient::Negative,
                "a step of r/(|e|^2+1) along a chart axis stays strictly inside the circle"
            );
            out.push(p);
        }
    }
    // The chords' points, one wall class each.
    let mut walls: Vec<usize> = Vec::new();
    for e in ring {
        if let Carrier::Plane { wall, .. } = e.carrier
            && !walls.contains(&wall)
        {
            walls.push(wall);
        }
    }
    let zero = Rat::from_int(0);
    let recip = |v: Rat| Rat::new(v.denom(), v.numer());
    for wall in walls {
        let Some(w) = class_coeffs_rat(jd, wall) else {
            continue;
        };
        let wn = [w[0], w[1], w[2]];
        let Some(chord) = (|| {
            let nn = nacre_exact::dot3_rat(&wn, &wn)?;
            let q = nacre_exact::dot3_rat(&wn, centre)?
                .checked_add(w[3])?
                .checked_mul(recip(nn)?)?;
            if q == zero {
                return None; // through the centre: the axis steps' case
            }
            let t_cap = r2.narrow()?.checked_mul(recip(nn)?)?;
            let q2 = q.checked_mul(q)?;
            let t_far = Rat::from_int(2)
                .checked_mul(q)?
                .checked_mul(t_cap)?
                .checked_mul(recip(q2.checked_add(t_cap)?)?)?;
            let t_near = q.checked_mul(Rat::new(1, 2)?)?;
            Some([t_far, t_near])
        })() else {
            continue;
        };
        for t in chord {
            let mut p = *centre;
            for i in 0..3 {
                p[i] = p[i].checked_sub(t.checked_mul(wn[i])?)?;
            }
            debug_assert_eq!(
                nacre_exact::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), r2),
                nacre_exact::Orient::Negative,
                "a chord's near and far points stay strictly inside the circle"
            );
            out.push(p);
        }
    }
    Some(out)
}

/// **Coordinate probes of a component whose faces carry no vertex.**
///
/// The witness is a **cap disk's centre**: the face is planar with a circular outer bound, so the
/// centre is the axis point at that plane's own axis parameter — the rule `bands.rs` already reads
/// a class's position along an axis with ([`axis_param_of_plane`](crate::planes::axis_param_of_plane)).
/// It is strictly inside the circle for any positive radius, so nothing needs to test that.
///
/// ★ **A holed cap names a point between its rims** ([`holed_cap_witness`]). The centre
/// of an annulus is in its hole, not on the face, and a witness that is not on the boundary is a
/// confidently wrong depth rather than an abstention — but passing a holed cap over
/// leaves a **tube** (two annular caps, two bands, no vertex anywhere) with no witness at all and
/// the multi-body fuse refused `RingHasNoWitness` (measured, the bushing). The remedy is the one
/// the
/// cut cap already uses: candidates derived from the face's own radii, and the **face asked** which
/// is on it. A face none of them is on is still passed over.
///
/// Several directions per point, because one ray can graze and the remedy is another direction;
/// the order is not load-bearing.
pub(crate) fn coord_probes(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[CompFace],
) -> Vec<Probe> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let mut out = Vec::new();
    for f in faces {
        let CompSurf::Plane(q) = &f.surf else {
            continue;
        };
        // ★★★★★ **A cut cap is a disk too, and it names its own interior the same way.** The
        // witness has always been the circle's centre; what was missing is that a face whose
        // boundary a wall has cut still *has* a circle ([`face_circle`]), and its centre may then
        // lie **on** the chord rather than inside the face — exactly what happens when the axis
        // rides the wall, which makes the chord a diameter. So the centre is offered as one
        // candidate among several and the **ring itself** says which one is inside.
        let def = match &f.outer {
            BoundEdges::Circle(def) => &**def,
            BoundEdges::Ring(r) => match face_circle(jd, *q, r) {
                Some(def) => def,
                None => continue,
            },
            _ => continue,
        };
        let Some(centre) = class_coeffs_rat(jd, *q)
            .and_then(|coeffs| crate::planes::axis_param_of_plane(&coeffs, def))
            .and_then(|t| {
                let (o, m) = (def.origin(), def.dir());
                let mut p = [zero; 3];
                for k in 0..3 {
                    p[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
                }
                Some(p)
            })
        else {
            continue;
        };
        let p = if !f.inner.is_empty() {
            match holed_cap_witness(jd, cyls, *q, def, &centre, f) {
                Some(p) => p,
                None => continue,
            }
        } else {
            match &f.outer {
                // A whole circle: the centre is strictly inside for any positive radius, and asking
                // would only add a road where none is needed. Today's answer, unchanged.
                BoundEdges::Circle(_) => centre,
                // ★ **Derived, not searched.** Each step is `centre ± λ·e` along the class chart's own
                // rational axes ([`Chart2dRat::axes`] — the one spelling for "a rational basis of this
                // plane"), with `λ = r / (|e|² + 1)`. Then `λ²|e|² = r²·x/(x+1)²` for `x = |e|²`, and
                // `x/(x+1)² ≤ 1/4` at its maximum, so every candidate is strictly inside the circle —
                // a rational inequality, no magic constant and no halving loop.
                //
                // ★★★★★ **Which one is inside the *face* is asked, not derived.** Deriving it would
                // mean spelling "the material side of the chord" in some frame, and this road has no
                // oracle for that sign; the ring already answers the question exactly
                // ([`point_in_mixed_ring`]), and an abstention just moves to the next candidate. The
                // centre goes first, so a cap the wall cuts off-centre still answers with it.
                //
                // ☑ Measured: the centre and the axis steps answer 40 of 40 faces of the
                // through-axis corpus; the offset wall's thin segment answers by a chord point (8
                // faces), and the fall-through below is the named residual — a segment cut again
                // along the chord's own normal line.
                _ => {
                    let BoundEdges::Ring(r) = &f.outer else {
                        continue;
                    };
                    let Some(cand) = ring_interior_candidates(jd, *q, def, &centre, r) else {
                        continue;
                    };
                    let Some(coeffs) = class_coeffs_rat(jd, *q) else {
                        continue;
                    };
                    match cand
                        .into_iter()
                        .find(|c| point_in_mixed_ring(jd, cyls, &coeffs, c, r) == Some(true))
                    {
                        Some(c) => c,
                        None => continue,
                    }
                }
            }
        };
        out.extend(
            probe_dirs(&def.dir())
                .into_iter()
                .map(|dir| Probe::Coord { p, dir }),
        );
    }
    out
}

/// The directions a coordinate probe casts along: the axis both ways and one perpendicular both
/// ways — several, because one ray can graze and the remedy is another direction; the order is not
/// load-bearing. One spelling for [`coord_probes`], [`corner_probes`], and the component
/// road's edge supply, which hands in a **plane's normal** where the other two hand in a cylinder's
/// axis: the argument is only "a nonzero direction to build a frame from".
pub(crate) fn probe_dirs(m: &[nacre_exact::Rat; 3]) -> Vec<[nacre_exact::Rat; 3]> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let nonzero = |v: &[Rat; 3]| v.iter().any(|c| *c != zero);
    let neg = |v: &[Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            zero.checked_sub(v[0])?,
            zero.checked_sub(v[1])?,
            zero.checked_sub(v[2])?,
        ])
    };
    let mut dirs = vec![*m];
    dirs.extend(neg(m));
    for k in 0..3 {
        let mut e = [zero; 3];
        e[k] = Rat::from_int(1);
        match nacre_exact::cross3_rat(&e, m) {
            Some(w) if nonzero(&w) => {
                dirs.push(w);
                dirs.extend(neg(&w));
                break;
            }
            _ => {}
        }
    }
    dirs
}

/// **Coordinate probes at a component's rational pierce corners** — the corners a
/// prism with arcs has where its walls meet its cylinders (a slot's, a fillet's, a D-prism's:
/// no three-plane name anywhere, and [`coord_probes`]' cap witness is not always on the face).
/// A corner whose root is rational ([`pierce_coords_rat`]) is a point of the boundary as exact as
/// a named vertex, cast along its own cylinder's directions. Each corner once.
pub(crate) fn corner_probes(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    nodes: impl Iterator<Item = NodeId>,
) -> Vec<Probe> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for n in nodes {
        if !seen.insert(n) {
            continue;
        }
        let Some((_, cyl, _)) = pierce_name(n) else {
            continue;
        };
        let Some(p) = pierce_coords_rat(jd, cyls, n) else {
            continue;
        };
        let Some(c) = cyls.get(cyl) else {
            continue;
        };
        out.extend(
            probe_dirs(&c.def.dir())
                .into_iter()
                .map(|dir| Probe::Coord { p, dir }),
        );
    }
    out
}

/// **A rational point on a holed planar cap** — a face whose outer bound is a circle (or a cut
/// cap that still has one, [`face_circle`]) and whose holes are circles or rings: the annular cap
/// of a tube, a boss's cap around a pin's hole.
///
/// The candidates are **derived from the face's own radii** and the cylinder's own rational
/// frame, then the **face is asked**: `centre + ρ·û` for `û` each of the four rational unit
/// directions of the circle's chart (`ref_dir/|ref_dir|`, `(m × ref_dir)/(|m||ref_dir|)` and
/// their negatives — rational exactly when both norms are, which every prism on a world frame
/// has; a rotated frame has none and the face is passed over) and `ρ` the half-radius and, per
/// circle hole, the **mid-radius** `(R + r_hole)/2` — the ring between two concentric rims. Which
/// candidate is *on the face* is an exact question: inside the outer (the circle's radial side,
/// or the mixed ring's parity) and outside every hole.
///
/// `None`: no rational frame, or no candidate on the face — the caller passes the face over, as
/// it did every holed cap before this was written.
fn holed_cap_witness(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    plane: usize,
    def: &nacre_topo::CylinderDef,
    centre: &[nacre_exact::Rat; 3],
    f: &CompFace,
) -> Option<[nacre_exact::Rat; 3]> {
    use nacre_exact::{Orient, Rat, inv_sqrt_exact, quad::cylinder_radial_side};
    let dot = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<Rat> {
        let mut acc = Rat::from_int(0);
        for k in 0..3 {
            acc = acc.checked_add(a[k].checked_mul(b[k])?)?;
        }
        Some(acc)
    };
    let scale = |v: &[Rat; 3], s: Rat| -> Option<[Rat; 3]> {
        Some([
            v[0].checked_mul(s)?,
            v[1].checked_mul(s)?,
            v[2].checked_mul(s)?,
        ])
    };
    let (m, e) = (def.dir(), def.ref_dir());
    // `1/|ref_dir|` and `1/(|m||ref_dir|)`, exact when the norms are (`inv_sqrt_exact`).
    let inv_e = inv_sqrt_exact(dot(&e, &e)?)?;
    let inv_me = inv_e.checked_mul(inv_sqrt_exact(dot(&m, &m)?)?)?;
    let u1 = scale(&e, inv_e)?;
    let u2 = scale(&nacre_exact::cross3_rat(&m, &e)?, inv_me)?;
    let neg = |v: &[Rat; 3]| scale(v, Rat::from_int(-1));
    let dirs = [u1, neg(&u1)?, u2, neg(&u2)?];
    // Probe circles midway between the outer rim and each hole's, at rational radii — which
    // needs the radii themselves. A face whose squared radii have no rational root is not probed
    // this way; `None` is the answer this function already gives when it cannot form a witness.
    let big_r = def.radius_exact()?;
    let half = Rat::new(1, 2)?;
    let mut radii = vec![big_r.checked_mul(half)?];
    for hole in &f.inner {
        if let BoundEdges::Circle(h) = hole {
            radii.push(big_r.checked_add(h.radius_exact()?)?.checked_mul(half)?);
        }
    }
    let coeffs = class_coeffs_rat(jd, plane);
    // On the face: inside the outer, outside every hole — each an exact question of the bound.
    let on_face = |p: &[Rat; 3]| -> Option<bool> {
        let inside_outer = match &f.outer {
            BoundEdges::Circle(_) => {
                cylinder_radial_side(p, &def.origin(), &def.dir(), def.r2()) == Orient::Negative
            }
            BoundEdges::Ring(r) => point_in_mixed_ring(jd, cyls, coeffs.as_ref()?, p, r)?,
            BoundEdges::Lateral(_) => return None,
        };
        if !inside_outer {
            return Some(false);
        }
        for hole in &f.inner {
            let inside_hole = match hole {
                BoundEdges::Circle(h) => {
                    cylinder_radial_side(p, &h.origin(), &h.dir(), h.r2()) != Orient::Positive
                }
                BoundEdges::Ring(r) => point_in_mixed_ring(jd, cyls, coeffs.as_ref()?, p, r)?,
                BoundEdges::Lateral(_) => return None,
            };
            if inside_hole {
                return Some(false);
            }
        }
        Some(true)
    };
    for rho in &radii {
        for u in &dirs {
            let step = scale(u, *rho)?;
            let mut p = *centre;
            for k in 0..3 {
                p[k] = p[k].checked_add(step[k])?;
            }
            if on_face(&p)? {
                return Some(p);
            }
        }
    }
    None
}

/// **Two rational planes whose meet is the line through `p` along `dir`.**
///
/// ★ The normals are **basis crosses** (`ê_k × dir`), whose components are a shuffle of `dir`'s —
/// no products, nothing to overflow. That is the same rule `SketchPlane::normal_def` states, and
/// for the same reason: the "obvious" second direction `n × u` squares the inputs' widths.
pub(super) fn planes_through_line(
    p: &[nacre_exact::Rat; 3],
    dir: &[nacre_exact::Rat; 3],
) -> Option<[[nacre_exact::Rat; 4]; 2]> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let basis = |k: usize| -> [Rat; 3] {
        let mut e = [zero; 3];
        e[k] = Rat::from_int(1);
        e
    };
    let mut ns: Vec<[Rat; 3]> = Vec::new();
    for k in 0..3 {
        if let Some(n) = nacre_exact::cross3_rat(&basis(k), dir) {
            if n.iter().any(|c| *c != zero)
                && (ns.is_empty()
                    || nacre_exact::cross3_rat(&ns[0], &n)
                        .is_some_and(|c| c.iter().any(|v| *v != zero)))
            {
                ns.push(n);
            }
        }
        if ns.len() == 2 {
            break;
        }
    }
    let [n0, n1] = <[[Rat; 3]; 2]>::try_from(ns).ok()?;
    let plane = |n: [Rat; 3]| -> Option<[Rat; 4]> {
        Some([
            n[0],
            n[1],
            n[2],
            zero.checked_sub(nacre_exact::dot3_rat(&n, p)?)?,
        ])
    };
    Some([plane(n0)?, plane(n1)?])
}
