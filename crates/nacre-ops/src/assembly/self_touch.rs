use super::*;
/// **What the tangency verdict needs, carried as one parameter** — the operation (the only thing
/// that decides `keep`) and the rows the gate stated. Empty for every boolean without a tangent
/// wall, which is almost all of them.
#[derive(Clone, Copy)]
pub(crate) struct Tangencies<'a> {
    pub(crate) kind: BoolKind,
    pub(crate) rows: &'a [crate::planes::Tangency],
}

impl Tangencies<'_> {
    /// For the roads that assemble without a gate run (the band-road fixtures). ★ The `kind` here
    /// is a placeholder, not a choice: with no rows the verdict returns before reading it, which is
    /// why naming an operation at this one site is not the casework the rule forbids.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Self {
            kind: BoolKind::Fuse,
            rows: &[],
        }
    }
}

/// **A tangency pinches the result when the material it leaves is two lumps that meet only along
/// the line — and both lumps end up in one solid.**
///
/// ★★★★★ **The operation enters here and only here, through [`crate::draft::keep`].** Near
/// the tangent line the material is three regions — the lens inside the cylinder, the **two**
/// wedges between the parabola and the plane, and the far half-space — and the adjacencies are
/// `L–W2` and `F–W2` only: `L` and `F` meet along the line itself and nowhere else, because at
/// `v = 0` the wedges are *inside* the cylinder. So the kept regions fall apart exactly when
///
/// ```text
/// (W2 ∧ ¬L ∧ ¬F)   the two wedges cannot reach each other
/// (L ∧ F ∧ ¬W2)    the lens and the far side meet only on the line
/// ```
///
/// Each region's `(in_A, in_B)` is local: the wall **face**'s material side and the lateral
/// **face**'s `orient_sign`, nothing else — which is why this is exact without touching the
/// arrangement. (It is exact only where no *other* plane holds the tangent line; such a plane is a
/// secant, so the picture there is six regions, and [`crate::planes::Tangency`] refuses to speak —
/// or, where that line is an edge both faces end on, hands the question to the structure.)
///
/// ★★ **Two lumps is not yet a defect** — and that was measured, not assumed. A boss tangent to a
/// wall *from outside* leaves `L ∧ F ∧ ¬W2`, and the honest answer is **two solids touching along
/// a line**, which the engine already produces (`Ok(2)`, `validate` clean, and it meshes). What
/// makes it a defect is the two lumps landing in **one** body, and that question is the grouping's,
/// not the geometry's: the tangency's wall class and cylinder class in one solid. The first case
/// above always answers yes (both wedges are bounded by the same wall face *and* the same lateral,
/// neither of which the tangency splits), and the second is where the grouping earns its keep — a
/// boss standing in a notch, tangent to the notch's wall and overlapping the block beside it, comes
/// back as one solid whose boundary touches itself, and nothing else in the kernel sees it.
/// **Does the material near a tangent line fall into more than one piece?** — the three regions
/// and their two adjacencies, asked through [`crate::draft::keep`] and nothing else.
///
/// `lens_in_wall_solid` says whether the cylinder's side of the wall plane is the wall **face**'s
/// material side; `cyl_orient` is `+1` for a boss (material inside) and `-1` for a bore. Those two
/// bits fix all six memberships, because near the line the only boundaries are those two surfaces.
pub(crate) fn lumps_fall_apart(
    kind: BoolKind,
    wall_solid: crate::planes::SolidSide,
    lens_in_wall_solid: bool,
    cyl_orient: i8,
) -> bool {
    use crate::planes::SolidSide;
    let region = |in_lens_side: bool, inside_cyl: bool| {
        let w = in_lens_side == lens_in_wall_solid; // in the wall face's solid
        let c = inside_cyl == (cyl_orient > 0); // in the lateral face's solid
        match wall_solid {
            SolidSide::A => crate::draft::keep(kind, w, c),
            SolidSide::B => crate::draft::keep(kind, c, w),
        }
    };
    let (l, w2, f) = (
        region(true, true),
        region(true, false),
        region(false, false),
    );
    (w2 && !l && !f) || (l && f && !w2)
}

pub(super) fn tangency_reject(
    tg: &Tangencies<'_>,
    faces: &[LocalFace],
    g: &Grouping,
) -> Result<(), BoolError> {
    use crate::planes::ClassIx;
    if tg.rows.is_empty() {
        return Ok(());
    }
    // The classes each *solid* carries — a solid is a material component plus its cavities.
    let bodies: Vec<Vec<crate::planes::ClassIx>> = g
        .positives
        .iter()
        .map(|m| {
            let comps = &g.comps_of[m];
            faces
                .iter()
                .enumerate()
                .filter(|(i, _)| comps.contains(&g.labels[*i]))
                .map(|(_, lf)| lf.surf)
                .collect()
        })
        .collect();
    for t in tg.rows {
        if t.undecided || (t.line_in_another_plane && !t.line_is_an_edge) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
        // ★ The line is an edge here (`Tangency::line_is_an_edge`): two lumps meeting on it meet
        // on an edge four faces use, and the shell guard names that — this row has nothing to
        // add, and the six regions around a line in a third plane are not its to judge.
        if t.line_is_an_edge {
            continue;
        }
        // Two lumps, and the contact is a *segment* — a corner grazing the line at a point is a
        // valid tangency (measured), and this must not convict it.
        if !lumps_fall_apart(tg.kind, t.wall_solid, t.lens_in_wall_solid, t.cyl_orient)
            || !t.straddles
        {
            continue;
        }
        let together = bodies
            .iter()
            .any(|cls| cls.contains(&ClassIx::Plane(t.wall)) && cls.contains(&ClassIx::Cyl(t.cyl)));
        if together {
            return Err(crate::reject_at(
                RejectReason::SelfTouchingResult,
                crate::RejectWhere::Point(t.witness),
            ));
        }
    }
    Ok(())
}

/// Push the reconstructed result and supersede the inputs.
///
/// ★ **The operands retire in exactly one place, and it is after the result is accepted.** Copied
/// into each arm that produces solids, every reject would have to know whether it stands before or
/// after its arm's copy; in one place, "a reject leaves the model alone" is true of the *shape* of
/// this function rather than of where the rejects happen to sit. [`reconstruct`] does the work and
/// never touches the live set.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_fuse_cut(
    model: &mut Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    seam: &[SeamVertex],
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    cut_rims: &crate::draft::CutRims,
    deferred: Option<BoolError>,
    tangencies: Tangencies<'_>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    let out = reconstruct(model, jd, seam, faces, cyls, cut_rims, deferred, tangencies)?;
    model.supersede_live(&[a, b]);
    Ok(out)
}

/// Per axis, the low and high ends of where a vertex's truth can be.
type Reach = ([f64; 3], [f64; 3]);

/// **Where a seam vertex's truth can be** — per axis, the interval its cache proves.
///
/// A `Bounded` cache holds the truth within `coord ± bound`; the radius is widened by a few ulps
/// of both operands (and the smallest normal), since `coord − bound` computed in `f64` may round
/// inward by half an ulp of the result — which the widening outgrows on either side of a binade.
/// Widening only over-keeps. Any other
/// variant proves nothing, and the answer is the whole line: the sieve then keeps every candidate
/// the vertex takes part in and leaves it to the exact test, which is the only safe reading of
/// "no bound". Measured: every seam point in the census and the perf folds is `Bounded`; the
/// suite's 72 others (a history past the cost cap, a carrier with no witness) are where this
/// arm runs.
fn reach(cache: &PointCache) -> Reach {
    match *cache {
        PointCache::Bounded { coord, bound } => {
            let c = coord.as_array();
            let w: [f64; 3] = core::array::from_fn(|k| {
                bound[k].upper_f64() * (1.0 + 4.0 * f64::EPSILON)
                    + c[k].abs() * (4.0 * f64::EPSILON)
                    + f64::MIN_POSITIVE
            });
            (
                core::array::from_fn(|k| c[k] - w[k]),
                core::array::from_fn(|k| c[k] + w[k]),
            )
        }
        PointCache::Ceiling { .. } | PointCache::Unrealized { .. } => {
            ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3])
        }
    }
}

/// **A solid whose surface touches itself is not a solid.**
///
/// The proposition: *an edge of this solid lies in the interior of one of this same solid's faces.*
/// An embedded boundary cannot do that — the surface would occupy the same points twice — so a
/// result that does is refused rather than returned. The unit is the **result solid**, outer shell
/// and cavities together: a wedge cut whose tip lands on the far wall touches its own cavity, and
/// asking the question per component would miss exactly that.
///
/// **Three sieves, cheapest first**, because the exact question is the expensive one:
///
/// 1. ★ **Which planes** — both endpoints must lie on the face's plane, and a vertex's plane triple
///    *is* the set of planes it lies on, so the candidates for an edge are `triple(u) ∩ triple(v)`,
///    usually two. Exact, no coordinates, a hash lookup per edge. Scanning every face for every
///    edge instead costs ~40% on the 80-fold star (measured).
/// 2. ★ **Which faces on that plane** — one plane can carry *eighty* faces in a folded star, so the
///    lookup above is not the end of it. A box per face rejects 99.99% of what survives step 1
///    (measured: 2,538,703 candidate pairs down to 203), and on the booleans that cost the most it
///    rejects all of them.
/// 3. The exact test, on what is left: [`combinatorics::segment_meets_face`] — **the open segment**,
///    not its endpoints. The endpoints alone miss the shape that motivated this: a cut taken all the
///    way *through* leaves a contact line whose ends sit on the face's ring while its middle crosses
///    the interior, and a point test reads both ends as "on the boundary, not inside" and passes a
///    body that cannot exist (measured: `Ok(n=1)`, volume exact, `validate` silent).
///
/// ★ **One judge, not two.** The endpoint test was here first and is *subsumed* — an endpoint
/// strictly inside puts that end in an inside stretch of the line, so the segment test answers it
/// too. Measured across the corpus: 1,696 candidate pairs, the two verdicts agreed on every one,
/// and no pair was undecidable. Keeping both would have left the same question answered in two
/// places.
///
/// ★ **The boxes hold where the truth can be, so they can only over-keep.** They are built from
/// realized `f64` coordinates, which are rounded; a box used to *reject* a candidate before an
/// exact test must therefore be conservative, or a real self-contact is dropped in silence. Each
/// vertex contributes its [`reach`] — the interval its cache proves the truth lies in — so the box
/// holds the true face, and a query asks whether the endpoint's reach meets it. The case this
/// protects is not hypothetical: the wedge that motivated this check touches at exactly `x = 1`,
/// which is a face of its own box.
///
/// The rings a face is trimmed by are built **lazily and once per face**: `Ring::edges` spends an
/// exact predicate per node, and building them for every result face is the same unconditional cost
/// that put 12x on the star when the void label did it.
pub(super) fn self_touch_reject(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    seam: &[SeamVertex],
    groups: &[Vec<usize>],
    by_comp_lf: &[Vec<&LocalFace>],
) -> Result<(), BoolError> {
    // Per node: the coordinate a reject names as its witness, and the reach the sieve reads.
    let pt: HashMap<NodeId, (Point3, Reach)> = seam
        .iter()
        .map(|sv| (sv.triple, (sv.cache.coord(), reach(&sv.cache))))
        .collect();
    for g in groups {
        // ★ **Planar faces only.** The proposition this sieve tests — "an edge of
        // this solid lies in the interior of one of its own faces" — is asked through plane
        // membership and in-plane boxes, and a lateral band has neither a plane nor a node to
        // offer. It abstains rather than pretending: the edge-use guard, the non-manifold vertex
        // check and `validate` all still stand behind it, so a curved self-touch surfaces there
        // instead of being missed silently. A band that pinches against a *plane* face is
        // therefore the case this does not see yet, and it is written down rather than assumed
        // away. ★ **An edge with a pierce endpoint joins that population**: it has no
        // three-plane name for the plane-membership question below, so it is skipped by the same
        // rule and for the same reason.
        let faces: Vec<&LocalFace> = g
            .iter()
            .flat_map(|&c| by_comp_lf[c].iter().copied())
            .filter(|lf| matches!(lf.surf, ClassIx::Plane(_)))
            .collect();
        // One pass: which faces sit on each plane, which two faces own each edge, and a box per
        // face holding its vertices' reach.
        let mut by_plane: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut owners: HashMap<[NodeId; 2], Vec<usize>> = HashMap::new();
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(faces.len());
        for (j, lf) in faces.iter().enumerate() {
            by_plane.entry(lf.surf.plane()).or_default().push(j);
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for r in rings_of(lf) {
                let k = r.nodes.len();
                for i in 0..k {
                    let (u, v) = (r.nodes[i], r.nodes[(i + 1) % k]);
                    owners
                        .entry(if u < v { [u, v] } else { [v, u] })
                        .or_default()
                        .push(j);
                    if let Some(&(_, (rl, rh))) = pt.get(&u) {
                        for k in 0..3 {
                            lo[k] = lo[k].min(rl[k]);
                            hi[k] = hi[k].max(rh[k]);
                        }
                    }
                }
            }
            boxes.push((lo, hi));
        }
        let mut rings_of_face: HashMap<usize, Vec<Vec<combinatorics::RingEdge>>> = HashMap::new();
        for (&[u, v], own) in &owners {
            // ★ **Abstain, exactly as the paragraph above already does for a lateral band.** An
            // edge with a pierce endpoint has no three-plane name to intersect, and this is a
            // *rejection guard*: declining here would turn "I cannot check this edge" into "this
            // model is invalid" — a reject on a possibly-sound solid. The nets named above (the
            // edge-use guard, the non-manifold vertex check, `validate`) still stand behind it,
            // and this population joins the one already written down there.
            let (Some(up), Some(vp)) = (
                combinatorics::three_plane_name(u),
                combinatorics::three_plane_name(v),
            ) else {
                continue;
            };
            for q in up.iter().copied().filter(|x| vp.contains(x)) {
                let Some(js) = by_plane.get(&q) else { continue };
                for &j in js.iter().filter(|j| !own.contains(j)) {
                    let (lo, hi) = boxes[j];
                    let in_box = [u, v].iter().all(|t| {
                        pt.get(t).is_some_and(|&(_, (rl, rh))| {
                            (0..3).all(|k| rh[k] >= lo[k] && rl[k] <= hi[k])
                        })
                    });
                    if !in_box {
                        continue;
                    }
                    // A face trimmed by a mixed ring joins the populations this sieve
                    // already passes over (a lateral face, a pierce-ended edge): the
                    // named walk cannot read it, and an `Err` here would turn "cannot
                    // check" into "model invalid". The net behind the skip stays what
                    // it is for those: the edge-use guard, the non-manifold vertex
                    // check, and `check_result_topology`'s plane-self edge guard.
                    if rings_of(faces[j]).any(Ring::is_mixed) {
                        continue;
                    }
                    if let std::collections::hash_map::Entry::Vacant(slot) = rings_of_face.entry(j)
                    {
                        slot.insert(
                            rings_of(faces[j])
                                .map(|r| r.edges(jd, cyls, Some(q)))
                                .collect::<Result<_, BoolError>>()?,
                        );
                    }
                    // ★ The line is `q ∩ w`, and `w` comes from **the edge's own two faces**, not
                    // from the endpoint names and not from the ring's recorded wall. Both of those
                    // fail, measured: an edge can lie on a pencil of three planes and the ring is
                    // free to call it by any of them — including `q` itself — and where four planes
                    // concur the endpoints' canonical triples need not share a second plane at all
                    // (31 pairs in the rotation sweep). An edge's two faces always give one, because
                    // that is what the edge *is*.
                    let Some(w2) = own.iter().map(|&o| faces[o].surf.plane()).find(|&x| x != q)
                    else {
                        continue; // an edge whose faces are both on `q` is not an edge
                    };
                    if combinatorics::segment_meets_face(jd, q, w2, up, vp, &rings_of_face[&j])? {
                        // The offending edge itself, as the witness: both endpoints are in `pt`
                        // (the `in_box` guard above already looked them up).
                        return Err(crate::reject_at(
                            RejectReason::SelfTouchingResult,
                            crate::RejectWhere::Segment([pt[&u].0, pt[&v].0]),
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// A face's rings, outer first.
pub(super) fn rings_of(lf: &LocalFace) -> impl Iterator<Item = &Ring> {
    lf.poly_rings()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_exact::Mag;

    /// **A vertex's reach holds its proven interval, and only a proven one is finite.**
    ///
    /// The sieve drops a candidate before the exact test whenever reaches miss a box, so a reach
    /// narrower than the truth's interval drops a real contact in silence. A `Bounded` reach must
    /// contain `coord ± bound` on every axis — an ordinary half-ulp, a zero bound (an exact
    /// coordinate), and a radius far below an `f64`'s — and a cache that proves nothing must
    /// keep every candidate.
    #[test]
    fn a_reach_holds_the_proven_interval_and_nothing_else_is_finite() {
        let coord = Point3::from_array([0.1, -3.0, 1.0e10]);
        let bound = [
            Mag::of(0.1).times(Mag::pow2(-53)),
            Mag::ZERO,
            Mag::pow2(-2000),
        ];
        let (lo, hi) = reach(&PointCache::Bounded { coord, bound });
        for k in 0..3 {
            let (c, r) = (coord.as_array()[k], bound[k].upper_f64());
            assert!(
                lo[k] < c && c < hi[k],
                "axis {k}: the coordinate itself is inside"
            );
            assert!(
                c - lo[k] >= r && hi[k] - c >= r,
                "axis {k}: the bound {r:e} is inside"
            );
            assert!(
                lo[k].is_finite() && hi[k].is_finite(),
                "axis {k}: a proven reach is finite"
            );
        }
        for cache in [
            PointCache::Ceiling { coord },
            PointCache::Unrealized { coord },
        ] {
            let (lo, hi) = reach(&cache);
            assert!(
                lo.iter().all(|x| *x == f64::NEG_INFINITY)
                    && hi.iter().all(|x| *x == f64::INFINITY),
                "{cache:?} proves no bound, so its reach is the whole line"
            );
        }
    }
}
