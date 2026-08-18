//! **Which parts of a cylinder's lateral surface survive a boolean** (M6-2a C4b).
//!
//! The plane arrangement decides one plane class at a time; a cylinder's wall is not a plane, so
//! it is decided here instead — and it is decided *coarsely*, because in this population it can
//! be. The M6-2a gate admits only ⊥ cuts and walls that stand clear of the lateral surface by
//! more than `r`, and those two facts together give the **uniform-slab theorem**:
//!
//! > Between two consecutive ⊥ cuts, the other operand's boundary does not meet the open
//! > cylinder slab at all — the ⊥ faces are the cuts themselves and every ∥ wall misses the
//! > lateral. So membership in the other operand is **uniform** over that slab, and one witness
//! > decides the whole band.
//!
//! ★ The theorem replaced a wrong one. "The band is in the other solid iff its z-range is" reads
//! plausibly and is **false**: an L-notch's inner corner is clear of every wall by more than `r`
//! while sitting outside the material. The witness query is what fixes it, and the counterexample
//! is a fixture of the road that answers it (`combinatorics::point_in_faces_rat`).
//!
//! The witness is a **rational point on the axis** — the band's midpoint — which is why the road
//! had to take coordinates rather than a plane triple (C4a).

use crate::boolean::{Bound, LocalFace};
use crate::combinatorics::{self, ComponentTriples};
use crate::planes::SolidSide;
use crate::planes::{ClassIx, FaceRow, WorkingPlane, axis_param_of_plane};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, reject};
use nacre_scalar::Rat;

/// One cylinder class, as the band pass reads it.
#[cfg_attr(not(test), allow(dead_code))] // the production caller arrives with C4b-3
pub(crate) struct CylRow {
    pub(crate) def: nacre_topo::CylinderDef,
    /// The lateral face's own extent in the axis parameter (its two rims).
    pub(crate) span: [Rat; 2],
    /// Which operand the cylinder belongs to.
    pub(crate) side: SolidSide,
}

/// One row per cylinder class, in `ClassIx::Cyl` numbering order.
///
/// `None` for a class whose lateral face could not state its rim span — the tracer declines that
/// face for the same reason (`DeclineKind::CylSpan`), so a missing span here is a class the
/// arrangement has already refused.
#[cfg_attr(not(test), allow(dead_code))] // the production caller arrives with C4b-3
pub(crate) fn cyl_rows(
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
) -> Result<Vec<CylRow>, BoolError> {
    let mut out: Vec<Option<CylRow>> = Vec::new();
    for (i, row) in faces.iter().enumerate() {
        let (ClassIx::Cyl(k), FaceRow::Cylinder(cf)) = (plane_ix[i], row) else {
            continue;
        };
        if out.len() <= k {
            out.resize_with(k + 1, || None);
        }
        if out[k].is_some() {
            continue; // one row per class; the first lateral face of it speaks
        }
        let Some(span) = cf.span else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        out[k] = Some(CylRow {
            def: cf.def.clone(),
            span,
            side: if i < n_a { SolidSide::A } else { SolidSide::B },
        });
    }
    out.into_iter()
        .map(|r| r.ok_or_else(|| reject(RejectReason::CylinderGateUndecided)))
        .collect()
}

/// **The surviving lateral bands, as result faces.**
///
/// `plane_faces` is the plane arrangement's output: its circle bounds are where the cuts are, and
/// taking the boundaries from *there* (rather than from every ⊥ class) is what makes each rim
/// shared by exactly the cap face and the band that meet on it.
///
/// Emitted in `(cylinder class, lower t)` order — `reconstruct` mints handles in face order, so
/// this ordering is the replay contract.
#[cfg_attr(not(test), allow(dead_code))] // the production caller arrives with C4b-3
/// `counterpart` is the **other operand's** planar faces — the only road the witness travels.
/// Taking one road rather than both is the API saying what the decline below enforces: every
/// cylinder in this population belongs to one side, and the side that is asked about is planar.
pub(crate) fn band_faces(
    kind: BoolKind,
    plane_faces: &[LocalFace],
    rows: &[CylRow],
    jd: &Judge<'_, WorkingPlane>,
    counterpart: &ComponentTriples,
) -> Result<Vec<LocalFace>, BoolError> {
    // The witness road reads **planar** faces. A counterpart that is itself a cylinder has no
    // road yet, and the gate's "parallel axes, clear of each other" proof is the *gate's*, not an
    // answer this pass derived — borrowing it here would be a shortcut resting on someone else's
    // assumption. Named decline instead.
    if rows.iter().any(|r| r.side == SolidSide::A) && rows.iter().any(|r| r.side == SolidSide::B) {
        return Err(reject(RejectReason::BandWitnessNotPlanar));
    }
    let mut out = Vec::new();
    for (k, row) in rows.iter().enumerate() {
        for (lo, hi) in bands_of(k, row, plane_faces, jd)? {
            let (lo_t, hi_t) = (param(jd, lo, row)?, param(jd, hi, row)?);
            let Some(mid) = midpoint(lo_t, hi_t) else {
                return Err(reject(RejectReason::CylinderGateUndecided));
            };
            let Some(witness) = axis_point(&row.def, mid) else {
                return Err(reject(RejectReason::CylinderGateUndecided));
            };
            let Some(in_other) = decide(jd, &witness, &row.def, counterpart)? else {
                return Err(reject(RejectReason::NoClearRay));
            };
            // Inside the wall the cylinder's solid is present, outside it is not; the other
            // operand is uniform across the slab (the theorem above). The face survives exactly
            // when those two chambers disagree under `keep`.
            let keep_side = |in_cyl: bool| match row.side {
                SolidSide::A => crate::arrangement::keep(kind, in_cyl, in_other),
                SolidSide::B => crate::arrangement::keep(kind, in_other, in_cyl),
            };
            let (keep_in, keep_out) = (keep_side(true), keep_side(false));
            if keep_in == keep_out {
                continue;
            }
            // The stored cylinder normal points **away from the axis**. It is the result face's
            // outward normal when the kept material is inside the wall (a boss), and must be
            // flipped when the material is outside it (a hole).
            out.push(LocalFace {
                surf: ClassIx::Cyl(k),
                outer: Bound::Band { lo, hi },
                inner: Vec::new(),
                flip: !keep_in,
            });
        }
    }
    Ok(out)
}

/// The plane classes bounding this cylinder's bands, as consecutive `(lo, hi)` pairs ordered by
/// axis parameter.
///
/// The boundaries are the classes that **emitted a circle** for this cylinder, plus the classes
/// the lateral's own rims sit on (found by matching the span's parameters — two ⊥ planes with the
/// same axis parameter *are* one plane, so the match is exact, not a tolerance).
fn bands_of(
    k: usize,
    row: &CylRow,
    plane_faces: &[LocalFace],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<(usize, usize)>, BoolError> {
    let mut classes: Vec<usize> = Vec::new();
    let push = |c: usize, classes: &mut Vec<usize>| {
        if !classes.contains(&c) {
            classes.push(c);
        }
    };
    for lf in plane_faces {
        let ClassIx::Plane(c) = lf.surf else { continue };
        let carries = std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, Bound::Circle { cyl } if *cyl == k));
        if carries {
            push(c, &mut classes);
        }
    }
    for c in 0..jd.planes.len() {
        if let Some(t) = param_opt(jd, c, row) {
            if t == row.span[0] || t == row.span[1] {
                push(c, &mut classes);
            }
        }
    }
    let mut keyed: Vec<(Rat, usize)> = classes
        .into_iter()
        .map(|c| Ok((param(jd, c, row)?, c)))
        .collect::<Result<_, BoolError>>()?;
    keyed.retain(|(t, _)| *t >= row.span[0] && *t <= row.span[1]);
    keyed.sort_by_key(|(t, _)| *t);
    Ok(keyed.windows(2).map(|w| (w[0].1, w[1].1)).collect())
}

/// The axis parameter of a plane class, or the honest refusal a class without an exact
/// description earns.
fn param(jd: &Judge<'_, WorkingPlane>, c: usize, row: &CylRow) -> Result<Rat, BoolError> {
    param_opt(jd, c, row).ok_or_else(|| reject(RejectReason::CylinderGateUndecided))
}

fn param_opt(jd: &Judge<'_, WorkingPlane>, c: usize, row: &CylRow) -> Option<Rat> {
    let coeffs = combinatorics::class_coeffs_rat(jd, c)?;
    axis_param_of_plane(&coeffs, &row.def)
}

/// The exact midpoint of two axis parameters.
fn midpoint(a: Rat, b: Rat) -> Option<Rat> {
    a.checked_add(b)?.checked_mul(Rat::new(1, 2)?)
}

/// `origin + t·dir`, exactly.
fn axis_point(def: &nacre_topo::CylinderDef, t: Rat) -> Option<[Rat; 3]> {
    let (o, m) = (def.origin(), def.dir());
    let mut p = o;
    for k in 0..3 {
        p[k] = p[k].checked_add(t.checked_mul(m[k])?)?;
    }
    Some(p)
}

/// Whether the witness lies inside the counterpart solid, asking along several directions until
/// one decides.
///
/// Several, because a single ray can graze a face boundary and abstain — the road says so in its
/// type, and the remedy is another direction. ★ **The order is not load-bearing** and this
/// comment used to claim it was: "the axis is the worst direction, a drilled cap's hole is
/// centred on it". Measured false for this pass — the road reads the **operand's** faces, not the
/// result's, and the operand's cap is not drilled. (Where a hole *does* exist in an operand — a
/// plate from an earlier boolean — the road's input cannot describe it at all, which is that
/// builder's contract, not this ordering's job; see `point_in_faces_rat`.)
fn decide(
    jd: &Judge<'_, WorkingPlane>,
    witness: &[Rat; 3],
    def: &nacre_topo::CylinderDef,
    other: &ComponentTriples,
) -> Result<Option<bool>, BoolError> {
    let (m, u) = (def.dir(), def.ref_dir());
    let neg = |v: &[Rat; 3]| -> Option<[Rat; 3]> {
        let zero = Rat::from_int(0);
        Some([
            zero.checked_sub(v[0])?,
            zero.checked_sub(v[1])?,
            zero.checked_sub(v[2])?,
        ])
    };
    let cross = |x: &[Rat; 3], y: &[Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            x[1].checked_mul(y[2])?
                .checked_sub(x[2].checked_mul(y[1])?)?,
            x[2].checked_mul(y[0])?
                .checked_sub(x[0].checked_mul(y[2])?)?,
            x[0].checked_mul(y[1])?
                .checked_sub(x[1].checked_mul(y[0])?)?,
        ])
    };
    let mut dirs: Vec<[Rat; 3]> = vec![u];
    if let Some(d) = neg(&u) {
        dirs.push(d);
    }
    if let Some(w) = cross(&m, &u) {
        dirs.push(w);
        if let Some(d) = neg(&w) {
            dirs.push(d);
        }
    }
    dirs.push(m);
    for d in dirs {
        if let Some(v) = combinatorics::point_in_faces_rat(jd, witness, &d, other)? {
            return Ok(Some(v));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planes::{PlaneSetup, plane_index_setup_past_stopper};
    use nacre_math::{Point3, Vector3};
    use nacre_store::Handle;
    use nacre_topo::{Model, Solid};

    /// Run the plane arrangement past the C2 stopper and hand the band pass what it needs.
    /// Returns `(band faces, the plane classes' axis parameters keyed by class)`.
    fn bands(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        kind: BoolKind,
    ) -> (Vec<LocalFace>, Vec<(usize, f64)>) {
        let setup = plane_index_setup_past_stopper(m, a, b).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            m,
            [(a, inc_a), (b, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
        );
        let plane_faces = crate::arrangement::trace_result_faces_for_test(
            m,
            kind,
            a,
            b,
            &jd,
            faces_tab,
            plane_ix,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the drill population traces");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("cylinder rows");
        // Only the planar operand gets a road — the cylinder-bearing one is never asked about
        // (and `component_triples` would rightly refuse to describe its lateral face).
        let counterpart = crate::tests::component_triples(m, a, &setup, &jd);
        let out = band_faces(kind, &plane_faces, &rows, &jd, &counterpart).expect("bands");
        // The classes' **z**, not their axis parameter: `t` is measured from the cylinder's own
        // origin along its raw `dir`, so a drill starting at z=−1 puts the box's cap at t=1. The
        // assertions read in world z, which is the vocabulary the fixtures are written in.
        let ts = (0..geom.len())
            .filter_map(|c| {
                let t = param_opt(&jd, c, &rows[0])?;
                let p = axis_point(&rows[0].def, t)?;
                Some((c, p[2].to_f64()))
            })
            .collect();
        (out, ts)
    }

    /// A band's two ends as world `z` — the assertion vocabulary.
    fn ends(lf: &LocalFace, ts: &[(usize, f64)]) -> (f64, f64) {
        let Bound::Band { lo, hi } = lf.outer else {
            panic!("a band face bounds a band");
        };
        let at = |c: usize| ts.iter().find(|(k, _)| *k == c).expect("a ⊥ class").1;
        (at(lo), at(hi))
    }

    fn box_and_drill(z0: f64, h: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cylinder(
            Point3::from_array([1.0, 1.0, z0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
        m.rebuild_adjacency();
        (m, a, b)
    }

    /// **A through hole.** The drill runs from z=−1 to z=3 through a `[0,2]³` box; `Cut` keeps
    /// exactly the wall between the box's two caps, and the material there is **outside** the
    /// cylinder, so the face is flipped against the stored (outward) normal.
    #[test]
    fn a_through_hole_keeps_the_band_between_the_caps_and_flips_it() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Cut);
        assert_eq!(out.len(), 1, "one surviving band: {out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
        assert!(out[0].flip, "a hole's wall faces its own axis");
        assert!(matches!(out[0].surf, ClassIx::Cyl(0)));
    }

    /// The same drill, fused instead: the material is **inside** the cylinder exactly where the
    /// box is not, so the two protruding stretches survive and the buried one does not — and
    /// neither survivor is flipped.
    #[test]
    fn a_fused_drill_keeps_the_two_protruding_bands_unflipped() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Fuse);
        let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
        spans.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        assert_eq!(spans, vec![(-1.0, 0.0), (2.0, 3.0)], "{out:?}");
        assert!(out.iter().all(|lf| !lf.flip), "a boss's wall faces outward");
    }

    /// `Common` keeps the buried stretch and nothing else — the mirror of the fuse case, and the
    /// one that would pass if the keep rule ignored the operand order.
    #[test]
    fn a_common_keeps_only_the_buried_band() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Common);
        assert_eq!(out.len(), 1, "{out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
        assert!(!out[0].flip, "the kept material is inside the wall");
    }

    /// **A blind hole.** The drill stops at z=1 inside the box: the wall survives from the box's
    /// cap down to the drill's own cap, which is where the bottom disk closes it.
    #[test]
    fn a_blind_hole_stops_at_the_drills_own_cap() {
        let (m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
        let (out, ts) = bands(&m, a, b, BoolKind::Cut);
        assert_eq!(out.len(), 1, "{out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 1.0));
        assert!(out[0].flip);
    }

    /// ★ **The uniform-slab theorem's counterexample, as a band case.** A cylinder standing in an
    /// L-prism's notch is clear of every wall by more than `r` — the refuted "z-range" rule would
    /// keep its wall — but it is outside the material, so `Cut` changes nothing and no band
    /// survives.
    #[test]
    fn a_cylinder_in_the_notch_keeps_no_band() {
        let profile = crate::ops::Profile2d::polygon(vec![
            nacre_math::Point2::from_array([0.0, 0.0]),
            nacre_math::Point2::from_array([2.0, 0.0]),
            nacre_math::Point2::from_array([2.0, 1.0]),
            nacre_math::Point2::from_array([1.0, 1.0]),
            nacre_math::Point2::from_array([1.0, 2.0]),
            nacre_math::Point2::from_array([0.0, 2.0]),
        ])
        .unwrap();
        let mut m = crate::ops::replay(&[crate::tests::extrude_log_op(profile, 1.0)]).unwrap();
        let a = m.live_solids[0];
        // In the notch (x,y ∈ [1,2]²), a slim drill standing clear of both notch walls.
        let b = m.add_cylinder(
            Point3::from_array([1.5, 1.5, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let (out, _) = bands(&m, a, b, BoolKind::Cut);
        assert!(
            out.is_empty(),
            "the notch is void — nothing to cut: {out:?}"
        );
    }
}
