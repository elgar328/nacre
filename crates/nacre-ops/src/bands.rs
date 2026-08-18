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
//! while sitting outside the material. Asking the arrangement is what fixes it, and the
//! counterexample is a fixture of this pass (`a_cylinder_in_the_notch_keeps_no_band`, with its
//! end-to-end twin: cutting with a drill standing in the notch removes nothing).
//!
//! ★★ **The chamber is read off the arrangement, not measured.** The plane arrangement already
//! makes a **disk cell** for every circle a cylinder leaves on a class, and `label_cells` writes
//! four bits on it: which solid's material lies immediately above and below that plane *inside the
//! circle*. That is the band's chamber, stated by the engine that decided it.
//!
//! ★ This replaced a witness ray (C4a's `point_in_faces_rat`): a rational point on the axis, cast
//! through the counterpart's faces, counting crossings. It gave the same answers, but it answered
//! a **3D containment** question — the shape the *component* probe asks — when the band's question
//! is the same shape as "does this face survive": two chambers either side of a boundary. Reading
//! the label needs no coordinates, no ray direction, no abstention retry, and has no width
//! ceiling, and it is why a cylinder may now stand on either side of the boolean.

use crate::boolean::{Bound, LocalFace};
use crate::combinatorics;
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
/// `labels` is the arrangement's own answer per `(cylinder class, plane class)` — see the module
/// docs. Cylinders may sit on **both** operands: nothing here asks a solid to describe itself.
pub(crate) fn band_faces(
    kind: BoolKind,
    plane_faces: &[LocalFace],
    rows: &[CylRow],
    jd: &Judge<'_, WorkingPlane>,
    labels: &crate::arrangement::DiskLabels,
) -> Result<Vec<LocalFace>, BoolError> {
    let mut out = Vec::new();
    for (k, row) in rows.iter().enumerate() {
        for (lo, hi) in bands_of(k, row, plane_faces, jd)? {
            let (in_own_inside, in_other) = chamber(jd, row, k, lo, hi, labels)?;
            // ★ **The wall is a boundary face of its own solid, so that solid's membership flips
            // across it** — read inside from the label, and outside is its negation. The
            // counterpart does *not* flip: the gate keeps its boundary clear of the lateral, so
            // the slab theorem covers both sides. The face survives exactly when the two chambers
            // disagree under `keep`.
            let keep_side = |in_own: bool| match row.side {
                SolidSide::A => crate::arrangement::keep(kind, in_own, in_other),
                SolidSide::B => crate::arrangement::keep(kind, in_other, in_own),
            };
            let (keep_in, keep_out) = (keep_side(in_own_inside), keep_side(!in_own_inside));
            if keep_in == keep_out {
                continue;
            }
            // The stored cylinder normal points **away from the axis**. It is the result face's
            // outward normal when the kept material is inside the wall (a boss), and must be
            // flipped when the material is outside it (a hole). `keep_in` is that question
            // directly, whichever solid the wall came from.
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
        // ★ A class ⊥ to the axis is a **potential band boundary**, so failing to place it is a
        // refusal, not a skip: silently missing one would merge two bands whose membership
        // differs and hand back a closed, wrong solid. Classes that are not ⊥ cannot bound a
        // band at all (the gate proved the ∥ ones clear), so those are skipped for cause.
        let Some(coeffs) = combinatorics::class_coeffs_rat(jd, c) else {
            continue;
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        if !nacre_scalar::parallel_rat(&n, &row.def.dir()) {
            continue;
        }
        let t = param(jd, c, row)?;
        if t == row.span[0] || t == row.span[1] {
            push(c, &mut classes);
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

/// The axis parameter of a plane class, or the honest refusal earned by a class with no exact
/// description — and by one whose parameter is a **value** too wide for `Rat` (the road's own
/// name: the gate's questions are signs and were made total, this one is not).
fn param(jd: &Judge<'_, WorkingPlane>, c: usize, row: &CylRow) -> Result<Rat, BoolError> {
    param_opt(jd, c, row).ok_or_else(|| reject(RejectReason::WitnessNotRational))
}

fn param_opt(jd: &Judge<'_, WorkingPlane>, c: usize, row: &CylRow) -> Option<Rat> {
    let coeffs = combinatorics::class_coeffs_rat(jd, c)?;
    axis_param_of_plane(&coeffs, &row.def)
}

/// **Which chamber the band sits in, read off the arrangement.**
///
/// The disk cell on a boundary class carries `[inA_above, inA_below, inB_above, inB_below]` — the
/// material immediately above and below that plane *inside the circle*. The band leaves `lo`
/// going up in the axis parameter, so which of "above"/"below" faces it is decided by the class
/// normal against the axis (`sign(n · m)`); at `hi` the band is on the other side.
///
/// ★ **Both ends are read, and they must agree.** The uniform-slab theorem says the counterpart's
/// membership is constant across the open slab; two ends disagreeing means the theorem's premise
/// broke (a ∥ wall crossed the slab — the gate's clearance proof failing), and that is refused
/// rather than resolved by picking one. Free, because both labels are already in hand.
///
/// ★★ **Both solids' bits are data, not a control — and assuming one of them was a real bug.**
/// The first spelling asserted "the band lies inside its own cylinder, so its own solid's bit is
/// true". That holds for a *tool* (a drill solid fills its own cylinder) and is **false for a wall
/// inherited from an earlier bore**: inside that circle the plate has no material — the bore is a
/// hole. Both are ordinary inputs the moment a cylinder may sit on either operand, so the own-solid
/// membership is read from the label too, and the wall's other side follows from what the wall
/// *is*: a boundary face of that solid, so its own membership flips across it while the
/// counterpart's (by the slab theorem) does not.
fn chamber(
    jd: &Judge<'_, WorkingPlane>,
    row: &CylRow,
    k: usize,
    lo: usize,
    hi: usize,
    labels: &crate::arrangement::DiskLabels,
) -> Result<(bool, bool), BoolError> {
    let (cyl_bit, other_bit) = match row.side {
        SolidSide::A => (0usize, 2usize), // [A above, A below, B above, B below]
        SolidSide::B => (2usize, 0usize),
    };
    // `above` in the label is the class normal's side; the band leaves `lo` toward `hi`, i.e.
    // toward increasing axis parameter.
    let toward_hi = |c: usize| -> bool { jd.planes[c].plane.normal().dot(dir_f64(&row.def)) > 0.0 };
    let read = |c: usize, band_is_above: bool| -> Option<(bool, bool)> {
        let l = labels.get(&(k, c))?;
        let i = usize::from(!band_is_above);
        Some((l[cyl_bit + i], l[other_bit + i]))
    };
    let ends = [(lo, toward_hi(lo)), (hi, !toward_hi(hi))];
    let mut answer: Option<(bool, bool)> = None;
    for (c, band_is_above) in ends {
        let Some(pair) = read(c, band_is_above) else {
            continue; // this end carries no disk cell — the other end speaks
        };
        match answer {
            None => answer = Some(pair),
            // ★ The slab is uniform by the gate's clearance proof — **both** bits must agree at
            // the two ends. Ends disagreeing means that proof failed (something crossed the open
            // slab), and guessing which end to believe is exactly what this kernel does not do.
            // This is also where a flipped side selection surfaces: the two ends read opposite
            // sides, so getting the sign wrong makes them contradict.
            Some(prev) if prev != pair => {
                return Err(reject(RejectReason::CylinderGateUndecided));
            }
            Some(_) => {}
        }
    }
    // Every boundary class of a band carries a disk cell (the emitted circles put them there, and
    // the cylinder's own caps are always arranged) — so this is a wiring failure, not an input.
    answer.ok_or_else(|| reject(RejectReason::CylinderGateUndecided))
}

/// The cylinder's axis direction as `f64` — used only for the class-normal sign, where the class
/// is ⊥ to the axis and the dot is a full magnitude from zero.
fn dir_f64(def: &nacre_topo::CylinderDef) -> nacre_math::Vector3 {
    let m = def.dir();
    nacre_math::Vector3::from_array([m[0].to_f64(), m[1].to_f64(), m[2].to_f64()])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planes::{PlaneSetup, plane_index_setup};
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
        let setup = plane_index_setup(m, a, b).unwrap();
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
            ..
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
        let (plane_faces, disk_labels) = crate::arrangement::trace_result_faces_full_for_test(
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
        let out = band_faces(kind, &plane_faces, &rows, &jd, &disk_labels).expect("bands");
        // The classes' **z**, not their axis parameter: `t` is measured from the cylinder's own
        // origin along its raw `dir`, so a drill starting at z=−1 puts the box's cap at t=1. The
        // assertions read in world z, which is the vocabulary the fixtures are written in.
        let ts = (0..geom.len())
            .filter_map(|c| {
                let t = param_opt(&jd, c, &rows[0])?;
                // The class's world z, via the axis point at that parameter.
                let (o, m) = (rows[0].def.origin(), rows[0].def.dir());
                let z = o[2].checked_add(t.checked_mul(m[2])?)?;
                Some((c, z.to_f64()))
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

    /// ★★ **The first end-to-end cylinder boolean.** A `[0,2]³` box drilled through by a
    /// radius-0.5 bore: seven faces (four walls, two drilled caps, one hole wall), genus 1, and a
    /// volume of `8 − π·0.25·2`. The volume is the **winding lock** — props integrates by the
    /// divergence theorem, so a hole loop wound the wrong way returns `8 + πr²h` instead.
    #[test]
    fn a_through_hole_is_built_and_measures_what_it_should() {
        let (mut m, a, b) = box_and_drill(-1.0, 4.0);
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the drill cuts");
        assert_eq!(out.len(), 1, "one body");
        let s = out[0];
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let faces = crate::planes::solid_shell_handles(&m, s)
            .into_iter()
            .map(|sh| m.shells.get(sh).faces.len())
            .sum::<usize>();
        assert_eq!(faces, 7, "4 walls + 2 drilled caps + the bore's wall");
        // ★ **Genus 1, re-derived rather than quoted.** An earlier plan wrote `V10−E15+F7−L2`
        // from memory; the counts are measured here and what is asserted is the relation
        // (`χ = V − E + F − L = 2(S − G)`, so one shell with one through hole gives `χ = 0`).
        let reach = m.reachable();
        let (v_n, e_n) = (reach.vertices.len() as i64, reach.edges.len() as i64);
        let f_n = reach.faces.len() as i64;
        let l_n: i64 = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len() as i64)
            .sum();
        assert_eq!(
            v_n - e_n + f_n - l_n,
            0,
            "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        let v = nacre_props::mass_props(&m, s)
            .expect("a closed solid has mass props")
            .volume;
        let want = 8.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (v - want).abs() < 1e-6,
            "volume {v} vs {want} — a bore removes material, it does not add it"
        );
        // The same box without the bore, for the sign of the correction rather than its value.
        let mut m2 = Model::new();
        let plain = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        m2.rebuild_adjacency();
        let plain_v = nacre_props::mass_props(&m2, plain)
            .expect("a box has mass props")
            .volume;
        assert!(v < plain_v, "the bore removes material: {v} vs {plain_v}");
    }

    /// ★★ **Two bores in one plate — the case the band pass was rebuilt for.**
    ///
    /// The second cut's counterpart already carries a cylinder face (the first bore's wall), and
    /// the old witness road could not describe such a body, so this used to decline by name. The
    /// band pass now reads the arrangement's own disk labels, and a curved counterpart is no
    /// longer a question anyone has to answer.
    ///
    /// ★ It is also where the wall's **own** membership stopped being assumed: the first bore's
    /// wall bounds the *plate*, and inside that circle the plate has no material — the opposite of
    /// a drill, which fills its own cylinder. Assuming the drill's case (as the pass did while
    /// only drills were tested) puts the second bore's `keep` on the wrong chamber.
    #[test]
    fn a_two_hole_plate_drills_both() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 2.0, 1.0]),
        );
        let d1 = m.add_cylinder(
            Point3::from_array([1.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        let d2 = m.add_cylinder(
            Point3::from_array([3.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let one = crate::boolean(&mut m, BoolKind::Cut, a, d1).expect("first bore");
        let v1 = nacre_props::mass_props(&m, one[0]).expect("props").volume;
        assert!(
            (v1 - (8.0 - std::f64::consts::PI * 0.0625)).abs() < 1e-9,
            "one bore: {v1}"
        );
        let two = crate::boolean(&mut m, BoolKind::Cut, one[0], d2).expect("second bore");
        assert_eq!(two.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v2 = nacre_props::mass_props(&m, two[0]).expect("props").volume;
        let want = 8.0 - 2.0 * std::f64::consts::PI * 0.0625;
        assert!((v2 - want).abs() < 1e-9, "two bores: {v2} vs {want}");
        // Genus 2: two through holes in one shell.
        let reach = m.reachable();
        let l: i64 = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len() as i64)
            .sum();
        let chi =
            reach.vertices.len() as i64 - reach.edges.len() as i64 + reach.faces.len() as i64 - l;
        assert_eq!(chi, -2, "two handles: χ = 2(1 − 2)");
    }

    /// **Three bores.** Two was the case the label route was built for; three is the check that
    /// nothing in it counts to two — each cut's counterpart carries one more cylinder face than
    /// the last.
    #[test]
    fn three_bores_in_a_row() {
        let mut m = Model::new();
        let mut solid = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([6.0, 2.0, 1.0]),
        );
        for x in [1.0, 3.0, 5.0] {
            let d = m.add_cylinder(
                Point3::from_array([x, 1.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.25,
                3.0,
            );
            m.rebuild_adjacency();
            solid = crate::boolean(&mut m, BoolKind::Cut, solid, d).expect("a bore")[0];
        }
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, solid).expect("props").volume;
        let want = 12.0 - 3.0 * std::f64::consts::PI * 0.0625;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **A drilled plate cut by a body flush with its faces** — the population the seated rule
    /// was refusing while its own sentence promised something narrower.
    ///
    /// `Common` of a drilled plate with a box covering half of it: the two operands share the
    /// plate's top and bottom planes, which used to read as "a cylinder cap lies flush" and
    /// decline. No cap lies anywhere near them — the bore's rims end on *holed* faces — and the
    /// engine gets the answer exactly right, which is what the guard now lets it do.
    #[test]
    fn a_drilled_plate_can_be_cut_by_a_body_flush_with_its_faces() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 2.0, 1.0]),
        );
        let d = m.add_cylinder(
            Point3::from_array([1.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let drilled = crate::boolean(&mut m, BoolKind::Cut, plate, d).expect("bore")[0];
        let half = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Common, drilled, half).expect("flush common");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4.0 - std::f64::consts::PI * 0.0625;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **A cylinder whose own cap is coplanar with the other body's face declines as
    /// "seated" — even standing far away.** Two coplanar planes are one class whichever way they
    /// sit, so a *disk* of this cylinder on a class the other operand also has faces on is the
    /// deferred flush-seating population, distance notwithstanding.
    ///
    /// ★★ The rule used to fire on **any** shared ⊥ class, which refused ordinary work: a drilled
    /// plate cut by a body flush with its faces (see the test above). It now asks for this
    /// cylinder's own **cap face**, which is what its user-facing sentence always claimed.
    #[test]
    fn a_distant_cap_on_the_boxs_own_plane_declines_as_seated() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        // Far from the box in x/y, but its caps land exactly on z = 0 and z = 2.
        let b = m.add_cylinder(
            Point3::from_array([10.0, 10.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let before = m.live_solids.clone();
        let err = crate::boolean(&mut m, BoolKind::Cut, a, b).expect_err("coplanar caps");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::SeatedCylinderCap,
                    ..
                }
            ),
            "{err:?}"
        );
        // ★ **A reject consumes nothing.** The operands are still live and still the same two, in
        // the same order — the property the arena residue note is about, checked on the cylinder
        // population where the refusals are new.
        assert_eq!(m.live_solids, before, "a refused boolean retires nothing");
    }

    /// **A drill that misses.** `Cut` returns the box untouched and `Fuse` returns two bodies —
    /// the population where a cylinder is present but no circle is ever emitted, so the whole
    /// curved path has to stay out of the way.
    #[test]
    fn a_cylinder_that_misses_changes_nothing_and_fuses_apart() {
        let (mut m, a, b) = {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cylinder(
                Point3::from_array([10.0, 10.0, 0.5]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            (m, a, b)
        };
        let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
        assert_eq!(cut.len(), 1);
        let v = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
        assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");

        // ★ **The fuse this test's name promises — and the case that used to abort the kernel.**
        // Two disjoint bodies means `n = 2`, and classifying them asks a ray probe that cannot be
        // handed a lateral face. Until the guard landed, that asked a band face for its plane
        // class and panicked; the honest answer is a named refusal, and this is where it is
        // measured. (The name said "fuses apart" for a while before the body did — the fixture
        // owed the population its own name promised.)
        let mut m2 = Model::new();
        let a2 = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b2 = m2.add_cylinder(
            Point3::from_array([10.0, 10.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m2.rebuild_adjacency();
        let live_before = m2.live_solids.clone();
        let err = crate::boolean(&mut m2, BoolKind::Fuse, a2, b2)
            .expect_err("two bodies, one of them curved");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CurvedComponentDepth,
                    ..
                }
            ),
            "{err:?}"
        );
        assert_eq!(
            m2.live_solids, live_before,
            "a refused boolean retires nothing"
        );
    }

    /// **A blind hole, end to end.** The bore stops inside the box, so its own cap closes the
    /// bottom — a disk face that only the seated-circle path can emit.
    #[test]
    fn a_blind_hole_is_built_and_measures_what_it_should() {
        let (mut m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the blind bore cuts");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 8.0 - std::f64::consts::PI * 0.25 * 1.0; // the bore reaches z = 1
        assert!((v - want).abs() < 1e-9, "volume {v} vs {want}");
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

        // ★ …and end to end: `Cut` returns the prism untouched. That is the counterexample's real
        // assertion — a band pass answering by z-range would run a wall through the notch, and
        // this volume would move.
        let before = nacre_props::mass_props(&m, a).expect("props").volume;
        let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that changes nothing");
        assert_eq!(cut.len(), 1);
        let after = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
        assert!(
            (after - before).abs() < 1e-12,
            "the notch drill removes nothing: {before} → {after}"
        );
    }
}
