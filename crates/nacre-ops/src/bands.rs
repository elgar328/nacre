//! **Which parts of a cylinder's lateral surface survive a boolean** (M6-2a C4b).
//!
//! The plane arrangement decides one plane class at a time; a cylinder's wall is not a plane, so
//! it is decided here instead — and it is decided *coarsely*, because in this population it can
//! be. The M6-2a gate admits only ⊥ cuts and ∥ walls whose **faces** stand clear of the lateral
//! surface, and those two facts together give the **uniform-slab theorem**:
//!
//! > Between two consecutive ⊥ cuts, the other operand's boundary does not meet the open
//! > cylinder slab at all — the ⊥ faces are the cuts themselves and every ∥ wall **face** misses
//! > the lateral. So membership in the other operand is **uniform** over that slab, and one
//! > witness decides the whole band.
//!
//! ★ **Faces, because that is what a boundary is made of.** The gate long tested each ∥ wall's
//! infinite *plane*, which is a cheaper sufficient condition and reads almost the same — but it
//! refused bodies standing well clear of the cylinder whose plane, extended, crossed it. The
//! theorem never asked for that; it asks about the boundary, and the gate now does too.
//!
//! ★★ **A wall face is asked about a rectangle, and the "consecutive cuts" above are per face.**
//! What a ∥ plane meets of the solid cylinder is a rectangle of that plane: the strip across, a
//! lateral face's axis-parameter span along. So the gate reads two separating axes, and a wall
//! clear on either one misses the rectangle (`planes.rs`' `face_clears_footprint`). With several
//! lateral faces there are several rectangles, and that is exactly the granularity the theorem
//! needs: `bands_of` clips every cut to its own row's span, so no band is ever built in the gap
//! between two of them, and a wall sitting in such a gap breaks no premise — there is none there
//! to break. Reading the spans as one `min..max` would invent both the band and the refusal.
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

/// **One lateral face**, as the band pass reads it.
///
/// ★★ **The row is a face, not a class.** A cylinder *class* is one lateral `Surface` handle
/// (`planes.rs`' `cyl_ix.entry(cf.surf)`), and one surface can carry **several faces of one
/// solid** — bore a plate through and then cut away the middle of the bore, and its wall becomes
/// two disjoint bands. "Which parts survive" is a question about a face; the class only says which
/// surface that face lies on.
///
/// This used to be one row per class, with the first face speaking for the rest. The others' spans
/// were then clipped away in [`bands_of`] and their bands never emitted, which left the result
/// shell open — the boolean came back `OpenResultShell`, naming a symptom of our own omission.
///
/// ★ **Merging the spans into one `min..max` is not the fix either**: it would invent a band across
/// the gap where the solid has no face at all.
pub(crate) struct CylRow {
    /// The [`ClassIx::Cyl`] payload — which lateral surface this face lies on. Several rows may
    /// share it.
    pub(crate) class: usize,
    pub(crate) def: nacre_topo::CylinderDef,
    /// **This face's** own extent in the axis parameter (its two rims).
    pub(crate) span: [Rat; 2],
    /// Which operand this face belongs to.
    pub(crate) side: SolidSide,
}

/// One row per lateral **face**, ordered by `(class, lower t)`.
///
/// ★ That order is [`band_faces`]' replay contract stated at the row level. A class's faces have
/// disjoint spans, so sorting by `(class, span[0])` keeps each class's bands in ascending `t` —
/// the contract generalizes rather than bends. With one face per class it *is* the old class-index
/// order, which is why existing results do not move.
///
/// ★★ **That disjointness holds because every lateral face is a full 2π band today.** M6-2b's
/// θ-partial faces will put two faces at the *same* `t`, and the tie then falls to the stable
/// sort's face order — still deterministic, but the sentence above stops being the reason. Naming
/// the premise here so the day it expires is a thing a reader can check, not a surprise.
///
/// A face whose rim span cannot be stated declines by name (`DeclineKind::CylSpan`) — the tracer
/// declines it for the same reason.
pub(crate) fn cyl_rows(
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
) -> Result<Vec<CylRow>, BoolError> {
    let mut out: Vec<CylRow> = Vec::new();
    let mut n_class = 0usize;
    for (i, row) in faces.iter().enumerate() {
        let (ClassIx::Cyl(k), FaceRow::Cylinder(cf)) = (plane_ix[i], row) else {
            continue;
        };
        n_class = n_class.max(k + 1);
        let Some(span) = cf.span else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        out.push(CylRow {
            class: k,
            def: cf.def.clone(),
            span,
            side: if i < n_a { SolidSide::A } else { SolidSide::B },
        });
    }
    out.sort_by(|a, b| (a.class, a.span[0]).cmp(&(b.class, b.span[0])));
    // ★ The guard the old per-class table carried, kept as its own sentence: a class exists
    // because a face made it, so a class with no row is a wiring failure rather than an input.
    // It has never fired; an unfired *guard* stays (an unfired *name* would not).
    for k in 0..n_class {
        if !out.iter().any(|r| r.class == k) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
    }
    Ok(out)
}

/// **The surviving lateral bands, as result faces.**
///
/// `plane_faces` is the plane arrangement's output: its circle bounds are where the cuts are, and
/// taking the boundaries from *there* (rather than from every ⊥ class) is what makes each rim
/// shared by exactly the cap face and the band that meet on it.
///
/// Emitted in `(cylinder class, lower t)` order — `reconstruct` mints handles in face order, so
/// this ordering is the replay contract.
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
    for row in rows {
        for (lo, hi) in bands_of(row, plane_faces, jd)? {
            let (in_own_inside, in_other) = chamber(jd, row, lo, hi, labels)?;
            // ★ **The wall is a boundary face of its own solid, so that solid's membership flips
            // across it** — read inside from the label, and outside is its negation. The
            // counterpart does *not* flip: the gate keeps its boundary faces clear of the lateral, so
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
                surf: ClassIx::Cyl(row.class),
                outer: Bound::Band { lo, hi },
                inner: Vec::new(),
                flip: !keep_in,
            });
        }
    }
    Ok(out)
}

/// The plane classes bounding **this face's** bands, as consecutive `(lo, hi)` pairs ordered by
/// axis parameter.
///
/// The boundaries are the classes that **emitted a circle** for this face's cylinder, plus the
/// classes this face's own rims sit on (found by matching the span's parameters — two ⊥ planes
/// with the same axis parameter *are* one plane, so the match is exact, not a tolerance).
///
/// ★ Circles are gathered per *class*, so a sibling face's boundaries are collected here too — and
/// then clipped away by this face's span. That clip is what keeps two faces of one surface from
/// borrowing each other's cuts.
fn bands_of(
    row: &CylRow,
    plane_faces: &[LocalFace],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<(usize, usize)>, BoolError> {
    let k = row.class;
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
        // band at all (the gate proved the ∥ ones clear of **every band** — its wall rule reads a
        // footprint rectangle per lateral face, so a wall it passed misses each of them), so
        // those are skipped for cause.
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
    lo: usize,
    hi: usize,
    labels: &crate::arrangement::DiskLabels,
) -> Result<(bool, bool), BoolError> {
    let k = row.class;
    let (cyl_bit, other_bit) = match row.side {
        SolidSide::A => (0usize, 2usize), // [A above, A below, B above, B below]
        SolidSide::B => (2usize, 0usize),
    };
    // `above` in the label is the class's **stored** normal side; the band leaves `lo` toward
    // `hi`, i.e. toward increasing axis parameter — which `plus_t_is_above` answers (and where the
    // f64 dot's exactness argument lives).
    let toward_hi = |c: usize| -> bool { crate::planes::plus_t_is_above(&jd.planes[c], &row.def) };
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

    /// **A drilled plate cut by a body flush with its faces.** `Common` of a drilled plate with a
    /// box covering half of it: the two operands share the plate's top and bottom planes, and the
    /// bore's rims end on *holed* faces rather than on disks.
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

    /// **A cylinder standing far away, whose caps happen to land on the box's own planes.** Two
    /// coplanar planes are one class whichever way they sit, so this cylinder's *disks* share a
    /// class with the box's faces while touching nothing — the shape that made the retired seated
    /// rule mistake a shared class for a contact. `Cut` takes nothing away.
    #[test]
    fn a_distant_cap_on_the_boxs_own_plane_removes_nothing() {
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
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
        assert_eq!(out.len(), 1);
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");
    }

    // ---- Seated caps: a cap flush on the other body's face (M6-2 finish) ----
    //
    // ★★ These five are the population the `SeatedCylinderCap` rule refused. What actually makes a
    // seated circle hard is its boundary meeting the counterpart's — and a boundary is either an
    // edge on a plane (parallel to the axis → `WallMeetsLateral`, oblique → `ObliqueCylinderCut`)
    // or another cylinder's rim (→ `CylinderPairContact`). The two fences at the end of this block
    // are those names still standing, so this file says both what was opened and what was not.
    //
    // ★★★ **Measured against the pre-deletion kernel, all seven of these came back
    // `SeatedCylinderCap` — the two fences included.** The rule stood before the wall and the
    // curved-depth rules and answered for them: an overhanging boss and a pair of coplanar-capped
    // cylinders are refused for reasons that have nothing to do with seating, and the sentence a
    // user got named the seating anyway. Deleting it does not only open the five; it lets the two
    // that stay shut say what is actually in the way.

    /// **(a) A through hole whose two caps are flush with the plate's own faces.** The drill neither
    /// overshoots nor stops short: both cap planes coincide with the plate's, and every ⊥ class in
    /// the operation carries faces of both operands.
    #[test]
    fn a_flush_through_drill_bores_the_plate() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a flush through hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(b) A blind hole seated on the face it is drilled from.** Only the lower cap is flush;
    /// the upper one stops inside the plate and closes the bore itself.
    #[test]
    fn a_blind_drill_seated_on_the_plates_base_bores_it() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a seated blind hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 - std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(c) A boss standing on the plate.** `Fuse` with the cylinder's base cap flush on the
    /// plate's top face — the seating a person draws first, and the one the retired rule refused.
    #[test]
    fn a_boss_seated_on_the_plate_fuses_to_it() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a seated boss");
        assert_eq!(out.len(), 1, "one body, not two");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(d) `Common` where the cylinder's caps sit on the box's own faces.** The intersection is
    /// the whole cylinder — every one of its faces is seated or shared.
    #[test]
    fn a_flush_cylinder_meets_the_box_in_itself() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cylinder(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Common, a, b).expect("a flush common");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(i) Drilling the floor of a pocket.** The seated face here is not an original operand
    /// face at all — it is a *floor the previous boolean made*, so the seating arrives through the
    /// arrangement rather than from the modeller. A blind bore down from that floor.
    #[test]
    fn a_drill_seated_on_a_pocket_floor_bores_it() {
        let mut m = Model::new();
        let block = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
        // Open at the top — a tool that stopped inside would leave a void, not a pocket.
        let pocket_tool = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 2.0]),
            Point3::from_array([5.0, 5.0, 7.0]),
        );
        m.rebuild_adjacency();
        let pocketed =
            crate::boolean(&mut m, BoolKind::Cut, block, pocket_tool).expect("pocket")[0];
        // The bore hangs from the pocket floor (z = 2) down into the material below it.
        let drill = m.add_cylinder(
            Point3::from_array([3.0, 3.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, pocketed, drill).expect("a floor bore");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 216.0 - 64.0 - std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **The fence: a boss hanging over the plate's edge is still refused — and the refusal moved
    /// to the name that is true.**
    ///
    /// It used to be `WallMeetsLateral`, because the plate's `x = 4` wall is parallel to the boss's
    /// axis and the wall rule judged the whole infinite band across that plane. But the wall lies
    /// *below* the boss (`t ∈ [−2, 0]` against the boss's span `[0, 1]`), touching only the plane
    /// its base sits in — "the wall meets the lateral surface" was simply false. Reading the span
    /// as the **open** interval the uniform-slab theorem asks for lets the wall through, and the
    /// real obstruction then names itself where it lives: the boss's rim circle crosses the plate
    /// top's boundary segment, which is [`RejectReason::CircleMeetsSegment`] (M6-2b's arcs).
    ///
    /// ★ **Measured, not reasoned:** with the span read as closed instead, this comes back to
    /// `WallMeetsLateral` — that is the whole difference the open reading makes here.
    #[test]
    fn a_boss_overhanging_the_plates_edge_is_still_refused() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let before = m.live_solids.clone();
        let err =
            crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect_err("overhanging boss");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CircleMeetsSegment,
                    ..
                }
            ),
            "{err:?}"
        );
        assert_eq!(m.live_solids, before, "a refused boolean retires nothing");
    }

    /// **The fence: two cylinders with coplanar caps are still refused** — by the name that
    /// describes what is actually hard about them (two curved bodies to classify), not by a rule
    /// about seating.
    #[test]
    fn two_cylinders_with_coplanar_caps_are_still_refused() {
        let mut m = Model::new();
        let a = m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        let b = m.add_cylinder(
            Point3::from_array([5.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let before = m.live_solids.clone();
        let err = crate::boolean(&mut m, BoolKind::Fuse, a, b).expect_err("two curved bodies");
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
        assert_eq!(m.live_solids, before, "a refused boolean retires nothing");
    }

    // ---- Wall faces: the gate asks the boundary, not the infinite plane (M6-2b preparation) ----
    //
    // ★★ The uniform-slab theorem's premise is about the other operand's **boundary**. The gate
    // used to test the wall's infinite *plane*, which is a cheaper sufficient condition — and it
    // refused a whole family the engine serves: a body standing well clear of a bore whose wall
    // plane, extended, happens to pass through it. These lock what that opened and what it did not.

    /// A plate with a bore, fused to a boss standing far away in `x` — whose `y = 12` wall **plane**
    /// stands only 2 from the bore's axis (`r = 3`). The boss's face is 20 away; nothing meets.
    fn plate_bore_and_boss(hole_y: f64, boss_z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, hole_y, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let boss = m.add_cuboid(
            Point3::from_array([28.0, 4.0, boss_z0]),
            Point3::from_array([36.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        (m, holed, boss)
    }

    #[test]
    fn a_boss_whose_wall_plane_crosses_a_distant_bore_still_fuses() {
        let (mut m, holed, boss) = plate_bore_and_boss(10.0, 5.0);
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("the boss fuses");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The `Cut` twin — the same wall planes reach the gate whichever way the operation runs.
    /// The tool is sunk into the plate (`z` from 3) so it removes material: a body merely *resting*
    /// on the top face is a coplanar contact and meets a different guard entirely, which would make
    /// this fixture measure that instead of the wall rule.
    #[test]
    fn the_same_boss_cuts_the_bored_plate() {
        let (mut m, holed, boss) = plate_bore_and_boss(10.0, 3.0);
        let out = crate::boolean(&mut m, BoolKind::Cut, holed, boss).expect("the boss cuts");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 - 8.0 * 8.0 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **Two bores, one wall plane crossing both.** The rule is per (cylinder, class), so a second
    /// cylinder is a second set of questions about the same face — and the face answers for each.
    #[test]
    fn one_wall_plane_may_cross_two_bores() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let mut holed = plate;
        for x in [8.0, 20.0] {
            let hole = m.add_cylinder(
                Point3::from_array([x, 10.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                3.0,
                7.0,
            );
            m.rebuild_adjacency();
            holed = crate::boolean(&mut m, BoolKind::Cut, holed, hole).expect("bore")[0];
        }
        let boss = m.add_cuboid(
            Point3::from_array([28.0, 4.0, 5.0]),
            Point3::from_array([36.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("two bores and a boss");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - 2.0 * std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **The fence: a wall face that really does cross the bore.** The plate's own `y = 20` wall
    /// would clear, so the tool here is a slab whose face runs right across the hole — the wall
    /// rule's true population, and M6-2b's.
    #[test]
    fn a_wall_face_that_really_crosses_the_bore_is_refused() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // A slab covering x ∈ [0, 40]: its y = 12 face passes straight through the bore.
        let slab = m.add_cuboid(
            Point3::from_array([0.0, 12.0, 0.0]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        m.rebuild_adjacency();
        let before = m.live_solids.clone();
        let err = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect_err("a real crossing");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
        assert_eq!(m.live_solids, before, "a refused boolean retires nothing");
    }

    /// **The fence: exact tangency.** The slab's face stands exactly `r` from the axis, so it
    /// touches the lateral surface along one ruling — the boundary case, and `Inside` includes it.
    #[test]
    fn a_wall_face_tangent_to_the_bore_is_refused() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let slab = m.add_cuboid(
            Point3::from_array([0.0, 13.0, 0.0]), // y = 13 is exactly r = 3 from y = 10
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect_err("tangency");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★ **The fence that measures the straight-edge barrier** — the only fixture that does.
    ///
    /// The plate carries a **crosswise** bore, so the plane `x = 30` holds that bore's circular cap
    /// face: an outer loop of one arc and one seam vertex. That plane is also parallel to the
    /// vertical drill's axis and passes within `r` of it, so the face test is reached — and a rule
    /// reading vertices alone would find "every vertex on one side" true of a **single point** and
    /// wave a disk straight through the strip. Asking the edge's curve kind is what stops it.
    #[test]
    fn a_disk_face_on_a_wall_plane_is_undecided_not_waved_through() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 10.0]),
        );
        // A crosswise bore along +X, ending inside the plate at x = 30: its cap lies on x = 30.
        let cross = m.add_cylinder(
            Point3::from_array([-1.0, 10.0, 5.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            2.0,
            31.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, cross).expect("crosswise bore")[0];
        // A vertical drill whose axis stands 1 from the plane x = 30 — inside its radius 3 — and
        // ★ **6 from the crosswise bore's axis**, clear of the radius sum 5, so the pair rule is
        // not what answers here. (At y = 10 the two axes actually meet and this fixture would be
        // measuring `CylinderPairContact` instead — the adjacent proposition.)
        let drill = m.add_cylinder(
            Point3::from_array([31.0, 4.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            12.0,
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, bored, drill)
            .expect_err("a disk face on the wall plane cannot be judged by its vertices");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CylinderGateUndecided,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    // ---- One surface, several lateral faces: the band belongs to the face ----

    /// The shape at the heart of it: a bored plate whose bore's **middle** is cut away, leaving two
    /// disjoint bands on one lateral surface. Returns the model and that solid.
    fn plate_with_a_split_bore(cut_x0: f64, cut_x1: f64) -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let hole = m.add_cylinder(
            Point3::from_array([5.0, 5.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            12.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // Walls 3 from the axis (r = 2), so the wall rule is not what this measures.
        let mid = m.add_cuboid(
            Point3::from_array([cut_x0, cut_x0, 4.0]),
            Point3::from_array([cut_x1, cut_x1, 6.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, bored, mid).expect("the middle cut");
        (m, out[0])
    }

    /// The lid tool the locks below cut with.
    fn add_lid(m: &mut Model) -> Handle<Solid> {
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        lid
    }

    /// How many faces of `s` lie on cylinder surfaces, grouped by surface.
    fn lateral_face_counts(m: &Model, s: Handle<Solid>) -> Vec<usize> {
        let mut counts: std::collections::HashMap<Handle<nacre_geom::Surface>, usize> =
            Default::default();
        let solid = m.solids.get(s);
        for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
            for &fh in &m.shells.get(sh).faces {
                let surf = m.faces.get(fh).surface;
                if matches!(
                    m.surface_truth(surf),
                    nacre_topo::SurfaceTruth::Cylinder { .. }
                ) {
                    *counts.entry(surf).or_default() += 1;
                }
            }
        }
        let mut v: Vec<usize> = counts.into_values().collect();
        v.sort_unstable();
        v
    }

    /// **Making it.** The split bore is an ordinary, correct solid — and nothing measured that
    /// until this defect was found, which is why it is its own test: if building it ever breaks,
    /// the test below (which *uses* it) must not be the one that goes red.
    #[test]
    fn a_bore_cut_across_the_middle_leaves_two_bands_on_one_surface() {
        let (m, s) = plate_with_a_split_bore(2.0, 8.0);
        assert_eq!(
            lateral_face_counts(&m, s),
            vec![2],
            "one lateral surface, two faces"
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, s).expect("props").volume;
        // 1000 − bore(π·4·10) − the middle cube outside the bore(6·6·2 − π·4·2)
        let want = 1000.0 - std::f64::consts::PI * 40.0 - (72.0 - std::f64::consts::PI * 8.0);
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **Using it.** The two-banded solid as an operand. `cyl_rows` used to let the first
    /// lateral face speak for the class, so the second band's span was clipped away and never
    /// emitted — the result shell came back open (`OpenResultShell`), naming a symptom of our own
    /// omission rather than anything about the input.
    #[test]
    fn a_solid_with_two_bands_on_one_surface_can_be_cut_again() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        // ★ The property under repair, counted directly — a volume alone does not say "two bands
        // came out", and this defect is exactly one band going missing.
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![2],
            "both bands survive the cut"
        );
        let p = nacre_props::mass_props(&m, out[0]).expect("props");
        let want = 1000.0
            - std::f64::consts::PI * 40.0
            - (72.0 - std::f64::consts::PI * 8.0)
            - (200.0 - std::f64::consts::PI * 8.0);
        assert!((p.volume - want).abs() < 1e-9, "{} vs {want}", p.volume);
        // ★★ **Area is the sharper oracle for this defect.** A missing band does not change the
        // volume — the boolean simply refuses — but it takes its `2πr·h` out of the surface. The
        // two bands here are `z ∈ [0,4]` and `[6,8]`, so `16π + 8π` of lateral area must be in
        // this number: 640 + 8π once every planar face is counted.
        let want_area = 640.0 + std::f64::consts::PI * 8.0;
        assert!(
            (p.area - want_area).abs() < 1e-9,
            "{} vs {want_area}",
            p.area
        );
    }

    /// ★ **The negative control.** Move the middle cut into a corner, away from the bore: the wall
    /// stays **one** face and the same three steps already worked before the fix. Without this, a
    /// change that merely made *any* third boolean succeed would look like a repair.
    #[test]
    fn a_middle_cut_that_misses_the_bore_leaves_one_band() {
        let (mut m, s) = plate_with_a_split_bore(0.0, 2.0);
        assert_eq!(lateral_face_counts(&m, s), vec![1], "the wall is untouched");
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want =
            1000.0 - std::f64::consts::PI * 40.0 - 8.0 - (200.0 - std::f64::consts::PI * 8.0);
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **The mechanism, read directly.** Two rows, one per lateral face, with the spans the two
    /// bands actually occupy — and the `t` axis here starts at the drill's origin `z = −1`, so the
    /// bands `z ∈ [0,4]` and `[6,10]` are `t ∈ [1,5]` and `[7,11]`.
    #[test]
    fn cyl_rows_gives_one_row_per_lateral_face() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = add_lid(&mut m);
        let setup = plane_index_setup(&m, s, lid).unwrap();
        let rows = cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).expect("rows");
        let spans: Vec<[f64; 2]> = rows
            .iter()
            .map(|r| [r.span[0].to_f64(), r.span[1].to_f64()])
            .collect();
        assert_eq!(
            spans,
            vec![[1.0, 5.0], [7.0, 11.0]],
            "one row per face, in t order"
        );
        assert!(rows.iter().all(|r| r.class == 0), "both on the one surface");
    }

    /// ★★ **The band pass makes nothing in the gap.** This is where the right fix parts company
    /// with the plausible one: merging the two faces' spans into a single `min..max` would put a
    /// band across `z ∈ [4,6]`, where the solid has no lateral face at all.
    #[test]
    fn no_band_is_invented_where_the_solid_has_no_face() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = add_lid(&mut m);
        let (out, ts) = bands(&m, s, lid, BoolKind::Cut);
        let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
        spans.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for (lo, hi) in &spans {
            // Disjoint from the open gap (4, 6): a band either ends at or below 4, or starts at or
            // above 6. ★ Measured: with the spans merged into one `min..max` this fires with
            // "a band 4..6 crosses the gap", which is the whole reason the clause is here.
            assert!(
                *hi <= 4.0 || *lo >= 6.0,
                "a band {lo}..{hi} crosses the gap z 4..6 where there is no face"
            );
        }
        assert!(
            spans.iter().any(|(lo, _)| (*lo - 0.0).abs() < 1e-9),
            "the lower band is emitted: {spans:?}"
        );
        assert!(
            spans.iter().any(|(lo, _)| (*lo - 6.0).abs() < 1e-9),
            "the upper band is emitted too: {spans:?}"
        );
    }

    // ---- The footprint's other axis: along the cylinder (M6-2b preparation) ----
    //
    // ★★ A wall parallel to the axis meets the cylinder in a **rectangle** of the wall's own
    // plane — the strip across, a lateral face's span along — so "does this face miss it" is one
    // question with two separating axes. The gate used to read only the first, and refused every
    // wall that stood clear of the cylinder along its length.
    //
    // ★ With several lateral faces there are several rectangles: the strip is shared, the spans
    // are not. That is why the gap between two bands is passable at all, and it is not a special
    // case bolted on — `bands_of` clips every cut to its own row's span, so no band is ever built
    // in a gap and there is no premise there for a wall to break.

    /// **The case this opened.** The two-banded bore of `plate_with_a_split_bore` has bands at
    /// `z ∈ [0,4]` and `[6,10]` — `t ∈ [1,5]` and `[7,11]`. A tool whose `y = 6` wall stands only 1
    /// from the axis (`r = 2`) crosses the strip, and its face runs the plate's whole width, so the
    /// axis across cannot clear it. Along the axis it is `t ∈ [5.5, 6.5]` — inside the gap, missing
    /// both rectangles.
    #[test]
    fn a_wall_face_in_the_gap_between_two_bands_clears() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let tool = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 4.5]),
            Point3::from_array([10.0, 6.0, 5.5]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the wall clears in t");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![2],
            "the tool passes between the bands and leaves both"
        );
        // 1000 − bore(π·4·10) − the middle cut(6·6·2 − π·4·2) − this tool, which meets the solid
        // over 10×6 minus the 6×4 the middle cut already took, one deep.
        let want = 1000.0
            - std::f64::consts::PI * 40.0
            - (72.0 - std::f64::consts::PI * 8.0)
            - (60.0 - 24.0);
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **The two axes are one judgement, measured as one.** The same tool three ways against
    /// the same two-banded bore, with only the numbers moved: clearing *either* axis is enough,
    /// and clearing neither is the refusal. A rule that had merely gained a second, independent
    /// test would pass (a) and (b) too — what this pins is (c), that the two are OR-ed rather
    /// than each able to wave a face through on its own terms.
    #[test]
    fn either_axis_clears_the_footprint_and_neither_does_not() {
        // (a) Across only: the `y = 6` face sits at `x ∈ [0, 2.5]`, clear of the strip
        //     `x ∈ [3.27, 6.73]`, while its plane still crosses the bore. Along the axis it runs
        //     `t ∈ [4,8]`, straddling both bands. ★ The `x = 2.5` wall must stand clear of the axis
        //     by more than `r`, or *it* becomes the face under test — at `x = 3` it is exactly
        //     tangent and this arm measures that instead (measured: it refuses).
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let across = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 3.0]),
            Point3::from_array([2.5, 6.0, 7.0]),
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, s, across).expect("clear across the strip");

        // (b) Along only: the face spans the full width, so only the gap saves it.
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let along = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 4.5]),
            Point3::from_array([10.0, 6.0, 5.5]),
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, s, along).expect("clear along the axis");

        // (c) Neither: full width *and* straddling both bands.
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let neither = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 3.0]),
            Point3::from_array([10.0, 6.0, 7.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, s, neither).expect_err("crosses both");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    // ---- A blind bore is usable, not just buildable ----
    //
    // ★★ A cylinder's lateral face touching a ⊥ class at its **rim** contributes a `Graze`, the
    // same as any planar wall meeting a class along an edge. It used to contribute nothing, and
    // `edge_mask`'s `Graze > Seated` precedence exists exactly for the corner where the two
    // disagree — a **blind bore's ceiling**. Missing it, the cap's seated rule flipped the wrong
    // label bit, the band's two ends contradicted each other, and the *next* boolean on that solid
    // came back `CylinderGateUndecided`. A through bore has no ceiling and was always fine.

    /// A plate bored to `hole_h` deep, ready to be operated on again.
    fn plate_with_a_bore(hole_h: f64) -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            hole_h,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore");
        (m, out[0])
    }

    /// A boss well clear of the bore (walls 4+ from the axis, r = 3), fused onto the plate.
    fn fuse_a_clear_boss(m: &mut Model, s: Handle<Solid>, z0: f64) -> Vec<Handle<Solid>> {
        let boss = m.add_cuboid(
            Point3::from_array([4.0, 4.0, z0]),
            Point3::from_array([12.0, 6.0, 8.0]),
        );
        m.rebuild_adjacency();
        crate::boolean(m, BoolKind::Fuse, s, boss).expect("the boss")
    }

    #[test]
    fn a_blind_bore_can_be_fused_onto_afterwards() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let out = fuse_a_clear_boss(&mut m, s, 5.0);
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![1],
            "the bore's wall survives"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 + 48.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The `Cut` twin. ★ The tool is sunk into the plate (`z` from 3): a body merely *resting* on
    /// the top face is a coplanar contact that meets `StraightAngle`, a different refusal
    /// entirely, and this fixture would then be measuring that instead of the bore.
    #[test]
    fn a_blind_bore_can_be_cut_afterwards() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let tool = m.add_cuboid(
            Point3::from_array([4.0, 4.0, 3.0]),
            Point3::from_array([12.0, 6.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the cut");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 - 8.0 * 2.0 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **The script a person actually writes**: a blind hole, then a second hole somewhere
    /// else. It goes through the cylinder-pair rule and two lateral surfaces at once — a road the
    /// boss fixtures above do not take.
    #[test]
    fn a_blind_bore_can_be_drilled_beside() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let second = m.add_cylinder(
            Point3::from_array([20.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, second).expect("the second hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![1, 1],
            "two bores, one wall each"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 - std::f64::consts::PI * 45.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **The control**: a through bore and an overshooting one already worked, and must keep
    /// working. Without this, a change that merely made *any* second boolean succeed would look
    /// like the repair.
    #[test]
    fn a_through_bore_is_unaffected() {
        for h in [5.0, 7.0] {
            let (mut m, s) = plate_with_a_bore(h);
            let out = fuse_a_clear_boss(&mut m, s, 5.0);
            assert_eq!(out.len(), 1);
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{:?}",
                nacre_validate::validate(&m)
            );
            let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
            let want = 4000.0 - std::f64::consts::PI * 45.0 + 48.0;
            assert!((v - want).abs() < 1e-9, "h={h}: {v} vs {want}");
        }
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
