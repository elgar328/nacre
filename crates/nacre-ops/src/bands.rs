//! **The lateral face's row, and the vocabulary the chart reads it with.**
//!
//! This file used to *decide* which parts of a cylinder's lateral surface survive a boolean (the
//! band road, M6-2a C4b → M6-2's rulings ladder). Since the D2b cutover (2026-08-30) that decision
//! is `cyl_chart::emit_lateral`'s, cell by cell; D3 deleted the band road (`band_faces`, `bands_of`,
//! `chamber`, `panel_faces`) once the census had held the two emissions equal face by face. What
//! stays here is what the chart reads, one spelling each: [`cyl_rows`]/[`CylRow`] (a face's span
//! is the existence truth at an uncut end), [`face_spans`] (existence at a cut end), [`read_bits`]
//! (which two bits of a label are a cell's chamber), [`keep_for`] (the wall is a boundary of its
//! own solid) and [`axis_param`]. The prose below is the theory those spellings state — the
//! uniform-slab theorem and why the chamber is *read* off the arrangement rather than measured.
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
//! needs: a row's span is the existence truth at an uncut end and its ends are band boundaries
//! (`cyl_chart::boundary_lines`), so no band is ever built in the gap between two of them, and a
//! wall sitting in such a gap breaks no premise — there is none there to break. Reading the spans
//! as one `min..max` would invent both the band and the refusal.
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

#[cfg(test)]
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
/// were then clipped away by the band road and their bands never emitted, which left the result
/// shell open — the boolean came back `OpenResultShell`, naming a symptom of our own omission.
/// The chart reads the same span as the cell's existence at an uncut end (`Chart::read_cell`).
///
/// ★ **Merging the spans into one `min..max` is not the fix either**: it would invent a band across
/// the gap where the solid has no face at all.
pub(crate) struct CylRow {
    /// The [`ClassIx::Cyl`] payload — which lateral surface this face lies on. Several rows may
    /// share it.
    pub(crate) class: usize,
    /// **This face's** own extent in the axis parameter (its two rims).
    pub(crate) span: [Rat; 2],
    /// Which operand this face belongs to.
    pub(crate) side: SolidSide,
}

/// One row per lateral **face**, ordered by `(class, lower t)`.
///
/// ★ That order was the band road's replay contract stated at the row level, and the chart's
/// emitter keeps it (classes ascending, then lines ascending). A class's faces have
/// disjoint spans, so sorting by `(class, span[0])` keeps each class's bands in ascending `t` —
/// the contract generalizes rather than bends. With one face per class it *is* the old class-index
/// order, which is why existing results do not move.
///
/// ★★ **That disjointness used to hold because every lateral face that reached here was a
/// band** — two whole rims and holes. Since E2-2 the tracer states a panel and a chain rim too
/// (`arrangement::lateral_shape`), so two faces of one class *can* share a `t`; the tie then
/// falls to the stable sort's face order — still deterministic, and the corpus so far carries one
/// lateral per class (the cleaning pass merged the pieces one operation earlier). Named here so
/// the day it is exercised is a thing a reader can check, not a surprise.
///
/// A face whose ⊥ range cannot be stated declines by name (`DeclineKind::CylSpan`) — the tracer
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
        let Some(span) = cf.footprint.span else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        // A lateral with no world statement declines the same way a missing span does — the
        // chart states every axis parameter against world planes (`axis_param`).
        if cf.def.is_none() {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        }
        out.push(CylRow {
            class: k,
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

/// The same question asked of a **cylinder** rather than one of its faces — the chart road
/// (`crate::cyl_chart`) has a class and a def but no row, and this rule may not be spelled twice.
pub(crate) fn axis_param(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    def: &nacre_topo::CylinderDef,
) -> Result<Rat, BoolError> {
    param_opt(jd, c, def).ok_or_else(|| reject(RejectReason::WitnessNotRational))
}

pub(crate) fn param_opt(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<Rat> {
    let coeffs = combinatorics::class_coeffs_rat(jd, c)?;
    axis_param_of_plane(&coeffs, def)
}

/// **Does this row's lateral face reach into the interval, at this arc?** — the *existence*
/// question, which no label answers.
///
/// A label states where material is. It is written about a **cell**, and it stays the same whether
/// the cylinder's own face bounds that cell or some other face does. That is enough while a
/// lateral marks a class over its whole circle; it stops being enough the moment a **holed**
/// lateral comes back in as an operand, because then two sectors of one rim can carry a literally
/// identical label and differ only in whether the face is there at all.
///
/// The trace already said which. One rule, no shapes counted:
///
/// | this face's mark on the arc | reaches into the interval |
/// |---|---|
/// | [`SegKind::Transversal`] | **yes** — the face passes through the plane here |
/// | [`SegKind::Graze`] | only when `body_above` names the interval's side |
/// | nothing | **no** — the face stops short of this arc |
///
/// `band_is_above` is the interval's side of this rim's plane in that plane's **stored** frame —
/// the very argument [`read_bits`] takes, and produced by the same `toward_hi` the caller already
/// holds. There is no second derivation of "which way is the band" here, deliberately: this
/// ladder has been bitten four times by a sign re-derived one call away from its twin.
///
/// ★ [`SegKind::Seated`] is skipped because it is a *planar* face's word — see [`ArcLabel::marks`],
/// where the type is what rules a lateral out, not a convention.
///
/// ★★ **The rule is not about holes.** An outer rim grazes too (`circle_on_class` says
/// `Grazes { body_above: up }` at the face's own end), and there the interval on the face's side
/// gets `true` from this same test — one sentence about every rim. What is out of reach is only
/// an **uncut** outer rim: it never becomes an [`ArcLabels`] entry at all (a whole circle goes to
/// `DiskLabels` and the band road), so today only cut rims ask here.
///
/// ☑ **How often each row actually fires** (whole binary, production calls only, measured on the
/// band road 2026-08-29): `Transversal` **265** · `Graze` **13**, of which **1** reaches and **12**
/// do not · nothing at all **0**. The chart reads it once per cut end (`cyl_chart::census`'s
/// `arcs_read`, 756 over the lib suite in D3).
/// The twelve are the six dropped sectors read at both rims, which is the cross-check that the
/// sector census and this one describe the same events.
pub(crate) fn face_spans(
    r: &crate::arrangement::ArcLabel,
    side: SolidSide,
    band_is_above: bool,
) -> Result<bool, BoolError> {
    use crate::arrangement::SegKind;
    let mut answer: Option<bool> = None;
    for (_, kind) in r.marks.iter().filter(|(s, _)| *s == side) {
        let reaches = match *kind {
            // A planar face's word — see `ArcLabel::marks`; a tangent ruling is a seated face's
            // too (cell ⑩), and never a rim arc's mark.
            SegKind::Seated { .. } | SegKind::Tangent { .. } => continue,
            SegKind::Transversal { .. } => true,
            SegKind::Graze { body_above } => body_above == band_is_above,
        };
        // ★★★ **Two of this solid's faces meeting at one arc and *disagreeing*: refused, and the
        // refusal is a placeholder.** It takes one cylinder class carrying several faces that
        // share a rim — reachable geometry (a split bore whose two bands touch), just not
        // reachable today. The answer that day is almost certainly `.any()`: if any of the
        // solid's faces reaches into the interval, a face is there. It is not written that way
        // now because it cannot be measured, and this ladder has twice shipped an unmeasured rule
        // that turned out wrong. Marks that **agree** decide nothing by themselves, so they are
        // taken: refusing there would be refusing a case with no guess in it.
        if answer.replace(reaches).is_some_and(|prev| prev != reaches) {
            return Err(reject(RejectReason::CylinderFaceUndecided));
        }
    }
    Ok(answer.unwrap_or(false))
}

/// The one spelling of "which two bits of a label are this cell's chamber": the row's own
/// solid's bit and the counterpart's, on the side of the plane the cell occupies. Read for a disk
/// end and an arc end alike (`Chart::read_cell`), so the two cannot drift.
pub(crate) fn read_bits(
    l: &crate::arrangement::Label,
    side: SolidSide,
    band_is_above: bool,
) -> (bool, bool) {
    let (cyl_bit, other_bit) = match side {
        SolidSide::A => (0usize, 2usize), // [A above, A below, B above, B below]
        SolidSide::B => (2usize, 0usize),
    };
    let i = usize::from(!band_is_above);
    (l[cyl_bit + i], l[other_bit + i])
}

/// The one spelling of the band's keep decision — the wall is a boundary face of its own solid,
/// so that solid's membership flips across it while the counterpart's does not.
pub(crate) fn keep_for(kind: BoolKind, side: SolidSide, in_own: bool, in_other: bool) -> bool {
    match side {
        SolidSide::A => crate::arrangement::keep(kind, in_own, in_other),
        SolidSide::B => crate::arrangement::keep(kind, in_other, in_own),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planes::{PlaneSetup, plane_index_setup};
    use nacre_math::{Point3, Vector3};
    use nacre_store::Handle;
    use nacre_topo::Surface;
    use nacre_topo::{Model, Solid};

    /// A minted vertex sits within `1e-12` of its definition — asserted from the cache where the
    /// cache proves it, and **measured here** where it does not.
    ///
    /// ★ The second half used to read a stored residual. The cache no longer stores one, and the
    /// honest replacement is not "skip those" but "measure the same distance the residual was":
    /// the point against the surfaces its definition names. Dropping to a bare `panic!` for the
    /// unproven variants would have turned this from a measurement into a restatement of which
    /// variant the funnel chose.
    fn knowledge_is_tight(m: &Model, h: Handle<nacre_topo::Vertex>) {
        match m.vertex_cache(h) {
            nacre_topo::PointCache::Bounded { bound, .. } => assert!(
                bound
                    .iter()
                    .all(|b| nacre_scalar::Mag::lt(*b, nacre_scalar::Mag::of(1e-12))),
                "realized within rounding: {bound:?}"
            ),
            nacre_topo::PointCache::Ceiling { coord }
            | nacre_topo::PointCache::Unrealized { coord } => {
                for sh in m.vertex(h).carriers() {
                    let d = m.surface_cache(sh).distance(*coord);
                    assert!(
                        d < 1e-12,
                        "unrealized, but {d:e} from surface {}",
                        sh.index()
                    );
                }
            }
        }
    }

    /// **The existence rule's whole truth table** — [`face_spans`] against marks written by hand.
    ///
    /// The production lock reaches this through an entire boolean, and a boolean exercises two of
    /// the rows below: a `Transversal` that reaches and a `Graze` pointing away. The rule is one
    /// sentence about every rim, not a description of that fixture, so the rest are stated here —
    /// including the two that can only be reached from geometry the road does not build yet.
    #[test]
    fn face_spans_reads_the_trace_not_the_label() {
        use crate::arrangement::{ArcLabel, SegKind};
        let n = combinatorics::NodeId::ThreePlane([0, 1, 2]);
        let row = |marks: Vec<(SolidSide, SegKind)>| ArcLabel {
            ends: [n, n],
            // Deliberately the label that says "material everywhere": nothing below may read it.
            label: [true; 4],
            marks,
        };
        let (a, b) = (SolidSide::A, SolidSide::B);
        let spans = |marks, above| face_spans(&row(marks), a, above);
        // A face running through the plane is on both sides of it.
        for above in [true, false] {
            assert!(spans(vec![(a, SegKind::Transversal { mat: 1 })], above).unwrap());
        }
        // A face whose boundary stops at the arc is on exactly the side it occupies.
        for body_above in [true, false] {
            for above in [true, false] {
                assert_eq!(
                    spans(vec![(a, SegKind::Graze { body_above })], above).unwrap(),
                    body_above == above,
                    "graze body_above={body_above} against band above={above}"
                );
            }
        }
        // The counterpart's marks are not this row's face, whatever they say.
        assert!(!spans(vec![(b, SegKind::Transversal { mat: 1 })], true).unwrap());
        // A planar face's seated rim says nothing about the lateral — see `ArcLabel::marks`.
        assert!(!spans(vec![(a, SegKind::Seated { body_above: true })], true).unwrap());
        // No mark at all: the face stops short of this arc.
        assert!(!spans(Vec::new(), true).unwrap());
        // Two of this solid's faces on one arc: agreeing decides, disagreeing refuses by name.
        assert!(
            spans(
                vec![
                    (a, SegKind::Graze { body_above: true }),
                    (a, SegKind::Transversal { mat: 1 }),
                ],
                true,
            )
            .unwrap()
        );
        assert!(matches!(
            spans(
                vec![
                    (a, SegKind::Graze { body_above: false }),
                    (a, SegKind::Transversal { mat: 1 }),
                ],
                true,
            ),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderFaceUndecided,
                ..
            })
        ));
    }

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
            cyls,
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
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
            m,
            kind,
            a,
            b,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the drill population traces");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("cylinder rows");
        let out = crate::cyl_chart::emit_lateral(kind, &jd, cyls, &plane_faces, &curved, &rows)
            .expect("the lateral faces emit");
        // The classes' **z**, not their axis parameter: `t` is measured from the cylinder's own
        // origin along its raw `dir`, so a drill starting at z=−1 puts the box's cap at t=1. The
        // assertions read in world z, which is the vocabulary the fixtures are written in.
        let ts = (0..geom.len())
            .filter_map(|c| {
                let def = &cyls[0].def;
                let t = param_opt(&jd, c, def)?;
                // The class's world z, via the axis point at that parameter.
                let (o, m) = (def.origin(), def.dir());
                let z = o[2].checked_add(t.checked_mul(m[2])?)?;
                Some((c, z.to_f64()))
            })
            .collect();
        (out, ts)
    }

    /// A band's two ends as world `z` — the assertion vocabulary.
    fn ends(lf: &LocalFace, ts: &[(usize, f64)]) -> (f64, f64) {
        let Bound::Band { lo, hi } = &lf.outer else {
            panic!("a band face bounds a band");
        };
        let (Some(lo), Some(hi)) = (lo.circle(), hi.circle()) else {
            panic!("a band's rims are whole circles here");
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
            .map(|sh| m.shell(sh).faces.len())
            .sum::<usize>();
        assert_eq!(faces, 7, "4 walls + 2 drilled caps + the bore's wall");
        // ★ **Genus 1, re-derived rather than quoted.** An earlier plan wrote `V10−E15+F7−L2`
        // from memory; the counts are measured here and what is asserted is the relation
        // (`χ = V − E + F − L = 2(S − G)`, so one shell with one through hole gives `χ = 0`).
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
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
            .map(|fh| m.face(*fh).inner.len() as i64)
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
    // edge on a plane (parallel to the axis → the wall rule, oblique → asked to miss every lateral
    // face since cell ⑩, else `ObliqueCylinderCut`) or another cylinder's rim (→ proved apart per
    // face pair since cell ⑩, else `CylinderPairContact`). ★ The wall rule is no longer a *fence*:
    // it records a crossing (cell ③) or a tangency (cell ⑥) and the roads behind it answer, so
    // what still stands at the end of this block is the oblique cut that does meet a face and the
    // cylinder pair whose faces do meet.
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

    /// One rulings-road boolean on the 40×40×20 plate, end to end on the production road:
    /// build, retire the operands, validate clean, answer the exact volume, tessellate
    /// watertight. The volume oracle is what finally *measures* the ladder's sign roster —
    /// a flipped panel winding, disk-side selector, ruling turn or chord sense moves it.
    fn through_boss_builds(kind: BoolKind, base: [f64; 3], want: f64) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            50.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, kind, plate, boss).expect("the rulings road builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the panels tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
    }

    /// ★ **The through-boss builds** — the rulings ladder's milestone: a boss standing through
    /// the plate's wall (axis exactly on the `x = 40` plane), the app-measured refusal that
    /// opened M6-2's remainder. Fuse keeps the outer half (plate + the boss outside it), cut
    /// carves the notch, common keeps the inner half-cylinder — each an exact closed form, and
    /// the three exercise both complementary θ-sectors.
    #[test]
    fn a_through_boss_fuses() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Fuse, [40.0, 20.0, -10.0], 32000.0 + 1000.0 * pi);
    }

    #[test]
    fn a_through_boss_cuts_a_notch() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Cut, [40.0, 20.0, -10.0], 32000.0 - 250.0 * pi);
    }

    #[test]
    fn a_through_boss_common_is_the_inner_half() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Common, [40.0, 20.0, -10.0], 250.0 * pi);
    }

    /// ★ A boss on the plate's **corner** — its axis on *two* wall planes, both recorded, each
    /// circle cut by both walls: the rim pairing (`(wall class, root)`) and the θ-panel
    /// machinery quarter the lateral. Measured to assemble in the lift probe; the volume pins
    /// it (plate + three quarters of the cylinder outside).
    #[test]
    fn a_corner_boss_fuses() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Fuse, [40.0, 40.0, -10.0], 32000.0 + 1125.0 * pi);
    }

    /// ★★★★★ **Both of this fence's exclusions are gone, and the tangent one took the longest
    /// because its old sentence was *true*.** The **offset** crossing (`0 <` distance `< r`) stood
    /// here on «walks to `OpenResultShell`», a measurement from before the region emitter, and
    /// builds exactly since cell ③. The **tangent** wall (distance exactly `r`) stood on «assembles
    /// a volume-correct zero-thickness pinch `validate` cannot see» — and cell ⑥ measured that this
    /// is *exactly right*: with the gate passing it, `Cut` returns `Ok`, `validate` is clean, and
    /// nothing in the kernel sees the contact. So the answer was never to keep refusing at the
    /// gate; it was to give that pinch a judge (`boolean::tangency_reject`), which is what this row
    /// now exercises — the operation decides, and the ones that do not pinch build.
    ///
    /// ★ A third row lived here until the D2b cutover: the **half-height** boss, whose upper cap
    /// sits inside the plate's material — the band road's `chamber` had no sector answer for that
    /// end and refused it `RulingBoundNotYet`. The chart reads it (a band below the plate, the
    /// outer sector beside it), so it builds now — [`Self::a_half_height_boss_builds`].
    #[test]
    fn one_chord_formula_covers_the_offset_wall_and_its_tangent_limit() {
        let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
        let pi = std::f64::consts::PI;
        // ★★★★★ **One expression, both rows — which is the evidence that the rule is general.**
        // The boss's part inside the plate is the disk less the `x > 40` circular segment, over the
        // plate's height; the tangent row is that same expression at `d = r`, where `seg(5, 5) = 0`
        // and the whole disk is inside. Nothing here is copied from an engine run: the tangent
        // volume is *derived* by walking the offset row's formula to its limit.
        let want = |d: f64| 32000.0 + 1250.0 * pi - (25.0 * pi - seg(d, 5.0)) * 20.0;
        for (base, d) in [([38.0, 20.0, -10.0], 2.0), ([35.0, 20.0, -10.0], 5.0)] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 20.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(base),
                Vector3::from_array([0.0, 0.0, 1.0]),
                5.0,
                50.0,
            );
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .unwrap_or_else(|e| panic!("{base:?}: {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), 1, "{base:?}");
            assert!(nacre_validate::validate(&m).is_empty(), "{base:?}");
            let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
            assert!((v - want(d)).abs() < 1e-9, "{base:?}: {v} vs {}", want(d));
        }
    }

    /// ★★★★★ **The half-height boss builds — the capability the D2b cutover opened.** A boss
    /// through the plate's wall whose upper cap (z = 10) sits inside the plate: below the plate
    /// the lateral is a whole band, beside it only the outer sector survives, and the cap is a
    /// half-disk. The band road refused this end (`RulingBoundNotYet`); the chart's cells read
    /// it off the same labels. Volumes derived, not copied: the boss outside the plate is the
    /// outer half-cylinder over the full height (`π·25·20/2 = 250π`) plus the inner half below
    /// the plate (`π·25·10/2 = 125π`); the boss inside the plate is the inner half over
    /// `z ∈ [0, 10]` (`125π`).
    ///
    /// The lateral is **one** face (D4): the band and the outer panel merge into a band whose
    /// upper rim is a **wrapping chain** (the z = 0 inner arc, a ruling, the z = 10 outer arc,
    /// a ruling) — the cleaning pass used to abstain here, having no `Bound` for it. ★ Both
    /// senses: the boss with its cap inside the plate (`z0 = −10`, a `hi` chain) and its mirror
    /// with its base inside (`z0 = 10`, a `lo` chain) — the walk's rotation rule is measured
    /// on each rather than assumed symmetric. Same volumes by symmetry.
    #[test]
    fn a_half_height_boss_builds() {
        let pi = std::f64::consts::PI;
        for (z0, kind, want, lateral) in [
            (-10.0, BoolKind::Fuse, 32000.0 + 375.0 * pi, vec![1]),
            (-10.0, BoolKind::Cut, 32000.0 - 125.0 * pi, vec![1]),
            (-10.0, BoolKind::Common, 125.0 * pi, vec![1]),
            (10.0, BoolKind::Fuse, 32000.0 + 375.0 * pi, vec![1]),
            (10.0, BoolKind::Cut, 32000.0 - 125.0 * pi, vec![1]),
            (10.0, BoolKind::Common, 125.0 * pi, vec![1]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 20.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array([40.0, 20.0, z0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                5.0,
                20.0,
            );
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, kind, plate, boss)
                .unwrap_or_else(|e| panic!("{kind:?} z0 {z0}: the half-height boss builds: {e:?}"));
            assert_eq!(out.len(), 1, "{kind:?}: one solid");
            m.rebuild_adjacency();
            assert_eq!(m.live_solids, out, "{kind:?}: the operands retired");
            let issues = nacre_validate::validate(&m);
            assert!(issues.is_empty(), "{kind:?}: {issues:?}");
            let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
            assert!((v - want).abs() <= 1e-9 * want, "{kind:?}: {v} vs {want}");
            assert_eq!(
                lateral_face_counts(&m, out[0]),
                lateral,
                "{kind:?}: lateral faces"
            );
            let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                .expect("the half-height boss tessellates");
            let mut uses: std::collections::HashMap<(u32, u32), usize> =
                std::collections::HashMap::new();
            for (_, tri) in mesh.triangles.iter() {
                for k in 0..3 {
                    let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                    *uses.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            assert_eq!(
                uses.values().filter(|&&n| n != 2).count(),
                0,
                "{kind:?}: the mesh is watertight"
            );
        }
    }

    /// **The straddling boss builds.** The M6-2b milestone fence: a boss hanging over the
    /// plate's edge — its rim circle cut by the plate top's boundary segment — fuses into one
    /// valid solid. The ladder this closes, in the names its rungs wore: `WallMeetsLateral`
    /// (a reason since retired — the wall rule read the closed span) → `CircleMeetsSegment` →
    /// `ArcBoundNotYet` (the
    /// stopper, walked from the class arrangement to the very end of the assembly) → built.
    ///
    /// The volume and the mixed-loop integrals are pinned by
    /// [`Self::a_cut_rim_boolean_builds_a_complete_solid`]; here the result answers the two
    /// whole-model judges: `validate` (whose winding check reads the arcs as witnesses now) and
    /// the tessellation (watertight, arcs sampled as sub-arcs).
    #[test]
    fn a_boss_overhanging_the_plates_edge_builds() {
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
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        assert_eq!(out.len(), 1, "one fused solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the arcs tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "the mesh is watertight");
    }

    /// ★★ **A boss standing over a bore, its whole outline *inside* the rim.** No segment crosses
    /// the circle, so the circle keeps its closed cell — what the shape really asks is that the
    /// **disk cell host a polygon**, and `nest_cells` refused to (a `debug_assert` said a disk
    /// hosts nothing, for a reason that was already wrong: the guard actually keeping the
    /// population out was `circles_meet_no_segment`).
    ///
    /// ★★★ **The population immediately found a live defect in the nesting predicate.** With the
    /// disk allowed to host, the *circle*'s contour came back as a hole of the **square** — the
    /// footprint `[7,9]×[9,11]` straddles the axis `(8,10)`, so "the circle's centre is inside the
    /// ring" is true in the nesting that does not hold. One witness point says *whether* two
    /// disjoint loops nest, never *which way*; the other direction is what settles it
    /// (`cell_in_cell`).
    #[test]
    fn a_segment_inside_the_rim_builds_the_boss_over_the_hole() {
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
        // Footprint `[7,9] × [9,11]` sits wholly inside the rim `(8,10)`, `r = 3`; `z ∈ [5,8]`
        // leaves the bore's span `t ∈ [1,6]` clear along the axis, so the wall rule passes it on.
        let boss = m.add_cuboid(
            Point3::from_array([7.0, 9.0, 5.0]),
            Point3::from_array([9.0, 11.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a boss over a bore");
        // ★ The boss's whole footprint is inside the rim, so it stands over the **hole** and
        // touches no material: two bodies, not one. That is forced by the *geometry* — the
        // footprint is strictly inside the rim, so there is no material to reach — and no longer
        // by a refusal: the wall rule used to make the argument for us, and since cells ③ and ⑥ it
        // serves that population instead of fencing it off.
        assert_eq!(out.len(), 2, "the boss touches nothing");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [
            2.0 * 2.0 * 3.0,
            40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0,
        ];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "the boss and the bored plate: {vols:?} vs {want:?}"
        );
        // The plate keeps its through hole (χ = 0), the boss is a box (χ = 2).
        let mut chi: Vec<i64> = out
            .iter()
            .map(|&s| {
                let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
                v_n - e_n + f_n - l_n
            })
            .collect();
        chi.sort_unstable();
        assert_eq!(chi, vec![0, 2], "genus 1 and genus 0");
    }

    /// **The tangency is a *slit*, not a cut — and that is why no combinatorial check can see it**
    /// (cell ⑤). Two clauses, because either alone passes vacuously:
    ///
    /// 1. **No vertex within `1e-9` of the touch.** The arrangement names the point exactly
    ///    (`NodeId::pierce(.., QuadRoot::Double)`) and then *discards* it — a tangency touches
    ///    without separating, so the circle keeps its closed cell. `nonmanifold_vertices` and the
    ///    Euler count read topology, and there is none here: **their silence is not evidence**
    ///    about this point, in either direction.
    /// 2. **The circle that touches is still whole** — a closed `[v, v]` rim edge passing through
    ///    the point. Without this the first clause cannot tell "the circle was never cut" from
    ///    "the circle was cut and its vertex landed elsewhere".
    ///
    /// What *is* the evidence that the solid is sound: the link of the boundary at the touch is a
    /// single circle (the pinched face's two lobes are joined around through the neighbouring
    /// curved face), the material is locally one piece, and a second kernel returns a body of the
    /// same volume and area (`nacre-oracle`'s `a_segment_tangent_to_a_rim_is_a_body_to_occt` and
    /// `a_rim_tangent_to_a_plate_top_is_a_body_to_occt`, measured against an ε-twin whose tangency
    /// is broken). ⇒ **the surface is a 2-manifold there; only the face is pinched.**
    fn the_touch_is_a_slit(m: &Model, s: Handle<Solid>, at: [f64; 3]) {
        let p = Point3::from_array(at);
        let sol = m.solid(s).clone();
        let (mut vertices_at, mut whole_circle_through) = (0usize, 0usize);
        let mut seen = std::collections::HashSet::new();
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let face = m.face(fh);
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        if !seen.insert(he.edge) {
                            continue;
                        }
                        let e = m.edge(he.edge);
                        for &vh in e.vertices.iter() {
                            if m.vertex_point(vh).distance(p) <= 1e-9 {
                                vertices_at += 1;
                            }
                        }
                        // A whole circle is the closed `[v, v]` rim spelling; a cut one is arcs.
                        if e.vertices[0] == e.vertices[1]
                            && m.edge_curve(he.edge).distance(p) <= 1e-9
                        {
                            whole_circle_through += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(vertices_at, 0, "a tangency mints no vertex: {at:?}");
        assert!(
            whole_circle_through > 0,
            "the touching circle is still whole (a closed rim edge through {at:?})"
        );
    }

    /// **A tangency builds.** The boss's `x = 11` edge is exactly tangent to the rim `(8,10)`,
    /// `r = 3`, so the locator's quadratic has a double root — the third arm of `CylinderMeet`
    /// that can reach here, and the only one whose point carries no radical (`(11, 10, 5)`,
    /// rational).
    ///
    /// ★★ **And a touch is not a break**: the circle keeps its closed cell, so nothing about the
    /// arrangement had to change for this to work. It was refused only because the guard's
    /// sentence was wider than its proposition — measured before the narrowing landed, with the
    /// whole guard off, this already produced exactly the solid asserted below.
    #[test]
    fn a_segment_tangent_to_the_rim_builds() {
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
        let boss = m.add_cuboid(
            Point3::from_array([11.0, 8.0, 5.0]),
            Point3::from_array([15.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a tangent edge");
        assert_eq!(out.len(), 1, "the boss sits on material");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0 + 4.0 * 4.0 * 3.0;
        assert!(
            (v - want).abs() < 1e-9,
            "plate less bore plus boss: {v} vs {want}"
        );
        // The through bore survives the fuse, so the one shell still has one handle.
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            0,
            "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        // ★ **This χ is the genus's lock, not the tangency's.** A pinch *that has a vertex* makes
        // χ odd — what `check_result_topology` reads — but a tangency mints no vertex at all
        // (below), so an even χ is compatible with a pinch and with none: it decides nothing
        // here. Deleting this assertion would lose the genus, not the manifold claim.
        m.rebuild_adjacency();
        the_touch_is_a_slit(&m, out[0], [11.0, 10.0, 5.0]);
        // ★★★★★ **And the solid meshes — across a bridge, since cell ⑦c.** The boss's footprint
        // touches the bore's rim at exactly one point, so the plate's top face has two inner
        // loops meeting there: its interior is pinched. Until ⑦c that was `SelfTouchingBoundary`
        // — a statement about the *tessellator* (`validate` is clean above and the volume is
        // exact). Now the touching sample is put into the straight edge it lies on, the two
        // holes are spliced into one at that point, and the sweep orders the twins symbolically.
        // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
        // census that only agreed with itself would be measuring when it looks rather than what
        // came out. Both spellings say the same thing here.
        m.rebuild_adjacency();
        crate::tests::mesh_covers_faces("segment tangent to the rim", &m, &out);
    }

    /// ★★ **The fence post: a crossing that lands exactly on a segment's endpoint.** The boss's
    /// corner `(8,13)` sits on the rim `(8,10)`, `r = 3`, so both edges leaving it meet the circle
    /// *at* their own end — the one place the in-segment test's inequality can be open or closed,
    /// and the two spellings give different answers.
    ///
    /// **Measured, both ways:** counting `Zero` as inside names the corner `(8, 13, 5)`; excluding
    /// it drops every root and falls back to the containment witness — a `Segment`, which would
    /// say "this edge lies inside the circle" about an edge that runs *outward* from a single
    /// touching point. Inclusive is the true sentence, and it is the one the arc split will need.
    #[test]
    fn a_crossing_on_a_segments_endpoint_is_inside_it() {
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
        let boss = m.add_cuboid(
            Point3::from_array([8.0, 13.0, 5.0]),
            Point3::from_array([12.0, 17.0, 8.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect_err("corner on rim");
        // ★★ **One point, two names** — and the split says so by name. `(8,13,5)` is the boss's
        // corner, so the arrangement already holds it as a three-plane vertex; the rim crossing
        // names the *same* point as a `NodeId::Pierce`. The DCEL keys vertices by name, so shipping
        // both would put two vertices where there is one — folding them is its own step, and until
        // then `CoincidentNodes` ("two names for one point") is the true sentence.
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CoincidentNodes,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★ **The population where the naming rule actually fires — and the first non-`+Z`
    /// cylinder in this repository.**
    ///
    /// Every `add_cylinder*` call in the suite states the axis `[0, 0, 1]` (measured: 87 of them),
    /// and that is not an accident of taste — `add_cuboid` pushes faces `[−Z, +Z, −Y, +Y, −X, +X]`
    /// and `build_prism` pushes "base cap, top cap, then walls", so a `+Z` cylinder's circle always
    /// lives on a class interned **before** the wall that crosses it. Measured over the whole ops
    /// suite: 463 of the 2060 `(class, wall)` pairs the gate examines are descending, but **all
    /// nine that reach the locator are ascending**. So `NodeId::pierce`'s canonicalization —
    /// re-sorting the pair and restating the root with it — has never once been exercised by a
    /// production path.
    ///
    /// Turning the boss onto `+X` inverts it: the rim lives on the box's `x = 4` face (class 5,
    /// interned last) and the segment it crosses is on the `z = 2` cap (class 1). The pair swaps.
    ///
    /// ★ **A tangency, so there is exactly one root** — and since 2026-08-21 that also means the
    /// circle is not separated, so this **builds** rather than refusing. The single root is still
    /// what makes the naming rule visible here (it is what the fixture was written for), and the
    /// witness it once asserted is now the *split point that never happens*. That matters: the
    /// straddling fixtures deliberately
    /// accept either of their two crossings, and a fixture copied from that template would be green
    /// whether or not the root rule is right. Here `disc = 0` — the rim (`x = 4`, centre
    /// `(4, 2, 1.5)`, `r = 0.5`) touches `z = 2` at the single point `(4, 2, 2)`, solved from the
    /// fixture's own numbers — and the name it must take is `QuadRoot::Double`, which a swap must
    /// **not** toggle.
    #[test]
    fn a_turned_boss_tangent_to_the_plate_top_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, 1.5]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a turned boss");
        assert_eq!(
            out.len(),
            1,
            "the boss's base disk sits on the plate's face"
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!(
            (v - want).abs() < 1e-9,
            "plate plus the whole boss: {v} vs {want}"
        );
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        // ★ The genus's lock, not the tangency's — see the note in the fixture above.
        m.rebuild_adjacency();
        the_touch_is_a_slit(&m, out[0], [4.0, 2.0, 2.0]);
        // ★★★★★ **And the solid meshes — the same pinch, spelled inner-to-outer.** The boss's
        // base circle is tangent to the plate's `z = 2` edge at `(4, 2, 2)`, so the `x = 4` face's
        // hole touches its own outer ring at one point and the face's interior is pinched there.
        // `validate` is clean above and the volume is exact; since cell ⑦c the touch is bridged
        // and the face's triangles are held to its exact area.
        // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
        // census that only agreed with itself would be measuring when it looks rather than what
        // came out. Both spellings say the same thing here.
        m.rebuild_adjacency();
        crate::tests::mesh_covers_faces("turned boss tangent to the plate top", &m, &out);
    }

    /// **The turned boss builds — and its crescent is the winding witness's red switch.**
    ///
    /// The near cap's circle is cut by two different plate edges (top and corner), so this
    /// result carries the arc-dominated crescent face whose chord Newell reads **backwards**
    /// (`cos = −1`, the measured wall) — `validate == []` here is what pins `loop_winding`'s
    /// segment witnesses. ★ The old fence's proposition — a root that fails to follow its pair
    /// through the sort names the wrong crossing — did not retire with the reject: the mint
    /// fence asserts the pierce vertices sit **on the derived crossings**, and a wrong root
    /// moves the minted point itself.
    #[test]
    fn a_turned_boss_over_the_plates_corner_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 0.25, 2.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        assert_eq!(out.len(), 1, "one fused solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-12, "{v} vs {want}");
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }

    /// **The audit and the boolean say the same thing about an arc class — and it says what the
    /// arrangement produced.**
    ///
    /// ★★★★ Three things are pinned here and none had a fence before.
    ///
    /// **One: the two copies of the per-class pipeline agree.** `frame_audit` re-runs pass A and B
    /// inline rather than calling `arrange`, so a split hoisted into one and not the other compiles
    /// perfectly and makes the audit report *no failure* for an input the boolean refuses —
    /// `decline_to_reject`'s doc calls that "the worst possible time to be lying". The existing
    /// guard (`the_audit_does_not_invent_failures`) cannot see it: its fixture is deliberately one
    /// that fails **outside** the class pipeline.
    ///
    /// **Two: what the arrangement produced on an arc class.** The stopper stands after every
    /// arrangement stage and swallows their answer, so without this the whole arc population is
    /// measured by a probe and the probe is deleted before the commit. The numbers are **derived,
    /// not read back**: the
    /// boss's rim crosses one edge of the plate's top face twice, cutting the circle into **2**
    /// arcs and the plane into **4** cells — `rect−disk`, `rect∩disk`, `disk−rect` and the outside.
    /// Only the outside winds `−1`, so there is **1** root group; and it shares nodes with all
    /// three `+1` cells (the plate's corners, the two crossings), so `cell_in_cell` answers "not
    /// comparable" every time and there are **0** holes. Exactly **one** class is cut: the boss
    /// stands *on* the top face, so only its base circle lies in a plane of the plate.
    ///
    /// **Three: the faces it emitted, as coordinates** — and this is the one that sees a wrong
    /// `sense`. Flipping the sense the split carries onto its sub-segments attaches the arcs to
    /// the wrong cells, and **every number in part Two stays put** (the cells still read four with
    /// one `−1`, the roots one, the holes zero, the sorted labels identical) — that blindness was
    /// derived before it was measured, and it held for four rungs. The emitted ring is where it
    /// finally shows, as the **exact reverse**, which is why rotation is normalized here and
    /// reversal is not. See `ClassAudit::outer_rings` for both arguments and for why the red probe
    /// is read on the **turned** boss: the straddling one is blind to the flip *and* its 2-node
    /// ring cannot express a reversal at all.
    ///
    /// ★ **Two earlier claims here were wrong and are recorded rather than quietly dropped**: that
    /// this fence is *green* under that flip (it is red — measured, three times), and that the
    /// flip's first reader is `label_cells`' keep decision (refuted a rung earlier: every
    /// order-independent summary of the labels is identical; it is `emit_faces`).
    #[test]
    fn the_audit_and_the_boolean_agree_about_an_arc_class() {
        // `y = 0` meets the turned boss's circle (centre `(y,z) = (0.25, 2)`, `r = 0.5`) at
        // `(z - 2)² = 0.1875`. The one irrational coordinate in either fixture, and the reason the
        // ring comparison is `near()` rather than `==`.
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        for (origin, axis, n_out, want_rings) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                // The plate's top face carries this class, so the ring is CCW about `+z`.
                [0.0, 0.0, 1.0],
                vec![
                    // `rect - disk`: the plate's top, three corners and the inner arc.
                    vec![
                        [0.0, 0.0, 2.0],
                        [4.0, 0.0, 2.0],
                        [4.0, 1.5, 2.0],
                        [4.0, 2.5, 2.0],
                        [4.0, 4.0, 2.0],
                        [0.0, 4.0, 2.0],
                    ],
                    // `disk - rect`: the overhang's underside — chord and outer arc.
                    vec![[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
                ],
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                // The plate's `x = 4` face carries this one, so CCW is about `+x`.
                [1.0, 0.0, 0.0],
                vec![
                    // `rect - disk`: the plate's right face, the disk biting its top-left corner.
                    vec![
                        [4.0, 0.0, 0.0],
                        [4.0, 4.0, 0.0],
                        [4.0, 4.0, 2.0],
                        [4.0, 0.75, 2.0],
                        [4.0, 0.0, s],
                    ],
                    // `disk - rect`: the boss's base cap outside the plate. The corner `(4,0,2)` is
                    // `0.25` from the circle's centre, so it sits *inside* the disk and is one of
                    // this ring's three nodes.
                    vec![[4.0, 0.0, s], [4.0, 0.0, 2.0], [4.0, 0.75, 2.0]],
                ],
            ),
        ] {
            // The boolean takes `&mut Model` and the audit `&Model`; a rejected boolean leaves
            // arena residue, so each gets its own build of the same fixture.
            let build = || {
                let mut m = Model::new();
                let plate = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([4.0, 4.0, 2.0]),
                );
                let boss = m.add_cylinder(
                    Point3::from_array(origin),
                    Vector3::from_array(axis),
                    0.5,
                    1.0,
                );
                m.rebuild_adjacency();
                (m, plate, boss)
            };
            let (mut m, plate, boss) = build();
            crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the arc class assembles");
            let (m, plate, boss) = build();
            let audits = crate::arrangement::frame_audit(&m, BoolKind::Fuse, plate, boss).unwrap();
            // The agreement, now that the population is green: the boolean built, and the audit's
            // replay of the same pipeline stops nowhere either.
            let stopped: Vec<_> = audits.iter().filter(|a| a.failed_at.is_some()).collect();
            assert!(stopped.is_empty(), "no class stops: {stopped:?}");
            // ★ The whole audit, not just `produced`: the ring coordinates are a sibling field, and
            // joining two filtered lists by index is the seam where they could come from different
            // classes. ★ By reference: `Produced` carries the labels now, so it is no longer `Copy`.
            let cut: Vec<_> = audits
                .iter()
                .filter(|a| a.produced.as_ref().is_some_and(|p| p.arcs > 0))
                .collect();
            assert_eq!(cut.len(), 1, "exactly one class has its circle cut");

            assert_eq!(
                *cut[0].produced.as_ref().unwrap(),
                crate::arrangement::Produced {
                    cells: 4,
                    arcs: 2,
                    roots: 1,
                    holes: 0,
                    // ★★ **Derived from the geometry, not read back.** `[A_above, A_below,
                    // B_above, B_below]` for the three `+1` regions, sorted: the boss stands *on*
                    // the plate, so plate-minus-disk has the plate below it and nothing above;
                    // disk-minus-plate — the overhang's underside — has the boss above it and
                    // nothing below; and their intersection has both. The outside cell winds `-1`
                    // and is not here.
                    //
                    // ★ **The same three come out for the turned boss, and the story is the same
                    // one in its own frame** — this assertion is inside the two-fixture loop, so
                    // reading it as a sentence about the straddling boss alone would be reading
                    // half of what it checks. There the class is `x = 4`: the plate lies on the
                    // `−x` side, the boss's base cap outside it on `+x`, and the overlap has both.
                    pos_labels: vec![
                        [false, false, true, false],
                        [false, true, false, false],
                        [false, true, true, false],
                    ],
                },
                "the arrangement walked the arcs, nested them, and labelled them"
            );
            // ★★★★★ **And the faces it emitted, as coordinates.** This is what a wrong `sense` on a
            // split segment moves and nothing before it does: the cell counts, the nesting and the
            // sorted labels above are all identical whichever way the rings run.
            //
            // ★★★ **The premise first.** `emit_faces` emits CCW about `n_out(wc)`, and which face
            // is the class root — hence which way `n_out` points — is a plane-table fact, not one
            // of the fixture's numbers. Pinning it here means a root flip fails *this* assertion,
            // with its own sentence, instead of silently reversing every ring below. Both come
            // from the plane table, which no arrangement stage can move.
            let got_n_out: Vec<f64> = cut[0]
                .root_normal
                .iter()
                .map(|c| c * cut[0].orient_sign as f64)
                .collect();
            assert!(
                got_n_out
                    .iter()
                    .zip(&n_out)
                    .all(|(a, b)| (a - b).abs() < 1e-9),
                "the rings below are derived CCW about {n_out:?}, but the class faces {got_n_out:?}"
            );
            let got = cut[0]
                .outer_rings
                .as_ref()
                .expect("the per-class road emits the rings");
            assert_eq!(
                got.len(),
                want_rings.len(),
                "one outer ring per emitted face"
            );
            for (g, w) in got.iter().zip(&want_rings) {
                assert_eq!(g.len(), w.len(), "ring length: got {g:?}, want {w:?}");
                assert!(
                    g.iter()
                        .zip(w)
                        .all(|(p, q)| (0..3).all(|i| (p[i] - q[i]).abs() < 1e-9)),
                    "ring: got {g:?}, want {w:?}"
                );
            }
        }
    }

    /// **The seam realizes a pierce vertex and measures it — seen through the door production
    /// uses, because nothing else can see it at all.**
    ///
    /// ★★★ The deferred stopper intercepts the whole seam stretch, so with the pierce arm
    /// deleted the seam fails `PierceVertexUnnamed`, the interception swallows it, and **every
    /// boolean-level fence stays green** (measured — that probe is what forced this test). The
    /// arm's only witness is a direct second consumer of `seam_table` on the very faces
    /// production feeds it, so this walks production's stretch step for step: trace, clean,
    /// append the bands, build the seam.
    ///
    /// ★ The tolerance is asserted as a **bound**, never a copied value; the pierce coordinates
    /// themselves are pinned by `ClassAudit::outer_rings` through the same `pierce_point` road,
    /// so re-asserting them here would be a second copy of an existing lock — and a `Lo`/`Hi`
    /// mix-up cannot hide behind the bound either, because both crossings lie on every defining
    /// surface and `outer_rings` is what tells them apart.
    #[test]
    fn the_seam_realizes_a_pierce_vertex_and_measures_it() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
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
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, deferred) =
                crate::arrangement::trace_result_faces_full_for_test(
                    &m,
                    BoolKind::Fuse,
                    plate,
                    boss,
                    &jd,
                    faces_tab,
                    plane_ix,
                    cyls,
                    *n_a,
                    class_owner,
                    &trace_in,
                )
                .expect("the arc population traces");
            // The stopper socket is empty since the population went green — nothing defers.
            assert!(deferred.is_none(), "{deferred:?}");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                    .expect("the lateral faces emit"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
                .expect("the seam realizes pierce nodes");
            let pierce: Vec<_> = seam
                .iter()
                .filter(|sv| crate::combinatorics::pierce_name(sv.triple).is_some())
                .collect();
            assert_eq!(pierce.len(), 2, "both crossings reach the seam, once each");
            for sv in pierce {
                assert!(
                    sv.tol < 1e-12,
                    "a pierce realization sits on everything that defines it: tol {}",
                    sv.tol
                );
            }
        }
    }

    /// **Every ring node of the arc population's result faces earns a definition — measured
    /// through the door production uses.**
    ///
    /// ★★★ Same instrument shape as the seam fence, same reason: the deferred stopper stands
    /// behind `name_result_vertices` (at the assembly's very end now) and intercepts
    /// everything, so no reject name can testify the naming completed — a direct second consumer
    /// is the only witness. This walks production's road (trace → clean → bands → seam → naming)
    /// and asserts on its product.
    ///
    /// ★★ Disabling the pierce arm reddens both fixtures (and the walls-fallback cannot fake a
    /// pierce def past its arc-carrier guard). ★ The fallback itself went **zero-population** when
    /// the split-twin subdivision landed — the bitten corner's twins match now, so its def comes
    /// down the far-plane road and no probe reddens on the fallback alone; the subdivision has
    /// its own fence (`a_subdivided_twin_matches_its_neighbour_edge_for_edge`).
    #[test]
    fn every_result_vertex_of_the_arc_population_is_named() {
        for (origin, axis, bites_corner) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], false),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], true),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
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
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                    .expect("the lateral faces emit"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming stages run");
            let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
            // Completeness: every ring node has a definition.
            let mut missing = Vec::new();
            for (fi, lf) in live.iter().enumerate() {
                for ring in lf.poly_rings() {
                    for &n in ring.iter() {
                        if !named.defs.contains_key(&(named.group_of[fi], n)) {
                            missing.push(n);
                        }
                    }
                }
            }
            assert!(missing.is_empty(), "def-less nodes: {missing:?}");
            let pierce_defs = named
                .defs
                .values()
                .filter(|d| matches!(d, crate::boolean::Def::Pierce { .. }))
                .count();
            assert_eq!(pierce_defs, 2, "both crossings are declared, once each");
            // The bitten corner's def names the right point: realize its three planes and land
            // on (4, 0, 2) — the fixture's own number, no class index copied.
            if bites_corner {
                let hit = named.defs.values().any(|d| {
                    let crate::boolean::Def::Three(t) = d else {
                        return false;
                    };
                    let t = t.planes();
                    nacre_geom::intersect::three_planes(
                        &geom[t[0]].plane,
                        &geom[t[1]].plane,
                        &geom[t[2]].plane,
                    )
                    .is_some_and(|p| {
                        let c = p.as_array();
                        (0..3).all(|i| (c[i] - [4.0, 0.0, 2.0][i]).abs() < 1e-9)
                    })
                });
                assert!(hit, "the bitten corner's def realizes to (4, 0, 2)");
            }
        }
    }

    /// **After the split-twin subdivision, every segment edge has exactly one twin.**
    ///
    /// ★★★ The subdivision cannot be locked through the naming fence — the walls-fallback
    /// rescues the bitten corner whether or not the twins match, so that fence is green either
    /// way. What the subdivision actually changes is the **edge-key census**: a neighbour's whole
    /// edge and the arc class's subdivided pieces share `norm_edge` keys only once the whole edge
    /// is cut at the same pierce nodes. So the proposition, with its one exception stated:
    ///
    /// > every segment ring edge (a `Wall::Plane` carrier) has its key used by exactly two
    /// > faces.
    ///
    /// Arc edges are excluded because their far side is the band, which contributes no ring.
    /// ★ Both-ends-pierce keys used to be excluded too — the chord and the two arcs between one
    /// pierce pair folded into a single `norm_edge` key — but the carrier gave arcs their own
    /// ordered key, so the chord's line key counts exactly its two coplanar faces now (measured:
    /// the exclusion removed, both fixtures stay green — the tightening the carrier cell's plan
    /// predicted).
    ///
    /// red: with the subdivision disabled, the pre-split single-use keys come back.
    #[test]
    fn a_subdivided_twin_matches_its_neighbour_edge_for_edge() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
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
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                    .expect("the lateral faces emit"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming stages run");
            let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
            let mut uses: std::collections::HashMap<
                (crate::combinatorics::NodeId, crate::combinatorics::NodeId),
                usize,
            > = std::collections::HashMap::new();
            let mut plain: Vec<(crate::combinatorics::NodeId, crate::combinatorics::NodeId)> =
                Vec::new();
            for lf in live {
                for ring in lf.poly_rings() {
                    let k = ring.nodes.len();
                    for t in 0..k {
                        let (a, b) = (ring.nodes[t], ring.nodes[(t + 1) % k]);
                        if matches!(ring.walls[t], crate::boolean::Wall::Arc { .. }) {
                            continue;
                        }
                        let key = crate::boolean::norm_edge(a, b);
                        *uses.entry(key).or_insert(0) += 1;
                        if !plain.contains(&key) {
                            plain.push(key);
                        }
                    }
                }
            }
            let odd: Vec<_> = plain
                .iter()
                .filter(|k| uses[*k] != 2)
                .map(|k| (*k, uses[k]))
                .collect();
            assert!(odd.is_empty(), "edges without exactly one twin: {odd:?}");
        }
    }

    /// **The pierce vertices are minted — canonical, measured, on the derived crossings.**
    ///
    /// Before the boolean no `Vertex::Pierce` exists anywhere in the model, so a whole-store
    /// filter is position-independent; a wrong `QuadRoot` canonicalization moves the minted
    /// point itself, which is what keeps the old toggle-lock alive now that the reject (whose
    /// witness once carried it) is gone.
    ///
    /// ★ The coordinates are the fixtures' own crossing derivations (the same numbers the ring
    /// and seam fences pin) — nothing here is copied from a run. The tolerance is a bound, and
    /// ascending handle order is `Vertex::Pierce`'s own contract, minted through
    /// `QuadRoot::canonical`'s second answer.
    #[test]
    fn a_pierce_vertex_is_minted_and_measured() {
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        for (origin, axis, crossings) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
            ),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let pierce: Vec<_> = (0..m.vertex_count() as u32)
                .filter_map(|i| m.vertex_handle_at(i))
                .map(|h| (h, m.vertex(h)))
                .filter(|(_, v)| matches!(**v, nacre_topo::Vertex::Pierce { .. }))
                .collect();
            assert_eq!(pierce.len(), 2, "both crossings minted, once each");
            for (h, v) in pierce {
                let nacre_topo::Vertex::Pierce { planes: [a, b], .. } = *v else {
                    unreachable!("filtered above");
                };
                assert!(a < b, "planes in ascending handle order: {a:?} vs {b:?}");
                let p = m.vertex_point(h);
                assert!(
                    crossings
                        .iter()
                        .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                    "a minted pierce vertex sits on a derived crossing: {p:?}"
                );
                knowledge_is_tight(&m, h);
            }
        }
    }

    /// **The two complementary arcs are two edges, and the cut rim is none.**
    ///
    /// ★★★ The observable form of the `[A, B]`-CCW convention (`derive_edge_curve`'s circle
    /// arm): between one pair of pierce vertices a circle offers two pieces, the endpoints
    /// alone cannot tell them apart, and the *vertex order* is the bit that does — so the store
    /// must hold **two** circle-carrier edges whose vertex pairs are each other's reverse.
    /// Erase the order from the welding key and they fold into one edge (the red probe this
    /// fence was built against).
    ///
    /// ★★ The rim skip's two sides, on one store: the **cut** circle mints no closed `[v, v]`
    /// edge (before the skip, both fixtures minted one that only the reject discarded), while
    /// the **uncut** far rim still mints exactly one — the skip's negative control, pinned to
    /// the far cap's axis coordinate so a skip that turned into "skip every rim" reddens here.
    ///
    /// ★ The two populations differ on the chord, deliberately: the straddling boss's pierce
    /// pair is joined by the plate-top chord (welded with the subdivided middle piece into
    /// **one** line edge used by both coplanar faces — the subdivision cell's promise realized
    /// in the store), while the turned boss's pair sits across the plate corner, joined through
    /// it by split boundary edges — no chord at all. Scoped to the edges the boolean minted
    /// (a snapshot, not a whole-store filter: the input cylinder's own rims are `[v, v]` too).
    #[test]
    fn an_arc_and_its_complement_are_minted_as_two_ordered_edges() {
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        // `seam_split`: where θ = 0 sits. `None` = a pierce vertex lies on the seam generator
        // (the straddling boss — the split's own `SeamIncident` case), so no piece splits;
        // `Some(p)` = the seam vertex S is minted at `p` (centre + ref_dir·r, derived) and the
        // wrap arc is cut there into two pieces.
        for (origin, axis, height, crossings, chord_edges, seam_split) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                1.0,
                [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
                1usize,
                None,
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                1.0,
                [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
                0usize,
                Some([4.0, 0.25, 1.5]),
            ),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                height,
            );
            m.rebuild_adjacency();
            let minted_from = m.edge_count();
            let vertices_from = m.vertex_count();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            let minted: Vec<_> = (minted_from as u32..m.edge_count() as u32)
                .filter_map(|i| m.edge_handle_at(i))
                .map(|h| (h, m.edge(h)))
                .collect();
            let is_cyl = |sh| matches!(m.surface_cache(sh), nacre_geom::Surface::Cylinder(_));
            // A circle-carrier edge is a **mixed** pair (the cap plane and the lateral); the
            // band's seam edge is `[lat, lat]` — both carriers the cylinder — and is not a
            // piece of any circle.
            let on_circle = |e: &nacre_topo::Edge| is_cyl(e.surfaces[0]) != is_cyl(e.surfaces[1]);

            // The cut circle's pieces: their directed vertex pairs chain into **one** cycle —
            // the observable of the `[A, B]`-CCW convention (for two pieces the cycle *is* the
            // mutually-reversed pair the first version of this fence asserted), and of the seam
            // split (three pieces when S stands, still one circle).
            let arcs: Vec<_> = minted
                .iter()
                .filter(|(_, e)| on_circle(e) && e.vertices[0] != e.vertices[1])
                .collect();
            let expect_pieces = 2 + usize::from(seam_split.is_some());
            assert_eq!(
                arcs.len(),
                expect_pieces,
                "the cut circle's pieces: {arcs:?}"
            );
            let mut succ: std::collections::HashMap<_, _> = std::collections::HashMap::new();
            for (_, e) in &arcs {
                let prev = succ.insert(e.vertices[0], e.vertices[1]);
                assert!(prev.is_none(), "two pieces leave one vertex CCW: {arcs:?}");
            }
            let start = arcs[0].1.vertices[0];
            let (mut cur, mut steps) = (start, 0usize);
            loop {
                cur = succ[&cur];
                steps += 1;
                if cur == start || steps > expect_pieces {
                    break;
                }
            }
            assert_eq!(steps, expect_pieces, "the pieces close one circle");
            // Endpoints: exactly two pierce vertices on the derived crossings, plus — when the
            // seam splits an arc — one OnSeam vertex at the derived seam point, with its
            // tolerance measured.
            let (mut pierce, mut on_seam) = (Vec::new(), Vec::new());
            for &v in succ.keys() {
                match *m.vertex(v) {
                    nacre_topo::Vertex::Pierce { .. } => pierce.push(v),
                    nacre_topo::Vertex::OnSeam(_) => on_seam.push(v),
                    ref d => panic!("an arc endpoint is neither pierce nor seam: {d:?}"),
                }
            }
            assert_eq!(pierce.len(), 2, "one pierce pair");
            for &v in &pierce {
                let p = m.vertex_point(v);
                assert!(
                    crossings
                        .iter()
                        .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                    "an arc endpoint sits on a derived crossing: {p:?}"
                );
            }
            match seam_split {
                None => assert!(on_seam.is_empty(), "the seam is the pierce vertex itself"),
                Some(sp) => {
                    assert_eq!(on_seam.len(), 1, "one seam vertex on the cut circle");
                    let p = m.vertex_point(on_seam[0]);
                    assert!(
                        (0..3).all(|i| (p.as_array()[i] - sp[i]).abs() < 1e-12),
                        "S sits on the derived seam point: {p:?}"
                    );
                    knowledge_is_tight(&m, on_seam[0]);
                }
            }
            // The minted OnSeam census: the uncut far rim's vertex, plus S when it stands —
            // and nothing else (a duplicate S at a seam-incident pierce vertex would show here).
            let minted_on_seam = (0..m.vertex_count() as u32)
                .filter_map(|i| m.vertex_handle_at(i))
                .map(|h| (h, m.vertex(h)))
                .skip(vertices_from)
                .filter(|(_, v)| matches!(**v, nacre_topo::Vertex::OnSeam(_)))
                .count();
            assert_eq!(
                minted_on_seam,
                1 + on_seam.len(),
                "far rim + S, nothing else"
            );

            // The rim skip's two sides: no closed edge on the cut circle, exactly one on the
            // uncut far rim (its axis coordinate is the far cap's, derived from the fixture).
            let closed: Vec<_> = minted
                .iter()
                .filter(|(_, e)| e.vertices[0] == e.vertices[1])
                .collect();
            assert_eq!(closed.len(), 1, "one uncut rim, no cut one: {closed:?}");
            let far = (0..3)
                .map(|i| (origin[i] + axis[i] * height) * axis[i])
                .sum::<f64>();
            let p = m.vertex_point(closed[0].1.vertices[0]).as_array();
            let along = (0..3).map(|i| p[i] * axis[i]).sum::<f64>();
            assert!(
                (along - far).abs() < 1e-12,
                "the closed rim is the far cap's: {along} vs {far}"
            );

            // The chord: welded into one line edge on the straddling boss, absent across the
            // turned boss's corner.
            let (ba, bb) = (pierce[0], pierce[1]);
            let chords = minted
                .iter()
                .filter(|(_, e)| {
                    !on_circle(e) && (e.vertices == [ba, bb] || e.vertices == [bb, ba])
                })
                .count();
            assert_eq!(chords, chord_edges, "the pierce pair's line edges");
        }
    }

    /// **The band's loop is one continuous cycle, and the shell closes over it.**
    ///
    /// ★★★ The band assembles its cut rim from the arc chain; this fence counts the closure
    /// directly on the store beside the production guard: every edge the boolean minted is used
    /// exactly twice across its faces, and the band face's outer loop is one vertex-continuous
    /// cycle of the derived length, with the seam edge traversed once in each sense.
    ///
    /// ★ Three fixtures: the straddling boss (lo rim cut, seam ≡ pierce), the turned boss
    /// (lo rim cut, seam splits the wrap arc — six half-edges), and the **hung** boss (the
    /// straddling boss mirrored under the plate — the cut circle is the band's **hi** end, so
    /// the chain is walked reversed; measured to pass the gate before this fence was written).
    /// The hi-cut *and* seam-split combination has no fixture yet — the chain logic is shared,
    /// and its population brings one when it arrives.
    #[test]
    fn the_bands_loop_is_one_continuous_cycle() {
        for (origin, axis, band_len) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], 5usize),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], 6usize),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0], 5usize),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let faces_from = m.face_count();
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let garbage: Vec<_> = (faces_from as u32..m.face_count() as u32)
                .filter_map(|i| m.face_handle_at(i))
                .map(|h| (h, m.face(h)))
                .collect();
            assert!(!garbage.is_empty(), "the face loop ran to completion");

            // Closure: every edge of the garbage faces is used exactly twice.
            let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
            for (_, f) in &garbage {
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for he in &lp.half_edges {
                        *uses.entry(he.edge).or_default() += 1;
                    }
                }
            }
            let odd: Vec<_> = uses.iter().filter(|&(_, &n)| n != 2).collect();
            assert!(
                odd.is_empty(),
                "a closed shell uses every edge twice: {odd:?}"
            );

            // The band face: one, on the cylinder, its outer loop a continuous cycle of the
            // derived length, the seam edge once in each sense.
            let bands: Vec<_> = garbage
                .iter()
                .filter(|(_, f)| {
                    matches!(m.surface_cache(f.surface), nacre_geom::Surface::Cylinder(_))
                })
                .collect();
            assert_eq!(bands.len(), 1, "one band face");
            let lp = &bands[0].1.outer;
            assert_eq!(lp.half_edges.len(), band_len, "the derived loop length");
            let ends = |he: &nacre_topo::HalfEdge| {
                let [a, b] = m.edge(he.edge).vertices;
                if he.forward { (a, b) } else { (b, a) }
            };
            for w in 0..lp.half_edges.len() {
                let (_, e0) = ends(&lp.half_edges[w]);
                let (s1, _) = ends(&lp.half_edges[(w + 1) % lp.half_edges.len()]);
                assert_eq!(e0, s1, "the loop chains vertex to vertex at step {w}");
            }
            let mut seen: std::collections::HashMap<_, Vec<bool>> =
                std::collections::HashMap::new();
            for he in &lp.half_edges {
                seen.entry(he.edge).or_default().push(he.forward);
            }
            let twice: Vec<_> = seen.values().filter(|v| v.len() == 2).collect();
            assert_eq!(twice.len(), 1, "exactly one edge is walked twice: the seam");
            assert_ne!(twice[0][0], twice[0][1], "once in each sense");
        }
    }

    /// **The grouping joins across a cut rim — through the door production uses.**
    ///
    /// ★★★ A cut circle bounds no whole disk, so the rim-key rule cannot see it, and the
    /// unordered node rule cannot either: on a 2-node circle the chord and both complementary
    /// arcs fold into one `norm_edge` pair (six users where each piece has two). The `JoinKey`'s
    /// ordered arc pairs — "a line is unordered, a circle is ordered", third appearance — give
    /// every piece exactly its cap and the band, so the result comes back as **one** component.
    ///
    /// ★ Asserted on `name_result_vertices`' own product (the subdivision must run first for the
    /// chord's line key to match), across all three fixtures. red: the band's registration
    /// removed → n == 2 and the held grouping is the old `PierceVertexUnnamed`.
    #[test]
    fn the_grouping_joins_across_a_cut_rim() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
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
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                crate::cyl_chart::emit_lateral(BoolKind::Fuse, &jd, cyls, &faces, &curved, &rows)
                    .expect("the lateral faces emit"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
                .expect("the seam realizes pierce nodes");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming runs");
            let g = named
                .grouping
                .as_ref()
                .expect("the grouping joins across the cut rim");
            assert_eq!(g.n, 1, "one component");
            assert_eq!(g.positives, vec![0], "one material piece, no cavity");
        }
    }

    /// **A cut-rim boolean builds a complete solid, and the integrals pin it exactly.**
    ///
    /// ★★★ The volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap — every fixture
    /// is a contact), so a wrong segment sign or scale moves an exact number; and because all
    /// three fixtures share that one number, the straddling boss's two z = 2 faces additionally
    /// split the segment **signs** (the plate top loses its half-disk bite, the digon *is* the
    /// bite). The structure is pinned beside it: one solid, no cavities, its outer shell exactly
    /// the boolean's minted faces.
    #[test]
    fn a_cut_rim_boolean_builds_a_complete_solid() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let (faces_from, solids_from) = (m.face_count(), m.solid_count());
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let pushed: Vec<_> = (solids_from as u32..m.solid_count() as u32)
                .filter_map(|i| m.solid_handle_at(i))
                .map(|h| (h, m.solid(h)))
                .collect();
            let [(sh, solid)] = pushed[..] else {
                panic!("one result solid, got {}", pushed.len());
            };
            assert_eq!(sh, out[0], "the pushed solid is the returned one");
            assert!(solid.cavities.is_empty(), "one material piece, no cavity");
            let shell_faces: std::collections::HashSet<_> =
                m.shell(solid.outer).faces.iter().copied().collect();
            let minted: std::collections::HashSet<_> = (faces_from as u32..m.face_count() as u32)
                .filter_map(|i| m.face_handle_at(i))
                .collect();
            assert_eq!(
                shell_faces, minted,
                "the outer shell is exactly the boolean's minted faces"
            );
            // ★★ **The arc integrals answer exactly** —
            // the volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap: every
            // fixture is a contact), so a wrong segment sign or scale moves an exact number.
            let props = nacre_props::mass_props(&m, sh).expect("the mixed loops integrate");
            let expect = 32.0 + std::f64::consts::PI / 4.0;
            assert!(
                (props.volume - expect).abs() < 1e-12,
                "volume {} vs derived {expect}",
                props.volume
            );
            // ★ All three fixtures share that one number, so a global scale error fools them
            // together — the straddling boss's two z = 2 faces split the segment **signs**: the
            // plate top loses its half-disk bite (16 − π/8), the overhang digon *is* the other
            // half-disk (π/8).
            if origin == [4.0, 2.0, 2.0] {
                let (mut top, mut digon) = (None, None);
                let mut i = faces_from as u32;
                while let Some(h) = m.face_handle_at(i) {
                    i += 1;
                    let f = m.face(h);
                    let nacre_geom::Surface::Plane(p) = m.surface_cache(f.surface) else {
                        continue;
                    };
                    if (p.normal().as_array()[2] - 1.0).abs() > 1e-9 {
                        continue;
                    }
                    match f.outer.half_edges.len() {
                        2 => digon = Some(h),
                        6 => top = Some(h),
                        1 => {} // the boss's closed top cap
                        n => panic!("an unexpected z-normal face with {n} half-edges"),
                    }
                }
                let area = |h| nacre_props::face_props(&m, h).expect("planar").area;
                let bite = std::f64::consts::PI / 8.0;
                assert!(
                    (area(top.expect("plate top")) - (16.0 - bite)).abs() < 1e-12,
                    "the plate top loses its bite"
                );
                assert!(
                    (area(digon.expect("overhang digon")) - bite).abs() < 1e-12,
                    "the digon is the bite"
                );
            }
        }
    }

    #[test]
    #[ignore = "scratch OBJ dump for eyeballing — run on demand"]
    fn dump_straddling_boss_obj() {
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
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("builds");
        m.rebuild_adjacency();
        let obj = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("tessellates")
            .to_obj();
        // Written only on request — the `--ignored` sweep runs this test too, and a test that
        // writes outside the workspace on every sweep is a side effect nobody asked for.
        match std::env::var("OBJ_OUT") {
            Ok(path) => {
                std::fs::write(&path, obj).expect("write");
                println!("wrote {path}");
            }
            Err(_) => println!("set OBJ_OUT=<path> to write the OBJ ({} bytes)", obj.len()),
        }
    }

    /// The chained bore-then-boss body: `cut` a through-bore at `(12,12)`, then `fuse` a boss
    /// whose rim the plate's edge cuts. The chain is what the milestone fence could not ask:
    /// nesting must place the bore's rim circle inside a **bitten** top ring — two pierce
    /// corners and an arc step — where the chart road's parity had no rational corners to read
    /// (`WitnessNotRational`, chaining wall 3) and the mixed parity answers in ℚ(√c).
    fn a_bored_plate_with_a_boss(boss_base: [f64; 3]) -> f64 {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([12.0, 12.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            4.0,
            20.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        let boss = m.add_cylinder(
            Point3::from_array(boss_base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            10.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, bored, boss).expect("the boss fuses");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the arcs tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "the mesh is watertight");
        nacre_props::mass_props(&m, out[0]).expect("props").volume
    }

    /// ★ **A bored plate takes a straddling boss — the first two-boolean chain through the arc
    /// population is green.** Volume `40·40·20 − π·4²·20 + π·5²·10 = 32000 − 70π` (the boss
    /// stands wholly above the plate top, so fuse adds its full cylinder).
    #[test]
    fn a_bored_plate_takes_a_straddling_boss() {
        let v = a_bored_plate_with_a_boss([40.0, 20.0, 20.0]);
        let want = 32000.0 - 70.0 * std::f64::consts::PI;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The mirror: the boss hangs under the plate's bottom edge — same chain, same volume, the
    /// cut circle at the band's hi end instead of lo.
    #[test]
    fn a_bored_plate_hangs_a_straddling_boss() {
        let v = a_bored_plate_with_a_boss([40.0, 20.0, -10.0]);
        let want = 32000.0 - 70.0 * std::f64::consts::PI;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **The chained contact-cut builds clean** (grouping-arm cell). Cut the same chain
    /// instead of fusing: the boss only *touches* the bored plate's top, so the cut removes
    /// nothing — and the result is the bored plate, exactly. It used to refuse
    /// `VertexNamesAbsentSurface` here: the top ring's rim seam could not merge (the arc pair
    /// collided in the merge's node-pair key), the unmerged ring kept pierce vertices whose
    /// definitions name the boss's cylinder, and the result has no face on it. The two-pass
    /// erase removes the seam — and the vertices with it — so the honest refusal became the
    /// honest build. Volume and validate lock that the build is *right*, not merely green.
    #[test]
    fn a_chained_contact_cut_builds_clean() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([12.0, 12.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            4.0,
            20.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        let boss = m.add_cylinder(
            Point3::from_array([40.0, 20.0, 20.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            10.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, bored, boss).expect("the contact cut");
        assert_eq!(out.len(), 1, "one body");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32000.0 - std::f64::consts::PI * 16.0 * 20.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **Two cylinders with coplanar caps fuse apart.** They stand `5` apart with their caps in
    /// the same two planes, which is what once made them look like a seating problem; what was
    /// actually hard was classifying **two curved bodies**, and both are now probed from a cap
    /// disk's centre — neither has a vertex to probe from.
    ///
    /// ★ Both volumes are the same number derived the same way (`π·0.5²·2`), so the assertion
    /// cannot pass by matching one body against the other.
    #[test]
    fn two_cylinders_with_coplanar_caps_fuse_apart() {
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
        let out = crate::boolean(&mut m, BoolKind::Fuse, a, b).expect("two curved bodies");
        assert_eq!(out.len(), 2, "they do not touch");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let want = std::f64::consts::PI * 0.25 * 2.0;
        for &s in &out {
            let v = nacre_props::mass_props(&m, s).expect("props").volume;
            assert!((v - want).abs() < 1e-9, "an untouched cylinder: {v}");
            let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
            assert_eq!(
                v_n - e_n + f_n - l_n,
                2,
                "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
            );
        }
    }

    /// `χ = V − E + F − L` over **one solid's** shells — the counts the genus relation
    /// `χ = 2(S − G)` is read from.
    ///
    /// ★ Per solid, not per model: `Model::reachable` spans every live solid, which coincides with
    /// this only while a fixture makes exactly one. The fixtures that make two need the sum split,
    /// and one spelling for all of them is what keeps the relation from being restated.
    fn euler_counts(m: &Model, s: Handle<Solid>) -> (i64, i64, i64, i64) {
        let faces: Vec<_> = crate::planes::solid_shell_handles(m, s)
            .into_iter()
            .flat_map(|sh| m.shell(sh).faces.clone())
            .collect();
        let mut verts = std::collections::HashSet::new();
        let mut edges = std::collections::HashSet::new();
        let mut loops = 0i64;
        for &fh in &faces {
            let f = m.face(fh);
            loops += f.inner.len() as i64;
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    edges.insert(he.edge);
                    verts.extend(m.edge(he.edge).vertices);
                }
            }
        }
        (
            verts.len() as i64,
            edges.len() as i64,
            faces.len() as i64,
            loops,
        )
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

    /// ★ **The d = 0 member of the family: the wall plane runs exactly through the bore's
    /// axis** (M6-2 rulings ladder). The gate passes it the same way (the boss's face clears the
    /// footprint along the span), and the rulings road stays **silent** — its trigger is the
    /// gate's carried `crossings` record, which is empty for every gate-passed pair. This
    /// geometry had no corpus row (its siblings are all d = 2), and it is the population an
    /// unconditionally-firing contribution arm broke — measured, 14 arc tests red — so it locks
    /// "an empty record changes nothing".
    #[test]
    fn a_boss_whose_wall_plane_holds_the_bores_axis_still_fuses() {
        let (mut m, holed, boss) = plate_bore_and_boss(12.0, 5.0);
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

    /// The d = 0 member with **caps in play**: a cylinder tool standing clear of the plate, its
    /// axis exactly on the plate's `x = 40` wall plane and its caps strictly inside the plate's
    /// height (no coplanar contact anywhere). Two disjoint bodies is the valid fuse answer, and
    /// the chord arm — like the ruling arm — stays silent behind the empty record.
    #[test]
    fn a_capped_tool_on_the_plates_wall_plane_fuses_as_two_bodies() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let tool = m.add_cylinder(
            Point3::from_array([40.0, 30.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            3.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, tool).expect("a distant fuse");
        assert_eq!(out.len(), 2, "two disjoint bodies");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [std::f64::consts::PI * 9.0 * 3.0, 4000.0];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "{vols:?} vs {want:?}"
        );
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

    /// **A wall face that really does cross the bore builds** (cell ③). The plate's own `y = 20`
    /// wall would clear, so the tool here is a slab whose face runs right across the hole — the
    /// wall rule's true population, and M6-2b's. It used to be the fence (`WallMeetsLateral`, a
    /// reason since retired); the gate records the pair now and the tracer cuts the bore's lateral
    /// along two rulings
    /// (2 from the axis, r = 3) and its caps along the chord: one body, the exact volume.
    #[test]
    fn a_wall_face_that_really_crosses_the_bore_builds() {
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
        let out =
            crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect("a real crossing builds");
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1);
        assert!(nacre_validate::validate(&m).is_empty());
        // The slab's box less the bore's `y ≥ 12` segment (d = 2, r = 3), taken from the bored plate.
        let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
        let pi = std::f64::consts::PI;
        let want = 4000.0 - 45.0 * pi - (1600.0 - 5.0 * seg(2.0, 3.0));
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **Exact tangency, and the operation is what refuses it.** The slab's face stands exactly
    /// `r` from the axis, so it touches the bore's lateral along one line. The gate passes that
    /// now (cell ⑥) — what convicts this shape is the *verdict*: cutting the slab away leaves the
    /// material as the **two wedges** between the parabola and the plane, which meet only on the
    /// line, and both are bounded by the same faces so they are one solid. `SelfTouchingResult`.
    ///
    /// ★ Its siblings are the control: the same wall with the boss *inside* the plate builds under
    /// `Fuse` and `Common` and only `Cut` pinches, and a boss tangent from **outside** comes back
    /// as two valid bodies. So this is not "the gate refuses tangencies"; it is one operation's
    /// answer about one shape.
    #[test]
    fn a_wall_face_tangent_to_the_bore_pinches_under_cut() {
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
                    reason: RejectReason::SelfTouchingResult,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★★★ **A blind stud tangent to the wall, and the three answers it must get.** A unit cube
    /// and a stud whose axis stands `0.3` from the origin with `r = 0.2`: the wall `x = 0.5` is
    /// **exactly** `r` away, which is what a round set of dimensions produces. The three operations
    /// differ, and they differ without any of them being named — `keep` is asked about the three
    /// regions beside the tangent line and the answers fall out.
    ///
    /// ★ This used to be called *the user's script*. It is not: the app's `cylinder({center})`
    /// anchors the cylinder at mid-height, so the user's stud runs `z ∈ [−1, 1]` **through** the
    /// cube (cell ⑦b caught the difference). That model is the next test; this one — base at
    /// `z = 0`, one cap pinched — is its blind neighbour, and the shape cell ⑥ and OCCT measured.
    ///
    /// Volumes **derived, not copied**: the stud's footprint is `x ∈ [0.1, 0.5] × y ∈ [−0.2, 0.2]`,
    /// inside the cube's, and only `z ∈ [0, 0.5]` overlaps — so the shared volume is `π r² · 0.5`.
    #[test]
    fn a_blind_stud_tangent_to_the_wall_gives_three_answers() {
        let pi = std::f64::consts::PI;
        let (whole, shared) = (pi * 0.04 * 2.0, pi * 0.04 * 0.5);
        let build = || {
            let mut m = Model::new();
            let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
            let stud = m.add_cylinder(
                Point3::from_array([0.3, 0.0, 0.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.2,
                2.0,
            );
            m.rebuild_adjacency();
            (m, cube, stud)
        };
        for (kind, want) in [
            (BoolKind::Fuse, 1.0 + whole - shared),
            (BoolKind::Common, shared),
        ] {
            let (mut m, a, b) = build();
            let out =
                crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), 1, "{kind:?}");
            assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
            let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
            assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
        }
        // Cutting the stud out leaves the material as the **two wedges** beside the tangent line,
        // which meet only on it — one solid whose surface touches itself.
        let (mut m, a, b) = build();
        let err = crate::boolean(&mut m, BoolKind::Cut, a, b).expect_err("the bore pinches");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::SelfTouchingResult,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★★★ **The user's script — `cuboid()` fused with `cylinder({r: 0.2, h: 2, center: [0.3,0,0]})`
    /// — and its three answers.** `center` is mid-height, so the stud runs `z ∈ [−1, 1]` through
    /// the cube and is tangent to the wall `x = 0.5` along the cube's whole height. Both caps are
    /// pinched (the stud's circle touches each square's edge at one point) — and since cell ⑦c
    /// both are bridged and drawn: this is the model the tessellator was opened for.
    ///
    /// Volumes derived: the stud is `π r² · 2`, the part inside the cube `π r² · 1`.
    #[test]
    fn the_users_through_stud_gives_three_answers() {
        let pi = std::f64::consts::PI;
        let (whole, shared) = (pi * 0.04 * 2.0, pi * 0.04 * 1.0);
        let build = || {
            let mut m = Model::new();
            let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
            let stud = m.add_cylinder(
                Point3::from_array([0.3, 0.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.2,
                2.0,
            );
            m.rebuild_adjacency();
            (m, cube, stud)
        };
        for (kind, want) in [
            (BoolKind::Fuse, 1.0 + whole - shared),
            (BoolKind::Common, shared),
        ] {
            let (mut m, a, b) = build();
            let out =
                crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), 1, "{kind:?}");
            assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
            let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
            assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
            if kind == BoolKind::Fuse {
                let faces: usize = std::iter::once(m.solid(out[0]).outer)
                    .map(|sh| m.shell(sh).faces.len())
                    .sum();
                assert_eq!(
                    faces, 10,
                    "four walls, two pinched caps, two bands, two stud caps"
                );
                crate::tests::mesh_covers_faces("the user's through stud", &m, &out);
            }
        }
        // Cutting the stud out leaves the two wedges beside the tangent line, meeting only on
        // it, the whole height of the cube — one solid whose surface touches itself.
        let (mut m, a, b) = build();
        let err = crate::boolean(&mut m, BoolKind::Cut, a, b).expect_err("the bore pinches");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::SelfTouchingResult,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★★★ **The user's cross studs — a stud along `z` and one along `y` through the cube,
    /// fused in either order.** After the first fuse the `z` stud's remaining lateral faces sit
    /// at `|z| ≥ 0.5` and the `y` stud's whole surface within `|z| ≤ 0.2`: the two cylinder
    /// classes share no face, and the arrangement builds the result — measured before the gate
    /// opened, with the pair let through by hand: 14 faces, volume `1 + 0.08π`, meshed. The
    /// gate used to refuse the pair on the distance between their *axes*, zero since they cross
    /// at the origin — a fact about two infinite surfaces, not about any face.
    ///
    /// Volumes: each stud is `0.08π`, half of it inside the cube; the second stud meets the
    /// first only inside the cube.
    #[test]
    fn the_users_cross_studs_build_in_either_order() {
        let pi = std::f64::consts::PI;
        let z_stud = |m: &mut Model| {
            m.add_cylinder(
                Point3::from_array([0.0, 0.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.2,
                2.0,
            )
        };
        let y_stud = |m: &mut Model| {
            m.add_cylinder(
                Point3::from_array([0.0, -1.0, 0.0]),
                Vector3::from_array([0.0, 1.0, 0.0]),
                0.2,
                2.0,
            )
        };
        for z_first in [true, false] {
            let build = || {
                let mut m = Model::new();
                let cube =
                    m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
                let (first, second) = if z_first {
                    (z_stud(&mut m), y_stud(&mut m))
                } else {
                    (y_stud(&mut m), z_stud(&mut m))
                };
                m.rebuild_adjacency();
                let studded = crate::boolean(&mut m, BoolKind::Fuse, cube, first)
                    .expect("the first stud fuses")[0];
                m.rebuild_adjacency();
                (m, studded, second)
            };
            for (kind, want, faces_want) in [
                (BoolKind::Fuse, 1.0 + 0.08 * pi, Some(14)),
                (BoolKind::Common, 0.04 * pi, Some(3)),
                (BoolKind::Cut, 1.0, None),
            ] {
                let (mut m, a, b) = build();
                let out = crate::boolean(&mut m, kind, a, b)
                    .unwrap_or_else(|e| panic!("{kind:?} (z first: {z_first}): {e:?}"));
                m.rebuild_adjacency();
                assert_eq!(out.len(), 1, "{kind:?} (z first: {z_first})");
                assert!(
                    nacre_validate::validate(&m).is_empty(),
                    "{kind:?}: {:?}",
                    nacre_validate::validate(&m)
                );
                let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
                assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
                if let Some(n) = faces_want {
                    let faces = m.shell(m.solid(out[0]).outer).faces.len();
                    assert_eq!(faces, n, "{kind:?} (z first: {z_first})");
                }
                crate::tests::mesh_covers_faces("the user's cross studs", &m, &out);
            }
        }
    }

    /// **The same two studs without the cube really cross**, and the gate still says so: their
    /// faces meet along a quartic curve (M6b). This is the negative control of the face-level
    /// clearance — letting perpendicular pairs through blindly makes this fixture come back as
    /// two separate, individually valid solids, which is wrong.
    #[test]
    fn crossing_studs_are_still_a_cylinder_pair_that_meets() {
        let mut m = Model::new();
        let z = m.add_cylinder(
            Point3::from_array([0.0, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        let y = m.add_cylinder(
            Point3::from_array([0.0, -1.0, 0.0]),
            Vector3::from_array([0.0, 1.0, 0.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Fuse, z, y).expect_err("the studs cross");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CylinderPairContact,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// **A cross above a stud**: a `y` cylinder at `z = 2.5` over a `z` cylinder ending at
    /// `z = 2`. Their axes cross (distance zero), so the surface rule refuses, but the upper
    /// one's reach along `z` is `[2.3, 2.7]` and the lower one's face spans `[0, 2]` — clear
    /// by the face rule: two solids that never touch, and the boolean says so.
    #[test]
    fn a_cross_above_a_stud_is_two_solids() {
        let mut m = Model::new();
        let low = m.add_cylinder(
            Point3::from_array([0.0; 3]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        let high = m.add_cylinder(
            Point3::from_array([0.0, -1.0, 2.5]),
            Vector3::from_array([0.0, 1.0, 0.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Fuse, low, high).unwrap_or_else(|e| panic!("{e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), 2, "two bodies that never touch");
        assert!(nacre_validate::validate(&m).is_empty());
        for s in &out {
            let v = nacre_props::mass_props(&m, *s).unwrap().volume;
            assert!((v - 0.08 * std::f64::consts::PI).abs() < 1e-9, "{v}");
        }
    }

    /// **An oblique cross above the stud** — the same, with the upper axis `(0, 1, 1)`. The
    /// face rule would clear it too (`lateral_reach`'s `d·m ≠ 0` arm), but the pair never
    /// reaches the pair loop: the upper cylinder's caps are planes oblique to the lower axis,
    /// and the plane–cylinder gate refuses those without asking whether they clear — the same
    /// proposition still spelled at surface level there. A lock on today's name, for the cell
    /// that opens that arm to flip.
    #[test]
    fn an_oblique_cross_above_the_stud_is_refused_by_its_caps() {
        let mut m = Model::new();
        let low = m.add_cylinder(
            Point3::from_array([0.0; 3]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        let high = m.add_cylinder(
            Point3::from_array([0.0, -1.0, 3.0]),
            Vector3::from_array([0.0, 1.0, 1.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Fuse, low, high).expect_err("oblique caps");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::ObliqueCylinderCut,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★★★ **The bridge pre-pass splits the wall's two straight edges under the through stud.**
    ///
    /// A stud through the cube (`center` anchoring, `z ∈ [−1, 1]`) pinches **both** caps: on
    /// each, the stud's circle touches the square's `x = 0.5` edge at one point, and that edge
    /// is shared with the `x = 0.5` wall. Before the pre-pass a straight edge carries exactly its
    /// two ends, so the touching sample — a vertex of the circle — is a vertex of neither the
    /// square ring nor the wall. The pre-pass inserts it into both edges' polylines, so the cap
    /// gains its twin and the wall gains the same vertex (no T-vertex, no crack).
    ///
    /// Read through the test-only `bridge_report`, because `tessellate` still refuses these caps
    /// (the bridge itself is the next commit) and discards everything on the way out.
    #[test]
    fn the_bridge_prepass_splits_the_walls_edges_under_both_caps() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let stud = m.add_cylinder(
            Point3::from_array([0.3, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Fuse, cube, stud).expect("the through stud fuses");
        assert_eq!(out.len(), 1);
        m.rebuild_adjacency();
        let (report, t) = nacre_tess::bridge_report(&m, &nacre_tess::TessConfig::default());
        assert!(
            report.declined.is_empty(),
            "declined: {:?}",
            report.declined
        );
        assert_eq!(report.splits.len(), 2, "splits: {:?}", report.splits);
        assert_ne!(
            report.splits[0].edge, report.splits[1].edge,
            "one split per cap"
        );
        for s in &report.splits {
            assert!(
                matches!(m.edge_curve(s.edge), nacre_geom::Curve::Line(_)),
                "the split edge is straight"
            );
            let poly = &t.by_edge[&s.edge];
            assert_eq!(poly.len(), 3, "two ends and the touching sample");
            assert_eq!(poly[s.at], s.vertex);
            assert!(
                matches!(
                    t.vertices.get(s.vertex).origin,
                    nacre_tess::TessOrigin::OnEdge { .. }
                ),
                "the inserted vertex is the circle's own sample, not a new one"
            );
        }
        // The blind stud (base anchoring) pinches one cap only.
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let stud = m.add_cylinder(
            Point3::from_array([0.3, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.2,
            2.0,
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Fuse, cube, stud).expect("the blind stud fuses");
        m.rebuild_adjacency();
        let (report, _) = nacre_tess::bridge_report(&m, &nacre_tess::TessConfig::default());
        assert!(report.declined.is_empty(), "{:?}", report.declined);
        assert_eq!(report.splits.len(), 1, "{:?}", report.splits);
    }
    /// ★★★★★ **A tangency from *outside* is two bodies, not a refusal** — the control that says
    /// the rule is not "reject every tangency", and the measurement that corrected it.
    ///
    /// The boss's axis stands at `x = −0.5` with `r = 0.5`, so it touches the plate's wall `x = 0`
    /// from the void side. `Fuse` keeps the lens (the boss) and the far side (the plate) but not
    /// the wedges between — two lumps meeting only on the line — and because they share no face
    /// they are **two valid solids**, which is what the kernel returns. `Cut` removes nothing and
    /// `Common` is empty.
    #[test]
    fn a_boss_tangent_from_outside_is_two_bodies() {
        let pi = std::f64::consts::PI;
        let build = || {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array([-0.5, 2.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            (m, plate, boss)
        };
        for (kind, bodies, want) in [
            (BoolKind::Fuse, 2, 32.0 + pi),
            (BoolKind::Cut, 1, 32.0),
            (BoolKind::Common, 0, 0.0),
        ] {
            let (mut m, a, b) = build();
            let out =
                crate::boolean(&mut m, kind, a, b).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
            m.rebuild_adjacency();
            assert_eq!(out.len(), bodies, "{kind:?}");
            assert!(nacre_validate::validate(&m).is_empty(), "{kind:?}");
            let v: f64 = out
                .iter()
                .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                .sum();
            assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
        }
    }

    /// ★★★★★ **The same two lumps, joined elsewhere — and this one *is* a pinch.** The block has a
    /// notch, the boss stands in it tangent to the notch's wall `x = 2` from the void side, and it
    /// overlaps the block in `y` (the notch is `0.8` wide, the boss `1.0`). So the lens and the far
    /// side are the same body, and the tangent line is where that body's surface meets itself.
    ///
    /// ★ This is why the verdict asks the **grouping** and not only the geometry: without it, the
    /// answer here would be the one the fixture above earns, and this solid shipped `Ok` with
    /// `validate` clean and a mesh — measured before the check existed.
    #[test]
    fn a_boss_tangent_in_a_notch_pinches_the_block_it_joins() {
        let mut m = Model::new();
        let outer = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([6.0, 6.0, 2.0]),
        );
        let notch = m.add_cuboid(
            Point3::from_array([2.0, 2.6, -0.5]),
            Point3::from_array([6.5, 3.4, 2.5]),
        );
        m.rebuild_adjacency();
        let block = crate::boolean(&mut m, BoolKind::Cut, outer, notch).expect("notch")[0];
        m.rebuild_adjacency();
        let boss = m.add_cylinder(
            Point3::from_array([2.5, 3.0, -0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            3.0,
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Fuse, block, boss).expect_err("a joined pinch");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::SelfTouchingResult,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// **The pinch formula, as a truth table** — four configurations × three operations × both
    /// owner orders, checked against the geometry by hand rather than against an engine run.
    ///
    /// Read the rows as: the cylinder's side of the wall plane is (or is not) the wall face's
    /// material side, and the cylinder keeps its material inside (a boss) or outside (a bore).
    /// ★ The third row's `Common` is the one that says this is not "a rule about `Cut`": a bore
    /// intersected with the wall's solid leaves the two wedges alone.
    #[test]
    fn the_pinch_formula_is_a_truth_table() {
        use crate::boolean::lumps_fall_apart;
        use crate::planes::SolidSide::{A, B};
        // (lens_in_wall_solid, cyl_orient, [Fuse, Cut, Common]) with the wall on side A.
        for (lens, orient, want) in [
            (true, 1i8, [false, true, false]), // an inside boss: only Cut pinches
            (false, 1, [true, false, false]),  // an outside boss: only Fuse splits into lumps
            (true, -1, [false, false, true]),  // a bore on the material side: only Common
            (false, -1, [false, false, false]), // a bore on the void side: never
        ] {
            for (i, kind) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
                .into_iter()
                .enumerate()
            {
                assert_eq!(
                    lumps_fall_apart(kind, A, lens, orient),
                    want[i],
                    "wall on A, lens {lens}, orient {orient}, {kind:?}"
                );
            }
        }
        // ★ Swapping the owners is not a symmetry: `Fuse` and `Common` are commutative but `Cut`
        // is not, so the same geometry judged with the wall on `B` reads `A − B` the other way.
        assert!(lumps_fall_apart(BoolKind::Cut, A, true, 1));
        assert!(!lumps_fall_apart(BoolKind::Cut, B, true, 1));
        assert!(lumps_fall_apart(BoolKind::Cut, B, false, -1));
        assert!(!lumps_fall_apart(BoolKind::Cut, A, false, -1));
    }

    /// ★★★★★ **The disk on a wall plane is read, and it clears** — cell ⑮.
    ///
    /// The plate carries a **crosswise** bore, so the plane `x = 30` holds that bore's circular cap
    /// face: an outer loop of one arc and one seam vertex. That plane is also parallel to the
    /// vertical drill's axis and passes within `r` of it, so the face test is reached. Until this
    /// cell the road answered «cannot read this shape» and the boolean stopped there — the barrier
    /// was sound (a rule reading vertices alone would have found "every vertex on one side" true
    /// of a **single point**) but it turned away the true answer along with the false one.
    ///
    /// ☑ The disk's own numbers: centre `(30, 10, 5)` radius `2`, against a strip centred on
    /// `y = 4` of half-width `√(9 − 1) ≈ 2.83`. Six apart, so it clears by more than its radius —
    /// and now says so. The volume is the oracle that the answer is not merely *an* answer.
    #[test]
    fn a_disk_face_on_a_wall_plane_is_read_and_clears() {
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
        let out = crate::boolean(&mut m, BoolKind::Cut, bored, drill).expect("the drill cuts");
        assert_eq!(out.len(), 1, "one body");
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        // 40 × 20 × 10, less the crosswise bore (r 2, thirty long inside) and the drill (r 3,
        // through the ten of thickness). Their axes pass six apart, clear of the radius sum.
        let want = 8000.0 - std::f64::consts::PI * (4.0 * 30.0 + 9.0 * 10.0);
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
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
        let mut counts: std::collections::HashMap<Handle<Surface>, usize> = Default::default();
        let solid = m.solid(s);
        for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
            for &fh in &m.shell(sh).faces {
                let surf = m.face(fh).surface;
                if matches!(m.surface(surf), nacre_topo::Surface::Cylinder { .. }) {
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
    // case bolted on — a row's span is the existence truth at an uncut end, so no band is ever
    // built in a gap and there is no premise there for a wall to break.

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
    /// and clearing neither is the **record** (cell ③ — it used to be the refusal). A rule that
    /// had merely gained a second, independent test would pass (a) and (b) too — what this pins
    /// is (c), that the two are OR-ed rather than each able to wave a face through on its own
    /// terms: a face waved through is an unlisted pair the tracer stays silent for, and the bore
    /// would be left uncut where the face crosses it — (c)'s exact volume is what says the pair
    /// was recorded.
    #[test]
    fn either_axis_clears_the_footprint_and_neither_does_not() {
        // (a) Across only: the `y = 6` face sits at `x ∈ [0, 2.5]`, clear of the strip
        //     `x ∈ [3.27, 6.73]`, while its plane still crosses the bore. Along the axis it runs
        //     `t ∈ [4,8]`, straddling both bands. ★ The `x = 2.5` wall must stand clear of the axis
        //     by more than `r`, or *it* becomes the face under test — at `x = 3` it is exactly
        //     tangent, and this arm would then measure the tangency instead of the clearance.
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

        // (c) Neither: full width *and* straddling both bands — a genuine crossing of both
        // laterals (`y = 6`, 1 from the axis, r = 2), which the gate used to refuse and records
        // now (cell ③): two ruling pieces per band at `x = 5 ± √3`, the chord on the caps and on
        // the middle cut's ceiling and floor. One body, the exact volume.
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let neither = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 3.0]),
            Point3::from_array([10.0, 6.0, 7.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, neither).expect("crosses both, builds");
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1);
        assert!(nacre_validate::validate(&m).is_empty());
        let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
        let pi = std::f64::consts::PI;
        // The plate less the bore, less the middle cut (less the bore inside it), less the tool's
        // box (less the bore's `y ≤ 6` part over its height, less the middle cut inside it — which
        // had already lost the bore's `y ≤ 6` part).
        let disk_le6 = 4.0 * pi - seg(1.0, 2.0);
        let want = 1000.0
            - 40.0 * pi
            - (72.0 - 8.0 * pi)
            - (240.0 - 4.0 * disk_le6 - (48.0 - 2.0 * disk_le6));
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
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

    /// **An enclosed cylindrical void is a cavity, not a second body.** A `[0,4]³` box `Cut` by a
    /// radius-`0.5`, height-`2` cylinder wholly inside it: the boolean's two components are the
    /// box's shell and the void's, and the void must come back **depth 1** — odd, so a hole.
    ///
    /// ★★★ **This is the population the coordinate probe exists for, and the only one that makes
    /// it say `true`.** The void's boundary is two disks and a band, so it carries no vertex at
    /// all and `nodes_of` came back empty — which is what `curved_component_depth`, now retired,
    /// used to refuse.
    /// Two *disjoint* bodies exercise the same road, but their answer is "outside" either way, so
    /// a road that always said `false` would pass them; here it would make the void a **second
    /// solid** and the box's volume whole.
    #[test]
    fn an_enclosed_cylindrical_void_is_a_cavity() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let void = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
        assert_eq!(out.len(), 1, "one body, hollow — not two");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(m.solid(out[0]).cavities.len(), 1, "the void is a cavity");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (v - want).abs() < 1e-9,
            "the box less the void: {v} vs {want}"
        );
        // Two shells, genus 0 each: `χ = V − E + F − L = 2(S − G) = 4`.
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            4,
            "two genus-0 shells: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }

    /// **A bore and a sealed void in one body.** `[0,8]×[0,4]×[0,4]` with a through bore at
    /// `(2,2)` and, `4` away, a wholly enclosed cylinder at `(6,2)` — `128 − π − π/2`.
    ///
    /// ★ It is the only fixture where the **coordinate** probe meets a face with a hole: the
    /// void has no vertex, so it is probed from a cap centre, and the body it asks about has
    /// annular caps. Two cylinders in one boolean also need `cylinders_clear`, which the `4`
    /// between the axes supplies.
    #[test]
    fn a_bore_and_a_sealed_void_live_in_one_body() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([8.0, 4.0, 4.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([2.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            6.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, a, bore).expect("the bore");
        let void = m.add_cylinder(
            Point3::from_array([6.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, bored[0], void).expect("the sealed void");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(m.solid(out[0]).cavities.len(), 1, "the void is a cavity");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 128.0 - std::f64::consts::PI * 0.25 * 4.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "the box less both: {v} vs {want}");
    }

    /// **A cavity finds its owner when there is more than one material.** The void fixture above
    /// takes a shortcut the code makes explicit — one material owns every cavity, no search — so
    /// the search itself is only reached when two materials stand apart. Here the hollow box is
    /// fused with a distant second box, and the void must still come back as the **hollow box's**
    /// cavity rather than a third body.
    ///
    /// ★ It is the only fixture that runs `first_deciding` over a **coordinate** probe in the
    /// cavity-owner search: the void has no vertex, and every other multi-material case in the
    /// suite has cavities with corners.
    #[test]
    fn a_cavity_with_no_vertex_still_finds_its_owner() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let void = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let hollow = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
        let b = m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([14.0, 14.0, 14.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, hollow[0], b).expect("two bodies apart");
        assert_eq!(out.len(), 2, "the void is a cavity, not a third body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut seen: Vec<(f64, usize)> = out
            .iter()
            .map(|&s| {
                (
                    nacre_props::mass_props(&m, s).expect("props").volume,
                    m.solid(s).cavities.len(),
                )
            })
            .collect();
        seen.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite"));
        let hollowed = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (seen[0].0 - hollowed).abs() < 1e-9 && seen[0].1 == 1,
            "the hollow box keeps its one cavity: {seen:?}"
        );
        assert!(
            (seen[1].0 - 64.0).abs() < 1e-9 && seen[1].1 == 0,
            "the second box is solid: {seen:?}"
        );
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

        // ★★ **The fuse this test's name promises — and the case that used to abort the kernel.**
        // Two disjoint bodies means `n = 2`, so each is classified by a ray probe. That probe
        // first panicked (it asked a band face for its plane class), then refused by name
        // (`curved_component_depth`, since retired), and now answers: the ray counts a
        // cylinder's crossings, and
        // the bare cylinder — whose boundary carries **no vertex at all** — is probed from a cap
        // disk's centre instead of from a corner.
        //
        // ★ The values are derived from the inputs, not read back from the result: a `[0,2]³` box
        // is `8`, and a radius-`0.5`, height-`1` cylinder is `π/4`. Two separate bodies, each a
        // sphere topologically (`χ = 2`), and `validate` clean.
        let mut m2 = Model::new();
        let a2 = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b2 = m2.add_cylinder(
            Point3::from_array([10.0, 10.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m2.rebuild_adjacency();
        let out = crate::boolean(&mut m2, BoolKind::Fuse, a2, b2).expect("two bodies apart");
        assert_eq!(out.len(), 2, "a fuse of two disjoint bodies is two bodies");
        assert!(
            nacre_validate::validate(&m2).is_empty(),
            "{:?}",
            nacre_validate::validate(&m2)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m2, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [std::f64::consts::PI * 0.25, 8.0];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "the cylinder and the box, each untouched: {vols:?}"
        );
        for &s in &out {
            let (v_n, e_n, f_n, l_n) = euler_counts(&m2, s);
            assert_eq!(
                v_n - e_n + f_n - l_n,
                2,
                "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
            );
        }
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
