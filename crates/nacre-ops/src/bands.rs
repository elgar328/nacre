//! **The lateral face's row, and the vocabulary the chart reads it with.**
//!
//! Which parts of a cylinder's lateral surface survive a boolean is
//! `cyl_chart::emit_lateral`'s decision, cell by cell. What
//! is here is what the chart reads, one spelling each: [`cyl_rows`]/[`CylRow`] (a face's span
//! is the existence truth at an uncut end), [`face_spans`] (existence at a cut end), [`read_bits`]
//! (which two bits of a label are a cell's chamber), [`keep_for`] (the wall is a boundary of its
//! own solid) and [`axis_param`]. The prose below is the theory those spellings state — the
//! uniform-slab theorem and why the chamber is *read* off the arrangement rather than measured.
//!
//! The plane arrangement decides one plane class at a time; a cylinder's wall is not a plane, so
//! it is decided here instead — and it is decided *coarsely*, because in this population it can
//! be. The gate admits only ⊥ cuts and ∥ walls whose **faces** stand clear of the lateral
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
//! clear on either one misses the rectangle (`planes`' `face_clears_footprint`). With several
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
//! ★ This replaced a witness ray (`point_in_faces_rat`): a rational point on the axis, cast
//! through the counterpart's faces, counting crossings. It gave the same answers, but it answered
//! a **3D containment** question — the shape the *component* probe asks — when the band's question
//! is the same shape as "does this face survive": two chambers either side of a boundary. Reading
//! the label needs no coordinates, no ray direction, no abstention retry, and has no width
//! ceiling, and it is why a cylinder may now stand on either side of the boolean.

use crate::combinatorics;
#[cfg(test)]
use crate::draft::{Bound, LocalFace};
use crate::planes::SolidSide;
use crate::planes::{ClassIx, FaceRow, WorkingPlane, axis_param_of_plane};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, reject};
use nacre_exact::Rat;

/// **One lateral face**, as the band pass reads it.
///
/// ★★ **The row is a face, not a class.** A cylinder *class* is one lateral `Surface` handle
/// (`planes`' `cyl_ix.entry(cf.surf)`), and one surface can carry **several faces of one
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
/// ★★ **That disjointness holds where every lateral face that reaches here is a
/// band** — two whole rims and holes. The tracer states a panel and a chain rim too
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
/// ☑ **How often each row actually fires** (whole binary, production calls only):
/// `Transversal` **265** · `Graze` **13**, of which **1** reaches and **12**
/// do not · nothing at all **0**. The chart reads it once per cut end (`cyl_chart::census`'s
/// `arcs_read`, 756 over the lib suite).
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
            // too, and never a rim arc's mark.
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
        SolidSide::A => crate::draft::keep(kind, in_own, in_other),
        SolidSide::B => crate::draft::keep(kind, in_other, in_own),
    }
}

#[cfg(test)]
#[path = "tests/bands/mod.rs"]
mod tests;
