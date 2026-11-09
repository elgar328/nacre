//! **The lateral face's row.**
//!
//! Which parts of a cylinder's lateral surface survive a boolean is
//! `cyl_chart::emit_lateral`'s decision, cell by cell, and the vocabulary it reads that decision
//! with lives beside it (`cyl_chart::read_cell`). What is here is the two facts the chart needs
//! *before* it can read anything: [`cyl_rows`]/[`CylRow`] (which lateral faces
//! there are, and the two boundary sections each ends at).
//!
//! The prose below is the theory [`cyl_rows`] rests on — the uniform-slab theorem, and the gate
//! that admits only the population it holds for.
//!
//! The plane arrangement decides one plane class at a time; a cylinder's wall is not a plane, so
//! it is decided here instead — and it is decided *coarsely*, because in this population it can
//! be. The gate admits only ⊥ cuts, ∥ walls whose **faces** stand clear of the lateral surface,
//! and oblique planes, and for ⊥ cuts and ∥ walls those facts give the **uniform-slab theorem**:
//!
//! > Between two consecutive ⊥ cuts, the other operand's boundary does not meet the open
//! > cylinder slab at all — the ⊥ faces are the cuts themselves, and every ∥ wall **face** misses
//! > the lateral. So membership in the other operand is **uniform** over that slab, and one witness
//! > decides the whole band.
//!
//! An oblique plane the gate passes touches no circle of a lateral's boundary, and where its faces
//! meet the lateral they meet it along the plane's section, which the chart reads as one more line
//! (`cyl_chart`'s module doc): the uniform ground is then the chart's cell between that line and
//! the cuts, and the cell is what one witness decides.
//!
//! **Faces, because that is what a boundary is made of.** Testing each ∥ wall's infinite *plane*
//! instead is a cheaper sufficient condition and reads almost the same — but it refuses bodies
//! standing well clear of the cylinder whose plane, extended, crosses it. The theorem asks about
//! the boundary, and so does the gate.
//!
//! **A wall face is asked about a rectangle, and the "consecutive cuts" above are per face.**
//! What a ∥ plane meets of the solid cylinder is a rectangle of that plane: the strip across, a
//! lateral face's axis-parameter span along. So the gate reads two separating axes, and a wall
//! clear on either one misses the rectangle (`planes`' `face_clears_footprint`). With several
//! lateral faces there are several rectangles, and that is exactly the granularity the theorem
//! needs: a row's ends are the existence truth at an uncut end and band boundaries
//! (`cyl_chart::boundary_lines`), so no band is ever built in the gap between two of them, and a
//! wall sitting in such a gap breaks no premise — there is none there to break. Reading the spans
//! as one `min..max` would invent both the band and the refusal.
//!
//! "The band is in the other solid iff its z-range is" reads plausibly and is **false**: an
//! L-notch's inner corner is clear of every wall by more than `r` while sitting outside the
//! material. Asking the arrangement is what answers it, and the counterexample is a fixture
//! (`a_cylinder_in_the_notch_keeps_no_band`, with its end-to-end twin: cutting with a drill
//! standing in the notch removes nothing).

#[cfg(test)]
use crate::BoolKind;
#[cfg(test)]
use crate::draft::{Bound, LocalFace};
use crate::planes::SolidSide;
use crate::planes::{ClassIx, FaceRow};
use crate::{BoolError, RejectReason, reject};
use nacre_exact::Rat;

/// **One lateral face**, as the band pass reads it.
///
/// **The row is a face, not a class.** A cylinder *class* is one lateral `Surface` handle
/// (`planes`' `cyl_ix.entry(cf.surf)`), and one surface can carry **several faces of one
/// solid** — bore a plate through and then cut away the middle of the bore, and its wall becomes
/// two disjoint bands. "Which parts survive" is a question about a face; the class only says which
/// surface that face lies on.
///
/// One row per class, with the first face speaking for the rest, would clip the others' spans away
/// on the band road and never emit their bands, leaving the result shell open — the boolean would
/// come back `OpenResultShell`, naming a symptom of the omission. The chart reads the same span as
/// the cell's existence at an uncut end (`Chart::read_cell`).
///
/// **Merging the spans into one `min..max` is not the fix either**: it would invent a band across
/// the gap where the solid has no face at all.
pub(crate) struct CylRow {
    /// The [`ClassIx::Cyl`] payload — which lateral surface this face lies on. Several rows may
    /// share it.
    pub(crate) class: usize,
    /// **This face's** own ends: the planes (world coefficients) of its lowest and highest boundary
    /// sections, low first — the lines a cell between which exists at an uncut end.
    pub(crate) ends: [[Rat; 4]; 2],
    /// How far this face reaches along the axis (its footprint's span) — the rows' order.
    pub(crate) reach: [Rat; 2],
    /// Which operand this face belongs to.
    pub(crate) side: SolidSide,
}

/// One row per lateral **face**, ordered by `(class, lowest reach)`.
///
/// That order is the replay contract stated at the row level, and the chart's emitter keeps it
/// (classes ascending, then lines ascending). A class's faces have disjoint reaches, so sorting by
/// `(class, reach[0])` keeps each class's bands in ascending `t`; with one face per class it *is*
/// class-index order.
///
/// **That disjointness holds where every lateral face that reaches here is a
/// band** — two whole rims and holes. The tracer states a panel and a chain rim too
/// (`arrangement::lateral_shape`), so two faces of one class *can* share a `t`; the tie then
/// falls to the stable sort's face order — still deterministic, and the corpus carries one
/// lateral per class (the cleaning pass merged the pieces one operation earlier). Named here so
/// the day it is exercised is a thing a reader can check, not a surprise.
///
/// A face whose ends cannot be stated declines by name (`DeclineKind::CylSpan`) — the tracer
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
        let (Some(ends), Some(reach)) = (cf.ends, cf.footprint.span) else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        // A lateral with no world statement declines the same way missing ends do — the chart
        // states every line against world planes.
        if cf.def.is_none() {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        }
        out.push(CylRow {
            class: k,
            ends,
            reach,
            side: if i < n_a { SolidSide::A } else { SolidSide::B },
        });
    }
    out.sort_by(|a, b| (a.class, a.reach[0]).cmp(&(b.class, b.reach[0])));
    // A class exists because a face made it, so a class with no row is the class table and the
    // face rows disagreeing rather than an input. It has never fired; an unfired guard stays,
    // under the name its whole class of guards shares (`CylinderStagesDisagree`).
    for k in 0..n_class {
        if !out.iter().any(|r| r.class == k) {
            return Err(reject(RejectReason::CylinderStagesDisagree));
        }
    }
    Ok(out)
}

#[cfg(test)]
#[path = "tests/bands/mod.rs"]
mod tests;
