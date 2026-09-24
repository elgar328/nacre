//! The b-rep side of the toleranced predicates: `nacre-ops`'s arrangement tables
//! ([`WorkingPlane`], [`crate::planes::FaceInfo`]) implement the [`Witness`]/[`PlaneWitness`]
//! ports so the
//! rotation-general sign predicates in [`nacre_judge::predicate`] can run over them.
//!
//! The tables are **pure description** — geometry and provenance, nothing about how this
//! operation judges. That belongs to [`Judge`], which the engine builds once per boolean
//! (`arrangement::plane_index_setup` supplies the standard and the collector) and passes down; the predicates
//! are its methods. The predicate logic itself (exact-vs-kernel routing, the frame3 judges) lives
//! in `nacre-judge`.

use crate::planes::{FaceRow, WorkingPlane};
use nacre_judge::WitnessPoint;
use nacre_judge::predicate::{PlaneWitness, Witness};
use nacre_math::Point3;

// The judging context and the helpers the engine reaches for by name. The predicates themselves
// are methods on `Judge`, so there is nothing else to re-export.
#[cfg(test)]
use nacre_judge::predicate::plane_def;
pub(crate) use nacre_judge::predicate::{ImplicitPoint, Judge};

impl Witness for WorkingPlane {
    fn base_coeffs_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        self.base_rat
    }
    fn tri(&self) -> [Point3; 3] {
        self.tri
    }
    fn chain_id(&self) -> u64 {
        self.base.chain_id
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        self.base.tri
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.rotated
    }
}

impl Witness for FaceRow {
    // ★ Every arm reads through [`FaceRow::plane`], which **panics on a cylinder row**: the
    // only predicate that runs on the face table is `Judge::planes_coplanar` during class
    // discovery, and that sweep filters to plane rows before asking (a cylinder's identity is
    // the cylinder class table's question). A panic here is an upstream filter bug made
    // loud, never a silently wrong plane answer.
    fn base_coeffs_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        self.plane().base_rat
    }
    fn tri(&self) -> [Point3; 3] {
        self.plane().tri
    }
    // A face table exists before plane classes do, and the only predicate that runs on it is
    // `Judge::planes_coplanar` during class discovery. Opting out here keeps that path unchanged.
    fn chain_id(&self) -> u64 {
        0
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        None
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.plane().tri_pt3
    }
    fn is_rotated(&self) -> bool {
        self.plane().rotated
    }
}

impl PlaneWitness for WorkingPlane {
    fn coeffs(&self) -> [f64; 4] {
        self.plane.coefficients()
    }
    // ★ `rotated` gates both: a moved plane's name speaks in its pre-motion frame, and the row
    // is a world description only when nothing moved the plane (the table sets `rotated` from
    // that — including a `Through` plane whose vertices meet in a frame its motion does not name).
    fn exact_coeffs(&self) -> Option<[f64; 4]> {
        self.name_ints.as_ref().filter(|_| !self.rotated)?.row
    }
    fn exact_normal(&self) -> Option<[f64; 3]> {
        self.name_ints.as_ref().filter(|_| !self.rotated)?.normal
    }
    fn frame_sign(&self) -> i8 {
        self.frame_sign
    }
    fn base_coeffs(&self) -> Option<[f64; 4]> {
        self.base.coeffs
    }
    fn name_ints(&self) -> Option<&nacre_judge::predicate::NameInts> {
        self.name_ints.as_ref()
    }
}

#[cfg(test)]
#[path = "tests/tolerant.rs"]
mod tests;
