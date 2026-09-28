//! Certified indirect predicates (CIP) — the toleranced sign layer for the nacre kernel.
//!
//! Two parts, both pure numeric (neither touches `nacre-topo`):
//! - [`kernel`] — the toleranced-sign kernel: a rotated point ([`WitnessPoint`]) carries a sound
//!   directional tol, and an orientation determinant is decided by an f64 filter → astro-float
//!   escalation → a proved coincidence or an honest "undecided" ([`Decision`]). The
//!   certified-toleranced twin of `nacre-predicates` (exact f64).
//! - [`predicate`] — the plane-arrangement geometric predicates (`Judge::orient3d`, `Judge::cmp_coord`,
//!   …) that route each query to the exact path (`nacre-predicates`) or the kernel by whether
//!   its planes are rotated. The b-rep supplies witnesses through the
//!   [`Witness`](predicate::Witness) / [`PlaneWitness`](predicate::PlaneWitness) ports.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
pub mod kernel;
pub mod predicate;

// The toleranced-sign kernel symbols consumers use directly (`nacre-ops`: `WitnessPoint`/`MoveNode` for
// vertex assembly, `dir_orient3d_judge` for a per-face normal sign).
pub use kernel::frame3::{
    Decision, FrameThrough, JudgedPoint, MoveNode, Standard, WideFrame, WitnessPoint, chain_parity,
    dir_orient3d_judge, dir_sign_judge, fold_suffix, indirect_cmp_coord_judge,
    indirect_orient3d_judge, judge_precision, meet_hp, normals_agree_judge, orient3d_filter,
    orient3d_judge, plane_hp, precision_for, trial_bound,
};
