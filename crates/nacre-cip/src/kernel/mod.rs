//! The toleranced-sign kernel — certified indirect predicates over `Pt2`/`Pt3` (design.md
//! §9 CIP). The 2D ([`frame2`]) and 3D ([`frame3`]) frames each carry a rotated point with a
//! sound directional tol and decide the sign of an orientation determinant: an f64 filter
//! with a sound error bound handles the easy cases, the ambiguous ones escalate to
//! astro-float from the exact point definitions, and one that still cannot separate from zero is
//! turned into the **distance it stands for** and compared against the operation's coincidence
//! limit ([`frame3::Standard`]) — proved coincident, or said out loud ([`frame3::Decision`]) rather
//! than assumed. Pure numeric layer — depends only on `nacre-scalar` (Rat/Angle/Orient) and
//! astro-float, never on `nacre-math`/`nacre-topo`.

use astro_float::RoundingMode;

/// Rounding mode for the high-precision (astro-float) realization layer.
pub(crate) const HP_RM: RoundingMode = RoundingMode::ToEven;

pub mod frame2;
pub mod frame3;
pub(crate) mod interval;
