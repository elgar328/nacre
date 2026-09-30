//! Boolean **capability coverage** suite — the kernel's acceptance matrix, driven
//! through the public API (front door), grouped by scenario. Each submodule is one
//! theme of "shape × operation → expected topology", asserted against the
//! independent `nacre_props`/`nacre_validate` oracles. Finer checks of private
//! helpers live in each source module's own `#[cfg(test)] mod tests`.
//!
//! One integration binary (submodules under `coverage/`, referenced by `#[path]`)
//! rather than one binary per theme, to avoid re-linking the crate per file.

#[path = "support/fixtures.rs"]
mod common;
#[path = "support/stated.rs"]
mod stated;

#[path = "coverage/convex.rs"]
mod convex;

#[path = "coverage/multisolid.rs"]
mod multisolid;

#[path = "coverage/rotation.rs"]
mod rotation;

#[path = "coverage/placement.rs"]
mod placement;

#[path = "coverage/edge_on_a_ruling.rs"]
mod edge_on_a_ruling;

#[path = "coverage/nonconvex.rs"]
mod nonconvex;

#[path = "coverage/features.rs"]
mod features;

#[path = "coverage/coplanar.rs"]
mod coplanar;

#[path = "coverage/invariants.rs"]
mod invariants;

#[path = "coverage/rejects.rs"]
mod rejects;

#[path = "coverage/copy.rs"]
mod copy;

#[path = "coverage/mirror.rs"]
mod mirror;

#[path = "coverage/sketch.rs"]
mod sketch;

#[path = "coverage/report.rs"]
mod report;

#[path = "coverage/fin_array.rs"]
mod fin_array;

#[path = "coverage/arc_crossings.rs"]
mod arc_crossings;
