//! Boolean **capability coverage** suite — the kernel's acceptance matrix, driven
//! through the public API (front door), grouped by scenario. Each submodule is one
//! theme of "shape × operation → expected topology", asserted against the
//! independent `nacre_props`/`nacre_validate` oracles. Finer checks of private
//! helpers live in each source module's own `#[cfg(test)] mod tests`.
//!
//! One integration binary (submodules under `coverage/`, referenced by `#[path]`)
//! rather than one binary per theme, to avoid re-linking the crate per file.

#[path = "coverage/common.rs"]
mod common;

#[path = "coverage/convex.rs"]
mod convex;

#[path = "coverage/multisolid.rs"]
mod multisolid;

#[path = "coverage/rotation.rs"]
mod rotation;

#[path = "coverage/nonconvex.rs"]
mod nonconvex;

#[path = "coverage/features.rs"]
mod features;
