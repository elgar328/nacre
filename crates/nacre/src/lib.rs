//! # nacre
//!
//! A pure-Rust CAD kernel that keeps an exact b-rep as the source of truth: analytic
//! surfaces and curves are never discretized into the model, meshes are derived, and a
//! tolerance is attached only to geometry the kernel *discovered* by intersection —
//! never to what a caller constructed.
//!
//! This is the **facade**. The kernel is built as one crate per layer so the layers
//! cannot reach into each other (geometry knows nothing of topology; topology refers to
//! geometry only by handle), and this crate re-exports them so a consumer depends on
//! one name instead of ten.
//!
//! ## The layers
//!
//! | module | what lives there |
//! |---|---|
//! | [`store`] | `Store<T>` / `Handle<T>` — append-only typed arenas |
//! | [`math`] | `Point<D>` / `Vector<D>` — plain f64 linear algebra |
//! | [`exact`] | `Rat` / `Angle` / `Isometry` / `PlaneName` — what is answered exactly in rationals |
//! | [`geom`] | `Plane` / `Line` / `Circle` / `Cylinder`, and the isolated `intersect` predicates |
//! | [`topo`] | `Model` and the b-rep cells — the truth aggregate |
//! | [`ops`] | `Operation` / `apply` / `replay`, sketch profiles, extrude, pad, pocket, boolean |
//! | [`props`] | volume, area, centroid, bounding box, per-face facts |
//! | [`validate`] | b-rep invariant checking |
//! | [`tess`] | triangulation and OBJ export |
//! | [`step`] | STEP (AP242) export |
//!
//! Two names need care:
//!
//! - **`Surface` exists twice.** `topo::Surface` is the *truth* — the exact definition a
//!   face's handle points at; `geom::Surface` is its f64 *realization*, the cache the
//!   numeric layers read. The [`prelude`] carries neither; name the module.
//! - **`validate` is both a module and a function.** [`validate`] is the layer, so the
//!   checker is `validate::validate`. The [`prelude`] carries the function.
//!
//! ## Example
//!
//! ```
//! use nacre::prelude::*;
//!
//! // A square profile, as a closure so the example can make two of them.
//! let square = |a: f64, b: f64| {
//!     Profile2d::polygon(vec![
//!         Point2::from_array([a, a]),
//!         Point2::from_array([b, a]),
//!         Point2::from_array([b, b]),
//!         Point2::from_array([a, b]),
//!     ]).unwrap()
//! };
//! // An extrude *names* the plane it sketches on — `SketchFrame::world` is the frame twin of
//! // `SketchPlane::world_xy`, reading the plane `Model::new` already seeded.
//! let extrude = |m: &mut Model, profile, dist| {
//!     let frame = SketchFrame::world(m, Axis::Z);
//!     match apply(m, &Operation::Extrude { frame, profile, dist }) {
//!         Ok(OpOutput::Extrude { solid, .. }) => solid,
//!         other => panic!("{other:?}"),
//!     }
//! };
//!
//! // Sweep a 4x4 plate 1 tall, then a 2x2 post 3 tall standing in the middle of it.
//! let mut model = Model::new();
//! let plate = extrude(&mut model, square(0.0, 4.0), 1.0);
//! let post = extrude(&mut model, square(1.0, 3.0), 3.0);
//!
//! // Union them. Both operands are consumed; the result is the new live solid.
//! let joined = boolean(&mut model, BoolKind::Fuse, plate, post).unwrap();
//! model.rebuild_adjacency();
//!
//! // The result is a valid closed b-rep, and its volume is the plate plus the part of
//! // the post standing above it: 4*4*1 + 2*2*2.
//! assert!(validate(&model).is_empty());
//! let volume = mass_props(&model, joined[0]).unwrap().volume;
//! assert!((volume - 24.0).abs() < 1e-9, "{volume}");
//! ```
//!
//! Errors are named, never silent: an operation the kernel cannot compute *correctly*
//! returns a reason (`ops::RejectReason`, classified by `ops::RejectClass`) rather than
//! a plausible-looking wrong answer.
#![doc(html_root_url = "https://docs.rs/nacre")]
#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]

pub use nacre_exact as exact;
pub use nacre_geom as geom;
pub use nacre_math as math;
pub use nacre_ops as ops;
pub use nacre_props as props;
pub use nacre_step as step;
pub use nacre_store as store;
pub use nacre_tess as tess;
pub use nacre_topo as topo;
pub use nacre_validate as validate;

/// The names a modelling script reaches for constantly, in one import.
///
/// Everything else is a module path away (`crate::geom`, `crate::step`, the error
/// detail types). The bar for being in here is empirical: it covers what a real
/// consumer imported, plus every name needed to **construct any `ops::Operation`
/// variant** — so `use nacre::prelude::*` is enough to drive the kernel.
pub mod prelude {
    pub use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    pub use nacre_math::{Point2, Point3, Vector3};
    pub use nacre_ops::{
        BoolError, BoolKind, DatumDef, Edge2d, OpError, OpOutput, Operation, Profile2d,
        RejectClass, RejectReason, RejectWhere, Ring2d, SketchError, SketchFrame, SketchPlane,
        apply, arc_to_rat, arc_turns, arc_turns_rat, boolean, face_plane, face_sketch_frame,
        from_paths, from_rings, replay,
    };
    pub use nacre_props::{
        FaceProps, MassProps, bounds, centroid, face_normal_at, face_props, mass_props,
    };
    pub use nacre_store::Handle;
    pub use nacre_tess::{TessConfig, Tessellation, tessellate};
    pub use nacre_topo::{Face, FramePlacement, Model, Solid};
    pub use nacre_validate::validate;
}
