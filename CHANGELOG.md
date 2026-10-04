# Changelog

All notable changes to the `nacre` facade are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
While the version is `0.0.z`, every release may break anything.

## [Unreleased]

### Changed

- An extruded circle's side face has no seam line: it is bounded by its two rims, one its outer loop and the other an inner loop. The solid has two edges instead of three, and its STEP output has no seam `LINE`.
- `Operation::Extrude` takes a signed distance: a negative one sweeps against the frame's normal, in the same frame, with the base cap still on the frame's plane. A distance of zero is `OpError::ZeroDistance`.
- **Breaking:** `topo::Model::push_edge`, `topo::Model::derive_edge_curve` and `topo::Model::rebuild_edge_cache` take a closure that gives, when asked, what the caller realized of the edge's truth (`topo::EdgeGiven`) — pass `|_| EdgeGiven::NONE` to give nothing.
- **Breaking:** `topo::Model::push_plane`, `topo::Model::push_plane_through` and `topo::Model::push_cylinder` take, after the figure, what that figure knows about itself (`topo::CacheStanding`) — pass `CacheStanding::Unrealized` for a figure nobody realized; a surface's cache now says whether it is the truth's realization (`topo::Model::surface_cache_standing`).
- **Breaking:** `ops::refine_vertex_cache` is now `ops::refine_caches`, the door to call before exporting a model with a long history: it also raises surface caches and re-derives edges, and its `RefineReport` counts, for vertices, surfaces and edges, what it raised and what it left — undecided at the top of the ladder, or with no way to realize it.

### Removed

- `Operation::PadOnFace` and `Operation::PocketOnFace`, their `OpOutput` variants, and `OpError::PadMissesFace`, `OpError::PocketNotBlind` and `OpError::NonPositiveDistance`. A pad or pocket is the face's sketch frame (`face_sketch_frame`), an `Extrude` with a positive or negative distance, and a `Boolean` with the face's solid; a pocket deeper than the body cuts through.
- `LogCell::Face`, which only the removed operations used.
- `topo::Model::make_live`, whose only callers were the removed operations.

### Fixed

- A cylinder bitten at both ends — a box over part of each rim, fused on or cut away — now tessellates. Its side face comes out with no seam line, bounded by its two cut rims, and the mesher refused it, so the whole model failed to mesh; the mesher now cuts such a face open along one line of the cylinder.
- Such a side face with windows whose spans together go all the way round the axis now tessellates too: the mesher cuts it open through one window rather than refusing the face.
- A cylinder crossing the side of such a cylinder, or a slanted wall cutting it, is now refused as `CylinderPairContact` or `ObliqueCylinderCut`; the side face's height was read from its lower rim alone, the contact was judged clear, and a later stage refused it as `LabelConflict` or `CylinderStagesDisagree`.
- A plane through three vertices that a transform moved (a datum plane, or a face built on one) is now exported at the nearest `f64` of its true position and orientation, unless its history is longer than 192 recorded motions behind a turn that is not a multiple of 90° or a sketch plane with no rational frame; before, it carried the previous position moved in floating point.
- A circular edge in STEP export (and in tessellation and mass properties) now carries its cylinder's axis, reference direction and radius bit for bit, and its centre is the nearest `f64` of the exact centre — including cylinders turned by an angle that is not a multiple of 90° or built on a sketch plane with no rational frame, unless that history is longer than 192 recorded motions; before, these values could be a few units in the last place off.
- A straight edge in STEP export is written with a unit `VECTOR` magnitude, and its direction is the nearest `f64` of the exact direction — including edges on faces turned by an angle that is not a multiple of 90° or built on a sketch plane with no rational frame, unless that history is longer than 192 recorded motions; before, the direction and the length were recomputed from the rounded end points.

## [0.0.1] - 2026-10-02

### Added

- First release of the kernel itself: exact b-rep solids built from sketches and modeling operations, all through the `nacre` crate.
- Sketches of lines and arcs on datum planes; the kernel works out which rings are holes and which are separate islands.
- Extrude, pad and pocket, translation and rotation by angles in rational degrees, mirroring and copying.
- Cut, fuse and common on planar solids, including coplanar contact, containment, cavities, multiple bodies and rotated operands.
- Booleans that mix planar faces with cylinders, where each plane is perpendicular or parallel to the cylinder's axis.
- Replay: the same operation log always rebuilds the same model.
- An operation that cannot be computed correctly fails with a named reason: not supported, impossible, or a suspected kernel defect.
- B-rep validation, exact mass properties (volume, area, centroid), tessellation with OBJ output, and compact shape-only STEP (AP242) export.
- Pure Rust with no C dependencies. The default `parallel` feature runs work on rayon; turn it off to build for `wasm32-unknown-unknown`.

## [0.0.0] - 2026-07-07

### Added

- Placeholder that reserves the `nacre` name on crates.io; it contains no kernel.

[Unreleased]: https://github.com/elgar328/nacre/compare/v0.0.1...HEAD
[0.0.1]: https://github.com/elgar328/nacre/releases/tag/v0.0.1
[0.0.0]: https://crates.io/crates/nacre/0.0.0
