# Changelog

All notable changes to the `nacre` facade are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
While the version is `0.0.z`, every release may break anything.

## [Unreleased]

### Changed

- `Operation::Extrude` takes a signed distance: a negative one sweeps against the frame's normal, in the same frame, with the base cap still on the frame's plane. A distance of zero is `OpError::ZeroDistance`.
- **Breaking:** `topo::Model::push_edge`, `topo::Model::derive_edge_curve` and `topo::Model::rebuild_edge_cache` take a closure that gives, when asked, what the caller realized of the edge's truth (`topo::EdgeGiven`) — pass `|_| EdgeGiven::NONE` to give nothing.

### Removed

- `Operation::PadOnFace` and `Operation::PocketOnFace`, their `OpOutput` variants, and `OpError::PadMissesFace`, `OpError::PocketNotBlind` and `OpError::NonPositiveDistance`. A pad or pocket is the face's sketch frame (`face_sketch_frame`), an `Extrude` with a positive or negative distance, and a `Boolean` with the face's solid; a pocket deeper than the body cuts through.
- `LogCell::Face`, which only the removed operations used.
- `topo::Model::make_live`, whose only callers were the removed operations.

### Fixed

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
