# Changelog

All notable changes to the `nacre` facade are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
While the version is `0.0.z`, every release may break anything.

## [Unreleased]

### Changed

- `Operation::Extrude` takes a signed distance: a negative one sweeps against the frame's normal, in the same frame, with the base cap still on the frame's plane. A distance of zero is `OpError::ZeroDistance`.

### Removed

- `Operation::PadOnFace` and `Operation::PocketOnFace`, their `OpOutput` variants, and `OpError::PadMissesFace`, `OpError::PocketNotBlind` and `OpError::NonPositiveDistance`. A pad or pocket is the face's sketch frame (`face_sketch_frame`), an `Extrude` with a positive or negative distance, and a `Boolean` with the face's solid; a pocket deeper than the body cuts through.
- `LogCell::Face`, which only the removed operations used.
- `topo::Model::make_live`, whose only callers were the removed operations.

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
