# Nacre

**Nacre is an exact b-rep CAD kernel.** Every topological decision is exact, and it never fails silently.

- **Exact by definition.** Dimensions keep the decimals you typed: `1.1` is 11/10, not the nearest binary fraction. A vertex is defined by the surfaces that meet there, and its coordinates are only a cache with a proven error bound, recomputable at any precision.

- **No global tolerance.** Many CAD kernels treat two points as the same whenever they fall within one global epsilon. Nacre makes every decision from the exact definitions instead. A fast f64 check settles most cases, and the rest, including those where a rotation has made coordinates irrational, are recomputed with as much precision as the geometry requires.

- **Named errors, not wrong solids.** When an operation cannot be computed correctly, Nacre refuses it and names the reason: the case is outside what the kernel supports, it has no valid result, or it points to a bug in the kernel itself.

- **Pure Rust, down to the arithmetic.** There are no C dependencies, so the same kernel runs natively and in the browser through WebAssembly.

> Named for *nacre*, mother-of-pearl, which grows one layer at a time and never rewrites the layers beneath. The kernel treats geometry the same way and, unlike most CAD kernels, its topology too.

## Try it

[![Parts built with nacre in the playground](https://github.com/elgar328/nacre-playground/releases/download/assets/showcase.png)](https://elgar328.github.io/nacre-playground/gallery/)

The [playground](https://github.com/elgar328/nacre-playground) runs the kernel right in your browser: write a short script and watch it build the exact solid. The [gallery](https://elgar328.github.io/nacre-playground/gallery/) collects example parts, rebuilt every day from the latest kernel. Click one to open it in the playground.

## Status

Nacre is at an early stage. While the version is `0.0.z`, any release may change anything. The kernel grows milestone by milestone, from planar solids to quadric surfaces to general NURBS, and each milestone leaves a working kernel behind.

`M1 ✓   M2 ✓   M3 ✓   M4 ✓   M5 ✓   M6 ◔   M7 ○`

| Milestone | Coverage | Status |
|---|---|---|
| **M1–M4**: foundation | Exact storage and handles, sketches with lines and arcs, extrude either way off a plane or a face, datum planes, rigid motions and mirroring, replay, validation, mass properties, tessellation, shape-only STEP export | **Done** |
| **M5**: polyhedral boolean | Booleans on planar solids, including coplanar contact, containment, cavities, multiple bodies and rotated operands | **Done** |
| **M6**: quadric boolean | Planes with cylinders, spheres and cones | **In progress**: cylinders work against planes perpendicular or parallel to their axis |
| **M7**: general NURBS | Intersections of free-form surfaces | Research track |

Anything outside the current coverage is refused, never approximated.

## Scope

Nacre is a geometry kernel, not a CAD application. It covers the geometry and topology of solids, modeling operations, tessellation and validation. Everything else a CAD application needs, such as GD&T, appearance or assemblies, is outside its scope and belongs to the application.

`nacre-step` exports compact STEP files that carry only the shape. An application that wants complete STEP files, with product data and the rest, should use its own exporter instead.

Two separate repositories put the kernel to work: [nacre-kit](https://github.com/elgar328/nacre-kit), a convenience layer on top of it, and [nacre-playground](https://github.com/elgar328/nacre-playground), a browser app built on the kit.

## Relation to Fornjot

Nacre began after [Fornjot](https://github.com/hannobraun/fornjot), Hanno Braun's b-rep CAD kernel in Rust, which was developed from 2020 until it was shut down in 2026. The two share much of their outlook: code-first mechanical CAD, immutable objects referenced by handle, and clear errors instead of quietly wrong results. No code is shared.

Where they differ is in what counts as the truth. Fornjot's [final experiment](https://github.com/hannobraun/fornjot/tree/main/experiments/2025-12-03) made approximated geometry the uniform representation and recorded topology alongside it as the geometry was built, deliberately trading exactness for simplicity. Nacre makes the opposite trade. Exact geometry is the truth, and approximations are a cache derived from it.

## License

Licensed under either MIT or Apache-2.0, at your option.
