# occt-helper — OCCT oracle backend

An **out-of-process** helper that scores nacre's geometry against OpenCASCADE
(OCCT), the dev oracle (design §7). OCCT is **never linked into the kernel** —
it is a dev-only test scorer, reached only through this helper + STEP files
(overview "OCCT / 외부 코드 규칙"). A crash on adversarial input kills the helper,
not the kernel.

## Protocol

```
occt-helper props <in.step>
```
Reads the STEP file and prints mass properties as key-value lines:
```
volume <v>
area <a>
faces <n>
bbox_min <x> <y> <z>
bbox_max <x> <y> <z>
```
Exit code: `0` success · `1` geometry failure (file missing, unreadable, no
shape) · `2` DRAWEXE crash.

(The design's `fuse|cut|common` boolean commands land with the M5 boolean
oracle; the protocol is fixed so the Rust side — `nacre-oracle` — is unchanged
when new commands or a new backend arrive.)

## Backend

`props` runs the Homebrew OCCT `DRAWEXE` Tcl shell in batch mode
(`stepread → vprops/sprops/nbshapes/bounding`) and parses its console output.
Requires `brew install opencascade` (provides `DRAWEXE` on `PATH`).

If the DRAWEXE console scraping ever proves fragile (version/output drift), swap
this script for a thin C++ helper linking OCCT (`STEPControl_Reader` +
`BRepGProp` + `BRepBndLib`) that emits the same lines — the protocol is fixed,
so nothing else changes.
