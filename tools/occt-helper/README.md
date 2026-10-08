# occt-helper — OCCT oracle backend

An **out-of-process** helper that scores nacre's geometry against OpenCASCADE
(OCCT), the dev oracle. OCCT is **never linked into the kernel** —
it is a dev-only test scorer, reached only through this helper + STEP files
(overview "OCCT / 외부 코드 규칙"). A crash on adversarial input kills the helper,
not the kernel.

## Protocol

```
occt-helper props <in.step>
occt-helper <fuse|cut|common> <a.step> <b.step> [out.step]
```
`props` reports one solid; the boolean commands run OCCT's `bfuse`/`bcut`/
`bcommon` on two STEP solids and report the **result's** properties. `cut` is
`A − B`. `out` optionally writes the boolean result as STEP (unused by the props
oracle). Every command prints the same key-value schema:
```
volume <v>
area <a>
faces <n>
bbox_min <x> <y> <z>
bbox_max <x> <y> <z>
centroid <x> <y> <z>
valid <1|0>
```
`centroid` is the volume's centre of mass. `valid` is DRAWEXE's `checkshape` on every shape
read (`a_1`, and `b_1` for a boolean) and on the boolean's result: `1` when each "seems to be
valid", `0` when any has faulty sub-shapes. A faulty shape adds one line per complaint,
`fault <a|b|r> <BRepCheck_code>`, and with `valid 0` the other keys appear only as far as they
were measured. Each `checkshape` is read in its own fenced section of the console output — one
grep over the whole output would read "input faulty, result valid" as valid. `nacre-oracle`
refuses a faulty shape (`OracleError::Faulty`).
Exit code: `0` success · `1` geometry failure (file missing, unreadable, no/empty
shape) · `2` DRAWEXE crash. The schema and exit codes are identical across
commands, so the Rust side (`nacre-oracle`) shares one parser.

Inputs must be **single-solid** STEP files (nacre's `to_step_solid`): a STEP with
multiple transferable roots reads as `x_1, x_2, …` and the boolean would use only
`x_1`. nacre's per-solid export always emits one root.

## Reading orientation as written

OCCT's STEP read repairs what it reads (`FromSTEP.exec.op : FixShape`), and the repair silently
turns a wrongly oriented face, edge, shell or void the right way round — such a file would score
valid with the right volume. The helper turns off only the three orientation fixes
(`param FromSTEP.FixShape.{FixOrientationMode,FixShellOrientationMode,FixFaceOrientationMode} 0`,
inside the DRAWEXE session). The rest of the repair builds the 2D parameter curves the files do
not carry and inserts the seam a seamless cylinder lateral needs; without it every face on a
cylinder reads invalid. A caller whose `CSF_STEPDefaults` empties `exec.op` runs no repair at all,
and the three settings then change nothing.

## Backend

Every command runs the Homebrew OCCT `DRAWEXE` Tcl shell in batch mode
(`stepread → [bfuse/bcut/bcommon →] vprops/sprops/nbshapes/bounding`) and parses
its console output. Requires `brew install opencascade` (provides `DRAWEXE` on
`PATH`).

If the DRAWEXE console scraping ever proves fragile (version/output drift), swap
this script for a thin C++ helper linking OCCT (`STEPControl_Reader` +
`BRepGProp` + `BRepBndLib`) that emits the same lines — the protocol is fixed,
so nothing else changes.
