"""Time FreeCAD on the fin model, feature by feature, headless.

    /Applications/FreeCAD.app/Contents/Resources/bin/freecadcmd tools/freecad-bench.py

**Why FreeCAD and not just DRAWEXE.** They are the same engine, so DRAWEXE answers "what does
OCCT cost". Only FreeCAD answers "what does a *user* of OCCT wait for" -- a PolarPattern is one
click that hides an 80-way boolean, and that is the thing a person compares us against.

Two things this measures that a stopwatch in the GUI cannot:

  - **Wall clock with the cache actually cold.** A recompute is a no-op unless the features are
    `touch()`ed first, which is why the Report view stays empty and the GUI feels instant.
  - **CPU time next to it.** OCCT runs its booleans on many threads, and the ratio says whether a
    number is work or waiting.

Measured (14-core M-series, FreeCAD 1.1.0):

    PolarPattern x80    wall 3.74 s    CPU 21.0 s    -> 5.6 cores average
    full recompute      wall 4.13 s

The 5.6 cores are mostly *spin*, not work: DRAWEXE does the same fuse in 3.01 s on one thread
(`brunparallel 0`, CPU 3.03 s == wall), and turning OCCT's parallelism on buys 1.05x while
spending 4x the CPU. So FreeCAD is not fast here because it is parallel. Its single-threaded
logic is simply quicker than ours.
"""

import time

import FreeCAD as App

N = 80


def timed(doc, label, feat, body):
    w0, c0 = time.perf_counter(), time.process_time()
    doc.recompute()
    w, c = time.perf_counter() - w0, time.process_time() - c0
    body.Tip = feat
    sh = feat.Shape
    print(
        f"  {label:<22}wall {w:7.3f} s   CPU {c:7.3f} s ({c / max(w, 1e-9):5.2f} cores)"
        f"   -> {len(sh.Faces):>5} faces, vol {sh.Volume:9.4f}",
        flush=True,
    )


def main():
    doc = App.newDocument("bench")
    body = doc.addObject("PartDesign::Body", "Body")
    doc.recompute()

    # The fin: a box reaching out from x=0.5, which puts every copy inside its neighbours near
    # the axis. That overlap is the whole cost -- 80 disjoint copies would be nearly free.
    fin = doc.addObject("PartDesign::AdditiveBox", "Fin")
    body.addObject(fin)
    fin.Length, fin.Width, fin.Height = 3.5, 0.4, 1.0
    fin.Placement = App.Placement(App.Vector(0.5, -0.2, 1.0), App.Rotation())
    timed(doc, "fin (pad)", fin, body)

    pat = doc.addObject("PartDesign::PolarPattern", "PolarPattern")
    body.addObject(pat)
    pat.Originals = [fin]
    pat.Axis = (body.Origin.OriginFeatures[2], [""])  # index 2 is Z_Axis
    pat.Angle, pat.Occurrences = 360.0, N
    timed(doc, f"PolarPattern x{N}", pat, body)

    hub = doc.addObject("PartDesign::AdditiveBox", "Hub")
    body.addObject(hub)
    hub.BaseFeature = pat  # without this the feature chain does not advance and `body.Shape`
    hub.Length, hub.Width, hub.Height = 2.0, 2.0, 3.0  # silently reports the previous tip
    hub.Placement = App.Placement(App.Vector(-1, -1, 0), App.Rotation())
    timed(doc, "hub (pad)", hub, body)

    bore = doc.addObject("PartDesign::SubtractiveBox", "Bore")
    body.addObject(bore)
    bore.BaseFeature = hub
    bore.Length, bore.Width, bore.Height = 1.2, 1.2, 3.0
    bore.Placement = App.Placement(App.Vector(-0.6, -0.6, 0), App.Rotation())
    timed(doc, "bore (pocket)", bore, body)

    for o in doc.Objects:
        o.touch()  # a recompute without this is a no-op -- the shapes are cached
    w0, c0 = time.perf_counter(), time.process_time()
    doc.recompute(None, True)
    w, c = time.perf_counter() - w0, time.process_time() - c0
    print(f"  {'full recompute':<22}wall {w:7.3f} s   CPU {c:7.3f} s", flush=True)


main()
