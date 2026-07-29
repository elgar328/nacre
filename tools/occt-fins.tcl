# Time OCCT on the playground's fin model, with the operand **order** as a variable.
#
#   DRAWEXE -f tools/occt-fins.tcl
#
# **Why a second script next to `occt-bench.tcl`.** That one folds a hub with fins reaching from
# r=2 to r=8. This one is the model a user actually wrote in the playground: fins from r=0.5, so
# every copy sits deep inside its neighbours near the axis. The two fixtures disagree about
# n-ary, which is exactly why both are kept -- one fixture is not a finding.
#
# Measured 2026-07-29 (14-core M-series, OCCT from Homebrew, `brunparallel 0`). Each pair is the
# same answer, so the times may be compared:
#
#                                    pairwise    n-ary    n-ary buys
#     fins only        (vol 49.5171)    8.27 s   3.03 s      2.7x
#     hub first + fins (vol 58.3029)    3.64 s   3.49 s      1.04x
#
# **The finding is the row difference, not the column difference.** Fusing the hub *first* makes
# the pairwise fold 2.3x faster (8.27 -> 3.64) all by itself, and once the order is right, n-ary
# buys 4%. So most of what an n-ary boolean is credited with here is really *avoiding a bad
# accumulation order*: a big central body absorbs the overlapping inner ends, and every later
# fuse meets a simpler accumulation.
#
# That matters for a code-CAD, where the script writes the fold as a `for` loop and the kernel
# never sees the whole set at once. The good order is the one a person writes anyway --
# `part = hub; for fin { part = fuse(part, fin) }` -- so the cheap structure is already what we
# get. nacre shows the same effect, smaller: 14.88 s fins-only vs 11.0 s hub-first.
#
# Threads, since it was assumed once and the assumption was wrong: DRAWEXE's default really is
# single-threaded (CPU == wall), and `brunparallel 1` buys 1.05x for 4x the CPU -- almost all of
# those cores spin. FreeCAD's PolarPattern shows the same shape (wall 3.74 s, CPU 21.0 s). So
# FreeCAD is not quick here because it is parallel; its single-threaded logic is. Like for like
# on one thread nacre is **4.9x** slower (14.88 s vs 3.03 s), while nacre's own parallelism is
# the outlier that works -- 4.25x against OCCT's 1.05x, closing the gap to 1.2x.

pload MODELING
brunparallel 0
set N 80
box hub -1 -1 0 2 2 3
for {set i 0} {$i < $N} {incr i} {
    box fin$i 0.5 -0.2 1 3.5 0.4 1
    trotate fin$i 0 0 0 0 0 1 [expr {360.0 * $i / $N}]
}

# ---- A: fins only, pairwise (what a code-CAD `for` loop produces) ----
copy fin0 a0
set t0 [clock milliseconds]
for {set i 1} {$i < $N} {incr i} { bfuse a$i a[expr {$i - 1}] fin$i }
set t1 [clock milliseconds]
puts "A fins pairwise ms [expr {$t1 - $t0}]"
puts "A [string map {\n { }} [vprops a[expr {$N - 1}]]]"

# ---- B: fins only, n-ary ----
bclearobjects
bcleartools
baddobjects fin0
for {set i 1} {$i < $N} {incr i} { baddtools fin$i }
set t2 [clock milliseconds]
bfillds
bbop bres 1
set t3 [clock milliseconds]
puts "B fins n-ary ms [expr {$t3 - $t2}]"
puts "B [string map {\n { }} [vprops bres]]"

# ---- C: hub FIRST, then the fins pairwise -- same answer, different order ----
copy hub c0
set t4 [clock milliseconds]
for {set i 0} {$i < $N} {incr i} { bfuse c[expr {$i + 1}] c$i fin$i }
set t5 [clock milliseconds]
puts "C hub-first pairwise ms [expr {$t5 - $t4}]"
puts "C [string map {\n { }} [vprops c$N]]"

# ---- D: hub + fins, n-ary ----
bclearobjects
bcleartools
baddobjects hub
for {set i 0} {$i < $N} {incr i} { baddtools fin$i }
set t6 [clock milliseconds]
bfillds
bbop dres 1
set t7 [clock milliseconds]
puts "D hub+fins n-ary ms [expr {$t7 - $t6}]"
puts "D [string map {\n { }} [vprops dres]]"
