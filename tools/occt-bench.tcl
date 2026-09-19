# Time OCCT on the same fold `tests/perf.rs` times, on the same machine.
#
#   DRAWEXE -f tools/occt-bench.tcl
#
# **Why a script and not the `occt-helper` oracle.** That helper spawns one DRAWEXE per
# boolean and round-trips STEP each time, which is right for *scoring* one result and wrong
# for timing eighty: it would measure process start-up and STEP parsing. Everything here
# happens inside one interpreter, so the number is the boolean's.
#
# Both structures are timed, because they are different algorithms and the comparison is
# meaningless without saying which one:
#   - **pairwise** — 80 successive fuses, which is what nacre's fold does today
#   - **n-ary**    — one `BOPAlgo` call with all 80 arguments, which nacre has no equivalent of
#
# Measured (14-core M-series, OCCT from Homebrew, DRAWEXE single-threaded):
#
#     nacre serial          11.51 s     407 faces
#     nacre parallel         2.67 s     407 faces
#     OCCT pairwise          4.65 s   1,974 faces
#     OCCT n-ary            *1.54 s*  1,974 faces
#
# Volumes agree (nacre 237.211567, OCCT 237.212), so this is the same answer and the times
# may be compared. Three things it says, and the third is the one worth acting on:
#
#   1. **nacre's face set is 4.8x smaller.** Coplanar merging is a mandatory post-op cleaning
#      here and OCCT leaves the arrangement's pieces, so nacre does strictly more per boolean.
#   2. **Like for like — same algorithm, one thread each — nacre is 2.5x slower.** That is
#      what certified arithmetic costs against a global 1e-7 tolerance. It is a factor, not
#      an order of magnitude, and parallelism already more than covers it.
#   3. **OCCT's n-ary beats its own pairwise by 3x.** That is a *realized* number from a
#      mature implementation, on this exact model — not an upper bound. nacre has no n-ary
#      path, and this is the measurement that says building one is worth it.

pload MODELING
set N 80

# Exactly `perf.rs::fin_fold`: a 6x6x2 hub with N fins arrayed around it. Keep the two in
# step -- a benchmark comparing two different shapes says nothing.
box hub -3 -3 0 6 6 2
for {set i 0} {$i < $N} {incr i} {
    box fin$i 2 -0.4 0 6 0.8 1
    trotate fin$i 0 0 0 0 0 1 [expr {360.0 * $i / $N}]
}

# ---- pairwise: the structure nacre's fold has ----
# Each result gets its own name; feeding a shape back in as its own result silently produces
# nothing measurable.
copy hub acc0
set t0 [clock milliseconds]
for {set i 0} {$i < $N} {incr i} {
    bfuse acc[expr {$i + 1}] acc$i fin$i
}
set t1 [clock milliseconds]
puts "OCCT pairwise ms [expr {$t1 - $t0}]"
# `nbshapes`/`vprops` write on Draw's own channel, which is lost when stdout is redirected --
# wrap them in `puts` or the numbers silently vanish from a captured log.
puts "OCCT pairwise [string map {\n { }} [nbshapes acc$N]]"
puts "OCCT pairwise [string map {\n { }} [vprops acc$N]]"

# ---- n-ary: one arrangement over all arguments ----
bclearobjects
bcleartools
baddobjects hub
for {set i 0} {$i < $N} {incr i} { baddtools fin$i }
set t2 [clock milliseconds]
bfillds
bbop res 1
set t3 [clock milliseconds]
puts "OCCT n-ary ms [expr {$t3 - $t2}]"
puts "OCCT n-ary [string map {\n { }} [nbshapes res]]"
puts "OCCT n-ary [string map {\n { }} [vprops res]]"
