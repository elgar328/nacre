//! **Where the boolean's wall clock is measured, so that speed claims come from a number.**
//!
//! `#[ignore]`d, and it asserts nothing about time: a timing assertion is flaky on a shared
//! machine and would eventually be relaxed into meaninglessness. The values belong in the
//! commit message of whatever change moved them, next to the census diff that says the
//! answers did not change.
//!
//! ★★★ **Run it `--release`, and it is skipped by the routine `--ignored` sweep**.
//! Measured: this file's two tests were **1673 of the sweep's 1802 seconds — 93%** — and they ran
//! in a *debug* build, so the numbers they printed were not the ones anyone wants. `--release`
//! turns 28 minutes into 4, and the gate skips them by name (`overview.md` carries the command).
//!
//! ★★ **But they are coverage, not only measurement, and that is what the split costs.** Nothing
//! here asserts a duration, yet every boolean is `.expect("fuse")`ed and every fold's volume is
//! taken — so a large *rotated* fold that stopped building would panic here and nowhere else in
//! the suite at this scale. What `--release` gives up is the **`debug_assert!`** net over exactly
//! those folds. So run this file in **debug** now and then (`cargo test -p nacre-ops --test perf
//! -- --ignored`, ~28 min) — not every cell, but before believing the rotated path is healthy.
//!
//! Two workloads, because parallelism pulls them in opposite directions:
//!
//! - **Rotated fold** — a hub with `n` fins arrayed around it. Rotation puts every judgement
//!   on the toleranced CIP path, plane classes grow into the dozens, and per-class work is
//!   milliseconds. This is where more cores help.
//! - **Axis-aligned tooth fold** — the same fold with exact coordinates, so the pair brackets
//!   what the certified path costs (measured 2-4x) against what the *algorithm* costs. Both
//!   grow at about `N^1.9`, which is the finding that matters: it is not a precision problem.
//! - **Axis-aligned small booleans** — a cube fused with a cube, many times. Exact
//!   coordinates, ~11 plane classes, per-class work in microseconds. This is where the
//!   scheduling overhead can cost more than it buys, so it is measured separately rather
//!   than averaged away.
//!
//! **The first run of anything here is discarded.** Under `parallel` the first boolean pays
//! rayon's pool start-up *and* — less obviously — `nacre-scalar`'s `HP_CONSTS` is a
//! `thread_local`, so every worker builds its own π at the model's precision the first time
//! it judges. Both are once-per-process; measuring them as if they were per-boolean is how
//! one concludes that parallelism is slow.

use nacre_math::{Point2, Point3};
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, Operation, Profile2d, apply, boolean};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

fn extrude(m: &mut Model, ring: Vec<Point2>, dist: f64) -> Handle<Solid> {
    let out = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile: Profile2d::polygon(ring).unwrap(),
            dist,
        },
    )
    .expect("extrude");
    m.rebuild_adjacency();
    match out {
        nacre_ops::OpOutput::Extrude { solid, .. } => solid,
        o => panic!("{o:?}"),
    }
}

/// The same fold, **axis-aligned**: a spine with `n` teeth fused on one at a time. Same shape
/// of work (n booleans against a growing solid, a few faces added each time) with every
/// coordinate exactly representable, so the pair of them separates *the certified path's cost*
/// from *the algorithm's cost*. Measured: the exact one is only 2-4x faster, which is how we
/// know the growth is algorithmic and not a precision problem.
fn tooth_fold(n: i128) -> std::time::Duration {
    let mut m = Model::new();
    let w = n as f64 * 0.5 + 1.0;
    let mut acc = extrude(
        &mut m,
        vec![p2(-1.0, -1.0), p2(w, -1.0), p2(w, 1.0), p2(-1.0, 1.0)],
        3.0,
    );
    let mut spent = std::time::Duration::ZERO;
    for i in 0..n {
        let x = i as f64 * 0.5;
        let tooth = extrude(
            &mut m,
            vec![p2(x, 0.5), p2(x + 0.2, 0.5), p2(x + 0.2, 4.0), p2(x, 4.0)],
            2.0,
        );
        let t = std::time::Instant::now();
        acc = boolean(&mut m, BoolKind::Fuse, acc, tooth).expect("fuse")[0];
        spent += t.elapsed();
        m.rebuild_adjacency();
    }
    spent
}

fn turn(m: &mut Model, s: Handle<Solid>, deg: Rat) -> Handle<Solid> {
    let out = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(deg).expect("angle"),
            }),
        },
    )
    .expect("rotate");
    m.rebuild_adjacency();
    match out {
        nacre_ops::OpOutput::Transform { solid } => solid,
        o => panic!("{o:?}"),
    }
}

/// `n` fins arrayed around a hub, folded one at a time. Returns the time spent in the
/// booleans alone — building the operands is not what is being measured.
fn fin_fold(n: i128) -> std::time::Duration {
    fin_fold_full(n).0
}

/// The fold, plus the answer it produced — face count and volume. Those are what license the
/// comparison in `tools/occt-bench.tcl`: two timings of *different* answers say nothing, and
/// OCCT leaves an unmerged face set where nacre merges coplanar faces, so the counts differ
/// even when the volumes agree.
fn fin_fold_full(n: i128) -> (std::time::Duration, usize, f64) {
    let mut m = Model::new();
    let mut acc = extrude(
        &mut m,
        vec![p2(-3.0, -3.0), p2(3.0, -3.0), p2(3.0, 3.0), p2(-3.0, 3.0)],
        2.0,
    );
    let mut spent = std::time::Duration::ZERO;
    for i in 0..n {
        let fin = extrude(
            &mut m,
            vec![p2(2.0, -0.4), p2(8.0, -0.4), p2(8.0, 0.4), p2(2.0, 0.4)],
            1.0,
        );
        let fin = turn(&mut m, fin, Rat::new(360 * i, n).expect("angle"));
        let t = std::time::Instant::now();
        acc = boolean(&mut m, BoolKind::Fuse, acc, fin).expect("fuse")[0];
        spent += t.elapsed();
        m.rebuild_adjacency();
    }
    let faces = m.shell(m.solid(acc).outer).faces.len();
    let volume = nacre_props::mass_props(&m, acc).expect("props").volume;
    (spent, faces, volume)
}

/// The playground's fin model, folded **two ways**, because the order is a variable and
/// `tools/occt-fins.tcl` found it to be a bigger one than n-ary.
///
/// - `hub_first = false` — the 80 fins fused to each other, nothing in the middle.
/// - `hub_first = true`  — the hub fused first, which is what the script a person writes does.
///
/// Both produce the same fins; the hub only adds itself. Measured on OCCT, fusing the hub
/// first makes the *pairwise* fold 2.3x faster and leaves n-ary buying 4% — so the cheap
/// structure is already the one a `for` loop produces, and n-ary is not the lever it looked
/// like. nacre shows the same effect, smaller.
fn fin_ring(n: i128, hub_first: bool) -> (std::time::Duration, usize, f64) {
    let mut m = Model::new();
    let fin_ring =
        |x0: f64, y0: f64, x1: f64, y1: f64| vec![p2(x0, y0), p2(x1, y0), p2(x1, y1), p2(x0, y1)];
    let mut spent = std::time::Duration::ZERO;
    let (mut acc, start) = if hub_first {
        (extrude(&mut m, fin_ring(-1.0, -1.0, 1.0, 1.0), 3.0), 0)
    } else {
        (extrude(&mut m, fin_ring(0.5, -0.2, 4.0, 0.2), 1.0), 1)
    };
    for i in start..n {
        let fin = extrude(&mut m, fin_ring(0.5, -0.2, 4.0, 0.2), 1.0);
        let fin = turn(&mut m, fin, Rat::new(360 * i, n).expect("angle"));
        let t = std::time::Instant::now();
        acc = boolean(&mut m, BoolKind::Fuse, acc, fin).expect("fuse")[0];
        spent += t.elapsed();
        m.rebuild_adjacency();
    }
    let faces = m.shell(m.solid(acc).outer).faces.len();
    let volume = nacre_props::mass_props(&m, acc).expect("props").volume;
    (spent, faces, volume)
}

/// `reps` independent small axis-aligned fuses — the shape the census is full of, and the
/// one where per-class work is smallest.
fn small_booleans(reps: usize) -> std::time::Duration {
    let mut spent = std::time::Duration::ZERO;
    for _ in 0..reps {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        m.rebuild_adjacency();
        let t = std::time::Instant::now();
        boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse");
        spent += t.elapsed();
    }
    spent
}

#[test]
#[ignore = "a measurement, not an assertion (run with --ignored --nocapture)"]
fn measure_boolean_wall_clock() {
    let threads = std::thread::available_parallelism().map_or(0, |n| n.get());
    let mode = if cfg!(feature = "parallel") {
        "parallel"
    } else {
        "serial"
    };
    println!("== nacre boolean wall clock ==  mode={mode}  cores={threads}");

    // Discard the first of each: rayon's pool and every worker's `HP_CONSTS` are built here.
    let _ = fin_fold(5);
    let _ = tooth_fold(5);
    let _ = small_booleans(20);

    for n in [7i128, 13, 25, 40, 60, 80] {
        let (r, a) = (fin_fold(n), tooth_fold(n));
        println!(
            "  fold {n:>2} : rotated {r:>10.2?}   axis-aligned {a:>10.2?}   certified-path tax {:.1}x",
            r.as_secs_f64() / a.as_secs_f64().max(1e-9)
        );
    }
    // The answer, so `tools/occt-bench.tcl`'s timings can be compared to these at all.
    let (_, faces, volume) = fin_fold_full(80);
    println!(
        "  fold 80 answer : {faces} faces, volume {volume:.6}   (OCCT, same fixture: 1974 faces, 237.212)"
    );
    // The order experiment, against `tools/occt-fins.tcl`'s OCCT numbers on the same shapes.
    for (label, hub_first, occt) in [
        ("fins only     ", false, "8.27s / n-ary 3.03s"),
        ("hub first+fins", true, "3.64s / n-ary 3.49s"),
    ] {
        let (d, faces, v) = fin_ring(80, hub_first);
        println!(
            "  ring 80 {label} : {d:>9.2?}   {faces} faces, vol {v:.4}   (OCCT 1 thread: {occt})"
        );
    }
    let reps = 200;
    let d = small_booleans(reps);
    println!(
        "  axis-aligned small   : {d:>12.2?} over {reps} booleans ({:.1?} each)",
        d / reps as u32
    );
}

/// Where `Profile2d::check`'s doc numbers come from. A convex ring is the worst case — nothing
/// short-circuits — and the predicates run on the rational truth (`orient2d_rat`,
/// narrow-first), so this is the price of exactness-on-the-truth over the old f64 signs.
#[test]
#[ignore]
fn measure_profile_check_wall_clock() {
    for n in [100usize, 1_000, 5_000] {
        let ring: Vec<Point2> = (0..n)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                p2(10.0 * a.cos(), 10.0 * a.sin())
            })
            .collect();
        let profile = Profile2d::polygon(ring).unwrap();
        let best = (0..3)
            .map(|_| {
                let t = std::time::Instant::now();
                profile.check().expect("a convex ring is clean");
                t.elapsed()
            })
            .min()
            .unwrap();
        println!("  check {n:>5} points : {best:>10.2?}");
    }
}
