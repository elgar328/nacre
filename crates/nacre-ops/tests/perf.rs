//! **Where the boolean's wall clock is measured, so that speed claims come from a number.**
//!
//! `#[ignore]`d, and it asserts nothing about time: a timing assertion is flaky on a shared
//! machine and would eventually be relaxed into meaninglessness. The values belong in the
//! commit message of whatever change moved them, next to the census diff that says the
//! answers did not change.
//!
//! Two workloads, because parallelism pulls them in opposite directions:
//!
//! - **Rotated fold** — a hub with `n` fins arrayed around it. Rotation puts every judgement
//!   on the toleranced CIP path, plane classes grow into the dozens, and per-class work is
//!   milliseconds. This is where more cores help.
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
use nacre_ops::{BoolKind, Operation, Profile2d, SketchPlane, apply, boolean};
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
            plane: SketchPlane::world_xy(),
            profile: Profile2d::polygon(ring),
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

fn turn(m: &mut Model, s: Handle<Solid>, deg: Rat) -> Handle<Solid> {
    let out = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
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
    spent
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
fn boolean_wall_clock() {
    let threads = std::thread::available_parallelism().map_or(0, |n| n.get());
    let mode = if cfg!(feature = "parallel") {
        "parallel"
    } else {
        "serial"
    };
    println!("== nacre boolean wall clock ==  mode={mode}  cores={threads}");

    // Discard the first of each: rayon's pool and every worker's `HP_CONSTS` are built here.
    let _ = fin_fold(5);
    let _ = small_booleans(20);

    for n in [7i128, 13, 25] {
        let d = fin_fold(n);
        println!("  rotated fold, {n:>2} fins : {d:>12.2?}");
    }
    let reps = 200;
    let d = small_booleans(reps);
    println!(
        "  axis-aligned small   : {d:>12.2?} over {reps} booleans ({:.1?} each)",
        d / reps as u32
    );
}
