//! The parallel build answers exactly as the sequential one (and the U-prism fixture it uses).

use super::*;

/// A signature that changes if `Store::push` order (hence handle identity) changes:
/// vertex points in handle order — exactly what `assemble_fuse_cut` assigns by first
/// appearance across `faces` — plus edge/face/solid counts and sorted volumes.
#[cfg(feature = "parallel")]
fn model_sig(m: &Model, solids: &[Handle<Solid>]) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    let _ = write!(s, "S{}", solids.len());
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let p = m.vertex_point(vh).as_array();
        let _ = write!(
            s,
            "|{:x},{:x},{:x}",
            p[0].to_bits(),
            p[1].to_bits(),
            p[2].to_bits()
        );
    }
    let _ = write!(s, "|E{}F{}", m.edge_count(), m.face_count());
    let mut vols: Vec<u64> = solids
        .iter()
        .map(|&sh| nacre_props::mass_props(m, sh).unwrap().volume.to_bits())
        .collect();
    vols.sort_unstable();
    let _ = write!(s, "|V{vols:?}");
    s
}

/// The parallel boolean must be bit-identical regardless of rayon thread count — replay
/// determinism (DNA) requires thread-order independence. Each fixture's result under a
/// 1-thread pool must equal the default many-thread result, run repeatedly so scheduling
/// jitter would show.
///
/// **★ Two things this test has to keep honest about itself.**
///
/// First, it passed for the whole stretch when there was *no* parallelism — the cutover
/// took the old engine's `par_iter` calls with it, and a one-thread pool is trivially
/// equal to a many-thread one when nothing forks. So the fixtures must have real
/// parallel width: the fin fold below reaches into the dozens of plane classes, where
/// the two-ngon fuse has barely a dozen.
///
/// Second, `model_sig` compares coordinates, counts and volumes — **not the report**.
/// `Notes::sorted` is a stable sort by site, so any two entries sharing a site keep
/// their *arrival* order, and under `parallel` that is the schedule's to decide. So the
/// signature includes `BoolReport`, which is the only thing that would catch it.
#[cfg(feature = "parallel")]
#[test]
fn parallel_boolean_is_thread_order_independent() {
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    // A hub with fins arrayed around it — every fin is turned by an angle with no exact
    // f64, so every judgement is on the toleranced path and the report is non-empty.
    let fin_fold = |n: i128| -> String {
        let mut m = replay(&[extrude_log_op(
            Profile2d::polygon(vec![
                p2(-3.0, -3.0),
                p2(3.0, -3.0),
                p2(3.0, 3.0),
                p2(-3.0, 3.0),
            ])
            .unwrap(),
            2.0,
        )])
        .unwrap();
        let mut acc = m.live_solids()[0];
        let mut sig = String::new();
        for i in 0..n {
            let fin = {
                let __w7 = SketchFrame::world(&m, Axis::Z);
                let out = ops::apply(
                    &mut m,
                    &Operation::Extrude {
                        frame: __w7,
                        profile: Profile2d::polygon(vec![
                            p2(2.0, -0.4),
                            p2(8.0, -0.4),
                            p2(8.0, 0.4),
                            p2(2.0, 0.4),
                        ])
                        .unwrap(),
                        dist: 1.0,
                    },
                )
                .unwrap();
                m.rebuild_adjacency();
                match out {
                    OpOutput::Extrude { solid, .. } => solid,
                    o => panic!("{o:?}"),
                }
            };
            let fin = transform(
                &mut m,
                fin,
                &Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::new(360 * i, n).unwrap()).unwrap(),
                }),
            )
            .unwrap();
            m.rebuild_adjacency();
            let (solids, report) = boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
            m.rebuild_adjacency();
            acc = solids[0];
            // **The report's decisions travel in the signature**, not just the geometry. Measured
            // on this fixture: up to 841 coincidences per boolean, so `loosest`'s `max_by_key` is
            // choosing among hundreds of candidates — which is the tie-break that a schedule could
            // otherwise decide.
            //
            // ★★★★★ **The count does not travel, because it is not an answer.**
            // `BoolReport::coincidences` is *"how many judgements were answered by a proved
            // coincidence"* — a measure of **work done**, and how much work the parallel phases do
            // is exactly what a schedule decides. The two decision fields carry their own rule
            // (`loosest`: *"Ties keep the first in sorted order, so the answer does not depend on
            // the schedule"*); the count carries none and never could.
            //
            // Measured with the tests optimized, which changes the interleaving: the same fold
            // reported 755 coincidences on one schedule and 735 on another while the 2,678
            // characters of decisions and geometry around it were **identical**. Comparing it
            // asserted something the kernel does not claim. What it is here for — that the report
            // is not empty, so comparing it proves something — is the assertion below, kept.
            assert!(
                i == 0 || report.coincidences > 0,
                "the report is empty, so comparing it proves nothing"
            );
            sig.push_str(&format!("{:?}{:?}", report.merges, report.loosest));
        }
        sig.push_str(&model_sig(&m, &[acc]));
        sig
    };

    let build = || {
        let ngon = |n: usize, r: f64, cx: f64, cy: f64| {
            Profile2d::polygon(
                (0..n)
                    .map(|i| {
                        let ang = std::f64::consts::TAU * (i as f64) / (n as f64);
                        p2(cx + r * ang.cos(), cy + r * ang.sin())
                    })
                    .collect(),
            )
            .unwrap()
        };
        let mut m = replay(&[
            extrude_log_op(ngon(16, 2.0, 0.0, 0.0), 3.0),
            extrude_log_op(ngon(16, 2.0, 2.5, 0.5), 3.0),
        ])
        .unwrap();
        let a = m.live_solids()[0];
        let b = m.live_solids()[1];
        let up = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
        let b = transform(&mut m, b, &up).unwrap();
        m.rebuild_adjacency();
        let tilt = rot_iso(Axis::X, 30);
        let a = transform(&mut m, a, &tilt).unwrap();
        m.rebuild_adjacency();
        let b = transform(&mut m, b, &tilt).unwrap();
        m.rebuild_adjacency();
        (m, a, b)
    };
    let two_ngons = || {
        let (mut m, a, b) = build();
        let (solids, report) = boolean_with_report(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        // The same rule as the fold above: decisions travel, the work count does not.
        format!(
            "{:?}{:?}{}",
            report.merges,
            report.loosest,
            model_sig(&m, &solids)
        )
    };
    let pool1 = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let seven_fins = || fin_fold(7);
    // **The only fixture here whose alias table is not empty.** `Aliases::record` returns
    // immediately below four planes, so every other model leaves the table at zero and
    // never exercises the snapshot-and-absorb a parallel round is built on. Measured on
    // this one: 36 aliases, settled over two rounds.
    let four_plane = || {
        let (mut m, target, bar) = four_plane_model(0.2, 45);
        let (solids, report) = boolean_with_report(&mut m, BoolKind::Cut, target, bar).unwrap();
        m.rebuild_adjacency();
        // The same rule as the fold above: decisions travel, the work count does not.
        format!(
            "{:?}{:?}{}",
            report.merges,
            report.loosest,
            model_sig(&m, &solids)
        )
    };
    for (name, run) in [
        (
            "two rotated ngons",
            &two_ngons as &(dyn Fn() -> String + Sync),
        ),
        (
            "a seven-fin fold",
            &seven_fins as &(dyn Fn() -> String + Sync),
        ),
        (
            "a four-plane concurrency",
            &four_plane as &(dyn Fn() -> String + Sync),
        ),
    ] {
        let reference = pool1.install(run);
        // A signature that came back empty would make this vacuous whatever the schedule.
        assert!(!reference.is_empty(), "{name}: nothing to compare");
        for _ in 0..4 {
            assert_eq!(
                run(),
                reference,
                "{name}: the result depends on thread order"
            );
        }
    }
}

#[test]
fn u_prism_is_valid() {
    // Pin the fixture itself: a mistyped profile could still trip `multichord`
    // below, for the wrong reason.
    let (m, u) = u_prism();
    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");
    let vol = nacre_props::mass_props(&m, u).unwrap().volume;
    assert!((vol - 5.3).abs() < 1e-9, "volume {vol}");
}
