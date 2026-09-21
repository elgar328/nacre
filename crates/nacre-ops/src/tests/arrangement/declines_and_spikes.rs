//! What the trace declines, and the measurement spikes (where the boolean spends its time).

use super::*;

/// Axis-aligned cubes never decline: every face is seated, a clean transversal chord, or a
/// parallel miss. This is falsifiable — a `declined` entry would mean the producer hit a
/// degeneracy it should not on this input.
#[test]
fn axis_aligned_cubes_decline_nothing() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    for wc in 0..planes.len() {
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        assert!(tr.declined.is_empty(), "class {wc} declined: {tr:?}");
    }
}

fn tilt_by(m: &mut Model, s: Handle<Solid>, deg: nacre_exact::Rat) -> Handle<Solid> {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let out = transform(
        m,
        s,
        &Isometry::rotation(Rotation {
            axis: Axis::Z,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(deg).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    out
}

/// **The crossing collector's direction families really are equivalence classes.**
///
/// ★ `split_at_crossings` replaced a predicate per wall pair with "same family?", which is sound
/// only because `wc ∩ w ∥ wc ∩ r` is transitive. The argument is in that function; this is the
/// **check**, over every wall pair of every class of a rotated and an axis-aligned fold: the
/// predicate's own answer must agree with the partition, everywhere.
///
/// It matters that this is a test and not a spike. A violation would not reject or panic — it
/// would invent a crossing point where the lines never meet, or lose one where they do, and the
/// census would drift by one vertex somewhere. **A silent wrong answer is exactly what a corpus
/// this size can hide**, so the invariant is asserted rather than inspected.
///
/// Separate from the timing spike on purpose: this asks the predicate about every pair, which
/// fills `WitnessPoint`'s realization cells and would make the phase timers read the collector 1.76x
/// cheaper than it is.
#[test]
#[ignore = "slow: every wall pair of every class, two folds"]
fn direction_families_partition_the_walls() {
    for rotated in [true, false] {
        let n = 24i128;
        let mut m = Model::new();
        let mut acc = if rotated {
            m.add_cuboid(
                Point3::from_array([-3.0, -3.0, 0.0]),
                Point3::from_array([3.0, 3.0, 2.0]),
            )
        } else {
            m.add_cuboid(
                Point3::from_array([-1.0, -1.0, 0.0]),
                Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
            )
        };
        m.rebuild_adjacency();
        let mut pairs = 0usize;
        for i in 0..n {
            let fin = if rotated {
                let f = m.add_cuboid(
                    Point3::from_array([2.0, -0.4, 0.0]),
                    Point3::from_array([8.0, 0.4, 1.0]),
                );
                m.rebuild_adjacency();
                tilt_by(&mut m, f, nacre_exact::Rat::new(360 * i, n).unwrap())
            } else {
                let x = i as f64 * 0.5;
                let f = m.add_cuboid(
                    Point3::from_array([x, 0.5, 0.0]),
                    Point3::from_array([x + 0.2, 4.0, 2.0]),
                );
                m.rebuild_adjacency();
                f
            };
            // The audit needs the *same* judging context the boolean uses, and the walls of each
            // class as the tracer finds them — so it re-runs the setup and the trace, then checks
            // the partition the collector would build.
            let setup = plane_index_setup(&m, acc, fin).expect("setup");
            let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
            let trace_in = combinatorics::trace_input(
                &m,
                [(acc, &setup.inc_a), (fin, &setup.inc_b)],
                &setup.surf_ix,
                setup.planes.len(),
                &jd,
                &setup.plane_ix,
                &setup.cyls,
                Default::default(),
            );
            for wc in 0..setup.geom.len() {
                let tr = trace_on_class(
                    &trace_in,
                    wc,
                    &jd,
                    &setup.cyls,
                    &setup.planes,
                    &setup.plane_ix,
                    &Aliases::default(),
                );
                let mut walls: Vec<usize> = Vec::new();
                for s in &tr.segs {
                    if !walls.contains(&s.wall) {
                        walls.push(s.wall);
                    }
                }
                let par = |a: usize, b: usize| jd.plane_pair_dir_sign(wc, a, b) == 0;
                // The partition, exactly as `split_at_crossings` builds it.
                let mut reps: Vec<usize> = Vec::new();
                let dir: Vec<usize> = walls
                    .iter()
                    .map(|&w| {
                        reps.iter().position(|&rep| par(rep, w)).unwrap_or_else(|| {
                            reps.push(w);
                            reps.len() - 1
                        })
                    })
                    .collect();
                for (i, &a) in walls.iter().enumerate() {
                    for (j, &b) in walls.iter().enumerate().skip(i + 1) {
                        pairs += 1;
                        assert_eq!(
                            par(a, b),
                            dir[i] == dir[j],
                            "class {wc}: walls {a} and {b} — the predicate and the partition \
                                 disagree, so parallelism is not transitive here"
                        );
                    }
                }
            }
            acc = super::super::boolean(&mut m, BoolKind::Fuse, acc, fin)
                .expect("fuse")
                .0[0];
            m.rebuild_adjacency();
        }
        println!(
            "  {} fold: {pairs} wall pairs audited, all agree",
            if rotated { "rotated" } else { "axis-aligned" }
        );
        assert!(
            pairs > 10_000,
            "{} fold audited only {pairs} pairs — too few to mean anything",
            if rotated { "rotated" } else { "axis-aligned" }
        );
    }
}

/// **Where does a whole boolean's time go?** Every earlier profile answered a share *of a
/// phase* — and the one that mattered was never taken.
///
/// ★ Two ways to be wrong that this is built against:
///
/// 1. **Mixed clocks.** The last breakdown summed per-thread timers (CPU) and read the result
///    against the fold's wall clock. The part came out larger than the whole — 3,719ms of a
///    2,540ms fold — and the ratio it implied was meaningless. **Run this with
///    `--no-default-features`**, where every phase and the total are the same clock.
/// 2. **Phases that do not add up.** The sum is printed against the measured whole, so anything
///    unaccounted for shows as a gap rather than hiding inside a phase's share.
/// 3. ★★★ **Another test adding to the same counters.** `phase::` are process-global atomics
///    and `cargo test` runs this binary's tests on parallel threads, so `reset()` here does not
///    fence anything: every concurrently running test that calls `boolean` — and in an
///    `--ignored` pass that is the *other* spikes and the two rotation stress tests — lands in
///    the same buckets, while `whole` below is this thread's wall clock alone. The share
///    percentages then exceed 100 for reasons that have nothing to do with the code.
///    **Run it alone**: `cargo test -p nacre-ops --no-default-features spike_where -- --ignored
///    --nocapture --test-threads=1`.
///
/// Timers live on the production path (`crate::phase`, charged from the stages themselves), not
/// in a replica of it.
#[test]
#[ignore = "spike"]
fn measure_spike_where_the_boolean_spends_it() {
    for rotated in [true, false] {
        spend(60, rotated);
    }
}

/// One fold's phase breakdown. `rotated` picks which predicate routes the judgements take:
/// an axis-aligned fold answers on the exact path, a rotated one mostly on the certified one,
/// and the difference between the two breakdowns is what the certification actually costs.
fn spend(n: i128, rotated: bool) {
    let mut m = Model::new();
    // ★ **The two folds must be the same *shape* of work, not the same model minus a rotation.**
    // Dropping the tilt stacks all 60 fins on one another — 7.6x fewer segments, a degenerate
    // model, and a comparison that says nothing. The axis-aligned arm places each fin at its
    // own x instead, so both arms fuse `n` distinct blades onto a growing solid and the only
    // difference is which predicate route the judgements take.
    let mut acc = if rotated {
        m.add_cuboid(
            Point3::from_array([-3.0, -3.0, 0.0]),
            Point3::from_array([3.0, 3.0, 2.0]),
        )
    } else {
        m.add_cuboid(
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([n as f64 * 0.5 + 1.0, 1.0, 3.0]),
        )
    };
    m.rebuild_adjacency();
    let blade = |m: &mut Model, i: i128| {
        let f = if rotated {
            m.add_cuboid(
                Point3::from_array([2.0, -0.4, 0.0]),
                Point3::from_array([8.0, 0.4, 1.0]),
            )
        } else {
            let x = i as f64 * 0.5;
            m.add_cuboid(
                Point3::from_array([x, 0.5, 0.0]),
                Point3::from_array([x + 0.2, 4.0, 2.0]),
            )
        };
        m.rebuild_adjacency();
        if rotated {
            tilt_by(m, f, nacre_exact::Rat::new(360 * i, n).unwrap())
        } else {
            f
        }
    };
    // Warm the code paths, then zero the counters: the first boolean pays for lazily-built
    // caches that the other n do not.
    {
        let f = blade(&mut m, if rotated { 1 } else { n });
        acc = super::super::boolean(&mut m, BoolKind::Fuse, acc, f)
            .expect("warm")
            .0[0];
        m.rebuild_adjacency();
    }
    phase::reset();
    phase::scale::reset();

    let mut whole = std::time::Duration::ZERO;
    for i in 0..n {
        let fin = blade(&mut m, i);
        let t = std::time::Instant::now();
        acc = super::super::boolean(&mut m, BoolKind::Fuse, acc, fin)
            .expect("fuse")
            .0[0];
        whole += t.elapsed();
        m.rebuild_adjacency();
    }

    let rows = phase::all();
    // ★ Indentation is nesting: depth 0 are `boolean`'s own phases, depth 1 the inside of
    // `trace_result_faces`, depth 2 the inside of `split_at_crossings`. Depth 1 partitions the
    // work depth 0 does not name, so 0 and 1 add to the whole — and adding depth 2 on top would
    // count `split_at_crossings` twice.
    let depth = |l: &str| (l.len() - l.trim_start().len()) / 2;
    let sum: u64 = rows
        .iter()
        .filter(|(l, _)| depth(l) <= 1)
        .map(|(_, ns)| ns)
        .sum();
    let total = whole.as_nanos() as u64;
    println!(
        "\n{n}-fin fold ({}), whole boolean, serial build:",
        if rotated { "rotated" } else { "axis-aligned" }
    );
    for (label, ns) in &rows {
        println!(
            "  {label:<32} {:>8.1?}  {:>5.1}%",
            std::time::Duration::from_nanos(*ns),
            100.0 * *ns as f64 / total as f64
        );
    }
    println!(
        "  {:<32} {:>8.1?}  {:>5.1}%   ← accounted",
        "sum of the above",
        std::time::Duration::from_nanos(sum),
        100.0 * sum as f64 / total as f64
    );
    println!("  {:<32} {whole:>8.1?}  100.0%   ← measured", "the fold");
    println!(
        "\n  ★ split_at_crossings is {:.1}% of the whole boolean. The plan continues at 25%.",
        100.0
            * rows
                .iter()
                .find(|(l, _)| l.trim() == "split_at_crossings")
                .map(|(_, ns)| *ns)
                .unwrap_or(0) as f64
            / total as f64
    );

    // ★ **Scale, which decides whether the structural fixes are worth building.** A hull test
    // over wall *pairs* replaces `2 × |segs on r|` predicate calls with `2` — worth nothing when
    // a wall carries one segment, and worth the ratio when it carries many.
    use phase::scale as sc;
    let (segs, walls) = (sc::get(&sc::SEGS), sc::get(&sc::WALLS));
    println!("\n  scale, summed over every class of every boolean:");
    println!(
        "    segments / walls        {segs:>10} / {walls:<10} = {:.2}   ← S3b lives or dies here",
        segs as f64 / walls.max(1) as f64
    );
    println!("    split points (reps)     {:>10}", sc::get(&sc::PTS));
    println!(
        "    (1) collect trips       {:>10}",
        sc::get(&sc::COLLECT_TRIPS)
    );
    println!(
        "    (3) cover trips         {:>10}   = {:.2}x the collecting loop",
        sc::get(&sc::COVER_TRIPS),
        sc::get(&sc::COVER_TRIPS) as f64 / sc::get(&sc::COLLECT_TRIPS).max(1) as f64
    );
}

/// The loop-invariant hoist's gate: **what did hoisting it do to the evidence?** `plane_pair_dir_sign`
/// records on the rotated path, and `BoolReport::coincidences` is a *count*, so asking the same
/// question fewer times moves it. The answer must not move; the count may.
#[test]
#[ignore = "spike"]
fn measure_spike_report_after_hoisting() {
    let n = 24i128;
    let mut m = Model::new();
    let mut acc = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let (mut total, mut loosest) = (0usize, String::new());
    for i in 0..n {
        let fin = m.add_cuboid(
            Point3::from_array([2.0, -0.4, 0.0]),
            Point3::from_array([8.0, 0.4, 1.0]),
        );
        m.rebuild_adjacency();
        let fin = tilt_by(&mut m, fin, nacre_exact::Rat::new(360 * i, n).unwrap());
        let (solids, report) =
            crate::boolean_with_report(&mut m, BoolKind::Fuse, acc, fin).unwrap();
        total += report.coincidences;
        if let Some(e) = &report.loosest {
            loosest = format!("{e:?}");
        }
        acc = solids[0];
        m.rebuild_adjacency();
    }
    let v = nacre_props::mass_props(&m, acc).unwrap().volume;
    println!("\n  coincidences over {n} booleans : {total}");
    println!("  loosest (last)                : {loosest}");
    println!("  volume                        : {v:.9}");
}

/// What is actually left to save, **in production's configuration**?
///
/// The earlier phase timing used `ClassReuse::Off`, so it counted work production never does —
/// `reuse.rs` skips most classes outright. This counts what survives that, how much of it a
/// per-class bounding box would cull, and how often the two things that would make the cull
/// unsound actually occur. Counting only; nothing is built and nothing is timed (the machine is
/// not quiet).
#[test]
#[ignore = "spike"]
fn measure_spike_cull_potential() {
    let n = 60i128;
    let mut m = Model::new();
    let mut acc = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    let mut tot = Counts::default();
    for i in 0..n {
        let fin = m.add_cuboid(
            Point3::from_array([2.0, -0.4, 0.0]),
            Point3::from_array([8.0, 0.4, 1.0]),
        );
        m.rebuild_adjacency();
        let fin = tilt_by(&mut m, fin, nacre_exact::Rat::new(360 * i, n).unwrap());
        tot.add(&count_cull(&m, acc, fin));
        acc = super::super::boolean(&mut m, BoolKind::Fuse, acc, fin)
            .expect("fuse")
            .0[0];
        m.rebuild_adjacency();
    }
    println!("over the fold, production config (ClassReuse::Proved):");
    println!("  classes                 {:>9}", tot.classes);
    println!(
        "  arranged after reuse    {:>9}  ({:.1}%)",
        tot.arranged,
        100.0 * tot.arranged as f64 / tot.classes as f64
    );
    println!(
        "  (face, class) pairs     {:>9}   in arranged classes",
        tot.pairs
    );
    println!(
        "  ★ cullable by box       {:>9}  ({:.1}%)",
        tot.cullable,
        100.0 * tot.cullable as f64 / tot.pairs.max(1) as f64
    );
    println!(
        "  seated footprint is one piece: {} of {} arranged classes",
        tot.one_piece, tot.arranged
    );
    // ★ Recorded, and **not** the operative number. The plan restricted the cull to classes
    // whose seated footprint is one connected piece, fearing that several pieces would need
    // several seeds. They do not: what the correction needs is that the culled chords' parity
    // be *constant over the footprint*, and the cull criterion (`box(F)` disjoint from the
    // footprint box) already guarantees no culled chord enters that box. A box is connected, so
    // the parity is one constant however many pieces the seated faces form.
    println!(
        "  (one-piece classes only:  {:>6}  of {} — not the limit, see the note)",
        tot.cullable_1p, tot.pairs_1p
    );
    println!(
        "  (class, solid) needing a seed: {}  of {}   — of those, {} are free (solid box misses)",
        tot.needs_seed,
        2 * tot.arranged,
        tot.seed_free
    );
}

fn count_cull(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> Counts {
    let mut c = Counts::default();
    let Ok(setup) = plane_index_setup(m, a, b) else {
        return c;
    };
    let PlaneSetup {
        planes: faces_tab,
        n_a,
        geom,
        plane_ix,
        class_owner,
        standard,
        notes,
        ..
    } = setup;
    let jd = Judge::new(&geom, standard, &notes);
    let _ = &jd;
    let plans = crate::reuse::class_plans(
        m,
        crate::reuse::ClassReuse::Proved,
        BoolKind::Fuse,
        a,
        b,
        &geom,
        &class_owner,
    );
    let boxes: Vec<[[f64; 2]; 3]> = faces_tab
        .iter()
        .map(|f| face_box(m, f.face().expect("real")))
        .collect();
    let overlap = |x: &[[f64; 2]; 3], y: &[[f64; 2]; 3]| {
        (0..3).all(|k| x[k][0] <= y[k][1] && y[k][0] <= x[k][1])
    };
    c.classes = geom.len();
    for (wc, plan) in plans.iter().enumerate() {
        if *plan != crate::reuse::ClassPlan::Arrange {
            continue;
        }
        c.arranged += 1;
        let seated: Vec<usize> = (0..faces_tab.len())
            .filter(|&f| plane_ix[f].plane() == wc)
            .collect();
        if seated.is_empty() {
            continue;
        }
        // The class's region of interest, and whether it is one connected piece (box graph).
        let mut foot = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
        for &s in &seated {
            for k in 0..3 {
                foot[k][0] = foot[k][0].min(boxes[s][k][0]);
                foot[k][1] = foot[k][1].max(boxes[s][k][1]);
            }
        }
        let mut parent: Vec<usize> = (0..seated.len()).collect();
        for i in 0..seated.len() {
            for j in (i + 1)..seated.len() {
                if overlap(&boxes[seated[i]], &boxes[seated[j]]) {
                    let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                    parent[ri] = rj;
                }
            }
        }
        let pieces = (0..seated.len())
            .map(|i| uf_find(&mut parent, i))
            .collect::<std::collections::HashSet<_>>()
            .len();
        if pieces == 1 {
            c.one_piece += 1;
        }
        // How many of this class's traced faces the box would cull, and whether a solid that
        // only crosses `W` has any culled face (that is the one needing a seed).
        let mut culled_side = [false; 2];
        let mut seated_side = [false; 2];
        for &s in &seated {
            seated_side[usize::from(s >= n_a)] = true;
        }
        for f in 0..faces_tab.len() {
            c.pairs += 1;
            if pieces == 1 {
                c.pairs_1p += 1;
            }
            if plane_ix[f].plane() == wc {
                continue; // seated: never culled
            }
            if !overlap(&boxes[f], &foot) {
                c.cullable += 1;
                if pieces == 1 {
                    c.cullable_1p += 1;
                }
                culled_side[usize::from(f >= n_a)] = true;
            }
        }
        // Each operand's whole box: if it misses the footprint entirely it cannot enclose it,
        // so its parity is false for free.
        let mut solid_box = [[[f64::INFINITY, f64::NEG_INFINITY]; 3]; 2];
        for (f, fb) in boxes.iter().enumerate() {
            let side = usize::from(f >= n_a);
            for k in 0..3 {
                solid_box[side][k][0] = solid_box[side][k][0].min(fb[k][0]);
                solid_box[side][k][1] = solid_box[side][k][1].max(fb[k][1]);
            }
        }
        for side in 0..2 {
            if culled_side[side] && !seated_side[side] {
                c.needs_seed += 1;
                if !overlap(&solid_box[side], &foot) {
                    c.seed_free += 1;
                }
            }
        }
    }
    c
}
