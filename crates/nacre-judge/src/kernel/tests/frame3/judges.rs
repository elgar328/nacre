//! The judges: orient3d, indirect orient3d, indirect cmp_coord, dir_sign.

use super::*;

/// A point rotated about one random axis (through a random pivot) by one **inexact**
/// rational angle — generic non-degenerate, heterogeneous provenance (each point its
/// own rotation), `det3_bound_soundness`'s corpus shape. Inexact angles only (an exact/near-coplanar
/// mix would expose the tol-0 GT noise floor — that is the declare-0 test's job).
fn rand_point(st: &mut u64) -> WitnessPoint {
    let base = rand_base(st);
    let axis = axis_of(rng(st, 0, 2));
    let ang = deg(rng(st, 0, 360_000), rng(st, 1, 9973)); // inexact
    let pivot = if rng(st, 0, 1) == 0 {
        [Rat::from_int(0); 3]
    } else {
        rand_base(st) // non-origin pivot: exercises the pivot tol → det3_bound path
    };
    WitnessPoint::at(base).rotate_about(axis, ang, pivot)
}

/// `det3_bound` soundness: over many random heterogeneous-rotation 4-point configs (inexact
/// angles, arbitrary pivots), the bound must upper-bound the real error of the f64 determinant vs
/// the astro-float truth — never exceeded. `#[ignore]`: slow (4× astro-float per sample).
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn det3_bound_soundness() {
    const GT: usize = 512;
    let mut st = 0x3D00_1234_ABCD_EF01u64;
    let (mut bad, mut worst) = (0usize, 0.0_f64);
    for _ in 0..5000 {
        let (pa, pb, pc, pd) = (
            rand_point(&mut st),
            rand_point(&mut st),
            rand_point(&mut st),
            rand_point(&mut st),
        );
        let det = det3_f64(pa.coord(), pb.coord(), pc.coord(), pd.coord());
        let truth = det3_hp(&pa, &pb, &pc, &pd, GT);
        let err = abs_err(det, &truth, GT);
        let bound = det3_bound(
            [pa.coord(), pb.coord(), pc.coord(), pd.coord()],
            [pa.tol(), pb.tol(), pc.tol(), pd.tol()],
        );
        if err > bound {
            bad += 1;
        }
        if bound > 0.0 {
            worst = worst.max(err / bound);
        }
    }
    eprintln!("[det3_bound N=5000] violations: {bad}; worst tightness: {worst:.3}");
    assert_eq!(
        bad, 0,
        "det3_bound must bound the determinant error on every sample"
    );
}

/// Fast port check (default suite): a handful of random configs must not violate
/// `det3_bound` (catches a transcription bug; the full statistical soundness is the
/// `#[ignore]`d `det3_bound_soundness`).
#[test]
fn det3_bound_port_check() {
    const GT: usize = 512;
    let mut st = 0x0BAD_F00D_1234_5678u64;
    for _ in 0..40 {
        let (pa, pb, pc, pd) = (
            rand_point(&mut st),
            rand_point(&mut st),
            rand_point(&mut st),
            rand_point(&mut st),
        );
        let det = det3_f64(pa.coord(), pb.coord(), pc.coord(), pd.coord());
        let err = abs_err(det, &det3_hp(&pa, &pb, &pc, &pd, GT), GT);
        let bound = det3_bound(
            [pa.coord(), pb.coord(), pc.coord(), pd.coord()],
            [pa.tol(), pb.tol(), pc.tol(), pd.tol()],
        );
        assert!(err <= bound, "det3_bound {bound:e} < err {err:e}");
    }
}

/// A known tetrahedron judges `Positive`, its mirror `Negative`, and the sign is
/// preserved under a shared rotation (rotating all four points does not flip
/// orientation). Exact (tol 0) points take the filter's fast path.
#[test]
fn orient3d_sanity_and_rotation_invariance() {
    let pt = |x, y, z| WitnessPoint::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
    // d at origin; a,b,c along +x,+y,+z → right-handed → Positive.
    let (a, b, c, d) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1), pt(0, 0, 0));
    assert_eq!(
        orient3d_judge(&a, &b, &c, &d, fixture()).orient(),
        Orient::Positive
    );
    assert_eq!(
        orient3d_judge(&b, &a, &c, &d, fixture()).orient(),
        Orient::Negative,
        "swap → mirror"
    );

    // Shared rotation of all four (37° about Z through a rational pivot) keeps the sign.
    let rot = |p: &WitnessPoint| {
        p.clone()
            .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-3, 1), ri(0, 1)])
    };
    assert_eq!(
        orient3d_judge(&rot(&a), &rot(&b), &rot(&c), &rot(&d), fixture()).orient(),
        Orient::Positive,
        "orientation is rotation-invariant"
    );
}

/// Declare-0: four **exactly coplanar** rational points (tol 0) → the determinant is
/// exactly 0 → `Orient::Zero` (the ask-the-user case).
#[test]
fn orient3d_coplanar_is_zero() {
    let pt = |x, y, z| WitnessPoint::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
    // all four in the plane z = 0.
    let (a, b, c, d) = (pt(0, 0, 0), pt(3, 0, 0), pt(0, 5, 0), pt(2, 7, 0));
    assert_eq!(
        orient3d_judge(&a, &b, &c, &d, fixture()).orient(),
        Orient::Zero
    );
}

/// `Bounded` arithmetic keeps a sound radius (worst-case interval), and `sign` is
/// definite exactly when the interval clears 0.
#[test]
fn iv_arithmetic_is_sound_and_sign_decides() {
    let (a, b) = (Bounded::new(3.0, 0.1), Bounded::new(-2.0, 0.2));
    // sub radius ≥ sum of input radii.
    assert!(a.sub(b).error >= 0.1 + 0.2);
    // mul radius ≥ |mid_a|·rad_b + |mid_b|·rad_a (+ error·error + rounding).
    assert!(a.mul(b).error >= 3.0 * 0.2 + 2.0 * 0.1);
    assert_eq!(Bounded::new(1.0, 0.5).sign(), Some(true));
    assert_eq!(Bounded::new(-1.0, 0.5).sign(), Some(false));
    assert_eq!(Bounded::new(0.3, 0.5).sign(), None); // straddles 0 → escalate
}

/// Sanity: three coordinate planes meet at the origin `V=(0,0,0)`; `orient3d(V,q,r,s)`
/// judges a known sign, flips on a q/r swap, is `Zero` when `s` is coplanar with
/// `V,q,r`, and is invariant under a shared rotation (rotations preserve orientation).
#[test]
fn indirect_sanity_and_rotation_invariance() {
    fn tri(
        t: &(WitnessPoint, WitnessPoint, WitnessPoint),
    ) -> (&WitnessPoint, &WitnessPoint, &WitnessPoint) {
        (&t.0, &t.1, &t.2)
    }
    let pt = |x, y, z| WitnessPoint::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
    let pa = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 1, 0)); // z = 0
    let pb = (pt(0, 0, 0), pt(1, 0, 0), pt(0, 0, 1)); // y = 0
    let pc = (pt(0, 0, 0), pt(0, 1, 0), pt(0, 0, 1)); // x = 0  → V = (0,0,0)
    let (q, r, s) = (pt(1, 0, 0), pt(0, 1, 0), pt(0, 0, 1));
    assert_eq!(
        indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s, fixture()).orient(),
        Orient::Negative
    );
    assert_eq!(
        indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &r, &q, &s, fixture()).orient(),
        Orient::Positive,
        "q/r swap flips the sign"
    );
    // s coplanar with V,q,r (all z = 0) → orient exactly 0 → declare-0.
    let s0 = pt(1, 1, 0);
    assert_eq!(
        indirect_orient3d_judge(tri(&pa), tri(&pb), tri(&pc), &q, &r, &s0, fixture()).orient(),
        Orient::Zero
    );
    // Shared rotation of all twelve points keeps the definite sign.
    let rot = |p: &WitnessPoint| {
        p.clone()
            .rotate_about(Axis::Z, deg(37, 1), [ri(2, 1), ri(-1, 1), ri(0, 1)])
    };
    let rp = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| (rot(t.0), rot(t.1), rot(t.2));
    let (ra, rb, rc) = (rp(tri(&pa)), rp(tri(&pb)), rp(tri(&pc)));
    let (rq, rr, rs) = (rot(&q), rot(&r), rot(&s));
    assert_eq!(
        indirect_orient3d_judge(tri(&ra), tri(&rb), tri(&rc), &rq, &rr, &rs, fixture()).orient(),
        Orient::Negative,
        "indirect orient is rotation-invariant"
    );
}

/// Rotated plane coefficient tol soundness. A plane's four coefficients
/// derive from three rotated points (subtraction/cross/dot); the interval `error` on
/// each must upper-bound the real f64 error vs the astro-float truth. Corpus mixes
/// origin and arbitrary pivots (`rand_point`). `#[ignore]`: slow astro-float GT.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn plane_coefficient_tol_soundness() {
    const GT: usize = 512;
    const N: usize = 10_000;
    let mut st = 0xB0B0_5555_1111_2222u64;
    let (mut bad, mut worst) = (0usize, 0.0_f64);
    for _ in 0..N {
        let (p0, p1, p2) = (
            rand_point(&mut st),
            rand_point(&mut st),
            rand_point(&mut st),
        );
        let iv = plane_iv(&p0, &p1, &p2);
        let hp = plane_hp(&p0, &p1, &p2, GT);
        for k in 0..4 {
            let err = abs_err(iv[k].value, &hp[k], GT);
            if err > iv[k].error {
                bad += 1;
            }
            if iv[k].error > 0.0 {
                worst = worst.max(err / iv[k].error);
            }
        }
    }
    eprintln!("[plane_coefficient_tol N={N}] violations: {bad}; worst tightness: {worst:.3}");
    assert_eq!(bad, 0, "plane coefficient tol must bound the error");
}

/// **Asking one question four ways must give one answer.**
///
/// Four planes through a common point are concurrent or they are not, and that fact does not
/// depend on which three of them you call "the point" and which one you call "the query". The
/// judge is free to *abstain* — but not to say `Positive` for one splitting and `Zero` for
/// another, because then two callers reading the same geometry disagree, and the arrangement
/// enters one vertex twice.
///
/// **The corpus is built, not sampled.** Random points never place a query vertex *on* the
/// implicit point, and that coincidence is the whole difficulty: it makes `row1 = D·(V − s)`
/// cancel to nothing. Here every plane is spanned by the shared point `p` and two others, so
/// `p` is both the meet of any three and a defining point of the fourth — the exact shape the
/// engine hits when a tool edge lands in a target plane.
#[test]
fn one_concurrency_read_four_ways_gives_one_answer() {
    let mut st = 0x5EED_1234_ABCD_9999u64;
    let mut checked = 0usize;
    for _ in 0..200 {
        // A rotated common point, and four planes each spanned by it and two more.
        let axis = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 1, 359_000), rng(&mut st, 1, 997)); // inexact ⇒ toleranced
        let pivot = rand_base(&mut st);
        let turn = |b: [Rat; 3], st: &mut u64| {
            let _ = st;
            WitnessPoint::at(b).rotate_about(axis, ang, pivot)
        };
        let p = turn(rand_base(&mut st), &mut st);
        let spans: Vec<[WitnessPoint; 2]> = (0..4)
            .map(|_| {
                [
                    turn(rand_base(&mut st), &mut st),
                    turn(rand_base(&mut st), &mut st),
                ]
            })
            .collect();
        let plane = |i: usize| (&p, &spans[i][0], &spans[i][1]);

        // Every way of choosing which three planes make the point and which one is queried.
        let mut verdicts = Vec::new();
        for q in 0..4 {
            let tri: Vec<usize> = (0..4).filter(|&i| i != q).collect();
            // Skip a splitting whose three planes do not meet at a single point at all.
            if dir_sign_judge(plane(tri[0]), plane(tri[1]), plane(tri[2]), fixture()).orient()
                == Orient::Zero
            {
                continue;
            }
            // `p` goes in the **third** slot: `indirect_hp` forms `row1 = Dvec − D·s` from
            // that one, so putting the shared point there is what makes the subtraction cancel
            // — the configuration the engine actually hit. With `p` first the term stays
            // healthy and the corpus measures nothing.
            let t = plane(q);
            verdicts.push(
                indirect_orient3d_judge(
                    plane(tri[0]),
                    plane(tri[1]),
                    plane(tri[2]),
                    t.1,
                    t.2,
                    t.0,
                    fixture(),
                )
                .orient(),
            );
        }
        if verdicts.len() < 2 {
            continue; // nothing to compare
        }
        checked += 1;
        assert!(
            verdicts.iter().all(|v| *v == verdicts[0]),
            "the same concurrency read four ways: {verdicts:?}"
        );
    }
    // A sweep that quietly stops constructing anything reads as agreement.
    assert!(
        checked > 20,
        "only {checked} configurations were comparable"
    );
    eprintln!("[form-invariance] {checked} concurrent configurations, all splittings agreed");
}

/// Indirect orient3d soundness over heterogeneous provenance (corpus A) and a
/// near-coplanar escalation-forcing corpus (corpus B). The interval filter →
/// astro-float judge must never claim a sign opposite to a GT-stable 512-bit truth
/// ([`truth::orient`]), and a constructed zero (corpus B's `ε = 0`, corpus C) must be judged
/// `Zero` — no oracle is asked there (see `truth`).
/// Asserts both paths are live: `escalated > 0` (filter defers) **and**
/// `filter_resolved > 0` (fast path resolves — a filter stuck at `None` would pass
/// wrong-sign 0 vacuously). `#[ignore]`: slow astro-float GT.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn indirect_orient3d_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let (mut zero_asserted, mut zero_signed) = (0usize, Vec::<Orient>::new());
    let mut check = |a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     q: &WitnessPoint,
                     r: &WitnessPoint,
                     s: &WitnessPoint,
                     zero: bool| {
        let judged = indirect_orient3d_judge(a, b, c, q, r, s, fixture()).orient();
        let planes = [
            plane_iv(a.0, a.1, a.2),
            plane_iv(b.0, b.1, b.2),
            plane_iv(c.0, c.1, c.2),
        ];
        if indirect_filter(planes, q.realized, r.realized, s.realized).is_none() {
            escalated += 1;
        } else {
            filter_resolved += 1;
        }
        if zero {
            zero_asserted += 1;
            if judged != Orient::Zero {
                zero_signed.push(judged);
            }
            return;
        }
        match truth::orient(a, b, c, q, r, s, GT) {
            None => skipped += 1,
            Some(truth) => {
                tested += 1;
                if judged == Orient::Zero {
                    declined += 1;
                } else if judged != truth {
                    wrong += 1;
                }
            }
        }
    };

    // Corpus A — heterogeneous provenance (each of the twelve points its own
    // rotation, origin or arbitrary pivot). Generic → mostly filter-resolved.
    let mut st = 0xC0C0_9999_ABAB_CDCDu64;
    for _ in 0..2000 {
        let p: Vec<WitnessPoint> = (0..12).map(|_| rand_point(&mut st)).collect();
        check(
            (&p[0], &p[1], &p[2]),
            (&p[3], &p[4], &p[5]),
            (&p[6], &p[7], &p[8]),
            &p[9],
            &p[10],
            &p[11],
            false,
        );
    }

    // Corpus B — one shared rotation about a shared pivot (affine → coplanarity
    // preserved), near-coplanar so the escalation's sign-resolution is exercised.
    let add3 = |a: [Rat; 3], b: [Rat; 3]| {
        [
            a[0].checked_add(b[0]).unwrap(),
            a[1].checked_add(b[1]).unwrap(),
            a[2].checked_add(b[2]).unwrap(),
        ]
    };
    let smul = |k: Rat, a: [Rat; 3]| {
        [
            a[0].checked_mul(k).unwrap(),
            a[1].checked_mul(k).unwrap(),
            a[2].checked_mul(k).unwrap(),
        ]
    };
    let cross = |a: [Rat; 3], b: [Rat; 3]| {
        [
            a[1].checked_mul(b[2])
                .unwrap()
                .checked_sub(a[2].checked_mul(b[1]).unwrap())
                .unwrap(),
            a[2].checked_mul(b[0])
                .unwrap()
                .checked_sub(a[0].checked_mul(b[2]).unwrap())
                .unwrap(),
            a[0].checked_mul(b[1])
                .unwrap()
                .checked_sub(a[1].checked_mul(b[0]).unwrap())
                .unwrap(),
        ]
    };
    for _ in 0..2000 {
        let off = |st: &mut u64| {
            [
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
            ]
        };
        let v = [
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
        ];
        let (u, w) = (off(&mut st), off(&mut st));
        let n = cross(u, w); // normal to the V-plane
        // ε off-plane push: 0 (exactly degenerate), or a wide window down to ~1e-15
        // where the interval filter is ambiguous but the true sign is definite.
        let eps = match rng(&mut st, 0, 2) {
            0 => ri(0, 1),
            1 => ri(1, rng(&mut st, 5_000, 200_000)),
            _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
        };
        let plane_pts = |st: &mut u64| [v, add3(v, off(st)), add3(v, off(st))];
        let a = plane_pts(&mut st);
        let b = plane_pts(&mut st);
        let c = plane_pts(&mut st);
        // Triangle q=v+u, r=v+w, s=v+u+w+ε·n (coplanar with V=v, s pushed ε off).
        let qb = add3(v, u);
        let rb = add3(v, w);
        let sb = add3(add3(add3(v, u), w), smul(eps, n));
        // One shared rotation (shared pivot) for all twelve points.
        let axis = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let rp = |p: [Rat; 3]| WitnessPoint::at(p).rotate_about(axis, ang, piv);
        let (a0, a1, a2) = (rp(a[0]), rp(a[1]), rp(a[2]));
        let (b0, b1, b2) = (rp(b[0]), rp(b[1]), rp(b[2]));
        let (c0, c1, c2) = (rp(c[0]), rp(c[1]), rp(c[2]));
        let (q, r, s) = (rp(qb), rp(rb), rp(sb));
        check(
            (&a0, &a1, &a2),
            (&b0, &b1, &b2),
            (&c0, &c1, &c2),
            &q,
            &r,
            &s,
            eps == ri(0, 1),
        );
    }

    // Corpus C — a **constructed** four-plane concurrency, because a random one never happens
    // and this is the configuration that matters: every plane is spanned by one shared point
    // `p` and two others, so `p` is both the meet of any three and a defining point of the
    // fourth. With `p` in the third slot, `indirect_hp`'s `row1 = Dvec − D·s` cancels to
    // nothing — the shape a tool edge lying in a target plane produces, and the one that let a
    // wrong sign through a suite of random corpora.
    for _ in 0..400 {
        let axis = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 1, 359_000), rng(&mut st, 1, 997)); // inexact
        let pivot = rand_base(&mut st);
        let turn = |b: [Rat; 3]| WitnessPoint::at(b).rotate_about(axis, ang, pivot);
        let p = turn(rand_base(&mut st));
        let sp: Vec<[WitnessPoint; 2]> = (0..4)
            .map(|_| [turn(rand_base(&mut st)), turn(rand_base(&mut st))])
            .collect();
        for qi in 0..4 {
            let t: Vec<usize> = (0..4).filter(|&i| i != qi).collect();
            check(
                (&p, &sp[t[0]][0], &sp[t[0]][1]),
                (&p, &sp[t[1]][0], &sp[t[1]][1]),
                (&p, &sp[t[2]][0], &sp[t[2]][1]),
                &sp[qi][0],
                &sp[qi][1],
                &p,
                true,
            );
        }
    }

    eprintln!(
        "[indirect_orient3d] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}; zero_asserted: {zero_asserted}"
    );
    assert!(
        zero_asserted > 0,
        "corpus must exercise the constructed-zero path"
    );
    assert!(
        zero_signed.is_empty(),
        "indirect_orient3d: a constructed zero was judged with a sign — {} of {zero_asserted}, first {:?}",
        zero_signed.len(),
        &zero_signed[..zero_signed.len().min(5)]
    );

    assert_eq!(
        wrong, 0,
        "indirect judge must never disagree with GT (soundness)"
    );
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path (else a stuck filter passes vacuously)"
    );
}

/// Fast port check (default suite): a handful of generic configs — the indirect
/// judge must match a moderate-precision truth (catches a transcription bug; full
/// statistical soundness is the `#[ignore]`d `plane_coefficient_tol_soundness` and
/// `indirect_orient3d_soundness`).
#[test]
fn indirect_port_check() {
    const GT: usize = 384;
    let mut st = 0x1DEA_2C11_9F00_5A5Au64;
    for _ in 0..16 {
        let p: Vec<WitnessPoint> = (0..12).map(|_| rand_point(&mut st)).collect();
        let a = (&p[0], &p[1], &p[2]);
        let b = (&p[3], &p[4], &p[5]);
        let c = (&p[6], &p[7], &p[8]);
        let judged = indirect_orient3d_judge(a, b, c, &p[9], &p[10], &p[11], fixture()).orient();
        if let Some(truth) = truth::orient(a, b, c, &p[9], &p[10], &p[11], GT) {
            assert!(
                judged == truth || judged == Orient::Zero,
                "indirect judge {judged:?} disagrees with truth {truth:?}"
            );
        }
    }
}

/// The three axis-perpendicular planes through integer point `p` (meet exactly at `p`,
/// tol 0) — an exact axis-aligned implicit point for the sanity oracle.
fn axis_planes(p: [i128; 3]) -> [[WitnessPoint; 3]; 3] {
    let pt = |x, y, z| WitnessPoint::at([ri(x, 1), ri(y, 1), ri(z, 1)]);
    let [x, y, z] = p;
    [
        [pt(x, y, z), pt(x, y + 1, z), pt(x, y, z + 1)], // ⊥ x
        [pt(x, y, z), pt(x + 1, y, z), pt(x, y, z + 1)], // ⊥ y
        [pt(x, y, z), pt(x + 1, y, z), pt(x, y + 1, z)], // ⊥ z
    ]
}

/// Sanity: two exact axis-aligned implicit points order by the compared axis; a swap
/// flips it; an equal coordinate is `Zero`.
#[test]
fn cmp_sanity_axis_aligned() {
    let a = axis_planes([1, 2, 3]);
    let b = axis_planes([1, 5, 3]);
    // y: 2 < 5 → a below b → Negative; swap → Positive.
    assert_eq!(
        indirect_cmp_coord_judge(tr(&a), tr(&b), 1, fixture()).orient(),
        Orient::Negative
    );
    assert_eq!(
        indirect_cmp_coord_judge(tr(&b), tr(&a), 1, fixture()).orient(),
        Orient::Positive
    );
    // x and z equal → Zero.
    assert_eq!(
        indirect_cmp_coord_judge(tr(&a), tr(&b), 0, fixture()).orient(),
        Orient::Zero
    );
    assert_eq!(
        indirect_cmp_coord_judge(tr(&a), tr(&b), 2, fixture()).orient(),
        Orient::Zero
    );
}

/// Fast port check: the cmp judge matches a moderate-precision truth on generic
/// rotated configs (catches a transcription bug; full soundness is the `#[ignore]`d
/// `indirect_cmp_coord_soundness`).
#[test]
fn cmp_port_check() {
    const GT: usize = 384;
    let mut st = 0xC301_7A5E_2266_9911u64;
    for _ in 0..16 {
        let a = rand_triple(&mut st);
        let b = rand_triple(&mut st);
        let axis = rng(&mut st, 0, 2) as usize;
        let judged = indirect_cmp_coord_judge(tr(&a), tr(&b), axis, fixture()).orient();
        if let Some(truth) = truth::cmp(tr(&a), tr(&b), axis, GT) {
            assert!(
                judged == truth || judged == Orient::Zero,
                "cmp judge {judged:?} disagrees with truth {truth:?}"
            );
        }
    }
}

fn rat_cross(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
    let m = |x: Rat, y: Rat| x.checked_mul(y).unwrap();
    let s = |x: Rat, y: Rat| x.checked_sub(y).unwrap();
    [
        s(m(a[1], b[2]), m(a[2], b[1])),
        s(m(a[2], b[0]), m(a[0], b[2])),
        s(m(a[0], b[1]), m(a[1], b[0])),
    ]
}

fn smul(k: Rat, a: [Rat; 3]) -> [Rat; 3] {
    [
        a[0].checked_mul(k).unwrap(),
        a[1].checked_mul(k).unwrap(),
        a[2].checked_mul(k).unwrap(),
    ]
}

/// Three base points defining a plane whose normal is parallel to `n` — `p0` and two
/// in-plane edges `n × e1`, `n × e2` (single crosses, small magnitude).
fn plane_norm(n: [Rat; 3], p0: [Rat; 3]) -> [[Rat; 3]; 3] {
    let e1 = [ri(1, 1), ri(2, 1), ri(3, 1)];
    let e2 = [ri(2, 1), ri(3, 1), ri(1, 1)];
    [p0, add3(p0, rat_cross(n, e1)), add3(p0, rat_cross(n, e2))]
}

fn rot_plane(b: [[Rat; 3]; 3], ax: Axis, ang: Angle, piv: [Rat; 3]) -> [WitnessPoint; 3] {
    b.map(|p| WitnessPoint::at(p).rotate_about(ax, ang, piv))
}

fn t3(p: &[WitnessPoint; 3]) -> (&WitnessPoint, &WitnessPoint, &WitnessPoint) {
    (&p[0], &p[1], &p[2])
}

/// The three points are far from collinear (`sin²` of the corner angle above a
/// threshold) — a well-formed plane. A degenerate plane (tiny normal) is not the
/// near-coplanar-*normals* regime under test, so the corpus skips it.
fn well_conditioned(p: &[WitnessPoint; 3]) -> bool {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let (e1, e2) = (
        sub(p[1].coord(), p[0].coord()),
        sub(p[2].coord(), p[0].coord()),
    );
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    dot(n, n) > 1e-6 * dot(e1, e1) * dot(e2, e2)
}

/// Sanity: the three axis-perpendicular planes have normals `+x, −y, +z`, so their
/// determinant is `−1`; a swap flips it, and three coplanar normals give `Zero`.
#[test]
fn dir_sign_judge_sanity() {
    let ap = axis_planes([0, 0, 0]);
    assert_eq!(
        dir_sign_judge(t3(&ap[0]), t3(&ap[1]), t3(&ap[2]), fixture()).orient(),
        Orient::Negative,
        "det[+x, -y, +z] = -1"
    );
    assert_eq!(
        dir_sign_judge(t3(&ap[0]), t3(&ap[2]), t3(&ap[1]), fixture()).orient(),
        Orient::Positive,
        "one swap flips the sign"
    );
    // Three normals in the plane z = 0 → coplanar → D = 0 → Zero.
    let mk = |n: [i128; 3]| {
        plane_norm([ri(n[0], 1), ri(n[1], 1), ri(n[2], 1)], [ri(0, 1); 3]).map(WitnessPoint::at)
    };
    let (a, b, c) = (mk([1, 0, 0]), mk([0, 1, 0]), mk([1, 1, 0]));
    assert_eq!(
        dir_sign_judge(t3(&a), t3(&b), t3(&c), fixture()).orient(),
        Orient::Zero,
        "coplanar normals → D = 0"
    );
}

/// `dir_sign` is a determinant of normals, invariant under a shared rotation
/// (`det(R·n) = det(R)·det(n) = det(n)`).
#[test]
fn dir_sign_rotation_invariant() {
    let mut st = 0x0D12_5157_ABCD_0007u64;
    for _ in 0..40 {
        let mk = |st: &mut u64| plane_norm(rand_base(st), rand_base(st));
        let (ba, bb, bc) = (mk(&mut st), mk(&mut st), mk(&mut st));
        let un = |b: [[Rat; 3]; 3]| b.map(WitnessPoint::at);
        let (ua, ub, uc) = (un(ba), un(bb), un(bc));
        let s = dir_sign_judge(t3(&ua), t3(&ub), t3(&uc), fixture()).orient();
        if s == Orient::Zero {
            continue;
        }
        let ax = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let (ra, rb, rc) = (
            rot_plane(ba, ax, ang, piv),
            rot_plane(bb, ax, ang, piv),
            rot_plane(bc, ax, ang, piv),
        );
        assert_eq!(
            dir_sign_judge(t3(&ra), t3(&rb), t3(&rc), fixture()).orient(),
            s,
            "rotation-invariant"
        );
    }
}

/// `dir_sign` soundness over a **near-coplanar-normals** corpus (which the indirect orient3d and
/// cmp_coord soundness tests do not stress: they force `M ≈ 0`, not `D ≈ 0`). Three plane normals `n0, n1,
/// n2 = α·n0 + β·n1 + ε·(n0×n1)` (ε tiny → `D ≈ ε` → escalation), shared rotation. The
/// judge must never disagree with a GT-stable 512-bit `D` sign ([`truth::dir_sign`]), and must
/// judge the constructed zero `ε = 0` `Zero`; both paths exercised.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn dir_sign_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let (mut zero_asserted, mut zero_signed) = (0usize, Vec::<Orient>::new());
    let mut st = 0x0D18_5160_C0C0_2323u64;
    let rand_n = |st: &mut u64| {
        [
            ri(rng(st, -20, 20), 1),
            ri(rng(st, -20, 20), 1),
            ri(rng(st, -20, 20), 1),
        ]
    };
    for _ in 0..2000 {
        let n0 = rand_n(&mut st);
        let n1 = rand_n(&mut st);
        let (alpha, beta) = (ri(rng(&mut st, -5, 5), 1), ri(rng(&mut st, -5, 5), 1));
        let eps = match rng(&mut st, 0, 2) {
            0 => ri(0, 1),
            1 => ri(1, rng(&mut st, 5_000, 200_000)),
            _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
        };
        // n2 = α·n0 + β·n1 + ε·(n0×n1) — near-coplanar with n0, n1.
        let n2 = add3(
            add3(smul(alpha, n0), smul(beta, n1)),
            smul(eps, rat_cross(n0, n1)),
        );
        let ax = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let mk = |n: [Rat; 3], st: &mut u64| rot_plane(plane_norm(n, rand_base(st)), ax, ang, piv);
        let pa = mk(n0, &mut st);
        let pb = mk(n1, &mut st);
        let pc = mk(n2, &mut st);
        if !(well_conditioned(&pa) && well_conditioned(&pb) && well_conditioned(&pc)) {
            continue; // a degenerate plane is not the regime under test
        }
        let judged = dir_sign_judge(t3(&pa), t3(&pb), t3(&pc), fixture()).orient();
        let (d, _) = cramer_iv([
            plane_iv(&pa[0], &pa[1], &pa[2]),
            plane_iv(&pb[0], &pb[1], &pb[2]),
            plane_iv(&pc[0], &pc[1], &pc[2]),
        ]);
        if d.sign().is_none() {
            escalated += 1;
        } else {
            filter_resolved += 1;
        }
        let zero = eps == ri(0, 1);
        if zero {
            zero_asserted += 1;
            if judged != Orient::Zero {
                zero_signed.push(judged);
            }
            continue;
        }
        match truth::dir_sign(t3(&pa), t3(&pb), t3(&pc), GT) {
            None => skipped += 1,
            Some(truth) => {
                tested += 1;
                if judged == Orient::Zero {
                    declined += 1;
                } else if judged != truth {
                    wrong += 1;
                }
            }
        }
    }
    eprintln!(
        "[dir_sign] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}; zero_asserted: {zero_asserted}"
    );
    assert!(
        zero_asserted > 0,
        "corpus must exercise the constructed-zero path"
    );
    assert!(
        zero_signed.is_empty(),
        "dir_sign: a constructed zero was judged with a sign — {} of {zero_asserted}, first {:?}",
        zero_signed.len(),
        &zero_signed[..zero_signed.len().min(5)]
    );

    assert_eq!(wrong, 0, "dir_sign must never disagree with GT (soundness)");
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path"
    );
}

/// Sanity: `det[d, x−base, y−base] = d·((x−base)×(y−base))`. With `base=0, x=e_x, y=e_y`
/// the edge normal is `+e_z`, so `d=+e_z → Positive`, `−e_z → Negative`, in-plane →
/// `Zero`; swapping the edge pair flips the sign.
#[test]
fn dir_orient3d_judge_sanity() {
    let o = WitnessPoint::at([ri(0, 1); 3]);
    let x = WitnessPoint::at([ri(1, 1), ri(0, 1), ri(0, 1)]);
    let y = WitnessPoint::at([ri(0, 1), ri(1, 1), ri(0, 1)]);
    let e = |a: i128, b: i128, c: i128| [ri(a, 1), ri(b, 1), ri(c, 1)];
    assert_eq!(
        dir_orient3d_judge(e(0, 0, 1), &o, &x, &y, FIXTURE_PREC),
        Orient::Positive
    );
    assert_eq!(
        dir_orient3d_judge(e(0, 0, -1), &o, &x, &y, FIXTURE_PREC),
        Orient::Negative
    );
    assert_eq!(
        dir_orient3d_judge(e(1, 0, 0), &o, &x, &y, FIXTURE_PREC),
        Orient::Zero,
        "d in the edge plane → det 0"
    );
    assert_eq!(
        dir_orient3d_judge(e(0, 0, 1), &o, &y, &x, FIXTURE_PREC),
        Orient::Negative,
        "swapping x,y flips the sign"
    );
}

/// `orient3d_ray(base, dir, x, y)` is exactly `orient3d(base, base+dir, x, y)`: on
/// unrotated points (where `base+dir` *is* an exact `WitnessPoint`) it must equal `orient3d_judge`
/// with the ideal point materialized — the argument-for-argument reduction the ops
/// ray-triangle caller relies on.
#[test]
fn orient3d_ray_matches_materialized() {
    let mut st = 0x0D1B_7A44_0F0F_9001u64;
    let at = |b: [Rat; 3]| WitnessPoint::at(b);
    for _ in 0..200 {
        let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
        let dir = [
            ri(rng(&mut st, -9, 9), 1),
            ri(rng(&mut st, -9, 9), 1),
            ri(rng(&mut st, -9, 9), 1),
        ];
        let q = [
            bb[0].checked_add(dir[0]).unwrap(),
            bb[1].checked_add(dir[1]).unwrap(),
            bb[2].checked_add(dir[2]).unwrap(),
        ];
        let (base, x, y) = (at(bb), at(xb), at(yb));
        assert_eq!(
            orient3d_ray(&base, dir, &x, &y, FIXTURE_PREC),
            orient3d_judge(&base, &at(q), &x, &y, fixture()).orient(),
            "orient3d_ray == orient3d(base, base+dir, x, y)"
        );
    }
}

/// Soundness — `dir_orient3d` over a rotated corpus, oracle = a GT-stable 512-bit
/// determinant (NOT rotation-invariance: `d` is fixed while the points rotate, so the
/// sign is not preserved). Two regimes: random directions, and **near-grazing** (`d`
/// almost in the edge plane, `det ≈ 0`) to force escalation. The judge must never claim
/// a sign opposite to the GT; both the fast filter and astro-float paths are exercised.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn dir_orient3d_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let mut st = 0x0D1A_5170_BEEF_4242u64;
    for _ in 0..1500 {
        let (bb, xb, yb) = (rand_base(&mut st), rand_base(&mut st), rand_base(&mut st));
        let grazing = rng(&mut st, 0, 1) == 1;
        let ax = axis_of(rng(&mut st, 0, 2));
        let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let rp = |b: [Rat; 3]| WitnessPoint::at(b).rotate_about(ax, ang, piv);
        let (base, x, y) = (rp(bb), rp(xb), rp(yb));
        let tri = [base.clone(), x.clone(), y.clone()];
        if !well_conditioned(&tri) {
            continue; // a degenerate edge fan is not the regime under test
        }
        // Two direction regimes. Random: any integer `d`. near-grazing: a large integer
        // multiple of the **rotated** in-plane edge `x−base` (built from the f64 coords),
        // so `d` is nearly ⊥ the rotated normal → `det ≈ 0` → the fast filter fails and
        // the astro-float path decides. Built after rotation, since the rotated normal is
        // what `d` must graze.
        let d = if !grazing {
            [
                ri(rng(&mut st, -19, 19), 1),
                ri(rng(&mut st, -19, 19), 1),
                ri(rng(&mut st, -19, 19), 1),
            ]
        } else {
            let big = rng(&mut st, 100_000_000, 9_000_000_000) as f64;
            let comp = |k: usize| {
                ri(
                    (big * (x.realized[k].value - base.realized[k].value)).round() as i128,
                    1,
                )
            };
            [comp(0), comp(1), comp(2)]
        };
        let judged = dir_orient3d_judge(d, &base, &x, &y, FIXTURE_PREC);
        // Recompute the Bounded filter to tally which path resolved (mirrors the judge).
        let dp = WitnessPoint::at(d);
        let (bi, xi, yi, di) = (base.realized, x.realized, y.realized, dp.realized);
        let subi =
            |u: [Bounded; 3], v: [Bounded; 3]| [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
        if det3_iv([di, subi(xi, bi), subi(yi, bi)]).sign().is_none() {
            escalated += 1;
        } else {
            filter_resolved += 1;
        }
        match truth::dir_orient(d, &base, &x, &y, GT) {
            None => skipped += 1,
            Some(truth) => {
                tested += 1;
                if judged == Orient::Zero {
                    declined += 1;
                } else if judged != truth {
                    wrong += 1;
                }
            }
        }
    }
    eprintln!(
        "[dir_orient3d] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
    );
    assert_eq!(
        wrong, 0,
        "dir_orient3d must never disagree with GT (soundness)"
    );
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path"
    );
}
