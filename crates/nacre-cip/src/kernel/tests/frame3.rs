use super::*;

/// A witness at an `f64`-representable point, stated as the rational it is — the fixture
/// spelling of what `WitnessPoint::exact` used to be. That door is retired: no production
/// caller, and a precondition ("exactly representable") no caller could check — the one
/// production site that handed it a rounded cache named a different point (nacre-ops
/// reuse). `at_nearest` states the same tol `0` here and stays honest elsewhere.
fn exact(c: [f64; 3]) -> Option<WitnessPoint> {
    let b = |x: f64| Rat::try_from_f64(x);
    Some(WitnessPoint::at_nearest([b(c[0])?, b(c[1])?, b(c[2])?]))
}

/// ★★★★ **One pass equals node by node — coordinate, tol and realization, bit for bit.**
///
/// [`WitnessPoint::apply_chain`] exists only to stop rebuilding the remembered chain per node
/// (`remember`'s doc has the measurement). It is allowed to be faster; it is not allowed to be
/// different, and "different" here means the last bit of a coordinate the cache will store.
#[test]
fn one_pass_equals_node_by_node() {
    let r = Rat::from_int;
    let deg = |d: i128| Angle::from_deg(r(d)).expect("a whole-degree angle");
    let chains: Vec<Vec<MoveNode>> = vec![
        // An irrational turn about an offset pivot, so the pivot arithmetic is charged too.
        (0..40)
            .map(|_| MoveNode::Rotate {
                axis: Axis::Z,
                angle: deg(37),
                pivot: [r(1), r(-2), r(0)],
            })
            .collect(),
        // Exact inputs: translations and a mirror, whose radii must stay where they were.
        (0..40)
            .map(|i| match i % 2 {
                0 => MoveNode::Translate {
                    offset: [Rat::new(1, 7).expect("1/7"), r(3), r(0)],
                },
                _ => MoveNode::Mirror {
                    axis: Axis::X,
                    offset: Rat::new(1, 2).expect("1/2"),
                },
            })
            .collect(),
        // Mixed, and quadrantal turns among them (exact cos/sin — the tol-0 arm).
        (0..40)
            .map(|i| match i % 3 {
                0 => MoveNode::Rotate {
                    axis: Axis::Y,
                    angle: deg(90),
                    pivot: [r(0); 3],
                },
                1 => MoveNode::Translate {
                    offset: [r(0), r(0), Rat::new(2, 5).expect("2/5")],
                },
                _ => MoveNode::Rotate {
                    axis: Axis::X,
                    angle: deg(11),
                    pivot: [r(0); 3],
                },
            })
            .collect(),
    ];
    for (i, chain) in chains.iter().enumerate() {
        let base = [r(2), r(3), r(4)];
        let one_pass = WitnessPoint::at(base).apply_chain(chain).expect("one pass");
        let node_by_node = chain.iter().fold(WitnessPoint::at(base), |q, n| match n {
            MoveNode::Rotate { axis, angle, pivot } => q.rotate_about(*axis, *angle, *pivot),
            MoveNode::Translate { offset } => q.translate(*offset),
            MoveNode::Mirror { axis, offset } => q.mirror(*axis, *offset),
            _ => unreachable!("the fixtures above hold no frame node"),
        });
        assert_eq!(
            one_pass.coord(),
            node_by_node.coord(),
            "chain {i}: coordinate"
        );
        assert_eq!(one_pass.tol(), node_by_node.tol(), "chain {i}: tol");
        assert_eq!(
            one_pass.chain.len(),
            node_by_node.chain.len(),
            "chain {i}: the definition is remembered whole"
        );
        for (a, b) in one_pass
            .realize(256)
            .iter()
            .zip(node_by_node.realize(256).iter())
        {
            assert_eq!(a.error, b.error, "chain {i}: realization radius");
            assert_eq!(
                nacre_scalar::round_to_f64(&a.value, a.error, 256),
                nacre_scalar::round_to_f64(&b.value, b.error, 256),
                "chain {i}: realized coordinate"
            );
        }
    }
}

/// The precision these fixtures judge at. Production chooses it per model
/// ([`judge_precision`]); a fixture pins one so its expectations stay fixed.
const FIXTURE_PREC: usize = 256;

/// The judgement context those fixtures use: that precision, plus the coincidence limit the
/// default derivation gives a unit-scale model (output resolution `2⁻⁵²`, two words further
/// down) and the production cap.
fn fixture() -> Standard {
    Standard {
        prec: FIXTURE_PREC,
        coincidence: Mag::pow2(-180),
        scale: Mag::of(1.0),
        cap: 4096,
    }
}

fn ri(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}

fn deg(n: i128, d: i128) -> Angle {
    Angle::from_deg(ri(n, d)).unwrap()
}

/// **The climb lands where the shortfall says, in one jump.**
///
/// The gap is `C · 2⁻ᵖʳᵉᶜ` over a cofactor and `C` does not move with the precision, so
/// `log₂(gap / limit)` is exactly the number of missing bits — the same reading
/// [`judge_precision`] takes to size the model. This pins that the loop uses it: the fixture's
/// gap shrinks bit-for-bit with the precision, it starts 100 bits short, and the next attempt
/// must arrive at 384 — `256 + 100` rounded up to a word — not at 512 (doubling), and not at
/// 320 (a fixed word step, which would need two more rounds).
#[test]
fn the_climb_lands_where_the_shortfall_names_in_one_jump() {
    let mut asked = Vec::new();
    let j = Standard {
        prec: 256,
        coincidence: Mag::pow2(-300),
        scale: Mag::of(1.0),
        cap: 4096,
    };
    let out = escalate(j, j.coincidence, |prec| {
        asked.push(prec);
        // `C = 2⁵⁶`: at 256 bits the gap is `2⁻²⁰⁰`, a hundred bits above the limit.
        Err(Gap::Of(Mag::pow2(56 - prec as i64)))
    });
    assert_eq!(asked, vec![256, 384], "the climb overshot or crept");
    assert!(
        matches!(out, Decision::Coincident { within } if within.lt(Mag::pow2(-299))),
        "{out:?} — the second attempt was inside the limit and had to be reported as proof"
    );
}

/// **At the cap the judgement says so; it does not quietly become a zero.**
///
/// This is the failure the whole cell exists to remove — an undecided determinant that reads
/// as "these coincide" and merges two things that are apart. A gap that never reaches the
/// limit must come back [`Decision::Exhausted`], and `orient()` may collapse it to `Zero` for
/// the geometry only because the outcome itself is still there to be reported.
#[test]
fn a_gap_that_never_reaches_the_limit_is_exhausted_not_coincident() {
    let j = Standard {
        prec: 256,
        coincidence: Mag::pow2(-300),
        scale: Mag::of(1.0),
        cap: 512,
    };
    let mut rounds = 0;
    // A gap that ignores the precision entirely — a cofactor collapsing as fast as the
    // radius shrinks. No depth settles it, which is what the cap is for.
    let out = escalate(j, j.coincidence, |_| {
        rounds += 1;
        Err(Gap::Of(Mag::pow2(-10)))
    });
    // **The bound it did establish rides out with it.** Without it the caller cannot tell a
    // judgement that stopped `2⁻¹⁰` short from one that stopped `2⁻³⁰⁰` short, and those are
    // not the same news about the model.
    assert_eq!(
        out,
        Decision::Exhausted {
            at: 512,
            within: Some(Mag::pow2(-10)),
        }
    );
    assert_eq!(out.orient(), Orient::Zero, "the geometry still gets a sign");
    assert!(rounds > 1 && rounds < 10, "climbed {rounds} times");

    // An unresolved cofactor that never resolves has no distance to quote at all — reporting
    // one would be inventing it.
    let out = escalate(j, j.coincidence, |_| Err(Gap::Unresolved { short: 64 }));
    assert_eq!(
        out,
        Decision::Exhausted {
            at: 512,
            within: None,
        }
    );
}

/// **A degenerate witness is a different answer from an exhausted one, and must not climb.**
///
/// No gap at all means the quantity that turns the determinant into a distance could not be
/// bounded away from zero — a collapsed witness triangle, or three planes with no meeting
/// point. More bits sharpen a distance; they do not conjure one. Telling the two apart is the
/// whole reason the outcome is not a single "could not decide": one cause is a budget, the
/// other is the geometry, and a caller that cannot tell them apart cannot act on either.
#[test]
fn a_degenerate_witness_is_told_apart_from_an_exhausted_one() {
    let j = Standard {
        prec: 256,
        coincidence: Mag::pow2(-300),
        scale: Mag::of(1.0),
        cap: 4096,
    };
    let mut rounds = 0;
    let out = escalate(j, j.coincidence, |_| {
        rounds += 1;
        Err(Gap::Vanished)
    });
    assert_eq!(out, Decision::Degenerate);
    assert_eq!(rounds, 1, "a vanished cofactor must not be re-realized");

    // …but a cofactor that is merely *unresolved* is the other cause, and it must climb —
    // reporting it as degenerate is how a judgement nobody looked at closely enough turns
    // into a silent merge. It names its own shortfall, so the climb is one jump here too.
    let mut asked = Vec::new();
    let out = escalate(j, j.coincidence, |prec| {
        asked.push(prec);
        if prec < 512 {
            Err(Gap::Unresolved { short: 200 })
        } else {
            Ok(Orient::Positive)
        }
    });
    assert_eq!(out, Decision::Sign(Orient::Positive));
    assert_eq!(
        asked,
        vec![256, 512],
        "an unresolved cofactor did not climb"
    );
}

/// **A structural zero comes back proved, not assumed.**
///
/// Four points that are coplanar *before* a rotation stay coplanar after it, but the judge
/// cannot see that: one base is `1/3`, so the pre-rotation shortcut cannot hand the question
/// to the exact predicate and the determinant is a transcendental zero no finite precision
/// separates. Before the gap was wired in, that came back `Orient::Zero` with nothing behind
/// it. Now it comes back with the distance it established — and that distance has to be
/// **below the coincidence limit**, which is what makes "these are the same plane" a proof
/// rather than a shrug.
///
/// The plane is **tilted** (`z = x`), which the first version of this fixture was not: a plane
/// perpendicular to the rotation axis keeps its `z` coordinate exactly, so every error term
/// meets an exactly-zero factor and the determinant comes back a *proved* zero instead — a
/// real outcome (and one this loop reports), but not the one under test here.
#[test]
fn a_structural_zero_comes_back_proved_not_assumed() {
    // The plane `z = x` through four points, one of them irrational, all turned together.
    let pt = |x: Rat, y: Rat| {
        WitnessPoint::at([x, y, x]).rotate_about(
            Axis::Z,
            deg(30, 1),
            [ri(1, 3), ri(1, 7), ri(0, 1)],
        )
    };
    let (a, b, c, d) = (
        pt(ri(1, 3), ri(1, 1)),
        pt(ri(0, 1), ri(0, 1)),
        pt(ri(1, 1), ri(0, 1)),
        pt(ri(0, 1), ri(1, 1)),
    );
    let j = fixture();
    let out = orient3d_judge(&a, &b, &c, &d, j);
    let Decision::Coincident { within } = out else {
        panic!("{out:?} — a structural zero must be proved, not assumed");
    };
    assert!(
        !j.coincidence.lt(within),
        "reported {:?} against a limit of {:?}",
        within.exp2(),
        j.coincidence.exp2()
    );
    assert_eq!(out.orient(), Orient::Zero);

    // …and the limit is genuinely consulted: ask for a coincidence a thousand bits finer than
    // the model can carry and the same judgement must refuse to call it one.
    let strict = Standard {
        coincidence: Mag::pow2(-2000),
        cap: 512,
        ..j
    };
    let strict_out = orient3d_judge(&a, &b, &c, &d, strict);
    let Decision::Exhausted {
        at: 512,
        within: Some(w),
    } = strict_out
    else {
        panic!("{strict_out:?} — the coincidence limit was ignored, any zero would pass");
    };
    // What it *could* show is still the honest number, and it is nowhere near the limit asked
    // for: the report says "within this much", not "these are the same".
    assert!(strict.coincidence.lt(w), "reported {:?}", w.exp2());
}

/// **The design principle, as a test: raising the precision must actually narrow the bound.**
///
/// The old declare-0 floor was a constant, so a deeper escalation bought a smaller threshold
/// only because `2⁻ᵖʳᵉᶜ` appeared in it — the *estimate* it multiplied never improved. A
/// computed radius has to do better than that, and if it did not, every "a deeper rung settles
/// nothing" measurement would be reporting a broken radius rather than a fact about the
/// geometry. So this pins the response directly: doubling the precision must take roughly a
/// factor of `2⁻ᵖʳᵉᶜ` off the radius, both for a realized coordinate and for a determinant
/// built out of several of them.
#[test]
fn a_deeper_precision_actually_narrows_the_radius() {
    // A rotation with an irrational cos/sin, about a non-origin pivot, so the radius has
    // every term in it: the trig bound, the pivot arithmetic, and the per-operation rounding.
    let pt = |x: i128, y: i128, z: i128| {
        WitnessPoint::at([ri(x, 10), ri(y, 10), ri(z, 10)]).rotate_about(
            Axis::Z,
            deg(37, 1),
            [ri(1, 3), ri(1, 7), ri(0, 1)],
        )
    };
    let (a, b, c, d) = (pt(3, 5, 0), pt(11, 2, 4), pt(-7, 9, 13), pt(1, -6, 2));

    let mut prev_coord: Option<i64> = None;
    let mut prev_det: Option<i64> = None;
    for prec in [256usize, 512, 1024, 2048] {
        let coord = a.hp_coord(prec)[0]
            .error
            .exp2()
            .expect("a rotated coordinate is not exact, so its radius is not zero");
        let det = det3_hp(&a, &b, &c, &d, prec)
            .error
            .exp2()
            .expect("a determinant over rotated points carries a radius");
        if let (Some(pc), Some(pd)) = (prev_coord, prev_det) {
            // Each step doubles `prec`, so the radius should drop by about that many binary
            // orders. Demand most of it — the seed terms and the operation count add a
            // constant offset that does not shrink, but it must not dominate.
            let want = (prec / 2) as i64 - 32;
            assert!(
                pc - coord >= want,
                "coordinate radius went 2^{pc} → 2^{coord} at {prec} bits: only {} orders, \
                     wanted {want}. A radius that ignores precision makes the ladder meaningless.",
                pc - coord
            );
            assert!(
                pd - det >= want,
                "determinant radius went 2^{pd} → 2^{det} at {prec} bits: only {} orders, \
                     wanted {want}",
                pd - det
            );
        }
        prev_coord = Some(coord);
        prev_det = Some(det);
    }
}

/// **Is the normalized determinant actually a distance?**
///
/// The coincidence limit is a length, so the quantity compared against it has to be one too.
/// `orient3d` is a signed volume; dividing by the area term is supposed to leave a height.
/// Nothing else in the suite would notice if that division were wrong by a factor of the
/// triangle's size — the sign would still be right, and only the *threshold* would silently
/// become a number that scales with the model.
///
/// So this builds a point a **known** height above a plane and checks the normalized value is
/// that height: across heights spanning 12 orders of magnitude, model scales spanning 9, and
/// a triangle deliberately made 1000× larger, which is exactly what an unnormalized
/// determinant would be fooled by.
#[test]
fn the_normalized_determinant_is_a_height_in_model_units() {
    let prec = 256;
    for scale in [1i128, 1_000, 1_000_000_000] {
        for (hn, hd) in [(1i128, 1i128), (1, 1_000), (1, 1_000_000_000_000)] {
            for tri_span in [1i128, 1_000] {
                let p = |x: i128, y: i128, z: (i128, i128)| {
                    WitnessPoint::at([ri(x * scale, 1), ri(y * scale, 1), ri(z.0, z.1)])
                };
                // Plane z = 0 through three points, and the query a height h above it.
                let (b, c, d) = (
                    p(tri_span, 0, (0, 1)),
                    p(0, tri_span, (0, 1)),
                    p(0, 0, (0, 1)),
                );
                let a = p(0, 0, (hn, hd));
                let det = det3_hp(&a, &b, &c, &d, prec);
                let sub = |x: &HpBounded, y: &HpBounded| x.sub(y, prec);
                let (bh, ch, dh) = (b.hp_coord(prec), c.hp_coord(prec), d.hp_coord(prec));
                let (u, v) = (
                    [
                        sub(&bh[0], &dh[0]),
                        sub(&bh[1], &dh[1]),
                        sub(&bh[2], &dh[2]),
                    ],
                    [
                        sub(&ch[0], &dh[0]),
                        sub(&ch[1], &dh[1]),
                        sub(&ch[2], &dh[2]),
                    ],
                );
                let cross = [
                    u[1].mul(&v[2], prec).sub(&u[2].mul(&v[1], prec), prec),
                    u[2].mul(&v[0], prec).sub(&u[0].mul(&v[2], prec), prec),
                    u[0].mul(&v[1], prec).sub(&u[1].mul(&v[0], prec), prec),
                ];
                // Squared, so no square root is needed: `det² / |cross|²` must be `h²`.
                let sq = |x: &BigFloat| x.mul(x, prec, HP_RM);
                let norm2 = cross.iter().fold(BigFloat::from_f64(0.0, prec), |acc, k| {
                    acc.add(&sq(&k.value), prec, HP_RM)
                });
                let got = sq(&det.value).div(&norm2, prec, HP_RM);
                let h =
                    BigFloat::from_i128(hn, prec).div(&BigFloat::from_i128(hd, prec), prec, HP_RM);
                let err = rel_err(&got, &sq(&h), prec);
                assert!(
                    err < 1e-6,
                    "scale {scale}, h {hn}/{hd}, triangle span {tri_span}: the normalized \
                         value is not the height (relative error {err:e})"
                );
            }
        }
    }
}

/// **Is `cmp_coord`'s normalized value actually a coordinate difference?**
///
/// The same trap as the orient3d height: `M` alone carries both Cramer denominators, so a
/// threshold applied to it moves when the planes are scaled — and the sign, which is all the
/// suite checks, does not move at all. Two implicit points a **known** distance apart, with
/// the plane coefficients deliberately scaled by 1000 (which multiplies `M` by 10⁹ and must
/// leave the gap untouched).
#[test]
fn the_normalized_cmp_is_a_coordinate_difference() {
    let prec = 256;
    for (gn, gd) in [(1i128, 1i128), (7, 100), (1, 1_000_000_000)] {
        for k in [1i128, 1_000] {
            // x = 0, y = 0, z = 0  and  x = 0, y = 0, z = g: two points on the z axis, `g`
            // apart. Coefficients scaled by `k` — the same planes, differently written.
            let pl = |c: [i128; 4]| {
                [
                    HpBounded::exact(BigFloat::from_i128(c[0] * k, prec)),
                    HpBounded::exact(BigFloat::from_i128(c[1] * k, prec)),
                    HpBounded::exact(BigFloat::from_i128(c[2] * k, prec)),
                    // `d` carries the offset, which must not be scaled away: `g = gn/gd`.
                    HpBounded::new(
                        BigFloat::from_i128(c[3] * k, prec).div(
                            &BigFloat::from_i128(gd, prec),
                            prec,
                            HP_RM,
                        ),
                        Mag::ZERO,
                    ),
                ]
            };
            let x0 = pl([1, 0, 0, 0]);
            let y0 = pl([0, 1, 0, 0]);
            let z0 = pl([0, 0, 1, 0]);
            let zg = pl([0, 0, 1, -gn]); // z = gn/gd
            let (da, dva) = cramer_hp(&[x0.clone(), y0.clone(), z0], prec);
            let (db, dvb) = cramer_hp(&[x0, y0, zg], prec);
            let m = dva[2].mul(&db, prec).sub(&dvb[2].mul(&da, prec), prec);
            // The midpoint version of `coord_gap`: |M| / |D_a·D_b| must be the separation.
            let got = m
                .value
                .div(&da.value.mul(&db.value, prec, HP_RM), prec, HP_RM)
                .abs();
            let want =
                BigFloat::from_i128(gn, prec).div(&BigFloat::from_i128(gd, prec), prec, HP_RM);
            let err = rel_err(&got, &want, prec);
            assert!(
                err < 1e-6,
                "gap {gn}/{gd}, coefficients scaled by {k}: the normalized value is not the \
                     separation (relative error {err:e})"
            );
        }
    }
}

/// **Is the indirect orient3d's normalized value a point-to-plane distance?**
///
/// This one carries *two* denominators — `|D|` because the implicit point is `Dvec/D`, and
/// `|cross|` because the dot product carries the triangle's area — so there are two separate
/// ways for it to stop being a length. The fixture scales each independently: the plane
/// coefficients by `k` (which moves `D`) and the query triangle by `t` (which moves `cross`),
/// while the true distance stays put.
#[test]
fn the_normalized_indirect_orient3d_is_a_distance() {
    let prec = 256;
    let big = |v: i128| HpBounded::exact(BigFloat::from_i128(v, prec));
    for (hn, hd) in [(1i128, 1i128), (3, 100), (1, 1_000_000)] {
        for k in [1i128, 1_000] {
            for t in [1i128, 1_000] {
                // V = ∩(x=0, y=0, z=h) sits `h` above the plane z = 0 through the triangle
                // (0,0,0), (t,0,0), (0,t,0).
                let pl = |c: [i128; 3], d: (i128, i128)| {
                    [
                        big(c[0] * k),
                        big(c[1] * k),
                        big(c[2] * k),
                        HpBounded::new(
                            BigFloat::from_i128(d.0 * k, prec).div(
                                &BigFloat::from_i128(d.1, prec),
                                prec,
                                HP_RM,
                            ),
                            Mag::ZERO,
                        ),
                    ]
                };
                let planes = [
                    pl([1, 0, 0], (0, 1)),
                    pl([0, 1, 0], (0, 1)),
                    pl([0, 0, 1], (-hn, hd)),
                ];
                let pt = |x: i128, y: i128| [big(x), big(y), big(0)];
                let (d, m, _) = indirect_hp(planes, pt(t, 0), pt(0, t), pt(0, 0), prec);
                // |M| / (|D| · |cross|); `cross` here is (t,0,0)×(0,t,0) = (0,0,t²).
                // `cross` here is (t,0,0)×(0,t,0) = (0,0,t²), so `|cross| = t²`.
                let norm = BigFloat::from_i128(t * t, prec);
                let got = m
                    .value
                    .div(&d.value.mul(&norm, prec, HP_RM), prec, HP_RM)
                    .abs();
                let want =
                    BigFloat::from_i128(hn, prec).div(&BigFloat::from_i128(hd, prec), prec, HP_RM);
                let err = rel_err(&got, &want, prec);
                assert!(
                    err < 1e-6,
                    "h {hn}/{hd}, coefficients x{k}, triangle x{t}: the normalized value is \
                         not the distance (relative error {err:e})"
                );
            }
        }
    }
}

/// **Is `dir_sign`'s normalized value an angle?**
///
/// This is the one judge whose answer is not a length: it asks whether three plane normals
/// are coplanar, and its determinant carries each normal's magnitude — which is arbitrary,
/// since a plane's coefficients scale freely. Normals `(1,0,0)`, `(0,1,0)`, `(0,1,ε)` miss
/// coplanarity by `ε/√(1+ε²)` in sine; scaling any coefficient set must leave that alone,
/// and without the normalization it does not.
#[test]
fn the_normalized_dir_sign_is_an_angle() {
    let prec = 256;
    let big = |v: i128| HpBounded::exact(BigFloat::from_i128(v, prec));
    for (en, ed) in [(1i128, 1i128), (1, 1_000), (1, 1_000_000_000)] {
        for k in [1i128, 1_000] {
            let eps = HpBounded::new(
                BigFloat::from_i128(en * k, prec).div(&BigFloat::from_i128(ed, prec), prec, HP_RM),
                Mag::ZERO,
            );
            let planes = [
                [big(1), big(0), big(0), big(0)],
                [big(0), big(1), big(0), big(0)],
                [big(0), big(k), eps, big(0)],
            ];
            let (d, _) = cramer_hp(&planes, prec);
            // |D| / (|n_a||n_b||n_c|), with each |n| taken as its largest component.
            // Squared again, so the three norms need no square root: `D² / Π|n|²` is `sin²`.
            let sq = |x: &BigFloat| x.mul(x, prec, HP_RM);
            let norm2 = planes.iter().fold(BigFloat::from_f64(1.0, prec), |acc, p| {
                let n2 = (0..3).fold(BigFloat::from_f64(0.0, prec), |a, i| {
                    a.add(&sq(&p[i].value), prec, HP_RM)
                });
                acc.mul(&n2, prec, HP_RM)
            });
            let got = sq(&d.value).div(&norm2, prec, HP_RM);
            let sine =
                BigFloat::from_i128(en, prec).div(&BigFloat::from_i128(ed, prec), prec, HP_RM);
            // The exact sine of the miss is `ε/√(1+ε²)` — `ε` only for small `ε`, and the
            // fixture spans up to `ε = 1` where the two differ by √2. Squared: `ε²/(1+ε²)`.
            let s2 = sq(&sine);
            let want2 = s2.div(
                &s2.add(&BigFloat::from_f64(1.0, prec), prec, HP_RM),
                prec,
                HP_RM,
            );
            let err = rel_err(&got, &want2, prec);
            assert!(
                err < 1e-6,
                "sine {en}/{ed}, coefficients x{k}: the normalized value is not the angle \
                     (relative error {err:e})"
            );
        }
    }
}

/// Relative error of `got` against `want`, as an `f64` magnitude — computed **in astro-float**
/// so the comparison is not limited by `bf_mag`'s power-of-two rounding. Used by the
/// normalization tests, where a factor-of-two slop would hide a real units error.
fn rel_err(got: &BigFloat, want: &BigFloat, prec: usize) -> f64 {
    let d = got.sub(want, prec, HP_RM);
    bf_mag(&d.abs()) / bf_mag(&want.abs())
}

/// Deterministic PRNG (splitmix64) for reproducible stress corpora.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}
fn rng(state: &mut u64, lo: i128, hi: i128) -> i128 {
    lo + (u128::from(splitmix64(state)) % (hi - lo + 1) as u128) as i128
}
fn axis_of(k: i128) -> Axis {
    match k.rem_euclid(3) {
        0 => Axis::X,
        1 => Axis::Y,
        _ => Axis::Z,
    }
}
fn rand_base(st: &mut u64) -> [Rat; 3] {
    [
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
        ri(rng(st, -100_000, 100_000), rng(st, 1, 100)),
    ]
}
/// The f64 value's distance from the high-precision realization's **midpoint** — what the
/// f64 tol has to bound. (The realization's own radius is a separate, far smaller quantity;
/// `GT` is deep enough that it does not enter these comparisons.)
fn abs_err(f: f64, truth: &HpBounded, gt: usize) -> f64 {
    bf_mag(&BigFloat::from_f64(f, gt).sub(&truth.value, gt, HP_RM).abs())
}

/// A frame's realized basis `[origin, û, v̂, ŵ]`, probed the way `frame_world_basis` probes —
/// apply to the world origin and the three unit points, subtract.
fn probe_through(f: &FrameThrough) -> [[f64; 3]; 4] {
    let ap = |c: [f64; 3]| {
        exact(c)
            .expect("probe coords are f64")
            .frame_through(f)
            .expect("a constructed node realizes")
            .coord()
    };
    let o = ap([0.0, 0.0, 0.0]);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    [
        o,
        sub(ap([1.0, 0.0, 0.0]), o),
        sub(ap([0.0, 1.0, 0.0]), o),
        sub(ap([0.0, 0.0, 1.0]), o),
    ]
}

fn probe_named(pf: nacre_scalar::PlaneFrame) -> [[f64; 3]; 4] {
    let ap = |c: [f64; 3]| {
        exact(c)
            .expect("probe coords are f64")
            .frame(pf)
            .expect("a named frame realizes")
            .coord()
    };
    let o = ap([0.0, 0.0, 0.0]);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    [
        o,
        sub(ap([1.0, 0.0, 0.0]), o),
        sub(ap([0.0, 1.0, 0.0]), o),
        sub(ap([0.0, 0.0, 1.0]), o),
    ]
}

/// ★ **The two-road differential that locks the judged frame** (`three_planes_big`'s
/// pattern): a plane both roads can describe — world-rational points, so the narrow road has
/// exact canonical coefficients — must get the same basis from `plane_frame_default` +
/// `plane_frame_named` (exact route) and from [`FrameThrough`] (judged route), within the
/// judged route's own stated error. The fixture's raw normal is canonical already (first
/// coefficient positive, content 1), so `flip: false` compares like for like.
#[test]
fn the_judged_frame_agrees_with_the_named_road() {
    // (0,0,1), (1,0,0), (0,1,1): n = (1,0,1), d = −1 — canonical letter for letter.
    let pts = || {
        [
            exact([0.0, 0.0, 1.0]).unwrap(),
            exact([1.0, 0.0, 0.0]).unwrap(),
            exact([0.0, 1.0, 1.0]).unwrap(),
        ]
    };
    let c = [ri(1, 1), ri(0, 1), ri(1, 1), ri(-1, 1)];
    let (origin, ref_dir) = nacre_scalar::plane_frame_default(c).expect("derivable");
    let pf = nacre_scalar::plane_frame_named(c, origin, ref_dir).expect("narrow frame");
    let narrow = probe_named(pf);

    let ft = FrameThrough::of(pts().map(JudgedPoint::Pure), false).expect("a healthy plane frames");
    assert!(!ft.vertical, "n = (1,0,1) is nowhere near vertical");
    let judged = probe_through(&ft);

    for (row, (n, j)) in narrow.iter().zip(&judged).enumerate() {
        for k in 0..3 {
            assert!(
                (n[k] - j[k]).abs() < 1e-9,
                "row {row} component {k}: narrow {} vs judged {}",
                n[k],
                j[k]
            );
        }
    }

    // ★ The negative control: `flip` negates all four coefficients, which must flip ŵ and û
    // together and leave v̂ — `Frame`'s own documented semantics, reproduced by the judge.
    let flipped = probe_through(
        &FrameThrough::of(pts().map(JudgedPoint::Pure), true).expect("flip frames too"),
    );
    for k in 0..3 {
        assert!((flipped[1][k] + judged[1][k]).abs() < 1e-9, "û flips");
        assert!((flipped[2][k] - judged[2][k]).abs() < 1e-9, "v̂ stays");
        assert!((flipped[3][k] + judged[3][k]).abs() < 1e-9, "ŵ flips");
    }
}

/// The vertical branch: `ẑ×n` vanishes **exactly** for a horizontal plane (exact inputs make
/// the interval a true zero), so the constructor must take `ŷ×n` — and agree with the narrow
/// road, whose branch condition is the exact `n₀ = n₁ = 0`.
#[test]
fn the_judged_frame_takes_the_vertical_branch_where_the_named_road_does() {
    let pts = [
        exact([0.0, 0.0, 5.0]).unwrap(),
        exact([1.0, 0.0, 5.0]).unwrap(),
        exact([0.0, 1.0, 5.0]).unwrap(),
    ];
    let ft = FrameThrough::of(pts.map(JudgedPoint::Pure), false).expect("z = 5 frames");
    assert!(ft.vertical, "a horizontal plane must take ŷ×n");
    let judged = probe_through(&ft);

    let c = [ri(0, 1), ri(0, 1), ri(1, 1), ri(-5, 1)];
    let (origin, ref_dir) = nacre_scalar::plane_frame_default(c).expect("derivable");
    let pf = nacre_scalar::plane_frame_named(c, origin, ref_dir).expect("narrow frame");
    let narrow = probe_named(pf);
    for (row, (n, j)) in narrow.iter().zip(&judged).enumerate() {
        for k in 0..3 {
            assert!(
                (n[k] - j[k]).abs() < 1e-9,
                "row {row} component {k}: narrow {} vs judged {}",
                n[k],
                j[k]
            );
        }
    }
}

/// A carrier triple whose meet is the rational point `(a, b, c)` — the axis planes
/// `x = a`, `y = b`, `z = c`, each witnessed by three exact points, each optionally turned
/// by `spin` about Z. The cheapest constructible meet whose affine value is *also* statable
/// as a `Pure` point, which is what the two-road differential needs.
fn axis_carriers(a: f64, b: f64, c: f64, spin: Option<i128>) -> [[WitnessPoint; 3]; 3] {
    let zero3 = [ri(0, 1), ri(0, 1), ri(0, 1)];
    let turn = |w: WitnessPoint| match spin {
        None => w,
        Some(deg) => w.rotate_about(Axis::Z, Angle::from_deg(ri(deg, 1)).unwrap(), zero3),
    };
    let e = |x: f64, y: f64, z: f64| turn(exact([x, y, z]).unwrap());
    [
        [e(a, 0.0, 0.0), e(a, 1.0, 0.0), e(a, 0.0, 1.0)], // x = a
        [e(0.0, b, 0.0), e(1.0, b, 0.0), e(0.0, b, 1.0)], // y = b
        [e(0.0, 0.0, c), e(1.0, 0.0, c), e(0.0, 1.0, c)], // z = c
    ]
}

/// ★★★★★ **The Pure-vs-Meet differential — the tooth of the meet route** (`three_planes_big`'s
/// pattern). A pure point can *also* be written as the meet of its three carriers, so the
/// same three vertices spelled `[Pure; 3]` and `[Meet; 3]` must derive one basis — that is
/// the only proof that the affine shortcut (the all-pure road) and the homogeneous join
/// describe the same plane. The rotated variant is the hard
/// half: the carriers turn as planes, the pure spelling turns as a point, and the two meet
/// only if the whole chain of machinery — carrier realization, Cramer, join, normalization —
/// is right.
#[test]
fn a_meet_spelling_derives_the_same_frame_as_the_pure_one() {
    let zero3 = [ri(0, 1), ri(0, 1), ri(0, 1)];
    for spin in [None, Some(37)] {
        let pure_pt = |a: i128, b: i128, c: i128| {
            let w = WitnessPoint::at([ri(a, 1), ri(b, 1), ri(c, 1)]);
            JudgedPoint::Pure(match spin {
                None => w,
                Some(d) => w.rotate_about(Axis::Z, Angle::from_deg(ri(d, 1)).unwrap(), zero3),
            })
        };
        let meet_pt = |a: i128, b: i128, c: i128| {
            JudgedPoint::Meet(Box::new(axis_carriers(a as f64, b as f64, c as f64, spin)))
        };
        // Three non-collinear points, none axis-degenerate.
        let coords = [(1, 2, 3), (4, 1, 2), (2, 5, 7)];
        let pure = FrameThrough::of(coords.map(|(a, b, c)| pure_pt(a, b, c)), false)
            .expect("the pure spelling frames");
        let meet = FrameThrough::of(coords.map(|(a, b, c)| meet_pt(a, b, c)), false)
            .expect("the meet spelling frames");
        assert_eq!(pure.vertical, meet.vertical, "same branch (spin {spin:?})");
        let (bp, bm) = (probe_through(&pure), probe_through(&meet));
        for (row, (p, m)) in bp.iter().zip(&bm).enumerate() {
            for k in 0..3 {
                assert!(
                    (p[k] - m[k]).abs() < 1e-9,
                    "spin {spin:?}, row {row}, component {k}: pure {} vs meet {}",
                    p[k],
                    m[k]
                );
            }
        }
        // And the anchor — the one place a meet is divided — agrees with the pure cache.
        let (ap, am) = (
            pure.anchor_coord().expect("pure anchor"),
            meet.anchor_coord().expect("meet anchor"),
        );
        for k in 0..3 {
            assert!(
                (ap[k] - am[k]).abs() < 1e-9,
                "anchor component {k}: {} vs {}",
                ap[k],
                am[k]
            );
        }
    }
}

/// A genuinely straddling definition — carriers that do **not** share a chain within one
/// meet — realizes, and its image under the frame carries honest positive tol that bounds
/// the high-precision realization. This is the population no coordinate can state.
#[test]
fn a_straddling_meet_realizes_within_its_stated_tol() {
    const GT: usize = 384;
    let zero3 = [ri(0, 1), ri(0, 1), ri(0, 1)];
    let deg = |d: i128| Angle::from_deg(ri(d, 1)).unwrap();
    // One carrier turned 37° about Z, the other two still: the meet straddles frames.
    let e = |x: f64, y: f64, z: f64| exact([x, y, z]).unwrap();
    let straddle = JudgedPoint::Meet(Box::new([
        [
            e(1.0, 0.0, 0.0).rotate_about(Axis::Z, deg(37), zero3),
            e(1.0, 1.0, 0.0).rotate_about(Axis::Z, deg(37), zero3),
            e(1.0, 0.0, 1.0).rotate_about(Axis::Z, deg(37), zero3),
        ],
        [e(0.0, 2.0, 0.0), e(1.0, 2.0, 0.0), e(0.0, 2.0, 1.0)], // y = 2, still
        [e(0.0, 0.0, 3.0), e(1.0, 0.0, 3.0), e(0.0, 1.0, 3.0)], // z = 3, still
    ]));
    let pts = [
        straddle,
        JudgedPoint::Pure(e(5.0, 1.0, 0.0)),
        JudgedPoint::Pure(e(0.0, 6.0, 1.0)),
    ];
    let ft = FrameThrough::of(pts.clone(), false).expect("a straddling meet frames");
    let q = WitnessPoint::at([ri(1, 2), ri(1, 3), ri(2, 1)])
        .frame_through(&ft)
        .expect("realizes");
    let hp = q.hp_coord(GT);
    for (k, h) in hp.iter().enumerate() {
        let err = abs_err(q.realized[k].value, h, GT);
        assert!(
            err <= q.realized[k].error.max(1e-300),
            "axis {k}: cache off by {err:e}, stated tol {:e}",
            q.realized[k].error
        );
    }
    // Determinism — same statement, same node, same realization (replay's requirement).
    let ft2 = FrameThrough::of(pts, false).expect("again");
    assert_eq!(ft, ft2);
    let q2 = WitnessPoint::at([ri(1, 2), ri(1, 3), ri(2, 1)])
        .frame_through(&ft2)
        .expect("realizes");
    assert_eq!(q.coord(), q2.coord());
    assert_eq!(q.tol(), q2.tol());
}

/// A meet whose carriers all pass through one line has no unique point — `D` straddles zero
/// at every precision, the join has nothing to normalize against, and the constructor must
/// refuse (the producer then rejects by name, without claiming degeneracy).
#[test]
fn a_meet_whose_carriers_share_a_line_is_refused() {
    let e = |x: f64, y: f64, z: f64| exact([x, y, z]).unwrap();
    // Three planes through the z-axis: x = 0, y = 0, x = y — D exactly 0.
    let sheaf = JudgedPoint::Meet(Box::new([
        [e(0.0, 0.0, 0.0), e(0.0, 1.0, 0.0), e(0.0, 0.0, 1.0)],
        [e(0.0, 0.0, 0.0), e(1.0, 0.0, 0.0), e(0.0, 0.0, 1.0)],
        [e(0.0, 0.0, 0.0), e(1.0, 1.0, 0.0), e(0.0, 0.0, 1.0)],
    ]));
    let pts = [
        sheaf,
        JudgedPoint::Pure(e(5.0, 0.0, 0.0)),
        JudgedPoint::Pure(e(0.0, 5.0, 0.0)),
    ];
    assert!(
        FrameThrough::of(pts, false).is_none(),
        "a meet with no unique point must not frame"
    );
}

/// A definition that cannot prove a basis is `None` from the constructor — collinear points
/// have a normal that is exactly zero, so neither branch's squared length clears zero. The
/// producer's duty (reject **by name**, and not as a proof of collinearity) is ops'; here
/// the contract is only that no frame comes back.
#[test]
fn a_definition_that_cannot_prove_a_basis_is_refused() {
    let pts = [
        exact([0.0, 0.0, 0.0]).unwrap(),
        exact([1.0, 1.0, 1.0]).unwrap(),
        exact([2.0, 2.0, 2.0]).unwrap(),
    ];
    assert!(
        FrameThrough::of(pts.map(JudgedPoint::Pure), false).is_none(),
        "collinear points name no plane and must not frame"
    );
}

/// ★ The point of the whole node: three defining points with **different** chains — the
/// population no exact frame can serve. The f64 cache must sit within its own stated tol of
/// the high-precision realization, and construction must be deterministic (same statement,
/// same branch, bit-identical realization — replay's requirement).
#[test]
fn a_heterogeneous_definition_realizes_within_its_stated_tol() {
    const GT: usize = 384;
    let deg = |d: i128| Angle::from_deg(ri(d, 1)).unwrap();
    let zero3 = [ri(0, 1), ri(0, 1), ri(0, 1)];
    let pts = || {
        [
            WitnessPoint::at([ri(2, 1), ri(0, 1), ri(0, 1)]).rotate_about(Axis::Z, deg(37), zero3),
            exact([5.0, 1.0, 0.0]).unwrap(),
            WitnessPoint::at([ri(0, 1), ri(3, 1), ri(1, 1)]).rotate_about(Axis::X, deg(22), zero3),
        ]
    };
    let ft =
        FrameThrough::of(pts().map(JudgedPoint::Pure), false).expect("a mixed-chain plane frames");
    let q = WitnessPoint::at([ri(1, 2), ri(1, 3), ri(2, 1)])
        .frame_through(&ft)
        .expect("realizes");
    let hp = q.hp_coord(GT);
    for (k, h) in hp.iter().enumerate() {
        let err = abs_err(q.realized[k].value, h, GT);
        assert!(
            err <= q.realized[k].error.max(1e-300),
            "axis {k}: cache off the realization by {err:e}, stated tol {:e}",
            q.realized[k].error
        );
        assert!(
            q.realized[k].error > 0.0,
            "a judged frame's image carries honest positive tol"
        );
    }

    // Determinism: the same statement builds the same node (definitional equality — the
    // `PartialEq` `shared_base` relies on) and realizes bit-identically.
    let ft2 = FrameThrough::of(pts().map(JudgedPoint::Pure), false).expect("again");
    assert_eq!(ft, ft2, "same statement, same node");
    let q2 = WitnessPoint::at([ri(1, 2), ri(1, 3), ri(2, 1)])
        .frame_through(&ft2)
        .expect("realizes");
    assert_eq!(q.coord(), q2.coord(), "same statement, same realization");
    assert_eq!(q.tol(), q2.tol(), "and the same stated error");
}

/// **The shared-motion shortcut must answer the same question a reflection is in the chain.**
///
/// `shared_base` cancels one shared motion out of `det[a−d, b−d, c−d]` and answers on the
/// pre-motion coordinates, *exactly* — on the strength of "a rigid motion preserves this
/// determinant". A **reflection does not**: it negates it. And the failure is not conservative,
/// because all four points carry the same chain, so the flip is uniform and the shortcut
/// returns a confidently wrong `Sign`.
///
/// So the claim is measured rather than argued: random chains that mix rotations, translations
/// and reflections, against a 512-bit realization of the same four definitions.
///
/// `#[ignore]`: astro-float ground truth is slow.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn the_shared_motion_shortcut_agrees_with_ground_truth() {
    const GT: usize = 512;
    let mut st = 0x51D2_7E44_0C13_A001u64;
    let (mut took_shortcut, mut mirror_seen) = (0usize, false);
    for _ in 0..2000 {
        // Four points sharing one chain — the shortcut's precondition. **Dyadic** bases, or
        // `shared_base` declines before the sign is ever in question (its own guard is that a
        // base must round-trip through f64) and the test would pass vacuously.
        let dyadic = |st: &mut u64| -> [Rat; 3] {
            std::array::from_fn(|_| ri(rng(st, -100_000, 100_000), 1 << rng(st, 0, 6)))
        };
        let bases: [[Rat; 3]; 4] = std::array::from_fn(|_| dyadic(&mut st));
        let mut pts: Vec<WitnessPoint> = bases.iter().map(|&b| WitnessPoint::at(b)).collect();
        for _ in 0..rng(&mut st, 1, 4) {
            match rng(&mut st, 0, 2) {
                0 => {
                    let off = dyadic(&mut st);
                    pts = pts.into_iter().map(|p| p.translate(off)).collect();
                }
                1 => {
                    mirror_seen = true;
                    let ax = axis_of(rng(&mut st, 0, 2));
                    let off = if rng(&mut st, 0, 2) == 0 {
                        ri(0, 1)
                    } else {
                        dyadic(&mut st)[0]
                    };
                    pts = pts.into_iter().map(|p| p.mirror(ax, off)).collect();
                }
                _ => {
                    let ax = axis_of(rng(&mut st, 0, 2));
                    let ang = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 997));
                    let piv = dyadic(&mut st);
                    pts = pts
                        .into_iter()
                        .map(|p| p.rotate_about(ax, ang, piv))
                        .collect();
                }
            }
        }
        let [pa, pb, pc, pd] = [&pts[0], &pts[1], &pts[2], &pts[3]];
        let Some(b) = shared_base(&[pa, pb, pc, pd]) else {
            continue; // a base that does not round-trip; the toleranced path answers
        };
        took_shortcut += 1;
        let shortcut = match nacre_predicates::orient3d(b[0], b[1], b[2], b[3]) {
            x if x > 0.0 => 1i8,
            x if x < 0.0 => -1,
            _ => 0,
        };
        // Ground truth: the same determinant, realized from the definitions at 512 bits. Deep
        // enough that only a genuinely near-zero case is ambiguous, and those are skipped.
        let det = det3_hp(pa, pb, pc, pd, GT);
        let Some(pos) = det.sign() else {
            continue;
        };
        let truth = if pos { 1i8 } else { -1 };
        assert_eq!(
            shortcut, truth,
            "the shortcut answered a different question than the definitions do"
        );
    }
    assert!(mirror_seen, "the sample must include a reflection");
    assert!(
        took_shortcut > 100,
        "the shortcut must actually fire ({took_shortcut} times)"
    );
}

/// Soundness: over random **motion** chains (mixed axes, arbitrary pivots, exact and inexact
/// angles, **and rational translations and reflections interleaved**), the direction-wise tol
/// must bound the true f64 error on every axis (astro-float 512-bit ground truth) — the
/// production mirror of exact3d H-d/H-f. (An `err` below 1e-100 is 512-bit GT noise.)
///
/// **Every node kind must appear here.** Each arrives with its own tol term, and nothing else
/// in the suite checks that term is an upper bound — the whole judgment layer is sound only if
/// it is. (The translate term was once added without this, and had to be back-filled.)
/// `#[ignore]`: astro-float ground truth is slow; run with `--ignored` (+ CI). The
/// per-predicate statistical validation is the `h_*` tests beside this one.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn tol_bounds_error_over_random_chains() {
    const GT: usize = 512;
    let mut st = 0x2A5C_1234_ABCD_9999u64;
    let (mut exact_seen, mut pivot_seen) = (false, false);
    let (mut translate_seen, mut mirror_seen) = (false, false);
    let mut worst_ratio = 0.0_f64;
    let mut worst_rot = 0.0_f64;
    // ★ **What the chain was made of, so a failure can be acted on.** A gate that only says
    // "the bound was exceeded" cannot be debugged: the fix for an exceeded bound is to find
    // the *missing term*, and that means knowing which links the offending chain had. The
    // per-shape bests below also say where the bound is already tight and where it is not —
    // which is what decides whether a given term is worth deriving more carefully.
    #[derive(Default, Clone)]
    struct Shape {
        piv: usize,
        tr: usize,
        mir: usize,
        desc: String,
    }
    // `(worst ratio, count, the chain that reached it)` per shape.
    let mut by_shape: Vec<(&str, f64, usize, String)> = vec![
        ("origin pivot only     ", 0.0, 0, String::new()),
        ("has a non-origin pivot", 0.0, 0, String::new()),
        ("no mirror             ", 0.0, 0, String::new()),
        ("has a mirror          ", 0.0, 0, String::new()),
        ("no translate          ", 0.0, 0, String::new()),
        ("has a translate       ", 0.0, 0, String::new()),
        // ★ **Pivot rotations and nothing else.** Added to isolate the `piv` term, and it
        // showed something else instead: this row sits at **exactly 0.5000** and does not
        // move when `piv` changes, because what dominates there is the *base seed* — `WitnessPoint::at`
        // measures the base's own rounding and charges twice it, so a chain that only rotates
        // reports half its bound and nothing else can shift that. **A row that cannot move is
        // reporting a different term than the one it was built for.**
        ("pivot rotations only  ", 0.0, 0, String::new()),
    ];
    // ★ And the **mean**, because a term that shrinks everywhere by a little moves no maximum.
    let (mut ratio_sum, mut ratio_n) = (0.0_f64, 0usize);
    let mut worst_desc = String::new();
    for _ in 0..2000 {
        let base = rand_base(&mut st);
        let seed_free = WitnessPoint::at(base).tol() == [0.0; 3];
        let mut shape = Shape::default();
        let mut p = WitnessPoint::at(base);
        // **From one node, not two.** A single origin rotation of an exact base is the case
        // the deleted 2D frame validated on its own; sampling it here is what makes this test
        // strictly cover that one, rather than merely resemble it.
        for _ in 0..rng(&mut st, 1, 5) {
            // Links of every kind, in both orders — which is where a tol term that only holds
            // "on its own" would show.
            match rng(&mut st, 0, 3) {
                0 => {
                    translate_seen = true;
                    shape.tr += 1;
                    shape.desc.push_str(" T");
                    p = p.translate(rand_base(&mut st));
                    continue;
                }
                1 => {
                    mirror_seen = true;
                    let off = if rng(&mut st, 0, 2) == 0 {
                        Rat::from_int(0) // the exact case: a pure sign flip
                    } else {
                        rand_base(&mut st)[0]
                    };
                    shape.mir += 1;
                    shape
                        .desc
                        .push_str(if off == Rat::from_int(0) { " M0" } else { " M" });
                    p = p.mirror(axis_of(rng(&mut st, 0, 2)), off);
                    continue;
                }
                _ => {}
            }
            let ax = axis_of(rng(&mut st, 0, 2));
            let pivot = match rng(&mut st, 0, 2) {
                0 => [Rat::from_int(0); 3],
                1 => {
                    pivot_seen = true;
                    rand_base(&mut st)
                }
                _ => {
                    pivot_seen = true;
                    base // point on the pivot
                }
            };
            let ang = if rng(&mut st, 0, 3) == 0 {
                exact_seen = true;
                deg(90 * rng(&mut st, 0, 3), 1)
            } else {
                deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973))
            };
            let at_origin = pivot.iter().all(|r| *r == Rat::from_int(0));
            if !at_origin {
                shape.piv += 1;
            }
            shape.desc.push_str(if at_origin { " R" } else { " Rp" });
            if ang.try_exact_cos_sin().is_some() {
                shape.desc.push('!'); // an exact angle: the rotation term is zero
            }
            p = p.rotate_about(ax, ang, pivot);
        }
        let hp = p.hp_coord(GT);
        for (axis, hp_a) in hp.iter().enumerate() {
            let err = abs_err(p.realized[axis].value, hp_a, GT);
            // ★ **The chain belongs in the failure, not only in the summary below.** An
            // exceeded bound means a term is missing, and the first thing needed is which
            // links were involved — but the panic aborts before any summary prints, so the
            // shape has to travel with the message. It is how the asymmetry in `rot` was
            // found: the offending chain was `R M`, and nothing else was.
            assert!(
                err <= p.realized[axis].error || err < 1e-100,
                "tol must bound the error: axis {axis}, err {err:e} > tol {:e}, links:{}",
                p.realized[axis].error,
                shape.desc
            );
            if p.realized[axis].error > 0.0 {
                let r = err / p.realized[axis].error;
                ratio_sum += r;
                ratio_n += 1;
                if r > worst_ratio {
                    worst_ratio = r;
                    worst_desc = format!("axis {axis}, links:{}", shape.desc);
                }
                if seed_free {
                    worst_rot = worst_rot.max(r);
                }
                for (i, hit) in [
                    shape.piv == 0,
                    shape.piv > 0,
                    shape.mir == 0,
                    shape.mir > 0,
                    shape.tr == 0,
                    shape.tr > 0,
                    shape.piv > 0 && shape.tr == 0 && shape.mir == 0,
                ]
                .into_iter()
                .enumerate()
                {
                    if hit && r > by_shape[i].1 {
                        by_shape[i].1 = r;
                        by_shape[i].3 = shape.desc.clone();
                    }
                }
            }
        }
        for (i, hit) in [
            shape.piv == 0,
            shape.piv > 0,
            shape.mir == 0,
            shape.mir > 0,
            shape.tr == 0,
            shape.tr > 0,
            shape.piv > 0 && shape.tr == 0 && shape.mir == 0,
        ]
        .into_iter()
        .enumerate()
        {
            if hit {
                by_shape[i].2 += 1;
            }
        }
    }
    assert!(
        exact_seen && pivot_seen && translate_seen && mirror_seen,
        "corpus must mix exact angles, pivots, translations and reflections"
    );
    // **How much of the bound the real error actually uses.** Every term feeding `tol` is now
    // either derived from round-to-nearest or measured against arbitrary precision — no
    // constant is left to justify — so this is what says the *assembly* of them holds.
    //
    // ★ It stays a figure nobody writes down: a number pasted into a doc went stale the moment
    // reflections joined this corpus, and was wrong for two commits.
    eprintln!("[tol tightness] worst err/tol: {worst_ratio:.4} all, {worst_rot:.4} rotation-only");
    // ★ **Which chain got closest, and which shapes are already tight.** The fix for an
    // exceeded bound is to find the *missing term*, so the gate has to say what the offending
    // chain was made of. The per-shape bests say the same thing in advance: a shape sitting
    // near 1.0 is a term already derived tightly, and one sitting low is where a lumped charge
    // still has slack. `R` a rotation about the origin, `Rp` about a pivot, `!` an exact angle,
    // `T` a translation, `M` a reflection (`M0` a pure sign flip).
    eprintln!(
        "[tol worst chain]  {worst_ratio:.4}  {worst_desc}   mean {:.4} over {ratio_n}",
        ratio_sum / ratio_n.max(1) as f64
    );
    for (label, r, n, desc) in &by_shape {
        eprintln!("[tol by shape]  {label}  worst {r:.4}  over {n:>4} chains  ←{desc}");
    }
    assert!(
        worst_ratio < 1.0,
        "the bound was reached exactly, which leaves nothing for a platform whose trig is \
             one ulp worse than this one's"
    );
}

/// A 90°-family chain about the origin stays exactly tol 0 (Niven exact realization).
/// **Why the tol is a coordinate *mixing* term and not a tangential one.**
///
/// The design's first directional tol was tangential — `tol_x = da·|y|`, `tol_y = da·|x|` —
/// on the reasoning that rotation moves a point along its circle. It is **not sound**: `cos`
/// and `sin` round independently, which puts a *radial* component into the error, and near an
/// axis the tangential prediction goes to zero while the real error does not. The adopted
/// `(|bx|+|by|)·da` per component covers both.
///
/// This is the standing guard on that shape. It moved here when the 2D frame it was written
/// in was deleted (nothing consumed that frame), because the *invariant* is about the tol
/// formula, which 3D uses verbatim — a rotation about Z is the 2D case with `z` carried along.
#[test]
fn the_tangential_tol_is_unsound_and_the_mixing_one_is_not() {
    const GT: usize = 256;
    let (mut tangential_sound, mut mixing_sound) = (true, true);
    let mut first_break = None;
    for (bxn, bxd, byn, byd) in [(1, 1, 0, 1), (3, 1, 4, 1), (5, 2, 7, 3), (1, 1, 1, 1)] {
        // Near-axis angles are where the tangential model is thinnest, so they are in here.
        for (an, ad) in [
            (30, 1),
            (45, 1),
            (60, 1),
            (37, 1),
            (899, 10),
            (1, 10),
            (3, 1),
        ] {
            let base = [ri(bxn, bxd), ri(byn, byd), ri(0, 1)];
            let angle = deg(an, ad);
            let p = WitnessPoint::at(base).rotate(Axis::Z, angle);
            let hp = p.compute_hp(GT);
            // The realization error this angle actually has — the same quantity `rotate_about`
            // charges. Using it makes the refutation stronger than a constant would: even the
            // *measured* trig error, applied tangentially, fails to bound.
            let (c, s) = angle.cos_sin_f64();
            let (dc, ds) = angle.realization_error_of(c, s);
            for (k, hp_k) in hp.iter().enumerate().take(2) {
                let err = abs_err(p.realized[k].value, hp_k, GT);
                // (1) the refuted tangential prediction: the *other* coordinate's magnitude.
                if err > [dc, ds][k] * p.realized[1 - k].value.abs() {
                    tangential_sound = false;
                    first_break.get_or_insert((bxn, bxd, byn, byd, an, ad, k));
                }
                // (2) the adopted mixing tol, which is what `rotate_about` computes.
                if err > p.realized[k].error {
                    mixing_sound = false;
                }
            }
        }
    }
    assert!(
        !tangential_sound,
        "the tangential formula bounded every sample — either the fixtures stopped reaching \
             near an axis, or the realization changed and this guard is now vacuous"
    );
    assert!(
        mixing_sound,
        "the adopted tol failed to bound a single origin rotation: {first_break:?}"
    );
}

/// Planes a face-based operation actually lands on: the world XY, a raised cap, a wall, the
/// Pythagorean `(3,4,0)` slope a user really draws, and two genuinely irrational tilts.
const FRAME_PLANES: [[i128; 4]; 6] = [
    [0, 0, 1, 0],
    [0, 0, 1, -5],
    [1, 0, 0, -2],
    [3, 4, 0, -10],
    [1, 1, 1, -3],
    [2, -3, 7, 11],
];

fn frame_of(c: [i128; 4]) -> nacre_scalar::PlaneFrame {
    let coeffs = c.map(Rat::from_int);
    nacre_scalar::plane_frame(coeffs).expect("these planes all have exact frames")
}

/// ★★★★ **The property a frame exists to have: a point drawn at `w = 0` is *on* the plane.**
///
/// Checked on the high-precision realization, against the plane's own exact coefficients — so
/// it catches a wrong origin, a `u` that is not in the plane, and a `v` that is not
/// perpendicular to both, none of which a "does it run" test would notice. The projection step
/// this design skips (`(n·n)·ref − (ref·n)·n`) is skipped on the strength of exactly this
/// being true, so it is the claim that has to be measured rather than argued.
#[test]
fn a_point_drawn_in_a_frame_lies_on_that_frame_s_plane() {
    const GT: usize = 512;
    for c in FRAME_PLANES {
        let fr = frame_of(c);
        for (u, v) in [(0, 0), (1, 0), (0, 1), (3, -7), (-2, 5)] {
            let base = [ri(u, 1), ri(v, 1), Rat::from_int(0)];
            let p = WitnessPoint::at(base)
                .frame(fr)
                .expect("a frame with positive lengths");
            let hp = p.hp_coord(GT);
            // a·x + b·y + c·z + d, realized — its interval must contain zero.
            let mut e = HpBounded::new(BigFloat::from_i128(c[3], GT), Mag::ZERO);
            for k in 0..3 {
                e = e.add(
                    &hp[k].mul(&HpBounded::of_rat(Rat::from_int(c[k]), GT), GT),
                    GT,
                );
            }
            let mag = bf_mag(&e.value.abs());
            let error = e.error.exp2().map_or(0.0, |x| 2f64.powi(x as i32));
            assert!(
                mag <= error || mag < 1e-100,
                "plane {c:?}, point ({u}, {v}): off the plane by {mag:e}, radius {error:e}"
            );
        }
    }
}

/// ★★★★ **A named frame is a frame the derivation cannot produce.**
///
/// `plane_frame_default` gives the convention a *face* takes — origin at the world origin's
/// projection, `+u` along `ẑ × n`. A plane a caller **named** may want neither: the script
/// layer's `ZX` plane has `+u = +ẑ`, while `ẑ × n` there is `−x̂`. So the node carries the
/// pair, and this checks that naming them actually moves the frame — and that the frame is
/// still orthonormal when it does.
///
/// ★★ **And two spellings of one direction are one frame.** `ref_dir` is reduced to its
/// primitive form, so `[0,0,1]` and `[0,0,5]` build the identical `PlaneFrame` — which is what
/// lets a motion node intern (the document's blocker 6).
#[test]
fn a_named_frame_takes_its_own_origin_and_u() {
    let coeffs = [0, 1, 0, 0].map(Rat::from_int); // the ZX plane, y = 0
    let (o_d, r_d) = nacre_scalar::plane_frame_default(coeffs).unwrap();
    assert_eq!(
        r_d.map(|r| r.to_f64()),
        [-1.0, 0.0, 0.0],
        "ẑ × n is −x̂ here"
    );

    let zhat = [0, 0, 1].map(Rat::from_int);
    let named = nacre_scalar::plane_frame_named(coeffs, o_d, zhat).unwrap();
    let derived = nacre_scalar::plane_frame(coeffs).unwrap();
    assert_ne!(
        named.u_raw, derived.u_raw,
        "naming +u must actually move it"
    );

    // The named frame is still a frame: `+u` realizes to `ẑ`, and a point at `w = 0` is on
    // the plane `y = 0`.
    let p = WitnessPoint::at([ri(3, 1), ri(-7, 1), Rat::from_int(0)])
        .frame(named)
        .unwrap();
    let q = WitnessPoint::at([Rat::from_int(0); 3])
        .frame(named)
        .unwrap();
    let u = [0, 1, 2].map(|k| p.realized[k].value - q.realized[k].value);
    assert!(p.realized[1].value.abs() < 1e-15, "on the plane y = 0");
    assert!((u[2] - 3.0).abs() < 1e-15, "+u ran along ẑ, got {u:?}");

    // ★ Two spellings of one direction, one frame.
    let long = [0, 0, 5].map(Rat::from_int);
    assert_eq!(
        nacre_scalar::plane_frame_named(coeffs, o_d, long),
        Some(named)
    );
}

/// ★★★ **The `v̂` fallback is taken, and the frame is still a frame.**
///
/// `|v_raw|² = |n|²·|u_raw|²` needs twice the width the lengths themselves do, so a plane with
/// wide coefficients cannot carry it — and the realization takes `v̂ = ŵ × û` instead. That
/// route is what this crate had before `v_raw` existed; what must not change is that the basis
/// is still a basis. Checked by the property that matters: a point drawn at `w = 0` is *on*
/// the plane, and `tol` still bounds the coordinate it wrote.
///
/// ★ **A fallback nothing takes is indistinguishable from one that is not there**, so the
/// first assertion is that these planes really do decline the exact route.
#[test]
fn a_frame_whose_v_does_not_fit_falls_back_and_still_works() {
    const GT: usize = 512;
    // Components near 2⁴⁰: the squared lengths fit, their product does not.
    const WIDE: [[i128; 4]; 3] = [
        [1_099_511_627_776, 1_099_511_627_791, 1_099_511_627_803, -7],
        [999_999_999_989, -1_000_000_000_039, 1_000_000_000_061, 13],
        [2_199_023_255_552, 3_298_534_883_329, -1_099_511_627_777, 0],
    ];
    for c in WIDE {
        let fr = frame_of(c);
        assert!(
            fr.v.is_none(),
            "plane {c:?} was expected to decline the exact v̂"
        );
        for (u, v) in [(0, 0), (1, 0), (0, 1), (3, -7)] {
            let p = WitnessPoint::at([ri(u, 1), ri(v, 1), Rat::from_int(0)])
                .frame(fr)
                .expect("a frame with positive lengths");
            let hp = p.hp_coord(GT);
            let mut e = HpBounded::new(BigFloat::from_i128(c[3], GT), Mag::ZERO);
            for k in 0..3 {
                e = e.add(
                    &hp[k].mul(&HpBounded::of_rat(Rat::from_int(c[k]), GT), GT),
                    GT,
                );
            }
            let mag = bf_mag(&e.value.abs());
            let error = e.error.exp2().map_or(0.0, |x| 2f64.powi(x as i32));
            assert!(
                mag <= error || mag < 1e-100,
                "plane {c:?}, point ({u}, {v}): off the plane by {mag:e}, radius {error:e}"
            );
            for (k, h) in hp.iter().enumerate() {
                let err = abs_err(p.realized[k].value, h, GT);
                assert!(
                    err <= p.realized[k].error || err < 1e-100,
                    "plane {c:?} axis {k}: err {err:e} > tol {:e}",
                    p.realized[k].error
                );
            }
        }
    }
    // …and a narrow plane still takes the exact route, so the decline above is the width's
    // doing and not something that turned the exact `v̂` off everywhere.
    assert!(frame_of([2, -3, 7, 11]).v.is_some());
}

/// **`tol` bounds the error in the `coord` the frame actually wrote** — the same contract
/// `a_rotation_chain_s_tol_bounds_its_realization` holds a rotation to, for the one motion
/// whose realization is a square root rather than a cosine.
#[test]
fn a_frame_s_tol_bounds_its_own_realization() {
    const GT: usize = 512;
    let (mut worst, mut worst_at) = (0.0f64, String::new());
    for c in FRAME_PLANES {
        let fr = frame_of(c);
        for (u, v, w) in [(0, 0, 0), (1, 0, 0), (3, -7, 2), (-2, 5, -1), (11, 13, 17)] {
            let base = [ri(u, 1), ri(v, 1), ri(w, 1)];
            let p = WitnessPoint::at(base).frame(fr).unwrap();
            let hp = p.hp_coord(GT);
            for (k, h) in hp.iter().enumerate() {
                let err = abs_err(p.realized[k].value, h, GT);
                assert!(
                    err <= p.realized[k].error || err < 1e-100,
                    "plane {c:?} at ({u},{v},{w}) axis {k}: err {err:e} > tol {:e}",
                    p.realized[k].error
                );
                if p.realized[k].error > 0.0 && err / p.realized[k].error > worst {
                    worst = err / p.realized[k].error;
                    worst_at = format!("plane {c:?} at ({u},{v},{w}) axis {k}");
                }
            }
        }
    }
    eprintln!("[cip] worst frame tol usage: {worst:.3} of the bound ({worst_at})");
    // A bound nothing ever approaches is a bound nobody derived — the rotation tests hold
    // themselves to the same reading.
    assert!(worst > 1e-6, "the frame tol is never approached: {worst:e}");
}

/// The wide triple's plane — a canonical answer past `i128` — as its exact
/// arbitrary-precision coefficients and the canonical [`WideFrame`] built on them.
/// The width is **asserted**, not assumed: qualify the fixture.
fn wide_fixture() -> ([num_bigint::BigInt; 4], WideFrame) {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let a = [q(big1, 3), q(big2, 7), q(0, 1)];
    let b = [q(-big2, 5), q(big1, 11), q(0, 1)];
    let c = [q(1, 13), q(1, 17), q(1, 19)];
    let name = nacre_scalar::plane_name_exact(a, b, c).expect("not collinear");
    assert!(
        name.narrow().is_none(),
        "the fixture must be genuinely wide"
    );
    let nacre_scalar::PlaneName::Wide(cs) = name else {
        unreachable!("narrow() said Wide")
    };
    let n = [cs[0].clone(), cs[1].clone(), cs[2].clone()];
    let fr = WideFrame::canonical(n, &cs[3]).expect("a nonzero normal");
    (cs, fr)
}

/// ★★★★ **The wide twin of `a_point_drawn_in_a_frame_lies_on_that_frame_s_plane`**:
/// a point drawn at `w = 0` in a [`WideFrame`] is *on* that plane, checked on the
/// high-precision realization against the exact arbitrary-precision coefficients. This is
/// the property that lets a sketch land on a plane whose data exceeds `i128` at all.
#[test]
fn a_point_drawn_in_a_wide_frame_lies_on_that_plane() {
    const GT: usize = 512;
    let (c, fr) = wide_fixture();
    for (u, v) in [(0, 0), (1, 0), (0, 1), (3, -7), (-2, 5)] {
        let base = [ri(u, 1), ri(v, 1), Rat::from_int(0)];
        let p = WitnessPoint::at(base)
            .frame_wide(&fr)
            .expect("a wide frame has positive lengths");
        let hp = p.hp_coord(GT);
        // a·x + b·y + c·z + d, realized — its interval must contain zero.
        let mut e = HpBounded::exact(nacre_scalar::bigint_to_bigfloat(&c[3], GT));
        for k in 0..3 {
            e = e.add(&hp[k].mul(&HpBounded::of_bigint(&c[k], GT), GT), GT);
        }
        let mag = bf_mag(&e.value.abs());
        let error = e.error.exp2().map_or(0.0, |x| 2f64.powi(x as i32));
        assert!(
            mag <= error,
            "wide plane, point ({u}, {v}): off the plane by {mag:e}, radius {error:e}"
        );
    }
}

/// **The wide twin of `a_frame_s_tol_bounds_its_own_realization`**: the f64 `coord`
/// a wide frame writes is within the `tol` it writes, measured against a far deeper
/// realization of the same definition.
#[test]
fn a_wide_frame_s_tol_bounds_its_own_realization() {
    const GT: usize = 512;
    let (_, fr) = wide_fixture();
    let (mut worst, mut worst_at) = (0.0f64, String::new());
    for (u, v, w) in [(0, 0, 0), (1, 0, 0), (3, -7, 2), (-2, 5, -1), (11, 13, 17)] {
        let base = [ri(u, 1), ri(v, 1), ri(w, 1)];
        let p = WitnessPoint::at(base).frame_wide(&fr).unwrap();
        let hp = p.hp_coord(GT);
        for (k, h) in hp.iter().enumerate() {
            let err = abs_err(p.realized[k].value, h, GT);
            assert!(
                err <= p.realized[k].error || err < 1e-100,
                "wide frame at ({u},{v},{w}) axis {k}: err {err:e} > tol {:e}",
                p.realized[k].error
            );
            if p.realized[k].error > 0.0 && err / p.realized[k].error > worst {
                worst = err / p.realized[k].error;
                worst_at = format!("({u},{v},{w}) axis {k}");
            }
        }
    }
    eprintln!("[cip] worst wide-frame tol usage: {worst:.3} of the bound ({worst_at})");
    assert!(
        worst > 1e-6,
        "the wide-frame tol is never approached: {worst:e}"
    );
}

/// ★★★ **A frame on an axis-aligned plane through the origin costs nothing at all** — the
/// basis is a signed permutation, so every product and every sum is exact and `tol` stays a
/// literal zero. That is what keeps such a sketch on the exact predicate path, and it is the
/// case the `perm` gate in [`WitnessPoint::frame`] exists for.
#[test]
fn an_axis_aligned_frame_through_the_origin_is_tol_zero() {
    for c in [[0i128, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0]] {
        let fr = frame_of(c);
        for (u, v, w) in [(0, 0, 0), (3, -7, 2), (11, 13, 17)] {
            let p = WitnessPoint::at([ri(u, 1), ri(v, 1), ri(w, 1)])
                .frame(fr)
                .unwrap();
            assert_eq!(p.tol(), [0.0; 3], "plane {c:?} at ({u},{v},{w})");
        }
    }
    // …and the world XY frame is the identity, which is what makes a frame on it harmless.
    let fr = frame_of([0, 0, 1, 0]);
    let p = WitnessPoint::at([ri(3, 1), ri(-7, 1), ri(2, 1)])
        .frame(fr)
        .unwrap();
    assert_eq!(p.coord(), [3.0, -7.0, 2.0]);
}

/// ★★★★ **The prize: two points sketched in one frame are judged *exactly*.**
///
/// `shared_base` cancels a motion every input carries and hands the exact predicate the
/// pre-motion coordinates. A frame is a rigid motion like any other, so a whole sketch on one
/// tilted face cancels down to its own rational `(u, v, w)` — and this asks for the answer only
/// an exact route can give: a **proved zero**. The f64 filter can prove a sign is nonzero; it
/// can never prove a determinant *is* zero, so `Some(Orient::Zero)` here is the cancellation
/// firing and nothing else.
///
/// ★★ **`shared_base` needed no change for this**, and that is worth pinning rather than
/// leaving to luck: it was written against a *property* (all chains structurally equal, and
/// the motion preserves the determinant) rather than against a list of variants, so a new
/// proper motion joined it for free. A test is what tells "it is right" from "it happens to
/// be right".
#[test]
fn four_points_sharing_one_frame_are_judged_exactly() {
    let fr = frame_of([2, -3, 7, 11]);
    let at = |u: i128, v: i128, w: i128| {
        WitnessPoint::at([ri(u, 1), ri(v, 1), ri(w, 1)])
            .frame(fr)
            .unwrap()
    };
    // Four points on the sketch plane itself (`w = 0`) — exactly coplanar, by construction.
    let (a, b, c, d) = (at(0, 0, 0), at(1, 0, 0), at(0, 1, 0), at(2, 3, 0));
    assert_eq!(
        orient3d_filter(&a, &b, &c, &d),
        Some(Orient::Zero),
        "a proved zero needs the exact route; the f64 filter cannot reach one"
    );
    // And a point off the plane still gets a sign, so the shortcut is not answering `Zero`
    // to everything.
    let e = at(2, 3, 1);
    assert!(matches!(
        orient3d_filter(&a, &b, &c, &e),
        Some(Orient::Positive | Orient::Negative)
    ));
}

/// ★★ **A frame is proper, so it must not touch the mirror parity** — the design says so
/// (`v̂ = ŵ × û` makes the basis right-handed) and `chain_parity` counts only reflections, so
/// this holds today by not being written. Pinned because the consequence of it changing is a
/// determinant handed back with the wrong sign, which is a confident wrong answer rather than
/// a slow one.
#[test]
fn a_frame_does_not_flip_the_chain_parity() {
    let fr = frame_of([3, 4, 0, -10]);
    let framed = WitnessPoint::at([ri(1, 1), ri(2, 1), ri(3, 1)])
        .frame(fr)
        .unwrap();
    assert_eq!(chain_parity(&framed.chain), 1);
    let mirrored = framed.clone().mirror(Axis::X, Rat::from_int(0));
    assert_eq!(
        chain_parity(&mirrored.chain),
        -1,
        "the reflection still counts"
    );
    let fr2 = frame_of([1, 1, 1, -3]);
    let twice = mirrored.frame(fr2).unwrap();
    assert_eq!(
        chain_parity(&twice.chain),
        -1,
        "a second frame leaves the odd reflection odd"
    );
}

#[test]
fn quadrantal_origin_chain_is_tol_zero() {
    let p = WitnessPoint::at([ri(3, 1), ri(5, 1), ri(7, 1)])
        .rotate(Axis::Z, deg(90, 1))
        .rotate(Axis::X, deg(180, 1))
        .rotate(Axis::Y, deg(270, 1));
    assert_eq!(p.tol(), [0.0; 3]);
    // and the realized coords are the exact permutation/negation (no spurious term).
    assert_eq!(p.coord(), [7.0, -3.0, -5.0]);
}

/// **A representable base measures exactly zero — `at` and `at_nearest` agree on it.**
///
/// `at_nearest` states tol `0` for a base every coordinate of which is an `f64`, without
/// measuring; this pins that the measurement would have said the same. It rests on three
/// links, and the third is the one worth a test: `try_from_f64` represents an f64 exactly,
/// realizing that base back yields the same f64, and `bf_mag` of an exact zero is `0.0` (not
/// a floor). If any link broke, `at` would report a nonzero tol here and the stated zero
/// would be silently changing geometry rather than skipping arithmetic.
#[test]
fn a_representable_base_measures_exactly_zero() {
    for c in [
        [0.0, 1.0, -1.0],         // integers, both signs
        [0.5, 0.25, -0.125],      // dyadic fractions
        [4.0, 0.2, 3.0],          // 0.2 is not decimal-exact but *is* an exact f64
        [1e-8, -2.5e-9, 7.5e-7],  // small
        [1e18, -4e17, 3.5e19],    // large, still inside Rat's exponent range
        [1e-20, -2.5e-21, 5e-18], // small, still inside it (the floor is 2^-74 ≈ 5.3e-23)
    ] {
        let fast = exact(c).expect("representable");
        let measured = WitnessPoint::at(c.map(|x| Rat::try_from_f64(x).expect("representable")));
        assert_eq!(fast.coord(), c, "the coordinates round-trip: {c:?}");
        assert_eq!(fast.coord(), measured.coord(), "same coord for {c:?}");
        assert_eq!(
            measured.tol(),
            [0.0; 3],
            "`at` must measure exactly zero for an f64-derived base: {c:?}"
        );
        assert_eq!(fast.tol(), measured.tol(), "same tol for {c:?}");
    }
}

/// **`at_nearest` states the tol from `Rat::to_f64`'s contract: zero where the coordinate is an
/// f64, `|x|·2⁻⁵³` where it is not — and never below what `at` measures.**
#[test]
fn at_nearest_states_zero_for_an_f64_and_the_contract_bound_otherwise() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // Representable throughout: identical to `exact`, and to what `at` measures.
    let base = [r(1, 2), r(-3, 1), r(1, 8)];
    let p = WitnessPoint::at_nearest(base);
    assert_eq!(p.tol(), [0.0; 3]);
    assert_eq!(p.coord(), [0.5, -3.0, 0.125]);
    assert_eq!(p.tol(), WitnessPoint::at(base).tol());
    // One coordinate `f64` cannot hold: the contract bound on every axis, never below the
    // measured rounding, and the coordinate is the nearest f64 of the base itself.
    let base = [r(3, 10), r(1, 1), r(7, 10)];
    let p = WitnessPoint::at_nearest(base);
    let measured = WitnessPoint::at(base);
    assert_eq!(p.coord(), measured.coord());
    assert_eq!(p.coord(), [0.3, 1.0, 0.7]);
    assert!(
        measured.tol()[0] > 0.0,
        "0.3 is not an f64 — the measurement must see it"
    );
    for k in 0..3 {
        let x = p.coord()[k];
        let want = (x.abs() * (f64::EPSILON * 0.5)).max(f64::MIN_POSITIVE);
        assert_eq!(p.tol()[k], want, "axis {k}");
        // The contract: `to_f64` is nearest, so the rounding is at most half an ulp of the
        // value it produced — and the stated bound covers that. (`at`'s measurement is
        // *inflated* by 2 for conservatism, so it is not the yardstick here.)
        let half_ulp = 0.5 * (f64::from_bits(x.to_bits() + 1) - x);
        assert!(
            p.tol()[k] >= half_ulp,
            "the contract bound covers half an ulp on axis {k}: {} < {half_ulp}",
            p.tol()[k]
        );
    }
}

/// **An exactly-representable point demands no precision at all** — its trial bound is `0`.
///
/// The twin of the test above, one level up: that one pins the `f64` tol, this one pins the
/// *realization* bound, which is what [`judge_precision`] reads to size a model. So a model
/// with no rotation history asks for nothing, and a caller that already knows a point is
/// `Constructed` can skip [`trial_bound`] rather than spend a realization computing a zero.
///
/// Not an accident of small numbers. The fixture builds its base with `Rat::try_from_f64` =
/// `mantissa · 2^exp`, so the **denominator is a power of two** and the **numerator fits
/// `i128`** — and `rat_to_hp` returns an exact interval exactly when those two hold at `prec`.
/// The corpus therefore reaches **both ends of `Rat`'s range**: where the numerator is widest
/// (a large integer, ~127 bits) and where the denominator is deepest (`2^126`).
///
/// ★ The coupling this rests on — `TRIAL_PREC ≥ 127` — is asserted at the constant itself,
/// where lowering it fails the build rather than one test.
#[test]
fn an_exact_point_demands_no_precision() {
    for c in [
        [0.0, 1.0, -1.0],         // integers, both signs
        [0.5, 0.25, -0.125],      // dyadic fractions
        [4.0, 0.2, 3.0],          // 0.2 is not decimal-exact but *is* an exact f64
        [1e18, -4e17, 3.5e19],    // large
        [1e38, -1.5e38, 1.0],     // ★ near Rat's ceiling: the widest numerator, ~127 bits
        [1e-20, -2.5e-21, 5e-18], // small
        [1e-22, -2e-22, 1.0],     // ★ near Rat's floor (2^-74): the deepest denominator
    ] {
        let p = exact(c).expect("representable");
        assert!(
            trial_bound(&p).is_zero(),
            "an exact point must realize exactly: {c:?} gave {:?}",
            trial_bound(&p)
        );
    }
}

/// Outside `Rat`'s exponent range there is no exact base, and `Rat::try_from_f64` says so
/// instead of panicking — the caller decides (an operation turns it into a named reject).
///
/// **The range is much narrower than f64's, at both ends** — easy to get wrong, and I did
/// on the first attempt. `Rat` is `Ratio<i128>`, so `mantissa · 2^exp` must fit `i128`
/// (`|x| ≲ 1.7e38`) *and* the denominator `2^k` must (`k = 1075 − exp_field ≤ 126`, i.e.
/// `|x| ≳ 2^-74 ≈ 5.3e-23`). Exact zero is special-cased and always representable.
/// A CAD model at either extreme is not real; the limit is.
#[test]
fn a_coordinate_outside_rats_range_has_no_exact_base() {
    assert!(exact([1e300, 0.0, 0.0]).is_none(), "too large");
    assert!(exact([1e-30, 0.0, 0.0]).is_none(), "below the 2^-74 floor");
    assert!(exact([f64::MIN_POSITIVE, 0.0, 0.0]).is_none(), "subnormal");
    // Zero is not a boundary case — it is special-cased and exact.
    assert_eq!(exact([0.0; 3]).expect("zero is exact").tol(), [0.0; 3]);
}

/// `at` seeds the base→f64 rounding: 0 for an integer base, positive for a base
/// that is not f64-representable (e.g. 1/3).
#[test]
fn at_seeds_base_rounding_tol() {
    assert_eq!(
        WitnessPoint::at([ri(2, 1), ri(3, 1), ri(4, 1)]).tol(),
        [0.0; 3]
    );
    let third = WitnessPoint::at([ri(1, 3), ri(0, 1), ri(0, 1)]);
    assert!(third.realized[0].error > 0.0 && third.realized[1].error == 0.0);
}

/// `at_with_tol` seeds a nonzero root tol (a Discovered seam) and the chain
/// transports it (`|R|·old`): a 90° rotation about Z swaps the x/y tol components.
#[test]
fn seeded_tol_transports_through_rotation() {
    let p = WitnessPoint::at_with_tol([ri(1, 1), ri(0, 1), ri(0, 1)], [1e-9, 2e-9, 3e-9])
        .rotate(Axis::Z, deg(90, 1));
    // 90° about Z: |R| swaps x,y → tol[0]=old tol[1], tol[1]=old tol[0]; z unchanged.
    assert_eq!(p.tol(), [2e-9, 1e-9, 3e-9]);
}

/// Same-axis bundling is tighter than incremental (H-e): K steps of θ amplify the
/// tol (each transported by |c|+|s| ≥ 1) vs one step of Kθ.
#[test]
fn bundling_is_tighter_than_incremental() {
    let base = [ri(11, 1), ri(-7, 1), ri(4, 1)];
    let (k, theta) = (30i128, ri(1, 7));
    let mut incr = WitnessPoint::at(base);
    for _ in 0..k {
        incr = incr.rotate(Axis::Z, Angle::from_deg(theta).unwrap());
    }
    let k_theta = theta.checked_mul(Rat::from_int(k)).unwrap();
    let bundled = WitnessPoint::at(base).rotate(Axis::Z, Angle::from_deg(k_theta).unwrap());
    assert!(
        bundled.realized[0].error < incr.realized[0].error
            && bundled.realized[1].error < incr.realized[1].error,
        "bundled {:?} must be tighter than incremental {:?}",
        bundled.tol(),
        incr.tol()
    );
}

// ---- orient3d_judge (2a-ii) ----

/// A point rotated about one random axis (through a random pivot) by one **inexact**
/// rational angle — generic non-degenerate, heterogeneous provenance (each point its
/// own rotation), the H-a corpus shape. Inexact angles only (an exact/near-coplanar
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

/// H-a — `det3_bound` soundness: over many random heterogeneous-rotation 4-point
/// configs (inexact angles, arbitrary pivots), the bound must upper-bound the real
/// error of the f64 determinant vs the astro-float truth — never exceeded. The
/// production mirror of exact3d H-a. `#[ignore]`: slow (4× astro-float per sample).
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn h_a_det3_bound_soundness() {
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
    eprintln!("[H-a N=5000] det3_bound violations: {bad}; worst tightness: {worst:.3}");
    assert_eq!(
        bad, 0,
        "det3_bound must bound the determinant error on every sample"
    );
}

/// Fast port check (default suite): a handful of random configs must not violate
/// `det3_bound` (catches a transcription bug; the full statistical soundness is the
/// `#[ignore]`d H-a + exact3d).
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

// ---- indirect orient3d (2c-i) ----

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

/// The high-precision indirect truth with a stability flag: `None` if even the
/// ground truth cannot resolve `D` or `M` above its floor (a genuine degeneracy).
#[allow(clippy::too_many_arguments)]
/// Ground truth for the indirect judge — and **not** by running the judge harder.
///
/// It used to call `sign_with_floor` with the judge's own `mag`, differing only in precision.
/// A floor that is wrong is then wrong identically in both, so the two agree and the test
/// passes: the oracle shared the defect it existed to find, which is how the `mag`-collapse
/// bug survived a soundness suite. Here the verdict comes from **comparing two precisions**
/// instead: a sign is trusted only when `prec` and `2·prec` produce the same nonzero sign, and
/// anything else is `None` (not asserted against). That borrows no formula from the judge.
fn indirect_truth(
    pa: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    pb: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    pc: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    q: &WitnessPoint,
    r: &WitnessPoint,
    s: &WitnessPoint,
    prec: usize,
) -> Option<Orient> {
    let at = |prec: usize| -> Orient {
        let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
        let (d, m, _gap) = indirect_hp(
            [ph(pa), ph(pb), ph(pc)],
            q.hp_coord(prec),
            r.hp_coord(prec),
            s.hp_coord(prec),
            prec,
        );
        // Raw signs, no floor: the agreement between two precisions is what filters noise.
        let raw = |x: &BigFloat| {
            if x.is_zero() {
                None
            } else {
                Some(x.is_positive())
            }
        };
        combine(raw(&d.value), raw(&m.value)).unwrap_or(Orient::Zero)
    };
    let (lo, hi) = (at(prec), at(2 * prec));
    (lo == hi && lo != Orient::Zero).then_some(lo)
}

/// H-b — rotated plane coefficient tol soundness. A plane's four coefficients
/// derive from three rotated points (subtraction/cross/dot); the interval `error` on
/// each must upper-bound the real f64 error vs the astro-float truth. Corpus mixes
/// origin and arbitrary pivots (`rand_point`). `#[ignore]`: slow astro-float GT.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn h_b_plane_coefficient_tol_soundness() {
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
    eprintln!("[H-b N={N}] coefficient tol violations: {bad}; worst tightness: {worst:.3}");
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

/// H-c — indirect orient3d soundness over heterogeneous provenance (corpus A) and a
/// near-coplanar escalation-forcing corpus (corpus B). The interval filter →
/// astro-float judge must never claim a sign opposite to a GT-stable 512-bit truth.
/// Asserts both paths are live: `escalated > 0` (filter defers) **and**
/// `filter_resolved > 0` (fast path resolves — a filter stuck at `None` would pass
/// wrong-sign 0 vacuously). `#[ignore]`: slow astro-float GT.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn h_c_indirect_orient3d_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let mut check = |a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
                     q: &WitnessPoint,
                     r: &WitnessPoint,
                     s: &WitnessPoint| {
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
        match indirect_truth(a, b, c, q, r, s, GT) {
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
            );
        }
    }

    eprintln!(
        "[H-c] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
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
/// statistical soundness is the `#[ignore]`d H-b/H-c + exact3d).
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
        if let Some(truth) = indirect_truth(a, b, c, &p[9], &p[10], &p[11], GT) {
            assert!(
                judged == truth || judged == Orient::Zero,
                "indirect judge {judged:?} disagrees with truth {truth:?}"
            );
        }
    }
}

// ---- indirect cmp_coord (cmp-i) ----

fn add3(a: [Rat; 3], b: [Rat; 3]) -> [Rat; 3] {
    [
        a[0].checked_add(b[0]).unwrap(),
        a[1].checked_add(b[1]).unwrap(),
        a[2].checked_add(b[2]).unwrap(),
    ]
}

/// Borrow an owned plane-triple as the `&WitnessPoint` tuples the judge takes.
fn tr(t: &[[WitnessPoint; 3]; 3]) -> [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3] {
    [
        (&t[0][0], &t[0][1], &t[0][2]),
        (&t[1][0], &t[1][1], &t[1][2]),
        (&t[2][0], &t[2][1], &t[2][2]),
    ]
}

/// Three planes meeting at `v` (each through `v` + two small offsets), all rotated by
/// `(ax, ang, piv)` — an implicit point at `rotate(v)` with heterogeneous provenance.
fn triple_pts(
    v: [Rat; 3],
    st: &mut u64,
    ax: Axis,
    ang: Angle,
    piv: [Rat; 3],
) -> [[WitnessPoint; 3]; 3] {
    let plane = |st: &mut u64| {
        let off = |st: &mut u64| {
            [
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
                ri(rng(st, -20, 20), rng(st, 1, 5)),
            ]
        };
        let (o1, o2) = (off(st), off(st));
        [
            WitnessPoint::at(v).rotate_about(ax, ang, piv),
            WitnessPoint::at(add3(v, o1)).rotate_about(ax, ang, piv),
            WitnessPoint::at(add3(v, o2)).rotate_about(ax, ang, piv),
        ]
    };
    [plane(st), plane(st), plane(st)]
}

/// A random rotated three-plane triple (own random center, axis, inexact angle,
/// pivot) — heterogeneous provenance. Sequences the RNG draws so each `&mut st`
/// borrow ends before the next.
fn rand_triple(st: &mut u64) -> [[WitnessPoint; 3]; 3] {
    let v = rand_base(st);
    let ax = axis_of(rng(st, 0, 2));
    let angle = deg(rng(st, 0, 360_000), rng(st, 1, 9973));
    let piv = rand_base(st);
    triple_pts(v, st, ax, angle, piv)
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

/// The high-precision cmp truth (`None` when the coordinates are equal or below the
/// GT floor — a genuine tie).
fn cmp_truth(
    a: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    b: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3],
    axis: usize,
    prec: usize,
) -> Option<Orient> {
    let hp = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
        [
            plane_hp(t[0].0, t[0].1, t[0].2, prec),
            plane_hp(t[1].0, t[1].1, t[1].2, prec),
            plane_hp(t[2].0, t[2].1, t[2].2, prec),
        ]
    };
    cmp_hp_with_gap(hp(a), hp(b), axis, prec).ok()
}

/// `cmp_combine` maps the parity of negative signs to the ordering.
#[test]
fn cmp_combine_counts_negatives() {
    assert_eq!(
        cmp_combine(Some(true), Some(true), Some(true)),
        Some(Orient::Positive)
    );
    assert_eq!(
        cmp_combine(Some(false), Some(true), Some(true)),
        Some(Orient::Negative)
    );
    assert_eq!(
        cmp_combine(Some(false), Some(false), Some(true)),
        Some(Orient::Positive)
    );
    assert_eq!(cmp_combine(None, Some(true), Some(true)), None);
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
/// rotated configs (catches a transcription bug; full soundness is `#[ignore]`d H-g).
#[test]
fn cmp_port_check() {
    const GT: usize = 384;
    let mut st = 0xC301_7A5E_2266_9911u64;
    for _ in 0..16 {
        let a = rand_triple(&mut st);
        let b = rand_triple(&mut st);
        let axis = rng(&mut st, 0, 2) as usize;
        let judged = indirect_cmp_coord_judge(tr(&a), tr(&b), axis, fixture()).orient();
        if let Some(truth) = cmp_truth(tr(&a), tr(&b), axis, GT) {
            assert!(
                judged == truth || judged == Orient::Zero,
                "cmp judge {judged:?} disagrees with truth {truth:?}"
            );
        }
    }
}

/// H-g — indirect cmp_coord soundness over heterogeneous provenance (corpus A) and a
/// near-tie corpus (corpus B: two points sharing an axis coordinate up to ε, rotated
/// about that same axis so the near-tie survives). The judge must never disagree with
/// a GT-stable 512-bit truth; both the fast filter and the escalation are exercised.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn h_g_indirect_cmp_coord_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
    let mut check = |a: &[[WitnessPoint; 3]; 3], b: &[[WitnessPoint; 3]; 3], axis: usize| {
        let (ta, tb) = (tr(a), tr(b));
        let judged = indirect_cmp_coord_judge(ta, tb, axis, fixture()).orient();
        let iv = |t: [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3]| {
            [
                plane_iv(t[0].0, t[0].1, t[0].2),
                plane_iv(t[1].0, t[1].1, t[1].2),
                plane_iv(t[2].0, t[2].1, t[2].2),
            ]
        };
        if cmp_filter(iv(ta), iv(tb), axis).is_none() {
            escalated += 1;
        } else {
            filter_resolved += 1;
        }
        match cmp_truth(ta, tb, axis, GT) {
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

    // Corpus A — heterogeneous provenance (each triple its own rotation and pivot).
    let mut st = 0x6A11_C0DE_5151_2323u64;
    for _ in 0..2000 {
        let a = rand_triple(&mut st);
        let b = rand_triple(&mut st);
        check(&a, &b, rng(&mut st, 0, 2) as usize);
    }

    // Corpus B — near-tie in z, shared Z-rotation (a Z-rotation leaves z unchanged, so
    // the pre-rotation z near-equality survives → M near 0 → escalation forced).
    for _ in 0..2000 {
        let z0 = rng(&mut st, -200, 200);
        let eps = match rng(&mut st, 0, 2) {
            0 => ri(0, 1),
            1 => ri(1, rng(&mut st, 5_000, 200_000)),
            _ => ri(1, rng(&mut st, 1_000_000_000, 1_000_000_000_000_000)),
        };
        let va = [
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
            ri(z0, 1),
        ];
        let vb = [
            ri(rng(&mut st, -200, 200), 1),
            ri(rng(&mut st, -200, 200), 1),
            ri(z0, 1).checked_add(eps).unwrap(),
        ];
        let angle = deg(rng(&mut st, 0, 360_000), rng(&mut st, 1, 9973));
        let piv = rand_base(&mut st);
        let a = triple_pts(va, &mut st, Axis::Z, angle, piv);
        let b = triple_pts(vb, &mut st, Axis::Z, angle, piv);
        check(&a, &b, 2);
    }

    eprintln!(
        "[H-g] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
    );
    assert_eq!(
        wrong, 0,
        "cmp judge must never disagree with GT (soundness)"
    );
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path"
    );
}

// ---- dir_sign_judge (3a-iii) ----

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

/// The **raw** sign of `D` (det of the three normals) at `prec` bits — no radius, no limit.
fn dir_d_sign_at(
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    prec: usize,
) -> Option<bool> {
    let ph = |t: (&WitnessPoint, &WitnessPoint, &WitnessPoint)| plane_hp(t.0, t.1, t.2, prec);
    // The **raw** sign, with no radius and no floor — see `dir_orient_at` for why the oracle
    // must not borrow the judge's bound. `dir_sign_truth` gets its confidence from two
    // precisions agreeing instead.
    let (dh, _) = cramer_hp(&[ph(a), ph(b), ph(c)], prec);
    (!dh.value.is_zero()).then(|| dh.value.is_positive())
}

/// The **GT-stable** `dir_sign` truth: `Some` only when `prec` and `prec + 128` agree
/// on a definite sign — otherwise the config is degenerate beyond what the ground
/// truth itself resolves, so it is `None` (skip), not a spurious "wrong". (Mirrors
/// exact3d's `indirect_truth` stability check.)
fn dir_sign_truth(
    a: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    b: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    c: (&WitnessPoint, &WitnessPoint, &WitnessPoint),
    prec: usize,
) -> Option<Orient> {
    match (
        dir_d_sign_at(a, b, c, prec),
        dir_d_sign_at(a, b, c, prec + 128),
    ) {
        (Some(x), Some(y)) if x == y => Some(orient_of(x)),
        _ => None,
    }
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

/// H-i — dir_sign soundness over a **near-coplanar-normals** corpus (which H-c/H-g do
/// not stress: they force `M ≈ 0`, not `D ≈ 0`). Three plane normals `n0, n1,
/// n2 = α·n0 + β·n1 + ε·(n0×n1)` (ε tiny → `D ≈ ε` → escalation), shared rotation. The
/// judge must never disagree with a GT-stable 512-bit `D` sign; both paths exercised.
#[test]
#[ignore = "slow astro-float ground truth (run with --ignored)"]
fn h_i_dir_sign_soundness() {
    const GT: usize = 512;
    let (mut wrong, mut declined, mut escalated, mut filter_resolved, mut skipped, mut tested) =
        (0usize, 0, 0, 0, 0, 0);
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
        match dir_sign_truth(t3(&pa), t3(&pb), t3(&pc), GT) {
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
        "[H-i] wrong: {wrong}/{tested}; declined: {declined}; escalated: {escalated}; filter_resolved: {filter_resolved}; skipped: {skipped}"
    );
    assert_eq!(wrong, 0, "dir_sign must never disagree with GT (soundness)");
    assert!(escalated > 0, "corpus must exercise the escalation path");
    assert!(
        filter_resolved > 0,
        "corpus must exercise the fast filter path"
    );
}

/// The `dir_orient3d` determinant `det[d, x−base, y−base]` at `prec` bits — the GT /
/// escalation realization (mirrors the judge's hp path).
fn dir_orient_at(
    d: [Rat; 3],
    base: &WitnessPoint,
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Option<bool> {
    let dp = WitnessPoint::at(d);
    let sub = |u: &HpBounded, v: &HpBounded| u.sub(v, prec);
    let (bh, xh, yh, dh) = (
        base.hp_coord(prec),
        x.hp_coord(prec),
        y.hp_coord(prec),
        dp.hp_coord(prec),
    );
    let rows = [
        dh,
        [
            sub(&xh[0], &bh[0]),
            sub(&xh[1], &bh[1]),
            sub(&xh[2], &bh[2]),
        ],
        [
            sub(&yh[0], &bh[0]),
            sub(&yh[1], &bh[1]),
            sub(&yh[2], &bh[2]),
        ],
    ];
    let det = det3_big_rows(&rows, prec);
    // The **raw** sign at `prec` bits, with no radius and no floor. Independence from the
    // judge is the whole point of an oracle: `dir_orient_truth` gets its confidence from two
    // precisions agreeing, not from any bound this file also ships to production.
    (!det.value.is_zero()).then(|| det.value.is_positive())
}

/// GT-stable truth: `Some` only when `prec` and `prec + 128` agree (else too degenerate
/// for the ground truth itself — skip, not a spurious "wrong"). Mirrors `dir_sign_truth`.
fn dir_orient_truth(
    d: [Rat; 3],
    base: &WitnessPoint,
    x: &WitnessPoint,
    y: &WitnessPoint,
    prec: usize,
) -> Option<Orient> {
    match (
        dir_orient_at(d, base, x, y, prec),
        dir_orient_at(d, base, x, y, prec + 128),
    ) {
        (Some(a), Some(b)) if a == b => Some(orient_of(a)),
        _ => None,
    }
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
fn h_dir_orient3d_soundness() {
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
        match dir_orient_truth(d, &base, &x, &y, GT) {
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
