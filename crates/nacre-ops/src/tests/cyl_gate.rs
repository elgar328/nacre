//! Unit tests of the population gate's own functions — the separating directions, a lateral
//! face's reach, the tangent wall's row. (`tests/suite/cylinder_gate.rs` drives the gate end to
//! end through booleans.)

use super::*;
use nacre_exact::Rat;

/// **Three directions, and each one is the only one that answers**.
///
/// `A` stands on the `z` axis through the origin, radius `1/5`, `z ∈ [0, 2]`. `B` lies along
/// `x`, same radius — and where it is put decides which direction sees it apart. The fourth
/// row is the one that matters most: a `B` that really runs through `A` is separated by
/// **none** of the three, which is what keeps a genuine crossing refused.
#[test]
fn each_separating_direction_is_the_only_one_for_some_pair() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let cyl = |o: [Rat; 3], m: [Rat; 3], e: [Rat; 3]| {
        nacre_topo::CylinderDef::new(o, m, e, nacre_exact::BigRat::from(q(1, 25))).unwrap() // radius 1/5, as r²
    };
    let a = cyl([z(0); 3], [z(0), z(0), z(1)], [z(1), z(0), z(0)]);
    let fp = |span: [i128; 2]| Footprint {
        span: Some([z(span[0]), z(span[1])]),
        theta: Theta::Whole,
    };
    let (x, dirs_of) = (fp([0, 2]), |b: &nacre_topo::CylinderDef| {
        separating_dirs(&a, &fp([0, 2]), b, &fp([0, 2]))
    });
    // `B` along `x`: its origin's `y`/`z` and its own span place it.
    let along_x =
        |y: i128, zz: i128| cyl([z(0), z(y), z(zz)], [z(1), z(0), z(0)], [z(0), z(0), z(1)]);
    let cases: [(nacre_topo::CylinderDef, [i128; 2], [bool; 3]); 4] = [
        // Clear of `A`'s span along `A`'s own axis, and nothing else.
        (along_x(0, 5), [-10, 10], [true, false, false]),
        // Off the end of `B`'s own span, and nothing else.
        (along_x(0, 1), [5, 10], [false, true, false]),
        // Beside it, across both axes — only the common perpendicular sees this.
        (along_x(5, 1), [-10, 10], [false, false, true]),
        // Straight through it: no direction separates them, and none may.
        (along_x(0, 1), [-10, 10], [false, false, false]),
    ];
    for (b, span, want) in cases {
        let dirs = dirs_of(&b);
        assert_eq!(dirs.len(), 3, "skew axes offer three directions");
        let got: Vec<bool> = dirs
            .iter()
            .map(|d| separated(&a, &x, &b, &fp(span), d) == Some(true))
            .collect();
        assert_eq!(got, want.to_vec(), "b at {:?} span {span:?}", b.origin());
    }
    // Parallel axes have no third direction to offer.
    let parallel = cyl([z(1), z(0), z(0)], [z(0), z(0), z(1)], [z(1), z(0), z(0)]);
    assert_eq!(separating_dirs(&a, &x, &parallel, &x).len(), 2);
    // A line on a cylinder names one more: its own radial direction.
    let line = Footprint {
        span: None,
        theta: Theta::At([z(1), z(0), z(0)]),
    };
    assert_eq!(
        separating_dirs(&a, &line, &parallel, &x),
        vec![a.dir(), parallel.dir(), [z(1), z(0), z(0)]]
    );
}

/// **Each of [`pieces_clear`]'s three witnesses is the only one that answers for some pair**, and a
/// pair that touches has none. Hand geometry, `[separated along some direction, one piece inside the
/// other's cylinder, parallel cross-section]` per pair:
/// - apart along `x`: a short crosswise cylinder beside a tall one;
/// - inside: a crosswise `r = 5` over `x ∈ [−5, 5]` within an `r = 20` about `z` — no direction
///   parts them;
/// - cross-section: two parallel unit cylinders `3/2` apart, each the half towards `+x`, whose
///   circles cross and whose halves do not;
/// - none: the crosswise `r = 3` over `x ∈ [−4, 4]` within `r = 5`, its rim on the outer lateral —
///   the box's corners land exactly on it, and strictly inside is what the witness asks.
#[test]
fn each_witness_of_pieces_clear_is_the_only_one_for_some_pair() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let cyl = |o: [Rat; 3], m: [Rat; 3], e: [Rat; 3], r2: Rat| {
        nacre_topo::CylinderDef::new(o, m, e, nacre_exact::BigRat::from(r2)).unwrap()
    };
    let (zx, xz) = (
        ([z(0), z(0), z(1)], [z(1), z(0), z(0)]),
        ([z(1), z(0), z(0)], [z(0), z(0), z(1)]),
    );
    let whole = |lo: i128, hi: i128| Footprint {
        span: Some([z(lo), z(hi)]),
        theta: Theta::Whole,
    };
    let right_half = Footprint {
        span: Some([z(0), z(1)]),
        theta: Theta::Arc(crate::planes::RimArc {
            from: [z(0), z(-1), z(0)],
            to: [z(0), z(1), z(0)],
        }),
    };
    let witnesses =
        |a: &nacre_topo::CylinderDef, x: &Footprint, b: &nacre_topo::CylinderDef, y: &Footprint| {
            [
                separating_dirs(a, x, b, y)
                    .iter()
                    .any(|d| separated(a, x, b, y, d) == Some(true)),
                piece_inside(a, x, b) || piece_inside(b, y, a),
                nacre_exact::parallel_rat(&a.dir(), &b.dir())
                    && cross_sections_clear(a, x, b, y) == Some(true),
            ]
        };
    let cases = [
        (
            cyl([z(0); 3], zx.0, zx.1, z(1)),
            whole(0, 3),
            cyl([z(2), z(0), q(1, 2)], xz.0, xz.1, q(1, 25)),
            whole(0, 1),
            [true, false, false],
        ),
        (
            cyl([z(0); 3], zx.0, zx.1, z(400)),
            whole(-20, 20),
            cyl([z(-5), z(0), z(0)], xz.0, xz.1, z(25)),
            whole(0, 10),
            [false, true, false],
        ),
        (
            cyl([z(0); 3], zx.0, zx.1, z(1)),
            right_half,
            cyl([q(3, 2), z(0), z(0)], zx.0, zx.1, z(1)),
            right_half,
            [false, false, true],
        ),
        (
            cyl([z(0); 3], zx.0, zx.1, z(25)),
            whole(-5, 5),
            cyl([z(-4), z(0), z(0)], xz.0, xz.1, z(9)),
            whole(0, 8),
            [false, false, false],
        ),
    ];
    for (a, x, b, y, want) in cases {
        let got = witnesses(&a, &x, &b, &y);
        assert_eq!(got, want, "a at {:?}, b at {:?}", a.origin(), b.origin());
        assert_eq!(pieces_clear(&a, &x, &b, &y), want.contains(&true));
    }
}

/// **The third direction is the face-level twin of the surface rung.** With a whole circle
/// and no span, `separated` along the common perpendicular says exactly what
/// [`nacre_exact::cylinders_clear`] says about the infinite surfaces — so the face test does
/// not disagree with the rung that runs before it, it only knows more when a face knows more.
#[test]
fn the_common_perpendicular_degenerates_to_the_surface_rule() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let cyl = |o: [Rat; 3], m: [Rat; 3], e: [Rat; 3]| {
        nacre_topo::CylinderDef::new(o, m, e, nacre_exact::BigRat::from(q(1, 25))).unwrap() // radius 1/5, as r²
    };
    let a = cyl([z(0); 3], [z(0), z(0), z(1)], [z(1), z(0), z(0)]);
    let unbounded = Footprint {
        span: None,
        theta: Theta::Whole,
    };
    // Offsets either side of the radius sum `2/5`, and exactly on it.
    for (n, d) in [(1i128, 2i128), (2, 5), (1, 3), (9, 10)] {
        let b = cyl(
            [z(0), q(n, d), z(0)],
            [z(1), z(0), z(0)],
            [z(0), z(0), z(1)],
        );
        let perp = separating_dirs(&a, &unbounded, &b, &unbounded)[2];
        let face = separated(&a, &unbounded, &b, &unbounded, &perp);
        let surface = nacre_exact::cylinders_clear(
            &a.origin(),
            &a.dir(),
            a.r2(),
            &b.origin(),
            &b.dir(),
            b.r2(),
        ) == nacre_exact::Orient::Positive;
        assert_eq!(face, Some(surface), "offset {n}/{d}");
    }
}

/// **The reach of a lateral face along a direction, and the clearance read off it** — hand
/// geometry, every arm of the two predicates.
///
/// The face: a cylinder on the `z` axis through the origin, radius `1/5`, span `t ∈ [0, 2]`
/// (`z ∈ [0, 2]`). Along `d = (0, 0, 1)` its reach is the span itself (`d·m ≠ 0`, no radial
/// part: `|d⊥| = 0`); along `d = (0, 1, 0)` it is `[−1/5, 1/5]` whatever the span — and with
/// no span at all, which is the only arm a production boolean reaches today.
#[test]
fn a_lateral_faces_reach_and_what_clears_it() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let def = nacre_topo::CylinderDef::new(
        [z(0); 3],
        [z(0), z(0), z(1)],
        [z(1), z(0), z(0)],
        nacre_exact::BigRat::from(q(1, 25)), // radius 1/5, as r²
    )
    .unwrap();
    let fp = |span: Option<[Rat; 2]>| Footprint {
        span,
        theta: Theta::Whole,
    };
    let span = Some([z(0), z(2)]);
    // A whole-circle reach is `lo − √ρ² .. hi + √ρ²` with one `ρ²` on both ends.
    let whole = |r: &Reach| {
        assert_eq!(
            r.rho2_lo, r.rho2_hi,
            "a whole circle reaches both ways alike"
        );
        (r.lo, r.hi, r.rho2_hi)
    };
    // Along its own axis: the projected span, no radial reach.
    let along = lateral_reach(&def, &fp(span), &[z(0), z(0), z(1)]).unwrap();
    assert_eq!(whole(&along), (z(0), z(2), z(0)));
    // Across it: a point widened by r² — with or without a span.
    for s in [span, None] {
        let across = lateral_reach(&def, &fp(s), &[z(0), z(1), z(0)]).unwrap();
        assert_eq!(whole(&across), (z(0), z(0), q(1, 25)));
    }
    // No span and a projection that needs one: nothing proved.
    assert!(lateral_reach(&def, &fp(None), &[z(0), z(0), z(1)]).is_none());
    // Slanted, `d = (0, 1, 1)`: `d·m = 1`, `|d⊥|² = 2 − 1 = 1` — span and radial part both.
    let slant = lateral_reach(&def, &fp(span), &[z(0), z(1), z(1)]).unwrap();
    assert_eq!(whole(&slant), (z(0), z(2), q(1, 25)));
    // Reversed, the projection's ends swap and the reach is still stated low end first.
    let back = lateral_reach(&def, &fp(span), &[z(0), z(0), z(-1)]).unwrap();
    assert_eq!(whole(&back), (z(-2), z(0), z(0)));
    // Overflow is `None`, never a verdict.
    assert!(lateral_reach(&def, &fp(span), &[q(i128::MAX / 2, 1); 3]).is_none());

    // Clearance against the across reach, `0 ± 1/5`:
    let reach = lateral_reach(&def, &fp(None), &[z(0), z(1), z(0)]).unwrap();
    assert_eq!(reach_clears(&reach, q(1, 4), z(1)), Some(true)); // past the far end
    assert_eq!(reach_clears(&reach, z(-1), q(-1, 4)), Some(true)); // past the near end
    assert_eq!(reach_clears(&reach, q(1, 5), z(1)), Some(false)); // touches the end: a point
    assert_eq!(reach_clears(&reach, q(-1, 10), q(1, 10)), Some(false)); // inside
    assert_eq!(reach_clears(&reach, q(-1, 10), z(1)), Some(false)); // straddles

    // **An arc reaches less than its circle.** The quarter arc from `(−r, 0)` to
    // `(0, −r)` — the third quadrant, counter-clockwise about `+z`.
    let m = [z(0), z(0), z(1)];
    let arc = RimArc {
        from: [q(-1, 5), z(0), z(0)],
        to: [z(0), q(-1, 5), z(0)],
    };
    assert_eq!(
        arc_contains(&arc, &rat_dir(&[z(-1), z(-1), z(0)]), &m),
        Some(true)
    ); // 225°
    assert_eq!(arc_contains(&arc, &rat_dir(&arc.from), &m), Some(true)); // an end is on it
    assert_eq!(
        arc_contains(&arc, &rat_dir(&[z(0), z(1), z(0)]), &m),
        Some(false)
    ); // 90°
    assert_eq!(
        arc_contains(&arc, &rat_dir(&[z(1), z(-1), z(0)]), &m),
        Some(false)
    ); // 315°
    // A direction carrying a radical — a root of the circle on a line — is read the same way.
    let root =
        |c: [(Rat, Rat); 2]| c.map(|(a, b)| nacre_exact::quad::QuadVal::new(a, b, z(3)).unwrap());
    let [x240, y240] = root([(q(-1, 2), z(0)), (z(0), q(-1, 2))]);
    let [x330, y330] = root([(z(0), q(1, 2)), (q(-1, 2), z(0))]);
    let zq = nacre_exact::quad::QuadVal::from_rat(z(0));
    assert_eq!(arc_contains(&arc, &[x240, y240, zq], &m), Some(true)); // 240°
    assert_eq!(arc_contains(&arc, &[x330, y330, zq], &m), Some(false)); // 330°
    let long = RimArc {
        from: arc.to,
        to: arc.from,
    }; // the other three quarters
    assert_eq!(
        arc_contains(&long, &rat_dir(&[z(1), z(1), z(0)]), &m),
        Some(true)
    );
    assert_eq!(
        arc_contains(&long, &rat_dir(&[z(-1), z(-1), z(0)]), &m),
        Some(false)
    );
    let half = RimArc {
        from: [z(1), z(0), z(0)],
        to: [z(-1), z(0), z(0)],
    }; // exactly a half turn, through +y
    assert_eq!(
        arc_contains(&half, &rat_dir(&[z(0), z(1), z(0)]), &m),
        Some(true)
    );
    assert_eq!(
        arc_contains(&half, &rat_dir(&[z(0), z(-1), z(0)]), &m),
        Some(false)
    );
    let fq = |span: Option<[Rat; 2]>| Footprint {
        span,
        theta: Theta::Arc(arc),
    };
    // `d = +y`: the peak direction `+y` (90°) is off the arc, so the far end is the larger of
    // the ends' `d·v` — `0` at `(−r, 0)` — with no radical; `−y` (270°) is the arc's own end,
    // so the near side keeps `ρ²`.
    let up = lateral_reach(&def, &fq(None), &[z(0), z(1), z(0)]).unwrap();
    assert_eq!(
        (up.lo, up.hi, up.rho2_lo, up.rho2_hi),
        (z(0), z(0), q(1, 25), z(0))
    );
    // `d = −y`: the mirror image.
    let down = lateral_reach(&def, &fq(None), &[z(0), z(-1), z(0)]).unwrap();
    assert_eq!(
        (down.lo, down.hi, down.rho2_lo, down.rho2_hi),
        (z(0), z(0), z(0), q(1, 25))
    );
    // `d = (1, 1, 0)`: `+d⊥` at 45° is off, `−d⊥` at 225° is on; both ends give `−1/5`.
    let diag = lateral_reach(&def, &fq(None), &[z(1), z(1), z(0)]).unwrap();
    assert_eq!(
        (diag.lo, diag.hi, diag.rho2_lo, diag.rho2_hi),
        (z(0), q(-1, 5), q(2, 25), z(0))
    );
    // And what the arc clears that the circle did not: a station at `+1/10` along `+y`.
    assert_eq!(reach_clears(&up, q(1, 10), q(1, 10)), Some(true));
    assert_eq!(reach_clears(&reach, q(1, 10), q(1, 10)), Some(false));
    // Along the axis the arc changes nothing.
    let along_arc = lateral_reach(&def, &fq(span), &[z(0), z(0), z(1)]).unwrap();
    assert_eq!(whole(&along_arc), (z(0), z(2), z(0)));
}

/// **The tangency is *stated*** — and every field of that
/// statement is asserted here, because the verdict downstream is only as good as this row.
///
/// The frozen census shape (`cylinder-wall-tangent`): a `2³` cube and a cylinder of radius
/// `0.5` whose axis stands at `x = 0.5`, so the wall `x = 0` is **exactly** `r` away. Hand
/// geometry throughout — nothing below is copied from an engine run.
///
/// Called directly rather than through the gate's ledger: a process-global read by tests
/// running in parallel attributes one fixture's rows to another, and the row builder is a pure
/// function of the setup, so the honest reading is to hand it that setup.
#[test]
fn the_tangent_wall_states_itself_exactly() {
    let mut m = Model::new();
    let cube = crate::fixtures::cuboid(&mut m, [0.0; 3], [2.0, 2.0, 2.0]);
    let cyl = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.5, 1.0, -1.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        4.0,
    )
    .solid;
    let (setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(&m, cube, cyl).expect("setup");
    let surf = cyl_surfs[0];
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let (o, r2) = (def.origin(), def.r2());
    // The tangent class, found the gate's own way: one clearance call per class.
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let c = (0..setup.geom.len())
        .find(|&k| {
            jd.world_rat(k).is_some_and(|w| {
                nacre_exact::point_plane_clearance_rat(&w, &o, r2) == nacre_exact::Orient::Zero
            })
        })
        .expect("the wall x = 0 is exactly r from the axis");
    let coeffs = jd.world_rat(c).expect("checked above");
    // No other plane of either operand holds the line `x = 0, y = 1` (a plane that did would
    // have to contain `+Z` and pass through `(0, 1)`).
    let holders = classes_holding_the_line(&jd, c, &coeffs, &def).expect("exact");
    assert!(holders.is_empty(), "{holders:?}");
    let rows = tangency_rows(
        &m,
        &setup.planes,
        &setup.plane_ix,
        setup.n_a,
        c,
        0,
        &coeffs,
        &def,
        surf,
        SolidSide::B, // the cylinder is the second operand; the wall face is the cube's
        &jd,
    );
    assert_eq!(rows.len(), 1, "one wall face, one lateral face: {rows:?}");
    let t = &rows[0];
    assert_eq!((t.wall, t.cyl), (c, 0));
    assert_eq!(
        (t.wall_solid, t.cyl_solid),
        (SolidSide::A, SolidSide::B),
        "the cube is the first operand"
    );
    let read = t.reading.as_ref().expect("the row is read");
    // The cube's material is `x > 0` and the whole cylinder stands there too.
    assert!(read.lens_in_wall_solid, "{t:?}");
    // A solid cylinder operand keeps its material inside.
    assert_eq!(t.cyl_orient, 1, "{t:?}");
    // The face `x = 0` runs from `y = 0` to `y = 2` across the tangent line `y = 1`, and its
    // corners sit at `t = 1` and `t = 3` inside the lateral's span `[0, 4]`, which runs on either
    // side of the line: both solids run through it over `[1, 3]`.
    let through: Vec<[Option<Rat>; 2]> =
        read.through.iter().map(|p| p.map(|q| q.as_rat())).collect();
    assert_eq!(
        through,
        vec![[Some(Rat::from_int(1)), Some(Rat::from_int(3))]],
        "{t:?}"
    );
    // The foot of the perpendicular is `(0, 1, −1)`; the span's middle carries it to `t = 2`.
    let w = read.witness.as_array();
    assert!(
        (w[0] - 0.0).abs() < 1e-12 && (w[1] - 1.0).abs() < 1e-12 && (w[2] - 1.0).abs() < 1e-12,
        "the witness sits on the contact, mid-span: {w:?}"
    );
}

/// **A line along the axis, against a lateral face's angular extent** — the half turn
/// `(1, 0) → (−1, 0)` counter-clockwise about `+z`, the upper half: a line through `(0, −1)` is shown
/// apart from the face, through `(0, 1)` is not, and neither is a line on an end (a line on a
/// face's end ruling is its edge). The witness is the line's own radial direction
/// ([`separating_dirs`]). The height the point is given at does not matter, a whole circle is never
/// apart from a line, and a second statement running its axis the other way reads the same arc as
/// the lower half — which holds `(0, −1)`.
#[test]
fn a_line_off_a_faces_angular_extent_is_apart_from_it_ends_excluded() {
    let z = Rat::from_int;
    let (o, m) = ([z(0); 3], [z(0), z(0), z(1)]);
    let def =
        nacre_topo::CylinderDef::new(o, m, [z(1), z(0), z(0)], nacre_exact::BigRat::from(z(1)))
            .unwrap();
    let upper = Footprint {
        span: None,
        theta: Theta::Arc(RimArc {
            from: [z(1), z(0), z(0)],
            to: [z(-1), z(0), z(0)],
        }),
    };
    let line = |p: [Rat; 3]| Footprint {
        span: None,
        theta: Theta::at_line(&p, &o, &m).unwrap(),
    };
    for (p, apart) in [
        ([z(0), z(1), z(5)], false),
        ([z(0), z(-1), z(0)], true),
        ([z(1), z(0), z(3)], false),
        ([z(-1), z(0), z(-2)], false),
    ] {
        assert_eq!(pieces_clear(&def, &line(p), &def, &upper), apart, "{p:?}");
        assert_eq!(
            pieces_clear(&def, &upper, &def, &line(p)),
            apart,
            "{p:?} swapped"
        );
        // The witness that decides a line, read three ways.
        assert_eq!(line_apart(&def, &line(p), &upper), Some(apart), "{p:?}");
    }
    let whole = Footprint {
        span: None,
        theta: Theta::Whole,
    };
    assert!(!pieces_clear(
        &def,
        &line([z(0), z(-1), z(0)]),
        &def,
        &whole
    ));
    let flipped = nacre_topo::CylinderDef::new(
        o,
        [z(0), z(0), z(-1)],
        [z(1), z(0), z(0)],
        nacre_exact::BigRat::from(z(1)),
    )
    .unwrap();
    assert!(!pieces_clear(
        &def,
        &line([z(0), z(-1), z(0)]),
        &flipped,
        &upper
    ));
}

/// **A line's angle is its point off the axis**, wherever the axis stands and whichever way it
/// runs: on the axis `(1, 2, 3) + s·(1, 1, 0)`, the point `(3, 2, 3)` is `(1, −1, 0)` off it (its
/// component along the axis, `(1, 1, 0)`, removed).
#[test]
fn a_lines_angle_is_its_point_off_the_axis() {
    let z = Rat::from_int;
    assert_eq!(
        Theta::at_line(
            &[z(3), z(2), z(3)],
            &[z(1), z(2), z(3)],
            &[z(1), z(1), z(0)]
        ),
        Some(Theta::At([z(1), z(-1), z(0)]))
    );
}

/// A half cylinder over the upper half of the unit disk, `z ∈ [0, 2]` (its lateral the half turn
/// `(1, 0) → (−1, 0)`), and a box `[lo, hi] × z ∈ [−0.5, 2.5]` — the box's height covers the
/// lateral's so the axis never clears a wall face first. The gate's tangency rows, or its refusal.
fn half_cylinder_and_box(lo: [f64; 2], hi: [f64; 2]) -> Result<Vec<Tangency>, BoolError> {
    let r = Rat::from_int;
    let mut m = Model::new();
    let half = crate::Ring2d::new(
        vec![[r(1), r(0)], [r(-1), r(0)]],
        vec![
            crate::Edge2d::Arc {
                center: [r(0), r(0)],
                r2: r(1),
                ccw: true,
            },
            crate::Edge2d::Line,
        ],
    )
    .expect("a half disk");
    let frame = crate::SketchFrame::world(&m, nacre_exact::Axis::Z);
    let Ok(crate::OpOutput::Extrude { solid: a, .. }) = crate::apply(
        &mut m,
        &crate::Operation::Extrude {
            frame,
            profile: crate::from_paths(vec![half]).unwrap().remove(0),
            dist: 2.0,
        },
    ) else {
        panic!("the half cylinder extrudes")
    };
    let b = crate::fixtures::cuboid(&mut m, [lo[0], lo[1], -0.5], [hi[0], hi[1], 2.5]);
    crate::arrangement::plane_index_setup(&m, a, b).map(|s| s.tangencies)
}

/// **A tangency row is about the lateral face, not its surface.** The half cylinder's surface
/// is tangent to three box walls; the gate writes a row only where the face is:
///
/// - `y = −1` touches the surface at `(0, −1)`, on the half the lateral does not have — no row;
/// - `y = 1` touches it at `(0, 1)`, on the lateral — one row, decided;
/// - `x = 1` touches it on the lateral's **end** ruling `(1, 0)` — a row (ends are the face's
///   edge), with nothing running through on both sides: the box's wall runs across the line, but
///   the half cylinder ends on it, its flat `y = 0` beside — an edge of the arrangement, the
///   structure's to judge.
#[test]
fn a_tangency_row_is_written_only_where_the_lateral_face_is() {
    let missing = half_cylinder_and_box([-0.5, -1.0], [0.5, 2.0]).expect("setup");
    assert!(missing.is_empty(), "{missing:?}");
    let present = half_cylinder_and_box([-0.5, 1.0], [0.5, 2.0]).expect("setup");
    assert_eq!(present.len(), 1, "{present:?}");
    assert!(!present[0].reading.through.is_empty(), "{present:?}");
    let end = half_cylinder_and_box([1.0, -1.0], [2.0, 2.0]).expect("setup");
    assert_eq!(end.len(), 1, "{end:?}");
    assert!(end[0].reading.through.is_empty(), "{end:?}");
}

/// **The same U face, read over two spans.** The block `[0, 6] × [−3, 3] × [0, 2]` with the slot
/// `y ∈ [−2, 2]`, `z ∈ [1, 3]` has a U-shaped face on `x = 0`, and a boss of radius `5/2` with
/// its axis at `(5/2, 0)` is tangent to that plane along `(0, 0)`. The face crosses the line
/// where `z ∈ [0, 1]` and nowhere above — the U's gap:
///
/// - the boss over `z ∈ [6/5, 9/5]` meets the line only in the gap, so the face touches nothing
///   on the lateral and writes no row (its corners lie on both sides of the line all the same);
/// - the boss over `z ∈ [1/2, 9/5]` meets it where the face is, `z ∈ (1/2, 1)`: one row, the line
///   through the face's interior, its witness the middle of that stretch, `z = 3/4`.
#[test]
fn a_u_face_writes_a_row_only_over_the_span_it_crosses() {
    let rows = |z0: f64, h: f64| {
        let mut m = Model::new();
        let block = crate::fixtures::cuboid(&mut m, [0.0, -3.0, 0.0], [6.0, 3.0, 2.0]);
        let slot = crate::fixtures::cuboid(&mut m, [-1.0, -2.0, 1.0], [7.0, 2.0, 3.0]);
        let a = crate::boolean(&mut m, crate::BoolKind::Cut, block, slot).expect("the slot")[0];
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([2.5, 0.0, z0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            nacre_math::Vector3::from_array([0.0, -1.0, 0.0]),
            2.5,
            h,
        )
        .solid;
        crate::arrangement::plane_index_setup(&m, a, b)
            .expect("the gate passes")
            .tangencies
    };
    let in_the_gap = rows(1.2, 0.6);
    assert!(in_the_gap.is_empty(), "{in_the_gap:?}");
    let across = rows(0.5, 1.3);
    assert_eq!(across.len(), 1, "{across:?}");
    assert!(!across[0].reading.through.is_empty(), "{across:?}");
    let w = across[0].reading.witness.as_array();
    assert!(
        w[0].abs() < 1e-12 && w[1].abs() < 1e-12 && (w[2] - 0.75).abs() < 1e-12,
        "the witness sits on the contact: {w:?}"
    );
}

/// **The same question from the lateral's side: a row only where the lateral face is on the
/// line.** [`crate::fixtures::windowed_boss`] has no lateral face on `x = 1, y = 0` over
/// `z ∈ [1.5, 2.5]`, and the slab's face `x = 1` is tangent to the cylinder along that line:
///
/// - the slab over `z ∈ [1.7, 2.3]` meets the line only in the window, so nothing touches the
///   lateral and the pair writes no row — though the surface is tangent there and the face's
///   angular extent and axial span both hold the line;
/// - the slab over `z ∈ [1, 3]` meets it where the face is, `z ∈ (1, 1.5) ∪ (2.5, 3)`: one row,
///   the line through the wall face's interior.
///
/// Both window widths (corners at `x = √0.91` and at `x = 4/5`) and both `seam` sides (the rims'
/// vertices at the window's angle, or opposite it; the window is an inner loop either way).
#[test]
fn a_lateral_window_writes_a_row_only_where_the_lateral_is() {
    let rows = |w: f64, seam: f64, z: [f64; 2]| {
        let mut m = Model::new();
        let boss = crate::fixtures::windowed_boss(&mut m, w, seam);
        let slab = crate::fixtures::cuboid(&mut m, [-2.0, -2.0, z[0]], [1.0, 2.0, z[1]]);
        crate::arrangement::plane_index_setup(&m, slab, boss)
            .expect("the gate passes")
            .tangencies
    };
    for w in [0.3, 0.6] {
        for seam in [1.0, -1.0] {
            let inside = rows(w, seam, [1.7, 2.3]);
            assert!(inside.is_empty(), "w {w} seam {seam}: {inside:?}");
            let longer = rows(w, seam, [1.0, 3.0]);
            assert_eq!(longer.len(), 1, "w {w} seam {seam}: {longer:?}");
            assert!(
                !longer[0].reading.through.is_empty(),
                "w {w} seam {seam}: {longer:?}"
            );
        }
    }
}

/// **A vertex's reach holds its proven interval, and only a proven one is finite.**
///
/// The sieve drops a candidate before the exact test whenever reaches miss a box, so a reach
/// narrower than the truth's interval drops a real contact in silence. A `Realized` reach must
/// contain `coord ± bound` on every axis — an ordinary half-ulp, a zero bound (an exact
/// coordinate), and a radius far below an `f64`'s — and a cache that proves nothing must
/// keep every candidate.
#[test]
fn a_reach_holds_the_proven_interval_and_nothing_else_is_finite() {
    let coord = nacre_math::Point3::from_array([0.1, -3.0, 1.0e10]);
    let bound = [
        nacre_exact::Mag::of(0.1).times(nacre_exact::Mag::pow2(-53)),
        nacre_exact::Mag::ZERO,
        nacre_exact::Mag::pow2(-2000),
    ];
    let (lo, hi) = cache_reach(&nacre_topo::PointCache::Realized { coord, bound });
    for k in 0..3 {
        let (c, r) = (coord.as_array()[k], bound[k].upper_f64());
        assert!(
            lo[k] < c && c < hi[k],
            "axis {k}: the coordinate itself is inside"
        );
        assert!(
            c - lo[k] >= r && hi[k] - c >= r,
            "axis {k}: the bound {r:e} is inside"
        );
        assert!(
            lo[k].is_finite() && hi[k].is_finite(),
            "axis {k}: a proven reach is finite"
        );
    }
    for cache in [
        nacre_topo::PointCache::Ceiling { coord },
        nacre_topo::PointCache::Unrealized { coord },
    ] {
        let (lo, hi) = cache_reach(&cache);
        assert!(
            lo.iter().all(|x| *x == f64::NEG_INFINITY) && hi.iter().all(|x| *x == f64::INFINITY),
            "{cache:?} proves no bound, so its reach is the whole line"
        );
    }
}

/// **A class whose faces cannot be listed still has a piece to ask about** — [`lateral_pieces`]
/// answers an unusable list (here none at all) with one piece of no span and the whole circle, so
/// «every piece is apart» is not an empty truth: a line at any angle is not shown apart from it.
#[test]
fn an_unusable_face_list_is_one_unbounded_piece() {
    let z = Rat::from_int;
    let mut m = Model::new();
    crate::fixtures::cylinder(
        &mut m,
        nacre_math::Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    let surf = (0u32..)
        .map_while(|i| m.face_handle_at(i))
        .map(|f| m.face(f).surface)
        .find(|&s| matches!(m.surface(s), Surface::Cylinder { .. }))
        .expect("a lateral");
    let pieces = lateral_pieces(&[], surf);
    assert_eq!(
        pieces,
        vec![Footprint {
            span: None,
            theta: Theta::Whole
        }]
    );
    let (o, axis) = ([z(0); 3], [z(0), z(0), z(1)]);
    let def =
        nacre_topo::CylinderDef::new(o, axis, [z(1), z(0), z(0)], nacre_exact::BigRat::from(z(1)))
            .unwrap();
    let line = Footprint {
        span: None,
        theta: Theta::at_line(&[z(0), z(-1), z(0)], &o, &axis).unwrap(),
    };
    assert!(!pieces.iter().all(|fp| pieces_clear(&def, &line, &def, fp)));
}

/// A prism over the triangle `(y, z)` = `tri` on the world `YZ` plane, swept along `+x`; its face on
/// `x = 0`, and that face's plane as the gate names it.
fn yz_prism_face(tri: [(f64, f64); 3]) -> (Model, Handle<Face>, [Rat; 4]) {
    let mut m = Model::new();
    let pts = tri.map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = crate::Ring2d::polygon_decimal(pts.to_vec()).expect("a triangle");
    let profile = crate::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = crate::SketchFrame::world(&m, nacre_exact::Axis::X);
    let Ok(crate::OpOutput::Extrude { solid, .. }) = crate::apply(
        &mut m,
        &crate::Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    ) else {
        panic!("the prism extrudes")
    };
    let zero = Rat::from_int(0);
    let (face, name) = m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find_map(|&f| {
            let name = *m.world_plane_name(m.face(f).surface)?.narrow()?;
            (name[1] == zero && name[2] == zero && name[3] == zero).then_some((f, name))
        })
        .expect("the face on x = 0");
    (m, face, name)
}

/// **A face meeting the line at a single corner does not touch it, and the contact is where the
/// face meets the line.** The line `x = 0, y = 1` along `z`, read over `0 < z < 2` on the plane `x
/// = 0`: a triangle with one corner on it at `z = 1.75` and the rest at `y > 1` meets it at that
/// point alone — no positive length, so it does not touch: its contact is empty; a triangle with an
/// edge on it over `1 ≤ z ≤ 1.5` touches there, along its boundary, never through its interior, and
/// its contact is that edge; a triangle the line crosses runs through it over `0.5 ≤ z ≤ 1.5`, its
/// contact. The corner's other two steps both cross the line at `1.75` under the reading that
/// counts a piece on it as above, so the parity stays even.
#[test]
fn a_face_meeting_the_line_at_a_corner_does_not_touch_it() {
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = Rat::from_int;
    let read = |tri: [(f64, f64); 3]| {
        let (m, face, coeffs) = yz_prism_face(tri);
        // The frame puts `u` on `y` and `v` on `z`: the triangle's first corner is `(0, y, z)`.
        let (y, zz) = tri[0];
        assert!(
            m.face(face).outer.half_edges.iter().any(|he| {
                m.vertex_point(m.edge(he.edge).vertices[0])
                    .distance(Point3::from_array([0.0, y, zz]))
                    < 1e-12
            }),
            "the corner (0, {y}, {zz})"
        );
        let c = line_runs_through_face(
            &m,
            m.face(face),
            &coeffs,
            &[z(0), z(1), z(0)],
            &[z(0), z(0), z(1)],
            [z(0), z(2)],
        )
        .expect("the face reads");
        let contact: Vec<[Rat; 2]> = c
            .contact
            .iter()
            .map(|p| p.map(|x| x.as_rat().expect("a rational end")))
            .collect();
        let at = row_witness(&[([z(0), z(2)], c.clone())]).and_then(|(_, at)| at);
        (!c.contact.is_empty(), c.runs_through(), at, contact)
    };
    assert_eq!(
        read([(1.0, 1.75), (2.5, 1.2), (2.5, 2.5)]),
        (false, false, None, vec![])
    );
    assert_eq!(
        read([(1.0, 1.0), (1.0, 1.5), (2.5, 1.25)]),
        (true, false, Some(r(5, 4)), vec![[z(1), r(3, 2)]])
    );
    assert_eq!(
        read([(1.0, 1.5), (0.5, 0.5), (1.5, 0.5)]),
        (true, true, Some(z(1)), vec![[r(1, 2), r(3, 2)]])
    );
}

/// **A row's witness is the first piece the line runs through, over its stretches in order, else
/// the first piece it meets** ([`row_witness`]). A through piece in a later stretch outranks a
/// boundary piece in an earlier one; a stretch the face clears has no piece and is passed over; and
/// a through piece with irrational ends keeps its stretch with no middle rather than falling back
/// to a rational piece before it.
#[test]
fn the_rows_witness_is_the_first_piece_the_line_runs_through() {
    use nacre_exact::QuadVal;
    let r = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = Rat::from_int;
    let q = |v: i128| QuadVal::from_rat(z(v));
    let (along, through) = ([true, false], [true, true]);
    let contact = |pieces: &[([QuadVal; 2], [bool; 2])]| LineContact {
        contact: pieces.iter().map(|p| p.0).collect(),
        sides: pieces.iter().map(|p| p.1).collect(),
    };
    let (first, second) = ([z(0), z(4)], [z(4), z(8)]);
    assert_eq!(
        row_witness(&[
            (first, contact(&[([q(0), q(1)], along)])),
            (
                second,
                contact(&[([q(4), q(5)], along), ([q(6), q(7)], through)])
            ),
        ]),
        Some((second, Some(r(13, 2)))),
        "a through piece in a later stretch"
    );
    assert_eq!(
        row_witness(&[
            (first, contact(&[])),
            (
                second,
                contact(&[([q(4), q(5)], along), ([q(6), q(7)], along)])
            ),
        ]),
        Some((second, Some(r(9, 2)))),
        "nothing through, an empty stretch passed over"
    );
    let irrational = QuadVal::new(z(2), r(1, 4), z(2)).expect("2 + √2/4");
    assert_eq!(
        row_witness(&[(
            first,
            contact(&[([q(0), q(1)], along), ([q(2), irrational], through)])
        )]),
        Some((first, None)),
        "an irrational through piece"
    );
    assert_eq!(row_witness(&[(first, contact(&[]))]), None, "no piece");
}

/// **A plane against a face's boundary sections: clear of them, reaching them, or touching one.**
/// The cylinder about the `z` axis, radius 2; each plane is `z = k·y + c`, its section's height
/// running over `c ± 2k`. Against a band bounded by the circles `z = 0` and `z = 10` (or `dir =
/// (0, 0, 2)`, which leaves the planes as they are): between them the section reaches the band,
/// beyond both it is clear, across a rim it reaches it, and ending exactly at a rim — one point,
/// `c − 2k = 0` — it touches. Against the one circle `z = 0` a section beside it is clear, across it
/// reaches, and ending on it from above or below touches. A boundary section that is the plane
/// itself reaches it (a face's own slanted cap), and two circles never touch.
#[test]
fn a_plane_stands_clear_of_reaching_or_touching_a_faces_boundary() {
    use crate::planes::BoundaryStand::{Clear, Sections, Touch};
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let def = |dz: i128| {
        nacre_topo::CylinderDef::new(
            [z(0); 3],
            [z(0), z(0), z(dz)],
            [z(1), z(0), z(0)],
            nacre_exact::BigRat::from(z(4)),
        )
        .unwrap()
    };
    // `z = k·y + c` is `−k·y + z − c = 0`.
    let plane = |k: Rat, c: Rat| {
        [
            z(0),
            z(0).checked_sub(k).unwrap(),
            z(1),
            z(0).checked_sub(c).unwrap(),
        ]
    };
    let circle = |h: i128| [z(0), z(0), z(1), z(-h)];
    for dz in [1, 2] {
        let band = [circle(0), circle(10)];
        let at = |k: Rat, c: Rat| boundary_stand(&def(dz), &band, &plane(k, c));
        assert_eq!(at(q(1, 4), q(5, 2)), Some(Sections), "z ∈ [2, 3]");
        assert_eq!(at(q(1, 4), z(-5)), Some(Clear), "below the band");
        assert_eq!(at(q(1, 4), z(13)), Some(Clear), "above it");
        assert_eq!(
            at(z(3), z(5)),
            Some(Sections),
            "z ∈ [−1, 11]: across both rims"
        );
        assert_eq!(
            at(q(1, 4), q(1, 2)),
            Some(Touch),
            "z ∈ [0, 1]: a rim's point"
        );
        assert_eq!(at(q(1, 4), q(19, 2)), Some(Touch), "z ∈ [9, 10]");
        assert_eq!(at(z(0), z(5)), Some(Sections), "the circle z = 5, between");
        assert_eq!(at(z(0), z(10)), Some(Sections), "the rim's own plane");
    }
    let rim = [circle(0)];
    let at = |k: Rat, c: Rat| boundary_stand(&def(1), &rim, &plane(k, c));
    assert_eq!(at(q(1, 4), z(1)), Some(Clear), "z ∈ [1/2, 3/2]");
    assert_eq!(at(q(1, 4), z(-1)), Some(Clear), "z ∈ [−3/2, −1/2]");
    assert_eq!(at(z(3), z(5)), Some(Sections), "z ∈ [−1, 11]");
    assert_eq!(at(q(1, 4), q(2, 5)), Some(Sections), "z ∈ [−1/10, 9/10]");
    assert_eq!(at(q(1, 4), q(1, 2)), Some(Touch), "z ∈ [0, 1]");
    assert_eq!(at(q(1, 4), q(-1, 2)), Some(Touch), "z ∈ [−1, 0]");
    // The circle's own plane, with no other boundary to stand between: reached by being it.
    assert_eq!(at(z(0), z(0)), Some(Sections), "the circle z = 0 itself");
}

/// **A class carries closed traces no road nests** when an oblique section stands beside a circle
/// or a section of an axis that is not parallel — and not when the axes agree, nor beside rulings,
/// nor with no section at all.
#[test]
fn a_section_beside_another_axis_carries_traces_no_road_nests() {
    let z = |v: i128| Rat::from_int(v);
    let (up, up_long, along) = ([z(0), z(0), z(1)], [z(0), z(0), z(3)], [z(1), z(0), z(0)]);
    let (sec, across, cuts) = (ClassMeet::Sections, ClassMeet::Across, ClassMeet::Cuts);
    assert!(
        carries_unnested(&[(up, &sec), (along, &across)]),
        "a section beside a circle"
    );
    assert!(
        carries_unnested(&[(up, &sec), (along, &sec)]),
        "two sections of crossing axes"
    );
    assert!(
        !carries_unnested(&[(up, &sec), (up_long, &sec)]),
        "two sections of parallel axes"
    );
    assert!(
        !carries_unnested(&[(up, &sec), (along, &cuts)]),
        "a section beside rulings"
    );
    assert!(
        !carries_unnested(&[(up, &across), (along, &across)]),
        "no section: the class is not asked"
    );
}

/// A slab swept along `+x` from `x = 0` to `x = 2` over the rings, sketched on `frame` (a frame on
/// `x = 0` facing `+x`), its face on `x = x0` (`0` or `2`), and that face's plane as the gate names
/// it.
fn yz_slab_face(
    rings: Vec<crate::Ring2d>,
    frame: impl FnOnce(&mut Model) -> crate::SketchFrame,
    x0: i128,
) -> (Model, Handle<Face>, [Rat; 4]) {
    let mut m = Model::new();
    let profile = crate::from_paths(rings).expect("one profile").remove(0);
    let frame = frame(&mut m);
    let out = crate::apply(
        &mut m,
        &crate::Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    );
    let Ok(crate::OpOutput::Extrude { solid, .. }) = out else {
        panic!("the slab extrudes: {out:?}")
    };
    let zero = Rat::from_int(0);
    let (face, name) = m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find_map(|&f| {
            let name = *m.world_plane_name(m.face(f).surface)?.narrow()?;
            (name[1] == zero && name[2] == zero && name[3] == Rat::from_int(-x0))
                .then_some((f, name))
        })
        .expect("the face on x = x0");
    (m, face, name)
}

/// The rectangle `[y0, y1] × [z0, z1]` with its corners rounded — radii lower-left, lower-right,
/// upper-right, upper-left — counter-clockwise.
fn rounded_ring(lo: [Rat; 2], hi: [Rat; 2], radii: [Rat; 4]) -> crate::Ring2d {
    let ([y0, z0], [y1, z1]) = (lo, hi);
    let [ll, lr, ur, ul] = radii;
    let add = |a: Rat, b: Rat| a.checked_add(b).unwrap();
    let sub = |a: Rat, b: Rat| a.checked_sub(b).unwrap();
    let arc = |c: [Rat; 2], r: Rat| crate::Edge2d::Arc {
        center: c,
        r2: r.checked_mul(r).unwrap(),
        ccw: true,
    };
    crate::Ring2d::new(
        vec![
            [add(y0, ll), z0],
            [sub(y1, lr), z0],
            [y1, add(z0, lr)],
            [y1, sub(z1, ur)],
            [sub(y1, ur), z1],
            [add(y0, ul), z1],
            [y0, sub(z1, ul)],
            [y0, add(z0, ll)],
        ],
        vec![
            crate::Edge2d::Line,
            arc([sub(y1, lr), add(z0, lr)], lr),
            crate::Edge2d::Line,
            arc([sub(y1, ur), sub(z1, ur)], ur),
            crate::Edge2d::Line,
            arc([add(y0, ul), sub(z1, ul)], ul),
            crate::Edge2d::Line,
            arc([add(y0, ll), add(z0, ll)], ll),
        ],
    )
    .expect("a rounded rectangle")
}

/// **A tangent line meets a filleted face where each piece of its boundary changes sides of it.**
/// The line `x = x0, y = c` along `z`, read over a span, on a slab face whose outline carries arcs
/// (`line_runs_through_face`). On the square `[0, 4]²` with every corner rounded `r = 1`, over
/// `0 < z < 4`, on both of the slab's faces (`x = 0` and `x = 2` walk the loop opposite ways):
///
/// - `c = 2` meets only the straight edges: the whole span, through the interior;
/// - `c = 1/2` crosses the two left arcs inside them, at `z = 1 − √3/2` and `z = 3 + √3/2` (the
///   circles' other roots, `1 + √3/2` and `3 − √3/2`, lie off the arcs);
/// - `c = 1` runs through the arcs' ends `(1, 0)` and `(1, 4)`, where each arc meets a straight
///   edge: each end changes sides once, the whole span through the interior;
/// - `c = 0` lies along the left edge, and the arcs leave its ends `(0, 1)`, `(0, 3)` along the
///   line — each end is its circle's extreme point across the line, and the arc bends to the
///   face's side: `[1, 3]`, on that side only (`w = n × m = −y`, so the `−w` reading).
///
/// An arc bending away from the face, stated clockwise on the outer loop: the same square less the
/// quarter disk `r = 1` about its corner `(0, 0)`, on both faces. `c = 1/2` crosses that arc inside
/// it, at `√3/2`: `[√3/2, 4]`. `c = 1` meets the arc only at its end `(1, 0)`, its circle's extreme
/// point across the line, where the arc bends to the centre's side — off the face — and the bottom
/// edge leaves to the face's side: one change, `[0, 4]` through the interior. `c = 0` runs along the
/// left edge, whose end `(0, 1)` the arc leaves across the line: `[1, 4]`, on the `−w` side only.
///
/// Two radicals on one line: lower-left `r = 1`, upper-left `r = 1/2`, `c = 1/4` crosses them at
/// `1 − √7/4` and `7/2 + √3/4`. An end off the line with the same `v` as a root is another point:
/// on a 3-4-5 frame, a quarter arc about `(3/5, 4)`, `r = 1`, from `(−1/5, 17/5)` to `(6/5, 16/5)`
/// meets `y = 0` at `(0, 16/5)` — the end's `z` — and the face runs from there to where the line
/// leaves its far edge, `z = 73/10`.
/// And a hole whose arcs end on the line: the square `[−2, 6]²` less the rounded square, `c = 1`,
/// is the face on `[−2, 0]` and `[4, 6]`.
#[test]
fn a_tangent_line_meets_a_filleted_face_where_each_piece_crosses() {
    use nacre_exact::QuadVal;
    // A stretch of the line the face meets, and its sides.
    type Piece = ([QuadVal; 2], [bool; 2]);
    let z = Rat::from_int;
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let qv = |a: Rat, b: Rat, c: Rat| QuadVal::new(a, b, c).unwrap();
    let r = |a: Rat| QuadVal::from_rat(a);
    // The world frame on `x = 0`: `u` on `y`, `v` on `z`.
    let world = |m: &mut Model| crate::SketchFrame::world(m, nacre_exact::Axis::X);
    let read_on = |rings: Vec<crate::Ring2d>,
                   frame: &dyn Fn(&mut Model) -> crate::SketchFrame,
                   x0: i128,
                   c: Rat,
                   span: [Rat; 2]| {
        let (m, face, coeffs) = yz_slab_face(rings, frame, x0);
        line_runs_through_face(
            &m,
            m.face(face),
            &coeffs,
            &[z(x0), c, z(0)],
            &[z(0), z(0), z(1)],
            span,
        )
        .ok()
        .map(|lc| lc.contact.into_iter().zip(lc.sides).collect::<Vec<Piece>>())
    };
    let read = |rings, x0, c, span| read_on(rings, &world, x0, c, span);
    let same = |got: &Option<Vec<Piece>>, want: &[Piece]| {
        let Some(got) = got else {
            return false;
        };
        got.len() == want.len()
            && got.iter().zip(want).all(|((p, s), (wp, ws))| {
                s == ws
                    && p.iter()
                        .zip(wp)
                        .all(|(a, b)| a.cmp_value(b) == Some(core::cmp::Ordering::Equal))
            })
    };
    let square = || rounded_ring([z(0), z(0)], [z(4), z(4)], [z(1); 4]);
    let both = [true, true];
    for x0 in [0, 2] {
        let cases: [(Rat, Vec<Piece>); 4] = [
            (z(2), vec![([r(z(0)), r(z(4))], both)]),
            (
                q(1, 2),
                vec![([qv(z(1), z(-1), q(3, 4)), qv(z(3), z(1), q(3, 4))], both)],
            ),
            (z(1), vec![([r(z(0)), r(z(4))], both)]),
            (z(0), vec![([r(z(1)), r(z(3))], [false, true])]),
        ];
        for (c, want) in cases {
            let got = read(vec![square()], x0, c, [z(0), z(4)]);
            assert!(same(&got, &want), "x = {x0}, y = {c:?}: {got:?}");
        }
        let notched = || {
            crate::Ring2d::new(
                vec![
                    [z(1), z(0)],
                    [z(4), z(0)],
                    [z(4), z(4)],
                    [z(0), z(4)],
                    [z(0), z(1)],
                ],
                vec![
                    crate::Edge2d::Line,
                    crate::Edge2d::Line,
                    crate::Edge2d::Line,
                    crate::Edge2d::Line,
                    crate::Edge2d::Arc {
                        center: [z(0), z(0)],
                        r2: z(1),
                        ccw: false,
                    },
                ],
            )
            .expect("a square with a notched corner")
        };
        let cases: [(Rat, Vec<Piece>); 3] = [
            (q(1, 2), vec![([qv(z(0), z(1), q(3, 4)), r(z(4))], both)]),
            (z(1), vec![([r(z(0)), r(z(4))], both)]),
            (z(0), vec![([r(z(1)), r(z(4))], [false, true])]),
        ];
        for (c, want) in cases {
            let got = read(vec![notched()], x0, c, [z(0), z(4)]);
            assert!(same(&got, &want), "notched, x = {x0}, y = {c:?}: {got:?}");
        }
    }
    let two_radii = rounded_ring([z(0), z(0)], [z(4), z(4)], [z(1), z(1), z(1), q(1, 2)]);
    let got = read(vec![two_radii], 0, q(1, 4), [z(0), z(4)]);
    let want = [(
        [qv(z(1), z(-1), q(7, 16)), qv(q(7, 2), z(1), q(3, 16))],
        both,
    )];
    assert!(same(&got, &want), "two radicals: {got:?}");
    // The quarter arc about the frame's origin from `−v` to `+u`, on the 3-4-5 frame at
    // `(y, z) = (3/5, 4)` with `u = (3/5, −4/5)`, `v = (4/5, 3/5)`: its end `+u` is `(6/5, 16/5)`, and
    // `y = 0` meets the circle at `(0, 16/5)` — on the arc, where the radial is `(7/25, −24/25)` in
    // the frame — and at `(0, 24/5)`, the radial `−u`, off it. The face is the frame's rectangle
    // `[−3, 1] × [−1, 2]` with that corner rounded; `y = 0` leaves it on the left edge, `u = −3`, at
    // `z = 73/10`.
    let tilted = |m: &mut Model| match crate::apply(
        m,
        &crate::Operation::DatumPlane {
            def: crate::DatumDef::Stated(crate::SketchPlane::from_axes(
                Point3::from_array([0.0, 0.6, 4.0]),
                nacre_math::Vector3::from_array([0.0, 0.6, -0.8]),
                nacre_math::Vector3::from_array([0.0, 0.8, 0.6]),
            )),
        },
    ) {
        Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("the 3-4-5 frame: {other:?}"),
    };
    let bitten = crate::Ring2d::new(
        vec![
            [z(0), z(-1)],
            [z(1), z(0)],
            [z(1), z(2)],
            [z(-3), z(2)],
            [z(-3), z(-1)],
        ],
        vec![
            crate::Edge2d::Arc {
                center: [z(0), z(0)],
                r2: z(1),
                ccw: true,
            },
            crate::Edge2d::Line,
            crate::Edge2d::Line,
            crate::Edge2d::Line,
            crate::Edge2d::Line,
        ],
    )
    .expect("a rectangle with a rounded corner");
    let got = read_on(vec![bitten], &tilted, 0, z(0), [z(0), z(8)]);
    assert!(
        same(&got, &[([r(q(16, 5)), r(q(73, 10))], both)]),
        "a root beside an end off the line: {got:?}"
    );
    let frame = crate::Ring2d::polygon_decimal(
        [(-2.0, -2.0), (6.0, -2.0), (6.0, 6.0), (-2.0, 6.0)]
            .map(|(u, v)| nacre_math::Point2::from_array([u, v]))
            .to_vec(),
    )
    .expect("a square");
    let got = read(vec![frame, square()], 0, z(1), [z(-2), z(6)]);
    assert!(
        same(
            &got,
            &[([r(z(-2)), r(z(0))], both), ([r(z(4)), r(z(6))], both)]
        ),
        "a hole whose arcs end on the line: {got:?}"
    );
}

/// **A wall face the footprint reader cannot spell rides its tangency row with its own refusal,
/// and the gate reaches that row only when the class's first face is another.** The D rod and the
/// bitten slab of `d_rod_on_a_bitten_slab`, against the class `x = 1` and the rod's cylinder: the
/// gate's two readers of the class raise the very same value — [`face_clears_footprint`]'s
/// `CurvedOperandBoundary`, with no place — whichever reaches the slab's face, so the boolean's
/// answer cannot say which did, and each is asked here directly (the row builder is a pure
/// function of the setup, as [`the_tangent_wall_states_itself_exactly`] says).
///
/// - [`wall_faces_clear`] stops at the class's first face that does not clear: with the rod first
///   that is the rod's own flat (`Ok(false)`, the class is read row by row), with the slab first it
///   is the slab's face, whose refusal it raises.
/// - [`tangency_rows`] writes one row for the slab's face against the rod's lateral in either
///   order: it carries that refusal, and its stretches through the line are not read (`None`) —
///   an unread face is unknown, not innocent.
#[test]
fn an_unread_wall_face_rides_its_tangency_row() {
    let z = Rat::from_int;
    for rod_first in [true, false] {
        let mut m = Model::new();
        let (rod, slab) = crate::tests::d_rod_on_a_bitten_slab(&mut m);
        let (a, b) = if rod_first { (rod, slab) } else { (slab, rod) };
        let (setup, cyl_surfs) =
            crate::arrangement::plane_index_setup_inner(&m, a, b).expect("setup");
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        // The class `x = 1` by its name: `y = −1` is tangent to the rod's cylinder too.
        let c = (0..setup.geom.len())
            .find(|&k| {
                jd.world_rat(k).is_some_and(|w| {
                    w[0] != z(0)
                        && w[1] == z(0)
                        && w[2] == z(0)
                        && w[0].checked_add(w[3]) == Some(z(0))
                })
            })
            .expect("the class x = 1");
        let coeffs = jd.world_rat(c).expect("checked above");
        // The rod's cylinder runs along `z`; the bore that bit the slab runs along `x`.
        let (ci, surf, def) = cyl_surfs
            .iter()
            .enumerate()
            .find_map(|(i, &s)| {
                let d = m.world_cylinder_def(s)?;
                (d.dir()[0] == z(0) && d.dir()[1] == z(0)).then_some((i, s, d))
            })
            .expect("the rod's cylinder");
        let (rod_side, slab_side) = if rod_first {
            (SolidSide::A, SolidSide::B)
        } else {
            (SolidSide::B, SolidSide::A)
        };
        let order = if rod_first { "rod first" } else { "slab first" };
        let spans = lateral_spans(&setup.planes, surf);
        let class_wide = wall_faces_clear(
            &m,
            &setup.planes,
            &setup.plane_ix,
            c,
            &coeffs,
            &def.origin(),
            &def.dir(),
            def.r2(),
            &spans,
        );
        if rod_first {
            assert!(matches!(class_wide, Ok(false)), "{order}: {class_wide:?}");
        } else {
            assert!(
                matches!(
                    class_wide,
                    Err(BoolError::Rejected {
                        reason: RejectReason::CurvedOperandBoundary,
                        at: None,
                    })
                ),
                "{order}: {class_wide:?}"
            );
        }
        let rows = tangency_rows(
            &m,
            &setup.planes,
            &setup.plane_ix,
            setup.n_a,
            c,
            ci,
            &coeffs,
            &def,
            surf,
            rod_side,
            &jd,
        );
        assert_eq!(rows.len(), 1, "{order}: {rows:?}");
        let t = &rows[0];
        assert_eq!((t.wall, t.cyl), (c, ci), "{order}");
        assert_eq!(
            (t.wall_solid, t.cyl_solid),
            (slab_side, rod_side),
            "{order}"
        );
        assert!(
            matches!(
                t.reading,
                Err(BoolError::Rejected {
                    reason: RejectReason::CurvedOperandBoundary,
                    at: None,
                })
            ),
            "{order}: {t:?}"
        );
    }
}

/// **A wall face's material side is read from its facing, whatever the face states of itself.**
/// The rounded plate's corner fillet (axis `x = 40, y = −20`, `r = 5`) is tangent to the class
/// `x = 45`, whose root is the plate's own side (outward `+x`); a block's face on `x = 45` is the
/// wall. Two placements give both signs and both answers: beside the plate (the block's face looks
/// `−x`, against the root, and the lens is outside the block) and over it (`+x`, the lens inside).
/// The wall face is never the root, so a reading that dropped the facing or the root's frame sign
/// would answer the first placement wrongly.
///
/// [`tangency_rows`] is asked with three tables that differ only in the wall face's row: as built
/// (one surface with the root — [`face_facing`]'s first tier), on another surface handle (its world
/// name is the root's — the second tier), and on another handle with no world name (the judge's
/// witnesses — the third). The rows are the same in all three.
#[test]
fn a_wall_face_reads_its_material_side_from_its_facing() {
    for (lo, hi, facing, lens_in) in [
        ([45.0, -30.0, 0.0], [50.0, -10.0, 12.0], -1i8, false),
        ([35.0, -30.0, 0.0], [45.0, -10.0, 12.0], 1, true),
    ] {
        let (m, setup, c, ci, surf, def) = fillet_against_block(lo, hi);
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let coeffs = jd.world_rat(c).expect("narrow");
        let in_class: Vec<usize> = (0..setup.planes.len())
            .filter(|&i| setup.plane_ix[i] == ClassIx::Plane(c))
            .collect();
        let wall = *in_class
            .iter()
            .find(|&&i| i >= setup.n_a)
            .expect("the block's face on x = 45");
        assert_ne!(wall, in_class[0], "the wall face is not the class root");
        let root_surf = jd.planes[c].surf;
        let other = setup
            .planes
            .iter()
            .find_map(|row| match row {
                FaceRow::Plane(p) if p.surf != root_surf => Some(p.surf),
                _ => None,
            })
            .expect("another plane");
        let mut named = setup.planes.clone();
        if let FaceRow::Plane(fi) = &mut named[wall] {
            fi.surf = other;
        }
        let mut judged = named.clone();
        if let FaceRow::Plane(fi) = &mut judged[wall] {
            fi.world = None;
        }
        let rows_of = |faces: &[FaceRow]| {
            assert_eq!(face_facing(&jd, faces, wall, c), Ok(facing), "{lo:?}");
            tangency_rows(
                &m,
                faces,
                &setup.plane_ix,
                setup.n_a,
                c,
                ci,
                &coeffs,
                &def,
                surf,
                SolidSide::A,
                &jd,
            )
        };
        let built = rows_of(&setup.planes);
        assert!(!built.is_empty(), "{lo:?}: the block's face writes a row");
        for t in &built {
            let read = t.reading.as_ref().expect("the row is read");
            assert_eq!(read.lens_in_wall_solid, lens_in, "{lo:?}: {t:?}");
        }
        assert_eq!(
            setup.planes[wall].plane().surf,
            root_surf,
            "{lo:?}: one surface with the root"
        );
        assert!(named[wall].plane().surf != root_surf && named[wall].plane().world.is_some());
        assert!(judged[wall].plane().surf != root_surf && judged[wall].plane().world.is_none());
        for (tier, faces) in [("named", &named), ("judged", &judged)] {
            assert_eq!(
                format!("{:?}", rows_of(faces)),
                format!("{built:?}"),
                "{lo:?}: the {tier} table"
            );
        }
    }
}

/// The rounded plate of [`a_wall_face_reads_its_material_side_from_its_facing`] against a block
/// `[lo, hi]`, and the gate's setup (the plate first): the class `x = 45` and the corner fillet's
/// cylinder, tangent to it along `x = 45, y = −20`.
fn fillet_against_block(
    lo: [f64; 3],
    hi: [f64; 3],
) -> (
    Model,
    crate::arrangement::PlaneSetup,
    usize,
    usize,
    Handle<Surface>,
    nacre_topo::CylinderDef,
) {
    let z = Rat::from_int;
    let mut m = Model::new();
    let plate = crate::tests::rounded_plate(&mut m, 1, 0);
    let block = crate::fixtures::cuboid(&mut m, lo, hi);
    let (setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(&m, plate, block).expect("setup");
    let (c, coeffs) = {
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        let c = (0..setup.geom.len())
            .find(|&k| {
                jd.world_rat(k).is_some_and(|w| {
                    w[0] != z(0)
                        && w[1] == z(0)
                        && w[2] == z(0)
                        && w[0].checked_mul(z(45)).and_then(|a| a.checked_add(w[3])) == Some(z(0))
                })
            })
            .expect("the class x = 45");
        (c, jd.world_rat(c).expect("checked above"))
    };
    let (ci, surf, def) = cyl_surfs
        .iter()
        .enumerate()
        .find_map(|(i, &s)| {
            let d = m.world_cylinder_def(s)?;
            (nacre_exact::point_plane_clearance_rat(&coeffs, &d.origin(), d.r2())
                == nacre_exact::Orient::Zero)
                .then_some((i, s, d))
        })
        .expect("the fillet is tangent to x = 45");
    (m, setup, c, ci, surf, def)
}

/// **A tangency row refuses by the first cause that stopped it, named and placed as
/// [`tangency_rows`]' reading orders them.** The block beside the rounded plate of
/// [`a_wall_face_reads_its_material_side_from_its_facing`] writes one row, read; each edit below
/// reaches a cause no row of the corpus carries (the row builder is a pure function of the setup,
/// as [`the_tangent_wall_states_itself_exactly`] says), and they are stacked from the last cause to
/// the first, so every step moves the answer to the earlier cause's name and place:
///
/// 1. the fillet's lateral row without its face — its sides of the line cannot be read, so the
///    stretches through the line cannot either: `CylinderGateUndecided` at the witness;
/// 2. and the lateral's stretch `[1/p, 12 − 1/q]` (`p`, `q` coprime near `2⁶⁴`), whose middle
///    leaves `Rat`: no witness, `WitnessNotRational` with no place;
/// 3. and the wall face on another handle with no world name and a degenerate witness triangle —
///    its facing is the judge's `Degenerate`: `DegenerateWitness`;
/// 4. and the class with no world statement — no lens' side: `CylinderGateUndecided`;
/// 5. and the class's coefficients scaled by `2⁶⁴` (the same plane), whose `n·n` leaves `Rat`: no
///    foot, `WitnessNotRational`.
///
/// Alone: the axis direction scaled by `2⁴⁰` stops the lateral's reading for its width, which the
/// stretches through the line carry to the witness; by `2⁶⁴ + 1` it stops the reading of whether
/// the line lies on the lateral at all. And the D rod's bitten slab face, which the footprint
/// reader cannot spell, carries its own refusal under the last two edits of the stack.
#[test]
fn a_tangency_row_refuses_by_its_first_cause() {
    let z = Rat::from_int;
    let (m, setup, c, ci, surf, def) =
        fillet_against_block([45.0, -30.0, 0.0], [50.0, -10.0, 12.0]);
    let coeffs = Judge::new(&setup.geom, setup.standard, &setup.notes)
        .world_rat(c)
        .expect("narrow");
    let lateral = (0..setup.planes.len())
        .find(|&i| matches!(&setup.planes[i], FaceRow::Cylinder(cf) if cf.surf == surf))
        .expect("the fillet's lateral");
    let wall = (setup.n_a..setup.planes.len())
        .find(|&i| setup.plane_ix[i] == ClassIx::Plane(c))
        .expect("the block's face on x = 45");
    let other = setup
        .planes
        .iter()
        .find_map(|row| match row {
            FaceRow::Plane(p) if p.surf != setup.planes[wall].plane().surf => Some(p.surf),
            _ => None,
        })
        .expect("another plane");
    let big = (1i128 << 64) + 1;
    let scaled = |by: i128| {
        nacre_topo::CylinderDef::new(
            def.origin(),
            def.dir().map(|v| v.checked_mul(z(by)).expect("fits")),
            def.ref_dir(),
            def.r2().clone(),
        )
        .expect("a cylinder")
    };
    // The rows of an edited setup: the first `steps` edits of the stack, and the axis `dir`.
    let rows = |steps: usize, dir: &nacre_topo::CylinderDef| {
        let mut faces = setup.planes.clone();
        let mut geom = setup.geom.clone();
        let mut coeffs = coeffs;
        if steps >= 1
            && let FaceRow::Cylinder(cf) = &mut faces[lateral]
        {
            cf.face = None;
        }
        if steps >= 2
            && let FaceRow::Cylinder(cf) = &mut faces[lateral]
        {
            let q = big + 2;
            cf.stretch = Some([
                Rat::new(1, big).expect("a ratio"),
                Rat::new(12 * q - 1, q).expect("a ratio"),
            ]);
        }
        if steps >= 3
            && let FaceRow::Plane(fi) = &mut faces[wall]
        {
            fi.surf = other;
            fi.world = None;
            let a = fi.tri_pt3[0].clone();
            fi.tri_pt3 = [a.clone(), a.clone(), a];
        }
        if steps >= 4 {
            geom[c].world = None;
        }
        if steps >= 5 {
            coeffs = coeffs.map(|v| v.checked_mul(z(1 << 64)).expect("fits"));
        }
        let jd = Judge::new(&geom, setup.standard, &setup.notes);
        tangency_rows(
            &m,
            &faces,
            &setup.plane_ix,
            setup.n_a,
            c,
            ci,
            &coeffs,
            dir,
            surf,
            SolidSide::A,
            &jd,
        )
        .into_iter()
        .map(|t| t.reading.map(|r| r.witness))
        .collect::<Vec<_>>()
    };
    let witness = Point3::from_array([45.0, -20.0, 6.0]);
    let refused = |reason, at: Option<Point3>| {
        vec![Err(BoolError::Rejected {
            reason,
            at: at.map(crate::RejectWhere::Point),
        })]
    };
    use RejectReason::{
        CylinderGateUndecided as Undecided, DegenerateWitness, WitnessNotRational as Width,
    };
    assert_eq!(rows(0, &def), vec![Ok(witness)], "as built");
    for (steps, want) in [
        (1, refused(Undecided, Some(witness))),
        (2, refused(Width, None)),
        (3, refused(DegenerateWitness, None)),
        (4, refused(Undecided, None)),
        (5, refused(Width, None)),
    ] {
        assert_eq!(rows(steps, &def), want, "the first {steps} edits");
    }
    assert_eq!(
        rows(0, &scaled(1 << 40)),
        refused(Width, Some(witness)),
        "the lateral's width"
    );
    assert_eq!(
        rows(0, &scaled(big)),
        refused(Width, None),
        "the line on the lateral"
    );

    let mut m = Model::new();
    let (rod, slab) = crate::tests::d_rod_on_a_bitten_slab(&mut m);
    let (setup, c, ci, surf, def) = rod_against(&m, rod, slab);
    let coeffs = Judge::new(&setup.geom, setup.standard, &setup.notes)
        .world_rat(c)
        .expect("narrow")
        .map(|v| v.checked_mul(z(1 << 64)).expect("fits"));
    let mut geom = setup.geom.clone();
    geom[c].world = None;
    let jd = Judge::new(&geom, setup.standard, &setup.notes);
    let got: Vec<_> = tangency_rows(
        &m,
        &setup.planes,
        &setup.plane_ix,
        setup.n_a,
        c,
        ci,
        &coeffs,
        &def,
        surf,
        SolidSide::A,
        &jd,
    )
    .into_iter()
    .map(|t| t.reading.map(|r| r.witness))
    .collect();
    assert_eq!(
        got,
        refused(RejectReason::CurvedOperandBoundary, None),
        "the unread face first"
    );
}

/// **A tangent partner's record gathers its rows' halves of the line, joined, and refuses a join it
/// cannot order.** The block standing on `y ∈ [−20, −10]` beside the rounded plate puts a face on
/// `y = −20`, through the fillet's axis, which cuts the cylinder and holds the tangent line of
/// `x = 45`: the gate records that pair ([`SharedRuling`]). [`shared_rulings`] is asked with the
/// two records the gate keeps for it — the secant `y = −20` and the tangent `x = 45` holding the
/// line with it, checked against the gate's own record — and with the rows edited:
///
/// - two rows, one per lateral as two faces on the line would write, with `+τ` halves `[0, 1]`
///   and `[1, 2]`: one piece `[0, 2]`;
/// - `[1/p, 1]` and `[1/(p + 2), 2]` (`p` near `2⁶⁴`), whose ends `Rat` cannot subtract: `None`,
///   not an empty half, which would read as no face on that side.
#[test]
fn a_shared_line_gathers_its_rows_halves() {
    let z = Rat::from_int;
    let mut m = Model::new();
    let plate = crate::tests::rounded_plate(&mut m, 1, 0);
    let block = crate::fixtures::cuboid(&mut m, [45.0, -20.0, 0.0], [50.0, -10.0, 12.0]);
    let setup =
        crate::arrangement::plane_index_setup(&m, plate, block).expect("the gate passes it");
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let t = (0..setup.geom.len())
        .find(|&k| {
            jd.world_rat(k).is_some_and(|w| {
                w[0] != z(0)
                    && w[1] == z(0)
                    && w[2] == z(0)
                    && w[0].checked_mul(z(45)).and_then(|a| a.checked_add(w[3])) == Some(z(0))
            })
        })
        .expect("the class x = 45");
    let coeffs = jd.world_rat(t).expect("checked above");
    let ci = (0..setup.cyls.len())
        .find(|&i| {
            let d = &setup.cyls[i].def;
            nacre_exact::point_plane_clearance_rat(&coeffs, &d.origin(), d.r2())
                == nacre_exact::Orient::Zero
        })
        .expect("the fillet is tangent to x = 45");
    let holders = classes_holding_the_line(&jd, t, &coeffs, &setup.cyls[ci].def).expect("rational");
    assert_eq!(holders.len(), 1, "y = −20 holds the line: {holders:?}");
    let meets: Vec<Vec<ClassMeet>> = (0..setup.cyls.len())
        .map(|i| {
            (0..setup.geom.len())
                .map(|k| match k {
                    _ if i != ci => ClassMeet::Across,
                    _ if k == holders[0] => ClassMeet::Cuts,
                    _ if k == t => ClassMeet::Tangent {
                        holders: holders.clone(),
                        other_owner: true,
                    },
                    _ => ClassMeet::Across,
                })
                .collect()
        })
        .collect();
    let records = |rows: &[Tangency]| {
        shared_rulings(&setup.planes, &jd, &setup.cyls, &meets, rows).map(|out| {
            out.into_iter()
                .filter(|(i, _)| *i == ci)
                .map(|(_, sr)| sr)
                .collect::<Vec<_>>()
        })
    };
    let built = records(&setup.tangencies).expect("rational");
    assert_eq!(built.len(), 1, "{built:?}");
    assert_eq!(
        format!("{built:?}"),
        format!("{:?}", setup.cyls[ci].shared),
        "the gate's own record"
    );
    let row = setup
        .tangencies
        .iter()
        .find(|row| row.wall == t && row.cyl == ci)
        .expect("the block's row");
    let q = |v: Rat| nacre_exact::QuadVal::from_rat(v);
    let two_rows = |a: [Rat; 2], b: [Rat; 2]| {
        let mut rows = vec![row.clone(), row.clone()];
        rows[0].reading.halves[0] = vec![a.map(q)];
        rows[1].reading.halves[0] = vec![b.map(q)];
        records(&rows)
    };
    let joined = two_rows([z(0), z(1)], [z(1), z(2)]).expect("rational");
    let halves = joined[0]
        .tangent_halves
        .as_ref()
        .expect("a tangent partner");
    let rat = |side: &[[nacre_exact::QuadVal; 2]]| -> Vec<[Rat; 2]> {
        side.iter()
            .map(|piece| piece.map(|v| v.as_rat().expect("a rational end")))
            .collect()
    };
    assert_eq!(rat(&halves[0]), vec![[z(0), z(2)]]);
    let p = (1i128 << 64) + 1;
    assert!(
        two_rows(
            [Rat::new(1, p).expect("a ratio"), z(1)],
            [Rat::new(1, p + 2).expect("a ratio"), z(2)]
        )
        .is_none()
    );
}

/// The class `x = 1`, the D rod's cylinder and the gate's setup for one of the bitten slabs —
/// the rod is the first operand.
fn rod_against(
    m: &Model,
    rod: Handle<Solid>,
    slab: Handle<Solid>,
) -> (
    crate::arrangement::PlaneSetup,
    usize,
    usize,
    Handle<Surface>,
    nacre_topo::CylinderDef,
) {
    let z = Rat::from_int;
    let (setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(m, rod, slab).expect("setup");
    let c = {
        let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
        (0..setup.geom.len())
            .find(|&k| {
                jd.world_rat(k).is_some_and(|w| {
                    w[0] != z(0)
                        && w[1] == z(0)
                        && w[2] == z(0)
                        && w[0].checked_add(w[3]) == Some(z(0))
                })
            })
            .expect("the class x = 1")
    };
    let (ci, surf, def) = cyl_surfs
        .iter()
        .enumerate()
        .find_map(|(i, &s)| {
            let d = m.world_cylinder_def(s)?;
            (d.dir()[0] == z(0) && d.dir()[1] == z(0)).then_some((i, s, d))
        })
        .expect("the rod's cylinder");
    (setup, c, ci, surf, def)
}

/// **The line reader names what stopped it.** On the bitten slab's face `x = 1` the corner reader
/// cannot spell the bite's irrational arc ends, and the line reader refuses by the very name the
/// footprint reader gives that face (`CurvedOperandBoundary`, its boundary running along a
/// cylinder); on the rod's own flat, readable, an axis direction whose squared length leaves `Rat`
/// stops it for its width.
#[test]
fn the_line_reader_names_what_stopped_it() {
    let z = Rat::from_int;
    let mut m = Model::new();
    let (rod, slab) = crate::tests::d_rod_on_a_bitten_slab(&mut m);
    let (setup, c, _, surf, def) = rod_against(&m, rod, slab);
    let coeffs = Judge::new(&setup.geom, setup.standard, &setup.notes)
        .world_rat(c)
        .expect("narrow");
    let foot = tangency_foot(&coeffs, &def.origin()).expect("rational");
    let face_of = |slab_side: bool| {
        (0..setup.planes.len())
            .find_map(|i| match &setup.planes[i] {
                FaceRow::Plane(p)
                    if setup.plane_ix[i] == ClassIx::Plane(c) && (i >= setup.n_a) == slab_side =>
                {
                    p.face
                }
                _ => None,
            })
            .expect("a face on x = 1")
    };
    let bitten = m.face(face_of(true));
    assert!(
        bitten.outer.half_edges.iter().any(|he| matches!(
            corner_of(&m, bitten, &bitten.outer, he),
            Err(CornerFail::Shape)
        )),
        "the corner reader cannot spell the bitten face"
    );
    let spans = lateral_spans(&setup.planes, surf);
    assert!(matches!(
        face_clears_footprint(
            &m,
            bitten,
            &coeffs,
            &def.origin(),
            &def.dir(),
            def.r2(),
            &spans
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::CurvedOperandBoundary,
            at: None,
        })
    ));
    let got = line_runs_through_face(&m, bitten, &coeffs, &foot, &def.dir(), [z(0), z(2)]);
    assert_eq!(got.err(), Some(RejectReason::CurvedOperandBoundary));
    let big = Rat::from_int(1 << 64);
    let flat = m.face(face_of(false));
    assert!(line_runs_through_face(&m, flat, &coeffs, &foot, &def.dir(), [z(0), z(2)]).is_ok());
    let got = line_runs_through_face(&m, flat, &coeffs, &foot, &[z(0), big, big], [z(0), z(2)]);
    assert_eq!(got.err(), Some(RejectReason::WitnessNotRational));
}

/// **A lateral's sides of a ruling stop for width where a value leaves `Rat`.** The cylinder of
/// [`the_tangent_wall_states_itself_exactly`] stands at `x = 1/2`; a foot at `x = (2¹²⁷ − 1)/3`
/// cannot be taken off the axis.
#[test]
fn a_lateral_read_stops_for_its_width() {
    let mut m = Model::new();
    let cyl = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.5, 1.0, -1.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        4.0,
    );
    let lateral = cyl.lateral;
    let def = m
        .world_cylinder_def(m.face(lateral).surface)
        .expect("a world cylinder");
    let z = Rat::from_int;
    let foot = [Rat::new(i128::MAX, 3).expect("a value"), z(1), z(0)];
    assert_eq!(
        lateral_sides_on_ruling(&m, m.face(lateral), &def, &foot).err(),
        Some(RejectReason::WitnessNotRational)
    );
}

/// **A wall face read through a sibling refuses by the sibling's own refusal.** The split bitten
/// slab's two faces on `x = 1`: the lower runs through the rod's line and reads; the upper carries
/// the bite, stands above the lateral, and cannot be read. The lower face's row reads the upper
/// to place the slab's sides of the line, and it carries the upper face's footprint refusal — the
/// very value the upper face's own row carries — so whichever row the gate reaches first, the
/// refusal is one.
#[test]
fn a_sibling_wall_face_refuses_a_row_by_its_own_refusal() {
    let mut m = Model::new();
    let (rod, slab) = crate::tests::d_rod_on_a_split_bitten_slab(&mut m);
    let (setup, c, ci, surf, def) = rod_against(&m, rod, slab);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let coeffs = jd.world_rat(c).expect("narrow");
    let spans = lateral_spans(&setup.planes, surf);
    let reads: Vec<_> = (setup.n_a..setup.planes.len())
        .filter(|&i| setup.plane_ix[i] == ClassIx::Plane(c))
        .map(|i| {
            let face = m.face(setup.planes[i].plane().face.expect("a face"));
            face_clears_footprint(
                &m,
                face,
                &coeffs,
                &def.origin(),
                &def.dir(),
                def.r2(),
                &spans,
            )
        })
        .collect();
    let unread = BoolError::Rejected {
        reason: RejectReason::CurvedOperandBoundary,
        at: None,
    };
    assert_eq!(reads.len(), 2, "{reads:?}");
    assert!(
        reads.contains(&Ok(false)) && reads.contains(&Err(unread)),
        "{reads:?}"
    );
    let rows = tangency_rows(
        &m,
        &setup.planes,
        &setup.plane_ix,
        setup.n_a,
        c,
        ci,
        &coeffs,
        &def,
        surf,
        SolidSide::A,
        &jd,
    );
    assert_eq!(rows.len(), 2, "one row per slab face: {rows:?}");
    for t in &rows {
        assert_eq!(t.reading.as_ref().err(), Some(&unread), "{t:?}");
    }
}
