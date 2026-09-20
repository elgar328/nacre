use super::*;
use nacre_scalar::{Angle, Axis, Rat};

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
        nacre_topo::CylinderDef::new(o, m, e, nacre_scalar::BigRat::from(q(1, 25))).unwrap() // radius 1/5, as r²
    };
    let a = cyl([z(0); 3], [z(0), z(0), z(1)], [z(1), z(0), z(0)]);
    let fp = |span: [i128; 2]| Footprint {
        span: Some([z(span[0]), z(span[1])]),
        theta: None,
    };
    let (x, dirs_of) = (fp([0, 2]), |b: &nacre_topo::CylinderDef| {
        separating_dirs(&a, b)
    });
    // `B` along `x`: its origin's `y`/`z` and its own span place it.
    let along_x =
        |y: i128, zz: i128| cyl([z(0), z(y), z(zz)], [z(1), z(0), z(0)], [z(0), z(0), z(1)]);
    let cases: [(nacre_topo::CylinderDef, [i128; 2], [bool; 3]); 4] = [
        // Clear of `A`'s span along `A`'s own axis, and nothing else.
        (along_x(0, 5), [-10, 10], [true, false, false]),
        // Off the end of `B`'s own span, and nothing else.
        (along_x(0, 1), [5, 10], [false, true, false]),
        // ★ Beside it, across both axes — only the common perpendicular sees this.
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
    assert_eq!(separating_dirs(&a, &parallel).len(), 2);
}

/// ★★ **The third direction is the face-level twin of the surface rung.** With a whole circle
/// and no span, `separated` along the common perpendicular says exactly what
/// [`nacre_scalar::cylinders_clear`] says about the infinite surfaces — so the face test does
/// not disagree with the rung that runs before it, it only knows more when a face knows more.
#[test]
fn the_common_perpendicular_degenerates_to_the_surface_rule() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let cyl = |o: [Rat; 3], m: [Rat; 3], e: [Rat; 3]| {
        nacre_topo::CylinderDef::new(o, m, e, nacre_scalar::BigRat::from(q(1, 25))).unwrap() // radius 1/5, as r²
    };
    let a = cyl([z(0); 3], [z(0), z(0), z(1)], [z(1), z(0), z(0)]);
    let unbounded = Footprint {
        span: None,
        theta: None,
    };
    // Offsets either side of the radius sum `2/5`, and exactly on it.
    for (n, d) in [(1i128, 2i128), (2, 5), (1, 3), (9, 10)] {
        let b = cyl(
            [z(0), q(n, d), z(0)],
            [z(1), z(0), z(0)],
            [z(0), z(0), z(1)],
        );
        let perp = separating_dirs(&a, &b)[2];
        let face = separated(&a, &unbounded, &b, &unbounded, &perp);
        let surface = nacre_scalar::cylinders_clear(
            &a.origin(),
            &a.dir(),
            a.r2(),
            &b.origin(),
            &b.dir(),
            b.r2(),
        ) == nacre_scalar::Orient::Positive;
        assert_eq!(face, Some(surface), "offset {n}/{d}");
    }
}

/// ★★★★★ **The same circle and the same ruling, and the arc decides.**
///
/// The unit-ish circle of radius `6/5` about the origin in `z = 0`, and a cylinder along `x`
/// through the origin with radius `1` — so that plane cuts it in the two rulings `y = ±1`, and
/// the extent every piece is measured along is `e = n × m = ŷ`.
///
/// The whole circle reaches `y ∈ [−6/5, 6/5]` and holds both rulings. The **left** arc — from
/// `(−24/25, 18/25)` to `(−24/25, −18/25)`, counter-clockwise about `+ẑ`, so it runs through
/// `−x̂` and holds neither `±ŷ` — reaches only `[−18/25, 18/25]` and holds neither. A **top**
/// arc holds `+ŷ`, so it reaches `+6/5` and holds the `+1` ruling but not the `−1`.
///
/// Asking the circle instead would refuse a plate whose
/// corner fillets never come near the drill's rulings.
#[test]
fn the_arc_decides_which_ruling_a_circle_holds() {
    use nacre_scalar::{MeetPoint, cylinder_ruling_reached_extent};
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let coeffs = [z(0), z(0), z(1), z(0)]; // the plane `z = 0`
    let axis = [z(0), z(0), z(1)]; // the circle's carrier
    let (o, m, r) = (
        [z(0); 3],
        [z(1), z(0), z(0)],
        nacre_scalar::BigRat::from(z(1)),
    ); // the ruling cylinder, along `x`
    let e = [z(0), z(1), z(0)]; // `n × m`
    let centre = [z(0); 3];
    let rho2 = nacre_scalar::BigRat::from(q(36, 25)); // ρ = 6/5, stated as its square
    let holds = |arc: Option<&RimArc>, side: i8| {
        let (lo, hi) = arc_ends_along(&centre, &rho2, &axis, arc, &e).expect("an extent");
        cylinder_ruling_reached_extent(
            &coeffs,
            &nacre_scalar::StripReach {
                lo: (&MeetPoint::Narrow(lo.0), &lo.1),
                hi: Some((&MeetPoint::Narrow(hi.0), &hi.1)),
            },
            &o,
            &m,
            &r,
            side,
            // Production's question: a crossing, not a touch. Every case below is
            // clear of the boundary either way, which is why the boundary has its own lock.
            false,
        )
    };
    // The whole circle spans `±6/5` and holds both.
    for side in [-1i8, 1] {
        assert!(holds(None, side), "the whole circle, side {side}");
    }
    // The left arc spans `±18/25` and holds neither — the arc's own reach, not its circle's.
    let left = RimArc {
        from: [q(-24, 25), q(18, 25), z(0)],
        to: [q(-24, 25), q(-18, 25), z(0)],
    };
    for side in [-1i8, 1] {
        assert!(!holds(Some(&left), side), "the left arc, side {side}");
    }
    // A top arc holds `+ŷ`, so it reaches `+6/5` — and still never the other ruling.
    let top = RimArc {
        from: [q(24, 25), q(18, 25), z(0)],
        to: [q(-24, 25), q(18, 25), z(0)],
    };
    assert!(holds(Some(&top), 1), "the top arc reaches the `+1` ruling");
    assert!(!holds(Some(&top), -1), "and never the `−1` one");
}

/// ★★★★★ **An arc's reach, by hand, and which arc it is.**
///
/// The unit circle about the origin in the plane `z = 0`, axis `+ẑ`, and the **first quadrant**
/// as `from = (1,0) → to = (0,1)` (counter-clockwise about that axis). Its points are
/// `(cos θ, sin θ)` for `θ ∈ [0, π/2]`, so along `d = (1,1,0)` it reaches
/// `cos θ + sin θ ∈ [1, √2]` — **one end rational, the other a radical**, which is the shape
/// no symmetric margin can state and the reason this rule exists.
///
/// ★ The last block is the lock that matters most: swapping the ends names the **complementary**
/// arc, whose reach is genuinely different. Reading the face's normal instead of the carrier's
/// axis would make exactly that swap.
#[test]
fn an_arcs_reach_is_its_own_and_not_its_complements() {
    let z = |v: i128| Rat::from_int(v);
    let one = nacre_scalar::BigRat::from(z(1)); // r² = 1
    let axis = [z(0), z(0), z(1)];
    let quadrant = RimArc {
        from: [z(1), z(0), z(0)],
        to: [z(0), z(1), z(0)],
    };
    // Along `+x̂`: the peak is on the arc (it is `from`), so that end carries the radical and
    // the other stops at an end of the arc — `[0, 1]`, the quadrant's own span in `x`.
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(0), z(0), z(1)))
    );
    // Along the diagonal: `[1, √2]`. `ρ² = r²|d⊥|² = 2`, and the low end is the rational `1`.
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &[z(1), z(1), z(0)]),
        Some((z(1), z(0), z(0), z(2)))
    );
    // A whole circle reaches alike both ways, and `d ∥ axis` has no radial term at all.
    assert_eq!(
        arc_extent(None, &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(1), z(0), z(1)))
    );
    assert_eq!(
        arc_extent(Some(&quadrant), &one, &axis, &axis),
        Some((z(0), z(0), z(0), z(0)))
    );
    // ★ The complement is a different set: `[−1, 1]`, not `[0, 1]`.
    let rest = RimArc {
        from: quadrant.to,
        to: quadrant.from,
    };
    assert_eq!(
        arc_extent(Some(&rest), &one, &axis, &[z(1), z(0), z(0)]),
        Some((z(0), z(1), z(0), z(1)))
    );
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
        nacre_scalar::BigRat::from(q(1, 25)), // radius 1/5, as r²
    )
    .unwrap();
    let fp = |span: Option<[Rat; 2]>| Footprint { span, theta: None };
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

    // ★ **An arc reaches less than its circle.** The quarter arc from `(−r, 0)` to
    // `(0, −r)` — the third quadrant, counter-clockwise about `+z`.
    let m = [z(0), z(0), z(1)];
    let arc = RimArc {
        from: [q(-1, 5), z(0), z(0)],
        to: [z(0), q(-1, 5), z(0)],
    };
    assert_eq!(arc_contains(&arc, &[z(-1), z(-1), z(0)], &m), Some(true)); // 225°
    assert_eq!(arc_contains(&arc, &arc.from, &m), Some(true)); // an end is on it
    assert_eq!(arc_contains(&arc, &[z(0), z(1), z(0)], &m), Some(false)); // 90°
    assert_eq!(arc_contains(&arc, &[z(1), z(-1), z(0)], &m), Some(false)); // 315°
    let long = RimArc {
        from: arc.to,
        to: arc.from,
    }; // the other three quarters
    assert_eq!(arc_contains(&long, &[z(1), z(1), z(0)], &m), Some(true));
    assert_eq!(arc_contains(&long, &[z(-1), z(-1), z(0)], &m), Some(false));
    let half = RimArc {
        from: [z(1), z(0), z(0)],
        to: [z(-1), z(0), z(0)],
    }; // exactly a half turn, through +y
    assert_eq!(arc_contains(&half, &[z(0), z(1), z(0)], &m), Some(true));
    assert_eq!(arc_contains(&half, &[z(0), z(-1), z(0)], &m), Some(false));
    let fq = |span: Option<[Rat; 2]>| Footprint {
        span,
        theta: Some(arc),
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

/// **The tangency the gate used to refuse is now *stated*** — and every field of that
/// statement is asserted here, because the verdict downstream is only as good as this row.
///
/// The frozen census shape (`cylinder-wall-tangent`): a `2³` cube and a cylinder of radius
/// `0.5` whose axis stands at `x = 0.5`, so the wall `x = 0` is **exactly** `r` away. Hand
/// geometry throughout — nothing below is copied from an engine run.
///
/// ★ Called directly rather than through the gate's ledger: a process-global read by tests
/// running in parallel attributes one fixture's rows to another
/// ([[nondeterministic-fixtures-and-instruments]]), and the row builder is a pure function of
/// the setup, so the honest reading is to hand it that setup.
#[test]
fn the_tangent_wall_states_itself_exactly() {
    let mut m = Model::new();
    let cube = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let cyl = m.add_cylinder(
        Point3::from_array([0.5, 1.0, -1.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        4.0,
    );
    m.rebuild_adjacency();
    let (setup, cyl_surfs) = plane_index_setup_inner(&m, cube, cyl).expect("setup");
    let surf = cyl_surfs[0];
    let def = world_cylinder_def(&m, surf).expect("a world cylinder");
    let (o, r2) = (def.origin(), def.r2());
    // The tangent class, found the gate's own way: one clearance call per class.
    let c = (0..setup.geom.len())
        .find(|&k| {
            setup.geom[k].world_rat.is_some_and(|w| {
                nacre_scalar::point_plane_clearance_rat(&w, &o, r2) == nacre_scalar::Orient::Zero
            })
        })
        .expect("the wall x = 0 is exactly r from the axis");
    let coeffs = setup.geom[c].world_rat.expect("checked above");
    let rows = tangency_rows(
        &m,
        &setup.planes,
        &setup.plane_ix,
        &setup.geom,
        setup.n_a,
        c,
        0,
        &coeffs,
        &def,
        surf,
        SolidSide::B, // the cylinder is the second operand; the wall face is the cube's
    );
    assert_eq!(rows.len(), 1, "one wall face, one lateral face: {rows:?}");
    let t = &rows[0];
    assert_eq!((t.wall, t.cyl), (c, 0));
    assert_eq!(
        (t.wall_solid, t.cyl_solid),
        (SolidSide::A, SolidSide::B),
        "the cube is the first operand"
    );
    // The cube's material is `x > 0` and the whole cylinder stands there too.
    assert!(t.lens_in_wall_solid, "{t:?}");
    // A solid cylinder operand keeps its material inside.
    assert_eq!(t.cyl_orient, 1, "{t:?}");
    // The face `x = 0` runs from `y = 0` to `y = 2` across the tangent line `y = 1`, and its
    // corners sit at `t = 1` and `t = 3` inside the lateral's span `[0, 4]`.
    assert!(t.straddles, "{t:?}");
    // No other plane of either operand holds the line `x = 0, y = 1` (a plane that did would
    // have to contain `+Z` and pass through `(0, 1)`).
    assert!(!t.line_in_another_plane, "{t:?}");
    assert!(!t.undecided, "{t:?}");
    // The foot of the perpendicular is `(0, 1, −1)`; the span's middle carries it to `t = 2`.
    let w = t.witness.as_array();
    assert!(
        (w[0] - 0.0).abs() < 1e-12 && (w[1] - 1.0).abs() < 1e-12 && (w[2] - 1.0).abs() < 1e-12,
        "the witness sits on the contact, mid-span: {w:?}"
    );
}

/// **The pierce tolerance is a measurement, and each of its terms can move.**
///
/// Hand geometry — no class indices, so nothing here is copied from an engine run: the turned
/// boss's second crossing, where the meet line of `y = 0` and `x = 4` pierces a cylinder of
/// radius `0.5` about the `+X` axis through `(y, z) = (0.25, 2)`. The true point is
/// `(4, 0, 2 − √3/4)` (irrational on purpose — every term sees real rounding, not exact zeros).
///
/// ★ The second probe moves **along the meet line**: both planes and the line stay at zero, so
/// only the cylinder-surface term can see it — which is how this asserts that term is
/// *exercised*, not merely present ([[instrument-decides-the-answer]]'s negative control).
#[test]
fn pierce_vertex_tol_measures_and_each_term_moves() {
    let pa = Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .unwrap();
    let pb = Plane::from_point_normal(
        Point3::from_array([4.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
    )
    .unwrap();
    let cyl = nacre_geom::Cylinder::from_axis(
        Point3::from_array([0.0, 0.25, 2.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        0.5,
    )
    .unwrap();
    let s = 2.0 - 3.0f64.sqrt() / 4.0;
    let p = Point3::from_array([4.0, 0.0, s]);
    assert!(
        pierce_vertex_tol(p, &pa, &pb, &cyl) < 1e-12,
        "the true crossing measures at rounding scale: {}",
        pierce_vertex_tol(p, &pa, &pb, &cyl)
    );
    // Along the meet line: planes and line stay zero, the cylinder term alone answers.
    let along = Point3::from_array([4.0, 0.0, s - 1e-6]);
    assert!(
        pierce_vertex_tol(along, &pa, &pb, &cyl) > 1e-7,
        "the cylinder-surface term is exercised"
    );
    // Off a plane: the instrument moves by the full offset.
    let off = Point3::from_array([4.0, 1e-6, s]);
    assert!(
        pierce_vertex_tol(off, &pa, &pb, &cyl) > 0.9e-6,
        "a plane term is exercised"
    );
}

/// **A loop with points that do not turn still states its own outward direction.**
///
/// Two faces sharing an edge list the same vertices along it, so a pad that covers part of a
/// face leaves the neighbouring walls with points strung along one line. Those points cannot
/// be removed — that would leave a T-vertex — so the population is every partially-covering
/// pad, and `outer_tri` has to survive it.
///
/// ★ **Turned coordinates are the whole difficulty.** Axis-aligned collinear points cancel to
/// exactly `0.0`, which the old "first corner that is not flat" rule skipped correctly. Turned,
/// the cancellation leaves ~2⁻⁵³ — not zero, so it was taken, and the direction that came back
/// was rounding: measured 90° off its own plane, which is `orient_sign`, the winding of the
/// witness triangle every predicate borrows, and which side of the plane holds material.
///
/// Asserted on `collect_planes`' own output rather than through `debug_assert`, so it is a
/// measurement in release too. `tests/probes/collinear_loop_points.rs` sweeps the same proposition
/// from outside the crate over a band of angles.
#[test]
fn a_loop_whose_points_do_not_all_turn_still_faces_the_right_way() {
    use crate::{OpOutput, Operation, Profile2d, SketchFrame, apply};
    use nacre_math::Point2;
    use nacre_scalar::{Isometry, Rotation};

    let rect = |a: f64, b: f64, c: f64, d: f64| {
        Profile2d::polygon(vec![
            Point2::from_array([a, b]),
            Point2::from_array([c, b]),
            Point2::from_array([c, d]),
            Point2::from_array([a, d]),
        ])
        .expect("a rectangle is a fair profile")
    };

    let mut m = Model::new();
    let world = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: world,
            profile: rect(0.0, 0.0, 2.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("the block") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let OpOutput::Transform { solid } = apply(
        &mut m,
        &Operation::Transform {
            solid,
            isometry: Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
            }),
        },
    )
    .expect("the turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    // A pad covering part of the face: the walls around it inherit the split edge.
    let face = m.shell(m.solid(solid).outer).faces[0];
    let OpOutput::PadOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PadOnFace {
            face,
            profile: rect(0.25, 0.25, 1.25, 1.25),
            dist: 1.0,
        },
    )
    .expect("the pad") else {
        unreachable!()
    };
    m.rebuild_adjacency();

    let faces = collect_planes(&m, solid).expect("the planes of a padded block");
    // The population has to be present, or this asserts over nothing.
    let straight = faces
        .iter()
        .filter(|f| {
            let n = m
                .face(f.face().expect("a model face"))
                .outer
                .half_edges
                .len();
            n > 4
        })
        .count();
    assert!(
        straight > 0,
        "no face carries a split edge — the fixture stopped reaching the population"
    );
    // The lock is on the *winding* side now: `n_out` is read off the stored
    // orientation, so `plane.normal()·n_out` is ±1 by construction and asserts
    // nothing. What the widest-corner rule still owns is the triangle — revert
    // `outer_tri` to "first non-flat corner" and the cross below is rounding
    // noise again, failing this in release.
    for row in &faces {
        let f = row.plane();
        let cos = (f.tri[1] - f.tri[0])
            .cross(f.tri[2] - f.tri[0])
            .normalize()
            .expect("a widest corner spans area")
            .dot(f.n_out);
        assert!(
            cos > 0.5,
            "a face's outer triangle is {cos:.3} of the way to perpendicular against its \
                 stated outward"
        );
    }
}

/// **A judgement's headroom is relative to its model, not carved out of a shared ceiling.**
///
/// The two budgets answer different questions — how deep the model is, and how thin a witness
/// it may still judge — and sharing one absolute number silently couples them: a solid turned
/// a thousand times would leave a hard judgement no room, so the same sliver would be judged
/// in a fresh model and abandoned in a turned one. Turning a model must not change what
/// counts as judgeable, so this pins the headroom to a constant *above the model's own
/// precision*, at both ends of the depth range.
#[test]
fn the_climbing_headroom_survives_a_deep_model() {
    let deg = Angle::from_deg(Rat::from_int(37)).expect("angle");
    let mut p = WitnessPoint::at([Rat::from_int(1), Rat::from_int(2), Rat::from_int(3)]);
    let mut seen = Vec::new();
    for turn in 0..=200 {
        if turn == 0 || turn == 20 || turn == 200 {
            let j = standard_for_points(std::slice::from_ref(&p));
            assert_eq!(
                j.cap,
                j.prec + CLIMB_HEADROOM,
                "turn {turn}: headroom is not the model's own precision plus a constant"
            );
            seen.push(j.prec);
        }
        p = p.rotate(Axis::Z, deg);
    }
    // …and the precision really did grow with the history, or the test above is vacuous.
    assert!(
        seen[0] < seen[2],
        "precision did not grow with the rotation history: {seen:?}"
    );
}
