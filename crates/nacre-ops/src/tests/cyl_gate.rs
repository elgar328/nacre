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

/// **The tangency is *stated*** — and every field of that
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
    let cube = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let cyl = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([0.5, 1.0, -1.0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        4.0,
    )
    .solid;
    m.rebuild_adjacency();
    let (setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(&m, cube, cyl).expect("setup");
    let surf = cyl_surfs[0];
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    let (o, r2) = (def.origin(), def.r2());
    // The tangent class, found the gate's own way: one clearance call per class.
    let c = (0..setup.geom.len())
        .find(|&k| {
            setup.geom[k].world_rat().is_some_and(|w| {
                nacre_exact::point_plane_clearance_rat(&w, &o, r2) == nacre_exact::Orient::Zero
            })
        })
        .expect("the wall x = 0 is exactly r from the axis");
    let coeffs = setup.geom[c].world_rat().expect("checked above");
    // No other plane of either operand holds the line `x = 0, y = 1` (a plane that did would
    // have to contain `+Z` and pass through `(0, 1)`).
    let holders = classes_holding_the_line(&setup.geom, c, &coeffs, &def).expect("exact");
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
        !holders.is_empty(),
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
    assert!(!t.line_in_another_plane, "{t:?}");
    assert!(t.undecided.is_none(), "{t:?}");
    // The foot of the perpendicular is `(0, 1, −1)`; the span's middle carries it to `t = 2`.
    let w = t.witness.as_array();
    assert!(
        (w[0] - 0.0).abs() < 1e-12 && (w[1] - 1.0).abs() < 1e-12 && (w[2] - 1.0).abs() < 1e-12,
        "the witness sits on the contact, mid-span: {w:?}"
    );
}

/// **Whether a line along the axis lies within a lateral face's angular extent** — the half turn
/// `(1, 0) → (−1, 0)` counter-clockwise about `+z`, the upper half: a line through `(0, 1)` is on
/// it, through `(0, −1)` is not, and both ends are (a line on a face's end ruling is its edge).
/// The height the point is given at does not matter, and a whole circle holds every line.
#[test]
fn a_faces_angular_extent_holds_the_lines_on_it_ends_included() {
    let z = Rat::from_int;
    let upper = Footprint {
        span: None,
        theta: Some(RimArc {
            from: [z(1), z(0), z(0)],
            to: [z(-1), z(0), z(0)],
        }),
    };
    let (o, m) = ([z(0); 3], [z(0), z(0), z(1)]);
    for (p, want) in [
        ([z(0), z(1), z(5)], true),
        ([z(0), z(-1), z(0)], false),
        ([z(1), z(0), z(3)], true),
        ([z(-1), z(0), z(-2)], true),
        ([z(1), z(-1), z(0)], false),
    ] {
        assert_eq!(upper.theta_holds_line(&p, &o, &m), Some(want), "{p:?}");
    }
    let whole = Footprint {
        span: None,
        theta: None,
    };
    assert_eq!(
        whole.theta_holds_line(&[z(0), z(-1), z(0)], &o, &m),
        Some(true)
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
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([lo[0], lo[1], -0.5]),
        Point3::from_array([hi[0], hi[1], 2.5]),
    );
    m.rebuild_adjacency();
    crate::arrangement::plane_index_setup(&m, a, b).map(|s| s.tangencies)
}

/// ★ **A tangency row is about the lateral face, not its surface.** The half cylinder's surface
/// is tangent to three box walls; the gate writes a row only where the face is:
///
/// - `y = −1` touches the surface at `(0, −1)`, on the half the lateral does not have — no row;
/// - `y = 1` touches it at `(0, 1)`, on the lateral — one row, decided;
/// - `x = 1` touches it on the lateral's **end** ruling `(1, 0)` — a row (ends are the face's
///   edge), held by a third plane too: the half cylinder's own flat `y = 0`. The box's wall runs
///   across that line rather than ending on it, so the line is no edge of both faces, the verdict
///   could not read the six regions there, and the gate refuses the pair by that name.
#[test]
fn a_tangency_row_is_written_only_where_the_lateral_face_is() {
    let missing = half_cylinder_and_box([-0.5, -1.0], [0.5, 2.0]).expect("setup");
    assert!(missing.is_empty(), "{missing:?}");
    let present = half_cylinder_and_box([-0.5, 1.0], [0.5, 2.0]).expect("setup");
    assert_eq!(present.len(), 1, "{present:?}");
    assert!(
        present[0].undecided.is_none() && !present[0].line_in_another_plane,
        "{present:?}"
    );
    match half_cylinder_and_box([1.0, -1.0], [2.0, 2.0]) {
        Err(BoolError::Rejected { reason, .. }) => {
            assert_eq!(reason, RejectReason::TangentLineInAnotherPlane);
        }
        other => panic!("the end ruling's row is refused at the gate: {other:?}"),
    }
}
