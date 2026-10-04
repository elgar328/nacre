//! The cylinder gates in scalar form: ruling reach, strip side, the cylinder-pair gate — and the
//! push funnel that realizes a vertex or keeps its fallback.

use super::*;

/// ★★★★★ **The negative control — the net can see a crossing.**
///
/// The audit reports zero circle–ruling crossings over the corpus, and that number means
/// «none there» only if the instrument can say «there». This drives the rational net directly,
/// where the geometry is written by hand rather than found: a class `x = 0`, a cylinder about
/// `z` at `x = 3` with radius `5` — so its two rulings stand at `y = ±4` — and a circle centred
/// on that plane whose radius is walked across each of them.
///
/// ★★★★★ **The extent form and the disk form are one door.**
///
/// [`nacre_exact::cylinder_ruling_reached_extent`] answers for a piece stated by its two ends;
/// a disk is that piece with both ends the same. Handing it the **same** point and margin twice
/// must reproduce the disk answer exactly — otherwise the general form has drifted from the case
/// it generalizes, and every arc reading rests on it. Swept over centres, radii and both sides so
/// the agreement is a population, not a sample.
///
/// ★ And the last block is what the generalization is *for*: two **different** ends read a reach
/// the symmetric form cannot state, and they answer differently.
#[test]
fn the_extent_form_agrees_with_the_disk_it_generalizes() {
    use nacre_exact::{MeetPoint, Rat, cylinder_ruling_reached, cylinder_ruling_reached_extent};
    let r = |v: i128| Rat::from_int(v);
    let b = |v: i128| nacre_exact::BigRat::from(Rat::from_int(v)); // a squared radius or margin
    let wall = [r(1), r(0), r(0), r(0)];
    let axis = [r(0), r(0), r(1)];
    for ox in [3i128, 5, 7] {
        for rad in [4i128, 5, 9] {
            for y in [-10i128, -6, -4, -1, 0, 1, 4, 6, 10] {
                for rho in [0i128, 1, 2, 4, 7] {
                    for side in [-1i8, 0, 1] {
                        let p = MeetPoint::Narrow([r(0), r(y), r(0)]);
                        let o = [r(ox), r(0), r(0)];
                        let disk = cylinder_ruling_reached(
                            &wall,
                            &p,
                            &b(rho * rho),
                            &o,
                            &axis,
                            &b(rad * rad),
                            side,
                        );
                        let both = cylinder_ruling_reached_extent(
                            &wall,
                            &nacre_exact::StripReach {
                                lo: (&p, &b(rho * rho)),
                                hi: Some((&p, &b(rho * rho))),
                            },
                            &o,
                            &axis,
                            &b(rad * rad),
                            side,
                            true,
                        );
                        assert_eq!(disk, both, "ox {ox} rad {rad} y {y} rho {rho} side {side}");
                        // ★★ **Crossing implies reaching, everywhere.** The open
                        // reading is the closed one minus its boundary, so it can never say
                        // *more*. This is what sweeps the hand-written arms of the sign table
                        // the named boundary threads through; nothing else exercises them in bulk.
                        let crossed = cylinder_ruling_reached_extent(
                            &wall,
                            &nacre_exact::StripReach {
                                lo: (&p, &b(rho * rho)),
                                hi: Some((&p, &b(rho * rho))),
                            },
                            &o,
                            &axis,
                            &b(rad * rad),
                            side,
                            false,
                        );
                        assert!(
                            !crossed || both,
                            "crossed but not reached: ox {ox} rad {rad} y {y} rho {rho} side {side}"
                        );
                    }
                }
            }
        }
    }
    // ★ A reach from `y = 4` to `y = 6` with no margin: it holds the `−` ruling at `y = 4`
    // (closed, a touch counts) and misses the `+` one at `y = −4` entirely — which no single
    // centre-and-radius could say, because a disk covering `[4, 6]` is centred at `5`.
    let lo = MeetPoint::Narrow([r(0), r(4), r(0)]);
    let hi = MeetPoint::Narrow([r(0), r(6), r(0)]);
    let o = [r(3), r(0), r(0)];
    let reach = |side, touch_counts| {
        cylinder_ruling_reached_extent(
            &wall,
            &nacre_exact::StripReach {
                lo: (&hi, &b(0)),
                hi: Some((&lo, &b(0))),
            },
            &o,
            &axis,
            &b(25), // radius 5, as r²
            side,
            touch_counts,
        )
    };
    // The scalar side runs along `−y` here, so the ruling at `y = +4` is the `−` one.
    assert!(
        reach(-1, true),
        "the reach ends on the ruling, and a touch counts"
    );
    assert!(!reach(1, true), "and never comes near the other one");
    // ★★ **And the same reach does not *cross* it.** The arrangement's net asks about a
    // crossing, because a node no road mints is one; an edge tangent to another divides nothing,
    // which is the same sentence as for a tangency at a vertex.
    assert!(!reach(-1, false), "ending on the ruling is not crossing it");
    assert!(!reach(1, false));
    // ★ And for a **disk**, where the two readings are the only thing separating a tangency from
    // an overlap: the same rulings at `y = ±4`, and a disk covering `[4, 8]` about `y = 6`.
    let tangent = MeetPoint::Narrow([r(0), r(6), r(0)]);
    let touching = |touch_counts| {
        cylinder_ruling_reached_extent(
            &wall,
            &nacre_exact::StripReach {
                lo: (&tangent, &b(4)), // ρ = 2, as ρ²
                hi: None,
            },
            &o,
            &axis,
            &b(25),
            -1,
            touch_counts,
        )
    };
    assert!(touching(true), "the disk ends exactly on the ruling");
    assert!(!touching(false), "and does not cross it");
}

/// It also pins the two answers a *side* separates (a circle reaching `y = +4` does not reach
/// `y = −4`), the tangent-wall spelling (`side == 0`, the ruling at `y = 0`), a plane that clears
/// the cylinder outright, and the precondition (a class not parallel to the axis is `None`, never
/// `false`).
#[test]
#[cfg(debug_assertions)]
fn a_circle_and_a_ruling_meet_or_clear() {
    use nacre_exact::{MeetPoint, Rat, cylinder_ruling_reached as reaches};
    let r = |v: i128| Rat::from_int(v);
    let b = |v: i128| nacre_exact::BigRat::from(Rat::from_int(v)); // a squared radius or margin
    // The class `x = 0`; a cylinder about `z` at `x = 3` with radius 5, so its two rulings stand
    // at `y = ±4`; and a circle centred on that plane whose radius is walked across each of them.
    let wall = [r(1), r(0), r(0), r(0)];
    let axis = [r(0), r(0), r(1)];
    let at = |x: i128| [r(x), r(0), r(0)];
    let centre = |y: i128| MeetPoint::Narrow([r(0), r(y), r(0)]);
    let meets = |ox: i128, rad: i128, y: i128, rho: i128, side: i8| {
        // The doors take the radii as squares.
        reaches(
            &wall,
            &centre(y),
            &b(rho * rho),
            &at(ox),
            &axis,
            &b(rad * rad),
            side,
        )
    };
    // ★ Sides here are the **scalar** family's: the sign of `(p − o)·(n × m)`, which for this
    // wall and axis points along `−y`. So a circle centred at `y = 5` sits on the `−` side, and
    // the ruling it reaches is the `−` one.
    assert!(meets(3, 5, 5, 2, -1));
    // …and not the `+` one, which is nine away.
    assert!(!meets(3, 5, 5, 2, 1));
    // The mirror: a centre at `y = −5` reaches the `+` ruling.
    assert!(meets(3, 5, -5, 2, 1));
    // Far enough away and it reaches neither.
    for side in [-1i8, 1] {
        assert!(!meets(3, 5, 10, 2, side), "side {side}");
    }
    // A tangent wall: one ruling at `y = 0`, `side == 0`.
    assert!(meets(5, 5, 1, 2, 0));
    assert!(!meets(5, 5, 5, 2, 0));
    // A plane that clears the cylinder carries no ruling to meet.
    assert!(!meets(10, 5, 0, 2, 1));
    // ★ The precondition (a class parallel to the axis) is the caller's to check — the arrangement
    // wrapper asks it and answers `None`; here it is a `debug_assert`, so this test does not drive
    // an off-axis class through the door.
}

/// ★★★★★ **A disk clears the strip, spans it, or floats inside it, and the three are
/// told apart.** The gate's disk arm rests on this, and a predicate that answered one of them
/// always would still pass every model in the corpus.
///
/// The class `x = 0`; a cylinder about `z` at `x = 3` with radius 5, so the strip runs `|U| ≤ 4`
/// in the `n × m` coordinate — which for this wall and axis points along `−y`.
///
/// ★ The last case is the one the vocabulary exists for: a **point** on a ruling spans nothing, so
/// `Crosses` needs width and never comes back from the spellings a vertex takes.
#[test]
fn a_disk_clears_a_strip_or_spans_it() {
    use nacre_exact::{MeetPoint, Rat, StripSide, cylinder_strip_side_margin as side};
    let r = |v: i128| Rat::from_int(v);
    let b = |v: i128| nacre_exact::BigRat::from(Rat::from_int(v)); // a squared radius or margin
    let wall = [r(1), r(0), r(0), r(0)];
    let axis = [r(0), r(0), r(1)];
    let at = [r(3), r(0), r(0)];
    let disk = |y: i128, rho: i128| {
        side(
            &wall,
            &MeetPoint::Narrow([r(0), r(y), r(0)]),
            &b(rho * rho), // the doors take the radii as squares
            &at,
            &axis,
            &b(25),
        )
    };
    // Ten away with a radius of two: clear, and on the `−(n × m)` side.
    assert_eq!(disk(10, 2), StripSide::Minus);
    assert_eq!(disk(-10, 2), StripSide::Plus);
    // Centred on a ruling: it spans that boundary.
    assert_eq!(disk(4, 2), StripSide::Crosses);
    assert_eq!(disk(-4, 2), StripSide::Crosses);
    // Reaching a ruling from outside, and stopping just short of one from inside.
    assert_eq!(disk(7, 4), StripSide::Crosses);
    assert_eq!(disk(0, 3), StripSide::Inside);
    // Touching is not spanning, on either side of the boundary.
    assert_eq!(disk(6, 2), StripSide::Inside);
    // ★ A point has no width, so it can never span: on the ruling it is `Inside`, as it was
    // before this cell and as `cylinder_strip_side` still answers.
    assert_eq!(disk(4, 0), StripSide::Inside);
    assert_eq!(
        disk(4, 0),
        nacre_exact::cylinder_strip_side(
            &wall,
            &MeetPoint::Narrow([r(0), r(4), r(0)]),
            &at,
            &axis,
            &b(25)
        )
    );
}

/// ★★★★★ **Why the net's population is empty today, frozen by name.**
///
/// A circle and a ruling of one class meet only if two lateral faces share a point, and the
/// cylinder-pair gate refuses that before any arrangement runs. This is that gate seen from
/// outside: two cylinders on crossing axes, one per operand, whose surfaces intersect. The
/// refusal is `CylinderPairContact` — *not* the arrangement's name — which is the whole content
/// of «the guard this cell deleted was standing where nothing could arrive».
#[test]
fn two_crossing_cylinders_are_refused_by_the_pair_gate() {
    let mut m = Model::new();
    let upright = {
        let profile = circle_profile([0.0, 0.0], 5.0);
        let frame = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 20.0,
            },
        )
        .expect("the upright extrudes") else {
            unreachable!()
        };
        solid
    };
    let across = turned_cylinder(&mut m, 5.0, 20.0, [-10.0, 0.0, 10.0]);
    m.rebuild_adjacency();
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        let out = boolean(&mut m, kind, upright, across);
        assert!(
            matches!(
                out,
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderPairContact,
                    ..
                })
            ),
            "{kind:?}: {out:?}"
        );
    }
}

/// ★★★★ **The one shape the walk has never seen, and why no fixture stands for it.**
///
/// A class carrying a circle **inside** a ruling cylinder's strip would put a disk afloat in a
/// ruling-bounded region: the arrangement's four edge kinds allow it, and nothing downstream has
/// ever been handed it. The audit counts that population and reports **zero** — so this asks
/// whether one could be *built*: the small cylinder crosses the big one's axis while lying wholly
/// **inside** it, which is exactly the placement whose cap plane would carry its own circle and
/// the big cylinder's rulings, the circle well within them.
///
/// ☑ Measured: it does not reach the arrangement at all — the operand gate speaks first, for both
/// kinds. So «zero inside» is not an instrument that failed to look; it is a shape the gate does
/// not admit yet. One placement is not a population claim, which is why this freezes the placement
/// by name rather than asserting anything wider.
///
/// ★ With the caps readable the plane-level questions
/// are answered instead of declined, and the pair rule gets to speak: these two cylinders really
/// do share points, so `CylinderPairContact` is the true name and `CylinderGateUndecided` would be
/// the road running out before reaching it.
#[test]
fn a_cylinder_nested_across_another_waits_at_the_gate() {
    let mut m = Model::new();
    let big = {
        let profile = circle_profile([0.0, 0.0], 20.0);
        let frame = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 40.0,
            },
        )
        .expect("the big one extrudes") else {
            unreachable!()
        };
        let z = |v: i128| nacre_exact::Rat::from_int(v);
        transform(
            &mut m,
            solid,
            &nacre_exact::Isometry::translation([z(0), z(0), z(-20)]),
        )
        .expect("centred on the origin")
    };
    let small = turned_cylinder(&mut m, 5.0, 10.0, [-5.0, 0.0, 0.0]);
    m.rebuild_adjacency();
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        let out = boolean(&mut m, kind, big, small);
        assert!(
            matches!(
                out,
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderPairContact,
                    ..
                })
            ),
            "{kind:?}: {out:?}"
        );
    }
}

/// ★★★★ **The shape that makes a disk *span* a tangent line, and why no verdict sees it
/// yet.**
///
/// A boss tangent to the box's wall `x = 0`, and a bar ending exactly on that wall whose cap is a
/// disk **containing** the tangent line. The gate writes the tangency row (measured:
/// `runs_through` comes back **true** — that wall face is one disk, and the line runs through it
/// along the disk's open chord, which only [`nacre_exact::StripSide::Crosses`] opens), and then the
/// cylinder-pair rule speaks first and the operation stops.
///
/// ☑ **The two cannot be separated today.** For the row to be written at all the disk must overlap
/// the boss's axis span; overlapping it puts the bar's own lateral inside the boss's reach, which
/// is exactly what `lateral_faces_clear` denies. So the spanning fold is exercised but no
/// *outcome* depends on it — the same shape as the «circle inside a strip». This freezes the
/// geometry so the work that opens `CylinderPairContact` finds it waiting.
///
/// ★ The pair rule can clear by a third separating direction and a reader stands in front of it,
/// so «cannot be separated today» is a statement about *this* geometry
/// and not about the rule in general. The freeze is what says so if that changes.
#[test]
fn a_disk_spanning_a_tangent_line_is_seen_before_the_pair_rule_speaks() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([40.0, 20.0, 10.0]),
    );
    // Tangent to `x = 0`: axis at x = 3, radius 3.
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([3.0, 10.0, 10.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        10.0,
    )
    .solid;
    m.rebuild_adjacency();
    let a = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss fuses")[0];
    // The bar's `+x` cap lands on `x = 0`, centred on the tangent line and overlapping the boss's
    // span, so its disk spans that line rather than clearing it.
    let bar = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([-20.0, 10.0, 9.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, -1.0]),
        3.0,
        20.0,
    )
    .solid;
    m.rebuild_adjacency();
    for kind in [BoolKind::Fuse, BoolKind::Cut] {
        let out = crate::boolean(&mut m, kind, a, bar);
        assert!(
            matches!(
                out,
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderPairContact,
                    ..
                })
            ),
            "{kind:?}: {out:?}"
        );
    }
}

/// ★ **The push funnel realizes; the fallback stands only where the realization declines** (cell
/// 52). Pushing a realizable definition again with a fallback that is wrong by three ulp gives a
/// `Bounded` cache holding the realization, not the fallback; a definition the road declines by
/// name — a pierce corner on a cylinder under an irrational turn, whose meet no road solves off
/// the world statement — keeps the fallback it was handed. The same turned cylinder's seam
/// vertices stand beside them, realized (met in the world).
#[test]
fn the_push_funnel_realizes_and_keeps_the_fallback_only_on_refusal() {
    let (mut m, l) = l_prism();
    let fh = m.shell(m.solid(l).outer).faces[0];
    let vh = m.edge(m.face(fh).outer.half_edges[0].edge).vertices[0];
    let def = *m.vertex(vh);
    let p = m.vertex_point(vh).as_array();
    let wrong = Point3::from_array([f64::from_bits(p[0].to_bits() + 3), p[1], p[2]]);
    let h = crate::realize::push_vertex_realized(
        &mut m,
        def,
        nacre_topo::PointCache::Unrealized { coord: wrong },
        crate::realize::ChainLink::Fresh,
    );
    assert!(matches!(
        m.vertex_cache(h),
        nacre_topo::PointCache::Bounded { .. }
    ));
    assert_eq!(
        m.vertex_point(h).as_array(),
        p,
        "the realization, not the fallback"
    );
    assert_ne!(m.vertex_point(h), wrong);

    let mut m = Model::new();
    let c = crate::fixtures::windowed_boss(&mut m, 0.6, -1.0);
    let turned =
        transform(&mut m, c, &rot_iso(Axis::Y, 37)).expect("an irrational turn records a motion");
    m.rebuild_adjacency();
    let (mut pierces, mut seams) = (0, 0);
    let vertices: std::collections::BTreeSet<_> = m
        .shell(m.solid(turned).outer)
        .faces
        .iter()
        .flat_map(|&fh| {
            let f = m.face(fh);
            std::iter::once(&f.outer)
                .chain(f.inner.iter())
                .flat_map(|lp| lp.half_edges.iter())
                .flat_map(|he| m.edge(he.edge).vertices)
                .collect::<Vec<_>>()
        })
        .collect();
    for vh in vertices {
        match *m.vertex(vh) {
            Vertex::Pierce { .. } => {
                pierces += 1;
                assert!(
                    matches!(
                        m.vertex_cache(vh),
                        nacre_topo::PointCache::Unrealized { coord: _ }
                    ),
                    "{:?}",
                    m.vertex_cache(vh)
                );
                assert_eq!(
                    crate::realize_vertex(&m, vh, crate::Precision::NearestF64).err(),
                    Some(crate::RealizeError::NoCurvedPoint)
                );
            }
            Vertex::OnSeam(_) => {
                seams += 1;
                assert!(
                    matches!(m.vertex_cache(vh), nacre_topo::PointCache::Bounded { .. }),
                    "{:?}",
                    m.vertex_cache(vh)
                );
            }
            Vertex::ThreePlane(_) => {}
        }
    }
    assert!(pierces > 0 && seams > 0, "{pierces} pierces, {seams} seams");
}

/// The slanted wall `z = 2 + x/6` over `x ∈ [−3, 3]`, `y ∈ [0, 4]`, everything above it to `z = 6`:
/// a trapezoid on the world `ZX` plane (`(u, v) = (z, x)`) extruded along `+y`. It cuts a lateral
/// on the `z` axis obliquely at about `z = 2`.
fn slanted_wedge(m: &mut Model) -> Handle<Solid> {
    let pts = [(1.5, -3.0), (2.5, 3.0), (6.0, 3.0), (6.0, -3.0)]
        .map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = crate::Ring2d::polygon_decimal(pts.to_vec()).expect("a trapezoid");
    let profile = crate::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(m, Axis::Y);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 4.0,
        },
    )
    .expect("the wedge extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// ★★★★ **A lateral bounded by two cut rims states its axis span from both.**
///
/// [`crate::fixtures::cut_at_both_rims`]'s lateral wraps the axis between two cut rims and has no
/// seam edge, so its upper rim is an **inner** loop. Its span is `[0, 4]`; read off the outer loop
/// alone it is `[0, 0.5]` — the lower rim's two stations — and the gates then prove a cylinder
/// crossing the lateral at `z = 2`, and a slanted wall at that height, clear of a face they cut:
/// the arrangement refuses them by the wrong names (`LabelConflict`, `CylinderStagesDisagree`),
/// and a debug build stops at `cyl_trace`'s span assertion. The crossing at `z = 0.25`, inside the
/// lower rim's stations, is the control the outer loop alone answers too.
#[test]
fn a_lateral_cut_at_both_rims_reads_its_span_from_both() {
    let refused_by = |out: Result<Vec<Handle<Solid>>, BoolError>, want: RejectReason| matches!(out, Err(BoolError::Rejected { reason, .. }) if reason == want);
    for built in [BoolKind::Fuse, BoolKind::Cut] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            for z in [2.0, 0.25] {
                let mut m = Model::new();
                let band = crate::fixtures::cut_at_both_rims(&mut m, built, -2.0);
                let across = turned_cylinder(&mut m, 0.2, 6.0, [-3.0, 0.0, z]);
                m.rebuild_adjacency();
                let out = boolean(&mut m, kind, band, across);
                assert!(
                    refused_by(out.clone(), RejectReason::CylinderPairContact),
                    "built by {built:?}, {kind:?} with a cylinder across at z = {z}: {out:?}"
                );
            }
            let mut m = Model::new();
            let band = crate::fixtures::cut_at_both_rims(&mut m, built, -2.0);
            let wall = slanted_wedge(&mut m);
            let out = boolean(&mut m, kind, band, wall);
            assert!(
                refused_by(out.clone(), RejectReason::ObliqueCylinderCut),
                "built by {built:?}, {kind:?} with a slanted wall: {out:?}"
            );
        }
    }
}

/// ★★★ **A lateral bounded by two cut rims is an operand like any other.** Fed back to a boolean
/// with a box slicing it through the middle (`x ∈ [−0.5, 3]`, `|y| ≤ 2`, `z ∈ [1.5, 2.5]`), every
/// kind builds a clean model of the analytic volume — the disk part with `x > −0.5`, one tall,
/// is `π − s` for `s` the segment past `|x| = 0.5`.
#[test]
fn a_lateral_cut_at_both_rims_is_an_operand_like_any_other() {
    let s = 0.5f64.acos() - 0.5 * 0.75f64.sqrt();
    let pi = std::f64::consts::PI;
    let slice = pi - s;
    for (built, body) in [
        (BoolKind::Fuse, 4.0 * pi + 30.0 - s),
        (BoolKind::Cut, 4.0 * pi - s),
    ] {
        for (kind, want) in [
            (BoolKind::Fuse, body + 14.0 - slice),
            (BoolKind::Cut, body - slice),
            (BoolKind::Common, slice),
        ] {
            let mut m = Model::new();
            let band = crate::fixtures::cut_at_both_rims(&mut m, built, -2.0);
            let tool = crate::fixtures::cuboid(
                &mut m,
                nacre_math::Point3::from_array([-0.5, -2.0, 1.5]),
                nacre_math::Point3::from_array([3.0, 2.0, 2.5]),
            );
            m.rebuild_adjacency();
            let out = boolean(&mut m, kind, band, tool)
                .unwrap_or_else(|e| panic!("built by {built:?}, {kind:?}: {e:?}"));
            m.rebuild_adjacency();
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{built:?} {kind:?}"
            );
            let v: f64 = out
                .iter()
                .map(|&b| nacre_props::mass_props(&m, b).unwrap().volume)
                .sum();
            assert!(
                (v - want).abs() < 1e-6,
                "built by {built:?}, {kind:?}: {v} vs {want}"
            );
        }
    }
}
