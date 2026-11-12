//! The cylinder gates in scalar form: ruling reach, strip side, the cylinder-pair gate — and the
//! push funnel that realizes a vertex or keeps its fallback.

use super::*;

/// **The negative control — the net can see a crossing.**
///
/// The audit reports zero circle–ruling crossings over the corpus, and that number means
/// «none there» only if the instrument can say «there». This drives the rational net directly,
/// where the geometry is written by hand rather than found: a class `x = 0`, a cylinder about
/// `z` at `x = 3` with radius `5` — so its two rulings stand at `y = ±4` — and a circle centred
/// on that plane whose radius is walked across each of them.
///
/// **The extent form and the disk form are one door.**
///
/// [`nacre_exact::cylinder_ruling_reached_extent`] answers for a piece stated by its two ends;
/// a disk is that piece with both ends the same. Handing it the **same** point and margin twice
/// must reproduce the disk answer exactly — otherwise the general form has drifted from the case
/// it generalizes, and every arc reading rests on it. Swept over centres, radii and both sides so
/// the agreement is a population, not a sample.
///
/// And the last block is what the generalization is *for*: two **different** ends read a reach
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
                        // **Crossing implies reaching, everywhere.** The open
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
    // A reach from `y = 4` to `y = 6` with no margin: it holds the `−` ruling at `y = 4`
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
    // **And the same reach does not *cross* it.** The arrangement's net asks about a
    // crossing, because a node no road mints is one; an edge tangent to another divides nothing,
    // which is the same sentence as for a tangency at a vertex.
    assert!(!reach(-1, false), "ending on the ruling is not crossing it");
    assert!(!reach(1, false));
    // And for a **disk**, where the two readings are the only thing separating a tangency from
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
    // Sides here are the **scalar** family's: the sign of `(p − o)·(n × m)`, which for this
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
    // The precondition (a class parallel to the axis) is the caller's to check — the arrangement
    // wrapper asks it and answers `None`; here it is a `debug_assert`, so this test does not drive
    // an off-axis class through the door.
}

/// **A disk clears the strip, spans it, or floats inside it, and the three are
/// told apart.** The gate's disk arm rests on this, and a predicate that answered one of them
/// always would still pass every model in the corpus.
///
/// The class `x = 0`; a cylinder about `z` at `x = 3` with radius 5, so the strip runs `|U| ≤ 4`
/// in the `n × m` coordinate — which for this wall and axis points along `−y`.
///
/// The last case is the one the vocabulary exists for: a **point** on a ruling spans nothing, so
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
    // A point has no width, so it can never span: on the ruling it is `Inside`, as
    // `cylinder_strip_side` answers.
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

/// **Why the net's population is empty today, frozen by name.**
///
/// A circle and a ruling of one class meet only if two lateral faces share a point, and the
/// cylinder-pair gate refuses that before any arrangement runs. This is that gate seen from
/// outside: two cylinders on crossing axes, one per operand, whose surfaces intersect. The refusal
/// is `CylinderPairContact` — *not* the arrangement's name — so a guard for this case inside the
/// arrangement would stand where nothing can arrive.
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

/// **Every kind, both orders, against the analytic volume** — the cylinder pairs whose faces do
/// not meet. `want` is `(bodies, volume)` for `[Fuse, Cut, Common]` with the operands as built,
/// then swapped ([`every_kind`] reads the rest).
fn builds_every_kind(
    build: impl Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>),
    want: [[(usize, f64); 3]; 2],
) {
    every_kind(build, want.map(|row| row.map(|(n, v)| Want::Built(n, v))));
}

fn axis_cylinder(m: &mut Model, base: [f64; 3], axis: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
    crate::fixtures::cylinder(m, Point3::from_array(base), Vector3::from_array(axis), r, h).solid
}

/// **A cylinder nested crosswise inside another clears the pair rule** — by containment, the
/// witness no separating direction can stand in for.
///
/// The small cylinder (`r = 5`, along `x` over `[−5, 5]`) lies wholly inside the big one (`r = 20`,
/// along `z`): its lateral stays within `√50` of the big axis. Its axes cross the big one's, and
/// along every rational direction the two reaches overlap, so the separation witness has nothing
/// to say; the corners of the small lateral's box all lie strictly inside the big cylinder, and that
/// proves the two laterals apart. Downstream, the small cylinder's caps `x = ±5` carry its circle
/// and the big cylinder's rulings at once, the circle inside the strip between them — a class the
/// arrangement draws like any other.
#[test]
fn a_cylinder_nested_across_another_clears_the_pair_rule() {
    use std::f64::consts::PI;
    builds_every_kind(
        |m| {
            // Built as a prism and moved, so the small one's lateral rides a motion chain.
            let profile = circle_profile([0.0, 0.0], 20.0);
            let frame = SketchFrame::world(m, Axis::Z);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
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
            let big = transform(
                m,
                solid,
                &nacre_exact::Isometry::translation([z(0), z(0), z(-20)]),
            )
            .expect("centred on the origin");
            (big, turned_cylinder(m, 5.0, 10.0, [-5.0, 0.0, 0.0]))
        },
        [
            [(1, 16000.0 * PI), (1, 15750.0 * PI), (1, 250.0 * PI)],
            [(1, 16000.0 * PI), (0, 0.0), (1, 250.0 * PI)],
        ],
    );
}

/// **A short crosswise cylinder beside another builds, though its ruling's line runs through the
/// other's cap circle.** The tall one (`r = 1`, cap on `z = 1/2`) and the short one (`r = 1/5`,
/// along `x` over `[2, 3]`, its axis at `z = 0.6`) are apart along `x`. The cap plane `z = 1/2`
/// cuts the short cylinder in two rulings at `y = ±√0.03` — irrational, so their ends are pierce
/// corners in a radical — and each ruling's **line** crosses the tall cylinder's cap circle near
/// `x = ±1`; its **segment**, `x ∈ [2, 3]`, never comes near it. The class-edge net asks the
/// segment.
#[test]
fn a_short_cylinder_beside_another_builds_past_its_cap_circle() {
    use std::f64::consts::PI;
    let (tall, short) = (2.5 * PI, 0.04 * PI);
    builds_every_kind(
        |m| {
            (
                axis_cylinder(m, [0.0, 0.0, 0.5], [0.0, 0.0, 1.0], 1.0, 2.5),
                axis_cylinder(m, [2.0, 0.0, 0.6], [1.0, 0.0, 0.0], 0.2, 1.0),
            )
        },
        [
            [(2, tall + short), (1, tall), (0, 0.0)],
            [(2, tall + short), (1, short), (0, 0.0)],
        ],
    );
}

/// **A crosswise stud standing proud of a big cylinder's cap.** The stud (`r = 5`, along `x` over
/// `[−5, 5]`, axis at `z = 18`) lies inside the big cylinder's surface (`r = 20`) and pokes out
/// through its cap `z = 20`. The cap's class carries the big circle and the stud's rulings
/// (`y = ±√21`): the rulings' lines cross the circle near `x = ±19.5`, their segments lie inside
/// its disk. Common is the stud less the segment above the cap, `10·(25π − (25·acos(2/5) − 2√21))`.
#[test]
fn a_crosswise_stud_through_a_cylinders_cap_builds() {
    use std::f64::consts::PI;
    let (big, stud) = (16000.0 * PI, 250.0 * PI);
    let common = 10.0 * (25.0 * PI - (25.0 * 0.4f64.acos() - 2.0 * 21f64.sqrt()));
    builds_every_kind(
        |m| {
            (
                axis_cylinder(m, [0.0, 0.0, -20.0], [0.0, 0.0, 1.0], 20.0, 40.0),
                axis_cylinder(m, [-5.0, 0.0, 18.0], [1.0, 0.0, 0.0], 5.0, 10.0),
            )
        },
        [
            [(1, big + stud - common), (1, big - common), (1, common)],
            [(1, big + stud - common), (1, stud - common), (1, common)],
        ],
    );
}

/// **A pin floating crosswise inside a bore** — the containing lateral is a hole's (material
/// outside it): a `40 × 40 × 8` plate bored `r = 10` about `z`, a pin `r = 2` along `x` over
/// `[−3, 3]`. Apart, so the Fuse is two bodies.
#[test]
fn a_pin_floating_across_a_bore_builds() {
    use std::f64::consts::PI;
    let (bored, pin) = (12800.0 - 800.0 * PI, 24.0 * PI);
    builds_every_kind(
        |m| {
            let plate = crate::fixtures::cuboid(m, [-20.0, -20.0, -4.0], [20.0, 20.0, 4.0]);
            let bore = axis_cylinder(m, [0.0, 0.0, -10.0], [0.0, 0.0, 1.0], 10.0, 20.0);
            let bored = boolean(m, BoolKind::Cut, plate, bore).expect("the bore")[0];
            (
                bored,
                axis_cylinder(m, [-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 6.0),
            )
        },
        [
            [(2, bored + pin), (1, bored), (0, 0.0)],
            [(2, bored + pin), (1, pin), (0, 0.0)],
        ],
    );
}

/// **One class carrying two cylinders' rulings.** A half cylinder (`r = 20` about `z`, `y ≤ 0`)
/// and a crosswise one (`r = 5` along `x`) inside its surface and through its flat wall `y = 0`:
/// that wall holds the half cylinder's rulings at `x = ±20` and the crosswise one's at `z = ±5`,
/// whose lines cross and whose segments do not.
#[test]
fn a_cylinder_across_a_half_cylinders_flat_wall_builds() {
    use std::f64::consts::PI;
    let (half, inside) = (8000.0 * PI, 125.0 * PI);
    builds_every_kind(
        |m| {
            let big = axis_cylinder(m, [0.0, 0.0, -20.0], [0.0, 0.0, 1.0], 20.0, 40.0);
            let cut = crate::fixtures::cuboid(m, [-30.0, 0.0, -30.0], [30.0, 30.0, 30.0]);
            let half = boolean(m, BoolKind::Cut, big, cut).expect("the half")[0];
            (
                half,
                axis_cylinder(m, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5.0, 10.0),
            )
        },
        [
            [(1, half + inside), (1, half - inside), (1, inside)],
            [(1, half + inside), (1, inside), (1, inside)],
        ],
    );
}

/// **Containment is strict — a rim touching the other lateral from inside is still a contact.**
/// A crosswise cylinder `r = 3` over `x ∈ [−4, 4]` inside one of `r = 5` about `z`: its rim points
/// `(±4, ±3, 0)` lie on the outer lateral (`16 + 9 = 25`), and so do its box's corners. Read closed,
/// the containment witness would wave the touching pair through.
#[test]
fn a_rim_touching_the_other_lateral_from_inside_is_refused() {
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        let mut m = Model::new();
        let outer = axis_cylinder(&mut m, [0.0, 0.0, -5.0], [0.0, 0.0, 1.0], 5.0, 10.0);
        let inner = axis_cylinder(&mut m, [-4.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 8.0);
        let out = boolean(&mut m, kind, outer, inner);
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

/// **The shape that makes a disk *span* a tangent line, and why no verdict sees it
/// yet.**
///
/// A boss tangent to the box's wall `x = 0`, and a bar ending exactly on that wall whose cap is a
/// disk **containing** the tangent line. The gate writes the tangency row (measured: the line runs
/// through that wall face — one disk — along the disk's open chord, which only
/// [`nacre_exact::StripSide::Crosses`] opens), and then the
/// cylinder-pair rule speaks first and the operation stops.
///
/// **The two cannot be separated today.** For the row to be written at all the disk must overlap
/// the boss's axis span; overlapping it puts the bar's rim on the tangent line — the point
/// `(0, 10, 12)` is on both laterals — so the two faces really touch, and `lateral_faces_clear`
/// rightly clears nothing. So the spanning fold is exercised but no *outcome* depends on it. This
/// freezes the geometry so the work that opens `CylinderPairContact` finds it waiting.
///
/// The pair rule can clear by a third separating direction and a reader stands in front of it,
/// so «cannot be separated today» is a statement about *this* geometry
/// and not about the rule in general. The freeze is what says so if that changes.
#[test]
fn a_disk_spanning_a_tangent_line_is_seen_before_the_pair_rule_speaks() {
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(&mut m, [0.0, 0.0, 0.0], [40.0, 20.0, 10.0]);
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

/// **The push funnel realizes; the fallback stands only where the realization declines**. Pushing a
/// realizable definition again with a fallback that is wrong by three ulp gives a
/// `Realized` cache holding the realization, not the fallback; a definition the road declines by
/// name — a pierce corner asked for a third crossing, which no line makes with a cylinder — keeps
/// the fallback it was handed. Beside it a cylinder under an irrational turn has its pierce
/// corners and its seam vertices realized (on their rulings and rims, met in the world).
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
        nacre_topo::PointCache::Realized { .. }
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
    let (mut pierces, mut seams) = (Vec::new(), 0);
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
            Vertex::Pierce {
                planes, cylinder, ..
            } => pierces.push((planes, cylinder)),
            Vertex::OnSeam(_) => seams += 1,
            Vertex::ThreePlane(_) | Vertex::CylinderPair { .. } => continue,
        }
        assert!(
            matches!(m.vertex_cache(vh), nacre_topo::PointCache::Realized { .. }),
            "{:?}: {:?}",
            m.vertex(vh),
            m.vertex_cache(vh)
        );
    }
    assert!(
        !pierces.is_empty() && seams > 0,
        "{pierces:?}, {seams} seams"
    );
    let (planes, cylinder) = pierces[0];
    let fallback = Point3::from_array([1.0, 2.0, 3.0]);
    let third = crate::realize::push_vertex_realized(
        &mut m,
        Vertex::Pierce {
            planes,
            cylinder,
            rank: 2,
        },
        nacre_topo::PointCache::Unrealized { coord: fallback },
        crate::realize::ChainLink::Fresh,
    );
    assert_eq!(
        *m.vertex_cache(third),
        nacre_topo::PointCache::Unrealized { coord: fallback },
        "the declined definition keeps the fallback it was handed"
    );
    assert_eq!(
        crate::realize_vertex(&m, third, crate::Precision::NearestF64).err(),
        Some(crate::RealizeError::NoCurvedPoint)
    );
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
    solid
}

/// **A lateral bounded by two cut rims states its axis span from both.**
///
/// [`crate::fixtures::cut_at_both_rims`]'s lateral wraps the axis between two cut rims and has no
/// seam edge, so its upper rim is an **inner** loop. Its span is `[0, 4]`; read off the outer loop
/// alone it is `[0, 0.5]` — the lower rim's two stations — and the gates then prove a cylinder
/// crossing the lateral at `z = 2`, and a slanted wall at that height, clear of a face they cut:
/// the arrangement refuses them by the wrong names (`LabelConflict`, `CylinderStagesDisagree`),
/// and a debug build stops at `cyl_trace`'s span assertion. The crossing at `z = 0.25`, inside the
/// lower rim's stations, is the control the outer loop alone answers too.
///
/// The slanted wall (`slanted_wedge`, `z ≥ 2 + x/6` over `y ∈ [0, 4]`) sections the lateral
/// between its two cut rims and crosses neither, so the whole ellipse lies in the face: every kind
/// builds at its derived volume. With `s` the disk's segment past `x = 0.5`, the wedge (volume 96)
/// takes from the cylinder the half `y ≥ 0` above the plane, `∬(2 − x/6) = π`, and from the upper
/// bite (fused) its part `y ≥ 0`, `7.5`, less what the two share, `s/4`.
#[test]
fn a_lateral_cut_at_both_rims_reads_its_span_from_both() {
    let refused_by = |out: Result<Vec<Handle<Solid>>, BoolError>, want: RejectReason| matches!(out, Err(BoolError::Rejected { reason, .. }) if reason == want);
    for built in [BoolKind::Fuse, BoolKind::Cut] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            for z in [2.0, 0.25] {
                let mut m = Model::new();
                let band = crate::fixtures::cut_at_both_rims(&mut m, built, -2.0);
                let across = turned_cylinder(&mut m, 0.2, 6.0, [-3.0, 0.0, z]);
                let out = boolean(&mut m, kind, band, across);
                assert!(
                    refused_by(out.clone(), RejectReason::CylinderPairContact),
                    "built by {built:?}, {kind:?} with a cylinder across at z = {z}: {out:?}"
                );
            }
            let s = 0.5f64.acos() - 0.5 * 0.75f64.sqrt();
            let pi = std::f64::consts::PI;
            let (band_volume, common) = match built {
                BoolKind::Fuse => (4.0 * pi + 30.0 - s, pi + 7.5 - s / 4.0),
                _ => (4.0 * pi - s, pi - s / 4.0),
            };
            let want = match kind {
                BoolKind::Fuse => band_volume + 96.0 - common,
                BoolKind::Cut => band_volume - common,
                BoolKind::Common => common,
            };
            let mut m = Model::new();
            let band = crate::fixtures::cut_at_both_rims(&mut m, built, -2.0);
            let wall = slanted_wedge(&mut m);
            let out = boolean(&mut m, kind, band, wall).unwrap_or_else(|e| {
                panic!("built by {built:?}, {kind:?} with a slanted wall: {e:?}")
            });
            let got = built_volume(
                &m,
                &out,
                1,
                &format!("{built:?} {kind:?} with a slanted wall"),
            );
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "built by {built:?}, {kind:?} with a slanted wall: {got} vs {want}"
            );
        }
    }
}

/// **A lateral bounded by two cut rims is an operand like any other.** Fed back to a boolean
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
            let tool = crate::fixtures::cuboid(&mut m, [-0.5, -2.0, 1.5], [3.0, 2.0, 2.5]);
            let out = boolean(&mut m, kind, band, tool)
                .unwrap_or_else(|e| panic!("built by {built:?}, {kind:?}: {e:?}"));
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

/// **The record of a tangent wall's line carries where that wall's faces lie on each side of it**
/// (`SharedRuling::tangent_halves`, in the lateral's `+τ`, `−τ` order). `A` is the notched block
/// (its notch's lateral through the ruling `(1, 0)`, `τ = +y` there), `B` a prism with its wall
/// `x = 1` tangent to the notch along that line and a face on the plane `y = 0`, which holds the
/// line and cuts the notch — the pair the record is written for. Over the lateral's span `z ∈ [0, 2]`:
///
/// - a box above the line (`y ∈ [0, 2]`) puts `x = 1` on the `+τ` side only — a face ending on the
///   line, the whole span; below it (`y ∈ [−2, 0]`), on the `−τ` side only — the sign of `±w → ±τ`;
/// - the same box standing in `z ∈ [1/2, 3/2]`, only there;
/// - a prism whose wall `x = 1` runs across the line (`y ∈ [−1, 1]`), with its face on `y = 0` off to
///   the side, on both sides — a contact, which divides nothing.
#[test]
fn a_tangent_walls_record_says_which_side_of_the_line_its_faces_lie_on() {
    let halves_of = |quad: &[[f64; 2]], lift: f64, height: f64| {
        let (m, a, b) = super::fixtures::a_solid_and_a_placed_prism(
            |m| crate::SketchFrame::world(m, nacre_exact::Axis::Z),
            super::fixtures::notched_block(),
            quad,
            lift,
            height,
        );
        let setup = plane_index_setup(&m, a, b).expect("the gate passes it");
        // The line through `(1, 0)` — the walls `y = ±1` of the prism across it are tangent to the
        // notch's end rulings too, with no face there.
        let at = [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)];
        let records: Vec<_> = setup.cyls[0]
            .shared
            .iter()
            .filter(|sr| sr.at == at)
            .filter_map(|sr| sr.tangent_halves.clone())
            .collect();
        assert_eq!(
            records.len(),
            1,
            "one tangent wall, one holder: {records:?}"
        );
        records[0]
            .iter()
            .map(|side| {
                side.iter()
                    .map(|p| p.map(|v| v.as_rat().expect("a rational end")))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a ratio");
    let whole = vec![[r(0, 1), r(2, 1)]];
    let above: &[[f64; 2]] = &[[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]];
    let below: &[[f64; 2]] = &[[1.0, -2.0], [3.0, -2.0], [3.0, 0.0], [1.0, 0.0]];
    let across: &[[f64; 2]] = &[
        [1.0, -1.0],
        [3.0, -1.0],
        [3.0, 0.0],
        [5.0, 0.0],
        [5.0, 1.0],
        [1.0, 1.0],
    ];
    assert_eq!(halves_of(above, -0.5, 3.0), vec![whole.clone(), vec![]]);
    assert_eq!(halves_of(below, -0.5, 3.0), vec![vec![], whole.clone()]);
    assert_eq!(
        halves_of(above, 0.5, 1.0),
        vec![vec![[r(1, 2), r(3, 2)]], vec![]]
    );
    assert_eq!(halves_of(across, -0.5, 3.0), vec![whole.clone(), whole]);
}

/// A rod — `r = 2` along `x` over `[−20, 20]`, its axis at height `14` — lying on [`rounded_plate`]'s
/// top face `z = 12`, tangent to it along the line `y = y0`.
fn rod_on_the_plate(m: &mut Model, fillets: usize, y0: f64) -> (Handle<Solid>, Handle<Solid>) {
    let plate = rounded_plate(m, fillets, 4);
    let rod = axis_cylinder(m, [-20.0, y0, 14.0], [1.0, 0.0, 0.0], 2.0, 40.0);
    (plate, rod)
}

/// **A rod lying on a plate whose corners are rounded.** The rod touches the plate's top face
/// along a line, so the gate writes a tangency row for that face, and the row reads where the line
/// runs through the face (`line_runs_through_face`) — a face whose outline carries the four fillet
/// arcs. Square-cornered (`rounded_plate(_, 0, _)`) or rounded, across the middle (`y = 0`, far
/// from every arc) or along the edge (`y = −22.5`, where the line runs on past the rod and meets
/// the outline only on two fillet arcs, at `x = ±(40 + √18.75)`), every kind builds: the fuse is
/// the two bodies that meet along the line alone, the cuts leave each operand whole, the common
/// part is empty.
#[test]
fn a_rod_lying_on_a_rounded_plate() {
    let pi = std::f64::consts::PI;
    let rod = pi * 4.0 * 40.0;
    let square = 90.0 * 50.0 * 12.0 - 4.0 * pi * 3.5 * 3.5 * 12.0;
    for (fillets, y0, plate) in [
        (0, 0.0, square),
        (4, 0.0, rounded_plate_volume()),
        (4, -22.5, rounded_plate_volume()),
    ] {
        every_kind(
            |m| rod_on_the_plate(m, fillets, y0),
            [
                [
                    Want::Built(2, plate + rod),
                    Want::Built(1, plate),
                    Want::Built(0, 0.0),
                ],
                [
                    Want::Built(2, plate + rod),
                    Want::Built(1, rod),
                    Want::Built(0, 0.0),
                ],
            ],
        );
    }
}

/// **A slab across a cylinder's tangent line, where only its fillet arcs meet that line.** `A` is
/// the unit cylinder along `z` over `[0, 2]`; `B` the slab `x ∈ [0, 1]` over the rectangle
/// `y ∈ [−4, 3/10]`, `z ∈ [−3, 5]`, square or with every corner rounded `r = 1/2`. Its face `x = 1` is
/// tangent to `A` along `(1, 0)` and lies on both sides of that line, as `A` does: both solids run
/// through it. The line meets the face's outline only on the two fillet arcs at `y1 = 3/10`, far
/// past `A`'s span — none of the straight edges reaches it — so where the line runs through the
/// face is read from the arcs alone (`line_runs_through_face`). `B − A` keeps two cusps of the slab
/// meeting along the line, a pinch refused as `SelfTouchingResult`; read as meeting nothing, the line
/// would write no row and the pinch would build. Every other kind builds: with `C = A ∩ B`, the part
/// of the disk at `x ≥ 0`, `y ≤ 3/10`, two tall — `π/2 + (3/10)√(91/100) + asin(3/10)` — `A ∪ B` is
/// `2π + |B| − C` and `A − B` is `2π − C`.
#[test]
fn a_rounded_slab_across_a_cylinders_tangent_line() {
    let pi = std::f64::consts::PI;
    let a = 2.0 * pi;
    let common = pi / 2.0 + 0.3 * 0.91f64.sqrt() + 0.3f64.asin();
    let rounded_rect = |r: f64| {
        let (y0, y1, z0, z1) = (-4.0, 0.3, -3.0, 5.0);
        if r == 0.0 {
            return crate::Profile2d::polygon(vec![p2(y0, z0), p2(y1, z0), p2(y1, z1), p2(y0, z1)])
                .expect("a rectangle");
        }
        stated(vec![
            line(p2(y0 + r, z0), p2(y1 - r, z0)),
            arc_turns(p2(y1 - r, z0 + r), p2(y1 - r, z0), 1),
            line(p2(y1, z0 + r), p2(y1, z1 - r)),
            arc_turns(p2(y1 - r, z1 - r), p2(y1, z1 - r), 1),
            line(p2(y1 - r, z1), p2(y0 + r, z1)),
            arc_turns(p2(y0 + r, z1 - r), p2(y0 + r, z1), 1),
            line(p2(y0, z1 - r), p2(y0, z0 + r)),
            arc_turns(p2(y0 + r, z0 + r), p2(y0, z0 + r), 1),
        ])
        .expect("a rounded rectangle")
        .remove(0)
    };
    for (r, b) in [(0.0, 34.4), (0.5, 33.4 + pi / 4.0)] {
        let build = |m: &mut Model| {
            let cyl = axis_cylinder(m, [0.0; 3], [0.0, 0.0, 1.0], 1.0, 2.0);
            let frame = match apply(
                m,
                &Operation::DatumPlane {
                    def: DatumDef::Stated(SketchPlane::from_axes(
                        Point3::from_array([0.0; 3]),
                        Vector3::from_array([0.0, 1.0, 0.0]),
                        Vector3::from_array([0.0, 0.0, 1.0]),
                    )),
                },
            ) {
                Ok(OpOutput::DatumPlane { frame, .. }) => frame,
                other => panic!("the frame on x = 0: {other:?}"),
            };
            let OpOutput::Extrude { solid: slab, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile: rounded_rect(r),
                    dist: 1.0,
                },
            )
            .expect("the slab extrudes") else {
                unreachable!()
            };
            (cyl, slab)
        };
        every_kind(
            build,
            [
                [
                    Want::Built(1, a + b - common),
                    Want::Built(1, a - common),
                    Want::Built(1, common),
                ],
                [
                    Want::Built(1, a + b - common),
                    Want::Refused(RejectReason::SelfTouchingResult),
                    Want::Built(1, common),
                ],
            ],
        );
    }
}

/// A D-shaped rod — the unit disk's three quarters from `(1, 0)` round to `(0, −1)`, closed by the
/// square corner `(1, −1)` — along `z` over `[0, 2]`, its flat `x = 1` tangent to its own lateral
/// along `(1, 0)`; and a slab `x ∈ [1, 2]` against that flat, bitten at its far edge `y = 4` by a
/// bore along `x` about `(y, z) = (4.5, 1)`, `r = 1`, so the slab's face `x = 1` carries an arc
/// whose ends `z = 1 ± √3/2` are irrational.
pub(crate) fn d_rod_on_a_bitten_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let rod = d_rod(m);
    let slab = crate::fixtures::cuboid(m, [1.0, -4.0, -3.0], [2.0, 4.0, 5.0]);
    let bore = axis_cylinder(m, [0.5, 4.5, 1.0], [1.0, 0.0, 0.0], 1.0, 2.0);
    let slab = boolean(m, BoolKind::Cut, slab, bore).expect("the bite")[0];
    (rod, slab)
}

/// [`d_rod_on_a_bitten_slab`]'s rod against the slab grooved across `x = 1` over `z ∈ [2.4, 2.6]`,
/// which splits that face in two: the lower one, plain, runs through the rod's line over the
/// lateral's `z ∈ [0, 2]`; the upper one carries the bite (about `z = 3.8`, ends
/// `z = 3.8 ± √3/2`) and stands above the lateral.
pub(crate) fn d_rod_on_a_split_bitten_slab(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let rod = d_rod(m);
    let slab = crate::fixtures::cuboid(m, [1.0, -4.0, -3.0], [2.0, 4.0, 5.0]);
    let groove = crate::fixtures::cuboid(m, [0.5, -5.0, 2.4], [1.5, 5.0, 2.6]);
    let slab = boolean(m, BoolKind::Cut, slab, groove).expect("the groove")[0];
    let bore = axis_cylinder(m, [0.5, 4.5, 3.8], [1.0, 0.0, 0.0], 1.0, 2.0);
    let slab = boolean(m, BoolKind::Cut, slab, bore).expect("the bite")[0];
    (rod, slab)
}

fn d_rod(m: &mut Model) -> Handle<Solid> {
    let profile = stated(vec![
        arc_turns(p2(0.0, 0.0), p2(1.0, 0.0), 3),
        line(p2(0.0, -1.0), p2(1.0, -1.0)),
        line(p2(1.0, -1.0), p2(1.0, 0.0)),
    ])
    .expect("the D")
    .remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid: rod, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    )
    .expect("the D rod extrudes") else {
        unreachable!()
    };
    rod
}

/// **A wall face the gate cannot read carries its own refusal on the tangency row.** The bitten
/// slab's face `x = 1` touches the D rod's lateral along `(1, 0)`, and its arc's ends are
/// irrational, which the footprint reader cannot spell (`face_clears_footprint`) — the rod's
/// solid's own flat on that plane is the class's first face and does not clear the footprint, so
/// the class is read row by row (`tangency_rows`), and the slab's face writes a row carrying that
/// reader's refusal, `CurvedOperandBoundary`, which the gate raises. With the slab given first, the
/// class's first face is the bitten one and the class-wide read raises the same name before any
/// row is written. Both raise the very value the footprint reader made, so this answer cannot say
/// which road refused; `an_unread_wall_face_rides_its_tangency_row` asks each road directly.
#[test]
fn an_unread_wall_face_refuses_through_its_tangency_row() {
    let refused = [Want::Refused(RejectReason::CurvedOperandBoundary); 3];
    every_kind(d_rod_on_a_bitten_slab, [refused, refused]);
}
