//! Profiles with holes — extruding a region that is not simply connected.
//!
//! The kernel could already *produce* a donut, by cutting a bar out of a box; what it could not do
//! was **construct** one. That distinction is the point: the boolean route makes every corner a
//! `Discovered` intersection with a measured tolerance, while extruding the profile directly makes
//! them `Constructed` — exact by construction (overview 절대원칙 4). So the strongest check here
//! is not a number typed by hand but the two producers agreeing.

use crate::common::*;
use crate::stated::*;
use nacre_math::{Point2, Point3};
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply, from_rings};
use nacre_scalar::Axis;
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, Solid};

/// `[0,4]²` with a `[1,3]²` hole, extruded 1 high: volume 16 − 4 = 12.
fn donut_profile() -> Profile2d {
    let sq = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    Profile2d::with_holes(sq(0.0, 4.0), vec![sq(1.0, 3.0)]).unwrap()
}

fn extrude(m: &mut Model, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile,
            dist,
        },
    )
    .unwrap() else {
        unreachable!("Extrude yields an Extrude output")
    };
    m.rebuild_adjacency();
    solid
}

/// **The check that cannot be fooled by my arithmetic.** The same donut, built two ways: swept
/// from a profile with a hole, and cut out of a solid box. Every measurable must agree — a hole
/// that failed to open shows up in the volume, walls that are missing or inside-out show up in
/// the area, and a mis-wound cap shows up in the face count.
#[test]
fn a_swept_hole_matches_the_same_shape_cut_out() {
    let mut m = Model::new();
    let swept = extrude(&mut m, donut_profile(), 1.0);
    let swept_props = nacre_props::mass_props(&m, swept).unwrap();
    let swept_faces = m.shells.get(m.solids.get(swept).outer).faces.len();

    let mut m2 = Model::new();
    let block = m2.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([4.0, 4.0, 1.0]),
    );
    let bar = m2.add_cuboid(
        Point3::from_array([1.0, 1.0, -1.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m2.rebuild_adjacency();
    let cut = boolean_one(&mut m2, BoolKind::Cut, block, bar).unwrap();
    m2.rebuild_adjacency();
    let cut_props = nacre_props::mass_props(&m2, cut).unwrap();
    let cut_faces = m2.shells.get(m2.solids.get(cut).outer).faces.len();

    assert!(
        (swept_props.volume - cut_props.volume).abs() < 1e-12,
        "volume {} vs {}",
        swept_props.volume,
        cut_props.volume
    );
    assert!(
        (swept_props.area - cut_props.area).abs() < 1e-12,
        "area {} vs {}",
        swept_props.area,
        cut_props.area
    );
    assert_eq!(swept_faces, cut_faces, "face count");
}

/// The construction is exact where the boolean's is not: sweeping a profile makes every vertex
/// `Constructed`, so no tolerance is recorded anywhere. This is the reason the kernel grows a
/// producer for a shape it could already cut.
#[test]
fn a_swept_hole_is_constructed_throughout() {
    let mut m = Model::new();
    let d = extrude(&mut m, donut_profile(), 1.0);
    for &fh in &m.shells.get(m.solids.get(d).outer).faces {
        let face = m.faces.get(fh).clone();
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for he in &lp.half_edges {
                for vh in m.edges.get(he.edge).vertices.iter() {
                    assert!(
                        m.vertex_tol(*vh).is_none(),
                        "a swept vertex carries no tolerance"
                    );
                    assert!(
                        matches!(m.vertex_cache(*vh), PointCache::Bounded { .. }),
                        "a swept vertex is realized from its definition (cell 52)"
                    );
                }
            }
        }
    }
}

/// Euler with the inner-loop term, on the first shape that actually needs it: a donut prism is
/// genus 1, and `χ = V − E + F − L_i = 16 − 24 + 10 − 2 = 0 = 2(S − G)` only if both cap holes
/// are counted. Drop the inner loops and the same solid reads as an impossible genus.
#[test]
fn a_donut_prism_is_a_valid_genus_one_solid() {
    let mut m = Model::new();
    let d = extrude(&mut m, donut_profile(), 1.0);

    let vs = nacre_validate::validate(&m);
    assert!(vs.is_empty(), "{vs:?}");

    let faces = &m.shells.get(m.solids.get(d).outer).faces;
    assert_eq!(faces.len(), 10, "4 outer walls + 4 hole walls + 2 caps");
    let inner_loops: usize = faces.iter().map(|&fh| m.faces.get(fh).inner.len()).sum();
    assert_eq!(inner_loops, 2, "one hole loop per cap");
    assert!((volume(&m, d) - 12.0).abs() < 1e-12, "{}", volume(&m, d));
}

/// A hole is only a hole if the mesh and the exchange format agree it is one. The tessellator
/// sweeps the cap with its hole as a sibling ring, and STEP writes it as an inner face bound.
///
/// ★ This used to say tessellation *"bridges"* cap holes. It does not, and has not since the
/// monotone sweep landed — bridging is exactly what that sweep exists to avoid.
#[test]
fn a_donut_prism_tessellates_and_exports() {
    let mut m = Model::new();
    extrude(&mut m, donut_profile(), 1.0);

    assert!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok(),
        "caps triangulate with holes"
    );
    let step = nacre_step::to_step(&m).unwrap();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
    assert!(step.contains("FACE_BOUND"), "the hole is an inner bound");
}

/// Several holes at once — the walls and cap loops are per-ring, so two is not a special case.
#[test]
fn a_profile_with_two_holes() {
    let sq = |ax: f64, ay: f64, bx: f64, by: f64| {
        vec![
            Point2::from_array([ax, ay]),
            Point2::from_array([bx, ay]),
            Point2::from_array([bx, by]),
            Point2::from_array([ax, by]),
        ]
    };
    let mut m = Model::new();
    let p = Profile2d::with_holes(
        sq(0.0, 0.0, 6.0, 3.0),
        vec![sq(1.0, 1.0, 2.0, 2.0), sq(4.0, 1.0, 5.0, 2.0)],
    )
    .unwrap();
    let d = extrude(&mut m, p, 1.0);

    assert!(nacre_validate::validate(&m).is_empty());
    assert!(
        (volume(&m, d) - (18.0 - 2.0)).abs() < 1e-12,
        "{}",
        volume(&m, d)
    );
    let inner_loops: usize = m
        .shells
        .get(m.solids.get(d).outer)
        .faces
        .iter()
        .map(|&fh| m.faces.get(fh).inner.len())
        .sum();
    assert_eq!(inner_loops, 4, "two holes on each of two caps");
}

/// A solid built with a hole is an ordinary operand — the boolean has always accepted faces with
/// inner loops, and this checks the constructed variety is no different.
#[test]
fn a_swept_hole_is_a_boolean_operand() {
    let mut m = Model::new();
    let d = extrude(&mut m, donut_profile(), 1.0);
    let knife = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.5]),
        Point3::from_array([5.0, 5.0, 2.0]),
    );
    m.rebuild_adjacency();

    let r = boolean_one(&mut m, BoolKind::Cut, d, knife).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    assert!((volume(&m, r) - 6.0).abs() < 1e-12, "{}", volume(&m, r));
}

/// The other sweep direction. A pocket sweeps *into* its face, so the outer ring is reversed on
/// the way in — the one path where "which way is counter-clockwise" differs from every extrude
/// above.
///
/// It was written expecting to be the only test that catches a mis-wound hole, and measurement
/// said otherwise: winding each hole opposite to the *normalized* outer ring and winding it
/// clockwise about the sweep are the same rule, because the outer ring is always normalized
/// counter-clockwise about the sweep. There is no second rule to get wrong. What remains real is
/// the hazard the design removed — deciding winding in two places (the frame's 2-D area *and* the
/// sweep's) — and this test is the end-to-end cover for the direction where those two disagree.
#[test]
fn a_pocket_with_a_hole_sweeps_the_other_way() {
    let mut m = Model::new();
    let block = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([10.0, 10.0, 4.0]),
    );
    m.rebuild_adjacency();
    let top = m
        .shells
        .get(m.solids.get(block).outer)
        .faces
        .iter()
        .copied()
        .find(|&fh| {
            has_face_on_plane(
                &m,
                block,
                Point3::from_array([5.0, 5.0, 4.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            ) && m.faces.get(fh).outer.half_edges.len() == 4
                && m.vertex_point(nacre_topo_first_vertex(&m, fh)).as_array()[2] == 4.0
        })
        .expect("top face");

    // An annular pocket: a square trench with an untouched island in the middle.
    //
    // `sq` takes the square's **world** bounds, which on this lid are also its frame coordinates:
    // the sketch origin is the world origin projected onto `z = 4` and the axes are `u = +x̂`,
    // `v = +ŷ`.
    let sq = |lo: f64, hi: f64| {
        vec![
            Point2::from_array([lo, hi]),
            Point2::from_array([lo, lo]),
            Point2::from_array([hi, lo]),
            Point2::from_array([hi, hi]),
        ]
    };
    let profile = Profile2d::with_holes(sq(1.0, 9.0), vec![sq(3.0, 7.0)]).unwrap();
    let OpOutput::PocketOnFace { solid, .. } = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: top,
            profile,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("PocketOnFace yields its own output")
    };
    m.rebuild_adjacency();

    // The trench is (8² − 4²) = 48 in plan, 1 deep. A hole wound the wrong way would remove the
    // whole 8² footprint (64) instead — or produce an invalid solid.
    assert!(nacre_validate::validate(&m).is_empty());
    assert!(
        (volume(&m, solid) - (400.0 - 48.0)).abs() < 1e-9,
        "{}",
        volume(&m, solid)
    );
}

/// The first vertex of a face's outer loop — a local helper, kept out of `common` because only
/// the pocket fixture above needs it.
fn nacre_topo_first_vertex(m: &Model, fh: Handle<nacre_topo::Face>) -> Handle<nacre_topo::Vertex> {
    let he = m.faces.get(fh).outer.half_edges[0];
    let [a, b] = m.edges.get(he.edge).vertices;
    if he.forward { a } else { b }
}

/// The whole point of the nesting step, through the front door: the caller hands over closed
/// rings in any order — no "this one is the hole" — and gets a solid with a hole in it.
#[test]
fn rings_in_any_order_become_a_donut_without_being_told_which_is_the_hole() {
    let sq = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    // Hole first, outer second — order carries no meaning.
    let profiles = from_rings(vec![sq(1.0, 3.0), sq(0.0, 4.0)]).unwrap();
    assert_eq!(profiles.len(), 1, "one body");

    let mut m = Model::new();
    let d = extrude(&mut m, profiles.into_iter().next().unwrap(), 1.0);

    assert!(nacre_validate::validate(&m).is_empty());
    assert!((volume(&m, d) - 12.0).abs() < 1e-12, "{}", volume(&m, d));
}

/// Rings inside a hole come back as their own body, and each extrudes on its own — this is the
/// "island" the syntax promises, assembled by the caller as a loop over the returned profiles.
#[test]
fn an_island_extrudes_as_a_second_body() {
    let sq = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    let profiles = from_rings(vec![sq(0.0, 9.0), sq(1.0, 8.0), sq(2.0, 7.0)]).unwrap();
    assert_eq!(profiles.len(), 2);

    let mut m = Model::new();
    let mut total = 0.0;
    for p in profiles {
        let s = extrude(&mut m, p, 1.0);
        total += volume(&m, s);
    }
    assert!(nacre_validate::validate(&m).is_empty());
    // ring 0..9 minus hole 1..8 = 81 − 49 = 32, plus the island 2..7 = 25.
    assert!((total - (32.0 + 25.0)).abs() < 1e-12, "{total}");
}

// --- the profile contract: what the kernel refuses to build from ---

fn sq(a: f64, b: f64) -> Vec<Point2> {
    vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ]
}

fn try_extrude(profile: Profile2d) -> Result<(), nacre_ops::OpError> {
    apply(
        &mut Model::new(),
        &Operation::Extrude {
            frame: SketchFrame::world(&Model::new(), Axis::Z),
            profile,
            dist: 1.0,
        },
    )
    .map(|_| ())
}

/// **The measurement this whole contract exists for.** Before the check, each of these built a
/// solid that `validate` called clean and `nacre-step` happily exported. The bowtie's volume came
/// out `NaN`; the two nesting mistakes came out as *plausible numbers* — 12.0 for a hole that
/// lies nowhere near the outline (it should be 16), and 20.0 where even-odd says 52. A wrong
/// number that looks right is the failure mode the kernel exists to prevent.
#[test]
fn profiles_the_kernel_used_to_build_silently_wrong_are_now_refused() {
    let bowtie = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([4.0, 4.0]),
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    assert!(matches!(
        try_extrude(Profile2d::polygon(bowtie).unwrap()),
        Err(nacre_ops::OpError::SelfIntersectingProfile { .. })
    ));
    // A hole that misses the outline entirely: it used to be subtracted anyway.
    assert!(matches!(
        try_extrude(Profile2d::with_holes(sq(0.0, 4.0), vec![sq(10.0, 12.0)]).unwrap()),
        Err(nacre_ops::OpError::HoleNotInsideOuter { hole: 0 })
    ));
    // A hole inside a hole is an island — material, not a second subtraction.
    assert!(matches!(
        try_extrude(
            Profile2d::with_holes(sq(0.0, 10.0), vec![sq(1.0, 9.0), sq(3.0, 7.0)]).unwrap()
        ),
        Err(nacre_ops::OpError::NestedHole { .. })
    ));
}

/// A ring that touches itself without crossing. There is no strict inside at the pinch, and the
/// prism it used to build had a believable volume (14.0) and a clean `validate`.
#[test]
fn a_pinched_ring_is_refused() {
    let pinch = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 2.0]),
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([4.0, 4.0]),
        Point2::from_array([2.0, 2.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    assert!(matches!(
        try_extrude(Profile2d::polygon(pinch).unwrap()),
        Err(nacre_ops::OpError::SelfIntersectingProfile { .. })
    ));
}

/// Holes touching the outline, or each other, have no unambiguous inside either — and the pad and
/// pocket path shares the same gate, so it refuses them too.
#[test]
fn touching_rings_are_refused_on_every_profile_entry_point() {
    let touching = Profile2d::with_holes(sq(0.0, 4.0), vec![sq(0.0, 2.0)]).unwrap();
    assert!(matches!(
        try_extrude(touching.clone()),
        Err(nacre_ops::OpError::ProfileRingsMeet { .. })
    ));

    // Same profile, arriving through `pad` on a face of an existing solid.
    let mut m = Model::new();
    let base = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([10.0, 10.0, 1.0]),
    );
    m.rebuild_adjacency();
    let top = *m.shells.get(m.solids.get(base).outer).faces.last().unwrap();
    assert!(matches!(
        apply(
            &mut m,
            &Operation::PadOnFace {
                face: top,
                profile: touching,
                dist: 1.0,
            },
        ),
        Err(nacre_ops::OpError::ProfileRingsMeet { .. })
    ));
}

/// The sketch layer refuses before it classifies, because it must: containment is decided by
/// even-odd parity, which only means "inside" on a simple ring. Reported by **point**, not edge
/// index — a point is what an editor can mark, and the enum's other rejections carry points too.
/// The bowtie here is a stated ring of four lines, exactly how it reaches the front door.
#[test]
fn a_self_crossing_outline_is_refused_before_the_rings_are_sorted() {
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let bowtie = vec![
        line(p(0.0, 0.0), p(4.0, 4.0)),
        line(p(4.0, 4.0), p(4.0, 0.0)),
        line(p(4.0, 0.0), p(0.0, 4.0)),
        line(p(0.0, 4.0), p(0.0, 0.0)),
    ];
    assert!(matches!(
        stated(bowtie),
        Err(nacre_ops::SketchError::RingSelfIntersects { .. })
    ));

    // With a second, perfectly good ring alongside it, the self-intersection still wins — the
    // nesting pass never runs on a ring whose inside is undefined.
    let crossing = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([4.0, 4.0]),
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    assert_eq!(
        from_rings(vec![sq(-10.0, -5.0), crossing]),
        Err(nacre_ops::SketchError::RingSelfIntersects {
            ring: 1,
            at: [[2.0, 2.0], [2.0, 2.0]],
        })
    );
}

/// The other half of the contract: everything legitimate still goes through. A reflex outline, a
/// donut, and a flat (collinear) corner are all simple polygons, and `check` must not reject
/// them. The flat corner does not survive as *data* — construction dissolves it (that is S3's
/// normalization, asserted below so the `Ok` is not vacuous) — but the author's input is legal.
#[test]
fn the_contract_does_not_bite_legitimate_profiles() {
    let l = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([4.0, 2.0]),
        Point2::from_array([2.0, 2.0]),
        Point2::from_array([2.0, 4.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    assert_eq!(Profile2d::polygon(l).unwrap().check(), Ok(()));
    assert_eq!(donut_profile().check(), Ok(()));
    let flat = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 0.0]), // mid-run on a straight edge
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([4.0, 4.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    let p = Profile2d::polygon(flat).unwrap();
    assert_eq!(p.check(), Ok(()));
    assert_eq!(
        p.outer().vertices().len(),
        4,
        "the flat corner is dissolved at construction, not merely tolerated"
    );
}

/// ★ **A coordinate outside the decimal window is a named error at construction** — S3's other
/// half. It used to fall silently to the f64 path: the prism still built, but recorded no exact
/// points, and its surfaces could not survive a motion undemoted. Both profile entry points name
/// it now.
#[test]
fn a_dimension_outside_the_decimal_window_is_refused_at_construction() {
    let ring = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([1e300, 0.0]), // no decimal this side of i128 spells it
        Point2::from_array([0.0, 4.0]),
    ];
    assert_eq!(
        Profile2d::polygon(ring.clone()).unwrap_err(),
        nacre_ops::OpError::ProfileOutsideDecimalWindow {
            ring: nacre_ops::ProfileRing::Outer,
            point: 1,
        }
    );
    assert_eq!(
        from_rings(vec![ring]).unwrap_err(),
        nacre_ops::SketchError::OutsideDecimalWindow { at: [1e300, 0.0] }
    );
}

/// ★ **The truth is the decimal the author wrote, and it wins over the binary carrier.**
///
/// `(0, 0.1) → (0.1, 0.2) → (0.2, 0.3)` is collinear in decimal (slope one) but *not* in the f64
/// binary values (`0.3` is not exactly `3 × 0.1bin`), so an f64-exact orient2d calls the corner
/// a real one. Construction judges the rational truth and dissolves it — the profile the kernel
/// keeps is what the author meant, not what the carrier rounded to.
#[test]
fn a_corner_flat_in_decimal_but_not_in_binary_is_dissolved() {
    let p = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.1]),
        Point2::from_array([0.1, 0.2]), // decimal-collinear midpoint, binary-bent corner
        Point2::from_array([0.2, 0.3]),
        Point2::from_array([0.2, 4.0]),
        Point2::from_array([-4.0, 4.0]),
    ])
    .unwrap();
    assert_eq!(
        p.outer().vertices().len(),
        4,
        "the decimal truth decided, so the midpoint is gone"
    );
}

/// **Circles and arcs both build; what the builder cannot stand is refused by name.** The
/// vocabulary accepts a `3-4-5` lens as a region, but its arcs are no quarter turns, so the
/// exact winding cannot read it; a leaf of two quarter arcs winds fine and fails on its corners
/// instead — a point on two cylinders is a definition the kernel does not have.
#[test]
fn a_circle_and_a_slot_extrude_and_a_lens_is_refused_by_name() {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let r = |n: i128| nacre_scalar::Rat::from_int(n);
    let circle = stated(vec![circle(p2(0.0, 0.0), 1.0)]).unwrap();
    let slot = stated(vec![
        line(p2(0.0, -1.0), p2(4.0, -1.0)),
        arc_turns(p2(4.0, 0.0), p2(4.0, -1.0), 2),
        line(p2(4.0, 1.0), p2(0.0, 1.0)),
        arc_turns(p2(0.0, 0.0), p2(0.0, 1.0), 2),
    ])
    .unwrap();
    let lens = stated(vec![
        arc_rat([r(3), r(4)], [r(0), r(0)], [r(6), r(0)], true),
        arc_rat([r(3), r(-4)], [r(6), r(0)], [r(0), r(0)], true),
    ])
    .unwrap();
    for profiles in [circle, slot] {
        let mut m = Model::new();
        let frame = SketchFrame::world(&m, Axis::Z);
        let out = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: profiles[0].clone(),
                dist: 1.0,
            },
        );
        assert!(matches!(out, Ok(OpOutput::Extrude { .. })), "{out:?}");
    }
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    assert_eq!(
        apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: lens[0].clone(),
                dist: 1.0
            }
        ),
        Err(nacre_ops::OpError::ArcSweepNotQuarterTurn)
    );
    // A leaf: two quarter arcs of different circles between (0,0) and (5,5).
    let leaf = stated(vec![
        arc_turns(p2(5.0, 0.0), p2(0.0, 0.0), -1), // → (5,5), bulging up-left
        arc_turns(p2(0.0, 5.0), p2(5.0, 5.0), -1), // → (0,0), bulging down-right
    ])
    .unwrap();
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    assert_eq!(
        apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: leaf[0].clone(),
                dist: 1.0
            }
        ),
        Err(nacre_ops::OpError::ArcsMeetAtVertex)
    );
}

/// **A slot prism exports**: its cylinder walls and arc edges reach STEP as the surfaces and
/// curves they are.
#[test]
fn a_slot_prism_tessellates_and_exports() {
    let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
    let slot = stated(vec![
        line(p2(0.0, -1.0), p2(4.0, -1.0)),
        arc_turns(p2(4.0, 0.0), p2(4.0, -1.0), 2),
        line(p2(4.0, 1.0), p2(0.0, 1.0)),
        arc_turns(p2(0.0, 0.0), p2(0.0, 1.0), 2),
    ])
    .unwrap();
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: slot[0].clone(),
            dist: 1.0,
        },
    )
    .expect("the slot extrudes");
    assert!(nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok());
    let step = nacre_step::to_step(&m).unwrap();
    for needle in ["MANIFOLD_SOLID_BREP", "CYLINDRICAL_SURFACE", "CIRCLE"] {
        assert!(step.contains(needle), "missing {needle}");
    }
}

// --- the radius is a square (open item 25, step 1) ---

/// ★ **The truth of a circle is its squared radius.** A quarter arc about the origin from
/// `(1, 1)` to `(−1, 1)` has `r² = 2` and no rational radius; until 2026-09-15 the sketch door
/// refused it by name (`ArcRadiusNotRational`) although no predicate ever needed the radius
/// unsquared. Now it is a ring — the circular segment above the chord `y = 1` — and extrudes to
/// a solid whose volume is `(π/2 − 1)·h`, its lateral cylinder's cache carrying `√2` correctly
/// rounded.
#[test]
fn a_non_pythagorean_arc_is_stated_and_extruded() {
    use nacre_ops::{Edge2d, Ring2d, arc_to_rat, arc_turns, from_paths};
    use nacre_scalar::Rat;
    let r = |n: i128| Rat::from_int(n);
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    // Both step doors state the same arc: the quarter turn lands where the stated end is.
    let (turned, end) = arc_turns(p(0.0, 0.0), p(1.0, 1.0), 1).expect("r² = 2 is a circle");
    assert_eq!(end, [r(-1), r(1)]);
    let stated =
        arc_to_rat([r(0), r(0)], [r(1), r(1)], [r(-1), r(1)], true).expect("on the circle");
    assert_eq!(turned, stated);
    assert!(matches!(stated, Edge2d::Arc { r2, .. } if r2 == r(2)));
    // The segment: the arc over the top, the chord back along `y = 1` — counter-clockwise.
    let ring = Ring2d::new(
        vec![[r(1), r(1)], [r(-1), r(1)]],
        vec![stated, Edge2d::Line],
    )
    .expect("a ring of one arc and one line");
    let mut profiles = from_paths(vec![ring]).expect("one region");
    let profile = profiles.pop().expect("one profile");
    assert!(profiles.is_empty());

    let mut m = Model::new();
    let s = extrude(&mut m, profile, 1.0);
    let violations = nacre_validate::validate(&m);
    assert!(violations.is_empty(), "{violations:?}");
    let want = std::f64::consts::FRAC_PI_2 - 1.0; // (r²/2)(θ − sin θ) with r² = 2, θ = π/2
    let got = volume(&m, s);
    assert!((got - want).abs() < 1e-9, "volume {got}, want {want}");
    // The lateral face's cache is the truth realized: √2, correctly rounded.
    let radii: Vec<f64> = m
        .faces
        .iter()
        .filter_map(|(_, f)| match m.surface(f.surface) {
            nacre_geom::Surface::Cylinder(c) => Some(c.radius()),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 1, "one lateral cylinder");
    assert_eq!(radii[0].to_bits(), 2f64.sqrt().to_bits());
}

/// ★ The boolean road carries `r² = 2` end to end: the segment prism above, cut by a box covering
/// `x ≥ 0`, leaves the left half of the segment — validate clean, volume `(π/2 − 1)/2`. Written as
/// the cell's probe («which named refusal is the first wall?») and promoted to a lock when it came
/// back green; a refusal is still reported in the panic so a regression names itself.
#[test]
fn a_non_pythagorean_prism_is_cut_by_a_box() {
    use nacre_ops::{Edge2d, Ring2d, arc_to_rat, from_paths};
    use nacre_scalar::Rat;
    let r = |n: i128| Rat::from_int(n);
    let stated =
        arc_to_rat([r(0), r(0)], [r(1), r(1)], [r(-1), r(1)], true).expect("on the circle");
    let ring = Ring2d::new(
        vec![[r(1), r(1)], [r(-1), r(1)]],
        vec![stated, Edge2d::Line],
    )
    .unwrap();
    let profile = from_paths(vec![ring]).unwrap().pop().unwrap();
    let mut m = Model::new();
    let prism = extrude(&mut m, profile, 1.0);
    // A box covering `x ≥ 0`: the cut leaves the left half of the segment.
    let cutter = m.add_cuboid(
        Point3::from_array([0.0, 0.0, -1.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    match apply(
        &mut m,
        &Operation::Boolean {
            kind: BoolKind::Cut,
            a: prism,
            b: cutter,
        },
    ) {
        Ok(OpOutput::Boolean { solids }) => {
            let violations = nacre_validate::validate(&m);
            assert!(violations.is_empty(), "{violations:?}");
            let total: f64 = solids.iter().map(|&s| volume(&m, s)).sum();
            let want = (std::f64::consts::FRAC_PI_2 - 1.0) / 2.0;
            assert!((total - want).abs() < 1e-9, "volume {total}, want {want}");
        }
        Ok(other) => panic!("unexpected output {other:?}"),
        Err(e) => panic!("the boolean road refused r² = 2: {e:?}"),
    }
}
