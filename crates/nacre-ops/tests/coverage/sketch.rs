//! Profiles with holes — extruding a region that is not simply connected.
//!
//! The kernel could already *produce* a donut, by cutting a bar out of a box; what it could not do
//! was **construct** one. That distinction is the point: the boolean route makes every corner a
//! `Discovered` intersection with a measured tolerance, while extruding the profile directly makes
//! them `Constructed` — exact by construction (overview 절대원칙 4). So the strongest check here
//! is not a number typed by hand but the two producers agreeing.

use crate::common::*;
use nacre_math::{Point2, Point3};
use nacre_ops::{
    BoolKind, Edge2d, OpOutput, Operation, Profile2d, SketchPlane, apply, from_edges, from_rings,
};
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Solid};

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
    Profile2d::with_holes(sq(0.0, 4.0), vec![sq(1.0, 3.0)])
}

fn extrude(m: &mut Model, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            plane: SketchPlane::world_xy(),
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
                for vh in m.edges.get(he.edge).bounds.iter().flatten() {
                    assert!(
                        matches!(m.vertices.get(*vh).origin, Origin::Constructed),
                        "a swept vertex carries no tolerance"
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

/// A hole is only a hole if the mesh and the exchange format agree it is one. Tessellation
/// bridges cap holes, and STEP writes them as inner face bounds.
#[test]
fn a_donut_prism_tessellates_and_exports() {
    let mut m = Model::new();
    extrude(&mut m, donut_profile(), 1.0);

    assert!(
        nacre_tess::to_obj(&m).is_ok(),
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
    );
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
                && m.vertices
                    .get(nacre_topo_first_vertex(&m, fh))
                    .point
                    .as_array()[2]
                    == 4.0
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
    let profile = Profile2d::with_holes(sq(1.0, 9.0), vec![sq(3.0, 7.0)]);
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
    let [a, b] = m.edges.get(he.edge).bounds.expect("bounded");
    if he.forward { a } else { b }
}

/// The whole point of the nesting step, through the front door: the caller hands over loose
/// closed rings — no "this one is the hole" — and gets a solid with a hole in it.
#[test]
fn loose_rings_become_a_donut_without_being_told_which_is_the_hole() {
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

/// The front door the syntax actually describes: hand over drawn segments, in any order, and get
/// the solid. Nothing declares the hole and nothing declares the order.
#[test]
fn drawn_segments_become_a_donut() {
    let ring = |a: f64, b: f64| {
        let p = |x: f64, y: f64| Point2::from_array([x, y]);
        vec![
            Edge2d::line(p(a, a), p(b, a)),
            Edge2d::line(p(b, b), p(b, a)), // backwards on purpose
            Edge2d::line(p(b, b), p(a, b)),
            Edge2d::line(p(a, b), p(a, a)),
        ]
    };
    let mut edges = ring(1.0, 3.0); // the hole, drawn first
    edges.extend(ring(0.0, 4.0));

    let profiles = from_edges(edges).unwrap();
    assert_eq!(profiles.len(), 1);

    let mut m = Model::new();
    let d = extrude(&mut m, profiles.into_iter().next().unwrap(), 1.0);
    assert!(nacre_validate::validate(&m).is_empty());
    assert!((volume(&m, d) - 12.0).abs() < 1e-12, "{}", volume(&m, d));
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
            plane: SketchPlane::world_xy(),
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
        try_extrude(Profile2d::polygon(bowtie)),
        Err(nacre_ops::OpError::SelfIntersectingProfile { .. })
    ));
    // A hole that misses the outline entirely: it used to be subtracted anyway.
    assert!(matches!(
        try_extrude(Profile2d::with_holes(sq(0.0, 4.0), vec![sq(10.0, 12.0)])),
        Err(nacre_ops::OpError::HoleNotInsideOuter { hole: 0 })
    ));
    // A hole inside a hole is an island — material, not a second subtraction.
    assert!(matches!(
        try_extrude(Profile2d::with_holes(
            sq(0.0, 10.0),
            vec![sq(1.0, 9.0), sq(3.0, 7.0)]
        )),
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
        try_extrude(Profile2d::polygon(pinch)),
        Err(nacre_ops::OpError::SelfIntersectingProfile { .. })
    ));
}

/// Holes touching the outline, or each other, have no unambiguous inside either — and the pad and
/// pocket path shares the same gate, so it refuses them too.
#[test]
fn touching_rings_are_refused_on_every_profile_entry_point() {
    let touching = Profile2d::with_holes(sq(0.0, 4.0), vec![sq(0.0, 2.0)]);
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
/// index — `from_edges` chains rings in walk order, so an index would name nothing the author
/// wrote. The bowtie here is drawn as four loose segments, exactly how it reaches the front door.
#[test]
fn a_self_crossing_outline_is_refused_before_the_rings_are_sorted() {
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let bowtie = vec![
        Edge2d::line(p(0.0, 0.0), p(4.0, 4.0)),
        Edge2d::line(p(4.0, 4.0), p(4.0, 0.0)),
        Edge2d::line(p(4.0, 0.0), p(0.0, 4.0)),
        Edge2d::line(p(0.0, 4.0), p(0.0, 0.0)),
    ];
    assert!(matches!(
        from_edges(bowtie),
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
/// donut, and a flat (collinear) corner are all simple polygons, and `check` must not touch them.
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
    assert_eq!(Profile2d::polygon(l).check(), Ok(()));
    assert_eq!(donut_profile().check(), Ok(()));
    let flat = vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([2.0, 0.0]), // mid-run on a straight edge
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([4.0, 4.0]),
        Point2::from_array([0.0, 4.0]),
    ];
    assert_eq!(Profile2d::polygon(flat).check(), Ok(()));
}
