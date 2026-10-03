//! The closed cylinder's b-rep — the seam model — as the product builds it: a whole circle
//! sketched and extruded (`nacre_ops::fixtures`). Two seam vertices, two full-circle rims, one
//! straight seam the lateral uses twice, two caps; and what the kernel keeps true of it (the
//! adjacency, the caps' names, the edge cache's regeneration). Validate-clean lives in
//! `nacre-validate`.

use nacre_geom::Curve;
use nacre_math::{Point3, Vector3};
use nacre_ops::fixtures::{CylinderSolid, cylinder};
use nacre_store::Handle;
use nacre_topo::{Adjacency, Model, Motion, Orientation, Surface, Vertex};
use proptest::prelude::*;

fn z_cylinder(r: f64, h: f64) -> (Model, CylinderSolid) {
    let mut m = Model::new();
    let c = cylinder(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        r,
        h,
    );
    m.rebuild_adjacency();
    (m, c)
}

/// Every edge used exactly twice, in opposite senses — the seam included, by one face.
fn every_edge_used_twice_opposite(m: &Model) -> Result<(), String> {
    let adj = Adjacency::rebuild(m);
    for (e, uses) in &adj.edge_uses {
        if uses.len() != 2 || uses[0].1 == uses[1].1 {
            return Err(format!("edge {} used {uses:?}", e.index()));
        }
    }
    Ok(())
}

#[test]
fn cylinder_counts_and_euler() {
    let (m, _) = z_cylinder(2.0, 5.0);
    assert_eq!(m.vertex_count(), 2);
    assert_eq!(m.edge_count(), 3);
    assert_eq!(m.face_count(), 3);
    assert_eq!(m.shell_count(), 1);
    assert_eq!(m.solid_count(), 1);
    // Euler χ = V − E + F = 2 (one shell, genus 0, no inner loops).
    assert_eq!(
        m.vertex_count() as i64 - m.edge_count() as i64 + m.face_count() as i64,
        2
    );
}

#[test]
fn cylinder_caps_face_outward() {
    let (m, c) = z_cylinder(2.0, 5.0);
    for (cap, want) in [(c.base, -1.0), (c.top, 1.0)] {
        let f = m.face(cap);
        let nacre_geom::Surface::Plane(p) = m.surface_cache(f.surface) else {
            panic!("a cap is planar")
        };
        let n = p
            .normal()
            .as_array()
            .map(|x| x * f64::from(f.orientation.sign()));
        assert_eq!(n, [0.0, 0.0, want], "the cap faces out along the axis");
    }
}

/// ★★ **A seam vertex states its two carriers** — `OnSeam([lateral, its own cap])`.
/// The pair pins the rim circle; the cylinder's `ref_dir` truth pins the point
/// (see `Vertex::OnSeam`). The lateral surface must be the cylinder, the other its cap.
#[test]
fn a_seam_vertex_states_its_rim_carriers() {
    let (m, _) = z_cylinder(2.0, 5.0);
    let vs: Vec<_> = (0..m.vertex_count() as u32)
        .filter_map(|i| m.vertex_handle_at(i))
        .collect();
    assert_eq!(vs.len(), 2, "a cylinder has exactly its two seam vertices");
    for vh in vs {
        let Vertex::OnSeam([a, b]) = *m.vertex(vh) else {
            panic!("a seam vertex carries OnSeam, got {:?}", m.vertex(vh))
        };
        assert!(
            matches!(m.surface_cache(a), nacre_geom::Surface::Cylinder(_)),
            "first carrier is the lateral cylinder"
        );
        let nacre_geom::Surface::Plane(p) = m.surface_cache(b) else {
            panic!("second carrier is the cap plane")
        };
        assert_eq!(
            p.distance(m.vertex_point(vh)),
            0.0,
            "the vertex's cap is the one it lies on"
        );
    }
}

/// The novel topology: the seam edge is used twice by the **same** lateral face with opposite
/// orientation. Every edge is still used exactly twice, opposite.
#[test]
fn cylinder_seam_edge_is_self_adjacent() {
    let (m, c) = z_cylinder(2.0, 5.0);
    every_edge_used_twice_opposite(&m).unwrap();
    let adj = Adjacency::rebuild(&m);
    let seam = (0..m.edge_count() as u32)
        .filter_map(|i| m.edge_handle_at(i))
        .find(|&eh| m.edge(eh).surfaces[0] == m.edge(eh).surfaces[1])
        .expect("a seam line edge exists");
    let uses = &adj.edge_uses[&seam];
    assert_eq!(uses[0].0, uses[1].0, "both uses are the same face");
    assert_eq!(uses[0].0, c.lateral, "and it is the lateral");
}

/// ★★ The «discard and regenerate» warrant: the edge-curve cache rebuilt from the carriers and
/// endpoints is bit-identical to the one `push_edge` filled eagerly — proof that nothing in it
/// was truth. A tilted cylinder beside a box, so circles, a seam and straight edges all regenerate.
#[test]
fn edge_cache_discard_and_regenerate_bit_identical() {
    let mut m = Model::new();
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-2.0, 1.0, 0.0]),
        Point3::from_array([3.0, 4.0, 10.0]),
    );
    cylinder(
        &mut m,
        Point3::from_array([8.0, 0.0, 0.0]),
        Vector3::from_array([0.3, -0.4, 1.0]),
        1.25,
        2.5,
    );
    let snapshot = |m: &Model| -> Vec<Curve> {
        (0..m.edge_count() as u32)
            .filter_map(|i| m.edge_handle_at(i))
            .map(|eh| m.edge_curve(eh).clone())
            .collect()
    };
    let before = snapshot(&m);
    m.rebuild_edge_cache();
    // ⚠ True of a model nothing has refined. `Model::refine_vertex_cache` moves coordinates,
    // and a rebuild after *that* is not a no-op — which is the whole reason the paying caller
    // re-derives. This asserts the derivation is stable, not that rebuilding is always free.
    assert_eq!(
        before,
        snapshot(&m),
        "regeneration must reproduce the cache bit for bit"
    );
}

/// ★★ **A rim's frame is its cylinder's, to the bit** — normal, `ref_dir` and radius are the
/// cylinder cache's, which is the truth's correctly rounded frame. Asked where re-deriving it in
/// `f64` moves it: a cylinder on a tilted frame node (axis `(1,1,1)`) turned 37° about `z` — its
/// seam direction has a component whose truth is exactly `0`, which a renormalization brings back
/// as a residue — and a `+z` cylinder turned 37° about `y` then 11° about `x`.
///
/// ★ The sweep also counts the rims whose frame a renormalization *would* move, and requires one:
/// a fixture whose frame survives `Circle::from_center_normal` bit for bit cannot tell the two
/// constructors apart, and this lock would pass with either.
#[test]
fn a_rims_frame_is_its_cylinders_cache() {
    use crate::fixtures::rot_iso;
    use nacre_exact::Axis;
    let cases: [(&str, [f64; 3], Vec<nacre_exact::Isometry>); 2] = [
        (
            "frame node (1,1,1), turned 37 about z",
            [1.0, 1.0, 1.0],
            vec![rot_iso(Axis::Z, 37)],
        ),
        (
            "turned 37 about y, 11 about x",
            [0.0, 0.0, 1.0],
            vec![rot_iso(Axis::Y, 37), rot_iso(Axis::X, 11)],
        ),
    ];
    let mut would_move = 0;
    for (name, axis, turns) in cases {
        let mut m = Model::new();
        let mut s = cylinder(
            &mut m,
            Point3::origin(),
            Vector3::from_array(axis),
            0.3,
            1.1,
        )
        .solid;
        for iso in turns {
            s = crate::fixtures::xf(&mut m, s, iso);
        }
        m.rebuild_adjacency();
        let bits = |v: Vector3| v.as_array().map(f64::to_bits);
        let mut rims = 0;
        for &fh in &m.shell(m.solid(s).outer).faces {
            for he in &m.face(fh).outer.half_edges {
                let Curve::Circle(c) = m.edge_curve(he.edge) else {
                    continue;
                };
                let cyl = m
                    .edge(he.edge)
                    .surfaces
                    .into_iter()
                    .find(|&h| matches!(m.surface(h), Surface::Cylinder { .. }))
                    .expect("a rim has a cylinder carrier");
                let nacre_geom::Surface::Cylinder(cc) = m.surface_cache(cyl) else {
                    unreachable!("a cylinder truth has a cylinder cache")
                };
                assert_eq!(
                    bits(c.normal()),
                    bits(cc.axis().direction()),
                    "{name}: normal"
                );
                assert_eq!(bits(c.ref_dir()), bits(cc.ref_dir()), "{name}: ref_dir");
                assert_eq!(
                    c.radius().to_bits(),
                    cc.radius().to_bits(),
                    "{name}: radius"
                );
                let again = nacre_geom::Circle::from_center_normal(
                    c.center(),
                    cc.axis().direction(),
                    cc.ref_dir(),
                    cc.radius(),
                )
                .expect("a nondegenerate frame");
                if bits(again.normal()) != bits(c.normal())
                    || bits(again.ref_dir()) != bits(c.ref_dir())
                {
                    would_move += 1;
                }
                rims += 1;
            }
        }
        assert!(rims >= 2, "{name}: the sweep met {rims} rim uses");
    }
    assert!(
        would_move > 0,
        "no rim here has a frame a renormalization moves — the lock cannot tell the constructors apart"
    );
}

/// ★★ **A stated rim's centre is the truth's nearest `f64`** — the cylinder's axis met with its
/// cap, exactly, then rounded. The oracle is this test's own: the fixture's decimals (the base and
/// the height, read as the decimals they spell) and an axis whose unit vector is rational by hand
/// (`3-4-0 / 5`, `0-3-4 / 5`, `12-0-5 / 13`, and `+z`), so the bottom rim's centre is the base and
/// the top's is `base + h·û`. Each of these axes has a rational frame across it, so the cylinder is
/// stated in the world — asserted, because a rational unit axis alone is not enough: `2-3-6 / 7`
/// has no rational perpendicular frame and is built on a frame node, where the centre is not yet
/// the truth's.
///
/// ★ The sweep also counts the rims where the caches' `f64` meet (`line_plane` of the cylinder
/// cache's axis and the cap cache) misses the truth, and requires one — otherwise this lock could
/// not tell the two roads apart.
#[test]
fn a_stated_rims_centre_is_the_truths_nearest() {
    use nacre_exact::Rat;
    let r = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    let cases: [([f64; 3], [Rat; 3]); 4] = [
        ([3.0, 4.0, 0.0], [r(3, 5), r(4, 5), r(0, 1)]),
        ([0.0, 3.0, 4.0], [r(0, 1), r(3, 5), r(4, 5)]),
        ([12.0, 0.0, 5.0], [r(12, 13), r(0, 1), r(5, 13)]),
        ([0.0, 0.0, 1.0], [r(0, 1), r(0, 1), r(1, 1)]),
    ];
    let bases: [[f64; 3]; 2] = [[0.1, 0.7, 0.3], [1.5, -2.25, 12.75]];
    let height = 1.1;
    let mut caches_miss = 0;
    for (axis, unit) in cases {
        for base in bases {
            let mut m = Model::new();
            let s = cylinder(
                &mut m,
                Point3::from_array(base),
                Vector3::from_array(axis),
                0.3,
                height,
            )
            .solid;
            m.rebuild_adjacency();
            let b = base.map(|x| Rat::from_decimal(x).expect("a decimal"));
            let h = Rat::from_decimal(height).expect("a decimal");
            let top: [Rat; 3] = core::array::from_fn(|k| {
                b[k].checked_add(h.checked_mul(unit[k]).unwrap()).unwrap()
            });
            let mut want: Vec<[u64; 3]> = [b, top]
                .iter()
                .map(|p| p.map(|x| x.to_f64().to_bits()))
                .collect();
            want.sort();
            let mut got = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for &fh in &m.shell(m.solid(s).outer).faces {
                for he in &m.face(fh).outer.half_edges {
                    let Curve::Circle(c) = m.edge_curve(he.edge) else {
                        continue;
                    };
                    if !seen.insert(he.edge) {
                        continue;
                    }
                    let centre = c.center().as_array().map(f64::to_bits);
                    got.push(centre);
                    let [s0, s1] = m.edge(he.edge).surfaces;
                    let (cyl, cap) = match m.surface(s0) {
                        Surface::Cylinder { .. } => (s0, s1),
                        Surface::Plane { .. } => (s1, s0),
                    };
                    assert!(
                        m.world_cylinder_def(cyl).is_some(),
                        "axis {axis:?}: the premise — this cylinder is stated in the world"
                    );
                    let (nacre_geom::Surface::Cylinder(cc), nacre_geom::Surface::Plane(pc)) =
                        (m.surface_cache(cyl), m.surface_cache(cap))
                    else {
                        unreachable!("a rim is a cylinder against a plane")
                    };
                    let meet = nacre_geom::intersect::line_plane(&cc.axis(), pc)
                        .expect("the axis crosses the cap");
                    if meet.as_array().map(f64::to_bits) != centre {
                        caches_miss += 1;
                    }
                }
            }
            got.sort();
            assert_eq!(got, want, "axis {axis:?} base {base:?}: rim centres");
        }
    }
    assert!(
        caches_miss > 0,
        "every rim here is one the caches' f64 meet already gets right — the lock cannot tell the roads apart"
    );
}

/// A cylinder's caps state points and a name, and a cap on the plane of another producer's face
/// **is** that face's plane — one plane, one handle, across two producers.
#[test]
fn a_cylinders_caps_record_points_and_intern_with_a_coplanar_face() {
    let mut m = Model::new();
    let cuboid = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([-1.0, -1.0, 2.0]),
    );
    let c = cylinder(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    m.rebuild_adjacency();
    let caps = [c.base, c.top].map(|f| m.face(f).surface);
    for s in &caps {
        assert!(
            matches!(m.surface(*s), Surface::Plane { .. }),
            "a cap without truth"
        );
        assert!(m.surface_name.contains_key(s), "a cap without a name");
    }
    // The top cap lies on `z = 2`, the same plane as the box's top face.
    let box_top = m
        .shell(m.solid(cuboid).outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .find(|s| {
            m.surface_name
                .get(s)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, -2.0])
        })
        .expect("the box top names z = 2");
    assert_eq!(
        caps[1], box_top,
        "the coplanar cap and box face must intern to one handle"
    );
}

/// ★★ **A nameless `Through` statement interns by its statement.** Two seam vertices and a box
/// corner span a plane no name road reaches (`plane_name_through` is `None`), so the door stores
/// the statement itself — once — and a second statement of it is the same handle: facing the
/// other way reports `flipped`, and a different motion is a different statement. The door's
/// contract ("store the statement, once") holds for every nameless reason identically, which is
/// what is pinned here.
#[test]
fn a_nameless_through_statement_interns_by_its_statement() {
    let mut m = Model::new();
    cylinder(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    );
    nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([4.0, 0.0, 0.0]),
        Point3::from_array([5.0, 1.0, 1.0]),
    );
    let mut vs: Vec<Handle<Vertex>> = Vec::new();
    let mut i = 0u32;
    while let Some(h) = m.vertex_handle_at(i) {
        i += 1;
        if matches!(*m.vertex(h), Vertex::OnSeam(_)) && vs.len() < 2 {
            vs.push(h);
        } else if matches!(*m.vertex(h), Vertex::ThreePlane(_)) && vs.len() == 2 {
            vs.push(h);
            break;
        }
    }
    let mut triple: [Handle<Vertex>; 3] = [vs[0], vs[1], vs[2]];
    triple.sort_by_key(|v| v.index());
    assert!(
        m.plane_name_through(triple).is_none(),
        "the fixture must be nameless, or this test measures the name road"
    );

    let cache = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.5]),
        Vector3::from_array([0.0, 0.0, -1.0]),
    )
    .unwrap();
    let before = m.surface_count();
    let (h, flipped) = m.push_plane_through(cache, triple, None, Orientation::Forward);
    assert!(!flipped);
    assert!(
        !m.surface_name.contains_key(&h),
        "a nameless statement must not invent a name"
    );
    assert_eq!(m.surface_count(), before + 1, "stored once");

    let (again, flipped_same) = m.push_plane_through(cache, triple, None, Orientation::Forward);
    assert_eq!(h, again, "one statement, one handle");
    assert!(!flipped_same);
    let reversed = nacre_geom::Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 0.5]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let (still, flipped_now) = m.push_plane_through(reversed, triple, None, Orientation::Reversed);
    assert_eq!(h, still, "direction is not part of the statement");
    assert!(
        flipped_now,
        "but a statement facing the other way is reported"
    );
    assert_eq!(m.surface_count(), before + 1, "and nothing new was stored");

    // A different motion is a different statement — the same rule the name key keeps.
    let node = m.push_motion(
        Motion::Translate {
            offset: [
                nacre_exact::Rat::from_int(1),
                nacre_exact::Rat::from_int(0),
                nacre_exact::Rat::from_int(0),
            ],
        },
        None,
    );
    let (moved, _) = m.push_plane_through(cache, triple, Some(node), Orientation::Forward);
    assert_ne!(h, moved, "the motion belongs in the statement key");
}

/// The kinds of a solid's edges: (circles, lines).
fn edge_kinds(m: &Model, s: Handle<nacre_topo::Solid>) -> (usize, usize) {
    let mut seen = std::collections::HashSet::new();
    let (mut circles, mut lines) = (0, 0);
    for &fh in &m.shell(m.solid(s).outer).faces {
        for he in &m.face(fh).outer.half_edges {
            if seen.insert(he.edge) {
                match m.edge_curve(he.edge) {
                    Curve::Circle(_) => circles += 1,
                    Curve::Line(_) => lines += 1,
                }
            }
        }
    }
    (circles, lines)
}

/// ★★ **A moved cylinder keeps its edges' kinds** — asked of two populations whose carriers do
/// not share one chain after the move: a cylinder on a tilted frame (a frame node) moved by a
/// non-dyadic translation, recorded behind that history, and a `+z` cylinder turned 30° about `z`
/// (its seam reference goes irrational, so the turn is recorded) then moved along its axis — the
/// lateral records the move behind its turn, the caps the turn fixed carry it. The edge kinds are
/// relations between two directions, which a rigid motion keeps; a move must neither refuse these
/// nor change them.
#[test]
fn a_moved_cylinder_keeps_its_edges_kinds() {
    use nacre_exact::{Isometry, Rat};
    let third = Rat::new(1, 3).expect("a third");
    let zero = Rat::from_int(0);
    let cases: [(&str, [f64; 3], Vec<Isometry>); 2] = [
        (
            "tilted frame, moved",
            [0.3, -0.4, 1.0],
            vec![Isometry::translation([third, zero, third])],
        ),
        (
            "turned, moved along the axis",
            [0.0, 0.0, 1.0],
            vec![
                crate::fixtures::rot_iso(nacre_exact::Axis::Z, 30),
                Isometry::translation([zero, zero, third]),
            ],
        ),
    ];
    for (name, axis, moves) in cases {
        let mut m = Model::new();
        let mut s = cylinder(
            &mut m,
            Point3::from_array([1.0, 2.0, 0.5]),
            Vector3::from_array(axis),
            0.7,
            3.0,
        )
        .solid;
        m.rebuild_adjacency();
        let v0 = nacre_props::mass_props(&m, s).unwrap().volume;
        for iso in moves {
            s = crate::fixtures::xf(&mut m, s, iso);
        }
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{name}: {:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(edge_kinds(&m, s), (2, 1), "{name}: two rims and a seam");
        let v = nacre_props::mass_props(&m, s).unwrap().volume;
        assert!((v - v0).abs() < 1e-9, "{name}: {v0} → {v}");
    }
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    #[test]
    fn prop_cylinder_structural(
        base in prop::array::uniform3(-1e3f64..1e3),
        axis in prop::array::uniform3(-1.0f64..1.0),
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        let axis = Vector3::from_array(axis);
        prop_assume!(axis.norm() > 0.1); // skip near-zero axes
        let mut m = Model::new();
        cylinder(&mut m, Point3::from_array(base), axis, r, h);
        m.rebuild_adjacency();
        prop_assert_eq!(m.vertex_count(), 2);
        prop_assert_eq!(m.edge_count(), 3);
        prop_assert_eq!(m.face_count(), 3);
        prop_assert!(every_edge_used_twice_opposite(&m).is_ok());
    }

    /// The same structure over the population that **actually stressed the exact
    /// arithmetic**: an axis a hair off `ẑ`, whose tiny components carry a full f64's worth
    /// of decimal digits and so lift to rationals with ~10²⁰ denominators — squaring one of
    /// those leaves `i128`. The sketch plane states itself from basis crosses of the written
    /// decimals (`SketchPlane::from_origin_normal`), which is what keeps this family built.
    ///
    /// ★ The sibling above cannot stand in for this: its uniform axis reaches this family
    /// about **0.01%** of the time (measured). Here it is ~80%.
    #[test]
    fn prop_cylinder_near_axis_aligned_is_built_not_refused(
        u in -1.0f64..1.0,
        v in -1.0f64..1.0,
        k in 1i32..9,
        j in 1i32..9,
        r in 0.5f64..10.0,
        h in 0.1f64..10.0,
    ) {
        let axis = Vector3::from_array([u * 10f64.powi(-k), v * 10f64.powi(-j), 1.0]);
        let mut m = Model::new();
        cylinder(&mut m, Point3::origin(), axis, r, h);
        m.rebuild_adjacency();
        prop_assert_eq!(m.vertex_count(), 2);
        prop_assert_eq!(m.edge_count(), 3);
        prop_assert_eq!(m.face_count(), 3);
        prop_assert!(every_edge_used_twice_opposite(&m).is_ok());
    }
}
