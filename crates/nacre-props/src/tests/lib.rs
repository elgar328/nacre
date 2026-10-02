use super::*;
use nacre_exact::Axis;
use nacre_math::{Point2, Vector3};
use nacre_ops::SketchFrame;
use nacre_ops::{OpOutput, Operation, Profile2d, apply};

/// Relative-or-absolute comparison sized for accumulated f64 error.
fn close(got: f64, expected: f64) -> bool {
    (got - expected).abs() <= 1e-9 * expected.abs().max(1.0)
}

// --- golden ---

#[test]
fn cube_mass_matches_analytic() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let props = mass_props(&m, s).unwrap();
    assert!(close(props.volume, 24.0), "volume {}", props.volume); // 2·3·4
    assert!(close(props.area, 52.0), "area {}", props.area); // 2(6+8+12)
}

#[test]
fn cylinder_mass_matches_analytic() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    )
    .solid;
    let props = mass_props(&m, s).unwrap();
    assert!(close(props.volume, 20.0 * PI), "volume {}", props.volume); // πr²h
    assert!(close(props.area, 28.0 * PI), "area {}", props.area); // 2πr² + 2πrh
}

/// ★ **Two doors, one fact.** A planar face's normal is answered twice — by
/// [`face_props`] for the whole face and by [`face_normal_at`] at a point — and the two
/// must say the same thing, or a caller would get a different answer depending on which
/// it happened to ask. (The point plays no part on a plane, which this also shows.)
#[test]
fn the_two_normal_doors_agree_on_a_plane() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    m.rebuild_adjacency();
    let shell = m.solid(s).outer;
    for &fh in &m.shell(shell).faces {
        let whole = face_props(&m, fh)
            .unwrap()
            .normal
            .expect("a box face is planar");
        // Two different points of the same face — its own centroid and a corner.
        let corner = m.vertex_point(he_start(&m, m.face(fh).outer.half_edges[0]).unwrap());
        for p in [face_props(&m, fh).unwrap().centroid, corner] {
            let at = face_normal_at(&m, fh, p).expect("a planar face has a normal anywhere");
            assert!((at - whole).norm() < 1e-12, "{at:?} vs {whole:?}");
        }
    }
}

/// The cylinder branch of the same door, where the point is the whole question: a free
/// cylinder's wall faces away from its axis, everywhere on it.
#[test]
fn a_free_cylinders_wall_faces_away_from_its_axis() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    )
    .solid;
    m.rebuild_adjacency();
    let shell = m.solid(s).outer;
    let wall = *m
        .shell(shell)
        .faces
        .iter()
        .find(|&&fh| matches!(m.surface_cache(m.face(fh).surface), Surface::Cylinder(_)))
        .expect("a lateral face");
    assert!(
        face_props(&m, wall).unwrap().normal.is_none(),
        "a curved face has no single normal — that is why the point door exists"
    );
    let cyl = match m.surface_cache(m.face(wall).surface) {
        Surface::Cylinder(c) => *c,
        _ => unreachable!(),
    };
    for k in 0..8 {
        let u = std::f64::consts::TAU * f64::from(k) / 8.0;
        let p = cyl.point_at(u, 2.5);
        let n = face_normal_at(&m, wall, p).expect("off the axis");
        let radial = p - cyl.axis().origin();
        let radial = radial - cyl.axis().direction() * radial.dot(cyl.axis().direction());
        assert!(
            n.dot(radial) > 0.0,
            "u={u}: the wall faces away from the axis"
        );
        assert!((n.norm() - 1.0).abs() < 1e-12, "unit");
    }
}

#[test]
fn concave_extrude_mass() {
    // An L-shaped profile: a 2×2 square with the top-right 1×1 corner removed.
    // Concave, so the area-weighted centroid differs from the vertex average —
    // a vertex-average bug in the flux would fail the volume assertion.
    let pts = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let profile = Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap();
    let dist = 3.0;
    let mut m = Model::new();
    let op = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile,
        dist,
    };
    let OpOutput::Extrude { solid: s, .. } = apply(&mut m, &op).unwrap() else {
        unreachable!("extrude yields Extrude output");
    };

    // Independent expected base area via the 2D shoelace (not mass_props).
    let a_l = shoelace(&pts);
    assert!((a_l - 3.0).abs() < 1e-12); // 4 − 1

    let props = mass_props(&m, s).unwrap();
    assert!(close(props.volume, a_l * dist), "volume {}", props.volume);
}

// --- bounds / face_props / centroid ---

/// **The trap this exists for.** A cylinder's lateral face bulges past its two
/// seam vertices, so a hull of the vertices reports a box that is too small in
/// the radial directions — a wrong answer with nothing to tip the caller off.
/// The axis is deliberately oblique so the analytic term is not axis-aligned, and the seam is
/// put on `(0, −0.8, 0.6)`, square to `x`, so the vertices sit at `x = 0` while the barrel
/// spans `x ∈ [−2, 2]`.
#[test]
fn bounds_follow_the_curve_not_the_vertices() {
    let axis = Vector3::from_array([0.0, 3.0, 4.0]); // unit (0, .6, .8)
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        axis,
        Vector3::from_array([0.0, -0.8, 0.6]),
        2.0,
        5.0,
    )
    .solid;
    let (lo, hi) = bounds(&m, s).unwrap();

    // x is fully perpendicular to the axis, so the barrel spans the diameter.
    assert!(close(lo[0], -2.0) && close(hi[0], 2.0), "{lo:?} {hi:?}");
    // Along y and z the circles foreshorten by √(1 − (n̂·e)²).
    let u = axis.normalize().unwrap().as_array();
    let end = Point3::origin() + axis.normalize().unwrap() * 5.0;
    for i in [1, 2] {
        let pad = 2.0 * (1.0 - u[i] * u[i]).sqrt();
        assert!(close(lo[i], -pad), "axis {i} lo {}", lo[i]);
        assert!(
            close(hi[i], end.as_array()[i] + pad),
            "axis {i} hi {}",
            hi[i]
        );
    }

    // And the vertex hull really is smaller — otherwise this test proves nothing.
    let mut vlo = [f64::INFINITY; 3];
    let mut vhi = [f64::NEG_INFINITY; 3];
    let mut i = 0u32;
    while let Some(vh) = m.vertex_handle_at(i) {
        i += 1;
        let p = m.vertex_point(vh).as_array();
        for i in 0..3 {
            vlo[i] = vlo[i].min(p[i]);
            vhi[i] = vhi[i].max(p[i]);
        }
    }
    assert!(
        vlo[0] > lo[0] + 1.0,
        "the vertex hull must be visibly wrong here: {vlo:?} vs {lo:?}"
    );
}

#[test]
fn bounds_of_a_box_are_its_corners() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, 2.0, 0.5]),
        Point3::from_array([3.0, 4.0, 9.0]),
    );
    let (lo, hi) = bounds(&m, s).unwrap();
    assert!(
        close(lo[0], -1.0) && close(lo[1], 2.0) && close(lo[2], 0.5),
        "{lo:?}"
    );
    assert!(
        close(hi[0], 3.0) && close(hi[1], 4.0) && close(hi[2], 9.0),
        "{hi:?}"
    );
}

/// The query a script needs to *name* a face: filter by normal, then by position.
#[test]
fn faces_can_be_picked_by_their_geometry() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::origin(),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let top = m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&f| (f, face_props(&m, f).unwrap()))
        .filter(|(_, p)| p.normal.is_some_and(|n| n.dot(up) > 0.5))
        .max_by(|a, b| a.1.centroid[2].partial_cmp(&b.1.centroid[2]).unwrap())
        .unwrap();
    assert!(close(top.1.area, 6.0), "area {}", top.1.area); // 2·3
    assert!(close(top.1.centroid[2], 4.0), "z {}", top.1.centroid[2]);
    assert!(close(top.1.centroid[0], 1.0) && close(top.1.centroid[1], 1.5));
}

/// A hole moves the face's centroid; using the outer ring alone would not.
#[test]
fn a_face_with_an_off_centre_hole_reports_the_region() {
    let sq = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    let mut m = Model::new();
    let op = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: Profile2d::with_holes(sq(0.0, 10.0), vec![sq(1.0, 3.0)]).unwrap(),
        dist: 1.0,
    };
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        unreachable!()
    };
    // faces[1] is the top cap — the one carrying the inner loop.
    let p = face_props(&m, faces[1]).unwrap();
    assert!(close(p.area, 96.0), "area {}", p.area); // 100 − 4
    let want = (100.0 * 5.0 - 4.0 * 2.0) / 96.0;
    assert!(
        close(p.centroid[0], want) && close(p.centroid[1], want),
        "{:?}",
        p.centroid
    );
}

/// Asymmetric on purpose: a symmetric solid puts the centroid at the middle
/// whatever the arithmetic does, so it proves nothing.
#[test]
fn centroid_of_an_l_prism_is_not_its_bounding_centre() {
    let pts = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let mut m = Model::new();
    let op = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap(),
        dist: 3.0,
    };
    let OpOutput::Extrude { solid: s, .. } = apply(&mut m, &op).unwrap() else {
        unreachable!()
    };
    let c = centroid(&m, s).unwrap();
    // Area 3 = unit squares at (.5,.5), (1.5,.5), (.5,1.5) ⇒ centroid (5/6, 5/6).
    assert!(close(c[0], 5.0 / 6.0) && close(c[1], 5.0 / 6.0), "{c:?}");
    assert!(close(c[2], 1.5), "{c:?}");
}

/// A cavity carries the opposite sign, so an off-centre void must push the
/// centroid away from it. Dropping that sign is invisible on a centred void.
#[test]
fn an_off_centre_void_pushes_the_centroid_away() {
    let mut m = Model::new();
    let outer =
        nacre_ops::fixtures::cuboid(&mut m, Point3::origin(), Point3::from_array([10.0; 3]));
    let inner = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0; 3]),
        Point3::from_array([3.0, 3.0, 3.0]),
    );
    let r = nacre_ops::boolean(&mut m, nacre_ops::BoolKind::Cut, outer, inner).unwrap();
    m.rebuild_adjacency();
    assert_eq!(r.len(), 1);
    let c = centroid(&m, r[0]).unwrap();
    // (1000·5 − 8·2)/992 on every axis: the void at (2,2,2) drags the rest up.
    let want = (1000.0 * 5.0 - 8.0 * 2.0) / 992.0;
    assert!(
        want > 5.0,
        "the void must be off-centre for this to prove anything"
    );
    for i in 0..3 {
        assert!(close(c[i], want), "axis {i}: {} vs {want}", c[i]);
    }
}

#[test]
fn a_curved_solid_refuses_a_centroid_but_still_reports_volume() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cylinder_with_seam(
        &mut m,
        Point3::origin(),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        5.0,
    )
    .solid;
    assert!(matches!(
        centroid(&m, s),
        Err(PropsError::CentroidOfCurvedFace)
    ));
    assert!(close(mass_props(&m, s).unwrap().volume, 20.0 * PI));
}

/// A `2·hw` square centred on a `size` cube's lid.
///
/// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
/// onto the face's plane and the axes are `u = +x_hat`, `v = +y_hat`, so a frame point `(a, b)`
/// is world `(a, b, size)`.
fn centred_on_the_lid(size: f64, hw: f64) -> Profile2d {
    let (cx, cy) = (0.5 * size, 0.5 * size);
    Profile2d::polygon(
        [
            [cx - hw, cy + hw],
            [cx - hw, cy - hw],
            [cx + hw, cy - hw],
            [cx + hw, cy + hw],
        ]
        .iter()
        .map(|&p| Point2::from_array(p))
        .collect(),
    )
    .unwrap()
}

/// Pad a `2·hw` square boss of height `dist` on a `size` cube's top face,
/// returning the padded solid's mass.
fn cube_then_pad(size: f64, hw: f64, dist: f64) -> MassProps {
    let sq = |s: f64| {
        Profile2d::polygon(
            [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut m = Model::new();
    let __w1 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w1,
            profile: sq(size),
            dist: size,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let boss = centred_on_the_lid(size, hw);
    let solid = nacre_ops::fixtures::pad(&mut m, faces[1], boss, dist)
        .unwrap()
        .solid();
    mass_props(&m, solid).unwrap()
}

#[test]
fn pad_boss_mass() {
    // Unit cube + a 0.4-square boss (area 0.16, perimeter 1.6) of height 0.5.
    // Volume = 1 + 0.16·0.5 = 1.08; area = 6 + 1.6·0.5 = 6.8. This is the
    // *absolute* test of the inner-loop subtraction (the hole is not filled).
    let m = cube_then_pad(1.0, 0.2, 0.5);
    assert!(close(m.volume, 1.08), "vol {}", m.volume);
    assert!(close(m.area, 6.8), "area {}", m.area);
}

/// Carve a `2·hw` square pocket of depth `dist` into a `size` cube's top face.
fn cube_then_pocket(size: f64, hw: f64, dist: f64) -> MassProps {
    let sq = |s: f64| {
        Profile2d::polygon(
            [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| Point2::from_array(p))
                .collect(),
        )
        .unwrap()
    };
    let mut m = Model::new();
    let __w0 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w0,
            profile: sq(size),
            dist: size,
        },
    )
    .unwrap() else {
        unreachable!()
    };
    let pocket = centred_on_the_lid(size, hw);
    let solid = nacre_ops::fixtures::pocket(&mut m, faces[1], pocket, dist)
        .unwrap()
        .solid();
    mass_props(&m, solid).unwrap()
}

#[test]
fn pocket_mass() {
    // Unit cube − a 0.4-square pocket of depth 0.5. Volume = 1 − 0.16·0.5 =
    // 0.92 (the inward walls contribute negatively); area = 6 + 1.6·0.5 = 6.8
    // (same as the boss). Exercises the hole subtraction with inward walls.
    let m = cube_then_pocket(1.0, 0.2, 0.5);
    assert!(close(m.volume, 0.92), "vol {}", m.volume);
    assert!(close(m.area, 6.8), "area {}", m.area);
}

/// Build a hollow solid: an `outer`-cube with a concentric `inner`-cube void,
/// the inner cube's shell reversed inward (M5 containment). Returns its mass.
fn cube_in_cube(min: Point3, outer: f64, inner: f64) -> MassProps {
    let mut m = Model::new();
    let ext = |s: f64| Vector3::from_array([s, s, s]);
    // Floor and height: a corner plus a random size differs from that corner by no decimal an
    // f64 need carry.
    let on = |m: &mut Model, lo: Point3, size: f64| {
        let [x, y, z] = lo.as_array();
        nacre_ops::fixtures::cuboid_on(m, [x, y], [x + size, y + size], z, size)
    };
    let a = on(&mut m, min, outer);
    let gap = 0.5 * (outer - inner); // centered ⇒ strictly interior on all sides
    let inner_min = min + ext(gap);
    let b = on(&mut m, inner_min, inner);

    let b_outer = m.solid(b).outer;
    let void = m.reversed_shell(b_outer);
    let a_outer = m.solid(a).outer;
    let hollow = m.push_solid(Solid {
        outer: a_outer,
        cavities: vec![void],
    });
    m.restore_live(vec![hollow]); // supersede the two source cubes
    mass_props(&m, hollow).unwrap()
}

#[test]
fn cube_in_cube_subtracts_the_void() {
    // 4-cube with a concentric 2-cube void: V = 4³ − 2³ = 56; total surface
    // = outer 6·4² + void 6·2² = 96 + 24 = 120 (both surfaces bound material).
    let props = cube_in_cube(Point3::origin(), 4.0, 2.0);
    assert!(close(props.volume, 56.0), "vol {}", props.volume);
    assert!(close(props.area, 120.0), "area {}", props.area);
}

/// Unsigned area of a 2D polygon (independent check for the concave test).
fn shoelace(pts: &[[f64; 2]]) -> f64 {
    let mut two_area = 0.0;
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        two_area += a[0] * b[1] - b[0] * a[1];
    }
    0.5 * two_area.abs()
}

// --- proptest ---

use proptest::prelude::*;

fn wide_point() -> impl Strategy<Value = Point3> {
    prop::array::uniform3(-1e6f64..1e6).prop_map(Point3::from_array)
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Random box, possibly far from the origin — exercises the R
    /// cancellation guard as well as the polygon path.
    #[test]
    fn prop_cuboid_mass(
        min in wide_point(),
        a in 0.1f64..1e3, b in 0.1f64..1e3, c in 0.1f64..1e3,
    ) {
        // Floor and height, as drawn: two random corners differ by no decimal an f64 carries.
        let [x, y, z] = min.as_array();
        let mut m = Model::new();
        let s = nacre_ops::fixtures::cuboid_on(&mut m, [x, y], [x + a, y + b], z, c);
        let props = mass_props(&m, s).unwrap();
        prop_assert!(close(props.volume, a * b * c), "vol {} vs {}", props.volume, a*b*c);
        prop_assert!(close(props.area, 2.0 * (a*b + b*c + c*a)), "area {}", props.area);
    }

    /// Random cylinder with an arbitrary axis direction — exercises the
    /// non-axis-aligned frame, the lateral sign, and the 2π band height.
    #[test]
    fn prop_cylinder_mass(
        base in wide_point(),
        dir in prop::array::uniform3(-1e3f64..1e3)
            .prop_map(Vector3::from_array)
            .prop_filter("axis must be clearly nonzero", |v| v.norm() > 1e-3),
        r in 0.5f64..10.0, h in 0.5f64..1e3,
    ) {
        let mut m = Model::new();
        let s = nacre_ops::fixtures::cylinder(&mut m, base, dir, r, h).solid;
        let props = mass_props(&m, s).unwrap();
        prop_assert!(close(props.volume, PI * r * r * h), "vol {} vs {}", props.volume, PI*r*r*h);
        prop_assert!(close(props.area, 2.0 * PI * r * (r + h)), "area {}", props.area);
    }

    /// A boss adds `A_p·dist` of volume and `P·dist` of surface (the top hole
    /// area cancels the boss cap). For a `2·hw` square: `A_p = 4hw²`, `P = 8hw`.
    /// This absolutely exercises the inner-loop subtraction (uncancelled hole).
    #[test]
    fn prop_pad_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..5.0) {
        let m = cube_then_pad(size, hw, dist);
        let vol = size * size * size + 4.0 * hw * hw * dist;
        let area = 6.0 * size * size + 8.0 * hw * dist;
        prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
        prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
    }

    /// A pocket *removes* `A_p·dist` of volume (inward walls) while adding the
    /// same `P·dist` of surface as a boss. `dist ≤ 0.9 < size` keeps it blind.
    #[test]
    fn prop_pocket_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..0.9) {
        let m = cube_then_pocket(size, hw, dist);
        let vol = size * size * size - 4.0 * hw * hw * dist;
        let area = 6.0 * size * size + 8.0 * hw * dist;
        prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
        prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
    }

    /// A concentric cube void of any interior size: V = outer³ − inner³, and
    /// the total surface is the sum of both. `min` ranges far from the origin
    /// to exercise the R-cancellation guard over the *full* (outer + void)
    /// boundary, not just the outer shell.
    #[test]
    fn prop_cube_in_cube(
        min in wide_point(),
        outer in 2.0f64..1e3,
        ratio in 0.1f64..0.85,
    ) {
        let inner = outer * ratio;
        let props = cube_in_cube(min, outer, inner);
        let vol = outer * outer * outer - inner * inner * inner;
        let area = 6.0 * (outer * outer + inner * inner);
        prop_assert!(close(props.volume, vol), "vol {} vs {}", props.volume, vol);
        prop_assert!(close(props.area, area), "area {} vs {}", props.area, area);
    }
}
