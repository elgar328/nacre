//! Tests for the boolean assembly — what it builds from the engine's draft faces.

/// **The half-height boss's lateral is one face** — a band with one whole rim and one
/// chain rim (the boss's cap sits inside the plate, so the lateral ends on a chain of arcs and
/// rulings). The region emitter builds it whole, and this locks the model.
mod region_tests {
    use crate::BoolKind;
    use nacre_math::{Point3, Vector3};
    use nacre_topo::Model;

    #[test]
    fn the_half_height_boss_lateral_is_one_face() {
        let mut m = Model::new();
        let plate = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([40.0, 20.0, -10.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            5.0,
            20.0,
        )
        .solid;
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the half-height boss builds");
        m.rebuild_adjacency();
        let laterals: Vec<_> = m
            .reachable()
            .faces
            .iter()
            .copied()
            .filter(|&f| {
                matches!(
                    m.surface_cache(m.face(f).surface),
                    nacre_geom::Surface::Cylinder(_)
                )
            })
            .collect();
        assert_eq!(
            laterals.len(),
            1,
            "one lateral face for the half-height boss"
        );
        let face = m.face(laterals[0]);
        assert_eq!(face.inner.len(), 1, "two loops: one rim each");
        let (mut closed, mut arcs) = (0usize, 0usize);
        for he in std::iter::once(&face.outer)
            .chain(&face.inner)
            .flat_map(|l| &l.half_edges)
        {
            let e = m.edge(he.edge);
            assert_ne!(e.surfaces[0], e.surfaces[1], "no seam edge");
            let circular = matches!(m.edge_curve(he.edge), nacre_geom::Curve::Circle(_));
            match (circular, e.vertices[0] == e.vertices[1]) {
                (true, true) => closed += 1,
                (true, false) => arcs += 1,
                _ => {}
            }
        }
        assert_eq!(closed, 1, "one whole rim");
        assert!(arcs >= 1, "a chain rim of arcs");
    }
}

/// **The result guard's carrier half** (`check_result_topology`, `EdgeCarrierMismatch`): a closed
/// cube, every count clean, but one edge stated on a pair its two faces do not lie on — Bottom ·
/// Front restated as Bottom · Right. No boolean in the suite produces this, so the guard is
/// planted: the untouched cube passes, the restated one is refused with the edge's end.
#[test]
fn an_edge_stated_on_a_pair_its_faces_do_not_keep_is_refused() {
    use crate::RejectReason;
    use nacre_math::Point3;
    use nacre_topo::{HalfEdge, Model, Shell, Solid};
    let mut m = Model::new();
    let cube = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    assert_eq!(
        super::check_result_topology(&m, &[cube]),
        None,
        "the cube itself"
    );
    let faces = m.shell(m.solid(cube).outer).faces.clone();
    // The extrude's face order: the floor, the top, then the walls from the rectangle's first
    // edge — `y = 0` first, `x = 0` last. The last wall's plane does not hold the floor–front edge.
    let (bottom, front, side) = (faces[0], faces[2], faces[5]);
    let edges_of = |m: &Model, f| -> Vec<_> {
        m.face(f)
            .outer
            .half_edges
            .iter()
            .map(|he| he.edge)
            .collect()
    };
    let shared = *edges_of(&m, bottom)
        .iter()
        .find(|e| edges_of(&m, front).contains(e))
        .expect("Bottom and Front share an edge");
    let wrong = m
        .push_edge(
            [m.face(bottom).surface, m.face(side).surface],
            m.edge(shared).vertices,
            |_| nacre_topo::EdgeGiven::NONE,
        )
        .expect("a line of two planes");
    let restated = |m: &mut Model, f| {
        let mut face = m.face(f).clone();
        for he in &mut face.outer.half_edges {
            if he.edge == shared {
                *he = HalfEdge {
                    edge: wrong,
                    forward: he.forward,
                };
            }
        }
        m.push_face(face)
    };
    let (b2, f2) = (restated(&mut m, bottom), restated(&mut m, front));
    let shell = m.push_shell(Shell {
        faces: faces
            .iter()
            .map(|&f| match f {
                f if f == bottom => b2,
                f if f == front => f2,
                f => f,
            })
            .collect(),
    });
    let solid = m.push_solid(Solid {
        outer: shell,
        cavities: vec![],
    });
    let got = super::check_result_topology(&m, &[solid]);
    assert!(
        matches!(got, Some((RejectReason::EdgeCarrierMismatch, Some(_)))),
        "{got:?}"
    );
}

/// The two pierce nodes where plane classes 4 and 5 meet plane class 3 on cylinder class 0's
/// lateral — a rim `(0, 3)` cut at two points, and a face of either kind walking it as two arcs.
fn two_arc_rim() -> (
    crate::combinatorics::NodeId,
    crate::combinatorics::NodeId,
    impl Fn(crate::planes::ClassIx, crate::combinatorics::Wall) -> crate::draft::LocalFace,
) {
    use crate::combinatorics::NodeId;
    use crate::draft::{Bound, LocalFace, Ring};
    use nacre_topo::QuadRoot;
    let n1 = NodeId::pierce(3, 4, 0, QuadRoot::Lo);
    let n2 = NodeId::pierce(3, 5, 0, QuadRoot::Lo);
    let face = move |surf, wall| LocalFace {
        surf,
        outer: Bound::Ring(Ring::new(vec![n1, n2], vec![wall, wall])),
        inner: Vec::new(),
        flip: false,
    };
    (n1, n2, face)
}

/// **A ring that dissolves to its whole circle closes in its face's own spelling**
/// (`dissolve_straight_angles`'s `closed`). The rim `(cylinder 0, plane 3)` is one edge: the cap
/// bounds it as `Circle { cyl: 0 }`, the lateral as `Rim { plane: 3, ccw }` — the cap's spelling
/// on a lateral names no plane and the assembly cannot build it. Both faces together is the
/// per-solid pass's real call; each alone separates "the face kind decides" from "a cap is there
/// too". An arc wall of another cylinder on the lateral closes nothing: the ring stays as its
/// producer wrote it.
#[test]
fn a_dissolved_rim_closes_in_its_faces_spelling() {
    use crate::combinatorics::Wall;
    use crate::draft::{Bound, LocalFace};
    use crate::planes::ClassIx;
    let (n1, n2, face) = two_arc_rim();
    let arc = |cyl, ccw| Wall::Arc { cyl, ccw, plane: 3 };
    let lateral = face(ClassIx::Cyl(0), arc(0, true));
    let cap = face(ClassIx::Plane(3), arc(0, false));
    let run = |mut v: Vec<LocalFace>| -> Vec<LocalFace> {
        let which: Vec<usize> = (0..v.len()).collect();
        super::coplanar::dissolve_straight_angles(&mut v, &which);
        v
    };
    let both = run(vec![lateral.clone(), cap.clone()]);
    assert!(
        matches!(
            both[0].outer,
            Bound::Rim {
                plane: 3,
                ccw: true
            }
        ),
        "the lateral closes as its rim: {:?}",
        both[0].outer
    );
    assert!(
        matches!(both[1].outer, Bound::Circle { cyl: 0 }),
        "the cap closes as its circle: {:?}",
        both[1].outer
    );
    let alone = run(vec![lateral]);
    assert!(
        matches!(
            alone[0].outer,
            Bound::Rim {
                plane: 3,
                ccw: true
            }
        ),
        "{:?}",
        alone[0].outer
    );
    let alone = run(vec![cap]);
    assert!(
        matches!(alone[0].outer, Bound::Circle { cyl: 0 }),
        "{:?}",
        alone[0].outer
    );
    let foreign = run(vec![face(ClassIx::Cyl(0), arc(1, true))]);
    assert!(
        matches!(&foreign[0].outer, Bound::Ring(r) if r.nodes == vec![n1, n2]),
        "another cylinder's arcs keep their nodes: {:?}",
        foreign[0].outer
    );
}

/// **A rim is cut for the solid that holds its node, and only for it**
/// (`HeldRims::cut_per_group`). Two result solids meet the rim `(0, 3)`: solid 0's face on plane 3
/// still holds both nodes, solid 1's per-solid pass dissolved them and closed its cap to the
/// circle. The table the chart reads, blind to solids, calls the rim cut; the assembly asks per
/// solid, so solid 1's whole circle is not refused for solid 0's corners.
#[test]
fn a_rim_is_cut_only_for_the_solid_that_holds_its_node() {
    use crate::combinatorics::Wall;
    use crate::draft::{Bound, CutRim, LocalFace, held_rims};
    use crate::planes::ClassIx;
    let (n1, n2, face) = two_arc_rim();
    let holds = face(ClassIx::Plane(3), Wall::Plane(4));
    let whole = LocalFace {
        surf: ClassIx::Plane(3),
        outer: Bound::Circle { cyl: 0 },
        inner: Vec::new(),
        flip: false,
    };
    let faces = vec![holds, whole];
    let split = std::collections::HashMap::from([(
        (0, 3),
        CutRim {
            nodes: vec![n1, n2],
        },
    )]);
    let held = held_rims(&faces, &split);
    assert!(
        held.get(&(0, 3)).is_some(),
        "premise: blind to solids, the rim is cut"
    );
    let cut = held.cut_per_group(&faces, &[0, 1]);
    assert_eq!(
        cut,
        std::collections::HashSet::from([(0, 0, 3)]),
        "cut for solid 0 only"
    );
}

/// **An arch standing on a cylinder's cap, touching its rim at two corners, fused** — two
/// result solids, and the cylinder's per-solid straight-angle pass dissolves both rim nodes:
/// for the cylinder they are straight runs of its rim, for the arch they are corners. The
/// cylinder's cap closes as its circle and its lateral as its rim, and the assembly asks per
/// solid whether that rim is cut, so the cylinder builds. Asked of the whole result, the arch's
/// corners on the cap's plane make the rim cut for the cylinder too, which
/// `CylinderStagesDisagree` refuses; the plant of a group-blind check measures that.
///
/// ★ **The arch's corner is named by the arch's planes.** The rim node is a pierce point — two
/// planes and the cylinder — and the arch has no face on the cylinder, so the naming takes the
/// canonical triple of the arch's own three planes there (`x = 7`, `y = 4`, `z = 0`); kept, the
/// pierce name would name a surface the arch lacks and the shipped fence would refuse the vertex
/// (`VertexNamesAbsentSurface`).
///
/// Two solids, `validate`-clean, every vertex named by its own solid's faces, volumes `1040` and
/// `250π`.
#[test]
fn an_arch_on_a_cylinders_rim_fuses_to_two_solids() {
    use crate::{BoolKind, OpOutput, Operation, Profile2d, SketchFrame, apply};
    use nacre_exact::Axis;
    use nacre_math::{Point2, Point3, Vector3};
    use nacre_topo::Model;
    let mut m = Model::new();
    // Axis `x = 4, y = 0`, radius 5, `z ∈ [−10, 0]`: the top rim passes through `(7, 4)` and
    // `(1, 4)` (`3² + 4² = 5²`).
    let cylinder = crate::fixtures::cylinder(
        &mut m,
        Point3::from_array([4.0, 0.0, -10.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        5.0,
        10.0,
    )
    .solid;
    // `z ∈ [0, 10]` over an arch whose two inner corners are those rim points. Near each, the
    // arch is a quadrant pointing away from the disk (`x ≥ 7, y ≥ 4` and `x ≤ 1, y ≥ 4`), and the
    // notch `[1, 7] × [4, 8]` keeps it off the disk everywhere else: the two solids touch at the
    // two corners only.
    let arch = Profile2d::polygon(
        [
            [7.0, 4.0],
            [12.0, 4.0],
            [12.0, 12.0],
            [-4.0, 12.0],
            [-4.0, 4.0],
            [1.0, 4.0],
            [1.0, 8.0],
            [7.0, 8.0],
        ]
        .iter()
        .map(|&[x, y]| Point2::from_array([x, y]))
        .collect(),
    )
    .expect("the arch");
    let op = Operation::Extrude {
        frame: SketchFrame::world(&m, Axis::Z),
        profile: arch,
        dist: 10.0,
    };
    let Ok(OpOutput::Extrude { solid: arch, .. }) = apply(&mut m, &op) else {
        panic!("the arch extrudes")
    };
    m.rebuild_adjacency();
    let out = crate::boolean(&mut m, BoolKind::Fuse, cylinder, arch)
        .unwrap_or_else(|e| panic!("the arch and the cylinder fuse: {e:?}"));
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2, "two solids touching at two points");
    assert_eq!(nacre_validate::validate(&m), vec![], "a valid model");
    let mut volumes: Vec<f64> = out
        .iter()
        .map(|&s| {
            assert_eq!(
                crate::transform::foreign_named_vertex(&m, s),
                None,
                "every vertex named by its own solid's faces"
            );
            nacre_props::mass_props(&m, s).expect("mass props").volume
        })
        .collect();
    volumes.sort_by(f64::total_cmp);
    let want = [250.0 * std::f64::consts::PI, 1040.0];
    for (got, want) in volumes.iter().zip(want) {
        assert!((got - want).abs() < 1e-9 * want, "{got} vs {want}");
    }
}
