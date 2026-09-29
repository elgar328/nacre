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
        let plate = m.add_cuboid(
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
        let outer = &m.face(laterals[0]).outer;
        let (mut closed, mut arcs) = (0usize, 0usize);
        for he in &outer.half_edges {
            let circular = matches!(m.edge_curve(he.edge), nacre_geom::Curve::Circle(_));
            let [a, b] = m.edge(he.edge).vertices;
            match (circular, a == b) {
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
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    assert_eq!(
        super::check_result_topology(&m, &[cube]),
        None,
        "the cube itself"
    );
    let faces = m.shell(m.solid(cube).outer).faces.clone();
    // `add_cuboid`'s face order: Bottom, Top, Front, Back, Left, Right.
    let (bottom, front, right) = (faces[0], faces[2], faces[5]);
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
            [m.face(bottom).surface, m.face(right).surface],
            m.edge(shared).vertices,
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
