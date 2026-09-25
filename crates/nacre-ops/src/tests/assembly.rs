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
