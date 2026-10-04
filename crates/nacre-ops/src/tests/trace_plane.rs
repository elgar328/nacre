//! Unit tests of the planar trace's own ring reader (`plane_ring`, `RingFail`).

use super::*;

/// ★★★ **Where a named curved ring stops, and under which name.**
///
/// ★★★★★ **The vessel carries a curved ring through, and only a *collapsed* name is refused.**
///
/// Refusing a pierce corner or a curved carrier would be the type running out — plane ids can
/// hold neither — rather than a decision; the ring holds names. What stops here is a name that is
/// degenerate *as a name* — two
/// of its three planes equal — which is a fact about the triple and not about cylinders.
#[test]
fn a_named_curved_ring_rides_through_and_only_a_collapsed_name_stops() {
    let three =
        |a, b, c| combinatorics::NodeId::three_planes(combinatorics::Canon3::three([a, b, c]));
    let pierce = combinatorics::NodeId::pierce(0, 1, 0, nacre_topo::QuadRoot::Lo);
    let plane = |c| crate::combinatorics::Wall::Plane(c);
    let ruling = crate::combinatorics::Wall::Ruling {
        cyl: 0,
        side: 1,
        up: true,
    };
    // A curved carrier rides through, carried as itself.
    let curved_carrier = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
        walls: vec![plane(1), ruling, plane(3)],
        concurrencies: vec![],
    };
    let (_, walls) = plane_ring(&curved_carrier).expect("a curved carrier is describable");
    assert_eq!(
        walls[1], ruling,
        "the carrier is the producer's, unflattened"
    );
    // So does a pierce corner, as its own name.
    let curved_corner = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), pierce, three(0, 1, 3)],
        walls: vec![plane(1), ruling, plane(3)],
        concurrencies: vec![],
    };
    let (ts, _) = plane_ring(&curved_corner).expect("a pierce corner is describable");
    assert_eq!(ts[1], pierce, "the corner is the producer's, unflattened");
    // ★ What is still refused, and for a reason that has nothing to do with cylinders.
    let collapsed = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 0, 3), three(0, 1, 3)],
        walls: vec![plane(1), plane(2), plane(3)],
        concurrencies: vec![],
    };
    assert!(matches!(plane_ring(&collapsed), Err(RingFail::Collapsed)));
    // And the plane-only ring still comes back with both halves.
    let plain = combinatorics::NamedRing {
        triples: vec![three(0, 1, 2), three(0, 2, 3), three(0, 1, 3)],
        walls: vec![plane(1), plane(2), plane(3)],
        concurrencies: vec![],
    };
    let (ts, walls) = plane_ring(&plain).expect("a plane ring");
    assert_eq!(ts.len(), 3);
    assert_eq!(walls, vec![plane(1), plane(2), plane(3)]);
}
