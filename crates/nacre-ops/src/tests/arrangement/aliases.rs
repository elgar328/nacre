//! The alias table's own bookkeeping.

use super::*;

/// **Joining two components whose roots are both keys is learning** — the fixed-point loop's
/// measure counts it although the table gains no key.
#[test]
fn joining_two_seated_components_is_learned_without_a_new_key() {
    let name = |k: usize| NodeId::three_planes(Canon3::three([0, 1, k]));
    let mut al = Aliases::default();
    al.union_point(name(2), name(3));
    al.union_point(name(4), name(5));
    let (keys, learned) = (al.point.len(), al.learned());
    al.union_point(name(3), name(5));
    assert_eq!(
        al.point.len(),
        keys,
        "both roots were keys: no key is added"
    );
    assert_eq!(al.learned(), learned + 1, "the partition changed");
    assert!(al.one_point(name(2), name(4)));
    al.union_point(name(2), name(5));
    assert_eq!(
        al.learned(),
        learned + 1,
        "one component already: nothing learned"
    );
}
