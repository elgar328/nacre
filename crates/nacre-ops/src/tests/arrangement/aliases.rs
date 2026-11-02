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

/// **A line learnt in pieces is one line on every class** — two pencils that share two classes
/// are one line, so a class that heard of only one piece still names the line by its least
/// holder. A table kept per class answers `6` on class `16` and leaves `16` out of class `3`'s
/// family: the suite and census carry that shape (the line `{3, 6, 7, 16}`).
#[test]
fn a_line_learnt_in_pieces_is_one_line_on_every_class() {
    let mut lines = Lines::default();
    assert!(lines.learn(&[3, 6, 7]));
    assert!(lines.learn(&[6, 7, 16]));
    assert!(!lines.learn(&[3, 7, 16]), "already one line");
    for class in [3, 6, 7, 16] {
        let others: Vec<usize> = [3, 6, 7, 16].into_iter().filter(|&x| x != class).collect();
        for &w in &others {
            assert_eq!(lines.family(class, w), others, "class {class}, wall {w}");
            assert_eq!(lines.wall(class, w), others[0], "class {class}, wall {w}");
        }
    }
    assert_eq!(
        lines.wall(3, 9),
        9,
        "a line no pencil holds keeps its pair's name"
    );
    assert_eq!(lines.held().count(), 1);
}
