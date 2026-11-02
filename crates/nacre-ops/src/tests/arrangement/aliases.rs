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
    let learn = |lines: &mut Lines, held: &[usize]| lines.learn(held, 0, |_| true, |_, _, _| 1);
    assert!(learn(&mut lines, &[3, 6, 7]));
    assert!(learn(&mut lines, &[6, 7, 16]));
    assert!(!learn(&mut lines, &[3, 7, 16]), "already one line");
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
    assert_eq!(lines.learnt().count(), 1);
}

/// **A component's representative is the same whichever came first, the union or the line.** Classes
/// `0`, `1`, `2` have no world statement, `3` has one; the line `{0, 1, 3}` is learnt, and a point
/// wears `Pierce{[0, 1]}` (on the line) and `Pierce{[0, 2]}` (not). Every pair of the line places
/// its names by intervals alike, so the line is named by its least pair by number, `[0, 1]`, and a
/// restated name is never greater than the name it restates. A line named by its one named class
/// first, `[0, 3]`, would make the representative `Pierce{[0, 1]}` restated (union first) or
/// `Pierce{[0, 2]}` (line first).
#[test]
fn the_representative_does_not_depend_on_what_was_learnt_first() {
    use crate::combinatorics::QuadRoot::Lo;
    let named = |c: usize| c == 3;
    let (on_line, off_line) = (NodeId::pierce(0, 1, 0, Lo), NodeId::pierce(0, 2, 0, Lo));
    let table = |union_first: bool| {
        let mut al = Aliases {
            nameless: (0..4).map(|c| !named(c)).collect(),
            ..Aliases::default()
        };
        if union_first {
            al.union_point(on_line, off_line);
        }
        al.learn_line_by(&[0, 1, 3], 2, named, |_, _, _| 1);
        if !union_first {
            al.union_point(on_line, off_line);
        }
        al.canon_point(off_line)
    };
    assert_eq!(table(true), table(false));
    assert_eq!(table(true), on_line, "the line's name is its least pair");
}

/// **A line with two named classes is named by them**, and a class calls the line by its least
/// partner in the same order: a named class by a named partner where there is one.
#[test]
fn a_line_is_named_by_its_least_pair_and_walls_follow() {
    let mut lines = Lines::default();
    // `0` has no world statement; `1` and `2` do.
    assert!(lines.learn(&[0, 1, 2], 9, |c| c != 0, |_, _, _| 1));
    let line = lines.line(0, 1).expect("learnt");
    assert_eq!(line.name(), [1, 2]);
    assert_eq!(
        lines.wall(1, 0),
        2,
        "a named class calls it by its named partner"
    );
    assert_eq!(lines.wall(2, 0), 1);
    assert_eq!(
        lines.wall(0, 2),
        1,
        "a class without one by its least partner"
    );
}
