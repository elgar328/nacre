use super::*;
use proptest::prelude::*;
use std::collections::{HashMap, HashSet};

#[test]
fn push_then_get_roundtrips() {
    let mut s = Store::new();
    let a = s.push("alpha");
    let b = s.push("beta");
    assert_eq!(*s.get(a), "alpha");
    assert_eq!(*s.get(b), "beta");
    assert_eq!(s.len(), 2);
}

#[test]
fn indices_are_sequential_and_distinct() {
    let mut s = Store::new();
    let a = s.push(10);
    let b = s.push(20);
    let c = s.push(30);
    assert_eq!((a.index(), b.index(), c.index()), (0, 1, 2));
    assert_ne!(a, b);
    assert_ne!(b, c);
}

#[test]
fn equality_is_by_index_only() {
    let mut s = Store::new();
    let a = s.push(1);
    let a_copy = a; // Handle is Copy
    assert_eq!(a, a_copy);
}

/// A `Handle` must be `Copy + Eq + Hash` even when the
/// stored type is neither `Eq` nor `Hash`. `NotHashable` holds an `f64`, so
/// if `Handle` derived its impls this would fail to compile — exactly the
/// bug that would break topo's `HashMap<Handle<Edge>, _>`.
#[test]
fn handle_is_hashable_even_when_t_is_not() {
    struct NotHashable(#[allow(dead_code)] f64);
    let mut s = Store::new();
    let h = s.push(NotHashable(0.1));
    let mut set: HashSet<Handle<NotHashable>> = HashSet::new();
    set.insert(h);
    assert!(set.contains(&h));
    let mut map: HashMap<Handle<NotHashable>, u8> = HashMap::new();
    map.insert(h, 7);
    assert_eq!(map.get(&h), Some(&7));
}

#[test]
fn handles_sort_by_insertion_order() {
    let mut s = Store::new();
    let mut hs: Vec<_> = (0..5).map(|v| s.push(v)).collect();
    hs.reverse();
    hs.sort();
    assert_eq!(
        hs.iter().map(|h| h.index()).collect::<Vec<_>>(),
        vec![0, 1, 2, 3, 4]
    );
}

#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "different Store")]
fn cross_store_handle_panics_in_debug() {
    let mut a: Store<u32> = Store::new();
    let b: Store<u32> = Store::new();
    let h = a.push(42);
    // Using store A's handle on store B is a bug — caught in debug builds.
    let _ = b.get(h);
}

/// `handle_at` gives back exactly what `push` returned for that slot.
#[test]
fn handle_at_is_the_handle_push_returned() {
    let mut s: Store<u32> = Store::new();
    let handles: Vec<_> = (0..4).map(|v| s.push(v * 10)).collect();
    for (i, &h) in handles.iter().enumerate() {
        assert_eq!(s.handle_at(i as u32), Some(h));
        assert_eq!(*s.get(s.handle_at(i as u32).unwrap()), (i as u32) * 10);
    }
}

/// A slot that does not exist has no handle — the boundary is `len()`, and an empty store
/// has none at all. This `None` is what lets a caller reject a log naming a cell that this
/// model does not have, instead of fabricating a handle into thin air.
#[test]
fn handle_at_past_the_end_is_none() {
    let empty: Store<u32> = Store::new();
    assert_eq!(empty.handle_at(0), None);
    let mut s: Store<u32> = Store::new();
    s.push(1);
    s.push(2);
    assert!(s.handle_at(1).is_some());
    assert_eq!(s.handle_at(2), None, "the boundary is len(), not len()+1");
    assert_eq!(s.handle_at(u32::MAX), None);
}

/// ★ **The new door grants no new power** — `iter().nth(i)` already minted this handle, so
/// `handle_at` is a name and an `O(1)` path, not an escape hatch. If these two ever
/// disagree, one of them is lying about which store it belongs to.
#[test]
fn handle_at_agrees_with_iter() {
    let mut s: Store<char> = Store::new();
    for c in "nacre".chars() {
        s.push(c);
    }
    for i in 0..s.len() {
        assert_eq!(
            s.handle_at(i as u32),
            s.iter().nth(i).map(|(h, _)| h),
            "slot {i}"
        );
    }
}

/// ★ **The positive control for re-anchoring**: an index taken from store A names *B's*
/// item when re-anchored on B — no panic, because the handle `handle_at` returns is B's by
/// construction. Its twin is [`cross_store_handle_panics_in_debug`] above: without the
/// re-anchoring the very same index still dies at the guard. (They cannot share a body —
/// `should_panic` needs the whole test to die.)
#[test]
fn handle_at_rebinds_to_this_store_not_the_other() {
    let mut a: Store<&str> = Store::new();
    let mut b: Store<&str> = Store::new();
    a.push("from A");
    b.push("from B");
    let h_a = a.handle_at(0).expect("A has slot 0");
    let rebound = b.handle_at(h_a.index()).expect("B has slot 0 too");
    assert_eq!(*b.get(rebound), "from B");
    assert_eq!(
        h_a, rebound,
        "handles compare by index — that is the vocabulary"
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// `handle_at` and `iter` agree over arbitrary stores — the property form of
    /// `handle_at_agrees_with_iter`, which is what makes «no new power» a claim about the
    /// API rather than about one fixture.
    #[test]
    fn handle_at_agrees_with_iter_everywhere(items in prop::collection::vec(any::<i64>(), 0..64)) {
        let mut s = Store::new();
        for &v in &items {
            s.push(v);
        }
        for i in 0..items.len() {
            prop_assert_eq!(s.handle_at(i as u32), s.iter().nth(i).map(|(h, _)| h));
        }
        prop_assert_eq!(s.handle_at(items.len() as u32), None);
    }

    /// Every handle retrieves exactly the item that was pushed at its slot,
    /// for an arbitrary sequence of pushes.
    #[test]
    fn all_handles_retrieve_their_own_item(items in prop::collection::vec(any::<i64>(), 0..256)) {
        let mut s = Store::new();
        let handles: Vec<_> = items.iter().map(|&v| s.push(v)).collect();
        prop_assert_eq!(s.len(), items.len());
        for (h, &expected) in handles.iter().zip(items.iter()) {
            prop_assert_eq!(*s.get(*h), expected);
        }
    }

    /// Distinct slots yield distinct (non-equal) handles.
    #[test]
    fn distinct_pushes_are_distinct_handles(n in 0usize..256) {
        let mut s = Store::new();
        let handles: Vec<_> = (0..n).map(|i| s.push(i)).collect();
        let unique: HashSet<_> = handles.iter().copied().collect();
        prop_assert_eq!(unique.len(), n);
    }

    /// `iter` visits items in insertion order and its handles round-trip.
    #[test]
    fn iter_matches_insertion_order(items in prop::collection::vec(any::<i32>(), 0..128)) {
        let mut s = Store::new();
        for &v in &items {
            s.push(v);
        }
        let collected: Vec<_> = s.iter().map(|(h, &v)| (h.index(), v)).collect();
        let expected: Vec<_> = items.iter().enumerate().map(|(i, &v)| (i as u32, v)).collect();
        prop_assert_eq!(collected, expected);
    }
}
