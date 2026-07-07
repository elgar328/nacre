//! Typed-index, append-only storage — the foundation of the nacre kernel.
//!
//! Every kernel object lives exactly once in a [`Store`] and is referenced only
//! by a [`Handle`]. There is no removal or mutation: handles stay valid forever,
//! and identity is an integer comparison (`h1 == h2`), never a floating-point
//! coordinate comparison. See `docs/design.md` §2.
//!
//! `Store`/`Handle` live in this lowest crate on purpose: `nacre-geom` already
//! needs `Handle` (its `Curve::Intersection` holds `Handle<Surface>`), so placing
//! them in `nacre-topo` would create a geom→topo→geom cycle.

use core::cmp::Ordering;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;

/// A permanent reference to an item in a [`Store<T>`].
///
/// This is just a `u32` index plus a type tag — it never holds a `T`. We do
/// **not** `#[derive]` the trait impls: derive would add a spurious `T: Trait`
/// bound (rust-lang/rust#26925), so e.g. `#[derive(Hash)]` would demand
/// `T: Hash`. That breaks `HashMap<Handle<Edge>, _>` (topo `Adjacency`), because
/// `Edge` carries an `f64` tolerance and thus implements neither `Hash` nor `Eq`.
/// The hand-written impls below compare/hash only `index`, so `Handle<T>` is
/// `Copy + Eq + Hash + Ord` for *any* `T`.
///
/// `PhantomData<fn() -> T>` (rather than `PhantomData<T>`) keeps `Handle<T>`
/// covariant in `T` and unconditionally `Send + Sync`, regardless of `T`.
pub struct Handle<T> {
    index: u32,
    /// Debug-only guard: which `Store` minted this handle. Used by
    /// [`Store::get`] to catch a handle applied to the wrong store/model.
    /// Absent in release builds, so `Handle` is `size_of::<u32>()` there.
    #[cfg(debug_assertions)]
    store: StoreId,
    _t: PhantomData<fn() -> T>,
}

impl<T> Handle<T> {
    /// The raw index. Exposed for debugging and OBJ/diagnostic dumps; identity
    /// should still go through `==`, not this.
    #[inline]
    pub fn index(&self) -> u32 {
        self.index
    }
}

// --- Manual trait impls: all key off `index` only, with no bound on `T`. ---

impl<T> Clone for Handle<T> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Handle<T> {}

impl<T> PartialEq for Handle<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}
impl<T> Eq for Handle<T> {}

impl<T> Hash for Handle<T> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}

impl<T> PartialOrd for Handle<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl<T> Ord for Handle<T> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.index.cmp(&other.index)
    }
}

impl<T> fmt::Debug for Handle<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Show the element type name so `Handle<Surface>` reads clearly in logs.
        write!(f, "Handle<{}>({})", short_type_name::<T>(), self.index)
    }
}

/// Append-only storage for a single object type. Push returns a permanent
/// [`Handle`]; there is no `remove`, so every handle stays valid for the life
/// of the store. See `docs/design.md` §2.
pub struct Store<T> {
    items: Vec<T>,
    /// Debug-only identity of this store, stamped into every handle it mints.
    #[cfg(debug_assertions)]
    id: StoreId,
}

impl<T> Store<T> {
    /// Create an empty store.
    #[inline]
    pub fn new() -> Self {
        Store {
            items: Vec::new(),
            #[cfg(debug_assertions)]
            id: StoreId::next(),
        }
    }

    /// Append `item` and return a handle to it. The handle's index is the item's
    /// position, so replaying the same pushes reproduces identical handles.
    #[inline]
    pub fn push(&mut self, item: T) -> Handle<T> {
        let index =
            u32::try_from(self.items.len()).expect("Store overflow: more than u32::MAX items");
        self.items.push(item);
        Handle {
            index,
            #[cfg(debug_assertions)]
            store: self.id,
            _t: PhantomData,
        }
    }

    /// Borrow the item a handle points to. Always valid — there is no removal.
    ///
    /// In debug builds this asserts the handle was minted by *this* store,
    /// catching the "used model A's handle on model B" bug (`docs/design.md` §2).
    #[inline]
    pub fn get(&self, h: Handle<T>) -> &T {
        #[cfg(debug_assertions)]
        assert_eq!(
            h.store, self.id,
            "Handle was minted by a different Store (cross-model/cross-store misuse)"
        );
        &self.items[h.index as usize]
    }

    /// Number of items stored.
    #[inline]
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the store is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Iterate over `(handle, &item)` pairs in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (Handle<T>, &T)> {
        self.items.iter().enumerate().map(|(i, item)| {
            let h = Handle {
                index: i as u32,
                #[cfg(debug_assertions)]
                store: self.id,
                _t: PhantomData,
            };
            (h, item)
        })
    }
}

impl<T> Default for Store<T> {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

/// Concise, `T`-agnostic debug output (item count only) — so `Store<T>` and
/// aggregates over it are `Debug` without requiring `T: Debug`, and printing a
/// big store doesn't dump every element.
impl<T> fmt::Debug for Store<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Store")
            .field("len", &self.items.len())
            .finish()
    }
}

/// Strip the module path from a type name (`nacre_geom::Surface` -> `Surface`)
/// for readable `Handle` debug output.
fn short_type_name<T>() -> &'static str {
    let full = core::any::type_name::<T>();
    full.rsplit("::").next().unwrap_or(full)
}

/// Process-local store identity, minted from a monotonic counter. Handed to
/// [`Store::new`] so distinct stores (and thus distinct models) get distinct
/// ids without any `Date`/random source. Debug builds only.
#[cfg(debug_assertions)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct StoreId(u64);

#[cfg(debug_assertions)]
impl StoreId {
    fn next() -> Self {
        use core::sync::atomic::{AtomicU64, Ordering};
        // Start at 1 so a zeroed handle is never mistaken for a valid one.
        static COUNTER: AtomicU64 = AtomicU64::new(1);
        StoreId(COUNTER.fetch_add(1, Ordering::Relaxed))
    }
}

#[cfg(test)]
mod tests {
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

    /// The F1 regression: a `Handle` must be `Copy + Eq + Hash` even when the
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

    proptest! {
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
}
