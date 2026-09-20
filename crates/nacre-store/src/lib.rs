//! Typed-index, append-only storage — the foundation of the nacre kernel.
//!
//! Every kernel object lives exactly once in a [`Store`] and is referenced only
//! by a [`Handle`]. There is no removal or mutation: handles stay valid forever,
//! and identity is an integer comparison (`h1 == h2`), never a floating-point
//! coordinate comparison.
//!
//! `Store`/`Handle` live in this lowest crate on purpose: a typed-index arena knows
//! nothing of geometry or topology, and everything above it — `nacre-topo`'s arenas,
//! `nacre-ops`, `nacre-validate`, `nacre-step` — names handles freely.
//!
//! (The older note here gave a different reason — that `nacre-geom` would need `Handle`
//! for a `Curve::Intersection` variant. It will not: a handle names an arena entry, and
//! the surface arena holds the **truth**, a `nacre-topo` type. A crate below topo cannot
//! name it, so geom stays `Handle`-free.)

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
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
/// of the store.
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
    /// catching the "used model A's handle on model B" bug.
    #[inline]
    pub fn get(&self, h: Handle<T>) -> &T {
        #[cfg(debug_assertions)]
        assert_eq!(
            h.store, self.id,
            "Handle was minted by a different Store (cross-model/cross-store misuse)"
        );
        &self.items[h.index as usize]
    }

    /// **This store's handle for a slot named by index** — `None` past the end.
    ///
    /// An operation log names cells by *index*, not by handle: a `Handle`'s identity is its
    /// index (the `Eq`/`Hash` impls below use nothing else), and a log outlives the model it
    /// was recorded against. So `replay`, which builds a fresh model from scratch, re-anchors
    /// each of the log's indices onto the model it is building — this is that re-anchoring.
    ///
    /// **It does not weaken [`Store::get`]'s guard.** The handle returned is this store's by
    /// construction, so nothing crosses models here; and the ability itself is not new —
    /// `iter().nth(i)` already mints exactly this handle. What this adds is a name, `O(1)`,
    /// an `Option` where `nth` gives one for a different reason, and somewhere to write the
    /// contract down.
    ///
    /// **The one legitimate shape is "re-anchor a log's index onto the model I am building".**
    /// Reaching for this because a cross-store panic fired is hiding the bug it caught. (The
    /// throwaway-`Store` trick in `nacre-validate`'s fixtures is a *different* job — it needs
    /// handles that point **past** the end, which this cannot make by construction.)
    #[inline]
    #[must_use]
    pub fn handle_at(&self, index: u32) -> Option<Handle<T>> {
        ((index as usize) < self.items.len()).then_some(Handle {
            index,
            #[cfg(debug_assertions)]
            store: self.id,
            _t: PhantomData,
        })
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
#[path = "tests/lib.rs"]
mod tests;
