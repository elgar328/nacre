//! H2/H3 — explicit (reference-based) sharing.
//!
//! Identity is by [`Handle`] (an append-only reference), **not** by coordinate.
//! Sharing happens only when a construction *references* an existing element;
//! never because two elements happen to land on the same coordinate (§5). That is
//! what makes a later coplanarity/coincidence test an O(1) `handle == handle`
//! instead of a coordinate comparison, and what keeps a dimension edit from
//! silently changing topology (the global-auto-merge / TNP conflict, §5/§7).

use crate::Rat;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

/// A typed, append-only reference. Identity is the index — two handles are equal
/// iff they refer to the same stored element, regardless of its coordinates.
/// Hand-written impls avoid spurious `T: Copy/Eq/...` bounds from `derive`.
pub struct Handle<T> {
    idx: usize,
    _t: PhantomData<fn() -> T>,
}

impl<T> Handle<T> {
    fn new(idx: usize) -> Self {
        Handle {
            idx,
            _t: PhantomData,
        }
    }

    pub fn index(self) -> usize {
        self.idx
    }
}

impl<T> Clone for Handle<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Handle<T> {}
impl<T> PartialEq for Handle<T> {
    fn eq(&self, other: &Self) -> bool {
        self.idx == other.idx
    }
}
impl<T> Eq for Handle<T> {}
impl<T> Hash for Handle<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.idx.hash(state);
    }
}
impl<T> std::fmt::Debug for Handle<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Handle({})", self.idx)
    }
}

/// Append-only store; `push` hands back a stable [`Handle`].
pub struct Store<T> {
    items: Vec<T>,
}

impl<T> Store<T> {
    pub fn new() -> Self {
        Store { items: Vec::new() }
    }

    pub fn push(&mut self, v: T) -> Handle<T> {
        let h = Handle::new(self.items.len());
        self.items.push(v);
        h
    }

    pub fn get(&self, h: Handle<T>) -> &T {
        &self.items[h.idx]
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

impl<T> Default for Store<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// A 2D sketch point (exact rational coordinate).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point2 {
    pub xy: (Rat, Rat),
}

/// A segment *referencing* two existing points by handle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Segment {
    pub a: Handle<Point2>,
    pub b: Handle<Point2>,
}

/// A closed loop *referencing* existing segments by handle.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Loop {
    pub edges: Vec<Handle<Segment>>,
}

/// A minimal 2D sketch: append-only stores of points and segments.
pub struct Sketch {
    pub points: Store<Point2>,
    pub segments: Store<Segment>,
}

impl Sketch {
    pub fn new() -> Self {
        Sketch {
            points: Store::new(),
            segments: Store::new(),
        }
    }

    /// Add a *new* point (fresh handle). Never merges by coordinate — two points
    /// at the same coordinate are distinct elements unless one *references* the
    /// other's handle.
    pub fn add_point(&mut self, x: Rat, y: Rat) -> Handle<Point2> {
        self.points.push(Point2 { xy: (x, y) })
    }

    /// Add a segment *referencing* existing point handles (sharing by reference).
    pub fn add_segment(&mut self, a: Handle<Point2>, b: Handle<Point2>) -> Handle<Segment> {
        self.segments.push(Segment { a, b })
    }
}

impl Default for Sketch {
    fn default() -> Self {
        Self::new()
    }
}

/// Coincidence by *reference*: O(1) handle equality, no coordinate comparison.
/// Shared elements are equal; independently-built ones are not, even at the same
/// coordinate — the boolean judges those at contact time (§5/§6), and a later
/// dimension edit can move them apart without a topology change (§7 TNP).
pub fn coincident<T>(a: Handle<T>, b: Handle<T>) -> bool {
    a == b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(n: i128) -> Rat {
        Rat::from_int(n)
    }

    /// Build a unit square and return its four corner handles (CCW) plus edges.
    fn unit_square(s: &mut Sketch) -> ([Handle<Point2>; 4], [Handle<Segment>; 4]) {
        let c = [
            s.add_point(r(0), r(0)),
            s.add_point(r(1), r(0)),
            s.add_point(r(1), r(1)),
            s.add_point(r(0), r(1)),
        ];
        let e = [
            s.add_segment(c[0], c[1]),
            s.add_segment(c[1], c[2]),
            s.add_segment(c[2], c[3]),
            s.add_segment(c[3], c[0]),
        ];
        (c, e)
    }

    /// H3: sketching "the same square" by *referencing* the existing corners and
    /// edges adds **zero** new points and segments — a second loop reuses every
    /// handle.
    #[test]
    fn sketch_on_existing_reuses_handles_zero_new() {
        let mut s = Sketch::new();
        let (_c, e) = unit_square(&mut s);
        assert_eq!(s.points.len(), 4);
        assert_eq!(s.segments.len(), 4);

        // A face drawn on the same boundary = a loop referencing the same edges.
        let loop_a = Loop { edges: e.to_vec() };
        let loop_b = Loop { edges: e.to_vec() };

        // Zero new elements: nothing was pushed to either store.
        assert_eq!(s.points.len(), 4, "no new points");
        assert_eq!(s.segments.len(), 4, "no new segments");
        // The two loops share every element by handle identity.
        assert!(
            loop_a
                .edges
                .iter()
                .zip(&loop_b.edges)
                .all(|(a, b)| coincident(*a, *b))
        );
    }

    /// H2: a second square built by an *independent* construction that happens to
    /// occupy the same coordinates is **not** merged — its points are distinct
    /// elements. Same coordinate, different handle; kept separate (§5).
    #[test]
    fn accidental_coincidence_is_not_merged() {
        let mut s = Sketch::new();
        let (c, _e) = unit_square(&mut s);

        // A different construction that lands on the same coordinates.
        let c2 = [
            s.add_point(r(0), r(0)),
            s.add_point(r(1), r(0)),
            s.add_point(r(1), r(1)),
            s.add_point(r(0), r(1)),
        ];

        // Eight distinct points now exist — nothing auto-merged.
        assert_eq!(s.points.len(), 8);
        for i in 0..4 {
            // coordinates are equal …
            assert_eq!(s.points.get(c[i]).xy, s.points.get(c2[i]).xy);
            // … but the elements are distinct (identity by reference, not coord).
            assert!(!coincident(c[i], c2[i]));
        }
    }

    /// The payoff (§5): a coincidence/coplanarity test is O(1) handle equality,
    /// never a coordinate comparison. Referencing an element makes it *the same*
    /// element; the check reads no coordinates.
    #[test]
    fn coincidence_is_handle_equality_not_coordinates() {
        let mut s = Sketch::new();
        let p = s.add_point(r(3), r(5));
        let shared = p; // referencing the same handle
        let independent = s.add_point(r(3), r(5)); // same coord, new element

        assert!(coincident(p, shared)); // shared reference → coincident
        assert!(!coincident(p, independent)); // accidental coord match → distinct
    }
}
