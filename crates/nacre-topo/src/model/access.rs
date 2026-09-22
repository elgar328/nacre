//! Reading the arenas: by handle, by index, and the counts.

use super::*;

impl Model {
    /// The surface's exact truth — what it *is*, beside the f64 realization
    /// [`Model::surface_cache`] returns. Total: a surface without a truth is unrepresentable — the
    /// truth **is** the arena entry a handle names.
    ///
    /// ★★★★ **Which door answers which question** (measured).
    /// **Classification, comparison and branching ask the truth**; display, tessellation,
    /// bounding and measurement ask [`Model::surface_cache`]. A *kind* — "is this face planar"
    /// — is a fact about what the surface **is**, so it is asked here even though the cache
    /// would answer it correctly (the two can never disagree: `push_plane_raw` and
    /// `push_cylinder_raw` pair them by type). The difference is when a wrong answer becomes
    /// possible: a new surface kind lands in *this* enum first, so a `match` here goes red a
    /// step before one on the cache does.
    #[inline]
    pub fn surface(&self, h: Handle<Surface>) -> &Surface {
        self.surfaces.get(h)
    }

    /// The seeded world plane whose **normal** runs along `axis` — `Z` names the XY plane
    /// (z = 0), `X` the YZ plane, `Y` the ZX plane. Deterministic (handles 0–2, pushed by
    /// [`Model::new`]), so a handle in an op log and one from a live session agree.
    pub fn world_plane(&self, axis: nacre_exact::Axis) -> Handle<Surface> {
        let ix = match axis {
            nacre_exact::Axis::Z => 0usize,
            nacre_exact::Axis::X => 1,
            nacre_exact::Axis::Y => 2,
        };
        let h = self.world_planes[ix];
        debug_assert!(
            matches!(self.surface(h), Surface::Plane { motion: None, .. }),
            "seed handles must stay the world planes"
        );
        h
    }

    /// **Debug-only: the cross-store guard for an index-parallel cache read.**
    ///
    /// A cache is a plain `Vec` indexed by `h.index()`, so reading it alone accepts a handle
    /// minted by *another* model and silently answers with the wrong cell. [`Store::get`] is
    /// where that guard lives (*"Handle was minted by a different Store"*), so a cache read asks
    /// its store first and throws the answer away.
    ///
    /// ★ Measured: a `Store::get` carries this guard for free; indexing a cache alone drops
    /// it, and a foreign handle reaches `surface_cache`/`vertex_point`/`edge_curve` without a
    /// sound. Release builds pay nothing — `Store::get`'s assertion is `cfg(debug_assertions)`
    /// and so is this call.
    #[inline]
    #[cfg(debug_assertions)]
    pub(super) fn debug_guard<T>(store: &Store<T>, h: Handle<T>) {
        let _ = store.get(h);
    }

    /// The surface a handle names, realized — the **f64 cache** of [`Model::surface`]'s
    /// answer.
    ///
    /// ★★★★ **What belongs here**: display, tessellation, bounding and
    /// measurement — every question whose answer is a number a rounded copy can carry, plus the
    /// one that compares a cached coordinate against its cached carrier. A *kind* question does
    /// not, even though it would answer correctly; [`Model::surface`] says why.
    ///
    /// Reading is open; **writing is not** — the store is private, so a surface can only
    /// enter through [`Model::push_plane`]/[`Model::push_cylinder`], which state its truth. The
    /// lock:
    ///
    /// ```compile_fail,E0616
    /// let m = nacre_topo::Model::new();
    /// let _ = m.surfaces.len(); // private field — read through `surface_cache`/`surface_count`
    /// ```
    #[inline]
    pub fn surface_cache(&self, h: Handle<Surface>) -> &nacre_geom::Surface {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.surfaces, h);
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "surface cache out of step with the surface store — push surfaces through \
             Model::push_plane / push_cylinder"
        );
        &self.surface_cache[h.index() as usize].realized
    }

    /// How many surfaces the arena holds — live and superseded alike. Handle-validity checks
    /// and the `a_boolean_mints_no_surface` lock read this; nothing iterates the store
    /// (superseded surfaces are still in it — consumers walk the live faces).
    #[inline]
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// The surface an **index** names — how a log's index vocabulary is re-anchored onto the model
    /// being built.
    ///
    /// A handle in an operation log carries only its index across models, so `replay` turns that
    /// index back into a handle of its own arena before applying the operation. `None` past the
    /// end: existence, not legality.
    ///
    /// **This does not open the store.** Reading was already open ([`Model::surface_cache`],
    /// [`Model::surface_count`]); writing still goes only through [`Model::push_plane`] /
    /// [`Model::push_cylinder`], which state the truth. The `compile_fail` lock on
    /// [`Model::surface_cache`] is untouched.
    ///
    /// The one legitimate shape is "re-anchor a log's index onto the model I am building" — using
    /// it to quiet a cross-model panic hides the bug instead of fixing it. Its consumer today is
    /// `nacre-ops`' `rebind`, for the operations that name a plane.
    #[inline]
    #[must_use]
    pub fn surface_handle_at(&self, index: u32) -> Option<Handle<Surface>> {
        self.surfaces.handle_at(index)
    }

    /// The vertex a log's index names — [`Model::surface_handle_at`]'s twin, for the same reason
    /// (`replay` re-anchors an operation's handles onto the model it is rebuilding).
    ///
    /// ★ A bounds check is all this can be, and for vertices that is a weaker guarantee than it
    /// looks: a rejected operation still leaves cells behind (measured 63–84), so a log recorded
    /// across a reject can re-anchor **in range and onto the wrong vertex**. The discipline that
    /// answers it — rebuild from the log before recording again — lives with `replay`, and the
    /// generated-session proptest is what holds it.
    #[inline]
    pub fn vertex_handle_at(&self, index: u32) -> Option<Handle<Vertex>> {
        self.vertices.handle_at(index)
    }

    /// The vertex a handle names — the **truth**, which for a vertex is its definition.
    ///
    /// ★ One door per entity per side:
    /// `x(h)` is the arena entry itself, `x_cache(h)` is what was realized from it, and the pieces
    /// underneath are reached by chaining on the returned type rather than by more doors here.
    #[inline]
    pub fn vertex(&self, h: Handle<Vertex>) -> &Vertex {
        self.vertices.get(h)
    }

    /// The edge a handle names — the truth beside [`Model::edge_curve`]'s cache.
    #[inline]
    pub fn edge(&self, h: Handle<Edge>) -> &Edge {
        self.edges.get(h)
    }

    /// The face a handle names. ★ There is no `face_cache`: a face has no realized twin, which is
    /// why [`Model::push_face`] is worth an invariant rather than a cache pairing.
    #[inline]
    pub fn face(&self, h: Handle<Face>) -> &Face {
        self.faces.get(h)
    }

    /// The shell a handle names — cache-less for the same reason as [`Model::face`].
    #[inline]
    pub fn shell(&self, h: Handle<Shell>) -> &Shell {
        self.shells.get(h)
    }

    /// The solid a handle names. Which solids are *live* is a separate question —
    /// [`Model::live_solids`] answers it; the arena keeps superseded ones forever.
    #[inline]
    pub fn solid(&self, h: Handle<Solid>) -> &Solid {
        self.solids.get(h)
    }

    /// How many vertices the arena holds — live and superseded alike, like
    /// [`Model::surface_count`].
    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// How many edges the arena holds — live and superseded alike.
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// How many faces the arena holds — live and superseded alike.
    #[inline]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// How many shells the arena holds — live and superseded alike.
    #[inline]
    pub fn shell_count(&self) -> usize {
        self.shells.len()
    }

    /// How many solids the arena holds — live and superseded alike. Most of them are dead:
    /// measured, a box moved 120 times leaves 120 superseded solids behind one live one.
    #[inline]
    pub fn solid_count(&self) -> usize {
        self.solids.len()
    }

    /// The edge an **index** names — [`Model::surface_handle_at`]'s twin for edges.
    ///
    /// ★★★★ **This family is also how the whole arena is walked, and that is deliberate.** There
    /// is no `iter()` door and there will not be one: exposing one invites treating the arena as
    /// the model, when most of what it holds is superseded (measured: 99% of the vertices after a
    /// box is moved 120 times). A consumer that wants the live model walks `live_solids` and the
    /// faces under it.
    ///
    /// But a consumer that wants the **arena** — `validate`'s reference-integrity and vertex-def
    /// checks — genuinely needs every cell, dead ones included: a torn page is torn whether or not
    /// anything still points at it. Those checks walk `0..x_count()` through these doors. Filtering
    /// them to the reachable set would make four planted-corruption tests pass while measuring
    /// nothing, because each plants its corruption on an unreachable cell.
    ///
    /// ☑ Walking by index is not a weaker `iter()`: `nacre-store` locks `handle_at` and `iter`
    /// to the same order, by unit test and by proptest, under the note that the door "grants no
    /// new power".
    #[inline]
    #[must_use]
    pub fn edge_handle_at(&self, index: u32) -> Option<Handle<Edge>> {
        self.edges.handle_at(index)
    }

    /// The face an index names — see [`Model::edge_handle_at`] for why this family exists.
    #[inline]
    #[must_use]
    pub fn face_handle_at(&self, index: u32) -> Option<Handle<Face>> {
        self.faces.handle_at(index)
    }

    /// The shell an index names — see [`Model::edge_handle_at`].
    #[inline]
    #[must_use]
    pub fn shell_handle_at(&self, index: u32) -> Option<Handle<Shell>> {
        self.shells.handle_at(index)
    }

    /// The solid an index names — see [`Model::edge_handle_at`].
    #[inline]
    #[must_use]
    pub fn solid_handle_at(&self, index: u32) -> Option<Handle<Solid>> {
        self.solids.handle_at(index)
    }
}
