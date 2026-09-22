//! Vertices, edges, faces and shells: the pushes, their caches, and reachability.

use super::*;

impl Model {
    /// A vertex's cache — **the one road to a coordinate from a vertex**, read from the
    /// index-parallel store [`Model::push_vertex`] fills. The coordinate piece is
    /// [`Model::vertex_point`] and the proven bound is [`PointCache::bound`]; read the whole when
    /// the variant itself — what the cache *knows* — is the question.
    #[inline]
    pub fn vertex_cache(&self, vh: Handle<Vertex>) -> &PointCache {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.vertices, vh);
        debug_assert_eq!(
            self.vertex_cache.len(),
            self.vertices.len(),
            "vertex cache out of step with the vertex store — push vertices through Model::push_vertex"
        );
        &self.vertex_cache[vh.index() as usize]
    }

    /// A vertex's realized coordinate — [`Model::vertex_cache`]'s coordinate piece.
    #[inline]
    pub fn vertex_point(&self, vh: Handle<Vertex>) -> Point3 {
        self.vertex_cache(vh).coord()
    }

    /// Push a vertex: its definition (the truth) plus its cache. **The road that builds a vertex,
    /// and the only one** — nothing else may add to the store.
    ///
    /// ★ There is no `rebuild_vertex_cache` that *replaces* what is here: an operation realizes the
    /// definition *before* it pushes (`nacre_ops`'s push funnel), so what arrives is already the
    /// realization's memo ([`PointCache::Bounded`]) or a named reason it is not
    /// ([`PointCache::Ceiling`], [`PointCache::Unrealized`]). The old warrant for the absence
    /// (*"238 of 1,992 differ from a naive re-solve"*) compared against a naive re-solve; the
    /// exact-rounding realization is not one, and the census says so vertex by vertex.
    ///
    /// ⚠ **A second writer does exist, and it only ever raises**: `nacre_ops::refine_vertex_cache`
    /// lifts a `Ceiling` to `Bounded` by paying for the realization this road would not. It cannot
    /// reach the other variants and cannot move a coordinate anywhere but closer to the truth, so
    /// what stays true of the cache behind this door is *append-only in accuracy*, not in bytes.
    pub fn push_vertex(&mut self, def: Vertex, cache: PointCache) -> Handle<Vertex> {
        match def {
            Vertex::ThreePlane([a, b, c]) => debug_assert!(
                a != b && b != c && a != c,
                "a three-plane definition needs three distinct planes"
            ),
            Vertex::OnSeam([a, b]) => {
                debug_assert!(a != b, "a seam vertex needs two distinct carriers")
            }
            Vertex::Pierce {
                planes: [a, b],
                cylinder,
                ..
            } => debug_assert!(
                a.index() < b.index() && cylinder != a && cylinder != b,
                "a pierce definition needs two sorted distinct planes and a distinct cylinder"
            ),
        }
        let h = self.vertices.push(def);
        self.vertex_cache.push(cache);
        h
    }

    /// **Raise a [`PointCache::Ceiling`] to [`PointCache::Bounded`]** — the second writer of the
    /// vertex cache, and the only thing it can do is make a coordinate more accurate.
    ///
    /// ★ Deliberately not `set_vertex_cache`. It cannot reach the other two variants: an
    /// `Unrealized` coordinate has no realization behind it (raising it would be inventing one) and
    /// a `Bounded` one is already the realization. The debug assert is the type saying so out loud
    /// — the caller that pays for the realization is `nacre_ops::refine_vertex_cache`, and its own
    /// contract is that it looks at nothing else.
    ///
    /// ⚠ **Anything derived from this coordinate is now stale**, starting with the edge curves the
    /// endpoints decide ([`Model::rebuild_edge_cache`]). The paying caller is responsible for that,
    /// because it is the one that knows whether it moved anything at all.
    pub fn refine_vertex_cache(&mut self, vh: Handle<Vertex>, coord: Point3, bound: [Mag; 3]) {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.vertices, vh);
        let slot = &mut self.vertex_cache[vh.index() as usize];
        debug_assert!(
            matches!(slot, PointCache::Ceiling { .. }),
            "only a Ceiling is raised, not {slot:?}"
        );
        *slot = PointCache::Bounded { coord, bound };
    }

    /// An edge's curve — **the one road to a curve from an edge**, read from the
    /// index-parallel cache [`Model::push_edge`] fills.
    #[inline]
    pub fn edge_curve(&self, e: Handle<Edge>) -> &Curve {
        #[cfg(debug_assertions)]
        Self::debug_guard(&self.edges, e);
        debug_assert_eq!(
            self.edge_cache.len(),
            self.edges.len(),
            "edge cache out of step with the edge store — push edges through Model::push_edge"
        );
        &self.edge_cache[e.index() as usize].curve
    }

    /// Push an edge: canonicalize the carrier pair, derive its curve, fill the cache — the one
    /// write road. `None` when the curve does not derive, which for the line arms means
    /// coincident endpoints (a zero-length edge); the caller maps that to its own reject
    /// (`DegenerateGeometry` / `ZeroLengthEdge`). ★ A rim is `[v, v]` and NOT degenerate — the
    /// circle arm never reads the endpoints (see [`Model::derive_edge_curve`]).
    pub fn push_edge(
        &mut self,
        surfaces: [Handle<Surface>; 2],
        vertices: [Handle<Vertex>; 2],
    ) -> Option<Handle<Edge>> {
        let surfaces = Edge::carrier_pair(surfaces[0], surfaces[1]);
        let curve = self.derive_edge_curve(surfaces, vertices)?;
        let h = self.edges.push(Edge { surfaces, vertices });
        self.edge_cache.push(EdgeCache { curve });
        Some(h)
    }

    /// Push a face — **the door that carries the face invariant**.
    ///
    /// ★ Unlike [`Model::push_vertex`] this pairs no cache with the truth, because a face has
    /// none (`surface_cache`/`edge_cache`/`vertex_cache` exist; a face cache does not). What this
    /// door is worth is therefore the **invariant**, not cache-sync: a face whose loops do not
    /// close is a face no consumer can walk, and this says so at the moment it is built.
    ///
    /// ☑ **Measured before it shipped** — every face the production road builds satisfies this,
    /// live or superseded, across the boolean / twice-cut / cylinder / tilted-frame fixtures. The
    /// one fixture that did not was a deliberately malformed «franken» face in a `transform` test,
    /// and it was reshaped rather than exempted, so the invariant holds with no hole behind it.
    pub fn push_face(&mut self, face: Face) -> Handle<Face> {
        debug_assert!(
            (face.surface.index() as usize) < self.surfaces.len(),
            "a face names a surface the arena does not hold"
        );
        debug_assert!(
            !face.outer.half_edges.is_empty(),
            "a face's outer loop has no half-edges"
        );
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            debug_assert!(
                !lp.half_edges.is_empty(),
                "a face's inner loop has no half-edges"
            );
            let n = lp.half_edges.len();
            for i in 0..n {
                let (he, nx) = (&lp.half_edges[i], &lp.half_edges[(i + 1) % n]);
                let [a, b] = self.edges.get(he.edge).vertices;
                let end = if he.forward { b } else { a };
                let [c, d] = self.edges.get(nx.edge).vertices;
                let start = if nx.forward { c } else { d };
                debug_assert_eq!(end, start, "a face loop does not close at half-edge {i}");
            }
        }
        self.faces.push(face)
    }

    /// Push a shell — [`Model::push_face`]'s twin, and a cache-less door for the same reason.
    ///
    /// ☑ **Measured before it shipped**: across the same fixtures no shell the production road
    /// builds is empty or names a face outside the arena, live or superseded (36 shells, 0 and 0).
    pub fn push_shell(&mut self, shell: Shell) -> Handle<Shell> {
        debug_assert!(
            !shell.faces.is_empty(),
            "a shell with no faces bounds nothing"
        );
        debug_assert!(
            shell
                .faces
                .iter()
                .all(|f| (f.index() as usize) < self.faces.len()),
            "a shell names a face the arena does not hold"
        );
        self.shells.push(shell)
    }

    /// Discard every edge-curve cache and derive it afresh — the «cache, not truth» warrant:
    /// nothing is lost, because nothing there was truth.
    /// ⚠★★★ **Only the reachable edges are re-derived, and a superseded one keeps what it has.**
    /// The arena is append-only, so most of what is in it is dead: measured, a boolean corner has
    /// 24 dead edges of 48, a twice-cut one 60 of 108, a thrice-moved box 36 of 48. Re-deriving
    /// those costs the work twice over and — once a caller can *move* a coordinate
    /// ([`Model::refine_vertex_cache`]) — risks a dead cell's endpoints becoming coincident, where
    /// the derivation answers `None` and this would die on the `expect`.
    ///
    /// ☑ That `None` does not happen today: measured over the same fixtures, **zero** stored edges
    /// fail to derive, dead or live. The filter is not a workaround for a live failure — it is what
    /// makes the failure structurally unreachable, because a superseded edge is never derived again.
    pub fn rebuild_edge_cache(&mut self) {
        let reach = self.reachable();
        self.edge_cache = self
            .edges
            .iter()
            .map(|(eh, e)| match reach.edges.contains(&eh) {
                true => EdgeCache {
                    curve: self
                        .derive_edge_curve(e.surfaces, e.vertices)
                        .expect("every live edge derives its curve"),
                },
                false => self.edge_cache[eh.index() as usize].clone(),
            })
            .collect();
    }

    /// The vertex a half-edge starts at: its edge's `vertices[0]` when the use runs
    /// forward, `vertices[1]` when it runs back.
    ///
    /// A traversal accessor, not an analysis — the same kind of thing as
    /// [`Model::reachable`], and the reason it lives here: walking a loop's
    /// corners is the first thing every consumer above does, and it was written
    /// twice (with two different failure policies) before this existed.
    ///
    /// Total: every edge has both endpoints by type, and a closed rim
    /// states that by repeating its seam vertex — `[v, v]`, so the start is `v` either way.
    /// A standalone full circle with no seam would have no start, but the type does not
    /// express that form: it is a wireframe/open-shell element and a v1 non-goal. Nothing
    /// guards it because nothing can build it.
    #[inline]
    pub fn he_start(&self, he: HalfEdge) -> Handle<Vertex> {
        let [a, b] = self.edges.get(he.edge).vertices;
        if he.forward { a } else { b }
    }

    /// **An edge's curve, derived from what the model already holds** — the edge's
    /// truth/cache split: the carriers and the endpoints decide the curve, so the stored one
    /// is a cache that can be discarded and regenerated.
    ///
    /// Dispatch by carrier type:
    /// * **Plane × Plane** (and the self-adjacent cylinder **seam**): the line through the two
    ///   endpoint coordinates — the very expression a producer would build the stored
    ///   curve with, so the derivation is bit-identical, and an endpoint pair that coincides is the
    ///   `None` (a degenerate line — the check lives in this arm only).
    /// * **Plane × Cylinder** (a rim, or an arc of one): the circle centred where the cylinder's
    ///   axis crosses the cap plane, with the **cylinder's** frame (`axis direction`, `ref_dir`,
    ///   `radius`) — the same parameters `add_cylinder` builds the stored rims from, so
    ///   tessellation's `θ` parameterization is preserved. The endpoints are not read: a full rim
    ///   is a closed edge (`[v, v]`), which is not a degeneracy.
    ///
    ///   ★★ **On a circle carrier, the vertex *order* says which arc**: two distinct
    ///   endpoints cut a circle into two pieces the endpoints alone cannot tell apart, so
    ///   `[A, B]` means the piece from A to B **counter-clockwise about the axis direction**,
    ///   and the two complementary arcs between one vertex pair are the two orders. Producers
    ///   uphold this (`boolean`'s edge welding keys arcs in CCW order); the curve stored here is
    ///   the whole circle either way. Both consumers read it the same way: tessellation's
    ///   `sample_edge` walks `θ(v0) → θ(v0) + Δθ` with `Circle::angle_of` as the one spelling of
    ///   θ, and `validate`'s `loop_winding` adds each arc's circular segment with Δθ from the
    ///   same order.
    /// * **Cylinder × Cylinder**: no producer builds one before M6 — `None`, honestly.
    ///
    /// ★ The M3 rim population is axis-perpendicular by construction; a *tilted* plane over a
    /// cylinder would cross in an ellipse, which this arm cannot express (M6). The debug
    /// assertion keeps that boundary visible.
    pub fn derive_edge_curve(
        &self,
        surfaces: [Handle<Surface>; 2],
        vertices: [Handle<Vertex>; 2],
    ) -> Option<Curve> {
        let endpoints_line = || -> Option<Curve> {
            let p0 = self.vertex_point(vertices[0]);
            let p1 = self.vertex_point(vertices[1]);
            Some(Curve::Line(Line::through_points(p0, p1)?))
        };
        // ★★ **Both kinds asked of the cache here, deliberately**. Every arm
        // reads cache *values* out of the very binding it matched — `p.normal()`, `c.axis()`,
        // `c.radius()` — so dispatching on the truth would double the lookups and leave two
        // matches whose agreement no reader could check. The rule sends *kind questions* to the
        // truth; this match's answer is a curve, and the kinds only choose how to derive it.
        match (
            self.surface_cache(surfaces[0]),
            self.surface_cache(surfaces[1]),
        ) {
            (nacre_geom::Surface::Plane(_), nacre_geom::Surface::Plane(_)) => endpoints_line(),
            (nacre_geom::Surface::Cylinder(_), nacre_geom::Surface::Cylinder(_))
                if surfaces[0] == surfaces[1] =>
            {
                endpoints_line() // the seam — a parameterization joint, straight along the axis
            }
            (nacre_geom::Surface::Plane(p), nacre_geom::Surface::Cylinder(c))
            | (nacre_geom::Surface::Cylinder(c), nacre_geom::Surface::Plane(p)) => {
                let axis = c.axis();
                // A plane **parallel** to the axis meets the lateral along rulings — straight,
                // so the endpoints decide, exactly like the seam arm above (the rulings
                // ladder). Same scale convention as the ⊥ assertion below, so the band between
                // the two tests is symmetric and only a genuinely tilted plane (an ellipse)
                // falls through to it.
                {
                    let n = p.normal();
                    let d = axis.direction();
                    if n.dot(d).powi(2) <= 1e-18 * n.norm_squared() * d.norm_squared() {
                        return endpoints_line();
                    }
                }
                debug_assert!(
                    {
                        let n = p.normal();
                        let d = axis.direction();
                        n.cross(d).norm_squared() <= 1e-18 * n.norm_squared()
                    },
                    "a tilted plane over a cylinder crosses in an ellipse — M6-3, no producer yet"
                );
                let center = nacre_geom::intersect::line_plane(&axis, p)?;
                Some(Curve::Circle(Circle::from_center_normal(
                    center,
                    axis.direction(),
                    c.ref_dir(),
                    c.radius(),
                )?))
            }
            (nacre_geom::Surface::Cylinder(_), nacre_geom::Surface::Cylinder(_)) => None, // two distinct cylinders: M6
        }
    }

    /// The handles reachable from the live solids — the live model.
    ///
    /// Superseded cells left in the append-only arena are excluded (nothing live
    /// references them). Every step is bounds-guarded, so this is safe even on a
    /// corrupt or partially-built model (a dangling handle simply prunes that
    /// branch; `validate`'s reference-integrity check reports it separately).
    pub fn reachable(&self) -> Reachable {
        let mut r = Reachable::default();
        for &solid_h in &self.live_solids {
            if !in_bounds(solid_h, &self.solids) {
                continue;
            }
            let solid = self.solids.get(solid_h);
            for &shell_h in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
                if !in_bounds(shell_h, &self.shells) || !r.shells.insert(shell_h) {
                    continue;
                }
                for &face_h in &self.shells.get(shell_h).faces {
                    if !in_bounds(face_h, &self.faces) || !r.faces.insert(face_h) {
                        continue;
                    }
                    let face = self.faces.get(face_h);
                    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                        for he in &lp.half_edges {
                            if !in_bounds(he.edge, &self.edges) || !r.edges.insert(he.edge) {
                                continue;
                            }
                            for v in self.edges.get(he.edge).vertices {
                                if in_bounds(v, &self.vertices) {
                                    r.vertices.insert(v);
                                }
                            }
                        }
                    }
                }
            }
        }
        r
    }
}
