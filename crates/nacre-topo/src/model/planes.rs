//! Stating a plane or a cylinder: interning, the through-points form, and the vertex meets.

use super::*;

impl Model {
    /// Push a **plane**, stating its truth outright: three exact points and the motion they are
    /// written before (`None` = the world), with the producer's own `f64` figure as the
    /// `fallback` cache. The truth is not optional: there is no point-less plane.
    ///
    /// ★ **`fallback` stands only where the plane has no world name** (a chain that does not fold,
    /// a mixed-frame `Through`) — everywhere else the cache is the truth's realization
    /// (`derive_surface_cache`), and the figure handed in is dropped.
    ///
    /// ★★★★★ **The points are the only thing a producer states.** The canonical name
    /// ([`Model::surface_name`]) is *derived* here, from those points, by
    /// [`nacre_exact::plane_name_exact`] — so a plane cannot be described two ways, because there
    /// is only one place to describe it.
    ///
    /// ★ **Stated in the frame `motion` names** — the world for `None`, the pre-motion frame
    /// otherwise. An interned plane keeps the first pusher's triple.
    ///
    /// ★ **The `bool` says the returned surface faces the *other* way** from the statement handed
    /// in (its points and `sense`), and a caller that meets it must record its face
    /// `Orientation::flipped()`. It exists because a plane's canonical form has no direction —
    /// `[0,0,1,−3]` and `[0,0,−1,3]` are one plane — so once identical planes share a handle the
    /// direction has to be reconciled somewhere, and the honest place is where the caller still
    /// knows what it asked for. It is read off the two truths, never the caches.
    ///
    /// ★★ **It is not a dormant path.** Measured over the census corpus, interning hits 3,125
    /// times and **246 of those report `flipped`** — a boss meeting the plate it sits on is one
    /// plane approached from both sides, which is as ordinary as it sounds.
    pub fn push_plane(
        &mut self,
        fallback: nacre_geom::Plane,
        points: [[nacre_exact::Rat; 3]; 3],
        motion: Option<Handle<MotionNode>>,
        sense: Orientation,
    ) -> (Handle<Surface>, bool) {
        // ★★★★★ **The name is derived, so it cannot disagree with the thing it names.**
        // `plane_name_exact` computes the canonical form at unbounded precision — `None` only
        // for collinear points (which no production `PlaneDef` can supply); the vessel
        // (`PlaneName::Narrow | Wide`) always holds the answer, so every plane interns, wide
        // ones included. [`WIDE_PLANES`] counts the names that took the wide vessel.
        let name = nacre_exact::plane_name_exact(points[0], points[1], points[2]);
        self.intern_plane(fallback, name, PlanePoints::Known(points), motion, sense)
    }

    /// **Interning, once — the half every plane producer shares.**
    ///
    /// A producer differs only in *how it derives the name* and *which `PlanePoints` it stores*.
    /// Everything after that — the key, the already-issued reply and its `flipped`, the arena
    /// push, the two side tables, the two counters — is the same, and was duplicated once, which
    /// promptly cost both counters on the new road ([`WIDE_PLANES`] and [`SEEDED_HITS`] were
    /// simply absent from it). Sharing the tail makes losing them structurally impossible rather
    /// than a thing to remember.
    fn intern_plane(
        &mut self,
        fallback: nacre_geom::Plane,
        name: Option<nacre_exact::PlaneName>,
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
        sense: Orientation,
    ) -> (Handle<Surface>, bool) {
        debug_assert!(
            self.cache_agrees_with_sense(&fallback, &points, motion, sense) != Some(false),
            "the stated sense {sense:?} disagrees with the figure the producer built"
        );
        if name.as_ref().is_some_and(|n| n.narrow().is_none()) {
            WIDE_PLANES.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
        let key = name.clone().map(|n| (n, motion));
        if let Some(k) = &key {
            if let Some(&h) = self.surface_ids.get(k) {
                if (h.index() as usize) < 3 {
                    SEEDED_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                // ★ The incoming statement is **dropped here** — the survivor's cache stands. Its
                // normal is the shared name's; its anchor is the first pusher's first point, and
                // [`Model::count_discarded_cache`] is the only measurement of how often the
                // incoming statement would have anchored elsewhere.
                self.count_discarded_cache(h, &points, motion, &fallback);
                // Same plane, already issued. The canonical form says nothing about direction, so
                // report whether the survivor points the other way and let the caller spell its
                // outward the other way round.
                let flipped = self
                    .flipped_by_truth(h, &points, sense)
                    .expect("a named plane's points span a direction");
                return (h, flipped);
            }
        }
        // One clone per push — the name is derived once here, never on a judging loop.
        let h = self.push_plane_raw(points, motion, sense, name, fallback);
        if let Some(k) = key {
            self.surface_ids.insert(k, h);
        }
        (h, false)
    }

    /// **`flipped`, from the two truths** — whether the survivor of an interning hit faces the
    /// other way from the statement that just arrived. The key is `(name, motion)`, so both
    /// statements ride one chain and their frame directions compare as their world ones do:
    /// flipped ⇔ the two points' directions agree exactly when the two senses differ.
    ///
    /// `None` only where a statement's points span no direction, which a named plane's cannot.
    fn flipped_by_truth(
        &self,
        survivor: Handle<Surface>,
        points: &PlanePoints,
        sense: Orientation,
    ) -> Option<bool> {
        let Surface::Plane {
            points: theirs,
            sense: their_sense,
            ..
        } = self.surface(survivor)
        else {
            return None;
        };
        let normal = |p: &PlanePoints| {
            let m = self.statement_points(p)?;
            nacre_exact::triple_normal([&m[0], &m[1], &m[2]])
        };
        let agree =
            nacre_exact::same_sense(&normal(points)?, &normal(theirs)?) == (sense == *their_sense);
        Some(!agree)
    }

    /// **A plane statement's three points, exactly, in the frame the statement is written in** —
    /// what its direction `(p₁ − p₀) × (p₂ − p₀)` and so its [`Surface::Plane::sense`] are read
    /// against. `Known` points are themselves; a `Through` statement's are its vertices' meets
    /// ([`Model::through_meets`], any width), in the order the statement holds them.
    ///
    /// `None` where a `Through` statement's meets cannot be placed in one frame.
    pub(super) fn statement_points(&self, p: &PlanePoints) -> Option<[nacre_exact::MeetPoint; 3]> {
        match p {
            PlanePoints::Known(k) => Some(k.map(nacre_exact::MeetPoint::Narrow)),
            PlanePoints::Through(vs) => self.through_meets(*vs),
        }
    }
    /// **The world direction a plane's points span**, before its sense — in `f64`, for the reads
    /// below that compare it with a cache (two normals of one plane, so the sign of their dot is
    /// not a rounding question).
    ///
    /// `None` where this crate cannot carry the points to the world: a chain that does not fold
    /// (a frame, a turn off the quarters — their realization lives in `nacre-ops`), or a
    /// `Through` statement under any motion (its vertices' caches are world points only when
    /// nothing moved the plane after them). Through a folded chain the points' cross travels as
    /// `det M · M·n` — the linear part carries it, and a reflection, carrying points, reverses
    /// it; the offset plays no part, so a chain whose offset overflowed still answers.
    pub(super) fn points_world_direction(
        &self,
        points: &PlanePoints,
        motion: Option<Handle<MotionNode>>,
    ) -> Option<[f64; 3]> {
        let w = match points {
            PlanePoints::Known(p) => {
                let f = |q: [nacre_exact::Rat; 3]| Point3::from_array(q.map(|x| x.to_f64()));
                (f(p[1]) - f(p[0])).cross(f(p[2]) - f(p[0])).as_array()
            }
            PlanePoints::Through(vs) => {
                motion.is_none().then_some(())?;
                let q = vs.map(|v| self.vertex_point(v));
                return Some((q[1] - q[0]).cross(q[2] - q[0]).as_array());
            }
        };
        match motion {
            None => Some(w),
            Some(leaf) => {
                let fold = self.chain_fold(leaf)?;
                let det = f64::from(fold.det());
                Some(fold.dir_f64(w).map(|c| c * det))
            }
        }
    }

    /// Whether the cache a pusher brought faces the way its stated `sense` says — the push door's
    /// cross-check between the figure a producer built (or `nacre-ops` realized) and the sense it
    /// stated. A read only.
    ///
    /// `Known` statements only: a `Through` datum's figure comes from its caller's `f64` cross,
    /// which a nearly collinear triple can turn, so a disagreement there is reachable from input
    /// and is not an assertion's business (a named plane's normal is derived from its name, an
    /// unnamed one's is turned by [`Model::align_cache_sense`]).
    fn cache_agrees_with_sense(
        &self,
        cache: &nacre_geom::Plane,
        points: &PlanePoints,
        motion: Option<Handle<MotionNode>>,
        sense: Orientation,
    ) -> Option<bool> {
        if matches!(points, PlanePoints::Through(_)) {
            return None;
        }
        let w = self.points_world_direction(points, motion)?;
        let n = cache.normal().as_array();
        let d = (0..3).map(|k| w[k] * n[k]).sum::<f64>() * f64::from(sense.sign());
        d.is_finite().then_some(d > 0.0)
    }
    /// **Push a plane stated as the three vertices it passes through** — [`push_plane`]'s twin
    /// for the datum vocabulary, with the same interning contract and the same `flipped` report.
    ///
    /// The name is derived the same way, one step further back: solve each vertex from its three
    /// carriers, then take the canonical form of the plane through those points. So a `Through`
    /// plane and a `Known` plane that *are* the same plane share one handle, which is the whole
    /// point of interning — the variant records how this plane's existence is grounded, not who
    /// asked for it first.
    ///
    /// ★ `vertices` is **sorted** by the caller before it gets here (the same three vertices are
    /// the same statement in any order). Direction is not lost: `sense` states it against the
    /// sorted order, which a caller with its own order reaches through that permutation's parity.
    ///
    /// ★★ **A statement the name key cannot hold still interns — by the statement itself.**
    /// A mixed-frame datum's exact world coefficients are irrational, so
    /// [`Model::plane_name_through`] answers `None`; such a plane takes the second key
    /// (`surface_through_ids`) — the sorted triple and the motion. That is *statement*
    /// identity: the same three vertices under the same motion are one handle, and geometric
    /// identity across different statements is the predicates' to answer per question. This is
    /// **not** the record-less population: the truth (handles + motion) is complete; what does
    /// not exist is a rational description of it.
    ///
    /// ★ The producer remains responsible for rejecting **before** pushing whatever it cannot
    /// frame — this door stores; it does not validate framability.
    ///
    /// [`push_plane`]: Model::push_plane
    pub fn push_plane_through(
        &mut self,
        fallback: nacre_geom::Plane,
        vertices: [Handle<Vertex>; 3],
        motion: Option<Handle<MotionNode>>,
        sense: Orientation,
    ) -> (Handle<Surface>, bool) {
        debug_assert!(
            vertices[0].index() < vertices[1].index() && vertices[1].index() < vertices[2].index(),
            "a Through statement must arrive sorted and duplicate-free"
        );
        let name = self.plane_name_through(vertices);
        if name.is_some() {
            return self.intern_plane(
                fallback,
                name,
                PlanePoints::Through(vertices),
                motion,
                sense,
            );
        }
        if let Some(&h) = self.surface_through_ids.get(&(vertices, motion)) {
            // One statement, one handle: the same sorted vertices span one direction, so the
            // senses alone say whether the survivor faces the other way.
            let flipped =
                matches!(self.surface(h), Surface::Plane { sense: theirs, .. } if *theirs != sense);
            return (h, flipped);
        }
        let h = self.push_plane_raw(
            PlanePoints::Through(vertices),
            motion,
            sense,
            None,
            fallback,
        );
        self.surface_through_ids.insert((vertices, motion), h);
        (h, false)
    }

    /// **The name a `Through` statement derives**, and the one place that derivation lives — the
    /// producer's check and [`Model::push_plane_through`] read the same answer, so "we rejected
    /// what we could not name" is structural rather than two functions agreeing by habit.
    ///
    /// `None` when any vertex is not a three-plane point, when the carriers do not share one
    /// motion (no frame holds a rational coordinate then), or when the three points are
    /// collinear. ★ A meet too wide for `Rat` is **not** on the list: the
    /// name is derived from the meets at whatever width they need
    /// ([`nacre_exact::plane_name_from_meets`]) — width was the arithmetic's problem, never
    /// the statement's.
    pub fn plane_name_through(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<nacre_exact::PlaneName> {
        let m = self.through_meets(vertices)?;
        nacre_exact::plane_name_from_meets([&m[0], &m[1], &m[2]])
    }

    /// **The three vertices' exact meeting points, in the one frame they share** — the single
    /// solve behind [`Model::plane_name_through`] (at push, width-free) and, through the
    /// all-narrow projection [`Model::through_points_rat`], the judging table's witness
    /// triangle. One spelling, so the name a plane interns under and the points a predicate
    /// reasons about cannot describe different planes.
    ///
    /// `None` on any of: a vertex that is not a three-plane point, a vertex [`Model::vertex_meet`]
    /// cannot place in one frame, a carrier with no recorded name, or three vertices whose frames
    /// disagree.
    pub fn through_meets(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<[nacre_exact::MeetPoint; 3]> {
        let mut pts: [Option<nacre_exact::MeetPoint>; 3] = [None, None, None];
        let mut frame = None;
        for (i, vh) in vertices.iter().enumerate() {
            let (p, mine) = self.vertex_meet(*vh)?;
            match frame {
                None if i == 0 => frame = Some(mine),
                f if f == Some(mine) => {}
                _ => return None, // the three vertices do not share one frame
            }
            pts[i] = Some(p);
        }
        Some([pts[0].take()?, pts[1].take()?, pts[2].take()?])
    }

    /// **One vertex's exact meeting point, and the frame it is stated in** — the per-vertex half
    /// of [`Model::through_meets`], which is its caller for the three-vertex case.
    ///
    /// The frame comes back beside the point because a coordinate means nothing without it: `None`
    /// is the world, `Some(node)` a pre-motion frame. A consumer comparing this against anything
    /// stated in world coordinates must **require `None`** — the three-vertex caller above only
    /// needs the three to *agree*, which is a weaker demand and would silently mix frames if it
    /// were copied.
    ///
    /// ★★★ **"Carriers in one frame" is not "carriers with one `motion` field".** A motion that
    /// **fixes** a plane restates nothing, so that plane stays world-stated beside carriers that
    /// moved; its world equation *is* its equation in their pre-motion frame, and
    /// [`Model::chain_fixes_plane`] is what proves it. Reading "no shared motion field" as "no
    /// common frame" costs every turned solid's corners their named datum road.
    ///
    /// `None` on any of: a vertex that is not a three-plane point; carriers carrying **two**
    /// motion histories; a world-stated carrier the shared chain does not fix (or whose name is
    /// `Wide`, since the licence reads narrow coefficients — conservative, and recorded); a
    /// carrier with no recorded name; or three carriers that meet in no point.
    pub fn vertex_meet(
        &self,
        v: Handle<Vertex>,
    ) -> Option<(nacre_exact::MeetPoint, Option<Handle<MotionNode>>)> {
        self.vertex_meet_of(self.vertices.get(v))
    }

    /// [`Model::vertex_meet`] on a definition that has not been pushed yet — what an operation
    /// asks before it states a vertex, so the cache it pushes is already the realization.
    pub fn vertex_meet_of(
        &self,
        def: &Vertex,
    ) -> Option<(nacre_exact::MeetPoint, Option<Handle<MotionNode>>)> {
        let tri = match *def {
            Vertex::ThreePlane(tri) => tri,
            // OnSeam pins a curve, not a point; a Pierce *is* a point but its coordinates
            // are quadratic-irrational — neither has the rational meet a datum statement
            // needs, so both decline here (honest, and spelled per variant so the next
            // variant is a compile error, not a silent fall-through).
            Vertex::OnSeam(_) | Vertex::Pierce { .. } => return None,
        };
        // ★★★ **The third door — solve in the world**. The two below want *one*
        // frame: the world (nothing moved) or one shared chain. A second-generation array breaks
        // both — an array fused in x and then moved in y puts carriers with chains `T2` and
        // `T1·T2` on one corner — and yet every one of those planes states the world exactly,
        // because a rational translation carries a rational name ([`Model::world_plane_name`]).
        // So when the shared-frame roads decline, the triple is solved from the **world names**
        // and the answer carries no leaf: the caller must not replay anything.
        //
        // ★ **It is a fallback, deliberately.** Running it first would re-spell points the two
        // doors already answer, and those spellings are what the corpus is pinned on. Reached
        // only where the two doors answer `None`, it can open a population and cannot move one.
        //
        // ★ It transports the **name**, never the point: `world_plane_name` moves each carrier's
        // equation into the world, and `three_planes_big` then meets three world planes. Nothing
        // here realizes a coordinate and shifts it.
        let world_road = || -> Option<nacre_exact::MeetPoint> {
            let (a, b, c) = (
                self.world_plane_name(tri[0])?,
                self.world_plane_name(tri[1])?,
                self.world_plane_name(tri[2])?,
            );
            nacre_exact::three_planes_big([&a, &b, &c])
        };
        let motions = tri.map(|h| self.plane_motion(h));
        // One shared leaf among the **moved** carriers; no moved carrier means the world.
        let mut leaf = None;
        for m in motions.iter().flatten() {
            match leaf {
                None => leaf = Some(*m),
                Some(l) if l == *m => {}
                // Two histories — no shared frame, so ask whether both state the world.
                Some(_) => return world_road().map(|p| (p, None)),
            }
        }
        let names = tri.map(|h| self.surface_name.get(&h));
        let [Some(a), Some(b), Some(c)] = names else {
            return None;
        };
        // ★★ **A world-stated carrier among chained ones is not a straddle if the chain fixes
        // it.** The invariant-plane restatement mints exactly that shape — a turned block's
        // corner is a restated cap × two chained walls — and a fixed plane's world equation *is*
        // its equation in the pre-motion frame, so the corner solves as if all three shared the
        // chain. Reading the mismatch as a straddle is what cost every turned solid's corners
        // their named datum road; `chain_fixes_plane` is the licence that tells them apart.
        //
        // ★ `narrow()` is asked **only of the carriers being licensed**, never of the other two:
        // this function's answer is a `MeetPoint` of any width, and refusing a `Wide` carrier
        // here would quietly switch off the capability `a_datum_through_wide_meets_keeps_its_name`
        // locks. A `Wide` *fixed* carrier declines the licence — the same conservatism the atoms
        // have, recorded rather than papered over.
        if let Some(leaf) = leaf {
            for (n, m) in [a, b, c].iter().zip(&motions) {
                if m.is_none() && !n.narrow().is_some_and(|c| self.chain_fixes_plane(leaf, c)) {
                    // The carriers straddle frames — the world road is the remaining question.
                    return world_road().map(|p| (p, None));
                }
            }
        }
        Some((nacre_exact::three_planes_big([a, b, c])?, leaf))
    }

    /// [`Model::through_meets`]' all-narrow projection — the form a **witness triangle** takes,
    /// since a witness base is a `[Rat; 3]` by type. `None` additionally when any meet is
    /// [`nacre_exact::MeetPoint::Wide`]; the judging table then builds its witness another way
    /// (the plane's own frame probes), so this is a road fork, not a refusal.
    pub fn through_points_rat(
        &self,
        vertices: [Handle<Vertex>; 3],
    ) -> Option<[[nacre_exact::Rat; 3]; 3]> {
        let meets = self.through_meets(vertices)?;
        let mut pts = [[nacre_exact::Rat::from_int(0); 3]; 3];
        for (o, m) in pts.iter_mut().zip(&meets) {
            *o = *m.narrow()?;
        }
        Some(pts)
    }

    /// The motion a surface's truth records, whichever variant it is.
    #[inline]
    pub fn plane_motion(&self, h: Handle<Surface>) -> Option<Handle<MotionNode>> {
        match self.surface(h) {
            Surface::Plane { motion, .. } | Surface::Cylinder { motion, .. } => *motion,
        }
    }

    /// Push a **cylinder** — the lateral surface, stating its exact truth, with
    /// interning by the whole statement (see [`CylinderKey`] for why the key is deliberately
    /// that literal: merging two `ref_dir`s would split the seam). No `flipped` report — a
    /// literal-identical statement realizes to a literal-identical cache, so there is no other
    /// way round to report.
    pub fn push_cylinder(
        &mut self,
        cache: nacre_geom::Cylinder,
        def: CylinderDef,
        motion: Option<Handle<MotionNode>>,
    ) -> Handle<Surface> {
        let key = (def.clone(), motion);
        if let Some(&h) = self.cylinder_ids.get(&key) {
            return h;
        }
        let h = self.push_cylinder_raw(def, motion, cache);
        self.cylinder_ids.insert(key, h);
        h
    }
}
