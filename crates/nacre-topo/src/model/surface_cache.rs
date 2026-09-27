//! The surface cache is *derived* from the truth: the raw pushes, the derivation, and the
//! counters that measure where it declines.

use super::*;

impl Model {
    /// A plane's push (private): the exact truth and its f64 cache, index-parallel, in one
    /// motion — so the two cannot come apart.
    ///
    /// ★★★★ **The truth comes first because the arena holds it, and the cache is derived from it
    /// here**: this door calls [`Model::apply_derivation`], so wherever the plane has a world name
    /// its cache is the truth's realization — anchor and unit normal — and `fallback`, the
    /// producer's own figure, stands only where the derivation declines (no world name: a chain
    /// that does not fold, a mixed-frame `Through`). Unlike a vertex's cache, a plane's does not
    /// say which of the two it holds; the census counters do (`surface_derive_counts`).
    ///
    /// ★★★ **Two doors split by kind, rather than one taking both enums.** A single
    /// `push_raw(truth: Surface, cache: nacre_geom::Surface)` could be handed a plane truth
    /// beside a cylinder cache, and **four** sites downstream would need an `unreachable!` to
    /// say the pairing holds. A typed door makes the mismatch unspellable, and those four cite
    /// the door instead of asserting the fact.
    ///
    /// ⚠ **No lock here, and this is why.** The proposition is the signature itself, and a
    /// `compile_fail` doc-test cannot reach a private function to demonstrate it. The public
    /// doors ([`Model::push_plane`], [`Model::push_cylinder`]) already took narrow types; what
    /// was wide was this crate-internal one. Recorded rather than locked.
    ///
    /// ★★★ **The name enters here too**, for the same reason the cache does: this is
    /// the one place a plane reaches the arena, so it is the one place that can derive a cache
    /// from the truth — and the derivation reads the name. Inserting it one line later
    /// (at interning) would leave a window in which a surface exists without the name that
    /// describes it.
    pub(super) fn push_plane_raw(
        &mut self,
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
        sense: Orientation,
        name: Option<nacre_exact::PlaneName>,
        fallback: nacre_geom::Plane,
    ) -> Handle<Surface> {
        let h = self.surfaces.push(Surface::Plane {
            points,
            motion,
            sense,
        });
        self.surface_cache.push(SurfaceCache {
            realized: nacre_geom::Surface::Plane(fallback),
        });
        if let Some(n) = name {
            self.surface_name.insert(h, n);
        }
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "the truth and its cache enter together or not at all"
        );
        if !self.apply_derivation(h) {
            self.align_cache_sense(h);
        }
        h
    }

    /// **Turn a plane's cache to face the way its truth says** — the sense is truth, so the cache
    /// follows it (truth → cache, the allowed direction). Returns whether it turned anything.
    ///
    /// For a plane whose normal the door **did not** derive — one with no world name (a
    /// mixed-frame `Through`, a chain that does not fold), whose cache is the producer's figure.
    /// A plane with a world name has nothing to turn: its normal is the name's times the name's
    /// sense ([`Model::derive_surface_cache`]), which is the truth's direction by construction.
    /// Here the `f64` cross of the plane's points, carried to the world where this crate can
    /// carry them, is compared with the cache; chains with a turn or a frame are out of reach —
    /// the census lock (`nacre_ops::audit_plane_senses`) holds them. For a `Known` statement the
    /// push door has already asserted the two agree, so this is the lock that the cache is a
    /// function of the truth rather than a correction.
    pub(crate) fn align_cache_sense(&mut self, h: Handle<Surface>) -> bool {
        let Surface::Plane {
            points,
            motion,
            sense,
        } = self.surface(h)
        else {
            return false;
        };
        let (points, motion, sense) = (points.clone(), *motion, *sense);
        let nacre_geom::Surface::Plane(cache) = *self.surface_cache(h) else {
            return false;
        };
        let n = cache.normal().as_array();
        let Some(w) = self.points_world_direction(&points, motion) else {
            return false;
        };
        let turn = (0..3).map(|k| w[k] * n[k]).sum::<f64>() * f64::from(sense.sign()) < 0.0;
        if turn {
            self.surface_cache[h.index() as usize] = SurfaceCache {
                realized: nacre_geom::Surface::Plane(cache.reversed()),
            };
        }
        turn
    }

    /// **Realize this surface's cache from its truth**, keeping the producer's value where the
    /// truth cannot say ([`Model::derive_surface_cache`] declines). Returns whether the cache
    /// now carries a derived normal — the one fact [`Model::push_plane_raw`] needs to know that
    /// the sense is already the truth's.
    ///
    /// ★ Counting reads the derivation computed here, against the value the producer stated, so
    /// the census `stat` rows keep meaning "how far the producer's value was from the truth's".
    ///
    /// ⚠ Planes only. [`Model::push_cylinder_raw`] calls [`Model::measure_derivation`]: the
    /// cylinder arm is derived and **thrown away** — whether applying it moves any cache bit is
    /// measured for the moved cylinders only (bit-identical), not for the unmoved ones.
    fn apply_derivation(&mut self, h: Handle<Surface>) -> bool {
        let derived = self.derive_surface_cache(h);
        self.count_derivation(h, derived.as_ref());
        match derived {
            Some(realized) => {
                self.surface_cache[h.index() as usize] = SurfaceCache { realized };
                true
            }
            None => false,
        }
    }

    /// A cylinder's push (private) — [`Model::push_plane_raw`]'s twin, and the other half of the
    /// reason neither takes a `nacre_geom::Surface`.
    pub(super) fn push_cylinder_raw(
        &mut self,
        def: CylinderDef,
        motion: Option<Handle<MotionNode>>,
        cache: nacre_geom::Cylinder,
    ) -> Handle<Surface> {
        let h = self.surfaces.push(Surface::Cylinder { def, motion });
        self.surface_cache.push(SurfaceCache {
            realized: nacre_geom::Surface::Cylinder(cache),
        });
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "the truth and its cache enter together or not at all"
        );
        self.measure_derivation(h);
        h
    }

    /// **The f64 cache this surface's truth realizes to** — the one realization road, for surfaces.
    ///
    /// A vertex's cache is its definition realized at birth (`push_vertex_realized`). A plane's is
    /// the same wherever the plane has a world name:
    /// * **Unit normal** — the world name's normal times the name's sense
    ///   ([`Model::world_plane_name_sense`]), each component the `f64` nearest the truth
    ///   ([`nacre_exact::unit_vector_f64`]). Correctly rounded, so it is unique: asking at more
    ///   bits cannot move it, and two statements of one plane cannot disagree about it. The sense
    ///   is multiplied into the integers, so a zero component is `+0.0`.
    /// * **Anchor** — the truth's first point, carried to the world and realized: a `Known`
    ///   statement's through a chain that folds ([`Model::chain_point_rat`]), an unmoved
    ///   `Through` statement's first meet when the meets are rational
    ///   ([`Model::through_points_rat`]). Where neither stands (an overflow, a moved `Through`)
    ///   the producer's anchor is kept beside the derived normal.
    /// * **Cylinder** — the world statement ([`Model::world_cylinder_def`]'s, without its
    ///   postcondition — this reader measures the disagreement it would assert away) realized: `origin` and
    ///   `radius` descend exactly; the axis direction and `ref_dir` do not
    ///   (`Cylinder::from_axis` normalizes both). ⚠ Nothing **applies** this arm today —
    ///   [`Model::push_cylinder_raw`] only measures it.
    ///
    /// ★ **What the anchor buys** (measured): a producer-stated anchor on a tilted plane
    /// varies by up to **22 ulps** with which face asked for the plane first, and for **20 of 29**
    /// moved planes it does not satisfy the plane's own coefficients. The anchor is
    /// a function of the truth, so neither varies with who asked.
    ///
    /// ★★★★★ **Why the anchor and not a canonical row** (measured). A canonical row is the
    /// tidier answer, but it moved 32 census result rows while the judge's exact shortcut still
    /// read the cache's coefficients (a rescaled row failed that road's check against the face
    /// corners — the road now reads the plane's name and not this cache), and it makes
    /// `plane_origin_projection` square the coefficients — which overflows for
    /// 8 of this corpus's planes and spends half the `i128` width budget. The anchor costs none
    /// of that: **0** census rows move, and the arithmetic is one
    /// addition. What it buys is the same thing: the cache becomes **recomputable from the
    /// arena**, so two statements of one plane cannot disagree about where it is anchored.
    ///
    /// ⚠ **This is not «the cache is a function of the geometry».** The anchor is the *first*
    /// point of the *first* pusher's triple; two files whose first pushers state the plane
    /// differently still differ. Inside one model interning makes that unreachable — one name,
    /// one handle, one truth.
    ///
    /// `None` — the caller keeps the cache it has — for a plane with no world name: an unnamed
    /// one, a chain that does not fold, a moved `Wide` name, or an overflow in the carry.
    ///
    /// ⚠★★★ **One door still bypasses this entirely** — [`Model::push_plane_unregistered`], which
    /// skips the name and therefore the derivation. It is the only remaining way for a surface's
    /// truth and its cache to disagree, and it is `cfg(test)`: every call site is a fixture that
    /// wants one geometric plane held as two handles.
    pub(crate) fn derive_surface_cache(&self, h: Handle<Surface>) -> Option<nacre_geom::Surface> {
        let rat3 = |v: [Rat; 3]| [v[0].to_f64(), v[1].to_f64(), v[2].to_f64()];
        match self.surface(h) {
            Surface::Plane { points, motion, .. } => {
                let stated = match self.surface_cache(h) {
                    nacre_geom::Surface::Plane(p) => *p,
                    nacre_geom::Surface::Cylinder(_) => return None,
                };
                let normal = self.world_unit_normal(h)?;
                let anchor = self
                    .truth_anchor(points, *motion)
                    .map_or(stated.origin(), Point3::from_array);
                Some(nacre_geom::Surface::Plane(Plane::from_point_unit_normal(
                    anchor, normal,
                )))
            }
            Surface::Cylinder { .. } => {
                let def = self.world_cylinder_statement(h)?;
                Cylinder::from_axis(
                    Point3::from_array(rat3(def.origin())),
                    Vector3::from_array(rat3(def.dir())),
                    Vector3::from_array(rat3(def.ref_dir())),
                    def.radius_f64(),
                )
                .map(nacre_geom::Surface::Cylinder)
            }
        }
    }

    /// The plane's unit normal realized from its world name and the name's sense — `None` exactly
    /// where the plane has no world name.
    fn world_unit_normal(&self, h: Handle<Surface>) -> Option<Vector3> {
        let name = self.world_plane_name(h)?;
        let sense = self.world_plane_name_sense(h)?;
        let [a, b, c, _] = name.coeff_ints();
        let toward = if sense.sign() < 0 {
            [-a, -b, -c]
        } else {
            [a, b, c]
        };
        nacre_exact::unit_vector_f64(&toward).map(Vector3::from_array)
    }

    /// A plane statement's first point in the world, realized — `None` where the truth does not
    /// place it rationally (a chain that does not fold or overflows, a moved `Through`, meets
    /// wider than `Rat`).
    fn truth_anchor(
        &self,
        points: &PlanePoints,
        motion: Option<Handle<MotionNode>>,
    ) -> Option<[f64; 3]> {
        let first = match (points, motion) {
            (PlanePoints::Known(pts), None) => pts[0],
            (PlanePoints::Known(pts), Some(leaf)) => self.chain_point_rat(leaf, pts[0])?,
            (PlanePoints::Through(vs), None) => self.through_points_rat(*vs)?[0],
            (PlanePoints::Through(_), Some(_)) => return None,
        };
        Some(first.map(|x| x.to_f64()))
    }

    /// Count what [`Model::derive_surface_cache`] would do at this push, and change nothing —
    /// the measuring half of [`Model::apply_derivation`], for the doors that do not apply it.
    fn measure_derivation(&self, h: Handle<Surface>) {
        self.count_derivation(h, self.derive_surface_cache(h).as_ref());
    }

    /// Count one derivation's outcome against the cache the producer stated.
    fn count_derivation(&self, h: Handle<Surface>, derived: Option<&nacre_geom::Surface>) {
        use std::sync::atomic::Ordering::Relaxed;
        match derived {
            None => {
                SURFACE_DECLINED.fetch_add(1, Relaxed);
                self.decline_reason(h).fetch_add(1, Relaxed);
            }
            Some(d) => {
                if surface_bits(d) != surface_bits(self.surface_cache(h)) {
                    SURFACE_DIFFERS.fetch_add(1, Relaxed);
                }
                SURFACE_DERIVED.fetch_add(1, Relaxed);
            }
        };
    }

    /// Which counter a decline belongs to — a **diagnosis of the road already taken**, not a
    /// second derivation. [`Model::derive_surface_cache`] asks [`Model::world_plane_name`],
    /// which folds every cause into one `None`; the populations have to be told apart because
    /// they need different work, and only one of them (`Wide`) is arithmetic at all.
    fn decline_reason(&self, h: Handle<Surface>) -> &'static std::sync::atomic::AtomicU64 {
        match self.surface(h) {
            Surface::Cylinder { .. } => &SURFACE_DECLINED_CYLINDER,
            Surface::Plane { .. } => match self.surface_name.get(&h) {
                None => &SURFACE_DECLINED_UNNAMED,
                Some(n) => match (self.plane_motion(h), n.narrow()) {
                    (_, None) => &SURFACE_DECLINED_WIDE,
                    (Some(leaf), Some(_)) if self.chain_fold(leaf).is_none() => {
                        &SURFACE_DECLINED_MOTION
                    }
                    _ => &SURFACE_DECLINED_ARITH,
                },
            },
        }
    }

    /// Count an interning hit whose incoming statement would have written a **different** cache
    /// from the survivor's — the product-population measurement of "the same geometry, described
    /// twice, writes the same file".
    ///
    /// The survivor's normal is its world name's, which the incoming statement shares (one name,
    /// one handle); only the anchor depends on which statement came first. So where the incoming
    /// truth places its own first point, that is what is compared; where it does not, the
    /// incoming producer's figure — the cache it would have kept.
    pub(super) fn count_discarded_cache(
        &self,
        h: Handle<Surface>,
        points: &PlanePoints,
        motion: Option<Handle<MotionNode>>,
        discarded: &nacre_geom::Plane,
    ) {
        let survivor = self.surface_cache(h);
        let differs = match (self.truth_anchor(points, motion), survivor) {
            (Some(a), nacre_geom::Surface::Plane(p)) => {
                a.map(f64::to_bits) != p.origin().as_array().map(f64::to_bits)
            }
            _ => surface_bits(survivor) != surface_bits(&nacre_geom::Surface::Plane(*discarded)),
        };
        if differs {
            SURFACE_DISCARDED_DIFFERING.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
}
