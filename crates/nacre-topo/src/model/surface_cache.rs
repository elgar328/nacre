//! The surface cache is *derived* from the truth: the raw pushes, the derivation, and the
//! counters that measure where it declines.

use super::*;

impl Model {
    /// A plane's push (private): the exact truth and its f64 cache, index-parallel, in one
    /// motion — so the two cannot come apart.
    ///
    /// ★★★★ **The truth comes first because the arena holds it, and the cache is derived from it
    /// here**: this door calls [`Model::apply_derivation`], so what the producer hands
    /// in survives only in the parts the truth does not decide — the row (`raw`) and with it the
    /// sense — and wholesale where [`Model::derive_surface_cache`] declines.
    /// ⚠ Only the anchor is derived; «the door takes only the truth» is **not** reached
    /// while `cache` is still a parameter.
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
        cache: nacre_geom::Plane,
    ) -> Handle<Surface> {
        let h = self.surfaces.push(Surface::Plane {
            points,
            motion,
            sense,
        });
        self.surface_cache.push(SurfaceCache {
            realized: nacre_geom::Surface::Plane(cache),
        });
        if let Some(n) = name {
            self.surface_name.insert(h, n);
        }
        debug_assert_eq!(
            self.surface_cache.len(),
            self.surfaces.len(),
            "the truth and its cache enter together or not at all"
        );
        self.apply_derivation(h);
        h
    }

    /// **Realize this surface's cache from its truth**, keeping the producer's value where the
    /// truth cannot say ([`Model::derive_surface_cache`] declines).
    ///
    /// ★ Counting happens **first**, against the value the producer stated, so the census `stat`
    /// rows keep meaning "how far the producer's value was from the truth's".
    ///
    /// ⚠ Planes only. [`Model::push_cylinder_raw`] still calls [`Model::measure_derivation`]:
    /// the cylinder arm is derived and **thrown away**, deliberately, because a moved cylinder's
    /// world statement lives in `nacre-ops` and this door cannot reach it.
    fn apply_derivation(&mut self, h: Handle<Surface>) {
        self.measure_derivation(h);
        if let Some(realized) = self.derive_surface_cache(h) {
            self.surface_cache[h.index() as usize] = SurfaceCache { realized };
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
    /// A vertex has a whole one (`push_vertex_realized` realizes the definition at
    /// birth). A surface has **half** of one: the anchor is derived here, the row and
    /// the sense are whatever the producer handed in.
    ///
    /// ★ **What that half buys** (measured): a producer-stated anchor on a tilted plane
    /// varies by up to **22 ulps** with which face asked for the plane first, and for **20 of 29**
    /// moved planes it does not satisfy the plane's own coefficients. The anchor is
    /// a function of the truth, so neither varies with who asked.
    ///
    /// What it derives:
    /// * **Plane** — the **anchor**, and nothing else: the truth's first point, carried to the
    ///   world and realized. The row (`raw`) and with it the sense are copied from the value the
    ///   producer stated.
    /// * **Cylinder** — `origin` and `radius` descend exactly from [`CylinderDef`]; the axis
    ///   direction and `ref_dir` do not (`Cylinder::from_axis` normalizes both). ⚠ Nothing
    ///   **applies** this arm today — [`Model::push_cylinder_raw`] only measures it.
    ///
    /// ★★★★★ **Why the anchor and not the row** (measured). A canonical row is the
    /// tidier answer, but it moves 32 census result rows, costs an exact-coefficient
    /// road (`WorkingPlane::reconcile` gates on `Plane::spans_exactly`, which a rescaled row
    /// fails), and makes `plane_origin_projection` square the coefficients — which overflows for
    /// 8 of this corpus's planes and spends half the `i128` width budget. The anchor costs none
    /// of that: **0** census rows move, the exact road is untouched, and the arithmetic is one
    /// addition. What it buys is the same thing: the cache becomes **recomputable from the
    /// arena**, so two statements of one plane cannot disagree about where it is anchored.
    ///
    /// ⚠ **This is not «the cache is a function of the geometry».** The anchor is the *first*
    /// point of the *first* pusher's triple; two files whose first pushers state the plane
    /// differently still differ. Inside one model interning makes that unreachable — one name,
    /// one handle, one truth.
    ///
    /// ☑ **The sense cannot come from the truth** (measured): 342 non-seed `Known` planes carry a
    /// cache normal opposing their own point order, so the point order does not name a direction.
    /// Copying the row sidesteps the question — `flipped` and every face's outward spelling are
    /// bit-unchanged by this derivation.
    ///
    /// `None` — the caller keeps the cache it has — for an unnamed plane, a `Wide` name, a
    /// motion that is not a rational translation chain, a `Through` truth, a moved cylinder, or
    /// an overflow. ⚠ The name is still required even though the anchor does not read it: every
    /// number above was measured with that gate on, and widening it is its own measurement.
    ///
    /// ⚠★★★ **One door still bypasses this entirely** — [`Model::push_plane_unregistered`], which
    /// skips the name and therefore the derivation. It is the only remaining way for a surface's
    /// truth and its cache to disagree about where the plane is anchored, and it is `cfg(test)`:
    /// every call site is a fixture that wants one geometric plane held as two handles.
    pub(crate) fn derive_surface_cache(&self, h: Handle<Surface>) -> Option<nacre_geom::Surface> {
        let rat3 = |v: [Rat; 3]| [v[0].to_f64(), v[1].to_f64(), v[2].to_f64()];
        match self.surface(h) {
            Surface::Plane { points, motion, .. } => {
                // The gate, kept verbatim: a plane the model cannot name in the world is one this
                // derivation declines, and every measured number above assumes that population.
                self.world_plane_name(h)?.narrow()?;
                let stated = match self.surface_cache(h) {
                    nacre_geom::Surface::Plane(p) => *p,
                    nacre_geom::Surface::Cylinder(_) => return None,
                };
                let PlanePoints::Known(pts) = points else {
                    // A `Through` truth names vertices, whose meet may not fit `Rat` at all.
                    return None;
                };
                let t = match motion {
                    None => [Rat::from_int(0); 3],
                    Some(leaf) => self.chain_translation(*leaf)?,
                };
                let mut anchor = [0.0f64; 3];
                for (k, a) in anchor.iter_mut().enumerate() {
                    *a = pts[0][k].checked_add(t[k])?.to_f64();
                }
                // The row verbatim — `coefficients()` is `[raw, −raw·origin]`, so its first three
                // are the `raw` the producer built, and copying them keeps the sense with it.
                let c = stated.coefficients();
                let raw = Vector3::from_array([c[0], c[1], c[2]]);
                Plane::from_point_normal(Point3::from_array(anchor), raw)
                    .map(nacre_geom::Surface::Plane)
            }
            Surface::Cylinder { def, motion } => {
                motion.is_none().then_some(())?;
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

    /// Count what [`Model::derive_surface_cache`] would do at this push, and change nothing —
    /// the measuring half of [`Model::apply_derivation`], for the doors that do not apply it.
    fn measure_derivation(&self, h: Handle<Surface>) {
        use std::sync::atomic::Ordering::Relaxed;
        match self.derive_surface_cache(h) {
            None => {
                SURFACE_DECLINED.fetch_add(1, Relaxed);
                self.decline_reason(h).fetch_add(1, Relaxed);
            }
            Some(d) => {
                if surface_bits(&d) != surface_bits(self.surface_cache(h)) {
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
                    (Some(leaf), Some(_)) if self.chain_translation(leaf).is_none() => {
                        &SURFACE_DECLINED_MOTION
                    }
                    _ => &SURFACE_DECLINED_ARITH,
                },
            },
        }
    }

    /// Count an interning hit whose incoming cache is **discarded** in favour of the survivor's,
    /// when the two are not the same bits — the product-population measurement of "the same
    /// geometry, described twice, writes the same file".
    pub(super) fn count_discarded_cache(&self, h: Handle<Surface>, discarded: &nacre_geom::Plane) {
        if surface_bits(self.surface_cache(h))
            != surface_bits(&nacre_geom::Surface::Plane(*discarded))
        {
            SURFACE_DISCARDED_DIFFERING.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
}
