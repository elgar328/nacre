//! Questions asked of a motion chain.

use super::*;

impl Model {
    /// The motion node a handle names. Same seal as [`Model::surface_cache`]: writing goes through
    /// [`Model::push_motion`] (interned), reading through here.
    #[inline]
    pub fn motion(&self, h: Handle<MotionNode>) -> &MotionNode {
        self.motions.get(h)
    }

    /// Whether more than `n` recorded motions stand between `leaf` and the world — the length of
    /// the chain a replay would walk (a `Frame` node counts once; its own short expansion is not a
    /// history). Walks at most `n + 1` nodes, so asking about a 4,000-deep history costs `n`.
    pub fn motion_deeper_than(&self, leaf: Handle<MotionNode>, n: usize) -> bool {
        let mut depth = 0;
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            depth += 1;
            if depth > n {
                return true;
            }
            cur = self.motion(h).parent;
        }
        false
    }

    /// Whether every node of `leaf`'s recorded chain fixes the plane `coeffs` **as a set** —
    /// the consumer-side twin of the producer's `Isometry::fixes_plane`, on the same scalar
    /// atoms (one rule per motion kind, two thin composers).
    ///
    /// What it licenses: a world-stated plane among motion-carrying carriers is usable in
    /// the carriers' *pre-motion* frame **iff** the chain fixes it — then its world equation
    /// is the same equation there, and the corner solves as if all three shared the chain.
    /// The invariant-plane restatement mints exactly this shape (a turned block's corner =
    /// restated cap × two chained walls), and this check is what keeps the licence honest:
    /// a frame-hosted datum's cap hits the `Frame` arm and the caller stays declined.
    ///
    /// ★ **It lives here, below `nacre-ops`, because [`Model::vertex_meet`] needs it.** Without
    /// the rescue that door reads a turned solid's corner as straddling frames — every carrier
    /// but the fixed cap moved — and the datum road it gates loses its named form. Keeping the
    /// rule in `nacre-ops` would mean a second walk down here, and this rule already cost the
    /// kernel a regression by existing in more than one spelling.
    ///
    /// Conservative by construction — `Frame` nodes are never fixed (their basis is
    /// irrational), and the scalar atoms answer `false` on overflow.
    pub fn chain_fixes_plane(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_exact::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, false)
    }

    /// **A plane's canonical name in the world**, whatever frame its truth is written in — the
    /// door between "how this surface got here" and "where it is".
    ///
    /// Unmoved: the name itself, which is already world. Moved by a chain that folds — translations,
    /// quarter turns and axis reflections, [`Model::chain_plane_coeffs`]: that name carried out
    /// exactly and canonicalized. `None` for anything else — a frame node, a turn off the
    /// quarters, a **moved** [`nacre_exact::PlaneName::Wide`] (no narrow vessel for the transport
    /// to take), or an overflow — and the caller declines rather than guessing. An *unmoved* `Wide` name comes back verbatim: it is already world,
    /// and refusing it would switch off a capability that is locked elsewhere.
    ///
    /// **Planes only.** A cylinder surface has no entry in `surface_name`, so it answers `None`;
    /// its own world statement is `nacre-ops`' `world_cylinder_def`.
    pub fn world_plane_name(&self, surf: Handle<Surface>) -> Option<nacre_exact::PlaneName> {
        let name = self.surface_name.get(&surf)?;
        match self.plane_motion(surf) {
            None => Some(name.clone()),
            Some(leaf) => Some(nacre_exact::PlaneName::Narrow(
                self.chain_plane_coeffs(leaf, *name.narrow()?)?,
            )),
        }
    }

    /// **A cylinder's exact statement in the world** — the one door between a cylinder's truth
    /// (written in the frame its motion names) and every consumer that compares it against world
    /// planes: the population gate's clearance arithmetic, the arrangement's circles and rulings,
    /// the band pass, the derived cache. [`Model::world_plane_name`]'s twin, one door per surface
    /// kind.
    ///
    /// Unmoved: the statement itself. Moved by a chain that folds: origin through
    /// [`Model::chain_point_rat`], axis and seam reference through [`Model::chain_dir_rat`], the
    /// squared radius unchanged — a folded chain is a signed permutation plus a rational offset,
    /// so every invariant [`CylinderDef::new`] checks survives it. Anything else (a frame, a turn
    /// off the quarters, overflow): `None`, and the caller refuses rather than measuring across
    /// two frames.
    ///
    /// ★ **The postcondition is checked, not assumed.** The surface's `f64` cache is already the
    /// *realized* world cylinder, so it is an independent second description of the very thing
    /// this door claims to produce: the origin must lie on the cache's axis, the axis must run
    /// along it, and the radius must match. A fold in the wrong order disagrees on the origin or
    /// the direction.
    pub fn world_cylinder_def(&self, surf: Handle<Surface>) -> Option<CylinderDef> {
        let out = self.world_cylinder_statement(surf)?;
        debug_assert!(
            {
                let nacre_geom::Surface::Cylinder(cache) = self.surface_cache(surf) else {
                    unreachable!("a cylinder truth is pushed beside a cylinder cache")
                };
                let o = Point3::from_array(out.origin().map(|r| r.to_f64()));
                let d = nacre_math::Vector3::from_array(out.dir().map(|r| r.to_f64()));
                let scale = 1.0 + o.as_array().iter().fold(0.0, |m: f64, c| m.max(c.abs()));
                cache.axis().distance(o) <= 1e-9 * scale
                    && cache.axis().direction().cross(d).norm() <= 1e-9 * d.norm()
                    && (cache.radius() - out.radius_f64()).abs() <= 1e-9 * scale
            },
            "the world statement and the realized cache describe one cylinder"
        );
        Some(out)
    }

    /// [`Model::world_cylinder_def`] without its postcondition — for the one reader whose job is
    /// to **measure** how far the cache stands from the truth (the derived cache), where a
    /// disagreement is the measurement, not a defect: a planted lie must reach `validate`.
    pub(super) fn world_cylinder_statement(&self, surf: Handle<Surface>) -> Option<CylinderDef> {
        let Surface::Cylinder { def, motion } = self.surface(surf) else {
            unreachable!("a cylinder surface carries a cylinder truth")
        };
        match motion {
            None => Some(def.clone()),
            Some(leaf) => CylinderDef::new(
                self.chain_point_rat(*leaf, def.origin())?,
                self.chain_dir_rat(*leaf, def.dir())?,
                self.chain_dir_rat(*leaf, def.ref_dir())?,
                def.r2().clone(),
            ),
        }
    }

    /// `leaf`'s whole chain, folded ([`Model::motion_folds`]) — read through the store first, so
    /// a handle from another model dies here as it does at [`Model::motion`].
    pub(super) fn chain_fold(&self, leaf: Handle<MotionNode>) -> Option<&nacre_exact::AxisAffine> {
        let _ = self.motion(leaf);
        self.motion_folds[leaf.index() as usize].as_ref()
    }

    /// **The plane `c`, stated before `leaf`'s chain, restated in the world** — canonicalized, so
    /// a plane reached by two routes lands on one array. `None` where the chain does not fold (a
    /// frame, a turn off the quarters) or the arithmetic overflows `i128`.
    pub fn chain_plane_coeffs(
        &self,
        leaf: Handle<MotionNode>,
        c: [nacre_exact::Rat; 4],
    ) -> Option<[nacre_exact::Rat; 4]> {
        self.chain_fold(leaf)?.plane_coeffs(c)
    }

    /// **The point `p`, stated before `leaf`'s chain, carried to the world** — exactly. `None` as
    /// [`Model::chain_plane_coeffs`].
    pub fn chain_point_rat(
        &self,
        leaf: Handle<MotionNode>,
        p: [nacre_exact::Rat; 3],
    ) -> Option<[nacre_exact::Rat; 3]> {
        self.chain_fold(leaf)?.point_rat(p)
    }

    /// **The direction `d`, stated before `leaf`'s chain, carried to the world** — the linear
    /// part only, so a chain whose offset overflowed still answers.
    pub fn chain_dir_rat(
        &self,
        leaf: Handle<MotionNode>,
        d: [nacre_exact::Rat; 3],
    ) -> Option<[nacre_exact::Rat; 3]> {
        self.chain_fold(leaf)?.dir_rat(d)
    }

    /// The strict twin: every node carries the plane's **coefficient row verbatim**, not merely
    /// the set. The one place they part is a mirror whose plane is the carrier itself (`n ∥ axis`,
    /// on-plane): the set maps to itself but the row comes back negated — harmless to a Cramer
    /// solve (both determinants negate, the ratio stands), fatal to a determinant *sign* read.
    /// So [`Model::chain_fixes_plane`] licenses solving, and this licenses judging-table
    /// descriptions (`nacre_ops`' mirror re-chaining).
    pub fn chain_preserves_plane_row(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_exact::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, true)
    }

    fn chain_fixes(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_exact::Rat; 4],
        rows_verbatim: bool,
    ) -> bool {
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            let n: &MotionNode = self.motion(h);
            let fixed = match &n.motion {
                Motion::Rotate { axis, .. } => {
                    nacre_exact::axis_rotation_fixes_plane(*axis, coeffs)
                }
                Motion::Translate { offset } => {
                    nacre_exact::translation_fixes_plane(offset, coeffs)
                }
                Motion::Mirror { axis, offset } => {
                    if rows_verbatim {
                        coeffs[match axis {
                            nacre_exact::Axis::X => 0,
                            nacre_exact::Axis::Y => 1,
                            nacre_exact::Axis::Z => 2,
                        }] == nacre_exact::Rat::from_int(0)
                    } else {
                        nacre_exact::mirror_fixes_plane(*axis, *offset, coeffs)
                    }
                }
                Motion::Frame { .. } => false,
            };
            if !fixed {
                return false;
            }
            cur = n.parent;
        }
        true
    }
}
