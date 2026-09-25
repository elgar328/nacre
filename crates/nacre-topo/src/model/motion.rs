//! Questions asked of a motion chain.

use super::*;
use nacre_exact::LineDir;

/// One step of a motion chain **as a line sees it** — the chain unrolled root first, a frame node
/// opened into its two parts: the frame's `ẑ` onto its plane's normal ([`LineStep::Frame`]), then
/// that plane's own chain. Compared by value: two nodes stating one motion are one step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LineStep<'a> {
    Motion(&'a Motion),
    Frame(Handle<Surface>),
}

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

    /// **Which way a plane's own name faces, against the plane** — `Forward` when the normal of
    /// its name ([`Model::surface_name`]) points the way the plane faces, in the frame its truth is
    /// written in; `Reversed` when it points the other way.
    ///
    /// The name is canonical and so carries no direction (its first nonzero component is
    /// positive); the plane's direction is its truth's — [`Surface::Plane::sense`] against its
    /// points' turn. This compares the two exactly, at any width
    /// ([`nacre_exact::name_along_points`]): the answer is a fact about the truth, never read off
    /// the `f64` cache. `None` for a surface with no name (a mixed-frame `Through` statement, a
    /// cylinder) or a `Through` statement whose meets cannot be placed in one frame.
    pub fn plane_name_sense(&self, surf: Handle<Surface>) -> Option<Orientation> {
        let Surface::Plane { points, sense, .. } = self.surface(surf) else {
            return None;
        };
        let name = self.surface_name.get(&surf)?;
        let m = self.statement_points(points)?;
        let along = nacre_exact::name_along_points(name, [&m[0], &m[1], &m[2]])?;
        Some(if along { *sense } else { sense.flipped() })
    }

    /// **Which way a plane's world name faces, against the plane** — [`Model::plane_name_sense`]
    /// for [`Model::world_plane_name`]: `Forward` when the world name's normal points the way the
    /// plane faces in the world.
    ///
    /// Unmoved, the world name is the name. Moved by a chain that folds (`M` a signed
    /// permutation, `det M = ±1`): the plane's world direction is `det M · M` of its own (a
    /// reflection reverses the points' turn), the name's normal travels as `M·n`, and
    /// canonicalizing that image may negate it — so the answer is the name's own sense times
    /// `det M` times that canonical sign. Every factor is exact, and every one exists wherever the
    /// world name does: `None` exactly where [`Model::world_plane_name`] is `None` (for a named
    /// plane), so a caller that holds the world name always holds its sense.
    pub fn world_plane_name_sense(&self, surf: Handle<Surface>) -> Option<Orientation> {
        let own = self.plane_name_sense(surf)?;
        let Some(leaf) = self.plane_motion(surf) else {
            return Some(own);
        };
        let fold = self.chain_fold(leaf)?;
        let n = self.surface_name.get(&surf)?.narrow()?;
        let carried = fold.dir_rat([n[0], n[1], n[2]])?;
        let world = self.world_plane_name(surf)?;
        let world = world.narrow()?;
        let zero = nacre_exact::Rat::from_int(0);
        let k = (0..3).find(|&k| world[k] != zero)?;
        let canonical_kept = (carried[k] > zero) == (world[k] > zero);
        Some(if canonical_kept == (fold.det() > 0) {
            own
        } else {
            own.flipped()
        })
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
    /// along it — **the same way round**, because a rim circle's edge curve is parameterized about
    /// the cache's axis while the loops that walk it are decided on this statement's — and the
    /// radius must match. A fold in the wrong order disagrees on the origin or the direction.
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
                    && cache.axis().direction().dot(d) > 0.0
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

    /// **The point `p`, stated in the world, carried back to before `leaf`'s chain** — the inverse
    /// of [`Model::chain_point_rat`]: a folded chain permutes the axes with signs and adds a
    /// rational offset, so its inverse is exact. `None` where the chain does not fold or the
    /// arithmetic overflows.
    pub fn chain_point_rat_inverse(
        &self,
        leaf: Handle<MotionNode>,
        p: [nacre_exact::Rat; 3],
    ) -> Option<[nacre_exact::Rat; 3]> {
        self.chain_fold(leaf)?.inverse().point_rat(p)
    }

    /// **The direction `d`, stated in the world, carried back to before `leaf`'s chain** — the
    /// inverse of [`Model::chain_dir_rat`], the linear part only.
    pub fn chain_dir_rat_inverse(
        &self,
        leaf: Handle<MotionNode>,
        d: [nacre_exact::Rat; 3],
    ) -> Option<[nacre_exact::Rat; 3]> {
        self.chain_fold(leaf)?.inverse().dir_rat(d)
    }

    /// **How a plane stands to a cylinder's axis, from their truths** — the relation that decides
    /// an edge's curve on the pair ([`Model::derive_edge_curve`]) and the population gate's arms,
    /// by the one rule ([`nacre_exact::LineDir::relation_to`]).
    ///
    /// ★ **A question about two directions, so only the chains' linear parts matter** — and a
    /// rigid motion keeps the relation. So the two statements need not reach the world, only one
    /// frame: each direction (the plane's name normal, the cylinder's `def.dir()`) is carried out
    /// along its own chain **as far as it goes exactly** — a translation leaves a line alone, a
    /// mirror or a quarter turn permutes it, any turn leaves its own axis alone, a frame takes its
    /// `ẑ` onto its plane's normal — and where the two stop, what is left must be **the same
    /// motions** (compared by value). When one went further, it is carried **back** through the
    /// other's extra steps (the same rules inverted): a plane a turn fixes is read in the turned
    /// cylinder's frame. A plane that is itself a frame's plane stands at `ẑ` just before that
    /// frame's step, which is how a cylinder sketched on a nameless datum still knows its base.
    ///
    /// `None` when the two cannot be placed in one frame — a relation the truth does not state,
    /// which a caller refuses by that name rather than reading the cache.
    pub fn plane_cylinder_relation(
        &self,
        plane: Handle<Surface>,
        cyl: Handle<Surface>,
    ) -> Option<nacre_exact::AxisRelation> {
        let Surface::Cylinder { def, motion: cm } = self.surface(cyl) else {
            return None;
        };
        let axis = LineDir::of_rat(&def.dir());
        let pm = self.plane_motion(plane);
        // Written in one frame already, or both carried whole to the world: the same rule, asked
        // where the walk below would have arrived without walking.
        if pm == *cm {
            if let Some(n) = self.name_normal(plane) {
                return Some(n.relation_to(&axis));
            }
        }
        if let (Some(name), Some(world)) = (
            self.world_plane_name(plane),
            self.world_cylinder_statement(cyl),
        ) {
            return Some(LineDir::normal_of(&name).relation_to(&LineDir::of_rat(&world.dir())));
        }
        let c_path = self.line_path(*cm);
        let (c_dir, c_at) = self.carry_line(axis, &c_path);
        let (p_dir, p_rest) = match c_path.iter().position(|s| *s == LineStep::Frame(plane)) {
            Some(i) => {
                let (d, at) = self.carry_line(LineDir::z(), &c_path[i..]);
                (d, c_path[i + at..].to_vec())
            }
            None => {
                let p_path = self.line_path(pm);
                let (d, at) = self.carry_line(self.name_normal(plane)?, &p_path);
                (d, p_path[at..].to_vec())
            }
        };
        let c_rest = &c_path[c_at..];
        if p_rest == c_rest {
            return Some(p_dir.relation_to(&c_dir));
        }
        // One went further: carry it back through the other's extra steps, last step first.
        let back = |d: LineDir, extra: &[LineStep]| {
            extra
                .iter()
                .rev()
                .try_fold(d, |d, step| self.line_back(&d, step))
        };
        if c_rest.ends_with(&p_rest) {
            let extra = &c_rest[..c_rest.len() - p_rest.len()];
            return Some(back(p_dir, extra)?.relation_to(&c_dir));
        }
        if p_rest.ends_with(c_rest) {
            let extra = &p_rest[..p_rest.len() - c_rest.len()];
            return Some(p_dir.relation_to(&back(c_dir, extra)?));
        }
        None
    }

    /// A plane's name normal, in the frame its truth is written in.
    fn name_normal(&self, plane: Handle<Surface>) -> Option<LineDir> {
        Some(LineDir::normal_of(self.surface_name.get(&plane)?))
    }

    /// `leaf`'s chain as [`LineStep`]s, root first.
    fn line_path(&self, leaf: Option<Handle<MotionNode>>) -> Vec<LineStep<'_>> {
        let mut nodes = Vec::new();
        let mut cur = leaf;
        while let Some(h) = cur {
            let n = self.motion(h);
            nodes.push(&n.motion);
            cur = n.parent;
        }
        let mut out = Vec::new();
        for m in nodes.into_iter().rev() {
            match *m {
                Motion::Frame { plane, .. } => {
                    out.push(LineStep::Frame(plane));
                    out.extend(self.line_path(self.plane_motion(plane)));
                }
                _ => out.push(LineStep::Motion(m)),
            }
        }
        out
    }

    /// Carry the line `d` along `path` as far as it goes exactly: the line reached, and how many
    /// steps it took.
    fn carry_line(&self, d: LineDir, path: &[LineStep]) -> (LineDir, usize) {
        let mut d = d;
        for (i, step) in path.iter().enumerate() {
            match self.line_step(&d, step, false) {
                Some(next) => d = next,
                None => return (d, i),
            }
        }
        (d, path.len())
    }

    /// One step on a line, forward or `back` — a translation leaves it, a mirror or a turn moves
    /// it as [`LineDir`] says, and a frame takes `ẑ` onto its plane's normal (and back). `None`
    /// where the step takes the line somewhere irrational.
    fn line_step(&self, d: &LineDir, step: &LineStep, back: bool) -> Option<LineDir> {
        match *step {
            LineStep::Motion(Motion::Translate { .. }) => Some(d.clone()),
            LineStep::Motion(&Motion::Mirror { axis, .. }) => Some(d.mirrored(axis)),
            LineStep::Motion(&Motion::Rotate { axis, pivot, angle }) => {
                d.turned(nacre_exact::Rotation { axis, pivot, angle }, back)
            }
            // Opened into `LineStep::Frame` and the plane's own chain by `line_path`.
            LineStep::Motion(Motion::Frame { .. }) => None,
            LineStep::Frame(p) if back => {
                let n = self.name_normal(p)?;
                (d.relation_to(&n) == nacre_exact::AxisRelation::Across).then(LineDir::z)
            }
            LineStep::Frame(p) => {
                if d.is_z() {
                    self.name_normal(p)
                } else {
                    None
                }
            }
        }
    }

    /// [`Model::line_step`] backwards.
    fn line_back(&self, d: &LineDir, step: &LineStep) -> Option<LineDir> {
        self.line_step(d, step, true)
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
