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

    /// Whether a replay through `leaf`'s chain walks more than `n` nodes
    /// ([`Model::chain_depth_up_to`]). Each step hands the budget that is left down to what it
    /// unfolds into, and the walk stops once past `n`, so asking about a 4,000-deep history costs
    /// about `n` steps.
    pub fn chain_deeper_than(&self, leaf: Handle<MotionNode>, n: usize) -> bool {
        self.chain_depth_up_to(leaf, n) > n
    }

    /// Whether the chain ending at `leaf` **continues** the one ending at `prefix` — `prefix` is
    /// `leaf` or one of its ancestors, and the world (`None`) is continued by every chain. Nodes are
    /// interned by `(motion, parent)`, so this is chain-prefix equality, asked of handles. Walks
    /// `leaf`'s ancestors only when the two differ.
    pub fn chain_continues(
        &self,
        leaf: Option<Handle<MotionNode>>,
        prefix: Option<Handle<MotionNode>>,
    ) -> bool {
        let Some(prefix) = prefix else {
            return true;
        };
        let mut cur = leaf;
        while let Some(h) = cur {
            if h == prefix {
                return true;
            }
            cur = self.motion(h).parent;
        }
        false
    }

    /// **How many nodes a replay through `leaf`'s chain walks**, counted until the count passes `n` —
    /// the same length [`Model::chain_deeper_than`] asks about, as a number, for a caller that sums
    /// the replays several chains will cost. Past `n` the number says only that: a frame step adds
    /// what its plane unfolds to at once, so it can overshoot.
    ///
    /// **A frame node counts what it unfolds to.** It names a plane, and a replay runs the frame
    /// and then that plane's own chain (`nacre_ops`' `frame_chain`), so the node is one step plus
    /// [`Model::plane_replay_depth`]. Counted once, a sketch on a face of a 7°-and-200-translations
    /// box replayed 201 nodes behind a chain of depth one. The walk stops once past `n`, through the
    /// frames too.
    pub fn chain_depth_up_to(&self, leaf: Handle<MotionNode>, n: usize) -> usize {
        let mut depth = 0;
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            depth += 1;
            if let Motion::Frame { plane, .. } = self.motion(h).motion {
                depth += self.plane_replay_depth(plane, n.saturating_sub(depth));
            }
            if depth > n {
                break;
            }
            cur = self.motion(h).parent;
        }
        depth
    }

    /// **How many nodes a replay walks to place `plane`'s statement in the world**, counted no
    /// further than past `n`: its own chain, and for a nameless `Through` statement what its witness
    /// replays for each of its three vertices too, summed — the walk `rotated_vertex::through_judged_points`
    /// takes. A vertex whose names meet in a narrow point in one frame is one point through that
    /// frame's chain; any other (it straddles frames, or its meet is wider than a witness base) is a
    /// meet of its carriers' witness triangles, three points through each carrier's replay — the
    /// count the mixed vertex road's cap keeps too. (A nameless statement's carriers are named: the
    /// datum door refuses one on a nameless carrier.)
    pub fn plane_replay_depth(&self, plane: Handle<Surface>, n: usize) -> usize {
        let own = self
            .surface_motion(plane)
            .map_or(0, |m| self.chain_depth_up_to(m, n));
        let Surface::Plane {
            points: PlanePoints::Through(vs),
            ..
        } = self.surface(plane)
        else {
            return own;
        };
        if own > n || self.surface_name(plane).is_some() {
            return own;
        }
        let mut depth = own;
        for &v in vs {
            if depth > n {
                break;
            }
            let def = self.vertex(v);
            depth += match self.vertex_named_meet(def) {
                Some(crate::NamedMeet::Point(nacre_exact::MeetPoint::Narrow(_), frame)) => {
                    frame.map_or(0, |m| self.chain_depth_up_to(m, n - depth))
                }
                // Straddling, wide, or naming no point: a meet of its carriers' witness triangles.
                _ => def
                    .carriers()
                    .map(|c| 3 * self.plane_replay_depth(c, n - depth))
                    .sum(),
            };
        }
        depth
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
    /// **It lives here, below `nacre-ops`, because [`Model::vertex_meet`] needs it.** Without
    /// the rescue that door reads a turned solid's corner as straddling frames — every carrier
    /// but the fixed cap moved — and the datum road it gates loses its named form. Keeping the rule
    /// in `nacre-ops` would mean a second walk down here — one rule in two spellings, which drift.
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

    /// **The line a plane's chain turns it about** — its hinge, in the world: `(point, dir)`, both
    /// rational.
    ///
    /// A plane with no world statement sits on a chain with a turn off the quarters, and such a
    /// turn leaves exactly one line in place: its own axis line, through the pivot along the
    /// axis. When the chain's first node is such a turn about a line on the plane's name, and every
    /// later turn off the quarters is about the same line (as the exact nodes between have carried
    /// it), every node leaves the line's points in place or carries them exactly — so the plane
    /// holds that line wherever the chain puts it. It comes out carried by the exact nodes after
    /// the turn.
    ///
    /// `None` where that is not shown: no turn off the quarters (the chain folds, and the plane has
    /// a world statement instead), a frame node anywhere on the chain, an exact first node (the
    /// transform carries an exact motion into a plane's statement instead while the plane has no
    /// history), two turns about different lines, a line off the name, an overflow. Nothing is read from an
    /// approximation.
    pub fn plane_hinge(
        &self,
        plane: Handle<Surface>,
    ) -> Option<([nacre_exact::Rat; 3], [nacre_exact::Rat; 3])> {
        use nacre_exact::Rat;
        let name = self.surface_name(plane)?;
        let mut nodes = Vec::new();
        let mut cur = self.surface_motion(plane);
        while let Some(h) = cur {
            nodes.push(self.motion(h).motion);
            cur = self.motion(h).parent;
        }
        let zero = Rat::from_int(0);
        let axis_line = |motion: Motion| -> Option<([Rat; 3], [Rat; 3])> {
            let Motion::Rotate { axis, pivot, .. } = motion else {
                return None;
            };
            let mut dir = [zero; 3];
            dir[match axis {
                nacre_exact::Axis::X => 0,
                nacre_exact::Axis::Y => 1,
                nacre_exact::Axis::Z => 2,
            }] = Rat::from_int(1);
            Some((pivot, dir))
        };
        let mut nodes = nodes.into_iter().rev();
        let first = nodes.next()?;
        if first.exact_affine().is_some() {
            return None;
        }
        let (mut p, mut d) = axis_line(first)?;
        let n = name.narrow()?;
        if nacre_exact::dot_sign_rat(&[n[0], n[1], n[2]], &d) != nacre_exact::Orient::Zero
            || nacre_exact::plane_residual_sign(name, p) != 0
        {
            return None;
        }
        for motion in nodes {
            if let Some(a) = motion.exact_affine() {
                (p, d) = (a.point_rat(p)?, a.dir_rat(d)?);
                continue;
            }
            // The same line: parallel, and the pivot on it.
            let (q, e) = axis_line(motion)?;
            let off = [
                q[0].checked_sub(p[0])?,
                q[1].checked_sub(p[1])?,
                q[2].checked_sub(p[2])?,
            ];
            let on_line = |v: [Rat; 3]| {
                nacre_exact::cross3_rat(&v, &d).is_some_and(|c| c.iter().all(|x| *x == zero))
            };
            if !on_line(e) || !on_line(off) {
                return None;
            }
        }
        Some((p, d))
    }

    /// **The cylinder `def`, stated in the world, carried back into the frame `plane`'s statement is
    /// written in** — where the plane's name and the cylinder are both rational, so a relation
    /// between them is asked exactly.
    ///
    /// A plane with no world statement sits on a chain with a turn off the quarters, and such a turn
    /// has no rational inverse. It needs none when it turns about the cylinder's own axis line: the
    /// turn then leaves the cylinder in place as a set, and the plane stands to it before the turn
    /// as it does after. So the chain is walked from its leaf back to its root — an exact node
    /// (translation, quarter turn, axis mirror) carries the cylinder back exactly, its squared
    /// radius unchanged under a rigid motion; a turn off the quarters must turn about the cylinder's
    /// axis line as it stands in that node's frame (parallel to the axis, its pivot on it).
    ///
    /// `None` where that is not shown: a turn about another line, a frame node, an overflow.
    /// Nothing is read from an approximation.
    pub fn cylinder_before_chain(
        &self,
        plane: Handle<Surface>,
        def: &CylinderDef,
    ) -> Option<CylinderDef> {
        use nacre_exact::Rat;
        let zero = Rat::from_int(0);
        let mut out = def.clone();
        let mut cur = self.surface_motion(plane);
        while let Some(h) = cur {
            let node = self.motion(h);
            match node.motion.exact_affine() {
                Some(a) => {
                    let back = a.inverse();
                    out = CylinderDef::new(
                        back.point_rat(out.origin())?,
                        back.dir_rat(out.dir())?,
                        back.dir_rat(out.ref_dir())?,
                        out.r2().clone(),
                    )?;
                }
                None => {
                    let Motion::Rotate { axis, pivot, .. } = node.motion else {
                        return None;
                    };
                    let mut along = [zero; 3];
                    along[match axis {
                        nacre_exact::Axis::X => 0,
                        nacre_exact::Axis::Y => 1,
                        nacre_exact::Axis::Z => 2,
                    }] = Rat::from_int(1);
                    let o = out.origin();
                    let off = [
                        pivot[0].checked_sub(o[0])?,
                        pivot[1].checked_sub(o[1])?,
                        pivot[2].checked_sub(o[2])?,
                    ];
                    let on_axis = |v: [Rat; 3]| {
                        nacre_exact::cross3_rat(&v, &out.dir())
                            .is_some_and(|c| c.iter().all(|x| *x == zero))
                    };
                    if !on_axis(along) || !on_axis(off) {
                        return None;
                    }
                }
            }
            cur = node.parent;
        }
        Some(out)
    }

    /// **A plane's canonical name in the world**, whatever frame its truth is written in — the
    /// door between "how this surface got here" and "where it is".
    ///
    /// Unmoved: the name itself, which is already world. Moved: that name carried out by the
    /// chain's fold ([`nacre_exact::ChainFold::plane_name`]) and canonicalized, wherever the fold carries it
    /// to the rationals — every plane under translations, quarter turns and axis reflections, and
    /// under a frame the planes its axis lengths allow (the far cap of a frame whose normal has a
    /// rational length, a wall through its origin). No name for anything else — a turn off the
    /// quarters, a plane a frame carries off the rationals, a **moved**
    /// [`nacre_exact::PlaneName::Wide`] (`Width` even where its image is narrow — see
    /// `moved_world_name`), or an overflow — and the caller declines rather than guessing, by the cause the door names
    /// ([`WorldFail`]: no name, no fold or an irrational image is `Unstated`, which is asked before
    /// a moved `Wide` name's `Width`). An *unmoved* `Wide` name comes back verbatim: it is already world,
    /// and refusing it would switch off a capability that is locked elsewhere.
    ///
    /// **Planes only.** A cylinder has no name ([`Model::surface_name`]), so it answers
    /// `Unstated`; its own world statement is [`Model::world_cylinder_def`].
    pub fn world_plane_name(
        &self,
        surf: Handle<Surface>,
    ) -> Result<nacre_exact::PlaneName, WorldFail> {
        let name = self.surface_name(surf).ok_or(WorldFail::Unstated)?;
        match self.surface_motion(surf) {
            None => Ok(name.clone()),
            Some(_) => Ok(self.moved_world_name(surf, name)?.0),
        }
    }

    /// **A moved plane's world name and which way it faces against its carried name** — one
    /// carry ([`nacre_exact::ChainFold::plane_name`]) answers both, so [`Model::world_plane_name`]
    /// and [`Model::world_plane_name_sense`] cannot disagree on whether a name exists. The image's
    /// geometry is asked first, a wide name's too; a moved `Wide` name is then `Width` even where
    /// its image is narrow — whether such a plane takes a world name is decided with the cylinder
    /// road's wide coefficients (todo), since naming it would move bodies that build.
    fn moved_world_name(
        &self,
        surf: Handle<Surface>,
        name: &nacre_exact::PlaneName,
    ) -> Result<(nacre_exact::PlaneName, i8), WorldFail> {
        let leaf = self.surface_motion(surf).ok_or(WorldFail::Unstated)?;
        let (image, turn) = self.chain_fold(leaf)?.plane_name(name)?;
        match (name, &image) {
            (nacre_exact::PlaneName::Narrow(_), nacre_exact::PlaneName::Narrow(_)) => {
                Ok((image, turn))
            }
            _ => Err(WorldFail::Width),
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
        let name = self.surface_name(surf)?;
        let m = self.statement_points(points)?;
        let along = nacre_exact::name_along_points(name, [&m[0], &m[1], &m[2]])?;
        Some(if along { *sense } else { sense.flipped() })
    }

    /// **Which way a plane's world name faces, against the plane** — [`Model::plane_name_sense`]
    /// for [`Model::world_plane_name`]: `Forward` when the world name's normal points the way the
    /// plane faces in the world.
    ///
    /// Unmoved, the world name is the name. Moved by a folded chain (a rigid map with linear part
    /// `M`, `det M = ±1`): the plane's world direction is `det M · M` of its own (a reflection
    /// reverses the points' turn), the name's normal travels along `M`, and canonicalizing that
    /// image may negate it — so the answer is the name's own sense times `det M` times that
    /// canonical sign, which the carry that makes the world name reports beside it. Every factor
    /// is exact, and every one exists wherever the world name does: `None` exactly where
    /// [`Model::world_plane_name`] has no name (for a named plane), so a caller that holds the
    /// world name always holds its sense.
    pub fn world_plane_name_sense(&self, surf: Handle<Surface>) -> Option<Orientation> {
        let own = self.plane_name_sense(surf)?;
        let Some(leaf) = self.surface_motion(surf) else {
            return Some(own);
        };
        let (_, turn) = self.moved_world_name(surf, self.surface_name(surf)?).ok()?;
        let det = self.chain_fold(leaf).ok()?.det();
        Some(if (turn > 0) == (det > 0) {
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
    /// Unmoved: the statement itself. Moved: carried by the chain's fold
    /// ([`nacre_exact::ChainFold::cylinder`]) — the origin as a point, the axis and the seam
    /// reference as lines (they are raw directions), the squared radius unchanged — a rigid map
    /// keeps every invariant [`CylinderDef::new`] checks. Where the fold carries a part off the
    /// rationals (a turn off the quarters, a circle off a frame's origin along an irrational axis
    /// length) or past `Rat`: the cause ([`WorldFail`]), every part's geometry asked before any
    /// width, and the caller refuses rather than measuring across two frames.
    ///
    /// **The postcondition is checked, not assumed.** The surface's `f64` cache is already the
    /// *realized* world cylinder, so it is an independent second description of the very thing
    /// this door claims to produce: the origin must lie on the cache's axis, the axis must run
    /// along it — **the same way round**, because a rim circle's edge curve is parameterized about
    /// the cache's axis while the loops that walk it are decided on this statement's — and the
    /// radius must match. A fold in the wrong order disagrees on the origin or the direction.
    pub fn world_cylinder_def(&self, surf: Handle<Surface>) -> Result<CylinderDef, WorldFail> {
        let out = self.world_cylinder_def_unchecked(surf)?;
        debug_assert!(
            {
                let Some(cache) = self.surface_cache(surf).cylinder() else {
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
        Ok(out)
    }

    /// [`Model::world_cylinder_def`] without its postcondition — for the one reader whose job is
    /// to **measure** how far the cache stands from the truth (the derived cache), where a
    /// disagreement is the measurement, not a defect: a planted lie must reach `validate`.
    pub(super) fn world_cylinder_def_unchecked(
        &self,
        surf: Handle<Surface>,
    ) -> Result<CylinderDef, WorldFail> {
        let Surface::Cylinder { def, motion } = self.surface(surf) else {
            unreachable!("a cylinder surface carries a cylinder truth")
        };
        match motion {
            None => Ok(def.clone()),
            Some(leaf) => Ok(self.chain_fold(*leaf)?.cylinder(def)?),
        }
    }

    /// `leaf`'s whole chain, folded — read through the store first, so a handle from another model
    /// dies here as it does at [`Model::motion`]. `Unstated` where the chain has no fold (a turn off
    /// the quarters, a frame with no fold of its own ([`Model::frame_fold`]), a frame the fold
    /// cannot compose with what carries it).
    pub(super) fn chain_fold(
        &self,
        leaf: Handle<MotionNode>,
    ) -> Result<&nacre_exact::ChainFold, WorldFail> {
        let _ = self.motion(leaf);
        self.chain_affines[leaf.index() as usize]
            .fold
            .as_ref()
            .ok_or(WorldFail::Unstated)
    }

    /// **What a frame node stating `placement` on `plane` (turned by `flip`) folds to**, without
    /// pushing one — the frame it states ([`FramePlacement::frame_on`], the one the judge replays),
    /// carried by the plane's own chain: coordinates read in that frame, carried to the world.
    /// `None` for a plane with no narrow name, a frame the narrow derivation cannot state, or a
    /// plane chain the frame does not compose with ([`nacre_exact::ChainFold::then`]) — every
    /// statement on such a chain is then `Unstated`.
    ///
    /// [`Model::push_motion`] folds a frame node by it once, as the node is born; a producer asks it
    /// to learn whether a frame's basis is rational in the world before deciding to push a node at
    /// all. The plane's name is read now: a name planted on the plane later (`plant_surface_name`,
    /// a test door) leaves a born node's fold describing the plane as it was.
    pub fn frame_fold(
        &self,
        plane: Handle<Surface>,
        placement: FramePlacement,
        flip: bool,
    ) -> Option<nacre_exact::ChainFold> {
        let name = self.surface_name(plane)?.narrow()?;
        let raw = nacre_exact::ChainFold::frame(&placement.frame_on(*name, flip)?);
        match self.surface_motion(plane) {
            None => Some(raw),
            Some(m) => raw.then(self.chain_fold(m).ok()?),
        }
    }

    /// **Whether `leaf`'s chain holds a frame node** — read as the node was born, whether or not
    /// the chain folds. The sketch frame of a face asks it: a plane carried through a frame keeps
    /// that frame's convention even where it has a world name.
    pub fn chain_through_frame(&self, leaf: Handle<MotionNode>) -> bool {
        let _ = self.motion(leaf);
        self.chain_affines[leaf.index() as usize].through_frame
    }

    /// **The plane `c`, stated before `leaf`'s chain, restated in the world** — canonicalized, so
    /// a plane reached by two routes lands on one array. [`WorldFail::Unstated`] where the chain
    /// has no fold or carries the plane off the rationals, `Width` where `i128` overflows.
    pub fn chain_plane_coeffs(
        &self,
        leaf: Handle<MotionNode>,
        c: [nacre_exact::Rat; 4],
    ) -> Result<[nacre_exact::Rat; 4], WorldFail> {
        Ok(self.chain_fold(leaf)?.plane_coeffs(c)?)
    }

    /// **The point `p`, stated before `leaf`'s chain, carried to the world** — exactly, or why not,
    /// as [`Model::chain_plane_coeffs`].
    pub fn chain_point_rat(
        &self,
        leaf: Handle<MotionNode>,
        p: [nacre_exact::Rat; 3],
    ) -> Result<[nacre_exact::Rat; 3], WorldFail> {
        Ok(self.chain_fold(leaf)?.point_rat(p)?)
    }

    /// **The vector `d`, stated before `leaf`'s chain, carried to the world** — its length kept,
    /// for a caller that measures along it; the linear part only, so a chain whose offset
    /// overflowed still answers.
    pub fn chain_dir_rat(
        &self,
        leaf: Handle<MotionNode>,
        d: [nacre_exact::Rat; 3],
    ) -> Result<[nacre_exact::Rat; 3], WorldFail> {
        Ok(self.chain_fold(leaf)?.dir_rat(d)?)
    }

    /// **The direction `d`, stated before `leaf`'s chain, carried to the world as a line** — some
    /// positive multiple of its image: where it points, for a caller that asks only that (a line's
    /// meet with a plane does not change with the length of its direction).
    pub fn chain_dir_along(
        &self,
        leaf: Handle<MotionNode>,
        d: [nacre_exact::Rat; 3],
    ) -> Result<[nacre_exact::Rat; 3], WorldFail> {
        Ok(self.chain_fold(leaf)?.dir_along(d)?)
    }

    /// **The point `p`, stated in the world, carried back to before `leaf`'s chain** — the inverse
    /// of [`Model::chain_point_rat`], exact or why not.
    pub fn chain_point_rat_inverse(
        &self,
        leaf: Handle<MotionNode>,
        p: [nacre_exact::Rat; 3],
    ) -> Result<[nacre_exact::Rat; 3], WorldFail> {
        Ok(self.chain_fold(leaf)?.point_rat_inverse(p)?)
    }

    /// **The direction `d`, stated in the world, carried back to before `leaf`'s chain** — some
    /// positive multiple of it: where it points, for a caller that states a direction by its line
    /// (a frame's `+u`).
    pub fn chain_dir_along_inverse(
        &self,
        leaf: Handle<MotionNode>,
        d: [nacre_exact::Rat; 3],
    ) -> Result<[nacre_exact::Rat; 3], WorldFail> {
        Ok(self.chain_fold(leaf)?.dir_along_inverse(d)?)
    }

    /// **How a plane stands to a cylinder's axis, from their truths** — the relation that decides
    /// an edge's curve on the pair ([`Model::derive_edge_curve`]) and the population gate's arms,
    /// by the one rule ([`nacre_exact::LineDir::relation_to`]).
    ///
    /// **A question about two directions, so only the chains' linear parts matter** — and a
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
        let pm = self.surface_motion(plane);
        // Written in one frame already, or both carried whole to the world: the same rule, asked
        // where the walk below would have arrived without walking.
        if pm == *cm {
            if let Some(n) = self.name_normal(plane) {
                return Some(n.relation_to(&axis));
            }
        }
        if let (Ok(name), Ok(world)) = (
            self.world_plane_name(plane),
            self.world_cylinder_def_unchecked(cyl),
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

    /// **Are two cylinders' axes parallel, as the truth states them** — `Some(true)`: where the two
    /// laterals meet, they meet on rulings.
    ///
    /// Asked in the frame the two statements share — one motion, or both carried whole to the
    /// world (a rigid motion keeps the relation). `None` when neither holds: a relation the truth
    /// does not state.
    pub fn cylinder_pair_relation(&self, a: Handle<Surface>, b: Handle<Surface>) -> Option<bool> {
        let (
            Surface::Cylinder {
                def: da,
                motion: ma,
            },
            Surface::Cylinder {
                def: db,
                motion: mb,
            },
        ) = (self.surface(a), self.surface(b))
        else {
            return None;
        };
        if ma == mb {
            return Some(nacre_exact::parallel_rat(&da.dir(), &db.dir()));
        }
        let (wa, wb) = (
            self.world_cylinder_def(a).ok()?,
            self.world_cylinder_def(b).ok()?,
        );
        Some(nacre_exact::parallel_rat(&wa.dir(), &wb.dir()))
    }

    /// A plane's name normal, in the frame its truth is written in.
    fn name_normal(&self, plane: Handle<Surface>) -> Option<LineDir> {
        Some(LineDir::normal_of(self.surface_name(plane)?))
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
                    out.extend(self.line_path(self.surface_motion(plane)));
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
