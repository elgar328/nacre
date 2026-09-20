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
        coeffs: &[nacre_scalar::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, false)
    }

    /// **The single rational translation `leaf`'s whole chain amounts to**, or `None` when the
    /// chain is anything else — the third question of this family, and the one that lets a
    /// *moved* statement be restated in the world exactly.
    ///
    /// A rational translation maps a rational statement to a rational statement, so a chain made
    /// only of [`Motion::Translate`] nodes loses nothing: a plane's `d` shifts by `−n·t`
    /// ([`nacre_scalar::Isometry::plane_coeffs`]), a cylinder's origin by `+t`. What such a move
    /// *does* lose is the exactness of the `f64` **cache** — which is why the producer still
    /// records the node (`transform`'s `carry_of`, and `nacre-ops`' reuse road reads a
    /// world-stated carrier's coordinate as the statement itself). So this answers a question
    /// about *descriptions*, for the per-operation mirrors that carry them; it does not license
    /// dropping the history.
    ///
    /// ★ **The parent chain is walked raw, so a [`Motion::Frame`] node refuses outright.** A
    /// frame's expansion can come back as translations, and a statement written *under* a frame
    /// is in that frame's coordinates — folding those as world translations is a different
    /// question. Translations commute and compose by addition, so no order is implied here.
    /// `None` on overflow too (checked throughout).
    pub fn chain_translation(&self, leaf: Handle<MotionNode>) -> Option<[nacre_scalar::Rat; 3]> {
        let mut total = [nacre_scalar::Rat::from_int(0); 3];
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            let node = self.motion(h);
            let Motion::Translate { offset } = node.motion else {
                return None;
            };
            for (o, t) in total.iter_mut().zip(offset) {
                *o = o.checked_add(t)?;
            }
            cur = node.parent;
        }
        Some(total)
    }

    /// **A plane's canonical name in the world**, whatever frame its truth is written in — the
    /// door between "how this surface got here" and "where it is".
    ///
    /// Unmoved: the name itself, which is already world. Moved by a chain that folds to a
    /// rational translation ([`Model::chain_translation`]): that name carried out exactly
    /// (`d' = d − n·t`, [`nacre_scalar::Isometry::plane_coeffs`], canonicalized). `None` for
    /// anything else — a rotation, a frame node, a **moved** [`nacre_scalar::PlaneName::Wide`]
    /// (no narrow vessel for the transport to take), or an overflow — and the caller declines
    /// rather than guessing. An *unmoved* `Wide` name comes back verbatim: it is already world,
    /// and refusing it would switch off a capability that is locked elsewhere.
    ///
    /// **Planes only.** A cylinder surface has no entry in `surface_name`, so it answers `None`;
    /// its own world statement is `nacre-ops`' `world_cylinder_def`.
    pub fn world_plane_name(&self, surf: Handle<Surface>) -> Option<nacre_scalar::PlaneName> {
        let name = self.surface_name.get(&surf)?;
        match self.plane_motion(surf) {
            None => Some(name.clone()),
            Some(leaf) => {
                let t = self.chain_translation(leaf)?;
                let c = name.narrow()?;
                Some(nacre_scalar::PlaneName::Narrow(
                    nacre_scalar::Isometry::translation(t).plane_coeffs(*c)?,
                ))
            }
        }
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
        coeffs: &[nacre_scalar::Rat; 4],
    ) -> bool {
        self.chain_fixes(leaf, coeffs, true)
    }

    fn chain_fixes(
        &self,
        leaf: Handle<MotionNode>,
        coeffs: &[nacre_scalar::Rat; 4],
        rows_verbatim: bool,
    ) -> bool {
        let mut cur = Some(leaf);
        while let Some(h) = cur {
            let n: &MotionNode = self.motion(h);
            let fixed = match &n.motion {
                Motion::Rotate { axis, .. } => {
                    nacre_scalar::axis_rotation_fixes_plane(*axis, coeffs)
                }
                Motion::Translate { offset } => {
                    nacre_scalar::translation_fixes_plane(offset, coeffs)
                }
                Motion::Mirror { axis, offset } => {
                    if rows_verbatim {
                        coeffs[match axis {
                            nacre_scalar::Axis::X => 0,
                            nacre_scalar::Axis::Y => 1,
                            nacre_scalar::Axis::Z => 2,
                        }] == nacre_scalar::Rat::from_int(0)
                    } else {
                        nacre_scalar::mirror_fixes_plane(*axis, *offset, coeffs)
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
