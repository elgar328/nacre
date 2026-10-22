use super::*;
/// The identity of the motion a plane was moved by — `0` for none, else the motion-history leaf's
/// index plus one — what [`PlaneWitness::chain_id`](nacre_judge::predicate::PlaneWitness::chain_id)
/// reports.
///
/// **Identity by handle, not by hash.** A 64-bit `DefaultHasher` digest of the chain would do worse
/// than slow a judgement on a collision: it would make the shared-motion shortcuts read two
/// incompatible pre-motion frames as one, answering a different question with full confidence.
/// The motion-history leaf is the canonical name of "which motion": equal handles are the same
/// chain by construction, and two structurally-equal chains under different handles are a
/// conservative miss. The three witnesses share one chain by construction — `collect_planes`
/// builds all three from the same surface truth — so there is nothing to cross-check.
pub(super) fn chain_id_of(motion: Option<Handle<nacre_topo::MotionNode>>) -> u64 {
    motion.map_or(0, |leaf| leaf.index() as u64 + 1)
}

/// One plane of the arrangement, indexed by a **dense** class id.
///
/// The face table cannot answer "which plane" without a convention: a class holds faces from both
/// operands, and two of them can face opposite ways, so there is no such thing as *the* plane's
/// outward normal. What a plane has is a **frame** — the class root's stored normal, which is the
/// way its truth faces (`Surface::Plane::sense` against its points; the `f64` cache follows it, and
/// no direction is read off that cache) — and the only
/// direction fact anyone needs from it is [`WorkingPlane::frame_sign`]. Everything else here is a
/// witness: three points known to lie on this plane, used to reconstruct it exactly.
#[derive(Clone)]
pub(crate) struct WorkingPlane {
    pub(crate) plane: Plane,
    /// The class root's exact rational coefficients — see [`FaceInfo::base_rat`].
    pub(crate) base_rat: Option<[nacre_exact::Rat; 4]>,
    /// The class root's statement **in the world** — see [`FaceInfo::world`]. Its narrow
    /// projection ([`WorkingPlane::world_rat`]) is the one description the cylinder roads may
    /// compare against a world axis, and the one [`crate::combinatorics::class_coeffs_rat`] hands
    /// out; [`WorkingPlane::stored_world_rat`] is the same turned to face the way the plane does.
    pub(crate) world: Option<WorldName>,
    /// The class's representative surface — what `assemble_fuse_cut` records in a
    /// `Vertex::ThreePlane`.
    pub(crate) surf: Handle<Surface>,
    /// The witness as exact `WitnessPoint` definitions (the root face's), borrowed by every predicate.
    pub(crate) tri_pt3: [WitnessPoint; 3],
    /// Whether the root face's solid is rotated — copied from it together with `tri_pt3` so the
    /// pair cannot disagree. See [`FaceInfo::rotated`].
    pub(crate) rotated: bool,
    /// `+1` when the plane's stored normal agrees with the root face's outward normal, `-1` when
    /// they oppose. This *is* the label frame: `[A_above, A_below, …]` is defined about the class
    /// root's stored normal, and this sign is what relates it to material. Precomputed here so the
    /// two `debug_assert`s that guard the convention run once, at construction.
    ///
    /// **And it is the chart's frame.** The arrangement's cell walk keeps a cell **on the
    /// left of its boundary's travel in the root face's outward frame**, `n_out = frame_sign ·
    /// n_P` — so any rule that turns a 3D fact stated against the stored normal into *which
    /// half-edge's cell* multiplies by this sign exactly once. Two rules do (the ⊥ road's
    /// cut-circle disk side in `emit_faces`, the ∥ road's `ruling_interior_is_even`), and this
    /// sentence is their one home; its population is a face lying on a **seed plane** with its
    /// outward along +axis (`Model::new` plants x = 0, y = 0, z = 0 with cache direction −axis),
    /// which a max-side face reaches by a translation or a rotation.
    pub(crate) frame_sign: i8,
    /// Which motion the root face's plane carries — see [`chain_id_of`]. Planes with one nonzero
    /// id share a pre-motion frame, which their names (`name_ints`) speak.
    pub(crate) chain_id: u64,
    /// The class root's name integers, folded to the stored orientation
    /// ([`nacre_judge::predicate::name_stored_ints`]), with their nearest `f64` spellings. `None`
    /// when the root's surface has no name.
    ///
    /// **The judge's exact answers read this and nothing else** — the name-integer road (the
    /// filters over the nearest spellings, then the integers): an unmoved plane's name for the
    /// world, a moved plane's for its pre-motion frame. The name is derived from the plane's
    /// defining points without rounding, so a predicate over it answers for the truth. The plane
    /// cache's coefficients (`plane`) and the face corners' caches are rounded images: a shortcut
    /// that read them and checked them only against each other answered for the rounded model —
    /// measured, a box whose corner lies on a wall exactly (`3·0.1 = 0.3`) was judged off it and
    /// the common refused.
    pub(crate) name_ints: Option<nacre_judge::predicate::NameInts>,
}

impl WorkingPlane {
    /// The world name's narrow projection — exact rational coefficients in the world, what the
    /// cylinder roads compare against a world axis. A projection, not a second record: it cannot
    /// disagree with [`WorkingPlane::world`].
    pub(crate) fn world_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        self.world.as_ref()?.name.narrow().copied()
    }

    /// The same coefficients **turned to face the way the plane does** — the frame every cell
    /// label is written in («above» is this normal's side). Read off the truth
    /// ([`WorldName::oriented_rat`]), never off the plane cache: the name alone carries no
    /// direction, and the cache is a rounded image of the plane.
    pub(crate) fn stored_world_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        self.world.as_ref()?.oriented_rat()
    }

    /// **The class's plane facing the way it does, as an interval** — [`Self::stored_world_rat`]
    /// exactly when the class has a world statement, else the witness triangle's plane
    /// ([`nacre_judge::plane_hp`]) turned by `frame_sign`: the triangle runs outward
    /// (`n_out = frame_sign · n_P`), so `frame_sign` carries it back to the stored normal.
    ///
    /// For a class with no world statement this **is** its canonical direction too: such a class
    /// has no name whose sign could disagree, so the one direction it has — the stored one — is the
    /// one a pierce name's roots are ordered along (`n₁ × n₂`). A named class orders them along its
    /// name instead ([`Self::world_rat`]), which is why [`Self::canon_hp`] is a separate door.
    pub(crate) fn stored_hp(&self, prec: usize) -> [nacre_exact::HpBounded; 4] {
        if let Some(w) = self.stored_world_rat() {
            return w.map(|x| nacre_exact::HpBounded::of_rat(x, prec));
        }
        self.witness_hp(prec).map(|x| {
            if self.frame_sign < 0 {
                negate(&x, prec)
            } else {
                x
            }
        })
    }

    /// The class's plane facing **outward** (`n_out = frame_sign · n_P`, the frame `orient3d`
    /// answers in) as an interval — exactly from the world statement where there is one, else the
    /// witness triangle, which runs outward by construction.
    pub(crate) fn outward_hp(&self, prec: usize) -> [nacre_exact::HpBounded; 4] {
        if self.world.is_some() {
            return self.stored_hp(prec).map(|x| {
                if self.frame_sign < 0 {
                    negate(&x, prec)
                } else {
                    x
                }
            });
        }
        self.witness_hp(prec)
    }

    /// The witness triangle's plane ([`nacre_judge::plane_hp`]) — outward.
    fn witness_hp(&self, prec: usize) -> [nacre_exact::HpBounded; 4] {
        let t = &self.tri_pt3;
        nacre_judge::plane_hp(&t[0], &t[1], &t[2], prec)
    }

    /// Whether the witness triangle and the world statement agree in sign on every component the
    /// interval decides — the premise [`Self::stored_hp`] stands on for a class with no world
    /// statement, checked where both exist. `true` when there is no world statement to check.
    #[cfg(debug_assertions)]
    pub(crate) fn witness_agrees_with_world(&self, prec: usize) -> bool {
        let Some(w) = self.stored_world_rat() else {
            return true;
        };
        let t = self.witness_hp(prec);
        let zero = nacre_exact::Rat::from_int(0);
        (0..3).all(|k| {
            let stored = if self.frame_sign < 0 {
                negate(&t[k], prec)
            } else {
                t[k].clone()
            };
            w[k] == zero || stored.sign().is_none_or(|pos| pos == (w[k] > zero))
        })
    }

    /// The coefficients a pierce name's roots are ordered by — the world name exactly where there is
    /// one (what [`crate::combinatorics::class_coeffs_rat`] hands the exact solve), else the stored
    /// facing ([`Self::stored_hp`]).
    pub(crate) fn canon_hp(&self, prec: usize) -> [nacre_exact::HpBounded; 4] {
        match self.world_rat() {
            Some(w) => w.map(|x| nacre_exact::HpBounded::of_rat(x, prec)),
            None => self.stored_hp(prec),
        }
    }

    /// **The side of this class every point of `q` provably lies on**, strictly — `None` when a
    /// point is on the plane or not proved off it, when the points straddle it, and for an empty
    /// `q`. A signed distance to a plane is affine, so all of `q` on one side puts its whole convex
    /// hull there.
    ///
    /// **It never escalates**: it reads `nacre_judge::orient3d_filter`, which answers only where
    /// its error bound certifies the sign. A missed proof sends the caller back to the work it
    /// hoped to skip; it never changes an answer.
    ///
    /// **The side is spelled as [`crate::combinatorics::side_of`] spells it** — the side
    /// `Judge::orient3d` gives a point named by three planes, about this class's witness
    /// triangle. `orient3d_filter(t₀, t₁, t₂, v)` is `det[t₀ − v, t₁ − v, t₂ − v]`; the judge's
    /// sign for the same point is that of `(v − t₂)·((t₀ − t₂) × (t₁ − t₂))`, the same
    /// determinant with `v` carried from the last row to the first — a 4-cycle, odd — so the two
    /// are opposite and the filter's sign is turned here, the one place a side is minted.
    pub(crate) fn misses(&self, q: &[WitnessPoint]) -> Option<i8> {
        use nacre_exact::Orient;
        let t = &self.tri_pt3;
        let mut side: Option<Orient> = None;
        for v in q {
            match nacre_judge::orient3d_filter(&t[0], &t[1], &t[2], v) {
                Some(o) if o != Orient::Zero && side.is_none_or(|s| s == o) => side = Some(o),
                _ => return None,
            }
        }
        side.map(|o| if o == Orient::Positive { -1 } else { 1 })
    }
}

/// **Which way a face's outward runs against its plane class's root's** — `+1` when they agree.
///
/// The tiers the class merge itself stands on (`Judge::planes_coplanar`): one surface is one
/// plane, so the two outwards relate by the faces' orientations alone; one world name is one plane,
/// and the two statements' senses ([`WorldName`]) relate its two facings; otherwise the two
/// witnesses' definitions answer ([`nacre_judge::normals_agree_judge`] — each witness is wound to
/// its own face's outward by the truth). The plane cache takes no part: its normal is a rounded
/// image, and a relation between two of them is a relation between two roundings.
///
/// `Err` only where the judge could not separate the two turns, with what it did establish — a
/// reader that answers for the result records it ([`face_facing_noted`]); one that can do without the
/// answer (reuse, which then arranges the class) drops it.
pub(crate) fn face_facing(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    fp: usize,
    c: usize,
) -> Result<i8, nacre_judge::Decision> {
    let (face, root) = (faces[fp].plane(), &jd.planes[c]);
    if face.surf == root.surf {
        return Ok(face.orient_sign * root.frame_sign);
    }
    if let (Some(fw), Some(rw)) = (&face.world, &root.world) {
        if fw.name == rw.name {
            return Ok(face.orient_sign * fw.sense.sign() * rw.sense.sign() * root.frame_sign);
        }
    }
    let (a, b) = (&face.tri_pt3, &root.tri_pt3);
    match nacre_judge::normals_agree_judge(
        (&a[0], &a[1], &a[2]),
        (&b[0], &b[1], &b[2]),
        jd.standard,
    ) {
        nacre_judge::Decision::Sign(nacre_exact::Orient::Positive) => Ok(1),
        nacre_judge::Decision::Sign(nacre_exact::Orient::Negative) => Ok(-1),
        outcome => Err(outcome),
    }
}

/// [`face_facing`] for a reader whose answer the result rests on: an undecided facing is recorded
/// as evidence, so the operation refuses by that judgement's name (`undecided_reject`) before
/// anything built on it is used — the rule every undecided judgement follows.
pub(crate) fn face_facing_noted(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    fp: usize,
    c: usize,
) -> Option<i8> {
    face_facing(jd, faces, fp, c)
        .map_err(|outcome| {
            jd.notes.push(nacre_judge::predicate::Evidence {
                site: nacre_judge::predicate::Site::FacesAgree { face: fp, class: c },
                outcome,
            })
        })
        .ok()
}

/// `−x` as an interval.
fn negate(x: &nacre_exact::HpBounded, prec: usize) -> nacre_exact::HpBounded {
    nacre_exact::HpBounded::of_rat(nacre_exact::Rat::from_int(0), prec).sub(x, prec)
}
