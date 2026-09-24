use super::*;
/// One plane of the arrangement, indexed by a **dense** class id.
///
/// The face table cannot answer "which plane" without a convention: a class holds faces from both
/// operands, and two of them can face opposite ways, so there is no such thing as *the* plane's
/// outward normal. What a plane has is a **frame** — the class root's stored normal — and the only
/// direction fact anyone needs from it is [`WorkingPlane::frame_sign`]. Everything else here is a
/// witness: three points known to lie on this plane, used to reconstruct it exactly.
/// The pre-rotation twin of a witness triangle: its `chain_id`, base points and base plane.
///
/// A rigid motion preserves the determinants the predicates take, so a judgement whose inputs all
/// carry **one** motion can be answered on these instead — exactly, off the toleranced path
/// entirely. `None` for the base data when a base coordinate is not `f64`-representable, since the
/// exact predicate takes `f64`; the judgement then stays toleranced (slower, never wrong).
#[derive(Clone, Copy, Debug)]
pub(crate) struct BaseFrame {
    /// `0` = no motion. Equal only for structurally identical chains.
    pub(crate) chain_id: u64,
    pub(crate) tri: Option<[Point3; 3]>,
    pub(crate) coeffs: Option<[f64; 4]>,
}

impl BaseFrame {
    /// No motion to cancel — for a plane with no history: a hand-built table in a test, or the
    /// synthetic split plane a subdivided boolean cuts with. Identical to what `of` returns for an
    /// unmoved face, so such a plane takes the same predicate routes an axis-aligned model does.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Self {
            chain_id: 0,
            tri: None,
            coeffs: None,
        }
    }

    pub(super) fn of(
        tri_pt3: &[WitnessPoint; 3],
        motion: Option<Handle<nacre_topo::MotionNode>>,
        frame_sign: i8,
        name: Option<&nacre_exact::PlaneName>,
    ) -> Self {
        // **Identity by handle, not by hash.** This used to fold the chain into a 64-bit
        // `DefaultHasher` digest and compare digests — and a collision does not make a judgement
        // slow, it makes `shared_base` hand the *exact* predicate two incompatible pre-motion
        // frames and answer a different question with full confidence. The motion-history leaf is
        // the canonical name of "which motion": equal handles are the same chain by construction,
        // and two structurally-equal chains under different handles are a conservative miss.
        //
        // The three witnesses share one chain by construction — `collect_planes` builds all three
        // from the same surface truth — so there is nothing to cross-check here either.
        let Some(leaf) = motion else {
            return Self {
                chain_id: 0,
                tri: None,
                coeffs: None,
            }; // no motion to cancel
        };
        // 0 is reserved for "no motion", so never hand it out as an id.
        let chain_id = leaf.index() as u64 + 1;
        let exact = |r: nacre_exact::Rat| nacre_exact::Rat::try_from_f64(r.to_f64()) == Some(r);
        if !tri_pt3.iter().all(|p| p.base.iter().all(|&r| exact(r))) {
            return Self {
                chain_id,
                tri: None,
                coeffs: None,
            };
        }
        let pt = |p: &WitnessPoint| {
            Point3::from_array([p.base[0].to_f64(), p.base[1].to_f64(), p.base[2].to_f64()])
        };
        let tri = [pt(&tri_pt3[0]), pt(&tri_pt3[1]), pt(&tri_pt3[2])];
        // ★ **An improper chain is corrected here, on the points, before anything is derived
        // from them.** An odd number of reflections leaves the base frame with the opposite
        // handedness, and the shortcuts' whole licence is that the motion preserves the
        // determinants they take. One more reflection puts the handedness back, and then the base
        // is related to the moved frame by a *proper* motion again — which is what every consumer
        // of this struct assumes. A sign flip on x is exact for every finite `f64`, and it is the
        // same convention `nacre_judge::frame3::shared_base` applies to its own points.
        let improper = nacre_judge::chain_parity(&tri_pt3[0].chain) < 0;
        let tri = if improper {
            tri.map(|p| {
                let [x, y, z] = p.as_array();
                Point3::from_array([-x, y, z])
            })
        } else {
            tri
        };
        // ★ The base plane must carry the **stored** orientation, not the triangle's. A class's
        // stored normal and its witness triangle's `cross` can oppose — that is exactly what
        // `frame_sign` records — and `through_points` gives the triangle's. A proper motion
        // preserves the cross product (`det = 1`), so multiplying by `frame_sign` reproduces the
        // same relation in the base frame. Without it the exact path answers with a flipped sign,
        // which the suite caught immediately.
        //
        // ★ **And it is derived from the corrected triangle, not corrected afterwards.** A plane
        // is not a bag of points: reflecting a triangle and re-deriving its normal is *not* the
        // same as reflecting the normal, because the cross product is a pseudovector — the two
        // differ by a global sign, and `plane_pair_dir_sign` reads exactly that sign. Deriving
        // last removes the question: the points are corrected once, and everything downstream is
        // the ordinary derivation from them. (The earlier spelling corrected the plane separately
        // and was off by that one sign; `a_reflected_spelling_takes_the_same_direction_signs`
        // is what found it.)
        let derived = Plane::through_points(tri[0], tri[1], tri[2]).map(|pl| {
            let c = pl.coefficients();
            let k = f64::from(frame_sign);
            [c[0] * k, c[1] * k, c[2] * k, c[3] * k]
        });
        // ★★★ **Take the plane from the record, not from the triangle.**
        //
        // The derivation above is the two-descriptions problem in miniature: `d` comes out of an
        // f64 dot product, so the plane it names is not quite the one `tri` lies on — measured, for
        // 27% of the census's rotated classes and 40% of the fin sweep's. The surface's recorded
        // pre-motion coefficients (`Model::surface_name`) *are* that plane, exactly, with no
        // triangle in the derivation at all.
        //
        // ★ **And its direction from σ**, which reads the same record against the witness exactly
        // ([`nacre_judge::predicate::witness_name_sense`] — the fold `name_ints` makes). The record
        // is canonical, so its sign is a normal form, not this face's orientation; σ is what
        // relates the two.
        let exact_coeffs = name.and_then(|name| {
            let c = *name.narrow()?;
            let sigma = nacre_judge::predicate::witness_name_sense(name, tri_pt3, frame_sign)?;
            // ★ **The same correction, applied to the plane.** When the chain is improper the
            // triangle above was reflected in `x`, so everything derived from it lives in the
            // reflected base frame — and the record does not. An unreflected plane and a
            // reflected one are *different planes*, not one plane spelled with the opposite sign,
            // so the record is reflected too: `x ↦ −x` negates its `x` coefficient (`d` stays —
            // the mirror is at `0`). Left uncanonical, because the sign is set next: the
            // reflection reverses the triangle's turn (the cross product is a pseudovector) and
            // not the record's, so the reflected record faces the reflected triangle's stored
            // direction exactly when σ says it does not.
            let zero = nacre_exact::Rat::from_int(0);
            let (v, s) = if improper {
                ([zero.checked_sub(c[0])?, c[1], c[2], c[3]], -sigma)
            } else {
                (c, sigma)
            };
            let v = if s < 0 {
                [
                    zero.checked_sub(v[0])?,
                    zero.checked_sub(v[1])?,
                    zero.checked_sub(v[2])?,
                    zero.checked_sub(v[3])?,
                ]
            } else {
                v
            };
            let f = v.map(|r| r.to_f64());
            // A canonicalized vector is integral; if it does not survive the round trip the
            // realization is a rounding and buys nothing over the derivation.
            v.iter()
                .zip(f)
                .all(|(&r, x)| nacre_exact::Rat::try_from_f64(x) == Some(r))
                .then_some(f)
        });
        let coeffs = exact_coeffs.or(derived);
        Self {
            chain_id,
            tri: Some(tri),
            coeffs,
        }
    }
}

#[derive(Clone)]
pub(crate) struct WorkingPlane {
    pub(crate) plane: Plane,
    /// The class root's exact rational coefficients — see [`FaceInfo::base_rat`].
    pub(crate) base_rat: Option<[nacre_exact::Rat; 4]>,
    /// The class root's canonical name **in the world** — see [`FaceInfo::world_name`]. Its narrow
    /// projection ([`WorkingPlane::world_rat`]) is the one description the cylinder roads may
    /// compare against a world axis, and the one [`crate::combinatorics::class_coeffs_rat`] hands
    /// out.
    pub(crate) world_name: Option<nacre_exact::PlaneName>,
    /// The class's representative surface — what `assemble_fuse_cut` records in a
    /// `Vertex::ThreePlane`.
    pub(crate) surf: Handle<Surface>,
    /// The root face's corner caches (`f64`), wound outward for that face — read by
    /// [`WorkingPlane::tri_n_out`] alone, never by a judgement.
    pub(crate) tri: [Point3; 3],
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
    /// ★★★ **And it is the chart's frame.** The arrangement's cell walk keeps a cell **on the
    /// left of its boundary's travel in the root face's outward frame**, `n_out = frame_sign ·
    /// n_P` — so any rule that turns a 3D fact stated against the stored normal into *which
    /// half-edge's cell* multiplies by this sign exactly once. Two rules do (the ⊥ road's
    /// cut-circle disk side in `emit_faces`, the ∥ road's `ruling_interior_is_even`), and this
    /// sentence is their one home; its population is a face lying on a **seed plane** with its
    /// outward along +axis (`Model::new` plants x = 0, y = 0, z = 0 with cache direction −axis),
    /// which a max-side face reaches by a translation or a rotation.
    pub(crate) frame_sign: i8,
    /// The pre-rotation twin — see [`BaseFrame`].
    pub(crate) base: BaseFrame,
    /// The class root's name integers, folded to the stored orientation
    /// ([`nacre_judge::predicate::name_stored_ints`]), and their `f64` row where it is exact.
    /// `None` when the root's surface has no name.
    ///
    /// ★★★ **An unmoved plane's exact shortcuts read this and nothing else** (the `PlaneWitness`
    /// impl in `tolerant`). The name is derived from the plane's defining points without
    /// rounding, so a predicate over its row answers for the truth. The plane cache's
    /// coefficients (`plane`) and the face corners' caches (`tri`) are rounded images: a shortcut
    /// that read them and checked them only against each other answered for the rounded model —
    /// measured, a box whose corner lies on a wall exactly (`3·0.1 = 0.3`) was judged off it and
    /// the common refused. A wide name has no `f64` row and keeps its BigInt rescue.
    pub(crate) name_ints: Option<nacre_judge::predicate::NameInts>,
}

impl WorkingPlane {
    /// The world name's narrow projection — exact rational coefficients in the world, what the
    /// cylinder roads compare against a world axis. A projection, not a second record: it cannot
    /// disagree with [`WorkingPlane::world_name`].
    pub(crate) fn world_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        self.world_name.as_ref()?.narrow().copied()
    }

    /// The class's outward normal — the root face's, which is what `tri` is wound for and what
    /// `emit_faces` winds its rings about. Not normalized: only its direction is ever read.
    pub(crate) fn tri_n_out(&self) -> Vector3 {
        (self.tri[1] - self.tri[0]).cross(self.tri[2] - self.tri[0])
    }
}
