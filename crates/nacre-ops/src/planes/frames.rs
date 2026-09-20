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
        base_rat: Option<[nacre_scalar::Rat; 4]>,
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
        let exact = |r: nacre_scalar::Rat| nacre_scalar::Rat::try_from_f64(r.to_f64()) == Some(r);
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
        // ★★★ **Take `d` from the record and the direction from the triangle.**
        //
        // The derivation above is the two-descriptions problem in miniature: `d` comes out of an
        // f64 dot product, so the plane it names is not quite the one `tri` lies on — measured, for
        // 27% of the census's rotated classes and 40% of the fin sweep's. The surface's recorded
        // pre-motion coefficients (`Model::surface_name`) *are* that plane, exactly, with no
        // triangle in the derivation at all.
        //
        // ★ Only the **direction** still comes from the triangle, and that is deliberate. The
        // record is canonicalized, so its sign is a normal form, not this face's outward sense;
        // and the reflection correction above cannot simply be applied to a normal, because the
        // cross product is a pseudovector and reflecting-then-deriving differs from
        // deriving-then-reflecting by a global sign (the comment above, and the test
        // `a_reflected_spelling_takes_the_same_direction_signs` that found it). Orienting the
        // exact plane to agree with the derived one reproduces whatever convention the derivation
        // had, without re-deriving the convention — and *direction* is the half where the two
        // descriptions do not part.
        let exact_coeffs = base_rat.and_then(|c| {
            // ★ **The same correction, applied to the plane.** When the chain is improper the
            // triangle above was reflected in `x`, so everything derived from it lives in the
            // reflected base frame — and the record does not. Orienting the normals afterwards
            // cannot repair that: an unreflected plane and a reflected one are *different planes*,
            // not the same plane spelled with the opposite sign, so the two mirror fixtures fail
            // outright. Reflect the plane, then let the orientation step below settle the sign
            // (which is where reflecting-then-deriving and deriving-then-reflecting differ).
            let c = if improper {
                nacre_scalar::mirror_plane_coeffs(
                    c,
                    nacre_scalar::Axis::X,
                    nacre_scalar::Rat::from_int(0),
                )?
            } else {
                c
            };
            let f = c.map(|r| r.to_f64());
            // A canonicalized vector is integral; if it does not survive the round trip the
            // realization is a rounding and buys nothing over the derivation.
            c.iter()
                .zip(f)
                .all(|(&r, x)| nacre_scalar::Rat::try_from_f64(x) == Some(r))
                .then_some(f)
        });
        let coeffs = match (exact_coeffs, derived) {
            (Some(e), Some(d)) => {
                let dot = e[0] * d[0] + e[1] * d[1] + e[2] * d[2];
                Some(if dot < 0.0 { e.map(|x| -x) } else { e })
            }
            _ => derived,
        };
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
    pub(crate) base_rat: Option<[nacre_scalar::Rat; 4]>,
    /// The class root's exact rational coefficients **in the world** — see
    /// [`FaceInfo::world_rat`]. The one description the cylinder roads may compare against a
    /// world axis, and the one [`crate::combinatorics::class_coeffs_rat`] hands out.
    pub(crate) world_rat: Option<[nacre_scalar::Rat; 4]>,
    /// The class's representative surface — what `assemble_fuse_cut` records in a
    /// `Vertex::ThreePlane`.
    pub(crate) surf: Handle<Surface>,
    /// Witness points on this plane (the root face's `tri`), outward-ordered for that face.
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
    /// ★★★ **The coefficients, but only where they describe the same plane as [`tri`](Self::tri).**
    ///
    /// A plane has two exact descriptions and they need not agree: `d` is the `f64` product
    /// `raw·origin`, so a face at `y = −0.2` gets a coefficient plane `2⁻⁵⁴` from the one its own
    /// witness spans (`Plane::coefficients` has the numbers). A predicate that describes one plane
    /// by its coefficients in one question and by its triangle in the next composes answers about
    /// **two different planes**, and what comes out is not even an order.
    ///
    /// So the disagreement is resolved here rather than guarded against at every call: when the
    /// two do not agree, this is `None` and there is nothing to describe the plane with except its
    /// triangle. **The same shape [`BaseFrame`] already uses** — a description that cannot be
    /// trusted is not carried, so no consumer has to remember to check it.
    ///
    /// `None` for a rotated plane too: there are no exact `f64` coefficients for one.
    pub(crate) exact_coeffs: Option<[f64; 4]>,
    /// The **normal** under the weaker agreement — parallel to what `tri` spans, direction not
    /// required (`frame_sign` records that separately).
    ///
    /// ★ `d` is where the two descriptions part, so a predicate that never reads it can keep its
    /// exact route on a plane [`exact_coeffs`](Self::exact_coeffs) has to refuse. Measured:
    /// demanding the full agreement for those cost 4.7x on the axis-aligned fold and bought
    /// nothing.
    pub(crate) exact_normal: Option<[f64; 3]>,
    /// The class root's name integers, folded to the stored orientation
    /// ([`nacre_judge::predicate::name_stored_ints`]) — what gives a **wide** name its exact
    /// judging shortcuts back. `None` when the root's surface
    /// has no name.
    pub(crate) name_ints: Option<nacre_judge::predicate::NameInts>,
}

impl WorkingPlane {
    /// **Reconcile a plane's two exact descriptions, once, at construction.**
    ///
    /// Returns what may be carried: the coefficients when they describe the same plane the witness
    /// spans, and the normal under the weaker agreement (parallel — the direction is `frame_sign`'s
    /// to record). `None` for a rotated plane, which has no exact `f64` coefficients at all.
    ///
    /// ★ **One function, so a test fixture cannot route differently from the arrangement.** Filling
    /// the two fields by hand at a second construction site is how the fixture and the engine come
    /// to disagree about which planes are describable — and this whole item exists because two
    /// descriptions of one plane disagreed.
    pub(crate) fn reconcile(
        plane: &Plane,
        tri: [Point3; 3],
        rotated: bool,
    ) -> (Option<[f64; 4]>, Option<[f64; 3]>) {
        if rotated {
            return (None, None);
        }
        let c = plane.coefficients();
        (
            plane.spans_exactly(tri).then_some(c),
            plane.normal_spans(tri).then(|| [c[0], c[1], c[2]]),
        )
    }

    /// The class's outward normal — the root face's, which is what `tri` is wound for and what
    /// `emit_faces` winds its rings about. Not normalized: only its direction is ever read.
    pub(crate) fn tri_n_out(&self) -> Vector3 {
        (self.tri[1] - self.tri[0]).cross(self.tri[2] - self.tri[0])
    }
}
