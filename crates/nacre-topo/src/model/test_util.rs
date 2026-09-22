//! Doors that exist only for tests (`test-util`): unchecked pushes and the primitive solids.

use super::*;

/// Everything a cylinder solid is assembled from, already derived: the **exact truth** (`def`,
/// the caps' points) and its **realization** (the caches, the seam coordinates).
///
/// ★ The two cylinder constructors differ only in *how they reach here* — a caller's f64
/// statement (`Model::add_cylinder`) or an exact orthonormal frame
/// ([`Model::add_cylinder_exact`]) — so the b-rep below is written once. The f64 road cannot
/// simply delegate to the exact one: its axis comes out of `normalize()`, and a normalized f64
/// direction has no rational form to hand over.
#[cfg(any(test, feature = "test-util"))]
struct CylinderParts {
    def: CylinderDef,
    lateral: Cylinder,
    /// Bottom then top: each cap's realized plane beside its three exact points.
    caps: [(Plane, [[Rat; 3]; 3]); 2],
    /// The rims' `θ = 0` points (bottom, top) — where the seam vertices sit.
    seam_pts: [Point3; 2],
    motion: Option<Handle<MotionNode>>,
}

impl Model {
    /// Push a plane with truth but **no name and no interning** — test-only.
    ///
    /// Two fixture populations need this door: hand-built merge fixtures that deliberately hold
    /// *one geometric plane as two handles* (interning would collapse them), and dummy planes
    /// whose handles are never dereferenced. The truth is still stated, so nothing point-less
    /// enters the arena even from tests.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_plane_unregistered(
        &mut self,
        cache: nacre_geom::Plane,
        points: [[nacre_exact::Rat; 3]; 3],
    ) -> Handle<Surface> {
        self.push_plane_raw(PlanePoints::Known(points), None, None, cache)
    }

    /// A new shell whose faces are copies of `src`'s with their outward normals
    /// flipped inward: every loop's winding is reversed and every face
    /// `orientation` is toggled. Pushes fresh [`Face`] cells and a fresh
    /// [`Shell`], but **reuses** `src`'s surfaces, edges, curves, and vertices —
    /// which stay valid handles after the source solid is superseded
    /// (append-only). This is the building block for a cavity (void) shell: the
    /// boundary of a solid whose interior becomes empty space (M5
    /// containment). The reversed winding keeps each shared edge used with
    /// opposed half-edges (a valid 2-manifold), and the toggled orientation makes
    /// the outward normal point into the void.
    #[cfg(any(test, feature = "test-util"))]
    pub fn reversed_shell(&mut self, src: Handle<Shell>) -> Handle<Shell> {
        let src_faces = self.shells.get(src).faces.clone();
        let faces: Vec<Handle<Face>> = src_faces
            .iter()
            .map(|&fh| {
                let face = self.faces.get(fh).clone();
                self.faces.push(Face {
                    surface: face.surface,
                    outer: face.outer.reversed(),
                    inner: face.inner.iter().map(Loop::reversed).collect(),
                    orientation: face.orientation.flipped(),
                })
            })
            .collect();
        self.shells.push(Shell { faces })
    }

    /// Push a face **without its invariant** — test-only.
    ///
    /// `validate`'s own tests have to plant models that are wrong: a loop that does not close, a
    /// face duplicated onto another's loop, a winding deliberately reversed. Those cannot go
    /// through [`Model::push_face`], whose assert dereferences the very edge being dangled. This
    /// is the same exception [`Model::push_plane_unregistered`] is, and the only one that is ever
    /// justified: something the product cannot express, kept for the tests that must express it.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_face_unchecked(&mut self, face: Face) -> Handle<Face> {
        self.faces.push(face)
    }

    /// Push a shell **without its invariant** — test-only, see [`Model::push_face_unchecked`].
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_shell_unchecked(&mut self, shell: Shell) -> Handle<Shell> {
        self.shells.push(shell)
    }

    /// Push a solid **without making it live** — test-only.
    ///
    /// ⚠ Not [`Model::push_solid`], and the difference is the whole point: that door marks the
    /// solid live, while a planted dangling reference has to stay **unreachable**, because what it
    /// proves is that the arena-wide checks see cells nothing points at.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_solid_unlisted(&mut self, solid: Solid) -> Handle<Solid> {
        self.solids.push(solid)
    }

    /// Add an axis-aligned box `min`..`max` to this model and return its solid.
    ///
    /// Requires `max[i] > min[i]` on every axis (a degenerate box is a caller
    /// bug → panic); each
    /// face is wound so its plane normal points outward, so all faces are
    /// [`Orientation::Forward`]. Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cuboid(&mut self, min: Point3, max: Point3) -> Handle<Solid> {
        let [x0, y0, z0] = min.as_array();
        let [x1, y1, z1] = max.as_array();
        debug_assert!(
            x1 > x0 && y1 > y0 && z1 > z0,
            "cuboid needs max > min on every axis"
        );

        // 8 corners, numbered by bits (i·x, j·y, k·z).
        let corners = [
            Point3::from_array([x0, y0, z0]), // V0
            Point3::from_array([x1, y0, z0]), // V1
            Point3::from_array([x1, y1, z0]), // V2
            Point3::from_array([x0, y1, z0]), // V3
            Point3::from_array([x0, y0, z1]), // V4
            Point3::from_array([x1, y0, z1]), // V5
            Point3::from_array([x1, y1, z1]), // V6
            Point3::from_array([x0, y1, z1]), // V7
        ];
        // 6 faces: (first 3 loop vertices for the outward plane, half-edges as
        // (edge index, forward)). Loops wound CCW seen from outside → outward
        // normal. See the design plan's winding table (hand-verified).
        type FaceDef = ([usize; 3], [(usize, bool); 4]);
        let faces_def: [FaceDef; 6] = [
            ([0, 3, 2], [(3, false), (2, false), (1, false), (0, false)]), // Bottom −Z
            ([4, 5, 6], [(4, true), (5, true), (6, true), (7, true)]),     // Top +Z
            ([0, 1, 5], [(0, true), (9, true), (4, false), (8, false)]),   // Front −Y
            ([2, 3, 7], [(2, true), (11, true), (6, false), (10, false)]), // Back +Y
            ([0, 4, 7], [(8, true), (7, false), (11, false), (3, true)]),  // Left −X
            ([1, 2, 6], [(1, true), (10, true), (5, false), (9, false)]),  // Right +X
        ];
        // ★ **Surfaces before vertices**, so each corner can name the three it lies on. Separate
        // arenas, so the interleaving does not move any handle; the surfaces' order among
        // themselves is what matters and it is unchanged.
        let surf: [(Handle<Surface>, bool); 6] = core::array::from_fn(|i| {
            let (tri, _) = &faces_def[i];
            // The same three corners in rationals. `from_decimal` because a corner is a value the
            // caller *wrote* — lifting the f64 bit pattern instead would carry its drift in and
            // defeat the whole point (see `Model::surface_name`).
            let rat_corner = |k: usize| -> Option<[nacre_exact::Rat; 3]> {
                let c = corners[k].as_array();
                Some([
                    nacre_exact::Rat::from_decimal(c[0])?,
                    nacre_exact::Rat::from_decimal(c[1])?,
                    nacre_exact::Rat::from_decimal(c[2])?,
                ])
            };
            self.push_plane(
                Plane::through_points(corners[tri[0]], corners[tri[1]], corners[tri[2]])
                    .expect("non-degenerate box"),
                // ★ The same three corners, exactly. The name is derived from them, so a corner
                // whose decimals are wide enough to overflow the narrow derivation still gets
                // one. A corner outside the decimal window is a caller bug (the radius/height
                // precedent in `add_cylinder`): the truth is not optional any more.
                core::array::from_fn(|k| {
                    rat_corner(tri[k]).expect("cuboid corners inside the decimal window")
                }),
                None,
            )
        });

        let vh: [Handle<Vertex>; 8] = core::array::from_fn(|i| {
            // Which three of the six faces meet at corner `i`. The corner order is
            // `V0..V3` round the bottom then `V4..V7` round the top, so the high bit picks the cap
            // and the position round the ring picks the two walls.
            let m = i % 4;
            let cap = if i < 4 { 0 } else { 1 }; // Bottom −Z / Top +Z
            let along_x = if m == 1 || m == 2 { 5 } else { 4 }; // Right +X / Left −X
            let along_y = if m >= 2 { 3 } else { 2 }; // Back +Y / Front −Y
            self.push_vertex(
                Vertex::ThreePlane([surf[cap].0, surf[along_y].0, surf[along_x].0]),
                PointCache::Unrealized { coord: corners[i] },
            )
        });

        // 12 edges as (start, end) vertex indices plus the two faces each edge runs between
        // (indices into `faces_def` — its carriers): bottom ring, top ring, verticals.
        const EDGES: [(usize, usize, [usize; 2]); 12] = [
            (0, 1, [0, 2]), // Bottom·Front
            (1, 2, [0, 5]), // Bottom·Right
            (2, 3, [0, 3]), // Bottom·Back
            (3, 0, [0, 4]), // Bottom·Left
            (4, 5, [1, 2]), // Top·Front
            (5, 6, [1, 5]), // Top·Right
            (6, 7, [1, 3]), // Top·Back
            (7, 4, [1, 4]), // Top·Left
            (0, 4, [2, 4]), // Front·Left
            (1, 5, [2, 5]), // Front·Right
            (2, 6, [3, 5]), // Back·Right
            (3, 7, [3, 4]), // Back·Left
        ];
        // The carrier columns restate what `faces_def` already says (which loops use which
        // edge) — keep the two tables from drifting apart.
        debug_assert!(EDGES.iter().enumerate().all(|(e, &(_, _, fs))| {
            faces_def
                .iter()
                .enumerate()
                .all(|(f, (_, hes))| hes.iter().any(|&(he, _)| he == e) == fs.contains(&f))
        }));
        let eh: [Handle<Edge>; 12] = core::array::from_fn(|i| {
            let (a, b, [fa, fb]) = EDGES[i];
            self.push_edge([surf[fa].0, surf[fb].0], [vh[a], vh[b]])
                .expect("non-degenerate box")
        });

        let fh: [Handle<Face>; 6] = core::array::from_fn(|i| {
            let (_, hes) = &faces_def[i];
            let (surface, flipped) = surf[i];
            let outer = Loop {
                half_edges: hes
                    .iter()
                    .map(|&(e, forward)| HalfEdge {
                        edge: eh[e],
                        forward,
                    })
                    .collect(),
            };
            self.faces.push(Face {
                surface,
                outer,
                inner: vec![],
                // The winding table below is written for a surface whose normal this face uses as-is;
                // a shared surface may point the other way, and then the same outward direction is
                // spelled `Reversed`.
                orientation: if flipped {
                    Orientation::Forward.flipped()
                } else {
                    Orientation::Forward
                },
            })
        });

        let shell = self.shells.push(Shell { faces: fh.to_vec() });
        self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        })
    }

    /// Add a closed cylinder solid and return it: the axis runs from `base` along
    /// `axis` for `height`, with the given `radius`. The b-rep is the shared seam form
    /// (`Model::cylinder_solid`); this entry's own work is **deriving the exact truth from a
    /// caller's f64 statement**.
    ///
    /// `radius`/`height` must be positive, `axis` nonzero, and every stated coordinate inside
    /// the decimal window (caller bug → panic). Does **not** rebuild adjacency — call
    /// [`Model::rebuild_adjacency`] once after all additions.
    ///
    /// ★ **The panics are why this is a test convenience.** A statement whose *computed*
    /// coordinates leave the decimal window is not a caller bug when the caller is an
    /// application — it is an input. The production road states a cylinder through
    /// [`Model::add_cylinder_exact`], whose frame makes those coordinates rational by
    /// construction and which names every refusal instead of panicking.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cylinder(
        &mut self,
        base: Point3,
        axis: Vector3,
        radius: f64,
        height: f64,
    ) -> Handle<Solid> {
        debug_assert!(
            radius > 0.0 && height > 0.0,
            "cylinder needs positive radius and height"
        );
        let d = axis.normalize().expect("cylinder axis must be nonzero");
        let u = d
            .any_perpendicular()
            .expect("a unit axis has a perpendicular");
        let c0 = base;
        let c1 = base + d * height;
        let p_bot = c0 + u * radius; // seam point on the bottom rim (angle 0)
        let p_top = c1 + u * radius; // seam point on the top rim

        // The exact truth: a direct lift of the caller's statement — origin and the raw,
        // unnormalized axis (normalizing would destroy the exact form; the `normal_def`
        // precedent). `ref_dir` replicates `any_perpendicular`'s own rule in rationals: cross
        // the axis with the basis axis of its smallest |component| (ties X→Y→Z).
        //
        // ★ **The basis choice reads `d` — the very components `any_perpendicular` reads — so
        // agreement is structural, not order-theoretic.** Comparing the *raw*
        // components on the argument "one positive scale preserves |·| order" is a real-number
        // argument, and f64 division rounds: a strict `|x| > |y|` can collapse to equality in
        // `d`, flipping which side of the `<=` tie-break each rule lands on (measured — axis
        // `[0.34, 0.33999999999999997, 1.0]`: raw picks Y, `d` picks X, seam ~90° apart, the
        // validate net fires on a healthy model; pinned in `tests/cylinder_truth.rs`).
        //
        // The *cross* still uses the raw exact components — `ê_k × raw` is positively parallel
        // to `ê_k × d` whichever values chose `k` — so the seam direction stays exact.
        // A statement outside the decimal window is a caller bug → panic (the radius/height
        // precedent above).
        let lift = |x: f64| -> nacre_exact::Rat {
            nacre_exact::Rat::from_decimal(x).expect("cylinder statement inside the decimal window")
        };
        let def = {
            let zero = nacre_exact::Rat::from_int(0);
            // Lift then negate (not lift the negated f64): `-0.0` has no decimal of its own.
            let neg = |x: f64| {
                zero.checked_sub(lift(x))
                    .expect("negating a lifted decimal cannot overflow")
            };
            let a = axis.as_array();
            let ax = d.as_array().map(f64::abs);
            // ê_k × axis, k = the smallest-|component| basis axis of `d` — the identical
            // comparison chain `any_perpendicular` runs on the identical inputs.
            let ref_dir = if ax[0] <= ax[1] && ax[0] <= ax[2] {
                [zero, neg(a[2]), lift(a[1])]
            } else if ax[1] <= ax[2] {
                [lift(a[2]), zero, neg(a[0])]
            } else {
                [neg(a[1]), lift(a[0]), zero]
            };
            // ★ **Unreachable, and now provably so.** `CylinderDef::new` refuses a zero axis,
            // a non-positive radius, and a `ref_dir` parallel to the axis. The first two are
            // the caller's debug_asserts above; the third cannot happen here: `ref_dir` is
            // `ê_k × a` for the basis axis `k` of *smallest* |component|, and `a ∥ ê_k` would
            // need `|a_k|` to be both the largest and the smallest component — true only for
            // the zero axis. (Positive scaling preserves parallelism, so picking `k` from the
            // normalized `d` while crossing the raw `a` does not disturb the argument.)
            // Before the parallelism test became total, this `expect` also fired on statements
            // it had no business rejecting — a long-decimal axis component whose square left
            // `i128`.
            CylinderDef::new(
                base.as_array().map(lift),
                a.map(lift),
                ref_dir,
                nacre_exact::BigRat::square_of(lift(radius)),
            )
            .expect("non-degenerate cylinder")
        };

        // A cap plane's three exact points: the decimal truth of the realized center and two
        // rim-direction offsets the construction already computed (the `add_cuboid` precedent —
        // the producer's own f64 is its statement). For an axis whose normalization is exact
        // (`ẑ`, a Pythagorean triple) these are exact by construction; for an irrational axis
        // they are the truth of what was *built*, which is all a `Constructed` surface ever
        // claims. A cap outside the decimal window is a caller bug (the radius/height
        // precedent above): the truth is not optional any more.
        let w = d.cross(u);
        let cap_points = |c: Point3| -> [[nacre_exact::Rat; 3]; 3] {
            let lift = |p: Point3| -> [nacre_exact::Rat; 3] {
                p.as_array().map(|x| {
                    nacre_exact::Rat::from_decimal(x)
                        .expect("cylinder caps inside the decimal window")
                })
            };
            [lift(c), lift(c + u * radius), lift(c + w * radius)]
        };
        self.cylinder_solid(CylinderParts {
            def,
            lateral: Cylinder::from_axis(c0, d, u, radius).expect("non-degenerate cylinder"),
            caps: [
                (
                    Plane::from_point_normal(c0, -d).expect("nonzero axis"),
                    cap_points(c0),
                ),
                (
                    Plane::from_point_normal(c1, d).expect("nonzero axis"),
                    cap_points(c1),
                ),
            ],
            seam_pts: [p_bot, p_top],
            motion: None,
        })
        .expect("a rim derives its circle from the cap and the cylinder")
        .0
    }

    /// **State a cylinder exactly** — the production road, and the one with no panics in it.
    ///
    /// `axis` and `ref_dir` must be **unit and perpendicular** (checked, exactly). That single
    /// precondition is what makes everything below exact rather than realized-then-lifted: the
    /// far centre `base + axis·height`, both seam points `c + ref_dir·radius`, and each cap
    /// plane's third point `c + (axis × ref_dir)·radius` are rational products of rational
    /// inputs. `Model::add_cylinder` (the test entry) cannot do this — it normalizes an f64 axis, so its caps
    /// and seams have to be lifted back out of computed floats, and a statement whose computed
    /// coordinates leave the decimal window panics there.
    ///
    /// ★ The precondition is not a restriction on *what can be built*: an application states a
    /// cylinder on a sketch frame, and a frame either has an exact orthonormal basis or is
    /// carried by `motion` — in which case the cylinder is stated in the frame's own coordinates,
    /// where the basis is `{0, ±1}`. So the tilted case is not the irrational case.
    ///
    /// Returns the solid beside its three faces (lateral, bottom cap, top cap). Does **not**
    /// rebuild adjacency.
    #[cfg(any(test, feature = "test-util"))]
    pub fn add_cylinder_exact(
        &mut self,
        base: [Rat; 3],
        axis: [Rat; 3],
        ref_dir: [Rat; 3],
        radius: Rat,
        height: Rat,
        motion: Option<Handle<MotionNode>>,
    ) -> Result<(Handle<Solid>, [Handle<Face>; 3]), CylinderError> {
        let zero = Rat::from_int(0);
        let one = Rat::from_int(1);
        if radius <= zero {
            return Err(CylinderError::NonPositiveRadius);
        }
        if height <= zero {
            return Err(CylinderError::NonPositiveHeight);
        }
        let dot = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<Rat> {
            let mut acc = zero;
            for k in 0..3 {
                acc = acc.checked_add(a[k].checked_mul(b[k])?)?;
            }
            Some(acc)
        };
        let (aa, rr, ar) = (
            dot(&axis, &axis).ok_or(CylinderError::Overflow)?,
            dot(&ref_dir, &ref_dir).ok_or(CylinderError::Overflow)?,
            dot(&axis, &ref_dir).ok_or(CylinderError::Overflow)?,
        );
        if aa != one || rr != one || ar != zero {
            return Err(CylinderError::FrameNotOrthonormal);
        }
        // The third direction of the frame, so a cap plane gets a second rim point rather than a
        // second statement of the same one — three points on a circle name its plane.
        let cross = |a: &[Rat; 3], b: &[Rat; 3]| -> Option<[Rat; 3]> {
            let term =
                |i: usize, j: usize| a[i].checked_mul(b[j])?.checked_sub(a[j].checked_mul(b[i])?);
            Some([term(1, 2)?, term(2, 0)?, term(0, 1)?])
        };
        let w = cross(&axis, &ref_dir).ok_or(CylinderError::Overflow)?;
        // `p + dir·s`, exactly — the only arithmetic this entry does, and the reason the
        // orthonormal precondition is worth checking.
        let step = |p: &[Rat; 3], dir: &[Rat; 3], s: Rat| -> Option<[Rat; 3]> {
            let mut out = *p;
            for k in 0..3 {
                out[k] = p[k].checked_add(dir[k].checked_mul(s)?)?;
            }
            Some(out)
        };
        let over = CylinderError::Overflow;
        let c0 = base;
        let c1 = step(&c0, &axis, height).ok_or(over)?;
        let cap_points = |c: &[Rat; 3]| -> Option<[[Rat; 3]; 3]> {
            Some([*c, step(c, &ref_dir, radius)?, step(c, &w, radius)?])
        };
        let (bottom_points, top_points) =
            (cap_points(&c0).ok_or(over)?, cap_points(&c1).ok_or(over)?);
        let (p_bot, p_top) = (bottom_points[1], top_points[1]);
        // The truth carries the squared radius; a stated `radius` is squared once, here — wide,
        // so no radius is refused for the width of its square.
        let def = CylinderDef::new(base, axis, ref_dir, nacre_exact::BigRat::square_of(radius))
            .ok_or(CylinderError::Degenerate)?;

        // The caches are the realization of exactly these statements — nothing here is measured
        // or re-derived, so a cache cannot disagree with the truth beside it.
        let point = |p: [Rat; 3]| Point3::from_array(p.map(Rat::to_f64));
        let vector = |p: [Rat; 3]| Vector3::from_array(p.map(Rat::to_f64));
        let (d, u) = (vector(axis), vector(ref_dir));
        let deg = CylinderError::Degenerate;
        self.cylinder_solid(CylinderParts {
            def,
            lateral: Cylinder::from_axis(point(c0), d, u, radius.to_f64()).ok_or(deg)?,
            caps: [
                (
                    Plane::from_point_normal(point(c0), -d).ok_or(deg)?,
                    bottom_points,
                ),
                (
                    Plane::from_point_normal(point(c1), d).ok_or(deg)?,
                    top_points,
                ),
            ],
            seam_pts: [point(p_bot), point(p_top)],
            motion,
        })
        .ok_or(deg)
    }

    /// Assemble a cylinder's b-rep from parts already derived — **the one spelling** of the
    /// V2/E3/F3 seam form, shared by both constructors.
    ///
    /// Two seam vertices, two full-circle rim edges (`bounds: Some([seam, seam])`,
    /// start == end), one straight seam edge, a cylindrical lateral face whose loop uses the
    /// seam edge twice (opposite orientation), and two planar caps (each a single-half-edge rim
    /// loop). This forms a valid CW-complex (Euler χ = 2) that
    /// [`validate`](../nacre_validate/fn.validate.html) accepts — a closed periodic surface
    /// needs a seam vertex, so the rims repeat it as `[v, v]` — an endpointless edge (the
    /// standalone full circle) is a form this type does not express.
    ///
    /// Returns the solid beside its three faces in push order (lateral, bottom cap, top cap).
    /// `None` if an edge cannot derive its curve from the carriers it states — each caller says
    /// what that means for it. Does **not** rebuild adjacency.
    #[cfg(any(test, feature = "test-util"))]
    fn cylinder_solid(
        &mut self,
        parts: CylinderParts,
    ) -> Option<(Handle<Solid>, [Handle<Face>; 3])> {
        let CylinderParts {
            def,
            lateral: lateral_cache,
            caps: [(bottom_plane, bottom_points), (top_plane, top_points)],
            seam_pts: [p_bot, p_top],
            motion,
        } = parts;

        // ★ **Surfaces before edges**: an edge states its two carriers, so the lateral
        // cylinder and both cap planes must exist first. Separate arenas — the interleaving
        // moves no handle; the surfaces' order among themselves (lateral → bottom cap → top
        // cap) is what matters (the `add_cuboid` precedent).
        let lateral_surface = self.push_cylinder(lateral_cache, def, motion);
        // ★★★ **A cap's plane may already be in the model, facing the other way** — and it is the
        // *rule*, not the exception, once cylinders come from operations: the bottom cap lies on
        // the very plane the sketch frame names. `push_plane` interns by the plane's canonical
        // name, which has no direction, so it hands back "the surface exists, its cache points
        // the other way" and the face must record `Orientation::flipped()` — the `add_cuboid`
        // spelling. Dropping the bit gave a cylinder standing on a face **two upward caps**
        // (measured: bottom cap outward `+Z` on a top-face frame), which is not a solid at all.
        let (bottom_cap_surface, bottom_flipped) =
            self.push_plane(bottom_plane, bottom_points, motion);
        let (top_cap_surface, top_flipped) = self.push_plane(top_plane, top_points, motion);
        let facing = |flipped: bool| {
            if flipped {
                Orientation::Forward.flipped()
            } else {
                Orientation::Forward
            }
        };

        // A seam vertex lies on two surfaces only — the rim circle's `θ = 0` point. `OnSeam`
        // states exactly that, and the designation is complete: the cylinder's
        // truth carries `ref_dir`, so the pair means "rim ∩ the `+ref_dir` ray" — one point
        // (see `Vertex::OnSeam`; regenerating the cached coordinate is not done here).
        let v_bot = self.push_vertex(
            Vertex::OnSeam([lateral_surface, bottom_cap_surface]),
            PointCache::Unrealized { coord: p_bot },
        );
        let v_top = self.push_vertex(
            Vertex::OnSeam([lateral_surface, top_cap_surface]),
            PointCache::Unrealized { coord: p_top },
        );

        // Rims are full circles seamed at their vertex (start == end); the seam is
        // a straight edge joining the two rim seam points.
        let bottom = self.push_edge([lateral_surface, bottom_cap_surface], [v_bot, v_bot])?;
        let top = self.push_edge([lateral_surface, top_cap_surface], [v_top, v_top])?;
        // Self-adjacent: a seam is a parameterization joint of ONE surface, not an
        // intersection of two — the confirmed spelling (see `Edge::surfaces`), guarded
        // by validate's "self-adjacent ⇔ cylinder" carrier rule.
        let seam = self.push_edge([lateral_surface, lateral_surface], [v_bot, v_top])?;

        // Lateral cylindrical face: one loop wrapping the seam twice (opposite).
        let lateral = {
            let outer = Loop {
                half_edges: vec![
                    HalfEdge {
                        edge: bottom,
                        forward: true,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: true,
                    },
                    HalfEdge {
                        edge: top,
                        forward: false,
                    },
                    HalfEdge {
                        edge: seam,
                        forward: false,
                    },
                ],
            };
            self.faces.push(Face {
                surface: lateral_surface,
                outer,
                inner: vec![],
                orientation: Orientation::Forward,
            })
        };
        // Bottom cap: outward normal −d, the bottom rim reversed.
        let bottom_cap = {
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: bottom,
                    forward: false,
                }],
            };
            self.faces.push(Face {
                surface: bottom_cap_surface,
                outer,
                inner: vec![],
                orientation: facing(bottom_flipped),
            })
        };
        // Top cap: outward normal +d, the top rim forward.
        let top_cap = {
            let outer = Loop {
                half_edges: vec![HalfEdge {
                    edge: top,
                    forward: true,
                }],
            };
            self.faces.push(Face {
                surface: top_cap_surface,
                outer,
                inner: vec![],
                orientation: facing(top_flipped),
            })
        };

        let shell = self.shells.push(Shell {
            faces: vec![lateral, bottom_cap, top_cap],
        });
        let solid = self.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        Some((solid, [lateral, bottom_cap, top_cap]))
    }
}
