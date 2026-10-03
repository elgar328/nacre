//! Doors that exist only for tests (`test-util`): unchecked pushes and `reversed_shell`.

use super::*;

impl Model {
    /// Push a plane with truth but **no name and no interning** — test-only.
    ///
    /// Two fixture populations need this door: hand-built merge fixtures that deliberately hold
    /// *one geometric plane as two handles* (interning would collapse them), and dummy planes
    /// whose handles are never dereferenced. The truth is still stated, so nothing point-less
    /// enters the arena even from tests.
    ///
    /// ⚠ The cache still follows the stated sense ([`Model::align_cache_sense`] runs at every raw
    /// push), so a model whose plane cache opposes its truth cannot be planted through here.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_plane_unregistered(
        &mut self,
        cache: nacre_geom::Plane,
        points: [[nacre_exact::Rat; 3]; 3],
        sense: Orientation,
    ) -> Handle<Surface> {
        self.push_plane_raw(
            PlanePoints::Known(points),
            None,
            sense,
            None,
            cache,
            CacheStanding::Unrealized,
        )
    }

    /// Push a cylinder whose cache is **not** derived from its statement, and not interned —
    /// test-only: the door a fixture uses to plant a cache its truth contradicts (the check that
    /// catches the lie is the proposition), which the public door would overwrite with the
    /// statement's own realization.
    #[cfg(any(test, feature = "test-util"))]
    pub fn push_cylinder_unregistered(
        &mut self,
        cache: nacre_geom::Cylinder,
        def: CylinderDef,
        motion: Option<Handle<MotionNode>>,
    ) -> Handle<Surface> {
        self.push_cylinder_raw(def, motion, cache, CacheStanding::Unrealized, false)
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
}
