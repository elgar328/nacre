//! [`Model`]'s methods, one file per concern. The types — and the fields these methods read —
//! are in the crate root; this file keeps construction, the prefix memo, the live-solid list and
//! the motion arena.

use super::*;

mod access;
mod cells;
mod motion;
mod planes;
mod surface_cache;
mod test_util;

impl Model {
    /// A model with the three **world axis planes pre-seeded** — surface handles 0 (XY, z = 0),
    /// 1 (YZ, x = 0), 2 (ZX, y = 0), deterministic so a replayed log and a live session name the
    /// same planes.
    ///
    /// Each seed's truth is the canonical triple `[0, u, v]` — the very points
    /// `SketchPlane::world_*` states — so any later producer of the same plane (a cuboid face on
    /// an axis, a `z = 0` sketch's base cap) **interns onto the seed**: one plane, one handle,
    /// stated once.
    ///
    /// ★ **The seed's f64 cache points down the −axis**, not up. This is the sense the dominant
    /// producers push — `extrude` builds base caps with `−normal` (its comment records that
    /// pushing `+normal` "flipped the stored normal on 781 base caps … for no reason"), and an
    /// origin cuboid's bottom/left/front faces cross to `−axis` raws — so a −axis seed keeps
    /// their `flipped` bits false and their `Orientation` spellings unchanged. The *sketch*
    /// convention ("`world_xy`'s normal is `+ẑ`") is about the frame, not the stored cache;
    /// frame derivation goes through the canonical name plus a measured `flip`, so the cache's
    /// direction never leaks to a caller.
    pub fn new() -> Self {
        let mut m = Model {
            surfaces: Store::default(),
            surface_cache: Vec::new(),
            edge_cache: Vec::new(),
            vertex_cache: Vec::new(),
            motions: Store::default(),
            motion_ids: HashMap::new(),
            surface_name: HashMap::new(),
            surface_ids: HashMap::new(),
            surface_through_ids: HashMap::new(),
            cylinder_ids: HashMap::new(),
            vertices: Store::default(),
            edges: Store::default(),
            faces: Store::default(),
            shells: Store::default(),
            solids: Store::default(),
            live_solids: Vec::new(),
            adj: Adjacency::default(),
            prefix_hp: HashMap::new(),
            world_planes: Vec::new(),
        };
        let r = nacre_exact::Rat::from_int;
        // (normal axis, +u, +v) for XY / YZ / ZX — the `SketchPlane::axis_plane` triples.
        let seeds: [([f64; 3], [i128; 3], [i128; 3]); 3] = [
            ([0.0, 0.0, 1.0], [1, 0, 0], [0, 1, 0]),
            ([1.0, 0.0, 0.0], [0, 1, 0], [0, 0, 1]),
            ([0.0, 1.0, 0.0], [0, 0, 1], [1, 0, 0]),
        ];
        let world_planes: Vec<Handle<Surface>> = seeds
            .into_iter()
            .map(|(n, u, v)| {
                let cache = nacre_geom::Plane::from_point_normal(
                    nacre_math::Point3::origin(),
                    nacre_math::Vector3::from_array(n.map(|c| -c)),
                )
                .expect("a unit axis");
                let points = [[r(0); 3], u.map(r), v.map(r)];
                let (h, flipped) = m.push_plane(cache, points, None);
                debug_assert!(!flipped, "an empty model cannot intern a seed");
                h
            })
            .collect();
        m.world_planes = world_planes;
        debug_assert_eq!(
            m.world_planes.iter().map(|h| h.index()).collect::<Vec<_>>(),
            [0, 1, 2],
            "seed handles are deterministic"
        );
        m
    }

    /// The prefix accelerator's value at one chain node — `(nodes folded, the point)` — if it is
    /// remembered. `None` is never an error: the caller folds from the base instead.
    #[inline]
    pub fn prefix_hp(
        &self,
        base: [Rat; 3],
        leaf: Handle<MotionNode>,
        prec: usize,
    ) -> Option<&PrefixValue> {
        self.prefix_hp.get(&(base, leaf, prec))
    }

    /// **Take one prefix value and leave another — "consumed" is how this table evicts.**
    ///
    /// A remembered prefix is read by exactly one successor (a transform maps one vertex to one
    /// vertex), so removing what was used and inserting what was produced keeps the table at a
    /// single live generation without tracking generations at all. `used` is the key a reader hit,
    /// `None` when it folded from the base.
    ///
    /// ⚠ **Only a caller that is extending a chain may insert**, which is why this is one door and
    /// not two: a vertex minted fresh by an arrangement is nobody's prefix, and remembering it
    /// would grow the table by one entry per boolean, forever.
    pub fn hand_over_prefix_hp(
        &mut self,
        used: Option<PrefixKey>,
        key: PrefixKey,
        value: PrefixValue,
    ) {
        if let Some(used) = used {
            self.prefix_hp.remove(&used);
        }
        self.prefix_hp.insert(key, value);
    }

    /// Drop every remembered prefix. Costs nothing but time: the next realization folds from the
    /// base and reaches the same bits (`Model::prefix_hp`'s contract).
    ///
    /// ☑ Production-facing on purpose even though nothing in this workspace calls it yet: it is
    /// the door that makes "empty is always correct" usable rather than merely true, and a
    /// consumer holding a long-lived model is the caller it is for.
    #[cfg(any(test, feature = "test-util"))]
    pub fn clear_prefix_hp(&mut self) {
        self.prefix_hp.clear();
    }

    /// How many prefixes are remembered — the instrument behind "the table stays bounded by the
    /// live generation, not by history length".
    ///
    /// ⚠ **Test-gated because it has no production consumer**, the same reason
    /// [`Model::push_plane_unregistered`] is: this counts an accelerator's internals, which is a
    /// thing to assert about, not a thing to build on. An ungated `pub` here would ship a
    /// permanent public API through the `nacre` facade for a caller that does not exist.
    #[cfg(any(test, feature = "test-util"))]
    #[inline]
    pub fn prefix_hp_len(&self) -> usize {
        self.prefix_hp.len()
    }

    /// Recompute the [`Adjacency`] cache from the current topology stores.
    /// Call once after a batch of additions (the cache is otherwise stale).
    pub fn rebuild_adjacency(&mut self) {
        // Build against an immutable borrow, then move into place — avoids
        // borrowing `self.adj` mutably while iterating the other stores.
        let adj = Adjacency::rebuild(&*self);
        self.adj = adj;
    }

    /// Push a solid into the store **and mark it live**. This is the
    /// blessed way for a producer to add a solid; the reachable closure
    /// ([`Model::reachable`]) grows to include it. Editing ops instead mutate
    /// [`Model::live_solids`] directly (drop the superseded solid, add the new).
    pub fn push_solid(&mut self, solid: Solid) -> Handle<Solid> {
        let h = self.solids.push(solid);
        self.live_solids.push(h);
        h
    }

    /// The live solids — the "current model" (supersede semantics).
    #[inline]
    pub fn live_solids(&self) -> &[Handle<Solid>] {
        &self.live_solids
    }

    /// **Drop these solids from the live set** — supersede, the editing ops' half.
    ///
    /// ⚠★★★ **Order-preserving on the survivors, and that is load-bearing, not incidental.**
    /// `nacre_step::to_step` exports the live set **in order**, so permuting it here would
    /// silently permute the exported STEP entities — and nothing would catch that: the census
    /// reads the arena (it never sees live order) and there is no golden STEP text anywhere.
    /// `Vec::retain` keeps relative order, which is why this is spelled with it.
    ///
    /// ★ Takes what to **drop**, not a keep-predicate: every caller reads "supersede these", and a
    /// predicate would make the call site say the opposite of the name.
    pub fn supersede_live(&mut self, drop: &[Handle<Solid>]) {
        self.live_solids.retain(|h| !drop.contains(h));
    }

    /// **Make a solid that is already in the arena live again.**
    ///
    /// ★ [`Model::push_solid`]'s doc has always described this move — *"editing ops instead mutate
    /// live_solids directly (drop the superseded solid, **add the new**)"* — but there was no door
    /// for the second half, so callers reached for the field. The population is the reject paths
    /// in `ops`, which retire the operands through a boolean and then have to put the original
    /// back when the op itself refuses.
    pub fn make_live(&mut self, h: Handle<Solid>) {
        debug_assert!(
            (h.index() as usize) < self.solids.len(),
            "a solid the arena does not hold cannot be live"
        );
        self.live_solids.push(h);
    }

    /// **Put the live set back** — the rollback half of a rejected operation.
    ///
    /// ★ This exists so the transaction has a name. The idiom it replaces (`let snapshot =
    /// …clone()` … `model.live_solids = snapshot`) is syntactically just an assignment and says
    /// nothing about what it is for; it is there because a rejected op once left the model changed.
    pub fn restore_live(&mut self, snapshot: Vec<Handle<Solid>>) {
        self.live_solids = snapshot;
    }

    /// Push a motion node, **interned**: the same `(motion, parent)` always yields the same
    /// handle.
    ///
    /// The handle is the canonical name of "which motion", and judgments use it to decide whether
    /// a set of points shares one rigid motion — which lets the whole judgement be answered
    /// exactly in the pre-motion frame. That identity has to be *both* collision-free and free of
    /// false misses:
    ///
    /// - a 64-bit hash of the chain's contents can collide, and a collision hands the exact
    ///   predicate two incompatible frames and answers a different question with full confidence;
    /// - a raw `Store::push` per transform is collision-free but *misses*: turning two solids by
    ///   the same 30° would make two nodes, and their shared motion would stop cancelling — so
    ///   rotating a model would turn its exact questions into assumed ones (measured: it breaks
    ///   `a_shared_rotation_still_assumes_nothing`).
    ///
    /// Interning gives both. It also keeps the forest small, since a chain shared by many solids
    /// is stored once.
    pub fn push_motion(
        &mut self,
        motion: Motion,
        parent: Option<Handle<MotionNode>>,
    ) -> Handle<MotionNode> {
        let node = MotionNode { motion, parent };
        if let Some(&h) = self.motion_ids.get(&node) {
            return h;
        }
        let h = self.motions.push(node);
        self.motion_ids.insert(node, h);
        h
    }
}
