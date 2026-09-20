use super::*;

use crate::SketchFrame;

use nacre_scalar::Axis;

/// The engine entry with the evidence dropped — these tests assert geometry, and the report
/// has its own tests. Shadows [`super::boolean`] so the call sites read as they always did.
fn boolean(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Vec<Handle<Solid>>, BoolError> {
    super::boolean(model, kind, a, b).map(|(solids, ..)| solids)
}

// --- Rotation-generality: the trace engine's decisions are coordinate-free. Rigidly rotating
// both operands by the same isometry must leave the result invariant. A non-90° angle flags the
// solid, whose faces record the motion, routing every predicate to the exact frame3 backend. `rot30` (lib.rs)
// is in a sibling test module and unreachable here, so the isometries are built inline.

/// 30° about `axis` through (1,1,0) — non-90°, so the motion is recorded (exact frame3 path).
fn rot_iso(axis: nacre_scalar::Axis) -> nacre_scalar::Isometry {
    use nacre_scalar::{Angle, Isometry, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// Rotate a solid by each axis in turn. A compound tilt needs `rebuild_adjacency` BETWEEN the
/// transforms (matching the oracle's `rotated_boolean_matches_occt`), or the second reads a
/// stale topology.
fn tilt(m: &mut Model, mut s: Handle<Solid>, axes: &[nacre_scalar::Axis]) -> Handle<Solid> {
    for &ax in axes {
        s = transform(m, s, &rot_iso(ax)).unwrap();
        m.rebuild_adjacency();
    }
    s
}

/// The z=1 class both solids seat a cap on.
fn shared_cap_class(
    m: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
) -> usize {
    let n_planes = plane_ix
        .iter()
        .map(|c| c.plane())
        .max()
        .map_or(0, |m| m + 1);
    (0..n_planes)
        .find(|&c| {
            let seats = |s: Handle<Solid>| {
                solid_shell_handles(m, s).into_iter().any(|sh| {
                    m.shell(sh).faces.iter().any(|fh| {
                        plane_ix[surf_ix[fh]].plane() == c && face_on_z1(*fh, surf_ix, faces)
                    })
                })
            };
            seats(a) && seats(b)
        })
        .expect("a shared z=1 cap class")
}

fn face_on_z1(fh: Handle<Face>, surf_ix: &HashMap<Handle<Face>, usize>, faces: &[FaceRow]) -> bool {
    let p = &faces[surf_ix[&fh]];
    // A cap in the z=1 plane: all three defining points at z=1.
    p.plane()
        .tri
        .iter()
        .all(|q| (q.as_array()[2] - 1.0).abs() < 1e-12)
}

#[derive(Default)]
struct Counts {
    classes: usize,
    arranged: usize,
    pairs: usize,
    cullable: usize,
    one_piece: usize,
    needs_seed: usize,
    /// Pairs in classes whose footprint is one piece — the realistic cull, since the
    /// multi-piece ones cannot take a single seed.
    pairs_1p: usize,
    cullable_1p: usize,
    /// Of `needs_seed`, the ones where the solid's whole box misses the footprint, so it
    /// **cannot** enclose it and the seed is false without any work.
    seed_free: usize,
}

impl Counts {
    fn add(&mut self, o: &Counts) {
        self.classes += o.classes;
        self.arranged += o.arranged;
        self.pairs += o.pairs;
        self.cullable += o.cullable;
        self.one_piece += o.one_piece;
        self.needs_seed += o.needs_seed;
        self.pairs_1p += o.pairs_1p;
        self.cullable_1p += o.cullable_1p;
        self.seed_free += o.seed_free;
    }
}

mod curved;
mod declines_and_spikes;
mod planar_cells;
mod rotated;
