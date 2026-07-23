//! Operations for the nacre kernel, plus a replayable operation log (design §6).
//!
//! [`Operation::Extrude`] (M2) sweeps a planar polygon profile into a prism;
//! [`Operation::PadOnFace`]/[`Operation::PocketOnFace`] (M4) consume a prior op's face by
//! `Handle` (exposed via [`OpOutput`]) and supersede a solid (design §2 live-solid
//! semantics) — each is a tool prism plus a boolean, not a direct face-split. Ops are
//! applied by [`apply`] and folded by [`replay`]; every result is a **closed** solid, so
//! `nacre-validate` applies fully.

use nacre_math::{Point2, Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Face, HalfEdge, Model, Solid, Vertex};

mod arrangement;
mod boolean;
mod combinatorics;
mod ops;
mod planes;
mod tolerant;
mod transform;

pub use boolean::boolean;
pub use ops::{BoolKind, OpError, OpOutput, Operation, Profile2d, SketchPlane, apply, replay};

impl SketchPlane {
    /// The world XY plane (normal +Z).
    pub fn world_xy() -> Self {
        Self {
            origin: Point3::origin(),
            x_axis: Vector3::from_array([1.0, 0.0, 0.0]),
            y_axis: Vector3::from_array([0.0, 1.0, 0.0]),
        }
    }

    /// A plane through `origin` with the given `normal`; `x`/`y` axes are
    /// synthesized (`x = n.any_perpendicular()`, `y = n × x`). `None` if `normal`
    /// is zero.
    pub fn from_origin_normal(origin: Point3, normal: Vector3) -> Option<Self> {
        let n = normal.normalize()?;
        let x = n.any_perpendicular()?;
        Some(Self {
            origin,
            x_axis: x,
            y_axis: n.cross(x),
        })
    }

    /// The 3-D point for sketch coordinates `p = (u, v)`.
    #[inline]
    pub fn point(&self, p: Point2) -> Point3 {
        self.origin + self.x_axis * p[0] + self.y_axis * p[1]
    }

    /// The plane normal `x × y` (unit when the axes are unit and orthogonal).
    #[inline]
    pub fn normal(&self) -> Vector3 {
        self.x_axis.cross(self.y_axis)
    }
}

/// Why a boolean could not be computed. The engine rejects out-of-coverage
/// input honestly rather than returning a plausibly-wrong solid (overview
/// 불리언 전략).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolError {
    /// Outside the current coverage: a non-`Common` kind, a non-planar face, a
    /// non-convex input, coplanar faces (across or within an input), a 4-plane
    /// concurrency, a tangential contact, or any other general-position
    /// violation.
    Unsupported,
    /// An input solid handle is not in `model.live_solids`.
    InputNotLive,
}

/// Names of the `Unsupported` reject sites, shared by the guard that raises one
/// and the test that asserts it — a renamed tag then cannot silently drift out
/// of a test's expectation. Not `#[cfg(test)]`: the guards name these in release
/// builds too. They are `const`, so they inline away where the tag is unused.
pub(crate) mod tag {
    /// An outer-shell edge used by other than two face loops. A backstop with no
    /// firing test: `validate` calls this `NonOpposedEdge` and every shell the
    /// operations build is manifold — but `boolean` never runs `validate` on its
    /// inputs, so a direct caller could still hand one in.
    pub const NON_MANIFOLD_EDGE: &str = "non_manifold_edge";
    /// A closed seam loop's edges disagree about which side its material lies on, or
    /// two of its nodes coincide, or two neighbours share no plane pair.
    ///
    /// A cross-check rather than a defence: the exact
    /// predicate (`order_along`) and the orientation bookkeeping (`FaceInfo::n_out`
    /// vs its plane's normal) must agree, edge by edge, and neither is assumed right.
    /// It survives release on purpose. Downstream, `validate` catches a wrong loop and
    /// `tessellate` catches a wrong *hole* — but a caller may run neither.
    ///
    /// It cannot catch a *globally* flipped loop: every edge would be wrong together.
    /// A globally flipped loop would need a golden on the loop producer itself.
    ///
    /// Unreachable today: "material on the left" is a global property, so a consistent
    /// loop makes every edge agree. Not dead code — relax `FOURPLANE` or cell 3c's
    /// node-identity argument and this is what speaks first.
    pub const LOOP_ORIENT_MISMATCH: &str = "loop_orient_mismatch";
    pub const COPLANAR_PAIR: &str = "coplanar_pair";
    /// The result severs into two or more material solids *and* at least one enclosed void
    /// (cavity) survives. Which outer shell owns which cavity needs a shell-scoped point-in-shell
    /// test we do not have yet, so this is honestly rejected and deferred to a follow-on cell.
    /// Reachable: `Cut` a hollow part with a cut that isolates the void into one severed piece.
    /// Born with its firing test (`severed_with_cavity_is_rejected`).
    ///
    /// ★ The starting point for that cell is the removed `point_in_solid_idx` (excised by the
    /// seam-engine cell, 2026-07-22 dev-log) once the arrangement stopped needing it (it seeds the unbounded cell and propagates inward
    /// instead). It answers *point in solid* by ray casting; what this needs is *point in shell*,
    /// so recover it from that commit and re-scope it rather than deriving one from scratch.
    pub const SEVERED_WITH_CAVITY: &str = "severed_with_cavity";
    /// No material-enclosing (outward) shell among the result components — every component is
    /// inward-oriented. Geometrically impossible for a real solid result; a defensive backstop
    /// with no firing test (cf. `FOURPLANE`).
    pub const NO_OUTWARD_SHELL: &str = "no_outward_shell";
    /// A boolean input is a rotated solid (overhaul stage 1b). Rotated planar geometry
    /// is representable but its predicates are not yet sound (no CIP until stage 3), so
    /// the boolean honestly rejects until then. Fired by a `Transform`-rotated operand.
    pub const ROTATED_UNSUPPORTED: &str = "rotated_unsupported";
    pub const THREE_PLANES: &str = "three_planes";
    /// Two **different** arrangement vertices (distinct plane triples) materialized to the same
    /// coordinate. The triple is the truth and the coordinate only its cache (overview §5), so this
    /// says the exact substrate and the f64 cache disagree about how many vertices exist — always a
    /// defect upstream, never a property of the input. Raised where the seam table is built, while
    /// both triples are still in hand; without it the disagreement surfaces much later as a
    /// zero-length edge. The known cause is a **split plane table** (one geometric plane carried by
    /// two classes); a genuine 4-plane concurrency would do the same.
    pub const SEAM_ALIAS: &str = "seam_alias";
    /// A result loop asked for an edge between two vertices at the same coordinate. Every ring node
    /// is a distinct arrangement vertex, so this cannot happen for well-named input — it is the
    /// backstop that keeps a degenerate one from aborting the kernel (`Line::through_points` used to
    /// `expect`). `SEAM_ALIAS` catches the known cause earlier, so this has no firing test.
    pub const ZERO_LENGTH_EDGE: &str = "zero_length_edge";
    pub const FOURPLANE: &str = "fourplane";
    pub const CYLINDER_FACE: &str = "cylinder_face";
    pub const DEGENERATE_FACE: &str = "degenerate_face";
    pub const DEGENERATE_NORMAL: &str = "degenerate_normal";
    /// Every candidate ray from a loop's nodes has a ring node on its line.
    ///
    /// `point_in_ring` casts along `P ∩ Q_a` for a node's own plane `Q_a`; a ring node on
    /// that line makes the crossing parity ambiguous. Candidates are `2 · |loop|` lines and
    /// two directions, and half of them can be spoiled at once — `l_and_staple`'s loop and
    /// arc share both `y` planes, so only the `x` lines are clear there. Unfired today.
    pub const NO_CLEAR_RAY: &str = "no_clear_ray";
    /// A loop's node lies *on* the ring it is being tested against.
    ///
    /// A hole ring never touches the outer ring it sits in, and `point_in_ring` checks that
    /// exactly: the ray's line meets an edge at `X`, and `X == v` strictly inside that edge means
    /// `v` is on the ring. Unfired.
    pub const POINT_ON_RING: &str = "point_on_ring";
    /// The trace arrangement on one plane class nested a hole whose containment depth exceeds one,
    /// or produced more than one unbounded contour (several disjoint bodies on the plane).
    /// `nest_cells` resolves any number of holes at depth one inside one outer loop; deeper nesting
    /// and multiple bodies are honestly rejected until the general nesting cell lands. Two distinct
    /// tags so a refactor cannot silently merge the conditions.
    pub const HOLE_DEPTH: &str = "hole_depth";
    pub const HOLE_ROOTS: &str = "hole_roots";
    /// A face whose boundary never crosses the seam, yet the seam lies on its plane — the
    /// convex path only.
    ///
    pub const MISSING_SEAM: &str = "missing_seam";
    /// A `Whole`-survival contact face whose footprint OVERLAPS the other's (∂P × ∂Q cross) rather
    /// than nesting, in the one such case still unbuilt. `Whole` has two entries: `Fuse`/same-normal,
    /// which the E1 union cell now builds, and `Cut`/opposite-normal, which is exact whenever the
    /// contact plane separates the two solids (nothing to remove). What is left is a `Cut` whose tool
    /// reaches back across that plane — a pin below its own contact face — where the cut owes a notch
    /// this path cannot yet cut. Honest reject rather than a whole cap that ignores the pin.
    pub const COPLANAR_MERGE: &str = "coplanar_merge";
}

#[cfg(test)]
thread_local! {
    static LAST_REJECT: std::cell::Cell<Option<&'static str>> =
        const { std::cell::Cell::new(None) };
}

/// Build an `Unsupported`, recording *which* guard raised it. `Err(Unsupported)`
/// alone cannot distinguish the guards, so a reject test whose fixture drifts
/// onto a different guard would still pass — [`assert_rejects`] closes that hole.
/// Every `Unsupported` site in this crate goes through here: two call sites
/// discard a `collect_planes` error (`detect_coincident_interface`), so only
/// exhaustive tagging makes "last tag written == the site that returned" hold.
#[inline]
#[cfg_attr(not(test), allow(unused_variables))]
pub(crate) fn reject(tag: &'static str) -> BoolError {
    #[cfg(test)]
    LAST_REJECT.with(|c| c.set(Some(tag)));
    BoolError::Unsupported
}

/// Assert that `f` rejects *through the intended guard*. Clears any stale tag
/// first, so a prior call in the same test cannot be mistaken for this one.
#[cfg(test)]
fn assert_rejects<T: std::fmt::Debug + PartialEq>(
    f: impl FnOnce() -> Result<T, BoolError>,
    expect: &'static str,
) {
    LAST_REJECT.with(|c| c.take());
    assert_eq!(f(), Err(BoolError::Unsupported));
    assert_eq!(LAST_REJECT.with(|c| c.take()), Some(expect));
}

/// The start vertex of a half-edge (`bounds[0]` if forward, else `bounds[1]`).
/// Every half-edge walked here belongs to a valid solid, so its edge is bounded.
pub(crate) fn he_start(model: &Model, he: HalfEdge) -> Handle<Vertex> {
    let [a, b] = model
        .edges
        .get(he.edge)
        .bounds
        .expect("a solid's loop edge is bounded");
    if he.forward { a } else { b }
}

use std::collections::HashMap;

fn unordered(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

#[cfg(test)]
pub mod tests {

    use super::*;
    use crate::transform::transform;
    use crate::{boolean::*, ops::*, planes::*};
    use nacre_geom::intersect::{planes_coplanar, three_planes};
    use nacre_geom::{Plane, Surface};
    use nacre_scalar::frame3::Pt3;
    use nacre_topo::{Loop, Orientation, Origin, Shell, VertexDef};
    use proptest::prelude::*;
    use std::collections::HashMap;

    /// Test shim: a boolean whose result is exactly one solid. Most tests operate on a single
    /// body; this asserts that and returns the lone handle, so call sites read as before while
    /// `boolean` itself returns the full `Vec` (cell 0.4 multi-solid).
    fn boolean_one(
        model: &mut Model,
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> Result<Handle<Solid>, BoolError> {
        let solids = boolean(model, kind, a, b)?;
        assert_eq!(
            solids.len(),
            1,
            "boolean_one: expected one solid, got {}",
            solids.len()
        );
        Ok(solids[0])
    }

    fn p2(x: f64, y: f64) -> Point2 {
        Point2::from_array([x, y])
    }

    /// Is there an outer-shell face on the plane through `pt` with normal `n`, oriented that way?
    /// The production path names a cap by the *face* that made it (`find_face_coplanar_with`); a
    /// test that wants to say "a face sits on z = 1.5 facing +z" has no such face in hand, and
    /// asserting geometry from coordinates is exactly what a test may do.
    fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
        let Some(target) = Plane::from_point_normal(pt, n) else {
            return false;
        };
        let shell = m.solids.get(solid).outer;
        m.shells.get(shell).faces.iter().any(|&fh| {
            let f = m.faces.get(fh);
            let Surface::Plane(plane) = m.surfaces.get(f.surface) else {
                return false;
            };
            let sign = match f.orientation {
                Orientation::Forward => 1.0,
                Orientation::Reversed => -1.0,
            };
            planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n) > 0.0
        })
    }

    fn square() -> Profile2d {
        Profile2d {
            points: vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)],
        }
    }

    fn extrude_op(profile: Profile2d, dist: f64) -> Operation {
        Operation::Extrude {
            plane: SketchPlane::world_xy(),
            profile,
            dist,
        }
    }

    fn regular_ngon(n: usize, r: f64) -> Profile2d {
        let points = (0..n)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                p2(r * a.cos(), r * a.sin())
            })
            .collect();
        Profile2d { points }
    }

    #[test]
    fn square_extrudes_to_a_cube() {
        let m = replay(&[extrude_op(square(), 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.edges.len(), 12);
        assert_eq!(m.faces.len(), 6);
        assert_eq!(m.solids.len(), 1);

        let mut got: Vec<[f64; 3]> = m.vertices.iter().map(|(_, v)| v.point.as_array()).collect();
        let mut want = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let key = |p: &[f64; 3]| (p[0].to_bits(), p[1].to_bits(), p[2].to_bits());
        got.sort_by_key(key);
        want.sort_by_key(key);
        assert_eq!(got, want);
    }

    #[test]
    fn triangle_extrudes_to_a_prism() {
        let tri = Profile2d {
            points: vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(1.0, 1.5)],
        };
        let m = replay(&[extrude_op(tri, 3.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 6);
        assert_eq!(m.edges.len(), 9);
        assert_eq!(m.faces.len(), 5);
    }

    #[test]
    fn pentagon_extrudes_clean() {
        let m = replay(&[extrude_op(regular_ngon(5, 2.0), 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 10);
        assert_eq!(m.faces.len(), 7);
    }

    #[test]
    fn concave_l_profile_is_valid() {
        // An L-shape (a reflex vertex) — a simple concave hexagon.
        let l = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.vertices.len(), 12);
        assert_eq!(m.faces.len(), 8);
    }

    /// The L-prism: profile `[(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)]` extruded to
    /// z ∈ [0,1]. Material = bottom bar (x∈[0,2],y∈[0,1]) ∪ left bar (x∈[0,1],
    /// y∈[1,2]); the notch (x∈[1,2],y∈[1,2]) is empty.
    fn l_prism() -> (Model, Handle<Solid>) {
        let l = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// The same L-prism, its profile started one vertex earlier so the reflex corner
    /// `(1,1)` lands at index 1 of the cap's loop. Geometrically identical.
    fn rotated_l_prism() -> (Model, Handle<Solid>) {
        let l = Profile2d {
            points: vec![
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
                p2(0.0, 0.0),
                p2(2.0, 0.0),
            ],
        };
        let m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// `FaceInfo::n_out` is documented as the single source of "outward". Two
    /// independent sources say which way that is: the ring, which winds CCW about the
    /// outward normal, and the b-rep's own `Surface` plus `Orientation`. They must agree
    /// on every face of every solid.
    ///
    /// `outer_tri` used to read the turn at the first non-collinear corner, which is the
    /// ring's winding **only when that corner is convex**. Nothing enforced that. The
    /// four fixtures below were safe by accident — none starts its cap loop one vertex
    /// before a reflex corner. `rotated_l_prism` does, and it is the same solid.
    #[test]
    fn outward_normals_agree_with_their_orientation() {
        // `collect_planes` debug_asserts, per face, that `sign(normal·n_out)` equals the topo
        // `face.orientation` — the invariant this test used to check by reading a `.orient` field
        // it kept alongside. That field is gone (its production role is `orient_sign`), so the
        // check lives at construction now; running collect_planes on shapes with reversed faces
        // (the L/U notch, a rotated solid) exercises it. Here we add the parallel invariant, which
        // is not debug_asserted the same way.
        let mut cube = Model::new();
        let c = cube.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let (ml, sl) = l_prism();
        let (mu, su) = u_prism();
        let (mr, sr) = rotated_l_prism();
        for (name, m, s) in [
            ("cube", &cube, c),
            ("l_prism", &ml, sl),
            ("u_prism", &mu, su),
            ("rotated_l_prism", &mr, sr),
        ] {
            for pi in &collect_planes(m, s).unwrap() {
                let dot = pi.plane.normal().dot(pi.n_out);
                assert!(
                    dot.abs() > 0.5,
                    "{name}: n_out is not parallel to its plane"
                );
            }
        }
    }


    /// `pt3_base_collinear` is exact on the pre-rotation rational bases: three genuinely
    /// collinear points stay collinear under rotation (→ skipped), and a real sliver (one point
    /// off the line) is never falsely called collinear (→ its crossing is kept, no silent-wrong).
    #[test]
    fn pt3_base_collinear_exact() {
        use nacre_scalar::{Angle, Axis, Rat};
        let ang = Angle::from_deg(Rat::from_int(37)).unwrap();
        let piv = [Rat::from_int(2), Rat::from_int(-1), Rat::from_int(0)];
        let rp = |x: i128, y: i128, z: i128| {
            Pt3::at([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)]).rotate_about(
                Axis::Z,
                ang,
                piv,
            )
        };
        // (0,0,0), (2,4,6), (1,2,3): all on the line t·(1,2,3) → collinear.
        assert!(pt3_base_collinear(&rp(0, 0, 0), &rp(2, 4, 6), &rp(1, 2, 3)));
        // (1,2,4) is off that line (z), a real nonzero-area triangle → not collinear.
        assert!(!pt3_base_collinear(
            &rp(0, 0, 0),
            &rp(2, 4, 6),
            &rp(1, 2, 4)
        ));
    }

    /// The L-prism with a `[0.1,0.9]³` box strictly inside its bottom bar
    /// (non-coplanar coordinates ⇒ no shared face planes). `V_L = 3`, `V_box =
    /// 0.512`.
    fn l_and_inner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
        (m, l, bx)
    }






    /// The L-prism with a box biting its convex corner `(2, 0)` — the first
    /// non-convex *overlap* (a real single-chord seam), M5-d2. The box spans
    /// `x∈[1.3,2.4]`, `y∈[-0.3,0.4]`, `z∈[0.2,1.4]`: it straddles the corner in x
    /// and y, and its z-range pokes above the L (`z=1`) while its floor `z=0.2`
    /// sits inside — so every crossing edge is a clean straddle (no edge tunnels
    /// fully through the other) and no box face is coplanar with an L face. The
    /// span is deliberately asymmetric so no seam point lands on a face centre
    /// (where both fan diagonals cross and every apex would graze).
    /// Overlap = `x∈[1.3,2]·y∈[0,0.4]·z∈[0.2,1]` = `0.224`;
    /// `V_L=3`, `V_box=1.1·0.7·1.2=0.924`.
    fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([1.3, -0.3, 0.2]),
            Point3::from_array([2.4, 0.4, 1.4]),
        );
        (m, l, bx)
    }





    /// The L with a box straddling its reflex corner (1,1): a *single* chord with one
    /// reflex bend. The box vertex `(1.6,1.6,·)` sits in the L's notch — inside the
    /// convex hull, outside the L — exactly where a convex half-space test would
    /// misclassify it `Inside`. Overlap = `xy(1.0 − notch 0.36) · z(0.8)` = `0.512`.
    fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.6, 0.6, 0.2]),
            Point3::from_array([1.6, 1.6, 1.4]),
        );
        (m, l, bx)
    }


    /// The same two solids the other way round: the bar severs the rod, and `Cut` answers with
    /// two solids (cell 0.4). Each severed stub is its own genus-0 box, so `validate` is clean —
    /// the pieces share no vertices or edges. (Before cell 0.4 this was `DISCONNECTED_RESULT`:
    /// one handle could not name two solids, and forcing both into one shell read as
    /// `NegativeGenus { genus: -1 }`. `pierced_multi` had been hiding it: severing A takes an
    /// edge of A through B.)
    #[test]
    fn cut_rod_by_l_severs_into_two() {
        let (mut m, l, rod) = l_and_rod();
        let solids = boolean(&mut m, BoolKind::Cut, rod, l).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((vol - 0.06).abs() < 1e-9, "total volume {vol}");
    }

    // --- Rotated booleans go live (overhaul 3d-i, `ROTATED_UNSUPPORTED` retired) ---
    // A boolean commutes with a rigid motion, so rotating both operands by the same
    // irrational-angle isometry must give the rigid image of the unrotated result — identical
    // volume, solid count, and cavity count, and still valid. These are the first live proof
    // that the CIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
    // sound end-to-end on rotated (rounded-irrational) geometry.






    /// Adversarial rotation stress (overhaul 3d-iii): many fixtures × kinds × rotations
    /// (single-axis, and Euler chains reaching arbitrary orientation) confirm the DNA
    /// invariant — a rotated boolean is *never silently wrong*: its result either equals the
    /// unrotated one (a boolean commutes with a rigid motion, so volume/solid-count/cavity-count
    /// are invariant) or is an honest reject. `#[ignore]`: each rotated boolean escalates its
    /// CIP predicates to astro-float and costs ~0.5–2.5 s, so this runs on demand, not per commit
    /// (the invariance regression guard is the fast `rotated_*` tests above).
    #[test]
    #[ignore = "slow: rotated booleans ~2s each (run with --ignored)"]
    fn rotation_invariance_stress() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        #[derive(PartialEq, Debug)]
        enum Out {
            Rej,
            Ok(f64, usize, usize),
        }
        let rot = |axis: Axis, deg: i128, piv: [i128; 3]| {
            Isometry::rotation(Rotation {
                axis,
                point: [
                    Rat::from_int(piv[0]),
                    Rat::from_int(piv[1]),
                    Rat::from_int(piv[2]),
                ],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            })
        };
        let run = |build: &dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>),
                   kind: BoolKind,
                   isos: &[Isometry]|
         -> Out {
            let (mut m, mut a, mut b) = build();
            for iso in isos {
                a = transform(&mut m, a, iso).unwrap();
                m.rebuild_adjacency();
                b = transform(&mut m, b, iso).unwrap();
                m.rebuild_adjacency();
            }
            match boolean(&mut m, kind, a, b) {
                Ok(solids) => {
                    m.rebuild_adjacency();
                    assert!(
                        nacre_validate::validate(&m).is_empty(),
                        "INVALID rotated result"
                    );
                    let vol: f64 = solids
                        .iter()
                        .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
                        .sum();
                    let cav: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
                    Out::Ok(vol, solids.len(), cav)
                }
                Err(_) => Out::Rej,
            }
        };
        let vclose =
            |x: f64, y: f64| (x - y).abs() <= 1e-6 || (x - y).abs() <= 1e-4 * x.abs().max(y.abs());
        let matches = |base: &Out, r: &Out| match (base, r) {
            (Out::Rej, Out::Rej) => true,
            (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2)) => vclose(*v1, *v2) && s1 == s2 && c1 == c2,
            _ => false,
        };
        let cube = |lo: [f64; 3], hi: [f64; 3]| {
            move || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
                let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
                (m, a, b)
            }
        };
        type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;
        let fixtures: Vec<(&str, Build)> = vec![
            ("corner", Box::new(l_and_corner_box)),
            (
                "rod",
                Box::new(|| {
                    let (m, l, r) = l_and_rod();
                    (m, r, l) // sever = Cut(rod, l): a=rod, b=l
                }),
            ),
            ("inner", Box::new(l_and_inner_box)),
            ("cube_corner", Box::new(cube([0.5; 3], [1.5; 3]))),
        ];
        let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
            ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
            ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
            (
                "Z50>X37",
                vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::X, 37, [0, 0, 1])],
            ),
            // A three-axis Euler chain reaches an arbitrary orientation (axes are X/Y/Z only).
            (
                "Z30>X30>Y73",
                vec![
                    rot(Axis::Z, 30, [1, 1, 0]),
                    rot(Axis::X, 30, [0, 0, 0]),
                    rot(Axis::Y, 73, [0, 2, 0]),
                ],
            ),
        ];
        let (mut success, mut reject, mut silent, mut skipped) = (0, 0, 0, 0);
        for (fname, build) in &fixtures {
            for kind in [BoolKind::Cut, BoolKind::Fuse, BoolKind::Common] {
                let base = run(build.as_ref(), kind, &[]);
                for (rname, isos) in &isos_list {
                    if matches!(base, Out::Rej) {
                        skipped += 1;
                        continue;
                    }
                    let r = run(build.as_ref(), kind, isos);
                    if matches!(r, Out::Rej) {
                        reject += 1; // an honest reject is acceptable, not a counterexample
                    } else if matches(&base, &r) {
                        success += 1;
                    } else {
                        silent += 1;
                        eprintln!("SILENT-WRONG {fname} {kind:?} {rname} base={base:?} rot={r:?}");
                    }
                }
            }
        }
        eprintln!(
            "ROTATION STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent} skipped(base_rej)={skipped}"
        );
        assert_eq!(
            silent, 0,
            "a rotated boolean was silently wrong (valid but != unrotated)"
        );
    }

    /// A signature that changes if `Store::push` order (hence handle identity) changes:
    /// vertex points in handle order — exactly what `assemble_fuse_cut` assigns by first
    /// appearance across `faces` — plus edge/face/solid counts and sorted volumes.
    #[cfg(feature = "parallel")]
    fn model_sig(m: &Model, solids: &[Handle<Solid>]) -> String {
        use std::fmt::Write;
        let mut s = String::new();
        let _ = write!(s, "S{}", solids.len());
        for (_, v) in m.vertices.iter() {
            let p = v.point.as_array();
            let _ = write!(
                s,
                "|{:x},{:x},{:x}",
                p[0].to_bits(),
                p[1].to_bits(),
                p[2].to_bits()
            );
        }
        let _ = write!(s, "|E{}F{}", m.edges.len(), m.faces.len());
        let mut vols: Vec<u64> = solids
            .iter()
            .map(|&sh| nacre_props::mass_props(m, sh).unwrap().volume.to_bits())
            .collect();
        vols.sort_unstable();
        let _ = write!(s, "|V{vols:?}");
        s
    }

    /// The parallel boolean must be bit-identical regardless of rayon thread count — replay
    /// determinism (DNA) requires thread-order independence. A rotated multi-face Fuse
    /// exercises the parallel reconstruction/classification; its result under a 1-thread
    /// pool (par_iter code, index order) must equal the default many-thread result, run
    /// repeatedly so scheduling jitter would show.
    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_boolean_is_thread_order_independent() {
        use nacre_scalar::{Axis, Isometry, Rat};
        let build = || {
            let ngon = |n: usize, r: f64, cx: f64, cy: f64| Profile2d {
                points: (0..n)
                    .map(|i| {
                        let ang = std::f64::consts::TAU * (i as f64) / (n as f64);
                        p2(cx + r * ang.cos(), cy + r * ang.sin())
                    })
                    .collect(),
            };
            let mut m = replay(&[
                extrude_op(ngon(16, 2.0, 0.0, 0.0), 3.0),
                extrude_op(ngon(16, 2.0, 2.5, 0.5), 3.0),
            ])
            .unwrap();
            let a = m.live_solids[0];
            let b = m.live_solids[1];
            let up = Isometry::translation([Rat::from_int(0), Rat::from_int(0), Rat::from_int(1)]);
            let b = transform(&mut m, b, &up).unwrap();
            m.rebuild_adjacency();
            let tilt = rot_iso(Axis::X, 30);
            let a = transform(&mut m, a, &tilt).unwrap();
            m.rebuild_adjacency();
            let b = transform(&mut m, b, &tilt).unwrap();
            m.rebuild_adjacency();
            (m, a, b)
        };
        let run = || {
            let (mut m, a, b) = build();
            let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
            m.rebuild_adjacency();
            model_sig(&m, &solids)
        };
        let pool1 = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let reference = pool1.install(run);
        for _ in 0..8 {
            assert_eq!(
                run(),
                reference,
                "parallel boolean result depends on thread order"
            );
        }
    }

    /// A U-prism: a bottom bar `y∈[0,1]` with two prongs rising from it. The prong
    /// tops sit at *different* heights (y=2.3 and y=2.0) on purpose — level tops
    /// would be coplanar faces and `has_coplanar_pair` would reject before the seam
    /// machinery ran. Area 3 + 1 + 1.3, extruded 1.0 ⇒ volume 5.3.
    fn u_prism() -> (Model, Handle<Solid>) {
        let u = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(3.0, 0.0),
                p2(3.0, 2.3),
                p2(2.0, 2.3),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let m = replay(&[extrude_op(u, 1.0)]).unwrap();
        let s = m.live_solids[0];
        (m, s)
    }

    /// The U with a slab shearing off both prong tops. The slab overhangs the U in
    /// x and z, so **every slab edge lies outside the U** (pierces nothing) and every
    /// U edge either straddles cleanly (one crossing of the slab's `y=1.5` face) or
    /// misses. No edge threads the other solid, so the arcs are the whole story.
    fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, u) = u_prism();
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.5, -0.5]),
            Point3::from_array([3.5, 2.5, 1.5]),
        );
        (m, u, slab)
    }

    #[test]
    fn u_prism_is_valid() {
        // Pin the fixture itself: a mistyped profile could still trip `multichord`
        // below, for the wrong reason.
        let (m, u) = u_prism();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, u).unwrap().volume;
        assert!((vol - 5.3).abs() < 1e-9, "volume {vol}");
    }

    // ---- arrangement: seam segment gathering (M5-d3 cell 3b) ----

    fn near(a: Point3, b: [f64; 3]) -> bool {
        (a - Point3::from_array(b)).norm() < 1e-9
    }

    proptest! {}

    /// A big cube whose `y=0, z=0` edge is crossed **twice** by the seam: the notch
    /// spans `x∈[3,7]` and hangs below both `y=0` and `z=0`, so that one edge enters
    /// and leaves it. Extents are asymmetric so no crossing lands on a face centre or a
    /// fan diagonal — a debt to the fan, which cell (5a) deleted. Kept: one variable at
    /// a time.
    ///
    /// `edge_seam` maps an edge to *one* seam triple. This input is what that map
    /// cannot represent — and `boolean` rejects it today (see
    /// `an_edge_crossed_twice_is_rejected`).
    fn cube_and_notch() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let y = m.add_cuboid(
            Point3::from_array([3.0, -1.0, -1.0]),
            Point3::from_array([7.0, 1.4, 1.2]),
        );
        (m, a, y)
    }

    /// A thin rod skewering the L's bottom bar in `z`, both ends outside. Each of its four
    /// vertical edges pierces the L's two caps, so the caps take a closed seam loop and the
    /// rod's walls take two chords apiece — and each wall's vertical edges are crossed
    /// **twice**, leaving runs with no vertex at all.
    fn l_and_rod() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let rod = m.add_cuboid(
            Point3::from_array([0.3, 0.3, -0.5]),
            Point3::from_array([0.5, 0.6, 1.5]),
        );
        (m, l, rod)
    }

    /// A slab over the pocketed cube, its underside at height `z0`. The rectangle is
    /// asymmetric so that the cube's four vertical edges, which pierce the underside at
    /// `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`, miss its fan diagonals; a square slab has all
    /// four apexes degenerate at once.
    ///
    /// `z0 = 0.3` runs below the pocket floor, so the slab's underside meets only the
    /// cube's outer walls. `z0 = 0.7` runs between the floor and the lid and meets the
    /// pocket walls as well — two loops, nested.
    fn pocket_and_slab(z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, pc) = pocketed_cube();
        let slab = m.add_cuboid(
            Point3::from_array([-0.2, -0.25, z0]),
            Point3::from_array([1.3, 1.2, 1.5]),
        );
        (m, slab, pc)
    }

    /// Nested loops, at last (cell 3f-7). Cell 3f-3 argued nesting was reachable and reached
    /// for a polyhedral torus; a pocket sliced between its floor and its lid is enough. The
    /// slab's underside carries the cube's cross-section as one loop and the pocket's inside
    /// it — a loop within a loop, which the pairwise guard over-rejected.
    ///
    /// `Cut` opens, both orders: `A ∩ B = B ∩ {z ≥ 0.7}` = `0.30 − 0.048 = 0.252`, slab `1.74`,
    /// pocketed cube `0.92`, so `Cut(slab,pc) = 1.488` and `Cut(pc,slab) = 0.668`. `Cut(pc,slab)`
    /// is the one that makes the cube cross-section an *island with a hole* (the slab is B, its
    /// underside dropped, so no region survives to own the loops); `Cut(slab,pc)` keeps the slab
    /// region and hangs the pocket loop in it as a hole beside an island.
    ///
    /// `Fuse` seals the pocket (`[0.3,0.7]² × [0.5,0.7]`, capped by the slab at `z = 0.7`) into
    /// an enclosed cavity — a second shell. Cell (5c) assembles it: the union material is
    /// `2.408` and the result carries one cavity of volume `0.032`, two shells, `validate`
    /// clean (the void's inward orientation). Both orders (Fuse is commutative).
    #[test]
    fn a_slab_between_the_lid_and_the_floor_nests_two_loops() {
        for (swap, expect) in [(false, 1.488), (true, 0.668)] {
            let (mut m, slab, pc) = pocket_and_slab(0.7);
            let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
            let r = boolean_one(&mut m, BoolKind::Cut, x, y).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "Cut swap={swap}: {vs:?}");
            let props = nacre_props::mass_props(&m, r).unwrap();
            assert!(
                (props.volume - expect).abs() < 1e-9,
                "Cut swap={swap}: {} vs {expect}",
                props.volume
            );
        }
        // Fuse seals the pocket into a cavity — a second shell, assembled by cell (5c).
        for swap in [false, true] {
            let (mut m, slab, pc) = pocket_and_slab(0.7);
            let (x, y) = if swap { (pc, slab) } else { (slab, pc) };
            let r = boolean_one(&mut m, BoolKind::Fuse, x, y).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "Fuse swap={swap}: {vs:?}");
            let props = nacre_props::mass_props(&m, r).unwrap();
            assert!(
                (props.volume - 2.408).abs() < 1e-9,
                "Fuse swap={swap}: {}",
                props.volume
            );
            assert_eq!(m.solids.get(r).cavities.len(), 1, "Fuse swap={swap}");
            assert_eq!(m.reachable().shells.len(), 2, "Fuse swap={swap}");
        }
    }

    /// `is_shell_outward` — the exact sign `assemble_fuse_cut` labels components by
    /// ((5d)#5, replacing the f64 signed-volume flux) — is true for an outward,
    /// material-enclosing shell and false for an inward void shell. `reversed_shell`
    /// flips one into the other, so the same faces read opposite orientations.
    #[test]
    fn is_shell_outward_true_for_outer_false_for_void() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let outer = m.solids.get(cube).outer;
        let out_faces = m.shells.get(outer).faces.clone();
        assert!(is_shell_outward(&m, &out_faces), "outer shell is outward");
        let void = m.reversed_shell(outer);
        let void_faces = m.shells.get(void).faces.clone();
        assert!(
            !is_shell_outward(&m, &void_faces),
            "reversed shell is a void"
        );
    }




    /// The L with a box biting its reflex corner and poking out the top. The box top
    /// (z=1.2) clears the L's z=1 **deliberately**: sunk inside the L's slab, the L's
    /// vertical edges at (2,1) and (1,1) would pierce the box's bottom *and* top face,
    /// which `pierced_multi` used to reject. Cell 3e-3 supports it; the fixture keeps its
    /// clearance so that it goes on testing one thing.
    fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.2]),
            Point3::from_array([2.5, 1.5, 1.2]),
        );
        (m, l, bx)
    }

    #[test]
    fn cut_staircase_seam_arc() {
        // On the box's bottom face the seam runs (2,0.5) → (2,1) → (1,1) → (1,1.5): a
        // staircase whose two bends turn opposite ways. `strict` used to reject it,
        // unable to tell a reflex turn from an arc folded back on itself.
        //
        // The reconstructed face there is (0.5,0.5) → (0.5,1.5) → (1,1.5) → (1,1) →
        // (2,1) → (2,0.5): the overlap footprint, area 1.0, reflex at (1,1). A correct
        // simple polygon. `strict` was rejecting a right answer.
        //
        // Overlap = footprint 1.0 × z∈[0.2,1] = 0.8.
        let (mut m, l, bx) = l_and_popup_box();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.8)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn fuse_staircase_seam_arc() {
        // The `Fuse` counterpart, closing the inclusion–exclusion: V_L + V_box − 0.8.
        // `validate` cannot see a self-intersecting face (it stays manifold, Euler
        // holds), so the volume is what pins the folded arc — with OCCT alongside.
        let (mut m, l, bx) = l_and_popup_box();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 2.0 - 0.8)).abs() < 1e-9, "volume {vol}");
    }

    /// The L-prism and an L-shaped bar lying in its notch, biting two convex corners of
    /// the L's top face. The bar spans `z ∈ [0.5, 1.5]`, so its body clears the cap.
    ///
    /// Each bite crosses **two different** edges of the cap, which is exactly why no edge
    /// is pierced twice — the bar takes corners, not edges. Two chords, no closed loop.
    fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let bar = Profile2d {
            points: vec![
                p2(1.8, 0.8),
                p2(2.1, 0.8),
                p2(2.1, 2.1),
                p2(0.8, 2.1),
                p2(0.8, 1.8),
                p2(1.8, 1.8),
            ],
        };
        let OpOutput::Extrude { solid: b, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 0.0, 0.5]),
                    ..SketchPlane::world_xy()
                },
                profile: bar,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, b)
    }



    /// The L-prism with an **L-shaped** stub standing wholly inside its top face,
    /// `z ∈ [0.5, 1.5]`. The seam on the cap is a closed loop with a reflex node — the
    /// suite's first non-convex inner loop, and the shape a winding must be read from.
    ///
    /// Its coordinates dodge the cap's fan diagonals from `(0,0)` (`y = x`, `y = x/2`,
    /// `y = 2x`), which `segment_crosses_face` would graze.
    fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let ell = Profile2d {
            points: vec![
                p2(0.2, 0.25),
                p2(0.85, 0.25),
                p2(0.85, 0.4),
                p2(0.35, 0.4), // reflex
                p2(0.35, 0.9),
                p2(0.2, 0.9),
            ],
        };
        let OpOutput::Extrude { solid: stub, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 0.0, 0.5]),
                    ..SketchPlane::world_xy()
                },
                profile: ell,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, stub)
    }

    /// A ring's nodes as coordinates — each triple is three planes, so its point is their meet.
    fn ring_points(planes: &[PlaneGeom], ring: &[[usize; 3]]) -> Vec<[f64; 3]> {
        ring.iter()
            .map(|t| {
                three_planes(
                    &planes[t[0]].plane,
                    &planes[t[1]].plane,
                    &planes[t[2]].plane,
                )
                .unwrap()
                .as_array()
            })
            .collect()
    }

    #[test]
    fn a_hole_winds_clockwise_and_an_island_counter_clockwise() {
        // The two rings of a holed face are stored with opposite windings — that is what makes one
        // a hole and the other its outer boundary — and `loop_winding` must read exactly that.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        assert_eq!(combinatorics::loop_winding(&planes, p, &outer).unwrap(), 1);
        assert_eq!(combinatorics::loop_winding(&planes, p, &hole).unwrap(), -1);

        // Nothing but the ring's direction went into that. Reversing it by hand agrees.
        for (name, ring, want) in [("outer", &outer, -1i8), ("hole", &hole, 1)] {
            let mut reversed = ring.clone();
            reversed.reverse();
            assert_eq!(
                combinatorics::loop_winding(&planes, p, &reversed).unwrap(),
                want,
                "{name} reversed"
            );
        }
    }

    #[test]
    fn a_reflex_node_turns_against_its_ring() {
        // The L cap's outer ring is a hexagon with exactly one reflex corner, at `(1, 1)`. A
        // convex node turns with the ring and the reflex one turns against it, so the turn signs
        // are *not* all equal — which is why the winding cannot be read off an arbitrary node.
        let (planes, p, outer, _) = holed_face_rings("dimple");
        let pts = ring_points(&planes, &outer);
        let reflex = pts
            .iter()
            .position(|q| near(Point3::from_array(*q), [1.0, 1.0, 1.0]))
            .expect("the reflex node");
        assert_eq!(
            combinatorics::turn_at(&planes, p, &outer, reflex).unwrap(),
            -1
        );
        let turns: Vec<i8> = (0..outer.len())
            .map(|i| combinatorics::turn_at(&planes, p, &outer, i).unwrap())
            .collect();
        assert_eq!(
            turns.iter().filter(|&&t| t == -1).count(),
            1,
            "one reflex corner: {turns:?}"
        );

        // The node `loop_winding` lands on is the lexicographically least — a hull vertex, where
        // the turn *is* the winding. The test finds it by reading coordinates; `loop_winding`
        // finds it with an exact predicate.
        let lo = (0..pts.len())
            .min_by(|&i, &j| pts[i].partial_cmp(&pts[j]).unwrap())
            .unwrap();
        assert_ne!(lo, reflex);
        assert_eq!(combinatorics::turn_at(&planes, p, &outer, lo).unwrap(), 1);

        // ★ The teeth. A ring is a cycle, so its winding cannot depend on where the walk began.
        // Start it at the reflex node and a `turn_at(ring[0])` implementation reads the reflex
        // sign — the exact fault `outer_tri` shipped. Measured: without this rotation, such an
        // implementation passes every assertion above.
        let mut rotated = outer.clone();
        rotated.rotate_left(reflex);
        assert_eq!(combinatorics::turn_at(&planes, p, &rotated, 0).unwrap(), -1);
        assert_eq!(
            combinatorics::loop_winding(&planes, p, &rotated).unwrap(),
            1
        );
    }

    #[test]
    fn cut_ell_dimple() {
        // The blind pocket is the stub's L-shaped section, `0.65·0.15 + 0.15·0.5 = 0.1725`,
        // half a unit deep. Area `14 − 0.1725 + 2.6·0.5 + 0.1725`: the lid gives up exactly
        // what the floor hands back, so only the walls move it.
        let (mut m, l, stub) = l_and_ell_stub();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (3.0 - 0.1725 * 0.5)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 15.3).abs() < 1e-9, "area {}", props.area);
    }

    /// The L-prism and a П-shaped staple straddling the L's reflex corner. The profile
    /// lives in the **XZ** sketch plane and extrudes along `−y`, so the L's cap (`z = 1`)
    /// is *parallel* to the extrusion axis and the staple's section there falls into two
    /// pieces: one wholly inside the cap, one wrapping the corner `(1,1)`.
    ///
    /// That parallelism is the whole point. A prism cut by a plane **perpendicular** to
    /// its axis meets a face in the profile, which is connected — so every component of
    /// `profile ∩ f` reaches `∂f`, and a face can never carry both an arc and a loop. Every
    /// earlier attempt at such a fixture died on that.
    ///
    /// Leg bottoms sit at `z = 0.5` and `z = 0.45`: two coplanar faces of *one* operand
    /// trip `coplanar_pair` at the door, exactly as `u_prism`'s staggered prongs avoid.
    /// And the legs span `y ∈ [0.65, 1.3]`, not `[0.7, 1.3]`, because `(1.4, 0.7)` lies on
    /// the cap's fan diagonal `y = x/2` and `segment_crosses_face` would graze it.
    fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let staple = Profile2d {
            points: vec![
                p2(0.1, 0.5),
                p2(0.6, 0.5),
                p2(0.6, 1.3),
                p2(0.8, 1.3),
                p2(0.8, 0.45),
                p2(1.4, 0.45),
                p2(1.4, 1.5),
                p2(0.1, 1.5),
            ],
        };
        let OpOutput::Extrude { solid: st, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane {
                    origin: Point3::from_array([0.0, 1.3, 0.0]),
                    x_axis: Vector3::from_array([1.0, 0.0, 0.0]),
                    y_axis: Vector3::from_array([0.0, 0.0, 1.0]),
                },
                profile: staple,
                dist: 0.65,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output")
        };
        (m, l, st)
    }

    /// One face ring, as plane triples.
    type Ring = Vec<[usize; 3]>;

    /// A **holed reflex face** from the live engine, as plane triples: `(planes, p, outer, hole)`.
    ///
    /// `Cut(L-prism, stub)` leaves the L's top cap carrying a hole — `"dimple"` a square one,
    /// `"ell"` an L-shaped one. Both rings come from [`combinatorics::face_vertex_triples`] and
    /// [`combinatorics::hole_rings`], which the boolean itself uses, so the fixture exercises only code
    /// the kernel runs.
    ///
    /// This replaces a helper that built its rings from `seam_paths_on`/`orient_seam_loop` — the
    /// retired seam engine. The *properties* below are about `point_in_ring`/`every_ray`, which are
    /// live and load-bearing (`nest_cells` picks a hole's host with them, `unify_coplanar_faces`
    /// groups by them), so they had to be re-homed rather than deleted with their old fixture.
    fn holed_face_rings(which: &str) -> (Vec<PlaneGeom>, usize, Ring, Ring) {
        let (mut m, l, stub) = if which == "dimple" {
            l_and_dimple()
        } else {
            l_and_ell_stub()
        };
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).expect("the cut");
        m.rebuild_adjacency();
        let faces_tab = collect_planes(&m, r).unwrap();
        let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
        for (i, pi) in faces_tab.iter().enumerate() {
            surf_ix.insert(pi.face, i);
        }
        let canon = plane_classes(&faces_tab);
        let (planes, plane_ix) = dense_planes(&faces_tab, &canon);
        let inc = combinatorics::edge_faces(&m, r, &surf_ix).unwrap();
        for &fh in &m.shells.get(m.solids.get(r).outer).faces {
            let fp = surf_ix[&fh];
            let holes = combinatorics::hole_rings(&m, fh, fp, &inc, &planes, &plane_ix).unwrap();
            if let Some(hole) = holes.into_iter().next() {
                let outer =
                    combinatorics::face_vertex_triples(&m, fh, fp, &inc, &planes, &plane_ix)
                        .unwrap();
                assert_eq!(outer.len(), 6, "{which}: the L's cap is a reflex hexagon");
                return (planes, plane_ix[fp], outer, hole);
            }
        }
        panic!("{which}: no holed face");
    }

    #[test]
    fn a_loop_is_inside_the_face_it_was_found_on() {
        // A hole ring never touches its face's outer ring: every node is *strictly* inside, so
        // `point_in_ring` must say so for all of them. The outer ring is the L's cap, a hexagon
        // with a reflex corner, so this is not a convex test.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        for t in &hole {
            assert!(combinatorics::point_in_ring(&planes, p, *t, &outer).unwrap());
        }
    }

    #[test]
    fn the_older_loops_are_inside_their_faces_too() {
        // Two hole shapes that must not move: a square one and an L-shaped one. Both sit
        // strictly inside the same reflex hexagon, and **every** clear ray agrees — the parity
        // cannot depend on which ray was cast, which is a second machine for free.
        for which in ["dimple", "ell"] {
            let (planes, p, outer, hole) = holed_face_rings(which);
            for t in &hole {
                let rays = combinatorics::every_ray(&planes, p, *t, &outer).unwrap();
                assert!(!rays.is_empty(), "{which}: no clear ray");
                assert!(rays.iter().all(|&x| x), "{which}: {rays:?}");
            }
        }
    }

    #[test]
    fn a_loop_is_placed_by_where_it_is_not_by_how_it_winds() {
        // Containment is about **where** a ring is, never which way it runs. A hole ring is stored
        // clockwise about the face normal and its outer ring counter-clockwise, and neither
        // direction may enter the answer: reversing either must change nothing.
        //
        // And containment is **not symmetric** — the classic way to get this wrong is a test that
        // only ever asks it one way round.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        let (mut rev_outer, mut rev_hole) = (outer.clone(), hole.clone());
        rev_outer.reverse();
        rev_hole.reverse();
        for t in &hole {
            assert!(combinatorics::point_in_ring(&planes, p, *t, &outer).unwrap());
            assert!(
                combinatorics::point_in_ring(&planes, p, *t, &rev_outer).unwrap(),
                "reversing the outer ring must not move the hole"
            );
        }
        for t in &outer {
            assert!(!combinatorics::point_in_ring(&planes, p, *t, &hole).unwrap());
            assert!(
                !combinatorics::point_in_ring(&planes, p, *t, &rev_hole).unwrap(),
                "nor may reversing the hole swallow the outer ring"
            );
        }
    }

    #[test]
    fn every_clear_ray_agrees() {
        // The ring is simple, so the parity cannot depend on the ray. `every_ray` returns one
        // answer per usable candidate and they must be unanimous; a disagreement means the ray
        // choice leaked into the result.
        //
        // (Its ancestor also pinned that *half* the candidates were unusable — that count came
        // from the retired staple fixture, whose loop and arc shared a plane. The holed L cap has
        // no such sharing, so only the unanimity survives the move.)
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        let rays = combinatorics::every_ray(&planes, p, hole[0], &outer).unwrap();
        assert!(!rays.is_empty(), "at least one candidate is clear");
        assert!(rays.iter().all(|&x| x), "and they agree: inside — {rays:?}");
    }

    #[test]
    fn a_ring_inside_a_ring_is_what_nesting_looks_like() {
        // `nested_loops` has no operand in the suite that produces it — a polyhedral torus would.
        // The detector can still be aimed at real geometry: a holed face *is* a ring inside a ring,
        // which fires exactly the condition the nesting brick asks about. It is the detector under
        // test, not the fixture.
        let (planes, p, outer, hole) = holed_face_rings("dimple");
        assert!(combinatorics::point_in_ring(&planes, p, hole[0], &outer).unwrap());
        assert!(!combinatorics::point_in_ring(&planes, p, outer[0], &hole).unwrap());
    }

    #[test]
    fn cut_l_staple() {
        // A loop beside an arc, resolved. The cap's kept region is the hexagon minus the
        // corner bite, and the near leg's rectangle sits inside it — so the loop is a
        // **hole** of that region. The winding says clockwise, which is only a cross-check
        // now: containment decided.
        //
        // `V_∩ = 0.5·0.5·0.65 + 0.27·0.55 = 0.311`, the two legs' parts inside the L.
        let (mut m, l, st) = l_and_staple();
        let r = boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 2.689).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 15.755).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn fuse_l_staple() {
        // `3 + 0.7605 − 0.311`. The staple measures `1.17` in section, `0.65` deep.
        let (mut m, l, st) = l_and_staple();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, st).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (3.0 + 0.7605 - 0.311)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 16.72).abs() < 1e-9, "area {}", props.area);
    }

    #[test]
    fn cut_staple_by_l() {
        // ★ The same loop, on the same face, is now an **island**. Swap the operands and the
        // cap's kept region becomes the corner bite alone; the loop lies outside it, in a
        // dropped region, so its interior is what survives.
        //
        // Cell 3f-3's counterexample made flesh: a hole and an island wind oppositely without
        // nesting, so no winding could have told these two apart. Position did.
        //
        // `0.7605 − 0.311`, and the three volumes close inclusion–exclusion exactly.
        let (mut m, l, st) = l_and_staple();
        let r = boolean_one(&mut m, BoolKind::Cut, st, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (0.7605 - 0.311)).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 4.68).abs() < 1e-9, "area {}", props.area);
    }

    /// The L with a stub rising out of its top face, footprint strictly inside that
    /// face. Unlike the rod of `drill_through_the_l`, the stub enters the L from within,
    /// so each of its vertical edges crosses exactly one face and its bottom ring stays
    /// inside — one chord, no tunnel, and the seam loop is the whole story.
    fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let stub = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 0.5]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        (m, l, stub)
    }

    #[test]
    fn cut_the_stub_by_the_l_leaves_an_island_face() {
        // Swap the operands of `cut_blind_dimple` and the same loop lands on a face whose
        // boundary is *all* dropped: the kept region is the loop's interior alone. That
        // face has no `∂f` at all — its outer loop *is* the seam ring, four `Discovered`
        // vertices and nothing else. The answer is the `0.4 × 0.4 × 0.5` box above `z = 1`.
        //
        // This is where `flip` first meets a discovered hole ring. It is wound CCW about the
        // L's `+z`, keeping the material (inside the stub) on
        // its left; `flip` reverses it and the face becomes the box's downward-facing
        // floor. Nothing but `validate` and the signed mesh volume can see that go wrong.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Cut, stub, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.08).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn fuse_the_stub_and_the_l() {
        // Not a new branch — the same hole, on the same face of the L, reached with the
        // operands the other way round. `Fuse` keeps both outsides and flips neither, so
        // what this pins is that the answer does not depend on which solid is `a`: the
        // hole now lands on B, and 3.08 is 3.08.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Fuse, stub, l).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn cut_blind_dimple() {
        // The stub's footprint never reaches the L's top-face boundary, so the seam is
        // a closed ring in the face interior. It is the face's inner loop, and the
        // result is a blind pocket: `3 − 0.4² × 0.5`.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 - 0.08)).abs() < 1e-9, "volume {vol}");
    }

    /// **`Origin` no longer tells result faces apart.** The arrangement names every vertex
    /// it emits by the three planes meeting there, so an operand corner the cut never touched
    /// comes back as `Discovered`, exactly like a seam vertex. Nothing carries over as
    /// `Constructed` (the arrangement builds no vertex from an original handle).
    ///
    /// This is a contract, not a curiosity: `pipeline.rs`'s island test selected a face by
    /// "all its vertices are `Discovered`", which was unique under the old engine and is
    /// now true of *every* face. It flipped the wrong face and only the last assertion
    /// noticed. Selecting a face by provenance is what this locks out.
    ///
    /// The subject is `Cut(l, stub)` — the **holed** result, so `face_half_edges` walks
    /// `inner` rings too (`count_discovered` walks only `outer` and would miss them).
    ///
    /// **Unrotated only.** Rotating a boolean result re-marks these vertices `Rotated` over
    /// a `Discovered` base — that is `transform_rotate_boolean_result_keeps_discovered_base`,
    /// and this lock must not be read as contradicting it.
    #[test]
    fn an_unrotated_boolean_names_every_vertex_by_its_plane_triple() {
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();

        let mut seen = std::collections::HashSet::new();
        let mut holed = 0;
        for sh in solid_shell_handles(&m, r) {
            for &fh in &m.shells.get(sh).faces {
                let face = m.faces.get(fh);
                holed += usize::from(!face.inner.is_empty());
                for he in face_half_edges(face) {
                    for vh in m.edges.get(he.edge).bounds.into_iter().flatten() {
                        if !seen.insert(vh) {
                            continue;
                        }
                        assert!(
                            matches!(
                                m.vertices.get(vh).origin,
                                Origin::Discovered {
                                    definition: VertexDef::ThreePlane(_),
                                    ..
                                }
                            ),
                            "vertex {:?} is {:?}, not a plane triple",
                            m.vertices.get(vh).point.as_array(),
                            m.vertices.get(vh).origin
                        );
                    }
                }
            }
        }
        assert_eq!(holed, 1, "the blind dimple leaves exactly one holed face");
        assert_eq!(seen.len(), 20, "the L's 12 corners + the dimple's 8");
    }

    #[test]
    fn a_flipped_hole_loop_is_caught() {
        // `loop_orient_mismatch` cannot see a loop flipped as a whole, and `f2`'s golden
        // pins the derivation — but the two downstream detectors must actually fire on
        // *this* shape, not merely exist. They are different in kind: `validate` sees a
        // rim edge used twice the same way, `tessellate` sees a hole wound like its
        // outer ring. Volume and OCCT see neither: `props` sums `|area|`.
        //
        // `Store` is append-only, so the face cannot be edited. Push a replacement with
        // the hole reversed, swap it into a fresh shell and solid, and move the live
        // handle: the old face falls out of `reachable()`.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        m.rebuild_adjacency();

        // The dimple is blind, so exactly one face carries a hole. Counted, not assumed:
        // its sibling in `pipeline.rs` said "the only face …" in prose, used `find`, and
        // silently flipped a different face once the engine stopped making the predicate
        // unique.
        let faces = m.shells.get(m.solids.get(r).outer).faces.clone();
        let holed_faces: Vec<_> = faces
            .iter()
            .copied()
            .filter(|&f| !m.faces.get(f).inner.is_empty())
            .collect();
        assert_eq!(holed_faces.len(), 1, "the L's top face carries the hole");
        let holed = holed_faces[0];
        let f = m.faces.get(holed).clone();
        let mut hole = f.inner[0].clone();
        hole.half_edges.reverse();
        for he in &mut hole.half_edges {
            he.forward = !he.forward;
        }
        let bad = m.faces.push(Face {
            inner: vec![hole],
            ..f
        });
        let swapped = faces
            .iter()
            .map(|&x| if x == holed { bad } else { x })
            .collect();
        let shell = m.shells.push(Shell { faces: swapped });
        let solid = m.push_solid(Solid {
            outer: shell,
            cavities: vec![],
        });
        m.live_solids = vec![solid];
        m.rebuild_adjacency();

        let vs = nacre_validate::validate(&m);
        assert!(
            vs.iter()
                .any(|v| matches!(v, nacre_validate::Violation::NonOpposedEdge { .. })),
            "validate stayed quiet: {vs:?}"
        );
        assert!(matches!(
            nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()),
            Err(nacre_tess::TessError::HoleWinding)
        ));
    }

    #[test]
    fn fuse_blind_dimple() {
        // The `Fuse` counterpart — a boss on the L — closing the inclusion–exclusion:
        // `V_L + V_stub − V_overlap`. Both put a hole in the same face of the L.
        let (mut m, l, stub) = l_and_dimple();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, stub).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (3.0 + 0.16 - 0.08)).abs() < 1e-9, "volume {vol}");
    }

    /// The unit cube with a 0.4-square pocket, 0.5 deep, in its top face: the void
    /// is `[0.3,0.7]² × [0.5,1]` and the solid measures `1 − 0.16·0.5 = 0.92`. Its
    /// lid is the only face in the suite that carries an inner loop.
    fn pocketed_cube() -> (Model, Handle<Solid>) {
        let (mut m, top) = cube_with_top();
        let OpOutput::PocketOnFace { solid, .. } =
            apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        (m, solid)
    }

    /// The first time a boolean result is fed back as an operand: the overlap box
    /// (Discovered corners) stacked on a third box merges through the coincident-
    /// interface path — which runs `is_convex` on that Discovered-cornered
    /// operand. Volume is the sum and the shell stays closed.
    #[test]
    fn a_boolean_result_stacks_as_an_operand() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let c = boolean_one(&mut m, BoolKind::Common, a, b).unwrap(); // [1,2]³
        m.rebuild_adjacency();
        let d = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 2.0]),
            Point3::from_array([2.0, 2.0, 3.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, c, d).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 2.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    /// A holed operand already survives the arrangement — it names holes as inner rings rather
    /// than guarding against them, so the pocket rides through untouched. Nothing tested it.
    /// Pin it before the guard comes down.
    #[test]
    fn containment_boolean_already_keeps_a_pocket() {
        for (kind, want) in [(BoolKind::Cut, 0.919), (BoolKind::Fuse, 0.92)] {
            let (mut m, pc) = pocketed_cube();
            let bx = m.add_cuboid(
                Point3::from_array([0.05, 0.05, 0.05]),
                Point3::from_array([0.15, 0.15, 0.15]),
            );
            let r = boolean_one(&mut m, kind, pc, bx).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(vs.is_empty(), "{kind:?} {vs:?}");
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!((vol - want).abs() < 1e-9, "{kind:?} volume {vol}");
        }
    }

    /// An edge of one convex solid, threading the other, severs it. The bar runs through the
    /// cube and out both ends, so `Cut(bar, cube)` leaves the bar in two 1×1×1 stubs — two solids
    /// (cell 0.4), each a clean genus-0 box (`validate` clean, pieces share nothing). The convex
    /// path rejected this as `poke_through`; the seam path returns both pieces. (Before cell 0.4
    /// this was `disconnected_result`; cell 3e-3's non-convex sibling was `Cut(rod, L)`.)
    #[test]
    fn a_convex_cut_severs_its_operand() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, bar, a).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        for &s in &solids {
            assert!((nacre_props::mass_props(&m, s).unwrap().volume - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn cut_by_a_box_inside_the_pocket_is_a_no_op() {
        // The box sits wholly in the void, so the solids are disjoint and `A − B = A`.
        // Reading the lid as filled instead classified the box's eight corners five
        // Inside and three Outside, and the seam-free path's own `debug_assert`
        // ("classification must be consistent per solid") caught it.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.6]),
            Point3::from_array([0.6, 0.6, 0.9]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!(m.solids.get(r).cavities.is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.92).abs() < 1e-9, "volume {vol}");
    }

    #[test]
    fn cut_with_a_hollow_operand_far_from_the_void() {
        // A cavitied operand whose seam misses the void: the corner cut is far from
        // the [1,2]³ void, so the void is carried through and preserved (cell (5c-in)).
        // The seam front-end walks all shells, so the void no longer silently vanishes
        // (it used to read as convex and return vol 26.875, cavities 0). Now: correct
        // 25.875 (27 − 1 void − 0.125 corner) with the cavity intact.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let cutter = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, cutter).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn blind_hole_drills_into_a_void() {
        // A cut reaching *into* a void: a stub drilled from below the hollow box up into
        // its void. The void loses its enclosure and merges with the outer shell (an
        // open pocket, cavities 0). The seam machinery (all-shells) + the (5c) component
        // split reconstruct it exactly: remove the channel [1.4,1.6]²×[0,1]=0.04 through
        // the floor, void interior removes nothing ⇒ 26 − 0.04 = 25.96.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        m.rebuild_adjacency();
        let stub = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, stub).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.96).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0); // the void opened to outside
    }

    #[test]
    fn a_tunnel_drilled_through_a_void() {
        // A cut passing all the way through the void (bottom to top). Removes the floor
        // and ceiling channels [1.4,1.6]²×([0,1]∪[2,3]) = 0.08; the void interior removes
        // nothing ⇒ 26 − 0.08 = 25.92. Result is a genus-1 solid (a straight tunnel),
        // one shell, no cavity.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let tunnel = m.add_cuboid(
            Point3::from_array([1.4, 1.4, -0.5]),
            Point3::from_array([1.6, 1.6, 3.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, tunnel).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.92).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    #[test]
    fn a_slab_splits_a_hollow_box_into_two() {
        // A slab cut through the whole box (and its void) severs it into two solids (cell 0.4).
        // The slab spans the full cross-section, so it opens the void — both pieces are
        // cavity-free. (Before cell 0.4 this was rejected as `DISCONNECTED_RESULT`.)
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.4, -0.5]),
            Point3::from_array([3.5, 1.6, 3.5]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        for &s in &solids {
            assert_eq!(m.solids.get(s).cavities.len(), 0);
        }
        // Hollow 26 (= 27 − 1 void); the slab removes 1.6 of material (8 area × 0.2 thick).
        let vol: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((vol - 24.4).abs() < 1e-9, "total volume {vol}");
    }

    /// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
    /// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
    /// cavity, the x>2 piece is solid — two outward shells *and* one inward. Which outer owns the
    /// cavity needs a containment test we do not have yet, so it is honestly rejected
    /// (`SEVERED_WITH_CAVITY`) rather than mis-assembled. This is the firing test the guard is
    /// born with (design.md); n0 measured the two-outward-plus-one-inward component split.
    #[test]
    fn severed_with_cavity_is_rejected() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        // Void near the x-low side, clear of the x=2 cut.
        let inner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 2.5]),
        );
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void)
        // and x>2 (solid).
        let slab = m.add_cuboid(
            Point3::from_array([2.0, -1.0, -1.0]),
            Point3::from_array([2.2, 4.0, 4.0]),
        );
        assert_rejects(
            || boolean(&mut m, BoolKind::Cut, hollow, slab),
            tag::SEVERED_WITH_CAVITY,
        );
    }

    // A hollow operand in a COPLANAR contact — the combination nothing covered until now. The
    // cavity goldens above all take the seam path (transversal cuts) and every coplanar golden uses
    // solid operands, so the intersection of the two was a blind spot, and the coplanar driver had
    // never received the all-shell patch the seam front-end got in cell (5c-in). It emitted only
    // outer-shell faces, so the void vanished: the fuse read 27.0625 — the *un-hollowed* cube plus
    // the boss — with `cavities: 0` and a clean `validate`, because what remained was still a
    // closed shell. Silent-wrong, invisible to every guard. Now the driver walks all shells.
    #[test]
    fn a_hollow_part_takes_a_coplanar_boss() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        // Top-flush boss on z=3: a genuine coplanar contact, clear of the void's planes.
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 3.0]),
            Point3::from_array([0.75, 0.75, 4.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // 27 − 1 void + 0.25·0.25·1 boss.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 26.0625).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
    }

    #[test]
    fn a_hollow_part_takes_a_coplanar_pocket() {
        // The Cut twin of the boss case: a top-flush pocket sunk into a hollow part. Same blind
        // spot, same silent-wrong before the fix (26.96875 with the void gone).
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let tool = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 2.5]),
            Point3::from_array([0.75, 0.75, 3.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, tool).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        // 27 − 1 void − 0.25·0.25·0.5 pocket.
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.96875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
    }

    #[test]
    fn a_coplanar_boss_over_a_void_plane_is_solved() {
        // The boss straddles the void's x=1 and y=1 planes (it spans [0.9,1.1]²), so classifying
        // the void's walls against it is not the clean whole-face case. This was pinned as an
        // honest reject with the standing instruction that a future change may "turn it into a
        // correct 26.04, never a silent answer" — and family #3 did: once one geometric plane is
        // one class whatever the two faces' sizes, the arrangement solves it. 26 (hollow) + 0.04
        // (boss); area 60 + 0.8 (boss sides) + 0.04 (its top) − 0.04 (its footprint); void intact.
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let boss = m.add_cuboid(
            Point3::from_array([0.9, 0.9, 3.0]),
            Point3::from_array([1.1, 1.1, 4.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Fuse, hollow, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - 26.04).abs() < 1e-9,
            "volume {}",
            props.volume
        );
        assert!((props.area - 60.8).abs() < 1e-9, "area {}", props.area);
        assert_eq!(m.solids.get(r).cavities.len(), 1, "the void survives");
    }

    #[test]
    fn a_hollow_part_takes_a_second_far_cut() {
        // The headline: keep cutting a part after it is hollow. A bore at the corner
        // opposite the void — the void survives, and the result is fed back as an
        // operand (chaining past the first cavity-producing op, which the door used to
        // block).
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        let bore = m.add_cuboid(Point3::from_array([-0.5; 3]), Point3::from_array([0.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, hollow, bore).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 25.875).abs() < 1e-9, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn clockwise_input_is_auto_corrected() {
        // The square wound CW; auto-CCW makes it a valid cube anyway.
        let cw = Profile2d {
            points: vec![p2(0.0, 1.0), p2(1.0, 1.0), p2(1.0, 0.0), p2(0.0, 0.0)],
        };
        let m = replay(&[extrude_op(cw, 1.0)]).unwrap();
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(m.faces.len(), 6);
    }

    #[test]
    fn replay_is_deterministic() {
        let log = vec![extrude_op(square(), 1.0)];
        let m1 = replay(&log).unwrap();
        let m2 = replay(&log).unwrap();
        let pts = |m: &Model| {
            m.vertices
                .iter()
                .map(|(_, v)| v.point.as_array())
                .collect::<Vec<_>>()
        };
        assert_eq!(pts(&m1), pts(&m2));
        assert_eq!(m1.edges.len(), m2.edges.len());
        assert_eq!(m1.faces.len(), m2.faces.len());
    }

    #[test]
    fn two_extrudes_make_two_solids() {
        let far = SketchPlane {
            origin: Point3::from_array([5.0, 0.0, 0.0]),
            ..SketchPlane::world_xy()
        };
        let log = vec![
            extrude_op(square(), 1.0),
            Operation::Extrude {
                plane: far,
                profile: square(),
                dist: 1.0,
            },
        ];
        let m = replay(&log).unwrap();
        assert_eq!(m.solids.len(), 2);
        assert!(nacre_validate::validate(&m).is_empty());
    }

    #[test]
    fn degenerate_inputs_are_rejected() {
        let plane = SketchPlane::world_xy();
        let two = Profile2d {
            points: vec![p2(0.0, 0.0), p2(1.0, 0.0)],
        };
        assert_eq!(
            apply(
                &mut Model::new(),
                &Operation::Extrude {
                    plane,
                    profile: two,
                    dist: 1.0
                }
            ),
            Err(OpError::DegenerateProfile)
        );
        assert_eq!(
            apply(&mut Model::new(), &extrude_op(square(), 0.0)),
            Err(OpError::NonPositiveDistance)
        );
        let dup = Profile2d {
            points: vec![p2(0.0, 0.0), p2(0.0, 0.0), p2(1.0, 1.0)],
        };
        assert_eq!(
            apply(
                &mut Model::new(),
                &Operation::Extrude {
                    plane,
                    profile: dup,
                    dist: 1.0
                }
            ),
            Err(OpError::DegenerateGeometry)
        );
    }

    /// Extrude a unit cube and return `(model, top face handle)`.
    fn cube_with_top() -> (Model, Handle<Face>) {
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude_op(square(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        let top = faces[1]; // base, top, sides…
        (m, top)
    }

    fn small_square() -> Profile2d {
        Profile2d {
            points: vec![p2(-0.2, -0.2), p2(0.2, -0.2), p2(0.2, 0.2), p2(-0.2, 0.2)],
        }
    }

    // A single-edge overhang footprint on the cube top: world x∈[0.25,0.75], y∈[-0.25,0.75]
    // (overhangs the y=0 edge), area 0.5. (Frame maps local [px,py] → world (0.5+py, 0.5−px).)
    fn edge_overhang_profile() -> Profile2d {
        Profile2d {
            points: vec![
                p2(-0.25, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(-0.25, 0.25),
            ],
        }
    }

    // A spanning slab: world x∈[0.25,0.75], y∈[-0.25,1.25] (crosses both y edges), area 0.75.
    fn spanning_slab_profile() -> Profile2d {
        Profile2d {
            points: vec![
                p2(-0.75, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(-0.75, 0.25),
            ],
        }
    }

    #[test]
    fn pad_an_overhanging_boss() {
        // The profile reaches past one face edge: part of the boss sits on the face, part
        // cantilevers into the air. Relaxing the containment gate routes it to the overhang Fuse
        // sidecar (Ok here proves the routing — a contained-only pad would reject). The boss lives
        // wholly above z=1, so vol = cube 1 + footprint 0.5 · dist 1 = 1.5.
        let (mut m, top) = cube_with_top();
        let OpOutput::PadOnFace { solid, top_face } =
            apply(&mut m, &pad_op(top, edge_overhang_profile(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.5).abs() < 1e-12);
        assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
    }

    #[test]
    fn pad_a_spanning_slab_boss() {
        // A slab crossing the whole face (overhangs two opposite edges). vol = 1 + 0.75 · 1 = 1.75.
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pad_op(top, spanning_slab_profile(), 1.0)).unwrap();
        let OpOutput::PadOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 1.75).abs() < 1e-12);
    }

    #[test]
    fn pocket_an_edge_slot() {
        // A blind pocket whose footprint overhangs one edge — an edge slot open to the side.
        // Only the on-face part (world x[0.25,0.75]×y[0,0.75] = 0.375) carves: 1 − 0.375·0.5 = 0.8125.
        let (mut m, top) = cube_with_top();
        let OpOutput::PocketOnFace { solid, bottom_face } =
            apply(&mut m, &pocket_op(top, edge_overhang_profile(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.8125).abs() < 1e-12);
        assert!(m.reachable().faces.contains(&bottom_face)); // slot floor recovered
    }

    #[test]
    fn pocket_a_slab_channel() {
        // A blind channel crossing the whole face (breaches two opposite walls). On-face carve
        // world x[0.25,0.75]×y[0,1] = 0.5: 1 − 0.5·0.5 = 0.75.
        let (mut m, top) = cube_with_top();
        let out = apply(&mut m, &pocket_op(top, spanning_slab_profile(), 0.5)).unwrap();
        let OpOutput::PocketOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        assert!((nacre_props::mass_props(&m, solid).unwrap().volume - 0.75).abs() < 1e-12);
    }

    #[test]
    fn pad_overhang_off_the_face_is_rejected() {
        // A footprint that does not touch the face at all. The boolean is not what fails here — it
        // fuses the two into a base plus a detached boss, which is the right answer (see
        // `a_touchless_boss_fuses_into_two_solids`). What breaks is the *pad's* premise, so the
        // error names that, and the model the caller is left holding is the one it started with.
        let (mut m, top) = cube_with_top();
        let far = Profile2d {
            points: vec![p2(1.8, 1.8), p2(2.2, 1.8), p2(2.2, 2.2), p2(1.8, 2.2)],
        };
        let before = m.live_solids.clone();
        assert_eq!(
            apply(&mut m, &pad_op(top, far, 0.3)),
            Err(OpError::PadMissesFace)
        );
        let (mut a, mut b) = (before, m.live_solids.clone());
        a.sort_by_key(|h| h.index());
        b.sort_by_key(|h| h.index());
        assert_eq!(a, b, "a rejected pad must leave the live model untouched");
    }

    /// **Chaining onto a fused boss.** The fuse leaves the base's `z=1` face a *ring* — a face with
    /// a hole where the boss sits — and the second boolean cuts through both. Every plane class the
    /// cut opens then meets that ring along the **hole's own edge**, which is the case that used to
    /// label inconsistently: the ring's neighbouring vertices there point *into* the hole, so
    /// reading the occupied side off a flank put the material on the wrong side of `W`. The side
    /// now comes from the ring's travel ([`arrangement::run_body_above`]), and the run leaves as its own
    /// homogeneous segment rather than being swallowed by the straddling stretch beside it.
    ///
    /// Hand volume: `1 + 0.5·0.5·1` fused, less the cutter's `0.2·0.2` column over `z ∈ [0.5, 2]`
    /// — `1.25 − 0.06 = 1.19`. The same shape is scored against OCCT by
    /// `boss_fuse_then_cut_matches_occt`, but that oracle is `#[ignore]`d, so this is the copy that
    /// runs on every `cargo test`.
    #[test]
    fn a_boss_fused_then_cut_through() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 1.0]),
            Point3::from_array([0.75, 0.75, 2.0]),
        );
        let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, bossed).unwrap().volume - 1.25).abs() < 1e-12,
            "the fused boss itself"
        );
        let cutter = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).expect("the chained cut");
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, r).unwrap().volume - 1.19).abs() < 1e-12,
            "base + boss less the drilled column: {}",
            nacre_props::mass_props(&m, r).unwrap().volume
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "a chained result is still a clean model"
        );
    }

    /// The kernel's answer for a boss that misses the face, stated on its own so nobody "fixes" the
    /// boolean to reject it: fusing two solids that do not touch **is** two solids, and both are
    /// whole. Only `pad` refuses that outcome, because a pad is defined as material joined to a face
    /// (`pad_overhang_off_the_face_is_rejected`).
    #[test]
    fn a_touchless_boss_fuses_into_two_solids() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // Shares the z = 1 plane class with the base's top, but sits far away in x/y — so the plane
        // carries two separate bodies, which is exactly what `hole_roots` used to refuse.
        let boss = m.add_cuboid(
            Point3::from_array([1.8, 1.8, 1.0]),
            Point3::from_array([2.2, 2.2, 1.3]),
        );
        let solids = boolean(&mut m, BoolKind::Fuse, base, boss).unwrap();
        assert_eq!(solids.len(), 2, "disjoint operands stay two solids");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let mut vols: Vec<f64> = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!(
            (vols[0] - 0.048).abs() < 1e-12 && (vols[1] - 1.0).abs() < 1e-12,
            "both pieces whole: {vols:?}"
        );
    }

    #[test]
    fn pocket_through_overhang_is_rejected() {
        // An overhang pocket deep enough to pierce the far side is not blind — no single floor.
        // Honest reject via whichever path fires (the overhang detector declines, the seam path
        // rejects the mixed contact), mirroring `pocket_through_the_solid_is_rejected`.
        let (mut m, top) = cube_with_top();
        let got = apply(&mut m, &pocket_op(top, edge_overhang_profile(), 1.5));
        assert!(
            matches!(got, Err(OpError::Boolean(_)) | Err(OpError::PocketNotBlind)),
            "through overhang must reject honestly, got {got:?}"
        );
    }

    #[test]
    fn a_non_convex_pad_cantilevers_and_runs_flush() {
        // **Three hard properties at once**, which is what makes this footprint worth keeping. The
        // frame maps local `[px, py]` to world `(0.5 + py, 0.5 − px)`, so the L below lands on
        // `(0.25,0.75) (0.25,−0.25) (0.75,−0.25) (0.75,0.25) (1.0,0.25) (1.0,0.75)`:
        //   1. **non-convex** — the L has a reflex corner at `(0.75, 0.25)`;
        //   2. **overhanging** — `y < 0` cantilevers past the cube's `y = 0` edge;
        //   3. **flush** — the edge `x = 1.0, y∈[0.25,0.75]` lies *exactly* on the face's `x = 1`
        //      boundary, the "profile rim shares the face rim" case.
        // It used to reject because the overhang sidecars gated on convexity; that gate is gone.
        //
        // Hand-checked shape: footprint `0.25 + 0.375 = 0.625`, prism wholly above `z = 1`, so
        //   volume 1 + 0.625 = 1.625
        //   area   5 (cube minus its top) + 0.5 (top left uncovered) + 3.5 (prism sides)
        //          + 0.625 (prism cap) + 0.125 (the cantilever's underside) = 9.75
        // The underside term is the cantilever: a contained pad would not have one.
        //
        // OCCT cannot score this directly — `pad` builds its tool prism internally, and rebuilding
        // it here would lean on the same frame mapping the assertion is testing.
        let (mut m, top) = cube_with_top();
        let l_over = Profile2d {
            points: vec![
                p2(-0.25, -0.25),
                p2(0.75, -0.25),
                p2(0.75, 0.25),
                p2(0.25, 0.25),
                p2(0.25, 0.5),
                p2(-0.25, 0.5),
            ],
        };
        let OpOutput::PadOnFace { solid, top_face } =
            apply(&mut m, &pad_op(top, l_over, 1.0)).expect("the cantilevered L pad")
        else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, solid).unwrap();
        assert!((p.volume - 1.625).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 9.75).abs() < 1e-12, "area {}", p.area);
        assert!(m.reachable().faces.contains(&top_face)); // boss top cap recovered
    }

    fn pad_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PadOnFace {
            face,
            profile,
            dist,
        }
    }


    /// explicit sharing (overhaul #3): a prism built with a shared base-cap
    /// surface reuses that `Surface` handle for its flush cap, and reconciles the
    /// cap's face orientation so the materialized outward normal stays `−sweep`.
    #[test]
    fn build_prism_base_cap_reuses_shared_surface() {
        let mut m = Model::new();
        // A face-plane surface with outward normal +z (as a face on the base solid).
        let sf = m.surfaces.push(Surface::Plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                .unwrap(),
        ));
        let base_pts = [
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([0.0, 1.0, 0.0]),
        ];
        let (_prism, faces) = build_prism(
            &mut m,
            &base_pts,
            Vector3::from_array([0.0, 0.0, 1.0]),
            Some(sf),
        )
        .unwrap();
        let cap = m.faces.get(faces[0]); // base cap is pushed first
        // Shared handle (was a fresh push before overhaul #3).
        assert_eq!(cap.surface, sf, "base cap reuses the shared surface handle");
        // Orientation reconciled: materialized outward normal is −sweep (−z).
        let Surface::Plane(p) = m.surfaces.get(cap.surface) else {
            unreachable!()
        };
        let sign = match cap.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
        let materialized = p.normal() * sign;
        assert!(
            (materialized - Vector3::from_array([0.0, 0.0, -1.0])).norm() < 1e-12,
            "materialized cap normal stays −z, got {materialized:?}"
        );
    }

    /// The handle branch of `shares_or_coplanar` is load-bearing: a shared
    /// `Surface` handle reports coplanar even when the stored `plane` values are
    /// *not* geometrically coplanar (so the fallback would not fire). This is the
    /// path a referenced coplanar contact takes; on axis-aligned M5 it is redundant
    /// with the geometric test, but the branch must work for rotated frames.
    #[test]
    fn shares_or_coplanar_uses_the_handle_branch() {
        let mut m = Model::new();
        let shared = m.surfaces.push(Surface::Plane(
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0]))
                .unwrap(),
        ));
        let fh = m.faces.push(Face {
            surface: shared,
            outer: Loop { half_edges: vec![] },
            inner: vec![],
            orientation: Orientation::Forward,
        });
        let plane_x0 =
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([1.0, 0.0, 0.0]))
                .unwrap();
        let plane_z0 =
            Plane::from_point_normal(Point3::origin(), Vector3::from_array([0.0, 0.0, 1.0]))
                .unwrap();
        // Each `tri` is three NON-collinear points of its own plane. A degenerate `tri` (three
        // equal points) would make every `orient3d` vanish, so the coordinate branch would report
        // coplanar and this test would pass without the handle branch ever mattering.
        let mk = |plane, tri| FaceInfo {
            surf: shared,
            face: fh,
            plane,
            tri,
            n_out: Vector3::from_array([0.0; 3]),
            // Unread: this table only ever reaches `t_planes_coplanar`, which decides on `tri`.
            orient_sign: 1,
            tri_pt3: None,
        };
        let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
        let planes = vec![
            mk(plane_x0, [p(0., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]), // in x = 0
            mk(plane_z0, [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.)]), // in z = 0
        ];
        // Neither fallback fires: the coefficients are not proportional, and the coordinates say
        // these really are two different planes.
        assert!(!planes_coplanar(&planes[0].plane, &planes[1].plane));
        assert!(!tolerant::t_planes_coplanar(&planes, 0, 1));
        // The shared handle alone makes them coplanar-by-reference.
        assert!(shares_or_coplanar(&planes, 0, 1));
    }




    #[test]
    fn pad_step_exports() {
        // The boss (holed outer face + walls + cap) exports without error.
        let (mut m, top) = cube_with_top();
        apply(&mut m, &pad_op(top, small_square(), 0.5)).unwrap();
        let step = nacre_step::to_step(&m).expect("boss exports");
        assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
    }

    fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        }
    }

    /// A prism raised on a `(1,1,1)`-slanted sketch plane, its far cap holding a blind pocket. The
    /// cap and its two anti-parallel side walls meet in a triple whose `raw` coefficients are
    /// exactly dependent (`det = 0`); before family #3's dir-sign fix the guard read `sqrt`-rounded
    /// unit normals, called that triple non-degenerate, and the consumer aborted on `D = 0`. Now
    /// the guard reads the same coefficients the consumer does, so the arrangement runs. Volume:
    /// a `2×2` base × `2` deep block is `8`, less the `0.4²×0.5` pocket.
    #[test]
    fn a_pocket_on_a_slanted_face() {
        let plane =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        let mut m = Model::new();
        let big = Profile2d {
            points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
        };
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane,
                profile: big,
                dist: 2.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let out = apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap();
        let OpOutput::PocketOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
        assert!((vol - (8.0 - 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
    }

    /// The same slanted cap, but a boss (pad, Fuse) instead of a pocket — the sweep runs the other
    /// way, a different code path. Volume: the `8` block plus a `0.4²×0.5` stub.
    #[test]
    fn a_pad_on_a_slanted_face() {
        let plane =
            SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array([1.0, 1.0, 1.0]))
                .unwrap();
        let mut m = Model::new();
        let big = Profile2d {
            points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
        };
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane,
                profile: big,
                dist: 2.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let out = apply(
            &mut m,
            &Operation::PadOnFace {
                face: faces[1],
                profile: small_square(),
                dist: 0.5,
            },
        )
        .unwrap();
        let OpOutput::PadOnFace { solid, .. } = out else {
            unreachable!()
        };
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
        assert!((vol - (8.0 + 0.16 * 0.5)).abs() < 1e-9, "volume {vol}");
    }






    #[test]
    fn pocket_step_exports() {
        let (mut m, top) = cube_with_top();
        apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap();
        let step = nacre_step::to_step(&m).expect("pocket exports");
        assert!(step.contains("FACE_BOUND("), "the hole emits a FACE_BOUND");
    }

    proptest! {
        #[test]
        fn prop_regular_ngon_on_xy_is_clean(
            n in 3usize..8,
            r in 0.5f64..10.0,
            dist in 0.1f64..10.0,
        ) {
            let m = replay(&[extrude_op(regular_ngon(n, r), dist)]).unwrap();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            prop_assert_eq!(m.vertices.len(), 2 * n);
            prop_assert_eq!(m.faces.len(), n + 2);
        }

        #[test]
        fn prop_ngon_on_arbitrary_plane_is_clean(
            n in 3usize..8,
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in -1.0f64..1.0,
            dist in 0.1f64..10.0,
        ) {
            let normal = Vector3::from_array([nx, ny, nz]);
            prop_assume!(normal.norm() > 0.1); // skip near-zero normals
            let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
            let m = replay(&[Operation::Extrude {
                plane,
                profile: regular_ngon(n, 2.0),
                dist,
            }])
            .unwrap();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            prop_assert_eq!(m.faces.len(), n + 2);
        }

        /// A blind pocket on a randomly-slanted face: the arrangement must give a valid solid of the
        /// right volume or reject honestly — **never panic**. Drives general (non-axis) plane normals
        /// through the dir-sign guard and the `angular_order`/`turn_at` consumers that read its zeros.
        ///
        /// Was `#[ignore]`d for a residual `D = 0` panic on general normals: a triple naming one
        /// geometric plane through two coincident faces. Symmetric normals like `(1,1,1)` cleared
        /// it, `(0.446, 0.737, 0.990)` did not. **Un-ignored 2026-07-22** — naming every plane by
        /// its class made those two faces one index, so the degenerate triple can no longer form.
        #[test]
        fn pocket_on_a_random_slanted_face_is_valid_or_rejects(
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in 0.2f64..1.0, // keep the normal clear of the sketch's degenerate zero
        ) {
            let normal = Vector3::from_array([nx, ny, nz]);
            prop_assume!(normal.norm() > 0.3);
            let plane = SketchPlane::from_origin_normal(Point3::origin(), normal).unwrap();
            let mut m = Model::new();
            let big = Profile2d {
                points: vec![p2(-1.0, -1.0), p2(1.0, -1.0), p2(1.0, 1.0), p2(-1.0, 1.0)],
            };
            let OpOutput::Extrude { faces, .. } =
                apply(&mut m, &Operation::Extrude { plane, profile: big, dist: 2.0 }).unwrap()
            else { unreachable!() };
            match apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)) {
                Ok(OpOutput::PocketOnFace { solid, .. }) => {
                    m.rebuild_adjacency();
                    prop_assert!(nacre_validate::validate(&m).is_empty());
                    let vol = nacre_props::mass_props(&m, solid).unwrap().volume;
                    prop_assert!((vol - (8.0 - 0.16 * 0.5)).abs() <= 1e-9 * 8.0, "volume {}", vol);
                }
                Ok(_) => prop_assert!(false, "unexpected op output"),
                Err(_) => {} // an honest reject is acceptable; a panic is not (and would fail the test)
            }
        }

        /// A random boss on a random box stays a valid b-rep (any interior
        /// profile, any positive height).
        #[test]
        fn prop_pad_stays_valid(
            sx in 0.5f64..5.0,
            sy in 0.5f64..5.0,
            sz in 0.5f64..5.0,
            h in 0.05f64..0.15,
            dist in 0.1f64..5.0,
        ) {
            let rect = Profile2d {
                points: vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)],
            };
            let mut m = Model::new();
            let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: rect,
                dist: sz,
            }).unwrap() else { unreachable!() };
            let hole = Profile2d {
                points: vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)],
            };
            apply(&mut m, &Operation::PadOnFace { face: faces[1], profile: hole, dist }).unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
        }

        /// A random blind pocket on a random box stays valid. `dist ≤ 0.8 < sz`
        /// keeps the pocket from punching through the box (height `sz ≥ 1`).
        #[test]
        fn prop_pocket_stays_valid(
            sx in 0.5f64..5.0,
            sy in 0.5f64..5.0,
            sz in 1.0f64..5.0,
            h in 0.05f64..0.15,
            dist in 0.1f64..0.8,
        ) {
            let rect = Profile2d {
                points: vec![p2(0.0, 0.0), p2(sx, 0.0), p2(sx, sy), p2(0.0, sy)],
            };
            let mut m = Model::new();
            let OpOutput::Extrude { faces, .. } = apply(&mut m, &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: rect,
                dist: sz,
            }).unwrap() else { unreachable!() };
            let hole = Profile2d {
                points: vec![p2(-h, -h), p2(h, -h), p2(h, h), p2(-h, h)],
            };
            apply(&mut m, &Operation::PocketOnFace { face: faces[1], profile: hole, dist }).unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
        }
    }

    // ---- boolean API (M5-c3 commit 1) ----

    fn two_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        (m, a, b)
    }

    /// Outer A = [0,3]³ (volume 27) with inner B = [1,2]³ (volume 1) strictly
    /// inside it — the containment fixture (returns `(model, outer, inner)`).
    fn nested_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        (m, a, b)
    }

    /// **Every producer states its side in the label frame** — the invariant family #2 restored,
    /// swept over the whole two-solid corpus (prints the interesting classes with `--nocapture`).
    ///
    /// The arrangement states its cell labels as `[*_above, *_below]` about one direction per plane
    /// class: the class root's **stored surface normal** (`Seated{body_above}` and `emit_faces`'
    /// `flip` are written against it). `combinatorics::side_of` answers in the root's **outward** frame
    /// instead, and the two are opposite exactly when the root face is `Reversed`
    /// (`orient_sign == -1`) — which no `add_cuboid` face ever is, but a face an earlier boolean
    /// re-emitted flipped is. `graze_above` read `side_of` raw, so on a pocket wall it flipped the
    /// wrong label bit. It no longer reads a point's side at all — [`arrangement::run_body_above`] derives
    /// the occupied side from the ring's travel, and the frame term cancels there because
    /// `order_along`'s direction and the label frame are defined by the same stored normal — but
    /// this class stays the corpus's only crossed-frame witness, so it is what would catch a
    /// producer that regresses to a raw `side_of`.
    ///
    /// What this pins, measured before the fix (2026-07-22):
    /// - the pocket fixture has 5 `orient_sign == -1` classes carrying seated *and* graze segments
    ///   (4 walls + the floor); every other fixture has **no** `Reversed` root at all, which is why
    ///   the whole corpus passed with the frames crossed and why converting cannot regress it;
    /// - the four wall classes stopped at `loop_orient_mismatch`; the floor class did **not** — its
    ///   four rim edges all carry a graze, so seated and graze were wrong *together*, consistently,
    ///   and the label survived verification while being inverted (a silent wrong, not a reject);
    /// - the hand-derived sides on those classes, so a re-crossed frame fails here first.
    #[test]
    fn every_producer_states_its_side_in_the_label_frame() {
        let mut boxed: Vec<(&str, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
        macro_rules! fixture {
            ($name:ident) => {{
                let (m, a, b) = $name();
                boxed.push((stringify!($name), m, a, b));
            }};
        }
        fixture!(two_boxes);
        fixture!(nested_boxes);
        fixture!(cube_and_notch);
        fixture!(stacked_cubes);
        fixture!(l_and_corner_box);
        fixture!(l_and_reflex_box);
        fixture!(l_and_inner_box);
        fixture!(l_and_popup_box);
        fixture!(l_and_notch_bar);
        fixture!(l_and_ell_stub);
        fixture!(l_and_staple);
        fixture!(l_and_dimple);
        fixture!(l_and_rod);
        fixture!(u_and_slab);
        {
            // The pocket family: `pocketed_cube` is itself a boolean result, so its pocket walls
            // are `Reversed` faces. Box coordinates are `pocket_corner_cut`'s.
            let (mut m, pc) = pocketed_cube();
            let bx = m.add_cuboid(
                Point3::from_array([0.85, 0.85, 0.85]),
                Point3::from_array([1.15, 1.15, 1.15]),
            );
            boxed.push(("pocket_corner_cut", m, pc, bx));
        }

        let mut reversed_with_graze: Vec<String> = Vec::new();
        let mut reversed_seated_only: Vec<String> = Vec::new();
        for (name, m, a, b) in &boxed {
            let audits = arrangement::frame_audit(m, BoolKind::Cut, *a, *b).unwrap();
            for au in &audits {
                let interesting = au.orient_sign < 0 || au.failed_at.is_some();
                if !interesting {
                    continue;
                }
                let where_ = format!(
                    "{name}: wc={} pt={:?} n={:?} orient_sign={} seated={:?} graze={:?} trans={} \
                     declined={:?} failed_at={:?}",
                    au.wc,
                    au.root_point,
                    au.root_normal,
                    au.orient_sign,
                    au.seated,
                    au.grazes,
                    au.transversals,
                    au.declined,
                    au.failed_at,
                );
                println!("{where_}");
                if au.orient_sign < 0 {
                    if au.grazes.is_empty() {
                        if !au.seated.is_empty() {
                            reversed_seated_only.push(where_);
                        }
                    } else {
                        reversed_with_graze.push(where_);
                    }
                }
            }
        }
        println!("--- reversed-root classes carrying a graze (the set the fix moves) ---");
        for r in &reversed_with_graze {
            println!("  {r}");
        }
        println!("--- reversed-root classes with seated but no graze (the alternative's risk) ---");
        for r in &reversed_seated_only {
            println!("  {r}");
        }
        // The pocket fixture must keep supplying such classes, or this test has stopped exercising
        // the crossed-frame configuration and would pass vacuously.
        assert_eq!(
            reversed_with_graze.len(),
            5,
            "the pocket's 4 walls + floor are the corpus's only reversed-root classes with a graze"
        );
        assert!(
            reversed_with_graze.iter().all(|r| r.starts_with("pocket")),
            "no other fixture may have one: {reversed_with_graze:?}"
        );

        // Hand-derived sides on the pocket wall class `x = 0.7` (root = the pocket's +x wall, whose
        // outward normal points into the void, so the stored normal `+x` makes "above" the material
        // side `x > 0.7`). The wall is seated with its body above; the two side walls and the floor
        // graze it from `x < 0.7`, i.e. below. Crossed frames invert the grazes.
        //
        // The **box's top face** grazes it too, from `x > 0.7`: the pocket's opening makes that face
        // a notched region whose edge rides this plane with the material outside the pocket. That is
        // a run whose flanks *differ*, which the engine used to read as a straddling transversal —
        // this class is the corpus's only `Reversed` root, so it is also the only place the frame
        // handling of `arrangement::run_body_above` is exercised against a crossed frame: the three
        // `false` entries below are the pre-existing answers, unchanged by the new rule.
        let (m, pc, bx) = boxed
            .iter()
            .find_map(|(n, m, a, b)| (*n == "pocket_corner_cut").then_some((m, *a, *b)))
            .unwrap();
        let wall = arrangement::frame_audit(m, BoolKind::Cut, pc, bx)
            .unwrap()
            .into_iter()
            .find(|au| au.root_point == [0.7, 0.7, 1.0] && au.root_normal == [1.0, 0.0, 0.0])
            .expect("the x=0.7 pocket wall class");
        assert_eq!(wall.orient_sign, -1, "the pocket wall is a Reversed face");
        assert_eq!(
            wall.seated,
            vec![true; 4],
            "body above = material at x > 0.7"
        );
        assert_eq!(
            wall.grazes,
            vec![true, false, false, false],
            "the top face grazes from x > 0.7 (its pocket-opening edge, material outside); \
             the two side walls and the floor graze from x < 0.7"
        );
        assert_eq!(wall.failed_at, None, "the class labels consistently");
    }

    /// `pocket_corner_cut` by hand, so the pocket family keeps a regression net that runs without
    /// OCCT: the unit cube less a `0.4²×0.5` pocket is `0.92`, and the corner box `[0.85,1.15]³`
    /// bites `0.15³` of solid (it clears the pocket, whose footprint stops at `x = 0.7`).
    #[test]
    fn a_corner_cut_off_a_pocketed_cube() {
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, pc, bx).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let props = nacre_props::mass_props(&m, r).unwrap();
        assert!(
            (props.volume - (0.92 - 0.15 * 0.15 * 0.15)).abs() < 1e-12,
            "volume {}",
            props.volume
        );
        // A corner bite replaces three 0.15² squares with three more: the area is unchanged at
        // 6 − 0.16 (the lid's hole) + 0.8 (four pocket walls) + 0.16 (its floor).
        assert!((props.area - 6.8).abs() < 1e-12, "area {}", props.area);
    }


    /// Lower corner of a solid's outer-shell vertex bounding box (for translation
    /// tests: a rigid move shifts it by exactly the offset).
    fn bbox_lo(m: &Model, s: Handle<Solid>) -> [f64; 3] {
        let mut lo = [f64::INFINITY; 3];
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        let p = m.vertices.get(vh).point.as_array();
                        for k in 0..3 {
                            lo[k] = lo[k].min(p[k]);
                        }
                    }
                }
            }
        }
        lo
    }

    fn count_discovered(m: &Model, s: Handle<Solid>) -> usize {
        let mut seen = std::collections::HashSet::new();
        let mut n = 0;
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh)
                            && matches!(m.vertices.get(vh).origin, Origin::Discovered { .. })
                        {
                            n += 1;
                        }
                    }
                }
            }
        }
        n
    }

    fn test_iso() -> (nacre_scalar::Isometry, [f64; 3]) {
        use nacre_scalar::Rat;
        (
            nacre_scalar::Isometry::translation([
                Rat::new(7, 2).unwrap(),
                Rat::from_int(-4),
                Rat::from_int(11),
            ]),
            [3.5, -4.0, 11.0],
        )
    }

    /// A rational translation supersedes a cuboid: rigid, so volume/area are
    /// invariant and the bounding box shifts by exactly the offset; validate/tess/
    /// STEP all accept the moved solid, and the input drops from `live_solids`.
    #[test]
    fn transform_translate_cuboid() {
        let (iso, off) = test_iso();
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let lo0 = bbox_lo(&m, c);

        let c2 = transform(&mut m, c, &iso).unwrap();
        m.rebuild_adjacency();

        assert_eq!(m.live_solids, vec![c2], "input superseded");
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-12,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-12, "area invariant");
        let lo1 = bbox_lo(&m, c2);
        for k in 0..3 {
            assert!(
                (lo1[k] - (lo0[k] + off[k])).abs() < 1e-12,
                "bbox shifted by offset"
            );
        }
        assert!(nacre_tess::to_obj(&m).is_ok(), "moved solid tessellates");
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP"),
            "moved solid exports to STEP"
        );
    }

    /// Transforming a boolean *result* (which carries `Discovered` seam vertices)
    /// preserves those vertices' `Origin` and remaps their `ThreePlane` definition
    /// onto the moved surfaces — the count survives and validate stays clean, so the
    /// definition was not silently downgraded to `Constructed`.
    #[test]
    fn transform_translate_preserves_discovered_definition() {
        let (iso, _) = test_iso();
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        let before = nacre_props::mass_props(&m, r).unwrap().volume;
        let disc = count_discovered(&m, r);
        assert!(
            disc > 0,
            "the Cut result must have Discovered seam vertices"
        );

        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(
            count_discovered(&m, r2),
            disc,
            "Discovered vertices preserved"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-12, "volume invariant");
    }

    /// Replay determinism (DNA 3): the same construction + transform reproduces the
    /// same geometry and the same handle down to the index.
    #[test]
    fn transform_is_deterministic() {
        let (iso, _) = test_iso();
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c2 = transform(&mut m, c, &iso).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }

    /// The `Transform` op flows through `apply`, superseding via the op dispatch.
    #[test]
    fn transform_op_applies() {
        let (iso, _) = test_iso();
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: iso,
            },
        )
        .unwrap();
        match out {
            OpOutput::Transform { solid } => assert_eq!(m.live_solids, vec![solid]),
            other => panic!("expected Transform output, got {other:?}"),
        }
    }

    /// A genuinely tilted rigid rotation: 30° about Z through the rational axis
    /// point (1,1,0). Non-90° and non-axis-aligned, so it exercises the Rotated
    /// origin and the boolean reject guard (unlike the 90° family, which stays exact).
    fn rot30() -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        })
    }

    /// Distinct outer-shell vertex points of a solid (dedup by handle).
    fn outer_points(m: &Model, s: Handle<Solid>) -> Vec<[f64; 3]> {
        let mut seen = std::collections::HashSet::new();
        let mut pts = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            pts.push(m.vertices.get(vh).point.as_array());
                        }
                    }
                }
            }
        }
        pts
    }

    /// A non-90° rotation genuinely tilts the solid: rigid (volume/area invariant),
    /// validate/tess/STEP clean, a known corner lands at its exact rotated image, the
    /// vertices carry `Origin::Rotated` (`solid_is_rotated`), and a boolean against it now
    /// runs (a *mixed*-rotation cut: rotated `c2` minus an axis-aligned `d` it contains, so
    /// `d` becomes a cavity — overhaul 3d-i retired the `ROTATED_UNSUPPORTED` entry guard).
    #[test]
    fn transform_rotate_cuboid_tilts_and_cuts() {
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let c2 = transform(&mut m, c, &rot30()).unwrap();
        m.rebuild_adjacency();

        assert_eq!(m.live_solids, vec![c2], "input superseded");
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-9,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-9, "area invariant");
        assert!(nacre_tess::to_obj(&m).is_ok(), "rotated solid tessellates");
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP"),
            "rotated solid exports to STEP"
        );

        assert!(solid_is_rotated(&m, c2), "vertices carry Origin::Rotated");

        // Corner (0,0,0) rotates about pivot (1,1) by 30°: dx=dy=-1, so
        // x' = 1 - cos30 + sin30, y' = 1 - sin30 - cos30, z' = 0.
        let (c30, s30) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let want = [1.0 - c30 + s30, 1.0 - s30 - c30, 0.0];
        let pts = outer_points(&m, c2);
        assert!(
            pts.iter()
                .any(|p| p.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9)),
            "corner (0,0,0) rotated to its exact image {want:?}; got {pts:?}"
        );

        // A mixed-rotation cut now runs: the axis-aligned `d` sits inside the rotated `c2`, so
        // `Cut(c2, d)` leaves `c2` with `d` carved out as a cavity (volume 24 − 1 = 23).
        let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "mixed cut is valid"
        );
        assert_eq!(
            m.solids.get(r).cavities.len(),
            1,
            "the contained box is a cavity"
        );
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 23.0).abs() < 1e-9, "volume {vol}");
    }

    /// A 90° rotation about Z is axis-aligned and exact: the solid stays
    /// `Constructed` (`solid_is_rotated` false), volume is exact, and boolean is
    /// still allowed — a following cut succeeds and validates.
    #[test]
    fn transform_rotate_90_is_exact_and_allows_boolean() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let rot90 = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
        });
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &rot90).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert!(!solid_is_rotated(&m, c2), "90° stays exact (Constructed)");
        assert_eq!(
            nacre_props::mass_props(&m, c2).unwrap().volume,
            24.0,
            "exact volume"
        );

        // c rotated 90° about origin occupies x∈[-3,0], y∈[0,2], z∈[0,4].
        // Cut with d = [-1,0.5,1]-[0.5,1.5,2]: overlap volume 1 → 24 − 1 = 23.
        let d = m.add_cuboid(
            Point3::from_array([-1.0, 0.5, 1.0]),
            Point3::from_array([0.5, 1.5, 2.0]),
        );
        let r = boolean(&mut m, BoolKind::Cut, c2, d).expect("exact rotation → boolean allowed");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r[0]).unwrap().volume;
        assert!((vol - 23.0).abs() < 1e-9, "cut volume {vol}");
    }

    /// Rotating a boolean *result* marks its `Discovered` seam vertices `Rotated`
    /// over the pre-rotation vertex as base: validate stays clean, volume is
    /// invariant, and at least one Rotated base is a Discovered vertex (the seam
    /// definition is preserved through the rotation, not downgraded).
    #[test]
    fn transform_rotate_boolean_result_keeps_discovered_base() {
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        assert!(
            count_discovered(&m, r) > 0,
            "Cut result has Discovered seams"
        );
        let before = nacre_props::mass_props(&m, r).unwrap().volume;

        let r2 = transform(&mut m, r, &rot30()).unwrap();
        m.rebuild_adjacency();

        assert!(nacre_validate::validate(&m).is_empty());
        assert!(
            solid_is_rotated(&m, r2),
            "rotated result carries Rotated origin"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");

        let sh = m.solids.get(r2).outer;
        let mut found_disc_base = false;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if let Origin::Rotated { base, .. } = m.vertices.get(vh).origin {
                            if matches!(m.vertices.get(base).origin, Origin::Discovered { .. }) {
                                found_disc_base = true;
                            }
                        }
                    }
                }
            }
        }
        assert!(
            found_disc_base,
            "a Rotated vertex's base is its Discovered seam vertex"
        );
    }

    /// Replay determinism (DNA 3): the same construction + rotation reproduces the
    /// same geometry and the same handle.
    #[test]
    fn transform_rotate_is_deterministic() {
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c2 = transform(&mut m, c, &rot30()).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }

    /// A rotation `Transform` flows through `apply` and marks the result Rotated.
    #[test]
    fn transform_rotate_op_applies() {
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: rot30(),
            },
        )
        .unwrap();
        let OpOutput::Transform { solid } = out else {
            panic!("expected Transform output, got {out:?}");
        };
        assert_eq!(m.live_solids, vec![solid]);
        assert!(solid_is_rotated(&m, solid));
    }

    fn boundary_verts(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
        let mut seen = std::collections::HashSet::new();
        let mut vs = Vec::new();
        let sh = m.solids.get(s).outer;
        for &fh in &m.shells.get(sh).faces {
            for he in &m.faces.get(fh).outer.half_edges {
                if let Some(bd) = m.edges.get(he.edge).bounds {
                    for vh in bd {
                        if seen.insert(vh) {
                            vs.push(vh);
                        }
                    }
                }
            }
        }
        vs
    }

    fn rot_iso(axis: nacre_scalar::Axis, deg: i128) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
        Isometry::rotation(SRot {
            axis,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    }

    /// Chain: (leaf, base_is_rotated, node_count, axes-root-to-leaf) for the first
    /// Rotated boundary vertex of `s`.
    fn forest_probe(m: &Model, s: Handle<Solid>) -> Option<(bool, usize, Vec<nacre_scalar::Axis>)> {
        let vh = *boundary_verts(m, s).first()?;
        let Origin::Rotated { base, rotation } = m.vertices.get(vh).origin else {
            return None;
        };
        let base_is_rotated = matches!(m.vertices.get(base).origin, Origin::Rotated { .. });
        let mut axes = Vec::new();
        let mut cur = Some(rotation);
        while let Some(h) = cur {
            let n = m.rotations.get(h);
            axes.push(n.axis);
            cur = n.parent;
        }
        axes.reverse();
        Some((base_is_rotated, axes.len(), axes))
    }

    fn translate_iso(off: [i128; 3]) -> nacre_scalar::Isometry {
        use nacre_scalar::{Isometry, Rat};
        Isometry::translation([
            Rat::from_int(off[0]),
            Rat::from_int(off[1]),
            Rat::from_int(off[2]),
        ])
    }

    fn rigid_iso(axis: nacre_scalar::Axis, deg: i128, off: [i128; 3]) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation as SRot};
        Isometry::rigid(
            SRot {
                axis,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            },
            [
                Rat::from_int(off[0]),
                Rat::from_int(off[1]),
                Rat::from_int(off[2]),
            ],
        )
    }

    /// A boolean produces `Discovered` seam vertices with tol 0 (exact axis-aligned
    /// intersections). Before exact quadrantal realization, rotating them exactly 90°
    /// left an ~8e-17 f64 residual that exceeded tol 0 → `VertexOffSurface`. Now the
    /// rotation is exact, so the residual stays 0 and validate is clean — both for a
    /// pure 90° rotation and for a rigid 90°+translation (the offset cancels in
    /// vertex−plane, so it does not reintroduce a residual).
    #[test]
    fn boolean_result_rotated_90_validates() {
        use nacre_scalar::Axis;
        for iso in [rot_iso(Axis::Z, 90), rigid_iso(Axis::Z, 90, [5, -3, 2])] {
            let (mut m, a, b) = two_boxes();
            let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
            let before = nacre_props::mass_props(&m, r).unwrap().volume;
            let r2 = transform(&mut m, r, &iso).unwrap();
            m.rebuild_adjacency();
            let vs = nacre_validate::validate(&m);
            assert!(
                vs.is_empty(),
                "exact 90° realization → validate clean: {vs:?}"
            );
            let after = nacre_props::mass_props(&m, r2).unwrap().volume;
            assert!((after - before).abs() < 1e-12, "volume invariant");
        }
    }

    /// A 90°-family rotation lands a cuboid's vertices exactly on the axis-aligned grid
    /// (no ~6e-17 spurious offset): the corner (2,3,4) rotated 90° about Z maps to
    /// exactly (-3,2,4).
    #[test]
    fn rotate_90_lands_vertices_exactly() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c2 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
        m.rebuild_adjacency();
        let pts = outer_points(&m, c2);
        // (2,3,4) about Z by 90°: (x,y)→(-y,x) → (-3,2,4). Bit-exact.
        assert!(
            pts.iter().any(|p| *p == [-3.0, 2.0, 4.0]),
            "corner lands exactly on the grid; got {pts:?}"
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// Re-rotating about the same axis chains a second forest node onto the first
    /// (this cell does not bundle): the leaf's parent is the earlier rotation, `base`
    /// stays the Constructed root, and the solid remains a rigid (volume/area-invariant)
    /// `Rotated` solid that validate/tess/STEP accept and a boolean now runs against.
    #[test]
    fn rerotate_same_axis_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap();
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::Z, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::Z])),
            "two Z nodes chained, base = the Constructed root"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap();
        assert!(
            (after.volume - before.volume).abs() < 1e-9,
            "volume invariant"
        );
        assert!((after.area - before.area).abs() < 1e-9, "area invariant");
        assert!(nacre_tess::to_obj(&m).is_ok());
        assert!(
            nacre_step::to_step(&m)
                .unwrap()
                .contains("MANIFOLD_SOLID_BREP")
        );
        assert!(solid_is_rotated(&m, c2), "re-rotated solid stays Rotated");
        // A boolean against the chain-rotated solid runs (guard retired, 3d-i): the axis-aligned
        // `d` inside the re-rotated `c2` is carved out, and the result is a valid solid.
        let d = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "chain-rotated cut is valid"
        );
    }

    /// Re-rotating about a different axis chains a node whose parent is the first
    /// rotation (root → Z → X), `base` still the root.
    #[test]
    fn rerotate_different_axis_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap().volume;
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::X])),
            "chain root→Z→X, base = root"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c2).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");
    }

    /// An **exact** (90°-family) rotation applied to an already-rotated solid is still
    /// recorded as a chain node — the composite is inexact (an ancestor is), so the
    /// forest must stay complete (1b silently dropped it). The solid stays Rotated and
    /// boolean-rejected.
    #[test]
    fn rerotate_exact_after_inexact_records_node() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 90)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c2),
            Some((false, 2, vec![Axis::Z, Axis::X])),
            "the exact 90°X is recorded as a chain node, not dropped"
        );
        assert!(
            solid_is_rotated(&m, c2),
            "composite is inexact → still Rotated"
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A translation between two same-axis rotations forces a chain (this cell never
    /// bundles anyway): the forest records both rotations, `base` stays the root, and
    /// the result is rigid and valid — sound with no adjacency guard (each rotation is
    /// its own node).
    #[test]
    fn rerotate_across_translation_chains() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let before = nacre_props::mass_props(&m, c).unwrap().volume;
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &translate_iso([5, -3, 2])).unwrap();
        m.rebuild_adjacency();
        let c3 = transform(&mut m, c2, &rot_iso(Axis::Z, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c3),
            Some((false, 2, vec![Axis::Z, Axis::Z])),
            "both rotations recorded; the intervening translation is not in the forest"
        );
        assert!(nacre_validate::validate(&m).is_empty());
        let after = nacre_props::mass_props(&m, c3).unwrap().volume;
        assert!((after - before).abs() < 1e-9, "volume invariant");
    }

    /// A three-axis chain records three nodes root→Z→X→Y.
    #[test]
    fn rerotate_deep_chain() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
        m.rebuild_adjacency();
        let c3 = transform(&mut m, c2, &rot_iso(Axis::Y, 20)).unwrap();
        m.rebuild_adjacency();

        assert_eq!(
            forest_probe(&m, c3),
            Some((false, 3, vec![Axis::Z, Axis::X, Axis::Y])),
        );
        assert!(nacre_validate::validate(&m).is_empty());
    }

    /// A fresh rotation of a Constructed solid is unchanged from 1b (B0): an inexact
    /// angle records a single root node; a 90°-family angle stays Constructed (no node,
    /// boolean allowed). Guards that the B0/B1 split preserves fresh-rotation behavior.
    #[test]
    fn fresh_rotation_of_constructed_unchanged() {
        use nacre_scalar::Axis;
        // inexact → one root node, base = the cuboid's Constructed vertices.
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        assert_eq!(forest_probe(&m, c1), Some((false, 1, vec![Axis::Z])));

        // exact 90° → Constructed (no node), boolean allowed.
        let mut m = Model::new();
        let c = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
        m.rebuild_adjacency();
        assert!(!solid_is_rotated(&m, c1), "fresh 90° stays Constructed");
        assert!(nacre_validate::validate(&m).is_empty());
        assert_eq!(forest_probe(&m, c1), None, "no rotation node");
    }

    /// Replay determinism (DNA 3): a re-rotation sequence reproduces the same forest
    /// and handles.
    #[test]
    fn rerotate_is_deterministic() {
        use nacre_scalar::Axis;
        let build = || {
            let mut m = Model::new();
            let c = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
            let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
            (bbox_lo(&m, c2), c2)
        };
        assert_eq!(build(), build(), "same ops → same geometry and handle");
    }


    #[test]
    fn cut_containment_makes_a_cavity() {
        // A = [0,3]³ (27) with B = [1,2]³ (1) strictly inside ⇒ A − B is a
        // hollow solid: volume 26, an outer + one void shell (V16/E24/F12/S2).
        let (mut m, a, b) = nested_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 26.0).abs() < 1e-12, "volume {vol}");
        assert_eq!(m.solids.get(r).cavities.len(), 1);
        let reach = m.reachable();
        assert_eq!(reach.shells.len(), 2);
        assert_eq!(reach.faces.len(), 12);
        assert_eq!(reach.vertices.len(), 16);
        assert_eq!(reach.edges.len(), 24);
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn cut_containment_off_center_cavity() {
        // The inner box need not be concentric — any strictly-interior B works.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 3.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - (64.0 - 6.0)).abs() < 1e-12, "volume {vol}"); // 4³ − 1·2·3
        assert_eq!(m.solids.get(r).cavities.len(), 1);
    }

    #[test]
    fn fuse_containment_is_the_container() {
        // A ∪ B with B ⊂ A is just A (no cavity).
        let (mut m, a, b) = nested_boxes();
        let vol_a = nacre_props::mass_props(&m, a).unwrap().volume;
        let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_a).abs() < 1e-12, "volume {vol}");
        assert!(m.solids.get(r).cavities.is_empty());
        assert_eq!(m.live_solids, vec![r]);
    }

    #[test]
    fn containment_symmetric_when_a_inside_b() {
        // Arguments swapped: A = inner ⊂ B = outer.
        // Cut(inner − outer): inner is wholly removed ⇒ empty, which is an answer, not an error —
        // and a successful boolean consumes its operands, so the Fuse below needs a fresh model
        // (it used to reuse this one only because the empty Cut was an error that consumed nothing).
        let (mut m, outer, inner) = nested_boxes();
        assert!(
            boolean(&mut m, BoolKind::Cut, inner, outer)
                .unwrap()
                .is_empty()
        );
        assert!(m.live_solids.is_empty(), "both operands are consumed");

        // Fuse(inner ∪ outer) = outer.
        let (mut m, outer, inner) = nested_boxes();
        let vol_outer = nacre_props::mass_props(&m, outer).unwrap().volume;
        let r = boolean_one(&mut m, BoolKind::Fuse, inner, outer).unwrap();
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - vol_outer).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn common_containment_is_the_inner_solid() {
        // Either argument order ⇒ the intersection is the inner solid (handled
        // by the existing half-space enumeration path, no cavity code).
        for swap in [false, true] {
            let (mut m, outer, inner) = nested_boxes();
            let vol_inner = nacre_props::mass_props(&m, inner).unwrap().volume;
            let (x, y) = if swap { (inner, outer) } else { (outer, inner) };
            let r = boolean_one(&mut m, BoolKind::Common, x, y).unwrap();
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            assert!((vol - vol_inner).abs() < 1e-9, "swap={swap} volume {vol}");
        }
    }


    // ---- coincident-coplanar merge (M5-c5) ----

    fn stacked_cubes() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        (m, a, b)
    }

    #[test]
    fn fuse_a_boss_onto_a_non_convex_solid() {
        // A contained boss on the top of an L-prism (non-convex kept `a`). The convexity gate
        // used to decline this to the seam path, which rejected the seamless contact; the
        // contained-coplanar Fuse now admits it. Volume 3 (L) + 0.4²·0.5 = 3.08.
        let (mut m, l) = l_prism(); // L footprint area 3, height 1
        let boss = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 1.0]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 3.08).abs() < 1e-12, "volume {vol}");
        // The boss top cap sits on the z = 1.5 plane, its outward normal +z.
        assert!(has_face_on_plane(
            &m,
            r,
            Point3::from_array([0.5, 0.5, 1.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        ));
    }

    #[test]
    fn fuse_a_non_convex_profile_boss() {
        // An L-shaped boss (non-convex cutter `b`) on a cube top. The gate used to decline the
        // non-convex prism; the contained-coplanar Fuse now carries the L footprint as a hole.
        // Volume 1 (cube) + 0.12 (L area) · 0.4 = 1.048.
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let l_base: Vec<Point3> = [
            [0.3, 0.3],
            [0.7, 0.3],
            [0.7, 0.5],
            [0.5, 0.5],
            [0.5, 0.7],
            [0.3, 0.7],
        ]
        .iter()
        .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
        .collect();
        let (boss, _) =
            build_prism(&mut m, &l_base, Vector3::from_array([0.0, 0.0, 0.4]), None).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, cube, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.048).abs() < 1e-12, "volume {vol}");
        // The L boss top cap sits on the z = 1.4 plane, its outward normal +z.
        assert!(has_face_on_plane(
            &m,
            r,
            Point3::from_array([0.4, 0.4, 1.4]),
            Vector3::from_array([0.0, 0.0, 1.0]),
        ));
    }

    #[test]
    fn an_edge_slot_through_the_bottom() {
        // The prism pokes out the base's bottom too, so the old convex/blind overhang-cut gate
        // declined it and the seam path could not build it either (honest reject). The F2 dispatch
        // collapse hands it to the unified coplanar driver, which carves the slot exactly:
        // base 1.0 − (x∈[0.5,1] · y∈[0.25,0.75] · z∈[0,1]) = 1 − 0.25 = 0.75.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.25, -0.5]),
            Point3::from_array([1.5, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_corner_cut_through_the_bottom() {
        // The corner prism pokes out the base's bottom, so the old blind gate declined it. The F2
        // collapse routes it to the unified driver: base 1.0 − corner column (x,y ∈ [0.5,1], full
        // height) = 1 − 0.25 = 0.75.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.5, -0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.75).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_boss_that_pierces_the_base_is_not_an_overhang() {
        // The boss dips below the base's top (its walls cross the base) — a transversal seam cut,
        // not a coplanar overhang, which the arrangement handles as a normal crossing.
        // Union = 1.0 + boss 0.75 − overlap 0.125 = 1.625.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 0.5]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.625).abs() < 1e-12, "volume {vol}");
    }

    // The Cut and Common twins of `fuse_a_corner_overhanging_boss` below — the same two solids,
    // the same shared z=1 plane. The boss sits entirely above it, so it removes nothing and shares
    // nothing: Cut is the base untouched and Common is empty. Both used to be rejected
    // (`coplanar_merge` for Cut; Common's every face dropped, which assembly reported as
    // `no_outward_shell`). The `Whole` survival cell now checks whether the contact plane actually
    // separates the solids, which is what makes the whole cap correct here.
    #[test]
    fn cut_by_a_corner_overhanging_boss_removes_nothing() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
        // Structure, not just volume: the base comes through as itself. Six faces means the cap was
        // not split along ∂Q and the boss contributed nothing.
        let s = m.solids.get(r);
        assert!(s.cavities.is_empty());
        assert_eq!(m.shells.get(s.outer).faces.len(), 6, "a clean cube");
    }

    #[test]
    fn cut_a_seated_block_by_the_part_below_it() {
        // The operands swapped: now the canonical contact face is the upper block's *lower* cap, so
        // the separation test runs with the plane's normal the other way round. Same answer — the
        // block keeps its volume.
        let mut m = Model::new();
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        m.rebuild_adjacency();
        let r = boolean_one(&mut m, BoolKind::Cut, block, base).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.0).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn common_with_a_corner_overhanging_boss_is_empty() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        m.rebuild_adjacency();
        // They meet only along the base's top face — a contact of zero volume.
        assert!(
            boolean(&mut m, BoolKind::Common, base, corner)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn cut_by_an_overhanging_boss_carrying_a_pin_owes_a_notch() {
        // Same seating, but the tool carries a pin reaching below the contact plane, so the plane
        // no longer separates the solids and the cut owes a real notch (1 − 0.2·0.2·0.5 = 0.98).
        //
        // This used to be an honest reject: the tool's z=1 cap is an annulus-like face whose
        // *inner* edge rides the pin's walls, and reading its occupancy off the ring's flank put
        // the material on the wrong side, so the class would not label. With the side read from
        // the ring's travel instead (`arrangement::run_body_above`), the notch comes out at the
        // hand-computed volume with a clean model.
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let block = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let pin = m.add_cuboid(
            Point3::from_array([0.55, 0.55, 0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        m.rebuild_adjacency();
        let tool = boolean_one(&mut m, BoolKind::Fuse, block, pin).unwrap();
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, tool).unwrap().volume - 1.02).abs() < 1e-12,
            "the pinned tool itself"
        );
        let notched =
            boolean_one(&mut m, BoolKind::Cut, base, tool).expect("the notch is buildable");
        m.rebuild_adjacency();
        assert!(
            (nacre_props::mass_props(&m, notched).unwrap().volume - 0.98).abs() < 1e-12,
            "the notch the tool owes: {}",
            nacre_props::mass_props(&m, notched).unwrap().volume
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "a notched result is still a clean model"
        );
    }

    // Build a unit cube with a blind pocket in its top — a non-convex solid whose side faces stay
    // convex. Returns `(model, pocketed solid)`.
    fn top_pocketed_cube() -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &extrude_op(square(), 1.0)).unwrap()
        else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid, .. } =
            apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap()
        else {
            unreachable!()
        };
        (m, solid)
    }

    // A1: plane-class canonicalization — coplanar walls of the two operands fold into one line.
    #[test]
    fn plane_classes_merge_a_shared_wall() {
        // Two unit cubes side by side share the plane x=1 (a's +x wall, b's -x wall — the same
        // plane, opposite normals). `plane_classes` must merge those two into one line class and
        // keep the far walls (a's x=0, b's x=2) distinct.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([1.0, 0.0, 0.0]),
            Point3::from_array([2.0, 1.0, 1.0]),
        );
        m.rebuild_adjacency();
        let planes_a = collect_planes(&m, a).unwrap();
        let na = planes_a.len();
        let mut planes = planes_a;
        planes.extend(collect_planes(&m, b).unwrap());
        // Find a plane by outward-normal x-sign and its x coordinate, within an index range.
        let find = |rng: std::ops::Range<usize>, nx: f64, x: f64| -> usize {
            rng.clone()
                .find(|&i| {
                    let n = planes[i].n_out.as_array();
                    n[0] * nx > 0.5 && (planes[i].tri[0].as_array()[0] - x).abs() < 1e-9
                })
                .expect("plane")
        };
        let a_xp = find(0..na, 1.0, 1.0); // a's +x wall at x=1
        let b_xm = find(na..planes.len(), -1.0, 1.0); // b's -x wall at x=1
        let a_xm = find(0..na, -1.0, 0.0); // a's -x wall at x=0
        let b_xp = find(na..planes.len(), 1.0, 2.0); // b's +x wall at x=2
        let canon = plane_classes(&planes);
        assert_eq!(canon[a_xp], canon[b_xm], "shared x=1 wall is one class");
        assert_ne!(canon[a_xm], canon[b_xp], "far walls stay distinct");
        assert_ne!(canon[a_xp], canon[a_xm], "x=1 and x=0 are different lines");
        // Every b face but its +x wall is coplanar with an a face, so 12 planes fold to 7 classes.
        let distinct: std::collections::HashSet<usize> = canon.iter().copied().collect();
        assert_eq!(distinct.len(), na + 1, "only b's far wall is a new class");
        // The class root is the smallest index in the class (deterministic canon).
        assert_eq!(canon[a_xp], a_xp.min(b_xm));
    }

    #[test]
    fn fuse_an_overhanging_boss_onto_a_non_convex_solid() {
        // A boss cantilevers off the +x side face of a top-pocketed cube (non-convex solid),
        // overhanging the bottom edge. The whole-solid gate used to block it; the contact face
        // (+x side) is a convex square, so the footprint gate admits it and the Fuse reconstruction
        // is local (the far pocket is verbatim-copied). Volume: pocketed 0.92 + boss 0.25 = 1.17.
        let (mut m, pc) = top_pocketed_cube();
        let boss = m.add_cuboid(
            Point3::from_array([1.0, 0.25, -0.25]),
            Point3::from_array([1.5, 0.75, 0.75]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, pc, boss).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        assert!((nacre_props::mass_props(&m, r).unwrap().volume - 1.17).abs() < 1e-12);
    }

    #[test]
    fn overhang_boss_with_a_non_convex_footprint() {
        // An L-shaped (non-convex) boss footprint overhanging a cube edge. The old convexity-gated
        // overhang detector declined this, so it used to be an honest reject; the F2 dispatch
        // collapse hands it to the unified coplanar driver, which builds it exactly. Volume =
        // cube 1.0 + L-prism (area 0.9·0.2 + 0.3·0.2 = 0.24) · height 0.4 = 1.096.
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let l_base: Vec<Point3> = [
            [0.3, 0.3],
            [1.2, 0.3],
            [1.2, 0.5],
            [0.6, 0.5],
            [0.6, 0.7],
            [0.3, 0.7],
        ]
        .iter()
        .map(|&[x, y]| Point3::from_array([x, y, 1.0]))
        .collect();
        let (l_tool, _) =
            build_prism(&mut m, &l_base, Vector3::from_array([0.0, 0.0, 0.4]), None).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, cube, l_tool).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.096).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn cut_a_blind_pocket_into_a_non_convex_solid() {
        // A blind pocket carved into an already-pocketed (non-convex) cube: a second contained
        // top-flush prism in a corner away from the first pocket. The kept solid `a` is non-convex,
        // which the pocket contact now admits (the gates are convexity-agnostic). Removed
        // 0.2·0.1·0.4 = 0.008 on top of the first pocket's 0.08 → 1 − 0.08 − 0.008 = 0.912.
        let (mut m, pc) = pocketed_cube();
        let corner = m.add_cuboid(
            Point3::from_array([0.05, 0.1, 0.6]),
            Point3::from_array([0.25, 0.2, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, pc, corner).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 0.912).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn cut_a_non_convex_blind_pocket() {
        // A blind pocket with a non-convex (L-shaped) footprint: the cutter prism is non-convex,
        // which the pocket contact now admits. The L extrudes to z∈[0,0.5], top-flush on the base's
        // z=0.5 face, blind. L area = 0.6² − 0.3² = 0.27, depth 0.5 → removed 0.135; base 3²·1.5 =
        // 13.5 → 13.365.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 0.5]),
        );
        let l = Profile2d {
            points: vec![
                p2(-0.3, -0.3),
                p2(0.3, -0.3),
                p2(0.3, 0.0),
                p2(0.0, 0.0),
                p2(0.0, 0.3),
                p2(-0.3, 0.3),
            ],
        };
        let OpOutput::Extrude { solid: lp, .. } = apply(&mut m, &extrude_op(l, 0.5)).unwrap()
        else {
            unreachable!()
        };
        let r = boolean_one(&mut m, BoolKind::Cut, base, lp).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 13.365).abs() < 1e-12, "volume {vol}");
    }

    #[test]
    fn a_pocket_that_punches_through_drills_a_bore() {
        // A top-flush tool that pokes out the base's bottom: the pocket becomes a through hole.
        // The exit face has to come out annular, and until the coplanar reconstruct learned to
        // emit a hole it came out whole instead, leaving the bore's walls nothing to close
        // against — an open shell the assembly guard rejected. (Honest reject, never a wrong
        // answer; the previous cell pinned it as such.)
        //
        // Area is the assertion that matters here: volume alone cannot tell a bore from a shape
        // that merely displaces the same material. 0.75 (top) + 0.75 (bottom) + 4 (sides) +
        // 2.0 (the bore's four inner walls) = 7.5, against 6.0 for the cube.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 0.75).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.5).abs() < 1e-12, "area {}", p.area);
        // A bore, not a void: no cavity shell, and both caps carry the hole (the top from the
        // coincident contact, the bottom from the section the tool cuts through it).
        let s = m.solids.get(r);
        assert!(s.cavities.is_empty(), "a through hole is not a cavity");
        let faces = &m.shells.get(s.outer).faces;
        assert_eq!(faces.len(), 10, "6 base faces + the bore's 4 walls");
        assert_eq!(
            faces
                .iter()
                .filter(|&&fh| !m.faces.get(fh).inner.is_empty())
                .count(),
            2,
            "both caps are annular"
        );
    }

    #[test]
    fn a_boss_that_punches_through_keeps_the_stub() {
        // The Fuse twin, and the same emission path: the tool's a-side face keeps `P∖Q`, so the
        // base's bottom needs the same hole for the stub below it to join on. Volume
        // 1 + 0.5·0.5·0.5 = 1.125; area 1.0 (top, the flush tool cap dissolves into it) + 0.75
        // (bottom) + 4 (sides) + 1.0 (stub walls) + 0.25 (stub floor) = 7.0.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let through = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, through).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 1.125).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
        let s = m.solids.get(r);
        assert_eq!(
            m.shells
                .get(s.outer)
                .faces
                .iter()
                .filter(|&&fh| !m.faces.get(fh).inner.is_empty())
                .count(),
            1,
            "only the bottom is annular — the flush top merges away"
        );
    }

    #[test]
    fn a_fused_stack_chains_through_a_cut() {
        // The dissolved 1×1×2 box (cell fuse-coplanar-merge) feeds a second boolean. Before
        // the merge/dissolve this rejected — first as COPLANAR_PAIR (the flat edges), then as
        // LOOP_ORIENT_MISMATCH (the straight-angle interface corners). A clean box cuts.
        let (mut m, a, b) = stacked_cubes();
        let stack = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        // A cutter straddling z=1 (the fused interface) — the seam runs where the split
        // vertical edges used to be. Result: 2 − 0.5·0.5·1.0.
        let cutter = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Cut, stack, cutter).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!((vol - 1.75).abs() < 1e-12, "volume {vol}");
    }

    // ---- `unify_coplanar_faces`: the general (interface-free) coplanar merge, on hand-built
    // `LocalFace` lists the coincident goldens above never reach — chains, opposite normals,
    // holes, seam edges, and the asymmetric T-junction the global dissolve exists to prevent. ----

    /// An axis-aligned `FaceInfo` at `d` along its normal, with a **non-degenerate `tri`** whose
    /// right-hand normal is `n_out`. The merge reads more than `n_out` now — `loop_winding` and
    /// `point_in_ring` name their arguments by plane and evaluate exact predicates on `tri` — so a
    /// dummy triangle would make those answers meaningless.
    fn mk_axis_plane(m: &mut Model, axis: usize, d: f64, positive: bool) -> PlaneGeom {
        let mut n = [0.0; 3];
        n[axis] = if positive { 1.0 } else { -1.0 };
        let normal = Vector3::from_array(n);
        let mut at = [0.0; 3];
        at[axis] = d;
        let origin = Point3::from_array(at);
        let plane = Plane::from_point_normal(origin, normal).unwrap();
        let surf = m.surfaces.push(Surface::Plane(plane));
        let face = m.faces.push(Face {
            surface: surf,
            outer: Loop { half_edges: vec![] },
            inner: vec![],
            orientation: Orientation::Forward,
        });
        // Two in-plane directions whose cross product is `+normal`, so `tri` winds outward.
        let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
        let (i, j) = if positive { (i, j) } else { (j, i) };
        let step = |k: usize| {
            let mut q = at;
            q[k] += 1.0;
            Point3::from_array(q)
        };
        let _ = face;
        PlaneGeom {
            surf,
            plane,
            tri: [origin, step(i), step(j)],
            tri_pt3: None,
            frame_sign: 1, // `plane` is built from `normal`, so the two agree
        }
    }

    /// A `FaceInfo` for the `unify_coplanar_faces` tests, which read none of its geometry; the rest is a
    /// valid-but-unreferenced dummy (`surf`/`face`/`plane` are never dereferenced there).
    fn face(plane_idx: usize, nodes: Vec<Node>, inner: Vec<Vec<Node>>) -> LocalFace {
        LocalFace {
            plane_idx,
            loop_nodes: nodes,
            inner,
            flip: false,
        }
    }

    #[test]
    fn unify_merges_a_coplanar_chain() {
        // Three unit squares on z=0 (+z), tiled in x, each sharing a vertical edge with the
        // next. One plane class, one normal ⇒ all fuse into a single face; the four
        // straight-angle mid-edge vertices dissolve, leaving one 4-corner rectangle.
        //
        // Named the way the arrangement names things — every vertex is the meeting of three
        // planes — because the straight-angle test reads those triples. (The old `Orig` fixture
        // exercised a path the engine stopped producing when it went all-`Seam`.)
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0, the shared class
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
            mk_axis_plane(&mut m, 0, 3.0, true),  // 4: x=3
            mk_axis_plane(&mut m, 1, 0.0, false), // 5: y=0
            mk_axis_plane(&mut m, 1, 1.0, true),  // 6: y=1
        ];
        let _canon: Vec<usize> = (0..p.len()).collect();
        let v = |x: usize, y: usize| Node::Seam([0, x, y]); // sorted: class, x-plane, y-plane
        let (c00, c10, c20, c30) = (v(1, 5), v(2, 5), v(3, 5), v(4, 5));
        let (c01, c11, c21, c31) = (v(1, 6), v(2, 6), v(3, 6), v(4, 6));
        let faces = vec![
            face(0, vec![c00, c10, c11, c01], vec![]),
            face(0, vec![c10, c20, c21, c11], vec![]),
            face(0, vec![c20, c30, c31, c21], vec![]),
        ];
        let out = unify_coplanar_faces(faces, &p).unwrap();
        assert_eq!(out.len(), 1, "three coplanar faces fuse into one");
        let l = &out[0].loop_nodes;
        assert_eq!(l.len(), 4, "straight-angle mid vertices dissolved: {l:?}");
        for c in [c00, c30, c31, c01] {
            assert!(l.contains(&c), "corner kept");
        }
        for c in [c10, c20, c11, c21] {
            assert!(!l.contains(&c), "mid vertex dropped");
        }
    }

    #[test]
    fn an_overhang_fuse_keeps_the_two_z1_caps_separate() {
        // An overhanging boss splits `z = 1` between two coplanar faces with **opposite** outward
        // normals — the base's exposed top (`+z`) and the boss underside (`-z`). They must not be
        // fused into one face: their `flip` differs, so `unify`'s `(plane_idx, flip)` group key
        // keeps them apart.
        //
        // This replaces two retired tests (`unify_keeps_opposite_normal_coplanar`,
        // `unify_keeps_holed_faces`) that built `Node::Orig` faces to exercise a passthrough the
        // arrangement never triggers — it emits all-`Seam`. The real invariant is exercised here on
        // the production path, in the default `cargo test` run: `overhang_fuse_then_cut_matches_occt`
        // proves it against OCCT but is `#[ignore]`, so this hand-computed volume is the non-ignored
        // guard. A wrong merge collapses the topology — the volume shifts or `validate` speaks.
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let r = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        assert!(
            (vol - 1.5).abs() < 1e-12,
            "base 1 + boss 0.5, no overlap: {vol}"
        );
    }

    #[test]
    fn a_hole_filled_by_two_faces_still_merges() {
        // Replaces `unify_skips_seam_shared_edges`, whose premise is gone twice over: the `Orig`
        // gate it pinned was removed when the engine went all-`Seam`, and `splice_along`, whose
        // panic it guarded against, no longer exists.
        //
        // What matters now is that the merge is not special-cased to "a hole filled by exactly one
        // neighbour". A [0,3]² face with a [1,2]² hole, and that hole filled by **two** pieces split
        // at x=1.5: every ring edge between them is carried in both directions, so erasing interior
        // boundary leaves only the outer square — one face, no hole, whatever the filling is cut
        // into. This is the case that separates a general rule from a bespoke one.
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 1.5, true),  // 3: x=1.5, where the filling is split
            mk_axis_plane(&mut m, 0, 2.0, true),  // 4: x=2
            mk_axis_plane(&mut m, 0, 3.0, true),  // 5: x=3
            mk_axis_plane(&mut m, 1, 0.0, false), // 6: y=0
            mk_axis_plane(&mut m, 1, 1.0, true),  // 7: y=1
            mk_axis_plane(&mut m, 1, 2.0, true),  // 8: y=2
            mk_axis_plane(&mut m, 1, 3.0, true),  // 9: y=3
        ];
        let _canon: Vec<usize> = (0..p.len()).collect();
        let v = |x: usize, y: usize| Node::Seam([0, x, y]);
        let (o00, o30, o33, o03) = (v(1, 6), v(5, 6), v(5, 9), v(1, 9));
        let (h11, h12, h22, h21) = (v(2, 7), v(2, 8), v(4, 8), v(4, 7));
        let (m12, m11) = (v(3, 8), v(3, 7)); // the split points on the hole's top and bottom
        let faces = vec![
            // Outer square with the hole, wound the way `emit_faces` states it: outer CCW, hole CW.
            face(
                0,
                vec![o00, o30, o33, o03],
                vec![vec![h11, h12, m12, h22, h21, m11]],
            ),
            face(0, vec![h11, m11, m12, h12], vec![]), // left filler
            face(0, vec![m11, h21, h22, m12], vec![]), // right filler
        ];
        let out = unify_coplanar_faces(faces, &p).unwrap();
        assert_eq!(out.len(), 1, "the hole is filled, so one face remains");
        assert!(out[0].inner.is_empty(), "and it has no hole left");
        assert_eq!(out[0].loop_nodes.len(), 4, "just the outer square");
        for c in [o00, o30, o33, o03] {
            assert!(out[0].loop_nodes.contains(&c), "outer corner kept");
        }
    }

    #[test]
    fn unify_keeps_a_vertex_that_is_a_corner_elsewhere() {
        // F0,F1 on z=0 merge; their shared-edge endpoint (1,0,0) is a straight angle on the
        // merged face but a real corner on a perpendicular face G (plane y=0). Global degree 3
        // ⇒ it is NOT dissolved — a per-face local rule would have, opening a T-junction. Its
        // twin (1,1,0), on the merged face only, IS dissolved.
        let mut m = Model::new();
        let p = vec![
            mk_axis_plane(&mut m, 2, 0.0, true),  // 0: z=0
            mk_axis_plane(&mut m, 0, 0.0, false), // 1: x=0
            mk_axis_plane(&mut m, 0, 1.0, true),  // 2: x=1
            mk_axis_plane(&mut m, 0, 2.0, true),  // 3: x=2
            mk_axis_plane(&mut m, 1, 0.0, false), // 4: y=0, the perpendicular face's plane
            mk_axis_plane(&mut m, 1, 1.0, true),  // 5: y=1
            mk_axis_plane(&mut m, 2, 1.0, true),  // 6: z=1
        ];
        let _canon: Vec<usize> = (0..p.len()).collect();
        let (v000, v100, v200) = (
            Node::Seam([0, 1, 4]),
            Node::Seam([0, 2, 4]),
            Node::Seam([0, 3, 4]),
        );
        let (v010, v110, v210) = (
            Node::Seam([0, 1, 5]),
            Node::Seam([0, 2, 5]),
            Node::Seam([0, 3, 5]),
        );
        let (v101, v201) = (Node::Seam([2, 4, 6]), Node::Seam([3, 4, 6]));
        let faces = vec![
            face(0, vec![v000, v100, v110, v010], vec![]),
            face(0, vec![v100, v200, v210, v110], vec![]),
            face(4, vec![v200, v100, v101, v201], vec![]), // perpendicular, not coplanar
        ];
        let out = unify_coplanar_faces(faces, &p).unwrap();
        assert_eq!(out.len(), 2, "z=0 pair merges; G stays");
        let merged = out.iter().find(|lf| lf.plane_idx == 0).unwrap();
        assert!(
            merged.loop_nodes.contains(&v100),
            "corner-elsewhere vertex kept (no T-junction)"
        );
        assert!(
            !merged.loop_nodes.contains(&v110),
            "pure straight-angle vertex dropped"
        );
    }



    proptest! {
        /// Diagonal corner overlaps (clean seam): fuse/cut volumes match the
        /// independent AABB formula (not nacre's own common).
        #[test]
        fn fuse_cut_diagonal_boxes_match_aabb(
            amin in prop::array::uniform3(-3.0f64..3.0),
            aext in prop::array::uniform3(1.0f64..3.0),
            t in prop::array::uniform3(0.15f64..0.6),
            s in prop::array::uniform3(0.3f64..2.0),
        ) {
            let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
            let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
            let bmax: [f64; 3] = std::array::from_fn(|i| amax[i] + s[i]);
            let ov: f64 = (0..3).map(|i| amax[i] - bmin[i]).product();
            let va: f64 = aext.iter().product();
            let vb: f64 = (0..3).map(|i| bmax[i] - bmin[i]).product();

            let build = || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
                let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
                (m, a, b)
            };

            let (mut m1, a1, b1) = build();
            let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1);
            prop_assume!(rf.is_ok()); // skip rare coplanar/degenerate configs
            let rf = rf.unwrap();
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb - ov)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            let rc = boolean_one(&mut m2, BoolKind::Cut, a2, b2).unwrap();
            m2.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m2).is_empty());
            let vc = nacre_props::mass_props(&m2, rc).unwrap().volume;
            prop_assert!((vc - (va - ov)).abs() <= 1e-9 * va, "cut {vc}");
        }

        /// Matched-footprint stacked boxes (coincident z-interface): fuse volume
        /// is the sum, cut is A, common is empty.
        #[test]
        fn stacked_boxes_merge_volumes(
            x0 in -3.0f64..3.0,
            y0 in -3.0f64..3.0,
            dx in 0.5f64..3.0,
            dy in 0.5f64..3.0,
            z0 in -3.0f64..3.0,
            h1 in 0.5f64..3.0,
            h2 in 0.5f64..3.0,
        ) {
            let (x1, y1) = (x0 + dx, y0 + dy);
            let (zm, z1) = (z0 + h1, z0 + h1 + h2);
            let build = || {
                let mut m = Model::new();
                let a = m.add_cuboid(Point3::from_array([x0, y0, z0]), Point3::from_array([x1, y1, zm]));
                let b = m.add_cuboid(Point3::from_array([x0, y0, zm]), Point3::from_array([x1, y1, z1]));
                (m, a, b)
            };
            let (va, vb) = (dx * dy * h1, dx * dy * h2);

            let (mut m1, a1, b1) = build();
            // Not `prop_assume!`: a matched-footprint stack is squarely in coverage whatever the
            // dimensions are, so a reject here is a defect, not an uninteresting sample. Assuming
            // it away is how this property went on passing while the kernel aborted on 2% of the
            // space and rejected 95% of it (family #3).
            let rf = boolean_one(&mut m1, BoolKind::Fuse, a1, b1)
                .expect("stacked boxes fuse at any dimensions");
            m1.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m1).is_empty());
            let vf = nacre_props::mass_props(&m1, rf).unwrap().volume;
            prop_assert!((vf - (va + vb)).abs() <= 1e-9 * (va + vb), "fuse {vf}");

            let (mut m2, a2, b2) = build();
            // The stack shares only its interface plane ⇒ no volume in common, at any dimensions.
            prop_assert!(boolean(&mut m2, BoolKind::Common, a2, b2).unwrap().is_empty());

            let (mut m3, a3, b3) = build();
            let rc = boolean_one(&mut m3, BoolKind::Cut, a3, b3).unwrap();
            let vc = nacre_props::mass_props(&m3, rc).unwrap().volume;
            prop_assert!((vc - va).abs() <= 1e-9 * va, "cut {vc}");
        }
    }

    /// **The topology of a boolean does not depend on whether the coordinates are
    /// f64-representable.** The same shape is built twice — once on tidy integers, once on
    /// dimensions that are not exact binary fractions — and both must give the same b-rep counts,
    /// with each volume matching its own formula.
    ///
    /// This is the invariant family #3 restored. Plane identity used to be read from the faces'
    /// *derived* coefficients, which are not exactly proportional for two differently-sized faces
    /// on one plane, so one plane became two classes and the arrangement named one point twice —
    /// but only when the arithmetic did not happen to cancel, which tidy coordinates hid
    /// (measured before the fix: 200/200 random stacked pairs under-merged, 12/600 ops aborted).
    #[test]
    fn boolean_topology_is_the_same_on_untidy_coordinates() {
        let counts = |dx: f64, dy: f64, z0: f64, h1: f64, h2: f64, kind: BoolKind| {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0, 0.0, z0]),
                Point3::from_array([dx, dy, z0 + h1]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.0, 0.0, z0 + h1]),
                Point3::from_array([dx, dy, z0 + h1 + h2]),
            );
            let r = boolean_one(&mut m, kind, a, b).expect("stacked boxes fuse/cut");
            m.rebuild_adjacency();
            assert!(nacre_validate::validate(&m).is_empty());
            let s = m.solids.get(r);
            let sh = m.shells.get(s.outer);
            let faces = sh.faces.len();
            let mut edges = std::collections::HashSet::new();
            let mut verts = std::collections::HashSet::new();
            for &fh in &sh.faces {
                let f = m.faces.get(fh);
                for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for he in &l.half_edges {
                        edges.insert(he.edge);
                        if let Some(bd) = m.edges.get(he.edge).bounds {
                            verts.extend(bd);
                        }
                    }
                }
            }
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            ((faces, edges.len(), verts.len(), s.cavities.len()), vol)
        };
        // The untidy dimensions are the minimal case the proptest shrank to when the kernel aborted.
        let (tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2) = (2.0, 0.5, 0.0, 0.5, 2.0);
        let (dx, dy, z0, h1, h2) = (
            1.628165457453874,
            0.5,
            0.11200046228159026,
            0.5,
            2.07926124157585,
        );
        for kind in [BoolKind::Fuse, BoolKind::Cut] {
            let (tidy_shape, tidy_vol) = counts(tidy_dx, tidy_dy, tidy_z0, tidy_h1, tidy_h2, kind);
            let (shape, vol) = counts(dx, dy, z0, h1, h2, kind);
            assert_eq!(tidy_shape, shape, "{kind:?}: same topology either way");
            let want = |a: f64, b: f64| match kind {
                BoolKind::Fuse => a + b,
                _ => a,
            };
            let (tw, w) = (
                want(tidy_dx * tidy_dy * tidy_h1, tidy_dx * tidy_dy * tidy_h2),
                want(dx * dy * h1, dx * dy * h2),
            );
            assert!((tidy_vol - tw).abs() < 1e-9, "{kind:?} tidy {tidy_vol}");
            assert!((vol - w).abs() < 1e-9 * w, "{kind:?} untidy {vol}");
        }
    }

    /// One geometric plane is one class **whatever the two faces' sizes**, and the coefficient test
    /// alone still cannot say so — the two walls' un-normalized coefficient 4-vectors are not
    /// exactly proportional. Pins that the coordinate branch is what earns the merge.
    #[test]
    fn one_plane_is_one_class_whatever_the_face_size() {
        let (dx, dy) = (1.628165457453874f64, 0.5f64);
        let (z0, h1, h2) = (0.11200046228159026f64, 0.5f64, 2.07926124157585f64);
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0]),
            Point3::from_array([dx, dy, z0 + h1]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, z0 + h1]),
            Point3::from_array([dx, dy, z0 + h1 + h2]),
        );
        let PlaneSetup {
            planes: faces_tab,
            geom: _planes,
            plane_ix,
            ..
        } = plane_index_setup(&m, a, b).unwrap();
        // The two `+X` walls: same plane x = dx, different face sizes (heights h1 vs h2). Two
        // *faces*, so this searches the face table — the plane table holds one entry for both,
        // which is the property under test.
        let x_walls: Vec<usize> = (0..faces_tab.len())
            .filter(|&i| {
                faces_tab[i].n_out.as_array() == [1.0, 0.0, 0.0]
                    && (faces_tab[i].tri[0].as_array()[0] - dx).abs() < 1e-12
            })
            .collect();
        assert_eq!(x_walls.len(), 2, "one wall from each box: {x_walls:?}");
        let (i, j) = (x_walls[0], x_walls[1]);
        assert!(
            !planes_coplanar(&faces_tab[i].plane, &faces_tab[j].plane),
            "the coefficient test still cannot prove these coplanar — that is the whole point"
        );
        assert!(
            tolerant::t_planes_coplanar(&faces_tab, i, j),
            "coordinates can"
        );
        assert_eq!(plane_ix[i], plane_ix[j], "so they are one plane-table row");
    }

    /// A vertex where one plane is split between two faces is named by **the planes that touch it**,
    /// not by the loop's two neighbours.
    ///
    /// The overhang chain puts the base's exposed top and the cantilever's underside on one plane
    /// (`z = 1`, opposite normals — `unify` rightly keeps them apart, `canon` rightly calls them one
    /// class). At `(1, 0.25, 1)` the side wall's loop runs straight through their shared line, so the
    /// old rule named that vertex with the same plane twice: a triple defining no point, which the
    /// exact predicates — whose precondition is `D ≠ 0` — aborted on. Measured 2026-07-22 as the only
    /// path a degenerate triple reached them (16 arrivals in the OCCT suite, now 0).
    #[test]
    fn a_vertex_is_named_by_the_planes_that_touch_it() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let cutter = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            ..
        } = plane_index_setup(&m, overhung, cutter).unwrap();
        let _ = &faces_tab;
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, overhung) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &planes, &plane_ix)
                        .unwrap();
                for t in &tris {
                    // The triple is already dense plane ids: distinct means three real planes.
                    assert!(
                        t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                        "vertex triple {t:?} names one plane twice"
                    );
                    // A name that denotes three distinct classes must denote a real point.
                    assert!(
                        three_planes(
                            &planes[t[0]].plane,
                            &planes[t[1]].plane,
                            &planes[t[2]].plane
                        )
                        .is_some(),
                        "triple {t:?} defines no point"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "the chained operand has vertices to name");
    }

    /// **Dense plane ids are order-isomorphic to the sparse roots.** `dense_planes` ranks the class
    /// roots, so any comparison, sort or lex-min over plane indices reads the same either way.
    ///
    /// This is a **migration gate, not a permanent invariant**: it exists so the claim is measured
    /// before the split rides on it, and it retires with `canon` — its subject, not its coverage,
    /// is what goes away.
    #[test]
    fn dense_plane_ids_are_monotone_in_canon() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // An overhanging boss splits `z = 1` between two faces, so classes really do merge and the
        // ranking really does compress — without that the map is the identity and proves nothing.
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let probe = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        // Rebuild the pieces `dense_planes` consumes, so this locks its contract without needing
        // `canon` to escape `plane_index_setup`. `plane_classes` is the same union-find the setup
        // runs; `dense_planes` the same ranking.
        let mut faces = collect_planes(&m, chained).unwrap();
        faces.extend(collect_planes(&m, probe).unwrap());
        let canon = plane_classes(&faces);
        let (geom, plane_ix) = dense_planes(&faces, &canon);
        assert!(
            canon.iter().enumerate().any(|(i, &c)| c != i),
            "fixture has no split plane — the invariant would be vacuous"
        );
        assert!(
            geom.len() < canon.len(),
            "the ranking must actually compress"
        );
        for i in 0..canon.len() {
            for j in 0..canon.len() {
                assert_eq!(
                    canon[i].cmp(&canon[j]),
                    plane_ix[i].cmp(&plane_ix[j]),
                    "faces {i}/{j}: canon {}/{} vs dense {}/{}",
                    canon[i],
                    canon[j],
                    plane_ix[i],
                    plane_ix[j]
                );
            }
        }
    }

    /// **A plane triple is always in class form.** `planes` is a per-face table, so the same
    /// `usize` could mean "face" or "plane"; producers settle it by emitting class roots, and a
    /// consumer's raw `==` then means "same plane". Four silent-wrong bugs on this branch came from
    /// the two meanings meeting in one comparison, so the invariant is asserted, not assumed.
    #[test]
    fn plane_triples_are_always_canon() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        // An *overhanging* boss splits `z = 1` between two faces with opposite normals (the base's
        // exposed top and the boss underside) — the shape that makes "face index" and "plane index"
        // differ at all. A boss sitting wholly inside the top merges into one holed face instead.
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let probe = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            ..
        } = plane_index_setup(&m, chained, probe).unwrap();
        // The fixture must actually merge two faces into one plane, or this proves nothing.
        assert!(
            planes.len() < faces_tab.len(),
            "fixture has no split plane — the invariant would be vacuous"
        );
        // A producer hands out dense plane ids (`loop_triples` maps face indices through
        // `plane_ix`), so "every element is a plane, not a face" is now the type, not a runtime
        // check. What remains testable is that the ids are in range and sorted-distinct.
        let mut checked = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let mut rings = vec![
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &planes, &plane_ix)
                        .unwrap(),
                ];
                rings.extend(
                    combinatorics::hole_rings(&m, fh, p, &inc_a, &planes, &plane_ix).unwrap(),
                );
                for t in rings.iter().flatten() {
                    for &k in t {
                        assert!(
                            k < planes.len(),
                            "triple {t:?} names {k}, out of the plane table"
                        );
                    }
                    assert!(
                        t[0] < t[1] && t[1] < t[2],
                        "triple {t:?} is not three distinct planes in sorted order"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "the chained operand has vertices to name");
    }

    /// **A point reads zero on each of its three defining planes.** `t_orient3d`'s on-plane
    /// shortcut is a raw `==` against the triple, so a vertex on the query plane must name it by the
    /// same id the query uses. The face/plane split makes that automatic — a plane has exactly one
    /// id now, so the old failure (a vertex named by face 6 of the `z = 1` class invisible to a
    /// query about face 1 of it) cannot be expressed. What is left to check is the identity itself.
    #[test]
    fn a_vertex_on_the_cut_plane_reads_zero_whichever_face_names_it() {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let chained = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        let probe = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        let PlaneSetup {
            planes: faces_tab,
            geom: planes,
            surf_ix,
            inc_a,
            plane_ix,
            ..
        } = plane_index_setup(&m, chained, probe).unwrap();
        assert!(
            planes.len() < faces_tab.len(),
            "fixture has no split plane — the sibling faces this used to distinguish"
        );
        let mut on_plane = 0usize;
        for sh in solid_shell_handles(&m, chained) {
            for &fh in &m.shells.get(sh).faces {
                let p = surf_ix[&fh];
                let tris =
                    combinatorics::face_vertex_triples(&m, fh, p, &inc_a, &planes, &plane_ix)
                        .unwrap();
                for t in &tris {
                    // The vertex lies on exactly its three defining planes; each must read 0.
                    for &q in t {
                        assert_eq!(
                            combinatorics::side_of(&planes, *t, q),
                            0,
                            "vertex {t:?} lies on plane {q} but does not read 0"
                        );
                        on_plane += 1;
                    }
                }
            }
        }
        assert!(on_plane > 0, "some vertex lies on some queried plane");
    }



    #[test]
    fn boolean_rejects_non_live_input() {
        let (mut m, a, b) = two_boxes();
        m.live_solids.retain(|&s| s != b); // as if superseded
        assert_eq!(
            boolean_one(&mut m, BoolKind::Common, a, b),
            Err(BoolError::InputNotLive)
        );
    }

    #[test]
    fn boolean_op_applies_and_wraps_error() {
        // A failing boolean's error is surfaced as `OpError::Boolean`. This used to be driven by a
        // disjoint `Common`, but that is no longer an error (it is an empty result, see
        // `boolean_op_passes_an_empty_result_through`), so the wrapping is exercised with a boolean
        // that genuinely fails: a handle that is not live.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        m.live_solids.retain(|&s| s != b); // retire `b` behind the op's back
        assert_eq!(
            apply(
                &mut m,
                &Operation::Boolean {
                    kind: BoolKind::Common,
                    a,
                    b
                }
            ),
            Err(OpError::Boolean(BoolError::InputNotLive))
        );
    }

    /// **An empty result is an answer.** Two solids that miss each other have no intersection, and
    /// that is what `Common` reports: `Ok` with no solids, both operands consumed like any other
    /// successful boolean. Stated on its own because the name is the contract — if someone makes
    /// this an error again, the failure points straight at what was decided (2026-07-22), and the
    /// `live_solids` assertion pins the retire that an early return would otherwise skip.
    #[test]
    fn a_disjoint_common_is_empty_not_an_error() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
        let solids = boolean(&mut m, BoolKind::Common, a, b).expect("empty is not a failure");
        assert!(solids.is_empty());
        assert!(
            m.live_solids.is_empty(),
            "a successful boolean consumes its operands"
        );
    }

    /// An empty boolean reaches the caller as an empty solid list, not an error — the op layer
    /// passes the kernel's answer through rather than reinterpreting it.
    #[test]
    fn boolean_op_passes_an_empty_result_through() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        // Offset in all axes so no faces are coplanar with A.
        let b = m.add_cuboid(Point3::from_array([10.0; 3]), Point3::from_array([11.0; 3]));
        let out = apply(
            &mut m,
            &Operation::Boolean {
                kind: BoolKind::Common,
                a,
                b,
            },
        )
        .unwrap();
        assert_eq!(out, OpOutput::Boolean { solids: vec![] });
        assert!(m.live_solids.is_empty(), "both operands are consumed");
    }

    // ---- boolean Common algorithm (M5-c3 commit 2) ----




    #[test]
    fn common_rejects_non_planar_input() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let cyl = m.add_cylinder(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        assert_rejects(
            || boolean_one(&mut m, BoolKind::Common, a, cyl),
            tag::CYLINDER_FACE,
        );
    }

    #[test]
    fn a_corner_flush_common_keeps_the_non_convex_overlap() {
        // A **corner-flush** `Common`: the L-prism and the box both start at the origin, so
        // **three** of their face planes coincide — `z = 0` (both floors), `x = 0`, `y = 0`. Every
        // vertex of the shared corner lies exactly on the other solid's face planes, which is what
        // the old reject tag said: `VERTEX_ON_FACE_PLANE`. The arrangement engine names such a
        // point by its plane triple like any other, so the configuration is no longer special.
        //
        // Not covered by the other two non-convex `Common` locks:
        // `common_non_convex_overlap_is_their_intersection` (l_and_corner_box) and
        // `common_non_convex_containment_is_inner` both meet transversally, with no coplanar pair.
        //
        // Hand-checked shape, not just volume. The overlap is the L
        // `x∈[0,1.5]×y∈[0,1]` (1.5) plus `x∈[0,1]×y∈[1,1.5]` (0.5) = 2.0, over `z∈[0,0.5]`:
        //   volume 2.0 · 0.5 = 1.0
        //   area   2 · 2.0 (caps) + 6.0 (the L's perimeter) · 0.5 = 7.0
        // and the L has six sides, so eight faces.
        let l = Profile2d {
            points: vec![
                p2(0.0, 0.0),
                p2(2.0, 0.0),
                p2(2.0, 1.0),
                p2(1.0, 1.0),
                p2(1.0, 2.0),
                p2(0.0, 2.0),
            ],
        };
        let mut m = replay(&[extrude_op(l, 1.0)]).unwrap();
        let lsolid = *m.live_solids.first().unwrap();
        let b = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.5, 1.5, 0.5]),
        );
        let r = boolean_one(&mut m, BoolKind::Common, lsolid, b).expect("corner-flush Common");
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{vs:?}");
        let p = nacre_props::mass_props(&m, r).unwrap();
        assert!((p.volume - 1.0).abs() < 1e-12, "volume {}", p.volume);
        assert!((p.area - 7.0).abs() < 1e-12, "area {}", p.area);
        assert_eq!(m.solids.get(r).cavities.len(), 0);
    }

    proptest! {
        /// Overlapping axis-aligned boxes: the intersection volume equals the
        /// independent AABB-overlap product (mixed A/B axis-aligned vertices).
        #[test]
        fn common_axis_boxes_volume_matches_aabb_overlap(
            amin in prop::array::uniform3(-5.0f64..5.0),
            aext in prop::array::uniform3(1.0f64..4.0),
            t in prop::array::uniform3(0.05f64..0.7),
            bext in prop::array::uniform3(1.0f64..4.0),
        ) {
            let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
            let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
            let bmax: [f64; 3] = std::array::from_fn(|i| bmin[i] + bext[i]);
            let expected: f64 = (0..3)
                .map(|i| (amax[i].min(bmax[i]) - bmin[i]).max(0.0))
                .product();
            prop_assume!(expected > 1e-3);

            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
            let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
            let res = boolean_one(&mut m, BoolKind::Common, a, b);
            prop_assume!(res.is_ok()); // skip rare coplanar/degenerate configs
            let r = res.unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            let vol = nacre_props::mass_props(&m, r).unwrap().volume;
            prop_assert!((vol - expected).abs() <= 1e-9 * expected.max(1.0), "{vol} vs {expected}");
        }

        /// A tilted square prism (oblique planes) intersected with a big enclosing
        /// box is the prism — exercises non-axis-aligned face normals (the in/out
        /// sign and CCW ordering) with an independent oracle (the prism's own mass).
        #[test]
        fn common_tilted_prism_with_enclosing_box_is_the_prism(
            nx in -0.5f64..0.5,
            ny in -0.5f64..0.5,
        ) {
            let plane = SketchPlane::from_origin_normal(
                Point3::origin(),
                Vector3::from_array([nx, ny, 1.0]),
            )
            .unwrap();
            let mut m = replay(&[Operation::Extrude {
                plane,
                profile: square(),
                dist: 1.0,
            }])
            .unwrap();
            let prism = *m.live_solids.first().unwrap();
            let vol_prism = nacre_props::mass_props(&m, prism).unwrap().volume;
            let c = m.add_cuboid(Point3::from_array([-10.0; 3]), Point3::from_array([10.0; 3]));
            let res = boolean_one(&mut m, BoolKind::Common, prism, c);
            prop_assume!(res.is_ok());
            let r = res.unwrap();
            m.rebuild_adjacency();
            prop_assert!(nacre_validate::validate(&m).is_empty());
            let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
            prop_assert!(
                (vol_r - vol_prism).abs() <= 1e-9 * vol_prism.max(1.0),
                "{vol_r} vs {vol_prism}"
            );
        }
    }
}
