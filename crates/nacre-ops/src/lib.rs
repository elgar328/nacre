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
mod rotated_vertex;
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
    /// A surviving cavity that no material component contains — geometrically impossible for a
    /// valid boolean result (a void lies inside exactly one piece). A defensive backstop; cavity
    /// ownership is otherwise decided exactly by [`combinatorics::point_in_component`] containment.
    pub const CAVITY_NO_OWNER: &str = "cavity_no_owner";
    /// No material-enclosing (outward) shell among the result components — every component is
    /// inward-oriented. Geometrically impossible for a real solid result; a defensive backstop
    /// with no firing test (cf. `FOURPLANE`).
    pub const NO_OUTWARD_SHELL: &str = "no_outward_shell";
    /// A rotated-result face's supporting plane could not be witnessed exactly: neither three
    /// of the face's own vertices assemble (a survivor wall) nor a unique operand plane `π` is
    /// recoverable from its seam corners' provenance ([`crate::rotated_vertex::face_plane_witness`]).
    /// The deep rotated-chain floor (e.g. a seam-only face whose `π` is itself a rotated surface).
    /// Honest reject, never a wrong result.
    pub const ROTATED_UNDERDETERMINED: &str = "rotated_underdetermined";
    /// The assembled result has an **odd Euler characteristic** (`V − E + F − L_i`), which no
    /// closed 2-manifold can have (it must equal the even `2(S − G)`) — so the arrangement produced
    /// a malformed solid and the boolean rejects rather than return it (DNA: never silently wrong).
    /// This is the Euler-parity backstop for malformity that is *not* a pinch (see
    /// `NON_MANIFOLD_VERTEX`); e.g. a dropped face. Rotation-independent. Checked post-assembly in
    /// `boolean`, per solid.
    pub const EULER_PARITY: &str = "euler_parity";
    /// The assembled result has a **non-manifold vertex** — a "pinch" where two or more face-fans
    /// meet at one point (a cutter's convex corner exactly on the target's concave corner; two
    /// solids touching only at a corner), even though every edge is manifold. No valid 2-manifold
    /// solid has one, so the boolean rejects with this clear reason rather than the incidental
    /// `EULER_PARITY` (which also misses an *even* number of pinches). Rotation-independent; checked
    /// per solid post-assembly via `nacre_topo::nonmanifold_vertices`.
    pub const NON_MANIFOLD_VERTEX: &str = "non_manifold_vertex";
    /// The assembled result has an even Euler characteristic but a **negative genus** (`S − χ/2 < 0`)
    /// — more handles than a solid can have, so it is not a valid closed 2-manifold. A count-based
    /// backstop below the pinch and parity checks; checked per solid post-assembly.
    pub const NEGATIVE_GENUS: &str = "negative_genus";
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
    use nacre_cip::Pt3;
    use nacre_geom::intersect::{planes_coplanar, three_planes};
    use nacre_geom::{Plane, Surface};
    use nacre_topo::{Loop, Orientation, Origin, VertexDef};
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

    // --- Rotated booleans go live (overhaul 3d-i, `ROTATED_UNSUPPORTED` retired) ---
    // A boolean commutes with a rigid motion, so rotating both operands by the same
    // irrational-angle isometry must give the rigid image of the unrotated result — identical
    // volume, solid count, and cavity count, and still valid. These are the first live proof
    // that the CIP-wired machinery (arrangement, seam, in/out, outer/cavity — 3a–3c-vi) is
    // sound end-to-end on rotated (rounded-irrational) geometry.

    /// A cutter's convex corner landing exactly on the target's concave corner — three planes
    /// (x=1,y=1,z=1) meeting at one point (1,1,1), six faces there — makes a **non-manifold pinch**
    /// (two face-fans touching at the point). No valid 2-manifold solid exists, so `boolean` rejects
    /// with the clear `NON_MANIFOLD_VERTEX` reason (not the incidental `EULER_PARITY`), leaving the
    /// live model untouched. **Rotation-independent**: both the axis-aligned and the rotated framings
    /// (exact and CIP-kernel paths) hit the same pinch and reject. R = a cube minus a far-corner
    /// octant; C = a cube whose +corner is that removed octant's inner corner. (Coplanar contact away
    /// from a corner works — [`a_rotated_boolean_result_can_be_cut_again`].)
    #[test]
    fn a_corner_coincident_cut_is_rejected_not_silently_wrong() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        // `rotate`: None = axis-aligned (exact path); Some = the result and cutter tilted (CIP path).
        let run = |rotate: bool| {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
            let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
            // C's +corner is (1,1,1) = R's concave corner → three shared planes meet there.
            let c = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, -1.0]),
                Point3::from_array([1.0; 3]),
            );
            m.rebuild_adjacency();
            let (r, c) = if rotate {
                let r = transform(&mut m, r, &iso).unwrap();
                m.rebuild_adjacency();
                let c = transform(&mut m, c, &iso).unwrap();
                m.rebuild_adjacency();
                (r, c)
            } else {
                (r, c)
            };
            let live = m.live_solids.clone();
            assert_rejects(
                || boolean(&mut m, BoolKind::Cut, r, c),
                tag::NON_MANIFOLD_VERTEX,
            );
            assert_eq!(m.live_solids, live, "reject must not mutate the live set");
        };
        run(false); // axis-aligned
        run(true); // rotated
    }

    /// Two cubes touching only at the corner (1,1,1): their Fuse pinches two solids at a single
    /// vertex (non-manifold), so it is rejected with the clear `NON_MANIFOLD_VERTEX` — the direct,
    /// minimal pinch (every edge is manifold; only the vertex is the defect).
    #[test]
    fn two_cubes_touching_at_a_corner_fuse_is_non_manifold() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        m.rebuild_adjacency();
        let live = m.live_solids.clone();
        assert_rejects(
            || boolean(&mut m, BoolKind::Fuse, a, b),
            tag::NON_MANIFOLD_VERTEX,
        );
        assert_eq!(m.live_solids, live, "reject must not mutate the live set");
    }

    /// The kernel's validator now catches the pinch too (it previously only exposed it indirectly as
    /// `EulerParity`). Bypass `boolean`'s reject via `arrangement::boolean` to obtain the malformed
    /// corner-touch Fuse solid, then confirm `validate` reports a `NonManifoldVertex`.
    #[test]
    fn validate_reports_the_non_manifold_pinch() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        m.rebuild_adjacency();
        crate::arrangement::boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(
            issues
                .iter()
                .any(|v| matches!(v, nacre_validate::Violation::NonManifoldVertex { .. })),
            "validate must flag the pinch: {issues:?}"
        );
    }

    /// Predicates over a rotated result's witness planes are rotation-invariant against the same
    /// result unrotated — the provenance witness (with outward winding) defines the exact plane,
    /// so `t_orient3d` agrees on every definite triple. A regression guard on the witness itself.
    #[test]
    fn rotated_result_witness_predicates_are_invariant() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        m.rebuild_adjacency();
        let table = |m: &Model, s: Handle<Solid>| {
            let f = collect_planes(m, s).unwrap();
            let c = plane_classes(&f);
            dense_planes(&f, &c).0
        };
        let pu = table(&m, r);
        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let pr = table(&m, r2);
        assert_eq!(pu.len(), pr.len(), "rotation preserves the plane count");
        let n = pu.len();
        let indep = |p: &[PlaneGeom], a: usize, b: usize, c: usize| {
            let nrm = |k: usize| p[k].plane.normal();
            nrm(a).dot(nrm(b).cross(nrm(c))).abs() > 0.3
        };
        let mut disagree = 0;
        for p in 0..n {
            for q in (p + 1)..n {
                for rr in (q + 1)..n {
                    if !indep(&pu, p, q, rr) {
                        continue;
                    }
                    for j in 0..n {
                        if j == p || j == q || j == rr {
                            continue;
                        }
                        let su = crate::tolerant::t_orient3d(&pu, p, q, rr, j);
                        let sr = crate::tolerant::t_orient3d(&pr, p, q, rr, j);
                        if su != 0 && sr != 0 && su != sr {
                            eprintln!("DISAGREE orient3d ({p},{q},{rr},{j}): u={su} r={sr}");
                            disagree += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(
            disagree, 0,
            "{disagree} predicate disagreements (witness wrong)"
        );
    }

    /// Rotating a boolean *result* and feeding it back into a boolean (was `ROTATED_UNSUPPORTED`):
    /// `collect_planes` now witnesses each rotated seam face's plane through provenance — its
    /// plane is `R(π)` for an operand plane `π`, recovered from the operand face still on `π` and
    /// rotated by the face's own chain. A boolean commutes with a rigid motion, so the rotated
    /// chain's result matches the unrotated chain's (volume, solid count) and stays valid.
    #[test]
    fn a_rotated_boolean_result_can_be_cut_again() {
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let iso = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });
        // Chain: R = Cut(A, B) removes a far-corner octant; then Cut(R, C) removes a near one.
        let build = |m: &mut Model| -> (Handle<Solid>, Handle<Solid>) {
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
            let r = boolean_one(m, BoolKind::Cut, a, b).unwrap();
            // A clean slab that severs R at x = 0.5 — no plane of C coincides with any of R's
            // (avoids the separate rotated-coplanar-contact gap; isolates the witness).
            let c = m.add_cuboid(
                Point3::from_array([-1.0, -1.0, -1.0]),
                Point3::from_array([0.5, 4.0, 4.0]),
            );
            (r, c)
        };
        // Unrotated reference (reuse already works when nothing is rotated).
        let mut m0 = Model::new();
        let (r0, c0) = build(&mut m0);
        m0.rebuild_adjacency();
        let ref_out = boolean(&mut m0, BoolKind::Cut, r0, c0).unwrap();
        m0.rebuild_adjacency();
        let ref_vol: f64 = ref_out
            .iter()
            .map(|&s| nacre_props::mass_props(&m0, s).unwrap().volume)
            .sum();
        // Rotated: turn the *result* R (and C) by the same isometry, then reuse R.
        let mut m = Model::new();
        let (r, c) = build(&mut m);
        m.rebuild_adjacency();
        let r = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let c = transform(&mut m, c, &iso).unwrap();
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Cut, r, c).unwrap();
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(
            issues.is_empty(),
            "rotated-result reuse must be valid: {issues:?}"
        );
        let vol: f64 = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert_eq!(
            out.len(),
            ref_out.len(),
            "solid count invariant under rotation"
        );
        assert!(
            (vol - ref_vol).abs() < 1e-6,
            "rotated reuse volume {vol} vs unrotated {ref_vol}"
        );
    }

    /// Result-reuse rotation stress: build R with a first boolean, then feed R into a second
    /// boolean with a fresh cutter C — once unrotated, once with R and C rotated by the same
    /// isometry. A boolean commutes with a rigid motion, so the rotated reuse must equal the
    /// unrotated one (volume, solid count, cavity count) or be an honest reject — never silently
    /// wrong. This is the invariant on the newly-enabled rotated-*Discovered* geometry (every
    /// vertex of a boolean result is `Discovered`, so its rotation exercises the provenance
    /// witness on every face). `#[ignore]`: rotated booleans escalate to astro-float (~1–3 s).
    #[test]
    #[ignore = "slow: rotated result-reuse booleans (run with --ignored)"]
    fn rotated_result_reuse_stress() {
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
        // (first kind, second kind, |m| -> (a, b, c)). R = kind1(a, b); out = kind2(R, c).
        type Build = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>, Handle<Solid>)>;
        let cuboid = |m: &mut Model, lo: [f64; 3], hi: [f64; 3]| {
            m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
        };
        let fixtures: Vec<(&str, BoolKind, BoolKind, Build)> = vec![
            (
                "corner_then_slab",
                BoolKind::Cut,
                BoolKind::Cut,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0; 3]);
                    let b = cuboid(m, [1.0; 3], [3.0; 3]);
                    let c = cuboid(m, [-1.0, -1.0, -1.0], [0.5, 4.0, 4.0]);
                    (a, b, c)
                }),
            ),
            (
                "fuse_then_bite",
                BoolKind::Fuse,
                BoolKind::Cut,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0, 1.0, 1.0]);
                    let b = cuboid(m, [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]);
                    let c = cuboid(m, [1.3, 1.3, -1.0], [3.0, 3.0, 2.0]);
                    (a, b, c)
                }),
            ),
            (
                "cut_then_fuse",
                BoolKind::Cut,
                BoolKind::Fuse,
                Box::new(move |m: &mut Model| {
                    let a = cuboid(m, [0.0; 3], [2.0; 3]);
                    let b = cuboid(m, [1.3, 1.3, 1.3], [3.0, 3.0, 3.0]);
                    let c = cuboid(m, [-0.7, 0.4, 0.4], [0.3, 1.4, 1.4]);
                    (a, b, c)
                }),
            ),
        ];
        let isos_list: Vec<(&str, Vec<Isometry>)> = vec![
            ("Z43", vec![rot(Axis::Z, 43, [1, 1, 0])]),
            ("X67", vec![rot(Axis::X, 67, [2, -1, 0])]),
            (
                "Z50>Y37",
                vec![rot(Axis::Z, 50, [1, 1, 0]), rot(Axis::Y, 37, [0, 0, 1])],
            ),
        ];
        let run = |k1: BoolKind, k2: BoolKind, build: &Build, isos: &[Isometry]| -> Out {
            let mut m = Model::new();
            let (a, b, c) = build(&mut m);
            m.rebuild_adjacency();
            let Ok(r) = boolean(&mut m, k1, a, b) else {
                return Out::Rej;
            };
            assert_eq!(r.len(), 1, "first boolean is a single solid");
            let mut r = r[0];
            let mut c = c;
            m.rebuild_adjacency();
            for iso in isos {
                r = transform(&mut m, r, iso).unwrap();
                m.rebuild_adjacency();
                c = transform(&mut m, c, iso).unwrap();
                m.rebuild_adjacency();
            }
            match boolean(&mut m, k2, r, c) {
                Ok(solids) => {
                    m.rebuild_adjacency();
                    assert!(
                        nacre_validate::validate(&m).is_empty(),
                        "INVALID rotated reuse result"
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
        let (mut success, mut reject, mut silent) = (0, 0, 0);
        for (fname, k1, k2, build) in &fixtures {
            let base = run(*k1, *k2, build, &[]);
            for (rname, isos) in &isos_list {
                let r = run(*k1, *k2, build, isos);
                match (&base, &r) {
                    (_, Out::Rej) => reject += 1,
                    (Out::Ok(v1, s1, c1), Out::Ok(v2, s2, c2))
                        if vclose(*v1, *v2) && s1 == s2 && c1 == c2 =>
                    {
                        success += 1
                    }
                    _ => {
                        silent += 1;
                        eprintln!("SILENT-WRONG {fname} {rname} base={base:?} rot={r:?}");
                    }
                }
            }
        }
        eprintln!("REUSE STRESS: success={success} honest_reject={reject} SILENT_WRONG={silent}");
        assert_eq!(
            silent, 0,
            "a rotated result-reuse boolean was silently wrong"
        );
        assert!(
            success >= 1,
            "at least one rotated reuse must actually succeed"
        );
    }

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

    /// A sever that also leaves a surviving cavity: a hollow box whose void sits to one side,
    /// cut by a slab that severs it without touching the void. The x<2 piece keeps the void as a
    /// cavity, the x>2 piece is solid — two outward shells *and* one inward. `point_in_component`
    /// assigns the void to the x<2 piece that nests it (a containment test), rather than rejecting.
    #[test]
    fn severed_with_cavity_assigns_the_void_to_its_piece() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        // Void near the x-low side (1×2×2 = 4), clear of the x=2 cut.
        let inner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 2.5]),
        );
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // A slab spanning full y,z, thin in x at x∈[2,2.2] — severs into x<2 (holds the void, vol
        // 2·3·3 − 4 = 14) and x>2 (solid, vol 0.8·3·3 = 7.2).
        let slab = m.add_cuboid(
            Point3::from_array([2.0, -1.0, -1.0]),
            Point3::from_array([2.2, 4.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(
            solids.len(),
            2,
            "the slab severs the hollow box into two pieces"
        );
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // Exactly one piece owns the void; volumes match the hand calculation.
        let with_cav: Vec<_> = solids
            .iter()
            .filter(|&&s| !m.solids.get(s).cavities.is_empty())
            .collect();
        assert_eq!(
            with_cav.len(),
            1,
            "the void is assigned to exactly one piece"
        );
        let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
        let hollow_piece = *with_cav[0];
        assert!(
            (vol(hollow_piece) - 14.0).abs() < 1e-9,
            "hollow piece {}",
            vol(hollow_piece)
        );
        let total: f64 = solids.iter().map(|&s| vol(s)).sum();
        assert!((total - 21.2).abs() < 1e-9, "total {total}");
    }

    /// The adjacent case: a cut that passes *through* the void opens it — the void wall becomes
    /// exterior boundary, so no cavity survives. Handled by the plain sever path (no cavity to
    /// assign), not the containment code, but pinned so a regression there is caught.
    #[test]
    fn a_cut_through_the_void_leaves_no_cavity() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let inner = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3])); // void 1³
        let hollow = boolean_one(&mut m, BoolKind::Cut, big, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(m.solids.get(hollow).cavities.len(), 1);
        // Slab x∈[1.4,1.6] passes through the void (x∈[1,2]) → severs AND opens the void.
        let slab = m.add_cuboid(
            Point3::from_array([1.4, -1.0, -1.0]),
            Point3::from_array([1.6, 4.0, 4.0]),
        );
        let solids = boolean(&mut m, BoolKind::Cut, hollow, slab).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // The void is opened, so neither piece keeps a cavity; material = 26 − (1.8 − 0.2) = 24.4.
        let total_cavities: usize = solids.iter().map(|&s| m.solids.get(s).cavities.len()).sum();
        assert_eq!(
            total_cavities, 0,
            "the cut opened the void — no surviving cavity"
        );
        let total: f64 = solids
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).unwrap().volume)
            .sum();
        assert!((total - 24.4).abs() < 1e-9, "total {total}");
    }

    /// Nested cavities: a hollow box A ([0,6]³ − [1,5]³ void) with a smaller hollow box B
    /// ([2,4]³ − [2.5,3.5]³ void) floating inside A's void. `Fuse(A,B)` is one arrangement with
    /// four components (two materials, two voids); B's void is contained by **both** A's outer
    /// shell and B's own, so the containment assignment must pick the **innermost** (B), not A.
    /// The result is two solids, each keeping its own void (A: 216−64 = 152, B: 8−1 = 7).
    #[test]
    fn a_void_nested_in_a_floating_island_goes_to_the_inner_solid() {
        let mut m = Model::new();
        let big = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
        let void = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([5.0; 3]));
        let a = boolean_one(&mut m, BoolKind::Cut, big, void).unwrap();
        m.rebuild_adjacency();
        let bbig = m.add_cuboid(Point3::from_array([2.0; 3]), Point3::from_array([4.0; 3]));
        let bvoid = m.add_cuboid(Point3::from_array([2.5; 3]), Point3::from_array([3.5; 3]));
        let b = boolean_one(&mut m, BoolKind::Cut, bbig, bvoid).unwrap();
        m.rebuild_adjacency();
        let solids = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        assert_eq!(solids.len(), 2);
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        // Each solid keeps exactly one void — B's void was assigned to B (innermost), not A.
        let vol = |s| nacre_props::mass_props(&m, s).unwrap().volume;
        for &s in &solids {
            assert_eq!(
                m.solids.get(s).cavities.len(),
                1,
                "each piece keeps its own void"
            );
        }
        let mut vols: Vec<f64> = solids.iter().map(|&s| vol(s)).collect();
        vols.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((vols[0] - 7.0).abs() < 1e-9, "inner {}", vols[0]);
        assert!((vols[1] - 152.0).abs() < 1e-9, "outer {}", vols[1]);
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

    fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
        Operation::PocketOnFace {
            face,
            profile,
            dist,
        }
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
