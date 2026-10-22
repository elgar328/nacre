/// **Nothing carries plane coefficients through a frame, and that is the design, not an
/// omission.**
///
/// A plane that passes through a motion has two descriptions — coefficients and a witness triangle
/// — and mixing them can produce intransitive comparisons and a confident wrong winding. The frame
/// design sidesteps the whole question: a plane built in a frame
/// states itself **in that frame**, its name and its base triangle speak the *same* frame, and
/// every judgement downstream is a determinant question a rigid motion
/// preserves. So no `d` transformation term is needed anywhere.
///
/// The one place that does transport coefficients is
/// [`coplanar_by_composed_rotation`](super::coplanar_by_composed_rotation), whose whole method
/// is composing two chains into a single coordinate-axis rotation and applying it to a
/// coefficient vector. A frame is not that, so it must decline — a conservative miss that
/// sends the pair to the interval path. This asks for that decline directly, because the
/// alternative is `Isometry::plane_coeffs` being handed a frame it cannot express and
/// answering anyway.
#[test]
fn a_chain_holding_a_frame_never_reaches_the_coefficient_transport() {
    use super::{WitnessPoint, single_rotation};
    use nacre_exact::{Angle, Axis, Rat};
    let ri = |n: i128, d: i128| Rat::new(n, d).unwrap();
    // A tilted plane's frame: 2x - 3y + 7z + 11 = 0.
    let fr = nacre_exact::plane_frame([2, -3, 7, 11].map(Rat::from_int)).unwrap();
    let framed = |u: i128, v: i128, w: i128| {
        WitnessPoint::at([ri(u, 1), ri(v, 1), ri(w, 1)])
            .frame(fr)
            .unwrap()
    };
    let def = [framed(0, 0, 0), framed(1, 0, 0), framed(0, 1, 0)];
    assert!(
        single_rotation(&def).is_none(),
        "a frame is not a rotation about a coordinate axis; transporting coefficients \
             through it would answer with full confidence and be wrong"
    );
    // The same points with an ordinary rotation *are* accepted, so the decline above is the
    // frame's doing and not a shortcut that never fires.
    let turned = [
        WitnessPoint::at([ri(0, 1), ri(0, 1), ri(0, 1)]),
        WitnessPoint::at([ri(1, 1), ri(0, 1), ri(0, 1)]),
        WitnessPoint::at([ri(0, 1), ri(1, 1), ri(0, 1)]),
    ]
    .map(|p| {
        p.rotate_about(
            Axis::Z,
            Angle::from_deg(Rat::from_int(30)).unwrap(),
            [Rat::from_int(0); 3],
        )
    });
    assert!(single_rotation(&turned).is_some());
}

use super::*;
use nacre_exact::Mag;
use nacre_math::Point3;

/// How these fixtures judge; production chooses both per model. The coincidence limit is the
/// derived default for a unit-scale model — output resolution (`2⁻⁵²`) two words further down.
fn fixture() -> ProofStandard {
    ProofStandard {
        prec: 256,
        coincidence: crate::coincidence_for(Mag::of(1.0)),
        scale: Mag::of(1.0),
        cap: 4096,
    }
}

/// A synthetic axis-aligned plane witness: three points on the plane, its exact
/// coefficients, the exact `WitnessPoint` definition of those points, and the name those points
/// span, folded as `collect_planes` folds it. `is_rotated` is `false`, so the predicates answer by
/// the name — the definition is there but unused, which is precisely the arrangement the
/// production tables have.
struct W {
    tri: [Point3; 3],
    coeffs: [f64; 4],
    def: [WitnessPoint; 3],
    ints: Option<NameInts>,
}
impl W {
    fn new(tri: [Point3; 3], coeffs: [f64; 4]) -> W {
        let def = tri.map(|p| {
            WitnessPoint::at(
                p.as_array()
                    .map(|x| nacre_exact::Rat::try_from_f64(x).expect("test coordinate")),
            )
        });
        let mut w = W {
            tri,
            coeffs,
            def,
            ints: None,
        };
        let b = w.def.each_ref().map(|p| p.base);
        let name =
            nacre_exact::plane_name_exact(b[0], b[1], b[2]).expect("three points span a plane");
        w.ints = name_stored_ints(Some(&name), &w.def, w.frame_sign());
        assert!(w.ints.is_some(), "the points lie on their own name");
        w
    }
}
impl W {
    /// The coefficients and the three points, exactly.
    fn exact_inputs(&self) -> ([nacre_exact::Rat; 4], [[nacre_exact::Rat; 3]; 3]) {
        let r = |x: f64| nacre_exact::Rat::try_from_f64(x).expect("test value");
        (self.coeffs.map(r), self.tri.map(|p| p.as_array().map(r)))
    }

    /// The triangle spans a plane and the coefficients' normal is square to both its edges.
    fn normal_spans(&self) -> bool {
        let ([a, b, c, _], t) = self.exact_inputs();
        let edge = |i: usize| -> [nacre_exact::Rat; 3] {
            core::array::from_fn(|k| t[i][k].checked_sub(t[0][k]).expect("small"))
        };
        let (e1, e2) = (edge(1), edge(2));
        let n = [a, b, c];
        !nacre_exact::parallel_rat(&e1, &e2)
            && nacre_exact::dot_sign_rat(&n, &e1) == nacre_exact::Orient::Zero
            && nacre_exact::dot_sign_rat(&n, &e2) == nacre_exact::Orient::Zero
    }

    /// Every point satisfies `a·x + b·y + c·z + d = 0`.
    fn points_satisfy(&self) -> bool {
        let (c, t) = self.exact_inputs();
        t.iter().all(|p| {
            let mut acc = c[3];
            for k in 0..3 {
                acc = acc
                    .checked_add(c[k].checked_mul(p[k]).expect("small"))
                    .expect("small");
            }
            acc == nacre_exact::Rat::from_int(0)
        })
    }
}
impl Witness for W {
    // A hand-built witness has no name: every same-plane question takes the definitions' road.
    fn world_name(&self) -> Option<&nacre_exact::PlaneName> {
        None
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.def
    }
    fn is_rotated(&self) -> bool {
        false
    }
}
impl PlaneWitness for W {
    // These witnesses carry no motion, so there is nothing to cancel.
    fn chain_id(&self) -> u64 {
        0
    }
    fn name_ints(&self) -> Option<&NameInts> {
        self.ints.as_ref()
    }
    fn frame_sign(&self) -> i8 {
        let t = self.tri;
        let x = (t[1] - t[0]).cross(t[2] - t[0]).as_array();
        let c = self.coeffs;
        if c[0] * x[0] + c[1] * x[1] + c[2] * x[2] > 0.0 {
            1
        } else {
            -1
        }
    }
    // A synthetic witness has no name: its coefficients stand in for the plane, and they are
    // handed out only where its triangle lies on them, so a test's two descriptions are one plane
    // — what `exact_coeffs` promises (the plane itself). Asked of the witness's own inputs, lifted
    // exactly (an `f64` is a binary rational).
    fn exact_coeffs(&self) -> Option<[f64; 4]> {
        (self.normal_spans() && self.points_satisfy()).then_some(self.coeffs)
    }
    fn exact_normal(&self) -> Option<[f64; 3]> {
        self.normal_spans()
            .then(|| [self.coeffs[0], self.coeffs[1], self.coeffs[2]])
    }
    fn base_coeffs(&self) -> Option<[f64; 4]> {
        None
    }
}

/// The three coordinate planes through `(1,1,1)`: `x=1`, `y=1`, `z=1`, plus a `z=0` plane.
fn cube_corner_planes() -> Vec<W> {
    let p = |a, b, c| Point3::from_array([a, b, c]);
    vec![
        // x = 1  →  1·x + 0 + 0 − 1 = 0
        W::new(
            [p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(1.0, 0.0, 1.0)],
            [1.0, 0.0, 0.0, -1.0],
        ),
        // y = 1
        W::new(
            [p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 1.0)],
            [0.0, 1.0, 0.0, -1.0],
        ),
        // z = 1
        W::new(
            [p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, 1.0)],
            [0.0, 0.0, 1.0, -1.0],
        ),
        // z = 0
        W::new(
            [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)],
            [0.0, 0.0, 1.0, 0.0],
        ),
    ]
}

/// A plane witness carrying one shared rotation — with its pre-motion name ([`RW::named`]), what
/// [`cancel_cmp_coord`] needs. The motion identity, the names and the world statement default to
/// "one shared motion, no names"; the builders below change them for fixtures that need a chain of
/// their own or a name.
struct RW {
    tri: [Point3; 3],
    coeffs: [f64; 4],
    def: [WitnessPoint; 3],
    base_coeffs: [f64; 4],
    id: u64,
    world: Option<nacre_exact::PlaneName>,
    ints: Option<NameInts>,
}

/// `[a,b,c,d]` of the plane through three points (`n = e1 × e2`, `d = −n·p0`).
fn plane_of(t: [[f64; 3]; 3]) -> [f64; 4] {
    let e = |i: usize| [t[i][0] - t[0][0], t[i][1] - t[0][1], t[i][2] - t[0][2]];
    let (u, v) = (e(1), e(2));
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    [
        n[0],
        n[1],
        n[2],
        -(n[0] * t[0][0] + n[1] * t[0][1] + n[2] * t[0][2]),
    ]
}

impl RW {
    /// Three integer points on a plane, turned `deg`° about `axis` through a non-origin
    /// pivot — so the pivot translation is present and has to cancel in the difference.
    fn rotated(pts: [[i128; 3]; 3], axis: Axis, deg: i128) -> RW {
        RW::turned(pts, &[(axis, deg)])
    }

    /// The same, through a chain of turns — one node per `(axis, deg)`.
    fn turned(pts: [[i128; 3]; 3], nodes: &[(Axis, i128)]) -> RW {
        RW::moved(pts, |p| {
            nodes.iter().fold(p, |p, &(axis, deg)| {
                p.rotate_about(axis, Angle::from_deg(ri(deg, 1)).unwrap(), RW::pivot())
            })
        })
    }

    /// Off the origin, so a pivot's translation is present and has to cancel in a difference.
    fn pivot() -> [nacre_exact::Rat; 3] {
        [ri(1, 3), ri(1, 7), ri(0, 1)]
    }

    /// Three integer points on a plane, each carried by `motion`.
    fn moved(pts: [[i128; 3]; 3], motion: impl Fn(WitnessPoint) -> WitnessPoint) -> RW {
        let def: [WitnessPoint; 3] =
            pts.map(|q| motion(WitnessPoint::at([ri(q[0], 1), ri(q[1], 1), ri(q[2], 1)])));
        let base = pts.map(|q| Point3::from_array([q[0] as f64, q[1] as f64, q[2] as f64]));
        let tri: [Point3; 3] = std::array::from_fn(|i| Point3::from_array(def[i].coord_approx()));
        RW {
            coeffs: plane_of(tri.map(|p| p.as_array())),
            base_coeffs: plane_of(base.map(|p| p.as_array())),
            tri,
            def,
            id: 1,
            world: None,
            ints: None,
        }
    }

    /// …carrying a motion no other plane of the table shares (`0` for none).
    fn chain(mut self, id: u64) -> RW {
        self.id = id;
        self
    }

    /// The name of the plane the bases span — before any motion.
    fn base_name(&self) -> nacre_exact::PlaneName {
        let b = self.def.each_ref().map(|w| w.base);
        nacre_exact::plane_name_exact(b[0], b[1], b[2]).expect("three bases span a plane")
    }

    /// …with its pre-motion name's integers, folded as `collect_planes` folds them.
    fn named(mut self) -> RW {
        let name = self.base_name();
        self.ints = name_stored_ints(Some(&name), &self.def, self.frame_sign());
        assert!(self.ints.is_some(), "the bases lie on their own name");
        self
    }

    /// …stated in the world, which an unmoved plane's name is.
    fn stated(mut self) -> RW {
        assert!(
            self.def[0].chain.is_empty(),
            "only an unmoved plane names the world"
        );
        self.world = Some(self.base_name());
        self
    }
}
impl Witness for RW {
    // A hand-built witness has no name unless [`RW::stated`] gave it one: every same-plane
    // question takes the definitions' road.
    fn world_name(&self) -> Option<&nacre_exact::PlaneName> {
        self.world.as_ref()
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.def
    }
    fn is_rotated(&self) -> bool {
        true
    }
}
impl PlaneWitness for RW {
    fn chain_id(&self) -> u64 {
        self.id
    }
    fn name_ints(&self) -> Option<&NameInts> {
        self.ints.as_ref()
    }
    fn frame_sign(&self) -> i8 {
        let t = self.tri;
        let x = (t[1] - t[0]).cross(t[2] - t[0]).as_array();
        let c = self.coeffs;
        if c[0] * x[0] + c[1] * x[1] + c[2] * x[2] > 0.0 {
            1
        } else {
            -1
        }
    }
    // Rotated witnesses: no exact `f64` coefficients exist, as in the arrangement.
    fn exact_coeffs(&self) -> Option<[f64; 4]> {
        None
    }
    fn exact_normal(&self) -> Option<[f64; 3]> {
        None
    }
    fn base_coeffs(&self) -> Option<[f64; 4]> {
        Some(self.base_coeffs)
    }
}

/// A judging context over a fixture table.
///
/// The collector is leaked so a fixture stays a one-liner — a `Vec` per call, in a test binary,
/// and nothing reads it. A fixture that *does* want the evidence builds its own [`Notes`] and
/// calls [`Judge::new`].
fn jd<W>(planes: &[W]) -> Judge<'_, W> {
    Judge::new(planes, fixture(), Box::leak(Box::new(Notes::new())))
}

fn ri(n: i128, d: i128) -> nacre_exact::Rat {
    nacre_exact::Rat::new(n, d).unwrap()
}
use nacre_exact::{Angle, Axis};

/// **The pre-rotation shortcut must give the answer the escalation gives.**
///
/// [`cancel_cmp_coord`] answers in the pre-rotation frame on the planes' pre-motion names;
/// [`indirect_cmp_coord_judge`] answers by realizing the rotation in astro-float and reading
/// an interval. They share no bound and no code, so agreement is a real cross-check — and it
/// is the only one available, because the shortcut's whole point is to *not* do what the
/// judge does. This asks through the front door, where the shared-plane identity answers the
/// zeros along the rotation axis first; the shortcut itself is the next test's.
///
/// The fixtures turn about a **non-origin pivot**: a rotation about `c` is `x ↦ R(x−c)+c`, and the
/// `+c` cancels in a difference. If it did not, every answer here would be wrong.
#[test]
fn the_pre_rotation_shortcut_agrees_with_the_escalation() {
    for deg in [30i128, 37, 120, 200] {
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            // ∩(x=1, y=2, z=3) = (1,2,3) and ∩(x=1, y=2, z=7) = (1,2,7): the difference is
            // along z, so an in-plane comparison is exactly 0 and a z comparison is definite.
            // Then a pair differing along x, which the shortcut must decline rather than
            // guess (`tan Θ` is irrational, so it is provably nonzero but not signed here).
            let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg).named();
            let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg).named();
            let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg).named();
            let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg).named();
            let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg).named();
            let ps = vec![px, py, z3, z7, x5];
            for (a, b, what) in [
                ([0usize, 1, 2], [0usize, 1, 3], "differ along z"),
                ([0, 1, 2], [4, 1, 2], "differ along x"),
            ] {
                for k in 0..3 {
                    let got = jd(&ps).cmp_coord(a, b, k);
                    let want = to_i8(
                        indirect_cmp_coord_judge(
                            borrow_triple(a.map(|i| plane_def(&ps, i))),
                            borrow_triple(b.map(|i| plane_def(&ps, i))),
                            k,
                            fixture(),
                        )
                        .orient(),
                    );
                    assert_eq!(
                        got, want,
                        "{deg}° about {axis:?}, {what}, axis {k}: shortcut {got} vs \
                             escalation {want}"
                    );
                }
            }
        }
    }
}

/// …and the shortcut must actually fire, or the test above only proves the escalation agrees
/// with itself. The rotation axis is answered exactly, and so is a difference lying along it —
/// asked of the shortcut directly, so the identity in front of it cannot answer for it.
#[test]
fn the_shortcut_fires_where_it_should_and_declines_where_it_cannot() {
    let (axis, deg) = (Axis::Z, 30i128);
    let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg).named();
    let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg).named();
    let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg).named();
    let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg).named();
    let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg).named();
    let ps = vec![px, py, z3, z7, x5];
    let (a, b) = ([0usize, 1, 2], [0usize, 1, 3]); // differ along z only
    // Z is the rotation axis: preserved, so the sign comes back exactly.
    assert_eq!(cancel_cmp_coord(&ps, a, b, 2), Some(-1), "z: 3 < 7");
    // x and y: the difference lies along the rotation axis, so both are exactly 0.
    assert_eq!(cancel_cmp_coord(&ps, a, b, 0), Some(0));
    assert_eq!(cancel_cmp_coord(&ps, a, b, 1), Some(0));
    // A difference in the rotation plane: provably nonzero, but its sign needs the rotation
    // realized, so the shortcut declines instead of guessing.
    let c = [4usize, 1, 2];
    assert_eq!(cancel_cmp_coord(&ps, a, c, 0), None);
    assert_eq!(cancel_cmp_coord(&ps, a, c, 1), None);
    // …while the rotation axis still answers for that pair (both points share z = 3).
    assert_eq!(cancel_cmp_coord(&ps, a, c, 2), Some(0));
}

/// **A chain that turns about more than one axis must be declined, not answered.**
///
/// The whole derivation rests on the product of the rotations being a rotation *about a
/// coordinate axis* — that is what makes one coordinate fixed and the other two a plane
/// rotation with a single angle. Compose an X turn with a Z turn and none of that holds:
/// there is no preserved coordinate, and the in-plane formula is about the wrong plane. The
/// guard is the only thing standing between that and a confidently wrong sign, and without a
/// mixed-axis fixture nothing else in the suite notices if it is removed. The planes are named,
/// so the decline is the axes' doing and not a missing name's.
#[test]
fn a_chain_about_two_axes_is_declined() {
    let nodes: &[(Axis, i128)] = &[(Axis::X, 30), (Axis::Z, 40)];
    let t = |pts| RW::turned(pts, nodes).named();
    let ps = vec![
        t([[1, 0, 0], [1, 1, 0], [1, 0, 1]]),
        t([[0, 2, 0], [1, 2, 0], [0, 2, 1]]),
        t([[0, 0, 3], [1, 0, 3], [0, 1, 3]]),
        t([[0, 0, 7], [1, 0, 7], [0, 1, 7]]),
    ];
    let (a, b) = ([0usize, 1, 2], [0usize, 1, 3]);
    for k in 0..3 {
        assert_eq!(
            cancel_cmp_coord(&ps, a, b, k),
            None,
            "axis {k}: a two-axis chain has no preserved coordinate, so nothing here is \
                 decidable in the pre-rotation frame"
        );
        // …and the toleranced path still answers it, so declining costs only speed.
        let want = to_i8(
            indirect_cmp_coord_judge(
                borrow_triple(a.map(|i| plane_def(&ps, i))),
                borrow_triple(b.map(|i| plane_def(&ps, i))),
                k,
                fixture(),
            )
            .orient(),
        );
        assert_eq!(jd(&ps).cmp_coord(a, b, k), want);
    }
}

/// The corner `∩(x=1, y=1, z=1) = (1,1,1)` sits above the `z=0` plane, so its orient3d
/// against `z=0` is definite (nonzero), and querying `z=1` (a defining plane) is exactly 0.
#[test]
fn t_orient3d_axis_definite_and_on_plane() {
    let ps = cube_corner_planes();
    // query plane j = 3 (z=0): definite.
    assert_ne!(jd(&ps).orient3d(0, 1, 2, 3), 0);
    // query plane j = 2 (z=1) is one of the defining planes → exactly 0.
    assert_eq!(jd(&ps).orient3d(0, 1, 2, 2), 0);
}

/// A plane is coplanar with itself; two distinct planes are not.
#[test]
fn t_planes_coplanar_reflexive_and_distinct() {
    let ps = cube_corner_planes();
    assert!(
        jd(&ps).planes_coplanar(0, 0),
        "a plane is coplanar with itself"
    );
    assert!(
        !jd(&ps).planes_coplanar(0, 1),
        "x=1 and y=1 are distinct planes"
    );
}

/// `any_rotated` is false for axis-aligned witnesses (none of them `is_rotated`).
#[test]
fn any_rotated_false_for_axis_aligned() {
    let ps = cube_corner_planes();
    assert!(!any_rotated(&ps, &[0, 1, 2, 3]));
}

/// `frame_sign` recomputes the stored-vs-outward sign from `coeffs` + `tri` alone: `+1`
/// when the coefficient-normal agrees with `cross(tri)`, `-1` when the winding is reversed.
#[test]
fn frame_sign_from_coeffs_and_tri() {
    let ps = cube_corner_planes();
    // x=1: tri wound so cross(tri) = +x, and the coeffs normal is +x → +1.
    assert_eq!(ps[0].frame_sign(), 1);
    // reversing the tri winding flips cross(tri) → -1 (coeffs unchanged).
    let flipped = W::new([ps[0].tri[0], ps[0].tri[2], ps[0].tri[1]], ps[0].coeffs);
    assert_eq!(flipped.frame_sign(), -1);
}

/// **A name's exact `f64` row stands where every integer it holds *is* an `f64` within 127
/// bits; its nearest spelling stands wherever the integers are in the `f64` range.** The exact
/// row is what an unmoved plane's exact shortcuts read, so a row that stood for an integer no
/// `f64` holds would hand the predicates a value that is not the name. A power of two past 53
/// bits *is* an `f64` and stands; the 127-bit cap keeps the expansions' degree-6 products in
/// range, so a wide name has no exact row even where its integers are `f64`s. The nearest
/// spelling is what the name-integer road's filters read, and it fails only past the `f64` range.
#[test]
fn a_names_exact_row_stands_only_where_every_integer_is_an_f64() {
    use super::{WitnessPoint, name_stored_ints};
    use num_bigint::BigInt;
    let r = |n: i128, d: i128| Rat::new(n, d).expect("rational");
    let ints = |pts: [[Rat; 3]; 3]| {
        let name = nacre_exact::plane_name_exact(pts[0], pts[1], pts[2]).expect("a plane");
        name_stored_ints(Some(&name), &pts.map(WitnessPoint::at), 1).expect("folded")
    };
    // `z = 1/10`: the name `[0, 0, 10, −1]` is small, so the row is the name exactly — the cap
    // whose `f64` image (`fl(0.1)`) is not on it.
    let tenth = ints([
        [r(0, 1), r(0, 1), r(1, 10)],
        [r(1, 1), r(0, 1), r(1, 10)],
        [r(0, 1), r(1, 1), r(1, 10)],
    ]);
    let row = tenth.row.expect("a small name has a row");
    assert_eq!(row.map(f64::abs), [0.0, 0.0, 10.0, 1.0]);
    assert_eq!(
        tenth.normal.map(|n| n.map(f64::abs)),
        Some([0.0, 0.0, 10.0])
    );
    assert_eq!(tenth.near, tenth.row);
    // `a·x − z = 0`: `a = 2⁵³ + 1` is no `f64` — no exact row, a nearest one; `a = 2⁶⁰` is.
    let slope = |a: i128| {
        ints([
            [r(0, 1), r(0, 1), r(0, 1)],
            [r(0, 1), r(1, 1), r(0, 1)],
            [r(1, 1), r(0, 1), r(a, 1)],
        ])
    };
    let odd = slope((1i128 << 53) + 1);
    assert_eq!((odd.row, odd.normal), (None, None), "2⁵³ + 1 is no f64");
    assert!(odd.near.is_some() && odd.near_normal.is_some());
    let pow = slope(1i128 << 60);
    assert!(pow.row.is_some() && pow.normal.is_some(), "2⁶⁰ is an f64");
    // `z = 2⁶⁰ + 1`: only `d` is no `f64`, so the normal stands and the row does not.
    let far = (1i128 << 60) + 1;
    let high = ints([
        [r(0, 1), r(0, 1), r(far, 1)],
        [r(1, 1), r(0, 1), r(far, 1)],
        [r(0, 1), r(1, 1), r(far, 1)],
    ]);
    assert_eq!(high.row, None);
    assert_eq!(high.normal.map(|n| n.map(f64::abs)), Some([0.0, 0.0, 1.0]));
    // Wide names: the plane `x = 1` spelled with a scale of `2¹³⁰` or `2¹¹⁰⁰` (the predicates are
    // row-linear, so a scale is still that plane's name). `2¹³⁰` is an `f64` but past the cap, so
    // it has a nearest spelling and no exact row; `2¹¹⁰⁰` is past the `f64` range and has neither.
    let scaled = |e: u32| {
        let k = BigInt::from(1) << e;
        let name = nacre_exact::PlaneName::Wide([k.clone(), BigInt::from(0), BigInt::from(0), -k]);
        let on_x1 = [
            [r(1, 1), r(0, 1), r(0, 1)],
            [r(1, 1), r(1, 1), r(0, 1)],
            [r(1, 1), r(0, 1), r(1, 1)],
        ];
        name_stored_ints(Some(&name), &on_x1.map(WitnessPoint::at), 1).expect("folded")
    };
    let wide = scaled(130);
    assert_eq!((wide.row, wide.normal), (None, None), "past the cap");
    assert!(
        wide.near.is_some() && wide.near_normal.is_some(),
        "inside the f64 range"
    );
    let huge = scaled(1100);
    assert_eq!((huge.near, huge.near_normal, huge.row), (None, None, None));
}

/// A judge too poor to climb: every question the exact routes and the identities cannot answer
/// comes back `Exhausted`, so an empty collector after a question means one of those answered it.
fn starved() -> ProofStandard {
    ProofStandard {
        prec: 8,
        coincidence: Mag::pow2(-4000),
        scale: Mag::of(16.0),
        cap: 8,
    }
}

/// **The table the identity locks read.** Walls turned about z by different angles, each on a
/// chain of its own, beside two caps stated in the world — so no question here has one motion
/// under all its planes, and the shared-motion shortcuts (`cancel_cmp_coord`, the pre-motion
/// determinant) never answer first. A table with one shared chain would leave every lock below
/// green with the identity switched off.
fn identity_table() -> Vec<RW> {
    identity_planes(true)
}

/// The same planes, the turned ones with or without their pre-motion names — without them the
/// identity knows none of their zeros.
fn identity_planes(named: bool) -> Vec<RW> {
    let name = |w: RW| if named { w.named() } else { w };
    let z = |pts, deg, id| name(RW::rotated(pts, Axis::Z, deg).chain(id));
    vec![
        z([[1, 0, 0], [1, 1, 0], [1, 0, 1]], 30, 1), // 0: x = 1 turned 30° about z
        z([[0, 2, 0], [1, 2, 0], [0, 2, 1]], 37, 2), // 1: y = 2 turned 37°
        RW::turned([[0, 0, 3], [1, 0, 3], [0, 1, 3]], &[])
            .chain(0)
            .stated(), // 2: z = 3
        RW::turned([[0, 0, 7], [1, 0, 7], [0, 1, 7]], &[])
            .chain(0)
            .stated(), // 3: z = 7
        z([[0, 5, 0], [1, 5, 0], [0, 5, 1]], 50, 3), // 4: y = 5 turned 50°
        name(RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], Axis::X, 30).chain(4)), // 5: y = 2 turned 30° about x — normal (0, c, s)
        z([[3, 0, 0], [3, 1, 0], [3, 0, 1]], 20, 5), // 6: x = 3 turned 20°
    ]
}

/// The escalation's answer — what the identity must agree with, on no shared code.
fn escalated(ps: &[RW], a: [usize; 3], b: [usize; 3], axis: usize) -> i8 {
    to_i8(
        indirect_cmp_coord_judge(
            borrow_triple(a.map(|i| plane_def(ps, i))),
            borrow_triple(b.map(|i| plane_def(ps, i))),
            axis,
            fixture(),
        )
        .orient(),
    )
}

/// **Where the shared planes hold a coordinate, it is equal with no arithmetic** — one shared
/// plane whose normal is the axis, two whose line runs across it, three that are one point. The
/// starved judge answers each `0` and records nothing, which only an answer before the
/// escalation can do; the escalation at full budget agrees.
#[test]
fn the_planes_two_points_share_answer_a_held_coordinate_outright() {
    let ps = identity_table();
    let held: [([usize; 3], [usize; 3], usize, &str); 5] = [
        ([0, 1, 2], [0, 1, 3], 0, "two walls' vertical line, along x"),
        ([0, 1, 2], [0, 1, 3], 1, "two walls' vertical line, along y"),
        ([0, 1, 2], [4, 6, 2], 2, "one cap z = 3, along z"),
        ([0, 1, 2], [0, 1, 2], 0, "one point"),
        (
            [5, 0, 2],
            [5, 6, 2],
            2,
            "one cap again, beside a plane turned about x",
        ),
    ];
    for (a, b, axis, what) in held {
        let notes = Notes::new();
        let jd = Judge::new(&ps, starved(), &notes);
        assert_eq!(jd.cmp_coord(a, b, axis), 0, "{what}");
        assert!(
            notes.sorted().is_empty(),
            "{what}: the starved judge climbed"
        );
        assert_eq!(
            escalated(&ps, a, b, axis),
            0,
            "{what}: the escalation disagrees"
        );
    }
    // The control: the same true zero, with the walls' zeros unknown, is left `Exhausted` — a
    // nonzero answer would not do, since the filter decides those at any budget.
    let unnamed = identity_planes(false);
    let notes = Notes::new();
    let jd = Judge::new(&unnamed, starved(), &notes);
    assert_eq!(jd.cmp_coord([0, 1, 2], [0, 1, 3], 0), 0);
    assert!(
        notes
            .sorted()
            .iter()
            .any(|e| matches!(e.outcome, Decision::Exhausted { .. })),
        "the starved judge answered without the identity, so the locks above prove nothing"
    );
}

/// **…and nowhere else.** Each question's true answer is not zero, and each is the one a
/// specific wrong rule would answer `0`:
/// - the vertical line along z — the line runs along it; a rule that kept `axis` among the
///   columns would find the two rows' only nonzero column shared, and call it held;
/// - two points sharing only the wall turned off `x = 1`, along x — its world normal is
///   `(cos 30°, sin 30°, 0)`; a rule that trusted every pre-motion zero, whatever axis the chain
///   turns about, would read the wall as `x = const`;
/// - the line of a plane turned about x (`(0, c, s)`) and a wall turned about z (`(a, b, 0)`),
///   along z — its direction's z component is `−c·a`; of the two ways to give the rows columns
///   `{x, y}`, one lands on the first row's zero, and a rule that tried only that one would call
///   it held.
#[test]
fn a_coordinate_the_shared_planes_do_not_hold_is_left_to_the_judges() {
    let ps = identity_table();
    let free: [([usize; 3], [usize; 3], usize, &str); 3] = [
        ([0, 1, 2], [0, 1, 3], 2, "the vertical line, along z"),
        ([0, 1, 2], [0, 4, 3], 0, "the turned wall alone, along x"),
        (
            [5, 0, 2],
            [5, 0, 3],
            2,
            "a turned-about-x plane and a wall, along z",
        ),
    ];
    for (a, b, axis, what) in free {
        let want = escalated(&ps, a, b, axis);
        assert_ne!(want, 0, "{what}: the fixture is meant to differ there");
        assert_eq!(jd(&ps).cmp_coord(a, b, axis), want, "{what}");
    }
}

/// **A chain the rule cannot read is not read.** A wall turned about z and then reflected stays
/// vertical, so its line with another wall is vertical and the x comparison really is `0` — but
/// the reflection ends `single_axis_motion`, the rule knows none of the wall's zeros, and the
/// starved judge is left `Exhausted` rather than handed an answer the rule did not prove.
#[test]
fn a_reflected_chain_lends_the_identity_nothing() {
    let mirrored = RW::moved([[1, 0, 0], [1, 1, 0], [1, 0, 1]], |p| {
        p.rotate_about(Axis::Z, Angle::from_deg(ri(30, 1)).unwrap(), RW::pivot())
            .mirror(Axis::X, ri(0, 1))
    })
    .chain(6)
    .named();
    let mut ps = identity_table();
    ps[0] = mirrored;
    let (a, b) = ([0usize, 1, 2], [0usize, 1, 3]);
    assert_eq!(
        escalated(&ps, a, b, 0),
        0,
        "the reflected wall is still vertical"
    );
    let notes = Notes::new();
    let jd = Judge::new(&ps, starved(), &notes);
    jd.cmp_coord(a, b, 0);
    assert!(
        notes
            .sorted()
            .iter()
            .any(|e| matches!(e.outcome, Decision::Exhausted { .. })),
        "the identity answered across a reflection it cannot read"
    );
}

/// The escalation's direction sign, bridged to the stored normals as the judge bridges it.
fn escalated_dir(ps: &[RW], p: usize, a: usize, b: usize) -> i8 {
    let d = |k: usize| plane_def(ps, k);
    ps[p].frame_sign()
        * ps[a].frame_sign()
        * ps[b].frame_sign()
        * to_i8(dir_sign_judge(borrow3(d(p)), borrow3(d(a)), borrow3(d(b)), fixture()).orient())
}

/// **Three normals dependent by their indices or their known zeros give `0` outright** — a
/// plane named twice, and three walls turned about z by three angles (every normal's z is `0`).
/// The starved judge answers and records nothing; the escalation agrees. Two walls and a cap, or
/// two walls and a plane turned about x, leave each normal a column of its own, and there the
/// judges answer the sign they find.
#[test]
fn normals_dependent_by_their_zeros_answer_a_direction_outright() {
    let ps = identity_table();
    for (p, a, b, what) in [
        (0, 0, 2, "a wall named twice"),
        (0, 2, 2, "a cap named twice"),
        (0, 1, 4, "three walls turned about z"),
    ] {
        let notes = Notes::new();
        let jd = Judge::new(&ps, starved(), &notes);
        assert_eq!(jd.plane_pair_dir_sign(p, a, b), 0, "{what}");
        assert!(
            notes.sorted().is_empty(),
            "{what}: the starved judge climbed"
        );
        assert_eq!(
            escalated_dir(&ps, p, a, b),
            0,
            "{what}: the escalation disagrees"
        );
    }
    for (p, a, b, what) in [
        (0, 1, 2, "two walls and a cap"),
        (0, 1, 5, "two walls and a plane turned about x"),
        // First in line, the turned-about-x plane's first column is its zero.
        (5, 0, 1, "a plane turned about x and two walls"),
    ] {
        let want = escalated_dir(&ps, p, a, b);
        assert_ne!(want, 0, "{what}: the fixture is meant to be independent");
        assert_eq!(jd(&ps).plane_pair_dir_sign(p, a, b), want, "{what}");
    }
    // The control: the walls' zeros unknown, the same true zero leaves the starved judge
    // `Exhausted`.
    let unnamed = identity_planes(false);
    let notes = Notes::new();
    let jd = Judge::new(&unnamed, starved(), &notes);
    assert_eq!(jd.plane_pair_dir_sign(0, 1, 4), 0);
    assert!(
        notes
            .sorted()
            .iter()
            .any(|e| matches!(e.outcome, Decision::Exhausted { .. })),
        "the starved judge answered without the identity, so the locks above prove nothing"
    );
}
