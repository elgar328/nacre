/// ★★★★ **Nothing carries plane coefficients through a frame, and that is the design, not an
/// omission.**
///
/// A plane that passes through a motion has two descriptions — coefficients and a witness
/// triangle — and mixing them is what once produced intransitive comparisons and a confident
/// wrong winding. The frame design sidesteps the whole question: a plane built in a frame
/// states itself **in that frame**, `BaseFrame` keeps its triangle and its coefficients in the
/// *same* frame, and every judgement downstream is a determinant question a rigid motion
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

/// How these fixtures judge; production chooses both per model. The coincidence limit is the
/// derived default for a unit-scale model — output resolution (`2⁻⁵²`) two words further down.
fn fixture() -> Standard {
    Standard {
        prec: 256,
        coincidence: Mag::pow2(-180),
        scale: Mag::of(1.0),
        cap: 4096,
    }
}

/// A synthetic axis-aligned plane witness: three points on the plane, its exact
/// coefficients, and the exact `WitnessPoint` definition of those points. `is_rotated` is `false`,
/// so predicates take the exact path — the definition is there but unused, which is
/// precisely the arrangement the production tables now have.
struct W {
    tri: [Point3; 3],
    coeffs: [f64; 4],
    def: [WitnessPoint; 3],
}
impl W {
    fn new(tri: [Point3; 3], coeffs: [f64; 4]) -> W {
        let def = tri.map(|p| {
            WitnessPoint::at_nearest(
                p.as_array()
                    .map(|x| nacre_exact::Rat::try_from_f64(x).expect("test coordinate")),
            )
        });
        W { tri, coeffs, def }
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
    // These witnesses carry no motion, so there is nothing to cancel.
    fn chain_id(&self) -> u64 {
        0
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        None
    }
}
impl PlaneWitness for W {
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

/// A plane witness carrying one shared rotation — what [`cancel_cmp_coord`] needs.
struct RW {
    tri: [Point3; 3],
    coeffs: [f64; 4],
    def: [WitnessPoint; 3],
    base: [Point3; 3],
    base_coeffs: [f64; 4],
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
        let pivot = [ri(1, 3), ri(1, 7), ri(0, 1)];
        let def: [WitnessPoint; 3] = pts.map(|q| {
            let mut p = WitnessPoint::at([ri(q[0], 1), ri(q[1], 1), ri(q[2], 1)]);
            for &(axis, deg) in nodes {
                p = p.rotate_about(axis, Angle::from_deg(ri(deg, 1)).unwrap(), pivot);
            }
            p
        });
        let base = pts.map(|q| Point3::from_array([q[0] as f64, q[1] as f64, q[2] as f64]));
        let tri: [Point3; 3] = std::array::from_fn(|i| Point3::from_array(def[i].coord()));
        RW {
            coeffs: plane_of(tri.map(|p| p.as_array())),
            base_coeffs: plane_of(base.map(|p| p.as_array())),
            tri,
            def,
            base,
        }
    }
}
impl Witness for RW {
    // A hand-built witness has no name: every same-plane question takes the definitions' road.
    fn world_name(&self) -> Option<&nacre_exact::PlaneName> {
        None
    }
    fn tri_pt3(&self) -> &[WitnessPoint; 3] {
        &self.def
    }
    fn is_rotated(&self) -> bool {
        true
    }
    fn chain_id(&self) -> u64 {
        1 // one shared motion for every witness in these fixtures
    }
    fn base_tri(&self) -> Option<[Point3; 3]> {
        Some(self.base)
    }
}
impl PlaneWitness for RW {
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
/// [`cancel_cmp_coord`] answers in the pre-rotation frame with exact rational arithmetic;
/// [`indirect_cmp_coord_judge`] answers by realizing the rotation in astro-float and reading
/// an interval. They share no bound and no code, so agreement is a real cross-check — and it
/// is the only one available, because the shortcut's whole point is to *not* do what the
/// judge does.
///
/// The fixtures turn about a **non-origin pivot**, which is the step the original comment
/// missed: a rotation about `c` is `x ↦ R(x−c)+c`, and the `+c` cancels in a difference. If
/// it did not, every answer here would be wrong.
#[test]
fn the_pre_rotation_shortcut_agrees_with_the_escalation() {
    for deg in [30i128, 37, 120, 200] {
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            // ∩(x=1, y=2, z=3) = (1,2,3) and ∩(x=1, y=2, z=7) = (1,2,7): the difference is
            // along z, so an in-plane comparison is exactly 0 and a z comparison is definite.
            // Then a pair differing along x, which the shortcut must decline rather than
            // guess (`tan Θ` is irrational, so it is provably nonzero but not signed here).
            let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg);
            let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg);
            let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg);
            let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg);
            let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg);
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
/// with itself. The rotation axis is answered exactly, and so is a difference lying along it.
#[test]
fn the_shortcut_fires_where_it_should_and_declines_where_it_cannot() {
    let (axis, deg) = (Axis::Z, 30i128);
    let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg);
    let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg);
    let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg);
    let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg);
    let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg);
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
/// mixed-axis fixture nothing else in the suite notices if it is removed.
#[test]
fn a_chain_about_two_axes_is_declined() {
    let nodes: &[(Axis, i128)] = &[(Axis::X, 30), (Axis::Z, 40)];
    let t = |pts| RW::turned(pts, nodes);
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

/// **A name's `f64` row stands exactly when every integer it holds fits 53 bits.** The row is
/// what an unmoved plane's exact shortcuts read, so a row that stood for a wider integer would
/// hand the predicates a value that is not the name — or, for a power of two that *is*
/// representable, a magnitude their `Expansion` products overflow on.
#[test]
fn a_names_f64_row_stands_only_where_every_integer_fits_53_bits() {
    use super::{WitnessPoint, name_stored_ints};
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
    // `a·x − z = 0` for `a = 2⁵³ + 1` (not an `f64`) and `a = 2⁶⁰` (an `f64`, too wide anyway).
    for a in [(1i128 << 53) + 1, 1i128 << 60] {
        let wide = ints([
            [r(0, 1), r(0, 1), r(0, 1)],
            [r(0, 1), r(1, 1), r(0, 1)],
            [r(1, 1), r(0, 1), r(a, 1)],
        ]);
        assert_eq!((wide.row, wide.normal), (None, None), "a = {a}");
    }
    // `z = 2⁶⁰ + 1`: only `d` is wide, so the normal stands and the row does not.
    let far = (1i128 << 60) + 1;
    let high = ints([
        [r(0, 1), r(0, 1), r(far, 1)],
        [r(1, 1), r(0, 1), r(far, 1)],
        [r(0, 1), r(1, 1), r(far, 1)],
    ]);
    assert_eq!(high.row, None);
    assert_eq!(high.normal.map(|n| n.map(f64::abs)), Some([0.0, 0.0, 1.0]));
}
