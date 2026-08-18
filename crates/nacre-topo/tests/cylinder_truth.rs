//! M6-0 gates for the cylinder's exact truth (`CylinderDef`).
//!
//! The bit literals below are **pre-change measurements**: the same four fixtures were probed
//! on the code *before* `CylinderDef` existed, and their seam-vertex coordinates and lateral
//! cache fields recorded as bit patterns. M6-0 adds truth beside the cache without touching how
//! the cache is realized, so every literal must still match bit-for-bit — this file is the eye
//! the existing gates do not have (the bit census has no cylinder models, STEP asserts
//! structure, tess is property-based; none of them would see an ulp drift here).
//!
//! Fixture algebra, not story (the fixtures-that-measure-nothing rule):
//! - **z-axis**: the axis a cylinder is usually built on — exact normalization, `ref_dir`
//!   realization must be bit-identical.
//! - **pythagorean (3,4,0)**: irrational-looking but exactly normalizable (|·| = 5) — still
//!   bit-identical, and the raw truth keeps the caller's `(3,4,0)`, not the unit vector.
//! - **tilted (1,2,3)**: genuinely irrational normalization — the def realizes through a
//!   different arithmetic order than the cache, so agreement is ≤ 1 ulp, not bitwise.
//! - **near-tie (1,1,2)**: |x| == |y| exercises the tie-break (X before Y) that the rational
//!   rule must copy from `any_perpendicular` — a wrong tie-break turns the seam 90°, an
//!   O(radius) error no epsilon gate would forgive.

use nacre_math::{Point3, Vector3};
use nacre_scalar::Rat;
use nacre_topo::{CylinderDef, Model, SurfaceTruth, VertexDef};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::from_array([x, y, z])
}

fn vec(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::from_array([x, y, z])
}

fn rat(x: f64) -> Rat {
    Rat::from_decimal(x).expect("fixture decimals are in the window")
}

/// Build one cylinder and return (seam-vertex coords bottom/top, lateral cache, lateral def).
fn build(
    base: Point3,
    axis: Vector3,
    radius: f64,
    height: f64,
) -> (Point3, Point3, nacre_geom::Cylinder, CylinderDef) {
    let mut m = Model::new();
    let s = m.add_cylinder(base, axis, radius, height);
    let shell = m.solids.get(s).outer;
    let lateral = m
        .shells
        .get(shell)
        .faces
        .iter()
        .map(|&f| m.faces.get(f).surface)
        .find(|&su| matches!(m.surface(su), nacre_geom::Surface::Cylinder(_)))
        .expect("a cylinder solid has a lateral face");
    let cache = match m.surface(lateral) {
        nacre_geom::Surface::Cylinder(c) => *c,
        _ => unreachable!(),
    };
    let def = match m.surface_truth(lateral) {
        SurfaceTruth::Cylinder { def, motion } => {
            assert!(motion.is_none(), "construction states the world");
            def.clone()
        }
        _ => panic!("a lateral face carries a cylinder truth"),
    };
    // The two seam vertices, bottom first (lower vertex-store index).
    let mut seams: Vec<_> = m
        .vertices
        .iter()
        .filter(|(_, v)| matches!(v.def, VertexDef::OnSeam(_)))
        .map(|(h, _)| h)
        .collect();
    seams.sort_by_key(|h| h.index());
    assert_eq!(seams.len(), 2, "one seam vertex per rim");
    (
        m.vertex_point(seams[0]),
        m.vertex_point(seams[1]),
        cache,
        def,
    )
}

fn bits3(p: [f64; 3]) -> [u64; 3] {
    p.map(f64::to_bits)
}

/// The def's `ref_dir` realized to f64: normalize the exact direction's own realization.
fn realized_ref_dir(def: &CylinderDef) -> Vector3 {
    let r = def.ref_dir();
    Vector3::from_array([r[0].to_f64(), r[1].to_f64(), r[2].to_f64()])
        .normalize()
        .expect("ref_dir is nonzero by construction")
}

fn ulp_apart(a: f64, b: f64) -> u64 {
    a.to_bits().abs_diff(b.to_bits())
}

#[test]
fn the_z_axis_cylinder_is_bit_identical_and_its_truth_is_the_statement() {
    let (v0, v1, cache, def) = build(pt(0.5, -1.25, 2.0), vec(0.0, 0.0, 1.0), 1.5, 2.5);
    // Pre-change bit literals (probe 2026-08-17).
    assert_eq!(
        bits3(v0.as_array()),
        [
            4602678819172646912,
            13836746905142427648,
            4611686018427387904
        ]
    );
    assert_eq!(
        bits3(v1.as_array()),
        [
            4602678819172646912,
            13836746905142427648,
            4616752568008179712
        ]
    );
    assert_eq!(
        bits3(cache.axis().origin().as_array()),
        [
            4602678819172646912,
            13831680355561635840,
            4611686018427387904
        ]
    );
    assert_eq!(
        bits3(cache.axis().direction().as_array()),
        [0, 0, 4607182418800017408]
    );
    assert_eq!(
        bits3(cache.ref_dir().as_array()),
        [0, 13830554455654793216, 0]
    );
    assert_eq!(cache.radius().to_bits(), 4609434218613702656);
    // The truth is the raw statement: ẑ has smallest |x| (ties X→Y→Z pick X), x̂ × ẑ = (0,−1,0).
    assert_eq!(def.origin(), [rat(0.5), rat(-1.25), rat(2.0)]);
    assert_eq!(def.dir(), [rat(0.0), rat(0.0), rat(1.0)]);
    assert_eq!(def.ref_dir(), [rat(0.0), rat(-1.0), rat(0.0)]);
    assert_eq!(def.radius(), rat(1.5));
    // Exactly normalizable → the def's realization is the cache's ref_dir, bit for bit.
    assert_eq!(
        bits3(realized_ref_dir(&def).as_array()),
        bits3(cache.ref_dir().as_array())
    );
}

#[test]
fn a_pythagorean_axis_keeps_the_raw_statement_and_the_exact_seam() {
    let (v0, v1, cache, def) = build(pt(1.0, 2.0, 3.0), vec(3.0, 4.0, 0.0), 2.0, 5.0);
    assert_eq!(
        bits3(v0.as_array()),
        [
            13826951575952896820,
            4614388178203810202,
            4613937818241073152
        ]
    );
    assert_eq!(
        bits3(v1.as_array()),
        [
            4612586738352862003,
            4619792497756654797,
            4613937818241073152
        ]
    );
    assert_eq!(
        bits3(cache.ref_dir().as_array()),
        [13828753015803845018, 4603579539098121011, 0]
    );
    // Raw truth: the caller's (3,4,0), never the unit (0.6,0.8,0). Smallest |component| is z;
    // ẑ × (3,4,0) = (−4,3,0).
    assert_eq!(def.dir(), [rat(3.0), rat(4.0), rat(0.0)]);
    assert_eq!(def.ref_dir(), [rat(-4.0), rat(3.0), rat(0.0)]);
    // |(−4,3,0)| = 5 exactly → realization is bit-identical to the cache.
    assert_eq!(
        bits3(realized_ref_dir(&def).as_array()),
        bits3(cache.ref_dir().as_array())
    );
}

#[test]
fn a_tilted_axis_realizes_the_same_seam_within_one_ulp() {
    let (v0, v1, cache, def) = build(pt(3.0, -1.0, 2.0), vec(1.0, 2.0, 3.0), 0.7, 4.0);
    assert_eq!(
        bits3(v0.as_array()),
        [
            4613937818241073152,
            13833177510631666613,
            4612560370086345703
        ]
    );
    assert_eq!(
        bits3(v1.as_array()),
        [
            4616267355777403146,
            4603180112408586558,
            4617985906959014164
        ]
    );
    assert_eq!(
        bits3(cache.ref_dir().as_array()),
        [0, 13829041699191119073, 4603171514739320982]
    );
    // Smallest |component| is x; x̂ × (1,2,3) = (0,−3,2) — the exact seam direction.
    assert_eq!(def.ref_dir(), [rat(0.0), rat(-3.0), rat(2.0)]);
    // Irrational normalization: the def realizes through a different arithmetic order than the
    // cache (Gram–Schmidt on the rounded unit axis), so agreement is ≤ 1 ulp per component —
    // and the direction itself is exact in the def, which is the point of M6-0.
    let r = realized_ref_dir(&def);
    let c = cache.ref_dir();
    for k in 0..3 {
        assert!(
            ulp_apart(r[k], c[k]) <= 1,
            "component {k}: def realization {} vs cache {} differ by more than 1 ulp",
            r[k],
            c[k]
        );
    }
    assert!(r.dot(c) > 0.0, "the seam direction must be preserved");
}

#[test]
fn a_near_tie_axis_takes_the_same_tie_break_as_the_cache() {
    let (v0, v1, cache, def) = build(pt(0.0, 0.0, 0.0), vec(1.0, 1.0, 2.0), 1.0, 1.0);
    assert_eq!(
        bits3(v0.as_array()),
        [0, 13829603540328246745, 4601727903846100441]
    );
    assert_eq!(
        bits3(v1.as_array()),
        [
            4601025967313136703,
            13825801877233839987,
            4608370063852310934
        ]
    );
    // Measured before the change: the cache's ref_dir and the seam vertex already disagree by
    // 1 ulp today (two normalization paths) — the literal records that honestly.
    assert_eq!(
        bits3(cache.ref_dir().as_array()),
        [0, 13829603540328246746, 4601727903846100442]
    );
    // |x| == |y|: the tie must break X before Y, exactly as `any_perpendicular` does —
    // x̂ × (1,1,2) = (0,−2,1). A Y pick would give (2,0,−1): a seam turned ~90°.
    assert_eq!(def.ref_dir(), [rat(0.0), rat(-2.0), rat(1.0)]);
    let r = realized_ref_dir(&def);
    let c = cache.ref_dir();
    for k in 0..3 {
        assert!(ulp_apart(r[k], c[k]) <= 1, "component {k} beyond 1 ulp");
    }
    assert!(r.dot(c) > 0.0, "the seam direction must be preserved");
}

#[test]
fn a_tie_born_of_rounding_picks_the_caches_basis() {
    // ★ The axis that refuted the first spelling's order argument ("one positive scale
    // preserves |·| order" — a real-number claim). Raw components: |x| > |y| strictly, so a
    // raw-reading rule picks basis Y; but f64 normalization rounds both to the SAME value, and
    // `any_perpendicular` (reading the normalized axis) tie-breaks to X. The rule must read
    // `d`'s components — the cache's actual inputs — or the seam lands ~90° from the cache's
    // and the validate net fires on a healthy model (measured, 2026-08-17).
    let (_, _, cache, def) = build(
        pt(0.0, 0.0, 0.0),
        vec(0.34, 0.33999999999999997, 1.0),
        0.5,
        1.0,
    );
    // X basis: x̂ × dir = (0, −dir_z, dir_y) — the zero slot names the chosen basis.
    assert_eq!(
        def.ref_dir(),
        [rat(0.0), rat(-1.0), rat(0.33999999999999997)],
        "the basis choice must be the cache's (X), not the raw order's (Y)"
    );
    let r = realized_ref_dir(&def);
    let c = cache.ref_dir();
    for k in 0..3 {
        assert!(ulp_apart(r[k], c[k]) <= 1, "component {k} beyond 1 ulp");
    }
    assert!(r.dot(c) > 0.0, "the seam direction must be preserved");
}

#[test]
fn the_same_statement_interns_and_a_different_ref_dir_does_not() {
    let mut m = Model::new();
    let cache = nacre_geom::Cylinder::from_axis(
        pt(0.0, 0.0, 0.0),
        vec(0.0, 0.0, 1.0),
        vec(0.0, -1.0, 0.0),
        1.0,
    )
    .expect("non-degenerate");
    let zero = || rat(0.0);
    let def = CylinderDef::new(
        [zero(), zero(), zero()],
        [zero(), zero(), rat(1.0)],
        [zero(), rat(-1.0), zero()],
        rat(1.0),
    )
    .expect("non-degenerate");
    let a = m.push_cylinder(cache, def.clone(), None);
    let b = m.push_cylinder(cache, def, None);
    assert_eq!(a, b, "one statement, one handle");
    // Same axis and radius, different ref_dir: a merge would split the seam, so the
    // conservative key deliberately keeps two handles (geometric identity is M6-1's,
    // per predicate).
    let other = CylinderDef::new(
        [zero(), zero(), zero()],
        [zero(), zero(), rat(1.0)],
        [rat(1.0), zero(), zero()],
        rat(1.0),
    )
    .expect("non-degenerate");
    let c = m.push_cylinder(cache, other, None);
    assert_ne!(a, c, "a different seam statement is a different surface");
}

#[test]
fn the_checked_constructor_refuses_what_means_no_cylinder() {
    let zero = || rat(0.0);
    let o = [zero(), zero(), zero()];
    let z = [zero(), zero(), rat(1.0)];
    let r = [rat(1.0), zero(), zero()];
    assert!(CylinderDef::new(o, o, r, rat(1.0)).is_none(), "zero axis");
    assert!(CylinderDef::new(o, z, r, zero()).is_none(), "zero radius");
    assert!(
        CylinderDef::new(o, z, r, rat(-1.0)).is_none(),
        "negative radius"
    );
    assert!(
        CylinderDef::new(o, z, [zero(), zero(), rat(2.0)], rat(1.0)).is_none(),
        "ref_dir parallel to the axis pins no seam"
    );
    assert!(
        CylinderDef::new(o, z, r, rat(1.0)).is_some(),
        "the sane statement stands"
    );
}

/// **A statement the arithmetic used to lose.** The axis carries a component that is small and
/// spelled with a full f64's digits, so its exact rational has a ~10²⁰ denominator; the
/// parallelism test squares it, and in `i128` that overflowed. The old constructor answered
/// `None` — "no cylinder" — and `add_cylinder`'s `expect` turned a width limit into a crash.
///
/// The three inputs are one family, chosen from the algebra rather than a story: a tiny
/// component at 1e-7 (the seed a proptest actually found), one at 1e-9, and one where **both**
/// off-axis components are wide, so the cross has no zero term to hide behind.
#[test]
fn a_wide_decimal_axis_is_a_cylinder_not_a_refusal() {
    let zero = || rat(0.0);
    let o = [zero(), zero(), zero()];
    let r = rat(1.0);
    for axis in [
        [rat(2.088798035473136e-7), zero(), rat(0.7055489621854671)],
        [rat(3.141592653589793e-9), zero(), rat(0.8414709848078965)],
        [
            rat(1.414213562373095e-8),
            rat(2.7182818284590453e-6),
            rat(0.9092974268256817),
        ],
    ] {
        // `ref_dir` as `add_cylinder` builds it for a nearly-`ẑ` axis: `ê_x × axis`.
        let neg = |x: Rat| {
            Rat::from_int(0)
                .checked_sub(x)
                .expect("negating a lifted decimal")
        };
        let ref_dir = [zero(), neg(axis[2]), axis[1]];
        assert!(
            CylinderDef::new(o, axis, ref_dir, r).is_some(),
            "a wide-decimal axis states a cylinder: {axis:?}"
        );
    }
}

/// The negative control for the totalization: making the test unable to overflow must not make
/// it unable to **refuse**. A `ref_dir` parallel to the axis still pins no seam when both are
/// written with wide decimals — where the old checked test could only shrug.
#[test]
fn a_wide_parallel_ref_dir_is_still_refused() {
    let zero = || rat(0.0);
    let o = [zero(), zero(), zero()];
    let axis = [rat(2.088798035473136e-7), zero(), rat(0.7055489621854671)];
    let doubled = [
        axis[0].checked_add(axis[0]).expect("small doubling"),
        zero(),
        axis[2].checked_add(axis[2]).expect("small doubling"),
    ];
    assert!(
        CylinderDef::new(o, axis, doubled, rat(1.0)).is_none(),
        "a parallel ref_dir pins no seam, however wide its spelling"
    );
    assert!(
        CylinderDef::new(o, axis, [zero(), zero(), zero()], rat(1.0)).is_none(),
        "a zero ref_dir pins no seam either"
    );
}

// ── K2: the exact entry ───────────────────────────────────────────────────────────────────────
//
// Everything above measures the f64 entry, which *derives* its truth by lifting computed floats.
// `add_cylinder_exact` is handed the truth instead, and these gates measure what that buys:
// statements no decimal window can hold, refusals with names rather than panics, and — where the
// two roads can say the same thing — the same model.

fn int3(v: [i128; 3]) -> [Rat; 3] {
    v.map(Rat::from_int)
}

/// The three faces `add_cylinder_exact` hands back, in push order.
type CylFaces = [nacre_store::Handle<nacre_topo::Face>; 3];

/// The lateral surface's truth and the two seam coordinates of a model built exactly.
fn read_exact(m: &Model, faces: CylFaces) -> (CylinderDef, Point3, Point3) {
    let lateral = m.faces.get(faces[0]).surface;
    let def = match m.surface_truth(lateral) {
        SurfaceTruth::Cylinder { def, motion } => {
            assert!(motion.is_none(), "this fixture states the world");
            def.clone()
        }
        _ => panic!("the first face of a cylinder solid is its lateral"),
    };
    let mut seams: Vec<_> = m
        .vertices
        .iter()
        .filter(|(_, v)| matches!(v.def, VertexDef::OnSeam(_)))
        .map(|(h, _)| h)
        .collect();
    seams.sort_by_key(|h| h.index());
    assert_eq!(seams.len(), 2, "one seam vertex per rim");
    (def, m.vertex_point(seams[0]), m.vertex_point(seams[1]))
}

/// ★ **The two roads are one spelling.** For a z-axis statement the f64 entry derives
/// `ref_dir = ê_x × axis = (0, −1, 0)`; handing the exact entry that same truth must produce the
/// same model — same def, same seam coordinates, bit for bit. If this ever parts, the shared
/// b-rep body grew a second version of a rule (the one thing extracting it was meant to prevent).
#[test]
fn both_roads_state_the_same_z_axis_cylinder() {
    let (bot_f64, top_f64, _, def_f64) = build(pt(0.0, 0.0, 0.0), vec(0.0, 0.0, 1.0), 2.0, 5.0);

    let mut m = Model::new();
    let (_, faces) = m
        .add_cylinder_exact(
            int3([0, 0, 0]),
            int3([0, 0, 1]),
            int3([0, -1, 0]),
            Rat::from_int(2),
            Rat::from_int(5),
            None,
        )
        .expect("an orthonormal statement builds");
    let (def, bot, top) = read_exact(&m, faces);

    assert_eq!(def.origin(), def_f64.origin(), "same axis point");
    assert_eq!(def.dir(), def_f64.dir(), "same axis direction");
    assert_eq!(def.ref_dir(), def_f64.ref_dir(), "same seam reference");
    assert_eq!(def.radius(), def_f64.radius(), "same radius");
    assert_eq!(bits3(bot.as_array()), bits3(bot_f64.as_array()));
    assert_eq!(bits3(top.as_array()), bits3(top_f64.as_array()));
}

/// ★★ **A third is a cylinder.** `1/3` has no decimal spelling, so the f64 entry cannot state
/// this base at all: it would lift `Rat::from_decimal(0.333…)`, a *different* point. The exact
/// entry keeps what it was given — and the second assertion is what gives the first its meaning
/// (without it this test would pass on a road that quietly rounded).
#[test]
fn a_base_with_no_decimal_form_is_stated_exactly() {
    let third = Rat::new(1, 3).expect("nonzero denominator");
    let mut m = Model::new();
    let (_, faces) = m
        .add_cylinder_exact(
            [third, Rat::from_int(0), Rat::from_int(0)],
            int3([0, 0, 1]),
            int3([1, 0, 0]),
            Rat::from_int(1),
            Rat::from_int(4),
            None,
        )
        .expect("a rational base is a statement, not a decimal");
    let (def, _, _) = read_exact(&m, faces);
    assert_eq!(def.origin()[0], third, "the truth is a third, exactly");
    assert_ne!(
        rat(third.to_f64()),
        third,
        "the decimal road's nearest statement is a different point — the window is why this \
         entry takes rationals"
    );
}

/// ★★ **Every refusal has a name, and none of them touches the arena.** The f64 entry answers
/// these four with `expect`/`debug_assert`; here they are values, because an application's
/// numbers are input and a panic in wasm is a dead session. The store-length assertion is the
/// other half: a late refusal would leave cells behind and shift every later log index.
#[test]
fn a_refused_statement_is_named_and_leaves_nothing_behind() {
    use nacre_topo::CylinderError;

    let mut m = Model::new();
    let before = (
        m.surface_count(),
        m.vertices.iter().count(),
        m.edges.iter().count(),
        m.live_solids.len(),
    );
    let two = Rat::from_int(2);
    let five = Rat::from_int(5);
    /// A refused statement: axis, seam reference, radius, height — and the name it earns.
    type Case = ([Rat; 3], [Rat; 3], Rat, Rat, CylinderError);
    let cases: [Case; 5] = [
        // A raw axis: the right direction, but not unit — every derived point would be scaled.
        (
            int3([0, 0, 2]),
            int3([1, 0, 0]),
            two,
            five,
            CylinderError::FrameNotOrthonormal,
        ),
        // Unit, but the seam reference lies along the axis: it pins no angle.
        (
            int3([0, 0, 1]),
            int3([0, 0, 1]),
            two,
            five,
            CylinderError::FrameNotOrthonormal,
        ),
        // Unit and perpendicular directions, but a non-unit reference length.
        (
            int3([0, 0, 1]),
            int3([3, 0, 0]),
            two,
            five,
            CylinderError::FrameNotOrthonormal,
        ),
        (
            int3([0, 0, 1]),
            int3([1, 0, 0]),
            Rat::from_int(0),
            five,
            CylinderError::NonPositiveRadius,
        ),
        (
            int3([0, 0, 1]),
            int3([1, 0, 0]),
            two,
            Rat::from_int(-5),
            CylinderError::NonPositiveHeight,
        ),
    ];
    for (axis, ref_dir, radius, height, want) in cases {
        assert_eq!(
            m.add_cylinder_exact(int3([0, 0, 0]), axis, ref_dir, radius, height, None)
                .unwrap_err(),
            want
        );
    }
    assert_eq!(
        (
            m.surface_count(),
            m.vertices.iter().count(),
            m.edges.iter().count(),
            m.live_solids.len()
        ),
        before,
        "a refusal is decided before anything is pushed"
    );
}
