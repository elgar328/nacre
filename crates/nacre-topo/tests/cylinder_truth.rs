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
