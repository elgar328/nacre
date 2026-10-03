//! Gates for the cylinder's exact truth (`CylinderDef`) at the kernel's own door: what a
//! statement is (the checked constructor refuses what names no cylinder, and only that — however
//! wide the decimals), and what one statement is in the arena (`push_cylinder` interns by the
//! whole statement, `ref_dir` included).

use nacre_exact::Rat;
use nacre_math::{Point3, Vector3};
use nacre_topo::{CylinderDef, Model};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::from_array([x, y, z])
}

fn vec(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::from_array([x, y, z])
}

fn rat(x: f64) -> Rat {
    Rat::from_decimal(x).expect("fixture decimals are in the window")
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
        nacre_exact::BigRat::from(rat(1.0)),
    )
    .expect("non-degenerate");
    let a = m.push_cylinder(
        cache,
        nacre_topo::CacheStanding::Unrealized,
        def.clone(),
        None,
    );
    let b = m.push_cylinder(cache, nacre_topo::CacheStanding::Unrealized, def, None);
    assert_eq!(a, b, "one statement, one handle");
    // Same axis and radius, different ref_dir: a merge would split the seam, so the
    // conservative key deliberately keeps two handles (geometric identity is decided
    // per predicate, not by the key).
    let other = CylinderDef::new(
        [zero(), zero(), zero()],
        [zero(), zero(), rat(1.0)],
        [rat(1.0), zero(), zero()],
        nacre_exact::BigRat::from(rat(1.0)),
    )
    .expect("non-degenerate");
    // Its own cache — the door holds a producer's cache to the statement it comes with.
    let other_cache = nacre_geom::Cylinder::from_axis(
        pt(0.0, 0.0, 0.0),
        vec(0.0, 0.0, 1.0),
        vec(1.0, 0.0, 0.0),
        1.0,
    )
    .expect("non-degenerate");
    let c = m.push_cylinder(
        other_cache,
        nacre_topo::CacheStanding::Unrealized,
        other,
        None,
    );
    assert_ne!(a, c, "a different seam statement is a different surface");
}

#[test]
fn the_checked_constructor_refuses_what_means_no_cylinder() {
    let zero = || rat(0.0);
    let o = [zero(), zero(), zero()];
    let z = [zero(), zero(), rat(1.0)];
    let r = [rat(1.0), zero(), zero()];
    assert!(
        CylinderDef::new(o, o, r, nacre_exact::BigRat::from(rat(1.0))).is_none(),
        "zero axis"
    );
    assert!(
        CylinderDef::new(o, z, r, nacre_exact::BigRat::from(zero())).is_none(),
        "zero radius"
    );
    assert!(
        CylinderDef::new(o, z, r, nacre_exact::BigRat::from(rat(-1.0))).is_none(),
        "negative radius"
    );
    assert!(
        CylinderDef::new(
            o,
            z,
            [zero(), zero(), rat(2.0)],
            nacre_exact::BigRat::from(rat(1.0))
        )
        .is_none(),
        "ref_dir parallel to the axis pins no seam"
    );
    assert!(
        CylinderDef::new(o, z, r, nacre_exact::BigRat::from(rat(1.0))).is_some(),
        "the sane statement stands"
    );
}

/// **A statement checked arithmetic loses.** The axis carries a component that is small and
/// spelled with a full f64's digits, so its exact rational has a ~10²⁰ denominator; the
/// parallelism test squares it, and in `i128` that overflows — a checked constructor would
/// answer `None` ("no cylinder"), turning a width limit into a refusal of a sane statement.
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
        // A seam reference square to a nearly-`ẑ` axis: the basis cross `ê_x × axis`.
        let neg = |x: Rat| {
            Rat::from_int(0)
                .checked_sub(x)
                .expect("negating a lifted decimal")
        };
        let ref_dir = [zero(), neg(axis[2]), axis[1]];
        assert!(
            CylinderDef::new(o, axis, ref_dir, nacre_exact::BigRat::from(r)).is_some(),
            "a wide-decimal axis states a cylinder: {axis:?}"
        );
    }
}

/// The negative control for the totalization: making the test unable to overflow must not make
/// it unable to **refuse**. A `ref_dir` parallel to the axis still pins no seam when both are
/// written with wide decimals — where a checked test could only shrug.
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
        CylinderDef::new(o, axis, doubled, nacre_exact::BigRat::from(rat(1.0))).is_none(),
        "a parallel ref_dir pins no seam, however wide its spelling"
    );
    assert!(
        CylinderDef::new(
            o,
            axis,
            [zero(), zero(), zero()],
            nacre_exact::BigRat::from(rat(1.0))
        )
        .is_none(),
        "a zero ref_dir pins no seam either"
    );
}

/// **The cache is the statement's correct rounding, not the producer's figure.** Two producers'
/// figures that miss it in the last place — a tilted axis normalized in `f64` (`(1, 1, 1)` and
/// `(1, −1, 0)` are two whose `f64` normalization is not the nearest unit vector)
/// ([`nacre_geom::Cylinder::from_axis`]), and an origin moved in `f64` beside a statement moved
/// exactly — are replaced at the door by the statement's origin rounded once and its unit frame
/// rounded once ([`nacre_exact::cyl_unit_frame_f64`]).
#[test]
fn the_door_realizes_a_cylinder_cache_from_its_statement() {
    let mut m = Model::new();
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let dir = [q(1, 1), q(1, 1), q(1, 1)];
    let ref_dir = [q(1, 1), q(-1, 1), q(0, 1)];
    let origin = [q(1, 3), q(2, 7), q(5, 11)];
    let def = CylinderDef::new(origin, dir, ref_dir, nacre_exact::BigRat::from(q(9, 4))).unwrap();
    // The producer's figure: every part computed in `f64`, and the origin a step away from the
    // statement's rounding, as an `f64` transport leaves it.
    let figure = nacre_geom::Cylinder::from_axis(
        pt(1.0 / 3.0 + 1e-15, 2.0 / 7.0, 5.0 / 11.0),
        vec(1.0, 1.0, 1.0),
        vec(1.0, -1.0, 0.0),
        1.5,
    )
    .unwrap();
    let h = m.push_cylinder(figure, nacre_topo::CacheStanding::Unrealized, def, None);
    let nacre_geom::Surface::Cylinder(c) = m.surface_cache(h) else {
        panic!("a cylinder")
    };
    let (axis, across) = nacre_exact::cyl_unit_frame_f64(&dir, &ref_dir).unwrap();
    assert_eq!(c.axis().origin().as_array(), origin.map(|r| r.to_f64()));
    assert_eq!(c.axis().direction().as_array(), axis);
    assert_eq!(c.ref_dir().as_array(), across);
    assert_eq!(c.radius(), 1.5);
    // The producer's own figure is not the rounding — the case the door exists for, on the
    // origin and on the frame alike.
    assert_ne!(
        figure.axis().origin().as_array(),
        origin.map(|r| r.to_f64())
    );
    assert_ne!(
        (
            figure.axis().direction().as_array(),
            figure.ref_dir().as_array()
        ),
        (axis, across)
    );
}
