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
    let a = m.push_cylinder(cache, def.clone(), None);
    let b = m.push_cylinder(cache, def, None);
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
    let c = m.push_cylinder(cache, other, None);
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
