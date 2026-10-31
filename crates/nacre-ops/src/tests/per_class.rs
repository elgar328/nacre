//! Unit tests of the per-class stage's own functions.

use super::*;
use nacre_exact::Rat;

/// **The band a closed trace sweeps along its own axis.** The cylinder stands on the `z` axis,
/// radius 2. On `z = 5` the trace is the circle at station 5. On `z = y/4 + 5/2` it is an ellipse
/// over `z ∈ [2, 3]`, `t ± Δ` with `Δ = 2·(1/4)/1`. With the raw axis `(0, 0, 2)` the stations
/// halve and so does `Δ` — `|n·m|` doubles. With radius `√2` and `z = y + 5` the band is
/// `5 ± √2`, rounded outward to rationals. A plane along the axis sweeps no band.
#[test]
fn a_traces_band_is_its_axis_meet_widened_by_the_ellipses_reach() {
    let q = |n: i128, d: i128| Rat::new(n, d).unwrap();
    let z = |v: i128| Rat::from_int(v);
    let cylinder = |dir: [Rat; 3], r2: Rat| {
        nacre_topo::CylinderDef::new(
            [z(0); 3],
            dir,
            [z(1), z(0), z(0)],
            nacre_exact::BigRat::from(r2),
        )
        .unwrap()
    };
    let up = cylinder([z(0), z(0), z(1)], z(4));
    assert_eq!(
        trace_band(&[z(0), z(0), z(1), z(-5)], &up),
        Some([z(5), z(5)])
    );
    let slant = [z(0), q(-1, 4), z(1), q(-5, 2)];
    assert_eq!(trace_band(&slant, &up), Some([z(2), z(3)]));
    let long = cylinder([z(0), z(0), z(2)], z(4));
    assert_eq!(trace_band(&slant, &long), Some([z(1), q(3, 2)]));
    let thin = cylinder([z(0), z(0), z(1)], z(2));
    let [lo, hi] = trace_band(&[z(0), z(-1), z(1), z(-5)], &thin).unwrap();
    let (below, above) = (z(5).checked_sub(lo).unwrap(), hi.checked_sub(z(5)).unwrap());
    for d in [below, above] {
        assert!(d.checked_mul(d).unwrap() >= z(2), "never narrower than √2");
        assert!(d < q(1415, 1000), "and near it");
    }
    assert_eq!(
        trace_band(&[z(1), z(0), z(0), z(-1)], &up),
        None,
        "along the axis"
    );
}

/// **The net cannot measure an oblique section its pieces do not part, and says so.** The bored
/// plate with a blind bore along `x` — radius ½ about `y = 10, z = 5`, `x ∈ [13, 40]` — and the sunk
/// gusset, whose slanted plane `z = y/4 + 5/2` sections the `z` bore and runs through the `x` bore's
/// axis: the class carries the first's ellipse and the second's rulings. With the ruling's own ends
/// the two pieces part along `x` (`[8, 12]` against `[13, 40]`) and the net answers «no crossing».
/// With its ends moved out to the plate's wall `x = 0` they overlap along every direction the pair
/// rule asks — the two axes meet at `(10, 10, 5)` — and the strip's arc ends are a circle's, so the
/// answer is «not measured», whether the circle states no contribution or one whole one: refused as
/// the oblique coverage's end, not as the net's backstop.
#[test]
fn an_oblique_section_its_pieces_do_not_part_is_not_measured() {
    let mut m = Model::new();
    let plate = crate::tests::plate_with_bore(&mut m);
    let across = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([13.0, 10.0, 5.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        28.0,
    )
    .solid;
    m.rebuild_adjacency();
    let a = crate::boolean(&mut m, BoolKind::Cut, plate, across).expect("the blind bore")[0];
    m.rebuild_adjacency();
    let b = crate::tests::gusset(&mut m, crate::tests::SUNK, 30.0);
    m.rebuild_adjacency();
    let setup = crate::arrangement::plane_index_setup(&m, a, b).expect("the gate passes");
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let z = |v: i128| Rat::from_int(v);
    let class = |want: [Rat; 4]| {
        (0..setup.geom.len())
            .find(|&c| {
                jd.world_rat(c).is_some_and(|w| {
                    // The same plane up to a positive or negative scale.
                    let k = (0..4).find(|&i| want[i] != z(0)).unwrap();
                    let s = w[k]
                        .checked_mul(Rat::new(want[k].denom(), want[k].numer()).unwrap())
                        .unwrap();
                    (0..4).all(|i| w[i] == want[i].checked_mul(s).unwrap())
                })
            })
            .expect("the class is in the table")
    };
    let wc = class([
        z(0),
        Rat::new(-1, 4).unwrap(),
        z(1),
        Rat::new(-5, 2).unwrap(),
    ]);
    let (x0, x13, x40) = (
        class([z(1), z(0), z(0), z(0)]),
        class([z(1), z(0), z(0), z(-13)]),
        class([z(1), z(0), z(0), z(-40)]),
    );
    let cyl_along = |dir: [Rat; 3]| {
        setup
            .cyls
            .iter()
            .position(|cy| nacre_exact::parallel_rat(&cy.def.dir(), &dir))
            .expect("the cylinder is in the table")
    };
    let (upright, lying) = (cyl_along([z(0), z(0), z(1)]), cyl_along([z(1), z(0), z(0)]));
    assert!(
        setup.cyls[upright].sections_oblique(wc),
        "the gate records the section"
    );
    assert!(
        setup.cyls[lying].cuts(wc),
        "and the cut along the other axis"
    );
    let end = |c: usize| {
        crate::combinatorics::NodeId::pierce(wc, c, lying, crate::combinatorics::QuadRoot::Lo)
    };
    let ruling = |from: usize| MergedRuling {
        cyl: lying,
        side: 1,
        end: [end(from), end(x40)],
        merged: Vec::new(),
        orient: -1,
    };
    let whole = (
        crate::planes::SolidSide::A,
        SegKind::Transversal { mat: -1 },
        None,
    );
    for merged in [Vec::new(), vec![whole]] {
        let circle = MergedCircle {
            cyl: upright,
            merged,
        };
        assert_eq!(
            circle_crosses_ruling(&jd, &setup.cyls, wc, &circle, &ruling(x13)),
            Ok(false),
            "parted along x"
        );
        assert_eq!(
            circle_crosses_ruling(&jd, &setup.cyls, wc, &circle, &ruling(x0)),
            Err(RejectReason::ObliqueCylinderCut),
            "not parted, and an ellipse is not measured by a circle's strip — the oblique coverage's end"
        );
    }
}
