//! The chain fold against hand-worked frames and against the exact-node map it wraps.

use super::*;

fn r(n: i128, d: i128) -> Rat {
    Rat::new(n, d).expect("a rational")
}

fn z(n: i128) -> Rat {
    Rat::from_int(n)
}

/// The datum an application states with normal `(2, 3, 6)` (length `7`) through `o`, under the
/// derived placement: `û = (−3, 2, 0)/√13`, `v̂ = n × u / (7√13)`, `ŵ = (2, 3, 6)/7`.
fn seven(o: [Rat; 3]) -> ChainFold {
    let d = z(0)
        .checked_sub(z(2).checked_mul(o[0]).unwrap())
        .unwrap()
        .checked_sub(z(3).checked_mul(o[1]).unwrap())
        .unwrap()
        .checked_sub(z(6).checked_mul(o[2]).unwrap())
        .unwrap();
    let pf = plane_frame_placed([z(2), z(3), z(6), d], None, false).expect("a frame");
    assert_eq!(pf.u_raw, [z(-3), z(2), z(0)], "premise: the derived +u");
    ChainFold::frame(&pf)
}

/// A frame whose lengths are all rational — the plane `z = 2`, named at `(1, −1, 2)` with
/// `+u = (3, 4, 0)`: `v = n × u = (−4, 3, 0)`, lengths `5, 5, 1` — so it carries every rational
/// statement rationally.
fn five() -> ChainFold {
    let pf = plane_frame_placed(
        [z(0), z(0), z(1), z(-2)],
        Some(([z(1), z(-1), z(2)], [z(3), z(4), z(0)])),
        false,
    )
    .expect("a frame");
    assert_eq!(pf.u_raw, [z(3), z(4), z(0)], "premise: +u as named");
    ChainFold::frame(&pf)
}

/// **Without a frame the fold is the exact-node map**, carry for carry.
#[test]
fn a_fold_of_exact_nodes_is_its_axis_map() {
    let a = AxisAffine::rotation(Rotation {
        axis: Axis::X,
        pivot: [z(1), z(0), r(1, 2)],
        angle: Angle::from_deg(z(90)).unwrap(),
    })
    .unwrap()
    .then(&AxisAffine::mirror(Axis::Y, r(3, 4)))
    .then(&AxisAffine::translation([r(1, 3), z(-2), z(5)]));
    let f = ChainFold::axis(a);
    let p = [r(5, 7), z(-1), r(2, 3)];
    assert_eq!(f.point_rat(p).ok(), a.point_rat(p));
    assert_eq!(f.dir_rat(p).ok(), a.dir_rat(p));
    assert_eq!(f.dir_along(p).ok(), a.dir_rat(p));
    assert_eq!(
        f.plane_coeffs([z(1), z(2), z(0), z(-3)]).ok(),
        a.plane_coeffs([z(1), z(2), z(0), z(-3)])
    );
    assert_eq!(f.point_rat_inverse(p).ok(), a.inverse().point_rat(p));
    assert_eq!(f.det(), a.det());
}

/// **A frame with one rational length carries exactly what lies along it.**
#[test]
fn a_frame_carries_the_statements_its_lengths_allow() {
    let o = [z(1), z(2), z(3)];
    let f = seven(o);
    // The origin of the frame is its projection of the world origin, plus nothing.
    let origin = f.point_rat([z(0); 3]).expect("the frame's origin");
    // Four along `ŵ` is `4·(2, 3, 6)/7` from there.
    let up = f.point_rat([z(0), z(0), z(4)]).expect("a point along ŵ");
    let want: [Rat; 3] = core::array::from_fn(|k| {
        origin[k]
            .checked_add(r(4 * [2, 3, 6][k], 7))
            .expect("small")
    });
    assert_eq!(up, want);
    assert_eq!(
        f.point_rat([z(1), z(0), z(0)]),
        Err(Carry::Irrational),
        "a step of one along û is u/√13"
    );
    assert_eq!(f.dir_rat([z(1), z(0), z(0)]), Err(Carry::Irrational));
    assert_eq!(f.dir_rat([z(0), z(0), z(7)]), Ok([z(2), z(3), z(6)]));
    assert_eq!(
        f.dir_along([z(1), z(0), z(0)]),
        Ok([z(-3), z(2), z(0)]),
        "along û"
    );
    // `|v|²/|u|² = 637/13 = 49`: `û` and `v̂` are in rational ratio, and `û + v̂` runs along
    // `7·u + v = 7·(−3, 2, 0) + (−12, −18, 13)`.
    assert_eq!(
        f.dir_along([z(1), z(1), z(0)]),
        Ok([z(-33), z(-4), z(13)]),
        "û + v̂"
    );
    assert_eq!(
        f.dir_along_inverse([z(-33), z(-4), z(13)]),
        Ok([z(1), z(1), z(0)]),
        "back along û + v̂"
    );
    assert_eq!(
        f.dir_along([z(1), z(0), z(1)]),
        Err(Carry::Irrational),
        "|u|²/|n|² = 13/49 is not a square"
    );
    // The far cap `w = 4`, across `ŵ`.
    let cap = f
        .plane_coeffs([z(0), z(0), z(1), z(-4)])
        .expect("the far cap");
    let dot = |c: &[Rat; 4], p: [Rat; 3]| {
        (0..3).fold(c[3], |s, k| {
            s.checked_add(c[k].checked_mul(p[k]).unwrap()).unwrap()
        })
    };
    assert_eq!(dot(&cap, up), z(0), "the far cap holds the point four up");
    assert_eq!(
        f.plane_coeffs([z(1), z(0), z(0), z(-1)]),
        Err(Carry::Irrational),
        "the plane u = 1 lies one along û, at u/√13, off the origin"
    );
    let wall = f.plane_coeffs([z(1), z(0), z(0), z(0)]).expect("u = 0");
    assert_eq!(dot(&wall, origin), z(0));
    assert_eq!(dot(&wall, up), z(0));
    // Back: the point four up is `(0, 0, 4)` in the frame, a point off `ŵ` is not rational there.
    assert_eq!(f.point_rat_inverse(up), Ok([z(0), z(0), z(4)]));
    // The derived origin is the world origin's projection, so the world origin is on `ŵ`'s line;
    // a step along `x` from it is not.
    assert!(f.point_rat_inverse([z(0), z(0), z(0)]).is_ok());
    assert_eq!(
        f.point_rat_inverse([z(1), z(0), z(0)]),
        Err(Carry::Irrational)
    );
    let back = f.dir_along_inverse([z(-3), z(2), z(0)]).expect("û's line");
    assert!(
        back[0] > z(0) && back[1] == z(0) && back[2] == z(0),
        "a positive multiple of +u: {back:?}"
    );
    assert_eq!(f.det(), 1);
}

/// **Composition: an exact node after, a rational frame inside, and what does not compose.**
#[test]
fn a_frame_composes_where_its_columns_carry_along() {
    let f = seven([z(0); 3]);
    // A translation after the frame moves every carried point by it.
    let t = AxisAffine::translation([r(3, 2), z(0), z(-1)]);
    let moved = f.then_axis(&t);
    let up = f.point_rat([z(0), z(0), z(1)]).unwrap();
    assert_eq!(
        moved.point_rat([z(0), z(0), z(1)]),
        Ok(t.point_rat(up).unwrap())
    );
    // A mirror after it reverses the handedness.
    assert_eq!(f.then_axis(&AxisAffine::mirror(Axis::Z, z(0))).det(), -1);
    // An axis-aligned frame inside it (the plane `z = 2` there, derived placement): every column
    // is one of the outer frame's axes, so the composite folds and carries as the outer one does
    // two up.
    let inner =
        ChainFold::frame(&plane_frame_placed([z(0), z(0), z(1), z(-2)], None, false).unwrap());
    let both = inner.then(&f).expect("an axis-aligned frame composes");
    assert_eq!(
        both.point_rat([z(0), z(0), z(1)]),
        f.point_rat([z(0), z(0), z(3)])
    );
    assert_eq!(
        both.dir_along([z(1), z(0), z(0)]),
        f.dir_along([z(1), z(0), z(0)])
    );
    // A frame inside whose origin the outer frame carries off the rationals does not compose.
    let off =
        ChainFold::frame(&plane_frame_placed([z(1), z(0), z(0), z(-1)], None, false).unwrap());
    assert!(
        off.then(&f).is_none(),
        "the plane x = 1 has its origin one along û — u/√13 — from the datum's"
    );
}

/// **Width is asked after geometry.** A point whose image overflows `Rat` is `Width`; one whose
/// image is irrational is `Irrational`, however wide.
#[test]
fn a_carry_names_geometry_before_width() {
    let near_max = z(i128::MAX - 1);
    let f = seven([z(0); 3]).then_axis(&AxisAffine::translation([z(0), z(0), near_max]));
    assert_eq!(
        f.point_rat([z(0), z(0), z(7)]),
        Err(Carry::Width),
        "6 past near-max"
    );
    assert_eq!(f.point_rat([z(1), z(0), z(7)]), Err(Carry::Irrational));
}

/// **An exact node before a frame folds into the frame's columns** — `F ∘ A` with `A` a quarter
/// turn about `z` and a step `7` up: a point, a direction and the handedness come out as the
/// frame alone gives them for `A`'s image, worked by hand.
#[test]
fn an_exact_node_inside_a_frame_folds_into_its_columns() {
    let f = seven([z(1), z(2), z(3)]);
    let quarter = AxisAffine::rotation(Rotation {
        axis: Axis::Z,
        pivot: [z(0); 3],
        angle: Angle::from_deg(z(90)).unwrap(),
    })
    .unwrap()
    .then(&AxisAffine::translation([z(0), z(0), z(7)]));
    let both = ChainFold::axis(quarter)
        .then(&f)
        .expect("the frame composes");
    // `(0, 0, 1)` turns to itself and steps to `(0, 0, 8)`.
    assert_eq!(
        both.point_rat([z(0), z(0), z(1)]),
        f.point_rat([z(0), z(0), z(8)])
    );
    // `x̂` turns to `ŷ`, which the frame carries along `v̂`.
    assert_eq!(
        both.dir_along([z(1), z(0), z(0)]),
        f.dir_along([z(0), z(1), z(0)])
    );
    assert_eq!(both.point_rat([z(1), z(0), z(0)]), Err(Carry::Irrational));
    assert_eq!(both.det(), 1);
    let mirrored = ChainFold::axis(AxisAffine::mirror(Axis::X, z(0)))
        .then(&f)
        .expect("a mirror composes");
    assert_eq!(mirrored.det(), -1);
    assert_eq!(
        mirrored.dir_along([z(1), z(0), z(0)]),
        f.dir_along([z(-1), z(0), z(0)])
    );
}

/// **Back through the frame, geometry is asked before width too.** The preimage of `(1, 0, −5)`
/// under the frame and a near-max step is irrational — its `v̂` coordinate divides by `7√13` —
/// and the step's size does not hide that.
#[test]
fn a_carry_back_names_geometry_before_width() {
    let f = seven([z(0); 3]).then_axis(&AxisAffine::translation([z(0), z(0), z(i128::MAX - 1)]));
    assert_eq!(
        f.point_rat_inverse([z(1), z(0), z(-5)]),
        Err(Carry::Irrational)
    );
    let small = seven([z(0); 3]).then_axis(&AxisAffine::translation([z(0), z(0), z(7)]));
    assert_eq!(
        small.point_rat_inverse([z(1), z(0), z(2)]),
        Err(Carry::Irrational)
    );
}

/// **An offset lost to overflow leaves the geometry asked.** Two near-max steps compose to an
/// offset `Rat` cannot hold: a statement the frame puts off the rationals is still `Irrational`,
/// one it keeps rational is `Width`, and directions — which do not read the offset — answer. The
/// same holds where the lost offset is inside the frame.
#[test]
fn a_lost_offset_still_asks_the_geometry() {
    let big = AxisAffine::translation([z(0), z(0), z(i128::MAX - 1)]);
    let lost = big.then(&big);
    let f = seven([z(0); 3]);
    for (what, fold) in [
        ("after the frame", f.then_axis(&lost)),
        (
            "inside the frame",
            ChainFold::axis(lost).then(&f).expect("composes"),
        ),
    ] {
        assert_eq!(
            fold.point_rat([z(0), z(0), z(7)]),
            Err(Carry::Width),
            "{what}"
        );
        assert_eq!(
            fold.point_rat([z(1), z(0), z(0)]),
            Err(Carry::Irrational),
            "{what}"
        );
        assert_eq!(
            fold.plane_coeffs([z(1), z(0), z(0), z(-1)]),
            Err(Carry::Irrational),
            "{what}: u = 1"
        );
        assert_eq!(
            fold.plane_coeffs([z(0), z(0), z(1), z(-4)]),
            Err(Carry::Width),
            "{what}: w = 4"
        );
        assert_eq!(
            fold.dir_along([z(1), z(0), z(0)]),
            f.dir_along([z(1), z(0), z(0)]),
            "{what}"
        );
    }
}

/// **The zero direction is its own image** — through a frame as without one.
#[test]
fn the_zero_direction_carries_to_zero() {
    let f = seven([z(1), z(2), z(3)]);
    assert_eq!(f.dir_along([z(0); 3]), Ok([z(0); 3]));
    assert_eq!(f.dir_along_inverse([z(0); 3]), Ok([z(0); 3]));
    assert_eq!(f.dir_rat([z(0); 3]), Ok([z(0); 3]));
}

/// **A carried plane's name and its facing come from one carry.** The turn
/// [`ChainFold::plane_name`] reports is checked against points, not against another carry of the
/// normal: a point on the plane's positive side lands where the image name, times the turn, is
/// positive. Through `seven` only points along `ŵ` carry rationally, so there the planes are
/// across `ŵ`; through a frame whose lengths are all rational, and without a frame, any plane.
#[test]
fn a_planes_name_and_facing_come_from_one_carry() {
    let f = seven([z(1), z(2), z(3)]);
    let g = five();
    let across = vec![[z(0), z(0), z(1), z(-4)], [z(0), z(0), z(1), z(4)]];
    let any = vec![
        [z(0), z(0), z(1), z(-4)],
        [z(1), z(0), z(0), z(0)],
        [z(1), z(1), z(0), z(-3)],
        [z(0), z(1), z(-2), z(5)],
    ];
    let mirror = AxisAffine::mirror(Axis::X, z(1));
    let folds = [
        ("seven", f.clone(), &across),
        ("seven, mirrored after", f.then_axis(&mirror), &across),
        ("five", g.clone(), &any),
        ("five, mirrored after", g.then_axis(&mirror), &any),
        (
            "no frame",
            ChainFold::axis(
                AxisAffine::mirror(Axis::Y, r(1, 2)).then(&AxisAffine::translation([
                    z(3),
                    z(-1),
                    r(1, 3),
                ])),
            ),
            &any,
        ),
    ];
    let value = |c: &[Rat; 4], p: &[Rat; 3]| {
        (0..3).fold(c[3], |s, k| {
            s.checked_add(c[k].checked_mul(p[k]).unwrap()).unwrap()
        })
    };
    let points = [
        [z(0), z(0), z(5)],
        [z(0), z(0), z(-5)],
        [z(2), z(1), z(0)],
        [z(-2), z(3), z(1)],
    ];
    let mut turns = [0, 0];
    for (what, fold, planes) in &folds {
        for c in planes.iter() {
            let name = canonical_plane_coeffs(*c).unwrap();
            let (image, turn) = fold
                .plane_name(&PlaneName::Narrow(name))
                .unwrap_or_else(|e| panic!("{what} {c:?}: {e:?}"));
            let image = *image.narrow().expect("narrow");
            assert_eq!(fold.plane_coeffs(name), Ok(image), "{what} {c:?}");
            let mut asked = 0;
            for p in &points {
                let side = value(&name, p);
                let Ok(q) = fold.point_rat(*p) else { continue };
                if side == z(0) {
                    continue;
                }
                let there = value(&image, &q);
                assert_eq!(
                    (there > z(0)) == (side > z(0)),
                    turn == 1,
                    "{what} {c:?} at {p:?}"
                );
                asked += 1;
            }
            assert!(asked > 0, "{what} {c:?}: no point asked");
            turns[usize::from(turn == 1)] += 1;
        }
    }
    assert!(
        turns[0] > 0 && turns[1] > 0,
        "both facings asked: {turns:?}"
    );
}

/// **A wide name is carried at any width, its geometry asked first.** `x = 2¹³⁰` and `z = 2¹³⁰`
/// as frame statements: across `ŵ` the image is rational (and wide), off it the frame's `√13`
/// makes it irrational; without a frame a wide name stays wide.
#[test]
fn a_wide_name_is_carried_and_its_geometry_asked() {
    use num_bigint::BigInt;
    let far: BigInt = -(BigInt::from(1) << 130_u32);
    let one = BigInt::from(1);
    let zero = BigInt::from(0);
    let across = PlaneName::Wide([zero.clone(), zero.clone(), one.clone(), far.clone()]);
    let off = PlaneName::Wide([one.clone(), zero.clone(), zero.clone(), far.clone()]);
    let f = seven([z(1), z(2), z(3)]);
    // `w = 2¹³⁰` is `n·(x − t) = 7·2¹³⁰` with `n·t = n·(1, 2, 3) = 26`.
    let seven_far = &far * 7 - 26;
    assert_eq!(
        f.plane_name(&across),
        Ok((
            PlaneName::Wide([BigInt::from(2), BigInt::from(3), BigInt::from(6), seven_far]),
            1
        ))
    );
    assert_eq!(f.plane_name(&off), Err(Carry::Irrational));
    let t = ChainFold::axis(AxisAffine::translation([z(5), z(0), z(0)]));
    assert_eq!(
        t.plane_name(&off),
        Ok((PlaneName::Wide([one, zero.clone(), zero, far - 5]), 1))
    );
}

/// **A frame inside a frame, an exact node between** — `G ∘ A ∘ F` with `A` a quarter turn, a
/// mirror and a step (`det A = −1`) and `G` a frame whose lengths are all rational, so `G` carries
/// whatever `A ∘ F` gives it: the composite answers each statement as the two carried in turn do,
/// the turns and the handedness multiplied.
#[test]
fn a_frame_inside_a_frame_carries_as_the_two_in_turn() {
    let a = AxisAffine::rotation(Rotation {
        axis: Axis::Z,
        pivot: [z(0); 3],
        angle: Angle::from_deg(z(90)).unwrap(),
    })
    .unwrap()
    .then(&AxisAffine::mirror(Axis::Y, z(1)))
    .then(&AxisAffine::translation([z(2), r(1, 3), z(-5)]));
    assert_eq!(a.det(), -1, "premise");
    let inner = seven([z(1), z(2), z(3)]).then_axis(&a);
    let outer = five();
    let both = inner.then(&outer).expect("a rational outer frame composes");
    assert_eq!(both.det(), inner.det() * outer.det());
    for p in [
        [z(0), z(0), z(0)],
        [z(0), z(0), z(7)],
        [z(0), z(0), r(-3, 2)],
    ] {
        let mid = inner.point_rat(p).unwrap();
        assert_eq!(both.point_rat(p), outer.point_rat(mid), "{p:?}");
    }
    for d in [[z(1), z(0), z(0)], [z(1), z(1), z(0)], [z(0), z(0), z(1)]] {
        let mid = inner.dir_along(d).unwrap();
        assert_eq!(both.dir_along(d), outer.dir_along(mid), "{d:?}");
    }
    for c in [
        [z(0), z(0), z(1), z(-4)],
        [z(1), z(0), z(0), z(0)],
        [z(1), z(1), z(0), z(0)],
    ] {
        let name = PlaneName::Narrow(canonical_plane_coeffs(c).unwrap());
        let (mid, t1) = inner.plane_name(&name).unwrap();
        let (out, t2) = outer.plane_name(&mid).unwrap();
        assert_eq!(both.plane_name(&name), Ok((out, t1 * t2)), "{c:?}");
    }
}

/// **A cylinder is one statement of three carries, and its geometry is asked first.** On `ŵ`
/// through `seven`'s origin: a seam reference along `û` carries, one mixing `û` and `ŵ` does not
/// (`√(13/49)`). With the offset lost to overflow the origin cannot be read, yet the mixed
/// reference still names the cylinder `Irrational`.
#[test]
fn a_cylinder_carries_its_geometry_before_its_width() {
    let cylinder = |reference: [Rat; 3]| {
        CylinderDef::new(
            [z(0); 3],
            [z(0), z(0), z(1)],
            reference,
            BigRat::square_of(z(1)),
        )
        .unwrap()
    };
    let (along_u, mixed) = (cylinder([z(1), z(0), z(0)]), cylinder([z(1), z(0), z(1)]));
    let f = seven([z(0); 3]);
    let carried = f.cylinder(&along_u).expect("rational");
    assert_eq!(carried.origin(), f.point_rat([z(0); 3]).unwrap());
    assert_eq!(carried.dir(), [z(2), z(3), z(6)]);
    assert_eq!(carried.ref_dir(), [z(-3), z(2), z(0)]);
    assert_eq!(f.cylinder(&mixed), Err(Carry::Irrational));
    let big = AxisAffine::translation([z(0), z(0), z(i128::MAX - 1)]);
    let lost = f.then_axis(&big.then(&big));
    assert_eq!(lost.cylinder(&along_u), Err(Carry::Width));
    assert_eq!(lost.cylinder(&mixed), Err(Carry::Irrational));
}
