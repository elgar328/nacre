//! **The bit census — a baseline to diff across a change that moves coordinates.**
//!
//! Prints one line per case: the exact bits of volume/area/centroid **and a hash of the result's
//! sorted vertex coordinates**. The derived quantities alone are not enough: a 1-ULP coordinate
//! change can cancel out of a volume integral, so a census of volumes can pass while the geometry
//! moved. The coordinate hash is what makes "nothing changed" checkable.
//!
//! Not an assertion — a dump. Run it on two commits and `diff`:
//!
//! ★ **Re-baselined once at S9** (world-plane seeding): the seeds intern with every
//! origin-touching producer, so the survivor's f64 plane cache — and with it the `in:` plane
//! digest — moved on the origin-touching lines. The switch was gated by an ε-equivalence check
//! (topology exact, volumes/areas/centroids within 2⁻⁴⁰ relative, vertex hashes byte-equal) and
//! the measured deviations are in `docs/dev-log.md`. `stat seeded_hits` below is the
//! falsifiability bridge for that population.
//!
//! ```text
//! cargo test -p nacre-ops --release --test census -- --ignored --nocapture | grep '^c '
//! ```
//!
//! ★★★ **Diff it across *profiles* too, not only across commits.** Drop `--release` and the same
//! 148 lines must come out — they do, measured. That is not a formality: `Angle`'s f64 route is
//! `(deg.to_f64() * PI / 180.0).cos()`, and LLVM evaluates that at compile time wherever it can see
//! the angle, one ulp away from what libm returns at run time. Two builds disagreeing here would
//! mean the coordinates a model stores depend on how it was compiled.
//!
//! ★★ **What keeps them equal is that an angle crosses the model store** before anything realizes
//! it — written into an `Operation`, pushed, read back — and no optimiser propagates a constant
//! through a heap structure. So this is a real check with a real way to fail, and no in-process
//! test can express it: a test runs in one profile. It belongs here, next to the diff it extends.
//!
//! ★★★★★ **What a census can and cannot see.** Every case is built from coordinates written here,
//! so it sees a change that moves a coordinate — and it is **blind by construction** to a change in
//! a population it does not contain. The first 130 lines are all *short* decimals, which is why a
//! regression that closed the exact path for every arbitrarily-tilted plane crossed this file
//! bit-identical, twice. The `fw` (full-width coordinates) and `tp` (tilted sketch plane) families
//! exist for that reason. **Read "bit-identical" as evidence only about the population present.**

use nacre_math::Point3;
use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
use nacre_ops::{DatumDef, SketchFrame, SketchPlane};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

fn xf(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    )
    .expect("transform") else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();
    solid
}

fn mirror(m: &mut Model, s: Handle<Solid>, axis: Axis, offset: Rat) -> Handle<Solid> {
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid: s,
            axis,
            offset,
        },
    )
    .expect("mirror") else {
        unreachable!("mirror yields Mirror output")
    };
    m.rebuild_adjacency();
    solid
}

/// A hash of the solid's **plane coefficients**, sorted by bits — the other half of what a
/// boolean actually reads.
///
/// **Measured, not assumed:** a boolean's result vertices are all recomputed from plane triples
/// (a measured `VertexDef::ThreePlane`), so the result carries *no* trace of the operands' vertex
/// coordinates. A change that moves operand vertices but leaves the planes alone is therefore
/// invisible in the result — which is exactly what happened the first time this census was used.
/// Recording the operands is what makes the census see the change it exists to see.
fn plane_digest(m: &Model, s: Handle<Solid>) -> (usize, u64) {
    use std::hash::{Hash, Hasher};
    let mut bits: Vec<[u64; 4]> = Vec::new();
    let src = m.solids.get(s).clone();
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            match m.surface(m.faces.get(fh).surface) {
                nacre_geom::Surface::Plane(pl) => {
                    bits.push(pl.coefficients().map(f64::to_bits));
                }
                nacre_geom::Surface::Cylinder(_) => {}
            }
        }
    }
    bits.sort_unstable();
    bits.dedup();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bits.hash(&mut h);
    (bits.len(), h.finish())
}

/// The operands' own state, before the boolean consumes them: vertices *and* planes.
fn operands(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> String {
    let f = |s| {
        let (vn, vh) = coord_digest(m, s);
        let (pn, ph) = plane_digest(m, s);
        format!("v{vn}h{vh:016x}p{pn}h{ph:016x}")
    };
    format!("{}+{}", f(a), f(b))
}

/// A hash of the solid's vertex coordinates, **sorted by bits** so it does not depend on traversal
/// or handle order — only on the set of coordinates the operation produced.
fn coord_digest(m: &Model, s: Handle<Solid>) -> (usize, u64) {
    use std::hash::{Hash, Hasher};
    let mut bits: Vec<[u64; 3]> = Vec::new();
    let src = m.solids.get(s).clone();
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            let face = m.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in m.edges.get(he.edge).vertices.iter() {
                        let p = m.vertex_point(vh).as_array();
                        bits.push([p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]);
                    }
                }
            }
        }
    }
    bits.sort_unstable();
    bits.dedup();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bits.hash(&mut h);
    (bits.len(), h.finish())
}

fn record(
    tag: &str,
    m: &Model,
    inputs: &str,
    out: &Result<Vec<Handle<Solid>>, nacre_ops::BoolError>,
) {
    print!("c {tag} in:{inputs} ");
    match out {
        Err(e) => println!("ERR {e:?}"),
        Ok(v) if v.is_empty() => println!("EMPTY"),
        Ok(v) => {
            // Every body, in the order the engine returned them — the count is part of the answer.
            let mut parts = Vec::new();
            for &s in v {
                let p = nacre_props::mass_props(m, s).expect("props");
                // ★ **A curved result has no centroid to print, and that is a record, not a
                // crash.** `centroid` refuses a face it cannot integrate a first moment over
                // (`CentroidOfCurvedFace`), which the census met the moment cylinder booleans
                // started returning solids. Volume and area still measure the whole body, so the
                // row keeps its comparable numbers and marks the missing field rather than
                // taking the process down.
                let (n, d) = coord_digest(m, s);
                let centroid = match nacre_props::centroid(m, s) {
                    Ok(c) => format!(
                        "{:016x},{:016x},{:016x}",
                        c[0].to_bits(),
                        c[1].to_bits(),
                        c[2].to_bits()
                    ),
                    Err(nacre_props::PropsError::CentroidOfCurvedFace) => "curved".to_string(),
                    Err(e) => panic!("centroid: {e:?}"),
                };
                parts.push(format!(
                    "{:016x}/{:016x}/{centroid}/v{n}h{d:016x}",
                    p.volume.to_bits(),
                    p.area.to_bits(),
                ));
            }
            println!("{} {}", v.len(), parts.join(" "));
        }
    }
}

const KINDS: [(&str, BoolKind); 3] = [
    ("fuse", BoolKind::Fuse),
    ("cut", BoolKind::Cut),
    ("common", BoolKind::Common),
];

#[test]
#[ignore = "census dump, not an assertion (run with --ignored --nocapture)"]
fn dump() {
    // ── Axis-aligned, no motion at all: the untouched baseline.
    let boxes: [([f64; 3], [f64; 3]); 6] = [
        ([0.5, 0.5, 0.5], [2.5, 2.5, 2.5]),
        ([0.25, 0.25, 0.25], [0.75, 0.75, 0.75]),
        ([1.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
        ([0.3, -1.0, 0.3], [0.7, 2.0, 0.7]),
        ([-1.0, 0.2, 0.2], [0.5, 0.8, 0.8]),
        ([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]),
    ];
    for (kn, k) in KINDS {
        for (i, (lo, hi)) in boxes.iter().enumerate() {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
            let b = m.add_cuboid(Point3::from_array(*lo), Point3::from_array(*hi));
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("aa {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── Rotated: one operand turned about each axis by several angles.
    for ax in [Axis::X, Axis::Y, Axis::Z] {
        for deg in [7i128, 30, 45, 90, 123] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([2.0, 2.0, 2.0]),
                );
                let b = m.add_cuboid(
                    Point3::from_array([0.4, 0.4, -1.0]),
                    Point3::from_array([1.6, 1.6, 3.0]),
                );
                m.rebuild_adjacency();
                let b = xf(
                    &mut m,
                    b,
                    Isometry::rotation(Rotation {
                        axis: ax,
                        point: [Rat::from_int(1); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                    }),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("rot {ax:?} {deg} {kn}"), &m, &inputs, &out);
            }
        }
    }

    // ── **Translated — what the previous census was missing entirely.**
    // Dyadic offsets (expected untouched by the motion work) and non-dyadic ones (expected to be
    // where every difference lands), on both operands so the *relative* placement varies too.
    let offsets: [(i128, i128); 10] = [
        (1, 1),
        (3, 1),
        (1, 2),
        (-5, 4),
        (7, 8),
        (1, 3),
        (2, 5),
        (3, 10),
        (7, 11),
        (13, 23),
    ];
    for (n, d) in offsets {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cuboid(
                Point3::from_array([1.0, 1.0, 1.0]),
                Point3::from_array([3.0, 3.0, 3.0]),
            );
            m.rebuild_adjacency();
            let t = Rat::new(n, d).expect("offset");
            let a = xf(
                &mut m,
                a,
                Isometry::translation([t, Rat::from_int(0), Rat::from_int(0)]),
            );
            let b = xf(&mut m, b, Isometry::translation([t, t, Rat::from_int(0)]));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("tr {n}/{d} {kn}"), &m, &inputs, &out);
        }
        // Chained: the same total move split in two, which f64 accumulation and exact folding
        // disagree about.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        m.rebuild_adjacency();
        let t = Rat::new(n, d).expect("offset");
        let mut a = a;
        for _ in 0..2 {
            a = xf(
                &mut m,
                a,
                Isometry::translation([t, Rat::from_int(0), Rat::from_int(0)]),
            );
        }
        let b = xf(
            &mut m,
            b,
            Isometry::translation([
                t.checked_add(t).expect("double"),
                Rat::from_int(0),
                Rat::from_int(0),
            ]),
        );
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, BoolKind::Fuse, a, b);
        m.rebuild_adjacency();
        record(&format!("tr2 {n}/{d} fuse"), &m, &inputs, &out);
    }

    // ── Mirrored, about dyadic and non-dyadic planes.
    for (n, d) in [(0i128, 1i128), (1, 2), (1, 3), (7, 22), (5, 7)] {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([1.0, 0.0, 0.0]),
                Point3::from_array([2.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([-1.0, 0.25, 0.25]),
                Point3::from_array([1.5, 0.75, 0.75]),
            );
            m.rebuild_adjacency();
            let a = mirror(&mut m, a, Axis::X, Rat::new(n, d).expect("plane"));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("mir {n}/{d} {kn}"), &m, &inputs, &out);
        }
    }

    // ── Rotate then place: refused today, expected to build after the motion work.
    for deg in [7i128, 30, 45] {
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([-2.0, -2.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        let tool = m.add_cuboid(
            Point3::from_array([-0.5, -0.5, -1.0]),
            Point3::from_array([0.5, 0.5, 2.0]),
        );
        m.rebuild_adjacency();
        let tool = xf(
            &mut m,
            tool,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
            }),
        );
        let place = Isometry::translation([Rat::from_int(5), Rat::from_int(-3), Rat::from_int(2)]);
        let tool = xf(&mut m, tool, place);
        let base = xf(&mut m, base, place);
        let inputs = operands(&m, base, tool);
        let out = boolean(&mut m, BoolKind::Cut, base, tool);
        m.rebuild_adjacency();
        record(&format!("rt {deg}"), &m, &inputs, &out);
    }

    // ── **Constructed by `extrude`, with the dimension split across steps.**
    //
    // Every family above builds its operands with `add_cuboid`, which takes literal corners and
    // does no arithmetic — so none of them can see a change to how construction *accumulates*.
    // `add_cuboid` is also `#[cfg(feature = "test-util")]`, i.e. the census was measuring a path
    // production does not use: the playground and the kit both go through `extrude`.
    //
    // Here one operand is raised in a single step and the other in two that should add to the
    // same height. In `f64` they do not (`1.1 + 6.6 != 7.7`), and the difference lands in the
    // operand digests below — which is the point of the family.
    for (kn, k) in KINDS {
        for (i, (whole, first, second)) in [(7.7, 1.1, 6.6), (3.0, 0.1, 2.9), (1.0, 0.3, 0.7)]
            .iter()
            .enumerate()
        {
            let mut m = Model::new();
            let a = ex(&mut m, 0.0, [0.0, 0.0], [2.0, 2.0], *whole);
            let b1 = ex(&mut m, 0.0, [1.0, 1.0], [3.0, 3.0], *first);
            let b2 = ex(&mut m, *first, [1.0, 1.0], [3.0, 3.0], *second);
            let b = boolean(&mut m, BoolKind::Fuse, b1, b2).expect("stack fuses")[0];
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("ex {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── ★★★★★ **Full-width coordinates — the population every family above is blind to.**
    //
    // Every corner and dimension above is a *short* decimal (`0.5`, `2.5`, `7.7`), and every
    // rational derived from one stays far inside `i128`. So a change to what the kernel can and
    // cannot **name** leaves all 130 lines bit-identical, and this census reports nothing.
    //
    // ★★★★★ That is not hypothetical: a regression that closed the exact path for every
    // arbitrarily-tilted plane shipped, and the census was bit-identical across it — **twice**,
    // because the same blindness had already been recorded once.
    //
    // ★★ **Measured, and it picks the right target.** A cuboid on seventeen-digit corners
    // overflows `plane_through_points`' *intermediates* on two of its six faces — the
    // narrow-route population the wide derivation (`plane_name_big`) names, since their
    // canonical answers are small. (An earlier reading here claimed those two faces went
    // *unnamed*; that stopped being true when the wide derivation landed — the genuinely
    // wide-vessel population lives in computed ring coordinates and is counted by
    // `WIDE_PLANES`.) The proptests are where this population lives today
    // (`stacked_boxes_merge_volumes` and friends generate arbitrary `f64`), and a proptest
    // cannot be a census line: it has no fixed coordinates to diff. These constants are those
    // coordinates, pinned.
    let fw: [([f64; 3], [f64; 3]); 3] = [
        (
            [-2.8374652839472, 1.0937465283947, -0.5837465283947],
            [1.4738264859372, 2.9384756293847, 0.8473625849372],
        ),
        // A stacked pair sharing one interface plane — the coplanar-contact route, on coordinates
        // whose interface plane cannot be named today.
        (
            [-2.8374652839472, 1.0937465283947, 0.8473625849372],
            [1.4738264859372, 2.9384756293847, 2.1937465283947],
        ),
        // Overlapping, so the result's planes come from both operands.
        (
            [-1.1937465283947, 1.9384756293847, -0.1837465283947],
            [2.8473625849372, 3.4738264859372, 1.4937465283947],
        ),
    ];
    for (kn, k) in KINDS {
        for (i, (lo, hi)) in fw.iter().enumerate() {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array(fw[0].0), Point3::from_array(fw[0].1));
            let b = m.add_cuboid(Point3::from_array(*lo), Point3::from_array(*hi));
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("fw {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── **A prism on a tilted sketch plane** — the `Motion::Frame` path, which no family above
    // reaches either. Its ring lives in the plane's own frame, so its walls *are* nameable; what
    // it exercises is the frame machinery, not the width limit.
    for (i, n) in [
        [0.3141592653589793, -0.2718281828459045, 1.0],
        [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
        [0.1, 0.2, 0.30000000000000004],
    ]
    .iter()
    .enumerate()
    {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = tilted_prism(&mut m, *n, 0.0, 1.0, 0.75);
            // A second prism on the *same* tilted plane, offset in the sketch and raised further:
            // its base cap and `a`'s are one plane, and its walls can be coplanar with `a`'s.
            let b = tilted_prism(&mut m, *n, 0.5, 1.6, 1.25);
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("tp {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── **A framed sketch on an `n·n`-overflow wall** (`wf`) — the population S4 opens, which
    // no family above contains: `fw` is all-narrow (measured at S2) and `tp`'s walls are
    // in-frame narrow. A prism on a fully tilted, exactly-orthonormal *decimal* frame has walls
    // whose names run ~110 bits — narrow, with squared lengths past `i128` — and a pad or
    // pocket on such a wall used to fall to the f64 path. Without this family a regression in
    // the wide-frame road would cross this file bit-identical (the 8b lesson, a third time).
    {
        use nacre_math::{Point2, Vector3};
        let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
        for (kn, pad) in [("pocket", false), ("pad", true)] {
            let mut m = Model::new();
            let plane = nacre_ops::SketchPlane::from_axes(
                Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                Vector3::from_array([0.6, 0.8, 0.0]),
                Vector3::from_array([-0.48, 0.36, 0.8]),
            );
            let __f105 = datum_frame(&mut m, plane);
            let OpOutput::Extrude { solid, .. } = apply(
                &mut m,
                &Operation::Extrude {
                    frame: __f105,
                    profile: nacre_ops::Profile2d::polygon(vec![
                        p2(0.1111111111111111, 0.1234567890123456),
                        p2(4.123456789012345, 0.2345678901234567),
                        p2(3.9876543210987654, 3.1234567890123459),
                        p2(0.2222222222222222, 2.765432109876543),
                    ])
                    .unwrap(),
                    dist: 2.5,
                },
            )
            .expect("the wf base prism") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            // The family qualifies itself: a wall whose name is narrow but whose squared
            // lengths overflow — the exact population the narrow frame derivation declines.
            let wall = *m
                .shells
                .get(m.solids.get(solid).outer)
                .faces
                .iter()
                .find(|&&f| {
                    let s = m.faces.get(f).surface;
                    m.surface_name
                        .get(&s)
                        .and_then(|n| n.narrow())
                        .is_some_and(|c| nacre_scalar::plane_frame_default(*c).is_none())
                })
                .expect("the wf population vanished — retune the constants");
            // A small square centred on the wall, in its own sketch frame.
            let sp = nacre_ops::face_plane(&m, wall).expect("planar");
            let d = nacre_props::face_props(&m, wall).unwrap().centroid - sp.origin();
            let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
            let profile = nacre_ops::Profile2d::polygon(vec![
                p2(cu - 0.3, cv - 0.3),
                p2(cu + 0.3, cv - 0.3),
                p2(cu + 0.3, cv + 0.3),
                p2(cu - 0.3, cv + 0.3),
            ])
            .unwrap();
            let inputs = operands(&m, solid, solid);
            let op = if pad {
                Operation::PadOnFace {
                    face: wall,
                    profile,
                    dist: 0.4,
                }
            } else {
                Operation::PocketOnFace {
                    face: wall,
                    profile,
                    dist: 0.4,
                }
            };
            let out = match apply(&mut m, &op).expect("the wf feature must build (S4)") {
                OpOutput::PadOnFace { solid, .. } | OpOutput::PocketOnFace { solid, .. } => solid,
                _ => unreachable!(),
            };
            m.rebuild_adjacency();
            record(&format!("wf {kn}"), &m, &inputs, &Ok(vec![out]));
        }
    }

    // ── **Contact — two bodies that meet along a line or at a point, and one pair that cannot part.**
    //
    // ★ The population this file was blind to. Before the contact work every one of these was a
    // reject, and census carried **no reject line at all** (measured: 152 lines, 7 `EMPTY`, zero
    // `ERR`), so the whole family was invisible to the gate that reads every commit. It records
    // *both* halves on purpose: the separable contacts coming back as the bodies they are, and the
    // pinch that stays a reject because the material loops around it.
    for (kn, k) in KINDS {
        // Two cubes sharing exactly the vertical line x=1, y=1.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        m.rebuild_adjacency();
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, k, a, b);
        m.rebuild_adjacency();
        record(&format!("ct edge {kn}"), &m, &inputs, &out);

        // Two cubes sharing exactly the point (1,1,1).
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        m.rebuild_adjacency();
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, k, a, b);
        m.rebuild_adjacency();
        record(&format!("ct corner {kn}"), &m, &inputs, &out);
    }
    // ★★ **The pinch that must stay a reject.** A and B meet only along the line x=2, y=2, but a
    // bridge overlaps both, so the material loops around the contact: cut it there and one piece
    // remains, not two. A line that ever stops saying `ERR` here is the separation going too far.
    {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        let bridge = m.add_cuboid(
            Point3::from_array([1.0, 0.3, 0.2]),
            Point3::from_array([3.0, 3.0, 0.8]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 2.0, 0.0]),
            Point3::from_array([4.0, 4.0, 1.0]),
        );
        m.rebuild_adjacency();
        let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
        m.rebuild_adjacency();
        let inputs = operands(&m, ab[0], b);
        let out = boolean(&mut m, BoolKind::Fuse, ab[0], b);
        m.rebuild_adjacency();
        record("ct ring-pinch fuse", &m, &inputs, &out);
    }
    // ── Cylinders (M6-0): the truth rides beside the cache. Solo and moved bodies are digest
    // lines (a boolean never runs); the boolean rows record the honest refusal — that reject
    // string is part of the corpus, so the day M6-2 admits cylinders, these lines change from
    // ERR to results *in the diff*, not silently.
    {
        use nacre_math::Vector3;
        let solo = |m: &Model, s: Handle<Solid>| {
            let (vn, vh) = coord_digest(m, s);
            let (pn, ph) = plane_digest(m, s);
            format!("v{vn}h{vh:016x}p{pn}h{ph:016x}")
        };
        let mut m = Model::new();
        let c = m.add_cylinder(
            Point3::from_array([0.5, -1.25, 2.0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            1.5,
            2.5,
        );
        m.rebuild_adjacency();
        println!("c cyl solo {}", solo(&m, c));
        // An exact 90° turn (truth transported, nothing recorded) and an inexact 31° turn
        // (node recorded, def carried verbatim) — both roads in the corpus.
        let turned = xf(
            &mut m,
            c,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(90)).expect("angle"),
            }),
        );
        println!("c cyl turn90 {}", solo(&m, turned));
        let leaned = xf(
            &mut m,
            turned,
            Isometry::rotation(Rotation {
                axis: Axis::X,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(31)).expect("angle"),
            }),
        );
        println!("c cyl turn31 {}", solo(&m, leaned));

        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 2.0, 2.0]),
            );
            let b = m.add_cylinder(
                Point3::from_array([1.0, 1.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("cyl {kn}"), &m, &inputs, &out);
        }
        // ── **Cut rims** (M6-2b green): a boss whose circle a boundary segment cuts — the arc
        // population. Three placements: straddling the plate's top edge (seam ≡ branch), turned
        // over the corner (the seam splits the wrap arc), and hung under the bottom edge (the
        // cut circle is the band's hi end).
        for (pn, origin, axis) in [
            ("straddle", [4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ("turned", [4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ("hung", [4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([4.0, 4.0, 2.0]),
                );
                let b = m.add_cylinder(
                    Point3::from_array(origin),
                    Vector3::from_array(axis),
                    0.5,
                    1.0,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arc {pn} {kn}"), &m, &inputs, &out);
            }
        }
        // ── **A bored plate, then a straddling boss** (chaining wall 3): the second boolean's
        // nesting reads a bitten top ring — branch corners and an arc step — so the containment
        // parity runs the mixed road instead of the rational chart. The first boolean is fixed
        // (`cut` the through-bore); the second varies by kind.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let bore = m.add_cylinder(
                Point3::from_array([1.0, 1.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            let bored = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            let boss = m.add_cylinder(
                Point3::from_array([4.0, 2.0, 2.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, bored, boss);
            let out = boolean(&mut m, k, bored, boss);
            m.rebuild_adjacency();
            record(&format!("arc bored straddle {kn}"), &m, &inputs, &out);
        }
        // ── **The rulings road** (cell 4, the gate's record-and-pass arm): a boss whose axis
        // lies exactly on the plate's wall plane builds (through, the corner's two walls, an
        // asymmetric station); the arm's deliberate exclusions record their refusals — offset
        // (`0 < d < r`) and tangent (`d = r`) keep `WallMeetsLateral` at the gate, the
        // half-height and flush-cap variants walk to the ladder's own refusals.
        for (pn, base, h) in [
            ("through", [40.0, 20.0, -10.0], 50.0),
            ("corner", [40.0, 40.0, -10.0], 50.0),
            ("offmid", [40.0, 10.0, -10.0], 50.0),
            ("offset", [38.0, 20.0, -10.0], 50.0),
            ("tangent", [35.0, 20.0, -10.0], 50.0),
            ("half", [40.0, 20.0, -10.0], 20.0),
            ("flush", [40.0, 20.0, 0.0], 20.0),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([40.0, 40.0, 20.0]),
                );
                let b = m.add_cylinder(
                    Point3::from_array(base),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    5.0,
                    h,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("rul {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **Face-to-face contact** (`ct2`): two cells meeting on a full wall, which is what a
    // pattern of parts is made of. The shared wall is interior, so the corners on it must
    // dissolve — and they only can once each cell's caps merge, which is what the coplanar
    // merge's circle carry opened. Bored, pocketed+bored, and a three-cell chain (the second
    // fuse meets a result that already carries merged caps).
    {
        let cell = |m: &mut Model, x0: f64, pocket: bool| {
            let plate = m.add_cuboid(
                Point3::from_array([x0, 0.0, 0.0]),
                Point3::from_array([x0 + 20.0, 20.0, 10.0]),
            );
            m.rebuild_adjacency();
            let mut out = plate;
            if pocket {
                let p = m.add_cuboid(
                    Point3::from_array([x0 + 2.0, 2.0, 4.0]),
                    Point3::from_array([x0 + 8.0, 8.0, 10.0]),
                );
                m.rebuild_adjacency();
                out = boolean(m, BoolKind::Cut, out, p).expect("the pocket cuts")[0];
                m.rebuild_adjacency();
            }
            let bore = m.add_cylinder(
                Point3::from_array([x0 + 14.0, 14.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.0,
                12.0,
            );
            m.rebuild_adjacency();
            out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            out
        };
        for (pn, pocket) in [("bored", false), ("pocketed", true)] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = cell(&mut m, 0.0, pocket);
                let b = cell(&mut m, 20.0, pocket);
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("ct2 {pn} {kn}"), &m, &inputs, &out);
            }
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = cell(&mut m, 0.0, false);
            let b = cell(&mut m, 20.0, false);
            let ab = boolean(&mut m, BoolKind::Fuse, a, b).expect("the first pair fuses")[0];
            m.rebuild_adjacency();
            let c = cell(&mut m, 40.0, false);
            let inputs = operands(&m, ab, c);
            let out = boolean(&mut m, k, ab, c);
            m.rebuild_adjacency();
            record(&format!("ct2 chain {kn}"), &m, &inputs, &out);
        }
    }
    // ── **Translated cylinders** (`trc`): the population the moved-cylinder road opened. A
    // non-dyadic offset records a chain, so each of these carries one (the dyadic twin records
    // nothing and is the same corpus row the `cyl` family already holds). Four placements: a
    // tool cutting, a boss fusing, a bored body fused onto its twin (the plane side), and a
    // bore translated **onto** another bore — whose axes then coincide, which the cylinder-pair
    // rule refuses by name (`CylinderPairContact`); before the road opened it never got that
    // far and said `CylinderGateUndecided` instead.
    {
        let off = |x: f64, y: f64, z: f64| {
            Isometry::translation([x, y, z].map(|c| Rat::try_from_f64(c).expect("rational")))
        };
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let tool = m.add_cylinder(
                Point3::from_array([7.3, 7.3, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let tool = xf(&mut m, tool, off(10.7, 0.0, 0.0));
            let inputs = operands(&m, plate, tool);
            let out = boolean(&mut m, k, plate, tool);
            m.rebuild_adjacency();
            record(&format!("trc tool {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array([10.0, 20.0, 5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                3.0,
                8.0,
            );
            m.rebuild_adjacency();
            let boss = xf(&mut m, boss, off(15.3, 0.1, 0.0));
            let inputs = operands(&m, plate, boss);
            let out = boolean(&mut m, k, plate, boss);
            m.rebuild_adjacency();
            record(&format!("trc boss {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let bored = |m: &mut Model| {
                let plate = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([40.0, 40.0, 10.0]),
                );
                let bore = m.add_cylinder(
                    Point3::from_array([7.3, 7.3, -5.0]),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    2.1,
                    30.0,
                );
                m.rebuild_adjacency();
                let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
                m.rebuild_adjacency();
                out
            };
            let a = bored(&mut m);
            let b = bored(&mut m);
            let b = xf(&mut m, b, off(25.7, 3.3, 2.0));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("trc body {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let bore = m.add_cylinder(
                Point3::from_array([18.0, 20.0, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let holed = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            let twin = m.add_cylinder(
                Point3::from_array([7.3, 20.0, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let twin = xf(&mut m, twin, off(10.7, 0.0, 0.0)); // lands on the bore's own axis
            let inputs = operands(&m, holed, twin);
            let out = boolean(&mut m, k, holed, twin);
            m.rebuild_adjacency();
            record(&format!("trc onaxis {kn}"), &m, &inputs, &out);
        }
    }
    // ── **A grid built in two generations** (`xy`): an array joined along x, then joined again
    // along y. What makes this its own population is not the second axis but the second
    // *generation* — a row is a body whose surfaces come from two provenances, so moving it puts
    // carriers with different chains on one corner, and the corner road answers those in the
    // world (no shared frame exists to answer them in). The `ct2` family above never reaches it:
    // there each cell is *built* at its place, so nothing carries a chain at all.
    //
    // The cells are deliberately feature-poor (one bore) — this table runs twice at every gate,
    // and the user-scale cell costs tens of seconds. What the rows have to pin is the road, and
    // one bore already puts a cylinder gate on a twice-moved wall.
    {
        let cell = |m: &mut Model| {
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([20.0, 20.0, 10.0]),
            );
            let bore = m.add_cylinder(
                Point3::from_array([6.3, 6.3, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                12.0,
            );
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            out
        };
        let row = |m: &mut Model, n: usize, step: [f64; 3]| {
            let mut acc = cell(m);
            for i in 1..n {
                let c = cell(m);
                let c = xf(
                    m,
                    c,
                    Isometry::translation(
                        step.map(|v| Rat::try_from_f64(v * i as f64).expect("rational")),
                    ),
                );
                acc = boolean(m, BoolKind::Fuse, acc, c).expect("the row fuses")[0];
                m.rebuild_adjacency();
            }
            acc
        };
        // x-then-y and y-then-x: the same grid by two routes. The second row is a *moved result*,
        // which is the shape the third door exists for.
        for (an, first, second) in [
            ("grid", [20.0, 0.0, 0.0], [0.0, 20.0, 0.0]),
            ("yx", [0.0, 20.0, 0.0], [20.0, 0.0, 0.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = row(&mut m, 2, first);
                let b = row(&mut m, 2, first);
                let b = xf(
                    &mut m,
                    b,
                    Isometry::translation(second.map(|v| Rat::try_from_f64(v).expect("rational"))),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("xy {an} {kn}"), &m, &inputs, &out);
            }
        }
        // Three in x, then y: the moved row carries three chains, not two.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = row(&mut m, 3, [20.0, 0.0, 0.0]);
            let b = row(&mut m, 3, [20.0, 0.0, 0.0]);
            let b = xf(
                &mut m,
                b,
                Isometry::translation(
                    [0.0, 20.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                ),
            );
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("xy three {kn}"), &m, &inputs, &out);
        }
        // The negative control: a turned row has no world statement, so it must keep taking the
        // old roads — a rotation does not fold into a translation and this door never opens for it.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = row(&mut m, 2, [20.0, 0.0, 0.0]);
            let b = row(&mut m, 2, [20.0, 0.0, 0.0]);
            let b = xf(
                &mut m,
                b,
                Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                }),
            );
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("xy turned {kn}"), &m, &inputs, &out);
        }
    }
    // ── **A datum through a turned solid's corners, and a boolean over it** (`dt`).
    //
    // ★★★★★ **This family exists because its absence hid a regression for 169 commits.** Every
    // datum above is `DatumDef::Stated`; nothing in this corpus stated one **through vertices**,
    // so when the invariant-plane restatement made `vertex_meet` read a turned solid's corners as
    // straddling — costing every one of them the named datum road — the census crossed it
    // bit-identical and said nothing. The rule (`docs/truth-and-cache.md`, gate 8b) is to put the
    // population a change touches into the ledger, and this is that population.
    //
    // ★★ **It does not stop at making the datum.** The plane a `ThroughVertices` datum mints is
    // where `collect_planes` has to choose between the exact witness road and the judged one, and
    // that fork is only reached when such a plane meets a boolean. So the datum hosts a prism and
    // the prism is an operand: a row that only *built* a datum would never walk the road worth
    // watching.
    for (kn, k) in KINDS {
        let mut m = Model::new();
        let block = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let turned = xf(
            &mut m,
            block,
            Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
            }),
        );
        // Three corners of the turned block, in handle order so the statement is replay-stable.
        let mut vs: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
        {
            let sol = m.solids.get(turned);
            'pick: for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
                for &fh in &m.shells.get(sh).faces {
                    for &he in &m.faces.get(fh).outer.half_edges {
                        let vh = m.he_start(he);
                        if !vs.contains(&vh) {
                            vs.push(vh);
                        }
                        if vs.len() == 3 {
                            break 'pick;
                        }
                    }
                }
            }
        }
        vs.sort_by_key(|v| v.index());
        let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
            &mut m,
            &Operation::DatumPlane {
                def: DatumDef::ThroughVertices([vs[0], vs[1], vs[2]]),
            },
        ) else {
            panic!("a datum through three corners of a turned block")
        };
        let OpOutput::Extrude { solid: tool, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: nacre_ops::Profile2d::polygon(vec![
                    nacre_math::Point2::from_array([-0.6, -0.6]),
                    nacre_math::Point2::from_array([0.6, -0.6]),
                    nacre_math::Point2::from_array([0.6, 0.6]),
                    nacre_math::Point2::from_array([-0.6, 0.6]),
                ])
                .unwrap(),
                dist: 1.4,
            },
        )
        .expect("a prism on the datum") else {
            unreachable!("extrude yields Extrude output")
        };
        m.rebuild_adjacency();
        let inputs = operands(&m, turned, tool);
        let out = boolean(&mut m, k, turned, tool);
        m.rebuild_adjacency();
        record(&format!("dt {kn}"), &m, &inputs, &out);
    }

    // ★★★★★ **The link that turns "interning explains it" into something falsifiable.** Since S2
    // every plane with points is named (wide ones in the arbitrary-precision vessel), so a
    // coordinate can move only when wide planes *merge* — a `c ` line that moves must come with
    // a nonzero count here. If the coordinates move and this stays zero, the cause is something
    // else and the diff is not explained.
    println!(
        "stat wide_planes {}",
        nacre_topo::WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed)
    );
    // S9: the falsifiability bridge for seeding-shaped changes. A `c ` diff that the wide
    // counter cannot explain (name-collision populations are narrow) must come with a nonzero
    // count here instead — pushes that interned onto a seeded world plane.
    println!(
        "stat seeded_hits {}",
        nacre_topo::SEEDED_HITS.load(std::sync::atomic::Ordering::Relaxed)
    );
}

/// A square prism on a plane through the origin with normal `n` — the tilted twin of [`ex`].
///
/// ★ The footprint is a **U** (a notched square), so two of its walls lie on one plane. Where the
/// kernel can name that plane the two share a `Surface` handle; where it cannot they are two
/// handles for one plane, which is the state this line of work removes.
fn tilted_prism(m: &mut Model, n: [f64; 3], off: f64, size: f64, dist: f64) -> Handle<Solid> {
    use nacre_math::{Point2, Vector3};
    let p = |x: f64, y: f64| Point2::from_array([off + x * size, off + y * size]);
    let plane =
        nacre_ops::SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array(n))
            .expect("a tilted plane");
    let __g208 = datum_frame(m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g208,
            // U-shaped: the notch's two side walls are the coplanar pair.
            profile: nacre_ops::Profile2d::polygon(vec![
                p(0.0, 0.0),
                p(3.0, 0.0),
                p(3.0, 2.0),
                p(2.0, 2.0),
                p(2.0, 1.0),
                p(1.0, 1.0),
                p(1.0, 2.0),
                p(0.0, 2.0),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("tilted extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    solid
}

/// A rectangular prism raised from `z` by `dist` — the construction path production uses.
fn ex(m: &mut Model, z: f64, lo: [f64; 2], hi: [f64; 2], dist: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let __g207 = datum_frame(
        m,
        nacre_ops::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
    );
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g207,
            profile: nacre_ops::Profile2d::polygon(vec![
                p(lo[0], lo[1]),
                p(hi[0], lo[1]),
                p(hi[0], hi[1]),
                p(lo[0], hi[1]),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    solid
}
