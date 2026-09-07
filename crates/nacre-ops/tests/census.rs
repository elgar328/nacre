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
//! exist for that reason, and so does `mot` (the boss corpus under rigid motion — cell ④):
//! until it was added, no row held a **rotated** cylinder boolean at all. **Read "bit-identical"
//! as evidence only about the population present.**

use nacre_math::Point3;
use nacre_ops::{BoolKind, Edge2d, OpOutput, Operation, apply, boolean, from_edges};
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
        // asymmetric station); the arm's deliberate exclusions are gone — offset (`0 < d < r`,
        // cell ③) and tangent (`d = r`, cell ⑥) both pass the gate now, and the tangent row's
        // three kinds are where the *verdict* answers (`Fuse` builds, `Cut` pinches).
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
    // ── **The same boss, on each of the four walls** (`wal`). A part does not care which side of
    // itself a boss sits on, and neither should the kernel — but the rulings road was built and
    // measured on the `+x` wall alone, and the corpus inherited that: every `rul` row above puts
    // the axis on a max-side wall. These rows put the same solid on all four, so the table can
    // see a rule that holds on one side and not the other.
    //
    // ★ **Recorded before the fix, deliberately.** As this family lands, **eleven** of these
    // twenty-two rows are refusals: six `RingOrientation` (the min-side walls, all three kinds),
    // two `MissingSeam` (`+y` cut and common, where `+y` fuse builds), and three
    // `OpenResultShell` (the min corner). Writing them down first is what makes the next commit's
    // diff the evidence: without it the rows would be born green and the corpus could not say
    // what changed.
    //
    // The box twins are the control that says what breaks needs the **cylinder's** edges: the
    // same straddle with a cuboid tool builds on every wall, volume 35 exactly. ★ Which *kind* of
    // cylinder edge it is, is not something this table can say — measured later, the wrong turn
    // is read on a **ruling**, and `turn` itself is correct. A population fact is not a code fact.
    // Fuse alone for the twins: the control only has to say that the planar straddle builds.
    {
        let plate = |m: &mut Model| {
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            m.rebuild_adjacency();
            a
        };
        for (pn, base) in [
            ("xlo", [0.0, 2.0, -1.0]),
            ("xhi", [4.0, 2.0, -1.0]),
            ("ylo", [2.0, 0.0, -1.0]),
            ("yhi", [2.0, 4.0, -1.0]),
            ("corner-lo", [0.0, 0.0, -1.0]),
            ("corner-hi", [4.0, 4.0, -1.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = plate(&mut m);
                let b = m.add_cylinder(
                    Point3::from_array(base),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    0.5,
                    4.0,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("wal {pn} {kn}"), &m, &inputs, &out);
            }
        }
        for (pn, lo, hi) in [
            ("box-xlo", [-0.5, 1.5, -1.0], [0.5, 2.5, 3.0]),
            ("box-xhi", [3.5, 1.5, -1.0], [4.5, 2.5, 3.0]),
            ("box-ylo", [1.5, -0.5, -1.0], [2.5, 0.5, 3.0]),
            ("box-yhi", [1.5, 3.5, -1.0], [2.5, 4.5, 3.0]),
        ] {
            let mut m = Model::new();
            let a = plate(&mut m);
            let b = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, BoolKind::Fuse, a, b);
            m.rebuild_adjacency();
            record(&format!("wal {pn} fuse"), &m, &inputs, &out);
        }
    }
    // ── **The boss corpus under rigid motion** (`mot`, cell ④): the production-side rows of the
    // commuting oracle (`tests.rs::the_boolean_commutes_with_rigid_motion`), one per sign class
    // the oracle names — ∥ wall classes with `frame_sign = −1` (a max-side wall put on a seed
    // plane by a translation; a wall turned onto one by rz90), ⊥ classes with `(axis_up, frame)`
    // `= (true, −1)` (corner-lo with the axis turned to −y) and `(false, −1)` (the top cap put on
    // z = 0), the transport's exactness boundary (the offset boss under a rigid motion) and the
    // recorded path (a non-dyadic translation). ★ Rows were added only where **both** profiles
    // dump them: the transport row (`offset-out` under `rz90 + t(5,−3,2)`) came with the
    // transport law (cell ④ stage 3) — before it, `world_cylinder_def`'s postcondition (a
    // `debug_assert`) took the dev census down there while release answered «disjoint» (Fuse 2
    // bodies, Cut the plate untouched, Common empty): the silent wrong answer the law closed.
    {
        let boss = |m: &mut Model, base: [f64; 3]| -> (Handle<Solid>, Handle<Solid>) {
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let b = m.add_cylinder(
                Point3::from_array(base),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            (plate, b)
        };
        let rot = |ax: Axis, deg: i128| {
            Isometry::rotation(Rotation {
                axis: ax,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
            })
        };
        let t = |x: i128, y: i128, z: i128| {
            Isometry::translation([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)])
        };
        let rows: Vec<(&str, [f64; 3], Isometry)> = vec![
            ("par wall+x t", [4.0, 2.0, -1.0], t(-4, -4, -2)),
            ("par wall-y rz90", [2.0, 0.0, -1.0], rot(Axis::Z, 90)),
            ("perp corner-lo rx90", [0.0, 0.0, -1.0], rot(Axis::X, 90)),
            ("perp wall+y t", [2.0, 4.0, -1.0], t(-4, -4, -2)),
            (
                "xport offset-out rz90+t",
                [4.3, 2.0, -1.0],
                Isometry::rigid(
                    Rotation {
                        axis: Axis::Z,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(90)).expect("angle"),
                    },
                    [Rat::from_int(5), Rat::from_int(-3), Rat::from_int(2)],
                ),
            ),
            (
                "rec through t",
                [2.0, 2.0, -1.0],
                Isometry::translation([
                    Rat::new(7, 11).expect("rat"),
                    Rat::new(3, 10).expect("rat"),
                    Rat::new(1, 4).expect("rat"),
                ]),
            ),
        ];
        for (name, base, iso) in &rows {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (p, b) = boss(&mut m, *base);
                let p = xf(&mut m, p, *iso);
                let b = xf(&mut m, b, *iso);
                let inputs = operands(&m, p, b);
                let out = boolean(&mut m, k, p, b);
                m.rebuild_adjacency();
                record(&format!("mot {name} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A cap that lies in another face's plane** (`cap`): the population where a *disk* is a
    // face of the result and the face around it holds the same circle as a hole. The two are
    // adjacent across that circle and nothing else — a disk has no nodes — so before the merge
    // learned to read circles, this table could not have held the answer: every row here either
    // was refused outright or came back with a plane split in two.
    //
    // ★ **Why the census had nothing to say about that change**: it contained no such fixture at
    // all (measured — the diff over the whole table was empty), which is the same lesson three
    // cells running: **the corpus does not contain the population a cell opens until the cell
    // adds it.** These rows are that addition. The bore row is the negative control: its circle
    // is shared with a *cylinder* face, which is not in the plane's group, so it must stay a hole.
    for (pn, base, h) in [
        ("ontop", [2.0, 2.0, 2.0], 1.0), // a boss standing on the top face (contact)
        ("sunk", [2.0, 2.0, 1.0], 1.0),  // a boss whose cap is flush with the top, inside
        ("through", [2.0, 2.0, 0.0], 3.0), // a boss whose cap is flush with the bottom
        ("bore", [2.0, 2.0, -1.0], 4.0), // the negative control: a real hole
    ] {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let b = m.add_cylinder(
                Point3::from_array(base),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                h,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("cap {pn} {kn}"), &m, &inputs, &out);
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
        // ★★ **The row that actually needs the world road.** Everything above builds whether or
        // not the third door exists — measured, by switching the door off and watching these rows
        // come back byte-identical. A grid only reaches a corner the population gate must judge
        // when the cell carries enough features, and the relation is **not monotone** (two pockets
        // with four bores build; one of those pockets with two *other* bores does not). So this row holds
        // the smallest configuration measured to refuse `CylinderGateUndecided` with the door off,
        // taken from the user's own part: one pocket and two through-bores in an 86 mm cell.
        // Treat the numbers as pinned — nudging one moves the row out of its population.
        {
            let featured = |m: &mut Model| {
                let plate = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([86.0, 86.0, 71.5]),
                );
                let pocket = m.add_cuboid(
                    Point3::from_array([1.0, 38.6, 0.0]),
                    Point3::from_array([33.2, 58.6, 68.5]),
                );
                m.rebuild_adjacency();
                let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
                m.rebuild_adjacency();
                let mut out = out;
                for (c, r) in [([17.1, 8.9], 2.22), ([46.75, 37.35], 2.34)] {
                    let bore = m.add_cylinder(
                        Point3::from_array([c[0], c[1], 0.0]),
                        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                        r,
                        143.0,
                    );
                    m.rebuild_adjacency();
                    out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
                    m.rebuild_adjacency();
                }
                out
            };
            let featured_row = |m: &mut Model| {
                let a = featured(m);
                let b = featured(m);
                let b = xf(
                    m,
                    b,
                    Isometry::translation(
                        [86.0, 0.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                    ),
                );
                let out = boolean(m, BoolKind::Fuse, a, b).expect("the row fuses")[0];
                m.rebuild_adjacency();
                out
            };
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = featured_row(&mut m);
                let b = featured_row(&mut m);
                let b = xf(
                    &mut m,
                    b,
                    Isometry::translation(
                        [0.0, 86.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                    ),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("xy needsworld {kn}"), &m, &inputs, &out);
            }
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
    // ── **Arc profiles** (cell ⑨): sketched circles and arcs extruded, then met by a box. The
    // half disk's chord wall *crosses* its cylinder (branch corners `Lo`/`Hi`), the annulus and
    // the bored plate have only whole circles; the slot's straight walls are *tangent* to its half
    // cylinders, the ruling the tracer has no side for — its rows record that refusal by name.
    {
        let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let sketch = |m: &mut Model, edges: Vec<Edge2d>, dist: f64| -> Handle<Solid> {
            let profile = from_edges(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, Axis::Z);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the arc profile extrudes") else {
                unreachable!()
            };
            solid
        };
        type Shape = Box<dyn Fn(&mut Model) -> Handle<Solid>>;
        let shapes: [(&str, Shape); 4] = [
            (
                "halfdisk",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            Edge2d::line(p2(2.0, 5.0), p2(2.0, -1.0)).unwrap(),
                            Edge2d::arc_turns(p2(2.0, 2.0), p2(2.0, -1.0), 2).unwrap(),
                        ],
                        3.0,
                    )
                }),
            ),
            (
                "annulus",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            Edge2d::circle(p2(2.0, 2.0), 3.0).unwrap(),
                            Edge2d::circle(p2(2.0, 2.0), 1.5).unwrap(),
                        ],
                        3.0,
                    )
                }),
            ),
            (
                "bored",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            Edge2d::line(p2(-1.0, -1.0), p2(5.0, -1.0)).unwrap(),
                            Edge2d::line(p2(5.0, -1.0), p2(5.0, 5.0)).unwrap(),
                            Edge2d::line(p2(5.0, 5.0), p2(-1.0, 5.0)).unwrap(),
                            Edge2d::line(p2(-1.0, 5.0), p2(-1.0, -1.0)).unwrap(),
                            Edge2d::circle(p2(2.0, 2.0), 1.0).unwrap(),
                        ],
                        3.0,
                    )
                }),
            ),
            (
                "slot",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            Edge2d::line(p2(0.0, 1.0), p2(4.0, 1.0)).unwrap(),
                            Edge2d::arc_turns(p2(4.0, 2.0), p2(4.0, 1.0), 2).unwrap(),
                            Edge2d::line(p2(4.0, 3.0), p2(0.0, 3.0)).unwrap(),
                            Edge2d::arc_turns(p2(0.0, 2.0), p2(0.0, 3.0), 2).unwrap(),
                        ],
                        3.0,
                    )
                }),
            ),
        ];
        for (sn, build) in &shapes {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = build(&mut m);
                let b = m.add_cuboid(
                    Point3::from_array([1.0, 0.0, 1.0]),
                    Point3::from_array([6.0, 4.0, 5.0]),
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arcprofile {sn} {kn}"), &m, &inputs, &out);
            }
        }
    }

    // ── **Arc walls** (cell ⑩): the user's first assembly and its controls. A filleted plate
    // with two holes (XY), a standing plate with a slot window (ZX, moved to y ∈ [3, 4]) and two
    // triangular gussets (YZ) met three walls at once — same-solid parallel cylinder pairs, the
    // oblique gusset plane, and the fillets' tangent rulings. Every row is a pair the assembly's
    // fold visits or a control that isolates one wall; the plan's prediction table names what
    // each stage flips.
    {
        // Helpers as items, not closures: the shape builders below are boxed and `move`d, and a
        // closure they borrowed would not live long enough.
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Edge2d>, dist: f64) -> Handle<Solid> {
            let profile = from_edges(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the arc-wall profile extrudes") else {
                unreachable!()
            };
            solid
        }
        // An exact translation by decimals — the anchor arithmetic a script hands the kit.
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn far_box(m: &mut Model) -> Handle<Solid> {
            m.add_cuboid(
                Point3::from_array([20.0, 20.0, 20.0]),
                Point3::from_array([21.0, 21.0, 21.0]),
            )
        }
        fn circle_prism(m: &mut Model, c: [f64; 2], r: f64, dist: f64) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![Edge2d::circle(p2(c[0], c[1]), r).unwrap()],
                dist,
            )
        }
        // p1: the plate, bottom corners filleted r 2 with the holes centred on the fillet axes.
        let p1 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-1.5, -4.0), p2(1.5, -4.0)).unwrap(),
                    Edge2d::arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1).unwrap(),
                    Edge2d::line(p2(3.5, -2.0), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -2.0)).unwrap(),
                    Edge2d::arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1).unwrap(),
                    Edge2d::circle(p2(-1.5, -2.0), 1.0).unwrap(),
                    Edge2d::circle(p2(1.5, -2.0), 1.0).unwrap(),
                ],
                1.0,
            )
        };
        // p2: the standing plate on ZX (sketch u = z, v = x) with a 2×1 slot window, moved to y = 3.
        let p2p = |m: &mut Model| {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    Edge2d::line(p2(1.0, -3.5), p2(6.0, -3.5)).unwrap(),
                    Edge2d::line(p2(6.0, -3.5), p2(6.0, 3.5)).unwrap(),
                    Edge2d::line(p2(6.0, 3.5), p2(1.0, 3.5)).unwrap(),
                    Edge2d::line(p2(1.0, 3.5), p2(1.0, -3.5)).unwrap(),
                    Edge2d::line(p2(3.0, -0.5), p2(4.0, -0.5)).unwrap(),
                    Edge2d::arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2).unwrap(),
                    Edge2d::line(p2(4.0, 0.5), p2(3.0, 0.5)).unwrap(),
                    Edge2d::arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        };
        // p3: a triangular gusset on YZ (sketch u = y, v = z), moved to x = 1.5 as the script does.
        let p3 = |m: &mut Model| {
            let s = prism(
                m,
                Axis::X,
                vec![
                    Edge2d::line(p2(0.0, 1.0), p2(3.0, 1.0)).unwrap(),
                    Edge2d::line(p2(3.0, 1.0), p2(3.0, 6.0)).unwrap(),
                    Edge2d::line(p2(3.0, 6.0), p2(0.0, 1.0)).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [1.5, 0.0, 0.0])
        };
        // Controls.
        let fillet15 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-2.0, -4.0), p2(2.0, -4.0)).unwrap(),
                    Edge2d::arc_turns(p2(2.0, -2.5), p2(2.0, -4.0), 1).unwrap(),
                    Edge2d::line(p2(3.5, -2.5), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -2.5)).unwrap(),
                    Edge2d::arc_turns(p2(-2.0, -2.5), p2(-3.5, -2.5), 1).unwrap(),
                ],
                1.0,
            )
        };
        let slot30 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(0.0, -5.0), p2(30.0, -5.0)).unwrap(),
                    Edge2d::arc_turns(p2(30.0, 0.0), p2(30.0, -5.0), 2).unwrap(),
                    Edge2d::line(p2(30.0, 5.0), p2(0.0, 5.0)).unwrap(),
                    Edge2d::arc_turns(p2(0.0, 0.0), p2(0.0, 5.0), 2).unwrap(),
                ],
                2.0,
            )
        };
        let holes_plate = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-3.5, -4.0), p2(3.5, -4.0)).unwrap(),
                    Edge2d::line(p2(3.5, -4.0), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -4.0)).unwrap(),
                    Edge2d::circle(p2(-1.5, -2.0), 1.0).unwrap(),
                    Edge2d::circle(p2(1.5, -2.0), 1.0).unwrap(),
                ],
                1.0,
            )
        };
        let rrect = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(2.0, 0.0), p2(8.0, 0.0)).unwrap(),
                    Edge2d::arc_turns(p2(8.0, 2.0), p2(8.0, 0.0), 1).unwrap(),
                    Edge2d::line(p2(10.0, 2.0), p2(10.0, 6.0)).unwrap(),
                    Edge2d::arc_turns(p2(8.0, 6.0), p2(10.0, 6.0), 1).unwrap(),
                    Edge2d::line(p2(8.0, 8.0), p2(2.0, 8.0)).unwrap(),
                    Edge2d::arc_turns(p2(2.0, 6.0), p2(2.0, 8.0), 1).unwrap(),
                    Edge2d::line(p2(0.0, 6.0), p2(0.0, 2.0)).unwrap(),
                    Edge2d::arc_turns(p2(2.0, 2.0), p2(0.0, 2.0), 1).unwrap(),
                ],
                3.0,
            )
        };
        let tube = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::circle(p2(2.0, 2.0), 3.0).unwrap(),
                    Edge2d::circle(p2(2.0, 2.0), 1.5).unwrap(),
                ],
                3.0,
            )
        };
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            ("p1p2", Box::new(move |m| (p1(m), p2p(m)))),
            ("p1p3", Box::new(move |m| (p1(m), p3(m)))),
            ("p2p3", Box::new(move |m| (p2p(m), p3(m)))),
            ("fillet15-far", Box::new(move |m| (fillet15(m), far_box(m)))),
            ("slot-far", Box::new(move |m| (slot30(m), far_box(m)))),
            ("holes-gusset", Box::new(move |m| (holes_plate(m), p3(m)))),
            (
                "rrect-box",
                Box::new(move |m| {
                    let a = rrect(m);
                    let b = m.add_cuboid(
                        Point3::from_array([3.0, 2.0, 1.0]),
                        Point3::from_array([7.0, 6.0, 5.0]),
                    );
                    (a, b)
                }),
            ),
            (
                "bushing",
                Box::new(move |m| {
                    let a = tube(m);
                    let pin = circle_prism(m, [2.0, 2.0], 1.0, 5.0);
                    (a, shift(m, pin, [0.0, 0.0, -1.0]))
                }),
            ),
            (
                "stacked",
                Box::new(move |m| {
                    let boss = circle_prism(m, [0.0, 0.0], 2.0, 2.0);
                    let pin = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    (boss, shift(m, pin, [0.0, 0.0, 2.0]))
                }),
            ),
            (
                "stacked-same",
                Box::new(move |m| {
                    let a = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    let b = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    (a, shift(m, b, [0.0, 0.0, 2.0]))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arcwalls {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **Four-plane operand vertices** (cell ⑪): a gusset whose apex lands exactly on a wall's top
    // edge makes a result vertex where **four faces** meet. Fed back as an operand, that vertex
    // used to be named once per face — four names, one of them a triple whose planes share a
    // line — and the next boolean refused by whichever symptom its build order met first
    // (`DegenerateWitness` in the parallel build, `FourPlane` sequentially; measured at S0). Since
    // S1 the vertex names itself from its incident classes and every pair builds. The rows hold
    // both operand orders, the user's four-part fold in both fold orders, the
    // near misses (apex above and below the edge), a box that shares only the top plane's class,
    // and the mirrored fold (class numbers permuted).
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Edge2d>, dist: f64) -> Handle<Solid> {
            let profile = from_edges(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the four-plane profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn mirror_x(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
            let OpOutput::Mirror { solid } = apply(
                m,
                &Operation::Mirror {
                    solid: s,
                    axis: Axis::X,
                    offset: Rat::from_int(0),
                },
            )
            .expect("the mirror applies") else {
                unreachable!()
            };
            solid
        }
        // A fuse inside a fixture: the operand a later row feeds back.
        fn fuse(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Fuse, a, b).expect("the fixture's own fuse builds");
            assert_eq!(out.len(), 1, "the fixture's own fuse is one body");
            out[0]
        }
        // The standing wall: 7 × 5 × 1 at y ∈ [3, 4], z ∈ [1, 6] (the slot plate without its slot).
        fn wall(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    Edge2d::line(p2(1.0, -3.5), p2(6.0, -3.5)).unwrap(),
                    Edge2d::line(p2(6.0, -3.5), p2(6.0, 3.5)).unwrap(),
                    Edge2d::line(p2(6.0, 3.5), p2(1.0, 3.5)).unwrap(),
                    Edge2d::line(p2(1.0, 3.5), p2(1.0, -3.5)).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        // The user's slot plate: the wall with its 2 × 1 slot window.
        fn slot_plate(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    Edge2d::line(p2(1.0, -3.5), p2(6.0, -3.5)).unwrap(),
                    Edge2d::line(p2(6.0, -3.5), p2(6.0, 3.5)).unwrap(),
                    Edge2d::line(p2(6.0, 3.5), p2(1.0, 3.5)).unwrap(),
                    Edge2d::line(p2(1.0, 3.5), p2(1.0, -3.5)).unwrap(),
                    Edge2d::line(p2(3.0, -0.5), p2(4.0, -0.5)).unwrap(),
                    Edge2d::arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2).unwrap(),
                    Edge2d::line(p2(4.0, 0.5), p2(3.0, 0.5)).unwrap(),
                    Edge2d::arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        // The user's filleted, bored plate.
        fn plate(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-1.5, -4.0), p2(1.5, -4.0)).unwrap(),
                    Edge2d::arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1).unwrap(),
                    Edge2d::line(p2(3.5, -2.0), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -2.0)).unwrap(),
                    Edge2d::arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1).unwrap(),
                    Edge2d::circle(p2(-1.5, -2.0), 1.0).unwrap(),
                    Edge2d::circle(p2(1.5, -2.0), 1.0).unwrap(),
                ],
                1.0,
            )
        }
        // The gusset on YZ (u = y, v = z), apex at `(3, apex)`, moved to `x`.
        fn gusset(m: &mut Model, x: f64, apex: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::X,
                vec![
                    Edge2d::line(p2(0.0, 1.0), p2(3.0, 1.0)).unwrap(),
                    Edge2d::line(p2(3.0, 1.0), p2(3.0, apex)).unwrap(),
                    Edge2d::line(p2(3.0, apex), p2(0.0, 1.0)).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [x, 0.0, 0.0])
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            (
                "wall+g",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 6.0))
                }),
            ),
            (
                "wall+g swapped",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let b = fuse(m, w, g1);
                    (gusset(m, -2.4, 6.0), b)
                }),
            ),
            (
                "user4 1234",
                Box::new(|m| {
                    let p1 = plate(m);
                    let p2s = slot_plate(m);
                    let ab = fuse(m, p1, p2s);
                    let g1 = gusset(m, 1.4, 6.0);
                    let abc = fuse(m, ab, g1);
                    (abc, gusset(m, -2.4, 6.0))
                }),
            ),
            (
                "user4 1342",
                Box::new(|m| {
                    let p1 = plate(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let ac = fuse(m, p1, g1);
                    let g2 = gusset(m, -2.4, 6.0);
                    let acd = fuse(m, ac, g2);
                    (acd, slot_plate(m))
                }),
            ),
            (
                "apex-5.9",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 5.9);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 5.9))
                }),
            ),
            (
                "apex-6.1",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.1);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 6.1))
                }),
            ),
            (
                "far-box",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    let b = m.add_cuboid(
                        Point3::from_array([-3.0, 0.0, 6.0]),
                        Point3::from_array([-1.0, 2.0, 7.0]),
                    );
                    (a, b)
                }),
            ),
            (
                "mirrored",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    let a = mirror_x(m, a);
                    let g2 = gusset(m, -2.4, 6.0);
                    (a, mirror_x(m, g2))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("fourplane {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A plane through a fillet's axis** (cell ⑫): a class plane that holds a fillet's axis
    // has two rulings on the fillet, and one of them can be the fillet's own tangent ruling with
    // the plate's wall. That line then carries names from two vocabularies — the operand's tangent
    // corner (`Branch … Double`) and the class's ruling crossing (`Branch … Lo/Hi`) — with nothing
    // that knows they are one point, and the lateral's ruling sweep declines by name. The rows:
    // the minimal plate-and-slab in both contact shapes, an off-axis control, and the user's
    // four-part fold at the script's own gusset positions (`1.5`/`−2.5`).
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Edge2d>, dist: f64) -> Handle<Solid> {
            let profile = from_edges(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the tangent-line profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn fuse(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Fuse, a, b).expect("the fixture's own fuse builds");
            assert_eq!(out.len(), 1, "the fixture's own fuse is one body");
            out[0]
        }
        // A 7 × 8 plate with one fillet (r 2) at its bottom-right corner: axis at (1.5, −2), tangent
        // to the bottom wall `y = −4` along `x = 1.5` and to the right wall `x = 3.5` along `y = −2`.
        fn plate1(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-3.5, -4.0), p2(1.5, -4.0)).unwrap(),
                    Edge2d::arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1).unwrap(),
                    Edge2d::line(p2(3.5, -2.0), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -4.0)).unwrap(),
                ],
                1.0,
            )
        }
        fn slab(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
            m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
        }
        // The user's parts: the filleted, bored plate, the slot plate, the gusset at `x`.
        fn plate(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    Edge2d::line(p2(-1.5, -4.0), p2(1.5, -4.0)).unwrap(),
                    Edge2d::arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1).unwrap(),
                    Edge2d::line(p2(3.5, -2.0), p2(3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(3.5, 4.0), p2(-3.5, 4.0)).unwrap(),
                    Edge2d::line(p2(-3.5, 4.0), p2(-3.5, -2.0)).unwrap(),
                    Edge2d::arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1).unwrap(),
                    Edge2d::circle(p2(-1.5, -2.0), 1.0).unwrap(),
                    Edge2d::circle(p2(1.5, -2.0), 1.0).unwrap(),
                ],
                1.0,
            )
        }
        fn slot_plate(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    Edge2d::line(p2(1.0, -3.5), p2(6.0, -3.5)).unwrap(),
                    Edge2d::line(p2(6.0, -3.5), p2(6.0, 3.5)).unwrap(),
                    Edge2d::line(p2(6.0, 3.5), p2(1.0, 3.5)).unwrap(),
                    Edge2d::line(p2(1.0, 3.5), p2(1.0, -3.5)).unwrap(),
                    Edge2d::line(p2(3.0, -0.5), p2(4.0, -0.5)).unwrap(),
                    Edge2d::arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2).unwrap(),
                    Edge2d::line(p2(4.0, 0.5), p2(3.0, 0.5)).unwrap(),
                    Edge2d::arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        fn gusset(m: &mut Model, x: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::X,
                vec![
                    Edge2d::line(p2(0.0, 1.0), p2(3.0, 1.0)).unwrap(),
                    Edge2d::line(p2(3.0, 1.0), p2(3.0, 6.0)).unwrap(),
                    Edge2d::line(p2(3.0, 6.0), p2(0.0, 1.0)).unwrap(),
                ],
                1.0,
            );
            shift(m, s, [x, 0.0, 0.0])
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            // The slab stands on the plate; its face `x = 1.5` holds the fillet's axis. It stops
            // at `y = −3`, short of the plate's bottom wall, so that plane is the only coincidence.
            (
                "fillet-slab",
                Box::new(|m| (plate1(m), slab(m, [1.5, -3.0, 1.0], [2.5, 0.0, 2.0]))),
            ),
            // Same plane, the slab far from the plate: the class still cuts the fillet.
            (
                "fillet-slab far",
                Box::new(|m| (plate1(m), slab(m, [1.5, 10.0, 1.0], [2.5, 12.0, 2.0]))),
            ),
            // Off the axis by 0.1: an ordinary offset wall — the control.
            (
                "slab 1.6",
                Box::new(|m| (plate1(m), slab(m, [1.6, -3.0, 1.0], [2.6, 0.0, 2.0]))),
            ),
            // ★ Measured on the way (S0): a slab reaching the plate's bottom wall plane `y = −4`
            // — its face coplanar with the wall, not overlapping it, its bottom edge on the
            // wall's line past the tangent point — meets *other* walls: on the axis plane
            // `CoincidentNodes`, off it a **cut** that should leave the plate untouched answers
            // `OpenResultShell` (an assembly defect, unfired before). Both frozen here by name
            // for the cell that takes them; not this cell's proposition.
            (
                "slab wall 1.5",
                Box::new(|m| (plate1(m), slab(m, [1.5, -4.0, 1.0], [2.5, 0.0, 2.0]))),
            ),
            (
                "slab wall 1.6",
                Box::new(|m| (plate1(m), slab(m, [1.6, -4.0, 1.0], [2.6, 0.0, 2.0]))),
            ),
            // The user's fold at the script's own gusset position: the third fuse of the script
            // (plate + slot plate, then the gusset at `1.5`). The fourth (`−2.5`) is the ops lock's.
            (
                "user 1.5",
                Box::new(|m| {
                    let p1 = plate(m);
                    let p2s = slot_plate(m);
                    let ab = fuse(m, p1, p2s);
                    (ab, gusset(m, 1.5))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("tangentline {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A fully rounded outline** (cell ⑬): fillet *every* corner of a plate and its cap ring
    // has no three-plane corner left — eight tangencies and nothing else. The predicate that asks
    // "is this ring inside that bore's disk" looks for a witness among three-plane names only, so
    // it finds none and refuses, and the plate cannot enter **any** boolean. The rows say exactly
    // that: three partners that share nothing but the plate (one standing on it, one overlapping,
    // one a hundred units away), and two controls that differ in one thing each — one corner left
    // sharp, and no bores. The last row carries the user's own script, whose *next* wall this cell
    // hands on to the one after it.
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Edge2d>, dist: f64) -> Handle<Solid> {
            let profile = from_edges(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the round-plate profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        // A 90 × 50 × 12 plate whose corners are filleted (r 5) and which carries `bores` holes
        // (d 7) at `(±38, ±18)`. `fillets` says how many of the four corners are rounded, walking
        // counter-clockwise from the bottom-right: with four the outline has **no** sharp corner.
        fn plate(m: &mut Model, fillets: usize, bores: usize) -> Handle<Solid> {
            // The corners in walk order, each with the direction in and the direction out.
            let corner = [
                ([45.0, -25.0], [1.0, 0.0], [0.0, 1.0]),
                ([45.0, 25.0], [0.0, 1.0], [-1.0, 0.0]),
                ([-45.0, 25.0], [-1.0, 0.0], [0.0, -1.0]),
                ([-45.0, -25.0], [0.0, -1.0], [1.0, 0.0]),
            ];
            let r = 5.0;
            let mut edges: Vec<Edge2d> = Vec::new();
            // Where the previous corner left the pen.
            let mut at = {
                let (c, _, d_out) = corner[3];
                if fillets > 3 {
                    [c[0] + r * d_out[0], c[1] + r * d_out[1]]
                } else {
                    c
                }
            };
            for (i, (c, d_in, d_out)) in corner.into_iter().enumerate() {
                if i < fillets {
                    let tin = [c[0] - r * d_in[0], c[1] - r * d_in[1]];
                    let centre = [tin[0] + r * d_out[0], tin[1] + r * d_out[1]];
                    edges.push(Edge2d::line(p2(at[0], at[1]), p2(tin[0], tin[1])).unwrap());
                    edges.push(
                        Edge2d::arc_turns(p2(centre[0], centre[1]), p2(tin[0], tin[1]), 1).unwrap(),
                    );
                    at = [c[0] + r * d_out[0], c[1] + r * d_out[1]];
                } else {
                    edges.push(Edge2d::line(p2(at[0], at[1]), p2(c[0], c[1])).unwrap());
                    at = c;
                }
            }
            for &[x, y] in [[38.0, 18.0], [38.0, -18.0], [-38.0, 18.0], [-38.0, -18.0]]
                .iter()
                .take(bores)
            {
                edges.push(Edge2d::circle(p2(x, y), 3.5).unwrap());
            }
            prism(m, Axis::Z, edges, 12.0)
        }
        // The user's rib: a stepped profile standing on the plate's top face, extruded ±20 in y.
        fn rib(m: &mut Model, x: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    Edge2d::line(p2(12.0, 0.0), p2(12.0, 7.5)).unwrap(),
                    Edge2d::line(p2(12.0, 7.5), p2(27.0, 7.5)).unwrap(),
                    Edge2d::line(p2(27.0, 7.5), p2(27.0, 5.5)).unwrap(),
                    Edge2d::line(p2(27.0, 5.5), p2(62.0, 5.5)).unwrap(),
                    Edge2d::line(p2(62.0, 5.5), p2(62.0, -5.5)).unwrap(),
                    Edge2d::line(p2(62.0, -5.5), p2(27.0, -5.5)).unwrap(),
                    Edge2d::line(p2(27.0, -5.5), p2(27.0, -7.5)).unwrap(),
                    Edge2d::line(p2(27.0, -7.5), p2(12.0, -7.5)).unwrap(),
                    Edge2d::line(p2(12.0, -7.5), p2(12.0, 0.0)).unwrap(),
                ],
                40.0,
            );
            shift(m, s, [x, -20.0, 0.0])
        }
        fn box_at(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
            m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi))
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            // Standing on the plate's top face.
            (
                "box on top",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // Overlapping it.
            (
                "box through",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [10.0, -20.0, 6.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // ★ A hundred units away — the partner has nothing to do with it. This row is what
            // says the wall is the plate's own.
            (
                "box far",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [200.0, -5.0, -5.0], [210.0, 5.0, 5.0]),
                    )
                }),
            ),
            // One corner left sharp — the control that differs in exactly one thing.
            (
                "one sharp corner",
                Box::new(|m| {
                    (
                        plate(m, 3, 4),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // Rounded but unbored — the other control.
            (
                "no bores",
                Box::new(|m| {
                    (
                        plate(m, 4, 0),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // The rib the user's script builds, on one side.
            ("rib", Box::new(|m| (plate(m, 4, 4), rib(m, 17.5)))),
            // ★ The user's own script — plate, both ribs, then the cylinder that bores across
            // them. It could not be a row until the fuse built; now it is one, and the name it
            // lands on is the **next** wall, a class carrying a circle and rulings at once. That
            // is the hand-off to the cell after this one, held in the corpus so it cannot be lost.
            (
                "user script",
                Box::new(|m| {
                    let p1 = plate(m, 4, 4);
                    let r1 = rib(m, 17.5);
                    m.rebuild_adjacency();
                    let ab = boolean(m, BoolKind::Fuse, p1, r1).expect("plate + rib")[0];
                    let r2 = rib(m, -17.5);
                    m.rebuild_adjacency();
                    let abc = boolean(m, BoolKind::Fuse, ab, r2).expect("+ the second rib")[0];
                    let c = prism(
                        m,
                        Axis::Z,
                        vec![Edge2d::circle(p2(0.0, 0.0), 10.0).unwrap()],
                        90.0,
                    );
                    let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
                    let turned = xf(
                        m,
                        c,
                        Isometry::rotation(nacre_scalar::Rotation {
                            axis: Axis::Y,
                            point: [r(0.0), r(0.0), r(0.0)],
                            angle: nacre_scalar::Angle::from_deg(Rat::from_int(90))
                                .expect("a right angle"),
                        }),
                    );
                    (abc, shift(m, turned, [45.0, 0.0, 47.0]))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("roundplate {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
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
