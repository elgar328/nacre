//! **The bit census — a baseline to diff across a change that moves coordinates.**
//!
//! Prints one line per case: the exact bits of volume/area/centroid **and a hash of the result's
//! sorted vertex coordinates**. The derived quantities alone are not enough: a 1-ULP coordinate
//! change can cancel out of a volume integral, so a census of volumes can pass while the geometry
//! moved. The coordinate hash is what makes "nothing changed" checkable.
//!
//! Not an assertion — a dump. Run it on two commits and `diff`:
//!
//! ```text
//! cargo test -p nacre-ops --release --test census -- --ignored --nocapture | grep '^c '
//! ```

use nacre_math::Point3;
use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

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
/// (`Origin::Discovered { ThreePlane }`), so the result carries *no* trace of the operands' vertex
/// coordinates. A change that moves operand vertices but leaves the planes alone is therefore
/// invisible in the result — which is exactly what happened the first time this census was used.
/// Recording the operands is what makes the census see the change it exists to see.
fn plane_digest(m: &Model, s: Handle<Solid>) -> (usize, u64) {
    use std::hash::{Hash, Hasher};
    let mut bits: Vec<[u64; 4]> = Vec::new();
    let src = m.solids.get(s).clone();
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            match m.surfaces.get(m.faces.get(fh).surface) {
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
                    for &vh in m.edges.get(he.edge).bounds.iter().flatten() {
                        let p = m.vertices.get(vh).point.as_array();
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
                let c = nacre_props::centroid(m, s).expect("centroid");
                let (n, d) = coord_digest(m, s);
                parts.push(format!(
                    "{:016x}/{:016x}/{:016x},{:016x},{:016x}/v{n}h{d:016x}",
                    p.volume.to_bits(),
                    p.area.to_bits(),
                    c[0].to_bits(),
                    c[1].to_bits(),
                    c[2].to_bits(),
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
}

/// A rectangular prism raised from `z` by `dist` — the construction path production uses.
fn ex(m: &mut Model, z: f64, lo: [f64; 2], hi: [f64; 2], dist: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            plane: nacre_ops::SketchPlane {
                origin: Point3::from_array([0.0, 0.0, z]),
                ..nacre_ops::SketchPlane::world_xy()
            },
            profile: nacre_ops::Profile2d::polygon(vec![
                p(lo[0], lo[1]),
                p(hi[0], lo[1]),
                p(hi[0], hi[1]),
                p(lo[0], hi[1]),
            ]),
            dist,
        },
    )
    .expect("extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    solid
}
