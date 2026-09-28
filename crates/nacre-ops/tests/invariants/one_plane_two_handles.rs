//! **One plane held as two surface handles booleans as it does as one.**
//!
//! A mixed-frame `Through` datum has no name — no one frame solves its three vertices
//! rationally — so it cannot intern onto a plane the model already names, and the same geometric
//! plane is two handles. The class discovery does not read handles for that: it asks
//! `Judge::planes_coplanar` of every pair of distinct surfaces (an exact zero, or a coincidence
//! proved within the coincidence precision and reported in `BoolReport::merges`), so the two
//! faces land in one plane class. What that leaves to measure is whether the result is the one
//! the single-handle model gets, and that is what these fixtures assert.
//!
//! ★ **Geometry, not model contents.** A result face takes its class representative's surface —
//! the lowest face row, which is the first operand's side — so where the first operand sits on
//! the nameless handle the merged face keeps it, and the named twin's model does not. That is a
//! known defect of its own (todo), so the comparison reads positions, windings, volumes and
//! counts, and not which handle a face carries.

use crate::fixtures::{live_vertices, rot_iso};
use nacre_exact::Axis;
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolKind, BoolReport, DatumDef, OpOutput, Operation, Profile2d, SketchFrame, SketchPlane,
    apply, boolean_with_report,
};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, Surface, Vertex};

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Profile2d {
    let p = |x, y| Point2::from_array([x, y]);
    Profile2d::polygon(vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]).unwrap()
}

fn extrude(m: &mut Model, frame: SketchFrame, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("extrude") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

fn datum(m: &mut Model, def: DatumDef) -> (Handle<Surface>, SketchFrame) {
    let OpOutput::DatumPlane { plane, frame } =
        apply(m, &Operation::DatumPlane { def }).expect("datum")
    else {
        unreachable!()
    };
    (plane, frame)
}

/// The live vertex at `want` — fixtures pick vertices by where they are, never by walk order.
fn at(m: &Model, want: [f64; 3]) -> Handle<Vertex> {
    *live_vertices(m)
        .iter()
        .find(|&&v| {
            let p = m.vertex_point(v).as_array();
            (0..3).all(|k| (p[k] - want[k]).abs() < 1e-9)
        })
        .unwrap_or_else(|| panic!("no vertex at {want:?}"))
}

/// A world box `B = [−2,2]² × [0,1]` (or, `below`, moved down to `[−1,0]` — a translation, carried
/// into its statements, so B is world-stated either way and its `z = 0` corners are the same four)
/// and a box `A = [5,6] × [0,1] × [0,1]` turned 37° about `z`: A's floor corners meet in A's turned
/// frame (its floor is fixed by the turn), B's in the world, so a triple mixing them has no one
/// frame. Returns B and A's floor corners that were `(5,0)` and `(6,0)` — both have `y > 2`, so
/// either makes a triple with B's front edge span `+z`.
fn floor_pair(m: &mut Model, below: bool) -> (Handle<Solid>, Handle<Vertex>, Handle<Vertex>) {
    let world = SketchFrame::world(m, Axis::Z);
    let mut b = extrude(m, world, rect(-2.0, -2.0, 2.0, 2.0), 1.0);
    if below {
        let down = nacre_exact::Isometry::translation([
            nacre_exact::Rat::from_int(0),
            nacre_exact::Rat::from_int(0),
            nacre_exact::Rat::from_int(-1),
        ]);
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: b,
                isometry: down,
            },
        )
        .expect("move B down") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        b = solid;
    }
    let world = SketchFrame::world(m, Axis::Z);
    let a = extrude(m, world, rect(5.0, 0.0, 6.0, 1.0), 1.0);
    let OpOutput::Transform { .. } = apply(
        m,
        &Operation::Transform {
            solid: a,
            isometry: rot_iso(Axis::Z, 37),
        },
    )
    .expect("turn") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let (c, s) = (37f64.to_radians().cos(), 37f64.to_radians().sin());
    (
        b,
        at(m, [5.0 * c, 5.0 * s, 0.0]),
        at(m, [6.0 * c, 6.0 * s, 0.0]),
    )
}

/// One model of one scene. `nameless` builds the plane as the mixed-frame datum, otherwise as a
/// named statement of the same plane; `flush` puts the other operand on the far side of the plane,
/// so the two only touch there — the case that leans hardest on the merge.
struct Scene {
    m: Model,
    lhs: Handle<Solid>,
    rhs: Handle<Solid>,
    /// The two handles of the one plane in the nameless model (checked distinct), or the one
    /// handle in the named model.
    planes: Vec<Handle<Surface>>,
    /// The prism(s) built on the datum, before the boolean.
    prisms: Vec<Handle<Solid>>,
}

/// S1 — B's floor `z = 0` stated again as a datum through two of B's floor corners and one of
/// A's; a prism on that datum against B. The merge is an exact zero: every point is on `z = 0`
/// exactly, since a turn about `z` keeps `z`.
fn s1(nameless: bool, flush: bool) -> Scene {
    let mut m = Model::new();
    let (b, a0, _) = floor_pair(&mut m, flush);
    let (p00, p10, p11) = (
        at(&m, [-2.0, -2.0, 0.0]),
        at(&m, [2.0, -2.0, 0.0]),
        at(&m, [2.0, 2.0, 0.0]),
    );
    let third = if nameless { a0 } else { p11 };
    let (t, frame) = datum(&mut m, DatumDef::ThroughVertices([p00, p10, third]));
    // B's floor: the named `z = 0` the datum states again.
    let floor = m
        .shell(m.solid(b).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .find(|h| {
            m.surface_name
                .get(h)
                .and_then(|n| n.narrow())
                .is_some_and(|c| c.map(|r| r.to_f64()) == [0.0, 0.0, 1.0, 0.0])
        })
        .expect("B's floor names z = 0");
    let p = extrude(&mut m, frame, rect(-1.0, -1.0, 1.0, 1.0), 0.5);
    let planes = if nameless { vec![t, floor] } else { vec![t] };
    Scene {
        m,
        lhs: b,
        rhs: p,
        planes,
        prisms: vec![p],
    }
}

/// S2 — the same `z = 0` stated as **two** mixed-frame datums through different triples; a prism
/// on each, one against the other. Two nameless handles: the merge is a coincidence proved
/// within the coincidence precision (A's turned corners are irrational in `x`, `y`).
fn s2(nameless: bool, _flush: bool) -> Scene {
    let mut m = Model::new();
    let (_, a0, a1) = floor_pair(&mut m, false);
    let (p00, p10, p11, p01) = (
        at(&m, [-2.0, -2.0, 0.0]),
        at(&m, [2.0, -2.0, 0.0]),
        at(&m, [2.0, 2.0, 0.0]),
        at(&m, [-2.0, 2.0, 0.0]),
    );
    let first = if nameless {
        [p00, p10, a0]
    } else {
        [p00, p10, p11]
    };
    let second = if nameless {
        [p01, p11, a1]
    } else {
        [p10, p11, p01]
    };
    let (t1, f1) = datum(&mut m, DatumDef::ThroughVertices(first));
    let (t2, f2) = datum(&mut m, DatumDef::ThroughVertices(second));
    let p1 = extrude(&mut m, f1, rect(-1.0, -1.0, 1.0, 1.0), 0.5);
    let p2 = extrude(&mut m, f2, rect(0.0, 0.0, 1.5, 1.5), 0.5);
    let planes = vec![t1, t2];
    Scene {
        m,
        lhs: p1,
        rhs: p2,
        planes,
        prisms: vec![p1, p2],
    }
}

/// S3 — a tilted plane (normal `(1,1,1)`) named by a stated datum D, a box B raised on it, and the
/// same plane stated again through three of B's floor corners. Those corners straddle — D is
/// world-stated, B's walls hang off D's frame node — so the datum is nameless and its witness
/// is the judged frame's; the merge with D is a coincidence proved within the coincidence
/// precision. The corners are picked by position and ordered so their cross faces `(1,1,1)`, the
/// way D's frame faces. `flush` raises B on D's other side (D stated facing `−(1,1,1)`, the same
/// handle — both models build B the same way) so the prism only touches it.
fn s3(nameless: bool, flush: bool) -> Scene {
    let mut m = Model::new();
    let n = [1.0, 1.0, 1.0];
    let stated = |sign: f64| {
        DatumDef::Stated(
            SketchPlane::from_origin_normal(
                Point3::from_array([0.0; 3]),
                Vector3::from_array(n.map(|x| sign * x)),
            )
            .expect("a tilted plane"),
        )
    };
    let (d, fd) = datum(&mut m, stated(1.0));
    let b = if flush {
        let (d_back, back) = datum(&mut m, stated(-1.0));
        assert_eq!(d_back, d, "one plane, stated either way, is one handle");
        extrude(&mut m, back, rect(-2.0, -2.0, 2.0, 2.0), 1.0)
    } else {
        extrude(&mut m, fd, rect(-2.0, -2.0, 2.0, 2.0), 1.0)
    };
    let mut floor: Vec<Handle<Vertex>> = live_vertices(&m)
        .into_iter()
        .filter(|&v| matches!(*m.vertex(v), Vertex::ThreePlane(t) if t.contains(&d)))
        .collect();
    assert_eq!(floor.len(), 4, "B's floor has four corners on D");
    assert!(
        floor.iter().all(|&v| m.vertex_meet(v).is_none()),
        "the fixture's premise: B's floor corners straddle D's world and its walls' frame"
    );
    let key = |v: &Handle<Vertex>| m.vertex_point(*v).as_array().map(f64::to_bits);
    floor.sort_by_key(key);
    let mut tri = [floor[0], floor[1], floor[2]];
    let w = tri.map(|v| m.vertex_point(v));
    let cross = (w[1] - w[0]).cross(w[2] - w[0]).as_array();
    if (0..3).map(|k| cross[k] * n[k]).sum::<f64>() < 0.0 {
        tri.swap(0, 1);
    }
    let (t, frame) = if nameless {
        datum(&mut m, DatumDef::ThroughVertices(tri))
    } else {
        (d, fd)
    };
    let p = extrude(&mut m, frame, rect(-1.0, -1.0, 1.0, 1.0), 0.5);
    let planes = if nameless { vec![t, d] } else { vec![t] };
    Scene {
        m,
        lhs: b,
        rhs: p,
        planes,
        prisms: vec![p],
    }
}

/// Sorted vertex positions of `s`, as bits.
fn positions(m: &Model, s: Handle<Solid>) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = Vec::new();
    for &f in &m.shell(m.solid(s).outer).faces {
        let face = m.face(f);
        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
            for &he in &lp.half_edges {
                out.push(m.vertex_point(m.he_start(he)).as_array().map(f64::to_bits));
            }
        }
    }
    out.sort();
    out.dedup();
    out
}

/// What a result is, geometrically: per face, its loops as cyclic sequences of vertex positions
/// (bits) — the winding says which way the face points — plus each solid's volume (bits) and the
/// counts of distinct vertex and edge handles.
#[derive(Debug, PartialEq)]
struct Shape {
    faces: Vec<Vec<Vec<[u64; 3]>>>,
    volumes: Vec<u64>,
    vertices: usize,
    edges: usize,
}

fn shape(m: &Model, solids: &[Handle<Solid>]) -> Shape {
    let mut faces = Vec::new();
    let mut volumes = Vec::new();
    let (mut vs, mut es) = (
        std::collections::HashSet::new(),
        std::collections::HashSet::new(),
    );
    for &s in solids {
        volumes.push(
            nacre_props::mass_props(m, s)
                .expect("measurable")
                .volume
                .to_bits(),
        );
        let sol = m.solid(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &f in &m.shell(sh).faces {
                let face = m.face(f);
                let mut loops: Vec<Vec<[u64; 3]>> = std::iter::once(&face.outer)
                    .chain(face.inner.iter())
                    .map(|lp| {
                        let ring: Vec<[u64; 3]> = lp
                            .half_edges
                            .iter()
                            .map(|he| {
                                vs.insert(m.he_start(*he));
                                es.insert(he.edge);
                                m.vertex_point(m.he_start(*he)).as_array().map(f64::to_bits)
                            })
                            .collect();
                        // A cycle has no first element: start it at its least one.
                        let k = (0..ring.len()).min_by_key(|&i| ring[i]).unwrap_or(0);
                        ring[k..].iter().chain(&ring[..k]).copied().collect()
                    })
                    .collect();
                loops.sort();
                faces.push(loops);
            }
        }
    }
    faces.sort();
    volumes.sort();
    Shape {
        faces,
        volumes,
        vertices: vs.len(),
        edges: es.len(),
    }
}

fn run(sc: &mut Scene, kind: BoolKind) -> (Shape, BoolReport) {
    let (out, report) = boolean_with_report(&mut sc.m, kind, sc.lhs, sc.rhs)
        .unwrap_or_else(|e| panic!("{kind:?} must answer: {e:?}"));
    sc.m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&sc.m).is_empty(),
        "{kind:?}: the result is a valid b-rep"
    );
    (shape(&sc.m, &out), report)
}

/// ★★★★ **A plane stated twice — once without a name — booleans as it does stated once.** Three
/// scenes (an exact merge, two nameless handles merged on a proved coincidence, a tilted plane's
/// judged datum merged with its named twin), each `Fuse`/`Cut`/`Common`, each overlapping and (S1,
/// S3) in flush contact. The nameless model and the named one must give the same shape.
///
/// ★ S3's prisms are **not** compared by bits before the boolean: corners on a nameless carrier
/// have no realization road (`Unrealized`, the `NoMeet` population), so they hold the
/// construction's figure until the boolean mints the result's vertices from their definitions.
/// They agree within `1e-12` there, and the results agree to the bit.
#[test]
fn a_plane_stated_twice_booleans_as_it_does_once() {
    type Build = fn(bool, bool) -> Scene;
    let scenes: [(&str, Build, bool, &[bool]); 3] = [
        ("S1", s1, false, &[false, true]),
        ("S2", s2, true, &[false]),
        ("S3", s3, true, &[false, true]),
    ];
    for (name, build, proved, contacts) in scenes {
        for &flush in contacts {
            for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                let (mut two, mut one) = (build(true, flush), build(false, flush));
                let at = format!("{name} flush={flush} {kind:?}");

                // Premises: two handles of one plane on one side, one named handle on the other.
                assert!(
                    two.planes
                        .iter()
                        .any(|h| !two.m.surface_name.contains_key(h)),
                    "{at}: the datum is nameless"
                );
                assert_ne!(two.planes[0], two.planes[1], "{at}: two handles");
                assert!(
                    one.planes
                        .iter()
                        .all(|h| one.m.surface_name.contains_key(h)),
                    "{at}: the control's plane is named"
                );
                for (p2, p1) in two.prisms.iter().zip(&one.prisms) {
                    let (a, b) = (positions(&two.m, *p2), positions(&one.m, *p1));
                    if name == "S3" {
                        let near = a.len() == b.len()
                            && a.iter().zip(&b).all(|(x, y)| {
                                (0..3).all(|k| {
                                    (f64::from_bits(x[k]) - f64::from_bits(y[k])).abs() < 1e-12
                                })
                            });
                        assert!(near, "{at}: the two prisms stand in one place");
                    } else {
                        assert_eq!(a, b, "{at}: the two prisms are the same prism");
                    }
                }

                let (shape_two, report) = run(&mut two, kind);
                let (shape_one, _) = run(&mut one, kind);
                // The merge took the road the scene is for.
                if proved {
                    assert!(
                        !report.merges.is_empty(),
                        "{at}: the class merge rested on a proved coincidence"
                    );
                } else {
                    assert!(
                        report.merges.is_empty(),
                        "{at}: the class merge was exact (B is the first operand)"
                    );
                }
                if flush && kind == BoolKind::Common {
                    assert!(
                        shape_two.faces.is_empty(),
                        "{at}: the fixture's premise — flush operands only touch"
                    );
                }
                assert_eq!(
                    shape_two, shape_one,
                    "{at}: two handles of one plane gave a different result"
                );
            }
        }
    }
}
