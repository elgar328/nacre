/// The shared solid fixtures (boxes, L/U prisms, pockets), one copy for unit and integration tests.
#[path = "../../../tests/support/fixtures.rs"]
mod fixtures;
/// The fixture vocabulary, shared with the integration tests (one copy).
#[path = "../../../tests/support/stated.rs"]
mod stated;
use crate::draft::LocalFace;
use fixtures::{
    boolean_one, cube_and_notch, extrude_op, l_and_corner_box, l_and_dimple, l_and_inner_box,
    l_and_popup_box, l_and_reflex_box, l_and_rod, l_prism, nested_boxes, outer_points, p2,
    pocket_op, regular_ngon, rotated_l_prism, small_square, square, stacked_cubes, two_boxes,
    u_and_slab, u_prism,
};
// Reached from outside this module (`bands`' tests) as `crate::tests::extrude_log_op`.
pub(crate) use fixtures::extrude_log_op;
use stated::{Stated, arc_turns, circle, line, stated};

use super::*;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: crate::SketchPlane) -> crate::SketchFrame {
    match crate::apply(
        m,
        &crate::Operation::DatumPlane {
            def: crate::DatumDef::Stated(plane),
        },
    ) {
        Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}
use crate::arrangement::{PlaneSetup, plane_index_setup};
use crate::combinatorics::{Canon3, NodeId};
use crate::tolerant::Judge;
use crate::transform::transform;
use crate::{assembly::*, ops::*, planes::*};
use nacre_exact::Axis;
use nacre_geom::Plane;
use nacre_geom::intersect::{planes_coplanar, three_planes};
use nacre_judge::WitnessPoint;
use nacre_topo::{Loop, Orientation, Surface, Vertex};
use proptest::prelude::*;
use std::collections::HashMap;

/// ★ **A fixture with no cylinders, said as a fact rather than left as a hole.** The table is how
/// a `NodeId::Pierce` reaches its definition, so an all-plane fixture has nothing to put in it —
/// and a bare `&[]` at a call site reads like something forgotten.
const NO_CYLS: &[crate::planes::WorkingCyl] = &[];

/// Is there an outer-shell face on the plane through `pt` with normal `n`, oriented that way?
/// The production path names a cap by the *face* that made it (`find_face_coplanar_with`); a
/// test that wants to say "a face sits on z = 1.5 facing +z" has no such face in hand, and
/// asserting geometry from coordinates is exactly what a test may do.
fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
    let Some(target) = Plane::from_point_normal(pt, n) else {
        return false;
    };
    let shell = m.solid(solid).outer;
    m.shell(shell).faces.iter().any(|&fh| {
        let f = m.face(fh);
        let nacre_geom::Surface::Plane(plane) = m.surface_cache(f.surface) else {
            return false;
        };
        let sign = f64::from(f.orientation.sign());
        planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n) > 0.0
    })
}

fn near(a: Point3, b: [f64; 3]) -> bool {
    (a - Point3::from_array(b)).norm() < 1e-9
}

/// The L-prism with an **L-shaped** stub standing wholly inside its top face,
/// `z ∈ [0.5, 1.5]`. The seam on the cap is a closed loop with a reflex node — the
/// suite's first non-convex inner loop, and the shape a winding must be read from.
///
/// Its coordinates dodge the cap's fan diagonals from `(0,0)` (`y = x`, `y = x/2`,
/// `y = 2x`), which `segment_crosses_face` would graze.
fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let ell = Profile2d::polygon(vec![
        p2(0.2, 0.25),
        p2(0.85, 0.25),
        p2(0.85, 0.4),
        p2(0.35, 0.4), // reflex
        p2(0.35, 0.9),
        p2(0.2, 0.9),
    ])
    .unwrap();
    let __f103 = datum_frame(
        &mut m,
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5])),
    );
    let OpOutput::Extrude { solid: stub, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f103,
            profile: ell,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, stub)
}

/// One face ring, as plane triples.
type Ring = Vec<[usize; 3]>;

/// A ring of vertex names as the plane triples these fixtures assert about.
fn names(ring: &[combinatorics::NodeId]) -> Ring {
    ring.iter()
        .map(|&n| combinatorics::three_plane_name(n).expect("a three-plane node"))
        .collect()
}

/// A **holed reflex face** from the live engine, as plane triples: `(planes, p, outer, hole)`.
///
/// `Cut(L-prism, stub)` leaves the L's top cap carrying a hole — `"dimple"` a square one,
/// `"ell"` an L-shaped one. Both rings come from [`combinatorics::face_vertex_triples`] and
/// [`combinatorics::hole_rings`], which the boolean itself uses, so the fixture exercises only code
/// the kernel runs.
///
/// This replaces a helper that built its rings from `seam_paths_on`/`orient_seam_loop` — the
/// retired seam engine. The *properties* below are about `point_in_ring`/`every_ray`, which are
/// live and load-bearing (`nest_cells` picks a hole's host with them, `unify_coplanar_faces`
/// groups by them), so they had to be re-homed rather than deleted with their old fixture.
fn holed_face_rings(which: &str) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let (m, l, stub) = if which == "dimple" {
        l_and_dimple()
    } else {
        l_and_ell_stub()
    };
    let (planes, p, outer, hole) = holed_face_rings_of(m, l, stub);
    assert_eq!(outer.len(), 6, "{which}: the L's cap is a reflex hexagon");
    (planes, p, outer, hole)
}

/// Cut `stub` out of `l` and hand back the first holed face's two rings, with the plane table
/// they are named in. The ring-length assertion belongs to the caller: a fixture built to put
/// a specific corner on a specific line **must** say how many corners it expects, because
/// `Profile2d` dissolves a vertex that sits mid-run on a straight edge and a silently
/// dissolved corner is a fixture that measures nothing.
fn holed_face_rings_of(
    mut m: Model,
    l: Handle<Solid>,
    stub: Handle<Solid>,
) -> (Vec<WorkingPlane>, usize, Ring, Ring) {
    let r = boolean_one(&mut m, BoolKind::Cut, l, stub).expect("the cut");
    m.rebuild_adjacency();
    let faces_tab = collect_planes(&m, r).unwrap();
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        surf_ix.insert(pi.face().expect("a real face table"), i);
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, _cyls) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).unwrap();
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        let holes = combinatorics::hole_rings(
            &m,
            fh,
            fp,
            &inc,
            &crate::planes::test_judge(&planes),
            &plane_ix,
            &[],
        )
        .unwrap();
        if let Some(hole) = holes.into_iter().next() {
            let outer = combinatorics::face_vertex_triples(
                &m,
                fh,
                fp,
                &inc,
                &crate::planes::test_judge(&planes),
                &plane_ix,
                &[],
            )
            .unwrap();
            return (
                planes,
                plane_ix[fp].plane(),
                names(&outer.poly().expect("a poly outer").triples),
                names(&hole.poly().expect("a poly hole").triples),
            );
        }
    }
    panic!("no holed face");
}

/// The unit cube with a 0.4-square pocket, 0.5 deep, in its top face: the void
/// is `[0.3,0.7]² × [0.5,1]` and the solid measures `1 − 0.16·0.5 = 0.92`. Its
/// lid is the only face in the suite that carries an inner loop.
fn pocketed_cube() -> (Model, Handle<Solid>) {
    let (mut m, top) = cube_with_top();
    let OpOutput::PocketOnFace { solid, .. } =
        apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
    else {
        unreachable!()
    };
    (m, solid)
}

/// Extrude a unit cube and return `(model, top face handle)`.
fn cube_with_top() -> (Model, Handle<Face>) {
    let mut m = Model::new();
    let op = extrude_op(&m, square(), 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
        unreachable!()
    };
    let top = faces[1]; // base, top, sides…
    (m, top)
}

/// An exact world-frame `Swept` from decimal f64 points — the truth-stating successor of the
/// retired `Swept::along` (S6b): the fallback is gone, so a test states its rings the way a
/// producer does. The f64 base is kept as handed in (decimals realize back bit-identically);
/// the top is the realization of the exact sum.
fn swept_world(base: Vec<Point3>, sweep: Vector3) -> crate::exact::Swept {
    let lift = |p: [f64; 3]| p.map(|x| nacre_exact::Rat::from_decimal(x).expect("decimal fixture"));
    let sv = lift(sweep.as_array());
    let rb: Vec<[nacre_exact::Rat; 3]> = base.iter().map(|p| lift(p.as_array())).collect();
    let rt: Vec<[nacre_exact::Rat; 3]> = rb
        .iter()
        .map(|b| core::array::from_fn(|i| b[i].checked_add(sv[i]).expect("fixture widths")))
        .collect();
    let top = crate::exact::realize(&rt);
    // The fixture's frame normal: the sweep's direction, which these fixtures keep axis-aligned.
    let normal: [nacre_exact::Rat; 3] = {
        let len2 = sv
            .iter()
            .map(|c| c.checked_mul(*c).unwrap())
            .fold(nacre_exact::Rat::from_int(0), |a, b| {
                a.checked_add(b).unwrap()
            });
        let inv = nacre_exact::inv_sqrt_exact(len2).expect("an axis-aligned fixture sweep");
        sv.map(|c| c.checked_mul(inv).unwrap())
    };
    let n = rb.len();
    // The fixture's winding, read exactly on its own (small) coordinates.
    let winding = {
        let cross =
            |a: &[nacre_exact::Rat; 3], b: &[nacre_exact::Rat; 3]| -> [nacre_exact::Rat; 3] {
                let t = |i: usize, j: usize| {
                    a[i].checked_mul(b[j])
                        .unwrap()
                        .checked_sub(a[j].checked_mul(b[i]).unwrap())
                        .unwrap()
                };
                [t(1, 2), t(2, 0), t(0, 1)]
            };
        let mut acc = nacre_exact::Rat::from_int(0);
        for i in 0..n {
            let c = cross(&rb[i], &rb[(i + 1) % n]);
            for k in 0..3 {
                acc = acc
                    .checked_add(c[k].checked_mul(normal[k]).unwrap())
                    .unwrap();
            }
        }
        match acc.cmp(&nacre_exact::Rat::from_int(0)) {
            core::cmp::Ordering::Greater => nacre_exact::Orient::Positive,
            core::cmp::Ordering::Less => nacre_exact::Orient::Negative,
            core::cmp::Ordering::Equal => nacre_exact::Orient::Zero,
        }
    };
    crate::exact::Swept {
        base,
        top,
        normal: Vector3::from_array(normal.map(|c| c.to_f64())),
        exact: crate::exact::SweptRat {
            base: rb,
            top: rt,
            segs: vec![crate::exact::Seg3::Line; n],
            normal,
            winding,
            motion: None,
        },
    }
}

/// The reported four-plane model: a unit cube with a block fused on its top, and a bar spun
/// `deg`° about an axis in the block's `x = 0.5` plane. At 45° with `half_z == 0.2` the bar's
/// half-width and its pivot-to-bottom offset are equal, so its bottom corner edge lands in that
/// plane and the y-planes cutting the edge become four-plane vertices.
fn four_plane_model(half_z: f64, deg: i128) -> (Model, Handle<Solid>, Handle<Solid>) {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let mut m = Model::new();
    let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let block = m.add_cuboid(
        Point3::from_array([0.5, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    m.rebuild_adjacency();
    let target = boolean(&mut m, BoolKind::Fuse, cube, block).expect("the block fuses on")[0];
    m.rebuild_adjacency();
    let bar = m.add_cuboid(
        Point3::from_array([0.3, -0.5, 1.0 - half_z]),
        Point3::from_array([0.7, 1.5, 1.0 + half_z]),
    );
    m.rebuild_adjacency();
    let bar = transform(
        &mut m,
        bar,
        &Isometry::rotation(Rotation {
            axis: Axis::Y,
            pivot: [Rat::new(1, 2).unwrap(), Rat::from_int(0), Rat::from_int(1)],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        }),
    )
    .unwrap();
    m.rebuild_adjacency();
    (m, target, bar)
}

/// Vertices the kernel vouches for beyond the construction's bare figure: realized from their
/// definition (`Bounded`) or carrying a measured residual. A boolean's vertices realize, so the
/// count is the population under the name the cache gives it.
fn count_discovered(m: &Model, s: Handle<Solid>) -> usize {
    let mut seen = std::collections::HashSet::new();
    let mut n = 0;
    let sh = m.solid(s).outer;
    for &fh in &m.shell(sh).faces {
        for he in &m.face(fh).outer.half_edges {
            {
                for vh in m.edge(he.edge).vertices {
                    // ★ Positive: "realized from its definition" is the claim, and the negative
                    // form would let `Ceiling` in — a vertex the cache road stopped on.
                    if seen.insert(vh)
                        && matches!(m.vertex_cache(vh), nacre_topo::PointCache::Bounded { .. })
                    {
                        n += 1;
                    }
                }
            }
        }
    }
    n
}

/// Does any surface of this solid record a motion? **The premise every fixture below asserts
/// first**: a translation whose `f64` landing happens to be exact records nothing, and a fixture
/// that silently took that road would prove the opposite of what it claims (measured — a boss
/// fixture did exactly that while the road under test was still shut).
fn carries_motion(m: &Model, s: Handle<Solid>) -> bool {
    crate::planes::solid_shell_handles(m, s)
        .into_iter()
        .flat_map(|sh| m.shell(sh).faces.clone())
        .any(|fh| m.plane_motion(m.face(fh).surface).is_some())
}

/// ★★★★★ **The eye that was missing: does the mesh actually cover the face it approximates?**
///
/// A boolean result was once shipped whose lateral mesh lost **16% of its area** — four triangles
/// spanned half a turn of the cylinder as flat chords, and the boss rendered as two cones. Every
/// standing instrument was green at the time: `validate` clean, the mesh watertight, the **exact**
/// volume right (`mass_props` integrates the b-rep, not the mesh), the face counts as predicted.
/// Nothing asked the one question that mattered — whether the triangles lie on the surface they
/// claim — so nothing answered it.
///
/// Two readings, and **neither is derived from the other**: a face's triangles against
/// [`nacre_props::face_props`] (an exact boundary integral) and the whole solid's triangles against
/// [`nacre_props::mass_props`] (the divergence theorem on the b-rep).
///
/// The mesh-covers-the-faces oracle, shared: a face's triangles against
/// [`nacre_props::face_props`] and the solid's against [`nacre_props::mass_props`], at the derived
/// relative budget (see [`the_mesh_covers_the_faces_it_approximates`]). Any test that meshes a
/// solid for the first time calls this — it is the only reading in the workspace that asks
/// whether the triangles lie on the face they claim.
pub(crate) fn mesh_covers_faces(name: &str, m: &Model, solids: &[Handle<Solid>]) {
    // 5 × the derived worst case (2.0e-4, a full disk at the angular budget).
    const BUDGET: f64 = 1e-3;
    let mesh = nacre_tess::tessellate(m, &nacre_tess::TessConfig::default()).expect("tess");
    for &s in solids {
        let sol = m.solid(s).clone();
        let mut mesh_volume = 0.0f64;
        for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
            for fh in m.shell(sh).faces.clone() {
                let exact = nacre_props::face_props(m, fh)
                    .unwrap_or_else(|e| panic!("{name}: face_props {e:?}"))
                    .area;
                let mut area = 0.0f64;
                for &th in mesh.by_face.get(&fh).map(|v| v.as_slice()).unwrap_or(&[]) {
                    let t = mesh.triangles.get(th);
                    let p: Vec<_> = t
                        .vertices
                        .iter()
                        .map(|&h| mesh.vertices.get(h).pos)
                        .collect();
                    area += (p[1] - p[0]).cross(p[2] - p[0]).norm() * 0.5;
                    // The signed volume of the tetrahedron on the origin; summed over an
                    // outward-oriented closed mesh it is the volume that mesh encloses.
                    let (a, b, c) = (p[0].as_array(), p[1].as_array(), p[2].as_array());
                    mesh_volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
                        - a[1] * (b[0] * c[2] - b[2] * c[0])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]))
                        / 6.0;
                }
                assert!(
                    (area - exact).abs() <= BUDGET * exact,
                    "{name}: a face's mesh area {area} is not its own {exact}"
                );
            }
        }
        let exact = nacre_props::mass_props(m, s).expect("mass_props").volume;
        assert!(
            (mesh_volume - exact).abs() <= BUDGET * exact,
            "{name}: the mesh encloses {mesh_volume}, the solid is {exact}"
        );
    }
}

/// A plate, one cylinder op, then a second — the population the corpus did not have.
pub(crate) fn chained(
    first: ([f64; 3], f64, BoolKind),
    second: ([f64; 3], f64, BoolKind),
) -> (Model, Result<Vec<Handle<Solid>>, BoolError>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let up = Vector3::from_array([0.0, 0.0, 1.0]);
    let a = m.add_cylinder(Point3::from_array(first.0), up, 0.5, first.1);
    m.rebuild_adjacency();
    let out = boolean(&mut m, first.2, plate, a).expect("the first op builds");
    assert_eq!(out.len(), 1, "the first op is one solid");
    m.rebuild_adjacency();
    let b = m.add_cylinder(Point3::from_array(second.0), up, 0.5, second.1);
    m.rebuild_adjacency();
    let r = boolean(&mut m, second.2, out[0], b);
    (m, r)
}

fn pinned_ends_ordered(at: [f64; 3], dir: [f64; 3], kind: BoolKind) -> usize {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([12.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(Point3::from_array(at), Vector3::from_array(dir), 0.5, 4.0);
    m.rebuild_adjacency();
    let out = boolean(&mut m, kind, plate, boss).expect("the boss builds");
    m.rebuild_adjacency();
    let r = out[0];
    let faces_tab = collect_planes(&m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(&m, r, &surf_ix).expect("edge incidence");
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = crate::planes::world_cylinder_def(&m, surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(surf) else {
                unreachable!()
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                realized: *cache,
                owner: crate::planes::SolidSide::A,
            }
        })
        .collect();
    let jd = crate::planes::test_judge(&planes);
    let mut pins = 0usize;
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        let crate::planes::ClassIx::Plane(p) = plane_ix[fp] else {
            continue;
        };
        let Ok(ring) = combinatorics::face_vertex_triples(&m, fh, fp, &inc, &jd, &plane_ix, &cyls)
        else {
            continue;
        };
        let Some(nr) = ring.poly() else { continue };
        let n = nr.triples.len();
        for i in 0..n {
            let crate::combinatorics::Wall::Plane(q) = nr.walls[i] else {
                continue;
            };
            let (a, b) = (nr.triples[i], nr.triples[(i + 1) % n]);
            let pin = |x: combinatorics::NodeId| match combinatorics::three_plane_name(x) {
                Some(t) => t
                    .iter()
                    .copied()
                    .find(|&c| c != p && c != q)
                    .map(combinatorics::EndPin::Class),
                None => Some(combinatorics::EndPin::Cylinder),
            };
            let (Some(pa), Some(pb)) = (pin(a), pin(b)) else {
                continue;
            };
            if matches!(pa, combinatorics::EndPin::Class(_))
                && matches!(pb, combinatorics::EndPin::Class(_))
            {
                continue;
            }
            pins += 1;
            // The f64 road: realize both points and dot the difference with `n_p × n_q` taken
            // from the raw coefficients — the same direction `plane_pair_dir_sign` reads.
            let xyz = |x: combinatorics::NodeId| -> [f64; 3] {
                match combinatorics::pierce_name(x) {
                    Some((_, cyl, _)) => {
                        combinatorics::pierce_point(&jd, cyl, &cyls[cyl].def, x).expect("realizes")
                    }
                    None => {
                        let c = combinatorics::node_coords_rat(&jd, x).expect("coords");
                        [c[0].to_f64(), c[1].to_f64(), c[2].to_f64()]
                    }
                }
            };
            let nv = |c: usize| {
                let k = jd.planes[c].plane.coefficients();
                Vector3::from_array([k[0], k[1], k[2]])
            };
            let d = nv(p).cross(nv(q));
            let (xa, xb) = (xyz(a), xyz(b));
            let t: f64 = (0..3).map(|k| (xa[k] - xb[k]) * d.as_array()[k]).sum();
            assert!(t.abs() > 1e-6, "the oracle decides in the noise: {t}");
            let want = if t > 0.0 { 1i8 } else { -1 };
            for (x, y, w) in [((a, pa), (b, pb), want), ((b, pb), (a, pa), -want)] {
                assert_eq!(
                    combinatorics::order_pinned(&jd, &cyls, p, q, x, y),
                    Some(w),
                    "p={p} q={q} {x:?} vs {y:?}"
                );
                // ★ `edge_dir` inverts the order to get the travel sense, and it is the carrier
                // it pairs with — not a bare `i8` — that the walk consumes.
                let dir = combinatorics::edge_dir(&jd, &cyls, p, q, x, y)
                    .unwrap_or_else(|e| panic!("p={p} q={q}: {e:?}"));
                assert!(
                    matches!(dir, combinatorics::EdgeDir::Line { carrier, sense }
                             if carrier == q && sense == -w),
                    "p={p} q={q}: {dir:?}"
                );
            }
        }
    }
    pins
}

fn rot_iso(axis: nacre_exact::Axis, deg: i128) -> nacre_exact::Isometry {
    use nacre_exact::{Angle, Isometry, Rat, Rotation as SRot};
    Isometry::rotation(SRot {
        axis,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
    })
}

fn translate_iso(off: [i128; 3]) -> nacre_exact::Isometry {
    use nacre_exact::{Isometry, Rat};
    Isometry::translation([
        Rat::from_int(off[0]),
        Rat::from_int(off[1]),
        Rat::from_int(off[2]),
    ])
}

fn rigid_iso(axis: nacre_exact::Axis, deg: i128, off: [i128; 3]) -> nacre_exact::Isometry {
    use nacre_exact::{Angle, Isometry, Rat, Rotation as SRot};
    Isometry::rigid(
        SRot {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        },
        [
            Rat::from_int(off[0]),
            Rat::from_int(off[1]),
            Rat::from_int(off[2]),
        ],
    )
}

/// A `FaceInfo` for the `unify_coplanar_faces` tests, which read none of its geometry; the rest is a
/// valid-but-unreferenced dummy (`surf`/`face`/`plane` are never dereferenced there).
fn face(plane_idx: usize, nodes: Vec<NodeId>, inner: Vec<Vec<NodeId>>) -> LocalFace {
    LocalFace {
        surf: crate::planes::ClassIx::Plane(plane_idx),
        outer: crate::draft::Bound::Ring(crate::draft::Ring::from_clean_names(plane_idx, nodes)),
        inner: inner
            .into_iter()
            .map(|r| crate::draft::Bound::Ring(crate::draft::Ring::from_clean_names(plane_idx, r)))
            .collect(),
        flip: false,
    }
}

/// A pentagonal prism whose fourth wall is slanted, and that wall's handle. Footprint area 15
/// (a 4×4 square less the 1×2 triangle the slant cuts off), swept 3.
fn prism_with_a_slanted_wall() -> (Model, Handle<Face>) {
    let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let __w1 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w1,
            profile: Profile2d::polygon(vec![
                p(0.0, 0.0),
                p(4.0, 0.0),
                p(4.0, 2.0),
                p(3.0, 4.0),
                p(0.0, 4.0),
            ])
            .unwrap(),
            dist: 3.0,
        },
    )
    .expect("extrude") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let wall = *m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .find(|&&f| {
            crate::ops::face_plane(&m, f).is_ok_and(|sp| {
                let n = sp.normal().as_array();
                n[0].abs() > 0.1 && n[1].abs() > 0.1 && n[2].abs() < 1e-12
            })
        })
        .expect("a slanted wall");
    (m, wall)
}

// ── a cylinder in the operation log is a circle extruded ────────────────────────────────
//
// `Model::add_cylinder_exact` is a test convenience behind a feature; the road an application
// takes is `Operation::Extrude` of a one-circle profile — `cylinder()` is that sugar, and the
// primitive operation retired once the two were measured identical. What the frame buys is that
// the statement never leaves the rationals: the axis is the frame's unit normal, the seam its
// `+u`. The extrude's faces come in push order: bottom cap, top cap, then the lateral.

/// The lateral surface's exact truth in a model holding exactly one cylinder.
fn lone_cylinder_def(m: &Model) -> nacre_topo::CylinderDef {
    let mut found = None;
    let mut i = 0u32;
    while let Some(h) = m.face_handle_at(i) {
        i += 1;
        let s = m.face(h).surface;
        if let nacre_topo::Surface::Cylinder { def, .. } = m.surface(s) {
            found = Some(def.clone());
        }
    }
    found.expect("a cylinder solid has a lateral face")
}

fn cylinder_op(m: &Model, center: [f64; 2], radius: f64, dist: f64) -> Operation {
    Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile: circle_profile(center, radius),
        dist,
    }
}

/// **The boss corpus** — a 4×4×2 plate at the origin and an r = 0.5 boss, `(name, base, height)`,
/// one row per way a boss can sit on a plate: through it, standing on it, flush with its floor, on
/// each wall, at a corner, off the wall's middle, and the half-height variants whose lateral is a
/// chain. Shared by the re-operation census and the crossing census so the two read one corpus.
const BOSS_FAMILIES: [(&str, [f64; 3], f64); 17] = [
    ("through", [2.0, 2.0, -1.0], 4.0),
    ("on top", [2.0, 2.0, 2.0], 1.0),
    ("flush", [2.0, 2.0, 0.0], 3.0),
    ("wall -y", [2.0, 0.0, -1.0], 4.0),
    ("wall +y", [2.0, 4.0, -1.0], 4.0),
    ("wall -x", [0.0, 2.0, -1.0], 4.0),
    ("wall +x", [4.0, 2.0, -1.0], 4.0),
    ("corner", [4.0, 4.0, -1.0], 4.0),
    ("corner-lo", [0.0, 0.0, -1.0], 4.0),
    ("offmid", [4.0, 1.0, -1.0], 5.0),
    ("half wall", [2.0, 0.0, -1.0], 2.0),
    ("half wall, cap below", [2.0, 0.0, 1.0], 2.0),
    ("half +x", [4.0, 2.0, -1.0], 2.0),
    ("half +x, cap below", [4.0, 2.0, 1.0], 2.0),
    // ★ The offset wall: the boss's axis stands off the plate's wall `x = 4` by less
    // than `r`, so the wall crosses the lateral in two rulings that are not a diameter's ends.
    // Outside by 0.3 (rational rulings, `y = 2 ± 0.4`; the part inside the plate a 0.2-deep
    // segment), inside by 0.3, outside by 0.2 (irrational rulings).
    ("offset-out", [4.3, 2.0, -1.0], 4.0),
    ("offset-in", [3.7, 2.0, -1.0], 4.0),
    ("offset-irr", [4.2, 2.0, -1.0], 4.0),
];

/// The plate and boss of a [`BOSS_FAMILIES`] row, adjacency rebuilt.
fn boss_family(base: [f64; 3], h: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let boss = m.add_cylinder(
        Point3::from_array(base),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        h,
    );
    m.rebuild_adjacency();
    (m, plate, boss)
}

// ═══════════════════════════════════════════════════════════════════════════════════════════
// **The boolean commutes with rigid motion.**
//
// One oracle: build the operands and take the boolean; build them again, move **both** by the
// same isometry, take the boolean; the two answers must be one answer up to that motion. The
// corpus is the boss corpus (`BOSS_FAMILIES`), one boss enclosed by the plate (a cavity), one
// planar pair (the control that the oracle is not cylinder-specific) and three
// second-operation families: a bored plate under a wall tool — the lateral's material lies
// *outside* the cylinder there — and a fused through-boss parted by a slab, a two-solid result.
//
// ★★★ **The motion group is chosen by algebra, not by story.** Every label the arrangement
// writes reads some product of `(frame_sign, axis_up, κ, exactness)`, and the unmoved corpus
// never produces `frame_sign = −1` on a wall class: a face sits with `frame_sign = −1` exactly
// when it lies on a **seed plane** (`Model::new` plants x = 0, y = 0, z = 0 with cache direction
// −axis) with its outward along +axis. So one dyadic translation that puts the plate's max faces
// onto the seeds reaches that class without any rotation; the quadrantal rotations turn the axis
// to ±x/±y (the ⊥ road's `axis_up` and the arc extremum's axis); a rigid motion of a non-dyadic
// origin exercises the transport's exactness boundary; a non-dyadic translation the recorded
// path itself. The same translation takes the min-side families *off* the seeds (frame +1), so
// the two classes are complementary.
//
// ★ `KNOWN` is the ledger of cells the kernel does not commute on **today**, each named by the
// sentence of the assertion that catches it — the crossing census's discipline: a fix flips its
// rows to commuting and the count lock says how many are left. A known cell is **run**, and must
// panic at one of its named sites; commuting instead is red, and so is a panic somewhere else.

/// What «the same answer» can mean under a motion.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum MotionClass {
    /// A pivot-0 quadrantal rotation: a signed permutation of coordinates, exact on any `f64`,
    /// so the moved result's vertex bits equal the unmoved result's bits moved.
    Quadrantal,
    /// An exact rigid motion or dyadic translation: rational vertices stay bit-exact; a pierce
    /// vertex (`a + b√c`, irrational) rounds once more on the moved road.
    ExactRigid,
    /// A non-dyadic translation: the recorded path — volumes, counts and validity only.
    Recorded,
}

fn motion_group() -> Vec<(String, nacre_exact::Isometry, MotionClass)> {
    use MotionClass::*;
    use nacre_exact::{Isometry, Rat};
    let mut g = Vec::new();
    for (ax, an) in [(Axis::X, "x"), (Axis::Y, "y"), (Axis::Z, "z")] {
        for deg in [90i128, 180, 270] {
            g.push((format!("r{an}{deg}"), rot_iso(ax, deg), Quadrantal));
        }
    }
    for (ax, an) in [(Axis::X, "x"), (Axis::Y, "y"), (Axis::Z, "z")] {
        g.push((
            format!("r{an}90+t"),
            rigid_iso(ax, 90, [5, -3, 2]),
            ExactRigid,
        ));
    }
    g.push(("t(5,-3,2)".into(), translate_iso([5, -3, 2]), ExactRigid));
    g.push((
        "t(-4,-4,-2)".into(),
        translate_iso([-4, -4, -2]),
        ExactRigid,
    ));
    g.push((
        "t(7/11,3/10,1/4)".into(),
        Isometry::translation([
            Rat::new(7, 11).unwrap(),
            Rat::new(3, 10).unwrap(),
            Rat::new(1, 4).unwrap(),
        ]),
        Recorded,
    ));
    g
}

type Build = Box<dyn Fn() -> (Model, Handle<Solid>, Handle<Solid>)>;

/// One boolean's answer, as the oracle compares it: the outcome's name, the sorted volumes, the
/// per-solid (faces, edges, vertices, cavities) sorted, and validity.
#[derive(Debug, Clone, PartialEq)]
struct Answer {
    name: String,
    volumes: Vec<f64>,
    counts: Vec<(usize, usize, usize, usize)>,
    valid: bool,
}

fn answer(m: &Model, r: &Result<Vec<Handle<Solid>>, BoolError>) -> Answer {
    let name = match r {
        Ok(v) => format!("Ok({})", v.len()),
        Err(BoolError::Rejected { reason, .. }) => format!("Rejected({reason:?})"),
        Err(e) => format!("Err({e:?})"),
    };
    let (mut volumes, mut counts) = (Vec::new(), Vec::new());
    if let Ok(v) = r {
        for &s in v {
            volumes.push(
                nacre_props::mass_props(m, s)
                    .map(|p| p.volume)
                    .unwrap_or(f64::NAN),
            );
            let sol = m.solid(s);
            let (mut faces, mut edges, mut verts) = (
                0usize,
                std::collections::HashSet::new(),
                std::collections::HashSet::new(),
            );
            for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
                for &fh in &m.shell(sh).faces {
                    faces += 1;
                    let face = m.face(fh);
                    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                        for he in &lp.half_edges {
                            edges.insert(he.edge);
                            for &vh in m.edge(he.edge).vertices.iter() {
                                verts.insert(vh);
                            }
                        }
                    }
                }
            }
            counts.push((faces, edges.len(), verts.len(), sol.cavities.len()));
        }
    }
    volumes.sort_by(f64::total_cmp);
    counts.sort_unstable();
    Answer {
        name,
        volumes,
        counts,
        valid: nacre_validate::validate(m).is_empty(),
    }
}

enum Outcome {
    Commutes,
    Diverged(String),
    Panicked(String),
    /// The moved *input* could not be built: `transform` refused an operand.
    InputUntransportable(OpError),
}

impl Outcome {
    fn text(&self) -> String {
        match self {
            Outcome::Commutes => "commutes".into(),
            Outcome::Diverged(why) => format!("diverged: {why}"),
            Outcome::Panicked(msg) => format!("panicked: {}", msg.lines().next().unwrap_or("")),
            Outcome::InputUntransportable(e) => format!("input untransportable: {e:?}"),
        }
    }
}

fn circle_profile(center: [f64; 2], radius: f64) -> Profile2d {
    stated(vec![circle(nacre_math::Point2::from_array(center), radius)])
        .unwrap()
        .remove(0)
}

/// The four-plane operand: a wall with a gusset whose apex lands on the wall's top edge, fused —
/// and the second gusset that will be fused onto it (the user's `1.4`/`−2.4`).
///
/// ★ It is **not** a rigid-motion oracle family, and the reason is measured: that oracle's digest
/// is bit-exact, and a vertex realized from three planes is exact only for axis-aligned triples.
/// With the slant `5y − 3z + 3 = 0` the moved-then-fused gusset tips came out at `z = −3.9e−16`
/// against `0`; with a 45° slant (`y − z + 1 = 0`) and dyadic positions every rotation agreed and
/// one translation still differed by one ULP in `x`. Same vertices, same names, different last
/// bits — so «one name under every motion» is asked of the audit instead
/// (`a_four_plane_operand_vertex_has_one_name_under_rigid_motion`).
fn wall_and_gusset_operand() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (x1, x2, top) = (1.4, -2.4, 6.0);
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist: 1.0,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let mut m = Model::new();
    let w = prism(
        &mut m,
        Axis::Y,
        vec![
            line(p2(1.0, -3.5), p2(top, -3.5)),
            line(p2(top, -3.5), p2(top, 3.5)),
            line(p2(top, 3.5), p2(1.0, 3.5)),
            line(p2(1.0, 3.5), p2(1.0, -3.5)),
        ],
    );
    let w = shift(&mut m, w, [0.0, 3.0, 0.0]);
    let gusset = |m: &mut Model, x: f64| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.0, 1.0), p2(3.0, 1.0)),
                line(p2(3.0, 1.0), p2(3.0, top)),
                line(p2(3.0, top), p2(0.0, 1.0)),
            ],
        );
        shift(m, g, [x, 0.0, 0.0])
    };
    let g1 = gusset(&mut m, x1);
    let g2 = gusset(&mut m, x2);
    m.rebuild_adjacency();
    let wg = boolean(&mut m, BoolKind::Fuse, w, g1).expect("wall + gusset")[0];
    m.rebuild_adjacency();
    (m, wg, g2)
}

/// The cell-13 operand: a 90 × 50 × 12 plate whose corners are filleted (r 5) and which carries
/// four bores (d 7). With every corner rounded its cap ring has **no three-plane corner at all** —
/// eight fillet tangencies and nothing else.
fn rounded_plate(m: &mut Model, fillets: usize, bores: usize) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let corner = [
        ([45.0, -25.0], [1.0, 0.0], [0.0, 1.0]),
        ([45.0, 25.0], [0.0, 1.0], [-1.0, 0.0]),
        ([-45.0, 25.0], [-1.0, 0.0], [0.0, -1.0]),
        ([-45.0, -25.0], [0.0, -1.0], [1.0, 0.0]),
    ];
    let r = 5.0;
    let mut edges: Vec<Stated> = Vec::new();
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
            edges.push(line(p2(at[0], at[1]), p2(tin[0], tin[1])));
            edges.push(arc_turns(p2(centre[0], centre[1]), p2(tin[0], tin[1]), 1));
            at = [c[0] + r * d_out[0], c[1] + r * d_out[1]];
        } else {
            edges.push(line(p2(at[0], at[1]), p2(c[0], c[1])));
            at = c;
        }
    }
    for &[x, y] in [[38.0, 18.0], [38.0, -18.0], [-38.0, 18.0], [-38.0, -18.0]]
        .iter()
        .take(bores)
    {
        edges.push(stated::circle(p2(x, y), 3.5));
    }
    let profile = stated(edges).unwrap().remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 12.0,
        },
    )
    .expect("the rounded plate extrudes") else {
        unreachable!()
    };
    solid
}

/// A cylinder of radius `r` and length `h` about the world Z axis, its base at the origin, then
/// turned a quarter about Y and translated — the shape a script's `cylinder()` states.
fn turned_cylinder(m: &mut Model, r: f64, h: f64, shift: [f64; 3]) -> Handle<Solid> {
    let profile = circle_profile([0.0, 0.0], r);
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("the cylinder extrudes") else {
        unreachable!()
    };
    let q = |x: f64| nacre_exact::Rat::from_decimal(x).expect("a short decimal");
    let turned = transform(
        m,
        solid,
        &nacre_exact::Isometry::rotation(nacre_exact::Rotation {
            axis: Axis::Y,
            pivot: [q(0.0), q(0.0), q(0.0)],
            angle: nacre_exact::Angle::from_deg(nacre_exact::Rat::from_int(90))
                .expect("a right angle"),
        }),
    )
    .expect("the cylinder turns");
    transform(
        m,
        turned,
        &nacre_exact::Isometry::translation([q(shift[0]), q(shift[1]), q(shift[2])]),
    )
    .expect("the cylinder moves")
}

mod common_kind;
mod commuting_oracle;
mod coplanar;
mod curved_nesting;
mod cylinder_classes;
mod cylinder_gate;
mod cylinder_ops;
mod extrude;
mod frames;
mod loop_nesting;
mod motion;
mod naming;
mod parallel;
mod rotation;
mod sketch_models;
mod user_models;
