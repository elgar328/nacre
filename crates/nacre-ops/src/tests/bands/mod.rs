use super::*;

use crate::arrangement::{PlaneSetup, plane_index_setup};

use nacre_math::{Point3, Vector3};

use nacre_store::Handle;

use nacre_topo::Surface;

use nacre_topo::{Model, Solid};

/// Run the plane arrangement past the stopper and hand the band pass what it needs.
/// Returns `(band faces, the plane classes' axis parameters keyed by class)`.
fn bands(
    m: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    kind: BoolKind,
) -> (Vec<LocalFace>, Vec<(usize, f64)>) {
    let setup = plane_index_setup(m, a, b).unwrap();
    let PlaneSetup {
        planes: faces_tab,
        geom,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        class_owner,
        n_a,
        standard,
        notes,
        cyls,
        ..
    } = &setup;
    let jd = Judge::new(geom, *standard, notes);
    let trace_in = crate::combinatorics::trace_input(
        m,
        [(a, inc_a), (b, inc_b)],
        surf_ix,
        faces_tab.len(),
        &jd,
        plane_ix,
        cyls,
        Default::default(),
    );
    let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
        m,
        kind,
        a,
        b,
        &jd,
        faces_tab,
        plane_ix,
        cyls,
        *n_a,
        class_owner,
        &trace_in,
    )
    .expect("the drill population traces");
    let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("cylinder rows");
    let out = crate::arrangement::cyl_chart::emit_lateral(
        kind,
        &jd,
        cyls,
        &plane_faces,
        &curved,
        &rows,
        &crate::draft::held_rims(&plane_faces, &curved.split_rims),
    )
    .expect("the lateral faces emit");
    // The classes' **z**, not their axis parameter: `t` is measured from the cylinder's own
    // origin along its raw `dir`, so a drill starting at z=−1 puts the box's cap at t=1. The
    // assertions read in world z, which is the vocabulary the fixtures are written in.
    let ts = (0..geom.len())
        .filter_map(|c| {
            let def = &cyls[0].def;
            let t = param_opt(&jd, c, def)?;
            // The class's world z, via the axis point at that parameter.
            let (o, m) = (def.origin(), def.dir());
            let z = o[2].checked_add(t.checked_mul(m[2])?)?;
            Some((c, z.to_f64()))
        })
        .collect();
    (out, ts)
}

/// A band's two ends as world `z` — the assertion vocabulary: its lower whole rim (the outer
/// bound) and its upper one (an inner bound).
fn ends(lf: &LocalFace, ts: &[(usize, f64)]) -> (f64, f64) {
    let Bound::Rim {
        plane: lo,
        ccw: true,
    } = lf.outer
    else {
        panic!(
            "a band's lower rim is a whole circle here, got {:?}",
            lf.outer
        );
    };
    let [hi] = lf
        .inner
        .iter()
        .filter_map(|b| match b {
            Bound::Rim { plane, ccw: false } => Some(*plane),
            _ => None,
        })
        .collect::<Vec<_>>()[..]
    else {
        panic!(
            "a band's upper rim is a whole circle here, got {:?}",
            lf.inner
        );
    };
    let at = |c: usize| ts.iter().find(|(k, _)| *k == c).expect("a ⊥ class").1;
    (at(lo), at(hi))
}

fn box_and_drill(z0: f64, h: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0; 3]),
    );
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([1.0, 1.0, z0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    (m, a, b)
}

/// `χ = V − E + F − L` over **one solid's** shells — the counts the genus relation
/// `χ = 2(S − G)` is read from.
///
/// ★ Per solid, not per model: `Model::reachable` spans every live solid, which coincides with
/// this only while a fixture makes exactly one. The fixtures that make two need the sum split,
/// and one spelling for all of them is what keeps the relation from being restated.
fn euler_counts(m: &Model, s: Handle<Solid>) -> (i64, i64, i64, i64) {
    let faces: Vec<_> = crate::planes::solid_shell_handles(m, s)
        .into_iter()
        .flat_map(|sh| m.shell(sh).faces.clone())
        .collect();
    let mut verts = std::collections::HashSet::new();
    let mut edges = std::collections::HashSet::new();
    let mut loops = 0i64;
    for &fh in &faces {
        let f = m.face(fh);
        loops += f.inner.len() as i64;
        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
            for he in &lp.half_edges {
                edges.insert(he.edge);
                verts.extend(m.edge(he.edge).vertices);
            }
        }
    }
    (
        verts.len() as i64,
        edges.len() as i64,
        faces.len() as i64,
        loops,
    )
}

/// How many faces of `s` lie on cylinder surfaces, grouped by surface.
fn lateral_face_counts(m: &Model, s: Handle<Solid>) -> Vec<usize> {
    let mut counts: std::collections::HashMap<Handle<Surface>, usize> = Default::default();
    let solid = m.solid(s);
    for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
        for &fh in &m.shell(sh).faces {
            let surf = m.face(fh).surface;
            if matches!(m.surface(surf), nacre_topo::Surface::Cylinder { .. }) {
                *counts.entry(surf).or_default() += 1;
            }
        }
    }
    let mut v: Vec<usize> = counts.into_values().collect();
    v.sort_unstable();
    v
}

mod blind_bores;
mod contacts_and_arcs;
mod lateral_faces;
mod seated_caps;
mod slab_bands;
mod wall_faces;
