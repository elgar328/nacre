//! Toleranced geometric predicates: the plane-arrangement sign predicates routed through
//! the CIP kernel so they stay exact under rotation.
//!
//! A rotated face's plane coefficients and `tri` coordinates are rounded irrationals, so the
//! exact predicates (`nacre-predicates`) over them are exact only w.r.t. the *rounded*
//! geometry. When a predicate's planes are rotated, these wrappers rebuild each plane from the
//! three exact [`Pt3`] its face carries ([`Witness::tri_pt3`], or — for the axis-aligned
//! operand of a *mixed*-rotation boolean, whose `tri_pt3` is `None` — exactly from its `tri`
//! coordinates, see [`plane_def`]) and decide the sign with the [`crate::kernel`] judges
//! instead.
//!
//! **Routing is per-predicate, derived — no `rotated` flag is threaded.** Each wrapper asks
//! [`any_rotated`] of just the planes it touches: all-axis-aligned → the exact hot path (never
//! builds a `Pt3`); any rotated → the kernel judge. So a mixed-rotation boolean keeps the
//! axis-aligned operand's own predicates on the fast path, finer than a per-boolean flag.
//!
//! The caller (`nacre-ops`) provides the witnesses by implementing [`Witness`] (a face) and
//! [`PlaneWitness`] (a plane class) on its own tables — the port keeps this crate free of the
//! b-rep and geometry types.

use crate::kernel::frame3::{
    Pt3, dir_sign_judge, indirect_cmp_coord_judge, indirect_orient3d_judge, orient3d_judge,
};
use nacre_math::Point3;
use nacre_predicates::{
    ThreePlane, det3_sign, indirect_cmp_coord, indirect_orient3d, orient2d, orient3d,
};
use nacre_scalar::{Orient, Rat};

/// A plane witnessed by three points known to lie on it, and — when the solid was rotated —
/// their exact [`Pt3`] definitions.
///
/// The rotation-general predicates need only this, which is why one implementation can serve
/// both index spaces (a plane class and a single face) without confusing them. Only
/// [`t_planes_coplanar`] uses the face form — it is the predicate that *defines* the classes,
/// so it necessarily runs before a plane table exists.
pub trait Witness {
    fn tri(&self) -> [Point3; 3];
    fn tri_pt3(&self) -> Option<&[Pt3; 3]>;
}

/// A witness that additionally carries its plane's exact coefficients — what the plane-class
/// predicates ([`t_orient3d`], [`t_cmp_coord`], [`t_plane_pair_dir_sign`]) need on the exact
/// path. (The stored↔outward `frame_sign` [`t_plane_pair_dir_sign`] also uses is *derived* from
/// `coeffs` + `tri` by [`frame_sign`], not required from the impl.) A single face (which never
/// plays a plane-class role) implements only [`Witness`].
pub trait PlaneWitness: Witness {
    /// The plane's exact (un-normalized) coefficients `[a, b, c, d]` (`n·x + d = 0`).
    fn coeffs(&self) -> [f64; 4];
}

/// Whether any of the named planes is rotated — the per-predicate routing signal. A plane
/// carries [`Witness::tri_pt3`] `Some` exactly when it came from a rotated solid. A predicate
/// must escalate to the kernel if **any** — not all — of its planes is irrational: a single
/// rounded coordinate can flip an f64 `orient3d`/`cmp`, whereas all-rational planes are exact.
pub fn any_rotated<W: Witness>(planes: &[W], idx: &[usize]) -> bool {
    idx.iter().any(|&k| planes[k].tri_pt3().is_some())
}

/// The three exact [`Pt3`] defining plane `k`. A rotated plane carries them cached in
/// [`Witness::tri_pt3`]; an axis-aligned plane (`None`) is built exactly from its `tri`
/// coordinates, which are already exact f64.
pub fn plane_def<W: Witness>(planes: &[W], k: usize) -> [Pt3; 3] {
    match planes[k].tri_pt3() {
        Some(t) => t.clone(),
        None => planes[k].tri().map(pt3_from_exact),
    }
}

/// An exact axis-aligned point as a tol-0 `Pt3` (its f64 coordinates are exact rationals).
fn pt3_from_exact(p: Point3) -> Pt3 {
    let a = p.as_array();
    let rat = |x: f64| Rat::try_from_f64(x).expect("axis-aligned coordinate is an exact rational");
    Pt3::at([rat(a[0]), rat(a[1]), rat(a[2])])
}

fn to_i8(o: Orient) -> i8 {
    match o {
        Orient::Positive => 1,
        Orient::Negative => -1,
        Orient::Zero => 0,
    }
}

/// Borrow an owned plane def as the `&Pt3` tuple the judges take.
fn borrow3(d: &[Pt3; 3]) -> (&Pt3, &Pt3, &Pt3) {
    (&d[0], &d[1], &d[2])
}

/// Borrow three owned plane defs as the tuples `indirect_cmp_coord_judge` takes.
fn borrow_triple(d: &[[Pt3; 3]; 3]) -> [(&Pt3, &Pt3, &Pt3); 3] {
    [borrow3(&d[0]), borrow3(&d[1]), borrow3(&d[2])]
}

/// The sign of `orient3d(V, tri_j)` where `V = ∩(planes p, q, r)` is an implicit point,
/// matching `order_along`'s shape (`+1`/`-1`/`0`).
///
/// - `!rotated`: the exact path — the implicit-point `orient3d` (Attene) on the stored plane
///   coefficients and `tri` coordinates.
/// - `rotated`: each of `p, q, r` and the explicit triangle `j` is taken as its three exact
///   [`Pt3`] ([`plane_def`]), and [`indirect_orient3d_judge`] decides the sign from the
///   definitions — never materializing `V` or reading the rounded `tri`.
///
/// Winding-invariant. A query plane `j` equal to one of `p, q, r` means the point lies on `j`,
/// so the sign is exactly `0` (a combinatorial identity) — decided before either numeric
/// branch, exact and cheap for both.
pub fn t_orient3d<W: PlaneWitness>(planes: &[W], p: usize, q: usize, r: usize, j: usize) -> i8 {
    if j == p || j == q || j == r {
        return 0;
    }
    if !any_rotated(planes, &[p, q, r, j]) {
        let tp = ThreePlane([planes[p].coeffs(), planes[q].coeffs(), planes[r].coeffs()]);
        let tj = planes[j].tri();
        return indirect_orient3d(&tp, tj[0].as_array(), tj[1].as_array(), tj[2].as_array());
    }
    let (dp, dq, dr, dj) = (
        plane_def(planes, p),
        plane_def(planes, q),
        plane_def(planes, r),
        plane_def(planes, j),
    );
    to_i8(indirect_orient3d_judge(
        borrow3(&dp),
        borrow3(&dq),
        borrow3(&dr),
        &dj[0],
        &dj[1],
        &dj[2],
    ))
}

/// The sign of `a[axis] − b[axis]` between the two implicit points `a = ∩(planes a…)` and
/// `b = ∩(planes b…)` (`+1` = `a[axis] > b[axis]`). `!rotated` → the exact implicit
/// `cmp_coord` on the stored coefficients; `rotated` → each triple's three planes as exact
/// `Pt3` → [`indirect_cmp_coord_judge`].
pub fn t_cmp_coord<W: PlaneWitness>(planes: &[W], a: [usize; 3], b: [usize; 3], axis: usize) -> i8 {
    if !any_rotated(planes, &[a[0], a[1], a[2], b[0], b[1], b[2]]) {
        let tp = |t: [usize; 3]| {
            ThreePlane([
                planes[t[0]].coeffs(),
                planes[t[1]].coeffs(),
                planes[t[2]].coeffs(),
            ])
        };
        return indirect_cmp_coord(&tp(a), &tp(b), axis);
    }
    let da = a.map(|k| plane_def(planes, k));
    let db = b.map(|k| plane_def(planes, k));
    to_i8(indirect_cmp_coord_judge(
        borrow_triple(&da),
        borrow_triple(&db),
        axis,
    ))
}

/// Whether the three points are **exactly collinear**, decided by the three coordinate-plane
/// projections of the cross product (each an exact `orient2d`). Non-collinearity is the
/// standing precondition of [`t_planes_coplanar`].
fn tri_collinear(t: [Point3; 3]) -> bool {
    let [a, b, c] = t.map(|p| p.as_array());
    let proj = |i: usize, j: usize| orient2d([a[i], a[j]], [b[i], b[j]], [c[i], c[j]]) == 0.0;
    proj(0, 1) && proj(1, 2) && proj(2, 0)
}

/// The exact side of triangle `tri`'s plane that the explicit point `p` lies on
/// (`+1`/`-1`/`0`), sharing [`t_orient3d`]'s convention (`p` takes `V`'s slot).
fn plane_side_exact(tri: [Point3; 3], p: Point3) -> i8 {
    let d = orient3d(
        p.as_array(),
        tri[0].as_array(),
        tri[1].as_array(),
        tri[2].as_array(),
    );
    match d.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => 1,
        Some(std::cmp::Ordering::Less) => -1,
        _ => 0,
    }
}

/// Whether planes `i` and `j` are the **same plane**, decided on the faces' original
/// coordinates instead of on their derived coefficients (three non-collinear points on a plane
/// determine it, so "every point of `tri_j` lies on `tri_i`'s plane" is conclusive — but only
/// under non-collinearity, so a degenerate `tri` answers `false`). `!rotated` → the exact
/// `orient3d` on `tri`; any rotated → the exact `Pt3` definitions and [`orient3d_judge`].
pub fn t_planes_coplanar<W: Witness>(planes: &[W], i: usize, j: usize) -> bool {
    if tri_collinear(planes[i].tri()) || tri_collinear(planes[j].tri()) {
        return false;
    }
    if !any_rotated(planes, &[i, j]) {
        return planes[j]
            .tri()
            .iter()
            .all(|&q| plane_side_exact(planes[i].tri(), q) == 0);
    }
    let (di, dj) = (plane_def(planes, i), plane_def(planes, j));
    dj.iter()
        .all(|q| to_i8(orient3d_judge(q, &di[0], &di[1], &di[2])) == 0)
}

/// `+1` if plane `w`'s stored coefficient-normal points the same way as its outward `tri`
/// normal, `-1` otherwise (`det(stored) = frame_sign · det(outward)`). Derived from the port's
/// `coeffs` + `tri` alone: `n_out` is *defined* as `cross(tri)` and the stored normal is a
/// positive multiple of `coeffs[0..3]`, so `sign(coeffs · cross(tri))` reproduces the
/// arrangement's stored `frame_sign` exactly (the two are co-sourced from one face).
fn frame_sign<W: PlaneWitness>(w: &W) -> i8 {
    let t = w.tri();
    let cross = (t[1] - t[0]).cross(t[2] - t[0]);
    let c = w.coeffs();
    let cx = cross.as_array();
    if c[0] * cx[0] + c[1] * cx[1] + c[2] * cx[2] > 0.0 {
        1
    } else {
        -1
    }
}

/// `sign(det[n_p; n_a; n_b])` over the three planes' stored normals — how the line `p ∩ a`
/// runs relative to plane `b`, matching `plane_pair_dir_sign`'s shape (`+1`/`-1`/`0`).
///
/// `!rotated` → `det3_sign` of the stored (un-normalized) normals. `rotated` → the kernel `D`
/// (det of the *outward* `tri` normals, [`dir_sign_judge`]) bridged to the *stored*-normal
/// convention by the per-plane [`frame_sign`]: `det(stored) =
/// frame_sign(p)·frame_sign(a)·frame_sign(b)·det(outward)`.
pub fn t_plane_pair_dir_sign<W: PlaneWitness>(planes: &[W], p: usize, a: usize, b: usize) -> i8 {
    if !any_rotated(planes, &[p, a, b]) {
        let row = |k: usize| {
            let [x, y, z, _] = planes[k].coeffs();
            [x, y, z]
        };
        return det3_sign([row(p), row(a), row(b)]);
    }
    let (dp, da, db) = (
        plane_def(planes, p),
        plane_def(planes, a),
        plane_def(planes, b),
    );
    frame_sign(&planes[p])
        * frame_sign(&planes[a])
        * frame_sign(&planes[b])
        * to_i8(dir_sign_judge(borrow3(&dp), borrow3(&da), borrow3(&db)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A synthetic axis-aligned plane witness: three points on the plane plus its exact
    /// coefficients. `tri_pt3` is `None` (axis-aligned), so predicates take the exact path.
    struct W {
        tri: [Point3; 3],
        coeffs: [f64; 4],
    }
    impl Witness for W {
        fn tri(&self) -> [Point3; 3] {
            self.tri
        }
        fn tri_pt3(&self) -> Option<&[Pt3; 3]> {
            None
        }
    }
    impl PlaneWitness for W {
        fn coeffs(&self) -> [f64; 4] {
            self.coeffs
        }
    }

    /// The three coordinate planes through `(1,1,1)`: `x=1`, `y=1`, `z=1`, plus a `z=0` plane.
    fn cube_corner_planes() -> Vec<W> {
        let p = |a, b, c| Point3::from_array([a, b, c]);
        vec![
            // x = 1  →  1·x + 0 + 0 − 1 = 0
            W {
                tri: [p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(1.0, 0.0, 1.0)],
                coeffs: [1.0, 0.0, 0.0, -1.0],
            },
            // y = 1
            W {
                tri: [p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 1.0)],
                coeffs: [0.0, 1.0, 0.0, -1.0],
            },
            // z = 1
            W {
                tri: [p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, 1.0)],
                coeffs: [0.0, 0.0, 1.0, -1.0],
            },
            // z = 0
            W {
                tri: [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)],
                coeffs: [0.0, 0.0, 1.0, 0.0],
            },
        ]
    }

    /// The corner `∩(x=1, y=1, z=1) = (1,1,1)` sits above the `z=0` plane, so its orient3d
    /// against `z=0` is definite (nonzero), and querying `z=1` (a defining plane) is exactly 0.
    #[test]
    fn t_orient3d_axis_definite_and_on_plane() {
        let ps = cube_corner_planes();
        // query plane j = 3 (z=0): definite.
        assert_ne!(t_orient3d(&ps, 0, 1, 2, 3), 0);
        // query plane j = 2 (z=1) is one of the defining planes → exactly 0.
        assert_eq!(t_orient3d(&ps, 0, 1, 2, 2), 0);
    }

    /// A plane is coplanar with itself; two distinct planes are not.
    #[test]
    fn t_planes_coplanar_reflexive_and_distinct() {
        let ps = cube_corner_planes();
        assert!(
            t_planes_coplanar(&ps, 0, 0),
            "a plane is coplanar with itself"
        );
        assert!(
            !t_planes_coplanar(&ps, 0, 1),
            "x=1 and y=1 are distinct planes"
        );
    }

    /// `any_rotated` is false for axis-aligned witnesses (`tri_pt3` is `None`).
    #[test]
    fn any_rotated_false_for_axis_aligned() {
        let ps = cube_corner_planes();
        assert!(!any_rotated(&ps, &[0, 1, 2, 3]));
    }

    /// `frame_sign` recomputes the stored-vs-outward sign from `coeffs` + `tri` alone: `+1`
    /// when the coefficient-normal agrees with `cross(tri)`, `-1` when the winding is reversed.
    #[test]
    fn frame_sign_from_coeffs_and_tri() {
        let ps = cube_corner_planes();
        // x=1: tri wound so cross(tri) = +x, and the coeffs normal is +x → +1.
        assert_eq!(frame_sign(&ps[0]), 1);
        // reversing the tri winding flips cross(tri) → -1 (coeffs unchanged).
        let flipped = W {
            tri: [ps[0].tri[0], ps[0].tri[2], ps[0].tri[1]],
            coeffs: ps[0].coeffs,
        };
        assert_eq!(frame_sign(&flipped), -1);
    }
}
