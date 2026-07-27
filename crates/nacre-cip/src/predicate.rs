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
use nacre_scalar::Orient;

/// A plane witnessed by three points known to lie on it, and — when the solid was rotated —
/// their exact [`Pt3`] definitions.
///
/// The rotation-general predicates need only this, which is why one implementation can serve
/// both index spaces (a plane class and a single face) without confusing them. Only
/// [`t_planes_coplanar`] uses the face form — it is the predicate that *defines* the classes,
/// so it necessarily runs before a plane table exists.
pub trait Witness {
    fn tri(&self) -> [Point3; 3];
    /// The three `tri` points as exact [`Pt3`] definitions, **always present**.
    ///
    /// This is a *cache*: the definition is built once, where the witness is, and every predicate
    /// borrows it. It used to be an `Option` whose emptiness *also* meant "not rotated" — one
    /// field answering two questions, which is why it could not be filled in advance without
    /// changing which predicate path runs. [`Witness::is_rotated`] is now that second question.
    fn tri_pt3(&self) -> &[Pt3; 3];
    /// Whether this plane came from a rotated solid — the predicate-routing signal, and a
    /// **different fact** from "does the definition carry a rotation chain". Neither
    /// `chain.is_empty()` nor `tol == 0` is equivalent to it (a rotated solid's face may witness
    /// its plane through points with no chain, and a 90°-family rotation has tol exactly 0), so
    /// the answer is carried, not derived.
    fn is_rotated(&self) -> bool;

    /// Identifies the rigid motion this witness's definition carries — `0` for none, and equal
    /// values **only** for structurally identical chains (same nodes, order and pivots).
    ///
    /// A rigid motion preserves every determinant these predicates take, so when all of a
    /// judgement's inputs carry one motion the answer is the answer on their pre-rotation
    /// coordinates — exactly, with no tolerance at all. This is what lets such a judgement leave
    /// the toleranced path entirely.
    fn chain_id(&self) -> u64;

    /// The witness triangle **before** that motion, or `None` when the pre-rotation coordinates
    /// are not `f64`-representable and so cannot be handed to the exact predicate.
    fn base_tri(&self) -> Option<[Point3; 3]>;

    /// The precision (bits) the escalating judges realize this operation's definitions at.
    ///
    /// **A property of the model, not a constant.** The error a realization carries grows with
    /// the rotation history — one bit per turn, measured — so a fixed precision decides the
    /// longest chain a model may have before its judgements stop separating. It is chosen once
    /// per boolean (see `nacre_ops::judge_precision`) and carried here because every predicate
    /// call site has the plane table in hand and nothing else. Uniform within an operation, which
    /// is what keeps [`Pt3`]'s realization cache warm.
    fn judge_prec(&self) -> usize;
}

/// A witness that additionally carries its plane's exact coefficients — what the plane-class
/// predicates ([`t_orient3d`], [`t_cmp_coord`], [`t_plane_pair_dir_sign`]) need on the exact
/// path. (The stored↔outward `frame_sign` [`t_plane_pair_dir_sign`] also uses is *derived* from
/// `coeffs` + `tri` by [`frame_sign`], not required from the impl.) A single face (which never
/// plays a plane-class role) implements only [`Witness`].
pub trait PlaneWitness: Witness {
    /// The plane's exact (un-normalized) coefficients `[a, b, c, d]` (`n·x + d = 0`).
    fn coeffs(&self) -> [f64; 4];

    /// The plane's coefficients before the rigid motion — see [`Witness::base_tri`].
    fn base_coeffs(&self) -> Option<[f64; 4]>;
}

/// Whether any of the named planes is rotated — the per-predicate routing signal. A predicate
/// must escalate to the kernel if **any** — not all — of its planes is irrational: a single
/// rounded coordinate can flip an f64 `orient3d`/`cmp`, whereas all-rational planes are exact.
pub fn any_rotated<W: Witness>(planes: &[W], idx: &[usize]) -> bool {
    idx.iter().any(|&k| planes[k].is_rotated())
}

/// The three exact [`Pt3`] defining plane `k` — **borrowed**, never rebuilt.
///
/// This used to construct them per call: cloning a rotated witness (a heap allocation each time)
/// or rebuilding an axis-aligned one from its `tri`. A boolean over 25 rotated fins called it a
/// million times, which was 77% of its runtime. The definitions are the same every call, so the
/// witness owns them and this is a pure accessor.
/// Can this judgement be answered exactly in the pre-rotation frame?
///
/// Yes when every input carries **one and the same** rigid motion (`chain_id`), that motion is not
/// the identity, and every pre-rotation witness is `f64`-representable. Then the rotation cancels
/// out of the determinant and the exact predicate answers on the bases. A mismatch is a
/// conservative miss — the toleranced path still answers, just more slowly.
fn shared_motion<W: Witness>(planes: &[W], idx: &[usize]) -> bool {
    let Some(&first) = idx.first() else {
        return false;
    };
    let id = planes[first].chain_id();
    id != 0 && idx.iter().all(|&k| planes[k].chain_id() == id)
}

pub fn plane_def<W: Witness>(planes: &[W], k: usize) -> &[Pt3; 3] {
    planes[k].tri_pt3()
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

/// Borrow three plane defs as the tuples `indirect_cmp_coord_judge` takes.
fn borrow_triple(d: [&[Pt3; 3]; 3]) -> [(&Pt3, &Pt3, &Pt3); 3] {
    [borrow3(d[0]), borrow3(d[1]), borrow3(d[2])]
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
    // One shared rigid motion ⇒ the same question, exactly, on the pre-rotation coordinates.
    if shared_motion(planes, &[p, q, r, j]) {
        if let (Some(cp), Some(cq), Some(cr), Some(tj)) = (
            planes[p].base_coeffs(),
            planes[q].base_coeffs(),
            planes[r].base_coeffs(),
            planes[j].base_tri(),
        ) {
            let tp = ThreePlane([cp, cq, cr]);
            return indirect_orient3d(&tp, tj[0].as_array(), tj[1].as_array(), tj[2].as_array());
        }
    }
    let (dp, dq, dr, dj) = (
        plane_def(planes, p),
        plane_def(planes, q),
        plane_def(planes, r),
        plane_def(planes, j),
    );
    to_i8(indirect_orient3d_judge(
        borrow3(dp),
        borrow3(dq),
        borrow3(dr),
        &dj[0],
        &dj[1],
        &dj[2],
        planes[p].judge_prec(),
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
    if let Some(s) = cancel_cmp_coord(planes, a, b, axis) {
        return s;
    }
    let da = a.map(|k| plane_def(planes, k));
    let db = b.map(|k| plane_def(planes, k));
    to_i8(indirect_cmp_coord_judge(
        borrow_triple(da),
        borrow_triple(db),
        axis,
        planes[a[0]].judge_prec(),
    ))
}

/// `cmp_coord` answered exactly in the pre-rotation frame, when it can be.
///
/// **Cancelling a rotation out of a one-coordinate comparison is not the same move as cancelling
/// it out of a determinant**, and the difference is why this used to be left undone: a rotation
/// mixes the axes, so `(Ra)[k] − (Rb)[k]` is not `(a − b)[k]` and the pre-rotation *order* along
/// an axis says nothing about the rotated one. That much is still true. What it misses is that
/// the comparison only ever sees the **difference** `d = a − b`, and a rotation about a pivot `c`
/// is `x ↦ R(x−c) + c`, so the pivots cancel and `Ra − Rb = R d` exactly — for a whole chain,
/// with a different pivot per node. Only the product of the rotations survives.
///
/// So when every plane carries one motion **about a single axis** `k` (total angle `Θ`), the
/// answer is decidable exactly:
///
/// - **asked axis `k`** — `(Rd)[k] = d[k]`, the rotation axis is fixed. The pre-rotation
///   comparison *is* the answer, sign included.
/// - **asked axis in the rotation plane** — with `(i, j)` the plane's axes,
///   `(Rd)[i] = d[i]·cos Θ − d[j]·sin Θ` and `(Rd)[j] = d[i]·sin Θ + d[j]·cos Θ`. Zero when
///   `d[i] = d[j] = 0` (the difference lies along the axis). Otherwise it would need
///   `tan Θ = ±d[i]/d[j]`, a **rational** tangent — which by Niven's theorem a rational-degree
///   angle has only at multiples of 45°. Outside that family the value is therefore **provably
///   nonzero**, so it must never be reported as a coincidence; it is left to escalate, where the
///   interval separates it. At a multiple of 45° `tan Θ ∈ {0, ±1}`, and the test is exact
///   rational arithmetic again.
///
/// Each `d[m] = 0` question is `indirect_cmp_coord(base_a, base_b, m) == 0` on the pre-rotation
/// coefficients — the exact predicate, unchanged.
///
/// Returns `None` when the shortcut does not apply (mixed axes, no shared motion, a base that is
/// not `f64`-representable, or a provably-nonzero in-plane case), and the toleranced path answers.
fn cancel_cmp_coord<W: PlaneWitness>(
    planes: &[W],
    a: [usize; 3],
    b: [usize; 3],
    axis: usize,
) -> Option<i8> {
    let all = [a[0], a[1], a[2], b[0], b[1], b[2]];
    if !shared_motion(planes, &all) {
        return None;
    }
    let (k, theta) = single_axis_motion(plane_def(planes, all[0]))?;
    let base = |t: [usize; 3]| -> Option<ThreePlane> {
        Some(ThreePlane([
            planes[t[0]].base_coeffs()?,
            planes[t[1]].base_coeffs()?,
            planes[t[2]].base_coeffs()?,
        ]))
    };
    let (ba, bb) = (base(a)?, base(b)?);
    // The pre-rotation comparison along `m`, whose sign is the sign of `d[m]`.
    let cmp = |m: usize| indirect_cmp_coord(&ba, &bb, m);

    let (i, j) = k.plane();
    if axis != i && axis != j {
        return Some(cmp(axis)); // the rotation axis is fixed: `(Rd)[k] = d[k]`
    }
    let (di, dj) = (cmp(i), cmp(j));
    if di == 0 && dj == 0 {
        return Some(0); // the difference lies along the rotation axis
    }
    let _ = (di, dj, theta);
    // Not both zero, so the value is nonzero **unless** `tan Θ` is rational — which by Niven a
    // rational-degree angle manages only on the 45° family. Either way the answer is left to the
    // toleranced path, but for opposite reasons, and neither can be settled here:
    //
    // - off the 45° family it is *provably nonzero*, so escalation is guaranteed to separate it —
    //   and it must never be reported as a coincidence;
    // - on the 45° family `(Rd)[i] = (√2/2)·(±d[i] ∓ d[j])`, whose sign needs the two differences
    //   **compared in magnitude**, not just their signs. `indirect_cmp_coord` returns a sign, so
    //   deciding it exactly would need a predicate for `sign(d[i] − d[j])` that does not exist
    //   yet. Escalation answers it correctly; only the *proof* is missing.
    None
}

/// The single rotation axis and total angle of a chain, or `None` if it turns about more than
/// one axis (then the product is not a rotation about a coordinate axis and this shortcut does
/// not apply — a conservative miss).
fn single_axis_motion(def: &[Pt3; 3]) -> Option<(nacre_scalar::Axis, nacre_scalar::Angle)> {
    let chain = &def[0].chain;
    let first = chain.first()?;
    let mut total = nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(0))?;
    for n in chain.iter() {
        if n.axis != first.axis {
            return None;
        }
        total = total.checked_add(n.angle.deg())?;
    }
    Some((first.axis, total))
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
    // "Same plane" is a statement about incidence, which a rigid motion preserves.
    if shared_motion(planes, &[i, j]) {
        if let (Some(ti), Some(tj)) = (planes[i].base_tri(), planes[j].base_tri()) {
            return tj.iter().all(|&q| plane_side_exact(ti, q) == 0);
        }
    }
    let (di, dj) = (plane_def(planes, i), plane_def(planes, j));
    dj.iter().all(|q| {
        to_i8(orient3d_judge(
            q,
            &di[0],
            &di[1],
            &di[2],
            planes[i].judge_prec(),
        )) == 0
    })
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
    // A determinant of normals: a rigid motion multiplies it by `det(R) = 1`, so one shared motion
    // means the pre-rotation normals give the same sign, exactly.
    if shared_motion(planes, &[p, a, b]) {
        let row = |k: usize| planes[k].base_coeffs().map(|[x, y, z, _]| [x, y, z]);
        if let (Some(rp), Some(ra), Some(rb)) = (row(p), row(a), row(b)) {
            return det3_sign([rp, ra, rb]);
        }
    }
    let (dp, da, db) = (
        plane_def(planes, p),
        plane_def(planes, a),
        plane_def(planes, b),
    );
    frame_sign(&planes[p])
        * frame_sign(&planes[a])
        * frame_sign(&planes[b])
        * to_i8(dir_sign_judge(
            borrow3(dp),
            borrow3(da),
            borrow3(db),
            planes[p].judge_prec(),
        ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The precision these fixtures judge at; production chooses it per model.
    const FIXTURE_PREC: usize = 256;

    /// A synthetic axis-aligned plane witness: three points on the plane, its exact
    /// coefficients, and the exact `Pt3` definition of those points. `is_rotated` is `false`,
    /// so predicates take the exact path — the definition is there but unused, which is
    /// precisely the arrangement the production tables now have.
    struct W {
        tri: [Point3; 3],
        coeffs: [f64; 4],
        def: [Pt3; 3],
    }
    impl W {
        fn new(tri: [Point3; 3], coeffs: [f64; 4]) -> W {
            let def = tri.map(|p| Pt3::exact(p.as_array()).expect("test coordinate"));
            W { tri, coeffs, def }
        }
    }
    impl Witness for W {
        fn tri(&self) -> [Point3; 3] {
            self.tri
        }
        fn tri_pt3(&self) -> &[Pt3; 3] {
            &self.def
        }
        fn is_rotated(&self) -> bool {
            false
        }
        // These witnesses carry no motion, so there is nothing to cancel.
        fn chain_id(&self) -> u64 {
            0
        }
        fn base_tri(&self) -> Option<[Point3; 3]> {
            None
        }
        fn judge_prec(&self) -> usize {
            FIXTURE_PREC
        }
    }
    impl PlaneWitness for W {
        fn coeffs(&self) -> [f64; 4] {
            self.coeffs
        }
        fn base_coeffs(&self) -> Option<[f64; 4]> {
            None
        }
    }

    /// The three coordinate planes through `(1,1,1)`: `x=1`, `y=1`, `z=1`, plus a `z=0` plane.
    fn cube_corner_planes() -> Vec<W> {
        let p = |a, b, c| Point3::from_array([a, b, c]);
        vec![
            // x = 1  →  1·x + 0 + 0 − 1 = 0
            W::new(
                [p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(1.0, 0.0, 1.0)],
                [1.0, 0.0, 0.0, -1.0],
            ),
            // y = 1
            W::new(
                [p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 1.0)],
                [0.0, 1.0, 0.0, -1.0],
            ),
            // z = 1
            W::new(
                [p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(0.0, 1.0, 1.0)],
                [0.0, 0.0, 1.0, -1.0],
            ),
            // z = 0
            W::new(
                [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)],
                [0.0, 0.0, 1.0, 0.0],
            ),
        ]
    }

    /// A plane witness carrying one shared rotation — what [`cancel_cmp_coord`] needs.
    struct RW {
        tri: [Point3; 3],
        coeffs: [f64; 4],
        def: [Pt3; 3],
        base: [Point3; 3],
        base_coeffs: [f64; 4],
    }

    /// `[a,b,c,d]` of the plane through three points (`n = e1 × e2`, `d = −n·p0`).
    fn plane_of(t: [[f64; 3]; 3]) -> [f64; 4] {
        let e = |i: usize| [t[i][0] - t[0][0], t[i][1] - t[0][1], t[i][2] - t[0][2]];
        let (u, v) = (e(1), e(2));
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        [
            n[0],
            n[1],
            n[2],
            -(n[0] * t[0][0] + n[1] * t[0][1] + n[2] * t[0][2]),
        ]
    }

    impl RW {
        /// Three integer points on a plane, turned `deg`° about `axis` through a non-origin
        /// pivot — so the pivot translation is present and has to cancel in the difference.
        fn rotated(pts: [[i128; 3]; 3], axis: Axis, deg: i128) -> RW {
            RW::turned(pts, &[(axis, deg)])
        }

        /// The same, through a chain of turns — one node per `(axis, deg)`.
        fn turned(pts: [[i128; 3]; 3], nodes: &[(Axis, i128)]) -> RW {
            let pivot = [ri(1, 3), ri(1, 7), ri(0, 1)];
            let def: [Pt3; 3] = pts.map(|q| {
                let mut p = Pt3::at([ri(q[0], 1), ri(q[1], 1), ri(q[2], 1)]);
                for &(axis, deg) in nodes {
                    p = p.rotate_about(axis, Angle::from_deg(ri(deg, 1)).unwrap(), pivot);
                }
                p
            });
            let base = pts.map(|q| Point3::from_array([q[0] as f64, q[1] as f64, q[2] as f64]));
            let tri: [Point3; 3] = std::array::from_fn(|i| Point3::from_array(def[i].coord));
            RW {
                coeffs: plane_of(tri.map(|p| p.as_array())),
                base_coeffs: plane_of(base.map(|p| p.as_array())),
                tri,
                def,
                base,
            }
        }
    }
    impl Witness for RW {
        fn tri(&self) -> [Point3; 3] {
            self.tri
        }
        fn tri_pt3(&self) -> &[Pt3; 3] {
            &self.def
        }
        fn is_rotated(&self) -> bool {
            true
        }
        fn chain_id(&self) -> u64 {
            1 // one shared motion for every witness in these fixtures
        }
        fn base_tri(&self) -> Option<[Point3; 3]> {
            Some(self.base)
        }
        fn judge_prec(&self) -> usize {
            FIXTURE_PREC
        }
    }
    impl PlaneWitness for RW {
        fn coeffs(&self) -> [f64; 4] {
            self.coeffs
        }
        fn base_coeffs(&self) -> Option<[f64; 4]> {
            Some(self.base_coeffs)
        }
    }

    fn ri(n: i128, d: i128) -> nacre_scalar::Rat {
        nacre_scalar::Rat::new(n, d).unwrap()
    }
    use nacre_scalar::{Angle, Axis};

    /// **The pre-rotation shortcut must give the answer the escalation gives.**
    ///
    /// [`cancel_cmp_coord`] answers in the pre-rotation frame with exact rational arithmetic;
    /// [`indirect_cmp_coord_judge`] answers by realizing the rotation in astro-float and reading
    /// an interval. They share no bound and no code, so agreement is a real cross-check — and it
    /// is the only one available, because the shortcut's whole point is to *not* do what the
    /// judge does.
    ///
    /// The fixtures turn about a **non-origin pivot**, which is the step the original comment
    /// missed: a rotation about `c` is `x ↦ R(x−c)+c`, and the `+c` cancels in a difference. If
    /// it did not, every answer here would be wrong.
    #[test]
    fn the_pre_rotation_shortcut_agrees_with_the_escalation() {
        for deg in [30i128, 37, 120, 200] {
            for axis in [Axis::X, Axis::Y, Axis::Z] {
                // ∩(x=1, y=2, z=3) = (1,2,3) and ∩(x=1, y=2, z=7) = (1,2,7): the difference is
                // along z, so an in-plane comparison is exactly 0 and a z comparison is definite.
                // Then a pair differing along x, which the shortcut must decline rather than
                // guess (`tan Θ` is irrational, so it is provably nonzero but not signed here).
                let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg);
                let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg);
                let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg);
                let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg);
                let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg);
                let ps = vec![px, py, z3, z7, x5];
                for (a, b, what) in [
                    ([0usize, 1, 2], [0usize, 1, 3], "differ along z"),
                    ([0, 1, 2], [4, 1, 2], "differ along x"),
                ] {
                    for k in 0..3 {
                        let got = t_cmp_coord(&ps, a, b, k);
                        let want = to_i8(indirect_cmp_coord_judge(
                            borrow_triple(a.map(|i| plane_def(&ps, i))),
                            borrow_triple(b.map(|i| plane_def(&ps, i))),
                            k,
                            FIXTURE_PREC,
                        ));
                        assert_eq!(
                            got, want,
                            "{deg}° about {axis:?}, {what}, axis {k}: shortcut {got} vs \
                             escalation {want}"
                        );
                    }
                }
            }
        }
    }

    /// …and the shortcut must actually fire, or the test above only proves the escalation agrees
    /// with itself. The rotation axis is answered exactly, and so is a difference lying along it.
    #[test]
    fn the_shortcut_fires_where_it_should_and_declines_where_it_cannot() {
        let (axis, deg) = (Axis::Z, 30i128);
        let px = RW::rotated([[1, 0, 0], [1, 1, 0], [1, 0, 1]], axis, deg);
        let py = RW::rotated([[0, 2, 0], [1, 2, 0], [0, 2, 1]], axis, deg);
        let z3 = RW::rotated([[0, 0, 3], [1, 0, 3], [0, 1, 3]], axis, deg);
        let z7 = RW::rotated([[0, 0, 7], [1, 0, 7], [0, 1, 7]], axis, deg);
        let x5 = RW::rotated([[5, 0, 0], [5, 1, 0], [5, 0, 1]], axis, deg);
        let ps = vec![px, py, z3, z7, x5];
        let (a, b) = ([0usize, 1, 2], [0usize, 1, 3]); // differ along z only
        // Z is the rotation axis: preserved, so the sign comes back exactly.
        assert_eq!(cancel_cmp_coord(&ps, a, b, 2), Some(-1), "z: 3 < 7");
        // x and y: the difference lies along the rotation axis, so both are exactly 0.
        assert_eq!(cancel_cmp_coord(&ps, a, b, 0), Some(0));
        assert_eq!(cancel_cmp_coord(&ps, a, b, 1), Some(0));
        // A difference in the rotation plane: provably nonzero, but its sign needs the rotation
        // realized, so the shortcut declines instead of guessing.
        let c = [4usize, 1, 2];
        assert_eq!(cancel_cmp_coord(&ps, a, c, 0), None);
        assert_eq!(cancel_cmp_coord(&ps, a, c, 1), None);
        // …while the rotation axis still answers for that pair (both points share z = 3).
        assert_eq!(cancel_cmp_coord(&ps, a, c, 2), Some(0));
    }

    /// **A chain that turns about more than one axis must be declined, not answered.**
    ///
    /// The whole derivation rests on the product of the rotations being a rotation *about a
    /// coordinate axis* — that is what makes one coordinate fixed and the other two a plane
    /// rotation with a single angle. Compose an X turn with a Z turn and none of that holds:
    /// there is no preserved coordinate, and the in-plane formula is about the wrong plane. The
    /// guard is the only thing standing between that and a confidently wrong sign, and without a
    /// mixed-axis fixture nothing else in the suite notices if it is removed.
    #[test]
    fn a_chain_about_two_axes_is_declined() {
        let nodes: &[(Axis, i128)] = &[(Axis::X, 30), (Axis::Z, 40)];
        let t = |pts| RW::turned(pts, nodes);
        let ps = vec![
            t([[1, 0, 0], [1, 1, 0], [1, 0, 1]]),
            t([[0, 2, 0], [1, 2, 0], [0, 2, 1]]),
            t([[0, 0, 3], [1, 0, 3], [0, 1, 3]]),
            t([[0, 0, 7], [1, 0, 7], [0, 1, 7]]),
        ];
        let (a, b) = ([0usize, 1, 2], [0usize, 1, 3]);
        for k in 0..3 {
            assert_eq!(
                cancel_cmp_coord(&ps, a, b, k),
                None,
                "axis {k}: a two-axis chain has no preserved coordinate, so nothing here is \
                 decidable in the pre-rotation frame"
            );
            // …and the toleranced path still answers it, so declining costs only speed.
            let want = to_i8(indirect_cmp_coord_judge(
                borrow_triple(a.map(|i| plane_def(&ps, i))),
                borrow_triple(b.map(|i| plane_def(&ps, i))),
                k,
                FIXTURE_PREC,
            ));
            assert_eq!(t_cmp_coord(&ps, a, b, k), want);
        }
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
        let flipped = W::new([ps[0].tri[0], ps[0].tri[2], ps[0].tri[1]], ps[0].coeffs);
        assert_eq!(frame_sign(&flipped), -1);
    }
}
