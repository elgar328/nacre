//! **What a motion chain folds into, and why a statement it carries can have no rational image.**

use super::*;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

type Q = Ratio<BigInt>;

/// **Why a statement carried along a [`ChainFold`] has no image in `Rat`.**
///
/// The two are asked in this order — whether the image is rational at all, then whether `Rat`
/// holds it — so a statement the chain puts off the rationals is never named for its width.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Carry {
    /// The chain carries this statement to an irrational image: no width would state it.
    Irrational,
    /// The image is rational and `Rat` cannot hold it, or it needs an offset whose value an
    /// overflow composing the chain lost.
    Width,
}

/// **A motion chain, folded** — the map a statement written against the chain's leaf is carried
/// to the world by, asked one statement at a time.
///
/// ```text
/// x ↦ A( F(x) )      A: the exact nodes' map (AxisAffine)
/// F(x) = t + Σ xᵢ·cᵢ/rᵢ   cᵢ rational, mutually orthogonal, rᵢ = |cᵢ| — none where no frame is folded in
/// ```
///
/// A frame's basis is `û = u/|u|, v̂ = v/|v|, ŵ = n/|n|` over rational `u, v, n`, so its map is
/// irrational as a whole whenever a length is — and yet a statement can land on the rationals:
/// a point on `ŵ` or a plane across it when `|n|` is rational (one through the frame's origin
/// whatever `|n|` is), a direction along one axis. So every carry answers for the statement it is
/// given: its exact image, or why there is none ([`Carry`]).
///
/// **Two parts, because the exact nodes compose in closed form.** Translations, quarter turns
/// and axis mirrors fold into one [`AxisAffine`], which never fails and keeps answering
/// directions when its offset overflows; a chain of them carries in `Rat` exactly as it always
/// has. The frame part exists only where a frame node is folded in, and a carry through it runs
/// at any width to the end — the frame part, then `A` — and is narrowed once: the geometry is asked
/// before any width can be.
///
/// A direction is asked two ways because two questions are asked of it: as a vector
/// ([`ChainFold::dir_rat`] — its length kept, for a caller that measures along it) and as a line
/// ([`ChainFold::dir_along`] — some positive multiple of its image, for a caller that only needs
/// where it points). Only the second can be rational where a length is not.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainFold {
    frame: Option<Box<FramePart>>,
    axis: AxisAffine,
}

/// `F(x) = t + Σ xᵢ·cᵢ/rᵢ` with `sᵢ = rᵢ² = cᵢ·cᵢ`, held at any width so a long frame axis does
/// not stop the fold; `det` is the sign the columns' order and directions give.
///
/// `t` is `None` where composing the fold needed an exact node's offset that overflowed `Rat`: the
/// origin is then unknown, and a point or a plane through it is not asked further (`Width`) —
/// directions do not read it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct FramePart {
    t: Option<[Q; 3]>,
    cols: [[Q; 3]; 3],
    s: [Q; 3],
    det: i8,
}

fn q(r: Rat) -> Q {
    BigRat::from_rat(r).0
}

fn narrow(x: &Q) -> Result<Rat, Carry> {
    BigRat(x.clone()).narrow().ok_or(Carry::Width)
}

fn narrow3(v: &[Q; 3]) -> Result<[Rat; 3], Carry> {
    Ok([narrow(&v[0])?, narrow(&v[1])?, narrow(&v[2])?])
}

/// `√x` when it is rational.
fn root(x: &Q) -> Result<Q, Carry> {
    let (n, d) = crate::sqrt::exact_root(x).ok_or(Carry::Irrational)?;
    Ok(Ratio::new_raw(n, d))
}

fn dot(a: &[Q; 3], b: &[Q; 3]) -> Q {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn cross(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

/// The integers of `v` over the lcm of its denominators — a positive multiple of `v`.
fn cleared<const N: usize>(v: &[Q; N]) -> [BigInt; N] {
    let l = v.iter().fold(BigInt::one(), |l, x| l.lcm(x.denom()));
    core::array::from_fn(|i| v[i].numer() * (&l / v[i].denom()))
}

/// The primitive integer direction along `v`, its sign kept — `v` itself when it is zero.
fn primitive(v: &[Q; 3]) -> [Q; 3] {
    let n = cleared(v);
    let g = n.iter().fold(BigInt::zero(), |g, x| g.gcd(x));
    if g.is_zero() {
        return v.clone();
    }
    n.map(|x| Q::from_integer(x / &g))
}

impl FramePart {
    /// The frame a [`PlaneFrame`] states: origin, `û`, `v̂ = ŵ × û`, `ŵ` — the basis a judge
    /// replays (`v` is `n × u`, as [`plane_frame_named`] derives it).
    fn of(pf: &PlaneFrame) -> FramePart {
        let (u, n) = (pf.u_raw.map(q), pf.n.map(q));
        let v = cross(&n, &u);
        FramePart {
            t: Some(pf.origin.map(q)),
            s: [dot(&u, &u), dot(&v, &v), dot(&n, &n)],
            cols: [u, v, n],
            det: 1,
        }
    }

    /// `Σ kᵢ·cᵢ`.
    fn comb(&self, k: &[Q; 3]) -> [Q; 3] {
        core::array::from_fn(|j| (0..3).map(|i| &k[i] * &self.cols[i][j]).sum())
    }

    /// The origin, where the fold kept it.
    fn origin(&self) -> Result<&[Q; 3], Carry> {
        self.t.as_ref().ok_or(Carry::Width)
    }

    fn point(&self, p: &[Q; 3]) -> Result<[Q; 3], Carry> {
        let v = self.vector(p)?;
        let t = self.origin()?;
        Ok(core::array::from_fn(|j| &t[j] + &v[j]))
    }

    /// `Σ dᵢ·cᵢ/rᵢ` — rational exactly when every axis `d` uses has a rational length.
    fn vector(&self, d: &[Q; 3]) -> Result<[Q; 3], Carry> {
        let mut k: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for i in 0..3 {
            if !d[i].is_zero() {
                k[i] = &d[i] / root(&self.s[i])?;
            }
        }
        Ok(self.comb(&k))
    }

    /// `rₖ·Σ dᵢ·cᵢ/rᵢ = Σ dᵢ·√(sₖ/sᵢ)·cᵢ` for `k` the first axis `d` uses — a positive multiple of
    /// the image, rational exactly when the axes `d` uses have lengths in rational ratio. No
    /// division by `d`, so a frame whose lengths are all `1` gives the image itself. The zero
    /// direction is its own image (`k` is then `0`, and reads nothing).
    fn along(&self, d: &[Q; 3]) -> Result<([Q; 3], usize), Carry> {
        let Some(k) = (0..3).find(|&i| !d[i].is_zero()) else {
            return Ok((d.clone(), 0));
        };
        let mut co: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for i in 0..3 {
            if !d[i].is_zero() {
                co[i] = &d[i] * root(&(&self.s[k] / &self.s[i]))?;
            }
        }
        Ok((self.comb(&co), k))
    }

    /// The image of the plane `α·x + δ = 0`, up to a positive factor: the normal along, and the
    /// constant `δ·rₖ − n'·t`.
    fn plane(&self, c: &[Q; 4]) -> Result<[Q; 4], Carry> {
        let a = [c[0].clone(), c[1].clone(), c[2].clone()];
        let (n, k) = self.along(&a)?;
        let lift = if c[3].is_zero() {
            Q::zero()
        } else {
            &c[3] * root(&self.s[k])?
        };
        let d = lift - dot(&n, self.origin()?);
        let [n0, n1, n2] = n;
        Ok([n0, n1, n2, d])
    }

    /// `xⱼ = cⱼ·(X − t)/rⱼ`.
    fn point_inverse(&self, x: &[Q; 3]) -> Result<[Q; 3], Carry> {
        let t = self.origin()?;
        let off: [Q; 3] = core::array::from_fn(|j| &x[j] - &t[j]);
        let mut out: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for (j, o) in out.iter_mut().enumerate() {
            let e = dot(&self.cols[j], &off);
            if !e.is_zero() {
                *o = e / root(&self.s[j])?;
            }
        }
        Ok(out)
    }

    /// A positive multiple of the direction whose image is `d`: `eⱼ·√(sₖ/sⱼ)` with `eⱼ = cⱼ·d`.
    /// The zero direction is its own preimage.
    fn along_inverse(&self, d: &[Q; 3]) -> Result<[Q; 3], Carry> {
        let e: [Q; 3] = core::array::from_fn(|j| dot(&self.cols[j], d));
        let Some(k) = (0..3).find(|&j| !e[j].is_zero()) else {
            return Ok(e);
        };
        let mut out: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for j in 0..3 {
            if !e[j].is_zero() {
                out[j] = &e[j] * root(&(&self.s[k] / &self.s[j]))?;
            }
        }
        Ok(out)
    }

    /// `G ∘ A ∘ self` as one frame part, `G` the outer one: each column goes out along `G` (a
    /// positive multiple, its square length read off the result) and the origin as a point.
    /// `Irrational` where `G` cannot carry a column along, or the origin, rationally — the only
    /// way it fails; an origin whose offset was lost stays unknown.
    fn under(&self, a: &AxisAffine, g: &FramePart) -> Result<FramePart, Carry> {
        let t = match (&self.t, a.offset_q()) {
            (Some(t), Some(o)) => {
                let moved = a.dir_q(t);
                known(g.point(&core::array::from_fn(|j| &moved[j] + &o[j])))?
            }
            _ => None,
        };
        let mut cols: [[Q; 3]; 3] = core::array::from_fn(|_| core::array::from_fn(|_| Q::zero()));
        let mut s: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for i in 0..3 {
            let (c, _) = g.along(&a.dir_q(&self.cols[i]))?;
            s[i] = dot(&c, &c);
            cols[i] = c;
        }
        Ok(FramePart {
            t,
            cols,
            s,
            det: self.det * a.det() * g.det,
        })
    }

    /// `self ∘ A` — the exact map `A` run first, its signed permutation folded into the columns.
    fn after(&self, a: &AxisAffine) -> Result<FramePart, Carry> {
        let mut cols: [[Q; 3]; 3] = core::array::from_fn(|_| core::array::from_fn(|_| Q::zero()));
        let mut s: [Q; 3] = core::array::from_fn(|_| Q::zero());
        for j in 0..3 {
            // `A·e_j` is `±e_π(j)`: column `π(j)` of this frame, with the sign.
            let e: [Q; 3] =
                core::array::from_fn(|i| Q::from_integer(BigInt::from(i32::from(i == j))));
            let (c, _) = self.along(&a.dir_q(&e))?;
            s[j] = dot(&c, &c);
            cols[j] = c;
        }
        let t = match a.offset_q() {
            Some(o) => known(self.point(&o))?,
            None => None,
        };
        Ok(FramePart {
            t,
            cols,
            s,
            det: self.det * a.det(),
        })
    }
}

/// A composed origin: `Irrational` stops the composition, `Width` (an origin already lost) keeps it
/// unknown.
fn known(p: Result<[Q; 3], Carry>) -> Result<Option<[Q; 3]>, Carry> {
    match p {
        Ok(p) => Ok(Some(p)),
        Err(Carry::Width) => Ok(None),
        Err(Carry::Irrational) => Err(Carry::Irrational),
    }
}

impl ChainFold {
    /// The fold of a chain of exact nodes — translations, quarter turns and axis mirrors.
    pub fn axis(axis: AxisAffine) -> ChainFold {
        ChainFold { frame: None, axis }
    }

    /// The map a frame node states: coordinates read in the frame `pf` names, carried into the
    /// frame its plane is written in.
    pub fn frame(pf: &PlaneFrame) -> ChainFold {
        ChainFold {
            frame: Some(Box::new(FramePart::of(pf))),
            axis: AxisAffine::translation([Rat::from_int(0); 3]),
        }
    }

    /// This fold with the exact node `outer` applied after it.
    pub fn then_axis(&self, outer: &AxisAffine) -> ChainFold {
        ChainFold {
            frame: self.frame.clone(),
            axis: self.axis.then(outer),
        }
    }

    /// **This fold, then `outer`, as one fold** — `None` where the composite has no form this
    /// type holds: `outer`'s frame cannot carry this fold's columns along it, or its origin, to
    /// the rationals. Width never stops it: the columns compose at any width, and an origin that
    /// needs an offset lost to overflow is kept unknown.
    pub fn then(&self, outer: &ChainFold) -> Option<ChainFold> {
        let Some(g) = &outer.frame else {
            return Some(self.then_axis(&outer.axis));
        };
        let inner = match &self.frame {
            Some(f) => f.under(&self.axis, g),
            None => g.after(&self.axis),
        };
        Some(ChainFold {
            frame: Some(Box::new(inner.ok()?)),
            axis: outer.axis,
        })
    }

    /// `+1` for a proper map, `−1` for one carrying an odd number of reflections.
    pub fn det(&self) -> i8 {
        self.axis.det() * self.frame.as_ref().map_or(1, |f| f.det)
    }

    /// `A` applied to a point at any width — `Width` where its offset was lost.
    fn axis_point(&self, p: &[Q; 3]) -> Result<[Q; 3], Carry> {
        let o = self.axis.offset_q().ok_or(Carry::Width)?;
        let m = self.axis.dir_q(p);
        Ok(core::array::from_fn(|j| &m[j] + &o[j]))
    }

    /// The image of the point `p`.
    pub fn point_rat(&self, p: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        match &self.frame {
            Some(f) => narrow3(&self.axis_point(&f.point(&p.map(q))?)?),
            None => self.axis.point_rat(p).ok_or(Carry::Width),
        }
    }

    /// The image of the vector `d`, its length kept.
    pub fn dir_rat(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        match &self.frame {
            Some(f) => narrow3(&self.axis.dir_q(&f.vector(&d.map(q))?)),
            None => self.axis.dir_rat(d).ok_or(Carry::Width),
        }
    }

    /// A positive multiple of the image of the direction `d` — where it points, not how long it
    /// is. Through a frame it is the primitive integer direction; without one, the image itself.
    pub fn dir_along(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        match &self.frame {
            Some(f) => narrow3(&primitive(&self.axis.dir_q(&f.along(&d.map(q))?.0))),
            None => self.axis.dir_rat(d).ok_or(Carry::Width),
        }
    }

    /// **The image of the cylinder `def`** — its origin carried as a point, its axis and seam
    /// reference as lines (where they point, not how long they are), its squared radius kept (the
    /// map is rigid).
    ///
    /// One statement, three carries: every one's geometry is asked before any width, so a
    /// cylinder whose seam reference the chain puts off the rationals is `Irrational` even where
    /// its origin needs an offset lost to overflow.
    pub fn cylinder(&self, def: &CylinderDef) -> Result<CylinderDef, Carry> {
        let carried = [
            self.point_rat(def.origin()),
            self.dir_along(def.dir()),
            self.dir_along(def.ref_dir()),
        ];
        if carried.contains(&Err(Carry::Irrational)) {
            return Err(Carry::Irrational);
        }
        let [origin, dir, ref_dir] = carried;
        CylinderDef::new(origin?, dir?, ref_dir?, def.r2().clone()).ok_or(Carry::Width)
    }

    /// **The image of the plane `c·(x, 1) = 0` as a name, and which way that name faces** — the
    /// canonical name of the image ([`PlaneName`], narrow or wide), with `+1` where it is the image
    /// carried with its sense and `−1` where canonicalizing turned it. One answer for both, so a
    /// reader that holds the name holds its facing.
    ///
    /// Through a frame, and for a wide name, the carry runs at any width, so the image's geometry
    /// is asked before its width; without a frame a narrow name is carried in `Rat`, where an
    /// overflow is `Width`.
    pub fn plane_name(&self, name: &PlaneName) -> Result<(PlaneName, i8), Carry> {
        if self.frame.is_none()
            && let Some(c) = name.narrow()
        {
            // The `Rat` road: `A` carries the normal to `M·n` with its sense, and canonicalizing
            // keeps that sense where the first nonzero components agree in sign. An overflow on it
            // is the answer — this road is not widened.
            let (Some(image), Some(n)) = (
                self.axis.plane_coeffs(*c),
                self.axis.dir_rat([c[0], c[1], c[2]]),
            ) else {
                return Err(Carry::Width);
            };
            let zero = Rat::from_int(0);
            let k = (0..3).find(|&k| image[k] != zero).unwrap_or(0);
            let kept = (n[k] > zero) == (image[k] > zero);
            return Ok((PlaneName::Narrow(image), if kept { 1 } else { -1 }));
        }
        self.image_name(&name.coeff_ints().map(Q::from_integer))
    }

    /// [`ChainFold::plane_name`] for any coefficients. The zero vector, which is no plane, is its
    /// own image.
    fn image_name(&self, c: &[Q; 4]) -> Result<(PlaneName, i8), Carry> {
        let image = match &self.frame {
            Some(f) => {
                let [n0, n1, n2, d] = f.plane(c)?;
                self.axis_plane(&[n0, n1, n2], d)?
            }
            None => {
                let n = [c[0].clone(), c[1].clone(), c[2].clone()];
                self.axis_plane(&n, c[3].clone())?
            }
        };
        Ok(canonical_name_ints(cleared(&image))
            .unwrap_or((PlaneName::Narrow([Rat::from_int(0); 4]), 1)))
    }

    /// `A` applied to the plane `n·x + d = 0`: `n' = M·n`, `d' = d − n'·t` — `Width` where the
    /// offset was lost.
    fn axis_plane(&self, n: &[Q; 3], d: Q) -> Result<[Q; 4], Carry> {
        let o = self.axis.offset_q().ok_or(Carry::Width)?;
        let [n0, n1, n2] = self.axis.dir_q(n);
        let d = d - (&n0 * &o[0] + &n1 * &o[1] + &n2 * &o[2]);
        Ok([n0, n1, n2, d])
    }

    /// The image of the plane `n·x + d = 0`, canonicalized ([`canonical_plane_coeffs`]) — the
    /// narrow name of [`ChainFold::plane_name`], or `Width` where it is wide.
    pub fn plane_coeffs(&self, c: [Rat; 4]) -> Result<[Rat; 4], Carry> {
        if self.frame.is_none() {
            return self.axis.plane_coeffs(c).ok_or(Carry::Width);
        }
        match self.image_name(&c.map(q))?.0 {
            PlaneName::Narrow(c) => Ok(c),
            PlaneName::Wide(_) => Err(Carry::Width),
        }
    }

    /// The point whose image is `p`.
    pub fn point_rat_inverse(&self, p: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        match &self.frame {
            Some(f) => {
                let o = self.axis.offset_q().ok_or(Carry::Width)?;
                let off: [Q; 3] = core::array::from_fn(|j| q(p[j]) - &o[j]);
                narrow3(&f.point_inverse(&self.axis.inverse().dir_q(&off))?)
            }
            None => self.axis.inverse().point_rat(p).ok_or(Carry::Width),
        }
    }

    /// A positive multiple of the direction whose image is `d` — through a frame, the primitive
    /// integer direction.
    pub fn dir_along_inverse(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        match &self.frame {
            Some(f) => narrow3(&primitive(
                &f.along_inverse(&self.axis.inverse().dir_q(&d.map(q)))?,
            )),
            None => self.axis.inverse().dir_rat(d).ok_or(Carry::Width),
        }
    }

    /// The image of an `f64` vector where the map only permutes and negates — `None` where a frame
    /// is folded in, whose irrational lengths an `f64` carry would round.
    pub fn dir_f64(&self, v: [f64; 3]) -> Option<[f64; 3]> {
        self.frame.is_none().then(|| self.axis.dir_f64(v))
    }
}

#[cfg(test)]
#[path = "tests/chain_fold.rs"]
mod tests;
