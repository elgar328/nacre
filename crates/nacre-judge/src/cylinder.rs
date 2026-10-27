//! **Where a point on a cylinder is, and which side of a plane it lies on** — the pierce point
//! `p ∩ q ∩ cylinder` and the questions asked of it, beside the plane predicates of
//! [`crate::predicate`].
//!
//! A cylinder that reaches a judgement is always stated exactly (`nacre-ops`' gate refuses one
//! with no world statement before the arrangement starts), so what can be irrational is a plane:
//! a class with no world statement is known by its witness triangle alone. A pierce point of two
//! stated planes is the exact `a + b√c` pair; one of a plane with no statement is an interval box,
//! solved from the witness's realization ([`nacre_exact::quad::plane_plane_cylinder_hp`]).
//!
//! The plane rows are read through the same ports as every other predicate
//! ([`PlaneWitness`]): the world name's narrow coefficients, its sense, the witness triangle and
//! `frame_sign`. Nothing here reads a cached coordinate.

use crate::kernel::frame3::{meet_hp, plane_hp};
use crate::predicate::{Judge, PlaneWitness};
use nacre_exact::quad::{MeetLine, QuadRoot, QuadVal};
use nacre_exact::{CylinderDef, HpBounded, Orient, Rat};

/// **Where a pierce point is** — what [`Judge::pierce`] solves a pierce name to, and the one vessel
/// every question about the point reads.
///
/// The questions are this type's methods, so a caller asks *what it wants to know* (a side, an
/// interval) and never which number answers it. A caller that needs the exact pair itself — an
/// `a + b√c` comparison, a rational root — takes [`PiercePoint::into_exact`] or
/// [`PiercePoint::exact`] and declines where there is none.
// The exact arm is the wide one and the common one; its holders already box the whole point, so
// boxing the arm again would buy nothing.
#[allow(clippy::large_enum_variant)]
#[derive(Clone)]
pub enum PiercePoint {
    /// The meet line of the name's two planes and the root along it: `base + s·dir` with
    /// `s = a + b√c`, both planes stated in the world.
    Exact(MeetLine, QuadVal),
    /// A point one of whose planes has no world statement: its definition — the name's two planes,
    /// its root and its cylinder — and the interval box that holds it at the operation's precision.
    /// The definition is what the questions below answer **by name** where an interval cannot (a
    /// point on one of its own planes; two points of one ruling), and what realizes the point
    /// deeper.
    Bounded {
        at: Box<[HpBounded; 3]>,
        planes: [usize; 2],
        root: QuadRoot,
        cyl: usize,
        def: Box<CylinderDef>,
    },
}

impl<W: PlaneWitness> Judge<'_, W> {
    /// **Solve the pierce name `(planes, cyl, root)`** on cylinder `def`.
    ///
    /// Both planes stated in the world: the exact meet ([`nacre_exact::quad::plane_plane_cylinder`])
    /// and the root by name ([`QuadRoot::crossing`]). Otherwise the interval meet at this
    /// operation's precision, which holds a crossing only where it proves two roots apart.
    ///
    /// `None` for:
    /// * **a stated pair the exact solve cannot carry** — a plane with no narrow rational
    ///   description (a `Wide` name), or the checked arithmetic leaving `Rat`;
    /// * **a meet without the name's root** — the name was minted from this same meet, so only a
    ///   naming defect produces it;
    /// * on the interval road, **a root the box cannot form**: a tangency (an interval proves two
    ///   roots apart, never one double root) or a discriminant it cannot separate from zero.
    pub fn pierce(
        &self,
        planes: [usize; 2],
        cyl: usize,
        def: &CylinderDef,
        root: QuadRoot,
    ) -> Option<PiercePoint> {
        let (o, m, r2) = (def.origin(), def.dir(), def.r2());
        // **The interval road is for a plane with no world statement, and only for it.** A named
        // pair whose exact solve fails (a wide name, an overflow) keeps declining here under its
        // own name: answering it by intervals would put a second number on a population the exact
        // road owns.
        if planes
            .iter()
            .any(|&c| self.planes[c].world_name().is_none())
        {
            let at = self.pierce_box(planes, def, root, self.standard.prec)?;
            return Some(PiercePoint::Bounded {
                at: Box::new(at),
                planes,
                root,
                cyl,
                def: Box::new(def.clone()),
            });
        }
        let (p1, p2) = (self.world_rat(planes[0])?, self.world_rat(planes[1])?);
        let (line, s) = root.crossing(&nacre_exact::quad::plane_plane_cylinder(
            &p1, &p2, &o, &m, r2,
        )?)?;
        #[cfg(debug_assertions)]
        self.debug_assert_bounded_holds_exact(planes, def, root, &line, &s);
        Some(PiercePoint::Exact(line, s))
    }

    /// **The pierce name `(planes, root)` on cylinder `def` as an interval box at `prec`** — the
    /// interval meet ([`nacre_exact::quad::plane_plane_cylinder_hp`]) of the two planes as the name
    /// orders them ([`Self::canon_hp`]). `None` for a tangency (an interval proves two roots apart,
    /// never one double root) or a discriminant it cannot separate from zero.
    fn pierce_box(
        &self,
        planes: [usize; 2],
        def: &CylinderDef,
        root: QuadRoot,
        prec: usize,
    ) -> Option<[HpBounded; 3]> {
        let line = nacre_exact::quad::plane_plane_cylinder_hp(
            &self.canon_hp(planes[0], prec),
            &self.canon_hp(planes[1], prec),
            &def.origin(),
            &def.dir(),
            def.r2(),
            prec,
        )?;
        match root {
            QuadRoot::Lo => line.root(false, prec),
            QuadRoot::Hi => line.root(true, prec),
            QuadRoot::Double => None,
        }
    }

    /// **The interval road answers for the named population too, where it can be checked.** Every
    /// exact pierce point is also solved the interval way
    /// ([`nacre_exact::quad::plane_plane_cylinder_hp`] on the same coefficients, the same root by
    /// name), and the two boxes must overlap — both hold the true point. A convention that
    /// differed between the roads (a root order, a line direction) would put the interval point on
    /// the other ruling and fail here on the first named crossing.
    #[cfg(debug_assertions)]
    fn debug_assert_bounded_holds_exact(
        &self,
        planes: [usize; 2],
        def: &CylinderDef,
        root: QuadRoot,
        line: &MeetLine,
        s: &QuadVal,
    ) {
        let prec = self.standard.prec;
        let Some(at) = self.pierce_box(planes, def, root, prec) else {
            return;
        };
        let Some(ex) = PiercePoint::Exact(line.clone(), *s).point_hp(self) else {
            return;
        };
        debug_assert!(
            (0..3).all(|k| ex[k].sub(&at[k], prec).sign().is_none()),
            "the interval road put a named pierce point somewhere the exact road does not"
        );
    }

    /// Plane `k`'s world name as narrow rational coefficients — the one description the cylinder
    /// roads compare against a world axis. `None` with no world statement, or a `Wide` one.
    pub fn world_rat(&self, k: usize) -> Option<[Rat; 4]> {
        self.planes[k].world_name()?.narrow().copied()
    }

    /// The same coefficients **turned to face the way the plane does** — the name's normal times
    /// its sense ([`PlaneWitness::world_sense`]).
    fn stored_world_rat(&self, k: usize) -> Option<[Rat; 4]> {
        let c = self.world_rat(k)?;
        if self.planes[k].world_sense()? > 0 {
            return Some(c);
        }
        let zero = Rat::from_int(0);
        Some([
            zero.checked_sub(c[0])?,
            zero.checked_sub(c[1])?,
            zero.checked_sub(c[2])?,
            zero.checked_sub(c[3])?,
        ])
    }

    /// **Canonical → outward**: the sign that carries plane `k`'s world name ([`Self::world_rat`])
    /// to the frame `orient3d` answers in — `κ · frame_sign`.
    ///
    /// **The name is the plane's *name*, not an oriented normal.** `orient3d` answers against the
    /// plane's **outward** normal, and the name may point either way along the plane's facing. The
    /// world name's sense (`κ`, [`PlaneWitness::world_sense`]) relates the two — read off the
    /// truth, never off a plane cache — and `frame_sign` carries the plane's facing to the root
    /// face's outward. `None` with no world statement.
    pub fn outward_fix(&self, k: usize) -> Option<i8> {
        let w = &self.planes[k];
        Some(w.world_sense()? * w.frame_sign())
    }

    /// **Plane `k` facing the way it is stored, as an interval** — the world statement turned by its
    /// sense exactly where there is one, else the witness triangle's plane turned by `frame_sign`:
    /// the triangle runs outward (`n_out = frame_sign · n_P`), so `frame_sign` carries it back to
    /// the stored normal.
    ///
    /// For a plane with no world statement this **is** its canonical direction too: such a plane
    /// has no name whose sign could disagree, so the one direction it has — the stored one — is the
    /// one a pierce name's roots are ordered along (`n₁ × n₂`).
    pub fn stored_hp(&self, k: usize, prec: usize) -> [HpBounded; 4] {
        if let Some(w) = self.stored_world_rat(k) {
            return w.map(|x| HpBounded::of_rat(x, prec));
        }
        self.witness_hp(k, prec).map(|x| {
            if self.planes[k].frame_sign() < 0 {
                negate(&x, prec)
            } else {
                x
            }
        })
    }

    /// Plane `k` facing **outward** (`n_out = frame_sign · n_P`, the frame `orient3d` answers in)
    /// as an interval — exactly from the world statement where there is one, else the witness
    /// triangle, which runs outward by construction.
    pub fn outward_hp(&self, k: usize, prec: usize) -> [HpBounded; 4] {
        if self.planes[k].world_name().is_some() {
            return self.stored_hp(k, prec).map(|x| {
                if self.planes[k].frame_sign() < 0 {
                    negate(&x, prec)
                } else {
                    x
                }
            });
        }
        self.witness_hp(k, prec)
    }

    /// The coefficients a pierce name's roots are ordered by — the world name exactly where there
    /// is one (what the exact solve reads), else the stored facing ([`Self::stored_hp`]).
    pub fn canon_hp(&self, k: usize, prec: usize) -> [HpBounded; 4] {
        match self.world_rat(k) {
            Some(w) => w.map(|x| HpBounded::of_rat(x, prec)),
            None => self.stored_hp(k, prec),
        }
    }

    /// The witness triangle's plane ([`plane_hp`]) — outward.
    fn witness_hp(&self, k: usize, prec: usize) -> [HpBounded; 4] {
        let t = self.planes[k].tri_pt3();
        plane_hp(&t[0], &t[1], &t[2], prec)
    }
}

/// **A point a cylinder question orders** — the meet of three planes, or a pierce point.
#[derive(Clone, Copy)]
pub enum CylPoint<'p> {
    /// The three planes' meet, by their table indices.
    Three([usize; 3]),
    /// A pierce point.
    Pierce(&'p PiercePoint),
}

/// **A coordinate along a cylinder's axis** (`m · x`) — stated exactly (`a + b√c`, for a point
/// with an exact meet or on a plane across the axis), or the point whose interval box gives it.
#[derive(Clone, Copy)]
pub enum AxisValue<'p> {
    /// The coordinate, exact.
    Stated(&'p QuadVal),
    /// A point with no exact meet.
    Point(&'p PiercePoint),
}

impl<W: PlaneWitness> Judge<'_, W> {
    /// **A three-plane point as an interval** — its three planes' witness triangles ([`plane_hp`])
    /// met at `prec` ([`meet_hp`]): a box that holds the true point whether or not the planes have
    /// a world statement. `None` for a meet the precision could not separate.
    pub fn three_hp(&self, t: [usize; 3], prec: usize) -> Option<[HpBounded; 3]> {
        meet_hp(&t.map(|k| self.witness_hp(k, prec)), prec)
    }

    /// Either kind of point as an interval at `prec`.
    fn point_box(&self, p: CylPoint<'_>, prec: usize) -> Option<[HpBounded; 3]> {
        match p {
            CylPoint::Three(t) => self.three_hp(t, prec),
            CylPoint::Pierce(m) => m.point_hp_at(self, prec),
        }
    }

    /// **Which of two points has the larger coordinate on world axis `axis`** — `+1` for `a`, read
    /// off their interval boxes. `None` where a box cannot be formed or the two cannot be told
    /// apart: these are two distinct names, so an unseparated pair is never read as a tie.
    pub fn world_order(&self, a: CylPoint<'_>, b: CylPoint<'_>, axis: usize) -> Option<i8> {
        sign_i8(&self.world_diff(a, b, axis, self.standard.prec)?)
    }

    /// `x_a − x_b` on world axis `axis`, at `prec`.
    fn world_diff(
        &self,
        a: CylPoint<'_>,
        b: CylPoint<'_>,
        axis: usize,
        prec: usize,
    ) -> Option<HpBounded> {
        let (x, y) = (self.point_box(a, prec)?, self.point_box(b, prec)?);
        Some(x[axis].sub(&y[axis], prec))
    }

    /// **`sign((a − b) · d)` for two points on the line `p ∩ q`**, `d` its stored direction
    /// (`n_p × n_q`, [`Self::stored_hp`]) — read off interval boxes. `None` where they cannot
    /// separate the two: distinct names, so never a tie.
    pub fn line_order(&self, p: usize, q: usize, a: &PiercePoint, b: &PiercePoint) -> Option<i8> {
        sign_i8(&self.line_diff(p, q, a, b, self.standard.prec)?)
    }

    /// `(a − b) · (n_p × n_q)` at `prec`.
    fn line_diff(
        &self,
        p: usize,
        q: usize,
        a: &PiercePoint,
        b: &PiercePoint,
        prec: usize,
    ) -> Option<HpBounded> {
        let (x, y) = (a.point_hp_at(self, prec)?, b.point_hp_at(self, prec)?);
        let (np, nq) = (self.stored_hp(p, prec), self.stored_hp(q, prec));
        let d = [
            np[1].mul(&nq[2], prec).sub(&np[2].mul(&nq[1], prec), prec),
            np[2].mul(&nq[0], prec).sub(&np[0].mul(&nq[2], prec), prec),
            np[0].mul(&nq[1], prec).sub(&np[1].mul(&nq[0], prec), prec),
        ];
        (0..3)
            .map(|k| x[k].sub(&y[k], prec).mul(&d[k], prec))
            .reduce(|s, t| s.add(&t, prec))
    }

    /// **The order of two coordinates along cylinder `def`'s axis.** Two stated ones compare
    /// exactly ([`QuadVal::cmp_value`]); a pair with a point compares by intervals — the point's
    /// box dotted with the axis, a stated value realized — and is never a tie: `None` where they
    /// cannot be separated.
    pub fn axis_order(
        &self,
        def: &CylinderDef,
        a: AxisValue<'_>,
        b: AxisValue<'_>,
    ) -> Option<core::cmp::Ordering> {
        use core::cmp::Ordering;
        if let (AxisValue::Stated(x), AxisValue::Stated(y)) = (a, b) {
            return x.cmp_value(y);
        }
        match self.axis_diff(def, a, b, self.standard.prec)?.sign()? {
            true => Some(Ordering::Greater),
            false => Some(Ordering::Less),
        }
    }

    /// `m · (x_a − x_b)` at `prec`, a stated value realized.
    fn axis_diff(
        &self,
        def: &CylinderDef,
        a: AxisValue<'_>,
        b: AxisValue<'_>,
        prec: usize,
    ) -> Option<HpBounded> {
        let m = def.dir().map(|x| HpBounded::of_rat(x, prec));
        let hp = |v: AxisValue<'_>| match v {
            AxisValue::Stated(q) => nacre_exact::quad_hp(q, prec),
            AxisValue::Point(p) => {
                let at = p.point_hp_at(self, prec)?;
                Some(
                    at[0]
                        .mul(&m[0], prec)
                        .add(&at[1].mul(&m[1], prec), prec)
                        .add(&at[2].mul(&m[2], prec), prec),
                )
            }
        };
        Some(hp(a)?.sub(&hp(b)?, prec))
    }
}

impl<W: PlaneWitness> Judge<'_, W> {
    /// The line `p ∩ q` against cylinder `def`, as intervals at this operation's precision — the
    /// two planes as a pierce name orders them ([`Self::canon_hp`]).
    fn meet_line_hp(
        &self,
        p: usize,
        q: usize,
        def: &CylinderDef,
        prec: usize,
    ) -> Option<nacre_exact::quad::LineQuadHp> {
        nacre_exact::quad::plane_plane_cylinder_hp(
            &self.canon_hp(p, prec),
            &self.canon_hp(q, prec),
            &def.origin(),
            &def.dir(),
            def.r2(),
            prec,
        )
    }

    /// **Does the line `p ∩ q` cross cylinder `def` in two points?** — `Some(true)` two crossings,
    /// `Some(false)` a miss, both read off the discriminant's sign once the line is proved not
    /// parallel to the axis (`A > 0`). `None` where the intervals cannot say: a tangency or an
    /// axis-parallel line is a statement an interval cannot make.
    pub fn line_crosses(&self, p: usize, q: usize, def: &CylinderDef) -> Option<bool> {
        let line = self.meet_line_hp(p, q, def, self.standard.prec)?;
        if line.a.sign() != Some(true) {
            return None;
        }
        line.disc(self.standard.prec).sign()
    }

    /// **Is the line `p ∩ q`, exactly parallel to the axis, off cylinder `def`'s surface?** — its
    /// constant term `C ≠ 0` proved (inside the solid or outside it, never on the surface). `false`
    /// where it is not proved; the caller states the parallelism.
    pub fn parallel_line_off(&self, p: usize, q: usize, def: &CylinderDef) -> bool {
        self.meet_line_hp(p, q, def, self.standard.prec)
            .is_some_and(|line| line.c.sign().is_some())
    }

    /// **Does plane `k`, along cylinder `def`'s axis, cut its lateral?** — the axis's offset from
    /// the plane `n·o + d`, squared, against `r²|n|²`, on the plane's outward interval
    /// ([`Self::outward_hp`]). `Some(true)` within the radius (two rulings), `Some(false)` beyond it
    /// (the whole cylinder on one side), `None` where the intervals cannot say.
    pub fn plane_cuts(&self, k: usize, def: &CylinderDef) -> Option<bool> {
        Some(!self.cut_value(k, def, self.standard.prec)?.sign()?)
    }

    /// `(n·o + d)² − r²|n|²` on plane `k`'s outward interval at `prec`.
    fn cut_value(&self, k: usize, def: &CylinderDef, prec: usize) -> Option<HpBounded> {
        let pl = self.outward_hp(k, prec);
        let o = def.origin().map(|x| HpBounded::of_rat(x, prec));
        let off = pl[0]
            .mul(&o[0], prec)
            .add(&pl[1].mul(&o[1], prec), prec)
            .add(&pl[2].mul(&o[2], prec), prec)
            .add(&pl[3], prec);
        let nn = pl[0]
            .mul(&pl[0], prec)
            .add(&pl[1].mul(&pl[1], prec), prec)
            .add(&pl[2].mul(&pl[2], prec), prec);
        let r2 = HpBounded::of_bigint(def.r2().numer(), prec)
            .div(&HpBounded::of_bigint(def.r2().denom(), prec), prec)?;
        Some(off.mul(&off, prec).sub(&r2.mul(&nn, prec), prec))
    }

    /// **Are two points proved apart?** — some coordinate of their interval boxes separated.
    /// `None` where a box cannot be formed; `Some(false)` is «not shown», never «one point».
    pub fn points_apart(&self, a: CylPoint<'_>, b: CylPoint<'_>) -> Option<bool> {
        let prec = self.standard.prec;
        let (x, y) = (self.point_box(a, prec)?, self.point_box(b, prec)?);
        Some((0..3).any(|k| y[k].sub(&x[k], prec).sign().is_some()))
    }

    /// **Is `n` strictly inside the counter-clockwise arc `lo → hi`**, three points of one circle of
    /// cylinder `def`, read from their interval boxes: travelling counter-clockwise about the axis
    /// `m`, `(lo, n, hi)` is in that order exactly when `m · ((n − lo) × (hi − lo)) > 0` — the sign
    /// of the triangle's turn about the axis. `None` where the boxes cannot say (two names for one
    /// point among them, or too little precision).
    pub fn circle_turn(
        &self,
        def: &CylinderDef,
        lo: CylPoint<'_>,
        n: CylPoint<'_>,
        hi: CylPoint<'_>,
    ) -> Option<bool> {
        self.turn_value(def, lo, n, hi, self.standard.prec)?.sign()
    }

    /// `m · ((n − lo) × (hi − lo))` at `prec`.
    fn turn_value(
        &self,
        def: &CylinderDef,
        lo: CylPoint<'_>,
        n: CylPoint<'_>,
        hi: CylPoint<'_>,
        prec: usize,
    ) -> Option<HpBounded> {
        let (a, b, c) = (
            self.point_box(lo, prec)?,
            self.point_box(n, prec)?,
            self.point_box(hi, prec)?,
        );
        let m = def.dir().map(|x| HpBounded::of_rat(x, prec));
        let e1 = [
            b[0].sub(&a[0], prec),
            b[1].sub(&a[1], prec),
            b[2].sub(&a[2], prec),
        ];
        let e2 = [
            c[0].sub(&a[0], prec),
            c[1].sub(&a[1], prec),
            c[2].sub(&a[2], prec),
        ];
        let x = [
            e1[1].mul(&e2[2], prec).sub(&e1[2].mul(&e2[1], prec), prec),
            e1[2].mul(&e2[0], prec).sub(&e1[0].mul(&e2[2], prec), prec),
            e1[0].mul(&e2[1], prec).sub(&e1[1].mul(&e2[0], prec), prec),
        ];
        Some(
            m[0].mul(&x[0], prec)
                .add(&m[1].mul(&x[1], prec), prec)
                .add(&m[2].mul(&x[2], prec), prec),
        )
    }
}

/// An interval's sign as `±1`, `None` where it straddles zero.
fn sign_i8(v: &HpBounded) -> Option<i8> {
    v.sign().map(|pos| if pos { 1 } else { -1 })
}

fn negate(x: &HpBounded, prec: usize) -> HpBounded {
    HpBounded::of_rat(Rat::from_int(0), prec).sub(x, prec)
}

impl PiercePoint {
    /// The exact pair, borrowed.
    pub fn exact(&self) -> Option<(&MeetLine, &QuadVal)> {
        match self {
            Self::Exact(line, s) => Some((line, s)),
            Self::Bounded { .. } => None,
        }
    }

    /// The exact pair, moved out.
    pub fn into_exact(self) -> Option<(MeetLine, QuadVal)> {
        match self {
            Self::Exact(line, s) => Some((line, s)),
            Self::Bounded { .. } => None,
        }
    }

    /// **The side of plane `q` this point lies on, from `q`'s statement** — one `a + b√c` sign
    /// ([`nacre_exact::quad::plane_side`]) against `q`'s world coefficients, bridged from the name's
    /// facing to the outward one ([`Judge::outward_fix`]), so it answers in the outward convention
    /// `Judge::orient3d` gives a three-plane point. `None` for a named plane with no narrow
    /// rational description. Where the point or `q` has no world statement the sign is read off
    /// interval boxes — `q`'s outward plane ([`Judge::outward_hp`]) at the point's box, `None` where
    /// the box straddles `q` — and a bounded point is on its own two planes by name.
    pub fn side<W: PlaneWitness>(&self, jd: &Judge<'_, W>, q: usize) -> Option<i8> {
        match self {
            Self::Bounded { planes, .. } if planes.contains(&q) => Some(0),
            // A plane with no world statement has no coefficients for the exact sign; the point's
            // own box is read against its interval plane instead.
            Self::Bounded { .. } => sign_i8(&self.side_value(jd, q, jd.standard.prec)?),
            Self::Exact(..) if jd.planes[q].world_name().is_none() => {
                sign_i8(&self.side_value(jd, q, jd.standard.prec)?)
            }
            Self::Exact(line, s) => {
                let co = jd.world_rat(q)?;
                let fix = jd.outward_fix(q)?;
                Some(
                    fix * match nacre_exact::quad::plane_side(&co, line, s) {
                        Orient::Positive => 1,
                        Orient::Negative => -1,
                        Orient::Zero => 0,
                    },
                )
            }
        }
    }

    /// `n_q · x + d_q` on `q`'s outward interval plane at `prec`.
    fn side_value<W: PlaneWitness>(
        &self,
        jd: &Judge<'_, W>,
        q: usize,
        prec: usize,
    ) -> Option<HpBounded> {
        let at = self.point_hp_at(jd, prec)?;
        let pl = jd.outward_hp(q, prec);
        Some(
            pl[0]
                .mul(&at[0], prec)
                .add(&pl[1].mul(&at[1], prec), prec)
                .add(&pl[2].mul(&at[2], prec), prec)
                .add(&pl[3], prec),
        )
    }

    /// **This point as an interval** at the operation's precision ([`Self::point_hp_at`]).
    pub fn point_hp<W: PlaneWitness>(&self, jd: &Judge<'_, W>) -> Option<[HpBounded; 3]> {
        self.point_hp_at(jd, jd.standard.prec)
    }

    /// **This point as an interval at `prec`** — `base + s·dir` along its exact meet line, the root
    /// realized at `prec` ([`nacre_exact::quad_hp`]); a point with no exact meet re-solved from its
    /// definition, its box at the operation's precision kept from its solve.
    pub fn point_hp_at<W: PlaneWitness>(
        &self,
        jd: &Judge<'_, W>,
        prec: usize,
    ) -> Option<[HpBounded; 3]> {
        match self {
            Self::Bounded { at, .. } if prec == jd.standard.prec => Some((**at).clone()),
            Self::Bounded {
                planes, root, def, ..
            } => jd.pierce_box(*planes, def, *root, prec),
            Self::Exact(line, s) => {
                let sb = nacre_exact::quad_hp(s, prec)?;
                let (b, d) = (line.base(), line.dir());
                let at = |k: usize| nacre_exact::affine_hp(b[k], d[k], &sb, prec);
                Some([at(0)?, at(1)?, at(2)?])
            }
        }
    }

    /// The name's two planes and its cylinder, for a point with no exact meet — what the by-name
    /// rules read.
    pub fn name_planes_and_cyl(&self) -> Option<([usize; 2], usize)> {
        match self {
            Self::Exact(..) => None,
            Self::Bounded { planes, cyl, .. } => Some((*planes, *cyl)),
        }
    }

    /// The name's two planes, for a point with no exact meet.
    pub fn name_planes(&self) -> Option<[usize; 2]> {
        self.name_planes_and_cyl().map(|(p, _)| p)
    }

    /// **Which of wall `wall`'s two rulings of `def` this point lies on** — the sign of
    /// `(x − o) · (m × n_wall)` with `n_wall` the wall's coefficients as a pierce name orders them
    /// ([`Judge::canon_hp`]). Exact for an exact point against a named wall
    /// ([`nacre_exact::quad::ruling_side_signed`]), `None` there for a wall with no narrow
    /// rational description (a `Wide` name — the exact road's population, which the interval road
    /// does not answer for, as [`Judge::pierce`] and [`Self::side`] rule); otherwise read off the
    /// interval box, `None` where it straddles the plane through the axis.
    pub fn ruling_side<W: PlaneWitness>(
        &self,
        jd: &Judge<'_, W>,
        wall: usize,
        def: &CylinderDef,
    ) -> Option<i8> {
        if let Some((line, s)) = self.exact()
            && jd.planes[wall].world_name().is_some()
        {
            return nacre_exact::quad::ruling_side_signed(&jd.world_rat(wall)?, def, (line, s));
        }
        sign_i8(&self.ruling_value(jd, wall, def, jd.standard.prec)?)
    }

    /// `(x − o) · (m × n_wall)` at `prec`.
    fn ruling_value<W: PlaneWitness>(
        &self,
        jd: &Judge<'_, W>,
        wall: usize,
        def: &CylinderDef,
        prec: usize,
    ) -> Option<HpBounded> {
        let at = self.point_hp_at(jd, prec)?;
        let n = jd.canon_hp(wall, prec);
        let m = def.dir().map(|x| HpBounded::of_rat(x, prec));
        let o = def.origin().map(|x| HpBounded::of_rat(x, prec));
        let c = [
            m[1].mul(&n[2], prec).sub(&m[2].mul(&n[1], prec), prec),
            m[2].mul(&n[0], prec).sub(&m[0].mul(&n[2], prec), prec),
            m[0].mul(&n[1], prec).sub(&m[1].mul(&n[0], prec), prec),
        ];
        let d = [
            at[0].sub(&o[0], prec),
            at[1].sub(&o[1], prec),
            at[2].sub(&o[2], prec),
        ];
        Some(
            d[0].mul(&c[0], prec)
                .add(&d[1].mul(&c[1], prec), prec)
                .add(&d[2].mul(&c[2], prec), prec),
        )
    }
}

#[cfg(test)]
#[path = "tests/cylinder.rs"]
mod tests;
