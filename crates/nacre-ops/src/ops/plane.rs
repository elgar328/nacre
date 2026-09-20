use super::*;
/// A sketch-plane frame: a 2-D point `(u, v)` maps to `origin + u·x + v·y`.
/// The axes are unit and orthogonal (the constructors ensure it); the normal is `x × y`.
///
/// ★★★★★ **The fields are private, and that is the whole point.** They used to be `pub`, so a
/// caller handed the kernel three *normalized* f64 vectors — and normalizing is where the
/// exactness dies: a plane with normal `(1, 1, 1)` has coefficients `[1, 1, 1, 0]`, three
/// integers, but its unit axes square to `0.9999999999999999…` and no exact form survives. The
/// kernel then had nothing to build on and dropped the whole prism to f64.
///
/// So a plane is built through a constructor that **keeps what the caller stated**
/// ([`PlaneDef`]), and the axes below are the *realization* of that. The two cannot describe
/// different planes because only one of them is written down.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchPlane {
    // `pub(crate)`, not `pub`: the constructors live in this crate's `lib.rs`, and what has to be
    // closed is the **public** surface — a caller outside must not be able to hand in three
    // normalized axes and call that a plane.
    pub(crate) origin: Point3,
    pub(crate) x_axis: Vector3,
    pub(crate) y_axis: Vector3,
    /// What the caller stated, exactly — `None` when they stated only axes (a frame the kernel
    /// cannot reconstruct, which then takes the f64 path it always took).
    pub(crate) def: Option<PlaneDef>,
}

impl SketchPlane {
    /// The world XY plane: `+u = x̂`, `+v = ŷ`, normal `+ẑ`.
    pub fn world_xy() -> Self {
        Self::axis_plane([0, 0, 1], [1, 0, 0], [0, 1, 0])
    }

    /// The world YZ plane: `+u = ŷ`, `+v = ẑ`, normal `+x̂`.
    pub fn world_yz() -> Self {
        Self::axis_plane([1, 0, 0], [0, 1, 0], [0, 0, 1])
    }

    /// The world ZX plane: `+u = ẑ`, `+v = x̂`, normal `+ŷ`.
    ///
    /// ★★ **The axes are named, not derived.** `ẑ × n` would give `−x̂` here; the convention a
    /// person expects (and the one the script layer documents) is `+u = ẑ`. A named plane gets to
    /// say, which is exactly what [`PlaneDef::ref_dir`] is for.
    pub fn world_zx() -> Self {
        Self::axis_plane([0, 1, 0], [0, 0, 1], [1, 0, 0])
    }

    /// One of the three world planes, stated exactly: normal, `+u`, `+v` as integer triples.
    /// The defining points are `[0, u, v]` — `u × v = n` for all three world planes, so the
    /// point order carries the same normal the coefficients used to state, and `points[1]`
    /// carries the *named* `+u` (the `world_zx` convention `+u = ẑ` included).
    fn axis_plane(n: [i128; 3], u: [i128; 3], v: [i128; 3]) -> Self {
        let r = |a: [i128; 3]| a.map(Rat::from_int);
        let f = |a: [i128; 3]| Vector3::from_array(a.map(|c| c as f64));
        debug_assert_eq!(
            {
                let (u, v) = (f(u), f(v));
                u.cross(v).as_array()
            },
            f(n).as_array(),
            "axis_plane point order must reproduce the stated normal"
        );
        Self {
            origin: Point3::origin(),
            x_axis: f(u),
            y_axis: f(v),
            def: Some(PlaneDef {
                points: [[Rat::from_int(0); 3], r(u), r(v)],
            }),
        }
    }

    /// A plane through `origin` with the given `normal`, its axes synthesized by the same
    /// convention a face's frame uses (`ops::frame_axes` — cross the world `ẑ` into the normal,
    /// or `ŷ` when the normal is vertical). `None` if `normal` is zero.
    ///
    /// ★ It has to be the same convention: this and [`face_plane`] answer the same question, and a
    /// caller that builds a frame here and compares it with one read off a face would otherwise
    /// find them ninety degrees apart.
    ///
    /// ★★★★ **`normal` is kept, not just consumed.** Normalizing it is what destroys the exact
    /// form — the caller's `(1, 1, 1)` is coefficients `[1, 1, 1, 0]`, while `normalize` of it
    /// squares to `0.9999999999999999…`. The unit axes below stay the f64 cache; the definition
    /// records what was handed in. A normal outside the decimal window simply leaves `def` empty,
    /// so this constructor is **never stricter than it was**.
    pub fn from_origin_normal(origin: Point3, normal: Vector3) -> Option<Self> {
        let n = normal.normalize()?;
        let (x_axis, y_axis) = frame_axes(n)?;
        Some(Self {
            origin,
            x_axis,
            y_axis,
            def: Self::normal_def(origin, normal),
        })
    }

    /// `from_origin_normal`'s exact half: the plane through `origin` with normal `normal`, and the
    /// same `ẑ × n` convention spelled in rationals (un-normalized — a cross product is already in
    /// the plane, so nothing needs projecting).
    ///
    /// The definition is the three points `[o, o + u, o + w]`, and **both in-plane directions are
    /// basis crosses** — `u = ẑ × n` (or `ŷ × n` for a vertical normal), `w = x̂ × n` (or a
    /// sibling), which are component *shuffles* of the written decimals: no products, no
    /// divisions, nothing to overflow. The only declines left are a coordinate outside the
    /// decimal window and a zero normal.
    ///
    /// ★★ **Polarity is algebraic**: `(a × n) × (b × n) = det[a, b, n] · n`, so `u × w` is a
    /// *known scalar* times `n` — `n₁` for the `(ẑ, x̂)` pair — and flipping `w`'s sign when that
    /// scalar is negative makes the point order face the caller's normal, both signs, exactly.
    ///
    /// ★★★ **Why not `w = n × u`, the "obvious" second direction: its components are products.**
    /// S6a shipped that and the checked arithmetic looked like a formality; the day the silent
    /// f64 fallback stopped absorbing failures (S6b), a proptest found the window's
    /// small-exponent corner — a `10²¹` denominator squares to `10⁴²`, and even the *primitive*
    /// direction of that cross needs 137 bits. The retired `named_plane_points` solved an axis
    /// for the same reason. Basis crosses stay inside the inputs' own widths.
    fn normal_def(origin: Point3, normal: Vector3) -> Option<PlaneDef> {
        let o = origin.as_array().map(Rat::from_decimal);
        let n = normal.as_array().map(Rat::from_decimal);
        let (o, n) = ([o[0]?, o[1]?, o[2]?], [n[0]?, n[1]?, n[2]?]);
        let zero = Rat::from_int(0);
        let neg = |x: Rat| {
            zero.checked_sub(x)
                .expect("negation cannot overflow a lifted decimal")
        };
        // (u, w0, s): two independent in-plane basis crosses and the scalar with
        // `u × w0 = s · n`, per the identity above (arms chosen so `s ≠ 0`).
        let (u, w0, s) = if n[0] == zero && n[1] == zero {
            if n[2] == zero {
                return None; // a zero normal names no plane
            }
            // Vertical: u = ŷ × n, w0 = x̂ × n; det[ŷ, x̂, n] = −n₂.
            ([n[2], zero, zero], [zero, neg(n[2]), zero], neg(n[2]))
        } else if n[1] != zero {
            // The stated convention: u = ẑ × n; w0 = x̂ × n; det[ẑ, x̂, n] = n₁.
            ([neg(n[1]), n[0], zero], [zero, neg(n[2]), n[1]], n[1])
        } else {
            // n₁ = 0, n₀ ≠ 0: u = ẑ × n; w0 = ŷ × n; det[ẑ, ŷ, n] = −n₀.
            ([neg(n[1]), n[0], zero], [n[2], zero, neg(n[0])], neg(n[0]))
        };
        let w = if s > zero { w0 } else { w0.map(neg) };
        let add = |a: [Rat; 3], b: [Rat; 3]| -> Option<[Rat; 3]> {
            Some([
                a[0].checked_add(b[0])?,
                a[1].checked_add(b[1])?,
                a[2].checked_add(b[2])?,
            ])
        };
        Some(PlaneDef {
            points: [o, add(o, u)?, add(o, w)?],
        })
    }

    /// **A plane through three written points**: `origin` is the sketch's `(0, 0)`, `+u` runs
    /// toward `x_point`, and `+v` leans toward `y_hint`.
    ///
    /// ★★★ **Everything here is exact by construction.** The written points, lifted to their
    /// decimal truth, **are** the definition — origin first, so the sketch `(0, 0)` and the `+u`
    /// direction (`x_point − origin`) fall out of the structure. `None` if the three are
    /// collinear or fall outside the decimal window.
    pub fn through_points(origin: Point3, x_point: Point3, y_hint: Point3) -> Option<Self> {
        let x = (x_point - origin).normalize()?;
        let v = y_hint - origin;
        let y = (v - x * v.dot(x)).normalize()?;
        let lift = |p: Point3| {
            let a = p.as_array().map(Rat::from_decimal);
            Some([a[0]?, a[1]?, a[2]?])
        };
        let def = (|| {
            let (o, xp, yh) = (lift(origin)?, lift(x_point)?, lift(y_hint)?);
            // The written points ARE the definition; the only thing to verify is that they name
            // a plane at all. `plane_name_exact` is total (Narrow | Wide), so `None` means
            // exactly one thing: collinear. There is no
            // "answer does not fit i128" failure class.
            nacre_scalar::plane_name_exact(o, xp, yh)?;
            Some(PlaneDef {
                points: [o, xp, yh],
            })
        })();
        Some(Self {
            origin,
            x_axis: x,
            y_axis: y,
            def,
        })
    }

    /// The same plane **moved to pass through `p`**, with `p` as the sketch's `(0, 0)`.
    ///
    /// ★★★★ **The plane travels with the origin** — `plane(ZX, { origin: … })` sets the position
    /// as well as the 2-D origin, and `SketchPlane { origin, ..world_xy() }` always meant that.
    /// Keeping the plane in place while moving the origin would leave the definition describing
    /// one plane and its origin sitting on another: measured, `world_xy().with_origin([0, 0, 0.5])`
    /// recorded `z = 0` for a cap at `z = 0.5`, and a boolean built on that lost 0.04 of volume.
    /// Since the definition is three points and the origin is the first of them, the move is a
    /// translation of the whole triple: differences (`ref_dir`) and the normal are untouched,
    /// exactly the "translation does not turn `+u`" the old form promised.
    pub fn with_origin(mut self, p: Point3) -> Self {
        self.origin = p;
        self.def = self.def.and_then(|d| {
            let a = p.as_array().map(Rat::from_decimal);
            let origin = [a[0]?, a[1]?, a[2]?];
            let shift = |q: [Rat; 3]| -> Option<[Rat; 3]> {
                Some([
                    q[0].checked_sub(d.points[0][0])?.checked_add(origin[0])?,
                    q[1].checked_sub(d.points[0][1])?.checked_add(origin[1])?,
                    q[2].checked_sub(d.points[0][2])?.checked_add(origin[2])?,
                ])
            };
            Some(PlaneDef {
                points: [origin, shift(d.points[1])?, shift(d.points[2])?],
            })
        });
        self
    }

    /// A frame from axes the caller already holds — **and the axes' decimal truth is its
    /// definition**. The boundary rule that `Profile2d` applies to coordinates applies to
    /// axes too: what the caller wrote *is* the statement, so the plane through
    /// `[o, o + x, o + y]` and the `+u` direction `x` are recorded exactly. A 45°-rotated frame
    /// — whose axes never lift to exact orthonormal rationals — extrudes through the frame
    /// road (wide names included) instead of falling silently to f64.
    ///
    /// ★ The **world lift** still comes first: axes whose decimals square and cross to exact
    /// `1`/`0` — a Pythagorean frame like `(0.6, 0.8, 0)`/`(−0.48, 0.36, 0.8)` — pass `exact()`
    /// and take the world-rational path, definition or not. That population is how the
    /// `n·n`-overflow walls (the census `wf` family) are built.
    ///
    /// ★★ **What the definition states — and what it does not.** Three points, `+u`, and the
    /// polarity (point order); the realized frame is *orthonormal*, exactly as for every other
    /// definition (`ref_dir` is any length, `+v` is derived on `y`'s side). Exact decimal
    /// orthogonality is deliberately **not** required — a rotated pair `[c, s, 0]/[−s, c, 0]`
    /// cancels to exactly zero, but two independently rounded axes need not, and requiring it
    /// would strand exactly the callers this lift exists for. A skewed (non-orthogonal) pair is
    /// outside [`SketchPlane`]'s contract; the frame realization drops the skew.
    ///
    /// `def` stays `None` only for axes outside the decimal window or a degenerate pair
    /// (`x × y = 0` in the decimals) — those keep the f64 path they had.
    pub fn from_axes(origin: Point3, x_axis: Vector3, y_axis: Vector3) -> Self {
        let def = (|| {
            let lift = |p: [f64; 3]| {
                let a = p.map(Rat::from_decimal);
                Some([a[0]?, a[1]?, a[2]?])
            };
            let o = lift(origin.as_array())?;
            let add = |a: [Rat; 3], b: [Rat; 3]| -> Option<[Rat; 3]> {
                Some([
                    a[0].checked_add(b[0])?,
                    a[1].checked_add(b[1])?,
                    a[2].checked_add(b[2])?,
                ])
            };
            let px = add(o, lift(x_axis.as_array())?)?;
            let py = add(o, lift(y_axis.as_array())?)?;
            // Total: `None` means exactly one thing — the axes are parallel (or zero) in their
            // decimal truth, and name no plane.
            nacre_scalar::plane_name_exact(o, px, py)?;
            Some(PlaneDef {
                points: [o, px, py],
            })
        })();
        Self {
            origin,
            x_axis,
            y_axis,
            def,
        }
    }

    /// The sketch's `(0, 0)` in space.
    #[inline]
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    /// The `+u` direction.
    #[inline]
    pub fn x_axis(&self) -> Vector3 {
        self.x_axis
    }

    /// The `+v` direction.
    #[inline]
    pub fn y_axis(&self) -> Vector3 {
        self.y_axis
    }

    /// The 3-D point for sketch coordinates `p = (u, v)`.
    #[inline]
    pub fn point(&self, p: Point2) -> Point3 {
        self.origin + self.x_axis * p[0] + self.y_axis * p[1]
    }

    /// The plane normal `x × y` (unit when the axes are unit and orthogonal).
    #[inline]
    pub fn normal(&self) -> Vector3 {
        self.x_axis.cross(self.y_axis)
    }
}

/// **A sketch plane as its author stated it** — the exact truth behind [`SketchPlane`]'s f64 axes.
///
/// ★★★ **One field, and the invariants are structural.** This used to carry coefficients, an
/// origin, and a reference direction as three halves that every constructor had to keep agreeing
/// ("an origin that is not on `coeffs` is a definition describing two different planes", which
/// cost 0.04 of volume the one time it happened). Now the definition is the three points alone:
///
/// - the sketch's `(0, 0)` **is** `points[0]`,
/// - `+u` **is** `points[1] − points[0]` — a difference of two points of the plane, so it lies in
///   the plane by definition,
/// - the normal's direction is `(p1 − p0) × (p2 − p0)` — the point order carries the polarity.
///
/// Nothing is left to check, and nothing can disagree. The canonical coefficients are *derived*
/// (`nacre_scalar::plane_name_exact` — total, `Narrow | Wide`), so there is no
/// failure class "the coefficients do not fit `i128`": three in-window points
/// always name their plane, however wide its canonical form.
///
/// ★ `ref_dir()` is **not** a unit vector and is not projected; the normalization a frame needs
/// is exactly one `1/√(rational)` at realization time — never something stored.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneDef {
    /// Three points of the plane, non-collinear (`plane_name_exact` is what verified it — every
    /// constructor rejects a collinear triple as "no plane"). `points[0]` is the sketch origin,
    /// `points[1] − points[0]` the `+u` direction, and the order fixes the normal's sign.
    pub(crate) points: [[nacre_scalar::Rat; 3]; 3],
}

impl PlaneDef {
    /// The three defining points — the sketch origin first, then the point `+u` runs toward,
    /// then the point fixing the normal's side.
    pub fn points(&self) -> [[nacre_scalar::Rat; 3]; 3] {
        self.points
    }

    /// Where the sketch's `(0, 0)` sits — the first defining point.
    pub fn origin(&self) -> [nacre_scalar::Rat; 3] {
        self.points[0]
    }

    /// The `+u` direction, in the plane, not unit length — `points[1] − points[0]`.
    ///
    /// The subtraction cannot overflow: both points passed through a constructor, and every
    /// constructor either lifted decimals (narrow) or added one lifted decimal to another —
    /// widths nowhere near `i128`'s ceiling.
    pub fn ref_dir(&self) -> [nacre_scalar::Rat; 3] {
        core::array::from_fn(|i| {
            self.points[1][i]
                .checked_sub(self.points[0][i])
                .expect("constructor-bounded widths")
        })
    }
}
