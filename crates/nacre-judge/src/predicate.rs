//! Toleranced geometric predicates: the plane-arrangement sign predicates routed through
//! the CIP kernel so they stay exact under rotation.
//!
//! A rotated face's plane coefficients and `tri` coordinates are rounded irrationals, so the
//! exact predicates (`nacre-predicates`) over them are exact only w.r.t. the *rounded*
//! geometry. When a predicate's planes are rotated, these wrappers rebuild each plane from the
//! three exact [`WitnessPoint`] its face carries ([`Witness::tri_pt3`], or — for the axis-aligned
//! operand of a *mixed*-rotation boolean, whose `tri_pt3` is `None` — exactly from its `tri`
//! coordinates, see [`plane_def`]) and decide the sign with the [`crate::kernel`] judges
//! instead.
//!
//! **Routing is per-predicate, derived — no `rotated` flag is threaded.** Each wrapper asks
//! [`any_rotated`] of just the planes it touches: all-axis-aligned → the exact hot path (never
//! builds a `WitnessPoint`); any rotated → the kernel judge. So a mixed-rotation boolean keeps the
//! axis-aligned operand's own predicates on the fast path, finer than a per-boolean flag.
//!
//! The caller (`nacre-ops`) provides the witnesses by implementing [`Witness`] (a face) and
//! [`PlaneWitness`] (a plane class) on its own tables — the port keeps this crate free of the
//! b-rep and geometry types — and pairs a table with a [`Judge`], which is where the three facts
//! that belong to the *operation* live: the standard of proof, the collector for what could not
//! be proved, and the table itself. The predicates are its methods.

use crate::kernel::frame3::{
    Decision, MoveNode, Standard, WitnessPoint, chain_parity, cramer_iv, dir_sign_judge,
    indirect_cmp_coord_judge, orient3d_from_cramer, orient3d_judge,
};
use nacre_exact::{Orient, Rat};
use nacre_math::Point3;
use nacre_predicates::{
    ThreePlane, det3_sign, indirect_cmp_coord, indirect_plane_side, orient2d, orient3d,
};
use num_bigint::BigInt;

/// **What a judgement was asked about**, in the only vocabulary this crate has: plane-table
/// indices.
///
/// Turning these into vertices, faces and solids is the caller's job and costs a reverse lookup
/// this crate cannot do — so the report starts with indices and distances, and grows names only
/// once it has proved its worth.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Site {
    /// Two faces judged to lie on **the same plane** — the merge that decides the plane classes,
    /// and therefore everything downstream. First in any report for that reason.
    PlanesCoplanar { i: usize, j: usize },
    /// Which side of plane `j` the implicit point `∩(p, q, r)` lies on.
    Orient3d {
        p: usize,
        q: usize,
        r: usize,
        j: usize,
    },
    /// The order of two implicit points along one axis.
    CmpCoord {
        a: [usize; 3],
        b: [usize; 3],
        axis: usize,
    },
    /// How the line `p ∩ a` runs relative to plane `b`.
    DirSign { p: usize, a: usize, b: usize },
}

/// One judgement that did **not** come back with a proved sign, and what it did establish.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Evidence {
    pub site: Site,
    pub outcome: Decision,
}

/// Where an operation's [`Evidence`] is collected — **diagnostic only**.
///
/// Nothing read from here changes a sign, a merge, or a coordinate: the geometry consumes
/// [`Decision::orient`]'s `i8` exactly as before, and this rides alongside. That is what lets the
/// port be added without touching a single geometric decision.
///
/// **It is not a global.** The collector belongs to the operation and reaches the predicates the
/// same way the judging precision does — through the witness table the caller already holds. A
/// thread-local would have put permanent mutable state in a pure numeric crate, and would have
/// been wrong the moment two operations ran on one thread.
#[derive(Clone, Debug, Default)]
pub struct Notes(NotesCell);

#[cfg(feature = "parallel")]
type NotesCell = std::sync::Arc<std::sync::Mutex<Vec<Evidence>>>;
#[cfg(not(feature = "parallel"))]
type NotesCell = std::rc::Rc<std::cell::RefCell<Vec<Evidence>>>;

impl Notes {
    pub fn new() -> Notes {
        Notes::default()
    }

    /// Record one inconclusive judgement. Clones share one list, so a collector may be handed
    /// around freely and still collect into one place.
    pub fn push(&self, e: Evidence) {
        #[cfg(feature = "parallel")]
        self.0.lock().expect("notes lock").push(e);
        #[cfg(not(feature = "parallel"))]
        self.0.borrow_mut().push(e);
    }

    /// Everything recorded, in a **deterministic** order.
    ///
    /// Sorted by site, not by arrival: predicates may run in any order (and, under `parallel`, on
    /// any thread), and a report that changed shape with the schedule would be a poor thing to
    /// hand a user who is trying to reproduce a result.
    pub fn sorted(&self) -> Vec<Evidence> {
        #[cfg(feature = "parallel")]
        let mut v = self.0.lock().expect("notes lock").clone();
        #[cfg(not(feature = "parallel"))]
        let mut v = self.0.borrow().clone();
        v.sort_by_key(|e| e.site);
        v
    }
}

/// A plane witnessed by three points known to lie on it, and — when the solid was rotated —
/// their exact [`WitnessPoint`] definitions.
///
/// The rotation-general predicates need only this, which is why one implementation can serve
/// both index spaces (a plane class and a single face) without confusing them. Only
/// [`Judge::planes_coplanar`] uses the face form — it is the predicate that *defines* the classes,
/// so it necessarily runs before a plane table exists.
pub trait Witness {
    fn tri(&self) -> [Point3; 3];
    /// The three `tri` points as exact [`WitnessPoint`] definitions, **always present**.
    ///
    /// This is a *cache*: the definition is built once, where the witness is, and every predicate
    /// borrows it. It is not an `Option` whose emptiness *also* means "not rotated": one field
    /// answering two questions could not be filled in advance without changing which predicate
    /// path runs. [`Witness::is_rotated`] is that second question.
    fn tri_pt3(&self) -> &[WitnessPoint; 3];
    /// Whether this plane came from a rotated solid — the predicate-routing signal, and a
    /// **different fact** from "does the definition carry a rotation chain". Neither
    /// `chain.is_empty()` nor `tol == 0` is equivalent to it (a rotated solid's face may witness
    /// its plane through points with no chain, and a 90°-family rotation has tol exactly 0), so
    /// the answer is carried, not derived.
    fn is_rotated(&self) -> bool;

    /// The plane's **exact rational coefficients in the frame its provenance names** — the world
    /// when unmoved, the pre-motion frame when moved. `None` when the producer recorded none.
    fn base_coeffs_rat(&self) -> Option<[nacre_exact::Rat; 4]> {
        None
    }

    /// Identifies the motion this witness's definition carries — `0` for none, and equal values
    /// **only** for structurally identical chains (same nodes, order and pivots).
    ///
    /// A motion preserves every determinant these predicates take *up to its own determinant*, so
    /// when all of a judgement's inputs carry one motion the answer is the answer on their
    /// pre-motion data — exactly, with no tolerance at all. This is what lets such a judgement
    /// leave the toleranced path entirely.
    fn chain_id(&self) -> u64;

    /// The witness triangle in the **pre-motion frame**, or `None` when those coordinates are not
    /// `f64`-representable and so cannot be handed to the exact predicate.
    ///
    /// **★ It is the pre-motion triangle *canonicalised for handedness*, not the raw one — and
    /// that is the implementor's job.** A chain containing an odd number of reflections is
    /// improper (`det = −1`), and the "up to its own determinant" above then bites: every
    /// judgement's inputs carry the *same* chain, so the determinant flips **uniformly** and the
    /// shortcut returns a confidently wrong sign rather than a conservative miss. An implementor
    /// therefore reflects the pre-motion data once more when the parity is odd
    /// ([`crate::chain_parity`]), which restores the moved frame's handedness and makes every
    /// determinant question transfer unchanged. `nacre-ops`' `BaseFrame::of` is the reference
    /// implementation; `crate::frame3`'s `shared_base` applies the same convention (a sign flip
    /// on x, exact for every finite `f64`) to the data it derives itself.
    ///
    /// **And a plane must be derived from the corrected points, never corrected separately** —
    /// reflecting a triangle and re-deriving its normal differs from reflecting the normal by a
    /// global sign, because the cross product is a pseudovector, and
    /// [`Judge::plane_pair_dir_sign`] reads exactly that sign.
    fn base_tri(&self) -> Option<[Point3; 3]>;
}

/// A witness that additionally carries its plane's exact coefficients — what the plane-class
/// predicates ([`Judge::orient3d`], [`Judge::cmp_coord`], [`Judge::plane_pair_dir_sign`]) need on the exact
/// path. (The stored↔outward `frame_sign` [`Judge::plane_pair_dir_sign`] also uses is *derived* from
/// `coeffs` + `tri` by [`PlaneWitness::frame_sign`], not required from the impl.) A single face
/// (which never
/// plays a plane-class role) implements only [`Witness`].
pub trait PlaneWitness: Witness {
    /// The plane cache's (un-normalized) coefficients `[a, b, c, d]` (`n·x + d = 0`) — **raw**: a
    /// rounded image of the plane, not a description a predicate may decide on.
    ///
    /// ★★★ **Raw means "not known to describe the plane".** `d` is an `f64` product, so a face at
    /// `y = −0.2` gets a coefficient plane `2⁻⁵⁴` from the one its own points span, and a
    /// predicate that answers one question from here and the next from another description is
    /// describing two planes — answers composed across them are not even an order. Predicates read
    /// [`exact_coeffs`](Self::exact_coeffs) / [`exact_normal`](Self::exact_normal); this is for a
    /// caller that means to *relate* descriptions (the `cfg(test)` oracles below).
    fn coeffs(&self) -> [f64; 4];

    /// **The plane itself, in `f64`** — its canonical name (derived from its defining points
    /// without rounding) as a row, for an unmoved plane whose name fits 53 bits a coefficient;
    /// `None` for a rotated plane, an unnamed one, or a wider name.
    ///
    /// ★ No other description is consulted: the witness triangle [`tri`](Witness::tri) is the face
    /// corners' `f64` caches and need not lie on this plane (`z = 0.1` is not `fl(0.1)`), so a
    /// predicate asks every plane of an exact question by this row.
    fn exact_coeffs(&self) -> Option<[f64; 4]>;

    /// The **normal** of the same row under the same rule — standing where only `d` is too wide,
    /// so a question that reads no `d` (a direction) keeps its exact route there.
    fn exact_normal(&self) -> Option<[f64; 3]>;

    /// `+1` when the stored normal agrees with the witness triangle's, `-1` when they oppose —
    /// **the relation between the two descriptions**, not a choice between them.
    ///
    /// This is the one thing a caller legitimately wants both descriptions for, and the arrangement
    /// already computes and stores it. It is exposed rather than re-derived here so there is one
    /// spelling of it.
    fn frame_sign(&self) -> i8;

    /// The plane's coefficients in the pre-motion frame — **derived from the canonicalised
    /// `base_tri`**, with all the handedness caveats there.
    fn base_coeffs(&self) -> Option<[f64; 4]>;

    /// The plane's **name integers, folded to the stored orientation** — built once by
    /// [`name_stored_ints`] where the witness is constructed, `None` when the plane has no name
    /// (or the fold declined). It serves twice. Its `f64` row is what an implementor hands out as
    /// [`Self::exact_coeffs`] for an unmoved plane. And it gives a *wide* name — no `f64` row —
    /// its exact shortcuts back through the BigInt rescue: [`Witness::base_coeffs_rat`] and
    /// [`Self::base_coeffs`] read `PlaneName::narrow()` or `f64`, so without it a wide name answers
    /// `None` everywhere and the judgement climbs — **2.0×** the
    /// escalations (240 against 475; the confound
    /// `wide_datum_cost.rs` states — the two arms pick different vertex triples — applies to
    /// either number, so read it as "several times", not as a coefficient).
    fn name_ints(&self) -> Option<&NameInts> {
        None
    }
}

/// A plane's canonical name integers, **oriented as the stored coefficients are** — the exact
/// integer stand-in for [`PlaneWitness::coeffs`], any width.
///
/// The canonical name deliberately carries no direction (first nonzero coefficient positive),
/// and the stored plane may hold it either way round — measured, 246 of the census corpus's 3,125
/// interning hits arrive `flipped`. So the raw name cannot be handed to a direction-sensitive
/// predicate; the σ that relates the two is folded in **here, once, at construction**
/// ([`name_stored_ints`]), so every consumer inherits the stored convention and the existing
/// `frame_sign` bridges apply verbatim.
#[derive(Clone, Debug)]
pub struct NameInts {
    /// The integer coefficients, in the frame the name speaks (the world when unmoved, the
    /// pre-motion frame when moved), sign-folded to the stored orientation.
    pub ints: [BigInt; 4],
    /// Whether the canonical name needed the wide vessel — the rescue gate reads this: a
    /// question whose planes are all narrow keeps its existing routes (behavior-identical),
    /// so the corpus (`wide_planes` 0) is untouched.
    pub wide: bool,
    /// The same row in `f64`, `Some` only when **every** coefficient fits 53 bits — then each `f64`
    /// *is* its integer, and a predicate over the row answers for the plane the defining points
    /// span (the name is derived from them without rounding). This is what an unmoved plane's
    /// exact shortcuts read ([`PlaneWitness::exact_coeffs`]). A wider integer keeps the BigInt
    /// rescue: an `Expansion` over it overflows in `cmp_coord`'s degree-6 products, and
    /// `nacre_exact::nearest_f64_big_exact` alone would call `2⁶⁰` exact.
    pub row: Option<[f64; 4]>,
    /// The first three under the same rule — what a direction question reads, which stands where
    /// `d` alone is too wide.
    pub normal: Option<[f64; 3]>,
}

impl NameInts {
    /// The sign of one coefficient: `+1`, `0`, `-1`.
    ///
    /// ★ A door, so a consumer above this crate can ask which way the stored-oriented name points
    /// without touching `BigInt`. `nacre-ops` reads it to check its own canonical→stored turn
    /// against this exact one; giving it the integers instead would put an arbitrary-precision
    /// type into a crate that deliberately has none.
    pub fn coeff_sign(&self, k: usize) -> i8 {
        match self.ints[k].sign() {
            num_bigint::Sign::Plus => 1,
            num_bigint::Sign::Minus => -1,
            num_bigint::Sign::NoSign => 0,
        }
    }
}

/// Fold a plane's canonical name to the **stored orientation** — the one-time σ computation
/// [`PlaneWitness::name_ints`] carries.
///
/// The witness triangle's `base` points lie exactly on the named plane (the table contract:
/// three exact points on the plane, wound to the face's outward normal), so the name's normal
/// and the triangle's cross product are exactly parallel and their dot's sign is σ against the
/// *triangle's* orientation — computed in integers after clearing all nine coordinate
/// denominators by one common positive factor (a global scale moves neither the cross's
/// direction nor the dot's sign; per-point scales would). `frame_sign` (stored vs. triangle)
/// then carries it the rest of the way: `σ = sign(name·cross) · frame_sign`.
///
/// `None` when the plane has no name, or the dot is zero (a degenerate witness — collinear
/// `base` points span no direction to compare against), in which case the rescue simply
/// declines and the judgement keeps its toleranced route: slower, never wrong.
pub fn name_stored_ints(
    name: Option<&nacre_exact::PlaneName>,
    tri_pt3: &[WitnessPoint; 3],
    frame_sign: i8,
) -> Option<NameInts> {
    use num_integer::Integer;
    let name = name?;
    let bases: [&[Rat; 3]; 3] = [&tri_pt3[0].base, &tri_pt3[1].base, &tri_pt3[2].base];
    // One common positive scale for all nine coordinates.
    let lcm = bases
        .iter()
        .flat_map(|p| p.iter())
        .fold(BigInt::from(1), |l, r| l.lcm(&BigInt::from(r.denom())));
    let lift = |p: &[Rat; 3]| -> [BigInt; 3] {
        core::array::from_fn(|i| BigInt::from(p[i].numer()) * (&lcm / BigInt::from(p[i].denom())))
    };
    let (p0, p1, p2) = (lift(bases[0]), lift(bases[1]), lift(bases[2]));
    let edge = |q: &[BigInt; 3]| -> [BigInt; 3] { core::array::from_fn(|i| &q[i] - &p0[i]) };
    let (u, v) = (edge(&p1), edge(&p2));
    let cross = [
        &u[1] * &v[2] - &u[2] * &v[1],
        &u[2] * &v[0] - &u[0] * &v[2],
        &u[0] * &v[1] - &u[1] * &v[0],
    ];
    let mut ints = name.coeff_ints();
    // The σ below is exact only when the witness bases lie exactly on the named plane — i.e.
    // when witness and name **speak the same frame**. A nonzero residual here is not a broken
    // table: a named plane whose meets are wider than any witness base is
    // witnessed by its own frame's probes, whose bases are frame-local coordinates — a
    // different frame from the name's. The fold then declines, `name_ints` stays `None`, and
    // the plane keeps the toleranced routes: slower, never wrong.
    let on_plane = [&p0, &p1, &p2].iter().all(|p| {
        let r: BigInt = (0..3).map(|i| &ints[i] * &p[i]).sum::<BigInt>() + &ints[3] * &lcm;
        r.sign() == num_bigint::Sign::NoSign
    });
    if !on_plane {
        return None;
    }
    let dot: BigInt = (0..3).map(|i| &ints[i] * &cross[i]).sum();
    let sigma = match dot.sign() {
        num_bigint::Sign::Plus => frame_sign,
        num_bigint::Sign::Minus => -frame_sign,
        num_bigint::Sign::NoSign => return None,
    };
    if sigma < 0 {
        for c in &mut ints {
            *c = -&*c;
        }
    }
    let exact = |x: &BigInt| -> Option<f64> {
        if x.bits() > 53 {
            return None;
        }
        let (v, exact) = nacre_exact::nearest_f64_big_exact(x, &BigInt::from(1))?;
        exact.then_some(v)
    };
    let normal = (|| Some([exact(&ints[0])?, exact(&ints[1])?, exact(&ints[2])?]))();
    let row = normal.and_then(|[a, b, c]| Some([a, b, c, exact(&ints[3])?]));
    Some(NameInts {
        ints,
        wide: name.narrow().is_none(),
        row,
        normal,
    })
}

/// **One operation's judging**: the witnesses it reasons over, the standard it holds them to, and
/// where the evidence goes.
///
/// These three are properties of the *operation*, not of a plane, and this is what says so.
/// Stamped on every table row, they would need a two-phase construction (build the rows, then
/// stamp them), a placeholder for the gap between, and a guard for forgetting — a whole class of
/// mistake that exists only when a fact is stored somewhere it does not belong.
///
/// The witness table stays a **pure description** of geometry, which is what makes one
/// implementation able to serve both index spaces (a plane class and a single face).
pub struct Judge<'a, W> {
    /// The witness table these indices name.
    pub planes: &'a [W],
    /// How deeply to realize, and how close counts as one thing.
    pub standard: Standard,
    /// Where an inconclusive judgement's evidence is collected — diagnostic only; nothing read
    /// from here changes a sign, a merge or a coordinate.
    pub notes: &'a Notes,
    /// ★ **Each plane's interval coefficients, built on first use and borrowed thereafter.**
    ///
    /// The same cache shape [`Witness::tri_pt3`] already documents, one level up: the filter every
    /// certified judgement runs takes the planes as intervals, and those are a function of the
    /// definitions alone. Rebuilding them per call is ~20-24% of a certified judgement, and the
    /// arrangement's crossing collector hands in the same three definitions hundreds of times.
    ///
    /// ★ **Why it is here and not on the witness, where [`Witness::tri_pt3`] says caches go.**
    /// [`PlaneWitness`] is a trait its consumers implement. A method on it returning `[Bounded; 4]`
    /// would ask every implementor to *produce* intervals — to supply radii that actually bound
    /// the coefficients — and the soundness contract that type carries would leave this crate
    /// with it. The type itself is public vocabulary (`nacre_exact::Bounded`, beside `Mag` and
    /// `Rat`); the obligation to mint one correctly stays behind this crate's own constructors.
    /// `frame_sign` and `coeffs` are `i8` and `[f64; 4]`, so they *do* live on the witness; this
    /// one cannot follow them.
    ///
    /// Two workers racing to fill one cell compute the same value, so the answer does not depend on
    /// who won — the same argument `HpCell` rests on.
    iv: Vec<BoundedCell>,
}

/// Lazily-filled cell for one plane's interval coefficients — `OnceLock` under `parallel` because
/// the boolean hands every worker the same `&Judge`, `OnceCell` otherwise.
#[cfg(feature = "parallel")]
type BoundedCell = std::sync::OnceLock<[nacre_exact::Bounded; 4]>;
#[cfg(not(feature = "parallel"))]
type BoundedCell = std::cell::OnceCell<[nacre_exact::Bounded; 4]>;

impl<'a, W> Judge<'a, W> {
    pub fn new(planes: &'a [W], standard: Standard, notes: &'a Notes) -> Judge<'a, W> {
        Judge {
            planes,
            standard,
            notes,
            iv: (0..planes.len()).map(|_| BoundedCell::new()).collect(),
        }
    }

    /// The sign the geometry consumes, **and** a note of what backed it when it was not a proved
    /// one.
    ///
    /// A proved sign — including a proved zero — is the ordinary case and says nothing worth
    /// reporting. A coincidence, an exhausted judgement or a degenerate witness is a statement
    /// about *this model* that the caller cannot recover afterwards, because by then it is just a
    /// `0`.
    fn record(&self, site: Site, d: Decision) -> i8 {
        if !matches!(d, Decision::Sign(_)) {
            self.notes.push(Evidence { site, outcome: d });
        }
        to_i8(d.orient())
    }
}

/// The plane-class predicates — the questions the arrangement asks of a plane table.
///
/// Each routes itself: all-axis-aligned planes take the exact path and never build a `WitnessPoint`;
/// any rotated plane goes to the kernel judges, under this operation's [`Standard`], with what
/// it could not prove recorded in [`Judge::notes`].
impl<W: Witness> Judge<'_, W> {
    /// Plane `k`'s interval coefficients, built once and copied thereafter.
    ///
    /// Returns a copy rather than a borrow because `[Bounded; 4]` is four pairs of `f64` — cheaper to
    /// move than to keep a reference alive across the judge call, and it keeps the cell's borrow
    /// from outliving the lookup.
    fn plane_iv(&self, k: usize) -> [nacre_exact::Bounded; 4] {
        *self.iv[k].get_or_init(|| {
            let d = plane_def(self.planes, k);
            crate::kernel::frame3::plane_iv(&d[0], &d[1], &d[2])
        })
    }
}

impl<W: PlaneWitness> Judge<'_, W> {
    /// The sign of `orient3d(V, tri_j)` where `V = ∩(planes p, q, r)` is an implicit point,
    /// matching `order_along`'s shape (`+1`/`-1`/`0`).
    ///
    /// - `!rotated`: the exact path — which side of plane `j` the implicit point lies on, all four
    ///   planes given by their names' `f64` rows ([`PlaneWitness::exact_coeffs`]), bridged to
    ///   `tri`'s winding by `frame_sign`.
    /// - `rotated`: each of `p, q, r` and the explicit triangle `j` is taken as its three exact
    ///   [`WitnessPoint`] ([`plane_def`]), and [`crate::indirect_orient3d_judge`] decides the sign
    ///   from the
    ///   definitions — never materializing `V` or reading the rounded `tri`.
    ///
    /// Winding-invariant. A query plane `j` equal to one of `p, q, r` means the point lies on `j`,
    /// so the sign is exactly `0` (a combinatorial identity) — decided before either numeric
    /// branch, exact and cheap for both.
    pub fn orient3d(&self, p: usize, q: usize, r: usize, j: usize) -> i8 {
        match self.orient3d_cheap(p, q, r, j) {
            Some(o) => o,
            None => self.orient3d_given(self.cramer_of(p, q, r), p, q, r, j),
        }
    }

    fn orient3d_cheap(&self, p: usize, q: usize, r: usize, j: usize) -> Option<i8> {
        if j == p || j == q || j == r {
            return Some(0);
        }
        // ★ **`j` is asked by its coefficients, not by its triangle.** The rows are the planes'
        // names in `f64`, derived from their defining points without rounding; `tri` is the face's
        // corners' `f64` caches, which a decimal plane does not pass through (`z = 0.1` is not
        // `fl(0.1)`). Asking `j` by `tri` answered for the rounded model — measured, a box whose
        // corner lies on a wall exactly (`3·0.1 = 0.3`) was judged off it and the common refused as
        // `ZeroLengthEdge`. `frame_sign` carries `j`'s stored normal to its triangle's winding, the
        // bridge the shared-motion arm below already uses.
        if let (Some(cp), Some(cq), Some(cr), Some(cj)) = (
            self.planes[p].exact_coeffs(),
            self.planes[q].exact_coeffs(),
            self.planes[r].exact_coeffs(),
            self.planes[j].exact_coeffs(),
        ) {
            let tp = ThreePlane([cp, cq, cr]);
            return Some(indirect_plane_side(&tp, cj) * self.planes[j].frame_sign());
        }
        // One shared motion ⇒ the same question, exactly, on the canonicalised pre-motion data.
        //
        // ★★★ **All four planes are asked for coefficients, `j` included** — the branch above's
        // rule. Asking `j` by its pre-motion *triangle* while the other three speak in
        // coefficients describes `j` twice, and the two part whenever its `d` is a rounded product
        // — measured at 27% of the census's rotated classes and 40% of the fin sweep's, answering
        // 30% and 75% of this branch's judgements under the mismatch. One vocabulary removes the
        // mismatch instead of detecting it.
        //
        // ★ `frame_sign` is what carries the convention across: `indirect_plane_side` answers
        // about the plane's own normal, while this predicate's contract is the *triangle's*
        // right-hand normal, and `frame_sign` is exactly the relation between the two.
        if shared_motion(self.planes, &[p, q, r, j]) {
            if let (Some(cp), Some(cq), Some(cr), Some(cj)) = (
                self.planes[p].base_coeffs(),
                self.planes[q].base_coeffs(),
                self.planes[r].base_coeffs(),
                self.planes[j].base_coeffs(),
            ) {
                let tp = ThreePlane([cp, cq, cr]);
                return Some(indirect_plane_side(&tp, cj) * self.planes[j].frame_sign());
            }
        }
        // A wide name reaches here (no f64 spelling exists for it), and its integers are still
        // an exact description — the same question the arm above asks, in the integer twin,
        // with the same `frame_sign` bridge on `j` (the one direction-sensitive slot: a negated
        // `p`/`q`/`r` row negates `D` and the dot together, so only `j`'s orientation matters).
        if let Some([bp, bq, br, bj]) = self.name_rescue([p, q, r, j], true) {
            return Some(
                nacre_exact::int_plane_side([&bp, &bq, &br], &bj) * self.planes[j].frame_sign(),
            );
        }
        None
    }

    /// **The name-integer rescue gate**: the rows for a question every one of whose planes is
    /// named, **at least one wide**, and whose frames the name can speak for — either no plane
    /// is moved (the names are world descriptions) or all carry one motion (the names are one
    /// shared pre-motion description, parity-corrected below). `None` keeps every existing
    /// route exactly as it was — in particular, an all-narrow question never takes this gate,
    /// so a corpus with no wide plane is untouched to the bit.
    ///
    /// ★ **The parity correction is on the coefficients, and it is not "negate x".** The
    /// canonicalised base frame reflects the *points* in `x` (see [`Witness::base_tri`]), and a
    /// plane derived from reflected points is `det(C)·(C·n, d)` — the cross product is a
    /// pseudovector, so the reflection `C` (negate the x *coefficient*) arrives with one more
    /// global sign `det(C) = −1`. The two together negate `y`, `z` and `d` and keep `x`:
    /// exactly what makes these rows a positive multiple of what `base_coeffs` would hold if
    /// the base were `f64`-representable, so the arms above transplant verbatim.
    fn name_rescue<const N: usize>(
        &self,
        idx: [usize; N],
        allow_shared: bool,
    ) -> Option<[[BigInt; 4]; N]> {
        let planes = self.planes;
        let rows = idx.map(|k| planes[k].name_ints());
        if rows.iter().any(|r| r.is_none()) || !rows.iter().flatten().any(|n| n.wide) {
            return None;
        }
        if !any_rotated(planes, &idx) {
            return Some(rows.map(|r| r.expect("checked above").ints.clone()));
        }
        if allow_shared && shared_motion(planes, &idx) {
            let odd = chain_parity(&plane_def(planes, idx[0])[0].chain) < 0;
            return Some(rows.map(|r| {
                let v = &r.expect("checked above").ints;
                if odd {
                    [v[0].clone(), -&v[1], -&v[2], -&v[3]]
                } else {
                    v.clone()
                }
            }));
        }
        None
    }

    /// **The one body of the certified `orient3d`**, given the implicit point's Cramer parts.
    ///
    /// ★ Both entry points end here — [`Judge::orient3d`] builds the parts and throws them away,
    /// [`ImplicitPoint::orient3d`] keeps them. Extracting it is not tidiness: a second copy of this
    /// body is a second place for the escalation policy to live, free to drift into a different
    /// answer for the same question. The route selection has the same hazard and the same answer,
    /// [`Judge::orient3d_cheap`].
    fn orient3d_given(&self, cr: CramerParts, p: usize, q: usize, r: usize, j: usize) -> i8 {
        let (dp, dq, dr, dj) = (
            plane_def(self.planes, p),
            plane_def(self.planes, q),
            plane_def(self.planes, r),
            plane_def(self.planes, j),
        );
        self.record(
            Site::Orient3d { p, q, r, j },
            orient3d_from_cramer(
                cr,
                borrow3(dp),
                borrow3(dq),
                borrow3(dr),
                &dj[0],
                &dj[1],
                &dj[2],
                self.standard,
            ),
        )
    }

    /// The interval planes of `∩(p, q, r)`, from the per-class cache.
    fn cramer_of(&self, p: usize, q: usize, r: usize) -> CramerParts {
        cramer_iv([self.plane_iv(p), self.plane_iv(q), self.plane_iv(r)])
    }

    /// **The implicit point `∩(p, q, r)`, as something you can ask questions of.**
    ///
    /// ★ On the certified route a point's Cramer parts are four determinants that do not mention the
    /// query at all — `(D, Dvec)` *is* the point. Naming the point makes sharing them a property of
    /// the value, so a caller with several questions about one point needs no "pair" entry point at
    /// every layer (there were three, one per crate boundary, and they capped the sharing at exactly
    /// two questions).
    ///
    /// **Hoist the handle as far as the point is constant.** In the arrangement's crossing collector
    /// that is the *wall pair*, so every segment on one wall shares it — not just one segment's two
    /// endpoints.
    pub fn point(&self, p: usize, q: usize, r: usize) -> ImplicitPoint<'_, W> {
        ImplicitPoint {
            jd: self,
            p,
            q,
            r,
            cramer: std::cell::OnceCell::new(),
        }
    }

    /// The sign of `a[axis] − b[axis]` between the two implicit points `a = ∩(planes a…)` and
    /// `b = ∩(planes b…)` (`+1` = `a[axis] > b[axis]`). `!rotated` → the exact implicit
    /// `cmp_coord` on the planes' name rows ([`PlaneWitness::exact_coeffs`]); `rotated` → each
    /// triple's three planes as exact
    /// `WitnessPoint` → [`indirect_cmp_coord_judge`].
    pub fn cmp_coord(&self, a: [usize; 3], b: [usize; 3], axis: usize) -> i8 {
        let planes = self.planes;
        let tp = |t: [usize; 3]| {
            Some(ThreePlane([
                planes[t[0]].exact_coeffs()?,
                planes[t[1]].exact_coeffs()?,
                planes[t[2]].exact_coeffs()?,
            ]))
        };
        if let (Some(ta), Some(tb)) = (tp(a), tp(b)) {
            return indirect_cmp_coord(&ta, &tb, axis);
        }
        if let Some(s) = cancel_cmp_coord(planes, a, b, axis) {
            return s;
        }
        // The world-gate only (`allow_shared: false`): a coordinate comparison is not
        // motion-invariant, and the rotation-axis reasoning that makes some moved cases exact
        // is `cancel_cmp_coord`'s above. Orientation-invariant, so the σ fold is harmless here.
        if let Some([a0, a1, a2, b0, b1, b2]) =
            self.name_rescue([a[0], a[1], a[2], b[0], b[1], b[2]], false)
        {
            return nacre_exact::int_cmp_coord([&a0, &a1, &a2], [&b0, &b1, &b2], axis);
        }
        let da = a.map(|k| plane_def(planes, k));
        let db = b.map(|k| plane_def(planes, k));
        self.record(
            Site::CmpCoord { a, b, axis },
            indirect_cmp_coord_judge(borrow_triple(da), borrow_triple(db), axis, self.standard),
        )
    }

    /// `sign(det[n_p; n_a; n_b])` over the three planes' stored normals — how the line `p ∩ a`
    /// runs relative to plane `b`, matching `plane_pair_dir_sign`'s shape (`+1`/`-1`/`0`).
    ///
    /// `!rotated` → `det3_sign` of the name rows' normals ([`PlaneWitness::exact_normal`], folded
    /// to the stored orientation). `rotated` → the kernel `D`
    /// (det of the *outward* `tri` normals, [`dir_sign_judge`]) bridged to the *stored*-normal
    /// convention by the per-plane [`PlaneWitness::frame_sign`]: `det(stored) =
    /// frame_sign(p)·frame_sign(a)·frame_sign(b)·det(outward)`.
    pub fn plane_pair_dir_sign(&self, p: usize, a: usize, b: usize) -> i8 {
        let planes = self.planes;
        // ★ **Only the normals are read**, so the row's normal half is enough — a plane whose `d`
        // is too wide for `f64` keeps this route.
        if let (Some(np), Some(na), Some(nb)) = (
            planes[p].exact_normal(),
            planes[a].exact_normal(),
            planes[b].exact_normal(),
        ) {
            return det3_sign([np, na, nb]);
        }
        // A determinant of normals: a motion multiplies it by `det(R)`, which the canonicalised
        // base frame has already made `+1` (see `Witness::base_tri`), so one shared motion means
        // the pre-motion normals give the same sign, exactly.
        if shared_motion(planes, &[p, a, b]) {
            let row = |k: usize| planes[k].base_coeffs().map(|[x, y, z, _]| [x, y, z]);
            if let (Some(rp), Some(ra), Some(rb)) = (row(p), row(a), row(b)) {
                return det3_sign([rp, ra, rb]);
            }
        }
        // The integer twin of both arms above, for a wide name: the rows carry the stored
        // orientation already (the σ fold), so no `frame_sign` bridge appears — exactly as in
        // the `det3_sign` arms, and unlike the toleranced route below, which reads the outward
        // triangles and bridges back.
        if let Some([bp, ba, bb]) = self.name_rescue([p, a, b], true) {
            return nacre_exact::int_dir_sign([&bp, &ba, &bb]);
        }
        let (dp, da, db) = (
            plane_def(planes, p),
            plane_def(planes, a),
            plane_def(planes, b),
        );
        planes[p].frame_sign()
            * planes[a].frame_sign()
            * planes[b].frame_sign()
            * self.record(
                Site::DirSign { p, a, b },
                dir_sign_judge(borrow3(dp), borrow3(da), borrow3(db), self.standard),
            )
    }
}

/// One implicit point of a plane table — `∩(p, q, r)` — and the questions asked of it.
///
/// Made by [`Judge::point`]. The Cramer parts are filled by the first question that reaches the
/// certified route and reused by every later one.
///
/// ★ **`OnceCell`, not `Cell<Option<_>>`**: "filled once, never changes" is said by the type, the way
/// [`Witness::tri_pt3`] and `WitnessPoint`'s realization cell already say it. And a **cell rather than a
/// lock** because the handle is a local value that never crosses a thread — the shared thing is the
/// `Judge` behind it, borrowed immutably. (On the `Judge` this cache would need a lock, and a lock
/// here would serialize the parallel arrangement.)
///
/// One lifetime suffices: `Judge`'s fields are `&'a [W]`, `Standard`, `&'a Notes` and an owned
/// `Vec`, all covariant in `'a`, so `&'p Judge<'a, W>` coerces to `&'p Judge<'p, W>`.
pub struct ImplicitPoint<'a, W> {
    jd: &'a Judge<'a, W>,
    p: usize,
    q: usize,
    r: usize,
    cramer: std::cell::OnceCell<CramerParts>,
}

/// `(D, Dvec)`. ★ Kept out of every public signature — an interval in a public signature asks its
/// caller to reason about radii; see [`Judge`]'s `iv` field for why that stays inside.
type CramerParts = (nacre_exact::Bounded, [nacre_exact::Bounded; 3]);

impl<W: PlaneWitness> ImplicitPoint<'_, W> {
    /// Which side of plane `j` this point lies on — [`Judge::orient3d`] for the same four planes,
    /// sharing this point's Cramer parts with every other question asked of this handle.
    ///
    /// A cheap route answers without touching the cell, so a handle whose questions all take the
    /// exact or shared-motion path costs nothing over calling [`Judge::orient3d`] directly.
    pub fn orient3d(&self, j: usize) -> i8 {
        let (p, q, r) = (self.p, self.q, self.r);
        match self.jd.orient3d_cheap(p, q, r, j) {
            Some(o) => o,
            None => {
                let cr = *self.cramer.get_or_init(|| self.jd.cramer_of(p, q, r));
                self.jd.orient3d_given(cr, p, q, r, j)
            }
        }
    }
}

/// The predicate that runs **before** a plane table exists — it is what *defines* the classes,
/// so it asks only for a [`Witness`], never a plane's coefficients.
impl<W: Witness> Judge<'_, W> {
    /// Whether planes `i` and `j` are the **same plane**, decided on the faces' original
    /// coordinates instead of on their derived coefficients (three non-collinear points on a plane
    /// determine it, so "every point of `tri_j` lies on `tri_i`'s plane" is conclusive — but only
    /// under non-collinearity, so a degenerate `tri` answers `false`). `!rotated` → the exact
    /// `orient3d` on `tri`; any rotated → the exact `WitnessPoint` definitions and [`orient3d_judge`].
    pub fn planes_coplanar(&self, i: usize, j: usize) -> bool {
        let planes = self.planes;
        if tri_collinear(planes[i].tri()) || tri_collinear(planes[j].tri()) {
            return false;
        }
        if !any_rotated(planes, &[i, j]) {
            return planes[j]
                .tri()
                .iter()
                .all(|&q| plane_side_exact(planes[i].tri(), q) == 0);
        }
        // "Same plane" is a statement about incidence, which *any* motion preserves — proper or
        // not — so this one shortcut would survive an improper chain even uncorrected. It reads
        // the canonicalised base anyway, because one convention is cheaper to keep than two.
        if shared_motion(planes, &[i, j]) {
            if let (Some(ti), Some(tj)) = (planes[i].base_tri(), planes[j].base_tri()) {
                return tj.iter().all(|&q| plane_side_exact(ti, q) == 0);
            }
        }
        // ★ Compose the two motions and compare the exact planes. This is the only route that can
        // *prove* two differently-turned planes identical: everything below decides by narrowing
        // intervals, which can refute equality but never establish it, so it lands on
        // `Coincident` and records that it did.
        if let Some(same) = coplanar_by_composed_rotation(planes, i, j) {
            return same;
        }
        let (di, dj) = (plane_def(planes, i), plane_def(planes, j));
        // Three point-on-plane judgements, and the answer is their conjunction. **The note belongs to
        // the merge, not to the points**: what a reader needs to know is "these two faces became one
        // plane, on this evidence", and a definite sign anywhere means no merge happened and there is
        // nothing to report. So the loosest of the three is recorded, and only once all three agreed.
        let mut loosest: Option<Decision> = None;
        for q in dj.iter() {
            let d = orient3d_judge(q, &di[0], &di[1], &di[2], self.standard);
            if d.orient() != Orient::Zero {
                return false;
            }
            if looser(d, loosest) {
                loosest = Some(d);
            }
        }
        if let Some(d) = loosest {
            self.notes.push(Evidence {
                site: Site::PlanesCoplanar { i, j },
                outcome: d,
            });
        }
        true
    }
}

/// Whether any of the named planes is rotated — the per-predicate routing signal. A predicate
/// must escalate to the kernel if **any** — not all — of its planes is irrational: a single
/// rounded coordinate can flip an f64 `orient3d`/`cmp`, whereas all-rational planes are exact.
/// **Does plane `k`'s stored *normal* point the same way as its witness triangle's?**
///
/// The weaker of the two agreements, and the one a normals-only predicate needs. `[a, b, c]` is
/// the triangle's normal exactly when it is orthogonal to both edges — two dot products, in
/// `Expansion` so the test is exact rather than a rounding of one. The direction is then read from
/// their dot product, whose sign is safe in `f64`: parallel non-zero vectors cannot cancel.
///
/// ★ **Worth separating from [`coeff_exact`] because `d` is where the disagreement lives.**
/// Measured, the failures are all of the shape `raw·origin` rounding — `3.5 × 0.2` landing on
/// `0.7000000000000001` — which moves the plane without turning it. Demanding the stronger
/// agreement here cost 4.7x on the axis-aligned fold for nothing.
///
/// ★ **`cfg(test)`: this is the oracle, not the road.** A plane now *carries* whether its
/// coefficients and normal can be trusted ([`PlaneWitness::exact_coeffs`],
/// [`PlaneWitness::exact_normal`]), so production reads the answer instead of deriving it here.
/// What is left is checking that a producer's claim is true, which is a test's question.
#[cfg(test)]
pub(crate) fn coeff_normal_ok<W: PlaneWitness>(planes: &[W], k: usize) -> bool {
    use nacre_predicates::Expansion;
    let [ca, cb, cc, _] = planes[k].coeffs();
    let t = planes[k].tri().map(|p| p.as_array());
    let ortho = |q: [f64; 3]| {
        let e = [q[0] - t[0][0], q[1] - t[0][1], q[2] - t[0][2]];
        Expansion::two_product(ca, e[0])
            .add(&Expansion::two_product(cb, e[1]))
            .add(&Expansion::two_product(cc, e[2]))
            .sign()
            == 0
    };
    // ★ **Parallel is the whole condition — the direction is not part of it.** The stored normal
    // is allowed to *oppose* the triangle's, and `WorkingPlane::frame_sign` exists to record exactly
    // that; both branches of `plane_pair_dir_sign` already carry the convention (the exact one
    // takes the determinant of stored normals, the toleranced one multiplies the outward
    // determinant by the three `frame_sign`s). An earlier spelling here also demanded
    // `coeffs · cross(tri) > 0`, which would refuse every `frame_sign == -1` plane for a
    // disagreement it does not have. Measured: no such plane exists in any model in the suite, so
    // it was costing nothing — but a guard that is wrong for a reason nobody has hit yet is still
    // wrong, and the next model to carry one would lose its fast route silently.
    ortho(t[1]) && ortho(t[2])
}

/// Is plane `k`'s witness triangle exactly on its own stored coefficients? — the **full**
/// agreement, `d` included, which the predicates that build implicit points need.
///
/// ★ **`cfg(test)`, for the reason [`coeff_normal_ok`] states**: the plane carries this answer
/// now, and what remains here is the independent check of it.
#[cfg(test)]
pub(crate) fn coeff_exact<W: PlaneWitness>(planes: &[W], k: usize) -> bool {
    use nacre_predicates::Expansion;
    let [ca, cb, cc, cd] = planes[k].coeffs();
    planes[k].tri().iter().all(|q| {
        let [x, y, z] = q.as_array();
        Expansion::two_product(ca, x)
            .add(&Expansion::two_product(cb, y))
            .add(&Expansion::two_product(cc, z))
            .add(&Expansion::two_product(cd, 1.0))
            .sign()
            == 0
    })
}

pub fn any_rotated<W: Witness>(planes: &[W], idx: &[usize]) -> bool {
    idx.iter().any(|&k| planes[k].is_rotated())
}

/// The three exact [`WitnessPoint`] defining plane `k` — **borrowed**, never rebuilt.
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

pub fn plane_def<W: Witness>(planes: &[W], k: usize) -> &[WitnessPoint; 3] {
    planes[k].tri_pt3()
}

fn to_i8(o: Orient) -> i8 {
    match o {
        Orient::Positive => 1,
        Orient::Negative => -1,
        Orient::Zero => 0,
    }
}

/// Borrow an owned plane def as the `&WitnessPoint` tuple the judges take.
fn borrow3(d: &[WitnessPoint; 3]) -> (&WitnessPoint, &WitnessPoint, &WitnessPoint) {
    (&d[0], &d[1], &d[2])
}

/// Borrow three plane defs as the tuples `indirect_cmp_coord_judge` takes.
fn borrow_triple(d: [&[WitnessPoint; 3]; 3]) -> [(&WitnessPoint, &WitnessPoint, &WitnessPoint); 3] {
    [borrow3(d[0]), borrow3(d[1]), borrow3(d[2])]
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

/// **The one rotation a chain amounts to** — axis, total angle and the pivot they share — or
/// `None` when it is not of that shape.
///
/// The sibling of [`single_axis_motion`], and the difference is the **pivot**. That one runs on a
/// *difference* of two points, where a pivot and a translation both cancel; this one is for
/// planes, and a plane's `d` is a statement about position, so nothing cancels. Hence: one axis,
/// one pivot, no translation, no reflection — anything else is `None` and the caller escalates.
///
/// `None` in the axis slot means the chain does not turn at all, which composes with anything.
type OneRotation = (
    Option<nacre_exact::Axis>,
    nacre_exact::Angle,
    Option<[Rat; 3]>,
);

fn single_rotation(def: &[WitnessPoint; 3]) -> Option<OneRotation> {
    let mut axis: Option<nacre_exact::Axis> = None;
    let mut pivot: Option<[Rat; 3]> = None;
    let mut total = nacre_exact::Angle::from_deg(Rat::from_int(0))?;
    for n in def[0].chain.iter() {
        match n {
            MoveNode::Rotate {
                axis: a,
                angle,
                pivot: p,
            } => {
                total = total.checked_add(angle.deg())?;
                if *axis.get_or_insert(*a) != *a {
                    return None;
                }
                if *pivot.get_or_insert(*p) != *p {
                    return None;
                }
            }
            // A translation moves a plane, so unlike `single_axis_motion` it cannot be skipped;
            // a reflection is improper; and a frame (narrow or wide) is a rotation this shortcut
            // cannot *state* — it is not about a coordinate axis, so there is no axis and angle
            // to accumulate. All are conservative misses rather than wrong answers: the caller
            // escalates to the general path, which realizes the chain whatever its shape.
            MoveNode::Translate { .. }
            | MoveNode::Mirror { .. }
            | MoveNode::Frame { .. }
            | MoveNode::FrameWide(_)
            | MoveNode::FrameThrough(_) => {
                return None;
            }
        }
    }
    Some((axis, total, pivot))
}

/// **Are planes `i` and `j` the same plane, decided by composing their motions?** `None` when the
/// two chains do not compose into something the rationals can state, which is the caller's cue to
/// escalate.
///
/// Both planes state themselves exactly in their own pre-motion frame
/// (`base_coeffs_rat`), so the question `M_A(P_A) = M_B(P_B)` becomes
/// `M_B⁻¹M_A(P_A) = P_B`. When both chains turn about **the same axis through the same pivot**,
/// rotations commute and their angles add, so `M_B⁻¹M_A` is exactly `Rotate(axis, θ_A − θ_B)` —
/// and if that angle is one the rationals can state (the 90° family), carrying `P_A` across is
/// exact. Both sides are canonical, so `==` *is* plane identity.
///
/// ★ **A `false` here is a proof too**, not a failure to prove: the transported coefficients are
/// exact, so differing means the planes differ. That is what lets a non-coplanar pair skip the
/// escalation entirely.
///
/// ★★ **The preconditions are what make the composition valid, not an optimisation.** Rotations
/// about *different* axes or pivots compose into a motion whose translation part contains
/// `R_B⁻¹p` — irrational — so dropping either check would answer confidently and wrongly.
fn coplanar_by_composed_rotation<W: Witness>(planes: &[W], i: usize, j: usize) -> Option<bool> {
    let (ca, cb) = (planes[i].base_coeffs_rat()?, planes[j].base_coeffs_rat()?);
    let (axis_a, theta_a, pivot_a) = single_rotation(plane_def(planes, i))?;
    let (axis_b, theta_b, pivot_b) = single_rotation(plane_def(planes, j))?;
    let axis = match (axis_a, axis_b) {
        (None, None) => return Some(ca == cb), // neither turns: the frames already coincide
        (Some(a), None) | (None, Some(a)) => a,
        (Some(a), Some(b)) if a == b => a,
        _ => return None,
    };
    let pivot = match (pivot_a, pivot_b) {
        (Some(p), Some(q)) if p != q => return None,
        (Some(p), _) | (_, Some(p)) => p,
        (None, None) => [Rat::from_int(0); 3],
    };
    let delta = theta_a.checked_add(Rat::from_int(0).checked_sub(theta_b.deg())?)?;
    let iso = nacre_exact::Isometry::rotation(nacre_exact::Rotation {
        axis,
        pivot,
        angle: delta,
    });
    Some(iso.plane_coeffs(ca)? == cb)
}

/// The single rotation axis and total angle of a chain, or `None` if it turns about more than
/// one axis (then the product is not a rotation about a coordinate axis and this shortcut does
/// not apply — a conservative miss).
///
/// **Translations are skipped, and pivots are ignored, for the same reason**: this runs on the
/// *difference* of two points that share one chain, and a translation is the identity on a
/// difference — `(p + t) − (q + t) = p − q`. A rotation about any pivot likewise acts on a
/// difference as the pure linear `R`, which is why only the axis and angle are read.
///
/// **A reflection ends it.** Unlike the determinant shortcuts, which a canonicalised base frame
/// repairs, this one asks about a *particular axis*, and a reflection negates the axis it fixes —
/// there is no handedness correction that answers an axis-wise question. Two reflections do not
/// rescue it either: their product is a half-turn, not a turn by the angles read here. So the
/// declaration is that this shortcut is defined for **axis-preserving** motions, and anything
/// else escalates. A conservative miss, never a wrong sign.
fn single_axis_motion(def: &[WitnessPoint; 3]) -> Option<(nacre_exact::Axis, nacre_exact::Angle)> {
    let chain = &def[0].chain;
    let mut axis: Option<nacre_exact::Axis> = None;
    let mut total = nacre_exact::Angle::from_deg(nacre_exact::Rat::from_int(0))?;
    for n in chain.iter() {
        let a = match n {
            MoveNode::Rotate { axis, angle, .. } => {
                total = total.checked_add(angle.deg())?;
                axis
            }
            MoveNode::Translate { .. } => continue, // cancels in the difference
            MoveNode::Mirror { .. } => return None,
            // A frame *is* a rotation, and it does act on a difference as its pure linear part —
            // but this shortcut's answer is "which coordinate axis is preserved", and a frame
            // (narrow or wide) preserves none. Escalating is the honest miss.
            MoveNode::Frame { .. } | MoveNode::FrameWide(_) | MoveNode::FrameThrough(_) => {
                return None; // a judged frame is a rotation this shortcut cannot state, like Frame
            }
        };
        if *axis.get_or_insert(*a) != *a {
            return None;
        }
    }
    Some((axis?, total))
}

/// Whether the three points are **exactly collinear**, decided by the three coordinate-plane
/// projections of the cross product (each an exact `orient2d`). Non-collinearity is the
/// standing precondition of [`Judge::planes_coplanar`].
fn tri_collinear(t: [Point3; 3]) -> bool {
    let [a, b, c] = t.map(|p| p.as_array());
    let proj = |i: usize, j: usize| orient2d([a[i], a[j]], [b[i], b[j]], [c[i], c[j]]) == 0.0;
    proj(0, 1) && proj(1, 2) && proj(2, 0)
}

/// The exact side of triangle `tri`'s plane that the explicit point `p` lies on
/// (`+1`/`-1`/`0`), sharing [`Judge::orient3d`]'s convention (`p` takes `V`'s slot).
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

/// Is `d` weaker evidence than `best` — the one a report should quote?
///
/// The order is by how much is left unsaid: a proved sign says everything, a coincidence names a
/// distance, and the two inconclusive outcomes say the least. Among coincidences the wider bound
/// wins, since that is the one closest to being wrong.
fn looser(d: Decision, best: Option<Decision>) -> bool {
    let rank = |x: Decision| match x {
        Decision::Sign(_) => 0,
        Decision::Coincident { .. } => 1,
        Decision::Exhausted { .. } => 2,
        Decision::Degenerate => 3,
    };
    match best {
        None => rank(d) > 0,
        Some(b) if rank(d) != rank(b) => rank(d) > rank(b),
        Some(Decision::Coincident { within: y }) => match d {
            Decision::Coincident { within: x } => y.lt(x),
            _ => false,
        },
        Some(_) => false,
    }
}

#[cfg(test)]
#[path = "tests/predicate.rs"]
mod tests;
