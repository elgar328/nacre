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
use nacre_math::Point3;
use nacre_predicates::{
    ThreePlane, det3_sign, indirect_cmp_coord, indirect_orient3d, indirect_plane_side, orient2d,
    orient3d,
};
use nacre_scalar::{Orient, Rat};
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
    /// borrows it. It used to be an `Option` whose emptiness *also* meant "not rotated" — one
    /// field answering two questions, which is why it could not be filled in advance without
    /// changing which predicate path runs. [`Witness::is_rotated`] is now that second question.
    fn tri_pt3(&self) -> &[WitnessPoint; 3];
    /// Whether this plane came from a rotated solid — the predicate-routing signal, and a
    /// **different fact** from "does the definition carry a rotation chain". Neither
    /// `chain.is_empty()` nor `tol == 0` is equivalent to it (a rotated solid's face may witness
    /// its plane through points with no chain, and a 90°-family rotation has tol exactly 0), so
    /// the answer is carried, not derived.
    fn is_rotated(&self) -> bool;

    /// The plane's **exact rational coefficients in the frame its provenance names** — the world
    /// when unmoved, the pre-motion frame when moved. `None` when the producer recorded none.
    fn base_coeffs_rat(&self) -> Option<[nacre_scalar::Rat; 4]> {
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
    /// The plane's (un-normalized) coefficients `[a, b, c, d]` (`n·x + d = 0`) — **raw**.
    ///
    /// ★★★ **Raw means "not known to describe the same plane as [`tri`](Witness::tri)".** The two
    /// are both exact descriptions and they need not agree: `d` is an `f64` product, so a face at
    /// `y = −0.2` gets a coefficient plane `2⁻⁵⁴` from the one its own witness spans. A predicate
    /// that answers one question from here and the next from `tri` is describing two planes, and
    /// answers composed across them are not even an order — which is a defect this kernel has
    /// already had.
    ///
    /// **So a predicate reads [`exact_coeffs`](Self::exact_coeffs) or
    /// [`exact_normal`](Self::exact_normal) instead**, which are `None` exactly when the
    /// descriptions part. This one is for a caller that means to *relate* the two rather than
    /// choose between them.
    fn coeffs(&self) -> [f64; 4];

    /// The coefficients **only when they describe the same plane as [`tri`](Witness::tri)** —
    /// `None` when they do not, and for a rotated plane, which has no exact `f64` coefficients.
    ///
    /// A predicate that reads `d` needs this one: `d` is where the two descriptions part.
    fn exact_coeffs(&self) -> Option<[f64; 4]>;

    /// The **normal** under the weaker agreement — parallel to what `tri` spans, direction not
    /// required (that is what `frame_sign` records).
    ///
    /// ★ A predicate that never reads `d` may use this and keep its exact route on a plane
    /// [`exact_coeffs`](Self::exact_coeffs) has to refuse. Measured: requiring the full agreement
    /// for those cost 4.7x on the axis-aligned fold and bought nothing.
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
    /// (or the fold declined). This is what gives a *wide* name its exact shortcuts back:
    /// [`Witness::base_coeffs_rat`] and
    /// [`Self::base_coeffs`]/[`Self::exact_coeffs`] all read `PlaneName::narrow()` or `f64`, so
    /// a wide name answers `None` to every one of them and the judgement climbs — **1.9×** the
    /// escalations (246 against 475; the confound
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
/// and the stored plane may hold it either way round — measured, 1,916 interned statements in
/// the suite arrive `flipped`. So the raw name cannot be handed to a direction-sensitive
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
    name: Option<&nacre_scalar::PlaneName>,
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
    Some(NameInts {
        ints,
        wide: name.narrow().is_none(),
    })
}

/// **One operation's judging**: the witnesses it reasons over, the standard it holds them to, and
/// where the evidence goes.
///
/// These three are properties of the *operation*, not of a plane, and this is what says so. They
/// used to be stamped on every table row — which meant a two-phase construction (build the rows,
/// then stamp them), a placeholder for the gap between, and a guard for forgetting; a whole class
/// of mistake that exists only when a fact is stored somewhere it does not belong.
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
    /// with it. The type itself is public vocabulary (`nacre_scalar::Bounded`, beside `Mag` and
    /// `Rat`); the obligation to mint one correctly stays behind this crate's own constructors.
    /// `frame_sign` and `coeffs` are `i8` and `[f64; 4]`, so they *do* live on the witness; this
    /// one cannot follow them.
    ///
    /// Two workers racing to fill one cell compute the same value, so the answer does not depend on
    /// who won — the same argument `HpCell` rests on.
    iv: Vec<BoundedCell>,
    /// ★★★ **Whether each plane's stored coefficients and its witness triangle describe the same
    /// plane** — the condition under which the exact route may be taken.
    ///
    /// The two are *both* exact descriptions, and they need not agree: `coefficients()` is exact
    /// integer arithmetic only when the defining vertices are integers, and a face inherited
    /// through a boolean has `Discovered` vertices that are not. Measured, the witness of such a
    /// plane sits `2⁻⁵⁴` off its own coefficients. Routing by the *question* then describes one
    /// plane two ways, and answers composed across the two are not even an order.
    ///
    /// Same cell shape and same race argument as [`Judge::iv`].
    coeff_ok: Vec<OkCell>,
    /// The weaker agreement — normals only. See [`coeff_normal_ok`].
    normal_ok: Vec<OkCell>,
}

/// Lazily-filled cell for one plane's interval coefficients — `OnceLock` under `parallel` because
/// the boolean hands every worker the same `&Judge`, `OnceCell` otherwise.
#[cfg(feature = "parallel")]
type BoundedCell = std::sync::OnceLock<[nacre_scalar::Bounded; 4]>;
#[cfg(not(feature = "parallel"))]
type BoundedCell = std::cell::OnceCell<[nacre_scalar::Bounded; 4]>;

#[cfg(feature = "parallel")]
type OkCell = std::sync::OnceLock<bool>;
#[cfg(not(feature = "parallel"))]
type OkCell = std::cell::OnceCell<bool>;

impl<'a, W> Judge<'a, W> {
    pub fn new(planes: &'a [W], standard: Standard, notes: &'a Notes) -> Judge<'a, W> {
        Judge {
            planes,
            standard,
            notes,
            iv: (0..planes.len()).map(|_| BoundedCell::new()).collect(),
            coeff_ok: (0..planes.len()).map(|_| OkCell::new()).collect(),
            normal_ok: (0..planes.len()).map(|_| OkCell::new()).collect(),
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
    fn plane_iv(&self, k: usize) -> [nacre_scalar::Bounded; 4] {
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
    /// - `!rotated`: the exact path — the implicit-point `orient3d` (Attene) on the stored plane
    ///   coefficients and `tri` coordinates.
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

    /// Every route of [`Judge::orient3d`] **except** the certified one, or `None` when only that one
    /// is left.
    ///
    /// ★ **The routing lives here and nowhere else.** `Judge::orient3d_pair` needs to know whether
    /// two questions will both reach the certified path — that is the only branch with anything to
    /// share — and asking it by re-testing `any_rotated`/`shared_motion` would put the route
    /// selection in two places, free to drift into two different answers for one question.
    /// **May the exact route describe these planes?** — no rotation *and* every one of them
    /// coefficient-exact, so both routes would be talking about the same geometry.
    ///
    /// The rotation test alone is what this replaced, and it is not enough: it says the toleranced
    /// route is *needed*, not that the exact route is *equivalent*.
    pub fn exact_route_ok(&self, idx: &[usize]) -> bool {
        !any_rotated(self.planes, idx) && idx.iter().all(|&k| self.coeff_exact(k))
    }

    /// Plane `k`'s witness triangle sits exactly on its own stored coefficients — computed once.
    fn coeff_exact(&self, k: usize) -> bool {
        *self.coeff_ok[k].get_or_init(|| coeff_exact(self.planes, k))
    }

    /// **The route test for a predicate that reads only normals** — `d` cannot reach it, so the
    /// weaker agreement is what has to hold.
    pub fn exact_normal_route_ok(&self, idx: &[usize]) -> bool {
        !any_rotated(self.planes, idx)
            && idx
                .iter()
                .all(|&k| *self.normal_ok[k].get_or_init(|| coeff_normal_ok(self.planes, k)))
    }

    fn orient3d_cheap(&self, p: usize, q: usize, r: usize, j: usize) -> Option<i8> {
        if j == p || j == q || j == r {
            return Some(0);
        }
        // ★ **`j` is in the list although its coefficients are never read.** The query point comes
        // from `j`'s *triangle*, so "asked with `j`'s triangle" and "asked with `j`'s coefficients"
        // have to be the same question — which is exactly what `exact_coeffs` being `Some` says.
        // Dropping `j` here because "its coefficients are unused" would reopen the defect.
        if let (Some(cp), Some(cq), Some(cr), Some(_)) = (
            self.planes[p].exact_coeffs(),
            self.planes[q].exact_coeffs(),
            self.planes[r].exact_coeffs(),
            self.planes[j].exact_coeffs(),
        ) {
            let tp = ThreePlane([cp, cq, cr]);
            let tj = self.planes[j].tri();
            return Some(indirect_orient3d(
                &tp,
                tj[0].as_array(),
                tj[1].as_array(),
                tj[2].as_array(),
            ));
        }
        // One shared motion ⇒ the same question, exactly, on the canonicalised pre-motion data.
        //
        // ★★★ **All four planes are asked for coefficients, `j` included** — the same discipline
        // the branch above states, and for a stronger reason: here `j` has no triangle in the
        // question at all. This used to read `j`'s pre-motion *triangle* while the other three
        // spoke in coefficients, and the two descriptions of `j` part whenever its `d` was a
        // rounded product — measured at 27% of the census's rotated classes and 40% of the fin
        // sweep's, answering 30% and 75% of this branch's judgements under the mismatch. Asking
        // `j` in the same vocabulary as the rest removes the mismatch instead of detecting it.
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
                nacre_scalar::int_plane_side([&bp, &bq, &br], &bj) * self.planes[j].frame_sign(),
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
    /// `cmp_coord` on the stored coefficients; `rotated` → each triple's three planes as exact
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
            return nacre_scalar::int_cmp_coord([&a0, &a1, &a2], [&b0, &b1, &b2], axis);
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
    /// `!rotated` → `det3_sign` of the stored (un-normalized) normals. `rotated` → the kernel `D`
    /// (det of the *outward* `tri` normals, [`dir_sign_judge`]) bridged to the *stored*-normal
    /// convention by the per-plane [`PlaneWitness::frame_sign`]: `det(stored) =
    /// frame_sign(p)·frame_sign(a)·frame_sign(b)·det(outward)`.
    pub fn plane_pair_dir_sign(&self, p: usize, a: usize, b: usize) -> i8 {
        let planes = self.planes;
        // ★ **Only the normals are read below, so only they have to agree** — `d`'s rounding,
        // which is where the two descriptions actually part, never reaches this determinant.
        // Only the normals are read, so only they have to agree — `d`, where the two descriptions
        // actually part, never reaches this determinant.
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
            return nacre_scalar::int_dir_sign([&bp, &ba, &bb]);
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
type CramerParts = (nacre_scalar::Bounded, [nacre_scalar::Bounded; 3]);

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
pub fn coeff_normal_ok<W: PlaneWitness>(planes: &[W], k: usize) -> bool {
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
pub fn coeff_exact<W: PlaneWitness>(planes: &[W], k: usize) -> bool {
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
    Option<nacre_scalar::Axis>,
    nacre_scalar::Angle,
    Option<[Rat; 3]>,
);

fn single_rotation(def: &[WitnessPoint; 3]) -> Option<OneRotation> {
    let mut axis: Option<nacre_scalar::Axis> = None;
    let mut pivot: Option<[Rat; 3]> = None;
    let mut total = nacre_scalar::Angle::from_deg(Rat::from_int(0))?;
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
    let iso = nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
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
fn single_axis_motion(
    def: &[WitnessPoint; 3],
) -> Option<(nacre_scalar::Axis, nacre_scalar::Angle)> {
    let chain = &def[0].chain;
    let mut axis: Option<nacre_scalar::Axis> = None;
    let mut total = nacre_scalar::Angle::from_deg(nacre_scalar::Rat::from_int(0))?;
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
mod tests {
    /// ★★★★ **Nothing carries plane coefficients through a frame, and that is the design, not an
    /// omission.**
    ///
    /// A plane that passes through a motion has two descriptions — coefficients and a witness
    /// triangle — and mixing them is what once produced intransitive comparisons and a confident
    /// wrong winding. The frame design sidesteps the whole question: a plane built in a frame
    /// states itself **in that frame**, `BaseFrame` keeps its triangle and its coefficients in the
    /// *same* frame, and every judgement downstream is a determinant question a rigid motion
    /// preserves. So no `d` transformation term is needed anywhere.
    ///
    /// The one place that does transport coefficients is
    /// [`coplanar_by_composed_rotation`](super::coplanar_by_composed_rotation), whose whole method
    /// is composing two chains into a single coordinate-axis rotation and applying it to a
    /// coefficient vector. A frame is not that, so it must decline — a conservative miss that
    /// sends the pair to the interval path. This asks for that decline directly, because the
    /// alternative is `Isometry::plane_coeffs` being handed a frame it cannot express and
    /// answering anyway.
    #[test]
    fn a_chain_holding_a_frame_never_reaches_the_coefficient_transport() {
        use super::{WitnessPoint, single_rotation};
        use nacre_scalar::{Angle, Axis, Rat};
        let ri = |n: i128, d: i128| Rat::new(n, d).unwrap();
        // A tilted plane's frame: 2x - 3y + 7z + 11 = 0.
        let fr = nacre_scalar::plane_frame([2, -3, 7, 11].map(Rat::from_int)).unwrap();
        let framed = |u: i128, v: i128, w: i128| {
            WitnessPoint::at([ri(u, 1), ri(v, 1), ri(w, 1)])
                .frame(fr)
                .unwrap()
        };
        let def = [framed(0, 0, 0), framed(1, 0, 0), framed(0, 1, 0)];
        assert!(
            single_rotation(&def).is_none(),
            "a frame is not a rotation about a coordinate axis; transporting coefficients \
             through it would answer with full confidence and be wrong"
        );
        // The same points with an ordinary rotation *are* accepted, so the decline above is the
        // frame's doing and not a shortcut that never fires.
        let turned = [
            WitnessPoint::at([ri(0, 1), ri(0, 1), ri(0, 1)]),
            WitnessPoint::at([ri(1, 1), ri(0, 1), ri(0, 1)]),
            WitnessPoint::at([ri(0, 1), ri(1, 1), ri(0, 1)]),
        ]
        .map(|p| {
            p.rotate_about(
                Axis::Z,
                Angle::from_deg(Rat::from_int(30)).unwrap(),
                [Rat::from_int(0); 3],
            )
        });
        assert!(single_rotation(&turned).is_some());
    }

    use super::*;
    use nacre_scalar::Mag;

    /// How these fixtures judge; production chooses both per model. The coincidence limit is the
    /// derived default for a unit-scale model — output resolution (`2⁻⁵²`) two words further down.
    fn fixture() -> Standard {
        Standard {
            prec: 256,
            coincidence: Mag::pow2(-180),
            scale: Mag::of(1.0),
            cap: 4096,
        }
    }

    /// A synthetic axis-aligned plane witness: three points on the plane, its exact
    /// coefficients, and the exact `WitnessPoint` definition of those points. `is_rotated` is `false`,
    /// so predicates take the exact path — the definition is there but unused, which is
    /// precisely the arrangement the production tables now have.
    struct W {
        tri: [Point3; 3],
        coeffs: [f64; 4],
        def: [WitnessPoint; 3],
    }
    impl W {
        fn new(tri: [Point3; 3], coeffs: [f64; 4]) -> W {
            let def = tri.map(|p| {
                WitnessPoint::at_nearest(
                    p.as_array()
                        .map(|x| nacre_scalar::Rat::try_from_f64(x).expect("test coordinate")),
                )
            });
            W { tri, coeffs, def }
        }
    }
    impl Witness for W {
        fn tri(&self) -> [Point3; 3] {
            self.tri
        }
        fn tri_pt3(&self) -> &[WitnessPoint; 3] {
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
    }
    impl PlaneWitness for W {
        fn coeffs(&self) -> [f64; 4] {
            self.coeffs
        }
        fn frame_sign(&self) -> i8 {
            let t = self.tri();
            let x = (t[1] - t[0]).cross(t[2] - t[0]).as_array();
            let c = self.coeffs;
            if c[0] * x[0] + c[1] * x[1] + c[2] * x[2] > 0.0 {
                1
            } else {
                -1
            }
        }
        // The same rule the arrangement applies at construction — one implementation, so a test
        // witness routes exactly as the real one would.
        fn exact_coeffs(&self) -> Option<[f64; 4]> {
            nacre_predicates::plane_spanned_by(self.coeffs, self.tri.map(|p| p.as_array()))
                .then_some(self.coeffs)
        }
        fn exact_normal(&self) -> Option<[f64; 3]> {
            nacre_predicates::plane_normal_spanned_by(self.coeffs, self.tri.map(|p| p.as_array()))
                .then(|| [self.coeffs[0], self.coeffs[1], self.coeffs[2]])
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
        def: [WitnessPoint; 3],
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
            let def: [WitnessPoint; 3] = pts.map(|q| {
                let mut p = WitnessPoint::at([ri(q[0], 1), ri(q[1], 1), ri(q[2], 1)]);
                for &(axis, deg) in nodes {
                    p = p.rotate_about(axis, Angle::from_deg(ri(deg, 1)).unwrap(), pivot);
                }
                p
            });
            let base = pts.map(|q| Point3::from_array([q[0] as f64, q[1] as f64, q[2] as f64]));
            let tri: [Point3; 3] = std::array::from_fn(|i| Point3::from_array(def[i].coord()));
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
        fn tri_pt3(&self) -> &[WitnessPoint; 3] {
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
    }
    impl PlaneWitness for RW {
        fn coeffs(&self) -> [f64; 4] {
            self.coeffs
        }
        fn frame_sign(&self) -> i8 {
            let t = self.tri();
            let x = (t[1] - t[0]).cross(t[2] - t[0]).as_array();
            let c = self.coeffs;
            if c[0] * x[0] + c[1] * x[1] + c[2] * x[2] > 0.0 {
                1
            } else {
                -1
            }
        }
        // Rotated witnesses: no exact `f64` coefficients exist, as in the arrangement.
        fn exact_coeffs(&self) -> Option<[f64; 4]> {
            None
        }
        fn exact_normal(&self) -> Option<[f64; 3]> {
            None
        }
        fn base_coeffs(&self) -> Option<[f64; 4]> {
            Some(self.base_coeffs)
        }
    }

    /// A judging context over a fixture table.
    ///
    /// The collector is leaked so a fixture stays a one-liner — a `Vec` per call, in a test binary,
    /// and nothing reads it. A fixture that *does* want the evidence builds its own [`Notes`] and
    /// calls [`Judge::new`].
    fn jd<W>(planes: &[W]) -> Judge<'_, W> {
        Judge::new(planes, fixture(), Box::leak(Box::new(Notes::new())))
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
                        let got = jd(&ps).cmp_coord(a, b, k);
                        let want = to_i8(
                            indirect_cmp_coord_judge(
                                borrow_triple(a.map(|i| plane_def(&ps, i))),
                                borrow_triple(b.map(|i| plane_def(&ps, i))),
                                k,
                                fixture(),
                            )
                            .orient(),
                        );
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
            let want = to_i8(
                indirect_cmp_coord_judge(
                    borrow_triple(a.map(|i| plane_def(&ps, i))),
                    borrow_triple(b.map(|i| plane_def(&ps, i))),
                    k,
                    fixture(),
                )
                .orient(),
            );
            assert_eq!(jd(&ps).cmp_coord(a, b, k), want);
        }
    }

    /// The corner `∩(x=1, y=1, z=1) = (1,1,1)` sits above the `z=0` plane, so its orient3d
    /// against `z=0` is definite (nonzero), and querying `z=1` (a defining plane) is exactly 0.
    #[test]
    fn t_orient3d_axis_definite_and_on_plane() {
        let ps = cube_corner_planes();
        // query plane j = 3 (z=0): definite.
        assert_ne!(jd(&ps).orient3d(0, 1, 2, 3), 0);
        // query plane j = 2 (z=1) is one of the defining planes → exactly 0.
        assert_eq!(jd(&ps).orient3d(0, 1, 2, 2), 0);
    }

    /// A plane is coplanar with itself; two distinct planes are not.
    #[test]
    fn t_planes_coplanar_reflexive_and_distinct() {
        let ps = cube_corner_planes();
        assert!(
            jd(&ps).planes_coplanar(0, 0),
            "a plane is coplanar with itself"
        );
        assert!(
            !jd(&ps).planes_coplanar(0, 1),
            "x=1 and y=1 are distinct planes"
        );
    }

    /// **The two exact descriptions of one plane must agree, or the exact route must not be taken.**
    ///
    /// A plane carries stored coefficients `[a, b, c, d]` *and* a witness triangle, and both are
    /// exact — of different planes. `d` is `−(raw·origin)`, an `f64` product: for a face at
    /// `y = −0.2` with `raw = [0, −3.5, 0]` it lands on `0.7000000000000001`, which is a plane
    /// `2⁻⁵⁴` away from the one the triangle spans.
    ///
    /// That is tolerable as long as one plane is never described *both* ways. It was not: the
    /// route was chosen by the question — "does any plane here rotate?" — so a plane appeared at
    /// one position in one comparison and another in the next, and answers composed across the two
    /// were **not transitive**. `A == B`, `B < C`, `A > C` is what came out, and a lexicographic
    /// scan over that lands on a node that is not extreme.
    ///
    /// So the invariant is: **`coeff_exact` ⟹ the two describe one plane**, and only then may the
    /// exact route run. This pins the implication on a plane built to fail it.
    #[test]
    fn the_exact_route_is_refused_when_the_two_descriptions_disagree() {
        // A face at y = −0.2 spanned by an integer-ish triangle: `raw·origin` cannot be exact.
        let tri = [
            Point3::from_array([-4.0, -0.2, 0.0]),
            Point3::from_array([-0.5, -0.2, 0.0]),
            Point3::from_array([-0.5, -0.2, 1.0]),
        ];
        let coeffs = [0.0, -3.5, 0.0, -0.7000000000000001];
        // The witness is not on the coefficient plane: 3.5 × 0.2 is not 0.7000000000000001 / 1.
        let ps = [W::new(tri, coeffs)];
        assert!(
            !coeff_exact(&ps, 0),
            "the stored coefficients and the witness must be seen to disagree"
        );
        // ★ But their *directions* do agree — the rounding moved the plane, it did not turn it —
        // so a predicate that reads only normals keeps its fast route.
        assert!(coeff_normal_ok(&ps, 0));
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
        assert_eq!(ps[0].frame_sign(), 1);
        // reversing the tri winding flips cross(tri) → -1 (coeffs unchanged).
        let flipped = W::new([ps[0].tri[0], ps[0].tri[2], ps[0].tri[1]], ps[0].coeffs);
        assert_eq!(flipped.frame_sign(), -1);
    }
}
