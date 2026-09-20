use super::*;
/// Whether the point named by plane triple `query` lies inside a connected component — the exact
/// 3D lift of [`point_in_ring`]. A winding-parity ray whose line `L = a ∩ b` is built from two of
/// the query's own planes (never an arbitrary direction): each crossing with a face on plane `q`
/// is the exact three-plane point `{a,b,q}`, so the whole test is on the index-plane substrate —
/// no coordinate read, no f64. Non-convex is native (parity, not a convex test).
///
/// `faces` is the component as `(plane, rings)` per face — `rings[0]` = outer, `rings[1..]` = holes.
/// Only this component's faces are summed, so testing a void's vertex against a material component
/// reads `true` iff that material's outer shell nests the void.
///
/// A ray grazing a face boundary, or an undecidable in-face containment, is abandoned for the
/// query's next plane pair. **`Ok(None)` means every pair of this query's planes was blocked** —
/// this *node* cannot decide the question, which is a fact about the node, not an error: the
/// callers hold other nodes to try, and "every node abstained" is *their* proposition to reject
/// (`NoClearRay`, raised where the retries actually run out). `Err` is reserved for the
/// judgement itself failing (`JudgeExhausted`, a ring that cannot be named, …) — those must
/// propagate, never be traded for the next node: an abstention has other nodes as its remedy,
/// a failed judgement does not, and retrying it would let a real cause masquerade as
/// "no clear ray" once every node hit it.
/// One boundary of a component's face, with a polygon's edges already derived.
///
/// ★★★ **It mirrors [`crate::boolean::Bound`] on purpose.** The probe used to read a *flattened*
/// projection — `Vec<Vec<RingEdge>>`, which can only spell a polygon — so every circular and
/// banded boundary was **dropped on the way in** (`LocalFace::poly_rings`) and the component the
/// ray counted was not the component. Reading the engine's own boundary vocabulary is what lets
/// the ray count what is actually there.
#[derive(Clone, Debug)]
pub(crate) enum BoundEdges {
    /// A polygon, its edges carrying their walls (the only shape the ray counts today).
    Ring(Vec<RingEdge>),
    /// A whole circle, by the cylinder whose surface it rides — a disk face's outer bound, or a
    /// bored face's hole. Boxed for the same reason [`CompSurf::Cylinder`] is.
    Circle(Box<nacre_topo::CylinderDef>),
    /// A lateral face's **whole** boundary as loops on the cylinder's chart — a band's two rims,
    /// a panel's ring, a chain rim, the holes — outer and holes alike in one list.
    ///
    /// ★ One list and not «outer, then holes» because the chart is an annulus: a loop that
    /// wraps the cylinder has no inside, and a face between two chain rims is emitted as
    /// `Ring(outer)` + `inner = [the other chain]` by the region walk. What is true for every
    /// shape is the parity: a point of the cylinder is on the face iff the ray up the axis from
    /// it crosses the loops an odd number of times ([`loop_parity`]).
    Lateral(Vec<LateralLoop>),
}

/// One boundary loop of a lateral face on the cylinder's chart `(θ, z)`.
#[derive(Clone, Debug)]
pub(crate) enum LateralLoop {
    /// A whole-circle rim: the plane class it rides (⊥ to the axis).
    Circle(usize),
    /// A ring of arcs (on ⊥ classes) and rulings (on ∥ classes), its corners pierce names.
    Ring(Vec<RingEdge>),
}

/// A component face's surface, as the ray needs it.
///
/// ★ A cylinder carries its **truth**, not its class index: the crossings are solved against
/// `origin/dir/radius`, and the producer (`boolean`) is the one holding the class table. The
/// consumer should not have to look anything up — the same rule `RingEdge` states about walls.
#[derive(Clone, Debug)]
pub(crate) enum CompSurf {
    Plane(usize),
    /// Boxed because a `CylinderDef` is four rationals wide and every *planar* face would
    /// otherwise carry that much dead space — and planar faces are nearly all of them.
    Cylinder(Box<nacre_topo::CylinderDef>),
}

/// One face of a component — its surface, and its boundaries.
#[derive(Clone, Debug)]
pub(crate) struct CompFace {
    pub(crate) surf: CompSurf,
    pub(crate) outer: BoundEdges,
    pub(crate) inner: Vec<BoundEdges>,
}

/// How a ray met one cylindrical face.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CurvedHit {
    /// This many crossings lie on the **counted half** of the ray. ★ Which half that is belongs to
    /// the caller, not here: it hands in the plane, and this counts that plane's negative side.
    /// The two roads pick opposite halves — the named one counts behind its query, the coordinate
    /// one ahead of its origin — and the parity is the same either way, so a name that said
    /// "behind" would be false for one of them.
    Counted(usize),
    /// The ray touched a boundary — tangent, on a ruling, or exactly through the query. Abandon
    /// this ray and let the caller try the next plane pair, the same policy `point_on_ring` sets
    /// for a polygon.
    Graze,
}

/// **Is a face's material there**, given a way to ask one of its bounds — inside the outer bound
/// and outside every hole.
///
/// ★★★ **The combination is one rule; only the *asking* is two.** A named point answers a polygon
/// by a ring walk and a circle by its radial side; a rational one projects through a chart. Those
/// genuinely differ — the point arrives differently. But "outer and not any hole, and abandon the
/// moment a bound says *on me*" is the same sentence for both roads, and a sentence written twice
/// is one that drifts ([[rule-lives-inline-next-door]] is this repository's most-repeated defect).
/// It was written twice for one commit; this is that commit's correction.
///
/// ★★ **The hole clause fires and no test would notice if it stopped.** Measured: 7 entries with
/// a hole across `nacre-ops`, **6 of which subtract** — and stubbing the loop away leaves every
/// target green. That is a property of today's population, not of the rule: a hole here is a
/// **through** bore, so a ray that passes through it pierces the *pair* of annular caps and a
/// wrongly-counted crossing is wrongly counted **twice**, leaving the parity alone. A blind bore
/// on a classified component breaks the pairing and there is no such fixture yet. Extraction is
/// what guards this, not a lock.
pub(super) fn material_of(
    f: &CompFace,
    mut ask: impl FnMut(&BoundEdges) -> Result<Option<bool>, BoolError>,
) -> Result<Option<bool>, BoolError> {
    let Some(mut here) = ask(&f.outer)? else {
        return Ok(None);
    };
    if here {
        for hole in &f.inner {
            match ask(hole)? {
                Some(true) => {
                    here = false;
                    break;
                }
                Some(false) => {}
                None => return Ok(None),
            }
        }
    }
    Ok(Some(here))
}

/// **A component's probe** — a point to ask "how deep is this component nested" from.
///
/// ★ The two variants are two *descriptions of a point*, not two policies: a three-plane name is
/// exact without coordinates at all (so it survives a rotated class, where no rational coordinate
/// exists), and rational coordinates are what a component describes when **none of its faces
/// carries a vertex** — a lone cylinder's boundary is two disks and a band.
///
/// ★★ **Both are points *on* the component's boundary**, which is what makes their depths
/// comparable: today's named probe is a face vertex, and [`coord_probes`] takes a cap disk's
/// **centre**, which lies on that face. An *interior* witness — the axis midpoint, say — would
/// also answer containment, but it is a different kind of point and could be separated from the
/// boundary by another component's wall.
pub(crate) enum Probe {
    Named([usize; 3]),
    Coord {
        p: [nacre_exact::Rat; 3],
        dir: [nacre_exact::Rat; 3],
    },
}

/// **Is this probe inside the component?** — the one door in front of the two roads.
pub(crate) fn probe_in_component(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    probe: &Probe,
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    match probe {
        Probe::Named(x) => point_in_component(jd, cyls, *x, faces),
        Probe::Coord { p, dir } => point_in_faces_rat(jd, cyls, p, dir, faces),
    }
}
