use super::*;
use std::collections::BTreeMap;
/// The first topological defect in a boolean's output, checked **per solid** so one bad piece is
/// never masked by another (a sum-of-all Euler parity misses two odd-χ solids; a global count
/// misses much), or `None` if every result solid is a valid closed 2-manifold. In one pass per
/// solid it builds each edge's using faces (`edge_uses`), the vertex set and the `V/E/F/L/S` counts,
/// then, most-specific first:
///  0. an edge whose uses lie in **two shells** — a shell is one closed boundary, so one split
///     between two shells leaves each open though the solid's edge counts are whole and its Euler
///     count reads the extra shell as genus ([`RejectReason::OpenResultShell`]);
///  1. a **non-manifold vertex** (pinch) — [`nacre_topo::nonmanifold_vertices`]; sound for any
///     number of pinches (unlike Euler parity, which a second pinch flips back to even);
///  2. an edge stated on **one plane twice** — an interior boundary the coplanar merge had to
///     erase ([`RejectReason::CoplanarMerge`]);
///  3. an edge whose stated carriers are **not the surfaces its two faces lie on**
///     ([`RejectReason::EdgeCarrierMismatch`]);
///  4. an **odd Euler characteristic** `χ = V − E + F − L_i` (a valid closed solid's is even);
///  5. a **negative genus** `S − χ/2 < 0` (an even χ that still cannot be a solid).
///
/// A boolean output is freshly built, so its cells never alias another live solid's — the scoped
/// maps are exact.
pub(crate) fn check_result_topology(
    model: &Model,
    solids: &[Handle<Solid>],
) -> Option<(RejectReason, Option<crate::RejectWhere>)> {
    for &sh in solids {
        let solid = model.solid(sh);
        let mut shells: HashSet<Handle<Shell>> = HashSet::new();
        let mut faces: HashSet<Handle<Face>> = HashSet::new();
        let mut edge_uses: HashMap<Handle<Edge>, Vec<Handle<Face>>> = HashMap::new();
        let mut vertices: HashSet<Handle<Vertex>> = HashSet::new();
        let mut inner_loops: i64 = 0;
        let mut edge_shell: HashMap<Handle<Edge>, Handle<Shell>> = HashMap::new();
        let mut across: Vec<Handle<Edge>> = Vec::new();
        for &shell_h in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
            if !shells.insert(shell_h) {
                continue;
            }
            for &fh in &model.shell(shell_h).faces {
                if !faces.insert(fh) {
                    continue;
                }
                let face = model.face(fh);
                inner_loops += face.inner.len() as i64;
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        edge_uses.entry(he.edge).or_default().push(fh);
                        vertices.extend(model.edge(he.edge).vertices);
                        if *edge_shell.entry(he.edge).or_insert(shell_h) != shell_h {
                            across.push(he.edge);
                        }
                    }
                }
            }
        }
        // **Each shell closes on its own.** Witness: the smallest such edge's first vertex, for
        // replay determinism.
        if let Some(&eh) = across.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::OpenResultShell,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        // The returned list is sorted by handle index, so the first entry is a deterministic
        // witness (one of possibly several pinch vertices).
        if let Some(&vh) = nacre_topo::nonmanifold_vertices(model, faces.iter().copied()).first() {
            return Some((
                RejectReason::NonManifoldVertex,
                Some(crate::RejectWhere::Point(model.vertex_point(vh))),
            ));
        }
        // **No edge may say it separates one plane from itself.** An edge whose *stated*
        // surface pair is one plane twice is an interior boundary the coplanar cleaning was
        // obliged to erase. When the merge abstains instead — a 2-node cap's chord and its
        // arcs collide in the node-pair key, so the group ships as-is — the pair would sail
        // through without this, and `validate` names the spelling a producer bug
        // (`EdgeCarrierMismatch`: no edge separates a surface from itself). Measured on the
        // straddling flush boss (`rul flush` Fuse). An *arc*-carried boundary between two
        // coplanar faces states a cylinder in its pair and ships by design, so only a stated
        // plane-self pair is asked here — a cylinder-self pair is refused one door earlier, by
        // the edge push (`TwoCylinders`). The reject is the cleaning's own name: an abstention that
        // leaves this shape is a merge that was mandatory and did not happen. Witness: the
        // offending edge's first vertex, smallest edge handle for replay determinism.
        let mut same_plane: Vec<Handle<Edge>> = Vec::new();
        for &eh in edge_uses.keys() {
            let e = model.edge(eh);
            let [s0, s1] = e.surfaces;
            if s0 == s1 && matches!(model.surface(s0), nacre_topo::Surface::Plane { .. }) {
                same_plane.push(eh);
            }
        }
        if let Some(&eh) = same_plane.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::CoplanarMerge,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        // **And the other half of that proposition: the stated pair is the pair of surfaces
        // the two using faces lie on** — `validate`'s rule, read the same way: two faces on two
        // surfaces must *be* the stated pair; two faces on one surface (two laterals of one
        // cylinder sharing an arc) must have that surface among the two. The
        // edge's curve derives from the stated pair, so a pair the faces do not keep is a wrong
        // curve on a shell every count above calls closed, and only `validate` would say so — the
        // shape is a straight edge on a cylinder's ruling stated `(cylinder, wall)` between two
        // plane faces (`an_edge_stated_on_a_pair_its_faces_do_not_keep_is_refused` plants one).
        // An edge with other than two uses is the shell guard's to name.
        let mut disagreeing: Vec<Handle<Edge>> = Vec::new();
        for (&eh, uses) in &edge_uses {
            let [f0, f1] = uses[..] else {
                continue;
            };
            let stated = model.edge(eh).surfaces;
            let mut observed = [model.face(f0).surface, model.face(f1).surface];
            if observed[1].index() < observed[0].index() {
                observed.swap(0, 1);
            }
            let agrees = if observed[0] == observed[1] {
                stated.contains(&observed[0])
            } else {
                stated == observed
            };
            if !agrees {
                disagreeing.push(eh);
            }
        }
        if let Some(&eh) = disagreeing.iter().min() {
            let [va, _] = model.edge(eh).vertices;
            return Some((
                RejectReason::EdgeCarrierMismatch,
                Some(crate::RejectWhere::Point(model.vertex_point(va))),
            ));
        }
        let (v, e, f) = (
            vertices.len() as i64,
            edge_uses.len() as i64,
            faces.len() as i64,
        );
        let chi = v - e + f - inner_loops;
        if chi % 2 != 0 {
            return Some((RejectReason::EulerParity, None));
        }
        if shells.len() as i64 - chi / 2 < 0 {
            return Some((RejectReason::NegativeGenus, None));
        }
    }
    None
}

/// Connected components of the reconstructed faces, joined **where two of them bound one material
/// across a boundary curve**: round each curve, a face joins the next face round on its material
/// side — the two faces bounding one kept sector of the curve's star
/// ([`LineStar::bounding_pairs`](crate::draft::LineStar::bounding_pairs)). Returns a component label
/// per face (dense `0..n` in order of first appearance, for replay determinism) and the component
/// count.
///
/// **A curve two faces use joins them without reading anything** — two spokes have one cyclic
/// order, and each is the other's next face either way round. That is nearly every curve. Where
/// more faces use one, the curve's kind supplies the spokes' order and each face's material side,
/// and the one rule pairs them: [`circle_star`] for a whole rim or a cut arc, [`line_star`] for a
/// straight edge. A star that cannot
/// be read is refused by its cause, never left unjoined — an unjoined face is what turns two bodies
/// into a false `NonManifoldResultEdge`.
///
/// **Touching is not joining.** Two bodies that meet along a line or a circle share it — four
/// faces use it, two of each — and each face's next on its material side is its own body's, so
/// they land in two components and come back as the two bodies they are, not as one pinched
/// pseudo-solid. Where a body touches **itself** the material still runs around the contact, and
/// the component stays one: the guard fires exactly where no pair of solids exists.
///
/// This is the 3D reading of a rule the assembly stage already applies one dimension down:
/// [`unify_coplanar_faces`] splits a coplanar group "into **edge-connected components**, because
/// faces that merely lie on the same plane without touching must each survive on their own".
///
/// An enclosed void is its own component — its boundary shares no edge with the outer.
pub(super) fn face_components(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
) -> Result<(Vec<usize>, usize), BoolError> {
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    // The faces round each ring edge, with where each meets it. One use is a dangling edge — the
    // edge-use guard's to name. Both tables are walked in key order, so a refusal names the same
    // curve run after run.
    let mut users: BTreeMap<NodeEdge, Vec<EdgeUse<'_>>> = BTreeMap::new();
    for (face, lf) in faces.iter().enumerate() {
        for ring in rings_of(lf) {
            for t in 0..ring.len() {
                users
                    .entry(NodeEdge::of(ring, t))
                    .or_default()
                    .push(EdgeUse { face, ring, t });
            }
        }
    }
    // **A whole rim** has no node at all — a circle has none — so a cap face and the lateral that
    // meet on it are keyed by the rim, exactly as `reconstruct` keys it, and where each meets it is
    // read off the kind of bound.
    let mut rims: BTreeMap<(usize, usize), Vec<(usize, RimSpoke)>> = BTreeMap::new();
    for (i, lf) in faces.iter().enumerate() {
        let bounds = std::iter::once((&lf.outer, true)).chain(lf.inner.iter().map(|b| (b, false)));
        for (b, outer) in bounds {
            match (b, lf.surf) {
                (Bound::Circle { cyl }, ClassIx::Plane(c)) => {
                    let spoke = if outer {
                        RimSpoke::CapInside
                    } else {
                        RimSpoke::CapOutside
                    };
                    rims.entry((*cyl, c)).or_default().push((i, spoke));
                }
                (Bound::Rim { plane, ccw }, ClassIx::Cyl(k)) => {
                    let spoke = if *ccw {
                        RimSpoke::LateralAbove
                    } else {
                        RimSpoke::LateralBelow
                    };
                    rims.entry((k, *plane)).or_default().push((i, spoke));
                }
                _ => {}
            }
        }
    }
    let mut join = |a: usize, b: usize| {
        let (ra, rb) = (uf_find(&mut parent, a), uf_find(&mut parent, b));
        parent[ra] = rb;
    };
    for (key, us) in &users {
        match (us.len(), key) {
            (0 | 1, _) => {}
            (2, _) => join(us[0].face, us[1].face),
            (_, NodeEdge::Arc { cyl, plane, .. }) => {
                let spokes = us
                    .iter()
                    .map(|e| Ok((e.face, arc_spoke(jd, cyls, faces, e.face, e.wall(), *cyl)?)))
                    .collect::<Result<Vec<_>, BoolError>>()?;
                for (a, b) in circle_star(jd, cyls, faces, *cyl, *plane, &spokes)? {
                    join(a, b);
                }
            }
            (_, NodeEdge::Line(_)) => {
                for (a, b) in line_star(jd, cyls, faces, us)? {
                    join(a, b);
                }
            }
        }
    }
    for (&(k, c), us) in &rims {
        match us.len() {
            0 | 1 => {}
            2 => join(us[0].0, us[1].0),
            _ => {
                for (a, b) in circle_star(jd, cyls, faces, k, c, us)? {
                    join(a, b);
                }
            }
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let mut labels = Vec::with_capacity(faces.len());
    let mut n = 0;
    for i in 0..faces.len() {
        let root = uf_find(&mut parent, i);
        let next = label.len();
        let l = *label.entry(root).or_insert(next);
        labels.push(l);
        n = n.max(l + 1);
    }
    Ok((labels, n))
}

/// **Where a face meets a circle, round it** — the four directions, in the cross-section through
/// the axis, a face can leave a cylinder's circle on a plane along: the cap inside the circle, the
/// lateral above the plane, the cap outside it, the lateral below. In that cyclic order whatever
/// the plane's tilt, and these are all there are — a circle lies on its plane and its cylinder
/// alone (an arc's key names both).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RimSpoke {
    CapInside = 0,
    LateralAbove = 1,
    CapOutside = 2,
    LateralBelow = 3,
}

impl RimSpoke {
    /// The four, in their cyclic order.
    const ROUND: [RimSpoke; 4] = [
        RimSpoke::CapInside,
        RimSpoke::LateralAbove,
        RimSpoke::CapOutside,
        RimSpoke::LateralBelow,
    ];
}

/// **Where face `i` meets the cut arc of cylinder class `k` it walks along `wall`** — the arc's
/// travel against the face's frame says which side of it the face lies on, the way a whole circle's
/// bound kind does. Its ring keeps the face on its left about the frame's normal before `flip`
/// (`arrangement::emit`), so travelling counter-clockwise about the axis a lateral lies above the
/// arc (the convention `Bound::Rim` states), and a cap lies inside it exactly when its frame faces
/// up the axis (the even half-edge borders the disk exactly then — `emit_faces`' arc labels).
///
/// A lateral on the arc is one of cylinder `k`'s, the cylinder the arc's key names: an arc on two
/// cylinders' laterals is a contact the gate refuses (`CylinderPairContact`), so another reads as
/// the stages disagreeing.
fn arc_spoke(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    i: usize,
    wall: Wall,
    k: usize,
) -> Result<RimSpoke, BoolError> {
    let Wall::Arc { ccw, plane, .. } = wall else {
        return Err(reject(RejectReason::CylinderStagesDisagree));
    };
    Ok(match faces[i].surf {
        ClassIx::Cyl(c) if c != k => return Err(reject(RejectReason::CylinderStagesDisagree)),
        ClassIx::Cyl(_) if ccw => RimSpoke::LateralAbove,
        ClassIx::Cyl(_) => RimSpoke::LateralBelow,
        ClassIx::Plane(_) => {
            let up = crate::planes::frame_faces_up_the_axis(&jd.planes[plane], &cyls[k].def)
                .ok_or_else(|| reject(jd.planes[plane].world_absent()))?;
            if up == ccw {
                RimSpoke::CapInside
            } else {
                RimSpoke::CapOutside
            }
        }
    })
}

/// **The faces that bound one kept sector round the circle of cylinder class `k` on plane class
/// `c`** — `users` are the faces meeting it and where. Each face's material side is read off its
/// truth:
///
/// - a cap's material lies below the plane along the axis exactly when its outward normal points up
///   the axis — the class's frame facing up it ([`crate::planes::frame_faces_up_the_axis`]) and the
///   face not flipped, or the frame facing down and the face flipped;
/// - a lateral's lies inside the cylinder exactly when it is not flipped (its outward normal looks
///   away from the axis).
///
/// Two faces leaving along one direction have no order (`UnorderedEdges`); a class whose facing
/// cannot be read (by the cause of its missing world statement,
/// [`crate::planes::WorkingPlane::world_absent`], as `reconstruct` reads the same facing — the gate
/// keeps a class with no world statement off a cylinder's circles) and neighbours whose sides
/// disagree (`LabelConflict`) are refused.
fn circle_star(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    k: usize,
    c: usize,
    users: &[(usize, RimSpoke)],
) -> Result<Vec<(usize, usize)>, BoolError> {
    let mut at: [Option<usize>; 4] = [None; 4];
    for &(i, spoke) in users {
        if at[spoke as usize].replace(i).is_some() {
            return Err(reject(RejectReason::UnorderedEdges));
        }
    }
    // Whether the face's material lies the next spoke on in `RimSpoke` order.
    let forward = |i: usize, spoke: RimSpoke| -> Result<bool, BoolError> {
        let flip = faces[i].flip;
        Ok(match spoke {
            RimSpoke::CapInside | RimSpoke::CapOutside => {
                let up = crate::planes::frame_faces_up_the_axis(&jd.planes[c], &cyls[k].def)
                    .ok_or_else(|| reject(jd.planes[c].world_absent()))?;
                // Material above the plane is toward `LateralAbove`: on from inside, back from
                // outside.
                let above = up == flip;
                above == (spoke == RimSpoke::CapInside)
            }
            RimSpoke::LateralAbove | RimSpoke::LateralBelow => {
                // Material inside is toward `CapInside`: back from above, on from below.
                let inside = !flip;
                inside == (spoke == RimSpoke::LateralBelow)
            }
        })
    };
    let mut star = Vec::with_capacity(users.len());
    for (spoke, slot) in RimSpoke::ROUND.into_iter().zip(at) {
        if let Some(i) = slot {
            star.push((i, forward(i, spoke)?));
        }
    }
    crate::draft::LineStar::new(star)
        .bounding_pairs()
        .ok_or_else(|| reject(RejectReason::LabelConflict))
}

/// **One face's use of a ring edge** — the face, the ring it walks the edge in, and the edge's
/// index there.
#[derive(Clone, Copy)]
struct EdgeUse<'a> {
    face: usize,
    ring: &'a Ring,
    t: usize,
}

impl EdgeUse<'_> {
    fn wall(&self) -> Wall {
        self.ring.walls[self.t]
    }

    /// `+1` where this use walks the edge the way `other` does, `−1` the other way: a straight
    /// edge's uses all run between one pair of nodes.
    fn along(&self, other: &EdgeUse<'_>) -> i8 {
        if self.ring.nodes[self.t] == other.ring.nodes[other.t] {
            1
        } else {
            -1
        }
    }
}

/// **The faces that bound one kept sector round a straight edge** more than two faces use — the
/// spokes are half-planes and half-laterals about the line, ordered by [`ruling_spokes`] where the
/// line is a ruling of a cylinder some face states it on, by [`plane_spokes`] where only planes
/// hold it.
fn line_star(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    us: &[EdgeUse<'_>],
) -> Result<Vec<(usize, usize)>, BoolError> {
    let ruling = us.iter().find_map(|e| match e.wall() {
        Wall::Ruling { cyl, up, .. } => {
            let travel: i8 = if up { 1 } else { -1 };
            Some((cyl, travel * e.along(&us[0])))
        }
        Wall::Plane(_) | Wall::Arc { .. } => None,
    });
    let star = match ruling {
        Some((k, up)) => ruling_spokes(jd, cyls, faces, us, k, up)?,
        None => plane_spokes(jd, cyls, faces, us)?,
    };
    crate::draft::LineStar::new(star)
        .bounding_pairs()
        .ok_or_else(|| reject(RejectReason::LabelConflict))
}

/// **The spokes round a line only planes hold, in order** — each face with whether its material
/// lies counter-clockwise after it.
///
/// Cut the star by a plane class `t` the line `L = n_p × n_q` crosses, and each face's half-plane
/// is a ray on `t` from the line's point: along `n_t × n_a` for a face on class `a`, the way that
/// enters the face's half — so the half-planes' cyclic order about `L` is the rays' about `t`'s
/// outward normal ([`combinatorics::angular_order`], which also names two faces leaving one way).
/// Which `t` serves does not matter — the cut keeps the order, and `sign(n_t·L)` is in each ray's
/// sense — so the first class in the table that crosses the line is taken.
///
/// The signs, stored normals throughout (`EdgeDir`'s convention). A face on `a` lies left of its
/// travel `τ·L` about its frame normal `f_a·n_a` (its ring before `flip`, `arrangement::emit`), so
/// its half runs along `f_a·τ·(n_a × L)`. Against that, `(n_t × n_a)·(n_a × L) = −(n_t·L)·|n_a|²`,
/// so the ray's sense is `s = −sign(n_t·L)·f_a·τ`. Its material is on the side `−f_a·φ·n_a` of
/// `a` (`φ = −1` flipped), and the side counter-clockwise after the ray about `f_t·n_t` is
/// `−f_t·s`'s — so the material lies counter-clockwise after it exactly when `f_t·s = f_a·φ`.
/// `τ` is the first face's edge direction ([`combinatorics::dir_at`], `+1` along `n_p × n_q`) and,
/// for the others, how they walk the edge against it.
fn plane_spokes(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    us: &[EdgeUse<'_>],
) -> Result<Vec<(usize, bool)>, BoolError> {
    let stages = || reject(RejectReason::CylinderStagesDisagree);
    let first = &us[0];
    let ClassIx::Plane(p) = faces[first.face].surf else {
        return Err(stages());
    };
    let edges = first.ring.edges(jd, cyls, Some(p))?;
    let e = &edges[first.t];
    let combinatorics::EdgeDir::Line { carrier: q, sense } =
        combinatorics::dir_at(jd, cyls, p, e, e.node)?
    else {
        return Err(stages());
    };
    // A table with no class across the line has no cut to read the order on — the order is not
    // read, and the refusal says so. Measured: no straight star in the suite, the ignored sweep,
    // the census or the kit reaches it.
    let t = (0..jd.planes.len())
        .find(|&t| jd.plane_pair_dir_sign(t, p, q) != 0)
        .ok_or_else(|| reject(RejectReason::UnorderedEdges))?;
    let across = jd.plane_pair_dir_sign(t, p, q);
    let ft = jd.planes[t].frame_sign;
    let mut dirs = Vec::with_capacity(us.len());
    let mut material = Vec::with_capacity(us.len());
    for u in us {
        let ClassIx::Plane(a) = faces[u.face].surf else {
            return Err(stages());
        };
        let fa = jd.planes[a].frame_sign;
        let s = -across * fa * sense * u.along(first);
        let phi = if faces[u.face].flip { -1 } else { 1 };
        dirs.push(combinatorics::EdgeDir::Line {
            carrier: a,
            sense: s,
        });
        material.push(ft * s == fa * phi);
    }
    let order = combinatorics::angular_order(jd, t, &dirs)?;
    Ok(order
        .into_iter()
        .map(|j| (us[j].face, material[j]))
        .collect())
}

/// Where a face leaves a ruling, about the axis `m` — counter-clockwise from the outward radial
/// `r`: a tangent plane's half on the `+τ` side (`τ = m × r`), the lateral's half there (it bends
/// inward, past the plane), the secant planes' chords into the cylinder, the lateral's `−τ` half,
/// the tangent plane's, and the chords' outer halves.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum RulingSlot {
    TangentPlus,
    LateralPlus,
    Inward,
    LateralMinus,
    TangentMinus,
    Outward,
}

/// **The spokes round a ruling of cylinder class `k`, in order** — each face with whether its
/// material lies counter-clockwise after it about the axis `m`. `up` is `+1` where the first use
/// walks the edge up the axis.
///
/// - A lateral face lies left of its travel `σ·m` about its outward normal `r` (before `flip`):
///   `r × σm = −σ·τ` — up the axis it is the `−τ` half. Its material is inside the cylinder
///   exactly when it is not flipped, and inside lies counter-clockwise after the `+τ` half and
///   before the `−τ` one.
/// - A plane face on `a` leaves along `h = f_a·σ·(n_a × m)`, and its material is on the side
///   `−f_a·φ·n_a`; the side counter-clockwise after `h` is `m × h = f_a·σ·n_a`'s, so the material
///   lies there exactly when `σ = −φ` — up the axis exactly when flipped. Which way `h` points:
///   `h·r = f_a·σ·(n_a·τ)`, the chord into the cylinder where that is negative — `n_a·τ`'s sign is
///   [`crate::planes::plus_theta_is_above`] of the class's ruling side; on a tangent plane
///   (`side == 0`) `h·τ = −f_a·σ·(n_a·r)`, the `+τ` half where that is positive. Its normal is
///   square to the axis and the line lies on it, so `n_a·r = −(n_a·o + d_a)` for any point `o` of
///   the axis: the half is read off which side of the plane the axis lies on
///   ([`crate::planes::axis_above`]), and needs no point of the line.
/// - The chords are ordered by [`crate::planes::chord_before`], the chart's order at a station,
///   and their outer halves the same way. It reads the line's radial `r`
///   ([`crate::planes::SharedRuling::radial`]) off the gate's record of the two chords' classes:
///   two plane classes cutting the cylinder in one ruling of a lateral face are recorded as a pair.
///
/// A plane face states its side of the ruling on its own wall where it names the line as a ruling;
/// where it names it as the meet with another plane class, the line is the two classes' shared
/// ruling and the gate's record of it says ([`combinatorics::shared_line`]).
///
/// Every lateral and ruling wall here is cylinder `k`'s: a line on two cylinders' laterals is a
/// contact the gate refuses (`CylinderPairContact`), and a ruling's `side` means nothing against
/// another cylinder, so a stranger reads as the stages disagreeing.
fn ruling_spokes(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[LocalFace],
    us: &[EdgeUse<'_>],
    k: usize,
    up: i8,
) -> Result<Vec<(usize, bool)>, BoolError> {
    let stages = || reject(RejectReason::CylinderStagesDisagree);
    let def = &cyls[k].def;
    // Each plane face's class and ruling side; `None` for a lateral.
    let mut sides: Vec<Option<(usize, i8)>> = Vec::with_capacity(us.len());
    for u in us {
        if let Wall::Ruling { cyl, .. } = u.wall()
            && cyl != k
        {
            return Err(stages());
        }
        sides.push(match (faces[u.face].surf, u.wall()) {
            (ClassIx::Cyl(c), _) if c == k => None,
            (ClassIx::Cyl(_), _) => return Err(stages()),
            (ClassIx::Plane(a), Wall::Ruling { plane, side, .. }) if plane == a => Some((a, side)),
            (ClassIx::Plane(a), Wall::Plane(w)) => {
                let sr = cyls[k].shared_of(a, w).ok_or_else(stages)?;
                let line = combinatorics::shared_line(jd, cyls, k, sr, a)
                    .ok_or_else(|| reject(jd.planes[a].world_absent()))?;
                Some((a, line.side))
            }
            (ClassIx::Plane(_), _) => return Err(stages()),
        });
    }
    let mut placed: Vec<(RulingSlot, usize, usize, bool)> = Vec::with_capacity(us.len());
    for (j, u) in us.iter().enumerate() {
        let sigma = up * u.along(&us[0]);
        let flip = faces[u.face].flip;
        let (slot, class, material) = match sides[j] {
            None => {
                if sigma > 0 {
                    (RulingSlot::LateralMinus, usize::MAX, flip)
                } else {
                    (RulingSlot::LateralPlus, usize::MAX, !flip)
                }
            }
            Some((a, side)) => {
                let fa = jd.planes[a].frame_sign;
                let material = (sigma > 0) == flip;
                let slot = if side == 0 {
                    // `sign(n_a·r)`: the line lies on the far side of the plane from the axis.
                    let axis_above = crate::planes::axis_above(jd, a, def)
                        .ok_or_else(|| reject(jd.planes[a].world_absent()))?;
                    let nr: i8 = if axis_above { -1 } else { 1 };
                    if -fa * sigma * nr > 0 {
                        RulingSlot::TangentPlus
                    } else {
                        RulingSlot::TangentMinus
                    }
                } else {
                    let above =
                        crate::planes::plus_theta_is_above(jd, a, side).ok_or_else(stages)?;
                    if (fa * sigma > 0) != above {
                        RulingSlot::Inward
                    } else {
                        RulingSlot::Outward
                    }
                };
                (slot, a, material)
            }
        };
        placed.push((slot, class, u.face, material));
    }
    // In slot order, the chords within a slot by the chart's order — an insertion sort, so a
    // comparison that cannot answer is a refusal rather than a comparator that fails.
    let mut out: Vec<(RulingSlot, usize, usize, bool)> = Vec::with_capacity(placed.len());
    for x in placed {
        let mut at = out.len();
        for (i, y) in out.iter().enumerate() {
            let before = match x.0.cmp(&y.0) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Greater => false,
                std::cmp::Ordering::Equal => {
                    if x.1 == y.1 || !matches!(x.0, RulingSlot::Inward | RulingSlot::Outward) {
                        return Err(reject(RejectReason::UnorderedEdges));
                    }
                    let r = cyls[k]
                        .shared_of(x.1, y.1)
                        .ok_or_else(stages)?
                        .radial(def)
                        .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
                    crate::planes::chord_before(jd, def, &r, x.1, y.1).map_err(reject)?
                }
            };
            if before {
                at = i;
                break;
            }
        }
        out.insert(at, x);
    }
    Ok(out.into_iter().map(|(_, _, f, m)| (f, m)).collect())
}
