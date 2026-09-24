use super::*;
/// Dense plane ids for a face table: `(geom, plane_ix)` where `plane_ix[face]` indexes `geom`.
///
/// **The numbering is monotone in `canon`.** Roots are ranked in increasing order, so
/// `canon[i] < canon[j]` iff `plane_ix[i] < plane_ix[j]` — every comparison, sort and lex-min over
/// plane indices is order-isomorphic to the sparse form. Nothing found in the engine turns out to
/// depend on that (the two candidates — `loop_winding`'s lex-min node and `crossings`' pre-dedup
/// sort — are by coordinate and by set, respectively), but the audit cannot be proved exhaustive
/// over ~175 sites, so the numbering removes the question instead of answering it.
pub(crate) fn dense_planes(
    planes: &[FaceRow],
    canon: &[usize],
) -> (Vec<WorkingPlane>, Vec<ClassIx>, Vec<Handle<Surface>>) {
    // Plane classes densify from the union-find roots; cylinder rows never entered the
    // union-find (their identity question is the cylinder class table's), so here they
    // number by first-seen surface — the index space exists before the table does.
    let mut roots: Vec<usize> = canon
        .iter()
        .enumerate()
        .filter(|&(i, _)| matches!(planes[i], FaceRow::Plane(_)))
        .map(|(_, &c)| c)
        .collect();
    roots.sort_unstable();
    roots.dedup();
    let mut cyl_ix: HashMap<Handle<Surface>, usize> = HashMap::new();
    let plane_ix = canon
        .iter()
        .enumerate()
        .map(|(i, c)| match &planes[i] {
            FaceRow::Plane(_) => ClassIx::Plane(
                roots
                    .binary_search(c)
                    .expect("a class root is in the root set"),
            ),
            FaceRow::Cylinder(cf) => {
                let n = cyl_ix.len();
                ClassIx::Cyl(*cyl_ix.entry(cf.surf).or_insert(n))
            }
        })
        .collect();
    let geom = roots
        .iter()
        .map(|&r| {
            let pi = planes[r].plane();
            WorkingPlane {
                base_rat: pi.base_rat,
                world_rat: pi.world_rat,
                base: BaseFrame::of(&pi.tri_pt3, pi.motion, pi.orient_sign, pi.base_rat),
                name_ints: nacre_judge::predicate::name_stored_ints(
                    pi.name.as_ref(),
                    &pi.tri_pt3,
                    pi.orient_sign,
                ),
                plane: pi.plane,
                surf: pi.surf,
                tri: pi.tri,
                tri_pt3: pi.tri_pt3.clone(),
                rotated: pi.rotated,
                frame_sign: pi.orient_sign,
            }
        })
        .collect();
    // The cylinder surfaces in `ClassIx::Cyl` numbering order (first seen) — rebuilt from the
    // same map so the two cannot drift.
    let mut cyls: Vec<(usize, Handle<Surface>)> = cyl_ix.into_iter().map(|(s, k)| (k, s)).collect();
    cyls.sort_unstable_by_key(|&(k, _)| k);
    (geom, plane_ix, cyls.into_iter().map(|(_, s)| s).collect())
}

/// Whether three `WitnessPoint` are **exactly collinear** (zero-area triangle), decided on their
/// pre-rotation rational `base` coordinates. A rigid rotation preserves collinearity, and three
/// vertices of one solid share a rotation chain, so their bases are comparable; all three
/// coordinate-plane projections of `(b−a)×(c−a)` must vanish (exact `Rat`, no tolerance). An
/// i128 overflow returns `false` (treat as non-collinear): a genuinely-collinear triangle then
/// stays and is at worst rejected, never falsely skipped (which would drop a
/// real crossing — silent-wrong). Unrotated vertices carry `base == coord`, so this is the exact
/// zero-area (collinear) test on the vertices' rotation definitions.
#[cfg(test)]
pub(crate) fn pt3_base_collinear(a: &WitnessPoint, b: &WitnessPoint, c: &WitnessPoint) -> bool {
    use nacre_exact::Rat;
    let (a, b, c) = (&a.base, &b.base, &c.base);
    let proj_zero = |i: usize, j: usize| -> Option<bool> {
        let det = b[i]
            .checked_sub(a[i])?
            .checked_mul(c[j].checked_sub(a[j])?)?
            .checked_sub(
                b[j].checked_sub(a[j])?
                    .checked_mul(c[i].checked_sub(a[i])?)?,
            )?;
        Some(det == Rat::from_int(0))
    };
    matches!(
        (proj_zero(1, 2), proj_zero(2, 0), proj_zero(0, 1)),
        (Some(true), Some(true), Some(true))
    )
}

/// Each outer-shell edge with its bound vertices and the two combined-plane
/// indices of its adjacent faces, in first-seen (deterministic) order.
///
/// Every loop of every face is walked, holes included: a hole-ring edge is used
/// once by the holed face's inner loop and once by the neighbouring wall's outer
/// loop, so it too has exactly two incident faces. Walking `outer` before `inner`
/// on each face leaves the order of a hole-free solid untouched.
///
/// The pair is returned as `[usize; 2]`, so no caller can index a third slot: an
/// edge with any other incidence count is a non-manifold shell and rejects here.
#[allow(clippy::type_complexity)]
pub(crate) fn edge_incidence(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<Vec<(Handle<Edge>, [Handle<Vertex>; 2], [usize; 2])>, BoolError> {
    let mut order: Vec<Handle<Edge>> = Vec::new();
    let mut map: HashMap<Handle<Edge>, ([Handle<Vertex>; 2], Vec<usize>)> = HashMap::new();
    for sh in solid_shell_handles(model, solid) {
        for &fh in &model.shell(sh).faces {
            let face = model.face(fh);
            let pidx = surf_ix[&fh];
            for he in face_half_edges(face) {
                let bounds = model.edge(he.edge).vertices;
                let entry = map.entry(he.edge).or_insert_with(|| {
                    order.push(he.edge);
                    (bounds, Vec::new())
                });
                entry.1.push(pidx);
            }
        }
    }
    order
        .into_iter()
        .map(|e| {
            let (b, p) = map.remove(&e).unwrap();
            match p[..] {
                [x, y] => Ok((e, b, [x, y])),
                // `validate` would call this `NonOpposedEdge`, but `boolean` never runs
                // `validate` on its inputs, so the guard stays. No firing test.
                _ => Err(reject(RejectReason::NonManifoldEdge)),
            }
        })
        .collect()
}

/// Every half-edge of a face: its outer loop first, then each hole ring in order.
pub(crate) fn face_half_edges(face: &Face) -> impl Iterator<Item = &HalfEdge> {
    face.outer
        .half_edges
        .iter()
        .chain(face.inner.iter().flat_map(|l| l.half_edges.iter()))
}

/// Two faces lie on the same plane — by a **shared `Surface` handle** (explicit
/// sharing: O(1) `Handle` identity, exact, rotation-independent) or, as a fallback,
/// by the geometric rank-1 `planes_coplanar` test. A referenced coplanar contact —
/// a pad/pocket cap that reuses its face's surface — is caught by the handle path
/// without any coordinate test. On the axis-aligned M5 corpus the handle path is
/// redundant with `planes_coplanar` (same handle ⇒ same plane), so the geometric
/// fallback is what keeps independently-built coplanar contacts working; the handle
/// path's real payoff is rotated frames, where the geometric test would need the
/// rotation-exact judgment.
pub(crate) fn shares_or_coplanar(jd: &Judge<'_, FaceRow>, i: usize, j: usize) -> bool {
    let (pa, pb) = (jd.planes[i].plane(), jd.planes[j].plane());
    // Three independent witnesses, OR-ed, so this can only ever merge *more* than before:
    //  1. the same `Surface` handle — coplanar by reference (what an ops-built tool's base cap and
    //     its target face share, and what a chained operand's split coplanar faces share);
    //  2. exactly proportional coefficients — the original test, kept;
    //  3. the faces' own coordinates, exactly (`Judge::planes_coplanar`) — the only one of the three
    //     that does not read a *derived* value, and the one that catches two independently built
    //     solids whose walls coincide (`add_cuboid` stacked on `add_cuboid`), where the rounded
    //     coefficients of differently-sized faces are not exactly proportional.
    pa.surf == pb.surf || planes_coplanar(&pa.plane, &pb.plane) || jd.planes_coplanar(i, j)
}

/// Union-find root of `x` in `parent` (with path compression). Roots are the smallest index
/// of their class, so the result is deterministic — the same log replays to the same model.
/// Drives component grouping in [`crate::assembly::unify_coplanar_faces`].
pub(crate) fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut c = x;
    while parent[c] != r {
        let next = parent[c];
        parent[c] = r;
        c = next;
    }
    r
}

/// Canonicalize the combined plane table by coplanarity: two planes that are the same plane
/// (shared `Surface` handle, or exact rank-1 [`planes_coplanar`]) are merged into one class, so
/// a wall of `a` coplanar with a wall of `b` names a **single line** in a shared plane π. Without
/// it a wall and its coplanar twin order against each other as one (`order_along(R, R) == 0`);
/// canonicalizing turns that self-comparison into a real order. Returns `canon` where `canon[i]`
/// is the class root (the smallest index in the class). Every decision is exact
/// (`shares_or_coplanar`) — no coordinate. O(n²) scan over the (small) face count.
pub(crate) fn plane_classes(jd: &Judge<'_, FaceRow>) -> Vec<usize> {
    let planes = jd.planes;
    let n = planes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    // Cylinder rows stay self-rooted: their identity question belongs to the cylinder class
    // table, not the coplanarity union-find — they simply never enter `reps` below.

    // **Merge by `Surface` handle first, then compare only one face per distinct surface.**
    //
    // Coplanarity is an equivalence, and a shared handle *is* the same plane — exactly, by
    // identity, in O(1). So faces that share one are already one class, and asking the
    // geometric question of each of them separately asks the same question many times.
    //
    // The answer cannot move: union-find's result is the transitive closure of the pairs it was
    // given, and a class's root stays its minimum index, so dropping *redundant* pairs leaves
    // the partition and its roots alone. No tolerance is involved — a shared handle is identity.
    //
    // **Measured, and smaller than the pair count suggests.** On the 80-fin fold's largest
    // boolean, 406 faces carry 172 distinct surfaces: 82,215 pairs become 14,706, and over the
    // fold 7.52M become 2.30M — 3.3x fewer. But the scan only got ~1.4x faster (0.96s → 0.69s
    // over the fold, ~1.03x end to end), because **the pairs this drops are the cheapest ones**:
    // they matched on the handle and returned at the first `||`. What is left is the
    // geometrically distinct pairs, which are the ones that were expensive all along. Counting
    // removed operations overstates the saving whenever the removed ones are the cheap ones.
    let mut rep: HashMap<Handle<Surface>, usize> = HashMap::new();
    let mut reps: Vec<usize> = Vec::new();
    for (i, row) in planes.iter().enumerate() {
        let FaceRow::Plane(p) = row else { continue };
        match rep.get(&p.surf) {
            Some(&r) => {
                let (ri, rj) = (uf_find(&mut parent, r), uf_find(&mut parent, i));
                if ri != rj {
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
            None => {
                rep.insert(p.surf, i);
                reps.push(i);
            }
        }
    }

    for a in 0..reps.len() {
        for b in (a + 1)..reps.len() {
            let (i, j) = (reps[a], reps[b]);
            if shares_or_coplanar(jd, i, j) {
                let (ri, rj) = (uf_find(&mut parent, i), uf_find(&mut parent, j));
                if ri != rj {
                    // Attach the larger root under the smaller so a class's root is its min index.
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
    }
    (0..n).map(|i| uf_find(&mut parent, i)).collect()
}
