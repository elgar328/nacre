use super::*;
use crate::combinatorics::{EdgeFaces, edge_faces};
use nacre_judge::Standard;
use nacre_topo::Surface;
/// The minimal per-op plane table two solids share: the
/// concatenated plane list (`a`'s then `b`'s), the face→index map, and each solid's
/// [`crate::combinatorics::EdgeFaces`]. Built once and shared: indices into the returned `planes`/`surf_ix`
/// are common to both solids, so a vertex of `a` and a face of `b` compose in one index space.
/// Destructure it with `..` (`let PlaneSetup { planes: faces_tab, geom: planes, plane_ix, .. } = …`):
/// the tables here grow as the arrangement learns to say "plane" and "face" in different index
/// spaces, and a positional tuple made every one of those steps touch all ~25 call sites.
///
/// The plane classes (`canon`) are computed here to build `geom`/`plane_ix` and then dropped — the
/// dense `plane_ix` is the only face→plane map anything downstream needs, so the sparse union-find
/// output does not escape.
pub(crate) struct PlaneSetup {
    pub(crate) planes: Vec<FaceRow>,
    pub(crate) surf_ix: HashMap<Handle<Face>, usize>,
    pub(crate) inc_a: EdgeFaces,
    pub(crate) inc_b: EdgeFaces,
    /// Where `a`'s faces end and `b`'s begin in `planes`. The concatenation always created this
    /// boundary; it was just never written down, so every later "whose face is this?" had to
    /// rebuild it.
    pub(crate) n_a: usize,
    /// The arrangement's planes, densely indexed — see [`dense_planes`].
    pub(crate) geom: Vec<WorkingPlane>,
    /// `plane_ix[face]` is that face's class — a plane index into `geom`, or a cylinder class
    /// ([`ClassIx`]).
    pub(crate) plane_ix: Vec<ClassIx>,
    /// Whose faces each plane class carries — see [`class_owners`].
    pub(crate) class_owner: Vec<Option<SolidSide>>,
    /// The cylinder classes, in [`ClassIx::Cyl`] numbering order — empty for an all-planar
    /// boolean. Filled by the population gate, which is also what refuses the interactions this
    /// milestone does not build.
    pub(crate) cyls: Vec<WorkingCyl>,
    /// The gate's carried answer for the rulings road: `(plane class, cylinder class)`
    /// pairs allowed through **without** a clearance proof — see
    /// [`crate::combinatorics::TraceInput::crossings`]. ★ It used to say "always empty while the wall rule
    /// refuses that population" — the gate-opening cell arrived, and the `crossings.insert` below
    /// fills it for a wall whose plane holds the axis exactly.
    pub(crate) crossings: std::collections::HashSet<(usize, usize)>,
    /// **The tangent `(wall face, lateral face)` pairs** — the graze twin of `crossings`, and the
    /// reason they are two sets rather than one: a recorded *crossing* says "this plane runs within
    /// the radius, so it cuts the lateral in two rulings", and a tangency says the opposite ("it
    /// touches along one line and cuts nothing"). Merging them would make the record's own
    /// proposition false and break the four `debug_assert`s that lean on it. See [`Tangency`].
    pub(crate) tangencies: Vec<Tangency>,
    /// How this operation judges, and where its evidence goes — the two facts that belong to the
    /// operation rather than to any one plane. The caller pairs them with a table to make a
    /// [`Judge`].
    pub(crate) standard: Standard,
    pub(crate) notes: Notes,
}

pub(crate) fn plane_index_setup(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<PlaneSetup, BoolError> {
    let (mut setup, cyl_surfs) = plane_index_setup_inner(model, a, b)?;
    // ★ The cylinder door: the population gate decides **by name** what stands in the way (an
    // oblique cut, a wall touching the lateral surface, an undecidable pair), and what passes now
    // goes on to be arranged.
    if !cyl_surfs.is_empty() {
        let (cyls, crossings, tangencies) = cylinder_gate(
            model,
            &cyl_surfs,
            &setup.geom,
            &setup.planes,
            &setup.plane_ix,
            setup.n_a,
        )?;
        setup.cyls = cyls;
        setup.crossings = crossings;
        setup.tangencies = tangencies;
    }
    Ok(setup)
}

pub(crate) fn plane_index_setup_inner(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(PlaneSetup, Vec<Handle<Surface>>), BoolError> {
    let t = Watch::new();
    let mut planes = collect_planes(model, a)?;
    let n_a = planes.len();
    planes.extend(collect_planes(model, b)?);
    t.charge(Sub::Collect);
    let t = Watch::new();
    let standard = standard_for(&planes);
    t.charge(Sub::Std);
    let notes = Notes::new();
    if standard.prec > JUDGE_PREC_CAP {
        return Err(reject(RejectReason::PrecisionBudget {
            needed: standard.prec,
            cap: JUDGE_PREC_CAP,
        }));
    }
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        // Synthetic faces are appended later, after this table is built; every entry here is real.
        surf_ix.insert(pi.face().expect("collect_planes yields real faces"), i);
    }
    let t = Watch::new();
    let inc_a = edge_faces(model, a, &surf_ix)?;
    let inc_b = edge_faces(model, b, &surf_ix)?;
    t.charge(Sub::Edges);
    // One judging context for the whole operation: the witnesses, the standard they are held to,
    // and where the evidence goes. The face table judges first (it is what *defines* the plane
    // classes), then the dense plane table inherits the same three.
    let t = Watch::new();
    let canon = plane_classes(&Judge::new(&planes, standard, &notes));
    t.charge(Sub::Classes);
    let t = Watch::new();
    let (geom, plane_ix, cyl_surfs) = dense_planes(&planes, &canon);
    let class_owner = class_owners(&plane_ix, n_a, geom.len());
    t.charge(Sub::Dense);
    Ok((
        PlaneSetup {
            planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            geom,
            plane_ix,
            class_owner,
            cyls: Vec::new(),
            crossings: std::collections::HashSet::new(),
            tangencies: Vec::new(),
            standard,
            notes,
        },
        cyl_surfs,
    ))
}
