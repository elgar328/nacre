use super::*;
/// Apply one operation to `model`, returning the handles it created. Does not
/// rebuild the adjacency cache (do that once after a batch — see [`replay`]).
pub fn apply(model: &mut Model, op: &Operation) -> Result<OpOutput, OpError> {
    match op {
        Operation::DatumPlane { def } => {
            let (plane, frame) = datum_plane(model, def)?;
            Ok(OpOutput::DatumPlane { plane, frame })
        }
        Operation::Extrude {
            frame,
            profile,
            dist,
        } => {
            let (solid, faces) = extrude_on_frame(model, frame, profile, *dist)?;
            Ok(OpOutput::Extrude { solid, faces })
        }
        Operation::Boolean { kind, a, b } => {
            let solids = boolean(model, *kind, *a, *b).map_err(OpError::Boolean)?;
            Ok(OpOutput::Boolean { solids })
        }
        Operation::Transform { solid, isometry } => {
            let out = transform(model, *solid, isometry)?;
            Ok(OpOutput::Transform { solid: out })
        }
        Operation::Copy { solid } => {
            let out = crate::transform::copy(model, *solid)?;
            Ok(OpOutput::Copy { solid: out })
        }
        Operation::Mirror {
            solid,
            axis,
            offset,
        } => {
            let out = crate::transform::mirror(model, *solid, *axis, *offset)?;
            Ok(OpOutput::Mirror { solid: out })
        }
    }
}

/// Replay an operation log into a fresh model. Deterministic: the same log
/// reproduces the same model down to handle indices.
///
/// ★★ **A log's handles are an index vocabulary, and this is where they are re-anchored.**
/// Every [`Operation`] but a `Stated` datum names a cell by `Handle`, and those handles belong to the
/// model the log was *recorded* against — a different arena from the one being built here. A
/// `Handle`'s identity is its index (`Store`'s manual `Eq`/`Hash` use nothing else), so the
/// index is the part that carries meaning across models, and [`nacre_store::Store::handle_at`]
/// turns it back
/// into a handle of *this* model, one operation at a time.
///
/// **Why one operation at a time, and not a pre-pass**: operation *N*'s handle names a cell
/// operation *N−1* created, so nothing can be checked before the walk. That also makes it
/// atomic for free — a re-anchoring failure happens before the operation pushes anything, and
/// the half-built `model` is local, so the caller never sees a partly-mutated arena.
///
/// ★ **The premise this rests on: the log is the whole history.** Cells put into a model
/// outside the log (a fixture's `apply` with no log beside it, a direct `push_*`) shift every
/// later index, and so
/// does a *late* reject — a boolean that pushed cells before declining leaves them in the
/// append-only arena while the log has no entry for them. [`OpError::LogHandleOutOfRange`] catches only the case where the index
/// runs off the end; an index that lands on a real-but-wrong cell cannot be detected here. So a
/// session that keeps recording after a late reject must rebuild from its log first.
///
/// **`apply` does not re-anchor**, deliberately: its model belongs to the caller, so its
/// handles do too, and quietly re-anchoring there would launder a genuinely foreign handle and
/// destroy the cross-model guard that catches it.
pub fn replay(ops: &[Operation]) -> Result<Model, OpError> {
    let mut model = Model::new();
    for op in ops {
        let op = rebind(&model, op)?;
        apply(&mut model, &op)?;
    }
    model.rebuild_adjacency();
    Ok(model)
}

/// The operation with its handles re-anchored onto `model` — borrowed when there is nothing to
/// re-anchor, which since the plane-handle vocabulary landed is only a `Stated` datum.
fn rebind<'a>(model: &Model, op: &'a Operation) -> Result<Cow<'a, Operation>, OpError> {
    let surface = |h: Handle<Surface>| {
        model
            .surface_handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Surface,
                index: h.index(),
            })
    };
    let vertex = |h: Handle<Vertex>| {
        model
            .vertex_handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Vertex,
                index: h.index(),
            })
    };
    let solid = |h: Handle<Solid>| {
        model
            .solid_handle_at(h.index())
            .ok_or(OpError::LogHandleOutOfRange {
                cell: LogCell::Solid,
                index: h.index(),
            })
    };
    Ok(match op {
        // ★ `Extrude` carries a handle too — its frame's plane — so it is remapped like the rest.
        Operation::Extrude {
            frame,
            profile,
            dist,
        } => Cow::Owned(Operation::Extrude {
            frame: frame.rebound(surface(frame.plane())?),
            profile: profile.clone(),
            dist: *dist,
        }),
        Operation::Boolean { kind, a, b } => Cow::Owned(Operation::Boolean {
            kind: *kind,
            a: solid(*a)?,
            b: solid(*b)?,
        }),
        Operation::Transform { solid: s, isometry } => Cow::Owned(Operation::Transform {
            solid: solid(*s)?,
            isometry: *isometry,
        }),
        Operation::Mirror {
            solid: s,
            axis,
            offset,
        } => Cow::Owned(Operation::Mirror {
            solid: solid(*s)?,
            axis: *axis,
            offset: *offset,
        }),
        Operation::Copy { solid: s } => Cow::Owned(Operation::Copy { solid: solid(*s)? }),
        // A `Stated` datum names no cell — it is three rational points and nothing else.
        Operation::DatumPlane {
            def: DatumDef::Stated(_),
        } => Cow::Borrowed(op),
        Operation::DatumPlane {
            def: DatumDef::ThroughVertices(vs),
        } => {
            let mut out = *vs;
            for v in out.iter_mut() {
                *v = vertex(*v)?;
            }
            Cow::Owned(Operation::DatumPlane {
                def: DatumDef::ThroughVertices(out),
            })
        }
        Operation::DatumPlane {
            def: DatumDef::Offset { frame, dist },
        } => Cow::Owned(Operation::DatumPlane {
            def: DatumDef::Offset {
                frame: frame.rebound(model.surface_handle_at(frame.plane.index()).ok_or(
                    OpError::LogHandleOutOfRange {
                        cell: LogCell::Surface,
                        index: frame.plane.index(),
                    },
                )?),
                dist: *dist,
            },
        }),
    })
}
