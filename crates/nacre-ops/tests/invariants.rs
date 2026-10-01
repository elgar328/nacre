//! Invariants that hold on every production path: replay determinism, the ops -> validate -> step/tess
//! pipeline, edge carriers, recorded plane points, sketch frames, datum planes, vertex realization,
//! one plane held as two handles, the topology layer's doors and live set on a product box.

#[path = "support/fixtures.rs"]
mod fixtures;
#[path = "support/stated.rs"]
mod stated;

#[path = "invariants/cylinder_seam_model.rs"]
mod cylinder_seam_model;
#[path = "invariants/datum_plane.rs"]
mod datum_plane;
#[path = "invariants/edge_carriers.rs"]
mod edge_carriers;
#[path = "invariants/one_plane_two_handles.rs"]
mod one_plane_two_handles;
#[path = "invariants/pipeline.rs"]
mod pipeline;
#[path = "invariants/points_coverage.rs"]
mod points_coverage;
#[path = "invariants/realize_vertex.rs"]
mod realize_vertex;
#[path = "invariants/replay.rs"]
mod replay;
#[path = "invariants/sketch_frame.rs"]
mod sketch_frame;
#[path = "invariants/sketch_frame_contract.rs"]
mod sketch_frame_contract;
#[path = "invariants/topo_model.rs"]
mod topo_model;
