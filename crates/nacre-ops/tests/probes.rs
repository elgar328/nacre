//! Regression probes: each file is one shape that once broke the kernel, kept as a named test.

#[path = "probes/collinear_loop_points.rs"]
mod collinear_loop_points;
#[path = "probes/concurrent_line.rs"]
mod concurrent_line;
#[path = "probes/contact_separates.rs"]
mod contact_separates;
#[path = "probes/invariant_planes.rs"]
mod invariant_planes;
#[path = "probes/knife_edge.rs"]
mod knife_edge;
#[path = "probes/rotation_sweep.rs"]
mod rotation_sweep;
#[path = "probes/self_touch.rs"]
mod self_touch;
#[path = "probes/vertex_tolerance.rs"]
mod vertex_tolerance;
#[path = "probes/void_label.rs"]
mod void_label;
