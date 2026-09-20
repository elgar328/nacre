//! Instruments that measure a population and pin what they find. None of these reads process-global
//! state; the ones that do (`census`, `reject_census`, `wide_datum_cost`) keep their own binary, because
//! every boolean in the same process moves those counters.

#[path = "support/fixtures.rs"]
mod fixtures;

#[path = "instruments/plane_anchor.rs"]
mod plane_anchor;
#[path = "instruments/point_width.rs"]
mod point_width;
