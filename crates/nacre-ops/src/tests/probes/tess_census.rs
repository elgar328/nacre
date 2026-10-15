//! **Can the kernel's own output be drawn?** — asked here because nowhere else asked it.
//!
//! Every other census in this crate watches a decision the boolean makes ([`crate::reject_census`],
//! `ruling_probe`, `cyl_chart::probe`); this one asks whether the solid that comes out can be
//! meshed, and asserts that every one does, with no exemption. It sits on the boolean's success
//! exit (`boolean_with_report`), the one place both public entry points pass through.
//!
//! ★ **Scope, plainly**: `#[cfg(test)]` is this crate's *unit* tests. `tests/*.rs` do **not** carry
//! this hook; the frozen census corpus asserts the same of its own rows (`tests/census.rs`'s
//! `record`), and a front-door family that needs the claim asserts it itself
//! (`edge_on_a_ruling`'s `drawn`). It cannot be put on every build — the product does not depend
//! on `nacre-tess`.
//!
//! ★ A panic inside `tessellate` is left to escape. Swallowing it would make the census say
//! "meshed" about a model that killed the mesher.
//!
//! ★★★ The **invariant** is checked in [`tess_census::record`] rather than in a test that reads
//! the ledger afterwards — see the note there. The test's job is only to say the census ran.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::assembly`] mounts it with `#[path]` as `tess_census`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use crate::ledger::Ledger;
use nacre_topo::Model;

/// One entry per boolean result the census meshed.
pub(crate) static MESHED: Ledger<()> = Ledger::new();

pub(crate) fn record(model: &Model) {
    let cfg = nacre_tess::TessConfig::default();
    let mesh = nacre_tess::tessellate(model, &cfg);
    // ★ **And it lies on its surface**: no lateral triangle leans off the cylinder by more than the
    // budget lets an edge turn (`leaning_laterals`) — a lattice point placed a few ulps off the
    // boundary makes a sliver leaning up to 90° off, with nothing else failing.
    if let Ok(t) = &mesh {
        let leaning = nacre_tess::leaning_laterals(model, t, cfg.max_angle_deg);
        assert!(
            leaning.is_empty(),
            "a boolean built a solid whose mesh leans off a lateral: {leaning:?}"
        );
    }
    // ★★★★★ **The claim is asserted here, where the fact exists — not in a later test.**
    // A `#[test]` that reads the ledger sees only the booleans that ran *before* it, so a solid
    // built by any later test would go unchecked. Asserting at the record makes the coverage
    // total and names the offending test in the panic.
    if let Err(e) = &mesh {
        panic!("a boolean built a solid the mesher refuses: {e:?}");
    }
    MESHED.push(());
}
