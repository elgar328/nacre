//! **Can the kernel's own output be drawn?** — asked here because nowhere else asked it.
//!
//! Every other census in this crate watches a decision the boolean makes ([`crate::reject_census`],
//! `ruling_probe`, `cyl_chart::probe`). None watched whether the solid that comes out can be
//! meshed, and unwatched the answer was **no, twice in 2129** — a tangency whose face interior
//! pinches — with nothing in the suite tessellating those two fixtures. With
//! that pinch bridged the answer is **every one meshes**, and the assertion says so with no
//! exemption.
//! This hook sits on the single success exit above, which is the one place both public entry
//! points pass through.
//!
//! ★★ **It runs before `rebuild_adjacency`, deliberately measured**: the failures reproduce on a
//! rebuilt model too, so this is watching the result and not an artifact of when it looks.
//!
//! ★ **Scope, plainly**: `#[cfg(test)]` is this crate's *unit* tests. `tests/*.rs` do **not** carry
//! this hook, so "2129 booleans" is the lib suite's number and not the workspace's. The frozen census
//! corpus asserts the same of its own rows (`tests/census.rs`'s `record`); `perf` and
//! `invariants/pipeline` do not.
//!
//! ★ A panic inside `tessellate` is left to escape. Swallowing it would make the census say
//! "meshed" about a model that killed the mesher.
//!
//! ★★★ The **invariant** is checked in [`tess_census::record`] rather than in a test that reads
//! the vector afterwards — see the note there. The test's job is only to say the census ran and
//! that the known population is still in it.
//!
//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::assembly`] mounts it with `#[path]` as `tess_census`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use crate::ledger::Ledger;
use nacre_tess::TessError;
use nacre_topo::Model;

/// One entry per boolean that returned a solid: the triangle count, or why not.
pub(crate) static MESHED: Ledger<Result<usize, TessError>> = Ledger::new();

pub(crate) fn record(model: &Model) {
    let r = nacre_tess::tessellate(model, &nacre_tess::TessConfig::default())
        .map(|t| t.triangles.len());
    // ★★★★★ **The claim is asserted here, where the fact exists — not in a later test.**
    // A `#[test]` that reads this vector sees only the booleans that ran *before* it (☑
    // measured: 764 of them, under `--test-threads=1`, because the suite runs in name order),
    // so a solid built by any later test would go unchecked. Asserting at the record makes the
    // coverage total and names the offending test in the panic instead of a distant census.
    // ★ No exemption: pinched caps
    // are bridged, so every solid a boolean returns
    // meshes, or the census says which test built the one that does not.
    assert!(
        r.is_ok(),
        "a boolean built a solid the mesher refuses: {r:?}"
    );
    MESHED.push(r);
}
