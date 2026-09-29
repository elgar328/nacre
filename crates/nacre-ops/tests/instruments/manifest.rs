//! **A self dev-dependency takes no default features.**
//!
//! A crate's integration tests link the non-test library, so a crate that wants its own
//! `test-util` there names itself as a dev-dependency (`nacre-ops`, `nacre-topo`). Written with the
//! defaults, that one line also turns the crate's default features back on under
//! `--no-default-features`: `nacre-ops`'s did, so `cargo test -p nacre-ops --no-default-features`
//! — the serial build the wasm playground takes — ran the parallel build instead, green, for as
//! long as the line stood. The property lives in `Cargo.toml` and the pre-commit hook builds only the
//! defaults, so it is asserted against the manifests themselves, the facade's shape
//! (`the_ops_dependency_stays_default_free_so_consumers_can_drop_rayon`). The `cargo tree` gate line
//! in `overview.md` is the other half: it asks cargo's own resolution, which also sees a path that
//! does not go through a manifest line.

use std::path::{Path, PathBuf};

/// Every workspace crate's manifest, `crates/*/Cargo.toml`.
fn manifests() -> Vec<PathBuf> {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/");
    let mut out: Vec<PathBuf> = std::fs::read_dir(crates)
        .expect("crates/")
        .map(|e| e.expect("entry").path().join("Cargo.toml"))
        .filter(|p| p.is_file())
        .collect();
    out.sort();
    out
}

/// The crate's own dev-dependency line, if it names itself (`<name> = { path = "." … }`).
fn self_dev_dependency(manifest: &str) -> Option<String> {
    let name = manifest
        .lines()
        .find_map(|l| l.strip_prefix("name = "))?
        .trim()
        .trim_matches('"');
    let mut in_dev = false;
    for line in manifest.lines() {
        if line.starts_with('[') {
            in_dev = line.trim() == "[dev-dependencies]";
            continue;
        }
        if in_dev && line.starts_with(&format!("{name} =")) && line.contains("path = \".\"") {
            return Some(line.to_string());
        }
    }
    None
}

#[test]
fn a_self_dev_dependency_takes_no_default_features() {
    let mut found = Vec::new();
    for path in manifests() {
        let text = std::fs::read_to_string(&path).expect("a manifest");
        let Some(line) = self_dev_dependency(&text) else {
            continue;
        };
        assert!(
            line.contains("default-features = false"),
            "{}: a self dev-dependency with its defaults turns them back on under \
             `--no-default-features`: {line}",
            path.display()
        );
        let features = line.split("features = [").nth(1).unwrap_or("");
        assert!(
            !features.contains("\"default\"") && !features.contains("\"parallel\""),
            "{}: the self dev-dependency asks for the defaults by name: {line}",
            path.display()
        );
        found.push(path);
    }
    // Not vacuous: the line this exists for is among those read.
    assert!(
        found.iter().any(|p| p.ends_with("nacre-ops/Cargo.toml")),
        "nacre-ops's self dev-dependency was not found: {found:?}"
    );
}
