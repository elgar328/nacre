//! OCCT oracle harness — nacre's primary defense against "plausible but wrong"
//! geometry.
//!
//! This crate exports a [`Model`] to STEP and asks OpenCASCADE (OCCT) to score
//! it: volume, surface area, face count, bounding box. OCCT is the *oracle*
//! (the answer key), reached **out-of-process** through `tools/occt-helper` +
//! STEP files. OCCT is **never linked into the kernel** — this is a dev-only
//! test harness (`version = "0.0.0"`, never published; overview OCCT rules). A
//! crash on adversarial input kills the helper subprocess, not the kernel.
//!
//! v0 is `props` only (volume/area/faces/bbox of one solid). It validates that
//! (1) our STEP is OCCT-valid, (2) a whole closed oriented solid was read (OCCT
//! reports a positive analytic volume only then), and (3) the transport works —
//! promoting the old manual FreeCAD check to an automated regression. The
//! nacre-side volume/area (`nacre-props`) is now diffed directly against OCCT
//! here (M4), and the boolean `fuse|cut|common` oracle beside it since M5.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_store::Handle;
use nacre_topo::{Model, Solid};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Mass properties of a solid, as computed by OCCT.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OcctProps {
    pub volume: f64,
    pub area: f64,
    pub faces: usize,
    pub bbox_min: [f64; 3],
    pub bbox_max: [f64; 3],
    /// OCCT's **volume** centre of mass (`vprops`), not the surface's.
    pub centroid: [f64; 3],
}

/// A failure while scoring a model against the OCCT oracle.
#[derive(Debug)]
pub enum OracleError {
    /// `tools/occt-helper` was not found or is not executable.
    HelperMissing,
    /// The helper read no shape from the STEP (missing file, empty, non-solid) —
    /// helper exit 1.
    GeometryFailed,
    /// DRAWEXE crashed on the input — helper exit 2. Carries stderr for triage.
    Crashed(String),
    /// [`nacre_step::to_step`] rejected the model (its `StepError`, stringified).
    Export(String),
    /// Filesystem or subprocess I/O failed.
    Io(String),
    /// The helper's stdout could not be parsed into [`OcctProps`].
    Parse(String),
}

impl OcctProps {
    /// Parse the helper's `props` stdout — key-value lines, one per field:
    /// `volume <v>` / `area <a>` / `faces <n>` / `bbox_min <x y z>` /
    /// `bbox_max <x y z>` / `centroid <x y z>`. Order-independent; every field
    /// is required.
    fn parse(stdout: &str) -> Result<OcctProps, OracleError> {
        let mut volume = None;
        let mut area = None;
        let mut faces = None;
        let mut bbox_min = None;
        let mut bbox_max = None;
        let mut centroid = None;

        for line in stdout.lines() {
            let mut it = line.split_whitespace();
            let Some(key) = it.next() else { continue };
            match key {
                "volume" => volume = Some(parse_f64(&mut it, "volume")?),
                "area" => area = Some(parse_f64(&mut it, "area")?),
                "faces" => {
                    let tok = it.next().ok_or_else(|| miss("faces"))?;
                    let n = tok
                        .parse::<usize>()
                        .map_err(|e| OracleError::Parse(format!("faces: {e}")))?;
                    faces = Some(n);
                }
                "bbox_min" => bbox_min = Some(parse_triple(&mut it, "bbox_min")?),
                "bbox_max" => bbox_max = Some(parse_triple(&mut it, "bbox_max")?),
                "centroid" => centroid = Some(parse_triple(&mut it, "centroid")?),
                _ => {} // ignore unknown keys — forward-compatible with new fields
            }
        }

        Ok(OcctProps {
            volume: volume.ok_or_else(|| miss("volume"))?,
            area: area.ok_or_else(|| miss("area"))?,
            faces: faces.ok_or_else(|| miss("faces"))?,
            bbox_min: bbox_min.ok_or_else(|| miss("bbox_min"))?,
            bbox_max: bbox_max.ok_or_else(|| miss("bbox_max"))?,
            centroid: centroid.ok_or_else(|| miss("centroid"))?,
        })
    }
}

fn miss(field: &str) -> OracleError {
    OracleError::Parse(format!("missing field: {field}"))
}

fn parse_f64<'a>(it: &mut impl Iterator<Item = &'a str>, field: &str) -> Result<f64, OracleError> {
    let tok = it.next().ok_or_else(|| miss(field))?;
    tok.parse::<f64>()
        .map_err(|e| OracleError::Parse(format!("{field}: {e}")))
}

fn parse_triple<'a>(
    it: &mut impl Iterator<Item = &'a str>,
    field: &str,
) -> Result<[f64; 3], OracleError> {
    let x = parse_f64(it, field)?;
    let y = parse_f64(it, field)?;
    let z = parse_f64(it, field)?;
    Ok([x, y, z])
}

/// Absolute path to the `occt-helper` script (`<crate>/../../tools/...`).
fn helper_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tools/occt-helper/occt-helper")
}

/// Score STEP text against the OCCT oracle: write it to a unique temp file, run
/// `occt-helper props`, and parse the result.
pub fn occt_props(step: &str) -> Result<OcctProps, OracleError> {
    let helper = helper_path();
    if !helper.is_file() {
        return Err(OracleError::HelperMissing);
    }

    // Unique temp path: pid + a process-lifetime counter (no external deps).
    static SEQ: AtomicUsize = AtomicUsize::new(0);
    let n = SEQ.fetch_add(1, Ordering::Relaxed);
    let mut path = std::env::temp_dir();
    path.push(format!("nacre-oracle-{}-{}.step", std::process::id(), n));

    std::fs::write(&path, step).map_err(|e| OracleError::Io(e.to_string()))?;

    let output = std::process::Command::new(&helper)
        .arg("props")
        .arg(&path)
        .output();

    // Best-effort cleanup; ignore removal errors.
    let _ = std::fs::remove_file(&path);

    let output = output.map_err(|e| OracleError::Io(e.to_string()))?;
    match output.status.code() {
        Some(0) => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            OcctProps::parse(&stdout)
        }
        Some(1) => Err(OracleError::GeometryFailed),
        Some(2) => Err(OracleError::Crashed(
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )),
        other => Err(OracleError::Io(format!(
            "occt-helper exited with unexpected status {other:?}"
        ))),
    }
}

/// Convenience: export a [`Model`] to STEP, then score it with [`occt_props`].
pub fn occt_props_of(model: &Model) -> Result<OcctProps, OracleError> {
    let step = nacre_step::to_step(model).map_err(|e| OracleError::Export(format!("{e:?}")))?;
    occt_props(&step)
}

/// Which binary boolean the OCCT oracle should run. Its own enum — `nacre-ops`
/// (which will define `BoolKind`) is only a dev-dependency here, so it cannot
/// appear in this crate's public API.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OcctBool {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

impl OcctBool {
    /// The `occt-helper` command name.
    fn as_str(self) -> &'static str {
        match self {
            OcctBool::Fuse => "fuse",
            OcctBool::Cut => "cut",
            OcctBool::Common => "common",
        }
    }
}

/// Score the OCCT boolean of two single-solid STEP texts: write each to a unique
/// temp file, run `occt-helper <fuse|cut|common>`, and parse the result's props.
///
/// Each input must be a single-solid STEP ([`nacre_step::to_step_solid`]); a
/// multi-root STEP would have OCCT operate on only its first shape. This is the
/// M5 boolean ground truth — the answer key for the nacre `PolyhedralBoolean`.
pub fn occt_boolean(kind: OcctBool, a_step: &str, b_step: &str) -> Result<OcctProps, OracleError> {
    let helper = helper_path();
    if !helper.is_file() {
        return Err(OracleError::HelperMissing);
    }

    // Two unique temp paths: pid + a process-lifetime counter (no external deps).
    static SEQ: AtomicUsize = AtomicUsize::new(0);
    let stamp = |tag: char| {
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!(
            "nacre-oracle-{}-{}{}.step",
            std::process::id(),
            tag,
            n
        ));
        path
    };
    let a_path = stamp('a');
    let b_path = stamp('b');
    std::fs::write(&a_path, a_step).map_err(|e| OracleError::Io(e.to_string()))?;
    std::fs::write(&b_path, b_step).map_err(|e| OracleError::Io(e.to_string()))?;

    let output = std::process::Command::new(&helper)
        .arg(kind.as_str())
        .arg(&a_path)
        .arg(&b_path)
        .output();

    // Best-effort cleanup; ignore removal errors.
    let _ = std::fs::remove_file(&a_path);
    let _ = std::fs::remove_file(&b_path);

    let output = output.map_err(|e| OracleError::Io(e.to_string()))?;
    match output.status.code() {
        Some(0) => OcctProps::parse(&String::from_utf8_lossy(&output.stdout)),
        Some(1) => Err(OracleError::GeometryFailed),
        Some(2) => Err(OracleError::Crashed(
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )),
        other => Err(OracleError::Io(format!(
            "occt-helper exited with unexpected status {other:?}"
        ))),
    }
}

/// Convenience: export solids `a` and `b` from `model` (each as a single-solid
/// STEP) and score their OCCT boolean.
pub fn occt_boolean_of(
    model: &Model,
    kind: OcctBool,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<OcctProps, OracleError> {
    let export =
        |h| nacre_step::to_step_solid(model, h).map_err(|e| OracleError::Export(format!("{e:?}")));
    occt_boolean(kind, &export(a)?, &export(b)?)
}

#[cfg(test)]
#[path = "tests/lib/mod.rs"]
mod tests;
