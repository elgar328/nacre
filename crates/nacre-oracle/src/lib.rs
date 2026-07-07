//! OCCT oracle harness — nacre's primary defense against "plausible but wrong"
//! geometry (design §7).
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
//! boolean `fuse|cut|common` oracle and a nacre-side volume for direct diff
//! arrive with M5.

use nacre_topo::Model;
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
    /// `bbox_max <x y z>`. Order-independent; every field is required.
    fn parse(stdout: &str) -> Result<OcctProps, OracleError> {
        let mut volume = None;
        let mut area = None;
        let mut faces = None;
        let mut bbox_min = None;
        let mut bbox_max = None;

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
                _ => {} // ignore unknown keys — forward-compatible with new fields
            }
        }

        Ok(OcctProps {
            volume: volume.ok_or_else(|| miss("volume"))?,
            area: area.ok_or_else(|| miss("area"))?,
            faces: faces.ok_or_else(|| miss("faces"))?,
            bbox_min: bbox_min.ok_or_else(|| miss("bbox_min"))?,
            bbox_max: bbox_max.ok_or_else(|| miss("bbox_max"))?,
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

#[cfg(test)]
mod tests {
    use super::*;
    use nacre_math::{Point3, Vector3};
    use std::f64::consts::PI;

    /// Combined relative-or-absolute float comparison, sized to DRAWEXE's output
    /// precision — it prints ~6 significant figures, so an exact value can land
    /// ~5e-7 relative away (e.g. 20π → `62.8319`). A 1e-4 relative band clears
    /// that rounding noise by 100× while still catching any real geometry error
    /// (those miss by percents, not parts-per-thousand).
    fn approx(a: f64, b: f64) -> bool {
        let diff = (a - b).abs();
        diff <= 1e-6 || diff <= 1e-4 * a.abs().max(b.abs())
    }

    #[test]
    fn parse_reads_all_fields_order_independent() {
        let stdout = "\
faces 6
volume 24
area 52
bbox_max 2 3 4
bbox_min 0 0 0
";
        let p = OcctProps::parse(stdout).unwrap();
        assert_eq!(p.volume, 24.0);
        assert_eq!(p.area, 52.0);
        assert_eq!(p.faces, 6);
        assert_eq!(p.bbox_min, [0.0, 0.0, 0.0]);
        assert_eq!(p.bbox_max, [2.0, 3.0, 4.0]);
    }

    #[test]
    fn parse_reports_missing_field() {
        let stdout = "volume 24\narea 52\nfaces 6\nbbox_min 0 0 0\n";
        match OcctProps::parse(stdout) {
            Err(OracleError::Parse(msg)) => assert!(msg.contains("bbox_max")),
            other => panic!("expected Parse error, got {other:?}"),
        }
    }

    #[test]
    fn parse_reports_non_numeric() {
        let stdout = "volume oops\narea 52\nfaces 6\nbbox_min 0 0 0\nbbox_max 2 3 4\n";
        assert!(matches!(
            OcctProps::parse(stdout),
            Err(OracleError::Parse(_))
        ));
    }

    // End-to-end oracle tests: they shell out to DRAWEXE, so they are #[ignore]d
    // (the pre-commit hook still compiles + clippy + fmt them). Run on a machine
    // with `brew install opencascade`:  cargo test -p nacre-oracle -- --ignored
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn cuboid_matches_occt() {
        let mut model = Model::new();
        model.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let p = occt_props_of(&model).unwrap();
        assert!(approx(p.volume, 24.0), "volume {}", p.volume); // 2·3·4
        assert!(approx(p.area, 52.0), "area {}", p.area); // 2(6+8+12)
        assert_eq!(p.faces, 6);
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn cylinder_matches_occt() {
        let mut model = Model::new();
        model.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let p = occt_props_of(&model).unwrap();
        // OCCT computes these analytically from the CYLINDRICAL_SURFACE, so a
        // match confirms our curved STEP really is a cylinder — and, since a
        // positive volume needs a closed outward-oriented solid, this automates
        // the orientation/validity check M3 deferred to manual FreeCAD.
        assert!(approx(p.volume, 20.0 * PI), "volume {}", p.volume); // π·r²·h
        assert!(approx(p.area, 28.0 * PI), "area {}", p.area); // 2πr² + 2πr·h
        assert_eq!(p.faces, 3); // lateral + 2 caps
    }
}
