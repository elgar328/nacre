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
//! nacre-side volume/area (`nacre-props`) is now diffed directly against OCCT
//! here (M4); the boolean `fuse|cut|common` oracle arrives with M5.

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
mod tests {
    use super::*;
    use nacre_math::{Point3, Vector3};
    use nacre_props::mass_props;
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

    // nacre-vs-OCCT diff: nacre computes volume/area analytically (nacre-props),
    // OCCT computes them independently from the same STEP. Agreement cross-checks
    // both — the two are wholly separate implementations. #[ignore]d like the
    // other DRAWEXE tests.

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn cube_props_diff_occt() {
        let mut model = Model::new();
        let solid = model.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let occt = occt_props_of(&model).unwrap();
        let nacre = mass_props(&model, solid).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "{} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn cylinder_props_diff_occt() {
        let mut model = Model::new();
        let solid = model.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let occt = occt_props_of(&model).unwrap();
        let nacre = mass_props(&model, solid).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "{} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pad_diff_occt() {
        use nacre_ops::{OpOutput, Operation, Profile2d, SketchPlane, apply};

        // Unit cube, then a 0.4-square boss of height 0.5 on the top face. The
        // padded solid's top face has a real hole (a FACE_BOUND in STEP); this
        // checks OCCT reads that holed boss as a closed solid and agrees on its
        // volume (1.08) and area (6.8) with nacre's analytic value.
        let sq = |s: f64| Profile2d {
            points: [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let mut model = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut model,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(1.0),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let boss = Profile2d {
            points: [[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut model,
            &Operation::PadOnFace {
                face: faces[1],
                profile: boss,
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };

        let occt = occt_props_of(&model).unwrap();
        let nacre = mass_props(&model, solid).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "{} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_diff_occt() {
        use nacre_ops::{OpOutput, Operation, Profile2d, SketchPlane, apply};

        // Unit cube, then a 0.4-square pocket of depth 0.5 in the top face. The
        // inward walls remove material; this checks OCCT reads the holed,
        // concave solid and agrees (volume 0.92, area 6.8) with nacre.
        let sq = |s: f64| Profile2d {
            points: [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let mut model = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut model,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(1.0),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let pocket = Profile2d {
            points: [[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut model,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: pocket,
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };

        let occt = occt_props_of(&model).unwrap();
        let nacre = mass_props(&model, solid).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "{} vs {}",
            nacre.area,
            occt.area
        );
    }

    // --- boolean oracle (M5) ---

    /// Two overlapping unit boxes A = [0,1]³, B = [0.5,1.5]³ (overlap [0.5,1]³ =
    /// 0.125). Hand-computable boolean volumes: union 1.875, difference A−B
    /// 0.875, intersection 0.125 — so they calibrate the OCCT boolean harness
    /// end to end (helper commands + single-solid export + parse).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn boolean_volumes_match_occt() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
        let cut = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
        let common = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
        assert!(approx(fuse.volume, 1.875), "fuse {}", fuse.volume);
        assert!(approx(cut.volume, 0.875), "cut {}", cut.volume);
        assert!(approx(common.volume, 0.125), "common {}", common.volume);
    }

    /// Disjoint boxes fuse to a compound whose total volume is the sum — a sanity
    /// check that the harness handles a non-overlapping (compound) result.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn disjoint_fuse_sums_volumes() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([2.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
        assert!(approx(fuse.volume, 2.0), "disjoint fuse {}", fuse.volume);
    }

    /// nacre's own `Common` result diffed against OCCT: build two overlapping
    /// cubes, ask OCCT for the intersection volume, and compare it to
    /// `mass_props` of the solid nacre's half-space enumeration produced.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn common_result_volume_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        // OCCT ground truth from the inputs (before nacre supersedes them).
        let occt = occt_boolean_of(&m, OcctBool::Common, a, b).unwrap();
        let r = boolean(&mut m, BoolKind::Common, a, b).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// nacre's `Fuse` and `Cut` results diffed against OCCT for two overlapping
    /// cubes (M5-c4 face clipping).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn fuse_cut_result_volume_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let boxes = || {
            let mut m = Model::new();
            let a = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([1.0, 1.0, 1.0]),
            );
            let b = m.add_cuboid(
                Point3::from_array([0.5, 0.5, 0.5]),
                Point3::from_array([1.5, 1.5, 1.5]),
            );
            (m, a, b)
        };
        for kind in [BoolKind::Fuse, BoolKind::Cut] {
            let occt_kind = match kind {
                BoolKind::Fuse => OcctBool::Fuse,
                BoolKind::Cut => OcctBool::Cut,
                BoolKind::Common => unreachable!(),
            };
            let (mut m, a, b) = boxes();
            let occt = occt_boolean_of(&m, occt_kind, a, b).unwrap();
            let r = boolean(&mut m, kind, a, b).unwrap();
            let nacre = mass_props(&m, r).unwrap();
            assert!(
                approx(nacre.volume, occt.volume),
                "{kind:?}: {} vs {}",
                nacre.volume,
                occt.volume
            );
        }
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn containment_cut_cavity_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        // A = [0,3]³ with B = [1,2]³ strictly inside ⇒ A − B is a hollow solid
        // (an internal void). OCCT diffs the two cavity-free inputs; nacre builds
        // the cavity independently, so the volume agreement is non-self-referential.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        // Capture OCCT's answer before the boolean supersedes the inputs.
        let occt = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, a, b).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "cavity cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn hollow_solid_step_volume_matches_occt() {
        // A = [0,4]³ (64) with a concentric B = [1,3]³ (8) void ⇒ material 56.
        // OCCT reads nacre's BREP_WITH_VOIDS export and computes the material
        // volume — the true gate on the exported void orientation. A flipped
        // void would read as 72 (= V_A + V_B), which the face-count round-trip
        // in nacre-step cannot catch.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([3.0; 3]));
        let b_outer = m.solids.get(b).outer;
        let void = m.reversed_shell(b_outer);
        let a_outer = m.solids.get(a).outer;
        let hollow = m.push_solid(Solid {
            outer: a_outer,
            cavities: vec![void],
        });
        m.live_solids.retain(|&s| s == hollow);

        let occt = occt_props_of(&m).unwrap();
        let nacre = mass_props(&m, hollow).unwrap();
        assert!(
            (nacre.volume - 56.0).abs() < 1e-9,
            "nacre volume {}",
            nacre.volume
        );
        assert!(
            approx(nacre.volume, occt.volume),
            "hollow volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn nonconvex_containment_cut_matches_occt() {
        use nacre_math::Point2;
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean};
        // A concave L-prism with a box strictly inside its bottom bar — a
        // non-convex containment cut ⇒ the L with an internal box void. OCCT reads
        // both (concave) inputs and cuts them independently of nacre's
        // point-in-polyhedron classification.
        let mut m = Model::new();
        let l_profile = Profile2d {
            points: [
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ]
            .iter()
            .map(|&p| Point2::from_array(p))
            .collect(),
        };
        let OpOutput::Extrude { solid: l, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: l_profile,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output");
        };
        let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "non-convex containment cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// A concave L-prism (volume 3) plus a box. Rebuilt per test because `boolean`
    /// supersedes its operands.
    fn l_prism_and_box(lo: [f64; 3], hi: [f64; 3]) -> (Model, Handle<Solid>, Handle<Solid>) {
        use nacre_math::Point2;
        use nacre_ops::{OpOutput, Operation, Profile2d, SketchPlane, apply};
        let mut m = Model::new();
        let l_profile = Profile2d {
            points: [
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ]
            .iter()
            .map(|&p| Point2::from_array(p))
            .collect(),
        };
        let OpOutput::Extrude { solid: l, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: l_profile,
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output");
        };
        let bx = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
        (m, l, bx)
    }

    /// The box bites the L's convex corner `(2, 0)`: a single chord, monotone bends.
    /// The non-convex overlap fixture of M5-d2.
    fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        l_prism_and_box([1.3, -0.3, 0.2], [2.4, 0.4, 1.4])
    }

    /// The box straddles the L's *reflex* corner `(1, 1)`: one chord with one reflex
    /// bend, and that bend projects outside its chord's endpoints (M5-d3 cell 3d).
    fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        l_prism_and_box([0.6, 0.6, 0.2], [1.6, 1.6, 1.4])
    }

    /// The box crosses the reflex corner and pops out the L's top: on its bottom face
    /// the seam is a staircase whose two bends turn opposite ways. `strict` rejected
    /// this until M5-d3 cell 3e-1; the reconstructed face is a correct simple polygon.
    fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
        l_prism_and_box([0.5, 0.5, 0.2], [2.5, 1.5, 1.2])
    }

    /// nacre's non-convex single-chord overlap (M5-d2) vs OCCT, `Cut`. The seam is
    /// a real boundary crossing, so this scores `overlap_fuse_cut`'s exact
    /// classification and seam reconstruction against an independent kernel.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn nonconvex_overlap_cut_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_corner_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "non-convex overlap cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The `Fuse` counterpart — the box protrudes past the L, so the union grows.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn nonconvex_overlap_fuse_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_corner_box();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "non-convex overlap fuse: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The reflex-corner bite (M5-d3 cell 3d). The arc's single bend projects *outside*
    /// its chord's endpoints, so the old projection sort was ordering it by luck. Fills
    /// the OCCT gap that cell left.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn reflex_bite_cut_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_reflex_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "reflex bite cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The folded (staircase) arc admitted by M5-d3 cell 3e-1, `Cut`.
    ///
    /// **This is the only gate on that cell.** A folded arc mis-ordered would build a
    /// self-intersecting face, and `validate` accepts one — still manifold, Euler holds.
    /// Only an independent kernel's volume says the face is right.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn folded_arc_cut_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_popup_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "folded arc cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The `Fuse` counterpart of the folded arc.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn folded_arc_fuse_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_popup_box();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "folded arc fuse: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// A stub standing in the L's top face, its footprint strictly inside that face.
    /// The seam is a closed loop in the face interior, so the result has a face with
    /// an inner loop — the shape M5-d3 cell 3f-1 opens.
    fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
        l_prism_and_box([0.3, 0.3, 0.5], [0.7, 0.7, 1.5])
    }

    /// nacre's first boolean result carrying a *hole* (M5-d3 cell 3f-1) vs OCCT: the
    /// L with a blind pocket. Area is scored alongside volume here — the hole is an
    /// area-visible feature, and nacre's own gates all read the same rings, so an
    /// independent kernel is what makes the hole's size a real claim.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn blind_dimple_cut_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_dimple();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "blind dimple cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "blind dimple cut area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The `Fuse` counterpart — a boss on the L. The same face gains the same hole,
    /// so the two agree on area while their volumes straddle the L's own 3.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn blind_dimple_fuse_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let (mut m, l, bx) = l_and_dimple();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Fuse, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "blind dimple fuse volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "blind dimple fuse area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// nacre's hole-aware classification (M5-d3 cell 3a) vs OCCT: a pocketed cube
    /// cut by a box sitting wholly inside the pocket void. The two solids are
    /// disjoint, so the answer is the pocketed cube untouched — but only if the
    /// classifier reads the lid's inner loop. Fanning the lid's outer ring alone
    /// put the box's corners on both sides of the boundary.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocketed_cut_in_the_void_matches_occt() {
        use nacre_math::Point2;
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean};
        let mut m = Model::new();
        let sq = |pts: [[f64; 2]; 4]| Profile2d {
            points: pts.iter().map(|&p| Point2::from_array(p)).collect(),
        };
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output");
        };
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1], // top cap
                profile: sq([[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]),
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!("pocket yields PocketOnFace output");
        };
        let bx = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.6]),
            Point3::from_array([0.6, 0.6, 0.9]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Cut, solid, bx).unwrap();
        let r = boolean(&mut m, BoolKind::Cut, solid, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "pocketed cut in the void: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// nacre's coincident-coplanar merge (M5-c5) vs OCCT: two cubes stacked on a
    /// shared z=1 face fuse to a 1×1×2 box.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn stacked_fuse_matches_occt() {
        use nacre_ops::{BoolKind, boolean};
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
        let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
    }
}
