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
    use nacre_ops::SketchPlane;
    use nacre_ops::boolean;
    use nacre_ops::{BoolError, BoolKind};
    use nacre_props::mass_props;
    use std::f64::consts::PI;

    /// Test shim: a boolean whose result is exactly one solid (cell 0.4 multi-solid).
    fn boolean_one(
        model: &mut Model,
        kind: BoolKind,
        a: Handle<Solid>,
        b: Handle<Solid>,
    ) -> Result<Handle<Solid>, BoolError> {
        let solids = boolean(model, kind, a, b)?;
        assert_eq!(
            solids.len(),
            1,
            "boolean_one: expected one solid, got {}",
            solids.len()
        );
        Ok(solids[0])
    }

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
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};

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
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};

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

    /// A pad whose footprint overhangs one face edge (a boss cantilever). This is the first time
    /// the `PadOnFace` pipe produces an overhang — it routes to the overhang Fuse sidecar. OCCT
    /// scores the cantilevered boss.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn overhang_pad_matches_occt() {
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};
        let mut model = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut model,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d {
                    points: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
                        .iter()
                        .map(|&p| nacre_math::Point2::from_array(p))
                        .collect(),
                },
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        // Footprint world x in [0.25,0.75], y in [-0.25,0.75] - overhangs the y=0 edge.
        let boss = Profile2d {
            points: [[-0.25, -0.25], [0.75, -0.25], [0.75, 0.25], [-0.25, 0.25]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut model,
            &Operation::PadOnFace {
                face: faces[1],
                profile: boss,
                dist: 1.0,
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
    }

    /// A blind pocket whose footprint overhangs one face edge (an edge slot open to the side).
    /// First overhang through the `PocketOnFace` pipe - routes to the overhang Cut sidecar. OCCT
    /// scores the slotted solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn overhang_pocket_matches_occt() {
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};
        let mut model = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut model,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d {
                    points: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
                        .iter()
                        .map(|&p| nacre_math::Point2::from_array(p))
                        .collect(),
                },
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let slot = Profile2d {
            points: [[-0.25, -0.25], [0.75, -0.25], [0.75, 0.25], [-0.25, 0.25]]
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut model,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: slot,
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

    /// E1 same_ground union: A = [0,1]³ and B = [0.5,1.5]²×[0,1] overlap in volume and share the
    /// z=0 / z=1 planes with overlapping footprints. The union is an L-footprint prism: area
    /// (1 + 1 − 0.25) × height 1 = 1.75. OCCT confirms it independently — the oracle for nacre's
    /// general 2D coplanar merge (the driver's union-cell reconstruct + interpenetrating-wall clip).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn same_ground_union_matches_occt() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.0]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
        assert!(
            approx(fuse.volume, 1.75),
            "same_ground fuse {}",
            fuse.volume
        );
    }

    /// A flush-edge pocket (B4-R1b): the cutter sits flush on two adjacent base faces (top z=1 and
    /// front y=0), so its walls are coplanar with the part's walls along the shared boundary edge.
    /// OCCT confirms base − cutter = 0.92 independently — the oracle for nacre's flush handler.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn flush_edge_pocket_cut_matches_occt() {
        let mut m = Model::new();
        let base = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let cutter = m.add_cuboid(
            Point3::from_array([0.3, 0.0, 0.5]),
            Point3::from_array([0.7, 0.4, 1.0]),
        );
        let cut = occt_boolean_of(&m, OcctBool::Cut, base, cutter).unwrap();
        assert!(approx(cut.volume, 0.92), "flush cut {}", cut.volume);
    }

    /// A through-tunnel Cut (B4-R1b-part2c): a cutter spanning the bar's full height (top and bottom
    /// both flush) makes two parallel coplanar contacts, cut into a tunnel via the contained branch.
    /// OCCT confirms bar − tunnel = 0.84 independently.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn through_tunnel_cut_matches_occt() {
        let mut m = Model::new();
        let bar = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let cutter = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 0.0]),
            Point3::from_array([0.7, 0.7, 1.0]),
        );
        let cut = occt_boolean_of(&m, OcctBool::Cut, bar, cutter).unwrap();
        assert!(approx(cut.volume, 0.84), "tunnel cut {}", cut.volume);
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

    /// A severing `Cut` returns two nacre solids; OCCT returns a COMPOUND of two solids for the
    /// same inputs. Their aggregate volume and area agree — nacre's multi-solid output (cell 0.4)
    /// matches OCCT. The bar threads the cube and out both ends, leaving two 1×1×1 stubs.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn sever_cut_matches_occt_compound() {
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        // OCCT ground truth from the inputs, before nacre supersedes them.
        let occt = occt_boolean_of(&m, OcctBool::Cut, bar, cube).unwrap();
        let solids = boolean(&mut m, BoolKind::Cut, bar, cube).unwrap();
        assert_eq!(solids.len(), 2, "nacre severs into two solids");
        let vol: f64 = solids
            .iter()
            .map(|&s| mass_props(&m, s).unwrap().volume)
            .sum();
        let area: f64 = solids
            .iter()
            .map(|&s| mass_props(&m, s).unwrap().area)
            .sum();
        assert!(
            approx(vol, occt.volume),
            "volume {vol} vs occt {}",
            occt.volume
        );
        assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
    }

    /// A `Transform`-translated solid composes correctly into a boolean (overhaul
    /// stage 1a): translate a cube by a rational offset, then `Cut` an overlapping
    /// cube from it. The moved geometry's boolean matches OCCT on the same inputs —
    /// so the geometry-rewrite produced a boolean-valid solid, not just a
    /// volume-invariant one.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn translated_solid_cut_matches_occt() {
        use nacre_ops::{BoolKind, Operation, apply, boolean};
        use nacre_scalar::{Isometry, Rat};
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        // Translate A by [3/10, 2/5, 1/2] ⇒ A' = [0.3,1.3]×[0.4,1.4]×[0.5,1.5].
        let iso = Isometry::translation([
            Rat::new(3, 10).unwrap(),
            Rat::new(2, 5).unwrap(),
            Rat::new(1, 2).unwrap(),
        ]);
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: a,
                isometry: iso,
            },
        )
        .unwrap();
        let nacre_ops::OpOutput::Transform { solid: a2 } = out else {
            panic!("expected Transform output");
        };
        // B overlaps A' at a corner.
        let b = m.add_cuboid(Point3::from_array([0.8; 3]), Point3::from_array([1.8; 3]));
        // OCCT ground truth on the moved inputs, before nacre supersedes them.
        let occt = occt_boolean_of(&m, OcctBool::Cut, a2, b).unwrap();
        let solids = boolean(&mut m, BoolKind::Cut, a2, b).unwrap();
        let vol: f64 = solids
            .iter()
            .map(|&s| mass_props(&m, s).unwrap().volume)
            .sum();
        let area: f64 = solids
            .iter()
            .map(|&s| mass_props(&m, s).unwrap().area)
            .sum();
        assert!(
            approx(vol, occt.volume),
            "volume {vol} vs occt {}",
            occt.volume
        );
        assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
    }

    /// A `Transform`-rotated solid is a well-formed b-rep (overhaul stage 1b):
    /// rotate a cuboid 30° about Z through a rational axis point, then export just
    /// that solid and ask OCCT for its volume/area. A rigid rotation leaves both
    /// invariant, so OCCT must agree with nacre's `mass_props` — proving the
    /// rotation rewrite produced a valid solid, not merely a volume-preserving
    /// vertex shuffle. The same check on a rotated `Cut` result exercises the
    /// `Discovered`-vertex rotation path. Single-solid export (`to_step_solid`)
    /// avoids summing the superseded input still resident in the arena.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_solid_props_match_occt() {
        use nacre_ops::{Operation, apply};
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let rot30 = Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
        });

        // (a) rotate a plain cuboid.
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: c,
                isometry: rot30,
            },
        )
        .unwrap();
        let nacre_ops::OpOutput::Transform { solid: c2 } = out else {
            panic!("expected Transform output");
        };
        let occt = occt_props(&nacre_step::to_step_solid(&m, c2).unwrap()).unwrap();
        let nacre = mass_props(&m, c2).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "volume {} vs occt {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs occt {}",
            nacre.area,
            occt.area
        );

        // (b) rotate a Cut result (Discovered seam vertices → Rotated).
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        let cut = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        let out = apply(
            &mut m,
            &Operation::Transform {
                solid: cut,
                isometry: rot30,
            },
        )
        .unwrap();
        let nacre_ops::OpOutput::Transform { solid: cut2 } = out else {
            panic!("expected Transform output");
        };
        let occt = occt_props(&nacre_step::to_step_solid(&m, cut2).unwrap()).unwrap();
        let nacre = mass_props(&m, cut2).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "cut volume {} vs occt {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "cut area {} vs occt {}",
            nacre.area,
            occt.area
        );
    }

    // --- Rotated booleans vs OCCT (overhaul 3d-ii) ---
    // A boolean commutes with a rigid motion, so a rotated boolean is the rotated image of the
    // unrotated one — and OCCT, given the rotated STEP, is an *independent* ground truth for its
    // volume and area. 3d-i checked rotation-invariance (self-consistency); this cross-checks the
    // rotated seam/reconstruction against a mature kernel, catching a same-volume-wrong-topology
    // bug that invariance alone could miss.

    /// An axis-aligned rotation about the line through `(1,1,0)` by `deg` degrees (a non-90°
    /// degree makes cos/sin irrational, so the realization is `Origin::Rotated`).
    fn rot_about(axis: nacre_scalar::Axis, deg: i128) -> nacre_scalar::Isometry {
        use nacre_scalar::{Angle, Isometry, Rat, Rotation};
        Isometry::rotation(Rotation {
            axis,
            point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        })
    }

    /// Rotate both operands by every isometry in `isos` (a chain), then assert nacre's boolean
    /// matches OCCT on total volume and area. OCCT ground truth is taken *before* the nacre
    /// boolean supersedes the rotated inputs. Returns the result solids for further assertions.
    fn rotated_boolean_matches_occt(
        m: &mut Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        kind: BoolKind,
        isos: &[nacre_scalar::Isometry],
    ) -> Vec<Handle<Solid>> {
        use nacre_ops::{OpOutput, Operation, apply};
        let rot = |m: &mut Model, mut s: Handle<Solid>| -> Handle<Solid> {
            for iso in isos {
                let out = apply(
                    m,
                    &Operation::Transform {
                        solid: s,
                        isometry: *iso,
                    },
                )
                .unwrap();
                let OpOutput::Transform { solid } = out else {
                    panic!("expected Transform output");
                };
                s = solid;
                m.rebuild_adjacency();
            }
            s
        };
        let a2 = rot(m, a);
        let b2 = rot(m, b);
        let occt_kind = match kind {
            BoolKind::Fuse => OcctBool::Fuse,
            BoolKind::Cut => OcctBool::Cut,
            BoolKind::Common => OcctBool::Common,
        };
        // OCCT ground truth on the rotated inputs, before nacre supersedes them.
        let occt = occt_boolean_of(m, occt_kind, a2, b2).unwrap();
        let solids = boolean(m, kind, a2, b2).unwrap();
        let vol: f64 = solids
            .iter()
            .map(|&s| mass_props(m, s).unwrap().volume)
            .sum();
        let area: f64 = solids.iter().map(|&s| mass_props(m, s).unwrap().area).sum();
        assert!(
            approx(vol, occt.volume),
            "volume {vol} vs occt {}",
            occt.volume
        );
        assert!(approx(area, occt.area), "area {area} vs occt {}", occt.area);
        solids
    }

    /// A fully-tilted (Z then X rotation, every face normal irrational) corner-overlap `Cut`
    /// matches OCCT — the strongest cross-check of the rotated arrangement/seam.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_overlap_cut_matches_occt() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        let isos = [rot_about(Axis::Z, 30), rot_about(Axis::X, 30)];
        rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Cut, &isos);
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_overlap_fuse_matches_occt() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Fuse, &[rot_about(Axis::Z, 30)]);
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_overlap_common_matches_occt() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
        let b = m.add_cuboid(Point3::from_array([0.5; 3]), Point3::from_array([1.5; 3]));
        rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Common, &[rot_about(Axis::Z, 30)]);
    }

    /// The sever — a bar cut clean through a cube into two solids (OCCT COMPOUND) — under a full
    /// tilt. This is where a rotated `is_shell_outward` (3c-vi) bug would hide: both severed
    /// pieces must read outward, so OCCT's total volume/area confirms two material solids, not one
    /// with the other misjudged as a cavity.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_sever_cut_matches_occt() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let cube = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let isos = [rot_about(Axis::Z, 30), rot_about(Axis::X, 30)];
        let solids = rotated_boolean_matches_occt(&mut m, bar, cube, BoolKind::Cut, &isos);
        assert_eq!(solids.len(), 2, "rotated sever yields two solids");
    }

    /// A rotated `Cut` of a strictly-contained box leaves a cavity (OCCT BREP_WITH_VOIDS): the
    /// volume nets the void, and nacre records exactly one cavity.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rotated_containment_cut_makes_cavity_matches_occt() {
        use nacre_scalar::Axis;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        let solids =
            rotated_boolean_matches_occt(&mut m, a, b, BoolKind::Cut, &[rot_about(Axis::Z, 30)]);
        assert_eq!(solids.len(), 1, "one solid");
        assert_eq!(
            m.solids.get(solids[0]).cavities.len(),
            1,
            "the contained box is a cavity"
        );
    }

    /// A **re-rotated** solid (overhaul stage 1c) is still a valid b-rep: rotate a
    /// cuboid 30° about Z, then 45° about X, so its vertices carry a two-node rotation
    /// chain. A rigid re-rotation leaves volume/area invariant, so OCCT must agree with
    /// nacre's `mass_props` — confirming the re-rotation multi-pass clone (chained
    /// forest, base=root) produced a well-formed solid, not just an invariant-preserving
    /// vertex shuffle. Single-solid export avoids summing the superseded intermediates.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn rerotated_solid_props_match_occt() {
        use nacre_ops::{Operation, apply};
        use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
        let rot = |axis, deg: i128| {
            Isometry::rotation(Rotation {
                axis,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            })
        };
        let mut m = Model::new();
        let c = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let step = |m: &mut Model, s, iso| {
            let out = apply(
                m,
                &Operation::Transform {
                    solid: s,
                    isometry: iso,
                },
            )
            .unwrap();
            let nacre_ops::OpOutput::Transform { solid } = out else {
                panic!("expected Transform output");
            };
            solid
        };
        let c1 = step(&mut m, c, rot(Axis::Z, 30));
        let c2 = step(&mut m, c1, rot(Axis::X, 45));

        let occt = occt_props(&nacre_step::to_step_solid(&m, c2).unwrap()).unwrap();
        let nacre = mass_props(&m, c2).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "volume {} vs occt {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs occt {}",
            nacre.area,
            occt.area
        );
    }

    /// nacre's own `Common` result diffed against OCCT: build two overlapping
    /// cubes, ask OCCT for the intersection volume, and compare it to
    /// `mass_props` of the solid nacre's half-space enumeration produced.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn common_result_volume_matches_occt() {
        use nacre_ops::BoolKind;
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
        let r = boolean_one(&mut m, BoolKind::Common, a, b).unwrap();
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
        use nacre_ops::BoolKind;
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
            let r = boolean_one(&mut m, kind, a, b).unwrap();
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
        use nacre_ops::BoolKind;
        // A = [0,3]³ with B = [1,2]³ strictly inside ⇒ A − B is a hollow solid
        // (an internal void). OCCT diffs the two cavity-free inputs; nacre builds
        // the cavity independently, so the volume agreement is non-self-referential.
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
        // Capture OCCT's answer before the boolean supersedes the inputs.
        let occt = occt_boolean_of(&m, OcctBool::Cut, a, b).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
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
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
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
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        let (mut m, l) = l_prism();
        let bx = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
        (m, l, bx)
    }

    /// The concave L-prism alone (volume 3).
    fn l_prism() -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let l = extrude(
            &mut m,
            SketchPlane::world_xy(),
            &[
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
        );
        (m, l)
    }

    /// Extrude a closed profile 1.0 along `plane`'s normal.
    fn extrude(m: &mut Model, plane: SketchPlane, pts: &[[f64; 2]]) -> Handle<Solid> {
        extrude_dist(m, plane, pts, 1.0)
    }

    /// Extrude a closed profile `dist` along `plane`'s normal.
    fn extrude_dist(
        m: &mut Model,
        plane: SketchPlane,
        pts: &[[f64; 2]],
        dist: f64,
    ) -> Handle<Solid> {
        use nacre_math::Point2;
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};
        let profile = Profile2d {
            points: pts.iter().map(|&p| Point2::from_array(p)).collect(),
        };
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                plane,
                profile,
                dist,
            },
        )
        .unwrap() else {
            unreachable!("extrude yields Extrude output");
        };
        solid
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_corner_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_corner_box();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "non-convex overlap fuse: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The `Common` counterpart — cell 3g opened non-convex `Common` (the intersection is
    /// the corner bite, `0.224`) and retired the convex `common`. Same seam as the Cut/Fuse
    /// above, only the keep/flip table differs.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn nonconvex_overlap_common_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_corner_box();
        let occt = occt_boolean_of(&m, OcctBool::Common, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Common, l, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "vol {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// A slotted bar — a cuboid with a full-width groove — carries two coplanar top strips
    /// that share one Surface. Cell coplanar-narrow lets it chain a second boolean; OCCT
    /// scores the blind pocket cut into one strip on the exported b-rep.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slotted_bar_pocket_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let bar = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([3.0, 1.0, 1.0]),
        );
        let groove = m.add_cuboid(
            Point3::from_array([1.0, -0.5, 0.5]),
            Point3::from_array([2.0, 1.5, 1.5]),
        );
        let slotted = boolean_one(&mut m, BoolKind::Cut, bar, groove).unwrap();
        m.rebuild_adjacency();
        let pocket = m.add_cuboid(
            Point3::from_array([0.2, 0.2, 0.7]),
            Point3::from_array([0.5, 0.5, 1.5]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Cut, slotted, pocket).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, slotted, pocket).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "slotted bar pocket cut: {} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// A blind pocket cut into a face (cell coplanar-contact-cut): a prism inside the base with
    /// its top flush is carved out. OCCT scores the pocket and a transversal Cut chained onto
    /// the (non-convex) pocketed solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_cut_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let occt_pocket = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
        let pocketed = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, pocketed).unwrap().volume, occt_pocket.volume),
            "pocket cut: {} vs {}",
            mass_props(&m, pocketed).unwrap().volume,
            occt_pocket.volume
        );
        // Chain a transversal cut at a corner, away from the pocket and off its face planes
        // (the pocketed solid is non-convex → seam path).
        let cutter = m.add_cuboid(
            Point3::from_array([0.8, 0.8, 0.3]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, pocketed, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, pocketed, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "pocket then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A boss fuses onto a face inside its boundary (cell coplanar-contact-boss): the base's
    /// top gains the boss footprint as a hole and the boss rides on it. OCCT scores the fuse
    /// and a Cut chained onto the bossed solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn boss_fuse_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.25, 0.25, 1.0]),
            Point3::from_array([0.75, 0.75, 2.0]),
        );
        let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
        let bossed = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, bossed).unwrap().volume, occt_fuse.volume),
            "boss fuse: {} vs {}",
            mass_props(&m, bossed).unwrap().volume,
            occt_fuse.volume
        );
        // Chain a cut through the boss.
        let cutter = m.add_cuboid(
            Point3::from_array([0.4, 0.4, 0.5]),
            Point3::from_array([0.6, 0.6, 2.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, bossed, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, bossed, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "boss then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A boss overhangs a single edge of a face (cell coplanar-contact-overhang): part fuses,
    /// part cantilevers. OCCT scores the fuse and a transversal Cut chained onto the (non-convex)
    /// overhanging solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn overhang_fuse_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 1.0]),
            Point3::from_array([1.5, 0.75, 2.0]),
        );
        let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
        let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, overhung).unwrap().volume, occt_fuse.volume),
            "overhang fuse: {} vs {}",
            mass_props(&m, overhung).unwrap().volume,
            occt_fuse.volume
        );
        // Chain a transversal cut drilling straight through the cantilever (the overhung solid
        // is non-convex → seam path).
        let cutter = m.add_cuboid(
            Point3::from_array([1.1, 0.35, 0.5]),
            Point3::from_array([1.4, 0.65, 2.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, overhung, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, overhung, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "overhang then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A single-edge overhang Cut carves an edge-slot breaking out through a wall (cell
    /// coplanar-contact-overhang-cut). OCCT scores the slot and a transversal Cut chained onto
    /// the (non-convex) slotted solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn edge_slot_cut_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.25, 0.5]),
            Point3::from_array([1.5, 0.75, 1.0]),
        );
        let occt_slot = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
        let slotted = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, slotted).unwrap().volume, occt_slot.volume),
            "edge slot: {} vs {}",
            mass_props(&m, slotted).unwrap().volume,
            occt_slot.volume
        );
        // Chain a transversal cut drilling through the base away from the slot (the slotted
        // solid is non-convex → seam path).
        let cutter = m.add_cuboid(
            Point3::from_array([0.15, 0.8, -0.5]),
            Point3::from_array([0.35, 0.95, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slotted, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, slotted, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "slot then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A corner-overhanging boss fuses onto a face (cell coplanar-contact-overhang-corner): the
    /// boss footprint swallows a base-top corner, crossing two edges. OCCT scores the fuse and a
    /// transversal Cut chained onto the (non-convex, L-cantilever) result.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn corner_overhang_fuse_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let boss = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 1.0]),
            Point3::from_array([1.5, 1.5, 2.0]),
        );
        let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, boss).unwrap();
        let overhung = boolean_one(&mut m, BoolKind::Fuse, base, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, overhung).unwrap().volume, occt_fuse.volume),
            "corner fuse: {} vs {}",
            mass_props(&m, overhung).unwrap().volume,
            occt_fuse.volume
        );
        // Chain a transversal cut drilling through the L cantilever's outer corner (x>1, y>1).
        let cutter = m.add_cuboid(
            Point3::from_array([1.1, 1.1, 0.5]),
            Point3::from_array([1.4, 1.4, 2.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, overhung, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, overhung, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "corner then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A spanning-slab boss fuses onto a face (cell coplanar-contact-overhang-multi): the slab
    /// crosses the base top, splitting each contact face into two pieces. OCCT scores the fuse and
    /// a transversal Cut chained onto the (non-convex, two-cantilever) result.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn spanning_slab_fuse_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 0.4, 1.0]),
            Point3::from_array([1.5, 0.6, 2.0]),
        );
        let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, base, slab).unwrap();
        let slabbed = boolean_one(&mut m, BoolKind::Fuse, base, slab).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, slabbed).unwrap().volume, occt_fuse.volume),
            "slab fuse: {} vs {}",
            mass_props(&m, slabbed).unwrap().volume,
            occt_fuse.volume
        );
        // Chain a transversal cut drilling through the x>1 cantilever piece.
        let cutter = m.add_cuboid(
            Point3::from_array([1.1, 0.45, 0.5]),
            Point3::from_array([1.4, 0.55, 2.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slabbed, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, slabbed, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "slab then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A spanning-slab Cut carves a channel breaking out two opposite walls (cell
    /// coplanar-contact-overhang-slab-cut). OCCT scores the channel and a transversal Cut chained
    /// onto the (non-convex) channelled solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_channel_cut_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 0.4, 0.5]),
            Point3::from_array([1.5, 0.6, 1.0]),
        );
        let occt_channel = occt_boolean_of(&m, OcctBool::Cut, base, slab).unwrap();
        let channelled = boolean_one(&mut m, BoolKind::Cut, base, slab).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(
                mass_props(&m, channelled).unwrap().volume,
                occt_channel.volume
            ),
            "slab channel: {} vs {}",
            mass_props(&m, channelled).unwrap().volume,
            occt_channel.volume
        );
        // Chain a transversal drill through the base away from the channel (y > 0.6).
        let cutter = m.add_cuboid(
            Point3::from_array([0.2, 0.75, -0.5]),
            Point3::from_array([0.4, 0.95, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, channelled, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, channelled, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "channel then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A corner slot Cut breaks out two adjacent walls (cell coplanar-contact-overhang-corner-cut).
    /// OCCT scores the corner slot and a transversal Cut chained onto the (non-convex) result.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn corner_slot_cut_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let corner = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let occt_slot = occt_boolean_of(&m, OcctBool::Cut, base, corner).unwrap();
        let slotted = boolean_one(&mut m, BoolKind::Cut, base, corner).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, slotted).unwrap().volume, occt_slot.volume),
            "corner slot: {} vs {}",
            mass_props(&m, slotted).unwrap().volume,
            occt_slot.volume
        );
        // Chain a transversal drill through the base at the opposite (0,0) corner.
        let cutter = m.add_cuboid(
            Point3::from_array([0.05, 0.05, -0.5]),
            Point3::from_array([0.25, 0.25, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, slotted, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, slotted, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "corner then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// An L-step Cut breaks out three walls at once — two corners (x=0, x=1) and a fully covered
    /// wall (y=1) — via the general N-wall path (cell coplanar-contact-overhang-cut-general). OCCT
    /// scores the L rebate and a transversal Cut chained onto the (non-convex) result.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn l_step_cut_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([-0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let occt_step = occt_boolean_of(&m, OcctBool::Cut, base, prism).unwrap();
        let stepped = boolean_one(&mut m, BoolKind::Cut, base, prism).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, stepped).unwrap().volume, occt_step.volume),
            "l step: {} vs {}",
            mass_props(&m, stepped).unwrap().volume,
            occt_step.volume
        );
        // Chain a transversal drill through the full-height front strip (y < 0.5).
        let cutter = m.add_cuboid(
            Point3::from_array([0.1, 0.1, -0.5]),
            Point3::from_array([0.3, 0.3, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, stepped, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, stepped, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "l step then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// A `Common` of a top-flush overhang pair is the convex overlap R = a ∩ b (cell
    /// coplanar-contact-overhang-common). The prism hangs past one base edge; OCCT scores R.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn edge_overhang_common_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.3, 0.5, 0.5]),
            Point3::from_array([0.7, 1.5, 1.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Common, base, prism).unwrap();
        let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "edge overhang common: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A `Common` whose small box is contained in the base's top face (same-normal coplanar cap) yet
    /// pokes out the bottom (E0). OCCT cross-checks the contained-InterQ-island result (vol 0.25).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn contained_common_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let box_ = m.add_cuboid(
            Point3::from_array([0.25, 0.25, -0.5]),
            Point3::from_array([0.75, 0.75, 1.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Common, base, box_).unwrap();
        let r = boolean_one(&mut m, BoolKind::Common, base, box_).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "contained common: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A corner overhang `Common`: the prism swallows the base's (1,1) corner, so R meets at a
    /// corner column (cell coplanar-contact-overhang-common). OCCT cross-checks that cc topology.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn corner_overhang_common_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let prism = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Common, base, prism).unwrap();
        let r = boolean_one(&mut m, BoolKind::Common, base, prism).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "corner overhang common: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A blind pocket carved by a non-convex (L-shaped) cutter — the contained-coplanar Cut now
    /// admits non-convex operands. OCCT scores the L pocket.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn non_convex_profile_pocket_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 0.5]),
        );
        let l = Profile2d {
            points: [
                [-0.3, -0.3],
                [0.3, -0.3],
                [0.3, 0.0],
                [0.0, 0.0],
                [0.0, 0.3],
                [-0.3, 0.3],
            ]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
        };
        let OpOutput::Extrude { solid: lp, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: l,
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let occt = occt_boolean_of(&m, OcctBool::Cut, base, lp).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, base, lp).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "non-convex profile pocket: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A blind pocket carved into an already-pocketed (non-convex) cube. OCCT scores the second
    /// pocket cut on the concave kept solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_into_non_convex_solid_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d {
                    points: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
                        .iter()
                        .map(|&p| nacre_math::Point2::from_array(p))
                        .collect(),
                },
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid: pc, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: Profile2d {
                    points: [[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]
                        .iter()
                        .map(|&p| nacre_math::Point2::from_array(p))
                        .collect(),
                },
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let corner = m.add_cuboid(
            Point3::from_array([0.05, 0.1, 0.6]),
            Point3::from_array([0.25, 0.2, 1.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Cut, pc, corner).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, pc, corner).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "pocket into non-convex solid: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A boss raised by a non-convex (L-shaped) prism — the contained-coplanar Fuse now admits
    /// non-convex operands. The base's top sits at z = 0, the L boss extrudes onto it flush.
    /// OCCT scores the L boss.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn non_convex_profile_boss_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let mut m = Model::new();
        let base = m.add_cuboid(
            Point3::from_array([-1.0, -1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 0.0]),
        );
        let l = Profile2d {
            points: [
                [-0.3, -0.3],
                [0.3, -0.3],
                [0.3, 0.0],
                [0.0, 0.0],
                [0.0, 0.3],
                [-0.3, 0.3],
            ]
            .iter()
            .map(|&p| nacre_math::Point2::from_array(p))
            .collect(),
        };
        let OpOutput::Extrude { solid: lb, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: l,
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let occt = occt_boolean_of(&m, OcctBool::Fuse, base, lb).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, base, lb).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "non-convex profile boss: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A boss cantilevers off the side face of a top-pocketed cube (a non-convex solid). The
    /// overhang Fuse now admits a non-convex solid when the contact face is convex; OCCT scores the
    /// cantilever on the concave part.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn overhang_boss_on_non_convex_solid_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d {
                    points: vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)],
                },
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid: pc, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: Profile2d {
                    points: vec![p2(-0.2, -0.2), p2(0.2, -0.2), p2(0.2, 0.2), p2(-0.2, 0.2)],
                },
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        // Boss on the +x side face, overhanging the bottom edge.
        let boss = m.add_cuboid(
            Point3::from_array([1.0, 0.25, -0.25]),
            Point3::from_array([1.5, 0.75, 0.75]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, pc, boss).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, pc, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "overhang boss on non-convex solid: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// A boss raised on the top of a non-convex (L-prism) solid. OCCT scores the boss fused onto
    /// the concave base.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn boss_onto_non_convex_solid_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let mut m = Model::new();
        let OpOutput::Extrude { solid: l, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: Profile2d {
                    points: [
                        [0.0, 0.0],
                        [2.0, 0.0],
                        [2.0, 1.0],
                        [1.0, 1.0],
                        [1.0, 2.0],
                        [0.0, 2.0],
                    ]
                    .iter()
                    .map(|&p| nacre_math::Point2::from_array(p))
                    .collect(),
                },
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let boss = m.add_cuboid(
            Point3::from_array([0.3, 0.3, 1.0]),
            Point3::from_array([0.7, 0.7, 1.5]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, boss).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, boss).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt.volume),
            "boss onto non-convex solid: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt.volume
        );
    }

    /// Two face-to-face cubes fuse into a clean box (cell fuse-coplanar-merge merges the
    /// coplanar side faces and dissolves the interface corners), which then chains a Cut.
    /// OCCT scores both the fuse and the chained cut on the exported b-rep.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn stacked_fuse_then_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0, 1.0, 1.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        let occt_fuse = occt_boolean_of(&m, OcctBool::Fuse, a, b).unwrap();
        let stack = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
        m.rebuild_adjacency();
        assert!(
            approx(mass_props(&m, stack).unwrap().volume, occt_fuse.volume),
            "stacked fuse: {} vs {}",
            mass_props(&m, stack).unwrap().volume,
            occt_fuse.volume
        );
        // Chain a cut straddling the fused interface.
        let cutter = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 1.5, 1.5]),
        );
        let occt_cut = occt_boolean_of(&m, OcctBool::Cut, stack, cutter).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, stack, cutter).unwrap();
        assert!(
            approx(mass_props(&m, r).unwrap().volume, occt_cut.volume),
            "stacked fuse then cut: {} vs {}",
            mass_props(&m, r).unwrap().volume,
            occt_cut.volume
        );
    }

    /// The two convex pokes cell (5b) opened, against OCCT. Both operands are convex; the
    /// convex path rejected these as `poke_through` and cell (5b) routes them to the seam
    /// path. The notch is an edge crossed twice; the drill is a genus-1 solid.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn notch_cube_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let y = m.add_cuboid(
            Point3::from_array([3.0, -1.0, -1.0]),
            Point3::from_array([7.0, 1.4, 1.2]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Cut, a, y).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, a, y).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "vol {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn notch_cube_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let y = m.add_cuboid(
            Point3::from_array([3.0, -1.0, -1.0]),
            Point3::from_array([7.0, 1.4, 1.2]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, y).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, a, y).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "vol {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn drilled_cube_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Cut, a, bar).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, a, bar).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "vol {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn drilled_cube_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
        let bar = m.add_cuboid(
            Point3::from_array([1.0, 1.0, -1.0]),
            Point3::from_array([2.0, 2.0, 4.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, a, bar).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, a, bar).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "vol {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The reflex-corner bite (M5-d3 cell 3d). The arc's single bend projects *outside*
    /// its chord's endpoints, so the old projection sort was ordering it by luck. Fills
    /// the OCCT gap that cell left.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn reflex_bite_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_reflex_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_popup_box();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_popup_box();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_dimple();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bx).unwrap();
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
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_dimple();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bx).unwrap();
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

    /// The L-prism and an L-shaped bar lying in its notch, `z ∈ [0.5, 1.5]`, its two arm
    /// ends biting the cap's convex corners. Two chords on one face (M5-d3 cell 3e-2).
    fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let raised = SketchPlane {
            origin: Point3::from_array([0.0, 0.0, 0.5]),
            ..SketchPlane::world_xy()
        };
        let bar = extrude(
            &mut m,
            raised,
            &[
                [1.8, 0.8],
                [2.1, 0.8],
                [2.1, 2.1],
                [0.8, 2.1],
                [0.8, 1.8],
                [1.8, 1.8],
            ],
        );
        (m, l, bar)
    }

    /// The L-prism and a П-shaped staple drawn in the XZ plane and extruded along `−y`, so
    /// the L's cap is parallel to the extrusion axis and the staple's section there falls in
    /// two: a loop wholly inside the cap, and an arc wrapping the reflex corner (cell 3f-4).
    fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
        use nacre_math::Vector3;
        let (mut m, l) = l_prism();
        let xz = SketchPlane {
            origin: Point3::from_array([0.0, 1.3, 0.0]),
            x_axis: Vector3::from_array([1.0, 0.0, 0.0]),
            y_axis: Vector3::from_array([0.0, 0.0, 1.0]),
        };
        let st = extrude_dist(
            &mut m,
            xz,
            &[
                [0.1, 0.5],
                [0.6, 0.5],
                [0.6, 1.3],
                [0.8, 1.3],
                [0.8, 0.45],
                [1.4, 0.45],
                [1.4, 1.5],
                [0.1, 1.5],
            ],
            0.65,
        );
        (m, l, st)
    }

    /// The same loop is a **hole** here — inside the cap's kept region — and an **island**
    /// under `Cut(staple, L)`, where the kept region is the corner bite alone. Containment
    /// alone tells them apart; the three volumes close inclusion–exclusion on `V_∩ = 0.311`.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn l_staple_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, st) = l_and_staple();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, st).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, st).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "l staple cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "l staple cut area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn l_staple_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, st) = l_and_staple();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, st).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, st).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "l staple fuse volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "l staple fuse area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The island case, scored by a kernel that never heard of our containment test.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn staple_cut_by_l_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, st) = l_and_staple();
        let occt = occt_boolean_of(&m, OcctBool::Cut, st, l).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, st, l).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "staple cut by l volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "staple cut by l area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The U-prism (volume 5.3) and a slab shearing off both prong tops. The slab's
    /// `y = 1.5` face carries **two** closed loops — one per prong — wholly inside it.
    fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let u = extrude(
            &mut m,
            SketchPlane::world_xy(),
            &[
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 2.3],
                [2.0, 2.3],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
        );
        let slab = m.add_cuboid(
            Point3::from_array([-0.5, 1.5, -0.5]),
            Point3::from_array([3.5, 2.5, 1.5]),
        );
        (m, u, slab)
    }

    /// Three diffs on one fixture, because the three reconstruct different things: two
    /// **islands** from one face (`Cut(u, slab)`), the same face's two **holes** on B
    /// (`Fuse`), and two holes on A with the U's caps split into cycles (`Cut(slab, u)`).
    /// Cell 3f-3.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn u_slab_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, u, slab) = u_and_slab();
        let occt = occt_boolean_of(&m, OcctBool::Cut, u, slab).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, u, slab).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "u slab cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "u slab cut area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn u_slab_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, u, slab) = u_and_slab();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, u, slab).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, u, slab).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "u slab fuse volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "u slab fuse area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_cut_by_u_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, u, slab) = u_and_slab();
        let occt = occt_boolean_of(&m, OcctBool::Cut, slab, u).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, slab, u).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "slab cut by u volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "slab cut by u area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The L-prism with an L-shaped stub standing wholly inside its cap. The blind pocket's
    /// lid carries the suite's first **non-convex** inner loop (M5-d3 cell 3h).
    fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
        let (mut m, l) = l_prism();
        let raised = SketchPlane {
            origin: Point3::from_array([0.0, 0.0, 0.5]),
            ..SketchPlane::world_xy()
        };
        let stub = extrude(
            &mut m,
            raised,
            &[
                [0.2, 0.25],
                [0.85, 0.25],
                [0.85, 0.4],
                [0.35, 0.4],
                [0.35, 0.9],
                [0.2, 0.9],
            ],
        );
        (m, l, stub)
    }

    /// An independent kernel on a hole that is not a rectangle. Volume scores the pocket;
    /// area scores its six walls, since the lid gives up exactly what the floor hands back.
    /// Neither can see the loop's *winding* — both kernels sum unsigned areas — which is why
    /// cell 3h added a second source for that and did not lean on this diff.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn ell_dimple_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, stub) = l_and_ell_stub();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, stub).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, stub).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "ell dimple cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "ell dimple cut area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// nacre's first face carrying more than one chord (M5-d3 cell 3e-2) vs OCCT, `Cut`.
    /// The L's cap keeps a single ring that uses both arcs, while the bar's floor splits
    /// into two faces — the two bites' floors. Volume scores the bites; area cannot, since
    /// a corner cut hands back exactly the faces it removes (14.0 either way).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn notch_bar_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, bar) = l_and_notch_bar();
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, bar).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, bar).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "notch bar cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "notch bar cut area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The `Fuse` counterpart, and not a symmetry re-ask: it reconstructs different faces.
    /// Keeping both outsides, the bar's floor stays one ring spanning both arcs instead of
    /// splitting, and the area — `14 + 6.58 − 0.96` — finally has something to score.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn notch_bar_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, bar) = l_and_notch_bar();
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, bar).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, bar).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "notch bar fuse volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "notch bar fuse area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// nacre's first face whose outer loop is *all* seam (M5-d3 cell 3f-2) vs OCCT: the
    /// stub cut by the L, leaving the `0.4 × 0.4 × 0.5` box above `z = 1`. Its floor is
    /// that island face.
    ///
    /// Both kernels take a face's normal from their own bookkeeping and sum unsigned
    /// areas, so neither this diff nor any other can see the island wound backwards —
    /// `validate` and the signed mesh volume do that. What an independent kernel scores
    /// here is the face's *existence and extent*: get the ring's nodes wrong and the
    /// polygon's area moves with it.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn island_cut_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, bx) = l_and_dimple();
        let occt = occt_boolean_of(&m, OcctBool::Cut, bx, l).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, bx, l).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "island cut volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "island cut area: {} vs {}",
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
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
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
        let r = boolean_one(&mut m, BoolKind::Cut, solid, bx).unwrap();
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
        use nacre_ops::BoolKind;
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
        let r = boolean_one(&mut m, BoolKind::Fuse, a, b).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{} vs {}",
            nacre.volume,
            occt.volume
        );
    }

    /// The unit cube with a `0.4`-square pocket `0.5` deep in its top face — the fixture
    /// whose lid carries an inner loop. OCCT reads the same solid from STEP; nothing here
    /// asks it to reproduce `PocketOnFace`.
    fn pocketed_cube() -> (Model, Handle<Solid>) {
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};
        let prof = |pts: &[[f64; 2]]| Profile2d {
            points: pts
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: prof(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: prof(&[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]),
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        (m, solid)
    }

    /// Score one boolean of a holed operand against OCCT, on volume and on area. Area
    /// matters here: nacre's own gates all read the same rings, so only an independent
    /// kernel makes the surviving hole's size a real claim.
    fn diff_holed(name: &str, kind: OcctBool, boxes: [[f64; 3]; 2], swap: bool) {
        use nacre_ops::BoolKind;
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(Point3::from_array(boxes[0]), Point3::from_array(boxes[1]));
        let (x, y) = if swap { (bx, pc) } else { (pc, bx) };
        let occt = occt_boolean_of(&m, kind, x, y).unwrap();
        let bk = match kind {
            OcctBool::Cut => BoolKind::Cut,
            OcctBool::Fuse => BoolKind::Fuse,
            OcctBool::Common => BoolKind::Common,
        };
        let r = boolean_one(&mut m, bk, x, y).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "{name} volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "{name} area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// A boolean composing on a shape a boolean can make (M5-d3 cell 3f-5). The seam bites
    /// the holed lid's corner, and the hole is placed inside the region left behind.
    ///
    /// The box is the symmetric one cell (5a) gave back: its vertical edge pierces the lid
    /// at `(0.85, 0.85)`, on a fan diagonal from every apex. OCCT is asked the same question
    /// on the same coordinates, and it does not fan.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_corner_cut_matches_occt() {
        diff_holed(
            "pocket corner cut",
            OcctBool::Cut,
            [[0.85, 0.85, 0.85], [1.15, 1.15, 1.15]],
            false,
        );
    }

    /// The same bite mirrored in `z`: the seam misses the lid, which rides out whole.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_bottom_corner_cut_matches_occt() {
        diff_holed(
            "pocket bottom corner cut",
            OcctBool::Cut,
            [[0.85, 0.85, -0.15], [1.15, 1.15, 0.15]],
            false,
        );
    }

    /// A corner box whose footprint overlaps the pocket, so its walls cross the lid's hole
    /// rim (M5-d3 cell 3f-6). `∂(lid)` is two rings the seam threads into one notch, opening
    /// the pocket to the outside. Only an independent kernel makes the absorbed hole's area a
    /// real claim; OCCT is asked on the same coordinates. `Cut(pc, box) = 0.893`.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_rim_corner_cut_matches_occt() {
        diff_holed(
            "pocket rim corner cut",
            OcctBool::Cut,
            [[0.55, 0.55, 0.85], [1.15, 1.15, 1.15]],
            false,
        );
    }

    /// The same two solids named the other way — the box's face is arranged first and the
    /// pocket's rim pierces it. `Cut(box, pc) = 0.081`.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_rim_corner_cut_either_way_matches_occt() {
        diff_holed(
            "pocket rim corner cut either way",
            OcctBool::Cut,
            [[0.55, 0.55, 0.85], [1.15, 1.15, 1.15]],
            true,
        );
    }

    /// A slab through the pocket between its floor and its lid (M5-d3 cell 3f-7). Its underside
    /// carries the cube's cross-section as one loop with the pocket's nested inside it — a loop
    /// within a loop. `Cut(slab, pc) = 1.488`, the pocket loop hung in the slab region as a hole
    /// beside an island; the area is the independent claim that the nesting placed both rings.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_nests_pocket_cut_matches_occt() {
        diff_holed(
            "slab nests pocket cut",
            OcctBool::Cut,
            [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
            false,
        );
    }

    /// The other order: `Cut(pc, slab) = 0.668`, where the slab's dropped underside leaves no
    /// region, so the cube section becomes an island carrying the pocket loop as its own hole.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_nests_pocket_cut_either_way_matches_occt() {
        diff_holed(
            "slab nests pocket cut either way",
            OcctBool::Cut,
            [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
            true,
        );
    }

    /// `Fuse` of the same two seals the pocket into an enclosed cavity (cell 5c): material
    /// `2.408`, a void of `0.032`. OCCT builds the same hollow solid and its volume subtracts
    /// the void while its area adds both surfaces — the independent claim that the seam path
    /// assembled the void inward, not as a phantom outer piece.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_seals_pocket_fuse_matches_occt() {
        diff_holed(
            "slab seals pocket fuse",
            OcctBool::Fuse,
            [[-0.2, -0.25, 0.7], [1.3, 1.2, 1.5]],
            false,
        );
    }

    /// A slab over the pocket, its underside below the pocket floor. `Cut` keeps the lid
    /// as a reversed inside-B piece, hole and all — `flip` meeting `inner` for the first
    /// time. `Fuse` drops it. The third scores the complement.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_cut_by_pocket_matches_occt() {
        diff_holed(
            "slab cut by pocket",
            OcctBool::Cut,
            [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
            true,
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn slab_and_pocket_fuse_matches_occt() {
        diff_holed(
            "slab and pocket fuse",
            OcctBool::Fuse,
            [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
            true,
        );
    }

    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn pocket_cut_by_slab_matches_occt() {
        diff_holed(
            "pocket cut by slab",
            OcctBool::Cut,
            [[-0.2, -0.25, 0.3], [1.3, 1.2, 1.5]],
            false,
        );
    }

    /// `coincident_merge` with a holed operand. The imprint splits the `x = 1` face into a
    /// holed remainder and a coplanar region face — the one holed shape the convex and seam
    /// paths never see, since `has_coplanar_pair` stops it at their door. OCCT reads the
    /// STEP as one cube either way, so this scores the merge, not the imprint.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn imprinted_merge_matches_occt() {
        use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, apply};
        let prof = |pts: &[[f64; 2]]| Profile2d {
            points: pts
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: prof(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let OpOutput::ImprintSketch { solid, .. } = apply(
            &mut m,
            &Operation::ImprintSketch {
                face: faces[3],
                profile: prof(&[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]),
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let bx = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 1.0]),
            Point3::from_array([1.0, 1.0, 2.0]),
        );
        let occt = occt_boolean_of(&m, OcctBool::Fuse, solid, bx).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, solid, bx).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "imprinted merge volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "imprinted merge area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// A rod drilled clean through the L's bar: the first genus-1 solid this kernel makes
    /// (M5-d3 cell 3e-3). Area is scored alongside volume — a tunnel's walls are area, and
    /// nacre's own gates all read the same rings, so only an independent kernel makes the
    /// hole's size a claim rather than a restatement.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn drilled_l_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, rod) = l_prism_and_box([0.3, 0.3, -0.5], [0.5, 0.6, 1.5]);
        let occt = occt_boolean_of(&m, OcctBool::Cut, l, rod).unwrap();
        let r = boolean_one(&mut m, BoolKind::Cut, l, rod).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "drilled L volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "drilled L area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The `Fuse`: the rod stands proud on both faces of the bar, and each of its walls
    /// splits into the stub above and the stub below.
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn l_and_rod_fuse_matches_occt() {
        use nacre_ops::BoolKind;
        let (mut m, l, rod) = l_prism_and_box([0.3, 0.3, -0.5], [0.5, 0.6, 1.5]);
        let occt = occt_boolean_of(&m, OcctBool::Fuse, l, rod).unwrap();
        let r = boolean_one(&mut m, BoolKind::Fuse, l, rod).unwrap();
        let nacre = mass_props(&m, r).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "l and rod fuse volume: {} vs {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "l and rod fuse area: {} vs {}",
            nacre.area,
            occt.area
        );
    }

    /// The non-convex overhang Cut deliverable (unified coplanar handler): a slot cut from a
    /// top-pocketed cube, flush on the +x wall, breaking out the bottom, its top coplanar-disjoint
    /// with the pocket floor. OCCT confirms the accepted volume (hand estimate 0.8575).
    #[test]
    #[ignore = "requires OCCT DRAWEXE (run with --ignored)"]
    fn non_convex_overhang_cut_matches_occt() {
        use nacre_ops::{OpOutput, Operation, Profile2d, apply};
        let sq = |pts: &[[f64; 2]]| Profile2d {
            points: pts
                .iter()
                .map(|&p| nacre_math::Point2::from_array(p))
                .collect(),
        };
        let mut m = Model::new();
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                plane: SketchPlane::world_xy(),
                profile: sq(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
                dist: 1.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let OpOutput::PocketOnFace { solid: pc, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: sq(&[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]),
                dist: 0.5,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let slot = m.add_cuboid(
            Point3::from_array([0.75, 0.25, -0.25]),
            Point3::from_array([1.0, 0.75, 0.5]),
        );
        // OCCT ground truth from the inputs, before nacre supersedes them.
        let occt = occt_boolean_of(&m, OcctBool::Cut, pc, slot).unwrap();
        let solids = boolean(&mut m, BoolKind::Cut, pc, slot).unwrap();
        assert_eq!(solids.len(), 1);
        let nacre = mass_props(&m, solids[0]).unwrap();
        assert!(
            approx(nacre.volume, occt.volume),
            "volume {} vs occt {}",
            nacre.volume,
            occt.volume
        );
        assert!(
            approx(nacre.volume, 0.8575),
            "volume {} vs hand 0.8575",
            nacre.volume
        );
        assert!(
            approx(nacre.area, occt.area),
            "area {} vs occt {}",
            nacre.area,
            occt.area
        );
    }
}
