use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn binary_path() -> &'static str {
    env!("CARGO_BIN_EXE_mzio")
}

fn demo_mzml_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("assets/demo.mzML")
}

fn unique_tsv_path(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzml-spinoff-{stem}-{}-{timestamp}.tsv",
        std::process::id()
    ))
}

#[test]
fn root_help_mentions_diagnostics() {
    let output = Command::new(binary_path())
        .arg("--help")
        .output()
        .expect("run root --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("diagnostics"));
    assert!(stdout.contains("diag"));
}

#[test]
fn diagnostics_help_mentions_scanning_flags() {
    let output = Command::new(binary_path())
        .args(["diagnostics", "--help"])
        .output()
        .expect("run diagnostics --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--panel <name>"));
    assert!(stdout.contains("--target <spec>"));
    assert!(stdout.contains("--delta <spec>"));
    assert!(stdout.contains("--top-n <n>"));
    assert!(stdout.contains("--max-delta-hits <n>"));
}

#[test]
fn diag_alias_writes_tsv_with_header() {
    let tsv_path = unique_tsv_path("diagnostics");
    let output = Command::new(binary_path())
        .args([
            "diag",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--target",
            "hexnac:204.08665",
            "--tsv",
            tsv_path.to_str().expect("tsv path"),
        ])
        .output()
        .expect("run diag alias");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let tsv = fs::read_to_string(&tsv_path).expect("read diagnostics tsv");
    assert!(tsv.starts_with(
        "file_stem\tscan_index\tscan_number\tscan_id\tms_level\trt_min\tprecursor_mz\tprecursor_charge\tkind\tpanel\tlabel\ttheoretical"
    ));

    let _ = fs::remove_file(tsv_path);
}

#[test]
fn diagnostics_rejects_conflicting_tolerance_flags() {
    let output = Command::new(binary_path())
        .args([
            "diagnostics",
            "--mzml",
            "x.mzML",
            "--target",
            "hexnac:204.08665",
            "--tol-ppm",
            "20",
            "--tol-da",
            "0.05",
        ])
        .output()
        .expect("run diagnostics with conflicting tolerances");
    assert!(!output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("only one of --tol-ppm or --tol-da"));
}
