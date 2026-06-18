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

fn unique_output_path(stem: &str, extension: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzio-extract-{stem}-{}-{timestamp}.{extension}",
        std::process::id()
    ))
}

#[test]
fn extract_help_mentions_multiple_output_formats() {
    let output = Command::new(binary_path())
        .args(["extract", "--help"])
        .output()
        .expect("run extract --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--mzml <file>"));
    assert!(stdout.contains("--scan <n>"));
    assert!(stdout.contains("--format <tsv|ms2>"));
    assert!(stdout.contains("--ms2 <file>"));
    assert!(stdout.contains("`.tsv` or `.ms2` is inferred"));
}

#[test]
fn scan_help_mentions_extract_alias() {
    let output = Command::new(binary_path())
        .args(["scan", "--help"])
        .output()
        .expect("run scan --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Alias of"));
    assert!(stdout.contains("extract"));
}

#[test]
fn extract_exports_selected_scan_to_stdout_as_tsv_by_default() {
    let output = Command::new(binary_path())
        .args([
            "extract",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--scan",
            "2",
        ])
        .output()
        .expect("run extract export");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("# scan_number\t2"));
    assert!(stdout.contains("# ms_level\t2"));
    assert!(stdout.contains("# precursor_mz\t500.200000"));
    assert!(stdout.contains("mz\tintensity"));
    assert!(stdout.contains("150.000000\t5.000000"));
}

#[test]
fn scan_alias_writes_tsv_export_to_file() {
    let output_path = unique_output_path("scan-alias", "tsv");
    let output = Command::new(binary_path())
        .args([
            "scan",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--scan",
            "scan=2",
            "--output",
            output_path.to_str().expect("output path"),
        ])
        .output()
        .expect("run scan export to file");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Wrote TSV export:"));

    let written = fs::read_to_string(&output_path).expect("read scan export");
    assert!(written.contains("# scan_id\tscan=2"));
    assert!(written.contains("350.000000\t80.000000"));

    let _ = fs::remove_file(output_path);
}

#[test]
fn extract_infers_ms2_output_from_extension() {
    let output_path = unique_output_path("ms2-extension", "ms2");
    let output = Command::new(binary_path())
        .args([
            "extract",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--scan",
            "2",
            "--output",
            output_path.to_str().expect("output path"),
        ])
        .output()
        .expect("run extract ms2 export");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Wrote MS2 export:"));

    let written = fs::read_to_string(&output_path).expect("read ms2 export");
    assert!(written.contains("H\tExtractor\tmzio"));
    assert!(written.contains("S\t2\t2\t500.200000"));
    assert!(written.contains("Z\t2\t998.385447"));
    assert!(written.contains("150.000000\t5.000000"));

    let _ = fs::remove_file(output_path);
}

#[test]
fn extract_rejects_conflicting_format_and_extension() {
    let output_path = unique_output_path("conflict", "ms2");
    let output = Command::new(binary_path())
        .args([
            "extract",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--scan",
            "2",
            "--output",
            output_path.to_str().expect("output path"),
            "--format",
            "tsv",
        ])
        .output()
        .expect("run conflicting extract export");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("suggests `ms2`"));

    let _ = fs::remove_file(output_path);
}
