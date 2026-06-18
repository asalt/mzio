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

fn unique_svg_path(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzio-plot-survey-{stem}-{}-{timestamp}.svg",
        std::process::id()
    ))
}

fn unique_out_dir(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzio-plot-survey-{stem}-{}-{timestamp}",
        std::process::id()
    ))
}

#[test]
fn plot_survey_help_mentions_views() {
    let output = Command::new(binary_path())
        .args(["plot-survey", "--help"])
        .output()
        .expect("run plot-survey --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--view <chrom|ms1-map|dda|default|all>"));
    assert!(stdout.contains("--chrom"));
    assert!(stdout.contains("--map"));
    assert!(stdout.contains("--out-dir <dir>"));
    assert!(stdout.contains("--prefix <text>"));
    assert!(stdout.contains("--svg <path>"));
    assert!(stdout.contains("--png [path]"));
}

#[test]
fn plot_survey_writes_default_svgs_to_out_dir() {
    let out_dir = unique_out_dir("default");
    let output = Command::new(binary_path())
        .args([
            "plot-survey",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--out-dir",
            out_dir.to_str().expect("out dir"),
            "--prefix",
            "demo_run",
        ])
        .output()
        .expect("run plot-survey");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let chrom = out_dir.join("demo_run_survey_chrom.svg");
    let map = out_dir.join("demo_run_survey_ms1_map.svg");
    let chrom_svg = fs::read_to_string(&chrom).expect("read chrom survey svg");
    assert!(chrom_svg.contains("<svg"));
    assert!(chrom_svg.contains("Run chromatograms"));
    assert!(chrom_svg.contains("Source: demo.mzML"));
    assert!(chrom_svg.contains("TIC (MS1)"));
    assert!(chrom_svg.contains("BPC (MS1)"));

    let map_svg = fs::read_to_string(&map).expect("read map survey svg");
    assert!(map_svg.contains("<svg"));
    assert!(map_svg.contains("MS1 base peak map"));
    assert!(map_svg.contains("Source: demo.mzML"));
    assert!(map_svg.contains("MS1 Base Peak Map"));

    let _ = fs::remove_dir_all(out_dir);
}

#[test]
fn plot_survey_writes_single_chrom_svg_when_requested() {
    let svg_path = unique_svg_path("chrom");
    let output = Command::new(binary_path())
        .args([
            "plot-survey",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--view",
            "chrom",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run plot-survey chrom");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let svg = fs::read_to_string(&svg_path).expect("read chrom svg");
    assert!(svg.contains("Run chromatograms"));
    assert!(svg.contains("TIC (MS1)"));
    assert!(svg.contains("BPC (MS1)"));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_survey_writes_dda_svgs() {
    let out_dir = unique_out_dir("dda");
    let output = Command::new(binary_path())
        .args([
            "plot-survey",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--view",
            "dda",
            "--out-dir",
            out_dir.to_str().expect("out dir"),
            "--prefix",
            "demo_run",
        ])
        .output()
        .expect("run plot-survey dda");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let density =
        fs::read_to_string(out_dir.join("demo_run_survey_dda_density.svg")).expect("read density");
    assert!(density.contains("DDA acquisition density"));
    assert!(density.contains("MS2 scans per RT bin"));
    assert!(density.contains("Cycle time / MS2 per MS1 cycle"));

    let precursors = fs::read_to_string(out_dir.join("demo_run_survey_dda_precursors.svg"))
        .expect("read precursors");
    assert!(precursors.contains("DDA precursor survey"));
    assert!(precursors.contains("MS2 precursor m/z map"));
    assert!(precursors.contains("MS2 charge distribution"));
    assert!(precursors.contains("charge"));

    let _ = fs::remove_dir_all(out_dir);
}

#[test]
fn plot_survey_rejects_svg_for_multiple_outputs() {
    let svg_path = unique_svg_path("ambiguous");
    let output = Command::new(binary_path())
        .args([
            "plot-survey",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run plot-survey ambiguous svg");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--svg requires exactly one output"));
}

#[test]
fn plot_survey_rejects_unknown_view() {
    let output = Command::new(binary_path())
        .args([
            "plot-survey",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--view",
            "nope",
        ])
        .output()
        .expect("run plot-survey bad view");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("unknown survey view"));
}
