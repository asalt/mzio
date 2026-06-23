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

fn demo_ms2_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("assets/demo.ms2")
}

fn unique_svg_path(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzml-spinoff-{stem}-{}-{timestamp}.svg",
        std::process::id()
    ))
}

#[test]
fn plot_help_mentions_annotation_flags() {
    let output = Command::new(binary_path())
        .args(["plot", "--help"])
        .output()
        .expect("run plot --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--scan <n>"));
    assert!(stdout.contains("--ms2 <file>"));
    assert!(stdout.contains("defaults to first spectrum"));
    assert!(stdout.contains("--peptide <SEQ>"));
    assert!(stdout.contains("M[+15.9949]"));
    assert!(stdout.contains("--svg-prefix <text>"));
    assert!(stdout.contains("--mod <position>:<delta>"));
    assert!(stdout.contains("--neutral-losses"));
    assert!(stdout.contains("--neutral-loss-min-frac <f>"));
    assert!(stdout.contains("--isotope-errors <list>"));
    assert!(stdout.contains("--remove-precursor"));
    assert!(stdout.contains("--tol-da <da>"));
}

#[test]
fn plot_rejects_conflicting_tolerance_flags() {
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            "x.mzML",
            "--index",
            "1",
            "--peptide",
            "PEPTIDE",
            "--tol-ppm",
            "20",
            "--tol-da",
            "0.5",
        ])
        .output()
        .expect("run plot with conflicting tolerances");
    assert!(!output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("only one of --tol-ppm or --tol-da"));
}

#[test]
fn plot_requires_selector_for_mzml_inputs() {
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
        ])
        .output()
        .expect("run plot without mzml selector");
    assert!(!output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("for mzML input"));
}

#[test]
fn plot_writes_annotated_svg() {
    let svg_path = unique_svg_path("annotated");
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--peptide",
            "PEPTIDEK/2",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run annotated plot");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Quality: SNR="));
    assert!(stdout.contains("log2_SNR="));
    assert!(stdout.contains("frag_error_mae_ppm="));

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("Matched fragment ladder"));
    assert!(svg.contains("Source: demo.mzML"));
    assert!(svg.contains("Peptide: PEPTIDEK"));
    assert!(svg.contains("tolerance: 20.0 ppm"));
    assert!(svg.contains("Quality: SNR="));
    assert!(svg.contains("log2_SNR="));
    assert!(svg.contains("frag_error_mae_ppm="));
    assert!(svg.contains("class=\"ladder-index\""));
    assert!(svg.contains("Full ion table"));
    assert!(svg.contains("matched colored, missing grey; p/w/n = H3PO4/H2O/NH3"));
    assert!(!svg.contains(">b++<"));
    assert!(!svg.contains(">y++<"));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_reports_fragment_error_in_da_for_wide_da_tolerance() {
    let svg_path = unique_svg_path("annotated-da-error");
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--peptide",
            "PEPTIDEK/2",
            "--tol-da",
            "0.5",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run annotated plot with Da tolerance");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("frag_error_mae_da="));
    assert!(!stdout.contains("frag_error_mae_ppm="));

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("frag_error_mae_da="));
    assert!(!svg.contains("frag_error_mae_ppm="));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_accepts_inline_mass_shift_syntax() {
    let svg_path = unique_svg_path("inline-mods");
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--peptide",
            "[+42.0106]PEPM[+15.9949]IDEK/2",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run inline-mod plot");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("Peptide: [+42.010600]PEPM[+15.994900]IDEK"));
    assert!(svg.contains("Applied mods: n-term:+42.010600, 4:+15.994900"));
    assert!(svg.contains("Precursor isotope errors: 0,1,2"));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_ms2_defaults_to_first_spectrum() {
    let svg_path = unique_svg_path("ms2-default");
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--ms2",
            demo_ms2_path().to_str().expect("demo ms2 path"),
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run ms2 plot without selector");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Scan scan=101 (index 0) | ms2"));

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("Source: demo.ms2"));
    assert!(svg.contains("Scan scan=101 | index 0 | ms2"));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_ms2_accepts_explicit_scan_selection() {
    let svg_path = unique_svg_path("ms2-scan");
    let output = Command::new(binary_path())
        .args([
            "plot",
            "--ms2",
            demo_ms2_path().to_str().expect("demo ms2 path"),
            "--scan",
            "205",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run ms2 plot with scan selector");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Scan scan=205 (index 1) | ms2"));

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("Scan scan=205 | index 1 | ms2"));
    assert!(svg.contains("Source: demo.ms2"));

    let _ = fs::remove_file(svg_path);
}

#[test]
fn plot_remove_precursor_hides_guide_marker_in_svg() {
    let visible_svg_path = unique_svg_path("precursor-visible");
    let hidden_svg_path = unique_svg_path("precursor-hidden");

    let visible_output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--mz-max",
            "600",
            "--svg",
            visible_svg_path.to_str().expect("visible svg path"),
        ])
        .output()
        .expect("run default precursor plot");

    assert!(
        visible_output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&visible_output.stdout),
        String::from_utf8_lossy(&visible_output.stderr),
    );

    let hidden_output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--mz-max",
            "600",
            "--remove-precursor",
            "--svg",
            hidden_svg_path.to_str().expect("hidden svg path"),
        ])
        .output()
        .expect("run remove-precursor plot");

    assert!(
        hidden_output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&hidden_output.stdout),
        String::from_utf8_lossy(&hidden_output.stderr),
    );

    let visible_svg = fs::read_to_string(&visible_svg_path).expect("read visible svg");
    let hidden_svg = fs::read_to_string(&hidden_svg_path).expect("read hidden svg");

    assert!(visible_svg.contains(">precursor "));
    assert!(!hidden_svg.contains(">precursor "));
    assert!(hidden_svg.contains("Precursor:"));

    let _ = fs::remove_file(visible_svg_path);
    let _ = fs::remove_file(hidden_svg_path);
}
