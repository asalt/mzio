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

fn local_slim_bruker_dia_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data/timstof2_hela_dia_slice/bruker_d/hela_dia_iildlisespik_rt21_60_21_78.d")
}

fn unique_prefix(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzio-dia-slice-{stem}-{}-{timestamp}",
        std::process::id()
    ))
}

fn output_path(prefix: &PathBuf, suffix: &str) -> PathBuf {
    let file_name = prefix
        .file_name()
        .and_then(|value| value.to_str())
        .expect("prefix file name");
    prefix
        .parent()
        .map(|parent| parent.join(format!("{file_name}.{suffix}")))
        .unwrap_or_else(|| PathBuf::from(format!("{file_name}.{suffix}")))
}

#[test]
fn main_help_mentions_dia_slice() {
    let output = Command::new(binary_path())
        .arg("--help")
        .output()
        .expect("run --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("dia-slice"));
    assert!(stdout.contains("Bruker .d"));
}

#[test]
fn dia_slice_help_mentions_runtime_and_slice_flags() {
    let output = Command::new(binary_path())
        .args(["dia-slice", "--help"])
        .output()
        .expect("run dia-slice --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--bruker <run.d>"));
    assert!(stdout.contains("--mzml <file>"));
    assert!(stdout.contains("--bruker-so <path>"));
    assert!(stdout.contains("--bruker-backend <name>"));
    assert!(stdout.contains("--python <exe>"));
    assert!(stdout.contains("--peptide <SEQ>"));
    assert!(stdout.contains("--fragment <ion>"));
    assert!(stdout.contains("--target <label:mz>"));
    assert!(stdout.contains("--neutral-losses"));
    assert!(stdout.contains("--charge <int>"));
    assert!(stdout.contains("--pseudo-ms2"));
    assert!(stdout.contains("--pseudo-ms2-rt-window <min>"));
    assert!(stdout.contains("--trace-peaks"));
    assert!(stdout.contains("--emit-trace"));
    assert!(stdout.contains("--trace-peak-smooth-window <n>"));
    assert!(stdout.contains("--trace-peak-snr <x>"));
    assert!(stdout.contains("--trace-peak-boundary-snr <x>"));
    assert!(stdout.contains("--trace-peak-boundary-fraction <x>"));
    assert!(stdout.contains("--trace-peak-min-points <n>"));
    assert!(stdout.contains("--trace-peak-min-nonzero <n>"));
    assert!(stdout.contains("--mz-da <da>"));
    assert!(stdout.contains("all RT"));
    assert!(stdout.contains("all mobility"));
    assert!(stdout.contains("all windows"));
    assert!(stdout.contains("--quad-min <mz>"));
    assert!(stdout.contains("--outdir <dir>"));
    assert!(stdout.contains("--rt-smooth"));
    assert!(stdout.contains("--rt-smooth-window <n>"));
    assert!(stdout.contains("run_tic.json"));
    assert!(stdout.contains("peaks.tsv"));
    assert!(stdout.contains("trace.tsv"));
    assert!(stdout.contains("-v, --verbose"));
    assert!(stdout.contains("-q, --quiet"));
    assert!(stdout.contains("MZIO_BRUKER_SO"));
}

#[test]
fn dia_slice_mzml_writes_profiles_and_svg() {
    let outdir = unique_prefix("mzml");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--mz",
            "500",
            "--mz-da",
            "500",
            "--outdir",
            outdir.to_str().expect("outdir"),
        ])
        .output()
        .expect("run dia-slice on mzml");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Wrote DIA slice outputs"));

    let out_prefix = outdir.join("demo.dia_slice_mz_500.0000");
    let summary =
        fs::read_to_string(output_path(&out_prefix, "summary.txt")).expect("read summary output");
    assert!(summary.contains("backend\tmzML"));
    assert!(summary.contains("spectra_considered"));

    let rt_profile =
        fs::read_to_string(output_path(&out_prefix, "rt_profile.tsv")).expect("read rt profile");
    assert!(rt_profile.contains("scan_index"));

    let mz_profile =
        fs::read_to_string(output_path(&out_prefix, "mz_profile.tsv")).expect("read mz profile");
    assert!(mz_profile.contains("mz_center"));

    let svg = fs::read_to_string(output_path(&out_prefix, "svg")).expect("read svg");
    assert!(svg.contains("Source:"));
    assert!(!svg.contains(&demo_mzml_path().display().to_string()));
    assert!(svg.contains("Target: m/z 500.0000 | all RT | mzML"));
    assert!(svg.contains("all RT | mzML"));
    assert!(!svg.contains("RT/m/z"));
    assert!(!svg.contains("Backend:"));
    assert!(!svg.contains("Capabilities:"));
    assert!(svg.contains("TIC "));
    assert!(svg.contains("RT profile"));
    assert!(svg.contains("m/z profile"));

    let _ = fs::remove_file(output_path(&out_prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&out_prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "mz_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "svg"));
    let _ = fs::remove_dir(&outdir);
}

#[test]
fn dia_slice_mzml_trace_peaks_writes_peak_and_trace_outputs() {
    let outdir = unique_prefix("trace");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--mz",
            "500",
            "--mz-da",
            "500",
            "--trace-peaks",
            "--emit-trace",
            "--outdir",
            outdir.to_str().expect("outdir"),
        ])
        .output()
        .expect("run dia-slice trace peaks on mzml");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let out_prefix = outdir.join("demo.dia_slice_mz_500.0000");
    let peaks = fs::read_to_string(output_path(&out_prefix, "peaks.tsv")).expect("read peaks");
    assert!(peaks.contains("# method\tselected-ion trace"));
    assert!(peaks.contains("signal_to_noise"));
    assert!(peaks.contains("peak_area_baseline_corrected"));

    let trace = fs::read_to_string(output_path(&out_prefix, "trace.tsv")).expect("read trace");
    assert!(trace.contains("# target\tm/z 500.0000"));
    assert!(trace.contains("target_intensity"));
    assert!(trace.contains("matched_events"));

    let _ = fs::remove_file(output_path(&out_prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&out_prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "mz_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "peaks.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "trace.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "svg"));
    let _ = fs::remove_dir(&outdir);
}

#[test]
fn dia_slice_mzml_honors_outdir_and_quiet() {
    let outdir = unique_prefix("outdir");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--mz",
            "500",
            "--mz-da",
            "500",
            "--outdir",
            outdir.to_str().expect("outdir"),
            "--quiet",
        ])
        .output()
        .expect("run dia-slice on mzml with outdir");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());

    let prefix = outdir.join("demo.dia_slice_mz_500.0000");
    let summary =
        fs::read_to_string(output_path(&prefix, "summary.txt")).expect("read summary output");
    assert!(summary.contains("backend\tmzML"));

    let _ = fs::remove_file(output_path(&prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&prefix, "mz_profile.tsv"));
    let _ = fs::remove_file(output_path(&prefix, "svg"));
    let _ = fs::remove_dir(&outdir);
}

#[test]
fn dia_slice_mzml_accepts_peptide_target_without_mz() {
    let outdir = unique_prefix("peptide");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--peptide",
            "GAIIGLMVGGVVIA",
            "--mz-da",
            "500",
            "--outdir",
            outdir.to_str().expect("outdir"),
        ])
        .output()
        .expect("run dia-slice on mzml with peptide");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let out_prefix = outdir.join("demo.dia_slice_GAIIGLMVGGVVIA_z2");
    let summary =
        fs::read_to_string(output_path(&out_prefix, "summary.txt")).expect("read summary output");
    assert!(summary.contains("peptide_sequence\tGAIIGLMVGGVVIA"));
    assert!(summary.contains("peptide_charge\t2"));
    assert!(summary.contains("peptide_precursor_mz\t"));

    let svg = fs::read_to_string(output_path(&out_prefix, "svg")).expect("read svg");
    assert!(svg.contains("Target: GAIIGLMVGGVVIA/2 precursor m/z"));

    let _ = fs::remove_file(output_path(&out_prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&out_prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "mz_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "svg"));
    let _ = fs::remove_dir(&outdir);
}

#[test]
fn dia_slice_mzml_accepts_peptide_fragment_target_without_mz() {
    let outdir = unique_prefix("fragment");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--peptide",
            "GAIIGLMVGGVVIA",
            "--fragment",
            "b8",
            "--mz-da",
            "500",
            "--outdir",
            outdir.to_str().expect("outdir"),
        ])
        .output()
        .expect("run dia-slice on mzml with peptide fragment");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let out_prefix = outdir.join("demo.dia_slice_GAIIGLMVGGVVIA_z2_b8");
    let summary =
        fs::read_to_string(output_path(&out_prefix, "summary.txt")).expect("read summary output");
    assert!(summary.contains("mz_center\t755.448408"));
    assert!(summary.contains("peptide_precursor_mz\t635.383591"));
    assert!(summary.contains("peptide_fragment_label\tb8"));
    assert!(summary.contains("peptide_fragment_mz\t755.448408"));

    let svg = fs::read_to_string(output_path(&out_prefix, "svg")).expect("read svg");
    assert!(svg.contains("fragment b8 m/z 755.4484"));
    assert!(svg.contains("precursor m/z 635.3836"));

    let _ = fs::remove_file(output_path(&out_prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&out_prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "mz_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "svg"));
    let _ = fs::remove_dir(&outdir);
}

#[test]
fn dia_slice_bruker_missing_runtime_fails_early() {
    let run_dir = unique_prefix("bruker-missing");
    fs::create_dir_all(&run_dir).expect("create fake bruker dir");
    let missing_so = unique_prefix("missing-runtime");

    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--bruker",
            run_dir.to_str().expect("run dir"),
            "--bruker-so",
            missing_so.to_str().expect("missing so path"),
            "--mz",
            "500",
        ])
        .output()
        .expect("run dia-slice with missing bruker runtime");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Bruker .d support requires timsdata.so"));
    assert!(stderr.contains("--bruker-so"));

    let _ = fs::remove_dir(&run_dir);
}

#[test]
fn dia_slice_bruker_slim_dia_fixture_when_available() {
    let bruker_path = local_slim_bruker_dia_path();
    if !bruker_path.exists() {
        return;
    }

    let outdir = unique_prefix("bruker-slim-dia");
    let output = Command::new(binary_path())
        .args([
            "dia-slice",
            "--bruker",
            bruker_path.to_str().expect("bruker fixture path"),
            "--bruker-backend",
            "native",
            "--peptide",
            "IILDLISESPIK/2",
            "--pseudo-ms2",
            "--neutral-losses",
            "--trace-peaks",
            "--emit-trace",
            "--mz-da",
            "0.025",
            "--rt-min",
            "21.60",
            "--rt-max",
            "21.78",
            "--outdir",
            outdir.to_str().expect("outdir"),
            "--quiet",
        ])
        .output()
        .expect("run dia-slice on slim Bruker DIA fixture");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());

    let out_prefix = outdir.join("hela_dia_iildlisespik_rt21_60_21_78.dia_slice_IILDLISESPIK_z2");
    let summary =
        fs::read_to_string(output_path(&out_prefix, "summary.txt")).expect("read summary output");
    assert!(summary.contains("backend\tBruker .d (timsrust native)"));
    assert!(summary.contains("acquisition_mode\tdiaPASEF"));
    assert!(summary.contains("capabilities\trt, mz, mobility, quad-window"));
    assert!(summary.contains("spectra_with_signal\t68"));

    let pseudo_ms2 =
        fs::read_to_string(output_path(&out_prefix, "pseudo_ms2.tsv")).expect("read pseudo-ms2");
    assert!(pseudo_ms2.contains("# peptide_input\tIILDLISESPIK/2"));
    assert!(pseudo_ms2.contains("frames_considered\t89"));
    let pseudo_ms2_svg = fs::read_to_string(output_path(&out_prefix, "pseudo_ms2.svg"))
        .expect("read pseudo-MS2 SVG");
    assert!(pseudo_ms2_svg.contains("N  sequence  C"));
    assert!(pseudo_ms2_svg.contains("sequence: IILDLISESPIK"));
    assert!(pseudo_ms2_svg.contains("phosphoric acid"));
    assert!(!pseudo_ms2_svg.contains("+w"));
    assert!(!pseudo_ms2_svg.contains("+n"));
    assert!(!pseudo_ms2_svg.contains("+p"));

    let trace = fs::read_to_string(output_path(&out_prefix, "trace.tsv")).expect("read trace");
    assert!(trace.contains("# target\tIILDLISESPIK/2 precursor m/z 670.9054"));

    let peaks = fs::read_to_string(output_path(&out_prefix, "peaks.tsv")).expect("read peaks");
    assert!(peaks.contains("peak_area_baseline_corrected"));

    let run_tic =
        fs::read_to_string(output_path(&out_prefix, "run_tic.json")).expect("read run TIC JSON");
    assert!(run_tic.contains("\"acquisition_mode\": \"diaPASEF\""));
    assert!(run_tic
        .contains("\"intensity_column\": \"analysis.tdf SummedIntensities / AccumulationTime\""));

    let _ = fs::remove_dir_all(&outdir);
}
