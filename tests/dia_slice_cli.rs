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
    assert!(stdout.contains("--mz-da <da>"));
    assert!(stdout.contains("all RT"));
    assert!(stdout.contains("all mobility"));
    assert!(stdout.contains("all windows"));
    assert!(stdout.contains("--quad-min <mz>"));
    assert!(stdout.contains("--outdir <dir>"));
    assert!(stdout.contains("--rt-smooth"));
    assert!(stdout.contains("--rt-smooth-window <n>"));
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
    assert!(svg.contains("Target: m/z 500.0000 | window 0.0000-1000.0000"));
    assert!(svg.contains("mzML | RT/m/z"));
    assert!(!svg.contains("Backend:"));
    assert!(!svg.contains("Capabilities:"));
    assert!(svg.contains("RT profile"));
    assert!(svg.contains("m/z profile"));

    let _ = fs::remove_file(output_path(&out_prefix, "summary.txt"));
    let _ = fs::remove_file(output_path(&out_prefix, "rt_profile.tsv"));
    let _ = fs::remove_file(output_path(&out_prefix, "mz_profile.tsv"));
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
