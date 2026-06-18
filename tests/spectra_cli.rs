use std::process::Command;

fn binary_path() -> &'static str {
    env!("CARGO_BIN_EXE_mzio")
}

#[test]
fn root_help_mentions_spectra() {
    let output = Command::new(binary_path())
        .arg("--help")
        .output()
        .expect("run root --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("spectra"));
    assert!(stdout.contains("Browse mzML spectra"));
    assert!(stdout.contains("plot-survey"));
}

#[test]
fn spectra_help_mentions_cache_flags() {
    let output = Command::new(binary_path())
        .args(["spectra", "--help"])
        .output()
        .expect("run spectra --help");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--mzml <file>"));
    assert!(stdout.contains("--no-mzml-cache"));
    assert!(stdout.contains("--reindex"));
    assert!(stdout.contains("--mzml-cache-path <p>"));
}
