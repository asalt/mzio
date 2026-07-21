use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn binary_path() -> &'static str {
    env!("CARGO_BIN_EXE_mzio")
}

fn demo_mzml_path() -> PathBuf {
    repo_path("assets/demo.mzML")
}

fn demo_ms2_path() -> PathBuf {
    repo_path("assets/demo.ms2")
}

fn repo_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(relative)
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

fn unique_temp_dir(stem: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "mzml-spinoff-{stem}-{}-{timestamp}",
        std::process::id()
    ))
}

fn output_paths_with_extension(dir: &std::path::Path, extension: &str) -> Vec<PathBuf> {
    let mut paths = fs::read_dir(dir)
        .unwrap_or_else(|err| panic!("read {}: {err}", dir.display()))
        .map(|entry| entry.expect("read export entry").path())
        .filter(|path| {
            path.extension()
                .and_then(|value| value.to_str())
                .is_some_and(|value| value == extension)
        })
        .collect::<Vec<_>>();
    paths.sort();
    paths
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
    assert!(stdout.contains("--pepxml <file>"));
    assert!(stdout.contains("--top-n <n>"));
    assert!(stdout.contains("--mod <position>:<delta>"));
    assert!(stdout.contains("--neutral-losses"));
    assert!(stdout.contains("--neutral-loss-min-frac <f>"));
    assert!(stdout.contains("--isotope-errors <list>"));
    assert!(stdout.contains("--remove-precursor"));
    assert!(stdout.contains("--tol-da <da>"));
}

#[test]
fn plot_pepxml_writes_svg_and_json_sidecar() {
    let temp_dir = unique_temp_dir("pepxml-plot");
    fs::create_dir_all(&temp_dir).expect("create temp dir");
    let pepxml_path = temp_dir.join("search.pep.xml");
    let svg_path = temp_dir.join("scan2.svg");
    fs::write(
        &pepxml_path,
        r#"
<msms_pipeline_analysis>
  <msms_run_summary>
    <spectrum_query spectrum="demo.2.2.2" start_scan="2" end_scan="2" assumed_charge="2">
      <search_result>
        <search_hit hit_rank="1" peptide="PEPMIDEK" protein="sp|P1" calc_neutral_pep_mass="900.1" massdiff="0.01">
          <modification_info mod_nterm_mass="42.0106">
            <mod_aminoacid_mass position="4" mass="147.035384645"/>
          </modification_info>
          <search_score name="hyperscore" value="51.2"/>
        </search_hit>
      </search_result>
    </spectrum_query>
  </msms_run_summary>
</msms_pipeline_analysis>
"#,
    )
    .expect("write pepXML");

    let output = Command::new(binary_path())
        .args([
            "plot",
            "--mzml",
            demo_mzml_path().to_str().expect("demo mzml path"),
            "--id",
            "scan=2",
            "--pepxml",
            pepxml_path.to_str().expect("pepxml path"),
            "--top-n",
            "3",
            "--svg",
            svg_path.to_str().expect("svg path"),
        ])
        .output()
        .expect("run pepXML plot");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Wrote SVG:"));
    assert!(stdout.contains("Wrote JSON:"));
    assert!(stdout.contains("pepXML hit rank 1"));

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("requested --top-n 3 but pepXML has 1 hit"));

    let svg = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg.contains("Peptide: [+42.010600]PEPM[+15.994900]IDEK"));
    assert!(svg.contains("Applied mods: n-term:+42.010600, 4:+15.994900"));

    let json_path = svg_path.with_extension("json");
    let json = fs::read_to_string(&json_path).expect("read json");
    assert!(json.contains("\"rank\": 1"));
    assert!(json.contains("\"peptide\": \"PEPMIDEK\""));
    assert!(json.contains("\"modified_sequence\": \"[+42.010600]PEPM[+15.994900]IDEK\""));
    assert!(json.contains("\"name\": \"hyperscore\""));
    assert!(json.contains("\"delta\": 15.9949"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn plot_committed_hela_fixture_writes_ranked_outputs() {
    let temp_dir = unique_temp_dir("hela-fixture-plot");
    fs::create_dir_all(&temp_dir).expect("create temp dir");
    let mzml_path = repo_path(
        "data/exploris480_hela_digest/mzml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065.mzML",
    );
    let pepxml_path = repo_path(
        "data/exploris480_hela_digest/pepxml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065_top_hits.pep.xml",
    );

    assert!(mzml_path.is_file(), "missing {}", mzml_path.display());
    assert!(pepxml_path.is_file(), "missing {}", pepxml_path.display());

    let output = Command::new(binary_path())
        .current_dir(&temp_dir)
        .args([
            "plot",
            "--mzml",
            mzml_path.to_str().expect("mzML fixture path"),
            "--scan",
            "9065",
            "--pepxml",
            pepxml_path.to_str().expect("pepXML fixture path"),
            "--top-n",
            "3",
            "--svg-prefix",
            "exploris480_hela_digest",
        ])
        .output()
        .expect("run HeLa fixture plot");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("pepXML hit rank 1"));
    assert!(stdout.contains("pepXML hit rank 2"));
    assert!(stdout.contains("pepXML hit rank 3"));

    let exports_dir = temp_dir.join("exports");
    let svg_paths = output_paths_with_extension(&exports_dir, "svg");
    let json_paths = output_paths_with_extension(&exports_dir, "json");
    assert_eq!(svg_paths.len(), 3);
    assert_eq!(json_paths.len(), 3);
    assert!(svg_paths.iter().all(|path| path
        .file_name()
        .unwrap()
        .to_string_lossy()
        .contains("rank")));

    let json = json_paths
        .iter()
        .map(|path| fs::read_to_string(path).expect("read json sidecar"))
        .collect::<Vec<_>>()
        .join("\n");
    assert!(json.contains("\"rank\": 1"));
    assert!(json.contains("\"rank\": 3"));
    assert!(json.contains("\"peptide\": \"MNTNPSR\""));
    assert!(json.contains("\"name\": \"hyperscore\""));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn plot_committed_phospho_fixture_writes_localization_outputs() {
    let temp_dir = unique_temp_dir("phospho-fixture-plot");
    fs::create_dir_all(&temp_dir).expect("create temp dir");
    let mzml_path = repo_path(
        "data/exploris480_phospho_localization/mzml/49108_1_EXP_802_PDX_1mg_phos_oneforth_F13_scan4877.mzML",
    );
    let pepxml_path = repo_path(
        "data/exploris480_phospho_localization/pepxml/49108_F13_scan4877_phospho_top_hits.pep.xml",
    );

    assert!(mzml_path.is_file(), "missing {}", mzml_path.display());
    assert!(pepxml_path.is_file(), "missing {}", pepxml_path.display());

    let output = Command::new(binary_path())
        .current_dir(&temp_dir)
        .args([
            "plot",
            "--mzml",
            mzml_path.to_str().expect("mzML fixture path"),
            "--scan",
            "4877",
            "--pepxml",
            pepxml_path.to_str().expect("pepXML fixture path"),
            "--top-n",
            "5",
            "--neutral-losses",
            "--svg-prefix",
            "exploris480_phospho_localization",
        ])
        .output()
        .expect("run phospho fixture plot");

    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("pepXML hit rank 1"));
    assert!(stdout.contains("pepXML hit rank 5"));
    assert!(stdout.contains("TGS[+79.966324]ESSQTGTSTTSSR"));

    let exports_dir = temp_dir.join("exports");
    let svg_paths = output_paths_with_extension(&exports_dir, "svg");
    let json_paths = output_paths_with_extension(&exports_dir, "json");
    assert_eq!(svg_paths.len(), 5);
    assert_eq!(json_paths.len(), 5);
    assert!(svg_paths.iter().all(|path| path
        .file_name()
        .unwrap()
        .to_string_lossy()
        .contains("nl-on")));
    let rank1_svg_path = svg_paths
        .iter()
        .find(|path| {
            path.file_name()
                .unwrap()
                .to_string_lossy()
                .contains("rank1")
        })
        .expect("rank 1 SVG");
    let rank1_svg = fs::read_to_string(rank1_svg_path).expect("read rank 1 SVG");
    assert!(rank1_svg.contains("N  sequence  C"));
    assert!(rank1_svg.contains("sequence: TGS[+79.966324]ESSQTGTSTTSSR"));
    assert!(rank1_svg.contains("phosphoric acid"));
    assert!(rank1_svg.contains("spectrum-peak-matched-b"));
    assert!(rank1_svg.contains("spectrum-peak-matched-y"));
    assert!(!rank1_svg.contains("+p"));

    let json = json_paths
        .iter()
        .map(|path| fs::read_to_string(path).expect("read json sidecar"))
        .collect::<Vec<_>>()
        .join("\n");
    assert!(json.contains("\"rank\": 5"));
    assert!(json.contains("\"peptide\": \"TGSESSQTGTSTTSSR\""));
    assert!(json.contains("\"modified_sequence\": \"TGS[+79.966324]ESSQTGTSTTSSR\""));
    assert!(json.contains("\"frag_error_mae_da\""));

    let _ = fs::remove_dir_all(temp_dir);
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
    assert!(svg.contains("colored = matched evidence; grey = unmatched theoretical base ion"));
    assert!(svg.contains("loss lines = matched evidence:"));
    assert!(svg.contains("N  sequence  C"));
    assert!(!svg.contains("+w"));
    assert!(!svg.contains("+n"));
    assert!(!svg.contains("+p"));

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
