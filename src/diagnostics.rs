use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use mzdata::io::SpectrumSource;
use mzdata::spectrum::SignalContinuity;

use crate::annotate::MassTolerance;
use crate::mzml::{extract_scan_number, load_spectrum_by_index, open_reader, LoadedSpectrum};

const DEFAULT_TOLERANCE: MassTolerance = MassTolerance::Da(0.05);
const DEFAULT_TOP_N: usize = 40;
const DEFAULT_MIN_REL: f64 = 0.01;
const DEFAULT_MAX_DELTA_HITS: usize = 3;

#[derive(Clone, Debug)]
struct DiagnosticsOptions {
    mzml_path: Option<PathBuf>,
    tsv_path: Option<PathBuf>,
    panels: Vec<String>,
    target_specs: Vec<String>,
    delta_specs: Vec<String>,
    top_n: usize,
    min_rel: f64,
    max_delta_hits: usize,
    tolerance: MassTolerance,
}

impl Default for DiagnosticsOptions {
    fn default() -> Self {
        Self {
            mzml_path: None,
            tsv_path: None,
            panels: Vec::new(),
            target_specs: Vec::new(),
            delta_specs: Vec::new(),
            top_n: DEFAULT_TOP_N,
            min_rel: DEFAULT_MIN_REL,
            max_delta_hits: DEFAULT_MAX_DELTA_HITS,
            tolerance: DEFAULT_TOLERANCE,
        }
    }
}

#[derive(Clone, Debug)]
struct DiagnosticTarget {
    panel: String,
    label: String,
    mz: f64,
}

#[derive(Clone, Debug)]
struct DiagnosticDelta {
    panel: String,
    label: String,
    delta_da: f64,
}

#[derive(Clone, Debug)]
struct FilteredPeak {
    mz: f64,
    intensity: f64,
    rel: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum HitKind {
    Target,
    Delta,
}

impl HitKind {
    fn label(self) -> &'static str {
        match self {
            Self::Target => "target",
            Self::Delta => "delta",
        }
    }
}

#[derive(Clone, Debug)]
struct DiagnosticHit {
    kind: HitKind,
    panel: String,
    label: String,
    theoretical: f64,
    observed_a_mz: f64,
    observed_b_mz: Option<f64>,
    observed_value: f64,
    error_da: f64,
    error_ppm: f64,
    intensity_a: f64,
    rel_a: f64,
    intensity_b: Option<f64>,
    rel_b: Option<f64>,
}

pub fn run(args: Vec<String>) -> anyhow::Result<()> {
    if args
        .iter()
        .any(|arg| matches!(arg.as_str(), "-h" | "--help" | "help"))
    {
        print_help();
        return Ok(());
    }

    let options = parse_args(args)?;
    let mzml_path = options
        .mzml_path
        .as_ref()
        .expect("parse_args validates mzML path");

    let targets = build_targets(&options)?;
    let deltas = build_deltas(&options)?;
    if targets.is_empty() && deltas.is_empty() {
        anyhow::bail!("diagnostics requires at least one --panel, --target, or --delta");
    }

    let mut reader = open_reader(mzml_path)?;
    let total_spectra = reader.len();

    let tsv_path = options
        .tsv_path
        .clone()
        .unwrap_or_else(|| default_tsv_path(mzml_path));
    if let Some(parent) = tsv_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }

    let mut file = File::create(&tsv_path)
        .with_context(|| format!("failed to create {}", tsv_path.display()))?;
    writeln!(
        file,
        "file_stem\tscan_index\tscan_number\tscan_id\tms_level\trt_min\tprecursor_mz\tprecursor_charge\tkind\tpanel\tlabel\ttheoretical\tobserved_a_mz\tobserved_b_mz\tobserved_value\terror_da\terror_ppm\tintensity_a\trel_a\tintensity_b\trel_b"
    )?;

    let mut rows_written = 0usize;
    let mut hit_scans = BTreeSet::<u32>::new();
    let mut label_counts = BTreeMap::<(HitKind, String), usize>::new();
    let mut ms2_scanned = 0usize;
    let mut profile_skipped = 0usize;

    for idx in 0..total_spectra {
        let spectrum = load_spectrum_by_index(&mut reader, idx as u32)?;
        if spectrum.meta.ms_level != 2 {
            continue;
        }
        if matches!(spectrum.meta.continuity, SignalContinuity::Profile) {
            profile_skipped += 1;
            continue;
        }
        ms2_scanned += 1;

        let filtered_peaks = filter_peaks(&spectrum, options.top_n, options.min_rel);
        if filtered_peaks.is_empty() {
            continue;
        }

        let mut scan_hits = Vec::<DiagnosticHit>::new();
        for target in &targets {
            if let Some(hit) = find_target_ion(&filtered_peaks, target, options.tolerance) {
                scan_hits.push(hit);
            }
        }
        for delta in &deltas {
            scan_hits.extend(find_peak_pairs_by_delta(
                &filtered_peaks,
                delta,
                options.tolerance,
                options.max_delta_hits,
            ));
        }

        if scan_hits.is_empty() {
            continue;
        }

        hit_scans.insert(spectrum.meta.idx);
        let file_stem = mzml_path
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("mzml");
        let scan_number = extract_scan_number(&spectrum.meta.scan_id);
        for hit in scan_hits {
            *label_counts
                .entry((hit.kind, hit.label.clone()))
                .or_insert(0) += 1;
            write_hit_row(&mut file, file_stem, &spectrum, scan_number, &hit)?;
            rows_written += 1;
        }
    }

    println!("Wrote TSV: {}", tsv_path.display());
    println!(
        "Scanned {} spectra | centroid/unknown MS2 {} | profile MS2 skipped {} | hits {} across {} scans",
        total_spectra, ms2_scanned, profile_skipped, rows_written, hit_scans.len()
    );
    if !label_counts.is_empty() {
        println!("Top labels:");
        let mut top_labels = label_counts.into_iter().collect::<Vec<_>>();
        top_labels.sort_by(|left, right| {
            right
                .1
                .cmp(&left.1)
                .then_with(|| left.0 .0.cmp(&right.0 .0))
                .then_with(|| left.0 .1.cmp(&right.0 .1))
        });
        for ((kind, label), count) in top_labels.into_iter().take(10) {
            println!("  {} {}: {}", kind.label(), label, count);
        }
    }

    Ok(())
}

fn parse_args(args: Vec<String>) -> anyhow::Result<DiagnosticsOptions> {
    let mut options = DiagnosticsOptions::default();
    let mut tolerance_override = None::<MassTolerance>;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                options.mzml_path =
                    Some(PathBuf::from(iter.next().context("--mzml expects a path")?));
            }
            "--tsv" | "--output" => {
                options.tsv_path = Some(PathBuf::from(
                    iter.next().context("--tsv expects a file path")?,
                ));
            }
            "--panel" => {
                options
                    .panels
                    .push(iter.next().context("--panel expects a panel name")?);
            }
            "--target" => {
                options.target_specs.push(
                    iter.next()
                        .context("--target expects <label>:<mz> or <mz>")?,
                );
            }
            "--delta" => {
                options.delta_specs.push(
                    iter.next()
                        .context("--delta expects <label>:<da> or <da>")?,
                );
            }
            "--top-n" => {
                let raw = iter.next().context("--top-n expects an integer")?;
                let value = raw.parse::<usize>().context("invalid --top-n")?;
                if value == 0 {
                    anyhow::bail!("--top-n must be at least 1");
                }
                options.top_n = value;
            }
            "--min-rel" => {
                let raw = iter
                    .next()
                    .context("--min-rel expects a float between 0 and 1")?;
                let value = raw.parse::<f64>().context("invalid --min-rel")?;
                if !(0.0..=1.0).contains(&value) || !value.is_finite() {
                    anyhow::bail!("--min-rel must be between 0 and 1");
                }
                options.min_rel = value;
            }
            "--max-delta-hits" => {
                let raw = iter.next().context("--max-delta-hits expects an integer")?;
                let value = raw.parse::<usize>().context("invalid --max-delta-hits")?;
                if value == 0 {
                    anyhow::bail!("--max-delta-hits must be at least 1");
                }
                options.max_delta_hits = value;
            }
            "--tol-ppm" => {
                let raw = iter.next().context("--tol-ppm expects a float")?;
                let ppm = raw.parse::<f64>().context("invalid --tol-ppm")?;
                if ppm <= 0.0 || !ppm.is_finite() {
                    anyhow::bail!("--tol-ppm must be a positive finite number");
                }
                set_tolerance(
                    &mut tolerance_override,
                    MassTolerance::Ppm(ppm),
                    "--tol-ppm",
                )?;
            }
            "--tol-da" => {
                let raw = iter.next().context("--tol-da expects a float")?;
                let da = raw.parse::<f64>().context("invalid --tol-da")?;
                if da <= 0.0 || !da.is_finite() {
                    anyhow::bail!("--tol-da must be a positive finite number");
                }
                set_tolerance(&mut tolerance_override, MassTolerance::Da(da), "--tol-da")?;
            }
            other => anyhow::bail!("unknown diagnostics option `{other}`"),
        }
    }

    options.tolerance = tolerance_override.unwrap_or(DEFAULT_TOLERANCE);
    if options.mzml_path.is_none() {
        anyhow::bail!("diagnostics requires --mzml <path>");
    }

    Ok(options)
}

fn set_tolerance(
    slot: &mut Option<MassTolerance>,
    tolerance: MassTolerance,
    flag: &str,
) -> anyhow::Result<()> {
    if slot.is_some() {
        anyhow::bail!("specify only one of --tol-ppm or --tol-da (conflict at {flag})");
    }
    *slot = Some(tolerance);
    Ok(())
}

fn print_help() {
    let program = crate::program_name();
    println!("{program} diagnostics");
    println!();
    println!("USAGE:");
    println!(
        "  {program} diagnostics --mzml <file> [--panel <name>] [--target <spec>] [--delta <spec>] [options]"
    );
    println!();
    println!("OPTIONS:");
    println!("  --mzml <file>               Input mzML file");
    println!(
        "  --tsv <file>                Output TSV path (default: exports/...diagnostics...tsv)"
    );
    println!("  --panel <name>              Built-in panel: oxonium, glyco, phospho (repeatable)");
    println!(
        "  --target <spec>             Custom target ion as <label>:<mz> or bare <mz> (repeatable)"
    );
    println!(
        "  --delta <spec>              Custom delta as <label>:<da> or bare <da> (repeatable)"
    );
    println!("  --top-n <n>                 Keep at most this many peaks per scan after filtering (default: 40)");
    println!("  --min-rel <f>               Minimum relative intensity before top-N filtering (default: 0.01)");
    println!("  --max-delta-hits <n>        Report up to this many pairwise matches per delta per scan (default: 3)");
    println!("  --tol-da <da>               Matching tolerance in Daltons (default: 0.05)");
    println!("  --tol-ppm <ppm>             Matching tolerance in ppm");
    println!("  --help                      Show this help");
    println!();
    println!("EXAMPLES:");
    println!("  {program} diagnostics --mzml sample.mzML --panel oxonium");
    println!("  {program} diagnostics --mzml sample.mzML --panel phospho --delta water:18.010565");
    println!(
        "  {program} diagnostics --mzml sample.mzML --target HexNAc:204.08665 --delta phospho_nl:97.976896 --top-n 60 --min-rel 0.02"
    );
}

fn build_targets(options: &DiagnosticsOptions) -> anyhow::Result<Vec<DiagnosticTarget>> {
    let mut out = Vec::<DiagnosticTarget>::new();
    for panel in &options.panels {
        out.extend(panel_targets(panel)?);
    }
    for raw in &options.target_specs {
        out.push(parse_target_spec(raw)?);
    }
    dedup_targets(&mut out);
    Ok(out)
}

fn build_deltas(options: &DiagnosticsOptions) -> anyhow::Result<Vec<DiagnosticDelta>> {
    let mut out = Vec::<DiagnosticDelta>::new();
    for panel in &options.panels {
        out.extend(panel_deltas(panel)?);
    }
    for raw in &options.delta_specs {
        out.push(parse_delta_spec(raw)?);
    }
    dedup_deltas(&mut out);
    Ok(out)
}

fn panel_targets(panel: &str) -> anyhow::Result<Vec<DiagnosticTarget>> {
    match panel {
        "oxonium" | "glyco" => Ok(vec![
            target("oxonium", "oxonium_126", 126.0550),
            target("oxonium", "oxonium_138", 138.0550),
            target("oxonium", "oxonium_144", 144.0655),
            target("oxonium", "oxonium_168", 168.0655),
            target("oxonium", "oxonium_186", 186.0761),
            target("oxonium", "hexnac_204", 204.0867),
            target("oxonium", "neuac_274", 274.0921),
            target("oxonium", "neuac_292", 292.1027),
            target("oxonium", "hexhexnac_366", 366.1395),
        ]),
        "phospho" => Ok(Vec::new()),
        other => anyhow::bail!("unknown --panel `{other}` (expected oxonium, glyco, or phospho)"),
    }
}

fn panel_deltas(panel: &str) -> anyhow::Result<Vec<DiagnosticDelta>> {
    match panel {
        "oxonium" | "glyco" => Ok(Vec::new()),
        "phospho" => Ok(vec![delta("phospho", "phospho_nl", 97.976_896)]),
        other => anyhow::bail!("unknown --panel `{other}` (expected oxonium, glyco, or phospho)"),
    }
}

fn target(panel: &str, label: &str, mz: f64) -> DiagnosticTarget {
    DiagnosticTarget {
        panel: panel.to_string(),
        label: label.to_string(),
        mz,
    }
}

fn delta(panel: &str, label: &str, delta_da: f64) -> DiagnosticDelta {
    DiagnosticDelta {
        panel: panel.to_string(),
        label: label.to_string(),
        delta_da,
    }
}

fn parse_target_spec(raw: &str) -> anyhow::Result<DiagnosticTarget> {
    if let Some((label, mz_raw)) = raw.rsplit_once(':') {
        let mz = mz_raw
            .parse::<f64>()
            .with_context(|| format!("invalid target mz in `{raw}`"))?;
        return Ok(target("custom", label, mz));
    }
    let mz = raw
        .parse::<f64>()
        .with_context(|| format!("invalid target mz `{raw}`"))?;
    Ok(target(
        "custom",
        &format!("mz{mz:.4}").replace('.', "_"),
        mz,
    ))
}

fn parse_delta_spec(raw: &str) -> anyhow::Result<DiagnosticDelta> {
    if let Some((label, delta_raw)) = raw.rsplit_once(':') {
        let delta_da = delta_raw
            .parse::<f64>()
            .with_context(|| format!("invalid delta mass in `{raw}`"))?;
        return Ok(delta("custom", label, delta_da));
    }
    let delta_da = raw
        .parse::<f64>()
        .with_context(|| format!("invalid delta mass `{raw}`"))?;
    Ok(delta(
        "custom",
        &format!("delta{delta_da:.4}").replace('.', "_"),
        delta_da,
    ))
}

fn dedup_targets(targets: &mut Vec<DiagnosticTarget>) {
    targets.sort_by(|left, right| {
        left.label
            .cmp(&right.label)
            .then_with(|| left.panel.cmp(&right.panel))
            .then_with(|| left.mz.partial_cmp(&right.mz).unwrap_or(Ordering::Equal))
    });
    targets.dedup_by(|left, right| {
        left.label == right.label && left.panel == right.panel && (left.mz - right.mz).abs() <= 1e-9
    });
}

fn dedup_deltas(deltas: &mut Vec<DiagnosticDelta>) {
    deltas.sort_by(|left, right| {
        left.label
            .cmp(&right.label)
            .then_with(|| left.panel.cmp(&right.panel))
            .then_with(|| {
                left.delta_da
                    .partial_cmp(&right.delta_da)
                    .unwrap_or(Ordering::Equal)
            })
    });
    deltas.dedup_by(|left, right| {
        left.label == right.label
            && left.panel == right.panel
            && (left.delta_da - right.delta_da).abs() <= 1e-9
    });
}

fn filter_peaks(spectrum: &LoadedSpectrum, top_n: usize, min_rel: f64) -> Vec<FilteredPeak> {
    let base_peak = spectrum.stats.base_peak_intensity.max(1.0e-6) as f64;
    let mut peaks = spectrum
        .mz
        .iter()
        .copied()
        .zip(spectrum.intensity.iter().copied())
        .filter_map(|(mz, intensity)| {
            let intensity = intensity as f64;
            if !mz.is_finite() || !intensity.is_finite() || intensity <= 0.0 {
                return None;
            }
            let rel = intensity / base_peak;
            if rel < min_rel {
                return None;
            }
            Some(FilteredPeak { mz, intensity, rel })
        })
        .collect::<Vec<_>>();

    peaks.sort_by(|left, right| {
        right
            .intensity
            .partial_cmp(&left.intensity)
            .unwrap_or(Ordering::Equal)
    });
    peaks.truncate(top_n.max(1));
    peaks.sort_by(|left, right| left.mz.partial_cmp(&right.mz).unwrap_or(Ordering::Equal));
    peaks
}

fn find_target_ion(
    peaks: &[FilteredPeak],
    target: &DiagnosticTarget,
    tolerance: MassTolerance,
) -> Option<DiagnosticHit> {
    let mut best: Option<&FilteredPeak> = None;
    for peak in peaks {
        if !tolerance.contains(target.mz, peak.mz) {
            continue;
        }
        match best {
            None => best = Some(peak),
            Some(current) => {
                let peak_error = (peak.mz - target.mz).abs();
                let current_error = (current.mz - target.mz).abs();
                if peak_error < current_error
                    || ((peak_error - current_error).abs() <= 1e-12
                        && peak.intensity > current.intensity)
                {
                    best = Some(peak);
                }
            }
        }
    }

    best.map(|peak| DiagnosticHit {
        kind: HitKind::Target,
        panel: target.panel.clone(),
        label: target.label.clone(),
        theoretical: target.mz,
        observed_a_mz: peak.mz,
        observed_b_mz: None,
        observed_value: peak.mz,
        error_da: peak.mz - target.mz,
        error_ppm: tolerance.error_ppm(target.mz, peak.mz),
        intensity_a: peak.intensity,
        rel_a: peak.rel,
        intensity_b: None,
        rel_b: None,
    })
}

fn find_peak_pairs_by_delta(
    peaks: &[FilteredPeak],
    delta: &DiagnosticDelta,
    tolerance: MassTolerance,
    max_hits: usize,
) -> Vec<DiagnosticHit> {
    let mut hits = Vec::<DiagnosticHit>::new();
    for left_index in 0..peaks.len() {
        for right_index in (left_index + 1)..peaks.len() {
            let left = &peaks[left_index];
            let right = &peaks[right_index];
            let observed_delta = right.mz - left.mz;
            if !tolerance.contains(delta.delta_da, observed_delta) {
                continue;
            }
            hits.push(DiagnosticHit {
                kind: HitKind::Delta,
                panel: delta.panel.clone(),
                label: delta.label.clone(),
                theoretical: delta.delta_da,
                observed_a_mz: left.mz,
                observed_b_mz: Some(right.mz),
                observed_value: observed_delta,
                error_da: observed_delta - delta.delta_da,
                error_ppm: tolerance.error_ppm(delta.delta_da, observed_delta),
                intensity_a: left.intensity,
                rel_a: left.rel,
                intensity_b: Some(right.intensity),
                rel_b: Some(right.rel),
            });
        }
    }

    hits.sort_by(|left, right| {
        let left_score = left.rel_a + left.rel_b.unwrap_or(0.0);
        let right_score = right.rel_a + right.rel_b.unwrap_or(0.0);
        right_score
            .partial_cmp(&left_score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                left.error_da
                    .abs()
                    .partial_cmp(&right.error_da.abs())
                    .unwrap_or(Ordering::Equal)
            })
    });
    hits.truncate(max_hits.max(1));
    hits
}

fn write_hit_row(
    file: &mut File,
    file_stem: &str,
    spectrum: &LoadedSpectrum,
    scan_number: Option<u64>,
    hit: &DiagnosticHit,
) -> anyhow::Result<()> {
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{:.6}\t{:+.6}\t{:+.2}\t{:.3}\t{:.4}\t{}\t{}",
        sanitize_tsv(file_stem),
        spectrum.meta.idx,
        scan_number
            .map(|value| value.to_string())
            .unwrap_or_default(),
        sanitize_tsv(&spectrum.meta.scan_id),
        spectrum.meta.ms_level,
        spectrum
            .meta
            .rt_minutes
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default(),
        spectrum
            .meta
            .precursor_mz
            .map(|value| format!("{value:.6}"))
            .unwrap_or_default(),
        spectrum
            .meta
            .precursor_charge
            .map(|value| value.to_string())
            .unwrap_or_default(),
        hit.kind.label(),
        sanitize_tsv(&hit.panel),
        sanitize_tsv(&hit.label),
        hit.theoretical,
        hit.observed_a_mz,
        hit.observed_b_mz
            .map(|value| format!("{value:.6}"))
            .unwrap_or_default(),
        hit.observed_value,
        hit.error_da,
        hit.error_ppm,
        hit.intensity_a,
        hit.rel_a,
        hit.intensity_b
            .map(|value| format!("{value:.3}"))
            .unwrap_or_default(),
        hit.rel_b
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default(),
    )?;
    Ok(())
}

fn sanitize_tsv(input: &str) -> String {
    input
        .replace('\t', " ")
        .replace('\n', " ")
        .replace('\r', " ")
}

fn default_tsv_path(mzml_path: &Path) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let stem = mzml_path
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("mzml");
    PathBuf::from("exports").join(format!("{stem}__diagnostics__{ts}.tsv"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_target_spec_with_implicit_label() {
        let target = parse_target_spec("204.08665").expect("target spec");
        assert_eq!(target.panel, "custom");
        assert!((target.mz - 204.08665).abs() < 1e-9);
        assert!(target.label.contains("204_0866"));
    }

    #[test]
    fn parses_delta_spec_with_explicit_label() {
        let delta = parse_delta_spec("phospho_nl:97.976896").expect("delta spec");
        assert_eq!(delta.label, "phospho_nl");
        assert!((delta.delta_da - 97.976896).abs() < 1e-9);
    }

    #[test]
    fn panel_builders_include_expected_entries() {
        let targets = panel_targets("oxonium").expect("oxonium panel");
        let deltas = panel_deltas("phospho").expect("phospho panel");
        assert!(targets.iter().any(|target| target.label == "hexnac_204"));
        assert!(deltas.iter().any(|delta| delta.label == "phospho_nl"));
    }

    #[test]
    fn find_target_ion_picks_best_match() {
        let peaks = vec![
            FilteredPeak {
                mz: 204.0400,
                intensity: 20.0,
                rel: 0.2,
            },
            FilteredPeak {
                mz: 204.0867,
                intensity: 10.0,
                rel: 0.1,
            },
        ];
        let hit = find_target_ion(
            &peaks,
            &target("oxonium", "hexnac_204", 204.08665),
            MassTolerance::Da(0.05),
        )
        .expect("target hit");
        assert!((hit.observed_a_mz - 204.0867).abs() < 1e-6);
    }

    #[test]
    fn find_peak_pairs_by_delta_reports_combined_intensity_sorted_hits() {
        let peaks = vec![
            FilteredPeak {
                mz: 500.0,
                intensity: 200.0,
                rel: 0.4,
            },
            FilteredPeak {
                mz: 597.9769,
                intensity: 180.0,
                rel: 0.36,
            },
            FilteredPeak {
                mz: 700.0,
                intensity: 50.0,
                rel: 0.1,
            },
            FilteredPeak {
                mz: 797.9770,
                intensity: 40.0,
                rel: 0.08,
            },
        ];
        let hits = find_peak_pairs_by_delta(
            &peaks,
            &delta("phospho", "phospho_nl", 97.976896),
            MassTolerance::Da(0.02),
            2,
        );
        assert_eq!(hits.len(), 2);
        assert!((hits[0].observed_a_mz - 500.0).abs() < 1e-6);
        assert!((hits[0].observed_b_mz.unwrap_or_default() - 597.9769).abs() < 1e-6);
    }

    #[test]
    fn filter_peaks_applies_relative_cutoff_before_top_n() {
        let spectrum = LoadedSpectrum {
            meta: crate::mzml::SpectrumMeta {
                idx: 1,
                scan_id: "scan=1".to_string(),
                ms_level: 2,
                rt_minutes: None,
                precursor_mz: None,
                precursor_charge: None,
                continuity: SignalContinuity::Centroid,
            },
            mz: vec![100.0, 200.0, 300.0],
            intensity: vec![1000.0, 40.0, 5.0],
            stats: crate::mzml::SpectrumStats {
                points: 3,
                mz_min: 100.0,
                mz_max: 300.0,
                base_peak_mz: 100.0,
                base_peak_intensity: 1000.0,
            },
        };
        let peaks = filter_peaks(&spectrum, 5, 0.03);
        assert_eq!(peaks.len(), 2);
        assert!(peaks.iter().all(|peak| peak.rel >= 0.03));
    }
}
