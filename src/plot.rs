use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use mzdata::spectrum::SignalContinuity;
use serde::Serialize;

use crate::annotate::{
    annotate_peaks, prepare_annotation, prepare_annotation_with_modifications, AnnotationContext,
    AnnotationQualityMetrics, AnnotationReport, FragmentIon, FragmentMatch, FragmentSeries,
    MassTolerance, NeutralLossKind, PrecursorCheck, DEFAULT_PRECURSOR_ISOTOPE_ERRORS,
};
use crate::ms2::load_selected_spectrum as load_selected_ms2_spectrum;
use crate::mzml::{
    extract_scan_number, load_selected_spectrum as load_selected_mzml_spectrum, open_reader,
    LoadedSpectrum, SpectrumSelector,
};
use crate::pepxml::{load_hits_for_scan, PepXmlHit, PepXmlModification, PepXmlScore};
use crate::scale::CoordinateRange;
use crate::svg_canvas::{AxisOrientation, AxisProps, AxisTickLabelStyle, SvgCanvas};

const SVG_WIDTH: u32 = 1480;
const SVG_HEIGHT: u32 = 940;
const SVG_BINS: usize = 4000;
const ION_TABLE_MIN_WIDTH: f64 = 620.0;
const ION_TABLE_MAX_EXTRA_WIDTH: f64 = 520.0;
const ION_TABLE_GAP: f64 = 18.0;
const PLOT_HEADER_TITLE_FONT: u32 = 20;
const PLOT_HEADER_META_FONT: u32 = 14;
const PLOT_HEADER_DETAIL_FONT: u32 = 15;
const PLOT_TICK_FONT: u32 = 18;
const PLOT_AXIS_LABEL_FONT: u32 = 20;
const PLOT_PEAK_LABEL_FONT: u32 = 13;
const PLOT_PRECURSOR_LABEL_FONT: u32 = 13;
const LADDER_TITLE_FONT: u32 = 14;
const LADDER_RESIDUE_FONT: u32 = 22;
const LADDER_INDEX_FONT: u32 = 12;
const LADDER_ION_FONT: u32 = 13;
const ION_TABLE_TITLE_FONT: u32 = 18;
const ION_TABLE_META_FONT: u32 = 13;
const ION_TABLE_HEADER_FONT: f64 = 13.5;
const COLOR_TEXT: &str = "#122033";
const COLOR_SUBTLE: &str = "#5b6775";
const COLOR_WARNING: &str = "#b45309";
const COLOR_SERIES_B: &str = "#2563eb";
const COLOR_SERIES_Y: &str = "#c2410c";
const COLOR_SERIES_NEUTRAL: &str = "#94a3b8";
const COLOR_PLOT: &str = "#0f766e";
const COLOR_CARD_BORDER: &str = "#d8e0ea";
const COLOR_AXIS: &str = "#334155";
const FONT_TABLE: &str = "Menlo, Consolas, Liberation Mono, monospace";
const DEFAULT_TOLERANCE: MassTolerance = MassTolerance::Ppm(20.0);
const COMMON_NEUTRAL_LOSSES: [NeutralLossKind; 3] = [
    NeutralLossKind::Water,
    NeutralLossKind::Ammonia,
    NeutralLossKind::PhosphoricAcid,
];
const DEFAULT_NEUTRAL_LOSS_LABEL_MIN_FRAC: f64 = 0.03;
const PRECURSOR_REMOVAL_HALF_WIDTH_DA: f64 = 1.5;

#[derive(Clone, Copy, Debug)]
enum PlotMode {
    Auto,
    Sticks,
    Line,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlotInputFormat {
    Mzml,
    Ms2,
}

#[derive(Clone, Copy, Debug)]
enum PlotRenderMode {
    Sticks,
    Line,
}

#[derive(Clone, Debug)]
struct PlotOptions {
    input_path: Option<PathBuf>,
    input_format: Option<PlotInputFormat>,
    selector: Option<SpectrumSelector>,
    svg_path: Option<PathBuf>,
    svg_prefix: Option<String>,
    peptide_input: Option<String>,
    pepxml_path: Option<PathBuf>,
    top_n: usize,
    top_n_explicit: bool,
    mod_inputs: Vec<String>,
    neutral_losses_enabled: bool,
    neutral_loss_label_min_frac: f64,
    isotope_errors: Vec<u8>,
    isotope_errors_explicit: bool,
    charge_override: Option<i32>,
    remove_precursor: bool,
    tolerance: MassTolerance,
    tolerance_explicit: bool,
    normalize: bool,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    mode: PlotMode,
}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            input_path: None,
            input_format: None,
            selector: None,
            svg_path: None,
            svg_prefix: None,
            peptide_input: None,
            pepxml_path: None,
            top_n: 1,
            top_n_explicit: false,
            mod_inputs: Vec::new(),
            neutral_losses_enabled: false,
            neutral_loss_label_min_frac: DEFAULT_NEUTRAL_LOSS_LABEL_MIN_FRAC,
            isotope_errors: DEFAULT_PRECURSOR_ISOTOPE_ERRORS.to_vec(),
            isotope_errors_explicit: false,
            charge_override: None,
            remove_precursor: false,
            tolerance: DEFAULT_TOLERANCE,
            tolerance_explicit: false,
            normalize: false,
            mz_min: None,
            mz_max: None,
            mode: PlotMode::Auto,
        }
    }
}

#[derive(Clone, Debug)]
struct SvgHeaderLine {
    text: String,
    size: u32,
    color: &'static str,
}

#[derive(Clone, Debug)]
struct PeakLabel {
    observed_mz: f64,
    display_intensity: f64,
    labels: Vec<PeakLabelText>,
}

#[derive(Clone, Debug)]
struct PeakLabelText {
    text: String,
    color: &'static str,
    title: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RulerTickKind {
    Minor,
    Medium,
    Major,
}

#[derive(Clone, Copy, Debug)]
struct RulerTick {
    value: f64,
    kind: RulerTickKind,
}

#[derive(Clone, Debug)]
struct IonTableRow {
    cleavage_index: usize,
    y_ordinal: usize,
    b1: IonTableCell,
    b2: IonTableCell,
    y2: IonTableCell,
    y1: IonTableCell,
}

#[derive(Clone, Debug, Default)]
struct IonTableCell {
    entries: Vec<IonTableEntry>,
}

#[derive(Clone, Debug)]
struct IonTableEntry {
    text: String,
    color: &'static str,
    title: String,
    matched: bool,
}

#[derive(Clone, Copy, Debug)]
struct IonTableLayout {
    width: f64,
    cut_width: f64,
    pad: f64,
    value_columns: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct IonKey {
    series: FragmentSeries,
    cleavage_index: usize,
    charge: u8,
    neutral_loss: Option<NeutralLossKind>,
}

#[derive(Clone, Debug)]
struct PlotPreparedData {
    render_mode: PlotRenderMode,
    bounds: CoordinateRange<f64>,
    y_max: f64,
    points: Vec<(f64, f64)>,
}

#[derive(Clone, Debug, Serialize)]
struct PsmPlotJson {
    mzml: String,
    pepxml: String,
    svg: String,
    scan: PsmScanJson,
    psm: PsmHitJson,
    annotation: PsmAnnotationJson,
    quality: Option<PsmQualityJson>,
    precursor_check: Option<PsmPrecursorCheckJson>,
    warnings: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
struct PsmScanJson {
    scan_id: String,
    index: u32,
    ms_level: u8,
    rt_minutes: Option<f64>,
    precursor_mz: Option<f64>,
    precursor_charge: Option<i32>,
    points: usize,
    base_peak_mz: f64,
    base_peak_intensity: f32,
}

#[derive(Clone, Debug, Serialize)]
struct PsmHitJson {
    rank: usize,
    peptide: String,
    assumed_charge: Option<i32>,
    spectrum: Option<String>,
    start_scan: Option<u64>,
    end_scan: Option<u64>,
    protein: Option<String>,
    calc_neutral_pep_mass: Option<f64>,
    massdiff: Option<f64>,
    scores: Vec<PepXmlScore>,
    modifications: Vec<PepXmlModification>,
}

#[derive(Clone, Debug, Serialize)]
struct PsmAnnotationJson {
    modified_sequence: String,
    charge_context: Option<i32>,
    tolerance: String,
    theoretical_ions: usize,
    matches: usize,
    matched_observed_peaks: usize,
}

#[derive(Clone, Debug, Serialize)]
struct PsmQualityJson {
    snr: f64,
    log2_snr: f64,
    cosine: f64,
    frag_error_mae_ppm: Option<f64>,
    frag_error_mae_da: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
struct PsmPrecursorCheckJson {
    charge: i32,
    monoisotopic_theoretical_mz: f64,
    theoretical_mz: f64,
    observed_mz: f64,
    isotope_error: u8,
    error_da: f64,
    error_ppm: f64,
    within_tolerance: bool,
}

pub fn run(args: Vec<String>) -> anyhow::Result<()> {
    if args
        .iter()
        .any(|arg| matches!(arg.as_str(), "-h" | "--help" | "help"))
    {
        print_plot_help();
        return Ok(());
    }

    let options = parse_plot_args(args)?;
    let input_path = options
        .input_path
        .as_ref()
        .expect("parse_plot_args validates input path");
    let input_format = options
        .input_format
        .expect("parse_plot_args validates input format");

    let spectrum = match input_format {
        PlotInputFormat::Mzml => {
            let selector = options
                .selector
                .as_ref()
                .expect("parse_plot_args validates mzML selector");
            let mut reader = open_reader(input_path.as_path())?;
            load_selected_mzml_spectrum(&mut reader, selector)?
        }
        PlotInputFormat::Ms2 => {
            load_selected_ms2_spectrum(input_path.as_path(), options.selector.as_ref())?
        }
    };
    let prepared = prepare_plot_data(&spectrum, &options);
    let display_charge = options.charge_override.or(spectrum.meta.precursor_charge);
    let neutral_losses = if options.neutral_losses_enabled {
        COMMON_NEUTRAL_LOSSES.as_slice()
    } else {
        &[]
    };

    if let Some(pepxml_path) = options.pepxml_path.as_ref() {
        run_pepxml_plot(
            &options,
            input_path,
            pepxml_path,
            &spectrum,
            &prepared,
            display_charge,
            neutral_losses,
        )?;
        return Ok(());
    }

    let annotation_context = if let Some(peptide_input) = options.peptide_input.as_deref() {
        let mut context = prepare_annotation(
            peptide_input,
            &options.mod_inputs,
            neutral_losses,
            display_charge,
            options.tolerance,
        )?;
        context.isotope_errors = options.isotope_errors.clone();
        Some(context)
    } else {
        None
    };

    let mut warnings = Vec::<String>::new();
    let annotation_report = if let Some(context) = annotation_context.as_ref() {
        if matches!(spectrum.meta.continuity, SignalContinuity::Profile) {
            warnings.push(
                "Fragment annotation skipped because this spectrum is profile-like; rendering observed spectrum only."
                    .to_string(),
            );
            None
        } else {
            let report = annotate_peaks(
                context,
                spectrum.meta.precursor_mz,
                &spectrum.mz,
                &spectrum.intensity,
            );
            if let Some(check) = report.precursor_check.as_ref() {
                if !check.within_tolerance {
                    warnings.push(format!(
                        "Precursor mismatch for {}+: observed {:.4}, theoretical {:.4} ({:+.4} Da, {:+.1} ppm, isotope error {}).",
                        check.charge,
                        check.observed_mz,
                        check.theoretical_mz,
                        check.error_da,
                        check.error_ppm,
                        check.isotope_error,
                    ));
                }
            }
            Some(report)
        }
    } else {
        None
    };

    let svg_path = options.svg_path.clone().unwrap_or_else(|| {
        default_output_path(
            input_path,
            &spectrum,
            annotation_context.as_ref(),
            display_charge,
            options.neutral_losses_enabled,
            options.svg_prefix.as_deref(),
        )
    });
    write_plot_svg(
        &svg_path,
        input_path,
        &spectrum,
        &prepared,
        &options,
        display_charge,
        annotation_context.as_ref(),
        annotation_report.as_ref(),
        &warnings,
    )?;

    println!("Wrote SVG: {}", svg_path.display());
    print_plot_summary(
        &spectrum,
        display_charge,
        annotation_context.as_ref(),
        annotation_report.as_ref(),
        &warnings,
    );

    Ok(())
}

fn prepare_plot_data(spectrum: &LoadedSpectrum, options: &PlotOptions) -> PlotPreparedData {
    let render_mode = resolve_plot_mode(options.mode, spectrum.meta.continuity, spectrum.mz.len());
    let bounds = resolve_bounds(spectrum, options.mz_min, options.mz_max);
    let precursor_exclusion_window = if options.remove_precursor {
        precursor_exclusion_window(spectrum.meta.precursor_mz)
    } else {
        None
    };
    let (points, y_max) = downsample_max_per_bin(
        &spectrum.mz,
        &spectrum.intensity,
        bounds,
        options.normalize,
        precursor_exclusion_window,
        SVG_BINS,
    );
    PlotPreparedData {
        render_mode,
        bounds,
        y_max,
        points,
    }
}

fn run_pepxml_plot(
    options: &PlotOptions,
    input_path: &Path,
    pepxml_path: &Path,
    spectrum: &LoadedSpectrum,
    prepared: &PlotPreparedData,
    display_charge: Option<i32>,
    neutral_losses: &[NeutralLossKind],
) -> anyhow::Result<()> {
    let scan_number = resolve_pepxml_scan_number(spectrum, options.selector.as_ref())?;
    let scan_hits = load_hits_for_scan(pepxml_path, scan_number, options.top_n)?;
    if scan_hits.available_hits == 0 {
        anyhow::bail!(
            "no pepXML search hits found for scan {} in {}",
            scan_number,
            pepxml_path.display()
        );
    }
    if scan_hits.requested_top_n > scan_hits.available_hits {
        eprintln!(
            "warning: requested --top-n {} but pepXML has {} hit(s) for scan {}; plotting available hit(s)",
            scan_hits.requested_top_n, scan_hits.available_hits, scan_number
        );
    }
    if options.svg_path.is_some() && scan_hits.hits.len() > 1 {
        anyhow::bail!("--svg can be used with only one pepXML hit; omit --svg or reduce --top-n");
    }

    for hit in &scan_hits.hits {
        let mut context = prepare_annotation_with_modifications(
            &hit.peptide,
            hit.explicit_modifications(),
            neutral_losses,
            options
                .charge_override
                .or(hit.assumed_charge)
                .or(display_charge),
            options.tolerance,
        )?;
        context.isotope_errors = options.isotope_errors.clone();

        let mut warnings = Vec::<String>::new();
        let annotation_report = build_annotation_report(spectrum, &context, &mut warnings);
        let svg_path = options.svg_path.clone().unwrap_or_else(|| {
            default_pepxml_output_path(
                input_path,
                spectrum,
                &context,
                display_charge.or(context.charge_context),
                options.neutral_losses_enabled,
                options.svg_prefix.as_deref(),
                hit,
            )
        });
        write_plot_svg(
            &svg_path,
            input_path,
            spectrum,
            prepared,
            options,
            display_charge.or(context.charge_context),
            Some(&context),
            annotation_report.as_ref(),
            &warnings,
        )?;
        let json_path = svg_path.with_extension("json");
        write_psm_json(
            &json_path,
            input_path,
            pepxml_path,
            &svg_path,
            spectrum,
            hit,
            &context,
            annotation_report.as_ref(),
            &warnings,
        )?;

        println!("Wrote SVG: {}", svg_path.display());
        println!("Wrote JSON: {}", json_path.display());
        println!(
            "pepXML hit rank {} | peptide {} | charge {}",
            hit.hit_rank,
            context.modified_sequence(),
            context
                .charge_context
                .map(|charge| format!("{charge}+"))
                .unwrap_or_else(|| "-".to_string())
        );
        print_plot_summary(
            spectrum,
            display_charge.or(context.charge_context),
            Some(&context),
            annotation_report.as_ref(),
            &warnings,
        );
    }

    Ok(())
}

fn resolve_pepxml_scan_number(
    spectrum: &LoadedSpectrum,
    selector: Option<&SpectrumSelector>,
) -> anyhow::Result<u64> {
    if let Some(scan) = extract_scan_number(&spectrum.meta.scan_id) {
        return Ok(scan);
    }
    if let Some(SpectrumSelector::ScanNumber(scan)) = selector {
        return Ok(*scan);
    }
    anyhow::bail!(
        "pepXML annotation requires a scan number; use --scan or an mzML native id containing scan=<n>"
    )
}

fn build_annotation_report(
    spectrum: &LoadedSpectrum,
    context: &AnnotationContext,
    warnings: &mut Vec<String>,
) -> Option<AnnotationReport> {
    if matches!(spectrum.meta.continuity, SignalContinuity::Profile) {
        warnings.push(
            "Fragment annotation skipped because this spectrum is profile-like; rendering observed spectrum only."
                .to_string(),
        );
        return None;
    }

    let report = annotate_peaks(
        context,
        spectrum.meta.precursor_mz,
        &spectrum.mz,
        &spectrum.intensity,
    );
    if let Some(check) = report.precursor_check.as_ref() {
        if !check.within_tolerance {
            warnings.push(format!(
                "Precursor mismatch for {}+: observed {:.4}, theoretical {:.4} ({:+.4} Da, {:+.1} ppm, isotope error {}).",
                check.charge,
                check.observed_mz,
                check.theoretical_mz,
                check.error_da,
                check.error_ppm,
                check.isotope_error,
            ));
        }
    }
    Some(report)
}

fn write_plot_svg(
    svg_path: &Path,
    input_path: &Path,
    spectrum: &LoadedSpectrum,
    prepared: &PlotPreparedData,
    options: &PlotOptions,
    display_charge: Option<i32>,
    annotation_context: Option<&AnnotationContext>,
    annotation_report: Option<&AnnotationReport>,
    warnings: &[String],
) -> anyhow::Result<()> {
    if let Some(parent) = svg_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }
    write_spectrum_svg(
        svg_path,
        input_path,
        SVG_WIDTH,
        SVG_HEIGHT,
        spectrum,
        prepared.render_mode,
        prepared.bounds,
        prepared.y_max,
        &prepared.points,
        options.normalize,
        display_charge,
        !options.remove_precursor,
        options.neutral_loss_label_min_frac,
        annotation_context,
        annotation_report,
        warnings,
    )
}

fn print_plot_summary(
    spectrum: &LoadedSpectrum,
    display_charge: Option<i32>,
    annotation_context: Option<&AnnotationContext>,
    annotation_report: Option<&AnnotationReport>,
    warnings: &[String],
) {
    println!(
        "Scan {} (index {}) | ms{} | points {} | precursor {} | base peak {:.4} @ {:.3e}",
        spectrum.meta.scan_id,
        spectrum.meta.idx,
        spectrum.meta.ms_level,
        spectrum.stats.points,
        format_precursor(spectrum.meta.precursor_mz, display_charge),
        spectrum.stats.base_peak_mz,
        spectrum.stats.base_peak_intensity
    );
    if let Some(context) = annotation_context {
        if let Some(report) = annotation_report {
            println!(
                "Annotation: {} theoretical ions, {} matches across {} observed peaks using {}",
                report.fragments.len(),
                report.matches.len(),
                report.matched_peak_count(),
                context.tolerance.label(),
            );
            println!(
                "Quality: {}",
                format_quality_metrics(&report.quality, context.tolerance)
            );
            if let Some(check) = report.precursor_check.as_ref() {
                print_precursor_check(check);
            }
        } else {
            println!(
                "Annotation requested for peptide {} but skipped for this spectrum.",
                context.peptide.sequence()
            );
        }
    }
    for warning in warnings {
        eprintln!("warning: {warning}");
    }
}

fn print_precursor_check(check: &PrecursorCheck) {
    if check.isotope_error == 0 {
        println!(
            "Precursor check: observed {:.4} vs theoretical {:.4} for {}+ ({:+.4} Da, {:+.1} ppm)",
            check.observed_mz, check.theoretical_mz, check.charge, check.error_da, check.error_ppm,
        );
    } else {
        println!(
            "Precursor check: observed {:.4} vs theoretical {:.4} for {}+ ({:+.4} Da, {:+.1} ppm) using isotope error {} (monoisotopic {:.4})",
            check.observed_mz,
            check.theoretical_mz,
            check.charge,
            check.error_da,
            check.error_ppm,
            check.isotope_error,
            check.monoisotopic_theoretical_mz,
        );
    }
}

fn write_psm_json(
    path: &Path,
    input_path: &Path,
    pepxml_path: &Path,
    svg_path: &Path,
    spectrum: &LoadedSpectrum,
    hit: &PepXmlHit,
    context: &AnnotationContext,
    report: Option<&AnnotationReport>,
    warnings: &[String],
) -> anyhow::Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }

    let payload = PsmPlotJson {
        mzml: input_path.display().to_string(),
        pepxml: pepxml_path.display().to_string(),
        svg: svg_path.display().to_string(),
        scan: PsmScanJson {
            scan_id: spectrum.meta.scan_id.clone(),
            index: spectrum.meta.idx,
            ms_level: spectrum.meta.ms_level,
            rt_minutes: spectrum.meta.rt_minutes.map(f64::from),
            precursor_mz: spectrum.meta.precursor_mz,
            precursor_charge: spectrum.meta.precursor_charge,
            points: spectrum.stats.points as usize,
            base_peak_mz: spectrum.stats.base_peak_mz,
            base_peak_intensity: spectrum.stats.base_peak_intensity,
        },
        psm: PsmHitJson {
            rank: hit.hit_rank,
            peptide: hit.peptide.clone(),
            assumed_charge: hit.assumed_charge,
            spectrum: hit.spectrum.clone(),
            start_scan: hit.start_scan,
            end_scan: hit.end_scan,
            protein: hit.protein.clone(),
            calc_neutral_pep_mass: hit.calc_neutral_pep_mass,
            massdiff: hit.massdiff,
            scores: hit.scores.clone(),
            modifications: hit.modifications.clone(),
        },
        annotation: PsmAnnotationJson {
            modified_sequence: context.modified_sequence(),
            charge_context: context.charge_context,
            tolerance: context.tolerance.label(),
            theoretical_ions: report.map(|value| value.fragments.len()).unwrap_or(0),
            matches: report.map(|value| value.matches.len()).unwrap_or(0),
            matched_observed_peaks: report.map(|value| value.matched_peak_count()).unwrap_or(0),
        },
        quality: report.map(|value| PsmQualityJson {
            snr: value.quality.snr_like,
            log2_snr: value.quality.log2_snr_like,
            cosine: value.quality.cosine,
            frag_error_mae_ppm: value.quality.frag_error_mae_ppm,
            frag_error_mae_da: value.quality.frag_error_mae_da,
        }),
        precursor_check: report
            .and_then(|value| value.precursor_check.as_ref())
            .map(|check| PsmPrecursorCheckJson {
                charge: check.charge,
                monoisotopic_theoretical_mz: check.monoisotopic_theoretical_mz,
                theoretical_mz: check.theoretical_mz,
                observed_mz: check.observed_mz,
                isotope_error: check.isotope_error,
                error_da: check.error_da,
                error_ppm: check.error_ppm,
                within_tolerance: check.within_tolerance,
            }),
        warnings: warnings.to_vec(),
    };

    let file =
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    serde_json::to_writer_pretty(file, &payload)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn parse_plot_args(args: Vec<String>) -> anyhow::Result<PlotOptions> {
    let mut options = PlotOptions::default();
    let mut tol_override = None::<MassTolerance>;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                let path = PathBuf::from(iter.next().context("--mzml expects a path")?);
                set_input_source(&mut options, PlotInputFormat::Mzml, path)?;
            }
            "--ms2" => {
                let path = PathBuf::from(iter.next().context("--ms2 expects a path")?);
                set_input_source(&mut options, PlotInputFormat::Ms2, path)?;
            }
            "--index" => {
                let raw = iter.next().context("--index expects an integer")?;
                let idx = raw.parse::<u32>().context("invalid --index")?;
                set_selector(&mut options, SpectrumSelector::Index(idx))?;
            }
            "--scan" => {
                let raw = iter.next().context("--scan expects a scan number")?;
                let scan_number = parse_scan_number_arg(&raw)?;
                set_selector(&mut options, SpectrumSelector::ScanNumber(scan_number))?;
            }
            "--id" | "--native-id" => {
                let id = iter.next().context("--id expects a native id")?;
                set_selector(&mut options, SpectrumSelector::NativeId(id))?;
            }
            "--svg" | "--output" => {
                options.svg_path = Some(PathBuf::from(
                    iter.next().context("--svg expects a file path")?,
                ));
            }
            "--svg-prefix" => {
                options.svg_prefix = Some(
                    iter.next()
                        .context("--svg-prefix expects a short text label")?,
                );
            }
            "--peptide" | "--sequence" | "--sequence-modi" | "--sequencemodi" => {
                let peptide = iter
                    .next()
                    .context("--peptide/--sequence expects a value")?;
                set_peptide_input(&mut options, peptide)?;
            }
            "--pepxml" | "--pep-xml" => {
                if options.pepxml_path.is_some() {
                    anyhow::bail!("specify --pepxml only once");
                }
                options.pepxml_path = Some(PathBuf::from(
                    iter.next().context("--pepxml expects a path")?,
                ));
            }
            "--top-n" => {
                let raw = iter.next().context("--top-n expects an integer")?;
                let value = raw.parse::<usize>().context("invalid --top-n")?;
                if value == 0 {
                    anyhow::bail!("--top-n must be at least 1");
                }
                options.top_n = value;
                options.top_n_explicit = true;
            }
            "--mod" => {
                options
                    .mod_inputs
                    .push(iter.next().context("--mod expects <position>:<delta>")?);
            }
            "--neutral-losses" => {
                options.neutral_losses_enabled = true;
            }
            "--neutral-loss-min-frac" => {
                let raw = iter
                    .next()
                    .context("--neutral-loss-min-frac expects a float between 0 and 1")?;
                let value = raw
                    .parse::<f64>()
                    .context("invalid --neutral-loss-min-frac")?;
                if !(0.0..=1.0).contains(&value) || !value.is_finite() {
                    anyhow::bail!("--neutral-loss-min-frac must be between 0 and 1");
                }
                options.neutral_loss_label_min_frac = value;
            }
            "--isotope-errors" => {
                let raw = iter
                    .next()
                    .context("--isotope-errors expects a comma-separated list like 0,1,2")?;
                options.isotope_errors = parse_isotope_errors(&raw)?;
                options.isotope_errors_explicit = true;
            }
            "--charge" => {
                let raw = iter.next().context("--charge expects an integer")?;
                options.charge_override = Some(raw.parse::<i32>().context("invalid --charge")?);
            }
            "--remove-precursor" => {
                options.remove_precursor = true;
            }
            "--tol-ppm" => {
                let raw = iter.next().context("--tol-ppm expects a float")?;
                let ppm = raw.parse::<f64>().context("invalid --tol-ppm")?;
                if ppm <= 0.0 || !ppm.is_finite() {
                    anyhow::bail!("--tol-ppm must be a positive finite number");
                }
                set_tolerance(&mut tol_override, MassTolerance::Ppm(ppm), "--tol-ppm")?;
            }
            "--tol-da" => {
                let raw = iter.next().context("--tol-da expects a float")?;
                let da = raw.parse::<f64>().context("invalid --tol-da")?;
                if da <= 0.0 || !da.is_finite() {
                    anyhow::bail!("--tol-da must be a positive finite number");
                }
                set_tolerance(&mut tol_override, MassTolerance::Da(da), "--tol-da")?;
            }
            "--normalize" => {
                options.normalize = true;
            }
            "--mz-min" => {
                let raw = iter.next().context("--mz-min expects a float")?;
                options.mz_min = Some(raw.parse::<f64>().context("invalid --mz-min")?);
            }
            "--mz-max" => {
                let raw = iter.next().context("--mz-max expects a float")?;
                options.mz_max = Some(raw.parse::<f64>().context("invalid --mz-max")?);
            }
            "--mode" => {
                let raw = iter.next().context("--mode expects auto|sticks|line")?;
                options.mode = parse_mode(&raw)?;
            }
            other => anyhow::bail!("unknown plot option `{other}`"),
        }
    }

    options.tolerance_explicit = tol_override.is_some();
    options.tolerance = tol_override.unwrap_or(DEFAULT_TOLERANCE);

    if options.input_path.is_none() || options.input_format.is_none() {
        anyhow::bail!("plot requires --mzml <path> or --ms2 <path>");
    }
    if matches!(options.input_format, Some(PlotInputFormat::Mzml)) && options.selector.is_none() {
        anyhow::bail!(
            "plot requires one of --index <n>, --scan <n>, or --id <native-id> for mzML input"
        );
    }
    if let (Some(min), Some(max)) = (options.mz_min, options.mz_max) {
        if min >= max {
            anyhow::bail!("--mz-min must be smaller than --mz-max");
        }
    }
    if options.peptide_input.is_some() && options.pepxml_path.is_some() {
        anyhow::bail!("specify only one of --peptide or --pepxml");
    }
    if options.pepxml_path.is_some() && !matches!(options.input_format, Some(PlotInputFormat::Mzml))
    {
        anyhow::bail!("--pepxml is currently supported only with --mzml input");
    }

    let has_annotation_source = options.peptide_input.is_some() || options.pepxml_path.is_some();
    if !has_annotation_source {
        if !options.mod_inputs.is_empty() {
            anyhow::bail!("--mod requires --peptide (or a sequence alias)");
        }
        if options.neutral_losses_enabled {
            anyhow::bail!("--neutral-losses requires --peptide/--pepxml annotation input");
        }
        if (options.neutral_loss_label_min_frac - DEFAULT_NEUTRAL_LOSS_LABEL_MIN_FRAC).abs()
            > f64::EPSILON
        {
            anyhow::bail!(
                "--neutral-loss-min-frac requires --peptide/--pepxml together with --neutral-losses"
            );
        }
        if options.tolerance_explicit {
            anyhow::bail!("--tol-ppm/--tol-da require --peptide/--pepxml annotation input");
        }
        if options.isotope_errors_explicit {
            anyhow::bail!("--isotope-errors requires --peptide/--pepxml annotation input");
        }
    } else if !options.neutral_losses_enabled
        && (options.neutral_loss_label_min_frac - DEFAULT_NEUTRAL_LOSS_LABEL_MIN_FRAC).abs()
            > f64::EPSILON
    {
        anyhow::bail!("--neutral-loss-min-frac requires --neutral-losses");
    }
    if options.pepxml_path.is_some() && !options.mod_inputs.is_empty() {
        anyhow::bail!("--mod cannot be combined with --pepxml; use pepXML modifications");
    }
    if options.top_n_explicit && options.pepxml_path.is_none() {
        anyhow::bail!("--top-n requires --pepxml");
    }

    Ok(options)
}

fn set_input_source(
    options: &mut PlotOptions,
    input_format: PlotInputFormat,
    input_path: PathBuf,
) -> anyhow::Result<()> {
    if options.input_path.is_some() || options.input_format.is_some() {
        anyhow::bail!("specify only one of --mzml or --ms2");
    }
    options.input_path = Some(input_path);
    options.input_format = Some(input_format);
    Ok(())
}

fn set_selector(options: &mut PlotOptions, selector: SpectrumSelector) -> anyhow::Result<()> {
    if options.selector.is_some() {
        anyhow::bail!("specify only one of --index, --scan, or --id");
    }
    options.selector = Some(selector);
    Ok(())
}

fn set_peptide_input(options: &mut PlotOptions, peptide: String) -> anyhow::Result<()> {
    if options.peptide_input.is_some() {
        anyhow::bail!("specify peptide input only once");
    }
    options.peptide_input = Some(peptide);
    Ok(())
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

fn parse_scan_number_arg(raw: &str) -> anyhow::Result<u64> {
    if let Ok(value) = raw.parse::<u64>() {
        return Ok(value);
    }
    extract_scan_number(raw)
        .ok_or_else(|| anyhow::anyhow!("invalid --scan `{raw}` (expected 107468 or scan=107468)"))
}

fn parse_isotope_errors(raw: &str) -> anyhow::Result<Vec<u8>> {
    let mut values = raw
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(|token| {
            token
                .parse::<u8>()
                .with_context(|| format!("invalid isotope error `{token}` in `{raw}`"))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    if values.is_empty() {
        anyhow::bail!("--isotope-errors requires at least one integer");
    }
    values.sort_unstable();
    values.dedup();
    Ok(values)
}

fn parse_mode(raw: &str) -> anyhow::Result<PlotMode> {
    match raw {
        "auto" => Ok(PlotMode::Auto),
        "sticks" | "stick" => Ok(PlotMode::Sticks),
        "line" => Ok(PlotMode::Line),
        other => anyhow::bail!("invalid --mode `{other}` (expected auto|sticks|line)"),
    }
}

fn print_plot_help() {
    let program = crate::program_name();
    println!("{program} plot");
    println!();
    println!("USAGE:");
    println!(
        "  {program} plot (--mzml <file> (--index <n> | --scan <n> | --id <native-id>) | --ms2 <file> [--index <n> | --scan <n> | --id <scan-id>]) [options]"
    );
    println!();
    println!("OPTIONS:");
    println!(
        "  --mzml <file>                Input mzML file
  --ms2 <file>                 Input plain-text MS2 file; defaults to first spectrum"
    );
    println!("  --index <n>                  Zero-based spectrum index");
    println!("  --scan <n>                   Scan number, e.g. 4821 or scan=4821");
    println!("  --id <native-id>             Full or partial native id");
    println!("  --svg <file>                 Output SVG path (default: exports/auto-name.svg)");
    println!("  --svg-prefix <text>          Prefix for autogenerated SVG names, e.g. calibrated");
    println!("  --peptide <SEQ>              Preferred peptide input; supports M[+15.9949], [+42.0106]PEPTIDE, /charge");
    println!("  --sequence <SEQ>             Alias for --peptide");
    println!("  --sequence-modi <SEQ>        Alias for --peptide");
    println!("  --sequencemodi <SEQ>         Alias for --peptide");
    println!("  --pepxml <file>              Annotate selected mzML scan from pepXML/pepXML.gz hits; mutually exclusive with --peptide");
    println!("  --top-n <n>                  Number of pepXML hits to plot for the selected scan (default: 1)");
    println!("  --mod <position>:<delta>     Repeatable explicit modification, 1-based");
    println!(
        "  --neutral-losses             Enable residue-aware -H2O / -NH3 and phospho -H3PO4 fragment variants"
    );
    println!("  --neutral-loss-min-frac <f>  Label neutral losses only above this base-peak fraction (default: 0.03)");
    println!("  --isotope-errors <list>      Allowed precursor isotope errors, e.g. 0,1,2 (default: 0,1,2)");
    println!("  --charge <int>               Optional precursor charge override");
    println!("  --remove-precursor           Hide the precursor guide and omit peaks within +/-1.5 Da of precursor m/z");
    println!("  --tol-ppm <ppm>              Fragment tolerance in ppm (default: 20)");
    println!("  --tol-da <da>                Fragment tolerance in Daltons");
    println!("  --normalize                  Normalize intensities to base peak = 1");
    println!("  --mode <auto|sticks|line>    Plot rendering mode");
    println!("  --mz-min <float>             Left x-bound");
    println!("  --mz-max <float>             Right x-bound");
    println!("  --help                       Show this help");
    println!();
    println!("EXAMPLES:");
    println!("  {program} plot --mzml sample.mzML --index 4821");
    println!("  {program} plot --mzml sample.mzML --scan 4821 --svg out.svg");
    println!(
        "  {program} plot --mzml sample.mzML --scan 4821 --peptide DSAVYFCARTKILDFD --tol-da 0.5"
    );
    println!(
        "  {program} plot --mzml sample.mzML --index 4821 --peptide DSAVYFCARTKILDFD --mod 7:+57.021464"
    );
    println!("  {program} plot --mzml sample.mzML --scan 4821 --pepxml search.pep.xml --top-n 3");
    println!(
        "  {program} plot --mzml sample.mzML --scan 4821 --peptide DSAVYFCARTKILDFD --neutral-losses --neutral-loss-min-frac 0.05 --svg-prefix calibrated"
    );
    println!(
        "  {program} plot --mzml sample.mzML --scan 3079 --peptide [+304.2071]T[+79.9663]S[+79.9663]SSSPSR/3"
    );
    println!(
        "  {program} plot --mzml sample.mzML --scan 3480 --peptide [+304.2071]T[+79.9663]S[+79.9663]SSSPSR/3 --isotope-errors 0,1,2"
    );
}

fn default_output_path(
    input_path: &Path,
    spectrum: &LoadedSpectrum,
    annotation_context: Option<&AnnotationContext>,
    display_charge: Option<i32>,
    neutral_losses_enabled: bool,
    svg_prefix: Option<&str>,
) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let mut parts = Vec::<String>::new();
    if let Some(prefix) = svg_prefix.and_then(sanitize_filename_label) {
        parts.push(prefix);
    }

    let source_stem = input_path
        .file_stem()
        .and_then(|value| value.to_str())
        .and_then(sanitize_filename_label)
        .unwrap_or_else(|| "spectrum".to_string());
    parts.push(source_stem);

    let scan_component = extract_scan_number(&spectrum.meta.scan_id)
        .map(|scan| format!("scan{scan}"))
        .unwrap_or_else(|| format!("index{}", spectrum.meta.idx));
    parts.push(scan_component);

    if let Some(context) = annotation_context {
        parts.push(default_peptide_label(
            context.peptide.sequence(),
            display_charge.or(context.charge_context),
        ));
    }

    parts.push(if neutral_losses_enabled {
        "nl-on".to_string()
    } else {
        "nl-off".to_string()
    });
    parts.push(format!("ms{}", spectrum.meta.ms_level));
    parts.push(ts.to_string());

    let filename = format!("{}.svg", parts.join("__"));
    PathBuf::from("exports").join(filename)
}

fn default_pepxml_output_path(
    input_path: &Path,
    spectrum: &LoadedSpectrum,
    annotation_context: &AnnotationContext,
    display_charge: Option<i32>,
    neutral_losses_enabled: bool,
    svg_prefix: Option<&str>,
    hit: &PepXmlHit,
) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let mut parts = Vec::<String>::new();
    if let Some(prefix) = svg_prefix.and_then(sanitize_filename_label) {
        parts.push(prefix);
    }

    let source_stem = input_path
        .file_stem()
        .and_then(|value| value.to_str())
        .and_then(sanitize_filename_label)
        .unwrap_or_else(|| "spectrum".to_string());
    parts.push(source_stem);

    let scan_component = extract_scan_number(&spectrum.meta.scan_id)
        .map(|scan| format!("scan{scan}"))
        .unwrap_or_else(|| format!("index{}", spectrum.meta.idx));
    parts.push(scan_component);
    parts.push(format!("rank{}", hit.hit_rank));
    parts.push(default_peptide_label(
        annotation_context.peptide.sequence(),
        display_charge.or(annotation_context.charge_context),
    ));
    parts.push(if neutral_losses_enabled {
        "nl-on".to_string()
    } else {
        "nl-off".to_string()
    });
    parts.push(format!("ms{}", spectrum.meta.ms_level));
    parts.push(ts.to_string());

    let filename = format!("{}.svg", parts.join("__"));
    PathBuf::from("exports").join(filename)
}

fn format_precursor(precursor_mz: Option<f64>, precursor_charge: Option<i32>) -> String {
    match (precursor_mz, precursor_charge) {
        (Some(mz), Some(charge)) => format!("{mz:.4} ({charge}+)"),
        (Some(mz), None) => format!("{mz:.4}"),
        (None, Some(charge)) => format!("? ({charge}+)"),
        (None, None) => "-".to_string(),
    }
}

fn format_quality_metrics(metrics: &AnnotationQualityMetrics, tolerance: MassTolerance) -> String {
    let frag_error = format_fragment_error_mae(metrics, tolerance);
    format!(
        "SNR={:.3} | log2_SNR={:.3} | cosine={:.3} | {}",
        metrics.snr_like, metrics.log2_snr_like, metrics.cosine, frag_error
    )
}

fn format_fragment_error_mae(
    metrics: &AnnotationQualityMetrics,
    tolerance: MassTolerance,
) -> String {
    match tolerance {
        MassTolerance::Da(da) if da >= 0.25 => metrics
            .frag_error_mae_da
            .map(|value| format!("frag_error_mae_da={value:.4}"))
            .unwrap_or_else(|| "frag_error_mae_da=NA".to_string()),
        _ => metrics
            .frag_error_mae_ppm
            .map(|value| format!("frag_error_mae_ppm={value:.2}"))
            .unwrap_or_else(|| "frag_error_mae_ppm=NA".to_string()),
    }
}

fn resolve_plot_mode(
    mode: PlotMode,
    continuity: SignalContinuity,
    points: usize,
) -> PlotRenderMode {
    match mode {
        PlotMode::Auto => match continuity {
            SignalContinuity::Profile => PlotRenderMode::Line,
            SignalContinuity::Unknown if points > 5_000 => PlotRenderMode::Line,
            _ => PlotRenderMode::Sticks,
        },
        PlotMode::Sticks => PlotRenderMode::Sticks,
        PlotMode::Line => PlotRenderMode::Line,
    }
}

fn resolve_bounds(
    spectrum: &LoadedSpectrum,
    mz_min_override: Option<f64>,
    mz_max_override: Option<f64>,
) -> CoordinateRange<f64> {
    let mut min_x = mz_min_override.unwrap_or(spectrum.stats.mz_min);
    let mut max_x = mz_max_override.unwrap_or(spectrum.stats.mz_max);
    if !min_x.is_finite() {
        min_x = 0.0;
    }
    if !max_x.is_finite() {
        max_x = min_x + 1.0;
    }
    if min_x >= max_x {
        max_x = min_x + 1.0;
    }
    CoordinateRange::new(min_x, max_x)
}

fn downsample_max_per_bin(
    mz: &[f64],
    intensity: &[f32],
    bounds: CoordinateRange<f64>,
    normalize: bool,
    excluded_window: Option<(f64, f64)>,
    bins: usize,
) -> (Vec<(f64, f64)>, f64) {
    let mut min_x = bounds.start;
    let mut max_x = bounds.end;
    if min_x > max_x {
        std::mem::swap(&mut min_x, &mut max_x);
    }

    let span = (max_x - min_x).max(1e-9);
    let bins = bins.clamp(16, 100_000);

    let mut best_y = vec![f64::NEG_INFINITY; bins];
    let mut best_x = vec![0.0; bins];
    let mut has = vec![false; bins];

    for (&mz, &inten) in mz.iter().zip(intensity.iter()) {
        if mz < min_x || mz > max_x {
            continue;
        }
        if let Some((excluded_min, excluded_max)) = excluded_window {
            if mz >= excluded_min && mz <= excluded_max {
                continue;
            }
        }
        let inten = inten as f64;
        if !inten.is_finite() {
            continue;
        }

        let frac = ((mz - min_x) / span).clamp(0.0, 1.0);
        let mut bin = (frac * bins as f64) as usize;
        if bin >= bins {
            bin = bins - 1;
        }

        if !has[bin] || inten > best_y[bin] {
            has[bin] = true;
            best_y[bin] = inten;
            best_x[bin] = mz;
        }
    }

    let mut points = has
        .iter()
        .enumerate()
        .filter_map(|(idx, present)| {
            if *present {
                Some((best_x[idx], best_y[idx]))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    if points.is_empty() {
        return (points, 1.0);
    }

    points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    let max_intensity = points
        .iter()
        .fold(0.0_f64, |acc, (_, value)| acc.max(*value))
        .max(1e-6_f64);

    if normalize && max_intensity > 0.0 {
        for (_, value) in &mut points {
            *value /= max_intensity;
        }
        (points, 1.1)
    } else {
        (points, max_intensity * 1.1)
    }
}

fn precursor_exclusion_window(precursor_mz: Option<f64>) -> Option<(f64, f64)> {
    precursor_mz.filter(|value| value.is_finite()).map(|value| {
        (
            value - PRECURSOR_REMOVAL_HALF_WIDTH_DA,
            value + PRECURSOR_REMOVAL_HALF_WIDTH_DA,
        )
    })
}

fn write_spectrum_svg(
    path: &Path,
    source_path: &Path,
    width: u32,
    height: u32,
    spectrum: &LoadedSpectrum,
    mode: PlotRenderMode,
    x_bounds: CoordinateRange<f64>,
    y_max: f64,
    points: &[(f64, f64)],
    normalize: bool,
    display_charge: Option<i32>,
    show_precursor_marker: bool,
    neutral_loss_label_min_frac: f64,
    annotation_context: Option<&AnnotationContext>,
    annotation_report: Option<&AnnotationReport>,
    warnings: &[String],
) -> anyhow::Result<()> {
    let mut file =
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?;

    let header_lines = build_header_lines(
        source_path,
        spectrum,
        normalize,
        display_charge,
        neutral_loss_label_min_frac,
        annotation_context,
        annotation_report,
        warnings,
    );
    let ladder_height = if annotation_report.is_some() {
        122.0
    } else {
        0.0
    };
    let ion_table_rows = annotation_report.map(build_ion_table_rows);

    let margin_left = 116.0;
    let margin_right = 28.0;
    let margin_bottom = 122.0;
    let line_height = 24.0;
    let header_top = 34.0;
    let header_height = header_lines.len() as f64 * line_height;
    let ladder_top = header_top + header_height + 14.0;
    let plot_top = ladder_top + ladder_height + if ladder_height > 0.0 { 16.0 } else { 0.0 };

    let base_w = width as f64;
    let base_h = height as f64;
    let plot_w = (base_w - margin_left - margin_right).max(1.0);
    let ion_table_width = ion_table_rows
        .as_ref()
        .map(|rows| estimate_ion_table_width(rows, base_w))
        .unwrap_or(0.0);
    let extra_width = if ion_table_rows.is_some() {
        ION_TABLE_GAP + ion_table_width
    } else {
        0.0
    };
    let table_required_height = ion_table_rows
        .as_ref()
        .map(|rows| ion_table_content_height(rows.len()))
        .unwrap_or(0.0);
    let plot_h = (base_h - plot_top - margin_bottom)
        .max(table_required_height)
        .max(1.0);
    let total_w = base_w + extra_width;
    let total_h = plot_top + plot_h + margin_bottom;
    let svg_w = total_w.ceil() as u32;
    let svg_h = total_h.ceil() as u32;

    let y_span = y_max.max(1e-9);
    let plot_canvas = SvgCanvas::new(
        margin_left,
        plot_top,
        plot_w,
        plot_h,
        x_bounds,
        CoordinateRange::new(0.0, y_span),
    );
    let x_axis = AxisProps::new(AxisOrientation::Bottom, "m/z");
    let y_axis =
        AxisProps::new(AxisOrientation::Left, "Intensity").with_tick_label_style(if normalize {
            AxisTickLabelStyle::Precision(2)
        } else {
            AxisTickLabelStyle::Scientific(2)
        });

    writeln!(
        file,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"##,
        width = svg_w,
        height = svg_h,
    )?;
    writeln!(
        file,
        r##"<rect x="0" y="0" width="{width}" height="{height}" fill="#fbfcfe"/>"##,
        width = svg_w,
        height = svg_h,
    )?;
    writeln!(
        file,
        r##"<rect x="10" y="10" width="{w}" height="{h}" rx="12" fill="white" stroke="{stroke}" stroke-width="1"/>"##,
        w = svg_w.saturating_sub(20),
        h = svg_h.saturating_sub(20),
        stroke = COLOR_CARD_BORDER,
    )?;

    for (idx, line) in header_lines.iter().enumerate() {
        writeln!(
            file,
            r##"<text x="{x}" y="{y}" font-family="Helvetica, Arial, sans-serif" font-size="{size}" fill="{fill}">{text}</text>"##,
            x = margin_left,
            y = header_top + idx as f64 * line_height,
            size = line.size,
            fill = line.color,
            text = escape_xml(&line.text),
        )?;
    }

    let table_left = margin_left + plot_w + (extra_width - ion_table_width);
    if let Some(report) = annotation_report {
        draw_ladder(
            &mut file,
            margin_left,
            ladder_top,
            plot_w,
            ladder_height,
            report,
            spectrum.stats.base_peak_intensity,
            neutral_loss_label_min_frac,
        )?;
    }

    writeln!(
        file,
        r##"<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="none" stroke="{stroke}" stroke-width="1"/>"##,
        x = margin_left,
        y = plot_top,
        w = plot_w,
        h = plot_h,
        stroke = COLOR_AXIS,
    )?;

    for tick in 0..=4 {
        let frac = tick as f64 / 4.0;
        let y = plot_canvas.top() + frac * plot_canvas.height();
        let value = (1.0 - frac) * y_span;
        writeln!(
            file,
            r##"<line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="#e5eaf0" stroke-width="1"/>"##,
            x1 = plot_canvas.left(),
            x2 = plot_canvas.right(),
            y = y,
        )?;
        writeln!(
            file,
            r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="end" dominant-baseline="middle">{text}</text>"##,
            x = plot_canvas.left() - 12.0,
            y = y,
            font = PLOT_TICK_FONT,
            fill = COLOR_SUBTLE,
            text = y_axis.format_tick(value),
        )?;
    }

    let ruler_ticks = build_ruler_ticks(x_bounds.min(), x_bounds.max(), plot_canvas.width());
    let major_ruler_step = ruler_major_step(x_bounds.min(), x_bounds.max(), plot_canvas.width());
    let x_label_y = plot_canvas.bottom() + 36.0;
    let x_axis_title_y = plot_canvas.bottom() + 78.0;
    for tick in &ruler_ticks {
        let px = plot_canvas.x(tick.value);
        let tick_len = match tick.kind {
            RulerTickKind::Minor => 4.0,
            RulerTickKind::Medium => 7.0,
            RulerTickKind::Major => 11.0,
        };
        let tick_stroke = match tick.kind {
            RulerTickKind::Minor => "#8fa1b5",
            RulerTickKind::Medium => "#64748b",
            RulerTickKind::Major => COLOR_AXIS,
        };
        let tick_width = match tick.kind {
            RulerTickKind::Minor => 0.8,
            RulerTickKind::Medium => 1.0,
            RulerTickKind::Major => 1.2,
        };
        if matches!(tick.kind, RulerTickKind::Major) {
            writeln!(
                file,
                r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="#edf1f5" stroke-width="1"/>"##,
                x = px,
                y1 = plot_canvas.top(),
                y2 = plot_canvas.bottom(),
            )?;
        }
        writeln!(
            file,
            r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="{width:.1}"/>"##,
            x = px,
            y1 = plot_canvas.bottom(),
            y2 = plot_canvas.bottom() + tick_len,
            stroke = tick_stroke,
            width = tick_width,
        )?;
        if matches!(tick.kind, RulerTickKind::Major) {
            writeln!(
                file,
                r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle">{text}</text>"##,
                x = px,
                y = x_label_y,
                font = PLOT_TICK_FONT,
                fill = COLOR_AXIS,
                text = format_ruler_label(tick.value, major_ruler_step),
            )?;
        }
    }

    if show_precursor_marker {
        if let Some(precursor_mz) = spectrum.meta.precursor_mz {
            if precursor_mz >= x_bounds.min() && precursor_mz <= x_bounds.max() {
                let px = plot_canvas.x(precursor_mz);
                writeln!(
                    file,
                    r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="1.5" stroke-opacity="0.45"/>"##,
                    x = px,
                    y1 = plot_top,
                    y2 = plot_top + plot_h,
                    stroke = COLOR_WARNING,
                )?;
                writeln!(
                    file,
                    r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle">precursor {text}</text>"##,
                    x = px,
                    y = plot_canvas.top() - 8.0,
                    font = PLOT_PRECURSOR_LABEL_FONT,
                    fill = COLOR_WARNING,
                    text = format!("{precursor_mz:.4}"),
                )?;
            }
        }
    }

    match mode {
        PlotRenderMode::Sticks => {
            for &(x, y) in points {
                let (px, py) = plot_canvas.transform(x, y);
                let (px0, py0) = plot_canvas.transform(x, 0.0);
                writeln!(
                    file,
                    r##"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="1"/>"##,
                    x1 = px0,
                    y1 = py0,
                    x2 = px,
                    y2 = py,
                    stroke = COLOR_PLOT,
                )?;
            }
        }
        PlotRenderMode::Line => {
            let mut d = String::new();
            for (idx, &(x, y)) in points.iter().enumerate() {
                let (px, py) = plot_canvas.transform(x, y);
                if idx == 0 {
                    d.push_str(&format!("M{px:.2},{py:.2}"));
                } else {
                    d.push_str(&format!(" L{px:.2},{py:.2}"));
                }
            }
            writeln!(
                file,
                r##"<path d="{d}" fill="none" stroke="{stroke}" stroke-width="1.2"/>"##,
                d = d,
                stroke = COLOR_PLOT,
            )?;
        }
    }

    if let Some(report) = annotation_report {
        for peak_label in collect_peak_labels(
            report,
            spectrum,
            normalize,
            x_bounds,
            neutral_loss_label_min_frac,
        ) {
            let (px, py) =
                plot_canvas.transform(peak_label.observed_mz, peak_label.display_intensity);
            for (rank, label) in peak_label.labels.iter().enumerate() {
                let text_y = (py - 10.0 - rank as f64 * 12.0).max(plot_canvas.top() + 12.0);
                let use_end_anchor = px > plot_canvas.right() - 90.0;
                let text_x = if use_end_anchor { px - 4.0 } else { px + 4.0 };
                let anchor = if use_end_anchor { "end" } else { "start" };
                writeln!(
                    file,
                    r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="0.8" stroke-opacity="0.6"/>"##,
                    x = px,
                    y1 = py,
                    y2 = text_y + 2.0,
                    stroke = label.color,
                )?;
                writeln!(
                    file,
                    r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="{anchor}"><title>{title}</title>{text}</text>"##,
                    x = text_x,
                    y = text_y,
                    font = PLOT_PEAK_LABEL_FONT,
                    fill = label.color,
                    anchor = anchor,
                    title = escape_xml(&label.title),
                    text = escape_xml(&label.text),
                )?;
            }
        }
    }

    if let (Some(report), Some(rows)) = (annotation_report, ion_table_rows.as_ref()) {
        draw_ion_table(
            &mut file,
            table_left,
            plot_top,
            ion_table_width,
            plot_h,
            report,
            rows,
        )?;
    }
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle">{text}</text>"##,
        x = plot_canvas.left() + plot_canvas.width() / 2.0,
        y = x_axis_title_y,
        font = PLOT_AXIS_LABEL_FONT,
        fill = COLOR_AXIS,
        text = x_axis.label(),
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle" transform="rotate(-90 {x:.2} {y:.2})">{text}</text>"##,
        x = plot_canvas.left() - 88.0,
        y = plot_canvas.top() + plot_canvas.height() / 2.0,
        font = PLOT_AXIS_LABEL_FONT,
        fill = COLOR_AXIS,
        text = y_axis.label(),
    )?;

    writeln!(file, "</svg>")?;
    Ok(())
}

fn build_ruler_ticks(min_x: f64, max_x: f64, plot_width: f64) -> Vec<RulerTick> {
    let minor_step = choose_ruler_minor_step(min_x, max_x, plot_width);
    let major_step = minor_step * 10.0;
    let medium_step = minor_step * 5.0;
    let start = (min_x / minor_step).ceil() * minor_step;
    let end = (max_x / minor_step).floor() * minor_step;

    let mut ticks = Vec::new();
    let mut value = start;
    let mut guard = 0usize;
    while value <= end + minor_step * 0.5 && guard < 4096 {
        if value >= min_x - minor_step * 1e-6 && value <= max_x + minor_step * 1e-6 {
            let kind = if is_multiple_of_step(value, major_step) {
                RulerTickKind::Major
            } else if is_multiple_of_step(value, medium_step) {
                RulerTickKind::Medium
            } else {
                RulerTickKind::Minor
            };
            ticks.push(RulerTick { value, kind });
        }
        value += minor_step;
        guard += 1;
    }
    ticks
}

fn choose_ruler_minor_step(min_x: f64, max_x: f64, plot_width: f64) -> f64 {
    let span = (max_x - min_x).abs();
    if !span.is_finite() || span <= f64::EPSILON || !plot_width.is_finite() || plot_width <= 0.0 {
        return 10.0;
    }

    let px_per_unit = plot_width / span;
    let candidates = [
        0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 25.0, 50.0, 100.0, 200.0, 250.0, 500.0,
    ];
    candidates
        .into_iter()
        .find(|step| *step * px_per_unit >= 8.0)
        .unwrap_or(1000.0)
}

fn ruler_major_step(min_x: f64, max_x: f64, plot_width: f64) -> f64 {
    choose_ruler_minor_step(min_x, max_x, plot_width) * 10.0
}

fn is_multiple_of_step(value: f64, step: f64) -> bool {
    if !step.is_finite() || step <= 0.0 {
        return false;
    }
    let nearest = (value / step).round() * step;
    (nearest - value).abs() <= step * 1e-6
}

fn format_ruler_label(value: f64, major_step: f64) -> String {
    let decimals = tick_decimals(major_step);
    format!("{value:.decimals$}")
}

fn tick_decimals(step: f64) -> usize {
    if !step.is_finite() || step <= 0.0 {
        return 0;
    }
    let exponent = step.abs().log10().floor() as i32;
    if exponent >= 0 {
        0
    } else {
        (-exponent) as usize
    }
}

fn build_header_lines(
    source_path: &Path,
    spectrum: &LoadedSpectrum,
    normalize: bool,
    display_charge: Option<i32>,
    neutral_loss_label_min_frac: f64,
    annotation_context: Option<&AnnotationContext>,
    annotation_report: Option<&AnnotationReport>,
    warnings: &[String],
) -> Vec<SvgHeaderLine> {
    let mut lines = Vec::new();
    let source_label = source_path
        .file_name()
        .and_then(|value| value.to_str())
        .map(|value| value.to_string())
        .unwrap_or_else(|| source_path.display().to_string());
    lines.push(SvgHeaderLine {
        text: format!(
            "Scan {} | index {} | ms{} | rt={} min | {} | {}",
            spectrum.meta.scan_id,
            spectrum.meta.idx,
            spectrum.meta.ms_level,
            spectrum
                .meta
                .rt_minutes
                .map(|value| format!("{value:.2}"))
                .unwrap_or_else(|| "-".to_string()),
            spectrum.meta.continuity,
            if normalize { "normalized" } else { "raw" },
        ),
        size: PLOT_HEADER_TITLE_FONT,
        color: COLOR_TEXT,
    });
    lines.push(SvgHeaderLine {
        text: format!("Source: {source_label}"),
        size: PLOT_HEADER_META_FONT,
        color: COLOR_SUBTLE,
    });
    lines.push(SvgHeaderLine {
        text: format!(
            "Precursor: {} | points: {} | base peak: {:.4} @ {:.3e}",
            format_precursor(spectrum.meta.precursor_mz, display_charge),
            spectrum.stats.points,
            spectrum.stats.base_peak_mz,
            spectrum.stats.base_peak_intensity
        ),
        size: PLOT_HEADER_DETAIL_FONT,
        color: COLOR_TEXT,
    });

    if let Some(context) = annotation_context {
        let mut text = format!(
            "Peptide: {} | tolerance: {}",
            context.modified_sequence(),
            context.tolerance.label()
        );
        if let Some(charge) = context.charge_context {
            text.push_str(&format!(" | charge context: {charge}+"));
        }
        if let Some(report) = annotation_report {
            text.push_str(&format!(
                " | matched {} / {} ions across {} peaks",
                report.matches.len(),
                report.fragments.len(),
                report.matched_peak_count()
            ));
        }
        lines.push(SvgHeaderLine {
            text,
            size: PLOT_HEADER_DETAIL_FONT,
            color: COLOR_TEXT,
        });

        if let Some(report) = annotation_report {
            lines.push(SvgHeaderLine {
                text: format!(
                    "Quality: {}",
                    format_quality_metrics(&report.quality, context.tolerance)
                ),
                size: PLOT_HEADER_META_FONT,
                color: COLOR_SUBTLE,
            });
        }

        if let Some(mod_label) = context.modifications_label() {
            lines.push(SvgHeaderLine {
                text: format!("Applied mods: {mod_label}"),
                size: PLOT_HEADER_META_FONT,
                color: COLOR_SUBTLE,
            });
        }
        if let Some(isotope_errors) = context.isotope_errors_label() {
            lines.push(SvgHeaderLine {
                text: format!("Precursor isotope errors: {isotope_errors}"),
                size: PLOT_HEADER_META_FONT,
                color: COLOR_SUBTLE,
            });
        }
        if let Some(loss_label) = context.neutral_losses_label() {
            lines.push(SvgHeaderLine {
                text: format!(
                    "Neutral losses: {loss_label} (residue-aware, labels >= {:.1}% bp)",
                    neutral_loss_label_min_frac * 100.0
                ),
                size: PLOT_HEADER_META_FONT,
                color: COLOR_SUBTLE,
            });
        }
    }

    for warning in warnings {
        lines.push(SvgHeaderLine {
            text: warning.clone(),
            size: PLOT_HEADER_META_FONT,
            color: COLOR_WARNING,
        });
    }

    lines
}

fn draw_ladder(
    file: &mut File,
    left: f64,
    top: f64,
    width: f64,
    height: f64,
    report: &AnnotationReport,
    base_peak_intensity: f32,
    neutral_loss_label_min_frac: f64,
) -> anyhow::Result<()> {
    writeln!(
        file,
        r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" rx="8" fill="#f8fafc" stroke="#dbe2ea" stroke-width="1"/>"##,
        x = left,
        y = top,
        w = width,
        h = height,
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}">Matched fragment ladder</text>"##,
        x = left + 12.0,
        y = top + 20.0,
        font = LADDER_TITLE_FONT,
        fill = COLOR_SUBTLE,
    )?;

    let residues = report.context.peptide.residue_chars();
    if residues.len() < 2 {
        return Ok(());
    }

    let y_mid = top + height / 2.0 + 8.0;
    let b_label_y = top + 42.0;
    let y_label_y = top + height - 18.0;
    let cleavage_index_y = y_mid + 4.0;

    for (idx, residue) in residues.iter().enumerate() {
        let x = left + width * ((idx as f64) + 0.5) / (residues.len() as f64);
        writeln!(
            file,
            r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle">{text}</text>"##,
            x = x,
            y = y_mid,
            font = LADDER_RESIDUE_FONT,
            fill = COLOR_TEXT,
            text = residue,
        )?;
    }

    let site_labels = build_site_labels(report, base_peak_intensity, neutral_loss_label_min_frac);
    for cleavage_index in 1..residues.len() {
        let x = left + width * (cleavage_index as f64) / (residues.len() as f64);
        let labels = site_labels
            .get(&cleavage_index)
            .cloned()
            .unwrap_or_else(SiteLabels::default);
        let top_color = if labels.b.is_empty() {
            COLOR_SERIES_NEUTRAL
        } else {
            COLOR_SERIES_B
        };
        let bottom_color = if labels.y.is_empty() {
            COLOR_SERIES_NEUTRAL
        } else {
            COLOR_SERIES_Y
        };
        writeln!(
            file,
            r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="1.6"/>"##,
            x = x,
            y1 = y_mid - 6.0,
            y2 = y_mid - 18.0,
            stroke = top_color,
        )?;
        writeln!(
            file,
            r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="1.6"/>"##,
            x = x,
            y1 = y_mid + 6.0,
            y2 = y_mid + 18.0,
            stroke = bottom_color,
        )?;
        if should_render_ladder_index(cleavage_index, residues.len()) {
            writeln!(
                file,
                r##"<text class="ladder-index" x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle"><title>{title}</title>{text}</text>"##,
                x = x,
                y = cleavage_index_y,
                font = LADDER_INDEX_FONT,
                fill = COLOR_SUBTLE,
                title = escape_xml(&format!(
                    "Cleavage after residue {} / before residue {}",
                    cleavage_index,
                    cleavage_index + 1
                )),
                text = cleavage_index,
            )?;
        }

        if !labels.b.is_empty() {
            writeln!(
                file,
                r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle"><title>{title}</title>{text}</text>"##,
                x = x,
                y = b_label_y,
                font = LADDER_ION_FONT,
                fill = COLOR_SERIES_B,
                title = escape_xml(&join_site_label_titles(&labels.b)),
                text = escape_xml(&join_site_label_text(&labels.b)),
            )?;
        }
        if !labels.y.is_empty() {
            writeln!(
                file,
                r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}" text-anchor="middle"><title>{title}</title>{text}</text>"##,
                x = x,
                y = y_label_y,
                font = LADDER_ION_FONT,
                fill = COLOR_SERIES_Y,
                title = escape_xml(&join_site_label_titles(&labels.y)),
                text = escape_xml(&join_site_label_text(&labels.y)),
            )?;
        }
    }

    Ok(())
}

#[derive(Clone, Debug, Default)]
struct SiteLabels {
    b: Vec<SiteLabel>,
    y: Vec<SiteLabel>,
}

#[derive(Clone, Debug)]
struct SiteLabel {
    text: String,
    title: String,
}

fn build_site_labels(
    report: &AnnotationReport,
    base_peak_intensity: f32,
    neutral_loss_label_min_frac: f64,
) -> BTreeMap<usize, SiteLabels> {
    let mut labels = BTreeMap::<usize, SiteLabels>::new();
    for cleavage_index in 1..report.context.peptide.len() {
        labels.entry(cleavage_index).or_default();
    }
    let mut seen = HashSet::<(usize, FragmentSeries, String)>::new();
    for matched in &report.matches {
        if !should_render_match_label(matched, base_peak_intensity, neutral_loss_label_min_frac) {
            continue;
        }
        let entry = labels.entry(matched.fragment.cleavage_index).or_default();
        let label = matched.fragment.label();
        let key = (
            matched.fragment.cleavage_index,
            matched.fragment.series,
            label.clone(),
        );
        if !seen.insert(key) {
            continue;
        }
        let site_label = SiteLabel {
            text: label,
            title: format_match_title(matched),
        };
        match matched.fragment.series {
            FragmentSeries::B => entry.b.push(site_label),
            FragmentSeries::Y => entry.y.push(site_label),
        }
    }
    for site in labels.values_mut() {
        site.b.sort_by(|left, right| left.text.cmp(&right.text));
        site.y.sort_by(|left, right| left.text.cmp(&right.text));
    }
    labels
}

fn should_render_ladder_index(cleavage_index: usize, residue_count: usize) -> bool {
    let step = if residue_count <= 15 {
        1
    } else if residue_count <= 25 {
        2
    } else {
        5
    };
    cleavage_index == residue_count.saturating_sub(1) || cleavage_index % step == 0
}

fn join_site_label_text(labels: &[SiteLabel]) -> String {
    labels
        .iter()
        .map(|label| label.text.as_str())
        .collect::<Vec<_>>()
        .join(" ")
}

fn join_site_label_titles(labels: &[SiteLabel]) -> String {
    labels
        .iter()
        .map(|label| label.title.as_str())
        .collect::<Vec<_>>()
        .join(" | ")
}

fn ion_table_content_height(row_count: usize) -> f64 {
    let row_height = ion_table_row_height(row_count);
    92.0 + row_count as f64 * row_height + 34.0
}

fn ion_table_row_height(row_count: usize) -> f64 {
    if row_count <= 18 {
        32.0
    } else if row_count <= 28 {
        27.0
    } else {
        23.0
    }
}

fn ion_table_cell_summary_text(entries: &[IonTableEntry]) -> String {
    entries
        .iter()
        .map(|entry| entry.text.as_str())
        .collect::<Vec<_>>()
        .join(" ")
}

fn estimate_text_width(text: &str, font_size: f64) -> f64 {
    let visible_chars = text.chars().count() as f64;
    visible_chars * font_size * 0.62
}

fn estimate_ion_table_width(rows: &[IonTableRow], base_width: f64) -> f64 {
    let layout = estimate_ion_table_layout(rows);
    let max_width = ION_TABLE_MIN_WIDTH + ION_TABLE_MAX_EXTRA_WIDTH.min(base_width * 0.35);
    layout.width.clamp(ION_TABLE_MIN_WIDTH, max_width)
}

fn estimate_ion_table_layout(rows: &[IonTableRow]) -> IonTableLayout {
    let font_size = ion_table_font_size(rows.len());
    let pad = 20.0;
    let cut_width = 58.0;
    let column_gap = 14.0;
    let show_charge2 = ion_table_has_charge2(rows);
    let value_columns = if show_charge2 { 4 } else { 2 };
    let header_labels: &[&str] = if show_charge2 {
        &["b++", "b+", "y+", "y++"]
    } else {
        &["b+", "y+"]
    };
    let mut max_value_width = header_labels
        .iter()
        .map(|label| estimate_text_width(label, ION_TABLE_HEADER_FONT))
        .fold(0.0_f64, f64::max);

    for row in rows {
        let cells: Vec<&IonTableCell> = if show_charge2 {
            vec![&row.b1, &row.b2, &row.y2, &row.y1]
        } else {
            vec![&row.b1, &row.y1]
        };
        for cell in cells {
            max_value_width = max_value_width.max(estimate_text_width(
                &ion_table_cell_summary_text(&cell.entries),
                font_size,
            ));
        }
    }

    let width = pad * 2.0
        + cut_width
        + max_value_width * value_columns as f64
        + column_gap * value_columns as f64;
    IonTableLayout {
        width,
        cut_width,
        pad,
        value_columns,
    }
}

fn ion_table_has_charge2(rows: &[IonTableRow]) -> bool {
    rows.iter()
        .any(|row| !row.b2.entries.is_empty() || !row.y2.entries.is_empty())
}

fn ion_table_font_size(row_count: usize) -> f64 {
    if row_count <= 18 {
        15.0
    } else if row_count <= 28 {
        13.0
    } else {
        11.5
    }
}

fn build_ion_table_rows(report: &AnnotationReport) -> Vec<IonTableRow> {
    let mut match_by_key = HashMap::<IonKey, &FragmentMatch>::new();
    for matched in &report.matches {
        match_by_key.insert(fragment_key(&matched.fragment), matched);
    }

    let residue_count = report.context.peptide.len();
    let mut rows = Vec::with_capacity(residue_count.saturating_sub(1));
    for cleavage_index in 1..residue_count {
        rows.push(IonTableRow {
            cleavage_index,
            y_ordinal: residue_count - cleavage_index,
            b1: build_ion_table_cell(report, &match_by_key, FragmentSeries::B, cleavage_index, 1),
            b2: build_ion_table_cell(report, &match_by_key, FragmentSeries::B, cleavage_index, 2),
            y2: build_ion_table_cell(report, &match_by_key, FragmentSeries::Y, cleavage_index, 2),
            y1: build_ion_table_cell(report, &match_by_key, FragmentSeries::Y, cleavage_index, 1),
        });
    }
    rows
}

fn build_ion_table_cell(
    report: &AnnotationReport,
    match_by_key: &HashMap<IonKey, &FragmentMatch>,
    series: FragmentSeries,
    cleavage_index: usize,
    charge: u8,
) -> IonTableCell {
    let mut fragments = report
        .fragments
        .iter()
        .filter(|fragment| {
            fragment.series == series
                && fragment.cleavage_index == cleavage_index
                && fragment.charge == charge
        })
        .collect::<Vec<_>>();
    fragments.sort_by(|left, right| {
        neutral_loss_rank(left.neutral_loss)
            .cmp(&neutral_loss_rank(right.neutral_loss))
            .then_with(|| {
                left.theoretical_mz
                    .partial_cmp(&right.theoretical_mz)
                    .unwrap_or(Ordering::Equal)
            })
    });

    let entries = fragments
        .into_iter()
        .map(|fragment| {
            let matched = match_by_key.get(&fragment_key(fragment)).copied();
            IonTableEntry {
                text: ion_table_entry_text(fragment),
                color: if matched.is_some() {
                    series_color(series)
                } else {
                    COLOR_SERIES_NEUTRAL
                },
                title: ion_table_entry_title(fragment, matched),
                matched: matched.is_some(),
            }
        })
        .collect();
    IonTableCell { entries }
}

fn fragment_key(fragment: &FragmentIon) -> IonKey {
    IonKey {
        series: fragment.series,
        cleavage_index: fragment.cleavage_index,
        charge: fragment.charge,
        neutral_loss: fragment.neutral_loss,
    }
}

fn neutral_loss_rank(loss: Option<NeutralLossKind>) -> u8 {
    match loss {
        None => 0,
        Some(NeutralLossKind::PhosphoricAcid) => 1,
        Some(NeutralLossKind::Water) => 2,
        Some(NeutralLossKind::Ammonia) => 3,
    }
}

fn ion_table_entry_text(fragment: &FragmentIon) -> String {
    match fragment.neutral_loss {
        None => format!("{:.2}", fragment.theoretical_mz),
        Some(loss) => format!(
            "{}{:.2}",
            neutral_loss_short_label(loss),
            fragment.theoretical_mz
        ),
    }
}

fn neutral_loss_short_label(loss: NeutralLossKind) -> &'static str {
    match loss {
        NeutralLossKind::Water => "w",
        NeutralLossKind::Ammonia => "n",
        NeutralLossKind::PhosphoricAcid => "p",
    }
}

fn ion_table_entry_title(fragment: &FragmentIon, matched: Option<&FragmentMatch>) -> String {
    if let Some(matched) = matched {
        format!(
            "{} theoretical {:.4} | observed {:.4} | error {:+.4} Da / {:+.1} ppm",
            fragment.label(),
            fragment.theoretical_mz,
            matched.observed_mz,
            matched.error_da,
            matched.error_ppm,
        )
    } else {
        format!(
            "{} theoretical {:.4} | no matched peak within tolerance",
            fragment.label(),
            fragment.theoretical_mz,
        )
    }
}

fn draw_ion_table(
    file: &mut File,
    left: f64,
    top: f64,
    width: f64,
    height: f64,
    report: &AnnotationReport,
    rows: &[IonTableRow],
) -> anyhow::Result<()> {
    let row_height = ion_table_row_height(rows.len());
    let font_size = ion_table_font_size(rows.len());
    let column_gap = 14.0;
    let layout = estimate_ion_table_layout(rows);
    let pad = layout.pad;
    let cut_width = layout.cut_width;
    let show_charge2 = layout.value_columns == 4;
    let value_width = (width - pad * 2.0 - cut_width - column_gap * layout.value_columns as f64)
        / layout.value_columns as f64;
    let table_bottom = top + height;
    let title_y = top + 24.0;
    let legend_y = top + 47.0;
    let header_y = top + 76.0;
    let row_start_y = top + 104.0;
    let mut cursor = left + pad;
    let b2_right = if show_charge2 {
        let right = cursor + value_width;
        cursor = right + column_gap;
        Some(right)
    } else {
        None
    };
    let b1_right = cursor + value_width;
    cursor = b1_right + column_gap;
    let cut_left = cursor;
    let cut_center = cut_left + cut_width / 2.0;
    cursor = cut_left + cut_width + column_gap;
    let y1_left = cursor;
    cursor = y1_left + value_width + column_gap;
    let y2_left = show_charge2.then_some(cursor);

    writeln!(
        file,
        r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" rx="8" fill="#f8fafc" stroke="#dbe2ea" stroke-width="1"/>"##,
        x = left,
        y = top,
        w = width,
        h = height,
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}">Full ion table</text>"##,
        x = left + pad,
        y = title_y,
        font = ION_TABLE_TITLE_FONT,
        fill = COLOR_SUBTLE,
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}">matched colored, missing grey; p/w/n = H3PO4/H2O/NH3</text>"##,
        x = left + pad,
        y = legend_y,
        font = ION_TABLE_META_FONT,
        fill = COLOR_SUBTLE,
    )?;

    let mut headers = Vec::new();
    if let Some(x) = b2_right {
        headers.push(("b++", x, "end"));
    }
    headers.push(("b+", b1_right, "end"));
    headers.push(("cut", cut_center, "middle"));
    headers.push(("y+", y1_left, "start"));
    if let Some(x) = y2_left {
        headers.push(("y++", x, "start"));
    }
    for (label, x, anchor) in headers {
        writeln!(
            file,
            r##"<text x="{x:.2}" y="{y:.2}" font-family="{font_family}" font-size="{font_size:.1}" fill="{fill}" text-anchor="{anchor}">{label}</text>"##,
            x = x,
            y = header_y,
            font_family = FONT_TABLE,
            font_size = ION_TABLE_HEADER_FONT,
            fill = COLOR_AXIS,
            anchor = anchor,
            label = label,
        )?;
    }
    writeln!(
        file,
        r##"<line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="#dbe2ea" stroke-width="1"/>"##,
        x1 = left + pad,
        x2 = left + width - pad,
        y = header_y + 6.0,
    )?;

    for (idx, row) in rows.iter().enumerate() {
        let y = row_start_y + idx as f64 * row_height;
        if y + row_height > table_bottom - 8.0 {
            break;
        }
        if idx > 0 {
            writeln!(
                file,
                r##"<line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="#eef2f6" stroke-width="1"/>"##,
                x1 = left + pad,
                x2 = left + width - pad,
                y = y - row_height / 2.0 + 2.0,
            )?;
        }
        if let Some(x) = b2_right {
            write_ion_table_cell(file, x, y, "end", font_size, value_width, &row.b2.entries)?;
        }
        write_ion_table_cell(
            file,
            b1_right,
            y,
            "end",
            font_size,
            value_width,
            &row.b1.entries,
        )?;
        writeln!(
            file,
            r##"<text x="{x:.2}" y="{y:.2}" font-family="{font}" font-size="{size:.1}" fill="{fill}" text-anchor="middle"><title>{title}</title>{text}</text>"##,
            x = cut_center,
            y = y,
            font = FONT_TABLE,
            size = font_size,
            fill = COLOR_SUBTLE,
            title = escape_xml(&format!("b{} / y{}", row.cleavage_index, row.y_ordinal)),
            text = format!("{}|{}", row.cleavage_index, row.y_ordinal),
        )?;
        write_ion_table_cell(
            file,
            y1_left,
            y,
            "start",
            font_size,
            value_width,
            &row.y1.entries,
        )?;
        if let Some(x) = y2_left {
            write_ion_table_cell(file, x, y, "start", font_size, value_width, &row.y2.entries)?;
        }
    }

    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="{font}" fill="{fill}">sequence: {sequence}</text>"##,
        x = left + pad,
        y = table_bottom - 14.0,
        font = ION_TABLE_META_FONT,
        fill = COLOR_SUBTLE,
        sequence = escape_xml(report.context.peptide.sequence()),
    )?;
    Ok(())
}

fn write_ion_table_cell(
    file: &mut File,
    x: f64,
    y: f64,
    anchor: &str,
    font_size: f64,
    max_width: f64,
    entries: &[IonTableEntry],
) -> anyhow::Result<()> {
    if entries.is_empty() {
        return Ok(());
    }

    let title = ion_table_cell_title(entries);
    let visible = ion_table_cell_visible_text(entries, font_size, max_width);
    let fallback_fill = visible
        .first()
        .map(|(entry, _)| entry.color)
        .unwrap_or(COLOR_SERIES_NEUTRAL);

    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="{font}" font-size="{size:.1}" fill="{fill}" text-anchor="{anchor}"><title>{title}</title>"##,
        x = x,
        y = y,
        font = FONT_TABLE,
        size = font_size,
        fill = fallback_fill,
        anchor = anchor,
        title = escape_xml(&title),
    )?;
    for (entry, text) in visible {
        writeln!(
            file,
            r##"<tspan fill="{fill}">{text}</tspan>"##,
            fill = entry.color,
            text = escape_xml(&text),
        )?;
    }
    writeln!(file, "</text>")?;
    Ok(())
}

fn ion_table_cell_title(entries: &[IonTableEntry]) -> String {
    entries
        .iter()
        .map(|entry| entry.title.as_str())
        .collect::<Vec<_>>()
        .join("\n")
}

fn ion_table_cell_visible_text<'a>(
    entries: &'a [IonTableEntry],
    font_size: f64,
    max_width: f64,
) -> Vec<(&'a IonTableEntry, String)> {
    let budget = estimated_text_budget(font_size, max_width);
    if budget == 0 {
        return Vec::new();
    }

    let inline_entries = entries
        .iter()
        .take_while(|entry| entry.matched)
        .count()
        .max(1);
    let mut remaining = budget;
    let mut visible = Vec::new();

    for (idx, entry) in entries.iter().enumerate().take(inline_entries) {
        let prefix = if idx == 0 { "" } else { " " };
        let available = remaining.saturating_sub(prefix.len());
        if available == 0 {
            break;
        }
        let Some(text) = truncate_ion_table_text(&entry.text, available) else {
            break;
        };
        if idx > 0 && text == "…" {
            break;
        }
        remaining = remaining.saturating_sub(prefix.len() + text.chars().count());
        visible.push((entry, format!("{prefix}{text}")));
        if text.ends_with('…') {
            return visible;
        }
    }

    if visible.len() < entries.len() {
        if let Some((_, last_text)) = visible.last_mut() {
            if !last_text.ends_with('…') {
                *last_text = append_ellipsis(last_text);
            }
        } else if let Some(entry) = entries.first() {
            visible.push((entry, "…".to_string()));
        }
    }

    visible
}

fn estimated_text_budget(font_size: f64, max_width: f64) -> usize {
    if !font_size.is_finite() || font_size <= 0.0 || !max_width.is_finite() || max_width <= 0.0 {
        return 0;
    }
    ((max_width / (font_size * 0.62)).floor() as usize).max(1)
}

fn truncate_ion_table_text(text: &str, available: usize) -> Option<String> {
    let len = text.chars().count();
    if len <= available {
        return Some(text.to_string());
    }
    if available == 0 {
        return None;
    }
    if available == 1 {
        return Some("…".to_string());
    }

    let keep = available - 1;
    let truncated = text.chars().take(keep).collect::<String>();
    Some(format!("{truncated}…"))
}

fn append_ellipsis(text: &str) -> String {
    let trimmed = text.trim_end();
    if trimmed.is_empty() || trimmed.ends_with('…') {
        trimmed.to_string()
    } else {
        format!("{trimmed}…")
    }
}

fn collect_peak_labels(
    report: &AnnotationReport,
    spectrum: &LoadedSpectrum,
    normalize: bool,
    x_bounds: CoordinateRange<f64>,
    neutral_loss_label_min_frac: f64,
) -> Vec<PeakLabel> {
    let base_peak = spectrum.stats.base_peak_intensity.max(1.0e-6) as f64;
    let mut grouped = BTreeMap::<usize, Vec<&FragmentMatch>>::new();
    for matched in &report.matches {
        if matched.observed_mz >= x_bounds.min()
            && matched.observed_mz <= x_bounds.max()
            && should_render_match_label(
                matched,
                spectrum.stats.base_peak_intensity,
                neutral_loss_label_min_frac,
            )
        {
            grouped.entry(matched.peak_index).or_default().push(matched);
        }
    }

    let mut out = Vec::new();
    for matches in grouped.into_values() {
        let mut sorted = matches;
        sorted.sort_by(|left, right| {
            left.error_da
                .abs()
                .partial_cmp(&right.error_da.abs())
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    right
                        .observed_intensity
                        .partial_cmp(&left.observed_intensity)
                        .unwrap_or(Ordering::Equal)
                })
        });

        let mut labels = Vec::new();
        let mut seen = HashSet::<String>::new();
        let representative = sorted
            .first()
            .copied()
            .expect("grouped peak labels should not be empty");
        for matched in sorted.into_iter() {
            let label_text = matched.fragment.label();
            if !seen.insert(label_text.clone()) {
                continue;
            }
            labels.push(PeakLabelText {
                text: label_text.clone(),
                color: series_color(matched.fragment.series),
                title: format_match_title(matched),
            });
            if labels.len() >= 2 {
                break;
            }
        }

        if labels.is_empty() {
            continue;
        }

        let display_intensity = if normalize {
            representative.observed_intensity as f64 / base_peak
        } else {
            representative.observed_intensity as f64
        };
        out.push(PeakLabel {
            observed_mz: representative.observed_mz,
            display_intensity,
            labels,
        });
    }

    out.sort_by(|left, right| {
        left.observed_mz
            .partial_cmp(&right.observed_mz)
            .unwrap_or(Ordering::Equal)
    });
    out
}

fn should_render_match_label(
    matched: &FragmentMatch,
    base_peak_intensity: f32,
    neutral_loss_label_min_frac: f64,
) -> bool {
    if matched.fragment.neutral_loss.is_none() {
        return true;
    }
    if base_peak_intensity <= 0.0 {
        return true;
    }
    (matched.observed_intensity as f64) / (base_peak_intensity as f64)
        >= neutral_loss_label_min_frac
}

fn format_match_title(matched: &FragmentMatch) -> String {
    format!(
        "{} | observed {:.4} | theoretical {:.4} | error {:+.4} Da / {:+.1} ppm",
        matched.fragment.label(),
        matched.observed_mz,
        matched.fragment.theoretical_mz,
        matched.error_da,
        matched.error_ppm,
    )
}

fn series_color(series: FragmentSeries) -> &'static str {
    match series {
        FragmentSeries::B => COLOR_SERIES_B,
        FragmentSeries::Y => COLOR_SERIES_Y,
    }
}

fn sanitize_filename_component(input: &str) -> String {
    let mut out = String::with_capacity(input.len().min(80));
    for ch in input.chars().take(80) {
        match ch {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' | '.' => out.push(ch),
            _ => out.push('_'),
        }
    }
    if out.is_empty() {
        "scan".to_string()
    } else {
        out
    }
}

fn sanitize_filename_label(input: &str) -> Option<String> {
    let sanitized = sanitize_filename_component(input);
    let trimmed = sanitized.trim_matches(|ch| ch == '_' || ch == '-' || ch == '.');
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_string())
    }
}

fn default_peptide_label(sequence: &str, charge: Option<i32>) -> String {
    let mut label = sanitize_filename_component(sequence);
    if let Some(charge) = charge.filter(|value| *value > 0) {
        label.push_str(&charge.to_string());
    }
    label
}

fn escape_xml(input: &str) -> String {
    input
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('\"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mzml::{native_id_matches_query, SpectrumMeta, SpectrumStats};

    #[test]
    fn parse_plot_args_rejects_conflicting_tolerances() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--index".to_string(),
            "1".to_string(),
            "--peptide".to_string(),
            "PEPTIDE".to_string(),
            "--tol-ppm".to_string(),
            "20".to_string(),
            "--tol-da".to_string(),
            "0.5".to_string(),
        ])
        .expect_err("conflicting tolerances should fail");
        assert!(err
            .to_string()
            .contains("only one of --tol-ppm or --tol-da"));
    }

    #[test]
    fn parse_plot_args_accepts_peptide_mods_and_tol_da() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--id".to_string(),
            "scan=2".to_string(),
            "--peptide".to_string(),
            "DSAVYFCARTKILDFD".to_string(),
            "--mod".to_string(),
            "7:+57.021464".to_string(),
            "--tol-da".to_string(),
            "0.5".to_string(),
        ])
        .expect("options parse");
        assert!(matches!(
            options.selector,
            Some(SpectrumSelector::NativeId(_))
        ));
        assert_eq!(options.mod_inputs, vec!["7:+57.021464"]);
        assert_eq!(options.tolerance, MassTolerance::Da(0.5));
    }

    #[test]
    fn parse_plot_args_accepts_pepxml_top_n_and_annotation_options() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "2".to_string(),
            "--pepxml".to_string(),
            "search.pep.xml".to_string(),
            "--top-n".to_string(),
            "3".to_string(),
            "--tol-da".to_string(),
            "0.5".to_string(),
            "--neutral-losses".to_string(),
        ])
        .expect("options parse");
        assert_eq!(
            options.pepxml_path.as_deref(),
            Some(Path::new("search.pep.xml"))
        );
        assert_eq!(options.top_n, 3);
        assert_eq!(options.tolerance, MassTolerance::Da(0.5));
        assert!(options.neutral_losses_enabled);
    }

    #[test]
    fn parse_plot_args_rejects_peptide_and_pepxml_together() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "2".to_string(),
            "--peptide".to_string(),
            "PEPTIDEK".to_string(),
            "--pepxml".to_string(),
            "search.pep.xml".to_string(),
        ])
        .expect_err("conflicting annotation inputs should fail");
        assert!(err
            .to_string()
            .contains("only one of --peptide or --pepxml"));
    }

    #[test]
    fn parse_plot_args_rejects_mod_with_pepxml() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "2".to_string(),
            "--pepxml".to_string(),
            "search.pep.xml".to_string(),
            "--mod".to_string(),
            "4:+15.9949".to_string(),
        ])
        .expect_err("manual mods with pepXML should fail");
        assert!(err
            .to_string()
            .contains("--mod cannot be combined with --pepxml"));
    }

    #[test]
    fn parse_plot_args_rejects_top_n_without_pepxml() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "2".to_string(),
            "--top-n".to_string(),
            "3".to_string(),
            "--peptide".to_string(),
            "PEPTIDEK".to_string(),
        ])
        .expect_err("top-n without pepXML should fail");
        assert!(err.to_string().contains("--top-n requires --pepxml"));
    }

    #[test]
    fn parse_plot_args_rejects_neutral_losses_without_peptide() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--id".to_string(),
            "scan=2".to_string(),
            "--neutral-losses".to_string(),
        ])
        .expect_err("neutral losses without peptide should fail");
        assert!(err
            .to_string()
            .contains("--neutral-losses requires --peptide"));
    }

    #[test]
    fn parse_plot_args_accepts_neutral_losses_with_peptide() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--id".to_string(),
            "scan=2".to_string(),
            "--peptide".to_string(),
            "DSAVYFCARTKILDFD".to_string(),
            "--neutral-losses".to_string(),
        ])
        .expect("options parse");
        assert!(options.neutral_losses_enabled);
    }

    #[test]
    fn parse_plot_args_accepts_remove_precursor() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--id".to_string(),
            "scan=2".to_string(),
            "--remove-precursor".to_string(),
        ])
        .expect("options parse");
        assert!(options.remove_precursor);
    }

    #[test]
    fn parse_plot_args_accepts_neutral_loss_label_threshold() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "107468".to_string(),
            "--peptide".to_string(),
            "DSAVYFCARTKILDF".to_string(),
            "--neutral-losses".to_string(),
            "--neutral-loss-min-frac".to_string(),
            "0.05".to_string(),
        ])
        .expect("options parse");
        assert!((options.neutral_loss_label_min_frac - 0.05).abs() < 1e-12);
    }

    #[test]
    fn parse_plot_args_accepts_isotope_errors() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "107468".to_string(),
            "--peptide".to_string(),
            "[+304.2071]THS[+79.9663]GS[+79.9663]SGGSGSR/3".to_string(),
            "--isotope-errors".to_string(),
            "0,1,2".to_string(),
        ])
        .expect("options parse");
        assert_eq!(options.isotope_errors, vec![0, 1, 2]);
    }

    #[test]
    fn parse_plot_args_rejects_neutral_loss_threshold_without_neutral_losses() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "107468".to_string(),
            "--peptide".to_string(),
            "DSAVYFCARTKILDF".to_string(),
            "--neutral-loss-min-frac".to_string(),
            "0.05".to_string(),
        ])
        .expect_err("threshold without neutral-losses should fail");
        assert!(err
            .to_string()
            .contains("--neutral-loss-min-frac requires --neutral-losses"));
    }

    #[test]
    fn parse_plot_args_rejects_isotope_errors_without_peptide() {
        let err = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "107468".to_string(),
            "--isotope-errors".to_string(),
            "0,1,2".to_string(),
        ])
        .expect_err("isotope errors without peptide should fail");
        assert!(err
            .to_string()
            .contains("--isotope-errors requires --peptide"));
    }

    #[test]
    fn parse_plot_args_accepts_scan_selector_and_svg_prefix() {
        let options = parse_plot_args(vec![
            "--mzml".to_string(),
            "x.mzML".to_string(),
            "--scan".to_string(),
            "107468".to_string(),
            "--svg-prefix".to_string(),
            "calibrated_".to_string(),
        ])
        .expect("options parse");
        assert!(matches!(
            options.selector,
            Some(SpectrumSelector::ScanNumber(107468))
        ));
        assert_eq!(options.svg_prefix.as_deref(), Some("calibrated_"));
    }

    #[test]
    fn build_ruler_ticks_prefers_dense_minor_ticks_and_major_labels() {
        let ticks = build_ruler_ticks(141.23, 1748.86, 1200.0);
        assert!(ticks.len() > 20);
        assert!(ticks
            .windows(2)
            .all(|window| window[1].value > window[0].value));
        assert!(ticks
            .iter()
            .any(|tick| matches!(tick.kind, RulerTickKind::Major)));
        assert!(ticks
            .iter()
            .any(|tick| matches!(tick.kind, RulerTickKind::Medium)));
        assert!(ticks
            .iter()
            .any(|tick| matches!(tick.kind, RulerTickKind::Minor)));
    }

    #[test]
    fn downsample_max_per_bin_can_exclude_precursor_window() {
        let mz = [95.0, 100.0, 105.0];
        let intensity = [10.0_f32, 100.0_f32, 20.0_f32];
        let (points, y_max) = downsample_max_per_bin(
            &mz,
            &intensity,
            CoordinateRange::new(90.0, 110.0),
            false,
            Some((99.5, 100.5)),
            64,
        );
        assert_eq!(points.len(), 2);
        assert!(points.iter().all(|(mz, _)| (*mz - 100.0).abs() > 1e-9));
        assert!((y_max - 22.0).abs() < 1e-9);
    }

    #[test]
    fn ion_table_rows_include_missing_and_neutral_loss_entries() {
        let context = prepare_annotation(
            "[+304.2071]SLES[+79.9663]DNEEK[+304.2071]/3",
            &[],
            &[
                NeutralLossKind::Water,
                NeutralLossKind::Ammonia,
                NeutralLossKind::PhosphoricAcid,
            ],
            Some(3),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        let rows = build_ion_table_rows(&report);
        assert_eq!(rows.len(), 8);
        assert!(rows[3]
            .b1
            .entries
            .iter()
            .any(|entry| entry.text.starts_with('p')));
        assert!(rows[0]
            .y1
            .entries
            .iter()
            .any(|entry| entry.color == COLOR_SERIES_NEUTRAL));
        assert!(ion_table_has_charge2(&rows));
        assert_eq!(estimate_ion_table_layout(&rows).value_columns, 4);
    }

    #[test]
    fn ion_table_hides_charge_two_columns_when_no_charge_two_fragments_exist() {
        let context = prepare_annotation("PEPTIDEK/2", &[], &[], Some(2), MassTolerance::Ppm(20.0))
            .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        let rows = build_ion_table_rows(&report);

        assert!(!ion_table_has_charge2(&rows));
        assert_eq!(estimate_ion_table_layout(&rows).value_columns, 2);
    }

    #[test]
    fn ion_table_cell_visible_text_prioritizes_matched_entries_and_preserves_full_title() {
        let entries = vec![
            IonTableEntry {
                text: "123.45".to_string(),
                color: COLOR_SERIES_B,
                title: "b4 theoretical 123.4500 | observed 123.4501".to_string(),
                matched: true,
            },
            IonTableEntry {
                text: "p98.76".to_string(),
                color: COLOR_SERIES_B,
                title: "b4-H3PO4 theoretical 98.7600 | observed 98.7602".to_string(),
                matched: true,
            },
            IonTableEntry {
                text: "w105.32".to_string(),
                color: COLOR_SERIES_NEUTRAL,
                title: "b4-H2O theoretical 105.3200 | no matched peak within tolerance".to_string(),
                matched: false,
            },
        ];

        let visible = ion_table_cell_visible_text(&entries, 9.5, 52.0);
        let rendered = visible
            .iter()
            .map(|(_, text)| text.as_str())
            .collect::<String>();

        assert_eq!(rendered, "123.45…");
        assert_eq!(ion_table_cell_title(&entries).lines().count(), 3);
    }

    #[test]
    fn ladder_indices_scale_down_for_longer_peptides() {
        assert!(should_render_ladder_index(1, 10));
        assert!(should_render_ladder_index(9, 10));
        assert!(!should_render_ladder_index(1, 20));
        assert!(should_render_ladder_index(2, 20));
        assert!(should_render_ladder_index(19, 20));
        assert!(!should_render_ladder_index(3, 30));
        assert!(should_render_ladder_index(5, 30));
        assert!(should_render_ladder_index(29, 30));
    }

    #[test]
    fn default_output_path_includes_scan_peptide_and_neutral_loss_state() {
        let spectrum = LoadedSpectrum {
            meta: SpectrumMeta {
                idx: 107467,
                scan_id: "controllerType=0 controllerNumber=1 scan=107468".to_string(),
                ms_level: 2,
                rt_minutes: None,
                precursor_mz: Some(874.9339),
                precursor_charge: Some(2),
                continuity: SignalContinuity::Centroid,
            },
            mz: Vec::new(),
            intensity: Vec::new(),
            stats: SpectrumStats {
                points: 0,
                mz_min: 0.0,
                mz_max: 0.0,
                base_peak_mz: 0.0,
                base_peak_intensity: 0.0,
            },
        };
        let context = prepare_annotation(
            "DSAVYFCARTKILDF/2",
            &[],
            &[],
            Some(2),
            MassTolerance::Da(0.5),
        )
        .expect("annotation context");
        let path = default_output_path(
            Path::new("/tmp/example_run.mzML"),
            &spectrum,
            Some(&context),
            Some(2),
            true,
            Some("calibrated_"),
        );
        let filename = path
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or("");
        assert!(filename.contains("calibrated"));
        assert!(filename.contains("example_run"));
        assert!(filename.contains("scan107468"));
        assert!(filename.contains("DSAVYFCARTKILDF2"));
        assert!(filename.contains("nl-on"));
    }

    #[test]
    fn should_render_match_label_filters_weak_neutral_losses_only() {
        let regular = FragmentMatch {
            fragment: crate::annotate::FragmentIon {
                series: FragmentSeries::Y,
                ordinal: 9,
                cleavage_index: 6,
                charge: 1,
                neutral_loss: None,
                theoretical_mz: 1078.5891,
            },
            peak_index: 0,
            observed_mz: 1078.5759,
            observed_intensity: 10.0,
            error_da: -0.0132,
            error_ppm: -12.2,
        };
        let neutral_loss = FragmentMatch {
            fragment: crate::annotate::FragmentIon {
                series: FragmentSeries::Y,
                ordinal: 9,
                cleavage_index: 6,
                charge: 1,
                neutral_loss: Some(NeutralLossKind::Water),
                theoretical_mz: 1060.5786,
            },
            peak_index: 1,
            observed_mz: 1060.6544,
            observed_intensity: 2.0,
            error_da: 0.0758,
            error_ppm: 71.5,
        };
        assert!(should_render_match_label(&regular, 100.0, 0.03));
        assert!(!should_render_match_label(&neutral_loss, 100.0, 0.03));
        assert!(should_render_match_label(&neutral_loss, 100.0, 0.02));
    }

    #[test]
    fn native_id_matches_thermo_scan_shorthand() {
        assert!(native_id_matches_query(
            "controllerType=0 controllerNumber=1 scan=107468",
            "scan=107468",
        ));
    }

    #[test]
    fn native_id_rejects_different_scan_numbers() {
        assert!(!native_id_matches_query(
            "controllerType=0 controllerNumber=1 scan=107468",
            "scan=107469",
        ));
    }
}
