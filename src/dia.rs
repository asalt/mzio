use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use libloading::Library;
use mzdata::io::SpectrumSource;
use serde::{Deserialize, Serialize};
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::{AcquisitionType, Frame, MSLevel, Metadata};

use crate::annotate::{
    fragment_charge_states, generate_fragments, prepare_annotation, FragmentIon, FragmentSeries,
    MassTolerance, NeutralLossKind,
};
use crate::mzml::{load_spectrum_by_index, open_reader};
use crate::scale::CoordinateRange;
use crate::svg_canvas::{AxisOrientation, AxisProps, AxisTickLabelStyle, SvgCanvas};

const BRUKER_BRIDGE_PY: &str = include_str!("dia_bruker_bridge.py");
const DEFAULT_MZ_PPM: f64 = 20.0;
const DEFAULT_MZ_PROFILE_BINS: usize = 160;
const DEFAULT_BRUKER_SO_PATHS: &[&str] =
    &["/opt/bruker/linux64/timsdata.so", "/opt/bruker/timsdata.so"];
const BRUKER_SO_ENV_VAR: &str = "MZIO_BRUKER_SO";
const BRUKER_PYTHON_ENV_VAR: &str = "MZIO_PYTHON";
const DEFAULT_PEPTIDE_CHARGE: i32 = 2;
const DEFAULT_PSEUDO_MS2_RT_WINDOW_MIN: f64 = 0.5;
const DEFAULT_RT_SMOOTH_WINDOW_POINTS: usize = 9;
const RT_SMOOTH_ZERO_WEIGHT: f64 = 0.25;
const PSEUDO_MS2_NEUTRAL_LOSS_MIN_RELATIVE_INTENSITY: f64 = 0.03;
const COMMON_NEUTRAL_LOSSES: [NeutralLossKind; 3] = [
    NeutralLossKind::Water,
    NeutralLossKind::Ammonia,
    NeutralLossKind::PhosphoricAcid,
];
const SVG_WIDTH: u32 = 1320;
const SVG_HEIGHT: u32 = 980;
const SVG_MARGIN_X: f64 = 60.0;
const SVG_TITLE_FONT: f64 = 28.0;
const SVG_META_FONT: f64 = 15.0;
const SVG_PANEL_TITLE_FONT: f64 = 22.0;
const SVG_TICK_FONT: f64 = 14.0;
const SVG_AXIS_LABEL_FONT: f64 = 16.0;
const COLOR_TEXT: &str = "#122033";
const COLOR_SUBTLE: &str = "#5b6775";
const COLOR_CARD_BORDER: &str = "#d8e0ea";
const COLOR_AXIS: &str = "#334155";
const COLOR_RT: &str = "#1d4ed8";
const COLOR_MZ: &str = "#b45309";
const COLOR_GRID: &str = "#e5ebf2";
const COLOR_BG: &str = "#ffffff";

#[derive(Clone, Debug)]
enum DiaInput {
    Bruker(PathBuf),
    Mzml(PathBuf),
}

impl DiaInput {
    fn path(&self) -> &Path {
        match self {
            Self::Bruker(path) | Self::Mzml(path) => path.as_path(),
        }
    }

    fn default_stem(&self) -> String {
        let path = self.path();
        path.file_stem()
            .or_else(|| path.file_name())
            .and_then(|value| value.to_str())
            .filter(|value| !value.is_empty())
            .unwrap_or("dia_slice")
            .to_string()
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct DiaCapabilities {
    pub(crate) has_mobility: bool,
    pub(crate) has_isolation_window: bool,
    pub(crate) requires_vendor_runtime: bool,
}

impl DiaCapabilities {
    fn labels(self) -> String {
        let mut labels = vec!["rt", "mz"];
        if self.has_mobility {
            labels.push("mobility");
        }
        if self.has_isolation_window {
            labels.push("quad-window");
        }
        if self.requires_vendor_runtime {
            labels.push("vendor-runtime");
        }
        labels.join(", ")
    }
}

fn compact_backend_label(label: &str, acquisition_mode: Option<&str>) -> String {
    let backend = match label {
        "Bruker .d (timsrust native)" => "Bruker native",
        "Bruker .d (alphaTims bridge)" => "Bruker alphaTims",
        other => other,
    };
    match acquisition_mode.filter(|mode| !mode.is_empty()) {
        Some(mode) if backend.starts_with("Bruker") => format!("{backend} {mode}"),
        _ => backend.to_string(),
    }
}

fn compact_capability_label(capabilities: DiaCapabilities) -> String {
    let mut dimensions = vec!["RT", "m/z"];
    if capabilities.has_mobility {
        dimensions.push("mobility");
    }
    let mut label = dimensions.join("/");
    if capabilities.has_isolation_window {
        label.push_str("; quad windows");
    }
    if capabilities.requires_vendor_runtime {
        label.push_str("; vendor runtime");
    }
    label
}

#[derive(Clone, Debug)]
pub(crate) struct DiaSliceRequest {
    pub(crate) mz: f64,
    pub(crate) mz_ppm: f64,
    pub(crate) mz_da: Option<f64>,
    pub(crate) peptide_target: Option<DiaPeptideTarget>,
    pub(crate) rt_min: Option<f64>,
    pub(crate) rt_max: Option<f64>,
    pub(crate) im_min: Option<f64>,
    pub(crate) im_max: Option<f64>,
    pub(crate) quad_min: Option<f64>,
    pub(crate) quad_max: Option<f64>,
}

impl DiaSliceRequest {
    fn mz_bounds(&self) -> (f64, f64) {
        let delta = self
            .mz_da
            .unwrap_or_else(|| self.mz * self.mz_ppm / 1_000_000.0);
        (self.mz - delta, self.mz + delta)
    }

    fn rt_window_label(&self) -> String {
        match (self.rt_min, self.rt_max) {
            (Some(lo), Some(hi)) => format!("{lo:.3}-{hi:.3} min"),
            _ => "all RT".to_string(),
        }
    }
}

#[derive(Clone, Debug)]
struct DiaMzTarget {
    label: String,
    mz: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct DiaPeptideTarget {
    pub(crate) input: String,
    pub(crate) sequence: String,
    pub(crate) modified_sequence: String,
    pub(crate) charge: i32,
    pub(crate) precursor_mz: f64,
    pub(crate) fragment: Option<DiaFragmentTarget>,
    pub(crate) fragments: Vec<DiaFragmentTarget>,
}

#[derive(Clone, Debug)]
pub(crate) struct DiaFragmentTarget {
    pub(crate) input: String,
    pub(crate) label: String,
    pub(crate) series: String,
    pub(crate) cleavage_index: usize,
    pub(crate) neutral_loss: Option<NeutralLossKind>,
    pub(crate) charge: u8,
    pub(crate) mz: f64,
}

#[derive(Clone, Debug)]
struct DiaSliceOptions {
    input: Option<DiaInput>,
    bruker_so: Option<PathBuf>,
    python_bin: Option<PathBuf>,
    bruker_backend: BrukerBackend,
    request: Option<DiaSliceRequest>,
    out_prefix: Option<PathBuf>,
    outdir: Option<PathBuf>,
    mz_profile_bins: usize,
    mz_targets: Vec<DiaMzTarget>,
    peptide_input: Option<String>,
    fragment_input: Option<String>,
    mod_inputs: Vec<String>,
    charge_override: Option<i32>,
    neutral_losses_enabled: bool,
    pseudo_ms2: bool,
    pseudo_ms2_rt_window_min: f64,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
}

impl Default for DiaSliceOptions {
    fn default() -> Self {
        Self {
            input: None,
            bruker_so: None,
            python_bin: None,
            bruker_backend: BrukerBackend::Auto,
            request: None,
            out_prefix: None,
            outdir: None,
            mz_profile_bins: DEFAULT_MZ_PROFILE_BINS,
            mz_targets: Vec::new(),
            peptide_input: None,
            fragment_input: None,
            mod_inputs: Vec::new(),
            charge_override: None,
            neutral_losses_enabled: false,
            pseudo_ms2: false,
            pseudo_ms2_rt_window_min: DEFAULT_PSEUDO_MS2_RT_WINDOW_MIN,
            rt_smooth: false,
            rt_smooth_window: DEFAULT_RT_SMOOTH_WINDOW_POINTS,
            verbosity: Verbosity::Normal,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Verbosity {
    Quiet,
    Normal,
    Verbose,
}

impl Verbosity {
    fn status(self, args: std::fmt::Arguments<'_>) {
        if self != Self::Quiet {
            eprintln!("{args}");
        }
    }

    fn detail(self, args: std::fmt::Arguments<'_>) {
        if self == Self::Verbose {
            eprintln!("{args}");
        }
    }

    fn success(self, args: std::fmt::Arguments<'_>) {
        if self != Self::Quiet {
            println!("{args}");
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct TimsrustFrameCandidate {
    reader_index: usize,
    rt_seconds: f64,
}

#[derive(Clone, Debug)]
struct TimsrustExtraction {
    payloads: Vec<BrukerBridgePayload>,
    run_tic: BrukerRunTicSummary,
}

#[derive(Clone, Debug, Serialize)]
struct BrukerRunTicSummary {
    schema_version: u8,
    source: PathBuf,
    backend: &'static str,
    acquisition_mode: String,
    intensity_column: &'static str,
    total_frames: usize,
    rt_min_minutes: Option<f64>,
    rt_max_minutes: Option<f64>,
    ms1: BrukerRunTicLevelSummary,
    ms2: BrukerRunTicLevelSummary,
    unknown: BrukerRunTicLevelSummary,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
struct BrukerRunTicLevelSummary {
    frames: usize,
    frames_with_signal: usize,
    detector_events: usize,
    summed_tic: f64,
    max_tic: f64,
    max_tic_frame: Option<usize>,
    max_tic_rt_minutes: Option<f64>,
}

#[derive(Debug)]
struct BrukerRunTicAccumulator {
    total_frames: usize,
    rt_min_seconds: f64,
    rt_max_seconds: f64,
    ms1: BrukerRunTicLevelAccumulator,
    ms2: BrukerRunTicLevelAccumulator,
    unknown: BrukerRunTicLevelAccumulator,
}

#[derive(Debug, Default)]
struct BrukerRunTicLevelAccumulator {
    frames: usize,
    frames_with_signal: usize,
    detector_events: usize,
    summed_tic: f64,
    max_tic: f64,
    max_tic_frame: Option<usize>,
    max_tic_rt_minutes: Option<f64>,
}

impl BrukerRunTicAccumulator {
    fn new() -> Self {
        Self {
            total_frames: 0,
            rt_min_seconds: f64::INFINITY,
            rt_max_seconds: f64::NEG_INFINITY,
            ms1: BrukerRunTicLevelAccumulator::default(),
            ms2: BrukerRunTicLevelAccumulator::default(),
            unknown: BrukerRunTicLevelAccumulator::default(),
        }
    }

    fn add_frame(&mut self, frame: &Frame) {
        self.total_frames = self.total_frames.saturating_add(1);
        if frame.rt_in_seconds.is_finite() {
            self.rt_min_seconds = self.rt_min_seconds.min(frame.rt_in_seconds);
            self.rt_max_seconds = self.rt_max_seconds.max(frame.rt_in_seconds);
        }

        let tic = corrected_frame_tic(frame);
        let rt_minutes = frame
            .rt_in_seconds
            .is_finite()
            .then_some(frame.rt_in_seconds / 60.0);
        match frame.ms_level {
            MSLevel::MS1 => {
                self.ms1
                    .add_frame(frame.index, rt_minutes, frame.intensities.len(), tic)
            }
            MSLevel::MS2 => {
                self.ms2
                    .add_frame(frame.index, rt_minutes, frame.intensities.len(), tic)
            }
            MSLevel::Unknown => {
                self.unknown
                    .add_frame(frame.index, rt_minutes, frame.intensities.len(), tic)
            }
        }
    }

    fn finalize(self, source: &Path, acquisition_mode: String) -> BrukerRunTicSummary {
        BrukerRunTicSummary {
            schema_version: 1,
            source: source.to_path_buf(),
            backend: "Bruker .d (timsrust native)",
            acquisition_mode,
            intensity_column: "timsrust_corrected_intensity_values",
            total_frames: self.total_frames,
            rt_min_minutes: finite_seconds_to_minutes(self.rt_min_seconds),
            rt_max_minutes: finite_seconds_to_minutes(self.rt_max_seconds),
            ms1: self.ms1.finalize(),
            ms2: self.ms2.finalize(),
            unknown: self.unknown.finalize(),
        }
    }
}

impl BrukerRunTicLevelAccumulator {
    fn add_frame(
        &mut self,
        frame_index: usize,
        rt_minutes: Option<f64>,
        detector_events: usize,
        tic: f64,
    ) {
        self.frames = self.frames.saturating_add(1);
        self.detector_events = self.detector_events.saturating_add(detector_events);
        if tic.is_finite() && tic > 0.0 {
            self.frames_with_signal = self.frames_with_signal.saturating_add(1);
            self.summed_tic += tic;
            if tic > self.max_tic {
                self.max_tic = tic;
                self.max_tic_frame = Some(frame_index);
                self.max_tic_rt_minutes = rt_minutes;
            }
        }
    }

    fn finalize(self) -> BrukerRunTicLevelSummary {
        BrukerRunTicLevelSummary {
            frames: self.frames,
            frames_with_signal: self.frames_with_signal,
            detector_events: self.detector_events,
            summed_tic: self.summed_tic,
            max_tic: self.max_tic,
            max_tic_frame: self.max_tic_frame,
            max_tic_rt_minutes: self.max_tic_rt_minutes,
        }
    }
}

fn corrected_frame_tic(frame: &Frame) -> f64 {
    let factor = frame.intensity_correction_factor;
    if !factor.is_finite() || factor <= 0.0 {
        return 0.0;
    }
    frame
        .intensities
        .iter()
        .map(|value| *value as f64 * factor)
        .filter(|value| value.is_finite() && *value > 0.0)
        .sum()
}

fn finite_seconds_to_minutes(value: f64) -> Option<f64> {
    value.is_finite().then_some(value / 60.0)
}

#[derive(Debug)]
struct ProgressReporter {
    enabled: bool,
    interactive: bool,
    label: String,
    total: usize,
    current: usize,
    last_percent: Option<usize>,
    next_log_percent: usize,
    finished: bool,
}

impl ProgressReporter {
    fn new(label: impl Into<String>, total: usize, verbosity: Verbosity) -> Self {
        Self {
            enabled: verbosity != Verbosity::Quiet && total > 0,
            interactive: io::stderr().is_terminal(),
            label: label.into(),
            total,
            current: 0,
            last_percent: None,
            next_log_percent: 10,
            finished: false,
        }
    }

    fn advance(&mut self) {
        if !self.enabled {
            return;
        }
        self.current = (self.current + 1).min(self.total);
        let percent = progress_percent(self.current, self.total);
        if self.interactive {
            if self.last_percent != Some(percent) || self.current == self.total {
                self.render_inline(percent);
                self.last_percent = Some(percent);
            }
        } else if percent >= self.next_log_percent || self.current == self.total {
            eprintln!(
                "{}: {}/{} frames ({}%)",
                self.label, self.current, self.total, percent
            );
            while self.next_log_percent <= percent {
                self.next_log_percent += 10;
            }
            self.last_percent = Some(percent);
        }
    }

    fn finish(&mut self) {
        if !self.enabled || self.finished {
            return;
        }
        if self.current < self.total {
            self.current = self.total;
            let percent = progress_percent(self.current, self.total);
            if self.interactive {
                self.render_inline(percent);
            } else if self.last_percent != Some(percent) {
                eprintln!(
                    "{}: {}/{} frames ({}%)",
                    self.label, self.current, self.total, percent
                );
            }
        }
        if self.interactive {
            eprintln!();
        }
        self.finished = true;
    }

    fn render_inline(&self, percent: usize) {
        let bar = progress_bar_body(self.current, self.total, 28);
        eprint!(
            "\r{}: [{}] {}/{} frames ({:>3}%)",
            self.label, bar, self.current, self.total, percent
        );
        let _ = io::stderr().flush();
    }
}

impl Drop for ProgressReporter {
    fn drop(&mut self) {
        if self.enabled && self.interactive && !self.finished {
            eprintln!();
        }
    }
}

fn progress_percent(current: usize, total: usize) -> usize {
    if total == 0 {
        100
    } else {
        ((current.min(total) * 100) / total).min(100)
    }
}

fn progress_bar_body(current: usize, total: usize, width: usize) -> String {
    if width == 0 {
        return String::new();
    }
    let filled = if total == 0 {
        width
    } else {
        (current.min(total) * width) / total
    };
    let mut body = String::with_capacity(width);
    body.push_str(&"=".repeat(filled));
    body.push_str(&".".repeat(width.saturating_sub(filled)));
    body
}

#[derive(Clone, Debug)]
struct PseudoMs2Options {
    enabled: bool,
    rt_window_min: f64,
}

#[derive(Clone, Debug)]
struct PseudoMs2Report {
    input_path: PathBuf,
    out_prefix: PathBuf,
    peptide: DiaPeptideTarget,
    rt_window: PseudoMs2RtWindow,
    frames_considered: usize,
    frames_with_signal: usize,
    matched_events: usize,
    precursor_frames_with_signal: usize,
    precursor_apex_rt: Option<f64>,
    precursor_apex_intensity: f64,
    fragments: Vec<PseudoMs2FragmentEvidence>,
}

#[derive(Clone, Debug)]
struct PseudoMs2RtWindow {
    min: f64,
    max: f64,
    source: &'static str,
}

impl PseudoMs2RtWindow {
    fn label(&self) -> String {
        format!("{:.3}-{:.3} min ({})", self.min, self.max, self.source)
    }
}

#[derive(Clone, Debug)]
struct PseudoMs2FragmentEvidence {
    label: String,
    series: String,
    cleavage_index: usize,
    neutral_loss: Option<NeutralLossKind>,
    charge: u8,
    mz: f64,
    summed_intensity: f64,
    matched_events: usize,
    frames_with_signal: usize,
    apex_rt: Option<f64>,
    apex_intensity: f64,
}

#[derive(Clone, Debug)]
struct PseudoIonTableRow {
    cleavage_index: usize,
    y_ordinal: usize,
    b1: Vec<PseudoIonTableCell>,
    b2: Vec<PseudoIonTableCell>,
    y1: Vec<PseudoIonTableCell>,
    y2: Vec<PseudoIonTableCell>,
}

#[derive(Clone, Debug)]
struct PseudoIonTableCell {
    label: String,
    series: String,
    neutral_loss: Option<NeutralLossKind>,
    mz: f64,
    summed_intensity: f64,
    matched_events: usize,
    frames_with_signal: usize,
    apex_rt: Option<f64>,
}

impl PseudoIonTableCell {
    fn detected(&self) -> bool {
        self.summed_intensity > 0.0
    }
}

#[derive(Clone, Debug)]
struct RtProfileRow {
    scan_index: u32,
    scan_id: String,
    rt_minutes: Option<f64>,
    summed_intensity: f64,
    matched_peaks: usize,
    precursor_mz: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
struct MzProfileBin {
    mz_center: f64,
    summed_intensity: f64,
}

#[derive(Clone, Debug)]
struct DiaSliceSummary {
    backend_label: &'static str,
    input_path: PathBuf,
    out_prefix: PathBuf,
    mz_min: f64,
    mz_max: f64,
    spectra_considered: usize,
    spectra_with_signal: usize,
    matched_peaks: usize,
    capabilities: DiaCapabilities,
    acquisition_mode: Option<String>,
    intensity_column: Option<String>,
    vendor_runtime_path: Option<PathBuf>,
}

#[derive(Clone, Debug)]
struct BrukerRuntime {
    path: PathBuf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BrukerBackend {
    Auto,
    Native,
    AlphaTims,
}

impl BrukerBackend {
    fn parse(value: &str) -> anyhow::Result<Self> {
        match value {
            "auto" => Ok(Self::Auto),
            "native" | "timsrust" => Ok(Self::Native),
            "alphatims" | "python" => Ok(Self::AlphaTims),
            other => anyhow::bail!(
                "unknown --bruker-backend `{other}`; expected auto, native, or alphatims"
            ),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct MobilityProfileRow {
    mobility: f64,
    summed_intensity: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct BrukerBridgePayload {
    acquisition_mode: String,
    intensity_column: String,
    mz_min: f64,
    mz_max: f64,
    frames_considered: usize,
    frames_with_signal: usize,
    matched_events: usize,
    rt_profile: Vec<BrukerBridgeRtRow>,
    mz_profile: Vec<BrukerBridgeMzRow>,
    im_profile: Vec<BrukerBridgeImRow>,
}

#[derive(Clone, Debug, Deserialize)]
struct BrukerBridgeRtRow {
    frame_index: u32,
    rt_minutes: f64,
    summed_intensity: f64,
    matched_events: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct BrukerBridgeMzRow {
    mz_center: f64,
    summed_intensity: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct BrukerBridgeImRow {
    mobility: f64,
    summed_intensity: f64,
}

struct TimsrustSliceAccumulator {
    mz_min: f64,
    mz_max: f64,
    tof_min: u32,
    tof_max: u32,
    mz_bins: Vec<MzProfileBin>,
    frames_with_signal: usize,
    matched_events: usize,
    rt_profile: Vec<BrukerBridgeRtRow>,
    im_accumulator: BTreeMap<usize, (f64, f64)>,
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
    let input = options
        .input
        .as_ref()
        .expect("parse_args validates DIA input");
    let request = options
        .request
        .as_ref()
        .expect("parse_args validates DIA request");
    if !options.mz_targets.is_empty() {
        return run_multi_target_slices(input, request, &options);
    }
    let default_prefix = default_out_prefix(input, request);
    let file_prefix = options.out_prefix.clone().unwrap_or(default_prefix);
    let out_prefix = options
        .outdir
        .as_ref()
        .map(|outdir| outdir.join(&file_prefix))
        .unwrap_or(file_prefix);
    options.verbosity.status(format_args!(
        "DIA slice: input={} target={} output-prefix={}",
        compact_display_path(input.path(), 72),
        request_target_label(request),
        out_prefix.display()
    ));
    options.verbosity.detail(format_args!(
        "DIA slice: m/z {:.6}-{:.6}, RT {}, IM {:?}-{:?}, quad {:?}-{:?}",
        request.mz_bounds().0,
        request.mz_bounds().1,
        request.rt_window_label(),
        request.im_min,
        request.im_max,
        request.quad_min,
        request.quad_max,
    ));

    match input {
        DiaInput::Mzml(path) => run_mzml_slice(
            path,
            request,
            &out_prefix,
            options.mz_profile_bins,
            &PseudoMs2Options {
                enabled: options.pseudo_ms2,
                rt_window_min: options.pseudo_ms2_rt_window_min,
            },
            options.rt_smooth,
            options.rt_smooth_window,
            options.verbosity,
            DiaCapabilities {
                has_mobility: false,
                has_isolation_window: false,
                requires_vendor_runtime: false,
            },
        ),
        DiaInput::Bruker(path) => run_bruker_slice(
            path,
            request,
            options.bruker_so.as_deref(),
            options.python_bin.as_deref(),
            options.bruker_backend,
            options.mz_profile_bins,
            &PseudoMs2Options {
                enabled: options.pseudo_ms2,
                rt_window_min: options.pseudo_ms2_rt_window_min,
            },
            options.rt_smooth,
            options.rt_smooth_window,
            options.verbosity,
            &out_prefix,
        ),
    }
}

fn run_multi_target_slices(
    input: &DiaInput,
    base_request: &DiaSliceRequest,
    options: &DiaSliceOptions,
) -> anyhow::Result<()> {
    let DiaInput::Bruker(path) = input else {
        anyhow::bail!("repeatable --target currently supports native Bruker .d input only");
    };
    if options.bruker_backend == BrukerBackend::AlphaTims {
        anyhow::bail!("repeatable --target requires --bruker-backend native or auto");
    }

    let requests = options
        .mz_targets
        .iter()
        .map(|target| request_for_mz_target(base_request, target))
        .collect::<Vec<_>>();
    let out_prefixes = options
        .mz_targets
        .iter()
        .map(|target| {
            multi_target_out_prefix(
                input,
                options.outdir.as_ref(),
                options.out_prefix.as_ref(),
                target,
            )
        })
        .collect::<Vec<_>>();
    let target_labels = options
        .mz_targets
        .iter()
        .map(|target| format!("{}={:.4}", target.label, target.mz))
        .collect::<Vec<_>>()
        .join(", ");
    options.verbosity.status(format_args!(
        "DIA slice: input={} targets={} output-count={}",
        compact_display_path(input.path(), 72),
        target_labels,
        requests.len()
    ));
    run_bruker_multi_target_timsrust(
        path,
        &requests,
        &out_prefixes,
        &multi_target_run_tic_prefix(input, options.outdir.as_ref(), options.out_prefix.as_ref()),
        options.mz_profile_bins,
        options.rt_smooth,
        options.rt_smooth_window,
        options.verbosity,
    )
    .with_context(|| "native timsrust multi-target extraction failed")
}

fn print_help() {
    let program = crate::program_name();
    println!("{program} dia-slice");
    println!();
    println!("USAGE:");
    println!(
        "  {program} dia-slice (--mzml <file> | --bruker <run.d>) (--mz <center> | --peptide <SEQ> | --target <label:mz>...) [options]"
    );
    println!();
    println!("OPTIONS:");
    println!("  --mzml <file>            Input DIA mzML file");
    println!("  --bruker <run.d>         Input Bruker .d folder");
    println!("  --bruker-so <path>       Override timsdata.so location for Bruker .d input");
    println!("  --bruker-backend <name>  Bruker backend: auto, native, alphatims [auto]");
    println!("  --python <exe>           Python interpreter for the alphaTims Bruker bridge");
    println!("  --mz <center>            Target fragment/signal m/z center");
    println!(
        "  --target <label:mz>      Repeatable raw m/z target for one-pass native Bruker extraction"
    );
    println!("  --peptide <SEQ>          Peptide precursor target; supports inline mods and /charge suffix");
    println!("  --sequence <SEQ>         Alias for --peptide");
    println!("  --sequence-modi <SEQ>    Alias for --peptide");
    println!("  --sequencemodi <SEQ>     Alias for --peptide");
    println!(
        "  --fragment <ion>         Peptide fragment target, e.g. b8, y8, b8++ [default: precursor]"
    );
    println!(
        "  --neutral-losses       Enable residue-aware -H2O / -NH3 and phospho -H3PO4 fragments"
    );
    println!(
        "  --pseudo-ms2             Write aggregated DIA fragment evidence as pseudo-MS2 TSV/SVG"
    );
    println!("  --ms2                    Alias for --pseudo-ms2");
    println!(
        "  --pseudo-ms2-rt-window <min>  Inferred RT window width [{}]",
        DEFAULT_PSEUDO_MS2_RT_WINDOW_MIN
    );
    println!("  --mod <pos:delta>        Peptide mass shift, repeatable; requires --peptide");
    println!(
        "  --charge <int>           Peptide precursor charge override [default: peptide /charge or {}]",
        DEFAULT_PEPTIDE_CHARGE
    );
    println!(
        "  --mz-ppm <ppm>           m/z extraction tolerance in ppm [{}]",
        DEFAULT_MZ_PPM
    );
    println!("  --mz-da <da>             Absolute m/z tolerance in Th; overrides --mz-ppm [default: none]");
    println!("  --rt-min <min>           RT lower bound in minutes [default: none; all RT]");
    println!("  --rt-max <max>           RT upper bound in minutes [default: none; all RT]");
    println!("  --im-min <1/K0>          Ion mobility lower bound [default: none; all mobility]");
    println!("  --im-max <1/K0>          Ion mobility upper bound [default: none; all mobility]");
    println!(
        "  --quad-min <mz>          Quadrupole isolation lower bound [default: none; all windows]"
    );
    println!(
        "  --quad-max <mz>          Quadrupole isolation upper bound [default: none; all windows]"
    );
    println!("  --out-prefix <name>      Output file stem/name only; no directories");
    println!("  --outdir <dir>           Output directory using the default generated file prefix");
    println!(
        "  --mz-bins <n>            Number of bins for the m/z profile [{}]",
        DEFAULT_MZ_PROFILE_BINS
    );
    println!("  --rt-smooth              Draw a smoothed overlay on the RT profile");
    println!(
        "  --rt-smooth-window <n>   Gaussian smoothing window in profile points [{}]; implies --rt-smooth",
        DEFAULT_RT_SMOOTH_WINDOW_POINTS
    );
    println!("  -v, --verbose            Print detailed progress/status messages");
    println!("  -q, --quiet              Suppress non-error terminal output");
    println!("  --help                   Show this help");
    println!();
    println!("OUTPUTS:");
    println!("  <prefix>.summary.txt");
    println!("  <prefix>.rt_profile.tsv");
    println!("  <prefix>.mz_profile.tsv");
    println!("  <prefix>.im_profile.tsv  Mobility-capable backends only");
    println!("  <prefix>.run_tic.json    Native Bruker run-level MS1/MS2 TIC summary");
    println!("  <prefix>.svg");
    println!();
    println!("NOTES:");
    println!(
        "  Native Bruker .d support uses timsrust. The alphaTims fallback requires timsdata.so and checks {}",
        DEFAULT_BRUKER_SO_PATHS.join(", ")
    );
    println!("  Use --bruker-so or {BRUKER_SO_ENV_VAR} to override the default runtime path");
}

fn parse_args(args: Vec<String>) -> anyhow::Result<DiaSliceOptions> {
    let mut options = DiaSliceOptions::default();
    let mut mz = None::<f64>;
    let mut mz_ppm = DEFAULT_MZ_PPM;
    let mut mz_da = None::<f64>;
    let mut rt_min = None::<f64>;
    let mut rt_max = None::<f64>;
    let mut im_min = None::<f64>;
    let mut im_max = None::<f64>;
    let mut quad_min = None::<f64>;
    let mut quad_max = None::<f64>;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                let path = PathBuf::from(iter.next().context("--mzml expects a path")?);
                set_input(&mut options.input, DiaInput::Mzml(path), "--mzml")?;
            }
            "--bruker" => {
                let path = PathBuf::from(iter.next().context("--bruker expects a path")?);
                set_input(&mut options.input, DiaInput::Bruker(path), "--bruker")?;
            }
            "--bruker-so" => {
                options.bruker_so = Some(PathBuf::from(
                    iter.next().context("--bruker-so expects a path")?,
                ));
            }
            "--bruker-backend" => {
                let raw = iter.next().context("--bruker-backend expects a name")?;
                options.bruker_backend = BrukerBackend::parse(&raw)?;
            }
            "--python" => {
                options.python_bin = Some(PathBuf::from(
                    iter.next()
                        .context("--python expects a path or executable name")?,
                ));
            }
            "--mz" => {
                mz = Some(parse_f64_flag("--mz", iter.next())?);
            }
            "--target" | "--mz-target" => {
                let raw = iter
                    .next()
                    .with_context(|| format!("{arg} expects <label:mz>"))?;
                options.mz_targets.push(parse_mz_target(&raw)?);
            }
            "--peptide" | "--sequence" | "--sequence-modi" | "--sequencemodi" => {
                let peptide = iter
                    .next()
                    .with_context(|| format!("{arg} expects a peptide sequence"))?;
                set_peptide_input(&mut options, peptide)?;
            }
            "--fragment" => {
                let fragment = iter
                    .next()
                    .context("--fragment expects a fragment ion label")?;
                set_fragment_input(&mut options, fragment)?;
            }
            "--neutral-losses" => {
                options.neutral_losses_enabled = true;
            }
            "--pseudo-ms2" | "--ms2" => {
                options.pseudo_ms2 = true;
            }
            "--pseudo-ms2-rt-window" => {
                options.pseudo_ms2_rt_window_min =
                    parse_f64_flag("--pseudo-ms2-rt-window", iter.next())?;
            }
            "--mod" => {
                options
                    .mod_inputs
                    .push(iter.next().context("--mod expects <position>:<delta>")?);
            }
            "--charge" => {
                let raw = iter.next().context("--charge expects an integer")?;
                let charge = raw.parse::<i32>().context("invalid --charge")?;
                if charge <= 0 {
                    anyhow::bail!("--charge must be a positive integer");
                }
                options.charge_override = Some(charge);
            }
            "--mz-ppm" => {
                mz_ppm = parse_f64_flag("--mz-ppm", iter.next())?;
            }
            "--mz-da" => {
                mz_da = Some(parse_f64_flag("--mz-da", iter.next())?);
            }
            "--rt-min" => {
                rt_min = Some(parse_f64_flag("--rt-min", iter.next())?);
            }
            "--rt-max" => {
                rt_max = Some(parse_f64_flag("--rt-max", iter.next())?);
            }
            "--im-min" => {
                im_min = Some(parse_f64_flag("--im-min", iter.next())?);
            }
            "--im-max" => {
                im_max = Some(parse_f64_flag("--im-max", iter.next())?);
            }
            "--quad-min" => {
                quad_min = Some(parse_f64_flag("--quad-min", iter.next())?);
            }
            "--quad-max" => {
                quad_max = Some(parse_f64_flag("--quad-max", iter.next())?);
            }
            "--out-prefix" => {
                let raw = iter.next().context("--out-prefix expects a name")?;
                options.out_prefix = Some(parse_out_prefix(&raw)?);
            }
            "--outdir" => {
                options.outdir = Some(PathBuf::from(
                    iter.next().context("--outdir expects a path")?,
                ));
            }
            "--mz-bins" => {
                let raw = iter.next().context("--mz-bins expects an integer")?;
                let bins = raw.parse::<usize>().context("invalid --mz-bins")?;
                if bins == 0 {
                    anyhow::bail!("--mz-bins must be at least 1");
                }
                options.mz_profile_bins = bins;
            }
            "--rt-smooth" => {
                options.rt_smooth = true;
            }
            "--rt-smooth-window" => {
                let raw = iter
                    .next()
                    .context("--rt-smooth-window expects an integer")?;
                let window = raw.parse::<usize>().context("invalid --rt-smooth-window")?;
                if window == 0 {
                    anyhow::bail!("--rt-smooth-window must be at least 1");
                }
                options.rt_smooth = true;
                options.rt_smooth_window = window;
            }
            "-v" | "--verbose" => {
                options.verbosity = Verbosity::Verbose;
            }
            "-q" | "--quiet" => {
                options.verbosity = Verbosity::Quiet;
            }
            other => anyhow::bail!("unknown dia-slice option `{other}`"),
        }
    }

    if mz_ppm <= 0.0 || !mz_ppm.is_finite() {
        anyhow::bail!("--mz-ppm must be a positive finite number");
    }
    if let Some(value) = mz_da {
        if value <= 0.0 || !value.is_finite() {
            anyhow::bail!("--mz-da must be a positive finite number");
        }
    }
    if options.pseudo_ms2_rt_window_min <= 0.0 || !options.pseudo_ms2_rt_window_min.is_finite() {
        anyhow::bail!("--pseudo-ms2-rt-window must be a positive finite number");
    }
    let input = options.input.as_ref().ok_or_else(|| {
        anyhow::anyhow!("dia-slice requires exactly one of --mzml <file> or --bruker <run.d>")
    })?;
    if !options.mz_targets.is_empty() {
        if mz.is_some() {
            anyhow::bail!("specify only one of --mz <center> or repeatable --target <label:mz>");
        }
        if options.peptide_input.is_some() {
            anyhow::bail!("--target cannot be combined with --peptide");
        }
        if matches!(input, DiaInput::Mzml(_)) {
            anyhow::bail!("repeatable --target currently supports native Bruker .d input only");
        }
        if options.bruker_backend == BrukerBackend::AlphaTims {
            anyhow::bail!("repeatable --target requires --bruker-backend native or auto");
        }
        let mut seen = BTreeSet::new();
        for target in &options.mz_targets {
            if !seen.insert(target.label.clone()) {
                anyhow::bail!("duplicate --target label `{}`", target.label);
            }
        }
    }
    if options.peptide_input.is_none() && !options.mod_inputs.is_empty() {
        anyhow::bail!("--mod requires --peptide (or a sequence alias)");
    }
    if options.peptide_input.is_none() && options.charge_override.is_some() {
        anyhow::bail!("--charge requires --peptide (or a sequence alias)");
    }
    if options.peptide_input.is_none() && options.fragment_input.is_some() {
        anyhow::bail!("--fragment requires --peptide (or a sequence alias)");
    }
    if options.peptide_input.is_none() && options.neutral_losses_enabled {
        anyhow::bail!("--neutral-losses requires --peptide (or a sequence alias)");
    }
    if options.peptide_input.is_none() && options.pseudo_ms2 {
        anyhow::bail!("--pseudo-ms2 requires --peptide (or a sequence alias)");
    }
    if options.pseudo_ms2 && options.fragment_input.is_some() {
        anyhow::bail!("--pseudo-ms2 aggregates all peptide fragments; omit --fragment");
    }
    let neutral_losses: &[NeutralLossKind] = if options.neutral_losses_enabled {
        COMMON_NEUTRAL_LOSSES.as_slice()
    } else {
        &[]
    };
    let peptide_target = options
        .peptide_input
        .as_deref()
        .map(|peptide| {
            resolve_peptide_target(
                peptide,
                options.fragment_input.as_deref(),
                &options.mod_inputs,
                options.charge_override,
                neutral_losses,
            )
        })
        .transpose()?;
    let mz = match (mz, peptide_target.as_ref(), options.mz_targets.first()) {
        (Some(_), Some(_), _) => {
            anyhow::bail!("specify only one of --mz <center> or --peptide <SEQ>")
        }
        (Some(_), None, Some(_)) => {
            anyhow::bail!("specify only one of --mz <center> or repeatable --target <label:mz>")
        }
        (None, Some(_), Some(_)) => anyhow::bail!("--target cannot be combined with --peptide"),
        (Some(value), None, None) => value,
        (None, Some(target), None) => target
            .fragment
            .as_ref()
            .map(|fragment| fragment.mz)
            .unwrap_or(target.precursor_mz),
        (None, None, Some(target)) => target.mz,
        (None, None, None) => {
            anyhow::bail!(
                "dia-slice requires --mz <center>, --peptide <SEQ>, or --target <label:mz>"
            )
        }
    };
    if !mz.is_finite() {
        anyhow::bail!("--mz must be finite");
    }

    validate_bounds(rt_min, rt_max, "RT")?;
    validate_bounds(im_min, im_max, "ion mobility")?;
    validate_bounds(quad_min, quad_max, "quadrupole")?;

    if matches!(input, DiaInput::Mzml(_)) && (im_min.is_some() || quad_min.is_some()) {
        anyhow::bail!(
            "--mzml input does not support --im-* or --quad-* yet; those filters are reserved for backends with mobility/isolation-window dimensions"
        );
    }
    if options.bruker_so.is_some() && !matches!(input, DiaInput::Bruker(_)) {
        anyhow::bail!("--bruker-so is only valid together with --bruker <run.d>");
    }
    if options.python_bin.is_some() && !matches!(input, DiaInput::Bruker(_)) {
        anyhow::bail!("--python is only valid together with --bruker <run.d>");
    }
    if options.pseudo_ms2 && matches!(input, DiaInput::Mzml(_)) {
        anyhow::bail!("--pseudo-ms2 currently supports native Bruker .d input only");
    }
    if options.pseudo_ms2 && options.bruker_backend == BrukerBackend::AlphaTims {
        anyhow::bail!("--pseudo-ms2 currently requires --bruker-backend native or auto");
    }

    options.request = Some(DiaSliceRequest {
        mz,
        mz_ppm,
        mz_da,
        peptide_target,
        rt_min,
        rt_max,
        im_min,
        im_max,
        quad_min,
        quad_max,
    });
    Ok(options)
}

fn set_peptide_input(options: &mut DiaSliceOptions, peptide: String) -> anyhow::Result<()> {
    if options.peptide_input.is_some() {
        anyhow::bail!("specify peptide input only once");
    }
    options.peptide_input = Some(peptide);
    Ok(())
}

fn set_fragment_input(options: &mut DiaSliceOptions, fragment: String) -> anyhow::Result<()> {
    if options.fragment_input.is_some() {
        anyhow::bail!("specify fragment input only once");
    }
    options.fragment_input = Some(fragment);
    Ok(())
}

fn resolve_peptide_target(
    peptide_input: &str,
    fragment_input: Option<&str>,
    mod_inputs: &[String],
    charge_override: Option<i32>,
    neutral_losses: &[NeutralLossKind],
) -> anyhow::Result<DiaPeptideTarget> {
    let context = prepare_annotation(
        peptide_input,
        mod_inputs,
        neutral_losses,
        charge_override,
        MassTolerance::Ppm(DEFAULT_MZ_PPM),
    )?;
    let charge = context.charge_context.unwrap_or(DEFAULT_PEPTIDE_CHARGE);
    if charge <= 0 {
        anyhow::bail!("peptide target charge must be positive");
    }
    let precursor_mz = context
        .peptide
        .precursor_mz(charge)
        .ok_or_else(|| anyhow::anyhow!("failed to compute peptide precursor m/z"))?;
    let charges = fragment_charge_states(Some(charge));
    let fragments = generate_fragments(&context.peptide, &charges, neutral_losses)
        .into_iter()
        .map(dia_fragment_target)
        .collect::<Vec<_>>();
    let fragment = fragment_input
        .map(|fragment| resolve_fragment_target(fragment, &fragments, context.peptide.sequence()))
        .transpose()?;
    Ok(DiaPeptideTarget {
        input: peptide_input.to_string(),
        sequence: context.peptide.sequence().to_string(),
        modified_sequence: context.modified_sequence(),
        charge,
        precursor_mz,
        fragment,
        fragments,
    })
}

fn resolve_fragment_target(
    fragment_input: &str,
    fragments: &[DiaFragmentTarget],
    peptide_sequence: &str,
) -> anyhow::Result<DiaFragmentTarget> {
    let requested = fragment_input.trim();
    if requested.is_empty() {
        anyhow::bail!("--fragment cannot be empty");
    }
    let matches = fragments
        .iter()
        .filter(|fragment| fragment_label_matches(fragment, requested))
        .collect::<Vec<_>>();

    match matches.as_slice() {
        [fragment] => {
            let mut out = (*fragment).clone();
            out.input = requested.to_string();
            Ok(out)
        }
        [] => {
            let available = fragments
                .iter()
                .map(|fragment| fragment.label.clone())
                .collect::<Vec<_>>()
                .join(", ");
            anyhow::bail!(
                "unknown --fragment `{requested}` for peptide {}; available fragments: {}",
                peptide_sequence,
                available
            )
        }
        _ => anyhow::bail!("ambiguous --fragment `{requested}`"),
    }
}

fn dia_fragment_target(fragment: FragmentIon) -> DiaFragmentTarget {
    let label = fragment.label();
    DiaFragmentTarget {
        input: label.clone(),
        label,
        series: fragment_series_label(fragment.series).to_string(),
        cleavage_index: fragment.cleavage_index,
        neutral_loss: fragment.neutral_loss,
        charge: fragment.charge,
        mz: fragment.theoretical_mz,
    }
}

fn fragment_series_label(series: FragmentSeries) -> &'static str {
    match series {
        FragmentSeries::B => "b",
        FragmentSeries::Y => "y",
    }
}

fn fragment_label_matches(fragment: &DiaFragmentTarget, requested: &str) -> bool {
    if fragment.label.eq_ignore_ascii_case(requested) {
        return true;
    }
    if fragment.charge == 1 && requested.ends_with('+') {
        let trimmed = requested.trim_end_matches('+');
        return fragment.label.eq_ignore_ascii_case(trimmed);
    }
    false
}

fn set_input(slot: &mut Option<DiaInput>, value: DiaInput, flag: &str) -> anyhow::Result<()> {
    if slot.is_some() {
        anyhow::bail!("dia-slice accepts only one input source; `{flag}` conflicts with the earlier input flag");
    }
    *slot = Some(value);
    Ok(())
}

fn parse_f64_flag(flag: &str, value: Option<String>) -> anyhow::Result<f64> {
    let raw = value.with_context(|| format!("{flag} expects a value"))?;
    raw.parse::<f64>()
        .with_context(|| format!("invalid {flag} value `{raw}`"))
}

fn parse_mz_target(raw: &str) -> anyhow::Result<DiaMzTarget> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        anyhow::bail!("--target cannot be empty");
    }
    let (label_raw, mz_raw) = match trimmed.rsplit_once(':') {
        Some((label, mz)) => {
            if label.trim().is_empty() {
                anyhow::bail!("--target label cannot be empty in `{raw}`");
            }
            (Some(label.trim()), mz.trim())
        }
        None => (None, trimmed),
    };
    let mz = mz_raw
        .parse::<f64>()
        .with_context(|| format!("invalid --target m/z value `{mz_raw}`"))?;
    if mz <= 0.0 || !mz.is_finite() {
        anyhow::bail!("--target m/z must be a positive finite number");
    }
    let label = label_raw
        .map(sanitize_filename_component)
        .unwrap_or_else(|| default_mz_target_label(mz));
    Ok(DiaMzTarget { label, mz })
}

fn default_mz_target_label(mz: f64) -> String {
    format!("mz_{mz:.4}")
}

fn parse_out_prefix(raw: &str) -> anyhow::Result<PathBuf> {
    let path = Path::new(raw);
    let has_one_normal_component = path.components().count() == 1 && path.file_name().is_some();
    if raw.ends_with('/') || raw.ends_with('\\') || !has_one_normal_component {
        anyhow::bail!(
            "--out-prefix expects a file stem/name, not a path; use --outdir for output directories"
        );
    }
    Ok(PathBuf::from(raw))
}

fn validate_bounds(min: Option<f64>, max: Option<f64>, label: &str) -> anyhow::Result<()> {
    if min.is_some() != max.is_some() {
        anyhow::bail!("{label} bounds must specify both min and max");
    }
    if let (Some(lo), Some(hi)) = (min, max) {
        if !lo.is_finite() || !hi.is_finite() {
            anyhow::bail!("{label} bounds must be finite");
        }
        if lo >= hi {
            anyhow::bail!("{label} min must be smaller than max");
        }
    }
    Ok(())
}

fn default_out_prefix(input: &DiaInput, request: &DiaSliceRequest) -> PathBuf {
    let target = request
        .peptide_target
        .as_ref()
        .map(|target| {
            let fragment = target
                .fragment
                .as_ref()
                .map(|fragment| format!("_{}", sanitize_filename_component(&fragment.label)))
                .unwrap_or_default();
            format!(
                "{}_z{}{}",
                sanitize_filename_component(&target.sequence),
                target.charge,
                fragment
            )
        })
        .unwrap_or_else(|| format!("mz_{:.4}", request.mz));
    PathBuf::from(format!("{}.dia_slice_{}", input.default_stem(), target))
}

fn request_for_mz_target(base_request: &DiaSliceRequest, target: &DiaMzTarget) -> DiaSliceRequest {
    let mut request = base_request.clone();
    request.mz = target.mz;
    request.peptide_target = None;
    request
}

fn multi_target_out_prefix(
    input: &DiaInput,
    outdir: Option<&PathBuf>,
    out_prefix: Option<&PathBuf>,
    target: &DiaMzTarget,
) -> PathBuf {
    let file_prefix = out_prefix
        .map(|base| {
            let base = base
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("dia_slice");
            PathBuf::from(format!("{base}__{}", target.label))
        })
        .unwrap_or_else(|| {
            PathBuf::from(format!(
                "{}.dia_slice_{}",
                input.default_stem(),
                target.label
            ))
        });
    outdir
        .map(|outdir| outdir.join(&file_prefix))
        .unwrap_or(file_prefix)
}

fn multi_target_run_tic_prefix(
    input: &DiaInput,
    outdir: Option<&PathBuf>,
    out_prefix: Option<&PathBuf>,
) -> PathBuf {
    let file_prefix = out_prefix
        .map(|base| {
            base.file_name()
                .and_then(|value| value.to_str())
                .filter(|value| !value.is_empty())
                .unwrap_or("dia_slice")
                .to_string()
        })
        .unwrap_or_else(|| input.default_stem());
    outdir
        .map(|outdir| outdir.join(&file_prefix))
        .unwrap_or_else(|| PathBuf::from(file_prefix))
}

fn request_target_label(request: &DiaSliceRequest) -> String {
    request
        .peptide_target
        .as_ref()
        .map(|target| match target.fragment.as_ref() {
            Some(fragment) => format!(
                "{}/{} fragment {} m/z {:.4}",
                target.modified_sequence, target.charge, fragment.label, fragment.mz
            ),
            None => format!(
                "{}/{} precursor m/z {:.4}",
                target.modified_sequence, target.charge, target.precursor_mz
            ),
        })
        .unwrap_or_else(|| format!("m/z {:.4}", request.mz))
}

fn run_mzml_slice(
    path: &Path,
    request: &DiaSliceRequest,
    out_prefix: &Path,
    mz_profile_bins: usize,
    pseudo_ms2: &PseudoMs2Options,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
    capabilities: DiaCapabilities,
) -> anyhow::Result<()> {
    if pseudo_ms2.enabled {
        anyhow::bail!("--pseudo-ms2 currently supports native Bruker .d input only");
    }
    if !path.is_file() {
        anyhow::bail!("mzML input does not exist: {}", path.display());
    }

    let (mz_min, mz_max) = request.mz_bounds();
    verbosity.status(format_args!(
        "DIA slice: reading mzML and extracting MS2 signal in m/z {:.6}-{:.6}",
        mz_min, mz_max
    ));
    let mut reader = open_reader(path)?;
    let total = reader.len();
    verbosity.detail(format_args!("DIA slice: mzML spectrum count={total}"));
    let mut rt_rows = Vec::<RtProfileRow>::new();
    let mut mz_bins = vec![
        MzProfileBin {
            mz_center: 0.0,
            summed_intensity: 0.0,
        };
        mz_profile_bins
    ];
    initialize_mz_bins(&mut mz_bins, mz_min, mz_max);

    let mut spectra_considered = 0usize;
    let mut spectra_with_signal = 0usize;
    let mut matched_peaks = 0usize;

    for idx in 0..total {
        let spectrum = load_spectrum_by_index(&mut reader, idx as u32)?;
        if spectrum.meta.ms_level != 2 {
            continue;
        }

        let rt = spectrum.meta.rt_minutes.map(f64::from);
        if !rt_in_window(rt, request.rt_min, request.rt_max) {
            continue;
        }

        spectra_considered += 1;
        let mut sum_intensity = 0.0f64;
        let mut peaks_in_window = 0usize;
        for (&mz, &intensity) in spectrum.mz.iter().zip(spectrum.intensity.iter()) {
            if mz < mz_min || mz > mz_max {
                continue;
            }
            let intensity = intensity as f64;
            sum_intensity += intensity;
            peaks_in_window += 1;
            accumulate_mz_bin(&mut mz_bins, mz_min, mz_max, mz, intensity);
        }

        if peaks_in_window > 0 {
            spectra_with_signal += 1;
            matched_peaks += peaks_in_window;
        }

        rt_rows.push(RtProfileRow {
            scan_index: spectrum.meta.idx,
            scan_id: spectrum.meta.scan_id,
            rt_minutes: rt,
            summed_intensity: sum_intensity,
            matched_peaks: peaks_in_window,
            precursor_mz: spectrum.meta.precursor_mz,
        });
    }

    if rt_rows.is_empty() {
        anyhow::bail!(
            "no MS2 spectra matched the requested mzML DIA slice window (mz {:.4}-{:.4}, RT {})",
            mz_min,
            mz_max,
            request.rt_window_label(),
        );
    }

    let summary = DiaSliceSummary {
        backend_label: "mzML",
        input_path: path.to_path_buf(),
        out_prefix: out_prefix.to_path_buf(),
        mz_min,
        mz_max,
        spectra_considered,
        spectra_with_signal,
        matched_peaks,
        capabilities,
        acquisition_mode: None,
        intensity_column: None,
        vendor_runtime_path: None,
    };
    write_outputs(
        &summary,
        request,
        &rt_rows,
        &mz_bins,
        None,
        rt_smooth,
        rt_smooth_window,
    )?;

    verbosity.success(format_args!(
        "Wrote DIA slice outputs: {} ({} spectra considered, {} with signal, {} matched peaks)",
        out_prefix.display(),
        summary.spectra_considered,
        summary.spectra_with_signal,
        summary.matched_peaks,
    ));
    Ok(())
}

fn run_bruker_slice(
    path: &Path,
    request: &DiaSliceRequest,
    bruker_so_override: Option<&Path>,
    python_override: Option<&Path>,
    backend: BrukerBackend,
    mz_profile_bins: usize,
    pseudo_ms2: &PseudoMs2Options,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
    out_prefix: &Path,
) -> anyhow::Result<()> {
    if !path.is_dir() {
        anyhow::bail!(
            "Bruker input does not exist or is not a directory: {}",
            path.display()
        );
    }

    match backend {
        BrukerBackend::Native => run_bruker_slice_timsrust(
            path,
            request,
            mz_profile_bins,
            pseudo_ms2,
            rt_smooth,
            rt_smooth_window,
            verbosity,
            out_prefix,
        ),
        BrukerBackend::AlphaTims => run_bruker_slice_alphatims(
            path,
            request,
            bruker_so_override,
            python_override,
            mz_profile_bins,
            rt_smooth,
            rt_smooth_window,
            verbosity,
            out_prefix,
        ),
        BrukerBackend::Auto => {
            match run_bruker_slice_timsrust(
                path,
                request,
                mz_profile_bins,
                pseudo_ms2,
                rt_smooth,
                rt_smooth_window,
                verbosity,
                out_prefix,
            ) {
                Ok(()) => Ok(()),
                Err(native_err) => {
                    if pseudo_ms2.enabled {
                        return Err(native_err)
                            .context("--pseudo-ms2 requires the native timsrust Bruker backend");
                    }
                    verbosity.status(format_args!(
                        "Native timsrust Bruker backend failed; falling back to alphaTims bridge ({native_err:#})"
                    ));
                    run_bruker_slice_alphatims(
                        path,
                        request,
                        bruker_so_override,
                        python_override,
                        mz_profile_bins,
                        rt_smooth,
                        rt_smooth_window,
                        verbosity,
                        out_prefix,
                    )
                    .with_context(|| {
                        format!(
                            "native timsrust backend failed first ({native_err:#}); alphaTims fallback also failed"
                        )
                    })
                }
            }
        }
    }
}

fn run_bruker_slice_timsrust(
    path: &Path,
    request: &DiaSliceRequest,
    mz_profile_bins: usize,
    pseudo_ms2: &PseudoMs2Options,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
    out_prefix: &Path,
) -> anyhow::Result<()> {
    let (mz_min, mz_max) = request.mz_bounds();
    verbosity.status(format_args!(
        "DIA slice: reading Bruker .d with timsrust and extracting m/z {:.6}-{:.6}",
        mz_min, mz_max
    ));
    let (payload, run_tic) =
        build_timsrust_bruker_payload(path, request, mz_profile_bins, verbosity)?;
    verbosity.detail(format_args!(
        "DIA slice: native extraction yielded {} frames with signal from {} considered frames",
        payload.frames_with_signal, payload.frames_considered
    ));
    let precursor_rt_points = payload
        .rt_profile
        .iter()
        .map(|row| (row.rt_minutes, row.summed_intensity))
        .collect::<Vec<_>>();
    let has_isolation_window = payload.acquisition_mode == "diaPASEF";
    let summary = finish_bruker_slice_outputs(
        path,
        request,
        out_prefix,
        "Bruker .d (timsrust native)",
        has_isolation_window,
        false,
        None,
        payload,
        false,
        rt_smooth,
        rt_smooth_window,
    )?;
    write_bruker_run_tic_json(out_prefix, &run_tic)?;
    if pseudo_ms2.enabled {
        verbosity.status(format_args!(
            "Pseudo-MS2: inferring RT window from the extracted precursor trace"
        ));
        let report = build_timsrust_pseudo_ms2(
            path,
            request,
            out_prefix,
            pseudo_ms2,
            &precursor_rt_points,
            verbosity,
        )?;
        write_pseudo_ms2_outputs(&report)?;
        verbosity.success(format_args!(
            "Wrote pseudo-MS2 outputs: {} ({} frames considered, {} with signal, {} matched detector events)",
            out_prefix.display(),
            report.frames_considered,
            report.frames_with_signal,
            report.matched_events,
        ));
    }
    print_bruker_success(out_prefix, &summary, verbosity);
    Ok(())
}

fn run_bruker_multi_target_timsrust(
    path: &Path,
    requests: &[DiaSliceRequest],
    out_prefixes: &[PathBuf],
    run_tic_prefix: &Path,
    mz_profile_bins: usize,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
) -> anyhow::Result<()> {
    if !path.is_dir() {
        anyhow::bail!(
            "Bruker input does not exist or is not a directory: {}",
            path.display()
        );
    }
    if requests.len() != out_prefixes.len() {
        anyhow::bail!("internal error: target/output count mismatch");
    }
    let ranges = requests
        .iter()
        .map(|request| {
            let (mz_min, mz_max) = request.mz_bounds();
            format!("{mz_min:.4}-{mz_max:.4}")
        })
        .collect::<Vec<_>>()
        .join(", ");
    verbosity.status(format_args!(
        "DIA slice: reading Bruker .d with timsrust and extracting {} m/z target(s): {}",
        requests.len(),
        ranges
    ));
    let extraction = build_timsrust_bruker_payloads(path, requests, mz_profile_bins, verbosity)?;
    write_bruker_run_tic_json(run_tic_prefix, &extraction.run_tic)?;
    let mut summaries = Vec::with_capacity(extraction.payloads.len());
    for ((request, out_prefix), payload) in
        requests.iter().zip(out_prefixes).zip(extraction.payloads)
    {
        let has_isolation_window = payload.acquisition_mode == "diaPASEF";
        let summary = finish_bruker_slice_outputs(
            path,
            request,
            out_prefix,
            "Bruker .d (timsrust native)",
            has_isolation_window,
            false,
            None,
            payload,
            true,
            rt_smooth,
            rt_smooth_window,
        )?;
        summaries.push(summary);
    }
    let matched_events = summaries
        .iter()
        .map(|summary| summary.matched_peaks)
        .sum::<usize>();
    verbosity.success(format_args!(
        "Wrote DIA slice outputs for {} targets ({} total matched detector events)",
        summaries.len(),
        matched_events
    ));
    Ok(())
}

fn run_bruker_slice_alphatims(
    path: &Path,
    request: &DiaSliceRequest,
    bruker_so_override: Option<&Path>,
    python_override: Option<&Path>,
    mz_profile_bins: usize,
    rt_smooth: bool,
    rt_smooth_window: usize,
    verbosity: Verbosity,
    out_prefix: &Path,
) -> anyhow::Result<()> {
    let runtime = resolve_bruker_runtime(bruker_so_override)?;
    let python = resolve_python(python_override)?;
    let (mz_min, mz_max) = request.mz_bounds();
    verbosity.status(format_args!(
        "DIA slice: reading Bruker .d through alphaTims bridge and extracting m/z {:.6}-{:.6}",
        mz_min, mz_max
    ));
    verbosity.detail(format_args!(
        "DIA slice: alphaTims runtime={} python={}",
        runtime.path.display(),
        python.display()
    ));
    let payload = run_bruker_bridge(
        path,
        request,
        runtime.path.as_path(),
        python.as_path(),
        mz_profile_bins,
    )?;

    let summary = finish_bruker_slice_outputs(
        path,
        request,
        out_prefix,
        "Bruker .d (alphaTims bridge)",
        true,
        true,
        Some(runtime.path),
        payload,
        false,
        rt_smooth,
        rt_smooth_window,
    )?;
    print_bruker_success(out_prefix, &summary, verbosity);
    Ok(())
}

fn finish_bruker_slice_outputs(
    path: &Path,
    request: &DiaSliceRequest,
    out_prefix: &Path,
    backend_label: &'static str,
    has_isolation_window: bool,
    requires_vendor_runtime: bool,
    vendor_runtime_path: Option<PathBuf>,
    payload: BrukerBridgePayload,
    allow_empty: bool,
    rt_smooth: bool,
    rt_smooth_window: usize,
) -> anyhow::Result<DiaSliceSummary> {
    if payload.matched_events == 0 && !allow_empty {
        anyhow::bail!(
            "no detector events matched the requested Bruker slice window (mz {:.4}-{:.4}, RT {}, IM {:?}-{:?}, quad {:?}-{:?})",
            payload.mz_min,
            payload.mz_max,
            request.rt_window_label(),
            request.im_min,
            request.im_max,
            request.quad_min,
            request.quad_max,
        );
    }

    let rt_rows = payload
        .rt_profile
        .iter()
        .map(|row| RtProfileRow {
            scan_index: row.frame_index,
            scan_id: format!("frame={}", row.frame_index),
            rt_minutes: Some(row.rt_minutes),
            summed_intensity: row.summed_intensity,
            matched_peaks: row.matched_events,
            precursor_mz: None,
        })
        .collect::<Vec<_>>();
    let mz_bins = payload
        .mz_profile
        .iter()
        .map(|row| MzProfileBin {
            mz_center: row.mz_center,
            summed_intensity: row.summed_intensity,
        })
        .collect::<Vec<_>>();
    let im_rows = payload
        .im_profile
        .iter()
        .map(|row| MobilityProfileRow {
            mobility: row.mobility,
            summed_intensity: row.summed_intensity,
        })
        .collect::<Vec<_>>();

    let summary = DiaSliceSummary {
        backend_label,
        input_path: path.to_path_buf(),
        out_prefix: out_prefix.to_path_buf(),
        mz_min: payload.mz_min,
        mz_max: payload.mz_max,
        spectra_considered: payload.frames_considered,
        spectra_with_signal: payload.frames_with_signal,
        matched_peaks: payload.matched_events,
        capabilities: DiaCapabilities {
            has_mobility: true,
            has_isolation_window,
            requires_vendor_runtime,
        },
        acquisition_mode: Some(payload.acquisition_mode),
        intensity_column: Some(payload.intensity_column),
        vendor_runtime_path,
    };
    write_outputs(
        &summary,
        request,
        &rt_rows,
        &mz_bins,
        Some(&im_rows),
        rt_smooth,
        rt_smooth_window,
    )?;
    Ok(summary)
}

fn print_bruker_success(out_prefix: &Path, summary: &DiaSliceSummary, verbosity: Verbosity) {
    verbosity.success(format_args!(
        "Wrote DIA slice outputs: {} ({} frames considered, {} with signal, {} matched detector events)",
        out_prefix.display(),
        summary.spectra_considered,
        summary.spectra_with_signal,
        summary.matched_peaks,
    ));
}

fn collect_timsrust_ms2_frame_candidates<F>(
    frame_reader: &FrameReader,
    mut include_rt_seconds: F,
) -> anyhow::Result<Vec<TimsrustFrameCandidate>>
where
    F: FnMut(f64) -> bool,
{
    let mut candidates = Vec::new();
    for index in 0..frame_reader.len() {
        let frame_meta = frame_reader
            .get_frame_without_coordinates(index)
            .with_context(|| format!("failed to inspect Bruker frame {}", index + 1))?;
        if frame_meta.ms_level != MSLevel::MS2 {
            continue;
        }
        if !include_rt_seconds(frame_meta.rt_in_seconds) {
            continue;
        }
        candidates.push(TimsrustFrameCandidate {
            reader_index: index,
            rt_seconds: frame_meta.rt_in_seconds,
        });
    }
    Ok(candidates)
}

fn build_timsrust_bruker_payload(
    path: &Path,
    request: &DiaSliceRequest,
    mz_profile_bins: usize,
    verbosity: Verbosity,
) -> anyhow::Result<(BrukerBridgePayload, BrukerRunTicSummary)> {
    let mut extraction = build_timsrust_bruker_payloads(
        path,
        std::slice::from_ref(request),
        mz_profile_bins,
        verbosity,
    )?;
    let payload = extraction
        .payloads
        .pop()
        .ok_or_else(|| anyhow::anyhow!("native timsrust extraction produced no payload"))?;
    Ok((payload, extraction.run_tic))
}

fn build_timsrust_bruker_payloads(
    path: &Path,
    requests: &[DiaSliceRequest],
    mz_profile_bins: usize,
    verbosity: Verbosity,
) -> anyhow::Result<TimsrustExtraction> {
    let Some(base_request) = requests.first() else {
        anyhow::bail!("native timsrust extraction requires at least one m/z target");
    };
    let metadata = MetadataReader::new(path)
        .with_context(|| "failed to read Bruker metadata with timsrust")?;
    let frame_reader =
        FrameReader::new(path).with_context(|| "failed to open Bruker .d with timsrust")?;
    let acquisition = frame_reader.get_acquisition();
    if acquisition != AcquisitionType::DIAPASEF
        && (base_request.quad_min.is_some() || base_request.quad_max.is_some())
    {
        anyhow::bail!(
            "native timsrust Bruker backend supports --quad-* filters only for diaPASEF data"
        );
    }

    let mut accumulators = requests
        .iter()
        .map(|request| prepare_timsrust_slice_accumulator(&metadata, request, mz_profile_bins))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let total_frames = frame_reader.len();
    verbosity.detail(format_args!(
        "DIA slice: scanning {} total frames for run TIC and {} target(s)",
        total_frames,
        accumulators.len()
    ));
    let mut progress = ProgressReporter::new("DIA slice native frames", total_frames, verbosity);

    let mut run_tic_accumulator = BrukerRunTicAccumulator::new();
    let mut frames_considered = 0usize;
    for frame_index in 0..total_frames {
        let frame = frame_reader
            .get(frame_index)
            .with_context(|| format!("failed to read Bruker frame {}", frame_index + 1))?;
        run_tic_accumulator.add_frame(&frame);
        if frame.ms_level != MSLevel::MS2
            || !rt_seconds_in_window(
                frame.rt_in_seconds,
                base_request.rt_min,
                base_request.rt_max,
            )
        {
            progress.advance();
            continue;
        }
        frames_considered += 1;
        let mut frame_intensities = vec![0.0_f64; accumulators.len()];
        let mut frame_events = vec![0usize; accumulators.len()];

        for scan_index in 0..frame.scan_offsets.len().saturating_sub(1) {
            let mobility = metadata.im_converter.convert(scan_index as u32);
            if !mobility_in_window(mobility, base_request.im_min, base_request.im_max) {
                continue;
            }
            if !quad_in_window(
                &frame,
                scan_index,
                base_request.quad_min,
                base_request.quad_max,
            ) {
                continue;
            }

            let start = frame.scan_offsets[scan_index];
            let end = frame.scan_offsets[scan_index + 1];
            if start >= end {
                continue;
            }
            let scan_tofs = &frame.tof_indices[start..end];
            for (target_idx, accumulator) in accumulators.iter_mut().enumerate() {
                let peak_start = scan_tofs.partition_point(|tof| *tof < accumulator.tof_min);
                let peak_end = scan_tofs.partition_point(|tof| *tof <= accumulator.tof_max);

                for relative_peak_index in peak_start..peak_end {
                    let peak_index = start + relative_peak_index;
                    let mz = metadata.mz_converter.convert(frame.tof_indices[peak_index]);
                    if mz < accumulator.mz_min || mz > accumulator.mz_max {
                        continue;
                    }
                    let intensity = frame.get_corrected_intensity(peak_index);
                    if !intensity.is_finite() || intensity <= 0.0 {
                        continue;
                    }
                    accumulator.matched_events += 1;
                    frame_events[target_idx] += 1;
                    frame_intensities[target_idx] += intensity;
                    accumulate_mz_histogram_bin(
                        &mut accumulator.mz_bins,
                        accumulator.mz_min,
                        accumulator.mz_max,
                        mz,
                        intensity,
                    );
                    accumulator
                        .im_accumulator
                        .entry(scan_index)
                        .and_modify(|(_, summed_intensity)| *summed_intensity += intensity)
                        .or_insert((mobility, intensity));
                }
            }
        }

        for (target_idx, accumulator) in accumulators.iter_mut().enumerate() {
            if frame_events[target_idx] > 0 {
                accumulator.frames_with_signal += 1;
                accumulator.rt_profile.push(BrukerBridgeRtRow {
                    frame_index: frame.index as u32,
                    rt_minutes: frame.rt_in_seconds / 60.0,
                    summed_intensity: frame_intensities[target_idx],
                    matched_events: frame_events[target_idx],
                });
            }
        }
        progress.advance();
    }
    progress.finish();

    let acquisition_mode = timsrust_acquisition_label(acquisition).to_string();
    let payloads = accumulators
        .into_iter()
        .map(|accumulator| {
            finalize_timsrust_slice_payload(accumulator, &acquisition_mode, frames_considered)
        })
        .collect();
    Ok(TimsrustExtraction {
        payloads,
        run_tic: run_tic_accumulator.finalize(path, acquisition_mode),
    })
}

fn prepare_timsrust_slice_accumulator(
    metadata: &Metadata,
    request: &DiaSliceRequest,
    mz_profile_bins: usize,
) -> anyhow::Result<TimsrustSliceAccumulator> {
    let (mz_min, mz_max) = request.mz_bounds();
    let tof_min = metadata.mz_converter.invert(mz_min).floor().max(0.0);
    let tof_max = metadata.mz_converter.invert(mz_max).ceil().max(0.0);
    if !tof_min.is_finite() || !tof_max.is_finite() || tof_min > tof_max {
        anyhow::bail!("timsrust produced invalid TOF bounds for requested m/z slice");
    }
    let mut mz_bins = vec![
        MzProfileBin {
            mz_center: 0.0,
            summed_intensity: 0.0,
        };
        mz_profile_bins
    ];
    initialize_mz_histogram_bins(&mut mz_bins, mz_min, mz_max);
    Ok(TimsrustSliceAccumulator {
        mz_min,
        mz_max,
        tof_min: tof_min as u32,
        tof_max: tof_max as u32,
        mz_bins,
        frames_with_signal: 0,
        matched_events: 0,
        rt_profile: Vec::new(),
        im_accumulator: BTreeMap::new(),
    })
}

fn finalize_timsrust_slice_payload(
    accumulator: TimsrustSliceAccumulator,
    acquisition_mode: &str,
    frames_considered: usize,
) -> BrukerBridgePayload {
    let im_profile = accumulator
        .im_accumulator
        .into_iter()
        .map(|(_, (mobility, summed_intensity))| BrukerBridgeImRow {
            mobility,
            summed_intensity,
        })
        .collect::<Vec<_>>();
    let mz_profile = accumulator
        .mz_bins
        .into_iter()
        .map(|bin| BrukerBridgeMzRow {
            mz_center: bin.mz_center,
            summed_intensity: bin.summed_intensity,
        })
        .collect::<Vec<_>>();

    BrukerBridgePayload {
        acquisition_mode: acquisition_mode.to_string(),
        intensity_column: "timsrust_corrected_intensity_values".to_string(),
        mz_min: accumulator.mz_min,
        mz_max: accumulator.mz_max,
        frames_considered,
        frames_with_signal: accumulator.frames_with_signal,
        matched_events: accumulator.matched_events,
        rt_profile: accumulator.rt_profile,
        mz_profile,
        im_profile,
    }
}

fn build_timsrust_pseudo_ms2(
    path: &Path,
    request: &DiaSliceRequest,
    out_prefix: &Path,
    options: &PseudoMs2Options,
    precursor_rt_points: &[(f64, f64)],
    verbosity: Verbosity,
) -> anyhow::Result<PseudoMs2Report> {
    let peptide = request
        .peptide_target
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--pseudo-ms2 requires --peptide"))?
        .clone();
    if peptide.fragments.is_empty() {
        anyhow::bail!(
            "peptide {} produced no theoretical fragments",
            peptide.sequence
        );
    }

    let metadata = MetadataReader::new(path)
        .with_context(|| "failed to read Bruker metadata with timsrust")?;
    let frame_reader =
        FrameReader::new(path).with_context(|| "failed to open Bruker .d with timsrust")?;
    let acquisition = frame_reader.get_acquisition();
    let precursor_evidence = match (request.rt_min, request.rt_max) {
        (Some(lo), Some(hi)) => explicit_precursor_rt_evidence(precursor_rt_points, lo, hi),
        _ => infer_precursor_rt_window_from_points(
            precursor_rt_points,
            peptide.precursor_mz,
            options.rt_window_min,
        )?,
    };
    verbosity.status(format_args!(
        "Pseudo-MS2: using RT {} and aggregating {} theoretical fragments",
        precursor_evidence.rt_window.label(),
        peptide.fragments.len()
    ));

    let mut accumulators = Vec::<PseudoFragmentAccumulator>::new();
    for fragment in &peptide.fragments {
        let (mz_min, mz_max) = mz_window(fragment.mz, request);
        let (tof_min, tof_max) = tof_bounds(&metadata, mz_min, mz_max)
            .with_context(|| format!("failed to compute TOF bounds for {}", fragment.label))?;
        accumulators.push(PseudoFragmentAccumulator {
            tof_min,
            tof_max,
            evidence: PseudoMs2FragmentEvidence {
                label: fragment.label.clone(),
                series: fragment.series.clone(),
                cleavage_index: fragment.cleavage_index,
                neutral_loss: fragment.neutral_loss,
                charge: fragment.charge,
                mz: fragment.mz,
                summed_intensity: 0.0,
                matched_events: 0,
                frames_with_signal: 0,
                apex_rt: None,
                apex_intensity: 0.0,
            },
        });
    }

    let (quad_min, quad_max) = pseudo_ms2_quad_bounds(request, &peptide, acquisition);
    verbosity.detail(format_args!(
        "Pseudo-MS2: acquisition={:?} quadrupole filter {:?}-{:?}",
        acquisition, quad_min, quad_max
    ));
    let mut frames_considered = 0usize;
    let mut frames_with_signal = 0usize;
    let mut matched_events = 0usize;
    let total_frames = frame_reader.len();
    let candidate_frames = collect_timsrust_ms2_frame_candidates(&frame_reader, |rt_seconds| {
        let rt_minutes = rt_seconds / 60.0;
        rt_minutes >= precursor_evidence.rt_window.min
            && rt_minutes <= precursor_evidence.rt_window.max
    })?;
    verbosity.detail(format_args!(
        "Pseudo-MS2: scanning {} MS2 frames after RT filtering from {} total frames",
        candidate_frames.len(),
        total_frames
    ));
    let mut progress = ProgressReporter::new(
        "Pseudo-MS2 native MS2 frames",
        candidate_frames.len(),
        verbosity,
    );

    for candidate in candidate_frames {
        let rt_minutes = candidate.rt_seconds / 60.0;
        frames_considered += 1;
        let frame = frame_reader.get(candidate.reader_index).with_context(|| {
            format!("failed to read Bruker frame {}", candidate.reader_index + 1)
        })?;
        let mut frame_intensities = vec![0.0_f64; accumulators.len()];
        let mut frame_events = vec![0usize; accumulators.len()];

        for scan_index in 0..frame.scan_offsets.len().saturating_sub(1) {
            let mobility = metadata.im_converter.convert(scan_index as u32);
            if !mobility_in_window(mobility, request.im_min, request.im_max) {
                continue;
            }
            if !quad_in_window(&frame, scan_index, quad_min, quad_max) {
                continue;
            }

            let start = frame.scan_offsets[scan_index];
            let end = frame.scan_offsets[scan_index + 1];
            if start >= end {
                continue;
            }
            let scan_tofs = &frame.tof_indices[start..end];
            for frag_idx in 0..accumulators.len() {
                let peak_start =
                    scan_tofs.partition_point(|tof| *tof < accumulators[frag_idx].tof_min);
                let peak_end =
                    scan_tofs.partition_point(|tof| *tof <= accumulators[frag_idx].tof_max);
                for relative_peak_index in peak_start..peak_end {
                    let peak_index = start + relative_peak_index;
                    let intensity = frame.get_corrected_intensity(peak_index);
                    if !intensity.is_finite() || intensity <= 0.0 {
                        continue;
                    }
                    accumulators[frag_idx].evidence.summed_intensity += intensity;
                    accumulators[frag_idx].evidence.matched_events += 1;
                    frame_intensities[frag_idx] += intensity;
                    frame_events[frag_idx] += 1;
                }
            }
        }

        let mut frame_has_signal = false;
        for (frag_idx, frame_intensity) in frame_intensities.iter().copied().enumerate() {
            let events = frame_events[frag_idx];
            if events == 0 {
                continue;
            }
            frame_has_signal = true;
            matched_events += events;
            let evidence = &mut accumulators[frag_idx].evidence;
            evidence.frames_with_signal += 1;
            if frame_intensity > evidence.apex_intensity {
                evidence.apex_intensity = frame_intensity;
                evidence.apex_rt = Some(rt_minutes);
            }
        }
        if frame_has_signal {
            frames_with_signal += 1;
        }
        progress.advance();
    }
    progress.finish();

    Ok(PseudoMs2Report {
        input_path: path.to_path_buf(),
        out_prefix: out_prefix.to_path_buf(),
        peptide,
        rt_window: precursor_evidence.rt_window,
        frames_considered,
        frames_with_signal,
        matched_events,
        precursor_frames_with_signal: precursor_evidence.frames_with_signal,
        precursor_apex_rt: precursor_evidence.apex_rt,
        precursor_apex_intensity: precursor_evidence.apex_intensity,
        fragments: accumulators
            .into_iter()
            .map(|accumulator| accumulator.evidence)
            .collect(),
    })
}

#[derive(Clone, Debug)]
struct PseudoFragmentAccumulator {
    tof_min: u32,
    tof_max: u32,
    evidence: PseudoMs2FragmentEvidence,
}

#[derive(Clone, Debug)]
struct PrecursorRtEvidence {
    rt_window: PseudoMs2RtWindow,
    frames_with_signal: usize,
    apex_rt: Option<f64>,
    apex_intensity: f64,
}

fn explicit_precursor_rt_evidence(
    points: &[(f64, f64)],
    rt_min: f64,
    rt_max: f64,
) -> PrecursorRtEvidence {
    let mut filtered = points
        .iter()
        .copied()
        .filter(|(rt, intensity)| {
            rt.is_finite()
                && *rt >= rt_min
                && *rt <= rt_max
                && intensity.is_finite()
                && *intensity > 0.0
        })
        .collect::<Vec<_>>();
    filtered.sort_by(|left, right| {
        left.0
            .partial_cmp(&right.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let (apex_rt, apex_intensity) = filtered
        .iter()
        .copied()
        .max_by(|left, right| {
            left.1
                .partial_cmp(&right.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(rt, intensity)| (Some(rt), intensity))
        .unwrap_or((None, 0.0));
    PrecursorRtEvidence {
        rt_window: PseudoMs2RtWindow {
            min: rt_min,
            max: rt_max,
            source: "explicit",
        },
        frames_with_signal: filtered.len(),
        apex_rt,
        apex_intensity,
    }
}

fn infer_precursor_rt_window_from_points(
    points: &[(f64, f64)],
    precursor_mz: f64,
    window_min: f64,
) -> anyhow::Result<PrecursorRtEvidence> {
    let mut points = points
        .iter()
        .copied()
        .filter(|(rt, intensity)| rt.is_finite() && intensity.is_finite() && *intensity > 0.0)
        .collect::<Vec<_>>();
    if points.is_empty() {
        anyhow::bail!(
            "could not infer pseudo-MS2 RT window: no extracted precursor signal for m/z {:.4}; pass --rt-min/--rt-max explicitly",
            precursor_mz
        );
    }
    points.sort_by(|left, right| {
        left.0
            .partial_cmp(&right.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let (rt_window, apex_rt, apex_intensity) = best_precursor_rt_window(&points, window_min);
    Ok(PrecursorRtEvidence {
        rt_window,
        frames_with_signal: points.len(),
        apex_rt: Some(apex_rt),
        apex_intensity,
    })
}

fn best_precursor_rt_window(
    points: &[(f64, f64)],
    window_min: f64,
) -> (PseudoMs2RtWindow, f64, f64) {
    let mut best_start = points[0].0;
    let mut best_sum = f64::NEG_INFINITY;
    let mut best_left = 0usize;
    let mut best_right = 0usize;
    let mut right = 0usize;
    let mut running = 0.0;

    for left in 0..points.len() {
        let start = points[left].0;
        while right < points.len() && points[right].0 <= start + window_min {
            running += points[right].1;
            right += 1;
        }
        if running > best_sum {
            best_sum = running;
            best_start = start;
            best_left = left;
            best_right = right;
        }
        running -= points[left].1;
    }

    let (apex_rt, apex_intensity) = points[best_left..best_right]
        .iter()
        .copied()
        .max_by(|left, right| {
            left.1
                .partial_cmp(&right.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(points[best_left]);
    (
        PseudoMs2RtWindow {
            min: best_start,
            max: best_start + window_min,
            source: "precursor-inferred",
        },
        apex_rt,
        apex_intensity,
    )
}

fn mz_window(center: f64, request: &DiaSliceRequest) -> (f64, f64) {
    let delta = request
        .mz_da
        .unwrap_or_else(|| center * request.mz_ppm / 1_000_000.0);
    (center - delta, center + delta)
}

fn tof_bounds(metadata: &Metadata, mz_min: f64, mz_max: f64) -> anyhow::Result<(u32, u32)> {
    let tof_min = metadata.mz_converter.invert(mz_min).floor().max(0.0);
    let tof_max = metadata.mz_converter.invert(mz_max).ceil().max(0.0);
    if !tof_min.is_finite() || !tof_max.is_finite() || tof_min > tof_max {
        anyhow::bail!(
            "timsrust produced invalid TOF bounds for m/z {:.4}-{:.4}",
            mz_min,
            mz_max
        );
    }
    Ok((tof_min as u32, tof_max as u32))
}

fn pseudo_ms2_quad_bounds(
    request: &DiaSliceRequest,
    peptide: &DiaPeptideTarget,
    acquisition: AcquisitionType,
) -> (Option<f64>, Option<f64>) {
    match (request.quad_min, request.quad_max) {
        (Some(lo), Some(hi)) => (Some(lo), Some(hi)),
        _ if acquisition == AcquisitionType::DIAPASEF => {
            (Some(peptide.precursor_mz), Some(peptide.precursor_mz))
        }
        _ => (None, None),
    }
}

fn rt_seconds_in_window(rt_seconds: f64, rt_min: Option<f64>, rt_max: Option<f64>) -> bool {
    match (rt_min, rt_max) {
        (Some(lo), Some(hi)) => rt_seconds >= lo * 60.0 && rt_seconds <= hi * 60.0,
        _ => true,
    }
}

fn mobility_in_window(mobility: f64, im_min: Option<f64>, im_max: Option<f64>) -> bool {
    match (im_min, im_max) {
        (Some(lo), Some(hi)) => mobility >= lo && mobility <= hi,
        _ => true,
    }
}

fn quad_in_window(
    frame: &Frame,
    scan_index: usize,
    quad_min: Option<f64>,
    quad_max: Option<f64>,
) -> bool {
    let (Some(lo), Some(hi)) = (quad_min, quad_max) else {
        return true;
    };
    let settings = frame.quadrupole_settings.as_ref();
    if settings.len() == 0 {
        return false;
    }
    settings
        .scan_starts
        .iter()
        .zip(settings.scan_ends.iter())
        .zip(
            settings
                .isolation_mz
                .iter()
                .zip(settings.isolation_width.iter()),
        )
        .any(
            |((scan_start, scan_end), (isolation_mz, isolation_width))| {
                if scan_index < *scan_start || scan_index >= *scan_end {
                    return false;
                }
                let half_width = *isolation_width / 2.0;
                let window_min = *isolation_mz - half_width;
                let window_max = *isolation_mz + half_width;
                window_min <= hi && window_max >= lo
            },
        )
}

fn timsrust_acquisition_label(acquisition: AcquisitionType) -> &'static str {
    match acquisition {
        AcquisitionType::DDAPASEF => "ddaPASEF",
        AcquisitionType::DIAPASEF => "diaPASEF",
        AcquisitionType::DiagonalDIAPASEF => "diagonal diaPASEF",
        AcquisitionType::Unknown => "unknown",
    }
}

fn resolve_python(python_override: Option<&Path>) -> anyhow::Result<PathBuf> {
    if let Some(path) = python_override {
        return Ok(path.to_path_buf());
    }
    if let Some(path) = env::var_os(BRUKER_PYTHON_ENV_VAR).map(PathBuf::from) {
        return Ok(path);
    }

    for candidate in ["python", "python3"] {
        if Command::new(candidate).arg("--version").output().is_ok() {
            return Ok(PathBuf::from(candidate));
        }
    }

    anyhow::bail!(
        "Bruker .d support via alphaTims requires a Python interpreter. Pass --python <exe> or set {}",
        BRUKER_PYTHON_ENV_VAR
    );
}

fn run_bruker_bridge(
    path: &Path,
    request: &DiaSliceRequest,
    runtime_path: &Path,
    python: &Path,
    mz_profile_bins: usize,
) -> anyhow::Result<BrukerBridgePayload> {
    let script_path = write_temp_bridge_script()?;
    let mut command = Command::new(python);
    command
        .arg(script_path.as_path())
        .arg("--bruker")
        .arg(path)
        .arg("--bruker-so")
        .arg(runtime_path)
        .arg("--mz")
        .arg(request.mz.to_string())
        .arg("--mz-ppm")
        .arg(request.mz_ppm.to_string())
        .arg("--mz-bins")
        .arg(mz_profile_bins.to_string())
        .env("PYTHONUNBUFFERED", "1");

    if let Some(value) = request.mz_da {
        command.arg("--mz-da").arg(value.to_string());
    }
    if let (Some(lo), Some(hi)) = (request.rt_min, request.rt_max) {
        command
            .arg("--rt-min")
            .arg(lo.to_string())
            .arg("--rt-max")
            .arg(hi.to_string());
    }
    if let (Some(lo), Some(hi)) = (request.im_min, request.im_max) {
        command
            .arg("--im-min")
            .arg(lo.to_string())
            .arg("--im-max")
            .arg(hi.to_string());
    }
    if let (Some(lo), Some(hi)) = (request.quad_min, request.quad_max) {
        command
            .arg("--quad-min")
            .arg(lo.to_string())
            .arg("--quad-max")
            .arg(hi.to_string());
    }

    let output = command
        .output()
        .with_context(|| format!("failed to launch Python bridge `{}`", python.display()))?;
    let _ = fs::remove_file(&script_path);

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        anyhow::bail!(
            "Bruker alphaTims bridge failed.\nstdout:\n{}\n\nstderr:\n{}",
            stdout.trim(),
            stderr.trim()
        );
    }

    serde_json::from_slice::<BrukerBridgePayload>(&output.stdout)
        .context("failed to parse JSON from Bruker alphaTims bridge")
}

fn resolve_bruker_runtime(bruker_so_override: Option<&Path>) -> anyhow::Result<BrukerRuntime> {
    if let Some(path) = bruker_so_override {
        return load_bruker_runtime_candidate(path).with_context(|| {
            format!(
                "Bruker .d support requires timsdata.so. The explicit --bruker-so path did not work. Try a valid library path or place timsdata.so at {}",
                DEFAULT_BRUKER_SO_PATHS[0]
            )
        });
    }

    if let Some(path) = env::var_os(BRUKER_SO_ENV_VAR).map(PathBuf::from) {
        return load_bruker_runtime_candidate(path.as_path()).with_context(|| {
            format!(
                "Bruker .d support requires timsdata.so. The {BRUKER_SO_ENV_VAR} override did not work. Try a valid library path or place timsdata.so at {}",
                DEFAULT_BRUKER_SO_PATHS[0]
            )
        });
    }

    let mut checked = Vec::<PathBuf>::new();
    let mut failures = Vec::<String>::new();
    for candidate in DEFAULT_BRUKER_SO_PATHS {
        let path = PathBuf::from(candidate);
        checked.push(path.clone());
        if !path.exists() {
            continue;
        }
        match load_bruker_runtime_candidate(path.as_path()) {
            Ok(runtime) => return Ok(runtime),
            Err(err) => failures.push(format!("{} ({err:#})", path.display())),
        }
    }

    let checked_list = checked
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    if failures.is_empty() {
        anyhow::bail!(
            "Bruker .d support requires timsdata.so. Checked: {}. Pass --bruker-so <path>, set {}, or place timsdata.so at {}",
            checked_list,
            BRUKER_SO_ENV_VAR,
            DEFAULT_BRUKER_SO_PATHS[0]
        );
    }
    anyhow::bail!(
        "Bruker .d support found candidate timsdata.so files, but none could be loaded. Checked: {}. Load failures: {}",
        checked_list,
        failures.join("; ")
    );
}

fn load_bruker_runtime_candidate(path: &Path) -> anyhow::Result<BrukerRuntime> {
    if !path.exists() {
        anyhow::bail!("missing {}", path.display());
    }
    if !path.is_file() {
        anyhow::bail!("{} is not a file", path.display());
    }
    unsafe { Library::new(path) }.with_context(|| format!("failed to load {}", path.display()))?;
    Ok(BrukerRuntime {
        path: path.to_path_buf(),
    })
}

fn initialize_mz_bins(bins: &mut [MzProfileBin], mz_min: f64, mz_max: f64) {
    if bins.is_empty() {
        return;
    }
    if bins.len() == 1 {
        bins[0].mz_center = (mz_min + mz_max) / 2.0;
        return;
    }
    let span = mz_max - mz_min;
    let len = bins.len();
    for (idx, bin) in bins.iter_mut().enumerate() {
        let frac = idx as f64 / (len - 1) as f64;
        bin.mz_center = mz_min + span * frac;
    }
}

fn accumulate_mz_bin(bins: &mut [MzProfileBin], mz_min: f64, mz_max: f64, mz: f64, intensity: f64) {
    if bins.is_empty() {
        return;
    }
    if bins.len() == 1 || (mz_max - mz_min).abs() <= f64::EPSILON {
        bins[0].summed_intensity += intensity;
        return;
    }
    let frac = ((mz - mz_min) / (mz_max - mz_min)).clamp(0.0, 1.0);
    let idx = ((bins.len() - 1) as f64 * frac).round() as usize;
    bins[idx.min(bins.len() - 1)].summed_intensity += intensity;
}

fn initialize_mz_histogram_bins(bins: &mut [MzProfileBin], mz_min: f64, mz_max: f64) {
    if bins.is_empty() {
        return;
    }
    if bins.len() == 1 {
        bins[0].mz_center = (mz_min + mz_max) / 2.0;
        return;
    }
    let width = (mz_max - mz_min) / bins.len() as f64;
    for (idx, bin) in bins.iter_mut().enumerate() {
        bin.mz_center = mz_min + (idx as f64 + 0.5) * width;
    }
}

fn accumulate_mz_histogram_bin(
    bins: &mut [MzProfileBin],
    mz_min: f64,
    mz_max: f64,
    mz: f64,
    intensity: f64,
) {
    if bins.is_empty() {
        return;
    }
    if bins.len() == 1 || (mz_max - mz_min).abs() <= f64::EPSILON {
        bins[0].summed_intensity += intensity;
        return;
    }
    let width = (mz_max - mz_min) / bins.len() as f64;
    if width <= 0.0 || !width.is_finite() {
        return;
    }
    let idx = ((mz - mz_min) / width).floor() as isize;
    let idx = idx.clamp(0, bins.len() as isize - 1) as usize;
    bins[idx].summed_intensity += intensity;
}

fn rt_in_window(rt_minutes: Option<f64>, rt_min: Option<f64>, rt_max: Option<f64>) -> bool {
    match (rt_minutes, rt_min, rt_max) {
        (_, None, None) => true,
        (Some(rt), Some(lo), Some(hi)) => rt >= lo && rt <= hi,
        _ => false,
    }
}

fn write_outputs(
    summary: &DiaSliceSummary,
    request: &DiaSliceRequest,
    rt_rows: &[RtProfileRow],
    mz_bins: &[MzProfileBin],
    im_rows: Option<&[MobilityProfileRow]>,
    rt_smooth: bool,
    rt_smooth_window: usize,
) -> anyhow::Result<()> {
    if let Some(parent) = summary.out_prefix.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }

    write_summary_text(summary, request)?;
    write_rt_profile(summary, rt_rows)?;
    write_mz_profile(summary, mz_bins)?;
    write_im_profile(summary, im_rows)?;
    write_summary_svg(
        summary,
        request,
        rt_rows,
        mz_bins,
        im_rows,
        rt_smooth,
        rt_smooth_window,
    )?;
    Ok(())
}

fn write_bruker_run_tic_json(
    out_prefix: &Path,
    run_tic: &BrukerRunTicSummary,
) -> anyhow::Result<()> {
    let path = output_path(out_prefix, "run_tic.json");
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }
    let file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    serde_json::to_writer_pretty(file, run_tic)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn write_pseudo_ms2_outputs(report: &PseudoMs2Report) -> anyhow::Result<()> {
    if let Some(parent) = report.out_prefix.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }
    write_pseudo_ms2_tsv(report)?;
    write_pseudo_ms2_svg(report)?;
    Ok(())
}

fn write_pseudo_ms2_tsv(report: &PseudoMs2Report) -> anyhow::Result<()> {
    let path = output_path(&report.out_prefix, "pseudo_ms2.tsv");
    let mut file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    writeln!(file, "# source\t{}", report.input_path.display())?;
    writeln!(
        file,
        "# peptide_input\t{}",
        sanitize_tsv(&report.peptide.input)
    )?;
    writeln!(
        file,
        "# peptide_modified_sequence\t{}",
        sanitize_tsv(&report.peptide.modified_sequence)
    )?;
    writeln!(file, "# peptide_charge\t{}", report.peptide.charge)?;
    writeln!(file, "# precursor_mz\t{:.6}", report.peptide.precursor_mz)?;
    writeln!(file, "# rt_window\t{}", report.rt_window.label())?;
    writeln!(
        file,
        "# precursor_rt_apex\t{}\t{:.6}",
        format_optional_f64(report.precursor_apex_rt, 6),
        report.precursor_apex_intensity
    )?;
    writeln!(
        file,
        "# precursor_frames_with_signal\t{}",
        report.precursor_frames_with_signal
    )?;
    writeln!(
        file,
        "# frames_considered\t{}\tframes_with_signal\t{}\tmatched_events\t{}",
        report.frames_considered, report.frames_with_signal, report.matched_events
    )?;
    writeln!(
        file,
        "fragment\tseries\tcharge\tmz\tsummed_intensity\tmatched_events\tframes_with_signal\tapex_rt\tapex_intensity"
    )?;
    for row in &report.fragments {
        writeln!(
            file,
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{:.6}",
            row.label,
            row.series,
            row.charge,
            row.mz,
            row.summed_intensity,
            row.matched_events,
            row.frames_with_signal,
            format_optional_f64(row.apex_rt, 6),
            row.apex_intensity,
        )?;
    }
    Ok(())
}

fn write_pseudo_ms2_svg(report: &PseudoMs2Report) -> anyhow::Result<()> {
    let path = output_path(&report.out_prefix, "pseudo_ms2.svg");
    let mut svg = String::new();
    let table_rows = build_pseudo_ion_table_rows(report);
    let neutral_loss_cutoff = pseudo_ms2_neutral_loss_display_cutoff(&report.fragments);
    let table_width = pseudo_ion_table_width(&table_rows);
    let table_height = pseudo_ion_table_height(table_rows.len());
    let plot_left = 96.0;
    let plot_top = 210.0;
    let plot_width = 1094.0;
    let plot_height = 410.0;
    let table_left = SVG_MARGIN_X;
    let table_top = plot_top + plot_height + 104.0;
    let width = SVG_WIDTH;
    let height = (table_top + table_height + 72.0).max(820.0).ceil() as u32;
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">",
        w = width,
        h = height
    );
    let _ = writeln!(
        svg,
        "<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"{}\"/>",
        COLOR_BG
    );

    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        44.0,
        SVG_TITLE_FONT,
        30.0,
        COLOR_TEXT,
        Some("700"),
        &["Pseudo-MS2 fragment evidence".to_string()],
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        70.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &format!("Source: {}", compact_display_path(&report.input_path, 78)),
            SVG_META_FONT,
            width as f64 - 2.0 * SVG_MARGIN_X,
            1,
        ),
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        96.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &format!(
                "Target: {}/{} precursor m/z {:.4} | RT {} | {} frames, {} with signal",
                report.peptide.modified_sequence,
                report.peptide.charge,
                report.peptide.precursor_mz,
                report.rt_window.label(),
                report.frames_considered,
                report.frames_with_signal
            ),
            SVG_META_FONT,
            width as f64 - 2.0 * SVG_MARGIN_X,
            2,
        ),
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        136.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &format!(
                "Precursor RT evidence: {} frames with signal, apex RT {} min, apex intensity {:.2e}",
                report.precursor_frames_with_signal,
                format_optional_f64(report.precursor_apex_rt, 3),
                report.precursor_apex_intensity,
            ),
            SVG_META_FONT,
            width as f64 - 2.0 * SVG_MARGIN_X,
            1,
        ),
    );

    let x_domain = padded_range(report.fragments.iter().map(|row| row.mz));
    let y_domain = padded_y_range(report.fragments.iter().map(|row| row.summed_intensity));
    let canvas = SvgCanvas::new(
        plot_left,
        plot_top,
        plot_width,
        plot_height,
        x_domain,
        y_domain,
    );
    draw_pseudo_ms2_axes(&mut svg, canvas, x_domain, y_domain);
    draw_pseudo_ms2_sticks(&mut svg, canvas, &report.fragments, neutral_loss_cutoff);
    draw_pseudo_ms2_legend(&mut svg, canvas.right() - 210.0, canvas.top() - 26.0);
    draw_pseudo_ms2_ion_table(
        &mut svg,
        table_left,
        table_top,
        table_width,
        table_height,
        &table_rows,
        &report.peptide.sequence,
        neutral_loss_cutoff,
    );
    svg.push_str("</svg>\n");
    fs::write(&path, svg).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn build_pseudo_ion_table_rows(report: &PseudoMs2Report) -> Vec<PseudoIonTableRow> {
    let mut by_key = BTreeMap::<(&str, usize, u8), Vec<&PseudoMs2FragmentEvidence>>::new();
    for fragment in &report.fragments {
        by_key
            .entry((
                fragment.series.as_str(),
                fragment.cleavage_index,
                fragment.charge,
            ))
            .or_default()
            .push(fragment);
    }

    let residue_count = report.peptide.sequence.chars().count();
    let mut rows = Vec::with_capacity(residue_count.saturating_sub(1));
    for cleavage_index in 1..residue_count {
        rows.push(PseudoIonTableRow {
            cleavage_index,
            y_ordinal: residue_count - cleavage_index,
            b1: build_pseudo_ion_table_cell_entries(&by_key, "b", cleavage_index, 1),
            b2: build_pseudo_ion_table_cell_entries(&by_key, "b", cleavage_index, 2),
            y1: build_pseudo_ion_table_cell_entries(&by_key, "y", cleavage_index, 1),
            y2: build_pseudo_ion_table_cell_entries(&by_key, "y", cleavage_index, 2),
        });
    }
    rows
}

fn build_pseudo_ion_table_cell_entries(
    by_key: &BTreeMap<(&str, usize, u8), Vec<&PseudoMs2FragmentEvidence>>,
    series: &'static str,
    cleavage_index: usize,
    charge: u8,
) -> Vec<PseudoIonTableCell> {
    let mut fragments = by_key
        .get(&(series, cleavage_index, charge))
        .cloned()
        .unwrap_or_default();
    fragments.sort_by(|left, right| {
        neutral_loss_rank(left.neutral_loss)
            .cmp(&neutral_loss_rank(right.neutral_loss))
            .then_with(|| {
                right
                    .summed_intensity
                    .partial_cmp(&left.summed_intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                left.mz
                    .partial_cmp(&right.mz)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    fragments
        .into_iter()
        .map(pseudo_ion_table_cell)
        .collect::<Vec<_>>()
}

fn pseudo_ion_table_cell(fragment: &PseudoMs2FragmentEvidence) -> PseudoIonTableCell {
    PseudoIonTableCell {
        label: fragment.label.clone(),
        series: fragment.series.clone(),
        neutral_loss: fragment.neutral_loss,
        mz: fragment.mz,
        summed_intensity: fragment.summed_intensity,
        matched_events: fragment.matched_events,
        frames_with_signal: fragment.frames_with_signal,
        apex_rt: fragment.apex_rt,
    }
}

fn pseudo_ms2_neutral_loss_display_cutoff(fragments: &[PseudoMs2FragmentEvidence]) -> f64 {
    let max_base_intensity = fragments
        .iter()
        .filter(|fragment| fragment.neutral_loss.is_none())
        .map(|fragment| fragment.summed_intensity)
        .filter(|intensity| intensity.is_finite() && *intensity > 0.0)
        .fold(0.0, f64::max);
    let max_intensity = if max_base_intensity > 0.0 {
        max_base_intensity
    } else {
        fragments
            .iter()
            .map(|fragment| fragment.summed_intensity)
            .filter(|intensity| intensity.is_finite() && *intensity > 0.0)
            .fold(0.0, f64::max)
    };
    max_intensity * PSEUDO_MS2_NEUTRAL_LOSS_MIN_RELATIVE_INTENSITY
}

fn pseudo_ms2_neutral_loss_signal_is_visible(
    summed_intensity: f64,
    neutral_loss: Option<NeutralLossKind>,
    neutral_loss_cutoff: f64,
) -> bool {
    neutral_loss.is_none()
        || (summed_intensity.is_finite()
            && summed_intensity > 0.0
            && summed_intensity >= neutral_loss_cutoff)
}

fn pseudo_ion_table_has_charge2(rows: &[PseudoIonTableRow]) -> bool {
    rows.iter()
        .any(|row| !row.b2.is_empty() || !row.y2.is_empty())
}

fn pseudo_ion_table_width(_rows: &[PseudoIonTableRow]) -> f64 {
    SVG_WIDTH as f64 - 2.0 * SVG_MARGIN_X
}

fn pseudo_ion_table_height(row_count: usize) -> f64 {
    112.0 + row_count as f64 * pseudo_ion_table_row_height(row_count) + 68.0
}

fn pseudo_ion_table_row_height(row_count: usize) -> f64 {
    if row_count <= 18 {
        31.0
    } else if row_count <= 32 {
        28.0
    } else {
        24.0
    }
}

fn pseudo_ion_table_font_size(row_count: usize) -> f64 {
    if row_count <= 18 {
        16.0
    } else if row_count <= 32 {
        14.0
    } else {
        12.0
    }
}

fn draw_pseudo_ms2_ion_table(
    svg: &mut String,
    left: f64,
    top: f64,
    width: f64,
    height: f64,
    rows: &[PseudoIonTableRow],
    sequence: &str,
    neutral_loss_cutoff: f64,
) {
    let row_height = pseudo_ion_table_row_height(rows.len());
    let font_size = pseudo_ion_table_font_size(rows.len());
    let show_charge2 = pseudo_ion_table_has_charge2(rows);
    let pad = 18.0;
    let title_y = top + 24.0;
    let meta_y = top + 47.0;
    let meta2_y = top + 66.0;
    let header_y = top + 96.0;
    let row_start_y = top + 125.0;
    let table_bottom = top + height;

    let _ = writeln!(
        svg,
        "<rect x=\"{left:.2}\" y=\"{top:.2}\" width=\"{width:.2}\" height=\"{height:.2}\" rx=\"12\" fill=\"#fbfdff\" stroke=\"{}\" stroke-width=\"1\"/>",
        COLOR_CARD_BORDER
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"18\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">Pseudo-MS2 ion table</text>",
        left + pad,
        title_y,
        COLOR_TEXT
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"13\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">color = DIA evidence; grey = theoretical base ion</text>",
        left + pad,
        meta_y,
        COLOR_SUBTLE
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"13\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">suffixes: +p = -H3PO4 phosphoric acid loss; +w = -H2O water loss; +n = -NH3 ammonia loss</text>",
        left + pad,
        meta2_y,
        COLOR_SUBTLE
    );

    draw_pseudo_ion_table_block(
        svg,
        left + pad,
        width - pad * 2.0,
        header_y,
        row_start_y,
        row_height,
        font_size,
        show_charge2,
        rows,
        neutral_loss_cutoff,
    );

    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"13\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">sequence: {}</text>",
        left + pad,
        table_bottom - 40.0,
        COLOR_SUBTLE,
        escape_xml(sequence)
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"13\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">neutral losses with at least {:.0}% of max fragment signal</text>",
        left + pad,
        table_bottom - 16.0,
        COLOR_SUBTLE,
        PSEUDO_MS2_NEUTRAL_LOSS_MIN_RELATIVE_INTENSITY * 100.0
    );
}

fn draw_pseudo_ion_table_block(
    svg: &mut String,
    left: f64,
    width: f64,
    header_y: f64,
    row_start_y: f64,
    row_height: f64,
    font_size: f64,
    show_charge2: bool,
    rows: &[PseudoIonTableRow],
    neutral_loss_cutoff: f64,
) {
    let cut_width = 58.0;
    let column_gap = 14.0;
    let value_columns = if show_charge2 { 4.0 } else { 2.0 };
    let value_width = (width - cut_width - column_gap * value_columns) / value_columns;

    let mut cursor = left;
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
        let _ = writeln!(
            svg,
            "<text x=\"{x:.2}\" y=\"{header_y:.2}\" font-size=\"12\" fill=\"{}\" text-anchor=\"{anchor}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{label}</text>",
            COLOR_AXIS
        );
    }
    let _ = writeln!(
        svg,
        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>",
        left,
        header_y + 8.0,
        left + width,
        header_y + 8.0,
        COLOR_GRID
    );

    for (idx, row) in rows.iter().enumerate() {
        let y = row_start_y + idx as f64 * row_height;
        if idx > 0 {
            let _ = writeln!(
                svg,
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#eef2f6\" stroke-width=\"1\"/>",
                left,
                y - row_height / 2.0 + 3.0,
                left + width,
                y - row_height / 2.0 + 3.0
            );
        }
        if let Some(x) = b2_right {
            write_pseudo_ion_cell(
                svg,
                x,
                y,
                "end",
                font_size,
                value_width,
                &row.b2,
                neutral_loss_cutoff,
            );
        }
        write_pseudo_ion_cell(
            svg,
            b1_right,
            y,
            "end",
            font_size,
            value_width,
            &row.b1,
            neutral_loss_cutoff,
        );
        let _ = writeln!(
            svg,
            "<text x=\"{cut_center:.2}\" y=\"{y:.2}\" font-size=\"{font_size:.1}\" fill=\"{}\" text-anchor=\"middle\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\"><title>{}</title>{}|{}</text>",
            COLOR_SUBTLE,
            escape_xml(&format!(
                "b{} / y{} cleavage",
                row.cleavage_index, row.y_ordinal
            )),
            row.cleavage_index,
            row.y_ordinal
        );
        write_pseudo_ion_cell(
            svg,
            y1_left,
            y,
            "start",
            font_size,
            value_width,
            &row.y1,
            neutral_loss_cutoff,
        );
        if let Some(x) = y2_left {
            write_pseudo_ion_cell(
                svg,
                x,
                y,
                "start",
                font_size,
                value_width,
                &row.y2,
                neutral_loss_cutoff,
            );
        }
    }
}

fn write_pseudo_ion_cell(
    svg: &mut String,
    x: f64,
    y: f64,
    anchor: &str,
    font_size: f64,
    max_width: f64,
    entries: &[PseudoIonTableCell],
    neutral_loss_cutoff: f64,
) {
    if entries.is_empty() {
        return;
    }
    let color = entries
        .iter()
        .find(|entry| {
            entry.detected()
                && pseudo_ms2_neutral_loss_signal_is_visible(
                    entry.summed_intensity,
                    entry.neutral_loss,
                    neutral_loss_cutoff,
                )
        })
        .or_else(|| entries.first())
        .map(|entry| {
            if entry.detected()
                && pseudo_ms2_neutral_loss_signal_is_visible(
                    entry.summed_intensity,
                    entry.neutral_loss,
                    neutral_loss_cutoff,
                )
            {
                pseudo_ms2_series_color(&entry.series)
            } else {
                "#aeb8c4"
            }
        })
        .unwrap_or("#aeb8c4");
    let title = entries
        .iter()
        .map(|entry| {
            format!(
                "{} theoretical {:.4} | summed {:.3e} | events {} | frames {} | apex RT {}",
                entry.label,
                entry.mz,
                entry.summed_intensity,
                entry.matched_events,
                entry.frames_with_signal,
                format_optional_f64(entry.apex_rt, 3)
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    let budget = mono_char_budget(font_size, max_width);
    let text = pseudo_ion_table_visible_text(entries, budget, neutral_loss_cutoff);
    let _ = writeln!(
        svg,
        "<text x=\"{x:.2}\" y=\"{y:.2}\" font-size=\"{font_size:.1}\" fill=\"{color}\" text-anchor=\"{anchor}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\"><title>{}</title>{}</text>",
        escape_xml(&title),
        escape_xml(&text)
    );
}

fn pseudo_ion_table_visible_text(
    entries: &[PseudoIonTableCell],
    budget: usize,
    neutral_loss_cutoff: f64,
) -> String {
    if budget == 0 {
        return String::new();
    }
    let mut parts = Vec::new();
    if let Some(base) = entries.iter().find(|entry| entry.neutral_loss.is_none()) {
        parts.push(format!("{:.2}", base.mz));
        parts.extend(
            entries
                .iter()
                .filter(|entry| {
                    entry.neutral_loss.is_some()
                        && entry.detected()
                        && pseudo_ms2_neutral_loss_signal_is_visible(
                            entry.summed_intensity,
                            entry.neutral_loss,
                            neutral_loss_cutoff,
                        )
                })
                .filter_map(|entry| {
                    entry
                        .neutral_loss
                        .map(|loss| format!("+{}", neutral_loss_short_label(loss)))
                }),
        );
    } else {
        parts.extend(
            entries
                .iter()
                .filter(|entry| {
                    entry.neutral_loss.is_none()
                        || pseudo_ms2_neutral_loss_signal_is_visible(
                            entry.summed_intensity,
                            entry.neutral_loss,
                            neutral_loss_cutoff,
                        )
                })
                .map(pseudo_ion_table_entry_text),
        );
    }
    let text = parts.join(" ");
    if entries.iter().any(|entry| {
        entry.neutral_loss.is_some()
            && entry.detected()
            && !pseudo_ms2_neutral_loss_signal_is_visible(
                entry.summed_intensity,
                entry.neutral_loss,
                neutral_loss_cutoff,
            )
    }) {
        truncate_end(&format!("{text} ..."), budget)
    } else {
        truncate_end(&text, budget)
    }
}

fn pseudo_ion_table_entry_text(entry: &PseudoIonTableCell) -> String {
    match entry.neutral_loss {
        Some(loss) => format!("{}{:.2}", neutral_loss_short_label(loss), entry.mz),
        None => format!("{:.2}", entry.mz),
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

fn neutral_loss_short_label(loss: NeutralLossKind) -> &'static str {
    match loss {
        NeutralLossKind::Water => "w",
        NeutralLossKind::Ammonia => "n",
        NeutralLossKind::PhosphoricAcid => "p",
    }
}

fn draw_pseudo_ms2_axes(
    svg: &mut String,
    canvas: SvgCanvas,
    x_domain: CoordinateRange<f64>,
    y_domain: CoordinateRange<f64>,
) {
    let _ = writeln!(
        svg,
        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"12\" fill=\"#fbfdff\" stroke=\"{}\" stroke-width=\"1\"/>",
        canvas.left() - 24.0,
        canvas.top() - 46.0,
        canvas.width() + 48.0,
        canvas.height() + 120.0,
        COLOR_CARD_BORDER
    );
    let x_ticks = linear_ticks(x_domain.min(), x_domain.max(), 6);
    let y_ticks = linear_ticks(y_domain.min(), y_domain.max(), 5);
    for tick in x_ticks {
        let px = canvas.x(tick);
        let _ = writeln!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>",
            px,
            canvas.top(),
            px,
            canvas.bottom(),
            COLOR_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{:.1}</text>",
            px,
            canvas.bottom() + 28.0,
            SVG_TICK_FONT,
            COLOR_SUBTLE,
            tick
        );
    }
    for tick in y_ticks {
        let py = canvas.y(tick);
        let _ = writeln!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>",
            canvas.left(),
            py,
            canvas.right(),
            py,
            COLOR_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"end\" dominant-baseline=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{:.2e}</text>",
            canvas.left() - 12.0,
            py,
            SVG_TICK_FONT,
            COLOR_SUBTLE,
            tick
        );
    }
    let _ = writeln!(
        svg,
        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.2\"/>",
        canvas.left(),
        canvas.top(),
        canvas.width(),
        canvas.height(),
        COLOR_AXIS
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">fragment m/z</text>",
        canvas.left() + canvas.width() / 2.0,
        canvas.bottom() + 58.0,
        SVG_AXIS_LABEL_FONT,
        COLOR_TEXT
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\" transform=\"rotate(-90 {:.2} {:.2})\">summed extracted intensity</text>",
        canvas.left() - 66.0,
        canvas.top() + canvas.height() / 2.0,
        SVG_AXIS_LABEL_FONT,
        COLOR_SUBTLE,
        canvas.left() - 66.0,
        canvas.top() + canvas.height() / 2.0
    );
}

fn draw_pseudo_ms2_sticks(
    svg: &mut String,
    canvas: SvgCanvas,
    fragments: &[PseudoMs2FragmentEvidence],
    neutral_loss_cutoff: f64,
) {
    let baseline = canvas.y(0.0);

    for (idx, row) in fragments.iter().enumerate() {
        if !pseudo_ms2_neutral_loss_signal_is_visible(
            row.summed_intensity,
            row.neutral_loss,
            neutral_loss_cutoff,
        ) {
            continue;
        }
        let x = canvas.x(row.mz);
        if row.summed_intensity <= 0.0 {
            continue;
        }
        let y = canvas.y(row.summed_intensity);
        let color = pseudo_ms2_series_color(&row.series);
        let _ = writeln!(
            svg,
            "<line x1=\"{x:.2}\" y1=\"{baseline:.2}\" x2=\"{x:.2}\" y2=\"{y:.2}\" stroke=\"{color}\" stroke-width=\"2.4\" stroke-linecap=\"round\"/>"
        );
        let _ = writeln!(
            svg,
            "<circle cx=\"{x:.2}\" cy=\"{y:.2}\" r=\"3.4\" fill=\"{color}\"><title>{}: {:.4} m/z, {:.3e}</title></circle>",
            escape_xml(&row.label),
            row.mz,
            row.summed_intensity
        );
        let label_y = (y - 10.0 - (idx % 4) as f64 * 13.0).max(canvas.top() + 14.0);
        let _ = writeln!(
            svg,
            "<text x=\"{x:.2}\" y=\"{label_y:.2}\" font-size=\"12\" text-anchor=\"middle\" fill=\"{color}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\"><title>{}: {:.4} m/z, {:.3e}</title>{}</text>",
            escape_xml(&row.label),
            row.mz,
            row.summed_intensity,
            escape_xml(&row.label),
        );
    }
}

fn draw_pseudo_ms2_legend(svg: &mut String, x: f64, y: f64) {
    for (idx, (label, color)) in [("b ions", "#1d4ed8"), ("y ions", "#b45309")]
        .iter()
        .enumerate()
    {
        let y = y + idx as f64 * 18.0;
        let _ = writeln!(
            svg,
            "<line x1=\"{x:.2}\" y1=\"{y:.2}\" x2=\"{x2:.2}\" y2=\"{y:.2}\" stroke=\"{color}\" stroke-width=\"3\"/>",
            x2 = x + 22.0
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"14\" dominant-baseline=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
            x + 30.0,
            y,
            COLOR_SUBTLE,
            label
        );
    }
}

fn pseudo_ms2_series_color(series: &str) -> &'static str {
    match series {
        "b" => "#1d4ed8",
        "y" => "#b45309",
        _ => "#475569",
    }
}

fn write_summary_text(summary: &DiaSliceSummary, request: &DiaSliceRequest) -> anyhow::Result<()> {
    let path = output_path(&summary.out_prefix, "summary.txt");
    let mut file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    writeln!(file, "backend\t{}", summary.backend_label)?;
    writeln!(file, "input\t{}", summary.input_path.display())?;
    writeln!(file, "mz_center\t{:.6}", request.mz)?;
    writeln!(file, "mz_min\t{:.6}", summary.mz_min)?;
    writeln!(file, "mz_max\t{:.6}", summary.mz_max)?;
    if let Some(target) = &request.peptide_target {
        writeln!(file, "peptide_input\t{}", sanitize_tsv(&target.input))?;
        writeln!(file, "peptide_sequence\t{}", target.sequence)?;
        writeln!(
            file,
            "peptide_modified_sequence\t{}",
            sanitize_tsv(&target.modified_sequence)
        )?;
        writeln!(file, "peptide_charge\t{}", target.charge)?;
        writeln!(file, "peptide_precursor_mz\t{:.6}", target.precursor_mz)?;
        if let Some(fragment) = &target.fragment {
            writeln!(
                file,
                "peptide_fragment_input\t{}",
                sanitize_tsv(&fragment.input)
            )?;
            writeln!(file, "peptide_fragment_label\t{}", fragment.label)?;
            writeln!(file, "peptide_fragment_mz\t{:.6}", fragment.mz)?;
        }
    }
    writeln!(file, "rt_window\t{}", request.rt_window_label())?;
    writeln!(file, "spectra_considered\t{}", summary.spectra_considered)?;
    writeln!(file, "spectra_with_signal\t{}", summary.spectra_with_signal)?;
    writeln!(file, "matched_peaks\t{}", summary.matched_peaks)?;
    writeln!(file, "capabilities\t{}", summary.capabilities.labels())?;
    if let Some(mode) = &summary.acquisition_mode {
        writeln!(file, "acquisition_mode\t{}", mode)?;
    }
    if let Some(column) = &summary.intensity_column {
        writeln!(file, "intensity_column\t{}", column)?;
    }
    if let Some(path) = &summary.vendor_runtime_path {
        writeln!(file, "vendor_runtime\t{}", path.display())?;
    }
    Ok(())
}

fn write_rt_profile(summary: &DiaSliceSummary, rt_rows: &[RtProfileRow]) -> anyhow::Result<()> {
    let path = output_path(&summary.out_prefix, "rt_profile.tsv");
    let mut file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    writeln!(
        file,
        "scan_index\tscan_id\trt_minutes\tsummed_intensity\tmatched_peaks\tprecursor_mz"
    )?;
    for row in rt_rows {
        writeln!(
            file,
            "{}\t{}\t{}\t{:.6}\t{}\t{}",
            row.scan_index,
            sanitize_tsv(&row.scan_id),
            format_optional_f64(row.rt_minutes, 6),
            row.summed_intensity,
            row.matched_peaks,
            format_optional_f64(row.precursor_mz, 6),
        )?;
    }
    Ok(())
}

fn write_mz_profile(summary: &DiaSliceSummary, mz_bins: &[MzProfileBin]) -> anyhow::Result<()> {
    let path = output_path(&summary.out_prefix, "mz_profile.tsv");
    let mut file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    writeln!(file, "mz_center\tsummed_intensity")?;
    for bin in mz_bins {
        writeln!(file, "{:.6}\t{:.6}", bin.mz_center, bin.summed_intensity)?;
    }
    Ok(())
}

fn write_im_profile(
    summary: &DiaSliceSummary,
    im_rows: Option<&[MobilityProfileRow]>,
) -> anyhow::Result<()> {
    let Some(im_rows) = im_rows else {
        return Ok(());
    };
    if im_rows.is_empty() {
        return Ok(());
    }

    let path = output_path(&summary.out_prefix, "im_profile.tsv");
    let mut file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    writeln!(file, "mobility\tsummed_intensity")?;
    for row in im_rows {
        writeln!(file, "{:.6}\t{:.6}", row.mobility, row.summed_intensity)?;
    }
    Ok(())
}

fn write_summary_svg(
    summary: &DiaSliceSummary,
    request: &DiaSliceRequest,
    rt_rows: &[RtProfileRow],
    mz_bins: &[MzProfileBin],
    im_rows: Option<&[MobilityProfileRow]>,
    rt_smooth: bool,
    rt_smooth_window: usize,
) -> anyhow::Result<()> {
    let path = output_path(&summary.out_prefix, "svg");
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">",
        w = SVG_WIDTH,
        h = SVG_HEIGHT
    );
    let _ = writeln!(
        svg,
        "<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"{}\"/>",
        COLOR_BG
    );

    let source_label = compact_display_path(&summary.input_path, 78);
    let title = "DIA slice summary";
    let source_line = format!("Source: {source_label}");
    let filter_line = match &request.peptide_target {
        Some(target) => match &target.fragment {
            Some(fragment) => format!(
                "Target: {}/{} fragment {} m/z {:.4} | precursor m/z {:.4} | window {:.4}-{:.4} | RT {}",
                target.modified_sequence,
                target.charge,
                fragment.label,
                fragment.mz,
                target.precursor_mz,
                summary.mz_min,
                summary.mz_max,
                request.rt_window_label(),
            ),
            None => format!(
                "Target: {}/{} precursor m/z {:.4} | window {:.4}-{:.4} | RT {}",
                target.modified_sequence,
                target.charge,
                target.precursor_mz,
                summary.mz_min,
                summary.mz_max,
                request.rt_window_label(),
            ),
        },
        None => format!(
            "Target: m/z {:.4} | window {:.4}-{:.4} | RT {}",
            request.mz,
            summary.mz_min,
            summary.mz_max,
            request.rt_window_label(),
        ),
    };
    let data_context =
        compact_backend_label(summary.backend_label, summary.acquisition_mode.as_deref());
    let capability_context = compact_capability_label(summary.capabilities);
    let acquisition_line = format!("{filter_line} | {data_context} | {capability_context}",);
    let counts_line = format!(
        "Spectra considered: {} | With signal: {} | Matched peaks: {}",
        summary.spectra_considered, summary.spectra_with_signal, summary.matched_peaks,
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        44.0,
        SVG_TITLE_FONT,
        30.0,
        COLOR_TEXT,
        Some("700"),
        &[title.to_string()],
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        70.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &source_line,
            SVG_META_FONT,
            (SVG_WIDTH as f64) - 2.0 * SVG_MARGIN_X,
            1,
        ),
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        96.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &acquisition_line,
            SVG_META_FONT,
            (SVG_WIDTH as f64) - 2.0 * SVG_MARGIN_X,
            1,
        ),
    );
    append_mono_text_lines(
        &mut svg,
        SVG_MARGIN_X,
        120.0,
        SVG_META_FONT,
        20.0,
        COLOR_SUBTLE,
        None,
        &wrap_mono_text(
            &counts_line,
            SVG_META_FONT,
            (SVG_WIDTH as f64) - 2.0 * SVG_MARGIN_X,
            1,
        ),
    );

    let margin_left = 88.0;
    let panel_top = 174.0;
    let panel_width = 1100.0;
    let chart_count = if im_rows.is_some_and(|rows| !rows.is_empty()) {
        3
    } else {
        2
    };
    let panel_height = if chart_count == 3 { 176.0 } else { 280.0 };
    let panel_gap = if chart_count == 3 { 104.0 } else { 112.0 };
    let rt_x_domain = rt_x_domain(rt_rows);
    let rt_y_domain = padded_y_range(rt_rows.iter().map(|row| row.summed_intensity));
    let rt_points = rt_points(rt_rows);
    let rt_smoothed_points =
        rt_smooth.then(|| gaussian_smooth_points(&rt_points, rt_smooth_window));
    let rt_tic = rt_rows.iter().map(|row| row.summed_intensity).sum::<f64>();
    let rt_annotation = if rt_smoothed_points.is_some() {
        format!("TIC {} | smoothed", format_compact_intensity(rt_tic))
    } else {
        format!("TIC {}", format_compact_intensity(rt_tic))
    };
    let rt_canvas = SvgCanvas::new(
        margin_left,
        panel_top,
        panel_width,
        panel_height,
        rt_x_domain,
        rt_y_domain,
    );
    let mz_x_domain = mz_x_domain(mz_bins);
    let mz_y_domain = padded_y_range(mz_bins.iter().map(|bin| bin.summed_intensity));
    let mz_top = if chart_count == 3 {
        panel_top + 2.0 * (panel_height + panel_gap)
    } else {
        panel_top + panel_height + panel_gap
    };
    let mz_canvas = SvgCanvas::new(
        margin_left,
        mz_top,
        panel_width,
        panel_height,
        mz_x_domain,
        mz_y_domain,
    );

    append_chart(
        &mut svg,
        rt_canvas,
        rt_x_domain,
        rt_y_domain,
        "RT profile",
        if rt_rows.iter().all(|row| row.rt_minutes.is_some()) {
            AxisProps::new(AxisOrientation::Bottom, "RT (min)")
                .with_tick_label_style(AxisTickLabelStyle::Precision(2))
        } else {
            AxisProps::new(AxisOrientation::Bottom, "scan index")
                .with_tick_label_style(AxisTickLabelStyle::Precision(0))
        },
        AxisProps::new(AxisOrientation::Left, "summed intensity")
            .with_tick_label_style(AxisTickLabelStyle::Scientific(2)),
        &rt_points,
        rt_smoothed_points.as_deref(),
        Some(&rt_annotation),
        COLOR_RT,
    );
    if let Some(im_rows) = im_rows.filter(|rows| !rows.is_empty()) {
        let im_x_domain = mobility_x_domain(im_rows);
        let im_y_domain = padded_y_range(im_rows.iter().map(|row| row.summed_intensity));
        let im_canvas = SvgCanvas::new(
            margin_left,
            panel_top + (panel_height + panel_gap),
            panel_width,
            panel_height,
            im_x_domain,
            im_y_domain,
        );
        append_chart(
            &mut svg,
            im_canvas,
            im_x_domain,
            im_y_domain,
            "Mobility profile",
            AxisProps::new(AxisOrientation::Bottom, "1/K0")
                .with_tick_label_style(AxisTickLabelStyle::Precision(4)),
            AxisProps::new(AxisOrientation::Left, "summed intensity")
                .with_tick_label_style(AxisTickLabelStyle::Scientific(2)),
            &mobility_points(im_rows),
            None,
            None,
            "#0f766e",
        );
    }
    append_chart(
        &mut svg,
        mz_canvas,
        mz_x_domain,
        mz_y_domain,
        "m/z profile",
        AxisProps::new(AxisOrientation::Bottom, "m/z")
            .with_tick_label_style(AxisTickLabelStyle::Precision(4)),
        AxisProps::new(AxisOrientation::Left, "summed intensity")
            .with_tick_label_style(AxisTickLabelStyle::Scientific(2)),
        &mz_points(mz_bins),
        None,
        None,
        COLOR_MZ,
    );

    svg.push_str("</svg>\n");
    fs::write(&path, svg).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn append_chart(
    svg: &mut String,
    canvas: SvgCanvas,
    x_domain: CoordinateRange<f64>,
    y_domain: CoordinateRange<f64>,
    title: &str,
    x_axis: AxisProps,
    y_axis: AxisProps,
    points: &[(f64, f64)],
    overlay_points: Option<&[(f64, f64)]>,
    top_right_label: Option<&str>,
    color: &str,
) {
    let _ = writeln!(
        svg,
        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"12\" fill=\"#fbfdff\" stroke=\"{}\" stroke-width=\"1\"/>",
        canvas.left() - 14.0,
        canvas.top() - 38.0,
        canvas.width() + 28.0,
        canvas.height() + 112.0,
        COLOR_CARD_BORDER
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
        canvas.left(),
        canvas.top() - 16.0,
        SVG_PANEL_TITLE_FONT,
        COLOR_TEXT,
        escape_xml(title)
    );
    if let Some(label) = top_right_label {
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"12\" text-anchor=\"end\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
            canvas.right(),
            canvas.top() - 16.0,
            COLOR_SUBTLE,
            escape_xml(label)
        );
    } else if overlay_points.is_some_and(|points| points.len() > 1) {
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"12\" text-anchor=\"end\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">smoothed overlay</text>",
            canvas.right(),
            canvas.top() - 16.0,
            COLOR_SUBTLE
        );
    }

    let x_ticks = linear_ticks(x_domain.min(), x_domain.max(), 5);
    let y_ticks = linear_ticks(y_domain.min(), y_domain.max(), 5);
    for tick in x_ticks {
        let px = canvas.x(tick);
        let _ = writeln!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>",
            px,
            canvas.top(),
            px,
            canvas.bottom(),
            COLOR_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
            px,
            canvas.bottom() + 28.0,
            SVG_TICK_FONT,
            COLOR_SUBTLE,
            x_axis.format_tick(tick)
        );
    }
    for tick in y_ticks {
        let py = canvas.y(tick);
        let _ = writeln!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"1\"/>",
            canvas.left(),
            py,
            canvas.right(),
            py,
            COLOR_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"end\" dominant-baseline=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
            canvas.left() - 12.0,
            py,
            SVG_TICK_FONT,
            COLOR_SUBTLE,
            y_axis.format_tick(tick)
        );
    }

    let _ = writeln!(
        svg,
        "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.2\"/>",
        canvas.left(),
        canvas.top(),
        canvas.width(),
        canvas.height(),
        COLOR_AXIS
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\">{}</text>",
        canvas.left() + canvas.width() / 2.0,
        canvas.bottom() + 58.0,
        SVG_AXIS_LABEL_FONT,
        COLOR_TEXT,
        escape_xml(x_axis.label())
    );
    let _ = writeln!(
        svg,
        "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{:.1}\" text-anchor=\"middle\" fill=\"{}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\" transform=\"rotate(-90 {:.2} {:.2})\">{}</text>",
        canvas.left() - 64.0,
        canvas.top() + canvas.height() / 2.0,
        SVG_AXIS_LABEL_FONT,
        COLOR_SUBTLE,
        canvas.left() - 64.0,
        canvas.top() + canvas.height() / 2.0,
        escape_xml(y_axis.label())
    );

    if points.is_empty() {
        return;
    }
    if points.len() == 1 {
        let (x, y) = canvas.transform(points[0].0, points[0].1);
        let _ = writeln!(
            svg,
            "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"4\" fill=\"{}\"/>",
            x, y, color
        );
        return;
    }

    let path_data = chart_path_data(canvas, points);
    let raw_opacity = if overlay_points.is_some() { 0.42 } else { 1.0 };
    let raw_width = if overlay_points.is_some() { 1.4 } else { 2.2 };
    let _ = writeln!(
        svg,
        "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"{:.1}\" stroke-opacity=\"{:.2}\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>",
        path_data,
        color,
        raw_width,
        raw_opacity,
    );
    if let Some(overlay_points) = overlay_points.filter(|points| points.len() > 1) {
        let overlay_path = chart_path_data(canvas, overlay_points);
        let _ = writeln!(
            svg,
            "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"3.2\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>",
            overlay_path,
            color
        );
    }
}

fn chart_path_data(canvas: SvgCanvas, points: &[(f64, f64)]) -> String {
    let mut path_data = String::new();
    for (idx, (x, y)) in points.iter().enumerate() {
        let (px, py) = canvas.transform(*x, *y);
        let command = if idx == 0 { 'M' } else { 'L' };
        let _ = write!(path_data, "{command}{px:.2},{py:.2} ");
    }
    path_data.trim_end().to_string()
}

fn format_compact_intensity(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.3e}")
    } else {
        "NA".to_string()
    }
}

fn append_mono_text_lines(
    svg: &mut String,
    x: f64,
    first_y: f64,
    font_size: f64,
    line_height: f64,
    fill: &str,
    font_weight: Option<&str>,
    lines: &[String],
) {
    let weight_attr = font_weight
        .map(|weight| format!(" font-weight=\"{weight}\""))
        .unwrap_or_default();
    for (idx, line) in lines.iter().enumerate() {
        let _ = writeln!(
            svg,
            "<text x=\"{x:.2}\" y=\"{y:.2}\" font-size=\"{font_size:.1}\" fill=\"{fill}\" font-family=\"Menlo, Consolas, Liberation Mono, monospace\"{weight_attr}>{text}</text>",
            y = first_y + idx as f64 * line_height,
            text = escape_xml(line)
        );
    }
}

fn wrap_mono_text(text: &str, font_size: f64, max_width: f64, max_lines: usize) -> Vec<String> {
    let max_chars = mono_char_budget(font_size, max_width);
    if max_lines == 0 || max_chars == 0 {
        return Vec::new();
    }

    let mut lines = Vec::<String>::new();
    let mut current = String::new();
    for word in text.split_whitespace() {
        let candidate_len = if current.is_empty() {
            word.chars().count()
        } else {
            current.chars().count() + 1 + word.chars().count()
        };
        if candidate_len <= max_chars {
            if !current.is_empty() {
                current.push(' ');
            }
            current.push_str(word);
            continue;
        }

        if !current.is_empty() {
            lines.push(current);
            current = String::new();
            if lines.len() == max_lines {
                break;
            }
        }

        if word.chars().count() > max_chars {
            lines.push(truncate_middle(word, max_chars));
            if lines.len() == max_lines {
                break;
            }
        } else {
            current.push_str(word);
        }
    }

    if !current.is_empty() && lines.len() < max_lines {
        lines.push(current);
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    if lines.len() == max_lines && text_is_not_fully_represented(text, &lines) {
        if let Some(last) = lines.last_mut() {
            *last = truncate_middle(last, max_chars.saturating_sub(3));
            last.push_str("...");
        }
    }
    lines
}

fn text_is_not_fully_represented(original: &str, lines: &[String]) -> bool {
    let rendered = lines.join(" ");
    rendered != original && !rendered.ends_with("...")
}

fn mono_char_budget(font_size: f64, max_width: f64) -> usize {
    (max_width / (font_size * 0.62)).floor().max(1.0) as usize
}

fn compact_display_path(path: &Path, max_chars: usize) -> String {
    let cwd = env::current_dir().ok();
    let home = env::var_os("HOME").map(PathBuf::from);
    compact_display_path_with(path, cwd.as_deref(), home.as_deref(), max_chars)
}

fn compact_display_path_with(
    path: &Path,
    cwd: Option<&Path>,
    home: Option<&Path>,
    max_chars: usize,
) -> String {
    let display = if path.is_relative() {
        path.to_string_lossy().to_string()
    } else if let Some(cwd) = cwd.and_then(|cwd| path.strip_prefix(cwd).ok()) {
        cwd.to_string_lossy().to_string()
    } else if let Some(home_relative) = home.and_then(|home| path.strip_prefix(home).ok()) {
        let rest = home_relative.to_string_lossy();
        if rest.is_empty() {
            "~".to_string()
        } else {
            format!("~/{rest}")
        }
    } else {
        path.to_string_lossy().to_string()
    };
    collapse_middle_path(&display, max_chars)
}

fn collapse_middle_path(value: &str, max_chars: usize) -> String {
    if value.chars().count() <= max_chars {
        return value.to_string();
    }
    let parts = value.split('/').collect::<Vec<_>>();
    if parts.len() >= 4 {
        let prefix_count = if parts.first() == Some(&"~") {
            3.min(parts.len().saturating_sub(1))
        } else {
            2.min(parts.len().saturating_sub(1))
        };
        let prefix = parts[..prefix_count].join("/");
        let suffix = parts.last().copied().unwrap_or_default();
        let candidate = format!("{prefix}/.../{suffix}");
        if candidate.chars().count() <= max_chars {
            return candidate;
        }

        let reserved = prefix.chars().count() + "/.../".chars().count();
        let suffix_budget = max_chars.saturating_sub(reserved).max(8);
        return format!("{prefix}/.../{}", truncate_middle(suffix, suffix_budget));
    }
    truncate_middle(value, max_chars)
}

fn truncate_middle(value: &str, max_chars: usize) -> String {
    let count = value.chars().count();
    if count <= max_chars {
        return value.to_string();
    }
    if max_chars <= 3 {
        return ".".repeat(max_chars);
    }
    let remaining = max_chars - 3;
    let left_count = (remaining + 1) / 2;
    let right_count = remaining / 2;
    let left = value.chars().take(left_count).collect::<String>();
    let right = value
        .chars()
        .rev()
        .take(right_count)
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect::<String>();
    format!("{left}...{right}")
}

fn truncate_end(value: &str, max_chars: usize) -> String {
    let count = value.chars().count();
    if count <= max_chars {
        return value.to_string();
    }
    if max_chars <= 3 {
        return ".".repeat(max_chars);
    }
    let keep = max_chars - 3;
    let prefix = value.chars().take(keep).collect::<String>();
    format!("{prefix}...")
}

fn gaussian_smooth_points(points: &[(f64, f64)], window_points: usize) -> Vec<(f64, f64)> {
    if points.len() < 3 || window_points <= 1 {
        return points.to_vec();
    }
    let window = if window_points % 2 == 0 {
        window_points + 1
    } else {
        window_points
    }
    .min(points.len().saturating_mul(2).saturating_sub(1));
    let half = window / 2;
    let sigma = (window as f64 / 3.0).max(1.0);
    let two_sigma_sq = 2.0 * sigma * sigma;

    points
        .iter()
        .enumerate()
        .map(|(idx, (x, _))| {
            let start = idx.saturating_sub(half);
            let end = (idx + half + 1).min(points.len());
            let mut weighted_sum = 0.0;
            let mut weight_total = 0.0;
            for (neighbor_idx, (_, y)) in points.iter().enumerate().take(end).skip(start) {
                let distance = neighbor_idx.abs_diff(idx) as f64;
                let mut weight = (-distance * distance / two_sigma_sq).exp();
                if *y <= 0.0 {
                    weight *= RT_SMOOTH_ZERO_WEIGHT;
                }
                weighted_sum += *y * weight;
                weight_total += weight;
            }
            let smoothed = if weight_total > 0.0 {
                weighted_sum / weight_total
            } else {
                0.0
            };
            (*x, smoothed)
        })
        .collect()
}

fn rt_x_domain(rt_rows: &[RtProfileRow]) -> CoordinateRange<f64> {
    let values = rt_rows
        .iter()
        .map(|row| row.rt_minutes.unwrap_or(row.scan_index as f64));
    padded_range(values)
}

fn mz_x_domain(mz_bins: &[MzProfileBin]) -> CoordinateRange<f64> {
    padded_range(mz_bins.iter().map(|bin| bin.mz_center))
}

fn mobility_x_domain(im_rows: &[MobilityProfileRow]) -> CoordinateRange<f64> {
    padded_range(im_rows.iter().map(|row| row.mobility))
}

fn padded_y_range(values: impl Iterator<Item = f64>) -> CoordinateRange<f64> {
    let max = values.fold(0.0f64, f64::max).max(1.0);
    CoordinateRange::new(0.0, max * 1.05)
}

fn padded_range(values: impl Iterator<Item = f64>) -> CoordinateRange<f64> {
    let values = values.filter(|value| value.is_finite()).collect::<Vec<_>>();
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !min.is_finite() || !max.is_finite() {
        return CoordinateRange::new(0.0, 1.0);
    }
    if (max - min).abs() <= f64::EPSILON {
        let pad = max.abs().max(1.0) * 0.01;
        return CoordinateRange::new(min - pad, max + pad);
    }
    CoordinateRange::new(min, max)
}

fn rt_points(rt_rows: &[RtProfileRow]) -> Vec<(f64, f64)> {
    rt_rows
        .iter()
        .map(|row| {
            (
                row.rt_minutes.unwrap_or(row.scan_index as f64),
                row.summed_intensity,
            )
        })
        .collect()
}

fn mz_points(mz_bins: &[MzProfileBin]) -> Vec<(f64, f64)> {
    mz_bins
        .iter()
        .map(|bin| (bin.mz_center, bin.summed_intensity))
        .collect()
}

fn mobility_points(im_rows: &[MobilityProfileRow]) -> Vec<(f64, f64)> {
    im_rows
        .iter()
        .map(|row| (row.mobility, row.summed_intensity))
        .collect()
}

fn linear_ticks(start: f64, end: f64, count: usize) -> Vec<f64> {
    if count <= 1 {
        return vec![start];
    }
    let step = (end - start) / (count - 1) as f64;
    (0..count).map(|idx| start + step * idx as f64).collect()
}

fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    let file_name = prefix
        .file_name()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .unwrap_or("dia_slice");
    let path = prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(format!("{file_name}.{suffix}")));
    path.unwrap_or_else(|| PathBuf::from(format!("{file_name}.{suffix}")))
}

fn write_temp_bridge_script() -> anyhow::Result<PathBuf> {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let path = env::temp_dir().join(format!(
        "mzio-bruker-bridge-{}-{timestamp}.py",
        std::process::id()
    ));
    fs::write(&path, BRUKER_BRIDGE_PY)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(path)
}

fn sanitize_tsv(value: &str) -> String {
    value.replace('\t', " ").replace('\n', " ")
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
        "target".to_string()
    } else {
        out
    }
}

fn format_optional_f64(value: Option<f64>, precision: usize) -> String {
    match value {
        Some(value) => format!("{value:.precision$}"),
        None => String::new(),
    }
}

fn escape_xml(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use crate::annotate::NeutralLossKind;
    use crate::scale::CoordinateRange;
    use crate::svg_canvas::{AxisOrientation, AxisProps, SvgCanvas};
    use timsrust::{Frame, MSLevel};

    use super::{
        accumulate_mz_bin, compact_display_path_with, initialize_mz_bins, padded_range, parse_args,
    };

    #[test]
    fn parse_args_rejects_mobility_flags_for_mzml() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--im-min".into(),
            "0.7".into(),
            "--im-max".into(),
            "1.1".into(),
        ])
        .expect_err("expected mobility flags to be rejected for mzML");
        assert!(err
            .to_string()
            .contains("--mzml input does not support --im-* or --quad-* yet"));
    }

    #[test]
    fn parse_args_rejects_python_for_mzml() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--python".into(),
            "python3".into(),
        ])
        .expect_err("expected --python to be rejected for mzML");
        assert!(err
            .to_string()
            .contains("--python is only valid together with --bruker <run.d>"));
    }

    #[test]
    fn initialize_and_fill_mz_bins() {
        let mut bins = vec![
            super::MzProfileBin {
                mz_center: 0.0,
                summed_intensity: 0.0,
            };
            3
        ];
        initialize_mz_bins(&mut bins, 100.0, 106.0);
        assert!((bins[0].mz_center - 100.0).abs() < 1e-9);
        assert!((bins[1].mz_center - 103.0).abs() < 1e-9);
        assert!((bins[2].mz_center - 106.0).abs() < 1e-9);

        accumulate_mz_bin(&mut bins, 100.0, 106.0, 105.8, 42.0);
        assert!((bins[2].summed_intensity - 42.0).abs() < 1e-9);
    }

    #[test]
    fn padded_range_handles_single_value() {
        let range = padded_range([42.0].into_iter());
        assert!(range.start < 42.0);
        assert!(range.end > 42.0);
    }

    #[test]
    fn progress_bar_body_scales_to_requested_width() {
        assert_eq!(super::progress_bar_body(0, 10, 10), "..........");
        assert_eq!(super::progress_bar_body(5, 10, 10), "=====.....");
        assert_eq!(super::progress_bar_body(10, 10, 10), "==========");
        assert_eq!(super::progress_percent(3, 4), 75);
    }

    #[test]
    fn bruker_run_tic_accumulator_splits_ms_levels() {
        let mut acc = super::BrukerRunTicAccumulator::new();
        let ms1 = Frame {
            index: 7,
            rt_in_seconds: 30.0,
            ms_level: MSLevel::MS1,
            intensity_correction_factor: 2.0,
            intensities: vec![10, 20],
            ..Frame::default()
        };
        let ms2 = Frame {
            index: 8,
            rt_in_seconds: 90.0,
            ms_level: MSLevel::MS2,
            intensity_correction_factor: 0.5,
            intensities: vec![8, 12],
            ..Frame::default()
        };

        acc.add_frame(&ms1);
        acc.add_frame(&ms2);
        let summary = acc.finalize(Path::new("run.d"), "diaPASEF".to_string());

        assert_eq!(summary.schema_version, 1);
        assert_eq!(summary.total_frames, 2);
        assert_eq!(summary.acquisition_mode, "diaPASEF");
        assert_eq!(summary.ms1.frames, 1);
        assert_eq!(summary.ms1.detector_events, 2);
        assert_eq!(summary.ms1.max_tic_frame, Some(7));
        assert!((summary.ms1.summed_tic - 60.0).abs() < 1e-9);
        assert_eq!(summary.ms2.frames, 1);
        assert_eq!(summary.ms2.max_tic_frame, Some(8));
        assert!((summary.ms2.summed_tic - 10.0).abs() < 1e-9);
        assert_eq!(summary.unknown.frames, 0);
        assert!((summary.rt_min_minutes.expect("rt min") - 0.5).abs() < 1e-9);
        assert!((summary.rt_max_minutes.expect("rt max") - 1.5).abs() < 1e-9);

        let json = serde_json::to_string(&summary).expect("serialize run TIC summary");
        assert!(json.contains("\"ms1\""));
        assert!(json.contains("\"ms2\""));
    }

    #[test]
    fn parse_args_accepts_peptide_target_without_mz() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
        ])
        .expect("options parse");
        let request = options.request.expect("request");
        let target = request.peptide_target.expect("peptide target");
        assert_eq!(target.sequence, "GAIIGLMVGGVVIA");
        assert_eq!(target.charge, 2);
        assert!((request.mz - 635.3836).abs() < 0.0001);
        assert!((target.precursor_mz - request.mz).abs() < 1e-9);
    }

    #[test]
    fn parse_args_accepts_peptide_fragment_target_without_mz() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
            "--fragment".into(),
            "b8".into(),
        ])
        .expect("options parse");
        let request = options.request.expect("request");
        let target = request.peptide_target.expect("peptide target");
        let fragment = target.fragment.expect("fragment target");
        assert_eq!(fragment.label, "b8");
        assert!((target.precursor_mz - 635.3836).abs() < 0.0001);
        assert!((fragment.mz - 755.4484).abs() < 0.0001);
        assert!((request.mz - fragment.mz).abs() < 1e-9);
    }

    #[test]
    fn parse_args_accepts_pseudo_ms2_for_bruker_peptide() {
        let options = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
            "--ms2".into(),
            "--pseudo-ms2-rt-window".into(),
            "0.75".into(),
        ])
        .expect("pseudo-ms2 options parse");
        assert!(options.pseudo_ms2);
        assert!((options.pseudo_ms2_rt_window_min - 0.75).abs() < 1e-9);
        let target = options
            .request
            .expect("request")
            .peptide_target
            .expect("peptide target");
        assert!(target.fragment.is_none());
        assert!(target
            .fragments
            .iter()
            .any(|fragment| fragment.label == "b8"));
    }

    #[test]
    fn parse_args_accepts_pseudo_ms2_with_neutral_loss_fragments() {
        let options = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--peptide".into(),
            "S[+79.9663]PEPTIDE/2".into(),
            "--pseudo-ms2".into(),
            "--neutral-losses".into(),
        ])
        .expect("pseudo-ms2 neutral-loss options parse");
        assert!(options.neutral_losses_enabled);
        let target = options
            .request
            .expect("request")
            .peptide_target
            .expect("peptide target");
        assert!(target
            .fragments
            .iter()
            .any(|fragment| fragment.label.contains("-H3PO4")));
    }

    #[test]
    fn parse_args_rejects_neutral_losses_without_peptide() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--neutral-losses".into(),
        ])
        .expect_err("neutral losses without peptide should fail");
        assert!(err
            .to_string()
            .contains("--neutral-losses requires --peptide"));
    }

    #[test]
    fn parse_args_accepts_rt_smooth_window() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--rt-smooth-window".into(),
            "5".into(),
        ])
        .expect("rt smoothing options parse");
        assert!(options.rt_smooth);
        assert_eq!(options.rt_smooth_window, 5);
    }

    #[test]
    fn gaussian_smooth_points_preserves_x_and_reduces_spike() {
        let points = vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0)];
        let smoothed = super::gaussian_smooth_points(&points, 3);
        assert_eq!(smoothed.len(), points.len());
        assert_eq!(smoothed[1].0, 1.0);
        assert!(smoothed[1].1 > 6.0);
        assert!(smoothed[1].1 < 10.0);
        assert!(smoothed[0].1 > 0.0);
    }

    #[test]
    fn append_chart_draws_smoothed_overlay_when_points_are_available() {
        let canvas = SvgCanvas::new(
            0.0,
            0.0,
            100.0,
            80.0,
            CoordinateRange::new(0.0, 2.0),
            CoordinateRange::new(0.0, 10.0),
        );
        let points = vec![(0.0, 0.0), (1.0, 10.0), (2.0, 0.0)];
        let smoothed = super::gaussian_smooth_points(&points, 3);
        let mut svg = String::new();
        super::append_chart(
            &mut svg,
            canvas,
            CoordinateRange::new(0.0, 2.0),
            CoordinateRange::new(0.0, 10.0),
            "RT profile",
            AxisProps::new(AxisOrientation::Bottom, "RT (min)"),
            AxisProps::new(AxisOrientation::Left, "summed intensity"),
            &points,
            Some(&smoothed),
            Some("TIC 1.000e1 | smoothed"),
            "#1d4ed8",
        );
        assert!(svg.contains("TIC 1.000e1 | smoothed"));
        assert!(!svg.contains("smoothed overlay"));
        assert!(svg.contains("stroke-opacity=\"0.42\""));
    }

    #[test]
    fn pseudo_ms2_ion_table_marks_y_ion_evidence() {
        let peptide = super::resolve_peptide_target("PEPTIDE/3", None, &[], None, &[])
            .expect("peptide target");
        let fragments = peptide
            .fragments
            .iter()
            .map(|fragment| super::PseudoMs2FragmentEvidence {
                label: fragment.label.clone(),
                series: fragment.series.clone(),
                cleavage_index: fragment.cleavage_index,
                neutral_loss: fragment.neutral_loss,
                charge: fragment.charge,
                mz: fragment.mz,
                summed_intensity: if fragment.label == "y3" { 42.0 } else { 0.0 },
                matched_events: if fragment.label == "y3" { 2 } else { 0 },
                frames_with_signal: if fragment.label == "y3" { 1 } else { 0 },
                apex_rt: if fragment.label == "y3" {
                    Some(12.5)
                } else {
                    None
                },
                apex_intensity: if fragment.label == "y3" { 42.0 } else { 0.0 },
            })
            .collect::<Vec<_>>();
        let report = super::PseudoMs2Report {
            input_path: PathBuf::from("run.d"),
            out_prefix: PathBuf::from("out"),
            peptide,
            rt_window: super::PseudoMs2RtWindow {
                min: 12.0,
                max: 13.0,
                source: "test",
            },
            frames_considered: 1,
            frames_with_signal: 1,
            matched_events: 2,
            precursor_frames_with_signal: 1,
            precursor_apex_rt: Some(12.5),
            precursor_apex_intensity: 100.0,
            fragments,
        };

        let rows = super::build_pseudo_ion_table_rows(&report);
        let row = rows.iter().find(|row| row.y_ordinal == 3).expect("y3 row");
        assert!(row.y1.iter().any(super::PseudoIonTableCell::detected));
        assert!(!row.b1.iter().any(super::PseudoIonTableCell::detected));

        let mut svg = String::new();
        let width = super::pseudo_ion_table_width(&rows);
        let height = super::pseudo_ion_table_height(rows.len());
        super::draw_pseudo_ms2_ion_table(
            &mut svg,
            0.0,
            0.0,
            width,
            height,
            &rows,
            &report.peptide.sequence,
            0.0,
        );
        assert!(svg.contains("Pseudo-MS2 ion table"));
        assert!(svg.contains("color = DIA evidence"));
        assert!(svg.contains("y3 theoretical"));
        assert!(svg.contains("neutral losses with at least 3% of max fragment signal"));
    }

    #[test]
    fn pseudo_ms2_ion_table_height_keeps_single_long_ladder() {
        let row_height = super::pseudo_ion_table_row_height(13);
        assert!((super::pseudo_ion_table_height(13) - (180.0 + 13.0 * row_height)).abs() < 1e-9);
    }

    #[test]
    fn pseudo_ms2_stick_labels_include_low_intensity_fragments() {
        let fragments = vec![
            super::PseudoMs2FragmentEvidence {
                label: "b1".to_string(),
                series: "b".to_string(),
                cleavage_index: 1,
                neutral_loss: None,
                charge: 1,
                mz: 100.0,
                summed_intensity: 1000.0,
                matched_events: 10,
                frames_with_signal: 5,
                apex_rt: Some(1.0),
                apex_intensity: 1000.0,
            },
            super::PseudoMs2FragmentEvidence {
                label: "y1".to_string(),
                series: "y".to_string(),
                cleavage_index: 1,
                neutral_loss: None,
                charge: 1,
                mz: 300.0,
                summed_intensity: 1.0,
                matched_events: 1,
                frames_with_signal: 1,
                apex_rt: Some(1.1),
                apex_intensity: 1.0,
            },
            super::PseudoMs2FragmentEvidence {
                label: "b2".to_string(),
                series: "b".to_string(),
                cleavage_index: 2,
                neutral_loss: None,
                charge: 1,
                mz: 320.0,
                summed_intensity: 0.0,
                matched_events: 0,
                frames_with_signal: 0,
                apex_rt: None,
                apex_intensity: 0.0,
            },
        ];
        let canvas = SvgCanvas::new(
            0.0,
            0.0,
            400.0,
            200.0,
            CoordinateRange::new(50.0, 350.0),
            CoordinateRange::new(0.0, 1000.0),
        );
        let mut svg = String::new();
        super::draw_pseudo_ms2_sticks(&mut svg, canvas, &fragments, 0.0);

        assert!(svg.contains(">b1</text>"));
        assert!(svg.contains(">y1</text>"));
        assert!(!svg.contains(">b2</text>"));
        assert!(!svg.contains("#aeb8c4"));
    }

    #[test]
    fn pseudo_ms2_svg_hides_weak_neutral_loss_labels() {
        let fragments = vec![
            super::PseudoMs2FragmentEvidence {
                label: "b1".to_string(),
                series: "b".to_string(),
                cleavage_index: 1,
                neutral_loss: None,
                charge: 1,
                mz: 100.0,
                summed_intensity: 1000.0,
                matched_events: 10,
                frames_with_signal: 5,
                apex_rt: Some(1.0),
                apex_intensity: 1000.0,
            },
            super::PseudoMs2FragmentEvidence {
                label: "b1-H3PO4".to_string(),
                series: "b".to_string(),
                cleavage_index: 1,
                neutral_loss: Some(NeutralLossKind::PhosphoricAcid),
                charge: 1,
                mz: 2.0,
                summed_intensity: 10.0,
                matched_events: 1,
                frames_with_signal: 1,
                apex_rt: Some(1.1),
                apex_intensity: 10.0,
            },
            super::PseudoMs2FragmentEvidence {
                label: "y1-H3PO4".to_string(),
                series: "y".to_string(),
                cleavage_index: 1,
                neutral_loss: Some(NeutralLossKind::PhosphoricAcid),
                charge: 1,
                mz: 300.0,
                summed_intensity: 50.0,
                matched_events: 3,
                frames_with_signal: 2,
                apex_rt: Some(1.2),
                apex_intensity: 50.0,
            },
        ];
        let cutoff = super::pseudo_ms2_neutral_loss_display_cutoff(&fragments);
        assert!((cutoff - 30.0).abs() < 1e-9);
        let canvas = SvgCanvas::new(
            0.0,
            0.0,
            400.0,
            200.0,
            CoordinateRange::new(0.0, 350.0),
            CoordinateRange::new(0.0, 1000.0),
        );
        let mut svg = String::new();
        super::draw_pseudo_ms2_sticks(&mut svg, canvas, &fragments, cutoff);

        assert!(svg.contains(">b1</text>"));
        assert!(svg.contains(">y1-H3PO4</text>"));
        assert!(!svg.contains(">b1-H3PO4</text>"));
    }

    #[test]
    fn parse_args_rejects_pseudo_ms2_without_peptide() {
        let err = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--mz".into(),
            "500".into(),
            "--pseudo-ms2".into(),
        ])
        .expect_err("pseudo-ms2 without peptide should fail");
        assert!(err.to_string().contains("--pseudo-ms2 requires --peptide"));
    }

    #[test]
    fn parse_args_rejects_pseudo_ms2_with_fragment() {
        let err = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
            "--fragment".into(),
            "b8".into(),
            "--pseudo-ms2".into(),
        ])
        .expect_err("pseudo-ms2 with a selected fragment should fail");
        assert!(err
            .to_string()
            .contains("--pseudo-ms2 aggregates all peptide fragments"));
    }

    #[test]
    fn parse_args_rejects_pseudo_ms2_for_mzml() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
            "--pseudo-ms2".into(),
        ])
        .expect_err("pseudo-ms2 mzML mode should fail");
        assert!(err
            .to_string()
            .contains("--pseudo-ms2 currently supports native Bruker .d input only"));
    }

    #[test]
    fn parse_args_accepts_outdir_and_verbosity() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--outdir".into(),
            "plots".into(),
            "--verbose".into(),
        ])
        .expect("options parse");
        assert_eq!(options.outdir.as_deref(), Some(Path::new("plots")));
        assert_eq!(options.verbosity, super::Verbosity::Verbose);
    }

    #[test]
    fn parse_args_accepts_repeatable_mz_targets_for_bruker() {
        let options = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--target".into(),
            "hexnac_204:204.0867".into(),
            "--target".into(),
            "292.1027".into(),
            "--mz-da".into(),
            "0.02".into(),
        ])
        .expect("target options parse");
        assert_eq!(options.mz_targets.len(), 2);
        assert_eq!(options.mz_targets[0].label, "hexnac_204");
        assert!((options.mz_targets[0].mz - 204.0867).abs() < 1e-9);
        assert_eq!(options.mz_targets[1].label, "mz_292.1027");
        assert_eq!(
            options.request.as_ref().expect("request").mz,
            options.mz_targets[0].mz
        );
    }

    #[test]
    fn parse_args_rejects_repeatable_targets_for_mzml() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--target".into(),
            "hexnac_204:204.0867".into(),
        ])
        .expect_err("mzML target batch should fail");
        assert!(err
            .to_string()
            .contains("repeatable --target currently supports native Bruker .d"));
    }

    #[test]
    fn parse_args_rejects_duplicate_target_labels() {
        let err = parse_args(vec![
            "--bruker".into(),
            "run.d".into(),
            "--target".into(),
            "hexnac:204.0867".into(),
            "--target".into(),
            "hexnac:204.0868".into(),
        ])
        .expect_err("duplicate target labels should fail");
        assert!(err.to_string().contains("duplicate --target label"));
    }

    #[test]
    fn multi_target_out_prefix_uses_out_prefix_as_base() {
        let input = super::DiaInput::Bruker(PathBuf::from("run.d"));
        let target = super::DiaMzTarget {
            label: "hexnac_204".to_string(),
            mz: 204.0867,
        };
        let path = super::multi_target_out_prefix(
            &input,
            Some(&PathBuf::from("plots")),
            Some(&PathBuf::from("glyco")),
            &target,
        );
        assert_eq!(path, PathBuf::from("plots/glyco__hexnac_204"));
    }

    #[test]
    fn parse_args_accepts_outdir_and_out_prefix_together() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--outdir".into(),
            "plots".into(),
            "--out-prefix".into(),
            "custom".into(),
        ])
        .expect("outdir and out-prefix should compose");
        assert_eq!(options.outdir.as_deref(), Some(Path::new("plots")));
        assert_eq!(options.out_prefix.as_deref(), Some(Path::new("custom")));
    }

    #[test]
    fn parse_args_rejects_path_like_out_prefix() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--out-prefix".into(),
            "~/Documents/".into(),
        ])
        .expect_err("out-prefix should reject directories");
        assert!(err
            .to_string()
            .contains("--out-prefix expects a file stem/name"));
    }

    #[test]
    fn precursor_rt_window_uses_strongest_fixed_width_region() {
        let evidence = super::infer_precursor_rt_window_from_points(
            &[(0.0, 1.0), (0.2, 2.0), (5.0, 10.0), (5.2, 8.0), (9.0, 7.0)],
            500.0,
            0.5,
        )
        .expect("window evidence");
        assert_eq!(evidence.rt_window.source, "precursor-inferred");
        assert!((evidence.rt_window.min - 5.0).abs() < 1e-9);
        assert!((evidence.rt_window.max - 5.5).abs() < 1e-9);
        assert_eq!(evidence.frames_with_signal, 5);
        assert!((evidence.apex_rt.expect("apex rt") - 5.0).abs() < 1e-9);
        assert!((evidence.apex_intensity - 10.0).abs() < 1e-9);
    }

    #[test]
    fn parse_args_uses_peptide_charge_suffix_unless_overridden() {
        let suffix = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "PEPTIDE/3".into(),
        ])
        .expect("suffix options")
        .request
        .expect("suffix request")
        .peptide_target
        .expect("suffix target");
        assert_eq!(suffix.charge, 3);

        let overridden = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "PEPTIDE/3".into(),
            "--charge".into(),
            "2".into(),
        ])
        .expect("override options")
        .request
        .expect("override request")
        .peptide_target
        .expect("override target");
        assert_eq!(overridden.charge, 2);
    }

    #[test]
    fn parse_args_rejects_mz_and_peptide_together() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--peptide".into(),
            "PEPTIDE".into(),
        ])
        .expect_err("conflicting target modes should fail");
        assert!(err
            .to_string()
            .contains("specify only one of --mz <center> or --peptide <SEQ>"));
    }

    #[test]
    fn parse_args_rejects_fragment_without_peptide() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--mz".into(),
            "500".into(),
            "--fragment".into(),
            "b8".into(),
        ])
        .expect_err("fragment without peptide should fail");
        assert!(err.to_string().contains("--fragment requires --peptide"));
    }

    #[test]
    fn parse_args_rejects_unknown_fragment() {
        let err = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--peptide".into(),
            "GAIIGLMVGGVVIA".into(),
            "--fragment".into(),
            "b99".into(),
        ])
        .expect_err("unknown fragment should fail");
        assert!(err.to_string().contains("unknown --fragment `b99`"));
    }

    #[test]
    fn compact_display_path_prefers_current_directory_relative_path() {
        let rendered = compact_display_path_with(
            Path::new("/work/mzio/assets/demo.mzML"),
            Some(Path::new("/work/mzio")),
            Some(Path::new("/home/alex")),
            80,
        );
        assert_eq!(rendered, "assets/demo.mzML");
    }

    #[test]
    fn compact_display_path_collapses_long_home_path() {
        let rendered = compact_display_path_with(
            Path::new(
                "/home/alex/windows/tims-ultra-0006/dshare/58087_58093/58093_1_ECL_1608_m_TMT10_prof_250ng_F01.d",
            ),
            Some(Path::new("/home/alex/amms06/mnt/e/MSPC001546")),
            Some(Path::new("/home/alex")),
            72,
        );
        assert!(rendered.starts_with("~/windows/tims-ultra-0006/.../"));
        assert!(rendered.ends_with(".d"));
        assert!(rendered.chars().count() <= 72);
    }
}
