use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs;
use std::io::{ErrorKind, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::Context;
use flate2::read::GzDecoder;

const SURVEY_TITLE_FONT: f64 = 28.0;
const SURVEY_SUBTITLE_FONT: f64 = 14.0;
const SURVEY_PANEL_TITLE_FONT: f64 = 18.0;
const SURVEY_TICK_FONT: f64 = 13.0;
const SURVEY_AXIS_LABEL_FONT: f64 = 15.0;
const SURVEY_MESSAGE_FONT: f64 = 14.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SurveyView {
    Chrom,
    Ms1Map,
    Default,
    Dda,
    All,
}

impl SurveyView {
    fn parse(value: &str) -> anyhow::Result<Self> {
        match value {
            "chrom" | "chromatogram" | "chromatograms" => Ok(Self::Chrom),
            "map" | "heatmap" | "base-peak-map" | "ms1-map" => Ok(Self::Ms1Map),
            "default" | "survey" => Ok(Self::Default),
            "dda" => Ok(Self::Dda),
            "all" => Ok(Self::All),
            other => anyhow::bail!(
                "unknown survey view `{other}`; expected chrom, ms1-map, dda, default, or all"
            ),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlotKind {
    Chrom,
    Ms1Map,
    DdaDensity,
    DdaPrecursors,
}

impl PlotKind {
    fn suffix(self) -> &'static str {
        match self {
            Self::Chrom => "survey_chrom",
            Self::Ms1Map => "survey_ms1_map",
            Self::DdaDensity => "survey_dda_density",
            Self::DdaPrecursors => "survey_dda_precursors",
        }
    }

    fn default_height(self) -> u32 {
        match self {
            Self::Chrom => 760,
            Self::Ms1Map => 760,
            Self::DdaDensity => 980,
            Self::DdaPrecursors => 980,
        }
    }

    fn title(self) -> &'static str {
        match self {
            Self::Chrom => "Run chromatograms",
            Self::Ms1Map => "MS1 base peak map",
            Self::DdaDensity => "DDA acquisition density",
            Self::DdaPrecursors => "DDA precursor survey",
        }
    }
}

#[derive(Clone, Debug)]
struct SurveyOptions {
    mzml_path: PathBuf,
    view: SurveyView,
    out_dir: PathBuf,
    prefix: Option<String>,
    svg_path: Option<PathBuf>,
    png_requested: bool,
    png_path: Option<PathBuf>,
    width: u32,
    height: Option<u32>,
    bins: usize,
    normalize: bool,
    title: Option<String>,
}

#[derive(Clone, Debug)]
struct SurveyPoint {
    ms_level: u8,
    rt: f64,
    tic: Option<f64>,
    base_peak_mz: Option<f64>,
    base_peak_intensity: Option<f64>,
    injection_time_ms: Option<f64>,
    selected_ion_mz: Option<f64>,
    selected_ion_charge: Option<u8>,
    selected_ion_intensity: Option<f64>,
}

#[derive(Clone, Debug)]
struct SurveyData {
    file_name: String,
    total_spectra: usize,
    ms1_spectra: usize,
    msn_spectra: usize,
    rt_bounds: Option<(f64, f64)>,
    ms1_mz_bounds: Option<(f64, f64)>,
    precursor_mz_bounds: Option<(f64, f64)>,
    tic_max: Option<f64>,
    base_peak_intensity_max: Option<f64>,
    selected_ion_intensity_max: Option<f64>,
    points: Vec<SurveyPoint>,
}

#[derive(Clone, Copy, Debug)]
struct Panel {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
}

pub(crate) fn run(args: Vec<String>) -> anyhow::Result<()> {
    if args
        .iter()
        .any(|arg| matches!(arg.as_str(), "-h" | "--help" | "help"))
    {
        print_help();
        return Ok(());
    }

    let options = parse_args(args)?;
    let plot_kinds = resolve_plot_kinds(options.view);
    if options.svg_path.is_some() && plot_kinds.len() != 1 {
        anyhow::bail!("--svg requires exactly one output; select a single --view");
    }
    if options.png_path.is_some() && plot_kinds.len() != 1 {
        anyhow::bail!("--png <path> requires exactly one output; select a single --view");
    }

    let data = collect_survey_data(&options.mzml_path)?;
    for kind in plot_kinds {
        let svg_path = options
            .svg_path
            .clone()
            .unwrap_or_else(|| default_svg_path(&options, kind));
        write_survey_svg(&svg_path, &data, &options, kind)
            .with_context(|| format!("failed to write {}", svg_path.display()))?;
        println!("wrote {}", svg_path.display());

        if options.png_requested {
            let png_path = options
                .png_path
                .clone()
                .unwrap_or_else(|| svg_path.with_extension("png"));
            convert_svg_to_png(&svg_path, &png_path)?;
            println!("wrote {}", png_path.display());
        }
    }

    Ok(())
}

fn resolve_plot_kinds(view: SurveyView) -> Vec<PlotKind> {
    match view {
        SurveyView::Chrom => vec![PlotKind::Chrom],
        SurveyView::Ms1Map => vec![PlotKind::Ms1Map],
        SurveyView::Default => vec![PlotKind::Chrom, PlotKind::Ms1Map],
        SurveyView::Dda => vec![PlotKind::DdaDensity, PlotKind::DdaPrecursors],
        SurveyView::All => vec![
            PlotKind::Chrom,
            PlotKind::Ms1Map,
            PlotKind::DdaDensity,
            PlotKind::DdaPrecursors,
        ],
    }
}

fn parse_args(args: Vec<String>) -> anyhow::Result<SurveyOptions> {
    let mut mzml_path = None::<PathBuf>;
    let mut view = SurveyView::Default;
    let mut out_dir = PathBuf::from("exports");
    let mut prefix = None::<String>;
    let mut svg_path = None::<PathBuf>;
    let mut png_requested = false;
    let mut png_path = None::<PathBuf>;
    let mut width = 1400_u32;
    let mut height = None::<u32>;
    let mut bins = 4000_usize;
    let mut normalize = false;
    let mut title = None::<String>;

    let mut iter = args.into_iter().peekable();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                mzml_path = Some(PathBuf::from(next_value(&mut iter, "--mzml")?));
            }
            "--view" => {
                view = SurveyView::parse(&next_value(&mut iter, "--view")?)?;
            }
            "--chrom" => view = SurveyView::Chrom,
            "--map" => view = SurveyView::Ms1Map,
            "--all" => view = SurveyView::All,
            "--out-dir" => {
                out_dir = PathBuf::from(next_value(&mut iter, "--out-dir")?);
            }
            "--prefix" => {
                prefix = Some(sanitize_filename_component(&next_value(
                    &mut iter, "--prefix",
                )?));
            }
            "--svg" => {
                svg_path = Some(PathBuf::from(next_value(&mut iter, "--svg")?));
            }
            "--png" => {
                png_requested = true;
                if let Some(next) = iter.peek() {
                    if !next.starts_with('-') {
                        png_path = Some(PathBuf::from(iter.next().expect("peeked value exists")));
                    }
                }
            }
            "--width" => {
                width = next_value(&mut iter, "--width")?
                    .parse()
                    .context("--width expects an integer pixel width")?;
            }
            "--height" => {
                height = Some(
                    next_value(&mut iter, "--height")?
                        .parse()
                        .context("--height expects an integer pixel height")?,
                );
            }
            "--bins" => {
                bins = next_value(&mut iter, "--bins")?
                    .parse()
                    .context("--bins expects an integer")?;
            }
            "--normalize" | "--norm" => normalize = true,
            "--title" => {
                title = Some(next_value(&mut iter, "--title")?);
            }
            other if other.starts_with('-') => {
                anyhow::bail!("unknown plot-survey option `{other}`")
            }
            other => {
                if mzml_path.is_some() {
                    anyhow::bail!("plot-survey accepts only one mzML input path");
                }
                mzml_path = Some(PathBuf::from(other));
            }
        }
    }

    let mzml_path =
        mzml_path.ok_or_else(|| anyhow::anyhow!("plot-survey requires an mzML path"))?;
    if !mzml_path.exists() {
        anyhow::bail!("mzML input does not exist: {}", mzml_path.display());
    }
    if width < 600 {
        anyhow::bail!("--width must be at least 600");
    }
    if matches!(height, Some(value) if value < 420) {
        anyhow::bail!("--height must be at least 420");
    }
    if bins < 32 {
        anyhow::bail!("--bins must be at least 32");
    }

    Ok(SurveyOptions {
        mzml_path,
        view,
        out_dir,
        prefix,
        svg_path,
        png_requested,
        png_path,
        width,
        height,
        bins,
        normalize,
        title,
    })
}

fn next_value(
    iter: &mut std::iter::Peekable<std::vec::IntoIter<String>>,
    flag: &str,
) -> anyhow::Result<String> {
    iter.next()
        .ok_or_else(|| anyhow::anyhow!("{flag} expects a value"))
}

fn print_help() {
    let program = crate::program_name();
    println!("{program} plot-survey");
    println!();
    println!("USAGE:");
    println!("  {program} plot-survey <file.mzML> [options]");
    println!("  {program} plot-survey --mzml <file.mzML> [options]");
    println!();
    println!("OPTIONS:");
    println!(
        "  --view <chrom|ms1-map|dda|default|all>  Select one view or view pack (default: default)"
    );
    println!("  --chrom                   Shortcut for --view chrom");
    println!("  --map                     Shortcut for --view ms1-map");
    println!("  --all                     Shortcut for --view all");
    println!("  --out-dir <dir>           Output directory for generated SVGs (default: exports)");
    println!("  --prefix <text>           Output filename prefix (default: mzML stem)");
    println!("  --svg <path>              Output SVG path; only valid for a single concrete view");
    println!(
        "  --png [path]              Also render PNGs; explicit path requires one concrete view"
    );
    println!("  --width <px>              SVG width (default: 1400)");
    println!("  --height <px>             SVG height (default depends on --view)");
    println!(
        "  --bins <n>                Max trace/map bins for compact SVG output (default: 4000)"
    );
    println!("  --normalize               Plot TIC/BPC as relative intensity");
    println!("  --title <text>            Override the SVG title");
    println!("  --help                    Show this help");
    println!();
    println!("EXAMPLES:");
    println!("  {program} plot-survey sample.mzML");
    println!("  {program} plot-survey sample.mzML --view dda --png");
    println!("  {program} plot-survey sample.mzML --chrom --svg sample_chrom.svg");
}

#[derive(Debug, Default)]
struct CurrentSpectrum {
    ms_level: Option<u8>,
    rt: Option<f64>,
    tic: Option<f64>,
    base_peak_mz: Option<f64>,
    base_peak_intensity: Option<f64>,
    injection_time_ms: Option<f64>,
    selected_ion_mz: Option<f64>,
    selected_ion_charge: Option<u8>,
    selected_ion_intensity: Option<f64>,
}

#[derive(Debug, Default)]
struct SurveyAccumulator {
    total_spectra: usize,
    ms1_spectra: usize,
    msn_spectra: usize,
    rt_min: f64,
    rt_max: f64,
    ms1_mz_min: f64,
    ms1_mz_max: f64,
    precursor_mz_min: f64,
    precursor_mz_max: f64,
    tic_max: f64,
    bpi_max: f64,
    selected_ion_intensity_max: f64,
    points: Vec<SurveyPoint>,
}

impl SurveyAccumulator {
    fn new() -> Self {
        Self {
            rt_min: f64::INFINITY,
            rt_max: -f64::INFINITY,
            ms1_mz_min: f64::INFINITY,
            ms1_mz_max: -f64::INFINITY,
            precursor_mz_min: f64::INFINITY,
            precursor_mz_max: -f64::INFINITY,
            ..Self::default()
        }
    }

    fn push(&mut self, spectrum: CurrentSpectrum) {
        self.total_spectra = self.total_spectra.saturating_add(1);
        let ms_level = spectrum.ms_level.unwrap_or(1);
        let Some(rt) = spectrum
            .rt
            .filter(|value| value.is_finite() && *value >= 0.0)
        else {
            if ms_level <= 1 {
                self.ms1_spectra = self.ms1_spectra.saturating_add(1);
            } else {
                self.msn_spectra = self.msn_spectra.saturating_add(1);
            }
            return;
        };
        self.rt_min = self.rt_min.min(rt);
        self.rt_max = self.rt_max.max(rt);

        if ms_level <= 1 {
            self.ms1_spectra = self.ms1_spectra.saturating_add(1);
            if let Some(value) = spectrum
                .tic
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                self.tic_max = self.tic_max.max(value);
            }
            if let Some(value) = spectrum
                .base_peak_intensity
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                self.bpi_max = self.bpi_max.max(value);
            }
            if let Some(value) = spectrum
                .base_peak_mz
                .filter(|value| value.is_finite() && *value > 0.0)
            {
                self.ms1_mz_min = self.ms1_mz_min.min(value);
                self.ms1_mz_max = self.ms1_mz_max.max(value);
            }
        } else {
            self.msn_spectra = self.msn_spectra.saturating_add(1);
            if let Some(value) = spectrum
                .selected_ion_mz
                .filter(|value| value.is_finite() && *value > 0.0)
            {
                self.precursor_mz_min = self.precursor_mz_min.min(value);
                self.precursor_mz_max = self.precursor_mz_max.max(value);
            }
            if let Some(value) = spectrum
                .selected_ion_intensity
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                self.selected_ion_intensity_max = self.selected_ion_intensity_max.max(value);
            }
        }

        self.points.push(SurveyPoint {
            ms_level,
            rt,
            tic: spectrum.tic,
            base_peak_mz: spectrum.base_peak_mz,
            base_peak_intensity: spectrum.base_peak_intensity,
            injection_time_ms: spectrum.injection_time_ms,
            selected_ion_mz: spectrum.selected_ion_mz,
            selected_ion_charge: spectrum.selected_ion_charge,
            selected_ion_intensity: spectrum.selected_ion_intensity,
        });
    }
}

fn collect_survey_data(path: &Path) -> anyhow::Result<SurveyData> {
    let mut reader = open_maybe_gzip(path)?;
    let mut acc = SurveyAccumulator::new();
    let mut current = None::<CurrentSpectrum>;
    let mut in_tag = false;
    let mut tag = Vec::<u8>::with_capacity(512);
    let mut buffer = [0_u8; 64 * 1024];

    loop {
        let n = reader
            .read(&mut buffer)
            .with_context(|| format!("failed to read mzML at {}", path.display()))?;
        if n == 0 {
            break;
        }
        for &byte in &buffer[..n] {
            if in_tag {
                if byte == b'>' {
                    process_xml_tag(&tag, &mut current, &mut acc);
                    tag.clear();
                    in_tag = false;
                } else if tag.len() < 16 * 1024 {
                    tag.push(byte);
                }
            } else if byte == b'<' {
                in_tag = true;
                tag.clear();
            }
        }
    }

    if let Some(spectrum) = current.take() {
        acc.push(spectrum);
    }

    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .map(ToString::to_string)
        .unwrap_or_else(|| path.to_string_lossy().to_string());

    Ok(SurveyData {
        file_name,
        total_spectra: acc.total_spectra,
        ms1_spectra: acc.ms1_spectra,
        msn_spectra: acc.msn_spectra,
        rt_bounds: finite_bounds(acc.rt_min, acc.rt_max),
        ms1_mz_bounds: finite_bounds(acc.ms1_mz_min, acc.ms1_mz_max),
        precursor_mz_bounds: finite_bounds(acc.precursor_mz_min, acc.precursor_mz_max),
        tic_max: (acc.tic_max > 0.0).then_some(acc.tic_max),
        base_peak_intensity_max: (acc.bpi_max > 0.0).then_some(acc.bpi_max),
        selected_ion_intensity_max: (acc.selected_ion_intensity_max > 0.0)
            .then_some(acc.selected_ion_intensity_max),
        points: acc.points,
    })
}

fn open_maybe_gzip(path: &Path) -> anyhow::Result<Box<dyn Read>> {
    let mut file = fs::File::open(path)
        .with_context(|| format!("failed to open mzML at {}", path.display()))?;
    let mut magic = [0_u8; 2];
    let n = file
        .read(&mut magic)
        .with_context(|| format!("failed to read mzML magic bytes at {}", path.display()))?;
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("failed to rewind mzML at {}", path.display()))?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn process_xml_tag(
    tag_bytes: &[u8],
    current: &mut Option<CurrentSpectrum>,
    acc: &mut SurveyAccumulator,
) {
    let tag = String::from_utf8_lossy(tag_bytes);
    let tag = tag.trim();
    if is_start_tag(tag, "spectrum") {
        if let Some(spectrum) = current.take() {
            acc.push(spectrum);
        }
        *current = Some(CurrentSpectrum::default());
        return;
    }
    if is_end_tag(tag, "spectrum") {
        if let Some(spectrum) = current.take() {
            acc.push(spectrum);
        }
        return;
    }

    let Some(spectrum) = current.as_mut() else {
        return;
    };
    if !is_start_tag(tag, "cvParam") {
        return;
    }

    let accession = attr_value(tag, "accession");
    let value = attr_value(tag, "value");
    let Some(accession) = accession.as_deref() else {
        return;
    };

    match accession {
        "MS:1000511" => {
            spectrum.ms_level = value.and_then(|raw| raw.parse::<u8>().ok());
        }
        "MS:1000016" => {
            spectrum.rt = value.and_then(|raw| raw.parse::<f64>().ok()).map(|rt| {
                convert_time_value(
                    rt,
                    attr_value(tag, "unitName").as_deref(),
                    attr_value(tag, "unitAccession").as_deref(),
                    TimeUnit::Minute,
                )
            });
        }
        "MS:1000927" => {
            spectrum.injection_time_ms =
                value.and_then(|raw| raw.parse::<f64>().ok()).map(|time| {
                    convert_time_value(
                        time,
                        attr_value(tag, "unitName").as_deref(),
                        attr_value(tag, "unitAccession").as_deref(),
                        TimeUnit::Millisecond,
                    )
                });
        }
        "MS:1000285" => {
            spectrum.tic = value.and_then(|raw| raw.parse::<f64>().ok());
        }
        "MS:1000504" => {
            spectrum.base_peak_mz = value.and_then(|raw| raw.parse::<f64>().ok());
        }
        "MS:1000505" => {
            spectrum.base_peak_intensity = value.and_then(|raw| raw.parse::<f64>().ok());
        }
        "MS:1000744" => {
            spectrum.selected_ion_mz = value.and_then(|raw| raw.parse::<f64>().ok());
        }
        "MS:1000041" => {
            spectrum.selected_ion_charge = value.and_then(|raw| raw.parse::<u8>().ok());
        }
        "MS:1000042" => {
            spectrum.selected_ion_intensity = value.and_then(|raw| raw.parse::<f64>().ok());
        }
        _ => {}
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TimeUnit {
    Minute,
    Millisecond,
}

fn convert_time_value(
    value: f64,
    unit_name: Option<&str>,
    unit_accession: Option<&str>,
    target: TimeUnit,
) -> f64 {
    let unit_name = unit_name.unwrap_or_default().to_ascii_lowercase();
    let unit_accession = unit_accession.unwrap_or_default();
    let seconds = if unit_name.contains("millisecond") || unit_accession == "UO:0000028" {
        value / 1000.0
    } else if unit_name.contains("second") || unit_accession == "UO:0000010" {
        value
    } else if unit_name.contains("minute") || unit_accession == "UO:0000031" {
        value * 60.0
    } else {
        match target {
            TimeUnit::Minute => value * 60.0,
            TimeUnit::Millisecond => value / 1000.0,
        }
    };

    match target {
        TimeUnit::Minute => seconds / 60.0,
        TimeUnit::Millisecond => seconds * 1000.0,
    }
}

fn is_start_tag(tag: &str, name: &str) -> bool {
    let Some(rest) = tag.strip_prefix(name) else {
        return false;
    };
    rest.is_empty()
        || rest
            .chars()
            .next()
            .is_some_and(|ch| ch.is_ascii_whitespace() || ch == '/')
}

fn is_end_tag(tag: &str, name: &str) -> bool {
    tag.strip_prefix('/')
        .is_some_and(|rest| rest == name || rest.starts_with(&format!("{name}>")))
}

fn attr_value(tag: &str, name: &str) -> Option<String> {
    let needle = format!("{name}=");
    let idx = tag.find(&needle)?;
    let rest = &tag[idx + needle.len()..];
    let mut chars = rest.chars();
    let quote = chars.next()?;
    if quote != '"' && quote != '\'' {
        return None;
    }
    let body = chars.as_str();
    let end = body.find(quote)?;
    Some(body[..end].to_string())
}

fn finite_bounds(min: f64, max: f64) -> Option<(f64, f64)> {
    if !min.is_finite() || !max.is_finite() {
        return None;
    }
    if min >= max {
        Some((min, min + 1.0))
    } else {
        Some((min, max))
    }
}

fn write_survey_svg(
    path: &Path,
    data: &SurveyData,
    options: &SurveyOptions,
    kind: PlotKind,
) -> anyhow::Result<()> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let width = options.width;
    let height = options.height.unwrap_or_else(|| kind.default_height());
    let mut file =
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?;

    writeln!(
        file,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"##
    )?;
    writeln!(
        file,
        r##"<rect width="100%" height="100%" fill="#ffffff"/>"##
    )?;
    writeln!(
        file,
        r##"<style>
text {{ font-family: Helvetica, Arial, sans-serif; fill: #17202a; }}
.muted {{ fill: #667085; }}
.axis {{ stroke: #2f3945; stroke-width: 1; fill: none; }}
.grid {{ stroke: #d7dde5; stroke-width: 1; }}
.panel-title {{ font-size: {panel_title_font:.1}px; font-weight: 700; }}
.tick {{ font-size: {tick_font:.1}px; fill: #4b5563; }}
.label {{ font-size: {axis_label_font:.1}px; fill: #344054; }}
</style>"##,
        panel_title_font = SURVEY_PANEL_TITLE_FONT,
        tick_font = SURVEY_TICK_FONT,
        axis_label_font = SURVEY_AXIS_LABEL_FONT
    )?;

    draw_header(&mut file, data, options, kind, width)?;

    let margin_left = 108.0;
    let margin_right = 42.0;
    let plot_w = width as f64 - margin_left - margin_right;
    let top = 110.0;
    let gap = 68.0;
    let bottom = 84.0;
    match kind {
        PlotKind::Chrom => {
            let available = height as f64 - top - bottom - gap;
            let panel_h = (available / 2.0).max(120.0);
            draw_tic_bpc(
                &mut file,
                data,
                options,
                Panel {
                    x: margin_left,
                    y: top,
                    w: plot_w,
                    h: panel_h,
                },
                Panel {
                    x: margin_left,
                    y: top + panel_h + gap,
                    w: plot_w,
                    h: panel_h,
                },
            )?;
        }
        PlotKind::Ms1Map => {
            draw_base_peak_map(
                &mut file,
                data,
                options,
                Panel {
                    x: margin_left,
                    y: top,
                    w: plot_w,
                    h: (height as f64 - top - bottom).max(240.0),
                },
            )?;
        }
        PlotKind::DdaDensity => {
            draw_dda_density(
                &mut file,
                data,
                options,
                Panel {
                    x: margin_left,
                    y: top,
                    w: plot_w,
                    h: (height as f64 - top - bottom).max(420.0),
                },
            )?;
        }
        PlotKind::DdaPrecursors => {
            draw_dda_precursors(
                &mut file,
                data,
                options,
                Panel {
                    x: margin_left,
                    y: top,
                    w: plot_w,
                    h: (height as f64 - top - bottom).max(420.0),
                },
            )?;
        }
    }

    writeln!(file, "</svg>")?;
    Ok(())
}

fn draw_header<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    options: &SurveyOptions,
    kind: PlotKind,
    width: u32,
) -> anyhow::Result<()> {
    let title = options
        .title
        .clone()
        .unwrap_or_else(|| format!("{}: {}", kind.title(), data.file_name));
    let rt = data
        .rt_bounds
        .map(|(lo, hi)| format!("{lo:.2}-{hi:.2} min"))
        .unwrap_or_else(|| "RT unavailable".to_string());
    let tic = data
        .tic_max
        .map(format_intensity)
        .unwrap_or_else(|| "unavailable".to_string());
    let bpc = data
        .base_peak_intensity_max
        .map(format_intensity)
        .unwrap_or_else(|| "unavailable".to_string());
    let mut subtitle = format!(
        "Source: {} | spectra {} (MS1 {}, MS2+ {}) | RT {} | TIC max {} | BPC max {}",
        data.file_name, data.total_spectra, data.ms1_spectra, data.msn_spectra, rt, tic, bpc
    );
    if matches!(kind, PlotKind::DdaDensity | PlotKind::DdaPrecursors) {
        if let Some(value) = data.selected_ion_intensity_max {
            subtitle.push_str(&format!(
                " | precursor intensity max {}",
                format_intensity(value)
            ));
        }
        if let Some(value) = median_ms2_injection_time_ms(data) {
            subtitle.push_str(&format!(" | MS2 IT median {value:.1} ms"));
        }
    }

    writeln!(
        writer,
        r##"<text x="40" y="42" font-size="{:.1}" font-weight="700">{}</text>"##,
        SURVEY_TITLE_FONT,
        escape_xml(&title)
    )?;
    writeln!(
        writer,
        r##"<text x="40" y="70" font-size="{:.1}" class="muted">{}</text>"##,
        SURVEY_SUBTITLE_FONT,
        escape_xml(&subtitle)
    )?;
    writeln!(
        writer,
        r##"<line x1="40" y1="88" x2="{x2}" y2="88" stroke="#d0d5dd" stroke-width="1"/>"##,
        x2 = width.saturating_sub(40)
    )?;
    Ok(())
}

fn median_ms2_injection_time_ms(data: &SurveyData) -> Option<f64> {
    let mut values = data
        .points
        .iter()
        .filter(|point| point.ms_level > 1)
        .filter_map(|point| point.injection_time_ms)
        .filter(|value| value.is_finite() && *value >= 0.0)
        .collect::<Vec<_>>();
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        Some((values[mid - 1] + values[mid]) / 2.0)
    } else {
        Some(values[mid])
    }
}

fn draw_tic_bpc<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    options: &SurveyOptions,
    tic_panel: Panel,
    bpc_panel: Panel,
) -> anyhow::Result<()> {
    let rt_bounds = match data.rt_bounds {
        Some(bounds) => bounds,
        None => {
            draw_empty_panel(writer, tic_panel, "TIC (MS1)", "No retention time metadata")?;
            draw_empty_panel(writer, bpc_panel, "BPC (MS1)", "No retention time metadata")?;
            return Ok(());
        }
    };

    let tic_points = downsample_trace(
        data.points
            .iter()
            .filter(|point| point.ms_level <= 1)
            .filter_map(|point| Some((point.rt, point.tic?))),
        rt_bounds,
        options.bins,
    );
    draw_trace_panel(
        writer,
        tic_panel,
        "TIC (MS1)",
        "total ion current",
        rt_bounds,
        &tic_points,
        "#248f5f",
        options.normalize,
    )?;

    let bpc_points = downsample_trace(
        data.points
            .iter()
            .filter(|point| point.ms_level <= 1)
            .filter_map(|point| Some((point.rt, point.base_peak_intensity?))),
        rt_bounds,
        options.bins,
    );
    draw_trace_panel(
        writer,
        bpc_panel,
        "BPC (MS1)",
        "base peak intensity",
        rt_bounds,
        &bpc_points,
        "#0089a7",
        options.normalize,
    )?;

    Ok(())
}

fn draw_trace_panel<W: Write>(
    writer: &mut W,
    panel: Panel,
    title: &str,
    y_label: &str,
    rt_bounds: (f64, f64),
    points: &[(f64, f64)],
    color: &str,
    normalize: bool,
) -> anyhow::Result<()> {
    draw_panel_axes(
        writer,
        panel,
        title,
        "RT minutes",
        y_label,
        rt_bounds,
        (0.0, 1.0),
    )?;

    if points.is_empty() {
        draw_panel_message(writer, panel, "No metadata values found")?;
        return Ok(());
    }

    let raw_y_max = points
        .iter()
        .fold(0.0_f64, |acc, (_, y)| acc.max(*y))
        .max(1e-12);
    let y_max = if normalize { 1.05 } else { raw_y_max * 1.08 };
    draw_y_grid(writer, panel, (0.0, y_max), normalize)?;
    draw_x_ticks(writer, panel, rt_bounds)?;

    let x_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let mut d = String::new();
    for (idx, &(rt, y)) in points.iter().enumerate() {
        let y = if normalize { y / raw_y_max } else { y };
        let px = panel.x + ((rt - rt_bounds.0) / x_span).clamp(0.0, 1.0) * panel.w;
        let py = panel.y + (1.0 - (y / y_max).clamp(0.0, 1.0)) * panel.h;
        if idx == 0 {
            d.push_str(&format!("M{px:.2},{py:.2}"));
        } else {
            d.push_str(&format!(" L{px:.2},{py:.2}"));
        }
    }
    writeln!(
        writer,
        r##"<path d="{d}" fill="none" stroke="{color}" stroke-width="1.8" stroke-linejoin="round" stroke-linecap="round"/>"##
    )?;
    Ok(())
}

fn draw_base_peak_map<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    options: &SurveyOptions,
    panel: Panel,
) -> anyhow::Result<()> {
    let Some(rt_bounds) = data.rt_bounds else {
        draw_empty_panel(
            writer,
            panel,
            "MS1 Base Peak Map",
            "No retention time metadata",
        )?;
        return Ok(());
    };
    let Some(mz_bounds) = data.ms1_mz_bounds else {
        draw_empty_panel(
            writer,
            panel,
            "MS1 Base Peak Map",
            "No base peak m/z metadata",
        )?;
        return Ok(());
    };
    draw_panel_axes(
        writer,
        panel,
        "MS1 Base Peak Map",
        "RT minutes",
        "base peak m/z",
        rt_bounds,
        mz_bounds,
    )?;
    draw_y_grid(writer, panel, mz_bounds, false)?;
    draw_x_ticks(writer, panel, rt_bounds)?;

    let map_points = data
        .points
        .iter()
        .filter(|point| point.ms_level <= 1)
        .filter_map(|point| Some((point.rt, point.base_peak_mz?, point.base_peak_intensity?)))
        .filter(|(_, _, intensity)| intensity.is_finite() && *intensity > 0.0)
        .collect::<Vec<_>>();
    if map_points.is_empty() {
        draw_panel_message(writer, panel, "No base peak intensity metadata found")?;
        return Ok(());
    }

    let x_bins = options.bins.clamp(48, 1800);
    let y_bins = ((x_bins as f64 * panel.h / panel.w).round() as usize).clamp(32, 360);
    let cell_w = panel.w / x_bins as f64;
    let cell_h = panel.h / y_bins as f64;
    let x_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let y_span = (mz_bounds.1 - mz_bounds.0).max(1e-9);

    let mut cells = vec![f64::NEG_INFINITY; x_bins * y_bins];
    for (rt, mz, intensity) in map_points {
        let mut xi = (((rt - rt_bounds.0) / x_span).clamp(0.0, 1.0) * x_bins as f64) as usize;
        let mut yi = (((mz - mz_bounds.0) / y_span).clamp(0.0, 1.0) * y_bins as f64) as usize;
        if xi >= x_bins {
            xi = x_bins - 1;
        }
        if yi >= y_bins {
            yi = y_bins - 1;
        }
        let value = intensity.log10();
        let idx = yi * x_bins + xi;
        cells[idx] = cells[idx].max(value);
    }

    let filled = cells
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    if filled.is_empty() {
        draw_panel_message(writer, panel, "No finite base peak intensities found")?;
        return Ok(());
    }
    let min_v = filled
        .iter()
        .copied()
        .fold(f64::INFINITY, |acc, value| acc.min(value));
    let max_v = filled
        .iter()
        .copied()
        .fold(-f64::INFINITY, |acc, value| acc.max(value));
    let span = (max_v - min_v).max(1e-9);

    for yi in 0..y_bins {
        for xi in 0..x_bins {
            let value = cells[yi * x_bins + xi];
            if !value.is_finite() {
                continue;
            }
            let frac = ((value - min_v) / span).clamp(0.0, 1.0);
            let fill = viridis(frac);
            let x = panel.x + xi as f64 * cell_w;
            let y = panel.y + panel.h - (yi + 1) as f64 * cell_h;
            writeln!(
                writer,
                r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}" opacity="0.92"/>"##,
                w = cell_w.max(0.6),
                h = cell_h.max(0.6),
            )?;
        }
    }

    draw_map_legend(writer, panel, min_v, max_v)?;
    Ok(())
}

fn draw_dda_density<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    options: &SurveyOptions,
    area: Panel,
) -> anyhow::Result<()> {
    let gap = 62.0;
    let panel_h = ((area.h - gap * 3.0) / 4.0).max(95.0);
    let tic_panel = Panel {
        x: area.x,
        y: area.y,
        w: area.w,
        h: panel_h,
    };
    let bpc_panel = Panel {
        x: area.x,
        y: area.y + panel_h + gap,
        w: area.w,
        h: panel_h,
    };
    let density_panel = Panel {
        x: area.x,
        y: area.y + (panel_h + gap) * 2.0,
        w: area.w,
        h: panel_h,
    };
    let cycle_panel = Panel {
        x: area.x,
        y: area.y + (panel_h + gap) * 3.0,
        w: area.w,
        h: panel_h,
    };

    draw_tic_bpc(writer, data, options, tic_panel, bpc_panel)?;

    let Some(rt_bounds) = data.rt_bounds else {
        draw_empty_panel(
            writer,
            density_panel,
            "MS2 scans per RT bin",
            "No retention time metadata",
        )?;
        draw_empty_panel(
            writer,
            cycle_panel,
            "Cycle time / MS2 per MS1 cycle",
            "No retention time metadata",
        )?;
        return Ok(());
    };

    let density = ms2_density_points(data, rt_bounds, options.bins);
    draw_bar_trace_panel(
        writer,
        density_panel,
        "MS2 scans per RT bin",
        "MS2 scans/bin",
        rt_bounds,
        &density,
        "#7659a6",
    )?;

    let (cycle_time, ms2_counts) = cycle_summary_points(data);
    let series = [
        TraceSeries {
            label: "cycle time (s)",
            points: &cycle_time,
            color: "#2f6fb0",
        },
        TraceSeries {
            label: "MS2 per cycle",
            points: &ms2_counts,
            color: "#c65f2e",
        },
    ];
    draw_multi_trace_panel(
        writer,
        cycle_panel,
        "Cycle time / MS2 per MS1 cycle",
        "seconds / scans",
        rt_bounds,
        &series,
    )?;

    Ok(())
}

fn draw_dda_precursors<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    options: &SurveyOptions,
    area: Panel,
) -> anyhow::Result<()> {
    let gap = 78.0;
    let map_h = (area.h * 0.72).max(300.0);
    let charge_h = (area.h - map_h - gap).max(150.0);
    let map_panel = Panel {
        x: area.x,
        y: area.y,
        w: area.w,
        h: map_h,
    };
    let charge_panel = Panel {
        x: area.x,
        y: area.y + map_h + gap,
        w: area.w,
        h: charge_h,
    };

    draw_precursor_map(writer, data, options, map_panel)?;
    draw_charge_distribution(writer, data, charge_panel)?;
    Ok(())
}

fn draw_bar_trace_panel<W: Write>(
    writer: &mut W,
    panel: Panel,
    title: &str,
    y_label: &str,
    rt_bounds: (f64, f64),
    points: &[(f64, f64)],
    color: &str,
) -> anyhow::Result<()> {
    draw_panel_axes(
        writer,
        panel,
        title,
        "RT minutes",
        y_label,
        rt_bounds,
        (0.0, 1.0),
    )?;
    if points.is_empty() {
        draw_panel_message(writer, panel, "No MS2 scans found")?;
        return Ok(());
    }
    let y_max = points
        .iter()
        .fold(0.0_f64, |acc, (_, y)| acc.max(*y))
        .max(1.0)
        * 1.12;
    draw_y_grid(writer, panel, (0.0, y_max), false)?;
    draw_x_ticks(writer, panel, rt_bounds)?;

    let x_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let bar_w = (panel.w / points.len().max(1) as f64 * 0.82).max(1.0);
    for &(rt, value) in points {
        let x = panel.x + ((rt - rt_bounds.0) / x_span).clamp(0.0, 1.0) * panel.w - bar_w / 2.0;
        let h = (value / y_max).clamp(0.0, 1.0) * panel.h;
        let y = panel.y + panel.h - h;
        writeln!(
            writer,
            r##"<rect x="{x:.2}" y="{y:.2}" width="{bar_w:.2}" height="{h:.2}" fill="{color}" opacity="0.72"/>"##
        )?;
    }
    Ok(())
}

struct TraceSeries<'a> {
    label: &'a str,
    points: &'a [(f64, f64)],
    color: &'a str,
}

fn draw_multi_trace_panel<W: Write>(
    writer: &mut W,
    panel: Panel,
    title: &str,
    y_label: &str,
    rt_bounds: (f64, f64),
    series: &[TraceSeries<'_>],
) -> anyhow::Result<()> {
    draw_panel_axes(
        writer,
        panel,
        title,
        "RT minutes",
        y_label,
        rt_bounds,
        (0.0, 1.0),
    )?;

    let y_max = series
        .iter()
        .flat_map(|entry| entry.points.iter().map(|(_, y)| *y))
        .fold(0.0_f64, |acc, y| acc.max(y))
        .max(1.0)
        * 1.12;
    if !y_max.is_finite() || series.iter().all(|entry| entry.points.is_empty()) {
        draw_panel_message(writer, panel, "No MS1 cycle structure found")?;
        return Ok(());
    }

    draw_y_grid(writer, panel, (0.0, y_max), false)?;
    draw_x_ticks(writer, panel, rt_bounds)?;
    for entry in series {
        draw_trace_path(
            writer,
            panel,
            rt_bounds,
            y_max,
            entry.points,
            entry.color,
            1.9,
        )?;
    }
    draw_inline_legend(
        writer,
        panel,
        &series
            .iter()
            .map(|entry| (entry.label, entry.color))
            .collect::<Vec<_>>(),
    )?;
    Ok(())
}

fn draw_trace_path<W: Write>(
    writer: &mut W,
    panel: Panel,
    rt_bounds: (f64, f64),
    y_max: f64,
    points: &[(f64, f64)],
    color: &str,
    stroke_width: f64,
) -> anyhow::Result<()> {
    if points.is_empty() {
        return Ok(());
    }
    let x_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let mut d = String::new();
    for (idx, &(rt, y)) in points.iter().enumerate() {
        let px = panel.x + ((rt - rt_bounds.0) / x_span).clamp(0.0, 1.0) * panel.w;
        let py = panel.y + (1.0 - (y / y_max).clamp(0.0, 1.0)) * panel.h;
        if idx == 0 {
            d.push_str(&format!("M{px:.2},{py:.2}"));
        } else {
            d.push_str(&format!(" L{px:.2},{py:.2}"));
        }
    }
    writeln!(
        writer,
        r##"<path d="{d}" fill="none" stroke="{color}" stroke-width="{stroke_width:.2}" stroke-linejoin="round" stroke-linecap="round"/>"##
    )?;
    Ok(())
}

fn draw_precursor_map<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    _options: &SurveyOptions,
    panel: Panel,
) -> anyhow::Result<()> {
    let Some(rt_bounds) = data.rt_bounds else {
        draw_empty_panel(
            writer,
            panel,
            "MS2 precursor m/z map",
            "No retention time metadata",
        )?;
        return Ok(());
    };
    let Some(mz_bounds) = data.precursor_mz_bounds else {
        draw_empty_panel(
            writer,
            panel,
            "MS2 precursor m/z map",
            "No selected-ion m/z metadata",
        )?;
        return Ok(());
    };

    draw_panel_axes(
        writer,
        panel,
        "MS2 precursor m/z map",
        "RT minutes",
        "MS2 precursor m/z",
        rt_bounds,
        mz_bounds,
    )?;
    draw_y_grid(writer, panel, mz_bounds, false)?;
    draw_x_ticks(writer, panel, rt_bounds)?;

    let mut points = data
        .points
        .iter()
        .filter(|point| point.ms_level > 1)
        .filter_map(|point| {
            Some((
                point.rt,
                point.selected_ion_mz?,
                point.selected_ion_charge,
                point.selected_ion_intensity,
            ))
        })
        .filter(|(_, mz, _, _)| mz.is_finite() && *mz > 0.0)
        .collect::<Vec<_>>();
    if points.is_empty() {
        draw_panel_message(writer, panel, "No MS2 precursor metadata found")?;
        return Ok(());
    }

    points.sort_by(|left, right| {
        let l = left.3.unwrap_or(0.0);
        let r = right.3.unwrap_or(0.0);
        l.partial_cmp(&r).unwrap_or(Ordering::Equal)
    });

    let logs = points
        .iter()
        .filter_map(|(_, _, _, intensity)| intensity.filter(|value| *value > 0.0))
        .map(|value| value.log10())
        .collect::<Vec<_>>();
    let min_log = logs
        .iter()
        .copied()
        .fold(f64::INFINITY, |acc, value| acc.min(value));
    let max_log = logs
        .iter()
        .copied()
        .fold(-f64::INFINITY, |acc, value| acc.max(value));
    let log_span = (max_log - min_log).max(1e-9);
    let x_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let y_span = (mz_bounds.1 - mz_bounds.0).max(1e-9);

    for (rt, mz, charge, intensity) in points {
        let frac = intensity
            .filter(|value| *value > 0.0)
            .map(|value| ((value.log10() - min_log) / log_span).clamp(0.0, 1.0))
            .unwrap_or(0.45);
        let px = panel.x + ((rt - rt_bounds.0) / x_span).clamp(0.0, 1.0) * panel.w;
        let py = panel.y + (1.0 - ((mz - mz_bounds.0) / y_span).clamp(0.0, 1.0)) * panel.h;
        let radius = 1.25 + 2.35 * frac;
        let opacity = 0.28 + 0.58 * frac;
        let fill = charge_color(charge);
        writeln!(
            writer,
            r##"<circle cx="{px:.2}" cy="{py:.2}" r="{radius:.2}" fill="{fill}" opacity="{opacity:.3}"/>"##
        )?;
    }
    draw_inline_legend(
        writer,
        panel,
        &[
            ("charge 2", charge_color(Some(2))),
            ("charge 3", charge_color(Some(3))),
            ("charge 4+", charge_color(Some(4))),
            ("unknown", charge_color(None)),
        ],
    )?;
    Ok(())
}

fn draw_charge_distribution<W: Write>(
    writer: &mut W,
    data: &SurveyData,
    panel: Panel,
) -> anyhow::Result<()> {
    let mut counts = BTreeMap::<u8, usize>::new();
    let mut unknown = 0_usize;
    for point in data.points.iter().filter(|point| point.ms_level > 1) {
        match point.selected_ion_charge {
            Some(charge) if charge > 0 => {
                *counts.entry(charge).or_insert(0) += 1;
            }
            _ => unknown = unknown.saturating_add(1),
        }
    }

    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="panel-title">MS2 charge distribution</text>"##,
        x = panel.x,
        y = panel.y - 14.0
    )?;
    writeln!(
        writer,
        r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="#fbfcfe" stroke="#2f3945" stroke-width="1"/>"##,
        x = panel.x,
        y = panel.y,
        w = panel.w,
        h = panel.h
    )?;

    let mut entries = counts
        .into_iter()
        .map(|(charge, count)| (format!("{charge}+"), count, charge_color(Some(charge))))
        .collect::<Vec<_>>();
    if unknown > 0 {
        entries.push(("unknown".to_string(), unknown, charge_color(None)));
    }
    if entries.is_empty() {
        draw_panel_message(writer, panel, "No MS2 charge metadata found")?;
        return Ok(());
    }

    let y_max = entries
        .iter()
        .fold(0_usize, |acc, (_, count, _)| acc.max(*count))
        .max(1) as f64
        * 1.15;
    draw_y_grid(writer, panel, (0.0, y_max), false)?;
    let slot_w = panel.w / entries.len() as f64;
    let bar_w = (slot_w * 0.62).max(4.0);
    for (idx, (label, count, color)) in entries.iter().enumerate() {
        let x = panel.x + idx as f64 * slot_w + (slot_w - bar_w) / 2.0;
        let h = (*count as f64 / y_max).clamp(0.0, 1.0) * panel.h;
        let y = panel.y + panel.h - h;
        writeln!(
            writer,
            r##"<rect x="{x:.2}" y="{y:.2}" width="{bar_w:.2}" height="{h:.2}" fill="{color}" opacity="0.82"/>"##
        )?;
        writeln!(
            writer,
            r##"<text x="{tx:.2}" y="{ty:.2}" class="tick" text-anchor="middle">{}</text>"##,
            escape_xml(label),
            tx = panel.x + idx as f64 * slot_w + slot_w / 2.0,
            ty = panel.y + panel.h + 26.0
        )?;
    }
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="label" text-anchor="middle">charge</text>"##,
        x = panel.x + panel.w / 2.0,
        y = panel.y + panel.h + 58.0
    )?;
    Ok(())
}

fn draw_inline_legend<W: Write>(
    writer: &mut W,
    panel: Panel,
    entries: &[(&str, &str)],
) -> anyhow::Result<()> {
    let mut x = panel.x + panel.w - 12.0;
    let y = panel.y - 15.0;
    for (label, color) in entries.iter().rev() {
        let text_w = label.len() as f64 * 7.2 + 28.0;
        x -= text_w;
        writeln!(
            writer,
            r##"<line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="{color}" stroke-width="3"/>"##,
            x1 = x,
            x2 = x + 15.0,
            y = y - 4.0
        )?;
        writeln!(
            writer,
            r##"<text x="{tx:.2}" y="{ty:.2}" class="tick">{}</text>"##,
            escape_xml(label),
            tx = x + 20.0,
            ty = y
        )?;
    }
    Ok(())
}

fn ms2_density_points(
    data: &SurveyData,
    rt_bounds: (f64, f64),
    bins_hint: usize,
) -> Vec<(f64, f64)> {
    let bins = bins_hint.clamp(32, 320);
    let span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let mut counts = vec![0_u32; bins];
    for point in data.points.iter().filter(|point| point.ms_level > 1) {
        let mut bin = (((point.rt - rt_bounds.0) / span).clamp(0.0, 1.0) * bins as f64) as usize;
        if bin >= bins {
            bin = bins - 1;
        }
        counts[bin] = counts[bin].saturating_add(1);
    }
    counts
        .iter()
        .enumerate()
        .map(|(idx, count)| {
            let rt = rt_bounds.0 + (idx as f64 + 0.5) * span / bins as f64;
            (rt, *count as f64)
        })
        .collect()
}

fn cycle_summary_points(data: &SurveyData) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let mut points = data
        .points
        .iter()
        .map(|point| (point.rt, point.ms_level))
        .collect::<Vec<_>>();
    points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    let mut current_ms1_rt = None::<f64>;
    let mut current_ms2_count = 0_u32;
    let mut cycle_time = Vec::<(f64, f64)>::new();
    let mut ms2_counts = Vec::<(f64, f64)>::new();

    for (rt, ms_level) in points {
        if ms_level <= 1 {
            if let Some(prev_rt) = current_ms1_rt {
                if rt > prev_rt {
                    cycle_time.push((prev_rt, (rt - prev_rt) * 60.0));
                    ms2_counts.push((prev_rt, current_ms2_count as f64));
                }
            }
            current_ms1_rt = Some(rt);
            current_ms2_count = 0;
        } else if current_ms1_rt.is_some() {
            current_ms2_count = current_ms2_count.saturating_add(1);
        }
    }

    (cycle_time, ms2_counts)
}

fn charge_color(charge: Option<u8>) -> &'static str {
    match charge {
        Some(1) => "#64748b",
        Some(2) => "#2f6fb0",
        Some(3) => "#248f5f",
        Some(4) => "#c65f2e",
        Some(5) => "#8b5fbf",
        Some(_) => "#b7791f",
        None => "#8a94a6",
    }
}

fn draw_panel_axes<W: Write>(
    writer: &mut W,
    panel: Panel,
    title: &str,
    x_label: &str,
    y_label: &str,
    _x_bounds: (f64, f64),
    _y_bounds: (f64, f64),
) -> anyhow::Result<()> {
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="panel-title">{}</text>"##,
        escape_xml(title),
        x = panel.x,
        y = panel.y - 14.0
    )?;
    writeln!(
        writer,
        r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="#fbfcfe" stroke="#2f3945" stroke-width="1"/>"##,
        x = panel.x,
        y = panel.y,
        w = panel.w,
        h = panel.h
    )?;
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="label" text-anchor="middle">{}</text>"##,
        escape_xml(x_label),
        x = panel.x + panel.w / 2.0,
        y = panel.y + panel.h + 58.0
    )?;
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="label" transform="rotate(-90 {x:.2} {y:.2})" text-anchor="middle">{}</text>"##,
        escape_xml(y_label),
        x = panel.x - 78.0,
        y = panel.y + panel.h / 2.0
    )?;
    Ok(())
}

fn draw_y_grid<W: Write>(
    writer: &mut W,
    panel: Panel,
    y_bounds: (f64, f64),
    relative: bool,
) -> anyhow::Result<()> {
    for tick in linear_ticks(y_bounds.0, y_bounds.1, 5) {
        let frac = ((tick - y_bounds.0) / (y_bounds.1 - y_bounds.0).max(1e-9)).clamp(0.0, 1.0);
        let y = panel.y + (1.0 - frac) * panel.h;
        writeln!(
            writer,
            r##"<line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" class="grid"/>"##,
            x1 = panel.x,
            x2 = panel.x + panel.w
        )?;
        let label = if relative {
            format!("{:.0}%", tick * 100.0)
        } else {
            format_axis_value(tick)
        };
        writeln!(
            writer,
            r##"<text x="{x:.2}" y="{y:.2}" class="tick" text-anchor="end" dominant-baseline="middle">{}</text>"##,
            escape_xml(&label),
            x = panel.x - 10.0,
            y = y
        )?;
    }
    Ok(())
}

fn draw_x_ticks<W: Write>(
    writer: &mut W,
    panel: Panel,
    x_bounds: (f64, f64),
) -> anyhow::Result<()> {
    for tick in linear_ticks(x_bounds.0, x_bounds.1, 6) {
        let frac = ((tick - x_bounds.0) / (x_bounds.1 - x_bounds.0).max(1e-9)).clamp(0.0, 1.0);
        let x = panel.x + frac * panel.w;
        writeln!(
            writer,
            r##"<line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="#2f3945" stroke-width="1"/>"##,
            y1 = panel.y + panel.h,
            y2 = panel.y + panel.h + 5.0
        )?;
        writeln!(
            writer,
            r##"<text x="{x:.2}" y="{y:.2}" class="tick" text-anchor="middle">{:.2}</text>"##,
            tick,
            y = panel.y + panel.h + 26.0
        )?;
    }
    Ok(())
}

fn draw_empty_panel<W: Write>(
    writer: &mut W,
    panel: Panel,
    title: &str,
    message: &str,
) -> anyhow::Result<()> {
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="panel-title">{}</text>"##,
        escape_xml(title),
        x = panel.x,
        y = panel.y - 14.0
    )?;
    writeln!(
        writer,
        r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="#fbfcfe" stroke="#2f3945" stroke-width="1"/>"##,
        x = panel.x,
        y = panel.y,
        w = panel.w,
        h = panel.h
    )?;
    draw_panel_message(writer, panel, message)?;
    Ok(())
}

fn draw_panel_message<W: Write>(writer: &mut W, panel: Panel, message: &str) -> anyhow::Result<()> {
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="muted" font-size="{:.1}" text-anchor="middle" dominant-baseline="middle">{}</text>"##,
        SURVEY_MESSAGE_FONT,
        escape_xml(message),
        x = panel.x + panel.w / 2.0,
        y = panel.y + panel.h / 2.0
    )?;
    Ok(())
}

fn draw_map_legend<W: Write>(
    writer: &mut W,
    panel: Panel,
    min_v: f64,
    max_v: f64,
) -> anyhow::Result<()> {
    let legend_w = 180.0;
    let legend_h = 10.0;
    let x = panel.x + panel.w - legend_w;
    let y = panel.y - 28.0;
    let steps = 36;
    for i in 0..steps {
        let frac0 = i as f64 / steps as f64;
        let fill = viridis(frac0);
        writeln!(
            writer,
            r##"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}"/>"##,
            x = x + frac0 * legend_w,
            y = y,
            w = legend_w / steps as f64 + 0.4,
            h = legend_h
        )?;
    }
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="tick" text-anchor="start">log10 BPC {:.1}</text>"##,
        min_v,
        x = x,
        y = y - 4.0
    )?;
    writeln!(
        writer,
        r##"<text x="{x:.2}" y="{y:.2}" class="tick" text-anchor="end">{:.1}</text>"##,
        max_v,
        x = x + legend_w,
        y = y - 4.0
    )?;
    Ok(())
}

fn downsample_trace<I>(samples: I, x_bounds: (f64, f64), bins: usize) -> Vec<(f64, f64)>
where
    I: IntoIterator<Item = (f64, f64)>,
{
    let mut raw = samples
        .into_iter()
        .filter(|(x, y)| x.is_finite() && y.is_finite() && *y >= 0.0)
        .collect::<Vec<_>>();
    raw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    if raw.len() <= bins {
        return raw;
    }

    let x_min = x_bounds.0.min(x_bounds.1);
    let x_max = x_bounds.0.max(x_bounds.1);
    let span = (x_max - x_min).max(1e-9);
    let bins = bins.clamp(32, 100_000);
    let mut best_y = vec![f64::NEG_INFINITY; bins];
    let mut best_x = vec![0.0_f64; bins];
    let mut has = vec![false; bins];

    for (x, y) in raw {
        let mut bin = (((x - x_min) / span).clamp(0.0, 1.0) * bins as f64) as usize;
        if bin >= bins {
            bin = bins - 1;
        }
        if !has[bin] || y > best_y[bin] {
            has[bin] = true;
            best_y[bin] = y;
            best_x[bin] = x;
        }
    }

    has.iter()
        .enumerate()
        .filter_map(|(idx, present)| present.then_some((best_x[idx], best_y[idx])))
        .collect()
}

fn linear_ticks(min: f64, max: f64, n: usize) -> Vec<f64> {
    if n <= 1 || !min.is_finite() || !max.is_finite() {
        return Vec::new();
    }
    let span = (max - min).max(1e-9);
    (0..n)
        .map(|idx| min + span * idx as f64 / (n - 1) as f64)
        .collect()
}

fn viridis(frac: f64) -> String {
    const STOPS: &[(f64, (u8, u8, u8))] = &[
        (0.0, (68, 1, 84)),
        (0.25, (59, 82, 139)),
        (0.5, (33, 145, 140)),
        (0.75, (94, 201, 98)),
        (1.0, (253, 231, 37)),
    ];
    let frac = frac.clamp(0.0, 1.0);
    for window in STOPS.windows(2) {
        let (left_f, left) = window[0];
        let (right_f, right) = window[1];
        if frac <= right_f {
            let local = ((frac - left_f) / (right_f - left_f).max(1e-9)).clamp(0.0, 1.0);
            let r = lerp_u8(left.0, right.0, local);
            let g = lerp_u8(left.1, right.1, local);
            let b = lerp_u8(left.2, right.2, local);
            return format!("#{r:02x}{g:02x}{b:02x}");
        }
    }
    "#fde725".to_string()
}

fn lerp_u8(a: u8, b: u8, frac: f64) -> u8 {
    (a as f64 + (b as f64 - a as f64) * frac).round() as u8
}

fn default_svg_path(options: &SurveyOptions, kind: PlotKind) -> PathBuf {
    let stem = options
        .prefix
        .clone()
        .unwrap_or_else(|| default_output_prefix(&options.mzml_path));
    options
        .out_dir
        .join(format!("{stem}_{}.svg", kind.suffix()))
}

fn default_output_prefix(input: &Path) -> String {
    input
        .file_stem()
        .and_then(|value| value.to_str())
        .map(sanitize_filename_component)
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "run".to_string())
}

fn convert_svg_to_png(svg_path: &Path, png_path: &Path) -> anyhow::Result<()> {
    if let Some(parent) = png_path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let inkscape = Command::new("inkscape")
        .arg(svg_path)
        .arg("--export-type=png")
        .arg("--export-filename")
        .arg(png_path)
        .arg("--export-area-page")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    match inkscape {
        Ok(status) if status.success() => return Ok(()),
        Ok(_) => {}
        Err(err) if err.kind() == ErrorKind::NotFound => {}
        Err(err) => return Err(err).context("failed to run inkscape for PNG export"),
    }

    let convert = Command::new("convert")
        .arg(svg_path)
        .arg(png_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    match convert {
        Ok(status) if status.success() => Ok(()),
        Ok(_) => Err(anyhow::anyhow!(
            "PNG conversion failed via inkscape/convert (SVG saved at {})",
            svg_path.display()
        )),
        Err(err) if err.kind() == ErrorKind::NotFound => Err(anyhow::anyhow!(
            "no SVG-to-PNG converter found (tried inkscape/convert); SVG saved at {}",
            svg_path.display()
        )),
        Err(err) => Err(err).context("failed to run convert for PNG export"),
    }
}

fn format_intensity(value: f64) -> String {
    if value == 0.0 {
        "0".to_string()
    } else {
        format!("{value:.3e}")
    }
}

fn format_axis_value(value: f64) -> String {
    let abs = value.abs();
    if abs > 0.0 && !(0.01..10000.0).contains(&abs) {
        format!("{value:.1e}")
    } else if abs >= 100.0 {
        format!("{value:.0}")
    } else if abs >= 10.0 {
        format!("{value:.1}")
    } else {
        format!("{value:.2}")
    }
}

fn sanitize_filename_component(input: &str) -> String {
    let mut out = String::with_capacity(input.len().min(90));
    for ch in input.chars().take(90) {
        match ch {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' | '.' => out.push(ch),
            _ => out.push('_'),
        }
    }
    out.trim_matches('_').to_string()
}

fn escape_xml(input: &str) -> String {
    input
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('\"', "&quot;")
        .replace('\'', "&apos;")
}
