use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Context;
use mzdata::spectrum::SignalContinuity;

use crate::mzml::{
    extract_scan_number, load_selected_spectrum, open_reader, LoadedSpectrum, SpectrumSelector,
};

const PROTON_MASS: f64 = 1.007_276_466_812;

#[derive(Clone, Debug, Default)]
struct ScanOptions {
    mzml_path: Option<PathBuf>,
    selector: Option<SpectrumSelector>,
    output_path: Option<PathBuf>,
    output_format: Option<ScanOutputFormat>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ScanOutputFormat {
    Tsv,
    Ms2,
}

pub fn run(command_name: &str, args: Vec<String>) -> anyhow::Result<()> {
    if args
        .iter()
        .any(|arg| matches!(arg.as_str(), "-h" | "--help" | "help"))
    {
        print_help(command_name);
        return Ok(());
    }

    let options = parse_args(args)?;
    let output_format = resolve_output_format(&options)?;
    let mzml_path = options
        .mzml_path
        .as_ref()
        .expect("parse_args validates mzML path");
    let selector = options
        .selector
        .as_ref()
        .expect("parse_args validates selector");

    let mut reader = open_reader(mzml_path)?;
    let spectrum = load_selected_spectrum(&mut reader, selector)?;
    let rendered = render_scan_export(mzml_path, &spectrum, output_format);

    if let Some(output_path) = options.output_path {
        if let Some(parent) = output_path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }
        }
        fs::write(&output_path, rendered)
            .with_context(|| format!("failed to write {}", output_path.display()))?;
        println!(
            "Wrote {} export: {}",
            output_format.label(),
            output_path.display()
        );
    } else {
        print!("{rendered}");
    }

    Ok(())
}

fn print_help(command_name: &str) {
    let program = crate::program_name();
    println!("{program} {command_name}");
    println!();
    if command_name == "scan" {
        println!("Alias of `{program} extract`.");
        println!();
    }
    println!("USAGE:");
    println!("  {program} {command_name} --mzml <file> --scan <n> [options]");
    println!();
    println!("OPTIONS:");
    println!("  --mzml <file>            Input mzML file");
    println!("  --scan <n>               Scan number, e.g. 4821 or scan=4821");
    println!("  --output <file>          Write the export to a file instead of stdout");
    println!("  --format <tsv|ms2>       Force the output format");
    println!("  --tsv <file>             Write tabular output");
    println!("  --ms2 <file>             Write MS2 output");
    println!("  --help                   Show this help");
    println!();
    println!("FORMAT:");
    println!("  `.tsv` or `.ms2` is inferred from `--output` when possible; otherwise TSV is the default.");
    println!("  TSV emits `# key<TAB>value` metadata lines followed by `mz<TAB>intensity`.");
    println!(
        "  MS2 emits a minimal single-spectrum record with `S`, `I`, and `Z` lines when available."
    );
    println!();
    println!("EXAMPLES:");
    println!("  {program} {command_name} --mzml sample.mzML --scan 4821");
    println!("  {program} extract --mzml sample.mzML --scan 4821 --output scan4821.tsv");
    println!("  {program} extract --mzml sample.mzML --scan 4821 --output scan4821.ms2");
}

fn parse_args(args: Vec<String>) -> anyhow::Result<ScanOptions> {
    let mut options = ScanOptions::default();
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                options.mzml_path =
                    Some(PathBuf::from(iter.next().context("--mzml expects a path")?));
            }
            "--scan" => {
                let raw = iter.next().context("--scan expects a scan number")?;
                let scan_number = parse_scan_number_arg(&raw)?;
                set_selector(&mut options, SpectrumSelector::ScanNumber(scan_number))?;
            }
            "--output" => {
                let output_path =
                    PathBuf::from(iter.next().context("--output expects a file path")?);
                set_output_path(&mut options, output_path)?;
            }
            "--tsv" => {
                let output_path = PathBuf::from(iter.next().context("--tsv expects a file path")?);
                set_output_path(&mut options, output_path)?;
                set_output_format(&mut options, ScanOutputFormat::Tsv)?;
            }
            "--ms2" => {
                let output_path = PathBuf::from(iter.next().context("--ms2 expects a file path")?);
                set_output_path(&mut options, output_path)?;
                set_output_format(&mut options, ScanOutputFormat::Ms2)?;
            }
            "--format" => {
                let value = iter.next().context("--format expects tsv or ms2")?;
                let output_format = ScanOutputFormat::parse(&value)?;
                set_output_format(&mut options, output_format)?;
            }
            other => anyhow::bail!("unknown extract option `{other}`"),
        }
    }

    if options.mzml_path.is_none() {
        anyhow::bail!("extract requires --mzml <file>");
    }
    if options.selector.is_none() {
        anyhow::bail!("extract requires --scan <n>");
    }

    Ok(options)
}

fn set_selector(options: &mut ScanOptions, selector: SpectrumSelector) -> anyhow::Result<()> {
    if options.selector.is_some() {
        anyhow::bail!("specify only one --scan selector");
    }
    options.selector = Some(selector);
    Ok(())
}

fn set_output_path(options: &mut ScanOptions, output_path: PathBuf) -> anyhow::Result<()> {
    if options.output_path.is_some() {
        anyhow::bail!("specify only one output path");
    }
    options.output_path = Some(output_path);
    Ok(())
}

fn set_output_format(
    options: &mut ScanOptions,
    output_format: ScanOutputFormat,
) -> anyhow::Result<()> {
    if let Some(existing) = options.output_format {
        if existing != output_format {
            anyhow::bail!(
                "conflicting output format selections: `{}` and `{}`",
                existing.as_str(),
                output_format.as_str()
            );
        }
        return Ok(());
    }
    options.output_format = Some(output_format);
    Ok(())
}

fn parse_scan_number_arg(raw: &str) -> anyhow::Result<u64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        anyhow::bail!("--scan expects a non-empty value");
    }
    extract_scan_number(trimmed)
        .ok_or_else(|| anyhow::anyhow!("invalid --scan `{raw}` (expected 107468 or scan=107468)"))
}

fn resolve_output_format(options: &ScanOptions) -> anyhow::Result<ScanOutputFormat> {
    let inferred = options
        .output_path
        .as_deref()
        .and_then(ScanOutputFormat::from_path_extension);

    if let (Some(explicit), Some(inferred)) = (options.output_format, inferred) {
        if explicit != inferred {
            let path = options
                .output_path
                .as_ref()
                .expect("inferred extension requires output path");
            anyhow::bail!(
                "output path `{}` suggests `{}`, but the requested format is `{}`",
                path.display(),
                inferred.as_str(),
                explicit.as_str()
            );
        }
    }

    Ok(options
        .output_format
        .or(inferred)
        .unwrap_or(ScanOutputFormat::Tsv))
}

fn render_scan_export(
    path: &Path,
    spectrum: &LoadedSpectrum,
    output_format: ScanOutputFormat,
) -> String {
    match output_format {
        ScanOutputFormat::Tsv => render_scan_tsv(path, spectrum),
        ScanOutputFormat::Ms2 => render_scan_ms2(path, spectrum),
    }
}

fn render_scan_tsv(path: &Path, spectrum: &LoadedSpectrum) -> String {
    let mut out = String::new();
    out.push_str(&format!("# source\t{}\n", path.display()));
    if let Some(scan_number) = extract_scan_number(&spectrum.meta.scan_id) {
        out.push_str(&format!("# scan_number\t{scan_number}\n"));
    }
    out.push_str(&format!("# scan_id\t{}\n", spectrum.meta.scan_id));
    out.push_str(&format!("# scan_index\t{}\n", spectrum.meta.idx));
    out.push_str(&format!("# ms_level\t{}\n", spectrum.meta.ms_level));
    out.push_str(&format!(
        "# continuity\t{}\n",
        continuity_label(spectrum.meta.continuity)
    ));
    out.push_str(&format!(
        "# rt_minutes\t{}\n",
        format_optional_f64(spectrum.meta.rt_minutes.map(f64::from), 6)
    ));
    out.push_str(&format!(
        "# precursor_mz\t{}\n",
        format_optional_f64(spectrum.meta.precursor_mz, 6)
    ));
    out.push_str(&format!(
        "# precursor_charge\t{}\n",
        spectrum
            .meta
            .precursor_charge
            .map(|value| value.to_string())
            .unwrap_or_default()
    ));
    out.push_str(&format!("# points\t{}\n", spectrum.stats.points));
    out.push_str(&format!("# mz_min\t{:.6}\n", spectrum.stats.mz_min));
    out.push_str(&format!("# mz_max\t{:.6}\n", spectrum.stats.mz_max));
    out.push_str(&format!(
        "# base_peak_mz\t{:.6}\n",
        spectrum.stats.base_peak_mz
    ));
    out.push_str(&format!(
        "# base_peak_intensity\t{:.6}\n",
        spectrum.stats.base_peak_intensity
    ));
    out.push_str("mz\tintensity\n");
    for (&mz, &intensity) in spectrum.mz.iter().zip(spectrum.intensity.iter()) {
        out.push_str(&format!("{mz:.6}\t{intensity:.6}\n"));
    }
    out
}

fn render_scan_ms2(path: &Path, spectrum: &LoadedSpectrum) -> String {
    let scan_number = spectrum_scan_number(spectrum);
    let mut out = String::new();
    out.push_str("H\tExtractor\tmzio\n");
    out.push_str(&format!("H\tSource\t{}\n", path.display()));
    out.push_str(&format!("H\tNativeID\t{}\n", spectrum.meta.scan_id));
    match spectrum.meta.precursor_mz {
        Some(precursor_mz) => {
            out.push_str(&format!(
                "S\t{scan_number}\t{scan_number}\t{precursor_mz:.6}\n"
            ));
        }
        None => {
            out.push_str(&format!("S\t{scan_number}\t{scan_number}\n"));
        }
    }
    if let Some(rt_minutes) = spectrum.meta.rt_minutes {
        out.push_str(&format!("I\tRTime\t{:.6}\n", rt_minutes));
    }
    out.push_str(&format!("I\tMSLevel\t{}\n", spectrum.meta.ms_level));
    out.push_str(&format!(
        "I\tSignalContinuity\t{}\n",
        continuity_label(spectrum.meta.continuity)
    ));
    if let Some((charge, neutral_mass)) = precursor_charge_and_neutral_mass(spectrum) {
        out.push_str(&format!("Z\t{charge}\t{neutral_mass:.6}\n"));
    }
    for (&mz, &intensity) in spectrum.mz.iter().zip(spectrum.intensity.iter()) {
        out.push_str(&format!("{mz:.6}\t{intensity:.6}\n"));
    }
    out
}

fn precursor_charge_and_neutral_mass(spectrum: &LoadedSpectrum) -> Option<(i32, f64)> {
    let charge = spectrum.meta.precursor_charge?;
    if charge == 0 {
        return None;
    }
    let precursor_mz = spectrum.meta.precursor_mz?;
    let charge_f64 = charge as f64;
    Some((charge, charge_f64 * (precursor_mz - PROTON_MASS)))
}

fn spectrum_scan_number(spectrum: &LoadedSpectrum) -> u64 {
    extract_scan_number(&spectrum.meta.scan_id).unwrap_or(spectrum.meta.idx as u64)
}

fn continuity_label(continuity: SignalContinuity) -> &'static str {
    match continuity {
        SignalContinuity::Centroid => "centroid",
        SignalContinuity::Profile => "profile",
        SignalContinuity::Unknown => "unknown",
    }
}

fn format_optional_f64(value: Option<f64>, precision: usize) -> String {
    match value {
        Some(value) => format!("{value:.precision$}"),
        None => String::new(),
    }
}

impl ScanOutputFormat {
    fn parse(value: &str) -> anyhow::Result<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "tsv" | "tab" | "table" | "txt" => Ok(Self::Tsv),
            "ms2" => Ok(Self::Ms2),
            other => anyhow::bail!("unknown output format `{other}` (expected `tsv` or `ms2`)"),
        }
    }

    fn from_path_extension(path: &Path) -> Option<Self> {
        let extension = path.extension()?.to_str()?.trim();
        if extension.is_empty() {
            return None;
        }
        match extension.to_ascii_lowercase().as_str() {
            "ms2" => Some(Self::Ms2),
            "tsv" | "tab" | "txt" => Some(Self::Tsv),
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Tsv => "tsv",
            Self::Ms2 => "ms2",
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Tsv => "TSV",
            Self::Ms2 => "MS2",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        parse_args, parse_scan_number_arg, precursor_charge_and_neutral_mass, render_scan_export,
        resolve_output_format, ScanOutputFormat,
    };
    use crate::mzml::{LoadedSpectrum, SpectrumMeta, SpectrumStats};
    use mzdata::spectrum::SignalContinuity;
    use std::path::{Path, PathBuf};

    fn demo_spectrum() -> LoadedSpectrum {
        LoadedSpectrum {
            meta: SpectrumMeta {
                idx: 3,
                scan_id: "scan=42".to_string(),
                ms_level: 2,
                rt_minutes: Some(1.25),
                precursor_mz: Some(500.2),
                precursor_charge: Some(2),
                continuity: SignalContinuity::Centroid,
            },
            mz: vec![150.0, 250.0],
            intensity: vec![5.0, 20.0],
            stats: SpectrumStats {
                points: 2,
                mz_min: 150.0,
                mz_max: 250.0,
                base_peak_mz: 250.0,
                base_peak_intensity: 20.0,
            },
        }
    }

    #[test]
    fn parse_args_accepts_scan_and_output() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--scan".into(),
            "scan=4821".into(),
            "--output".into(),
            "scan4821.tsv".into(),
        ])
        .expect("parse valid scan args");
        assert!(options.mzml_path.is_some());
        assert!(options.selector.is_some());
        assert!(options.output_path.is_some());
        assert_eq!(
            resolve_output_format(&options).expect("resolve format"),
            ScanOutputFormat::Tsv
        );
    }

    #[test]
    fn parse_args_accepts_explicit_ms2_output() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--scan".into(),
            "4821".into(),
            "--ms2".into(),
            "scan4821.export".into(),
        ])
        .expect("parse valid ms2 args");
        assert_eq!(options.output_path, Some(PathBuf::from("scan4821.export")));
        assert_eq!(options.output_format, Some(ScanOutputFormat::Ms2));
    }

    #[test]
    fn parse_scan_number_arg_accepts_scan_prefix() {
        let scan = parse_scan_number_arg("scan=107468").expect("parse scan");
        assert_eq!(scan, 107468);
    }

    #[test]
    fn parse_scan_number_arg_accepts_plain_integer() {
        let scan = parse_scan_number_arg("107468").expect("parse scan");
        assert_eq!(scan, 107468);
    }

    #[test]
    fn resolve_output_format_detects_ms2_extension() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--scan".into(),
            "42".into(),
            "--output".into(),
            "scan42.ms2".into(),
        ])
        .expect("parse args");
        assert_eq!(
            resolve_output_format(&options).expect("resolve format"),
            ScanOutputFormat::Ms2
        );
    }

    #[test]
    fn resolve_output_format_rejects_conflicting_extension_and_flag() {
        let options = parse_args(vec![
            "--mzml".into(),
            "demo.mzML".into(),
            "--scan".into(),
            "42".into(),
            "--output".into(),
            "scan42.ms2".into(),
            "--format".into(),
            "tsv".into(),
        ])
        .expect("parse args");
        let error = resolve_output_format(&options).expect_err("expected conflict");
        assert!(error.to_string().contains("suggests `ms2`"));
    }

    #[test]
    fn render_scan_export_includes_metadata_and_peaks_for_tsv() {
        let spectrum = demo_spectrum();
        let text = render_scan_export(Path::new("demo.mzML"), &spectrum, ScanOutputFormat::Tsv);
        assert!(text.contains("# scan_number\t42"));
        assert!(text.contains("# ms_level\t2"));
        assert!(text.contains("mz\tintensity"));
        assert!(text.contains("150.000000\t5.000000"));
    }

    #[test]
    fn render_scan_export_writes_ms2_record() {
        let spectrum = demo_spectrum();
        let text = render_scan_export(Path::new("demo.mzML"), &spectrum, ScanOutputFormat::Ms2);
        assert!(text.contains("H\tExtractor\tmzio"));
        assert!(text.contains("S\t42\t42\t500.200000"));
        assert!(text.contains("I\tRTime\t1.250000"));
        assert!(text.contains("Z\t2\t998.385447"));
        assert!(text.contains("150.000000\t5.000000"));
    }

    #[test]
    fn precursor_charge_and_neutral_mass_matches_ms2_formula() {
        let spectrum = demo_spectrum();
        let (charge, neutral_mass) =
            precursor_charge_and_neutral_mass(&spectrum).expect("neutral mass");
        assert_eq!(charge, 2);
        assert!((neutral_mass - 998.385447066376).abs() < 1e-9);
    }
}
