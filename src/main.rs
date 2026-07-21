mod annotate;
mod dia;
mod diagnostics;
mod ion_table;
mod ms2;
mod mzml;
mod pepxml;
mod plot;
mod scale;
mod scan;
mod spectra;
mod survey;
mod svg_canvas;

use std::env;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::Context;

const PLANNED_COMMANDS: &[&str] = &["summary", "find", "serve"];

pub(crate) fn program_name() -> String {
    env::args()
        .next()
        .as_deref()
        .and_then(|arg0| Path::new(arg0).file_name().and_then(|value| value.to_str()))
        .filter(|name| !name.is_empty())
        .unwrap_or("mzio")
        .to_string()
}

fn print_help() {
    let program = program_name();
    println!("{program}");
    println!();
    println!("Standalone mzML and MS/MS tooling extracted from universal-tui.");
    println!();
    println!("USAGE:");
    println!("  {program} <command> [options]");
    println!("  {program} --help");
    println!();
    println!("AVAILABLE COMMANDS:");
    println!("  dia-slice    Export DIA slice summaries from mzML or Bruker .d");
    println!("  diagnostics  Scan centroided MS2 spectra for diagnostic ions and deltas");
    println!("  diag         Alias for diagnostics");
    println!("  extract      Export one mzML spectrum by scan number");
    println!("  plot         Export an mzML or ms2 spectrum as SVG");
    println!("  plot-survey  Export run-level survey and DDA QC SVGs");
    println!("  scan         Alias for extract");
    println!("  spectra      Browse mzML spectra in a local TUI with caching");
    for command in PLANNED_COMMANDS {
        println!("  {command}      Planned");
    }
    println!();
    println!("RUN:");
    println!("  {program} dia-slice --help");
    println!("  {program} extract --help");
    println!("  {program} plot --help");
    println!("  {program} plot-survey --help");
    println!("  {program} scan --help");
    println!("  {program} spectra --help");
}

fn print_spectra_help() {
    let program = program_name();
    println!("{program} spectra");
    println!();
    println!("USAGE:");
    println!("  {program} spectra --mzml <file> [options]");
    println!();
    println!("OPTIONS:");
    println!("  --mzml <file>             Input mzML file to browse");
    println!("  --no-mzml-cache          Disable on-disk index caching");
    println!("  --reindex                Rebuild the cache instead of reusing it");
    println!("  --mzml-cache-path <p>    Override the cache file or cache directory path");
    println!("  --help                   Show this help");
    println!();
    println!("KEYS:");
    println!("  j/k, arrows, PageUp/PageDown, Home/End   Move through spectra");
    println!("  Enter                                    Load selected spectrum");
    println!(
        "  o                                        Cycle right pane (spectrum/chrom/map/overview)"
    );
    println!("  n                                        Toggle normalization");
    println!("  z or Tab                                 Cycle zoom preset");
    println!("  m                                        Cycle plot mode");
    println!("  /, s, S                                  Search scan ids forward/backward");
    println!("  p / P                                    Export SVG / PNG");
    println!("  q or Esc                                 Quit");
}

fn print_stub(command: &str) {
    println!("`{command}` is not implemented yet.");
    println!("This package is set up as the standalone home for upcoming mzML work.");
}

fn parse_spectra_args(args: Vec<String>) -> anyhow::Result<(PathBuf, spectra::SpectraOptions)> {
    let mut mzml_path = None::<PathBuf>;
    let mut options = spectra::SpectraOptions::default();
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--mzml" => {
                mzml_path = Some(PathBuf::from(iter.next().context("--mzml expects a path")?));
            }
            "--no-mzml-cache" => {
                options.index_cache.enabled = false;
            }
            "--reindex" => {
                options.index_cache.refresh = true;
            }
            "--mzml-cache-path" => {
                options.index_cache.path = Some(PathBuf::from(
                    iter.next().context("--mzml-cache-path expects a value")?,
                ));
            }
            other => anyhow::bail!("unknown spectra option `{other}`"),
        }
    }

    let mzml_path = mzml_path.ok_or_else(|| anyhow::anyhow!("spectra requires --mzml <file>"))?;
    Ok((mzml_path, options))
}

fn run() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        None | Some("-h") | Some("--help") | Some("help") => {
            print_help();
            Ok(())
        }
        Some("diagnostics") | Some("diag") => diagnostics::run(args.collect()),
        Some("dia-slice") | Some("dia") => dia::run(args.collect()),
        Some("extract") => scan::run("extract", args.collect()),
        Some("plot") => plot::run(args.collect()),
        Some("plot-survey") | Some("survey") => survey::run(args.collect()),
        Some("scan") => scan::run("scan", args.collect()),
        Some("spectra") | Some("browse") => {
            let sub_args = args.collect::<Vec<_>>();
            if sub_args
                .iter()
                .any(|arg| matches!(arg.as_str(), "-h" | "--help" | "help"))
            {
                print_spectra_help();
                return Ok(());
            }
            let (path, options) = parse_spectra_args(sub_args)?;
            let runtime = tokio::runtime::Builder::new_multi_thread()
                .enable_all()
                .build()
                .context("failed to build tokio runtime for spectra browser")?;
            runtime.block_on(spectra::run_spectra_demo(path, options))
        }
        Some(command) if PLANNED_COMMANDS.contains(&command) => {
            print_stub(command);
            Ok(())
        }
        Some(other) => anyhow::bail!("unknown command `{other}`; run `{} --help`", program_name()),
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err:#}");
            ExitCode::from(1)
        }
    }
}
