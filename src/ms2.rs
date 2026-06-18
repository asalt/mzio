use std::fs;
use std::path::Path;

use anyhow::Context;
use mzdata::spectrum::SignalContinuity;

use crate::mzml::{extract_scan_number, LoadedSpectrum, SpectrumMeta, SpectrumSelector};

const PROTON_MASS: f64 = 1.007_276_466_812;

#[derive(Clone, Debug)]
struct Ms2SpectrumRecord {
    idx: u32,
    scan_start: u64,
    scan_end: u64,
    precursor_mz: Option<f64>,
    rt_minutes: Option<f32>,
    precursor_charge: Option<i32>,
    mz: Vec<f64>,
    intensity: Vec<f32>,
}

impl Ms2SpectrumRecord {
    fn new(idx: u32, scan_start: u64, scan_end: u64, precursor_mz: Option<f64>) -> Self {
        Self {
            idx,
            scan_start,
            scan_end,
            precursor_mz,
            rt_minutes: None,
            precursor_charge: None,
            mz: Vec::new(),
            intensity: Vec::new(),
        }
    }

    fn scan_id(&self) -> String {
        format!("scan={}", self.scan_start)
    }

    fn matches_scan_number(&self, scan_number: u64) -> bool {
        self.scan_start == scan_number || self.scan_end == scan_number
    }

    fn matches_native_id(&self, query: &str) -> bool {
        let query = query.trim();
        if query.is_empty() {
            return false;
        }
        if let Some(scan_number) = extract_scan_number(query) {
            return self.matches_scan_number(scan_number);
        }

        let scan_id = self.scan_id();
        scan_id == query || scan_id.contains(query)
    }

    fn into_loaded_spectrum(self) -> LoadedSpectrum {
        LoadedSpectrum::new(
            SpectrumMeta {
                idx: self.idx,
                scan_id: self.scan_id(),
                ms_level: 2,
                rt_minutes: self.rt_minutes,
                precursor_mz: self.precursor_mz,
                precursor_charge: self.precursor_charge,
                continuity: SignalContinuity::Centroid,
            },
            self.mz,
            self.intensity,
        )
    }
}

pub(crate) fn load_selected_spectrum(
    path: &Path,
    selector: Option<&SpectrumSelector>,
) -> anyhow::Result<LoadedSpectrum> {
    let spectra = parse_ms2_file(path)?;
    if spectra.is_empty() {
        anyhow::bail!("no spectra found in {}", path.display());
    }

    let selected = match selector {
        Some(SpectrumSelector::Index(idx)) => spectra
            .into_iter()
            .find(|record| record.idx == *idx)
            .ok_or_else(|| anyhow::anyhow!("index {idx} out of range for {}", path.display()))?,
        Some(SpectrumSelector::ScanNumber(scan_number)) => spectra
            .into_iter()
            .find(|record| record.matches_scan_number(*scan_number))
            .ok_or_else(|| anyhow::anyhow!("scan {scan_number} not found in {}", path.display()))?,
        Some(SpectrumSelector::NativeId(id)) => spectra
            .into_iter()
            .find(|record| record.matches_native_id(id))
            .ok_or_else(|| anyhow::anyhow!("id `{id}` not found in {}", path.display()))?,
        None => spectra
            .into_iter()
            .next()
            .expect("checked non-empty spectra"),
    };

    Ok(selected.into_loaded_spectrum())
}

fn parse_ms2_file(path: &Path) -> anyhow::Result<Vec<Ms2SpectrumRecord>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read ms2 file at {}", path.display()))?;

    let mut spectra = Vec::<Ms2SpectrumRecord>::new();
    let mut current = None::<Ms2SpectrumRecord>;

    for (line_idx, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line.starts_with("H\t") {
            continue;
        }

        if line.starts_with("S\t") {
            if let Some(record) = current.take() {
                spectra.push(record);
            }
            current = Some(parse_s_line(line, spectra.len() as u32, line_idx + 1)?);
            continue;
        }

        if line.starts_with("I\t") {
            if let Some(record) = current.as_mut() {
                parse_i_line(record, line);
            }
            continue;
        }

        if line.starts_with("Z\t") {
            if let Some(record) = current.as_mut() {
                parse_z_line(record, line);
            }
            continue;
        }

        if line.starts_with('D') {
            continue;
        }

        if let Some(record) = current.as_mut() {
            parse_peak_line(record, line, line_idx + 1)?;
        } else {
            anyhow::bail!(
                "encountered peak data before first spectrum header at {} line {}",
                path.display(),
                line_idx + 1
            );
        }
    }

    if let Some(record) = current.take() {
        spectra.push(record);
    }

    Ok(spectra)
}

fn parse_s_line(line: &str, idx: u32, line_number: usize) -> anyhow::Result<Ms2SpectrumRecord> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 3 {
        anyhow::bail!("line {line_number}: invalid S record `{line}`");
    }

    let scan_start = fields[1]
        .parse::<u64>()
        .with_context(|| format!("line {line_number}: invalid S scan start `{}`", fields[1]))?;
    let scan_end = fields[2]
        .parse::<u64>()
        .with_context(|| format!("line {line_number}: invalid S scan end `{}`", fields[2]))?;
    let precursor_mz = fields
        .get(3)
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite());

    Ok(Ms2SpectrumRecord::new(
        idx,
        scan_start,
        scan_end,
        precursor_mz,
    ))
}

fn parse_i_line(record: &mut Ms2SpectrumRecord, line: &str) {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 3 {
        return;
    }

    match fields[1] {
        "RTime" | "RetTime" | "RetentionTime" => {
            record.rt_minutes = fields[2]
                .parse::<f32>()
                .ok()
                .filter(|value| value.is_finite());
        }
        _ => {}
    }
}

fn parse_z_line(record: &mut Ms2SpectrumRecord, line: &str) {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 2 {
        return;
    }

    let charge = match fields[1].parse::<i32>() {
        Ok(value) if value != 0 => value,
        _ => return,
    };
    record.precursor_charge = Some(charge);

    if record.precursor_mz.is_some() || fields.len() < 3 {
        return;
    }

    if let Ok(neutral_mass) = fields[2].parse::<f64>() {
        if neutral_mass.is_finite() {
            let charge_f64 = charge as f64;
            record.precursor_mz = Some((neutral_mass + charge_f64 * PROTON_MASS) / charge_f64);
        }
    }
}

fn parse_peak_line(
    record: &mut Ms2SpectrumRecord,
    line: &str,
    line_number: usize,
) -> anyhow::Result<()> {
    let mut fields = line.split_whitespace();
    let mz = fields
        .next()
        .ok_or_else(|| anyhow::anyhow!("line {line_number}: missing peak m/z"))?
        .parse::<f64>()
        .with_context(|| format!("line {line_number}: invalid peak m/z in `{line}`"))?;
    let intensity = fields
        .next()
        .ok_or_else(|| anyhow::anyhow!("line {line_number}: missing peak intensity"))?
        .parse::<f32>()
        .with_context(|| format!("line {line_number}: invalid peak intensity in `{line}`"))?;

    record.mz.push(mz);
    record.intensity.push(intensity);
    Ok(())
}
