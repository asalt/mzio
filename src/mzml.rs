use std::fs::File;
use std::path::Path;

use anyhow::Context;
use mzdata::{
    io::{DetailLevel, MzMLReader, SpectrumSource},
    prelude::*,
    spectrum::{
        bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray},
        IonProperties, PrecursorSelection, RefPeakDataLevel, SignalContinuity,
    },
};
use mzpeaks::{CentroidLike, DeconvolutedCentroidLike};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum SpectrumSelector {
    Index(u32),
    ScanNumber(u64),
    NativeId(String),
}

#[derive(Clone, Debug)]
pub(crate) struct SpectrumMeta {
    pub(crate) idx: u32,
    pub(crate) scan_id: String,
    pub(crate) ms_level: u8,
    pub(crate) rt_minutes: Option<f32>,
    pub(crate) precursor_mz: Option<f64>,
    pub(crate) precursor_charge: Option<i32>,
    pub(crate) continuity: SignalContinuity,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct SpectrumStats {
    pub(crate) points: u32,
    pub(crate) mz_min: f64,
    pub(crate) mz_max: f64,
    pub(crate) base_peak_mz: f64,
    pub(crate) base_peak_intensity: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct LoadedSpectrum {
    pub(crate) meta: SpectrumMeta,
    pub(crate) mz: Vec<f64>,
    pub(crate) intensity: Vec<f32>,
    pub(crate) stats: SpectrumStats,
}

impl LoadedSpectrum {
    pub(crate) fn new(meta: SpectrumMeta, mz: Vec<f64>, intensity: Vec<f32>) -> Self {
        let stats = compute_stats(&mz, &intensity);
        Self {
            meta,
            mz,
            intensity,
            stats,
        }
    }
}

pub(crate) fn open_reader(path: &Path) -> anyhow::Result<MzMLReader<File>> {
    let mut reader = MzMLReader::open_path(path)
        .with_context(|| format!("failed to open mzML at {}", path.display()))?;
    reader.detail_level = DetailLevel::Lazy;
    Ok(reader)
}

pub(crate) fn load_selected_spectrum(
    reader: &mut MzMLReader<File>,
    selector: &SpectrumSelector,
) -> anyhow::Result<LoadedSpectrum> {
    let idx = match selector {
        SpectrumSelector::Index(idx) => *idx,
        SpectrumSelector::ScanNumber(scan_number) => {
            find_index_by_scan_number(reader, *scan_number)?
        }
        SpectrumSelector::NativeId(id) => find_index_by_id(reader, id)?,
    };
    load_spectrum_by_index(reader, idx)
}

pub(crate) fn load_spectrum_by_index(
    reader: &mut MzMLReader<File>,
    idx: u32,
) -> anyhow::Result<LoadedSpectrum> {
    let spec = reader
        .get_spectrum_by_index(idx as usize)
        .ok_or_else(|| anyhow::anyhow!("index {idx} out of range"))?;

    let (precursor_mz, precursor_charge) = selected_precursor_details(&spec);
    let meta = SpectrumMeta {
        idx,
        scan_id: spec.id().to_string(),
        ms_level: spec.ms_level(),
        rt_minutes: {
            let value = spec.start_time() as f32;
            if value > 0.0 {
                Some(value)
            } else {
                None
            }
        },
        precursor_mz,
        precursor_charge,
        continuity: spec.signal_continuity(),
    };

    let (mz, intensity) = load_peak_arrays(spec.raw_arrays(), spec.peaks(), idx)?;
    if mz.len() != intensity.len() {
        anyhow::bail!(
            "spectrum {idx}: mz/intensity length mismatch (mz.len={} intensity.len={})",
            mz.len(),
            intensity.len()
        );
    }

    Ok(LoadedSpectrum::new(meta, mz, intensity))
}

fn selected_precursor_details<C, D, S>(spec: &S) -> (Option<f64>, Option<i32>)
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    let Some(ion) = spec.precursor().and_then(|precursor| precursor.ion()) else {
        return (None, None);
    };
    (Some(ion.mz()), ion.charge())
}

pub(crate) fn find_index_by_id(
    reader: &mut MzMLReader<File>,
    target_id: &str,
) -> anyhow::Result<u32> {
    let offsets = reader.get_index();
    for (idx, (native_id, _)) in offsets.iter().enumerate() {
        if native_id_matches_query(native_id, target_id) {
            return Ok(idx as u32);
        }
    }

    let total = reader.len();
    for idx in 0..total {
        let spec = reader
            .get_spectrum_by_index(idx)
            .ok_or_else(|| anyhow::anyhow!("index {idx} missing while scanning for native id"))?;
        if native_id_matches_query(spec.id(), target_id) {
            return Ok(idx as u32);
        }
    }

    anyhow::bail!("native id `{target_id}` not found")
}

pub(crate) fn find_index_by_scan_number(
    reader: &mut MzMLReader<File>,
    target_scan: u64,
) -> anyhow::Result<u32> {
    let offsets = reader.get_index();
    for (idx, (native_id, _)) in offsets.iter().enumerate() {
        if extract_scan_number(native_id) == Some(target_scan) {
            return Ok(idx as u32);
        }
    }

    let total = reader.len();
    for idx in 0..total {
        let spec = reader
            .get_spectrum_by_index(idx)
            .ok_or_else(|| anyhow::anyhow!("index {idx} missing while scanning for scan number"))?;
        if extract_scan_number(spec.id()) == Some(target_scan) {
            return Ok(idx as u32);
        }
    }

    anyhow::bail!("scan `{target_scan}` not found")
}

pub(crate) fn native_id_matches_query(candidate: &str, query: &str) -> bool {
    if candidate == query {
        return true;
    }

    if candidate
        .split_ascii_whitespace()
        .any(|token| token == query || scan_token_matches(token, query))
    {
        return true;
    }

    scan_token_matches(candidate, query)
}

pub(crate) fn extract_scan_number(input: &str) -> Option<u64> {
    let trimmed = input.trim();
    if let Ok(value) = trimmed.parse::<u64>() {
        return Some(value);
    }

    let token = trimmed
        .split_ascii_whitespace()
        .find(|piece| piece.starts_with("scan="))
        .unwrap_or(trimmed);
    let value = token.strip_prefix("scan=")?;
    value.parse::<u64>().ok()
}

fn scan_token_matches(candidate: &str, query: &str) -> bool {
    match (extract_scan_number(candidate), extract_scan_number(query)) {
        (Some(left), Some(right)) => left == right,
        _ => false,
    }
}

fn load_peak_arrays<C, D>(
    raw_arrays: Option<&BinaryArrayMap>,
    peaks: RefPeakDataLevel<'_, C, D>,
    idx: u32,
) -> anyhow::Result<(Vec<f64>, Vec<f32>)>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
{
    if let Some(arrays) = raw_arrays {
        if let Ok((mz, intensity)) = decode_mz_intensity_from_arrays(arrays, idx) {
            return Ok((mz, intensity));
        }
    }

    let mut mz = Vec::<f64>::new();
    let mut intensity = Vec::<f32>::new();
    match peaks {
        RefPeakDataLevel::Missing => {}
        RefPeakDataLevel::RawData(arrays) => {
            if let Ok((decoded_mz, decoded_intensity)) =
                decode_mz_intensity_from_arrays(arrays, idx)
            {
                mz = decoded_mz;
                intensity = decoded_intensity;
            }
        }
        RefPeakDataLevel::Centroid(centroids) => {
            mz.reserve(centroids.len());
            intensity.reserve(centroids.len());
            for peak in centroids.iter() {
                let centroid = peak.as_centroid();
                mz.push(centroid.mz);
                intensity.push(centroid.intensity);
            }
        }
        RefPeakDataLevel::Deconvoluted(peaks) => {
            mz.reserve(peaks.len());
            intensity.reserve(peaks.len());
            for peak in peaks.iter() {
                let centroid = peak.as_centroid();
                let charge = centroid.charge as f64;
                let mz_value = if charge.abs() < f64::EPSILON {
                    centroid.neutral_mass
                } else {
                    (centroid.neutral_mass + 1.007_276 * charge) / charge
                };
                mz.push(mz_value);
                intensity.push(centroid.intensity);
            }
        }
    }
    Ok((mz, intensity))
}

fn decode_mz_intensity_from_arrays(
    arrays: &BinaryArrayMap,
    idx: u32,
) -> anyhow::Result<(Vec<f64>, Vec<f32>)> {
    let mz_array = arrays
        .get(&ArrayType::MZArray)
        .ok_or_else(|| anyhow::anyhow!("spectrum {idx}: missing m/z array"))?;
    let intensity_array = arrays
        .get(&ArrayType::IntensityArray)
        .ok_or_else(|| anyhow::anyhow!("spectrum {idx}: missing intensity array"))?;

    let mz_bytes = decode_data_array_bytes(mz_array, idx, "m/z")?;
    let intensity_bytes = decode_data_array_bytes(intensity_array, idx, "intensity")?;

    let mz = decode_bytes_as_f64(&mz_bytes, mz_array.dtype, idx, "m/z")?;
    let intensity = decode_bytes_as_f32(&intensity_bytes, intensity_array.dtype, idx, "intensity")?;

    Ok((mz, intensity))
}

fn decode_data_array_bytes(
    array: &DataArray,
    idx: u32,
    label: &'static str,
) -> anyhow::Result<Vec<u8>> {
    array
        .decode()
        .map(|cow| cow.into_owned())
        .with_context(|| format!("failed to decode {label} bytes for spectrum {idx}"))
}

fn decode_bytes_as_f64(
    bytes: &[u8],
    dtype: BinaryDataArrayType,
    idx: u32,
    label: &'static str,
) -> anyhow::Result<Vec<f64>> {
    match dtype {
        BinaryDataArrayType::Unknown => {
            anyhow::bail!("spectrum {idx}: {label} array has unknown dtype")
        }
        BinaryDataArrayType::Float64 => {
            if bytes.len() % 8 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 8",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 8);
            for chunk in bytes.chunks_exact(8) {
                let raw: [u8; 8] = chunk
                    .try_into()
                    .expect("chunks_exact(8) yields 8-byte slices");
                out.push(f64::from_le_bytes(raw));
            }
            Ok(out)
        }
        BinaryDataArrayType::Float32 => {
            if bytes.len() % 4 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 4",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 4);
            for chunk in bytes.chunks_exact(4) {
                let raw: [u8; 4] = chunk
                    .try_into()
                    .expect("chunks_exact(4) yields 4-byte slices");
                out.push(f32::from_le_bytes(raw) as f64);
            }
            Ok(out)
        }
        BinaryDataArrayType::Int32 => {
            if bytes.len() % 4 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 4",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 4);
            for chunk in bytes.chunks_exact(4) {
                let raw: [u8; 4] = chunk
                    .try_into()
                    .expect("chunks_exact(4) yields 4-byte slices");
                out.push(i32::from_le_bytes(raw) as f64);
            }
            Ok(out)
        }
        BinaryDataArrayType::Int64 => {
            if bytes.len() % 8 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 8",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 8);
            for chunk in bytes.chunks_exact(8) {
                let raw: [u8; 8] = chunk
                    .try_into()
                    .expect("chunks_exact(8) yields 8-byte slices");
                out.push(i64::from_le_bytes(raw) as f64);
            }
            Ok(out)
        }
        BinaryDataArrayType::ASCII => {
            anyhow::bail!("spectrum {idx}: {label} array is ASCII-encoded, not supported")
        }
    }
}

fn decode_bytes_as_f32(
    bytes: &[u8],
    dtype: BinaryDataArrayType,
    idx: u32,
    label: &'static str,
) -> anyhow::Result<Vec<f32>> {
    match dtype {
        BinaryDataArrayType::Unknown => {
            anyhow::bail!("spectrum {idx}: {label} array has unknown dtype")
        }
        BinaryDataArrayType::Float32 => {
            if bytes.len() % 4 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 4",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 4);
            for chunk in bytes.chunks_exact(4) {
                let raw: [u8; 4] = chunk
                    .try_into()
                    .expect("chunks_exact(4) yields 4-byte slices");
                out.push(f32::from_le_bytes(raw));
            }
            Ok(out)
        }
        BinaryDataArrayType::Float64 => {
            if bytes.len() % 8 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 8",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 8);
            for chunk in bytes.chunks_exact(8) {
                let raw: [u8; 8] = chunk
                    .try_into()
                    .expect("chunks_exact(8) yields 8-byte slices");
                out.push(f64::from_le_bytes(raw) as f32);
            }
            Ok(out)
        }
        BinaryDataArrayType::Int32 => {
            if bytes.len() % 4 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 4",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 4);
            for chunk in bytes.chunks_exact(4) {
                let raw: [u8; 4] = chunk
                    .try_into()
                    .expect("chunks_exact(4) yields 4-byte slices");
                out.push(i32::from_le_bytes(raw) as f32);
            }
            Ok(out)
        }
        BinaryDataArrayType::Int64 => {
            if bytes.len() % 8 != 0 {
                anyhow::bail!(
                    "spectrum {idx}: {label} byte length {} not divisible by 8",
                    bytes.len()
                );
            }
            let mut out = Vec::with_capacity(bytes.len() / 8);
            for chunk in bytes.chunks_exact(8) {
                let raw: [u8; 8] = chunk
                    .try_into()
                    .expect("chunks_exact(8) yields 8-byte slices");
                out.push(i64::from_le_bytes(raw) as f32);
            }
            Ok(out)
        }
        BinaryDataArrayType::ASCII => {
            anyhow::bail!("spectrum {idx}: {label} array is ASCII-encoded, not supported")
        }
    }
}

fn compute_stats(mz: &[f64], intensity: &[f32]) -> SpectrumStats {
    let mut points: u32 = 0;
    let mut mz_min = f64::INFINITY;
    let mut mz_max = -f64::INFINITY;
    let mut base_peak_mz = 0.0_f64;
    let mut base_peak_intensity = -f32::INFINITY;

    for (&m, &i) in mz.iter().zip(intensity.iter()) {
        points = points.saturating_add(1);
        if m.is_finite() {
            mz_min = mz_min.min(m);
            mz_max = mz_max.max(m);
        }
        if i.is_finite() && i > base_peak_intensity {
            base_peak_intensity = i;
            base_peak_mz = m;
        }
    }

    if points == 0 {
        SpectrumStats {
            points,
            mz_min: 0.0,
            mz_max: 0.0,
            base_peak_mz: 0.0,
            base_peak_intensity: 0.0,
        }
    } else {
        SpectrumStats {
            points,
            mz_min: if mz_min.is_finite() { mz_min } else { 0.0 },
            mz_max: if mz_max.is_finite() { mz_max } else { 0.0 },
            base_peak_mz,
            base_peak_intensity: if base_peak_intensity.is_finite() {
                base_peak_intensity
            } else {
                0.0
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{extract_scan_number, native_id_matches_query};

    #[test]
    fn extract_scan_number_accepts_plain_integer() {
        assert_eq!(extract_scan_number("107468"), Some(107468));
    }

    #[test]
    fn native_id_matches_plain_numeric_query() {
        assert!(native_id_matches_query(
            "controllerType=0 controllerNumber=1 scan=107468",
            "107468",
        ));
    }
}
