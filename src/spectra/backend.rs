use std::borrow::Cow;
use std::fs;
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::UNIX_EPOCH;

use anyhow::Context;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use mzdata::{
    io::{DetailLevel, MzMLReader, SpectrumSource},
    params::{Param, ParamDescribed},
    prelude::*,
    spectrum::{
        bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray},
        IonProperties, PrecursorSelection, RefPeakDataLevel, SignalContinuity,
    },
};

use super::{spectra_debug_log, spectra_debug_target, IndexCacheOptions};

#[derive(Debug, Clone)]
pub(super) struct SpectrumMeta {
    pub idx: u32,
    pub offset: u64,
    pub scan_id: String,
    pub ms_level: u8,
    pub rt_minutes: Option<f32>,
    pub tic: Option<f32>,
    pub precursor_mz: Option<f64>,
    pub charge: Option<u8>,
    pub continuity: SignalContinuity,
    pub points: Option<u32>,
    pub mz_min: Option<f64>,
    pub mz_max: Option<f64>,
    pub base_peak_mz: Option<f64>,
    pub base_peak_intensity: Option<f32>,
}

#[derive(Debug, Clone, Default)]
pub(super) struct MzmlFileInfo {
    pub run_id: Option<String>,
    pub start_time: Option<String>,
    pub default_instrument_id: Option<u32>,
    pub spectrum_count_hint: Option<u64>,
    pub contents: Vec<String>,
    pub source_files: Vec<String>,
    pub instrument_summaries: Vec<String>,
    pub software_summaries: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SpectrumStats {
    pub points: u32,
    pub mz_min: f64,
    pub mz_max: f64,
    pub base_peak_mz: f64,
    pub base_peak_intensity: f32,
}

#[derive(Debug, Clone)]
pub(super) struct SpectrumOwned {
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f32]>,
    pub stats: SpectrumStats,
}

impl SpectrumOwned {
    pub(super) fn approx_bytes(&self) -> usize {
        self.mz.len() * std::mem::size_of::<f64>()
            + self.intensity.len() * std::mem::size_of::<f32>()
    }
}

#[derive(Debug, Clone)]
pub(super) struct SpectrumData<'a> {
    pub mz: Cow<'a, [f64]>,
    pub intensity: Cow<'a, [f32]>,
}

#[derive(Debug)]
pub(super) struct IndexResult {
    pub metas: Vec<SpectrumMeta>,
    pub file_info: MzmlFileInfo,
    pub cached: bool,
    pub cache_path: Option<PathBuf>,
    pub cache_warning: Option<String>,
}

fn short_param_label(param: &Param) -> String {
    const MAX_VAL: usize = 48;

    let name = param.name.as_str();
    let value = param.value.to_string();
    let value = value.trim();
    if value.is_empty() {
        return name.to_string();
    }

    let mut clipped = value.to_string();
    if clipped.len() > MAX_VAL {
        clipped.truncate(MAX_VAL);
        clipped.push('…');
    }
    format!("{name}={clipped}")
}

fn summarize_params(params: &[Param], max: usize) -> String {
    params
        .iter()
        .filter(|p| !p.name.trim().is_empty())
        .take(max.max(1))
        .map(short_param_label)
        .collect::<Vec<_>>()
        .join(", ")
}

fn selected_precursor_details<C, D, S>(spec: &S) -> (Option<f64>, Option<u8>)
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    let Some(ion) = spec.precursor().and_then(|precursor| precursor.ion()) else {
        return (None, None);
    };
    let charge = ion.charge().and_then(|value| u8::try_from(value).ok());
    (Some(ion.mz()), charge)
}

fn mzml_file_info(reader: &MzMLReader<std::fs::File>) -> MzmlFileInfo {
    let run = reader.run_description();
    let (run_id, start_time, default_instrument_id) = match run {
        Some(run) => (
            run.id.clone(),
            run.start_time.as_ref().map(ToString::to_string),
            run.default_instrument_id,
        ),
        None => (None, None, None),
    };

    let file_description = reader.file_description();
    let contents = file_description
        .contents
        .iter()
        .filter(|p| !p.name.trim().is_empty())
        .take(12)
        .map(|p| p.name.clone())
        .collect::<Vec<_>>();

    let source_files = file_description
        .source_files
        .iter()
        .take(8)
        .map(|sf| {
            let name = sf.name.trim();
            let location = sf.location.trim();
            if !name.is_empty() && !location.is_empty() {
                format!("{name} @ {location}")
            } else if !name.is_empty() {
                name.to_string()
            } else if !location.is_empty() {
                location.to_string()
            } else {
                sf.id.clone()
            }
        })
        .filter(|s| !s.trim().is_empty())
        .collect::<Vec<_>>();

    let mut instrument_summaries = Vec::new();
    let mut instruments: Vec<_> = reader.instrument_configurations().iter().collect();
    instruments.sort_by_key(|(id, _)| **id);
    for (id, inst) in instruments.into_iter().take(6) {
        let mut parts = Vec::new();
        for component in inst.components.iter() {
            let label = component.component_type.to_string();
            let detail = summarize_params(&component.params, 2);
            if detail.is_empty() {
                parts.push(label);
            } else {
                parts.push(format!("{label}: {detail}"));
            }
        }
        if parts.is_empty() {
            let detail = summarize_params(&inst.params, 4);
            if !detail.is_empty() {
                parts.push(detail);
            }
        }
        let body = if parts.is_empty() {
            "(no instrument details)".to_string()
        } else {
            parts.join(" | ")
        };
        instrument_summaries.push(format!("Instrument {id}: {body}"));
    }

    let software_summaries = reader
        .softwares()
        .iter()
        .take(8)
        .map(|sw| {
            let mut base = sw.id.clone();
            if !sw.version.trim().is_empty() {
                base.push_str(&format!(" v{}", sw.version.trim()));
            }
            let detail = summarize_params(&sw.params, 2);
            if detail.is_empty() {
                base
            } else {
                format!("{base} ({detail})")
            }
        })
        .collect::<Vec<_>>();

    MzmlFileInfo {
        run_id,
        start_time,
        default_instrument_id,
        spectrum_count_hint: reader.spectrum_count_hint(),
        contents,
        source_files,
        instrument_summaries,
        software_summaries,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct FileFingerprint {
    len: u64,
    modified_secs: u64,
    modified_nanos: u32,
    signature: u64,
}

impl FileFingerprint {
    fn from_path(path: &Path) -> anyhow::Result<Self> {
        let meta = fs::metadata(path)
            .with_context(|| format!("failed to stat mzML at {}", path.display()))?;
        let len = meta.len();
        let modified = meta.modified().unwrap_or(UNIX_EPOCH);
        let (modified_secs, modified_nanos) = match modified.duration_since(UNIX_EPOCH) {
            Ok(dur) => (dur.as_secs(), dur.subsec_nanos()),
            Err(_) => (0, 0),
        };
        let signature = compute_content_signature(path, len).with_context(|| {
            format!("failed to compute content signature for {}", path.display())
        })?;
        Ok(Self {
            len,
            modified_secs,
            modified_nanos,
            signature,
        })
    }

    fn matches_len_mtime(&self, other: &Self) -> bool {
        self.len == other.len
            && self.modified_secs == other.modified_secs
            && self.modified_nanos == other.modified_nanos
    }
}

const INDEX_CACHE_MAGIC: &[u8; 8] = b"UTUIIDX\0";
const INDEX_CACHE_VERSION: u32 = 3;
const INDEX_CACHE_SUFFIX: &str = ".mzio.idx.gz";
const MAX_CACHE_STRING_LEN: usize = 64 * 1024 * 1024;
const CONTENT_SIGNATURE_WINDOW: usize = 64 * 1024;

fn fnv1a64_update(mut hash: u64, bytes: &[u8]) -> u64 {
    const FNV_PRIME: u64 = 0x100000001b3;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

fn compute_content_signature(path: &Path, len: u64) -> anyhow::Result<u64> {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;

    let mut file = fs::File::open(path)
        .with_context(|| format!("failed to open mzML for signature: {}", path.display()))?;

    let window = CONTENT_SIGNATURE_WINDOW as u64;
    let mut buf = vec![0u8; CONTENT_SIGNATURE_WINDOW];
    let mut hash = FNV_OFFSET;

    let start_len = len.min(window) as usize;
    file.seek(SeekFrom::Start(0))?;
    file.read_exact(&mut buf[..start_len])?;
    hash = fnv1a64_update(hash, &buf[..start_len]);

    if len > window {
        file.seek(SeekFrom::Start(len - window))?;
        file.read_exact(&mut buf)?;
        hash = fnv1a64_update(hash, &buf);
    }

    if len > 3 * window {
        let mid_start = (len / 2).saturating_sub(window / 2);
        file.seek(SeekFrom::Start(mid_start))?;
        file.read_exact(&mut buf)?;
        hash = fnv1a64_update(hash, &buf);
    }

    Ok(hash)
}

fn default_index_cache_path(mzml_path: &Path) -> PathBuf {
    let mut out = mzml_path.to_path_buf();
    let file_name = mzml_path
        .file_name()
        .unwrap_or_else(|| std::ffi::OsStr::new("mzml"));
    let mut name = file_name.to_os_string();
    name.push(INDEX_CACHE_SUFFIX);
    out.set_file_name(name);
    out
}

fn resolve_index_cache_path(mzml_path: &Path, cache: &IndexCacheOptions) -> Option<PathBuf> {
    if !cache.enabled {
        return None;
    }
    let Some(custom) = cache.path.as_ref() else {
        return Some(default_index_cache_path(mzml_path));
    };

    if custom.is_dir() {
        let file_name = mzml_path
            .file_name()
            .unwrap_or_else(|| std::ffi::OsStr::new("mzml"));
        let mut name = file_name.to_os_string();
        name.push(INDEX_CACHE_SUFFIX);
        return Some(custom.join(name));
    }

    Some(custom.clone())
}

fn is_permission_denied(err: &anyhow::Error) -> bool {
    err.chain().any(|cause| {
        cause
            .downcast_ref::<std::io::Error>()
            .is_some_and(|io| io.kind() == std::io::ErrorKind::PermissionDenied)
    })
}

fn cache_write_warning(err: &anyhow::Error) -> String {
    if is_permission_denied(err) {
        "cache write failed (permission denied; use --mzml-cache-path)".to_string()
    } else {
        "cache write failed (see logs)".to_string()
    }
}

fn write_u8(w: &mut impl Write, v: u8) -> anyhow::Result<()> {
    w.write_all(&[v])?;
    Ok(())
}

fn write_u32(w: &mut impl Write, v: u32) -> anyhow::Result<()> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_u64(w: &mut impl Write, v: u64) -> anyhow::Result<()> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_f32(w: &mut impl Write, v: f32) -> anyhow::Result<()> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_f64(w: &mut impl Write, v: f64) -> anyhow::Result<()> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_string(w: &mut impl Write, value: &str) -> anyhow::Result<()> {
    let bytes = value.as_bytes();
    let len = u32::try_from(bytes.len())
        .with_context(|| format!("string too long ({} bytes)", bytes.len()))?;
    write_u32(w, len)?;
    w.write_all(bytes)?;
    Ok(())
}

fn write_opt_string(w: &mut impl Write, value: &Option<String>) -> anyhow::Result<()> {
    match value {
        Some(s) => {
            write_u8(w, 1)?;
            write_string(w, s)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn write_vec_strings(w: &mut impl Write, values: &[String]) -> anyhow::Result<()> {
    let len = u32::try_from(values.len())
        .with_context(|| format!("too many strings ({})", values.len()))?;
    write_u32(w, len)?;
    for s in values {
        write_string(w, s)?;
    }
    Ok(())
}

fn write_opt_u32(w: &mut impl Write, value: Option<u32>) -> anyhow::Result<()> {
    match value {
        Some(v) => {
            write_u8(w, 1)?;
            write_u32(w, v)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn write_opt_u64(w: &mut impl Write, value: Option<u64>) -> anyhow::Result<()> {
    match value {
        Some(v) => {
            write_u8(w, 1)?;
            write_u64(w, v)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn write_opt_f32(w: &mut impl Write, value: Option<f32>) -> anyhow::Result<()> {
    match value {
        Some(v) => {
            write_u8(w, 1)?;
            write_f32(w, v)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn write_opt_f64(w: &mut impl Write, value: Option<f64>) -> anyhow::Result<()> {
    match value {
        Some(v) => {
            write_u8(w, 1)?;
            write_f64(w, v)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn write_opt_u8(w: &mut impl Write, value: Option<u8>) -> anyhow::Result<()> {
    match value {
        Some(v) => {
            write_u8(w, 1)?;
            write_u8(w, v)?;
        }
        None => write_u8(w, 0)?,
    }
    Ok(())
}

fn read_u8(r: &mut impl Read) -> anyhow::Result<u8> {
    let mut buf = [0u8; 1];
    r.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u32(r: &mut impl Read) -> anyhow::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(r: &mut impl Read) -> anyhow::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_f32(r: &mut impl Read) -> anyhow::Result<f32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(f32::from_le_bytes(buf))
}

fn read_f64(r: &mut impl Read) -> anyhow::Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

fn read_string(r: &mut impl Read) -> anyhow::Result<String> {
    let len = read_u32(r)? as usize;
    if len > MAX_CACHE_STRING_LEN {
        anyhow::bail!("refusing to allocate string of {len} bytes from cache");
    }
    let mut buf = vec![0u8; len];
    r.read_exact(&mut buf)?;
    String::from_utf8(buf).context("cache string was not valid UTF-8")
}

fn read_opt_string(r: &mut impl Read) -> anyhow::Result<Option<String>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_string(r)?)),
        other => anyhow::bail!("invalid option tag {other} for string"),
    }
}

fn read_vec_strings(r: &mut impl Read) -> anyhow::Result<Vec<String>> {
    let len = read_u32(r)? as usize;
    let mut out = Vec::with_capacity(len.min(16));
    for _ in 0..len {
        out.push(read_string(r)?);
    }
    Ok(out)
}

fn read_opt_u32(r: &mut impl Read) -> anyhow::Result<Option<u32>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_u32(r)?)),
        other => anyhow::bail!("invalid option tag {other} for u32"),
    }
}

fn read_opt_u64(r: &mut impl Read) -> anyhow::Result<Option<u64>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_u64(r)?)),
        other => anyhow::bail!("invalid option tag {other} for u64"),
    }
}

fn read_opt_f32(r: &mut impl Read) -> anyhow::Result<Option<f32>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_f32(r)?)),
        other => anyhow::bail!("invalid option tag {other} for f32"),
    }
}

fn read_opt_f64(r: &mut impl Read) -> anyhow::Result<Option<f64>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_f64(r)?)),
        other => anyhow::bail!("invalid option tag {other} for f64"),
    }
}

fn read_opt_u8(r: &mut impl Read) -> anyhow::Result<Option<u8>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(read_u8(r)?)),
        other => anyhow::bail!("invalid option tag {other} for u8"),
    }
}

fn write_mzml_file_info(w: &mut impl Write, info: &MzmlFileInfo) -> anyhow::Result<()> {
    write_opt_string(w, &info.run_id)?;
    write_opt_string(w, &info.start_time)?;
    write_opt_u32(w, info.default_instrument_id)?;
    write_opt_u64(w, info.spectrum_count_hint)?;
    write_vec_strings(w, &info.contents)?;
    write_vec_strings(w, &info.source_files)?;
    write_vec_strings(w, &info.instrument_summaries)?;
    write_vec_strings(w, &info.software_summaries)?;
    Ok(())
}

fn read_mzml_file_info(r: &mut impl Read) -> anyhow::Result<MzmlFileInfo> {
    Ok(MzmlFileInfo {
        run_id: read_opt_string(r)?,
        start_time: read_opt_string(r)?,
        default_instrument_id: read_opt_u32(r)?,
        spectrum_count_hint: read_opt_u64(r)?,
        contents: read_vec_strings(r)?,
        source_files: read_vec_strings(r)?,
        instrument_summaries: read_vec_strings(r)?,
        software_summaries: read_vec_strings(r)?,
    })
}

fn write_index_cache(
    cache_path: &Path,
    fingerprint: FileFingerprint,
    file_info: &MzmlFileInfo,
    metas: &[SpectrumMeta],
) -> anyhow::Result<()> {
    let Some(parent) = cache_path.parent() else {
        anyhow::bail!(
            "cache path {} has no parent directory",
            cache_path.display()
        );
    };
    if !parent.as_os_str().is_empty() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create cache directory {}", parent.display()))?;
    }

    let tmp_path = {
        let file_name = cache_path
            .file_name()
            .unwrap_or_else(|| std::ffi::OsStr::new("mzml-cache"));
        let mut name = file_name.to_os_string();
        name.push(".tmp");
        cache_path.with_file_name(name)
    };

    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&tmp_path)
        .with_context(|| format!("failed to create cache file {}", tmp_path.display()))?;

    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(INDEX_CACHE_MAGIC)?;
    write_u32(&mut encoder, INDEX_CACHE_VERSION)?;
    write_u64(&mut encoder, fingerprint.len)?;
    write_u64(&mut encoder, fingerprint.modified_secs)?;
    write_u32(&mut encoder, fingerprint.modified_nanos)?;
    write_u64(&mut encoder, fingerprint.signature)?;
    write_mzml_file_info(&mut encoder, file_info)?;

    let count = u32::try_from(metas.len())
        .with_context(|| format!("too many spectra to cache ({})", metas.len()))?;
    write_u32(&mut encoder, count)?;
    for meta in metas {
        write_u64(&mut encoder, meta.offset)?;
        write_string(&mut encoder, &meta.scan_id)?;
        write_u8(&mut encoder, meta.ms_level)?;
        write_u8(&mut encoder, meta.continuity as u8)?;
        write_opt_f32(&mut encoder, meta.rt_minutes)?;
        write_opt_f32(&mut encoder, meta.tic)?;
        write_opt_f64(&mut encoder, meta.precursor_mz)?;
        write_opt_u8(&mut encoder, meta.charge)?;
        write_opt_f64(&mut encoder, meta.base_peak_mz)?;
        write_opt_f32(&mut encoder, meta.base_peak_intensity)?;
    }

    let _file = encoder.finish()?;

    if fs::rename(&tmp_path, cache_path).is_err() {
        let _ = fs::remove_file(cache_path);
        fs::rename(&tmp_path, cache_path).with_context(|| {
            format!(
                "failed to move cache file {} into place at {}",
                tmp_path.display(),
                cache_path.display()
            )
        })?;
    }

    Ok(())
}

fn read_index_cache(
    cache_path: &Path,
    expected: &FileFingerprint,
) -> anyhow::Result<(Vec<SpectrumMeta>, MzmlFileInfo, u32)> {
    let file = fs::File::open(cache_path)
        .with_context(|| format!("failed to open cache file {}", cache_path.display()))?;
    let mut decoder = GzDecoder::new(file);

    let mut magic = [0u8; 8];
    decoder.read_exact(&mut magic)?;
    if &magic != INDEX_CACHE_MAGIC {
        anyhow::bail!("invalid cache magic in {}", cache_path.display());
    }
    let version = read_u32(&mut decoder)?;
    match version {
        1 => {
            let fingerprint = FileFingerprint {
                len: read_u64(&mut decoder)?,
                modified_secs: read_u64(&mut decoder)?,
                modified_nanos: read_u32(&mut decoder)?,
                signature: 0,
            };
            if !fingerprint.matches_len_mtime(expected) {
                anyhow::bail!(
                    "cache fingerprint mismatch (file len={} mtime={}.{}; cache len={} mtime={}.{})",
                    expected.len,
                    expected.modified_secs,
                    expected.modified_nanos,
                    fingerprint.len,
                    fingerprint.modified_secs,
                    fingerprint.modified_nanos
                );
            }
        }
        2 | 3 => {
            let fingerprint = FileFingerprint {
                len: read_u64(&mut decoder)?,
                modified_secs: read_u64(&mut decoder)?,
                modified_nanos: read_u32(&mut decoder)?,
                signature: read_u64(&mut decoder)?,
            };
            if fingerprint != *expected {
                anyhow::bail!(
                    "cache fingerprint mismatch (file len={} mtime={}.{} sig={:#016x}; cache len={} mtime={}.{} sig={:#016x})",
                    expected.len,
                    expected.modified_secs,
                    expected.modified_nanos,
                    expected.signature,
                    fingerprint.len,
                    fingerprint.modified_secs,
                    fingerprint.modified_nanos,
                    fingerprint.signature
                );
            }
        }
        other => {
            anyhow::bail!("unsupported cache version {other} (expected {INDEX_CACHE_VERSION})");
        }
    }

    let file_info = read_mzml_file_info(&mut decoder)?;
    let count = read_u32(&mut decoder)? as usize;
    let mut metas = Vec::with_capacity(count.min(50_000));
    for idx in 0..count {
        let offset = read_u64(&mut decoder)?;
        let scan_id = read_string(&mut decoder)?;
        let ms_level = read_u8(&mut decoder)?;
        let continuity_raw = read_u8(&mut decoder)?;
        let continuity = match continuity_raw {
            0 => SignalContinuity::Unknown,
            3 => SignalContinuity::Centroid,
            5 => SignalContinuity::Profile,
            _ => SignalContinuity::Unknown,
        };
        let rt_minutes = read_opt_f32(&mut decoder)?;
        let tic = if version >= 3 {
            read_opt_f32(&mut decoder)?
        } else {
            None
        };
        let precursor_mz = read_opt_f64(&mut decoder)?;
        let charge = read_opt_u8(&mut decoder)?;
        let base_peak_mz = if version >= 3 {
            read_opt_f64(&mut decoder)?
        } else {
            None
        };
        let base_peak_intensity = if version >= 3 {
            read_opt_f32(&mut decoder)?
        } else {
            None
        };

        metas.push(SpectrumMeta {
            idx: idx as u32,
            offset,
            scan_id,
            ms_level,
            rt_minutes,
            tic,
            precursor_mz,
            charge,
            continuity,
            points: None,
            mz_min: None,
            mz_max: None,
            base_peak_mz,
            base_peak_intensity,
        });
    }

    Ok((metas, file_info, version))
}

fn looks_like_spectrum_tag(buf: &[u8]) -> bool {
    let mut i = 0usize;
    while i < buf.len() && buf[i].is_ascii_whitespace() {
        i = i.saturating_add(1);
    }
    buf.get(i..)
        .is_some_and(|tail| tail.starts_with(b"<spectrum"))
}

fn validate_cached_offsets(path: &Path, metas: &[SpectrumMeta], len: u64) -> anyhow::Result<()> {
    if metas.is_empty() {
        anyhow::bail!("cache contained 0 spectra");
    }

    let mut indices = Vec::with_capacity(3);
    indices.push(0usize);
    indices.push(metas.len() / 2);
    indices.push(metas.len().saturating_sub(1));
    indices.sort_unstable();
    indices.dedup();

    let mut file = fs::File::open(path).with_context(|| {
        format!(
            "failed to open mzML for cache validation: {}",
            path.display()
        )
    })?;
    let mut buf = [0u8; 96];

    for &i in &indices {
        let offset = metas
            .get(i)
            .map(|m| m.offset)
            .ok_or_else(|| anyhow::anyhow!("cache meta index {i} missing"))?;
        if offset == 0 || offset >= len {
            anyhow::bail!("cache offset {offset} for meta[{i}] out of range (len={len})");
        }
        file.seek(SeekFrom::Start(offset))?;
        let n = file.read(&mut buf)?;
        if n == 0 || !looks_like_spectrum_tag(&buf[..n]) {
            anyhow::bail!("cache offset {offset} for meta[{i}] did not point to <spectrum>");
        }
    }

    Ok(())
}

pub(super) fn index_mzml(
    path: &Path,
    cache: &IndexCacheOptions,
    mut status: impl FnMut(String),
) -> anyhow::Result<IndexResult> {
    let fingerprint = FileFingerprint::from_path(path)?;
    let cache_path = resolve_index_cache_path(path, cache);
    let mut cache_warning = None::<String>;

    if let Some(cache_path) = cache_path.as_ref() {
        if cache.enabled && !cache.refresh && cache_path.exists() {
            status(format!(
                "Loading cached index... ({})",
                cache_path.display()
            ));
            match read_index_cache(cache_path.as_path(), &fingerprint) {
                Ok((metas, file_info, version)) => {
                    match validate_cached_offsets(path, &metas, fingerprint.len) {
                        Ok(()) => {
                            if version < INDEX_CACHE_VERSION {
                                status(format!(
                                    "Cached index outdated (v{version}); rebuilding..."
                                ));
                            } else {
                                return Ok(IndexResult {
                                    metas,
                                    file_info,
                                    cached: true,
                                    cache_path: Some(cache_path.clone()),
                                    cache_warning,
                                });
                            }
                        }
                        Err(err) => {
                            spectra_debug_log(format!(
                                "cache validation failed ({}): {err:?}",
                                cache_path.display()
                            ));
                            status("Cached index invalid; rebuilding...".to_string());
                        }
                    }
                }
                Err(err) => {
                    spectra_debug_log(format!(
                        "cache load failed ({}): {err:?}",
                        cache_path.display()
                    ));
                    status("Cached index invalid; rebuilding...".to_string());
                }
            }
        }
    }

    let mut reader = MzMLReader::open_path(path)
        .with_context(|| format!("failed to open mzML at {}", path.display()))?;
    reader.detail_level = DetailLevel::MetadataOnly;
    let total = reader.len();
    let file_info = mzml_file_info(&reader);
    let offsets: Vec<u64> = reader
        .get_index()
        .iter()
        .map(|(_id, offset)| *offset)
        .collect();

    let mut metas = Vec::with_capacity(50_000);
    for (i, spec) in reader.into_iter().enumerate() {
        if i % 5000 == 0 && total > 0 {
            status(format!("Indexing spectra... {}/{}", i, total));
        }

        let rt_minutes = {
            let t = spec.start_time() as f32;
            if t > 0.0 {
                Some(t)
            } else {
                None
            }
        };

        let (precursor_mz, charge) = selected_precursor_details(&spec);

        let ms_level = spec.ms_level();
        let continuity = spec.signal_continuity();
        let tic = spec
            .description
            .get_param_by_accession("MS:1000285")
            .and_then(|p| p.parse::<f32>().ok())
            .filter(|v| v.is_finite() && *v >= 0.0);
        let base_peak_mz = spec
            .description
            .get_param_by_accession("MS:1000504")
            .and_then(|p| p.parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v > 0.0);
        let base_peak_intensity = spec
            .description
            .get_param_by_accession("MS:1000505")
            .and_then(|p| p.parse::<f32>().ok())
            .filter(|v| v.is_finite() && *v >= 0.0);

        metas.push(SpectrumMeta {
            idx: i as u32,
            offset: offsets.get(i).copied().unwrap_or(0),
            scan_id: spec.id().to_string(),
            ms_level,
            rt_minutes,
            tic,
            precursor_mz,
            charge,
            continuity,
            points: None,
            mz_min: None,
            mz_max: None,
            base_peak_mz,
            base_peak_intensity,
        });
    }

    if metas.is_empty() && total > 0 {
        spectra_debug_log(format!(
            "indexer fallback: mzdata iterator returned 0 spectra (expected {total}); retrying with random access"
        ));
        let mut reader = MzMLReader::open_path(path)
            .with_context(|| format!("failed to reopen mzML at {}", path.display()))?;
        reader.detail_level = DetailLevel::MetadataOnly;
        metas = Vec::with_capacity(total.min(50_000));
        let offsets: Vec<u64> = reader
            .get_index()
            .iter()
            .map(|(_id, offset)| *offset)
            .collect();

        for i in 0..total {
            if i % 5000 == 0 && total > 0 {
                status(format!("Indexing spectra... {}/{}", i, total));
            }

            let Some(spec) = reader.get_spectrum_by_index(i) else {
                continue;
            };

            let rt_minutes = {
                let t = spec.start_time() as f32;
                if t > 0.0 {
                    Some(t)
                } else {
                    None
                }
            };

            let (precursor_mz, charge) = selected_precursor_details(&spec);
            let ms_level = spec.ms_level();
            let continuity = spec.signal_continuity();
            let tic = spec
                .description
                .get_param_by_accession("MS:1000285")
                .and_then(|p| p.parse::<f32>().ok())
                .filter(|v| v.is_finite() && *v >= 0.0);
            let base_peak_mz = spec
                .description
                .get_param_by_accession("MS:1000504")
                .and_then(|p| p.parse::<f64>().ok())
                .filter(|v| v.is_finite() && *v > 0.0);
            let base_peak_intensity = spec
                .description
                .get_param_by_accession("MS:1000505")
                .and_then(|p| p.parse::<f32>().ok())
                .filter(|v| v.is_finite() && *v >= 0.0);

            metas.push(SpectrumMeta {
                idx: i as u32,
                offset: offsets.get(i).copied().unwrap_or(0),
                scan_id: spec.id().to_string(),
                ms_level,
                rt_minutes,
                tic,
                precursor_mz,
                charge,
                continuity,
                points: None,
                mz_min: None,
                mz_max: None,
                base_peak_mz,
                base_peak_intensity,
            });
        }
    }

    if let Some(cache_path) = cache_path.as_ref() {
        if cache.enabled {
            status(format!("Saving index cache... ({})", cache_path.display()));
            if let Err(err) =
                write_index_cache(cache_path.as_path(), fingerprint, &file_info, &metas)
            {
                spectra_debug_log(format!(
                    "cache write failed ({}): {err:?}",
                    cache_path.display()
                ));
                cache_warning = Some(cache_write_warning(&err));
            }
        }
    }

    Ok(IndexResult {
        metas,
        file_info,
        cached: false,
        cache_path,
        cache_warning,
    })
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

    if spectra_debug_target(idx) {
        spectra_debug_log(format!(
            "spectrum {idx}: arrays mz[dtype={:?},compression={:?},bytes={}] intensity[dtype={:?},compression={:?},bytes={}]",
            mz_array.dtype,
            mz_array.compression,
            mz_bytes.len(),
            intensity_array.dtype,
            intensity_array.compression,
            intensity_bytes.len()
        ));
    }

    let mz = decode_bytes_as_f64(&mz_bytes, mz_array.dtype, idx, "m/z")?;
    let intensity = decode_bytes_as_f32(&intensity_bytes, intensity_array.dtype, idx, "intensity")?;

    Ok((mz, intensity))
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
        mz_min = 0.0;
        mz_max = 0.0;
        base_peak_mz = 0.0;
        base_peak_intensity = 0.0;
    } else {
        if !mz_min.is_finite() {
            mz_min = 0.0;
        }
        if !mz_max.is_finite() {
            mz_max = 0.0;
        }
        if !base_peak_intensity.is_finite() {
            base_peak_intensity = 0.0;
        }
    }

    SpectrumStats {
        points,
        mz_min,
        mz_max,
        base_peak_mz,
        base_peak_intensity,
    }
}

pub(super) fn load_spectrum(
    reader: &mut MzMLReader<std::fs::File>,
    idx: u32,
) -> anyhow::Result<SpectrumOwned> {
    let spec = reader
        .get_spectrum_by_index(idx as usize)
        .ok_or_else(|| anyhow::anyhow!("index {idx} out of range"))?;

    // Prefer raw arrays if present, falling back to peak sets.
    if let Some(arrays) = spec.raw_arrays() {
        match decode_mz_intensity_from_arrays(arrays, idx) {
            Ok((mz, intensity)) => {
                if spectra_debug_target(idx) {
                    spectra_debug_log(format!(
                        "decoded spectrum {idx}: mz.len={} intensity.len={}",
                        mz.len(),
                        intensity.len()
                    ));
                    let mz_head = mz.len().min(5);
                    let inten_head = intensity.len().min(5);
                    if mz_head > 0 {
                        spectra_debug_log(format!(
                            "mz head={:?} tail={:?}",
                            &mz[..mz_head],
                            &mz[mz.len().saturating_sub(mz_head)..]
                        ));
                    }
                    if inten_head > 0 {
                        spectra_debug_log(format!(
                            "intensity head={:?} tail={:?}",
                            &intensity[..inten_head],
                            &intensity[intensity.len().saturating_sub(inten_head)..]
                        ));
                    }

                    let mut mz_min = f64::INFINITY;
                    let mut mz_max = -f64::INFINITY;
                    let mut mz_nonfinite = 0usize;
                    for &v in mz.iter() {
                        if v.is_finite() {
                            mz_min = mz_min.min(v);
                            mz_max = mz_max.max(v);
                        } else {
                            mz_nonfinite = mz_nonfinite.saturating_add(1);
                        }
                    }
                    let mut i_min = f32::INFINITY;
                    let mut i_max = -f32::INFINITY;
                    let mut i_nonfinite = 0usize;
                    for &v in intensity.iter() {
                        if v.is_finite() {
                            i_min = i_min.min(v);
                            i_max = i_max.max(v);
                        } else {
                            i_nonfinite = i_nonfinite.saturating_add(1);
                        }
                    }

                    spectra_debug_log(format!(
                        "mz[min,max]=[{mz_min},{mz_max}] nonfinite={mz_nonfinite} intensity[min,max]=[{i_min},{i_max}] nonfinite={i_nonfinite}"
                    ));
                }

                if mz.len() != intensity.len() {
                    spectra_debug_log(format!(
                        "spectrum {idx}: length mismatch mz.len={} intensity.len={}",
                        mz.len(),
                        intensity.len()
                    ));
                    debug_assert_eq!(
                        mz.len(),
                        intensity.len(),
                        "spectrum {idx}: mz/intensity length mismatch"
                    );
                } else {
                    #[cfg(debug_assertions)]
                    {
                        let mut decreases = 0usize;
                        let mut comparisons = 0usize;
                        for w in mz.windows(2) {
                            comparisons = comparisons.saturating_add(1);
                            if w[1] < w[0] {
                                decreases = decreases.saturating_add(1);
                            }
                        }
                        if comparisons > 0 {
                            let frac = decreases as f64 / comparisons as f64;
                            debug_assert!(
                                frac <= 0.05,
                                "spectrum {idx}: mz not mostly sorted (decreases {decreases}/{comparisons})"
                            );
                        }
                    }

                    let stats = compute_stats(&mz, &intensity);
                    return Ok(SpectrumOwned {
                        mz: Arc::from(mz),
                        intensity: Arc::from(intensity),
                        stats,
                    });
                }
            }
            Err(err) => {
                spectra_debug_log(format!("spectrum {idx}: raw array decode failed: {err:?}"));
            }
        }
    }

    let peaks = spec.peaks();
    let mut mz = Vec::<f64>::new();
    let mut intensity = Vec::<f32>::new();

    match peaks {
        RefPeakDataLevel::Missing => {}
        RefPeakDataLevel::RawData(arrays) => {
            if let Ok((mzs, intensities)) = decode_mz_intensity_from_arrays(arrays, idx) {
                mz = mzs;
                intensity = intensities;
            };
        }
        RefPeakDataLevel::Centroid(peaks) => {
            mz.reserve(peaks.len());
            intensity.reserve(peaks.len());
            for peak in peaks.iter() {
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

    if mz.len() != intensity.len() {
        spectra_debug_log(format!(
            "spectrum {idx}: peak length mismatch mz.len={} intensity.len={}",
            mz.len(),
            intensity.len()
        ));
        debug_assert_eq!(
            mz.len(),
            intensity.len(),
            "spectrum {idx}: peak length mismatch"
        );
        anyhow::bail!(
            "spectrum {idx}: mz/intensity length mismatch (mz.len={} intensity.len={})",
            mz.len(),
            intensity.len()
        );
    }
    let stats = compute_stats(&mz, &intensity);

    Ok(SpectrumOwned {
        mz: Arc::from(mz),
        intensity: Arc::from(intensity),
        stats,
    })
}
