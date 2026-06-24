use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Context;
use flate2::read::MultiGzDecoder;
use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;
use serde::Serialize;

use crate::annotate::{residue_mass, ExplicitModification, ModificationSite};
use crate::mzml::extract_scan_number;

#[derive(Clone, Debug)]
pub(crate) struct PepXmlScanHits {
    pub(crate) requested_top_n: usize,
    pub(crate) available_hits: usize,
    pub(crate) hits: Vec<PepXmlHit>,
}

#[derive(Clone, Debug)]
pub(crate) struct PepXmlHit {
    pub(crate) hit_rank: usize,
    pub(crate) peptide: String,
    pub(crate) assumed_charge: Option<i32>,
    pub(crate) spectrum: Option<String>,
    pub(crate) start_scan: Option<u64>,
    pub(crate) end_scan: Option<u64>,
    pub(crate) protein: Option<String>,
    pub(crate) calc_neutral_pep_mass: Option<f64>,
    pub(crate) massdiff: Option<f64>,
    pub(crate) modifications: Vec<PepXmlModification>,
    pub(crate) scores: Vec<PepXmlScore>,
}

impl PepXmlHit {
    pub(crate) fn explicit_modifications(&self) -> Vec<ExplicitModification> {
        self.modifications
            .iter()
            .map(|modification| ExplicitModification {
                site: match modification.site {
                    PepXmlModificationSite::NTerm => ModificationSite::NTerm,
                    PepXmlModificationSite::CTerm => ModificationSite::CTerm,
                    PepXmlModificationSite::Residue(position) => {
                        ModificationSite::Residue(position)
                    }
                },
                delta: modification.delta,
            })
            .collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub(crate) enum PepXmlModificationSite {
    NTerm,
    CTerm,
    Residue(usize),
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct PepXmlModification {
    pub(crate) site: PepXmlModificationSite,
    pub(crate) residue: Option<char>,
    pub(crate) reported_mass: f64,
    pub(crate) base_mass: Option<f64>,
    pub(crate) delta: f64,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct PepXmlScore {
    pub(crate) name: String,
    pub(crate) value: String,
}

#[derive(Clone, Debug)]
struct QueryContext {
    matches_scan: bool,
    spectrum: Option<String>,
    start_scan: Option<u64>,
    end_scan: Option<u64>,
    assumed_charge: Option<i32>,
}

#[derive(Clone, Debug)]
struct HitBuilder {
    file_order: usize,
    hit_rank: Option<usize>,
    peptide: Option<String>,
    assumed_charge: Option<i32>,
    spectrum: Option<String>,
    start_scan: Option<u64>,
    end_scan: Option<u64>,
    protein: Option<String>,
    calc_neutral_pep_mass: Option<f64>,
    massdiff: Option<f64>,
    modifications: Vec<PepXmlModification>,
    scores: Vec<PepXmlScore>,
}

pub(crate) fn load_hits_for_scan(
    path: &Path,
    scan: u64,
    top_n: usize,
) -> anyhow::Result<PepXmlScanHits> {
    if top_n == 0 {
        anyhow::bail!("--top-n must be at least 1");
    }
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
    {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    parse_hits_for_scan(reader, scan, top_n)
        .with_context(|| format!("failed to parse {}", path.display()))
}

fn parse_hits_for_scan<R: BufRead>(
    reader: R,
    scan: u64,
    top_n: usize,
) -> anyhow::Result<PepXmlScanHits> {
    let mut xml = Reader::from_reader(reader);
    xml.trim_text(true);
    let mut buf = Vec::new();
    let mut query = None::<QueryContext>;
    let mut current_hit = None::<HitBuilder>;
    let mut hits = Vec::<PepXmlHit>::new();
    let mut file_order = 0_usize;

    loop {
        match xml.read_event_into(&mut buf)? {
            Event::Start(start) => {
                handle_start(&start, scan, &mut query, &mut current_hit, &mut file_order)?;
            }
            Event::Empty(start) => {
                handle_start(&start, scan, &mut query, &mut current_hit, &mut file_order)?;
                handle_end(
                    start.local_name().as_ref(),
                    &mut query,
                    &mut current_hit,
                    &mut hits,
                )?;
            }
            Event::End(end) => {
                handle_end(
                    end.local_name().as_ref(),
                    &mut query,
                    &mut current_hit,
                    &mut hits,
                )?;
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    hits.sort_by_key(|hit| hit.hit_rank);
    let available_hits = hits.len();
    let hits = hits.into_iter().take(top_n).collect::<Vec<_>>();
    Ok(PepXmlScanHits {
        requested_top_n: top_n,
        available_hits,
        hits,
    })
}

fn handle_start(
    start: &BytesStart<'_>,
    scan: u64,
    query: &mut Option<QueryContext>,
    current_hit: &mut Option<HitBuilder>,
    file_order: &mut usize,
) -> anyhow::Result<()> {
    match start.local_name().as_ref() {
        b"spectrum_query" => {
            let attrs = attributes(start)?;
            let spectrum = attr_value(&attrs, "spectrum");
            let start_scan = parse_optional_u64(attr_value(&attrs, "start_scan"))?;
            let end_scan = parse_optional_u64(attr_value(&attrs, "end_scan"))?;
            let assumed_charge = parse_optional_i32(attr_value(&attrs, "assumed_charge"))?;
            let matches_scan = query_matches_scan(scan, spectrum.as_deref(), start_scan, end_scan);
            *query = Some(QueryContext {
                matches_scan,
                spectrum,
                start_scan,
                end_scan,
                assumed_charge,
            });
        }
        b"search_hit" => {
            if query.as_ref().is_some_and(|ctx| ctx.matches_scan) {
                let attrs = attributes(start)?;
                let ctx = query
                    .as_ref()
                    .expect("query is checked above for matching scan");
                *file_order += 1;
                *current_hit = Some(HitBuilder {
                    file_order: *file_order,
                    hit_rank: parse_optional_usize(attr_value(&attrs, "hit_rank"))?,
                    peptide: attr_value(&attrs, "peptide"),
                    assumed_charge: ctx.assumed_charge,
                    spectrum: ctx.spectrum.clone(),
                    start_scan: ctx.start_scan,
                    end_scan: ctx.end_scan,
                    protein: attr_value(&attrs, "protein"),
                    calc_neutral_pep_mass: parse_optional_f64(attr_value(
                        &attrs,
                        "calc_neutral_pep_mass",
                    ))?,
                    massdiff: parse_optional_f64(attr_value(&attrs, "massdiff"))?,
                    modifications: Vec::new(),
                    scores: Vec::new(),
                });
            }
        }
        b"modification_info" => {
            if let Some(hit) = current_hit.as_mut() {
                let attrs = attributes(start)?;
                if let Some(value) = parse_optional_f64(attr_value(&attrs, "mod_nterm_mass"))? {
                    hit.modifications.push(PepXmlModification {
                        site: PepXmlModificationSite::NTerm,
                        residue: None,
                        reported_mass: value,
                        base_mass: None,
                        delta: value,
                    });
                }
                if let Some(value) = parse_optional_f64(attr_value(&attrs, "mod_cterm_mass"))? {
                    hit.modifications.push(PepXmlModification {
                        site: PepXmlModificationSite::CTerm,
                        residue: None,
                        reported_mass: value,
                        base_mass: None,
                        delta: value,
                    });
                }
            }
        }
        b"mod_aminoacid_mass" => {
            if let Some(hit) = current_hit.as_mut() {
                let attrs = attributes(start)?;
                let position = parse_required_usize(attr_value(&attrs, "position"), "position")?;
                let reported_mass = parse_required_f64(attr_value(&attrs, "mass"), "mass")?;
                let peptide = hit
                    .peptide
                    .as_deref()
                    .context("pepXML search_hit missing peptide before mod_aminoacid_mass")?;
                let residue = peptide
                    .chars()
                    .nth(position.saturating_sub(1))
                    .ok_or_else(|| {
                        anyhow::anyhow!(
                            "pepXML modification position {} is out of range for peptide {}",
                            position,
                            peptide
                        )
                    })?;
                let base_mass = residue_mass(residue).ok_or_else(|| {
                    anyhow::anyhow!(
                        "unsupported amino-acid code `{}` in pepXML peptide {}",
                        residue,
                        peptide
                    )
                })?;
                hit.modifications.push(PepXmlModification {
                    site: PepXmlModificationSite::Residue(position),
                    residue: Some(residue),
                    reported_mass,
                    base_mass: Some(base_mass),
                    delta: reported_mass - base_mass,
                });
            }
        }
        b"search_score" => {
            if let Some(hit) = current_hit.as_mut() {
                let attrs = attributes(start)?;
                if let Some(name) = attr_value(&attrs, "name") {
                    let value = attr_value(&attrs, "value").unwrap_or_default();
                    hit.scores.push(PepXmlScore { name, value });
                }
            }
        }
        _ => {}
    }
    Ok(())
}

fn handle_end(
    local_name: &[u8],
    query: &mut Option<QueryContext>,
    current_hit: &mut Option<HitBuilder>,
    hits: &mut Vec<PepXmlHit>,
) -> anyhow::Result<()> {
    match local_name {
        b"search_hit" => {
            if let Some(builder) = current_hit.take() {
                hits.push(builder.finish()?);
            }
        }
        b"spectrum_query" => {
            *query = None;
        }
        _ => {}
    }
    Ok(())
}

impl HitBuilder {
    fn finish(self) -> anyhow::Result<PepXmlHit> {
        let peptide = self
            .peptide
            .ok_or_else(|| anyhow::anyhow!("pepXML search_hit missing peptide"))?;
        Ok(PepXmlHit {
            hit_rank: self.hit_rank.unwrap_or(self.file_order),
            peptide,
            assumed_charge: self.assumed_charge,
            spectrum: self.spectrum,
            start_scan: self.start_scan,
            end_scan: self.end_scan,
            protein: self.protein,
            calc_neutral_pep_mass: self.calc_neutral_pep_mass,
            massdiff: self.massdiff,
            modifications: self.modifications,
            scores: self.scores,
        })
    }
}

fn attributes(start: &BytesStart<'_>) -> anyhow::Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for attr in start.attributes().with_checks(false) {
        let attr = attr?;
        let key = std::str::from_utf8(attr.key.as_ref())?.to_string();
        let value = attr.unescape_value()?.into_owned();
        out.push((key, value));
    }
    Ok(out)
}

fn attr_value(attrs: &[(String, String)], key: &str) -> Option<String> {
    attrs
        .iter()
        .find(|(candidate, _)| candidate == key)
        .map(|(_, value)| value.clone())
}

fn parse_optional_usize(raw: Option<String>) -> anyhow::Result<Option<usize>> {
    raw.map(|value| {
        value
            .parse::<usize>()
            .with_context(|| format!("invalid pepXML integer `{value}`"))
    })
    .transpose()
}

fn parse_optional_u64(raw: Option<String>) -> anyhow::Result<Option<u64>> {
    raw.map(|value| {
        value
            .parse::<u64>()
            .with_context(|| format!("invalid pepXML scan number `{value}`"))
    })
    .transpose()
}

fn parse_optional_i32(raw: Option<String>) -> anyhow::Result<Option<i32>> {
    raw.map(|value| {
        value
            .parse::<i32>()
            .with_context(|| format!("invalid pepXML charge `{value}`"))
    })
    .transpose()
}

fn parse_optional_f64(raw: Option<String>) -> anyhow::Result<Option<f64>> {
    raw.map(|value| {
        value
            .parse::<f64>()
            .with_context(|| format!("invalid pepXML float `{value}`"))
    })
    .transpose()
}

fn parse_required_usize(raw: Option<String>, attr: &str) -> anyhow::Result<usize> {
    raw.ok_or_else(|| anyhow::anyhow!("pepXML missing required `{attr}` attribute"))?
        .parse::<usize>()
        .with_context(|| format!("invalid pepXML `{attr}` integer"))
}

fn parse_required_f64(raw: Option<String>, attr: &str) -> anyhow::Result<f64> {
    raw.ok_or_else(|| anyhow::anyhow!("pepXML missing required `{attr}` attribute"))?
        .parse::<f64>()
        .with_context(|| format!("invalid pepXML `{attr}` float"))
}

fn query_matches_scan(
    scan: u64,
    spectrum: Option<&str>,
    start_scan: Option<u64>,
    end_scan: Option<u64>,
) -> bool {
    if let Some(start) = start_scan {
        let end = end_scan.unwrap_or(start);
        if start <= scan && scan <= end {
            return true;
        }
    }
    spectrum.is_some_and(|value| spectrum_string_matches_scan(value, scan))
}

fn spectrum_string_matches_scan(spectrum: &str, scan: u64) -> bool {
    if extract_scan_number(spectrum) == Some(scan) {
        return true;
    }

    let numeric_tokens = spectrum
        .split(|ch: char| !ch.is_ascii_digit())
        .filter(|token| !token.is_empty())
        .filter_map(|token| token.parse::<u64>().ok())
        .collect::<Vec<_>>();
    if numeric_tokens.len() >= 2 {
        numeric_tokens[..numeric_tokens.len() - 1].contains(&scan)
    } else {
        numeric_tokens.contains(&scan)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    const PEPXML: &str = r#"
<msms_pipeline_analysis>
  <msms_run_summary>
    <spectrum_query spectrum="demo.2.2.2" start_scan="2" end_scan="2" assumed_charge="2">
      <search_result>
        <search_hit hit_rank="1" peptide="PEPMUDEK" protein="sp|P1" calc_neutral_pep_mass="1000.1" massdiff="0.01">
          <modification_info mod_nterm_mass="42.0106" mod_cterm_mass="17.0030">
            <mod_aminoacid_mass position="4" mass="147.035384645"/>
            <mod_aminoacid_mass position="5" mass="150.953633405"/>
          </modification_info>
          <search_score name="hyperscore" value="51.2"/>
        </search_hit>
        <search_hit hit_rank="2" peptide="PEPTIDEK" protein="sp|P2">
          <search_score name="hyperscore" value="19.7"/>
        </search_hit>
      </search_result>
    </spectrum_query>
</msms_run_summary>
</msms_pipeline_analysis>
"#;

    #[test]
    fn parses_top_hits_and_converts_residue_masses_to_deltas() {
        let parsed = parse_hits_for_scan(PEPXML.as_bytes(), 2, 3).expect("parse pepXML");
        assert_eq!(parsed.requested_top_n, 3);
        assert_eq!(parsed.available_hits, 2);
        assert_eq!(parsed.hits.len(), 2);

        let top = &parsed.hits[0];
        assert_eq!(top.hit_rank, 1);
        assert_eq!(top.peptide, "PEPMUDEK");
        assert_eq!(top.assumed_charge, Some(2));
        assert_eq!(top.scores[0].name, "hyperscore");
        assert_eq!(top.scores[0].value, "51.2");

        assert_eq!(top.modifications.len(), 4);
        assert!(top.modifications.iter().any(|modification| matches!(
            modification.site,
            PepXmlModificationSite::NTerm
        ) && (modification.delta - 42.0106)
            .abs()
            < 1e-9));
        assert!(top.modifications.iter().any(|modification| matches!(
            modification.site,
            PepXmlModificationSite::CTerm
        ) && (modification.delta - 17.0030)
            .abs()
            < 1e-9));
        assert!(top.modifications.iter().any(|modification| {
            matches!(modification.site, PepXmlModificationSite::Residue(4))
                && modification.residue == Some('M')
                && (modification.delta - 15.9949).abs() < 1e-6
        }));
        assert!(top.modifications.iter().any(|modification| {
            matches!(modification.site, PepXmlModificationSite::Residue(5))
                && modification.residue == Some('U')
                && modification.delta.abs() < 1e-9
        }));
    }

    #[test]
    fn matches_common_spectrum_string_without_scan_attributes() {
        let xml = r#"
<msms_pipeline_analysis>
  <msms_run_summary>
    <spectrum_query spectrum="demo.44533.44533.2" assumed_charge="2">
      <search_result>
        <search_hit hit_rank="1" peptide="PEPTIDEK"/>
      </search_result>
    </spectrum_query>
  </msms_run_summary>
</msms_pipeline_analysis>
"#;
        let parsed = parse_hits_for_scan(xml.as_bytes(), 44533, 1).expect("parse pepXML");
        assert_eq!(parsed.available_hits, 1);
        assert_eq!(parsed.hits[0].peptide, "PEPTIDEK");
    }

    #[test]
    fn loads_gzipped_pepxml_path() {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "mzio-pepxml-test-{}-{timestamp}.pepXML.gz",
            std::process::id()
        ));
        let file = fs::File::create(&path).expect("create gz fixture");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder
            .write_all(PEPXML.as_bytes())
            .expect("write gz fixture");
        encoder.finish().expect("finish gz fixture");

        let parsed = load_hits_for_scan(&path, 2, 1).expect("load gzipped pepXML");
        assert_eq!(parsed.available_hits, 2);
        assert_eq!(parsed.hits.len(), 1);
        assert_eq!(parsed.hits[0].peptide, "PEPMUDEK");

        let _ = fs::remove_file(path);
    }
}
