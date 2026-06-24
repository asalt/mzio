use std::cmp::Ordering;

use anyhow::Context;

const PROTON_MASS: f64 = 1.007_276_466_812;
const WATER_MASS: f64 = 18.010_564_683_7;
const AMMONIA_MASS: f64 = 17.026_549_101;
const PHOSPHORIC_ACID_MASS: f64 = 97.976_895_573;
const NEUTRON_MASS_DIFF: f64 = 1.003_354_835_07;
const PHOSPHO_DELTA_MASS: f64 = 79.966_331;
const PHOSPHO_DELTA_EPSILON: f64 = 0.01;
const QUALITY_MAX_HEAVY_ISOTOPE: usize = 2;

pub(crate) const DEFAULT_PRECURSOR_ISOTOPE_ERRORS: [u8; 3] = [0, 1, 2];

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum MassTolerance {
    Ppm(f64),
    Da(f64),
}

impl MassTolerance {
    pub(crate) fn window_da(self, theoretical_mz: f64) -> f64 {
        match self {
            Self::Ppm(ppm) => theoretical_mz.abs() * ppm * 1e-6,
            Self::Da(da) => da.abs(),
        }
    }

    pub(crate) fn contains(self, theoretical_mz: f64, observed_mz: f64) -> bool {
        (observed_mz - theoretical_mz).abs() <= self.window_da(theoretical_mz)
    }

    pub(crate) fn error_ppm(self, theoretical_mz: f64, observed_mz: f64) -> f64 {
        if theoretical_mz.abs() <= f64::EPSILON {
            0.0
        } else {
            (observed_mz - theoretical_mz) / theoretical_mz * 1_000_000.0
        }
    }

    pub(crate) fn label(self) -> String {
        match self {
            Self::Ppm(ppm) => format!("{ppm:.1} ppm"),
            Self::Da(da) => format!("{da:.3} Da"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ModificationSite {
    NTerm,
    CTerm,
    Residue(usize),
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ExplicitModification {
    pub(crate) site: ModificationSite,
    pub(crate) delta: f64,
}

impl ExplicitModification {
    pub(crate) fn label(&self) -> String {
        match self.site {
            ModificationSite::NTerm => format!("n-term:{:+.6}", self.delta),
            ModificationSite::CTerm => format!("c-term:{:+.6}", self.delta),
            ModificationSite::Residue(position) => format!("{position}:{:+.6}", self.delta),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum FragmentSeries {
    B,
    Y,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum NeutralLossKind {
    Water,
    Ammonia,
    PhosphoricAcid,
}

impl NeutralLossKind {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Water => "H2O",
            Self::Ammonia => "NH3",
            Self::PhosphoricAcid => "H3PO4",
        }
    }

    fn mass(self) -> f64 {
        match self {
            Self::Water => WATER_MASS,
            Self::Ammonia => AMMONIA_MASS,
            Self::PhosphoricAcid => PHOSPHORIC_ACID_MASS,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct FragmentIon {
    pub(crate) series: FragmentSeries,
    pub(crate) ordinal: usize,
    pub(crate) cleavage_index: usize,
    pub(crate) charge: u8,
    pub(crate) neutral_loss: Option<NeutralLossKind>,
    pub(crate) theoretical_mz: f64,
}

impl FragmentIon {
    pub(crate) fn label(&self) -> String {
        let prefix = match self.series {
            FragmentSeries::B => "b",
            FragmentSeries::Y => "y",
        };
        let charge_suffix = if self.charge <= 1 {
            String::new()
        } else {
            "+".repeat(self.charge as usize)
        };
        let neutral_loss_suffix = self
            .neutral_loss
            .map(|loss| format!("-{}", loss.label()))
            .unwrap_or_default();
        format!(
            "{prefix}{}{neutral_loss_suffix}{charge_suffix}",
            self.ordinal
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct FragmentMatch {
    pub(crate) fragment: FragmentIon,
    pub(crate) peak_index: usize,
    pub(crate) observed_mz: f64,
    pub(crate) observed_intensity: f32,
    pub(crate) error_da: f64,
    pub(crate) error_ppm: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct AnnotationQualityMetrics {
    pub(crate) snr_like: f64,
    pub(crate) log2_snr_like: f64,
    pub(crate) cosine: f64,
    pub(crate) frag_error_mae_ppm: Option<f64>,
    pub(crate) frag_error_mae_da: Option<f64>,
}

#[derive(Clone, Debug)]
pub(crate) struct PrecursorCheck {
    pub(crate) charge: i32,
    pub(crate) monoisotopic_theoretical_mz: f64,
    pub(crate) theoretical_mz: f64,
    pub(crate) observed_mz: f64,
    pub(crate) isotope_error: u8,
    pub(crate) error_da: f64,
    pub(crate) error_ppm: f64,
    pub(crate) within_tolerance: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct Peptide {
    sequence: String,
    residue_chars: Vec<char>,
    residue_deltas: Vec<f64>,
    residue_masses: Vec<f64>,
    n_term_delta: f64,
    c_term_delta: f64,
    charge_hint: Option<i32>,
}

impl Peptide {
    pub(crate) fn sequence(&self) -> &str {
        &self.sequence
    }

    pub(crate) fn len(&self) -> usize {
        self.residue_chars.len()
    }

    pub(crate) fn residue_chars(&self) -> &[char] {
        &self.residue_chars
    }

    pub(crate) fn charge_hint(&self) -> Option<i32> {
        self.charge_hint
    }

    pub(crate) fn neutral_mass(&self) -> f64 {
        self.total_residue_mass() + self.n_term_delta + self.c_term_delta + WATER_MASS
    }

    pub(crate) fn precursor_mz(&self, charge: i32) -> Option<f64> {
        if charge <= 0 {
            return None;
        }
        let z = charge as f64;
        Some((self.neutral_mass() + z * PROTON_MASS) / z)
    }

    fn total_residue_mass(&self) -> f64 {
        self.residue_masses.iter().copied().sum()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct AnnotationContext {
    pub(crate) peptide: Peptide,
    pub(crate) modifications: Vec<ExplicitModification>,
    pub(crate) neutral_losses: Vec<NeutralLossKind>,
    pub(crate) charge_context: Option<i32>,
    pub(crate) isotope_errors: Vec<u8>,
    pub(crate) tolerance: MassTolerance,
}

impl AnnotationContext {
    pub(crate) fn modifications_label(&self) -> Option<String> {
        if self.modifications.is_empty() {
            None
        } else {
            Some(
                self.modifications
                    .iter()
                    .map(ExplicitModification::label)
                    .collect::<Vec<_>>()
                    .join(", "),
            )
        }
    }

    pub(crate) fn neutral_losses_label(&self) -> Option<String> {
        if self.neutral_losses.is_empty() {
            None
        } else {
            Some(
                self.neutral_losses
                    .iter()
                    .map(|loss| loss.label())
                    .collect::<Vec<_>>()
                    .join(", "),
            )
        }
    }

    pub(crate) fn isotope_errors_label(&self) -> Option<String> {
        if self.isotope_errors.is_empty() {
            None
        } else {
            Some(
                self.isotope_errors
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            )
        }
    }

    pub(crate) fn modified_sequence(&self) -> String {
        let mut residue_mods = vec![Vec::<String>::new(); self.peptide.len()];
        let mut n_term_mods = Vec::<String>::new();
        let mut c_term_mods = Vec::<String>::new();
        for modification in &self.modifications {
            match modification.site {
                ModificationSite::NTerm => {
                    n_term_mods.push(format!("[{:+.6}]", modification.delta))
                }
                ModificationSite::CTerm => {
                    c_term_mods.push(format!("[{:+.6}]", modification.delta))
                }
                ModificationSite::Residue(position) => {
                    if position > 0 && position <= residue_mods.len() {
                        residue_mods[position - 1].push(format!("[{:+.6}]", modification.delta));
                    }
                }
            }
        }

        let mut out = String::new();
        for value in n_term_mods {
            out.push_str(&value);
        }
        for (idx, aa) in self.peptide.residue_chars().iter().copied().enumerate() {
            out.push(aa);
            for value in &residue_mods[idx] {
                out.push_str(value);
            }
        }
        for value in c_term_mods {
            out.push_str(&value);
        }
        out
    }
}

#[derive(Clone, Debug)]
pub(crate) struct AnnotationReport {
    pub(crate) context: AnnotationContext,
    pub(crate) fragments: Vec<FragmentIon>,
    pub(crate) matches: Vec<FragmentMatch>,
    pub(crate) precursor_check: Option<PrecursorCheck>,
    pub(crate) quality: AnnotationQualityMetrics,
}

impl AnnotationReport {
    pub(crate) fn matched_peak_count(&self) -> usize {
        let mut peak_ids = self
            .matches
            .iter()
            .map(|m| m.peak_index)
            .collect::<Vec<_>>();
        peak_ids.sort_unstable();
        peak_ids.dedup();
        peak_ids.len()
    }
}

#[derive(Clone, Debug)]
struct ParsedPeptideInput {
    sequence: String,
    inline_modifications: Vec<ExplicitModification>,
    charge_hint: Option<i32>,
}

#[derive(Clone, Debug)]
struct ObservedPeak {
    index: usize,
    mz: f64,
    intensity: f32,
}

pub(crate) fn prepare_annotation(
    peptide_input: &str,
    mod_inputs: &[String],
    neutral_losses: &[NeutralLossKind],
    preferred_charge_context: Option<i32>,
    tolerance: MassTolerance,
) -> anyhow::Result<AnnotationContext> {
    let parsed = parse_peptide_input(peptide_input)?;
    let mut modifications = parsed.inline_modifications;
    modifications.reserve(mod_inputs.len());
    for raw in mod_inputs {
        modifications.push(parse_mod_spec(raw)?);
    }
    prepare_annotation_from_parts(
        parsed.sequence,
        parsed.charge_hint,
        modifications,
        neutral_losses,
        preferred_charge_context,
        tolerance,
    )
}

pub(crate) fn prepare_annotation_with_modifications(
    peptide_input: &str,
    explicit_modifications: Vec<ExplicitModification>,
    neutral_losses: &[NeutralLossKind],
    preferred_charge_context: Option<i32>,
    tolerance: MassTolerance,
) -> anyhow::Result<AnnotationContext> {
    let parsed = parse_peptide_input(peptide_input)?;
    let mut modifications = parsed.inline_modifications;
    modifications.extend(explicit_modifications);
    prepare_annotation_from_parts(
        parsed.sequence,
        parsed.charge_hint,
        modifications,
        neutral_losses,
        preferred_charge_context,
        tolerance,
    )
}

fn prepare_annotation_from_parts(
    sequence: String,
    charge_hint: Option<i32>,
    modifications: Vec<ExplicitModification>,
    neutral_losses: &[NeutralLossKind],
    preferred_charge_context: Option<i32>,
    tolerance: MassTolerance,
) -> anyhow::Result<AnnotationContext> {
    let peptide = build_peptide(&sequence, charge_hint, &modifications)?;
    let charge_context = preferred_charge_context.or(peptide.charge_hint());
    Ok(AnnotationContext {
        peptide,
        modifications,
        neutral_losses: neutral_losses.to_vec(),
        charge_context,
        isotope_errors: DEFAULT_PRECURSOR_ISOTOPE_ERRORS.to_vec(),
        tolerance,
    })
}

pub(crate) fn annotate_peaks(
    context: &AnnotationContext,
    observed_precursor_mz: Option<f64>,
    mz: &[f64],
    intensity: &[f32],
) -> AnnotationReport {
    let fragment_charges = fragment_charge_states(context.charge_context);
    let fragments =
        generate_fragments(&context.peptide, &fragment_charges, &context.neutral_losses);
    let matches = match_fragments(&fragments, mz, intensity, context.tolerance);
    let quality = calculate_quality_metrics(&matches, mz, intensity, context.tolerance);
    let precursor_check = match (context.charge_context, observed_precursor_mz) {
        (Some(charge), Some(observed_mz)) if charge > 0 => context
            .peptide
            .precursor_mz(charge)
            .map(|monoisotopic_theoretical_mz| {
                let z = charge as f64;
                let mut candidates = context
                    .isotope_errors
                    .iter()
                    .copied()
                    .map(|isotope_error| {
                        let theoretical_mz = monoisotopic_theoretical_mz
                            - isotope_error as f64 * NEUTRON_MASS_DIFF / z;
                        let error_da = observed_mz - theoretical_mz;
                        let error_ppm = context.tolerance.error_ppm(theoretical_mz, observed_mz);
                        PrecursorCheck {
                            charge,
                            monoisotopic_theoretical_mz,
                            theoretical_mz,
                            observed_mz,
                            isotope_error,
                            error_da,
                            error_ppm,
                            within_tolerance: context
                                .tolerance
                                .contains(theoretical_mz, observed_mz),
                        }
                    })
                    .collect::<Vec<_>>();
                candidates.sort_by(|left, right| {
                    right
                        .within_tolerance
                        .cmp(&left.within_tolerance)
                        .then_with(|| {
                            left.error_da
                                .abs()
                                .partial_cmp(&right.error_da.abs())
                                .unwrap_or(Ordering::Equal)
                        })
                        .then_with(|| left.isotope_error.cmp(&right.isotope_error))
                });
                candidates
                    .into_iter()
                    .next()
                    .expect("default isotope errors are non-empty")
            }),
        _ => None,
    };

    AnnotationReport {
        context: context.clone(),
        fragments,
        matches,
        precursor_check,
        quality,
    }
}

fn calculate_quality_metrics(
    matches: &[FragmentMatch],
    mz: &[f64],
    intensity: &[f32],
    tolerance: MassTolerance,
) -> AnnotationQualityMetrics {
    const LOG_RATIO_EPSILON: f64 = 1.0;

    let observed = sorted_observed_peaks(mz, intensity);
    let mut matched_peak_indices = matches
        .iter()
        .map(|matched| matched.peak_index)
        .collect::<Vec<_>>();

    for matched in matches {
        if !should_integrate_fragment_isotopes(&matched.fragment, tolerance) {
            continue;
        }
        let isotope_spacing = NEUTRON_MASS_DIFF / matched.fragment.charge as f64;
        for isotope_order in 1..=QUALITY_MAX_HEAVY_ISOTOPE {
            let target_mz =
                matched.fragment.theoretical_mz + isotope_order as f64 * isotope_spacing;
            if let Some(peak) = best_peak_for_target(&observed, target_mz, tolerance) {
                matched_peak_indices.push(peak.index);
            }
        }
    }

    matched_peak_indices.sort_unstable();
    matched_peak_indices.dedup();

    let mut matched_intensity = 0.0;
    let mut unmatched_intensity = 0.0;
    let mut matched_sum_squares = 0.0;
    let mut observed_sum_squares = 0.0;

    for (idx, value) in intensity.iter().copied().enumerate() {
        if !value.is_finite() {
            continue;
        }
        let value = (value as f64).max(0.0);
        observed_sum_squares += value * value;
        if matched_peak_indices.binary_search(&idx).is_ok() {
            matched_intensity += value;
            matched_sum_squares += value * value;
        } else {
            unmatched_intensity += value;
        }
    }

    let snr_like = if unmatched_intensity > 0.0 {
        matched_intensity / unmatched_intensity
    } else if matched_intensity > 0.0 {
        f64::INFINITY
    } else {
        0.0
    };
    let log2_snr_like = ((matched_intensity + LOG_RATIO_EPSILON)
        / (unmatched_intensity + LOG_RATIO_EPSILON))
        .log2();
    let cosine_denominator = matched_sum_squares.sqrt() * observed_sum_squares.sqrt();
    let cosine = if cosine_denominator > 0.0 {
        (matched_sum_squares / cosine_denominator).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (frag_error_mae_ppm, frag_error_mae_da) = if matches.is_empty() {
        (None, None)
    } else {
        (
            Some(
                matches
                    .iter()
                    .map(|matched| matched.error_ppm.abs())
                    .sum::<f64>()
                    / matches.len() as f64,
            ),
            Some(
                matches
                    .iter()
                    .map(|matched| matched.error_da.abs())
                    .sum::<f64>()
                    / matches.len() as f64,
            ),
        )
    };

    AnnotationQualityMetrics {
        snr_like,
        log2_snr_like,
        cosine,
        frag_error_mae_ppm,
        frag_error_mae_da,
    }
}

fn should_integrate_fragment_isotopes(fragment: &FragmentIon, tolerance: MassTolerance) -> bool {
    let isotope_spacing = NEUTRON_MASS_DIFF / fragment.charge as f64;
    isotope_spacing > 2.0 * tolerance.window_da(fragment.theoretical_mz)
}

pub(crate) fn fragment_charge_states(charge_context: Option<i32>) -> Vec<u8> {
    match charge_context {
        Some(charge) if charge > 1 => (1..charge)
            .filter(|value| *value <= u8::MAX as i32)
            .map(|value| value as u8)
            .collect(),
        Some(_) => vec![1],
        None => vec![1, 2],
    }
}

fn parse_peptide_input(input: &str) -> anyhow::Result<ParsedPeptideInput> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        anyhow::bail!("peptide input cannot be empty");
    }

    let (sequence_input, charge_hint) = if let Some((head, tail)) = trimmed.rsplit_once('/') {
        if !tail.is_empty() && tail.chars().all(|ch| ch.is_ascii_digit()) {
            let charge_hint = tail
                .parse::<i32>()
                .with_context(|| format!("invalid charge suffix in peptide `{trimmed}`"))?;
            (head, Some(charge_hint))
        } else {
            (trimmed, None)
        }
    } else {
        (trimmed, None)
    };

    let mut sequence = String::new();
    let mut inline_modifications = Vec::<ExplicitModification>::new();
    let bytes = sequence_input.as_bytes();
    let mut idx = 0usize;

    if bytes.first().copied() == Some(b'[') {
        let (delta, next_idx) = parse_inline_delta(sequence_input, 0)?;
        inline_modifications.push(ExplicitModification {
            site: ModificationSite::NTerm,
            delta,
        });
        idx = next_idx;
    }

    while idx < bytes.len() {
        let aa = bytes[idx] as char;
        if !aa.is_ascii_alphabetic() {
            anyhow::bail!(
                "peptide input must contain amino-acid letters with optional inline mass shifts like M[+15.9949], S[+79.9663], [+42.0106]PEPTIDE, and optional `/charge` suffix"
            );
        }
        sequence.push(aa.to_ascii_uppercase());
        idx += 1;

        if idx < bytes.len() && bytes[idx] == b'[' {
            let (delta, next_idx) = parse_inline_delta(sequence_input, idx)?;
            inline_modifications.push(ExplicitModification {
                site: ModificationSite::Residue(sequence.len()),
                delta,
            });
            idx = next_idx;
        }
    }

    if sequence.is_empty() {
        anyhow::bail!("peptide sequence cannot be empty");
    }

    Ok(ParsedPeptideInput {
        sequence,
        inline_modifications,
        charge_hint,
    })
}

fn parse_mod_spec(input: &str) -> anyhow::Result<ExplicitModification> {
    let (position_raw, delta_raw) = input
        .split_once(':')
        .ok_or_else(|| anyhow::anyhow!("invalid --mod `{input}` (expected <position>:<delta>)"))?;
    let position = position_raw
        .parse::<usize>()
        .with_context(|| format!("invalid modification position in `{input}`"))?;
    let delta = delta_raw
        .parse::<f64>()
        .with_context(|| format!("invalid modification delta in `{input}`"))?;
    Ok(ExplicitModification {
        site: ModificationSite::Residue(position),
        delta,
    })
}

fn build_peptide(
    sequence: &str,
    charge_hint: Option<i32>,
    modifications: &[ExplicitModification],
) -> anyhow::Result<Peptide> {
    let residue_chars = sequence.chars().collect::<Vec<_>>();
    if residue_chars.is_empty() {
        anyhow::bail!("peptide sequence cannot be empty");
    }

    let mut delta_by_position = vec![0.0_f64; residue_chars.len()];
    let mut n_term_delta = 0.0_f64;
    let mut c_term_delta = 0.0_f64;
    for modification in modifications {
        match modification.site {
            ModificationSite::NTerm => {
                n_term_delta += modification.delta;
            }
            ModificationSite::CTerm => {
                c_term_delta += modification.delta;
            }
            ModificationSite::Residue(position) => {
                if position == 0 || position > residue_chars.len() {
                    anyhow::bail!(
                        "modification position {} out of range for peptide length {}",
                        position,
                        residue_chars.len()
                    );
                }
                delta_by_position[position - 1] += modification.delta;
            }
        }
    }

    let mut residue_masses = Vec::with_capacity(residue_chars.len());
    for (idx, aa) in residue_chars.iter().copied().enumerate() {
        let base = residue_mass(aa)
            .ok_or_else(|| anyhow::anyhow!("unsupported amino-acid code `{aa}`"))?;
        residue_masses.push(base + delta_by_position[idx]);
    }

    Ok(Peptide {
        sequence: sequence.to_string(),
        residue_chars,
        residue_deltas: delta_by_position,
        residue_masses,
        n_term_delta,
        c_term_delta,
        charge_hint,
    })
}

fn parse_inline_delta(input: &str, start_idx: usize) -> anyhow::Result<(f64, usize)> {
    let rest = &input[start_idx + 1..];
    let end_rel = rest
        .find(']')
        .ok_or_else(|| anyhow::anyhow!("missing closing `]` in peptide `{input}`"))?;
    let raw = rest[..end_rel].trim();
    if raw.is_empty() {
        anyhow::bail!("empty mass shift in peptide `{input}`");
    }
    let normalized = if raw.starts_with('+') || raw.starts_with('-') {
        raw.to_string()
    } else {
        format!("+{raw}")
    };
    let delta = normalized
        .parse::<f64>()
        .with_context(|| format!("invalid mass shift `{raw}` in peptide `{input}`"))?;
    Ok((delta, start_idx + end_rel + 2))
}

pub(crate) fn residue_mass(aa: char) -> Option<f64> {
    match aa {
        'A' => Some(71.037_113_805),
        'R' => Some(156.101_111_05),
        'N' => Some(114.042_927_47),
        'D' => Some(115.026_943_065),
        'C' => Some(103.009_184_505),
        'E' => Some(129.042_593_135),
        'Q' => Some(128.058_577_54),
        'G' => Some(57.021_463_735),
        'H' => Some(137.058_911_875),
        'I' => Some(113.084_064_015),
        'L' => Some(113.084_064_015),
        'K' => Some(128.094_963_05),
        'M' => Some(131.040_484_645),
        'F' => Some(147.068_413_945),
        'P' => Some(97.052_763_875),
        'S' => Some(87.032_028_435),
        'T' => Some(101.047_678_505),
        'W' => Some(186.079_312_98),
        'Y' => Some(163.063_328_575),
        'U' => Some(150.953_633_405),
        'V' => Some(99.068_413_945),
        _ => None,
    }
}

pub(crate) fn generate_fragments(
    peptide: &Peptide,
    charges: &[u8],
    neutral_losses: &[NeutralLossKind],
) -> Vec<FragmentIon> {
    if peptide.len() < 2 {
        return Vec::new();
    }

    let mut prefix_mass = Vec::with_capacity(peptide.len());
    let mut running = 0.0_f64;
    for mass in &peptide.residue_masses {
        running += *mass;
        prefix_mass.push(running);
    }
    let total_residue_mass = peptide.total_residue_mass();

    let mut fragments = Vec::with_capacity((peptide.len() - 1) * charges.len() * 2);
    for cleavage_index in 1..peptide.len() {
        let prefix_residue_mass = prefix_mass[cleavage_index - 1];
        let b_residue_mass = prefix_residue_mass + peptide.n_term_delta;
        let b_residues = &peptide.residue_chars[..cleavage_index];
        let b_residue_deltas = &peptide.residue_deltas[..cleavage_index];
        let y_ordinal = peptide.len() - cleavage_index;
        let y_residue_mass = total_residue_mass - prefix_residue_mass + peptide.c_term_delta;
        let y_residues = &peptide.residue_chars[cleavage_index..];
        let y_residue_deltas = &peptide.residue_deltas[cleavage_index..];
        push_fragment_variants(
            &mut fragments,
            FragmentSeries::B,
            cleavage_index,
            cleavage_index,
            b_residue_mass,
            b_residues,
            b_residue_deltas,
            charges,
            neutral_losses,
        );
        push_fragment_variants(
            &mut fragments,
            FragmentSeries::Y,
            y_ordinal,
            cleavage_index,
            y_residue_mass + WATER_MASS,
            y_residues,
            y_residue_deltas,
            charges,
            neutral_losses,
        );
    }
    fragments.sort_by(|left, right| {
        left.theoretical_mz
            .partial_cmp(&right.theoretical_mz)
            .unwrap_or(Ordering::Equal)
            .then_with(|| left.series.cmp(&right.series))
            .then_with(|| left.ordinal.cmp(&right.ordinal))
            .then_with(|| left.neutral_loss.cmp(&right.neutral_loss))
            .then_with(|| left.charge.cmp(&right.charge))
    });
    fragments
}

fn push_fragment_variants(
    fragments: &mut Vec<FragmentIon>,
    series: FragmentSeries,
    ordinal: usize,
    cleavage_index: usize,
    neutral_mass: f64,
    residues: &[char],
    residue_deltas: &[f64],
    charges: &[u8],
    neutral_losses: &[NeutralLossKind],
) {
    for charge in charges.iter().copied() {
        let z = charge as f64;
        fragments.push(FragmentIon {
            series,
            ordinal,
            cleavage_index,
            charge,
            neutral_loss: None,
            theoretical_mz: (neutral_mass + z * PROTON_MASS) / z,
        });

        for loss in neutral_losses.iter().copied() {
            if !fragment_supports_neutral_loss(residues, residue_deltas, loss) {
                continue;
            }
            let shifted_neutral_mass = neutral_mass - loss.mass();
            if shifted_neutral_mass <= 0.0 {
                continue;
            }
            fragments.push(FragmentIon {
                series,
                ordinal,
                cleavage_index,
                charge,
                neutral_loss: Some(loss),
                theoretical_mz: (shifted_neutral_mass + z * PROTON_MASS) / z,
            });
        }
    }
}

fn fragment_supports_neutral_loss(
    residues: &[char],
    residue_deltas: &[f64],
    loss: NeutralLossKind,
) -> bool {
    residues
        .iter()
        .copied()
        .zip(residue_deltas.iter().copied())
        .any(|(aa, delta)| match loss {
            NeutralLossKind::Water => matches!(aa, 'S' | 'T' | 'D' | 'E'),
            NeutralLossKind::Ammonia => matches!(aa, 'K' | 'N' | 'Q' | 'R'),
            NeutralLossKind::PhosphoricAcid => {
                matches!(aa, 'S' | 'T') && is_phospho_like_delta(delta)
            }
        })
}

fn is_phospho_like_delta(delta: f64) -> bool {
    (delta - PHOSPHO_DELTA_MASS).abs() <= PHOSPHO_DELTA_EPSILON
}

fn match_fragments(
    fragments: &[FragmentIon],
    mz: &[f64],
    intensity: &[f32],
    tolerance: MassTolerance,
) -> Vec<FragmentMatch> {
    let observed = sorted_observed_peaks(mz, intensity);

    let mut out = Vec::new();
    for fragment in fragments {
        if let Some(peak) = best_peak_for_target(&observed, fragment.theoretical_mz, tolerance) {
            let error_da = peak.mz - fragment.theoretical_mz;
            out.push(FragmentMatch {
                fragment: fragment.clone(),
                peak_index: peak.index,
                observed_mz: peak.mz,
                observed_intensity: peak.intensity,
                error_da,
                error_ppm: tolerance.error_ppm(fragment.theoretical_mz, peak.mz),
            });
        }
    }

    out.sort_by(|left, right| {
        left.observed_mz
            .partial_cmp(&right.observed_mz)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                left.error_da
                    .abs()
                    .partial_cmp(&right.error_da.abs())
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| {
                right
                    .observed_intensity
                    .partial_cmp(&left.observed_intensity)
                    .unwrap_or(Ordering::Equal)
            })
    });
    out
}

fn sorted_observed_peaks(mz: &[f64], intensity: &[f32]) -> Vec<ObservedPeak> {
    let mut observed = mz
        .iter()
        .copied()
        .zip(intensity.iter().copied())
        .enumerate()
        .filter_map(|(index, (mz, intensity))| {
            if mz.is_finite() && intensity.is_finite() {
                Some(ObservedPeak {
                    index,
                    mz,
                    intensity,
                })
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    observed.sort_by(|left, right| left.mz.partial_cmp(&right.mz).unwrap_or(Ordering::Equal));
    observed
}

fn best_peak_for_target(
    observed: &[ObservedPeak],
    theoretical_mz: f64,
    tolerance: MassTolerance,
) -> Option<&ObservedPeak> {
    let tolerance_da = tolerance.window_da(theoretical_mz);
    let insertion = observed.partition_point(|peak| peak.mz < theoretical_mz);
    let mut best: Option<&ObservedPeak> = None;

    let mut left = insertion;
    while left > 0 {
        left -= 1;
        let candidate = &observed[left];
        let error_da = (candidate.mz - theoretical_mz).abs();
        if error_da > tolerance_da && candidate.mz < theoretical_mz {
            break;
        }
        if !tolerance.contains(theoretical_mz, candidate.mz) {
            continue;
        }
        if is_better_peak(candidate, best, theoretical_mz) {
            best = Some(candidate);
        }
    }

    let mut right = insertion;
    while right < observed.len() {
        let candidate = &observed[right];
        let error_da = (candidate.mz - theoretical_mz).abs();
        if error_da > tolerance_da && candidate.mz > theoretical_mz {
            break;
        }
        if tolerance.contains(theoretical_mz, candidate.mz)
            && is_better_peak(candidate, best, theoretical_mz)
        {
            best = Some(candidate);
        }
        right += 1;
    }

    best
}

fn is_better_peak<'a>(
    candidate: &'a ObservedPeak,
    best: Option<&'a ObservedPeak>,
    theoretical_mz: f64,
) -> bool {
    let Some(best) = best else {
        return true;
    };

    let candidate_error = (candidate.mz - theoretical_mz).abs();
    let best_error = (best.mz - theoretical_mz).abs();
    if candidate_error < best_error {
        true
    } else if (candidate_error - best_error).abs() <= 1e-12 {
        if candidate.intensity > best.intensity {
            true
        } else if (candidate.intensity - best.intensity).abs() <= f32::EPSILON {
            candidate.index < best.index
        } else {
            false
        }
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_mod_spec() {
        let parsed = parse_mod_spec("7:+57.021464").expect("mod parses");
        assert_eq!(parsed.site, ModificationSite::Residue(7));
        assert!((parsed.delta - 57.021_464).abs() < 1e-9);
    }

    #[test]
    fn parses_plain_peptide_with_charge_suffix() {
        let context = prepare_annotation(
            "PEPTIDEK/2",
            &["3:+15.994915".to_string()],
            &[],
            None,
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.peptide.sequence(), "PEPTIDEK");
        assert_eq!(context.charge_context, Some(2));
        assert_eq!(context.modifications.len(), 1);
    }

    #[test]
    fn parses_inline_residue_shift_with_charge_suffix() {
        let context = prepare_annotation(
            "PEPM[+15.9949]IDE/2",
            &[],
            &[],
            None,
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.peptide.sequence(), "PEPMIDE");
        assert_eq!(context.charge_context, Some(2));
        assert_eq!(context.modifications.len(), 1);
        assert_eq!(context.modifications[0].site, ModificationSite::Residue(4));
        assert!((context.modifications[0].delta - 15.9949).abs() < 1e-9);
    }

    #[test]
    fn parses_inline_positive_shift_without_sign() {
        let context = prepare_annotation(
            "PEPS[79.9663]TIDE",
            &[],
            &[],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.modifications.len(), 1);
        assert_eq!(context.modifications[0].site, ModificationSite::Residue(4));
        assert!((context.modifications[0].delta - 79.9663).abs() < 1e-9);
    }

    #[test]
    fn parses_inline_negative_shift() {
        let context = prepare_annotation(
            "PEPE[-18.0106]K",
            &[],
            &[],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.modifications.len(), 1);
        assert_eq!(context.modifications[0].site, ModificationSite::Residue(4));
        assert!((context.modifications[0].delta + 18.0106).abs() < 1e-9);
    }

    #[test]
    fn parses_n_term_inline_shift() {
        let context = prepare_annotation(
            "[+42.0106]PEPTIDE",
            &[],
            &[],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.modifications.len(), 1);
        assert_eq!(context.modifications[0].site, ModificationSite::NTerm);
        assert!((context.modifications[0].delta - 42.0106).abs() < 1e-9);
        assert_eq!(
            context.modifications_label().as_deref(),
            Some("n-term:+42.010600")
        );
    }

    #[test]
    fn explicit_c_term_shift_changes_precursor_and_y_ions() {
        let plain = prepare_annotation("PEPUK/2", &[], &[], Some(2), MassTolerance::Ppm(20.0))
            .expect("plain annotation");
        let shifted = prepare_annotation_with_modifications(
            "PEPUK",
            vec![ExplicitModification {
                site: ModificationSite::CTerm,
                delta: 17.003,
            }],
            &[],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("shifted annotation");

        assert_eq!(shifted.peptide.sequence(), "PEPUK");
        assert_eq!(shifted.charge_context, Some(2));
        assert_eq!(shifted.modified_sequence(), "PEPUK[+17.003000]");

        let plain_precursor = plain.peptide.precursor_mz(2).expect("plain precursor");
        let shifted_precursor = shifted.peptide.precursor_mz(2).expect("shifted precursor");
        assert!((shifted_precursor - plain_precursor - 17.003 / 2.0).abs() < 1e-9);

        let plain_y1 = generate_fragments(&plain.peptide, &[1], &[])
            .into_iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("plain y1");
        let shifted_y1 = generate_fragments(&shifted.peptide, &[1], &[])
            .into_iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("shifted y1");
        assert!((shifted_y1.theoretical_mz - plain_y1.theoretical_mz - 17.003).abs() < 1e-9);
    }

    #[test]
    fn merges_inline_and_cli_modifications() {
        let context = prepare_annotation(
            "[+42.0106]PEPM[+15.9949]IDE",
            &["2:+79.9663".to_string()],
            &[],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        assert_eq!(context.modifications.len(), 3);
        assert_eq!(context.modifications[0].site, ModificationSite::NTerm);
        assert_eq!(context.modifications[1].site, ModificationSite::Residue(4));
        assert_eq!(context.modifications[2].site, ModificationSite::Residue(2));
    }

    #[test]
    fn inline_n_term_shift_changes_b_ions_and_precursor_mass() {
        let plain = prepare_annotation("PEP", &[], &[], Some(2), MassTolerance::Ppm(20.0))
            .expect("plain context");
        let shifted =
            prepare_annotation("[+42.0106]PEP", &[], &[], Some(2), MassTolerance::Ppm(20.0))
                .expect("shifted context");
        let plain_fragments = generate_fragments(&plain.peptide, &[1], &[]);
        let shifted_fragments = generate_fragments(&shifted.peptide, &[1], &[]);
        let plain_b1 = plain_fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::B && fragment.ordinal == 1)
            .expect("plain b1")
            .theoretical_mz;
        let shifted_b1 = shifted_fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::B && fragment.ordinal == 1)
            .expect("shifted b1")
            .theoretical_mz;
        let plain_y1 = plain_fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("plain y1")
            .theoretical_mz;
        let shifted_y1 = shifted_fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("shifted y1")
            .theoretical_mz;
        assert!((shifted_b1 - plain_b1 - 42.0106).abs() < 1e-6);
        assert!((shifted_y1 - plain_y1).abs() < 1e-6);

        let plain_precursor = plain.peptide.precursor_mz(2).expect("plain precursor");
        let shifted_precursor = shifted.peptide.precursor_mz(2).expect("shifted precursor");
        assert!((shifted_precursor - plain_precursor - 21.0053).abs() < 1e-4);
    }

    #[test]
    fn precursor_check_prefers_allowed_isotope_error_match() {
        let context = prepare_annotation(
            "[+304.2071]T[+79.9663]S[+79.9663]SSSPSR/3",
            &[],
            &[],
            Some(3),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let monoisotopic = context.peptide.precursor_mz(3).expect("precursor mz");
        let observed = monoisotopic - 1.003_354_835_07 / 3.0 + 0.0005;
        let report = annotate_peaks(&context, Some(observed), &[], &[]);
        let check = report.precursor_check.expect("precursor check");
        assert_eq!(check.isotope_error, 1);
        assert!(check.within_tolerance);
        assert!((check.monoisotopic_theoretical_mz - monoisotopic).abs() < 1e-9);
    }

    #[test]
    fn rejects_malformed_inline_mass_shift() {
        let err = prepare_annotation("PEPM[+15.9949", &[], &[], Some(2), MassTolerance::Ppm(20.0))
            .expect_err("missing bracket should fail");
        assert!(err.to_string().contains("missing closing `]`"));
    }

    #[test]
    fn generates_expected_fragment_count() {
        let context = prepare_annotation("PEPTIDE", &[], &[], Some(3), MassTolerance::Ppm(20.0))
            .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert_eq!(report.fragments.len(), (context.peptide.len() - 1) * 4);
    }

    #[test]
    fn two_plus_precursors_only_generate_singly_charged_fragments() {
        let context = prepare_annotation("PEPTIDE", &[], &[], Some(2), MassTolerance::Ppm(20.0))
            .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert!(report.fragments.iter().all(|fragment| fragment.charge == 1));
        assert_eq!(report.fragments.len(), (context.peptide.len() - 1) * 2);
    }

    #[test]
    fn matches_fragments_in_da_space() {
        let context = prepare_annotation("PEP", &[], &[], Some(1), MassTolerance::Da(0.02))
            .expect("annotation context");
        let fragments = generate_fragments(&context.peptide, &[1], &[]);
        let target = fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("y1 fragment")
            .theoretical_mz;
        let mz = [50.0, target + 0.01, 400.0];
        let intensity = [10.0, 250.0, 5.0];
        let matches = match_fragments(&fragments, &mz, &intensity, MassTolerance::Da(0.02));
        assert!(matches.iter().any(|matched| {
            matched.fragment.series == FragmentSeries::Y
                && matched.fragment.ordinal == 1
                && (matched.error_da - 0.01).abs() < 1e-6
        }));
    }

    #[test]
    fn annotation_quality_metrics_use_unique_matched_peaks() {
        let context = prepare_annotation("PEP", &[], &[], Some(1), MassTolerance::Da(0.02))
            .expect("annotation context");
        let fragments = generate_fragments(&context.peptide, &[1], &[]);
        let target = fragments
            .iter()
            .find(|fragment| fragment.series == FragmentSeries::Y && fragment.ordinal == 1)
            .expect("y1 fragment")
            .theoretical_mz;
        let mz = [target + 0.01, target + NEUTRON_MASS_DIFF + 0.01, 1000.0];
        let intensity = [100.0, 25.0, 50.0];
        let report = annotate_peaks(&context, None, &mz, &intensity);

        assert_eq!(report.matched_peak_count(), 1);
        assert!((report.quality.snr_like - 2.5).abs() < 1e-9);
        assert!((report.quality.log2_snr_like - (126.0_f64 / 51.0).log2()).abs() < 1e-9);
        let matched_sum_squares = 100.0_f64 * 100.0 + 25.0 * 25.0;
        let observed_sum_squares = matched_sum_squares + 50.0 * 50.0;
        assert!(
            (report.quality.cosine
                - (matched_sum_squares
                    / (matched_sum_squares.sqrt() * observed_sum_squares.sqrt())))
            .abs()
                < 1e-9
        );
        assert!(
            report
                .quality
                .frag_error_mae_ppm
                .expect("fragment error mae")
                > 0.0
        );
    }

    #[test]
    fn annotation_quality_skips_isotopes_when_spacing_is_not_two_tolerances() {
        let context = prepare_annotation("PEP/3", &[], &[], Some(3), MassTolerance::Da(0.5))
            .expect("annotation context");
        let fragments = generate_fragments(&context.peptide, &[1, 2], &[]);
        let target = fragments
            .iter()
            .find(|fragment| {
                fragment.series == FragmentSeries::Y
                    && fragment.ordinal == 1
                    && fragment.charge == 2
            })
            .expect("y1++ fragment")
            .theoretical_mz;
        let mz = [
            target + 0.01,
            target + NEUTRON_MASS_DIFF / 2.0 + 0.01,
            1000.0,
        ];
        let intensity = [100.0, 25.0, 50.0];
        let report = annotate_peaks(&context, None, &mz, &intensity);

        assert_eq!(report.matched_peak_count(), 1);
        assert!((report.quality.snr_like - (100.0_f64 / 75.0)).abs() < 1e-9);
        assert!((report.quality.log2_snr_like - (101.0_f64 / 76.0).log2()).abs() < 1e-9);
    }

    #[test]
    fn fragment_isotope_gate_requires_spacing_above_two_tolerances() {
        let base = FragmentIon {
            series: FragmentSeries::Y,
            ordinal: 1,
            cleavage_index: 1,
            charge: 2,
            neutral_loss: None,
            theoretical_mz: 500.0,
        };
        assert!(should_integrate_fragment_isotopes(
            &base,
            MassTolerance::Da(0.25)
        ));
        assert!(!should_integrate_fragment_isotopes(
            &base,
            MassTolerance::Da(0.5)
        ));

        let mut charge_three = base.clone();
        charge_three.charge = 3;
        assert!(!should_integrate_fragment_isotopes(
            &charge_three,
            MassTolerance::Da(0.25)
        ));
    }

    #[test]
    fn generates_residue_aware_neutral_losses_when_enabled() {
        let context = prepare_annotation(
            "DSATK",
            &[],
            &[NeutralLossKind::Water, NeutralLossKind::Ammonia],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert!(report
            .fragments
            .iter()
            .any(|fragment| fragment.label() == "b1-H2O"));
        assert!(report
            .fragments
            .iter()
            .any(|fragment| fragment.label() == "y1-NH3"));
    }

    #[test]
    fn generates_phospho_neutral_losses_only_for_phospho_fragments() {
        let context = prepare_annotation(
            "AS[+79.9663]TK",
            &[],
            &[NeutralLossKind::PhosphoricAcid],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert!(report
            .fragments
            .iter()
            .any(|fragment| fragment.label() == "b2-H3PO4"));
        assert!(report
            .fragments
            .iter()
            .any(|fragment| fragment.label() == "y3-H3PO4"));
        assert!(!report
            .fragments
            .iter()
            .any(|fragment| fragment.label() == "b1-H3PO4"));
    }

    #[test]
    fn skips_phospho_neutral_losses_for_unmodified_serine_threonine() {
        let context = prepare_annotation(
            "ASTK",
            &[],
            &[NeutralLossKind::PhosphoricAcid],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert!(!report
            .fragments
            .iter()
            .any(|fragment| fragment.neutral_loss == Some(NeutralLossKind::PhosphoricAcid)));
    }

    #[test]
    fn skips_unsupported_neutral_losses() {
        let context = prepare_annotation(
            "PFGI",
            &[],
            &[NeutralLossKind::Water, NeutralLossKind::Ammonia],
            Some(2),
            MassTolerance::Ppm(20.0),
        )
        .expect("annotation context");
        let report = annotate_peaks(&context, None, &[], &[]);
        assert!(!report
            .fragments
            .iter()
            .any(|fragment| fragment.neutral_loss == Some(NeutralLossKind::Water)));
        assert!(!report
            .fragments
            .iter()
            .any(|fragment| fragment.neutral_loss == Some(NeutralLossKind::Ammonia)));
    }
}
