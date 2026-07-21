use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::annotate::{FragmentSeries, NeutralLossKind};

const FONT_FAMILY: &str = "Menlo, Consolas, Liberation Mono, monospace";
const COLOR_TEXT: &str = "#122033";
const COLOR_SUBTLE: &str = "#5b6775";
const COLOR_AXIS: &str = "#334155";
const COLOR_GRID: &str = "#e5ebf2";
const COLOR_ROW_GRID: &str = "#eef2f6";
const COLOR_BORDER: &str = "#d8e0ea";
const COLOR_B: &str = "#1d4ed8";
const COLOR_Y: &str = "#b45309";
const COLOR_MISSING: &str = "#aeb8c4";

const PAD: f64 = 18.0;
const COLUMN_GAP: f64 = 14.0;
const MIN_VALUE_WIDTH: f64 = 154.0;
const MIN_CENTER_WIDTH: f64 = 176.0;
const HEADER_HEIGHT: f64 = 126.0;
const FOOTER_HEIGHT: f64 = 72.0;
const BASE_FONT_SIZE: f64 = 15.0;
const LOSS_FONT_SIZE: f64 = 13.0;
const HEADER_FONT_SIZE: f64 = 13.5;
const POSITION_FONT_SIZE: f64 = 15.0;
const RESIDUE_FONT_SIZE: f64 = 16.0;
const CELL_LINE_HEIGHT: f64 = 18.0;

#[derive(Clone, Debug)]
pub(crate) struct SvgIonTableEntry {
    pub(crate) series: FragmentSeries,
    pub(crate) ordinal: usize,
    pub(crate) charge: u8,
    pub(crate) neutral_loss: Option<NeutralLossKind>,
    pub(crate) mz: f64,
    pub(crate) detected: bool,
    pub(crate) title: String,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct SvgIonTableCell {
    pub(crate) entries: Vec<SvgIonTableEntry>,
}

#[derive(Clone, Debug)]
pub(crate) struct SvgIonTableRow {
    pub(crate) n_position: usize,
    pub(crate) c_position: usize,
    pub(crate) residue_label: String,
    pub(crate) b: BTreeMap<u8, SvgIonTableCell>,
    pub(crate) y: BTreeMap<u8, SvgIonTableCell>,
}

#[derive(Clone, Debug)]
pub(crate) struct SvgIonTable {
    pub(crate) title: String,
    pub(crate) evidence_legend: String,
    pub(crate) loss_legend: String,
    pub(crate) sequence: String,
    pub(crate) footer_note: Option<String>,
    pub(crate) charges: Vec<u8>,
    pub(crate) rows: Vec<SvgIonTableRow>,
}

#[derive(Clone, Debug)]
pub(crate) struct SvgIonTableLayout {
    pub(crate) width: f64,
    pub(crate) height: f64,
    row_heights: Vec<f64>,
    center_width: f64,
    value_width: f64,
    footer_height: f64,
}

impl SvgIonTable {
    pub(crate) fn layout(&self, minimum_width: f64) -> SvgIonTableLayout {
        let charges = self.charges.len().max(1);
        let value_columns = charges * 2;
        let max_cell_width = self
            .rows
            .iter()
            .flat_map(|row| row.b.values().chain(row.y.values()))
            .flat_map(|cell| cell.entries.iter())
            .map(estimated_entry_width)
            .fold(MIN_VALUE_WIDTH, f64::max);
        let center_width = self
            .rows
            .iter()
            .map(|row| estimate_text_width(&row.residue_label, RESIDUE_FONT_SIZE) + 82.0)
            .fold(MIN_CENTER_WIDTH, f64::max);
        let gaps = COLUMN_GAP * value_columns as f64;
        let natural_width = PAD * 2.0 + center_width + max_cell_width * value_columns as f64 + gaps;
        let width = natural_width.max(minimum_width);
        let value_width =
            ((width - PAD * 2.0 - center_width - gaps) / value_columns as f64).max(MIN_VALUE_WIDTH);
        let row_heights = self
            .rows
            .iter()
            .map(|row| {
                let lines = row
                    .b
                    .values()
                    .chain(row.y.values())
                    .map(|cell| cell.entries.len())
                    .max()
                    .unwrap_or(1)
                    .max(1);
                (lines as f64 * CELL_LINE_HEIGHT + 14.0).max(34.0)
            })
            .collect::<Vec<_>>();
        let footer_height = if self.footer_note.is_some() {
            FOOTER_HEIGHT
        } else {
            48.0
        };
        let height = HEADER_HEIGHT + row_heights.iter().sum::<f64>() + footer_height;
        SvgIonTableLayout {
            width,
            height,
            row_heights,
            center_width,
            value_width,
            footer_height,
        }
    }

    pub(crate) fn render(&self, svg: &mut String, left: f64, top: f64, layout: &SvgIonTableLayout) {
        let _ = writeln!(
            svg,
            "<rect x=\"{left:.2}\" y=\"{top:.2}\" width=\"{width:.2}\" height=\"{height:.2}\" rx=\"8\" fill=\"#fbfdff\" stroke=\"{COLOR_BORDER}\" stroke-width=\"1\"/>",
            width = layout.width,
            height = layout.height,
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"18\" fill=\"{COLOR_TEXT}\">{}</text>",
            left + PAD,
            top + 26.0,
            escape_xml(&self.title),
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"13\" fill=\"{COLOR_SUBTLE}\">{}</text>",
            left + PAD,
            top + 49.0,
            escape_xml(&self.evidence_legend),
        );
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"13\" fill=\"{COLOR_SUBTLE}\">{}</text>",
            left + PAD,
            top + 70.0,
            escape_xml(&self.loss_legend),
        );

        let charges = normalized_charges(&self.charges);
        let table_left = left + PAD;
        let mut cursor = table_left;
        let mut b_positions = BTreeMap::<u8, f64>::new();
        for charge in charges.iter().rev() {
            let right = cursor + layout.value_width;
            b_positions.insert(*charge, right);
            cursor = right + COLUMN_GAP;
        }
        let center_left = cursor;
        let center_right = center_left + layout.center_width;
        cursor = center_right + COLUMN_GAP;
        let mut y_positions = BTreeMap::<u8, f64>::new();
        for (idx, charge) in charges.iter().enumerate() {
            y_positions.insert(*charge, cursor);
            cursor += layout.value_width;
            if idx + 1 < charges.len() {
                cursor += COLUMN_GAP;
            }
        }

        let header_y = top + 103.0;
        for charge in charges.iter().rev() {
            if let Some(x) = b_positions.get(charge) {
                write_charge_header(svg, *x, header_y, "end", FragmentSeries::B, *charge);
            }
        }
        let center = (center_left + center_right) / 2.0;
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{header_y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{HEADER_FONT_SIZE:.1}\" fill=\"{COLOR_AXIS}\" text-anchor=\"middle\">N  sequence  C</text>",
            center,
        );
        for charge in &charges {
            if let Some(x) = y_positions.get(charge) {
                write_charge_header(svg, *x, header_y, "start", FragmentSeries::Y, *charge);
            }
        }
        let _ = writeln!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{COLOR_GRID}\" stroke-width=\"1\"/>",
            table_left,
            header_y + 8.0,
            left + layout.width - PAD,
            header_y + 8.0,
        );

        let mut row_top = top + HEADER_HEIGHT;
        for (idx, row) in self.rows.iter().enumerate() {
            let row_height = layout.row_heights[idx];
            if idx > 0 {
                let _ = writeln!(
                    svg,
                    "<line x1=\"{:.2}\" y1=\"{row_top:.2}\" x2=\"{:.2}\" y2=\"{row_top:.2}\" stroke=\"{COLOR_ROW_GRID}\" stroke-width=\"1\"/>",
                    table_left,
                    left + layout.width - PAD,
                );
            }
            let center_y = row_top + row_height / 2.0 + 4.0;
            for charge in charges.iter().rev() {
                if let (Some(x), Some(cell)) = (b_positions.get(charge), row.b.get(charge)) {
                    write_cell(svg, *x, row_top, row_height, "end", cell);
                }
            }
            let _ = writeln!(
                svg,
                "<text x=\"{:.2}\" y=\"{center_y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{POSITION_FONT_SIZE:.1}\" fill=\"{COLOR_SUBTLE}\" text-anchor=\"end\">{}</text>",
                center_left + 28.0,
                row.n_position,
            );
            let _ = writeln!(
                svg,
                "<text x=\"{center:.2}\" y=\"{center_y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{RESIDUE_FONT_SIZE:.1}\" font-weight=\"700\" fill=\"{COLOR_TEXT}\" text-anchor=\"middle\">{}</text>",
                escape_xml(&row.residue_label),
            );
            let _ = writeln!(
                svg,
                "<text x=\"{:.2}\" y=\"{center_y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{POSITION_FONT_SIZE:.1}\" fill=\"{COLOR_SUBTLE}\" text-anchor=\"start\">{}</text>",
                center_right - 28.0,
                row.c_position,
            );
            for charge in &charges {
                if let (Some(x), Some(cell)) = (y_positions.get(charge), row.y.get(charge)) {
                    write_cell(svg, *x, row_top, row_height, "start", cell);
                }
            }
            row_top += row_height;
        }

        let footer_top = top + layout.height - layout.footer_height;
        let _ = writeln!(
            svg,
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"13\" fill=\"{COLOR_SUBTLE}\">sequence: {}</text>",
            left + PAD,
            footer_top + 28.0,
            escape_xml(&self.sequence),
        );
        if let Some(note) = &self.footer_note {
            let _ = writeln!(
                svg,
                "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"13\" fill=\"{COLOR_SUBTLE}\">{}</text>",
                left + PAD,
                footer_top + 52.0,
                escape_xml(note),
            );
        }
    }
}

pub(crate) fn fragment_label_markup(
    series: FragmentSeries,
    ordinal: usize,
    charge: u8,
    neutral_loss: Option<NeutralLossKind>,
    include_charge: bool,
) -> String {
    let mut out = String::new();
    out.push(match series {
        FragmentSeries::B => 'b',
        FragmentSeries::Y => 'y',
    });
    let _ = write!(
        out,
        "<tspan baseline-shift=\"sub\" font-size=\"75%\">{ordinal}</tspan>"
    );
    if let Some(loss) = neutral_loss {
        out.push_str(&neutral_loss_markup(loss));
    }
    if include_charge && charge > 1 {
        let _ = write!(
            out,
            "<tspan baseline-shift=\"super\" font-size=\"75%\">{charge}+</tspan>"
        );
    }
    out
}

pub(crate) fn series_color(series: FragmentSeries) -> &'static str {
    match series {
        FragmentSeries::B => COLOR_B,
        FragmentSeries::Y => COLOR_Y,
    }
}

fn normalized_charges(charges: &[u8]) -> Vec<u8> {
    let mut out = charges.to_vec();
    out.sort_unstable();
    out.dedup();
    if out.is_empty() {
        out.push(1);
    }
    out
}

fn write_charge_header(
    svg: &mut String,
    x: f64,
    y: f64,
    anchor: &str,
    series: FragmentSeries,
    charge: u8,
) {
    let prefix = match series {
        FragmentSeries::B => "b",
        FragmentSeries::Y => "y",
    };
    let charge_text = if charge <= 1 {
        "+".to_string()
    } else {
        format!("{charge}+")
    };
    let _ = writeln!(
        svg,
        "<text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{HEADER_FONT_SIZE:.1}\" fill=\"{COLOR_AXIS}\" text-anchor=\"{anchor}\">{prefix}<tspan baseline-shift=\"super\" font-size=\"75%\">{charge_text}</tspan></text>"
    );
}

fn write_cell(
    svg: &mut String,
    x: f64,
    row_top: f64,
    row_height: f64,
    anchor: &str,
    cell: &SvgIonTableCell,
) {
    if cell.entries.is_empty() {
        return;
    }
    let content_height = cell.entries.len() as f64 * CELL_LINE_HEIGHT;
    let first_y = row_top + (row_height - content_height) / 2.0 + BASE_FONT_SIZE;
    for (line_idx, entry) in cell.entries.iter().enumerate() {
        let y = first_y + line_idx as f64 * CELL_LINE_HEIGHT;
        let color = if entry.detected {
            series_color(entry.series)
        } else {
            COLOR_MISSING
        };
        let (font_size, weight, opacity, markup) = if let Some(loss) = entry.neutral_loss {
            (
                LOSS_FONT_SIZE,
                "400",
                "0.78",
                format!("{}&#160;&#160;{:.2}", neutral_loss_markup(loss), entry.mz),
            )
        } else {
            (
                BASE_FONT_SIZE,
                if entry.detected { "700" } else { "500" },
                "1",
                format!(
                    "{}&#160;&#160;{:.2}",
                    fragment_label_markup(entry.series, entry.ordinal, entry.charge, None, false,),
                    entry.mz
                ),
            )
        };
        let _ = writeln!(
            svg,
            "<text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"{FONT_FAMILY}\" font-size=\"{font_size:.1}\" font-weight=\"{weight}\" fill=\"{color}\" fill-opacity=\"{opacity}\" text-anchor=\"{anchor}\"><title>{}</title>{markup}</text>",
            escape_xml(&entry.title),
        );
    }
}

fn neutral_loss_markup(loss: NeutralLossKind) -> String {
    match loss {
        NeutralLossKind::Water => {
            "&#8722;H<tspan baseline-shift=\"sub\" font-size=\"75%\">2</tspan>O".to_string()
        }
        NeutralLossKind::Ammonia => {
            "&#8722;NH<tspan baseline-shift=\"sub\" font-size=\"75%\">3</tspan>".to_string()
        }
        NeutralLossKind::PhosphoricAcid => "&#8722;H<tspan baseline-shift=\"sub\" font-size=\"75%\">3</tspan>PO<tspan baseline-shift=\"sub\" font-size=\"75%\">4</tspan>".to_string(),
    }
}

fn estimated_entry_width(entry: &SvgIonTableEntry) -> f64 {
    let text = if let Some(loss) = entry.neutral_loss {
        format!("-{}  {:.2}", loss.label(), entry.mz)
    } else {
        let series = match entry.series {
            FragmentSeries::B => 'b',
            FragmentSeries::Y => 'y',
        };
        format!("{series}{}  {:.2}", entry.ordinal, entry.mz)
    };
    let font = if entry.neutral_loss.is_some() {
        LOSS_FONT_SIZE
    } else {
        BASE_FONT_SIZE
    };
    estimate_text_width(&text, font) + 12.0
}

fn estimate_text_width(text: &str, font_size: f64) -> f64 {
    text.chars().count() as f64 * font_size * 0.62
}

fn escape_xml(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragment_label_uses_subscript_and_higher_charge_superscript() {
        let label = fragment_label_markup(
            FragmentSeries::Y,
            10,
            2,
            Some(NeutralLossKind::Ammonia),
            true,
        );
        assert!(label.starts_with("y<tspan baseline-shift=\"sub\""));
        assert!(label.contains("&#8722;NH"));
        assert!(label.contains(">2+</tspan>"));

        let singly_charged = fragment_label_markup(FragmentSeries::B, 9, 1, None, true);
        assert!(!singly_charged.contains("baseline-shift=\"super\""));
    }

    #[test]
    fn layout_expands_for_more_charge_columns_and_loss_lines() {
        let entry = |charge, neutral_loss| SvgIonTableEntry {
            series: FragmentSeries::B,
            ordinal: 3,
            charge,
            neutral_loss,
            mz: 345.123,
            detected: true,
            title: "test".to_string(),
        };
        let mut b = BTreeMap::new();
        b.insert(
            1,
            SvgIonTableCell {
                entries: vec![
                    entry(1, None),
                    entry(1, Some(NeutralLossKind::Water)),
                    entry(1, Some(NeutralLossKind::Ammonia)),
                ],
            },
        );
        let table = SvgIonTable {
            title: "Ion table".to_string(),
            evidence_legend: "legend".to_string(),
            loss_legend: "losses".to_string(),
            sequence: "PEPTIDE".to_string(),
            footer_note: None,
            charges: vec![1, 2, 3],
            rows: vec![SvgIonTableRow {
                n_position: 1,
                c_position: 7,
                residue_label: "P".to_string(),
                b,
                y: BTreeMap::new(),
            }],
        };
        let layout = table.layout(100.0);
        assert!(layout.width > 1100.0);
        assert!(layout.row_heights[0] >= 68.0);
    }
}
