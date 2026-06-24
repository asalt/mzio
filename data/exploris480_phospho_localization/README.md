# Orbitrap Exploris 480 phospho-localization fixture

Stable slim fixture data for pepXML-backed phospho localization comparisons.

Contents:

- `mzml/49108_1_EXP_802_PDX_1mg_phos_oneforth_F13_scan4877.mzML`
  - One-scan mzML made from a full phospho-enriched DDA mzML.
  - Scan `4877`, charge `2`, 98 peaks.
  - About 552 KB.
- `pepxml/49108_F13_scan4877_phospho_top_hits.pep.xml`
  - One `spectrum_query` from the matching pepXML with five `search_hit`
    entries.
- `manifest.tsv`
  - Minimal source and scan metadata.

Why this scan:

- Same bare peptide in ranks 1-5: `TGSESSQTGTSTTSSR`.
- Alternate phospho localization across S/T sites.
- Useful for residue mass deltas, neutral-loss rendering, rank-specific output
  names, JSON sidecars, and localization comparison plots.

Top hits:

- rank 1: `TGS[+79.966324]ESSQTGTSTTSSR`
- rank 2: `TGSES[+79.966324]SQTGTSTTSSR`
- rank 3: `T[+79.966329]GSESSQTGTSTTSSR`
- rank 4: `TGSESS[+79.966324]QTGTSTTSSR`
- rank 5: `TGSESSQT[+79.966329]GTSTTSSR`

Try it:

```bash
cargo run -- plot \
  --mzml data/exploris480_phospho_localization/mzml/49108_1_EXP_802_PDX_1mg_phos_oneforth_F13_scan4877.mzML \
  --scan 4877 \
  --pepxml data/exploris480_phospho_localization/pepxml/49108_F13_scan4877_phospho_top_hits.pep.xml \
  --top-n 5 \
  --neutral-losses \
  --svg-prefix exploris480_phospho_localization
```

To regenerate the slim mzML, place the full source mzML here:

```text
data/_sources/exploris480_phospho_localization/49108_1_EXP_802_PDX_1mg_phos_oneforth_F13.mzML
```

Then run:

```bash
bash data/regenerate_slim_fixtures.sh
```
