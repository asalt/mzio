# Orbitrap Exploris 480 HeLa digest fixture

Stable slim fixture data for pepXML-backed `mzio plot` checks.

Contents:

- `mzml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065.mzML`
  - One-scan mzML made from a full HeLa digest DDA mzML.
  - Scan `9065`, charge `2`, 142 peaks.
  - About 1.3 MB.
- `pepxml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065_top_hits.pep.xml`
  - One `spectrum_query` from the matching pepXML with five `search_hit`
    entries.
- `manifest.tsv`
  - Minimal source and scan metadata.

This fixture is useful for `--pepxml`, `--top-n`, rank-specific plot filenames,
JSON sidecars, ordinary residue mass deltas, and lower-confidence alternative
hits.

Try it:

```bash
cargo run -- plot \
  --mzml data/exploris480_hela_digest/mzml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065.mzML \
  --scan 9065 \
  --pepxml data/exploris480_hela_digest/pepxml/99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065_top_hits.pep.xml \
  --top-n 3 \
  --svg-prefix exploris480_hela_digest
```

To regenerate the slim mzML, place the full source mzML here:

```text
data/_sources/exploris480_hela_digest/99990_236_EXP_Hela_100ng_200nl_Jan8.mzML
```

Then run:

```bash
bash data/regenerate_slim_fixtures.sh
```
