# mzio fixture data

This directory contains stable, slim real-data fixtures for `mzio plot`
integration checks. They complement the tiny synthetic files in `assets/`.

The committed fixtures are intentionally small one-scan mzML files plus matching
single-query pepXML excerpts. They are meant to stay stable; add new examples
when they cover a distinct behavior rather than regenerating these routinely.

Committed fixture sets:

- `exploris480_hela_digest`
  - Orbitrap Exploris 480 HeLa digest DDA MS/MS scan.
  - Five ranked pepXML hits for one charge-2 scan.
  - Useful for `--pepxml`, `--top-n`, rank-specific output names, JSON sidecars,
    and lower-confidence alternative hits.
- `exploris480_phospho_localization`
  - Orbitrap Exploris 480 phospho-enriched DDA MS/MS scan.
  - Five ranked localizations for the same bare peptide.
  - Useful for residue mass deltas, neutral-loss rendering, and localization
    comparison plots.
- `timstof2_hela_dia_slice`
  - Slim Bruker timsTOF Ultra HeLa diaPASEF `.d` slice generated from a local
    full raw run.
  - Keeps 100 frames around the DIA-NN precursor `IILDLISESPIK/2`.
  - Useful for opt-in native Bruker `mzio dia-slice` testing, pseudo-MS2 output,
    trace/peak sidecars, and run-level TIC JSON.

Local regeneration inputs, if needed, should be staged under ignored repo-local
paths:

- `data/_sources/exploris480_hela_digest/99990_236_EXP_Hela_100ng_200nl_Jan8.mzML`
- `data/_sources/exploris480_phospho_localization/49108_1_EXP_802_PDX_1mg_phos_oneforth_F13.mzML`
- `data/timstof2_hela/DIA/99992_111_TOF_HeLa_DIA_20ng_Apr18_calib_610.d`

Run `bash data/regenerate_slim_fixtures.sh` after staging those source mzMLs.
The script only writes the committed `data/*/mzml/*.mzML` fixture outputs.
Run `python3 data/make_slim_bruker_dia_fixture.py --help` for the Bruker DIA
slice fixture generator.
