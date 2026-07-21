# timsTOF2 HeLa diaPASEF slice fixture

This directory contains a slim Bruker `.d` fixture generated from a local
timsTOF Ultra HeLa diaPASEF run. It is intended for opt-in native Bruker
`mzio dia-slice` integration checks.

The current fixture keeps 100 source frames around a DIA-NN precursor:

- peptide: `IILDLISESPIK/2`
- precursor m/z: `670.905411`
- source RT window: `21.6007-21.7765 min`
- source frames: `12156-12255`
- output `.d`: `bruker_d/hela_dia_iildlisespik_rt21_60_21_78.d`

The `.d` frame IDs are renumbered to `1..N`; the source frame IDs are recorded
in the JSON manifest next to the fixture. Regenerate from a local full raw run
with:

```bash
python3 data/make_slim_bruker_dia_fixture.py \
  --source-d data/timstof2_hela/DIA/99992_111_TOF_HeLa_DIA_20ng_Apr18_calib_610.d \
  --out-d data/timstof2_hela_dia_slice/bruker_d/hela_dia_iildlisespik_rt21_60_21_78.d \
  --frame-start 12156 \
  --frame-end 12255 \
  --target-label IILDLISESPIK/2 \
  --overwrite
```

Validation command:

```bash
mzio dia-slice \
  --bruker data/timstof2_hela_dia_slice/bruker_d/hela_dia_iildlisespik_rt21_60_21_78.d \
  --bruker-backend native \
  --peptide IILDLISESPIK/2 \
  --pseudo-ms2 \
  --trace-peaks \
  --emit-trace \
  --mz-da 0.025 \
  --rt-min 21.60 \
  --rt-max 21.78 \
  --outdir /tmp/mzio-slim-dia-smoke \
  --out-prefix slim_iildlisespik
```
