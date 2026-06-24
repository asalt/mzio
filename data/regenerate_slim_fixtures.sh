#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${MZIO_MSCONVERT_IMAGE:-chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest}"

run_msconvert() {
  local fixture="$1"
  local source_name="$2"
  local scan="$3"
  local output_name="$4"

  local source_dir="${repo_root}/data/_sources/${fixture}"
  local output_dir="${repo_root}/data/${fixture}/mzml"

  if [[ ! -f "${source_dir}/${source_name}" ]]; then
    printf 'missing source mzML: %s\n' "${source_dir}/${source_name}" >&2
    exit 1
  fi

  mkdir -p "${output_dir}"

  docker run --rm \
    -e WINEDEBUG=-all \
    -v "${source_dir}:/source:ro" \
    -v "${output_dir}:/out" \
    "${image}" \
    wine msconvert "/source/${source_name}" \
      --mzML \
      --zlib \
      --stripLocationFromSourceFiles \
      --stripVersionFromSoftware \
      --filter "scanNumber ${scan}" \
      --outfile "${output_name}" \
      --outdir /out

  perl -0pi -e \
    "s|/source/${source_name}|data/_sources/${fixture}/${source_name}|g; s|--outdir /out|--outdir data/${fixture}/mzml|g" \
    "${output_dir}/${output_name}"
}

run_msconvert \
  "exploris480_hela_digest" \
  "99990_236_EXP_Hela_100ng_200nl_Jan8.mzML" \
  "9065" \
  "99990_236_EXP_Hela_100ng_200nl_Jan8_scan9065.mzML"

run_msconvert \
  "exploris480_phospho_localization" \
  "49108_1_EXP_802_PDX_1mg_phos_oneforth_F13.mzML" \
  "4877" \
  "49108_1_EXP_802_PDX_1mg_phos_oneforth_F13_scan4877.mzML"
