#!/usr/bin/env bash
set -euo pipefail

# =========================
# Config
# =========================
BASE_DIR="/data/supeng/env/env/horse/lfmm_vcf"
MAP_DIR="${BASE_DIR}"                         # chr*.map
ZSCORE_DIR="${BASE_DIR}/lfmm_bio1_results"    # chr*_s1.6.zscore
OUT_DIR="${ZSCORE_DIR}/merged"
OUT_FILE="${OUT_DIR}/all_chromosomes_LFMM.tsv"

SUFFIX="_s1.6.zscore"
CHR_START=1
CHR_END=31

mkdir -p "${OUT_DIR}"

# Header:
# PLINK .map format: CHR SNP_ID GD BP
# LFMM .zscore columns: Z  -log10(P)  P
echo -e "CHR\tSNP\tBP\tZ\tNEGLOG10P\tP" > "${OUT_FILE}"

# =========================
# Merge
# =========================
for chr in $(seq ${CHR_START} ${CHR_END}); do
  map_file="${MAP_DIR}/chr${chr}.map"
  z_file="${ZSCORE_DIR}/chr${chr}${SUFFIX}"

  if [[ ! -s "${map_file}" ]]; then
    echo "[ERROR] Missing or empty map: ${map_file}" >&2
    exit 1
  fi
  if [[ ! -s "${z_file}" ]]; then
    echo "[ERROR] Missing or empty zscore: ${z_file}" >&2
    exit 1
  fi

  # Line count check
  n_map=$(wc -l < "${map_file}" | tr -d ' ')
  n_z=$(wc -l < "${z_file}" | tr -d ' ')

  if [[ "${n_map}" -ne "${n_z}" ]]; then
    echo "[ERROR] Line count mismatch chr${chr}: map=${n_map}, zscore=${n_z}" >&2
    echo "        map: ${map_file}" >&2
    echo "        z  : ${z_file}" >&2
    exit 1
  fi

  # Sanity: ensure zscore has at least 3 columns on first line
  z_cols=$(awk 'NR==1{print NF}' "${z_file}")
  if [[ "${z_cols}" -lt 3 ]]; then
    echo "[ERROR] chr${chr} zscore seems to have <3 columns (NF=${z_cols}): ${z_file}" >&2
    echo "        First line:" >&2
    head -n 1 "${z_file}" >&2
    exit 1
  fi

  # Merge:
  # map -> CHR SNP BP (fields 1,2,4)
  # zscore -> Z NEGLOG10P P (fields 1,2,3)
  paste \
    <(awk 'BEGIN{OFS="\t"}{print $1,$2,$4}' "${map_file}") \
    <(awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "${z_file}") \
    >> "${OUT_FILE}"

  echo "[OK] chr${chr} merged: ${n_map} SNPs"
done

echo "=========================================="
echo "[DONE] Merged file:"
echo "  ${OUT_FILE}"
echo "Rows (including header): $(wc -l < "${OUT_FILE}")"
