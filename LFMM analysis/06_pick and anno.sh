#!/usr/bin/env bash
set -euo pipefail

# =========================
# Input / Output
# =========================
INPUT_FILE="/data/supeng/env/env/horse/lfmm_vcf/lfmm_bio1_results/merged/all_chromosomes_LFMM.tsv"
OUT_DIR="/data/supeng/env/env/horse/lfmm_vcf/lfmm_bio1_results/merged/significant_bonferroni"
ANNOT_BED="/data/supeng/bodysize/horse/gwas/852/825_maf005/output/horse3-gene_nochr.bed"

ALPHA="0.05"   # Bonferroni family-wise error rate

mkdir -p "${OUT_DIR}"

OUT_BED="${OUT_DIR}/significant_snps.bed"
OUT_ANN="${OUT_DIR}/significant_snps_annotated.tsv"
OUT_SORTED="${OUT_DIR}/sorted_significant_snps_annotated.tsv"
GENE_LIST="${OUT_DIR}/top_significant_genes.txt"
META="${OUT_DIR}/run_info.txt"

# =========================
# Sanity checks
# =========================
if [[ ! -s "${INPUT_FILE}" ]]; then
  echo "[ERROR] Missing input: ${INPUT_FILE}" >&2
  exit 1
fi
if [[ ! -s "${ANNOT_BED}" ]]; then
  echo "[ERROR] Missing annotation bed: ${ANNOT_BED}" >&2
  exit 1
fi
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found in PATH" >&2; exit 1; }
command -v bc >/dev/null 2>&1 || { echo "[ERROR] bc not found in PATH" >&2; exit 1; }

# =========================
# 1) Compute SNP count (exclude header) and Bonferroni threshold
# =========================
total_lines=$(wc -l < "${INPUT_FILE}" | tr -d ' ')
if [[ "${total_lines}" -le 1 ]]; then
  echo "[ERROR] Input seems empty (only header?)" >&2
  exit 1
fi

snp_count=$(( total_lines - 1 ))
bonf=$(echo "scale=20; ${ALPHA} / ${snp_count}" | bc -l)

{
  echo "INPUT_FILE=${INPUT_FILE}"
  echo "ANNOT_BED=${ANNOT_BED}"
  echo "ALPHA=${ALPHA}"
  echo "SNP_COUNT=${snp_count}"
  echo "BONFERRONI_THRESHOLD=${bonf}"
  echo "DATE=$(date '+%Y-%m-%d %H:%M:%S')"
} > "${META}"

echo "[INFO] SNP count = ${snp_count}"
echo "[INFO] Bonferroni threshold = ${bonf}"

# =========================
# 2) Harmonize chromosome naming with annotation bed (chr prefix or not)
#    Detect whether ANNOT_BED uses 'chr' prefix.
# =========================
annot_has_chr=0
if awk 'BEGIN{has=0} NR<=1000 {if($1 ~ /^chr/){has=1; exit}} END{exit !has}' "${ANNOT_BED}"; then
  annot_has_chr=1
fi

# =========================
# 3) Filter significant SNPs using P (col6) and write BED
#    INPUT columns: CHR SNP BP Z NEGLOG10P P
#    BED: chrom  start  end  SNP  Z  NEGLOG10P  P
#    start = BP-1 (0-based), end = BP
# =========================
awk -v thr="${bonf}" -v annot_chr="${annot_has_chr}" '
BEGIN{
  OFS="\t"
}
NR==1{next}  # skip header
{
  chr=$1; snp=$2; bp=$3; z=$4; nlogp=$5; p=$6;

  # add/remove "chr" to match annotation
  if(annot_chr==1 && chr !~ /^chr/) chr="chr"chr;
  if(annot_chr==0 && chr ~ /^chr/) sub(/^chr/,"",chr);

  # filter
  if(p+0 < thr+0){
    start=bp-1; end=bp;
    if(start<0) start=0;
    print chr, start, end, snp, z, nlogp, p
  }
}
' "${INPUT_FILE}" > "${OUT_BED}"

sig_n=$(wc -l < "${OUT_BED}" | tr -d ' ')
echo "[INFO] Significant SNPs (Bonferroni) = ${sig_n}"
echo "[INFO] BED written: ${OUT_BED}"

if [[ "${sig_n}" -eq 0 ]]; then
  echo "[WARN] No significant SNPs under Bonferroni. You may consider FDR (BH) or a relaxed threshold." >&2
  exit 0
fi

# =========================
# 4) Annotate with bedtools
# =========================
# -a: SNP bed (7 cols)
# -b: gene bed
# output: SNP cols + gene-bed cols
bedtools intersect -a "${OUT_BED}" -b "${ANNOT_BED}" -wa -wb > "${OUT_ANN}"
echo "[INFO] Annotated output: ${OUT_ANN}"

# =========================
# 5) Sort by P (SNP P is col7 in OUT_BED, so in OUT_ANN it's still col7)
# =========================
sort -k7,7g "${OUT_ANN}" > "${OUT_SORTED}"
echo "[INFO] Sorted by P: ${OUT_SORTED}"

# =========================
# 6) Extract gene list (assume gene name/id is last column of annotation bed)
# =========================
awk '{print $NF}' "${OUT_ANN}" | sort -u > "${GENE_LIST}"
echo "[INFO] Gene list: ${GENE_LIST}"

echo "=========================================="
echo "[DONE] All outputs in: ${OUT_DIR}"
echo " - ${OUT_BED}"
echo " - ${OUT_ANN}"
echo " - ${OUT_SORTED}"
echo " - ${GENE_LIST}"
echo " - ${META}"
