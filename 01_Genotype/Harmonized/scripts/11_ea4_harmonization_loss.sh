#!/usr/bin/env bash
# =============================================================================
# 11_ea4_harmonization_loss.sh — Quantify EA4 SNPs lost to harmonization
# =============================================================================
# Run on the cluster from PROJ_ROOT. Produces a table with, per cohort:
#
#   (A) EA4 weight file total SNPs        — all candidates SBayesR produced
#   (B) In cohort's pre-intersection bfile — what would score without harm.
#   (C) In 4-cohort intersection (no Twin) — what would score w/o TwinLife
#   (D) In 5-cohort shared intersection    — what actually scored (our pipeline)
#
# Exclusions due to:
#   (A - B) = this cohort's imputation/platform missingness
#   (B - C) = SOEP mildQC bottleneck (small)
#   (C - D) = TwinLife 1000G bottleneck
#   (A - D) = total harmonization loss
#
# USAGE:
#   bash scripts/11_ea4_harmonization_loss.sh
#   (Tardis-side, on-demand diagnostic; not part of the default DAG)
#
# INPUTS:
#   ${WEIGHTS_OUT}/EA4_excl_BASEII_sbayesR_chrpos.txt
#   ${WEIGHTS_OUT}/EA4_excl_SHIP_sbayesR_chrpos.txt
#   ${OUT_ROOT}/geno/<COHORT>_harmonized.bim         — pre-intersection bfile
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.bim  — 5-cohort intersection
#
# OUTPUTS:
#   ${PROJ_ROOT}/output/ea4_harmonization_loss.tsv   — per-cohort funnel table
#
# DEPENDS ON: steps 02 (weights) and 03 (harmonized + shared bfiles).
# =============================================================================
set -euo pipefail

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SUBMIT_DIR}/config.sh"

OUTFILE="${PROJ_ROOT}/output/ea4_harmonization_loss.tsv"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

WEIGHTS_EA4_BASEII="${WEIGHTS_OUT}/EA4_excl_BASEII_sbayesR_chrpos.txt"
WEIGHTS_EA4_SHIP="${WEIGHTS_OUT}/EA4_excl_SHIP_sbayesR_chrpos.txt"

# ---- helper: extract SNP ID column (col 1, header "SNP") -------------------
weight_snps () {
  awk 'NR>1 {print $1}' "$1" | sort -u
}

# ---- (A) Total weight SNPs -------------------------------------------------
A_BASEII=$(weight_snps "$WEIGHTS_EA4_BASEII" | wc -l)
A_SHIP=$(  weight_snps "$WEIGHTS_EA4_SHIP"   | wc -l)

# ---- Build per-cohort SNP lists at each funnel stage -----------------------
# Pre-intersection bfile: {cohort}_harmonized (no _shared suffix).
# Post-R2 is already applied here; only the 5-cohort intersection differs.
declare -A PRE_BIM
PRE_BIM[BASEII]="${GENO_OUT}/BASEII_harmonized.bim"
PRE_BIM[SHIP0]="${GENO_OUT}/SHIP0_harmonized.bim"
PRE_BIM[SHIPTD]="${GENO_OUT}/SHIPTD_harmonized.bim"
PRE_BIM[SOEP]="${GENO_OUT}/SOEP_harmonized.bim"
PRE_BIM[TWINLIFE]="${GENO_OUT}/TWINLIFE_harmonized.bim"

declare -A POST_BIM
POST_BIM[BASEII]="${GENO_OUT}/BASEII_harmonized_shared.bim"
POST_BIM[SHIP0]="${GENO_OUT}/SHIP0_harmonized_shared.bim"
POST_BIM[SHIPTD]="${GENO_OUT}/SHIPTD_harmonized_shared.bim"
POST_BIM[SOEP]="${GENO_OUT}/SOEP_harmonized_shared.bim"
POST_BIM[TWINLIFE]="${GENO_OUT}/TWINLIFE_harmonized_shared.bim"

bim_snps () { awk '{print $2}' "$1" | sort -u; }

# 4-cohort intersection (no TwinLife):
#   intersect(BASEII_pre, SHIP0_pre, SHIPTD_pre, SOEP_pre)
FOUR="${TMP}/four_cohort_shared.txt"
bim_snps "${PRE_BIM[BASEII]}" > "${TMP}/B.txt"
bim_snps "${PRE_BIM[SHIP0]}"  > "${TMP}/S0.txt"
bim_snps "${PRE_BIM[SHIPTD]}" > "${TMP}/ST.txt"
bim_snps "${PRE_BIM[SOEP]}"   > "${TMP}/SO.txt"
comm -12 "${TMP}/B.txt"   "${TMP}/S0.txt" > "${TMP}/x1.txt"
comm -12 "${TMP}/x1.txt"  "${TMP}/ST.txt" > "${TMP}/x2.txt"
comm -12 "${TMP}/x2.txt"  "${TMP}/SO.txt" > "$FOUR"
FOUR_N=$(wc -l < "$FOUR")

# ---- scoring intersection per cohort ---------------------------------------
# For a given weight file W and cohort bim B, scoreable SNPs =
# SNPs that appear in both (sort-unique, comm -12).
count_intersect () {
  local w="$1"; local b="$2"
  comm -12 <(weight_snps "$w") <(bim_snps "$b") | wc -l
}

# weight × 4-cohort intersection file
count_intersect_file () {
  local w="$1"; local f="$2"
  comm -12 <(weight_snps "$w") <(sort -u "$f") | wc -l
}

echo "Computing intersections (this takes 1-2 min per cohort)..."

# pick the EA4 weight file whose SNP set is attributed for each cohort
declare -A WEIGHT
WEIGHT[BASEII]="$WEIGHTS_EA4_BASEII"
WEIGHT[SHIP0]="$WEIGHTS_EA4_SHIP"
WEIGHT[SHIPTD]="$WEIGHTS_EA4_SHIP"
WEIGHT[SOEP]="$WEIGHTS_EA4_SHIP"
WEIGHT[TWINLIFE]="$WEIGHTS_EA4_SHIP"

declare -A A_TOT B_PRE C_4 D_5
for c in BASEII SHIP0 SHIPTD SOEP TWINLIFE; do
  w="${WEIGHT[$c]}"
  A_TOT[$c]=$(weight_snps "$w" | wc -l)
  B_PRE[$c]=$(count_intersect "$w" "${PRE_BIM[$c]}")
  # TwinLife cannot gain from a 4-cohort intersection that excludes itself
  if [[ "$c" == "TWINLIFE" ]]; then
    C_4[$c]="NA"
  else
    C_4[$c]=$(count_intersect_file "$w" "$FOUR")
  fi
  D_5[$c]=$(count_intersect "$w" "${POST_BIM[$c]}")
  printf "  %-8s A=%d  B=%d  C=%s  D=%d\n" \
    "$c" "${A_TOT[$c]}" "${B_PRE[$c]}" "${C_4[$c]}" "${D_5[$c]}"
done

# ---- write result ----------------------------------------------------------
{
  echo -e "cohort\tweight_file\tA_weight_total\tB_pre_intersect\tC_4cohort_no_twin\tD_5cohort_shared\tloss_B_to_D\tloss_C_to_D_twin_bottleneck"
  for c in BASEII SHIP0 SHIPTD SOEP TWINLIFE; do
    w=$(basename "${WEIGHT[$c]}")
    a="${A_TOT[$c]}"; b="${B_PRE[$c]}"; cc="${C_4[$c]}"; d="${D_5[$c]}"
    loss_bd=$(( b - d ))
    if [[ "$cc" == "NA" ]]; then loss_cd="NA"; else loss_cd=$(( cc - d )); fi
    echo -e "$c\t$w\t$a\t$b\t$cc\t$d\t$loss_bd\t$loss_cd"
  done
} > "$OUTFILE"

echo
echo "Result written: $OUTFILE"
column -s $'\t' -t "$OUTFILE"
echo
echo "4-cohort intersection size (no TwinLife): $FOUR_N"
echo "5-cohort shared intersection size:         $(wc -l < "${POST_BIM[BASEII]}")"
