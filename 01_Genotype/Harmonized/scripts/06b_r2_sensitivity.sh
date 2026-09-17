#!/usr/bin/env bash
# =============================================================================
# 06b_r2_sensitivity.sh — Score PGIs with strict R² > 0.8 SNP set
# =============================================================================
# Produces a second set of PGI scores using only variants that pass the
# strict imputation R² threshold (default 0.8) in ALL cohorts.
#
# This is a sensitivity analysis: the main pipeline uses R² >= 0.3.
# Comparing the two sets shows whether poorly-imputed variants drive
# cross-cohort PGI differences.
#
# USAGE:
#   bash scripts/06b_r2_sensitivity.sh
#   (typically via `bash run_pipeline.sh --step 06b`; skipped
#   when R2_FILTER=no or the strict shared SNP set is absent)
#
# INPUTS:
#   ${OUT_ROOT}/weights/_shared_snp_set_r2_${R2_THRESHOLD_STRICT}.txt
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}
#   ${WEIGHTS_OUT}/<TAG>_sbayesR_chrpos.txt
#
# OUTPUTS:
#   ${OUT_ROOT}/scores/<COHORT>_<trait>_r2strict.sscore
#
# DEPENDS ON: steps 02 (weights), 03 (harmonized bfiles + strict shared SNP
# set), 03a (strict R² pass lists).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

if [[ "${R2_FILTER:-yes}" != "yes" ]]; then
  log_info "=== Step 06b: SKIPPED (R2_FILTER=no) ==="
  exit 0
fi

STRICT_SNP_SET="${OUT_ROOT}/weights/_shared_snp_set_r2_${R2_THRESHOLD_STRICT}.txt"

if [[ ! -f "$STRICT_SNP_SET" ]]; then
  log_warn "=== Step 06b: SKIPPED (strict SNP set not found at $STRICT_SNP_SET) ==="
  log_warn "  This can happen if strict R2 pass lists were missing for some cohorts."
  exit 0
fi

N_STRICT=$(wc -l < "$STRICT_SNP_SET")
log_info "=== Step 06b: R² Sensitivity Scoring (R2 >= ${R2_THRESHOLD_STRICT}) ==="
log_info "  Strict shared SNP set: $N_STRICT SNPs"

DECISION_LOG="${PROJ_ROOT}/decision_log.md"

# --- Helper: score one cohort with one weight file, restricted to strict set --
score_strict() {
  local cohort="$1"
  local trait_tag="$2"
  local weight_file="$3"

  local bfile
  bfile=$(get_harmonized_bfile "$cohort")
  local out_prefix="${SCORES_OUT}/${cohort}_${trait_tag}_r2strict"
  local sscore="${out_prefix}.sscore"

  # Skip if already done
  if [[ -f "$sscore" ]]; then
    log_info "  SKIP: $cohort $trait_tag r2strict (already scored)"
    return 0
  fi

  if [[ ! -f "$weight_file" ]]; then
    log_error "  Weight file not found: $weight_file"
    return 1
  fi

  log_info "  Scoring: $cohort x $trait_tag (R2 strict)"

  "$PLINK2" \
    --bfile "$bfile" \
    --extract "$STRICT_SNP_SET" \
    --score "$weight_file" 1 2 3 header-read cols=+scoresums \
    --threads "$PLINK2_THREADS" \
    --out "$out_prefix" \
    2>&1 | tee "${LOG_DIR}/score_${cohort}_${trait_tag}_r2strict.log"

  if [[ -f "$sscore" ]]; then
    log_info "    -> Score written to $sscore"
  else
    log_error "    -> Scoring FAILED for $cohort $trait_tag (r2strict)"
  fi
}

# --- Score all cohorts --------------------------------------------------------
for cohort in "${COHORTS[@]}"; do
  log_info ""
  log_info "--- Scoring cohort: $cohort (R2 strict) ---"

  bfile=$(get_harmonized_bfile "$cohort")
  if [[ ! -f "${bfile}.bed" ]]; then
    log_warn "Harmonized genotype not found for $cohort — skipping"
    continue
  fi

  # EA4 (single common weight file for all cohorts)
  ea4_trait=$(get_ea4_trait_for_cohort "$cohort")
  ea4_tag="${TRAIT_SBAYESR_TAG[$ea4_trait]}"
  ea4_weight="${WEIGHTS_OUT}/${ea4_tag}_sbayesR_chrpos.txt"
  score_strict "$cohort" "EA4" "$ea4_weight"

  # Cognitive
  cog_weight="${WEIGHTS_OUT}/Cog_sbayesR_chrpos.txt"
  score_strict "$cohort" "Cog" "$cog_weight"

  # Non-Cognitive (published weights, all cohorts)
  ncog_trait=$(get_noncog_trait_for_cohort "$cohort")
  ncog_tag="${TRAIT_SBAYESR_TAG[$ncog_trait]}"
  ncog_weight="${WEIGHTS_OUT}/${ncog_tag}_sbayesR_chrpos.txt"
  score_strict "$cohort" "NonCog" "$ncog_weight"

  # Height
  height_weight="${WEIGHTS_OUT}/Height_sbayesR_chrpos.txt"
  score_strict "$cohort" "Height" "$height_weight"
done

# --- Log decisions -----------------------------------------------------------
{
  echo ""
  echo "## Step 06b: R² Sensitivity Scoring ($(date '+%Y-%m-%d %H:%M'))"
  echo ""
  echo "- Scored all cohorts x traits with strict R2 >= ${R2_THRESHOLD_STRICT} SNP set"
  echo "- Strict shared SNP set: $N_STRICT SNPs"
  echo "- Output: output/scores/{COHORT}_{trait}_r2strict.sscore"
  echo "- Compare with main scores (R2 >= ${R2_THRESHOLD}) to assess robustness"
  echo ""
} >> "$DECISION_LOG"

log_info ""
log_info "=== Step 06b complete ==="
