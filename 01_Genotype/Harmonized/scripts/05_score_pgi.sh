#!/usr/bin/env bash
# =============================================================================
# 05_score_pgi.sh — Compute PGIs using PLINK2 --score
# =============================================================================
# For each cohort, score all relevant traits using SBayesR weight files.
# Weight-to-cohort mapping (see config.sh helpers — authoritative):
#   EA4:    all cohorts → EA4_excl_SHIP        (get_ea4_trait_for_cohort)
#   NonCog: all cohorts → NCog (Malanchini)    (get_noncog_trait_for_cohort)
#   Cog:    same for all
#   Height: same for all
#
# USAGE:
#   bash scripts/05_score_pgi.sh
#   (typically via `bash run_pipeline.sh --step 05`)
#
# INPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}
#   ${WEIGHTS_OUT}/<TAG>_sbayesR_chrpos.txt  — scoring files from step 02
#
# OUTPUTS:
#   ${OUT_ROOT}/scores/<COHORT>_<trait>.sscore  — PLINK2 BETA_SUM scores
#
# DEPENDS ON: steps 02 (weights) and 03 (harmonized bfiles). Independent of
# step 04/04b — PC computation can run in parallel.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 05: PGI Scoring with PLINK2 ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
SCORE_SUMMARY=""

# --- Helper: score one cohort with one weight file ---------------------------
score_one() {
    local cohort="$1"
    local trait_tag="$2"
    local weight_file="$3"
    local extra_args="$4"  # e.g., "--maf 0.01" or ""

    local bfile
    bfile=$(get_harmonized_bfile "$cohort")
    local out_prefix="${SCORES_OUT}/${cohort}_${trait_tag}"
    local sscore="${out_prefix}.sscore"

    # Skip if already done
    if [[ -f "$sscore" ]]; then
        log_info "  SKIP: $cohort $trait_tag (already scored)"
        return 0
    fi

    if [[ ! -f "$weight_file" ]]; then
        log_error "  Weight file not found: $weight_file"
        return 1
    fi

    log_info "  Scoring: $cohort x $trait_tag"

    # Build PLINK2 command
    local cmd=("$PLINK2"
        --bfile "$bfile"
        --score "$weight_file" 1 2 3 header-read cols=+scoresums
        --threads "$PLINK2_THREADS"
        --out "$out_prefix"
    )

    # Add extra args (e.g., MAF filter)
    if [[ -n "$extra_args" ]]; then
        # shellcheck disable=SC2206
        cmd+=($extra_args)
    fi

    "${cmd[@]}" 2>&1 | tee "${LOG_DIR}/score_${cohort}_${trait_tag}.log"

    if [[ -f "$sscore" ]]; then
        local n_scored
        n_scored=$(grep -c "variants processed" "${out_prefix}.log" 2>/dev/null || echo "?")
        log_info "    -> Score written to $sscore"
    else
        log_error "    -> Scoring FAILED for $cohort $trait_tag"
    fi
}

# --- Score all cohorts --------------------------------------------------------
for cohort in "${COHORTS[@]}"; do
    log_info ""
    log_info "--- Scoring cohort: $cohort ---"

    bfile=$(get_harmonized_bfile "$cohort")
    if [[ ! -f "${bfile}.bed" ]]; then
        log_warn "Harmonized genotype not found for $cohort — skipping"
        continue
    fi

    # EA4 (single common weight file for all cohorts)
    ea4_trait=$(get_ea4_trait_for_cohort "$cohort")
    ea4_tag="${TRAIT_SBAYESR_TAG[$ea4_trait]}"
    ea4_weight="${WEIGHTS_OUT}/${ea4_tag}_sbayesR_chrpos.txt"
    score_one "$cohort" "EA4" "$ea4_weight" ""

    # Cognitive (same for all)
    cog_weight="${WEIGHTS_OUT}/Cog_sbayesR_chrpos.txt"
    score_one "$cohort" "Cog" "$cog_weight" ""

    # Non-Cognitive (published weights, all cohorts)
    ncog_trait=$(get_noncog_trait_for_cohort "$cohort")
    ncog_tag="${TRAIT_SBAYESR_TAG[$ncog_trait]}"
    ncog_weight="${WEIGHTS_OUT}/${ncog_tag}_sbayesR_chrpos.txt"
    score_one "$cohort" "NonCog" "$ncog_weight" ""

    # Height (same for all)
    height_weight="${WEIGHTS_OUT}/Height_sbayesR_chrpos.txt"
    score_one "$cohort" "Height" "$height_weight" ""

    SCORE_SUMMARY="${SCORE_SUMMARY}\n- $cohort: EA4=${ea4_trait}, Cog=Cog, NonCog=${ncog_trait}, Height=Height"
done

# --- Log decisions -----------------------------------------------------------
{
    echo ""
    echo "## Step 05: PGI Scoring ($(date '+%Y-%m-%d %H:%M'))"
    echo ""
    echo "- Scored with PLINK2 --score using chr:pos matched weights"
    echo "- Score output: cols=+scoresums (A1Effect_SUM column)"
    echo "- Variant IDs: chr:pos format (matching step 03 harmonization)"
    echo ""
    echo "### Weight-to-cohort mapping"
    echo -e "$SCORE_SUMMARY"
    echo ""
} >> "$DECISION_LOG"

log_info "=== Step 05 complete ==="
