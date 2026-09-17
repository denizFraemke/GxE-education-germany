#!/usr/bin/env bash
# =============================================================================
# 03_harmonize_geno.sh — Harmonize genotype variant IDs to chr:pos format
#                     + restrict to unified SNP set for cross-cohort comparability
# =============================================================================
# PHASE 1: Basic harmonization
#   For each cohort:
#     1. Convert variant IDs to chr{CHR}:{POS} using PLINK2 --set-all-var-ids
#     2. Remove duplicate positions (--rm-dup force-first)
#     3. Output: {COHORT}_harmonized.bed/bim/fam
#
# PHASE 2: Unified SNP set restriction (if SHARED_SNP_SET=yes)
#   After all cohorts are harmonized:
#     1. Compute intersection of SNPs across all 5 cohorts' genotypes
#        (SNPs must be present in ALL cohorts to be scorable for every cohort)
#     2. Optional: report which of those also have weight SNPs available
#     3. Restrict each cohort to the shared SNP set
#     4. Output: {COHORT}_harmonized_shared.bed/bim/fam
#
# NOTE: We do NOT restrict to SNPs present in ALL weight files.
#       Each trait is scored independently — only the cross-cohort genotype
#       intersection matters for comparability.
#
# USAGE:
#   bash scripts/03_harmonize_geno.sh
#   (typically via `bash run_pipeline.sh --step 03`)
#
# INPUTS:
#   $(get_geno_bfile <cohort>).{bed,bim,fam}                   — raw bfiles
#   ${OUT_ROOT}/geno/<COHORT>_r2_pass_${R2_THRESHOLD}.txt       — R² pass list
#   ${OUT_ROOT}/geno/<COHORT>_r2_pass_${R2_THRESHOLD_STRICT}.txt — strict R² list
#
# OUTPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized.{bed,bim,fam}         — phase 1 output
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}  — phase 2 output
#                                                                (used downstream)
#   ${OUT_ROOT}/weights/_shared_snp_set_r2_${R2_THRESHOLD}.txt
#   ${OUT_ROOT}/weights/_shared_snp_set_r2_${R2_THRESHOLD_STRICT}.txt
#
# DEPENDS ON: step 03a (R² pass lists) when R2_FILTER=yes; otherwise step 00.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

# Harmonization settings (from config or defaults)
SHARED_SNP_SET="${SHARED_SNP_SET:-yes}"

log_info "=== Step 03: Harmonizing genotype data ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
SUMMARY=""
PHASE1_DONE=""

# =============================================================================
# PHASE 1: Basic harmonization (chr:pos IDs, dedup)
# =============================================================================
for cohort in "${COHORTS[@]}"; do
    bfile=$(get_geno_bfile "$cohort")
    out_prefix=$(get_harmonized_bfile "$cohort")

    # Skip if already done
    if [[ -f "${out_prefix}.bed" ]]; then
        n_var=$(wc -l < "${out_prefix}.bim")
        log_info "SKIP: $cohort already harmonized ($n_var variants)"
        SUMMARY="${SUMMARY}\n- $cohort [phase1]: SKIPPED (already done, $n_var variants)"
        PHASE1_DONE="${PHASE1_DONE} ${cohort}"
        continue
    fi

    log_info "Harmonizing $cohort: $bfile -> $out_prefix"

    # Count input
    n_in=$(wc -l < "${bfile}.bim")
    log_info "  Input: $n_in variants"

    # PLINK2: set variant IDs to chr:pos, remove duplicates, make bed
    "$PLINK2" \
        --bfile "$bfile" \
        --set-all-var-ids "${PLINK2_SET_VAR_IDS}" \
        --rm-dup "$PLINK2_RM_DUP" \
        --make-bed \
        --threads "$PLINK2_THREADS" \
        --out "$out_prefix" \
        2>&1 | tee "${LOG_DIR}/harmonize_${cohort}.log"

    # Count output
    n_out=$(wc -l < "${out_prefix}.bim")
    n_dup=$((n_in - n_out))
    log_info "  Output: $n_out variants ($n_dup duplicates removed)"
    SUMMARY="${SUMMARY}\n- $cohort [phase1]: $n_in -> $n_out variants ($n_dup duplicates removed)"
    PHASE1_DONE="${PHASE1_DONE} ${cohort}"
done

# =============================================================================
# PHASE 2: Unified SNP set restriction (cross-cohort comparability)
# =============================================================================
if [[ "$SHARED_SNP_SET" == "yes" ]]; then
    log_info ""
    log_info "=== PHASE 2: Computing cross-cohort SNP intersection ==="

    # Step 2a: Compute global intersection of SNPs across all 5 cohorts' genotypes.
    # Uses sort+comm for fast set intersection (much faster than awk for large files).
    log_info "Computing SNP intersection across all 5 cohorts' genotypes..."

    # Extract sorted SNP IDs from each cohort's bim
    FIRST_COH="${COHORTS[0]}"
    first_bfile=$(get_harmonized_bfile "$FIRST_COH")
    if [[ ! -f "${first_bfile}.bim" ]]; then
        log_error "  ${FIRST_COH} harmonized bim not found — cannot compute intersection"
        exit 1
    fi
    awk '{print $2}' "${first_bfile}.bim" | sort > "${OUT_ROOT}/weights/_tmp_sorted.txt"
    log_info "  ${FIRST_COH}: $(wc -l < "${first_bfile}.bim") SNPs (starting set)"

    # Intersect with each remaining cohort using sort+comm
    for cohort in "${COHORTS[@]:1}"; do
        bfile=$(get_harmonized_bfile "$cohort")
        if [[ ! -f "${bfile}.bim" ]]; then
            log_warn "  $cohort: harmonized bim not found — skipping"
            continue
        fi
        awk '{print $2}' "${bfile}.bim" | sort > "${OUT_ROOT}/weights/_tmp_cohort_sorted.txt"
        comm -12 "${OUT_ROOT}/weights/_tmp_sorted.txt" "${OUT_ROOT}/weights/_tmp_cohort_sorted.txt" \
            > "${OUT_ROOT}/weights/_tmp_intersect.txt"
        mv "${OUT_ROOT}/weights/_tmp_intersect.txt" "${OUT_ROOT}/weights/_tmp_sorted.txt"
        rm -f "${OUT_ROOT}/weights/_tmp_cohort_sorted.txt"
        log_info "  ${cohort}: $(wc -l < "${OUT_ROOT}/weights/_tmp_sorted.txt") SNPs after intersection"
    done

    SHARED_SNP_SET_RAW="${OUT_ROOT}/weights/_shared_snp_set_raw.txt"
    mv "${OUT_ROOT}/weights/_tmp_sorted.txt" "$SHARED_SNP_SET_RAW"
    N_SHARED_RAW=$(wc -l < "$SHARED_SNP_SET_RAW")
    log_info "  Global SNP intersection (genotype only): $N_SHARED_RAW SNPs"

    # Step 2a-R2: Apply imputation R² filtering to the shared SNP set
    # Uses sort+comm for fast set intersection
    if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
        log_info ""
        log_info "Applying imputation R² filtering to shared SNP set..."

        # --- Lenient threshold (main pipeline) ---
        SHARED_SNP_SET="${OUT_ROOT}/weights/_shared_snp_set.txt"
        # shared set is already sorted from above
        cp "$SHARED_SNP_SET_RAW" "${OUT_ROOT}/weights/_tmp_r2.txt"
        for cohort in "${COHORTS[@]}"; do
            r2_pass="${GENO_OUT}/${cohort}_r2_pass_${R2_THRESHOLD}.txt"
            if [[ ! -f "$r2_pass" ]]; then
                log_error "  R2 pass list not found for $cohort: $r2_pass"
                log_error "  Run step 03a first."
                exit 1
            fi
            sort "$r2_pass" > "${OUT_ROOT}/weights/_tmp_r2_sorted.txt"
            comm -12 "${OUT_ROOT}/weights/_tmp_r2.txt" "${OUT_ROOT}/weights/_tmp_r2_sorted.txt" \
                > "${OUT_ROOT}/weights/_tmp_r2b.txt"
            mv "${OUT_ROOT}/weights/_tmp_r2b.txt" "${OUT_ROOT}/weights/_tmp_r2.txt"
            rm -f "${OUT_ROOT}/weights/_tmp_r2_sorted.txt"
            log_info "    $cohort: $(wc -l < "${OUT_ROOT}/weights/_tmp_r2.txt") SNPs remaining"
        done
        mv "${OUT_ROOT}/weights/_tmp_r2.txt" "$SHARED_SNP_SET"
        N_SHARED=$(wc -l < "$SHARED_SNP_SET")
        N_R2_REMOVED=$((N_SHARED_RAW - N_SHARED))
        log_info "  After R2 >= ${R2_THRESHOLD} filter: $N_SHARED SNPs ($N_R2_REMOVED removed)"

        # --- Strict threshold (sensitivity analysis) ---
        SHARED_SNP_SET_STRICT="${OUT_ROOT}/weights/_shared_snp_set_r2_${R2_THRESHOLD_STRICT}.txt"
        cp "$SHARED_SNP_SET_RAW" "${OUT_ROOT}/weights/_tmp_r2s.txt"
        for cohort in "${COHORTS[@]}"; do
            r2_pass="${GENO_OUT}/${cohort}_r2_pass_${R2_THRESHOLD_STRICT}.txt"
            if [[ ! -f "$r2_pass" ]]; then
                log_warn "  Strict R2 pass list not found for $cohort — skipping sensitivity set"
                rm -f "${OUT_ROOT}/weights/_tmp_r2s.txt"
                SHARED_SNP_SET_STRICT=""
                break
            fi
            sort "$r2_pass" > "${OUT_ROOT}/weights/_tmp_r2s_sorted.txt"
            comm -12 "${OUT_ROOT}/weights/_tmp_r2s.txt" "${OUT_ROOT}/weights/_tmp_r2s_sorted.txt" \
                > "${OUT_ROOT}/weights/_tmp_r2sb.txt"
            mv "${OUT_ROOT}/weights/_tmp_r2sb.txt" "${OUT_ROOT}/weights/_tmp_r2s.txt"
            rm -f "${OUT_ROOT}/weights/_tmp_r2s_sorted.txt"
        done
        if [[ -n "$SHARED_SNP_SET_STRICT" && -f "${OUT_ROOT}/weights/_tmp_r2s.txt" ]]; then
            mv "${OUT_ROOT}/weights/_tmp_r2s.txt" "$SHARED_SNP_SET_STRICT"
            N_SHARED_STRICT=$(wc -l < "$SHARED_SNP_SET_STRICT")
            log_info "  After R2 >= ${R2_THRESHOLD_STRICT} filter: $N_SHARED_STRICT SNPs (sensitivity)"
        fi
    else
        SHARED_SNP_SET="${OUT_ROOT}/weights/_shared_snp_set.txt"
        cp "$SHARED_SNP_SET_RAW" "$SHARED_SNP_SET"
        N_SHARED=$(wc -l < "$SHARED_SNP_SET")
        log_info "  R2 filtering disabled — using raw intersection: $N_SHARED SNPs"
    fi

    # Step 2b: Report weight file coverage of the shared SNP set
    log_info ""
    log_info "Checking weight file coverage of shared SNPs..."
    ALL_WEIGHT_SNPS="${OUT_ROOT}/weights/_all_weight_snps.txt"
    : > "$ALL_WEIGHT_SNPS"
    for trait in "${TRAIT_NAMES[@]}"; do
        tag="${TRAIT_SBAYESR_TAG[$trait]}"
        weight_file="${WEIGHTS_OUT}/${tag}_sbayesR_chrpos.txt"
        if [[ -f "$weight_file" ]]; then
            n_snps=$(awk 'NR==FNR{a[$1];next} $1 in a' "$weight_file" "$SHARED_SNP_SET" | wc -l)
            log_info "  ${tag}: $n_snps of $N_SHARED shared SNPs have weights"
        fi
    done

    # Step 2c: QC check
    if [[ "$N_SHARED" -lt 100000 ]]; then
        log_warn "  WARNING: Shared SNP set has only $N_SHARED SNPs (< 100k recommended)"
        log_warn "  Consider imputation or relaxing SNP set requirements"
    fi

    # Step 2d: Restrict each cohort to shared SNP set
    log_info ""
    log_info "Restricting each cohort to shared SNP set..."
    for cohort in "${COHORTS[@]}"; do
        bfile=$(get_harmonized_bfile "$cohort")
        out_prefix="${GENO_OUT}/${cohort}_harmonized_shared"
        n_before=$(wc -l < "${bfile}.bim")

        if [[ -f "${out_prefix}.bed" ]]; then
            n_after=$(wc -l < "${out_prefix}.bim")
            log_info "  SKIP: $cohort already restricted ($n_after variants)"
            SUMMARY="${SUMMARY}\n- $cohort [phase2]: SKIPPED (already done, $n_after shared SNPs)"
            continue
        fi

        "$PLINK2" \
            --bfile "$bfile" \
            --extract "$SHARED_SNP_SET" \
            --make-bed \
            --threads "$PLINK2_THREADS" \
            --out "$out_prefix" \
            2>&1 | tee "${LOG_DIR}/harmonize_shared_${cohort}.log"

        n_after=$(wc -l < "${out_prefix}.bim")
        log_info "  $cohort: $n_before -> $n_after shared SNPs"
        SUMMARY="${SUMMARY}\n- $cohort [phase2]: $n_before -> $n_after shared SNPs (from $N_SHARED global set)"
    done

    # Log the shared SNP set path
    {
        echo ""
        echo "## Step 03 Phase 2: Unified SNP Set ($(date '+%Y-%m-%d %H:%M'))"
        echo ""
        echo "- SHARED_SNP_SET=yes: restricted all cohorts to SNPs present in ALL 5 cohorts"
        echo "- Genotype-only intersection: $N_SHARED_RAW SNPs"
        if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
            echo "- R2_FILTER=yes: additionally filtered by imputation R²"
            echo "- After R2 >= ${R2_THRESHOLD}: $N_SHARED SNPs (main pipeline)"
            if [[ -n "${SHARED_SNP_SET_STRICT:-}" && -f "${SHARED_SNP_SET_STRICT:-}" ]]; then
                echo "- After R2 >= ${R2_THRESHOLD_STRICT}: $N_SHARED_STRICT SNPs (sensitivity)"
            fi
        fi
        echo "- Shared SNP set file: $SHARED_SNP_SET"
        echo "- Restricted bfiles: {COHORT}_harmonized_shared.bed/bim/fam"
        echo ""
    } >> "$DECISION_LOG"

    log_info ""
    log_info "=== Step 03 complete (PHASE 1 + PHASE 2) ==="
else
    log_info ""
    log_info "=== Step 03 complete (PHASE 1 only — SHARED_SNP_SET=no) ==="
fi

# --- Log decisions -----------------------------------------------------------
{
    echo ""
    echo "## Step 03: Genotype Harmonization ($(date '+%Y-%m-%d %H:%M'))"
    echo ""
    echo "- Variant IDs converted to chr:pos format using PLINK2 --set-all-var-ids '${PLINK2_SET_VAR_IDS}'"
    echo "- Duplicate positions resolved with --rm-dup ${PLINK2_RM_DUP}"
    echo "- Decision: Using chr:pos (without alleles) for simpler matching with SBayesR weights"
    echo "- Decision: Keeping first variant at duplicate positions (force-first)"
    echo "- SHARED_SNP_SET=${SHARED_SNP_SET}"
    if [[ "$SHARED_SNP_SET" == "yes" ]]; then
        echo "- Phase 2: Restricted to $N_SHARED shared SNPs across all cohorts"
        echo "- Restricted bfiles: {COHORT}_harmonized_shared.bed/bim/fam"
        echo "- Shared SNP set: $SHARED_SNP_SET"
    fi
    echo ""
    echo "### Per-cohort summary"
    echo -e "$SUMMARY"
    echo ""
} >> "$DECISION_LOG"
