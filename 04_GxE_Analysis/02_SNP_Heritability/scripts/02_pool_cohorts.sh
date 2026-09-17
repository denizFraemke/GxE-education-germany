#!/usr/bin/env bash
# =============================================================================
# 02_pool_cohorts.sh — Merge the 5 INFO>=.90 cohort bfiles into one pooled set
# =============================================================================
# plink 1.9 --merge-list with a missnp retry — same pattern as the upstream
# 04b_compute_cross_cohort_pcs.sh. Strand-ambiguous (A/T, C/G) variants can
# appear as 3+ alleles after harmonization; on the first merge failure we
# exclude the missnp from ALL cohorts and retry.
#
# Output: ${GENO_OUT}/pooled.{bed,bim,fam}
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 02: Pool ${#COHORTS[@]} cohort bfiles ==="

POOLED="${GENO_OUT}/pooled"
MERGE_LIST="${GENO_OUT}/_pool_merge_list.txt"

if [[ -f "${POOLED}.bed" ]]; then
    n_samples=$(wc -l < "${POOLED}.fam")
    n_variants=$(wc -l < "${POOLED}.bim")
    log_info "SKIP: pooled bfile already present (${n_samples} samples × ${n_variants} variants)"
    exit 0
fi

# Build merge list (every cohort except the first; first is passed via --bfile).
: > "$MERGE_LIST"
FIRST_BFILE=""
for cohort in "${COHORTS[@]}"; do
    bf="${GENO_OUT}/${cohort}_harmonized_r2_${R2_THRESHOLD_H2}"
    if [[ ! -f "${bf}.bed" ]]; then
        log_error "Missing INFO>=${R2_THRESHOLD_H2} bfile for ${cohort}: ${bf}.bed"
        log_error "Run step 01 first."
        exit 1
    fi
    if [[ -z "$FIRST_BFILE" ]]; then
        FIRST_BFILE="$bf"
    else
        echo "${bf}.bed ${bf}.bim ${bf}.fam" >> "$MERGE_LIST"
    fi
done

log_info "Merge anchor: $(basename "$FIRST_BFILE")"
log_info "Adding $(wc -l < "$MERGE_LIST") cohort(s) via merge-list."

MISSNP="${POOLED}-merge.missnp"
rm -f "${POOLED}".{bed,bim,fam,log,nosex} "$MISSNP"

log_info ""
log_info "Merge attempt 1..."
set +e
"$PLINK19" \
    --bfile "$FIRST_BFILE" \
    --merge-list "$MERGE_LIST" \
    --make-bed \
    --out "$POOLED" \
    2>&1 | tee "${LOG_DIR}/02_merge_attempt1.log"
exit1=$?
set -e

if [[ $exit1 -ne 0 ]] && [[ -f "$MISSNP" ]]; then
    n_missnp=$(wc -l < "$MISSNP")
    log_warn "Merge failed: $n_missnp 3+-allele variants. Excluding from all cohorts and retrying."

    CLEAN_LIST="${GENO_OUT}/_pool_merge_list_clean.txt"
    : > "$CLEAN_LIST"
    FIRST_CLEAN=""
    for cohort in "${COHORTS[@]}"; do
        bf="${GENO_OUT}/${cohort}_harmonized_r2_${R2_THRESHOLD_H2}"
        clean="${GENO_OUT}/_${cohort}_clean"
        "$PLINK2" \
            --bfile "$bf" \
            --exclude "$MISSNP" \
            --make-bed \
            --threads 4 \
            --out "$clean" \
            2>&1 | tee "${LOG_DIR}/02_clean_${cohort}.log"
        if [[ -z "$FIRST_CLEAN" ]]; then
            FIRST_CLEAN="$clean"
        else
            echo "${clean}.bed ${clean}.bim ${clean}.fam" >> "$CLEAN_LIST"
        fi
    done

    rm -f "${POOLED}".{bed,bim,fam,log,nosex}
    log_info "Merge attempt 2 (after excluding ${n_missnp} missnp variants)..."
    "$PLINK19" \
        --bfile "$FIRST_CLEAN" \
        --merge-list "$CLEAN_LIST" \
        --make-bed \
        --out "$POOLED" \
        2>&1 | tee "${LOG_DIR}/02_merge_attempt2.log"

    # Clean up per-cohort intermediates.
    for cohort in "${COHORTS[@]}"; do
        rm -f "${GENO_OUT}/_${cohort}_clean".{bed,bim,fam,log,nosex}
    done
    rm -f "$CLEAN_LIST"
fi

if [[ ! -f "${POOLED}.bed" ]]; then
    log_error "Pooled merge failed. See ${LOG_DIR}/02_merge_attempt*.log"
    exit 1
fi

n_samples=$(wc -l < "${POOLED}.fam")
n_variants=$(wc -l < "${POOLED}.bim")
log_info ""
log_info "Pooled bfile: ${n_samples} samples × ${n_variants} variants → ${POOLED}.{bed,bim,fam}"
log_info "=== Step 02 complete ==="
