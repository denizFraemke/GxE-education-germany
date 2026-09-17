#!/usr/bin/env bash
# =============================================================================
# 01_rebuild_bfiles_info090.sh — Per-cohort INFO>=.90 harmonized bfiles
# =============================================================================
# Plan S5 stipulates harmonized variants at imputation INFO >= 0.90 (stricter
# than the PGI pipeline's 0.30). We reuse the upstream PGI pipeline's
# pass-list machinery rather than re-implementing it:
#
#   1. Invoke the upstream 03a_filter_imputation_quality.sh from the genotype
#      project with R2_THRESHOLD_STRICT=0.9. That step is idempotent — it
#      writes ${GENOTYPE_OUT_GENO}/{COHORT}_r2_pass_0.9.txt if not present.
#
#   2. Compute the cross-cohort intersection of those pass lists (this
#      mirrors the SHARED_SNP_SET phase 2 logic in upstream 03 but only for
#      the strict 0.9 list; the upstream step does this already and writes
#      ${GENOTYPE_OUT_GENO}/_shared_snp_set_r2_0.9.txt — we just reuse it
#      if present; otherwise rebuild it here).
#
#   3. For each cohort, plink2 --extract that shared list from the cohort's
#      Phase-1 {COHORT}_harmonized.{bed,bim,fam} (NOT from the lenient
#      0.3-restricted {COHORT}_harmonized_shared.* — we want a separate,
#      strict set for h²).
#
# Outputs land in ${GENO_OUT}/{COHORT}_harmonized_r2_0.9.{bed,bim,fam}.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 01: Rebuild harmonized bfiles at R² >= ${R2_THRESHOLD_H2} ==="

UPSTREAM_03A="${GENOTYPE_PROJ_ROOT}/scripts/03a_filter_imputation_quality.sh"
if [[ ! -x "$UPSTREAM_03A" ]] && [[ ! -f "$UPSTREAM_03A" ]]; then
    log_error "Upstream R²-filter script missing: $UPSTREAM_03A"
    log_error "Make sure the genotype pipeline is staged at GENOTYPE_PROJ_ROOT."
    exit 1
fi

# Force the strict pass-list at 0.9 by overriding the env var the upstream
# script honours. Lenient threshold left at upstream default.
log_info "Invoking upstream 03a with R2_THRESHOLD_STRICT=${R2_THRESHOLD_H2}..."
R2_FILTER=yes R2_THRESHOLD_STRICT="${R2_THRESHOLD_H2}" \
    bash "$UPSTREAM_03A" 2>&1 | tee "${LOG_DIR}/01a_upstream_03a.log"

# --- Compute (or reuse) cross-cohort intersection of the 0.9 pass lists ----
SHARED_R2_LIST="${GENOTYPE_OUT_GENO}/_shared_snp_set_r2_${R2_THRESHOLD_H2}.txt"

if [[ -s "$SHARED_R2_LIST" ]]; then
    log_info "Reusing existing shared list: $SHARED_R2_LIST ($(wc -l < "$SHARED_R2_LIST") SNPs)"
else
    log_info "Computing cross-cohort intersection of *_r2_pass_${R2_THRESHOLD_H2}.txt..."
    tmp="${GENO_OUT}/_intersect_tmp.txt"
    first="${COHORTS[0]}"
    sort "${GENOTYPE_OUT_GENO}/${first}_r2_pass_${R2_THRESHOLD_H2}.txt" > "$tmp"
    log_info "  ${first}: $(wc -l < "$tmp") SNPs"
    for c in "${COHORTS[@]:1}"; do
        pass="${GENOTYPE_OUT_GENO}/${c}_r2_pass_${R2_THRESHOLD_H2}.txt"
        if [[ ! -s "$pass" ]]; then
            log_error "  Missing or empty: $pass — re-run upstream 03a."
            exit 1
        fi
        sort "$pass" | comm -12 - "$tmp" > "${tmp}.next"
        mv "${tmp}.next" "$tmp"
        log_info "  ${c}: $(wc -l < "$tmp") SNPs remaining"
    done
    mv "$tmp" "$SHARED_R2_LIST"
    log_info "Wrote: $SHARED_R2_LIST"
fi

N_SHARED=$(wc -l < "$SHARED_R2_LIST")
log_info "Cross-cohort INFO>=${R2_THRESHOLD_H2} shared SNP set: ${N_SHARED} SNPs"

if (( N_SHARED < 100000 )); then
    log_warn "Shared SNP set has only ${N_SHARED} variants (< 100k). GREML should"
    log_warn "still run but variance estimation may be less stable."
fi

# --- Per-cohort --extract from the Phase-1 (unrestricted) bfile -------------
for cohort in "${COHORTS[@]}"; do
    phase1_bfile="${GENOTYPE_OUT_GENO}/${cohort}_harmonized"
    out_prefix="${GENO_OUT}/${cohort}_harmonized_r2_${R2_THRESHOLD_H2}"

    if [[ -f "${out_prefix}.bed" ]]; then
        n=$(wc -l < "${out_prefix}.bim")
        log_info "SKIP: $cohort already restricted (${n} variants)"
        continue
    fi

    if [[ ! -f "${phase1_bfile}.bed" ]]; then
        log_error "Phase-1 bfile missing for $cohort: ${phase1_bfile}.bed"
        log_error "Re-run upstream step 03 (Phase 1) in the genotype pipeline."
        exit 1
    fi

    log_info ""
    log_info "Restricting $cohort to INFO>=${R2_THRESHOLD_H2} shared SNPs..."
    "$PLINK2" \
        --bfile "$phase1_bfile" \
        --extract "$SHARED_R2_LIST" \
        --make-bed \
        --threads 4 \
        --out "$out_prefix" \
        2>&1 | tee "${LOG_DIR}/01_extract_${cohort}.log"

    n_after=$(wc -l < "${out_prefix}.bim")
    log_info "  ${cohort}: ${n_after} variants in ${out_prefix}.{bed,bim,fam}"
done

log_info ""
log_info "=== Step 01 complete ==="
