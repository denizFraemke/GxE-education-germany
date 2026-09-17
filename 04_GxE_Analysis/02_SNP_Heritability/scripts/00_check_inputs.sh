#!/usr/bin/env bash
# =============================================================================
# 00_check_inputs.sh — Verify all required inputs are present (Tardis side)
# =============================================================================
# Fails loud if any prerequisite is missing. Run this once after a fresh
# upload to confirm the Tardis environment is wired up before kicking off
# the long-running steps.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 00: Checking inputs ==="

missing=0

check_file() {
    if [[ -e "$1" ]]; then
        log_info "  OK    $1"
    else
        log_error "  MISS  $1  ($2)"
        missing=$((missing + 1))
    fi
}

# --- Binaries ---------------------------------------------------------------
log_info ""
log_info "Tool binaries:"
for tool in "$GCTA" "$PLINK19" "$PLINK2"; do
    if [[ -x "$tool" ]]; then
        log_info "  OK    $tool"
    else
        log_error "  MISS  $tool  (not executable or not present)"
        missing=$((missing + 1))
    fi
done

# --- Upstream genotype outputs ---------------------------------------------
log_info ""
log_info "Upstream Phase-1 harmonized bfiles (from genotype pipeline):"
for cohort in "${COHORTS[@]}"; do
    bfile="${GENOTYPE_OUT_GENO}/${cohort}_harmonized"
    for ext in bed bim fam; do
        check_file "${bfile}.${ext}" "${cohort} Phase-1 ${ext}"
    done
done

log_info ""
log_info "Upstream imputation info files:"
check_file "${GENOTYPE_INFO_DIR}" "info-file directory"

# --- Phenotype + covariate + keep files (from Mac side) --------------------
log_info ""
log_info "Phenotype / covariates / strata (uploaded from Mac):"
check_file "${PHENO_FILE}"   "education.phen"
check_file "${QCOVAR_FILE}"  "crosspcs.qcovar"
check_file "${COVAR_FILE}"   "gender.covar"
check_file "${PROJ_ROOT}/data/manifest.tsv" "phenotype manifest"

log_info ""
log_info "Per-stratum keep files:"
for stratum in "${STRATA_MAIN[@]}"; do
    check_file "$(get_keep_file "$stratum")" "${stratum} keep list"
done

# --- Result ----------------------------------------------------------------
log_info ""
if [[ $missing -eq 0 ]]; then
    log_info "=== Step 00 OK — all inputs present ==="
else
    log_error "=== Step 00 FAILED — ${missing} missing input(s) ==="
    log_error "Resolve missing items above and re-run."
    exit 1
fi
