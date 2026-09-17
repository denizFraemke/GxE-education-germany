#!/usr/bin/env bash
# =============================================================================
# 00_check_inputs.sh — Verify all required input files exist before running
# =============================================================================
# USAGE:
#   bash scripts/00_check_inputs.sh
#   (run via `bash run_pipeline.sh --step 00`)
#
# INPUTS (existence + readability only):
#   ${GCTB}, ${PLINK19}, ${PLINK2}                       — tool binaries
#   ${LD_REF_DIR}/ukb50k_shrunk_chr{1..22}_mafpt01.ldm.sparse[.info]
#                                                        — UKB 50k LD matrices
#   ${SUMSTATS_DIR}/<trait subdir>/<per-chr file>        — per-trait sumstats
#   $(get_geno_bfile <cohort>).{bed,bim,fam}             — per-cohort bfiles
#
# OUTPUTS:
#   None — exits 0 if all checks pass, non-zero with a log line per missing
#   file otherwise.
#
# DEPENDS ON: nothing (first step of the DAG).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"

log_info "=== Step 00: Checking all required input files ==="

ERRORS=0

# --- Check tool binaries ----------------------------------------------------
for tool_var in GCTB PLINK19 PLINK2; do
    tool_path="${!tool_var}"
    if [[ ! -x "$tool_path" ]]; then
        log_error "Tool not found or not executable: $tool_var = $tool_path"
        ERRORS=$((ERRORS + 1))
    else
        log_info "OK: $tool_var = $tool_path"
    fi
done

# --- Check LD reference matrices --------------------------------------------
for CHR in $(seq 1 22); do
    ldm="${LD_REF_DIR}/ukb50k_shrunk_chr${CHR}_mafpt01.ldm.sparse"
    if [[ ! -f "${ldm}.info" ]] && [[ ! -f "$ldm" ]]; then
        log_error "LD matrix not found for chr${CHR}: $ldm"
        ERRORS=$((ERRORS + 1))
    fi
done
if [[ $ERRORS -eq 0 ]]; then
    log_info "OK: All 22 LD reference matrices found"
fi

# --- Check formatted summary statistics per trait ----------------------------
for trait in "${TRAIT_NAMES[@]}"; do
    subdir="${SUMSTATS_DIR}/${TRAIT_SUMSTATS_SUBDIR[$trait]}"
    pattern="${TRAIT_CHR_FILE_PATTERN[$trait]}"
    missing_chr=0
    for CHR in $(seq 1 22); do
        fname="${pattern//\{CHR\}/$CHR}"
        if [[ ! -f "${subdir}/${fname}" ]]; then
            missing_chr=$((missing_chr + 1))
        fi
    done
    if [[ $missing_chr -gt 0 ]]; then
        log_error "Trait $trait: $missing_chr of 22 per-chromosome sumstat files missing in $subdir"
        ERRORS=$((ERRORS + 1))
    else
        # Sanity-check the .ma column count. SBayesR expects exactly 8 columns
        # (SNP A1 A2 freq b se p N) and reads them positionally — a malformed
        # file (e.g. an extra leading CHR column) shifts every column, giving
        # "0 matched SNPs" and, on an NA row, a GCTB segfault at step 01.
        chr1_file="${subdir}/${pattern//\{CHR\}/1}"
        ncol=$(awk 'NR==1{print NF; exit}' "$chr1_file" 2>/dev/null || echo 0)
        if [[ "$ncol" != "8" ]]; then
            log_error "Trait $trait: $chr1_file has $ncol columns, expected 8 (SNP A1 A2 freq b se p N)."
            if [[ "$trait" == "Height" ]]; then
                log_error "  Height's pre-formatted share copy is malformed — regenerate from raw GIANT:"
                log_error "    bash scripts/format_height_sumstats.sh   (set GIANT_HEIGHT_RAW in user_config.sh)"
            fi
            ERRORS=$((ERRORS + 1))
        else
            log_info "OK: Trait $trait — all 22 chromosome files found in $subdir (8-col .ma)"
        fi
    fi
done

# --- Check genotype data (bed/bim/fam) per cohort ---------------------------
for cohort in "${COHORTS[@]}"; do
    bfile=$(get_geno_bfile "$cohort")
    for ext in .bed .bim .fam; do
        if [[ ! -f "${bfile}${ext}" ]]; then
            log_error "Genotype file not found: ${bfile}${ext}"
            ERRORS=$((ERRORS + 1))
        fi
    done
    if [[ -f "${bfile}.bed" ]]; then
        n_samples=$(wc -l < "${bfile}.fam" 2>/dev/null || echo 0)
        n_variants=$(wc -l < "${bfile}.bim" 2>/dev/null || echo 0)
        log_info "OK: $cohort genotype — $n_samples samples, $n_variants variants"
    fi
done

# --- Summary -----------------------------------------------------------------
if [[ $ERRORS -gt 0 ]]; then
    log_error "=== INPUT CHECK FAILED: $ERRORS errors found ==="
    log_error "Fix the issues above before running the pipeline."
    log_error "Missing sumstats or genotypes? Run:  bash scripts/setup_tardis_data.sh"
    log_error "  It stages them onto the cluster from the data share. The formatted"
    log_error "  per-chromosome sumstats live at Projects/04_data_analysis/013-PGS/data/formatted/"
    log_error "  (subdir per trait); see that script's header for every source path."
    exit 1
else
    log_info "=== All input checks passed ==="
fi
