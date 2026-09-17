#!/usr/bin/env bash
#SBATCH --job-name=h2_greml_blockdiag
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/reml_blockdiag_%x_%a.%j.out
#SBATCH --error=logs/reml_blockdiag_%x_%a.%j.err
# =============================================================================
# 05d_fit_greml_block_diagonal.slurm — One block-diagonal GREML fit per task
# =============================================================================
# Submitted by 05d_fit_greml_block_diagonal.sh. Each task reads one line of
# FITS_LIST and runs gcta --reml on a per-cell block-diagonal GRM.
#
# Covar: ${BLOCKDIAG_COVAR_FILE_OVERRIDE:-$BLOCKDIAG_COVAR_FILE}.
#   Cell primary + R1 → gender_cohort.covar (gender + cohort indicator)
#   G1              → region_cohort.covar  (region + cohort indicator)
#   RG1             → cohort_only.covar    (cohort indicator only)
# The cohort indicator absorbs between-cohort phenotype-mean differences
# that would otherwise inflate V_e under the block-diagonal design. See
# Plan_deviations.md §7g (cohort FE rationale) and §7h (new tracks).
#
# FITS_LIST line format (tab-separated):
#   <label>\t<grm_prefix>\t<extra_flags>
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${FITS_LIST:?FITS_LIST not exported by the submitter}"
: "${OUT_DIR:?OUT_DIR not exported by the submitter}"
IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"

ROW="$(awk -v idx="$IDX" -F '\t' 'NR == idx { print; exit }' "$FITS_LIST")"
if [[ -z "$ROW" ]]; then
    log_error "No fit at line $IDX of $FITS_LIST"
    exit 1
fi

LABEL=$(   echo "$ROW" | cut -f1)
GRM_PATH=$(echo "$ROW" | cut -f2)
EXTRA=$(   echo "$ROW" | cut -f3)

OUT_PREFIX="${OUT_DIR}/${LABEL}"

log_info "================================================================"
log_info "Block-diagonal GREML fit: ${LABEL} (array task ${IDX})"
log_info "  GRM source:  ${GRM_PATH}"
log_info "  Extra flags: ${EXTRA:-<none>}"
log_info "  Out prefix:  ${OUT_PREFIX}"
log_info "================================================================"

if [[ -f "${OUT_PREFIX}.hsq" ]]; then
    log_info "SKIP: ${LABEL}.hsq already exists"
    exit 0
fi

if [[ ! -f "${GRM_PATH}.grm.bin" ]]; then
    log_error "GRM missing: ${GRM_PATH}.grm.bin"
    exit 1
fi

__COVAR_FILE="${BLOCKDIAG_COVAR_FILE_OVERRIDE:-$BLOCKDIAG_COVAR_FILE}"
log_info "  Covar file:  ${__COVAR_FILE}"

# shellcheck disable=SC2086
"$GCTA" \
    --reml \
    --grm-bin "$GRM_PATH" \
    --pheno "$PHENO_FILE" \
    --qcovar "$QCOVAR_FILE" \
    --covar "$__COVAR_FILE" \
    --out "$OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-4}" \
    $EXTRA

if [[ -f "${OUT_PREFIX}.hsq" ]]; then
    log_info "  → ${OUT_PREFIX}.hsq"
else
    log_warn "  No .hsq produced (REML did not converge or input invalid)"
fi
