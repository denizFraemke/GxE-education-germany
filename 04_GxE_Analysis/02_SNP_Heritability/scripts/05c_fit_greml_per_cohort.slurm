#!/usr/bin/env bash
#SBATCH --job-name=h2_greml_cohort
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/reml_cohort_%x_%a.%j.out
#SBATCH --error=logs/reml_cohort_%x_%a.%j.err
# =============================================================================
# 05c_fit_greml_per_cohort.slurm — One per-cohort GREML fit per array task
# =============================================================================
# Each task reads one line of FITS_LIST and runs gcta --reml against a
# per-cohort unrel GRM, optionally with --keep to subset to a (cohort ×
# region) or (cohort × stratum) sub-sample.
#
# FITS_LIST line format (tab-separated):
#   <label>\t<grm_path>\t<keep_path_or_NONE>\t<extra_flags>
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
KEEP_PATH=$(echo "$ROW" | cut -f3)
EXTRA=$(   echo "$ROW" | cut -f4)

OUT_PREFIX="${OUT_DIR}/${LABEL}"

log_info "================================================================"
log_info "Per-cohort GREML fit: ${LABEL} (array task ${IDX})"
log_info "  GRM source:  ${GRM_PATH}"
log_info "  Keep filter: ${KEEP_PATH}"
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

# Build optional --keep flag.
KEEP_FLAG=()
if [[ "$KEEP_PATH" != "NONE" ]]; then
    if [[ ! -f "$KEEP_PATH" ]]; then
        log_error "Keep file missing: $KEEP_PATH"
        exit 1
    fi
    KEEP_FLAG=(--keep "$KEEP_PATH")
fi

# shellcheck disable=SC2086
"$GCTA" \
    --reml \
    --grm-bin "$GRM_PATH" \
    "${KEEP_FLAG[@]}" \
    --pheno "$PHENO_FILE" \
    --qcovar "$QCOVAR_FILE" \
    --covar "$COVAR_FILE" \
    --out "$OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-4}" \
    $EXTRA

if [[ -f "${OUT_PREFIX}.hsq" ]]; then
    log_info "  → ${OUT_PREFIX}.hsq"
else
    log_warn "  No .hsq produced (REML did not converge or sample too small)"
    # Don't fail the array task — downstream parsing handles missing .hsq.
fi
