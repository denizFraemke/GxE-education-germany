#!/usr/bin/env bash
#SBATCH --job-name=h2_greml_main
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/reml_%x_%a.%j.out
#SBATCH --error=logs/reml_%x_%a.%j.err
# =============================================================================
# 05_fit_greml_main.slurm — One GREML fit per SLURM array task (4-cell main)
# =============================================================================
# Each task reads one line of FITS_LIST and runs gcta --reml. Three fit
# flavours: joint mGRM, single-GRM standalone, both ± --reml-no-constrain.
#
# FITS_LIST line format (tab-separated):
#   <fit_label>\t<fit_type>\t<grm_or_mgrm_path>\t<extra_flags>
# where fit_type ∈ {mgrm, grm}.
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

# Parse the tab-separated row.
LABEL=$(   echo "$ROW" | cut -f1)
FIT_TYPE=$(echo "$ROW" | cut -f2)
GRM_PATH=$(echo "$ROW" | cut -f3)
EXTRA=$(   echo "$ROW" | cut -f4)

OUT_PREFIX="${OUT_DIR}/${LABEL}"

log_info "================================================================"
log_info "GREML fit: ${LABEL} (array task ${IDX})"
log_info "  Fit type:    ${FIT_TYPE}"
log_info "  GRM source:  ${GRM_PATH}"
log_info "  Extra flags: ${EXTRA:-<none>}"
log_info "  Out prefix:  ${OUT_PREFIX}"
log_info "================================================================"

if [[ -f "${OUT_PREFIX}.hsq" ]]; then
    log_info "SKIP: ${LABEL}.hsq already exists"
    exit 0
fi

# Build the gcta invocation.
case "$FIT_TYPE" in
    mgrm)
        if [[ ! -f "$GRM_PATH" ]]; then
            log_error "mGRM list file missing: $GRM_PATH"
            exit 1
        fi
        GCTA_GRM_FLAG=(--mgrm "$GRM_PATH")
        ;;
    grm)
        if [[ ! -f "${GRM_PATH}.grm.bin" ]]; then
            log_error "GRM missing: ${GRM_PATH}.grm.bin"
            exit 1
        fi
        GCTA_GRM_FLAG=(--grm-bin "$GRM_PATH")
        ;;
    *)
        log_error "Unknown fit type: $FIT_TYPE (expected mgrm or grm)"
        exit 1
        ;;
esac

# shellcheck disable=SC2086
"$GCTA" \
    --reml \
    "${GCTA_GRM_FLAG[@]}" \
    --pheno "$PHENO_FILE" \
    --qcovar "$QCOVAR_FILE" \
    --covar "$COVAR_FILE" \
    --out "$OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-4}" \
    $EXTRA

if [[ -f "${OUT_PREFIX}.hsq" ]]; then
    log_info "  → ${OUT_PREFIX}.hsq"
else
    log_warn "  No .hsq produced (REML did not converge or input invalid)"
    # Don't fail the array task — downstream parsing handles missing .hsq
    # gracefully (07_lrt.R, 09_assemble_output.R both check presence).
fi
