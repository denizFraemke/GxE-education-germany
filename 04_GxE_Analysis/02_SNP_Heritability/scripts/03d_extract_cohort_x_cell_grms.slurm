#!/usr/bin/env bash
#SBATCH --job-name=h2_blockcomp
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=15:00
#SBATCH --output=logs/blockcomp_%x_%a.%j.out
#SBATCH --error=logs/blockcomp_%x_%a.%j.err
# =============================================================================
# 03d_extract_cohort_x_cell_grms.slurm — Extract one block component
# =============================================================================
# Submitted by 03d_extract_cohort_x_cell_grms.sh. Each task reads one line
# of TASKS_FILE and runs gcta --grm-bin --keep --make-grm-bin to subset a
# per-cohort GRM to a (cohort × cell) sub-sample. Output is a "block
# component" GRM that step 04d combines block-diagonally into the cell GRM.
#
# Resource sizing: GRM subsetting is cheap (no genotype I/O), but the
# largest input GRM (SHIP, N≈6,673) is ~85 MB packed lower-triangular;
# allow 2 CPUs × 4 GB and 15 min generous ceiling.
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${TASKS_FILE:?TASKS_FILE not exported by the submitter}"
IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"

ROW="$(awk -v idx="$IDX" -F '\t' 'NR == idx { print; exit }' "$TASKS_FILE")"
if [[ -z "$ROW" ]]; then
    log_error "No task at line $IDX of $TASKS_FILE"
    exit 1
fi

COHORT=$(    echo "$ROW" | cut -f1)
STRATUM=$(   echo "$ROW" | cut -f2)
IN_PREFIX=$( echo "$ROW" | cut -f3)
KEEP_PATH=$( echo "$ROW" | cut -f4)
OUT_PREFIX=$(echo "$ROW" | cut -f5)

log_info "================================================================"
log_info "Block component: ${COHORT} × ${STRATUM} (array task ${IDX})"
log_info "  Input GRM:  ${IN_PREFIX}"
log_info "  Keep file:  ${KEEP_PATH}"
log_info "  Out prefix: ${OUT_PREFIX}"
log_info "================================================================"

if [[ -f "${OUT_PREFIX}.grm.bin" ]]; then
    N_DONE=$(wc -l < "${OUT_PREFIX}.grm.id")
    log_info "SKIP: block component already built (N=${N_DONE})"
    exit 0
fi

if [[ ! -f "${IN_PREFIX}.grm.bin" ]]; then
    log_error "Input GRM missing: ${IN_PREFIX}.grm.bin"
    exit 1
fi
if [[ ! -f "$KEEP_PATH" ]]; then
    log_error "Keep file missing: $KEEP_PATH"
    exit 1
fi

"$GCTA" \
    --grm-bin "$IN_PREFIX" \
    --keep "$KEEP_PATH" \
    --make-grm-bin \
    --out "$OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-2}"

if [[ -f "${OUT_PREFIX}.grm.bin" ]]; then
    N_OUT=$(wc -l < "${OUT_PREFIX}.grm.id")
    log_info "  ${COHORT} × ${STRATUM}: N=${N_OUT} → ${OUT_PREFIX}.grm.{bin,id,N.bin}"
else
    log_error "Block component not written for ${COHORT} × ${STRATUM}"
    exit 1
fi
