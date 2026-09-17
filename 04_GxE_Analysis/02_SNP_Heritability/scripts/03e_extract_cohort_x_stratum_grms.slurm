#!/usr/bin/env bash
#SBATCH --job-name=h2_blockcomp_subprune
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=20:00
#SBATCH --output=logs/blockcomp_subprune_%x_%a.%j.out
#SBATCH --error=logs/blockcomp_subprune_%x_%a.%j.err
# =============================================================================
# 03e_extract_cohort_x_stratum_grms.slurm — Subset-then-prune, one block
# =============================================================================
# Submitted by 03e_extract_cohort_x_stratum_grms.sh. Each task does two
# gcta calls:
#
#   (1) Subset the pre-cutoff per-cohort GRM to the (cohort × stratum)
#       sub-sample defined by <STRATUM>.iids
#         gcta --grm-bin grm_cohort_<COHORT> --keep <STRATUM>.iids \
#              --make-grm-bin --out <intermediate>
#
#   (2) Apply --grm-cutoff to the intermediate sub-block to remove within-
#       block relatives:
#         gcta --grm-bin <intermediate> --grm-cutoff $GRM_CUTOFF_BLOCK \
#              --make-grm-bin --out <out_prefix>
#
# Together these implement subset-then-prune: cross-stratum within-cohort
# relatives are retained as block-diagonal information; within each block,
# the conventional 0.05 cutoff is applied.
#
# The intermediate GRM is written under TMPDIR or output/grm/ depending on
# size and cleaned up at the end.
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${TASKS_FILE:?TASKS_FILE not exported by the submitter}"
: "${GRM_CUTOFF_BLOCK:?GRM_CUTOFF_BLOCK not exported by the submitter}"
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

INTERMEDIATE="${OUT_PREFIX}.precut"

log_info "================================================================"
log_info "Block component (subset-then-prune): ${COHORT} × ${STRATUM} (task ${IDX})"
log_info "  Source GRM (pre-cutoff): ${IN_PREFIX}"
log_info "  Keep file:               ${KEEP_PATH}"
log_info "  Intermediate prefix:     ${INTERMEDIATE}"
log_info "  Out prefix:              ${OUT_PREFIX}"
log_info "  Within-block --grm-cutoff: ${GRM_CUTOFF_BLOCK}"
log_info "================================================================"

if [[ -f "${OUT_PREFIX}.grm.bin" ]]; then
    N_DONE=$(wc -l < "${OUT_PREFIX}.grm.id")
    log_info "SKIP: block component already built (N=${N_DONE})"
    exit 0
fi

if [[ ! -f "${IN_PREFIX}.grm.bin" ]]; then
    log_error "Source GRM missing: ${IN_PREFIX}.grm.bin"
    exit 1
fi
if [[ ! -f "$KEEP_PATH" ]]; then
    log_error "Keep file missing: $KEEP_PATH"
    exit 1
fi

cleanup_intermediate() {
    rm -f "${INTERMEDIATE}.grm.bin" \
          "${INTERMEDIATE}.grm.id" \
          "${INTERMEDIATE}.grm.N.bin" \
          "${INTERMEDIATE}.log"
}
trap cleanup_intermediate EXIT

# ---- (1) Subset to (cohort × stratum) sub-sample (pre-cutoff) -----------
log_info ""
log_info "Step (1): subset pre-cutoff cohort GRM to ${COHORT} × ${STRATUM}..."
"$GCTA" \
    --grm-bin "$IN_PREFIX" \
    --keep "$KEEP_PATH" \
    --make-grm-bin \
    --out "$INTERMEDIATE" \
    --thread-num "${SLURM_CPUS_PER_TASK:-2}"

if [[ ! -f "${INTERMEDIATE}.grm.bin" ]]; then
    log_error "Intermediate (subset) GRM not written"
    exit 1
fi
N_PRECUT=$(wc -l < "${INTERMEDIATE}.grm.id")
log_info "  pre-cutoff N = ${N_PRECUT}"

# ---- (2) Within-block --grm-cutoff --------------------------------------
log_info ""
log_info "Step (2): apply --grm-cutoff ${GRM_CUTOFF_BLOCK} within (${COHORT} × ${STRATUM})..."
"$GCTA" \
    --grm-bin "$INTERMEDIATE" \
    --grm-cutoff "$GRM_CUTOFF_BLOCK" \
    --make-grm-bin \
    --out "$OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-2}"

if [[ -f "${OUT_PREFIX}.grm.bin" ]]; then
    N_OUT=$(wc -l < "${OUT_PREFIX}.grm.id")
    LOST=$((N_PRECUT - N_OUT))
    log_info "  ${COHORT} × ${STRATUM}: pre=${N_PRECUT}, post=${N_OUT} (lost ${LOST}) → ${OUT_PREFIX}.grm.{bin,id,N.bin}"
else
    log_error "Post-cutoff GRM not written for ${COHORT} × ${STRATUM}"
    exit 1
fi
