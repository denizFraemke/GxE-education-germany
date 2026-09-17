#!/usr/bin/env bash
#SBATCH --job-name=h2_unrel
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=30:00
#SBATCH --output=logs/unrel_%x_%a.%j.out
#SBATCH --error=logs/unrel_%x_%a.%j.err
# =============================================================================
# 04_unrelated_filter.slurm — Per-stratum --grm-cutoff (SLURM array task)
# =============================================================================
# Submitted by 04_unrelated_filter.sh with --array=1-<N>. Each task reads
# the corresponding stratum's GRM from step 03 and writes a kinship-filtered
# GRM via gcta --grm-cutoff at $GRM_CUTOFF (config.sh).
#
# Resource sizing (verified small):
#   cpus-per-task=2 — modest threading; this step is IO + pairwise compare
#   mem-per-cpu=4G  — largest cell GRM (N=8007) is ~256MB; 8G total
#   time=30:00      — even East_pre1990 finishes in a few minutes
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${STRATA_LIST:?STRATA_LIST not exported by the submitter}"
IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"
STRATUM="$(awk -v idx="$IDX" 'NR == idx { print; exit }' "$STRATA_LIST")"
if [[ -z "$STRATUM" ]]; then
    log_error "No stratum at line $IDX of $STRATA_LIST"
    exit 1
fi

GRM_IN="${GRM_OUT}/grm_${STRATUM}"
GRM_OUT_PREFIX="${GRM_OUT}/grm_${STRATUM}_unrel"

log_info "================================================================"
log_info "Unrelated filter for stratum: ${STRATUM} (array task ${IDX})"
log_info "  Input GRM:  ${GRM_IN}"
log_info "  Output GRM: ${GRM_OUT_PREFIX}"
log_info "  Cutoff:     ${GRM_CUTOFF}"
log_info "================================================================"

if [[ -f "${GRM_OUT_PREFIX}.grm.bin" ]]; then
    N_AFTER=$(wc -l < "${GRM_OUT_PREFIX}.grm.id")
    log_info "SKIP: already filtered (N=${N_AFTER})"
    exit 0
fi

if [[ ! -f "${GRM_IN}.grm.bin" ]]; then
    log_error "Input GRM missing: ${GRM_IN}.grm.bin"
    log_error "Run step 03 first."
    exit 1
fi

N_BEFORE=$(wc -l < "${GRM_IN}.grm.id")

"$GCTA" \
    --grm-bin "$GRM_IN" \
    --grm-cutoff "$GRM_CUTOFF" \
    --make-grm-bin \
    --out "$GRM_OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-2}"

if [[ -f "${GRM_OUT_PREFIX}.grm.bin" ]]; then
    N_AFTER=$(wc -l < "${GRM_OUT_PREFIX}.grm.id")
    log_info "  ${STRATUM}: ${N_BEFORE} -> ${N_AFTER}"
else
    log_error "Output GRM not written for ${STRATUM}"
    exit 1
fi
