#!/usr/bin/env bash
#SBATCH --job-name=h2_unrel_cohort
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=30:00
#SBATCH --output=logs/unrel_cohort_%x_%a.%j.out
#SBATCH --error=logs/unrel_cohort_%x_%a.%j.err
# =============================================================================
# 04c_unrelated_filter_per_cohort.slurm — Per-cohort --grm-cutoff (SLURM array)
# =============================================================================
# Submitted by 04c_unrelated_filter_per_cohort.sh. Each task reads one
# cohort's GRM from step 03c and writes a within-cohort unrelated-set GRM
# via gcta --grm-cutoff. Cutoff = GRM_CUTOFF_PER_COHORT (0.05, published
# convention; no cross-cohort fingerprint to absorb inside a single cohort).
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${COHORTS_LIST:?COHORTS_LIST not exported by the submitter}"
IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"
COHORT="$(awk -v idx="$IDX" 'NR == idx { print; exit }' "$COHORTS_LIST")"
if [[ -z "$COHORT" ]]; then
    log_error "No cohort at line $IDX of $COHORTS_LIST"
    exit 1
fi

GRM_IN="${GRM_OUT}/grm_cohort_${COHORT}"
GRM_OUT_PREFIX="${GRM_OUT}/grm_cohort_${COHORT}_unrel"

log_info "================================================================"
log_info "Unrelated filter for cohort: ${COHORT} (array task ${IDX})"
log_info "  Input GRM:  ${GRM_IN}"
log_info "  Output GRM: ${GRM_OUT_PREFIX}"
log_info "  Cutoff:     ${GRM_CUTOFF_PER_COHORT}"
log_info "================================================================"

if [[ -f "${GRM_OUT_PREFIX}.grm.bin" ]]; then
    N_AFTER=$(wc -l < "${GRM_OUT_PREFIX}.grm.id")
    log_info "SKIP: already filtered (N=${N_AFTER})"
    exit 0
fi

if [[ ! -f "${GRM_IN}.grm.bin" ]]; then
    log_error "Input GRM missing: ${GRM_IN}.grm.bin"
    log_error "Run step 03c first."
    exit 1
fi

N_BEFORE=$(wc -l < "${GRM_IN}.grm.id")

"$GCTA" \
    --grm-bin "$GRM_IN" \
    --grm-cutoff "$GRM_CUTOFF_PER_COHORT" \
    --make-grm-bin \
    --out "$GRM_OUT_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-2}"

if [[ -f "${GRM_OUT_PREFIX}.grm.bin" ]]; then
    N_AFTER=$(wc -l < "${GRM_OUT_PREFIX}.grm.id")
    log_info "  ${COHORT}: ${N_BEFORE} -> ${N_AFTER}"
else
    log_error "Output GRM not written for ${COHORT}"
    exit 1
fi
