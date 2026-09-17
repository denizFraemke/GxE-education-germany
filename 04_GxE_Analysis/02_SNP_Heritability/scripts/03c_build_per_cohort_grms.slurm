#!/usr/bin/env bash
#SBATCH --job-name=h2_grm_cohort
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=8:00:00
#SBATCH --output=logs/grm_cohort_%x_%a.%j.out
#SBATCH --error=logs/grm_cohort_%x_%a.%j.err
# =============================================================================
# 03c_build_per_cohort_grms.slurm — One SLURM array task = one cohort's GRM
# =============================================================================
# Submitted by 03c_build_per_cohort_grms.sh with --array=1-<N>%2 (same
# beegfs/FUSE concurrency cap as step 03). Each task runs gcta64
# --make-grm-bin against pooled.bed, masked to one cohort via --keep.
#
# Resource sizing is identical to step 03 — SHIP has ~7,770 individuals,
# similar to the East × pre-1990 cell. See 03_build_stratum_grms.slurm for
# the full rationale (exclusive, 200G, 16 threads, 8h ceiling).
# =============================================================================
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${COHORTS_LIST:?COHORTS_LIST not exported by the submitter}"
if [[ ! -f "$COHORTS_LIST" ]]; then
    log_error "COHORTS_LIST not found: $COHORTS_LIST"
    exit 1
fi

IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"
COHORT="$(awk -v idx="$IDX" 'NR == idx { print; exit }' "$COHORTS_LIST")"
if [[ -z "$COHORT" ]]; then
    log_error "No cohort at line $IDX of $COHORTS_LIST"
    exit 1
fi

GCTA_KEEP="${PROJ_ROOT}/data/keep/cohort_${COHORT}.gcta.iids"
GRM_PREFIX="${GRM_OUT}/grm_cohort_${COHORT}"

log_info "================================================================"
log_info "Building per-cohort GRM: ${COHORT} (array task ${IDX})"
log_info "  Pooled bfile: ${GENO_OUT}/pooled"
log_info "  Keep list:    ${GCTA_KEEP}"
log_info "  Out prefix:   ${GRM_PREFIX}"
log_info "================================================================"

if [[ -f "${GRM_PREFIX}.grm.bin" ]]; then
    N_DONE=$(wc -l < "${GRM_PREFIX}.grm.id")
    log_info "SKIP: ${COHORT} GRM already built (N=${N_DONE})"
    exit 0
fi

if [[ ! -f "$GCTA_KEEP" ]]; then
    log_error "Resolved keep file missing: $GCTA_KEEP"
    log_error "The submitter (03c_build_per_cohort_grms.sh) should have produced this."
    exit 1
fi
N_KEEP=$(wc -l < "$GCTA_KEEP")
log_info "  N in keep file: ${N_KEEP}"

"$GCTA" \
    --bfile "${GENO_OUT}/pooled" \
    --keep "$GCTA_KEEP" \
    --autosome \
    --maf "$GREML_MAF" \
    --make-grm-bin \
    --out "$GRM_PREFIX" \
    --thread-num "${SLURM_CPUS_PER_TASK:-4}"

if [[ -f "${GRM_PREFIX}.grm.bin" ]]; then
    N_GRM=$(wc -l < "${GRM_PREFIX}.grm.id")
    log_info "  GRM N=${N_GRM} → ${GRM_PREFIX}.grm.{bin,id,N.bin}"
else
    log_error "GRM output missing for ${COHORT}"
    exit 1
fi

log_info "Step 03c task for ${COHORT} complete."
