#!/usr/bin/env bash
#SBATCH --job-name=h2_grm
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=8:00:00
#SBATCH --output=logs/grm_%x_%a.%j.out
#SBATCH --error=logs/grm_%x_%a.%j.err
# =============================================================================
# 03_build_stratum_grms.slurm — One SLURM array task = one stratum's GRM
# =============================================================================
# Mirrors 01_Genotype/Harmonized/scripts/01_run_sbayesr.slurm.
#
# Submitted by 03_build_stratum_grms.sh with --array=1-<N>%2 (throttled to
# at most 2 concurrent tasks across the array). Each task runs one
# gcta64 --make-grm-bin call against the pooled bfile, masked to that
# stratum's individuals via --keep.
#
# Resource sizing rationale:
#   partition=short  — 1d time limit (short has the idle ood-* nodes).
#
#   --exclusive      — each task claims the entire node. Several gcta64
#                      processes per node, all mmap-reading the same 21 GB
#                      beegfs-hosted pooled.bed simultaneously, block in
#                      fuse_d kernel waits; FUSE/beegfs cannot serve that
#                      many concurrent readers from one mount, and
#                      whole-node isolation is the simplest fix.
#
#   mem=200G         — explicit absolute ceiling. --mem=0 is not honoured
#                      as "no cap" for this account/partition. 200G sits
#                      well above GCTA's peak (~115G for the largest cell,
#                      East_pre1990 at N=8007) and below the node's 256G
#                      capacity, leaving OS / page-cache headroom.
#
#   cpus-per-task=16 — matrix math is the bottleneck (N^2 × n_SNPs
#                      operations). With --exclusive the whole node's CPUs
#                      are reserved anyway, so 16 of them are passed to
#                      gcta64 --thread-num.
#
#   time=8:00:00     — East_pre1990 at N=8007 runs roughly 30 min matrix
#                      math + 30 min IO/recoding at 16 threads; the 8h
#                      ceiling covers any pathological beegfs slowdown.
#
# Why this is a SLURM job rather than a login-node loop: login-node cgroup
# CPU throttling stalls gcta64 (the same reason 01_run_sbayesr.slurm exists
# in the upstream genotype pipeline). Compute nodes get full cores.
#
# USAGE (submitted automatically):
#   sbatch --array=1-<N> --wait \
#          --export=ALL,STRATA_LIST=<path> \
#          scripts/03_build_stratum_grms.slurm
#
# INPUTS (passed via --export):
#   STRATA_LIST — path to a file with one stratum name per line; the array
#                 task ID indexes into it (1-based).
#
# OUTPUTS:
#   ${GRM_OUT}/grm_<stratum>.grm.{bin,id,N.bin}
# =============================================================================

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

: "${STRATA_LIST:?STRATA_LIST not exported by the submitter}"
if [[ ! -f "$STRATA_LIST" ]]; then
    log_error "STRATA_LIST not found: $STRATA_LIST"
    exit 1
fi

# Resolve this task's stratum from the list (array task IDs are 1-based).
IDX="${SLURM_ARRAY_TASK_ID:?must be running inside a SLURM array}"
STRATUM="$(awk -v idx="$IDX" 'NR == idx { print; exit }' "$STRATA_LIST")"
if [[ -z "$STRATUM" ]]; then
    log_error "No stratum at line $IDX of $STRATA_LIST"
    exit 1
fi

GCTA_KEEP="${PROJ_ROOT}/data/keep/${STRATUM}.gcta.iids"
GRM_PREFIX="${GRM_OUT}/grm_${STRATUM}"

log_info "================================================================"
log_info "Building GRM for stratum: ${STRATUM} (array task ${IDX})"
log_info "  Pooled bfile: ${GENO_OUT}/pooled"
log_info "  Keep list:    ${GCTA_KEEP}"
log_info "  Out prefix:   ${GRM_PREFIX}"
log_info "================================================================"

# Resume mode — skip if the GRM is already on disk.
if [[ -f "${GRM_PREFIX}.grm.bin" ]]; then
    N_DONE=$(wc -l < "${GRM_PREFIX}.grm.id")
    log_info "SKIP: ${STRATUM} GRM already built (N=${N_DONE})"
    exit 0
fi

if [[ ! -f "$GCTA_KEEP" ]]; then
    log_error "Resolved keep file missing: $GCTA_KEEP"
    log_error "The submitter (03_build_stratum_grms.sh) should have produced this."
    exit 1
fi
N_KEEP=$(wc -l < "$GCTA_KEEP")
log_info "  N in keep file: ${N_KEEP}"

# NOTE: GCTA does not accept --hwe (PLINK does). See Plan_deviations.md §7d
# for the rationale.
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
    log_error "GRM output missing for ${STRATUM}"
    exit 1
fi

log_info "Step 03 task for ${STRATUM} complete."
