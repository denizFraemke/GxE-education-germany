#!/usr/bin/env bash
# =============================================================================
# 03d_extract_cohort_x_cell_grms.sh — Extract per-(cohort × cell) GRM blocks
# =============================================================================
# First step of the block-diagonal cell-level track (Plan_deviations.md §7e).
# For each (cohort, cell) sub-sample with N >= MIN_BLOCKCOMP_N, run
#   gcta --grm-bin grm_cohort_<COHORT>_unrel --keep <cohort×cell>.iids \
#        --make-grm-bin --out grm_blockcomp_<COHORT>_<CELL>
# This is a cheap subsetting operation (no genotype I/O — GCTA reads the
# pre-built cohort GRM and writes a sub-GRM). Each task takes seconds.
#
# Reuses the per-cohort unrel GRMs from step 04c AND the per-(cohort × cell)
# keep files from step 05c (they're produced as a side-effect of 05c's
# subset-keep generation, independent of the MIN_SUBCELL_N gate on REML
# fits — we re-gate here at the lower MIN_BLOCKCOMP_N threshold).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 03d: Extract per-(cohort × cell) GRM blocks (SLURM array) ==="

SLURM_SCRIPT="${PROJ_ROOT}/scripts/03d_extract_cohort_x_cell_grms.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# Sanity: per-cohort unrel GRMs from step 04c must exist.
for c in "${COHORTS_REML[@]}"; do
    if [[ ! -f "${GRM_OUT}/grm_cohort_${c}_unrel.grm.bin" ]]; then
        log_error "Missing per-cohort unrel GRM for ${c} — run step 04c first."
        exit 1
    fi
done

KEEP_DIR="${PROJ_ROOT}/data/keep"

# Enumerate all (cohort × cell) combinations. Emit one task per combination
# with N >= MIN_BLOCKCOMP_N.
TASKS_FILE="${LOG_DIR}/_step03d_blockcomps.txt"
: > "$TASKS_FILE"

log_info ""
log_info "Per-(cohort × cell) sub-sample sizes (analytic ∩ cohort-unrel; threshold N >= ${MIN_BLOCKCOMP_N}):"
for c in "${COHORTS_REML[@]}"; do
    for s in "${STRATA_MAIN[@]}"; do
        keep="${KEEP_DIR}/cohort_${c}_${s}.gcta.iids"
        if [[ ! -f "$keep" ]]; then
            log_info "  ${c} × ${s}: keep file missing (run step 05c first) — skipping"
            continue
        fi
        n=$(wc -l < "$keep")
        if (( n < MIN_BLOCKCOMP_N )); then
            log_info "  ${c} × ${s}: N=${n} < ${MIN_BLOCKCOMP_N} — skipping"
            continue
        fi
        out_prefix="${GRM_OUT}/grm_blockcomp_${c}_${s}"
        if [[ -f "${out_prefix}.grm.bin" ]]; then
            log_info "  ${c} × ${s}: N=${n} — SKIP (block already extracted)"
            continue
        fi
        log_info "  ${c} × ${s}: N=${n} — pending"
        # Tab-separated: cohort \t stratum \t in_prefix \t keep_path \t out_prefix
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$c" "$s" \
            "${GRM_OUT}/grm_cohort_${c}_unrel" \
            "$keep" \
            "$out_prefix" \
            >> "$TASKS_FILE"
    done
done

N=$(wc -l < "$TASKS_FILE")
if (( N == 0 )); then
    log_info ""
    log_info "All block components already extracted — nothing to submit."
    log_info "=== Step 03d complete (no-op) ==="
    exit 0
fi

log_info ""
log_info "Submitting SLURM array (1-${N}) for ${N} block component(s):"
i=0
while IFS=$'\t' read -r c s _; do
    i=$((i + 1))
    log_info "  array task ${i}: ${c} × ${s}"
done < "$TASKS_FILE"
log_info "Tasks file: $TASKS_FILE"

JOB_OUTPUT=$(sbatch \
    --wait \
    --array="1-${N}" \
    --export="ALL,TASKS_FILE=${TASKS_FILE}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

log_info ""
log_info "Built block components:"
while IFS=$'\t' read -r c s _ _ out_prefix; do
    if [[ -f "${out_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${out_prefix}.grm.id")
        log_info "  ${c} × ${s}: N=${n_id}"
    else
        log_warn "  ${c} × ${s}: missing (task may have failed — check logs/)"
    fi
done < "$TASKS_FILE"

log_info ""
log_info "=== Step 03d complete ==="
