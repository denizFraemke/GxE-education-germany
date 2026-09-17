#!/usr/bin/env bash
# =============================================================================
# 04c_unrelated_filter_per_cohort.sh — Submit per-cohort --grm-cutoff
# =============================================================================
# Per-cohort sensitivity track (Plan_deviations.md §7e closing subsection).
# Mirrors 04_unrelated_filter.sh but operates on the cohort GRMs built by
# step 03c, using GRM_CUTOFF_PER_COHORT (0.05) — back to the published
# convention since the cross-cohort fingerprint that drove the cell-level
# cutoff to 0.20 doesn't exist inside a single cohort.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04c: Submit per-cohort --grm-cutoff (SLURM array) ==="
log_info "Cutoff: ${GRM_CUTOFF_PER_COHORT} (within-cohort; published convention)"

SLURM_SCRIPT="${PROJ_ROOT}/scripts/04c_unrelated_filter_per_cohort.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

PENDING=()
for cohort in "${COHORTS_REML[@]}"; do
    grm_in="${GRM_OUT}/grm_cohort_${cohort}"
    grm_out="${GRM_OUT}/grm_cohort_${cohort}_unrel"

    if [[ ! -f "${grm_in}.grm.bin" ]]; then
        log_warn "  ${cohort}: input GRM missing (${grm_in}.grm.bin) — run step 03c first; skipping"
        continue
    fi
    if [[ -f "${grm_out}.grm.bin" ]]; then
        n_after=$(wc -l < "${grm_out}.grm.id")
        log_info "  SKIP: ${cohort} already filtered (N=${n_after})"
        continue
    fi
    PENDING+=("$cohort")
done

if (( ${#PENDING[@]} == 0 )); then
    log_info ""
    log_info "All per-cohort GRMs already filtered — nothing to submit."
    log_info "=== Step 04c complete (no-op) ==="
    exit 0
fi

COHORTS_LIST="${LOG_DIR}/_step04c_cohorts.txt"
: > "$COHORTS_LIST"
for c in "${PENDING[@]}"; do echo "$c" >> "$COHORTS_LIST"; done
N=${#PENDING[@]}

log_info ""
log_info "Submitting SLURM array (1-${N}) for ${N} cohort(s):"
for ((i=1; i<=N; i++)); do log_info "  array task ${i}: ${PENDING[$((i-1))]}"; done
log_info "Cohorts list: $COHORTS_LIST"

JOB_OUTPUT=$(sbatch \
    --wait \
    --array="1-${N}" \
    --export="ALL,COHORTS_LIST=${COHORTS_LIST}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

log_info ""
log_info "Post-cutoff per-cohort N:"
for c in "${COHORTS_REML[@]}"; do
    f="${GRM_OUT}/grm_cohort_${c}_unrel.grm.id"
    [[ -f "$f" ]] && log_info "  ${c}: $(wc -l < "$f")"
done

log_info ""
log_info "=== Step 04c complete ==="
