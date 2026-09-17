#!/usr/bin/env bash
# =============================================================================
# 04_unrelated_filter.sh — Submit per-stratum --grm-cutoff as a SLURM array
# =============================================================================
# Mirrors 03_build_stratum_grms.sh. One array task = one stratum's
# unrelated-set filter via gcta --grm-cutoff.
#
# RESUME: tasks short-circuit if the *_unrel GRM already exists.
#         If --wait is interrupted, the array keeps running. Re-run this
#         script to pick up only what's missing.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04: Submit per-stratum --grm-cutoff (SLURM array) ==="

SLURM_SCRIPT="${PROJ_ROOT}/scripts/04_unrelated_filter.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# Stratum list (4-cell main only).
STRATA_ALL=("${STRATA_MAIN[@]}")

# Collect pending strata (those without an existing *_unrel GRM).
PENDING=()
for stratum in "${STRATA_ALL[@]}"; do
    grm_in="${GRM_OUT}/grm_${stratum}"
    grm_out="${GRM_OUT}/grm_${stratum}_unrel"

    if [[ ! -f "${grm_in}.grm.bin" ]]; then
        log_warn "  ${stratum}: input GRM missing (${grm_in}.grm.bin) — run step 03 first; skipping"
        continue
    fi
    if [[ -f "${grm_out}.grm.bin" ]]; then
        n_after=$(wc -l < "${grm_out}.grm.id")
        log_info "  SKIP: ${stratum} already filtered (N=${n_after})"
        continue
    fi
    PENDING+=("$stratum")
done

if (( ${#PENDING[@]} == 0 )); then
    log_info ""
    log_info "All strata already filtered — nothing to submit."
    log_info "=== Step 04 complete (no-op) ==="
    exit 0
fi

STRATA_LIST="${LOG_DIR}/_step04_strata.txt"
: > "$STRATA_LIST"
for s in "${PENDING[@]}"; do echo "$s" >> "$STRATA_LIST"; done
N=${#PENDING[@]}

log_info ""
log_info "Submitting SLURM array (1-${N}) for ${N} stratum/strata:"
for ((i=1; i<=N; i++)); do log_info "  array task ${i}: ${PENDING[$((i-1))]}"; done
log_info "Strata list: $STRATA_LIST"

JOB_OUTPUT=$(sbatch \
    --wait \
    --array="1-${N}" \
    --export="ALL,STRATA_LIST=${STRATA_LIST}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

# --- Per-stratum cell-N report ---
log_info ""
log_info "Post-cutoff cell N:"
for s in "${STRATA_MAIN[@]}"; do
    f="${GRM_OUT}/grm_${s}_unrel.grm.id"
    [[ -f "$f" ]] && log_info "  ${s}: $(wc -l < "$f")"
done

log_info ""
log_info "=== Step 04 complete ==="
