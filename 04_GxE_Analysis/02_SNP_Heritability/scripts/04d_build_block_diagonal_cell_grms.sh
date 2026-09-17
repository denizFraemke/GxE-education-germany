#!/usr/bin/env bash
# =============================================================================
# 04d_build_block_diagonal_cell_grms.sh — Combine block components into cells
# =============================================================================
# For each of the 4 main cells, take the per-(cohort × cell) block components
# built by step 03d and combine them into a single block-diagonal cell GRM
# using R/build_block_diagonal_grm.R. The cross-cohort entries of the
# resulting GRM are zero by construction → no cross-cohort allele-frequency
# fingerprint can drive REML or relatedness pruning.
#
# Runs on the login node (each cell's combine is a single R process that
# allocates one packed lower-triangular vector of size up to ~25M float32
# entries ≈ 100 MB, then writes it out). No SLURM array — the largest cell
# takes a couple of minutes; the four cells together complete in well under
# 10 minutes.
#
# Idempotent: skips cells whose blockdiag GRM is already on disk.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04d: Combine block components into per-cell block-diagonal GRMs ==="

COMBINER="${PROJ_ROOT}/R/build_block_diagonal_grm.R"
if [[ ! -f "$COMBINER" ]]; then
    log_error "Combiner R script missing: $COMBINER"
    exit 1
fi

# Parameterised stratum list. Default = STRATA_MAIN (validated cell-level
# primary). Secondary tracks (R1/G1/RG1) override this via STRATA_LIST env
# var; see scripts/run_blockdiag_*.sh. Existing behaviour is unchanged
# when the variable is unset.
if [[ -n "${STRATA_LIST_NAMES:-}" ]]; then
    # Whitespace-separated name list (set by wrappers via "${STRATA_LIST[*]}")
    read -r -a __STRATA <<< "${STRATA_LIST_NAMES}"
else
    __STRATA=("${STRATA_MAIN[@]}")
fi
log_info "  Stratum list (${#__STRATA[@]}): ${__STRATA[*]}"

for stratum in "${__STRATA[@]}"; do
    out_prefix="${GRM_OUT}/grm_blockdiag_${stratum}"

    if [[ -f "${out_prefix}.grm.bin" ]]; then
        n_done=$(wc -l < "${out_prefix}.grm.id")
        log_info "  SKIP: ${stratum} block-diagonal GRM already built (N=${n_done})"
        continue
    fi

    # Collect the block components for this cell.
    components=()
    log_info ""
    log_info "Building block-diagonal GRM for ${stratum}..."
    for cohort in "${COHORTS_REML[@]}"; do
        comp="${GRM_OUT}/grm_blockcomp_${cohort}_${stratum}"
        if [[ -f "${comp}.grm.bin" ]]; then
            n_comp=$(wc -l < "${comp}.grm.id")
            components+=("$comp")
            log_info "  + ${cohort} block: N=${n_comp}"
        fi
    done

    if (( ${#components[@]} == 0 )); then
        log_warn "  ${stratum}: no block components present — skipping"
        log_warn "    Run step 03d first."
        continue
    fi

    log_info "  Combining ${#components[@]} block(s) into ${out_prefix}..."
    Rscript "$COMBINER" "$out_prefix" "${components[@]}"

    if [[ -f "${out_prefix}.grm.bin" ]]; then
        n_total=$(wc -l < "${out_prefix}.grm.id")
        log_info "  ${stratum}: combined N=${n_total} → ${out_prefix}.grm.{bin,id,N.bin}"
    else
        log_error "  ${stratum}: block-diagonal GRM not written"
        exit 1
    fi
done

log_info ""
log_info "Per-cell block-diagonal GRM N (after combine):"
for s in "${__STRATA[@]}"; do
    f="${GRM_OUT}/grm_blockdiag_${s}.grm.id"
    [[ -f "$f" ]] && log_info "  ${s}: $(wc -l < "$f")"
done

log_info ""
log_info "=== Step 04d complete ==="
