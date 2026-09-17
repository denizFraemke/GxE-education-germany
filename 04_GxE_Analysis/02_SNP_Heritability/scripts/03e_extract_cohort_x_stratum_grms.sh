#!/usr/bin/env bash
# =============================================================================
# 03e_extract_cohort_x_stratum_grms.sh — Subset-then-prune per-(cohort × stratum)
# =============================================================================
# Block-component extractor for the secondary/exploratory tracks (R1, G1, RG1
# — see config.sh STRATA_REGION_MAIN / STRATA_GENDER_MAIN / STRATA_REGION_GENDER
# and Plan_deviations.md §7h). Differs from 03d in two ways:
#
#   1. Starts from the PRE-cutoff per-cohort GRM `grm_cohort_<COHORT>` written
#      by step 03c (rather than the POST-cutoff `..._unrel` GRM used by 03d).
#   2. Applies `--grm-cutoff 0.05` AFTER subsetting to the (cohort × stratum)
#      sub-block.
#
# Together these implement *subset-then-prune*: within-cohort relatives whose
# partner falls into a different comparison stratum (e.g. a cross-sex sib pair
# in TwinLife under G1) are not removed by the cell-level primary's earlier
# unrelated filter — they're retained inside the block they end up in. Each
# block then has its own within-block unrelated filter applied so REML still
# operates on an unrelated sub-sample within each block.
#
# The cell-level primary track (03d) deliberately stays on the
# *prune-then-subset* order, which is what its validated h² estimates are
# based on.
#
# Parameterisation (env vars set by scripts/run_blockdiag_*.sh):
#   STRATA_LIST_NAMES  whitespace-separated stratum labels (no default; if
#                      empty, this script is a no-op and exits)
#
# Output: grm_blockcomp_<COHORT>_<STRATUM>.{grm.bin,grm.id,grm.N.bin}
# Idempotent: skips any (cohort × stratum) whose output already exists.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

if [[ -z "${STRATA_LIST_NAMES:-}" ]]; then
    log_error "03e: STRATA_LIST_NAMES is empty — run via scripts/run_blockdiag_*.sh."
    log_error "  (03e is intended for secondary tracks; for the cell-level primary, use 03d.)"
    exit 1
fi
read -r -a __STRATA <<< "${STRATA_LIST_NAMES}"

log_info "=== Step 03e: Subset-then-prune per-(cohort × stratum) GRM blocks ==="
log_info "  Strata (${#__STRATA[@]}): ${__STRATA[*]}"
log_info "  Source GRMs:   grm_cohort_<COHORT> (pre-cutoff, from step 03c)"
log_info "  Within-block:  --grm-cutoff 0.05"

SLURM_SCRIPT="${PROJ_ROOT}/scripts/03e_extract_cohort_x_stratum_grms.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# Sanity: pre-cutoff per-cohort GRMs from step 03c must exist.
for c in "${COHORTS_REML[@]}"; do
    if [[ ! -f "${GRM_OUT}/grm_cohort_${c}.grm.bin" ]]; then
        log_error "Missing pre-cutoff per-cohort GRM for ${c} — run step 03c first."
        exit 1
    fi
done

KEEP_DIR="${PROJ_ROOT}/data/keep"

# Helper: build an FID-aligned per-(cohort × stratum) keep file. We can't pass
# the raw whole-stratum keep file straight to `gcta --keep` because GCTA
# matches on (FID, IID) pairs, and `export_stratum_phenotype.R` writes FIDs
# that don't match the cohort-pooled `pooled.fam` convention for BASE-II and
# TwinLife (their source FIDs differ; see 05d's align_input_fids for the same
# mismatch).
#
# For each (cohort × stratum):
#   - read the cohort's per-cohort GRM .grm.id (the canonical FID source)
#   - read the stratum keep file
#   - emit (FID_from_GRM, IID) for every IID that appears in both
# Returns the count of overlapping IIDs by printing the line count.
build_aligned_keep() {
    local cohort_grm_id="$1" stratum_keep="$2" out_path="$3"
    awk '
        NR==FNR { iid_to_fid[$2]=$1; next }
        ($2 in iid_to_fid) { print iid_to_fid[$2]"\t"$2 }
    ' "$cohort_grm_id" "$stratum_keep" > "$out_path"
    wc -l < "$out_path"
}

ALIGNED_KEEP_DIR="${LOG_DIR}/_step03e_aligned_keep"
mkdir -p "$ALIGNED_KEEP_DIR"

TASKS_FILE="${LOG_DIR}/_step03e_blockcomps.txt"
: > "$TASKS_FILE"

log_info ""
log_info "Per-(cohort × stratum) sub-sample sizes (FID-aligned; gate N >= ${MIN_BLOCKCOMP_N}):"
for c in "${COHORTS_REML[@]}"; do
    cohort_grm_id="${GRM_OUT}/grm_cohort_${c}.grm.id"
    for s in "${__STRATA[@]}"; do
        keep="${KEEP_DIR}/${s}.iids"
        if [[ ! -f "$keep" ]]; then
            log_info "  ${c} × ${s}: keep file missing (${keep}) — skipping"
            continue
        fi
        aligned_keep="${ALIGNED_KEEP_DIR}/keep_${c}_${s}.iids"
        n=$(build_aligned_keep "$cohort_grm_id" "$keep" "$aligned_keep")
        n=$((n + 0))
        if (( n < MIN_BLOCKCOMP_N )); then
            log_info "  ${c} × ${s}: pre-cutoff N=${n} < ${MIN_BLOCKCOMP_N} — skipping"
            continue
        fi
        out_prefix="${GRM_OUT}/grm_blockcomp_${c}_${s}"
        if [[ -f "${out_prefix}.grm.bin" ]]; then
            # An existing block component with N=0 is unusable; re-extract it.
            n_done=$(wc -l < "${out_prefix}.grm.id" 2>/dev/null || echo 0)
            if (( n_done > 0 )); then
                log_info "  ${c} × ${s}: pre-cutoff N=${n} — SKIP (already extracted, N=${n_done})"
                continue
            else
                log_warn "  ${c} × ${s}: existing block component has N=0 — forcing re-extract"
                rm -f "${out_prefix}.grm.bin" "${out_prefix}.grm.id" "${out_prefix}.grm.N.bin"
            fi
        fi
        log_info "  ${c} × ${s}: pre-cutoff N=${n} — pending"
        # Tab-separated: cohort \t stratum \t in_prefix \t keep_path \t out_prefix
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$c" "$s" \
            "${GRM_OUT}/grm_cohort_${c}" \
            "$aligned_keep" \
            "$out_prefix" \
            >> "$TASKS_FILE"
    done
done

N=$(wc -l < "$TASKS_FILE")
if (( N == 0 )); then
    log_info ""
    log_info "All block components already extracted (or below threshold) — nothing to submit."
    log_info "=== Step 03e complete (no-op) ==="
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
    --export="ALL,TASKS_FILE=${TASKS_FILE},GRM_CUTOFF_BLOCK=${GRM_CUTOFF_PER_COHORT}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

log_info ""
log_info "Built block components (post-cutoff N):"
while IFS=$'\t' read -r c s _ _ out_prefix; do
    if [[ -f "${out_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${out_prefix}.grm.id")
        log_info "  ${c} × ${s}: N=${n_id}"
    else
        log_warn "  ${c} × ${s}: missing (task may have failed — check logs/)"
    fi
done < "$TASKS_FILE"

log_info ""
log_info "=== Step 03e complete ==="
