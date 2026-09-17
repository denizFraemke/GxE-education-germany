#!/usr/bin/env bash
# =============================================================================
# 03c_build_per_cohort_grms.sh — Submit per-cohort GRM builds as a SLURM array
# =============================================================================
# Per-cohort sensitivity track answering "is one cohort driving the pooled
# result?" (Plan_deviations.md §7e closing subsection). Mirrors step 03 but
# operates per-cohort, not per-stratum.
#
# For each COHORTS_REML entry (BASEII, SHIP = SHIP-0 + SHIP-Td, SOEP,
# TWINLIFE) we:
#   1. Build a (FID, IID) keep file from the cohort's harmonized fam(s),
#      intersected with the analytic sample (union of the 4 stratum keep
#      files), joined with pooled.fam for the FID column.
#   2. Submit a SLURM array task that runs gcta64 --bfile pooled --keep
#      --make-grm-bin. GCTA's within-keep allele-frequency standardization
#      means the resulting GRM has no cross-cohort fingerprint (the whole
#      point of this track).
#
# RESUME: tasks short-circuit if their grm_cohort_<cohort>.grm.bin already
#         exists.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 03c: Submit per-cohort GRM builds (SLURM array) ==="

POOLED="${GENO_OUT}/pooled"
if [[ ! -f "${POOLED}.bed" ]]; then
    log_error "Pooled bfile missing: ${POOLED}.bed"
    log_error "Run step 02 first."
    exit 1
fi

SLURM_SCRIPT="${PROJ_ROOT}/scripts/03c_build_per_cohort_grms.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# ---- 1. IID -> FID lookup from pooled.fam (reuse from step 03) ----------
FAM_INDEX="${GENO_OUT}/_pooled_iid_to_fid.tsv"
if [[ ! -s "$FAM_INDEX" ]]; then
    awk '{ print $2"\t"$1 }' "${POOLED}.fam" > "$FAM_INDEX"
    log_info "Wrote IID->FID index: $FAM_INDEX ($(wc -l < "$FAM_INDEX") rows)"
else
    log_info "Reusing IID->FID index: $FAM_INDEX ($(wc -l < "$FAM_INDEX") rows)"
fi

# ---- 2. Analytic-sample IID union from the 4 stratum keep files ----------
ANALYTIC_IIDS="${GENO_OUT}/_analytic_iids.tsv"
: > "$ANALYTIC_IIDS"
for s in "${STRATA_MAIN[@]}"; do
    keep="${PROJ_ROOT}/data/keep/${s}.iids"
    if [[ -f "$keep" ]]; then
        awk '{ print $2 }' "$keep" >> "$ANALYTIC_IIDS"
    fi
done
sort -u "$ANALYTIC_IIDS" -o "$ANALYTIC_IIDS"
log_info "Analytic-sample IIDs: $(wc -l < "$ANALYTIC_IIDS") (union of 4 cells)"

# ---- 3. Build per-cohort keep files -------------------------------------
KEEP_DIR="${PROJ_ROOT}/data/keep"
mkdir -p "$KEEP_DIR"

build_cohort_keep() {
    local cohort="$1"
    local components
    components=$(cohort_reml_components "$cohort")
    local out="${KEEP_DIR}/cohort_${cohort}.gcta.iids"

    # 3a. Union of components' IIDs from harmonized_r2_0.9 fam files.
    local tmp_cohort_iids; tmp_cohort_iids="$(mktemp)"
    : > "$tmp_cohort_iids"
    for c in $components; do
        local fam="${GENO_OUT}/${c}_harmonized_r2_0.9.fam"
        if [[ ! -f "$fam" ]]; then
            log_error "  Missing harmonized fam: $fam (run step 01 first)"
            rm -f "$tmp_cohort_iids"
            return 1
        fi
        awk '{ print $2 }' "$fam" >> "$tmp_cohort_iids"
    done
    sort -u "$tmp_cohort_iids" -o "$tmp_cohort_iids"

    # 3b. Intersect with analytic IIDs.
    local tmp_intersect; tmp_intersect="$(mktemp)"
    comm -12 "$tmp_cohort_iids" "$ANALYTIC_IIDS" > "$tmp_intersect"

    # 3c. Join with FAM_INDEX → (FID, IID) two-column file for GCTA --keep.
    join -t $'\t' -1 1 -2 1 \
        <(sort -t $'\t' -k1,1 "$tmp_intersect") \
        <(sort -t $'\t' -k1,1 "$FAM_INDEX") \
        | awk -F'\t' 'BEGIN{OFS="\t"} { print $2, $1 }' \
        > "$out"

    local n_out
    n_out=$(wc -l < "$out")
    rm -f "$tmp_cohort_iids" "$tmp_intersect"
    echo "$n_out"
}

log_info ""
log_info "Cohorts to build (${#COHORTS_REML[@]} entries):"
for c in "${COHORTS_REML[@]}"; do log_info "  - $c ($(cohort_reml_components "$c"))"; done

PENDING=()
for cohort in "${COHORTS_REML[@]}"; do
    grm_prefix="${GRM_OUT}/grm_cohort_${cohort}"
    out_keep="${KEEP_DIR}/cohort_${cohort}.gcta.iids"

    if [[ -f "${grm_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${grm_prefix}.grm.id")
        log_info "  SKIP: cohort_${cohort} GRM already built (N=${n_id})"
        continue
    fi

    log_info "  Building keep file for ${cohort}..."
    n_keep=$(build_cohort_keep "$cohort")
    log_info "    N in keep (analytic ∩ cohort, joined with pooled.fam): ${n_keep}"
    if (( n_keep < 1 )); then
        log_warn "    Empty keep for ${cohort}; will not submit."
        continue
    fi
    PENDING+=("$cohort")
done

if (( ${#PENDING[@]} == 0 )); then
    log_info ""
    log_info "All cohort GRMs already built — nothing to submit."
    log_info "=== Step 03c complete (no-op) ==="
    exit 0
fi

# ---- 4. Submit SLURM array ----------------------------------------------
COHORTS_LIST="${LOG_DIR}/_step03c_cohorts.txt"
: > "$COHORTS_LIST"
for c in "${PENDING[@]}"; do echo "$c" >> "$COHORTS_LIST"; done
N=${#PENDING[@]}

# Same beegfs/FUSE concurrency cap as step 03. SHIP is the largest cohort
# (~7,770 individuals) and reads the full 21 GB pooled.bed, so the same
# %2 throttle applies.
MAX_CONCURRENT="${MAX_CONCURRENT:-2}"

log_info ""
log_info "Submitting SLURM array (1-${N}%${MAX_CONCURRENT}) for ${N} cohort(s):"
for ((i=1; i<=N; i++)); do log_info "  array task ${i}: ${PENDING[$((i-1))]}"; done
log_info "Cohorts list: $COHORTS_LIST"
log_info "SLURM worker: $SLURM_SCRIPT"
log_info "Max concurrent tasks: ${MAX_CONCURRENT}"

JOB_OUTPUT=$(sbatch \
    --wait \
    --array="1-${N}%${MAX_CONCURRENT}" \
    --export="ALL,COHORTS_LIST=${COHORTS_LIST}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    log_error ""
    log_error "Re-running this script picks up only cohorts without a"
    log_error "grm_cohort_<cohort>.grm.bin file."
    exit 1
}
log_info "$JOB_OUTPUT"

# ---- 5. Report ----------------------------------------------------------
log_info ""
log_info "Built per-cohort GRMs:"
for cohort in "${COHORTS_REML[@]}"; do
    grm_prefix="${GRM_OUT}/grm_cohort_${cohort}"
    if [[ -f "${grm_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${grm_prefix}.grm.id")
        log_info "  ${cohort}: N=${n_id}"
    else
        log_warn "  ${cohort}: GRM missing (task may have failed — check logs/)"
    fi
done

log_info ""
log_info "=== Step 03c complete ==="
