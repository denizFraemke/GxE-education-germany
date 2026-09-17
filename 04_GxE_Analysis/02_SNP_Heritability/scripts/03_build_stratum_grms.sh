#!/usr/bin/env bash
# =============================================================================
# 03_build_stratum_grms.sh — Submit per-stratum GRM builds as a SLURM array
# =============================================================================
# Mirrors 01_Genotype/Harmonized/scripts/01_submit_all_sbayesr.sh.
#
# For each of the 4 main cells (STRATA_MAIN) this script:
#   1. Builds the IID -> FID index from pooled.fam once.
#   2. Resolves each stratum's exporter-written "FID_placeholder IID" keep
#      file into a "FID_from_pooled IID" file (data/keep/<stratum>.gcta.iids).
#   3. Writes the strata list to a stable file in logs/.
#   4. Submits the SLURM array in scripts/03_build_stratum_grms.slurm with
#      one task per stratum, --wait, and blocks until all tasks finish.
#
# The actual gcta64 --make-grm-bin call runs on a compute node (sized by the
# slurm script's #SBATCH lines), not on the login node, whose cgroup CPU
# throttling stalls it.
#
# RESUME: tasks short-circuit if their grm_<stratum>.grm.bin already exists.
#         If --wait is interrupted (e.g. the ssh session drops), the SLURM
#         array keeps running — re-run this script later; already-done
#         strata will skip.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 03: Submit per-stratum GRM builds (SLURM array) ==="

POOLED="${GENO_OUT}/pooled"
if [[ ! -f "${POOLED}.bed" ]]; then
    log_error "Pooled bfile missing: ${POOLED}.bed"
    log_error "Run step 02 first."
    exit 1
fi

SLURM_SCRIPT="${PROJ_ROOT}/scripts/03_build_stratum_grms.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker script missing: $SLURM_SCRIPT"
    exit 1
fi

# ---- 1. IID -> FID index from pooled.fam (once) -------------------------
FAM_INDEX="${GENO_OUT}/_pooled_iid_to_fid.tsv"
if [[ ! -s "$FAM_INDEX" ]]; then
    awk '{ print $2"\t"$1 }' "${POOLED}.fam" > "$FAM_INDEX"
    log_info "Wrote IID->FID index: $FAM_INDEX ($(wc -l < "$FAM_INDEX") rows)"
else
    log_info "Reusing IID->FID index: $FAM_INDEX ($(wc -l < "$FAM_INDEX") rows)"
fi

# ---- 2. Resolve per-stratum keep files (cheap; do once on login node) ---
resolve_keep_to_gcta() {
    local in_keep="$1"
    local out_keep="$2"

    awk 'NF >= 2 { print $2 }' "$in_keep" | sort -u > "${out_keep}.iids"
    join -t $'\t' -1 1 -2 1 \
        <(sort -t $'\t' -k1,1 "${out_keep}.iids") \
        <(sort -t $'\t' -k1,1 "$FAM_INDEX") \
        | awk -F'\t' 'BEGIN{OFS="\t"} { print $2, $1 }' \
        > "$out_keep"

    local n_in n_out
    n_in=$(wc -l < "${out_keep}.iids")
    n_out=$(wc -l < "$out_keep")
    rm -f "${out_keep}.iids"

    if (( n_out < n_in )); then
        log_warn "  $(basename "$in_keep"): $((n_in - n_out)) IID(s) not in pooled.fam, dropped"
    fi
    echo "$n_out"
}

# ---- 3. Stratum list (4-cell main only) ---------------------------------
STRATA_ALL=("${STRATA_MAIN[@]}")

log_info ""
log_info "Strata to build (${#STRATA_ALL[@]} cells, 4-cell main analysis):"
for s in "${STRATA_ALL[@]}"; do log_info "  - $s"; done

# Resolve keep files + collect strata that still need building.
PENDING=()
for stratum in "${STRATA_ALL[@]}"; do
    in_keep="$(get_keep_file "$stratum")"
    gcta_keep="${PROJ_ROOT}/data/keep/${stratum}.gcta.iids"
    grm_prefix="${GRM_OUT}/grm_${stratum}"

    if [[ ! -f "$in_keep" ]]; then
        log_warn "  ${stratum}: no exporter-written keep file at $in_keep — skipping"
        continue
    fi

    if [[ -f "${grm_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${grm_prefix}.grm.id")
        log_info "  SKIP: ${stratum} GRM already built (N=${n_id})"
        continue
    fi

    log_info "  Resolving keep for ${stratum}..."
    n_keep=$(resolve_keep_to_gcta "$in_keep" "$gcta_keep")
    log_info "    N in keep (after pooled.fam join): ${n_keep}"
    if (( n_keep < 1 )); then
        log_warn "    Empty keep for ${stratum}; will not submit."
        continue
    fi
    PENDING+=("$stratum")
done

if (( ${#PENDING[@]} == 0 )); then
    log_info ""
    log_info "All strata already built — nothing to submit."
    log_info "=== Step 03 complete (no-op) ==="
    exit 0
fi

# ---- 4. Write the SLURM-array input list --------------------------------
STRATA_LIST="${LOG_DIR}/_step03_strata.txt"
: > "$STRATA_LIST"
for s in "${PENDING[@]}"; do echo "$s" >> "$STRATA_LIST"; done
N=${#PENDING[@]}

# Concurrency throttle: at most this many tasks running at once across the
# array. Combined with --exclusive in the .slurm worker, this caps how many
# nodes concurrently mmap-read the beegfs-hosted pooled.bed. Without the
# throttle, 3-4 gcta64 processes per node block in fuse_d kernel waits —
# the FUSE/beegfs mount cannot serve that many concurrent readers.
MAX_CONCURRENT="${MAX_CONCURRENT:-2}"

log_info ""
log_info "Submitting SLURM array (1-${N}%${MAX_CONCURRENT}) for ${N} stratum/strata:"
for ((i=1; i<=N; i++)); do log_info "  array task ${i}: ${PENDING[$((i-1))]}"; done
log_info "Strata list: $STRATA_LIST"
log_info "SLURM worker: $SLURM_SCRIPT"
log_info "Max concurrent tasks: ${MAX_CONCURRENT} (export MAX_CONCURRENT=K to override)"

# ---- 5. sbatch --wait ---------------------------------------------------
# --wait blocks this script until every array task has completed. If the
# session dies (e.g. the ssh connection drops), the SLURM array keeps
# running on the cluster; re-running this script picks up where it left off
# because completed strata short-circuit at the top.
JOB_OUTPUT=$(sbatch \
    --wait \
    --array="1-${N}%${MAX_CONCURRENT}" \
    --export="ALL,STRATA_LIST=${STRATA_LIST}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    log_error ""
    log_error "Re-running this script will resubmit only the unfinished strata"
    log_error "(those without a grm_<stratum>.grm.bin file)."
    exit 1
}

log_info "$JOB_OUTPUT"

# ---- 6. Report what got built -------------------------------------------
log_info ""
log_info "Built GRMs:"
for stratum in "${STRATA_ALL[@]}"; do
    grm_prefix="${GRM_OUT}/grm_${stratum}"
    if [[ -f "${grm_prefix}.grm.bin" ]]; then
        n_id=$(wc -l < "${grm_prefix}.grm.id")
        log_info "  ${stratum}: N=${n_id}"
    else
        log_warn "  ${stratum}: GRM missing (task may have failed — check logs/)"
    fi
done

log_info ""
log_info "=== Step 03 complete ==="
