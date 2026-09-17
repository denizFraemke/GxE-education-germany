#!/usr/bin/env bash
# =============================================================================
# 01_submit_all_sbayesr.sh — Submit SBayesR SLURM array jobs for all traits
# =============================================================================
# Submits one SLURM array job (chr1-22) per trait in TRAIT_NAMES. Skips traits
# whose 22 .snpRes files already exist on disk (resume mode).
#
# USAGE:
#   bash scripts/01_submit_all_sbayesr.sh
#   (typically via `bash run_pipeline.sh --step 01`)
#
# INPUTS:
#   ${PROJ_ROOT}/scripts/01_run_sbayesr.slurm
#   ${SUMSTATS_DIR}/<trait subdir>/<per-chr .ma>
#   ${LD_REF_DIR}/ukb50k_shrunk_chr{1..22}_mafpt01.ldm.sparse
#
# OUTPUTS:
#   ${LOG_DIR}/sbayesr_job_ids.txt — submitted SLURM job IDs
#   ${SBAYESR_OUT}/<TAG>/<TAG>_chr{1..22}_sbayesR.snpRes  (written async by the
#     SLURM array tasks; this script returns once submission is done)
#
# DEPENDS ON: step 00 (input checks).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

SLURM_SCRIPT="${PROJ_ROOT}/scripts/01_run_sbayesr.slurm"
JOB_IDS_FILE="${LOG_DIR}/sbayesr_job_ids.txt"
> "$JOB_IDS_FILE"

log_info "=== Step 01: Submitting SBayesR jobs for all traits ==="

# --- Submit standard traits ---------------------------------------------------
for trait in "${TRAIT_NAMES[@]}"; do
    tag="${TRAIT_SBAYESR_TAG[$trait]}"
    out_dir="${SBAYESR_OUT}/${tag}"

    # Check if already complete
    n_done=0
    for chr in $(seq 1 22); do
        [[ -f "${out_dir}/${tag}_chr${chr}_sbayesR.snpRes" ]] && n_done=$((n_done + 1))
    done
    if [[ $n_done -eq 22 ]]; then
        log_info "SKIP: $trait — all 22 chromosomes already complete"
        continue
    fi

    log_info "Submitting SBayesR for trait: $trait (tag: $tag, ${n_done}/22 done)"

    JOB_ID=$(sbatch \
        --job-name="sbr_${tag}" \
        --output="${LOG_DIR}/sbayesr_${tag}_chr%a.%j.out" \
        --error="${LOG_DIR}/sbayesr_${tag}_chr%a.%j.err" \
        --export="ALL,TRAIT=${trait},OUT_TAG=${tag}" \
        "$SLURM_SCRIPT" | awk '{print $NF}')

    echo "${tag}:${JOB_ID}" >> "$JOB_IDS_FILE"
    log_info "  -> Submitted job array $JOB_ID"
done

log_info ""
log_info "=== All SBayesR jobs submitted ==="
log_info "Job IDs saved to: $JOB_IDS_FILE"
log_info "Monitor with: squeue -u \$USER | grep sbr_"
log_info "When all complete, run: bash scripts/02_concat_weights.sh"
