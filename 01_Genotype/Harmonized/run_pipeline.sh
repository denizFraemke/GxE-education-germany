#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master orchestrator for PGI Computation Pipeline
#                  WITH Cross-Cohort Harmonization
# =============================================================================
# Runs all pipeline steps in order, including harmonization steps.
#
# USAGE:
#   # Run everything:
#   bash run_pipeline.sh
#
#   # Run a specific step:
#   bash run_pipeline.sh --step 00
#
#   # Run steps 03 through 07 (after SBayesR is done):
#   bash run_pipeline.sh --from 03 --to 07
#
# PIPELINE STEPS:
#   00   Check inputs         — verify all required files exist
#   01   SBayesR              — submit SLURM array jobs (per chr, per trait)
#   02   Concat weights       — concatenate per-chr SBayesR outputs
#   03a  R2 filter            — imputation R² quality filtering (creates pass lists)
#   03   Harmonize geno       — chr:pos IDs + unified SNP set restriction (uses R2 pass lists)
#   03b  Fix strand           — fix allele orientation for NonCog/Cog
#   03c  UKB align            — align to UKB reference (optional)
#   04   Compute PCs           — standard ancestry PCs per cohort
#   04b  Cross-cohort PCs      — ancestry PCs on the pooled sample
#   05   Score PGIs           — PLINK2 --score for all traits
#   06b  R2 sensitivity       — score PGIs with strict R2 > 0.8 SNP set
#   07   Assemble output       — cross-cohort harmonization
#   07b  Compare approaches    — compare harmonization approaches
#
# NOT DRIVEN BY THIS ORCHESTRATOR:
#   04c  Ancestry keep-lists   — scripts/04c_ancestry_1kg.sh followed by
#        scripts/04c_classify_eur.R. Projects each cohort onto 1000G phase-3
#        PCs and writes the European-ancestry keep-lists that 03_Merge
#        consumes. Depends on step 03; run it by hand after step 03/04.
#        See scripts/README_04c_ancestry.md.
#
# NOTES:
#   - Step 01 submits SLURM jobs that run asynchronously.
#   - Steps 02–07b run sequentially on the login node.
#   - Each step is idempotent: it skips work already done (resume mode).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
ensure_dirs

# Harmonization settings (from config.sh)
FIX_STRAND="${FIX_STRAND:-yes}"
SHARED_SNP_SET="${SHARED_SNP_SET:-yes}"
UKB_ALIGN="${UKB_ALIGN:-no}"

# --- Parse arguments ---------------------------------------------------------
STEP_FROM=""
STEP_TO=""
SINGLE_STEP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step)      SINGLE_STEP="$2"; shift 2 ;;
        --from)      STEP_FROM="$2"; shift 2 ;;
        --to)        STEP_TO="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: bash run_pipeline.sh [--step N] [--from N] [--to N]"
            echo "Steps: 00 01 02 03a 03 03b 03c 04 04b 05 06b 07 07b"
            echo ""
            echo "Harmonization settings (set in config.sh):"
            echo "  FIX_STRAND=$FIX_STRAND"
            echo "  SHARED_SNP_SET=$SHARED_SNP_SET"
            echo "  R2_FILTER=${R2_FILTER:-yes}"
            echo "  R2_THRESHOLD=${R2_THRESHOLD:-0.3}"
            echo "  R2_THRESHOLD_STRICT=${R2_THRESHOLD_STRICT:-0.8}"
            echo "  UKB_ALIGN=$UKB_ALIGN"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# All steps in order
ALL_STEPS=(00 01 02)
if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
    ALL_STEPS+=(03a)  # R2 quality filtering — MUST run before step 03 Phase 2
fi
ALL_STEPS+=(03)
if [[ "$FIX_STRAND" == "yes" ]]; then
    ALL_STEPS+=(03b)
fi
if [[ "$UKB_ALIGN" == "yes" ]]; then
    ALL_STEPS+=(03c)
fi
ALL_STEPS+=(04 04b 05)
if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
    ALL_STEPS+=(06b)
fi
ALL_STEPS+=(07 07b)

# Determine which steps to run
if [[ -n "$SINGLE_STEP" ]]; then
    STEPS=("$SINGLE_STEP")
elif [[ -n "$STEP_FROM" ]] || [[ -n "$STEP_TO" ]]; then
    STEP_FROM="${STEP_FROM:-00}"
    STEP_TO="${STEP_TO:-07b}"
    STEPS=()
    in_range=false
    for s in "${ALL_STEPS[@]}"; do
        [[ "$s" == "$STEP_FROM" ]] && in_range=true
        $in_range && STEPS+=("$s")
        [[ "$s" == "$STEP_TO" ]] && break
    done
else
    STEPS=("${ALL_STEPS[@]}")
fi

# --- Pipeline log ------------------------------------------------------------
PIPELINE_LOG="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$(dirname "$PIPELINE_LOG")"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

log_info "============================================================"
log_info "PGI Computation Pipeline — Cross-Cohort Harmonized"
log_info "Started: $(date)"
log_info "Project root: $PROJ_ROOT"
log_info "Steps to run: ${STEPS[*]}"
log_info "Harmonization settings:"
log_info "  FIX_STRAND=$FIX_STRAND"
log_info "  SHARED_SNP_SET=$SHARED_SNP_SET"
log_info "  R2_FILTER=${R2_FILTER:-yes}"
log_info "  R2_THRESHOLD=${R2_THRESHOLD:-0.3}"
log_info "  R2_THRESHOLD_STRICT=${R2_THRESHOLD_STRICT:-0.8}"
log_info "  UKB_ALIGN=$UKB_ALIGN"
log_info "Log file: $PIPELINE_LOG"
log_info "============================================================"

# --- Run steps ---------------------------------------------------------------
run_step() {
    local step_num="$1"
    local script="$2"
    local desc="$3"

    log_info ""
    log_info "============================================================"
    log_info "STEP ${step_num}: ${desc}"
    log_info "============================================================"

    local start_time
    start_time=$(date +%s)

    bash "$script"
    local exit_code=$?

    local end_time
    end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))

    if [[ $exit_code -eq 0 ]]; then
        log_info "STEP ${step_num} completed in ${elapsed}s"
    else
        log_error "STEP ${step_num} FAILED (exit code $exit_code) after ${elapsed}s"
        exit $exit_code
    fi
}

run_step_rscript() {
    local step_num="$1"
    local script="$2"
    local desc="$3"

    log_info ""
    log_info "============================================================"
    log_info "STEP ${step_num}: ${desc}"
    log_info "============================================================"

    local start_time
    start_time=$(date +%s)

    Rscript "$script"
    local exit_code=$?

    local end_time
    end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))

    if [[ $exit_code -eq 0 ]]; then
        log_info "STEP ${step_num} completed in ${elapsed}s"
    else
        log_error "STEP ${step_num} FAILED (exit code $exit_code) after ${elapsed}s"
        exit $exit_code
    fi
}

for step in "${STEPS[@]}"; do
    case "$step" in
        00)
            run_step 00 "${SCRIPT_DIR}/scripts/00_check_inputs.sh" \
                "Check all required input files" ;;
        01)
            log_info ""
            log_info "============================================================"
            log_info "STEP 01: Submit SBayesR SLURM jobs"
            log_info "============================================================"
            bash "${SCRIPT_DIR}/scripts/01_submit_all_sbayesr.sh"
            log_info ""
            log_info "SBayesR jobs submitted. They will run asynchronously."
            log_info "Monitor with: squeue -u \$USER | grep sbr_"
            log_info ""
            log_info "IMPORTANT: Wait for ALL SBayesR jobs to complete before"
            log_info "running step 02. You can check with:"
            log_info "  squeue -u \$USER | grep sbr_ | wc -l"
            log_info ""
            log_info "Then continue with:"
            log_info "  bash run_pipeline.sh --from 02"
            # If running the full pipeline, pause here
            if [[ ${#STEPS[@]} -gt 1 ]] && [[ "${STEPS[0]}" == "00" || "${STEPS[0]}" == "01" ]]; then
                log_info "Stopping here to wait for SLURM jobs."
                log_info "Resume with: bash run_pipeline.sh --from 02"
                exit 0
            fi
            ;;
        02)
            run_step 02 "${SCRIPT_DIR}/scripts/02_concat_weights.sh" \
                "Concatenate SBayesR outputs and prepare scoring weights" ;;
        03)
            run_step 03 "${SCRIPT_DIR}/scripts/03_harmonize_geno.sh" \
                "Harmonize genotype variant IDs + unified SNP set" ;;
        03a)
            if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
                run_step 03a "${SCRIPT_DIR}/scripts/03a_filter_imputation_quality.sh" \
                    "Imputation R² quality filtering"
            else
                log_info "STEP 03a skipped (R2_FILTER=no)"
            fi
            ;;
        03b)
            if [[ "$FIX_STRAND" == "yes" ]]; then
                run_step 03b "${SCRIPT_DIR}/scripts/03b_fix_strand.sh" \
                    "Fix strand/allele orientation issues"
            else
                log_info "STEP 03b skipped (FIX_STRAND=no)"
            fi
            ;;
        03c)
            if [[ "$UKB_ALIGN" == "yes" ]]; then
                run_step 03c "${SCRIPT_DIR}/scripts/03c_align_to_ukb.sh" \
                    "Align genotypes to UKB reference panel"
            else
                log_info "STEP 03c skipped (UKB_ALIGN=no)"
            fi
            ;;
        04)
            run_step 04 "${SCRIPT_DIR}/scripts/04_compute_pcs.sh" \
                "Compute standard ancestry PCs" ;;
        04b)
            run_step 04b "${SCRIPT_DIR}/scripts/04b_compute_cross_cohort_pcs.sh" \
                "Compute cross-cohort PCs (pooled sample)" ;;
        05)
            run_step 05 "${SCRIPT_DIR}/scripts/05_score_pgi.sh" \
                "Score PGIs with PLINK2" ;;
        06b)
            if [[ "${R2_FILTER:-yes}" == "yes" ]]; then
                run_step 06b "${SCRIPT_DIR}/scripts/06b_r2_sensitivity.sh" \
                    "R² sensitivity scoring (R2 >= ${R2_THRESHOLD_STRICT:-0.8})"
            else
                log_info "STEP 06b skipped (R2_FILTER=no)"
            fi
            ;;
        07)
            run_step_rscript 07 "${SCRIPT_DIR}/scripts/07_assemble_output.R" \
                "Assemble output with cross-cohort harmonization" ;;
        07b)
            run_step_rscript 07b "${SCRIPT_DIR}/scripts/07b_compare_harmonization.R" \
                "Compare harmonization approaches" ;;
        *)
            log_error "Unknown step: $step"
            exit 1
            ;;
    esac
done

log_info ""
log_info "============================================================"
log_info "Pipeline complete: $(date)"
log_info "Final output: ${FINAL_OUT}/"
log_info "Diagnostic output: ${OUT_ROOT}/diagnostics/"
log_info "Pipeline log: $PIPELINE_LOG"
log_info "============================================================"
