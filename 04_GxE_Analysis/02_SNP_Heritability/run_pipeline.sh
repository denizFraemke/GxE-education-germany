#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master orchestrator for the SNP-heritability pipeline
# =============================================================================
# RUNS ON TARDIS, AFTER the Mac side has uploaded code (via
# transfer_pipeline_to_tardis.sh) and phenotype/strata (via
# transfer_phenotype_to_tardis.sh).
#
# USAGE on Tardis:
#   cd ${TARDIS_PROJECT_ROOT}
#   bash run_pipeline.sh                # all steps
#   bash run_pipeline.sh --step 03      # one step
#   bash run_pipeline.sh --from 03 --to 06
#
# PIPELINE STEPS:
#   00   Check inputs            — gcta64/plink2 present; bfiles + phenotype reachable
#   01   Rebuild INFO>=.90 bfiles — re-run upstream 03a at R²>=0.9 + plink2 --extract
#   02   Pool cohorts            — plink 1.9 --merge-list with missnp retry
#   03   Build stratum GRMs      — per stratum gcta --make-grm-bin --keep
#   04   Unrelated filter        — per stratum gcta --grm-cutoff $GRM_CUTOFF
#   05   Fit GREML (4-cell)      — 8 per-stratum standalone REML fits
#   07   Heterogeneity tests     — Cochran's Q + pairwise Z on the standalones
#   08   Power calc              — Visscher 2014 SE(h²) + power on h² grid
#   09   Assemble output         — h2 + heterogeneity + power TSVs + dashboard row
#   (There is no step 06; the 05 → 07 numbering gap is intentional.)
#
# PER-COHORT SENSITIVITY TRACK (Plan_deviations.md §7e closing subsection):
#   03c  Build per-cohort GRMs   — one GRM per (BASE-II, SHIP, SOEP, TwinLife)
#   04c  Unrelated filter        — per cohort, cutoff 0.05 (within-cohort)
#   05c  Fit GREML per cohort    — cohort-only + cohort × region + cohort × stratum
#   06c  Assemble per-cohort h²  — h2_per_cohort.tsv + summary report
#   These steps are independent of 03/04/05 (no shared GRMs) and answer the
#   "is one cohort driving the pooled result?" question. They run after 02
#   and don't need 05 (FID alignment is duplicated in 05c).
#
# BLOCK-DIAGONAL CELL-LEVEL TRACK (preferred cell-level answer):
#   03d  Extract (cohort × cell) GRMs — sub-GRMs from each per-cohort unrel
#   04d  Combine into per-cell GRMs   — block-diagonal stacker (login-node R)
#   05d  Fit GREML per cell           — REML on block-diagonal cell GRMs
#   06d  Assemble cell-level h²       — h2_blockdiag.tsv + Q test + summary
#   Builds on 03c/04c/05c output (per-cohort unrel GRMs and the
#   per-(cohort × cell) keep files). The block-diagonal GRM sets cross-cohort
#   entries to zero, so V_G is identified entirely from within-cohort
#   relatedness signal — the methodologically clean cell-level estimate.
#
# SECONDARY/EXPLORATORY COMPARISON TRACKS (Plan_deviations.md §7h):
#   xtra Run R1 + G1 + RG1 (and the regression check that guards the primary)
#        Implemented by three wrapper scripts that set TRACK + STRATA_LIST_NAMES
#        + a track-specific covar override, then dispatch the same 03e → 04d →
#        05d → 06d chain. 03e (subset-then-prune) is used instead of 03d for
#        the secondary tracks. The primary track (03d → 04d → 05d → 06d) is
#        unchanged; a regression check verifies its h²/SE values before the
#        secondary tracks run.
#
# NOTES:
#   - Steps run sequentially on the Tardis login node (no SLURM array here;
#     GCTA REML is single-machine).
#   - Each step is idempotent: it skips work already done (resume mode).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
ensure_dirs

STEP_FROM=""
STEP_TO=""
SINGLE_STEP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step) SINGLE_STEP="$2"; shift 2 ;;
        --from) STEP_FROM="$2";   shift 2 ;;
        --to)   STEP_TO="$2";     shift 2 ;;
        --help|-h)
            echo "Usage: bash run_pipeline.sh [--step N] [--from N] [--to N]"
            echo "Steps: 00 01 02 03 04 05 07 08 09  (there is no step 06)"
            echo ""
            echo "Per-cohort sensitivity track (§7e closing subsection):"
            echo "  03c  build per-cohort GRMs"
            echo "  04c  per-cohort --grm-cutoff ${GRM_CUTOFF_PER_COHORT}"
            echo "  05c  per-cohort REML (cohort / region / stratum sub-fits)"
            echo "  06c  assemble per-cohort h² output"
            echo ""
            echo "Block-diagonal cell-level track (preferred cell-level answer):"
            echo "  03d  extract per-(cohort × cell) GRM blocks"
            echo "  04d  combine into per-cell block-diagonal GRMs"
            echo "  05d  per-cell block-diagonal REML"
            echo "  06d  assemble cell-level h² output"
            echo ""
            echo "Secondary / exploratory comparison tracks (§7h):"
            echo "  xtra regression-check primary + run R1 + G1 + RG1 wrappers"
            echo "       (regression_check_blockdiag_primary.sh + run_blockdiag_*.sh)"
            echo ""
            echo "Config (set in config.sh):"
            echo "  R2_THRESHOLD_H2=${R2_THRESHOLD_H2}"
            echo "  GRM_CUTOFF=${GRM_CUTOFF}"
            echo "  GRM_CUTOFF_PER_COHORT=${GRM_CUTOFF_PER_COHORT}"
            echo "  MIN_BLOCKCOMP_N=${MIN_BLOCKCOMP_N}"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

ALL_STEPS=(00 01 02 03 04 05 07 08 09 03c 04c 05c 06c 03d 04d 05d 06d xtra)

if [[ -n "$SINGLE_STEP" ]]; then
    STEPS=("$SINGLE_STEP")
elif [[ -n "$STEP_FROM" ]] || [[ -n "$STEP_TO" ]]; then
    STEP_FROM="${STEP_FROM:-00}"
    STEP_TO="${STEP_TO:-xtra}"
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

PIPELINE_LOG="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$(dirname "$PIPELINE_LOG")"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

log_info "============================================================"
log_info "SNP-heritability pipeline (Tardis)"
log_info "Started: $(date)"
log_info "PROJ_ROOT:          $PROJ_ROOT"
log_info "GENOTYPE_PROJ_ROOT: $GENOTYPE_PROJ_ROOT"
log_info "R2_THRESHOLD_H2:    $R2_THRESHOLD_H2"
log_info "GRM_CUTOFF:         $GRM_CUTOFF (pooled-GRM cell track)"
log_info "GRM_CUTOFF_PER_COHORT: $GRM_CUTOFF_PER_COHORT (per-cohort sensitivity track)"
log_info "Steps to run:       ${STEPS[*]}"
log_info "Log file:           $PIPELINE_LOG"
log_info "============================================================"

run_step() {
    local step_num="$1"
    local script="$2"
    local desc="$3"

    log_info ""
    log_info "============================================================"
    log_info "STEP ${step_num}: ${desc}"
    log_info "============================================================"

    local start_time end_time elapsed exit_code
    start_time=$(date +%s)

    if [[ "$script" == *.R ]]; then
        Rscript "$script"
        exit_code=$?
    else
        bash "$script"
        exit_code=$?
    fi

    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))

    if [[ $exit_code -eq 0 ]]; then
        log_info "STEP ${step_num} completed in ${elapsed}s"
    else
        log_error "STEP ${step_num} FAILED (exit code $exit_code) after ${elapsed}s"
        exit $exit_code
    fi
}

for step in "${STEPS[@]}"; do
    case "$step" in
        00) run_step 00 "${SCRIPT_DIR}/scripts/00_check_inputs.sh"           "Check all required inputs" ;;
        01) run_step 01 "${SCRIPT_DIR}/scripts/01_rebuild_bfiles_info090.sh" "Rebuild per-cohort INFO>=.90 bfiles" ;;
        02) run_step 02 "${SCRIPT_DIR}/scripts/02_pool_cohorts.sh"           "Pool the 5 cohort bfiles into one" ;;
        03) run_step 03 "${SCRIPT_DIR}/scripts/03_build_stratum_grms.sh"     "Build per-stratum GRMs" ;;
        04) run_step 04 "${SCRIPT_DIR}/scripts/04_unrelated_filter.sh"       "Apply --grm-cutoff per stratum" ;;
        05) run_step 05 "${SCRIPT_DIR}/scripts/05_fit_greml_main.sh"         "Fit GREML — 4-cell main" ;;
        07) run_step 07 "${SCRIPT_DIR}/scripts/07_lrt.R"                     "Cochran's Q heterogeneity tests" ;;
        08) run_step 08 "${SCRIPT_DIR}/scripts/08_power_calc.R"              "Visscher 2014 power per stratum" ;;
        09) run_step 09 "${SCRIPT_DIR}/scripts/09_assemble_output.R"         "Assemble h2/heterogeneity/power TSVs + dashboard" ;;
        03c) run_step 03c "${SCRIPT_DIR}/scripts/03c_build_per_cohort_grms.sh"        "Build per-cohort GRMs (sensitivity track)" ;;
        04c) run_step 04c "${SCRIPT_DIR}/scripts/04c_unrelated_filter_per_cohort.sh"  "Apply --grm-cutoff per cohort (sensitivity track)" ;;
        05c) run_step 05c "${SCRIPT_DIR}/scripts/05c_fit_greml_per_cohort.sh"         "Fit GREML per cohort + sub-fits (sensitivity track)" ;;
        06c) run_step 06c "${SCRIPT_DIR}/scripts/06c_assemble_per_cohort.R"           "Assemble per-cohort h² output (sensitivity track)" ;;
        03d) run_step 03d "${SCRIPT_DIR}/scripts/03d_extract_cohort_x_cell_grms.sh"           "Extract (cohort × cell) GRM blocks (block-diagonal track)" ;;
        04d) run_step 04d "${SCRIPT_DIR}/scripts/04d_build_block_diagonal_cell_grms.sh"       "Combine into per-cell block-diagonal GRMs (block-diagonal track)" ;;
        05d) run_step 05d "${SCRIPT_DIR}/scripts/05d_fit_greml_block_diagonal.sh"             "Fit GREML on block-diagonal cell GRMs (block-diagonal track)" ;;
        06d) run_step 06d "${SCRIPT_DIR}/scripts/06d_assemble_block_diagonal.R"               "Assemble cell-level h² output (block-diagonal track)" ;;
        xtra)
            run_step xtra-regcheck "${SCRIPT_DIR}/scripts/regression_check_blockdiag_primary.sh" "Regression check: primary cell-level h² unchanged"
            run_step xtra-R1       "${SCRIPT_DIR}/scripts/run_blockdiag_region_main.sh"          "Secondary track R1: Region main (East vs West)"
            run_step xtra-G1       "${SCRIPT_DIR}/scripts/run_blockdiag_gender_main.sh"          "Secondary track G1: recorded gender/sex main"
            run_step xtra-RG1      "${SCRIPT_DIR}/scripts/run_blockdiag_region_gender.sh"        "Secondary track RG1: Region × recorded gender/sex"
            ;;
        *)  log_error "Unknown step: $step"; exit 1 ;;
    esac
done

log_info ""
log_info "============================================================"
log_info "Pipeline complete: $(date)"
log_info "Final output: ${FINAL_OUT}/"
log_info "Pipeline log: $PIPELINE_LOG"
log_info "============================================================"
