#!/usr/bin/env bash
# =============================================================================
# run_blockdiag_region_main.sh — Track R1: Region main (East vs West)
# =============================================================================
# Secondary/exploratory block-diagonal GREML track. Compares East vs West
# regions collapsed across pre/post-1990. Covariates: gender_int + cohort_int.
#
# Pipeline order driven by env vars:
#   03e_extract_cohort_x_stratum_grms.sh   (subset-then-prune block components)
#   04d_build_block_diagonal_cell_grms.sh  (combine block components)
#   05d_fit_greml_block_diagonal.sh        (REML fit, with override covar file)
#   06d_assemble_block_diagonal.R          (TRACK-suffixed outputs)
#
# Outputs:
#   output/grm/grm_blockcomp_<COHORT>_region_{east,west}.*
#   output/grm/grm_blockdiag_region_{east,west}.*
#   output/reml/blockdiag_region_main/blockdiag_region_{east,west}.hsq
#   output/final/h2_region_main.{tsv,md}
#   output/final/h2_region_main_heterogeneity.tsv
#
# See Plan_deviations.md §7h for design rationale.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

export TRACK="region_main"
export STRATA_LIST_NAMES="${STRATA_REGION_MAIN[*]}"
export BLOCKDIAG_OUT_SUBDIR="blockdiag_region_main"
export BLOCKDIAG_COVAR_FILE_OVERRIDE="${BLOCKDIAG_COVAR_FILE}"

log_info "================================================================"
log_info "Track: ${TRACK} — Region main (East vs West)"
log_info "  Strata:    ${STRATA_LIST_NAMES}"
log_info "  Covar:     ${BLOCKDIAG_COVAR_FILE_OVERRIDE}"
log_info "  Out subdir: output/reml/${BLOCKDIAG_OUT_SUBDIR}/"
log_info "================================================================"

bash "${HERE}/03e_extract_cohort_x_stratum_grms.sh"
bash "${HERE}/04d_build_block_diagonal_cell_grms.sh"
bash "${HERE}/05d_fit_greml_block_diagonal.sh"
Rscript "${HERE}/06d_assemble_block_diagonal.R"

log_info ""
log_info "=== Track ${TRACK} complete ==="
