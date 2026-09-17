#!/usr/bin/env bash
# =============================================================================
# run_blockdiag_region_gender.sh — Track RG1: Region × recorded gender/sex
# =============================================================================
# Secondary/exploratory block-diagonal GREML track. Compares the four
# Region × gender/sex groups: East × Male, East × Female, West × Male,
# West × Female. Covariates: cohort_int only.
#
# Pipeline order driven by env vars:
#   03e_extract_cohort_x_stratum_grms.sh   (subset-then-prune block components)
#   04d_build_block_diagonal_cell_grms.sh  (combine block components)
#   05d_fit_greml_block_diagonal.sh        (REML fit, with override covar file)
#   06d_assemble_block_diagonal.R          (TRACK-suffixed outputs)
#
# Outputs:
#   output/grm/grm_blockcomp_<COHORT>_region_gender_<region>_<sex>.*
#   output/grm/grm_blockdiag_region_gender_<region>_<sex>.*
#   output/reml/blockdiag_region_gender/blockdiag_region_gender_<region>_<sex>.hsq
#   output/final/h2_region_gender.{tsv,md}
#   output/final/h2_region_gender_heterogeneity.tsv
#
# See Plan_deviations.md §7h for design rationale.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

export TRACK="region_gender"
export STRATA_LIST_NAMES="${STRATA_REGION_GENDER[*]}"
export BLOCKDIAG_OUT_SUBDIR="blockdiag_region_gender"
export BLOCKDIAG_COVAR_FILE_OVERRIDE="${COHORT_ONLY_COVAR_FILE}"

log_info "================================================================"
log_info "Track: ${TRACK} — Region × recorded gender/sex"
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
