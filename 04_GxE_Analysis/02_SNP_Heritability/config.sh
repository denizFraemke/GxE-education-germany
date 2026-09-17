#!/usr/bin/env bash
# =============================================================================
# config.sh — Central configuration for the SNP-heritability pipeline (Tardis)
# =============================================================================
# All paths, tool locations, threshold parameters, and stratum definitions.
# Sourced by every Tardis-side script under scripts/.
#
# ENVIRONMENT: MPIB Tardis HPC. Paths assume the SNP-h² project has been
# rsynced to ${PROJ_ROOT} via transfer_pipeline_to_tardis.sh, and the
# upstream PGI pipeline already lives at ${GENOTYPE_PROJ_ROOT}.
#
# This file is sourced, not executed: it deliberately does not set shell
# options, so that sourcing it from an interactive shell leaves errexit /
# nounset untouched. Every calling script sets `set -euo pipefail` itself.
# =============================================================================

# --- Project root (auto-detected from this file's location) -----------------
PROJ_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ_ROOT

# --- Upstream genotype pipeline location -------------------------------------
# Reused: harmonized bfiles, info files. Set in user_config.sh on the Mac and
# rsynced over; the fallback below keeps a Tardis-only checkout working when
# user_config.sh has not been transferred.
export GENOTYPE_PROJ_ROOT="${TARDIS_GENOTYPE_PROJECT_ROOT:-/home/mpib/${USER}/80_years_GxE_on_Education_in_Germany/genotype}"
export GENOTYPE_OUT_GENO="${GENOTYPE_PROJ_ROOT}/output/geno"
export GENOTYPE_INFO_DIR="${GENOTYPE_PROJ_ROOT}/data/info_files"

# --- Output / log directories (created automatically) ------------------------
export OUT_ROOT="${PROJ_ROOT}/output"
export GENO_OUT="${OUT_ROOT}/geno"
export GRM_OUT="${OUT_ROOT}/grm"
export REML_OUT="${OUT_ROOT}/reml"
export FINAL_OUT="${OUT_ROOT}/final"
export LOG_DIR="${PROJ_ROOT}/logs"

# --- Tool binaries -----------------------------------------------------------
# GCTA: REML / GRM. PLINK 1.9: cohort merge (handles missnp retry). PLINK 2.0:
# variant extraction. Override paths in user_config.sh if your install differs.
export GCTA="${TARDIS_BIN_DIR:-/data/home/${USER}/bin}/gcta64"
export PLINK19="${TARDIS_BIN_DIR:-/data/home/${USER}/bin}/plink"
export PLINK2="${TARDIS_BIN_DIR:-/data/home/${USER}/bin}/plink2"

# --- Cohort list -------------------------------------------------------------
export COHORTS=(BASEII SHIP0 SHIPTD SOEP TWINLIFE)

# --- Imputation R² threshold for this pipeline -------------------------------
# Plan S5: INFO >= 0.90 (stricter than the PGI pipeline's 0.30).
export R2_THRESHOLD_H2="0.9"

# --- GRM / unrelated-filter parameters --------------------------------------
# MAF / autosome filters per plan S5. Plan also calls for HWE p > 1e-4, but
# we skip that here — see scripts/03_build_stratum_grms.sh and
# Plan_deviations.md §7d for the rationale. (Briefly: GCTA's GRM step has
# no --hwe flag, and pooled-cohort HWE is confounded with cross-cohort
# frequency structure.)
export GREML_MAF="0.01"
# gcta --grm-cutoff; pairwise GRM threshold for the within-stratum unrelated
# filter on the pooled cross-cohort GRM. Plan S5 (line 247) specifies 0.025;
# we use 0.20.
#
# On the pooled GRM, cross-cohort allele-frequency structure inflates the
# off-diagonal tail, and GCTA's iterative removal drops the individuals who
# carry the most such pairs. At thresholds of 0.025-0.10 that removal takes
# out the *minority cohort* of each Region × Reunification cell almost
# entirely (>99% loss for TwinLife-East, SOEP-East and BASE-II-West at 0.05;
# the pre-1990 cells are still minority-depleted at 0.10 because ~7,000 SHIP
# individuals in East × pre-1990 give each minority individual a large tail
# of above-cutoff cross-cohort pairs). The resulting h^2 would describe
# "SHIP vs SOEP+TwinLife" rather than East vs West.
#
# At 0.20 the cutoff admits 2nd-cousin-level pairs (relatedness ~0.03-0.06)
# but preserves all four cells uniformly across cohorts; see
# Plan_deviations.md §7e for the per-cell retention tables. The retained
# sample therefore carries substantial relatedness, which is why this track
# is reported as a sensitivity rather than as a conventional
# unrelated-sample GREML.
export GRM_CUTOFF="0.20"

# --- Stratum definitions -----------------------------------------------------
# 4 cells = East/West × pre/post-reunification, where reunification cutoff =
# "turned 15 in or after 1990" → birth_year >= 1975.
export STRATA_MAIN=(
    East_pre1990
    East_post1990
    West_pre1990
    West_post1990
)

# --- Power-calc h² grid -----------------------------------------------------
# Visscher 2014 SE(h²) and LRT power are evaluated at each of these h² values
# per stratum, against the empirical Var(off-diag-GRM).
export H2_GRID=(0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40)

# --- Logging helpers ---------------------------------------------------------
log_msg()   { local level="$1"; shift; echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $*"; }
log_info()  { log_msg "INFO"  "$@"; }
log_warn()  { log_msg "WARN"  "$@"; }
log_error() { log_msg "ERROR" "$@"; }
export -f log_msg log_info log_warn log_error

ensure_dirs() {
    mkdir -p "$GENO_OUT" "$GRM_OUT" "$REML_OUT" "$FINAL_OUT" "$LOG_DIR"
}
export -f ensure_dirs

# --- Stratum → keep-list path helper ----------------------------------------
# Per-stratum two-column "FID IID" files are produced by the Mac-side
# export_stratum_phenotype.R and rsynced into ${PROJ_ROOT}/data/keep/.
get_keep_file() {
    local stratum="$1"
    echo "${PROJ_ROOT}/data/keep/${stratum}.iids"
}
export -f get_keep_file

# --- Per-cohort h² sensitivity track (steps 03c / 04c / 05c / 06c) -----------
# COHORTS_REML is the set of "cohorts" we fit per-cohort GREML on, treating
# SHIP-0 and SHIP-Td as a single combined SHIP (they're the same source
# population sampled in two waves; pooling them is the natural choice and
# also gives us a usefully large SHIP-East sample). See Plan_deviations.md
# §7e for the substantive question this track answers (whether one specific
# cohort is responsible for the pooled cell-level estimates).
export COHORTS_REML=(BASEII SHIP SOEP TWINLIFE)

# Map each REML cohort to the harmonized-fam IDs it comprises. Whitespace
# separated; consumed by 03c when building the (FID, IID) keep file.
cohort_reml_components() {
    case "$1" in
        BASEII)   echo "BASEII" ;;
        SHIP)     echo "SHIP0 SHIPTD" ;;
        SOEP)     echo "SOEP" ;;
        TWINLIFE) echo "TWINLIFE" ;;
        *) return 1 ;;
    esac
}
export -f cohort_reml_components

# Within-cohort GRM cutoff. Inside a single cohort, the
# pooled-allele-frequency fingerprint that drove the cell-level analysis to
# 0.20 doesn't exist (§7e). Back to the published convention.
export GRM_CUTOFF_PER_COHORT="0.05"

# --- Block-diagonal cell-level h² track (steps 03d / 04d / 05d / 06d) --------
# Builds the methodologically cleaner cell-level h² answer: per-(cohort × cell)
# GRMs combined block-diagonally so cross-cohort pairs contribute zero to V_G.
# Within-cohort allele-frequency standardization avoids the cross-cohort
# allele-frequency fingerprint that drove the cell-level pooled track to
# cutoff 0.20. Reuses the per-cohort GRMs from step 03c (cutoff 0.05) and
# the per-(cohort × cell) keep files from step 05c. See Plan_deviations.md §7e.

# Minimum (cohort × cell) sub-sample N to include as a block component. 30
# includes all sub-cells with meaningful information; smaller sub-cells are
# skipped to avoid log clutter (their per-pair signal is negligible and
# inclusion doesn't change V_G estimation).
export MIN_BLOCKCOMP_N="30"

# --- Phenotype + covariate file paths ---------------------------------------
# Produced on the Mac, rsynced into ${PROJ_ROOT}/data/. GCTA file formats:
#   *.phen   = "FID IID y"          (--pheno)
#   *.qcovar = "FID IID q1 q2 ..."  (--qcovar; numeric covariates: 20 crossPCs)
#   *.covar  = "FID IID c1 c2 ..."  (--covar;  categorical covariates: gender)
export PHENO_FILE="${PROJ_ROOT}/data/education.phen"
export QCOVAR_FILE="${PROJ_ROOT}/data/crosspcs.qcovar"
export COVAR_FILE="${PROJ_ROOT}/data/gender.covar"

# Block-diagonal track only (step 05d). Adds a categorical cohort indicator
# alongside gender to absorb between-cohort mean differences in education
# that would otherwise inflate V_e under the block-diagonal design (because
# cross-cohort pairs contribute zero to V_G by construction). See
# Plan_deviations.md §7g.
export BLOCKDIAG_COVAR_FILE="${PROJ_ROOT}/data/gender_cohort.covar"

# --- Additional GREML comparison tracks (Plan_deviations.md §7h) -----------
# Secondary/exploratory tracks built on the same block-diagonal design as
# the cell-level primary, with track-specific strata and covariates. The
# new tracks use the *subset-then-prune* block extractor (03e), which
# subsets PRE-cutoff per-cohort GRMs to (cohort × new stratum) and then
# applies --grm-cutoff 0.05 within that sub-block. This avoids removing
# within-cohort relatives who fall into different comparison strata (e.g.
# cross-gender sib pairs in TwinLife).
#
# The cell-level primary track (03d/04d/05d/06d) uses the opposite order:
# it subsets POST-cutoff per-cohort GRMs, which is what its validated h²
# estimates are based on. See Plan_deviations.md §7h for the deliberate
# pruning-order difference between the two designs.

# R1: Region main comparison (East vs West, collapsed across pre/post-1990)
export STRATA_REGION_MAIN=(region_east region_west)
# G1: recorded gender/sex main comparison
export STRATA_GENDER_MAIN=(gender_male gender_female)
# RG1: Region × recorded gender/sex
export STRATA_REGION_GENDER=(region_gender_east_male
                             region_gender_west_male
                             region_gender_east_female
                             region_gender_west_female)

# Track-specific covar files (written by export_stratum_phenotype.R):
#   region_cohort.covar:  region_int + cohort_int (used by G1)
#   cohort_only.covar:    cohort_int only         (used by RG1)
#   gender_cohort.covar:  gender_int + cohort_int (used by cell-level primary AND R1)
export REGION_COHORT_COVAR_FILE="${PROJ_ROOT}/data/region_cohort.covar"
export COHORT_ONLY_COVAR_FILE="${PROJ_ROOT}/data/cohort_only.covar"
