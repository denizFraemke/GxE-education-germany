#!/usr/bin/env bash
# =============================================================================
# config.sh — Central configuration for PGI Computation Pipeline
# =============================================================================
# All paths, tool locations, trait definitions, and parameters.
# This file is sourced by every pipeline script.
#
# ENVIRONMENT: SLURM cluster.
# Code, data, output, and logs all live under ${PROJ_ROOT} (auto-detected).
# Per-user settings (cluster username, project root used by the Mac→cluster
# transfer scripts, bin/LD-ref locations, local analysis-project path) come
# from user_config.sh, which is gitignored and written by
# _bootstrap_user_config.sh on first run.
# =============================================================================

set -euo pipefail

# --- Project root (auto-detected from this file's location) -----------------
PROJ_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ_ROOT

# --- Per-user overrides (optional) ------------------------------------------
# user_config.sh is gitignored and may not exist on the cluster if a new user
# hasn't run any transfer_*.sh yet. When absent we fall through to the
# `${VAR:-default}` defaults below.
if [[ -f "${PROJ_ROOT}/user_config.sh" ]]; then
    # shellcheck source=/dev/null
    source "${PROJ_ROOT}/user_config.sh"
fi

# --- Output directories (created automatically) -----------------------------
export OUT_ROOT="${PROJ_ROOT}/output"
export SBAYESR_OUT="${OUT_ROOT}/sbayesr"
export WEIGHTS_OUT="${OUT_ROOT}/weights"
export GENO_OUT="${OUT_ROOT}/geno"
export PCS_OUT="${OUT_ROOT}/pcs"
export SCORES_OUT="${OUT_ROOT}/scores"
export FINAL_OUT="${OUT_ROOT}/final"
export LOG_DIR="${PROJ_ROOT}/logs"

# --- Tool binaries ----------------------------------------------------------
# GCTB (SBayesR), PLINK 1.9, PLINK 2.0. Anchored at ${TARDIS_BIN_DIR}
# (set in user_config.sh on the Mac, rsync'd to the cluster); falls back to
# ${CLUSTER_BIN_DIR}/bin if user_config.sh is absent.
_BIN_DIR="${TARDIS_BIN_DIR:-${CLUSTER_BIN_DIR:-/path/to/cluster_bin_dir}/bin}"
export GCTB="${_BIN_DIR}/gctb"
export PLINK19="${_BIN_DIR}/plink"
export PLINK2="${_BIN_DIR}/plink2"
unset _BIN_DIR

# --- LD reference for SBayesR -----------------------------------------------
# UKB 50k shrunk sparse LD matrices (Zeng et al.).
# File pattern: ${LD_REF_DIR}/ukb50k_shrunk_chr${CHR}_mafpt01.ldm.sparse
export LD_REF_DIR="${TARDIS_LD_REF_DIR:-${CLUSTER_BIN_DIR:-/path/to/cluster_bin_dir}/reference/ukb_50k_bigset_2.8M}"

# --- Formatted summary statistics for SBayesR -------------------------------
# Pre-formatted per-chromosome GWAS summary statistics.
# Copy or symlink these from the data share (or wherever you keep your
# formatted sumstats) to the cluster at the path below.
export SUMSTATS_DIR="${PROJ_ROOT}/data/formatted_sumstats"

# Height is the one trait whose pre-formatted share copy is malformed (extra
# CHR column + unfiltered NA rows); regenerate it from the RAW GIANT/Yengo 2022
# file with scripts/format_height_sumstats.sh. Point this at that raw file (on
# the data share: Projects/03_data/009_SUMSTATS/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz),
# or stage the raw file to the cluster and set the path here / in user_config.sh.
export GIANT_HEIGHT_RAW="${GIANT_HEIGHT_RAW:-}"

# --- Raw genotype data (bfile prefixes) --------------------------------------
# These are the imputed/QC'd genotype files for each cohort, expected under
# ${PROJ_ROOT}/data/geno/. Stage them there once with setup_tardis_data.sh.
export GENO_BASEII="${PROJ_ROOT}/data/geno/LIFEBRAIN_BASEII_2015_AFFY"
export GENO_SHIP0="${PROJ_ROOT}/data/geno/Final-SHIP-0_R4a.chrall.dose"
export GENO_SHIPTD="${PROJ_ROOT}/data/geno/SHIP-Td-SHIP-Td_B2_merged.chrall"
# NOTE: the SOEP genotypes are not available on the share, so these paths do
# not resolve on a fresh transfer. Two releases exist: the provider-merged file
# is named SOEP-G.b37.chrall (no "ext_"), while this path points at the larger
# extended share fileset — reconcile the prefix if SOEP is ever re-staged.
export GENO_SOEP="${PROJ_ROOT}/data/geno/SOEP-G.b37.ext_chrall"
export GENO_TWINLIFE="${PROJ_ROOT}/data/geno/TwinLife.b37.chrall"

# --- Cohort list -------------------------------------------------------------
export COHORTS=("BASEII" "SHIP0" "SHIPTD" "SOEP" "TWINLIFE")

# --- SBayesR default parameters ----------------------------------------------
export SBAYESR_GAMMA="0,0.01,0.1,1"
export SBAYESR_PI="0.95,0.02,0.02,0.01"
# Chain length 10,000 with 2,000 burn-in are the settings that produced the
# reported polygenic indices; changing them requires recomputing the PGIs and
# re-importing every downstream result. GCTB's MCMC is stochastic and no RNG
# seed is set, so the per-SNP posterior weights are not bit-reproducible. The
# effect on the score is negligible: the seed reshuffles effects among
# correlated SNPs within an LD block and the PGI is a weighted sum over those
# same dosages, so the aggregation absorbs it. Measured cross-seed PGI
# correlation r = 1.000 (5 seeds, TwinLife, N = 5,861).
export SBAYESR_CHAIN_LENGTH=10000
export SBAYESR_BURNIN=2000
export SBAYESR_OUT_FREQ=10
export SBAYESR_EXCLUDE_MHC="yes"  # "yes" or "no"

# --- PLINK2 scoring defaults -------------------------------------------------
export PLINK2_THREADS=4
export PLINK2_SET_VAR_IDS="chr@:#"  # chr:pos format
export PLINK2_RM_DUP="force-first"

# =============================================================================
# CROSS-COHORT HARMONIZATION SETTINGS
# =============================================================================
# These settings control the harmonization pipeline.
# Adjust based on diagnostic results.

# Fix strand/allele orientation issues for NonCog/Cog (step 03b)?
# DISABLED BY DESIGN — there are no non-palindromic complement mismatches in
# any cohort against any weight file, so there is nothing for an allele-based
# strand fix to do on this dataset. The only SNPs the current script would
# touch are palindromes (A/T, C/G), whose strand cannot be resolved from
# alleles alone; flipping them inverts the dosage interpretation (a +13
# raw-unit additive shift in the NonCog PGI for SOEP/TwinLife) without
# justification. Leave at "no" unless a future cohort introduces real
# non-palindromic complement mismatches AND the script is rewritten to
# (a) exclude palindromes and (b) loop data-drivenly across cohorts.
# See README §Cross-cohort harmonization #4 for full context.
export FIX_STRAND="no"

# Restrict all cohorts to a unified SNP set (step 03, phase 2)?
# Set to "yes" to compute the intersection of SNPs across all cohorts
# and restrict genotypes to that set before scoring. This eliminates
# platform-specific SNP content as a source of cross-cohort PGI differences.
export SHARED_SNP_SET="yes"

# Align genotypes to UKB reference panel (step 03c)?
# Set to "yes" to additionally align alleles to the UKB 50k European reference
# panel (which the SBayesR weights were derived from). Requires the UKB
# reference to be accessible.
export UKB_ALIGN="no"

# --- Imputation quality (R²) filtering ----------------------------------------
# Filter variants by imputation R² before SNP intersection.
# Info files must be in data/info_files/ (transfer with scripts/transfer_info_files.sh).
export R2_FILTER="${R2_FILTER:-yes}"
export R2_THRESHOLD="${R2_THRESHOLD:-0.3}"          # default (lenient)
export R2_THRESHOLD_STRICT="${R2_THRESHOLD_STRICT:-0.8}"  # sensitivity analysis

# Info file directory on Tardis (created by transfer_info_files.sh)
export INFO_DIR="${PROJ_ROOT}/data/info_files"

# Per-cohort info file locations (Minimac format: col 1=SNP, col 7=Rsq)
# SNP IDs in info files use bare chr:pos (e.g., "22:16050435"); the pipeline
# prepends "chr" to match harmonized bim IDs ("chr22:16050435").
export INFO_BASEII="${INFO_DIR}/BASEII/LIFEBRAIN_BASEII_2015_AFFY.info.gz"
export INFO_SHIP0_PATTERN="${INFO_DIR}/SHIP0/SHIP-0_R4a.chr{CHR}.info.gz"
export INFO_SHIPTD_PATTERN="${INFO_DIR}/SHIPTD/SHIP-Td.chr{CHR}.info.gz"
export INFO_SHIPTD_B2_PATTERN="${INFO_DIR}/SHIPTD/SHIP-Td_B2.chr{CHR}.info.gz"
export INFO_SOEP_PATTERN="${INFO_DIR}/SOEP/SOEP-G.b37.mildQC.hrc1-1_imp.chr{CHR}.info.gz"
export INFO_TWINLIFE="${INFO_DIR}/TWINLIFE/TwinLife_imputation_r2.tsv.gz"

# =============================================================================
# TRAIT DEFINITIONS
# =============================================================================
# Each trait specifies:
#   - SUMSTATS_SUBDIR: subdirectory under SUMSTATS_DIR
#   - CHR_FILE_PATTERN: per-chromosome file pattern ({CHR} is replaced)
#   - SBAYESR_TAG: output tag for SBayesR results
#   - WEIGHT_COLS: "id_col a1_col beta_col" for PLINK2 --score
#   - TARGET_COHORTS: which cohorts to score (space-separated)
#
# Weight file format after SBayesR (.snpRes):
#   Col 1: Index
#   Col 2: Name (rsID or chr:pos:ref:alt)
#   Col 3: Chrom
#   Col 4: Position
#   Col 5: A1
#   Col 6: A2
#   Col 7: A1Frq
#   Col 8: A1Effect  <-- the SBayesR posterior mean effect
#   ...
# After conversion (step 02), scoring files have 3 columns:
#   Col 1: SNP (chr:pos)  Col 2: A1  Col 3: BETA
# For PLINK2 scoring: --score file 1 2 3 header-read
# =============================================================================

# Number of traits (set after defining arrays)
declare -a TRAIT_NAMES
declare -A TRAIT_SUMSTATS_SUBDIR
declare -A TRAIT_CHR_FILE_PATTERN
declare -A TRAIT_SBAYESR_TAG

# --- EA4 leave-one-out weights ----------------------------------------------
# Both cohort-specific LOO GWAS are run through SBayesR; scoring uses
# EA4_excl_SHIP for every cohort (see get_ea4_trait_for_cohort below).
TRAIT_NAMES+=("EA4_excl_BASEII")
TRAIT_SUMSTATS_SUBDIR["EA4_excl_BASEII"]="EA4_excl_BASEII"
TRAIT_CHR_FILE_PATTERN["EA4_excl_BASEII"]="EA4_excl_BASEII_chr{CHR}_SBayesRformat.txt"
TRAIT_SBAYESR_TAG["EA4_excl_BASEII"]="EA4_excl_BASEII"

TRAIT_NAMES+=("EA4_excl_SHIP")
TRAIT_SUMSTATS_SUBDIR["EA4_excl_SHIP"]="EA4_excl_SHIP"
TRAIT_CHR_FILE_PATTERN["EA4_excl_SHIP"]="EA4_excl_SHIP_chr{CHR}_SBayesRformat.txt"
TRAIT_SBAYESR_TAG["EA4_excl_SHIP"]="EA4_excl_SHIP"

# --- Cognitive (Malanchini et al., 2024) — all cohorts -----------------------
TRAIT_NAMES+=("Cog")
TRAIT_SUMSTATS_SUBDIR["Cog"]="Cog_Malanchini"
TRAIT_CHR_FILE_PATTERN["Cog"]="Cog_ext_chr{CHR}_SBayesRformat.txt"
TRAIT_SBAYESR_TAG["Cog"]="Cog"

# --- Non-Cognitive general (Malanchini et al.) — all cohorts -----------------
TRAIT_NAMES+=("NCog")
TRAIT_SUMSTATS_SUBDIR["NCog"]="NCog_Malanchini"
TRAIT_CHR_FILE_PATTERN["NCog"]="NCog_ext_chr{CHR}_SBayesRformat.txt"
TRAIT_SBAYESR_TAG["NCog"]="NCog"

# --- Height (Yengo et al., 2022 / GIANT) — all cohorts ----------------------
TRAIT_NAMES+=("Height")
TRAIT_SUMSTATS_SUBDIR["Height"]="Height_GIANT"
TRAIT_CHR_FILE_PATTERN["Height"]="HEIGHT_SBayesRformat.chr{CHR}.txt"
TRAIT_SBAYESR_TAG["Height"]="Height"

export TRAIT_NAMES

# =============================================================================
# HELPER: get genotype bfile for a cohort
# =============================================================================
get_geno_bfile() {
    local cohort="$1"
    case "$cohort" in
        BASEII)   echo "$GENO_BASEII" ;;
        SHIP0)    echo "$GENO_SHIP0" ;;
        SHIPTD)   echo "$GENO_SHIPTD" ;;
        SOEP)     echo "$GENO_SOEP" ;;
        TWINLIFE) echo "$GENO_TWINLIFE" ;;
        *) echo "ERROR: unknown cohort $cohort" >&2; return 1 ;;
    esac
}
export -f get_geno_bfile

# =============================================================================
# HELPER: get harmonized genotype prefix for a cohort
# =============================================================================
get_harmonized_bfile() {
    local cohort="$1"
    local shared_prefix="${GENO_OUT}/${cohort}_harmonized_shared"
    local base_prefix="${GENO_OUT}/${cohort}_harmonized"

    if [[ "${SHARED_SNP_SET:-yes}" == "yes" && -f "${shared_prefix}.bed" ]]; then
        echo "${shared_prefix}"
    else
        echo "${base_prefix}"
    fi
}
export -f get_harmonized_bfile

# =============================================================================
# HELPER: get the correct EA4 trait name for a cohort (LOO selection)
# =============================================================================
# All five cohorts use EA4_excl_SHIP. Rationale:
#   - BASEII contributes ~0.07% of the EA4 GWAS sample, so the in-sample bias
#     from scoring BASEII against EA4_excl_SHIP is statistically undetectable
#     (r² inflation ~ 1.0007).
#   - Using a single weight file across all five cohorts puts every PGI on the
#     same scale without requiring an analytical or empirical calibration
#     step. Cross-cohort raw-mean differences then reflect genuine population
#     variation (e.g. BASEII as the oldest cohort), not LOO methodology.
get_ea4_trait_for_cohort() {
    local cohort="$1"
    case "$cohort" in
        BASEII)   echo "EA4_excl_SHIP" ;;
        SHIP0)    echo "EA4_excl_SHIP" ;;
        SHIPTD)   echo "EA4_excl_SHIP" ;;
        SOEP)     echo "EA4_excl_SHIP" ;;
        TWINLIFE) echo "EA4_excl_SHIP" ;;
        *) echo "ERROR: unknown cohort $cohort" >&2; return 1 ;;
    esac
}
export -f get_ea4_trait_for_cohort

# =============================================================================
# HELPER: get the correct NonCog trait name for a cohort
# =============================================================================
# All five cohorts use the Malanchini PUBLISHED NonCog (NCog), whose LDSC
# genetic correlation with the Cog GWAS used here is rg = +0.16.
#
# NCog carries a different estimand scale (Σβ +3.40) and shifts the
# cross-cohort raw mean. That is a SCALING property, NOT a validity problem,
# and is neutralised by within-cohort standardisation (PGI_NonCog_within_z,
# computed in 03_Merge/merge.R). DOWNSTREAM: use the within-cohort z (or
# global z) for NonCog, never the raw cross-cohort mean.
get_noncog_trait_for_cohort() {
    local cohort="$1"
    case "$cohort" in
        BASEII|SHIP0|SHIPTD|SOEP|TWINLIFE) echo "NCog" ;;
        *) echo "ERROR: unknown cohort $cohort" >&2; return 1 ;;
    esac
}
export -f get_noncog_trait_for_cohort

# =============================================================================
# LOGGING
# =============================================================================
log_msg() {
    local level="$1"; shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $*"
}
export -f log_msg

log_info()  { log_msg "INFO"  "$@"; }
log_warn()  { log_msg "WARN"  "$@"; }
log_error() { log_msg "ERROR" "$@"; }
export -f log_info log_warn log_error

# =============================================================================
# Create output directories
# =============================================================================
ensure_dirs() {
    mkdir -p "$SBAYESR_OUT" "$WEIGHTS_OUT" "$GENO_OUT" "$PCS_OUT" \
             "$SCORES_OUT" "$FINAL_OUT" "$LOG_DIR"
}
export -f ensure_dirs
