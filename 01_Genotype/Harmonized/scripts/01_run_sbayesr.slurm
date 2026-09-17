#!/usr/bin/env bash
#SBATCH --job-name=sbayesr
#SBATCH --array=1-22
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=6:00:00
#SBATCH --output=logs/sbayesr_%x_chr%a.%j.out
#SBATCH --error=logs/sbayesr_%x_chr%a.%j.err
# =============================================================================
# 01_run_sbayesr.slurm — Run SBayesR for one trait, one chromosome
# =============================================================================
# USAGE (submitted by 01_submit_all_sbayesr.sh):
#   sbatch --job-name=sbayesr_${TRAIT} \
#          --export=TRAIT=...,GAMMA=...,PI=...,CHAIN=...,BURNIN=...,EXCLUDE_MHC=...,OUT_TAG=... \
#          scripts/01_run_sbayesr.slurm
#
# INPUTS (with CHR = $SLURM_ARRAY_TASK_ID):
#   ${SUMSTATS_DIR}/<trait subdir>/<per-chr .ma file> — SBayesR-format sumstats
#   ${LD_REF_DIR}/ukb50k_shrunk_chr${CHR}_mafpt01.ldm.sparse — UKB 50k LD matrix
#
# OUTPUTS:
#   ${SBAYESR_OUT}/${OUT_TAG}/${OUT_TAG}_chr${CHR}_sbayesR.{snpRes,parRes,mcmcsamples}
#
# DEPENDS ON: step 00 (input checks). Submitted by 01_submit_all_sbayesr.sh.
# =============================================================================
set -euo pipefail

# Under SLURM, $0 points to a temporary spool copy of this script. Use the
# original submission directory so config and project-relative paths resolve.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

CHR=${SLURM_ARRAY_TASK_ID}

# These are passed via --export from the submission script
: "${TRAIT:?TRAIT not set}"
: "${OUT_TAG:?OUT_TAG not set}"
: "${GAMMA:=$SBAYESR_GAMMA}"
: "${PI:=$SBAYESR_PI}"
: "${CHAIN:=$SBAYESR_CHAIN_LENGTH}"
: "${BURNIN:=$SBAYESR_BURNIN}"
: "${EXCLUDE_MHC:=$SBAYESR_EXCLUDE_MHC}"

# Resolve input file
SUBDIR="${SUMSTATS_DIR}/${TRAIT_SUMSTATS_SUBDIR[$TRAIT]}"
PATTERN="${TRAIT_CHR_FILE_PATTERN[$TRAIT]}"
SUMSTAT_FILE="${SUBDIR}/${PATTERN//\{CHR\}/$CHR}"

if [[ ! -f "$SUMSTAT_FILE" ]]; then
    echo "ERROR: Summary statistics file not found: $SUMSTAT_FILE"
    exit 1
fi

# Output directory for this trait
TRAIT_OUT_DIR="${SBAYESR_OUT}/${OUT_TAG}"
mkdir -p "$TRAIT_OUT_DIR"

OUT_PREFIX="${TRAIT_OUT_DIR}/${OUT_TAG}_chr${CHR}_sbayesR"

# Skip if output already exists
if [[ -f "${OUT_PREFIX}.snpRes" ]]; then
    echo "SKIP: Output already exists: ${OUT_PREFIX}.snpRes"
    exit 0
fi

# LD matrix
LDM="${LD_REF_DIR}/ukb50k_shrunk_chr${CHR}_mafpt01.ldm.sparse"
if [[ ! -f "${LDM}.info" ]] && [[ ! -f "$LDM" ]]; then
    echo "ERROR: LD matrix not found: $LDM"
    exit 1
fi

# Build GCTB command
CMD=("$GCTB" --sbayes R
    --ldm "$LDM"
    --gwas-summary "$SUMSTAT_FILE"
    --gamma "$GAMMA"
    --pi "$PI"
    --chain-length "$CHAIN"
    --burn-in "$BURNIN"
    --out-freq "$SBAYESR_OUT_FREQ"
    --out "$OUT_PREFIX"
)

# Conditionally exclude MHC
if [[ "$EXCLUDE_MHC" == "yes" ]]; then
    CMD+=(--exclude-mhc)
fi

echo "================================================================"
echo "SBayesR: ${OUT_TAG} chromosome ${CHR}"
echo "Input:   $SUMSTAT_FILE"
echo "LD:      $LDM"
echo "Output:  $OUT_PREFIX"
echo "Params:  gamma=$GAMMA pi=$PI chain=$CHAIN burnin=$BURNIN mhc_excl=$EXCLUDE_MHC"
echo "================================================================"

"${CMD[@]}"

echo "SBayesR complete for ${OUT_TAG} chr${CHR}"
