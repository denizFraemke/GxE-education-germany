#!/usr/bin/env bash
#SBATCH --job-name=pgi_pipeline
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pgi_pipeline_%j.out
#SBATCH --error=logs/pgi_pipeline_%j.err
# =============================================================================
# submit_post_sbayesr.slurm — Run steps 02-07 after SBayesR jobs complete
# =============================================================================
# USAGE:
#   # After all SBayesR jobs finish:
#   sbatch submit_post_sbayesr.slurm
#
#   # Or with dependency on SBayesR jobs:
#   sbatch --dependency=afterok:<job1>:<job2>:... submit_post_sbayesr.slurm
# =============================================================================

set -euo pipefail

# Use SLURM CPU count if available
export PLINK2_THREADS=${SLURM_CPUS_PER_TASK:-8}

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
bash run_pipeline.sh --from 02 --to 07
