#!/usr/bin/env bash
# =============================================================================
# transfer_results_from_tardis.sh — Pull final S5 TSVs back to the Mac
# =============================================================================
# rsync the contents of ${TARDIS_PROJECT_ROOT}/output/final/ on Tardis to
# ${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/output/final/ on the
# Mac (the shared-volume mirror). Run this AFTER the Tardis pipeline
# (steps 00–09) finishes and BEFORE running summarize_results.R.
#
# The bundle is small (TSVs + a markdown row); takes seconds.
#
# Usage:
#   bash 04_GxE_Analysis/02_SNP_Heritability/transfer_results_from_tardis.sh
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

SRC="${TARDIS_SSH}:${TARDIS_PROJECT_ROOT}/output/final/"
DEST="${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/output/final/"
LOG="${PIPELINE_DIR}/transfer_results.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "============================================================"
log "Pulling final S5 results from Tardis"
log "  Source: $SRC"
log "  Dest:   $DEST"
log "============================================================"

mkdir -p "$DEST"

# Pull only the deliverables (TSVs + markdown). Skip any stale large files.
rsync -avz \
    --include '*.tsv' \
    --include '*.md' \
    --include '*.json' \
    --exclude '*' \
    "$SRC" "$DEST" 2>&1 | tee -a "$LOG"

echo
log "Files now on Mac:"
ls -lh "$DEST" | tee -a "$LOG"

log "=== DONE ==="
log ""
log "Next: Rscript 04_GxE_Analysis/02_SNP_Heritability/summarize_results.R"
