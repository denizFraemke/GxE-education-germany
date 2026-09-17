#!/usr/bin/env bash
# =============================================================================
# transfer_pipeline_to_tardis.sh — Sync this section's CODE to Tardis (Mac side)
# =============================================================================
# Pushes the monorepo's `04_GxE_Analysis/02_SNP_Heritability/` to
# ${TARDIS_PROJECT_ROOT} on Tardis. Excludes runtime artefacts (data/,
# output/, logs/, *.bak*, *.log). user_config.sh IS rsynced — that's how the
# Tardis side learns its paths and TARDIS_GENOTYPE_PROJECT_ROOT.
#
# Run this from the Mac. After it lands, SSH into Tardis and:
#   cd ${TARDIS_PROJECT_ROOT}
#   bash run_pipeline.sh
#
# Usage:
#   bash 04_GxE_Analysis/02_SNP_Heritability/transfer_pipeline_to_tardis.sh
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT, ... from user_config.sh (interactive
# setup on first run).
# shellcheck source=_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

TARDIS="${TARDIS_SSH}"
TARDIS_DEST="${TARDIS_PROJECT_ROOT}"
LOG="${PIPELINE_DIR}/transfer_pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "============================================================"
log "Pipeline source:    $PIPELINE_DIR"
log "Tardis destination: ${TARDIS}:${TARDIS_DEST}/"
log "============================================================"

if [[ ! -f "${PIPELINE_DIR}/config.sh" ]] || [[ ! -f "${PIPELINE_DIR}/run_pipeline.sh" ]]; then
    log "ERROR: PIPELINE_DIR does not look like the SNP-h² section"
    log "       (missing config.sh or run_pipeline.sh under $PIPELINE_DIR)"
    exit 1
fi

# Ensure the destination root exists on Tardis (--mkpath needs newer rsync; do
# it via ssh for portability).
ssh "$TARDIS" "mkdir -p '${TARDIS_DEST}'"

RSYNC_EXCLUDES=(
    --exclude 'data/'
    --exclude 'output/'
    --exclude 'logs/'
    --exclude '*.log'
    --exclude '*.bak*'
    --exclude '.DS_Store'
)

echo ""
log "=== DRY RUN (preview only — nothing transferred yet) ==="
rsync -avzn --delete \
    "${RSYNC_EXCLUDES[@]}" \
    "${PIPELINE_DIR}/" \
    "${TARDIS}:${TARDIS_DEST}/"

echo ""
read -r -p "Proceed with the upload above? (y/n) " confirm
if [[ "$confirm" != "y" ]]; then
    log "Aborted by user."
    exit 0
fi

echo ""
log "=== UPLOADING ==="

if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -d -i -s -w $$ &
    CAFFEINATE_PID=$!
    log "Sleep prevention active (caffeinate PID: $CAFFEINATE_PID)"
fi

rsync -avz --delete \
    "${RSYNC_EXCLUDES[@]}" \
    "${PIPELINE_DIR}/" \
    "${TARDIS}:${TARDIS_DEST}/" \
    2>&1 | tee -a "$LOG"

if command -v caffeinate >/dev/null 2>&1 && [[ -n "${CAFFEINATE_PID:-}" ]]; then
    kill "$CAFFEINATE_PID" 2>/dev/null || true
fi

log "=== DONE ==="
log ""
log "Next: SSH into Tardis and run the pipeline:"
log "  ssh ${TARDIS}"
log "  cd ${TARDIS_DEST}"
log "  bash run_pipeline.sh"
