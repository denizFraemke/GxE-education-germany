#!/usr/bin/env bash
# =============================================================================
# transfer_pipeline_to_tardis.sh — Sync the Harmonized pipeline code to Tardis
# =============================================================================
# Pushes the monorepo's `01_Genotype/Harmonized/` to ${TARDIS_PROJECT_ROOT} on
# Tardis. Excludes everything that's gitignored or runtime-only (output/,
# logs/, data/, *.bak*, decision_log.md). caffeinate keeps the Mac awake
# during the rsync. user_config.sh IS rsynced — that's how your Tardis copy
# of the pipeline picks up your settings.
#
# Usage:
#   bash 01_Genotype/Harmonized/transfer_pipeline_to_tardis.sh
# =============================================================================
set -euo pipefail

# Resolve PIPELINE_DIR from this script's own location so the script is
# location-agnostic — works whether you're in the worktree, a normal clone,
# or a temp checkout.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT, ... from user_config.sh (or walk
# the user through interactive setup on first run).
# shellcheck source=_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

TARDIS="${TARDIS_SSH}"
TARDIS_DEST="${TARDIS_PROJECT_ROOT}"
LOG="${PIPELINE_DIR}/transfer_pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# Sanity: confirm we're rsyncing the right thing
log "============================================================"
log "Pipeline source:   $PIPELINE_DIR"
log "Tardis destination: ${TARDIS}:${TARDIS_DEST}/"
log "============================================================"

if [[ ! -f "${PIPELINE_DIR}/config.sh" ]] || [[ ! -f "${PIPELINE_DIR}/run_pipeline.sh" ]]; then
  log "ERROR: PIPELINE_DIR does not look like the Harmonized pipeline"
  log "       (missing config.sh or run_pipeline.sh under $PIPELINE_DIR)"
  exit 1
fi

# Common rsync excludes — keep code, drop runtime artifacts.
RSYNC_EXCLUDES=(
  --exclude 'output/'
  --exclude 'logs/'
  --exclude 'data/'
  --exclude '*.log'
  --exclude '*.bak*'
  --exclude '.DS_Store'
  --exclude '_attic_*/'
  --exclude '_archive_*/'
  --exclude 'transfer_geno.completed'
  --exclude 'transfer_geno_workbench.completed'
  --exclude 'decision_log.md'
)

# Dry-run preview first so the user can eyeball what's about to transfer.
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

# Real upload.
echo ""
log "=== UPLOADING ==="

# Prevent the Mac from sleeping during the transfer (no-op on Linux).
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
