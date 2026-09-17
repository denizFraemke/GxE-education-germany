#!/usr/bin/env bash
# =============================================================================
# transfer_final_from_tardis.sh — Pull genotype deliverables back from Tardis
# =============================================================================
# Phase 3 of the pipeline.
# Pulls ${TARDIS_PROJECT_ROOT}/output/final/ (and the diagnostics tree) from
# the cluster to the analysis project at ${ANALYSIS_ROOT}/data/final/, where the
# 03_Merge step consumes them:
#   - {COHORT}_PGI_PCs.tsv      — merge.R reads these
#   - {COHORT}_PC_eigenval.tsv  — validate.R reads these for the ancestry-PC
#                                 scree (variance explained per PC)
#   - PGI_summary_statistics.tsv — summary table
# The PGI tables and the eigenvalue tables both live in output/final/, so this
# single rsync pulls them together into the same folder.
#
# USAGE:
#   bash 01_Genotype/Harmonized/transfer_final_from_tardis.sh
#   (run on the Mac, SMB share mounted, after the cluster run finished)
# =============================================================================
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Loads TARDIS_SSH, TARDIS_PROJECT_ROOT, ANALYSIS_ROOT (writes user_config.sh
# interactively on first run).
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

: "${TARDIS_SSH:?set TARDIS_SSH in user_config.sh}"
: "${TARDIS_PROJECT_ROOT:?set TARDIS_PROJECT_ROOT in user_config.sh}"
: "${ANALYSIS_ROOT:?set ANALYSIS_ROOT in user_config.sh}"

SRC="${TARDIS_SSH}:${TARDIS_PROJECT_ROOT}/output"
DEST="${ANALYSIS_ROOT}/data/final"
LOG="${PIPELINE_DIR}/transfer_final.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "============================================================"
log "Pulling genotype deliverables from Tardis"
log "  Source: ${SRC}/final/  (+ diagnostics/)"
log "  Dest:   ${DEST}/"
log "============================================================"
mkdir -p "${DEST}" "${DEST}/diagnostics"

# Deliverables (TSV/MD/JSON only — the full output/ tree is large + regenerable).
rsync -avz \
    --include '*.tsv' --include '*.md' --include '*.json' --exclude '*' \
    "${SRC}/final/" "${DEST}/" 2>&1 | tee -a "$LOG"

# Diagnostics tree (recurse into subdirs, keep TSV/MD/JSON).
rsync -avz \
    --include '*/' --include '*.tsv' --include '*.md' --include '*.json' --exclude '*' \
    "${SRC}/diagnostics/" "${DEST}/diagnostics/" 2>&1 | tee -a "$LOG"

log ""
log "PC eigenvalue tables pulled (for the validate.R scree):"
if ls "${DEST}"/*_PC_eigenval.tsv >/dev/null 2>&1; then
    ls -1 "${DEST}"/*_PC_eigenval.tsv | tee -a "$LOG"
else
    log "  NONE found — re-run step 04 on Tardis (it emits them), then re-run this."
fi
log "=== DONE ==="
log "Next: Rscript 03_Merge/merge.R   then   Rscript 03_Merge/validate.R"
