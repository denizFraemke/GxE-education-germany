#!/usr/bin/env bash
# =============================================================================
# transfer_phenotype_to_tardis.sh — Sync phenotype + stratum keep-lists (Mac side)
# =============================================================================
# Pushes the files produced by export_stratum_phenotype.R from the shared
# volume to ${TARDIS_PROJECT_ROOT}/data/ on Tardis.
#
# Run order:
#   1) Rscript export_stratum_phenotype.R              # writes to ${DATA_ROOT}
#   2) bash    transfer_phenotype_to_tardis.sh         # uploads to Tardis
#   3) bash    transfer_pipeline_to_tardis.sh          # uploads code (or first-time)
#   4) ssh into Tardis, run `bash run_pipeline.sh`
#
# Usage:
#   bash 04_GxE_Analysis/02_SNP_Heritability/transfer_phenotype_to_tardis.sh
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# shellcheck source=_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

SRC_DIR="${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/data"
DEST="${TARDIS_SSH}:${TARDIS_PROJECT_ROOT}/data/"
LOG="${PIPELINE_DIR}/transfer_phenotype.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "============================================================"
log "Phenotype source: $SRC_DIR"
log "Tardis dest:      $DEST"
log "============================================================"

if [[ ! -d "$SRC_DIR" ]]; then
    log "ERROR: $SRC_DIR does not exist."
    log "Run Rscript ${PIPELINE_DIR}/export_stratum_phenotype.R first."
    exit 1
fi

for required in education.phen crosspcs.qcovar gender.covar manifest.tsv; do
    if [[ ! -f "${SRC_DIR}/${required}" ]]; then
        log "ERROR: missing ${SRC_DIR}/${required}"
        log "Re-run export_stratum_phenotype.R."
        exit 1
    fi
done

if [[ ! -d "${SRC_DIR}/keep" ]] || [[ -z "$(ls -A "${SRC_DIR}/keep" 2>/dev/null)" ]]; then
    log "ERROR: ${SRC_DIR}/keep/ is empty."
    log "Re-run export_stratum_phenotype.R."
    exit 1
fi

ssh "${TARDIS_SSH}" "mkdir -p '${TARDIS_PROJECT_ROOT}/data'"

echo ""
log "=== DRY RUN (preview only — nothing transferred yet) ==="
rsync -avzn --delete \
    --exclude '.DS_Store' \
    "${SRC_DIR}/" "${DEST}"

echo ""
read -r -p "Proceed with the upload above? (y/n) " confirm
if [[ "$confirm" != "y" ]]; then
    log "Aborted by user."
    exit 0
fi

echo ""
log "=== UPLOADING ==="

rsync -avz --delete \
    --exclude '.DS_Store' \
    "${SRC_DIR}/" "${DEST}" \
    2>&1 | tee -a "$LOG"

log "=== DONE ==="
log ""
log "Tardis now has: ${TARDIS_PROJECT_ROOT}/data/{education.phen,crosspcs.qcovar,gender.covar,manifest.tsv,keep/*.iids}"
