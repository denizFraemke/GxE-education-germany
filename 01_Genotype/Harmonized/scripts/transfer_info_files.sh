#!/usr/bin/env bash
# =============================================================================
# transfer_info_files.sh — Transfer imputation info files to Tardis HPC
# =============================================================================
# Copies Minimac .info.gz files for BASEII / SHIP-0 / SHIP-Td / SOEP from
# the project share volumes to the cluster under data/info_files/{COHORT}/.
#
# TwinLife's R² TSV (which substitutes for a missing .info.gz) is NOT handled
# here — it rides along with `transfer_geno_to_tardis.sh` because both come
# off the same TwinLife share.
#
# PREREQUISITES:
#   Mount the relevant share volumes and point ${BASE2_SHARE}, ${SHIP_SHARE}
#   and ${SOEP_SHARE} at them.
#
# USAGE:
#   bash scripts/transfer_info_files.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT from user_config.sh (or run the
# interactive bootstrap on first invocation).
# shellcheck source=../_bootstrap_user_config.sh
source "${SCRIPT_DIR}/_bootstrap_user_config.sh"

TARDIS="${TARDIS_SSH}"
TARDIS_INFO_DIR="${TARDIS_PROJECT_ROOT}/data/info_files"
LOG="${SCRIPT_DIR}/transfer_info_files.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# =============================================================================
# SOURCE PATHS (local volumes)
# =============================================================================
# Mount points of the per-cohort project shares. Override any of these in the
# environment; the defaults are placeholders.
BASE2_SHARE="${BASE2_SHARE:-/path/to/base2_share}"
SHIP_SHARE="${SHIP_SHARE:-/path/to/ship_share}"
SOEP_SHARE="${SOEP_SHARE:-/path/to/soep_share}"

BASEII_INFO="${BASE2_SHARE}/private/data/001_BASEII_vcf/BASEII_AFFY6.0_imputed/LIFEBRAIN_BASEII_2015_AFFY.info.gz"
SHIP_INFO_DIR="${SHIP_SHARE}/private/data/003_SHIP_info_files"
SOEP_INFO_DIR="${SOEP_SHARE}/private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/info_files"

# =============================================================================
# CHECK SOURCES
# =============================================================================
log "Checking source files..."
missing=0

# BASEII
if [[ -f "$BASEII_INFO" ]]; then
  log "  OK: BASEII info file"
else
  log "  MISSING: $BASEII_INFO"
  log "    Mount the BASE-II share and set \${BASE2_SHARE} first"
  missing=1
fi

# SHIP-0
ship0_test="${SHIP_INFO_DIR}/SHIP-0_R4a.chr1.info.gz"
if [[ -f "$ship0_test" ]]; then
  log "  OK: SHIP-0 info files (22 chr)"
else
  log "  MISSING: $ship0_test"
  log "    Mount the SHIP share and set \${SHIP_SHARE} first"
  missing=1
fi

# SHIP-Td (batch 1 + batch 2)
shiptd_test="${SHIP_INFO_DIR}/SHIP-Td.chr1.info.gz"
shiptd_b2_test="${SHIP_INFO_DIR}/SHIP-Td_B2.chr1.info.gz"
if [[ -f "$shiptd_test" && -f "$shiptd_b2_test" ]]; then
  log "  OK: SHIP-Td info files (batch 1 + batch 2, 22 chr each)"
else
  log "  MISSING: SHIP-Td info files"
  log "    Mount the SHIP share and set \${SHIP_SHARE} first"
  missing=1
fi

# SOEP — the chr1 info file is not part of the provider delivery (see README
# §"Cohort notes"). Probe chr2 instead so the script still runs.
soep_test="${SOEP_INFO_DIR}/SOEP-G.b37.mildQC.hrc1-1_imp.chr2.info.gz"
if [[ -f "$soep_test" ]]; then
  log "  OK: SOEP info files (chr2..22; chr1 expected absent — provider-missing)"
else
  log "  MISSING: $soep_test"
  log "    Mount the SOEP share and set \${SOEP_SHARE} first"
  missing=1
fi

if [[ $missing -eq 1 ]]; then
  log "ERROR: Some files missing. Fix the issues above and re-run."
  exit 1
fi

# =============================================================================
# CREATE REMOTE DIRECTORIES
# =============================================================================
log ""
log "Creating remote directories on Tardis..."
ssh "$TARDIS" "mkdir -p ${TARDIS_INFO_DIR}/{BASEII,SHIP0,SHIPTD,SOEP}"

# =============================================================================
# DRY RUN PREVIEW
# =============================================================================
log ""
log "=== FILES TO TRANSFER ==="
log ""
log "BASEII:   1 file (single genome-wide info)"
log "SHIP0:    22 files (per-chr info)"
log "SHIPTD:   44 files (22 chr x 2 batches)"
log "SOEP:     21 files (chr2..22; chr1 not delivered — see README, Cohort notes)"
log ""
log "Total: ~88 files"
log "(TwinLife R² rides along with transfer_geno_to_tardis.sh.)"
log "Destination: ${TARDIS}:${TARDIS_INFO_DIR}/"
log ""

read -p "Proceed with transfer? (y/n) " confirm
if [[ "$confirm" != "y" ]]; then
  echo "Aborted."
  exit 0
fi

# =============================================================================
# TRANSFER
# =============================================================================
log ""
log "=== TRANSFERRING ==="
caffeinate -d -i -s -w $$ &
CAFFEINATE_PID=$!

# BASEII (single file)
log "Transferring BASEII..."
rsync -avP "$BASEII_INFO" "${TARDIS}:${TARDIS_INFO_DIR}/BASEII/" 2>&1 | tee -a "$LOG"

# SHIP-0 (22 chr)
log "Transferring SHIP-0..."
rsync -avP "${SHIP_INFO_DIR}"/SHIP-0_R4a.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.info.gz \
  "${TARDIS}:${TARDIS_INFO_DIR}/SHIP0/" 2>&1 | tee -a "$LOG"

# SHIP-Td batch 1 + batch 2 (44 files)
log "Transferring SHIP-Td (both batches)..."
rsync -avP "${SHIP_INFO_DIR}"/SHIP-Td.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.info.gz \
  "${TARDIS}:${TARDIS_INFO_DIR}/SHIPTD/" 2>&1 | tee -a "$LOG"
rsync -avP "${SHIP_INFO_DIR}"/SHIP-Td_B2.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.info.gz \
  "${TARDIS}:${TARDIS_INFO_DIR}/SHIPTD/" 2>&1 | tee -a "$LOG"

# SOEP (chr2..22; chr1 is not part of the provider delivery — see README
# §"Cohort notes"). Glob picks up whatever chromosomes are actually on the
# share, so the missing chr1 doesn't break rsync.
log "Transferring SOEP (chr2..22; chr1 expected absent)..."
shopt -s nullglob
soep_files=( "${SOEP_INFO_DIR}"/SOEP-G.b37.mildQC.hrc1-1_imp.chr*.info.gz )
shopt -u nullglob
if [[ ${#soep_files[@]} -eq 0 ]]; then
  log "  ERROR: no SOEP info files found at ${SOEP_INFO_DIR}/"
  exit 1
fi
log "  Found ${#soep_files[@]} SOEP info files."
rsync -avP "${soep_files[@]}" "${TARDIS}:${TARDIS_INFO_DIR}/SOEP/" 2>&1 | tee -a "$LOG"

kill "$CAFFEINATE_PID" 2>/dev/null || true

log ""
log "=== TRANSFER COMPLETE ==="
log "Info files are now at: ${TARDIS}:${TARDIS_INFO_DIR}/"
log ""
log "Next: upload updated pipeline scripts with:"
log "  bash transfer_pipeline_to_tardis.sh"
