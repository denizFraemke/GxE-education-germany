#!/usr/bin/env bash
# =============================================================================
# extract_and_transfer_info_workbench.sh
# =============================================================================
# Run this on the R workbench / ARC.
# Transfers .info.gz files for BASEII / SHIP-0 / SHIP-Td / SOEP directly
# from the project shares to the cluster. Faster than the Mac/SMB route.
#
# TwinLife's R² TSV (which substitutes for a missing .info.gz) is NOT handled
# here — it rides along with `transfer_geno_to_tardis_workbench.sh` because
# both come off the same TwinLife share.
#
# USAGE (on workbench):
#   bash extract_and_transfer_info_workbench.sh
#
# PREREQUISITES:
#   - SSH access to the cluster (key-based auth; see user_config.sh)
#   - project shares mounted at ${SHARE_MOUNT_ROOT}/
# =============================================================================
set -euo pipefail

# --- Configuration -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT from user_config.sh (or run the
# interactive bootstrap on first invocation).
# shellcheck source=../_bootstrap_user_config.sh
source "${SCRIPT_DIR}/_bootstrap_user_config.sh"

# Mount root of the project shares. Override in the environment; the default
# is a placeholder.
SRT_BASE="${SHARE_MOUNT_ROOT:-/path/to/share_mount_root}"
TARDIS="${TARDIS_SSH}"
TARDIS_INFO_DIR="${TARDIS_PROJECT_ROOT}/data/info_files"
TMP_DIR="/tmp/pgi_info_files_$$"

# Source paths on workbench
BASEII_INFO="${SRT_BASE}/<base2-share>/private/data/001_BASEII_vcf/BASEII_AFFY6.0_imputed/LIFEBRAIN_BASEII_2015_AFFY.info.gz"
SHIP_INFO_DIR="${SRT_BASE}/<ship-share>/private/data/003_SHIP_info_files"
SOEP_INFO_DIR="${SRT_BASE}/<soep-share>/private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/info_files"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# =============================================================================
# 1. CHECK ALL SOURCE PATHS
# =============================================================================
log "Checking source paths on workbench..."
error=0

for label_path in \
  "BASEII info:$(dirname "$BASEII_INFO")" \
  "SHIP info:${SHIP_INFO_DIR}" \
  "SOEP info:${SOEP_INFO_DIR}"; do
  label="${label_path%%:*}"
  path="${label_path#*:}"
  if [[ -d "$path" ]]; then
    log "  OK: $label ($path)"
  else
    log "  MISSING: $label ($path)"
    error=1
  fi
done

if [[ $error -eq 1 ]]; then
  log "ERROR: Some shares not accessible. Check the \${SHARE_MOUNT_ROOT} mounts."
  exit 1
fi

# =============================================================================
# 2. STAGE ALL INFO FILES INTO TMP DIRECTORY
# =============================================================================
log ""
log "=== Staging info files ==="
mkdir -p "$TMP_DIR"

# BASEII (single file)
mkdir -p "${TMP_DIR}/BASEII"
cp "$BASEII_INFO" "${TMP_DIR}/BASEII/"
log "  BASEII: 1 file copied"

# SHIP-0 (per-chr)
mkdir -p "${TMP_DIR}/SHIP0"
n=0
for CHR in $(seq 1 22); do
  src="${SHIP_INFO_DIR}/SHIP-0_R4a.chr${CHR}.info.gz"
  if [[ -f "$src" ]]; then
    cp "$src" "${TMP_DIR}/SHIP0/"
    n=$((n + 1))
  else
    log "  WARNING: Missing SHIP-0 chr${CHR}"
  fi
done
log "  SHIP-0: $n files copied"

# SHIP-Td (batch 1 + batch 2)
mkdir -p "${TMP_DIR}/SHIPTD"
n=0
for CHR in $(seq 1 22); do
  for prefix in "SHIP-Td" "SHIP-Td_B2"; do
    src="${SHIP_INFO_DIR}/${prefix}.chr${CHR}.info.gz"
    if [[ -f "$src" ]]; then
      cp "$src" "${TMP_DIR}/SHIPTD/"
      n=$((n + 1))
    else
      log "  WARNING: Missing ${prefix} chr${CHR}"
    fi
  done
done
log "  SHIP-Td: $n files copied (batch 1 + batch 2)"

# SOEP (per-chr)
mkdir -p "${TMP_DIR}/SOEP"
n=0
for CHR in $(seq 1 22); do
  src="${SOEP_INFO_DIR}/SOEP-G.b37.mildQC.hrc1-1_imp.chr${CHR}.info.gz"
  if [[ -f "$src" ]]; then
    cp "$src" "${TMP_DIR}/SOEP/"
    n=$((n + 1))
  else
    log "  WARNING: Missing SOEP chr${CHR}"
  fi
done
log "  SOEP: $n files copied"

# =============================================================================
# 3. TRANSFER TO TARDIS
# =============================================================================
log ""
log "=== Transferring to Tardis ==="
log "  Destination: ${TARDIS}:${TARDIS_INFO_DIR}/"

# Create remote directories
ssh "$TARDIS" "mkdir -p ${TARDIS_INFO_DIR}/{BASEII,SHIP0,SHIPTD,SOEP}"

# Transfer all at once
rsync -avP "${TMP_DIR}/" "${TARDIS}:${TARDIS_INFO_DIR}/"

# =============================================================================
# 4. CLEANUP
# =============================================================================
rm -rf "$TMP_DIR"

log ""
log "=== DONE ==="
log "Info files are now at: ${TARDIS}:${TARDIS_INFO_DIR}/"
log ""
log "Contents:"
ssh "$TARDIS" "for d in ${TARDIS_INFO_DIR}/*/; do echo \"  \$(basename \$d): \$(ls \$d | wc -l) files\"; done"
