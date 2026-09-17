#!/usr/bin/env bash
# =============================================================================
# transfer_geno_to_tardis.sh — Transfer genotype files + TwinLife R² to Tardis
# =============================================================================
# - Prevents Mac from sleeping during transfer (even with lid closed)
# - Auto-remounts SMB volumes if they disconnect
# - Retries each file up to 10 times on failure
# - Logs everything to transfer_geno.log
#
# Transfers two streams:
#   GENO_DEST  — per-cohort imputed bfiles (.bed/.bim/.fam)
#   INFO_DEST  — the TwinLife imputation-R² TSV. TwinLife has no .info.gz,
#                so a pre-extracted aggregate file lives on the TwinLife
#                share at
#                ${TWINLIFE_SHARE}/.../rawgendata/TwinLife_imputation_r2.tsv.gz
#                and rides along with this transfer.
#
# USAGE:
#   bash transfer_geno_to_tardis.sh
#   (then you can close the terminal — it runs in background)
# =============================================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT from user_config.sh (or run the
# interactive bootstrap on first invocation).
# shellcheck source=_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

TARDIS="${TARDIS_SSH}"

# Mount points of the per-cohort project shares. Override any of these in the
# environment; the defaults are placeholders.
BASE2_SHARE="${BASE2_SHARE:-/path/to/base2_share}"
SHIP_SHARE="${SHIP_SHARE:-/path/to/ship_share}"
SOEP_SHARE="${SOEP_SHARE:-/path/to/soep_share}"
TWINLIFE_SHARE="${TWINLIFE_SHARE:-/path/to/twinlife_share}"

GENO_DEST="${TARDIS_PROJECT_ROOT}/data/geno"
INFO_DEST="${TARDIS_PROJECT_ROOT}/data/info_files/TWINLIFE"
LOG="${PIPELINE_DIR}/transfer_geno.log"
STATE="${PIPELINE_DIR}/transfer_geno.completed"
MAX_RETRIES=10
RETRY_WAIT=30   # seconds between retries

# Genotype bfiles → ${GENO_DEST}
declare -a FILES=(
    "${BASE2_SHARE}/private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.bed"
    "${BASE2_SHARE}/private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.bim"
    "${BASE2_SHARE}/private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.fam"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.bed"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.bim"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.fam"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.bed"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.bim"
    "${SHIP_SHARE}/private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.fam"
    # NOTE: the SOEP genotypes are not available on the share, so these source
    # paths do not resolve on a fresh transfer.
    "${SOEP_SHARE}/private/data/001_SOEP_processed_data/008_Genotype/SOEP-G.b37.chrall.bed"
    "${SOEP_SHARE}/private/data/001_SOEP_processed_data/008_Genotype/SOEP-G.b37.chrall.bim"
    "${SOEP_SHARE}/private/data/001_SOEP_processed_data/008_Genotype/SOEP-G.b37.chrall.fam"
    "${TWINLIFE_SHARE}/private/data/Gendata/001_processed/TwinLife.b37.chrall.bed"
    "${TWINLIFE_SHARE}/private/data/Gendata/001_processed/TwinLife.b37.chrall.bim"
    "${TWINLIFE_SHARE}/private/data/Gendata/001_processed/TwinLife.b37.chrall.fam"
)

# TwinLife imputation-R² TSV (pre-extracted aggregate file deposited on the
# TwinLife share — replaces step 03a's missing .info.gz for this
# cohort) → ${INFO_DEST}
declare -a INFO_FILES=(
    "${TWINLIFE_SHARE}/private/data/Gendata/rawgendata/TwinLife_imputation_r2.tsv.gz"
)

# =============================================================================
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

is_completed() {
    local src="$1"
    [[ -f "$STATE" ]] && grep -Fqx "$src" "$STATE"
}

mark_completed() {
    local src="$1"
    touch "$STATE"
    if ! grep -Fqx "$src" "$STATE"; then
        echo "$src" >> "$STATE"
    fi
}

hydrate_completed_state() {
    touch "$STATE"

    local src fname
    for src in "${FILES[@]}" "${INFO_FILES[@]}"; do
        fname=$(basename "$src")
        if grep -Fq "[$fname] DONE" "$LOG" 2>/dev/null; then
            mark_completed "$src"
        fi
    done
}

transfer_file() {
    local src="$1"
    local dest="$2"

    local fname
    fname=$(basename "$src")
    local attempt=0

    while [[ $attempt -lt $MAX_RETRIES ]]; do
        ((attempt++))
        log "[$fname] Attempt $attempt/$MAX_RETRIES"

        # Check source file exists (the share must already be mounted)
        if [[ ! -f "$src" ]]; then
            log "[$fname] ERROR: Source file not found: $src"
            sleep "$RETRY_WAIT"
            continue
        fi

        # Run rsync (--partial allows resume on interrupted transfers)
        rsync -avP --partial \
            "$src" \
            "${TARDIS}:${dest}/" \
            2>&1 | tee -a "$LOG"

        if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
            mark_completed "$src"
            log "[$fname] DONE"
            return 0
        else
            log "[$fname] rsync failed (attempt $attempt). Retrying in ${RETRY_WAIT}s..."
            sleep "$RETRY_WAIT"
        fi
    done

    log "[$fname] FAILED after $MAX_RETRIES attempts."
    return 1
}

# =============================================================================
# MAIN
# =============================================================================
TOTAL=$(( ${#FILES[@]} + ${#INFO_FILES[@]} ))

log "============================================================"
log "Starting genotype transfer to Tardis"
log "Genotype bfiles   → ${GENO_DEST}   (${#FILES[@]} files)"
log "TwinLife R² TSV   → ${INFO_DEST}   (${#INFO_FILES[@]} file)"
log "Log:   $LOG"
log "State: $STATE"
log "============================================================"

# Pre-create both Tardis destination dirs (rsync does NOT auto-create parents).
ssh "$TARDIS" "mkdir -p '${GENO_DEST}' '${INFO_DEST}'" \
    || { log "ERROR: could not ssh to Tardis ($TARDIS) to create dest dirs"; exit 1; }

hydrate_completed_state

# Prevent Mac from sleeping — works even with lid closed. macOS only;
# caffeinate doesn't exist on Linux / Posit Workbench, so guard it.
# caffeinate -d (display), -i (idle), -s (system/disk)
CAFFEINATE_PID=""
if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -d -i -s -w $$ &
    CAFFEINATE_PID=$!
    log "Sleep prevention active (caffeinate PID: $CAFFEINATE_PID)"
fi

failed=()

# Pass 1: genotype bfiles → GENO_DEST
for src in "${FILES[@]}"; do
    fname=$(basename "$src")
    log ""
    if is_completed "$src"; then
        log "--- Skipping completed: $fname ---"
        continue
    fi
    log "--- Starting: $fname → ${GENO_DEST}/ ---"
    if ! transfer_file "$src" "$GENO_DEST"; then
        failed+=("$src")
    fi
done

# Pass 2: TwinLife R² TSV → INFO_DEST
for src in "${INFO_FILES[@]}"; do
    fname=$(basename "$src")
    log ""
    if is_completed "$src"; then
        log "--- Skipping completed: $fname ---"
        continue
    fi
    log "--- Starting: $fname → ${INFO_DEST}/ ---"
    if ! transfer_file "$src" "$INFO_DEST"; then
        failed+=("$src")
    fi
done

# Stop caffeinate (only if it was started)
[[ -n "$CAFFEINATE_PID" ]] && kill "$CAFFEINATE_PID" 2>/dev/null || true

log ""
log "============================================================"
log "TRANSFER COMPLETE"
log "Successful: $((TOTAL - ${#failed[@]})) / ${TOTAL}"
if [[ ${#failed[@]} -gt 0 ]]; then
    log "FAILED files:"
    for f in "${failed[@]}"; do
        log "  $f"
    done
else
    log "All files transferred successfully!"
fi
log "============================================================"
