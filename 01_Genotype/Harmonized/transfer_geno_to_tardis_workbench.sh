#!/usr/bin/env bash

set -u

# =============================================================================
# transfer_geno_to_tardis_workbench.sh — Transfer genotype files + TwinLife R²
# from Posit Workbench / ARC to Tardis.
# =============================================================================
# - Linux/Workbench friendly: no /Volumes, no open, no caffeinate
# - Skips files already present on Tardis when byte sizes match
# - Resumes interrupted transfers via rsync --partial
# - Keeps a local completed-state file so reruns continue cleanly
#
# Transfers two streams (mirrors transfer_geno_to_tardis.sh):
#   GENO_DEST  — per-cohort imputed bfiles (.bed/.bim/.fam)
#   INFO_DEST  — TwinLife imputation-R² TSV deposited on the TwinLife share
#
# USAGE:
#   1. Edit SOURCE_ROOTS below to match the mount points on Workbench.
#   2. Run a dry check:
#        bash transfer_geno_to_tardis_workbench.sh --probe
#   3. Run the transfer:
#        nohup bash transfer_geno_to_tardis_workbench.sh \
#          > transfer_geno_workbench.nohup.out 2>&1 &
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load TARDIS_SSH, TARDIS_PROJECT_ROOT from user_config.sh (or run the
# interactive bootstrap on first invocation).
# shellcheck source=_bootstrap_user_config.sh
source "${SCRIPT_DIR}/_bootstrap_user_config.sh"

TARDIS="${TARDIS_SSH}"
GENO_DEST="${TARDIS_PROJECT_ROOT}/data/geno"
INFO_DEST="${TARDIS_PROJECT_ROOT}/data/info_files/TWINLIFE"
LOG="${SCRIPT_DIR}/transfer_geno_workbench.log"
STATE="${SCRIPT_DIR}/transfer_geno_workbench.completed"
MAX_RETRIES=10
RETRY_WAIT=30
PROBE_ONLY=0

if [[ "${1:-}" == "--probe" ]]; then
    PROBE_ONLY=1
fi

# Source roots visible after `kinit` + `biomount` on Posit Workbench / ARC.
# Override SHARE_MOUNT_ROOT in the environment; the default is a placeholder.
SHARE_MOUNT_ROOT="${SHARE_MOUNT_ROOT:-/path/to/share_mount_root}"
declare -A SOURCE_ROOTS=(
    ["<ship-share>"]="${SHARE_MOUNT_ROOT}/<ship-share>"
    ["<soep-share>"]="${SHARE_MOUNT_ROOT}/<soep-share>"
    ["<twinlife-share>"]="${SHARE_MOUNT_ROOT}/<twinlife-share>"
    ["<base2-share>"]="${SHARE_MOUNT_ROOT}/<base2-share>"
)

declare -a FILE_SPECS=(
    "<ship-share>|private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.bed"
    "<ship-share>|private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.bim"
    "<ship-share>|private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose.fam"
    "<ship-share>|private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.bed"
    "<ship-share>|private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.bim"
    "<ship-share>|private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall.fam"
    # SOEP-G extended dataset (encrypted per-chromosome; decrypt+merge on Tardis with setup_soep_geno.sh)
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr1.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr1.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr1.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr2.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr2.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr2.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr3.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr3.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr3.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr4.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr4.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr4.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr5.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr5.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr5.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr6.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr6.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr6.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr7.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr7.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr7.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr8.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr8.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr8.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr9.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr9.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr9.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr10.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr10.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr10.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr11.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr11.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr11.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr12.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr12.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr12.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr13.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr13.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr13.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr14.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr14.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr14.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr15.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr15.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr15.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr16.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr16.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr16.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr17.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr17.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr17.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr18.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr18.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr18.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr19.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr19.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr19.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr20.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr20.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr20.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr21.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr21.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr21.dose.vcf.fam.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr22.dose.vcf.bed.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr22.dose.vcf.bim.gpg"
    "<soep-share>|private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr22.dose.vcf.fam.gpg"
    "<twinlife-share>|private/data/Gendata/001_processed/TwinLife.b37.chrall.bed"
    "<twinlife-share>|private/data/Gendata/001_processed/TwinLife.b37.chrall.bim"
    "<twinlife-share>|private/data/Gendata/001_processed/TwinLife.b37.chrall.fam"
    "<base2-share>|private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.bed"
    "<base2-share>|private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.bim"
    "<base2-share>|private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY.fam"
)

# TwinLife imputation-R² TSV (deposited on the TwinLife share) →
# ${INFO_DEST}. Replaces the missing Minimac .info.gz for TwinLife.
declare -a INFO_FILE_SPECS=(
    "<twinlife-share>|private/data/Gendata/rawgendata/TwinLife_imputation_r2.tsv.gz"
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

file_size() {
    local path="$1"
    stat -c '%s' "$path"
}

remote_file_size() {
    local fname="$1"
    local dest="$2"
    ssh "$TARDIS" "if [ -f '${dest}/${fname}' ]; then stat -c '%s' '${dest}/${fname}'; else echo MISSING; fi"
}

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

    local spec dataset rel_path src
    for spec in "${FILE_SPECS[@]}" "${INFO_FILE_SPECS[@]}"; do
        dataset="${spec%%|*}"
        rel_path="${spec#*|}"
        src="${SOURCE_ROOTS[$dataset]}/${rel_path}"
        if grep -Fq "[$(basename "$src")] DONE" "$LOG" 2>/dev/null; then
            mark_completed "$src"
        fi
    done
}

resolve_source_path() {
    local dataset="$1"
    local rel_path="$2"
    local root="${SOURCE_ROOTS[$dataset]:-}"

    if [[ -z "$root" ]]; then
        return 1
    fi

    local src="${root}/${rel_path}"
    if [[ -f "$src" ]]; then
        echo "$src"
        return 0
    fi

    return 1
}

probe_sources() {
    local missing=0
    local spec dataset rel_path

    log "Probing configured source paths"
    for spec in "${FILE_SPECS[@]}" "${INFO_FILE_SPECS[@]}"; do
        dataset="${spec%%|*}"
        rel_path="${spec#*|}"
        if src=$(resolve_source_path "$dataset" "$rel_path"); then
            log "[OK] $src"
        else
            log "[MISSING] dataset=${dataset} rel_path=${rel_path} root=${SOURCE_ROOTS[$dataset]}"
            missing=1
        fi
    done

    return "$missing"
}

transfer_file() {
    local src="$1"
    local dest="$2"
    local fname
    fname=$(basename "$src")
    local local_size
    local_size=$(file_size "$src")
    local attempt=0

    while [[ $attempt -lt $MAX_RETRIES ]]; do
        ((attempt++))
        log "[$fname] Attempt $attempt/$MAX_RETRIES"

        if [[ ! -f "$src" ]]; then
            log "[$fname] ERROR: Source file not found: $src"
            sleep "$RETRY_WAIT"
            continue
        fi

        local remote_size
        remote_size=$(remote_file_size "$fname" "$dest")

        if [[ "$remote_size" == "$local_size" ]]; then
            log "[$fname] Already present on Tardis with matching size ($local_size bytes); skipping."
            mark_completed "$src"
            return 0
        fi

        rsync -avP --partial \
            "$src" \
            "${TARDIS}:${dest}/" \
            2>&1 | tee -a "$LOG"

        if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
            mark_completed "$src"
            log "[$fname] DONE"
            return 0
        fi

        log "[$fname] rsync failed (attempt $attempt). Retrying in ${RETRY_WAIT}s..."
        sleep "$RETRY_WAIT"
    done

    log "[$fname] FAILED after $MAX_RETRIES attempts."
    return 1
}

TOTAL=$(( ${#FILE_SPECS[@]} + ${#INFO_FILE_SPECS[@]} ))

log "============================================================"
log "Starting genotype transfer to Tardis from Workbench"
log "Genotype bfiles   → ${GENO_DEST}   (${#FILE_SPECS[@]} files)"
log "TwinLife R² TSV   → ${INFO_DEST}   (${#INFO_FILE_SPECS[@]} file)"
log "Log:   $LOG"
log "State: $STATE"
log "============================================================"

# Pre-create both Tardis destination dirs (rsync does NOT auto-create parents).
ssh "$TARDIS" "mkdir -p '${GENO_DEST}' '${INFO_DEST}'" \
    || { log "ERROR: could not ssh to Tardis ($TARDIS) to create dest dirs"; exit 1; }

hydrate_completed_state

if ! probe_sources; then
    if [[ $PROBE_ONLY -eq 1 ]]; then
        exit 1
    fi
    log "ERROR: Some source paths are not reachable from this machine."
    log "Edit SOURCE_ROOTS in this script, then rerun with --probe."
    exit 1
fi

if [[ $PROBE_ONLY -eq 1 ]]; then
    log "Probe completed successfully."
    exit 0
fi

failed=()

# Pass 1: genotype bfiles → GENO_DEST
for spec in "${FILE_SPECS[@]}"; do
    dataset="${spec%%|*}"
    rel_path="${spec#*|}"
    src="${SOURCE_ROOTS[$dataset]}/${rel_path}"
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
for spec in "${INFO_FILE_SPECS[@]}"; do
    dataset="${spec%%|*}"
    rel_path="${spec#*|}"
    src="${SOURCE_ROOTS[$dataset]}/${rel_path}"
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
