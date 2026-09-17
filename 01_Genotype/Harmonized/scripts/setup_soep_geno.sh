#!/usr/bin/env bash
# =============================================================================
# setup_soep_geno.sh — Decrypt and merge the per-chromosome SOEP-G ext dataset
# =============================================================================
# Run this ONCE on the cluster BEFORE the main pipeline (after transferring files).
# Decrypts the 66 per-chromosome .gpg files and merges them into a single
# PLINK BED file: SOEP-G.b37.ext_chrall.{bed,bim,fam}
#
# USAGE:
#   bash scripts/setup_soep_geno.sh --passphrase-file ~/soep_passphrase.txt
#
# Create the passphrase file first (chmod 600 to protect it):
#   echo 'YOUR_PASSPHRASE' > ~/soep_passphrase.txt
#   chmod 600 ~/soep_passphrase.txt
#   (delete it after setup is complete)
# =============================================================================
set -euo pipefail

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SUBMIT_DIR}/config.sh"
ensure_dirs

PASSPHRASE_FILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --passphrase-file) PASSPHRASE_FILE="$2"; shift 2 ;;
        *) log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$PASSPHRASE_FILE" ]] || [[ ! -f "$PASSPHRASE_FILE" ]]; then
    log_error "ERROR: Provide --passphrase-file <path>"
    log_error "Create it with:"
    log_error "  echo 'PASSPHRASE' > ~/soep_passphrase.txt"
    log_error "  chmod 600 ~/soep_passphrase.txt"
    exit 1
fi

PASSPHRASE=$(head -1 "$PASSPHRASE_FILE")
GENO_DIR="$(dirname "$GENO_SOEP")"
EXT_PREFIX="SOEP-G.b37.mildQC.hrc1-1_imp"
DECRYPT_DIR="${GENO_DIR}/soep_ext_decrypted"
MERGE_LIST="${DECRYPT_DIR}/merge_list.txt"

log_info "=== setup_soep_geno.sh: Decrypt and merge SOEP-G ext per-chromosome files ==="
log_info "Input:   ${GENO_DIR}/${EXT_PREFIX}.chr*.dose.vcf.{bed,bim,fam}.gpg"
log_info "Output:  ${GENO_SOEP}.{bed,bim,fam}"

# --- Check that the encrypted source files are present -----------------------
log_info "Checking encrypted source files..."
missing=0
for CHR in $(seq 1 22); do
    for ext in bed bim fam; do
        f="${GENO_DIR}/${EXT_PREFIX}.chr${CHR}.dose.vcf.${ext}.gpg"
        if [[ ! -f "$f" ]]; then
            log_error "  Missing: $f"
            missing=$((missing + 1))
        fi
    done
done
if [[ $missing -gt 0 ]]; then
    log_error "$missing encrypted files not found. Run transfer_geno_to_tardis_workbench.sh first."
    exit 1
fi
log_info "All 66 encrypted source files present."

# Skip if already done
if [[ -f "${GENO_SOEP}.bed" ]] && [[ -f "${GENO_SOEP}.bim" ]] && [[ -f "${GENO_SOEP}.fam" ]]; then
    n_samples=$(wc -l < "${GENO_SOEP}.fam")
    n_variants=$(wc -l < "${GENO_SOEP}.bim")
    log_info "SKIP: Merged file already exists — $n_samples samples, $n_variants variants"
    log_info "Delete ${GENO_SOEP}.{bed,bim,fam} to force re-run."
    exit 0
fi

# --- Step 1: Decrypt per-chromosome files ------------------------------------
mkdir -p "$DECRYPT_DIR"
> "$MERGE_LIST"

log_info "Decrypting 22 chromosomes (66 files)..."
for CHR in $(seq 1 22); do
    prefix="${EXT_PREFIX}.chr${CHR}.dose.vcf"
    for ext in bed bim fam; do
        src="${GENO_DIR}/${prefix}.${ext}.gpg"
        dst="${DECRYPT_DIR}/${prefix}.${ext}"
        if [[ -f "$dst" ]]; then
            log_info "  SKIP decrypt: chr${CHR} .${ext} (exists)"
        else
            gpg --batch --passphrase "$PASSPHRASE" --decrypt "$src" > "$dst" 2>/dev/null
            log_info "  Decrypted: chr${CHR} .${ext}"
        fi
    done
    echo "${DECRYPT_DIR}/${prefix}" >> "$MERGE_LIST"
done
log_info "All 22 chromosomes decrypted."

# --- Step 2: Merge with PLINK2 -----------------------------------------------
log_info "Merging 22 chromosomes into ${GENO_SOEP}..."

# Use chr1 as the base bfile, remaining chrs as the merge list
HEAD=$(head -1 "$MERGE_LIST")
TAIL_LIST="${DECRYPT_DIR}/merge_tail.txt"
tail -n +2 "$MERGE_LIST" > "$TAIL_LIST"

"$PLINK2" \
    --bfile "$HEAD" \
    --pmerge-list "$TAIL_LIST" bfile \
    --make-bed \
    --set-all-var-ids '@:#:$r:$a' \
    --threads "$PLINK2_THREADS" \
    --out "$GENO_SOEP" \
    2>&1 | tee "${LOG_DIR}/setup_soep_merge.log"

if [[ ! -f "${GENO_SOEP}.bed" ]]; then
    log_error "Merge failed — output not found. Check ${LOG_DIR}/setup_soep_merge.log"
    exit 1
fi

n_samples=$(wc -l < "${GENO_SOEP}.fam")
n_variants=$(wc -l < "${GENO_SOEP}.bim")
log_info "Merged: $n_samples samples, $n_variants variants"

# --- Step 3: Clean up decrypted per-chr files --------------------------------
log_info "Cleaning up decrypted per-chromosome files..."
rm -rf "$DECRYPT_DIR"

log_info "=== setup_soep_geno.sh complete ==="
log_info "Merged genotype: ${GENO_SOEP}.{bed,bim,fam}"
log_info "You can now delete the passphrase file: rm -f $PASSPHRASE_FILE"
