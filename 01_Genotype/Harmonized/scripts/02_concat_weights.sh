#!/usr/bin/env bash
# =============================================================================
# 02_concat_weights.sh — Concatenate per-chromosome SBayesR outputs and
#                         convert to PLINK2 scoring format (chr:pos IDs)
# =============================================================================
# USAGE:
#   bash scripts/02_concat_weights.sh
#   (typically via `bash run_pipeline.sh --step 02`)
#
# INPUTS:
#   ${SBAYESR_OUT}/<TAG>/<TAG>_chr{1..22}_sbayesR.snpRes
#                                            — per-chr SBayesR posteriors
#
# OUTPUTS:
#   ${SBAYESR_OUT}/<TAG>/<TAG>_allchr_sbayesR.snpRes — concatenated .snpRes
#   ${WEIGHTS_OUT}/<TAG>_sbayesR_chrpos.txt          — 3-col scoring file
#                                                       (chr:pos / A1 / BETA)
#
# DEPENDS ON: step 01 (SBayesR per-chr .snpRes files must exist).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

# Helper functions in this script are used in command substitution. Route their
# progress logging to stderr so stdout can be reserved for returned values.
log_info() {
    printf '[%s] [INFO] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

log_error() {
    printf '[%s] [ERROR] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

log_info "=== Step 02: Concatenating SBayesR outputs and preparing scoring weights ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"

# --- Helper: concatenate per-chr .snpRes into one file -----------------------
concat_snpres() {
    local tag="$1"
    local out_dir="${SBAYESR_OUT}/${tag}"
    local concat_file="${out_dir}/${tag}_allchr_sbayesR.snpRes"

    if [[ -f "$concat_file" ]]; then
        log_info "SKIP: Concatenated file exists: $concat_file"
        echo "$concat_file"
        return 0
    fi

    # Check all 22 chromosomes are present
    local missing=0
    for chr in $(seq 1 22); do
        if [[ ! -f "${out_dir}/${tag}_chr${chr}_sbayesR.snpRes" ]]; then
            log_error "Missing: ${out_dir}/${tag}_chr${chr}_sbayesR.snpRes"
            missing=$((missing + 1))
        fi
    done
    if [[ $missing -gt 0 ]]; then
        log_error "$tag: $missing chromosomes missing. Cannot concatenate."
        return 1
    fi

    # Concatenate: header from chr1, then data from all chromosomes
    head -1 "${out_dir}/${tag}_chr1_sbayesR.snpRes" > "$concat_file"
    for chr in $(seq 1 22); do
        tail -n +2 "${out_dir}/${tag}_chr${chr}_sbayesR.snpRes" >> "$concat_file"
    done

    local n_snps
    n_snps=$(( $(wc -l < "$concat_file") - 1 ))
    log_info "Concatenated $tag: $n_snps SNPs across 22 chromosomes"
    echo "$concat_file"
}

# --- Helper: convert .snpRes to chr:pos scoring format -----------------------
# The .snpRes Name column can be rsID or chr:pos:ref:alt.
# We create a scoring file with chr:pos as variant ID (col 2 of .snpRes
# is replaced with chr{Chrom}:{Position}).
convert_to_chrpos() {
    local input_snpres="$1"
    local output_weight="$2"

    if [[ -f "$output_weight" ]]; then
        log_info "SKIP: Weight file exists: $output_weight"
        return 0
    fi

    # Read .snpRes, create scoring file with columns: chrpos A1 A1Effect
    # .snpRes columns: Index Name Chrom Position A1 A2 A1Frq A1Effect SE VarExplained PIP LastSampleEff
    # We want: chr{Chrom}:{Position}  A1  A1Effect
    awk 'BEGIN {OFS="\t"}
         NR==1 {print "SNP", "A1", "BETA"; next}
         {print "chr" $3 ":" $4, $5, $8}' "$input_snpres" > "$output_weight"

    local n_snps
    n_snps=$(( $(wc -l < "$output_weight") - 1 ))
    log_info "Converted to chrpos scoring format: $output_weight ($n_snps SNPs)"
}

# --- Process standard traits -------------------------------------------------
N_PROCESSED=0
for trait in "${TRAIT_NAMES[@]}"; do
    tag="${TRAIT_SBAYESR_TAG[$trait]}"
    log_info "Processing trait: $trait (tag: $tag)"

    concat_file=$(concat_snpres "$tag") || continue

    # Create chr:pos scoring weight
    chrpos_weight="${WEIGHTS_OUT}/${tag}_sbayesR_chrpos.txt"
    convert_to_chrpos "$concat_file" "$chrpos_weight"

    N_PROCESSED=$((N_PROCESSED + 1))
done

# --- Log decisions -----------------------------------------------------------
{
    echo ""
    echo "## Step 02: Weight Preparation ($(date '+%Y-%m-%d %H:%M'))"
    echo ""
    echo "- Concatenated per-chromosome SBayesR .snpRes files into genome-wide files"
    echo "- Converted variant IDs to chr:pos format for PLINK2 scoring"
    echo "- Standard traits processed: ${N_PROCESSED}"
    echo "- Scoring file format: SNP(chr:pos) A1 BETA"
    echo ""
} >> "$DECISION_LOG"

log_info "=== Step 02 complete: $N_PROCESSED trait weight files prepared ==="
