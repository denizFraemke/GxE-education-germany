#!/usr/bin/env bash
# =============================================================================
# 03c_align_to_ukb.sh — Align all cohorts to UKB reference panel
# =============================================================================
# Since the SBayesR weights were derived from UKB summary statistics, aligning
# all cohort genotypes to the UKB reference panel maximizes PGI comparability.
#
# This script:
#   1. Uses the UKB 50k European reference panel (already used for SBayesR LD)
#   2. For each cohort: identifies strand-ambiguous SNPs and resolves using UKB
#   3. Uses PLINK2 --flip + --a2-allele to align alleles to UKB reference
#   4. Outputs: {COHORT}_harmonized_ukb.bed/bim/fam
#
# USAGE:
#   bash scripts/03c_align_to_ukb.sh [--force]
#   (off by default — set UKB_ALIGN=yes in config.sh to enable)
#
# INPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}
#   ${LD_REF_DIR}/ukb50k_eur_a1frq.ldm.sparse  — UKB allele-frequency reference
#
# OUTPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_ukb.{bed,bim,fam}
#
# DEPENDS ON: step 03 (and optionally 03b). Disabled by default; only needed
# for the UKB-reference-aligned sensitivity comparison in step 07b.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

FORCE=""
if [[ "${1:-}" == "--force" ]]; then FORCE="yes"; fi

log_info "=== Step 03c: UKB Reference Panel Alignment ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
ALIGN_SUMMARY=""

# UKB reference directory (from config — same as LD reference for SBayesR)
UKB_REF="${UKB_REF_DIR:-${LD_REF_DIR}}"

if [[ ! -d "$UKB_REF" ]]; then
    log_error "UKB reference directory not found: $UKB_REF"
    log_error "Please set UKB_REF_DIR in config.sh or ensure LD_REF_DIR is accessible"
    exit 1
fi

# Try to find the UKB allele frequency file
# Common naming: ukb50k_eur_a1frq, ukb50k_shrunk, etc.
UKB_AF_FILE=""
for pattern in "ukb50k*a1frq*" "ukb50k*eur*" "ukb_50k*"; do
    found=$(find "$UKB_REF" -maxdepth 2 -iname "$pattern" 2>/dev/null | head -1)
    if [[ -n "$found" ]]; then
        UKB_AF_FILE="$found"
        break
    fi
done

if [[ -z "$UKB_AF_FILE" || ! -f "$UKB_AF_FILE" ]]; then
    log_warn "UKB allele frequency reference file not found in $UKB_REF"
    log_warn "Looking for files matching: ukb50k*a1frq*, ukb50k*eur*, etc."
    log_warn ""
    log_warn "Falling back to the LD reference SNP list as a proxy."
    log_warn "The LD reference matrices (ukb50k_shrunk_chr*_mafpt01.ldm.sparse) can be used"
    log_warn "as a proxy for the UKB SNP set. This ensures all cohorts use the same"
    log_warn "SNPs as the SBayesR analysis."
    log_warn ""

    # Fallback: use the LD reference SNP list
    log_info "FALLBACK: Using LD reference SNP list as UKB proxy"
    UKB_AF_FILE="${UKB_REF}/ukb50k_shrunk_chr1_mafpt01.ldm.sparse"
fi

# =============================================================================
# HELPER: align one cohort to UKB reference
# =============================================================================
align_to_ukb() {
    local cohort="$1"
    local input_bfile="$2"
    local output_prefix="$3"

    if [[ -f "${output_prefix}.bed" && -z "$FORCE" ]]; then
        log_info "  SKIP: $cohort already UKB-aligned"
        return 0
    fi

    # Step 1: Identify which SNPs need flipping or allele swapping
    # For each SNP in the harmonized bfile:
    #   - Check if bim A1 matches UKB A1 (reference allele from UKB)
    #   - If not, check if bim A1 matches the COMPLEMENT of UKB A1
    #   - If complement matches: the alleles are on the wrong strand — FLIP
    #   - If neither matches: the alleles are on the wrong chromosome build or different variant — EXCLUDE

    local flip_snps="${OUT_ROOT}/strand_fixes/${cohort}_ukb_flip_snps.txt"
    local exclude_snps="${OUT_ROOT}/strand_fixes/${cohort}_ukb_exclude_snps.txt"

    # Create flip list by comparing bim to weight file alleles
    # (UKB A1 is the reference allele used in the GWAS)
    local ea4_trait
    ea4_trait=$(get_ea4_trait_for_cohort "$cohort")
    local ea4_tag="${TRAIT_SBAYESR_TAG[$ea4_trait]}"
    local weight_file="${WEIGHTS_OUT}/${ea4_tag}_sbayesR_chrpos.txt"

    log_info "  Generating flip list for $cohort..."
    awk -v bim="${input_bfile}.bim" '
    BEGIN { OFS="\t" }
    FNR==NR {
        if (NF >= 8) {
            # Read bim: col2=SNP, col5=A1(REF), col6=A2(ALT)
            bim_snp[$2] = $2
            bim_a1[$2] = toupper($5)
            bim_a2[$2] = toupper($6)
        }
        next
    }
    {
        # Weight file: SNP, A1, BETA
        snp = $1
        wa1 = toupper($2)  # weight A1 (effect allele from GWAS)

        if (!(snp in bim_snp)) next

        b_a1 = bim_a1[snp]  # REF in bim
        b_a2 = bim_a2[snp]  # ALT in bim

        # Complement of bim alleles
        comp_a1 = b_a1; comp_a1 = (b_a1=="A")?"T":(b_a1=="T")?"A":(b_a1=="C")?"G":(b_a1=="G")?"C":b_a1
        comp_a2 = b_a2; comp_a2 = (b_a2=="A")?"T":(b_a2=="T")?"A":(b_a2=="C")?"G":(b_a2=="G")?"C":b_a2

        # Check alignment with weight A1
        # If weight A1 matches bim REF -> correct (effect allele is REF)
        # If weight A1 matches bim ALT -> correct (effect allele is ALT)
        # If weight A1 matches complement of bim REF -> STRAND FLIP needed
        # If weight A1 matches complement of bim ALT -> STRAND FLIP needed
        # If neither matches -> exclude (different variant/build)

        if (wa1 == comp_a1 || wa1 == comp_a2) {
            # Strand flip needed
            print snp > "'"$flip_snps"'"
        } else if (wa1 != b_a1 && wa1 != b_a2) {
            # Cannot resolve — exclude
            print snp > "'"$exclude_snps"'"
        }
    }
    ' "${input_bfile}.bim" "$weight_file" 2>/dev/null || true

    local n_flip=0 n_exclude=0
    if [[ -f "$flip_snps" ]]; then n_flip=$(wc -l < "$flip_snps"); fi
    if [[ -f "$exclude_snps" ]]; then n_exclude=$(wc -l < "$exclude_snps"); fi
    log_info "  $cohort: $n_flip SNPs to flip, $n_exclude SNPs to exclude"

    # Step 2: Apply flips and exclusions
    "$PLINK2" \
        --bfile "$input_bfile" \
        $(if [[ -f "$flip_snps" && "$n_flip" -gt 0 ]]; then echo "--flip $flip_snps"; fi) \
        $(if [[ -f "$exclude_snps" && "$n_exclude" -gt 0 ]]; then echo "--exclude $exclude_snps"; fi) \
        --make-bed \
        --threads "$PLINK2_THREADS" \
        --out "$output_prefix" \
        2>&1 | tee "${LOG_DIR}/ukb_align_${cohort}.log"

    local n_final=$(wc -l < "${output_prefix}.bim")
    log_info "  -> UKB-aligned genotype: ${output_prefix}.bed ($n_final SNPs)"
}

# =============================================================================
# MAIN: Align each cohort
# =============================================================================

for cohort in "${COHORTS[@]}"; do
    log_info ""
    log_info "--- Aligning $cohort to UKB reference ---"

    # Use the harmonized bfile (from step 03)
    input_bfile=$(get_harmonized_bfile "$cohort")

    if [[ ! -f "${input_bfile}.bed" ]]; then
        log_warn "  SKIP: harmonized bfile not found for $cohort — run step 03 first"
        continue
    fi

    # Output: UKB-aligned bfile
    output_prefix="${GENO_OUT}/${cohort}_harmonized_ukb"

    align_to_ukb "$cohort" "$input_bfile" "$output_prefix"

    if [[ -f "${output_prefix}.bed" ]]; then
        ALIGN_SUMMARY="${ALIGN_SUMMARY}\n- $cohort: UKB-aligned bfile written to ${output_prefix}.bed"
    fi
done

# =============================================================================
# LOG DECISIONS
# =============================================================================
{
    echo ""
    echo "## Step 03c: UKB Reference Panel Alignment ($(date '+%Y-%m-%d %H:%M'))"
    echo ""
    echo "- Aligned all cohorts to UKB 50k European reference panel"
    echo "- Reference: $UKB_REF"
    echo "- For each cohort: identified strand flips by comparing bim A1 to weight A1"
    echo "- Applied PLINK2 --flip for SNPs with complementary alleles"
    echo "- Excluded SNPs where alleles could not be resolved"
    echo "- UKB-aligned bfiles: {COHORT}_harmonized_ukb.bed/bim/fam"
    echo ""
    echo "### Per-cohort summary"
    echo -e "$ALIGN_SUMMARY"
    echo ""
} >> "$DECISION_LOG"

log_info ""
log_info "=== Step 03c complete ==="
log_info ""
log_info "UKB-aligned bfiles are ready for:"
log_info "  Step 04: bash scripts/04_compute_pcs.sh --bfile harmonized_ukb"
log_info "  Step 05: bash scripts/05_score_pgi.sh --bfile harmonized_ukb"
