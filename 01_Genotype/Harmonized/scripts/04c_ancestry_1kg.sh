#!/usr/bin/env bash
# =============================================================================
# 04c_ancestry_1kg.sh — In-house European-ancestry keep-lists via 1000G projection
# =============================================================================
# Reference-based ancestry classification: project each cohort onto principal
# components trained on the 1000 Genomes phase-3 reference, then keep
# individuals classified as European (see 04c_classify_eur.R).
#
# Runs for the cohorts whose harmonized genotypes are available:
#   BASEII, SHIP0, SHIPTD, TWINLIFE.
# SOEP is not projected here — its genotypes are not available to project, so
# SOEP uses a keep-list computed elsewhere
# (see data/1Kg_exclusion_IDs/EUR_keep_IDs_SOEP.tsv and its README).
#
# STAGES
#   A. Build the 1000G reference PLINK set (once): VCF -> bed, chr:pos IDs,
#      biallelic ACGT SNPs, drop strand-ambiguous SNPs, MAF filter, exclude
#      long-range LD, LD-prune.
#   B. Reference PCA with allele weights (PLINK2 --pca allele-wts).
#   C. Project the reference samples AND each cohort onto those PCs via --score
#      (variance-standardize), so all sit on one common projected scale.
#   D. Hand off to 04c_classify_eur.R for EUR classification + keep-lists.
#
# USAGE
#   REF_1KG_VCF_DIR=/path/to/1000G_phase3_b37 bash scripts/04c_ancestry_1kg.sh
#
# INPUTS
#   $REF_1KG_VCF_DIR/ALL.chr{1..22}.*.genotypes.vcf.gz  (from the EBI release)
#   $REF_1KG_VCF_DIR/integrated_call_samples_v3.20130502.ALL.panel
#   ${GENO_OUT}/<COHORT>_harmonized_shared.{bed,bim,fam}
#
# OUTPUTS (under ${OUT_ROOT}/ancestry/)
#   ref_1kg_pruned.{bed,bim,fam}          reference SNP set used for PCA
#   ref_1kg.eigenvec.allele / .acount     PC loadings + freqs for projection
#   ref_projected_pcs.tsv                 reference sample PCs + super_pop label
#   <KEY>_projected_pcs.tsv               per-cohort projected PCs (KEY = merge key)
#
# DEPENDS ON: step 03 (harmonized bfiles). Reference download is manual (public).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04c: 1000G-EUR ancestry projection ==="

# --- Parameters --------------------------------------------------------------
REF_1KG_VCF_DIR="${REF_1KG_VCF_DIR:-${PROJ_ROOT}/reference/1000G_phase3_b37}"
PANEL="${REF_1KG_VCF_DIR}/integrated_call_samples_v3.20130502.ALL.panel"
ANCESTRY_OUT="${OUT_ROOT}/ancestry"
REF_MAF=0.01
N_PCS=20
mkdir -p "$ANCESTRY_OUT"

# Cohorts to project (SOEP excluded — see header). Map COHORTS-array name ->
# merge key / keep-list filename (merge.R uses SHIPTd / TwinLife, not the
# all-caps array spelling).
PROJECT_COHORTS=("BASEII" "SHIP0" "SHIPTD" "TWINLIFE")
declare -A MERGE_KEY=( [BASEII]=BASEII [SHIP0]=SHIP0 [SHIPTD]=SHIPTd [TWINLIFE]=TwinLife )

# Long-range LD regions (hg19; Price et al. 2008 / Anderson et al. 2010).
# Reuse step 04's file if present, else create it here so 04c is self-contained.
EXCLUDE_REGIONS="${PCS_OUT}/long_range_ld_hg19.txt"
if [[ ! -f "$EXCLUDE_REGIONS" ]]; then
    mkdir -p "$PCS_OUT"
    cat > "$EXCLUDE_REGIONS" << 'EOF'
1 48000000 52000000
2 86000000 100500000
2 134500000 138000000
2 183000000 190000000
3 47500000 50000000
3 83500000 87000000
3 89000000 97500000
5 44500000 50500000
5 98000000 100500000
5 129000000 132000000
5 135500000 138500000
6 25500000 33500000
6 57000000 64000000
6 140000000 142500000
7 55000000 66000000
8 8000000 12000000
8 43000000 50000000
8 112000000 115000000
10 37000000 43000000
11 46000000 57000000
11 87500000 90500000
12 33000000 40000000
12 109500000 112000000
20 32000000 34500000
EOF
    log_info "Created long-range LD exclusion regions file"
fi

# =============================================================================
# STAGE A — build the 1000G reference PLINK set (skipped if already present)
# =============================================================================
REF_PRUNED="${ANCESTRY_OUT}/ref_1kg_pruned"
if [[ -f "${REF_PRUNED}.bed" ]]; then
    log_info "Reference set already built: ${REF_PRUNED}.bed — skipping stage A"
else
    [[ -f "$PANEL" ]] || { log_error "1000G panel not found: $PANEL"; exit 1; }

    log_info "Converting 1000G per-chromosome VCFs -> PLINK (chr:pos IDs, biallelic ACGT SNPs)..."
    : > "${ANCESTRY_OUT}/merge_list.txt"
    for chr in $(seq 1 22); do
        if [[ -f "${ANCESTRY_OUT}/ref_chr${chr}.bed" ]]; then
            log_info "  chr${chr} already converted — skipping"
        else
            vcf=$(ls "${REF_1KG_VCF_DIR}"/ALL.chr${chr}.*.genotypes.vcf.gz 2>/dev/null | head -1 || true)
            [[ -n "$vcf" && -f "$vcf" ]] || { log_error "Missing VCF for chr${chr} in $REF_1KG_VCF_DIR"; exit 1; }
            "$PLINK2" \
                --vcf "$vcf" \
                --set-all-var-ids "$PLINK2_SET_VAR_IDS" \
                --new-id-max-allele-len 100 missing \
                --max-alleles 2 --min-alleles 2 \
                --snps-only just-acgt \
                --autosome \
                --maf "$REF_MAF" \
                --rm-dup "$PLINK2_RM_DUP" \
                --make-bed \
                --threads "$PLINK2_THREADS" \
                --out "${ANCESTRY_OUT}/ref_chr${chr}" \
                2>&1 | tee "${LOG_DIR}/anc_ref_chr${chr}.log"
        fi
        echo "${ANCESTRY_OUT}/ref_chr${chr}" >> "${ANCESTRY_OUT}/merge_list.txt"
    done

    if [[ -f "${ANCESTRY_OUT}/ref_1kg_all.bed" ]]; then
        log_info "Merged reference already present — skipping merge"
    else
        log_info "Merging chromosomes..."
        "$PLINK19" \
            --merge-list "${ANCESTRY_OUT}/merge_list.txt" \
            --make-bed \
            --out "${ANCESTRY_OUT}/ref_1kg_all" \
            2>&1 | tee "${LOG_DIR}/anc_ref_merge.log"
    fi

    # Drop strand-ambiguous SNPs (A/T, C/G) — they cannot be reliably aligned
    # across datasets by allele code alone, so they corrupt projection.
    if [[ ! -f "${ANCESTRY_OUT}/ambiguous.snps" ]]; then
        awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C"){print $2}' \
            "${ANCESTRY_OUT}/ref_1kg_all.bim" > "${ANCESTRY_OUT}/ambiguous.snps"
    fi
    log_info "  Strand-ambiguous SNPs to drop: $(wc -l < "${ANCESTRY_OUT}/ambiguous.snps")"

    # Exclude long-range LD regions. Use PLINK2 here: PLINK 1.9's --exclude range
    # requires a 4th (range-ID) column, whereas PLINK2 accepts the 3-column
    # chr/start/end file (matching step 04's usage).
    log_info "Excluding long-range LD + ambiguous, then LD-pruning..."
    "$PLINK2" \
        --bfile "${ANCESTRY_OUT}/ref_1kg_all" \
        --exclude range "$EXCLUDE_REGIONS" \
        --make-bed --out "${ANCESTRY_OUT}/ref_1kg_lrld" \
        --threads "$PLINK2_THREADS" \
        2>&1 | tee "${LOG_DIR}/anc_ref_lrld.log"
    "$PLINK19" \
        --bfile "${ANCESTRY_OUT}/ref_1kg_lrld" \
        --exclude "${ANCESTRY_OUT}/ambiguous.snps" \
        --indep-pairwise 1000 50 0.05 \
        --out "${ANCESTRY_OUT}/ref_prune" \
        2>&1 | tee "${LOG_DIR}/anc_ref_prune.log"
    "$PLINK19" \
        --bfile "${ANCESTRY_OUT}/ref_1kg_lrld" \
        --extract "${ANCESTRY_OUT}/ref_prune.prune.in" \
        --make-bed --out "$REF_PRUNED" \
        2>&1 | tee "${LOG_DIR}/anc_ref_pruned.log"
    log_info "  Reference SNPs for PCA: $(wc -l < "${REF_PRUNED}.bim")"

    rm -f "${ANCESTRY_OUT}"/ref_chr*.{bed,bim,fam,log} \
          "${ANCESTRY_OUT}"/ref_1kg_all.{bed,bim,fam,log} \
          "${ANCESTRY_OUT}"/ref_1kg_lrld.{bed,bim,fam,log} \
          "${ANCESTRY_OUT}"/ref_prune.{prune.in,prune.out,log}
fi

# Clean variant-ID list for --extract (PLINK expects IDs, not a .bim).
[[ -f "${REF_PRUNED}.snps" ]] || awk '{print $2}' "${REF_PRUNED}.bim" > "${REF_PRUNED}.snps"

# =============================================================================
# STAGE B — reference PCA with allele weights
# =============================================================================
REF_PCA="${ANCESTRY_OUT}/ref_1kg"
if [[ ! -f "${REF_PCA}.eigenvec.allele" ]]; then
    log_info "Computing reference PCA (allele weights)..."
    "$PLINK2" \
        --bfile "$REF_PRUNED" \
        --freq counts \
        --pca "$N_PCS" allele-wts \
        --threads "$PLINK2_THREADS" \
        --out "$REF_PCA" \
        2>&1 | tee "${LOG_DIR}/anc_ref_pca.log"
fi

# Resolve ID / A1 / PC1 column positions by NAME (robust across PLINK2 builds;
# same logic as step 04).
eav_id=""; eav_a1=""; eav_pc1=""
IFS=$'\t' read -r -a eav_cols < "${REF_PCA}.eigenvec.allele"
for eav_i in "${!eav_cols[@]}"; do
    case "${eav_cols[$eav_i]}" in
        \#ID|ID) eav_id=$((eav_i + 1)) ;;
        A1)      eav_a1=$((eav_i + 1)) ;;
        PC1)     eav_pc1=$((eav_i + 1)) ;;
    esac
done
[[ -n "$eav_id" && -n "$eav_a1" && -n "$eav_pc1" ]] || { log_error "Cannot locate ID/A1/PC1 in eigenvec.allele"; exit 1; }

# --- Reusable projection: --score onto the reference axes ---------------------
project_onto_ref() {  # $1 = source bfile prefix, $2 = output PC-table path
    # Projects onto the reference PCs using the SHARED common SNP set
    # ($COMMON_SNPS), so the reference and every cohort sum over the SAME
    # variants and land on a comparable PC scale.
    local src="$1" out_tsv="$2" tmp="${2%.tsv}"
    "$PLINK2" \
        --bfile "$src" \
        --extract "$COMMON_SNPS" \
        --read-freq "${REF_PCA}.acount" \
        --score "${REF_PCA}.eigenvec.allele" "$eav_id" "$eav_a1" header-read \
                no-mean-imputation variance-standardize \
                cols=fid,scoresums \
        --score-col-nums "${eav_pc1}-$((eav_pc1 + N_PCS - 1))" \
        --threads "$PLINK2_THREADS" \
        --out "$tmp" \
        2>&1 | tee "${LOG_DIR}/anc_proj_$(basename "$tmp").log"
    # .sscore -> tidy PC table (SCORE<n>_SUM or PC<n>_SUM -> PC<n>)
    awk 'BEGIN{FS=OFS="\t"}
         NR==1{for(i=1;i<=NF;i++){if($i~/_SUM$/){sub(/_SUM$/,"",$i); sub(/^SCORE/,"PC",$i)}}}
         {print}' "${tmp}.sscore" > "$out_tsv"
    rm -f "${tmp}.sscore" "${tmp}.log"
}

# Build ONE common allele-matched SNP set (reference ∩ every cohort). Each
# member must share the same chr:pos ID AND allele pair (swap allowed) with the
# reference in ALL cohorts — this both removes chr:pos collisions and guarantees
# every projection sums the identical variant set (comparable PC scale).
COMMON_SNPS="${ANCESTRY_OUT}/common_matched.snps"
cp "${REF_PRUNED}.snps" "$COMMON_SNPS"
for cohort in "${PROJECT_COHORTS[@]}"; do
    bfile=$(get_harmonized_bfile "$cohort")
    [[ -f "${bfile}.bed" ]] || continue
    awk 'NR==FNR{r1[$2]=$5; r2[$2]=$6; next}
         ($2 in r1){ if(($5==r1[$2]&&$6==r2[$2])||($5==r2[$2]&&$6==r1[$2])) print $2 }' \
        "${REF_PRUNED}.bim" "${bfile}.bim" > "${ANCESTRY_OUT}/_match_${cohort}.tmp"
    awk 'NR==FNR{k[$1]=1; next} ($1 in k)' \
        "${ANCESTRY_OUT}/_match_${cohort}.tmp" "$COMMON_SNPS" > "${COMMON_SNPS}.tmp" \
        && mv "${COMMON_SNPS}.tmp" "$COMMON_SNPS"
    rm -f "${ANCESTRY_OUT}/_match_${cohort}.tmp"
done
log_info "Common allele-matched SNPs (reference and all cohorts): $(wc -l < "$COMMON_SNPS")"

# =============================================================================
# STAGE C — project reference samples + each cohort
# =============================================================================
log_info "Projecting reference samples..."
project_onto_ref "$REF_PRUNED" "${ANCESTRY_OUT}/ref_projected_pcs_raw.tsv"
# attach super_pop label from the panel (sample IID -> super_pop)
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{ if(FNR>1) sp[$1]=$3; next }                    # panel: sample pop super_pop gender
     FNR==1{ print $0, "super_pop"; next }
     { iid=$2; print $0, (iid in sp ? sp[iid] : "NA") }' \
    "$PANEL" "${ANCESTRY_OUT}/ref_projected_pcs_raw.tsv" \
    > "${ANCESTRY_OUT}/ref_projected_pcs.tsv"
rm -f "${ANCESTRY_OUT}/ref_projected_pcs_raw.tsv"

for cohort in "${PROJECT_COHORTS[@]}"; do
    key="${MERGE_KEY[$cohort]}"
    bfile=$(get_harmonized_bfile "$cohort")
    if [[ ! -f "${bfile}.bed" ]]; then
        log_warn "Harmonized genotype not found for $cohort — skipping"; continue
    fi
    log_info "Projecting $cohort (-> key $key)..."
    project_onto_ref "$bfile" "${ANCESTRY_OUT}/${key}_projected_pcs.tsv"
    log_info "  $(($(wc -l < "${ANCESTRY_OUT}/${key}_projected_pcs.tsv") - 1)) samples projected"
done

log_info "=== Step 04c projection complete. Next: Rscript scripts/04c_classify_eur.R ==="
log_info "    Reference PCs + labels : ${ANCESTRY_OUT}/ref_projected_pcs.tsv"
log_info "    Cohort PCs             : ${ANCESTRY_OUT}/<KEY>_projected_pcs.tsv"
