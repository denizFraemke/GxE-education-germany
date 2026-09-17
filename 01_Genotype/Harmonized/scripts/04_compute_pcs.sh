#!/usr/bin/env bash
# =============================================================================
# 04_compute_pcs.sh — Compute within-cohort ancestry PCs
# =============================================================================
# For each cohort:
#   1. Filter to common, autosomal, well-genotyped SNPs
#   2. Exclude long-range LD regions (MHC, inversions)
#   3. LD-prune to get independent SNPs
#   4. Identify a maximum unrelated subset (KING-robust kinship < 0.0884,
#      i.e. less related than 3rd-degree relatives)
#   5. Compute PCs on the unrelated subset only (avoids family-structure bias;
#      especially relevant for TwinLife)
#   6. Project ALL participants onto those PCs so everyone receives PC values
#
# Output: ${cohort}_pcs.eigenvec  with columns #FID IID PC1 ... PC20
#
# USAGE:
#   bash scripts/04_compute_pcs.sh
#   (typically via `bash run_pipeline.sh --step 04`)
#
# INPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}
#
# OUTPUTS:
#   ${OUT_ROOT}/pcs/<COHORT>_pcs.eigenvec   — #FID IID PC1..PC20
#   ${OUT_ROOT}/pcs/<COHORT>_pcs.eigenval   — eigenvalues
#
# DEPENDS ON: step 03 (harmonized bfiles).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04: Computing within-cohort ancestry PCs ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
SUMMARY=""

# --- Long-range LD regions to exclude (hg19/GRCh37) -------------------------
# From Price et al. (2008) and Anderson et al. (2010)
EXCLUDE_REGIONS="${PCS_OUT}/long_range_ld_hg19.txt"
if [[ ! -f "$EXCLUDE_REGIONS" ]]; then
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

# --- PCA settings -----------------------------------------------------------
N_PCS=20
KING_CUTOFF=0.0884   # ~3rd-degree relatives (Manichaikul et al., 2010)

for cohort in "${COHORTS[@]}"; do
    bfile=$(get_harmonized_bfile "$cohort")
    pc_prefix="${PCS_OUT}/${cohort}_pcs"
    pc_out="${pc_prefix}.eigenvec"

    # Skip if already done
    if [[ -f "$pc_out" ]]; then
        n_pcs_done=$(head -1 "$pc_out" | awk '{print NF - 2}')
        log_info "SKIP: $cohort PCs already computed ($n_pcs_done PCs)"
        SUMMARY="${SUMMARY}\n- $cohort: SKIPPED (already done)"
        continue
    fi

    if [[ ! -f "${bfile}.bed" ]]; then
        log_warn "Harmonized genotype not found for $cohort — run step 03 first"
        continue
    fi

    log_info "Computing PCs for $cohort..."

    # --- Step 1: SNP filter — autosomal, MAF > 0.05, geno > 0.98, HWE > 1e-6 -
    FILTERED="${PCS_OUT}/${cohort}_pc_filtered"
    "$PLINK2" \
        --bfile "$bfile" \
        --autosome \
        --maf 0.05 \
        --geno 0.02 \
        --hwe 1e-6 \
        --exclude range "$EXCLUDE_REGIONS" \
        --make-bed \
        --threads "$PLINK2_THREADS" \
        --out "$FILTERED" \
        2>&1 | tee "${LOG_DIR}/pcs_filter_${cohort}.log"

    n_filtered=$(wc -l < "${FILTERED}.bim")
    log_info "  After filtering: $n_filtered variants"

    # --- Step 2: LD pruning (window 1000 kb, step 50, r² < 0.05) -------------
    PRUNED="${PCS_OUT}/${cohort}_pc_pruned"
    "$PLINK19" \
        --bfile "$FILTERED" \
        --indep-pairwise 1000 50 0.05 \
        --out "$PRUNED" \
        2>&1 | tee "${LOG_DIR}/pcs_prune_${cohort}.log"

    n_pruned=$(wc -l < "${PRUNED}.prune.in")
    log_info "  After LD pruning: $n_pruned independent SNPs"

    # --- Step 3: Identify maximum unrelated subset (KING) --------------------
    # PLINK2 --king-cutoff: greedy algorithm to find a maximum subset such
    # that no pair has KING-robust kinship >= cutoff (0.0884 ≈ 3rd-degree).
    # Requires at least 2 samples; works on small samples too.
    UNREL_PREFIX="${pc_prefix}_unrel_set"
    "$PLINK2" \
        --bfile "$FILTERED" \
        --extract "${PRUNED}.prune.in" \
        --king-cutoff "$KING_CUTOFF" \
        --threads "$PLINK2_THREADS" \
        --out "$UNREL_PREFIX" \
        2>&1 | tee "${LOG_DIR}/pcs_king_${cohort}.log"

    UNREL_KEEP="${UNREL_PREFIX}.king.cutoff.in.id"
    n_total=$(wc -l < "${FILTERED}.fam")
    # The .king.cutoff.in.id file has a header line — subtract 1 for true count
    n_unrel=$(( $(wc -l < "$UNREL_KEEP") - 1 ))
    n_excluded=$(( n_total - n_unrel ))
    log_info "  Unrelated subset (KING < ${KING_CUTOFF}): ${n_unrel} of ${n_total} samples (${n_excluded} relateds)"

    # --- Step 4: Compute PCs on the unrelated subset (with allele weights) ---
    # `allele-wts` writes ${prefix}.eigenvec.allele used for projection.
    # `--freq counts` writes ${prefix}.acount, the allele frequencies needed
    # by --read-freq during projection (variance-standardize).
    TRAIN_PREFIX="${pc_prefix}_train"
    "$PLINK2" \
        --bfile "$FILTERED" \
        --extract "${PRUNED}.prune.in" \
        --keep "$UNREL_KEEP" \
        --freq counts \
        --pca "$N_PCS" allele-wts \
        --threads "$PLINK2_THREADS" \
        --out "$TRAIN_PREFIX" \
        2>&1 | tee "${LOG_DIR}/pcs_train_${cohort}.log"

    # --- Step 5: Project ALL samples onto those PC axes ----------------------
    # `cols=fid,scoresums` keeps FID, IID (always emitted by --score), and the
    # SCORE*_SUM columns. (`iid` is NOT a valid --score cols= token — PLINK2
    # errors "Unrecognized ID 'iid'"; IID is mandatory output, not requestable.)
    # SCORE*_SUM (not SCORE*_AVG) matches the scale of the original eigenvecs.
    # `variance-standardize` rescales dosages using the training-set frequencies.
    PROJ_PREFIX="${pc_prefix}_proj"
    # Resolve the ID / A1 / PC1 column POSITIONS from the .eigenvec.allele header
    # by NAME. PLINK2's column order varies by build: a-7LM wrote
    # '#CHROM ID REF ALT A1 PC1..' (A1 in col 5), whereas a-7.1LM inserts a
    # PROVISIONAL_REF? column before A1 (A1 in col 6). Hardcoding "2 5 ... 6-N"
    # then makes --score read the wrong column as the allele -> "mismatching
    # allele codes / No valid variants". Column NAMES are stable across builds.
    eav_id=""; eav_a1=""; eav_pc1=""
    IFS=$'\t' read -r -a eav_cols < "${TRAIN_PREFIX}.eigenvec.allele"
    for eav_i in "${!eav_cols[@]}"; do
        case "${eav_cols[$eav_i]}" in
            \#ID|ID) eav_id=$((eav_i + 1)) ;;
            A1)      eav_a1=$((eav_i + 1)) ;;
            PC1)     eav_pc1=$((eav_i + 1)) ;;
        esac
    done
    if [[ -z "$eav_id" || -z "$eav_a1" || -z "$eav_pc1" ]]; then
        log_error "Cannot locate ID/A1/PC1 columns in ${TRAIN_PREFIX}.eigenvec.allele"
        exit 1
    fi
    "$PLINK2" \
        --bfile "$FILTERED" \
        --extract "${PRUNED}.prune.in" \
        --read-freq "${TRAIN_PREFIX}.acount" \
        --score "${TRAIN_PREFIX}.eigenvec.allele" "$eav_id" "$eav_a1" header-read \
                no-mean-imputation variance-standardize \
                cols=fid,scoresums \
        --score-col-nums "${eav_pc1}-$((eav_pc1 + N_PCS - 1))" \
        --threads "$PLINK2_THREADS" \
        --out "$PROJ_PREFIX" \
        2>&1 | tee "${LOG_DIR}/pcs_project_${cohort}.log"

    # --- Step 6: Reformat .sscore -> .eigenvec (rename SCORE*_SUM -> PC*) ----
    if [[ ! -f "${PROJ_PREFIX}.sscore" ]]; then
        log_error "  Projection output missing for $cohort"
        SUMMARY="${SUMMARY}\n- $cohort: FAILED (projection)"
        continue
    fi

    awk 'BEGIN{FS=OFS="\t"}
         NR==1 {
             for (i=1; i<=NF; i++) {
                 if ($i ~ /^SCORE[0-9]+_SUM$/) {
                     n=$i; sub("SCORE","",n); sub("_SUM","",n);
                     $i = "PC" n
                 }
             }
         }
         { print }' \
        "${PROJ_PREFIX}.sscore" > "$pc_out"

    if [[ -f "$pc_out" ]]; then
        n_samples=$(( $(wc -l < "$pc_out") - 1 ))
        log_info "  PCs written: $n_samples samples x $N_PCS PCs"
        SUMMARY="${SUMMARY}\n- $cohort: $n_samples projected (${n_unrel} train / ${n_excluded} relateds), $n_pruned SNPs, $N_PCS PCs"
    else
        log_error "  PCA output not found for $cohort"
        SUMMARY="${SUMMARY}\n- $cohort: FAILED"
    fi

    # --- Persist eigenvalues for the validation scree ------------------------
    # The raw .eigenval is deleted in cleanup below; keep a tidy per-cohort copy
    # in the deliverable dir so 03_Merge/validate.R can draw a real
    # variance-explained scree (the projected PC scores do NOT carry eigenvalue
    # magnitude). Columns: cohort, PC, eigenval.
    if [[ -f "${TRAIN_PREFIX}.eigenval" ]]; then
        mkdir -p "$FINAL_OUT"
        awk -v co="$cohort" 'BEGIN{OFS="\t"; print "cohort","PC","eigenval"}
                             {print co, "PC"NR, $1}' \
            "${TRAIN_PREFIX}.eigenval" > "${FINAL_OUT}/${cohort}_PC_eigenval.tsv"
        log_info "  Eigenvalues saved: ${cohort}_PC_eigenval.tsv"
    fi

    # --- Clean up intermediate files -----------------------------------------
    rm -f "${FILTERED}.bed" "${FILTERED}.bim" "${FILTERED}.fam" "${FILTERED}.log"
    rm -f "${PRUNED}.prune.in" "${PRUNED}.prune.out" "${PRUNED}.log" "${PRUNED}.nosex"
    rm -f "${UNREL_PREFIX}.king.cutoff.in.id" "${UNREL_PREFIX}.king.cutoff.out.id" \
          "${UNREL_PREFIX}.log"
    rm -f "${TRAIN_PREFIX}.eigenvec" "${TRAIN_PREFIX}.eigenvec.allele" \
          "${TRAIN_PREFIX}.eigenval" "${TRAIN_PREFIX}.acount" "${TRAIN_PREFIX}.log"
    rm -f "${PROJ_PREFIX}.sscore" "${PROJ_PREFIX}.log"
done

# --- Log decisions -----------------------------------------------------------
{
    echo ""
    echo "## Step 04: Within-Cohort Ancestry PCs ($(date '+%Y-%m-%d %H:%M'))"
    echo ""
    echo "- PCs computed within each cohort using PLINK2 --pca"
    echo "- QC filters: autosomal, MAF > 0.05, genotyping rate > 98%, HWE p > 1e-6"
    echo "- Excluded long-range LD regions (Price et al., 2008)"
    echo "- LD pruning: window 1000kb, step 50, r² < 0.05"
    echo "- Relatedness: PCs trained on KING-unrelated subset (kinship < ${KING_CUTOFF},"
    echo "  ~3rd-degree); related individuals are projected onto those axes"
    echo "  (Manichaikul et al., 2010 for the kinship cutoff)"
    echo "- Number of PCs: $N_PCS"
    echo ""
    echo "### Per-cohort summary"
    echo -e "$SUMMARY"
    echo ""
} >> "$DECISION_LOG"

log_info "=== Step 04 complete ==="
