#!/usr/bin/env bash
# =============================================================================
# 04b_compute_cross_cohort_pcs.sh — Compute PCs on pooled cross-cohort sample
# =============================================================================
# Merges all cohorts' harmonized_shared genotypes into one dataset, then
# computes PCs on the combined sample. These cross-cohort PCs capture
# ancestry differences BETWEEN cohorts (e.g., regional German substructure,
# platform batch effects) that within-cohort PCs miss.
#
# Essential for pooled GxE regressions across multiple cohorts.
#
# To avoid family structure leaking into the top PCs (TwinLife is a twin
# study), PCs are trained on a maximum unrelated subset of the merged
# sample (KING-robust kinship < 0.0884, ≈ 3rd-degree relatives) and then
# projected onto everyone, so all participants receive cross-cohort PC values.
#
# NOTE: After strand-flipping (step 03b), some cohorts may have flipped
# alleles for strand-ambiguous SNPs (A/T, C/G). PLINK 1.9 --merge-list
# sees these as 3+ allele variants. We handle this by:
#   1. Attempting merge
#   2. If .missnp file is produced, excluding those SNPs from ALL cohorts
#   3. Re-attempting merge with cleaned files
#
# Output: output/pcs/cross_cohort_pcs.eigenvec  (#FID IID PC1 ... PC20)
#
# USAGE:
#   bash scripts/04b_compute_cross_cohort_pcs.sh
#   (typically via `bash run_pipeline.sh --step 04b`)
#
# INPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_harmonized_shared.{bed,bim,fam}  — all 5 cohorts
#
# OUTPUTS:
#   ${OUT_ROOT}/pcs/cross_cohort_pcs.eigenvec  — #FID IID PC1..PC20 (pooled)
#   ${OUT_ROOT}/pcs/cross_cohort_pcs.eigenval
#
# DEPENDS ON: step 03 (harmonized bfiles for all cohorts).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== Step 04b: Computing cross-cohort PCs ==="

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
N_PCS=20
KING_CUTOFF=0.0884   # ~3rd-degree relatives (Manichaikul et al., 2010)
CROSS_PC_PREFIX="${PCS_OUT}/cross_cohort_pcs"
CROSS_PC_OUT="${CROSS_PC_PREFIX}.eigenvec"

# Skip if already done
if [[ -f "$CROSS_PC_OUT" ]]; then
  n_samples=$(( $(wc -l < "$CROSS_PC_OUT") - 1 ))
  n_pcs_done=$(head -1 "$CROSS_PC_OUT" | awk '{print NF - 2}')
  log_info "SKIP: Cross-cohort PCs already computed ($n_samples samples, $n_pcs_done PCs)"
  log_info "  Delete $CROSS_PC_OUT to force recomputation"
  exit 0
fi

# --- Step 1: Create merge list -----------------------------------------------
MERGE_LIST="${PCS_OUT}/_cross_cohort_merge_list.txt"
: > "$MERGE_LIST"

FIRST_BFILE=""
N_COHORTS=0
for cohort in "${COHORTS[@]}"; do
  bfile=$(get_harmonized_bfile "$cohort")
  if [[ ! -f "${bfile}.bed" ]]; then
    log_error "  $cohort harmonized_shared bfile not found: ${bfile}.bed"
    exit 1
  fi
  if [[ -z "$FIRST_BFILE" ]]; then
    FIRST_BFILE="$bfile"
  else
    echo "${bfile}.bed ${bfile}.bim ${bfile}.fam" >> "$MERGE_LIST"
  fi
  N_COHORTS=$((N_COHORTS + 1))
done

log_info "  Merging $N_COHORTS cohorts (reference: $(basename "$FIRST_BFILE"))..."

# --- Step 2: Merge all cohorts ------------------------------------------------
# Use cleaned per-cohort files if available (from a previous retry), else originals
MERGED="${PCS_OUT}/_cross_cohort_merged"
MISSNP="${MERGED}-merge.missnp"

# Clean up any leftover files from a previous failed attempt
rm -f "${MERGED}.bed" "${MERGED}.bim" "${MERGED}.fam" "${MERGED}.log"
rm -f "$MISSNP"

# First merge attempt
log_info "  Merge attempt 1..."
set +e
"$PLINK19" \
  --bfile "$FIRST_BFILE" \
  --merge-list "$MERGE_LIST" \
  --make-bed \
  --out "$MERGED" \
  2>&1 | tee "${LOG_DIR}/cross_cohort_merge.log"
merge_exit=$?
set -e

# If merge failed with .missnp (3+ allele variants from strand-ambiguous SNPs)
if [[ $merge_exit -ne 0 ]] && [[ -f "$MISSNP" ]]; then
  N_MISSNP=$(wc -l < "$MISSNP")
  log_warn "  Merge failed: $N_MISSNP variants with 3+ alleles (strand-ambiguous SNPs)"
  log_info "  Excluding these from ALL cohorts and retrying..."

  # Create cleaned versions of each cohort
  CLEAN_MERGE_LIST="${PCS_OUT}/_cross_cohort_clean_merge_list.txt"
  : > "$CLEAN_MERGE_LIST"
  FIRST_CLEAN=""

  for cohort in "${COHORTS[@]}"; do
    bfile=$(get_harmonized_bfile "$cohort")
    clean="${PCS_OUT}/_${cohort}_clean"

    "$PLINK2" \
      --bfile "$bfile" \
      --exclude "$MISSNP" \
      --make-bed \
      --threads "$PLINK2_THREADS" \
      --out "$clean" \
      2>&1 | tee "${LOG_DIR}/cross_cohort_clean_${cohort}.log"

    n_clean=$(wc -l < "${clean}.bim")
    log_info "    $cohort: $n_clean variants after excluding missnp"

    if [[ -z "$FIRST_CLEAN" ]]; then
      FIRST_CLEAN="$clean"
    else
      echo "${clean}.bed ${clean}.bim ${clean}.fam" >> "$CLEAN_MERGE_LIST"
    fi
  done

  # Retry merge with cleaned files
  rm -f "${MERGED}.bed" "${MERGED}.bim" "${MERGED}.fam" "${MERGED}.log"
  log_info "  Merge attempt 2 (after excluding $N_MISSNP problematic SNPs)..."

  "$PLINK19" \
    --bfile "$FIRST_CLEAN" \
    --merge-list "$CLEAN_MERGE_LIST" \
    --make-bed \
    --out "$MERGED" \
    2>&1 | tee "${LOG_DIR}/cross_cohort_merge_retry.log"

  # Clean up per-cohort cleaned files
  for cohort in "${COHORTS[@]}"; do
    clean="${PCS_OUT}/_${cohort}_clean"
    rm -f "${clean}.bed" "${clean}.bim" "${clean}.fam" "${clean}.log" "${clean}.nosex"
  done
  rm -f "$CLEAN_MERGE_LIST"
fi

if [[ ! -f "${MERGED}.bed" ]]; then
  log_error "  Merge failed — check ${LOG_DIR}/cross_cohort_merge.log"
  exit 1
fi

N_SAMPLES=$(wc -l < "${MERGED}.fam")
N_VARIANTS=$(wc -l < "${MERGED}.bim")
log_info "  Merged: $N_SAMPLES samples x $N_VARIANTS variants"

# --- Step 3: Filter variants -------------------------------------------------
# Same QC as within-cohort PCs: autosomal, MAF > 5%, geno > 98%, HWE > 1e-6
# Exclude long-range LD regions
EXCLUDE_REGIONS="${PCS_OUT}/long_range_ld_hg19.txt"

# Create long-range LD exclusion file if it doesn't exist
if [[ ! -f "$EXCLUDE_REGIONS" ]]; then
  log_info "  Creating long-range LD exclusion regions file..."
  cat > "$EXCLUDE_REGIONS" << 'LDEOF'
1 48000000 52000000 1p11
2 86000000 100500000 2p11
2 134500000 138000000 2q21
2 183000000 190000000 2q33
3 47500000 50000000 3p21
3 83500000 87000000 3q12
5 44500000 50500000 5p13-14
5 98000000 100500000 5q15
5 129000000 132000000 5q31
5 136000000 139500000 5q33
6 25000000 35000000 6p22-MHC
6 57000000 64000000 6p12
6 140000000 142500000 6q25
7 55000000 66000000 7p12-q11
8 7000000 13000000 8p23
8 43000000 50000000 8p11
8 112000000 115000000 8q23
10 37000000 43000000 10p11
11 46000000 57000000 11p11-q13
11 87500000 90500000 11q14
12 33000000 40000000 12p11
12 109500000 112000000 12q23
20 32000000 34500000 20p11
LDEOF
fi

FILTERED="${PCS_OUT}/_cross_cohort_filtered"
"$PLINK2" \
  --bfile "$MERGED" \
  --autosome \
  --maf 0.05 \
  --geno 0.02 \
  --hwe 1e-6 0 \
  --exclude range "$EXCLUDE_REGIONS" \
  --make-bed \
  --threads "$PLINK2_THREADS" \
  --out "$FILTERED" \
  2>&1 | tee "${LOG_DIR}/cross_cohort_filter.log"

N_FILTERED=$(wc -l < "${FILTERED}.bim")
log_info "  After QC filtering: $N_FILTERED variants"

# --- Step 4: LD pruning ------------------------------------------------------
PRUNED="${PCS_OUT}/_cross_cohort_pruned"
"$PLINK19" \
  --bfile "$FILTERED" \
  --indep-pairwise 1000 50 0.05 \
  --out "$PRUNED" \
  2>&1 | tee "${LOG_DIR}/cross_cohort_prune.log"

N_PRUNED=$(wc -l < "${PRUNED}.prune.in")
log_info "  After LD pruning: $N_PRUNED independent SNPs"

# --- Step 5: Identify maximum unrelated subset across the pooled sample ------
# Captures both within-cohort relateds (TwinLife twins/sibs) and any cryptic
# cross-cohort relateds. PCs are trained on the unrelated set and projected
# onto everyone.
UNREL_PREFIX="${PCS_OUT}/_cross_cohort_unrel_set"
"$PLINK2" \
  --bfile "$FILTERED" \
  --extract "${PRUNED}.prune.in" \
  --king-cutoff "$KING_CUTOFF" \
  --threads "$PLINK2_THREADS" \
  --out "$UNREL_PREFIX" \
  2>&1 | tee "${LOG_DIR}/cross_cohort_king.log"

UNREL_KEEP="${UNREL_PREFIX}.king.cutoff.in.id"
N_TOTAL=$(wc -l < "${FILTERED}.fam")
N_UNREL=$(( $(wc -l < "$UNREL_KEEP") - 1 ))
N_RELATED=$(( N_TOTAL - N_UNREL ))
log_info "  Unrelated subset (KING < ${KING_CUTOFF}): ${N_UNREL} of ${N_TOTAL} pooled samples (${N_RELATED} relateds)"

# --- Step 6: Train PCs on unrelated subset (write allele weights + freqs) ----
TRAIN_PREFIX="${PCS_OUT}/_cross_cohort_train"
"$PLINK2" \
  --bfile "$FILTERED" \
  --extract "${PRUNED}.prune.in" \
  --keep "$UNREL_KEEP" \
  --freq counts \
  --pca "$N_PCS" allele-wts \
  --threads "$PLINK2_THREADS" \
  --out "$TRAIN_PREFIX" \
  2>&1 | tee "${LOG_DIR}/cross_cohort_train.log"

# --- Step 7: Project ALL pooled samples onto the trained axes ----------------
# SCORE*_SUM (not _AVG) keeps the eigenvec scale; variance-standardize uses
# training-set frequencies via --read-freq.
PROJ_PREFIX="${PCS_OUT}/_cross_cohort_proj"
# Resolve ID / A1 / PC1 column POSITIONS from the .eigenvec.allele header by
# NAME — PLINK2's column order varies by build (a-7.1LM inserts a
# PROVISIONAL_REF? column before A1), so hardcoding "2 5 ... 6-N" breaks
# --score with "mismatching allele codes / No valid variants". See
# 04_compute_pcs.sh for the full rationale.
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
  2>&1 | tee "${LOG_DIR}/cross_cohort_project.log"

if [[ ! -f "${PROJ_PREFIX}.sscore" ]]; then
  log_error "  Cross-cohort projection failed — check ${LOG_DIR}/cross_cohort_project.log"
  exit 1
fi

# --- Step 8: Reformat .sscore -> .eigenvec (rename SCORE*_SUM -> PC*) --------
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
    "${PROJ_PREFIX}.sscore" > "$CROSS_PC_OUT"

if [[ -f "$CROSS_PC_OUT" ]]; then
  N_PC_SAMPLES=$(( $(wc -l < "$CROSS_PC_OUT") - 1 ))
  log_info "  Cross-cohort PCs written: $N_PC_SAMPLES samples x $N_PCS PCs"
else
  log_error "  Reformat to eigenvec failed"
  exit 1
fi

# --- Step 9: Clean up intermediate files --------------------------------------
rm -f "${MERGED}.bed" "${MERGED}.bim" "${MERGED}.fam" "${MERGED}.log" "${MERGED}.nosex"
rm -f "${MERGED}-merge.missnp"
rm -f "${FILTERED}.bed" "${FILTERED}.bim" "${FILTERED}.fam" "${FILTERED}.log"
rm -f "${PRUNED}.prune.in" "${PRUNED}.prune.out" "${PRUNED}.log" "${PRUNED}.nosex"
rm -f "${UNREL_PREFIX}.king.cutoff.in.id" "${UNREL_PREFIX}.king.cutoff.out.id" \
      "${UNREL_PREFIX}.log"
rm -f "${TRAIN_PREFIX}.eigenvec" "${TRAIN_PREFIX}.eigenvec.allele" \
      "${TRAIN_PREFIX}.eigenval" "${TRAIN_PREFIX}.acount" "${TRAIN_PREFIX}.log"
rm -f "${PROJ_PREFIX}.sscore" "${PROJ_PREFIX}.log"
rm -f "$MERGE_LIST"

# --- Log decisions ------------------------------------------------------------
N_MISSNP_LOG="${N_MISSNP:-0}"
{
  echo ""
  echo "## Step 04b: Cross-Cohort PCs ($(date '+%Y-%m-%d %H:%M'))"
  echo ""
  echo "- PCs computed on POOLED sample across all $N_COHORTS cohorts"
  echo "- These capture between-cohort ancestry/batch differences"
  echo "- Used harmonized_shared genotypes (same SNP set in all cohorts)"
  if [[ "${N_MISSNP_LOG}" -gt 0 ]]; then
    echo "- Excluded $N_MISSNP_LOG strand-ambiguous SNPs (3+ alleles after merge)"
  fi
  echo "- QC: autosomal, MAF > 0.05, genotyping rate > 98%, HWE p > 1e-6"
  echo "- Excluded long-range LD regions (Price et al., 2008)"
  echo "- LD pruning: window 1000kb, step 50, r² < 0.05"
  echo "- $N_PRUNED independent SNPs used"
  echo "- Relatedness: trained on KING-unrelated subset (kinship < ${KING_CUTOFF},"
  echo "  ~3rd-degree, Manichaikul et al., 2010): ${N_UNREL} of ${N_TOTAL}"
  echo "  pooled samples (${N_RELATED} relateds projected onto trained axes)"
  echo "- $N_PC_SAMPLES samples projected, $N_PCS PCs"
  echo ""
} >> "$DECISION_LOG"

log_info "=== Step 04b complete ==="
