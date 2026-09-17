#!/usr/bin/env bash
# =============================================================================
# 03a_filter_imputation_quality.sh — Create per-cohort R² pass lists
# =============================================================================
# For each cohort, reads Minimac imputation info files and outputs lists of
# SNPs that pass the R² threshold. These lists are used by step 03 (Phase 2)
# to restrict the cross-cohort shared SNP set to well-imputed variants.
#
# Two thresholds are applied simultaneously:
#   - R2_THRESHOLD (default 0.3): lenient, used for main pipeline
#   - R2_THRESHOLD_STRICT (default 0.8): stringent, used for sensitivity
#
# OUTPUT (per cohort):
#   output/geno/{COHORT}_r2_pass_0.3.txt  — SNP IDs passing lenient threshold
#   output/geno/{COHORT}_r2_pass_0.8.txt  — SNP IDs passing strict threshold
#
# SNP ID FORMAT:
#   Info files use bare "chr:pos" (e.g., "22:16050435").
#   Pipeline harmonized bim files use "chr22:16050435".
#   This script prepends "chr" to match.
#   SOEP info files use "chr:pos:ref:alt" — the ":ref:alt" suffix is stripped.
#
# SHIP-Td SPECIAL HANDLING:
#   SHIP-Td is a merge of two genotyping batches (Omni 2.5 + GSA), each
#   imputed separately. We take min(R2_batch1, R2_batch2) per variant —
#   a variant must be well-imputed in BOTH batches to pass.
#
# USAGE:
#   bash scripts/03a_filter_imputation_quality.sh
#   (typically via `bash run_pipeline.sh --step 03a`)
#
# INPUTS:
#   ${PROJ_ROOT}/data/info_files/<COHORT>/*.info.gz       — Minimac info files
#   ${PROJ_ROOT}/data/info_files/TWINLIFE/TwinLife_imputation_r2.tsv.gz
#                                                         — TwinLife custom TSV
#
# OUTPUTS:
#   ${OUT_ROOT}/geno/<COHORT>_r2_pass_${R2_THRESHOLD}.txt        — lenient list
#   ${OUT_ROOT}/geno/<COHORT>_r2_pass_${R2_THRESHOLD_STRICT}.txt — strict list
#
# DEPENDS ON: step 00 (input checks); info files must be on Tardis already.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

if [[ "${R2_FILTER:-yes}" != "yes" ]]; then
  log_info "=== Step 03a: SKIPPED (R2_FILTER=no) ==="
  exit 0
fi

log_info "=== Step 03a: Imputation R² quality filtering ==="
log_info "  Lenient threshold:  R2 >= ${R2_THRESHOLD}"
log_info "  Strict threshold:   R2 >= ${R2_THRESHOLD_STRICT}"

DECISION_LOG="${PROJ_ROOT}/decision_log.md"
SUMMARY=""

# =============================================================================
# HELPER: Extract passing SNPs from a single Minimac info file
# Reads col 1 (SNP) and col 7 (Rsq), prepends "chr" to SNP ID.
# Writes to two output files (lenient + strict) in a single pass.
# Args: $1=info_file  $2=out_lenient  $3=out_strict  $4=id_transform
#   id_transform: "standard" for chr:pos, "soep" for chr:pos:ref:alt
# =============================================================================
extract_passing_snps() {
  local info_file="$1"
  local out_lenient="$2"
  local out_strict="$3"
  local id_transform="${4:-standard}"
  local r2_len="$R2_THRESHOLD"
  local r2_str="$R2_THRESHOLD_STRICT"

  if [[ "$id_transform" == "soep" ]]; then
    # SOEP: SNP ID = "22:16050822:G:A" → strip to "22:16050822" → prepend "chr"
    zcat "$info_file" | awk -F'\t' -v r2_len="$r2_len" -v r2_str="$r2_str" '
      NR == 1 { next }
      {
        # Strip :ref:alt suffix
        split($1, parts, ":")
        snp_id = "chr" parts[1] ":" parts[2]
        rsq = $7
        if (rsq + 0 >= r2_len + 0) print snp_id >> out_l
        if (rsq + 0 >= r2_str + 0) print snp_id >> out_s
      }
    ' out_l="$out_lenient" out_s="$out_strict"
  elif [[ "$id_transform" == "twinlife" ]]; then
    # TwinLife: 2-column format (SNP, Rsq), SNP = "22:16050435"
    zcat "$info_file" | awk -F'\t' -v r2_len="$r2_len" -v r2_str="$r2_str" '
      NR == 1 { next }
      {
        snp_id = "chr" $1
        rsq = $2
        if (rsq + 0 >= r2_len + 0) print snp_id >> out_l
        if (rsq + 0 >= r2_str + 0) print snp_id >> out_s
      }
    ' out_l="$out_lenient" out_s="$out_strict"
  else
    # Standard Minimac: col 1 = "22:16050435", col 7 = Rsq
    zcat "$info_file" | awk -F'\t' -v r2_len="$r2_len" -v r2_str="$r2_str" '
      NR == 1 { next }
      {
        snp_id = "chr" $1
        rsq = $7
        if (rsq + 0 >= r2_len + 0) print snp_id >> out_l
        if (rsq + 0 >= r2_str + 0) print snp_id >> out_s
      }
    ' out_l="$out_lenient" out_s="$out_strict"
  fi
}

# =============================================================================
# PROCESS EACH COHORT
# =============================================================================

for cohort in "${COHORTS[@]}"; do
  log_info ""
  log_info "--- Processing $cohort ---"

  out_lenient="${GENO_OUT}/${cohort}_r2_pass_${R2_THRESHOLD}.txt"
  out_strict="${GENO_OUT}/${cohort}_r2_pass_${R2_THRESHOLD_STRICT}.txt"

  # Skip if already done
  if [[ -f "$out_lenient" && -f "$out_strict" ]]; then
    n_len=$(wc -l < "$out_lenient" | tr -d ' ')
    n_str=$(wc -l < "$out_strict" | tr -d ' ')
    log_info "  SKIP: already done (lenient=$n_len, strict=$n_str)"
    SUMMARY="${SUMMARY}\n- $cohort: SKIPPED (lenient=$n_len, strict=$n_str)"
    continue
  fi

  # Clear output files
  : > "$out_lenient"
  : > "$out_strict"

  case "$cohort" in
    BASEII)
      if [[ ! -f "$INFO_BASEII" ]]; then
        log_error "  Info file not found: $INFO_BASEII"
        exit 1
      fi
      log_info "  Reading: $INFO_BASEII (single genome-wide file)"
      extract_passing_snps "$INFO_BASEII" "$out_lenient" "$out_strict" "standard"
      ;;

    SHIP0)
      log_info "  Reading: SHIP-0 per-chr info files"
      for CHR in $(seq 1 22); do
        info_file="${INFO_SHIP0_PATTERN//\{CHR\}/$CHR}"
        if [[ ! -f "$info_file" ]]; then
          log_warn "  Missing: $info_file"
          continue
        fi
        extract_passing_snps "$info_file" "$out_lenient" "$out_strict" "standard"
      done
      ;;

    SHIPTD)
      # Special handling: merge two batches, take min(R2) per variant
      log_info "  Reading: SHIP-Td batch 1 + batch 2 info files (min R2 strategy)"
      tmp_b1="${GENO_OUT}/_tmp_shiptd_b1.tsv"
      tmp_b2="${GENO_OUT}/_tmp_shiptd_b2.tsv"
      : > "$tmp_b1"
      : > "$tmp_b2"

      for CHR in $(seq 1 22); do
        b1_file="${INFO_SHIPTD_PATTERN//\{CHR\}/$CHR}"
        b2_file="${INFO_SHIPTD_B2_PATTERN//\{CHR\}/$CHR}"

        if [[ -f "$b1_file" ]]; then
          zcat "$b1_file" | awk -F'\t' 'NR > 1 { print $1 "\t" $7 }' >> "$tmp_b1"
        else
          log_warn "  Missing batch 1: $b1_file"
        fi

        if [[ -f "$b2_file" ]]; then
          zcat "$b2_file" | awk -F'\t' 'NR > 1 { print $1 "\t" $7 }' >> "$tmp_b2"
        else
          log_warn "  Missing batch 2: $b2_file"
        fi
      done

      # Join on SNP ID, take min(R2), then filter
      log_info "  Computing min(R2) across batches..."
      sort -t$'\t' -k1,1 "$tmp_b1" > "${tmp_b1}.sorted"
      sort -t$'\t' -k1,1 "$tmp_b2" > "${tmp_b2}.sorted"

      join -t$'\t' -j 1 "${tmp_b1}.sorted" "${tmp_b2}.sorted" \
        | awk -F'\t' -v r2_len="$R2_THRESHOLD" -v r2_str="$R2_THRESHOLD_STRICT" '
          {
            # $1=SNP, $2=R2_batch1, $3=R2_batch2
            min_r2 = ($2 + 0 < $3 + 0) ? $2 + 0 : $3 + 0
            snp_id = "chr" $1
            if (min_r2 >= r2_len + 0) print snp_id >> out_l
            if (min_r2 >= r2_str + 0) print snp_id >> out_s
          }
        ' out_l="$out_lenient" out_s="$out_strict"

      rm -f "$tmp_b1" "$tmp_b2" "${tmp_b1}.sorted" "${tmp_b2}.sorted"
      ;;

    SOEP)
      log_info "  Reading: SOEP per-chr info files (stripping :ref:alt suffix)"
      soep_missing_chrs=""
      for CHR in $(seq 1 22); do
        info_file="${INFO_SOEP_PATTERN//\{CHR\}/$CHR}"
        if [[ ! -f "$info_file" ]]; then
          log_warn "  Missing info: chr${CHR} — will fall back to bim (mildQC already applied)"
          soep_missing_chrs="${soep_missing_chrs} ${CHR}"
          continue
        fi
        extract_passing_snps "$info_file" "$out_lenient" "$out_strict" "soep"
      done

      # Fallback for missing chromosomes: the SOEP bim was built from mildQC
      # data (R2 > 0.1 already applied). Include those variants as "passing"
      # since they already cleared basic imputation quality filtering.
      if [[ -n "$soep_missing_chrs" ]]; then
        soep_bfile=$(get_harmonized_bfile "SOEP")
        if [[ -f "${soep_bfile}.bim" ]]; then
          for CHR in $soep_missing_chrs; do
            n_fallback=$(awk -v chr="$CHR" '$1 == chr { print "chr" $1 ":" $4 }' "${soep_bfile}.bim" \
              | tee -a "$out_lenient" \
              | tee -a "$out_strict" \
              | wc -l | tr -d ' ')
            log_info "  Fallback chr${CHR}: $n_fallback variants from bim (mildQC, R2>0.1)"
          done
        else
          log_warn "  No harmonized bim found for SOEP fallback — chr${soep_missing_chrs} excluded"
        fi
      fi
      ;;

    TWINLIFE)
      if [[ ! -f "$INFO_TWINLIFE" ]]; then
        log_error "  Info file not found: $INFO_TWINLIFE"
        log_error "  Expected: the TwinLife imputation-R² TSV at this path."
        log_error "  Re-run \`bash transfer_geno_to_tardis.sh\` from the Mac"
        log_error "  (or transfer_geno_to_tardis_workbench.sh from Workbench)"
        log_error "  to ship it from the TwinLife share."
        exit 1
      fi
      log_info "  Reading: $INFO_TWINLIFE (extracted from VCFs)"
      extract_passing_snps "$INFO_TWINLIFE" "$out_lenient" "$out_strict" "twinlife"
      ;;

    *)
      log_error "Unknown cohort: $cohort"
      exit 1
      ;;
  esac

  # Report counts
  n_len=$(wc -l < "$out_lenient" | tr -d ' ')
  n_str=$(wc -l < "$out_strict" | tr -d ' ')
  log_info "  Variants passing R2 >= ${R2_THRESHOLD}: $n_len"
  log_info "  Variants passing R2 >= ${R2_THRESHOLD_STRICT}: $n_str"
  SUMMARY="${SUMMARY}\n- $cohort: lenient(>=${R2_THRESHOLD})=$n_len, strict(>=${R2_THRESHOLD_STRICT})=$n_str"
done

# =============================================================================
# LOG DECISIONS
# =============================================================================
{
  echo ""
  echo "## Step 03a: Imputation R² Filtering ($(date '+%Y-%m-%d %H:%M'))"
  echo ""
  echo "- R2_FILTER=yes"
  echo "- Lenient threshold: R2 >= ${R2_THRESHOLD} (used for main scoring)"
  echo "- Strict threshold: R2 >= ${R2_THRESHOLD_STRICT} (used for sensitivity analysis)"
  echo "- SHIP-Td: min(R2_batch1, R2_batch2) strategy (conservative)"
  echo "- SOEP: stripped :ref:alt suffix from SNP IDs"
  echo "- TwinLife: R2 extracted from dosage VCF INFO field"
  echo "- All SNP IDs prefixed with 'chr' to match harmonized bim format"
  echo ""
  echo "### Per-cohort variant counts"
  echo -e "$SUMMARY"
  echo ""
} >> "$DECISION_LOG"

log_info ""
log_info "=== Step 03a complete ==="
