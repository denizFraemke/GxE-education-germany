#!/usr/bin/env bash
# =============================================================================
# 05c_fit_greml_per_cohort.sh — Submit per-cohort GREML fits (SLURM array)
# =============================================================================
# Per-cohort sensitivity track (Plan_deviations.md §7e closing subsection).
# Three groups of fits, all using the per-cohort unrel GRMs from step 04c:
#
#   1. Cohort-only      — 4 fits (one per cohort, full cohort N)
#   2. Cohort × region  — up to 7 fits (cohorts with both regions)
#   3. Cohort × stratum — up to 12 fits (cells with N >= 100)
#
# Each fit is submitted ± --reml-no-constrain → up to ~46 array tasks. All
# small (no fit > a few minutes), all run as one SLURM array.
#
# Subset filtering is done via GCTA's --keep at REML time: the cohort's
# unrel GRM stays whole, and we pass a subset keep file to fit on the
# (cohort × region) or (cohort × stratum) sub-sample. The keep files are
# generated below by intersecting the cohort's unrel-GRM IIDs with the
# region / stratum keep files from the main pipeline.
#
# FID alignment: phen/qcovar/covar are aligned to pooled.fam by IID
# (same logic as step 05). Idempotent — re-running is a no-op once aligned.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

COHORT_OUT="${REML_OUT}/per_cohort"
mkdir -p "$COHORT_OUT"

log_info "=== Step 05c: Submit per-cohort GREML fits (SLURM array) ==="

SLURM_SCRIPT="${PROJ_ROOT}/scripts/05c_fit_greml_per_cohort.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# ---- 0. Sanity: per-cohort unrel GRMs exist -----------------------------
for c in "${COHORTS_REML[@]}"; do
    if [[ ! -f "${GRM_OUT}/grm_cohort_${c}_unrel.grm.bin" ]]; then
        log_error "Missing unrel GRM for cohort ${c} — run step 04c first."
        exit 1
    fi
done

# ---- 1. Align phen/qcovar/covar FIDs to pooled.fam (same as step 05) ---
# IID is the person-level key; GCTA matches by (FID, IID). Mac exporter
# writes FID=IID for everyone; pooled.fam has cohort-specific FIDs.
# Idempotent.
align_input_fids() {
    local pooled_fam="${GENO_OUT}/pooled.fam"
    if [[ ! -f "$pooled_fam" ]]; then
        log_error "pooled.fam not found: $pooled_fam (run step 02 first)"
        return 1
    fi
    local fam_index="${GENO_OUT}/_pooled_iid_to_fid.tsv"
    if [[ ! -s "$fam_index" ]]; then
        awk '{ print $2"\t"$1 }' "$pooled_fam" > "$fam_index"
    fi

    local in_file in_name tmp_file
    for in_file in "$PHENO_FILE" "$QCOVAR_FILE" "$COVAR_FILE"; do
        if [[ ! -f "$in_file" ]]; then
            log_error "  Input file missing: $in_file"
            return 1
        fi
        in_name="$(basename "$in_file")"
        tmp_file="${in_file}.fid_aligned.tmp"
        awk -F'\t' -v OFS='\t' -v famidx="$fam_index" -v label="$in_name" '
          BEGIN {
            while ((getline line < famidx) > 0) {
              n = split(line, a, "\t")
              if (n >= 2) fid_of[a[1]] = a[2]
            }
            close(famidx)
          }
          {
            total++
            if ($2 in fid_of) {
              aligned = fid_of[$2]
              if ($1 != aligned) changed++
              $1 = aligned
            } else {
              unmatched++
            }
            print
          }
          END {
            printf "  %-22s total=%d  changed=%d  unmatched=%d\n",
                   label, total, (changed+0), (unmatched+0) > "/dev/stderr"
          }
        ' "$in_file" > "$tmp_file"
        mv "$tmp_file" "$in_file"
    done
}

log_info ""
log_info "Matching phen / qcovar / covar to pooled.fam by IID..."
align_input_fids
log_info "Done."

# ---- 2. Build cohort × region and cohort × stratum keep files ------------
# Each subset keep file is the intersection of:
#   (cohort's unrel-GRM IIDs)  ∩  (region or stratum analytic-sample IIDs)
# joined with pooled.fam for the FID column.
KEEP_DIR="${PROJ_ROOT}/data/keep"

# Minimum sub-cell N at which we actually fit. Below this the REML SE is
# wider than the parameter space (h² in [0, 1]) and the fit is uninformative.
MIN_SUBCELL_N=100

build_subset_keep() {
    # Args: cohort, subset_label, list_of_stratum_keep_files
    local cohort="$1"; shift
    local subset_label="$1"; shift
    local out="${KEEP_DIR}/cohort_${cohort}_${subset_label}.gcta.iids"

    local cohort_unrel="${GRM_OUT}/grm_cohort_${cohort}_unrel.grm.id"
    if [[ ! -f "$cohort_unrel" ]]; then
        echo "0"
        return 0
    fi

    local tmp_cohort_iids; tmp_cohort_iids="$(mktemp)"
    awk '{ print $2 }' "$cohort_unrel" | sort -u > "$tmp_cohort_iids"

    local tmp_subset_iids; tmp_subset_iids="$(mktemp)"
    : > "$tmp_subset_iids"
    for f in "$@"; do
        [[ -f "$f" ]] && awk '{ print $2 }' "$f" >> "$tmp_subset_iids"
    done
    sort -u "$tmp_subset_iids" -o "$tmp_subset_iids"

    local tmp_int; tmp_int="$(mktemp)"
    comm -12 "$tmp_cohort_iids" "$tmp_subset_iids" > "$tmp_int"

    local fam_index="${GENO_OUT}/_pooled_iid_to_fid.tsv"
    join -t $'\t' -1 1 -2 1 \
        <(sort -t $'\t' -k1,1 "$tmp_int") \
        <(sort -t $'\t' -k1,1 "$fam_index") \
        | awk -F'\t' 'BEGIN{OFS="\t"} { print $2, $1 }' \
        > "$out"

    rm -f "$tmp_cohort_iids" "$tmp_subset_iids" "$tmp_int"
    wc -l < "$out"
}

# Region keep files (union of pre- and post-1990 keep files for each region).
EAST_KEEPS=( "$(get_keep_file East_pre1990)" "$(get_keep_file East_post1990)" )
WEST_KEEPS=( "$(get_keep_file West_pre1990)" "$(get_keep_file West_post1990)" )

log_info ""
log_info "Building subset keep files (cohort × region, cohort × stratum)..."
for c in "${COHORTS_REML[@]}"; do
    # cohort × region
    for region in East West; do
        if [[ "$region" == "East" ]]; then
            n=$(build_subset_keep "$c" "East" "${EAST_KEEPS[@]}")
        else
            n=$(build_subset_keep "$c" "West" "${WEST_KEEPS[@]}")
        fi
        log_info "  ${c} × ${region}: N=${n}"
    done
    # cohort × stratum
    for s in "${STRATA_MAIN[@]}"; do
        n=$(build_subset_keep "$c" "$s" "$(get_keep_file "$s")")
        log_info "  ${c} × ${s}: N=${n}"
    done
done

# ---- 3. Build fits-list file --------------------------------------------
# Tab-separated: <label>\t<grm_path>\t<keep_path_or_NONE>\t<extra_flags>
FITS_LIST="${LOG_DIR}/_step05c_fits.txt"
: > "$FITS_LIST"

emit_fit() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$FITS_LIST"
}

for c in "${COHORTS_REML[@]}"; do
    grm="${GRM_OUT}/grm_cohort_${c}_unrel"

    # 1. Cohort-only (no --keep needed; uses the whole cohort unrel GRM).
    emit_fit "cohort_${c}"               "$grm" "NONE" ""
    emit_fit "cohort_${c}_noconstrain"   "$grm" "NONE" "--reml-no-constrain"

    # 2. Cohort × region (one per region with N >= MIN_SUBCELL_N).
    for region in East West; do
        keep="${KEEP_DIR}/cohort_${c}_${region}.gcta.iids"
        if [[ -f "$keep" ]]; then
            n=$(wc -l < "$keep")
            if (( n >= MIN_SUBCELL_N )); then
                emit_fit "cohort_${c}_${region}"              "$grm" "$keep" ""
                emit_fit "cohort_${c}_${region}_noconstrain"  "$grm" "$keep" "--reml-no-constrain"
            else
                log_info "  Skip ${c} × ${region}: N=${n} < ${MIN_SUBCELL_N}"
            fi
        fi
    done

    # 3. Cohort × stratum (one per stratum with N >= MIN_SUBCELL_N).
    for s in "${STRATA_MAIN[@]}"; do
        keep="${KEEP_DIR}/cohort_${c}_${s}.gcta.iids"
        if [[ -f "$keep" ]]; then
            n=$(wc -l < "$keep")
            if (( n >= MIN_SUBCELL_N )); then
                emit_fit "cohort_${c}_${s}"              "$grm" "$keep" ""
                emit_fit "cohort_${c}_${s}_noconstrain"  "$grm" "$keep" "--reml-no-constrain"
            else
                log_info "  Skip ${c} × ${s}: N=${n} < ${MIN_SUBCELL_N}"
            fi
        fi
    done
done

# Filter to pending (no existing .hsq).
PENDING_LINES=()
i=0
while IFS=$'\t' read -r label rest; do
    i=$((i + 1))
    if [[ -f "${COHORT_OUT}/${label}.hsq" ]]; then
        log_info "  SKIP: ${label}.hsq already exists"
        continue
    fi
    PENDING_LINES+=("$i")
done < "$FITS_LIST"

N=${#PENDING_LINES[@]}
if (( N == 0 )); then
    log_info ""
    log_info "All ${i} fits already complete — nothing to submit."
    log_info "=== Step 05c complete (no-op) ==="
    exit 0
fi

ARRAY_SPEC=$(IFS=','; echo "${PENDING_LINES[*]}")

log_info ""
log_info "Submitting SLURM array for ${N} REML fit(s) (of ${i} total):"
for ln in "${PENDING_LINES[@]}"; do
    log_info "  array task ${ln}: $(awk -v idx="$ln" 'NR == idx { print $1; exit }' "$FITS_LIST")"
done
log_info "Fits list: $FITS_LIST"

JOB_OUTPUT=$(sbatch \
    --wait \
    --array="${ARRAY_SPEC}" \
    --export="ALL,FITS_LIST=${FITS_LIST},OUT_DIR=${COHORT_OUT}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

log_info ""
log_info "Fit outcomes:"
while IFS=$'\t' read -r label rest; do
    if [[ -f "${COHORT_OUT}/${label}.hsq" ]]; then
        log_info "  ${label}: OK"
    else
        log_warn "  ${label}: NO .hsq (REML may have failed to converge)"
    fi
done < "$FITS_LIST"

log_info ""
log_info "=== Step 05c complete (.hsq files in ${COHORT_OUT}/) ==="
