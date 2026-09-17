#!/usr/bin/env bash
# =============================================================================
# 05_fit_greml_main.sh — Submit per-stratum REML fits for the 4-cell main (SLURM)
# =============================================================================
# Writes a fits-list file with one row per REML fit, then submits a SLURM
# array. Only the 4 × 2 = 8 per-stratum standalone fits are emitted
# (constrained + --reml-no-constrain per cell).
#
# The originally-planned joint mGRM framework (full + nested) is NOT used:
# gcta --reml --mgrm segfaults on disjoint-sample GRMs. Cross-cell h²
# equality is instead tested via Cochran's Q in step 07. See
# Plan_deviations.md §7b for the full diagnosis.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

MAIN_OUT="${REML_OUT}/main"
mkdir -p "$MAIN_OUT"

log_info "=== Step 05: Submit GREML fits — 4-cell main (SLURM array) ==="

SLURM_SCRIPT="${PROJ_ROOT}/scripts/05_fit_greml_main.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# Sanity check: per-stratum _unrel GRMs must exist (step 04 must have run).
for s in East_pre1990 East_post1990 West_pre1990 West_post1990; do
    if [[ ! -f "${GRM_OUT}/grm_${s}_unrel.grm.bin" ]]; then
        log_error "Missing _unrel GRM for ${s} — run step 04 first."
        exit 1
    fi
done

# --- Match phen/qcovar/covar to pooled.fam by IID -------------------
# IID is the person-level key. GCTA happens to match by the (FID, IID)
# pair, so we copy pooled.fam's FID onto each input row before launching
# REML. (The Mac exporter writes FID=IID for everyone; pooled.fam
# inherits cohort-specific FID values from the source bfiles. Without
# this step the strict-pair match would drop the rows where the two
# disagree.) Idempotent — re-running is a no-op once the FIDs line up.
align_input_fids() {
    local pooled_fam="${GENO_OUT}/pooled.fam"
    if [[ ! -f "$pooled_fam" ]]; then
        log_error "pooled.fam not found: $pooled_fam (run step 02 first)"
        return 1
    fi

    local fam_index="${GENO_OUT}/_pooled_iid_to_fid.tsv"
    if [[ ! -s "$fam_index" ]]; then
        awk '{ print $2"\t"$1 }' "$pooled_fam" > "$fam_index"
        log_info "  Built IID → FID lookup: $fam_index ($(wc -l < "$fam_index") rows)"
    else
        log_info "  Reusing IID → FID lookup: $fam_index ($(wc -l < "$fam_index") rows)"
    fi

    local in_file in_name tmp_file
    for in_file in "$PHENO_FILE" "$QCOVAR_FILE" "$COVAR_FILE"; do
        if [[ ! -f "$in_file" ]]; then
            log_error "  Input file missing: $in_file"
            log_error "  Re-upload phenotype: transfer_phenotype_to_tardis.sh"
            return 1
        fi
        in_name="$(basename "$in_file")"
        tmp_file="${in_file}.fid_aligned.tmp"

        # awk: build IID → FID map from fam_index, then rewrite col 1 of
        # the input file with the looked-up FID. OFS set BEFORE the
        # assignment so $0 is rebuilt with tabs. Report aggregate counts
        # (total, changed, unmatched) to stderr.
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

# Build fits-list file (tab-separated: label, type, path, flags). Only
# per-stratum standalone fits — the joint mGRM framework was removed
# because gcta --reml --mgrm segfaults on disjoint-sample inputs
# (Plan_deviations.md §7b). Cross-cell heterogeneity is tested instead
# by step 07's Cochran's Q on the standalone estimates.
FITS_LIST="${LOG_DIR}/_step05_fits.txt"
: > "$FITS_LIST"

emit_fit() {
    # $1=label  $2=type (mgrm|grm)  $3=path  $4=extra-flags-or-empty
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$FITS_LIST"
}

# Per-stratum standalones (constrained + --reml-no-constrain).
for s in East_pre1990 East_post1990 West_pre1990 West_post1990; do
    emit_fit "standalone_${s}"             grm "${GRM_OUT}/grm_${s}_unrel" ""
    emit_fit "standalone_${s}_noconstrain" grm "${GRM_OUT}/grm_${s}_unrel" "--reml-no-constrain"
done

# Filter out fits whose .hsq already exists (resume mode).
PENDING_LINES=()
i=0
while IFS=$'\t' read -r label rest; do
    i=$((i + 1))
    if [[ -f "${MAIN_OUT}/${label}.hsq" ]]; then
        log_info "  SKIP: ${label}.hsq already exists"
        continue
    fi
    PENDING_LINES+=("$i")
done < "$FITS_LIST"

N=${#PENDING_LINES[@]}
if (( N == 0 )); then
    log_info ""
    log_info "All ${i} fits already complete — nothing to submit."
    log_info "=== Step 05 complete (no-op) ==="
    exit 0
fi

# Convert pending line numbers to a SLURM array index list.
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
    --export="ALL,FITS_LIST=${FITS_LIST},OUT_DIR=${MAIN_OUT}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

# --- Per-fit completion report ---
log_info ""
log_info "Fit outcomes:"
while IFS=$'\t' read -r label rest; do
    if [[ -f "${MAIN_OUT}/${label}.hsq" ]]; then
        log_info "  ${label}: OK"
    else
        log_warn "  ${label}: NO .hsq (REML may have failed to converge)"
    fi
done < "$FITS_LIST"

log_info ""
log_info "=== Step 05 complete (.hsq files in ${MAIN_OUT}/) ==="
