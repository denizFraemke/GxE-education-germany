#!/usr/bin/env bash
# =============================================================================
# 05d_fit_greml_block_diagonal.sh — Submit REML on per-cell block-diagonal GRMs
# =============================================================================
# Per-cell GREML fits using the block-diagonal cell GRMs from step 04d.
# Cross-cohort entries of these GRMs are zero by construction, so REML
# estimates V_G entirely from within-cohort relatedness signal.
#
# Covariates: $BLOCKDIAG_COVAR_FILE = data/gender_cohort.covar
#   (gender_int + cohort_int categorical). The cohort indicator absorbs
#   between-cohort mean differences in education that would otherwise load
#   into V_e under the block-diagonal design and deflate h². See
#   Plan_deviations.md §7g for the diagnostic that established this.
#
# 4 cells × 2 variants (constrained + --reml-no-constrain) = 8 fits, all
# small. Reuses the FID-alignment block from step 05, extended to align the
# block-diagonal covar file as well (idempotent).
#
# NOTE on resume mode: a fit is skipped whenever its .hsq already exists, and
# the filename does not encode the covar file or GRM. Delete the output
# directory before re-running after a covariate or GRM change.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

# Parameterised stratum list, output subdir, covar file. Defaults reproduce
# the validated cell-level primary track (R×R cells, gender_cohort.covar,
# output/reml/blockdiag/). Secondary tracks (R1/G1/RG1) override these via
# the wrappers in scripts/run_blockdiag_*.sh.
if [[ -n "${STRATA_LIST_NAMES:-}" ]]; then
    read -r -a __STRATA <<< "${STRATA_LIST_NAMES}"
else
    __STRATA=("${STRATA_MAIN[@]}")
fi
BLOCKDIAG_OUT="${REML_OUT}/${BLOCKDIAG_OUT_SUBDIR:-blockdiag}"
__COVAR_FILE="${BLOCKDIAG_COVAR_FILE_OVERRIDE:-$BLOCKDIAG_COVAR_FILE}"
mkdir -p "$BLOCKDIAG_OUT"

log_info "=== Step 05d: Submit per-cell block-diagonal REML fits (SLURM array) ==="
log_info "  Stratum list:  ${__STRATA[*]}"
log_info "  Output dir:    ${BLOCKDIAG_OUT}"
log_info "  Covar file:    ${__COVAR_FILE}"

SLURM_SCRIPT="${PROJ_ROOT}/scripts/05d_fit_greml_block_diagonal.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

# Sanity: per-cell block-diagonal GRMs from step 04d must exist.
for s in "${__STRATA[@]}"; do
    if [[ ! -f "${GRM_OUT}/grm_blockdiag_${s}.grm.bin" ]]; then
        log_error "Missing block-diagonal GRM for ${s} — run step 04d first."
        exit 1
    fi
done

# ---- FID alignment (same as step 05; idempotent) -----------------------
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
    # 05d uses a track-specific covar file (default: gender_cohort.covar; R1
    # also uses gender_cohort.covar; G1 uses region_cohort.covar; RG1 uses
    # cohort_only.covar). We align the chosen $__COVAR_FILE here plus the
    # plain gender.covar so anything else in this run that picks it up has
    # aligned FIDs — idempotent.
    for in_file in "$PHENO_FILE" "$QCOVAR_FILE" "$COVAR_FILE" "$__COVAR_FILE"; do
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

# ---- Build fits list ----------------------------------------------------
# Tab-separated: <label>\t<grm_prefix>\t<extra_flags>
FITS_LIST="${LOG_DIR}/_step05d_fits.txt"
: > "$FITS_LIST"

emit_fit() {
    printf '%s\t%s\t%s\n' "$1" "$2" "${3:-}" >> "$FITS_LIST"
}

for s in "${__STRATA[@]}"; do
    grm="${GRM_OUT}/grm_blockdiag_${s}"
    emit_fit "blockdiag_${s}"             "$grm" ""
    emit_fit "blockdiag_${s}_noconstrain" "$grm" "--reml-no-constrain"
done

# Filter to pending (no existing .hsq).
PENDING_LINES=()
i=0
while IFS=$'\t' read -r label rest; do
    i=$((i + 1))
    if [[ -f "${BLOCKDIAG_OUT}/${label}.hsq" ]]; then
        log_info "  SKIP: ${label}.hsq already exists"
        continue
    fi
    PENDING_LINES+=("$i")
done < "$FITS_LIST"

N=${#PENDING_LINES[@]}
if (( N == 0 )); then
    log_info ""
    log_info "All ${i} fits already complete — nothing to submit."
    log_info "=== Step 05d complete (no-op) ==="
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
    --export="ALL,FITS_LIST=${FITS_LIST},OUT_DIR=${BLOCKDIAG_OUT},BLOCKDIAG_COVAR_FILE_OVERRIDE=${__COVAR_FILE}" \
    "$SLURM_SCRIPT" 2>&1) || {
    log_error "sbatch returned non-zero. Output:"
    log_error "$JOB_OUTPUT"
    exit 1
}
log_info "$JOB_OUTPUT"

log_info ""
log_info "Fit outcomes:"
while IFS=$'\t' read -r label rest; do
    if [[ -f "${BLOCKDIAG_OUT}/${label}.hsq" ]]; then
        log_info "  ${label}: OK"
    else
        log_warn "  ${label}: NO .hsq (REML may have failed to converge)"
    fi
done < "$FITS_LIST"

log_info ""
log_info "=== Step 05d complete (.hsq files in ${BLOCKDIAG_OUT}/) ==="
