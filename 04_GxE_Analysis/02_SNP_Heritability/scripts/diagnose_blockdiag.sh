#!/usr/bin/env bash
# =============================================================================
# diagnose_blockdiag.sh — Diagnose the block-diagonal h² collapse
# =============================================================================
# Validates the block-diagonal cell-level specification (steps 03d/04d/05d)
# against the per-cohort × stratum fit from step 05c, and establishes whether
# a cohort fixed effect is required in --covar. Adding ~1,000 non-SHIP
# individuals block-diagonally to the ~6,000 SHIP individuals of the
# East × pre-1990 cell should not move the joint h² away from the SHIP-only
# estimate unless the model is misspecified.
#
# Two checks:
#
#   (1) SINGLE-BLOCK REPLICATION
#       Combine SHIP × East-pre1990 alone via R/build_block_diagonal_grm.R
#       (a 1-block "block-diagonal" GRM), then fit REML on it using the
#       same phen/qcovar/covar setup as step 05d. The result should
#       match SHIP × East-pre1990 from step 05c (h² ≈ 0.287). If it
#       doesn't, the block-diagonal *implementation* is buggy
#       (GRM writing, .grm.id ordering, .grm.N.bin, or covariate
#       alignment).
#
#   (2) COHORT FIXED-EFFECTS REFIT
#       Build a gender_cohort.covar that adds a categorical cohort
#       indicator to gender.covar. Fit the four full block-diagonal cell
#       GRMs with this covar as well as with the plain gender.covar. If
#       the East × pre-1990 h² rises substantially (toward 0.20–0.30)
#       when the cohort indicator is added, the gender-only specification
#       is misspecified — it absorbs between-cohort education-mean
#       differences into V_e because no cohort fixed effect is in the
#       model. If h² stays low (~0.05–0.15), the single-V_G joint REML
#       across heterogeneous cohorts is intrinsically too restrictive.
#
# Output:
#   output/reml/blockdiag_diagnostic/*.hsq   — 10 .hsq files (5 fits × 2 variants)
#   output/final/h2_blockdiag_diagnostic.tsv — tidy table
#   output/final/h2_blockdiag_diagnostic.md  — human-readable report
#
# The existing block-diagonal track output (output/reml/blockdiag/,
# output/final/h2_blockdiag*) is NOT modified.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"
ensure_dirs

log_info "=== diagnose_blockdiag: block-diagonal validation checks ==="

DIAG_OUT="${REML_OUT}/blockdiag_diagnostic"
mkdir -p "$DIAG_OUT"

SLURM_SCRIPT="${PROJ_ROOT}/scripts/diagnose_blockdiag.slurm"
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    log_error "SLURM worker missing: $SLURM_SCRIPT"
    exit 1
fi

COMBINER="${PROJ_ROOT}/R/build_block_diagonal_grm.R"
if [[ ! -f "$COMBINER" ]]; then
    log_error "Block-diagonal combiner missing: $COMBINER"
    exit 1
fi

# --------------------------------------------------------------------
# 1. Generate gender_cohort.covar (gender + categorical cohort indicator)
# --------------------------------------------------------------------
# Cohort coding: 1=BASEII, 2=SHIP (=SHIP-0 ∪ SHIP-Td), 3=SOEP, 4=TWINLIFE.
# Source-of-truth for IID→cohort is the per-cohort harmonized fam files
# in ${GENO_OUT}/<COHORT>_harmonized_r2_0.9.fam.
COVAR_COHORT_FILE="${PROJ_ROOT}/data/gender_cohort.covar"

log_info ""
log_info "Building IID → cohort lookup and writing $COVAR_COHORT_FILE..."

Rscript - <<RSCRIPT_EOF
geno_dir <- "${GENO_OUT}"
cohorts  <- c(BASEII = 1L, SHIP0 = 2L, SHIPTD = 2L, SOEP = 3L, TWINLIFE = 4L)

iid_cohort <- list()  # IID -> cohort_int; first-cohort-wins
for (c in names(cohorts)) {
  fam <- file.path(geno_dir, sprintf("%s_harmonized_r2_0.9.fam", c))
  if (!file.exists(fam)) stop("Missing fam: ", fam)
  iids <- read.table(fam, header = FALSE, stringsAsFactors = FALSE,
                     colClasses = c("character", "character",
                                    "NULL", "NULL", "NULL", "NULL"))[, 2]
  for (i in iids) if (is.null(iid_cohort[[i]])) iid_cohort[[i]] <- unname(cohorts[c])
}
cat(sprintf("  IIDs in cohort lookup: %d\n", length(iid_cohort)))

covar <- read.table("${COVAR_FILE}", header = FALSE, stringsAsFactors = FALSE,
                    colClasses = c("character","character","integer"))
names(covar) <- c("FID", "IID", "gender_int")
cat(sprintf("  IIDs in gender.covar: %d\n", nrow(covar)))

covar\$cohort_int <- unlist(iid_cohort[covar\$IID])
n_missing <- sum(is.na(covar\$cohort_int))
if (n_missing > 0L) {
  stop(sprintf("%d IIDs in gender.covar have no cohort assignment", n_missing))
}

write.table(covar, "${COVAR_COHORT_FILE}",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("  Wrote: ${COVAR_COHORT_FILE} (N=%d, cols: FID IID gender_int cohort_int)\n",
            nrow(covar)))
cat("  Cohort distribution in gender_cohort.covar:\n")
print(table(covar\$cohort_int, dnn = "cohort"))
RSCRIPT_EOF

if [[ ! -f "$COVAR_COHORT_FILE" ]]; then
    log_error "Failed to write $COVAR_COHORT_FILE"
    exit 1
fi
log_info "Done."

# --------------------------------------------------------------------
# 2. Single-block replication GRM: SHIP × East_pre1990 alone via R combiner
# --------------------------------------------------------------------
SOLO_BLOCKCOMP="${GRM_OUT}/grm_blockcomp_SHIP_East_pre1990"
SOLO_BLOCKDIAG="${GRM_OUT}/grm_blockdiag_SHIP_East_pre1990_solo"

log_info ""
log_info "Building single-block-via-combiner GRM (SHIP × East_pre1990 alone)..."
if [[ ! -f "${SOLO_BLOCKCOMP}.grm.bin" ]]; then
    log_error "Missing block component: ${SOLO_BLOCKCOMP}.grm.bin"
    log_error "Run step 03d first."
    exit 1
fi

if [[ -f "${SOLO_BLOCKDIAG}.grm.bin" ]]; then
    n=$(wc -l < "${SOLO_BLOCKDIAG}.grm.id")
    log_info "  SKIP: solo block-diagonal already built (N=${n})"
else
    Rscript "$COMBINER" "$SOLO_BLOCKDIAG" "$SOLO_BLOCKCOMP"
    if [[ -f "${SOLO_BLOCKDIAG}.grm.bin" ]]; then
        n=$(wc -l < "${SOLO_BLOCKDIAG}.grm.id")
        log_info "  Wrote: ${SOLO_BLOCKDIAG}.grm.{bin,id,N.bin} (N=${n})"
    else
        log_error "Failed to build solo block-diagonal GRM"
        exit 1
    fi
fi

# --------------------------------------------------------------------
# 3. Build FITS_LIST and submit SLURM array
# --------------------------------------------------------------------
# Format (tab-separated): label  grm_prefix  covar_file  extra_flags
#
# Fit catalogue:
#   - solo blockdiag for SHIP × East_pre (1-block validation): 2 fits
#   - cohort-dummy refit for each of 4 cells: 4 × 2 = 8 fits
# Total: 10.
#
# The gender-only block-diagonal fits are not repeated here — they are
# already in output/reml/blockdiag/ from step 05d — and neither is
# SHIP × East-pre via the cohort GRM + --keep, which is in
# output/reml/per_cohort/ from step 05c. The diagnostic report reads all
# four sources side-by-side.
FITS_LIST="${LOG_DIR}/_diagnose_blockdiag_fits.txt"
: > "$FITS_LIST"

emit_fit() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$FITS_LIST"
}

# (1) Single-block replication
emit_fit "solo_SHIP_East_pre1990"             "$SOLO_BLOCKDIAG" "$COVAR_FILE" ""
emit_fit "solo_SHIP_East_pre1990_noconstrain" "$SOLO_BLOCKDIAG" "$COVAR_FILE" "--reml-no-constrain"

# (2) Cohort-dummy refit of the 4 existing block-diagonal cells
for s in "${STRATA_MAIN[@]}"; do
    grm="${GRM_OUT}/grm_blockdiag_${s}"
    if [[ ! -f "${grm}.grm.bin" ]]; then
        log_warn "  Missing block-diagonal GRM for ${s} (run step 04d first) — skipping cohort-dummy refit"
        continue
    fi
    emit_fit "blockdiag_${s}_cohortcovar"             "$grm" "$COVAR_COHORT_FILE" ""
    emit_fit "blockdiag_${s}_cohortcovar_noconstrain" "$grm" "$COVAR_COHORT_FILE" "--reml-no-constrain"
done

PENDING_LINES=()
i=0
while IFS=$'\t' read -r label rest; do
    i=$((i + 1))
    if [[ -f "${DIAG_OUT}/${label}.hsq" ]]; then
        log_info "  SKIP: ${label}.hsq already exists"
        continue
    fi
    PENDING_LINES+=("$i")
done < "$FITS_LIST"

N=${#PENDING_LINES[@]}
if (( N == 0 )); then
    log_info ""
    log_info "All ${i} diagnostic fits already complete — running report only."
else
    ARRAY_SPEC=$(IFS=','; echo "${PENDING_LINES[*]}")
    log_info ""
    log_info "Submitting SLURM array for ${N} diagnostic fit(s) (of ${i} total):"
    for ln in "${PENDING_LINES[@]}"; do
        log_info "  array task ${ln}: $(awk -v idx="$ln" 'NR == idx { print $1; exit }' "$FITS_LIST")"
    done
    log_info "Fits list: $FITS_LIST"

    JOB_OUTPUT=$(sbatch \
        --wait \
        --array="${ARRAY_SPEC}" \
        --export="ALL,FITS_LIST=${FITS_LIST},OUT_DIR=${DIAG_OUT}" \
        "$SLURM_SCRIPT" 2>&1) || {
        log_error "sbatch returned non-zero. Output:"
        log_error "$JOB_OUTPUT"
        exit 1
    }
    log_info "$JOB_OUTPUT"

    log_info ""
    log_info "Fit outcomes:"
    while IFS=$'\t' read -r label rest; do
        if [[ -f "${DIAG_OUT}/${label}.hsq" ]]; then
            log_info "  ${label}: OK"
        else
            log_warn "  ${label}: NO .hsq"
        fi
    done < "$FITS_LIST"
fi

# --------------------------------------------------------------------
# 4. Run the report assembler
# --------------------------------------------------------------------
log_info ""
log_info "Running diagnostic report assembler..."
Rscript "${PROJ_ROOT}/scripts/diagnose_blockdiag_report.R"

log_info ""
log_info "=== diagnose_blockdiag complete ==="
log_info "Report: ${FINAL_OUT}/h2_blockdiag_diagnostic.md"
log_info "Table:  ${FINAL_OUT}/h2_blockdiag_diagnostic.tsv"
