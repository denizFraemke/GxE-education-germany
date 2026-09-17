#!/usr/bin/env bash
# =============================================================================
# regression_check_blockdiag_primary.sh — Guard primary cell-level h² values
# =============================================================================
# Before running the secondary tracks (R1, G1, RG1) the cell-level primary
# (East/West × pre/post-1990) must reproduce its validated h²/SE values
# bit-for-bit. This script reads the existing h2_blockdiag.tsv and exits
# non-zero on any deviation > tolerance.
#
# Reference values are the validated constrained-fit (model = "standalone")
# h² and REML SE of the four primary cells. See Plan_deviations.md §7g (the
# block-diagonal + cohort FE diagnostic that established them) and §7h (the
# guard rationale).
#
# Exit codes:
#   0 — all four cells reproduce within tolerance
#   1 — at least one cell deviates beyond tolerance (the secondary tracks
#       should NOT be run; investigate first)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"

TSV="${FINAL_OUT}/h2_blockdiag.tsv"
if [[ ! -f "$TSV" ]]; then
    log_error "Regression check: h2_blockdiag.tsv missing at ${TSV}."
    log_error "  Run the cell-level primary first (step 06d, with no TRACK override)."
    exit 1
fi

# Numeric tolerance: 1e-3 on h² and on SE. GCTA's REML output is reported to
# 4-6 decimal places, so this catches any non-trivial drift while tolerating
# fp-noise from a rerun.
TOL="0.001"

# Reference: stratum, model, h2, SE (all four constrained-model cells)
REFERENCE=$(cat <<'EOF'
East_pre1990	standalone	0.282	0.050
East_post1990	standalone	0.040	0.428
West_pre1990	standalone	0.281	0.172
West_post1990	standalone	0.000	0.347
EOF
)

log_info "=== Regression check: cell-level primary h² ==="
log_info "  TSV:       ${TSV}"
log_info "  Tolerance: ±${TOL} on h² and on REML SE"
log_info ""

FAIL=0
while IFS=$'\t' read -r ref_stratum ref_model ref_h2 ref_se; do
    [[ -z "${ref_stratum:-}" ]] && continue
    row=$(awk -F'\t' -v s="$ref_stratum" -v m="$ref_model" '
        NR==1 { for (i=1;i<=NF;i++) col[$i]=i; next }
        $col["stratum"]==s && $col["model"]==m {
            printf "%s\t%s\n", $col["h2_SNP"], $col["h2_SNP_SE"]; exit
        }
    ' "$TSV")
    if [[ -z "$row" ]]; then
        log_error "  ${ref_stratum} / ${ref_model}: row missing in ${TSV}"
        FAIL=1
        continue
    fi
    obs_h2=$(echo "$row" | cut -f1)
    obs_se=$(echo "$row" | cut -f2)
    d_h2=$(awk -v a="$obs_h2" -v b="$ref_h2" 'BEGIN { d=a-b; if (d<0) d=-d; print d }')
    d_se=$(awk -v a="$obs_se" -v b="$ref_se" 'BEGIN { d=a-b; if (d<0) d=-d; print d }')
    bad_h2=$(awk -v d="$d_h2" -v t="$TOL" 'BEGIN { print (d > t ? "1" : "0") }')
    bad_se=$(awk -v d="$d_se" -v t="$TOL" 'BEGIN { print (d > t ? "1" : "0") }')
    if [[ "$bad_h2" == "1" || "$bad_se" == "1" ]]; then
        log_error "  FAIL ${ref_stratum}: h²=${obs_h2} (ref ${ref_h2}, Δ=${d_h2}); SE=${obs_se} (ref ${ref_se}, Δ=${d_se})"
        FAIL=1
    else
        log_info "  OK   ${ref_stratum}: h²=${obs_h2} (Δ=${d_h2}); SE=${obs_se} (Δ=${d_se})"
    fi
done <<< "$REFERENCE"

if (( FAIL == 1 )); then
    log_error ""
    log_error "Regression check FAILED. The validated cell-level primary has changed."
    log_error "  Do NOT run the R1/G1/RG1 secondary tracks until this is investigated."
    log_error "  Likely causes: change in covariate file, GRM cutoff, --grm-cutoff order,"
    log_error "  or genotype/phenotype data inputs."
    exit 1
fi

log_info ""
log_info "=== Regression check passed: primary track unchanged ==="
exit 0
