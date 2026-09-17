#!/usr/bin/env bash
# Run Section B (Mobility) of the PGI workflow.
#
# Stages:
#   1. B_analysis.R  — load data, prepare samples (dat_cluster_mob),
#                      fit B0/B0_pre/B1/B2/B3 (ParEdu-controlled) + gender
#                      B4/B5/B6/B0g (family RE) + per-cell nonlinearity LRTs,
#                      LOO; save B_Mobility_Models.rds
#   2. B_paredu_moderation.R — parental education as a moderator of the
#                      PGI-Education association, by region (exploratory,
#                      deviations §17), plus the regional parental-education
#                      gradient and its leave-one-study-out check. Writes
#                      S7_ParEdu_*.csv, ParEdu_Gradient.csv and
#                      ParEdu_Gradient_LOO.csv into the run's manuscript_export/
#                      tree. Must run before the report stages: B_report.R
#                      consumes the gradient (Supplementary Data 5) and
#                      B_sensitivities_report.R consumes its leave-one-study-out
#                      table (Supplementary Data 6).
#   3. B_report.R    — build Supplementary_Data_5_Mobility_results.xlsx,
#                      B_Mobility_Plots.pdf, B_Mobility_LOO_Plots.pdf,
#                      and B_Mobility_Overview.md
#   4. B_sensitivities_analysis.R — S1-B/S2-B/S3-B; save B_Sensitivities_Models.rds
#   5. B_sensitivities_report.R   — build Supplementary_Data_6_Mobility_sensitivities_results.xlsx
#   6. B_appendix_analysis.R      — per-study E/E-RxG/F-linear; save B_Appendix_Models.rds
#   7. B_appendix_report.R        — build Supplementary_Data_7_Mobility_appendix_results.xlsx
#
# Output dir: ${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/B_Mobility/
# RUN_TS defaults to YYYYMMDD_HHMM; export RUN_TS before invoking to
# share a timestamp across sections.
#
# --dry-run: print the planned step sequence and exit without running R.

set -euo pipefail

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -z "${RUN_TS:-}" ]]; then
  export RUN_TS="$(date +%Y%m%d_%H%M)"
fi

DATA_ROOT="${DATA_ROOT:-/path/to/data_root}"
OUT_DIR="$DATA_ROOT/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_$RUN_TS/B_Mobility"
LOG_FILE="$OUT_DIR/B_run.log"

ANALYSIS_R="$SCRIPT_DIR/B_Mobility/B_analysis.R"
REPORT_R="$SCRIPT_DIR/B_Mobility/B_report.R"
SENS_ANALYSIS_R="$SCRIPT_DIR/B_Mobility/B_sensitivities_analysis.R"
SENS_REPORT_R="$SCRIPT_DIR/B_Mobility/B_sensitivities_report.R"
APP_ANALYSIS_R="$SCRIPT_DIR/B_Mobility/B_appendix_analysis.R"
APP_REPORT_R="$SCRIPT_DIR/B_Mobility/B_appendix_report.R"
PAREDU_R="$SCRIPT_DIR/B_Mobility/B_paredu_moderation.R"

echo "run_B.command"
echo "  RUN_TS    = $RUN_TS"
echo "  DATA_ROOT = $DATA_ROOT"
echo "  OUT_DIR   = $OUT_DIR"
echo "  Stage 1:  Rscript $ANALYSIS_R          (B0-B3 + gender B4/B5/B6/B0g)"
echo "  Stage 2:  Rscript $PAREDU_R         (ParEdu moderation + regional gradient + LOCO)"
echo "  Stage 3:  Rscript $REPORT_R"
echo "  Stage 4:  Rscript $SENS_ANALYSIS_R     (S1-B/S2-B/S3-B)"
echo "  Stage 5:  Rscript $SENS_REPORT_R"
echo "  Stage 6:  Rscript $APP_ANALYSIS_R      (per-study E/E-RxG/F-linear)"
echo "  Stage 7:  Rscript $APP_REPORT_R"
echo "  LOG       = $LOG_FILE"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[dry-run] no R execution, no writes to the data share."
  exit 0
fi

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "ERROR: DATA_ROOT not reachable: $DATA_ROOT" >&2
  echo "       Set DATA_ROOT to the mounted data share and retry." >&2
  exit 2
fi

# Create the output dir with a relative mkdir from inside the 04_GxE_Analysis
# parent. macOS `mkdir -p` walks the absolute path back to the share root,
# which on DFS-style automounts (where the share root is read-only and the
# project subtree is a separate writable submount) errors with
# "Read-only file system" at the boundary. A relative mkdir from inside
# the writable parent avoids the traversal.
GXE_DIR="$DATA_ROOT/04_GxE_Analysis"
REL_OUT_DIR="01_PGI_Analysis/runs/RUN_$RUN_TS/B_Mobility"
( cd "$GXE_DIR" && mkdir -p "$REL_OUT_DIR" )

# Stage 1: analysis
echo
echo "=== Stage 1: B_analysis.R ==="
Rscript "$ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 2: parental-education moderation + regional gradient + LOCO.
# Runs before the report stages because they consume its CSVs: B_report.R reads
# ParEdu_Gradient.csv for Supplementary Data 5, B_sensitivities_report.R reads
# ParEdu_Gradient_LOO.csv for Supplementary Data 6.
echo
echo "=== Stage 2: B_paredu_moderation.R ==="
Rscript "$PAREDU_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 3: report
echo
echo "=== Stage 3: B_report.R ==="
Rscript "$REPORT_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 4: sensitivities analysis (reads B_Mobility_Models.rds for gates)
echo
echo "=== Stage 4: B_sensitivities_analysis.R ==="
Rscript "$SENS_ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 5: sensitivities report
echo
echo "=== Stage 5: B_sensitivities_report.R ==="
Rscript "$SENS_REPORT_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 6: per-study appendix analysis
echo
echo "=== Stage 6: B_appendix_analysis.R ==="
Rscript "$APP_ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 7: per-study appendix report
echo
echo "=== Stage 7: B_appendix_report.R ==="
Rscript "$APP_REPORT_R" 2>&1 | tee -a "$LOG_FILE"

echo
echo "run_B.command: complete. Outputs under $OUT_DIR"
