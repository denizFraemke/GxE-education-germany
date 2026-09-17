#!/usr/bin/env bash
# Run Section A (Education) of the PGI workflow.
#
# Stages:
#   1. A_analysis.R               — load data, prepare samples, fit A1/A2/A3 and
#                                   the gender models; save A_Education_Models.rds
#   2. A_report.R                 — build Supplementary_Data_2_Education_results.xlsx,
#                                   A_Education_Plots.pdf, A_Education_LOO_Plots.pdf,
#                                   A_Education_Overview.md
#   3. A_sensitivities_analysis.R — S1-A/S2-A/S4/Method/S3-A; save
#                                   A_Sensitivities_Models.rds
#   4. A_sensitivities_report.R   — build
#                                   Supplementary_Data_3_Education_sensitivities_results.xlsx
#   5. A_appendix_analysis.R      — per-study D/F/F-linear; save A_Appendix_Models.rds
#   6. A_appendix_report.R        — build
#                                   Supplementary_Data_4_Education_appendix_results.xlsx
#
# Output dir: ${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/A_Education/
# RUN_TS defaults to YYYYMMDD_HHMM; export RUN_TS before invoking to
# share a timestamp across sections. The sensitivities and appendix stages
# read A_Education_Models.rds from this same run.
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
OUT_DIR="$DATA_ROOT/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_$RUN_TS/A_Education"
LOG_FILE="$OUT_DIR/A_run.log"

ANALYSIS_R="$SCRIPT_DIR/A_Education/A_analysis.R"
REPORT_R="$SCRIPT_DIR/A_Education/A_report.R"
SENS_ANALYSIS_R="$SCRIPT_DIR/A_Education/A_sensitivities_analysis.R"
SENS_REPORT_R="$SCRIPT_DIR/A_Education/A_sensitivities_report.R"
APP_ANALYSIS_R="$SCRIPT_DIR/A_Education/A_appendix_analysis.R"
APP_REPORT_R="$SCRIPT_DIR/A_Education/A_appendix_report.R"

echo "run_A.command"
echo "  RUN_TS    = $RUN_TS"
echo "  DATA_ROOT = $DATA_ROOT"
echo "  OUT_DIR   = $OUT_DIR"
echo "  Stage 1:  Rscript $ANALYSIS_R          (A1-A3 + gender A4/A5/A6/A0g/A-spec)"
echo "  Stage 2:  Rscript $REPORT_R"
echo "  Stage 3:  Rscript $SENS_ANALYSIS_R     (S1-A/S2-A/S4/Method/S3-A)"
echo "  Stage 4:  Rscript $SENS_REPORT_R"
echo "  Stage 5:  Rscript $APP_ANALYSIS_R      (per-study D/F/F-linear)"
echo "  Stage 6:  Rscript $APP_REPORT_R"
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
REL_OUT_DIR="01_PGI_Analysis/runs/RUN_$RUN_TS/A_Education"
( cd "$GXE_DIR" && mkdir -p "$REL_OUT_DIR" )

# Stage 1: analysis
echo
echo "=== Stage 1: A_analysis.R ==="
Rscript "$ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 2: report
echo
echo "=== Stage 2: A_report.R ==="
Rscript "$REPORT_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 3: sensitivities analysis (reads A_Education_Models.rds for gates)
echo
echo "=== Stage 3: A_sensitivities_analysis.R ==="
Rscript "$SENS_ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 4: sensitivities report
echo
echo "=== Stage 4: A_sensitivities_report.R ==="
Rscript "$SENS_REPORT_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 5: per-study appendix analysis
echo
echo "=== Stage 5: A_appendix_analysis.R ==="
Rscript "$APP_ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 6: per-study appendix report
echo
echo "=== Stage 6: A_appendix_report.R ==="
Rscript "$APP_REPORT_R" 2>&1 | tee -a "$LOG_FILE"

echo
echo "run_A.command: complete. Outputs under $OUT_DIR"
