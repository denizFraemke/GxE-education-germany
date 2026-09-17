#!/usr/bin/env bash
# Run the 00_Descriptives module of the PGI workflow.
#
# Stages:
#   1. Descriptives_analysis.R — load data, build the full-family and dedup
#                                analytic samples, compute aggregate descriptive
#                                tables, save Descriptives_Tables.rds
#   2. Descriptives_report.R   — build
#                                Supplementary_Data_1_Descriptives_results.xlsx
#                                and Descriptives_Overview.md
#
# Output dir: ${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/00_Descriptives/
# RUN_TS defaults to YYYYMMDD_HHMM; export RUN_TS before invoking to share a
# timestamp across sections.
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
OUT_DIR="$DATA_ROOT/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_$RUN_TS/00_Descriptives"
LOG_FILE="$OUT_DIR/Descriptives_run.log"

ANALYSIS_R="$SCRIPT_DIR/00_Descriptives/Descriptives_analysis.R"
REPORT_R="$SCRIPT_DIR/00_Descriptives/Descriptives_report.R"

echo "run_Descriptives.command"
echo "  RUN_TS    = $RUN_TS"
echo "  DATA_ROOT = $DATA_ROOT"
echo "  OUT_DIR   = $OUT_DIR"
echo "  Stage 1:  Rscript $ANALYSIS_R"
echo "  Stage 2:  Rscript $REPORT_R"
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

# Create the output dir with a relative mkdir from inside the writable
# 04_GxE_Analysis parent — avoids the DFS-submount "Read-only file system"
# boundary error that an absolute mkdir -p hits on a DFS-style automount.
GXE_DIR="$DATA_ROOT/04_GxE_Analysis"
REL_OUT_DIR="01_PGI_Analysis/runs/RUN_$RUN_TS/00_Descriptives"
( cd "$GXE_DIR" && mkdir -p "$REL_OUT_DIR" )

# Stage 1: analysis
echo
echo "=== Stage 1: Descriptives_analysis.R ==="
Rscript "$ANALYSIS_R" 2>&1 | tee -a "$LOG_FILE"

# Stage 2: report
echo
echo "=== Stage 2: Descriptives_report.R ==="
Rscript "$REPORT_R" 2>&1 | tee -a "$LOG_FILE"

echo
echo "run_Descriptives.command: complete. Outputs under $OUT_DIR"
