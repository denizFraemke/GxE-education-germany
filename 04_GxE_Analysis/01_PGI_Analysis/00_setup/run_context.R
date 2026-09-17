# Run-context setup for the 01_PGI_Analysis workflow.
# Sourced by every section's analysis and report script.
#
# Provides:
#   get_run_ts()             RUN_TS for this invocation (cached in env)
#   data_root()              DATA_ROOT (env override + default below)
#   gxe_dir()                the 04_GxE_Analysis/ root on the data share
#   output_dir(section_id)   per-section run dir under
#                            01_PGI_Analysis/runs/RUN_TS/; created if missing
#
# Outputs live under
#   ${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/<section>/

get_run_ts <- function() {
  ts <- Sys.getenv("RUN_TS", unset = "")
  if (!nzchar(ts)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M")
    Sys.setenv(RUN_TS = ts)
  }
  ts
}

data_root <- function() {
  Sys.getenv(
    "DATA_ROOT",
    unset = "/path/to/data_root"
  )
}

gxe_dir <- function() {
  file.path(data_root(), "04_GxE_Analysis")
}

output_dir <- function(section_id) {
  d <- file.path(gxe_dir(), "01_PGI_Analysis", "runs",
                 paste0("RUN_", get_run_ts()), section_id)
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  d
}
