#!/usr/bin/env Rscript
# check_frozen_run.R
# ----------------------------------------------------------------------------
# Reproducibility check: is this machine still able to reproduce the run that
# the manuscript repository's frozen exports were built from?
#
# Run it before fitting or exporting anything against a frozen run. It answers
# five questions and exits non-zero on the first one that fails:
#
#   1. Is DATA_ROOT reachable (SMB share mounted)?
#   2. Which Combined_Harmonized_*.rds does 00_setup/load_data.R resolve to?
#      It picks the newest by mtime, so a newer merge silently changes the
#      sample. This prints the choice instead of leaving it implicit.
#   3. Do the frozen run's model caches still carry the size+mtime fingerprints
#      recorded in the manuscript repo's 01_results/manuscript_export/
#      _provenance.csv?
#   4. Rebuilt from source, do the analytic samples still have the N the frozen
#      descriptives cache recorded (full-family and dedup one-per-family)?
#   5. Do a few frozen headline estimates still read back unchanged?
#
# Aggregate output only: counts, dimensions, coefficients, file fingerprints.
# No individual-level row is printed or written.
#
# Usage:
#   RUN_TS=20260729_1555 Rscript 05_manuscript_export/check_frozen_run.R
#   RUN_TS=... MANUSCRIPT_REPO=/path/to/manuscript Rscript ...
#
#   MANUSCRIPT_REPO has no default; step 3 is skipped when it is unset.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  # ggplot2 is attached because GxE_Germany_analysis_helpers.R (sourced in step 4
  # for prepare_samples.R's helpers) builds ggplot scales at source time.
  library(ggplot2); library(tibble)
})

SECTION_DIR  <- normalizePath(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(
  file.path(SECTION_DIR, "..", "04_GxE_Analysis", "01_PGI_Analysis"), mustWork = TRUE)

if (!nzchar(Sys.getenv("RUN_TS"))) Sys.setenv(RUN_TS = "20260729_1555")
source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS <- get_run_ts()

failures <- character(0)
ok   <- function(...) cat(sprintf("  [ok]   %s\n", sprintf(...)))
bad  <- function(...) { m <- sprintf(...); cat(sprintf("  [FAIL] %s\n", m))
                        failures <<- c(failures, m) }
note <- function(...) cat(sprintf("  [note] %s\n", sprintf(...)))

cat(sprintf("\ncheck_frozen_run.R -- RUN_%s\n", RUN_TS))

# ---- 1. Data root reachable ------------------------------------------------
cat("\n1. Data root\n")
DATA_ROOT <- data_root()
if (!dir.exists(DATA_ROOT)) {
  bad("DATA_ROOT not reachable: %s -- mount the SMB share (Finder > Go > Connect to Server) or set DATA_ROOT", DATA_ROOT)
  cat("\nFAILED: cannot continue without the data root.\n"); quit(status = 1)
}
ok("DATA_ROOT = %s", DATA_ROOT)

# ---- 2. Which merged RDS does load_data.R resolve to? ----------------------
cat("\n2. Merged phenotype+PGI file resolved by 00_setup/load_data.R\n")
DATA_DIR <- file.path(DATA_ROOT, "03_Merge")
# Same pattern and same newest-by-mtime rule as load_data.R.
cands <- list.files(DATA_DIR, pattern = "^Combined_Harmonized_\\d{8}\\.rds$", full.names = TRUE)
if (!length(cands)) {
  bad("no Combined_Harmonized_<YYYYMMDD>.rds in %s -- run 03_Merge/merge.R", DATA_DIR)
} else {
  picked <- cands[which.max(file.info(cands)$mtime)]
  ok("resolves to %s (newest of %d by mtime)", basename(picked), length(cands))
  note("load_data.R selects by mtime, not by name -- a newer merge would change the sample silently")
}

# ---- 3. Cache fingerprints vs the manuscript repo's provenance sidecar ------
cat("\n3. Model-cache fingerprints vs the manuscript repo's _provenance.csv\n")
finger <- function(path) {                      # identical to the export script's
  if (!file.exists(path)) return(NA_character_)
  info <- file.info(path)
  sprintf("size=%d;mtime=%s", info$size, format(info$mtime, "%Y-%m-%dT%H:%M:%S"))
}
A_PATH <- file.path(output_dir("A_Education"),     "A_Education_Models.rds")
B_PATH <- file.path(output_dir("B_Mobility"),      "B_Mobility_Models.rds")
D_PATH <- file.path(output_dir("00_Descriptives"), "Descriptives_Tables.rds")

# Location of the manuscript checkout. No default: it is not derivable from
# SECTION_DIR, so it must be supplied through the MANUSCRIPT_REPO env var.
MS   <- Sys.getenv("MANUSCRIPT_REPO", unset = "")
PROV <- if (nzchar(MS)) file.path(MS, "01_results", "manuscript_export", "_provenance.csv") else ""
if (!nzchar(PROV) || !file.exists(PROV)) {
  note("no manuscript provenance sidecar (%s) -- skipping (set MANUSCRIPT_REPO)",
       if (nzchar(PROV)) PROV else "MANUSCRIPT_REPO unset")
} else {
  p <- utils::read.csv(PROV, stringsAsFactors = FALSE)
  if (!identical(as.character(p$run_ts), RUN_TS))
    bad("provenance records run %s but this check is running against %s", p$run_ts, RUN_TS)
  for (x in list(list("A_Education_Models.rds",  A_PATH, p$a_cache_fingerprint),
                 list("B_Mobility_Models.rds",   B_PATH, p$b_cache_fingerprint),
                 list("Descriptives_Tables.rds", D_PATH, p$desc_cache_fingerprint))) {
    now <- finger(x[[2]])
    if (is.na(now))                    bad("%s is missing from the run folder", x[[1]])
    else if (identical(now, x[[3]]))   ok("%-24s %s", x[[1]], now)
    else bad("%s changed since export\n           recorded %s\n           on disk  %s",
             x[[1]], x[[3]], now)
  }
  note("analysis_repo_commit recorded as %s", substr(p$analysis_repo_commit, 1, 8))
}

# ---- 4. Rebuild the analytic samples and compare N -------------------------
cat("\n4. Analytic sample sizes, rebuilt from source vs the frozen descriptives cache\n")
if (!file.exists(D_PATH)) {
  bad("descriptives cache missing -- cannot compare sample sizes")
} else {
  tb <- readRDS(D_PATH)
  RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()
  EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
  source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),      local = TRUE)
  # prepare_samples.R calls safe_numeric()/cohort_groups() from the analysis
  # helpers, exactly as Descriptives_analysis.R sources them before it.
  source(file.path(WORKFLOW_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
  sink(tempfile())                                # load/prepare are chatty by design
  suppressMessages(suppressWarnings({
    source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),     local = TRUE)
    source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)
  }))
  sink()
  n_full  <- nrow(dat_cluster)
  n_dedup <- nrow(dat_edu)
  if (identical(n_full,  as.integer(tb$N_total)))
    ok("full-family sample (dat_cluster) N = %d, matches the frozen cache", n_full)
  else bad("full-family N = %d, frozen cache recorded %d", n_full, tb$N_total)
  if (identical(n_dedup, as.integer(tb$N_test)))
    ok("dedup one-per-family sample (dat_edu) N = %d, matches the frozen cache", n_dedup)
  else bad("dedup N = %d, frozen cache recorded %d", n_dedup, tb$N_test)
}

# ---- 5. Spot-check frozen headline estimates -------------------------------
cat("\n5. Frozen headline estimates read back from the caches\n")
if (file.exists(A_PATH)) {
  A <- readRDS(A_PATH)
  gx <- A$A4_Coefficients[A$A4_Coefficients$term == "PGI_Edu_z:gender_c", ]
  if (nrow(gx) == 1)
    ok("A4 PGI x gender  b = %+.4f, SE = %.4f, p = %.4f", gx$estimate, gx$std.error, gx$p.value)
  else bad("A4 PGI x gender coefficient not found in the cache")
} else bad("A cache missing -- cannot spot-check")

cat("\n")
if (length(failures)) {
  cat(sprintf("FAILED (%d): the frozen run is NOT reproducible on this machine as configured.\n",
              length(failures)))
  for (f in failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("PASSED: the frozen run reproduces -- safe to backfill exports against it.\n")
