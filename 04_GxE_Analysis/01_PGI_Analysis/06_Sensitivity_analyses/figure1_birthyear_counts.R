# figure1_birthyear_counts.R
# One aggregate export, which fits no model: per-study x birth-year Ns for
# Figure 1.
#
# It does not refit the frozen run. It reads the merged RDS and rebuilds the two
# analytic samples through the same prepare_samples.R the frozen run used.
#
# Usage:
#   RUN_TS=<ts> DATA_ROOT=<...> \
#     OUT_DIR=<...> Rscript 06_Sensitivity_analyses/figure1_birthyear_counts.R

# This script sources 00_setup/prepare_samples.R, which builds the run's
# descriptive plots on the way to the analytic samples, so it needs the same
# library set A_analysis.R and B_analysis.R load.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(tibble)
  library(mgcv)
  library(openxlsx)
})

# ---------------------------------------------------------------------------
# Locations
# ---------------------------------------------------------------------------

.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) {
      return(normalizePath(f$ofile, mustWork = TRUE))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fa) > 0 && nzchar(fa[1])) {
    return(normalizePath(fa[1], mustWork = TRUE))
  }
  stop("Cannot determine figure1_birthyear_counts.R location.")
}
SECTION_DIR <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
REPO_ROOT <- normalizePath(file.path(WORKFLOW_DIR, "..", ".."), mustWork = TRUE)

# The source commit stamped onto every export (rule: every export carries it).
source_commit <- function() {
  out <- suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ))
  if (length(out) == 1L && nzchar(out)) out else NA_character_
}

step_counts <- function(out_dir) {
  # Rebuild the two analytic samples exactly as the frozen run did, through the
  # same setup A_analysis.R and B_analysis.R use, then count. Counts only: no
  # value belonging to any individual leaves this function.
  source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
  source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
  source(file.path(WORKFLOW_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

  RESULTS <- list()
  PLOTS <- list()
  PRIMARY_RESULTS <- list()
  DATA_ROOT <- data_root()
  DATA_DIR <- file.path(DATA_ROOT, "03_Merge")
  EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
  source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"), local = TRUE)
  source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

  # dat_cluster / dat_cluster_mob are the primary random-effects analytic samples
  # the frozen focal models were fitted on - the Ns the manuscript reports.
  stopifnot(exists("dat_cluster"), exists("dat_cluster_mob"))

  # birth_year is stored as a double and is not integer-valued for every study,
  # so it is rounded to a whole year before tabulating - the same round() the
  # frozen fig02_birthyear_counts.csv export applies, so the two files bin
  # identically.
  cnt <- function(d, nm) {
    d <- as.data.frame(d)
    n_all <- nrow(d)
    keep <- !is.na(d$cohort) & !is.na(d$birth_year)
    if (any(!keep)) {
      cat(sprintf(
        "counts: %s - dropped %d of %d rows with missing study or birth year\n",
        nm, sum(!keep), n_all
      ))
    }
    d <- d[keep, , drop = FALSE]
    by_bin <- as.integer(round(d$birth_year))
    tab <- as.data.frame(table(study = as.character(d$cohort), birth_year = by_bin),
      stringsAsFactors = FALSE
    )
    tab <- tab[tab$Freq > 0, , drop = FALSE]
    tab$birth_year <- as.integer(as.character(tab$birth_year))
    names(tab)[names(tab) == "Freq"] <- nm
    stopifnot(
      sum(tab[[nm]]) == nrow(d),
      !anyDuplicated(paste(tab$study, tab$birth_year))
    )
    tab
  }
  a <- cnt(dat_cluster, "n_attainment")
  b <- cnt(dat_cluster_mob, "n_mobility")

  out <- merge(a, b, by = c("study", "birth_year"), all = TRUE)
  out$n_attainment[is.na(out$n_attainment)] <- 0L
  out$n_mobility[is.na(out$n_mobility)] <- 0L
  out <- out[order(out$study, out$birth_year), , drop = FALSE]

  out$source_run <- RUN_TS_VAL
  out$source_commit <- source_commit()
  out$note <- paste(
    "Counts only. n_attainment is the educational attainment analytic sample",
    "(dat_cluster); n_mobility the educational mobility subsample",
    "(dat_cluster_mob; SHIP absent - it carries no parental education). The",
    "mobility subsample is nested within the attainment sample in every",
    "study x birth-year cell, so n_mobility <= n_attainment throughout. Birth",
    "year is rounded to a whole year, matching fig02_birthyear_counts.csv.",
    "Additional sensitivity outside the prespecified analysis plan."
  )

  # Nesting is asserted, not assumed: a future rebuild that broke it would stop
  # here rather than quietly produce a figure with impossible cells.
  over <- which(out$n_mobility > out$n_attainment)
  if (length(over)) {
    stop(sprintf(
      "%d of %d study x birth-year cells have n_mobility > n_attainment",
      length(over), nrow(out)
    ))
  }

  stopifnot(sum(out$n_mobility) <= sum(out$n_attainment))
  # Totals must reproduce the frozen run's Ns, minus only rows the counting step
  # itself reported as dropped for a missing study or birth year.
  cat(sprintf(
    "counts: totals - attainment %d (sample N = %d), mobility %d (sample N = %d)\n",
    sum(out$n_attainment), nrow(dat_cluster),
    sum(out$n_mobility), nrow(dat_cluster_mob)
  ))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "figure1_study_birthyear_counts.csv")
  utils::write.csv(out, path, row.names = FALSE, na = "")
  cat(sprintf(
    "counts: wrote %s (%d study x birth-year rows; N att = %d, N mob = %d)\n",
    path, nrow(out), sum(out$n_attainment), sum(out$n_mobility)
  ))
  invisible(out)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

RUN_TS_VAL <- Sys.getenv("RUN_TS", "20260729_1555")
# Every script in this folder writes to manuscript_export/sensitivity_analyses/
# of the run, so the manuscript repository imports from one folder. This
# script's default sits under the section's own output/, which is gitignored but
# inside the repository, so import_results.R resolves a real source commit and
# branch rather than NA.
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(
  SECTION_DIR, "output", paste0("RUN_", RUN_TS_VAL),
  "manuscript_export", "sensitivity_analyses"
))
step_counts(OUT_DIR)
