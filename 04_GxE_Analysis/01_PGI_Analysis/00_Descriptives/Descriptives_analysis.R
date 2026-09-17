# Descriptives_analysis.R
# Phenotype-descriptives module — analysis stage.
# Loads data, builds the full-family (dat_cluster) and dedup one-adult-per-family
# (dat_edu) analytic samples, and computes the descriptive tables (overall +
# three birth-year bins, full sample + per study). All outputs are
# AGGREGATE ONLY (N, M(SD), Welch p, Hedges' g) — the saved RDS holds the
# computed tables, never individual-level rows.
#
# Usage (from run_Descriptives.command):
#   Rscript 01_PGI_Analysis/00_Descriptives/Descriptives_analysis.R
#
# Outputs (under runs/RUN_<TS>/00_Descriptives/):
#   Descriptives_Tables.rds  — coverage + main + per-study table list
#
# Descriptives_report.R then reads this RDS to build
# Supplementary_Data_1_Descriptives_results.xlsx and Descriptives_Overview.md.
#
# Sample / test choices (see ../../../Plan_deviations.md §5):
#   * N and M(SD) are reported on the full-family RE sample (dat_cluster /
#     dat_cluster_mob), matching the primary A/B analyses; the two-group tests
#     use the dedup one-adult-per-family sample (dat_edu) for clean
#     independence, so the test Ns are smaller.
#   * Region difference = East - West; Gender difference = Female - Male.
#   * Welch t-test p + Hedges' g; computed only when both groups have
#     >= DESC_MIN_CELL_N finite observations.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
  library(openxlsx)
})

# ---- Resolve paths (works under both Rscript and source()) ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine Descriptives_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR  # flattened: 01_PGI_Analysis is the workflow + script root
rm(.find_this_file)

# ---- Run context: RUN_TS and output dirs ----
source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL <- get_run_ts()
OUT_DIR    <- output_dir("00_Descriptives")

# ---- Constants and helpers ----
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "helpers_descriptive_tables.R"),   local = TRUE)

# ---- Result collectors (required by 00_setup/ scripts) ----
RESULTS         <- list()
PLOTS           <- list()
PRIMARY_RESULTS <- list()

# ---- Load + prepare samples ----
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var
DATA_DIR         <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")

source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

# ---- Samples ----
# N and M(SD) are reported on the full-family RE sample (dat_cluster /
# dat_cluster_mob, which the A/B analyses use); the East-West / Female-Male
# tests are computed on the dedup one-adult-per-family sample (dat_edu) for
# clean independence. See ../../../Plan_deviations.md §5.
full_sample  <- if (exists("dat_cluster")) dat_cluster else dat_edu
dedup_sample <- dat_edu
cat(sprintf("Descriptives_analysis: full-family N = %d; dedup (test) N = %d\n",
            nrow(full_sample), nrow(dedup_sample)))

# Ceiling / floor indicators (as percentages): the share of each birth-cohort
# group at the top code (18y = university) and the lower end (<= 9y, i.e. at
# most lower-secondary). Together with the within-cohort education SD (the
# "Years of education" M(SD) row), these separate genuine variance compression
# from mechanical saturation of the years scale: edu_z_kernel divides by the
# within-birth-year education SD (Plan_deviations.md §2/§3), so a stable SD and
# a stable ceiling share mean the cohort trend in the PGI association is not a
# denominator artifact.
full_sample$edu_ceiling  <- 100 * as.integer(full_sample$education  >= 18)
full_sample$edu_floor    <- 100 * as.integer(full_sample$education  <= 9)
dedup_sample$edu_ceiling <- 100 * as.integer(dedup_sample$education >= 18)
dedup_sample$edu_floor   <- 100 * as.integer(dedup_sample$education <= 9)

# ---- Measures to describe ----
# Years of education (raw, interpretable) plus the analysis outcome
# (edu_z_kernel), the ceiling/floor shares (scale-saturation check), mobility
# (kernel-z difference), PGI-Education (z) and parental education (raw years).
# Mobility / parental-education cells restrict automatically to non-missing
# observations (SHIP has neither).
measures <- list(
  "Years of education"        = "education",
  "Education (edu_z_kernel)"  = "edu_z_kernel",
  "% at ceiling (edu = 18y)"  = "edu_ceiling",
  "% at floor (edu <= 9y)"    = "edu_floor",
  "Mobility (kernel-z diff)"  = "mobility",
  "PGI-Education (z)"         = "PGI_Edu_z",
  "Parental education (yrs)"  = "parental_education"
)

# Studies whose birth-year span is too narrow for a meaningful by-year
# split (BASE-II is ~1927-1951) skip the per-study Table 2.
BYYEAR_SKIP <- c("BASE-II")

# ---- Coverage (per study) ----
coverage <- build_coverage_table(full_sample, dedup_sample)

# ---- Full-sample tables ----
main_overall <- build_desc_table(full_sample, dedup_sample, measures)
main_by_year <- build_desc_table_by_group(full_sample, dedup_sample, measures)

# ---- Per-study tables ----
cohorts <- sort(as.character(unique(full_sample$cohort)))
per_cohort <- lapply(cohorts, function(ch) {
  fs <- full_sample[full_sample$cohort == ch, , drop = FALSE]
  ds <- dedup_sample[dedup_sample$cohort == ch, , drop = FALSE]
  list(
    overall     = build_desc_table(fs, ds, measures),
    by_year     = build_desc_table_by_group(fs, ds, measures),
    show_byyear = !(ch %in% BYYEAR_SKIP),
    N           = nrow(fs)
  )
})
names(per_cohort) <- cohorts

# ---- Analysis-sample descriptives ----
# Attrition, PGI intercorrelations and leave-one-study-out sample composition,
# all on the full-family RE sample the A/B analyses use.
attrition <- if (!is.null(PRIMARY_RESULTS[["Attrition_Core"]]))
  PRIMARY_RESULTS[["Attrition_Core"]] else RESULTS[["Attrition_Core"]]

pgi_cols <- intersect(c("PGI_Edu", "PGI_Cog", "PGI_nonCog", "PGI_Height"), names(full_sample))
.pm  <- apply(as.matrix(full_sample[, pgi_cols, drop = FALSE]), 2,
              function(z) suppressWarnings(as.numeric(z)))
.cor <- cor(.pm, use = "pairwise.complete.obs")
pgi_correlations <- data.frame(PGI = rownames(.cor), round(as.data.frame(.cor), 3),
                               check.names = FALSE, stringsAsFactors = FALSE)

loo_folds <- c(list("Full sample" = character(0)), cohort_groups(full_sample))
loo_composition <- dplyr::bind_rows(lapply(names(loo_folds), function(fl) {
  dl  <- loo_folds[[fl]]
  sub <- if (length(dl) == 0) full_sample else
    full_sample[!(as.character(full_sample$cohort) %in% dl), , drop = FALSE]
  data.frame(
    fold       = if (length(dl) == 0) "Full sample" else paste0("-", fl),
    N          = nrow(sub),
    BY_mean    = round(mean(sub$birth_year, na.rm = TRUE), 1),
    Edu_mean   = round(mean(sub$education,  na.rm = TRUE), 2),
    PGI_mean   = round(mean(suppressWarnings(as.numeric(sub$PGI_Edu_z)), na.rm = TRUE), 3),
    pct_East   = round(100 * mean(sub$east_west == "East",   na.rm = TRUE), 1),
    pct_female = round(100 * mean(sub$gender    == "female", na.rm = TRUE), 1),
    stringsAsFactors = FALSE)
}))

# ---- Bundle (aggregate tables only) ----
tables <- list(
  RUN_TS       = RUN_TS_VAL,
  sample_label = "Full-family RE sample (dat_cluster); East-West/Female-Male tests on dedup one-per-family subset (dat_edu)",
  N_total      = nrow(full_sample),
  N_test       = nrow(dedup_sample),
  min_cell_n   = DESC_MIN_CELL_N,
  measures     = names(measures),
  coverage     = coverage,
  main_overall = main_overall,
  main_by_year = main_by_year,
  per_cohort   = per_cohort,
  attrition        = attrition,
  pgi_correlations = pgi_correlations,
  loo_composition  = loo_composition
)

RDS_PATH <- file.path(OUT_DIR, "Descriptives_Tables.rds")
saveRDS(tables, RDS_PATH)
cat(sprintf("Descriptives_analysis: wrote %s (%d study sheets)\n",
            RDS_PATH, length(per_cohort)))
