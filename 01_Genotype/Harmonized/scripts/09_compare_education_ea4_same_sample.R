#!/usr/bin/env Rscript
# =============================================================================
# 09_compare_education_ea4_same_sample.R — same-sample cohort means / t-tests
# =============================================================================
# Cohort mean/SD and mean-difference tests for education and EA4 on the exact
# same complete-case sample. Local-machine cross-check against the analysis
# project; not part of the production DAG.
#
# USAGE:
#   Rscript scripts/09_compare_education_ea4_same_sample.R
#
# INPUTS:
#   <analysis project>/outputs/GxE_Germany_Combined_<date>.rds
#     (most recent — script picks the newest matching file)
#
# OUTPUTS:
#   Console summary tables only.
#
# DEPENDS ON: nothing in this pipeline directly. Consumes the analysis-side
# merged RDS that itself ingests `output/final/<COHORT>_PGI_PCs.tsv` from
# step 07.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Resolve ANALYSIS_ROOT from env or ../user_config.sh.
source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "_user_config.R"))
analysis_dir <- get_user_config_value("ANALYSIS_ROOT")
outputs_dir <- file.path(analysis_dir, "outputs")
fig_dir <- file.path(outputs_dir, "figures")

rds_files <- list.files(
  outputs_dir,
  pattern = "^GxE_Germany_Combined_.*\\.rds$",
  full.names = TRUE
)
if (length(rds_files) == 0) {
  stop("No GxE_Germany_Combined_*.rds file found in: ", outputs_dir)
}

input_rds <- rds_files[which.max(file.info(rds_files)$mtime)]
d <- as.data.table(readRDS(input_rds))

required_cols <- c("cohort", "birth_year", "education", "PGI_Edu")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

age_reference_year <- as.integer(format(Sys.Date(), "%Y"))
age_min <- 25
pgi_outlier_cutoff <- 3
d[, age_years := age_reference_year - birth_year]
d[, cohort := factor(cohort, levels = c("BASE-II", "SHIP", "SOEP", "TwinLife"))]

same_sample_before_outlier_filter <- d[
  is.finite(age_years) &
    age_years > age_min &
    is.finite(education) &
    is.finite(PGI_Edu)
]
same_sample <- same_sample_before_outlier_filter[abs(PGI_Edu) <= pgi_outlier_cutoff]

vars <- data.table(
  variable = c("education", "PGI_Edu"),
  label = c("Education", "EA4 PGI global_z")
)

sample_counts <- d[, .(
  n_total = .N,
  n_age_gt25 = sum(is.finite(age_years) & age_years > age_min),
  n_education = sum(is.finite(education)),
  n_pgi_edu = sum(is.finite(PGI_Edu)),
  n_same_sample_before_outlier_filter = sum(
    is.finite(age_years) &
      age_years > age_min &
      is.finite(education) &
      is.finite(PGI_Edu)
  ),
  pgi_edu_outlier_cutoff = pgi_outlier_cutoff,
  n_pgi_edu_outliers_same_sample = sum(
    is.finite(age_years) &
      age_years > age_min &
      is.finite(education) &
      is.finite(PGI_Edu) &
      abs(PGI_Edu) > pgi_outlier_cutoff
  ),
  n_same_sample = sum(
    is.finite(age_years) &
      age_years > age_min &
      is.finite(education) &
      is.finite(PGI_Edu) &
      abs(PGI_Edu) <= pgi_outlier_cutoff
  )
), by = cohort]

summary_rows <- rbindlist(lapply(seq_len(nrow(vars)), function(i) {
  v <- vars$variable[i]
  lab <- vars$label[i]
  same_sample[, .(
    variable = lab,
    n = .N,
    mean = mean(get(v)),
    sd = sd(get(v)),
    se = sd(get(v)) / sqrt(.N),
    min = min(get(v)),
    p25 = quantile(get(v), 0.25),
    median = median(get(v)),
    p75 = quantile(get(v), 0.75),
    max = max(get(v))
  ), by = cohort]
}))
setcolorder(summary_rows, c("variable", "cohort", "n", "mean", "sd", "se", "min", "p25", "median", "p75", "max"))

pairwise_rows <- rbindlist(lapply(seq_len(nrow(vars)), function(i) {
  v <- vars$variable[i]
  lab <- vars$label[i]
  dd <- same_sample[, .(cohort, value = get(v))]
  cohorts <- levels(droplevels(dd$cohort))
  out <- list()
  k <- 1L
  for (a in seq_along(cohorts)) {
    for (b in seq_along(cohorts)) {
      if (b <= a) next
      ca <- cohorts[a]
      cb <- cohorts[b]
      xa <- dd[cohort == ca, value]
      xb <- dd[cohort == cb, value]
      tt <- t.test(xa, xb)
      out[[k]] <- data.table(
        variable = lab,
        cohort_1 = ca,
        cohort_2 = cb,
        n_1 = length(xa),
        n_2 = length(xb),
        mean_1 = mean(xa),
        mean_2 = mean(xb),
        diff_1_minus_2 = mean(xa) - mean(xb),
        t = unname(tt$statistic),
        df = unname(tt$parameter),
        p = tt$p.value
      )
      k <- k + 1L
    }
  }
  res <- rbindlist(out)
  res[, p_holm := p.adjust(p, method = "holm")]
  res[, significant_holm_0.05 := p_holm < 0.05]
  res[]
}))

overall_rows <- rbindlist(lapply(seq_len(nrow(vars)), function(i) {
  v <- vars$variable[i]
  lab <- vars$label[i]
  dd <- same_sample[, .(cohort, value = get(v))]
  ow <- oneway.test(value ~ cohort, data = dd, var.equal = FALSE)
  data.table(
    variable = lab,
    test = "Welch one-way ANOVA",
    statistic = unname(ow$statistic),
    df_num = unname(ow$parameter[1]),
    df_den = unname(ow$parameter[2]),
    p = ow$p.value,
    n = nrow(dd)
  )
}))

date_stamp <- format(Sys.Date(), "%Y%m%d")
suffix <- paste0("same_sample_age_gt25_no_pgi_outliers_abs", pgi_outlier_cutoff, "_", date_stamp)
sample_file <- file.path(fig_dir, paste0("education_ea4_", suffix, "_counts.tsv"))
summary_file <- file.path(fig_dir, paste0("education_ea4_", suffix, "_means_sds.tsv"))
pairwise_file <- file.path(fig_dir, paste0("education_ea4_", suffix, "_pairwise_mean_tests.tsv"))
overall_file <- file.path(fig_dir, paste0("education_ea4_", suffix, "_overall_mean_tests.tsv"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(sample_counts, sample_file, sep = "\t")
fwrite(summary_rows, summary_file, sep = "\t")
fwrite(pairwise_rows, pairwise_file, sep = "\t")
fwrite(overall_rows, overall_file, sep = "\t")

cat("Input dataset:\n  ", input_rds, "\n", sep = "")
cat(sprintf("Same-sample filter:\n  age = %d - birth_year > %d; education observed; PGI_Edu observed; abs(PGI_Edu) <= %s\n",
            age_reference_year, age_min, pgi_outlier_cutoff))
cat("Saved:\n")
cat("  ", sample_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
cat("  ", pairwise_file, "\n", sep = "")
cat("  ", overall_file, "\n", sep = "")
cat("\nSample counts:\n")
print(sample_counts)
cat("\nMeans/SDs on same sample:\n")
print(summary_rows)
cat("\nOverall Welch tests:\n")
print(overall_rows)
cat("\nPairwise Welch tests with Holm correction:\n")
print(pairwise_rows[, .(variable, cohort_1, cohort_2, diff_1_minus_2, p, p_holm, significant_holm_0.05)])
