#!/usr/bin/env Rscript
# =============================================================================
# 08_plot_pgi_qc.R — PGI QC plots from the analysis-side merged dataset
# =============================================================================
# Plots PGI distributions and education-vs-EA4 scatters from the merged
# harmonized analysis dataset. Local-machine diagnostic; not part of the
# production DAG.
#
# USAGE:
#   Rscript scripts/08_plot_pgi_qc.R
#   (runs on the Mac with the downstream analysis project mounted)
#
# INPUTS:
#   <analysis project>/outputs/GxE_Germany_Combined_<date>.rds
#     — the most recent merged phenotype + PGI dataset built by the analysis
#       project's `GxE_Germany_Data_Merge_harmonized.R`.
#
# OUTPUTS:
#   <analysis project>/outputs/figures/*.pdf  — distribution & scatter plots
#
# DEPENDS ON: nothing in this pipeline directly. Consumes the analysis-side
# merged RDS that itself ingests `output/final/<COHORT>_PGI_PCs.tsv` from
# step 07.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

d <- as.data.table(readRDS(input_rds))

required_cols <- c("cohort", "education", "PGI_Edu", "PGI_Cog", "PGI_nonCog", "PGI_Height")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) {
  stop("Missing required columns in merged dataset: ", paste(missing_cols, collapse = ", "))
}

age_reference_year <- as.integer(format(Sys.Date(), "%Y"))
age_min <- 25
pgi_outlier_cutoff <- 3
if (!"birth_year" %in% names(d)) {
  stop("Missing required column for age filtering: birth_year")
}
d[, age_years := age_reference_year - birth_year]
d_all <- copy(d)
d_age <- d[is.finite(age_years) & age_years > age_min]
d <- d_age[is.finite(PGI_Edu) & abs(PGI_Edu) <= pgi_outlier_cutoff]

plot_stamp <- format(Sys.Date(), "%Y%m%d")
age_suffix <- paste0("age_gt", age_min, "_no_pgi_outliers_abs", pgi_outlier_cutoff)

pgi_cols <- c("PGI_Edu", "PGI_Cog", "PGI_Height")
pgi_labels <- data.table(
  variable = pgi_cols,
  PGI = c("EA4 / Education", "Cognition", "Height")
)

pgi_long <- melt(
  d,
  id.vars = "cohort",
  measure.vars = pgi_cols,
  variable.name = "variable",
  value.name = "score"
)
pgi_long <- merge(pgi_long, pgi_labels, by = "variable", all.x = TRUE)
pgi_long <- pgi_long[is.finite(score)]
pgi_long[, cohort := factor(cohort, levels = c("BASE-II", "SHIP", "SOEP", "TwinLife"))]
pgi_long[, PGI := factor(PGI, levels = pgi_labels$PGI)]

hist_base <- ggplot(pgi_long, aes(x = score)) +
  geom_vline(xintercept = 0, color = "#B22222", linewidth = 0.5) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 45,
    fill = "#3E6B89",
    color = "white",
    linewidth = 0.15
  ) +
  geom_density(color = "#C75B39", linewidth = 0.45, na.rm = TRUE) +
  facet_grid(cohort ~ PGI, scales = "free_y") +
  labs(
    title = "Harmonized PGI Distributions by Cohort, Age > 25",
    subtitle = sprintf(
      "%s | age = %d - birth_year > %d | abs(PGI_Edu) <= %s",
      basename(input_rds), age_reference_year, age_min, pgi_outlier_cutoff
    ),
    x = "Global harmonized PGI z-score",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8)
  )

hist_full_png <- file.path(fig_dir, paste0("pgi_histograms_by_cohort_full_range_", age_suffix, "_", plot_stamp, ".png"))
hist_full_pdf <- file.path(fig_dir, paste0("pgi_histograms_by_cohort_full_range_", age_suffix, "_", plot_stamp, ".pdf"))
ggsave(hist_full_png, hist_base, width = 14, height = 9, dpi = 300)
ggsave(hist_full_pdf, hist_base, width = 14, height = 9)

hist_zoom <- hist_base +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(
    title = "Harmonized PGI Distributions by Cohort, Age > 25, Zoomed to +/-3 SD",
    x = "Global harmonized PGI z-score (zoomed)"
  )

hist_zoom_png <- file.path(fig_dir, paste0("pgi_histograms_by_cohort_zoom_pm3_", age_suffix, "_", plot_stamp, ".png"))
hist_zoom_pdf <- file.path(fig_dir, paste0("pgi_histograms_by_cohort_zoom_pm3_", age_suffix, "_", plot_stamp, ".pdf"))
ggsave(hist_zoom_png, hist_zoom, width = 14, height = 9, dpi = 300)
ggsave(hist_zoom_pdf, hist_zoom, width = 14, height = 9)

scatter_dt <- d[is.finite(education) & is.finite(PGI_Edu)]
scatter_dt[, cohort := factor(cohort, levels = c("BASE-II", "SHIP", "SOEP", "TwinLife"))]
scatter_dt[, PGI_Edu_within_cohort_z := {
  s <- sd(PGI_Edu)
  if (is.na(s) || s == 0) rep(NA_real_, .N) else (PGI_Edu - mean(PGI_Edu)) / s
}, by = cohort]

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

global_cor_dt <- scatter_dt[, {
  fit <- lm(education ~ PGI_Edu)
  coef_tab <- coef(summary(fit))
  r_val <- cor(PGI_Edu, education)
  .(
    n = .N,
    r = r_val,
    r2 = r_val^2,
    slope = unname(coef_tab["PGI_Edu", "Estimate"]),
    se = unname(coef_tab["PGI_Edu", "Std. Error"]),
    p = unname(coef_tab["PGI_Edu", "Pr(>|t|)"])
  )
}, by = cohort]
global_cor_dt[, label := sprintf("r = %.3f\nR2 = %.3f\nn = %s", r, r2, format(n, big.mark = ","))]

scatter_global_plot <- ggplot(scatter_dt, aes(x = PGI_Edu, y = education)) +
  geom_point(color = "#2B4C5C", alpha = 0.16, size = 0.55) +
  geom_smooth(method = "lm", se = FALSE, color = "#C75B39", linewidth = 0.8) +
  geom_text(
    data = global_cor_dt,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    size = 3.5
  ) +
  facet_wrap(~ cohort, ncol = 2, scales = "free_y") +
  labs(
    title = "Education by EA4 PGI",
    subtitle = sprintf(
      "%s | global harmonized EA4 PGI | age = %d - birth_year > %d | abs(PGI_Edu) <= %s",
      basename(input_rds), age_reference_year, age_min, pgi_outlier_cutoff
    ),
    x = "PGI_Edu (EA4, global harmonized z-score)",
    y = "Education"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

scatter_global_png <- file.path(fig_dir, paste0("education_vs_pgi_ea4_by_cohort_", age_suffix, "_", plot_stamp, ".png"))
scatter_global_pdf <- file.path(fig_dir, paste0("education_vs_pgi_ea4_by_cohort_", age_suffix, "_", plot_stamp, ".pdf"))
ggsave(scatter_global_png, scatter_global_plot, width = 10.5, height = 8, dpi = 300)
ggsave(scatter_global_pdf, scatter_global_plot, width = 10.5, height = 8)

cor_dt <- scatter_dt[, {
  fit <- lm(education ~ PGI_Edu_within_cohort_z)
  coef_tab <- coef(summary(fit))
  r_val <- cor(PGI_Edu_within_cohort_z, education)
  .(
  n = .N,
    r = r_val,
    r2 = r_val^2,
    slope = unname(coef_tab["PGI_Edu_within_cohort_z", "Estimate"]),
    se = unname(coef_tab["PGI_Edu_within_cohort_z", "Std. Error"]),
    p = unname(coef_tab["PGI_Edu_within_cohort_z", "Pr(>|t|)"])
  )
}, by = cohort]
cor_dt[, label := sprintf(
  "b = %.3f (SE %.3f)\np %s\nR2 = %.3f\nn = %s",
  slope, se, vapply(p, format_p, character(1)), r2, format(n, big.mark = ",")
)]

scatter_plot <- ggplot(scatter_dt, aes(x = PGI_Edu_within_cohort_z, y = education)) +
  geom_point(color = "#2B4C5C", alpha = 0.16, size = 0.55) +
  geom_smooth(method = "lm", se = FALSE, color = "#C75B39", linewidth = 0.8) +
  geom_text(
    data = cor_dt,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    size = 3.5
  ) +
  facet_wrap(~ cohort, ncol = 2, scales = "free_y") +
  labs(
    title = "Education by Within-Cohort Standardized EA4 PGI",
    subtitle = sprintf(
      "%s | global EA4 PGI standardized within cohort after filters | age = %d - birth_year > %d | abs(PGI_Edu) <= %s",
      basename(input_rds), age_reference_year, age_min, pgi_outlier_cutoff
    ),
    x = "EA4 PGI: global harmonized score, z-standardized within cohort",
    y = "Education"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

scatter_png <- file.path(fig_dir, paste0("education_vs_pgi_ea4_within_cohort_z_by_cohort_", age_suffix, "_", plot_stamp, ".png"))
scatter_pdf <- file.path(fig_dir, paste0("education_vs_pgi_ea4_within_cohort_z_by_cohort_", age_suffix, "_", plot_stamp, ".pdf"))
ggsave(scatter_png, scatter_plot, width = 10.5, height = 8, dpi = 300)
ggsave(scatter_pdf, scatter_plot, width = 10.5, height = 8)

age_counts <- data.table(
  cohort = sort(unique(d_all$cohort))
)[, {
  dc <- d_all[cohort == .BY$cohort]
  dc_age <- dc[is.finite(age_years) & age_years > age_min]
  .(
    n_total = nrow(dc),
    n_age_observed = sum(is.finite(dc$age_years)),
    n_age_gt25 = sum(is.finite(dc$age_years) & dc$age_years > age_min),
    n_age_le25 = sum(is.finite(dc$age_years) & dc$age_years <= age_min),
    n_age_missing = sum(!is.finite(dc$age_years)),
    n_pgi_edu = sum(is.finite(dc$PGI_Edu)),
    n_pgi_edu_age_gt25 = sum(is.finite(dc$PGI_Edu) & is.finite(dc$age_years) & dc$age_years > age_min),
    pgi_edu_outlier_cutoff = pgi_outlier_cutoff,
    n_pgi_edu_outliers_age_gt25 = sum(is.finite(dc_age$PGI_Edu) & abs(dc_age$PGI_Edu) > pgi_outlier_cutoff),
    n_pgi_edu_after_filters = sum(is.finite(dc_age$PGI_Edu) & abs(dc_age$PGI_Edu) <= pgi_outlier_cutoff),
    n_edu_pgi = sum(is.finite(dc$education) & is.finite(dc$PGI_Edu)),
    n_edu_pgi_age_gt25 = sum(is.finite(dc$education) & is.finite(dc$PGI_Edu) & is.finite(dc$age_years) & dc$age_years > age_min),
    n_edu_pgi_after_filters = sum(is.finite(dc_age$education) & is.finite(dc_age$PGI_Edu) & abs(dc_age$PGI_Edu) <= pgi_outlier_cutoff)
  )
}, by = cohort]

pgi_summary <- pgi_long[, .(
  n = .N,
  mean = mean(score),
  sd = sd(score),
  min = min(score),
  p01 = quantile(score, 0.01),
  p05 = quantile(score, 0.05),
  median = median(score),
  p95 = quantile(score, 0.95),
  p99 = quantile(score, 0.99),
  max = max(score)
), by = .(cohort, PGI)]

summary_file <- file.path(fig_dir, paste0("pgi_distribution_summary_", age_suffix, "_", plot_stamp, ".tsv"))
global_cor_file <- file.path(fig_dir, paste0("education_vs_pgi_ea4_global_correlations_", age_suffix, "_", plot_stamp, ".tsv"))
cor_file <- file.path(fig_dir, paste0("education_vs_pgi_ea4_within_cohort_z_correlations_", age_suffix, "_", plot_stamp, ".tsv"))
age_counts_file <- file.path(fig_dir, paste0("age_filter_counts_", age_suffix, "_", plot_stamp, ".tsv"))
fwrite(pgi_summary, summary_file, sep = "\t")
fwrite(global_cor_dt, global_cor_file, sep = "\t")
fwrite(cor_dt, cor_file, sep = "\t")
fwrite(age_counts, age_counts_file, sep = "\t")

cat("Input dataset:\n  ", input_rds, "\n", sep = "")
cat(sprintf(
  "Filters:\n  age = %d - birth_year > %d\n  abs(PGI_Edu) <= %s\n",
  age_reference_year, age_min, pgi_outlier_cutoff
))
cat("Saved plots:\n")
cat("  ", hist_full_png, "\n", sep = "")
cat("  ", hist_zoom_png, "\n", sep = "")
cat("  ", scatter_global_png, "\n", sep = "")
cat("  ", scatter_png, "\n", sep = "")
cat("Saved tables:\n")
cat("  ", summary_file, "\n", sep = "")
cat("  ", global_cor_file, "\n", sep = "")
cat("  ", cor_file, "\n", sep = "")
cat("  ", age_counts_file, "\n", sep = "")
cat("\nAge filter counts:\n")
print(age_counts)
cat("\nEA4 global correlations:\n")
print(global_cor_dt)
cat("\nEA4 within-cohort-z correlations and slopes:\n")
print(cor_dt)
