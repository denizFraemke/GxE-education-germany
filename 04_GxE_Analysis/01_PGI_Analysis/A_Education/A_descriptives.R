# A_descriptives.R
# Section A descriptives: educational attainment over strata, no gender split.
# Sourced by A_analysis.R after 00_setup/ has produced dat_cluster.
#
# Produces in the calling frame:
#   desc_cohort      data.frame by cohort
#   desc_cohort_region   data.frame by cohort x east_west
#   pgi_correlations data.frame: PGI intercorrelation matrix

stopifnot(exists("dat_cluster"),
          exists("safe_numeric"))

# ---- 4a. By cohort ----
desc_cohort <- dat_cluster |>
  dplyr::group_by(cohort) |>
  dplyr::summarise(
    N          = dplyr::n(),
    BY_min     = min(birth_year, na.rm = TRUE),
    BY_max     = max(birth_year, na.rm = TRUE),
    BY_mean    = round(mean(birth_year, na.rm = TRUE), 1),
    BY_sd      = round(sd(birth_year, na.rm = TRUE), 1),
    Edu_mean   = round(mean(education, na.rm = TRUE), 2),
    Edu_sd     = round(sd(education,   na.rm = TRUE), 2),
    PGI_mean   = round(mean(PGI_Edu_z, na.rm = TRUE), 3),
    PGI_sd     = round(sd(PGI_Edu_z,   na.rm = TRUE), 3),
    pct_East   = round(100 * mean(east_west == "East", na.rm = TRUE), 1),
    .groups    = "drop"
  )

# ---- 4b. By cohort x region (no gender split for Section A) ----
desc_cohort_region <- dat_cluster |>
  dplyr::group_by(cohort, east_west) |>
  dplyr::summarise(
    N        = dplyr::n(),
    BY_mean  = round(mean(birth_year, na.rm = TRUE), 1),
    BY_sd    = round(sd(birth_year,   na.rm = TRUE), 1),
    Edu_mean = round(mean(education,  na.rm = TRUE), 2),
    Edu_sd   = round(sd(education,    na.rm = TRUE), 2),
    PGI_mean = round(mean(PGI_Edu_z,  na.rm = TRUE), 3),
    PGI_sd   = round(sd(PGI_Edu_z,    na.rm = TRUE), 3),
    .groups  = "drop"
  )

# ---- 4c. PGI intercorrelation matrix ----
pgi_cols <- intersect(c("PGI_Edu", "PGI_Cog", "PGI_nonCog", "PGI_Height"),
                      names(dat_cluster))
pgi_mat <- as.matrix(dat_cluster[, pgi_cols, drop = FALSE])
pgi_mat <- apply(pgi_mat, 2, safe_numeric)
cor_mat <- cor(pgi_mat, use = "pairwise.complete.obs")
pgi_correlations <- data.frame(
  pgi = rownames(cor_mat),
  cor_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(sprintf("A_descriptives: %d cohorts, dat_cluster N = %d\n",
            length(unique(dat_cluster$cohort)), nrow(dat_cluster)))
