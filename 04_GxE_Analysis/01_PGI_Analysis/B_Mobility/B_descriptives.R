# B_descriptives.R
# Section B descriptives: intergenerational mobility over strata.
# Sourced by B_analysis.R after 00_setup/ has produced dat_cluster_mob.
#
# Produces in the calling frame:
#   desc_cohort         data.frame by cohort
#   desc_cohort_region  data.frame by cohort x east_west
#   pgi_correlations    PGI intercorrelation matrix (mobility sample)

stopifnot(exists("dat_cluster_mob"), exists("safe_numeric"))

# ---- By cohort ----
desc_cohort <- dat_cluster_mob |>
  dplyr::group_by(cohort) |>
  dplyr::summarise(
    N          = dplyr::n(),
    BY_min     = min(birth_year, na.rm = TRUE),
    BY_max     = max(birth_year, na.rm = TRUE),
    BY_mean    = round(mean(birth_year, na.rm = TRUE), 1),
    Mob_mean   = round(mean(mobility, na.rm = TRUE), 3),
    Mob_sd     = round(sd(mobility,   na.rm = TRUE), 3),
    ParEdu_mean= round(mean(parental_edu_z_kernel, na.rm = TRUE), 3),
    PGI_mean   = round(mean(PGI_Edu_z, na.rm = TRUE), 3),
    pct_East   = round(100 * mean(east_west == "East", na.rm = TRUE), 1),
    .groups    = "drop"
  )

# ---- By cohort x region ----
desc_cohort_region <- dat_cluster_mob |>
  dplyr::group_by(cohort, east_west) |>
  dplyr::summarise(
    N           = dplyr::n(),
    Mob_mean    = round(mean(mobility, na.rm = TRUE), 3),
    Mob_sd      = round(sd(mobility,   na.rm = TRUE), 3),
    ParEdu_mean = round(mean(parental_edu_z_kernel, na.rm = TRUE), 3),
    PGI_mean    = round(mean(PGI_Edu_z, na.rm = TRUE), 3),
    .groups     = "drop"
  )

# ---- PGI intercorrelation matrix (on the mobility sample) ----
pgi_cols <- intersect(c("PGI_Edu", "PGI_Cog", "PGI_nonCog", "PGI_Height"),
                      names(dat_cluster_mob))
pgi_mat <- apply(as.matrix(dat_cluster_mob[, pgi_cols, drop = FALSE]), 2, safe_numeric)
cor_mat <- cor(pgi_mat, use = "pairwise.complete.obs")
pgi_correlations <- data.frame(
  pgi = rownames(cor_mat), cor_mat, check.names = FALSE, stringsAsFactors = FALSE
)

cat(sprintf("B_descriptives: %d cohorts, dat_cluster_mob N = %d\n",
            length(unique(droplevels(dat_cluster_mob$cohort))), nrow(dat_cluster_mob)))
