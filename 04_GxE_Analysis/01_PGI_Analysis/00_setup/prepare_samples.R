# prepare_samples.R
# Build the analytic samples from df_raw.
#
# After source()ing, the calling environment contains:
#   df              : exclusions applied, all derived columns
#   dat_edu, dat_mob: supplementary (dedup-OLS) analytic samples, used by the
#                     attrition / mobility / descriptive plots and by the
#                     sensitivities that stay on the dedup track
#   dat_cluster     : primary RE analytic sample (full-family TwinLife,
#                     no family dedup; ready for + s(fid_re, bs = "re"))
#   by_mean, EDU_SD : centring constant and global education SD
#
# Pre-conditions in the caller's frame:
#   df_raw     : loaded data frame from load_data.R
#   gender_col : "gender" or "Gender"
#   RESULTS    : list collector (for TwinLife_FamilyStructure etc.)
#   PLOTS      : list collector (for Kernel_OwnEdu and Kernel_ParEdu)
#   PRIMARY_RESULTS : list collector (for primary-track Attrition_Core)
# Plus the analytic constants from constants.R (MAX_ABS_PGI_EDU,
# KERNEL_BW, PARENT_GENERATION_GAP, COL_KERNEL_EAST/WEST, ...).

stopifnot(exists("df_raw"),
          exists("gender_col"),
          exists("RESULTS"),
          exists("PLOTS"),
          exists("PRIMARY_RESULTS"))

# =====================================================================
# (a) Variable preparation: build df, dat_edu, dat_mob; attrition; kernel-z.
# =====================================================================

df <- df_raw

# ---- 3a. Harmonize column names ----
if (gender_col == "Gender" && !"gender" %in% names(df)) {
  names(df)[names(df) == "Gender"] <- "gender"
}
if ("Height" %in% names(df) && !"height" %in% names(df)) {
  names(df)[names(df) == "Height"] <- "height"
}

# ---- 3b. Numeric types ----
df$birth_year  <- safe_numeric(df$birth_year)
df$education   <- safe_numeric(df$education)
df$PGI_Edu     <- safe_numeric(df$PGI_Edu)
df$PGI_Cog     <- safe_numeric(df$PGI_Cog)
df$PGI_nonCog  <- safe_numeric(df$PGI_nonCog)
df$PGI_Height  <- safe_numeric(df$PGI_Height)
if ("parental_education" %in% names(df))
  df$parental_education <- safe_numeric(df$parental_education)
if ("height" %in% names(df)) df$height <- safe_numeric(df$height)

# ---- 3b2. TwinLife family structure ----
# fid: family ID (only TwinLife has it; NA for other cohorts)
# ptyp: person type in TwinLife
#   1, 2 = twins
#   200, 201 = non-twin children/siblings
#   110, 120, 300+ = parents or older relatives
# Both tracks exclude TwinLife twins AND non-twin children/siblings
# (`ptyp` 1, 2, 200, 201) so TwinLife contributes only its parent/adult
# generation. This is a design-scope choice (the adults span the historical
# cohorts of interest; keeping the children would concentrate the youngest
# cohort almost entirely in TwinLife), NOT an age-completeness filter --
# completeness is enforced separately by the `birth_year <= 1996` cutoff
# below, which some excluded twins already pass. See Plan_deviations.md
# section 5. The supplementary track additionally keeps at most one
# adult per family (random draw within `fid`); the primary track keeps
# all family adults and absorbs their within-family correlation via a
# `s(fid_re, bs = "re")` random intercept. See Random_effects_family_clustering.md.
if ("ptyp" %in% names(df)) {
  df$ptyp <- safe_numeric(df$ptyp)
  df$excluded_twin_sibling <- ifelse(df$cohort == "TwinLife" & df$ptyp %in% c(1, 2, 200, 201), TRUE, FALSE)
} else {
  df$excluded_twin_sibling <- FALSE
}
# fid: keep for TwinLife, set to NA for other cohorts (for clustering)
if (!"fid" %in% names(df)) df$fid <- NA_character_
df$fid <- ifelse(df$cohort == "TwinLife", as.character(df$fid), NA_character_)
if (!"pid" %in% names(df)) df$pid <- seq_len(nrow(df))

df$excluded_twinlife_family_duplicate <- FALSE
tw_adult_idx <- which(df$cohort == "TwinLife" &
                        !df$excluded_twin_sibling &
                        !is.na(df$fid))
if (length(tw_adult_idx) > 0) {
  tw_keys <- data.frame(
    idx = tw_adult_idx,
    fid = df$fid[tw_adult_idx],
    pid = as.character(df$pid[tw_adult_idx]),
    stringsAsFactors = FALSE
  )
  set.seed(20260410)
  keep_idx <- tw_keys |>
    dplyr::group_by(fid) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::pull(idx)
  dup_idx <- setdiff(tw_adult_idx, keep_idx)
  if (length(dup_idx) > 0) df$excluded_twinlife_family_duplicate[dup_idx] <- TRUE
}

tw_all_adults <- df$cohort == "TwinLife" & !df$excluded_twin_sibling & !is.na(df$fid)
tw_adult_family_sizes <- table(df$fid[tw_all_adults])
RESULTS[["TwinLife_FamilyStructure"]] <- data.frame(
  metric = c(
    "TwinLife twins/non-twin children excluded from primary",
    "TwinLife adult records eligible before family dedup",
    "TwinLife unique adult families eligible",
    "TwinLife families with >1 adult retained candidate",
    "TwinLife duplicate adults excluded from primary"
  ),
  value = c(
    sum(df$excluded_twin_sibling, na.rm = TRUE),
    sum(tw_all_adults, na.rm = TRUE),
    length(tw_adult_family_sizes),
    sum(tw_adult_family_sizes > 1),
    sum(df$excluded_twinlife_family_duplicate, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat(sprintf("TwinLife twins+siblings (excluded from primary): N = %d\n",
            sum(df$excluded_twin_sibling, na.rm = TRUE)))
cat(sprintf("TwinLife duplicate adults excluded to keep one per family: N = %d\n",
            sum(df$excluded_twinlife_family_duplicate, na.rm = TRUE)))

# ---- 3c. east_west as factor ----
ew_lower <- tolower(trimws(as.character(df$east_west)))
df$east_west <- ifelse(ew_lower == "west", "West",
                       ifelse(ew_lower == "east", "East", NA_character_))
df$east_west <- factor(df$east_west, levels = c("West", "East"))

# ---- 3d. Gender as factor ----
gen_raw <- tolower(trimws(as.character(df$gender)))
# Handle numeric coding (0=female, 1=male) or character
df$gender <- dplyr::case_when(
  gen_raw %in% c("female", "0") ~ "female",
  gen_raw %in% c("male",   "1") ~ "male",
  TRUE ~ NA_character_
)
df$gender <- factor(df$gender, levels = c("female", "male"))

# ---- 3e. Cohort as factor ----
df$cohort <- factor(df$cohort)

# ---- 3f. Numeric helpers ----
# Contrast-coded predictors for model formulas (ANALYSIS_PLAN.md, "Note on coding"):
#   east_west_c = -0.5 (West) / +0.5 (East)
#   gender_c    = +0.5 (female) / -0.5 (male)
# 0/1 helpers (east, gender_01) are retained ONLY for cell-indicator construction
# of region- and gender-specific PGI columns (PGI_W, PGI_E, PGI_WF, ..., PGI_EM).
df$east_west_c <- ifelse(df$east_west == "East",   0.5, -0.5)
df$gender_c    <- ifelse(df$gender    == "female", 0.5, -0.5)
df$east        <- as.integer(df$east_west == "East")
df$gender_01   <- as.integer(df$gender    == "male")

# ---- 3g. Attrition tracking ----
# Per-factor attrition steps. The TwinLife child-exclusion (twins +
# non-twin siblings) is shared by both tracks. The family-dedup step
# (one-adult-per-family) is ONLY applied in the supplementary dedup
# track; the primary RE track keeps all family adults and applies a
# random intercept instead, so its attrition list omits that final step.
core_attrition_steps_shared <- list(
  list(label = "Birth year observed", keep = function(d) !is.na(d$birth_year)),
  list(label = "Birth year <= 1996",  keep = function(d) d$birth_year <= 1996),
  list(label = "Education observed",  keep = function(d) !is.na(d$education)),
  list(label = "Region observed",     keep = function(d) !is.na(d$east_west)),
  list(label = "Gender observed",     keep = function(d) !is.na(d$gender)),
  list(label = "PGI-Edu observed",    keep = function(d) !is.na(d$PGI_Edu)),
  list(label = sprintf("PGI-Edu within +/- %s", MAX_ABS_PGI_EDU),
       keep = function(d) abs(d$PGI_Edu) <= MAX_ABS_PGI_EDU),
  list(label = "TwinLife children excluded",
       keep = function(d) !dplyr::coalesce(d$excluded_twin_sibling, FALSE))
)
# Supplementary (dedup) track: shared steps + family deduplication.
core_attrition_steps <- c(
  core_attrition_steps_shared,
  list(list(label = "TwinLife one adult per family",
            keep = function(d) !dplyr::coalesce(d$excluded_twinlife_family_duplicate, FALSE)))
)
# Primary (RE) track: shared steps only - no family dedup step.
core_attrition_steps_primary <- core_attrition_steps_shared

attrition_core <- compute_attrition_flow(df, core_attrition_steps, pipeline = "Core sample")
RESULTS[["Attrition_Core"]] <- attrition_core
attrition_core_primary <- compute_attrition_flow(df, core_attrition_steps_primary,
                                                 pipeline = "Core sample")
PRIMARY_RESULTS[["Attrition_Core"]] <- attrition_core_primary

# ---- 3h. Apply exclusions ----
n_twin_sib_excl <- sum(df$cohort == "TwinLife" & df$ptyp %in% c(1, 2, 200, 201), na.rm = TRUE)
n_twinlife_dup_excl <- sum(df$excluded_twinlife_family_duplicate, na.rm = TRUE)
n_outlier_excl <- sum(!is.na(df$PGI_Edu) & abs(df$PGI_Edu) > MAX_ABS_PGI_EDU, na.rm = TRUE)
n_ship_west_excl <- sum(df$cohort == "SHIP" & df$east_west == "West", na.rm = TRUE)

df <- df |>
  dplyr::filter(
    !is.na(birth_year),
    birth_year <= 1996,
    !is.na(education),
    !is.na(east_west),
    !is.na(gender),
    !is.na(PGI_Edu),
    abs(PGI_Edu) <= MAX_ABS_PGI_EDU,
    # Exclude twins AND non-twin siblings from the main analysis
    !excluded_twin_sibling,
    # Keep one adult per TwinLife family in the primary sample
    !excluded_twinlife_family_duplicate,
    # SHIP was recruited in the East (Greifswald, Western Pomerania) — exclude
    # the few West cases
    !(cohort == "SHIP" & east_west == "West")
  )
cat(sprintf(
  "After exclusions: N = %d (excluded %d EA4 outliers, %d twins+siblings, %d TwinLife duplicate adults, %d SHIP West)\n",
  nrow(df), n_outlier_excl, n_twin_sib_excl, n_twinlife_dup_excl, n_ship_west_excl
))

tw_primary_family_sizes <- table(df$fid[df$cohort == "TwinLife" & !is.na(df$fid)])
RESULTS[["TwinLife_PrimaryFamilyCheck"]] <- data.frame(
  metric = c(
    "TwinLife adults in primary sample",
    "TwinLife unique families in primary sample",
    "TwinLife families with >1 adult in primary sample",
    "TwinLife max adults per family in primary sample"
  ),
  value = c(
    sum(df$cohort == "TwinLife", na.rm = TRUE),
    length(tw_primary_family_sizes),
    sum(tw_primary_family_sizes > 1),
    if (length(tw_primary_family_sizes) > 0) max(tw_primary_family_sizes) else 0
  ),
  stringsAsFactors = FALSE
)

# ---- 3i. Center and standardize ----
by_mean <- mean(df$birth_year, na.rm = TRUE)
by_sd   <- sd(df$birth_year,   na.rm = TRUE)

df$BYc   <- df$birth_year - by_mean          # centered (for GAMs)
df$BYc_z <- z_scale(df$birth_year)           # z-standardized (for linear models)

# PGI: loaded as the global harmonized PGI from the merge. Residualise each PGI on
# cross-cohort ancestry PCs (crossPC1-10, pooled) to control for residual
# population stratification, then re-z-score once across the pooled analysis sample
# for model coefficients (do not standardize within cohort or within region).
# Genotype batch is not adjusted (no batch variable in the harmonized data). All
# PGI-present individuals have complete crossPC, so this drops no one.
# resid_pgi_on_xpcs(dat, pgi): residuals of pgi ~ crossPC1-10 within `dat`
# (NA-preserving). Applied at EVERY PGI z-scoring site — df here and dat_cluster
# below — because each sample re-z-scores its own PGI to sd = 1.
.XPCS <- paste0("crossPC", 1:10)
resid_pgi_on_xpcs <- function(dat, pgi) {
  stopifnot(all(.XPCS %in% names(dat)))
  ok <- stats::complete.cases(dat[, c(pgi, .XPCS)])
  r  <- rep(NA_real_, nrow(dat))
  r[ok] <- stats::residuals(
    stats::lm(stats::as.formula(paste(pgi, "~", paste(.XPCS, collapse = " + "))),
              data = dat[ok, , drop = FALSE]))
  r
}
df$PGI_Edu_z    <- z_scale(resid_pgi_on_xpcs(df, "PGI_Edu"))
df$PGI_Cog_z    <- z_scale(resid_pgi_on_xpcs(df, "PGI_Cog"))
df$PGI_nonCog_z <- z_scale(resid_pgi_on_xpcs(df, "PGI_nonCog"))
df$PGI_Height_z <- z_scale(resid_pgi_on_xpcs(df, "PGI_Height"))

# edu_z_kernel is the PRIMARY attainment outcome (A1/A2/A3 and the S1/S2/S4
# attainment sensitivities). It is years of education standardized WITHIN
# BIRTH YEAR via a Gaussian kernel over birth_year (bandwidth = KERNEL_BW
# years; see Plan_deviations.md §2/§3): each person is expressed as their
# education relative to others born in the same birth-year neighbourhood.
# This is the cohort-relative-standing construct the GxE design targets, and
# it puts attainment on the same metric as the mobility components
# (mobility = edu_z_kernel - parental_edu_z_kernel).
df$edu_z_kernel <- kernel_z_score(df$education, df$birth_year,
                                  bandwidth = KERNEL_BW)
cat(sprintf("edu_z_kernel computed: N non-missing = %d / %d\n",
            sum(!is.na(df$edu_z_kernel)), nrow(df)))
# height_z, the outcome of the S2 height negative control, is standardized within
# birth year like edu_z_kernel above - same Gaussian kernel, same bandwidth, same
# effective-N floor - so the sanity model sits on the same kind of scale as the
# attainment models it is read against. RUN_20260729_1555 predates this and its
# cached height tables are on a global z.
if ("height" %in% names(df))
  df$height_z <- kernel_z_score(df$height, df$birth_year, bandwidth = KERNEL_BW)

# ---- 3j. Reunification variable (factor) ----
# Turned 15 in/after 1990 -> born >= 1975 -> "post"
df$reunif <- factor(
  ifelse(df$birth_year >= 1975, "post", "pre"),
  levels = c("pre", "post")
)

# Store by_mean as attribute for extraction functions
attr(df, "by_mean") <- by_mean

# ---- 3k. Kernel-standardized parental education ----
# Kernel regression: weight closer birth_year individuals more.
# Bandwidth = KERNEL_BW (default 5 years), pooled across all studies.

HAS_parental_education <- "parental_education" %in% names(df) &&
                          !all(is.na(df$parental_education))

if (HAS_parental_education) {

  # Parental education covariate: standardize using CHILD's birth_year
  # This captures "relative family advantage within cohort"
  cat("\nKernel-standardizing parental education (bandwidth =", KERNEL_BW, "years)...\n")
  df$ParEdu_z <- kernel_z_score(df$parental_education, df$birth_year,
                                bandwidth = KERNEL_BW)

  # Diagnostic: kernel-smoothed moments
  pe_moments <- kernel_moments(df$parental_education, df$birth_year,
                               bandwidth = KERNEL_BW)
  cat(sprintf("  Kernel standardization: %d unique birth years, ",
              sum(!is.na(pe_moments$mu))))
  cat(sprintf("mean range [%.1f, %.1f], SD range [%.1f, %.1f]\n",
              min(pe_moments$mu, na.rm = TRUE), max(pe_moments$mu, na.rm = TRUE),
              min(pe_moments$sd, na.rm = TRUE), max(pe_moments$sd, na.rm = TRUE)))

  # For mobility: z-score parental education within PARENT's birth cohort
  # Use actual parent birth years if available, otherwise approximate
  if (all(c("mother_birth_year", "father_birth_year") %in% names(df))) {
    df$parent_by <- rowMeans(
      cbind(safe_numeric(df$mother_birth_year),
            safe_numeric(df$father_birth_year)),
      na.rm = TRUE
    )
    # Fall back to approximation where parent birth years are missing
    df$parent_by <- ifelse(is.finite(df$parent_by), df$parent_by,
                           df$birth_year - PARENT_GENERATION_GAP)
  } else {
    df$parent_by <- df$birth_year - PARENT_GENERATION_GAP
  }
  df$parental_edu_z_kernel <- kernel_z_score(df$parental_education,
                                              df$parent_by,
                                              bandwidth = KERNEL_BW)

  # Mobility = child cohort-z - parental cohort-z
  df$mobility <- df$edu_z_kernel - df$parental_edu_z_kernel

  cat(sprintf("Mobility computed for N = %d (non-missing)\n", sum(!is.na(df$mobility))))
  cat(sprintf("parental_edu_z_kernel computed for N = %d (non-missing)\n", sum(!is.na(df$parental_edu_z_kernel))))

  # Kernel diagnostic plot, split by region (non-red/blue palette so it does
  # not clash with the East/West colour mapping used in the slope plots).
  pe_moments_reg <- bind_rows(
    kernel_moments(df$parental_education[df$east_west == "West"],
                   df$birth_year[df$east_west == "West"],
                   bandwidth = KERNEL_BW) |> mutate(Region = "West"),
    kernel_moments(df$parental_education[df$east_west == "East"],
                   df$birth_year[df$east_west == "East"],
                   bandwidth = KERNEL_BW) |> mutate(Region = "East")
  )
  PLOTS[["Kernel_ParEdu"]] <- ggplot(pe_moments_reg, aes(x = t, colour = Region, fill = Region)) +
    geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mu), linewidth = 0.9) +
    geom_point(aes(y = mu, size = n_raw), alpha = 0.30) +
    scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_fill_manual  (values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_size_continuous(range = c(0.3, 3), guide = "none") +
    labs(x = "Child's birth year",
         y = "Parental education (years)",
         title = "Kernel-smoothed parental education by birth year",
         subtitle = sprintf("Bandwidth = %d years, shaded = +/- 1 SD; one smooth per region.", KERNEL_BW)) +
    .theme

} else {
  df$mobility <- NA_real_; df$parental_edu_z_kernel <- NA_real_
  cat("Parental education not found - mobility analyses will be skipped.\n")
}

# Kernel-smoothed own (individual) education by birth year, split by region.
own_edu_moments_reg <- bind_rows(
  kernel_moments(df$education[df$east_west == "West"],
                 df$birth_year[df$east_west == "West"],
                 bandwidth = KERNEL_BW) |> mutate(Region = "West"),
  kernel_moments(df$education[df$east_west == "East"],
                 df$birth_year[df$east_west == "East"],
                 bandwidth = KERNEL_BW) |> mutate(Region = "East")
)
PLOTS[["Kernel_OwnEdu"]] <- ggplot(own_edu_moments_reg, aes(x = t, colour = Region, fill = Region)) +
  geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
  geom_line(aes(y = mu), linewidth = 0.9) +
  geom_point(aes(y = mu, size = n_raw), alpha = 0.30) +
  scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_fill_manual  (values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_size_continuous(range = c(0.3, 3), guide = "none") +
  labs(x = "Birth year",
       y = "Years of education",
       title = "Kernel-smoothed individual education by birth year",
       subtitle = sprintf("Bandwidth = %d years, shaded = +/- 1 SD; one smooth per region.", KERNEL_BW)) +
  .theme

# ---- 3l. Region-specific PGI interaction terms (for GAMs) ----
# Built from PGI_Edu_z (sd = 1 on the analytic sample) so the GAM-implied
# slopes are reported per 1 analytic-sample SD of PGI - the same scale as
# the linear models (A1/A2/B1/B2) and the np-overlay slopes. PGI_Edu (the
# pooled-merge-z) has sd ~ 0.71 on the analytic sample after exclusions,
# so using it in the GAMs would put the GAM slopes on a ~1/0.71 larger
# scale than the linear-model slopes.
df$PGI_W  <- df$PGI_Edu_z * (1L - df$east)
df$PGI_E  <- df$PGI_Edu_z *        df$east
df$PGI_WF <- df$PGI_Edu_z * (1L - df$east) * (1L - df$gender_01)
df$PGI_WM <- df$PGI_Edu_z * (1L - df$east) *        df$gender_01
df$PGI_EF <- df$PGI_Edu_z *        df$east  * (1L - df$gender_01)
df$PGI_EM <- df$PGI_Edu_z *        df$east  *        df$gender_01

# ---- 3m. Analysis subsets ----
dat_edu <- df   # full analytic sample
dat_mob <- df |> filter(!is.na(mobility), !is.na(parental_edu_z_kernel))
attr(dat_mob, "by_mean") <- by_mean

EDU_SD <- sd(dat_edu$education, na.rm = TRUE)
MOB_SD <- if (nrow(dat_mob) > 1) sd(dat_mob$mobility, na.rm = TRUE) else NA_real_

# edu_std: education rescaled by a single global SD, available as an optional
# global-SD sensitivity. It is NOT the focal attainment outcome (that is
# edu_z_kernel; see Plan_deviations.md §2/§3): a global multiplicative rescale
# does not remove the cross-birth-year/region differences in education SD that
# the kernel z absorbs.
dat_edu$edu_std <- dat_edu$education / EDU_SD
dat_mob$edu_std <- dat_mob$education / EDU_SD

mobility_attrition_steps <- list(
  list(label = "Parental mobility observed", keep = function(d) !is.na(d$mobility)),
  list(label = "Parental covariate observed", keep = function(d) !is.na(d$parental_edu_z_kernel))
)
attrition_mob <- compute_attrition_flow(dat_edu, mobility_attrition_steps, pipeline = "Mobility sample")
RESULTS[["Attrition_Mobility"]] <- attrition_mob
RESULTS[["Attrition_All"]] <- bind_rows(attrition_core, attrition_mob)
# Split into two plots so each pipeline can have its own y-axis cap.
# The Core sample shares one y range across studies, capped at 12500;
# TwinLife's first one-two pre-exclusion columns sit above that cap and are
# clipped on purpose so the smaller-step changes in the analytic-sample range
# stay readable.
PLOTS[["Attrition_flow_Core"]] <- plot_attrition_flow(
  attrition_core,
  title = "Sample attrition by exclusion step and study - Core sample",
  y_max = 12500
)
PLOTS[["Attrition_flow_Mob"]] <- plot_attrition_flow(
  attrition_mob,
  title = "Sample attrition by exclusion step and study - Mobility sample"
)

cat(sprintf("\nSection A sample (dat_edu): N = %d\n", nrow(dat_edu)))
cat(sprintf("Section B sample (dat_mob): N = %d\n", nrow(dat_mob)))


# =====================================================================
# (b) Full-family analytic sample (dat_cluster) for the primary RE track.
#     Same data prep on the sample that keeps all family adults, plus fid_re.
#     Every focal model is fit on dat_cluster with + s(fid_re, bs = 're') in
#     mgcv; see Random_effects_family_clustering.md. The dedup'd dat_edu /
#     dat_mob samples built above are the supplementary sensitivity track.
# =====================================================================
if (any(df_raw$cohort == "TwinLife")) {
  cat("Building dat_cluster (full-family TwinLife) for the primary RE pipeline...\n")

  df_with_twins <- df_raw
  if (gender_col == "Gender" && !"gender" %in% names(df_with_twins)) {
    names(df_with_twins)[names(df_with_twins) == "Gender"] <- "gender"
  }
  if ("Height" %in% names(df_with_twins) && !"height" %in% names(df_with_twins)) {
    names(df_with_twins)[names(df_with_twins) == "Height"] <- "height"
  }
  df_with_twins$birth_year   <- safe_numeric(df_with_twins$birth_year)
  df_with_twins$education    <- safe_numeric(df_with_twins$education)
  df_with_twins$PGI_Edu     <- safe_numeric(df_with_twins$PGI_Edu)
  df_with_twins$PGI_Cog     <- safe_numeric(df_with_twins$PGI_Cog)
  df_with_twins$PGI_nonCog  <- safe_numeric(df_with_twins$PGI_nonCog)
  df_with_twins$PGI_Height  <- safe_numeric(df_with_twins$PGI_Height)
  if ("parental_education" %in% names(df_with_twins))
    df_with_twins$parental_education <- safe_numeric(df_with_twins$parental_education)

  ew_lower_tw <- tolower(trimws(as.character(df_with_twins$east_west)))
  df_with_twins$east_west <- ifelse(ew_lower_tw == "west", "West",
                                    ifelse(ew_lower_tw == "east", "East", NA_character_))
  df_with_twins$east_west <- factor(df_with_twins$east_west, levels = c("West", "East"))
  gen_raw_tw <- tolower(trimws(as.character(df_with_twins$gender)))
  df_with_twins$gender <- dplyr::case_when(
    gen_raw_tw %in% c("female", "0") ~ "female",
    gen_raw_tw %in% c("male",   "1") ~ "male",
    TRUE ~ NA_character_)
  df_with_twins$gender <- factor(df_with_twins$gender, levels = c("female", "male"))
  df_with_twins$cohort <- factor(df_with_twins$cohort)

  if ("ptyp" %in% names(df_with_twins)) {
    df_with_twins$ptyp <- safe_numeric(df_with_twins$ptyp)
    df_with_twins$is_twin_family <- ifelse(
      df_with_twins$cohort == "TwinLife" & df_with_twins$ptyp %in% c(1, 2, 200, 201),
      TRUE, FALSE)
  } else {
    df_with_twins$is_twin_family <- FALSE
  }
  if (!"fid" %in% names(df_with_twins)) df_with_twins$fid <- NA_character_
  df_with_twins$fid <- ifelse(df_with_twins$cohort == "TwinLife",
                              as.character(df_with_twins$fid), NA_character_)

  # Drop twins / non-twin siblings (ptyp 1, 2, 200, 201) from both tracks
  # (children, mostly under age-25 cutoff). Spousal pairs and parents are
  # retained and handled by the RE penalty.
  dat_cluster <- df_with_twins |>
    dplyr::filter(
      !is.na(birth_year), birth_year <= 1996,
      !is.na(education), !is.na(east_west),
      !is.na(gender), !is.na(PGI_Edu),
      abs(PGI_Edu) <= MAX_ABS_PGI_EDU,
      !(cohort == "SHIP" & east_west == "West"),
      !is_twin_family
    )
  cat(sprintf("  dat_cluster (full-family adults): N = %d\n", nrow(dat_cluster)))

  # Variable prep on dat_cluster: same scaling and helper columns as on
  # dat_edu, plus the fid_re grouping factor for the RE term.
  dat_cluster$BYc      <- dat_cluster$birth_year - mean(dat_cluster$birth_year, na.rm = TRUE)
  dat_cluster$BYc_z    <- z_scale(dat_cluster$birth_year)
  # Residualise on cross-cohort PCs before z-scoring, as for df above.
  # dat_cluster_mob inherits these columns, so Section B is covered too.
  dat_cluster$PGI_Edu_z   <- z_scale(resid_pgi_on_xpcs(dat_cluster, "PGI_Edu"))
  dat_cluster$PGI_Cog_z   <- z_scale(resid_pgi_on_xpcs(dat_cluster, "PGI_Cog"))
  dat_cluster$PGI_nonCog_z<- z_scale(resid_pgi_on_xpcs(dat_cluster, "PGI_nonCog"))
  dat_cluster$PGI_Height_z<- z_scale(resid_pgi_on_xpcs(dat_cluster, "PGI_Height"))
  dat_cluster$edu_z_kernel <- kernel_z_score(dat_cluster$education, dat_cluster$birth_year,
                                             bandwidth = KERNEL_BW)
  dat_cluster$east_west_c <- ifelse(dat_cluster$east_west == "East",   0.5, -0.5)
  dat_cluster$gender_c    <- ifelse(dat_cluster$gender    == "female", 0.5, -0.5)
  dat_cluster$east        <- as.integer(dat_cluster$east_west == "East")
  dat_cluster$gender_01   <- as.integer(dat_cluster$gender == "male")
  dat_cluster$edu_std     <- dat_cluster$education / EDU_SD

  # Random-effects grouping factor. TwinLife rows with a family ID get
  # their family-level intercept; every other row (BASE-II / SHIP / SOEP,
  # plus rare TwinLife rows without fid) is pooled into a single "non_TL"
  # level so the random-effects penalty does not have to estimate ~10000
  # singleton intercepts. The "non_TL" BLUP absorbs only an overall mean
  # shift, which does not affect the focal slopes.
  # See Random_effects_family_clustering.md section 3 for the rationale.
  dat_cluster$fid_re <- factor(ifelse(
    dat_cluster$cohort == "TwinLife" & !is.na(dat_cluster$fid),
    paste0("TL_", dat_cluster$fid),
    "non_TL"
  ))

  # ----------------------------------------------------------------------
  # Derived columns that make dat_cluster self-sufficient for A_analysis.R /
  # B_analysis.R, including their gender and sensitivity stages.
  # ----------------------------------------------------------------------

  # Reunification step-function predictor (used by A2 RE).
  dat_cluster$reunif <- factor(
    ifelse(dat_cluster$birth_year >= 1975, "post", "pre"),
    levels = c("pre", "post")
  )

  # Region- and gender-specific PGI interaction columns. Built from
  # PGI_Edu_z (sd = 1 on dat_cluster) so the GAM-implied slopes are
  # reported per 1 analytic-sample SD of PGI.
  dat_cluster$PGI_W  <- dat_cluster$PGI_Edu_z * (1L - dat_cluster$east)
  dat_cluster$PGI_E  <- dat_cluster$PGI_Edu_z *        dat_cluster$east
  dat_cluster$PGI_WF <- dat_cluster$PGI_Edu_z * (1L - dat_cluster$east) * (1L - dat_cluster$gender_01)
  dat_cluster$PGI_WM <- dat_cluster$PGI_Edu_z * (1L - dat_cluster$east) *        dat_cluster$gender_01
  dat_cluster$PGI_EF <- dat_cluster$PGI_Edu_z *        dat_cluster$east  * (1L - dat_cluster$gender_01)
  dat_cluster$PGI_EM <- dat_cluster$PGI_Edu_z *        dat_cluster$east  *        dat_cluster$gender_01

  # Centring constant as an attribute so extract_pgi_slope() /
  # diff_smooth_simul() can recover original birth years from BYc.
  attr(dat_cluster, "by_mean") <- mean(dat_cluster$birth_year, na.rm = TRUE)

  # ----------------------------------------------------------------------
  # dat_cluster_mob: the full-family RE mobility sample (Section B).
  # mobility = child cohort-z education - parental cohort-z education;
  # parental education is z-scored within the PARENT birth cohort (kernel-z
  # over parent_by). SHIP has no parental education, so it drops out here.
  # ----------------------------------------------------------------------
  if ("parental_education" %in% names(dat_cluster)) {
    dat_cluster$parental_education <- safe_numeric(dat_cluster$parental_education)
    if (all(c("mother_birth_year", "father_birth_year") %in% names(dat_cluster))) {
      .pby <- rowMeans(cbind(safe_numeric(dat_cluster$mother_birth_year),
                             safe_numeric(dat_cluster$father_birth_year)),
                       na.rm = TRUE)
      dat_cluster$parent_by <- ifelse(is.finite(.pby), .pby,
                                      dat_cluster$birth_year - PARENT_GENERATION_GAP)
    } else {
      dat_cluster$parent_by <- dat_cluster$birth_year - PARENT_GENERATION_GAP
    }
    dat_cluster$parental_edu_z_kernel <- kernel_z_score(
      dat_cluster$parental_education, dat_cluster$parent_by, bandwidth = KERNEL_BW)
    dat_cluster$mobility <- dat_cluster$edu_z_kernel - dat_cluster$parental_edu_z_kernel
  } else {
    dat_cluster$parental_edu_z_kernel <- NA_real_
    dat_cluster$mobility <- NA_real_
  }

  dat_cluster_mob <- dat_cluster |>
    dplyr::filter(!is.na(mobility), !is.na(parental_edu_z_kernel))
  # Re-residualise + z-score the PGIs WITHIN the mobility sample (SHIP excluded):
  # the cross-cohort ancestry axis shifts once SHIP is dropped, so inheriting
  # dat_cluster's full-sample residual leaves East-West PGI variance uncontrolled.
  # Residualise on this sample's own crossPC, then rebuild the region/gender PGI
  # interaction columns (PGI_W/E/WF/WM/EF/EM) from the re-standardised PGI.
  dat_cluster_mob$PGI_Edu_z    <- z_scale(resid_pgi_on_xpcs(dat_cluster_mob, "PGI_Edu"))
  dat_cluster_mob$PGI_Cog_z    <- z_scale(resid_pgi_on_xpcs(dat_cluster_mob, "PGI_Cog"))
  dat_cluster_mob$PGI_nonCog_z <- z_scale(resid_pgi_on_xpcs(dat_cluster_mob, "PGI_nonCog"))
  dat_cluster_mob$PGI_Height_z <- z_scale(resid_pgi_on_xpcs(dat_cluster_mob, "PGI_Height"))
  dat_cluster_mob$PGI_W  <- dat_cluster_mob$PGI_Edu_z * (1L - dat_cluster_mob$east)
  dat_cluster_mob$PGI_E  <- dat_cluster_mob$PGI_Edu_z *        dat_cluster_mob$east
  dat_cluster_mob$PGI_WF <- dat_cluster_mob$PGI_Edu_z * (1L - dat_cluster_mob$east) * (1L - dat_cluster_mob$gender_01)
  dat_cluster_mob$PGI_WM <- dat_cluster_mob$PGI_Edu_z * (1L - dat_cluster_mob$east) *        dat_cluster_mob$gender_01
  dat_cluster_mob$PGI_EF <- dat_cluster_mob$PGI_Edu_z *        dat_cluster_mob$east  * (1L - dat_cluster_mob$gender_01)
  dat_cluster_mob$PGI_EM <- dat_cluster_mob$PGI_Edu_z *        dat_cluster_mob$east  *        dat_cluster_mob$gender_01
  # Region-specific ParEdu cell columns (parallel to PGI_W/PGI_E) so the
  # B3 GAM can give parental education the same region-specific smooth
  # flexibility as PGI (see Plan_deviations.md — B3 ParEdu control).
  dat_cluster_mob$ParEdu_W <- dat_cluster_mob$parental_edu_z_kernel * (1L - dat_cluster_mob$east)
  dat_cluster_mob$ParEdu_E <- dat_cluster_mob$parental_edu_z_kernel *        dat_cluster_mob$east
  attr(dat_cluster_mob, "by_mean") <- attr(dat_cluster, "by_mean")
  cat(sprintf("dat_cluster_mob (full-family mobility RE sample): N = %d\n",
              nrow(dat_cluster_mob)))
}

