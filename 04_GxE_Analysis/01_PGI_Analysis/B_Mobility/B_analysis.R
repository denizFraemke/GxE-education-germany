# B_analysis.R
# Section B entry point. Loads data, builds samples, computes mobility
# descriptives, fits B0/B0_pre/B1/B2/B3 and the gender models B4/B5/B6/B0g
# (family RE on TwinLife families via bam), runs leave-one-study-out, saves
# the models RDS.
#
# ParEdu control (see ../../../Plan_deviations.md §9):
#   B1/B2 -> full preregistered Eq-4 ParEdu interaction set
#   B3    -> ParEdu main + ParEdu:region + region-specific ParEdu smooths
#            over birth year (parallel to the PGI smooths)
#
# Usage (from run_B.command): Rscript B_Mobility/B_analysis.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
})

# ---- Resolve paths (works under Rscript and source()) ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine B_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR  # flattened: 01_PGI_Analysis is the workflow + script root
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL <- get_run_ts()
OUT_DIR_B  <- output_dir("B_Mobility")

source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()

DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var
DATA_DIR        <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")

source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

if (!exists("dat_cluster_mob") || nrow(dat_cluster_mob) < 200) {
  stop("dat_cluster_mob missing or < 200 rows — cannot run Section B.")
}
cat(sprintf("B_analysis: dat_cluster_mob N = %d\n", nrow(dat_cluster_mob)))

source(file.path(SECTION_DIR, "B_descriptives.R"), local = TRUE)

re_term <- "+ s(fid_re, bs = 're')"

# ---- ParEdu interaction strings (preregistered Eq-4) ----
paredu_b1 <- paste(
  "+ parental_edu_z_kernel",
  "+ parental_edu_z_kernel:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z",
  "+ parental_edu_z_kernel:BYc_z",
  "+ parental_edu_z_kernel:BYc_z:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z:BYc_z",
  "+ parental_edu_z_kernel:PGI_Edu_z:east_west_c")
paredu_b2 <- paste(
  "+ parental_edu_z_kernel",
  "+ parental_edu_z_kernel:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z",
  "+ parental_edu_z_kernel:reunif",
  "+ parental_edu_z_kernel:reunif:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z:reunif",
  "+ parental_edu_z_kernel:PGI_Edu_z:east_west_c")

# ---- B1: linear focal (Eq-4 ParEdu) ----
cat("B1 (RE)...\n")
m_B1_re <- mgcv::bam(
  as.formula(paste(
    "mobility ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3",
    paredu_b1, re_term)),
  data = dat_cluster_mob, method = "fREML", discrete = TRUE)

# ---- B2: step-function focal (Eq-4 ParEdu, reunif) ----
cat("B2 (RE)...\n")
m_B2_re <- mgcv::bam(
  as.formula(paste(
    "mobility ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3",
    paredu_b2, re_term)),
  data = dat_cluster_mob, method = "fREML", discrete = TRUE)

# ---- B3: GAM null / linear / nonlinear, with ParEdu control + smooths ----
# ParEdu control common to all three B3 levels: main + region + the
# region-specific ParEdu smooths over birth year (../../../Plan_deviations.md §9).
b3_formula <- function(level, basis, rt = re_term) {
  pe_smooth <- sprintf(
    "s(BYc, by = ParEdu_W, k = %d, bs = '%s') + s(BYc, by = ParEdu_E, k = %d, bs = '%s')",
    K_DEFAULT, basis, K_DEFAULT, basis)
  # Baseline carries the PGI main-effect cells (PGI_W + PGI_E, constant
  # slope) so baseline -> linear isolates the region-specific cohort slopes.
  rhs <- paste("mobility ~ BYc + east_west_c + gender_c +",
               "parental_edu_z_kernel + parental_edu_z_kernel:east_west_c +",
               pe_smooth, "+ PGI_W + PGI_E")
  if (level %in% c("lin", "nonlin"))
    rhs <- paste(rhs, "+ PGI_W:BYc + PGI_E:BYc")
  if (level == "nonlin")
    rhs <- paste(rhs, sprintf(
      "+ s(BYc, by = PGI_W, k = %d, bs = '%s') + s(BYc, by = PGI_E, k = %d, bs = '%s')",
      K_DEFAULT, basis, K_DEFAULT, basis))
  paste(rhs, rt)
}

cat("B3 null (RE)...\n")
m_B3_null_re <- mgcv::bam(as.formula(b3_formula("null", BS_DEFAULT)),
                         data = dat_cluster_mob, method = "fREML", discrete = TRUE)
cat("B3 linear (RE)...\n")
m_B3_lin_re  <- mgcv::bam(as.formula(b3_formula("lin", BS_DEFAULT)),
                         data = dat_cluster_mob, method = "fREML", discrete = TRUE)
cat("B3 nonlinear tp (RE)...\n")
m_B3_nonlin_re <- mgcv::bam(as.formula(b3_formula("nonlin", BS_DEFAULT)),
                         data = dat_cluster_mob, method = "fREML", discrete = TRUE)
cat("B3 nonlinear cr sensitivity (RE)...\n")
m_B3_cr_re   <- mgcv::bam(as.formula(b3_formula("nonlin", BS_SENS)),
                         data = dat_cluster_mob, method = "fREML", discrete = TRUE)

# Guard: PGI main effect in the baseline (constant slope) so the AIC/BIC
# ladder is honest; B3's focal test is non-linearity only (linear ->
# nonlinear). The linear PGI cohort trend is tested in B1 (Wald).
assert_pgi_in_baseline(m_B3_null_re, m_B3_lin_re, c("PGI_W", "PGI_E"), "B3")

# ---- B0 / B0_pre: overall + pre-reunification East-West contrast (exploratory) ----
b0_f_base <- paste(
  "mobility ~ PGI_Edu_z * east_west_c +",
  "east_west_c * (gender_c + BYc_z + parental_edu_z_kernel) +",
  "gender_c * BYc_z + gender_c * parental_edu_z_kernel +",
  "BYc_z * parental_edu_z_kernel")
b0_f <- paste(b0_f_base, re_term)
cat("B0 (RE)...\n")
m_B0_re <- mgcv::bam(as.formula(b0_f), data = dat_cluster_mob,
                     method = "fREML", discrete = TRUE)
dat_cluster_mob_pre <- dat_cluster_mob |> dplyr::filter(reunif == "pre")
m_B0_pre_re <- NULL
if (nrow(dat_cluster_mob_pre) >= 100 &&
    length(unique(dat_cluster_mob_pre$east_west)) == 2) {
  cat("B0 pre-reunification (RE)...\n")
  m_B0_pre_re <- mgcv::bam(as.formula(b0_f), data = dat_cluster_mob_pre,
                           method = "fREML", discrete = TRUE)
}

# ---- Coefficient / smooth / fit tables ----
.drop_re_smooths <- function(sm) {
  if (is.null(sm) || !nrow(sm)) return(sm)
  sm[!grepl("fid_re", sm$smooth_term), , drop = FALSE]
}
B1_Coefficients <- tidy_coefs(m_B1_re, "B1_RE")
B1_ModelFit     <- model_fit_row(m_B1_re, "B1_RE")
B2_Coefficients <- tidy_coefs(m_B2_re, "B2_RE")
B2_ModelFit     <- model_fit_row(m_B2_re, "B2_RE")

# B3 nested LRT: null->linear (ML), linear->nonlinear (fREML).
.b3_ln <- lrt_pair(m_B3_lin_re,  m_B3_nonlin_re)   # linear -> nonlinear (fREML)
.gf <- function(m) c(AIC = AIC(m), BIC = BIC(m),
                     logLik = as.numeric(logLik(m)), dev = summary(m)$dev.expl)
.f0 <- .gf(m_B3_null_re); .f1 <- .gf(m_B3_lin_re); .f2 <- .gf(m_B3_nonlin_re)
B3_NestedLRT <- data.frame(
  model      = c("B3_baseline", "B3_lin", "B3_nonlin"),
  AIC        = c(.f0["AIC"], .f1["AIC"], .f2["AIC"]),
  BIC        = c(.f0["BIC"], .f1["BIC"], .f2["BIC"]),
  logLik     = c(.f0["logLik"], .f1["logLik"], .f2["logLik"]),
  dev_expl   = c(.f0["dev"], .f1["dev"], .f2["dev"]),
  LRT_chisq  = c(NA_real_, NA_real_, .b3_ln[["chisq"]]),
  LRT_df     = c(NA_real_, NA_real_, .b3_ln[["df"]]),
  LRT_p      = c(NA_real_, NA_real_, .b3_ln[["p"]]),
  comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE
)

# Smooth tables (drop RE). Keep PGI + ParEdu smooths for the sheet; the
# GAM-overview nonlinearity_conclusion uses the PGI smooths only.
B3_Smooths_tp <- .drop_re_smooths(tidy_smooths(m_B3_nonlin_re, "B3_RE_nonlin_tp"))
B3_Smooths_cr <- .drop_re_smooths(tidy_smooths(m_B3_cr_re,     "B3_RE_nonlin_cr"))
B3_Smooths_pgi_tp <- B3_Smooths_tp[grepl("PGI", B3_Smooths_tp$smooth_term), , drop = FALSE]

B_ModelFit <- dplyr::bind_rows(
  B1_ModelFit, B2_ModelFit,
  model_fit_row(m_B3_null_re,   "B3_baseline_RE"),
  model_fit_row(m_B3_lin_re,    "B3_lin_RE"),
  model_fit_row(m_B3_nonlin_re, "B3_nonlin_RE_tp"),
  model_fit_row(m_B3_cr_re,     "B3_nonlin_RE_cr"),
  model_fit_row(m_B0_re,        "B0_RE")
)
if (!is.null(m_B0_pre_re))
  B_ModelFit <- dplyr::bind_rows(B_ModelFit, model_fit_row(m_B0_pre_re, "B0_pre_RE"))

# k.check on the B3 GAM smooths (drop RE row).
k_index_rows <- list()
for (lbl in c("B3_nonlin_RE_tp", "B3_nonlin_RE_cr")) {
  mod <- if (lbl == "B3_nonlin_RE_tp") m_B3_nonlin_re else m_B3_cr_re
  kc <- tryCatch(mgcv::k.check(mod), error = function(e) NULL)
  if (!is.null(kc)) {
    df_kc <- data.frame(model = lbl, smooth_term = rownames(kc),
                        k_prime = kc[, "k'"], edf = kc[, "edf"],
                        k_index = kc[, "k-index"], p_value = kc[, "p-value"],
                        stringsAsFactors = FALSE)
    k_index_rows[[lbl]] <- df_kc[!grepl("fid_re", df_kc$smooth_term), , drop = FALSE]
  }
}
B3_K_Index <- if (length(k_index_rows) > 0) dplyr::bind_rows(k_index_rows) else NULL

# B0 / B0_pre region-specific PGI-mobility slopes + interaction coef.
B0_RegionSlopes     <- extract_lm_region_pgi_effects(m_B0_re, dat_cluster_mob)
B0_Coefficients     <- tidy_coefs(m_B0_re, "B0_RE")
B0_pre_RegionSlopes <- if (!is.null(m_B0_pre_re))
  extract_lm_region_pgi_effects(m_B0_pre_re, dat_cluster_mob_pre) else NULL
B0_pre_Coefficients <- if (!is.null(m_B0_pre_re))
  tidy_coefs(m_B0_pre_re, "B0_pre_RE") else NULL

# ---- Leave-one-study-out (folds: studies present in the mobility sample) ----
cat("LOO: leave-one-study-out refits (B0, B1, B2, B3)...\n")
.loo_re_term <- function(sub) {
  if (!"fid_re" %in% names(sub)) return("")
  sub$fid_re <- droplevels(factor(sub$fid_re))
  if (nlevels(sub$fid_re) >= 2) "+ s(fid_re, bs = 're')" else ""
}

loo_fit_B1 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(as.formula(paste(
    "mobility ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3", paredu_b1, rt)),
    data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "B1"); r <- r[r$term == "PGI_Edu_z:BYc_z:east_west_c", ]
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5), p = round(r$p.value, 5))
}
B1_LOO <- build_loo_table(leave_one_cohort_out(dat_cluster_mob, loo_fit_B1),
                          estimate_col = "estimate")

loo_fit_B2 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(as.formula(paste(
    "mobility ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3", paredu_b2, rt)),
    data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "B2"); r <- r[r$term == "PGI_Edu_z:reunifpost:east_west_c", ]
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5), p = round(r$p.value, 5))
}
B2_LOO <- build_loo_table(leave_one_cohort_out(dat_cluster_mob, loo_fit_B2),
                          estimate_col = "estimate")

# B0 / B0_pre LOO: overall (and pre-reunification) East-West PGI contrast.
loo_fit_B0 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- tryCatch(mgcv::bam(as.formula(paste(b0_f_base, rt)),
                          data = sub, method = "fREML", discrete = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(c(estimate = NA_real_, SE = NA_real_, p = NA_real_))
  r <- tidy_coefs(m, "B0"); r <- r[r$term == "PGI_Edu_z:east_west_c", ]
  if (!nrow(r)) return(c(estimate = NA_real_, SE = NA_real_, p = NA_real_))
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5), p = round(r$p.value, 5))
}
B0_LOO <- build_loo_table(leave_one_cohort_out(dat_cluster_mob, loo_fit_B0),
                          estimate_col = "estimate")
B0_pre_LOO <- if (!is.null(m_B0_pre_re))
  build_loo_table(leave_one_cohort_out(dat_cluster_mob_pre, loo_fit_B0),
                  estimate_col = "estimate") else NULL

# B0 LOO per-fold region slopes: West/East PGI-mobility slope per held-out
# study (for the per-fold LOO plot, the leave-one-study-out analogue of
# the main B0 region-slope figure).
.b0_loo_slopes <- function(dat, analysis_label) {
  folds <- c(list("Full sample" = character(0)), cohort_groups(dat))
  out <- list()
  for (fl in names(folds)) {
    dl  <- folds[[fl]]
    sub <- if (length(dl) == 0) dat else
      dat[!(as.character(dat$cohort) %in% dl), , drop = FALSE]
    if (length(unique(sub$east_west)) < 2) next
    sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
    m <- tryCatch(mgcv::bam(as.formula(paste(b0_f_base, rt)), data = sub,
                            method = "fREML", discrete = TRUE), error = function(e) NULL)
    if (is.null(m)) next
    rs <- extract_lm_region_pgi_effects(m, sub)
    rs$fold     <- if (length(dl) == 0) "Full sample" else paste0("-", fl)
    rs$analysis <- analysis_label
    out[[fl]] <- rs
  }
  dplyr::bind_rows(out)
}
B0_LOO_slopes <- dplyr::bind_rows(
  .b0_loo_slopes(dat_cluster_mob, "Overall"),
  if (!is.null(m_B0_pre_re)) .b0_loo_slopes(dat_cluster_mob_pre, "Pre-reunification") else NULL)

# B3 LOO: nonlinearity LRT + per-fold slope curves (one pass).
b3_loo_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster_mob))
b3_loo_rows <- list(); b3_loo_curves <- list()
for (fl in names(b3_loo_folds)) {
  drop_levs <- b3_loo_folds[[fl]]
  sub <- if (length(drop_levs) == 0) dat_cluster_mob else
    dat_cluster_mob[!(as.character(dat_cluster_mob$cohort) %in% drop_levs), , drop = FALSE]
  attr(sub, "by_mean") <- attr(dat_cluster_mob, "by_mean")
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  fold_lab <- if (length(drop_levs) == 0) "Full sample" else paste0("-", fl)
  held     <- if (length(drop_levs) == 0) "(none)" else fl
  m_lin <- mgcv::bam(as.formula(b3_formula("lin",    BS_DEFAULT, rt)), data = sub, method = "fREML", discrete = TRUE)
  m_nl  <- mgcv::bam(as.formula(b3_formula("nonlin", BS_DEFAULT, rt)), data = sub, method = "fREML", discrete = TRUE)
  v <- lrt_pair(m_lin, m_nl)
  concl <- if (is.na(v[["p"]])) "Not tested"
           else if (v[["p"]] < 0.05) "Nonlinearity supported"
           else "No evidence that nonlinear specification is needed"
  b3_loo_rows[[fl]] <- data.frame(
    fold = fold_lab, held_out = held, N = nrow(sub),
    LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1),
    LRT_p = signif(v[["p"]], 4), nonlinearity_conclusion = concl,
    status = "OK", stringsAsFactors = FALSE)
  cv <- dplyr::bind_rows(
    extract_pgi_slope(m_nl, sub, east_west = "West"),
    extract_pgi_slope(m_nl, sub, east_west = "East"))
  cv$fold <- fold_lab
  b3_loo_curves[[fl]] <- cv
}
B3_LOO        <- dplyr::bind_rows(b3_loo_rows)
B3_LOO_curves <- dplyr::bind_rows(b3_loo_curves)
B3_LOO_curves$fold <- factor(B3_LOO_curves$fold, levels = unique(B3_LOO_curves$fold))
cat("LOO done.\n")

# =====================================================================
# GENDER as a main moderator of mobility (Section B, family random effects).
# B4/B6/B0g are the Section C mobility analyses of the preregistration
# (C1/C3/C0), estimated here with family random effects; B5 is the per-gender
# mobility GAM (analogue of A5). Parental education is controlled at the same
# flexibility as the focal PGI cells (../../../Plan_deviations.md §10).
# =====================================================================
.sm <- function(by, basis = BS_DEFAULT) sprintf("s(BYc, by = %s, k = %d, bs = '%s')", by, K_DEFAULT, basis)

# Gender / cell columns on the mobility sample (PGI_W/E/WF/WM/EF/EM and
# ParEdu_W/E are inherited from dat_cluster_mob; add the gender + cell ParEdu).
dat_cluster_mob$PGI_F    <- dat_cluster_mob$PGI_Edu_z * (1L - dat_cluster_mob$gender_01)
dat_cluster_mob$PGI_M    <- dat_cluster_mob$PGI_Edu_z *        dat_cluster_mob$gender_01
dat_cluster_mob$ParEdu_F <- dat_cluster_mob$parental_edu_z_kernel * (1L - dat_cluster_mob$gender_01)
dat_cluster_mob$ParEdu_M <- dat_cluster_mob$parental_edu_z_kernel *        dat_cluster_mob$gender_01
dat_cluster_mob$ParEdu_WF <- dat_cluster_mob$parental_edu_z_kernel * (1L - dat_cluster_mob$east) * (1L - dat_cluster_mob$gender_01)
dat_cluster_mob$ParEdu_WM <- dat_cluster_mob$parental_edu_z_kernel * (1L - dat_cluster_mob$east) *        dat_cluster_mob$gender_01
dat_cluster_mob$ParEdu_EF <- dat_cluster_mob$parental_edu_z_kernel *        dat_cluster_mob$east  * (1L - dat_cluster_mob$gender_01)
dat_cluster_mob$ParEdu_EM <- dat_cluster_mob$parental_edu_z_kernel *        dat_cluster_mob$east  *        dat_cluster_mob$gender_01

# ---- B4: linear four-way (RE), four-way ParEdu control ----
cat("B4 four-way (RE)...\n")
B4_focal <- "PGI_Edu_z:BYc_z:east_west_c:gender_c"
.b4_rhs  <- paste("mobility ~ PGI_Edu_z * BYc_z * east_west_c * gender_c +",
                  "parental_edu_z_kernel * BYc_z * east_west_c * gender_c")
m_B4_re <- mgcv::bam(as.formula(paste(.b4_rhs, re_term)),
                     data = dat_cluster_mob, method = "fREML", discrete = TRUE)
B4_Coefficients <- tidy_coefs(m_B4_re, "B4_RE")
B4_ModelFit     <- model_fit_row(m_B4_re, "B4_RE")
B4_LOO <- build_loo_table(leave_one_cohort_out(dat_cluster_mob, function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(as.formula(paste(.b4_rhs, rt)), data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "B4_RE"); r <- r[r$term == B4_focal, ]
  if (!nrow(r)) return(c(estimate = NA_real_, SE = NA_real_, p = NA_real_))
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5), p = round(r$p.value, 5))
}), estimate_col = "estimate")
# B4_LOO_gender: leave-one-study-out of the POOLED PGI x gender mobility
# interaction (PGI_Edu_z:gender_c) from the same B4 model (four-way ParEdu
# control preserved); B4_LOO tracks only the four-way term. Finding-matched to
# the reported PGI x gender mobility result; additional sensitivity
# (../../../Plan_deviations.md §14). Estimate, SE, 95% CI, exact p, supported
# birth-year range, delta-vs-full, convergence.
B4_LOO_gender <- loco_pooled_coef(
  dat_cluster_mob, .b4_rhs, "PGI_Edu_z:gender_c", "B4_gender_RE")

# ---- B5: per-gender PGI-mobility GAM, regions pooled (RE) ----
# Gender analog of A3/A5 for mobility; gender-specific ParEdu control at the
# same flexibility as the focal PGI cells.
cat("B5 per-gender pooled GAM (RE)...\n")
.b5_rhs <- function(level, k, basis) {
  pe  <- paste("parental_edu_z_kernel + parental_edu_z_kernel:gender_c +",
               .sm("ParEdu_F", basis), "+", .sm("ParEdu_M", basis))
  # Baseline carries the PGI main-effect cells (PGI_F + PGI_M) so
  # baseline -> linear isolates the per-gender cohort slopes.
  rhs <- paste("mobility ~ BYc + east_west_c + gender_c +", pe, "+ PGI_F + PGI_M")
  if (level %in% c("lin", "nonlin")) rhs <- paste(rhs, "+ PGI_F:BYc + PGI_M:BYc")
  if (level == "nonlin") rhs <- paste(rhs, "+", .sm("PGI_F", basis), "+", .sm("PGI_M", basis))
  rhs
}
m_B5_null_re   <- mgcv::bam(as.formula(paste(.b5_rhs("null",   K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B5_lin_re    <- mgcv::bam(as.formula(paste(.b5_rhs("lin",    K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B5_nonlin_re <- mgcv::bam(as.formula(paste(.b5_rhs("nonlin", K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B5_cr_re     <- mgcv::bam(as.formula(paste(.b5_rhs("nonlin", K_SENS,    BS_SENS),    re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
assert_pgi_in_baseline(m_B5_null_re, m_B5_lin_re, c("PGI_F", "PGI_M"), "B5")
.b5_ln <- lrt_pair(m_B5_lin_re, m_B5_nonlin_re)   # linear -> nonlinear (fREML)
.b5_0 <- .gf(m_B5_null_re); .b5_1 <- .gf(m_B5_lin_re); .b5_2 <- .gf(m_B5_nonlin_re)
B5_NestedLRT <- data.frame(model = c("B5_baseline", "B5_lin", "B5_nonlin"),
  AIC = c(.b5_0["AIC"], .b5_1["AIC"], .b5_2["AIC"]), BIC = c(.b5_0["BIC"], .b5_1["BIC"], .b5_2["BIC"]),
  logLik = c(.b5_0["logLik"], .b5_1["logLik"], .b5_2["logLik"]), dev_expl = c(.b5_0["dev"], .b5_1["dev"], .b5_2["dev"]),
  LRT_chisq = c(NA_real_, NA_real_, .b5_ln[["chisq"]]), LRT_df = c(NA_real_, NA_real_, .b5_ln[["df"]]),
  LRT_p = c(NA_real_, NA_real_, .b5_ln[["p"]]), comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE)
B5_Smooths_tp <- .drop_re_smooths(tidy_smooths(m_B5_nonlin_re, "B5_nonlin_tp"))
B5_Smooths_cr <- .drop_re_smooths(tidy_smooths(m_B5_cr_re,     "B5_nonlin_cr"))
B5_ModelFit   <- dplyr::bind_rows(
  model_fit_row(m_B5_null_re, "B5_baseline_RE"), model_fit_row(m_B5_lin_re, "B5_lin_RE"),
  model_fit_row(m_B5_nonlin_re, "B5_nonlin_RE_tp"), model_fit_row(m_B5_cr_re, "B5_nonlin_RE_cr"))
.b5_slope <- function(mod, dat, gender) {
  bym <- attr(dat, "by_mean"); if (is.null(bym)) bym <- mean(dat$birth_year, na.rm = TRUE)
  gr  <- seq(min(dat$BYc, na.rm = TRUE), max(dat$BYc, na.rm = TRUE), length.out = 200)
  is_m <- gender == "male"
  nd0 <- data.frame(BYc = gr, east_west_c = 0, gender_c = if (is_m) -0.5 else 0.5,
                    parental_edu_z_kernel = 0, ParEdu_F = 0, ParEdu_M = 0, PGI_F = 0, PGI_M = 0)
  if ("fid_re" %in% names(mod$model)) { fl <- levels(mod$model$fid_re); nd0$fid_re <- factor(fl[1], levels = fl) }
  nd1 <- nd0; if (is_m) nd1$PGI_M <- 1 else nd1$PGI_F <- 1
  D <- predict(mod, nd1, type = "lpmatrix") - predict(mod, nd0, type = "lpmatrix")
  b <- coef(mod); V <- vcov(mod); fit <- as.vector(D %*% b); se <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))
  data.frame(birth_year = gr + bym, slope = fit, se = se, ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se, gender = gender)
}
B5_curves <- dplyr::bind_rows(.b5_slope(m_B5_nonlin_re, dat_cluster_mob, "female"),
                              .b5_slope(m_B5_nonlin_re, dat_cluster_mob, "male"))

# ---- B5 leave-one-study-out: nonlinearity + difference smooth ----
# B5_LOO: per fold, the omnibus linear-vs-nonlinear LRT plus female/male
# curvature LRTs, with the supported birth-year range — checks that the
# nonlinear cohort-varying gender pattern (male-driven) is not carried by a
# single study. Additional sensitivity (../../../Plan_deviations.md §14).
B5_LOO <- loco_pergender_gam(
  dat_cluster_mob, rhs_without_re(m_B5_lin_re), rhs_without_re(m_B5_nonlin_re),
  cells = list(female = "PGI_F", male = "PGI_M"), label = "B5")
# B5_LOO_curves: per-fold region-pooled per-gender PGI-slope curves across birth
# year (leave-one-study-out analogue of B5_curves / Fig-4C). The curve is read
# from each fold's nonlinear per-gender fit, matching B3_LOO_curves /
# B6_LOO_curves.
B5_LOO_curves <- loco_pergender_curves(
  dat_cluster_mob, rhs_without_re(m_B5_nonlin_re), .b5_slope)
# B5_LOO_diffsmooth / B5_LOO_diff_intervals: per fold, the female-minus-male
# PGI-slope difference smooth over birth year with BOTH pointwise and
# simultaneous 95% bands (same orientation and Monte-Carlo procedure as the
# primary B5 export, diff_smooth_pooled_gender_simul), restricted to the fold's
# supported birth-year range, and the contiguous birth-year intervals where the
# simultaneous band excludes zero.
.b5_loco_ds        <- loco_b5_diffsmooth(dat_cluster_mob, rhs_without_re(m_B5_nonlin_re))
B5_LOO_diffsmooth      <- .b5_loco_ds$smooth
B5_LOO_diff_intervals  <- .b5_loco_ds$intervals

# ---- B6: Region x Gender PGI-mobility GAM (RE), cell-specific ParEdu ----
cat("B6 Region x Gender GAM (RE)...\n")
.b6_rhs <- function(level, k, basis) {
  pe  <- paste("parental_edu_z_kernel + parental_edu_z_kernel:east_west_c + parental_edu_z_kernel:gender_c +",
               .sm("ParEdu_WF", basis), "+", .sm("ParEdu_WM", basis), "+", .sm("ParEdu_EF", basis), "+", .sm("ParEdu_EM", basis))
  rhs <- paste("mobility ~ BYc + east_west_c + gender_c +", pe, "+ PGI_WF + PGI_WM + PGI_EF + PGI_EM")
  if (level %in% c("lin", "nonlin")) rhs <- paste(rhs, "+ PGI_WF:BYc + PGI_WM:BYc + PGI_EF:BYc + PGI_EM:BYc")
  if (level == "nonlin") rhs <- paste(rhs, "+", .sm("PGI_WF", basis), "+", .sm("PGI_WM", basis), "+", .sm("PGI_EF", basis), "+", .sm("PGI_EM", basis))
  rhs
}
m_B6_null_re   <- mgcv::bam(as.formula(paste(.b6_rhs("null",   K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B6_lin_re    <- mgcv::bam(as.formula(paste(.b6_rhs("lin",    K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B6_nonlin_re <- mgcv::bam(as.formula(paste(.b6_rhs("nonlin", K_DEFAULT, BS_DEFAULT), re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
m_B6_cr_re     <- mgcv::bam(as.formula(paste(.b6_rhs("nonlin", K_SENS,    BS_SENS),    re_term)), data = dat_cluster_mob, method = "fREML", discrete = TRUE)
assert_pgi_in_baseline(m_B6_null_re, m_B6_lin_re, c("PGI_WF", "PGI_WM", "PGI_EF", "PGI_EM"), "B6")
.b6_ln <- lrt_pair(m_B6_lin_re, m_B6_nonlin_re)   # linear -> nonlinear (fREML)
.b6_0 <- .gf(m_B6_null_re); .b6_1 <- .gf(m_B6_lin_re); .b6_2 <- .gf(m_B6_nonlin_re)
B6_NestedLRT <- data.frame(model = c("B6_baseline", "B6_lin", "B6_nonlin"),
  AIC = c(.b6_0["AIC"], .b6_1["AIC"], .b6_2["AIC"]), BIC = c(.b6_0["BIC"], .b6_1["BIC"], .b6_2["BIC"]),
  logLik = c(.b6_0["logLik"], .b6_1["logLik"], .b6_2["logLik"]), dev_expl = c(.b6_0["dev"], .b6_1["dev"], .b6_2["dev"]),
  LRT_chisq = c(NA_real_, NA_real_, .b6_ln[["chisq"]]), LRT_df = c(NA_real_, NA_real_, .b6_ln[["df"]]),
  LRT_p = c(NA_real_, NA_real_, .b6_ln[["p"]]), comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE)
B6_Smooths_tp  <- .drop_re_smooths(tidy_smooths(m_B6_nonlin_re, "B6_nonlin_tp"))
B6_Smooths_pgi <- B6_Smooths_tp[grepl("PGI", B6_Smooths_tp$smooth_term), , drop = FALSE]
B6_Smooths_cr  <- .drop_re_smooths(tidy_smooths(m_B6_cr_re, "B6_nonlin_cr"))
B6_ModelFit    <- dplyr::bind_rows(
  model_fit_row(m_B6_null_re, "B6_baseline_RE"), model_fit_row(m_B6_lin_re, "B6_lin_RE"),
  model_fit_row(m_B6_nonlin_re, "B6_nonlin_RE_tp"), model_fit_row(m_B6_cr_re, "B6_nonlin_RE_cr"))
B6_OverallSlopes <- extract_region_gender_pgi_effects(m_B6_null_re, dat_cluster_mob)
.b6_kc <- tryCatch(mgcv::k.check(m_B6_nonlin_re), error = function(e) NULL)
B6_K_Index <- if (!is.null(.b6_kc)) {
  d <- data.frame(model = "B6_nonlin_RE_tp", smooth_term = rownames(.b6_kc), k_prime = .b6_kc[, "k'"],
    edf = .b6_kc[, "edf"], k_index = .b6_kc[, "k-index"], p_value = .b6_kc[, "p-value"], stringsAsFactors = FALSE)
  d[!grepl("fid_re", d$smooth_term), , drop = FALSE]
} else NULL
cat("B6 LOO...\n")
.b6_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster_mob))
.b6_rows <- list(); .b6_curves <- list()
for (fl in names(.b6_folds)) {
  dl  <- .b6_folds[[fl]]
  sub <- if (length(dl) == 0) dat_cluster_mob else dat_cluster_mob[!(as.character(dat_cluster_mob$cohort) %in% dl), , drop = FALSE]
  attr(sub, "by_mean") <- attr(dat_cluster_mob, "by_mean")
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  fold_lab <- if (length(dl) == 0) "Full sample" else paste0("-", fl); held <- if (length(dl) == 0) "(none)" else fl
  m_lin <- mgcv::bam(as.formula(paste(.b6_rhs("lin",    K_DEFAULT, BS_DEFAULT), rt)), data = sub, method = "fREML", discrete = TRUE)
  m_nl  <- mgcv::bam(as.formula(paste(.b6_rhs("nonlin", K_DEFAULT, BS_DEFAULT), rt)), data = sub, method = "fREML", discrete = TRUE)
  v <- lrt_pair(m_lin, m_nl)
  concl <- if (is.na(v[["p"]])) "Not tested" else if (v[["p"]] < 0.05) "Nonlinearity supported" else "No evidence that nonlinear specification is needed"
  .b6_rows[[fl]] <- data.frame(fold = fold_lab, held_out = held, N = nrow(sub),
    LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1), LRT_p = signif(v[["p"]], 4),
    nonlinearity_conclusion = concl, status = "OK", stringsAsFactors = FALSE)
  cv <- dplyr::bind_rows(
    extract_pgi_slope(m_nl, sub, east_west = "West", gender = "female"),
    extract_pgi_slope(m_nl, sub, east_west = "West", gender = "male"),
    extract_pgi_slope(m_nl, sub, east_west = "East", gender = "female"),
    extract_pgi_slope(m_nl, sub, east_west = "East", gender = "male"))
  cv$fold <- fold_lab; .b6_curves[[fl]] <- cv
}
B6_LOO <- dplyr::bind_rows(.b6_rows)
B6_LOO_curves <- dplyr::bind_rows(.b6_curves)
B6_LOO_curves$fold  <- factor(B6_LOO_curves$fold, levels = unique(B6_LOO_curves$fold))
B6_LOO_curves$group <- paste(B6_LOO_curves$east_west, B6_LOO_curves$gender)

# ---- B0g: overall Region x Gender PGI-mobility contrast (RE), full + pre ----
cat("B0g overall Region x Gender contrast (RE; full + pre-reunif)...\n")
B0g_focal <- "PGI_Edu_z:east_west_c:gender_c"
.b0g_cell_slopes <- function(mod, dat) {
  cells <- list(c("West", "female", -0.5, 0.5), c("West", "male", -0.5, -0.5),
                c("East", "female", 0.5, 0.5), c("East", "male", 0.5, -0.5))
  fl <- if ("fid_re" %in% names(mod$model)) levels(mod$model$fid_re) else NULL
  b <- coef(mod); V <- vcov(mod)
  do.call(rbind, lapply(cells, function(cl) {
    rc <- as.numeric(cl[3]); gc <- as.numeric(cl[4])
    mk <- function(pg) {
      nd <- data.frame(PGI_Edu_z = pg, east_west_c = rc, gender_c = gc, parental_edu_z_kernel = 0, BYc_z = 0)
      if (!is.null(fl)) nd$fid_re <- factor(fl[1], levels = fl)
      predict(mod, newdata = nd, type = "lpmatrix")
    }
    D <- mk(1) - mk(0); e <- as.numeric(D %*% b); s <- sqrt(as.numeric(D %*% V %*% t(D)))
    data.frame(east_west = cl[1], gender = cl[2], group = paste(cl[1], cl[2]),
      slope = e, se = s, ci_lo = e - 1.96 * s, ci_hi = e + 1.96 * s, p = 2 * pnorm(-abs(e / s)), stringsAsFactors = FALSE)
  }))
}
.fit_b0g <- function(dat, tag) {
  m <- mgcv::bam(as.formula(paste(
    "mobility ~ PGI_Edu_z * east_west_c * gender_c + parental_edu_z_kernel * east_west_c * gender_c + BYc_z", re_term)),
    data = dat, method = "fREML", discrete = TRUE)
  list(coefs = tidy_coefs(m, tag), fit = model_fit_row(m, tag), cell_slopes = .b0g_cell_slopes(m, dat), N = nrow(dat))
}
B0g <- .fit_b0g(dat_cluster_mob, "B0g_Mob")
.dcm_pre <- dat_cluster_mob[dat_cluster_mob$reunif == "pre", , drop = FALSE]
.dcm_pre$fid_re <- droplevels(factor(.dcm_pre$fid_re))
B0g_pre <- if (nrow(.dcm_pre) >= 100 && length(unique(.dcm_pre$east_west)) == 2 && nlevels(.dcm_pre$fid_re) >= 2)
  .fit_b0g(.dcm_pre, "B0g_Mob_pre") else NULL

# ---- B0g LOO: leave-one-study-out of the Region x Gender per-cell slopes ----
# Leave-one-study-out analogue of the Fig-4A per-cell region x gender slopes
# (source: B0g$cell_slopes). Per fold, refit the constant-slope B0g model on the
# reduced mobility sample and re-extract the four Region x Gender cell slopes.
# Fold that drops TwinLife uses the no-family-RE term (.loo_re_term). Four rows
# per fold (Full sample first).
.b0g_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster_mob))
.b0g_loo   <- list()
for (fl in names(.b0g_folds)) {
  dl   <- .b0g_folds[[fl]]; full <- length(dl) == 0
  sub  <- if (full) dat_cluster_mob else
    dat_cluster_mob[!(as.character(dat_cluster_mob$cohort) %in% dl), , drop = FALSE]
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(as.formula(paste(
    "mobility ~ PGI_Edu_z * east_west_c * gender_c + parental_edu_z_kernel * east_west_c * gender_c + BYc_z", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  cs <- .b0g_cell_slopes(m, sub)
  cs$fold     <- if (full) "Full sample" else paste0("-", fl)
  cs$held_out <- if (full) "(none)" else fl
  cs$N        <- nrow(sub)
  .b0g_loo[[fl]] <- cs
}
B0g_LOO_cells <- dplyr::bind_rows(.b0g_loo)
B0g_LOO_cells$fold <- factor(B0g_LOO_cells$fold, levels = unique(B0g_LOO_cells$fold))

# ---- Per-cell nonlinearity: clean nested LRTs + display curves (B3, B5, B6) ----
cat("Per-cell nonlinearity LRTs + display curves...\n")
.cell_lrt_table <- function(base_rhs, lin_mod, cell_smooths, dat, rt, tag) {
  do.call(rbind, lapply(names(cell_smooths), function(cell) {
    m <- mgcv::bam(as.formula(paste(base_rhs, "+", cell_smooths[[cell]], rt)), data = dat, method = "fREML", discrete = TRUE)
    v <- lrt_pair(lin_mod, m)
    data.frame(analysis = tag, cell = cell, LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1),
      LRT_p = signif(v[["p"]], 4), nonlinear_sig = isTRUE(!is.na(v[["p"]]) && v[["p"]] < 0.05), stringsAsFactors = FALSE)
  }))
}
B3_cell_LRT <- .cell_lrt_table(b3_formula("lin", BS_DEFAULT, ""), m_B3_lin_re,
  list(West = .sm("PGI_W"), East = .sm("PGI_E")), dat_cluster_mob, re_term, "B3")
B3_display_curves <- dplyr::bind_rows(lapply(c("West", "East"), function(rg) {
  sig <- B3_cell_LRT$nonlinear_sig[B3_cell_LRT$cell == rg]
  m <- if (isTRUE(sig)) m_B3_nonlin_re else m_B3_lin_re
  cv <- extract_pgi_slope(m, dat_cluster_mob, east_west = rg); cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv }))
B5_cell_LRT <- .cell_lrt_table(.b5_rhs("lin", K_DEFAULT, BS_DEFAULT), m_B5_lin_re,
  list(female = .sm("PGI_F"), male = .sm("PGI_M")), dat_cluster_mob, re_term, "B5")
B5_curves <- dplyr::bind_rows(lapply(c("female", "male"), function(g) {
  sig <- B5_cell_LRT$nonlinear_sig[B5_cell_LRT$cell == g]
  m <- if (isTRUE(sig)) m_B5_nonlin_re else m_B5_lin_re
  cv <- .b5_slope(m, dat_cluster_mob, g); cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv }))
B6_cell_LRT <- .cell_lrt_table(.b6_rhs("lin", K_DEFAULT, BS_DEFAULT), m_B6_lin_re,
  list(`West female` = .sm("PGI_WF"), `West male` = .sm("PGI_WM"), `East female` = .sm("PGI_EF"), `East male` = .sm("PGI_EM")),
  dat_cluster_mob, re_term, "B6")
B6_display_curves <- dplyr::bind_rows(lapply(
  list(c("West", "female"), c("West", "male"), c("East", "female"), c("East", "male")), function(cl) {
    sig <- B6_cell_LRT$nonlinear_sig[B6_cell_LRT$cell == paste(cl[1], cl[2])]
    m <- if (isTRUE(sig)) m_B6_nonlin_re else m_B6_lin_re
    cv <- extract_pgi_slope(m, dat_cluster_mob, east_west = cl[1], gender = cl[2])
    cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv }))
cat("Gender block done.\n")

# ---- Save ----
model_cache <- list(
  RUN_TS = RUN_TS_VAL, by_mean = attr(dat_cluster_mob, "by_mean"),
  dat_cluster_mob = dat_cluster_mob, dat_cluster_mob_pre = dat_cluster_mob_pre,
  m_B1_re = m_B1_re, m_B2_re = m_B2_re,
  m_B3_null_re = m_B3_null_re, m_B3_lin_re = m_B3_lin_re,
  m_B3_nonlin_re = m_B3_nonlin_re, m_B3_cr_re = m_B3_cr_re,
  m_B0_re = m_B0_re, m_B0_pre_re = m_B0_pre_re,
  B1_Coefficients = B1_Coefficients, B1_ModelFit = B1_ModelFit,
  B2_Coefficients = B2_Coefficients, B2_ModelFit = B2_ModelFit,
  B3_NestedLRT = B3_NestedLRT, B3_Smooths_tp = B3_Smooths_tp,
  B3_Smooths_cr = B3_Smooths_cr, B3_Smooths_pgi_tp = B3_Smooths_pgi_tp,
  B3_K_Index = B3_K_Index, B_ModelFit = B_ModelFit,
  B0_RegionSlopes = B0_RegionSlopes, B0_Coefficients = B0_Coefficients,
  B0_pre_RegionSlopes = B0_pre_RegionSlopes, B0_pre_Coefficients = B0_pre_Coefficients,
  B1_LOO = B1_LOO, B2_LOO = B2_LOO, B3_LOO = B3_LOO, B3_LOO_curves = B3_LOO_curves,
  B0_LOO = B0_LOO, B0_pre_LOO = B0_pre_LOO, B0_LOO_slopes = B0_LOO_slopes,
  desc_cohort = desc_cohort, desc_cohort_region = desc_cohort_region,
  pgi_correlations = pgi_correlations,
  # ---- Gender block (Section C mobility analyses + per-gender GAM, family RE) ----
  B4_focal = B4_focal, m_B4_re = m_B4_re, B4_Coefficients = B4_Coefficients,
  B4_ModelFit = B4_ModelFit, B4_LOO = B4_LOO, B4_LOO_gender = B4_LOO_gender,
  m_B5_lin_re = m_B5_lin_re, m_B5_nonlin_re = m_B5_nonlin_re, m_B5_cr_re = m_B5_cr_re,
  B5_NestedLRT = B5_NestedLRT, B5_Smooths_tp = B5_Smooths_tp, B5_Smooths_cr = B5_Smooths_cr,
  B5_ModelFit = B5_ModelFit, B5_curves = B5_curves, B5_cell_LRT = B5_cell_LRT,
  B5_LOO = B5_LOO, B5_LOO_curves = B5_LOO_curves, B5_LOO_diffsmooth = B5_LOO_diffsmooth,
  B5_LOO_diff_intervals = B5_LOO_diff_intervals,
  B3_cell_LRT = B3_cell_LRT, B3_display_curves = B3_display_curves,
  m_B6_lin_re = m_B6_lin_re, m_B6_nonlin_re = m_B6_nonlin_re, m_B6_cr_re = m_B6_cr_re,
  B6_NestedLRT = B6_NestedLRT, B6_Smooths_tp = B6_Smooths_tp, B6_Smooths_pgi = B6_Smooths_pgi,
  B6_Smooths_cr = B6_Smooths_cr, B6_ModelFit = B6_ModelFit, B6_K_Index = B6_K_Index,
  B6_OverallSlopes = B6_OverallSlopes, B6_LOO = B6_LOO, B6_LOO_curves = B6_LOO_curves,
  B6_cell_LRT = B6_cell_LRT, B6_display_curves = B6_display_curves,
  B0g_focal = B0g_focal, B0g = B0g, B0g_pre = B0g_pre, B0g_LOO_cells = B0g_LOO_cells
)
OUT_RDS <- file.path(OUT_DIR_B, "B_Mobility_Models.rds")
saveRDS(model_cache, OUT_RDS)
cat(sprintf("B_analysis: models RDS saved to %s\n", OUT_RDS))
