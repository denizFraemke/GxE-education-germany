# A_analysis.R
# Section A entry point. Loads data, builds samples, computes descriptives,
# fits A1/A2/A3 and the gender models A4/A5/A6/A0g/A-spec (family RE on
# TwinLife families), saves the models RDS.
#
# Usage (from run_A.command):
#   Rscript 01_PGI_Analysis/A_Education/A_analysis.R
#
# Inputs (env vars; defaults via run_context.R):
#   DATA_ROOT (default: /path/to/data_root)
#   RUN_TS    (default: format(Sys.time(), "%Y%m%d_%H%M"))
#
# Outputs (under ${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/A_Education/):
#   A_Education_Models.rds  — fitted models + analytic frames + key tables
#
# A_report.R then reads this RDS to build
# Supplementary_Data_2_Education_results.xlsx, A_Education_Plots.pdf,
# A_Education_LOO_Plots.pdf and A_Education_Overview.md.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
})

# ---- Resolve paths (works under both Rscript and source()) ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  # Rscript replaces spaces in the --file= argument with `~+~`; undo that
  # before normalizePath so paths like ".../015_Data analysis/..." resolve.
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine A_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR  # flattened: 01_PGI_Analysis is the workflow + script root
rm(.find_this_file)

# ---- Run context: RUN_TS and output dirs ----
source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL  <- get_run_ts()
OUT_DIR_A   <- output_dir("A_Education")

# ---- Constants and helpers ----
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"),
       local = TRUE)

# ---- Result collectors required by the 00_setup/ scripts ----
RESULTS         <- list()
PLOTS           <- list()
PRIMARY_RESULTS <- list()

# ---- Load + prepare samples ----
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var
DATA_DIR        <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")

source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

cat(sprintf("A_analysis: dat_cluster N = %d\n", nrow(dat_cluster)))

# ---- Descriptives (no gender split) ----
source(file.path(SECTION_DIR, "A_descriptives.R"), local = TRUE)

# All RE fits use mgcv::bam() rather than gam(): bam is built for large
# samples with big random effects (here s(fid_re, bs="re") has ~1,323
# levels). Fits use method="fREML" + discrete=TRUE (fast REML with covariate
# discretization). bam objects inherit class "gam", so all downstream helpers
# apply unchanged.

# ---- A1: linear three-way interaction (RE) ----
cat("A1 (RE)...\n")
re_term <- "+ s(fid_re, bs = 're')"
m_A1_re <- mgcv::bam(
  as.formula(paste(
    "edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3", re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

# ---- A2: reunification step-function (RE) ----
cat("A2 (RE)...\n")
m_A2_re <- mgcv::bam(
  as.formula(paste(
    "edu_z_kernel ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3", re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

# ---- A3: GAM varying-coefficient (RE) — baseline / linear / nonlinear ----
# Three nested formulas so build_gam_overview_table can compute the nested
# LRT ladder on the RE sample.
cat("A3 baseline (RE)...\n")
# Baseline carries the PGI main-effect cells (PGI_W + PGI_E, constant slope),
# so the baseline -> linear step isolates the region-specific cohort slopes
# (PGI_W:BYc + PGI_E:BYc, 2 df) and is not confounded with the PGI main
# effect. Enforced by assert_pgi_in_baseline() below.
m_A3_null_re <- mgcv::bam(
  as.formula(paste("edu_z_kernel ~ BYc + east_west_c + gender_c + PGI_W + PGI_E", re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

cat("A3 linear (RE)...\n")
m_A3_lin_re <- mgcv::bam(
  as.formula(paste(
    "edu_z_kernel ~ BYc + east_west_c + gender_c +
       PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc", re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

cat("A3 nonlinear tp (RE)...\n")
m_A3_nonlin_re <- mgcv::bam(
  as.formula(sprintf(
    "edu_z_kernel ~ BYc + east_west_c + gender_c +
       PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc +
       s(BYc, by = PGI_W, k = %d, bs = '%s') +
       s(BYc, by = PGI_E, k = %d, bs = '%s') %s",
    K_DEFAULT, BS_DEFAULT, K_DEFAULT, BS_DEFAULT, re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

# Sensitivity: cubic regression spline basis (cr).
cat("A3 nonlinear cr sensitivity (RE)...\n")
m_A3_cr_re <- mgcv::bam(
  as.formula(sprintf(
    "edu_z_kernel ~ BYc + east_west_c + gender_c +
       PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc +
       s(BYc, by = PGI_W, k = %d, bs = '%s') +
       s(BYc, by = PGI_E, k = %d, bs = '%s') %s",
    K_SENS, BS_SENS, K_SENS, BS_SENS, re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)

# Guard: the PGI main effect lives in the baseline, so the baseline is the
# plan's "model without cohort-varying terms" (constant PGI slope) and the
# AIC/BIC ladder below is an honest information-criteria comparison. A3's
# focal test is non-linearity only (linear -> nonlinear); the linear PGI
# cohort trend is tested in A1 (Wald), so there is no baseline -> linear LRT.
assert_pgi_in_baseline(m_A3_null_re, m_A3_lin_re, c("PGI_W", "PGI_E"), "A3")

# ---- Coefficient + smooth tables ----
# Helper: drop random-effects smooths (s(fid_re), bs="re") from a
# tidy_smooths() table so the GAM nonlinearity reporting reflects only
# the PGI x birth-year smooths, not the family random intercept.
.drop_re_smooths <- function(sm) {
  if (is.null(sm) || !nrow(sm)) return(sm)
  sm[!grepl("fid_re", sm$smooth_term), , drop = FALSE]
}

A1_Coefficients <- tidy_coefs(m_A1_re,        "A1_RE")
A1_ModelFit     <- model_fit_row(m_A1_re,     "A1_RE")
A2_Coefficients <- tidy_coefs(m_A2_re,        "A2_RE")
A2_ModelFit     <- model_fit_row(m_A2_re,     "A2_RE")
# A3 focal test = non-linearity: linear -> nonlinear (shared fixed effects,
# fREML LRT). The baseline/linear/nonlinear rows carry AIC/BIC/logLik/dev_expl
# as the plan's information-criteria ladder (level a); only the nonlinear row
# carries an LRT. The linear cohort trend is A1's (Wald) - no baseline->linear.
.a3_lrt_ln <- lrt_pair(m_A3_lin_re,  m_A3_nonlin_re)     # linear -> nonlinear (fREML)
.a3_gf <- function(mod) c(
  AIC = AIC(mod), BIC = BIC(mod),
  logLik = as.numeric(logLik(mod)), dev = summary(mod)$dev.expl)
.f0 <- .a3_gf(m_A3_null_re); .f1 <- .a3_gf(m_A3_lin_re); .f2 <- .a3_gf(m_A3_nonlin_re)
A3_NestedLRT <- data.frame(
  model      = c("A3_baseline", "A3_lin", "A3_nonlin"),
  AIC        = c(.f0["AIC"],    .f1["AIC"],    .f2["AIC"]),
  BIC        = c(.f0["BIC"],    .f1["BIC"],    .f2["BIC"]),
  logLik     = c(.f0["logLik"], .f1["logLik"], .f2["logLik"]),
  dev_expl   = c(.f0["dev"],    .f1["dev"],    .f2["dev"]),
  LRT_chisq  = c(NA_real_, NA_real_, .a3_lrt_ln[["chisq"]]),
  LRT_df     = c(NA_real_, NA_real_, .a3_lrt_ln[["df"]]),
  LRT_p      = c(NA_real_, NA_real_, .a3_lrt_ln[["p"]]),
  comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE
)
# Smooth tables from the REML fits, with the RE smooth filtered out.
A3_Smooths_tp   <- .drop_re_smooths(tidy_smooths(m_A3_nonlin_re, "A3_RE_nonlin_tp"))
A3_Smooths_cr   <- .drop_re_smooths(tidy_smooths(m_A3_cr_re,     "A3_RE_nonlin_cr"))
A_ModelFit      <- dplyr::bind_rows(
  A1_ModelFit, A2_ModelFit,
  model_fit_row(m_A3_null_re,   "A3_baseline_RE"),
  model_fit_row(m_A3_lin_re,    "A3_lin_RE"),
  model_fit_row(m_A3_nonlin_re, "A3_nonlin_RE_tp"),
  model_fit_row(m_A3_cr_re,     "A3_nonlin_RE_cr")
)

# k.check diagnostics for the GAM smooths.
k_index_rows <- list()
for (lbl in c("A3_nonlin_RE_tp", "A3_nonlin_RE_cr")) {
  mod <- if (lbl == "A3_nonlin_RE_tp") m_A3_nonlin_re else m_A3_cr_re
  kc <- tryCatch(mgcv::k.check(mod), error = function(e) NULL)
  if (!is.null(kc)) {
    df_kc <- data.frame(
      model       = lbl,
      smooth_term = rownames(kc),
      k_prime     = kc[, "k'"],
      edf         = kc[, "edf"],
      k_index     = kc[, "k-index"],
      p_value     = kc[, "p-value"],
      stringsAsFactors = FALSE
    )
    # Drop the random-effects smooth row (not a nonlinearity diagnostic).
    df_kc <- df_kc[!grepl("fid_re", df_kc$smooth_term), , drop = FALSE]
    k_index_rows[[lbl]] <- df_kc
  }
}
A3_K_Index <- if (length(k_index_rows) > 0) dplyr::bind_rows(k_index_rows) else NULL

# ---- Leave-one-study-out (LOO) sensitivity ----
# Refit each focal model on the analytic sample minus one study fold
# (BASE-II / SHIP / SOEP / TwinLife; SHIP is one fold) and track the
# focal result. Uses the same bam fitting as the main models. Only the
# small result tables are cached, not the refit models.
cat("LOO: leave-one-study-out refits (A1, A2, A3)...\n")

# Per-fold random-effects term. When a fold drops TwinLife, no families
# remain and s(fid_re, bs="re") would be a single-level (degenerate) RE
# that mgcv cannot fit — so that fold is fit WITHOUT the family RE (there
# is no within-family clustering to model). droplevels() also clears the
# stale TL_<fid> levels left by subsetting.
.loo_re_term <- function(sub) {
  if (!"fid_re" %in% names(sub)) return("")
  sub$fid_re <- droplevels(factor(sub$fid_re))
  if (nlevels(sub$fid_re) >= 2) "+ s(fid_re, bs = 're')" else ""
}

# A1 — linear focal coefficient PGI_Edu_z:BYc_z:east_west_c
loo_fit_A1 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re))
  rt <- .loo_re_term(sub)
  m <- mgcv::bam(
    as.formula(paste(
      "edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "A1_RE")
  r <- r[r$term == "PGI_Edu_z:BYc_z:east_west_c", ]
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5),
    p = round(r$p.value, 5))
}
A1_LOO <- build_loo_table(
  leave_one_cohort_out(dat_cluster, loo_fit_A1), estimate_col = "estimate")

# A2 — step-function focal coefficient PGI_Edu_z:reunifpost:east_west_c
loo_fit_A2 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re))
  rt <- .loo_re_term(sub)
  m <- mgcv::bam(
    as.formula(paste(
      "edu_z_kernel ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "A2_RE")
  r <- r[r$term == "PGI_Edu_z:reunifpost:east_west_c", ]
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5),
    p = round(r$p.value, 5))
}
A2_LOO <- build_loo_table(
  leave_one_cohort_out(dat_cluster, loo_fit_A2), estimate_col = "estimate")

# A3 — nonlinearity via the linear -> nonlinear LRT, plus per-fold
# cohort-varying slope curves for the LOO PDF. One pass over the folds
# fits the linear + nonlinear GAM and produces both the LRT row and the
# West/East slope curves, so no model is fit twice. Both fits share the
# same fixed effects, so fREML log-likelihoods are comparable (matching
# the main A3 linear->nonlinear test). The fold that drops TwinLife is
# fit without the family RE (see .loo_re_term).
a3_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster))
a3_loo_rows   <- list()
a3_loo_curves <- list()
for (fl in names(a3_folds)) {
  drop_levs <- a3_folds[[fl]]
  sub <- if (length(drop_levs) == 0) dat_cluster else
    dat_cluster[!(as.character(dat_cluster$cohort) %in% drop_levs), , drop = FALSE]
  attr(sub, "by_mean") <- attr(dat_cluster, "by_mean")
  sub$fid_re <- droplevels(factor(sub$fid_re))
  rt       <- .loo_re_term(sub)
  fold_lab <- if (length(drop_levs) == 0) "Full sample" else paste0("-", fl)
  held     <- if (length(drop_levs) == 0) "(none)" else fl

  m_lin <- mgcv::bam(
    as.formula(paste(
      "edu_z_kernel ~ BYc + east_west_c + gender_c +
         PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  m_nl <- mgcv::bam(
    as.formula(sprintf(
      "edu_z_kernel ~ BYc + east_west_c + gender_c +
         PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc +
         s(BYc, by = PGI_W, k = %d, bs = '%s') +
         s(BYc, by = PGI_E, k = %d, bs = '%s') %s",
      K_DEFAULT, BS_DEFAULT, K_DEFAULT, BS_DEFAULT, rt)),
    data = sub, method = "fREML", discrete = TRUE)

  v <- lrt_pair(m_lin, m_nl)
  concl <- if (is.na(v[["p"]])) "Not tested"
           else if (v[["p"]] < 0.05) "Nonlinearity supported"
           else "No evidence that nonlinear specification is needed"
  a3_loo_rows[[fl]] <- data.frame(
    fold = fold_lab, held_out = held, N = nrow(sub),
    LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1),
    LRT_p = signif(v[["p"]], 4), nonlinearity_conclusion = concl,
    status = "OK", stringsAsFactors = FALSE)

  cv <- dplyr::bind_rows(
    extract_pgi_slope(m_nl, sub, east_west = "West"),
    extract_pgi_slope(m_nl, sub, east_west = "East"))
  cv$fold <- fold_lab
  a3_loo_curves[[fl]] <- cv
}
A3_LOO        <- dplyr::bind_rows(a3_loo_rows)
A3_LOO_curves <- dplyr::bind_rows(a3_loo_curves)
# Preserve fold order (Full sample first) for faceting.
A3_LOO_curves$fold <- factor(A3_LOO_curves$fold, levels = unique(A3_LOO_curves$fold))

cat("LOO done.\n")

# =====================================================================
# GENDER as a main moderator of attainment (Section A, family random effects).
# A4/A6/A0g/A-spec are the Section C analyses of the preregistration
# (C1/C3/C0/C2), estimated here with family random effects rather than the
# preregistered dedup-OLS (see ../../../Plan_deviations.md §10). A5 is the
# per-gender attainment GAM.
#   A4   linear four-way PGI x BY x Region x Gender         (C1 in the plan)
#   A5   per-gender PGI-over-birth-year GAM, regions pooled
#   A6   Region x Gender PGI cohort GAM                     (C3 in the plan)
#   A0g  overall Region x Gender PGI contrast (+ pre-reunif)(C0 in the plan)
#   A-spec  logistic-transition specification check         (C2 in the plan)
# A4/A5/A6/A0g use the full-family RE sample (dat_cluster). A-spec is an
# NLS specification check with no RE analogue, kept on the dedup sample
# (dat_edu), reported as a convergence/fit diagnostic only.
# =====================================================================

# ---- A4: linear four-way interaction (RE) ----
cat("A4 four-way (RE)...\n")
A4_focal <- "PGI_Edu_z:BYc_z:east_west_c:gender_c"
m_A4_re <- mgcv::bam(
  as.formula(paste(
    "edu_z_kernel ~ PGI_Edu_z * BYc_z * east_west_c * gender_c", re_term)),
  data = dat_cluster, method = "fREML", discrete = TRUE)
A4_Coefficients <- tidy_coefs(m_A4_re, "A4_RE")
A4_ModelFit     <- model_fit_row(m_A4_re, "A4_RE")
loo_fit_A4 <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(as.formula(paste(
    "edu_z_kernel ~ PGI_Edu_z * BYc_z * east_west_c * gender_c", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  r <- tidy_coefs(m, "A4_RE"); r <- r[r$term == A4_focal, ]
  if (!nrow(r)) return(c(estimate = NA_real_, SE = NA_real_, p = NA_real_))
  c(estimate = round(r$estimate, 5), SE = round(r$std.error, 5), p = round(r$p.value, 5))
}
A4_LOO <- build_loo_table(
  leave_one_cohort_out(dat_cluster, loo_fit_A4), estimate_col = "estimate")
# A4_LOO_gender: leave-one-study-out of the POOLED PGI x gender interaction
# (PGI_Edu_z:gender_c) from the same A4 model; A4_LOO tracks only the four-way
# term. Finding-matched to the reported PGI x gender attainment result;
# additional sensitivity, not preregistered (../../../Plan_deviations.md §14).
# Carries estimate, SE, 95% CI, exact p, supported birth-year range,
# delta-vs-full and a convergence flag.
A4_LOO_gender <- loco_pooled_coef(
  dat_cluster, rhs_without_re(m_A4_re), "PGI_Edu_z:gender_c", "A4_gender_RE")

# ---- A5: per-gender PGI-over-birth-year GAM, regions pooled (RE) ----
# Gender analog of A3: null/linear/nonlinear nested LRT plus the per-gender
# PGI->education slope over birth year (the headline Female-vs-Male curve).
# Regions enter only as a fixed-effect main term.
cat("A5 per-gender pooled GAM (RE)...\n")
dat_cluster$PGI_F <- dat_cluster$PGI_Edu_z * (1L - dat_cluster$gender_01)
dat_cluster$PGI_M <- dat_cluster$PGI_Edu_z *        dat_cluster$gender_01
.a5_rhs <- function(level, k, basis) {
  # Baseline carries the PGI main-effect cells (PGI_F + PGI_M) so
  # baseline -> linear isolates the per-gender cohort slopes (2 df).
  rhs <- "edu_z_kernel ~ BYc + east_west_c + gender_c + PGI_F + PGI_M"
  if (level %in% c("lin", "nonlin"))
    rhs <- paste(rhs, "+ PGI_F:BYc + PGI_M:BYc")
  if (level == "nonlin")
    rhs <- paste(rhs, sprintf(
      "+ s(BYc, by = PGI_F, k = %d, bs = '%s') + s(BYc, by = PGI_M, k = %d, bs = '%s')",
      k, basis, k, basis))
  rhs
}
m_A5_null_re   <- mgcv::bam(as.formula(paste(.a5_rhs("null",   K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A5_lin_re    <- mgcv::bam(as.formula(paste(.a5_rhs("lin",    K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A5_nonlin_re <- mgcv::bam(as.formula(paste(.a5_rhs("nonlin", K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A5_cr_re     <- mgcv::bam(as.formula(paste(.a5_rhs("nonlin", K_SENS,    BS_SENS),    re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
assert_pgi_in_baseline(m_A5_null_re, m_A5_lin_re, c("PGI_F", "PGI_M"), "A5")
.a5_ln <- lrt_pair(m_A5_lin_re,  m_A5_nonlin_re)   # linear -> nonlinear (fREML)
.a5_g0 <- .a3_gf(m_A5_null_re); .a5_g1 <- .a3_gf(m_A5_lin_re); .a5_g2 <- .a3_gf(m_A5_nonlin_re)
A5_NestedLRT <- data.frame(
  model = c("A5_baseline", "A5_lin", "A5_nonlin"),
  AIC = c(.a5_g0["AIC"], .a5_g1["AIC"], .a5_g2["AIC"]),
  BIC = c(.a5_g0["BIC"], .a5_g1["BIC"], .a5_g2["BIC"]),
  logLik = c(.a5_g0["logLik"], .a5_g1["logLik"], .a5_g2["logLik"]),
  dev_expl = c(.a5_g0["dev"], .a5_g1["dev"], .a5_g2["dev"]),
  LRT_chisq = c(NA_real_, NA_real_, .a5_ln[["chisq"]]),
  LRT_df = c(NA_real_, NA_real_, .a5_ln[["df"]]),
  LRT_p = c(NA_real_, NA_real_, .a5_ln[["p"]]),
  comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE)
A5_Smooths_tp <- .drop_re_smooths(tidy_smooths(m_A5_nonlin_re, "A5_nonlin_tp"))
A5_Smooths_cr <- .drop_re_smooths(tidy_smooths(m_A5_cr_re,     "A5_nonlin_cr"))
A5_ModelFit   <- dplyr::bind_rows(
  model_fit_row(m_A5_null_re,   "A5_baseline_RE"),
  model_fit_row(m_A5_lin_re,    "A5_lin_RE"),
  model_fit_row(m_A5_nonlin_re, "A5_nonlin_RE_tp"),
  model_fit_row(m_A5_cr_re,     "A5_nonlin_RE_cr"))
# Per-gender PGI->education slope over birth year (RE-aware; the PGI_F/PGI_M
# pooled-gender cells are not handled by extract_pgi_slope, so use a bespoke
# lpmatrix contrast with an fid_re placeholder that cancels in nd1 - nd0).
.a5_slope <- function(mod, dat, gender) {
  bym <- attr(dat, "by_mean"); if (is.null(bym)) bym <- mean(dat$birth_year, na.rm = TRUE)
  gr  <- seq(min(dat$BYc, na.rm = TRUE), max(dat$BYc, na.rm = TRUE), length.out = 200)
  is_m <- gender == "male"
  nd0 <- data.frame(BYc = gr, east_west_c = 0,
                    gender_c = if (is_m) -0.5 else 0.5, PGI_F = 0, PGI_M = 0)
  if ("fid_re" %in% names(mod$model)) {
    fl <- levels(mod$model$fid_re); nd0$fid_re <- factor(fl[1], levels = fl)
  }
  nd1 <- nd0; if (is_m) nd1$PGI_M <- 1 else nd1$PGI_F <- 1
  D <- predict(mod, nd1, type = "lpmatrix") - predict(mod, nd0, type = "lpmatrix")
  b <- coef(mod); V <- vcov(mod)
  fit <- as.vector(D %*% b); se <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))
  data.frame(birth_year = gr + bym, slope = fit, se = se,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se, gender = gender)
}
A5_curves <- dplyr::bind_rows(.a5_slope(m_A5_nonlin_re, dat_cluster, "female"),
                              .a5_slope(m_A5_nonlin_re, dat_cluster, "male"))

# ---- A6: Region x Gender PGI cohort GAM (RE) ----
# Baseline/linear/nonlinear nested LRT with the four Region x Gender PGI
# cells, plus per-cell cohort-varying slope curves.
cat("A6 Region x Gender GAM (RE)...\n")
.a6_rhs <- function(level, k, basis) {
  rhs <- "edu_z_kernel ~ BYc + east_west_c + gender_c + PGI_WF + PGI_WM + PGI_EF + PGI_EM"
  if (level %in% c("lin", "nonlin"))
    rhs <- paste(rhs, "+ PGI_WF:BYc + PGI_WM:BYc + PGI_EF:BYc + PGI_EM:BYc")
  if (level == "nonlin")
    rhs <- paste(rhs, sprintf(
      "+ s(BYc, by = PGI_WF, k = %d, bs = '%s') + s(BYc, by = PGI_WM, k = %d, bs = '%s') + s(BYc, by = PGI_EF, k = %d, bs = '%s') + s(BYc, by = PGI_EM, k = %d, bs = '%s')",
      k, basis, k, basis, k, basis, k, basis))
  rhs
}
m_A6_null_re   <- mgcv::bam(as.formula(paste(.a6_rhs("null",   K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A6_lin_re    <- mgcv::bam(as.formula(paste(.a6_rhs("lin",    K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A6_nonlin_re <- mgcv::bam(as.formula(paste(.a6_rhs("nonlin", K_DEFAULT, BS_DEFAULT), re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
m_A6_cr_re     <- mgcv::bam(as.formula(paste(.a6_rhs("nonlin", K_SENS,    BS_SENS),    re_term)),
                            data = dat_cluster, method = "fREML", discrete = TRUE)
assert_pgi_in_baseline(m_A6_null_re, m_A6_lin_re,
                       c("PGI_WF", "PGI_WM", "PGI_EF", "PGI_EM"), "A6")
.a6_ln <- lrt_pair(m_A6_lin_re,  m_A6_nonlin_re)   # linear -> nonlinear (fREML)
.a6_g0 <- .a3_gf(m_A6_null_re); .a6_g1 <- .a3_gf(m_A6_lin_re); .a6_g2 <- .a3_gf(m_A6_nonlin_re)
A6_NestedLRT <- data.frame(
  model = c("A6_baseline", "A6_lin", "A6_nonlin"),
  AIC = c(.a6_g0["AIC"], .a6_g1["AIC"], .a6_g2["AIC"]),
  BIC = c(.a6_g0["BIC"], .a6_g1["BIC"], .a6_g2["BIC"]),
  logLik = c(.a6_g0["logLik"], .a6_g1["logLik"], .a6_g2["logLik"]),
  dev_expl = c(.a6_g0["dev"], .a6_g1["dev"], .a6_g2["dev"]),
  LRT_chisq = c(NA_real_, NA_real_, .a6_ln[["chisq"]]),
  LRT_df = c(NA_real_, NA_real_, .a6_ln[["df"]]),
  LRT_p = c(NA_real_, NA_real_, .a6_ln[["p"]]),
  comparison = c(NA_character_, NA_character_, "linear vs nonlinear"),
  stringsAsFactors = FALSE)
A6_Smooths_tp  <- .drop_re_smooths(tidy_smooths(m_A6_nonlin_re, "A6_nonlin_tp"))
A6_Smooths_pgi <- A6_Smooths_tp[grepl("PGI", A6_Smooths_tp$smooth_term), , drop = FALSE]
A6_Smooths_cr  <- .drop_re_smooths(tidy_smooths(m_A6_cr_re, "A6_nonlin_cr"))
A6_ModelFit    <- dplyr::bind_rows(
  model_fit_row(m_A6_null_re,   "A6_baseline_RE"),
  model_fit_row(m_A6_lin_re,    "A6_lin_RE"),
  model_fit_row(m_A6_nonlin_re, "A6_nonlin_RE_tp"),
  model_fit_row(m_A6_cr_re,     "A6_nonlin_RE_cr"))
A6_OverallSlopes <- extract_region_gender_pgi_effects(m_A6_null_re, dat_cluster)
.a6_kc <- tryCatch(mgcv::k.check(m_A6_nonlin_re), error = function(e) NULL)
A6_K_Index <- if (!is.null(.a6_kc)) {
  d <- data.frame(model = "A6_nonlin_RE_tp", smooth_term = rownames(.a6_kc),
    k_prime = .a6_kc[, "k'"], edf = .a6_kc[, "edf"], k_index = .a6_kc[, "k-index"],
    p_value = .a6_kc[, "p-value"], stringsAsFactors = FALSE)
  d[!grepl("fid_re", d$smooth_term), , drop = FALSE]
} else NULL
# A6 LOO: nonlinearity LRT + per-fold 4-cell slope curves (one pass).
cat("A6 LOO...\n")
.a6_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster))
.a6_rows  <- list(); .a6_curves <- list()
for (fl in names(.a6_folds)) {
  dl  <- .a6_folds[[fl]]
  sub <- if (length(dl) == 0) dat_cluster else
    dat_cluster[!(as.character(dat_cluster$cohort) %in% dl), , drop = FALSE]
  attr(sub, "by_mean") <- attr(dat_cluster, "by_mean")
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  fold_lab <- if (length(dl) == 0) "Full sample" else paste0("-", fl)
  held     <- if (length(dl) == 0) "(none)" else fl
  m_lin <- mgcv::bam(as.formula(paste(.a6_rhs("lin",    K_DEFAULT, BS_DEFAULT), rt)),
                     data = sub, method = "fREML", discrete = TRUE)
  m_nl  <- mgcv::bam(as.formula(paste(.a6_rhs("nonlin", K_DEFAULT, BS_DEFAULT), rt)),
                     data = sub, method = "fREML", discrete = TRUE)
  v <- lrt_pair(m_lin, m_nl)
  concl <- if (is.na(v[["p"]])) "Not tested"
           else if (v[["p"]] < 0.05) "Nonlinearity supported"
           else "No evidence that nonlinear specification is needed"
  .a6_rows[[fl]] <- data.frame(fold = fold_lab, held_out = held, N = nrow(sub),
    LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1),
    LRT_p = signif(v[["p"]], 4), nonlinearity_conclusion = concl,
    status = "OK", stringsAsFactors = FALSE)
  cv <- dplyr::bind_rows(
    extract_pgi_slope(m_nl, sub, east_west = "West", gender = "female"),
    extract_pgi_slope(m_nl, sub, east_west = "West", gender = "male"),
    extract_pgi_slope(m_nl, sub, east_west = "East", gender = "female"),
    extract_pgi_slope(m_nl, sub, east_west = "East", gender = "male"))
  cv$fold <- fold_lab; .a6_curves[[fl]] <- cv
}
A6_LOO        <- dplyr::bind_rows(.a6_rows)
A6_LOO_curves <- dplyr::bind_rows(.a6_curves)
A6_LOO_curves$fold  <- factor(A6_LOO_curves$fold, levels = unique(A6_LOO_curves$fold))
A6_LOO_curves$group <- paste(A6_LOO_curves$east_west, A6_LOO_curves$gender)

# ---- A0g: overall Region x Gender PGI contrast (RE), full + pre-reunif ----
# No BY x PGI interaction (birth year kept as a covariate). The focal
# PGI_Edu_z:east_west_c:gender_c asks whether the PGI-education slope's gender
# gap differs by region; cell slopes per Region x Gender are reported
# additively.
cat("A0g overall Region x Gender contrast (RE; full + pre-reunif)...\n")
A0g_focal <- "PGI_Edu_z:east_west_c:gender_c"
.a0g_cell_slopes <- function(mod, dat) {
  cells <- list(c("West", "female", -0.5, 0.5), c("West", "male", -0.5, -0.5),
                c("East", "female", 0.5, 0.5), c("East", "male", 0.5, -0.5))
  fl <- if ("fid_re" %in% names(mod$model)) levels(mod$model$fid_re) else NULL
  b <- coef(mod); V <- vcov(mod)
  do.call(rbind, lapply(cells, function(cl) {
    rc <- as.numeric(cl[3]); gc <- as.numeric(cl[4])
    mk <- function(pg) {
      nd <- data.frame(PGI_Edu_z = pg, east_west_c = rc, gender_c = gc, BYc_z = 0)
      if (!is.null(fl)) nd$fid_re <- factor(fl[1], levels = fl)
      predict(mod, newdata = nd, type = "lpmatrix")
    }
    D <- mk(1) - mk(0); e <- as.numeric(D %*% b); s <- sqrt(as.numeric(D %*% V %*% t(D)))
    data.frame(east_west = cl[1], gender = cl[2], group = paste(cl[1], cl[2]),
      slope = e, se = s, ci_lo = e - 1.96 * s, ci_hi = e + 1.96 * s,
      p = 2 * pnorm(-abs(e / s)), stringsAsFactors = FALSE)
  }))
}
.fit_a0g <- function(dat, tag) {
  m <- mgcv::bam(
    as.formula(paste("edu_z_kernel ~ PGI_Edu_z * east_west_c * gender_c + BYc_z", re_term)),
    data = dat, method = "fREML", discrete = TRUE)
  list(coefs = tidy_coefs(m, tag), fit = model_fit_row(m, tag),
       cell_slopes = .a0g_cell_slopes(m, dat), N = nrow(dat))
}
A0g <- .fit_a0g(dat_cluster, "A0g_Edu")
.dc_pre <- dat_cluster[dat_cluster$reunif == "pre", , drop = FALSE]
.dc_pre$fid_re <- droplevels(factor(.dc_pre$fid_re))
A0g_pre <- if (nrow(.dc_pre) >= 100 && length(unique(.dc_pre$east_west)) == 2 &&
               nlevels(.dc_pre$fid_re) >= 2) .fit_a0g(.dc_pre, "A0g_Edu_pre") else NULL

# ---- A0g LOO: leave-one-study-out of the Region x Gender per-cell slopes ----
# Leave-one-study-out analogue of the Fig-4A per-cell region x gender slopes
# (source: A0g$cell_slopes). Per fold, refit the constant-slope A0g model on the
# reduced sample and re-extract the four Region x Gender cell slopes. The fold
# that drops TwinLife uses the no-family-RE term (.loo_re_term), matching A3/A6
# LOO. Emits four rows per fold (Full sample first).
.a0g_folds <- c(list("Full sample" = character(0)), cohort_groups(dat_cluster))
.a0g_loo   <- list()
for (fl in names(.a0g_folds)) {
  dl   <- .a0g_folds[[fl]]; full <- length(dl) == 0
  sub  <- if (full) dat_cluster else
    dat_cluster[!(as.character(dat_cluster$cohort) %in% dl), , drop = FALSE]
  sub$fid_re <- droplevels(factor(sub$fid_re)); rt <- .loo_re_term(sub)
  m <- mgcv::bam(
    as.formula(paste("edu_z_kernel ~ PGI_Edu_z * east_west_c * gender_c + BYc_z", rt)),
    data = sub, method = "fREML", discrete = TRUE)
  cs <- .a0g_cell_slopes(m, sub)
  cs$fold     <- if (full) "Full sample" else paste0("-", fl)
  cs$held_out <- if (full) "(none)" else fl
  cs$N        <- nrow(sub)
  .a0g_loo[[fl]] <- cs
}
A0g_LOO_cells <- dplyr::bind_rows(.a0g_loo)
A0g_LOO_cells$fold <- factor(A0g_LOO_cells$fold, levels = unique(A0g_LOO_cells$fold))

# ---- A-spec: logistic-transition specification check (dedup sample) ----
# The preregistered logistic transition (Eqs. 8-10) is reported as a
# convergence / fit diagnostic; the NLS has no RE analogue, so it is fit on
# dat_edu.
cat("A-spec logistic specification check (NLS, dedup)...\n")
m_Aspec_linear <- lm(education ~ BYc + PGI_Edu + PGI_Edu:BYc + east_west_c + gender_c,
                     data = dat_edu)
.aspec_ll_lin <- as.numeric(logLik(m_Aspec_linear)); .aspec_k_lin <- length(coef(m_Aspec_linear))
.aspec_fit_nls <- function(label, form, start, maxiter) {
  res <- tryCatch(nls(form, data = dat_edu, start = start,
    control = nls.control(maxiter = maxiter, warnOnly = TRUE)), error = function(e) e)
  if (inherits(res, "error"))
    return(data.frame(model = label, converged = FALSE, n_iter = NA_integer_,
      message = paste("error:", conditionMessage(res)), AIC = NA_real_, BIC = NA_real_,
      U_minus_L = NA_real_, LRT_chisq = NA_real_, LRT_df = NA_real_,
      LRT_p_boundary = NA_real_, stringsAsFactors = FALSE))
  ci <- res$convInfo; cf <- coef(res)
  UmL <- if (all(c("U", "L") %in% names(cf))) unname(cf["U"] - cf["L"]) else NA_real_
  ll <- as.numeric(logLik(res)); kk <- length(cf)
  chi <- max(0, 2 * (ll - .aspec_ll_lin)); df <- kk - .aspec_k_lin
  p_b <- if (df > 0) 0.5 * pchisq(chi, df, lower.tail = FALSE) else NA_real_
  data.frame(model = label, converged = isTRUE(ci$isConv),
    n_iter = if (!is.null(ci$finIter)) ci$finIter else NA_integer_,
    message = if (!is.null(ci$stopMessage)) ci$stopMessage else "(none)",
    AIC = round(AIC(res), 2), BIC = round(BIC(res), 2), U_minus_L = round(UmL, 5),
    LRT_chisq = round(chi, 3), LRT_df = df, LRT_p_boundary = signif(p_b, 4),
    stringsAsFactors = FALSE)
}
Aspec_Diagnostic <- dplyr::bind_rows(
  data.frame(model = "Aspec_linear_null", converged = TRUE, n_iter = NA_integer_,
    message = "(OLS reference for the boundary LRT)", AIC = round(AIC(m_Aspec_linear), 2),
    BIC = round(BIC(m_Aspec_linear), 2), U_minus_L = NA_real_, LRT_chisq = NA_real_,
    LRT_df = NA_real_, LRT_p_boundary = NA_real_, stringsAsFactors = FALSE),
  .aspec_fit_nls("Aspec_logistic_simple",
    education ~ alpha + delta_BY * BYc +
      (L + (U - L) / (1 + exp(-kappa * (BYc - tau)))) * PGI_Edu +
      gamma_R * east_west_c + gamma_G * gender_c,
    list(alpha = 0, delta_BY = 0, L = 0.05, U = 0.3, kappa = 0.05, tau = 0,
         gamma_R = 0, gamma_G = 0), 300),
  .aspec_fit_nls("Aspec_logistic_RxG",
    education ~ alpha + delta_BY * BYc +
      (L + (U - L) / (1 + exp(
        -(kappa_0 + kappa_R * east_west_c + kappa_G * gender_c + kappa_RG * east_west_c * gender_c) *
          (BYc - (tau_0 + tau_R * east_west_c + tau_G * gender_c + tau_RG * east_west_c * gender_c))))) * PGI_Edu +
      gamma_R * east_west_c + gamma_G * gender_c,
    list(alpha = 0, delta_BY = 0, L = 0.05, U = 0.3,
         kappa_0 = 0.05, kappa_R = 0, kappa_G = 0, kappa_RG = 0,
         tau_0 = 0, tau_R = 0, tau_G = 0, tau_RG = 0, gamma_R = 0, gamma_G = 0), 500))
.aspec_nls   <- Aspec_Diagnostic[Aspec_Diagnostic$model != "Aspec_linear_null", , drop = FALSE]
.aspec_conv  <- any(.aspec_nls$converged, na.rm = TRUE)
.aspec_beats <- any(!is.na(.aspec_nls$LRT_p_boundary) & .aspec_nls$LRT_p_boundary < 0.05, na.rm = TRUE)
.aspec_degen <- any(!is.na(.aspec_nls$U_minus_L) & abs(.aspec_nls$U_minus_L) < 1e-3, na.rm = TRUE)
Aspec_note <- if (!.aspec_conv) {
  "Logistic-transition NLS did not converge; specification not interpreted."
} else if (.aspec_degen) {
  "Logistic-transition NLS converged onto a degenerate boundary (U = L, no transition); not interpreted."
} else if (!.aspec_beats) {
  "Logistic transition does not improve on the linear model (boundary LRT n.s.); not interpreted."
} else {
  "Logistic transition converged and improves on the linear model (see Aspec_Diagnostic)."
}
cat(sprintf("  A-spec: %s\n", Aspec_note))

# =====================================================================
# Per-cell nonlinearity: clean nested LRTs (add ONLY that cell's smooth
# to the linear model) for A3 (region), A5 (gender), A6 (Region x Gender).
# The summary.gam per-smooth F in the smooths tables is confounded with
# the parametric PGI:BYc linear slope (its null space includes the linear
# term), so it overstates "nonlinearity"; these added-smooth LRTs isolate
# curvature. The display slope curves then draw the nonlinear smooth ONLY
# for cells with significant curvature (LRT p < .05) and the linear slope
# otherwise, so plots never show wiggle the data do not support. B3/B5/B6
# follow the same convention. See ../../../Plan_deviations.md §6.
# =====================================================================
cat("Per-cell nonlinearity LRTs + display curves...\n")
.cell_lrt_table <- function(base_rhs, lin_mod, cell_smooths, dat, rt, tag) {
  do.call(rbind, lapply(names(cell_smooths), function(cell) {
    m <- mgcv::bam(as.formula(paste(base_rhs, "+", cell_smooths[[cell]], rt)),
                   data = dat, method = "fREML", discrete = TRUE)
    v <- lrt_pair(lin_mod, m)
    data.frame(analysis = tag, cell = cell,
      LRT_chisq = round(v[["chisq"]], 3), LRT_df = round(v[["df"]], 1),
      LRT_p = signif(v[["p"]], 4),
      nonlinear_sig = isTRUE(!is.na(v[["p"]]) && v[["p"]] < 0.05),
      stringsAsFactors = FALSE)
  }))
}
.sm <- function(by) sprintf("s(BYc, by = %s, k = %d, bs = '%s')", by, K_DEFAULT, BS_DEFAULT)

# A3 (region): West / East curvature
A3_cell_LRT <- .cell_lrt_table(
  "edu_z_kernel ~ BYc + east_west_c + gender_c + PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc",
  m_A3_lin_re, list(West = .sm("PGI_W"), East = .sm("PGI_E")), dat_cluster, re_term, "A3")
A3_display_curves <- dplyr::bind_rows(lapply(c("West", "East"), function(rg) {
  sig <- A3_cell_LRT$nonlinear_sig[A3_cell_LRT$cell == rg]
  m   <- if (isTRUE(sig)) m_A3_nonlin_re else m_A3_lin_re
  cv  <- extract_pgi_slope(m, dat_cluster, east_west = rg)
  cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv
}))

# A5 (gender, regions pooled): female / male curvature
A5_cell_LRT <- .cell_lrt_table(
  .a5_rhs("lin", K_DEFAULT, BS_DEFAULT), m_A5_lin_re,
  list(female = .sm("PGI_F"), male = .sm("PGI_M")), dat_cluster, re_term, "A5")
A5_curves <- dplyr::bind_rows(lapply(c("female", "male"), function(g) {
  sig <- A5_cell_LRT$nonlinear_sig[A5_cell_LRT$cell == g]
  m   <- if (isTRUE(sig)) m_A5_nonlin_re else m_A5_lin_re
  cv  <- .a5_slope(m, dat_cluster, g)
  cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv
}))

# A5_LOO: leave-one-study-out of the per-gender attainment GAM — per fold, the
# omnibus linear-vs-nonlinear LRT plus the female- and male-specific curvature
# LRTs, with the supported birth-year range. The full-sample A5 result is a
# stable null (no cohort-varying gender pattern in attainment); this checks the
# null is neither masked nor created by any single study. Additional
# sensitivity (../../../Plan_deviations.md §14), reported for symmetry with B5.
A5_LOO <- loco_pergender_gam(
  dat_cluster, rhs_without_re(m_A5_lin_re), rhs_without_re(m_A5_nonlin_re),
  cells = list(female = "PGI_F", male = "PGI_M"), label = "A5")

# A5_LOO_curves: per-fold region-pooled per-gender PGI-slope curves across birth
# year (leave-one-study-out analogue of A5_curves / Fig-4B). The curve is read
# from each fold's nonlinear per-gender fit, matching A3_LOO_curves /
# A6_LOO_curves.
A5_LOO_curves <- loco_pergender_curves(
  dat_cluster, rhs_without_re(m_A5_nonlin_re), .a5_slope)

# A6 (Region x Gender): WF / WM / EF / EM curvature
A6_cell_LRT <- .cell_lrt_table(
  .a6_rhs("lin", K_DEFAULT, BS_DEFAULT), m_A6_lin_re,
  list(`West female` = .sm("PGI_WF"), `West male` = .sm("PGI_WM"),
       `East female` = .sm("PGI_EF"), `East male` = .sm("PGI_EM")),
  dat_cluster, re_term, "A6")
A6_display_curves <- dplyr::bind_rows(lapply(
  list(c("West", "female"), c("West", "male"), c("East", "female"), c("East", "male")),
  function(cl) {
    cellname <- paste(cl[1], cl[2])
    sig <- A6_cell_LRT$nonlinear_sig[A6_cell_LRT$cell == cellname]
    m   <- if (isTRUE(sig)) m_A6_nonlin_re else m_A6_lin_re
    cv  <- extract_pgi_slope(m, dat_cluster, east_west = cl[1], gender = cl[2])
    cv$shape <- if (isTRUE(sig)) "nonlinear" else "linear"; cv
  }))
cat(sprintf("  A3 cells nonlinear: %s | A5: %s | A6: %s\n",
            paste(A3_cell_LRT$cell[A3_cell_LRT$nonlinear_sig], collapse = ",") ,
            paste(A5_cell_LRT$cell[A5_cell_LRT$nonlinear_sig], collapse = ",") ,
            paste(A6_cell_LRT$cell[A6_cell_LRT$nonlinear_sig], collapse = ",")))

cat("Gender block done.\n")

# ---- Save models + intermediates ----
model_cache <- list(
  RUN_TS              = RUN_TS_VAL,
  by_mean             = by_mean,
  EDU_SD              = EDU_SD,
  dat_cluster         = dat_cluster,
  m_A1_re             = m_A1_re,
  m_A2_re             = m_A2_re,
  m_A3_null_re        = m_A3_null_re,
  m_A3_lin_re         = m_A3_lin_re,
  m_A3_nonlin_re      = m_A3_nonlin_re,
  m_A3_cr_re          = m_A3_cr_re,
  # Tables (so A_report.R can rebuild xlsx without recomputing)
  A1_Coefficients     = A1_Coefficients,
  A1_ModelFit         = A1_ModelFit,
  A2_Coefficients     = A2_Coefficients,
  A2_ModelFit         = A2_ModelFit,
  A3_NestedLRT        = A3_NestedLRT,
  A3_Smooths_tp       = A3_Smooths_tp,
  A3_Smooths_cr       = A3_Smooths_cr,
  A_ModelFit          = A_ModelFit,
  A3_K_Index          = A3_K_Index,
  A1_LOO              = A1_LOO,
  A2_LOO              = A2_LOO,
  A3_LOO              = A3_LOO,
  A3_LOO_curves       = A3_LOO_curves,
  desc_cohort         = desc_cohort,
  desc_cohort_region  = desc_cohort_region,
  pgi_correlations    = pgi_correlations,
  attrition_core_primary = attrition_core_primary,
  # ---- Gender block (Section C analyses + per-gender attainment GAM, family RE) ----
  A4_focal            = A4_focal,
  m_A4_re             = m_A4_re,
  A4_Coefficients     = A4_Coefficients,
  A4_ModelFit         = A4_ModelFit,
  A4_LOO              = A4_LOO,
  A4_LOO_gender       = A4_LOO_gender,
  m_A5_lin_re         = m_A5_lin_re,
  m_A5_nonlin_re      = m_A5_nonlin_re,
  m_A5_cr_re          = m_A5_cr_re,
  A5_NestedLRT        = A5_NestedLRT,
  A5_Smooths_tp       = A5_Smooths_tp,
  A5_Smooths_cr       = A5_Smooths_cr,
  A5_ModelFit         = A5_ModelFit,
  A5_curves           = A5_curves,
  A5_cell_LRT         = A5_cell_LRT,
  A5_LOO              = A5_LOO,
  A5_LOO_curves       = A5_LOO_curves,
  A3_cell_LRT         = A3_cell_LRT,
  A3_display_curves   = A3_display_curves,
  m_A6_lin_re         = m_A6_lin_re,
  m_A6_nonlin_re      = m_A6_nonlin_re,
  m_A6_cr_re          = m_A6_cr_re,
  A6_NestedLRT        = A6_NestedLRT,
  A6_Smooths_tp       = A6_Smooths_tp,
  A6_Smooths_pgi      = A6_Smooths_pgi,
  A6_Smooths_cr       = A6_Smooths_cr,
  A6_ModelFit         = A6_ModelFit,
  A6_K_Index          = A6_K_Index,
  A6_OverallSlopes    = A6_OverallSlopes,
  A6_cell_LRT         = A6_cell_LRT,
  A6_display_curves   = A6_display_curves,
  A6_LOO              = A6_LOO,
  A6_LOO_curves       = A6_LOO_curves,
  A0g_focal           = A0g_focal,
  A0g                 = A0g,
  A0g_pre             = A0g_pre,
  A0g_LOO_cells       = A0g_LOO_cells,
  Aspec_Diagnostic    = Aspec_Diagnostic,
  Aspec_note          = Aspec_note
)

OUT_RDS <- file.path(OUT_DIR_A, "A_Education_Models.rds")
saveRDS(model_cache, OUT_RDS)
cat(sprintf("A_analysis: models RDS saved to %s\n", OUT_RDS))
