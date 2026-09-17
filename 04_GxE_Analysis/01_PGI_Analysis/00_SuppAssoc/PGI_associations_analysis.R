# PGI_associations_analysis.R
# Supplementary associations module — analysis stage.
#
# Computes, per PGI (Education / Cognitive / Non-cognitive) x outcome
# (educational attainment, educational mobility):
#   * the standardised main association (slope, 95% CI, p), and
#   * the PGI's INCREMENTAL R^2 (delta R^2 of adding that PGI to the
#     covariate-only model).
#
# Same primary analytic samples and covariate structure as the reported
# main-effect models (region, gender, birth year; + parental education for
# mobility), fitted with OLS so a clean R^2 is available (the RE point
# estimates match to ~1e-2). AGGREGATE OUTPUT ONLY — no individual rows.
#
# The Education row is the "R^2 of the PGI" reported in the manuscript; the
# Cognitive / Non-cognitive rows give their main associations with educational
# attainment and mobility. Moderation of the alternative PGIs is the S1
# alt-PGI focal analysis; this module reports main associations only.
#
# Writes: <run>/manuscript_export/pgi_associations.csv
#
# Usage:
#   RUN_TS=<ts> Rscript 01_PGI_Analysis/00_SuppAssoc/PGI_associations_analysis.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})

# ---- Resolve paths (works under Rscript and source()) ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i); if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine PGI_associations_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL <- get_run_ts()
RUN_DIR    <- file.path(gxe_dir(), "01_PGI_Analysis", "runs", paste0("RUN_", RUN_TS_VAL))
OUT_DIR    <- file.path(RUN_DIR, "manuscript_export")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()

DATA_ROOT        <- data_root()
DATA_DIR         <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")

source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

# ---- association + incremental R^2 (with bootstrap CI) ---------------------
set.seed(20260718)                                   # reproducible bootstrap
BOOT_B <- as.integer(Sys.getenv("BOOT_B", "2000"))   # bootstrap resamples

.dr2 <- function(dat, base_f, full_f)
  summary(stats::lm(full_f, dat))$r.squared - summary(stats::lm(base_f, dat))$r.squared

assoc_row <- function(dat, outcome_col, outcome_label, pgi, covars) {
  need <- c(outcome_col, pgi, covars)
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("missing columns in ", outcome_label, ": ", paste(miss, collapse = ", "))
  dat <- dat[stats::complete.cases(dat[, need]), ]
  base_f <- stats::as.formula(paste(outcome_col, "~", paste(covars, collapse = " + ")))
  full_f <- stats::as.formula(paste(outcome_col, "~", pgi, "+", paste(covars, collapse = " + ")))
  m0 <- stats::lm(base_f, data = dat)
  m1 <- stats::lm(full_f, data = dat)
  co <- summary(m1)$coefficients[pgi, ]
  ci <- stats::confint(m1, pgi)
  dr2 <- unname(summary(m1)$r.squared - summary(m0)$r.squared)
  # Bias-corrected (BC) percentile bootstrap CI for the incremental R^2.
  # Resampling with replacement inflates in-sample R^2, so the plain percentile
  # interval is biased upward; the BC step (Efron) recentres it on the estimate.
  n <- nrow(dat)
  boots <- vapply(seq_len(BOOT_B), function(b)
    .dr2(dat[sample.int(n, n, replace = TRUE), , drop = FALSE], base_f, full_f),
    numeric(1))
  z0 <- stats::qnorm(mean(boots < dr2)); if (!is.finite(z0)) z0 <- 0
  probs <- stats::pnorm(2 * z0 + stats::qnorm(c(.025, .975)))
  r2ci <- stats::quantile(boots, probs, names = FALSE, na.rm = TRUE)
  data.frame(
    pgi = sub("_z$", "", sub("^PGI_", "", pgi)), outcome = outcome_label, n = nrow(dat),
    estimate = unname(co["Estimate"]), se = unname(co["Std. Error"]),
    ci_lo = unname(ci[1]), ci_hi = unname(ci[2]), p = unname(co["Pr(>|t|)"]),
    incremental_r2 = dr2, r2_ci_lo = r2ci[1], r2_ci_hi = r2ci[2],
    stringsAsFactors = FALSE)
}

pgis    <- c("PGI_Edu_z", "PGI_Cog_z", "PGI_nonCog_z")
att_cov <- c("BYc_z", "east_west_c", "gender_c")
mob_cov <- c("BYc_z", "east_west_c", "gender_c", "parental_edu_z_kernel")

att <- do.call(rbind, lapply(pgis, function(p) assoc_row(dat_cluster,     "edu_z_kernel", "attainment", p, att_cov)))
mob <- do.call(rbind, lapply(pgis, function(p) assoc_row(dat_cluster_mob, "mobility",     "mobility",   p, mob_cov)))
out <- rbind(att, mob)

utils::write.csv(out, file.path(OUT_DIR, "pgi_associations.csv"), row.names = FALSE)
cat("Wrote", file.path(OUT_DIR, "pgi_associations.csv"), "\n\n")

# aggregate console summary (safe: model-level statistics only)
disp <- transform(out, estimate = round(estimate, 3), incremental_r2 = round(incremental_r2, 3),
                  r2_ci_lo = round(r2_ci_lo, 3), r2_ci_hi = round(r2_ci_hi, 3), p = signif(p, 3))
print(disp[, c("outcome", "pgi", "n", "estimate", "p", "incremental_r2", "r2_ci_lo", "r2_ci_hi")],
      row.names = FALSE)
