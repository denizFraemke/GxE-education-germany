#!/usr/bin/env Rscript
# export_manuscript_aggregates.R
# ----------------------------------------------------------------------------
# Produce tidy, AGGREGATE-ONLY figure-source and descriptive files for the
# manuscript repository, from the frozen model caches of a completed run.
#
# Design contract (mirrors the manuscript repo's frozen-results workflow):
#   - Reads the run's *_Models.rds / Descriptives_Tables.rds caches and the
#     SNP-heritability TSVs.
#   - Mostly *predicts* fitted curves/slopes from already-fitted model objects
#     (aggregate reporting, not estimation) and reshapes pre-computed aggregate
#     tables. The one exception is Fig 3B: it re-estimates the B0 East-West
#     contrast on the pre-/post-1990 subsamples (mirrors the cached B0_pre fit).
#   - NEVER writes individual-level rows. Every output is a curve, a slope, a
#     smooth, an M(SD) table, or a count. No row of one person's values is
#     ever written or printed.
#
# Outputs -> <run>/manuscript_export/ :
#   fig01_*.csv .. fig04_*.csv, figS_*.csv   figure-source data (see plan)
#   descriptives_table1_overall.csv          Table 1 (overall)
#   descriptives_table2_by_year.csv          Table 2 (by birth-year group)
#   focal_estimates.csv                       consolidated focal + robustness
#   snp_h2*.csv                               copies of the S5 TSVs
#   _provenance.csv                           run id, commit, checksums, date
#   SUPPLEMENTARY_DATA_TRANSFER.md            which workbooks go where
#
# Usage:
#   RUN_TS=20260703_1304 Rscript export_manuscript_aggregates.R
#   Rscript export_manuscript_aggregates.R --run 20260703_1304
#   (DATA_ROOT env respected exactly as in 00_setup/run_context.R)
#
#   EXPORT_ONLY="figS06 gender x PGI,pooled SNP" RUN_TS=... Rscript ...
#     Backfill path: run only the steps whose label contains one of the given
#     substrings, leaving every other already-exported CSV untouched. Use when
#     adding an output to a run whose other files are already frozen upstream.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(mgcv)
  # ggplot2/broom/tibble are attached because GxE_Germany_analysis_helpers.R
  # (sourced below for the fitted-curve builders) also loads the plotting
  # helpers, which build ggplot scales at source time.
  library(ggplot2); library(broom); library(tibble)
})

# ---- Resolve script location (works under Rscript and source()) ------------
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine export_manuscript_aggregates.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
# This module lives at the repo top level (05_manuscript_export/); it reuses the
# GxE analysis setup + helpers (run_context, constants, fitted-curve builders) to
# regenerate figure-source data, so WORKFLOW_DIR points into the 04 analysis stage.
WORKFLOW_DIR <- normalizePath(
  file.path(SECTION_DIR, "..", "04_GxE_Analysis", "01_PGI_Analysis"), mustWork = TRUE)

# ---- Run id: --run <ts> | RUN_TS env | default the frozen run --------------
argv <- commandArgs(trailingOnly = TRUE)
ri <- match("--run", argv)
if (!is.na(ri) && length(argv) >= ri + 1) Sys.setenv(RUN_TS = argv[ri + 1])
if (!nzchar(Sys.getenv("RUN_TS"))) Sys.setenv(RUN_TS = "20260703_1304")

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(WORKFLOW_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

RUN_TS  <- get_run_ts()
EXP_DIR <- output_dir("manuscript_export")   # <run>/manuscript_export/, created
message(sprintf("Exporting manuscript aggregates for RUN_%s\n  -> %s", RUN_TS, EXP_DIR))

# ---- Small IO helpers (write only; report dims, never rows) ----------------
written <- character(0)
w <- function(df, name) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    message("  [skip] ", name, " (empty/NULL)"); return(invisible(NULL))
  }
  path <- file.path(EXP_DIR, name)
  utils::write.csv(df, path, row.names = FALSE)
  written <<- c(written, name)
  message(sprintf("  [ok]   %-42s %d x %d", name, nrow(df), ncol(df)))
  invisible(path)
}
# EXPORT_ONLY: comma-separated substrings. When set, only steps whose label
# contains one of them run; every other step is skipped and its existing CSV in
# the run's manuscript_export/ is left untouched. This is the cheap backfill
# path for adding one output to an already-frozen run without rewriting (and
# re-checksumming) the other hundred-odd files. Unset = run everything.
EXPORT_ONLY <- trimws(strsplit(Sys.getenv("EXPORT_ONLY", ""), ",")[[1]])
EXPORT_ONLY <- EXPORT_ONLY[nzchar(EXPORT_ONLY)]
if (length(EXPORT_ONLY))
  message(sprintf("EXPORT_ONLY set: running only steps matching %s",
                  paste(sQuote(EXPORT_ONLY), collapse = ", ")))
try_step <- function(label, expr) {
  if (length(EXPORT_ONLY) && !any(vapply(EXPORT_ONLY, grepl, logical(1),
                                         x = label, fixed = TRUE))) {
    message("  [skip] ", label, " (EXPORT_ONLY)"); return(invisible(NULL))
  }
  tryCatch(force(expr),
           error = function(e) message("  [FAIL] ", label, ": ", conditionMessage(e)))
}

# ---- Load caches -----------------------------------------------------------
A_PATH <- file.path(output_dir("A_Education"),  "A_Education_Models.rds")
B_PATH <- file.path(output_dir("B_Mobility"),   "B_Mobility_Models.rds")
D_PATH <- file.path(output_dir("00_Descriptives"), "Descriptives_Tables.rds")
for (p in c(A_PATH, B_PATH, D_PATH))
  if (!file.exists(p)) stop("Missing cache: ", p, " (run the analysis first).")

A  <- readRDS(A_PATH)
B  <- readRDS(B_PATH)
tb <- readRDS(D_PATH)
datA <- A$dat_cluster
datB <- B$dat_cluster_mob
attr(datB, "by_mean") <- B$by_mean

# ---- Dedup one-per-family sample, rebuilt on demand -------------------------
# Table 1/2 report N and M(SD) on the full-family RE sample but compute their
# TESTS on the dedup one-adult-per-family sample, for clean independence
# (Descriptives_analysis.R; Plan_deviations §5). The companion tests below are
# tests, so they use the same dedup sample -- but no cache stores it
# (Descriptives_Tables.rds holds only the finished tables), so it is rebuilt
# from the merge. `check_frozen_run.R` verifies that this rebuild still yields
# the N the frozen cache recorded; run it before trusting these outputs.
# Lazy + memoised: only the steps that need it pay the load.
.dedup <- NULL
dedup_sample <- function() {
  if (!is.null(.dedup)) return(.dedup)
  RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()
  DATA_DIR         <- file.path(data_root(), "03_Merge")
  EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
  e <- environment()
  sink(tempfile())                     # load/prepare print diagnostics by design
  suppressMessages(suppressWarnings({
    source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = e)
    source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = e)
  }))
  sink()
  if (!identical(nrow(e$dat_edu), as.integer(tb$N_test)))
    stop(sprintf("dedup sample rebuilt with N = %d but the frozen cache recorded %d -- refusing to export against a different sample",
                 nrow(e$dat_edu), tb$N_test))
  .dedup <<- e$dat_edu
  .dedup
}
source(file.path(WORKFLOW_DIR, "R", "helpers_descriptive_tables.R"), local = TRUE)

# ---- Descriptive kernel-mean smoother (Fig 1): fitted here, exported as -----
#      coordinates only. Groups by `group_var` (gender or east_west); trimmed
#      to the birth-year range where BOTH group levels have >= min_cell
#      observations within +/- half_win years. Band = outcome mean +/- 1 SD.
fit_desc_curves <- function(dat, outcome, label, group_var = "gender",
                            k = 10, n_grid = 200,
                            min_cell = 30, half_win = 2.5, min_year = NA_real_) {
  need <- c(outcome, "birth_year", group_var)
  d <- dat[stats::complete.cases(dat[, need]), need]
  gv <- as.character(d[[group_var]])
  levs <- sort(unique(gv))
  supp_range <- function(g) {
    dg <- d[gv == g, ]
    yrs <- seq(floor(min(dg$birth_year)), ceiling(max(dg$birth_year)))
    ok  <- yrs[vapply(yrs, function(y)
      sum(abs(dg$birth_year - y) <= half_win) >= min_cell, logical(1))]
    if (length(ok) < 2) NULL else range(ok)
  }
  ranges <- Filter(Negate(is.null), lapply(levs, supp_range))
  if (!length(ranges)) return(NULL)
  lo <- max(vapply(ranges, `[`, numeric(1), 1))
  hi <- min(vapply(ranges, `[`, numeric(1), 2))
  if (!is.na(min_year)) lo <- max(lo, min_year)   # optional hard floor
  if (!(hi > lo)) return(NULL)
  grid <- seq(lo, hi, length.out = n_grid)
  out <- lapply(levs, function(g) {
    dg <- d[gv == g, ]
    if (nrow(dg) < 50) return(NULL)
    # mean smoother + an SD smoother on squared residuals: the band is the
    # outcome spread (mean +/- 1 SD), not a CI of the mean.
    mean_fit <- mgcv::gam(
      stats::as.formula(paste0(outcome, " ~ s(birth_year, k = ", k, ")")), data = dg)
    dg$.resid2 <- (dg[[outcome]] - as.numeric(fitted(mean_fit)))^2
    sd_fit <- mgcv::gam(.resid2 ~ s(birth_year, k = k), data = dg)
    mu <- as.numeric(predict(mean_fit, newdata = data.frame(birth_year = grid)))
    s  <- sqrt(pmax(0, as.numeric(predict(sd_fit, newdata = data.frame(birth_year = grid)))))
    data.frame(outcome = label, group_type = group_var, group = g,
               birth_year = grid, fit = mu, sd = s, lo = mu - s, hi = mu + s)
  })
  dplyr::bind_rows(out)
}

# ---- Trim slope curves to each group's supported birth-year range ----------
#      (>= min_cell observations within +/- half_win years; same shared rule
#      as the Fig 1 kernel curves). Keeps GAM slope curves from displaying
#      wide, data-sparse boundary wiggles.
supported_ranges <- function(dat, group_var, min_cell = 30, half_win = 2.5) {
  d <- dat[!is.na(dat$birth_year) & !is.na(dat[[group_var]]),
           c("birth_year", group_var), drop = FALSE]
  gv <- as.character(d[[group_var]])
  levs <- sort(unique(gv))
  setNames(lapply(levs, function(g) {
    dg <- d[gv == g, ]
    yrs <- seq(floor(min(dg$birth_year)), ceiling(max(dg$birth_year)))
    ok  <- yrs[vapply(yrs, function(y)
      sum(abs(dg$birth_year - y) <= half_win) >= min_cell, logical(1))]
    if (length(ok) < 2) c(-Inf, Inf) else range(ok)
  }), levs)
}
trim_curve <- function(curve, dat, group_var, ...) {
  if (is.null(curve) || !all(c(group_var, "birth_year") %in% names(curve))) return(curve)
  rng <- supported_ranges(dat, group_var, ...)
  keep <- vapply(seq_len(nrow(curve)), function(i) {
    r <- rng[[as.character(curve[[group_var]][i])]]
    if (is.null(r)) TRUE else curve$birth_year[i] >= r[1] && curve$birth_year[i] <= r[2]
  }, logical(1))
  curve[keep, , drop = FALSE]
}

# ===========================================================================
# FINDING 0 - Descriptives (Fig 1 + Table 1/2)
# ===========================================================================
try_step("descriptives tables", {
  w(tb$main_overall, "descriptives_table1_overall.csv")
  w(tb$main_by_year, "descriptives_table2_by_year.csv")
  if (!is.null(tb$measures))
    w(data.frame(measure = tb$measures), "descriptives_measure_order.csv")
})

# ---------------------------------------------------------------------------
# Companions to Table 1 / Table 2.
#
# Exploratory and not preregistered. These are SAMPLE-DESCRIPTION statistics,
# not tests of any hypothesis in the analysis plan. With N = 13,049 (11,707 in
# the dedup test sample) almost any contrast reaches significance, and a
# birth-year difference in mean education is a property of the sampling frame --
# the four contributing studies cover different birth-year windows by design --
# not a finding about Germany. The `label` and `caveat` columns carry that
# statement into the manuscript repo so it cannot be dropped on the way to the
# text.
#
# Table 1 and Table 2 already carry Hedges' g, a Welch p and a significance
# marker for BOTH the East-West and the Female-Male contrast, for all seven
# measures, overall and in three birth-year groups. Two descriptive quantities
# are not in them, and only those two are computed here:
#   (a) a test of the BIRTH-YEAR TREND itself (Table 2 bins birth year and
#       compares region/gender within a bin; it never tests across bins);
#   (b) variances -- every existing cell is about means.
# Same sample and conventions as the tables they accompany: dedup
# one-per-family, East - West and Female - Male sign convention, cells gated on
# DESC_MIN_CELL_N.
# ---------------------------------------------------------------------------
DESC_LABEL  <- "exploratory; not preregistered"
DESC_CAVEAT <- paste(
  "Describes the composition of the analytic sample; it does not test a study",
  "hypothesis. At this N nearly any contrast is significant, and the birth-year",
  "trend in particular reflects which studies cover which birth years.")

# ---- (a) Birth-year trend, and whether it differs by gender or region ------
#      Tests the trend behind "higher in more recent birth years", and whether
#      the descriptive gender x birth-year crossover survives an interaction
#      test. One model family answers both.
#      Birth year is scaled per DECADE so b reads in measure-units per decade.
try_step("descriptive birth-year trend tests (exploratory)", {
  d <- dedup_sample()
  d$edu_ceiling <- 100 * as.integer(d$education >= 18)
  d$edu_floor   <- 100 * as.integer(d$education <= 9)
  d$by_dec      <- (suppressWarnings(as.numeric(d$birth_year)) - 1950) / 10
  d$gender_c    <- ifelse(as.character(d$gender)    == "female", 0.5, -0.5)
  d$region_c    <- ifelse(as.character(d$east_west) == "East",   0.5, -0.5)
  measures <- list(
    "Years of education"       = "education",
    "Education (edu_z_kernel)" = "edu_z_kernel",
    "% at ceiling (edu = 18y)" = "edu_ceiling",
    "% at floor (edu <= 9y)"   = "edu_floor",
    "Mobility (kernel-z diff)" = "mobility",
    "PGI-Education (z)"        = "PGI_Edu_z",
    "Parental education (yrs)" = "parental_education")
  pull <- function(m, term, label, model, n) {
    co <- summary(m)$coefficients
    if (!term %in% rownames(co)) return(NULL)
    b <- co[term, "Estimate"]; se <- co[term, "Std. Error"]; p <- co[term, "Pr(>|t|)"]
    data.frame(model = model, term = label, n = n,
               b = b, se = se, ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se,
               p = p,
               sig = if (is.na(p)) "" else if (p < .001) "***" else
                     if (p < .01) "**" else if (p < .05) "*" else "",
               stringsAsFactors = FALSE)
  }
  rows <- lapply(names(measures), function(lab) {
    col <- measures[[lab]]
    if (is.null(d[[col]])) return(NULL)
    dd <- d[is.finite(suppressWarnings(as.numeric(d[[col]]))) & is.finite(d$by_dec), ]
    if (nrow(dd) < DESC_MIN_CELL_N) return(NULL)
    dd$y <- suppressWarnings(as.numeric(dd[[col]]))
    n <- nrow(dd)
    # Marginal trend: the descriptive claim, unadjusted.
    m1 <- stats::lm(y ~ by_dec, data = dd)
    out <- list(pull(m1, "by_dec", "birth year (per decade)", "marginal", n))
    # Moderated: does that trend differ by gender / by region? Both moderators
    # are centred +/-0.5, so the by_dec row stays the average trend.
    if (length(unique(dd$gender_c)) == 2 && length(unique(dd$region_c)) == 2) {
      m2 <- stats::lm(y ~ by_dec * gender_c + by_dec * region_c, data = dd)
      out <- c(out, list(
        pull(m2, "by_dec",          "birth year (per decade), adjusted", "moderated", n),
        pull(m2, "by_dec:gender_c", "birth year x gender (female - male)", "moderated", n),
        pull(m2, "by_dec:region_c", "birth year x region (East - West)",  "moderated", n)))
    }
    out <- Filter(Negate(is.null), out)
    if (!length(out)) return(NULL)
    cbind(Measure = lab, do.call(rbind, out), stringsAsFactors = FALSE)
  })
  trend <- do.call(rbind, Filter(Negate(is.null), rows))
  trend$unit   <- "measure units per decade of birth year"
  trend$label  <- DESC_LABEL
  trend$caveat <- DESC_CAVEAT
  w(trend, "descriptives_birthyear_trend.csv")
})

# ---- (b) Variances by region and by gender --------------------------------
#      Every existing cell is a mean; these add variances.
#      Fligner-Killeen is the repo's variance test (it is what
#      figS06_gender_pgi_scalefree.csv already uses) -- rank-based, so it
#      survives the bounded 7-18 education scale and the skew at the ceiling.
#      Sign convention matches the tables: East - West, Female - Male, so the
#      ratio is var(East)/var(West) and var(female)/var(male).
try_step("descriptive variance tests (exploratory)", {
  d <- dedup_sample()
  d$edu_ceiling <- 100 * as.integer(d$education >= 18)
  d$edu_floor   <- 100 * as.integer(d$education <= 9)
  measures <- list(
    "Years of education"       = "education",
    "Education (edu_z_kernel)" = "edu_z_kernel",
    "% at ceiling (edu = 18y)" = "edu_ceiling",
    "% at floor (edu <= 9y)"   = "edu_floor",
    "Mobility (kernel-z diff)" = "mobility",
    "PGI-Education (z)"        = "PGI_Edu_z",
    "Parental education (yrs)" = "parental_education")
  # A 0/100 indicator's variance is a deterministic function of its mean, so a
  # variance test on it restates the mean test. SDs are still reported (they
  # are read as the share), but no p-value is emitted for one.
  BINARY <- c("edu_ceiling", "edu_floor")
  cell <- function(x, g, lvl1, lvl2, col) {
    keep <- is.finite(x)
    x <- x[keep]; g <- as.character(g)[keep]
    x1 <- x[g == lvl1]; x2 <- x[g == lvl2]
    if (length(x1) < DESC_MIN_CELL_N || length(x2) < DESC_MIN_CELL_N) return(NULL)
    v1 <- stats::var(x1); v2 <- stats::var(x2)
    binary <- col %in% BINARY
    p <- if (binary) NA_real_ else tryCatch(
      stats::fligner.test(c(x1, x2), factor(rep(c(lvl1, lvl2), c(length(x1), length(x2)))))$p.value,
      error = function(e) NA_real_)
    data.frame(
      group1 = lvl1, group2 = lvl2,
      n_group1 = length(x1), sd_group1 = sqrt(v1),
      n_group2 = length(x2), sd_group2 = sqrt(v2),
      var_ratio_g1_over_g2 = if (is.finite(v2) && v2 > 0) v1 / v2 else NA_real_,
      fligner_p = p,
      sig = if (is.na(p)) "" else if (p < .001) "***" else
            if (p < .01) "**" else if (p < .05) "*" else "",
      note = if (binary)
        "no test: a 0/100 indicator's variance is a function of its mean, so this would restate the mean contrast" else "",
      stringsAsFactors = FALSE)
  }
  strata <- c("Overall", levels(by_group3(d$birth_year)))
  bg <- by_group3(d$birth_year)
  rows <- list()
  for (st in strata) {
    ds <- if (st == "Overall") d else d[which(bg == st), , drop = FALSE]
    for (lab in names(measures)) {
      col <- measures[[lab]]
      if (is.null(ds[[col]])) next
      x <- suppressWarnings(as.numeric(ds[[col]]))
      for (cn in list(list("East-West",   ds$east_west, "East",   "West"),
                      list("Female-Male", ds$gender,    "female", "male"))) {
        r <- cell(x, cn[[2]], cn[[3]], cn[[4]], col)
        if (is.null(r)) next
        key <- data.frame(`Birth-year group` = st, Measure = lab, contrast = cn[[1]],
                          check.names = FALSE, stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- cbind(key, r)
      }
    }
  }
  vt <- do.call(rbind, rows)
  # Two of these cells have a frozen counterpart computed on a DIFFERENT scale.
  # Say so in the file: a reader comparing them side by side would otherwise
  # read agreement as duplication and disagreement as an error.
  vt$cross_reference <- ""
  is_ov <- vt$`Birth-year group` == "Overall" & vt$contrast == "East-West"
  vt$cross_reference[is_ov & vt$Measure == "Years of education"] <- paste(
    "S4c (attainment/sensitivities/S4_Migration_Tests.csv) tests the same",
    "East-West education-variance contrast on the same sample as a TOST",
    "equivalence test; its descriptive ratio agrees with this row.")
  vt$cross_reference[is_ov & vt$Measure == "PGI-Education (z)"] <- paste(
    "DIFFERENT SCALE from S4b in attainment/sensitivities/S4_Migration_Tests.csv,",
    "which finds East PGI variance BELOW West (ratio 0.88, p = 2.1e-05). S4b uses",
    "the RAW PGI_Edu; this row uses PGI_Edu_z, residualised on the ten pooled",
    "cross-cohort ancestry PCs and then re-z-scored. The East-West variance gap in",
    "the PGI is therefore largely absorbed by the ancestry-PC adjustment. Both are",
    "correct on their own scale; they are not interchangeable, and the migration",
    "diagnostic in the plan (S4 Test 2) is the raw-PGI one.")
  vt$label  <- DESC_LABEL
  vt$caveat <- DESC_CAVEAT
  w(vt, "descriptives_variance_tests.csv")
})
try_step("fig01 kernel curves (gender + region)", {
  curves <- dplyr::bind_rows(
    fit_desc_curves(datA, "education", "Years of education", "gender"),
    fit_desc_curves(datA, "education", "Years of education", "east_west"),
    fit_desc_curves(datB, "mobility",  "Educational mobility", "gender",    min_year = 1935),
    fit_desc_curves(datB, "mobility",  "Educational mobility", "east_west", min_year = 1935)
  )
  w(curves, "fig01_kernel_curves.csv")
})

# ===========================================================================
# FINDING 1 - Attainment: regionally stable PGI slope (Fig 2)
# ===========================================================================
try_step("fig02 attainment slope curves", {
  curves <- if (!is.null(A$A3_display_curves)) A$A3_display_curves else
    dplyr::bind_rows(
      extract_pgi_slope(A$m_A3_nonlin_re, datA, east_west = "West"),
      extract_pgi_slope(A$m_A3_nonlin_re, datA, east_west = "East"))
  w(trim_curve(curves, datA, "east_west"), "fig02_attainment_slope_curves.csv")
})
try_step("fig02 attainment np overlay", {
  np <- compute_binned_slopes(datA, outcome = "edu_z_kernel", pgi = "PGI_Edu_z",
                              group_vars = "east_west")
  w(np, "fig02_attainment_np_slopes.csv")
})
try_step("fig02 birth-year histogram (counts)", {
  h <- datA %>% filter(!is.na(birth_year), !is.na(east_west)) %>%
    mutate(by_bin = round(birth_year)) %>%
    count(east_west, by_bin, name = "n")
  w(as.data.frame(h), "fig02_birthyear_counts.csv")
})

# ---------------------------------------------------------------------------
# SUPPLEMENT (supports Finding 1) - SOEP-only East/West PGI-Education attainment
# analysis using the WITHIN-cohort PC-residualised PGI (Fig S5).
# ---------------------------------------------------------------------------
# SOEP-G is the only non-SHIP cohort with both East and West participants above
# the per-region threshold (n >= 30), so it is the natural single-cohort test of
# the East-West pattern originally reported in SOEP-G alone (Fraemke et al., 2025).
# That paper adjusted the PGI on WITHIN-SOEP ancestry PCs; the pooled multi-cohort
# pipeline here instead residualises on the shared cross-cohort PCs (crossPC1-10).
# For the single-cohort comparison we reconstruct the 2025-style PGI: residualise
# the raw PGI on the within-cohort (SOEP) PCs PC1-10 - same PC count as the
# pipeline, only the basis changes - and standardise within SOEP. This is a
# documented refit (like the Fig 3B B0 re-estimation), from the frozen dat_cluster
# so the sample is identical. Aggregate output only.
try_step("figS05 SOEP-only within-cohort-PC region slopes + three-way", {
  dSO <- datA[as.character(datA$cohort) == "SOEP", , drop = FALSE]
  pcs <- paste0("PC", 1:10)
  need <- c("PGI_Edu_raw", pcs, "edu_z_kernel", "birth_year", "east",
            "east_west", "east_west_c", "gender_c", "BYc_z")
  if (!all(need %in% names(dSO)))
    stop("dat_cluster lacks within-PC inputs (need PGI_Edu_raw + PC1-10): ",
         paste(setdiff(need, names(dSO)), collapse = ", "))
  # keep the canonical cross-cohort-PC PGI (for the specification contrast) before
  # overwriting PGI_Edu_z with the within-cohort-PC version used by the fits below.
  dSO$PGI_cross <- dSO$PGI_Edu_z
  # within-SOEP PC-residualised, within-SOEP-standardised PGI (overwrite PGI_Edu_z,
  # which the GAM/linear fits below read as the focal PGI).
  ok <- stats::complete.cases(dSO[, c("PGI_Edu_raw", pcs)])
  r  <- rep(NA_real_, nrow(dSO))
  r[ok] <- stats::residuals(
    lm(stats::as.formula(paste("PGI_Edu_raw ~", paste(pcs, collapse = "+"))),
       data = dSO[ok, , drop = FALSE]))
  dSO$PGI_Edu_z <- as.numeric(scale(r))

  # --- region varying-coefficient GAM for the slope curves (mirrors the appendix
  #     pc_region_gam 2-region branch, but on the within-PC PGI) ---
  d <- dSO[!is.na(dSO$edu_z_kernel) & !is.na(dSO$PGI_Edu_z) &
           !is.na(dSO$birth_year)   & !is.na(dSO$east), , drop = FALSE]
  nby <- length(unique(d$birth_year)); k <- min(K_DEFAULT, max(4L, floor(nby / 3)))
  bym <- mean(d$birth_year); d$BYc <- d$birth_year - bym; attr(d, "by_mean") <- bym
  d$PGI_W <- d$PGI_Edu_z * (1L - d$east); d$PGI_E <- d$PGI_Edu_z * d$east
  rhs_lin <- "edu_z_kernel ~ BYc + east_west_c + gender_c + PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc"
  rhs_nl  <- sprintf("%s + s(BYc, by = PGI_W, k = %d, bs = '%s') + s(BYc, by = PGI_E, k = %d, bs = '%s')",
                     rhs_lin, k, BS_DEFAULT, k, BS_DEFAULT)
  m_lin <- mgcv::gam(stats::as.formula(rhs_lin), data = d, method = "REML")
  m_nl  <- mgcv::gam(stats::as.formula(rhs_nl),  data = d, method = "REML")
  st <- as.data.frame(summary(m_nl)$s.table)
  smp <- function(pat) { rr <- st[grepl(pat, rownames(st), fixed = TRUE), , drop = FALSE]
    if (nrow(rr)) signif(rr[1, "p-value"], 3) else NA_real_ }
  clip <- function(crv, by_obs, lo = 0.02, hi = 0.98) {
    qs <- quantile(by_obs, c(lo, hi), na.rm = TRUE)
    crv[crv$birth_year >= qs[1] & crv$birth_year <= qs[2], , drop = FALSE] }
  curves <- dplyr::bind_rows(lapply(c("West", "East"), function(reg) {
    crv <- extract_pgi_slope(m_nl, d, east_west = reg, n_grid = 150)
    clip(crv, d$birth_year[if (reg == "East") d$east == 1 else d$east == 0]) }))
  reg_shape <- c(West = if (isTRUE(smp("PGI_W") < .05)) "nonlinear" else "linear",
                 East = if (isTRUE(smp("PGI_E") < .05)) "nonlinear" else "linear")
  curves$shape <- unname(reg_shape[as.character(curves$east_west)])
  cc <- intersect(c("birth_year","BYc","east_west","gender","slope","se","ci_lo","ci_hi","shape"),
                  names(curves))
  curves <- curves[order(as.character(curves$east_west), curves$birth_year), cc, drop = FALSE]
  w(curves, "figS05_soep_region_slope_curves.csv")

  # --- linear focal tests on the within-PC PGI: three-way + per-region amplifications ---
  cw <- function(m, term) { pt <- summary(m)$p.table; e <- pt[term,1]; s <- pt[term,2]
    c(b = e, lo = e - 1.96*s, hi = e + 1.96*s, p = pt[term,4]) }
  m3 <- mgcv::bam(edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3,
                  data = dSO, method = "REML")
  tw <- cw(m3, "PGI_Edu_z:BYc_z:east_west_c")
  # canonical cross-cohort-PC PGI, same three-way (the specification contrast)
  m3c <- mgcv::bam(edu_z_kernel ~ (PGI_cross + BYc_z + east_west_c + gender_c)^3,
                   data = dSO, method = "REML")
  twc <- cw(m3c, "PGI_cross:BYc_z:east_west_c")
  dSO$reunif <- droplevels(factor(dSO$reunif))
  twr <- tryCatch({
    mr <- mgcv::bam(edu_z_kernel ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3,
                    data = dSO, method = "REML")
    cw(mr, grep("^PGI_Edu_z:reunif.*:east_west_c$", rownames(summary(mr)$p.table), value = TRUE)[1])
  }, error = function(e) c(b = NA, lo = NA, hi = NA, p = NA))
  amp <- function(reg) { dr <- dSO[dSO$east_west == reg, , drop = FALSE]
    cw(mgcv::bam(edu_z_kernel ~ PGI_Edu_z * BYc_z + gender_c, data = dr, method = "REML"),
       "PGI_Edu_z:BYc_z") }
  aE <- amp("East"); aW <- amp("West")

  endpt <- function(reg, a) { s <- curves[as.character(curves$east_west) == reg, , drop = FALSE]
    s <- s[order(s$birth_year), , drop = FALSE]; e <- s[1, ]; l <- s[nrow(s), ]
    data.frame(region = reg,
      n = if (reg == "East") sum(d$east == 1) else sum(d$east == 0), n_total = nrow(d),
      shape = unname(reg_shape[reg]), smooth_p = if (reg == "East") smp("PGI_E") else smp("PGI_W"),
      by_early = e$birth_year, slope_early = e$slope, ci_lo_early = e$ci_lo, ci_hi_early = e$ci_hi,
      by_late = l$birth_year, slope_late = l$slope, ci_lo_late = l$ci_lo, ci_hi_late = l$ci_hi,
      ampl_b = a["b"], ampl_ci_lo = a["lo"], ampl_ci_hi = a["hi"], ampl_p = a["p"],
      stringsAsFactors = FALSE, row.names = NULL) }
  w(dplyr::bind_rows(endpt("East", aE), endpt("West", aW)), "figS05_soep_region_stats.csv")

  within_lab <- "within-cohort PC-residualised (10 SOEP PCs), within-SOEP z"
  cross_lab  <- "cross-cohort PC-residualised (canonical pooled pipeline)"
  w(data.frame(
      contrast = c("PGI x birth-year x region",
                   "PGI x birth-year x region",
                   "PGI x reunification-step x region"),
      construction = c(within_lab, cross_lab, within_lab),
      b     = c(tw["b"],  twc["b"],  twr["b"]),
      ci_lo = c(tw["lo"], twc["lo"], twr["lo"]),
      ci_hi = c(tw["hi"], twc["hi"], twr["hi"]),
      p     = c(tw["p"],  twc["p"],  twr["p"]),
      stringsAsFactors = FALSE, row.names = NULL),
    "figS05_soep_threeway.csv")
})

# ===========================================================================
# FINDING 2 - Mobility: PGI slope stronger in the East (Fig 3)
# ===========================================================================
try_step("fig03 mobility slope curves", {
  # Per-region display curves: the nonlinear smooth is drawn ONLY where the
  # per-region curvature LRT is significant (West, p<.001); the other region
  # is shown as its linear trend (East curvature ns, p=.60). This is what
  # B3_display_curves encodes (shape column = nonlinear/linear).
  curves <- if (!is.null(B$B3_display_curves)) as.data.frame(B$B3_display_curves) else
    dplyr::bind_rows(
      extract_pgi_slope(B$m_B3_nonlin_re, datB, east_west = "West"),
      extract_pgi_slope(B$m_B3_nonlin_re, datB, east_west = "East"))
  w(trim_curve(curves, datB, "east_west"), "fig03_mobility_slope_curves.csv")
})
try_step("fig03 mobility East-West difference smooth", {
  w(diff_smooth_simul(B$m_B3_nonlin_re, datB), "fig03_mobility_diffsmooth.csv")
})
try_step("fig03 mobility per-region slopes (B0 pre vs post)", {
  # B0 East-West contrast re-estimated on the pre- and post-reunification
  # subsamples (mirrors the cached B0_pre analysis) so Fig 3B shows the gap
  # both predates and persists across 1990. This is the one place the export
  # re-fits rather than only predicting from cached models.
  b0_f <- paste(
    "mobility ~ PGI_Edu_z * east_west_c +",
    "east_west_c * (gender_c + BYc_z + parental_edu_z_kernel) +",
    "gender_c * BYc_z + gender_c * parental_edu_z_kernel +",
    "BYc_z * parental_edu_z_kernel + s(fid_re, bs = 're')")
  fit_period <- function(period, label) {
    sub <- datB[!is.na(datB$reunif) & datB$reunif == period, , drop = FALSE]
    if (nrow(sub) < 100 || length(unique(sub$east_west)) < 2) return(NULL)
    m <- mgcv::bam(stats::as.formula(b0_f), data = sub, method = "fREML", discrete = TRUE)
    dplyr::mutate(as.data.frame(extract_lm_region_pgi_effects(m, sub)), period = label)
  }
  rs <- dplyr::bind_rows(fit_period("pre", "Pre-1990"), fit_period("post", "Post-1990"))
  w(rs, "fig03_mobility_region_slopes.csv")
})
try_step("fig03 mobility per-cell slopes (B0g)", {
  cs <- dplyr::bind_rows(
    dplyr::mutate(as.data.frame(B$B0g$cell_slopes), sample = "full"),
    if (!is.null(B$B0g_pre))
      dplyr::mutate(as.data.frame(B$B0g_pre$cell_slopes), sample = "pre1990"))
  w(cs, "fig03_mobility_percell_slopes.csv")
})

# ===========================================================================
# FINDING 3 - Gender: women < men PGI association (Fig 4)
# ===========================================================================
try_step("fig04 attainment per-gender PGI slope curves (A5)", {
  w(trim_curve(as.data.frame(A$A5_curves), datA, "gender"),
    "fig04_attainment_gender_slopes.csv")
})
try_step("fig04 mobility per-gender PGI slope curves (B5)", {
  w(trim_curve(as.data.frame(B$B5_curves), datB, "gender"),
    "fig04_mobility_gender_slopes.csv")
})
try_step("fig04 mobility female-male difference smooth (B5, simultaneous CI)", {
  # Planned female-male contrast (B5): diff = female - male PGI-mobility slope with a
  # 95% simultaneous CI. Mirrors the East-West fig03_mobility_diffsmooth export.
  w(diff_smooth_pooled_gender_simul(B$m_B5_nonlin_re, datB),
    "fig04_mobility_gender_diffsmooth.csv")
})
try_step("fig04 per-cell region x gender slopes (A0g + B0g)", {
  cells <- dplyr::bind_rows(
    dplyr::mutate(as.data.frame(A$A0g$cell_slopes), outcome = "Years of education"),
    dplyr::mutate(as.data.frame(B$B0g$cell_slopes), outcome = "Educational mobility"))
  w(cells, "fig04_percell_slopes.csv")
})

# ---- figS06: is the gender gap driven by low- or high-PGI individuals? ------
#      Exploratory; not preregistered. Descriptive raw means only -- NO model:
#      split the attainment sample into PGI-Education quintiles (pooled across
#      region/cohort) and report, per quintile, the male and female mean of
#      relative years of education (edu_z_kernel) with SEs, and the
#      male-minus-female gap. A gap that grows toward Q5 = the gender difference
#      is concentrated at high PGI. Aggregate cells only.
#
#      The same quintile means are also reported on the RAW years-of-education
#      scale (`education`, 7-18), so the gap can be stated in years rather than
#      in SD units of the within-birth-year kernel z.
try_step("figS06 gender x PGI-quintile attainment gap (raw, descriptive)", {
  d <- datA[!is.na(datA$PGI_Edu_z) & !is.na(datA$edu_z_kernel) &
              !is.na(datA$education) & !is.na(datA$gender), ]
  br <- stats::quantile(d$PGI_Edu_z, probs = seq(0, 1, 0.2), na.rm = TRUE)
  d$q <- cut(d$PGI_Edu_z, breaks = br, include.lowest = TRUE, labels = 1:5)
  agg <- do.call(rbind, lapply(split(d, list(d$q, as.character(d$gender))), function(s) {
    if (!nrow(s)) return(NULL)
    data.frame(quintile = as.integer(as.character(s$q[1])),
               gender   = as.character(s$gender[1]),
               n        = nrow(s),
               mean     = mean(s$edu_z_kernel),
               se       = stats::sd(s$edu_z_kernel) / sqrt(nrow(s)),
               # same cell on the raw years-of-education scale (7-18)
               mean_years = mean(s$education),
               se_years   = stats::sd(s$education) / sqrt(nrow(s)),
               pgi_mid  = mean(s$PGI_Edu_z),
               # share at the education ceiling (18 yrs) / floor (<= 9 yrs), raw years
               pct_ceiling = 100 * mean(s$education >= 18),
               pct_floor   = 100 * mean(s$education <= 9))
  }))
  wide <- stats::reshape(
    agg[, c("quintile", "gender", "n", "mean", "se", "mean_years", "se_years",
            "pct_ceiling", "pct_floor")],
    idvar = "quintile", timevar = "gender", direction = "wide")
  wide <- merge(wide, stats::aggregate(pgi_mid ~ quintile, agg, mean), by = "quintile")
  wide$gap    <- wide$mean.male - wide$mean.female
  wide$se_gap <- sqrt(wide$se.male^2 + wide$se.female^2)
  wide$gap_years    <- wide$mean_years.male - wide$mean_years.female
  wide$se_gap_years <- sqrt(wide$se_years.male^2 + wide$se_years.female^2)
  wide <- wide[order(wide$quintile),
               c("quintile", "pgi_mid", "n.female", "mean.female", "se.female",
                 "n.male", "mean.male", "se.male", "gap", "se_gap",
                 "pct_ceiling.female", "pct_ceiling.male",
                 "pct_floor.female", "pct_floor.male",
                 "mean_years.female", "se_years.female",
                 "mean_years.male", "se_years.male", "gap_years", "se_gap_years")]
  names(wide) <- c("quintile", "pgi_mid", "n_female", "mean_female", "se_female",
                   "n_male", "mean_male", "se_male", "gap", "se_gap",
                   "pct_ceiling_female", "pct_ceiling_male",
                   "pct_floor_female", "pct_floor_male",
                   "mean_years_female", "se_years_female",
                   "mean_years_male", "se_years_male", "gap_years", "se_gap_years")
  w(wide, "figS06_gender_pgi_quintiles.csv")
})

# ---- figS06 top-tail band: the same descriptive gap collapsed to PGI bands --
#      The gap is stated for the top 40% (Q4 and Q5), which is a two-quintile
#      band, not a quintile.
#      Reported for the bottom 40% (Q1+Q2), middle 20% (Q3) and top 40% (Q4+Q5)
#      so that "the gap is larger at higher PGI" can be read off the same table
#      rather than eyeballed across the five quintile rows. Descriptive raw
#      means only -- no model and no test. Aggregate only.
try_step("figS06 gender x PGI top-tail band gap (raw, descriptive)", {
  d <- datA[!is.na(datA$PGI_Edu_z) & !is.na(datA$edu_z_kernel) &
              !is.na(datA$education) & !is.na(datA$gender), ]
  br <- stats::quantile(d$PGI_Edu_z, probs = seq(0, 1, 0.2), na.rm = TRUE)
  q  <- cut(d$PGI_Edu_z, breaks = br, include.lowest = TRUE, labels = 1:5)
  d$band <- factor(c("bottom40", "bottom40", "middle20", "top40", "top40")[as.integer(q)],
                   levels = c("bottom40", "middle20", "top40"))
  one <- function(s, band) {
    f <- s[s$gender == "female", ]; m <- s[s$gender == "male", ]
    se <- function(x) stats::sd(x) / sqrt(length(x))
    data.frame(
      band        = band,
      quintiles   = c(bottom40 = "Q1-Q2", middle20 = "Q3", top40 = "Q4-Q5")[[band]],
      pgi_mid     = mean(s$PGI_Edu_z),
      n_female    = nrow(f), n_male = nrow(m),
      mean_female = mean(f$edu_z_kernel),       se_female = se(f$edu_z_kernel),
      mean_male   = mean(m$edu_z_kernel),       se_male   = se(m$edu_z_kernel),
      gap         = mean(m$edu_z_kernel) - mean(f$edu_z_kernel),
      se_gap      = sqrt(se(m$edu_z_kernel)^2 + se(f$edu_z_kernel)^2),
      mean_years_female = mean(f$education),    se_years_female = se(f$education),
      mean_years_male   = mean(m$education),    se_years_male   = se(m$education),
      gap_years   = mean(m$education) - mean(f$education),
      se_gap_years = sqrt(se(m$education)^2 + se(f$education)^2),
      pct_ceiling_female = 100 * mean(f$education >= 18),
      pct_ceiling_male   = 100 * mean(m$education >= 18),
      stringsAsFactors = FALSE)
  }
  bands <- do.call(rbind, lapply(levels(d$band), function(b) one(d[d$band == b, ], b)))
  w(bands, "figS06_gender_pgi_toptail.csv")
})

# ---- figS06 scale-free robustness + gender EA variances. Descriptive gender
#      differences in the variance of educational attainment (male:female variance
#      ratio on raw years and within-cohort-z scales + Fligner-Killeen tests, and
#      per-gender SDs), then two measures that are invariant to that variance and to
#      the bounded 18-year scale: (1) the Spearman rank correlation PGI x raw
#      education per gender (+ Fisher-z test), and (2) an ordered-probit PGI x gender
#      interaction (ordered education; birth year + region adjusted; no family RE --
#      polr limitation). Aggregate scalars only.
try_step("figS06 gender x PGI scale-free robustness + EA variances", {
  d <- datA[!is.na(datA$PGI_Edu_z) & !is.na(datA$education) &
              !is.na(datA$edu_z_kernel) & !is.na(datA$gender), ]
  sp <- function(g) {
    s <- d[d$gender == g, ]; n <- nrow(s)
    rho <- suppressWarnings(stats::cor(s$PGI_Edu_z, s$education, method = "spearman"))
    ci <- tanh(atanh(rho) + c(-1.96, 1.96) / sqrt(n - 3))
    c(rho = rho, lo = ci[1], hi = ci[2], n = n)
  }
  m <- sp("male"); f <- sp("female")
  zf <- (atanh(m["rho"]) - atanh(f["rho"])) / sqrt(1 / (m["n"] - 3) + 1 / (f["n"] - 3))
  dd <- d; dd$edu_ord <- factor(dd$education, ordered = TRUE)
  dd$gender_c <- ifelse(dd$gender == "male", 0.5, -0.5)
  om <- MASS::polr(edu_ord ~ PGI_Edu_z * gender_c + BYc_z + east_west_c,
                   data = dd, method = "probit", Hess = TRUE)
  r <- summary(om)$coefficients["PGI_Edu_z:gender_c", ]
  # gender EA variances (male:female variance ratio + robust Fligner-Killeen test)
  gv <- function(v) {
    mm <- d[[v]][d$gender == "male"]; ff <- d[[v]][d$gender == "female"]
    c(ratio = stats::var(mm) / stats::var(ff),
      sd_m = stats::sd(mm), sd_f = stats::sd(ff),
      flig = stats::fligner.test(d[[v]], factor(d$gender))$p.value)
  }
  vr <- gv("education"); vz <- gv("edu_z_kernel")
  kv <- data.frame(
    key = c("rho_male", "rho_male_lo", "rho_male_hi", "n_male",
            "rho_female", "rho_female_lo", "rho_female_hi", "n_female",
            "fisherz", "fisherz_p", "ord_b", "ord_se", "ord_p", "ord_n",
            "var_ratio_raw", "sd_male_raw", "sd_female_raw", "fligner_p_raw",
            "var_ratio_z", "sd_male_z", "sd_female_z", "fligner_p_z"),
    value = c(m["rho"], m["lo"], m["hi"], m["n"],
              f["rho"], f["lo"], f["hi"], f["n"],
              zf, 2 * pnorm(-abs(zf)),
              r["Value"], r["Std. Error"], 2 * pnorm(-abs(r["t value"])), nrow(d),
              vr["ratio"], vr["sd_m"], vr["sd_f"], vr["flig"],
              vz["ratio"], vz["sd_m"], vz["sd_f"], vz["flig"]),
    row.names = NULL)
  w(kv, "figS06_gender_pgi_scalefree.csv")
})

# ===========================================================================
# SUPPLEMENT - SNP-heritability (Fig S1) + focal forest (Fig S2)
# ===========================================================================
try_step("SNP-heritability copies", {
  h2_dir <- file.path(gxe_dir(), "02_SNP_Heritability", "output", "final")
  snp_map <- c(h2_snp = "snp_h2.csv", heterogeneity = "snp_heterogeneity.csv",
               power = "snp_power.csv",
               h2_region_main = "snp_h2_region.csv",
               h2_region_main_heterogeneity = "snp_h2_region_heterogeneity.csv",
               h2_gender_main = "snp_h2_gender.csv",
               h2_gender_main_heterogeneity = "snp_h2_gender_heterogeneity.csv")
  for (nm in names(snp_map)) {
    src <- file.path(h2_dir, paste0(nm, ".tsv"))
    if (file.exists(src))
      w(utils::read.delim(src, check.names = FALSE), snp_map[[nm]])
    else message("  [skip] ", nm, ".tsv not found at ", h2_dir)
  }
})

# ---- Pooled SNP-heritability of educational attainment ---------------------
#      heterogeneity.tsv carries `h2_pooled` (the inverse-variance combination
#      of the four region x era cell estimates, computed in
#      02_SNP_Heritability/scripts/07_lrt.R) but no standard error for it, so
#      the estimate cannot be reported with an interval from that file alone.
#      The pooled SE follows from the same weights:
#      SE_pooled = 1/sqrt(sum(w)), w_k = 1/SE_k^2.
#
#      This is computed here rather than in 07_lrt.R because 07_lrt.R is a
#      Tardis-side script that reads the per-cell .hsq fits; those files are not
#      transferred back, so it cannot be re-run locally. The formula below is a
#      deliberate copy of 07_lrt.R lines 122-124 and is checked against the
#      frozen `h2_pooled` value on every run (see the stopifnot).
#
#      NOT a whole-sample GREML fit: no such fit exists (the strata are the four
#      cells plus the secondary R1/G1/RG1 tracks, Plan_deviations.md §7h). The
#      cells differ greatly in precision, so this is a precision-weighted
#      combination and not a population heritability for Germany.
try_step("pooled SNP-heritability (inverse-variance over the four cells)", {
  h2_dir <- file.path(gxe_dir(), "02_SNP_Heritability", "output", "final")
  src <- file.path(h2_dir, "h2_snp.tsv")
  het <- file.path(h2_dir, "heterogeneity.tsv")
  if (!file.exists(src)) { message("  [skip] h2_snp.tsv not found at ", h2_dir); return(invisible()) }
  h2 <- utils::read.delim(src, check.names = FALSE)
  hd <- if (file.exists(het)) utils::read.delim(het, check.names = FALSE) else NULL

  pool_one <- function(model, suffix) {
    r <- h2[h2$model == model & !is.na(h2$h2_SNP) & !is.na(h2$h2_SNP_SE) &
              h2$h2_SNP_SE > 0, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    wts <- 1 / r$h2_SNP_SE^2
    pooled <- sum(wts * r$h2_SNP) / sum(wts)
    se     <- sqrt(1 / sum(wts))
    data.frame(
      analysis      = r$analysis[1],
      model_suffix  = suffix,
      K             = nrow(r),
      cells         = paste(r$stratum, collapse = ","),
      N_cells_total = sum(r$N),
      h2_pooled     = pooled,
      h2_pooled_SE  = se,
      ci_lo         = pooled - 1.96 * se,
      ci_hi         = pooled + 1.96 * se,
      note          = paste("inverse-variance combination of the K per-stratum",
                            "standalone GREML h2 estimates; not a whole-sample fit"),
      stringsAsFactors = FALSE)
  }
  out <- rbind(pool_one("per_stratum_standalone", "constrained"),
               pool_one("per_stratum_standalone_noconstrain", "noconstrain"))
  if (is.null(out)) { message("  [skip] no usable per-stratum h2 rows"); return(invisible()) }

  # Reproduction check against the value 07_lrt.R already wrote, so the two
  # implementations cannot silently drift apart.
  if (!is.null(hd) && "h2_pooled" %in% names(hd)) {
    ref <- hd[hd$test == "cochrans_Q", ]
    for (i in seq_len(nrow(out))) {
      r <- ref[ref$model_suffix == out$model_suffix[i], ]
      if (nrow(r))
        stopifnot(isTRUE(all.equal(out$h2_pooled[i], r$h2_pooled[1], tolerance = 1e-8)))
    }
    message("  [ok]   pooled h2 reproduces heterogeneity.tsv h2_pooled exactly")
  }
  w(out, "snp_h2_pooled.csv")
})
# ===========================================================================
# Gender simple slopes - decomposition of the prespecified PGI x gender
# interaction
# ===========================================================================
# The pooled PGI-Education slope WITHIN women and WITHIN men, collapsing region,
# for both outcomes. This is a linear combination of coefficients already in the
# frozen A4 / B4 models -- the PGI main effect plus the PGI x gender term at each
# level of gender -- with the SE from the model's own coefficient covariance. It
# is NOT a fresh fit and NOT a re-derived sample: the estimates come from the
# same fitted objects as the interaction the Results already report, so the three
# numbers cohere and a reader who adds them up gets the interaction back.
#
# gender_c is coded +0.5 female / -0.5 male, so the interaction is the change in
# slope per unit of gender_c and the simple slopes are main effect +/- 0.5 x
# interaction. The female-minus-male contrast is therefore IDENTICAL to the
# interaction row, in estimate, SE and p -- which is what the recovery assertion
# below checks. It runs on every export and stops rather than warns.
#
# PRESPECIFIED. The interaction is Analysis C of the plan; these slopes
# reparameterise that same fitted model at fixed levels of a contrast-coded
# moderator. There is exactly one such combination and no analytic choice was
# made after seeing data, so the inferential status is the interaction's: this
# is a reparameterisation of the Analysis C interaction, not a new test.
try_step("gender simple slopes (decomposition of the prespecified PGI x gender interaction)", {
  MAIN <- "PGI_Edu_z"; INT <- "PGI_Edu_z:gender_c"
  GENDERS <- c(female = 0.5, male = -0.5)   # gender_c coding, 00_setup/prepare_samples.R

  # Contrast on the fitted model: estimate, SE from vcov, CI and p.
  # broom::tidy.gam(parametric = TRUE, conf.int = TRUE) builds the frozen
  # coefficient tables with a qnorm-based CI; the p-value comes from
  # summary.gam()'s p.table, which uses a t with the model's residual df. Rather
  # than assume which, both are computed and the one that reproduces the frozen
  # interaction p is used -- and if neither does, the step stops.
  contrast_row <- function(mod, k, rdf) {
    b <- coef(mod); V <- vcov(mod)
    est <- sum(k * b)
    se  <- sqrt(as.numeric(t(k) %*% V %*% k))
    tst <- est / se
    list(estimate = est, se = se, statistic = tst,
         p_norm = 2 * stats::pnorm(-abs(tst)),
         p_t    = 2 * stats::pt(-abs(tst), rdf),
         ci_lo  = est - stats::qnorm(.975) * se,
         ci_hi  = est + stats::qnorm(.975) * se)
  }
  kvec <- function(mod, gc) {
    k <- setNames(rep(0, length(coef(mod))), names(coef(mod)))
    stopifnot(all(c(MAIN, INT) %in% names(k)))
    k[MAIN] <- 1; k[INT] <- gc; k
  }

  slopes_for <- function(mod, frozen, outcome, model_label) {
    rdf <- tryCatch(summary(mod)$residual.df, error = function(e) stats::df.residual(mod))
    n   <- tryCatch(stats::nobs(mod), error = function(e) NA_integer_)
    fz  <- frozen[frozen$term == INT, ]
    if (nrow(fz) != 1) stop("frozen ", INT, " row not found for ", model_label)

    # --- which p-value convention reproduces the frozen interaction? ---------
    dif <- contrast_row(mod, kvec(mod, 0.5) - kvec(mod, -0.5), rdf)   # = the INT row
    use_t <- abs(dif$p_t - fz$p.value) <= abs(dif$p_norm - fz$p.value)
    pick  <- function(r) if (use_t) r$p_t else r$p_norm

    # --- RECOVERY ASSERTION -------------------------------------------------
    # The female-minus-male contrast IS the interaction coefficient. If the
    # linear combination were formed wrongly this is where it shows.
    for (chk in list(list("estimate", dif$estimate, fz$estimate, 1e-10),
                     list("SE",       dif$se,       fz$std.error, 1e-10),
                     list("p",        pick(dif),    fz$p.value,   1e-8))) {
      if (!isTRUE(all.equal(chk[[2]], chk[[3]], tolerance = chk[[4]])))
        stop(sprintf("%s: female-minus-male %s = %.15g does not recover the frozen interaction %s = %.15g",
                     model_label, chk[[1]], chk[[2]], chk[[1]], chk[[3]]))
    }
    message(sprintf("  [ok]   %s: women - men recovers the frozen %s exactly (b = %+.6f, SE = %.6f, p = %.6f; %s p-values)",
                    model_label, INT, dif$estimate, dif$se, pick(dif),
                    if (use_t) sprintf("t(%.0f)", rdf) else "normal"))

    rows <- lapply(names(GENDERS), function(g) {
      r <- contrast_row(mod, kvec(mod, GENDERS[[g]]), rdf)
      data.frame(outcome = outcome, model = model_label, term = paste0("PGI-Education slope within ", g),
                 gender = g, gender_c = GENDERS[[g]], n = n,
                 estimate = r$estimate, se = r$se, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
                 statistic = r$statistic, p = pick(r), stringsAsFactors = FALSE)
    })
    # The contrast is carried in the file so the recovery is visible to a reader,
    # not only asserted at export time.
    rows[[length(rows) + 1L]] <- data.frame(
      outcome = outcome, model = model_label,
      term = "difference (women - men) = PGI_Edu_z:gender_c",
      gender = "difference", gender_c = NA_real_, n = n,
      estimate = dif$estimate, se = dif$se, ci_lo = dif$ci_lo, ci_hi = dif$ci_hi,
      statistic = dif$statistic, p = pick(dif), stringsAsFactors = FALSE)
    do.call(rbind, rows)
  }

  ss <- rbind(
    slopes_for(A$m_A4_re, A$A4_Coefficients, "Educational attainment", "A4_RE"),
    slopes_for(B$m_B4_re, B$B4_Coefficients, "Intergenerational educational mobility", "B4_RE"))

  # --- soft check: the pooled slopes must sit inside the per-cell slopes -----
  # fig04_percell_slopes.csv is a separate per-cell model (A0g/B0g), so this is
  # a sanity bracket on the sign convention, not an identity.
  pc <- tryCatch(utils::read.csv(file.path(EXP_DIR, "fig04_percell_slopes.csv"),
                                 stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(pc)) {
    map <- c("Educational attainment" = "Years of education",
             "Intergenerational educational mobility" = "Educational mobility")
    for (i in which(ss$gender %in% c("female", "male"))) {
      s <- pc$slope[pc$outcome == map[[ss$outcome[i]]] & pc$gender == ss$gender[i]]
      if (!length(s)) next
      inside <- ss$estimate[i] > min(s) && ss$estimate[i] < max(s)
      message(sprintf("  [%s] %s %s: pooled %.4f vs per-cell [%.4f, %.4f]",
                      if (inside) "ok" else "WARN", ss$model[i], ss$gender[i],
                      ss$estimate[i], min(s), max(s)))
    }
  }

  ss$scale  <- "SD of birth-year-standardised outcome per SD of PGI-Education"
  ss$label  <- "prespecified (decomposes the Analysis C PGI x gender interaction)"
  ss$caveat <- paste(
    "A linear combination of coefficients from the SAME frozen model as the",
    "interaction the Results report (A4 attainment, B4 mobility) -- not a new fit,",
    "not a new sample, and not an independent test. The two slopes differ by exactly",
    "the interaction coefficient, so they carry no evidence beyond it. Region is",
    "collapsed at the unweighted midpoint of the +/-0.5 contrast coding, not the",
    "sample-weighted regional mean.")
  w(ss, "gender_pgi_simple_slopes.csv")
})

try_step("focal-estimate table", {
  # Consolidated focal terms as a table.
  # Linear/step focals carry an estimate + CI; GAM omnibus carries an LRT.
  # Column pickers tolerant to broom vs custom naming.
  pick <- function(df, cands) { m <- intersect(cands, names(df)); if (length(m)) df[[m[1]]] else NULL }
  grab <- function(coef_df, term, id, label) {
    if (is.null(coef_df) || !("term" %in% names(coef_df))) return(NULL)
    r <- coef_df[coef_df$term == term, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    g <- function(c) { v <- pick(r, c); if (is.null(v)) NA_real_ else v[1] }
    data.frame(analysis = id, label = label, kind = "coef",
               estimate = g(c("estimate", "Estimate")),
               ci_lo    = g(c("conf.low", "ci_lo", "CI lower")),
               ci_hi    = g(c("conf.high", "ci_hi", "CI upper")),
               p        = g(c("p.value", "p_value", "p")))
  }
  lrt <- function(lrt_df, id, label) {
    if (is.null(lrt_df)) return(NULL)
    # omnibus = the "linear vs nonlinear" row of nested_lrt_table
    row <- if ("comparison" %in% names(lrt_df))
      lrt_df[grepl("nonlinear", lrt_df$comparison), , drop = FALSE] else lrt_df
    if (!nrow(row)) row <- lrt_df[nrow(lrt_df), , drop = FALSE]
    g <- function(c) { v <- pick(row, c); if (is.null(v)) NA_real_ else v[1] }
    data.frame(analysis = id, label = label, kind = "lrt",
               estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
               p = g(c("LRT_p", "lrt_p", "p.value", "p_value", "p")))
  }
  forest <- dplyr::bind_rows(
    grab(A$A1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c", "A1", "Attainment: PGI x BY x Region"),
    grab(A$A2_Coefficients, "PGI_Edu_z:reunifpost:east_west_c", "A2", "Attainment: PGI x Post x Region"),
    lrt(A$A3_NestedLRT, "A3", "Attainment: nonlinearity (omnibus)"),
    grab(B$B1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c", "B1", "Mobility: PGI x BY x Region"),
    grab(B$B2_Coefficients, "PGI_Edu_z:reunifpost:east_west_c", "B2", "Mobility: PGI x Post x Region"),
    lrt(B$B3_NestedLRT, "B3", "Mobility: nonlinearity (omnibus)"))
  w(forest, "focal_estimates.csv")
})

# ===========================================================================
# Provenance sidecar + supplementary-workbook transfer manifest
# ===========================================================================
# Cheap cache fingerprint: size + mtime. Full SHA-256 of the caches is avoided
# on purpose — they are ~0.6-0.7 GB each and live on the SMB share, so hashing
# them would read gigabytes over the network. import_results.R still SHA-256s
# each imported CSV in the manuscript repo; this sidecar only needs to pin which
# run/commit produced them.
finger <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  info <- file.info(path)
  sprintf("size=%d;mtime=%s", info$size, format(info$mtime, "%Y-%m-%dT%H:%M:%S"))
}
git1 <- function(...) tryCatch(
  system2("git", c("-C", shQuote(WORKFLOW_DIR), ...), stdout = TRUE)[1],
  error = function(e) NA_character_)

try_step("provenance sidecar", {
  prov <- data.frame(
    run_ts               = RUN_TS,
    exported             = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    analysis_repo_commit = git1("rev-parse", "HEAD"),
    analysis_repo_branch = git1("rev-parse", "--abbrev-ref", "HEAD"),
    a_cache_fingerprint  = finger(A_PATH),
    b_cache_fingerprint  = finger(B_PATH),
    desc_cache_fingerprint = finger(D_PATH),
    stringsAsFactors = FALSE)
  # write directly (single-row provenance is aggregate metadata, not data)
  utils::write.csv(prov, file.path(EXP_DIR, "_provenance.csv"), row.names = FALSE)
  message("  [ok]   _provenance.csv")
})

try_step("supplementary-data transfer manifest", {
  wbs <- data.frame(
    supplement_id = paste0("S", 1:7),
    source_file = c(
      "00_Descriptives/Supplementary_Data_1_Descriptives_results.xlsx",
      "A_Education/Supplementary_Data_2_Education_results.xlsx",
      "A_Education/Supplementary_Data_3_Education_sensitivities_results.xlsx",
      "A_Education/Supplementary_Data_4_Education_appendix_results.xlsx",
      "B_Mobility/Supplementary_Data_5_Mobility_results.xlsx",
      "B_Mobility/Supplementary_Data_6_Mobility_sensitivities_results.xlsx",
      "B_Mobility/Supplementary_Data_7_Mobility_appendix_results.xlsx"),
    stringsAsFactors = FALSE)
  lines <- c(
    "# Supplementary-data transfer manifest",
    "",
    sprintf("Run: RUN_%s. Exported: %s.", RUN_TS, format(Sys.time(), "%Y-%m-%d")),
    "",
    "The seven Supplementary Data workbooks are publication-ready and transfer",
    "**directly** to the manuscript repo's `03_supplement/` folder (renamed S1-S7).",
    "They are NOT imported through `import_results.R` (that channel targets",
    "`01_results/` and is for the tidy figure-source CSVs in this same folder).",
    "",
    "| Supplement | Source (relative to run folder) |",
    "|---|---|",
    paste0("| ", wbs$supplement_id, " | `", wbs$source_file, "` |"),
    "",
    "Tidy figure-source CSVs in this folder -> `01_results/` via:",
    "```",
    "Rscript import_results.R --from <run_folder> \\",
    "  manuscript_export/fig01_kernel_curves.csv \\",
    "  manuscript_export/fig02_attainment_slope_curves.csv  [...]",
    "```")
  writeLines(lines, file.path(EXP_DIR, "SUPPLEMENTARY_DATA_TRANSFER.md"))
  utils::write.csv(wbs, file.path(EXP_DIR, "supplementary_data_manifest.csv"),
                   row.names = FALSE)
  message("  [ok]   SUPPLEMENTARY_DATA_TRANSFER.md + supplementary_data_manifest.csv")
})

message(sprintf("\nDone. %d figure-source CSV(s) + metadata written to:\n  %s",
                length(written), EXP_DIR))
