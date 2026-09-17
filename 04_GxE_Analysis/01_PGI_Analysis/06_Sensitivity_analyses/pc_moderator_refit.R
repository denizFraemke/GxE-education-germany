# pc_moderator_refit.R
# Sensitivity of the four published focal interactions to ancestry adjustment
# entering the outcome model.
#
# The published models residualise each PGI on the pooled cross-cohort PCs
# before z-scoring it (00_setup/prepare_samples.R), which removes PC variance
# from the predictor but leaves PC main effects on the outcome and PC x moderator
# associations. A focal PGI x moderator interaction can absorb either.
#
# Three rungs, one factor apart:
#
#   M0_published       as published
#   pc_main          + crossPC1-10 main effects in the outcome model
#   pc_x_moderators  + crossPC1-10 main effects + PC x {BYc_z, east_west_c, gender_c}
#
# pc_main is not optional: a PC x moderator interaction cannot be fitted without
# PC main effects, so the third rung against M0 alone would bundle two changes.
# Within-study residualisation scope is covered by ancestry_sensitivity.R.
#
# Each focal is refitted on the model it is reported from, with that model's rhs
# verbatim from the published script, so M0's term names match the frozen
# coefficient tables:
#
#   PGI-Education x birth year, attainment   A1   A_Education/A_analysis.R:85,236
#   PGI-Education x region, mobility         B0   B_Mobility/B_analysis.R:134-139
#   PGI-Education x gender, attainment       A4   A_Education/A_analysis.R:339,346
#   PGI-Education x gender, mobility         B4   B_Mobility/B_analysis.R:349-352
#
# The gender focals use A4/B4 rather than the per-cell slope models A0g/B0g,
# which give four slopes and no test of the difference. The run stops unless M0
# reproduces the frozen estimate for every focal.
#
# Every model is linear in its parameters; the only smooth is the family random
# intercept the published focals already carry. The GAM ladder is not touched.
# The rungs add coefficients without adding observations, so intervals are
# expected to widen; each row carries its SE, CI, coefficient count and the two
# counterfactual p-values computed below.
#
# Nothing here overwrites the frozen run and its model caches are not loaded
# (~650 MB each). The analytic samples are rebuilt once from the merged data
# through 00_setup/, as 05_manuscript_export/check_frozen_run.R step 4 does.
# Aggregate output only; no individual row is printed or written.
#
# Usage:
#   RUN_TS=<ts> Rscript 05_manuscript_export/check_frozen_run.R
#   RUN_TS=<ts> Rscript 06_Sensitivity_analyses/pc_moderator_refit.R
#
#   OUT_DIR=<dir>            write the CSVs elsewhere (default: the run's
#                            manuscript_export/sensitivity_analyses/)
#   MANUSCRIPT_REPO=<path>   where 01_results/ lives
#   DATA_ROOT=<path>         data root

suppressPackageStartupMessages({
  library(dplyr); library(mgcv)
})

# ---------------------------------------------------------------------------
# Locations
# ---------------------------------------------------------------------------
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine pc_moderator_refit.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
REPO_ROOT    <- normalizePath(file.path(WORKFLOW_DIR, "..", ".."), mustWork = TRUE)
rm(.find_this_file)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
FROZEN_RUN_TS <- "20260729_1555"
Z975          <- stats::qnorm(0.975)

RE_TERM   <- "+ s(fid_re, bs = 're')"          # verbatim from A_/B_analysis.R
CROSS_PCS <- paste0("crossPC", 1:10)           # the canonical shared PC basis
MODS      <- c("BYc_z", "east_west_c", "gender_c")   # spec_bases()$mods, verbatim

# Reproduction tolerances for the M0 anchor. The samples are rebuilt from the
# same merge through the same 00_setup/ code, so the slack is for BLAS and mgcv
# numerical drift only.
TOL_EST <- 1e-6
TOL_SE  <- 1e-6
TOL_P   <- 1e-5

# The four bases. Each rhs and focal list is verbatim from the published
# analysis script.
BASES <- list(
  A1 = list(
    sample   = "attainment", outcome = "edu_z_kernel",
    rhs      = "edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3",
    focals   = c("PGI_Edu_z", "PGI_Edu_z:BYc_z", "PGI_Edu_z:east_west_c",
                 "PGI_Edu_z:gender_c", "PGI_Edu_z:BYc_z:east_west_c"),
    headline = "PGI_Edu_z:BYc_z",
    finding  = "PGI-Education x birth year (educational attainment)",
    frozen_file  = file.path("attainment", "main", "A1_Coefficients.csv"),
    frozen_model = "A1_RE",
    note     = "A1 as published: three-way linear attainment model (A_analysis.R:85)"),
  B0 = list(
    sample   = "mobility", outcome = "mobility",
    rhs      = paste("mobility ~ PGI_Edu_z * east_west_c +",
                     "east_west_c * (gender_c + BYc_z + parental_edu_z_kernel) +",
                     "gender_c * BYc_z + gender_c * parental_edu_z_kernel +",
                     "BYc_z * parental_edu_z_kernel"),
    focals   = c("PGI_Edu_z", "PGI_Edu_z:east_west_c"),
    headline = "PGI_Edu_z:east_west_c",
    finding  = "PGI-Education x region (educational mobility)",
    frozen_file  = file.path("mobility", "main", "B0_Coefficients.csv"),
    frozen_model = "B0_RE",
    note     = "B0 as published: overall PGI x region mobility contrast (B_analysis.R:134)"),
  A4 = list(
    sample   = "attainment", outcome = "edu_z_kernel",
    rhs      = "edu_z_kernel ~ PGI_Edu_z * BYc_z * east_west_c * gender_c",
    focals   = c("PGI_Edu_z:gender_c", "PGI_Edu_z:BYc_z:gender_c",
                 "PGI_Edu_z:east_west_c:gender_c",
                 "PGI_Edu_z:BYc_z:east_west_c:gender_c"),
    headline = "PGI_Edu_z:gender_c",
    finding  = "PGI-Education x gender (educational attainment)",
    frozen_file  = file.path("attainment", "main", "A4_Coefficients.csv"),
    frozen_model = "A4_RE",
    note     = "A4 as published: four-way linear attainment model (A_analysis.R:339)"),
  B4 = list(
    sample   = "mobility", outcome = "mobility",
    rhs      = paste("mobility ~ PGI_Edu_z * BYc_z * east_west_c * gender_c +",
                     "parental_edu_z_kernel * BYc_z * east_west_c * gender_c"),
    focals   = c("PGI_Edu_z:gender_c", "PGI_Edu_z:BYc_z:gender_c",
                 "PGI_Edu_z:BYc_z:east_west_c:gender_c"),
    headline = "PGI_Edu_z:gender_c",
    finding  = "PGI-Education x gender (educational mobility)",
    frozen_file  = file.path("mobility", "main", "B4_Coefficients.csv"),
    frozen_model = "B4_RE",
    note     = "B4 as published: four-way linear mobility model (B_analysis.R:349)")
)

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
source_commit <- function() {
  out <- suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), "rev-parse", "HEAD"),
                                  stdout = TRUE, stderr = FALSE))
  if (length(out) == 1L && nzchar(out)) out else NA_character_
}
source_dirty <- function() {
  out <- suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), "status", "--porcelain",
                                           "--", shQuote(SECTION_DIR)),
                                  stdout = TRUE, stderr = FALSE))
  length(out) > 0
}

# Checks collector, verbatim from ancestry_sensitivity.R. kind = "check" gates the
# run; kind = "diagnostic" records an aggregate for the record and never fails.
new_checks <- function() {
  e <- new.env(); e$rows <- list()
  push <- function(kind, name, expected, observed, tol, note) {
    d <- abs(expected - observed)
    e$rows[[length(e$rows) + 1]] <- data.frame(
      kind = kind, check = name, expected = expected, observed = observed,
      abs_diff = d, tol = tol,
      pass = if (identical(kind, "diagnostic")) NA else (is.finite(d) && d <= tol),
      note = note, stringsAsFactors = FALSE)
    invisible(NULL)
  }
  e$add  <- function(name, expected, observed, tol, note = "")
    push("check", name, expected, observed, tol, note)
  e$diag <- function(name, observed, note = "")
    push("diagnostic", name, NA_real_, observed, NA_real_, note)
  e$table  <- function() do.call(rbind, e$rows)
  e$failed <- function() { t <- e$table(); any(!is.na(t$pass) & !t$pass) }
  e
}

# ---------------------------------------------------------------------------
# The ladder
# ---------------------------------------------------------------------------
# `pc_main` is the PC main-effect block and `pcx` crosses that whole block with
# each moderator, so the third rung adds 10 main effects and 10 x length(mods)
# interactions to the published rhs.
# `pc_terms_added` is cumulative against M0; `wald_*` describe the block this
# rung adds against the rung above it, which is what the joint test is about.
# R labels an interaction by the order its variables first appear in the
# formula, and the moderators appear in the published rhs before the PC block,
# so the terms come back as `BYc_z:crossPC1`, not `crossPC1:BYc_z`. Both orders
# are matched so the test cannot silently select nothing if that changes, and
# `wald_n` is checked against the number of terms matched.
build_rungs <- function(base, pc_main = paste(CROSS_PCS, collapse = " + "), mods = MODS) {
  pcx  <- paste(sprintf("(%s):%s", pc_main, mods), collapse = " + ")
  mpat <- paste(mods, collapse = "|")
  list(
    list(id = "M0_published", order = 1L, cluster = "-", rhs = base$rhs,
         pc_terms_added = 0L, wald_pattern = NA_character_, wald_n = 0L,
         wald_block = NA_character_,
         desc = "as published (PGI pre-residualised on pooled crossPC1-10)"),
    list(id = "pc_main", order = 2L, cluster = "A",
         rhs = paste(base$rhs, "+", pc_main),
         pc_terms_added = length(CROSS_PCS),
         wald_pattern = "^crossPC[0-9]+$",
         wald_n = length(CROSS_PCS),
         wald_block = "crossPC1-10 main effects",
         desc = "+ crossPC1-10 main effects in the outcome model"),
    list(id = "pc_x_moderators", order = 3L, cluster = "A",
         rhs = paste(base$rhs, "+", pc_main, "+", pcx),
         pc_terms_added = length(CROSS_PCS) * (1L + length(mods)),
         wald_pattern = sprintf("^((%s):crossPC[0-9]+|crossPC[0-9]+:(%s))$", mpat, mpat),
         wald_n = length(CROSS_PCS) * length(mods),
         wald_block = sprintf("crossPC1-10 x {%s}", paste(mods, collapse = ", ")),
         desc = sprintf("+ crossPC1-10 main effects + PC x {%s}", paste(mods, collapse = ", ")))
  )
}

# Fit one rung. The RE term is appended exactly as A_/B_analysis.R append it, and
# the engine (discrete or the discrete = FALSE fallback) is recorded on every row
# so a fallback is never silent.
fit_rung <- function(dat, rhs, focals) {
  d <- dat
  d$fid_re <- droplevels(factor(d$fid_re))
  fml <- stats::as.formula(paste(rhs, RE_TERM))

  engine <- "bam fREML discrete"
  m <- tryCatch(mgcv::bam(fml, data = d, method = "fREML", discrete = TRUE),
                error = function(e) NULL)
  if (is.null(m)) {
    engine <- "bam fREML (discrete=FALSE fallback)"
    m <- tryCatch(mgcv::bam(fml, data = d, method = "fREML", discrete = FALSE),
                  error = function(e) e)
  }
  if (is.null(m) || inherits(m, "error"))
    stop("rung failed to fit: ",
         if (inherits(m, "error")) conditionMessage(m) else "fit returned NULL")

  sm   <- summary(m)
  cs   <- sm$p.table
  pcol <- if ("Pr(>|t|)" %in% colnames(cs)) "Pr(>|t|)" else colnames(cs)[ncol(cs)]
  n    <- as.integer(sm$n)
  rows <- lapply(focals, function(tm) {
    if (!tm %in% rownames(cs))
      stop("focal term absent from the fit: ", tm)
    est <- unname(cs[tm, "Estimate"]); se <- unname(cs[tm, "Std. Error"])
    data.frame(term = tm, estimate = est, se = se,
               ci_lo = est - Z975 * se, ci_hi = est + Z975 * se,
               p = unname(cs[tm, pcol]), n = n, n_coef = nrow(cs),
               engine = engine, stringsAsFactors = FALSE)
  })
  list(coefs = dplyr::bind_rows(rows), mod = m, n = n, n_coef = nrow(cs))
}

# Joint Wald chi-square that all coefficients matching `pattern` are zero. A Wald test and
# not a likelihood-ratio test because fREML log-likelihoods are not comparable
# across models with different fixed effects.
wald_joint <- function(mod, pattern) {
  b <- stats::coef(mod); V <- stats::vcov(mod)
  idx <- grep(pattern, names(b))
  if (!length(idx)) return(list(df = 0L, chisq = NA_real_, p = NA_real_, n_terms = 0L))
  Vi <- tryCatch(solve(V[idx, idx, drop = FALSE]), error = function(e) NULL)
  if (is.null(Vi)) return(list(df = length(idx), chisq = NA_real_, p = NA_real_,
                               n_terms = length(idx)))
  chi <- as.numeric(t(b[idx]) %*% Vi %*% b[idx])
  list(df = length(idx), chisq = chi,
       p = stats::pchisq(chi, length(idx), lower.tail = FALSE), n_terms = length(idx))
}

# ---------------------------------------------------------------------------
# Attenuation vs precision loss
# ---------------------------------------------------------------------------
# Two counterfactuals: hold the SE at its M0 value and let only the estimate
# move, then hold the estimate at M0 and let only the SE move. Both are on the
# normal scale, the scale the reported CIs use, so they are comparable to each
# other rather than to bam's own p-value, which they bracket.
p_normal <- function(est, se) 2 * stats::pnorm(-abs(est / se))

# `pct_change_estimate_vs_M0` is signed against M0's own sign, so a POSITIVE
# value always means the coefficient moved AWAY from zero, whatever its sign.
# The driver label names the channel and the direction of the change.
decompose <- function(est, se, p, est0, se0, p0, floor_pct = 5) {
  d_est <- 100 * (est  - est0) / est0
  d_se  <- 100 * (se   - se0)  / se0
  z     <- abs(est / se); z0 <- abs(est0 / se0)
  driver <- if (!is.finite(d_est) || !is.finite(d_se)) NA_character_
    else if (abs(d_est) < floor_pct && abs(d_se) < floor_pct) "no material change"
    else if (abs(d_est) >= 2 * abs(d_se))
      sprintf("estimate-driven (%s)",
              if (abs(est) < abs(est0)) "attenuated" else "amplified")
    else if (abs(d_se)  >= 2 * abs(d_est))
      sprintf("precision-driven (%s)", if (se > se0) "SE widened" else "SE narrowed")
    else "both"
  data.frame(
    estimate_M0 = est0, se_M0 = se0, p_M0 = p0,
    delta_estimate_vs_M0    = est - est0,
    pct_change_estimate_vs_M0 = d_est,
    se_ratio_vs_M0          = se / se0,
    pct_change_se_vs_M0     = d_se,
    abs_z = z, abs_z_M0 = z0,
    pct_change_abs_z_vs_M0  = 100 * (z - z0) / z0,
    p_if_only_estimate_moved = p_normal(est,  se0),
    p_if_only_se_moved       = p_normal(est0, se),
    driver = driver,
    sign_flip = is.finite(est) & is.finite(est0) & sign(est) != sign(est0),
    crosses_05_boundary = is.finite(p) & is.finite(p0) & ((p < .05) != (p0 < .05)),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Samples: rebuilt once from the merge through 00_setup/ (check_frozen_run.R step 4)
# ---------------------------------------------------------------------------
rebuild_samples <- function() {
  suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(tibble) })
  RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()
  EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
  source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
  source(file.path(WORKFLOW_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
  DATA_DIR <- file.path(Sys.getenv("DATA_ROOT"), "03_Merge")
  sink(tempfile())                                   # load/prepare are chatty by design
  on.exit(sink(), add = TRUE)
  suppressMessages(suppressWarnings({
    source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
    source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)
  }))
  sink(); on.exit()

  need_A <- c("cohort", "east_west", "east_west_c", "gender", "gender_c", "birth_year",
              "BYc_z", "PGI_Edu_z", "fid_re", "edu_z_kernel", CROSS_PCS)
  need_B <- c(need_A, "parental_edu_z_kernel", "mobility")
  datA <- as.data.frame(dat_cluster)[, need_A]
  datB <- as.data.frame(dat_cluster_mob)[, need_B]
  rm(list = intersect(c("dat_cluster", "dat_cluster_mob", "df_raw", "df", "df_with_twins",
                        "dat_edu", "dat_mob"), ls()))
  invisible(gc())
  list(attainment = datA, mobility = datB, data_file = DATA_FILE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_pc_moderator_refit <- function(out_dir, results_dir, samples, run_ts) {
  run_ts <- sub("^RUN_", "", run_ts)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  commit <- source_commit(); dirty <- source_dirty()
  checks <- new_checks()
  t0 <- Sys.time()
  message(sprintf("A5 PC x moderator ladder: RUN_%s\n  exports  -> %s\n  commit      %s%s",
                  run_ts, out_dir, substr(commit, 1, 8),
                  if (dirty) " (working tree dirty)" else ""))

  prov <- utils::read.csv(file.path(results_dir, "_provenance.csv"), stringsAsFactors = FALSE)
  frozen_commit <- as.character(prov$analysis_repo_commit[1])

  # ---- PC completeness: the ladder must not change the sample -------------
  # prepare_samples.R states that every PGI-present individual has complete
  # crossPC, so adding PC terms drops no one. Verified rather than assumed: if it
  # failed, a moved coefficient could be a sample change rather than an
  # adjustment effect, and the ladder would not be a ladder.
  for (sm in c("attainment", "mobility")) {
    d <- samples[[sm]]
    checks$add(sprintf("%s: crossPC1-10 complete on every row (ladder cannot change the sample)", sm),
               nrow(d), sum(stats::complete.cases(d[, CROSS_PCS, drop = FALSE])), 0,
               note = "0 tolerance: any missing PC would make A1/A2 a different sample from M0")
  }

  ladder <- list(); wald <- list()
  for (bn in names(BASES)) {
    base <- BASES[[bn]]
    dat  <- samples[[base$sample]]
    message(sprintf("\nBase %s (%s) - %s", bn, base$sample, base$finding))

    frozen <- utils::read.csv(file.path(results_dir, base$frozen_file), stringsAsFactors = FALSE)
    frozen <- frozen[frozen$model == base$frozen_model, , drop = FALSE]
    if (!nrow(frozen)) stop("frozen model rows not found: ", base$frozen_model)

    m0 <- NULL
    for (rg in build_rungs(base)) {
      message(sprintf("  - %-20s", rg$id), appendLF = FALSE)
      r <- fit_rung(dat, rg$rhs, base$focals)
      message(sprintf(" n = %5d, %2d coefficients", r$n, r$n_coef))

      # M0 is the anchor: it must reproduce the frozen published estimate for
      # every focal of this base, or the whole base is meaningless.
      if (rg$id == "M0_published") {
        m0 <- r$coefs
        for (tm in base$focals) {
          fr <- frozen[frozen$term == tm, , drop = FALSE]
          if (!nrow(fr)) stop("frozen row not found: ", base$frozen_model, " / ", tm)
          rr <- r$coefs[r$coefs$term == tm, ]
          tag <- sprintf("M0 reproduces frozen %s %s", base$frozen_model, tm)
          checks$add(paste(tag, "- estimate"), fr$estimate[1], rr$estimate, TOL_EST,
                     note = base$frozen_file)
          checks$add(paste(tag, "- SE"),       fr$std.error[1], rr$se,       TOL_SE,
                     note = base$frozen_file)
          checks$add(paste(tag, "- p"),        fr$p.value[1],   rr$p,        TOL_P,
                     note = base$frozen_file)
        }
        checks$diag(sprintf("%s: published model coefficient count", bn), r$n_coef,
                    note = sprintf("%s, n = %d, %.1f observations per coefficient",
                                   base$outcome, r$n, r$n / r$n_coef))
      } else {
        # Every rung must be fitted on the same rows as M0, and must add exactly
        # the coefficients the rung says it adds.
        checks$add(sprintf("%s / %s: n unchanged from M0", bn, rg$id),
                   m0$n[1], r$n, 0,
                   note = "the rungs differ by adjustment only, never by sample")
        checks$add(sprintf("%s / %s: adds exactly %d PC coefficients", bn, rg$id,
                           rg$pc_terms_added),
                   m0$n_coef[1] + rg$pc_terms_added, r$n_coef, 0,
                   note = sprintf("published %d + %d PC terms", m0$n_coef[1], rg$pc_terms_added))

        # Does the added PC block carry any signal at all?
        # Gated on matching the whole block: a pattern that silently matches
        # nothing would report df = 0 and look like an absent effect.
        w <- wald_joint(r$mod, rg$wald_pattern)
        checks$add(sprintf("%s / %s: joint Wald matches all %d coefficients of the added block",
                           bn, rg$id, rg$wald_n),
                   rg$wald_n, w$n_terms, 0,
                   note = sprintf("pattern %s", rg$wald_pattern))
        wald[[paste(bn, rg$id)]] <- data.frame(
          base = bn, sample = base$sample, rung = rg$id, block = rg$wald_block,
          df = w$df, chisq = w$chisq, p = w$p, stringsAsFactors = FALSE)
        checks$diag(sprintf("%s / %s: joint Wald chi-square on the added PC block", bn, rg$id),
                    w$chisq, note = sprintf("%s, df = %d", rg$wald_block, w$df))
        checks$diag(sprintf("%s / %s: joint Wald p on the added PC block", bn, rg$id),
                    w$p, note = "H0: every coefficient of the block this rung adds is zero")
      }

      dec <- do.call(rbind, lapply(base$focals, function(tm) {
        rr <- r$coefs[r$coefs$term == tm, ]
        ss <- m0[m0$term == tm, ]
        decompose(rr$estimate, rr$se, rr$p, ss$estimate, ss$se, ss$p)
      }))
      ladder[[paste(bn, rg$id)]] <- cbind(
        data.frame(
          finding   = ifelse(base$focals == base$headline, base$finding, ""),
          base      = bn, sample = base$sample, outcome = base$outcome,
          rung      = rg$id, rung_order = rg$order, cluster = rg$cluster,
          rung_desc = rg$desc, pc_terms_added = rg$pc_terms_added,
          is_headline_focal = base$focals == base$headline,
          stringsAsFactors = FALSE),
        r$coefs,
        data.frame(obs_per_coef = round(r$coefs$n / r$coefs$n_coef, 1),
                   stringsAsFactors = FALSE),
        dec)
      rm(r); invisible(gc())
    }
  }

  tab <- dplyr::bind_rows(ladder)
  num <- c("estimate", "se", "ci_lo", "ci_hi", "estimate_M0", "se_M0",
           "delta_estimate_vs_M0", "se_ratio_vs_M0", "abs_z", "abs_z_M0")
  tab[num] <- lapply(tab[num], function(x) round(x, 6))
  pct <- c("pct_change_estimate_vs_M0", "pct_change_se_vs_M0", "pct_change_abs_z_vs_M0")
  tab[pct] <- lapply(tab[pct], function(x) round(x, 1))
  tab <- tab[, c("finding", "base", "sample", "outcome", "term", "is_headline_focal",
                 "rung", "rung_order", "cluster", "rung_desc", "pc_terms_added",
                 "estimate", "se", "ci_lo", "ci_hi", "p", "n", "n_coef", "obs_per_coef",
                 "estimate_M0", "se_M0", "p_M0",
                 "delta_estimate_vs_M0", "pct_change_estimate_vs_M0",
                 "se_ratio_vs_M0", "pct_change_se_vs_M0",
                 "abs_z", "abs_z_M0", "pct_change_abs_z_vs_M0",
                 "p_if_only_estimate_moved", "p_if_only_se_moved", "driver",
                 "sign_flip", "crosses_05_boundary", "engine")]
  tab$ci_rule <- "normal, qnorm(0.975) x SE (as the frozen exports write it)"
  tab$anchor  <- ifelse(tab$rung == "M0_published",
                        "verified against the frozen published coefficient table (pc_moderator_checks.csv)", "")

  # ---- write -------------------------------------------------------------
  w <- function(df, name) {
    utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
    message(sprintf("  [ok]   %-34s %d x %d", name, nrow(df), ncol(df)))
  }
  w(tab, "pc_moderator_ladder.csv")
  chk <- checks$table(); chk$source_commit <- commit; chk$run_ts <- run_ts
  w(chk, "pc_moderator_checks.csv")
  provo <- data.frame(
    run_ts = run_ts, frozen_analysis_repo_commit = frozen_commit, script_commit = commit,
    script_working_tree_dirty = dirty, reference_source = results_dir,
    data_source = samples$data_file,
    n_attainment = nrow(samples$attainment), n_mobility = nrow(samples$mobility),
    bases_fitted = paste(names(BASES), collapse = ", "),
    rungs_per_base = 3L, models_fitted = length(BASES) * 3L,
    focal_rows = nrow(tab),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R_version = R.version.string,
    mgcv_version = as.character(utils::packageVersion("mgcv")),
    checks_total = sum(chk$kind == "check"),
    checks_passed = sum(chk$kind == "check" & chk$pass),
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    stringsAsFactors = FALSE)
  w(provo, "pc_moderator_provenance.csv")

  # ---- console summary (aggregates only) ---------------------------------
  cat("\n== M0 anchor ==\n")
  anch <- chk[chk$kind == "check" & grepl("^M0 reproduces", chk$check), ]
  cat(sprintf("  %d of %d frozen-reproduction checks pass\n", sum(anch$pass), nrow(anch)))
  bad <- chk[chk$kind == "check" & !chk$pass, ]
  if (nrow(bad)) for (i in seq_len(nrow(bad)))
    cat(sprintf("  [FAIL] %-80s diff %.3g (tol %.3g)\n",
                substr(bad$check[i], 1, 80), bad$abs_diff[i], bad$tol[i]))

  cat("\n== the ladder (headline focals) ==\n")
  hd <- tab[tab$is_headline_focal, ]
  for (bn in unique(hd$base)) {
    ht <- hd[hd$base == bn, ]
    cat(sprintf("\n  %s  [%s]  %s\n", bn, ht$term[1], ht$finding[1]))
    for (i in seq_len(nrow(ht)))
      cat(sprintf("    %-20s b = %+.4f  SE %.4f  95%% CI [%+.4f, %+.4f]  p = %.4f  n = %5d  k = %2d\n",
                  ht$rung[i], ht$estimate[i], ht$se[i], ht$ci_lo[i], ht$ci_hi[i],
                  ht$p[i], ht$n[i], ht$n_coef[i]))
    for (i in which(ht$rung != "M0_published"))
      cat(sprintf("      %-20s estimate %+.1f%%, SE %+.1f%%, |z| %+.1f%%  -> %s\n",
                  ht$rung[i], ht$pct_change_estimate_vs_M0[i], ht$pct_change_se_vs_M0[i],
                  ht$pct_change_abs_z_vs_M0[i], ht$driver[i]))
  }

  if (length(wald)) {
    cat("\n== joint Wald on the added PC block ==\n")
    wt <- dplyr::bind_rows(wald)
    for (i in seq_len(nrow(wt)))
      cat(sprintf("  %-3s %-20s %-38s df = %2d  chi2 = %8.2f  p = %.4g\n",
                  wt$base[i], wt$rung[i], substr(wt$block[i], 1, 38), wt$df[i],
                  wt$chisq[i], wt$p[i]))
  }

  if (checks$failed())
    cat("\nWARNING: some checks failed - read pc_moderator_checks.csv before using the exports.\n")
  cat(sprintf("\nDone in %.1f min.\n", provo$minutes))
  invisible(list(table = tab, checks = chk, provenance = provo))
}

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
if (!nzchar(Sys.getenv("RUN_TS"))) Sys.setenv(RUN_TS = FROZEN_RUN_TS)
RUN_TS <- get_run_ts()
if (!nzchar(Sys.getenv("DATA_ROOT"))) Sys.setenv(DATA_ROOT = data_root())
RUN_DIR <- file.path(gxe_dir(), "01_PGI_Analysis", "runs", paste0("RUN_", RUN_TS))
if (!dir.exists(RUN_DIR)) stop("Run folder not reachable: ", RUN_DIR)
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(RUN_DIR, "manuscript_export", "sensitivity_analyses"))
MS <- Sys.getenv("MANUSCRIPT_REPO", unset = file.path(
  path.expand("~"), "Documents", "017_manuscript", "80_years_GxE_on_Education_in_Germany-manuscript"))
RESULTS_DIR <- file.path(MS, "01_results", "manuscript_export")
if (!file.exists(file.path(RESULTS_DIR, "_provenance.csv"))) {
  message("01_results/ not found at ", MS, " - reading the frozen exports from the run folder")
  RESULTS_DIR <- file.path(RUN_DIR, "manuscript_export")
}
message("Rebuilding the analytic samples from the merge (00_setup/) ...")
SAMPLES <- rebuild_samples(); invisible(gc())
invisible(run_pc_moderator_refit(OUT_DIR, RESULTS_DIR, SAMPLES, RUN_TS))
