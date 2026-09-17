# helpers_gender_loco_altpgi.R
# ---------------------------------------------------------------------------
# Finding-matched sensitivity infrastructure for the gender results and the
# alternative-PGI substitutions. Two independent groups of helpers:
#
#   LOCO    Gender leave-one-study-out refits of the primary A4/B4/A5/B5
#           gender models: drop one contributing-study fold and refit. "LOCO"
#           is named to distinguish it from the leave-one-out (LOO) GWAS
#           weight construction upstream. These cover the pooled PGI x gender
#           interaction and the per-gender cohort GAMs (A5/B5), alongside the
#           A4_LOO/B4_LOO tables, which track only the four-way
#           PGI x BY x region x gender term.
#
#   AltPGI  Finding-matched random-effects alternative-PGI families
#           (PGI-Education, PGI-Cognition, PGI-Noncognitive) on a common
#           complete-case sample per outcome, re-residualised on ancestry PCs
#           and re-z-scored within that sample. Alongside the dedup-OLS S1
#           (which tests only the A1/B1 three-way and A4/B4 four-way), these
#           cover the reported focal findings: the pooled PGI x BY increase
#           (A1), the East-West mobility contrast (B0) and the pooled
#           PGI x gender interactions (A4/B4), plus the B5 nonlinear gender
#           pattern. Benjamini-Hochberg is applied within each focal
#           hypothesis across the three PGIs.
#
# All model fits mirror the primary specification exactly: mgcv::bam(...,
# method = "fREML", discrete = TRUE) with the family random intercept
# s(fid_re, bs = "re"), dropped for a fold only when < 2 family levels remain.
#
# Depends (defined in sibling helpers_*.R, all loaded by the wrapper before
# any of these run): z_scale, tidy_coefs, lrt_pair, cohort_groups,
# leave_one_cohort_out, build_loo_table, diff_smooth_pooled_gender_simul,
# diff_smooth_significant_intervals, and the constants K_DEFAULT / BS_DEFAULT /
# K_SENS / BS_SENS / N_SIM.
# ---------------------------------------------------------------------------

# ==== shared plumbing ======================================================

# loco_re_term — the per-fold family RE term. Drops the random intercept when
# a fold leaves < 2 distinct family levels (e.g. the fold that holds out
# TwinLife, the only study carrying real family ids). Mirrors the inline
# .loo_re_term used in A_analysis.R / B_analysis.R.
loco_re_term <- function(sub) {
  if (!"fid_re" %in% names(sub)) return("")
  sub$fid_re <- droplevels(factor(sub$fid_re))
  if (nlevels(sub$fid_re) >= 2) "+ s(fid_re, bs = 're')" else ""
}

# rhs_without_re — return a model's full "lhs ~ rhs" string with the family
# RE smooth s(fid_re, ...) removed, so the LOCO engine can re-append the
# fold-specific RE term. Deriving the spec from the fitted primary model (not
# a copied string) guarantees the LOCO refits use the identical specification.
rhs_without_re <- function(mod) {
  f <- paste(deparse(stats::formula(mod)), collapse = " ")
  f <- gsub("\\+\\s*s\\(fid_re[^)]*\\)", "", f)   # drop the family RE smooth
  trimws(gsub("\\s+", " ", f))
}

# .loco_finite — convergence flag for an mgcv bam fit. bam exposes no single
# canonical convergence slot across versions, so a finite focal estimate with
# a positive SE is the operational convergence proxy; the engine status
# ("OK"/"fit_failed") is surfaced separately upstream.
.loco_finite <- function(est, se) is.finite(est) && is.finite(se) && se > 0

# ==== pooled-coefficient LOCO (A4 / B4 PGI x gender) ==================

# loco_pooled_coef — leave-one-study-out refits of a full-family RE model,
# extracting ONE focal coefficient per fold with its 95% CI, supported
# birth-year range, exact p and a convergence flag. `base_rhs` is the model
# formula WITHOUT the family RE term (use rhs_without_re(m_A4_re) etc.).
# Returns the leave_one_cohort_out() table (Full sample row first) with
# delta_vs_full appended by build_loo_table.
loco_pooled_coef <- function(dat, base_rhs, focal_term, label,
                             cohort_var = "cohort") {
  fit_extract <- function(sub) {
    sub$fid_re <- droplevels(factor(sub$fid_re))
    rt  <- loco_re_term(sub)
    m   <- mgcv::bam(stats::as.formula(paste(base_rhs, rt)),
                     data = sub, method = "fREML", discrete = TRUE)
    r   <- tidy_coefs(m, label)
    r   <- r[r$term == focal_term, , drop = FALSE]
    byr <- range(sub$birth_year, na.rm = TRUE)
    if (!nrow(r)) {
      return(c(estimate = NA_real_, SE = NA_real_, ci_lo = NA_real_,
               ci_hi = NA_real_, p = NA_real_, by_min = byr[1], by_max = byr[2],
               converged = 0))
    }
    c(estimate  = round(r$estimate, 5),  SE    = round(r$std.error, 5),
      ci_lo     = round(r$conf.low, 5),  ci_hi = round(r$conf.high, 5),
      p         = round(r$p.value, 6),
      by_min    = byr[1], by_max = byr[2],
      converged = as.numeric(.loco_finite(r$estimate, r$std.error)))
  }
  tab <- leave_one_cohort_out(dat, fit_extract, cohort_var = cohort_var)
  tab <- build_loo_table(tab, estimate_col = "estimate")
  tab$focal_term <- focal_term
  tab$model      <- label
  tab
}

# ==== per-gender GAM LOCO (A5 / B5 nonlinearity) =====================

# .pergender_tests_one — for ONE analytic sample, the omnibus linear-vs-
# nonlinear LRT plus each cell's (female / male) added-smooth curvature LRT
# for the pooled per-gender GAM. lin_rhs / nonlin_rhs are formulas WITHOUT the
# family RE term; cells maps a cell label -> its PGI by-variable (PGI_F/PGI_M).
.pergender_tests_one <- function(sub, lin_rhs, nonlin_rhs, cells,
                                 k = K_DEFAULT, basis = BS_DEFAULT) {
  sub$fid_re <- droplevels(factor(sub$fid_re))
  rt    <- loco_re_term(sub)
  m_lin <- mgcv::bam(stats::as.formula(paste(lin_rhs, rt)),
                     data = sub, method = "fREML", discrete = TRUE)
  m_nl  <- mgcv::bam(stats::as.formula(paste(nonlin_rhs, rt)),
                     data = sub, method = "fREML", discrete = TRUE)
  v_omni <- lrt_pair(m_lin, m_nl)
  mk <- function(test, v) data.frame(
    test = test, LRT_chisq = round(v[["chisq"]], 3),
    LRT_df = round(v[["df"]], 1), LRT_p = signif(v[["p"]], 4),
    nonlinear_sig = isTRUE(!is.na(v[["p"]]) && v[["p"]] < 0.05),
    stringsAsFactors = FALSE)
  rows <- list(mk("omnibus linear vs nonlinear", v_omni))
  for (cn in names(cells)) {
    sm  <- sprintf("s(BYc, by = %s, k = %d, bs = '%s')", cells[[cn]], k, basis)
    m_c <- mgcv::bam(stats::as.formula(paste(lin_rhs, "+", sm, rt)),
                     data = sub, method = "fREML", discrete = TRUE)
    rows[[length(rows) + 1]] <- mk(paste0(cn, " curvature"), lrt_pair(m_lin, m_c))
  }
  dplyr::bind_rows(rows)
}

# loco_pergender_gam — .pergender_tests_one across every leave-one-study-out
# fold (Full sample first). Emits one row per (fold x test) with the supported
# birth-year range, N and status.
loco_pergender_gam <- function(dat, lin_rhs, nonlin_rhs,
                               cells = list(female = "PGI_F", male = "PGI_M"),
                               label, k = K_DEFAULT, basis = BS_DEFAULT,
                               cohort_var = "cohort") {
  folds <- c(list("Full sample" = character(0)), cohort_groups(dat, cohort_var))
  out <- list()
  for (fl in names(folds)) {
    drop_levs <- folds[[fl]]
    full      <- length(drop_levs) == 0
    present   <- full || any(as.character(dat[[cohort_var]]) %in% drop_levs)
    fold_lab  <- if (full) "Full sample" else paste0("-", fl)
    held      <- if (full) "(none)" else fl
    if (!present) {
      out[[fl]] <- data.frame(model = label, fold = fold_lab, held_out = held,
        N = 0L, by_min = NA_real_, by_max = NA_real_, test = NA_character_,
        LRT_chisq = NA_real_, LRT_df = NA_real_, LRT_p = NA_real_,
        nonlinear_sig = NA, status = "N/A - study absent",
        stringsAsFactors = FALSE)
      next
    }
    sub <- if (full) dat else
      dat[!(as.character(dat[[cohort_var]]) %in% drop_levs), , drop = FALSE]
    attr(sub, "by_mean") <- attr(dat, "by_mean")
    byr <- range(sub$birth_year, na.rm = TRUE)
    res <- tryCatch(.pergender_tests_one(sub, lin_rhs, nonlin_rhs, cells, k, basis),
                    error = function(e) NULL)
    if (is.null(res)) {
      out[[fl]] <- data.frame(model = label, fold = fold_lab, held_out = held,
        N = nrow(sub), by_min = byr[1], by_max = byr[2], test = NA_character_,
        LRT_chisq = NA_real_, LRT_df = NA_real_, LRT_p = NA_real_,
        nonlinear_sig = NA, status = "fit_failed", stringsAsFactors = FALSE)
      next
    }
    out[[fl]] <- cbind(
      data.frame(model = label, fold = fold_lab, held_out = held,
                 N = nrow(sub), by_min = byr[1], by_max = byr[2],
                 stringsAsFactors = FALSE),
      res, data.frame(status = "OK", stringsAsFactors = FALSE))
  }
  dplyr::bind_rows(out)
}

# loco_pergender_curves — per-gender PGI-slope curves across birth year for every
# leave-one-study-out fold (Full sample first); the gender analogue of
# A3_LOO_curves / A6_LOO_curves. Per fold, refit the nonlinear per-gender GAM on
# the reduced sample and read off each gender's slope curve via the section
# extractor slope_fn(mod, sub, gender) (.a5_slope / .b5_slope), so the mobility
# model's extra covariates are handled by that section's own extractor.
# nonlin_rhs is the formula WITHOUT the family RE term (the fold RE term is
# appended here).
loco_pergender_curves <- function(dat, nonlin_rhs, slope_fn,
                                  genders = c("female", "male"),
                                  k = K_DEFAULT, basis = BS_DEFAULT,
                                  cohort_var = "cohort") {
  folds <- c(list("Full sample" = character(0)), cohort_groups(dat, cohort_var))
  out <- list()
  for (fl in names(folds)) {
    drop_levs <- folds[[fl]]
    full      <- length(drop_levs) == 0
    if (!full && !any(as.character(dat[[cohort_var]]) %in% drop_levs)) next
    sub <- if (full) dat else
      dat[!(as.character(dat[[cohort_var]]) %in% drop_levs), , drop = FALSE]
    attr(sub, "by_mean") <- attr(dat, "by_mean")
    sub$fid_re <- droplevels(factor(sub$fid_re))
    rt   <- loco_re_term(sub)
    m_nl <- tryCatch(
      mgcv::bam(stats::as.formula(paste(nonlin_rhs, rt)),
                data = sub, method = "fREML", discrete = TRUE),
      error = function(e) NULL)
    if (is.null(m_nl)) next
    cv <- dplyr::bind_rows(lapply(genders, function(g) slope_fn(m_nl, sub, g)))
    cv$fold <- if (full) "Full sample" else paste0("-", fl)
    out[[fl]] <- cv
  }
  res <- dplyr::bind_rows(out)
  res$fold <- factor(res$fold, levels = unique(res$fold))
  res
}

# ==== per-fold female-minus-male difference smooth (B5) ==============

# loco_b5_diffsmooth — per fold, the female-minus-male PGI-slope difference
# smooth over birth year (same orientation and inferential procedure as the
# primary B5 export, diff_smooth_pooled_gender_simul), reported with BOTH a
# pointwise (+/- 1.96 SE) and a simultaneous (N_SIM Monte-Carlo) 95% band,
# trimmed to each fold's supported birth-year range (the smooth is built over
# range(sub$BYc), so trimming is automatic). Also returns, per fold, the
# contiguous birth-year intervals where the SIMULTANEOUS band excludes zero.
# nonlin_rhs = rhs_without_re(m_B5_nonlin_re).
loco_b5_diffsmooth <- function(dat, nonlin_rhs, cohort_var = "cohort",
                               contrast = "Female - Male",
                               dir_pos = "Female > Male",
                               dir_neg = "Female < Male",
                               n_grid = 200, n_sim = N_SIM, seed = 42) {
  folds <- c(list("Full sample" = character(0)), cohort_groups(dat, cohort_var))
  smooths <- list(); intervals <- list()
  for (fl in names(folds)) {
    drop_levs <- folds[[fl]]
    full      <- length(drop_levs) == 0
    present   <- full || any(as.character(dat[[cohort_var]]) %in% drop_levs)
    fold_lab  <- if (full) "Full sample" else paste0("-", fl)
    held      <- if (full) "(none)" else fl
    if (!present) {
      intervals[[fl]] <- data.frame(model = "B5", fold = fold_lab,
        held_out = held, contrast = contrast, ci = "95% simultaneous CI",
        excludes_zero = "n/a", birth_year_interval = "study absent",
        direction = "n/a", max_abs_difference = NA_real_, status = "N/A - study absent",
        stringsAsFactors = FALSE)
      next
    }
    sub <- if (full) dat else
      dat[!(as.character(dat[[cohort_var]]) %in% drop_levs), , drop = FALSE]
    attr(sub, "by_mean") <- attr(dat, "by_mean")
    ds <- tryCatch({
      sub2 <- sub; sub2$fid_re <- droplevels(factor(sub2$fid_re))
      rt   <- loco_re_term(sub2)
      m    <- mgcv::bam(stats::as.formula(paste(nonlin_rhs, rt)),
                        data = sub2, method = "fREML", discrete = TRUE)
      diff_smooth_pooled_gender_simul(m, sub2, n_grid = n_grid,
                                      n_sim = n_sim, seed = seed)
    }, error = function(e) NULL)
    if (is.null(ds)) {
      intervals[[fl]] <- data.frame(model = "B5", fold = fold_lab,
        held_out = held, contrast = contrast, ci = "95% simultaneous CI",
        excludes_zero = "n/a", birth_year_interval = "fit failed",
        direction = "n/a", max_abs_difference = NA_real_, status = "fit_failed",
        stringsAsFactors = FALSE)
      next
    }
    ds$ci_lo_pointwise <- ds$diff - 1.96 * ds$se
    ds$ci_hi_pointwise <- ds$diff + 1.96 * ds$se
    ds$model <- "B5"; ds$fold <- fold_lab; ds$held_out <- held
    smooths[[fl]] <- ds
    iv <- diff_smooth_significant_intervals(
      ds, contrast = contrast, ci_label = "95% simultaneous CI",
      dir_pos = dir_pos, dir_neg = dir_neg)
    iv <- cbind(data.frame(model = "B5", fold = fold_lab, held_out = held,
                           stringsAsFactors = FALSE),
                iv, data.frame(status = "OK", stringsAsFactors = FALSE))
    intervals[[fl]] <- iv
  }
  sm <- dplyr::bind_rows(smooths)
  if (nrow(sm)) sm$fold <- factor(sm$fold, levels = unique(sm$fold))
  list(smooth = sm, intervals = dplyr::bind_rows(intervals))
}

# ==== common complete-case alternative-PGI sample ====================

# altpgi_common_sample — restrict `dat` to rows complete on all three raw PGIs
# and the ancestry PCs, then re-residualise each raw PGI on crossPC1-10 and
# re-z-score WITHIN this common sample (so PGI-Education, PGI-Cognition and
# PGI-Noncognitive are all standardised on the identical N). Overwrites the
# *_z columns. Errors clearly if the raw PGIs / PCs are absent (the caller
# should then rebuild the sample via prepare_samples.R).
altpgi_common_sample <- function(dat,
                                 pgis_raw = c("PGI_Edu", "PGI_Cog", "PGI_nonCog"),
                                 also_complete = character(0)) {
  xpcs <- paste0("crossPC", 1:10)
  need <- c(pgis_raw, xpcs)
  miss <- setdiff(need, names(dat))
  if (length(miss))
    stop("altpgi_common_sample: missing column(s) ",
         paste(miss, collapse = ", "),
         " - rebuild the sample via prepare_samples.R first.")
  keep_cols <- c(pgis_raw, xpcs, intersect(also_complete, names(dat)))
  cc  <- stats::complete.cases(dat[, keep_cols, drop = FALSE])
  sub <- dat[cc, , drop = FALSE]
  resid_on_pcs <- function(d, pgi) {
    stats::residuals(stats::lm(
      stats::as.formula(paste(pgi, "~", paste(xpcs, collapse = " + "))),
      data = d))
  }
  for (p in pgis_raw) sub[[paste0(p, "_z")]] <- z_scale(resid_on_pcs(sub, p))
  sub
}

# ==== finding-matched RE alternative-PGI families ====================

# .sub_pgi — substitute the placeholder token "PGI" in a formula/term string
# with a concrete PGI z-variable name. The templates use "PGI" only as the
# placeholder (no other literal "PGI" substring), so a single gsub is safe.
.sub_pgi <- function(template, pgi) gsub("PGI", pgi, template, fixed = TRUE)

# altpgi_re_family — refit one RE focal model per PGI on the common sample
# (substituting the PGI variable), extract one or more focal coefficients, and
# Benjamini-Hochberg-adjust p WITHIN each focal-term family across the three
# PGIs. `rhs_template` and each entry of `focal_map` (family label -> focal
# coefficient template) use the "PGI" placeholder.
altpgi_re_family <- function(dat, rhs_template, focal_map, label,
                             pgis = c("PGI_Edu_z", "PGI_Cog_z", "PGI_nonCog_z"),
                             re = "+ s(fid_re, bs = 're')") {
  rows <- list()
  for (pgi in pgis) {
    rhs <- .sub_pgi(rhs_template, pgi)
    m   <- mgcv::bam(stats::as.formula(paste(rhs, re)),
                     data = dat, method = "fREML", discrete = TRUE)
    co  <- tidy_coefs(m, label)
    n_m <- tryCatch(stats::nobs(m), error = function(e) NA_integer_)
    for (fam in names(focal_map)) {
      ft <- .sub_pgi(focal_map[[fam]], pgi)
      r  <- co[co$term == ft, , drop = FALSE]
      rows[[paste(fam, pgi)]] <- data.frame(
        analysis = label, family = fam, focal_term = ft,
        PGI = sub("_z$", "", pgi), PGI_z = pgi,
        estimate  = if (nrow(r)) round(r$estimate, 5)  else NA_real_,
        se        = if (nrow(r)) round(r$std.error, 5) else NA_real_,
        statistic = if (nrow(r)) round(r$statistic, 4) else NA_real_,
        ci_lo     = if (nrow(r)) round(r$conf.low, 5)  else NA_real_,
        ci_hi     = if (nrow(r)) round(r$conf.high, 5) else NA_real_,
        p_raw     = if (nrow(r)) r$p.value             else NA_real_,
        N = n_m, stringsAsFactors = FALSE)
    }
  }
  out <- dplyr::bind_rows(rows)
  out <- do.call(rbind, lapply(split(out, out$family), function(d) {
    d$p_BH <- stats::p.adjust(d$p_raw, method = "BH"); d
  }))
  rownames(out) <- NULL
  out[order(match(out$family, names(focal_map)),
            match(out$PGI_z, pgis)), , drop = FALSE]
}

# .mk_nd0 — a one-row newdata with every predictor column of `mod` set to 0
# (fid_re to its first level), then the named `overrides` applied. Used to
# build lpmatrix slope contrasts at specified region/gender cell codings.
.mk_nd0 <- function(mod, overrides = list()) {
  mv <- names(mod$model)[-1]              # drop the response
  nd <- as.data.frame(as.list(stats::setNames(rep(0, length(mv)), mv)))
  if ("fid_re" %in% mv) {
    fl <- levels(mod$model$fid_re); nd$fid_re <- factor(fl[1], levels = fl)
  }
  for (k in names(overrides)) nd[[k]] <- overrides[[k]]
  nd
}

# re_pgi_slopes — per-cell standardised PGI slope (per 1 SD) from a fitted RE
# model, via an lpmatrix contrast (PGI = 1 vs 0) at each cell's covariate
# coding. `cells` maps a cell label -> a list of covariate overrides, e.g.
# list(West = list(east_west_c = -0.5), East = list(east_west_c = 0.5)).
re_pgi_slopes <- function(mod, pgi_z, cells) {
  do.call(rbind, lapply(names(cells), function(cn) {
    nd0 <- .mk_nd0(mod, cells[[cn]]); nd0[[pgi_z]] <- 0
    nd1 <- nd0; nd1[[pgi_z]] <- 1
    D   <- stats::predict(mod, nd1, type = "lpmatrix") -
           stats::predict(mod, nd0, type = "lpmatrix")
    b <- stats::coef(mod); V <- stats::vcov(mod)
    est <- as.numeric(D %*% b); se <- sqrt(as.numeric(D %*% V %*% t(D)))
    data.frame(cell = cn, PGI = sub("_z$", "", pgi_z), PGI_z = pgi_z,
      slope = round(est, 5), se = round(se, 5),
      ci_lo = round(est - 1.96 * se, 5), ci_hi = round(est + 1.96 * se, 5),
      p = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)
  }))
}

# altpgi_slopes_family — the per-cell PGI slopes for every PGI (one RE fit per
# PGI on the common sample), stacked into a tidy table. Used for the East/West
# mobility slopes (Finding 2) and the female/male gender slopes (Finding 3).
altpgi_slopes_family <- function(dat, rhs_template, cells, label,
                                 pgis = c("PGI_Edu_z", "PGI_Cog_z", "PGI_nonCog_z"),
                                 re = "+ s(fid_re, bs = 're')") {
  do.call(rbind, lapply(pgis, function(pgi) {
    rhs <- .sub_pgi(rhs_template, pgi)
    m   <- mgcv::bam(stats::as.formula(paste(rhs, re)),
                     data = dat, method = "fREML", discrete = TRUE)
    cbind(data.frame(analysis = label, stringsAsFactors = FALSE),
          re_pgi_slopes(m, pgi, cells))
  }))
}

# ==== alternative-PGI per-gender GAM (A5 / B5 reproduction) ==========

# altpgi_pergender_gam — for each PGI on the common sample, the pooled per-
# gender GAM omnibus + female/male curvature LRTs (BH-adjusted across the three
# PGIs on the omnibus), and — when want_diffsmooth = TRUE (B5) — the female-
# minus-male difference smooth with a simultaneous 95% band plus its
# excludes-zero birth-year intervals. lin_rhs / nonlin_rhs reference the
# PGI_F / PGI_M cells, which this function rebuilds from each PGI in turn.
altpgi_pergender_gam <- function(dat, lin_rhs, nonlin_rhs, label,
                                 pgis = c("PGI_Edu_z", "PGI_Cog_z", "PGI_nonCog_z"),
                                 cells = list(female = "PGI_F", male = "PGI_M"),
                                 want_diffsmooth = FALSE,
                                 k = K_DEFAULT, basis = BS_DEFAULT,
                                 n_sim = N_SIM, seed = 42) {
  tests <- list(); smooths <- list(); intervals <- list()
  for (pgi in pgis) {
    sub <- dat
    sub$PGI_F <- sub[[pgi]] * (1L - sub$gender_01)
    sub$PGI_M <- sub[[pgi]] *        sub$gender_01
    tr <- tryCatch(.pergender_tests_one(sub, lin_rhs, nonlin_rhs, cells, k, basis),
                   error = function(e) NULL)
    if (!is.null(tr)) {
      tr <- cbind(data.frame(analysis = label, PGI = sub("_z$", "", pgi),
                             PGI_z = pgi, stringsAsFactors = FALSE), tr)
      tests[[pgi]] <- tr
    }
    if (want_diffsmooth) {
      ds <- tryCatch({
        s2 <- sub; s2$fid_re <- droplevels(factor(s2$fid_re))
        rt <- loco_re_term(s2)
        m  <- mgcv::bam(stats::as.formula(paste(nonlin_rhs, rt)),
                        data = s2, method = "fREML", discrete = TRUE)
        diff_smooth_pooled_gender_simul(m, s2, n_sim = n_sim, seed = seed)
      }, error = function(e) NULL)
      if (!is.null(ds)) {
        ds$ci_lo_pointwise <- ds$diff - 1.96 * ds$se
        ds$ci_hi_pointwise <- ds$diff + 1.96 * ds$se
        ds$analysis <- label; ds$PGI <- sub("_z$", "", pgi); ds$PGI_z <- pgi
        smooths[[pgi]] <- ds
        iv <- diff_smooth_significant_intervals(ds, contrast = "Female - Male",
                ci_label = "95% simultaneous CI",
                dir_pos = "Female > Male", dir_neg = "Female < Male")
        intervals[[pgi]] <- cbind(
          data.frame(analysis = label, PGI = sub("_z$", "", pgi), PGI_z = pgi,
                     stringsAsFactors = FALSE), iv)
      }
    }
  }
  tests_df <- dplyr::bind_rows(tests)
  # BH across the three PGIs on the omnibus nonlinearity test only.
  if (nrow(tests_df)) {
    om <- tests_df$test == "omnibus linear vs nonlinear"
    tests_df$p_BH <- NA_real_
    tests_df$p_BH[om] <- stats::p.adjust(tests_df$LRT_p[om], method = "BH")
  }
  list(tests = tests_df,
       diffsmooth = dplyr::bind_rows(smooths),
       diff_intervals = dplyr::bind_rows(intervals))
}
