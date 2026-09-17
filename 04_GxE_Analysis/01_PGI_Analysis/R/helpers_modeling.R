# helpers_modeling.R
# Model extraction, nested LRT, slope curves, GAM difference smooths,
# cluster-robust SEs.

# ---- Model-extraction helpers -------------------------------------------

tidy_coefs <- function(mod, label) {
  # broom::tidy.gam() defaults to parametric = FALSE, which returns smooth-
  # term tests (edf, F, p) instead of the fixed-effect coefficients. For
  # GAMs we want the parametric coefficients here (Estimate, SE, t, p,
  # 95% CI) so headline tables and focal-row lookups work. lm/glm dispatch
  # already returns coefficients by default; no flag needed.
  td <- if (inherits(mod, "gam")) {
    broom::tidy(mod, parametric = TRUE, conf.int = TRUE)
  } else {
    broom::tidy(mod, conf.int = TRUE)
  }
  td |> mutate(model = label, .before = 1)
}

tidy_smooths <- function(mod, label) {
  stab <- summary(mod)$s.table
  if (is.null(stab) || nrow(stab) == 0) return(NULL)
  as.data.frame(stab) |>
    tibble::rownames_to_column("smooth_term") |>
    mutate(model = label, .before = 1)
}

model_fit_row <- function(mod, label, n = NULL) {
  if (is.null(n)) n <- tryCatch(nobs(mod), error = function(e) NA_integer_)
  data.frame(
    model         = label, n = n,
    AIC           = tryCatch(AIC(mod),             error = function(e) NA_real_),
    BIC           = tryCatch(BIC(mod),             error = function(e) NA_real_),
    logLik        = tryCatch(as.numeric(logLik(mod)), error = function(e) NA_real_),
    R2_or_devExpl = if (inherits(mod, "gam")) summary(mod)$dev.expl
                    else if (inherits(mod, "lm")) summary(mod)$r.squared
                    else NA_real_,
    stringsAsFactors = FALSE
  )
}

compare_nested_models <- function(reduced_mod, full_mod, reduced_label, full_label) {
  lrt <- anova(reduced_mod, full_mod, test = "Chisq")

  get_val <- function(cols) {
    for (cn in cols) {
      if (cn %in% colnames(lrt)) return(lrt[2, cn])
    }
    NA_real_
  }

  data.frame(
    reduced_model = reduced_label,
    full_model    = full_label,
    n             = tryCatch(nobs(full_mod), error = function(e) NA_integer_),
    AIC_reduced   = tryCatch(AIC(reduced_mod), error = function(e) NA_real_),
    AIC_full      = tryCatch(AIC(full_mod),    error = function(e) NA_real_),
    BIC_reduced   = tryCatch(BIC(reduced_mod), error = function(e) NA_real_),
    BIC_full      = tryCatch(BIC(full_mod),    error = function(e) NA_real_),
    logLik_reduced = tryCatch(as.numeric(logLik(reduced_mod)), error = function(e) NA_real_),
    logLik_full    = tryCatch(as.numeric(logLik(full_mod)),    error = function(e) NA_real_),
    LRT_df        = get_val(c("Df","df","Chi Df")),
    LRT_chisq     = get_val(c("Deviance","X2","Chi.sq","Chisq","LRT")),
    LRT_p         = get_val(c("Pr(>Chi)","p-value","P(>|Chi|)","Pr(>Chisq)")),
    stringsAsFactors = FALSE
  )
}

extract_region_pgi_effects <- function(mod, dat, pgi_var = "PGI_Edu",
                                       regions = c("West", "East"),
                                       gender = "female",
                                       byc_value = 0,
                                       par_edu_value = 0) {
  make_nd <- function(region, pgi_value) {
    is_east <- region == "East"
    is_male <- gender == "male"
    nd <- data.frame(
      BYc = byc_value,
      PGI_Edu = 0,
      PGI_Cog = 0,
      PGI_nonCog = 0,
      PGI_Height = 0,
      stringsAsFactors = FALSE
    )
    nd[[pgi_var]] <- pgi_value

    if ("east_west" %in% names(mod$model))
      nd$east_west <- factor(region, levels = levels(dat$east_west))
    if ("gender" %in% names(mod$model))
      nd$gender <- factor(gender, levels = levels(dat$gender))
    if ("east_west_c" %in% names(mod$model))
      nd$east_west_c <- if (is_east) 0.5 else -0.5
    if ("gender_c" %in% names(mod$model))
      nd$gender_c <- if (is_male) -0.5 else 0.5
    if ("east" %in% names(mod$model)) nd$east <- as.integer(is_east)
    if ("gender_01" %in% names(mod$model)) nd$gender_01 <- as.integer(is_male)
    if ("parental_edu_z_kernel" %in% names(mod$model))
      nd$parental_edu_z_kernel <- par_edu_value
    if ("ParEdu_z" %in% names(mod$model)) nd$ParEdu_z <- par_edu_value
    if ("fid_re" %in% names(mod$model)) {
      fid_levels <- levels(mod$model$fid_re)
      nd$fid_re <- factor(fid_levels[1], levels = fid_levels)
    }

    for (v in c("PGI_W","PGI_E","PGI_WF","PGI_WM","PGI_EF","PGI_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0
    # Region/gender-specific ParEdu control cells (e.g. the B6 by-smooths) –
    # held at 0; identical in the PGI=0 and PGI=1 newdata, so they cancel.
    for (v in c("ParEdu_W","ParEdu_E","ParEdu_WF","ParEdu_WM","ParEdu_EF","ParEdu_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0

    if ("PGI_W"  %in% names(nd)) nd$PGI_W  <- pgi_value * as.numeric(!is_east)
    if ("PGI_E"  %in% names(nd)) nd$PGI_E  <- pgi_value * as.numeric( is_east)
    if ("PGI_WF" %in% names(nd)) nd$PGI_WF <- pgi_value * as.numeric(!is_east & gender == "female")
    if ("PGI_WM" %in% names(nd)) nd$PGI_WM <- pgi_value * as.numeric(!is_east & gender == "male")
    if ("PGI_EF" %in% names(nd)) nd$PGI_EF <- pgi_value * as.numeric( is_east & gender == "female")
    if ("PGI_EM" %in% names(nd)) nd$PGI_EM <- pgi_value * as.numeric( is_east & gender == "male")

    nd
  }

  V <- vcov(mod)
  out <- lapply(regions, function(region) {
    nd0 <- make_nd(region, 0)
    nd1 <- make_nd(region, 1)
    X0 <- predict(mod, newdata = nd0, type = "lpmatrix")
    X1 <- predict(mod, newdata = nd1, type = "lpmatrix")
    D  <- X1 - X0
    fit <- as.vector(D %*% coef(mod))
    se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))
    data.frame(
      east_west = region,
      slope = fit,
      se = se,
      ci_lo = fit - 1.96 * se,
      ci_hi = fit + 1.96 * se,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

extract_region_gender_pgi_effects <- function(mod, dat, pgi_var = "PGI_Edu",
                                              regions = c("West", "East"),
                                              genders = c("female", "male"),
                                              byc_value = 0,
                                              par_edu_value = 0) {
  out <- lapply(regions, function(region) {
    bind_rows(lapply(genders, function(gen) {
      extract_region_pgi_effects(
        mod = mod,
        dat = dat,
        pgi_var = pgi_var,
        regions = region,
        gender = gen,
        byc_value = byc_value,
        par_edu_value = par_edu_value
      ) |>
        mutate(gender = gen,
               group = paste(region, gen))
    }))
  })
  bind_rows(out)
}

extract_lm_region_pgi_effects <- function(mod, dat,
                                          regions = c("West", "East"),
                                          gender = "female",
                                          byc_z_value = 0,
                                          par_edu_value = 0) {
  mm_terms <- delete.response(terms(mod))
  V <- vcov(mod)
  is_male <- gender == "male"

  bind_rows(lapply(regions, function(region) {
    is_east <- region == "East"
    nd0 <- data.frame(
      PGI_Edu_z   = 0,
      BYc_z       = byc_z_value,
      east_west_c = if (is_east) 0.5 else -0.5,
      gender_c    = if (is_male) -0.5 else 0.5,
      parental_edu_z_kernel = par_edu_value,
      stringsAsFactors = FALSE
    )
    # Family random-effects placeholder (cancels in nd1 - nd0)
    if ("fid_re" %in% names(mod$model)) {
      fid_levels <- levels(mod$model$fid_re)
      nd0$fid_re <- factor(fid_levels[1], levels = fid_levels)
    }
    nd1 <- nd0
    nd1$PGI_Edu_z <- 1

    if (inherits(mod, "gam")) {
      X0 <- predict(mod, newdata = nd0, type = "lpmatrix")
      X1 <- predict(mod, newdata = nd1, type = "lpmatrix")
    } else {
      X0 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd0, xlev = mod$xlevels))
      X1 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd1, xlev = mod$xlevels))
    }
    D  <- X1 - X0

    fit <- as.vector(D %*% coef(mod))
    se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

    data.frame(
      east_west = region,
      slope = fit,
      se = se,
      ci_lo = fit - 1.96 * se,
      ci_hi = fit + 1.96 * se,
      stringsAsFactors = FALSE
    )
  }))
}

# ---- Nested LRT table (null -> linear -> nonlinear) -----------------
# Works for both lm and gam objects.
#
# The LRT is computed manually rather than via anova() because anova.gam
# derives the comparison df from the difference in *effective* degrees of
# freedom. For models containing a random effect (s(fid_re, bs="re")) the
# RE soaks up different amounts of edf between the two models, which can
# produce a negative or NA df (and therefore no p-value) for comparisons
# that differ in fixed effects — exactly the null -> linear case here.
#
# Manual LRT:
#   chisq = 2 * (logLik_full - logLik_reduced)
#   df    = length(coef(full)) - length(coef(reduced))
# The coefficient-count df is exact for unpenalized fixed-effect additions
# (null -> linear: 4) and the basis dimension — a conservative upper bound
# on the effective df — for penalized smooth additions (linear ->
# nonlinear). For a valid fixed-effects LRT the inputs must be ML fits;
# REML log-likelihoods are not comparable across models with different
# fixed effects.

# lrt_pair — manual likelihood-ratio test between a reduced and a full
# (nested) model. Returns c(chisq, df, p). Robust to the effective-df
# accounting of anova.gam for RE/penalized models: df is the difference
# in coefficient counts (exact for unpenalized additions; basis dimension
# / conservative upper bound for penalized smooths). For a valid
# fixed-effects LRT the inputs must be ML fits.
lrt_pair <- function(reduced, full) {
  ll_r <- tryCatch(as.numeric(logLik(reduced)), error = function(e) NA_real_)
  ll_f <- tryCatch(as.numeric(logLik(full)),    error = function(e) NA_real_)
  nc_r <- tryCatch(length(stats::coef(reduced)), error = function(e) NA_integer_)
  nc_f <- tryCatch(length(stats::coef(full)),    error = function(e) NA_integer_)
  chisq <- if (is.na(ll_r) || is.na(ll_f)) NA_real_ else max(0, 2 * (ll_f - ll_r))
  df    <- if (is.na(nc_r) || is.na(nc_f)) NA_real_ else as.numeric(nc_f - nc_r)
  p     <- if (!is.na(chisq) && !is.na(df) && df > 0)
             stats::pchisq(chisq, df, lower.tail = FALSE) else NA_real_
  c(chisq = chisq, df = df, p = p)
}

nested_lrt_table <- function(null_mod, lin_mod, nonlin_mod, label = "A3") {
  v1 <- lrt_pair(null_mod, lin_mod)
  v2 <- lrt_pair(lin_mod,  nonlin_mod)

  get_fit <- function(mod) {
    dev_expl <- if (inherits(mod, "gam")) summary(mod)$dev.expl
                else if (inherits(mod, "lm")) summary(mod)$r.squared
                else NA_real_
    c(AIC     = tryCatch(AIC(mod),                error=function(e) NA_real_),
      BIC     = tryCatch(BIC(mod),                error=function(e) NA_real_),
      logLik  = tryCatch(as.numeric(logLik(mod)), error=function(e) NA_real_),
      devExpl = dev_expl)
  }

  f0 <- get_fit(null_mod); f1 <- get_fit(lin_mod); f2 <- get_fit(nonlin_mod)

  data.frame(
    model      = paste0(label, c("_null","_lin","_nonlin")),
    AIC        = c(f0["AIC"],    f1["AIC"],    f2["AIC"]),
    BIC        = c(f0["BIC"],    f1["BIC"],    f2["BIC"]),
    logLik     = c(f0["logLik"], f1["logLik"], f2["logLik"]),
    dev_expl   = c(f0["devExpl"],f1["devExpl"],f2["devExpl"]),
    LRT_chisq  = c(NA_real_,     v1["chisq"],  v2["chisq"]),
    LRT_df     = c(NA_real_,     v1["df"],     v2["df"]),
    LRT_p      = c(NA_real_,     v1["p"],      v2["p"]),
    comparison = c(NA_character_,"PGI baseline vs linear (cohort slope)","linear vs nonlinear"),
    stringsAsFactors = FALSE
  )
}

# ---- Baseline-convention guard for the nested ladders ------------
# The nested LRT ladder tests the SHAPE of the PGI slope over birth cohort.
# For the baseline -> linear rung to isolate the cohort SLOPE (rather than
# conflate it with the PGI main effect), the PGI main-effect cell terms must
# be present in the baseline model already. This guard fails the run if any
# ladder's baseline omits a PGI main-effect cell, so the convention is
# enforced identically across every ladder (A3/A5/A6 and B3/B5/B6) and cannot
# silently regress. It inspects the fitted models' coefficient names, so it
# stays in lockstep with the actual formulas rather than a copied string.
assert_pgi_in_baseline <- function(baseline_mod, linear_mod, pgi_cells, tag) {
  bn <- names(stats::coef(baseline_mod))
  ln <- names(stats::coef(linear_mod))
  for (cell in pgi_cells) {
    if (!(cell %in% bn))
      stop(sprintf(paste0("[%s] PGI main-effect cell '%s' is missing from the ",
        "baseline model: the baseline -> linear LRT would conflate the PGI ",
        "main effect with the cohort slope. Add '%s' to the baseline formula."),
        tag, cell, cell))
    if (!(cell %in% ln))
      stop(sprintf("[%s] PGI main-effect cell '%s' is missing from the linear model.",
                   tag, cell))
  }
  invisible(TRUE)
}

# ---- Non-parametric kernel-regression PGI slopes ---------------------
# Estimates beta_PGI(birth_year) using kernel-weighted local linear regression
# within each group cell.

compute_binned_slopes <- function(dat, outcome, pgi = "PGI_Edu",
                                   covariates = character(0),
                                   group_vars = "east_west",
                                   bandwidth  = NP_BANDWIDTH,
                                   n_grid     = 200,
                                   min_n      = 20,
                                   min_eff_n  = 10) {
  needed <- c(outcome, pgi, "birth_year", covariates)
  dat2 <- dat |>
    filter(if_all(all_of(needed), ~ !is.na(.x)))

  if (length(group_vars) > 0) {
    dat2 <- dat2 |>
      filter(if_all(all_of(group_vars), ~ !is.na(.x)))
  }

  rhs <- if (length(covariates) > 0)
    paste(c(pgi, covariates), collapse = " + ") else pgi
  fmla_str <- paste(outcome, "~", rhs)

  kernel_w <- function(x, x0, h) {
    dnorm((x - x0) / h)
  }

  fit_local_slope <- function(d, x0) {
    w <- kernel_w(d$birth_year, x0, bandwidth)
    n_eff <- sum(w > max(w) * exp(-2), na.rm = TRUE)
    if (sum(w > 0, na.rm = TRUE) < min_n || n_eff < min_eff_n) return(NULL)

    # Build the formula INSIDE this scope so its environment captures `w`;
    # otherwise lm(weights = w) errors with "object 'w' not found" because
    # the formula's enclosing environment doesn't contain w.
    fmla <- as.formula(fmla_str)
    m <- tryCatch(
      lm(fmla, data = d, weights = w),
      error = function(e) NULL
    )
    if (is.null(m)) return(NULL)

    b <- coef(summary(m))
    if (!pgi %in% rownames(b)) return(NULL)

    data.frame(
      birth_year = x0,
      slope = b[pgi, "Estimate"],
      se = b[pgi, "Std. Error"],
      n = nrow(d),
      n_eff = n_eff,
      stringsAsFactors = FALSE
    )
  }

  dat2 |>
    group_by(across(all_of(group_vars))) |>
    group_modify(~ {
      d <- .x
      if (nrow(d) < min_n) return(data.frame())
      x_grid <- seq(min(d$birth_year), max(d$birth_year), length.out = n_grid)
      out <- lapply(x_grid, function(x0) fit_local_slope(d, x0))
      bind_rows(out)
    }) |>
    ungroup()
}

# ---- Linear-model PGI curve (for A1 / B1 plots) ---------------------
# Returns the PGI slope as a function of BYc_z, by region, holding gender fixed.
# Difference predict(PGI=1) - predict(PGI=0) at a grid of BYc_z values.

extract_lm_A1_curves <- function(mod, dat, n_grid = 200, gender_value = "female") {
  BYc_z_range <- range(dat$BYc_z, na.rm = TRUE)
  BYc_z_grid  <- seq(BYc_z_range[1], BYc_z_range[2], length.out = n_grid)
  by_sd   <- sd(dat$birth_year,  na.rm = TRUE)
  by_mean <- mean(dat$birth_year, na.rm = TRUE)

  V <- vcov(mod)
  mm_terms <- delete.response(terms(mod))

  is_male <- gender_value == "male"
  results <- list()
  for (reg in c("West", "East")) {
    is_east <- reg == "East"
    nd0 <- data.frame(
      PGI_Edu_z    = 0,
      BYc_z        = BYc_z_grid,
      east_west_c  = if (is_east) 0.5 else -0.5,
      gender_c     = if (is_male) -0.5 else 0.5,
      PGI_Cog_z    = 0,
      PGI_nonCog_z = 0,
      PGI_Height_z = 0,
      stringsAsFactors = FALSE
    )
    if ("parental_edu_z_kernel" %in% names(mod$model)) nd0$parental_edu_z_kernel <- 0
    if ("ParEdu_z" %in% names(mod$model)) nd0$ParEdu_z <- 0
    nd1 <- nd0; nd1$PGI_Edu_z <- 1

    slope <- tryCatch(
      as.numeric(predict(mod, newdata = nd1) - predict(mod, newdata = nd0)),
      error = function(e) rep(NA_real_, n_grid)
    )

    se <- tryCatch({
      X0 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd0, xlev = mod$xlevels))
      X1 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd1, xlev = mod$xlevels))
      D  <- X1 - X0
      sqrt(pmax(1e-12, rowSums((D %*% V) * D)))
    }, error = function(e) rep(NA_real_, n_grid))

    results[[reg]] <- data.frame(
      birth_year = BYc_z_grid * by_sd + by_mean,
      BYc_z      = BYc_z_grid,
      east_west  = reg,
      gender     = gender_value,
      slope      = slope,
      se         = se,
      ci_lo      = slope - 1.96 * se,
      ci_hi      = slope + 1.96 * se
    )
  }
  bind_rows(results)
}

# ---- A2 PGI effects per east_west x reunif cell -------------------------
# Returns the PGI slope for each east_west x period cell (female reference).

extract_a2_pgi_effects <- function(mod, dat) {
  V <- vcov(mod)
  mm_terms <- delete.response(terms(mod))

  results <- list()
  for (reg in c("West", "East")) {
    is_east <- reg == "East"
    for (per in c("pre", "post")) {
      nd0 <- data.frame(
        PGI_Edu_z    = 0,
        reunif       = factor(per, levels = levels(dat$reunif)),
        east_west_c  = if (is_east) 0.5 else -0.5,
        gender_c     = 0.5,                       # female reference
        PGI_Cog_z    = 0,
        PGI_nonCog_z = 0,
        PGI_Height_z = 0,
        stringsAsFactors = FALSE
      )
      if ("parental_edu_z_kernel" %in% names(mod$model)) nd0$parental_edu_z_kernel <- 0
      if ("ParEdu_z" %in% names(mod$model)) nd0$ParEdu_z <- 0
      if ("fid_re" %in% names(mod$model)) {
        fid_levels <- levels(mod$model$fid_re)
        nd0$fid_re <- factor(fid_levels[1], levels = fid_levels)
      }
      nd1 <- nd0; nd1$PGI_Edu_z <- 1

      slope <- tryCatch(
        as.numeric(predict(mod, newdata = nd1) - predict(mod, newdata = nd0)),
        error = function(e) NA_real_
      )

      se <- tryCatch({
        if (inherits(mod, "gam")) {
          X0 <- predict(mod, newdata = nd0, type = "lpmatrix")
          X1 <- predict(mod, newdata = nd1, type = "lpmatrix")
        } else {
          X0 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd0, xlev = mod$xlevels))
          X1 <- model.matrix(mm_terms, data = model.frame(mm_terms, nd1, xlev = mod$xlevels))
        }
        D  <- X1 - X0
        sqrt(max(1e-12, sum((D %*% V) * D)))
      }, error = function(e) NA_real_)

      results[[paste(reg, per)]] <- data.frame(
        east_west = reg,
        period    = factor(per, levels = c("pre", "post")),
        slope     = slope,
        se        = se,
        ci_lo     = slope - 1.96 * se,
        ci_hi     = slope + 1.96 * se,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(results)
}

# ---- GAM slope extraction -----------------------------------------------
# Predicts at PGI=1 vs PGI=0 across a BYc grid. Captures linear + smooth parts.
# east_west: "West" or "East" | gender: "female" or "male"

extract_pgi_slope <- function(mod, dat, east_west, gender = "female",
                              pgi_var = "PGI_Edu", n_grid = 200) {
  BYc_range <- range(dat$BYc, na.rm = TRUE)
  BYc_grid  <- seq(BYc_range[1], BYc_range[2], length.out = n_grid)
  by_mean   <- attr(dat, "by_mean")
  if (is.null(by_mean)) by_mean <- mean(dat$birth_year, na.rm = TRUE)

  is_east <- east_west == "East"
  is_male <- gender == "male"

  # Core newdata - include both raw and z-scored PGI columns so the
  # function works regardless of which one the model formula uses.
  by_sd_local <- sd(dat$birth_year, na.rm = TRUE)
  nd <- data.frame(
    BYc          = BYc_grid,
    BYc_z        = if (is.finite(by_sd_local) && by_sd_local > 0)
                     BYc_grid / by_sd_local else BYc_grid,
    PGI_Edu      = 0,
    PGI_Edu_z    = 0,
    PGI_Cog      = 0,
    PGI_Cog_z    = 0,
    PGI_nonCog   = 0,
    PGI_nonCog_z = 0,
    PGI_Height   = 0,
    PGI_Height_z = 0,
    stringsAsFactors = FALSE
  )

  # Factor predictors (if present in model)
  if ("east_west" %in% names(mod$model))
    nd$east_west <- factor(east_west, levels = levels(dat$east_west))
  if ("gender" %in% names(mod$model))
    nd$gender <- factor(gender, levels = levels(dat$gender))
  # Family random-effects placeholder (any valid level - the random-effects
  # contribution to predict(nd1) - predict(nd0) cancels because nd1 and
  # nd0 share the same fid_re value)
  if ("fid_re" %in% names(mod$model)) {
    fid_levels <- levels(mod$model$fid_re)
    nd$fid_re <- factor(fid_levels[1], levels = fid_levels)
  }

  # Contrast-coded numeric predictors (ANALYSIS_PLAN.md "Note on coding")
  if ("east_west_c" %in% names(mod$model))
    nd$east_west_c <- if (is_east) 0.5 else -0.5
  if ("gender_c" %in% names(mod$model))
    nd$gender_c <- if (is_male) -0.5 else 0.5

  # 0/1 cell-indicator helpers (only if present in model)
  if ("east"      %in% names(mod$model)) nd$east      <- as.integer(is_east)
  if ("gender_01" %in% names(mod$model)) nd$gender_01 <- as.integer(is_male)

  # Optional covariates
  if ("parental_edu_z_kernel" %in% names(mod$model)) nd$parental_edu_z_kernel <- 0
  if ("ParEdu_z" %in% names(mod$model)) nd$ParEdu_z <- 0

  # east_west x gender PGI interaction terms – start at 0
  for (v in c("PGI_W","PGI_E","PGI_WF","PGI_WM","PGI_EF","PGI_EM"))
    if (v %in% names(mod$model)) nd[[v]] <- 0

  # Region/gender-specific ParEdu control cells (e.g. B3 by-smooths
  # s(BYc, by = ParEdu_W/ParEdu_E)) – held at 0; they are identical in
  # nd0 and nd1, so they cancel in the X1 - X0 contrast below.
  for (v in c("ParEdu_W","ParEdu_E","ParEdu_WF","ParEdu_WM","ParEdu_EF","ParEdu_EM"))
    if (v %in% names(mod$model)) nd[[v]] <- 0

  nd0 <- nd; nd1 <- nd
  nd0[[pgi_var]] <- 0; nd1[[pgi_var]] <- 1

  # Set the correct PGI_* interaction term for this east_west x gender combo
  if ("PGI_W"  %in% names(nd1)) nd1$PGI_W  <- as.numeric(!is_east)
  if ("PGI_E"  %in% names(nd1)) nd1$PGI_E  <- as.numeric( is_east)
  if ("PGI_WF" %in% names(nd1)) nd1$PGI_WF <- as.numeric(!is_east & !is_male)
  if ("PGI_WM" %in% names(nd1)) nd1$PGI_WM <- as.numeric(!is_east &  is_male)
  if ("PGI_EF" %in% names(nd1)) nd1$PGI_EF <- as.numeric( is_east & !is_male)
  if ("PGI_EM" %in% names(nd1)) nd1$PGI_EM <- as.numeric( is_east &  is_male)

  X0 <- predict(mod, newdata = nd0, type = "lpmatrix")
  X1 <- predict(mod, newdata = nd1, type = "lpmatrix")
  D  <- X1 - X0

  b   <- coef(mod); V <- vcov(mod)
  fit <- as.vector(D %*% b)
  se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

  data.frame(
    birth_year = BYc_grid + by_mean,
    BYc        = BYc_grid,
    east_west     = east_west,
    gender     = gender,
    slope      = fit,
    se         = se,
    ci_lo      = fit - 1.96 * se,
    ci_hi      = fit + 1.96 * se
  )
}

# ---- East-West difference smooth with simultaneous 95% CI -----------

diff_smooth_simul <- function(mod, dat, gender = "female",
                              pgi_var = "PGI_Edu",
                              n_grid = 200, n_sim = N_SIM, seed = 42) {
  BYc_range <- range(dat$BYc, na.rm = TRUE)
  BYc_grid  <- seq(BYc_range[1], BYc_range[2], length.out = n_grid)
  by_mean   <- attr(dat, "by_mean")
  if (is.null(by_mean)) by_mean <- mean(dat$birth_year, na.rm = TRUE)

  is_male <- gender == "male"

  by_sd_local <- sd(dat$birth_year, na.rm = TRUE)
  make_nd <- function(east_west, pgi_val) {
    is_east <- east_west == "East"
    nd <- data.frame(
      BYc          = BYc_grid,
      BYc_z        = if (is.finite(by_sd_local) && by_sd_local > 0)
                       BYc_grid / by_sd_local else BYc_grid,
      PGI_Edu      = pgi_val, PGI_Edu_z = pgi_val,
      PGI_Cog      = 0,       PGI_Cog_z = 0,
      PGI_nonCog   = 0,       PGI_nonCog_z = 0,
      PGI_Height   = 0,       PGI_Height_z = 0,
      stringsAsFactors = FALSE
    )
    if ("east_west" %in% names(mod$model))
      nd$east_west <- factor(east_west, levels = levels(dat$east_west))
    if ("gender" %in% names(mod$model))
      nd$gender <- factor(gender, levels = levels(dat$gender))
    if ("east_west_c" %in% names(mod$model))
      nd$east_west_c <- if (is_east) 0.5 else -0.5
    if ("gender_c" %in% names(mod$model))
      nd$gender_c <- if (is_male) -0.5 else 0.5
    if ("east"      %in% names(mod$model)) nd$east      <- as.integer(is_east)
    if ("gender_01" %in% names(mod$model)) nd$gender_01 <- as.integer(is_male)
    if ("parental_edu_z_kernel" %in% names(mod$model)) nd$parental_edu_z_kernel <- 0
    if ("ParEdu_z"  %in% names(mod$model)) nd$ParEdu_z  <- 0
    if ("fid_re" %in% names(mod$model)) {
      fid_levels <- levels(mod$model$fid_re)
      nd$fid_re <- factor(fid_levels[1], levels = fid_levels)
    }

    for (v in c("PGI_W","PGI_E","PGI_WF","PGI_WM","PGI_EF","PGI_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0
    # Region/gender-specific ParEdu control cells (B3 by-smooths) – held at 0
    for (v in c("ParEdu_W","ParEdu_E","ParEdu_WF","ParEdu_WM","ParEdu_EF","ParEdu_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0

    if ("PGI_W"  %in% names(nd)) nd$PGI_W  <- pgi_val * as.numeric(!is_east)
    if ("PGI_E"  %in% names(nd)) nd$PGI_E  <- pgi_val * as.numeric( is_east)
    if ("PGI_WF" %in% names(nd)) nd$PGI_WF <- pgi_val * as.numeric(!is_east & !is_male)
    if ("PGI_WM" %in% names(nd)) nd$PGI_WM <- pgi_val * as.numeric(!is_east &  is_male)
    if ("PGI_EF" %in% names(nd)) nd$PGI_EF <- pgi_val * as.numeric( is_east & !is_male)
    if ("PGI_EM" %in% names(nd)) nd$PGI_EM <- pgi_val * as.numeric( is_east &  is_male)
    nd
  }

  X_E1 <- predict(mod, newdata = make_nd("East", 1), type = "lpmatrix")
  X_E0 <- predict(mod, newdata = make_nd("East", 0), type = "lpmatrix")
  X_W1 <- predict(mod, newdata = make_nd("West", 1), type = "lpmatrix")
  X_W0 <- predict(mod, newdata = make_nd("West", 0), type = "lpmatrix")
  D <- (X_E1 - X_E0) - (X_W1 - X_W0)

  b <- coef(mod); V <- vcov(mod)
  fit <- as.vector(D %*% b)
  se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

  set.seed(seed)
  L <- tryCatch(chol(V), error = function(e) {
    Ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
    t(Ee$vectors %*% diag(sqrt(pmax(Ee$values, 0))) %*% t(Ee$vectors))
  })
  Zmax <- replicate(n_sim, {
    z <- rnorm(ncol(V))
    f_draw <- as.vector(D %*% (b + L %*% z))
    max(abs((f_draw - fit) / se), na.rm = TRUE)
  })
  crit <- as.numeric(quantile(Zmax, 0.95, na.rm = TRUE))

  data.frame(
    birth_year = BYc_grid + by_mean,
    BYc = BYc_grid,
    diff  = fit, se = se,
    ci_lo = fit - crit * se,
    ci_hi = fit + crit * se
  )
}

# Female-male difference smooth with simultaneous 95 % CI, holding region
# fixed. Mirrors diff_smooth_simul but contrasts gender within a region
# rather than region within a gender.
diff_smooth_gender_simul <- function(mod, dat, east_west,
                                     pgi_var = "PGI_Edu",
                                     n_grid = 200, n_sim = N_SIM, seed = 42) {
  BYc_range <- range(dat$BYc, na.rm = TRUE)
  BYc_grid  <- seq(BYc_range[1], BYc_range[2], length.out = n_grid)
  by_mean   <- attr(dat, "by_mean")
  if (is.null(by_mean)) by_mean <- mean(dat$birth_year, na.rm = TRUE)

  is_east <- east_west == "East"
  by_sd_local <- sd(dat$birth_year, na.rm = TRUE)

  make_nd <- function(gender, pgi_val) {
    is_male <- gender == "male"
    nd <- data.frame(
      BYc          = BYc_grid,
      BYc_z        = if (is.finite(by_sd_local) && by_sd_local > 0)
                       BYc_grid / by_sd_local else BYc_grid,
      PGI_Edu      = pgi_val, PGI_Edu_z = pgi_val,
      PGI_Cog      = 0,       PGI_Cog_z = 0,
      PGI_nonCog   = 0,       PGI_nonCog_z = 0,
      PGI_Height   = 0,       PGI_Height_z = 0,
      stringsAsFactors = FALSE
    )
    if ("east_west" %in% names(mod$model))
      nd$east_west <- factor(east_west, levels = levels(dat$east_west))
    if ("gender" %in% names(mod$model))
      nd$gender <- factor(gender, levels = levels(dat$gender))
    if ("east_west_c" %in% names(mod$model))
      nd$east_west_c <- if (is_east) 0.5 else -0.5
    if ("gender_c" %in% names(mod$model))
      nd$gender_c <- if (is_male) -0.5 else 0.5
    if ("east"      %in% names(mod$model)) nd$east      <- as.integer(is_east)
    if ("gender_01" %in% names(mod$model)) nd$gender_01 <- as.integer(is_male)
    if ("parental_edu_z_kernel" %in% names(mod$model)) nd$parental_edu_z_kernel <- 0
    if ("ParEdu_z"  %in% names(mod$model)) nd$ParEdu_z  <- 0
    if ("fid_re" %in% names(mod$model)) {
      fid_levels <- levels(mod$model$fid_re)
      nd$fid_re <- factor(fid_levels[1], levels = fid_levels)
    }

    for (v in c("PGI_W","PGI_E","PGI_WF","PGI_WM","PGI_EF","PGI_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0
    # Region/gender-specific ParEdu control cells (B3 by-smooths) – held at 0
    for (v in c("ParEdu_W","ParEdu_E","ParEdu_WF","ParEdu_WM","ParEdu_EF","ParEdu_EM"))
      if (v %in% names(mod$model)) nd[[v]] <- 0

    if ("PGI_W"  %in% names(nd)) nd$PGI_W  <- pgi_val * as.numeric(!is_east)
    if ("PGI_E"  %in% names(nd)) nd$PGI_E  <- pgi_val * as.numeric( is_east)
    if ("PGI_WF" %in% names(nd)) nd$PGI_WF <- pgi_val * as.numeric(!is_east & !is_male)
    if ("PGI_WM" %in% names(nd)) nd$PGI_WM <- pgi_val * as.numeric(!is_east &  is_male)
    if ("PGI_EF" %in% names(nd)) nd$PGI_EF <- pgi_val * as.numeric( is_east & !is_male)
    if ("PGI_EM" %in% names(nd)) nd$PGI_EM <- pgi_val * as.numeric( is_east &  is_male)
    nd
  }

  X_F1 <- predict(mod, newdata = make_nd("female", 1), type = "lpmatrix")
  X_F0 <- predict(mod, newdata = make_nd("female", 0), type = "lpmatrix")
  X_M1 <- predict(mod, newdata = make_nd("male",   1), type = "lpmatrix")
  X_M0 <- predict(mod, newdata = make_nd("male",   0), type = "lpmatrix")
  D <- (X_F1 - X_F0) - (X_M1 - X_M0)

  b <- coef(mod); V <- vcov(mod)
  fit <- as.vector(D %*% b)
  se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

  set.seed(seed)
  L <- tryCatch(chol(V), error = function(e) {
    Ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
    t(Ee$vectors %*% diag(sqrt(pmax(Ee$values, 0))) %*% t(Ee$vectors))
  })
  Zmax <- replicate(n_sim, {
    z <- rnorm(ncol(V))
    f_draw <- as.vector(D %*% (b + L %*% z))
    max(abs((f_draw - fit) / se), na.rm = TRUE)
  })
  crit <- as.numeric(quantile(Zmax, 0.95, na.rm = TRUE))

  data.frame(
    birth_year = BYc_grid + by_mean,
    BYc = BYc_grid,
    diff  = fit, se = se,
    ci_lo = fit - crit * se,
    ci_hi = fit + crit * se
  )
}

# Female - male difference smooth for the per-gender GAMs that pool region
# and use PGI_F / PGI_M cells (A5 / B5). Mirrors diff_smooth_gender_simul but
# for the pooled-region parameterisation (no PGI_WF/... cells). Returns the
# female-minus-male PGI-slope difference over birth year with a simultaneous
# 95% band. ParEdu_F / ParEdu_M (B5) are held at 0 and cancel in the contrast.
diff_smooth_pooled_gender_simul <- function(mod, dat, n_grid = 200,
                                            n_sim = N_SIM, seed = 42) {
  BYc_range <- range(dat$BYc, na.rm = TRUE)
  BYc_grid  <- seq(BYc_range[1], BYc_range[2], length.out = n_grid)
  by_mean   <- attr(dat, "by_mean")
  if (is.null(by_mean)) by_mean <- mean(dat$birth_year, na.rm = TRUE)

  make_nd <- function(is_female, pgi_val) {
    nd <- data.frame(BYc = BYc_grid, east_west_c = 0, gender_c = 0,
                     PGI_F = 0, PGI_M = 0, stringsAsFactors = FALSE)
    for (v in c("parental_edu_z_kernel", "ParEdu_F", "ParEdu_M"))
      if (v %in% names(mod$model)) nd[[v]] <- 0
    if ("fid_re" %in% names(mod$model)) {
      fl <- levels(mod$model$fid_re); nd$fid_re <- factor(fl[1], levels = fl)
    }
    if (is_female) nd$PGI_F <- pgi_val else nd$PGI_M <- pgi_val
    nd
  }
  X_F1 <- predict(mod, make_nd(TRUE,  1), type = "lpmatrix")
  X_F0 <- predict(mod, make_nd(TRUE,  0), type = "lpmatrix")
  X_M1 <- predict(mod, make_nd(FALSE, 1), type = "lpmatrix")
  X_M0 <- predict(mod, make_nd(FALSE, 0), type = "lpmatrix")
  D <- (X_F1 - X_F0) - (X_M1 - X_M0)

  b <- coef(mod); V <- vcov(mod)
  fit <- as.vector(D %*% b)
  se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

  set.seed(seed)
  L <- tryCatch(chol(V), error = function(e) {
    Ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
    t(Ee$vectors %*% diag(sqrt(pmax(Ee$values, 0))) %*% t(Ee$vectors))
  })
  Zmax <- replicate(n_sim, {
    z <- rnorm(ncol(V))
    f_draw <- as.vector(D %*% (b + L %*% z))
    max(abs((f_draw - fit) / se), na.rm = TRUE)
  })
  crit <- as.numeric(quantile(Zmax, 0.95, na.rm = TRUE))

  data.frame(
    birth_year = BYc_grid + by_mean, BYc = BYc_grid,
    diff = fit, se = se, ci_lo = fit - crit * se, ci_hi = fit + crit * se)
}

# ---- Cluster-robust helpers -----------------------------------------

cluster_vcov <- function(mod, cluster) {
  sandwich::vcovCL(mod, cluster = cluster, type = "HC1")
}

tidy_coefs_cluster <- function(mod, label, cluster_vec) {
  V  <- cluster_vcov(mod, cluster_vec)
  ct <- lmtest::coeftest(mod, vcov. = V)
  ci <- lmtest::coefci(mod, vcov. = V)
  out <- data.frame(
    model     = label,
    term      = rownames(ct),
    estimate  = ct[, "Estimate"],
    std.error = ct[, "Std. Error"],
    statistic = ct[, "t value"],
    p.value   = ct[, "Pr(>|t|)"],
    conf.low  = ci[, 1],
    conf.high = ci[, 2],
    vcov_type = "cluster-robust (HC1, by cohort_fid)",
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

# Returns (estimate, se_naive, p_naive, se_clust, p_clust) for a focal term.
focal_clust_vs_naive <- function(mod, term, cluster_vec) {
  cs <- coef(summary(mod))
  V  <- cluster_vcov(mod, cluster_vec)
  ct <- lmtest::coeftest(mod, vcov. = V)
  data.frame(
    term       = term,
    estimate   = if (term %in% rownames(cs)) cs[term, "Estimate"] else NA_real_,
    se_naive   = if (term %in% rownames(cs)) cs[term, "Std. Error"] else NA_real_,
    p_naive    = if (term %in% rownames(cs)) cs[term, "Pr(>|t|)"] else NA_real_,
    se_cluster = if (term %in% rownames(ct)) ct[term, "Std. Error"] else NA_real_,
    p_cluster  = if (term %in% rownames(ct)) ct[term, "Pr(>|t|)"] else NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---- Within-study slope curves ------------------------------------------
# Fit a varying-coefficient GAM separately in each contributing study and
# return per-study slope curves over birth year, optionally stratified by
# region and/or gender. For TwinLife, an fid random effect is added
# automatically.

within_cohort_slope_curves <- function(dat, outcome,
                                       pgi_var = "PGI_Edu",
                                       include_parental_edu = FALSE,
                                       gender_stratified    = FALSE,
                                       k = K_DEFAULT, bs = BS_DEFAULT,
                                       min_n = 200,
                                       min_cell_n = 50) {
  cohorts <- levels(droplevels(factor(dat$cohort)))
  all_curves <- list()
  fits <- list()

  for (coh in cohorts) {
    sub <- dat[as.character(dat$cohort) == coh, , drop = FALSE]
    sub <- sub[!is.na(sub[[outcome]]) & !is.na(sub[[pgi_var]]) & !is.na(sub$birth_year), ]
    if (nrow(sub) < min_n) next

    attr(sub, "by_mean") <- attr(dat, "by_mean")

    has_east  <- length(unique(na.omit(sub$east_west))) > 1
    has_male  <- length(unique(na.omit(sub$gender)))   > 1
    use_fid_re <- coh == "TwinLife" && "fid" %in% names(sub) &&
                  length(unique(sub$fid)) >= 5 &&
                  length(unique(sub$fid)) < nrow(sub)

    pgi_vec <- as.numeric(sub[[pgi_var]])
    east_i  <- as.integer(as.character(sub$east_west) == "East")
    male_i  <- as.integer(as.character(sub$gender)    == "male")

    # Region- (and optionally gender-) specific PGI columns so the PGI slope
    # is free to vary by cell, mirroring the main A3/B3/A6/B6 parameterization.
    if (gender_stratified && has_male && has_east) {
      sub$PGI_WF <- pgi_vec * (1L - east_i) * (1L - male_i)
      sub$PGI_WM <- pgi_vec * (1L - east_i) *        male_i
      sub$PGI_EF <- pgi_vec *        east_i * (1L - male_i)
      sub$PGI_EM <- pgi_vec *        east_i *        male_i
      cells <- list(
        list(reg = "West", gn = "female", col = "PGI_WF"),
        list(reg = "West", gn = "male",   col = "PGI_WM"),
        list(reg = "East", gn = "female", col = "PGI_EF"),
        list(reg = "East", gn = "male",   col = "PGI_EM")
      )
    } else if (has_east) {
      sub$PGI_W <- pgi_vec * (1L - east_i)
      sub$PGI_E <- pgi_vec *        east_i
      cells <- list(
        list(reg = "West", gn = NA_character_, col = "PGI_W"),
        list(reg = "East", gn = NA_character_, col = "PGI_E")
      )
    } else {
      cells <- list(list(reg = NA_character_, gn = NA_character_, col = pgi_var))
    }

    cells <- Filter(function(c) {
      m_rows <- rep(TRUE, nrow(sub))
      if (!is.na(c$reg)) m_rows <- m_rows & as.character(sub$east_west) == c$reg
      if (!is.na(c$gn))  m_rows <- m_rows & as.character(sub$gender)    == c$gn
      sum(m_rows) >= min_cell_n
    }, cells)
    if (length(cells) == 0) next

    cell_cols <- unique(vapply(cells, function(c) c$col, character(1)))

    rhs_parts <- c("birth_year")
    if (has_east) rhs_parts <- c(rhs_parts, "east_west_c")
    if (has_male) rhs_parts <- c(rhs_parts, "gender_c")
    have_pared <- include_parental_edu && "parental_edu_z_kernel" %in% names(sub) &&
                  sum(!is.na(sub$parental_edu_z_kernel)) > min_n
    if (have_pared) rhs_parts <- c(rhs_parts, "parental_edu_z_kernel")
    rhs_parts <- c(rhs_parts, cell_cols,
                   paste0(cell_cols, ":birth_year"),
                   sprintf("s(birth_year, by = %s, k = %d, bs = '%s')",
                           cell_cols, k, bs))
    if (have_pared)
      rhs_parts <- c(rhs_parts,
                     sprintf("s(birth_year, by = parental_edu_z_kernel, k = %d, bs = '%s')", k, bs))
    if (use_fid_re) {
      sub$fid <- factor(sub$fid)
      rhs_parts <- c(rhs_parts, "s(fid, bs = 're')")
    }

    f <- as.formula(paste(outcome, "~", paste(rhs_parts, collapse = " + ")))
    m <- tryCatch(mgcv::gam(f, data = sub, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) next
    fits[[coh]] <- m

    n_grid <- 200
    by_min <- min(sub$birth_year); by_max <- max(sub$birth_year)
    by_grid <- seq(by_min, by_max, length.out = n_grid)

    for (c in cells) {
      nd_base <- data.frame(birth_year = by_grid, stringsAsFactors = FALSE)
      if (has_east) {
        cell_east <- if (is.na(c$reg)) "West" else c$reg
        nd_base$east_west_c <- if (cell_east == "East") 0.5 else -0.5
      }
      if (has_male) {
        cell_male <- !is.na(c$gn) && c$gn == "male"
        nd_base$gender_c <- if (cell_male) -0.5 else 0.5
      }
      if (include_parental_edu && "parental_edu_z_kernel" %in% names(m$model))
        nd_base$parental_edu_z_kernel <- 0
      if (use_fid_re) nd_base$fid <- factor(levels(sub$fid)[1], levels = levels(sub$fid))
      for (cc in cell_cols) nd_base[[cc]] <- 0

      nd0 <- nd_base
      nd1 <- nd_base; nd1[[c$col]] <- 1

      X0 <- predict(m, newdata = nd0, type = "lpmatrix")
      X1 <- predict(m, newdata = nd1, type = "lpmatrix")
      D  <- X1 - X0
      b  <- coef(m); V <- vcov(m)
      fit <- as.vector(D %*% b)
      se  <- sqrt(pmax(1e-12, rowSums((D %*% V) * D)))

      all_curves[[length(all_curves) + 1]] <- data.frame(
        cohort     = coh,
        east_west  = if (has_east) c$reg else NA_character_,
        gender     = if (has_male && gender_stratified) c$gn else NA_character_,
        birth_year = by_grid,
        slope      = fit,
        se         = se,
        ci_lo      = fit - 1.96 * se,
        ci_hi      = fit + 1.96 * se,
        stringsAsFactors = FALSE
      )
    }
  }

  list(curves = if (length(all_curves)) bind_rows(all_curves) else NULL,
       fits = fits)
}

# ---------------------------------------------------------------------------
# focal_dispersion_check(): per-focal heteroscedasticity robustness, gated on
# the focal being significant in the main random-effects analysis.
#
# For each candidate focal interaction we (a) only run the check if the focal
# was significant in the main RE model (per-focal gate via each spec's
# `main_p`), then (b) re-estimate the mean focal allowing the residual variance
# to follow the SAME structure as the focal (Gaussian location-scale,
# mgcv::gaulss): homoscedastic baseline (scale ~ 1) vs heteroscedastic
# (scale ~ scale_rhs). gaulss carries no family random effect, so the
# homo-vs-het comparison is RE-free; the RE p-value is used only for gating.
#
# specs: list of lists, each with $label, $mean_rhs, $focal, $scale_rhs, $main_p.
# Returns one tidy row per spec (non-significant focals flagged Tested = FALSE).
focal_dispersion_check <- function(dat, outcome, specs, alpha = 0.05) {
  .fp <- function(p) .fmt_p_vec(p)   # workflow-wide p convention
  rows <- lapply(specs, function(s) {
    base <- data.frame(
      Focal = s$label,
      `Main RE p` = if (is.null(s$main_p)) NA_character_ else .fp(s$main_p),
      Tested = FALSE, `Homosc. est` = NA_real_, `Homosc. p` = NA_character_,
      `Het est` = NA_real_, `Het p` = NA_character_, `Scale LRT p` = NA_character_,
      Conclusion = NA_character_, check.names = FALSE, stringsAsFactors = FALSE)
    if (is.null(s$main_p) || is.na(s$main_p) || s$main_p >= alpha) {
      base$Conclusion <- "Not significant in the main analysis; dispersion check not run."
      return(base)
    }
    meanf <- stats::as.formula(paste(outcome, "~", s$mean_rhs))
    m_homo <- tryCatch(mgcv::gam(list(meanf, ~ 1), family = mgcv::gaulss(), data = dat, method = "REML"),
                       error = function(e) NULL)
    m_het  <- tryCatch(mgcv::gam(list(meanf, stats::as.formula(paste("~", s$scale_rhs))),
                                 family = mgcv::gaulss(), data = dat, method = "REML"),
                       error = function(e) NULL)
    gp <- function(m) { if (is.null(m)) return(c(NA_real_, NA_real_))
      p <- summary(m)$p.table
      if (s$focal %in% rownames(p)) c(p[s$focal, 1], p[s$focal, ncol(p)]) else c(NA_real_, NA_real_) }
    vh <- gp(m_homo); ve <- gp(m_het)
    lrt_p <- if (!is.null(m_homo) && !is.null(m_het))
      tryCatch(stats::anova(m_homo, m_het, test = "Chisq")[2, "Pr(>Chi)"], error = function(e) NA_real_) else NA_real_
    base$Tested <- TRUE
    base$`Homosc. est` <- round(vh[1], 4); base$`Homosc. p` <- .fp(vh[2])
    base$`Het est`     <- round(ve[1], 4); base$`Het p`     <- .fp(ve[2])
    base$`Scale LRT p` <- .fp(lrt_p)
    base$Conclusion <- if (is.na(ve[2])) "Location-scale model failed to fit." else
      if (ve[2] < alpha) "Robust: focal remains significant under heteroscedasticity." else
      "Attenuated: focal drops below alpha once dispersion is modelled."
    base
  })
  do.call(rbind, rows)
}
