# scale_and_censoring.R
# Sensitivity of two focal terms to the outcome's scale and to censoring:
#
#   attainment  PGI-Education x birth year   (A1_RE, term PGI_Edu_z:BYc_z)
#   mobility    PGI-Education x region       (B0_RE, term PGI_Edu_z:east_west_c)
#
# Each row is the published focal model with one thing changed, reported beside
# the published estimate and the row above it.
#
#   row 0   no model: ceiling (>= 18 years) and floor (<= 9 years) shares in the
#           mobility subsample and in the attainment sample (where they must
#           reproduce Table 1), by region, study and birth-year group; parental
#           education likewise.
#   row 1   outcome on raw years of education (attainment) and raw years above
#           the parents (mobility). A diagnostic; the standardised outcome stays
#           focal.
#   row 2   outcome replaced by its rank-based inverse normal transform, model
#           and covariates unchanged.
#   row 3   ordered probit on the raw education categories, same rhs as the
#           published model, no family random effect (as in the supplement's
#           gender check), SEs cluster-robust by family. Fitted by maximum
#           likelihood in base R so the cluster-robust sandwich is available;
#           agreement with MASS::polr is asserted in scale_censoring_checks.csv.
#
# Nothing here overwrites the frozen run and its model caches are not loaded
# (~650 MB each). The analytic samples are rebuilt once from the merged data
# through 00_setup/, as 05_manuscript_export/check_frozen_run.R step 4 does, and
# the script stops unless the published A1 and B0 focal coefficients reproduce.
# Reference estimates come from the manuscript repository's 01_results/. Models
# are released (rm + gc) between steps. Aggregate output only; no individual row
# is printed or written.
#
# Usage:
#   RUN_TS=<ts> Rscript 06_Sensitivity_analyses/scale_and_censoring.R
#   OUT_DIR=<dir>              write the CSVs elsewhere (default: the run's
#                              manuscript_export/sensitivity_analyses/)
#   MANUSCRIPT_REPO=<path>     where 01_results/ lives
#   DATA_ROOT=<path>           data root

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
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine scale_and_censoring.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
REPO_ROOT    <- normalizePath(file.path(WORKFLOW_DIR, "..", ".."), mustWork = TRUE)
rm(.find_this_file)


# ---------------------------------------------------------------------------
# Constants: the published specifications, verbatim from A_analysis.R / B_analysis.R
# ---------------------------------------------------------------------------
FROZEN_RUN_TS <- "20260729_1555"
CEIL_YEARS    <- 18      # Descriptives_analysis.R: education >= 18
FLOOR_YEARS   <- 9       # Descriptives_analysis.R: education <= 9
Z975          <- stats::qnorm(0.975)

A1_RHS  <- "(PGI_Edu_z + BYc_z + east_west_c + gender_c)^3"
B0_RHS  <- paste("PGI_Edu_z * east_west_c +",
                 "east_west_c * (gender_c + BYc_z + parental_edu_z_kernel) +",
                 "gender_c * BYc_z + gender_c * parental_edu_z_kernel +",
                 "BYc_z * parental_edu_z_kernel")
RE_TERM <- "s(fid_re, bs = 're')"
A_FOCAL <- "PGI_Edu_z:BYc_z"
B_FOCAL <- "PGI_Edu_z:east_west_c"
A_LABEL <- "attainment_pgi_x_birthyear"
B_LABEL <- "mobility_pgi_x_region"

# Row-1 mobility: the same B0 right-hand side with parental education in raw
# years, so outcome and covariate sit on one scale.
B0_RHS_RAW <- gsub("parental_edu_z_kernel", "parental_education", B0_RHS, fixed = TRUE)

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

# by_group3 - the study's three birth-year bins (helpers_descriptive_tables.R).
by_group3 <- function(by) {
  factor(ifelse(by < 1950, "< 1950", ifelse(by < 1975, "1950-1974", ">= 1975")),
         levels = c("< 1950", "1950-1974", ">= 1975"))
}

# Family cluster id: TwinLife families keep their id; everyone else is a singleton.
family_cluster <- function(d) {
  f <- as.character(d$fid_re)
  ifelse(!is.na(f) & f != "non_TL", f, paste0("i", seq_len(nrow(d))))
}

# Rank-based inverse normal transform (Blom offsets, mid-ranks for ties).
rank_int <- function(y) {
  r <- rank(y, ties.method = "average", na.last = "keep")
  stats::qnorm((r - 3 / 8) / (sum(!is.na(y)) + 1 / 4))
}

fit_bam <- function(lhs, rhs, data) {
  mgcv::bam(stats::as.formula(paste(lhs, "~", rhs, "+", RE_TERM)),
            data = data, method = "fREML", discrete = TRUE)
}

# One parametric coefficient from a bam/gam: estimate, SE, statistic, p from the
# p.table (as tidy_coefs() / broom read them), normal-based CI (as the frozen
# export carries: conf.low = estimate - qnorm(.975) * SE).
p_row <- function(m, term) {
  pt <- summary(m)$p.table
  if (!term %in% rownames(pt)) stop("term not in model: ", term)
  est <- unname(pt[term, 1]); se <- unname(pt[term, 2])
  data.frame(term = term, estimate = est, se = se, statistic = unname(pt[term, 3]),
             p = unname(pt[term, 4]), ci_lo = est - Z975 * se, ci_hi = est + Z975 * se,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Ordered probit, maximum likelihood in base R
# ---------------------------------------------------------------------------
# y in 1..K (ordered categories), X = design without intercept. Latent
# y* = X b + e, e ~ N(0, 1), thresholds tau_1 < ... < tau_{K-1} parameterised as
# tau_1 = zeta_1, tau_k = tau_{k-1} + exp(zeta_k) (MASS::polr's parameterisation).
# Same likelihood as MASS::polr(method = "probit"); written out so that
# per-observation scores are available for a cluster-robust sandwich covariance
# (sandwich::vcovCL's HC1 small-sample adjustment). Returns coefficients, the
# observed-information covariance (numerically differentiated analytic gradient),
# the cluster-robust covariance, log-likelihood and convergence diagnostics. The
# agreement with MASS::polr is asserted in scale_censoring_checks.csv.
fit_oprobit <- function(y, X, cluster = NULL, start = NULL, reltol = 1e-12, maxit = 5000) {
  y <- as.integer(y); K <- max(y)
  stopifnot(min(y) == 1L, K >= 3L, all(seq_len(K) %in% y))
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  nt <- K - 1L; P <- p + nt
  pnames <- c(colnames(X), paste0("zeta", seq_len(nt)))

  unpack <- function(par) {
    beta <- par[seq_len(p)]; zeta <- par[p + seq_len(nt)]
    list(beta = beta, zeta = zeta, tau = cumsum(c(zeta[1], exp(zeta[-1]))))
  }
  parts <- function(par) {
    u <- unpack(par)
    eta  <- drop(X %*% u$beta)
    u_hi <- c(u$tau, Inf)[y] - eta
    u_lo <- c(-Inf, u$tau)[y] - eta
    pr <- ifelse(u_lo > 0,
                 stats::pnorm(u_lo, lower.tail = FALSE) - stats::pnorm(u_hi, lower.tail = FALSE),
                 stats::pnorm(u_hi) - stats::pnorm(u_lo))
    list(u = u, u_hi = u_hi, u_lo = u_lo, pr = pmax(pr, 1e-300))
  }
  nll <- function(par) -sum(log(parts(par)$pr))
  scores <- function(par) {                       # n x P matrix of per-observation scores
    s <- parts(par); u <- s$u
    d_hi <- stats::dnorm(s$u_hi); d_lo <- stats::dnorm(s$u_lo)   # dnorm(+-Inf) = 0
    S_beta <- X * (-(d_hi - d_lo) / s$pr)
    A <- matrix(0, n, nt)                         # d loglik / d tau_k
    a_hi <-  d_hi / s$pr; a_lo <- -d_lo / s$pr
    i_hi <- which(y <= nt); A[cbind(i_hi, y[i_hi])] <- a_hi[i_hi]
    i_lo <- which(y >= 2L); A[cbind(i_lo, y[i_lo] - 1L)] <- A[cbind(i_lo, y[i_lo] - 1L)] + a_lo[i_lo]
    J <- matrix(0, nt, nt); J[, 1] <- 1           # d tau / d zeta
    if (nt >= 2) for (m in 2:nt) J[m:nt, m] <- exp(u$zeta[m])
    cbind(S_beta, A %*% J)
  }
  gr <- function(par) -colSums(scores(par))

  if (is.null(start)) {                           # marginal-quantile thresholds, zero slopes
    cp  <- cumsum(tabulate(y, K) / n)[seq_len(nt)]
    tau0 <- stats::qnorm(pmin(pmax(cp, 1e-6), 1 - 1e-6))
    start <- c(rep(0, p), tau0[1], log(pmax(diff(tau0), 1e-6)))
  }
  opt <- stats::optim(start, nll, gr, method = "BFGS",
                      control = list(maxit = maxit, reltol = reltol))
  par <- opt$par
  # Newton polish to a tight gradient (BFGS stops on the function value).
  H <- stats::optimHess(par, nll, gr)
  for (it in seq_len(25)) {
    g <- gr(par); if (sqrt(sum(g^2)) < 1e-6) break
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    cand <- par - step; f_old <- nll(par); f_new <- nll(cand); h <- 1
    while (!is.finite(f_new) || f_new > f_old + 1e-10) {      # backtrack
      h <- h / 2; if (h < 1e-6) break
      cand <- par - h * step; f_new <- nll(cand)
    }
    if (f_new > f_old + 1e-10) break
    par <- cand; H <- stats::optimHess(par, nll, gr)
  }
  g <- gr(par); grad_norm <- sqrt(sum(g^2))
  V_naive <- tryCatch(solve(H), error = function(e) matrix(NA_real_, P, P))
  S <- scores(par)
  if (is.null(cluster)) cluster <- seq_len(n)
  Sc <- rowsum(S, as.character(cluster)); G <- nrow(Sc)
  adj <- (G / (G - 1)) * ((n - 1) / (n - P))     # sandwich::vcovCL type = "HC1"
  V_cl <- adj * V_naive %*% crossprod(Sc) %*% V_naive
  dimnames(V_naive) <- dimnames(V_cl) <- list(pnames, pnames)
  names(par) <- pnames
  u <- unpack(par)
  list(coef = par, beta = setNames(u$beta, colnames(X)), tau = u$tau,
       vcov_naive = V_naive, vcov_cluster = V_cl,
       se_naive = sqrt(diag(V_naive)), se_cluster = sqrt(diag(V_cl)),
       logLik = -nll(par), n = n, K = K, n_clusters = G, P = P,
       converged = (opt$convergence == 0L) && grad_norm < 1e-3, grad_norm = grad_norm,
       optim_message = opt$message, nll = nll, gr = gr)
}

# Design matrix for the ordered probit from a model formula string (no intercept).
design_no_intercept <- function(rhs, data) {
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = data)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

# One tidy row for a focal term of a fit_oprobit() object, cluster-robust
# inference primary, model-based alongside.
oprobit_row <- function(fit, term) {
  est <- unname(fit$coef[term]); se_c <- unname(fit$se_cluster[term]); se_n <- unname(fit$se_naive[term])
  data.frame(term = term, estimate = est, se = se_c, statistic = est / se_c,
             p = 2 * stats::pnorm(-abs(est / se_c)),
             ci_lo = est - Z975 * se_c, ci_hi = est + Z975 * se_c,
             se_naive = se_n, p_naive = 2 * stats::pnorm(-abs(est / se_n)),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Row 0: the crosstab (no model)
# ---------------------------------------------------------------------------
censoring_crosstab <- function(d, sample_label) {
  one <- function(v, grp_type, grp, n_total, measure) {
    v <- v[!is.na(v)]
    data.frame(sample = sample_label, measure = measure, group_type = grp_type, group = grp,
               n = length(v), share_of_sample_pct = 100 * length(v) / n_total,
               pct_ceiling = 100 * mean(v >= CEIL_YEARS), pct_floor = 100 * mean(v <= FLOOR_YEARS),
               pct_censored = 100 * mean(v >= CEIL_YEARS | v <= FLOOR_YEARS),
               mean_years = mean(v), sd_years = stats::sd(v), stringsAsFactors = FALSE)
  }
  measures <- list(own_education = "education")
  if ("parental_education" %in% names(d) && any(!is.na(d$parental_education)))
    measures$parental_education <- "parental_education"
  rows <- list()
  for (mname in names(measures)) {
    v <- d[[measures[[mname]]]]; n_total <- sum(!is.na(v))
    rows[[length(rows) + 1]] <- one(v, "overall", "all", n_total, mname)
    for (r in c("West", "East"))
      rows[[length(rows) + 1]] <- one(v[d$east_west == r], "region", r, n_total, mname)
    for (s in sort(unique(as.character(d$cohort))))
      rows[[length(rows) + 1]] <- one(v[as.character(d$cohort) == s], "study", s, n_total, mname)
    bg <- by_group3(d$birth_year)
    for (b in levels(bg))
      rows[[length(rows) + 1]] <- one(v[bg == b], "birth_year_group", b, n_total, mname)
    for (r in c("West", "East")) for (b in levels(bg))
      rows[[length(rows) + 1]] <- one(v[d$east_west == r & bg == b], "region_x_birth_year_group",
                                      paste(r, b, sep = " | "), n_total, mname)
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# The ladder for one focal term
# ---------------------------------------------------------------------------
# spec: list(label, focal, rhs, rhs_raw, outcome_z, outcome_raw, main_term,
#            ref_file, ref_model, sample)
ladder_for <- function(d, spec, ref_row, commit, run_ts, checks) {
  cl <- family_cluster(d)
  n  <- nrow(d); G <- length(unique(cl))
  base <- function(row, factor_changed, specification, outcome, estimator, r,
                   se_type, main_est = NA_real_, note = "") {
    data.frame(focal = spec$label, row = as.character(row), factor_changed = factor_changed,
               specification = specification, outcome = outcome, estimator = estimator,
               sample = spec$sample, n = n, n_clusters = G,
               term = r$term, estimate = r$estimate, se = r$se, statistic = r$statistic, p = r$p,
               ci_lo = r$ci_lo, ci_hi = r$ci_hi, se_type = se_type,
               main_effect_estimate = main_est,
               ratio_focal_to_main = if (is.na(main_est)) NA_real_ else r$estimate / main_est,
               sign_matches_published = sign(r$estimate) == sign(ref_row$estimate),
               significant_05 = r$p < 0.05, note = note,
               source_commit = commit, run_ts = run_ts, stringsAsFactors = FALSE)
  }
  out <- list()

  # ---- reference row: the published estimate, read from the frozen export ----
  out$ref <- base("ref", "none (published focal model)",
                  paste0(spec$ref_model, ": ", spec$outcome_z, " ~ ", spec$rhs, " + ", RE_TERM),
                  spec$outcome_z, "linear, family random effect (mgcv::bam fREML)", ref_row,
                  "model-based (frozen export)", main_est = ref_row$main_estimate,
                  note = paste0("read from ", spec$ref_file))
  out$ref$source_commit <- ref_row$source_commit

  # refit of the published model on the rebuilt sample: must reproduce the frozen row
  m_ref <- fit_bam(spec$outcome_z, spec$rhs, d)
  rr <- p_row(m_ref, spec$focal)
  checks$add(paste0(spec$label, ": published focal reproduced on the rebuilt sample"),
             ref_row$estimate, rr$estimate, tol = 1e-6)
  checks$add(paste0(spec$label, ": published focal SE reproduced"), ref_row$se, rr$se, tol = 1e-6)
  rm(m_ref); invisible(gc())

  # ---- row 1: outcome scale (raw years) ----
  m1 <- fit_bam(spec$outcome_raw, spec$rhs_raw, d)
  out$row1 <- base(1L, "outcome scale: raw years instead of the birth-year-standardised score",
                   paste0(spec$outcome_raw, " ~ ", spec$rhs_raw, " + ", RE_TERM), spec$outcome_raw,
                   "linear, family random effect (mgcv::bam fREML)", p_row(m1, spec$focal),
                   "model-based", main_est = p_row(m1, spec$main_term)$estimate,
                   note = "diagnostic, not a replacement outcome; years per SD of PGI")
  rm(m1); invisible(gc())

  # ---- row 2: outcome transform, rank-based inverse normal transform ----
  d$.y_rint <- rank_int(d[[spec$outcome_z]])
  m2 <- fit_bam(".y_rint", spec$rhs, d)
  out$row2 <- base(2L, "outcome transform: rank-based inverse normal transform of the outcome",
                   paste0("rank_int(", spec$outcome_z, ") ~ ", spec$rhs, " + ", RE_TERM),
                   paste0("rank_int(", spec$outcome_z, ")"),
                   "linear, family random effect (mgcv::bam fREML)", p_row(m2, spec$focal),
                   "model-based", main_est = p_row(m2, spec$main_term)$estimate,
                   note = "Blom scores of the outcome; covariates as published")
  rm(m2); invisible(gc())

  # ---- row 3: estimator, ordered probit on raw education categories ----
  y_ord <- as.integer(factor(d$education, ordered = TRUE))
  X <- design_no_intercept(spec$rhs, d)
  f3 <- fit_oprobit(y_ord, X, cluster = cl)
  if (!f3$converged) warning(spec$label, ": row-3 ordered probit did not converge (grad ", f3$grad_norm, ")")
  # MASS::polr on the same design: the estimator check (same likelihood, independent code)
  dfX <- as.data.frame(X); names(dfX) <- make.names(colnames(X)); dfX$.y <- factor(y_ord, ordered = TRUE)
  pol <- tryCatch(MASS::polr(.y ~ ., data = dfX, method = "probit", Hess = TRUE), error = function(e) NULL)
  if (is.null(pol))
    pol <- tryCatch(MASS::polr(.y ~ ., data = dfX, method = "probit", Hess = TRUE,
                               start = c(f3$beta, f3$tau)), error = function(e) NULL)
  if (!is.null(pol)) {
    checks$add(paste0(spec$label, ": row 3 ML slopes reproduce MASS::polr (max abs diff)"),
               0, max(abs(unname(coef(pol)) - unname(f3$beta))), tol = 1e-4)
    checks$add(paste0(spec$label, ": row 3 ML thresholds reproduce MASS::polr (max abs diff)"),
               0, max(abs(unname(pol$zeta) - f3$tau)), tol = 1e-4)
    pse <- sqrt(diag(vcov(pol)))[make.names(spec$focal)]
    checks$add(paste0(spec$label, ": row 3 focal model-based SE vs MASS::polr (ratio)"),
               1, unname(f3$se_naive[spec$focal]) / unname(pse), tol = 0.01)
  } else {
    checks$add(paste0(spec$label, ": MASS::polr comparison fit"), 1, 0, tol = 0,
               note = "polr failed to fit; no comparison")
  }
  r3 <- oprobit_row(f3, spec$focal)
  out$row3 <- base(3L, "estimator: ordered probit on the raw education categories (no family RE)",
                   paste0("ordered probit: edu_ord ~ ", spec$rhs),
                   sprintf("education as %d ordered categories", f3$K),
                   "ordered probit, ML (base R), cluster-robust SE by family", r3,
                   "cluster-robust by family (HC1); model-based SE in se_naive",
                   main_est = unname(f3$coef[spec$main_term]),
                   note = sprintf("latent-scale units; model-based p = %.4g; log-lik %.2f; %d categories",
                                  r3$p_naive, f3$logLik, f3$K))
  out$row3$se_naive <- r3$se_naive; out$row3$p_naive <- r3$p_naive
  rm(f3, X, pol, dfX); invisible(gc())

  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# Checks collector
# ---------------------------------------------------------------------------
new_checks <- function() {
  e <- new.env(); e$rows <- list()
  e$add <- function(name, expected, observed, tol, note = "") {
    d <- abs(expected - observed)
    e$rows[[length(e$rows) + 1]] <- data.frame(
      check = name, expected = expected, observed = observed, abs_diff = d, tol = tol,
      pass = is.finite(d) && d <= tol, note = note, stringsAsFactors = FALSE)
    invisible(NULL)
  }
  e$table <- function() do.call(rbind, e$rows)
  e
}

# Parse "12.97 (2.74)" -> 12.97 (Table 1 / Table 2 cells).
msd_num <- function(x) as.numeric(sub("\\s.*$", "", as.character(x)))

# ---------------------------------------------------------------------------
# Samples: rebuilt once from the merge through 00_setup/ (check_frozen_run.R step 4)
# ---------------------------------------------------------------------------
rebuild_samples <- function() {
  suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(tibble) })
  RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()
  EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
  source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
  # prepare_samples.R calls safe_numeric()/cohort_groups() from the analysis helpers.
  source(file.path(WORKFLOW_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
  DATA_DIR <- file.path(Sys.getenv("DATA_ROOT"), "03_Merge")
  sink(tempfile())                                     # load/prepare are chatty by design
  on.exit(sink(), add = TRUE)
  suppressMessages(suppressWarnings({
    source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
    source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)
  }))
  sink(); on.exit()
  need_A <- c("education", "edu_z_kernel", "PGI_Edu_z", "BYc_z", "east_west_c", "gender_c",
              "east_west", "gender", "cohort", "fid_re", "birth_year")
  need_B <- c(need_A, "mobility", "parental_edu_z_kernel", "parental_education")
  datA <- as.data.frame(dat_cluster)[, need_A]
  datB <- as.data.frame(dat_cluster_mob)[, need_B]
  list(datA = datA, datB = datB, data_file = DATA_FILE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_ladder <- function(out_dir, results_dir, samples, run_ts) {
  # results_dir: the frozen exports (01_results/manuscript_export in the manuscript
  # repo); samples: list(datA, datB, data_file) from rebuild_samples()
  run_ts <- sub("^RUN_", "", run_ts)
  exp_dir <- results_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  commit <- source_commit(); dirty <- source_dirty()
  checks <- new_checks()
  t0 <- Sys.time()
  message(sprintf("A2 ladder: RUN_%s\n  exports  -> %s\n  commit      %s%s", run_ts, out_dir,
                  substr(commit, 1, 8), if (dirty) " (working tree dirty)" else ""))

  # frozen references (never typed)
  prov <- utils::read.csv(file.path(exp_dir, "_provenance.csv"), stringsAsFactors = FALSE)
  frozen_commit <- as.character(prov$analysis_repo_commit[1])
  A1c <- utils::read.csv(file.path(exp_dir, "attainment", "main", "A1_Coefficients.csv"), stringsAsFactors = FALSE)
  B0c <- utils::read.csv(file.path(exp_dir, "mobility",   "main", "B0_Coefficients.csv"), stringsAsFactors = FALSE)
  t1  <- utils::read.csv(file.path(exp_dir, "descriptives_table1_overall.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  sf  <- utils::read.csv(file.path(exp_dir, "figS06_gender_pgi_scalefree.csv"), stringsAsFactors = FALSE)
  ref_from <- function(cf, term, main_term, file) {
    i <- which(cf$term == term)[1]; j <- which(cf$term == main_term)[1]
    if (is.na(i)) stop("frozen term not found: ", term, " in ", file)
    data.frame(term = term, estimate = cf$estimate[i], se = cf$std.error[i], statistic = cf$statistic[i],
               p = cf$p.value[i], ci_lo = cf$conf.low[i], ci_hi = cf$conf.high[i],
               main_estimate = cf$estimate[j], source_commit = frozen_commit, stringsAsFactors = FALSE)
  }
  refA <- ref_from(A1c, A_FOCAL, "PGI_Edu_z", "attainment/main/A1_Coefficients.csv")
  refB <- ref_from(B0c, B_FOCAL, "PGI_Edu_z", "mobility/main/B0_Coefficients.csv")

  # samples (rebuilt once by the caller; the merge is not re-read here)
  s <- samples; datA <- s$datA; datB <- s$datB; rm(samples)
  datA$east_west <- factor(as.character(datA$east_west), levels = c("West", "East"))
  datB$east_west <- factor(as.character(datB$east_west), levels = c("West", "East"))
  datA$fid_re <- factor(as.character(datA$fid_re)); datB$fid_re <- factor(as.character(datB$fid_re))
  message(sprintf("  samples: attainment N = %d, mobility N = %d", nrow(datA), nrow(datB)))
  t1n <- function(measure, col) msd_num(t1[[col]][t1$Measure == measure])
  checks$add("attainment sample N equals Table 1", t1$N_overall[t1$Measure == "Years of education"], nrow(datA), 0)
  checks$add("mobility sample N equals Table 1",   t1$N_overall[t1$Measure == "Mobility (kernel-z diff)"], nrow(datB), 0)
  if (!all(checks$table()$pass)) stop("Sample sizes differ from the frozen Table 1 - wrong run or wrong sample source.")

  # ---- row 0 ----
  message("row 0: censoring crosstab ...")
  ct <- rbind(censoring_crosstab(datA, "attainment (Table 1 sample)"),
              censoring_crosstab(datB, "mobility subsample (B0 sample)"))
  ct$source_commit <- commit; ct$run_ts <- run_ts
  g <- function(sample, grp, col) ct[[col]][ct$sample == sample & ct$measure == "own_education" &
                                              ct$group_type %in% c("overall", "region") & ct$group == grp][1]
  att <- "attainment (Table 1 sample)"
  checks$add("Table 1 % at ceiling, overall, reproduced", t1n("% at ceiling (edu = 18y)", "MSD_overall"), round(g(att, "all", "pct_ceiling"), 2), 0.005)
  checks$add("Table 1 % at ceiling, East, reproduced",    t1n("% at ceiling (edu = 18y)", "MSD_East"),    round(g(att, "East", "pct_ceiling"), 2), 0.005)
  checks$add("Table 1 % at ceiling, West, reproduced",    t1n("% at ceiling (edu = 18y)", "MSD_West"),    round(g(att, "West", "pct_ceiling"), 2), 0.005)
  checks$add("Table 1 % at floor, overall, reproduced",   t1n("% at floor (edu <= 9y)", "MSD_overall"),   round(g(att, "all", "pct_floor"), 2), 0.005)
  checks$add("Table 1 % at floor, East, reproduced",      t1n("% at floor (edu <= 9y)", "MSD_East"),      round(g(att, "East", "pct_floor"), 2), 0.005)
  checks$add("Table 1 % at floor, West, reproduced",      t1n("% at floor (edu <= 9y)", "MSD_West"),      round(g(att, "West", "pct_floor"), 2), 0.005)
  checks$add("Table 1 mobility N East reproduced", t1$N_East[t1$Measure == "Mobility (kernel-z diff)"],
             ct$n[ct$sample == "mobility subsample (B0 sample)" & ct$measure == "own_education" & ct$group == "East"][1], 0)

  # ---- the supplement's gender ordered probit, reproduced on the rebuilt sample ----
  # (proof that the sample and the estimator are the ones the supplement reports)
  message("reproducing the supplement's gender ordered-probit check ...")
  dd <- datA[!is.na(datA$PGI_Edu_z) & !is.na(datA$education) & !is.na(datA$edu_z_kernel) & !is.na(datA$gender), ]
  dd$edu_ord <- factor(dd$education, ordered = TRUE); dd$gender_c <- ifelse(dd$gender == "male", 0.5, -0.5)
  sfv <- function(k) sf$value[sf$key == k]
  om <- tryCatch(MASS::polr(edu_ord ~ PGI_Edu_z * gender_c + BYc_z + east_west_c, data = dd,
                            method = "probit", Hess = TRUE), error = function(e) NULL)
  if (!is.null(om)) {
    checks$add("figS06 ordered-probit PGI x gender b reproduced (MASS::polr)", sfv("ord_b"),
               unname(coef(om)["PGI_Edu_z:gender_c"]), 1e-6)
    Xg <- design_no_intercept("PGI_Edu_z * gender_c + BYc_z + east_west_c", dd)
    fg <- fit_oprobit(as.integer(dd$edu_ord), Xg)
    checks$add("figS06 ordered-probit PGI x gender b reproduced (base-R ML)",
               sfv("ord_b"), unname(fg$beta["PGI_Edu_z:gender_c"]), 1e-4)
    rm(Xg, fg)
  }
  rm(dd, om); invisible(gc())

  # ---- levels vs difference: the note's algebra, on the data ----
  # Exact under OLS (same X, outcome shifted by a column of X). Under the fREML
  # random-effects fit the smoothing parameter is re-estimated numerically, so the
  # two fits agree only to the optimiser's tolerance (observed ~5e-5 on 0.073).
  message("checking the levels/difference equivalence for the mobility model ...")
  m_dif_lm <- stats::lm(stats::as.formula(paste("mobility ~", B0_RHS)), data = datB)
  m_lev_lm <- stats::lm(stats::as.formula(paste("edu_z_kernel ~", B0_RHS)), data = datB)
  checks$add("mobility (OLS, exact algebra): PGI x region identical when Y* replaces M = Y* - P* as outcome",
             unname(coef(m_dif_lm)[B_FOCAL]), unname(coef(m_lev_lm)[B_FOCAL]), 1e-9)
  checks$add("mobility (OLS, exact algebra): P* coefficient shifts by exactly +1 between the two",
             1, unname(coef(m_lev_lm)["parental_edu_z_kernel"] - coef(m_dif_lm)["parental_edu_z_kernel"]), 1e-9)
  rm(m_dif_lm, m_lev_lm)
  m_lev <- fit_bam("edu_z_kernel", B0_RHS, datB)
  checks$add("mobility (bam fREML, optimiser tolerance): PGI x region identical when Y* replaces M as outcome",
             refB$estimate, p_row(m_lev, B_FOCAL)$estimate, 1e-4)
  checks$add("mobility (bam fREML, optimiser tolerance): P* coefficient shifts by +1 between the two",
             1, p_row(m_lev, "parental_edu_z_kernel")$estimate -
               B0c$estimate[B0c$term == "parental_edu_z_kernel"][1], 1e-4)
  rm(m_lev); invisible(gc())

  # ---- the two ladders ----
  specA <- list(label = A_LABEL, focal = A_FOCAL, rhs = A1_RHS, rhs_raw = A1_RHS,
                outcome_z = "edu_z_kernel", outcome_raw = "education",
                main_term = "PGI_Edu_z", ref_file = "attainment/main/A1_Coefficients.csv", ref_model = "A1_RE",
                sample = "attainment, full-family RE sample (dat_cluster)")
  specB <- list(label = B_LABEL, focal = B_FOCAL, rhs = B0_RHS, rhs_raw = B0_RHS_RAW,
                outcome_z = "mobility", outcome_raw = ".mob_raw",
                main_term = "PGI_Edu_z", ref_file = "mobility/main/B0_Coefficients.csv", ref_model = "B0_RE",
                sample = "mobility, full-family RE subsample (dat_cluster_mob)")
  datB$.mob_raw <- datB$education - datB$parental_education

  message("ladder: attainment, PGI x birth year ...")
  LA <- ladder_for(datA, specA, refA, commit, run_ts, checks)
  message("ladder: mobility, PGI x region ...")
  LB <- ladder_for(datB, specB, refB, commit, run_ts, checks)

  ladder <- dplyr::bind_rows(LA, LB)
  ladder$outcome[ladder$outcome == ".mob_raw"] <- "education - parental_education (raw years)"
  ladder$specification <- gsub(".mob_raw", "(education - parental_education)", ladder$specification, fixed = TRUE)

  # ---- write ----
  w <- function(df, name) {
    path <- file.path(out_dir, name)
    utils::write.csv(df, path, row.names = FALSE)
    message(sprintf("  [ok]   %-38s %d x %d", name, nrow(df), ncol(df)))
  }
  w(ladder, "scale_censoring_ladder.csv")
  w(ct,     "censoring_crosstab.csv")
  chk <- checks$table(); chk$source_commit <- commit; chk$run_ts <- run_ts
  w(chk, "scale_censoring_checks.csv")
  provo <- data.frame(
    run_ts = run_ts, frozen_analysis_repo_commit = frozen_commit, script_commit = commit,
    script_working_tree_dirty = dirty, reference_source = exp_dir, data_source = s$data_file,
    n_attainment = nrow(datA), n_mobility = nrow(datB),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R_version = R.version.string,
    mgcv_version = as.character(utils::packageVersion("mgcv")),
    MASS_version = as.character(utils::packageVersion("MASS")),
    checks_total = nrow(chk), checks_passed = sum(chk$pass),
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    stringsAsFactors = FALSE)
  w(provo, "scale_censoring_provenance.csv")

  # ---- console summary (aggregates only) ----
  cat("\n== checks ==\n")
  for (i in seq_len(nrow(chk)))
    cat(sprintf("  [%s] %-95s diff %.3g (tol %.3g)\n", if (chk$pass[i]) "ok" else "FAIL",
                chk$check[i], chk$abs_diff[i], chk$tol[i]))
  cat("\n== ladder ==\n")
  for (i in seq_len(nrow(ladder)))
    cat(sprintf("  %-28s row %-4s %-45s b = %+.4f  SE %.4f  95%% CI [%+.4f, %+.4f]  p = %.4g  %s\n",
                ladder$focal[i], ladder$row[i], substr(ladder$factor_changed[i], 1, 45),
                ladder$estimate[i], ladder$se[i], ladder$ci_lo[i], ladder$ci_hi[i], ladder$p[i],
                if (ladder$sign_matches_published[i]) "same sign" else "SIGN FLIP"))
  cat("\n== row 0: regional censoring, own education ==\n")
  r0 <- ct[ct$measure == "own_education" & ct$group_type == "region", ]
  for (i in seq_len(nrow(r0)))
    cat(sprintf("  %-32s %-5s n = %5d (%.1f%%)  ceiling %.2f%%  floor %.2f%%  censored %.2f%%\n",
                r0$sample[i], r0$group[i], r0$n[i], r0$share_of_sample_pct[i],
                r0$pct_ceiling[i], r0$pct_floor[i], r0$pct_censored[i]))
  if (!all(chk$pass)) cat("\nWARNING: some checks failed - read scale_censoring_checks.csv before using the exports.\n")
  cat(sprintf("\nDone in %.1f min.\n", provo$minutes))
  invisible(list(ladder = ladder, crosstab = ct, checks = chk, provenance = provo))
}

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
if (!nzchar(Sys.getenv("RUN_TS"))) Sys.setenv(RUN_TS = FROZEN_RUN_TS)
RUN_TS  <- get_run_ts()
if (!nzchar(Sys.getenv("DATA_ROOT"))) Sys.setenv(DATA_ROOT = data_root())
RUN_DIR <- file.path(gxe_dir(), "01_PGI_Analysis", "runs", paste0("RUN_", RUN_TS))
if (!dir.exists(RUN_DIR)) stop("Run folder not reachable: ", RUN_DIR)
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(RUN_DIR, "manuscript_export", "sensitivity_analyses"))
# Published references come from the manuscript repo's 01_results/ (the run's own
# manuscript_export/ holds the same files; used only if the repo is not here).
MS <- Sys.getenv("MANUSCRIPT_REPO", unset = file.path(
  path.expand("~"), "Documents", "017_manuscript", "80_years_GxE_on_Education_in_Germany-manuscript"))
RESULTS_DIR <- file.path(MS, "01_results", "manuscript_export")
if (!file.exists(file.path(RESULTS_DIR, "_provenance.csv"))) {
  message("01_results/ not found at ", MS, " - reading the frozen exports from the run folder")
  RESULTS_DIR <- file.path(RUN_DIR, "manuscript_export")
}
message("Rebuilding the analytic samples from the merge (00_setup/) ...")
SAMPLES <- rebuild_samples(); invisible(gc())
invisible(run_ladder(OUT_DIR, RESULTS_DIR, SAMPLES, run_ts = RUN_TS))
