# B_paredu_moderation.R
#
# Exploratory (post-plan). Tests the Morris et al. (2026, PNAS) G x E result —
# amplification of polygenic prediction by family socioeconomic advantage — in the
# German mobility sample, and asks whether it survives completing the Keller (2014)
# covariate-interaction set.
#
# Keller (2014) requires, for a G x E test with covariates C, that the model carry
# G x C and E x C for every C. Here G = PGI_Edu_z, E = parental_edu_z_kernel, and
# C = {BYc_z, east_west_c, gender_c}. Auditing the published B1 right-hand side:
#
#   G x C : PGI x BYc_z          present (^3 block)
#           PGI x east_west_c    present (^3 block)
#           PGI x gender_c       present (^3 block)
#   E x C : ParEdu x BYc_z       present
#           ParEdu x east_west_c present
#           ParEdu x gender_c    MISSING  <- the only genuinely absent term
#
# Two specifications are therefore fitted:
#   published : the B1 right-hand side exactly as published.
#   keller    : published + parental_edu_z_kernel:gender_c.
#
# Deliberately not included:
#   - contributing study as a covariate or in interactions. Study is near-collinear
#     with the design variables (SHIP is entirely East; BASE-II is almost entirely the
#     oldest birth years), so conditioning on it removes the contrasts of interest
#     rather than confounding. It was not prespecified and is out of scope.
#   - squared covariate terms (BYc_z^2 and its G/E interactions).
#   - ancestry principal components on the environment side.
#
# NOTE ON THE OUTCOME. mobility = edu_z_kernel - parental_edu_z_kernel, and
# parental_edu_z_kernel is on the right-hand side, so every coefficient here except the
# ParEdu main effect is numerically identical to the same model fitted to edu_z_kernel
# (educational attainment). The exports apply to BOTH outcomes; only the ParEdu main
# effect differs (by exactly +1).
#
# OUTPUTS. Five CSVs into the run's manuscript_export tree:
#   mobility/sensitivities/S7_ParEdu_Keller.csv        focal G x E terms, both specs
#   mobility/sensitivities/S7_ParEdu_SimpleSlopes.csv  PGI slope at -1/0/+1 SD ParEdu
#   mobility/sensitivities/S7_ParEdu_Lines.csv         fitted prediction grid
#   mobility/main/ParEdu_Gradient.csv                  regional social-origin gradient (S5)
#   mobility/sensitivities/ParEdu_Gradient_LOO.csv     its leave-one-study-out check (S6)
#
# The last two are NOT part of the exploratory Keller re-analysis and are NOT taken from
# B1. They are the social-origin gradient — how strongly educational attainment tracks
# parental education in each region — from a dedicated model fitted directly on
# edu_z_kernel: parental education x region, adjusted for birth year and gender, and
# NOT adjusted for PGI-Education. Conditioning on the PGI would strip the genetically
# mediated part out of parental transmission, which belongs in this quantity. See the
# block where m_grad is fitted. B1 is not modified; Findings 1 and 2 are untouched.
#
# DO NOT ADD A COHORT TREND HERE. The B1 right-hand side makes the gradient a function of
# birth year, and the trend looks compelling (p = 3e-25 linear; real curvature under the
# B3 ParEdu smooths). It is not reportable: it vanishes without TwinLife, and on a common
# birth-year window it is zero within SOEP alone and zero within TwinLife alone — it only
# appears when the two are combined, because the birth-year axis is stratified by study
# (BASE-II and TwinLife share no 5-year band).
#
# Fits models; nothing else.
# RUNS ON: Mac with the shared SMB volume mounted.

.find_this_file <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  # Rscript encodes spaces in --file= as "~+~", and this repository's path contains
  # one ("015_Data analysis"). Decode before normalizePath(), exactly as the other
  # section scripts do — without this the script only works when invoked with a
  # relative path, and fails from run_B.command.
  f <- gsub("~\\+~", " ", f, fixed = FALSE)
  if (length(f) && nzchar(f[1])) return(normalizePath(f[1], mustWork = TRUE))
  normalizePath(sys.frames()[[1]]$ofile)
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

# constants.R defines ggplot scales at source time, so ggplot2 must be attached first.
suppressPackageStartupMessages({ library(ggplot2); library(mgcv); library(dplyr) })

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL <- get_run_ts()
# Write into the run's manuscript_export tree, mirroring export_result_tables.R, so
# import_results.R can pick these up on the standard relative path
# (manuscript_export/mobility/sensitivities/...).
RUN_DIR <- dirname(output_dir("B_Mobility"))
OUT_DIR <- file.path(RUN_DIR, "manuscript_export", "mobility", "sensitivities")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# The regional parental-education gradient is a derived quantity of the B1 model, not a
# sensitivity, so it exports alongside the other main mobility tables (and lands in
# Supplementary Data 5). Its leave-one-study-out check stays on the sensitivities path.
OUT_DIR_MAIN <- file.path(RUN_DIR, "manuscript_export", "mobility", "main")
dir.create(OUT_DIR_MAIN, recursive = TRUE, showWarnings = FALSE)

source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()

DATA_ROOT <- data_root()
DATA_DIR  <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)

stopifnot(exists("dat_cluster_mob"), nrow(dat_cluster_mob) > 200)
d <- dat_cluster_mob
cat(sprintf("\nB_ParEdu_keller: N = %d\n", nrow(d)))

re_term <- "+ s(fid_re, bs = 're')"

paredu_b1 <- paste(
  "+ parental_edu_z_kernel",
  "+ parental_edu_z_kernel:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z",
  "+ parental_edu_z_kernel:BYc_z",
  "+ parental_edu_z_kernel:BYc_z:east_west_c",
  "+ parental_edu_z_kernel:PGI_Edu_z:BYc_z",
  "+ parental_edu_z_kernel:PGI_Edu_z:east_west_c")

cov3_block    <- "(PGI_Edu_z + BYc_z + east_west_c + gender_c)^3"
rhs_published <- paste(cov3_block, paredu_b1)
rhs_keller    <- paste(rhs_published, "+ parental_edu_z_kernel:gender_c")

fit <- function(rhs) mgcv::bam(as.formula(paste("mobility ~", rhs, re_term)),
                               data = d, method = "fREML", discrete = TRUE)

cat("fitting published spec...\n"); m_pub <- fit(rhs_published)
cat("fitting keller spec...\n");    m_kel <- fit(rhs_keller)

# ---- focal terms + region simple slopes, with the covariance term ----------------
T_MAIN <- "parental_edu_z_kernel:PGI_Edu_z"
T_XREG <- "PGI_Edu_z:parental_edu_z_kernel:east_west_c"

pick <- function(m, want) {
  nm  <- names(coef(m))
  hit <- nm[vapply(strsplit(nm, ":"), function(p)
    setequal(p, strsplit(want, ":")[[1]]), logical(1))]
  if (length(hit) != 1) stop("could not uniquely locate term: ", want)
  hit
}

summarise_spec <- function(m, label) {
  b  <- coef(m); V <- vcov(m)
  ia <- pick(m, T_MAIN); ib <- pick(m, T_XREG)
  # east_west_c is coded East +0.5 / West -0.5
  lc <- function(w) c(est = unname(b[ia] + w * b[ib]),
                      se  = unname(sqrt(V[ia, ia] + w^2 * V[ib, ib] + 2 * w * V[ia, ib])))
  rows <- rbind(
    c(term = "PGI x ParEdu (pooled)", lc(0)),
    c(term = "PGI x ParEdu x region", c(est = unname(b[ib]), se = unname(sqrt(V[ib, ib])))),
    c(term = "PGI x ParEdu | East",   lc(+0.5)),
    c(term = "PGI x ParEdu | West",   lc(-0.5)))
  out <- data.frame(spec = label, term = rows[, "term"],
                    estimate = as.numeric(rows[, "est"]),
                    se       = as.numeric(rows[, "se"]),
                    stringsAsFactors = FALSE)
  out$ci_lo <- out$estimate - 1.96 * out$se
  out$ci_hi <- out$estimate + 1.96 * out$se
  out$z     <- out$estimate / out$se
  out$p     <- 2 * pnorm(-abs(out$z))
  out$N     <- nrow(d)
  out
}

focal <- rbind(summarise_spec(m_pub, "published"), summarise_spec(m_kel, "keller"))
cat("\n-- focal terms and region simple slopes --\n")
print(focal, row.names = FALSE, digits = 4)

# ---- PGI slope at -1 / 0 / +1 SD parental education, by region -------------------
# Uses only terms involving PGI_Edu_z; birth year and gender held at 0 (their means).
simple_slopes <- function(m, label) {
  b <- coef(m); V <- vcov(m)
  i_pgi <- pick(m, "PGI_Edu_z"); i_reg <- pick(m, "PGI_Edu_z:east_west_c")
  i_pe  <- pick(m, T_MAIN);      i_x   <- pick(m, T_XREG)
  grid <- expand.grid(east_west = c("East", "West"), paredu_sd = c(-1, 0, 1),
                      stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    w  <- if (grid$east_west[k] == "East") 0.5 else -0.5
    pe <- grid$paredu_sd[k]
    cv <- setNames(numeric(length(b)), names(b))
    cv[i_pgi] <- 1; cv[i_reg] <- w; cv[i_pe] <- pe; cv[i_x] <- w * pe
    est <- sum(cv * b); se <- sqrt(drop(t(cv) %*% V %*% cv))
    data.frame(spec = label, east_west = grid$east_west[k], paredu_sd = pe,
               slope = est, se = se, ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se)
  }))
}

slopes <- rbind(simple_slopes(m_pub, "published"), simple_slopes(m_kel, "keller"))
cat("\n-- PGI slope across parental education, by region --\n")
print(slopes, row.names = FALSE, digits = 4)

# ---- fitted regression lines: education on PGI at low vs high parental education --
# Fitted on edu_z_kernel (educational attainment) rather than mobility so the lines are
# directly readable as "predicted education". Identical model, and identical
# coefficients apart from the ParEdu main effect (see header note); refitting simply
# gives the intercept structure the plot needs.
cat("\nfitting attainment-outcome version for the prediction grid...\n")
m_att <- mgcv::bam(as.formula(paste("edu_z_kernel ~", rhs_published, re_term)),
                   data = d, method = "fREML", discrete = TRUE)

# Predicted education at birth year = 0 and gender = 0 (both centred), so every term
# involving BYc_z or gender_c drops out. Built from the coefficient vector rather than
# predict() to avoid supplying random-effect levels for the family term.
pred_grid <- function(m) {
  b <- coef(m); V <- vcov(m)
  i0    <- "(Intercept)"
  i_pgi <- pick(m, "PGI_Edu_z")
  i_reg <- pick(m, "east_west_c")
  i_pe  <- pick(m, "parental_edu_z_kernel")
  i_pr  <- pick(m, "PGI_Edu_z:east_west_c")
  i_er  <- pick(m, "parental_edu_z_kernel:east_west_c")
  i_pe_pgi <- pick(m, T_MAIN)
  i_x      <- pick(m, T_XREG)

  # -1 SD / mean / +1 SD: the conventional simple-slopes probe of a linear moderator.
  # These are the same three anchors whose slopes S7_ParEdu_SimpleSlopes.csv reports.
  g <- expand.grid(pgi = seq(-2, 2, by = 0.25),
                   paredu_sd = c(-1, 0, 1),
                   east_west = c("East", "West"),
                   stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(g)), function(k) {
    w  <- if (g$east_west[k] == "East") 0.5 else -0.5
    pe <- g$paredu_sd[k]; x <- g$pgi[k]
    cv <- setNames(numeric(length(b)), names(b))
    cv[i0] <- 1; cv[i_pgi] <- x; cv[i_reg] <- w; cv[i_pe] <- pe
    cv[i_pr] <- x * w; cv[i_er] <- pe * w
    cv[i_pe_pgi] <- pe * x; cv[i_x] <- x * pe * w
    est <- sum(cv * b); se <- sqrt(drop(t(cv) %*% V %*% cv))
    data.frame(east_west = g$east_west[k], paredu_sd = pe, pgi = x,
               fit = est, se = se, ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se)
  }))
}

# Parental-education gradient by region, read directly off both fits. The region
# interaction is parameterisation-invariant; only the main effect shifts by +1 when the
# outcome changes from mobility to attainment. Printed so the claim in the supplement
# text is verified against the model rather than derived by hand.
# `shift` implements the mobility -> attainment reparameterisation exactly. Because
# mobility = edu_z_kernel - parental_edu_z_kernel and parental_edu_z_kernel is on the
# right-hand side, refitting on edu_z_kernel changes ONLY the ParEdu main effect, by
# exactly +1, and leaves the covariance matrix unchanged. Applying that algebraically to
# the mobility fit is preferable to a separate bam() refit on edu_z_kernel: bam
# re-estimates the family random-effect smoothing parameter per outcome, so an
# independent refit drifts by ~1e-4 on a coefficient that is invariant in theory. Taking
# the shift instead keeps these gradients exactly reconcilable with the published
# B1_Coefficients.csv, which is fitted on mobility.
paredu_gradient <- function(m, label, dat = d, shift = 0) {
  b <- coef(m); V <- vcov(m)
  ia <- pick(m, "parental_edu_z_kernel")
  ib <- pick(m, "parental_edu_z_kernel:east_west_c")
  b[ia] <- b[ia] + shift
  lc <- function(w) c(est = unname(b[ia] + w * b[ib]),
                      se  = unname(sqrt(V[ia, ia] + w^2 * V[ib, ib] + 2 * w * V[ia, ib])))
  rows <- rbind(c(term = "ParEdu (pooled)",     lc(0)),
                c(term = "ParEdu x region",     c(est = unname(b[ib]), se = unname(sqrt(V[ib, ib])))),
                c(term = "ParEdu | East",       lc(+0.5)),
                c(term = "ParEdu | West",       lc(-0.5)))
  out <- data.frame(outcome = label, term = rows[, "term"],
                    estimate = as.numeric(rows[, "est"]), se = as.numeric(rows[, "se"]),
                    stringsAsFactors = FALSE)
  out$ci_lo <- out$estimate - 1.96 * out$se
  out$ci_hi <- out$estimate + 1.96 * out$se
  out$z <- out$estimate / out$se; out$p <- 2 * pnorm(-abs(out$z))
  out$N <- nrow(dat)
  out
}

# ---- the reported social-origin gradient --------------------------------------
# HOW STRONGLY DOES PARENTAL EDUCATION PREDICT EDUCATIONAL ATTAINMENT, BY REGION.
#
# Fitted directly on edu_z_kernel. Not from B1, not on the mobility metric, and
# NOT ADJUSTED FOR THE PGI.
#
# Why no PGI. Parental education and the offspring PGI are correlated through genetic
# transmission, so conditioning on the PGI strips the genetically mediated part out of
# the social-origin gradient. The reported quantity is meant to describe how closely
# attainment tracks social origin, so that part belongs in it. Adjusting for the PGI
# removed about 9% of the pooled association and about 15% of the East gradient:
#   PGI-adjusted   pooled 0.456  East 0.341  West 0.571  E-W -0.229
#   no PGI         pooled 0.501  East 0.393  West 0.608  E-W -0.216
# The contrast is barely affected; the levels are.
#
# Birth year and gender are kept as covariates for consistency with the rest of the
# paper, but they change almost nothing here — edu_z_kernel and parental_edu_z_kernel
# are both already kernel-standardised WITHIN birth year, so that adjustment is baked
# in. The bare ParEdu * region model gives pooled 0.504 / East 0.393 / West 0.615.
#
# B1 IS NOT MODIFIED. This is a separate fit; Findings 1 and 2 are untouched.
rhs_grad <- "parental_edu_z_kernel * east_west_c + BYc_z + gender_c"
cat("fitting the social-origin gradient model...\n")
m_grad <- mgcv::bam(as.formula(paste("edu_z_kernel ~", rhs_grad, re_term)),
                    data = d, method = "fREML", discrete = TRUE)

.stamp_spec <- function(df) {
  df$spec <- paste("Educational attainment on parental education and its region",
                   "interaction, adjusted for birth year and gender. NOT adjusted for",
                   "PGI-Education, and not conditioned on any covariate value.")
  df
}

gradient_df <- .stamp_spec(paredu_gradient(m_grad, "attainment", dat = d))
cat("\n-- social-origin gradient in educational attainment --\n")
print(gradient_df[, c("outcome", "term", "estimate", "se", "ci_lo", "ci_hi", "p")],
      row.names = FALSE, digits = 4)

lines_df <- pred_grid(m_att)
cat("\n-- fitted slopes implied by the prediction grid (per 1 SD PGI) --\n")
print(do.call(rbind, lapply(split(lines_df, list(lines_df$east_west, lines_df$paredu_sd)),
  function(s) data.frame(east_west = s$east_west[1], paredu_sd = s$paredu_sd[1],
                         slope = coef(lm(fit ~ pgi, data = s))[["pgi"]]))),
  row.names = FALSE, digits = 4)

# ---- leave-one-study-out on the parental-education x region contrast -------------
# One refit per fold on mobility, shifted to the attainment metric, yields both the
# contrast (which is parameterisation-invariant) and the region-specific gradients.
# Folds come from cohort_groups(), which builds them from the
# cohorts actually present — SHIP has no parental education and is absent from the
# mobility sample, so the folds are -BASE-II, -SOEP, -TwinLife.
#
# This is the study-sensitivity check for the gradient the manuscript leans on: if the
# East-West difference only survives with a particular study in the sample, that has to
# be visible before any framing rests on it.
cat("\nLOCO: leave-one-study-out on the ParEdu x region contrast...\n")

loco_fit_paredu <- function(sub) {
  sub$fid_re <- droplevels(factor(sub$fid_re))
  rt <- if (nlevels(sub$fid_re) >= 2) re_term else ""
  # Same specification and outcome as the reported gradient, refit per fold.
  m  <- mgcv::bam(as.formula(paste("edu_z_kernel ~", rhs_grad, rt)),
                  data = sub, method = "fREML", discrete = TRUE)
  b <- coef(m); V <- vcov(m)
  ia <- pick(m, "parental_edu_z_kernel")
  ib <- pick(m, "parental_edu_z_kernel:east_west_c")
  lc <- function(w) c(unname(b[ia] + w * b[ib]),
                      unname(sqrt(V[ia, ia] + w^2 * V[ib, ib] + 2 * w * V[ia, ib])))
  ew <- c(unname(b[ib]), unname(sqrt(V[ib, ib])))
  po <- lc(0); ea <- lc(+0.5); we <- lc(-0.5)
  c(estimate = round(ew[1], 5), SE = round(ew[2], 5),
    ci_lo = round(ew[1] - 1.96 * ew[2], 5), ci_hi = round(ew[1] + 1.96 * ew[2], 5),
    p = signif(2 * pnorm(-abs(ew[1] / ew[2])), 4),
    gradient_pooled = round(po[1], 5), gradient_pooled_SE = round(po[2], 5),
    gradient_East   = round(ea[1], 5), gradient_East_SE   = round(ea[2], 5),
    gradient_West   = round(we[1], 5), gradient_West_SE   = round(we[2], 5))
}

paredu_loo <- build_loo_table(
  leave_one_cohort_out(d, loco_fit_paredu, folds = cohort_groups(d)),
  estimate_col = "estimate")
paredu_loo$outcome <- "attainment"
paredu_loo$focal_term <- "parental_edu_z_kernel:east_west_c"
paredu_loo <- .stamp_spec(paredu_loo)

cat("\n-- LOCO: ParEdu x region contrast, attainment metric --\n")
print(paredu_loo[, c("fold", "N", "estimate", "SE", "ci_lo", "ci_hi", "p",
                     "gradient_East", "gradient_West", "delta_vs_full", "status")],
      row.names = FALSE, digits = 4)

f1 <- file.path(OUT_DIR, "S7_ParEdu_Keller.csv")
f2 <- file.path(OUT_DIR, "S7_ParEdu_SimpleSlopes.csv")
f3 <- file.path(OUT_DIR, "S7_ParEdu_Lines.csv")
f4 <- file.path(OUT_DIR_MAIN, "ParEdu_Gradient.csv")
f5 <- file.path(OUT_DIR, "ParEdu_Gradient_LOO.csv")
write.csv(focal,       f1, row.names = FALSE)
write.csv(slopes,      f2, row.names = FALSE)
write.csv(lines_df,    f3, row.names = FALSE)
write.csv(gradient_df, f4, row.names = FALSE)
write.csv(paredu_loo,  f5, row.names = FALSE)
cat("\nwrote:\n  ", f1, "\n  ", f2, "\n  ", f3, "\n  ", f4, "\n  ", f5, "\n")
