# A_sensitivities_analysis.R
# Section A (Educational attainment) sensitivities — analysis stage.
# All sensitivities use the dedup-OLS sample (dat_edu) except the
# finding-matched S1-A, CR2/RE (full-family dat_cluster) and S3-A (gaulss on
# dat_cluster without the family RE). See ../../../Plan_deviations.md §11/§12.
#
#   S1-A   Alternative-PGI family (PGI-Edu + PGI-Cog + PGI-NonCog), BH per
#          family, on the A1 focal and the A4 four-way gender focal; plus the
#          finding-matched family-RE version on a common re-residualised sample.
#   S2-A   Height negative control (PGI-Height) on edu_z_kernel (A1 spec) and
#          on height_z (sanity).
#   S4     Migration diagnostics (A-only): S4a/b/c + gated S4d. Outer gate
#          on A focal significance (read from A_Education_Models.rds).
#   Method A1 focal under dedup-OLS / CR2 / RE (RE read from the A RDS).
#   S3-A   Heteroscedasticity (Gaussian location-scale gaulss), per region,
#          gated on A focal significance.
#
# Usage (from run_A.command, after A_analysis.R): Rscript A_Education/A_sensitivities_analysis.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
})

# ---- Resolve paths ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine A_sensitivities_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
RUN_TS_VAL <- get_run_ts()
OUT_DIR_A  <- output_dir("A_Education")

source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
RESULTS <- list(); PLOTS <- list(); PRIMARY_RESULTS <- list()

DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var
DATA_DIR        <- file.path(DATA_ROOT, "03_Merge")
EXCLUDE_TWINLIFE <- identical(Sys.getenv("EXCLUDE_TWINLIFE"), "1")
source(file.path(WORKFLOW_DIR, "00_setup", "load_data.R"),       local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "prepare_samples.R"), local = TRUE)
cat(sprintf("A_sensitivities: dat_edu N=%d, dat_cluster N=%d\n",
            nrow(dat_edu), if (exists("dat_cluster")) nrow(dat_cluster) else 0L))

A1_focal <- "PGI_Edu_z:BYc_z:east_west_c"
A4_focal <- "PGI_Edu_z:BYc_z:east_west_c:gender_c"

# Alternative PGIs available on the analytic sample.
alt_pgis <- c("PGI_Cog_z", "PGI_nonCog_z")
alt_pgis <- alt_pgis[vapply(alt_pgis, function(v)
  v %in% names(dat_edu) && !all(is.na(dat_edu[[v]])), logical(1))]
s1_pgis <- c("PGI_Edu_z", alt_pgis)

# =====================================================================
# S1-A: alternative-PGI family (BH per family) — A1 focal + A4 four-way
# =====================================================================
cat("S1-A: alternative-PGI families...\n")
.s1_family <- function(dat, formula_fun, focal_fun, family_label) {
  rows <- list()
  for (pgi in s1_pgis) {
    lbl   <- sub("_z$", "", pgi); focal <- focal_fun(pgi)
    m <- lm(as.formula(formula_fun(pgi)), data = dat)
    cs <- coef(summary(m))
    if (focal %in% rownames(cs))
      rows[[lbl]] <- data.frame(family = family_label, PGI = lbl, focal_term = focal,
        estimate = cs[focal, "Estimate"], se = cs[focal, "Std. Error"],
        t = cs[focal, "t value"], p_raw = cs[focal, 4], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (!is.null(out)) out$p_BH <- p.adjust(out$p_raw, method = "BH")
  out
}
S1_A1 <- .s1_family(dat_edu,
  function(pgi) paste("edu_z_kernel ~ (", pgi, "+ BYc_z + east_west_c + gender_c)^3"),
  function(pgi) paste0(pgi, ":BYc_z:east_west_c"), "A1 (PGI x BY x region)")
S1_A4 <- .s1_family(dat_edu,
  function(pgi) paste0("edu_z_kernel ~ ", pgi, " * BYc_z * east_west_c * gender_c"),
  function(pgi) paste0(pgi, ":BYc_z:east_west_c:gender_c"), "A4 (4-way PGI x BY x region x gender)")
S1_AltPGIs <- dplyr::bind_rows(S1_A1, S1_A4)

# =====================================================================
# S2-A: height negative control
# =====================================================================
cat("S2-A: height negative control...\n")
# Height negative control on the A1 specification, plus a height -> height
# sanity check.
S2_Height_Education <- S2_Height_Height <- NULL
if ("PGI_Height_z" %in% names(dat_edu) && !all(is.na(dat_edu$PGI_Height_z))) {
  m_s2 <- lm(edu_z_kernel ~ (PGI_Height_z + BYc_z + east_west_c + gender_c)^3, data = dat_edu)
  S2_Height_Education <- tidy_coefs(m_s2, "S2_Height_Edu (A1 spec)")
  if ("height_z" %in% names(dat_edu) && !all(is.na(dat_edu$height_z))) {
    m_s2h <- lm(height_z ~ (PGI_Height_z + BYc_z + east_west_c + gender_c)^3,
                data = dat_edu |> dplyr::filter(!is.na(height_z)))
    S2_Height_Height <- tidy_coefs(m_s2h, "S2_Height_Height (sanity)")
  }
}

# =====================================================================
# Read A_Education_Models.rds (S4 outer gate + S3-A gate + Method RE row)
# =====================================================================
# GATE_RUN_TS re-points the gate cache at another run, for the case where only
# this stage is re-run into a fresh run folder: the gates then read the run the
# gated decisions were made in, instead of degrading to "unavailable".
A_GATE_TS  <- Sys.getenv("GATE_RUN_TS", RUN_TS_VAL)
A_GATE_DIR <- file.path(dirname(dirname(OUT_DIR_A)), paste0("RUN_", A_GATE_TS), "A_Education")
A_RDS <- file.path(A_GATE_DIR, "A_Education_Models.rds")
A_mc  <- if (file.exists(A_RDS)) tryCatch(readRDS(A_RDS), error = function(e) NULL) else NULL
gxe_a_sig <- NA; a_gate_basis <- "unavailable"
if (!is.null(A_mc)) {
  a1_p <- { r <- A_mc$A1_Coefficients[A_mc$A1_Coefficients$term == A1_focal, ]
            if (nrow(r)) r$p.value[1] else NA_real_ }
  a3_p <- { sm <- A_mc$A3_Smooths_tp
            if (!is.null(sm) && nrow(sm)) suppressWarnings(min(sm[["p-value"]], na.rm = TRUE)) else NA_real_ }
  gxe_a_sig <- isTRUE(a1_p < 0.05) || isTRUE(a3_p < 0.05)
  a_gate_basis <- sprintf("RE A1 focal p=%s; min A3 smooth p=%s", signif(a1_p, 3), signif(a3_p, 3))
} else {
  m_a1_ols <- lm(edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3, data = dat_edu)
  a1_p_ols <- coef(summary(m_a1_ols))[A1_focal, 4]
  gxe_a_sig <- isTRUE(a1_p_ols < 0.05)
  a_gate_basis <- sprintf("A_RDS absent; fallback dedup-OLS A1 focal p=%s (A3 not checked)", signif(a1_p_ols, 3))
}
cat(sprintf("A gate: GxE_A_sig = %s [%s]\n", gxe_a_sig, a_gate_basis))

# =====================================================================
# S1-A (finding-matched, random effects) — ../../../Plan_deviations.md §15
# Alternative PGIs (Education / Cognition / Noncognitive) on the reported focal
# attainment findings, fit with the primary family-RE specification (bam,
# fREML, s(fid_re)) on a COMMON complete-case sample across the three PGIs;
# each PGI is re-residualised on crossPC1-10 and re-z-scored WITHIN that sample
# so the comparison is not confounded by sample or standardisation differences,
# and PGI-Education is fit on the same common sample. Covers, alongside the
# dedup-OLS S1 above (A1 three-way + A4 four-way):
#   Finding 1  pooled PGI x BY (+ contextual PGI x BY x region), A1-equivalent.
#   Finding 3  pooled PGI x gender + gender-specific slopes, A4-equivalent.
#   A5 per-gender GAM, reported for symmetry with B5 (a stable null).
# BH is applied within each focal hypothesis across the three PGIs. Model
# specifications are derived from the cached primary fits (A_mc) so the
# alternative-PGI models are structurally identical bar the PGI variable.
# =====================================================================
cat("S1-A (finding-matched, RE): common-sample alternative PGIs...\n")
S1fm_raw <- c("PGI_Edu", "PGI_Cog", "PGI_nonCog")
S1fm_raw <- S1fm_raw[vapply(S1fm_raw, function(v)
  exists("dat_cluster") && v %in% names(dat_cluster) &&
    !all(is.na(dat_cluster[[v]])), logical(1))]
S1fm_ok <- !is.null(A_mc) && exists("dat_cluster") &&
  all(paste0("crossPC", 1:10) %in% names(dat_cluster)) && length(S1fm_raw) >= 2
S1_finding1 <- S1_finding3_att <- S1_finding3_att_slopes <- S1_A5_nonlin <- NULL
S1fm_att_N <- NA_integer_; S1fm_status <- "not_run"
if (S1fm_ok) {
  S1fm_att   <- altpgi_common_sample(dat_cluster, pgis_raw = S1fm_raw)
  S1fm_pgis  <- paste0(S1fm_raw, "_z")
  S1fm_att_N <- nrow(S1fm_att)
  a1_tmpl <- gsub("PGI_Edu_z", "PGI", rhs_without_re(A_mc$m_A1_re), fixed = TRUE)
  a4_tmpl <- gsub("PGI_Edu_z", "PGI", rhs_without_re(A_mc$m_A4_re), fixed = TRUE)
  S1_finding1 <- altpgi_re_family(S1fm_att, a1_tmpl,
    list(`PGI x BY (pooled)` = "PGI:BYc_z",
         `PGI x BY x region`  = "PGI:BYc_z:east_west_c"),
    "S1_finding1 (attainment across cohorts)", pgis = S1fm_pgis)
  S1_finding3_att <- altpgi_re_family(S1fm_att, a4_tmpl,
    list(`PGI x gender` = "PGI:gender_c"),
    "S1_finding3 (attainment gender)", pgis = S1fm_pgis)
  S1_finding3_att_slopes <- altpgi_slopes_family(S1fm_att, a4_tmpl,
    list(Female = list(gender_c = 0.5), Male = list(gender_c = -0.5)),
    "S1_finding3 (attainment gender) slopes", pgis = S1fm_pgis)
  if (!is.null(A_mc$m_A5_lin_re) && !is.null(A_mc$m_A5_nonlin_re)) {
    .a5fm <- altpgi_pergender_gam(S1fm_att,
      rhs_without_re(A_mc$m_A5_lin_re), rhs_without_re(A_mc$m_A5_nonlin_re),
      "S1_A5_nonlin (attainment per-gender GAM; null check)",
      pgis = S1fm_pgis, want_diffsmooth = FALSE)
    S1_A5_nonlin <- .a5fm$tests
  }
  S1fm_status <- sprintf("run on common N=%d; PGIs=%s", S1fm_att_N,
                         paste(sub("_z$", "", S1fm_pgis), collapse = "/"))
} else {
  S1fm_status <- "skipped (A_mc, raw PGIs or crossPCs unavailable)"
}
S1_FindingMatched_Status <- data.frame(
  section = "A", status = S1fm_status, common_N = S1fm_att_N,
  bh_families = "BH within each focal hypothesis across PGI-Education/Cognition/Noncognitive",
  stringsAsFactors = FALSE)
cat(sprintf("  S1-A finding-matched: %s\n", S1fm_status))

# =====================================================================
# S4: migration (A-only; outer gate on A focal significance)
# =====================================================================
S4_Status <- data.frame(gate = "A focal significance", GxE_A_significant = gxe_a_sig,
  basis = a_gate_basis, triggered = gxe_a_sig,
  note = if (isTRUE(gxe_a_sig)) "S4 migration diagnostics run." else
    "S4 not triggered (no A focal significant at alpha=.05).", stringsAsFactors = FALSE)
S4_Migration_Tests <- S4d_GatingDecision <- S4d_Compare <- S4d_Coefficients <- NULL
if (isTRUE(gxe_a_sig)) {
  cat("S4: migration diagnostics (triggered)...\n")
  pe <- dat_edu$PGI_Edu[dat_edu$east == 1]; pw <- dat_edu$PGI_Edu[dat_edu$east == 0]
  ee <- dat_edu$education[dat_edu$east == 1]; ew <- dat_edu$education[dat_edu$east == 0]
  tt <- t.test(pe, pw, alternative = "less"); fp <- var.test(pe, pw, alternative = "less")
  # Test 3 is an EQUIVALENCE claim ("education variance does not differ to a
  # practically important degree between regions"), so it is tested with two
  # one-sided F-tests (TOST) against a variance-ratio bound - NOT by failing to
  # reject a difference (which an under-powered test passes trivially).
  # Equivalence is declared only if BOTH one-sided tests reject, i.e. the ratio
  # var(Edu_East)/var(Edu_West) sits inside [1/bound, bound].
  EDU_VAR_EQUIV_BOUND <- 1.5
  fe      <- var.test(ee, ew)                                             # descriptive ratio
  tost_lo <- var.test(ee, ew, ratio = 1 / EDU_VAR_EQUIV_BOUND, alternative = "greater")
  tost_hi <- var.test(ee, ew, ratio = EDU_VAR_EQUIV_BOUND,     alternative = "less")
  tost_p  <- max(tost_lo$p.value, tost_hi$p.value)                        # binding one-sided p
  S4_Migration_Tests <- data.frame(
    test = c("S4a: mean(PGI) East<West (Welch t)", "S4b: var(PGI) East<West (F)",
             sprintf("S4c: var(Edu) East~=West (TOST equivalence, ratio bound %.2f)", EDU_VAR_EQUIV_BOUND)),
    statistic = c(tt$statistic, fp$statistic, fe$statistic),
    p_value = c(tt$p.value, fp$p.value, tost_p),
    East_stat = c(mean(pe, na.rm = TRUE), var(pe, na.rm = TRUE), var(ee, na.rm = TRUE)),
    West_stat = c(mean(pw, na.rm = TRUE), var(pw, na.rm = TRUE), var(ew, na.rm = TRUE)),
    stringsAsFactors = FALSE)
  t1 <- isTRUE(tt$p.value < 0.05); t2 <- isTRUE(fp$p.value < 0.05)
  t3 <- isTRUE(tost_lo$p.value < 0.05 && tost_hi$p.value < 0.05)          # equivalence established
  s4d_run <- t1 && t2 && t3
  S4d_GatingDecision <- data.frame(test1_East_mean_lt_West_p_lt_05 = t1,
    test2_East_var_lt_West_p_lt_05 = t2, test3_var_edu_equivalent_TOST = t3,
    run_S4d = s4d_run, stringsAsFactors = FALSE)
  if (s4d_run) {
    cat("S4d: inner gate passes; within-region PGI refit.\n")
    # Standardize PGI within region x birth cohort (local standing), not region
    # only: the rank-restriction mechanism operates within region x birth
    # cohort. Region-specific kernel-z of PGI over birth year - the same
    # within-birth-year metric as edu_z_kernel. See ../../../ANALYSIS_PLAN.md.
    dat_edu$PGI_Edu_z_within <- NA_real_
    for (rg in c("West", "East")) {
      ix <- which(dat_edu$east_west == rg)
      dat_edu$PGI_Edu_z_within[ix] <- kernel_z_score(
        dat_edu$PGI_Edu[ix], dat_edu$birth_year[ix], bandwidth = KERNEL_BW)
    }
    # Compare on the common non-NA subset so the attenuation is apples-to-apples
    # (kernel-z is undefined at the birth-year boundaries).
    sub_s4d <- dat_edu[!is.na(dat_edu$PGI_Edu_z_within), , drop = FALSE]
    m_s4d <- lm(edu_z_kernel ~ (PGI_Edu_z_within + BYc_z + east_west_c + gender_c)^3, data = sub_s4d)
    m_a1  <- lm(edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3, data = sub_s4d)
    f_loc <- "PGI_Edu_z_within:BYc_z:east_west_c"
    cs_loc <- coef(summary(m_s4d)); cs_abs <- coef(summary(m_a1))
    S4d_Compare <- dplyr::bind_rows(
      data.frame(scope = "absolute PGI (A1, dedup-OLS)", term = A1_focal,
        estimate = cs_abs[A1_focal, 1], se = cs_abs[A1_focal, 2], p = cs_abs[A1_focal, 4], stringsAsFactors = FALSE),
      data.frame(scope = "within-region PGI (S4d)", term = f_loc,
        estimate = cs_loc[f_loc, 1], se = cs_loc[f_loc, 2], p = cs_loc[f_loc, 4], stringsAsFactors = FALSE))
    S4d_Compare$attenuation_pct <- c(NA_real_,
      round(100 * (S4d_Compare$estimate[1] - S4d_Compare$estimate[2]) / S4d_Compare$estimate[1], 1))
    S4d_Coefficients <- tidy_coefs(m_s4d, "S4d_WithinRegionPGI")
  } else {
    S4d_Compare <- data.frame(note = sprintf(
      "S4d skipped: inner gate failed (S4a p=%.3g, S4b p=%.3g, S4c TOST p=%.3g).",
      tt$p.value, fp$p.value, tost_p), stringsAsFactors = FALSE)
  }
}

# =====================================================================
# Method comparison: A1 focal under dedup-OLS / CR2 / RE
# =====================================================================
cat("Method comparison: dedup-OLS / CR2 / RE for the A1 focal...\n")
.cr2_available <- requireNamespace("clubSandwich", quietly = TRUE)
m_dedup <- lm(edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3, data = dat_edu)
cs_d <- coef(summary(m_dedup))
row_dedup <- data.frame(estimator = "dedup-OLS (one-per-family)", sample = "dat_edu",
  N = nrow(dat_edu), estimate = cs_d[A1_focal, 1], se = cs_d[A1_focal, 2],
  p = cs_d[A1_focal, 4], df = m_dedup$df.residual, stringsAsFactors = FALSE)
row_cr2 <- NULL
if (exists("dat_cluster") && "fid_re" %in% names(dat_cluster) && .cr2_available) {
  m_full <- lm(edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3, data = dat_cluster)
  ct <- tryCatch(clubSandwich::coef_test(m_full, vcov = "CR2", cluster = dat_cluster$fid_re,
                                         test = "Satterthwaite"), error = function(e) NULL)
  if (!is.null(ct) && A1_focal %in% rownames(ct)) {
    r <- ct[A1_focal, ]
    row_cr2 <- data.frame(estimator = "CR2 cluster-robust (by family)", sample = "dat_cluster",
      N = nrow(dat_cluster), estimate = r[["beta"]], se = r[["SE"]],
      p = r[["p_Satt"]], df = round(r[["df_Satt"]], 1), stringsAsFactors = FALSE)
  }
}
row_re <- NULL
if (!is.null(A_mc) && !is.null(A_mc$A1_Coefficients)) {
  r <- A_mc$A1_Coefficients[A_mc$A1_Coefficients$term == A1_focal, ]
  if (nrow(r)) row_re <- data.frame(estimator = "RE s(fid_re) (read from A RDS; not refit)",
    sample = "dat_cluster", N = nrow(A_mc$dat_cluster), estimate = r$estimate[1],
    se = r$std.error[1], p = r$p.value[1], df = NA_real_, stringsAsFactors = FALSE)
}
Method_Comparison <- dplyr::bind_rows(row_dedup, row_cr2, row_re)
Method_Status <- data.frame(focal_term = A1_focal, cr2_available = .cr2_available,
  a_rds_found = !is.null(A_mc),
  note = paste0(if (!.cr2_available) "clubSandwich not installed -> CR2 row omitted. " else "",
    if (is.null(A_mc)) "A RDS absent -> RE row omitted. " else "RE row read from A RDS (A1 not refit). "),
  stringsAsFactors = FALSE)

# =====================================================================
# S3-A: heteroscedasticity (gaulss), per region, gated on A focal sig.
# gaulss cannot carry the family RE, so dat_cluster is used WITHOUT it
# (documented caveat; ../../../Plan_deviations.md §12).
# =====================================================================
cat("S3-A: heteroscedasticity (gaulss)...\n")
.sm <- function(by) sprintf("s(BYc, by = %s, k = %d, bs = '%s')", by, K_DEFAULT, BS_DEFAULT)
.disp_table <- function(mod, region, kind) {
  if (is.null(mod)) return(NULL)
  cs <- summary(mod)$p.table; if (is.null(cs) || !nrow(cs)) return(NULL)
  sr <- grep("\\.1$", rownames(cs)); if (!length(sr)) return(NULL)
  zc <- if ("z value" %in% colnames(cs)) "z value" else "t value"
  pc <- if ("Pr(>|z|)" %in% colnames(cs)) "Pr(>|z|)" else "Pr(>|t|)"
  data.frame(region = region, model = kind, term = sub("\\.1$", "", rownames(cs)[sr]),
             estimate = cs[sr, "Estimate"], std.error = cs[sr, "Std. Error"],
             z = cs[sr, zc], p_value = cs[sr, pc], stringsAsFactors = FALSE, row.names = NULL)
}
.focal_smooth_row <- function(mod, region, variance_kind) {
  if (is.null(mod)) return(NULL)
  st <- summary(mod)$s.table; if (is.null(st) || !nrow(st)) return(NULL)
  rn <- rownames(st); idx <- which(grepl("PGI_Edu_z", rn) & !grepl("parental_edu", rn))
  if (!length(idx)) return(NULL)
  statcol <- intersect(c("F", "Chi.sq"), colnames(st))[1]
  data.frame(region = region, variance = variance_kind, smooth = rn[idx[1]],
             edf = st[idx[1], "edf"], statistic = st[idx[1], statcol],
             p_value = st[idx[1], "p-value"], stringsAsFactors = FALSE, row.names = NULL)
}
.fit_region <- function(d, region, outcome, mean_lin, mean_va, disp_env, disp_gen, disp_va) {
  if (nrow(d) < 200) return(list(note = sprintf("Skipped (%s N=%d < 200)", region, nrow(d))))
  gl <- function(mf, df, kind) tryCatch(
    mgcv::gam(list(as.formula(paste(outcome, "~", mf)), df), family = mgcv::gaulss(), data = d, method = "REML"),
    error = function(e) { cat(sprintf("  %s (%s) failed: %s\n", kind, region, conditionMessage(e))); NULL })
  va_homo <- tryCatch(mgcv::gam(as.formula(paste(outcome, "~", mean_va)), data = d, method = "REML"),
                      error = function(e) NULL)
  list(env = gl(mean_lin, disp_env, "env"), gen = gl(mean_lin, disp_gen, "gen"),
       va_het = gl(mean_va, disp_va, "var_aware"), va_homo = va_homo)
}
S3A <- NULL
if (isTRUE(gxe_a_sig)) {
  mean_lin <- "BYc + gender_c + PGI_Edu_z + PGI_Edu_z:BYc"
  mean_va  <- paste("BYc + gender_c + PGI_Edu_z +", .sm("PGI_Edu_z"))
  fits <- lapply(c("East", "West"), function(reg)
    .fit_region(dat_cluster[as.character(dat_cluster$east_west) == reg, , drop = FALSE], reg,
                "edu_z_kernel", mean_lin, mean_va,
                ~ BYc, ~ PGI_Edu_z, ~ BYc + PGI_Edu_z + BYc:PGI_Edu_z))
  names(fits) <- c("East", "West")
  disp <- list(); foc <- list()
  for (reg in names(fits)) {
    f <- fits[[reg]]; if (!is.null(f$note)) next
    disp[[paste(reg, "env")]] <- .disp_table(f$env, reg, "env (sigma~BYc)")
    disp[[paste(reg, "gen")]] <- .disp_table(f$gen, reg, "gen (sigma~PGI)")
    disp[[paste(reg, "va")]]  <- .disp_table(f$va_het, reg, "var_aware (sigma~BYc+PGI+BYc:PGI)")
    foc[[paste(reg, "homo")]] <- .focal_smooth_row(f$va_homo, reg, "homoscedastic")
    foc[[paste(reg, "het")]]  <- .focal_smooth_row(f$va_het,  reg, "heteroscedastic (gaulss)")
  }
  disp <- disp[!vapply(disp, is.null, logical(1))]; foc <- foc[!vapply(foc, is.null, logical(1))]
  notes <- vapply(fits, function(f) if (!is.null(f$note)) f$note else "OK", character(1))
  S3A <- list(
    dispersion = if (length(disp)) do.call(rbind, disp) else data.frame(note = "No dispersion summaries."),
    focal_under_het = if (length(foc)) do.call(rbind, foc) else data.frame(note = "No focal smooth extracted."),
    region_status = data.frame(region = names(notes), status = unname(notes), stringsAsFactors = FALSE))
}
S3_Status <- data.frame(section = "A-S3 (education)", GxE_significant = gxe_a_sig,
  basis = a_gate_basis, triggered = isTRUE(gxe_a_sig),
  note = if (isTRUE(gxe_a_sig)) "S3-A heteroscedasticity GAMs fitted." else "Not triggered; S3-A not fitted.",
  stringsAsFactors = FALSE)

# =====================================================================
# Focal dispersion robustness (per-focal gated): for each Section A focal
# that is significant in the main RE analysis, re-estimate the mean focal
# allowing the residual variance to follow the same structure (gaulss, no RE).
# =====================================================================
cat("Focal dispersion robustness (gated on main-analysis significance)...\n")
.cp <- function(tbl, term) { if (is.null(tbl)) return(NA_real_)
  r <- tbl[tbl$term == term, ]; if (nrow(r)) r$p.value[1] else NA_real_ }
Focal_Dispersion <- if (!is.null(A_mc) && exists("dat_cluster")) focal_dispersion_check(
  dat_cluster, "edu_z_kernel", list(
    list(label = "PGI x birth year x region (A1)",
         mean_rhs = "PGI_Edu_z*BYc_z*east_west_c + gender_c", focal = "PGI_Edu_z:BYc_z:east_west_c",
         scale_rhs = "PGI_Edu_z*BYc_z*east_west_c", main_p = .cp(A_mc$A1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c")),
    list(label = "PGI x gender (pooled across region)",
         mean_rhs = "PGI_Edu_z*gender_c + east_west_c + BYc_z", focal = "PGI_Edu_z:gender_c",
         scale_rhs = "PGI_Edu_z*gender_c", main_p = .cp(A_mc$A4_Coefficients, "PGI_Edu_z:gender_c")),
    list(label = "PGI x region x gender",
         mean_rhs = "PGI_Edu_z*east_west_c*gender_c + BYc_z", focal = "PGI_Edu_z:east_west_c:gender_c",
         scale_rhs = "PGI_Edu_z*east_west_c*gender_c", main_p = .cp(A_mc$A4_Coefficients, "PGI_Edu_z:east_west_c:gender_c")),
    list(label = "PGI x birth year x region x gender (A4)",
         mean_rhs = "PGI_Edu_z*BYc_z*east_west_c*gender_c", focal = "PGI_Edu_z:BYc_z:east_west_c:gender_c",
         scale_rhs = "PGI_Edu_z*east_west_c*gender_c", main_p = .cp(A_mc$A4_Coefficients, "PGI_Edu_z:BYc_z:east_west_c:gender_c")))) else NULL

# =====================================================================
# Save
# =====================================================================
model_cache <- list(
  RUN_TS = RUN_TS_VAL, s1_pgis = s1_pgis,
  S1_AltPGIs = S1_AltPGIs,
  S1_FindingMatched_Status = S1_FindingMatched_Status,
  S1_finding1 = S1_finding1,
  S1_finding3_att = S1_finding3_att, S1_finding3_att_slopes = S1_finding3_att_slopes,
  S1_A5_nonlin = S1_A5_nonlin,
  S2_Height_Education = S2_Height_Education, S2_Height_Height = S2_Height_Height,
  S4_Status = S4_Status, S4_Migration_Tests = S4_Migration_Tests,
  S4d_GatingDecision = S4d_GatingDecision, S4d_Compare = S4d_Compare, S4d_Coefficients = S4d_Coefficients,
  Method_Comparison = Method_Comparison, Method_Status = Method_Status,
  S3_Status = S3_Status, A_S3 = S3A, Focal_Dispersion = Focal_Dispersion,
  N_edu = nrow(dat_edu), N_cluster = if (exists("dat_cluster")) nrow(dat_cluster) else 0L)
OUT_RDS <- file.path(OUT_DIR_A, "A_Sensitivities_Models.rds")
saveRDS(model_cache, OUT_RDS)
cat(sprintf("A_sensitivities_analysis: models RDS saved to %s\n", OUT_RDS))
