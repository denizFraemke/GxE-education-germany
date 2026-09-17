# B_sensitivities_analysis.R
# Section B (Educational mobility) sensitivities — analysis stage.
# Mirrors A_Education/A_sensitivities_analysis.R for the mobility phenotype.
# All sensitivities carry the parental-education control (mobility = own minus
# parental education; ../../Plan_deviations.md §9/§11/§12).
#
#   S1-B   Alternative-PGI family (PGI-Edu + PGI-Cog + PGI-NonCog), BH per
#          family, on the B1 focal and the B4 four-way gender focal (dat_mob).
#   S2-B   Height negative control (PGI-Height) on mobility (B1 spec), plus a
#          PGI-Height -> height sanity check.
#   S3-B   Heteroscedasticity (Gaussian location-scale gaulss), per region, on
#          mobility, gated on B focal significance; also the B0 East-West
#          PGI-mobility contrast under heteroscedasticity (OLS vs gaulss).
#
# S4 migration and the dedup/CR2/RE method comparison are Section A only.
#
# Usage (from run_B.command, after B_analysis.R): Rscript B_Mobility/B_sensitivities_analysis.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
})

.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine B_sensitivities_analysis.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
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
cat(sprintf("B_sensitivities: dat_mob N=%d, dat_cluster_mob N=%d\n",
            nrow(dat_mob), if (exists("dat_cluster_mob")) nrow(dat_cluster_mob) else 0L))

B1_focal <- "PGI_Edu_z:BYc_z:east_west_c"
B4_focal <- "PGI_Edu_z:BYc_z:east_west_c:gender_c"

# Parental-education control (PGI-agnostic, so the focal PGI interaction is
# isolated when the PGI is swapped): ParEdu main + ParEdu x region + ParEdu x
# birth year (+ their three-way). Mirrors the B1 ParEdu adjustment minus the
# ParEdu x PGI terms (which would not be comparable across swapped PGIs).
paredu_ctrl <- paste("+ parental_edu_z_kernel + parental_edu_z_kernel:east_west_c",
                     "+ parental_edu_z_kernel:BYc_z + parental_edu_z_kernel:BYc_z:east_west_c")

# Alternative PGIs available on the mobility sample.
alt_pgis <- c("PGI_Cog_z", "PGI_nonCog_z")
alt_pgis <- alt_pgis[vapply(alt_pgis, function(v)
  v %in% names(dat_mob) && !all(is.na(dat_mob[[v]])), logical(1))]
s1_pgis <- c("PGI_Edu_z", alt_pgis)

# =====================================================================
# S1-B: alternative-PGI family (BH per family) — B1 focal + B4 four-way
# =====================================================================
cat("S1-B: alternative-PGI families...\n")
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
S1_B1 <- .s1_family(dat_mob,
  function(pgi) paste("mobility ~ (", pgi, "+ BYc_z + east_west_c + gender_c)^3", paredu_ctrl),
  function(pgi) paste0(pgi, ":BYc_z:east_west_c"), "B1 (PGI x BY x region)")
S1_B4 <- .s1_family(dat_mob,
  function(pgi) paste0("mobility ~ ", pgi, " * BYc_z * east_west_c * gender_c",
                       " + parental_edu_z_kernel * BYc_z * east_west_c * gender_c"),
  function(pgi) paste0(pgi, ":BYc_z:east_west_c:gender_c"), "B4 (4-way PGI x BY x region x gender)")
S1_AltPGIs <- dplyr::bind_rows(S1_B1, S1_B4)

# =====================================================================
# S2-B: height negative control on mobility (B1 spec) + sanity
# =====================================================================
cat("S2-B: height negative control on mobility...\n")
S2_Height_Mobility <- S2_Height_Height <- NULL
if ("PGI_Height_z" %in% names(dat_mob) && !all(is.na(dat_mob$PGI_Height_z))) {
  m_s2 <- lm(as.formula(paste("mobility ~ (PGI_Height_z + BYc_z + east_west_c + gender_c)^3", paredu_ctrl)),
             data = dat_mob)
  S2_Height_Mobility <- tidy_coefs(m_s2, "S2_Height_Mobility (B1 spec)")
  if ("height_z" %in% names(dat_mob) && !all(is.na(dat_mob$height_z))) {
    m_s2h <- lm(height_z ~ (PGI_Height_z + BYc_z + east_west_c + gender_c)^3,
                data = dat_mob |> dplyr::filter(!is.na(height_z)))
    S2_Height_Height <- tidy_coefs(m_s2h, "S2_Height_Height (sanity)")
  }
}

# =====================================================================
# Read B_Mobility_Models.rds (S3-B gate)
# =====================================================================
# GATE_RUN_TS re-points the gate cache at another run, for the case where only
# this stage is re-run into a fresh run folder: the gate then reads the run the
# gated decision was made in, instead of degrading to "unavailable".
B_GATE_TS  <- Sys.getenv("GATE_RUN_TS", RUN_TS_VAL)
B_GATE_DIR <- file.path(dirname(dirname(OUT_DIR_B)), paste0("RUN_", B_GATE_TS), "B_Mobility")
B_RDS <- file.path(B_GATE_DIR, "B_Mobility_Models.rds")
B_mc  <- if (file.exists(B_RDS)) tryCatch(readRDS(B_RDS), error = function(e) NULL) else NULL
gxe_b_sig <- NA; b_gate_basis <- "unavailable"
if (!is.null(B_mc)) {
  b1_p <- { r <- B_mc$B1_Coefficients[B_mc$B1_Coefficients$term == B1_focal, ]
            if (nrow(r)) r$p.value[1] else NA_real_ }
  b3_p <- { sm <- B_mc$B3_Smooths_pgi_tp
            if (!is.null(sm) && nrow(sm)) suppressWarnings(min(sm[["p-value"]], na.rm = TRUE)) else NA_real_ }
  gxe_b_sig <- isTRUE(b1_p < 0.05) || isTRUE(b3_p < 0.05)
  b_gate_basis <- sprintf("RE B1 focal p=%s; min B3 smooth p=%s", signif(b1_p, 3), signif(b3_p, 3))
} else {
  m_b1_ols <- lm(as.formula(paste("mobility ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3", paredu_ctrl)),
                 data = dat_mob)
  b1_p_ols <- coef(summary(m_b1_ols))[B1_focal, 4]
  gxe_b_sig <- isTRUE(b1_p_ols < 0.05)
  b_gate_basis <- sprintf("B_RDS absent; fallback dedup-OLS B1 focal p=%s (B3 not checked)", signif(b1_p_ols, 3))
}
cat(sprintf("B gate: GxE_B_sig = %s [%s]\n", gxe_b_sig, b_gate_basis))

# =====================================================================
# S1-B (finding-matched, random effects) — Plan_deviations.md §15
# Alternative PGIs (Education / Cognition / Noncognitive) on the ACTUAL focal
# mobility findings, refit with the primary family-RE specification on a COMMON
# complete-case sample across the three PGIs (each re-residualised on
# crossPC1-10 and re-z-scored within that sample; PGI-Education refit on the
# same sample). Extends the dedup-OLS S1 above (B1 three-way + B4 four-way) to:
#   Finding 2  East-West mobility contrast PGI x region,
#              with East/West standardised slopes — B0-equivalent.
#   Finding 3  pooled PGI x gender + gender-specific slopes — B4-equivalent.
#   B5 nonlinear gender: omnibus + female/male curvature + the female-minus-male
#              difference smooth (simultaneous 95% band), same procedure as the
#              primary PGI-Education analysis.
# BH within each focal hypothesis across the three PGIs. Specs derived from the
# cached primary fits (B_mc) so alt-PGI models are structurally identical bar
# the PGI variable and the parental-education control is preserved.
# =====================================================================
cat("S1-B (finding-matched, RE): common-sample alternative PGIs...\n")
S1fm_raw_b <- c("PGI_Edu", "PGI_Cog", "PGI_nonCog")
S1fm_raw_b <- S1fm_raw_b[vapply(S1fm_raw_b, function(v)
  exists("dat_cluster_mob") && v %in% names(dat_cluster_mob) &&
    !all(is.na(dat_cluster_mob[[v]])), logical(1))]
S1fm_ok_b <- !is.null(B_mc) && exists("dat_cluster_mob") &&
  all(paste0("crossPC", 1:10) %in% names(dat_cluster_mob)) && length(S1fm_raw_b) >= 2
S1_finding2 <- S1_finding2_slopes <- S1_finding3_mob <- S1_finding3_mob_slopes <-
  S1_B5_nonlin <- S1_B5_nonlin_diffsmooth <- S1_B5_nonlin_diff_intervals <- NULL
S1fm_mob_N <- NA_integer_; S1fm_status_b <- "not_run"
if (S1fm_ok_b) {
  S1fm_mob <- altpgi_common_sample(dat_cluster_mob, pgis_raw = S1fm_raw_b)
  # gender-specific ParEdu cells for the per-gender GAM (the linear/interaction
  # models use parental_edu_z_kernel directly, already present).
  S1fm_mob$ParEdu_F <- S1fm_mob$parental_edu_z_kernel * (1L - S1fm_mob$gender_01)
  S1fm_mob$ParEdu_M <- S1fm_mob$parental_edu_z_kernel *        S1fm_mob$gender_01
  S1fm_pgis_b <- paste0(S1fm_raw_b, "_z")
  S1fm_mob_N  <- nrow(S1fm_mob)
  b0_tmpl <- gsub("PGI_Edu_z", "PGI", rhs_without_re(B_mc$m_B0_re), fixed = TRUE)
  b4_tmpl <- gsub("PGI_Edu_z", "PGI", rhs_without_re(B_mc$m_B4_re), fixed = TRUE)
  # Finding 2: East-West mobility contrast + region slopes.
  S1_finding2 <- altpgi_re_family(S1fm_mob, b0_tmpl,
    list(`PGI x region` = "PGI:east_west_c"),
    "S1_finding2 (East-West mobility contrast)", pgis = S1fm_pgis_b)
  S1_finding2_slopes <- altpgi_slopes_family(S1fm_mob, b0_tmpl,
    list(West = list(east_west_c = -0.5), East = list(east_west_c = 0.5)),
    "S1_finding2 (East-West mobility contrast) slopes", pgis = S1fm_pgis_b)
  # Finding 3: mobility gender + gender-specific slopes.
  S1_finding3_mob <- altpgi_re_family(S1fm_mob, b4_tmpl,
    list(`PGI x gender` = "PGI:gender_c"),
    "S1_finding3 (mobility gender)", pgis = S1fm_pgis_b)
  S1_finding3_mob_slopes <- altpgi_slopes_family(S1fm_mob, b4_tmpl,
    list(Female = list(gender_c = 0.5), Male = list(gender_c = -0.5)),
    "S1_finding3 (mobility gender) slopes", pgis = S1fm_pgis_b)
  # B5 nonlinear gender reproduction (emphasis): omnibus + curvature + smooth.
  if (!is.null(B_mc$m_B5_lin_re) && !is.null(B_mc$m_B5_nonlin_re)) {
    .b5fm <- altpgi_pergender_gam(S1fm_mob,
      rhs_without_re(B_mc$m_B5_lin_re), rhs_without_re(B_mc$m_B5_nonlin_re),
      "S1_B5_nonlin (mobility per-gender GAM)",
      pgis = S1fm_pgis_b, want_diffsmooth = TRUE)
    S1_B5_nonlin                <- .b5fm$tests
    S1_B5_nonlin_diffsmooth     <- .b5fm$diffsmooth
    S1_B5_nonlin_diff_intervals <- .b5fm$diff_intervals
  }
  S1fm_status_b <- sprintf("run on common N=%d; PGIs=%s", S1fm_mob_N,
                           paste(sub("_z$", "", S1fm_pgis_b), collapse = "/"))
} else {
  S1fm_status_b <- "skipped (B_mc, raw PGIs or crossPCs unavailable)"
}
S1_FindingMatched_Status <- data.frame(
  section = "B", status = S1fm_status_b, common_N = S1fm_mob_N,
  bh_families = "BH within each focal hypothesis across PGI-Education/Cognition/Noncognitive",
  stringsAsFactors = FALSE)
cat(sprintf("  S1-B finding-matched: %s\n", S1fm_status_b))

# =====================================================================
# S3-B: heteroscedasticity (gaulss), per region, on mobility.
# gaulss cannot carry the family RE, so dat_cluster_mob is used WITHOUT it
# (documented caveat; ../../Plan_deviations.md §12).
# =====================================================================
cat("S3-B: heteroscedasticity (gaulss)...\n")
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
# B0 East-West PGI-mobility contrast under heteroscedasticity (OLS vs gaulss).
.b0_contrast_under_het <- function(dat) {
  b0_mean <- paste("mobility ~ PGI_Edu_z * east_west_c +",
    "east_west_c * (gender_c + BYc_z + parental_edu_z_kernel) +",
    "gender_c * BYc_z + gender_c * parental_edu_z_kernel + BYc_z * parental_edu_z_kernel")
  foc <- "PGI_Edu_z:east_west_c"; rows <- list()
  for (sc in list(c("full", "all"), c("pre-reunification", "pre"))) {
    d <- if (sc[2] == "all") dat else dat[dat$reunif == "pre", , drop = FALSE]
    if (nrow(d) < 100 || length(unique(d$east_west)) < 2) next
    m_ols <- lm(as.formula(b0_mean), data = d)
    m_het <- tryCatch(mgcv::gam(list(as.formula(b0_mean), ~ PGI_Edu_z * east_west_c + BYc_z),
                                family = mgcv::gaulss(), data = d, method = "REML"), error = function(e) NULL)
    co <- coef(summary(m_ols))
    rows[[paste0(sc[1], "_ols")]] <- data.frame(scope = sc[1], model = "OLS (homoscedastic)", N = nrow(d),
      estimate = co[foc, 1], se = co[foc, 2], p = co[foc, 4], stringsAsFactors = FALSE)
    if (!is.null(m_het)) { pt <- summary(m_het)$p.table
      rows[[paste0(sc[1], "_het")]] <- data.frame(scope = sc[1], model = "gaulss (heteroscedastic)", N = nrow(d),
        estimate = pt[foc, 1], se = pt[foc, 2], p = pt[foc, "Pr(>|z|)"], stringsAsFactors = FALSE) }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}
B_S3 <- NULL
if (isTRUE(gxe_b_sig)) {
  pe_smooth <- .sm("parental_edu_z_kernel")
  mean_lin  <- paste("BYc + gender_c + PGI_Edu_z + PGI_Edu_z:BYc + parental_edu_z_kernel +", pe_smooth)
  mean_va   <- paste("BYc + gender_c + PGI_Edu_z + parental_edu_z_kernel +", .sm("PGI_Edu_z"), "+", pe_smooth)
  fits <- lapply(c("East", "West"), function(reg)
    .fit_region(dat_cluster_mob[as.character(dat_cluster_mob$east_west) == reg, , drop = FALSE], reg,
                "mobility", mean_lin, mean_va,
                ~ BYc + PGI_Edu_z + BYc:PGI_Edu_z, ~ BYc + PGI_Edu_z + BYc:PGI_Edu_z,
                ~ BYc + PGI_Edu_z + BYc:PGI_Edu_z))
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
  B_S3 <- list(
    dispersion = if (length(disp)) do.call(rbind, disp) else data.frame(note = "No dispersion summaries."),
    focal_under_het = if (length(foc)) do.call(rbind, foc) else data.frame(note = "No focal smooth extracted."),
    b0_contrast = .b0_contrast_under_het(dat_cluster_mob),
    region_status = data.frame(region = names(notes), status = unname(notes), stringsAsFactors = FALSE))
}
S3_Status <- data.frame(section = "B-S3 (mobility)", GxE_significant = gxe_b_sig,
  basis = b_gate_basis, triggered = isTRUE(gxe_b_sig),
  note = if (isTRUE(gxe_b_sig)) "S3-B heteroscedasticity GAMs fitted." else "Not triggered; S3-B not fitted.",
  stringsAsFactors = FALSE)

# =====================================================================
# Focal dispersion robustness (per-focal gated): for each Section B focal
# that is significant in the main RE analysis, re-estimate the mean focal
# allowing the residual variance to follow the same structure (gaulss, no RE),
# parental education controlled. (The B0 Region x PGI focal is also covered in
# more detail, with the pre-reunification split, by B_S3$b0_contrast above.)
# =====================================================================
cat("Focal dispersion robustness (gated on main-analysis significance)...\n")
.cp <- function(tbl, term) { if (is.null(tbl)) return(NA_real_)
  r <- tbl[tbl$term == term, ]; if (nrow(r)) r$p.value[1] else NA_real_ }
.pe <- "parental_edu_z_kernel"
Focal_Dispersion <- if (!is.null(B_mc) && exists("dat_cluster_mob")) focal_dispersion_check(
  dat_cluster_mob, "mobility", list(
    list(label = "PGI x region (overall, B0)",
         mean_rhs = paste("PGI_Edu_z*east_west_c + gender_c + BYc_z +", .pe), focal = "PGI_Edu_z:east_west_c",
         scale_rhs = "PGI_Edu_z*east_west_c", main_p = .cp(B_mc$B0_Coefficients, "PGI_Edu_z:east_west_c")),
    list(label = "PGI x birth year x region (B1)",
         mean_rhs = paste("PGI_Edu_z*BYc_z*east_west_c + gender_c +", .pe), focal = "PGI_Edu_z:BYc_z:east_west_c",
         scale_rhs = "PGI_Edu_z*BYc_z*east_west_c", main_p = .cp(B_mc$B1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c")),
    list(label = "PGI x gender (pooled across region)",
         mean_rhs = paste("PGI_Edu_z*gender_c + east_west_c + BYc_z +", .pe), focal = "PGI_Edu_z:gender_c",
         scale_rhs = "PGI_Edu_z*gender_c", main_p = .cp(B_mc$B4_Coefficients, "PGI_Edu_z:gender_c")),
    list(label = "PGI x region x gender",
         mean_rhs = paste("PGI_Edu_z*east_west_c*gender_c + BYc_z +", .pe), focal = "PGI_Edu_z:east_west_c:gender_c",
         scale_rhs = "PGI_Edu_z*east_west_c*gender_c", main_p = .cp(B_mc$B4_Coefficients, "PGI_Edu_z:east_west_c:gender_c")),
    list(label = "PGI x birth year x region x gender (B4)",
         mean_rhs = paste("PGI_Edu_z*BYc_z*east_west_c*gender_c +", .pe), focal = "PGI_Edu_z:BYc_z:east_west_c:gender_c",
         scale_rhs = "PGI_Edu_z*east_west_c*gender_c", main_p = .cp(B_mc$B4_Coefficients, "PGI_Edu_z:BYc_z:east_west_c:gender_c")))) else NULL

# =====================================================================
# Save
# =====================================================================
model_cache <- list(
  RUN_TS = RUN_TS_VAL, s1_pgis = s1_pgis,
  S1_AltPGIs = S1_AltPGIs,
  S1_FindingMatched_Status = S1_FindingMatched_Status,
  S1_finding2 = S1_finding2, S1_finding2_slopes = S1_finding2_slopes,
  S1_finding3_mob = S1_finding3_mob, S1_finding3_mob_slopes = S1_finding3_mob_slopes,
  S1_B5_nonlin = S1_B5_nonlin, S1_B5_nonlin_diffsmooth = S1_B5_nonlin_diffsmooth,
  S1_B5_nonlin_diff_intervals = S1_B5_nonlin_diff_intervals,
  S2_Height_Mobility = S2_Height_Mobility, S2_Height_Height = S2_Height_Height,
  S3_Status = S3_Status, B_S3 = B_S3, Focal_Dispersion = Focal_Dispersion,
  N_mob = nrow(dat_mob), N_cluster_mob = if (exists("dat_cluster_mob")) nrow(dat_cluster_mob) else 0L)
OUT_RDS <- file.path(OUT_DIR_B, "B_Sensitivities_Models.rds")
saveRDS(model_cache, OUT_RDS)
cat(sprintf("B_sensitivities_analysis: models RDS saved to %s\n", OUT_RDS))
