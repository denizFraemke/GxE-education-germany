# A_appendix_analysis.R
# Section A per-study appendix (attainment) — analysis stage.
# Appendix-grade per-study education analyses owned by Section A (D, F,
# F-linear; the mobility per-study analysis E lives in Section B). All
# dedup-OLS / descriptive; within-study birth-year windows are narrow, so
# several study x cell combos are skipped for lack of N. The better-powered
# heterogeneity check is leave-one-study-out in the main Section A analyses.
# See ../../../Plan_deviations.md §13.
#
#   D  Per-study A3-style education GAM (region-specific PGI cohort smooths).
#   F  Per-study A6-style Region x Gender GAM on edu_z_kernel.
#   F-linear  Per-study robust linear gender model (PGI main, gender x PGI,
#             Region x Gender x PGI).
#
# Usage (from run_A.command, after A_analysis.R): Rscript A_Education/A_appendix_analysis.R

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
  stop("Cannot determine A_appendix_analysis.R location.")
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
cat(sprintf("A_appendix_analysis: dat_edu N=%d\n", nrow(dat_edu)))

.cohorts_in <- function(dat) sort(unique(as.character(dat$cohort)))
.clip_curve <- function(crv, by_obs, lo = 0.02, hi = 0.98) {
  if (is.null(crv) || !length(by_obs) || all(is.na(by_obs))) return(crv)
  qs <- quantile(by_obs, c(lo, hi), na.rm = TRUE)
  crv[crv$birth_year >= qs[1] & crv$birth_year <= qs[2], , drop = FALSE]
}

# ---- D: per-study region GAM (A3-style education) ----
pc_region_gam <- function(dat, outcome, covars = "", min_n = 100, min_reg = 30) {
  rows <- list(); curves <- list(); smooths <- list()
  for (coh in .cohorts_in(dat)) {
    d <- dat[as.character(dat$cohort) == coh, , drop = FALSE]
    keep <- !is.na(d[[outcome]]) & !is.na(d$PGI_Edu_z) & !is.na(d$birth_year) & !is.na(d$east)
    d <- d[keep, , drop = FALSE]
    n <- nrow(d); ne <- sum(d$east == 1); nw <- sum(d$east == 0)
    base <- data.frame(dataset = coh, N = n, N_west = nw, N_east = ne, stringsAsFactors = FALSE)
    mkrow <- function(fit, wp = NA, ep = NA, nl = NA)
      cbind(base, data.frame(fit = fit, West_smooth_p = wp, East_smooth_p = ep, nonlinear_LRT_p = nl))
    if (n < min_n) { rows[[coh]] <- mkrow("skip (N<100)"); next }
    nby <- length(unique(d$birth_year)); k <- min(K_DEFAULT, max(4L, floor(nby / 3)))
    bym <- mean(d$birth_year); d$BYc <- d$birth_year - bym; attr(d, "by_mean") <- bym
    d$PGI_W <- d$PGI_Edu_z * (1L - d$east); d$PGI_E <- d$PGI_Edu_z * d$east
    cov <- if (nzchar(covars)) paste("+", covars) else ""
    both <- nw >= min_reg && ne >= min_reg
    if (both) {
      rhs_lin <- sprintf("%s ~ BYc + east_west_c + gender_c %s + PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc", outcome, cov)
      rhs_nl  <- sprintf("%s + s(BYc, by = PGI_W, k = %d, bs = '%s') + s(BYc, by = PGI_E, k = %d, bs = '%s')", rhs_lin, k, BS_DEFAULT, k, BS_DEFAULT)
      regs <- c("West", "East")
    } else {
      reg <- if (nw >= min_reg) "West" else if (ne >= min_reg) "East" else NA_character_
      if (is.na(reg)) { rows[[coh]] <- mkrow("skip (no region >= 30)"); next }
      pc <- if (reg == "West") "PGI_W" else "PGI_E"
      rhs_lin <- sprintf("%s ~ BYc + gender_c %s + %s + %s:BYc", outcome, cov, pc, pc)
      rhs_nl  <- sprintf("%s + s(BYc, by = %s, k = %d, bs = '%s')", rhs_lin, pc, k, BS_DEFAULT)
      regs <- reg
    }
    m_lin <- tryCatch(mgcv::gam(as.formula(rhs_lin), data = d, method = "REML"), error = function(e) NULL)
    m_nl  <- tryCatch(mgcv::gam(as.formula(rhs_nl),  data = d, method = "REML"), error = function(e) NULL)
    if (is.null(m_nl)) { rows[[coh]] <- mkrow("GAM failed"); next }
    st <- as.data.frame(summary(m_nl)$s.table)
    smp <- function(pat) { r <- st[grepl(pat, rownames(st), fixed = TRUE), , drop = FALSE]; if (nrow(r)) signif(r[1, "p-value"], 3) else NA_real_ }
    nl_p <- tryCatch({ v <- lrt_pair(m_lin, m_nl); signif(v[["p"]], 3) }, error = function(e) NA_real_)
    smooths[[coh]] <- dplyr::mutate(tibble::rownames_to_column(st, "smooth_term"), dataset = coh, .before = 1)
    for (reg in regs) {
      crv <- tryCatch(extract_pgi_slope(m_nl, d, east_west = reg, n_grid = 150), error = function(e) NULL)
      if (!is.null(crv)) {
        crv <- .clip_curve(crv, d$birth_year[if (reg == "East") d$east == 1 else d$east == 0])
        crv$dataset <- coh; curves[[paste(coh, reg)]] <- crv
      }
    }
    rows[[coh]] <- mkrow(if (both) "2-region" else paste0(regs, "-only"),
                         wp = smp("PGI_W"), ep = smp("PGI_E"), nl = nl_p)
  }
  list(table = dplyr::bind_rows(rows),
       curves = if (length(curves)) dplyr::bind_rows(curves) else NULL,
       smooths = if (length(smooths)) dplyr::bind_rows(smooths) else NULL)
}

# ---- F: per-study Region x Gender GAM on edu_z_kernel ----
pc_rxg_gam <- function(dat, outcome, min_cell = 50) {
  rows <- list(); curves <- list(); smooths <- list()
  for (coh in .cohorts_in(dat)) {
    d <- dat[as.character(dat$cohort) == coh, , drop = FALSE]
    keep <- !is.na(d[[outcome]]) & !is.na(d$PGI_Edu_z) & !is.na(d$east) & !is.na(d$gender_01)
    d <- d[keep, , drop = FALSE]
    cn <- table(factor(ifelse(d$east == 1, "East", "West"), c("West", "East")),
                factor(ifelse(d$gender_01 == 1, "male", "female"), c("female", "male")))
    base <- data.frame(dataset = coh, N = nrow(d),
                       Wf = cn["West", "female"], Wm = cn["West", "male"],
                       Ef = cn["East", "female"], Em = cn["East", "male"], stringsAsFactors = FALSE)
    mkrow <- function(fit) cbind(base, data.frame(fit = fit))
    west_ok <- all(cn["West", ] >= min_cell); east_ok <- all(cn["East", ] >= min_cell)
    four <- west_ok && east_ok
    single <- (west_ok && all(cn["East", ] == 0)) || (east_ok && all(cn["West", ] == 0))
    if (!four && !single) { rows[[coh]] <- mkrow("skip (cell N)"); next }
    nby <- length(unique(d$birth_year)); k <- min(K_DEFAULT, max(4L, floor(nby / 4)))
    bym <- mean(d$birth_year); d$BYc <- d$birth_year - bym; attr(d, "by_mean") <- bym
    d$PGI_WF <- d$PGI_Edu_z * (1L - d$east) * (1L - d$gender_01)
    d$PGI_WM <- d$PGI_Edu_z * (1L - d$east) *        d$gender_01
    d$PGI_EF <- d$PGI_Edu_z *        d$east  * (1L - d$gender_01)
    d$PGI_EM <- d$PGI_Edu_z *        d$east  *        d$gender_01
    if (four) {
      f <- sprintf("%s ~ BYc + east_west_c + gender_c + PGI_WF + PGI_WM + PGI_EF + PGI_EM + PGI_WF:BYc + PGI_WM:BYc + PGI_EF:BYc + PGI_EM:BYc + s(BYc, by = PGI_WF, k = %d, bs = '%s') + s(BYc, by = PGI_WM, k = %d, bs = '%s') + s(BYc, by = PGI_EF, k = %d, bs = '%s') + s(BYc, by = PGI_EM, k = %d, bs = '%s')",
                   outcome, k, BS_DEFAULT, k, BS_DEFAULT, k, BS_DEFAULT, k, BS_DEFAULT)
      cells <- list(c("West", "female"), c("West", "male"), c("East", "female"), c("East", "male"))
    } else {
      reg <- if (west_ok) "West" else "East"
      cc  <- if (reg == "West") c("PGI_WF", "PGI_WM") else c("PGI_EF", "PGI_EM")
      f <- sprintf("%s ~ BYc + gender_c + %s + %s + %s:BYc + %s:BYc + s(BYc, by = %s, k = %d, bs = '%s') + s(BYc, by = %s, k = %d, bs = '%s')",
                   outcome, cc[1], cc[2], cc[1], cc[2], cc[1], k, BS_DEFAULT, cc[2], k, BS_DEFAULT)
      cells <- if (reg == "West") list(c("West", "female"), c("West", "male")) else list(c("East", "female"), c("East", "male"))
    }
    m <- tryCatch(mgcv::gam(as.formula(f), data = d, method = "REML"), error = function(e) NULL)
    if (is.null(m)) { rows[[coh]] <- mkrow("GAM failed"); next }
    smooths[[coh]] <- dplyr::mutate(tibble::rownames_to_column(as.data.frame(summary(m)$s.table), "smooth_term"), dataset = coh, .before = 1)
    for (cs in cells) {
      crv <- tryCatch(extract_pgi_slope(m, d, east_west = cs[1], gender = cs[2], n_grid = 150), error = function(e) NULL)
      if (!is.null(crv)) {
        mask <- (if (cs[1] == "East") d$east == 1 else d$east == 0) &
                (if (cs[2] == "male") d$gender_01 == 1 else d$gender_01 == 0)
        crv <- .clip_curve(crv, d$birth_year[mask])
        crv$dataset <- coh; curves[[paste(coh, cs[1], cs[2])]] <- crv
      }
    }
    rows[[coh]] <- mkrow(if (four) "4-cell" else paste0(reg, "-2cell"))
  }
  list(table = dplyr::bind_rows(rows),
       curves = if (length(curves)) dplyr::bind_rows(curves) else NULL,
       smooths = if (length(smooths)) dplyr::bind_rows(smooths) else NULL)
}

# ---- F-linear: per-study robust linear gender model ----
pc_gender_linear <- function(dat, outcome, paredu = FALSE, min_n = 80) {
  rows <- list()
  for (coh in .cohorts_in(dat)) {
    d <- dat[as.character(dat$cohort) == coh, , drop = FALSE]
    keep <- !is.na(d[[outcome]]) & !is.na(d$PGI_Edu_z) & !is.na(d$gender_c) & !is.na(d$east_west_c) & !is.na(d$BYc_z)
    if (paredu) keep <- keep & !is.na(d$parental_edu_z_kernel)
    d <- d[keep, , drop = FALSE]
    both <- length(unique(d$east_west)) == 2
    g <- function(m, t) { cf <- summary(m)$coefficients; if (t %in% rownames(cf)) sprintf("%.4f (p=%.3g)%s", cf[t, 1], cf[t, 4], ifelse(cf[t, 4] < .05, "*", "")) else "NA" }
    if (nrow(d) < min_n) {
      rows[[coh]] <- data.frame(dataset = coh, N = nrow(d), regions = "skip",
        PGI_main = NA, gender_x_PGI = NA, RegionxGenderxPGI = NA, stringsAsFactors = FALSE); next
    }
    rhs <- if (both) "PGI_Edu_z * east_west_c * gender_c + BYc_z" else "PGI_Edu_z * gender_c + BYc_z"
    if (paredu) rhs <- paste(rhs, "+ parental_edu_z_kernel")
    m <- lm(as.formula(paste(outcome, "~", rhs)), data = d)
    rows[[coh]] <- data.frame(dataset = coh, N = nrow(d),
      regions = if (both) "both" else paste(unique(as.character(d$east_west)), collapse = ""),
      PGI_main = g(m, "PGI_Edu_z"), gender_x_PGI = g(m, "PGI_Edu_z:gender_c"),
      RegionxGenderxPGI = if (both) g(m, "PGI_Edu_z:east_west_c:gender_c") else NA,
      stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}

cat("D: per-study education GAMs...\n");           D <- pc_region_gam(dat_edu, "edu_z_kernel")
cat("F: per-study Region x Gender GAMs...\n");      Fres <- pc_rxg_gam(dat_edu, "edu_z_kernel")
cat("F-linear: per-study linear gender model (edu)...\n"); F_linear_edu <- pc_gender_linear(dat_edu, "edu_z_kernel")

model_cache <- list(RUN_TS = RUN_TS_VAL, N_edu = nrow(dat_edu),
                    D = D, F = Fres, F_linear_edu = F_linear_edu)
OUT_RDS <- file.path(OUT_DIR_A, "A_Appendix_Models.rds")
saveRDS(model_cache, OUT_RDS)
cat(sprintf("A_appendix_analysis: models RDS saved to %s\n", OUT_RDS))
