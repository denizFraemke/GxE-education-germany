# ancestry_sensitivity.R
# Sample crossed with the scope of the ancestry-PC residualisation, a 2x2 for the
# attainment PGI-Education x birth year x region term.
#
# Sample and ancestry adjustment are otherwise confounded: the pooled estimate
# uses across-study PCs, while the SOEP-G estimate that reproduces the earlier
# pattern reconstructs within-SOEP residualisation. Crossing SOEP-G against the
# pooled sample with within-study against across-study PCs separates them.
#
# Three cells are read from the frozen run:
#
#   pooled / across-study PCs   attainment/main/A1_Coefficients.csv (A1_RE),
#                               cross-checked against focal_estimates.csv
#   SOEP-G / within-study PCs   figS05_soep_threeway.csv, row 1
#   SOEP-G / across-study PCs   figS05_soep_threeway.csv, row 2
#
# Each is verified by refitting on the rebuilt sample; the script stops unless
# all three reproduce.
#
# The fourth cell, pooled sample with within-study residualisation, is the one
# refit performed here: the PGI is residualised on the within-cohort PCs (PC1-10)
# separately within each study, z-scored once across the pooled attainment
# sample, and the published A1 specification fitted otherwise unchanged.
#
# "Within each study" means per genotyping batch. SHIP-0 and SHIP-Td carry
# separate PCAs (separate *_PC_eigenval.tsv under
# 01_Genotype/Harmonized/data/final/), so PC1 differs between them and pooling
# them in one regression would residualise on a mixture. Five groups: BASEII,
# SHIP0, SHIPTd, SOEP, TwinLife. The four-group alternative (SHIP pooled) is
# fitted too and recorded in ancestry_checks.csv.
#
# Nothing here overwrites the frozen run and its model caches are not loaded
# (~650 MB each). The analytic sample is rebuilt once from the merged data
# through 00_setup/, as 05_manuscript_export/check_frozen_run.R step 4 does. No
# GAM is refit. Aggregate output only; no individual row is printed or written.
#
# Usage:
#   RUN_TS=<ts> Rscript 05_manuscript_export/check_frozen_run.R
#   RUN_TS=<ts> Rscript 06_Sensitivity_analyses/ancestry_sensitivity.R
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
  stop("Cannot determine ancestry_sensitivity.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
REPO_ROOT    <- normalizePath(file.path(WORKFLOW_DIR, "..", ".."), mustWork = TRUE)
rm(.find_this_file)


# ---------------------------------------------------------------------------
# Constants: the published specification, verbatim from A_analysis.R
# ---------------------------------------------------------------------------
FROZEN_RUN_TS <- "20260729_1555"
Z975          <- stats::qnorm(0.975)
Z_FIGS05      <- 1.96          # figS05 wrote its CIs on 1.96, not qnorm(.975)

A1_RHS  <- "(PGI_Edu_z + BYc_z + east_west_c + gender_c)^3"
RE_TERM <- "s(fid_re, bs = 're')"
FOCAL   <- "PGI_Edu_z:BYc_z:east_west_c"

WITHIN_PCS <- paste0("PC", 1:10)        # per-cohort PCA basis (figS05 uses these)
CROSS_PCS  <- paste0("crossPC", 1:10)   # shared cross-cohort basis (canonical)

# The genotyping batches the within-cohort PCs were computed in. SHIP is two.
BATCH_FILES <- c(SHIP0 = "SHIP0_PGI_PCs.tsv", SHIPTd = "SHIPTD_PGI_PCs.tsv")

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

# z_scale(), verbatim from R/helpers_formatting.R (NA-preserving, sd == 0 -> NA).
z_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

fit_bam <- function(lhs, rhs, data) {
  mgcv::bam(stats::as.formula(paste(lhs, "~", rhs, "+", RE_TERM)),
            data = data, method = "fREML", discrete = TRUE)
}

# One parametric coefficient: estimate, SE, statistic, p from the p.table (as
# tidy_coefs()/broom read them), normal CI on qnorm(.975) as the frozen export
# carries it.
p_row <- function(m, term) {
  pt <- summary(m)$p.table
  if (!term %in% rownames(pt)) stop("term not in model: ", term)
  est <- unname(pt[term, 1]); se <- unname(pt[term, 2])
  data.frame(term = term, estimate = est, se = se, statistic = unname(pt[term, 3]),
             p = unname(pt[term, 4]), ci_lo = est - Z975 * se, ci_hi = est + Z975 * se,
             stringsAsFactors = FALSE)
}

# Number of rows a formula actually uses (complete cases over its variables).
n_used <- function(fml_string, data) {
  v <- intersect(all.vars(stats::as.formula(fml_string)), names(data))
  sum(stats::complete.cases(data[, v, drop = FALSE]))
}

# Checks collector. kind = "check" gates the run; kind = "diagnostic" records an
# aggregate for the record and never fails.
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
# Ancestry residualisation
# ---------------------------------------------------------------------------
# Residuals of `pgi` on `pcs`, fitted SEPARATELY within each level of `group`.
# NA where the PGI or any PC is missing, or where a group is too small to fit.
# The per-group intercept absorbs each group's mean PGI, so this also demeans
# the PGI by study, a property of within-study residualisation that is recorded
# in the table's note column.
resid_pgi_within_groups <- function(dat, pgi, pcs, group) {
  stopifnot(all(c(pgi, pcs) %in% names(dat)), length(group) == nrow(dat))
  fml <- stats::as.formula(paste(pgi, "~", paste(pcs, collapse = " + ")))
  r <- rep(NA_real_, nrow(dat))
  for (g in sort(unique(group[!is.na(group)]))) {
    idx <- which(group == g)
    ok  <- idx[stats::complete.cases(dat[idx, c(pgi, pcs), drop = FALSE])]
    if (length(ok) <= length(pcs) + 1L) next
    r[ok] <- stats::residuals(stats::lm(fml, data = dat[ok, , drop = FALSE]))
  }
  r
}

# Residuals of `pgi` on `pcs` over the whole frame (the canonical pooled form;
# resid_pgi_on_xpcs() in 00_setup/prepare_samples.R with pcs = crossPC1-10).
resid_pgi_pooled <- function(dat, pgi, pcs) {
  resid_pgi_within_groups(dat, pgi, pcs, rep("all", nrow(dat)))
}

# Genotyping-batch labels: cohort, except SHIP, split into SHIP0 / SHIPTd by
# matching the preserved genotype-side IID against the two per-batch PGI+PC
# tables. IIDs are pseudonymous identifiers and are never printed or written.
genotyping_batch <- function(dat, geno_final_dir) {
  grp <- as.character(dat$cohort)
  info <- list(n_ship = sum(grp == "SHIP", na.rm = TRUE), n_matched = 0L,
               n_unmatched = 0L, n_ambiguous = 0L, available = TRUE)
  if (info$n_ship == 0L) return(list(group = grp, info = info))
  read_iids <- function(f) {
    p <- file.path(geno_final_dir, f)
    if (!file.exists(p)) return(NULL)
    x <- utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
    v <- as.character(x$IID); rm(x); v
  }
  i0  <- read_iids(BATCH_FILES[["SHIP0"]])
  itd <- read_iids(BATCH_FILES[["SHIPTd"]])
  if (is.null(i0) || is.null(itd)) {
    info$available <- FALSE
    return(list(group = grp, info = info))
  }
  iid <- as.character(dat$IID)
  is_ship <- !is.na(grp) & grp == "SHIP"
  in0  <- is_ship & iid %in% i0
  intd <- is_ship & iid %in% itd
  info$n_ambiguous <- sum(in0 & intd)
  grp[in0  & !intd] <- "SHIP0"
  grp[intd & !in0]  <- "SHIPTd"
  info$n_matched   <- sum(in0 | intd)
  info$n_unmatched <- sum(is_ship & !(in0 | intd))
  rm(i0, itd, iid)
  list(group = grp, info = info)
}

# ---------------------------------------------------------------------------
# The frozen SOEP cells, reconstructed exactly as
# 05_manuscript_export/export_manuscript_aggregates.R writes figS05.
# ---------------------------------------------------------------------------
# Properties of that block, unchanged here: plain mgcv::bam with
# method = "REML" and no family random effect; CIs on 1.96; the "within" PGI
# residualised from PGI_Edu_raw on the SOEP PC1-10 and z-scored within SOEP; the
# "cross" PGI the pooled-residualised, pooled-z PGI_Edu_z restricted to SOEP.
figs05_soep_fits <- function(datA) {
  dSO <- datA[as.character(datA$cohort) == "SOEP", , drop = FALSE]
  pcs <- WITHIN_PCS
  dSO$PGI_cross <- dSO$PGI_Edu_z
  ok <- stats::complete.cases(dSO[, c("PGI_Edu_raw", pcs)])
  r  <- rep(NA_real_, nrow(dSO))
  r[ok] <- stats::residuals(
    stats::lm(stats::as.formula(paste("PGI_Edu_raw ~", paste(pcs, collapse = "+"))),
              data = dSO[ok, , drop = FALSE]))
  dSO$PGI_Edu_z <- as.numeric(scale(r))
  cw <- function(m, term) {
    pt <- summary(m)$p.table
    e <- unname(pt[term, 1]); s <- unname(pt[term, 2])
    data.frame(term = term, estimate = e, se = s, statistic = unname(pt[term, 3]),
               p = unname(pt[term, 4]), ci_lo = e - Z_FIGS05 * s, ci_hi = e + Z_FIGS05 * s,
               stringsAsFactors = FALSE)
  }
  f_within <- "edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3"
  f_cross  <- "edu_z_kernel ~ (PGI_cross + BYc_z + east_west_c + gender_c)^3"
  m3 <- mgcv::bam(stats::as.formula(f_within), data = dSO, method = "REML")
  tw <- cw(m3, "PGI_Edu_z:BYc_z:east_west_c"); n_within <- n_used(f_within, dSO)
  rm(m3); invisible(gc())
  m3c <- mgcv::bam(stats::as.formula(f_cross), data = dSO, method = "REML")
  twc <- cw(m3c, "PGI_cross:BYc_z:east_west_c"); n_cross <- n_used(f_cross, dSO)
  rm(m3c); invisible(gc())
  dSO$reunif <- droplevels(factor(dSO$reunif))
  f_reunif <- "edu_z_kernel ~ (PGI_Edu_z + reunif + east_west_c + gender_c)^3"
  mr <- mgcv::bam(stats::as.formula(f_reunif), data = dSO, method = "REML")
  rt <- grep("^PGI_Edu_z:reunif.*:east_west_c$", rownames(summary(mr)$p.table), value = TRUE)[1]
  twr <- cw(mr, rt); n_reunif <- n_used(f_reunif, dSO)
  rm(mr); invisible(gc())
  # scale diagnostic: the two SOEP cells are not on one PGI scale by construction
  sd_cross_in_soep <- stats::sd(dSO$PGI_cross, na.rm = TRUE)
  n_soep <- nrow(dSO)
  rm(dSO); invisible(gc())
  list(within = tw, cross = twc, reunif = twr,
       n = c(within = n_within, cross = n_cross, reunif = n_reunif, soep_rows = n_soep),
       sd_cross_in_soep = sd_cross_in_soep, reunif_term = rt)
}

# ---------------------------------------------------------------------------
# Sample: rebuilt once from the merge through 00_setup/ (check_frozen_run.R step 4)
# ---------------------------------------------------------------------------
rebuild_attainment_sample <- function() {
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
  need <- c("education", "edu_z_kernel", "PGI_Edu_z", "BYc_z", "east_west_c", "gender_c",
            "east_west", "east", "gender", "cohort", "fid_re", "birth_year", "reunif",
            "IID", "PGI_Edu", "PGI_Edu_raw", WITHIN_PCS, CROSS_PCS)
  datA <- as.data.frame(dat_cluster)[, need]
  rm(list = intersect(c("dat_cluster", "dat_cluster_mob", "df_raw", "df", "df_with_twins",
                        "dat_edu", "dat_mob"), ls()))
  invisible(gc())
  list(datA = datA, data_file = DATA_FILE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_ancestry_sensitivity <- function(out_dir, results_dir, samples, run_ts, geno_final_dir) {
  run_ts <- sub("^RUN_", "", run_ts)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  commit <- source_commit(); dirty <- source_dirty()
  checks <- new_checks()
  t0 <- Sys.time()
  message(sprintf("A4 ancestry 2x2: RUN_%s\n  exports  -> %s\n  commit      %s%s",
                  run_ts, out_dir, substr(commit, 1, 8),
                  if (dirty) " (working tree dirty)" else ""))

  # ---- frozen references (read, never typed) ----
  prov <- utils::read.csv(file.path(results_dir, "_provenance.csv"), stringsAsFactors = FALSE)
  frozen_commit <- as.character(prov$analysis_repo_commit[1])
  A1c <- utils::read.csv(file.path(results_dir, "attainment", "main", "A1_Coefficients.csv"),
                         stringsAsFactors = FALSE)
  foc <- utils::read.csv(file.path(results_dir, "focal_estimates.csv"), stringsAsFactors = FALSE)
  s05 <- utils::read.csv(file.path(results_dir, "figS05_soep_threeway.csv"), stringsAsFactors = FALSE)

  i <- which(A1c$model == "A1_RE" & A1c$term == FOCAL)[1]
  if (is.na(i)) stop("frozen A1_RE ", FOCAL, " not found in A1_Coefficients.csv")
  ref_pooled <- data.frame(estimate = A1c$estimate[i], se = A1c$std.error[i],
                           statistic = A1c$statistic[i], p = A1c$p.value[i],
                           ci_lo = A1c$conf.low[i], ci_hi = A1c$conf.high[i],
                           stringsAsFactors = FALSE)
  # the same cell as focal_estimates.csv reports it: the two frozen files must agree
  j <- which(foc$analysis == "A1" & grepl("PGI x BY x Region", foc$label, fixed = TRUE))[1]
  if (is.na(j)) stop("frozen focal_estimates.csv row for A1 not found")
  checks$add("frozen: focal_estimates.csv agrees with A1_Coefficients.csv on the pooled cell",
             ref_pooled$estimate, foc$estimate[j], 1e-12,
             note = paste0("focal_estimates.csv row: ", foc$analysis[j], " / ", foc$label[j]))
  checks$add("frozen: focal_estimates.csv agrees on the pooled cell p",
             ref_pooled$p, foc$p[j], 1e-12)
  checks$add("frozen: focal_estimates.csv agrees on the pooled cell CI (lower)",
             ref_pooled$ci_lo, foc$ci_lo[j], 1e-12)

  s05_row <- function(contrast, construction_grepl) {
    k <- which(s05$contrast == contrast & grepl(construction_grepl, s05$construction))[1]
    if (is.na(k)) stop("frozen figS05 row not found: ", contrast, " / ", construction_grepl)
    data.frame(estimate = s05$b[k], ci_lo = s05$ci_lo[k], ci_hi = s05$ci_hi[k], p = s05$p[k],
               construction = s05$construction[k], stringsAsFactors = FALSE)
  }
  ref_soep_within <- s05_row("PGI x birth-year x region",       "^within-cohort")
  ref_soep_cross  <- s05_row("PGI x birth-year x region",       "^cross-cohort")
  ref_soep_reunif <- s05_row("PGI x reunification-step x region", "^within-cohort")
  # figS05 wrote no SE; invert its own CI rule (b +- 1.96 * SE) to recover it
  se_from_ci <- function(r) (r$ci_hi - r$ci_lo) / (2 * Z_FIGS05)

  # ---- sample ----
  datA <- samples$datA
  datA$east_west <- factor(as.character(datA$east_west), levels = c("West", "East"))
  datA$fid_re    <- factor(as.character(datA$fid_re))
  message(sprintf("  attainment sample: N = %d", nrow(datA)))
  t1 <- utils::read.csv(file.path(results_dir, "descriptives_table1_overall.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  checks$add("attainment sample N equals Table 1",
             t1$N_overall[t1$Measure == "Years of education"], nrow(datA), 0)
  if (checks$failed())
    stop("Sample size differs from the frozen Table 1 - wrong run or wrong sample source.")

  # ---- PC completeness: the two bases must cover the same rows, or the 2x2's
  #      adjustment axis would be confounded with a sample change ----
  full_fml <- paste("edu_z_kernel ~", A1_RHS, "+", RE_TERM)
  ok_cross  <- stats::complete.cases(datA[, c("PGI_Edu", CROSS_PCS)])
  ok_within <- stats::complete.cases(datA[, c("PGI_Edu", WITHIN_PCS)])
  checks$add("within-cohort PCs cover exactly the rows the cross-cohort PCs cover",
             sum(ok_cross), sum(ok_within), 0,
             note = "if this fails the adjustment axis is confounded with a sample change")
  checks$diag("rows with complete PGI_Edu and all 10 cross-cohort PCs", sum(ok_cross))
  checks$diag("rows with complete PGI_Edu and all 10 within-cohort PCs", sum(ok_within))

  # ---- genotyping batches ----
  gb <- genotyping_batch(datA, geno_final_dir)
  datA$.batch  <- gb$group
  datA$.cohort <- as.character(datA$cohort)
  checks$add("every SHIP row is assigned to a genotyping batch (SHIP-0 or SHIP-Td)",
             0, gb$info$n_unmatched, 0,
             note = sprintf("SHIP rows %d, matched %d, ambiguous %d",
                            gb$info$n_ship, gb$info$n_matched, gb$info$n_ambiguous))
  checks$add("no SHIP row matches both batch tables", 0, gb$info$n_ambiguous, 0)
  bt <- table(datA$.batch)
  for (b in names(bt)) checks$diag(paste0("group size, genotyping batch: ", b), as.numeric(bt[[b]]))

  # ---- cell 1: the published pooled / across-study cell, verified by refit ----
  message("verifying the published pooled cell (across-study PCs) ...")
  m_pub <- fit_bam("edu_z_kernel", A1_RHS, datA)
  r_pub <- p_row(m_pub, FOCAL)
  n_pub <- n_used(full_fml, datA)
  checks$add("frozen pooled cell reproduces on the rebuilt sample (estimate)",
             ref_pooled$estimate, r_pub$estimate, 1e-6)
  checks$add("frozen pooled cell reproduces on the rebuilt sample (SE)",
             ref_pooled$se, r_pub$se, 1e-6)
  rm(m_pub); invisible(gc())

  # ---- cell 2: THE REFIT. pooled sample, within-study PC residualisation ----
  message("fitting the missing cell: pooled sample, within-study PC residualisation ...")
  r_within_raw <- resid_pgi_within_groups(datA, "PGI_Edu", WITHIN_PCS, datA$.batch)

  # What each residualisation removes from the PGI, before the z-score puts both
  # back on SD 1. Recorded because per-study residualisation necessarily removes
  # the between-study mean differences as well as the ancestry axes, and the
  # Discussion sentence should be able to say how much of each.
  sd_pgi     <- stats::sd(datA$PGI_Edu, na.rm = TRUE)
  r_cross    <- resid_pgi_pooled(datA, "PGI_Edu", CROSS_PCS)
  grp_mean   <- stats::ave(datA$PGI_Edu, datA$.batch,
                           FUN = function(v) mean(v, na.rm = TRUE))
  r_demeaned <- datA$PGI_Edu - grp_mean
  var_share <- function(r) 100 * (1 - stats::var(r, na.rm = TRUE) / sd_pgi^2)
  checks$diag("SD of PGI_Edu in the analytic sample (before any residualisation)", sd_pgi)
  checks$diag("SD of the across-study (crossPC1-10) residual before pooled z-scoring",
              stats::sd(r_cross, na.rm = TRUE),
              note = "the published construction")
  checks$diag("SD of the within-study PC residual before pooled z-scoring",
              stats::sd(r_within_raw, na.rm = TRUE),
              note = "the refit's construction; z-scored to 1 in the fitted model")
  checks$diag("SD of PGI_Edu after study demeaning alone (no PCs)",
              stats::sd(r_demeaned, na.rm = TRUE),
              note = "isolates the part of the refit's extra shrinkage that is study demeaning")
  checks$diag("% of PGI variance removed by the across-study residualisation", var_share(r_cross))
  checks$diag("% of PGI variance removed by the within-study residualisation",
              var_share(r_within_raw))
  checks$diag("% of PGI variance removed by study demeaning alone", var_share(r_demeaned))
  rm(r_cross, r_demeaned, grp_mean); invisible(gc())

  dat2 <- datA
  dat2$PGI_Edu_z <- z_scale(r_within_raw)            # the ONLY change vs the published model
  n_fit <- n_used(full_fml, dat2)
  checks$add("the refit uses the same N as the published model",
             n_pub, n_fit, 0,
             note = "same rows, same formula; only the PGI's ancestry residualisation differs")
  m_fit <- fit_bam("edu_z_kernel", A1_RHS, dat2)
  r_fit <- p_row(m_fit, FOCAL)
  r_fit_main <- p_row(m_fit, "PGI_Edu_z")
  rm(m_fit); invisible(gc())
  checks$diag("PGI-Education main effect in the refit (published: read A1_Coefficients.csv)",
              r_fit_main$estimate,
              note = "for comparison with the published main effect; not a 2x2 cell")
  checks$diag("correlation of the within-study-residualised PGI with the canonical PGI_Edu_z",
              stats::cor(dat2$PGI_Edu_z, datA$PGI_Edu_z, use = "complete.obs"))
  for (b in sort(unique(datA$.batch))) {
    ii <- which(datA$.batch == b)
    checks$diag(paste0("correlation of the two PGIs within ", b),
                stats::cor(dat2$PGI_Edu_z[ii], datA$PGI_Edu_z[ii], use = "complete.obs"))
  }
  rm(dat2); invisible(gc())

  # ---- construction robustness for the refit (recorded, not a 2x2 cell) ----
  message("robustness of the refit: source variable and grouping ...")
  r_raw_src <- resid_pgi_within_groups(datA, "PGI_Edu_raw", WITHIN_PCS, datA$.batch)
  checks$add("refit PGI is invariant to starting from PGI_Edu_raw instead of PGI_Edu (max abs diff)",
             0, max(abs(z_scale(r_raw_src) - z_scale(r_within_raw)), na.rm = TRUE), 1e-9,
             note = paste("PGI_Edu is a single global affine transform of PGI_Edu_raw, so the",
                          "per-group residuals differ by one common factor that the pooled",
                          "z-score removes; figS05 residualises PGI_Edu_raw, the pipeline PGI_Edu"))
  rm(r_raw_src); invisible(gc())
  d4 <- datA
  d4$PGI_Edu_z <- z_scale(resid_pgi_within_groups(datA, "PGI_Edu", WITHIN_PCS, datA$.cohort))
  m4 <- fit_bam("edu_z_kernel", A1_RHS, d4)
  r_4grp <- p_row(m4, FOCAL)
  rm(m4, d4); invisible(gc())
  checks$add("refit is robust to grouping by cohort label (4 groups, SHIP pooled) rather than by genotyping batch (5)",
             r_fit$estimate, r_4grp$estimate, 0.01,
             note = sprintf("4-group estimate %.5f, p = %.4g; the 5-group fit is primary because SHIP-0 and SHIP-Td carry separate PCAs",
                            r_4grp$estimate, r_4grp$p))
  rm(r_within_raw); invisible(gc())

  # ---- cells 3-5: the frozen SOEP rows, verified by refit ----
  message("verifying the frozen SOEP cells (figS05) ...")
  s <- figs05_soep_fits(datA)
  checks$add("frozen SOEP / within-study cell reproduces (estimate)",
             ref_soep_within$estimate, s$within$estimate, 1e-6)
  checks$add("frozen SOEP / within-study cell reproduces (p)",
             ref_soep_within$p, s$within$p, 1e-6)
  checks$add("frozen SOEP / across-study cell reproduces (estimate)",
             ref_soep_cross$estimate, s$cross$estimate, 1e-6)
  checks$add("frozen SOEP / across-study cell reproduces (p)",
             ref_soep_cross$p, s$cross$p, 1e-6)
  checks$add("frozen SOEP / reunification-step row reproduces (estimate)",
             ref_soep_reunif$estimate, s$reunif$estimate, 1e-6)
  checks$add("frozen SOEP / reunification-step row reproduces (p)",
             ref_soep_reunif$p, s$reunif$p, 1e-6)
  checks$diag("SD within SOEP of the canonical pooled-z PGI (the across-study SOEP cell's PGI)",
              s$sd_cross_in_soep,
              note = paste("the within-study SOEP cell is z-scored WITHIN SOEP (SD 1) while the",
                           "across-study SOEP cell keeps the pooled z, so part of the .097-vs-.081",
                           "gap inside SOEP is PGI rescaling rather than adjustment; the pooled",
                           "row of the 2x2 has no such asymmetry"))

  # ---- the deliverable ----
  soep_model  <- "edu_z_kernel ~ (PGI x BYc_z x east_west_c x gender_c, all 3-way) - no family random effect (mgcv::bam, REML)"
  pooled_model <- paste0("A1_RE: edu_z_kernel ~ ", A1_RHS, " + ", RE_TERM,
                         " (mgcv::bam, fREML, discrete)")
  mk <- function(cell_id, cell, role, contrast, sample, pc_scope, pgi_std, term,
                 est, se, lo, hi, p, n, n_clusters, model, ci_basis, src, verified, note)
    data.frame(cell_id = cell_id, cell = cell, role = role, contrast = contrast,
               sample = sample, pc_scope = pc_scope, pgi_standardisation = pgi_std,
               term = term, estimate = est, se = se, ci_lo = lo, ci_hi = hi, p = p,
               n = n, n_clusters = n_clusters, model = model, ci_basis = ci_basis,
               source = src, verified = verified, note = note,
               source_commit = commit, run_ts = run_ts, stringsAsFactors = FALSE)

  used <- stats::complete.cases(datA[, intersect(all.vars(stats::as.formula(full_fml)),
                                                 names(datA)), drop = FALSE])
  n_clusters_pooled <- length(unique(as.character(datA$fid_re[used])))
  tab <- dplyr::bind_rows(
    mk(1L, "pooled sample x across-study PCs",
       "published primary specification",
       "PGI-Education x birth year x region",
       "pooled: BASE-II, SHIP (SHIP-0 + SHIP-Td), SOEP-G, TwinLife; attainment RE sample (dat_cluster)",
       "across-study (crossPC1-10, residualised once over the pooled sample)",
       "z-scored once over the pooled attainment analytic sample",
       FOCAL, ref_pooled$estimate, ref_pooled$se, ref_pooled$ci_lo, ref_pooled$ci_hi,
       ref_pooled$p, n_pub, n_clusters_pooled, pooled_model, "normal, qnorm(.975) x SE",
       "frozen: manuscript_export/attainment/main/A1_Coefficients.csv (model A1_RE); identical in focal_estimates.csv",
       "yes - refit on the rebuilt sample reproduces estimate and SE (ancestry_checks.csv)",
       "The published three-way. N and n_clusters are from the verification refit; the estimate, CI and p are read from the frozen file."),
    mk(2L, "pooled sample x within-study PCs",
       "the missing cell (fitted for this comment)",
       "PGI-Education x birth year x region",
       "pooled: BASE-II, SHIP (SHIP-0 + SHIP-Td), SOEP-G, TwinLife; attainment RE sample (dat_cluster)",
       "within-study (PC1-10, residualised separately within each genotyping batch: BASEII, SHIP0, SHIPTd, SOEP, TwinLife)",
       "z-scored once over the pooled attainment analytic sample",
       FOCAL, r_fit$estimate, r_fit$se, r_fit$ci_lo, r_fit$ci_hi, r_fit$p,
       n_fit, n_clusters_pooled, pooled_model, "normal, qnorm(.975) x SE",
       "fitted here (ancestry_sensitivity.R)",
       "n/a - this is the refit",
       paste("The published A1 specification with exactly one factor changed: the scope of the",
             "ancestry residualisation. Same rows, same outcome, same covariates, same estimator.",
             "Per-study residualisation also demeans the PGI by study (each group's intercept),",
             "which is a property of the construction the comment asks for.")),
    mk(3L, "SOEP-G only x across-study PCs",
       "2x2 cell",
       "PGI-Education x birth year x region",
       "SOEP-G only (the SOEP rows of the same attainment RE sample)",
       "across-study (crossPC1-10, canonical pooled pipeline)",
       "pooled z inherited from the pooled pipeline (NOT re-standardised within SOEP)",
       "PGI_cross:BYc_z:east_west_c", ref_soep_cross$estimate, se_from_ci(ref_soep_cross),
       ref_soep_cross$ci_lo, ref_soep_cross$ci_hi, ref_soep_cross$p,
       unname(s$n["cross"]), NA_integer_, soep_model, "normal, 1.96 x SE (as figS05 wrote it)",
       "frozen: manuscript_export/figS05_soep_threeway.csv, row 2",
       "yes - refit on the rebuilt sample reproduces estimate and p (ancestry_checks.csv)",
       paste("SE derived by inverting figS05's own CI rule (b +/- 1.96 x SE); figS05 wrote no SE column.",
             "The SOEP fits carry no family random effect and use method = 'REML', so the SOEP pair is",
             "internally comparable and the pooled pair is internally comparable, but the two pairs are",
             "not fit-identical.")),
    mk(4L, "SOEP-G only x within-study PCs",
       "2x2 cell",
       "PGI-Education x birth year x region",
       "SOEP-G only (the SOEP rows of the same attainment RE sample)",
       "within-study (SOEP PC1-10; the 2025-style reconstruction)",
       "z-scored within SOEP",
       "PGI_Edu_z:BYc_z:east_west_c", ref_soep_within$estimate, se_from_ci(ref_soep_within),
       ref_soep_within$ci_lo, ref_soep_within$ci_hi, ref_soep_within$p,
       unname(s$n["within"]), NA_integer_, soep_model, "normal, 1.96 x SE (as figS05 wrote it)",
       "frozen: manuscript_export/figS05_soep_threeway.csv, row 1",
       "yes - refit on the rebuilt sample reproduces estimate and p (ancestry_checks.csv)",
       paste("The cell the comment quotes at p = .048. Residualised from PGI_Edu_raw on the SOEP PC1-10",
             "and z-scored within SOEP; because cell 3 keeps the pooled z instead, part of the difference",
             "between cells 3 and 4 is PGI rescaling and not adjustment. SE derived by inverting figS05's",
             "CI rule. No family random effect.")),
    mk(5L, "SOEP-G only x within-study PCs, reunification-step contrast",
       "the other contrast (not a cell of the 2x2)",
       "PGI-Education x reunification step (born >= 1975) x region",
       "SOEP-G only (the SOEP rows of the same attainment RE sample)",
       "within-study (SOEP PC1-10; the 2025-style reconstruction)",
       "z-scored within SOEP",
       s$reunif_term, ref_soep_reunif$estimate, se_from_ci(ref_soep_reunif),
       ref_soep_reunif$ci_lo, ref_soep_reunif$ci_hi, ref_soep_reunif$p,
       unname(s$n["reunif"]), NA_integer_,
       "edu_z_kernel ~ (PGI x reunif x east_west_c x gender_c, all 3-way) - no family random effect (mgcv::bam, REML)",
       "normal, 1.96 x SE (as figS05 wrote it)",
       "frozen: manuscript_export/figS05_soep_threeway.csv, row 3",
       "yes - refit on the rebuilt sample reproduces estimate and p (ancestry_checks.csv)",
       paste("Carried alongside the 2x2 because it is the same SOEP-G within-PC construction asking a",
             "different question: a pre/post reunification step rather than a continuous birth-year",
             "gradient. It is not a cell of the 2x2 and must not be read as one.")))

  # ---- write ----
  w <- function(df, name) {
    utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
    message(sprintf("  [ok]   %-34s %d x %d", name, nrow(df), ncol(df)))
  }
  w(tab, "ancestry_sensitivity.csv")
  chk <- checks$table(); chk$source_commit <- commit; chk$run_ts <- run_ts
  w(chk, "ancestry_checks.csv")
  provo <- data.frame(
    run_ts = run_ts, frozen_analysis_repo_commit = frozen_commit, script_commit = commit,
    script_working_tree_dirty = dirty, reference_source = results_dir,
    data_source = samples$data_file, genotype_source = geno_final_dir,
    n_attainment = nrow(datA), n_fitted_cell = n_fit,
    cells_read_frozen = 4L, cells_fitted_here = 1L,
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R_version = R.version.string,
    mgcv_version = as.character(utils::packageVersion("mgcv")),
    checks_total = sum(chk$kind == "check"), checks_passed = sum(chk$kind == "check" & chk$pass),
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    stringsAsFactors = FALSE)
  w(provo, "ancestry_provenance.csv")

  # ---- console summary (aggregates only) ----
  cat("\n== checks ==\n")
  for (i in which(chk$kind == "check"))
    cat(sprintf("  [%s] %-92s diff %.3g (tol %.3g)\n", if (chk$pass[i]) "ok" else "FAIL",
                substr(chk$check[i], 1, 92), chk$abs_diff[i], chk$tol[i]))
  cat("\n== diagnostics ==\n")
  for (i in which(chk$kind == "diagnostic"))
    cat(sprintf("  %-92s %.5g\n", substr(chk$check[i], 1, 92), chk$observed[i]))
  cat("\n== the 2x2 ==\n")
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %d %-52s b = %+.4f  95%% CI [%+.4f, %+.4f]  p = %.4g  n = %5d  %s\n",
                tab$cell_id[i], substr(tab$cell[i], 1, 52), tab$estimate[i],
                tab$ci_lo[i], tab$ci_hi[i], tab$p[i], tab$n[i],
                if (grepl("^fitted", tab$source[i])) "FITTED HERE" else "frozen"))
  if (checks$failed())
    cat("\nWARNING: some checks failed - read ancestry_checks.csv before using the exports.\n")
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
GENO_FINAL_DIR <- file.path(Sys.getenv("DATA_ROOT"), "01_Genotype", "Harmonized", "data", "final")
MS <- Sys.getenv("MANUSCRIPT_REPO", unset = file.path(
  path.expand("~"), "Documents", "017_manuscript", "80_years_GxE_on_Education_in_Germany-manuscript"))
RESULTS_DIR <- file.path(MS, "01_results", "manuscript_export")
if (!file.exists(file.path(RESULTS_DIR, "_provenance.csv"))) {
  message("01_results/ not found at ", MS, " - reading the frozen exports from the run folder")
  RESULTS_DIR <- file.path(RUN_DIR, "manuscript_export")
}
message("Rebuilding the attainment analytic sample from the merge (00_setup/) ...")
SAMPLES <- rebuild_attainment_sample(); invisible(gc())
invisible(run_ancestry_sensitivity(OUT_DIR, RESULTS_DIR, SAMPLES, RUN_TS, GENO_FINAL_DIR))
