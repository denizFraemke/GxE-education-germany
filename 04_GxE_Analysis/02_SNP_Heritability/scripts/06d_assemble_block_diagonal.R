#!/usr/bin/env Rscript
########################################################################
## 06d_assemble_block_diagonal.R — Assemble per-stratum block-diagonal h²
##
## Reads the .hsq files written by step 05d and writes:
##   output/final/h2_<TRACK>.tsv               — one row per (stratum, model)
##   output/final/h2_<TRACK>_heterogeneity.tsv — Q + pairwise contrasts
##   output/final/h2_<TRACK>_summary.md        — human-readable report
##
## Parameterised by env vars (set by run_pipeline.sh or the wrappers in
## scripts/run_blockdiag_*.sh). Defaults reproduce the validated cell-level
## primary track exactly (TRACK="blockdiag", STRATA = STRATA_MAIN).
##
##   STRATA_LIST_NAMES         space-separated stratum labels
##   TRACK                     short label used in output filenames
##   BLOCKDIAG_OUT_SUBDIR      subdir of output/reml/ containing .hsq files
##   BLOCKDIAG_COVAR_FILE_OVERRIDE  covar file actually used by 05d
##
## Aggregate output only (counts, h², SE, CI, group means/SDs, fractions) —
## no individual-level data values are read or written.
########################################################################
rm(list = ls())

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))
source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

# ---- Track parameterisation -------------------------------------------
TRACK     <- Sys.getenv("TRACK", unset = "blockdiag")
SUBDIR    <- Sys.getenv("BLOCKDIAG_OUT_SUBDIR", unset = "blockdiag")
STRATA_E  <- Sys.getenv("STRATA_LIST_NAMES", unset = "")
COVAR_USED<- Sys.getenv("BLOCKDIAG_COVAR_FILE_OVERRIDE",
                        unset = file.path(PROJ_ROOT, "data", "gender_cohort.covar"))

STRATA <- if (nzchar(STRATA_E)) {
  strsplit(trimws(STRATA_E), "\\s+")[[1]]
} else {
  c("East_pre1990", "East_post1990", "West_pre1990", "West_post1990")
}

BLOCKDIAG_DIR <- file.path(PROJ_ROOT, "output", "reml", SUBDIR)
GRM_DIR       <- file.path(PROJ_ROOT, "output", "grm")
FINAL_OUT     <- file.path(PROJ_ROOT, "output", "final")
PHENO_FILE    <- file.path(PROJ_ROOT, "data", "education.phen")
QCOVAR_FILE   <- file.path(PROJ_ROOT, "data", "crosspcs.qcovar")
KEEP_DIR      <- file.path(PROJ_ROOT, "data", "keep")
dir.create(FINAL_OUT, recursive = TRUE, showWarnings = FALSE)

OUT_TSV    <- file.path(FINAL_OUT, sprintf("h2_%s.tsv", TRACK))
OUT_HET    <- file.path(FINAL_OUT, sprintf("h2_%s_heterogeneity.tsv", TRACK))
OUT_REPORT <- file.path(FINAL_OUT, sprintf("h2_%s_summary.md", TRACK))

cat("=== Step 06d: Assemble per-stratum block-diagonal h² outputs ===\n")
cat(sprintf("  Track:       %s\n", TRACK))
cat(sprintf("  REML dir:    %s\n", BLOCKDIAG_DIR))
cat(sprintf("  Strata (%d): %s\n", length(STRATA), paste(STRATA, collapse = ", ")))
cat(sprintf("  Covar file:  %s\n", COVAR_USED))

if (!dir.exists(BLOCKDIAG_DIR)) {
  cat("\nNo block-diagonal REML directory — run step 05d first.\n")
  quit(status = 1)
}

# ---- Cohort label inference for block-component file lookup ----------
COHORTS_REML <- c("BASEII", "SHIP", "SOEP", "TWINLIFE")

# ---- Parse one stratum's .hsq into a tidy row ------------------------
parse_stratum <- function(stratum,
                          model = c("standalone", "standalone_noconstrain")) {
  model <- match.arg(model)
  label <- if (model == "standalone") {
    paste0("blockdiag_", stratum)
  } else {
    paste0("blockdiag_", stratum, "_noconstrain")
  }
  f <- file.path(BLOCKDIAG_DIR, paste0(label, ".hsq"))
  hsq <- parse_hsq(f)
  if (!isTRUE(hsq$.present)) {
    cat(sprintf("  WARN: %s.hsq missing or unparseable\n", label))
    return(NULL)
  }
  h2    <- hsq[["V(G)/Vp"]];   h2_se <- hsq[["V(G)/Vp_SE"]]
  v_g   <- hsq[["V(G)"]];      v_g_se<- hsq[["V(G)_SE"]]
  v_e   <- hsq[["V(e)"]];      v_e_se<- hsq[["V(e)_SE"]]
  v_p   <- hsq[["Vp"]]
  n     <- hsq[["n"]]
  if (is.null(h2) || is.null(h2_se) || is.null(n)) {
    cat(sprintf("  WARN: %s missing required fields\n", label))
    return(NULL)
  }
  data.frame(
    label   = label, stratum = stratum, model = model,
    N       = as.integer(n),
    V_G     = v_g, V_G_SE = v_g_se,
    V_e     = v_e, V_e_SE = v_e_se,
    Vp      = v_p,
    h2_SNP  = h2,  h2_SNP_SE = h2_se,
    ci_lower = h2 - 1.96 * h2_se,
    ci_upper = h2 + 1.96 * h2_se,
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (s in STRATA) {
  for (m in c("standalone", "standalone_noconstrain")) {
    r <- parse_stratum(s, m)
    if (!is.null(r)) rows[[length(rows) + 1L]] <- r
  }
}
if (length(rows) == 0L) stop("No parseable .hsq rows found.")
h2_tab <- do.call(rbind, rows)
h2_tab$stratum <- factor(h2_tab$stratum, levels = STRATA)
h2_tab <- h2_tab[order(h2_tab$stratum, h2_tab$model), ]
h2_tab$stratum <- as.character(h2_tab$stratum)

write.table(h2_tab, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote: %s (%d rows)\n", OUT_TSV, nrow(h2_tab)))

# ---- Heterogeneity: Cochran's Q + pairwise Z on constrained strata ---
main <- h2_tab[h2_tab$model == "standalone", ]
Q <- NA_real_; df <- NA_integer_; pQ <- NA_real_; I2 <- NA_real_
h2_pool <- NA_real_; pool_se <- NA_real_

if (nrow(main) >= 2L) {
  w  <- 1 / main$h2_SNP_SE^2
  h2_pool <- sum(w * main$h2_SNP) / sum(w)
  Q  <- sum(w * (main$h2_SNP - h2_pool)^2)
  df <- nrow(main) - 1L
  pQ <- pchisq(Q, df = df, lower.tail = FALSE)
  I2 <- max(0, (Q - df) / Q)
  pool_se <- 1 / sqrt(sum(w))

  cat(sprintf("\nOmnibus Cochran's Q (%s, constrained):\n", TRACK))
  cat(sprintf("  Q = %.3f, df = %d, p = %.3f, I² = %.3f\n", Q, df, pQ, I2))
  cat(sprintf("  Inverse-variance pooled h² = %.3f (pooled SE = %.3f, 95%% CI [%.3f, %.3f])\n",
              h2_pool, pool_se, h2_pool - 1.96 * pool_se, h2_pool + 1.96 * pool_se))

  het_rows <- list(data.frame(
    test = "cochrans_Q", analysis = TRACK, K = nrow(main), df = df,
    statistic = Q, p_value = pQ, I2 = I2, h2_pooled = h2_pool,
    pooled_SE = pool_se, cells = NA_character_,
    note = "Q on constrained stratum fits", stringsAsFactors = FALSE
  ))

  for (i in seq_len(nrow(main) - 1L)) {
    for (j in (i + 1L):nrow(main)) {
      z <- (main$h2_SNP[i] - main$h2_SNP[j]) /
           sqrt(main$h2_SNP_SE[i]^2 + main$h2_SNP_SE[j]^2)
      pz <- 2 * pnorm(-abs(z))
      het_rows[[length(het_rows) + 1L]] <- data.frame(
        test = "pairwise_Z", analysis = TRACK,
        K = NA_integer_, df = NA_integer_,
        statistic = z, p_value = pz, I2 = NA_real_, h2_pooled = NA_real_,
        pooled_SE = NA_real_,
        cells = paste(main$stratum[i], main$stratum[j], sep = ","),
        note = sprintf("h2[%s]=%.3f vs h2[%s]=%.3f",
                       main$stratum[i], main$h2_SNP[i],
                       main$stratum[j], main$h2_SNP[j]),
        stringsAsFactors = FALSE
      )
    }
  }

  het <- do.call(rbind, het_rows)
  write.table(het, OUT_HET, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Wrote: %s (%d rows)\n", OUT_HET, nrow(het)))
}

# ---- Diagnostics: cohort composition + phenotype/birth-year summaries -
# All outputs are aggregate (counts, means, SDs, fractions) — no individual
# rows are printed.

# Cohort composition per stratum, from per-(cohort × stratum) block-component
# .grm.id files written by step 03d / 03e.
cohort_composition <- function(stratum) {
  out <- data.frame(
    stratum = stratum, cohort = COHORTS_REML, N = 0L,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(COHORTS_REML)) {
    f <- file.path(GRM_DIR,
                   sprintf("grm_blockcomp_%s_%s.grm.id",
                           COHORTS_REML[k], stratum))
    if (file.exists(f)) {
      out$N[k] <- as.integer(length(readLines(f)))
    }
  }
  out$share <- out$N / max(1L, sum(out$N))
  out
}

comp_rows <- lapply(STRATA, cohort_composition)
comp_tab  <- do.call(rbind, comp_rows)

# Phenotype + birth-year aggregates per stratum.
#
# Matching is on IID only, not (FID, IID). The keep files retain source-data
# FIDs from export_stratum_phenotype.R, while the phen/qcovar files have been
# rewritten in place by 05d's align_input_fids() to use pooled.fam canonical
# FIDs. So (FID, IID) matching silently undercounts BASE-II and TwinLife
# rows in these diagnostic tables. IIDs are unique across the analytic
# sample (the export step guarantees this), so IID-only matching is correct
# for joining people across files.
read_keep_iids <- function(stratum) {
  f <- file.path(KEEP_DIR, paste0(stratum, ".iids"))
  if (!file.exists(f)) return(character(0))
  k <- read.table(f, header = FALSE, stringsAsFactors = FALSE,
                  colClasses = c("character", "character"))
  k$V2
}

read_pheno <- function() {
  if (!file.exists(PHENO_FILE)) return(NULL)
  d <- read.table(PHENO_FILE, header = FALSE, stringsAsFactors = FALSE,
                  colClasses = c("character", "character", "numeric"))
  names(d) <- c("FID", "IID", "y")
  d
}

read_qcovar_birthyear <- function() {
  if (!file.exists(QCOVAR_FILE)) return(NULL)
  hdr <- tryCatch(readLines(QCOVAR_FILE, n = 1L), error = function(e) NA_character_)
  if (is.na(hdr) || !nzchar(hdr)) return(NULL)
  # qcovar layout (export_stratum_phenotype.R): FID IID birth_year_c
  # birth_year_c_sq <PC1..PC20>. Read all columns and pick by position.
  d <- read.table(QCOVAR_FILE, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 3L) return(NULL)
  names(d)[1:3] <- c("FID", "IID", "birth_year_c")
  d[, c("IID", "birth_year_c")]
}

pheno <- read_pheno()
qcov  <- read_qcovar_birthyear()

pheno_summary <- data.frame(
  stratum = STRATA,
  N_pheno = NA_integer_,
  edu_mean = NA_real_, edu_sd = NA_real_, edu_var = NA_real_,
  by_centered_mean = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(STRATA)) {
  iids <- read_keep_iids(STRATA[i])
  if (length(iids) == 0L) next
  if (!is.null(pheno)) {
    y <- pheno$y[match(iids, pheno$IID)]
    y <- y[!is.na(y)]
    if (length(y) > 0L) {
      pheno_summary$N_pheno[i]  <- length(y)
      pheno_summary$edu_mean[i] <- mean(y)
      pheno_summary$edu_sd[i]   <- sd(y)
      pheno_summary$edu_var[i]  <- var(y)
    }
  }
  if (!is.null(qcov)) {
    by_c <- qcov$birth_year_c[match(iids, qcov$IID)]
    by_c <- by_c[!is.na(by_c)]
    if (length(by_c) > 0L) {
      pheno_summary$by_centered_mean[i] <- mean(by_c)
    }
  }
}

# ---- Human-readable report ----------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f"), x))
fmt_n <- function(n) formatC(n, big.mark = ",", format = "d")

covar_short <- basename(COVAR_USED)
covar_desc <- switch(
  covar_short,
  "gender_cohort.covar" = "gender_int + cohort_int (categorical)",
  "region_cohort.covar" = "region_int + cohort_int (categorical)",
  "cohort_only.covar"   = "cohort_int (categorical)",
  covar_short
)

md <- c(
  sprintf("# Block-diagonal stratum-level h²_SNP — track `%s` (step 06d)", TRACK),
  "",
  sprintf("_Generated %s by `06d_assemble_block_diagonal.R`._",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  sprintf(paste0(
    "**Specification.** Per-stratum GREML on block-diagonal GRMs built ",
    "from per-(cohort × stratum) sub-GRMs (within-cohort allele-frequency ",
    "standardization; cross-cohort entries set to zero). V_G is identified ",
    "entirely from within-cohort relatedness signal. Covariates: 20 cross-",
    "cohort PCs + birth_year_c + birth_year_c_sq in `--qcovar`; %s in ",
    "`--covar` (`%s`). The cohort indicator absorbs cohort-level mean ",
    "differences in years of education that would otherwise load into V_e ",
    "(see `Plan_deviations.md` §7g for the diagnostic that established ",
    "this; §7h documents the secondary tracks built on the same design)."),
    covar_desc, covar_short),
  "",
  "## Stratum-level h²",
  "",
  "| Stratum | N | h²_SNP | REML SE | 95% Wald CI | flag |",
  "|---|---:|---:|---:|:---:|:---|"
)
boundary_flagged <- character(0)
for (i in seq_len(nrow(main))) {
  r <- main[i, ]
  flag <- character(0)
  if (!is.na(r$h2_SNP) && r$h2_SNP <= 0.005) flag <- c(flag, "h²≈0 boundary")
  if (!is.na(r$h2_SNP) && r$h2_SNP >= 0.995) flag <- c(flag, "h²≈1 boundary")
  if (!is.na(r$h2_SNP_SE) && r$h2_SNP_SE >= 0.3) flag <- c(flag, "imprecise")
  if (length(flag) > 0L) boundary_flagged <- c(boundary_flagged, r$stratum)
  md <- c(md, sprintf("| %s | %s | **%s** | %s | [%s, %s] | %s |",
                      r$stratum, fmt_n(r$N), fmt(r$h2_SNP, 3),
                      fmt(r$h2_SNP_SE, 3), fmt(r$ci_lower, 3), fmt(r$ci_upper, 3),
                      paste(flag, collapse = "; ")))
}

# Variance components per stratum (ANALYSIS_PLAN §S5 line 254: V_G = h²_SNP × V_P
# reported per stratum so that differences in h² can be decomposed into
# differences in V_G vs differences in V_P). No formal Q test on V_G across
# strata is reported: V_G SEs are large in the imprecise cells (e.g. ≥ 2.6
# in the post-1990 primary cells), so a Q test would be uninformative. See
# Plan_deviations.md §7i for the rationale.
md <- c(md, "", "## Variance components per stratum",
        "",
        "Per the analysis plan (§S5): V_G is the additive SNP variance ",
        "(years² of education); V_P is the residual phenotypic variance on ",
        "the observed scale; h² = V_G / V_P.",
        "")
md <- c(md, "| Stratum | N | V_G | V_G SE | V_e | V_e SE | V_P | h² = V_G / V_P |",
        "|---|---:|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(main))) {
  r <- main[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                      r$stratum, fmt_n(r$N),
                      fmt(r$V_G, 3), fmt(r$V_G_SE, 3),
                      fmt(r$V_e, 3), fmt(r$V_e_SE, 3),
                      fmt(r$Vp, 3),
                      fmt(r$h2_SNP, 3)))
}

if (nrow(main) >= 2L) {
  md <- c(md, "",
          sprintf("**Cross-stratum heterogeneity (Cochran's Q):** Q = %.3f, df = %d, p = %.3f, I² = %.3f.",
                  Q, df, pQ, I2),
          sprintf("**Inverse-variance pooled h² (fixed-effects):** %.3f (pooled SE = %.3f, 95%% Wald CI [%.3f, %.3f]).",
                  h2_pool, pool_se, h2_pool - 1.96 * pool_se, h2_pool + 1.96 * pool_se),
          "",
          "Non-significant Q is consistent with **no evidence of heterogeneity** across strata; with the precision available in these data the Q test has limited power, so non-significance does not by itself establish that h² is the same across strata.")

  md <- c(md, "", "## Pairwise contrasts (Z = (h²_i − h²_j) / √(SE_i² + SE_j²))",
          "", "| Comparison | Z | p (uncorrected) | h² estimates |",
          "|---|---:|---:|---|")
  for (i in seq_len(nrow(main) - 1L)) {
    for (j in (i + 1L):nrow(main)) {
      z <- (main$h2_SNP[i] - main$h2_SNP[j]) /
           sqrt(main$h2_SNP_SE[i]^2 + main$h2_SNP_SE[j]^2)
      pz <- 2 * pnorm(-abs(z))
      md <- c(md, sprintf("| %s vs %s | %s | %s | %s vs %s |",
                          main$stratum[i], main$stratum[j],
                          fmt(z, 3), fmt(pz, 3),
                          fmt(main$h2_SNP[i], 3), fmt(main$h2_SNP[j], 3)))
    }
  }
}

# Cohort composition diagnostic.
md <- c(md, "", "## Cohort composition per stratum",
        "",
        "Counts are the number of individuals contributing to each ",
        "block component (cohort × stratum) of the block-diagonal GRM.",
        "")
md <- c(md, "| Stratum | BASE-II | SHIP | SOEP | TwinLife | Total | max share |",
        "|---|---:|---:|---:|---:|---:|---:|")
single_cohort_warn <- character(0)
for (s in STRATA) {
  cs <- comp_tab[comp_tab$stratum == s, ]
  totalN <- sum(cs$N)
  max_share <- if (totalN > 0) max(cs$share) else NA_real_
  if (!is.na(max_share) && max_share > 0.80) {
    single_cohort_warn <- c(single_cohort_warn,
                            sprintf("%s (%s, %s%%)", s,
                                    cs$cohort[which.max(cs$share)],
                                    fmt(100 * max_share, 1)))
  }
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                      s,
                      fmt_n(cs$N[cs$cohort == "BASEII"]),
                      fmt_n(cs$N[cs$cohort == "SHIP"]),
                      fmt_n(cs$N[cs$cohort == "SOEP"]),
                      fmt_n(cs$N[cs$cohort == "TWINLIFE"]),
                      fmt_n(totalN),
                      fmt(max_share, 3)))
}

if (length(single_cohort_warn) > 0L) {
  md <- c(md, "",
          sprintf("> **Single-cohort dominance:** one cohort contributes > 80%% of: %s. Treat these stratum estimates as cohort-weighted rather than population-representative.",
                  paste(single_cohort_warn, collapse = "; ")))
}

# Phenotype + birth-year aggregates.
md <- c(md, "", "## Phenotype + birth-year summary per stratum",
        "",
        "Aggregate group statistics on the analysis sample actually fit ",
        "by GREML. Education is years; birth_year_c is mean-centered ",
        "(value 0 = overall mean birth year).",
        "")
md <- c(md, "| Stratum | N (pheno) | edu mean | edu SD | edu Var | mean birth_year_c |",
        "|---|---:|---:|---:|---:|---:|")
for (i in seq_along(STRATA)) {
  ps <- pheno_summary[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |",
                      ps$stratum, fmt_n(ps$N_pheno),
                      fmt(ps$edu_mean, 3), fmt(ps$edu_sd, 3), fmt(ps$edu_var, 3),
                      fmt(ps$by_centered_mean, 3)))
}

# Track-specific interpretation: the default cell-level primary gets the
# numbers-anchored text, other tracks get conservative comparison-track
# language.
md <- c(md, "",
        "## Interpretation",
        "",
        sprintf("Track: **%s**. Covariates: %s (file `%s`). %d strata: %s.",
                TRACK, covar_desc, covar_short,
                length(STRATA), paste(STRATA, collapse = ", ")),
        "",
        "These stratum-level h² estimates use **within-cohort genomic relatedness only** (cross-cohort GRM entries set to zero), with cohort fixed effects in the model to absorb cohort-level mean differences in years of education.")

if (TRACK == "blockdiag") {
  md <- c(md, "",
          "Estimates converged across specifications in the well-powered cells (REML SE ≤ 0.10): East × pre-1990 and West × pre-1990 both around h² ≈ 0.28. The post-1990 cells remained imprecise (REML SE ~0.3–0.5; 95% CIs cover most of the [0, 1] range) and do not support meaningful cell-specific inference under any of the GRM designs run in this section.",
          "",
          "Across the three specifications run on this sample (pooled-GRM at `--grm-cutoff 0.20`, per-cohort, and block-diagonal with cohort fixed effects), there was no evidence of East/West or pre/post-Reunification heterogeneity in the well-powered cells. Power to detect heterogeneity involving the post-1990 cells was limited.",
          "",
          "Side-by-side comparison with the pooled-GRM track (cutoff 0.20, treated as sensitivity) and the per-cohort track is in the section README.")
} else {
  md <- c(md, "",
          "This is a **secondary / exploratory** comparison track, run on the same block-diagonal design as the cell-level primary. It uses the *subset-then-prune* GRM order (subset the pre-cutoff per-cohort GRM to the new stratum, then `--grm-cutoff 0.05` within each cohort × stratum block) so that within-cohort relatives split across strata are retained as block-diagonal information.",
          "",
          "Interpret the per-stratum h² and the heterogeneity test in light of the diagnostics above: if one cohort dominates a stratum (>80% share) the estimate is cohort-weighted rather than population-representative; the Cochran's Q test has limited power with REML SEs in this range.",
          "",
          "Sex labels are **recorded gender/sex** as harmonised across cohorts. Treat results as exploratory and conservative.")
}

writeLines(md, OUT_REPORT)
cat(sprintf("Wrote: %s\n", OUT_REPORT))

cat("\n=== Step 06d complete ===\n")
