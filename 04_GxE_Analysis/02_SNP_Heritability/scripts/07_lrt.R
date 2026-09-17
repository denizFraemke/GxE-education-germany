#!/usr/bin/env Rscript
########################################################################
## 07_lrt.R — Heterogeneity tests on per-stratum h²_SNP (Tardis)
##
## The analysis plan (Plan_deviations.md §7b) specified GCTA's
## `--reml --mgrm` joint multi-component framework to compute LRTs of
## "all sigma_G equal" vs the full free model. That approach requires
## either overlapping individuals across GRMs (GCTA's documented use
## case) or a single GRM combining all strata. Neither holds for this
## stratum design: per-stratum GRMs have **disjoint** individuals, and
## GCTA 1.95.1 segfaults (exit code 11) on every joint fit when
## --reml --mgrm is passed disjoint-sample lists.
##
## This script therefore computes heterogeneity tests directly from
## the per-stratum standalone REML fits:
##
##   - Omnibus Cochran's Q on the K cell h² point estimates (the
##     analogue of the plan's H0: all h²_SNP equal across strata).
##     Q = Σ_k (h²_k − h²_pooled)² / SE(h²_k)²,
##     h²_pooled = Σ_k w_k h²_k / Σ_k w_k,  w_k = 1 / SE(h²_k)²,
##     Q ~ χ²_{K-1} under H0.
##
##   - Pairwise Z-tests between each (i, j) pair of cells:
##     Z_ij = (h²_i − h²_j) / sqrt(SE_i² + SE_j²),
##     P_ij = 2 · (1 − Φ(|Z_ij|)).
##
##   - I² heterogeneity statistic (Higgins & Thompson 2002):
##     I² = max(0, (Q − (K-1)) / Q) — fraction of total variation due
##     to between-cell heterogeneity rather than within-cell sampling.
##
## Caveats — versus the joint LRT the plan specified:
##   - Cochran's Q on inverse-variance-weighted ratio estimates is an
##     asymptotic χ² test; for small per-cell SEs and modest K it
##     behaves like the LRT but loses power when SEs are large
##     relative to between-cell differences.
##   - This test treats each cell's h² as an independent estimate
##     with its own σ²_e, whereas the plan's joint mgrm framework
##     would have constrained σ²_e to be shared across cells
##     (Plan_deviations §7b). The null tested here is therefore h²
##     equality under per-cell residual variance — a slightly
##     different null than the plan's, but closer to the question
##     being asked ("is heritability equal across these contexts?"),
##     since fixing σ²_e a priori is an additional modelling
##     assumption.
##
## Output: ${REML_OUT}/heterogeneity.tsv, plus a compatibility symlink
## lrt.tsv -> heterogeneity.tsv for step 09 and other consumers that
## read the lrt.tsv name.
########################################################################

rm(list = ls())

# ---- Resolve project root via this script's path --------------------
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

REML_OUT  <- file.path(PROJ_ROOT, "output", "reml")
MAIN_DIR  <- file.path(REML_OUT, "main")
OUT_FILE  <- file.path(REML_OUT, "heterogeneity.tsv")
LRT_LINK  <- file.path(REML_OUT, "lrt.tsv")  # compatibility symlink

cat("=== Step 07: Heterogeneity tests (Cochran's Q + pairwise Z) ===\n")

# ---- Pull h² + SE per cell from each standalone .hsq -----------------
extract_h2 <- function(hsq_path) {
  h <- parse_hsq(hsq_path)
  if (!isTRUE(h$.present)) return(NULL)
  # parse_hsq stores keys exactly as they appear in the .hsq file.
  # For a single-component standalone fit those are:
  #   V(G)        -> point + SE
  #   V(e)        -> point + SE
  #   Vp          -> point + SE
  #   V(G)/Vp     -> point + SE  (this IS h² and its SE)
  list(
    h2     = h[["V(G)/Vp"]],
    h2_SE  = h[["V(G)/Vp_SE"]],
    n      = h[["n"]],
    V_G    = h[["V(G)"]],
    V_G_SE = h[["V(G)_SE"]],
    V_P    = h[["Vp"]],
    logL   = h[["logL"]]
  )
}

# ---- Run the heterogeneity tests for one analysis --------------------
heterogeneity_block <- function(analysis_label, dir, strata, model_suffix) {
  rows <- list()
  for (s in strata) {
    f <- file.path(dir, sprintf("standalone_%s%s.hsq", s, model_suffix))
    r <- extract_h2(f)
    if (is.null(r)) {
      cat(sprintf("  %s%s: missing .hsq — skipping\n", s, model_suffix))
      next
    }
    rows[[s]] <- c(stratum = s, r)
  }
  if (length(rows) < 2L) {
    cat(sprintf("  %s: fewer than 2 fits available — no test computed.\n",
                analysis_label))
    return(NULL)
  }

  h2_vec  <- sapply(rows, function(r) r$h2)
  se_vec  <- sapply(rows, function(r) r$h2_SE)
  ok      <- !is.na(h2_vec) & !is.na(se_vec) & is.finite(se_vec) & se_vec > 0
  if (sum(ok) < 2L) {
    cat(sprintf("  %s: <2 cells with usable SE — skipping.\n", analysis_label))
    return(NULL)
  }
  h2_vec <- h2_vec[ok]
  se_vec <- se_vec[ok]
  strata_used <- names(rows)[ok]
  K <- length(h2_vec)

  # Omnibus Cochran's Q
  w        <- 1 / se_vec^2
  h2_pool  <- sum(w * h2_vec) / sum(w)
  Q        <- sum(w * (h2_vec - h2_pool)^2)
  df       <- K - 1
  pQ       <- pchisq(Q, df = df, lower.tail = FALSE)
  I2       <- max(0, (Q - df) / Q)

  out <- list()

  out[[length(out) + 1L]] <- data.frame(
    analysis      = analysis_label,
    model_suffix  = ifelse(nzchar(model_suffix), sub("^_", "", model_suffix), "constrained"),
    test          = "cochrans_Q",
    cells         = paste(strata_used, collapse = ","),
    K             = K,
    df            = df,
    statistic     = Q,
    p_value       = pQ,
    I2            = I2,
    h2_pooled     = h2_pool,
    note          = "omnibus heterogeneity test of per-cell h^2",
    stringsAsFactors = FALSE
  )

  # Pairwise Z-tests
  for (i in seq_len(K - 1L)) {
    for (j in (i + 1L):K) {
      d  <- h2_vec[i] - h2_vec[j]
      se <- sqrt(se_vec[i]^2 + se_vec[j]^2)
      z  <- d / se
      p  <- 2 * pnorm(-abs(z))
      out[[length(out) + 1L]] <- data.frame(
        analysis      = analysis_label,
        model_suffix  = ifelse(nzchar(model_suffix), sub("^_", "", model_suffix), "constrained"),
        test          = sprintf("pairwise_Z[%s_vs_%s]", strata_used[i], strata_used[j]),
        cells         = paste(strata_used[c(i, j)], collapse = ","),
        K             = 2L,
        df            = 1L,
        statistic     = z,
        p_value       = p,
        I2            = NA_real_,
        h2_pooled     = NA_real_,
        note          = sprintf("h2[%s]=%.3f vs h2[%s]=%.3f",
                                strata_used[i], h2_vec[i],
                                strata_used[j], h2_vec[j]),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out)
}

# ---- 4-cell main analysis (constrained + noconstrain variants) -------
MAIN_STRATA <- c("East_pre1990", "East_post1990", "West_pre1990", "West_post1990")

results <- list()
for (suf in c("", "_noconstrain")) {
  block <- heterogeneity_block(
    analysis_label = "main_4cell",
    dir            = MAIN_DIR,
    strata         = MAIN_STRATA,
    model_suffix   = suf
  )
  if (!is.null(block)) results[[length(results) + 1L]] <- block
}

if (length(results) == 0L) {
  cat("No heterogeneity tests produced — were any standalone fits successful?\n")
  quit(status = 0)
}

final <- do.call(rbind, results)
final$timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")

write.table(final, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWrote: %s (%d rows)\n", OUT_FILE, nrow(final)))

# Compatibility symlink lrt.tsv -> heterogeneity.tsv for step 09's read
# path.
if (file.exists(LRT_LINK) && !identical(Sys.readlink(LRT_LINK), basename(OUT_FILE))) {
  unlink(LRT_LINK)
}
if (!file.exists(LRT_LINK)) {
  file.symlink(basename(OUT_FILE), LRT_LINK)
}

# Compact summary to stdout
cat("\nOmnibus Cochran's Q (constrained fits):\n")
omnibus <- final[final$test == "cochrans_Q" & final$model_suffix == "constrained",
                 c("analysis", "K", "df", "statistic", "p_value", "I2", "h2_pooled")]
print(omnibus, row.names = FALSE)

cat("\n=== Step 07 complete ===\n")
