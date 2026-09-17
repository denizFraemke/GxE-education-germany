#!/usr/bin/env Rscript
########################################################################
## diagnose_blockdiag_report.R — assemble the block-diagonal diagnostic
##
## Compares four cell-level h² estimates for East × pre-1990 (the case
## where the discrepancy is largest):
##
##   1. SHIP × East_pre1990, per-cohort × stratum fit
##      (output/reml/per_cohort/cohort_SHIP_East_pre1990.hsq)
##      — direct REML on the cohort SHIP unrel GRM with --keep
##        SHIP-East-pre individuals. This is the reference value
##        (h² ≈ 0.287).
##
##   2. SHIP × East_pre1990, one-block blockdiag replication
##      (output/reml/blockdiag_diagnostic/solo_SHIP_East_pre1990.hsq)
##      — REML on a 1-block "block-diagonal" GRM containing the SAME
##        individuals as (1), built by passing the SHIP × East-pre block
##        component through R/build_block_diagonal_grm.R. Tests whether
##        the combiner produces a valid GRM (it should reproduce (1)).
##
##   3. Full block-diagonal East × pre-1990, no cohort covariate
##      (output/reml/blockdiag/blockdiag_East_pre1990.hsq)
##      — REML on the 4-block diagonal GRM (SHIP + BASE-II + SOEP +
##        TwinLife sub-blocks), using only gender.covar.
##
##   4. Full block-diagonal East × pre-1990, with cohort covariate
##      (output/reml/blockdiag_diagnostic/blockdiag_East_pre1990_cohortcovar.hsq)
##      — REML on the same 4-block diagonal GRM but now with a
##        cohort-indicator dummy added via gender_cohort.covar.
##
## Plus the same (3 vs 4) comparison for the other three cells.
##
## Output:
##   output/final/h2_blockdiag_diagnostic.tsv  — tidy table
##   output/final/h2_blockdiag_diagnostic.md   — human-readable report
########################################################################
rm(list = ls())
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

REML_OUT   <- file.path(PROJ_ROOT, "output", "reml")
FINAL_OUT  <- file.path(PROJ_ROOT, "output", "final")
dir.create(FINAL_OUT, recursive = TRUE, showWarnings = FALSE)

OUT_TSV    <- file.path(FINAL_OUT, "h2_blockdiag_diagnostic.tsv")
OUT_REPORT <- file.path(FINAL_OUT, "h2_blockdiag_diagnostic.md")

cat("=== diagnose_blockdiag_report ===\n")

# ---- helper: parse one .hsq into a tidy row -----------------------------
parse_row <- function(hsq_path, model_label, model_description, source_dir,
                      covar_used, constrained) {
  hsq <- parse_hsq(hsq_path)
  if (!isTRUE(hsq$.present)) {
    return(data.frame(
      model_label = model_label, model_description = model_description,
      source_dir = source_dir, covar_used = covar_used,
      constrained = constrained,
      N = NA_integer_,
      V_G = NA_real_, V_G_SE = NA_real_,
      V_e = NA_real_, V_e_SE = NA_real_,
      Vp  = NA_real_,
      h2_SNP = NA_real_, h2_SNP_SE = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      converged = "missing",
      hsq_path = hsq_path,
      stringsAsFactors = FALSE
    ))
  }
  h2    <- hsq[["V(G)/Vp"]];   h2_se <- hsq[["V(G)/Vp_SE"]]
  v_g   <- hsq[["V(G)"]];      v_g_se<- hsq[["V(G)_SE"]]
  v_e   <- hsq[["V(e)"]];      v_e_se<- hsq[["V(e)_SE"]]
  v_p   <- hsq[["Vp"]]
  n     <- hsq[["n"]]
  data.frame(
    model_label = model_label, model_description = model_description,
    source_dir = source_dir, covar_used = covar_used,
    constrained = constrained,
    N = if (is.null(n)) NA_integer_ else as.integer(n),
    V_G = v_g, V_G_SE = v_g_se,
    V_e = v_e, V_e_SE = v_e_se,
    Vp  = v_p,
    h2_SNP = h2, h2_SNP_SE = h2_se,
    ci_lower = h2 - 1.96 * h2_se,
    ci_upper = h2 + 1.96 * h2_se,
    converged = "yes",   # parse_hsq presence implies GCTA wrote the .hsq
    hsq_path = hsq_path,
    stringsAsFactors = FALSE
  )
}

STRATA <- c("East_pre1990", "East_post1990", "West_pre1990", "West_post1990")

rows <- list()

# 1. Reference: SHIP × East-pre1990 per-cohort × stratum fit (from 05c)
for (suffix in c("", "_noconstrain")) {
  rows[[length(rows) + 1L]] <- parse_row(
    hsq_path = file.path(REML_OUT, "per_cohort",
                         paste0("cohort_SHIP_East_pre1990", suffix, ".hsq")),
    model_label = paste0("ref_SHIP_East_pre1990", suffix),
    model_description = "reference: SHIP × East-pre, per-cohort × stratum (cohort GRM + --keep)",
    source_dir = "output/reml/per_cohort",
    covar_used = "gender.covar",
    constrained = (suffix == "")
  )
}

# 2. Solo blockdiag replication (1-block via R combiner)
for (suffix in c("", "_noconstrain")) {
  rows[[length(rows) + 1L]] <- parse_row(
    hsq_path = file.path(REML_OUT, "blockdiag_diagnostic",
                         paste0("solo_SHIP_East_pre1990", suffix, ".hsq")),
    model_label = paste0("solo_SHIP_East_pre1990", suffix),
    model_description = "validation: 1-block blockdiag (SHIP × East-pre only, via R combiner)",
    source_dir = "output/reml/blockdiag_diagnostic",
    covar_used = "gender.covar",
    constrained = (suffix == "")
  )
}

# 3. Original full blockdiag (no cohort covar) — from step 05d
for (s in STRATA) {
  for (suffix in c("", "_noconstrain")) {
    rows[[length(rows) + 1L]] <- parse_row(
      hsq_path = file.path(REML_OUT, "blockdiag",
                           paste0("blockdiag_", s, suffix, ".hsq")),
      model_label = paste0("blockdiag_", s, suffix),
      model_description = sprintf("original full blockdiag: %s (gender only)", s),
      source_dir = "output/reml/blockdiag",
      covar_used = "gender.covar",
      constrained = (suffix == "")
    )
  }
}

# 4. Cohort-dummy refit of the full blockdiag (gender + cohort covar)
for (s in STRATA) {
  for (suffix in c("", "_noconstrain")) {
    rows[[length(rows) + 1L]] <- parse_row(
      hsq_path = file.path(REML_OUT, "blockdiag_diagnostic",
                           paste0("blockdiag_", s, "_cohortcovar", suffix, ".hsq")),
      model_label = paste0("blockdiag_", s, "_cohortcovar", suffix),
      model_description = sprintf("cohort-dummy refit: %s (gender + cohort)", s),
      source_dir = "output/reml/blockdiag_diagnostic",
      covar_used = "gender_cohort.covar",
      constrained = (suffix == "")
    )
  }
}

tab <- do.call(rbind, rows)
write.table(tab, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote: %s (%d rows)\n", OUT_TSV, nrow(tab)))

# ---- Markdown report ----------------------------------------------------
fmt   <- function(x, d = 3) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f"), x))
fmt_n <- function(n) if (is.na(n)) "NA" else formatC(n, big.mark = ",", format = "d")
fmt_ci <- function(lo, hi) {
  if (is.na(lo) || is.na(hi)) return("NA")
  sprintf("[%s, %s]", fmt(lo, 3), fmt(hi, 3))
}

# East-pre comparison (the headline of the diagnostic) — constrained rows only.
east_pre_rows <- tab[grepl("(SHIP_East_pre1990|blockdiag_East_pre1990)$", tab$model_label), ]

md <- c(
  "# Block-diagonal diagnostic report",
  "",
  sprintf("_Generated %s by `diagnose_blockdiag_report.R`._",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Compares four East × pre-1990 estimates (the cell where the gender-only block-diagonal h² is most discrepant from the SHIP-only estimate), plus a cohort-dummy refit for all four cells.",
  "",
  "## East × pre-1990: four-way comparison (constrained fits)",
  "",
  "| Model | N | V_G | V_e | V_P | h²_SNP | SE | 95% CI | Covariates | hsq |",
  "|---|---:|---:|---:|---:|---:|---:|:---:|---|---|"
)
for (i in seq_len(nrow(east_pre_rows))) {
  r <- east_pre_rows[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | **%s** | %s | %s | `%s` | `%s` |",
                      r$model_description, fmt_n(r$N),
                      fmt(r$V_G, 3), fmt(r$V_e, 3), fmt(r$Vp, 3),
                      fmt(r$h2_SNP, 3), fmt(r$h2_SNP_SE, 3),
                      fmt_ci(r$ci_lower, r$ci_upper),
                      r$covar_used,
                      basename(r$hsq_path)))
}

md <- c(md, "",
        "## Cohort-dummy refit for all four cells",
        "",
        "| Cell | N | h² (no cohort) | h² (with cohort) | Δ |",
        "|---|---:|---:|---:|---:|")
for (s in STRATA) {
  r_no  <- tab[tab$model_label == paste0("blockdiag_", s), ]
  r_yes <- tab[tab$model_label == paste0("blockdiag_", s, "_cohortcovar"), ]
  if (nrow(r_no) == 1L && nrow(r_yes) == 1L) {
    delta <- r_yes$h2_SNP - r_no$h2_SNP
    md <- c(md, sprintf("| %s | %s | %s ± %s | %s ± %s | %+s |",
                        s, fmt_n(r_yes$N),
                        fmt(r_no$h2_SNP, 3), fmt(r_no$h2_SNP_SE, 3),
                        fmt(r_yes$h2_SNP, 3), fmt(r_yes$h2_SNP_SE, 3),
                        fmt(delta, 3)))
  }
}

# ---- Interpretation block ------------------------------------------------
md <- c(md, "",
        "## Interpretation",
        "",
        "Use the decision rules in `diagnose_blockdiag.sh` to interpret:",
        "",
        "1. **Did the 1-block blockdiag replication reproduce the SHIP × East-pre reference?**",
        "   - If `solo_SHIP_East_pre1990` h² ≈ `ref_SHIP_East_pre1990` h² (both ≈ 0.287), the block-diagonal builder and the GCTA command are correctly implemented.",
        "   - If they disagree, the implementation has a bug (GRM packing, .grm.id ordering, .grm.N.bin, or covariate alignment). **Stop interpreting blockdiag results** until that's fixed.",
        "",
        "2. **Did adding a cohort dummy rescue the full-blockdiag East × pre-1990 h²?**",
        "   - If `blockdiag_East_pre1990_cohortcovar` h² rises substantially toward 0.20–0.30, the gender-only blockdiag specification is misspecified — it absorbs between-cohort education-mean differences into V_e because no cohort fixed effect is in the model.",
        "   - If it stays in 0.05–0.15, the one-component blockdiag GREML is intrinsically too restrictive for these heterogeneous cohort blocks — shared V_G and V_e across blocks with different V_P, measurement, and reliability is implausible. In that case the per-cohort × stratum estimates are the more defensible cell-level statement.",
        "",
        "## Full diagnostic table",
        "",
        "Tidy machine-readable rows are in `h2_blockdiag_diagnostic.tsv` (the constrained rows match the tables above; the `_noconstrain` rows show whether the constrained values hit the 0-boundary).")

writeLines(md, OUT_REPORT)
cat(sprintf("Wrote: %s\n", OUT_REPORT))

cat("\n=== Done ===\n")
