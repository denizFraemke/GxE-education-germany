#!/usr/bin/env Rscript
########################################################################
## 06c_assemble_per_cohort.R — Per-cohort h² output (Tardis)
##
## Walks output/reml/per_cohort/ for the .hsq files written by step 05c
## and produces a tidy TSV in output/final/h2_per_cohort.tsv with one row
## per (cohort, subset, model) combination.
##
## "subset" ∈ {whole, East, West, East_pre1990, East_post1990,
##             West_pre1990, West_post1990}
## "model"  ∈ {standalone, standalone_noconstrain}
##
## Also writes output/final/h2_per_cohort_summary.md — a small human-
## readable markdown report.
##
## Aggregate output only (counts, h², SE, CI) — no individual data.
########################################################################
rm(list = ls())

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))

source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

COHORT_DIR <- file.path(PROJ_ROOT, "output", "reml", "per_cohort")
FINAL_OUT  <- file.path(PROJ_ROOT, "output", "final")
dir.create(FINAL_OUT, recursive = TRUE, showWarnings = FALSE)

OUT_TSV     <- file.path(FINAL_OUT, "h2_per_cohort.tsv")
OUT_REPORT  <- file.path(FINAL_OUT, "h2_per_cohort_summary.md")

cat("=== Step 06c: Assemble per-cohort h² outputs ===\n")
cat(sprintf("  REML dir: %s\n", COHORT_DIR))
cat(sprintf("  Out TSV : %s\n", OUT_TSV))

if (!dir.exists(COHORT_DIR)) {
  cat("\nNo per-cohort REML directory found — run step 05c first.\n")
  quit(status = 1)
}

hsq_files <- list.files(COHORT_DIR, pattern = "\\.hsq$", full.names = TRUE)
if (length(hsq_files) == 0L) {
  cat("\nNo .hsq files in", COHORT_DIR, "— run step 05c first.\n")
  quit(status = 1)
}
cat(sprintf("  Found %d .hsq files.\n\n", length(hsq_files)))

# ---- Parse each .hsq into a tidy row ---------------------------------
# Label format from 05c:
#   cohort_<COHORT>[_<SUBSET>][_noconstrain]
# where <COHORT> ∈ {BASEII, SHIP, SOEP, TWINLIFE}
#   and <SUBSET> ∈ {<nothing>, East, West, East_pre1990, ...}
COHORTS <- c("BASEII", "SHIP", "SOEP", "TWINLIFE")
parse_label <- function(label) {
  noconstrain <- grepl("_noconstrain$", label)
  base        <- sub("_noconstrain$", "", label)
  # base = cohort_<COHORT>[_<SUBSET>]
  stopifnot(grepl("^cohort_", base))
  rest <- sub("^cohort_", "", base)
  # First token = COHORT
  m <- regexpr(paste0("^(", paste(COHORTS, collapse = "|"), ")"), rest)
  if (m[1] != 1L)
    stop("Could not parse cohort from label: ", label)
  cohort <- regmatches(rest, m)
  rest   <- sub(paste0("^", cohort), "", rest)
  subset <- if (nzchar(rest)) sub("^_", "", rest) else "whole"
  list(cohort = cohort, subset = subset,
       model  = if (noconstrain) "standalone_noconstrain" else "standalone")
}

rows <- list()
for (f in hsq_files) {
  label <- sub("\\.hsq$", "", basename(f))
  meta  <- parse_label(label)
  hsq   <- parse_hsq(f)
  if (!isTRUE(hsq$.present)) {
    cat(sprintf("  WARN: %s present but parseable as no fit\n", label))
    next
  }
  v_g  <- hsq[["V(G)"]];   v_g_se  <- hsq[["V(G)_SE"]]
  v_e  <- hsq[["V(e)"]];   v_e_se  <- hsq[["V(e)_SE"]]
  v_p  <- hsq[["Vp"]]
  h2   <- hsq[["V(G)/Vp"]]
  h2_se <- hsq[["V(G)/Vp_SE"]]
  n    <- hsq[["n"]]
  if (is.null(h2) || is.null(h2_se) || is.null(n)) {
    cat(sprintf("  WARN: %s missing required fields (h2/SE/n)\n", label))
    next
  }
  rows[[length(rows) + 1L]] <- data.frame(
    label    = label,
    cohort   = meta$cohort,
    subset   = meta$subset,
    model    = meta$model,
    N        = as.integer(n),
    V_G      = v_g,    V_G_SE   = v_g_se,
    V_e      = v_e,    V_e_SE   = v_e_se,
    Vp       = v_p,
    h2_SNP   = h2,     h2_SNP_SE = h2_se,
    ci_lower = h2 - 1.96 * h2_se,
    ci_upper = h2 + 1.96 * h2_se,
    stringsAsFactors = FALSE
  )
}

if (length(rows) == 0L) {
  stop("No parseable .hsq rows found.")
}

h2_tab <- do.call(rbind, rows)

# Stable ordering: cohort × subset × model.
SUBSET_ORDER <- c("whole",
                  "East", "West",
                  "East_pre1990", "East_post1990",
                  "West_pre1990", "West_post1990")
h2_tab$subset <- factor(h2_tab$subset, levels = SUBSET_ORDER)
h2_tab$cohort <- factor(h2_tab$cohort, levels = COHORTS)
h2_tab <- h2_tab[order(h2_tab$cohort, h2_tab$subset, h2_tab$model), ]
h2_tab$subset <- as.character(h2_tab$subset)
h2_tab$cohort <- as.character(h2_tab$cohort)

write.table(h2_tab, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote: %s (%d rows)\n", OUT_TSV, nrow(h2_tab)))

# ---- Human-readable summary ------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f"), x))
fmt_n <- function(n) formatC(n, big.mark = ",", format = "d")

md <- c(
  "# Per-cohort h²_SNP sensitivity (step 06c)",
  "",
  sprintf("_Generated %s by `06c_assemble_per_cohort.R`._",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Per-cohort GREML using cohort-specific GRMs (no cross-cohort fingerprint by construction).",
  "Each cohort's GRM was filtered at `--grm-cutoff 0.05` within-cohort (the published convention).",
  "",
  "## Cohort-only h² (full cohort)",
  "",
  "| Cohort | N | h²_SNP | SE | 95% CI |",
  "|---|---:|---:|---:|:---:|"
)

cohort_only <- h2_tab[h2_tab$subset == "whole" & h2_tab$model == "standalone", ]
for (i in seq_len(nrow(cohort_only))) {
  r <- cohort_only[i, ]
  md <- c(md, sprintf("| %s | %s | **%s** | %s | [%s, %s] |",
                     r$cohort, fmt_n(r$N), fmt(r$h2_SNP, 3),
                     fmt(r$h2_SNP_SE, 3), fmt(r$ci_lower, 3), fmt(r$ci_upper, 3)))
}

md <- c(md, "",
        "## Cohort × region h² (where N ≥ 100 in the unrel GRM × region subset)",
        "",
        "| Cohort | Region | N | h²_SNP | SE | 95% CI |",
        "|---|---|---:|---:|---:|:---:|")
sub_region <- h2_tab[h2_tab$subset %in% c("East", "West") & h2_tab$model == "standalone", ]
if (nrow(sub_region) > 0L) {
  for (i in seq_len(nrow(sub_region))) {
    r <- sub_region[i, ]
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | [%s, %s] |",
                       r$cohort, r$subset, fmt_n(r$N),
                       fmt(r$h2_SNP, 3), fmt(r$h2_SNP_SE, 3),
                       fmt(r$ci_lower, 3), fmt(r$ci_upper, 3)))
  }
} else {
  md <- c(md, "_None._")
}

md <- c(md, "",
        "## Cohort × stratum h² (where N ≥ 100 in the unrel GRM × stratum subset)",
        "",
        "| Cohort | Stratum | N | h²_SNP | SE | 95% CI |",
        "|---|---|---:|---:|---:|:---:|")
sub_strat <- h2_tab[h2_tab$subset %in% c("East_pre1990", "East_post1990",
                                         "West_pre1990", "West_post1990") &
                    h2_tab$model == "standalone", ]
if (nrow(sub_strat) > 0L) {
  for (i in seq_len(nrow(sub_strat))) {
    r <- sub_strat[i, ]
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | [%s, %s] |",
                       r$cohort, r$subset, fmt_n(r$N),
                       fmt(r$h2_SNP, 3), fmt(r$h2_SNP_SE, 3),
                       fmt(r$ci_lower, 3), fmt(r$ci_upper, 3)))
  }
} else {
  md <- c(md, "_None._")
}

md <- c(md, "",
        "## Interpretation",
        "",
        "If all four cohorts produce overlapping h²_SNP confidence intervals around the pooled headline (~0.29 from the cell-level analysis), no single cohort is driving the result.",
        "",
        "If one cohort's h² is markedly different (CIs do not overlap with the others), that cohort contributes differently to the pooled cell-level estimate and the pooled value should be read as cohort-weighted.",
        "",
        "All `_noconstrain` fits are in `h2_per_cohort.tsv` for inspection but are not displayed here — they're typically identical to the constrained fits when the constrained estimate stays comfortably away from the [0, 1] boundary.")

writeLines(md, OUT_REPORT)
cat(sprintf("Wrote: %s\n", OUT_REPORT))

cat("\n=== Step 06c complete ===\n")
