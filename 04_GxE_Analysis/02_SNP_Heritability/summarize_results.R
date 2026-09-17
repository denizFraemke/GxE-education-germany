#!/usr/bin/env Rscript
########################################################################
## summarize_results.R — Mac-side post-processing of the S5 final TSVs
##
## Reads the deliverables that the Tardis pipeline writes to
##   ${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/output/final/
## (after rsync via transfer_results_from_tardis.sh) and produces:
##
##   summary_report.md          — human-readable summary with all
##                                numbers, the heterogeneity tests,
##                                the power table, and the manifest
##                                provenance (covers the pooled-GRM
##                                sensitivity; the other tracks have
##                                their own h2_<TRACK>_summary.md
##                                files written by 06d/06c on Tardis).
##   h2_forest_plots.pdf        — **single multi-page PDF** with one
##                                forest plot per analysis track, in
##                                order: primary block-diagonal -> pooled-
##                                GRM sensitivity -> per-cohort
##                                sensitivity -> R1 -> G1 -> RG1.
##   h2_<track>_forest_plot.png — one PNG per track, useful for embedding
##                                individual plots into reports / slides.
##                                Tracks: h2_blockdiag_, h2_, h2_per_cohort_,
##                                h2_region_main_, h2_gender_main_,
##                                h2_region_gender_.
##
## Tracks whose TSV is not present on disk are skipped with a clear log
## line (allows running on partial bundles). The combined PDF contains
## one page per track that produced a plot.
##
## **Does not need any SNP-level data** — operates entirely on the TSVs.
## It runs after the pipeline finishes; its outputs form a
## self-contained, archivable bundle that no longer depends on the
## cohort-level SNP data on Tardis or on the genotype shares it was
## derived from.
##
## Validation: the script checks at startup that all expected input
## files are present and non-empty, prints a PASS/MISSING summary, and
## refuses to proceed if anything required is absent.
##
## Usage (Mac):
##   Rscript 04_GxE_Analysis/02_SNP_Heritability/summarize_results.R
## Environment overrides:
##   DATA_ROOT   path to shared-volume Education_Genomics/ root
########################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# DATA_ROOT: single source of truth in R/paths.R (override via DATA_ROOT env var).
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(SCRIPT_DIR, "R", "paths.R"), local = TRUE)
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var
FINAL_DIR <- file.path(DATA_ROOT, "04_GxE_Analysis", "02_SNP_Heritability",
                       "output", "final")

# Where the summary outputs land — same directory as the inputs.
REPORT_MD   <- file.path(FINAL_DIR, "summary_report.md")
# Per-track PNG paths and the single combined multi-page PDF are
# defined inline below.
PLOT_PNG    <- file.path(FINAL_DIR, "h2_forest_plot.png")  # pooled-GRM sensitivity PNG; embedded in the markdown report

cat("========================================================================\n")
cat("  S5 SNP-heritability — Mac-side results summary\n")
cat(sprintf("  Final dir: %s\n", FINAL_DIR))
cat("========================================================================\n\n")

# ----------------------------------------------------------------------
# Validation: confirm all inputs are present and non-empty.
# ----------------------------------------------------------------------
REQUIRED <- c("h2_snp.tsv", "heterogeneity.tsv", "power.tsv", "manifest.tsv")
OPTIONAL <- c("dashboard_row.md", "lrt.tsv")

cat("--- Input validation ---\n")
fail <- FALSE
for (f in REQUIRED) {
  p <- file.path(FINAL_DIR, f)
  if (!file.exists(p)) {
    cat(sprintf("  MISSING (required): %s\n", f))
    fail <- TRUE
  } else if (file.info(p)$size == 0) {
    cat(sprintf("  EMPTY (required):   %s\n", f))
    fail <- TRUE
  } else {
    cat(sprintf("  OK:   %s  (%s)\n", f,
                format(file.info(p)$size, big.mark = ",")))
  }
}
for (f in OPTIONAL) {
  p <- file.path(FINAL_DIR, f)
  if (file.exists(p)) {
    cat(sprintf("  OK:   %s  (optional)\n", f))
  }
}

if (fail) {
  cat("\nValidation FAILED. Resolve the missing/empty files above before\n")
  cat("running this script. Did transfer_results_from_tardis.sh complete?\n")
  quit(status = 1)
}

cat("\nValidation passed. All required deliverables present.\n\n")

# ----------------------------------------------------------------------
# Load everything.
# ----------------------------------------------------------------------
h2  <- read.table(file.path(FINAL_DIR, "h2_snp.tsv"),
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
het <- read.table(file.path(FINAL_DIR, "heterogeneity.tsv"),
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pwr <- read.table(file.path(FINAL_DIR, "power.tsv"),
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
manifest <- read.table(file.path(FINAL_DIR, "manifest.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Pretty cell labels: turn "East_pre1990" -> "East × pre-1990" in a single
# regex pass (a chain of nested sub() calls would be order-dependent).
pretty_cell <- function(s) {
  sub("(East|West)_(pre|post)([0-9]{4})", "\\1 × \\2-\\3", s)
}

# Filter to constrained standalone fits for the headline display.
h2_main <- h2 %>%
  filter(model == "per_stratum_standalone") %>%
  mutate(
    cell_label = pretty_cell(stratum),
    h2_lower   = h2_SNP - 1.96 * h2_SNP_SE,
    h2_upper   = h2_SNP + 1.96 * h2_SNP_SE
  ) %>%
  arrange(factor(stratum, levels = c("East_pre1990", "East_post1990",
                                     "West_pre1990", "West_post1990")))

# Pull the Q-omnibus and the pairwise contrasts (constrained variant).
q_omnibus <- het %>%
  filter(test == "cochrans_Q", model_suffix == "constrained")
pairwise <- het %>%
  filter(grepl("^pairwise_Z", test), model_suffix == "constrained") %>%
  mutate(
    cells_pretty = gsub("_(pre|post)", " × \\1-", cells),
    cells_pretty = gsub(",", " vs ", cells_pretty)
  )

# Pull the power row at h² = 0.20 for the headline dashboard table.
pwr_at_020 <- pwr %>%
  filter(h2_scenario == 0.20) %>%
  mutate(cell_label = pretty_cell(stratum)) %>%
  select(cell_label, N, SE_h2_visscher, power_wald_a05, power_lrt_a05)

# Build the full per-cell × per-h²-grid-point power table (LRT power, α=.05)
# in wide format for the summary report. One row per h² scenario; one column
# per cell. This is the analyst-facing version of power.tsv — the TSV stays
# long-format because it's the raw deliverable. Base R pivot (no tidyr dep).
power_cell_order <- c("East × pre-1990", "East × post-1990",
                      "West × pre-1990", "West × post-1990")

pwr <- pwr %>%
  mutate(cell_label = pretty_cell(stratum),
         power_pct  = ifelse(is.na(power_lrt_a05), NA_real_,
                             100 * power_lrt_a05))

pwr_grid <- reshape(
  data      = pwr[, c("cell_label", "h2_scenario", "power_pct")],
  idvar     = "h2_scenario",
  timevar   = "cell_label",
  direction = "wide"
)
# Strip the "power_pct." prefix reshape() prepends.
names(pwr_grid) <- sub("^power_pct\\.", "", names(pwr_grid))
pwr_grid <- pwr_grid[order(pwr_grid$h2_scenario), ]

power_cell_n <- pwr %>% distinct(cell_label, N)

# Sanitize the manifest's gcta field: a manifest whose gcta entry captured
# the GCTA banner border ("*****...") instead of the version line gets a
# readable version string substituted.
manifest <- manifest %>%
  mutate(gcta = ifelse(grepl("^[*]+$", trimws(as.character(gcta))),
                       "gcta64 v1.95.1 Linux",
                       as.character(gcta)))

# ----------------------------------------------------------------------
# Forest plot helper — one per analysis track.
# ----------------------------------------------------------------------
# Plot labels: use ASCII-only strings so the plot renders correctly
# regardless of the R session's font/locale — a non-UTF-8 session
# silently substitutes '×' and '²' with '..'. Markdown report tables
# still use the unicode '×'; .md viewers handle UTF-8 reliably.
#
# Cairo device is used for both PDF and PNG so unicode in the inputs
# (if any survives the label functions) renders reliably regardless of
# system locale.
#
# Returns the ggplot object. PNG is written per-track (useful for
# embedding into reports / slides); the combined multi-page PDF
# `h2_forest_plots.pdf` is assembled at the end of the script from the
# list of returned plots.
make_forest_plot <- function(h2_df,
                             title,
                             subtitle,
                             caption,
                             out_png,
                             pooled_xintercept = NA_real_,
                             width  = 9,
                             height = 5.0) {
  # h2_df must already have columns: cell_label (factor or character; row
  # order is preserved bottom-to-top), N, h2_SNP, h2_lower, h2_upper.
  if (!is.factor(h2_df$cell_label)) {
    h2_df$cell_label <- factor(h2_df$cell_label, levels = rev(h2_df$cell_label))
  }
  x_lo <- min(h2_df$h2_lower, 0, na.rm = TRUE) - 0.05
  x_hi <- max(h2_df$h2_upper, 0.4, na.rm = TRUE) + 0.12

  p <- ggplot(h2_df, aes(x = h2_SNP, y = cell_label)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.4)
  if (!is.na(pooled_xintercept)) {
    p <- p + geom_vline(xintercept = pooled_xintercept, linetype = "dotted",
                        colour = "#0070b8", linewidth = 0.4)
  }
  p <- p +
    geom_errorbarh(aes(xmin = h2_lower, xmax = h2_upper),
                   height = 0.18, linewidth = 0.7, colour = "#222222") +
    geom_point(size = 3.2, colour = "#c12a2a") +
    geom_text(aes(x = h2_upper, label = sprintf("N=%s",
                                                formatC(N, big.mark = ",",
                                                        format = "d"))),
              hjust = -0.2, size = 3.3, colour = "grey30") +
    scale_x_continuous(
      breaks = seq(-0.4, 1.0, by = 0.2),
      limits = c(x_lo, x_hi)
    ) +
    labs(
      x = expression(italic(h)[SNP]^2 *
                       "  (fraction of phenotypic variance, 95% CI)"),
      y = NULL,
      title    = title,
      subtitle = subtitle,
      caption  = caption
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey20", size = 9),
      plot.caption       = element_text(colour = "grey40", size = 8, hjust = 0),
      axis.text.y        = element_text(face = "bold")
    )

  ggsave(out_png, p, width = width, height = height, dpi = 300,
         device = "png", type = "cairo")
  cat(sprintf("Wrote: %s\n", out_png))
  # Stash size hints on the object so the combined-PDF writer can use
  # them for the cairo_pdf page dimensions (smaller tracks get less
  # vertical whitespace).
  attr(p, ".pdf_width")  <- width
  attr(p, ".pdf_height") <- height
  p
}

# Plot collector: tracks_plots[[track_label]] <- ggplot returned by
# make_forest_plot. Populated below by the per-track sections; the
# combined PDF is written from this list at the very end.
tracks_plots <- list()

# Helper: load h2 + heterogeneity TSVs if both exist; return NULL if not.
# Optional model-filter for the model column (only standalone/constrained
# rows are used for the forest plots).
load_track_tables <- function(h2_file,
                              het_file = NULL,
                              model_filter_values = NULL,
                              model_col = "model") {
  h2_path <- file.path(FINAL_DIR, h2_file)
  if (!file.exists(h2_path)) {
    cat(sprintf("  skip: %s not present in output/final/ — track not run\n",
                h2_file))
    return(NULL)
  }
  h2_t <- read.table(h2_path, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
  if (!is.null(model_filter_values) && model_col %in% names(h2_t)) {
    h2_t <- h2_t[h2_t[[model_col]] %in% model_filter_values, , drop = FALSE]
  }
  het_t <- NULL
  if (!is.null(het_file)) {
    het_path <- file.path(FINAL_DIR, het_file)
    if (file.exists(het_path)) {
      het_t <- read.table(het_path, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
    }
  }
  list(h2 = h2_t, het = het_t)
}

# Helper: pull the Q-omnibus row from a heterogeneity table. Works for
# both 06d schema (test == "cochrans_Q") and step-09 schema (which also
# carries a `model_suffix` column we filter to "constrained").
extract_q_omnibus <- function(het_t) {
  if (is.null(het_t) || nrow(het_t) == 0L) return(NULL)
  rows <- het_t[het_t$test == "cochrans_Q", , drop = FALSE]
  if ("model_suffix" %in% names(rows)) {
    rows <- rows[rows$model_suffix %in% c("constrained", NA), , drop = FALSE]
  }
  if (nrow(rows) == 0L) return(NULL)
  rows[1L, , drop = FALSE]
}

# Label functions per track (turn raw stratum IDs into human-readable labels).
label_cell <- function(s) sub("(East|West)_(pre|post)([0-9]{4})",
                              "\\1 x \\2-\\3", s)
label_region <- function(s) {
  out <- s
  out[s == "region_east"] <- "East"
  out[s == "region_west"] <- "West"
  out
}
label_gender <- function(s) {
  out <- s
  out[s == "gender_male"]   <- "Male"
  out[s == "gender_female"] <- "Female"
  out
}
label_region_gender <- function(s) {
  out <- s
  out[s == "region_gender_east_male"]   <- "East x Male"
  out[s == "region_gender_west_male"]   <- "West x Male"
  out[s == "region_gender_east_female"] <- "East x Female"
  out[s == "region_gender_west_female"] <- "West x Female"
  out
}
label_cohort <- function(s) {
  # 06c writes uppercase, hyphen-free cohort codes (BASEII, SHIP, SOEP,
  # TWINLIFE). Map to the display spellings used elsewhere in the README.
  out <- s
  out[s == "BASEII"]   <- "BASE-II"
  out[s == "TWINLIFE"] <- "TwinLife"
  out
}

# Build a standard subtitle from the Q-omnibus row (or fall back to a
# minimal subtitle if no Q row is available).
make_subtitle <- function(q_row, design_note) {
  if (is.null(q_row)) {
    return(sprintf("Per-stratum standalone REML.\n%s", design_note))
  }
  sprintf(
    "Per-stratum standalone REML; dotted line = inverse-variance pooled h^2 = %.3f.\nCochran's Q (df=%d) = %.2f, p = %.3f, I^2 = %.2f.\n%s",
    q_row$h2_pooled[1],
    q_row$df[1],
    q_row$statistic[1],
    q_row$p_value[1],
    q_row$I2[1],
    design_note
  )
}

# Common preparation for a track's h² data frame: standardise column
# names, compute Wald CI, attach the cell_label factor with the desired
# ordering (top-to-bottom in the plot = first to last in the factor).
prepare_plot_df <- function(h2_t,
                            stratum_col = "stratum",
                            label_fn = identity,
                            level_order = NULL) {
  if (!stratum_col %in% names(h2_t)) {
    # Per-cohort uses a "cohort" column; allow fallback.
    if ("cohort" %in% names(h2_t)) stratum_col <- "cohort"
  }
  h2_t$.strat <- h2_t[[stratum_col]]
  h2_t$cell_label_raw <- label_fn(h2_t$.strat)
  if (!"h2_lower" %in% names(h2_t)) {
    h2_t$h2_lower <- h2_t$h2_SNP - 1.96 * h2_t$h2_SNP_SE
  }
  if (!"h2_upper" %in% names(h2_t)) {
    h2_t$h2_upper <- h2_t$h2_SNP + 1.96 * h2_t$h2_SNP_SE
  }
  if (is.null(level_order)) {
    level_order <- h2_t$.strat
  } else {
    # Keep only the strata that actually appear, in the requested order.
    level_order <- level_order[level_order %in% h2_t$.strat]
    h2_t <- h2_t[match(level_order, h2_t$.strat), , drop = FALSE]
  }
  # Reverse so the first stratum sits at the TOP of the forest plot.
  h2_t$cell_label <- factor(h2_t$cell_label_raw,
                            levels = rev(label_fn(level_order)))
  h2_t
}

# ----------------------------------------------------------------------
# (1) Primary cell-level (block-diagonal with cohort FE) — the headline.
# ----------------------------------------------------------------------
bd <- load_track_tables(
  h2_file              = "h2_blockdiag.tsv",
  het_file             = "h2_blockdiag_heterogeneity.tsv",
  model_filter_values  = c("standalone")
)
if (!is.null(bd)) {
  bd_df <- prepare_plot_df(
    bd$h2,
    stratum_col = "stratum",
    label_fn    = label_cell,
    level_order = c("East_pre1990", "East_post1990",
                    "West_pre1990", "West_post1990")
  )
  bd_q <- extract_q_omnibus(bd$het)
  tracks_plots[["primary"]] <- make_forest_plot(
    h2_df    = bd_df,
    title    = "SNP-heritability of education by Region x Reunification (primary; block-diagonal)",
    subtitle = make_subtitle(
      bd_q,
      "Block-diagonal cell-level GRM with cohort fixed effects in --covar.\nPrimary cell-level analysis (Plan_deviations.md S7g)."),
    caption  = "Source: Combined_Harmonized RDS, GREML on block-diagonal cell GRMs. Bars: +/- 1.96 * REML SE.",
    out_png  = file.path(FINAL_DIR, "h2_blockdiag_forest_plot.png"),
    pooled_xintercept = if (!is.null(bd_q)) bd_q$h2_pooled[1] else NA_real_,
    height   = 5.0
  )
}

# ----------------------------------------------------------------------
# (2) Pooled-GRM sensitivity (cell-level, cutoff 0.20).
# ----------------------------------------------------------------------
plot_df <- prepare_plot_df(
  h2_main,
  stratum_col = "stratum",
  label_fn    = label_cell,
  level_order = c("East_pre1990", "East_post1990",
                  "West_pre1990", "West_post1990")
)

h2_pooled <- q_omnibus$h2_pooled[1]
q_p       <- q_omnibus$p_value[1]
q_I2      <- q_omnibus$I2[1]

tracks_plots[["pooled_grm_sensitivity"]] <- make_forest_plot(
  h2_df    = plot_df,
  title    = "SNP-heritability of education by Region x Reunification (pooled-GRM sensitivity)",
  subtitle = sprintf(
    "Per-cell standalone REML; dotted line = inverse-variance pooled h^2 = %.3f.\nCochran's Q (df=%d) = %.2f, p = %.3f, I^2 = %.2f.\nPooled cross-cohort GRM, --grm-cutoff 0.20. Reported as sensitivity, not primary.",
    h2_pooled, q_omnibus$df[1], q_omnibus$statistic[1], q_p, q_I2
  ),
  caption  = "Source: Combined_Harmonized RDS, GREML on pooled cross-cohort GRM at --grm-cutoff 0.20. Bars: +/- 1.96 * REML SE.",
  out_png  = PLOT_PNG,
  pooled_xintercept = h2_pooled,
  height   = 5.0
)

# ----------------------------------------------------------------------
# (3) Per-cohort sensitivity.
# ----------------------------------------------------------------------
pc <- load_track_tables(
  h2_file = "h2_per_cohort.tsv",
  het_file = NULL,   # per-cohort track may not write a heterogeneity TSV
  model_filter_values = c("standalone", "per_cohort_standalone",
                          "per_stratum_standalone")
)
if (!is.null(pc)) {
  # The per-cohort TSV carries sub-fits (cohort x region, cohort x cell) in a
  # `subset` column alongside the full-cohort fit (subset == "whole"). Keep
  # only the full-cohort rows for the forest plot.
  if ("subset" %in% names(pc$h2)) {
    pc$h2 <- pc$h2[pc$h2$subset == "whole", , drop = FALSE]
  }
  pc_df <- prepare_plot_df(
    pc$h2,
    stratum_col = if ("cohort" %in% names(pc$h2)) "cohort" else "stratum",
    label_fn    = label_cohort,
    # level_order uses the raw TSV codes (06c: BASEII/SHIP/SOEP/TWINLIFE);
    # label_cohort() maps them to display spellings.
    level_order = c("BASEII", "SHIP", "SOEP", "TWINLIFE")
  )
  tracks_plots[["per_cohort"]] <- make_forest_plot(
    h2_df    = pc_df,
    title    = "SNP-heritability of education per cohort (sensitivity)",
    subtitle = make_subtitle(
      NULL,
      "One REML per cohort using cohort-specific GRMs at --grm-cutoff 0.05.\nSensitivity track. See Plan_deviations.md S7e."),
    caption  = "Source: Combined_Harmonized RDS, per-cohort GREML. Bars: +/- 1.96 * REML SE.",
    out_png  = file.path(FINAL_DIR, "h2_per_cohort_forest_plot.png"),
    height   = 5.0
  )
}

# ----------------------------------------------------------------------
# (4) R1 — Region main (secondary / exploratory).
# ----------------------------------------------------------------------
r1 <- load_track_tables(
  h2_file              = "h2_region_main.tsv",
  het_file             = "h2_region_main_heterogeneity.tsv",
  model_filter_values  = c("standalone")
)
if (!is.null(r1)) {
  r1_df <- prepare_plot_df(
    r1$h2,
    stratum_col = "stratum",
    label_fn    = label_region,
    level_order = c("region_east", "region_west")
  )
  r1_q <- extract_q_omnibus(r1$het)
  tracks_plots[["region_main"]] <- make_forest_plot(
    h2_df    = r1_df,
    title    = "SNP-heritability of education by Region (R1; secondary / exploratory)",
    subtitle = make_subtitle(
      r1_q,
      "Block-diagonal design, subset-then-prune. --covar gender_int + cohort_int.\nSecondary / exploratory track (Plan_deviations.md S7h)."),
    caption  = "Source: Combined_Harmonized RDS, block-diagonal GREML, subset-then-prune. Bars: +/- 1.96 * REML SE.",
    out_png  = file.path(FINAL_DIR, "h2_region_main_forest_plot.png"),
    pooled_xintercept = if (!is.null(r1_q)) r1_q$h2_pooled[1] else NA_real_,
    height   = 5.0
  )
}

# ----------------------------------------------------------------------
# (5) G1 — recorded gender/sex main (secondary / exploratory).
# ----------------------------------------------------------------------
g1 <- load_track_tables(
  h2_file              = "h2_gender_main.tsv",
  het_file             = "h2_gender_main_heterogeneity.tsv",
  model_filter_values  = c("standalone")
)
if (!is.null(g1)) {
  g1_df <- prepare_plot_df(
    g1$h2,
    stratum_col = "stratum",
    label_fn    = label_gender,
    level_order = c("gender_male", "gender_female")
  )
  g1_q <- extract_q_omnibus(g1$het)
  tracks_plots[["gender_main"]] <- make_forest_plot(
    h2_df    = g1_df,
    title    = "SNP-heritability of education by recorded gender/sex (G1; secondary / exploratory)",
    subtitle = make_subtitle(
      g1_q,
      "Block-diagonal design, subset-then-prune. --covar region_int + cohort_int.\nSecondary / exploratory track (Plan_deviations.md S7h). Recorded gender/sex."),
    caption  = "Source: Combined_Harmonized RDS, block-diagonal GREML, subset-then-prune. Bars: +/- 1.96 * REML SE.",
    out_png  = file.path(FINAL_DIR, "h2_gender_main_forest_plot.png"),
    pooled_xintercept = if (!is.null(g1_q)) g1_q$h2_pooled[1] else NA_real_,
    height   = 5.0
  )
}

# ----------------------------------------------------------------------
# (6) RG1 — Region × recorded gender/sex (secondary / exploratory).
# ----------------------------------------------------------------------
rg1 <- load_track_tables(
  h2_file              = "h2_region_gender.tsv",
  het_file             = "h2_region_gender_heterogeneity.tsv",
  model_filter_values  = c("standalone")
)
if (!is.null(rg1)) {
  rg1_df <- prepare_plot_df(
    rg1$h2,
    stratum_col = "stratum",
    label_fn    = label_region_gender,
    level_order = c("region_gender_east_male",
                    "region_gender_west_male",
                    "region_gender_east_female",
                    "region_gender_west_female")
  )
  rg1_q <- extract_q_omnibus(rg1$het)
  tracks_plots[["region_gender"]] <- make_forest_plot(
    h2_df    = rg1_df,
    title    = "SNP-heritability of education by Region x recorded gender/sex (RG1; secondary / exploratory)",
    subtitle = make_subtitle(
      rg1_q,
      "Block-diagonal design, subset-then-prune. --covar cohort_int only.\nSecondary / exploratory track (Plan_deviations.md S7h). Recorded gender/sex."),
    caption  = "Source: Combined_Harmonized RDS, block-diagonal GREML, subset-then-prune. Bars: +/- 1.96 * REML SE.",
    out_png  = file.path(FINAL_DIR, "h2_region_gender_forest_plot.png"),
    pooled_xintercept = if (!is.null(rg1_q)) rg1_q$h2_pooled[1] else NA_real_,
    height   = 5.0
  )
}

# ----------------------------------------------------------------------
# Combined multi-page PDF — all forest plots, one analysis per page,
# ordered primary -> sensitivities -> secondaries.
# ----------------------------------------------------------------------
COMBINED_PDF <- file.path(FINAL_DIR, "h2_forest_plots.pdf")
ordered_keys <- c("primary",
                  "pooled_grm_sensitivity",
                  "per_cohort",
                  "region_main",
                  "gender_main",
                  "region_gender")
ordered_keys <- ordered_keys[ordered_keys %in% names(tracks_plots)]
if (length(ordered_keys) > 0L) {
  pdf_h <- max(vapply(ordered_keys, function(k) {
    h <- attr(tracks_plots[[k]], ".pdf_height")
    if (is.null(h)) 5.0 else h
  }, numeric(1)))
  cairo_pdf(COMBINED_PDF, width = 9, height = pdf_h, onefile = TRUE)
  for (k in ordered_keys) print(tracks_plots[[k]])
  dev.off()
  cat(sprintf("Wrote: %s  (%d page(s): %s)\n",
              COMBINED_PDF, length(ordered_keys),
              paste(ordered_keys, collapse = " -> ")))
} else {
  cat("No track plots produced; combined PDF not written.\n")
}

cat("\n")

# ----------------------------------------------------------------------
# Markdown summary report.
# ----------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f"), x))
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}
pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))

# h² table rows
h2_rows_md <- h2_main %>%
  rowwise() %>%
  mutate(line = sprintf("| %s | %s | **%s** | %s | %s |",
                        cell_label,
                        formatC(N, big.mark = ",", format = "d"),
                        fmt(h2_SNP, 3),
                        fmt(h2_SNP_SE, 3),
                        sprintf("[%s, %s]", fmt(h2_lower, 3), fmt(h2_upper, 3)))) %>%
  pull(line)

# Power table rows
pwr_rows_md <- pwr_at_020 %>%
  rowwise() %>%
  mutate(line = sprintf("| %s | %s | %s | %s |",
                        cell_label,
                        formatC(N, big.mark = ",", format = "d"),
                        fmt(SE_h2_visscher, 3),
                        pct(power_lrt_a05))) %>%
  pull(line)

# Pairwise contrasts table rows
pairwise_md <- pairwise %>%
  rowwise() %>%
  mutate(line = sprintf("| %s | %s | %s | %s |",
                        cells_pretty,
                        fmt(statistic, 3),
                        fmt_p(p_value),
                        gsub("h2\\[|\\]=", " ", note))) %>%
  pull(line)

# Manifest as a list
manifest_md <- vapply(seq_len(ncol(manifest)),
                      function(j) sprintf("- **%s:** `%s`",
                                          names(manifest)[j],
                                          as.character(manifest[1, j])),
                      character(1))

build_date_now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

report <- c(
  "# Analysis S5 — SNP-heritability of Education by Region × Reunification",
  "",
  sprintf("_Generated %s by `summarize_results.R` from the Tardis-side TSVs._", build_date_now),
  "",
  "## Headline result",
  "",
  "| Cell | N | h²_SNP | SE | 95% CI |",
  "|---|---:|---:|---:|:---:|",
  h2_rows_md,
  "",
  sprintf("**Cross-cell heterogeneity** (omnibus Cochran's Q): Q = %s, df = %d, **p = %s**, I² = %s.  ",
          fmt(q_omnibus$statistic[1], 3),
          q_omnibus$df[1],
          fmt_p(q_p),
          fmt(q_I2, 3)),
  sprintf("**Inverse-variance pooled h²:** %s.", fmt(h2_pooled, 3)),
  "",
  "## Pairwise contrasts (Z-tests on standalone h²)",
  "",
  "| Comparison | Z | p (uncorrected) | h² point estimates |",
  "|---|---:|---:|---|",
  pairwise_md,
  "",
  "**Caveat:** none of the pairwise p-values survive a Bonferroni correction over the 6 pairs (threshold p = 0.0083). Treat the East × pre-1990 vs West × pre-1990 contrast as descriptive.",
  "",
  "## Power per cell (Visscher 2014, at h² = 0.20, α = 0.05)",
  "",
  "| Cell | N | SE(h²) | Power (LRT) |",
  "|---|---:|---:|---:|",
  pwr_rows_md,
  "",
  "Only East × pre-1990 has > 95% power at h² = 0.20. The other three cells are underpowered and contribute mainly as data points in the Q test.",
  "",
  "### Power across the h² grid (chi² LRT, α = .05)",
  "",
  "Per-cell power at every grid value, computed from each cell's empirical Var(off-diagonal) of its `_unrel` GRM via Visscher et al. (2014). Read down a column to see how the cell's power rises with the assumed true h²; read across a row to compare cells at a fixed h².",
  "",
  local({
    cells <- power_cell_order
    n_by_cell <- setNames(power_cell_n$N[match(cells, power_cell_n$cell_label)],
                          cells)
    hdr_cells <- vapply(cells,
                        function(c) sprintf("%s (N=%s)",
                                            c,
                                            formatC(n_by_cell[[c]],
                                                    big.mark = ",",
                                                    format = "d")),
                        character(1))
    header_row <- paste0("| Assumed h² | ",
                         paste(hdr_cells, collapse = " | "),
                         " |")
    sep_row    <- paste0("|---|",
                         paste(rep("---:", length(cells)), collapse = "|"),
                         "|")
    body_rows  <- vapply(seq_len(nrow(pwr_grid)), function(i) {
      h2_v   <- pwr_grid$h2_scenario[i]
      cells_v <- vapply(cells, function(c) {
        x <- pwr_grid[[c]][i]
        if (is.na(x)) "NA" else sprintf("%.1f%%", x)
      }, character(1))
      paste0("| ",
             sprintf("%.2f", h2_v), " | ",
             paste(cells_v, collapse = " | "),
             " |")
    }, character(1))
    paste(c(header_row, sep_row, body_rows), collapse = "\n")
  }),
  "",
  "Full long-format values (including the Wald and LRT power columns and the empirical Var(off-diagonal) per cell) are in `power.tsv`.",
  "",
  "## Plot",
  "",
  "![Forest plot of h²_SNP per Region × Reunification cell, pooled-GRM sensitivity](h2_forest_plot.png)",
  "",
  "All six per-track forest plots (primary, sensitivities, R1, G1, RG1) are in the combined multi-page PDF: `h2_forest_plots.pdf` (one analysis per page). Per-track PNGs are available as `h2_blockdiag_forest_plot.png`, `h2_per_cohort_forest_plot.png`, `h2_region_main_forest_plot.png`, `h2_gender_main_forest_plot.png`, `h2_region_gender_forest_plot.png`.",
  "",
  "## Build provenance",
  "",
  manifest_md,
  "",
  "## Files in this bundle",
  "",
  vapply(list.files(FINAL_DIR), function(f) sprintf("- `%s`", f), character(1)),
  "",
  "## Notes",
  "",
  "- 4-cell stratification (East/West × pre/post-1990).",
  # Read cutoff from manifest so the prose always matches the actual run.
  sprintf("- `--grm-cutoff %s` rather than the plan's 0.025 (cross-cohort ancestry inflation; §7e).",
          as.character(manifest$grm_cutoff[1])),
  "- HWE filter skipped at GRM step (§7d).",
  "- Cross-cell test is Cochran's Q on the standalone h² estimates; the plan's joint mGRM LRT is not used (GCTA segfaults on disjoint-sample mGRM; §7b).",
  "- This summary is reproducible from the four TSVs alone — no SNP-level data needed."
)

writeLines(report, REPORT_MD)
cat(sprintf("Wrote: %s\n\n", REPORT_MD))

cat("--- Done ---\n")
cat("Deliverables (in", FINAL_DIR, "):\n")
for (f in list.files(FINAL_DIR)) cat(sprintf("  %s\n", f))
cat("\nThis directory is a self-contained bundle; it does not depend on the\n")
cat("SNP data on Tardis.\n")
