# B_sensitivities_report.R
# Section B sensitivities — report stage. Reads B_Sensitivities_Models.rds and
# builds Supplementary_Data_6_Mobility_sensitivities_results.xlsx,
# B_Sensitivities_Plots.pdf and B_Sensitivities_Overview.md.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble); library(openxlsx)
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
  stop("Cannot determine B_sensitivities_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

OUT_DIR_B <- output_dir("B_Mobility")
RDS_PATH  <- file.path(OUT_DIR_B, "B_Sensitivities_Models.rds")
if (!file.exists(RDS_PATH)) stop("B_Sensitivities_Models.rds not found — run B_sensitivities_analysis.R first.")
mc <- readRDS(RDS_PATH)
cat(sprintf("B_sensitivities_report: loaded %s\n", RDS_PATH))

run_metadata <- build_run_metadata("B_Mobility (sensitivities)", mc$RUN_TS, SCRIPT_DIR)
nz <- function(x, fallback = data.frame(note = "not run / not available", stringsAsFactors = FALSE))
  if (is.null(x)) fallback else x

# ---- leave-one-study-out inputs -------------------------------------------------
# The LOO tables are computed in B_analysis.R and cached in B_Mobility_Models.rds, not in
# the sensitivities cache this script otherwise reads. They were exported as CSVs but
# never surfaced in a reader-facing workbook; sheet 05 below closes that gap.
.b_models_rds <- file.path(OUT_DIR_B, "B_Mobility_Models.rds")
mb <- if (file.exists(.b_models_rds)) readRDS(.b_models_rds) else NULL
if (is.null(mb))
  cat("B_sensitivities_report: NOTE — B_Mobility_Models.rds absent; LOO sheet degrades to notes.\n")

# ParEdu x region LOCO, written by B_paredu_moderation.R (run_B.command stage 2).
# The gradient is fitted directly on edu_z_kernel rather than derived from B1, so the
# file carries no B1_ prefix.
.paredu_loo_csv <- file.path(dirname(OUT_DIR_B), "manuscript_export", "mobility",
                             "sensitivities", "ParEdu_Gradient_LOO.csv")
ParEdu_LOO <- if (file.exists(.paredu_loo_csv))
  utils::read.csv(.paredu_loo_csv, stringsAsFactors = FALSE) else NULL
if (is.null(ParEdu_LOO))
  cat("B_sensitivities_report: NOTE — ParEdu_Gradient_LOO.csv absent.\n")
# The per-row `spec` column repeats the block note on every fold and appears in no
# other sheet in the file; drop it so the block reads like its neighbours.
.paredu_loo_tbl <- nz(ParEdu_LOO)
.paredu_loo_tbl <- .paredu_loo_tbl[, setdiff(names(.paredu_loo_tbl), "spec"), drop = FALSE]

DATA_NO <- 6L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)

readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Contents", "Samples", "Significance stars"),
  value = c("6 - Educational mobility, sensitivity analyses",
            "Robustness checks for the Section B (educational mobility) results",
            "Alternative PGIs (S1-B), height negative control (S2-B), heteroscedasticity (S3-B), leave-one-study-out",
            "Dedup one-adult-per-family mobility sample (S1/S2), plus the full-family mobility sample for the location-scale (gaulss) heteroscedasticity model",
            "* p<.05, ** p<.01, *** p<.001"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_Alternative_PGIs", "02_Height_control", "03_Heteroscedasticity", "04_FindingMatched_AltPGIs",
            "05_Leave_one_study_out"),
  analysis_id = c("S1-B", "S2-B", "S3-B", "S1-B (finding-matched)", "LOO"),
  description = c("Cog / NonCog PGIs on the B1 and B4 focals (Benjamini-Hochberg per family)",
                  "Height PGI placebo on mobility (B1 spec) and on height (sanity)",
                  "Genetic/environmental heteroscedasticity (gaulss), per region, plus the B0 East-West contrast",
                  "Cog / NonCog / Edu PGIs on the ACTUAL focal findings (B0 East-West contrast, pooled PGI x gender, B5 nonlinear), family-RE, common re-residualised sample, BH per hypothesis",
                  paste("Every Section B focal estimate refitted with one contributing study held out",
                        "(BASE-II / SOEP / TwinLife), including the regional parental-education",
                        "gradient. Distinct from Supplementary Data 7, which fits WITHIN each cohort")),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Sensitivity analyses for educational mobility. Parental education is controlled",
                    "throughout. S1-B and S2-B use the one-adult-per-family mobility sample; S3-B is gated",
                    "on the Section B focal results being significant and is fitted on the full-family",
                    "mobility sample without the family random effect (gaulss limitation).")),
  list(header = "Provenance", df = run_metadata, note = "Publication-safe run metadata."),
  list(header = "Sheet index", df = readme_index, note = "Internal IDs (S1-B, S2-B, S3-B) are kept inside each sheet."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

sheet_S1 <- multiblock_sheet(
  list(header = "Alternative-PGI family", df = nz(mc$S1_AltPGIs),
       note = paste("Re-runs the B1 (PGI x birth year x region) and B4 (four-way) focal tests substituting",
                    "the cognitive and non-cognitive PGIs for the education PGI, with Benjamini-Hochberg",
                    "correction within each family. Parental education is controlled. Tests whether the",
                    "focal mobility pattern is education-PGI-specific.")),
  title = .t(1, "Alternative-PGI sensitivity for educational mobility"))

sheet_S1fm <- multiblock_sheet(
  list(header = "Method / multiple-testing families", df = nz(mc$S1_FindingMatched_Status),
       note = paste("Finding-matched extension of S1: PGI-Education / PGI-Cognition / PGI-Noncognitive",
                    "substituted into the ACTUAL focal mobility findings, refit with the primary family",
                    "random-effects specification (parental education controlled) on one common",
                    "complete-case sample - each PGI re-residualised on ancestry PCs and re-z-scored within",
                    "that sample. Benjamini-Hochberg within each focal hypothesis across the three PGIs.",
                    "Differences are in PGI predictive association, not genetic effects.",
                    "See Plan_deviations.md section 15.")),
  list(header = "Finding 2 - East-West mobility contrast (PGI x region)", df = nz(mc$S1_finding2),
       note = "The highest-priority finding-matched test: whether the East-West difference in the PGI-mobility slope holds under the alternative PGIs."),
  list(header = "Finding 2 - East / West standardised PGI-mobility slopes", df = nz(mc$S1_finding2_slopes)),
  list(header = "Finding 3 - mobility gender (pooled PGI x gender)", df = nz(mc$S1_finding3_mob),
       note = "Pooled PGI x gender interaction (not only the four-way term) under each PGI."),
  list(header = "Finding 3 - mobility gender: per-cell PGI slopes (female / male)",
       df = nz(mc$S1_finding3_mob_slopes)),
  list(header = "B5 nonlinear gender reproduction (omnibus + per-gender curvature)", df = nz(mc$S1_B5_nonlin),
       note = "Whether the male-driven nonlinear cohort pattern reproduces under Cognition / Noncognitive; BH on the omnibus across the three PGIs."),
  list(header = "B5 female-male difference smooth: excludes-zero intervals per PGI", df = nz(mc$S1_B5_nonlin_diff_intervals)),
  title = .t(4, "Finding-matched alternative-PGI sensitivity for educational mobility"))

sheet_S2 <- multiblock_sheet(
  list(header = "Height PGI on mobility (placebo; expect null)", df = nz(mc$S2_Height_Mobility),
       note = "Negative control: the height PGI substituted into the B1 specification on mobility (parental education controlled) - should be null if the focal pattern is not a generic genetic artefact."),
  list(header = "Height PGI on height (sanity; expect strong)", df = nz(mc$S2_Height_Height),
       note = "Sanity check that the height PGI does predict measured height in the mobility sample."),
  title = .t(2, "Height negative-control sensitivity for educational mobility"))

B3 <- mc$B_S3
sheet_S3 <- multiblock_sheet(
  list(header = "Gate", df = nz(mc$S3_Status),
       note = "S3-B runs only if a Section B focal test is significant (read from the Section B model cache)."),
  list(header = "Focal dispersion robustness (significant focals only)", df = nz(mc$Focal_Dispersion),
       note = paste("Per-focal heteroscedasticity check, GATED on the main random-effects result:",
                    "for each focal significant in the main analysis, the mean effect is re-estimated",
                    "allowing the residual variance to follow the same structure (Gaussian location-scale,",
                    "no family RE; parental education controlled). Tested = FALSE rows were not significant",
                    "in the main analysis and are not run. A focal is robust if it stays significant under",
                    "heteroscedasticity (Het p); Scale LRT p tests whether the dispersion structure improves",
                    "fit. The B0 Region x PGI focal is shown here and, with the pre-reunification split, below.")),
  list(header = "Dispersion structure (per region)", df = nz(if (!is.null(B3)) B3$dispersion),
       note = paste("Scale-model coefficients from Gaussian location-scale (gaulss) fits per region:",
                    "environmental (variance ~ birth year), genetic (variance ~ PGI) and a combined model.",
                    "Parental education enters both the mean and is smoothed over birth year.")),
  list(header = "Focal cohort smooth: homoscedastic vs location-scale", df = nz(if (!is.null(B3)) B3$focal_under_het),
       note = "The focal PGI cohort smooth re-estimated allowing heteroscedasticity, to check the focal mobility result is not a variance artefact."),
  list(header = "B0 East-West contrast under heteroscedasticity", df = nz(if (!is.null(B3)) B3$b0_contrast),
       note = "The overall East-West PGI-mobility contrast (full and pre-reunification) under OLS vs the gaulss location-scale model."),
  list(header = "Region fit status", df = nz(if (!is.null(B3)) B3$region_status),
       note = paste("Per-region fit status. Note: S3-B is an approximation - the gaulss model carries no",
                    "family random effect and pools the cohort smooth across cells.")),
  title = .t(3, "Heteroscedasticity sensitivity for educational mobility"))

# ---- 05: leave-one-study-out ----
# Each focal estimate refitted with one contributing study held out. Folds are BASE-II,
# SOEP and TwinLife; SHIP carries no parental education and is absent from the mobility
# sample entirely, so it is not a fold. Leaving a study OUT is a different question from
# fitting WITHIN a study — the latter is Supplementary Data 7 (per-cohort appendix).
.loo_note <- paste("Focal estimate refitted with the named study held out. 'Full sample' is the",
  "published estimate; delta_vs_full is the fold estimate minus it. Folds: BASE-II, SOEP,",
  "TwinLife (SHIP is absent from the mobility sample — no parental education).")
sheet_LOO <- multiblock_sheet(
  list(header = "(a) B0 - overall East-West PGI-mobility contrast", df = nz(mb$B0_LOO), note = .loo_note),
  list(header = "(b) B0 - per-fold region-specific PGI-mobility slopes", df = nz(mb$B0_LOO_slopes),
       note = "West and East PGI-mobility slopes per held-out study, for the overall (B0) contrast."),
  list(header = "(c) B1 - PGI x birth year x region", df = nz(mb$B1_LOO), note = .loo_note),
  list(header = "(d) B3 - nonlinearity of the PGI-mobility cohort trajectory", df = nz(mb$B3_LOO),
       note = "Omnibus linear-vs-nonlinear LRT per held-out study."),
  list(header = "(e) B4 - PGI x birth year x region x gender", df = nz(mb$B4_LOO), note = .loo_note),
  list(header = "(f) B4 - pooled PGI x gender", df = nz(mb$B4_LOO_gender),
       note = "Pooled PGI x gender term (region control preserved) per held-out study."),
  list(header = "(g) B5 - per-gender nonlinearity", df = nz(mb$B5_LOO),
       note = "Omnibus LRT plus per-gender curvature per held-out study, regions pooled."),
  list(header = "(h) B6 - region x gender cohort trajectory", df = nz(mb$B6_LOO), note = .loo_note),
  list(header = "(i) B0g - overall region x gender contrast, per cell", df = nz(mb$B0g_LOO_cells),
       note = "Per region x gender cell, per held-out study."),
  list(header = "(j) ParEdu gradient - parental-education x region contrast", df = .paredu_loo_tbl,
       note = paste("Leave-one-study-out on the regional social-origin gradient reported in",
                    "Supplementary Data 5, sheet 09. Estimates are on the educational-attainment",
                    "metric, adjusted for birth year and gender. gradient_East / gradient_West",
                    "are the region-specific gradients per fold; 'estimate' is the",
                    "East-minus-West difference.")),
  title = .t(5, "Leave-one-study-out sensitivity of the Section B focal estimates"))

sheets <- list(
  "00_Index"             = sheet_README,
  "01_Alternative_PGIs"  = sheet_S1,
  "02_Height_control"    = sheet_S2,
  "03_Heteroscedasticity"= sheet_S3,
  "04_FindingMatched_AltPGIs" = sheet_S1fm,
  "05_Leave_one_study_out"    = sheet_LOO)

sheets <- humanize_sheets(sheets)
OUT_XLSX <- file.path(OUT_DIR_B, "Supplementary_Data_6_Mobility_sensitivities_results.xlsx")
write_section_xlsx("B_sensitivities", sheets, OUT_XLSX)
cat(sprintf("B_sensitivities_report: xlsx written (%d sheets)\n", length(sheets)))

# ---- Plots: S1 alt-PGI forest ----
OUT_PDF <- file.path(OUT_DIR_B, "B_Sensitivities_Plots.pdf")
plots <- list()
if (!is.null(mc$S1_AltPGIs) && nrow(mc$S1_AltPGIs)) {
  s1 <- mc$S1_AltPGIs
  plots[["S1"]] <- tryCatch(
    ggplot(s1, aes(x = estimate, y = PGI)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), height = 0.2) +
      geom_point(size = 2.5) + facet_wrap(~ family, scales = "free_x") +
      labs(x = "Focal estimate (95% CI)", y = "PGI",
           title = "S1-B: alternative-PGI focal estimates (mobility)") + theme_minimal(base_size = 12),
    error = function(e) NULL)
}
plots <- Filter(Negate(is.null), plots)
if (length(plots)) {
  pdf(OUT_PDF, width = 10, height = 6)
  for (p in plots) tryCatch(print(p), error = function(e) cat("plot failed:", conditionMessage(e), "\n"))
  dev.off()
  cat(sprintf("B_sensitivities_report: pdf written (%d pages)\n", length(plots)))
}

# ---- Overview markdown ----
analyses_md <- data.frame(
  Sensitivity = c("S1-B", "S2-B", "S3-B"),
  What = c("Alternative PGIs (Cog / NonCog) on B1 + B4 focals (BH per family)",
           "Height PGI negative control on mobility (B1 spec) + height sanity",
           "Heteroscedasticity (gaulss) per region + B0 contrast, gated on B focal significance"),
  Sample = c("dedup-OLS (dat_mob)", "dedup-OLS (dat_mob)", "dat_cluster_mob (no RE)"),
  stringsAsFactors = FALSE)
outputs_md <- data.frame(
  Output = c("Workbook", "Figures", "Model cache"),
  File = c("Supplementary_Data_6_Mobility_sensitivities_results.xlsx", "B_Sensitivities_Plots.pdf", "B_Sensitivities_Models.rds"),
  stringsAsFactors = FALSE)
generate_section_overview_md(
  section_id = "B_sensitivities",
  section_title = "Section B — Sensitivities (mobility)",
  purpose = paste("Mobility-phenotype sensitivity checks owned by Section B: alternative PGIs (S1-B),",
                  "height negative control on mobility (S2-B, new), and heteroscedasticity (S3-B) with",
                  "the B0 East-West contrast. Parental education is controlled throughout."),
  analyses = analyses_md,
  findings = nz(mc$S1_AltPGIs),
  interpretation = paste("S1-B and S2-B use the dedup-OLS mobility sample (relative comparison across PGIs",
                         "unaffected by the estimator); S3-B is gated on the Section B focal significance",
                         "read from B_Mobility_Models.rds. See Plan_deviations.md §9/§11/§12."),
  outputs = outputs_md,
  path = file.path(OUT_DIR_B, "B_Sensitivities_Overview.md"))
cat("B_sensitivities_report: DONE.\n")
