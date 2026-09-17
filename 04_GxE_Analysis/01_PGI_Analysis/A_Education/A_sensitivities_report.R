# A_sensitivities_report.R
# Section A sensitivities — report stage. Reads A_Sensitivities_Models.rds and
# builds Supplementary_Data_3_Education_sensitivities_results.xlsx,
# A_Sensitivities_Plots.pdf and A_Sensitivities_Overview.md.

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
  stop("Cannot determine A_sensitivities_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

OUT_DIR_A <- output_dir("A_Education")
RDS_PATH  <- file.path(OUT_DIR_A, "A_Sensitivities_Models.rds")
if (!file.exists(RDS_PATH)) stop("A_Sensitivities_Models.rds not found — run A_sensitivities_analysis.R first.")
mc <- readRDS(RDS_PATH)
cat(sprintf("A_sensitivities_report: loaded %s\n", RDS_PATH))

run_metadata <- build_run_metadata("A_Education (sensitivities)", mc$RUN_TS, SCRIPT_DIR)
nz <- function(x, fallback = data.frame(note = "not run / not available", stringsAsFactors = FALSE))
  if (is.null(x)) fallback else x

# ---- leave-one-study-out inputs -------------------------------------------------
# Computed in A_analysis.R and cached in A_Education_Models.rds, not in the
# sensitivities cache this script otherwise reads.
.a_models_rds <- file.path(OUT_DIR_A, "A_Education_Models.rds")
ma <- if (file.exists(.a_models_rds)) readRDS(.a_models_rds) else NULL
if (is.null(ma))
  cat("A_sensitivities_report: NOTE — A_Education_Models.rds absent; LOO sheet degrades to notes.\n")

DATA_NO <- 3L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)

readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Contents", "Samples", "Significance stars"),
  value = c("3 - Educational attainment, sensitivity analyses",
            "Robustness checks for the Section A (educational attainment) results",
            "Alternative PGIs (S1), height negative control (S2), migration diagnostics (S4), estimator comparison, heteroscedasticity (S3), leave-one-study-out",
            "Dedup one-adult-per-family sample (S1/S2/S4), plus the full-family sample for CR2/RE and the location-scale (gaulss) heteroscedasticity model",
            "* p<.05, ** p<.01, *** p<.001"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_Alternative_PGIs", "02_Height_control", "03_Migration", "04_Estimator_comparison", "05_Heteroscedasticity", "06_FindingMatched_AltPGIs",
            "07_Leave_one_study_out"),
  analysis_id = c("S1-A", "S2-A", "S4", "Method", "S3-A", "S1-A (finding-matched)", "LOO"),
  description = c("Cog / NonCog PGIs on the A1 and A4 focals (Benjamini-Hochberg per family)",
                  "Height PGI placebo on education (A1 spec) and on height (sanity)",
                  "Migration-composition diagnostics, gated on A focal significance",
                  "A1 focal under dedup-OLS / cluster-robust (CR2) / random-effects estimators",
                  "Genetic/environmental heteroscedasticity (Gaussian location-scale model), per region",
                  "Cog / NonCog / Edu PGIs on the ACTUAL focal findings (pooled PGI x BY, pooled PGI x gender, per-gender GAM), family-RE, common re-residualised sample, BH per hypothesis",
                  paste("Every Section A focal estimate refitted with one contributing study held out.",
                        "Distinct from the per-study appendix, which fits WITHIN each study")),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Sensitivity analyses for educational attainment. Most are estimated on the",
                    "one-adult-per-family sample (relative comparison across PGIs is unaffected by the",
                    "estimator); S4 and S3 are gated on the Section A focal results being significant.")),
  list(header = "Provenance", df = run_metadata, note = "Publication-safe run metadata."),
  list(header = "Sheet index", df = readme_index, note = "Internal IDs (S1-A, S2-A, S4, Method, S3-A) are kept inside each sheet."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

sheet_S1 <- multiblock_sheet(
  list(header = "Alternative-PGI family", df = nz(mc$S1_AltPGIs),
       note = paste("Re-runs the A1 (PGI x birth year x region) and A4 (four-way) focal tests substituting",
                    "the cognitive and non-cognitive PGIs for the education PGI, with Benjamini-Hochberg",
                    "correction within each family. Tests whether the focal pattern is education-specific.")),
  title = .t(1, "Alternative-PGI sensitivity for educational attainment"))

sheet_S1fm <- multiblock_sheet(
  list(header = "Method / multiple-testing families", df = nz(mc$S1_FindingMatched_Status),
       note = paste("Finding-matched extension of S1: PGI-Education / PGI-Cognition / PGI-Noncognitive",
                    "substituted into the ACTUAL focal attainment models, refit with the primary family",
                    "random-effects specification on one common complete-case sample - each PGI",
                    "re-residualised on ancestry PCs and re-z-scored within that sample so all three are",
                    "on equal footing. Benjamini-Hochberg within each focal hypothesis across the three",
                    "PGIs. Differences are in PGI predictive association, not genetic effects.",
                    "See Plan_deviations.md section 15.")),
  list(header = "Finding 1 - attainment across cohorts (pooled PGI x birth year; + PGI x BY x region context)",
       df = nz(mc$S1_finding1),
       note = "Pooled PGI x birth-year interaction (the reported cohort increase) and the contextual three-way."),
  list(header = "Finding 3 - attainment gender (pooled PGI x gender)", df = nz(mc$S1_finding3_att),
       note = "Pooled PGI x gender interaction (not only the four-way term) under each PGI."),
  list(header = "Finding 3 - attainment gender: per-cell PGI slopes (female / male)",
       df = nz(mc$S1_finding3_att_slopes)),
  list(header = "Per-gender attainment GAM (A5-equivalent; omnibus + curvature, symmetry null-check)",
       df = nz(mc$S1_A5_nonlin),
       note = "Reported for symmetry with the mobility B5 reproduction; the attainment result is a stable null."),
  title = .t(6, "Finding-matched alternative-PGI sensitivity for educational attainment"))

sheet_S2 <- multiblock_sheet(
  list(header = "Height PGI on education (placebo; expect null)", df = nz(mc$S2_Height_Education),
       note = "Negative control: the height PGI substituted into the A1 specification on years of education - should be null if the focal pattern is not a generic genetic artefact."),
  list(header = "Height PGI on height (sanity; expect strong)", df = nz(mc$S2_Height_Height),
       note = "Sanity check that the height PGI does predict measured height."),
  title = .t(2, "Height negative-control sensitivity for educational attainment"))

sheet_S4 <- multiblock_sheet(
  list(header = "Outer gate", df = nz(mc$S4_Status),
       note = "S4 runs only if a Section A focal test is significant (read from the Section A model cache)."),
  list(header = "Migration-composition diagnostics", df = nz(mc$S4_Migration_Tests),
       note = paste("Three necessary conditions for selective East-West migration to drive the regional",
                    "PGI difference: lower East mean PGI (S4a), lower East PGI variance (S4b), and",
                    "practically equivalent education variance (S4c, TOST equivalence against a",
                    "variance-ratio bound - not merely a non-significant difference).")),
  list(header = "Inner gate decision", df = nz(mc$S4d_GatingDecision),
       note = "The within-region refit (S4d) runs only if all three conditions hold (S4c requires established equivalence, not absence of evidence)."),
  list(header = "Absolute vs within-region PGI", df = nz(mc$S4d_Compare),
       note = "If triggered, compares the focal estimate using the pooled PGI vs a PGI standardized within region x cohort (region-specific kernel-z over birth year); attenuation would suggest migration-driven composition."),
  list(header = "Within-region PGI coefficients", df = nz(mc$S4d_Coefficients),
       note = "Full coefficients from the within-region-PGI refit (if run)."),
  title = .t(3, "Migration-composition diagnostics for educational attainment"))

sheet_M <- multiblock_sheet(
  list(header = "Estimator comparison", df = nz(mc$Method_Comparison),
       note = paste("The A1 focal coefficient under three estimators: dedup one-adult-per-family OLS,",
                    "cluster-robust CR2 standard errors (clustered by family), and the random-effects",
                    "model. Agreement across estimators indicates the focal result is not a clustering artefact.")),
  list(header = "Status", df = nz(mc$Method_Status),
       note = "Availability of the CR2 package and the Section A model cache used for the RE row."),
  title = .t(4, "Estimator comparison for the educational-attainment focal test"))

A3 <- mc$A_S3
sheet_S3 <- multiblock_sheet(
  list(header = "Gate", df = nz(mc$S3_Status),
       note = "S3 runs only if a Section A focal test is significant."),
  list(header = "Focal dispersion robustness (significant focals only)", df = nz(mc$Focal_Dispersion),
       note = paste("Per-focal heteroscedasticity check, GATED on the main random-effects result:",
                    "for each focal significant in the main analysis, the mean effect is re-estimated",
                    "allowing the residual variance to follow the same structure (Gaussian location-scale,",
                    "no family RE). Tested = FALSE rows were not significant in the main analysis and are",
                    "not run. A focal is robust if it stays significant under heteroscedasticity (Het p);",
                    "Scale LRT p tests whether the dispersion structure itself improves fit.")),
  list(header = "Dispersion structure (per region)", df = nz(if (!is.null(A3)) A3$dispersion),
       note = paste("Scale-model coefficients from Gaussian location-scale (gaulss) fits per region:",
                    "environmental (variance ~ birth year), genetic (variance ~ PGI) and a combined model.")),
  list(header = "Focal cohort smooth: homoscedastic vs location-scale", df = nz(if (!is.null(A3)) A3$focal_under_het),
       note = "The focal PGI cohort smooth re-estimated allowing heteroscedasticity, to check the focal result is not a variance artefact."),
  list(header = "Region fit status", df = nz(if (!is.null(A3)) A3$region_status),
       note = paste("Per-region fit status. Note: S3-A is an approximation - the gaulss model carries no",
                    "family random effect and pools the cohort smooth across cells.")),
  title = .t(5, "Heteroscedasticity sensitivity for educational attainment"))

# ---- 07: leave-one-study-out ----
# Sits beside 04_Estimator_comparison: both ask whether a focal estimate survives a change
# in what is being fitted rather than a change in the predictor. Leaving a study OUT is a
# different question from fitting WITHIN a study — the latter is the per-study appendix.
.loo_note_a <- paste("Focal estimate refitted with the named study held out. 'Full sample' is the",
  "published estimate; delta_vs_full is the fold estimate minus it. Folds are the",
  "contributing studies present in the attainment sample (BASE-II, SHIP, SOEP, TwinLife).")
sheet_LOO_A <- multiblock_sheet(
  list(header = "(a) A1 - PGI x birth year x region", df = nz(ma$A1_LOO), note = .loo_note_a),
  list(header = "(b) A2 - PGI x reunification x region", df = nz(ma$A2_LOO), note = .loo_note_a),
  list(header = "(c) A3 - nonlinearity of the PGI-attainment cohort trajectory", df = nz(ma$A3_LOO),
       note = "Omnibus linear-vs-nonlinear LRT per held-out study."),
  list(header = "(d) A4 - PGI x birth year x region x gender", df = nz(ma$A4_LOO), note = .loo_note_a),
  list(header = "(e) A4 - pooled PGI x gender", df = nz(ma$A4_LOO_gender),
       note = "Pooled PGI x gender term (region control preserved) per held-out study."),
  list(header = "(f) A5 - per-gender nonlinearity", df = nz(ma$A5_LOO),
       note = "Omnibus LRT plus per-gender curvature per held-out study, regions pooled."),
  list(header = "(g) A6 - region x gender cohort trajectory", df = nz(ma$A6_LOO), note = .loo_note_a),
  list(header = "(h) A0g - overall region x gender contrast, per cell", df = nz(ma$A0g_LOO_cells),
       note = "Per region x gender cell, per held-out study."),
  title = .t(7, "Leave-one-study-out sensitivity of the Section A focal estimates"))

sheets <- list(
  "00_Index"               = sheet_README,
  "01_Alternative_PGIs"     = sheet_S1,
  "02_Height_control"       = sheet_S2,
  "03_Migration"            = sheet_S4,
  "04_Estimator_comparison" = sheet_M,
  "05_Heteroscedasticity"   = sheet_S3,
  "06_FindingMatched_AltPGIs" = sheet_S1fm,
  "07_Leave_one_study_out"    = sheet_LOO_A)

sheets <- humanize_sheets(sheets)
OUT_XLSX <- file.path(OUT_DIR_A, "Supplementary_Data_3_Education_sensitivities_results.xlsx")
write_section_xlsx("A_sensitivities", sheets, OUT_XLSX)
cat(sprintf("A_sensitivities_report: xlsx written (%d sheets)\n", length(sheets)))

# ---- Plots: S1 alt-PGI forest + Method-comparison forest ----
OUT_PDF <- file.path(OUT_DIR_A, "A_Sensitivities_Plots.pdf")
plots <- list()
if (!is.null(mc$S1_AltPGIs) && nrow(mc$S1_AltPGIs)) {
  s1 <- mc$S1_AltPGIs
  plots[["S1"]] <- tryCatch(
    ggplot(s1, aes(x = estimate, y = PGI)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), height = 0.2) +
      geom_point(size = 2.5) + facet_wrap(~ family, scales = "free_x") +
      labs(x = "Focal estimate (95% CI)", y = "PGI",
           title = "S1-A: alternative-PGI focal estimates") + theme_minimal(base_size = 12),
    error = function(e) NULL)
}
if (!is.null(mc$Method_Comparison) && nrow(mc$Method_Comparison)) {
  mcomp <- mc$Method_Comparison
  plots[["M"]] <- tryCatch(
    ggplot(mcomp, aes(x = estimate, y = estimator)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), height = 0.2) +
      geom_point(size = 2.5) +
      labs(x = "A1 focal estimate (95% CI)", y = NULL,
           title = "Method comparison: A1 focal under dedup-OLS / CR2 / RE") + theme_minimal(base_size = 12),
    error = function(e) NULL)
}
plots <- Filter(Negate(is.null), plots)
if (length(plots)) {
  pdf(OUT_PDF, width = 10, height = 6)
  for (p in plots) tryCatch(print(p), error = function(e) cat("plot failed:", conditionMessage(e), "\n"))
  dev.off()
  cat(sprintf("A_sensitivities_report: pdf written (%d pages)\n", length(plots)))
}

# ---- Overview markdown ----
analyses_md <- data.frame(
  Sensitivity = c("S1-A", "S2-A", "S4", "Method", "S3-A"),
  What = c("Alternative PGIs (Cog / NonCog) on A1 + A4 focals (BH per family)",
           "Height PGI negative control (education / height / four-way)",
           "Migration diagnostics (A-only), gated on A focal significance",
           "A1 focal under dedup-OLS / CR2 / RE",
           "Heteroscedasticity (gaulss) per region, gated on A focal significance"),
  Sample = c("dedup-OLS (dat_edu)", "dedup-OLS (dat_edu)", "dedup-OLS (dat_edu)",
             "dat_edu / dat_cluster", "dat_cluster (no RE)"),
  stringsAsFactors = FALSE)
outputs_md <- data.frame(
  Output = c("Workbook", "Figures", "Model cache"),
  File = c("Supplementary_Data_3_Education_sensitivities_results.xlsx", "A_Sensitivities_Plots.pdf", "A_Sensitivities_Models.rds"),
  stringsAsFactors = FALSE)
generate_section_overview_md(
  section_id = "A_sensitivities",
  section_title = "Section A — Sensitivities (attainment)",
  purpose = paste("Attainment-phenotype sensitivity checks owned by Section A: alternative PGIs (S1-A),",
                  "height negative control (S2-A), migration diagnostics (S4), estimator comparison",
                  "(Method), and heteroscedasticity (S3-A). All dedup-OLS except CR2/RE and the gaulss",
                  "S3-A (full-family, no RE)."),
  analyses = analyses_md,
  findings = nz(mc$S1_AltPGIs),
  interpretation = paste("S1-A and S2-A use the dedup-OLS sample (relative comparison across PGIs",
                         "unaffected by the estimator); S4 and S3-A are gated on A focal significance",
                         "read from A_Education_Models.rds. See Plan_deviations.md §11/§12."),
  outputs = outputs_md,
  path = file.path(OUT_DIR_A, "A_Sensitivities_Overview.md"))
cat("A_sensitivities_report: DONE.\n")
