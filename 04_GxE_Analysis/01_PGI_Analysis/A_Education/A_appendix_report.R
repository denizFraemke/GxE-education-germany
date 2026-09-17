# A_appendix_report.R
# Section A per-study appendix — report stage. Reads A_Appendix_Models.rds and
# builds Supplementary_Data_4_Education_appendix_results.xlsx,
# A_Appendix_Plots.pdf and A_Appendix_Overview.md.

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
  stop("Cannot determine A_appendix_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

OUT_DIR_A <- output_dir("A_Education")
RDS_PATH  <- file.path(OUT_DIR_A, "A_Appendix_Models.rds")
if (!file.exists(RDS_PATH)) stop("A_Appendix_Models.rds not found — run A_appendix_analysis.R first.")
mc <- readRDS(RDS_PATH)
cat(sprintf("A_appendix_report: loaded %s\n", RDS_PATH))

run_metadata <- build_run_metadata("A_Education (appendix)", mc$RUN_TS, SCRIPT_DIR)
nz <- function(x) if (is.null(x)) data.frame(note = "not available", stringsAsFactors = FALSE) else x

DATA_NO <- 4L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)

readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Contents", "Sample", "Status"),
  value = c("4 - Educational attainment, per-study appendix",
            "Per-study (within-dataset) education analyses, appendix-grade",
            "Per-study region GAM (D), per-study region x gender GAM (F), per-study linear gender model (F-linear)",
            "One-adult-per-family sample, fit separately within each study (BASE-II / SHIP / SOEP / TwinLife)",
            "Appendix-grade: per-study birth-year windows are narrow; cells with too few observations are skipped"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_PerStudy_region", "02_PerStudy_region_gender", "03_PerStudy_gender_linear"),
  analysis_id = c("D", "F", "F-linear"),
  description = c("Per-study region-varying PGI cohort GAM",
                  "Per-study region x gender PGI cohort GAM",
                  "Per-study robust linear gender model (PGI main / gender x PGI / region x gender x PGI)"),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Single-study complements to the pooled Section A analyses. The better-powered",
                    "heterogeneity check is the leave-one-study-out sensitivity in the main Section A",
                    "workbook; these per-study fits are appendix-grade because the within-study",
                    "birth-year windows are narrow.")),
  list(header = "Provenance", df = run_metadata, note = "Publication-safe run metadata."),
  list(header = "Sheet index", df = readme_index, note = "Internal IDs (D, F, F-linear) are kept inside each sheet."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

sheet_D <- multiblock_sheet(
  list(header = "Per-study fit summary", df = nz(mc$D$table),
       note = paste("Within each study, a region-varying PGI cohort GAM is fit (two regions where both",
                    "have enough observations, otherwise single-region). Columns report per-region smooth",
                    "p-values and the linear-vs-nonlinear LRT; 'fit' notes which specification was used.")),
  list(header = "Per-study smooth tables", df = nz(mc$D$smooths),
       note = "summary.gam smooth statistics for each study's region-specific PGI cohort smooths."),
  title = .t(1, "Per-study region-varying results for educational attainment"))

sheet_F <- multiblock_sheet(
  list(header = "Per-study fit summary", df = nz(mc$F$table),
       note = paste("Within each study, a region x gender PGI cohort GAM is fit (four cells where each",
                    "has enough observations, otherwise a single-region two-cell fit, otherwise skipped).")),
  list(header = "Per-study smooth tables", df = nz(mc$F$smooths),
       note = "summary.gam smooth statistics for each study's region x gender PGI cohort smooths."),
  title = .t(2, "Per-study region x gender results for educational attainment"))

sheet_Fl <- multiblock_sheet(
  list(header = "Per-study linear gender model", df = nz(mc$F_linear_edu),
       note = paste("Robust per-study linear model (education ~ PGI x region x gender + birth year),",
                    "reporting the PGI main effect, the gender x PGI interaction (men vs women) and the",
                    "region x gender x PGI interaction where both regions are present. Directly comparable",
                    "across studies; * marks p < .05.")),
  title = .t(3, "Per-study linear gender model for educational attainment"))

sheets <- list("00_Index" = sheet_README, "01_PerStudy_region" = sheet_D,
               "02_PerStudy_region_gender" = sheet_F, "03_PerStudy_gender_linear" = sheet_Fl)
sheets <- humanize_sheets(sheets)
OUT_XLSX <- file.path(OUT_DIR_A, "Supplementary_Data_4_Education_appendix_results.xlsx")
write_section_xlsx("A_appendix", sheets, OUT_XLSX)
cat(sprintf("A_appendix_report: xlsx written (%d sheets)\n", length(sheets)))

# ---- Plots: per-study slope curves ----
OUT_PDF <- file.path(OUT_DIR_A, "A_Appendix_Plots.pdf")
plots <- list()
if (!is.null(mc$D$curves) && nrow(mc$D$curves)) {
  plots[["D"]] <- tryCatch(
    ggplot(mc$D$curves, aes(x = birth_year, y = slope, colour = east_west, fill = east_west)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.9) + facet_wrap(~ dataset, scales = "free_x") +
      scale_colour_manual(values = c(East = COL_EAST, West = COL_WEST)) +
      scale_fill_manual(  values = c(East = COL_EAST, West = COL_WEST)) +
      labs(x = "Birth year", y = "PGI->education slope",
           title = "D: per-study cohort-varying PGI slope by region", colour = "Region", fill = "Region") +
      theme_minimal(base_size = 12),
    error = function(e) NULL)
}
if (!is.null(mc$F$curves) && nrow(mc$F$curves)) {
  fc <- mc$F$curves; fc$group <- paste(fc$east_west, fc$gender)
  plots[["F"]] <- tryCatch(
    ggplot(fc, aes(x = birth_year, y = slope, colour = group, fill = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.8) + facet_wrap(~ dataset, scales = "free_x") +
      labs(x = "Birth year", y = "PGI->education slope",
           title = "F: per-study cohort-varying PGI slope by Region x Gender", colour = "Cell", fill = "Cell") +
      theme_minimal(base_size = 12),
    error = function(e) NULL)
}
plots <- Filter(Negate(is.null), plots)
if (length(plots)) {
  pdf(OUT_PDF, width = 11, height = 7)
  for (p in plots) tryCatch(print(p), error = function(e) cat("plot failed:", conditionMessage(e), "\n"))
  dev.off()
  cat(sprintf("A_appendix_report: pdf written (%d pages)\n", length(plots)))
}

# ---- Overview ----
analyses_md <- data.frame(
  Analysis = c("D", "F", "F-linear"),
  What = c("Per-study A3-style education GAM (region-specific PGI cohort smooths)",
           "Per-study A6-style Region x Gender GAM on edu_z_kernel",
           "Per-study robust linear gender model (PGI main / gender x PGI / Region x Gender x PGI)"),
  Sample = rep("dedup-OLS (dat_edu), within each study", 3),
  stringsAsFactors = FALSE)
outputs_md <- data.frame(
  Output = c("Workbook", "Figures", "Model cache"),
  File = c("Supplementary_Data_4_Education_appendix_results.xlsx", "A_Appendix_Plots.pdf", "A_Appendix_Models.rds"),
  stringsAsFactors = FALSE)
generate_section_overview_md(
  section_id = "A_appendix",
  section_title = "Section A — Per-study appendix (attainment)",
  purpose = paste("Appendix-grade per-study education analyses (D, F, F-linear), within each study",
                  "(BASE-II / SHIP / SOEP / TwinLife). Narrow per-study birth-year windows mean several",
                  "study x cell combos are skipped; leave-one-study-out in the main Section A analyses",
                  "remains the better-powered heterogeneity check."),
  analyses = analyses_md, findings = nz(mc$F_linear_edu),
  interpretation = paste("Single-study complements to the pooled Section A analyses; see",
                         "Plan_deviations.md section 13. The mobility per-study analysis (E) lives in Section B."),
  outputs = outputs_md, path = file.path(OUT_DIR_A, "A_Appendix_Overview.md"))
cat("A_appendix_report: DONE.\n")
