# B_appendix_report.R
# Section B per-study appendix — report stage. Reads B_Appendix_Models.rds and
# builds Supplementary_Data_7_Mobility_appendix_results.xlsx,
# B_Appendix_Plots.pdf, B_Appendix_Overview.md.

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
  stop("Cannot determine B_appendix_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

OUT_DIR_B <- output_dir("B_Mobility")
RDS_PATH  <- file.path(OUT_DIR_B, "B_Appendix_Models.rds")
if (!file.exists(RDS_PATH)) stop("B_Appendix_Models.rds not found — run B_appendix_analysis.R first.")
mc <- readRDS(RDS_PATH)
cat(sprintf("B_appendix_report: loaded %s\n", RDS_PATH))

run_metadata <- build_run_metadata("B_Mobility (appendix)", mc$RUN_TS, SCRIPT_DIR)
nz <- function(x) if (is.null(x)) data.frame(note = "not available", stringsAsFactors = FALSE) else x

DATA_NO <- 7L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)

readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Contents", "Sample", "Status"),
  value = c("7 - Educational mobility, per-study appendix",
            "Per-study (within-dataset) mobility analyses, appendix-grade",
            "Per-study region GAM (E), per-study region x gender GAM (E-RxG), per-study linear gender model (F-linear)",
            "One-adult-per-family mobility sample, fit separately within each study (BASE-II / SOEP / TwinLife; SHIP absent)",
            "Appendix-grade: per-study birth-year windows are narrow; cells with too few observations are skipped"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_PerStudy_region", "02_PerStudy_region_gender", "03_PerStudy_gender_linear"),
  analysis_id = c("E", "E-RxG", "F-linear"),
  description = c("Per-study region-varying PGI mobility cohort GAM (ParEdu-controlled)",
                  "Per-study region x gender PGI mobility cohort GAM (ParEdu-controlled)",
                  "Per-study robust linear gender model on mobility (PGI main / gender x PGI / region x gender x PGI)"),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Single-study complements to the pooled Section B analyses (parental education",
                    "controlled). The better-powered heterogeneity check is the leave-one-study-out",
                    "sensitivity in the main Section B workbook; these per-study fits are appendix-grade",
                    "because the within-study birth-year windows are narrow.")),
  list(header = "Provenance", df = run_metadata, note = "Publication-safe run metadata."),
  list(header = "Sheet index", df = readme_index, note = "Internal IDs (E, E-RxG, F-linear) are kept inside each sheet."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

sheet_E <- multiblock_sheet(
  list(header = "Per-study fit summary", df = nz(mc$E$table),
       note = paste("Within each study, a region-varying PGI mobility cohort GAM is fit (parental",
                    "education controlled; two regions where both have enough observations, otherwise",
                    "single-region). Columns report per-region smooth p-values and the linear-vs-nonlinear LRT.")),
  list(header = "Per-study smooth tables", df = nz(mc$E$smooths),
       note = "summary.gam smooth statistics for each study's region-specific PGI mobility cohort smooths."),
  title = .t(1, "Per-study region-varying results for educational mobility"))

sheet_F <- multiblock_sheet(
  list(header = "Per-study fit summary", df = nz(mc$E_rxg$table),
       note = paste("Within each study, a region x gender PGI mobility cohort GAM is fit (parental",
                    "education controlled; four cells where each has enough observations, otherwise a",
                    "single-region two-cell fit, otherwise skipped).")),
  list(header = "Per-study smooth tables", df = nz(mc$E_rxg$smooths),
       note = "summary.gam smooth statistics for each study's region x gender PGI mobility cohort smooths."),
  title = .t(2, "Per-study region x gender results for educational mobility"))

sheet_Fl <- multiblock_sheet(
  list(header = "Per-study linear gender model", df = nz(mc$F_linear_mob),
       note = paste("Robust per-study linear model (mobility ~ PGI x region x gender + birth year +",
                    "parental education), reporting the PGI main effect, the gender x PGI interaction",
                    "(men vs women) and the region x gender x PGI interaction where both regions are",
                    "present. Directly comparable across studies; * marks p < .05.")),
  title = .t(3, "Per-study linear gender model for educational mobility"))

sheets <- list("00_Index" = sheet_README, "01_PerStudy_region" = sheet_E,
               "02_PerStudy_region_gender" = sheet_F, "03_PerStudy_gender_linear" = sheet_Fl)
sheets <- humanize_sheets(sheets)
OUT_XLSX <- file.path(OUT_DIR_B, "Supplementary_Data_7_Mobility_appendix_results.xlsx")
write_section_xlsx("B_appendix", sheets, OUT_XLSX)
cat(sprintf("B_appendix_report: xlsx written (%d sheets)\n", length(sheets)))

# ---- Plots: per-study slope curves ----
OUT_PDF <- file.path(OUT_DIR_B, "B_Appendix_Plots.pdf")
plots <- list()
if (!is.null(mc$E$curves) && nrow(mc$E$curves)) {
  plots[["E"]] <- tryCatch(
    ggplot(mc$E$curves, aes(x = birth_year, y = slope, colour = east_west, fill = east_west)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.9) + facet_wrap(~ dataset, scales = "free_x") +
      scale_colour_manual(values = c(East = COL_EAST, West = COL_WEST)) +
      scale_fill_manual(  values = c(East = COL_EAST, West = COL_WEST)) +
      labs(x = "Birth year", y = "PGI->mobility slope",
           title = "E: per-study cohort-varying PGI-mobility slope by region", colour = "Region", fill = "Region") +
      theme_minimal(base_size = 12),
    error = function(e) NULL)
}
if (!is.null(mc$E_rxg$curves) && nrow(mc$E_rxg$curves)) {
  fc <- mc$E_rxg$curves; fc$group <- paste(fc$east_west, fc$gender)
  plots[["F"]] <- tryCatch(
    ggplot(fc, aes(x = birth_year, y = slope, colour = group, fill = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.8) + facet_wrap(~ dataset, scales = "free_x") +
      labs(x = "Birth year", y = "PGI->mobility slope",
           title = "E-RxG: per-study cohort-varying PGI-mobility slope by Region x Gender", colour = "Cell", fill = "Cell") +
      theme_minimal(base_size = 12),
    error = function(e) NULL)
}
plots <- Filter(Negate(is.null), plots)
if (length(plots)) {
  pdf(OUT_PDF, width = 11, height = 7)
  for (p in plots) tryCatch(print(p), error = function(e) cat("plot failed:", conditionMessage(e), "\n"))
  dev.off()
  cat(sprintf("B_appendix_report: pdf written (%d pages)\n", length(plots)))
}

# ---- Overview ----
analyses_md <- data.frame(
  Analysis = c("E", "E-RxG", "F-linear"),
  What = c("Per-study B3-style mobility GAM (region-specific PGI cohort smooths)",
           "Per-study B6-style Region x Gender GAM on mobility (ParEdu-controlled)",
           "Per-study robust linear gender model on mobility (PGI main / gender x PGI / Region x Gender x PGI)"),
  Sample = rep("dedup-OLS (dat_mob), within each study", 3),
  stringsAsFactors = FALSE)
outputs_md <- data.frame(
  Output = c("Workbook", "Figures", "Model cache"),
  File = c("Supplementary_Data_7_Mobility_appendix_results.xlsx", "B_Appendix_Plots.pdf", "B_Appendix_Models.rds"),
  stringsAsFactors = FALSE)
generate_section_overview_md(
  section_id = "B_appendix",
  section_title = "Section B — Per-study appendix (mobility)",
  purpose = paste("Appendix-grade per-study mobility analyses (E, E-RxG, F-linear), within each study",
                  "(BASE-II / SOEP / TwinLife; SHIP absent). Parental education controlled. Narrow",
                  "per-study birth-year windows mean several study x cell combos are skipped;",
                  "leave-one-study-out in the main Section B analyses remains the better-powered check."),
  analyses = analyses_md, findings = nz(mc$F_linear_mob),
  interpretation = paste("Single-study complements to the pooled Section B analyses; see",
                         "Plan_deviations.md sections 9 and 13."),
  outputs = outputs_md, path = file.path(OUT_DIR_B, "B_Appendix_Overview.md"))
cat("B_appendix_report: DONE.\n")
