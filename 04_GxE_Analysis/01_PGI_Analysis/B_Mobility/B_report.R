# B_report.R
# Build Supplementary_Data_5_Mobility_results.xlsx + B_Mobility_Plots.pdf +
# B_Mobility_LOO_Plots.pdf + B_Mobility_Overview.md from the B model cache
# produced by B_analysis.R. Gender is a main moderator (B4/B5/B6/B0g).
# Set REPORT_SKIP_PLOTS=1 to write the workbook and overview without the PDFs.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
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
  stop("Cannot determine B_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)

OUT_DIR_B <- output_dir("B_Mobility")
RDS_PATH  <- file.path(OUT_DIR_B, "B_Mobility_Models.rds")
if (!file.exists(RDS_PATH)) stop("B_Mobility_Models.rds not found — run B_analysis.R first.")
mc <- readRDS(RDS_PATH)
# Optional tables degrade to a "not available" note rather than crashing the
# workbook build.
nz <- function(x, fallback = data.frame(note = "not run / not available", stringsAsFactors = FALSE))
  if (is.null(x)) fallback else x
dcm <- mc$dat_cluster_mob
attr(dcm, "by_mean") <- mc$by_mean
m_B5_nonlin <- mc$m_B5_nonlin_re
m_B6_nonlin <- mc$m_B6_nonlin_re
cohorts_b <- sort(as.character(unique(droplevels(dcm$cohort))))
N_b <- nrow(dcm)
cat(sprintf("B_report: loaded models (dat_cluster_mob N = %d)\n", N_b))

run_metadata <- build_run_metadata("B_Mobility", mc$RUN_TS, SCRIPT_DIR)

# Regional social-origin gradient — fitted by B_paredu_moderation.R (run_B.command
# stage 2, which is ordered before this script so the CSV exists here). Read rather
# than recomputed: it comes from a dedicated model that is not in this cache. A
# missing file degrades the sheet to a note rather than failing the workbook.
.grad_csv <- file.path(dirname(OUT_DIR_B), "manuscript_export", "mobility", "main",
                       "ParEdu_Gradient.csv")
ParEdu_Gradient <- if (file.exists(.grad_csv))
  utils::read.csv(.grad_csv, stringsAsFactors = FALSE) else NULL
if (is.null(ParEdu_Gradient))
  cat("B_report: NOTE — ParEdu_Gradient.csv absent; sheet 09 degrades to a note.\n")

DATA_NO <- 5L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)
.note_smoothF <- paste("summary.gam F-tests for the region/gender-specific PGI x birth-year smooths.",
  "These overstate curvature because the smooth's linear part is confounded with the parametric",
  "PGI x birth-year slope; use the clean per-cell nested-LRT curvature test above.")

# Workflow-wide p convention (.fmt_p_vec: "<.001" or three decimals) so the
# overview markdown and the workbook render the same p identically.
.fmt_p <- function(p) {
  if (length(p) == 0) return("NA"); if (length(p) > 1) p <- p[1]
  if (is.na(p)) return("NA")
  .fmt_p_vec(p)
}
.headline_from_coef <- function(coef_tbl, term, analysis_id, model_family) {
  r <- if (is.null(coef_tbl)) NULL else coef_tbl[coef_tbl$term == term, ]
  if (is.null(r) || !nrow(r)) return(build_headline_table(analysis_id = analysis_id, outcome = "mobility",
    sample_name = "dat_cluster_mob (full-family RE)", N = N_b, cohorts_included = cohorts_b,
    model_family = model_family, focal_term = term, estimate = NA, SE_or_CI = NA, p = NA))
  build_headline_table(analysis_id = analysis_id, outcome = "mobility",
    sample_name = "dat_cluster_mob (full-family RE)", N = N_b, cohorts_included = cohorts_b,
    model_family = model_family, focal_term = term, estimate = round(r$estimate, 4),
    SE_or_CI = sprintf("SE=%.4f (95%% CI [%.4f, %.4f])", r$std.error, r$conf.low, r$conf.high),
    p = round(r$p.value, 4))
}
.headline_lrt <- function(nested, analysis_id, model_family) {
  r <- nested[which(nested$comparison == "linear vs nonlinear"), ]
  build_headline_table(analysis_id = analysis_id, outcome = "mobility",
    sample_name = "dat_cluster_mob (full-family RE)", N = N_b, cohorts_included = cohorts_b,
    model_family = model_family, focal_term = "omnibus nested LRT (linear vs nonlinear)",
    estimate = if (nrow(r)) sprintf("chi2 = %.2f, df = %.2f", r$LRT_chisq, r$LRT_df) else NA,
    SE_or_CI = "n/a (see GAM overview)", p = if (nrow(r)) round(r$LRT_p, 4) else NA)
}

B1_Headline  <- .headline_from_coef(mc$B1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c", "B1", "Linear mixed model (family RE)")
B2_Headline  <- .headline_from_coef(mc$B2_Coefficients, "PGI_Edu_z:reunifpost:east_west_c", "B2", "Step-function mixed model (family RE)")
B3_Headline  <- .headline_lrt(mc$B3_NestedLRT, "B3", "Varying-coefficient GAM (family RE)")
B0_Headline  <- .headline_from_coef(mc$B0_Coefficients, "PGI_Edu_z:east_west_c", "B0", "Overall East-West contrast (family RE)")
B0_pre_Headline <- if (!is.null(mc$B0_pre_Coefficients))
  .headline_from_coef(mc$B0_pre_Coefficients, "PGI_Edu_z:east_west_c", "B0_pre", "Overall East-West contrast, pre-reunification (family RE)") else NULL
B4_Headline  <- .headline_from_coef(mc$B4_Coefficients, mc$B4_focal, "B4", "Linear four-way mixed model (family RE)")
B5_Headline  <- .headline_lrt(mc$B5_NestedLRT, "B5", "Per-gender varying-coefficient GAM (regions pooled, family RE)")
B6_Headline  <- .headline_lrt(mc$B6_NestedLRT, "B6", "Region x Gender varying-coefficient GAM (family RE)")
B0g_Headline <- .headline_from_coef(mc$B0g$coefs, mc$B0g_focal, "B0g", "Overall Region x Gender contrast (family RE)")

B5_tp_pgi <- mc$B5_Smooths_tp[grepl("PGI", mc$B5_Smooths_tp$smooth_term), , drop = FALSE]
B3_GAM_Overview <- build_gam_overview_table(mc$B3_NestedLRT, mc$B3_Smooths_pgi_tp, "dat_cluster_mob", N_b, "B3", alpha = 0.05, cell_lrt = mc$B3_cell_LRT)
B5_GAM_Overview <- build_gam_overview_table(mc$B5_NestedLRT, B5_tp_pgi,            "dat_cluster_mob", N_b, "B5", alpha = 0.05, cell_lrt = mc$B5_cell_LRT)
B6_GAM_Overview <- build_gam_overview_table(mc$B6_NestedLRT, mc$B6_Smooths_pgi,    "dat_cluster_mob", N_b, "B6", alpha = 0.05, cell_lrt = mc$B6_cell_LRT)

.headline_to_chr <- function(df) { for (col in c("estimate","SE_or_CI","p","N")) if (col %in% names(df)) df[[col]] <- as.character(df[[col]]); df }
# ---- 00_Index ----
readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Outcome variable", "Analytic sample", "Estimator", "Significance stars"),
  value = c("5 - Educational mobility",
            "Educational mobility (gene-environment interaction with the education PGI)",
            "Mobility = own minus parental education, each standardized within birth year (kernel-z difference)",
            "", # placeholder replaced below
            "Linear and generalized additive mixed models (mgcv::bam, REML); parental education is controlled throughout",
            "* p<.05, ** p<.01, *** p<.001"),
  stringsAsFactors = FALSE)
readme_file$value[4] <- sprintf("Full-family random-effects mobility sample, N = %d (SHIP absent: no parental education)", N_b)
readme_index <- data.frame(
  sheet = c("01_Linear_birth_year", "02_Step_reunification", "03_GAM_region", "04_Overall_EastWest",
            "05_Linear_gender", "06_GAM_gender", "07_GAM_region_gender", "08_Overall_RxG",
            "09_ParEdu_gradient"),
  analysis_id = c("B1", "B2", "B3", "B0 / B0_pre", "B4", "B5", "B6", "B0g", "ParEdu gradient"),
  description = c("Linear PGI x birth-year x region interaction",
                  "Step-function PGI x reunification x region interaction",
                  "Region-varying nonlinear (GAM) PGI-mobility cohort trajectory",
                  "Overall (and pre-reunification) East-West PGI-mobility contrast",
                  "Four-way PGI x birth-year x region x gender interaction",
                  "Per-gender PGI-mobility cohort trajectory, regions pooled (GAM)",
                  "Region x gender PGI-mobility cohort trajectory (GAM)",
                  "Overall additive region x gender PGI-mobility contrast",
                  paste("Regional social-origin gradient in educational attainment",
                        "(pooled, East-West difference, East, West), adjusted for birth",
                        "year and gender")),
  model = c("Linear mixed model", "Step-function mixed model", "Varying-coefficient GAM", "Linear contrast",
            "Linear four-way mixed model", "Per-gender GAM", "Region x gender GAM", "Linear contrast",
            "Linear mixed model (dedicated fit)"),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Section B reports how the education PGI predicts intergenerational educational",
                    "mobility across birth cohorts, East/West Germany, and gender, controlling for",
                    "parental education throughout (Plan_deviations.md §9/§10). Region and gender are",
                    "contrast-coded (East = +0.5 / West = -0.5; Female = +0.5 / Male = -0.5); coefficients",
                    "are per 1 SD of the PGI or birth year. Gender (sheets 5-8) is a main moderator,",
                    "corresponding to Section C of the preregistration plus a per-gender GAM.")),
  list(header = "Provenance", df = run_metadata, note = "Publication-safe run metadata."),
  list(header = "Sheet index", df = readme_index, note = "Internal analysis IDs (B0-B6, B0g) are kept inside each sheet."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

# ---- Planned East-West difference smooth (B3) ----
diff_B <- diff_smooth_simul(mc$m_B3_nonlin_re, dcm)
B3_diff_intervals <- diff_smooth_significant_intervals(diff_B, contrast = "East - West PGI->mobility slope")

# (c) contrasts for the gender GAMs: B5 female-male (regions pooled);
# B6 female-male within each region.
B5_diff_intervals <- diff_smooth_significant_intervals(
  diff_smooth_pooled_gender_simul(m_B5_nonlin, dcm),
  contrast = "Female - Male PGI->mobility slope",
  dir_pos = "Female > Male", dir_neg = "Female < Male")
B6_diff_intervals <- rbind(
  diff_smooth_significant_intervals(
    diff_smooth_gender_simul(m_B6_nonlin, dcm, east_west = "West"),
    contrast = "Female - Male PGI->mobility slope (West)",
    dir_pos = "Female > Male", dir_neg = "Female < Male"),
  diff_smooth_significant_intervals(
    diff_smooth_gender_simul(m_B6_nonlin, dcm, east_west = "East"),
    contrast = "Female - Male PGI->mobility slope (East)",
    dir_pos = "Female > Male", dir_neg = "Female < Male"))

# ---- Per-analysis fit rows ----
fit_B1 <- mc$B_ModelFit[mc$B_ModelFit$model == "B1_RE", , drop = FALSE]
fit_B2 <- mc$B_ModelFit[mc$B_ModelFit$model == "B2_RE", , drop = FALSE]
fit_B3 <- mc$B_ModelFit[grepl("^B3", mc$B_ModelFit$model), , drop = FALSE]
fit_B0 <- mc$B_ModelFit[grepl("^B0_", mc$B_ModelFit$model) | mc$B_ModelFit$model == "B0_RE", , drop = FALSE]
fit_B4 <- mc$B4_ModelFit

# B3 fragility caveat, reused by the sheet note and the overview markdown.
.b3_caveat <- paste("FRAGILE / secondary: the B3 focal nonlinearity is significant only under the",
  "region-specific ParEdu-smooth adjustment (Plan_deviations.md section 9); without that adjustment",
  "the preregistered focal is null, and the West cohort bend drops below significance under",
  "leave-one-study-out. Treat as suggestive.")

sheet_B1 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = B1_Headline,
       note = paste("B1: linear mixed model of mobility on the education PGI, birth year, region and gender",
                    "(all 2-/3-way interactions) plus the preregistered Eq-4 parental-education interaction set,",
                    "with a per-family random intercept. Focal: PGI x birth year x region.")),
  list(header = "Model coefficients", df = mc$B1_Coefficients,
       note = "Parametric coefficients (per 1 SD of the PGI/birth year); ParEdu Eq-4 interactions included; family random intercept not shown."),
  list(header = "Model fit", df = fit_B1, note = "Fit statistics for the B1 model."),
  list(header = "Leave-one-study-out sensitivity", df = mc$B1_LOO,
       note = "Focal coefficient re-estimated with each study held out (BASE-II / SOEP / TwinLife; SHIP absent from the mobility sample)."),
  title = .t(1, "Linear model results for educational mobility across birth years and region"))

sheet_B2 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = B2_Headline,
       note = paste("B2: as B1 but with a reunification step (born >= 1975) plus the Eq-4 ParEdu interactions.",
                    "Focal: PGI x post-reunification x region.")),
  list(header = "Model coefficients", df = mc$B2_Coefficients, note = "Parametric coefficients; ParEdu Eq-4 interactions included."),
  list(header = "Model fit", df = fit_B2, note = "Fit statistics for the B2 model."),
  list(header = "Leave-one-study-out sensitivity", df = mc$B2_LOO, note = "Focal coefficient re-estimated with each study held out."),
  title = .t(2, "Step-function (reunification) results for educational mobility"))

sheet_B3 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = B3_Headline,
       note = paste("B3: region-varying GAM of the PGI-mobility slope over birth year, controlling",
                    "parental education with region-specific ParEdu smooths (Plan_deviations.md section 9).",
                    "The linear PGI x cohort(x region) trend is tested in B1 (Wald); B3 adds only the",
                    "non-linear departure.", .b3_caveat)),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$B3_NestedLRT,
       note = paste("Baseline = PGI main effect (constant slope); linear adds the cohort slope;",
                    "nonlinear adds the smooth. The LRT column carries only linear -> nonlinear.")),
  # (b) Per-region curvature --------------------------------------------
  list(header = "(b) Per-region curvature test (clean nested LRT)", df = mc$B3_cell_LRT,
       note = paste("Adds only that region's PGI smooth to the linear model - the clean per-region",
                    "non-linearity test (supersedes the confounded summary.gam F; see",
                    "Plan_deviations.md section 6).")),
  # (c) East-West contrast ----------------------------------------------
  list(header = "(c) East-West difference smooth (planned contrast)", df = B3_diff_intervals,
       note = "Birth-year ranges where the East-West difference in the PGI-mobility slope has a 95% simultaneous CI excluding zero (see the difference-smooth figure)."),
  # Robustness + diagnostics --------------------------------------------
  list(header = "Robustness: leave-one-study-out", df = mc$B3_LOO, note = "Linear-vs-nonlinear LRT re-run with each study held out."),
  list(header = "Diagnostics: k-check", df = mc$B3_K_Index, note = "Basis-dimension adequacy for each smooth."),
  title = .t(3, "Region-varying nonlinear (GAM) results for educational mobility"))

.b0_blocks <- list(
  list(header = "B0 focal analysis summary", df = B0_Headline,
       note = paste("B0: overall East-West contrast in the PGI-mobility slope (no birth-year x PGI term),",
                    "controlling parental education, family random intercept. Focal: PGI x region.")),
  list(header = "B0 per-region PGI-mobility slopes", df = mc$B0_RegionSlopes, note = "PGI-mobility slope per region, with 95% CI."),
  list(header = "B0 coefficients", df = mc$B0_Coefficients, note = "Full B0 coefficients."),
  list(header = "B0 leave-one-study-out (focal PGI x region)", df = mc$B0_LOO, note = "Focal contrast re-estimated with each study held out."))
if (!is.null(mc$B0_pre_RegionSlopes)) {
  .b0_blocks <- c(.b0_blocks, list(
    list(header = "B0 (pre-reunification) per-region slopes", df = mc$B0_pre_RegionSlopes, note = "As above, born < 1975."),
    list(header = "B0 (pre-reunification) coefficients", df = mc$B0_pre_Coefficients, note = "Pre-reunification B0 coefficients.")))
  if (!is.null(mc$B0_pre_LOO))
    .b0_blocks <- c(.b0_blocks, list(list(header = "B0 (pre-reunification) leave-one-study-out", df = mc$B0_pre_LOO, note = "Pre-reunification focal contrast, study held out.")))
}
.b0_blocks <- c(.b0_blocks, list(list(header = "Model fit (B0 / B0_pre)", df = fit_B0, note = "Fit statistics.")))
sheet_B0 <- do.call(multiblock_sheet, c(.b0_blocks, list(title = .t(4, "Overall East-West PGI-mobility contrast"))))

sheet_B4 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = B4_Headline,
       note = paste("B4: linear four-way PGI x birth-year x region x gender mixed model with four-way ParEdu",
                    "control (Plan_deviations.md §10). Focal: the four-way interaction.")),
  list(header = "Model coefficients", df = mc$B4_Coefficients, note = "Full factorial parametric coefficients (PGI and ParEdu four-way sets)."),
  list(header = "Model fit", df = fit_B4, note = "Fit statistics for the B4 model."),
  list(header = "Leave-one-study-out sensitivity", df = mc$B4_LOO, note = "Four-way focal coefficient re-estimated with each study held out."),
  list(header = "Leave-one-study-out sensitivity (pooled PGI x gender)", df = nz(mc$B4_LOO_gender),
       note = paste("Pooled PGI x gender interaction (not the four-way term above) re-estimated with each",
                    "study held out; estimate, SE, 95% CI, exact p, supported birth-year range and change",
                    "vs full. Additional sensitivity - see Plan_deviations.md section 14.")),
  title = .t(5, "Linear gender-interaction results for educational mobility"))

sheet_B5 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = B5_Headline,
       note = paste("B5: per-gender GAM (regions pooled, family random intercept) of the PGI-mobility",
                    "slope over birth year, with gender-specific ParEdu control. The linear gender x",
                    "cohort trend is tested in B4 (Wald); B5 adds only the non-linear departure.")),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$B5_NestedLRT,
       note = paste("Baseline = PGI main effect (constant slope); linear adds the per-gender cohort",
                    "slopes; nonlinear adds the smooths. The LRT column carries only linear -> nonlinear.")),
  # (b) Per-gender curvature --------------------------------------------
  list(header = "(b) Per-gender curvature test (clean nested LRT)", df = mc$B5_cell_LRT,
       note = paste("Adds only that gender's PGI smooth to the linear model - the clean per-gender",
                    "non-linearity test (supersedes the confounded summary.gam F).")),
  # (c) Female-Male contrast --------------------------------------------
  list(header = "(c) Female-Male difference smooth (planned contrast)", df = B5_diff_intervals,
       note = "Birth-year ranges where the female-male difference in the PGI-mobility slope has a 95% simultaneous CI excluding zero."),
  # (d) Leave-one-study-out --------------------------------------------
  list(header = "(d) Leave-one-study-out (omnibus + per-gender curvature)", df = nz(mc$B5_LOO),
       note = paste("Per fold: the omnibus linear-vs-nonlinear LRT plus the female/male curvature tests,",
                    "with the supported birth-year range. Checks the (male-driven) nonlinear cohort pattern",
                    "is not carried by any one study. Additional sensitivity - see Plan_deviations.md section 14.")),
  list(header = "(e) Female-Male difference smooth, leave-one-study-out", df = nz(mc$B5_LOO_diff_intervals),
       note = paste("The birth-year intervals where the female-male difference smooth's 95% simultaneous CI",
                    "excludes zero, recomputed with each study held out (restricted to each fold's",
                    "supported birth-year range).")),
  title = .t(6, "Per-gender nonlinear (GAM) results for educational mobility"))

sheet_B6 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = B6_Headline,
       note = paste("B6: region x gender GAM of the PGI-mobility slope over birth year, with cell-specific",
                    "ParEdu control (Plan_deviations.md section 10). The linear four-way trend is tested",
                    "in B4 (Wald); B6 adds only the non-linear departure.")),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$B6_NestedLRT,
       note = paste("Baseline = PGI main effects (constant slopes); linear adds the four cell cohort",
                    "slopes; nonlinear adds the smooths. The LRT column carries only linear -> nonlinear.")),
  # (b) Per-cell curvature ----------------------------------------------
  list(header = "(b) Per-cell curvature test (clean nested LRT)", df = mc$B6_cell_LRT,
       note = paste("Adds only that region x gender cell's PGI smooth to the linear model - the clean",
                    "per-cell non-linearity test (supersedes the confounded summary.gam F).")),
  # (c) Female-Male contrast within region ------------------------------
  list(header = "(c) Female-Male difference smooth within region (planned contrast)", df = B6_diff_intervals,
       note = "Birth-year ranges where the female-male difference in the PGI-mobility slope has a 95% simultaneous CI excluding zero, within West and within East."),
  # Descriptive + robustness + diagnostics ------------------------------
  list(header = "Overall additive region x gender PGI slopes", df = mc$B6_OverallSlopes, note = "PGI-mobility slope per region x gender cell from the additive model."),
  list(header = "Robustness: leave-one-study-out", df = mc$B6_LOO, note = "Linear-vs-nonlinear LRT re-run with each study held out."),
  list(header = "Diagnostics: k-check", df = mc$B6_K_Index, note = "Basis-dimension adequacy for each smooth."),
  title = .t(7, "Region x gender nonlinear (GAM) results for educational mobility"))

.b0g_blocks <- list(
  list(header = "Focal analysis summary", df = B0g_Headline,
       note = paste("B0g: overall region x gender PGI-mobility contrast (no birth-year x PGI term), with",
                    "region x gender ParEdu control. Focal: PGI x region x gender.")),
  list(header = "Model coefficients (full sample)", df = mc$B0g$coefs, note = "Parametric coefficients on the full mobility sample."),
  list(header = "Per-cell region x gender PGI-mobility slopes (full sample)", df = mc$B0g$cell_slopes, note = "PGI-mobility slope per region x gender cell, with 95% CI and p."),
  list(header = "Model fit (full sample)", df = mc$B0g$fit, note = "Fit statistics."))
if (!is.null(mc$B0g_pre)) {
  .b0g_blocks <- c(.b0g_blocks, list(
    list(header = "Model coefficients (pre-reunification)", df = mc$B0g_pre$coefs, note = "Restricted to born < 1975."),
    list(header = "Per-cell PGI-mobility slopes (pre-reunification)", df = mc$B0g_pre$cell_slopes, note = "Per-cell slopes, pre-reunification."),
    list(header = "Model fit (pre-reunification)", df = mc$B0g_pre$fit, note = "Fit statistics.")))
}
sheet_B0g <- do.call(multiblock_sheet, c(.b0g_blocks, list(title = .t(8, "Overall region x gender contrast for educational mobility"))))

# ---- 09: regional social-origin gradient ----
# Its own sheet rather than a block on 01_Linear_birth_year: a different model
# from B1, reported on the educational-attainment metric while sheet 01 is on
# the mobility metric.
.grad_note <- paste(
  "Gradient of educational attainment on parental education, per 1 SD of kernel-standardised",
  "parental education. Fitted directly on educational attainment as parental education x",
  "region, adjusted for birth year and gender. east_west_c is coded East +0.5 / West -0.5,",
  "so 'ParEdu x region' is the East-minus-West difference and the pooled row is the grand mean.")
# The per-row `spec` column repeats .grad_note on every row; drop it so the table
# reads like the rest of the workbook.
.grad_tbl <- nz(ParEdu_Gradient)
.grad_tbl <- .grad_tbl[, setdiff(names(.grad_tbl), "spec"), drop = FALSE]
sheet_B1grad <- multiblock_sheet(
  list(header = "Regional social-origin gradient in educational attainment",
       df = .grad_tbl, note = .grad_note),
  title = .t(9, "Regional parental-education gradient in educational attainment"))

sheets <- list(
  "00_Index"             = sheet_README,
  "01_Linear_birth_year" = sheet_B1,
  "02_Step_reunification"= sheet_B2,
  "03_GAM_region"        = sheet_B3,
  "04_Overall_EastWest"  = sheet_B0,
  "05_Linear_gender"     = sheet_B4,
  "06_GAM_gender"        = sheet_B5,
  "07_GAM_region_gender" = sheet_B6,
  "08_Overall_RxG"       = sheet_B0g,
  "09_ParEdu_gradient"   = sheet_B1grad)
sheets <- humanize_sheets(sheets)
OUT_XLSX <- file.path(OUT_DIR_B, "Supplementary_Data_5_Mobility_results.xlsx")
write_section_xlsx("B", sheets, OUT_XLSX)
cat(sprintf("B_report: xlsx written to %s (%d sheets)\n", OUT_XLSX, length(sheets)))

# ---- Plots ----
# The PDFs are the expensive stage: the workbook above is already written, and
# everything below is prediction plus ggplot rendering. A wording-only change to
# the workbook does not need them, so REPORT_SKIP_PLOTS=1 stops here. The xlsx and
# the Overview markdown are unaffected either way.
if (identical(Sys.getenv("REPORT_SKIP_PLOTS"), "1")) {
  cat("B_report: REPORT_SKIP_PLOTS=1 -- skipping the PDF stage; xlsx and overview still written.\n")
} else {

cat("B_report: building plots...\n")
kern_mob <- dplyr::bind_rows(
  kernel_moments(dcm$mobility[dcm$east_west == "West"], dcm$birth_year[dcm$east_west == "West"], bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "West"),
  kernel_moments(dcm$mobility[dcm$east_west == "East"], dcm$birth_year[dcm$east_west == "East"], bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "East"))
p_kern <- ggplot(kern_mob, aes(x = t, colour = Region, fill = Region)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
  geom_line(aes(y = mu), linewidth = 0.9) + geom_point(aes(y = mu, size = n_raw), alpha = 0.3) +
  scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_fill_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_size_continuous(range = c(0.3, 3), guide = "none") +
  labs(x = "Birth year", y = "Mobility (cohort-z difference)",
       title = "Kernel-smoothed intergenerational mobility by birth year",
       subtitle = sprintf("Bandwidth = %d years, shaded = +/- 1 SD; one smooth per region.", KERNEL_BW)) + .theme

np_B <- compute_binned_slopes(dcm, outcome = "mobility", pgi = "PGI_Edu_z", covariates = "parental_edu_z_kernel", group_vars = "east_west")
np_B_gender <- tryCatch(compute_binned_slopes(dcm, outcome = "mobility", pgi = "PGI_Edu_z",
  covariates = "parental_edu_z_kernel", group_vars = c("east_west", "gender")), error = function(e) NULL)

curves_B1 <- dplyr::bind_rows(
  extract_pgi_slope(mc$m_B1_re, dcm, east_west = "West", pgi_var = "PGI_Edu_z"),
  extract_pgi_slope(mc$m_B1_re, dcm, east_west = "East", pgi_var = "PGI_Edu_z"))
p_B1 <- plot_combined_east_west(curves_B1, np_B, title = "B1: Linear PGI-mobility slope by Region",
  ylab = "Standardized beta (SD mobility per 1 SD PGI)", hist_dat = dcm, hist_alpha = 0.10, hist_bins = 50)

b2_eff <- extract_a2_pgi_effects(mc$m_B2_re, dcm)
p_B2 <- ggplot(b2_eff, aes(x = period, y = slope, colour = east_west, group = east_west)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.08, linewidth = 0.7) +
  geom_line(linewidth = 1) + geom_point(size = 3.5) + .scale_col +
  scale_x_discrete(labels = c(pre = "Pre-reunification", post = "Post-reunification")) +
  labs(x = "Period", y = "Standardized beta (SD mobility per 1 SD PGI)",
       title = "B2: PGI-mobility by reunification period", colour = "Region") + .theme

# B3 display curves: nonlinear smooth only where the per-region LRT is significant.
p_B3 <- plot_combined_east_west(mc$B3_display_curves, np_B,
  title = "B3: Cohort-varying PGI-mobility slope by Region (smooth shown only where curvature is significant)",
  ylab = "Standardized beta (SD mobility per 1 SD PGI)", hist_dat = dcm, hist_alpha = 0.10, hist_bins = 50)
p_B3_diff <- plot_diff_smooth(diff_B, title = "B3 Mobility: East-West standardized PGI difference",
  hist_dat = dcm, hist_alpha = 0.10, hist_bins = 50)

# B0 overall/pre region slopes.
b0_eff <- dplyr::bind_rows(
  dplyr::mutate(mc$B0_RegionSlopes, analysis = "Overall"),
  if (!is.null(mc$B0_pre_RegionSlopes)) dplyr::mutate(mc$B0_pre_RegionSlopes, analysis = "Pre-reunification") else NULL)
b0_eff$analysis <- factor(b0_eff$analysis, levels = c("Overall", "Pre-reunification")[c(TRUE, !is.null(mc$B0_pre_RegionSlopes))])
.dodge <- position_dodge(width = 0.4)
p_B0 <- ggplot(b0_eff, aes(x = analysis, y = slope, colour = east_west, group = east_west)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, linewidth = 0.7, position = .dodge) +
  geom_point(size = 3.5, position = .dodge) + .scale_col +
  labs(x = NULL, y = "Standardized beta (SD mobility per 1 SD PGI)",
       title = "B0: PGI-mobility slope by region", colour = "Region") + .theme

# B5 per-gender curve; B6 4-cell + diff smooths.
p_B5 <- tryCatch(
  ggplot(mc$B5_curves, aes(x = birth_year, y = slope, colour = gender, fill = gender, linetype = gender)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 1) +
    scale_colour_manual(values = c(female = "#7B3294", male = "#008837")) +
    scale_fill_manual(values = c(female = "#7B3294", male = "#008837")) +
    scale_linetype_manual(values = c(female = LT_FEMALE, male = LT_MALE)) +
    labs(x = "Birth year", y = "Standardized beta (SD mobility per 1 SD PGI)",
         title = "B5: Per-gender PGI-mobility slope over birth year (regions pooled; smooth shown only where curvature is significant)",
         colour = "Gender", fill = "Gender", linetype = "Gender") + .theme,
  error = function(e) { cat("p_B5 failed:", conditionMessage(e), "\n"); NULL })
p_B6 <- tryCatch(plot_combined_east_west_gender(mc$B6_display_curves, np_B_gender,
  title = "B6: Cohort-varying PGI-mobility slope by Region x Gender (smooth shown only where curvature is significant)",
  ylab = "Standardized beta (SD mobility per 1 SD PGI)", hist_dat = dcm, hist_alpha = 0.10, hist_bins = 50),
  error = function(e) { cat("p_B6 failed:", conditionMessage(e), "\n"); NULL })
p_B6_diff_W <- tryCatch(plot_diff_smooth(diff_smooth_gender_simul(m_B6_nonlin, dcm, east_west = "West"),
  title = "B6 West: Female-Male standardized PGI-mobility difference", hist_dat = dcm[dcm$east_west == "West", ], hist_alpha = 0.10),
  error = function(e) NULL)
p_B6_diff_E <- tryCatch(plot_diff_smooth(diff_smooth_gender_simul(m_B6_nonlin, dcm, east_west = "East"),
  title = "B6 East: Female-Male standardized PGI-mobility difference", hist_dat = dcm[dcm$east_west == "East", ], hist_alpha = 0.10),
  error = function(e) NULL)

OUT_PDF <- file.path(OUT_DIR_B, "B_Mobility_Plots.pdf")
pdf(OUT_PDF, width = 10, height = 7)
for (p in Filter(Negate(is.null), list(p_kern, p_B1, p_B2, p_B3, p_B3_diff, p_B0, p_B5, p_B6, p_B6_diff_W, p_B6_diff_E)))
  tryCatch(print(p), error = function(e) cat("Plot failed:", conditionMessage(e), "\n"))
dev.off()
cat(sprintf("B_report: pdf written to %s\n", OUT_PDF))

# ---- Leave-one-study-out PDF: B3 per-fold curves + kernel-mobility per fold + B0 forest/slopes ----
.folds_map <- cohort_groups(dcm)
if (!is.null(mc$B3_LOO_curves) && nrow(mc$B3_LOO_curves) > 0) {
  loo_curves <- mc$B3_LOO_curves
  fold_levels <- if (is.factor(loo_curves$fold)) levels(loo_curves$fold) else unique(loo_curves$fold)
  .trim_one <- function(fl) {
    fc <- loo_curves[as.character(loo_curves$fold) == fl, , drop = FALSE]; if (!nrow(fc)) return(fc)
    sub <- if (fl == "Full sample") dcm else dcm[!(as.character(dcm$cohort) %in% .folds_map[[sub("^-", "", fl)]]), , drop = FALSE]
    .trim_to_data_range(fc, sub, group_vars = "east_west")
  }
  loo_curves <- dplyr::bind_rows(lapply(fold_levels, .trim_one)); loo_curves$fold <- factor(loo_curves$fold, levels = fold_levels)
  p_loo_B3 <- ggplot(loo_curves, aes(x = birth_year, y = slope, colour = east_west, fill = east_west)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + facet_wrap(~ fold) + .scale_col + .scale_fill +
    labs(x = "Birth year", y = "Standardized beta (SD mobility per 1 SD PGI)",
         title = "B3 leave-one-study-out: cohort-varying PGI-mobility slope by region",
         subtitle = "Each panel refits the B3 GAM with one study held out.", colour = "Region", fill = "Region") + .theme
  .kern_one <- function(fl) {
    sub <- if (fl == "Full sample") dcm else dcm[!(as.character(dcm$cohort) %in% .folds_map[[sub("^-", "", fl)]]), , drop = FALSE]
    dplyr::bind_rows(
      kernel_moments(sub$mobility[sub$east_west == "West"], sub$birth_year[sub$east_west == "West"], bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "West", fold = fl),
      kernel_moments(sub$mobility[sub$east_west == "East"], sub$birth_year[sub$east_west == "East"], bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "East", fold = fl))
  }
  kern_loo <- dplyr::bind_rows(lapply(fold_levels, .kern_one)); kern_loo$fold <- factor(kern_loo$fold, levels = fold_levels)
  p_loo_kern <- ggplot(kern_loo, aes(x = t, colour = Region, fill = Region)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mu), linewidth = 0.8) + geom_point(aes(y = mu, size = n_raw), alpha = 0.3) + facet_wrap(~ fold) +
    scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_fill_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_size_continuous(range = c(0.2, 2.5), guide = "none") +
    labs(x = "Birth year", y = "Mobility (cohort-z difference)",
         title = "Leave-one-study-out: kernel-smoothed mobility by birth year",
         subtitle = sprintf("Bandwidth = %d years; one study held out per panel.", KERNEL_BW), colour = "Region", fill = "Region") + .theme
  loo_pages <- list(p_loo_B3, p_loo_kern)
  .b0_forest_df <- function(tbl, lbl) { if (is.null(tbl) || !nrow(tbl)) return(NULL)
    data.frame(fold = tbl$fold, estimate = tbl$estimate, ci_lo = tbl$estimate - 1.96 * tbl$SE, ci_hi = tbl$estimate + 1.96 * tbl$SE, analysis = lbl, stringsAsFactors = FALSE) }
  b0_forest <- dplyr::bind_rows(.b0_forest_df(mc$B0_LOO, "Overall"), .b0_forest_df(mc$B0_pre_LOO, "Pre-reunification"))
  if (!is.null(b0_forest) && nrow(b0_forest) > 0) {
    flev <- unique(b0_forest$fold); b0_forest$fold <- factor(b0_forest$fold, levels = rev(flev))
    b0_forest$analysis <- factor(b0_forest$analysis, levels = intersect(c("Overall", "Pre-reunification"), unique(b0_forest$analysis)))
    p_loo_B0 <- ggplot(b0_forest, aes(x = estimate, y = fold)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.2, linewidth = 0.7) +
      geom_point(size = 3) + facet_wrap(~ analysis) +
      labs(x = "B0 focal: East-West PGI-mobility contrast (standardized beta)", y = NULL,
           title = "B0 leave-one-study-out: East-West PGI-mobility contrast",
           subtitle = "Each row refits B0 with one study held out; bars = 95% CI.") + .theme
    loo_pages <- c(loo_pages, list(p_loo_B0))
  }
  if (!is.null(mc$B0_LOO_slopes) && nrow(mc$B0_LOO_slopes) > 0) {
    sl <- mc$B0_LOO_slopes; sl$fold <- factor(sl$fold, levels = unique(sl$fold))
    sl$analysis <- factor(sl$analysis, levels = intersect(c("Overall", "Pre-reunification"), unique(sl$analysis)))
    p_loo_B0_slopes <- ggplot(sl, aes(x = analysis, y = slope, colour = east_west, group = east_west)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, linewidth = 0.6, position = position_dodge(width = 0.4)) +
      geom_point(size = 2.8, position = position_dodge(width = 0.4)) + facet_wrap(~ fold) + .scale_col +
      labs(x = NULL, y = "Standardized beta (SD mobility per 1 SD PGI)",
           title = "B0 leave-one-study-out: PGI-mobility slope by region",
           subtitle = "Each panel refits B0 (and B0_pre) with one study held out.", colour = "Region") + .theme
    loo_pages <- c(loo_pages, list(p_loo_B0_slopes))
  }
  OUT_LOO_PDF <- file.path(OUT_DIR_B, "B_Mobility_LOO_Plots.pdf")
  pdf(OUT_LOO_PDF, width = 11, height = 7)
  for (p in loo_pages) tryCatch(print(p), error = function(e) cat("LOO plot failed:", conditionMessage(e), "\n"))
  dev.off()
  cat(sprintf("B_report: LOO pdf written to %s (%d pages)\n", OUT_LOO_PDF, length(loo_pages)))
}
}

# ---- Section overview markdown ----
b3_concl <- unique(B3_GAM_Overview$nonlinearity_conclusion)
b3_concl <- if (length(b3_concl) == 1) b3_concl else paste0("Mixed: ", paste(b3_concl, collapse = "; "))
get_focal <- function(tbl, term) { r <- if (is.null(tbl)) NULL else tbl[tbl$term == term, ]; if (is.null(r) || !nrow(r)) c(NA, NA) else c(r$estimate, r$p.value) }
b1f <- get_focal(mc$B1_Coefficients, "PGI_Edu_z:BYc_z:east_west_c")
b2f <- get_focal(mc$B2_Coefficients, "PGI_Edu_z:reunifpost:east_west_c")
b0f <- get_focal(mc$B0_Coefficients, "PGI_Edu_z:east_west_c")
b4f <- get_focal(mc$B4_Coefficients, mc$B4_focal)
b0gf <- get_focal(mc$B0g$coefs, mc$B0g_focal)
analyses_tbl <- data.frame(
  Analysis = c("B1","B2","B3","B0/B0_pre","B4","B5","B6","B0g"),
  Question = c("Does the PGI-mobility association change linearly across cohorts, by region?",
               "Does it change before vs. after the reunification cutoff, by region?",
               "Is a nonlinear birth-cohort pattern needed (controlling parental education)?",
               "Overall / pre-reunification East-West difference in the PGI-mobility slope",
               "Does the cohort x region pattern differ by gender (four-way)?",
               "Does the PGI-mobility slope's birth-cohort trajectory differ for women vs. men?",
               "Does the cohort-varying PGI-mobility slope differ across Region x Gender cells?",
               "Does the overall PGI-mobility slope's gender gap differ by region?"),
  Model = c("Linear RE", "Step-function RE", "GAM RE (ParEdu-smoothed)", "Linear RE",
            "Linear four-way RE", "Per-gender GAM RE", "Region x Gender GAM RE", "Overall R x G RE"),
  stringsAsFactors = FALSE)
findings_tbl <- data.frame(
  Analysis = c("B1","B2","B3","B0","B4","B6","B0g"),
  `Focal test` = c("PGI x birth year x region", "PGI x post-1975 x region",
                   "Nonlinear GAM vs. linear (ParEdu-controlled)", "PGI x region (overall)",
                   "PGI x birth year x region x gender", "Region x Gender nonlinear GAM vs. linear", "PGI x region x gender"),
  `Estimate / statistic` = c(
    if (is.na(b1f[1])) "n/a" else sprintf("beta = %.4f", b1f[1]),
    if (is.na(b2f[1])) "n/a" else sprintf("beta = %.4f", b2f[1]),
    { r <- mc$B3_NestedLRT[which(mc$B3_NestedLRT$comparison == "linear vs nonlinear"), ]; if (nrow(r)) sprintf("chi2 = %.2f", r$LRT_chisq) else "n/a" },
    if (is.na(b0f[1])) "n/a" else sprintf("beta = %.4f", b0f[1]),
    if (is.na(b4f[1])) "n/a" else sprintf("beta = %.4f", b4f[1]),
    { r <- mc$B6_NestedLRT[which(mc$B6_NestedLRT$comparison == "linear vs nonlinear"), ]; if (nrow(r)) sprintf("chi2 = %.2f", r$LRT_chisq) else "n/a" },
    if (is.na(b0gf[1])) "n/a" else sprintf("beta = %.4f", b0gf[1])),
  p = c(.fmt_p(b1f[2]), .fmt_p(b2f[2]),
        { r <- mc$B3_NestedLRT[which(mc$B3_NestedLRT$comparison == "linear vs nonlinear"), ]; if (nrow(r)) .fmt_p(r$LRT_p) else "NA" },
        .fmt_p(b0f[2]), .fmt_p(b4f[2]),
        { r <- mc$B6_NestedLRT[which(mc$B6_NestedLRT$comparison == "linear vs nonlinear"), ]; if (nrow(r)) .fmt_p(r$LRT_p) else "NA" },
        .fmt_p(b0gf[2])),
  Conclusion = c(
    if (!is.na(b1f[2]) && b1f[2] < .05) "Focal interaction significant" else "No focal interaction at alpha=.05",
    if (!is.na(b2f[2]) && b2f[2] < .05) "Focal interaction significant" else "No focal interaction at alpha=.05",
    paste0(b3_concl, " -- ", .b3_caveat),
    if (!is.na(b0f[2]) && b0f[2] < .05) "Significant East-West difference" else "No East-West difference at alpha=.05",
    if (!is.na(b4f[2]) && b4f[2] < .05) "Four-way interaction significant" else "No four-way interaction at alpha=.05",
    { cc <- unique(B6_GAM_Overview$nonlinearity_conclusion); if (length(cc) == 1) cc else paste(cc, collapse = "; ") },
    if (!is.na(b0gf[2]) && b0gf[2] < .05) "Region x Gender interaction significant" else "No Region x Gender interaction at alpha=.05"),
  check.names = FALSE, stringsAsFactors = FALSE)
outputs_tbl <- data.frame(
  Output = c("Excel workbook","Figure PDF","LOO figure PDF","Model cache"),
  File = c("Supplementary_Data_5_Mobility_results.xlsx","B_Mobility_Plots.pdf","B_Mobility_LOO_Plots.pdf","B_Mobility_Models.rds"),
  stringsAsFactors = FALSE)
interp <- paste0(
  "Section B tests whether the PGI-mobility association varies by birth cohort, region and gender, ",
  "controlling for parental education throughout (B1/B2 with the preregistered Eq-4 ParEdu interactions; ",
  "B3/B5/B6 with ParEdu smooths at the focal cell flexibility — Plan_deviations.md §9/§10). Gender is a ",
  "main moderator: B4 (four-way), B5 (per-gender GAM), B6 (Region x Gender GAM), B0g (overall additive ",
  "contrast); B4/B6/B0g are the Section C mobility analyses of the preregistration, estimated with family ",
  "random effects. B0/B0_pre are the overall (and pre-reunification) East-West contrasts. ", .b3_caveat,
  " Each focal result is checked with leave-one-study-out (BASE-II / SOEP / TwinLife; SHIP is absent from ",
  "the mobility sample for lack of parental education). Where higher-order interactions are null, the ",
  "lower-order / additive moderation is the reported effect.")

OUT_MD <- file.path(OUT_DIR_B, "B_Mobility_Overview.md")
generate_section_overview_md("B", "Section B — Educational mobility (cohort, region, gender)",
  paste0("Section B tests whether the association between the education PGI and relative intergenerational ",
         "mobility differs across birth cohorts, between East and West Germany, and by gender, controlling ",
         "for parental education. RE sample dat_cluster_mob (N = ", N_b, "; BASE-II, SOEP, TwinLife — SHIP ",
         "absent, no parental education). Gender (B4/B5/B6/B0g) is a main moderator."),
  analyses_tbl, findings_tbl, interp, outputs_tbl, OUT_MD)
cat(sprintf("B_report: overview written to %s\n", OUT_MD))
cat("\nB_report: DONE.\n")
