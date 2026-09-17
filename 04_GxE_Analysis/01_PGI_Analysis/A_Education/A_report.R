# A_report.R
# Build Supplementary_Data_2_Education_results.xlsx, A_Education_Plots.pdf,
# A_Education_LOO_Plots.pdf and A_Education_Overview.md from the model cache
# produced by A_analysis.R.
#
# Reads A_Education_Models.rds from runs/RUN_<TS>/A_Education/.
# Set REPORT_SKIP_PLOTS=1 to write the workbook and overview without the PDFs.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(tibble); library(mgcv); library(openxlsx)
})

# ---- Resolve paths ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  # Rscript replaces spaces in the --file= argument with `~+~`; undo that
  # before normalizePath so paths like ".../015_Data analysis/..." resolve.
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine A_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR  # flattened: 01_PGI_Analysis is the workflow + script root
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"),
       local = TRUE)

# ---- Load model cache ----
OUT_DIR_A <- output_dir("A_Education")
RDS_PATH  <- file.path(OUT_DIR_A, "A_Education_Models.rds")
if (!file.exists(RDS_PATH)) {
  stop("A_Education_Models.rds not found at ", RDS_PATH,
       " — run A_analysis.R first.")
}
mc <- readRDS(RDS_PATH)
# Optional tables degrade to a "not available" note rather than crashing the
# workbook build.
nz <- function(x, fallback = data.frame(note = "not run / not available", stringsAsFactors = FALSE))
  if (is.null(x)) fallback else x
dat_cluster  <- mc$dat_cluster
m_A1_re      <- mc$m_A1_re
m_A2_re      <- mc$m_A2_re
m_A3_nonlin  <- mc$m_A3_nonlin_re
m_A3_cr      <- mc$m_A3_cr_re
# Gender block (Section C analyses + per-gender attainment GAM, RE)
m_A4         <- mc$m_A4_re
m_A5_nonlin  <- mc$m_A5_nonlin_re
m_A6_nonlin  <- mc$m_A6_nonlin_re

cat(sprintf("A_report: loaded models from %s (dat_cluster N = %d)\n",
            RDS_PATH, nrow(dat_cluster)))

# ---- Run metadata sheet ----
run_metadata <- build_run_metadata("A_Education", mc$RUN_TS, SCRIPT_DIR)

# ---- Headline tables ----
cohorts_a <- sort(as.character(unique(dat_cluster$cohort)))
N_a       <- nrow(dat_cluster)

# A1 headline
A1_focal_term <- "PGI_Edu_z:BYc_z:east_west_c"
a1_row <- mc$A1_Coefficients[mc$A1_Coefficients$term == A1_focal_term, ]
A1_Headline <- if (nrow(a1_row) > 0) {
  build_headline_table(
    analysis_id = "A1", status = "preregistered", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Linear mixed model (family RE)",
    focal_term = A1_focal_term,
    estimate = round(a1_row$estimate, 4),
    SE_or_CI = sprintf("SE=%.4f (95%% CI [%.4f, %.4f])",
                       a1_row$std.error, a1_row$conf.low, a1_row$conf.high),
    p = round(a1_row$p.value, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#01_Linear_birth_year",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else {
  data.frame(note = "A1 focal term not found in coefficient table.",
             stringsAsFactors = FALSE)
}

# A2 headline (post coefficient on the interaction with reunifpost)
A2_focal_term <- "PGI_Edu_z:reunifpost:east_west_c"
a2_row <- mc$A2_Coefficients[mc$A2_Coefficients$term == A2_focal_term, ]
A2_Headline <- if (nrow(a2_row) > 0) {
  build_headline_table(
    analysis_id = "A2", status = "preregistered", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Step-function mixed model (family RE)",
    focal_term = A2_focal_term,
    estimate = round(a2_row$estimate, 4),
    SE_or_CI = sprintf("SE=%.4f (95%% CI [%.4f, %.4f])",
                       a2_row$std.error, a2_row$conf.low, a2_row$conf.high),
    p = round(a2_row$p.value, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#02_Step_reunification",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else {
  data.frame(note = "A2 focal term not found in coefficient table.",
             stringsAsFactors = FALSE)
}

# A3 headline (omnibus nested LRT, linear vs nonlinear).
# Use which() to drop the NA-comparison row (the null model has no
# comparison string and `NA == "linear vs nonlinear"` returns NA, which
# would otherwise produce a phantom all-NA row in the subset.)
a3_lrt_row <- mc$A3_NestedLRT[which(mc$A3_NestedLRT$comparison == "linear vs nonlinear"), ]
A3_Headline <- if (nrow(a3_lrt_row) > 0) {
  build_headline_table(
    analysis_id = "A3", status = "preregistered", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Varying-coefficient GAM (family RE)",
    focal_term = "omnibus nested LRT (linear vs nonlinear)",
    estimate = sprintf("chi2 = %.2f, df = %.2f",
                       a3_lrt_row$LRT_chisq, a3_lrt_row$LRT_df),
    SE_or_CI = "n/a (see A3_GAM_Overview for per-smooth tests)",
    p = round(a3_lrt_row$LRT_p, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#03_GAM_region",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else {
  data.frame(note = "A3 nested LRT row not found.", stringsAsFactors = FALSE)
}

# ---- A3 GAM overview (nonlinearity_conclusion schema) ----
A3_GAM_Overview <- build_gam_overview_table(
  nested_lrt   = mc$A3_NestedLRT,
  smooth_table = mc$A3_Smooths_tp,
  sample_name  = "dat_cluster",
  N            = N_a,
  analysis_id  = "A3",
  alpha        = 0.05,
  cell_lrt     = mc$A3_cell_LRT
)

# ---- Gender headlines + GAM overviews (A4, A5, A6, A0g) ----
a4_row <- mc$A4_Coefficients[mc$A4_Coefficients$term == mc$A4_focal, ]
A4_Headline <- if (nrow(a4_row) > 0) {
  build_headline_table(
    analysis_id = "A4", status = "preregistered", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Linear four-way mixed model (family RE)",
    focal_term = mc$A4_focal, estimate = round(a4_row$estimate, 4),
    SE_or_CI = sprintf("SE=%.4f (95%% CI [%.4f, %.4f])",
                       a4_row$std.error, a4_row$conf.low, a4_row$conf.high),
    p = round(a4_row$p.value, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#04_Linear_gender",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else data.frame(note = "A4 focal term not found.", stringsAsFactors = FALSE)

a5_lrt_row <- mc$A5_NestedLRT[which(mc$A5_NestedLRT$comparison == "linear vs nonlinear"), ]
A5_Headline <- if (nrow(a5_lrt_row) > 0) {
  build_headline_table(
    analysis_id = "A5", status = "main", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Per-gender varying-coefficient GAM (regions pooled, family RE)",
    focal_term = "omnibus nested LRT (linear vs nonlinear)",
    estimate = sprintf("chi2 = %.2f, df = %.2f", a5_lrt_row$LRT_chisq, a5_lrt_row$LRT_df),
    SE_or_CI = "n/a (see A5 per-gender smooth tests)",
    p = round(a5_lrt_row$LRT_p, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#05_GAM_gender",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else data.frame(note = "A5 nested LRT row not found.", stringsAsFactors = FALSE)
A5_GAM_Overview <- build_gam_overview_table(
  nested_lrt = mc$A5_NestedLRT, smooth_table = mc$A5_Smooths_tp,
  sample_name = "dat_cluster", N = N_a, analysis_id = "A5", alpha = 0.05,
  cell_lrt = mc$A5_cell_LRT)

a6_lrt_row <- mc$A6_NestedLRT[which(mc$A6_NestedLRT$comparison == "linear vs nonlinear"), ]
A6_Headline <- if (nrow(a6_lrt_row) > 0) {
  build_headline_table(
    analysis_id = "A6", status = "preregistered", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = N_a,
    cohorts_included = cohorts_a, model_family = "Region x Gender varying-coefficient GAM (family RE)",
    focal_term = "omnibus nested LRT (linear vs nonlinear)",
    estimate = sprintf("chi2 = %.2f, df = %.2f", a6_lrt_row$LRT_chisq, a6_lrt_row$LRT_df),
    SE_or_CI = "n/a (see A6 per-cell smooth tests)",
    p = round(a6_lrt_row$LRT_p, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#06_GAM_region_gender",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else data.frame(note = "A6 nested LRT row not found.", stringsAsFactors = FALSE)
A6_GAM_Overview <- build_gam_overview_table(
  nested_lrt = mc$A6_NestedLRT, smooth_table = mc$A6_Smooths_pgi,
  sample_name = "dat_cluster", N = N_a, analysis_id = "A6", alpha = 0.05,
  cell_lrt = mc$A6_cell_LRT)

a0g_row <- mc$A0g$coefs[mc$A0g$coefs$term == mc$A0g_focal, ]
A0g_Headline <- if (nrow(a0g_row) > 0) {
  build_headline_table(
    analysis_id = "A0g", status = "main", outcome = "edu_z_kernel",
    sample_name = "dat_cluster (full-family TwinLife, RE)", N = mc$A0g$N,
    cohorts_included = cohorts_a, model_family = "Overall Region x Gender PGI contrast (family RE)",
    focal_term = mc$A0g_focal, estimate = round(a0g_row$estimate, 4),
    SE_or_CI = sprintf("SE=%.4f (95%% CI [%.4f, %.4f])",
                       a0g_row$std.error, a0g_row$conf.low, a0g_row$conf.high),
    p = round(a0g_row$p.value, 4),
    output_source = "Supplementary_Data_2_Education_results.xlsx#07_Overall_RxG",
    code_source = "01_PGI_Analysis/A_Education/A_analysis.R")
} else data.frame(note = "A0g focal term not found.", stringsAsFactors = FALSE)

# ---- Assemble xlsx sheets ----
# Layout: an index sheet, then one self-contained sheet per analysis. Each
# analysis sheet stacks its publication headline, its detail tables and its
# model-fit row(s) as labelled blocks (multiblock_sheet).

# Coerce estimate / SE_or_CI / p / N to character so A3 (which reports
# the omnibus chi2 as a string) can bind_rows with A1/A2 (numeric
# estimates). Publication headlines are for human reading, not
# computation; uniform character columns avoid type conflicts.
.headline_to_chr <- function(df) {
  for (col in c("estimate", "SE_or_CI", "p", "N")) {
    if (col %in% names(df)) df[[col]] <- as.character(df[[col]])
  }
  df
}

# DATA_NO = this Supplementary Data file's number; .t() builds the per-sheet
# title row.
DATA_NO <- 2L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)

# 00_Index — file description, provenance, sheet index.
readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Outcome variable", "Analytic sample",
           "Estimator", "Significance stars"),
  value = c("2 - Educational attainment",
            "Educational attainment (gene-environment interaction with the education PGI)",
            "Years of education standardized within birth year via a Gaussian kernel (cohort-relative standing)",
            "Full-family random-effects sample, N = 12,763 (one random intercept per family; TwinLife twins/children excluded)",
            "Linear and generalized additive mixed models (mgcv::bam, REML); the specification check is a logistic-transition NLS on the one-adult-per-family sample",
            "* p<.05, ** p<.01, *** p<.001"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_Linear_birth_year", "02_Step_reunification", "03_GAM_region",
            "04_Linear_gender", "05_GAM_gender", "06_GAM_region_gender",
            "07_Overall_RxG", "08_Logistic_check"),
  analysis_id = c("A1", "A2", "A3", "A4", "A5", "A6", "A0g", "A-spec"),
  description = c("Linear PGI x birth-year x region interaction",
                  "Step-function PGI x reunification x region interaction",
                  "Region-varying nonlinear (GAM) PGI cohort trajectory",
                  "Four-way PGI x birth-year x region x gender interaction",
                  "Per-gender PGI cohort trajectory, regions pooled (GAM)",
                  "Region x gender PGI cohort trajectory (GAM)",
                  "Overall additive region x gender PGI contrast",
                  "Logistic-transition specification check"),
  model = c("Linear mixed model", "Step-function mixed model", "Varying-coefficient GAM",
            "Linear four-way mixed model", "Per-gender GAM", "Region x gender GAM",
            "Linear contrast", "Logistic NLS (specification check)"),
  stringsAsFactors = FALSE)
sheet_README <- multiblock_sheet(
  list(header = "File description", df = readme_file,
       note = paste("Section A reports how the education polygenic index (PGI) predicts educational",
                    "attainment across birth cohorts, East/West Germany, and gender. Region and gender",
                    "are contrast-coded (East = +0.5 / West = -0.5; Female = +0.5 / Male = -0.5);",
                    "coefficients are per 1 SD of the PGI or birth year. Gender (sheets 4-8) is a main",
                    "moderator, corresponding to Section C of the preregistration plus a per-gender GAM.")),
  list(header = "Provenance", df = run_metadata,
       note = "Publication-safe run metadata (no local paths, usernames, or volume names)."),
  list(header = "Sheet index", df = readme_index,
       note = "Internal analysis IDs (A1-A6, A0g, A-spec) are kept inside each sheet for cross-reference."),
  title = sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))

# Descriptives (attrition, by-study, study x region, PGI intercorrelations,
# leave-one-study-out sample composition) live in the 00_Descriptives module,
# not in this analysis workbook.

# Per-analysis model-fit rows, split out of the section-wide A_ModelFit.
fit_A1 <- mc$A_ModelFit[mc$A_ModelFit$model == "A1_RE", , drop = FALSE]
fit_A2 <- mc$A_ModelFit[mc$A_ModelFit$model == "A2_RE", , drop = FALSE]
fit_A3 <- mc$A_ModelFit[grepl("^A3", mc$A_ModelFit$model), , drop = FALSE]

sheet_A1 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = A1_Headline,
       note = paste("A1: linear mixed model of cohort-relative years of education on the education PGI,",
                    "birth year, region and gender, with all two- and three-way interactions and a",
                    "per-family random intercept. Focal term: PGI x birth year x region - whether the",
                    "PGI-education association changes across cohorts differently in East vs West.")),
  list(header = "Model coefficients", df = mc$A1_Coefficients,
       note = "Parametric coefficients (per 1 SD of the PGI/birth year); the family random intercept is not shown."),
  list(header = "Model fit", df = fit_A1, note = "Fit statistics (sample size, AIC, BIC, log-likelihood, deviance explained)."),
  list(header = "Leave-one-study-out sensitivity", df = mc$A1_LOO,
       note = "Focal coefficient re-estimated with each study held out in turn; change_vs_full is the shift from the full-sample estimate."),
  title = .t(1, "Linear model results for educational attainment across birth years and region"))

sheet_A2 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = A2_Headline,
       note = paste("A2: as A1 but replacing continuous birth year with a reunification step (born >= 1975,",
                    "i.e. turned 15 in/after 1990). Focal term: PGI x post-reunification x region - whether",
                    "the PGI-education association shifts after reunification differently in East vs West.")),
  list(header = "Model coefficients", df = mc$A2_Coefficients,
       note = "Parametric coefficients; region/gender contrast-coded; family random intercept not shown."),
  list(header = "Model fit", df = fit_A2, note = "Fit statistics for the A2 model."),
  list(header = "Leave-one-study-out sensitivity", df = mc$A2_LOO,
       note = "Focal coefficient re-estimated with each study held out in turn."),
  title = .t(2, "Step-function (reunification) results for educational attainment"))

# East-West difference smooth (the planned contrast). Computed here so
# both the A3 sheet and the A3_diff PDF page reuse the same object.
attr(dat_cluster, "by_mean") <- mc$by_mean
diff_A <- diff_smooth_simul(m_A3_nonlin, dat_cluster)
A3_diff_intervals <- diff_smooth_significant_intervals(
  diff_A, contrast = "East - West PGI->education slope")

# (c) contrasts for the gender GAMs: A5 female-male (regions pooled);
# A6 female-male within each region.
A5_diff_intervals <- diff_smooth_significant_intervals(
  diff_smooth_pooled_gender_simul(m_A5_nonlin, dat_cluster),
  contrast = "Female - Male PGI->education slope",
  dir_pos = "Female > Male", dir_neg = "Female < Male")
A6_diff_intervals <- rbind(
  diff_smooth_significant_intervals(
    diff_smooth_gender_simul(m_A6_nonlin, dat_cluster, east_west = "West"),
    contrast = "Female - Male PGI->education slope (West)",
    dir_pos = "Female > Male", dir_neg = "Female < Male"),
  diff_smooth_significant_intervals(
    diff_smooth_gender_simul(m_A6_nonlin, dat_cluster, east_west = "East"),
    contrast = "Female - Male PGI->education slope (East)",
    dir_pos = "Female > Male", dir_neg = "Female < Male"))

.note_smoothF <- paste("summary.gam F-tests for the region-specific PGI x birth-year smooths. These",
  "overstate curvature because the smooth's linear part is confounded with the parametric",
  "PGI x birth-year slope; use the clean nested-LRT curvature test above for per-cell nonlinearity.")

# A3 reports the plan's three A3 evaluation levels in order:
#   (a) omnibus - is any (non-linear) cohort variation needed?
#   (b) per-region curvature (clean nested LRT)
#   (c) East-West contrast (difference smooth, simultaneous CI)
# The LINEAR cohort trend is tested in A1 (Wald), not re-tested here.
sheet_A3 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = A3_Headline,
       note = paste("A3 tests whether the PGI-education slope departs from linear across birth",
                    "cohorts, separately by region (family random intercept). The linear PGI x",
                    "birth-year(x region) trend itself is tested in A1 (sheet 01, Wald); A3 adds",
                    "only the non-linear departure.")),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$A3_NestedLRT,
       note = paste("Baseline = PGI main effect, constant slope (the plan's model WITHOUT",
                    "cohort-varying terms); linear adds the cohort slope; nonlinear adds the",
                    "smooth. AIC/BIC/log-likelihood compare the three; the LRT column carries",
                    "only the linear -> nonlinear (non-linearity) test.")),
  # (b) Per-region curvature --------------------------------------------
  list(header = "(b) Per-region curvature test (clean nested LRT)", df = mc$A3_cell_LRT,
       note = paste("Adds only that region's PGI x birth-year smooth to the linear model and",
                    "LRT-tests it - the clean per-region non-linearity test (supersedes the",
                    "confounded summary.gam F; see Plan_deviations.md section 6). Curves in the",
                    "PDF draw the smooth only where this is significant.")),
  # (c) East-West contrast ----------------------------------------------
  list(header = "(c) East-West difference smooth (planned contrast)", df = A3_diff_intervals,
       note = paste("Birth-year ranges where the East-West difference in the PGI-education slope",
                    "has a 95% simultaneous CI excluding zero (see the difference-smooth figure",
                    "in the PDF).")),
  # Robustness + diagnostics --------------------------------------------
  list(header = "Robustness: leave-one-study-out", df = mc$A3_LOO,
       note = "Linear-vs-nonlinear LRT re-run with each study held out in turn."),
  list(header = "Diagnostics: k-check", df = mc$A3_K_Index,
       note = "Basis-dimension adequacy for each smooth (low p / k-index < 1 suggests k too small)."),
  title = .t(3, "Region-varying nonlinear (GAM) results for educational attainment"))

# ---- Gender sheets ----
fit_A4 <- mc$A4_ModelFit
sheet_A4 <- multiblock_sheet(
  list(header = "Focal analysis summary", df = A4_Headline,
       note = paste("A4: linear four-way PGI x birth-year x region x gender mixed model (family random",
                    "intercept). Focal term: the four-way interaction - whether the cohort x region",
                    "pattern of the PGI-education association differs between women and men.")),
  list(header = "Model coefficients", df = mc$A4_Coefficients,
       note = "Full factorial parametric coefficients (per 1 SD of the PGI/birth year)."),
  list(header = "Model fit", df = fit_A4, note = "Fit statistics for the A4 model."),
  list(header = "Leave-one-study-out sensitivity", df = mc$A4_LOO,
       note = "Four-way focal coefficient re-estimated with each study held out in turn."),
  list(header = "Leave-one-study-out sensitivity (pooled PGI x gender)", df = nz(mc$A4_LOO_gender),
       note = paste("Pooled PGI x gender interaction (not the four-way term above) re-estimated with",
                    "each study held out; carries estimate, SE, 95% CI, exact p, supported birth-year",
                    "range and change vs full. Additional sensitivity - see Plan_deviations.md section 14.")),
  title = .t(4, "Linear gender-interaction results for educational attainment"))

sheet_A5 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = A5_Headline,
       note = paste("A5: per-gender GAM (regions pooled, family random intercept) of how the",
                    "PGI-education slope changes across birth cohorts for women vs men. The linear",
                    "gender x cohort trend is tested in A4 (Wald); A5 adds only the non-linear departure.")),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$A5_NestedLRT,
       note = paste("Baseline = PGI main effect (constant slope); linear adds the per-gender cohort",
                    "slopes; nonlinear adds the smooths. The LRT column carries only the",
                    "linear -> nonlinear (non-linearity) test.")),
  # (b) Per-gender curvature --------------------------------------------
  list(header = "(b) Per-gender curvature test (clean nested LRT)", df = mc$A5_cell_LRT,
       note = paste("Adds only that gender's PGI x birth-year smooth to the linear model - the clean",
                    "per-gender non-linearity test (supersedes the confounded summary.gam F; see",
                    "Plan_deviations.md section 6).")),
  # (c) Female-Male contrast --------------------------------------------
  list(header = "(c) Female-Male difference smooth (planned contrast)", df = A5_diff_intervals,
       note = paste("Birth-year ranges where the female-male difference in the PGI-education slope",
                    "has a 95% simultaneous CI excluding zero (see the difference-smooth figure).")),
  # (d) Leave-one-study-out --------------------------------------------
  list(header = "(d) Leave-one-study-out (omnibus + per-gender curvature)", df = nz(mc$A5_LOO),
       note = paste("Per fold: the omnibus linear-vs-nonlinear LRT plus the female- and male-specific",
                    "curvature LRTs, with the supported birth-year range. Checks the (null) per-gender",
                    "cohort nonlinearity is not created or masked by any one study. Additional",
                    "sensitivity - see Plan_deviations.md section 14.")),
  title = .t(5, "Per-gender nonlinear (GAM) results for educational attainment"))

sheet_A6 <- multiblock_sheet(
  # (a) Omnibus ----------------------------------------------------------
  list(header = "(a) Omnibus nonlinearity test (linear vs nonlinear GAM)", df = A6_Headline,
       note = paste("A6: region x gender GAM (family random intercept) with a PGI-education cohort",
                    "smooth in each of the four region x gender cells. The linear four-way trend is",
                    "tested in A4 (Wald); A6 adds only the non-linear departure.")),
  list(header = "(a) Model-fit ladder (information criteria)", df = mc$A6_NestedLRT,
       note = paste("Baseline = PGI main effects (constant slopes); linear adds the four cell cohort",
                    "slopes; nonlinear adds the smooths. The LRT column carries only the",
                    "linear -> nonlinear (non-linearity) test.")),
  # (b) Per-cell curvature ----------------------------------------------
  list(header = "(b) Per-cell curvature test (clean nested LRT)", df = mc$A6_cell_LRT,
       note = paste("Adds only that region x gender cell's PGI smooth to the linear model - the clean",
                    "per-cell non-linearity test (supersedes the confounded summary.gam F; see",
                    "Plan_deviations.md section 6).")),
  # (c) Female-Male contrast within region ------------------------------
  list(header = "(c) Female-Male difference smooth within region (planned contrast)", df = A6_diff_intervals,
       note = paste("Birth-year ranges where the female-male difference in the PGI-education slope",
                    "has a 95% simultaneous CI excluding zero, within West and within East.")),
  # Descriptive + robustness + diagnostics ------------------------------
  list(header = "Overall additive region x gender PGI slopes", df = mc$A6_OverallSlopes,
       note = "PGI-education slope per region x gender cell from the additive (no-interaction) model."),
  list(header = "Robustness: leave-one-study-out", df = mc$A6_LOO,
       note = "Linear-vs-nonlinear LRT re-run with each study held out in turn."),
  list(header = "Diagnostics: k-check", df = mc$A6_K_Index, note = "Basis-dimension adequacy for each smooth."),
  title = .t(6, "Region x gender nonlinear (GAM) results for educational attainment"))

.a0g_blocks <- list(
  list(header = "Focal analysis summary", df = A0g_Headline,
       note = paste("A0g: overall region x gender PGI-education contrast (no birth-year x PGI term; birth",
                    "year kept as a covariate), family random intercept. Focal term: PGI x region x gender",
                    "- whether the gender gap in the PGI-education slope differs by region.")),
  list(header = "Model coefficients (full sample)", df = mc$A0g$coefs,
       note = "Parametric coefficients on the full-family sample."),
  list(header = "Per-cell region x gender PGI slopes (full sample)", df = mc$A0g$cell_slopes,
       note = "PGI-education slope in each region x gender cell, with 95% CI and p."),
  list(header = "Model fit (full sample)", df = mc$A0g$fit, note = "Fit statistics for the full-sample A0g model."))
if (!is.null(mc$A0g_pre)) {
  .a0g_blocks <- c(.a0g_blocks, list(
    list(header = "Model coefficients (pre-reunification subsample)", df = mc$A0g_pre$coefs,
         note = "As above, restricted to the pre-reunification (born < 1975) subsample."),
    list(header = "Per-cell PGI slopes (pre-reunification)", df = mc$A0g_pre$cell_slopes,
         note = "Per-cell PGI-education slopes in the pre-reunification subsample."),
    list(header = "Model fit (pre-reunification)", df = mc$A0g_pre$fit, note = "Fit statistics for the pre-reunification A0g model.")))
}
sheet_A0g <- do.call(multiblock_sheet, c(.a0g_blocks,
  list(title = .t(7, "Overall region x gender contrast for educational attainment"))))

sheet_Aspec <- multiblock_sheet(
  list(header = "Specification-check conclusion", df = data.frame(conclusion = mc$Aspec_note, stringsAsFactors = FALSE),
       note = paste("A-spec: the preregistered logistic-transition specification (Section C, C2) is demoted",
                    "to a convergence/fit check - it does not converge or improve on the linear model, so it",
                    "is reported as a diagnostic only (one-adult-per-family sample, no random effect).")),
  list(header = "Convergence and fit diagnostic", df = mc$Aspec_Diagnostic,
       note = "Convergence status, AIC/BIC, U - L (transition amplitude) and boundary LRT for each attempted specification."),
  title = .t(8, "Logistic-transition specification check for educational attainment"))

sheets <- list(
  "00_Index"            = sheet_README,
  "01_Linear_birth_year" = sheet_A1,
  "02_Step_reunification"= sheet_A2,
  "03_GAM_region"        = sheet_A3,
  "04_Linear_gender"     = sheet_A4,
  "05_GAM_gender"        = sheet_A5,
  "06_GAM_region_gender" = sheet_A6,
  "07_Overall_RxG"       = sheet_A0g,
  "08_Logistic_check"    = sheet_Aspec
)

# Final pass: human-readable terms, publication column names, significance stars.
sheets <- humanize_sheets(sheets)

OUT_XLSX <- file.path(OUT_DIR_A, "Supplementary_Data_2_Education_results.xlsx")
write_section_xlsx("A", sheets, OUT_XLSX)
cat(sprintf("A_report: xlsx written to %s (%d sheets)\n", OUT_XLSX, length(sheets)))

# ---- Plots ----
# The PDFs are the expensive stage: the workbook above is already written, and
# everything below is prediction plus ggplot rendering. A wording-only change to
# the workbook does not need them, so REPORT_SKIP_PLOTS=1 stops here. The xlsx and
# the Overview markdown are unaffected either way.
if (identical(Sys.getenv("REPORT_SKIP_PLOTS"), "1")) {
  cat("A_report: REPORT_SKIP_PLOTS=1 -- skipping the PDF stage; xlsx and overview still written.\n")
} else {

cat("A_report: building plots...\n")

# Kernel-smoothed individual education by birth year, by region
own_edu_moments_reg_re <- dplyr::bind_rows(
  kernel_moments(dat_cluster$education[dat_cluster$east_west == "West"],
                 dat_cluster$birth_year[dat_cluster$east_west == "West"],
                 bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "West"),
  kernel_moments(dat_cluster$education[dat_cluster$east_west == "East"],
                 dat_cluster$birth_year[dat_cluster$east_west == "East"],
                 bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "East")
)
p_kernel <- ggplot(own_edu_moments_reg_re,
                   aes(x = t, colour = Region, fill = Region)) +
  geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
  geom_line(aes(y = mu), linewidth = 0.9) +
  geom_point(aes(y = mu, size = n_raw), alpha = 0.30) +
  scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_fill_manual  (values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
  scale_size_continuous(range = c(0.3, 3), guide = "none") +
  labs(x = "Birth year", y = "Years of education",
       title = "Kernel-smoothed individual education by birth year",
       subtitle = sprintf("Bandwidth = %d years, shaded = +/- 1 SD; one smooth per region.", KERNEL_BW)) +
  .theme

p_cohort_cov <- plot_cohort_coverage(dat_cluster,
                                     title = "Birth-year distribution by region")

# Attrition (core, primary track)
p_attrition <- plot_attrition_flow(mc$attrition_core_primary,
                                   title = "Sample attrition — Section A (RE primary track)",
                                   y_max = 12500)

# A1 slopes — linear cohort-varying slope by region, with non-parametric overlay
attr(dat_cluster, "by_mean") <- mc$by_mean
curves_A1 <- dplyr::bind_rows(
  extract_pgi_slope(m_A1_re, dat_cluster, east_west = "West",
                    pgi_var = "PGI_Edu_z"),
  extract_pgi_slope(m_A1_re, dat_cluster, east_west = "East",
                    pgi_var = "PGI_Edu_z"))
np_slopes_A <- compute_binned_slopes(
  dat_cluster, outcome = "edu_z_kernel", pgi = "PGI_Edu_z",
  group_vars = "east_west")
p_A1 <- plot_combined_east_west(
  gam_curves = curves_A1, np_slopes = np_slopes_A,
  title = "A1: Linear PGI slope x Birth Year x Region",
  ylab = "Standardized beta (SD education per 1 SD PGI)",
  hist_dat = dat_cluster, hist_alpha = 0.10, hist_bins = 50)

# A2 effects — pre/post reunification slopes by region
a2_eff <- extract_a2_pgi_effects(m_A2_re, dat_cluster)
p_A2 <- ggplot(a2_eff,
               aes(x = period, y = slope, colour = east_west, group = east_west)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.08, linewidth = 0.7) +
  geom_line(linewidth = 1) + geom_point(size = 3.5) +
  .scale_col +
  scale_x_discrete(labels = c(pre = "Pre-reunification", post = "Post-reunification")) +
  labs(x = "Period", y = "Standardized beta (SD education per 1 SD PGI)",
       title = "A2: PGI-Education by reunification period",
       colour = "Region") +
  .theme

# A3 slopes + East–West difference smooth. Display curves draw the nonlinear
# smooth only for regions with significant curvature (per-region nested LRT);
# linear slope otherwise.
curves_A3 <- mc$A3_display_curves
p_A3 <- plot_combined_east_west(
  gam_curves = curves_A3, np_slopes = np_slopes_A,
  title = "A3: Cohort-varying PGI slope by Region (smooth shown only where curvature is significant)",
  ylab = "Standardized beta (SD education per 1 SD PGI)",
  hist_dat = dat_cluster, hist_alpha = 0.10, hist_bins = 50)

# diff_A already computed above (reused for both the 04_A3 table and
# this PDF page) so the simulation runs once.
p_A3_diff <- plot_diff_smooth(
  diff_A,
  title = "A3 Education: East-West standardized PGI difference",
  hist_dat = dat_cluster, hist_alpha = 0.10, hist_bins = 50)

# ---- Gender plots (A5 per-gender pooled; A6 Region x Gender 4-cell + diffs) ----
np_slopes_gender <- tryCatch(compute_binned_slopes(
  dat_cluster, outcome = "edu_z_kernel", pgi = "PGI_Edu_z",
  group_vars = c("east_west", "gender")), error = function(e) NULL)

# A5: Female vs Male PGI->education slope over birth year (regions pooled).
p_A5 <- tryCatch(
  ggplot(mc$A5_curves, aes(x = birth_year, y = slope, colour = gender, fill = gender,
                           linetype = gender)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = c(female = "#7B3294", male = "#008837")) +
    scale_fill_manual(  values = c(female = "#7B3294", male = "#008837")) +
    scale_linetype_manual(values = c(female = LT_FEMALE, male = LT_MALE)) +
    labs(x = "Birth year", y = "Standardized beta (SD education per 1 SD PGI)",
         title = "A5: Per-gender PGI-education slope over birth year (regions pooled; smooth shown only where curvature is significant)",
         colour = "Gender", fill = "Gender", linetype = "Gender") + .theme,
  error = function(e) { cat("p_A5 failed:", conditionMessage(e), "\n"); NULL })

# A6: 4-cell Region x Gender cohort-varying slopes (display curves: nonlinear
# smooth only for cells with significant curvature, linear otherwise) +
# Male-Female diff smooths.
curves_A6 <- mc$A6_display_curves
p_A6 <- if (!is.null(curves_A6)) tryCatch(plot_combined_east_west_gender(
  gam_curves = curves_A6, np_slopes = np_slopes_gender,
  title = "A6: Cohort-varying PGI slope by Region x Gender (smooth shown only where curvature is significant)",
  ylab = "Standardized beta (SD education per 1 SD PGI)",
  hist_dat = dat_cluster, hist_alpha = 0.10, hist_bins = 50),
  error = function(e) { cat("p_A6 failed:", conditionMessage(e), "\n"); NULL }) else NULL
p_A6_diff_W <- tryCatch(plot_diff_smooth(
  diff_smooth_gender_simul(m_A6_nonlin, dat_cluster, east_west = "West"),
  title = "A6 West: Female-Male standardized PGI difference",
  hist_dat = dat_cluster[dat_cluster$east_west == "West", ], hist_alpha = 0.10),
  error = function(e) { cat("p_A6_diff_W failed:", conditionMessage(e), "\n"); NULL })
p_A6_diff_E <- tryCatch(plot_diff_smooth(
  diff_smooth_gender_simul(m_A6_nonlin, dat_cluster, east_west = "East"),
  title = "A6 East: Female-Male standardized PGI difference",
  hist_dat = dat_cluster[dat_cluster$east_west == "East", ], hist_alpha = 0.10),
  error = function(e) { cat("p_A6_diff_E failed:", conditionMessage(e), "\n"); NULL })

# ---- PDF ----
OUT_PDF <- file.path(OUT_DIR_A, "A_Education_Plots.pdf")
pdf(OUT_PDF, width = 10, height = 7)
.plots_main <- Filter(Negate(is.null), list(
  p_attrition, p_kernel, p_cohort_cov, p_A1, p_A2, p_A3, p_A3_diff,
  p_A5, p_A6, p_A6_diff_W, p_A6_diff_E))
for (p in .plots_main) {
  tryCatch(print(p), error = function(e) cat("Plot failed:", conditionMessage(e), "\n"))
}
dev.off()
cat(sprintf("A_report: pdf written to %s\n", OUT_PDF))

# ---- Extra LOO PDF: A3 per-fold cohort-varying slope curves ----
# One panel per leave-one-study-out fold (Full sample + drop each
# cohort), showing how the A3 cohort-varying PGI slope by region changes
# when each cohort is removed. Makes the A3 LOO LRT findings visual.
if (!is.null(mc$A3_LOO_curves) && nrow(mc$A3_LOO_curves) > 0) {
  loo_curves <- mc$A3_LOO_curves
  fold_levels <- if (is.factor(loo_curves$fold)) levels(loo_curves$fold)
                 else unique(loo_curves$fold)

  # Trim each fold's curves to that fold's own per-region birth-year
  # coverage (2nd-98th percentile), the same .trim_to_data_range approach
  # the main A3 plot uses — otherwise the curves run far into sparse
  # early/late birth years. Each fold's subsample is reconstructed from
  # dat_cluster by dropping the held-out cohort.
  .folds_map <- cohort_groups(dat_cluster)
  .trim_one_fold <- function(fold_label) {
    fc <- loo_curves[as.character(loo_curves$fold) == fold_label, , drop = FALSE]
    if (!nrow(fc)) return(fc)
    if (fold_label == "Full sample") {
      sub <- dat_cluster
    } else {
      drop_levs <- .folds_map[[sub("^-", "", fold_label)]]
      sub <- dat_cluster[!(as.character(dat_cluster$cohort) %in% drop_levs), , drop = FALSE]
    }
    .trim_to_data_range(fc, sub, group_vars = "east_west")
  }
  loo_curves <- dplyr::bind_rows(lapply(fold_levels, .trim_one_fold))
  loo_curves$fold <- factor(loo_curves$fold, levels = fold_levels)

  p_loo_A3 <- ggplot(loo_curves,
      aes(x = birth_year, y = slope, colour = east_west, fill = east_west)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ fold) +
    .scale_col + .scale_fill +
    labs(x = "Birth year",
         y = "Standardized beta (SD education per 1 SD PGI)",
         title = "A3 leave-one-study-out: cohort-varying PGI slope by region",
         subtitle = "Each panel refits the A3 GAM with one study held out (Full sample = no study dropped).",
         colour = "Region", fill = "Region") +
    .theme

  # Second LOO page: kernel-smoothed individual education by birth year,
  # recomputed per fold (descriptive; no model fit). Mirrors the main
  # Kernel_OwnEdu plot but faceted by fold so you can see how the raw
  # education-by-cohort pattern shifts when each cohort is dropped.
  .kern_loo_one <- function(fold_label) {
    if (fold_label == "Full sample") {
      sub <- dat_cluster
    } else {
      drop_levs <- .folds_map[[sub("^-", "", fold_label)]]
      sub <- dat_cluster[!(as.character(dat_cluster$cohort) %in% drop_levs), , drop = FALSE]
    }
    dplyr::bind_rows(
      kernel_moments(sub$education[sub$east_west == "West"],
                     sub$birth_year[sub$east_west == "West"],
                     bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "West", fold = fold_label),
      kernel_moments(sub$education[sub$east_west == "East"],
                     sub$birth_year[sub$east_west == "East"],
                     bandwidth = KERNEL_BW) |> dplyr::mutate(Region = "East", fold = fold_label)
    )
  }
  kern_loo <- dplyr::bind_rows(lapply(fold_levels, .kern_loo_one))
  kern_loo$fold <- factor(kern_loo$fold, levels = fold_levels)
  p_loo_kernel <- ggplot(kern_loo, aes(x = t, colour = Region, fill = Region)) +
    geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mu), linewidth = 0.8) +
    geom_point(aes(y = mu, size = n_raw), alpha = 0.30) +
    facet_wrap(~ fold) +
    scale_colour_manual(values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_fill_manual  (values = c(West = COL_KERNEL_WEST, East = COL_KERNEL_EAST)) +
    scale_size_continuous(range = c(0.2, 2.5), guide = "none") +
    labs(x = "Birth year", y = "Years of education",
         title = "Leave-one-study-out: kernel-smoothed individual education by birth year",
         subtitle = sprintf("Bandwidth = %d years, shaded = +/- 1 SD; one study held out per panel.", KERNEL_BW),
         colour = "Region", fill = "Region") +
    .theme

  OUT_LOO_PDF <- file.path(OUT_DIR_A, "A_Education_LOO_Plots.pdf")
  pdf(OUT_LOO_PDF, width = 11, height = 7)
  for (p in list(p_loo_A3, p_loo_kernel)) {
    tryCatch(print(p),
             error = function(e) cat("LOO plot failed:", conditionMessage(e), "\n"))
  }
  dev.off()
  cat(sprintf("A_report: LOO pdf written to %s (2 pages)\n", OUT_LOO_PDF))
}
}

# ---- Section overview markdown ----
# Uses the workflow-wide p convention (.fmt_p_vec: "<.001" or three decimals)
# so the overview and the workbook render the same p identically.
.fmt_p <- function(p) {
  if (length(p) == 0) return("NA")
  if (length(p) > 1)  p <- p[1]    # scalar-safe
  if (is.na(p)) return("NA")
  .fmt_p_vec(p)
}

a1_focal_est <- if (nrow(a1_row) > 0)
  sprintf("beta = %.4f (SE = %.4f)", a1_row$estimate, a1_row$std.error) else "n/a"
a1_focal_p   <- if (nrow(a1_row) > 0) .fmt_p(a1_row$p.value) else "NA"
a2_focal_est <- if (nrow(a2_row) > 0)
  sprintf("beta = %.4f (SE = %.4f)", a2_row$estimate, a2_row$std.error) else "n/a"
a2_focal_p   <- if (nrow(a2_row) > 0) .fmt_p(a2_row$p.value) else "NA"
a3_focal_est <- if (nrow(a3_lrt_row) > 0)
  sprintf("chi2 = %.2f (df = %.2f)", a3_lrt_row$LRT_chisq, a3_lrt_row$LRT_df) else "n/a"
a3_focal_p   <- if (nrow(a3_lrt_row) > 0) .fmt_p(a3_lrt_row$LRT_p) else "NA"

a3_conclusion <- unique(A3_GAM_Overview$nonlinearity_conclusion)
a3_conclusion <- if (length(a3_conclusion) == 1) a3_conclusion else
  paste0("Mixed across smooths: ", paste(a3_conclusion, collapse = "; "))

a1_concl <- if (nrow(a1_row) > 0 && !is.na(a1_row$p.value) && a1_row$p.value < 0.05) {
  "Focal interaction significant"
} else {
  "No focal interaction at alpha=.05"
}
a2_concl <- if (nrow(a2_row) > 0 && !is.na(a2_row$p.value) && a2_row$p.value < 0.05) {
  "Focal interaction significant"
} else {
  "No focal interaction at alpha=.05"
}

# Gender focal strings (A4/A6 = focal coefficient; A5/A6 nonlinearity = LRT).
a4_focal_est <- if (nrow(a4_row) > 0)
  sprintf("beta = %.4f (SE = %.4f)", a4_row$estimate, a4_row$std.error) else "n/a"
a4_focal_p   <- if (nrow(a4_row) > 0) .fmt_p(a4_row$p.value) else "NA"
a4_concl     <- if (nrow(a4_row) > 0 && !is.na(a4_row$p.value) && a4_row$p.value < 0.05)
  "Four-way interaction significant" else "No four-way interaction at alpha=.05"
a5_focal_est <- if (nrow(a5_lrt_row) > 0)
  sprintf("chi2 = %.2f (df = %.2f)", a5_lrt_row$LRT_chisq, a5_lrt_row$LRT_df) else "n/a"
a5_focal_p   <- if (nrow(a5_lrt_row) > 0) .fmt_p(a5_lrt_row$LRT_p) else "NA"
a5_concl     <- unique(A5_GAM_Overview$nonlinearity_conclusion)
a5_concl     <- if (length(a5_concl) == 1) a5_concl else paste(a5_concl, collapse = "; ")
a6_focal_est <- if (nrow(a6_lrt_row) > 0)
  sprintf("chi2 = %.2f (df = %.2f)", a6_lrt_row$LRT_chisq, a6_lrt_row$LRT_df) else "n/a"
a6_focal_p   <- if (nrow(a6_lrt_row) > 0) .fmt_p(a6_lrt_row$LRT_p) else "NA"
a6_concl     <- unique(A6_GAM_Overview$nonlinearity_conclusion)
a6_concl     <- if (length(a6_concl) == 1) a6_concl else paste(a6_concl, collapse = "; ")
a0g_focal_est <- if (nrow(a0g_row) > 0)
  sprintf("beta = %.4f (SE = %.4f)", a0g_row$estimate, a0g_row$std.error) else "n/a"
a0g_focal_p   <- if (nrow(a0g_row) > 0) .fmt_p(a0g_row$p.value) else "NA"
a0g_concl     <- if (nrow(a0g_row) > 0 && !is.na(a0g_row$p.value) && a0g_row$p.value < 0.05)
  "Region x Gender interaction significant" else "No Region x Gender interaction at alpha=.05"

analyses_tbl <- data.frame(
  Analysis = c("A1", "A2", "A3", "A4", "A5", "A6", "A0g"),
  Question = c("Does the PGI x birth-year association differ by region?",
               "Does the PGI association differ before vs. after the reunification-related cutoff, by region?",
               "Is a nonlinear birth-cohort pattern needed?",
               "Does the PGI x birth-year x region pattern differ by gender (four-way)?",
               "Does the PGI-education slope's birth-cohort trajectory differ for women vs. men (regions pooled)?",
               "Does the cohort-varying PGI slope differ across Region x Gender cells?",
               "Does the overall PGI-education slope's gender gap differ by region?"),
  Model    = c("Linear RE model", "Step-function RE model", "GAM RE model",
               "Linear four-way RE model", "Per-gender GAM RE model (regions pooled)",
               "Region x Gender GAM RE model", "Overall Region x Gender RE contrast"),
  stringsAsFactors = FALSE
)
findings_tbl <- data.frame(
  Analysis = c("A1", "A2", "A3", "A4", "A5", "A6", "A0g"),
  `Focal test`         = c("PGI x birth year x region",
                           "PGI x post-1975 x region",
                           "Nonlinear GAM vs. linear model",
                           "PGI x birth year x region x gender",
                           "Per-gender nonlinear GAM vs. linear",
                           "Region x Gender nonlinear GAM vs. linear",
                           "PGI x region x gender"),
  `Estimate / statistic` = c(a1_focal_est, a2_focal_est, a3_focal_est,
                             a4_focal_est, a5_focal_est, a6_focal_est, a0g_focal_est),
  p          = c(a1_focal_p, a2_focal_p, a3_focal_p,
                 a4_focal_p, a5_focal_p, a6_focal_p, a0g_focal_p),
  Conclusion = c(a1_concl, a2_concl, a3_conclusion,
                 a4_concl, a5_concl, a6_concl, a0g_concl),
  check.names = FALSE, stringsAsFactors = FALSE
)
outputs_tbl <- data.frame(
  Output = c("Excel workbook", "Figure PDF", "LOO figure PDF", "Model cache"),
  File   = c("Supplementary_Data_2_Education_results.xlsx",
             "A_Education_Plots.pdf",
             "A_Education_LOO_Plots.pdf",
             "A_Education_Models.rds"),
  stringsAsFactors = FALSE
)

interp <- paste0(
  "A1 and A2 test whether the PGI-education association differs across ",
  "birth cohorts and East/West Germany using a linear interaction and a ",
  "step-function around the reunification-related birth-year cutoff. A3 ",
  "asks whether a nonlinear birth-cohort pattern is needed; the GAM result ",
  "is interpreted only as evidence about nonlinearity (`",
  a3_conclusion, "`). The shape of the function is shown in the PDF ",
  "figures rather than summarized through selected birth-year anchors. ",
  "Each focal result is checked for robustness with leave-one-study-out ",
  "refits (BASE-II / SHIP / SOEP / TwinLife held out in turn); see the ",
  "LOO block at the bottom of each analysis sheet in the workbook. ",
  "Gender is reported here as a main moderator: A4 tests the four-way ",
  "PGI x birth-year x region x gender interaction; A5 traces the per-gender ",
  "PGI-education slope over birth year (regions pooled); A6 fits a ",
  "Region x Gender varying-coefficient GAM; and A0g reports the overall ",
  "additive Region x Gender PGI contrast. A4, A6, A0g and the A-spec ",
  "logistic-transition specification check correspond to the Section C ",
  "analyses of the preregistration (C1, C3, C0, C2), estimated here with ",
  "family random effects. Where the higher-order interactions are null, the ",
  "lower-order / additive gender moderation is the reported effect."
)

OUT_MD <- file.path(OUT_DIR_A, "A_Education_Overview.md")
generate_section_overview_md(
  section_id     = "A",
  section_title  = "Section A — Educational attainment (cohort, region, gender)",
  purpose        = paste0(
    "Section A tests whether the association between the education PGI ",
    "and educational attainment differs across birth cohorts, between ",
    "East and West Germany, and by gender. Analyses use the primary RE ",
    "sample (`dat_cluster`, N = ", N_a, ") with a per-family random ",
    "intercept for TwinLife families. Gender is a main moderator here ",
    "(A4/A5/A6/A0g): the Section C analyses of the preregistration estimated ",
    "with family random effects, plus a per-gender attainment GAM (A5)."),
  analyses       = analyses_tbl,
  findings       = findings_tbl,
  interpretation = interp,
  outputs        = outputs_tbl,
  path           = OUT_MD
)
cat(sprintf("A_report: overview written to %s\n", OUT_MD))

cat("\nA_report: DONE.\n")
