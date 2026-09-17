# Descriptives_report.R
# Phenotype-descriptives module — report stage.
# Reads Descriptives_Tables.rds (aggregate tables) and assembles
# Supplementary_Data_1_Descriptives_results.xlsx (built directly with openxlsx
# so the tables carry two-tier grouped headers) plus Descriptives_Overview.md.
#
# Sheets:
#   00_Index                file description + provenance + sheet index
#   01_Coverage             per-study coverage + reading notes
#   02_Overall_and_by_year  Table 1 (overall) + Table 2 (by birth-year group)
#   03..06_<study>          same two tables, one sheet per study
#   07_Analysis_sample      attrition, PGI intercorrelations, leave-one-study-out
#                           sample composition
#
# Each Table 1/2 uses a two-row header: a merged group row (overall / East
# / West / East-West Difference / Female / Male / Gender Difference) over a
# sub-header row (N, M (SD), …, Hedge's g, p). Effect size (Hedge's g) is
# reported before the p-value.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
  library(openxlsx)
})

# ---- Resolve paths ----
.find_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  fa <- gsub("~\\+~", " ", fa, fixed = FALSE)
  if (length(fa) > 0 && nzchar(fa[1])) return(normalizePath(fa[1], mustWork = TRUE))
  stop("Cannot determine Descriptives_report.R location.")
}
SECTION_DIR  <- normalizePath(dirname(.find_this_file()), mustWork = TRUE)
WORKFLOW_DIR <- normalizePath(file.path(SECTION_DIR, ".."), mustWork = TRUE)
SCRIPT_DIR   <- WORKFLOW_DIR
rm(.find_this_file)

source(file.path(WORKFLOW_DIR, "00_setup", "run_context.R"), local = TRUE)
source(file.path(WORKFLOW_DIR, "00_setup", "constants.R"),   local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "GxE_Germany_analysis_helpers.R"), local = TRUE)
source(file.path(SCRIPT_DIR,  "R", "helpers_descriptive_tables.R"),   local = TRUE)

# ---- Load table cache ----
OUT_DIR  <- output_dir("00_Descriptives")
RDS_PATH <- file.path(OUT_DIR, "Descriptives_Tables.rds")
if (!file.exists(RDS_PATH)) {
  stop("Descriptives_Tables.rds not found at ", RDS_PATH,
       " — run Descriptives_analysis.R first.")
}
tb <- readRDS(RDS_PATH)
cat(sprintf("Descriptives_report: loaded %s (N = %d)\n", RDS_PATH, tb$N_total))

# ---- Shared note text (folded under the tables they describe, not as a
#      separate legend block) ----
# Each note is a character VECTOR: one element per line. write_note() renders
# each element on its own row so long notes read as a short paragraph rather
# than one endless cell.
.main_note <- c(
  sprintf(paste("Unadjusted descriptive comparisons. Each cell shows N and M(SD) on the",
                "full-family RE sample (N = %d)."), tb$N_total),
  sprintf(paste("The East-West and Female-Male columns give Hedge's g (bias-corrected Cohen's d)",
                "and a Welch t-test p on the one-adult-per-family subsample (N = %d); a group needs",
                ">= %d observations, else 'n/a'."), tb$N_test, tb$min_cell_n),
  "Stars: * p<.05, ** p<.01, *** p<.001. Region difference = East - West; gender difference = Female - Male. Blank cells = a structurally absent group.",
  paste("For the '% at ceiling/floor' rows the M(SD) cell is the percentage of the group at the top",
        "code (>= 18y, university) / bottom of the scale (<= 9y, at most lower-secondary), with its",
        "binomial SD."))
.note_overall <- c(
  paste("Overall region/gender contrasts mix studies of different selectivity (East is mostly the",
        "population sample SHIP plus the smaller high-SES BASE-II); the per-study sheets give the",
        "cleaner within-study contrasts."),
  paste("Raw years-of-education contrasts conflate the secular trend, so the analyses use",
        "edu_z_kernel (years standardized within birth year)."),
  .main_note)
.note_byyear <- c(
  "Birth-year groups: < 1950 / 1950-1974 / >= 1975 (the 1975 cut = turned 15 in/after 1990).",
  paste("Raw years-of-education contrasts across groups conflate the secular trend, which is why the",
        "analyses use edu_z_kernel (years standardized within birth year)."),
  paste("A within-study education SD (the 'Years of education' M(SD)) that is stable across birth-year",
        "groups, together with a stable '% at ceiling' share, indicates the cohort trend in PGI",
        "prediction on edu_z_kernel reflects a genuine change in association rather than compression",
        "of the outcome scale (Plan_deviations.md §2/§3)."),
  .main_note)

run_metadata <- build_run_metadata("00_Descriptives", tb$RUN_TS, SCRIPT_DIR)

# ---- Build the Supplementary Data 1 workbook ----
DATA_NO <- 1L
.t <- function(n, txt) sprintf("Supplementary Data %d, Sheet %d. %s", DATA_NO, n, txt)
wb <- openxlsx::createWorkbook()
LEFT <- 2L                                          # one empty column A (margin)
st_title <- openxlsx::createStyle(textDecoration = "bold", fontSize = 13)
st_head  <- openxlsx::createStyle(textDecoration = c("bold", "underline"),
                                  fontColour = "#1F4E79", fontSize = 11)
st_cols  <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
st_note  <- openxlsx::createStyle(fontColour = "#595959", textDecoration = "italic")  # no wrap -> compact rows

write_title <- function(sheet, txt) {
  openxlsx::writeData(wb, sheet, txt, startRow = 2L, startCol = LEFT, colNames = FALSE)
  openxlsx::addStyle(wb, sheet, st_title, rows = 2L, cols = LEFT)
}
# Labelled block: blue underlined heading, table, optional grey note. Returns next row.
write_block <- function(sheet, row, header, df, note = NULL) {
  openxlsx::writeData(wb, sheet, header, startRow = row, startCol = LEFT, colNames = FALSE)
  openxlsx::addStyle(wb, sheet, st_head, rows = row, cols = LEFT)
  openxlsx::writeData(wb, sheet, df, startRow = row + 1L, startCol = LEFT)
  openxlsx::addStyle(wb, sheet, st_cols, rows = row + 1L,
                     cols = LEFT:(LEFT + ncol(df) - 1L), gridExpand = TRUE)
  row <- row + nrow(df) + 2L
  if (!is.null(note)) {
    openxlsx::writeData(wb, sheet, paste("Note.", note), startRow = row, startCol = LEFT, colNames = FALSE)
    openxlsx::addStyle(wb, sheet, st_note, rows = row, cols = LEFT)
    row <- row + 1L
  }
  row + 1L
}
# Renders a note as one row per vector element (the first prefixed "Note."),
# so long notes break into readable lines instead of one endless cell.
# Returns the next free row (+1 blank spacer).
write_note <- function(sheet, row, note) {
  note <- as.character(note)
  note[1] <- paste("Note.", note[1])
  for (k in seq_along(note)) {
    openxlsx::writeData(wb, sheet, note[k], startRow = row + k - 1L, startCol = LEFT, colNames = FALSE)
    openxlsx::addStyle(wb, sheet, st_note, rows = row + k - 1L, cols = LEFT)
  }
  row + length(note) + 1L
}

# --- 00_Index ---
readme_file <- data.frame(
  item = c("Supplementary Data file", "Domain", "Outcomes", "Test sample", "Reporting sample", "Significance stars"),
  value = c("1 - Descriptive statistics",
            "Descriptive statistics for the GxE analytic sample (education and mobility)",
            "Years of education, edu_z_kernel, % at ceiling/floor, mobility, PGI-Education, parental education",
            "Tests (Welch t / Hedge's g) on the one-adult-per-family sample for independence",
            "N and M(SD) on the full-family random-effects sample (matches the A/B analyses)",
            "* p<.05, ** p<.01, *** p<.001"),
  stringsAsFactors = FALSE)
readme_index <- data.frame(
  sheet = c("01_Coverage", "02_Overall_and_by_year", "03_BASE_II", "04_SHIP", "05_SOEP", "06_TwinLife", "07_Analysis_sample"),
  contents = c("Per-study coverage and the legend / caveats",
               "Full-sample descriptives: overall (Table 1) and by birth-year group (Table 2)",
               "BASE-II descriptives", "SHIP descriptives", "SOEP descriptives", "TwinLife descriptives",
               "Attrition, PGI intercorrelations and leave-one-study-out sample composition"),
  stringsAsFactors = FALSE)
openxlsx::addWorksheet(wb, "00_Index")
write_title("00_Index", sprintf("Supplementary Data %d. Index - contents and methods", DATA_NO))
r <- 4L
r <- write_block("00_Index", r, "File description", readme_file,
  note = paste("Region difference = East - West; gender difference = Female - Male. N and M(SD) use the",
               "full-family sample; tests use the one-adult-per-family subsample (clean independence)."))
r <- write_block("00_Index", r, "Provenance", run_metadata, note = "Publication-safe run metadata.")
r <- write_block("00_Index", r, "Sheet index", readme_index)
openxlsx::setColWidths(wb, "00_Index", cols = 1, widths = 3)
openxlsx::setColWidths(wb, "00_Index", cols = 2, widths = 26)
openxlsx::setColWidths(wb, "00_Index", cols = 3:9, widths = 34)

# --- 01_Coverage ---
openxlsx::addWorksheet(wb, "01_Coverage")
write_title("01_Coverage", .t(1, "Study coverage and reading notes"))
r <- 4L
r <- write_block("01_Coverage", r, "Study coverage", humanize_df(setNames(tb$coverage, sub("^Cohort$", "Study", names(tb$coverage)))),
  note = paste("N = full-family sample size (the size used in the analyses); N_test = one-adult-per-family",
               "subsample used for the difference tests. Single-region studies (SHIP is East-only) leave",
               "their West and East-West cells blank; studies without parental education (SHIP) leave their",
               "mobility and parental rows blank."))
openxlsx::setColWidths(wb, "01_Coverage", cols = 1, widths = 3)
openxlsx::setColWidths(wb, "01_Coverage", cols = 2, widths = 16)
openxlsx::setColWidths(wb, "01_Coverage", cols = 3:10, widths = 13)

# --- 02_Overall_and_by_year ---
openxlsx::addWorksheet(wb, "02_Overall_and_by_year")
write_title("02_Overall_and_by_year", .t(2, "Full-sample descriptives, overall and by birth-year group"))
r <- 4L
openxlsx::writeData(wb, "02_Overall_and_by_year", "Table 1. Overall", startRow = r, startCol = LEFT, colNames = FALSE)
openxlsx::addStyle(wb, "02_Overall_and_by_year", st_head, rows = r, cols = LEFT)
r <- write_desc_table1(wb, "02_Overall_and_by_year", r + 1L, tb$main_overall)
r <- write_note("02_Overall_and_by_year", r, .note_overall)
openxlsx::writeData(wb, "02_Overall_and_by_year", "Table 2. By birth-year group", startRow = r, startCol = LEFT, colNames = FALSE)
openxlsx::addStyle(wb, "02_Overall_and_by_year", st_head, rows = r, cols = LEFT)
r <- write_desc_table2(wb, "02_Overall_and_by_year", r + 1L, tb$main_by_year, tb$measures)
r <- write_note("02_Overall_and_by_year", r, .note_byyear)
openxlsx::setColWidths(wb, "02_Overall_and_by_year", cols = LEFT, widths = 24)
openxlsx::setColWidths(wb, "02_Overall_and_by_year", cols = (LEFT + 1L):(LEFT + 15L), widths = 12)

# --- Per-study sheets ---
cohorts <- names(tb$per_cohort)
.coh_sheet <- c("BASE-II" = "03_BASE_II", "SHIP" = "04_SHIP", "SOEP" = "05_SOEP", "TwinLife" = "06_TwinLife")
.coh_no    <- c("BASE-II" = 3L, "SHIP" = 4L, "SOEP" = 5L, "TwinLife" = 6L)
for (ch in cohorts) {
  pc <- tb$per_cohort[[ch]]
  nm <- if (!is.na(.coh_sheet[ch])) .coh_sheet[[ch]] else substr(gsub("[^A-Za-z0-9]", "_", ch), 1, 28)
  openxlsx::addWorksheet(wb, nm)
  write_title(nm, .t(if (!is.na(.coh_no[ch])) .coh_no[[ch]] else 9L,
                     sprintf("%s descriptives", ch)))
  r <- 4L
  openxlsx::writeData(wb, nm, sprintf("Table 1. %s, overall (N = %d)", ch, pc$N),
                      startRow = r, startCol = LEFT, colNames = FALSE)
  openxlsx::addStyle(wb, nm, st_head, rows = r, cols = LEFT)
  r <- write_desc_table1(wb, nm, r + 1L, pc$overall)
  r <- write_note(nm, r, .main_note)
  if (isTRUE(pc$show_byyear)) {
    openxlsx::writeData(wb, nm, sprintf("Table 2. %s, by birth-year group", ch),
                        startRow = r, startCol = LEFT, colNames = FALSE)
    openxlsx::addStyle(wb, nm, st_head, rows = r, cols = LEFT)
    r <- write_desc_table2(wb, nm, r + 1L, pc$by_year, tb$measures)
    r <- write_note(nm, r, .note_byyear)
  } else {
    r <- write_note(nm, r, "By-birth-year table omitted: this study's birth-year span is too narrow for a meaningful split.")
  }
  openxlsx::setColWidths(wb, nm, cols = LEFT, widths = 24)
  openxlsx::setColWidths(wb, nm, cols = (LEFT + 1L):(LEFT + 15L), widths = 12)
}

# --- 07_Analysis_sample ---
openxlsx::addWorksheet(wb, "07_Analysis_sample")
write_title("07_Analysis_sample", .t(7, "Analysis-sample descriptives (full-family RE sample)"))
r <- 4L
r <- write_block("07_Analysis_sample", r, "Sample attrition (full-family RE track)", humanize_df(tb$attrition),
  note = "Exclusions applied to reach the primary random-effects analytic sample, by study.")
r <- write_block("07_Analysis_sample", r, "PGI intercorrelations (full-family sample)", humanize_df(tb$pgi_correlations),
  note = "Pearson correlations among the education, cognitive, non-cognitive and height PGIs.")
r <- write_block("07_Analysis_sample", r, "Leave-one-study-out sample composition", humanize_df(tb$loo_composition),
  note = "How the full-family sample's size and composition shift when each study is held out.")
openxlsx::setColWidths(wb, "07_Analysis_sample", cols = 1, widths = 3)
openxlsx::setColWidths(wb, "07_Analysis_sample", cols = 2:9, widths = 15)

XLSX_PATH <- file.path(OUT_DIR, "Supplementary_Data_1_Descriptives_results.xlsx")
saved <- save_desc_workbook(wb, XLSX_PATH)
cat(sprintf("Descriptives_report: wrote %s\n", saved))

# ---- Overview markdown ----
analyses_md <- data.frame(
  Table  = c("Table 1", "Table 2", "Per-study sheets", "Coverage"),
  Scope  = c("Full analytic sample, overall",
             "Full analytic sample x 3 birth-year bins (<1950 / 1950-1974 / >=1975)",
             "Each study (BASE-II, SHIP, SOEP, TwinLife): Tables 1 and 2",
             "Per-study N, birth-year range, regions, mobility availability"),
  Sample = rep(tb$sample_label, 4),
  stringsAsFactors = FALSE
)
interpretation_md <- paste(
  "Descriptive, unadjusted comparisons. Region difference = East - West and gender",
  "difference = Female - Male, with Welch t-test p-values and Hedge's g effect sizes",
  sprintf("computed only when both groups reach N >= %d.", tb$min_cell_n),
  "Overall region/gender contrasts mix studies of very different selectivity (East is",
  "mostly the population sample SHIP plus the smaller high-SES BASE-II), so the per-study",
  "sheets give the cleaner within-study read. Single-region (SHIP) and mobility-absent",
  "(SHIP) cells, and sparse study x birth-year cells, read '—' or 'n/a'. Raw",
  "years-of-education contrasts conflate the secular trend, which is why the analyses use",
  "edu_z_kernel."
)
outputs_md <- data.frame(
  Output = c("Workbook", "Table cache", "This overview"),
  File   = c("Supplementary_Data_1_Descriptives_results.xlsx",
             "Descriptives_Tables.rds", "Descriptives_Overview.md"),
  stringsAsFactors = FALSE
)
MD_PATH <- file.path(OUT_DIR, "Descriptives_Overview.md")
generate_section_overview_md(
  section_id     = "00_Descriptives",
  section_title  = "Descriptives — Educational attainment and mobility",
  purpose        = paste(
    "Publishable phenotype descriptives for the GxE analytic sample: years of education,",
    "the analysis outcome edu_z_kernel, educational mobility, PGI-Education and parental",
    "education, broken down by region and gender (overall and by birth-year group), for the",
    "full sample and each study."),
  analyses       = analyses_md,
  findings       = tb$coverage,
  interpretation = interpretation_md,
  outputs        = outputs_md,
  path           = MD_PATH
)
cat(sprintf("Descriptives_report: wrote %s\n", MD_PATH))
