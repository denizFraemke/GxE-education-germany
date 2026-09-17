# helpers_descriptive_tables.R
# Building blocks for the publishable phenotype-descriptives workbook
# (00_Descriptives/). All outputs are AGGREGATE ONLY (counts, means, SDs,
# test p-values, effect sizes) — no individual-level rows are produced or
# stored.
#
# Sign conventions for the two-group comparisons:
#   Region difference  = East  - West
#   Gender difference  = Female - Male
# p-values are Welch's t-test (unequal variances); effect sizes are
# Hedges' g (small-sample bias-corrected Cohen's d, pooled SD). A
# comparison is only computed when BOTH groups have at least
# DESC_MIN_CELL_N finite observations; otherwise the test cells read
# "n/a" and an empty group's M(SD) cell reads "—".
#
# Requires DESC_MIN_CELL_N from 00_setup/constants.R.

# --------------------------------------------------------------------------
# fmt_msd — "mean (sd)" to 2 dp, or "—" if no finite values.
# --------------------------------------------------------------------------
fmt_msd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("—")          # em dash
  sprintf("%.2f (%.2f)", mean(x), stats::sd(x))
}

# --------------------------------------------------------------------------
# fmt_p — the workflow-wide p-value convention (.fmt_p_vec: "<.001" or three
# decimals), with an untestable cell rendered "n/a" rather than blank.
# --------------------------------------------------------------------------
fmt_p <- function(p) {
  if (length(p) != 1 || is.na(p)) return("n/a")
  .fmt_p_vec(p)
}

# --------------------------------------------------------------------------
# hedges_g — bias-corrected standardized mean difference (x1 - x2), pooled
# SD. Returns NA if either group has < 2 finite values or pooled SD is 0.
# --------------------------------------------------------------------------
hedges_g <- function(x1, x2) {
  x1 <- x1[is.finite(x1)]; x2 <- x2[is.finite(x2)]
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  d <- (mean(x1) - mean(x2)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)              # Hedges' small-sample correction
  d * J
}

# --------------------------------------------------------------------------
# welch_p — Welch two-sample t-test p-value (x1 vs x2), NA on failure.
# --------------------------------------------------------------------------
welch_p <- function(x1, x2) {
  x1 <- x1[is.finite(x1)]; x2 <- x2[is.finite(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  tryCatch(stats::t.test(x1, x2)$p.value, error = function(e) NA_real_)
}

# --------------------------------------------------------------------------
# group_cells — per-group N and M(SD) for one measure split by a two-level
# grouping vector. Used on the FULL-FAMILY sample (N / M(SD) reporting).
# --------------------------------------------------------------------------
# An absent group (N = 0, e.g. SHIP West) renders blank (NA / "") rather
# than "—", so single-region cohorts simply leave those cells empty.
group_cells <- function(x, g, lvl1, lvl2) {
  g  <- as.character(g)
  x1 <- x[which(g == lvl1 & is.finite(x))]
  x2 <- x[which(g == lvl2 & is.finite(x))]
  list(
    n1 = if (length(x1) == 0) NA_integer_ else length(x1),
    cell1 = if (length(x1) == 0) "" else fmt_msd(x1),
    n2 = if (length(x2) == 0) NA_integer_ else length(x2),
    cell2 = if (length(x2) == 0) "" else fmt_msd(x2)
  )
}

# --------------------------------------------------------------------------
# group_test — (x1 - x2) Welch p and Hedges' g for one measure split by a
# two-level grouping vector. Used on the DEDUP one-per-family sample (clean
# independence). Gated on min_n in both groups.
#   lvl1 : "minuend" (East, female); lvl2 : "subtrahend" (West, male)
# --------------------------------------------------------------------------
group_test <- function(x, g, lvl1, lvl2, min_n = DESC_MIN_CELL_N) {
  g  <- as.character(g)
  x1 <- x[which(g == lvl1 & is.finite(x))]
  x2 <- x[which(g == lvl2 & is.finite(x))]
  # A structurally absent group (one side N = 0, e.g. SHIP West) leaves the
  # comparison blank; "n/a" is reserved for a present-but-too-small group.
  if (length(x1) == 0 || length(x2) == 0) return(list(g = "", p = "", sig = ""))
  testable <- length(x1) >= min_n && length(x2) >= min_n
  if (!testable) return(list(g = "n/a", p = "n/a", sig = ""))
  pv  <- welch_p(x1, x2)
  sig <- if (is.na(pv)) "" else if (pv < .001) "***" else if (pv < .01) "**" else if (pv < .05) "*" else ""
  list(g = sprintf("%.2f", hedges_g(x1, x2)), p = fmt_p(pv), sig = sig)
}

# --------------------------------------------------------------------------
# by_group3 — the study's three birth-year bins. The 1975 boundary is the
# reunification cut (born >= 1975 = turned 15 in/after 1990), shared with
# the A/B analyses; the lower cut is 1950.
# --------------------------------------------------------------------------
by_group3 <- function(birth_year) {
  by <- suppressWarnings(as.numeric(birth_year))
  factor(
    ifelse(by < 1950, "< 1950",
           ifelse(by < 1975, "1950-1974", ">= 1975")),
    levels = c("< 1950", "1950-1974", ">= 1975")
  )
}

# Column order is fixed so the report-stage writers can place grouped
# headers positionally: overall | East | West | East-West Difference
# (Hedge's g, p) | Female | Male | Gender Difference (Hedge's g, p).
DESC_STAT_COLS <- c("N_overall", "MSD_overall",
                    "N_East", "MSD_East", "N_West", "MSD_West",
                    "g_EW", "p_EW", "sig_EW",
                    "N_Female", "MSD_Female", "N_Male", "MSD_Male",
                    "g_FM", "p_FM", "sig_FM")

# --------------------------------------------------------------------------
# desc_measure_row — one descriptive row for a single measure.
#   N and M(SD) (overall + per group) come from `full`  (full-family sample);
#   the East-West / Female-Male tests (Hedges' g, p) come from `dedup`
#   (one-per-family sample) for clean independence.
# Both frames must carry `east_west` (West/East) and `gender` (female/male).
# --------------------------------------------------------------------------
desc_measure_row <- function(full, dedup, label, col, min_n = DESC_MIN_CELL_N) {
  xf <- if (is.null(full[[col]]))  rep(NA_real_, nrow(full))  else suppressWarnings(as.numeric(full[[col]]))
  xd <- if (is.null(dedup[[col]])) rep(NA_real_, nrow(dedup)) else suppressWarnings(as.numeric(dedup[[col]]))
  n_all     <- sum(is.finite(xf))
  overall   <- if (n_all == 0) "" else fmt_msd(xf)
  n_overall <- if (n_all == 0) NA_integer_ else n_all
  ew_c <- group_cells(xf, full$east_west, "East",   "West")
  fm_c <- group_cells(xf, full$gender,    "female", "male")
  ew_t <- group_test(xd, dedup$east_west, "East",   "West", min_n)
  fm_t <- group_test(xd, dedup$gender,    "female", "male", min_n)
  data.frame(
    Measure     = label,
    N_overall   = n_overall, MSD_overall = overall,
    N_East      = ew_c$n1,  MSD_East    = ew_c$cell1,
    N_West      = ew_c$n2,  MSD_West    = ew_c$cell2,
    g_EW        = ew_t$g,   p_EW        = ew_t$p,   sig_EW = ew_t$sig,
    N_Female    = fm_c$n1,  MSD_Female  = fm_c$cell1,
    N_Male      = fm_c$n2,  MSD_Male    = fm_c$cell2,
    g_FM        = fm_t$g,   p_FM        = fm_t$p,   sig_FM = fm_t$sig,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# build_desc_table — overall descriptive table: one row per measure.
#   full / dedup : full-family (N/M(SD)) and dedup (tests) samples
#   measures     : named list  label -> column name
# --------------------------------------------------------------------------
build_desc_table <- function(full, dedup, measures, min_n = DESC_MIN_CELL_N) {
  rows <- lapply(names(measures), function(lab)
    desc_measure_row(full, dedup, lab, measures[[lab]], min_n))
  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# build_desc_table_by_group — descriptive table split by the three
# birth-year bins: rows = measure x bin, with a leading "Birth-year group"
# column. Both samples are split by the same bins. Empty/sparse cells
# degrade gracefully via desc_measure_row.
# --------------------------------------------------------------------------
build_desc_table_by_group <- function(full, dedup, measures, min_n = DESC_MIN_CELL_N) {
  bgf    <- by_group3(full$birth_year)
  bgd    <- by_group3(dedup$birth_year)
  groups <- levels(bgf)
  rows <- list()
  for (lab in names(measures)) {
    for (gg in groups) {
      fs   <- full[which(bgf == gg), , drop = FALSE]
      ds   <- dedup[which(bgd == gg), , drop = FALSE]
      r    <- desc_measure_row(fs, ds, lab, measures[[lab]], min_n)
      gcol <- data.frame(`Birth-year group` = gg,
                         check.names = FALSE, stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- cbind(gcol, r)
    }
  }
  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# build_coverage_table — per-cohort coverage so empty cells elsewhere are
# self-explanatory. N is the full-family count (matches the analyses);
# N_test is the dedup one-per-family count (the t-test sample).
# --------------------------------------------------------------------------
build_coverage_table <- function(full, dedup) {
  cohorts <- sort(as.character(unique(full$cohort)))
  rows <- lapply(cohorts, function(ch) {
    sub <- full[full$cohort == ch, , drop = FALSE]
    by  <- suppressWarnings(as.numeric(sub$birth_year))
    regs <- sort(unique(as.character(sub$east_west)))
    nmob <- if (is.null(sub$mobility)) 0L else sum(is.finite(suppressWarnings(as.numeric(sub$mobility))))
    data.frame(
      Cohort        = ch,
      N             = nrow(sub),
      N_test        = sum(dedup$cohort == ch),
      BY_min        = if (all(is.na(by))) NA_integer_ else floor(min(by, na.rm = TRUE)),
      BY_max        = if (all(is.na(by))) NA_integer_ else floor(max(by, na.rm = TRUE)),
      Regions       = paste(regs, collapse = "/"),
      Single_region = length(regs) == 1L,
      N_mobility    = nmob,
      Has_mobility  = nmob > 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ==========================================================================
# Report-stage writers — render the publishable two-tier grouped tables
# directly to an openxlsx worksheet (group-header row + sub-header row +
# data). The 14 stat columns (DESC_STAT_COLS) are written starting at the
# first data column; a leading label column ("Measure" for Table 1,
# "Birth-year group" for Table 2) sits to their left.
# ==========================================================================

# Group spans relative to the first stat column (offset 0 = first stat col).
.desc_group_spans <- function() list(
  list(label = "overall",              off = c(0, 1)),
  list(label = "East",                 off = c(2, 3)),
  list(label = "West",                 off = c(4, 5)),
  list(label = "East-West Difference", off = c(6, 8)),
  list(label = "Female",               off = c(9, 10)),
  list(label = "Male",                 off = c(11, 12)),
  list(label = "Gender Difference",    off = c(13, 15))
)
.desc_subheaders <- c("N", "M (SD)", "N", "M (SD)", "N", "M (SD)",
                      "Hedge's g", "p", "sig.", "N", "M (SD)", "N", "M (SD)",
                      "Hedge's g", "p", "sig.")

# Left margin: one empty column (A); the leading label column sits at B.
DESC_LEFT <- 2L

.desc_styles <- function() list(
  group = openxlsx::createStyle(textDecoration = "bold", halign = "center",
                                fontColour = "#1F4E79", fgFill = "#D9E1F2"),
  sub   = openxlsx::createStyle(textDecoration = "bold", halign = "center",
                                fgFill = "#EEEEEE"),
  title = openxlsx::createStyle(textDecoration = "bold", fgFill = "#F2F2F2"),
  lab   = openxlsx::createStyle(textDecoration = "bold")
)

# Write the two header rows (group + sub) starting at `row`, with the first
# stat column at `c0`. `label_header` goes in the label column (c0 - 1) on
# the sub-header row. Returns the first data row.
.desc_write_header <- function(wb, sheet, row, c0, label_header, st) {
  for (g in .desc_group_spans()) {
    c1 <- c0 + g$off[1]; c2 <- c0 + g$off[2]
    openxlsx::writeData(wb, sheet, g$label, startRow = row, startCol = c1, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = c1:c2, rows = row)
    openxlsx::addStyle(wb, sheet, st$group, rows = row, cols = c1:c2, gridExpand = TRUE)
  }
  sub_row <- row + 1L
  openxlsx::writeData(wb, sheet, label_header, startRow = sub_row, startCol = c0 - 1L,
                      colNames = FALSE)
  openxlsx::addStyle(wb, sheet, st$sub, rows = sub_row, cols = c0 - 1L)
  for (j in seq_along(.desc_subheaders)) {
    openxlsx::writeData(wb, sheet, .desc_subheaders[j], startRow = sub_row,
                        startCol = c0 + j - 1L, colNames = FALSE)
  }
  openxlsx::addStyle(wb, sheet, st$sub, rows = sub_row,
                     cols = c0:(c0 + length(.desc_subheaders) - 1L), gridExpand = TRUE)
  sub_row + 1L
}

# Table 1 (overall): one row per measure, label column = "Measure" (at B,
# i.e. one empty column A for margin); stat columns start at C.
write_desc_table1 <- function(wb, sheet, start_row, df) {
  st <- .desc_styles()
  c0 <- DESC_LEFT + 1L                        # stat cols start one right of the label
  data_row <- .desc_write_header(wb, sheet, start_row, c0, "Measure", st)
  openxlsx::writeData(wb, sheet, data.frame(df$Measure), startRow = data_row,
                      startCol = DESC_LEFT, colNames = FALSE)
  openxlsx::writeData(wb, sheet, df[, DESC_STAT_COLS], startRow = data_row,
                      startCol = c0, colNames = FALSE)
  data_row + nrow(df) + 1L                    # +1 blank spacer
}

# Table 2 (by birth-year group): shared header, then one bold section title
# per measure followed by its three birth-year rows. Label column =
# "Birth-year group". `measures_order` fixes the section order.
write_desc_table2 <- function(wb, sheet, start_row, df, measures_order) {
  st <- .desc_styles()
  c0 <- DESC_LEFT + 1L
  row <- .desc_write_header(wb, sheet, start_row, c0, "Birth-year group", st)
  ncol_total <- c0 + length(.desc_subheaders) - 1L
  for (i in seq_along(measures_order)) {
    m   <- measures_order[i]
    sub <- df[df$Measure == m, , drop = FALSE]
    if (!nrow(sub)) next
    openxlsx::writeData(wb, sheet, sprintf("%d. %s", i, m), startRow = row, startCol = DESC_LEFT,
                        colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = DESC_LEFT:ncol_total, rows = row)
    openxlsx::addStyle(wb, sheet, st$title, rows = row, cols = DESC_LEFT:ncol_total, gridExpand = TRUE)
    row <- row + 1L
    openxlsx::writeData(wb, sheet, data.frame(sub$`Birth-year group`), startRow = row,
                        startCol = DESC_LEFT, colNames = FALSE)
    openxlsx::writeData(wb, sheet, sub[, DESC_STAT_COLS], startRow = row,
                        startCol = c0, colNames = FALSE)
    row <- row + nrow(sub)
  }
  row + 1L
}

# --------------------------------------------------------------------------
# save_desc_workbook — saveWorkbook with the same locked-SMB-target
# overwrite fallback used by save_xlsx (helpers_reporting.R).
# --------------------------------------------------------------------------
save_desc_workbook <- function(wb, path, label = "00_Descriptives") {
  tryCatch({
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    path
  }, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("temporarily unavailable|Permission denied|cannot open|cannot create",
              msg, ignore.case = TRUE)) {
      fb <- sub("\\.xlsx$", sprintf("_%s.xlsx", format(Sys.time(), "%Y%m%d_%H%M%S")), path)
      cat(sprintf("\n*** WARNING: could not overwrite %s (%s)\n*** %s -> fallback %s\n",
                  path, msg, label, fb))
      openxlsx::saveWorkbook(wb, fb, overwrite = TRUE)
      return(fb)
    }
    warning(w)
  })
}
