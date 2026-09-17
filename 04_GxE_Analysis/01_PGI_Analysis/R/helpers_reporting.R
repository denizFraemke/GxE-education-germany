# helpers_reporting.R
# Workbook assembly and section-overview generation for the
# 01_PGI_Analysis outputs.
#
# GAM overview schema (no birth-year anchors):
#   columns: analysis_id, sample, N, nested_LRT_chi2, nested_LRT_df,
#            nested_LRT_p, nested_aic_delta, nested_bic_delta,
#            smooth_label, edf, ref_df, F, p, alpha,
#            nonlinearity_conclusion, interpretation
#   nonlinearity_conclusion ∈ {
#     "Nonlinearity supported",
#     "No evidence that nonlinear specification is needed",
#     "Inconclusive",
#     "Not tested"
#   }
# The conclusion only answers whether the nonlinear GAM is needed
# relative to the corresponding linear model; it does not comment on
# the substantive importance of the PGI association.

# --------------------------------------------------------------------------
# build_run_metadata — standard run-provenance table for a section.
# repo_dir is shQuote()d before being passed to git so a repo path
# containing spaces (e.g. ".../015_Data analysis/...") is not split when
# system2() builds the shell command for output capture.
# --------------------------------------------------------------------------

build_run_metadata <- function(section_id, run_ts, repo_dir, alpha = "0.05") {
  git_sha <- tryCatch(
    system2("git", c("-C", shQuote(repo_dir), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = TRUE),
    error = function(e) "unknown")
  git_sha <- paste(git_sha, collapse = " ")
  if (grepl("fatal|error|not a git", git_sha, ignore.case = TRUE)) git_sha <- "unknown"
  # Publication-safe metadata only: no local paths / volumes / usernames.
  data.frame(
    key   = c("Analysis run", "Software", "Git commit", "Repository",
              "Significance threshold (alpha)", "Section"),
    value = c(run_ts, R.version.string, git_sha,
              "https://github.com/denizFraemke/80_years_GxE_on_Education_in_Germany-analysis",
              alpha, section_id),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# diff_smooth_significant_intervals — summarize the planned East-West
# (or female-male) difference smooth without birth-year anchors.
#
# Input: the data frame returned by diff_smooth_simul() /
#        diff_smooth_gender_simul(), with columns birth_year, diff,
#        ci_lo, ci_hi (simultaneous CI bounds).
#
# Reports the birth-year RANGES where the simultaneous CI excludes zero
# (i.e. where the band no longer touches the y = 0 line), as contiguous
# intervals — not anchored point estimates. A row per contiguous
# significant interval, or a single "none" row if the band covers zero
# everywhere.
# --------------------------------------------------------------------------

diff_smooth_significant_intervals <- function(diff_df,
                                              contrast = "East - West",
                                              ci_label = "95% simultaneous CI",
                                              dir_pos = "East > West",
                                              dir_neg = "East < West") {
  empty_row <- function(intervals, direction) {
    data.frame(
      contrast               = contrast,
      ci                     = ci_label,
      excludes_zero          = if (intervals == "none") "no" else "yes",
      birth_year_interval    = intervals,
      direction              = direction,
      max_abs_difference     = round(max(abs(diff_df$diff), na.rm = TRUE), 4),
      stringsAsFactors       = FALSE
    )
  }

  ok <- !is.na(diff_df$ci_lo) & !is.na(diff_df$ci_hi)
  sig <- ok & (diff_df$ci_lo > 0 | diff_df$ci_hi < 0)
  if (!any(sig)) {
    return(empty_row("none", "n/a"))
  }

  # Contiguous runs of significant grid points.
  r        <- rle(sig)
  ends     <- cumsum(r$lengths)
  starts   <- ends - r$lengths + 1L
  sig_idx  <- which(r$values)

  rows <- lapply(sig_idx, function(i) {
    rng_rows <- starts[i]:ends[i]
    by_lo <- min(diff_df$birth_year[rng_rows])
    by_hi <- max(diff_df$birth_year[rng_rows])
    # Direction is consistent within a contiguous CI-excludes-zero run.
    dir <- if (all(diff_df$ci_lo[rng_rows] > 0, na.rm = TRUE)) {
      dir_pos
    } else if (all(diff_df$ci_hi[rng_rows] < 0, na.rm = TRUE)) {
      dir_neg
    } else {
      "mixed"
    }
    data.frame(
      contrast            = contrast,
      ci                  = ci_label,
      excludes_zero       = "yes",
      birth_year_interval = sprintf("%.0f-%.0f", by_lo, by_hi),
      direction           = dir,
      max_abs_difference  = round(max(abs(diff_df$diff), na.rm = TRUE), 4),
      stringsAsFactors    = FALSE
    )
  })
  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# save_xlsx — write a named list of data frames to one workbook, with
# overwrite-with-fallback handling for locked network-share targets.
# --------------------------------------------------------------------------

save_xlsx <- function(tables, path, label) {
  wb <- openxlsx::createWorkbook()
  for (sheet_name in names(tables)) {
    sn <- substr(sheet_name, 1, 31)   # Excel sheet-name limit
    openxlsx::addWorksheet(wb, sn)
    tbl <- tables[[sheet_name]]
    if (inherits(tbl, "multiblock_sheet")) {
      .write_multiblock_sheet(wb, sn, tbl)
    } else {
      if (is.null(tbl)) tbl <- data.frame()
      for (j in seq_along(tbl)) {
        if (is.list(tbl[[j]])) tbl[[j]] <- vapply(tbl[[j]], as.character, character(1))
      }
      openxlsx::writeData(wb, sn, tbl)
    }
  }
  tryCatch({
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    path
  }, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("temporarily unavailable|Permission denied|cannot open|cannot create",
              msg, ignore.case = TRUE)) {
      fallback <- sub("\\.xlsx$",
                      sprintf("_%s.xlsx", format(Sys.time(), "%Y%m%d_%H%M%S")),
                      path)
      cat(sprintf(
        "\n*** WARNING: could not overwrite %s (%s)\n*** %s xlsx -> fallback %s\n",
        path, msg, label, fallback))
      openxlsx::saveWorkbook(wb, fallback, overwrite = TRUE)
      return(fallback)
    }
    warning(w)
  })
}

# --------------------------------------------------------------------------
# multiblock_sheet — build a structured "sheet" object that save_xlsx can
# render as a single Excel sheet containing several stacked tables, each
# preceded by a bold header row and separated by a blank row. Use this
# when several related tables share the same audience but have
# heterogeneous shapes (e.g. the section overview = run metadata + 3
# headline rows; the A3 nonlinearity sheet = headline + GAM overview +
# nested LRT + smooths + k.check).
#
# Usage:
#   multiblock_sheet(
#     list(header = "=== Run metadata ===", df = run_metadata_df),
#     list(header = "=== Headlines ===",    df = headlines_df),
#     ...
#   )
# --------------------------------------------------------------------------

# A "sheet" = an optional publication title + a list of blocks. Each block is
# list(header = , df = , note = NULL). `title` (the "Supplementary Data N,
# Sheet M. ..." line) is carried as an attribute.
multiblock_sheet <- function(..., title = NULL) {
  blocks <- list(...)
  for (b in blocks) {
    stopifnot(is.list(b),
              "header" %in% names(b),
              "df" %in% names(b))
  }
  structure(blocks, class = "multiblock_sheet", title = title)
}

# Internal: render a multiblock_sheet to an openxlsx worksheet — one empty top
# row + one empty left column for margin; optional bold sheet title; each block
# = a blue underlined heading, the table (bold blue column-name row), and an
# optional grey italic "Note." beneath it.
.write_multiblock_sheet <- function(wb, sn, sheet) {
  col0 <- 2L                       # one empty left column (A)
  row  <- 2L                       # one empty top row (1)
  st_title <- openxlsx::createStyle(textDecoration = "bold", fontSize = 13)
  st_head  <- openxlsx::createStyle(textDecoration = c("bold", "underline"),
                                    fontColour = "#1F4E79", fontSize = 11)
  st_cols  <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                                    border = "bottom", borderColour = "#1F4E79")
  # Notes are not wrapped, so they don't inflate row height (the full text is
  # visible in the formula bar when the cell is selected).
  st_note  <- openxlsx::createStyle(fontColour = "#595959", textDecoration = "italic")
  ttl <- attr(sheet, "title")
  if (!is.null(ttl) && nzchar(ttl)) {
    openxlsx::writeData(wb, sn, ttl, startRow = row, startCol = col0, colNames = FALSE)
    openxlsx::addStyle(wb, sn, st_title, rows = row, cols = col0)
    row <- row + 2L
  }
  for (b in sheet) {
    openxlsx::writeData(wb, sn, b$header, startRow = row, startCol = col0, colNames = FALSE)
    openxlsx::addStyle(wb, sn, st_head, rows = row, cols = col0)
    row <- row + 1L
    df <- b$df
    if (is.null(df) || (is.data.frame(df) && !ncol(df))) {
      openxlsx::writeData(wb, sn, "(no rows)", startRow = row, startCol = col0, colNames = FALSE)
      row <- row + 1L
    } else {
      for (j in seq_along(df)) {
        if (is.list(df[[j]])) df[[j]] <- vapply(df[[j]], as.character, character(1))
      }
      openxlsx::writeData(wb, sn, df, startRow = row, startCol = col0)
      openxlsx::addStyle(wb, sn, st_cols, rows = row,
                         cols = col0:(col0 + ncol(df) - 1L), gridExpand = TRUE)
      row <- row + nrow(df) + 1L
    }
    if (!is.null(b$note) && nzchar(b$note)) {
      openxlsx::writeData(wb, sn, paste("Note.", b$note), startRow = row,
                          startCol = col0, colNames = FALSE)
      openxlsx::addStyle(wb, sn, st_note, rows = row, cols = col0)
      row <- row + 1L
    }
    row <- row + 1L                # blank row between blocks
  }
}

# --------------------------------------------------------------------------
# build_sheet_index — produce a one-row-per-sheet index for the workbook.
# --------------------------------------------------------------------------

build_sheet_index <- function(sheets) {
  shape <- function(t) {
    if (is.null(t)) return(c(rows = 0L, cols = 0L))
    if (inherits(t, "multiblock_sheet")) {
      # Sum rows across blocks (header + colnames + data + blank), max cols.
      rows <- sum(vapply(t, function(b) {
        if (is.null(b$df)) 2L
        else as.integer(nrow(b$df) + 3L)   # header + colnames + N + blank
      }, integer(1)))
      cols <- max(vapply(t, function(b) {
        if (is.null(b$df)) 1L else as.integer(ncol(b$df))
      }, integer(1)))
      return(c(rows = rows, cols = cols))
    }
    c(rows = as.integer(nrow(t)), cols = as.integer(ncol(t)))
  }
  shapes <- vapply(sheets, shape, integer(2))
  data.frame(
    order      = seq_along(sheets),
    sheet      = substr(names(sheets), 1, 31),
    n_rows     = shapes["rows", ],
    n_cols     = shapes["cols", ],
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}

# --------------------------------------------------------------------------
# build_headline_table — one publication-ready row per analysis.
# Schema:
#   analysis_id | outcome | sample | N | cohorts_included |
#   model_family | focal_term | estimate | SE_or_CI | p
# The workbook is the human-facing artefact, so provenance columns are not
# emitted; run provenance lives in the 00_Index sheet. The `status`,
# `output_source` and `code_source` arguments are accepted and ignored.
# --------------------------------------------------------------------------

build_headline_table <- function(analysis_id,
                                 status = NULL,
                                 outcome,
                                 sample_name,
                                 N,
                                 cohorts_included,
                                 model_family,
                                 focal_term,
                                 estimate,
                                 SE_or_CI,
                                 p,
                                 output_source = NULL,
                                 code_source = NULL) {
  data.frame(
    analysis_id      = analysis_id,
    outcome          = outcome,
    sample           = sample_name,
    N                = N,
    cohorts_included = paste(cohorts_included, collapse = ", "),
    model_family     = model_family,
    focal_term       = focal_term,
    estimate         = estimate,
    SE_or_CI         = SE_or_CI,
    p                = p,
    stringsAsFactors = FALSE
  )
}

# ==========================================================================
# Human-readable relabelling for the published workbooks.
# Applied as a final pass over the assembled `sheets` object (humanize_sheets)
# so all the upstream focal-term lookups keep working on the raw R names and
# only the OUTPUT is prettified: R term names -> readable labels, no
# snake_case, clear column headers.
# ==========================================================================

.HUMAN_TOKENS <- c(
  "PGI_Edu_z" = "PGI-Education", "PGI_Edu" = "PGI-Education",
  "PGI_Edu_raw" = "PGI-Education (raw)", "PGI_Edu_z_within" = "PGI-Education (within-region)",
  "PGI_Cog_z" = "PGI-Cognition", "PGI_Cog" = "PGI-Cognition",
  "PGI_nonCog_z" = "PGI-NonCog", "PGI_nonCog" = "PGI-NonCog",
  "PGI_Height_z" = "PGI-Height", "PGI_Height" = "PGI-Height",
  "PGI_W" = "PGI-Education (West)", "PGI_E" = "PGI-Education (East)",
  "PGI_WF" = "PGI-Education (West, female)", "PGI_WM" = "PGI-Education (West, male)",
  "PGI_EF" = "PGI-Education (East, female)", "PGI_EM" = "PGI-Education (East, male)",
  "PGI_F" = "PGI-Education (female)", "PGI_M" = "PGI-Education (male)",
  "BYc_z" = "Birth year", "BYc" = "Birth year", "birth_year" = "Birth year",
  "east_west_c" = "Region (East vs West)", "east_west" = "Region", "east" = "East",
  "gender_c" = "Gender (Female vs Male)", "gender" = "Gender", "gender_01" = "Male",
  "reunifpost" = "Post-reunification", "reunif" = "Reunification", "reunifpre" = "Pre-reunification",
  "parental_edu_z_kernel" = "Parental education", "ParEdu_z" = "Parental education",
  "ParEdu_W" = "Parental education (West)", "ParEdu_E" = "Parental education (East)",
  "edu_z_kernel" = "Years of Education (kernel std.)", "edu_std" = "Years of education (global std.)",
  "education" = "Years of education", "mobility" = "Educational mobility",
  "attainment" = "Educational attainment",
  "height_z" = "Height (z)", "height" = "Height",
  "(Intercept)" = "Intercept", "fid_re" = "Family (random effect)")

# Human-readable, publication-facing column names.
.HUMAN_COLS <- c(
  analysis_id = "Analysis", outcome = "Outcome", sample = "Sample", cohorts_included = "Studies",
  model_family = "Model", focal_term = "Focal term", estimate = "Estimate", SE_or_CI = "SE / 95% CI",
  p = "p", p.value = "p", p_raw = "p (raw)", p_BH = "p (BH-adjusted)",
  std.error = "SE", statistic = "t / z", conf.low = "CI lower", conf.high = "CI upper",
  term = "Term", model = "Model", smooth_term = "Smooth term", smooth_label = "Smooth",
  edf = "edf", ref_df = "Ref. df", "p-value" = "p", p_value = "p",
  nested_LRT_chi2 = "LRT chi-sq", nested_LRT_df = "LRT df", nested_LRT_p = "LRT p",
  nested_aic_delta = "AIC delta", nested_bic_delta = "BIC delta", LRT_chisq = "LRT chi-sq",
  LRT_df = "LRT df", LRT_p = "LRT p", LRT_p_boundary = "LRT p (boundary)",
  nonlinearity_conclusion = "Nonlinearity conclusion", interpretation = "Interpretation",
  cell = "Cell", cohort = "Study", nonlinear_sig = "Nonlinear (significant)", fold = "Fold",
  held_out = "Study held out", delta_vs_full = "Change vs full", group = "Group",
  east_west = "Region", gender = "Gender", slope = "Slope", se = "SE",
  ci_lo = "CI lower", ci_hi = "CI upper", dev_expl = "Deviance explained",
  logLik = "Log-likelihood", R2_or_devExpl = "R-sq / deviance explained",
  family = "Family", PGI = "PGI", estimator = "Estimator", variance = "Variance model",
  region = "Region", scope = "Scope", dataset = "Study", k_prime = "k'", k_index = "k-index",
  comparison = "Comparison", N = "N", n = "N")

.tok <- function(p) if (p %in% names(.HUMAN_TOKENS)) .HUMAN_TOKENS[[p]] else gsub("_", " ", p)

humanize_term <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    mapped <- vapply(parts, function(p) {
      p <- trimws(p)
      m <- regmatches(p, regexec("^s\\((.*)\\)$", p))[[1]]
      if (length(m) == 2) return(paste0("Smooth(", .tok(m[2]), ")"))
      .tok(p)
    }, character(1))
    paste(mapped, collapse = " × ")
  }, character(1), USE.NAMES = FALSE)
}

.tidy_value <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(s)
    if (s %in% names(.HUMAN_TOKENS)) return(.HUMAN_TOKENS[[s]])
    if (grepl("_", s) && !grepl("\\s", s)) return(gsub("_", " ", s))  # identifier-ish -> de-snake
    s
  }, character(1), USE.NAMES = FALSE)
}

.tidy_sample <- function(x) {
  x <- as.character(x)
  x <- gsub("dat_cluster_mob", "Full-family RE mobility sample", x, fixed = TRUE)
  x <- gsub("dat_cluster", "Full-family RE sample", x, fixed = TRUE)
  x <- gsub("dat_edu", "Dedup one-per-family sample", x, fixed = TRUE)
  x <- gsub("dat_mob", "Dedup mobility sample", x, fixed = TRUE)
  x
}

.human_col <- function(nm) {
  if (grepl("__star$", nm)) return("sig.")
  if (nm %in% names(.HUMAN_COLS)) return(.HUMAN_COLS[[nm]])
  out <- gsub("[_.]", " ", nm)
  paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
}

# Significance stars from a raw p-value column (* <.05, ** <.01, *** <.001).
.stars <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  vapply(p, function(x) if (is.na(x)) "" else if (x < .001) "***" else
    if (x < .01) "**" else if (x < .05) "*" else "", character(1))
}

# p-value formatting — the single convention used by every artefact in this
# workflow (workbooks, section-overview markdown, descriptive tables):
# "<.001" for p below .001, otherwise three decimals. NA in, NA out.
.fmt_p_vec <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(p) if (is.na(p)) NA_character_ else if (p < .001) "<.001" else sprintf("%.3f", p),
         character(1))
}
.is_int_col <- function(x) {
  x <- x[is.finite(x)]
  length(x) == 0 || all(abs(x - round(x)) < 1e-9)
}

humanize_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !ncol(df)) return(df)
  # Drop a constant "model" column (e.g. coefficient tables all "A4 RE") - the
  # model is already named in the sheet title / block heading.
  if ("model" %in% names(df)) {
    mv <- df$model[!is.na(df$model)]
    if (length(mv) > 0 && length(unique(mv)) == 1L) df$model <- NULL
  }
  p_like <- intersect(c("p", "p.value", "p_raw", "p_BH", "p-value", "p_value",
    "LRT_p", "nested_LRT_p", "LRT_p_boundary"), names(df))
  # Insert a significance-stars column immediately after each p-value column
  # (stars computed from the raw numeric p BEFORE it is formatted to "<.001").
  if (length(p_like)) {
    cols <- list(); nms <- character(0)
    for (nm in names(df)) {
      cols[[length(cols) + 1L]] <- df[[nm]]; nms <- c(nms, nm)
      if (nm %in% p_like) {
        cols[[length(cols) + 1L]] <- .stars(df[[nm]]); nms <- c(nms, paste0(nm, "__star"))
      }
    }
    df <- as.data.frame(cols, stringsAsFactors = FALSE); names(df) <- nms
  }
  onames <- names(df)
  # Numeric tidying: p-like columns -> "<.001"/3dp; other non-integer numeric
  # columns -> 4 significant digits (counts / df / years left intact).
  for (cc in p_like) df[[cc]] <- .fmt_p_vec(df[[cc]])
  keep_int <- c("N", "n", "n_iter", "BY_min", "BY_max", "LRT_df", "nested_LRT_df", "df",
    "N_west", "N_east", "N_test", "N_mobility", "N_overall", "N_East", "N_West",
    "N_Female", "N_Male", "Wf", "Wm", "Ef", "Em", "step_order")
  star_cols <- grep("__star$", onames, value = TRUE)
  for (cc in setdiff(onames, c(p_like, keep_int, star_cols))) {
    if (is.numeric(df[[cc]]) && !.is_int_col(df[[cc]])) df[[cc]] <- signif(df[[cc]], 4)
  }
  term_cols <- intersect(c("term", "focal_term", "smooth_term", "smooth", "smooth_label"), onames)
  for (cc in term_cols) df[[cc]] <- humanize_term(df[[cc]])
  val_cols <- setdiff(intersect(c("model", "analysis_id", "comparison", "scope", "family", "PGI",
    "estimator", "variance", "model_family", "cell", "group", "east_west", "gender", "fold",
    "held_out", "outcome", "region", "dataset", "section"), onames), term_cols)
  for (cc in val_cols) df[[cc]] <- .tidy_value(df[[cc]])
  if ("sample" %in% onames) df[["sample"]] <- .tidy_sample(df[["sample"]])
  names(df) <- vapply(onames, .human_col, character(1))
  df
}

humanize_sheets <- function(sheets) {
  out <- lapply(sheets, function(sh) {
    if (inherits(sh, "multiblock_sheet")) {
      structure(lapply(sh, function(b) { b$df <- humanize_df(b$df); b }),
                class = "multiblock_sheet", title = attr(sh, "title"))
    } else humanize_df(sh)
  })
  names(out) <- names(sheets)
  out
}

# --------------------------------------------------------------------------
# build_gam_overview_table — locked GAM overview schema, no birth-year
# anchors. Reports only what the section report needs to decide whether
# nonlinear modelling is required.
#
# Inputs:
#   nested_lrt   : output of nested_lrt_table(null_mod, lin_mod, nonlin_mod, label)
#                  — three rows (label_null, label_lin, label_nonlin) with
#                  AIC, BIC, logLik, dev_expl, LRT_chisq, LRT_df, LRT_p, comparison
#   smooth_table : output of tidy_smooths(nonlin_mod, label)
#                  — one row per smooth term with edf, ref_df, F, p
#   sample_name  : character(1)
#   N            : integer(1)
#   analysis_id  : "A3" / "B3" / "C3-edu" / "C3-mob"
#   alpha        : significance threshold (default .05)
# --------------------------------------------------------------------------

# Override each row's nonlinearity conclusion with the clean per-cell nested
# LRT (add only that cell's smooth vs the linear model), so the conclusion is
# not driven by the confounded summary.gam F. `cell_lrt` has columns
# cell / LRT_chisq / LRT_df / LRT_p / nonlinear_sig.
.apply_cell_lrt <- function(out, cell_lrt) {
  map <- c(WF = "West female", WM = "West male", EF = "East female", EM = "East male",
           W = "West", E = "East", F = "female", M = "male")
  for (i in seq_len(nrow(out))) {
    sl <- out$smooth_label[i]; if (is.na(sl)) next
    suf <- sub(".*PGI_", "", sl)
    cell <- if (suf %in% names(map)) map[[suf]] else NA_character_
    row <- if (!is.na(cell)) cell_lrt[as.character(cell_lrt$cell) == cell, , drop = FALSE] else cell_lrt[0, ]
    if (nrow(row)) {
      sig <- isTRUE(as.logical(row$nonlinear_sig[1]))
      out$nonlinearity_conclusion[i] <- if (sig) "Nonlinearity supported" else
        "No evidence that nonlinear specification is needed"
      out$interpretation[i] <- sprintf(
        "Clean per-cell nested LRT (adds only this cell's smooth): chi-sq=%.2f, df=%.0f, p=%.4g (%s).",
        as.numeric(row$LRT_chisq[1]), as.numeric(row$LRT_df[1]), as.numeric(row$LRT_p[1]),
        if (sig) "significant" else "n.s.")
    }
  }
  out
}

build_gam_overview_table <- function(nested_lrt, smooth_table, sample_name,
                                     N, analysis_id, alpha = .05, cell_lrt = NULL) {
  if (is.null(nested_lrt) || nrow(nested_lrt) < 3) {
    return(.gam_overview_row(
      analysis_id = analysis_id, sample_name = sample_name, N = N,
      alpha = alpha, smooth_label = NA_character_,
      smooth_row = NULL, lrt_chi2 = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_,
      aic_delta = NA_real_, bic_delta = NA_real_, fit_failed = TRUE))
  }

  # nested LRT row: linear vs nonlinear (third row of nested_lrt_table)
  lin_row     <- nested_lrt[2, , drop = FALSE]
  nonlin_row  <- nested_lrt[3, , drop = FALSE]
  lrt_chi2    <- nonlin_row$LRT_chisq
  lrt_df      <- nonlin_row$LRT_df
  lrt_p       <- nonlin_row$LRT_p
  aic_delta   <- nonlin_row$AIC - lin_row$AIC
  bic_delta   <- nonlin_row$BIC - lin_row$BIC

  # Emit one row per smooth term (so e.g. A3 produces two rows: West smooth, East smooth)
  if (is.null(smooth_table) || nrow(smooth_table) == 0) {
    return(.gam_overview_row(
      analysis_id = analysis_id, sample_name = sample_name, N = N,
      alpha = alpha, smooth_label = NA_character_,
      smooth_row = NULL, lrt_chi2 = lrt_chi2, lrt_df = lrt_df, lrt_p = lrt_p,
      aic_delta = aic_delta, bic_delta = bic_delta, fit_failed = FALSE))
  }

  out_rows <- list()
  for (i in seq_len(nrow(smooth_table))) {
    sr <- smooth_table[i, , drop = FALSE]
    out_rows[[i]] <- .gam_overview_row(
      analysis_id = analysis_id, sample_name = sample_name, N = N,
      alpha = alpha,
      smooth_label = sr$smooth_term,
      smooth_row = sr,
      lrt_chi2 = lrt_chi2, lrt_df = lrt_df, lrt_p = lrt_p,
      aic_delta = aic_delta, bic_delta = bic_delta,
      fit_failed = FALSE)
  }
  out <- do.call(rbind, out_rows)
  if (!is.null(cell_lrt) && !is.null(out) && nrow(out)) out <- .apply_cell_lrt(out, cell_lrt)
  out
}

# Internal: build a single row of the GAM overview table.
# The conclusion is deterministic — no model judgment, only test outcomes.
.gam_overview_row <- function(analysis_id, sample_name, N, alpha,
                              smooth_label, smooth_row,
                              lrt_chi2, lrt_df, lrt_p,
                              aic_delta, bic_delta, fit_failed = FALSE) {
  if (fit_failed) {
    conclusion <- "Not tested"
    interp     <- "Smooth fit unavailable; nonlinearity cannot be evaluated."
    edf <- NA_real_; ref_df <- NA_real_; F <- NA_real_; p_smooth <- NA_real_
  } else if (is.null(smooth_row)) {
    edf <- NA_real_; ref_df <- NA_real_; F <- NA_real_; p_smooth <- NA_real_
    if (is.na(lrt_p)) {
      conclusion <- "Not tested"
      interp     <- "LRT and smooth-test results both unavailable."
    } else if (lrt_p < alpha) {
      conclusion <- "Inconclusive"
      interp     <- sprintf("Nested LRT supports nonlinearity (chi2=%.2f, df=%.2f, p=%.4f) but no per-smooth test available.",
                            lrt_chi2, lrt_df, lrt_p)
    } else {
      conclusion <- "No evidence that nonlinear specification is needed"
      interp     <- sprintf("Nested LRT does not reject the linear specification (chi2=%.2f, df=%.2f, p=%.4f); per-smooth tests unavailable.",
                            lrt_chi2, lrt_df, lrt_p)
    }
  } else {
    edf      <- as.numeric(smooth_row$edf)
    ref_df   <- as.numeric(smooth_row$Ref.df)
    F        <- as.numeric(smooth_row$F)
    p_smooth <- as.numeric(smooth_row[["p-value"]])

    lrt_sig    <- !is.na(lrt_p)    && lrt_p    < alpha
    smooth_sig <- !is.na(p_smooth) && p_smooth < alpha
    lrt_null   <- !is.na(lrt_p)    && lrt_p    >= alpha
    edf_linear <- !is.na(edf)      && abs(edf - 1) < 0.1
    edf_boundary <- !is.na(edf)    && !is.na(ref_df) && edf > (ref_df - 0.1)

    if (lrt_sig && smooth_sig) {
      conclusion <- "Nonlinearity supported"
    } else if (lrt_null && edf_linear) {
      conclusion <- "No evidence that nonlinear specification is needed"
    } else if (edf_boundary) {
      conclusion <- "Inconclusive"
    } else if ((lrt_sig && !smooth_sig) || (!lrt_sig && smooth_sig)) {
      conclusion <- "Inconclusive"
    } else if (is.na(lrt_p) || is.na(p_smooth)) {
      conclusion <- "Not tested"
    } else {
      conclusion <- "Inconclusive"
    }

    interp <- sprintf(
      "Nested LRT chi2=%.2f, df=%.2f, p=%.4f. Smooth %s: edf=%.2f, F=%.2f, p=%.4f. See GAM curve in PDF.",
      ifelse(is.na(lrt_chi2), NA, lrt_chi2),
      ifelse(is.na(lrt_df), NA, lrt_df),
      ifelse(is.na(lrt_p), NA, lrt_p),
      smooth_label, edf, F, p_smooth
    )
  }

  data.frame(
    analysis_id             = analysis_id,
    sample                  = sample_name,
    N                       = N,
    nested_LRT_chi2         = lrt_chi2,
    nested_LRT_df           = lrt_df,
    nested_LRT_p            = lrt_p,
    nested_aic_delta        = aic_delta,
    nested_bic_delta        = bic_delta,
    smooth_label            = smooth_label,
    edf                     = edf,
    ref_df                  = ref_df,
    F                       = F,
    p                       = p_smooth,
    alpha                   = alpha,
    nonlinearity_conclusion = conclusion,
    interpretation          = interp,
    stringsAsFactors        = FALSE
  )
}

# --------------------------------------------------------------------------
# build_loo_table — format the output of leave_one_cohort_out() for the
# workbook. Takes the raw LOO data frame (fold / held_out / N / focal
# stats / status) and, if a numeric focal-estimate column is present,
# appends delta_vs_full = estimate(fold) - estimate(full sample).
#
#   loo_df       : data frame from leave_one_cohort_out()
#   estimate_col : name of the focal point-estimate column to difference
#                  against the full-sample row (NULL to skip, e.g. for the
#                  A3 nonlinearity-LRT table where the focal is a test,
#                  not a single estimate)
# --------------------------------------------------------------------------

build_loo_table <- function(loo_df, estimate_col = NULL) {
  if (is.null(loo_df) || !nrow(loo_df)) return(loo_df)
  if (!is.null(estimate_col) && estimate_col %in% names(loo_df)) {
    full_row <- loo_df[loo_df$fold == "Full sample", , drop = FALSE]
    full_est <- if (nrow(full_row)) suppressWarnings(as.numeric(full_row[[estimate_col]][1])) else NA_real_
    est_num  <- suppressWarnings(as.numeric(loo_df[[estimate_col]]))
    loo_df$delta_vs_full <- ifelse(loo_df$fold == "Full sample",
                                   NA_real_, round(est_num - full_est, 4))
  }
  loo_df
}

# --------------------------------------------------------------------------
# write_section_xlsx — emit a sheet-ordered workbook.
# `sheets` is a named list of data frames or multiblock_sheet objects.
# Set `include_index = TRUE` to prepend an auto-generated `00_Sheet_Index`
# sheet (default FALSE — with ≤ ~8 sheets per workbook the index is
# friction rather than navigation).
# --------------------------------------------------------------------------

write_section_xlsx <- function(section_id, sheets, path,
                               include_index = FALSE) {
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  ordered <- if (include_index) {
    c(list("00_Sheet_Index" = build_sheet_index(sheets)), sheets)
  } else {
    sheets
  }
  save_xlsx(ordered, path, label = sprintf("Section %s", section_id))
}

# --------------------------------------------------------------------------
# generate_section_overview_md — emit the GitHub-readable per-section
# overview markdown alongside the xlsx and PDF. Generated from the same
# data frames used for the workbook so the two artefacts cannot drift.
#
# Inputs:
#   section_id        : "A" / "B" / "C"
#   section_title     : human-readable title, e.g. "Section A — Educational attainment"
#   purpose           : single-paragraph string
#   analyses          : data.frame with columns (Analysis, Question, Model, Status)
#   findings          : data.frame with columns (Analysis, `Focal test`,
#                       `Estimate / statistic`, p, Conclusion)
#   interpretation    : single-paragraph string
#   outputs           : data.frame with columns (Output, File)
#   path              : output file path (.md)
# --------------------------------------------------------------------------

generate_section_overview_md <- function(section_id, section_title,
                                         purpose, analyses, findings,
                                         interpretation, outputs, path) {
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }

  con <- file(path, open = "w")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    sprintf("# %s", section_title),
    "",
    "## Purpose",
    "",
    purpose,
    "",
    "## Analyses",
    ""
  ), con)
  writeLines(.df_to_md_table(analyses), con)

  writeLines(c(
    "",
    "## Main findings",
    ""
  ), con)
  writeLines(.df_to_md_table(findings), con)

  writeLines(c(
    "",
    "## Interpretation",
    "",
    interpretation,
    "",
    "## Output locations",
    ""
  ), con)
  writeLines(.df_to_md_table(outputs), con)

  invisible(path)
}

# Internal: convert a data.frame into a GitHub-flavoured markdown table.
.df_to_md_table <- function(df) {
  if (is.null(df) || !nrow(df)) return("_(no rows)_")
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  align  <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows   <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(df[i, , drop = FALSE], function(x) {
      if (is.numeric(x) && !is.na(x)) format(x, trim = TRUE) else as.character(x)
    }, character(1))
    paste0("| ", paste(cells, collapse = " | "), " |")
  }, character(1))
  c(header, align, rows)
}
