# helpers_loo.R
# Leave-one-study-out (LOO) sensitivity infrastructure. The contributing
# study is carried in the `cohort` column, so the function and argument
# names below use that column name.
#
# Each focal model is refit on the analytic sample minus one study fold;
# the focal result is tracked across folds to show how robust it is to
# any single study. SHIP is held out as one fold: if the analytic sample
# carries separate SHIP-0 / SHIP-Td levels they are collapsed into a single
# "SHIP" fold here, per the preregistered plan, which treats SHIP as one
# study.

# cohort_groups — define the LOO folds from the studies present in `dat`.
# Returns a named list mapping fold label -> the study level(s) that
# fold holds out. SHIP-0 / SHIP-Td (if present) collapse into "SHIP".
cohort_groups <- function(dat, cohort_var = "cohort") {
  levs <- sort(unique(as.character(dat[[cohort_var]])))
  levs <- levs[!is.na(levs)]
  ship_levels <- levs[grepl("^SHIP", levs)]
  non_ship    <- levs[!grepl("^SHIP", levs)]
  folds <- stats::setNames(as.list(non_ship), non_ship)
  if (length(ship_levels) > 0) folds[["SHIP"]] <- ship_levels
  # Order: alphabetical non-SHIP, then SHIP last for readability.
  folds[order(names(folds) == "SHIP", names(folds))]
}

# leave_one_cohort_out — generic LOO engine.
#
#   dat             : full analytic sample (must contain `cohort_var`)
#   fit_and_extract : function(subsample) -> named numeric/character
#                     vector of focal statistics for that subsample
#                     (e.g. c(estimate=..., SE=..., p=...) for a linear
#                     focal, or c(LRT_chisq=..., LRT_df=..., LRT_p=...,
#                     nonlinearity_conclusion=...) for a GAM focal).
#   folds           : named list from cohort_groups(dat)
#   cohort_var      : name of the study column
#
# Returns a data frame with one "Full sample" row followed by one row per
# fold, columns: fold, held_out, N, <the names returned by
# fit_and_extract...>, status. A fold whose study is absent from `dat`
# (e.g. SHIP in the mobility sample) is reported with status
# "N/A - study absent" and no refit; a fit that errors is reported with
# status "fit_failed".
leave_one_cohort_out <- function(dat, fit_and_extract,
                                 folds = cohort_groups(dat),
                                 cohort_var = "cohort") {
  run_one <- function(sub, fold_label, held_out, status) {
    if (!identical(status, "OK")) {
      return(data.frame(fold = fold_label, held_out = held_out,
                        N = nrow(sub), status = status,
                        stringsAsFactors = FALSE))
    }
    stats <- tryCatch(fit_and_extract(sub), error = function(e) {
      attr(e, "loo_failed") <- TRUE; e
    })
    if (inherits(stats, "error")) {
      return(data.frame(fold = fold_label, held_out = held_out,
                        N = nrow(sub), status = "fit_failed",
                        stringsAsFactors = FALSE))
    }
    cbind(
      data.frame(fold = fold_label, held_out = held_out, N = nrow(sub),
                 stringsAsFactors = FALSE),
      as.data.frame(as.list(stats), stringsAsFactors = FALSE),
      data.frame(status = "OK", stringsAsFactors = FALSE)
    )
  }

  rows <- list()
  # Full-sample baseline.
  rows[["Full sample"]] <- run_one(dat, "Full sample", "(none)", "OK")
  # One row per fold.
  for (fl in names(folds)) {
    drop_levels <- folds[[fl]]
    present <- any(as.character(dat[[cohort_var]]) %in% drop_levels)
    if (!present) {
      rows[[fl]] <- run_one(dat[0, , drop = FALSE],
                            paste0("-", fl), fl, "N/A - study absent")
      next
    }
    sub <- dat[!(as.character(dat[[cohort_var]]) %in% drop_levels), , drop = FALSE]
    # Preserve the by_mean attribute used by extraction helpers.
    attr(sub, "by_mean") <- attr(dat, "by_mean")
    rows[[fl]] <- run_one(sub, paste0("-", fl), fl, "OK")
  }
  dplyr::bind_rows(rows)
}
