# =============================================================================
# Shared phenotype helpers
# -----------------------------------------------------------------------------
# Small utility functions reused across cohort phenotype scripts. Source from
# any phenotype script via:
#
#   source(file.path("..", "utils", "edu_helpers.R"))
#
# (assumes the script is run from its own cohort directory, e.g.
# 02_Phenotype/BASE-II/).
# =============================================================================


#' Replace negative codes with NA
#'
#' Many German survey datasets encode missing-data reasons as negative
#' integers (e.g. -1 = "no answer", -2 = "does not apply"). This helper
#' coerces the input to numeric and replaces any negative value with NA.
#'
#' @param x  a vector (numeric, integer, character, or labelled)
#' @return   numeric vector of the same length with negatives -> NA
neg_to_na <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (!is.numeric(x_num)) return(x)
  x_num[x_num < 0] <- NA
  x_num
}


#' First non-missing value in time order
#'
#' Given a value vector and a time vector of the same length, return the
#' value at the earliest time where both value and time are non-missing.
#' Useful for collapsing person-wave panels to one row per person and
#' picking the first observed value of a measure.
#'
#' @param x     value vector
#' @param time  time vector (numeric or coercible)
#' @return      scalar; NA_real_ if no non-missing pair exists
earliest_nonmissing <- function(x, time) {
  t_num <- as.numeric(time)
  idx   <- which(!is.na(x) & !is.na(t_num))
  if (length(idx) == 0) return(NA_real_)
  x[idx[which.min(t_num[idx])]]
}


#' Last non-missing value in time order
#'
#' Mirror of `earliest_nonmissing()`: returns the value at the latest
#' time where both value and time are non-missing.
#'
#' @param x     value vector
#' @param time  time vector (numeric or coercible)
#' @return      scalar; NA_real_ if no non-missing pair exists
latest_nonmissing <- function(x, time) {
  t_num <- as.numeric(time)
  idx   <- which(!is.na(x) & !is.na(t_num))
  if (length(idx) == 0) return(NA_real_)
  x[idx[which.max(t_num[idx])]]
}


#' Mean ignoring NAs, returning NA when all-NA
#'
#' `mean(x, na.rm = TRUE)` returns NaN when every element is NA. This
#' helper returns NA_real_ in that case for consistency with the rest of
#' the dataset.
#'
#' @param x  numeric vector
#' @return   scalar mean, or NA_real_ if all-NA
mean_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}
