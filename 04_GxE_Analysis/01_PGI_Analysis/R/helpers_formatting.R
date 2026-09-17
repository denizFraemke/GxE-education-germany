# helpers_formatting.R
# Numeric / standardization utilities.

z_scale <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

z_within_group <- function(x, g, min_n = 5) {
  ave(x, g, FUN = function(v) {
    v <- as.numeric(v); n_ok <- sum(!is.na(v))
    if (n_ok < min_n) return(rep(NA_real_, length(v)))
    s <- sd(v, na.rm = TRUE); m <- mean(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
    (v - m) / s
  })
}

cohort_bin <- function(birth_year, width = COHORT_BIN_WIDTH) {
  floor(birth_year / width) * width
}
