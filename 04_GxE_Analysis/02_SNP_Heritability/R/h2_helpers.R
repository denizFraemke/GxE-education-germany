########################################################################
## h2_helpers.R — Shared helpers for the SNP-heritability pipeline (Tardis)
##
## Used by scripts/07_lrt.R, scripts/08_power_calc.R, scripts/09_assemble_output.R.
##
## Functions provided:
##   parse_hsq(path)                   — read a single GCTA .hsq file into a
##                                       named list of estimates.
##   read_grm_offdiag_var(grm_prefix)  — empirical variance of the off-diagonal
##                                       of a binary GCTA GRM (Visscher 2014
##                                       power-calc input).
##   visscher_se_h2(N, var_offdiag)    — analytical SE(h²) per Visscher et al.
##                                       (2014), Eq. 2.
##   visscher_power(h2, se_h2)         — one-sided Wald + chi² LRT power at
##                                       α = .05 (default).
##   lrt_pvalue(logL_full, logL_nested, df) — chi² LRT p.
########################################################################

suppressPackageStartupMessages({
  # base R only — keep helpers dependency-light so they run on any Tardis R.
})

# ---------- .hsq parser -----------------------------------------------------
# A GCTA --reml .hsq file is a simple key/value table; the multi-GRM variant
# has rows V(G1), V(G2), ..., V(e), Vp, V(G)/Vp, logL, logL0, LRT, df, n_snps,
# n. We parse it to a named list keyed by the raw labels (with whitespace
# stripped). Rows are "name value SE?".
parse_hsq <- function(path) {
  if (!file.exists(path)) {
    return(list(.path = path, .present = FALSE))
  }
  raw <- readLines(path, warn = FALSE)
  raw <- raw[nzchar(trimws(raw))]
  if (length(raw) < 2L)
    return(list(.path = path, .present = FALSE))

  # First line is a header ("Source\tVariance\tSE"); skip it. Some GCTA
  # versions write tab-separated, others whitespace; split flexibly.
  rows <- strsplit(raw[-1], "\\s+", perl = TRUE)
  out  <- list(.path = path, .present = TRUE)
  for (r in rows) {
    if (length(r) < 2L) next
    key <- r[[1L]]
    val <- suppressWarnings(as.numeric(r[[2L]]))
    out[[key]] <- val
    if (length(r) >= 3L) {
      se <- suppressWarnings(as.numeric(r[[3L]]))
      out[[paste0(key, "_SE")]] <- se
    }
  }
  out
}

# ---------- GRM-binary reader -----------------------------------------------
# GCTA's GRM-bin format: <prefix>.grm.bin is a packed lower-triangular matrix
# (including the diagonal) of single-precision floats; <prefix>.grm.N.bin is
# the matching count; <prefix>.grm.id is "FID IID" per individual.
#
# We compute Var(off-diag) on the fly to avoid materialising the full N x N
# matrix (an unrelated cohort GRM has ~N^2 / 2 off-diagonal entries; with
# N = 10k that's 50M floats — fits in RAM but no reason to load it all when
# all we need is mean + variance).
read_grm_offdiag_var <- function(grm_prefix) {
  id_file  <- paste0(grm_prefix, ".grm.id")
  bin_file <- paste0(grm_prefix, ".grm.bin")
  if (!file.exists(id_file) || !file.exists(bin_file))
    stop("GRM files missing for prefix: ", grm_prefix)

  ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("FID", "IID"))
  n <- nrow(ids)
  total_entries <- as.numeric(n) * (as.numeric(n) + 1) / 2

  con <- file(bin_file, "rb")
  on.exit(close(con))

  # Streaming pass: read one row of the packed triangle at a time. Row i
  # has i entries; the last entry on each row is the diagonal — skip it.
  s   <- 0
  ss  <- 0
  cnt <- 0
  for (i in seq_len(n)) {
    v <- readBin(con, what = "numeric", n = i, size = 4L, endian = "little")
    if (i >= 2L) {
      off <- v[-i]                # drop the diagonal entry
      s   <- s  + sum(off)
      ss  <- ss + sum(off * off)
      cnt <- cnt + length(off)
    }
  }

  if (cnt < 2L)
    return(list(N = n, n_offdiag = cnt, mean = NA_real_, var = NA_real_))

  mean_off <- s / cnt
  var_off  <- (ss - cnt * mean_off^2) / (cnt - 1L)
  list(N = n, n_offdiag = cnt, mean = mean_off, var = var_off)
}

# ---------- Visscher 2014 SE(h²) + power ------------------------------------
# Visscher et al. (2014), "Statistical Power to Detect Genetic (Co)Variance of
# Complex Traits Using SNP Data in Unrelated Samples", PLOS Genet 10(4):
# SE(h²_SNP) ≈ sqrt(2 / (N^2 · Var(off-diag GRM))) for a single-component
# REML in an unrelated sample. Independent of the true h² to leading order.
visscher_se_h2 <- function(N, var_offdiag) {
  if (is.na(N) || is.na(var_offdiag) || var_offdiag <= 0)
    return(NA_real_)
  sqrt(2 / (as.numeric(N)^2 * var_offdiag))
}

# Wald and chi² LRT power at alpha (default .05, two-sided Wald, df=1 LRT).
visscher_power <- function(h2, se_h2, alpha = 0.05) {
  if (is.na(h2) || is.na(se_h2) || se_h2 <= 0)
    return(list(wald = NA_real_, lrt = NA_real_))
  z_crit <- qnorm(1 - alpha / 2)
  ncp_z  <- h2 / se_h2
  wald   <- pnorm(z_crit - ncp_z, lower.tail = FALSE) +
              pnorm(-z_crit - ncp_z, lower.tail = TRUE)
  chi2_crit <- qchisq(1 - alpha, df = 1)
  lrt <- pchisq(chi2_crit, df = 1, ncp = ncp_z^2, lower.tail = FALSE)
  list(wald = wald, lrt = lrt)
}

# ---------- LRT p-value ------------------------------------------------------
lrt_pvalue <- function(logL_full, logL_nested, df) {
  if (any(is.na(c(logL_full, logL_nested, df))))
    return(list(LRT = NA_real_, p = NA_real_))
  stat <- -2 * (logL_nested - logL_full)
  if (stat < 0) stat <- 0   # tiny numerical negatives from REML round-trip
  list(LRT = stat, p = pchisq(stat, df = df, lower.tail = FALSE))
}
