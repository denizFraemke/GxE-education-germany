#!/usr/bin/env Rscript
########################################################################
## 08_power_calc.R — Per-stratum power via Visscher 2014 (Tardis)
##
## For each of the 4 main strata:
##   1. read the unrelated GRM produced by step 04 and compute
##      Var(off-diag) empirically (one streaming pass — no full N x N
##      matrix in memory);
##   2. compute SE(h^2) = sqrt(2 / (N^2 * Var(off-diag))) per Visscher
##      et al. (2014);
##   3. evaluate one-sided Wald and chi^2 LRT power across the configured
##      h^2 grid.
##
## Output: ${REML_OUT}/power.tsv
########################################################################

rm(list = ls())

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

GRM_OUT   <- file.path(PROJ_ROOT, "output", "grm")
REML_OUT  <- file.path(PROJ_ROOT, "output", "reml")
OUT_FILE  <- file.path(REML_OUT, "power.tsv")
dir.create(REML_OUT, recursive = TRUE, showWarnings = FALSE)

# Read the same stratum + grid config the shell pipeline uses.
# Pure-R regex parser: handles `export NAME=(...)` and `NAME=(...)`,
# single-line and multi-line array bodies, and trailing `# comments`.
# Parsing the file directly avoids shelling out to bash just to read
# two array literals.
read_config_array <- function(name) {
  cfg <- readLines(file.path(PROJ_ROOT, "config.sh"), warn = FALSE)
  pat_start <- sprintf("^[[:space:]]*(export[[:space:]]+)?%s=\\(", name)
  i <- grep(pat_start, cfg)
  if (length(i) == 0L) return(character(0))
  start_line <- i[1L]

  # Find the line containing the closing ')' (may be the same line).
  end_line <- start_line
  if (!grepl("\\)", cfg[start_line])) {
    after <- which(grepl("\\)", cfg) & seq_along(cfg) > start_line)
    if (length(after) == 0L) return(character(0))
    end_line <- after[1L]
  }

  # Concatenate the relevant lines, then pull out the text between
  # `<NAME>=(` and the first `)`.
  raw <- paste(cfg[start_line:end_line], collapse = "\n")
  open_re <- sprintf("(export[[:space:]]+)?%s=\\(", name)
  open_pos <- regexpr(open_re, raw)
  if (open_pos == -1) return(character(0))
  after_open <- substring(raw, open_pos + attr(open_pos, "match.length"))
  close_pos <- regexpr("\\)", after_open)
  if (close_pos == -1) return(character(0))
  inner <- substring(after_open, 1, close_pos - 1)

  # Strip trailing comments on any line (config.sh sometimes has `# foo`
  # after a value), then split on whitespace and strip surrounding quotes.
  inner <- gsub("#[^\n]*", "", inner)
  out <- strsplit(trimws(inner), "[[:space:]]+", perl = TRUE)[[1L]]
  out <- gsub('^"|"$', "", out)
  out <- gsub("^'|'$", "", out)
  out[nzchar(out)]
}

STRATA_MAIN <- read_config_array("STRATA_MAIN")
H2_GRID     <- as.numeric(read_config_array("H2_GRID"))

if (length(STRATA_MAIN) == 0L || length(H2_GRID) == 0L)
  stop("Could not parse STRATA_MAIN / H2_GRID from config.sh.")

cat("=== Step 08: Power calc (Visscher 2014) ===\n")
cat(sprintf("Strata:    %s\n", paste(STRATA_MAIN, collapse = ", ")))
cat(sprintf("h^2 grid:  %s\n", paste(H2_GRID, collapse = ", ")))

build_rows <- function(strata, analysis_label) {
  rows <- list()
  for (s in strata) {
    grm_prefix <- file.path(GRM_OUT, paste0("grm_", s, "_unrel"))
    if (!file.exists(paste0(grm_prefix, ".grm.bin"))) {
      cat(sprintf("  %s: GRM missing (%s) — skipping\n", s, grm_prefix))
      next
    }
    cat(sprintf("  %s: reading GRM off-diagonal ...\n", s))
    stat <- read_grm_offdiag_var(grm_prefix)
    N <- stat$N
    var_off <- stat$var
    se_h2 <- visscher_se_h2(N, var_off)
    for (h2 in H2_GRID) {
      pw <- visscher_power(h2, se_h2, alpha = 0.05)
      rows[[length(rows) + 1L]] <- data.frame(
        analysis        = analysis_label,
        stratum         = s,
        N               = N,
        var_offdiag_grm = var_off,
        h2_scenario     = h2,
        SE_h2_visscher  = se_h2,
        z               = if (is.na(se_h2)) NA_real_ else h2 / se_h2,
        power_wald_a05  = pw$wald,
        power_lrt_a05   = pw$lrt,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

out <- build_rows(STRATA_MAIN, "main_4cell")

if (is.null(out) || nrow(out) == 0L)
  stop("No power rows produced — were any unrelated GRMs available?")

out$timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
write.table(out, OUT_FILE,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWrote: %s (%d rows)\n", OUT_FILE, nrow(out)))

# Compact summary of SE at the centre of the grid.
mid <- 0.20
summary_rows <- out[out$h2_scenario == mid, c("analysis", "stratum", "N",
                                              "SE_h2_visscher",
                                              "power_lrt_a05")]
cat("\nAt h^2 = 0.20:\n")
print(summary_rows)
cat("\n=== Step 08 complete ===\n")
