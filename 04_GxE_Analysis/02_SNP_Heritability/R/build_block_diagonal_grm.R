#!/usr/bin/env Rscript
########################################################################
## build_block_diagonal_grm.R — combine disjoint-sample GRMs as a
## block-diagonal union.
##
## GCTA's `--mgrm --make-grm-bin` is built for combining GRMs **of the
## same individuals** computed on different SNP sets (e.g. per-chr).
## It computes the intersection of IDs across input GRMs and sums their
## values. For disjoint individual sets (our stratum-by-stratum design)
## that intersection is empty and gcta segfaults.
##
## What this script produces instead is a **block-diagonal union**:
## given K disjoint sub-GRMs G_k (each n_k × n_k), produce an N × N
## matrix (N = sum n_k) with G_k on each diagonal block and zeros for
## cross-block pairs. The result is a valid GCTA GRM that can be passed
## to `gcta --reml --grm-bin <prefix>` as a single variance-component
## GRM in which V_G is identified from within-block relatedness only.
##
## Usage (called by scripts/04d_build_block_diagonal_cell_grms.sh):
##   Rscript build_block_diagonal_grm.R <out_prefix> <in_prefix_1> <in_prefix_2> [...]
##
## Output:
##   <out_prefix>.grm.bin   — packed lower-triangular float32, block-diagonal
##   <out_prefix>.grm.id    — concatenated FID/IID rows (stratum stacking order)
##   <out_prefix>.grm.N.bin — packed lower-triangular float32; SNP counts
##                            within-block, 0 for cross-block pairs
##
## Notes:
##   - GCTA's .grm.bin / .grm.N.bin are headerless: packed float32 values
##     in row-major order over the lower triangle (including the diagonal).
##     Row i has i entries. Total length = N*(N+1)/2.
##   - Memory: combined GRM at N=14k has N*(N+1)/2 ~ 98M entries × 8 bytes
##     in R doubles ≈ 750 MB temporarily; written out as float32 (~390 MB).
########################################################################

suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage: Rscript build_block_diagonal_grm.R <out_prefix> <in_prefix_1> [<in_prefix_2> ...]",
    "\nAt least one input GRM is required. With a single input the output is",
    "byte-identical to the input (useful for 1-block validation of the",
    "combiner — same I/O path, same indexing logic, no cross-block zeros to add)."
  ))
}
out_prefix  <- args[1]
in_prefixes <- args[-1]

# ---------- Read one GRM ---------------------------------------------------
read_grm <- function(prefix) {
  id_file  <- paste0(prefix, ".grm.id")
  bin_file <- paste0(prefix, ".grm.bin")
  n_file   <- paste0(prefix, ".grm.N.bin")
  for (f in c(id_file, bin_file, n_file)) {
    if (!file.exists(f)) stop("Missing file: ", f)
  }

  ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("FID", "IID"))
  N <- nrow(ids)
  n_lt <- as.numeric(N) * (as.numeric(N) + 1) / 2  # numeric to avoid int overflow

  con <- file(bin_file, "rb")
  vals <- readBin(con, what = "numeric", n = n_lt, size = 4L, endian = "little")
  close(con)
  if (length(vals) != n_lt)
    stop(sprintf("Bin length mismatch for %s: got %d, expected %d",
                 bin_file, length(vals), n_lt))

  con <- file(n_file, "rb")
  Nsnps <- readBin(con, what = "numeric", n = n_lt, size = 4L, endian = "little")
  close(con)

  list(ids = ids, vals = vals, Nsnps = Nsnps, N = N)
}

# ---------- Read all ---------------------------------------------------------
cat(sprintf("Reading %d input GRM(s)...\n", length(in_prefixes)))
grms <- lapply(in_prefixes, function(p) {
  cat(sprintf("  %s\n", p))
  read_grm(p)
})

total_N <- sum(sapply(grms, function(g) g$N))
cat(sprintf("Combined N = %d (sum of %s)\n",
            total_N,
            paste(sapply(grms, function(g) g$N), collapse = " + ")))

# ---------- Allocate combined buffers ----------------------------------------
n_total_lt <- as.numeric(total_N) * (as.numeric(total_N) + 1) / 2
cat(sprintf("Allocating combined lower-triangular: %s entries\n",
            format(n_total_lt, big.mark = ",")))
combined_vals  <- numeric(n_total_lt)   # zero-initialised; cross-block stays 0
combined_Nsnps <- numeric(n_total_lt)

# ---------- Fill in each stratum's block --------------------------------------
# For combined row r = row_offset + i (i in 1..N_k):
#   - Row r in the lower triangle has r entries (indices r*(r-1)/2+1 .. r*(r+1)/2).
#   - First row_offset entries (cross-stratum pairs): zero — already initialised.
#   - Last i entries (within-stratum row i of g): copied from g.
row_offset <- 0
for (g_idx in seq_along(grms)) {
  g <- grms[[g_idx]]
  cat(sprintf("Placing stratum %d (N=%d) at rows %d..%d\n",
              g_idx, g$N, row_offset + 1, row_offset + g$N))
  for (i in seq_len(g$N)) {
    combined_row <- row_offset + i
    # 1-based start index of row combined_row in the lower triangle:
    combined_start <- as.numeric(combined_row - 1) * as.numeric(combined_row) / 2 + 1
    # 1-based start of row i in stratum g's triangle:
    g_start <- as.numeric(i - 1) * as.numeric(i) / 2 + 1

    # Destination slice (last i entries of combined row; first row_offset are zero):
    dst_lo <- combined_start + row_offset
    dst_hi <- combined_start + combined_row - 1
    src_lo <- g_start
    src_hi <- g_start + i - 1

    combined_vals[dst_lo:dst_hi]  <- g$vals[src_lo:src_hi]
    combined_Nsnps[dst_lo:dst_hi] <- g$Nsnps[src_lo:src_hi]
  }
  row_offset <- row_offset + g$N
}

# ---------- Write outputs -----------------------------------------------------
combined_ids <- do.call(rbind, lapply(grms, function(g) g$ids))
cat(sprintf("Writing combined GRM: %s.grm.{bin,id,N.bin}\n", out_prefix))

write.table(combined_ids, paste0(out_prefix, ".grm.id"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

con <- file(paste0(out_prefix, ".grm.bin"), "wb")
writeBin(combined_vals, con, size = 4L, endian = "little")
close(con)

con <- file(paste0(out_prefix, ".grm.N.bin"), "wb")
writeBin(combined_Nsnps, con, size = 4L, endian = "little")
close(con)

cat(sprintf("Done. Combined N = %d, components = %d.\n",
            total_N, length(in_prefixes)))
