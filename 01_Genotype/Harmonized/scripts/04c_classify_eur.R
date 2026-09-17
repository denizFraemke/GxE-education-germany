#!/usr/bin/env Rscript
# =============================================================================
# 04c_classify_eur.R — classify European ancestry from 1000G-projected PCs
# =============================================================================
# Reads the projected PC tables written by 04c_ancestry_1kg.sh (reference
# samples with super_pop labels + each cohort), fits a transparent per-
# population Gaussian rule in top-K PC space, and writes one EUR keep-list per
# cohort in the format filter_to_eur() expects (a single "IID" column).
#
# Rule: for each individual, compute the Mahalanobis distance to every 1000G
# super-population centroid (each with that population's own covariance) in the
# first K projected PCs. Keep the individual iff the nearest centroid is EUR AND
# its distance to the EUR centroid is within the chi-square cutoff (df = K).
# Both K and the cutoff are parameters and are echoed to the log.
#
# SOEP is not classified here — its keep-list is computed elsewhere.
#
# USAGE
#   Rscript scripts/04c_classify_eur.R                       # uses defaults
#   ANC_DIR=... OUT_DIR=... EXISTING_KEEP_DIR=... N_PC=6 EUR_Q=0.999 \
#     Rscript scripts/04c_classify_eur.R
# =============================================================================
suppressWarnings(suppressPackageStartupMessages({}))

ANC_DIR  <- Sys.getenv("ANC_DIR",  file.path(Sys.getenv("OUT_ROOT", "output"), "ancestry"))
OUT_DIR  <- Sys.getenv("OUT_DIR",  file.path(ANC_DIR, "keep_lists"))
EXISTING <- Sys.getenv("EXISTING_KEEP_DIR", "")   # optional: for concordance check
N_PC     <- as.integer(Sys.getenv("N_PC", "6"))   # PCs used for classification
EUR_Q    <- as.numeric(Sys.getenv("EUR_Q", "0.999"))  # chi-square keep cutoff
COHORT_KEYS <- c("BASEII", "SHIP0", "SHIPTd", "TwinLife")
SUPERPOPS   <- c("EUR", "AFR", "EAS", "SAS", "AMR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- helpers -----------------------------------------------------------------
read_pcs <- function(path) {
  d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
  names(d) <- sub("_SUM$", "", names(d))            # PLINK --score emits PC<n>_SUM -> PC<n>
  names(d)[grepl("IID", names(d))][1] <- "IID"
  d
}
pc_cols <- function(d) paste0("PC", seq_len(N_PC))

# --- reference: centroids + covariances per super-population ------------------
ref <- read_pcs(file.path(ANC_DIR, "ref_projected_pcs.tsv"))
stopifnot("super_pop" %in% names(ref))
pcs <- pc_cols(ref)
if (!all(pcs %in% names(ref))) stop("Reference lacks PC1..PC", N_PC)

cutoff <- qchisq(EUR_Q, df = N_PC)
centroid <- list(); covmat <- list()
for (sp in SUPERPOPS) {
  M <- as.matrix(ref[ref$super_pop == sp, pcs, drop = FALSE])
  if (nrow(M) < N_PC + 2) stop("Too few reference samples for ", sp)
  centroid[[sp]] <- colMeans(M)
  covmat[[sp]]   <- cov(M)
}

# Non-admixed reference populations used to decide the "nearest" label. AMR
# (1000G Admixed American) is deliberately excluded: it is an admixed population
# whose broad covariance overlaps the European cluster, so a plain
# nearest-centroid rule spuriously assigns clearly-European individuals to AMR.
# European ancestry is therefore: within the EUR Mahalanobis cutoff AND with EUR
# the nearest of the non-admixed references (EUR/AFR/EAS/SAS).
NONADMIXED <- c("EUR", "AFR", "EAS", "SAS")
classify_eur <- function(d) {
  X <- as.matrix(d[, pcs, drop = FALSE])
  D <- sapply(SUPERPOPS, function(sp) mahalanobis(X, centroid[[sp]], covmat[[sp]]))
  nearest <- NONADMIXED[max.col(-D[, NONADMIXED, drop = FALSE], ties.method = "first")]
  keep <- nearest == "EUR" & D[, "EUR"] <= cutoff
  list(keep = keep, nearest = nearest, d_eur = D[, "EUR"])
}

cat(sprintf("Classification: K=%d PCs, EUR keep cutoff = chisq_%.3f(df=%d) = %.2f\n\n",
            N_PC, EUR_Q, N_PC, cutoff))

# --- classify each cohort, write keep-lists, and report ----------------------
summary_rows <- list()
for (key in COHORT_KEYS) {
  f <- file.path(ANC_DIR, paste0(key, "_projected_pcs.tsv"))
  if (!file.exists(f)) { cat(sprintf("[%s] projected PCs not found — skipped\n", key)); next }
  d <- read_pcs(f)
  cl <- classify_eur(d)
  keep_iid <- d$IID[cl$keep]

  out <- file.path(OUT_DIR, paste0("EUR_keep_IDs_", key, ".tsv"))
  write.table(data.frame(IID = keep_iid), out, sep = "\t",
              quote = FALSE, row.names = FALSE)

  line <- sprintf("[%s] projected=%d  EUR-kept=%d  dropped=%d",
                  key, nrow(d), sum(cl$keep), sum(!cl$keep))

  # optional concordance vs an existing keep-list
  if (nzchar(EXISTING)) {
    ef <- file.path(EXISTING, paste0("EUR_keep_IDs_", key, ".tsv"))
    if (file.exists(ef)) {
      old <- as.character(read.table(ef, header = TRUE, sep = "\t")$IID)
      new <- as.character(keep_iid)
      both <- length(intersect(new, old))
      line <- paste0(line, sprintf("  | vs existing: kept-both=%d  new-only=%d  dropped-vs-old=%d  agree=%.1f%%",
                     both, length(setdiff(new, old)), length(setdiff(old, new)),
                     100 * both / length(union(new, old))))
    }
  }
  cat(line, "\n")
  summary_rows[[key]] <- line
}

cat("\nWrote keep-lists to: ", OUT_DIR, "\n")
cat("SOEP is not classified here — its keep-list, EUR_keep_IDs_SOEP.tsv, comes from elsewhere.\n")
