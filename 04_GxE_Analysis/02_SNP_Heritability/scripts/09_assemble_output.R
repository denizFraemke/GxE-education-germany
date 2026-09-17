#!/usr/bin/env Rscript
########################################################################
## 09_assemble_output.R — Final outputs + dashboard row (Tardis)
##
## Reads everything step 05/07/08 produced, builds three TSVs in
## ${FINAL_OUT}/ ready for download to the Mac, and emits a one-line
## markdown dashboard row appendable to the section README.
##
## Outputs:
##   output/final/h2_snp.tsv      — one row per (analysis, model, stratum)
##   output/final/lrt.tsv         — copy of step 07 result for symmetry
##   output/final/power.tsv       — copy of step 08 result for symmetry
##   output/final/manifest.tsv    — provenance
##   output/final/dashboard_row.md — markdown row to paste into README
########################################################################

rm(list = ls())

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
PROJ_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

source(file.path(PROJ_ROOT, "R", "h2_helpers.R"))

REML_OUT  <- file.path(PROJ_ROOT, "output", "reml")
MAIN_DIR  <- file.path(REML_OUT, "main")
FINAL_OUT <- file.path(PROJ_ROOT, "output", "final")
dir.create(FINAL_OUT, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 09: Assemble final outputs ===\n")

# ---- helper: build one h2_snp row from a parsed .hsq ----------------
extract_h2_rows <- function(hsq, analysis, model, stratum_labels) {
  if (!isTRUE(hsq$.present)) return(NULL)

  # Per-cell genetic variance components V(G1), V(G2), ... — only present
  # in joint mgrm fits. Standalone single-GRM .hsq has V(G).
  v_g_names <- sort(grep("^V\\(G\\d*\\)$", names(hsq), value = TRUE))
  if (length(v_g_names) == 0L) return(NULL)

  v_e   <- hsq[["V(e)"]]
  v_e_se <- hsq[["V(e)_SE"]]
  v_p   <- hsq[["Vp"]]

  rows <- list()
  for (k in seq_along(v_g_names)) {
    g_name <- v_g_names[k]
    g_se_name <- paste0(g_name, "_SE")
    v_g   <- hsq[[g_name]]
    v_g_se <- hsq[[g_se_name]]

    stratum <- if (length(v_g_names) == 1L && length(stratum_labels) >= 1L) {
      stratum_labels[1L]
    } else if (length(stratum_labels) >= k) {
      stratum_labels[k]
    } else {
      paste0("cell_", k)
    }

    # h^2 for the k-th component. For mgrm fits gcta also writes
    # V(G_k)/Vp rows; prefer those when present (they include the
    # joint-fit SE rather than a naive ratio).
    ratio_name <- sprintf("V(G%d)/Vp", k)
    if (length(v_g_names) == 1L && "V(G)/Vp" %in% names(hsq)) {
      h2 <- hsq[["V(G)/Vp"]]
      h2_se <- hsq[["V(G)/Vp_SE"]]
    } else if (ratio_name %in% names(hsq)) {
      h2 <- hsq[[ratio_name]]
      h2_se <- hsq[[paste0(ratio_name, "_SE")]]
    } else {
      h2 <- if (!is.na(v_p) && v_p > 0) v_g / v_p else NA_real_
      h2_se <- NA_real_
    }

    # Safe-getter — GCTA standalone .hsq files don't include n_SNPs (only
    # the joint mGRM fits do); reading NULL into data.frame() errors with
    # "differing number of rows". Default to NA when a field is absent.
    sg <- function(key, default = NA_real_) {
      v <- hsq[[key]]
      if (is.null(v) || length(v) == 0L) default else v
    }

    rows[[k]] <- data.frame(
      analysis = analysis,
      model    = model,
      stratum  = stratum,
      N        = sg("n", NA_integer_),
      n_snps   = sg("n_SNPs", NA_integer_),
      V_G      = v_g,
      V_G_SE   = if (is.null(v_g_se)) NA_real_ else v_g_se,
      V_e      = sg("V(e)"),
      V_e_SE   = sg("V(e)_SE"),
      V_P      = v_p,
      h2_SNP   = h2,
      h2_SNP_SE = if (is.null(h2_se)) NA_real_ else h2_se,
      logL     = sg("logL"),
      converged = TRUE,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

ORDER_MAIN <- c("East_pre1990", "East_post1990", "West_pre1990", "West_post1990")

h2_rows <- list()

# Per-stratum standalones (4-cell main analysis).
for (s in ORDER_MAIN) {
  for (suffix in c("", "_noconstrain")) {
    f <- file.path(MAIN_DIR, paste0("standalone_", s, suffix, ".hsq"))
    h <- parse_hsq(f)
    model <- if (suffix == "_noconstrain") "per_stratum_standalone_noconstrain"
             else "per_stratum_standalone"
    r <- extract_h2_rows(h, "main_4cell", model, s)
    if (!is.null(r)) h2_rows[[length(h2_rows) + 1L]] <- r
  }
}

h2_df <- if (length(h2_rows) > 0L) do.call(rbind, h2_rows) else
         data.frame()
h2_df$timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")

write.table(h2_df,
            file.path(FINAL_OUT, "h2_snp.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote: %s (%d rows)\n",
            file.path(FINAL_OUT, "h2_snp.tsv"), nrow(h2_df)))

# Copy through the heterogeneity TSV + power for a one-stop deliverable.
# Both filenames travel downstream: heterogeneity.tsv is the canonical
# name, lrt.tsv is the compatibility symlink that step 07 creates.
for (f in c("heterogeneity.tsv", "lrt.tsv", "power.tsv")) {
  src <- file.path(REML_OUT, f)
  dst <- file.path(FINAL_OUT, f)
  if (file.exists(src) || (!is.na(Sys.readlink(src)) && nzchar(Sys.readlink(src)))) {
    file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE)
    cat(sprintf("Copied: %s -> %s\n", src, dst))
  } else if (f != "lrt.tsv") {
    cat(sprintf("WARN: missing %s — step 07/08 didn't run?\n", src))
  }
}

# ---- Manifest -------------------------------------------------------
upstream_manifest <- file.path(PROJ_ROOT, "data", "manifest.tsv")
upstream <- if (file.exists(upstream_manifest)) {
  read.table(upstream_manifest, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE)
} else {
  data.frame()
}

gcta_version <- tryCatch({
  # GCTA's banner is a multi-line ASCII frame; the actual version string
  # lives on the line containing "version" (e.g. "* version v1.95.1
  # Linux"), not on the first line, which is the asterisk border.
  banner <- system(sprintf("'%s' 2>&1 | head -5",
                           Sys.getenv("GCTA", "/path/to/cluster_bin_dir/gcta64")),
                   intern = TRUE)
  version_line <- grep("version", banner, ignore.case = TRUE, value = TRUE)
  if (length(version_line) >= 1L) {
    trimws(sub("^\\*\\s*", "", version_line[1]))
  } else {
    NA_character_
  }
}, error = function(e) NA_character_)

# Provenance flag: did the SOEP chr1 Minimac info file exist when this run
# happened, or did upstream 03a fall back to the mildQC bim list for chr1?
# The SOEP provider has never shipped chr1.info.gz (see
# 01_Genotype/Harmonized/README.md "Known issues"). If they ever do, this
# field self-updates to "present" on the next run.
soep_chr1_info_path <- file.path(
  Sys.getenv("GENOTYPE_PROJ_ROOT",
             "/path/to/cluster_project_root/genotype"),
  "data", "info_files", "SOEP",
  "SOEP-G.b37.mildQC.hrc1-1_imp.chr1.info.gz"
)
soep_chr1_info <- if (file.exists(soep_chr1_info_path)) {
  "present"
} else {
  "fallback_mildQC_R2_0.1"
}

manifest <- data.frame(
  build_date  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  gcta        = gcta_version,
  proj_root   = PROJ_ROOT,
  # Cutoff recorded here so it travels with the output. Set from the environment
  # by run_pipeline.sh; this script assembles the pooled-GRM track, while the
  # per-cohort and block-diagonal tracks use GRM_CUTOFF_PER_COHORT.
  # See Plan_deviations.md §7e.
  grm_cutoff  = Sys.getenv("GRM_CUTOFF", "0.05"),
  soep_chr1_info = soep_chr1_info,
  upstream_manifest = if (nrow(upstream) > 0L) {
    paste(colnames(upstream), unlist(upstream[1, , drop = TRUE]),
          sep = "=", collapse = "; ")
  } else {
    "missing"
  },
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(FINAL_OUT, "manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Dashboard row --------------------------------------------------
# Single markdown row summarising main h² point estimates. The full table
# lives in h2_snp.tsv; this is for the README status dashboard.
fmt <- function(x, digits = 3) {
  if (is.na(x)) "NA" else sprintf(paste0("%.", digits, "f"), x)
}
get_h2 <- function(stratum) {
  ix <- which(h2_df$analysis == "main_4cell" &
              h2_df$model    == "per_stratum_standalone" &
              h2_df$stratum  == stratum)
  if (length(ix) == 0L) return(c(NA_real_, NA_real_, NA_integer_))
  c(h2_df$h2_SNP[ix[1L]], h2_df$h2_SNP_SE[ix[1L]], h2_df$N[ix[1L]])
}
dash <- paste(
  "| h²_SNP (per-stratum standalone) |",
  paste0(fmt(get_h2("East_pre1990")[1]),  " ± ", fmt(get_h2("East_pre1990")[2]),  " (N=", get_h2("East_pre1990")[3],  ")"),
  "|",
  paste0(fmt(get_h2("East_post1990")[1]), " ± ", fmt(get_h2("East_post1990")[2]), " (N=", get_h2("East_post1990")[3], ")"),
  "|",
  paste0(fmt(get_h2("West_pre1990")[1]),  " ± ", fmt(get_h2("West_pre1990")[2]),  " (N=", get_h2("West_pre1990")[3],  ")"),
  "|",
  paste0(fmt(get_h2("West_post1990")[1]), " ± ", fmt(get_h2("West_post1990")[2]), " (N=", get_h2("West_post1990")[3], ")"),
  "|"
)
writeLines(dash, file.path(FINAL_OUT, "dashboard_row.md"))
cat(sprintf("Wrote: %s\n", file.path(FINAL_OUT, "dashboard_row.md")))

cat("\n=== Step 09 complete — final outputs in ", FINAL_OUT, " ===\n", sep = "")
