#!/usr/bin/env Rscript
########################################################################
## 80 Years of GxE on Education in Germany — Merge
##
## Loads phenotypes from raw cohort shares, PGI+PC tables from the
## DATA_ROOT volume, applies cohort-specific crosswalks and uniform ID
## prefixes, stacks the four cohorts into one analysis-ready data frame
## and writes DATA_ROOT/03_Merge/Combined_Harmonized_<YYYYMMDD>.rds.
##
## RAW PGIs and ancestry PCs pass through unchanged — no residualisation
## happens here. Residualisation belongs in 04_GxE_Analysis/.
##
## Usage:
##   Rscript 03_Merge/merge.R
## Environment overrides:
##   DATA_ROOT    path to shared-volume data root
##                (default: /path/to/data_root)
########################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(haven)
})

# ---- Resolve script dir & DATA_ROOT ---------------------------------------
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
# DATA_ROOT: single source of truth in R/paths.R (override via DATA_ROOT env var).
source(file.path(SCRIPT_DIR, "R", "paths.R"), local = TRUE)
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var

# ---- Source helpers --------------------------------------------------------
source(file.path(SCRIPT_DIR, "R", "merge_common.R"), local = TRUE)

# ---- PGI loader ------------------------------------------------------------
# Reads a per-cohort {COHORT}_PGI_PCs.tsv produced by 01_Genotype/Harmonized/.
# After normalize_pgi_columns() (in merge_common.R) every trait t in
# {Edu, Cog, nonCog, Height} is represented by exactly three columns:
#   PGI_<t>           — analysis-friendly cross-cohort z
#   PGI_<t>_raw       — raw PLINK2 score (untouched)
#   PGI_<t>_within_z  — within-cohort z computed at merge time
#
# Per-cohort TSV columns:
#   IID
#   PGI_EA4, PGI_Cog, PGI_NonCog, PGI_Height                 (raw PLINK score)
#   PGI_EA4_raw, PGI_Cog_raw, PGI_NonCog_raw, PGI_Height_raw
#   PGI_EA4_global_z, PGI_Cog_global_z,
#   PGI_NonCog_global_z, PGI_Height_global_z                 (cross-cohort z)
#   PC1..PC20, crossPC1..crossPC20
load_pgi <- function(path) {
  if (!file.exists(path)) {
    warning("PGI TSV not found: ", path)
    return(NULL)
  }
  df <- read.table(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  df <- normalize_pgi_columns(df)
  df$IID <- as.character(df$IID)
  df
}

cat("========================================================================\n")
cat("  GxE Germany Merge\n")
cat(sprintf("  DATA_ROOT = %s\n", DATA_ROOT))
cat("========================================================================\n\n")

# Prerequisite: the cohort phenotype shares and the project data share must be
# mounted before this script runs. It reads them directly and does not mount.

# ---- Locate PGI TSVs on the shared volume ---------------------------------
PGI_DIR <- file.path(DATA_ROOT, "01_Genotype", "Harmonized", "data", "final")
if (!dir.exists(PGI_DIR))
  stop("PGI directory not found: ", PGI_DIR,
       "\nTransfer per-cohort {COHORT}_PGI_PCs.tsv files from Tardis to this folder first.")

tsv_files <- list(
  BASE   = file.path(PGI_DIR, "BASEII_PGI_PCs.tsv"),
  SHIP0  = file.path(PGI_DIR, "SHIP0_PGI_PCs.tsv"),
  SHIPTD = file.path(PGI_DIR, "SHIPTD_PGI_PCs.tsv"),
  SOEP   = file.path(PGI_DIR, "SOEP_PGI_PCs.tsv"),
  TWIN   = file.path(PGI_DIR, "TWINLIFE_PGI_PCs.tsv")
)
missing <- tsv_files[!sapply(tsv_files, file.exists)]
if (length(missing) > 0)
  stop("Missing PGI TSV file(s):\n  ",
       paste(unlist(missing), collapse = "\n  "))
cat(sprintf("\nPGI dir: %s\n", PGI_DIR))

# ---- 1. Load phenotypes ----------------------------------------------------
cat("\n--- 1. Loading phenotype files ---\n")
df_base <- load_pheno_base();     cat(sprintf("  BASE-II:  N = %d\n", nrow(df_base)))
df_ship <- load_pheno_ship();     cat(sprintf("  SHIP:     N = %d\n", nrow(df_ship)))
df_soep <- load_pheno_soep();     cat(sprintf("  SOEP:     N = %d\n", nrow(df_soep)))
df_twin <- load_pheno_twinlife(); cat(sprintf("  TwinLife: N = %d\n", nrow(df_twin)))

# ---- 2. Load PGI+PC tables, ancestry-filter, crosswalk --------------------
cat("\n--- 2. Loading PGI+PC tables and merging ---\n")

cat("\n[BASE-II]\n")
pgi_base <- load_pgi(tsv_files$BASE) |> filter_to_eur("BASEII", DATA_ROOT)
df_base <- baseii_crosswalk(df_base, pgi_base)
report_merge(df_base, "BASE-II")

cat("\n[SHIP]\n")
pgi_s0  <- load_pgi(tsv_files$SHIP0)  |> filter_to_eur("SHIP0",  DATA_ROOT)
pgi_std <- load_pgi(tsv_files$SHIPTD) |> filter_to_eur("SHIPTd", DATA_ROOT)
df_ship <- ship_crosswalk(df_ship, pgi_s0, pgi_std)
report_merge(df_ship, "SHIP")

cat("\n[SOEP]\n")
pgi_soep <- load_pgi(tsv_files$SOEP) |> filter_to_eur("SOEP", DATA_ROOT)
df_soep <- soep_crosswalk(df_soep, pgi_soep)
report_merge(df_soep, "SOEP")

cat("\n[TwinLife]\n")
pgi_twin <- load_pgi(tsv_files$TWIN) |> filter_to_eur("TwinLife", DATA_ROOT)
df_twin <- twinlife_crosswalk(df_twin, pgi_twin)
report_merge(df_twin, "TwinLife")
twin_correlation_check(df_twin)

# ---- 3. Apply uniform prefixes & stack ------------------------------------
cat("\n--- 3. Harmonising column names + stacking cohorts ---\n")

# `across(where(haven::is.labelled), as.numeric)` applied to every cohort:
# any haven_labelled (SPSS/Stata) column that survived from the cohort's
# DataMining script gets coerced to plain numeric so bind_rows() doesn't
# complain about mixed labelled/non-labelled column types.
#
# pid prefix policy:
#   BASE-II  — pheno already carries `b_<nkidpz>` (added upstream); no edit.
#   SHIP     — pheno carries `s_<id>` (SHIP-0) or `t_<id>` (SHIP-TREND);
#              both are intentionally collapsed to `sh_<id>` here so the
#              two SHIP sub-samples appear as a single SHIP cohort.
#   SOEP     — pheno carries the bare numeric ID; we add `so_`.
#   TwinLife — pheno carries the bare numeric ID; we add `tl_`.

df_base_out <- df_base %>%
  mutate(cohort = "BASE-II",
         across(where(haven::is.labelled), as.numeric)) %>%
  select(any_of(COMMON_COLS))

df_ship_out <- df_ship %>%
  mutate(cohort = "SHIP",
         pid    = sub("^[st]_", "sh_", pid),
         across(where(haven::is.labelled), as.numeric)) %>%
  select(any_of(COMMON_COLS))

df_soep_out <- df_soep %>%
  mutate(cohort = "SOEP",
         pid    = paste0("so_", pid),
         across(where(haven::is.labelled), as.numeric)) %>%
  select(any_of(COMMON_COLS))

df_twin_out <- df_twin %>%
  mutate(cohort = "TwinLife",
         pid    = paste0("tl_", pid),
         across(where(haven::is.labelled), as.numeric)) %>%
  select(any_of(COMMON_COLS))

df_combined <- bind_rows(df_base_out, df_ship_out, df_soep_out, df_twin_out)

# ---- 4. Coverage diagnostics ----------------------------------------------
cat(sprintf("\nCombined: %d rows x %d columns\n",
            nrow(df_combined), ncol(df_combined)))
cat("Cohort breakdown:\n")
print(table(df_combined$cohort, useNA = "ifany"))
cat("\nPGI_Edu coverage per cohort:\n")
print(tapply(!is.na(df_combined$PGI_Edu), df_combined$cohort,
             function(x) sprintf("%d / %d (%.1f%%)",
                                 sum(x), length(x), 100 * mean(x))))

# ---- 5. Save ---------------------------------------------------------------
cat("\n--- 5. Saving combined RDS ---\n")
save_combined(df_combined, DATA_ROOT, pgi_dir = PGI_DIR)

cat("\nDone.\n")
