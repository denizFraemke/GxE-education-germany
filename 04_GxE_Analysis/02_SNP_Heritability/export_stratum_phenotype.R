#!/usr/bin/env Rscript
########################################################################
## export_stratum_phenotype.R — Mac-side phenotype + stratum exporter
##
## Reads the most recent Combined_Harmonized_<YYYYMMDD>.rds from the
## shared volume, derives the Region x Time stratification used in
## Analysis S5 (Region from `east_west`; Time from `birth_year` via the
## "turned-15-in-or-after-1990" rule, i.e. BY >= 1975 = post), and
## writes GCTA-formatted inputs to
##   ${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/data/
## ready to be rsynced to Tardis by transfer_phenotype_to_tardis.sh.
##
## RUNS ON: Mac with the shared SMB volume mounted.
##
## Outputs (all two-column "FID IID ..." files, GCTA-compatible):
##   data/education.phen         — FID IID education (raw years)
##   data/crosspcs.qcovar        — FID IID crossPC1..20 birth_year_c birth_year_c_sq
##                                 (20 cross-cohort PCs + age covariate; see
##                                 Plan_deviations.md §7f for the BY rationale)
##   data/gender.covar           — FID IID gender (categorical; pooled-GRM
##                                 and per-cohort tracks)
##   data/gender_cohort.covar    — FID IID gender + cohort_int (categorical;
##                                 cell-level block-diagonal primary + R1.
##                                 cohort_int: 1=BASE-II, 2=SHIP, 3=SOEP,
##                                 4=TwinLife; see §7g.)
##   data/region_cohort.covar    — FID IID region_int + cohort_int
##                                 (block-diagonal G1 only; region_int:
##                                 1=East, 2=West; see §7h.)
##   data/cohort_only.covar      — FID IID cohort_int
##                                 (block-diagonal RG1 only; see §7h.)
##   data/keep/<stratum>.iids    — FID IID for each main cell
##   data/manifest.tsv           — provenance (source RDS, build date, cell N,
##                                 BY centring constant, qcovar schema)
##
## Aggregate-only stdout: dimensions and per-cell counts. No row prints.
##
## Usage:
##   Rscript 04_GxE_Analysis/02_SNP_Heritability/export_stratum_phenotype.R
## Environment overrides:
##   DATA_ROOT   path to shared-volume Education_Genomics/ root
########################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
})

# ---- Resolve paths ---------------------------------------------------------
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})

# DATA_ROOT: single source of truth in R/paths.R (override via DATA_ROOT env var).
source(file.path(SCRIPT_DIR, "R", "paths.R"), local = TRUE)
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var

MERGE_DIR <- file.path(DATA_ROOT, "03_Merge")
OUT_DIR   <- file.path(DATA_ROOT, "04_GxE_Analysis", "02_SNP_Heritability", "data")
KEEP_DIR  <- file.path(OUT_DIR, "keep")
dir.create(KEEP_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================================================\n")
cat("  S5 SNP-heritability — stratum + phenotype export (Mac side)\n")
cat(sprintf("  DATA_ROOT: %s\n", DATA_ROOT))
cat(sprintf("  Out dir:   %s\n", OUT_DIR))
cat("========================================================================\n\n")

# ---- Find most recent Combined_Harmonized_<YYYYMMDD>.rds -------------------
rds_candidates <- list.files(MERGE_DIR,
                             pattern = "^Combined_Harmonized_\\d{8}\\.rds$",
                             full.names = TRUE)
if (length(rds_candidates) == 0L)
  stop("No Combined_Harmonized_<date>.rds files in ", MERGE_DIR)
rds_path <- rds_candidates[which.max(as.integer(regmatches(
  rds_candidates, regexpr("\\d{8}", rds_candidates))))]
cat(sprintf("Source RDS: %s\n", basename(rds_path)))

df <- readRDS(rds_path)
cat(sprintf("Loaded:     %d rows x %d cols\n\n", nrow(df), ncol(df)))

# ---- Required columns ------------------------------------------------------
# IID is the genotype-side identifier (e.g. "CO12345", "s_56789",
# "CP00042677_CP00042677", "DEXxxxxxxxxxxx") — the one that lives in the
# per-cohort harmonized bfile's .fam. Required for GCTA --keep to actually
# match rows in the pooled bfile (using `pid` here would silently produce
# 0-row keep files and an empty GRM). It is preserved through
# merge_common.R's per-cohort crosswalks.
required <- c("pid", "IID", "cohort", "east_west",
              "birth_year", "gender", "education")
missing_cols <- setdiff(required, names(df))
if (length(missing_cols) > 0)
  stop("Required columns missing from merged RDS: ",
       paste(missing_cols, collapse = ", "),
       "\nRe-run 03_Merge/merge.R against the latest origin/main to regenerate",
       " a Combined_Harmonized_<date>.rds that carries the IID column.")
# crossPCs: prefer crossPC1..20 (Harmonized mode); fail loudly if absent.
cross_pcs <- paste0("crossPC", 1:20)
missing_pcs <- setdiff(cross_pcs, names(df))
if (length(missing_pcs) > 0)
  stop("Cross-cohort PC columns missing: ",
       paste(missing_pcs, collapse = ", "))

# ---- Filter to analytic sample ---------------------------------------------
# Plan S5 needs phenotype + region + birth year + ancestry covariates +
# gender, plus the genotype-side IID so the row can actually enter REML.
analytic <- df %>%
  filter(!is.na(IID),
         !is.na(east_west),
         !is.na(birth_year),
         !is.na(education),
         !is.na(gender),
         if_all(all_of(cross_pcs), ~ !is.na(.x)))

cat(sprintf("Analytic sample (non-missing east_west/birth_year/education/gender/crossPC1..20): N = %d\n",
            nrow(analytic)))
cat("Cohort breakdown of analytic sample:\n")
print(table(analytic$cohort, useNA = "ifany"))
cat("\n")

# ---- Assemble FID/IID -------------------------------------------------------
# GCTA expects two-column "FID IID" files. The pooled.fam on Tardis carries
# the original cohort-genotype-side IID (e.g. "CO12345", "s_56789",
# "CP00042677_CP00042677", "DEXxxxxxxxxxxx") and a cohort-specific FID. The
# exporter writes the IID we have here and duplicates it into the FID slot
# as a placeholder; the Tardis-side scripts/03_build_stratum_grms.sh resolves
# the actual FID by joining each stratum's IIDs against pooled.fam before
# calling gcta64 --keep.
analytic <- analytic %>%
  mutate(FID = IID)   # IID column already comes from the merged RDS

# ---- Derive strata ---------------------------------------------------------
# Region: east_west is character {"east","west"} in the merged RDS.
# Time (4-cell):  BY < 1975  -> pre1990   (turned 15 before 1990)
#                 BY >= 1975 -> post1990  (turned 15 in/after 1990)
# Recorded gender/sex: harmonized across cohorts via the existing
# Male/Female mapping below. The variable likely reflects each source
# cohort's recorded sex at recruitment or self-reported gender — see
# README.md / Plan_deviations.md §7h for the data-dictionary caveat.
analytic$region <- ifelse(analytic$east_west == "east", "East", "West")
time_main       <- ifelse(analytic$birth_year < 1975, "pre1990", "post1990")
analytic$stratum_main <- paste(analytic$region, time_main, sep = "_")
analytic$gender_norm  <- dplyr::case_when(
  as.character(analytic$gender) %in% c("Male",   "male",   "M", "1", "m") ~ "Male",
  as.character(analytic$gender) %in% c("Female", "female", "F", "2", "f") ~ "Female",
  TRUE ~ NA_character_
)

cat("MAIN 4-cell N (post non-missing filter):\n")
print(addmargins(table(region = analytic$region, time_main = time_main)))
cat("\n")
cat("Region marginal:\n");       print(table(analytic$region))
cat("Recorded gender marginal:\n"); print(table(analytic$gender_norm, useNA = "ifany"))
cat("Region × recorded gender:\n")
print(addmargins(table(region = analytic$region, gender = analytic$gender_norm,
                       useNA = "ifany")))
cat("\n")

# ---- Age covariate --------------------------------------------------------
# Plan_deviations.md §7f: the four cells span 25–55 years of birth-year
# variation. Without controlling for BY inside each cell, secular cohort
# trends in education load onto V_e and bias h² downward. Add BY (centred
# on the analytic-sample mean) and its square to the qcovar so REML
# partials out a smooth quadratic age trend per stratum. h² is invariant
# to fixed-effect addition; what changes is the cell-level recalibration
# when within-cell cohort variation no longer hides in V_e.
#
# Centring matters only for numerical conditioning of the quadratic
# (raw BY vs BY² are highly collinear; centred BY vs centred BY² are not),
# not for the partial-out itself.
BY_MEAN <- mean(analytic$birth_year)
analytic <- analytic %>%
  mutate(birth_year_c    = birth_year - BY_MEAN,
         birth_year_c_sq = (birth_year - BY_MEAN)^2)

cat(sprintf("Age covariate: birth_year centred at analytic-sample mean = %.2f.\n",
            BY_MEAN))
cat("Within-cell BY summary (aggregate; for sanity-checking):\n")
by_summary <- aggregate(birth_year ~ stratum_main, data = analytic,
                        FUN = function(x) c(N = length(x),
                                            min = min(x),
                                            mean = round(mean(x), 1),
                                            max = max(x),
                                            sd = round(sd(x), 2)))
print(by_summary)
cat("\n")

# ---- Write GCTA phenotype + covariate files --------------------------------
write_tsv_nohdr <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# *.phen — FID IID y. Education is raw years per the plan; no z-standardization
# happens here (we want V_G interpretable in years^2).
phen_df <- analytic %>% select(FID, IID, education)
write_tsv_nohdr(phen_df, file.path(OUT_DIR, "education.phen"))

# *.qcovar — FID IID crossPC1..20 birth_year_c birth_year_c_sq.
# Plan_deviations.md §7f explains the BY addition. The filename says
# "crosspcs" but the schema is broader; it matches config.sh's QCOVAR_FILE.
qcovar_cols <- c(cross_pcs, "birth_year_c", "birth_year_c_sq")
qcovar_df <- analytic %>% select(FID, IID, all_of(qcovar_cols))
write_tsv_nohdr(qcovar_df, file.path(OUT_DIR, "crosspcs.qcovar"))

# *.covar — FID IID gender (categorical covariate; gcta treats integer
# levels as factors). The merged-RDS `gender` is a character column with
# values "Male" / "Female" (cf. 03_Merge crosswalk), so we map to the
# PLINK / GCTA .fam convention: 1 = Male, 2 = Female. Anything else
# becomes NA and is silently dropped by gcta --covar.
covar_df <- analytic %>%
  mutate(gender_int = dplyr::case_when(
    as.character(gender) %in% c("Male",   "male",   "M", "1", "m") ~ 1L,
    as.character(gender) %in% c("Female", "female", "F", "2", "f") ~ 2L,
    TRUE ~ NA_integer_
  )) %>%
  select(FID, IID, gender_int)

n_total <- nrow(covar_df)
n_na    <- sum(is.na(covar_df$gender_int))
if (n_na > 0L) {
  cat(sprintf("  WARN: %d/%d rows have unmappable gender — written as NA\n",
              n_na, n_total))
}
write_tsv_nohdr(covar_df, file.path(OUT_DIR, "gender.covar"))

# ---- region_int / cohort_int derived columns ------------------------------
# Stable categorical coding for the secondary GREML comparison tracks
# (Plan_deviations.md §7h). region_int and cohort_int are written into
# the dedicated covar files below.
#   region_int: 1 = East, 2 = West
#   cohort_int: 1 = BASE-II, 2 = SHIP, 3 = SOEP, 4 = TwinLife
analytic <- analytic %>%
  mutate(region_int = dplyr::case_when(
    region == "East" ~ 1L,
    region == "West" ~ 2L,
    TRUE             ~ NA_integer_
  ),
  cohort_int = dplyr::case_when(
    cohort == "BASE-II"  ~ 1L,
    cohort == "SHIP"     ~ 2L,
    cohort == "SOEP"     ~ 3L,
    cohort == "TwinLife" ~ 4L,
    TRUE                 ~ NA_integer_
  ))

# ---- gender_cohort.covar (block-diagonal track only; see Plan_deviations §7g)
# Adds a categorical `cohort_int` column alongside `gender_int`. GCTA --covar
# dummy-codes both. Used ONLY by step 05d (block-diagonal cell-level REML).
# The cell-level pooled (steps 03-05) and per-cohort sensitivity (03c-06c)
# tracks continue to use the plain gender.covar — cohort fixed effects don't
# add information there (per-cohort track has constant cohort within a fit;
# cell-level pooled GRM has cross-cohort entries that absorb the cohort means
# differently).
# Stable cohort coding:
#   BASE-II = 1, SHIP = 2, SOEP = 3, TwinLife = 4
covar_cohort_df <- covar_df %>%
  mutate(cohort_int = dplyr::case_when(
    analytic$cohort == "BASE-II"  ~ 1L,
    analytic$cohort == "SHIP"     ~ 2L,
    analytic$cohort == "SOEP"     ~ 3L,
    analytic$cohort == "TwinLife" ~ 4L,
    TRUE                          ~ NA_integer_
  )) %>%
  select(FID, IID, gender_int, cohort_int)

n_cohort_na <- sum(is.na(covar_cohort_df$cohort_int))
if (n_cohort_na > 0L) {
  cat(sprintf("  WARN: %d/%d rows have unmappable cohort — written as NA\n",
              n_cohort_na, nrow(covar_cohort_df)))
}
write_tsv_nohdr(covar_cohort_df, file.path(OUT_DIR, "gender_cohort.covar"))

# ---- region_cohort.covar (G1: recorded gender/sex main comparison) --------
# Block-diagonal track G1 uses this covar in --covar. region_int + cohort_int
# are non-constant within each gender stratum; gender_int is constant (and
# therefore excluded). See Plan_deviations.md §7h.
covar_region_cohort_df <- covar_df %>%
  left_join(analytic %>% select(FID, IID, region_int, cohort_int),
            by = c("FID", "IID")) %>%
  select(FID, IID, region_int, cohort_int)
write_tsv_nohdr(covar_region_cohort_df, file.path(OUT_DIR, "region_cohort.covar"))

# ---- cohort_only.covar (RG1: Region × recorded gender/sex comparison) ----
# Block-diagonal track RG1 uses this covar in --covar. cohort_int is
# non-constant within each region × gender stratum; both region_int and
# gender_int are constant and excluded. See Plan_deviations.md §7h.
covar_cohort_only_df <- covar_df %>%
  left_join(analytic %>% select(FID, IID, cohort_int),
            by = c("FID", "IID")) %>%
  select(FID, IID, cohort_int)
write_tsv_nohdr(covar_cohort_only_df, file.path(OUT_DIR, "cohort_only.covar"))

cat(sprintf("Wrote: %s  (N=%d)\n", "education.phen",    nrow(phen_df)))
cat(sprintf("Wrote: %s  (N=%d, cols: %s)\n", "crosspcs.qcovar",
            nrow(qcovar_df), paste(qcovar_cols, collapse = ", ")))
cat(sprintf("Wrote: %s  (N=%d, cols: gender_int)\n", "gender.covar",
            nrow(covar_df)))
cat(sprintf("Wrote: %s  (N=%d, cols: gender_int + cohort_int [1=BASE-II,2=SHIP,3=SOEP,4=TwinLife])\n",
            "gender_cohort.covar", nrow(covar_cohort_df)))
cat(sprintf("Wrote: %s  (N=%d, cols: region_int [1=East,2=West] + cohort_int)\n",
            "region_cohort.covar", nrow(covar_region_cohort_df)))
cat(sprintf("Wrote: %s  (N=%d, cols: cohort_int)\n",
            "cohort_only.covar", nrow(covar_cohort_only_df)))

# ---- Write per-stratum FID/IID keep files ----------------------------------
write_keep <- function(stratum_label, mask) {
  keep <- analytic[mask, c("FID", "IID")]
  out  <- file.path(KEEP_DIR, paste0(stratum_label, ".iids"))
  write_tsv_nohdr(keep, out)
  invisible(nrow(keep))
}

main_levels <- c("East_pre1990", "East_post1990",
                 "West_pre1990", "West_post1990")

cat("\nKeep files (data/keep/<stratum>.iids):\n")
cat("  Primary cell-level (Region × Reunification):\n")
for (lev in main_levels) {
  n <- write_keep(lev, analytic$stratum_main == lev)
  cat(sprintf("    %-22s N=%d\n", lev, n))
}

# Secondary/exploratory tracks (Plan_deviations.md §7h):
#   R1  region_main         — East vs West, collapsed across pre/post-1990
#   G1  gender_main         — Male vs Female (recorded gender/sex), collapsed
#                             across Region and time
#   RG1 region_gender       — 4-group Region × recorded gender/sex, collapsed
#                             across pre/post-1990
cat("  R1 Region main (collapsed across pre/post-1990):\n")
n <- write_keep("region_east", analytic$region == "East")
cat(sprintf("    %-22s N=%d\n", "region_east", n))
n <- write_keep("region_west", analytic$region == "West")
cat(sprintf("    %-22s N=%d\n", "region_west", n))

cat("  G1 recorded gender/sex main (collapsed across Region and time):\n")
n <- write_keep("gender_male",   analytic$gender_norm == "Male")
cat(sprintf("    %-22s N=%d\n", "gender_male", n))
n <- write_keep("gender_female", analytic$gender_norm == "Female")
cat(sprintf("    %-22s N=%d\n", "gender_female", n))

cat("  RG1 Region × recorded gender/sex (collapsed across pre/post-1990):\n")
for (r in c("East", "West")) {
  for (g in c("Male", "Female")) {
    lab <- sprintf("region_gender_%s_%s", tolower(r), tolower(g))
    n   <- write_keep(lab,
                      analytic$region == r & analytic$gender_norm == g)
    cat(sprintf("    %-30s N=%d\n", lab, n))
  }
}

# ---- Manifest --------------------------------------------------------------
manifest <- data.frame(
  source_rds      = basename(rds_path),
  source_rds_mtime = format(file.info(rds_path)$mtime, "%Y-%m-%dT%H:%M:%S"),
  build_date      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  n_analytic      = nrow(analytic),
  reunification_cutoff_BY = 1975,
  birth_year_mean = sprintf("%.4f", BY_MEAN),
  qcovar_columns  = paste(qcovar_cols, collapse = ","),
  covar_columns   = "gender_int",
  blockdiag_covar_columns = "gender_int,cohort_int",
  region_cohort_covar_columns = "region_int,cohort_int",
  cohort_only_covar_columns   = "cohort_int",
  region_int_coding = "1=East,2=West",
  gender_int_coding = "1=Male,2=Female",
  cohort_int_coding = "1=BASE-II,2=SHIP,3=SOEP,4=TwinLife",
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(OUT_DIR, "manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Next step: bash transfer_phenotype_to_tardis.sh\n")
