########################################################################
## 80 Years of GxE on Education in Germany
## merge_common.R — shared helpers for 03_Merge/merge.R
##
## Responsibilities:
##   - SMB auto-mount of required volumes (macOS `open smb://…`).
##   - Per-cohort phenotype loaders (from raw cohort shares).
##   - Cohort-specific crosswalks (BASE-II, SHIP, SOEP, TwinLife).
##   - Ancestry-outlier filter on ancestry PCs.
##   - Uniform ID prefixing and the combined-rds save routine.
##   - PGI-column normalisation (maps PGI_EA4/PGI_NonCog -> PGI_Edu/PGI_nonCog;
##     promotes the cross-cohort _global_z value to the bare PGI_<trait>
##     name; always computes a real within-cohort z and exposes it as
##     *_within_z; keeps the raw PLINK2 score as *_raw).
##   - Twin-pair PGI correlation diagnostic.
##   - COMMON_COLS vector for the final dataset.
##
## No residualisation happens here. PGIs and PCs flow through unchanged.
########################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(haven)
})

# ------------------------------------------------------------------
# Phenotype input paths (raw phenotype stays on cohort shares; the
# loaders below pick the file with the most recent YYYYMMDD in its
# filename so we don't have to bump a hardcoded date every delivery).
# ------------------------------------------------------------------
PHENO_BASE_DIR <- "${BASE2_SHARE}/private/data/004_BASEII_Socioeceonomic_Pheno/Final_data"
PHENO_SHIP_DIR <- "${SHIP_SHARE}/private/data/007_SHIP_processed_pheno"
PHENO_SOEP_DIR <- "${SOEP_SHARE}/private/data/001_SOEP_processed_data/007_GermanGenetics"
PHENO_TWIN_DIR <- "${TWINLIFE_SHARE}/private/data/2025_TwinLife_GermanGenetics"

BASEII_KEY <- "${BASE2_SHARE}/private/data/002_BASEII_pheno/baseII_ID_key.txt"
TWIN_XWALK <- "${TWINLIFE_SHARE}/private/data/2025_TwinLife_GermanGenetics/twinlife_pid_to_genetic_iid_crosswalk.csv"

# ------------------------------------------------------------------
# Phenotype loaders — each returns a data.frame
#
# We resolve the actual file at runtime by picking the most recent
# YYYYMMDD-stamped file matching the cohort's prefix in the cohort's
# share directory. The loaded .rda contains a single data.frame with
# a cohort-specific name (df_final / df_SHIP / soepis_final / df_final)
# — we just grab whatever object is in the file.
# ------------------------------------------------------------------
.most_recent_rda <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L)
    stop("No file matching '", pattern, "' in:\n  ", dir)
  date_str <- regmatches(files, regexpr("\\d{8}", files))
  if (length(date_str) != length(files) || any(date_str == "")) {
    return(files[which.max(file.info(files)$mtime)])  # mtime fallback
  }
  files[which.max(as.integer(date_str))]
}

.load_first_object <- function(path) {
  if (!file.exists(path)) stop("Phenotype file not found: ", path)
  e <- new.env()
  load(path, envir = e)
  obj_names <- ls(e)
  if (length(obj_names) == 0L) stop("Empty .rda: ", path)
  e[[obj_names[1L]]]
}

load_pheno_base <- function() .load_first_object(
  .most_recent_rda(PHENO_BASE_DIR, "^EduGxE_DE_BASEII_Phenotype_\\d{8}\\.rda$")
)
load_pheno_ship <- function() .load_first_object(
  .most_recent_rda(PHENO_SHIP_DIR, "^EduGxE_DE_SHIP_Phenotype_\\d{8}\\.rda$")
)
load_pheno_soep <- function() .load_first_object(
  .most_recent_rda(PHENO_SOEP_DIR, "^EduGxE_DE_SOEP_Phenotype_\\d{8}\\.rda$")
)
load_pheno_twinlife <- function() .load_first_object(
  .most_recent_rda(PHENO_TWIN_DIR, "^EduGxE_DE_TwinLife_Phenotype_\\d{8}\\.(rda|RData)$")
)

# ------------------------------------------------------------------
# PGI column normalisation
#
# Standard internal trait labels (TSV-name aliases):
#   PGI_EA4    -> PGI_Edu
#   PGI_NonCog -> PGI_nonCog
#
# After normalisation, every trait t in {Edu, Cog, nonCog, Height} is
# represented by exactly THREE columns:
#
#   PGI_<t>             — analysis-friendly cross-cohort z (was
#                         PGI_<t>_global_z in the input TSV).
#   PGI_<t>_raw         — raw PLINK2 score (untouched).
#   PGI_<t>_within_z    — TRUE within-cohort z computed here from the
#                         raw column: (raw - cohort_mean) / cohort_sd.
#                         Always populated.
#
# Note: this function runs PER COHORT (one PGI table at a time), which
# is why within-cohort z computed here is genuinely cohort-local.
# ------------------------------------------------------------------
normalize_pgi_columns <- function(df) {
  nm <- names(df)
  nm <- sub("PGI_EA4",    "PGI_Edu",    nm, fixed = TRUE)
  nm <- sub("PGI_NonCog", "PGI_nonCog", nm, fixed = TRUE)
  names(df) <- nm

  for (t in c("Edu", "Cog", "nonCog", "Height")) {
    plain <- paste0("PGI_", t)
    raw   <- paste0("PGI_", t, "_raw")
    wz    <- paste0("PGI_", t, "_within_z")
    gz    <- paste0("PGI_", t, "_global_z")

    if (!plain %in% names(df)) next  # trait absent from this cohort's TSV

    # Capture raw values: use _raw if present, otherwise the bare column.
    raw_vec <- if (raw %in% names(df)) df[[raw]] else df[[plain]]

    # Compute REAL within-cohort z from the raw column.
    x <- as.numeric(raw_vec)
    cohort_mean <- mean(x, na.rm = TRUE)
    cohort_sd   <- sd(x,   na.rm = TRUE)
    df[[wz]] <- if (is.na(cohort_sd) || cohort_sd == 0) {
      rep(NA_real_, length(x))
    } else {
      (x - cohort_mean) / cohort_sd
    }

    # Ensure _raw column always exists.
    df[[raw]] <- raw_vec

    # Promote cross-cohort global_z to the bare PGI_<t> name and drop
    # the now-redundant _global_z column.
    if (!gz %in% names(df))
      stop("PGI TSV missing required column '", gz,
           "'. Expected from 01_Genotype/Harmonized/ output.")
    df[[plain]] <- df[[gz]]
    df[[gz]] <- NULL
  }
  df
}

# ------------------------------------------------------------------
# Ancestry-outlier filter: drop rows where any |z| on PC1..PC10 > cutoff
# ------------------------------------------------------------------
filter_ancestry_outliers <- function(pgi_df, pc_cols = paste0("PC", 1:10),
                                     z_cutoff = 3) {
  if (is.null(pgi_df)) return(NULL)
  pc_cols <- intersect(pc_cols, names(pgi_df))
  if (length(pc_cols) == 0) return(pgi_df)
  keep <- rep(TRUE, nrow(pgi_df))
  for (pc in pc_cols) {
    v <- as.numeric(pgi_df[[pc]])
    z <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    keep <- keep & (is.na(z) | abs(z) <= z_cutoff)
  }
  dropped <- sum(!keep, na.rm = TRUE)
  if (dropped > 0)
    message(sprintf("    dropped %d ancestry outlier(s) (|z|>%d on PC1..10)",
                    dropped, z_cutoff))
  pgi_df[keep, , drop = FALSE]
}

# ------------------------------------------------------------------
# 1kG-EUR keep-list filter.
# Keeps only individuals on the per-cohort European keep-list, read from
# data/1Kg_exclusion_IDs/EUR_keep_IDs_<cohort_key>.tsv.
#
# Keep-list provenance:
#   BASE-II / SHIP / TwinLife — classified in-house by projecting each cohort
#     onto principal components trained on the 1000 Genomes phase-3 reference
#     and keeping European-assigned individuals (see
#     01_Genotype/Harmonized/scripts/04c_ancestry_1kg.sh + 04c_classify_eur.R).
#   SOEP — a European-ancestry keep-list computed for SOEP in a separate
#     project; the SOEP genotypes are not available to re-project, so this
#     list is used as-is.
#
# SOEP special case: that keep-list was built on the strict-QC (chrall) release,
# whereas our scored SOEP is the mild-QC (ext_chrall) release. We keep
# (on-list EUR) UNION (ext-only individuals absent from chrall, i.e. mild-QC
# only, unassessable and ~98% European per the SOEP-G cohort profile) and drop
# only the chrall individuals classified non-European.
# Adds a geno_qc column ("strict" / "mild_only").
# ------------------------------------------------------------------
EUR_KEEP_DIR <- function(data_root)
  file.path(data_root, "01_Genotype", "Harmonized", "data", "1Kg_exclusion_IDs")

filter_to_eur <- function(pgi_df, cohort_key, data_root) {
  if (is.null(pgi_df)) return(NULL)
  # NO_EXCL=1 -> pass-through (no ancestry/PC exclusion at all), for isolation tests.
  if (identical(Sys.getenv("NO_EXCL"), "1")) {
    pgi_df$geno_qc <- "none"
    message(sprintf("    [%s] NO_EXCL: keeping all %d (no ancestry exclusion)", cohort_key, nrow(pgi_df)))
    return(pgi_df)
  }
  eur_dir  <- EUR_KEEP_DIR(data_root)
  eur_file <- file.path(eur_dir, sprintf("EUR_keep_IDs_%s.tsv", cohort_key))
  if (!file.exists(eur_file)) stop("EUR keep-list not found: ", eur_file)
  eur_iid  <- as.character(read.table(eur_file, header = TRUE, sep = "\t")$IID)
  pgi_df$IID <- as.character(pgi_df$IID)

  if (cohort_key == "SOEP") {
    gid <- function(x) as.integer(sub("^CP0*", "", sub("_.*$", "", as.character(x))))
    strict_file <- file.path(eur_dir, "SOEP_strict_chrall_IIDs.tsv")
    if (!file.exists(strict_file)) stop("SOEP strict-QC id list not found: ", strict_file)
    strict_gid <- gid(read.table(strict_file, header = TRUE, sep = "\t")$IID)
    eur_gid    <- gid(eur_iid)
    g          <- gid(pgi_df$IID)
    in_chrall  <- g %in% strict_gid
    keep       <- (g %in% eur_gid) | !in_chrall        # EUR OR ext-only (mild-QC)
    pgi_df$geno_qc <- ifelse(in_chrall, "strict", "mild_only")
  } else {
    keep <- pgi_df$IID %in% eur_iid
    pgi_df$geno_qc <- "strict"
  }
  message(sprintf("    [%s] 1kG-EUR filter: kept %d, dropped %d (of %d)",
                  cohort_key, sum(keep), sum(!keep), nrow(pgi_df)))
  pgi_df[keep, , drop = FALSE]
}

# ------------------------------------------------------------------
# Crosswalks
# ------------------------------------------------------------------

# BASE-II: {CO,DA,SP,TS}XXXX -> {22,21,23,24}XXXX -> nkidpznum
baseii_crosswalk <- function(df_pheno, pgi_df) {
  if (is.null(pgi_df)) return(df_pheno)
  key <- read.table(BASEII_KEY, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,
                    col.names = c("genetic_iid", "nkidpz"))
  pfx_map <- c(CO = "22", DA = "21", SP = "23", TS = "24")
  key$nkidpznum <- NA_integer_
  for (pfx in names(pfx_map)) {
    idx <- startsWith(key$nkidpz, pfx)
    key$nkidpznum[idx] <- as.integer(paste0(pfx_map[[pfx]],
                                            substring(key$nkidpz[idx], 3)))
  }
  df_pheno$nkidpznum <- as.integer(sub("^b_", "", df_pheno$pid))
  df_pheno %>%
    left_join(key[, c("genetic_iid", "nkidpznum")], by = "nkidpznum") %>%
    left_join(pgi_df, by = c("genetic_iid" = "IID")) %>%
    # Preserve the genotype-side IID (the join consumed pgi_df$IID into
    # genetic_iid). Needed by downstream pipelines that touch the per-cohort
    # bfiles, e.g. 04_GxE_Analysis/02_SNP_Heritability/.
    rename(IID = genetic_iid) %>%
    select(-nkidpznum)
}

# SHIP: last 5 chars of IID, prefixed s_ / t_, matched to pheno pid
ship_crosswalk <- function(df_pheno, pgi_0, pgi_td) {
  parse_prefix <- function(pgi_df, pfx) {
    if (is.null(pgi_df)) return(NULL)
    pgi_df$IID <- as.character(pgi_df$IID)
    pgi_df$ship_pid <- paste0(pfx,
                              substr(pgi_df$IID,
                                     nchar(pgi_df$IID) - 4,
                                     nchar(pgi_df$IID)))
    pgi_df
  }
  pgi <- bind_rows(parse_prefix(pgi_0, "s_"),
                   parse_prefix(pgi_td, "t_")) %>%
    distinct(ship_pid, .keep_all = TRUE)
  # Preserve pgi_df$IID (s_<lastfive> / t_<lastfive>): needed by downstream
  # pipelines that touch the per-cohort bfiles.
  df_pheno %>% left_join(pgi, by = c("pid" = "ship_pid"))
}

# SOEP: IID is "CP00042677_CP00042677"; strip to 42677 -> matches pheno ID
soep_crosswalk <- function(df_pheno, pgi_df) {
  if (is.null(pgi_df)) return(df_pheno)
  pgi_df$genetic_id <- as.integer(sub("^CP0+", "",
                                      sub("_.*$", "", pgi_df$IID)))
  # Preserve pgi_df$IID ("CPxxxxxxxx_CPxxxxxxxx"): needed by downstream
  # pipelines that touch the per-cohort bfiles.
  df_pheno %>% left_join(pgi_df, by = c("ID" = "genetic_id"))
}

# TwinLife: pid -> genetic_iid via CSV crosswalk
twinlife_crosswalk <- function(df_pheno, pgi_df) {
  if (is.null(pgi_df)) return(df_pheno)
  xwalk <- read.csv(TWIN_XWALK, stringsAsFactors = FALSE)
  df_pheno %>%
    left_join(xwalk[, c("pid", "genetic_iid")], by = "pid") %>%
    left_join(pgi_df, by = c("genetic_iid" = "IID")) %>%
    # Preserve the genotype-side IID (the join consumed pgi_df$IID into
    # genetic_iid). Needed by downstream pipelines that touch the per-cohort
    # bfiles.
    rename(IID = genetic_iid)
}

# ------------------------------------------------------------------
# Coverage report & twin diagnostic
# ------------------------------------------------------------------
report_merge <- function(df, name) {
  n  <- nrow(df)
  ne <- sum(!is.na(df$education))
  np <- sum(!is.na(df$PGI_Edu))
  nb <- sum(!is.na(df$PGI_Edu) & !is.na(df$education))
  message(sprintf("  %s: N=%d | edu=%d | PGI_Edu=%d | both=%d",
                  name, n, ne, np, nb))
}

twin_correlation_check <- function(df_twin) {
  # Defensive: ptyp may arrive as labelled / character / numeric depending on
  # how the upstream phenotype script saved it. Coerce to integer before the
  # %in% check so we don't silently filter out everything when ptyp is e.g.
  # haven_labelled.
  tc <- df_twin %>%
    mutate(ptyp = suppressWarnings(as.integer(ptyp))) %>%
    filter(ptyp %in% c(1L, 2L), !is.na(PGI_Edu)) %>%
    select(fid, ptyp, PGI_Edu) %>%
    pivot_wider(names_from = ptyp, values_from = PGI_Edu,
                names_prefix = "PGI_ptyp") %>%
    filter(!is.na(PGI_ptyp1) & !is.na(PGI_ptyp2))
  if (nrow(tc) > 10) {
    r <- cor(tc$PGI_ptyp1, tc$PGI_ptyp2, use = "complete.obs")
    message(sprintf("  Twin PGI_Edu correlation (n=%d pairs): r = %.3f",
                    nrow(tc), r))
    if (abs(r) < 0.2)
      warning("[TwinLife] Twin PGI correlation near 0 — crosswalk may be incorrect!")
    invisible(data.frame(n_pairs = nrow(tc), r = r))
  } else {
    message("  Too few twin pairs to validate (< 10).")
    invisible(data.frame(n_pairs = nrow(tc), r = NA_real_))
  }
}

# ------------------------------------------------------------------
# Common column contract (keeps columns in a stable order; missing
# ones drop silently via select(any_of(...)))
# ------------------------------------------------------------------
COMMON_COLS_PHENO <- c(
  "pid",
  "IID",                                  # genotype-side IID (preserved
                                          # through the per-cohort
                                          # crosswalks); needed by any
                                          # downstream step that touches the
                                          # cohort bfiles (GRM construction,
                                          # GCTA --keep, etc.).
  "cohort", "gender", "birth_year", "education",
  "mother_education", "father_education",
  "mother_birth_year", "father_birth_year", "parental_education",
  "east_west", "east_west_mig", "east_west_mig_edu",
  "height", "bmi", "height_age", "bmi_age",
  "fid", "ptyp"
)

COMMON_COLS_PGI <- c(
  # bare PGI_<trait> = analysis-friendly cross-cohort z (was *_global_z
  # in the input TSV); see normalize_pgi_columns() above
  "PGI_Edu", "PGI_Cog", "PGI_nonCog", "PGI_Height",
  "PGI_Edu_raw",      "PGI_Cog_raw",      "PGI_nonCog_raw",      "PGI_Height_raw",
  "PGI_Edu_within_z", "PGI_Cog_within_z", "PGI_nonCog_within_z", "PGI_Height_within_z"
)

COMMON_COLS_PCS <- c(paste0("PC", 1:20), paste0("crossPC", 1:20))

COMMON_COLS <- c(COMMON_COLS_PHENO, COMMON_COLS_PGI, COMMON_COLS_PCS, "geno_qc")

# ------------------------------------------------------------------
# Save the combined RDS with metadata
# ------------------------------------------------------------------
save_combined <- function(df, data_root,
                          date_str = format(Sys.Date(), "%Y%m%d"),
                          pgi_dir = NA_character_) {
  out_dir  <- file.path(data_root, "03_Merge")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir,
                        sprintf("Combined_Harmonized_%s.rds", date_str))
  attr(df, "pgi_source_dir") <- pgi_dir
  attr(df, "build_date")     <- date_str
  saveRDS(df, out_file)
  message(sprintf("  Saved: %s (%d rows x %d cols)",
                  out_file, nrow(df), ncol(df)))
  invisible(out_file)
}
