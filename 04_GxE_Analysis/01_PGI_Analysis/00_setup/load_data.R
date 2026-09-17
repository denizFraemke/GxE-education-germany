# load_data.R
# Locate and load the latest Combined_Harmonized_*.rds.
#
# After source()ing, the calling environment contains:
#   DATA_FILE   : absolute path of the loaded RDS
#   df_raw      : the loaded data frame (with TwinLife filter applied if requested)
#   gender_col  : "gender" or "Gender" depending on which is present
#
# Pre-conditions in the caller's frame:
#   DATA_DIR         (chr): directory to search for Combined_Harmonized_*.rds
#   EXCLUDE_TWINLIFE (lgl): if TRUE, drop cohort == "TwinLife" before returning
#   RESULTS          (list): collector for the Analysis_Metadata sheet

stopifnot(exists("DATA_DIR"),
          exists("EXCLUDE_TWINLIFE"),
          exists("RESULTS"))

# ---- 2a. Find and load latest combined RDS ----
rds_files <- list.files(DATA_DIR,
                         pattern = "^Combined_Harmonized_\\d{8}\\.rds$",
                         full.names = TRUE)
if (length(rds_files) == 0)
  stop("No combined RDS found in DATA_DIR. Run the merge script first.")
DATA_FILE <- rds_files[which.max(file.info(rds_files)$mtime)]
cat("Loading:", DATA_FILE, "\n")
df_raw <- readRDS(DATA_FILE)

if (EXCLUDE_TWINLIFE) {
  df_raw <- df_raw |> dplyr::filter(cohort != "TwinLife")
  cat("Run option: TwinLife excluded from this analysis.\n")
}

# normalize_pgi_columns() in merge_common.R promotes the cross-cohort
# _global_z value to the bare PGI_<trait> column, so the bare column is
# the analysis-friendly z-score.
required_pgi <- paste0("PGI_", c("Edu", "Cog", "nonCog", "Height"))
missing_pgi <- setdiff(required_pgi, names(df_raw))
if (length(missing_pgi) > 0) {
  stop("Combined RDS is missing PGI columns: ",
       paste(missing_pgi, collapse = ", "),
       "\nRun 03_Merge/merge.R against the Tardis-final PGIs.")
}
cat("PGI source dir:", attr(df_raw, "pgi_source_dir", exact = TRUE), "\n")

# ---- 2b. Diagnostics ----
cat("\n--- DIAGNOSTICS ---\n")
cat("Dimensions:", nrow(df_raw), "rows x", ncol(df_raw), "columns\n")
cat("Columns:", paste(sort(names(df_raw)), collapse = ", "), "\n")
cat("\nCohort breakdown:\n"); print(table(df_raw$cohort, useNA = "ifany"))

# Check required columns
required_cols <- c("pid", "cohort", "birth_year", "education", "east_west",
                   "PGI_Edu")
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Check for Gender column (capital G from merge script)
gender_col <- if ("gender" %in% names(df_raw)) { "gender"
} else if ("Gender" %in% names(df_raw)) { "Gender"
} else { stop("No gender/Gender column found.") }
cat("Gender column found as:", gender_col, "\n")

RESULTS[["Analysis_Metadata"]] <- data.frame(
  key   = c("data_file", "exclude_twinlife"),
  value = c(DATA_FILE, as.character(EXCLUDE_TWINLIFE)),
  stringsAsFactors = FALSE
)
