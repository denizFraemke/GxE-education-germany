#!/usr/bin/env Rscript
# =============================================================================
# export_result_tables.R — dump every AGGREGATE result table from the model
# caches to tidy CSVs under manuscript_export/, in a human-browsable tree. This
# is the single interface to the manuscript repo: analysis produces tidy CSVs,
# the manuscript only reads them (no xlsx, no xlsx->CSV round-trip, no readxl).
#
# Layout produced under OUT_DIR:
#   descriptives/                    sample descriptives, PGI correlations, attrition
#   attainment/{main,sensitivities,appendix}/   Section A tables (A1..A6, A0g, ...)
#   mobility/{main,sensitivities,appendix}/      Section B tables (B0..B6, B0g, ...)
#   README.md                        auto-generated legend (codes + file types)
#
# SAFETY: only allow-listed aggregate tables are written. Individual-level data
# frames (dat_cluster*, dat_edu*, dat_mob*, df_*) are NEVER exported.
#
# Usage:
#   RUN_DIR=/…/runs/RUN_<TS> OUT_DIR=/…/manuscript_export Rscript export_result_tables.R
#   CACHES="a.rds,b.rds"     OUT_DIR=/tmp/proto            Rscript export_result_tables.R
# =============================================================================

OUT_DIR <- Sys.getenv("OUT_DIR", ""); if (!nzchar(OUT_DIR)) stop("Set OUT_DIR")

caches <- Sys.getenv("CACHES", "")
if (nzchar(caches)) {
  cache_files <- strsplit(caches, ",")[[1]]
} else {
  RUN_DIR <- Sys.getenv("RUN_DIR", ""); if (!nzchar(RUN_DIR)) stop("Set RUN_DIR or CACHES")
  cache_files <- Filter(file.exists, c(
    file.path(RUN_DIR, "A_Education",     "A_Education_Models.rds"),
    file.path(RUN_DIR, "A_Education",     "A_Sensitivities_Models.rds"),
    file.path(RUN_DIR, "A_Education",     "A_Appendix_Models.rds"),
    file.path(RUN_DIR, "B_Mobility",      "B_Mobility_Models.rds"),
    file.path(RUN_DIR, "B_Mobility",      "B_Sensitivities_Models.rds"),
    file.path(RUN_DIR, "B_Mobility",      "B_Appendix_Models.rds"),
    file.path(RUN_DIR, "00_Descriptives", "Descriptives_Tables.rds"),
    # Optional add-on caches. A full run folds these tables into the
    # *_Models.rds / *_Sensitivities_Models.rds above; they are read separately
    # only when a run directory carries them as standalone files.
    file.path(RUN_DIR, "A_Education",     "A_Education_GenderLOCO.rds"),
    file.path(RUN_DIR, "A_Education",     "A_Sensitivities_AltPGI.rds"),
    file.path(RUN_DIR, "B_Mobility",      "B_Mobility_GenderLOCO.rds"),
    file.path(RUN_DIR, "B_Mobility",      "B_Sensitivities_AltPGI.rds")))
}

# Individual-level guard: skip a table by NAME (dat_/df_) OR by CONTENT (it carries
# an individual-ID column). Every other data.frame/matrix is an aggregate result
# table and is exported — content-based, so it also captures the bespoke
# sensitivities/appendix tables (S1_AltPGIs, Method_Comparison, F_linear_edu, …)
# that a name allow-list would miss.
ID_COLS <- c("IID", "PID", "FID", "GENETIC_ID", "SAMPLE_ID")
is_individual <- function(nm, obj)
  grepl("^dat_|^dat$|^df_|^df$", nm) ||
  any(toupper(as.character(colnames(obj))) %in% ID_COLS)

# cache basename -> (section, subrun)
classify_cache <- function(path) {
  b <- basename(path)
  section <- if (grepl("^A", b)) "attainment" else if (grepl("^B", b)) "mobility" else "descriptives"
  subrun  <- if (grepl("Sensitiv", b)) "sensitivities" else if (grepl("Append", b)) "appendix" else "main"
  list(section = section, subrun = subrun)
}
# object -> relative folder
route <- function(section, subrun, nm) {
  if (grepl("^desc_|^pgi_correlations$|^attrition", nm)) return("descriptives")
  if (section == "descriptives") return("descriptives")
  file.path(section, subrun)
}

written <- 0; skipped_priv <- character(0); index <- list()
for (f in cache_files) {
  m <- readRDS(f); cc <- classify_cache(f)
  for (nm in names(m)) {
    obj <- m[[nm]]
    if (!(is.data.frame(obj) || is.matrix(obj))) next
    if (is_individual(nm, obj)) { skipped_priv <- c(skipped_priv, nm); next }
    rel <- route(cc$section, cc$subrun, nm)
    dir.create(file.path(OUT_DIR, rel), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(as.data.frame(obj), file.path(OUT_DIR, rel, paste0(nm, ".csv")), row.names = FALSE)
    index[[length(index) + 1]] <- data.frame(folder = rel, file = paste0(nm, ".csv"),
                                             rows = nrow(obj), cols = ncol(obj))
    written <- written + 1
  }
}

# ---- README legend --------------------------------------------------------
codes <- c(
  "A1 = PGI x birth-year x region (attainment)",
  "A2 = PGI x post-reunification x region (attainment)",
  "A3 = birth-cohort nonlinearity, GAM (attainment)",
  "A4 = four-way PGI x birth-year x region x gender (attainment)",
  "A5 = per-gender cohort-trajectory GAM (attainment)",
  "A6 = region x gender varying-coefficient GAM (attainment)",
  "A0g = overall region x gender PGI contrast (attainment)",
  "B0 = East-West PGI-mobility contrast; B0_pre = pre-reunification",
  "B1 = PGI x birth-year x region (mobility)",
  "B2 = PGI x post-reunification x region (mobility)",
  "B3 = birth-cohort nonlinearity, GAM (mobility)",
  "B4/B5/B6/B0g = mobility gender analogues of A4/A5/A6/A0g",
  "A4_LOO_gender / B4_LOO_gender = leave-one-cohort-out of the POOLED PGI x gender term",
  "A5_LOO / B5_LOO = leave-one-cohort-out of the per-gender GAM (omnibus + curvature)",
  "B5_LOO_diffsmooth / _diff_intervals = per-fold female-minus-male difference smooth + excludes-zero intervals",
  "S1_finding1 = alt-PGI attainment across cohorts (pooled PGI x BY, RE, common sample)",
  "S1_finding2 = alt-PGI East-West mobility contrast (PGI x region, RE, common sample)",
  "S1_finding3_att / _mob = alt-PGI pooled PGI x gender (RE, common sample)",
  "S1_A5_nonlin / S1_B5_nonlin = alt-PGI per-gender GAM (omnibus + curvature; B5 adds diff smooth)")
types <- c(
  "*_Coefficients  = model coefficient table (estimate, SE, t/z, p, 95% CI)",
  "*_ModelFit      = model fit statistics",
  "*_NestedLRT     = linear-vs-nonlinear GAM likelihood-ratio test",
  "*_cell_LRT      = per-cell nonlinearity tests",
  "*_Smooths / *_curves / *_display_curves = fitted GAM smooth predictions",
  "*_RegionSlopes / *_OverallSlopes = per-region / overall PGI slopes",
  "*_LOO / *_LOO_* = leave-one-cohort-out refits",
  "*_K_Index       = GAM basis-dimension check",
  "*_Diagnostic    = specification diagnostics")
idx <- if (length(index)) do.call(rbind, index) else data.frame()
readme <- c(
  "# manuscript_export — tidy result tables", "",
  "Machine-readable CSVs generated directly from the analysis model caches",
  "(`export_result_tables.R` + `export_manuscript_aggregates.R`). The manuscript",
  "repo reads only this folder. No xlsx and no round-trip conversion.", "",
  "## Folders", "",
  "- `descriptives/` — sample descriptives, PGI correlations, attrition",
  "- `attainment/{main,sensitivities,appendix}/` — Section A (educational attainment)",
  "- `mobility/{main,sensitivities,appendix}/` — Section B (educational mobility)",
  "- `fig*.csv` — figure-source data at the export root (from export_manuscript_aggregates.R)", "",
  "## Analysis codes", "", paste0("- ", codes), "",
  "## File-type suffixes", "", paste0("- ", types), "",
  "## Contents", "",
  if (nrow(idx)) c("| folder | file | rows | cols |", "|---|---|---|---|",
                   sprintf("| %s | %s | %d | %d |", idx$folder, idx$file, idx$rows, idx$cols)) else "(none)")
writeLines(readme, file.path(OUT_DIR, "README.md"))

cat(sprintf("Wrote %d aggregate result CSVs + README.md to %s\n", written, OUT_DIR))
if (length(skipped_priv))
  cat("Skipped individual-level frames (never exported): ",
      paste(unique(skipped_priv), collapse = ", "), "\n", sep = "")
