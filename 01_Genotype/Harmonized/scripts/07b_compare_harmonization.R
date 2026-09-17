#!/usr/bin/env Rscript
# =============================================================================
# 07b_compare_harmonization.R — Compare PGI harmonization approaches
# =============================================================================
# Compares PGI distributions and (if phenotype data available) association
# statistics across cohorts for different harmonization approaches.
#
# APPROACH A: UKB-aligned genotypes (from 03c_align_to_ukb.sh)
#   - Genotypes aligned to the UKB reference panel before scoring
#
# APPROACH B: Cross-cohort pooled standardization (from 07_assemble_output.R)
#   - Raw scores z-standardized using pooled mean/SD across all cohorts
#
# This script:
#   1. Prints distribution statistics for each approach
#   2. Computes cross-cohort heterogeneity (CV of cohort means)
#
# USAGE:
#   Rscript scripts/07b_compare_harmonization.R [--approach-a-dir PATH]
#     --approach-a-dir: Path to approach A final output directory (optional)
#   (typically via `bash run_pipeline.sh --step 07b`)
#
# INPUTS:
#   ${FINAL_OUT}/<COHORT>_PGI_PCs.tsv  — per-cohort final tables from step 07
#   (optional, with --approach-a-dir) a parallel `final/` tree produced by an
#   alternate run with step 03c (UKB alignment) enabled.
#
# OUTPUTS:
#   ${OUT_ROOT}/diagnostics/harmonization_comparison.tsv (+ related TSVs)
#   Console summary of cross-cohort heterogeneity statistics.
#
# DEPENDS ON: step 07 (assembled final TSVs). The UKB-aligned approach
# comparison additionally requires step 03c to have been run.
#
# NOTE: This script compares ONLY the cross-cohort standardization approach
#       (which is always available). The UKB-aligned approach comparison
#       requires running the UKB alignment step first.
# =============================================================================

suppressPackageStartupMessages({
    for (.pkg in c("data.table")) {
        if (!requireNamespace(.pkg, quietly = TRUE)) {
            lib_path <- Sys.getenv("R_LIBS_USER")
            if (!nzchar(lib_path)) lib_path <- file.path(Sys.getenv("HOME"), "R", "library")
            dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
            install.packages(.pkg, repos = "https://cloud.r-project.org",
                             lib = lib_path, quiet = TRUE)
            .libPaths(c(lib_path, .libPaths()))
        }
        library(.pkg, character.only = TRUE)
    }
    rm(.pkg)
})

cat("================================================================================\n")
cat("  Step 07b: Comparing PGI Harmonization Approaches\n")
cat("================================================================================\n\n")

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
default_root <- if (length(script_path) > 0) {
    normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
    getwd()
}

proj_root <- Sys.getenv("PROJ_ROOT", unset = default_root)
if (proj_root == "") proj_root <- default_root

# Parse arguments
approach_a_dir <- NULL
for (i in seq_along(args)) {
    if (args[i] == "--approach-a-dir" && i + 1 <= length(args)) {
        approach_a_dir <- args[i + 1]
    }
}

cohorts <- c("BASEII", "SHIP0", "SHIPTD", "SOEP", "TWINLIFE")
traits  <- c("EA4", "Cog", "NonCog", "Height")

# =============================================================================
# APPROACH B: Cross-cohort pooled standardization (always available)
# =============================================================================
cat("--- APPROACH B: Cross-Cohort Pooled Standardization ---\n\n")

approach_b_dir <- file.path(proj_root, "output", "final")

approach_b_stats <- list()
for (cohort in cohorts) {
    f <- file.path(approach_b_dir, paste0(cohort, "_PGI_PCs.tsv"))
    if (!file.exists(f)) {
        cat(sprintf("  WARN: %s not found — run 07_assemble_output.R first\n", basename(f)))
        next
    }
    d <- fread(f, nrows = 0)

    for (trait in traits) {
        raw_col <- paste0("PGI_", trait, "_raw")
        gz_col <- paste0("PGI_", trait, "_global_z")
        final_col <- paste0("PGI_", trait)

        if (!raw_col %in% names(d)) next

        select_cols <- unique(c(raw_col, gz_col, final_col))
        dd <- fread(f, select = intersect(select_cols, names(d)))
        for (missing_col in setdiff(select_cols, names(dd))) {
            dd[[missing_col]] <- NA_real_
        }
        approach_b_stats[[paste(cohort, trait)]] <- data.frame(
            cohort = cohort,
            trait = trait,
            approach = "B (cross-cohort pooled)",
            n = sum(!is.na(dd[[raw_col]])),
            mean_raw = mean(dd[[raw_col]], na.rm = TRUE),
            sd_raw = sd(dd[[raw_col]], na.rm = TRUE),
            mean_global_z = mean(dd[[gz_col]], na.rm = TRUE),
            sd_global_z = sd(dd[[gz_col]], na.rm = TRUE),
            mean_final = mean(dd[[final_col]], na.rm = TRUE),
            sd_final = sd(dd[[final_col]], na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }
}

if (length(approach_b_stats) > 0) {
    approach_b_df <- rbindlist(approach_b_stats)
    cat("Approach B — Distribution Statistics:\n\n")
    print(approach_b_df[, .(cohort, trait, n, mean_raw, sd_raw, mean_global_z, sd_global_z)])
}

# =============================================================================
# APPROACH A: UKB-aligned (if available)
# =============================================================================
approach_a_df <- NULL
if (!is.null(approach_a_dir) && dir.exists(approach_a_dir)) {
    cat("\n--- APPROACH A: UKB Reference Alignment ---\n\n")
    approach_a_stats <- list()
    for (cohort in cohorts) {
        f <- file.path(approach_a_dir, paste0(cohort, "_PGI_PCs.tsv"))
        if (!file.exists(f)) next

        d <- fread(f, nrows = 0)
        for (trait in traits) {
            raw_col <- paste0("PGI_", trait, "_raw")
            gz_col <- paste0("PGI_", trait, "_global_z")
            final_col <- paste0("PGI_", trait)

            if (!raw_col %in% names(d)) next

            select_cols <- unique(c(raw_col, gz_col, final_col))
            dd <- fread(f, select = intersect(select_cols, names(d)))
            for (missing_col in setdiff(select_cols, names(dd))) {
                dd[[missing_col]] <- NA_real_
            }
            approach_a_stats[[paste(cohort, trait)]] <- data.frame(
                cohort = cohort,
                trait = trait,
                approach = "A (UKB-aligned)",
                n = sum(!is.na(dd[[raw_col]])),
                mean_raw = mean(dd[[raw_col]], na.rm = TRUE),
                sd_raw = sd(dd[[raw_col]], na.rm = TRUE),
                mean_global_z = mean(dd[[gz_col]], na.rm = TRUE),
                sd_global_z = sd(dd[[gz_col]], na.rm = TRUE),
                mean_final = mean(dd[[final_col]], na.rm = TRUE),
                sd_final = sd(dd[[final_col]], na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(approach_a_stats) > 0) {
        approach_a_df <- rbindlist(approach_a_stats)
        cat("Approach A — Distribution Statistics:\n\n")
        print(approach_a_df[, .(cohort, trait, n, mean_raw, sd_raw, mean_global_z, sd_global_z)])
    }
} else {
    cat("\n--- APPROACH A: UKB Reference Alignment ---\n\n")
    cat("  NOTE: Approach A results not available.\n")
    cat("  To compare with UKB alignment, run:\n")
    cat("    bash scripts/03c_align_to_ukb.sh\n")
    cat("    Rscript scripts/07_assemble_output.R --approach-a\n")
    cat("  Then re-run with: Rscript scripts/07b_compare_harmonization.R --approach-a-dir output/final_ukb\n\n")
}

# =============================================================================
# COMPARE: Raw vs Global_z vs Final
# =============================================================================
cat("\n================================================================================\n")
cat("  COMPARISON: Raw vs Global_z vs Final (Approach B)\n")
cat("================================================================================\n\n")

if (length(approach_b_stats) > 0) {
    # Compare SDs across cohorts for each standardization stage
    # Lower CV = more comparable across cohorts
    cat("Coefficient of Variation (CV = SD/mean) of raw PGIs across cohorts:\n")
    for (trait_name in traits) {
        subset <- approach_b_df[trait == trait_name, ]
        raw_cvs <- subset[, sd_raw / abs(mean_raw)]
        overall_cv <- sd(subset$mean_raw) / abs(mean(subset$mean_raw))
        cat(sprintf("  %s: CV of cohort means = %.3f (SDs range: %.4f - %.4f)\n",
                    trait_name, overall_cv,
                    min(subset$sd_raw), max(subset$sd_raw)))
    }

    cat("\nCoefficient of Variation of global_z PGIs across cohorts:\n")
    for (trait_name in traits) {
        subset <- approach_b_df[trait == trait_name, ]
        # After global z-standardization, means should all be ~0
        overall_cv <- sd(subset$mean_global_z) / abs(mean(subset$mean_global_z))
        cat(sprintf("  %s: CV of cohort means = %.6f (SDs range: %.4f - %.4f)\n",
                    trait_name, overall_cv,
                    min(subset$sd_global_z), max(subset$sd_global_z)))
    }

    cat("\nCoefficient of Variation of final (canonical = raw) PGIs across cohorts:\n")
    for (trait_name in traits) {
        subset <- approach_b_df[trait == trait_name, ]
        overall_cv <- sd(subset$mean_final) / abs(mean(subset$mean_final))
        cat(sprintf("  %s: CV of cohort means = %.6f\n",
                    trait_name, overall_cv))
    }

    # NonCog sign check
    cat("\n*** NonCog Sign Consistency Check ***\n")
    ncog_subset <- approach_b_df[trait == "NonCog", ]
    signs <- sign(ncog_subset$mean_global_z)
    if (all(is.na(signs))) {
        cat("  WARN: NonCog_global_z is unavailable in the current output files\n")
    } else if (all(signs > 0, na.rm = TRUE)) {
        cat("  OK: All cohorts have positive NonCog_global_z\n")
    } else if (all(signs < 0, na.rm = TRUE)) {
        cat("  OK: All cohorts have negative NonCog_global_z\n")
    } else {
        cat("  FAIL: NonCog has MIXED signs across cohorts!\n")
        for (i in 1:nrow(ncog_subset)) {
            cat(sprintf("    %s: mean=%.4f (%s)\n",
                        ncog_subset$cohort[i],
                        ncog_subset$mean_global_z[i],
                        if (ncog_subset$mean_global_z[i] > 0) "POSITIVE" else "NEGATIVE"))
        }
        cat("  -> Run 03b_fix_strand.sh to fix strand orientation issues\n")
    }
}

# =============================================================================
# COLUMN USAGE
# =============================================================================
cat("\n================================================================================\n")
cat("  COLUMN USAGE\n")
cat("================================================================================\n\n")

cat("Approach B (cross-cohort pooled standardization) is the primary pipeline:\n")
cat("  - PGI_<trait>_global_z: cross-cohort z, pooled SD units\n")
cat("  - PGI_<trait>: canonical column, equal to PGI_<trait>_raw\n")
cat("  - Downstream merge residualizes on PCs and re-derives the analysis z\n")
cat("  - Verify that NonCog has consistent sign across cohorts before analysis\n\n")

cat("================================================================================\n")
