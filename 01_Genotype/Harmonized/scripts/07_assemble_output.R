#!/usr/bin/env Rscript
# =============================================================================
# 07_assemble_output.R — Assemble final PGI + PC tables with
# CROSS-COHORT HARMONIZATION
# =============================================================================
# For each cohort, produces a table with:
#   - IID (individual ID)
#   - PGI_EA4, PGI_Cog, PGI_NonCog, PGI_Height
#   - PC1 through PC20 (standard ancestry PCs)
#
# HARMONIZATION STRATEGY (cross-cohort):
#
#   All five cohorts are scored against a single LOO weight file per trait
#   (see config.sh: get_ea4_trait_for_cohort, get_noncog_trait_for_cohort),
#   so raw PGIs are already on a common scale and directly comparable across
#   cohorts. No within-cohort standardization is applied — that would erase
#   real cohort-level differences in genetic propensity that the analysis
#   wants to preserve. We do also write a cross-cohort z-standardized
#   `_global_z` column for downstream code that wants standardized inputs.
#
# OUTPUT COLUMNS per trait:
#   - PGI_{trait}_raw        : Raw PLINK2 BETA_SUM (preserved for transparency)
#   - PGI_{trait}_global_z   : (raw - pooled_mean) / pooled_sd across all cohorts
#   - PGI_{trait}            : Canonical PGI = PGI_{trait}_raw. This is the
#                              column you'd use directly in regressions and
#                              cross-cohort comparisons.
#
# USAGE:
#   Rscript scripts/07_assemble_output.R
#   (typically via `bash run_pipeline.sh --step 07`)
#
# INPUTS:
#   ${OUT_ROOT}/scores/<COHORT>_<trait>.sscore        — main PGI scores
#   ${OUT_ROOT}/scores/<COHORT>_<trait>_r2strict.sscore (when present)
#   ${OUT_ROOT}/pcs/<COHORT>_pcs.eigenvec             — within-cohort PCs
#   ${OUT_ROOT}/pcs/cross_cohort_pcs.eigenvec         — cross-cohort PCs
#
# OUTPUTS:
#   ${FINAL_OUT}/<COHORT>_PGI_PCs.tsv                 — per-cohort deliverable
#   ${FINAL_OUT}/PGI_summary_statistics.tsv           — descriptives
#   ${OUT_ROOT}/diagnostics/global_pgi_statistics.tsv — pooled means/SDs
#
# DEPENDS ON: steps 04, 04b, 05 (always); 06b when its sscores exist.
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
cat("  Step 07: Assembling final PGI + PC tables\n")
cat("  Cross-Cohort Harmonization\n")
cat("================================================================================\n\n")

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && a != "") a else b

# --- Configuration -----------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
default_root <- if (length(script_path) > 0) {
    normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
    getwd()
}

proj_root <- Sys.getenv("PROJ_ROOT", unset = default_root)
if (proj_root == "") proj_root <- default_root

scores_dir   <- file.path(proj_root, "output", "scores")
pcs_dir      <- file.path(proj_root, "output", "pcs")
final_dir    <- file.path(proj_root, "output", "final")
log_dir      <- file.path(proj_root, "logs")
decision_log <- file.path(proj_root, "decision_log.md")
diag_dir     <- file.path(proj_root, "output", "diagnostics")

dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

cohorts <- c("BASEII", "SHIP0", "SHIPTD", "SOEP", "TWINLIFE")
traits  <- c("EA4", "Cog", "NonCog", "Height")

# --- Helper: read PLINK2 .sscore file ----------------------------------------
read_sscore <- function(path) {
    if (!file.exists(path)) {
        warning("Score file not found: ", path)
        return(NULL)
    }
    dt <- fread(path, header = TRUE, nThread = 4,
                colClasses = list(character = c("#FID", "IID")))
    if ("#FID" %in% names(dt)) setnames(dt, "#FID", "FID")
    dt
}

pick_score_col <- function(dt) {
    preferred <- c("BETA_SUM", "SCORE1_SUM", "BETA_AVG", "SCORE1_AVG")
    hit <- preferred[preferred %in% names(dt)]
    if (length(hit) > 0) return(hit[1])

    fallback <- grep("SUM$|AVG$", names(dt), value = TRUE)
    fallback <- fallback[!grepl("ALLELE|DOSAGE|DENOM", fallback, ignore.case = TRUE)]
    if (length(fallback) > 0) return(fallback[1])

    score_like <- grep("SCORE|BETA", names(dt), value = TRUE)
    score_like <- score_like[!grepl("ALLELE|DOSAGE|DENOM", score_like, ignore.case = TRUE)]
    if (length(score_like) > 0) return(score_like[1])

    NULL
}

# --- Helper: z-standardize ---------------------------------------------------
z_scale <- function(x) {
    x <- as.numeric(x)
    mu <- mean(x, na.rm = TRUE)
    s  <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mu) / s
}

# =============================================================================
# PASS 1: Read all raw PGI scores from all cohorts
# =============================================================================
cat("--- PASS 1: Reading raw PGI scores from all cohorts ---\n\n")

all_raw <- list()

for (coh in cohorts) {
    cat(sprintf("  Reading %s...\n", coh))
    pgi_list <- list()

    for (trait in traits) {
        sscore_file <- file.path(scores_dir, paste0(coh, "_", trait, ".sscore"))
        dt <- read_sscore(sscore_file)
        if (is.null(dt)) {
            cat(sprintf("    WARN: %s %s score file not found\n", coh, trait))
            next
        }
        score_col <- pick_score_col(dt)
        if (is.null(score_col)) {
            cat(sprintf("    WARN: No score column found in %s\n", sscore_file))
            next
        }
        pgi_col <- paste0("PGI_", trait)
        pgi_dt <- dt[, .(IID = as.character(IID),
                         cohort = coh,
                         score = as.numeric(get(score_col)))]
        setnames(pgi_dt, "score", pgi_col)
        pgi_list[[trait]] <- pgi_dt
        cat(sprintf("    %s: %d individuals, score col = %s\n", trait, nrow(pgi_dt), score_col))
    }

    if (length(pgi_list) > 0) {
        merged <- Reduce(function(a, b) merge(a, b, by = c("IID", "cohort"), all = TRUE), pgi_list)
        all_raw[[coh]] <- merged
    }
}

if (length(all_raw) == 0) {
    stop("No PGI data found for any cohort!")
}

raw_pooled <- rbindlist(all_raw, fill = TRUE)
cat(sprintf("\n  Total individuals across all cohorts: %d\n", nrow(raw_pooled)))

# =============================================================================
# PASS 2: Compute GLOBAL mean and SD for each PGI trait (cross-cohort)
# =============================================================================
cat("\n--- PASS 2: Computing global mean and SD (cross-cohort) ---\n\n")

pgi_cols <- grep("^PGI_", names(raw_pooled), value = TRUE)
pgi_cols <- setdiff(pgi_cols, "cohort")

global_stats <- list()
for (col in pgi_cols) {
    vals <- raw_pooled[[col]]
    n_valid <- sum(!is.na(vals))
    mu <- mean(vals, na.rm = TRUE)
    s <- sd(vals, na.rm = TRUE)
    global_stats[[col]] <- list(mean = mu, sd = s, n = n_valid)
    cat(sprintf("  %s: global_mean = %.6f, global_sd = %.6f, n = %d\n",
                col, mu, s, n_valid))
}

# Save global stats for downstream use
global_stats_df <- rbindlist(lapply(names(global_stats), function(col) {
    gs <- global_stats[[col]]
    data.frame(
        trait = sub("^PGI_", "", col),
        global_mean = gs$mean,
        global_sd = gs$sd,
        n = gs$n,
        stringsAsFactors = FALSE
    )
}))
fwrite(global_stats_df,
       file.path(diag_dir, "global_pgi_statistics.tsv"),
       sep = "\t")

# =============================================================================
# PASS 3: Assemble per-cohort tables with harmonized scores
# =============================================================================
cat("\n--- PASS 3: Assembling per-cohort tables with cross-cohort harmonization ---\n\n")

all_summary <- list()

for (coh in cohorts) {
    cat(sprintf("\n--- Processing cohort: %s ---\n", coh))

    # NOTE: Use coh (not "cohort") as loop variable to avoid data.table
    # column-name shadowing. In data.table, bare `cohort` inside [i, ]
    # resolves to the COLUMN "cohort", not the loop variable.
    cohort_mask <- raw_pooled[["cohort"]] == coh
    cohort_data <- raw_pooled[cohort_mask, ]
    if (nrow(cohort_data) == 0) {
        cat(sprintf("  SKIP: No data for %s\n", coh))
        next
    }
    cat(sprintf("  Individuals: %d\n", nrow(cohort_data)))

    # --- Apply cross-cohort harmonization ------------------------------------
    for (col in pgi_cols) {
        raw_col <- paste0(col, "_raw")
        global_z_col <- paste0(col, "_global_z")
        final_col <- col

        raw_vals <- cohort_data[[col]]
        mu <- global_stats[[col]]$mean
        s <- global_stats[[col]]$sd

        # Preserve raw scores
        cohort_data[, (raw_col) := raw_vals]

        # STEP 1: Cross-cohort harmonized (global z-standardization)
        if (!is.na(s) && s > 0) {
            cohort_data[, (global_z_col) := (raw_vals - mu) / s]
        } else {
            cohort_data[, (global_z_col) := NA_real_]
        }

        # Canonical PGI_{trait} = PGI_{trait}_raw. Single weight file per
        # trait across all cohorts means raw PGIs are directly comparable;
        # no within-cohort standardization (that would erase real cohort-
        # level genetic-propensity differences).
        gz_vals <- cohort_data[[global_z_col]]
        cohort_data[, (final_col) := raw_vals]

        # Report statistics
        cat(sprintf("  %s:\n", col))
        cat(sprintf("    Raw:     mean=%.4f  sd=%.4f\n",
                    mean(raw_vals, na.rm = TRUE), sd(raw_vals, na.rm = TRUE)))
        cat(sprintf("    global_z: mean=%.4f  sd=%.4f\n",
                    mean(gz_vals, na.rm = TRUE), sd(gz_vals, na.rm = TRUE)))
        cat(sprintf("    Canonical (= raw): mean=%.4f  sd=%.4f\n",
                    mean(cohort_data[[final_col]], na.rm = TRUE),
                    sd(cohort_data[[final_col]], na.rm = TRUE)))
    }

    # --- Read within-cohort PCs ------------------------------------------------
    pc_file <- file.path(pcs_dir, paste0(coh, "_pcs.eigenvec"))
    if (file.exists(pc_file)) {
        pc_dt <- fread(pc_file, header = TRUE, nThread = 4,
                       colClasses = list(character = c("#FID", "IID")))
        if ("#FID" %in% names(pc_dt)) setnames(pc_dt, "#FID", "FID")

        pc_cols_orig <- grep("^PC", names(pc_dt), value = TRUE)
        if (length(pc_cols_orig) > 0) {
            new_pc_names <- paste0("PC", seq_along(pc_cols_orig))
            setnames(pc_dt, pc_cols_orig, new_pc_names)
            pc_dt <- pc_dt[, c("IID", new_pc_names), with = FALSE]
        }
        cohort_data <- merge(cohort_data, pc_dt, by = "IID", all.x = TRUE)
        n_with_pcs <- sum(!is.na(cohort_data$PC1))
        cat(sprintf("  Within-cohort PCs: %d / %d individuals matched\n",
                    n_with_pcs, nrow(cohort_data)))
    } else {
        cat(sprintf("  WARN: Within-cohort PC file not found: %s\n", pc_file))
    }

    # --- Read cross-cohort PCs ------------------------------------------------
    cross_pc_file <- file.path(pcs_dir, "cross_cohort_pcs.eigenvec")
    if (file.exists(cross_pc_file)) {
        cross_pc_dt <- fread(cross_pc_file, header = TRUE, nThread = 4,
                             colClasses = list(character = c("#FID", "IID")))
        if ("#FID" %in% names(cross_pc_dt)) setnames(cross_pc_dt, "#FID", "FID")

        cross_cols_orig <- grep("^PC", names(cross_pc_dt), value = TRUE)
        if (length(cross_cols_orig) > 0) {
            cross_pc_names <- paste0("crossPC", seq_along(cross_cols_orig))
            setnames(cross_pc_dt, cross_cols_orig, cross_pc_names)
            cross_pc_dt <- cross_pc_dt[, c("IID", cross_pc_names), with = FALSE]
        }
        cohort_data <- merge(cohort_data, cross_pc_dt, by = "IID", all.x = TRUE)
        n_with_cross_pcs <- sum(!is.na(cohort_data$crossPC1))
        cat(sprintf("  Cross-cohort PCs: %d / %d individuals matched\n",
                    n_with_cross_pcs, nrow(cohort_data)))
    } else {
        cat(sprintf("  NOTE: Cross-cohort PC file not found (run step 04b to generate)\n"))
    }

    # Remove cohort column from final output
    cohort_data[, cohort := NULL]

    # --- Summary statistics ---------------------------------------------------
    for (col in pgi_cols) {
        raw_col <- paste0(col, "_raw")
        global_z_col <- paste0(col, "_global_z")
        final_col <- col

        raw_vals <- cohort_data[[raw_col]]
        gz_vals <- cohort_data[[global_z_col]]
        final_vals <- cohort_data[[final_col]]

        all_summary <- c(all_summary, list(data.frame(
            cohort = coh,
            trait = sub("^PGI_", "", col),
            n = sum(!is.na(final_vals)),
            mean_raw = mean(raw_vals, na.rm = TRUE),
            sd_raw = sd(raw_vals, na.rm = TRUE),
            mean_global_z = mean(gz_vals, na.rm = TRUE),
            sd_global_z = sd(gz_vals, na.rm = TRUE),
            mean_final = mean(final_vals, na.rm = TRUE),
            sd_final = sd(final_vals, na.rm = TRUE),
            min_final = min(final_vals, na.rm = TRUE),
            max_final = max(final_vals, na.rm = TRUE),
            stringsAsFactors = FALSE
        )))
    }

    # --- Save per-cohort table -----------------------------------------------
    out_file <- file.path(final_dir, paste0(coh, "_PGI_PCs.tsv"))
    fwrite(cohort_data, out_file, sep = "\t")
    cat(sprintf("  Saved: %s (%d rows x %d cols)\n", out_file, nrow(cohort_data), ncol(cohort_data)))
}

# =============================================================================
# PASS 4: Save summary tables
# =============================================================================
if (length(all_summary) > 0) {
    summary_dt <- rbindlist(all_summary)
    summary_file <- file.path(final_dir, "PGI_summary_statistics.tsv")
    fwrite(summary_dt, summary_file, sep = "\t")
    cat(sprintf("\nPGI summary saved: %s\n", summary_file))
}

# =============================================================================
# PASS 5: Print cross-cohort comparison
# =============================================================================
cat("\n================================================================================\n")
cat("  CROSS-COHORT HARMONIZATION VERIFICATION\n")
cat("================================================================================\n\n")

cat("Global (pooled across all cohorts) PGI statistics:\n")
for (col in pgi_cols) {
    gs <- global_stats[[col]]
    cat(sprintf("  %s: mean=%.6f  sd=%.6f  n=%d\n", col, gs$mean, gs$sd, gs$n))
}

cat("\nPer-cohort final statistics (after harmonization):\n")
if (length(all_summary) > 0) {
    print(summary_dt[, .(cohort, trait, n, mean_final, sd_final, min_final, max_final)])
}

cat("\n*** KEY CHECKS ***\n")
if (length(all_summary) > 0) {
    # Check: cross-cohort raw-mean spread per trait. With a unified weight
    # file per trait, large between-cohort drifts in the raw mean would flag
    # something genuinely off (population structure, scoring code, etc.).
    # Small drifts (<= 1 SD on the pooled scale) are expected from real
    # cohort-level allele-frequency differences and are not a problem.
    for (col in pgi_cols) {
        trait_name <- sub("^PGI_", "", col)
        cohort_raw_means <- summary_dt[trait == trait_name, mean_raw]
        cohort_raw_sds   <- summary_dt[trait == trait_name, sd_raw]

        spread <- max(cohort_raw_means, na.rm = TRUE) -
                  min(cohort_raw_means, na.rm = TRUE)
        median_sd <- stats::median(cohort_raw_sds, na.rm = TRUE)
        spread_in_sd <- if (is.finite(median_sd) && median_sd > 0)
                            spread / median_sd else NA_real_

        if (!is.na(spread_in_sd) && spread_in_sd > 1.5) {
            cat(sprintf("  WARN: %s cross-cohort raw-mean spread is %.2f units (~%.2f cohort SDs)\n",
                        col, spread, spread_in_sd))
        } else {
            cat(sprintf("  OK: %s cross-cohort raw-mean spread = %.3f units (~%.2f cohort SDs)\n",
                        col, spread, spread_in_sd))
        }
    }

    # Check: NonCog sign consistency on the cross-cohort z scale. Opposite
    # signs would mean some cohorts' NonCog PGI moves the OPPOSITE direction
    # from others, which is the textbook strand-pathology signal.
    ncog_means <- summary_dt[trait == "NonCog", mean_global_z]
    ncog_signs <- sign(ncog_means)
    if (length(unique(ncog_signs)) > 1) {
        cat(sprintf("\n  CRITICAL: NonCog_global_z has OPPOSITE signs across cohorts!\n"))
        cat(sprintf("  Cohorts with positive NonCog: %s\n",
                    paste(summary_dt[trait == "NonCog" & mean_global_z > 0, cohort], collapse=", ")))
        cat(sprintf("  Cohorts with negative NonCog: %s\n",
                    paste(summary_dt[trait == "NonCog" & mean_global_z < 0, cohort], collapse=", ")))
        cat("  -> investigate strand handling before trusting NonCog cross-cohort comparisons.\n")
    } else {
        cat("\n  OK: NonCog has consistent sign across all cohorts.\n")
    }
}

# =============================================================================
# LOG DECISIONS
# =============================================================================
sink(decision_log, append = TRUE)
cat("\n## Step 07: Output Assembly — Cross-Cohort Harmonized (",
    format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n\n")
cat("- HARMONIZATION STRATEGY:\n")
cat("  1. Raw scores pooled across all cohorts to compute global mean & SD\n")
cat("  2. Cross-cohort z: (raw - global_mean) / global_sd -> PGI_trait_global_z\n")
cat("  3. Canonical PGI_trait := PGI_trait_raw (raw PLINK2 BETA_SUM).\n")
cat("     A single LOO weight file per trait makes raw PGIs comparable across\n")
cat("     cohorts; no within-cohort standardization is applied — that would\n")
cat("     erase real cohort-level genetic-propensity differences.\n")
cat("- OUTPUT COLUMNS:\n")
cat("  - PGI_trait_raw:      Raw PGI from PLINK2 (BETA_SUM)\n")
cat("  - PGI_trait_global_z: Cross-cohort z (pooled SD units)\n")
cat("  - PGI_trait:          Canonical = PGI_trait_raw (bare column is raw)\n")
cat("\n### Global statistics used for harmonization\n")
for (col in pgi_cols) {
    gs <- global_stats[[col]]
    cat(sprintf("  %s: mean=%.6f, sd=%.6f, n=%d\n", col, gs$mean, gs$sd, gs$n))
}
cat("\n")
sink()

cat("\n=== Step 07 complete (Cross-Cohort Harmonized) ===\n")
