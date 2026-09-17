#!/usr/bin/env Rscript
########################################################################
## 80 Years of GxE on Education in Germany — Merge validation report
##
## Reads the most recent Combined_Harmonized_<YYYYMMDD>.rds under
## DATA_ROOT/03_Merge/ and writes:
##   DATA_ROOT/03_Merge/validation_reports/validate_Harmonized_<date>.pdf
##   DATA_ROOT/03_Merge/validation_reports/validate_Harmonized_<date>_summary.tsv
##
## Pages:
##   1  Summary / cohort breakdown
##   2..5  Histograms of each PGI (raw) per cohort
##   6..7  Histograms of each outcome (education, height) per cohort
##   8  Scatterplots: PGI × education, faceted by trait, coloured by cohort
##   9  Scatterplots: PGI × height,    faceted by trait, coloured by cohort
##   9b Incremental R² (PGI added to covariates + PCs) per cohort × trait
##   10 PGI inter-correlation heatmaps per cohort
##   11 TwinLife parental-education reliability check
##   12 ID prefix uniqueness + duplicate-pid check
##   13 Ancestry PC scree — share of variance per PC, per cohort (scale-invariant)
##
## Usage:
##   Rscript 03_Merge/validate.R
########################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

# Defined up-front so it can be used anywhere below.
# (Base R provides %||% as of 4.4, but defining locally keeps the script
# portable to older R versions.)
`%||%` <- function(a, b) if (!is.null(a)) a else b

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
})
# DATA_ROOT: single source of truth in R/paths.R (override via DATA_ROOT env var).
source(file.path(SCRIPT_DIR, "R", "paths.R"), local = TRUE)
DATA_ROOT <- data_root()  # resolved once; override via DATA_ROOT env var

MERGE_DIR  <- file.path(DATA_ROOT, "03_Merge")
REPORT_DIR <- file.path(MERGE_DIR, "validation_reports")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Locate latest combined RDS --------------------------------------------
rds_pattern <- "^Combined_Harmonized_\\d{8}\\.rds$"
rds_files <- list.files(MERGE_DIR, pattern = rds_pattern, full.names = TRUE)
if (length(rds_files) == 0)
  stop("No Combined_Harmonized_<date>.rds found under ", MERGE_DIR,
       "\nRun 03_Merge/merge.R first.")
rds_file <- rds_files[which.max(file.info(rds_files)$mtime)]
today    <- format(Sys.Date(), "%Y%m%d")
out_pdf  <- file.path(REPORT_DIR,
                      sprintf("validate_Harmonized_%s.pdf", today))
out_tsv  <- file.path(REPORT_DIR,
                      sprintf("validate_Harmonized_%s_summary.tsv", today))

cat("Reading:", rds_file, "\n")
df <- readRDS(rds_file)
cat(sprintf("%d rows x %d columns\n", nrow(df), ncol(df)))

cohort_levels <- c("BASE-II", "SHIP", "SOEP", "TwinLife")
df$cohort <- factor(df$cohort,
                    levels = intersect(cohort_levels, unique(df$cohort)))

pgi_traits <- c("Edu", "Cog", "nonCog", "Height")

# ---- Summary TSV -----------------------------------------------------------
# Per-cohort diagnostics in wide format. Includes:
#   - n, demographic counts, phenotype means/SDs
#   - n / mean / SD for each PGI
#   - PGI x phenotype correlations + n (the scientific validation)
#   - inter-PGI correlations within cohort (sanity check)
#   - incremental R^2 of PGI added to a base model of covariates + PCs
safe_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

# Validation runs on the raw PLINK PGI scores (PGI_<trait>_raw), not the
# z-scored bare PGI_<trait> columns produced by normalize_pgi_columns().
# Linear transformations leave correlations and incremental R^2 unchanged,
# but raw means/SDs are the more useful sanity-check (they should match
# what the per-cohort PGI tables had).
pgi_cols   <- intersect(paste0("PGI_", pgi_traits, "_raw"), names(df))
pheno_cols <- intersect(c("education", "height"),           names(df))

# Covariates for the incremental R^2 base model. Whatever is present in df
# at runtime gets used; missing covariates are silently skipped.
covar_cols_all <- intersect(c("gender", "birth_year", paste0("PC", 1:20)),
                            names(df))

incremental_r2 <- function(d, pgi_col, outcome_col, covars) {
  # Drop covariates that are all-NA or constant within this subset
  covars <- covars[vapply(covars, function(c) {
    v <- d[[c]]
    !all(is.na(v)) && length(unique(v[!is.na(v)])) > 1
  }, logical(1))]
  use <- c(pgi_col, outcome_col, covars)
  ok  <- complete.cases(d[, use, drop = FALSE])
  if (sum(ok) < 30) return(NA_real_)
  d_ok <- d[ok, , drop = FALSE]
  cv_str <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
  full <- tryCatch(
    summary(lm(as.formula(sprintf("%s ~ %s + %s", outcome_col, pgi_col, cv_str)),
               data = d_ok))$r.squared,
    error = function(e) NA_real_)
  base <- tryCatch(
    summary(lm(as.formula(sprintf("%s ~ %s", outcome_col, cv_str)),
               data = d_ok))$r.squared,
    error = function(e) NA_real_)
  full - base
}

build_diagnostics <- function(d) {
  out <- list(
    n            = nrow(d),
    n_education  = sum(!is.na(d$education)),
    n_height     = sum(!is.na(d$height)),
    n_gender     = sum(!is.na(d$gender)),
    n_east_west  = sum(!is.na(d$east_west)),
    edu_mean     = mean(d$education, na.rm = TRUE),
    edu_sd       = sd(d$education,   na.rm = TRUE),
    height_mean  = mean(d$height,    na.rm = TRUE),
    height_sd    = sd(d$height,      na.rm = TRUE)
  )
  for (p in pgi_cols) {
    out[[paste0("n_",   p)]]     <- sum(!is.na(d[[p]]))
    out[[paste0(p, "_mean")]]    <- mean(d[[p]], na.rm = TRUE)
    out[[paste0(p, "_sd")]]      <- sd(d[[p]],   na.rm = TRUE)
  }
  for (p in pgi_cols) for (o in pheno_cols) {
    out[[sprintf("r_%s_%s", p, o)]] <- safe_cor(d[[p]], d[[o]])
    out[[sprintf("n_%s_%s", p, o)]] <- sum(!is.na(d[[p]]) & !is.na(d[[o]]))
  }
  if (length(pgi_cols) >= 2) {
    for (pair in combn(pgi_cols, 2, simplify = FALSE)) {
      out[[sprintf("r_%s_%s", pair[1], pair[2])]] <-
        safe_cor(d[[pair[1]]], d[[pair[2]]])
    }
  }
  for (p in pgi_cols) for (o in pheno_cols) {
    out[[sprintf("incr_r2_%s_%s", p, o)]] <-
      incremental_r2(d, p, o, covar_cols_all)
  }
  as.data.frame(out, check.names = FALSE)
}

summary_df <- df %>%
  group_by(cohort) %>%
  group_modify(~ build_diagnostics(.x)) %>%
  ungroup()

write.table(summary_df, out_tsv, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Summary TSV: ", out_tsv, "\n")

# ---- Helpers ---------------------------------------------------------------
text_page <- function(title, lines) {
  par(mar = c(0.5, 0.5, 2.5, 0.5))
  plot.new()
  title(main = title)
  text(0, 1, paste(lines, collapse = "\n"),
       adj = c(0, 1), family = "mono", cex = 0.85)
}

# ---- PDF -------------------------------------------------------------------
pdf(out_pdf, width = 10, height = 7)

## Page 1 — Summary
text_page(
  sprintf("Validation report — %s", today),
  c(
    sprintf("Source RDS:   %s", rds_file),
    sprintf("Dimensions:   %d rows x %d columns", nrow(df), ncol(df)),
    sprintf("Build date:   %s", attr(df, "build_date") %||% "unknown"),
    "",
    "Cohort summary:",
    capture.output(print(as.data.frame(summary_df), row.names = FALSE))
  )
)

## Pages 2..5 — PGI histograms per cohort (raw PLINK score)
for (t in pgi_traits) {
  col <- paste0("PGI_", t, "_raw")
  if (!col %in% names(df)) next
  sub <- df %>% filter(!is.na(.data[[col]]))
  if (nrow(sub) == 0) next
  p <- ggplot(sub, aes(x = .data[[col]])) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.85) +
    facet_wrap(~ cohort, scales = "free") +
    labs(title = sprintf("Histogram: %s (raw PLINK score)", col),
         x = col, y = "Count") +
    theme_minimal(base_size = 11)
  print(p)
}

## Pages 6..7 — Outcome histograms
for (outcome in c("education", "height")) {
  if (!outcome %in% names(df)) next
  sub <- df %>% filter(!is.na(.data[[outcome]]))
  if (nrow(sub) == 0) next
  p <- ggplot(sub, aes(x = .data[[outcome]])) +
    geom_histogram(bins = 40, fill = "darkorange", alpha = 0.85) +
    facet_wrap(~ cohort, scales = "free") +
    labs(title = sprintf("Histogram: %s", outcome),
         x = outcome, y = "Count") +
    theme_minimal(base_size = 11)
  print(p)
}

## Pages 8..9 — PGI × outcome scatters (raw PLINK score)
make_scatter <- function(df, outcome) {
  present <- intersect(paste0("PGI_", pgi_traits, "_raw"), names(df))
  long <- df %>%
    select(cohort, all_of(c(outcome, present))) %>%
    pivot_longer(all_of(present), names_to = "trait",
                 values_to = "pgi") %>%
    mutate(trait = sub("_raw$", "", trait)) %>%
    filter(!is.na(pgi), !is.na(.data[[outcome]]))
  if (nrow(long) == 0) return(NULL)

  rs <- long %>%
    group_by(trait, cohort) %>%
    summarise(r = suppressWarnings(cor(pgi, .data[[outcome]],
                                       use = "pairwise.complete.obs")),
              .groups = "drop") %>%
    group_by(trait) %>%
    summarise(label = paste(sprintf("%s: r=%.2f", cohort, r),
                            collapse = "   "),
              .groups = "drop")

  ggplot(long, aes(x = pgi, y = .data[[outcome]], colour = cohort)) +
    geom_point(alpha = 0.10, size = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    facet_wrap(~ trait, scales = "free_x") +
    geom_text(data = rs, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.4,
              size = 3, colour = "grey20") +
    labs(title    = sprintf("%s ~ PGIs by cohort", outcome),
         subtitle = "Linear fit per cohort; r annotated top-left",
         x = "PGI (raw PLINK score)",
         y = outcome,
         colour = "Cohort") +
    theme_minimal(base_size = 11)
}

if ("education" %in% names(df)) {
  p <- make_scatter(df, "education"); if (!is.null(p)) print(p)
}
if ("height" %in% names(df)) {
  p <- make_scatter(df, "height"); if (!is.null(p)) print(p)
}

## Page 9b — Incremental R^2 of PGI controlling for covariates + PCs
incr_grid <- expand.grid(cohort  = levels(df$cohort),
                         pgi     = pgi_cols,
                         outcome = pheno_cols,
                         stringsAsFactors = FALSE)
incr_grid$incr_r2 <- mapply(function(coh, p, o) {
  sub <- df[df$cohort == coh, , drop = FALSE]
  incremental_r2(sub, p, o, covar_cols_all)
}, incr_grid$cohort, incr_grid$pgi, incr_grid$outcome)
incr_grid$pgi_short <- sub("^PGI_|_raw$", "", incr_grid$pgi)
incr_grid$cohort    <- factor(incr_grid$cohort, levels = levels(df$cohort))

if (any(!is.na(incr_grid$incr_r2))) {
  p_incr <- ggplot(incr_grid,
                   aes(x = cohort, y = pmax(incr_r2, 0) * 100,
                       fill = pgi_short)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    geom_text(aes(label = sprintf("%.2f", incr_r2 * 100)),
              position = position_dodge(width = 0.85),
              vjust = -0.4, size = 2.7, na.rm = TRUE) +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(title    = "Incremental R^2 of PGI (added to covariates + PCs)",
         subtitle = sprintf("Base model: %s",
                            paste(covar_cols_all, collapse = " + ")),
         x = "", y = "Delta R^2 (%)", fill = "PGI") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 25, hjust = 1))
  print(p_incr)
}

## Page 10 — PGI inter-correlation matrices (on raw PLINK scores;
## linear scaling does not change correlations)
cor_rows <- list()
for (coh in levels(df$cohort)) {
  sub <- df %>% filter(cohort == coh) %>%
    select(any_of(paste0("PGI_", pgi_traits, "_raw")))
  if (nrow(sub) < 30 || ncol(sub) < 2) next
  names(sub) <- sub("_raw$", "", names(sub))
  cm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs"))
  cor_rows[[coh]] <- cm %>% as.data.frame() %>%
    tibble::rownames_to_column("trait") %>%
    pivot_longer(-trait, names_to = "other", values_to = "r") %>%
    mutate(cohort = coh)
}
if (length(cor_rows) > 0) {
  cor_long <- bind_rows(cor_rows)
  p <- ggplot(cor_long, aes(x = trait, y = other, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
    scale_fill_gradient2(low = "red", mid = "white", high = "steelblue",
                         midpoint = 0, limits = c(-1, 1)) +
    facet_wrap(~ cohort) +
    labs(title = "PGI inter-correlations by cohort",
         x = "", y = "", fill = "r") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  print(p)
}

## Page 11 — TwinLife parental-education reliability
if ("TwinLife" %in% levels(df$cohort) &&
    all(c("parental_education", "fid", "ptyp", "education") %in% names(df))) {
  tl <- df %>% filter(cohort == "TwinLife") %>%
    mutate(ptyp_chr = as.character(ptyp),
           fid_chr  = as.character(fid))
  parent_ptyp <- c("300", "400")
  child_ptyp  <- c("1", "2", "110", "120", "200", "201")
  parents_own <- tl %>% filter(ptyp_chr %in% parent_ptyp,
                               !is.na(education)) %>%
    group_by(fid_chr) %>%
    summarise(mean_parent_own_edu = mean(education, na.rm = TRUE),
              n_parents_with_edu  = sum(!is.na(education)),
              .groups = "drop")
  children <- tl %>% filter(ptyp_chr %in% child_ptyp,
                            !is.na(parental_education))
  cross <- inner_join(
    children %>% select(fid_chr, pid, ptyp_chr, parental_education),
    parents_own, by = "fid_chr"
  )
  if (nrow(cross) >= 30) {
    r <- suppressWarnings(cor(cross$parental_education,
                              cross$mean_parent_own_edu,
                              use = "pairwise.complete.obs"))
    d <- cross$parental_education - cross$mean_parent_own_edu
    text_page(
      "TwinLife parental_education reliability",
      c("(child-row parental_education vs mean of parents' own education",
        "within the same FID)",
        "",
        sprintf("n child rows with parents-own-edu:  %d", nrow(cross)),
        sprintf("cor(child_pared, mean_parent_own):  %.3f", r),
        sprintf("mean diff (child - parent):         %.3f",
                mean(d, na.rm = TRUE)),
        sprintf("sd   diff:                          %.3f",
                sd(d,   na.rm = TRUE)),
        sprintf("%% within 0.5 years:                 %.1f%%",
                100 * mean(abs(d) < 0.5, na.rm = TRUE)))
    )
  } else {
    text_page("TwinLife parental_education reliability",
              "Not enough paired child/parent rows in the analysis sample.")
  }
}

## Page 12 — ID prefix uniqueness
id_check <- df %>%
  mutate(has_prefix = grepl("^(b_|sh_|so_|tl_)", pid)) %>%
  summarise(n_total          = n(),
            n_missing_prefix = sum(!has_prefix),
            n_duplicate_pids = sum(duplicated(pid)))
text_page(
  "ID prefix & uniqueness check",
  c(sprintf("Total rows:               %d", id_check$n_total),
    sprintf("Missing b_/sh_/so_/tl_:   %d   (must be 0)",
            id_check$n_missing_prefix),
    sprintf("Duplicate pid:            %d   (must be 0)",
            id_check$n_duplicate_pids))
)

## Page 13 — Ancestry PC scree (variance explained per PC)
## TRUE scree: step 04 emits the PCA eigenvalues as {COHORT}_PC_eigenval.tsv in
## the deliverable dir; lambda_k / sum(lambda) is the real, reproducible,
## scale-invariant "% variance explained" per ordered PC. If those files are
## absent, fall back to the share of PROJECTED-score variance. That fallback is
## only a weak sanity shape: the projection emits ~equal-variance PCs, so
## var(PC) does not carry eigenvalue magnitude and its raw scale is an
## arbitrary, PLINK-version-dependent artifact (~1e9 different between builds).
pretty_cohort <- function(x) {
  m <- c(BASEII = "BASE-II", SHIP0 = "SHIP-0", SHIPTD = "SHIP-Td",
         SOEP = "SOEP", TWINLIFE = "TwinLife")
  out <- m[x]; out[is.na(out)] <- x[is.na(out)]; unname(out)
}
eig_dir   <- file.path(DATA_ROOT, "01_Genotype", "Harmonized", "data", "final")
eig_files <- list.files(eig_dir, pattern = "_PC_eigenval\\.tsv$", full.names = TRUE)
if (length(eig_files) > 0) {
  eig <- do.call(rbind, lapply(eig_files, read.delim, stringsAsFactors = FALSE)) %>%
    mutate(cohort = pretty_cohort(cohort),
           PC_num = as.integer(sub("^PC", "", PC))) %>%
    group_by(cohort) %>%
    mutate(pct = 100 * eigenval / sum(eigenval, na.rm = TRUE)) %>%
    ungroup()
  p <- ggplot(eig, aes(x = PC_num, y = pct, colour = cohort)) +
    geom_line() + geom_point(size = 1.5) +
    scale_x_continuous(breaks = seq(0, max(eig$PC_num, na.rm = TRUE), by = 5)) +
    labs(title    = "Ancestry PC scree — variance explained per PC",
         subtitle = "PCA eigenvalue share (lambda_k / sum lambda) per cohort — true, reproducible scree",
         x = "PC index (ordered)", y = "Variance explained (%)", colour = "Cohort") +
    theme_minimal(base_size = 11)
  print(p)
} else {
  pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
  pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
  if (length(pc_cols) > 0) {
    pc_scree <- df %>% group_by(cohort) %>%
      summarise(across(all_of(pc_cols), ~ var(.x, na.rm = TRUE)), .groups = "drop") %>%
      pivot_longer(-cohort, names_to = "PC", values_to = "variance") %>%
      group_by(cohort) %>%
      mutate(share  = 100 * variance / sum(variance, na.rm = TRUE),
             PC_num = as.integer(sub("^PC", "", PC))) %>%
      ungroup()
    p <- ggplot(pc_scree, aes(x = PC_num, y = share, colour = cohort)) +
      geom_line() + geom_point(size = 1.5) +
      scale_x_continuous(breaks = seq(0, length(pc_cols), by = 5)) +
      labs(title    = "Ancestry PC scree — share of PROJECTED-score variance (FALLBACK)",
           subtitle = "eigenvalues unavailable; re-run step 04 to emit {COHORT}_PC_eigenval.tsv",
           x = "PC index (ordered)", y = "Share of PC variance (%)", colour = "Cohort") +
      theme_minimal(base_size = 11)
    print(p)
  }
}

dev.off()
cat("PDF written: ", out_pdf, "\n")
