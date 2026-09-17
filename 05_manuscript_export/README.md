# 05_manuscript_export — bridge from the analysis run to the manuscript repo

Pipeline stage **05**: turn a finished GxE analysis run (`04_GxE_Analysis`, the
`RUN_<TS>/` model caches) into **tidy, aggregate-only CSVs** that the separate
manuscript repository consumes.

**Contract with the manuscript repo:** the manuscript repo never fits models and
never stores individual-level data. It reads *only* the CSVs this stage writes
into `RUN_<TS>/manuscript_export/`. There is no Excel in the loop — the CSVs are
generated straight from the analysis result tables (no `data → xlsx → CSV`
round-trip). If publication-ready `.xlsx` workbooks are needed for readers, they
are assembled *from* these same CSVs as a final step, never parsed back.

```
04_GxE_Analysis run (RUN_<TS>/*.rds model caches)
        │
        ├─ export_result_tables.R        → every aggregate result table, tidy CSV
        │                                   (coefficients, GAM tests, LOO, slopes, …)
        └─ export_manuscript_aggregates.R → figure-source data, fitted curves,
                                            descriptives, focal_estimates, provenance
        ▼
RUN_<TS>/manuscript_export/   ── import_results.R (--from <run>) ──►  manuscript repo (01_results/)
```

## Output layout (`RUN_<TS>/manuscript_export/`)

Organised so a human can scan it; every folder name says what it holds.

```
manuscript_export/
├── README.md                 auto-generated legend (analysis codes + file-type suffixes)
├── focal_estimates.csv       consolidated headline focal tests
├── _provenance.csv           run id / analysis-repo commit / date
├── snp_h2_pooled.csv         pooled SNP-h2 + SE/CI over the four cells
├── figS06_gender_pgi_*.csv   gender x PGI-quintile gap, SD-unit + years scales,
│                             and the collapsed PGI bands
├── gender_pgi_simple_slopes.csv       PGI slope within women / within men,
│                             region collapsed, both outcomes
├── descriptives_birthyear_trend.csv   birth-year trend + gender/region
│                             moderation of it
├── descriptives_variance_tests.csv    East-West and Female-Male variance
│                             contrasts by birth-year stratum
├── descriptives/             sample descriptives, PGI correlations, attrition
├── fig01_*.csv … figS06_*.csv         figure-source data, at the export root
├── attainment/               Section A — educational attainment
│   ├── main/                 A1–A6, A0g: *_Coefficients, *_ModelFit, *_NestedLRT, *_LOO, …
│   ├── sensitivities/
│   └── appendix/             per-cohort
└── mobility/                 Section B — educational mobility
    ├── main/                 B0–B6, B0g: *_Coefficients, *_RegionSlopes, *_NestedLRT, *_LOO, …
    ├── sensitivities/
    └── appendix/
```

File names carry both the analysis code and the table type, e.g.
`attainment/main/A1_Coefficients.csv`, `attainment/main/A3_NestedLRT.csv`,
`mobility/main/B0_RegionSlopes.csv`. The generated `README.md` in that folder
spells out the codes (A1 = PGI×BY×region, A4 = four-way gender, B0 = East–West
mobility contrast, …) and the suffixes (`_Coefficients`, `_NestedLRT` = GAM
nonlinearity LRT, `_LOO` = leave-one-cohort-out, …).

## Scripts

- **`check_frozen_run.R`** — reproducibility check. Run it before fitting or
  exporting anything against a frozen run. It verifies that this machine still
  reproduces the run the manuscript's frozen exports were built from: data root
  mounted, which `Combined_Harmonized_*.rds` `load_data.R` resolves to (it picks
  newest-by-mtime, so a newer merge would change the sample silently), the model
  caches' size+mtime fingerprints against the manuscript repo's
  `_provenance.csv`, the rebuilt analytic sample sizes against the frozen
  descriptives cache, and a headline coefficient. Exits non-zero on the first
  failure.
  ```
  RUN_TS=<ts> MANUSCRIPT_REPO=<path> Rscript 05_manuscript_export/check_frozen_run.R
  ```
- **`export_result_tables.R`** — dumps every aggregate table from the model
  caches (`A_Education_Models.rds`, sensitivities, appendix; the B equivalents;
  `Descriptives_Tables.rds`) into the layout above. Self-contained (only reads the
  caches). An allow-list restricts output to aggregate tables; individual-level
  frames (`dat_cluster*`, `df_*`) are **never** written.
  ```
  RUN_DIR=<…>/runs/RUN_<TS> OUT_DIR=<…>/runs/RUN_<TS>/manuscript_export \
    Rscript 05_manuscript_export/export_result_tables.R
  ```
- **`export_manuscript_aggregates.R`** — figure-source data, fitted GAM
  curves/slopes, descriptive tables, `focal_estimates.csv`, and `_provenance.csv`.
  Reuses the `04_GxE_Analysis/01_PGI_Analysis` setup + helpers (via `WORKFLOW_DIR`)
  to rebuild fitted curves, so it depends on that analysis stage. Never writes
  individual rows.
  ```
  RUN_TS=<ts> DATA_ROOT=<path> Rscript 05_manuscript_export/export_manuscript_aggregates.R
  ```
  **Selective re-export.** `EXPORT_ONLY="<substring>,<substring>"` runs only the
  steps whose label contains one of the substrings and leaves every other
  already-exported CSV untouched. Use this when adding one output to a run whose
  other files are already frozen and imported downstream — a full re-run would
  rewrite and re-checksum all of them for no reason. Example:
  ```
  EXPORT_ONLY="figS06 gender x PGI,pooled SNP-heritability" \
    RUN_TS=<ts> Rscript 05_manuscript_export/export_manuscript_aggregates.R
  ```
The figures themselves are built in the manuscript repository, from the CSVs
this stage exports; the plot scripts live there under
`02_manuscript/figures/scripts/`.

## Transfer to the manuscript repo

`import_results.R --from "<run folder>" manuscript_export/<files…>` copies the
CSVs into the manuscript repo's `01_results/` (with SHA-256 + provenance). Because
`--from` points at the SMB run folder (not a git repo), commit/branch come back
`NA`; `_provenance.csv` carries the run id + analysis-repo commit.

## Privacy

Every output is a coefficient, curve, slope, smooth, M(SD) table, or count — no
individual rows. The `.rds` caches (which embed the analytic data frames) stay in
the analysis repo and are never transferred.
