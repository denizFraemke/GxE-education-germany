# 80 Years of GxE on Education in Germany

Monorepo for the gene-by-environment analysis of educational attainment across
four German cohorts (BASE-II, SHIP, SOEP, TwinLife). Eligibility was defined for
birth years 1920–2000; the realised analytic sample spans birth years 1918–1995.

The pipeline uses cross-cohort harmonized SBayesR PGIs: all cohorts are scored
against a common SBayesR posterior so the resulting PGIs are directly
comparable across cohorts.

## Repository layout

Top-level folders are numbered to reflect pipeline order. Steps annotated
`(Tardis)` run on the HPC cluster; everything else runs on a Mac with the
shared volume mounted.

```
.
├── 01_Genotype/
│   └── Harmonized/           # Cross-cohort harmonized SBayesR PGI pipeline (Tardis)
├── 02_Phenotype/             # Per-cohort phenotype construction
│   ├── BASE-II/ SHIP/ SOEP/ TwinLife/
│   └── utils/                # edu_helpers.R (shared education coding)
├── 03_Merge/
│   ├── merge.R               # Orchestrator: phenotypes + PGIs + PCs → combined RDS
│   ├── validate.R            # Post-merge PDF/TSV report
│   └── R/                    # merge_common.R, paths.R
├── 04_GxE_Analysis/
│   ├── 01_PGI_Analysis/      # Mega-analytic GxE workflow (PGI × cohort × region × gender)
│   │   └── 06_Sensitivity_analyses/    # Further sensitivity analyses
│   └── 02_SNP_Heritability/  # GREML SNP-h² by Region × Time (Tardis); Analysis S5
├── 05_manuscript_export/     # Finished analysis run → tidy aggregate CSVs
│   ├── export_result_tables.R
│   ├── export_manuscript_aggregates.R
│   └── check_frozen_run.R
├── ANALYSIS_PLAN.md
├── Plan_deviations.md
├── README.md
└── .gitignore
```

`02_Phenotype/` is tracked directly in this repository; it is not a git
submodule and there is no `.gitmodules`.

## Data location

All data products live on the shared volume, mirroring the code layout. The
root is written `${DATA_ROOT}` below; the documented default is
`/path/to/data_root`.

```
${DATA_ROOT}/
├── 01_Genotype/Harmonized/data/final/{COHORT}_PGI_PCs.tsv
├── 02_Phenotype/                                          (processed phenotype outputs)
├── 03_Merge/
│   ├── Combined_Harmonized_<YYYYMMDD>.rds
│   └── validation_reports/
└── 04_GxE_Analysis/
    ├── 01_PGI_Analysis/runs/RUN_<TS>/                     (xlsx / pdf / rds per section)
    └── 02_SNP_Heritability/                               (h2_snp.tsv / lrt.tsv / power.tsv)
```

The Mac-side R stages — `03_Merge`, `04_GxE_Analysis/01_PGI_Analysis`, the R
helpers of `04_GxE_Analysis/02_SNP_Heritability`, and `05_manuscript_export` —
resolve this root through a `data_root()` helper defined once per stage
(`03_Merge/R/paths.R`,
`04_GxE_Analysis/01_PGI_Analysis/00_setup/run_context.R`,
`04_GxE_Analysis/02_SNP_Heritability/R/paths.R`). Set the `DATA_ROOT`
environment variable to point those stages at a different root; don't hard-code
the path in individual scripts. The genotype stage does **not** read
`DATA_ROOT`: `01_Genotype/Harmonized/config.sh` derives `PROJ_ROOT` from its own
location (overridable through `user_config.sh`), and the Tardis side of
`04_GxE_Analysis/02_SNP_Heritability` locates the upstream genotype outputs
through `GENOTYPE_PROJ_ROOT`.

No individual-level or source cohort data is committed to this repository. The
manuscript-export stage does produce aggregate figure-source CSVs; those are
coefficients, curves, M(SD) tables and counts.

## Pipeline overview

1. **Genotype** *(Tardis)* — run the PGI pipeline (see
   `01_Genotype/Harmonized/README.md` in this repo). Final per-cohort TSVs are
   transferred to `${DATA_ROOT}/01_Genotype/Harmonized/data/final/`.
2. **Phenotype** *(per-cohort sources)* — run the per-cohort scripts under
   `02_Phenotype/`. Processed phenotype outputs live under
   `${DATA_ROOT}/02_Phenotype/`.
3. **Merge** *(Mac + shared volume)* — `Rscript 03_Merge/merge.R` writes
   `${DATA_ROOT}/03_Merge/Combined_Harmonized_<date>.rds`. Raw PGIs and ancestry
   PCs flow through unchanged; no residualization happens at this stage.
4. **Validate** *(Mac + shared volume)* — `Rscript 03_Merge/validate.R`
   inspects the latest merged RDS and writes a PDF report plus a summary TSV
   to `${DATA_ROOT}/03_Merge/validation_reports/`.
5. **Analysis — GxE mega-analysis** *(Mac + shared volume)* — launch the
   per-section launchers under `04_GxE_Analysis/01_PGI_Analysis/`
   (`run_Descriptives.command`, `run_A.command`, `run_B.command`), each of which
   runs its analysis + report stages and writes a timestamped log alongside the
   outputs. The PGIs are residualized on ancestry PCs as part of sample
   preparation.
6. **Analysis — SNP-heritability (S5)** *(Mac → Tardis)* — Mac side derives
   GCTA inputs from the merged RDS and rsyncs both code and phenotype/strata
   to Tardis; Tardis runs the multi-component GREML pipeline. See
   [`04_GxE_Analysis/02_SNP_Heritability/README.md`](04_GxE_Analysis/02_SNP_Heritability/README.md)
   for the end-to-end workflow.
7. **Manuscript export** *(Mac + shared volume)* — `05_manuscript_export/` turns a
   finished run into tidy, aggregate-only CSVs under `RUN_<TS>/manuscript_export/`
   (`export_result_tables.R` for every result table, `export_manuscript_aggregates.R`
   for figure-source data + descriptives). The figure scripts that draw those CSVs
   are maintained in the separate manuscript repository, which imports them. No
   individual-level data leaves the analysis repo. See
   [`05_manuscript_export/README.md`](05_manuscript_export/README.md).

Documentation of implementation choices that depart from `ANALYSIS_PLAN.md`
lives in [`Plan_deviations.md`](Plan_deviations.md).

## The frozen run

Every result the manuscript reports comes from a single frozen analysis run,
`RUN_20260729_1555`. `05_manuscript_export/check_frozen_run.R` verifies a
checkout against it: that the merged input is the one the run used, that the
model caches still carry their recorded fingerprints, that the analytic samples
rebuild to the same Ns, and that a set of headline estimates reads back
unchanged.

[`04_GxE_Analysis/01_PGI_Analysis/06_Sensitivity_analyses/`](04_GxE_Analysis/01_PGI_Analysis/06_Sensitivity_analyses/README.md)
holds three further sensitivity analyses — a scale-invariance and censoring
ladder, sample crossed with the scope of the ancestry-PC residualisation, and a
refit of the published focal models with ancestry-PC main effects and PC ×
moderator terms — alongside the script that counts the per-study birth-year Ns
behind Figure 1. None of them refits or overwrites the frozen run; they read it
and write their own aggregate CSVs, so the reported estimates are untouched.

## Dependencies

R ≥ 4.3 with `broom`, `clubSandwich`, `data.table`, `dplyr`, `ggplot2`,
`haven`, `lmtest`, `lubridate`, `MASS`, `mgcv`, `openxlsx`, `patchwork`,
`purrr`, `sandwich`, `stringr`, `tibble`, `tidyr`. There is no lockfile; each
analysis run bundles a `sessionInfo()` into its output folder for an
exact-version record.

## Licence

MIT — see [`LICENSE`](LICENSE).

## Citation

Fraemke, D., Miller, A., Koellinger, P., Hertwig, R., Richter, D., Zinn, S.,
Kandler, C., Forstner, A. J., Mönkediek, B., Diewald, M., Teumer, A.,
Baumeister, S. E., Völzke, H., Völker, U., Grabe, H. J., Bertram, L.,
Lindenberger, U., Drewelies, J., Kühn, S., Demuth, I., Gerstorf, D., Okbay, A.,
Biroli, P., Fuchs-Schündeln, N., Harden, K. P., Abdellaoui, A., Malanchini, M.,
Tucker-Drob, E. M., & Raffington, L. *Analysis code for: Genomics of
Educational Attainment Across 80 Years of Social and Political Transformation
in Germany* (Version 1.0-preprint) [Computer software].
https://github.com/Biosocial/GxE-education-germany

The preprint DOI is not yet issued. [`CITATION.cff`](CITATION.cff) carries the
same byline in machine-readable form, with a commented-out `identifiers:` block
to be filled in once the DOI exists.

Questions about the code go to Deniz Fraemke (fraemke@mpib-berlin.mpg.de).
