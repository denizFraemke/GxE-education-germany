# 01_PGI_Analysis

PGI-level workflow for the manuscript *Genomics of Educational Attainment
Across 80 Years of Social and Political Transformation in Germany*.

The workflow is organized **by phenotype** (not by test type): each
phenotype section reports its birth-cohort, East/West and **gender**
moderators together, and owns its sensitivities and per-study appendix.

Methodological departures are catalogued in `../../Plan_deviations.md`.

## Layout

```
01_PGI_Analysis/
  run_A.command, run_B.command, run_Descriptives.command
                              # independent per-section launchers
  00_setup/                   # data loading, sample prep, run context
  R/                          # shared helpers + GxE_Germany_analysis_helpers.R wrapper
  00_Descriptives/            # sample descriptives (attainment + mobility)
  A_Education/                # Section A: educational attainment
                              #   A1/A2/A3 (cohort, region) + A4/A5/A6/A0g/A-spec (gender, RE)
                              #   A_sensitivities_* (S1-A/S2-A/S4/Method/S3-A)
                              #   A_appendix_*       (per-study D/F/F-linear)
  B_Mobility/                 # Section B: educational mobility
                              #   B0/B0_pre/B1/B2/B3 (cohort, region) + B4/B5/B6/B0g (gender, RE)
                              #   B_sensitivities_*  (S1-B/S2-B/S3-B)
                              #   B_appendix_*       (per-study E/E-RxG/F-linear)
  Random_effects_family_clustering.md   # RE vs dedup-OLS methods note
```

Gender is a main moderator inside Sections A and B (the preregistration's
Section C analyses, estimated with family random effects), so there is no
standalone gender folder. The cross-section sensitivities and the per-study
appendix likewise live inside each phenotype section.

## Run order

Each launcher is independent and runs all of its section's stages under one
`RUN_TS`:

```bash
./run_Descriptives.command          # Supplementary_Data_1_Descriptives_results.xlsx
./run_A.command                     # Supplementary_Data_2/3/4 (Education: main / sensitivities / appendix)
./run_B.command                     # Supplementary_Data_5/6/7 (Mobility: main / sensitivities / appendix)
```

### Rebuilding a workbook without the PDFs

The plot stage of `A_report.R` / `B_report.R` runs *after* the workbook is
written. For a wording or layout change to the xlsx, skip it:

```bash
REPORT_SKIP_PLOTS=1 RUN_TS=<ts> Rscript A_Education/A_report.R
```

The xlsx and the `*_Overview.md` are written either way and are unchanged by the
flag. The existing PDFs are left in place, so re-run without the flag before a
release if the figures also need refreshing.

`run_A.command` runs six stages (main analysis → report → sensitivities
analysis → sensitivities report → appendix analysis → appendix report);
`run_B.command` runs seven (the parental-education moderation stage is inserted
after the main analysis). The sensitivities read the section's `*_Models.rds`
from the same run for their significance gates, so the stages share one
`RUN_TS`.

By default each launcher creates its own timestamped run directory under
`runs/RUN_<YYYYMMDD_HHMM>/<section>/`. Export `RUN_TS` before invoking to run
multiple sections into a single timestamp:

```bash
export RUN_TS=$(date +%Y%m%d_%H%M)
./run_Descriptives.command
./run_A.command
./run_B.command
```

## Output convention

Outputs land at
`${DATA_ROOT}/04_GxE_Analysis/01_PGI_Analysis/runs/RUN_<TS>/`, where
`DATA_ROOT` is the data share (default `/path/to/data_root`; override with the
`DATA_ROOT` environment variable). Each run is preserved in its own timestamped
directory; nothing is overwritten.

## Deviations from the preregistered plan

See `../../Plan_deviations.md`. Key publication-readiness deviations:

- gender is a family-RE main moderator inside A/B (preregistration
  Section C; §10), with a per-gender attainment/mobility GAM (A5/B5)
- per-cell nonlinearity decided by a clean nested LRT; plots draw the smooth
  only where that LRT is significant
- sensitivities and the per-study appendix live inside each phenotype section
- S3 (heteroscedasticity) and S4 (migration) conditional on focal-test
  significance (§1); per-focal dispersion robustness gated likewise
- wf-PGI omitted from S1 (§8)
- SNP-heritability is outside this workflow (lives in `../02_SNP_Heritability/`)
