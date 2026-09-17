# A_Education — Section A: Educational attainment (cohort, region, gender)

Section A is organized **by phenotype** (educational attainment): it reports the birth-cohort,
East/West and **gender** moderators of the PGI–education association together. Gender is a main
dimension here — the Section C analyses of the preregistration (C1/C3/C0/C2) estimated with family
random effects, plus a per-gender attainment GAM (A5); see `../../../Plan_deviations.md` §10/§13.

- **Sample:** `dat_cluster` (full-family TwinLife, random effects via `s(fid_re, bs = "re")`; see
  `../Random_effects_family_clustering.md`). A-spec (below) is the one exception — an NLS specification
  check with no RE analogue, kept on the dedup `dat_edu`.
- **Outcome:** `edu_z_kernel` (years of education standardized within birth year via a Gaussian kernel,
  pooled across regions; see `../../../Plan_deviations.md` §2/§3).

## Analyses

| id | what | model | preregistration |
|---|---|---|---|
| **A1** | linear PGI × birth-year × region | linear RE | Section A |
| **A2** | step-function PGI × reunification × region | step RE | Section A |
| **A3** | region-specific PGI cohort smooths | GAM RE | Section A |
| **A4** | four-way PGI × birth-year × region × gender | linear RE | Section C (C1) |
| **A5** | per-gender PGI–education slope over birth year (regions pooled) | GAM RE | per-gender attainment GAM |
| **A6** | Region × Gender PGI cohort smooths | GAM RE | Section C (C3) |
| **A0g** | overall Region × Gender PGI contrast (+ pre-reunification) | RE contrast | Section C (C0) |
| **A-spec** | logistic-transition specification check | NLS (dedup) | Section C (C2) |

Where the higher-order interactions (A4 four-way, A6 Region × Gender, A0g three-way) are null, the
lower-order / additive gender moderation is the reported effect.

## Files

| File | Role |
|---|---|
| `A_descriptives.R` | attainment descriptives over strata (sourced by `A_analysis.R`) |
| `A_analysis.R` | A1–A3 + A4/A5/A6/A0g/A-spec; per-cell nonlinearity LRTs; leave-one-study-out (incl. the gender LOCO tables `A4_LOO_gender`, `A5_LOO`, §14); saves `A_Education_Models.rds` |
| `A_report.R` | builds `Supplementary_Data_2_Education_results.xlsx`, `A_Education_Plots.pdf`, `A_Education_LOO_Plots.pdf`, `A_Education_Overview.md` |

Leave-one-study-out (LOCO) for the **pooled** `PGI × gender` term (`A4_LOO_gender`) and the
per-gender GAM (`A5_LOO`) accompany the four-way `A4_LOO` as an additional sensitivity
(`../../../Plan_deviations.md` §14). The finding-matched random-effects alternative-PGI checks
(`S1_finding1` / `S1_finding3_att` / `S1_A5_nonlin`) live in the sensitivities script (§15).

Run via `../run_A.command`. (Descriptives for both phenotypes live in `../00_Descriptives`.)

## Outputs (under `runs/RUN_<TS>/A_Education/`)

- `Supplementary_Data_2_Education_results.xlsx` — sheets `00_Index`, `01_Linear_birth_year`,
  `02_Step_reunification`, `03_GAM_region`, `04_Linear_gender`, `05_GAM_gender`,
  `06_GAM_region_gender`, `07_Overall_RxG`, `08_Logistic_check` (human-readable term/column labels)
- `A_Education_Plots.pdf` (incl. A5 per-gender and A6 Region × Gender slope curves + diff smooths)
- `A_Education_LOO_Plots.pdf`
- `A_Education_Models.rds`, `A_Education_Overview.md`, `A_run.log`

## Sensitivities (Section A owns its attainment checks)

`A_sensitivities_analysis.R` → `A_Sensitivities_Models.rds` → `A_sensitivities_report.R`
(`Supplementary_Data_3_Education_sensitivities_results.xlsx` + `A_Sensitivities_Plots.pdf` +
`A_Sensitivities_Overview.md`):

| id | what | sample |
|---|---|---|
| **S1-A** | alternative PGIs (Cog / NonCog) on the A1 and A4 focals, BH per family | dedup-OLS (`dat_edu`) |
| **S1-A (finding-matched)** | alternative PGIs on the reported focals, family RE, common re-residualised sample | `dat_cluster` |
| **S2-A** | height-PGI negative control: → education (A1 spec), → height (sanity) | dedup-OLS |
| **S4** | migration diagnostics (A-only), gated on A focal significance + gated S4d within-region refit | dedup-OLS |
| **Method** | A1 focal under dedup-OLS / CR2 / RE (RE read from `A_Education_Models.rds`) | dedup-OLS + `dat_cluster` |
| **S3-A** | heteroscedasticity (Gaussian location-scale `gaulss`) per region, gated on A focal significance | `dat_cluster` (no RE) |

Sensitivities stay dedup-OLS (relative comparison across PGIs is unaffected by the estimator;
`../../../Plan_deviations.md` §11) except the finding-matched S1-A, CR2/RE and the `gaulss` S3-A
(full-family, no RE; §12).

## Per-study appendix (attainment)

`A_appendix_analysis.R` → `A_Appendix_Models.rds` → `A_appendix_report.R`
(`Supplementary_Data_4_Education_appendix_results.xlsx` + `A_Appendix_Plots.pdf` +
`A_Appendix_Overview.md`):

| id | what |
|---|---|
| **D** | per-study A3-style education GAM (region-specific PGI cohort smooths) |
| **F** | per-study A6-style Region × Gender GAM on `edu_z_kernel` |
| **F-linear** | per-study robust linear gender model (PGI main / gender × PGI / Region × Gender × PGI) |

Appendix-grade (narrow within-study birth-year windows; some study × cell combos skipped). The mobility
per-study analysis (E) lives in Section B. See `../../../Plan_deviations.md` §13.

## Run order

`../run_A.command` runs all six stages under one `RUN_TS`: A_analysis → A_report →
A_sensitivities_analysis → A_sensitivities_report → A_appendix_analysis → A_appendix_report. Each stage
is independently runnable (the sensitivities read `A_Education_Models.rds` for their gates).
