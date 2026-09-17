# B_Mobility — Section B: Educational mobility (cohort, region, gender)

Section B is organized **by phenotype** (intergenerational educational mobility): it reports the
birth-cohort, East/West and **gender** moderators of the PGI→mobility association together, controlling
for parental education throughout. Gender is a main dimension here — the Section C analyses of the
preregistration (C1/C3/C0) estimated with family random effects, plus a per-gender mobility GAM (B5);
see `../../../Plan_deviations.md` §9/§10/§13.

- **Sample:** `dat_cluster_mob` (full-family TwinLife, random effects via `s(fid_re, bs = "re")`; see
  `../Random_effects_family_clustering.md`). SHIP has no parental education and is **absent** from this
  sample, so SHIP is not a LOO fold; the LOO folds are BASE-II / SOEP / TwinLife.
- **Outcome:** `mobility = edu_z_kernel − parental_edu_z_kernel` (each standardized within birth year via
  a Gaussian kernel; see `../../../Plan_deviations.md` §2/§3).
- **Parental-education control** (see `../../../Plan_deviations.md` §9/§10):
  - **B1 / B2** use the full preregistered Eq-4 ParEdu interaction set.
  - **B3** controls ParEdu with a main effect, a `ParEdu:east_west_c` term, and region-specific ParEdu
    smooths over birth year (`s(BYc, by = ParEdu_W) + s(BYc, by = ParEdu_E)`), parallel to the PGI
    smooths. This flips the focal PGI nonlinearity from null to West-driven significant relative to an
    unadjusted B3 (**fragile / secondary** — also fails leave-one-cohort-out; see the B3 caveat in the
    workbook and overview).
  - **B4 / B5 / B6 / B0g** control ParEdu at the same cell flexibility as the focal PGI cells.

## Analyses

| id | what | model | preregistration |
|---|---|---|---|
| **B0 / B0_pre** | overall (and pre-reunification) East-West PGI-mobility contrast | linear RE | overall contrast |
| **B1** | linear PGI × birth-year × region | linear RE | Section B |
| **B2** | step-function PGI × reunification × region | step RE | Section B |
| **B3** | region-specific PGI cohort smooths | GAM RE | Section B |
| **B4** | four-way PGI × birth-year × region × gender | linear RE | Section C (C1) |
| **B5** | per-gender PGI→mobility slope over birth year (regions pooled) | GAM RE | per-gender mobility GAM |
| **B6** | Region × Gender PGI cohort smooths | GAM RE | Section C (C3) |
| **B0g** | overall Region × Gender PGI contrast (+ pre-reunification) | RE contrast | Section C (C0) |

There is no logistic specification check (A-spec/C2) in Section B — that was education-only. Where the
higher-order interactions (B4 four-way, B6 Region × Gender, B0g three-way) are null, the lower-order /
additive gender moderation is the reported effect.

## Files

| File | Role |
|---|---|
| `B_descriptives.R` | mobility descriptives over strata (sourced by `B_analysis.R`) |
| `B_analysis.R` | B0–B3 + B4/B5/B6/B0g; per-cell nonlinearity LRTs; LOO (incl. gender LOCO `B4_LOO_gender`, `B5_LOO`, `B5_LOO_diffsmooth`, §14); saves `B_Mobility_Models.rds` |
| `B_report.R` | builds `Supplementary_Data_5_Mobility_results.xlsx`, `B_Mobility_Plots.pdf`, LOO PDF, `B_Mobility_Overview.md` |

Gender leave-one-cohort-out (LOCO) for the **pooled** `PGI × gender` term
(`B4_LOO_gender`), the per-gender GAM (`B5_LOO`) and the per-fold female-minus-male
difference smooth (`B5_LOO_diffsmooth` / `B5_LOO_diff_intervals`) extend the
four-way `B4_LOO` — additional sensitivity (Plan_deviations §14). The
finding-matched random-effects alternative-PGI checks — including the
highest-priority B0 East–West mobility contrast (`S1_finding2`) and the B5
nonlinear reproduction (`S1_B5_nonlin`) — live in the sensitivities script (§15).

Run via `../run_B.command`. (Descriptives also live in the dedicated `../00_Descriptives` module, which
covers both education and mobility.)

## Outputs (under `runs/RUN_<TS>/B_Mobility/`)

- `Supplementary_Data_5_Mobility_results.xlsx` — sheets `00_Index`, `01_Linear_birth_year`,
  `02_Step_reunification`, `03_GAM_region`, `04_Overall_EastWest`, `05_Linear_gender`, `06_GAM_gender`,
  `07_GAM_region_gender`, `08_Overall_RxG` (human-readable term/column labels)
- `B_Mobility_Plots.pdf` (incl. B5 per-gender and B6 Region × Gender slope curves + diff smooths)
- `B_Mobility_LOO_Plots.pdf`
- `B_Mobility_Models.rds`, `B_Mobility_Overview.md`, `B_run.log`

## Sensitivities (Section B owns its mobility checks)

`B_sensitivities_analysis.R` → `B_Sensitivities_Models.rds` → `B_sensitivities_report.R`
(`Supplementary_Data_6_Mobility_sensitivities_results.xlsx` + `B_Sensitivities_Plots.pdf` +
`B_Sensitivities_Overview.md`):

| id | what | sample |
|---|---|---|
| **S1-B** | alternative PGIs (Cog / NonCog) on the B1 and B4 focals, BH per family | dedup-OLS (`dat_mob`) |
| **S2-B** | height-PGI negative control: → mobility (B1 spec), → height (sanity) — **new** | dedup-OLS |
| **S3-B** | heteroscedasticity (`gaulss`) per region + the B0 East-West contrast, gated on B focal significance | `dat_cluster_mob` (no RE) |

Migration diagnostics (S4) and the dedup/CR2/RE estimator comparison are **Section A only**.
Sensitivities stay dedup-OLS except the `gaulss` S3-B (full-family, no RE; `../../Plan_deviations.md`
§11/§12). Parental education is controlled throughout.

## Per-cohort appendix (mobility)

`B_appendix_analysis.R` → `B_Appendix_Models.rds` → `B_appendix_report.R`
(`Supplementary_Data_7_Mobility_appendix_results.xlsx` + `B_Appendix_Plots.pdf` +
`B_Appendix_Overview.md`):

| id | what |
|---|---|
| **E** | per-cohort B3-style mobility GAM (region-specific PGI cohort smooths, ParEdu-controlled) |
| **E-RxG** | per-cohort B6-style Region × Gender GAM on `mobility` (ParEdu-controlled) |
| **F-linear** | per-cohort robust linear gender model on mobility (PGI main / gender × PGI / Region × Gender × PGI) |

Appendix-grade (narrow per-cohort birth-year windows; some cohort × cell combos skipped). See
`../../Plan_deviations.md` §9/§13.

## Run order

`../run_B.command` runs all six stages under one `RUN_TS`: B_analysis → B_report →
B_sensitivities_analysis → B_sensitivities_report → B_appendix_analysis → B_appendix_report. Each stage
is independently runnable (the sensitivities read `B_Mobility_Models.rds` for their gates).
