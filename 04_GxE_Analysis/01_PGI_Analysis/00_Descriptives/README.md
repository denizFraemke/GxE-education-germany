# 00_Descriptives — Phenotype descriptives

Descriptive workbook for the GxE analytic sample, produced by the same two-stage
pattern as the analysis sections (`*_analysis.R` → `*_Tables.rds` →
`*_report.R`), launched by `../run_Descriptives.command`.

## What it reports

For each measure — **Years of education** (raw), **Education (`edu_z_kernel`)** (the analysis
outcome), **Mobility** (`edu_z_kernel − parental_edu_z_kernel`), **PGI-Education (z)**, and
**Parental education (years)** — the workbook gives:

- **Overall** N and M(SD);
- **By region**: East and West N + M(SD), plus the **East − West** Welch t-test p and Hedges' g;
- **By gender**: Female and Male N + M(SD), plus the **Female − Male** Welch t-test p and Hedges' g;
- the same breakdown split into **three birth-year groups**: `< 1950`, `1950–1974`, `≥ 1975`
  (the 1975 boundary is the reunification cut — born ≥ 1975 = turned 15 in/after 1990 — shared with
  the A/B analyses).

## Method (see `../../../Plan_deviations.md` §5)

- **Sample.** **N and M(SD)** are reported on the **full-family RE sample** (`dat_cluster` /
  `dat_cluster_mob`) so they match the A/B analyses (TwinLife contributes all family adults). The
  **East−West / Female−Male tests** (Welch p + Hedges' g) are computed on the **dedup one-adult-per-family**
  subset (`dat_edu`) for clean independence. The coverage sheet reports both `N` (full-family) and
  `N_test` (dedup).
- **Sign conventions.** Region = East − West; Gender = Female − Male.
- **Test / effect size.** Welch two-sample t-test (unequal variances) for the p-value; **Hedges' g**
  (small-sample bias-corrected Cohen's d, pooled SD) for the effect. Computed only when both groups
  have ≥ `DESC_MIN_CELL_N` (= 30) observations in the dedup sample.
- **Graceful empties.** A **structurally absent** group renders **blank** — SHIP is East-only, so its West
  and East−West cells are blank; SHIP has no parental education, so its mobility/parental rows are blank.
  A **present-but-too-small** group (< 30 in the dedup sample) reads **`n/a`** for the test.
- **BASE-II** skips the by-birth-year table (Table 2): its birth-year span (~1927–1951) is too narrow for
  a meaningful split.
- **Terminology.** Each contributing dataset (BASE-II, SHIP, SOEP, TwinLife) is a **study**; "cohort"
  refers to birth cohort throughout.
- **Caveat.** These are unadjusted descriptive comparisons (no covariates, no multiplicity correction);
  raw years-of-education contrasts conflate the secular education trend, which is why the analyses use
  `edu_z_kernel`.

## Run

```
./run_Descriptives.command            # from 01_PGI_Analysis/
./run_Descriptives.command --dry-run  # print steps, no R / no SMB writes
```

## Outputs (under `runs/RUN_<TS>/00_Descriptives/`)

| File | Role |
|---|---|
| `Supplementary_Data_1_Descriptives_results.xlsx` | Workbook: `00_Index`, `01_Coverage`, `02_Overall_and_by_year`, one sheet per study (`03_BASE_II` … `06_TwinLife`), `07_Analysis_sample` |
| `Descriptives_Tables.rds` | Cached aggregate tables (no individual-level rows) |
| `Descriptives_Overview.md` | GitHub-readable summary + coverage table |
| `Descriptives_run.log` | Combined stage log |
