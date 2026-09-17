# Sensitivity analyses

Three sensitivity analyses and the counting script behind Figure 1, each
reading the frozen analysis run and writing aggregate CSVs beside it. None
refits or overwrites the run, and none loads its model caches. Reference
estimates are read from the frozen export. Each sensitivity script stops if the
published estimate it anchors on fails to reproduce on the rebuilt analytic
sample.

| Script | What it computes | Writes |
|---|---|---|
| `scale_and_censoring.R` | Ceiling and floor shares; the two focal terms refitted on the raw outcome, on a rank-based inverse normal transform, and under an ordered probit with cluster-robust SEs. | `scale_censoring_ladder.csv`, `censoring_crosstab.csv`, `scale_censoring_checks.csv`, `scale_censoring_provenance.csv` |
| `ancestry_sensitivity.R` | Sample (SOEP-G vs pooled) crossed with the scope of the ancestry-PC residualisation (within-study vs across-study), for the attainment PGI × birth year × region term. | `ancestry_sensitivity.csv`, `ancestry_checks.csv`, `ancestry_provenance.csv` |
| `pc_moderator_refit.R` | The four published focal models refitted with ancestry-PC main effects, then with PC × moderator interactions, decomposing any movement into estimate and standard-error channels. | `pc_moderator_ladder.csv`, `pc_moderator_checks.csv`, `pc_moderator_provenance.csv` |
| `figure1_birthyear_counts.R` | Per-study × birth-year Ns for Figure 1. Fits no model. | `figure1_study_birthyear_counts.csv` |

Each script's header states its specification and its checks. Output is
aggregate only — coefficients, standard errors, intervals, *p*-values, counts
and shares. No individual-level row is printed or written.

## Environment

| Variable | Meaning |
|---|---|
| `RUN_TS` | timestamp of the run to read |
| `OUT_DIR` | where the CSVs go (default: the run's `manuscript_export/sensitivity_analyses/`) |
| `MANUSCRIPT_REPO` | checkout whose `01_results/` holds the frozen reference estimates |
| `DATA_ROOT` | data root, as elsewhere in this repository |
