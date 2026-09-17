# Analysis S5 — SNP-heritability of education

SNP-based heritability (h²_SNP) of years of education, estimated by
**Region × Reunification** (East/West × pre/post-1990) in the harmonized
pooled sample of four cohorts (BASE-II, SHIP, SOEP, TwinLife), per Analysis
S5 of `ANALYSIS_PLAN.md`. GCTA-GREML; education in raw years; covariates 20
cross-cohort PCs + birth-year (linear + quadratic) + gender, with a cohort
fixed effect in the primary specification. This page summarises the design;
everything operational lives in the technical reference. The estimates are
reported in the manuscript, not here.

> **More detail:**
> [`README_technical.md`](README_technical.md) — full pipeline, workflow
> commands, methods, diagnostics, output tree.
> [`../../Plan_deviations.md`](../../Plan_deviations.md) §7a–§7i —
> departures from the registered plan.

---

## Analysis hierarchy

- **Primary** — block-diagonal GREML with cohort fixed effects. Per-cohort
  GRMs are combined block-diagonally (cross-cohort entries zero), so V_G is
  identified from within-cohort genomic relatedness only while cohort-level
  mean differences in education are absorbed by a cohort indicator in
  `--covar`. This is the headline estimator.
- **Sensitivities** — (i) single pooled cross-cohort GRM at a liberal
  `--grm-cutoff 0.20`; (ii) per-cohort GREML at `--grm-cutoff 0.05`.
- **Secondary / exploratory** — R1 Region (East vs West), G1 recorded
  gender/sex, RG1 Region × gender/sex. Same block-diagonal design, reported
  as exploratory only.

## Power

Analytical Visscher-style power estimates specified in the analysis plan were
computed as diagnostics, but precision is interpreted using REML SEs because
the approximation is not appropriate for the final GREML specifications.

## Main limitations

- **SHIP dominance.** The East × pre-1990 cell (and the East secondary
  strata) are predominantly SHIP; these are **cohort-weighted estimands**,
  closest to a SHIP-Pomerania-weighted estimate, not balanced over all East
  Germans.
- **Post-1990 imprecision.** Both post-1990 cells are small and their
  estimates imprecise; a point estimate at the zero boundary is a constraint
  of the REML fit, not evidence of zero heritability.
- **Single variance components per cell.** The block-diagonal model fits one
  V_G and one V_e per cell across cohort blocks (cohort means are adjusted,
  cohort-specific variances are not). A diagnostic of within-cell
  education-variance ratios shows this simplification is not strongly
  contradicted by the data.
- **SOEP chr1 fallback.** The SOEP chr1 INFO file was never delivered;
  chr1 SOEP variants enter via a mildQC (R² > 0.1) fallback (see
  `manifest.tsv`).
- **No causal interpretation.** h²_SNP is a population variance-partition
  quantity in these samples; it is not a causal or transportable estimate.
