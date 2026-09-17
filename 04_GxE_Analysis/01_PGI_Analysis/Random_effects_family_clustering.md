# Random-effects handling of TwinLife family clustering

This document describes the **primary** method the Section-4 pipeline uses
to deal with the family structure in TwinLife: a mixed model with a
random intercept per family, fitted via `mgcv::gam(... + s(fid_re,
bs = "re"))` on the full-family analytic sample. The
`Supplementary_dedup` track in the same pipeline (`*_supplementary_dedup.pdf`)
fits the same focal models on a one-adult-per-family subsample with plain
OLS, as a conservative cross-check.

It is referenced from the analysis script header and from
`Plan_deviations.md` (section 5).

---

## 1. Why family clustering matters here

TwinLife is the only one of the four studies (BASE-II, SHIP, SOEP, TwinLife)
that recruits genetically and environmentally related individuals
together. A single TwinLife `fid` family typically contains:

- two parents (`ptyp` 110 / 120 / 300+),
- one or more twins (`ptyp` 1, 2),
- possibly additional non-twin siblings (`ptyp` 200, 201).

The other three studies are samples of unrelated individuals; for them
every row is effectively its own family.

Within a TwinLife family, both **genetic** and **environmental** sources
of correlation are present:

- Siblings share roughly 50 % of segregating variants on average;
  MZ-twin pairs share ~100 %.
- Parents do not share a recent common ancestor, but their PGIs are
  correlated through assortative mating (educational homogamy; the
  PGI-EA spousal correlation is ~0.18-0.45 in published estimates).
  They also share with their children the household, region, parental
  SES, and educational context that the analysis uses as a predictor
  or outcome.

If we treat rows from one such family as independent, the focal-coefficient
SE under-states the true uncertainty and the focal-coefficient *estimate*
is over-influenced by families that happen to have more sampled members.

---

## 2. The three clustering adjustments

Three adjustments for TwinLife families are available for each focal model:

| Option | Sample | Method | Adjusts point estimate? | Adjusts SE? |
|---|---|:---:|:---:|:---:|
| **(a) One-adult-per-family OLS** | dedup'd, one row per `fid` | `lm(...)` | yes, by construction | yes, by construction |
| **(b) Cluster-robust SE (CR2)** | full family | `lm(...)` + sandwich | no | yes |
| **(c) Random effects on `fid`** | full family | `mgcv::gam(... + s(fid_re, bs="re"))` | **yes** | **yes** |

Option **(a)** is statistically clean — every row is independent by
construction — but it discards roughly half of the TwinLife adult sample:
~700 rows out of the ~2 600 adult-with-PGI TwinLife rows.

Option **(b)** preserves all the rows but only re-derives the SE under a
within-family-correlated residual assumption (Liang-Zeger sandwich,
finite-sample CR2 correction from `clubSandwich`). The focal coefficient
*estimate* is unchanged from plain OLS on the bigger sample, so a family
with four sampled siblings still contributes four rows of *influence* to
the slope. Only the SE shrinks the family's effective contribution.

Option **(c)** also preserves all the rows, but additionally puts a
shared per-`fid` intercept into the model and shrinks within-family
deviations toward the family mean. Both the focal *estimate* and the
focal *SE* reflect the within-family correlation: a four-sibling
family contributes somewhere between 1 and 4 effective rows to *both*,
with the exact number determined by the estimated intra-class
correlation.

Option **(c)** is the primary method: it is the only one that down-weights
family-correlated information in both the slope and the SE. Option (a) is
reported as a *supplementary* sensitivity analysis (see
`*_supplementary_dedup.pdf`); option (b) appears in the Section A estimator
comparison.

---

## 3. How the random effect is constructed

The grouping factor `fid_re` is built once on the full analytic sample
`dat_cluster` (no twin/sibling exclusion, no family dedup, plus the
same harmonized data-quality filters as the supplementary sample):

```r
dat_cluster$fid_re <- factor(ifelse(
  dat_cluster$cohort == "TwinLife" & !is.na(dat_cluster$fid),
  paste0("TL_", dat_cluster$fid),
  "non_TL"
))
```

- **`"TL_<fid>"`** — one level per TwinLife family with a known
  family ID. Most contain 2–6 members.
- **`"non_TL"`** — a single shared level pooling **every** non-TwinLife
  row (BASE-II, SHIP, SOEP), plus the rare TwinLife rows that lack
  an `fid`.

Pooling all non-TwinLife rows into a single random-effects level is a
*deliberate simplification*. If we gave every non-TwinLife row its own
singleton family (~10 000 levels), the random-effects penalty matrix
would have ~10 000 rows and the fit would be very slow without
changing any focal coefficient. Concretely:

- Each `non_TL` member's residual contribution is independent of every
  other `non_TL` member's, so there is no within-`non_TL` variance to
  estimate as a random effect.
- The single shared `non_TL` BLUP can only absorb an *overall mean
  shift* between TwinLife and non-TwinLife rows. The focal slope
  (`PGI_Edu_z : BYc_z : east_west_c` and analogues) is invariant to
  any such constant shift, so the focal estimates are unaffected by the
  pooling.

For TwinLife families, each `TL_<fid>` level contributes the family's
shared deviation from the population intercept, weighted toward zero
by the random-effects penalty (BLUP shrinkage), with the amount of
shrinkage governed by the estimated `σ²_fid`.

---

## 4. The model formula

For each focal model, the random-effects term is appended to the
existing fixed-effect specification. For example, A1 (the linear three-
way interaction model from the analysis plan) becomes:

```r
m_A1 <- mgcv::gam(
  edu_z_kernel ~ (PGI_Edu_z + BYc_z + east_west_c + gender_c)^3
            + s(fid_re, bs = "re"),
  data = dat_cluster,
  method = "REML"
)
```

A3 (the cohort-varying spline model) is

```r
m_A3 <- mgcv::gam(
  edu_z_kernel ~ BYc + east_west_c + gender_c +
            PGI_W + PGI_E + PGI_W:BYc + PGI_E:BYc +
            s(BYc, by = PGI_W, k = 8, bs = "tp") +
            s(BYc, by = PGI_E, k = 8, bs = "tp") +
            s(fid_re, bs = "re"),
  data = dat_cluster,
  method = "REML"
)
```

The B-set, C-set, and the auxiliary A2 / B0 / B2 / C3-null specifications
follow the same pattern: the random-effects term is appended to whatever
fixed-effect formula was already on the model.

For the primary slope plots (`A1_slopes`, `A3_slopes`, `B1_slopes`,
`B3_slopes`, `C3_edu_slopes`, `C3_mob_slopes`, and their `_diff`
companions), the slope is read off as the *difference* in fixed-effect
predictions at PGI = 1 versus PGI = 0, holding every other predictor
constant. Because the random-effect contribution depends only on
`fid_re` (not on PGI), the random effects cancel in the slope
difference and the curves shown are the **fixed-effect partial slopes**,
correctly conditioned on the within-family correlation that was used
in fitting.

---

## 5. Variance-component diagnostics

The script reports the `mgcv::gam.vcomp()` output for the A1 fit in
`RESULTS$S7_A1_re_VarComp`:

```
component       std.dev   lower   upper
s(fid_re)        0.537    0.489   0.590     <- between-family SD on edu_z_kernel
scale            0.892    0.880   0.903     <- residual SD on edu_z_kernel
```

These imply

```
ICC = σ²_fid / (σ²_fid + σ²_resid)
    = 0.537² / (0.537² + 0.892²)
    ≈ 0.27
```

so **about 27 % of variance in standardized education is between
TwinLife families** in this analytic sample. That is a substantial
within-family correlation: a naive OLS that pooled all family members
without an adjustment would meaningfully over-state the precision of
PGI-related slopes.

---

## 6. Focal-coefficient comparison

For the A1 focal `PGI_Edu_z : BYc_z : east_west_c`, the three approaches
give:

| Method | N | est | SE | p |
|---|---:|---:|---:|---:|
| (a) Dedup OLS (supplementary) | 11 410 | -0.002 | 0.021 | 0.92 |
| (b) Cluster-robust SE | 13 357 | 0.003 | 0.019 | 0.90 |
| (c) **Primary: random effects** | 13 357 | 0.012 | 0.019 | 0.51 |

Both (b) and (c) yield essentially the same SE (the random-effects
penalty and the sandwich both encode the same within-family correlation
on this dataset). They differ on the point estimate because (c)
additionally down-weights large-family contributions to the slope. None
of the three rejects the A1 null at α = .05, so the qualitative
conclusion is invariant across the three methods.

---

## 7. Caveats / known limitations

- **TwinLife contributes only its parent generation.** Twin and non-twin
  sibling rows (`ptyp ∈ {1, 2, 200, 201}`) — the *children* in the
  families — are excluded from both tracks, so TwinLife enters as its
  parent / adult-relative generation. This is a design-scope choice: the
  adults span the historical cohorts of interest, and keeping the children
  would concentrate the youngest cohort almost entirely in TwinLife. It is
  *not* an age-completeness filter (that is the separate `birth_year ≤ 1996`
  cutoff) — ≈594 of the excluded twins/siblings are in fact complete-education
  genotyped adults (≈371 families, 1986–96 cohorts) who are set aside for
  the design reason above. Their within-pair relatedness would be absorbed
  by `s(fid_re)`, so re-including them is an available sensitivity check.
  See `../../Plan_deviations.md` §5.
- **Spousal pairs are not explicitly modelled.** Spouses share `fid`
  and therefore share a random intercept under our model. This is
  appropriate: their main correlated context is the household, which
  the random intercept absorbs. We do not model assortative mating in
  PGI separately.
- **Per-study exploratory analyses (Sections D, E, F) are run on the
  dedup'd sample only.** This is a known scope limit: the TwinLife panel
  would, in principle, also benefit from a within-study random-effects
  refit, but the per-study GAMs are exploratory appendix-style results,
  not focal hypothesis tests. They remain in the supplementary PDF.
- **`Sensitivity analyses S1-S5` (per the plan) are reported on the
  supplementary OLS sample.** This is also a known scope limit. Their
  multiplicity-corrected status is already governed by Benjamini-
  Hochberg within their own family; re-running them with random effects
  would not change their inferential structure.

---

## 8. References

- Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*
  (2nd ed.), §6.3 on simple random effects via `s(..., bs = "re")`.
- Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting
  linear mixed-effects models using lme4. *Journal of Statistical
  Software*, 67(1), 1-48 — for the partial-pooling interpretation
  underlying option (c).
- Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
  cluster-robust variance estimation. *Journal of Business & Economic
  Statistics*, 36(4), 672-683 — for the CR2 sandwich estimator, option (b).
