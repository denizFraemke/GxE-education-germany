# Plan deviations

This file catalogues where the implemented analysis in `04_GxE_Analysis/`
departs from `ANALYSIS_PLAN.md`. The plan stays as preregistered; deviations are
documented here so the differences are explicit and the reasoning is preserved.

Each entry states what the implementation does, the plan specification, the
rationale, and caveats where relevant.

Section numbers are stable anchors (referenced from code and other docs).
§1–§6, §8–§19 and §20–§21 cover `01_PGI_Analysis/` (the PGI gene–environment
analysis); §7 covers `02_SNP_Heritability/` (Analysis S5).

---

## 1. Sensitivity analyses — gating policy

**What.** S1 (alternative PGIs) and S2 (height negative control) run on every
analytic sample regardless of focal-test significance. S3 (heteroscedasticity)
and S4 (migration) run conditionally: only if at least one of the section's
focal tests is significant at α = .05 — the linear focal coefficient or the
smallest focal-GAM PGI-cohort smooth p-value (the per-smooth approximate
F-test, not the omnibus linear→nonlinear LRT). The gate is read from each
section's `*_Models.rds`:

- **A (attainment):** `gxe_a_sig` — A1's `PGI:BYc:east_west_c`, or the smallest
  A3 region-specific PGI smooth p-value.
- **B (mobility):** `gxe_b_sig` — B1's `PGI:BYc:east_west_c`, or the smallest
  B3 PGI smooth p-value. The B gate has two triggers only; the B2
  reunification-step focal is not read into it.

S4 is A-only: the migration-composition explanation applies to East–West
differences in PGI–attainment associations and is not informative for the
mobility focal. The per-focal dispersion robustness check (§12) is gated
analogously, per individual focal.

**Plan.** "Sensitivity analyses (run only if A1–A3 yields significant
gene–environment interaction parameters, e.g., β₁₄ in A1)" (`ANALYSIS_PLAN.md`,
line 205); the same conditional language applies to B.

**Rationale.** Keeps the conditional schedule plan-conform. S1 and S2 remain
unconditional because the BH-corrected alternative-PGI family (§8) and the
negative control are informative regardless of focal significance; S1's
multiplicity is controlled with Benjamini–Hochberg (q = .05) per the plan.

**Caveats.** Because the GAM trigger uses the per-smooth approximate F-test
rather than the omnibus LRT, a section can trigger on a significant
region-specific cohort smooth even when its omnibus nonlinearity test is
non-significant (this is the case for A: a West smooth is significant under the
confounded F-test while A3 omnibus is null — see §6). When focal p-values are
near .05 the trigger is sensitive to specification; each section's `S3_Status` /
`S4_Status` sheet records the exact p-values used. Within S4, the within-region
refit (S4d) is additionally gated on the three migration preconditions (§11).

---

## 2. Primary Educational attainment outcome is `edu_z_kernel` (within-birth-year kernel z)

**What.** Education-outcome models (A1/A2/A3 and the gender models A4/A5/A6/A0g)
regress on `edu_z_kernel` — education standardized within birth year via the
Gaussian-kernel construction of §3. This puts attainment on the same metric
as the mobility outcome (`mobility = edu_z_kernel − parental_edu_z_kernel`).

**Plan.** "all continuous variables, including covariates, will be
z-standardized" (`ANALYSIS_PLAN.md`, line 177); the within-cohort-z basis
(line 120) and its kernel realization are covered in §3.

**Rationale.** `edu_z_kernel` is the cohort-relative-standing construct the GxE
design targets, and it matches the metric of the mobility components. A steeper
PGI-on-relative-standing slope in one cohort is a substantive finding under the
study's hypothesis, not an artifact to be standardized away.

**Caveat.** A global-SD rescale `edu_std = education / EDU_SD` is still computed
and stored as an optional global-SD sensitivity column, but is not the focal
outcome. Because `edu_z_kernel` is undefined at the birth-year boundaries (§3), a
few boundary-cohort rows drop from the attainment models.

---

## 3. Mobility components use Gaussian kernel-weighted z-scoring

**What.** Child and parental education are z-scored with `kernel_z_score`
(Gaussian kernel, bandwidth = 5 years) over `birth_year` and `parent_by`;
`mobility = edu_z_kernel − parental_edu_z_kernel`. The child-education
`edu_z_kernel` this produces is also the attainment outcome (§2).

**Plan.** "we will z-standardize years of education within birth cohort
(subtracting the cohort-specific mean and dividing by the cohort-specific
standard deviation)" (`ANALYSIS_PLAN.md`, line 120).

**Rationale.** Kernel-z is a smooth-cohort version of within-cohort z: for each
birth year, local mean and SD are Gaussian-weighted moments over the birth-year
neighbourhood. Compared with hard 5-year bins (the literal plan), z-scores
transition smoothly across adjacent years and sparse cohorts borrow strength
from neighbours rather than hitting a small-N floor.

**Caveats.** The weighted variance uses the ML (denominator-N) form, which
slightly under-estimates SD (|z| over-estimated by ≤ 1 % at typical effective N).
At the oldest/youngest birth years the kernel window is one-sided; an
`n_eff ≥ 10` floor guards degenerate cases but cannot correct the boundary tilt
(≈ first/last 3–4 birth years at bandwidth 5).

---

## 4. `|PGI_Edu| > 3` outlier exclusion

**What.** Rows whose absolute global-z PGI-Education exceeds 3 are dropped.

**Plan.** Not specified.

**Rationale.** Conventional QC for harmonized PGI distributions to prevent
extreme-PGI outliers from dominating regression diagnostics, especially for
higher-order interactions. Applied after cross-cohort PGI harmonization
(z-scoring on the pooled sample) and symmetric.

---

## 5. TwinLife family handling (primary: random effects; supplementary: one-adult-per-family)

**What.** TwinLife twins and non-twin siblings (`ptyp ∈ {1, 2, 200, 201}`) are
dropped from both analytic samples, so TwinLife contributes only its
parent/adult-relative generation. For the remaining adults, every focal model —
including the gender models (§10) — is run on two tracks:

- **Primary:** the full-family sample (no `fid` deduplication) with a per-family
  random intercept, `mgcv::bam(... + s(fid_re, bs = "re"))`. The focal estimate
  and SE are partial-pooled toward the within-family mean. Method detail in
  `04_GxE_Analysis/01_PGI_Analysis/Random_effects_family_clustering.md`.
- **Supplementary (dedup):** one-adult-per-family subsample (seeded draw within
  each `fid`), plain OLS / mgcv without random effects — the conservative
  cross-check for row-level independence.

**Plan.** The preregistration does not specify a TwinLife family-clustering
method.

**Rationale.** *Twins/siblings dropped:* the analytic interest is the adult /
parent generation, which spans the full historical cohort range (including the
pre-reunification East/West cohorts of interest); keeping the TwinLife children
would confound the youngest cohort with a TwinLife study/genotyping effect and
introduce genetically related pairs. (Education completeness is enforced
study-wide by the separate `birth_year ≤ 1996` filter (§20), not by this
exclusion; the ≈594 complete-education adult twins it sets aside remain
re-includable under the `s(fid_re)` intercept.) *Random effects:* the intercept
absorbs the within-family outcome correlation from shared household / region /
parental-SES context and assortative mating; on the current sample the education
ICC ≈ 0.27, so the penalty does substantive work. The dedup track is the
conservative cross-check; a three-way comparison with cluster-robust (CR2) SEs is
in `Random_effects_family_clustering.md` (and §11).

---

## 6. Per-cell nonlinearity decided by a clean nested LRT

**What.** For every varying-coefficient GAM (A3/A5/A6, B3/B5/B6) the
nonlinearity decision for each cell (region, gender, or region × gender) is
a clean nested likelihood-ratio test: add *only that cell's* PGI × birth-year
smooth to the linear model and test linear vs nonlinear. This per-cell LRT
drives (i) the "nonlinearity conclusion" in each GAM-overview table and (ii) the
display curves, which draw the nonlinear smooth only where the per-cell LRT is
significant and the linear fit otherwise.

**Plan.** The plan specifies the varying-coefficient GAM but not the per-cell
nonlinearity test.

**Rationale.** The `summary.gam` per-smooth F-test for `s(BYc, by = PGI_cell)`
is confounded with the parametric `PGI:BYc` slope (it tests the whole smooth,
linear part included) and overstates curvature; the omnibus linear→nonlinear LRT
can be dragged either way by other cells. The per-cell nested LRT isolates each
cell's departure from linearity. Under this convention the attainment cohort
effect is linear (A3 both regions, A5 both genders), with the only curvature a
West-male bend (A6); for mobility the nonlinear signal is West-male-concentrated
(B5, B6), and the West B3 bend is fragile under leave-one-cohort-out (§9).

---

## 7. SNP-heritability (S5) — implementation choices

These deviations apply to `04_GxE_Analysis/02_SNP_Heritability/`, which
implements Analysis S5 (multi-component GREML on a Region × Time stratification,
with heterogeneity tests across cells).

### 7a. Power calc uses a scripted R implementation of Visscher 2014

**What.** `scripts/08_power_calc.R` computes
`SE(h²) = √(2 / (N² · Var(off-diag-GRM)))` per stratum over an h² grid.
`Var(off-diag-GRM)` is the empirical variance of each stratum's `*_unrel` GRM
off-diagonals, computed in a single streaming pass.

**Plan.** "we will conduct a power analysis using the GCTA-GREML Power
calculator (Visscher et al., 2014)" (`ANALYSIS_PLAN.md`, line 247).

**Rationale.** The GCTA web calculator is not reachable from the compute nodes
used for this analysis; a scripted implementation of the same formula reproduces
the SE/power numbers, is auditable and reproducible, and a grid output (one row
per stratum × h²) is more useful for sensitivity reporting than a single-point
query.

**Caveat.** The Visscher `SE(h²)` is leading-order independent of the true h² and
treats the sample as fully unrelated; the `--grm-cutoff` filter satisfies this
approximately.

### 7b. Cross-cell h²_SNP test: Cochran's Q heterogeneity (not a joint mGRM LRT)

**What.** Heterogeneity across the K cells is tested directly from per-stratum
standalone REML fits (`scripts/07_lrt.R`):

- **Omnibus Cochran's Q** on the K cell h² estimates (the primary test of plan
  H₀: all h²_SNP equal): `Q = Σ_k (h²_k − h²_pooled)² / SE(h²_k)²`,
  `h²_pooled = Σ w_k h²_k / Σ w_k`, `w_k = 1/SE(h²_k)²`, `Q ~ χ²_{K−1}`.
- **Pairwise Z-tests:** `Z_ij = (h²_i − h²_j) / √(SE_i² + SE_j²)`.
- **I²** heterogeneity (Higgins & Thompson 2002).

Output: `output/reml/heterogeneity.tsv`.

**Plan.** "H₀: h²_SNP,East,pre = h²_SNP,East,post = h²_SNP,West,pre =
h²_SNP,West,post" (`ANALYSIS_PLAN.md`, line 249).

**Rationale.** `gcta --reml --mgrm` requires overlapping individuals across GRMs
and does not support the disjoint-sample cell structure here, so the
plan-literal joint-fit LRT is unavailable. Cochran's Q on
inverse-variance-weighted ratio estimates is the standard meta-analytic
heterogeneity test and a direct ratio-level analogue of the variance-component
LRT — in fact a closer match to the plan's H₀ (equality of h²_SNP, a variance
ratio) than a σ²_G-equality LRT under shared σ²_e would be.

**Caveats.** Cochran's Q is asymptotic and loses power when SEs are large
relative to between-cell differences (the imprecise post-1990 cells contribute
little); the omnibus Q does not identify which cells differ — the pairwise
Z-tests fill that gap.

### 7c. No 6-cell temporal sensitivity

**What.** Only the preregistered 4-cell (Region × pre/post-1990) stratification
is fit; a finer 6-cell split is not produced.

**Plan.** The plan specifies 4 cells (line 249).

**Rationale.** A 6-cell split leaves two cells below the N floor for stable REML;
the substantive question (h² heterogeneity across Region × Time) is answered by
the 4-cell main plus Cochran's Q (§7b).

### 7d. HWE filter is not applied at the GRM-construction step

**What.** `scripts/03_build_stratum_grms.sh` builds GRMs with
`--maf 0.01 --autosome` and no HWE filter.

**Plan.** "harmonized intersecting autosomal SNP set … (imputation INFO ≥ .90,
MAF ≥ .01, HWE p > 1×10⁻⁴)" (`ANALYSIS_PLAN.md`, line 247).

**Rationale.** GCTA's GRM routine exposes no `--hwe` flag, and pooled-sample HWE
p-values are confounded with cross-cohort allele-frequency structure (East vs
West, genotyping batch, imputation panel) — filtering by pooled HWE would
preferentially remove ancestry-informative variants. INFO ≥ 0.9 plus the
cross-cohort SNP intersection already pre-screen for genotyping quality.

**Caveat.** Per-stratum HWE filtering would break the shared-SNP assumption of
the multi-component GREML (different strata → different SNP sets); not pursued.

### 7e. Unrelated-filter threshold and the cell-level estimator

**What.** Two cell-level tracks are produced:

- **Block-diagonal (primary).** Per-(cohort × cell) GRMs built from per-cohort
  GRMs (within-cohort allele frequencies, `--grm-cutoff 0.05` within cohort),
  combined block-diagonally so cross-cohort pairs contribute zero to V_G
  (steps 03c/03d/04d/05d/06d).
- **Pooled cross-cohort GRM at `--grm-cutoff 0.20` (sensitivity).**

**Plan.** "filtered to be approximately unrelated (pairwise genetic relatedness
< 0.025)" (`ANALYSIS_PLAN.md`, line 247).

**Rationale.** The pooled cross-cohort GRM carries ancestry / imputation-panel
structure in the 0.025–0.10 off-diagonal band that GCTA misreads as relatedness;
at the plan's 0.025 cutoff this wipes out the minority cohort within each cell,
breaking the intended East/West composition. The block-diagonal design
removes cross-cohort structure by construction (within-cohort AF
standardization), so the conventional 0.05 cutoff suffices *and* cohort
composition is preserved — it is the primary cell-level estimator. The pooled
track at 0.20 is retained as a liberal-kinship sensitivity (a documented
relaxation past the 0.05 convention, accepted in exchange for cohort-correct
cell composition).

**Caveats.** The 0.20 pooled cutoff admits up to ~2nd-cousin pairs, biasing V_G
upward; the block-diagonal design identifies V_G from within-cohort relatedness
only (cross-cohort A = 0 by construction). The per-stratum REML SE is the
primary SE; the Visscher SE (§7a) is a power-calc input.

### 7f. Birth-year added to the GREML qcovar (linear + quadratic)

**What.** The qcovar carries `crossPC1..20`, `birth_year_c`, and
`birth_year_c_sq` (centred on the analytic-sample mean); `gender` is in
`--covar`. GCTA partials these out before partitioning V_P.

**Plan.** S5 (line 247) lists GRM filters and the unrelated cutoff but does not
enumerate fixed-effect covariates beyond ancestry.

**Rationale.** The cells span wide birth-year windows (the pre-1990 cells ≈ 55
years) over which years of education trends strongly; without a within-cell
birth-year adjustment that between-cohort variance loads onto V_e and biases h²
downward. A linear-plus-quadratic term lets each stratum absorb its secular
trend before variance partitioning (standard GREML setup). Centring is for
numerical conditioning (raw BY and BY² are collinear at r ≈ 0.99).

**Caveat.** h²_SNP is invariant to additive linear partialling-out, so the ratio
is well-defined either way; what changes is the per-cell V_G / V_e recalibration.

### 7g. Cohort fixed effects in the block-diagonal cell-level track

**What.** The block-diagonal track (step 05d) adds a categorical `cohort_int` to
`--covar` (`data/gender_cohort.covar`); the pooled and per-cohort tracks use
plain `gender.covar`.

**Plan.** Not specified — an implementation choice driven by the block-diagonal
design.

**Rationale.** Under the block-diagonal design cross-cohort relatedness entries
are zero by construction, so cohort-level phenotype-mean differences cannot be
absorbed by V_G; without an explicit cohort fixed effect they load into V_e and
deflate h². Including the cohort categorical removes those mean differences
before variance partitioning. With it, the two informative pre-1990 cells give
h² ≈ 0.28 (East-pre 0.282 ± 0.050, West-pre 0.281 ± 0.172), converging with the
per-cohort SHIP estimate (≈ 0.26) and the pooled-GRM track.

**Status of the three tracks.** Block-diagonal with cohort fixed effects (05d)
is the primary cell-level estimator; the pooled cross-cohort GRM at 0.20
(05) is a sensitivity (stronger assumptions); the per-cohort fits (05c)
answer "is any single cohort driving the pooled estimate?" (none is; SHIP
dominates precision).

**Caveat.** The post-1990 cells are imprecise under every specification
(REML SE ≥ 0.35); their near-zero point estimates are boundary values, not
evidence of zero heritability. Convergence across tracks in the pre-1990 cells
supports a calibrated reading; equality across cells is not established, only
consistent with the data given the precision.

### 7h. Secondary / exploratory comparison tracks (R1, G1, RG1)

**What.** Three additional block-diagonal GREML tracks compare h² across simpler
groupings of the same sample: R1 Region (East vs West), G1 recorded
gender/sex (Male vs Female), RG1 Region × recorded gender/sex. Each uses the
20 cross-cohort PCs + birth_year_c + birth_year_c_sq in `--qcovar` and the
remaining categorical factors + `cohort_int` in `--covar`. These are
secondary / exploratory — they surface group-level heterogeneity patterns,
not calibrated per-group population estimates.

**GRM order — subset-then-prune.** The secondary tracks subset the pre-cutoff
per-cohort GRM to each (cohort × stratum) sub-block, then apply
`--grm-cutoff 0.05` within the block (`scripts/03e_…`). This differs from the
primary track's prune-then-subset order: tracks like G1 split close relatives
across strata (a cross-sex sib pair lands in different strata), and
subset-then-prune keeps both members within their respective blocks. A
regression check (`regression_check_blockdiag_primary.sh`) verifies the primary
cell h²/SE are unchanged before any secondary track runs.

**Caveats.** Recorded gender/sex labels are not measures of gender identity or
biological sex. Sample sizes are unbalanced across cohorts; strata where one
cohort contributes > 80 % are flagged as cohort-weighted rather than
population-representative. Non-significant heterogeneity is reported as "no
evidence of heterogeneity", not "evidence of equality".

### 7i. V_G across strata is reported descriptively; no formal V_G heterogeneity test

**What.** V_G, V_G_SE, V_e, V_e_SE and V_P are reported per stratum in every
per-stratum TSV and summary, and the README surfaces V_G/V_P alongside h² for
the four primary cells. A formal Cochran's Q test on V_G is not reported.

**Plan.** "we will report the implied additive SNP variance V_G = h²_SNP × V_P …
Comparing V_G across strata will clarify whether differences in h²_SNP reflect
changes in genetic variance versus phenotypic variance" (`ANALYSIS_PLAN.md`,
line 254).

**Rationale.** The plan asks for a descriptive comparison, which is implemented.
The two post-1990 cells have V_G_SE ≥ 2.6 with V_G near zero, so a Q test on V_G
would be dominated by their noise and carry essentially no power; the
point-estimate comparison is the informative reading. (Substantively: in the
pre-1990 cells h² is near-identical, 0.282 vs 0.281, while V_G and V_P are both
≈ 27–28 % higher in West-pre — the two shifts offset; suggestive but not
statistically distinguishable.)

---

## 8. wf-PGI-Education omitted from S1

**What.** Within-family PGI-Education (Tan et al. 2024) is not included in S1.
The alternative-PGI Benjamini–Hochberg family per focal analysis is three
tests, not four: PGI-Education (reference), PGI-Cognition, PGI-Non-Cognitive.

**Plan.** "substituting the primary predictor (PGI-Education) with … PGI-Cognitive,
PGI-Non-Cognitive, and wf-PGI-Education … a family of four tests per analysis …
Benjamini–Hochberg at q = .05" (`ANALYSIS_PLAN.md`, line 209).

**Rationale.** wf-PGI summary statistics are not available for the cohorts
scored here, and carrying empty wf-PGI columns would imply a test that was not
performed. The BH threshold (q = .05) is unchanged.

**Caveat.** The remaining family is correlated, so BH on three correlated tests
is conservative; the family can be restored to four with no change to the
procedure if wf-PGI weights become available.

---

## 9. Section B parental-education control at focal-cell flexibility

**What.** The mobility models control parental education
(`parental_edu_z_kernel`) at the same functional flexibility as the focal PGI:

- **B1 / B2** implement the preregistered Eq-(4) ParEdu interaction set
  (`ParEdu` main + `ParEdu × {Region, PGI, BY}` and the three-way terms; B2
  substitutes the reunification step for BY). This is plan conformance, not a
  deviation.
- **B3 (GAM)** adds, beyond Eq-(5)'s flat covariate, a region-specific smooth
  of ParEdu over birth year `s(BYc, by = ParEdu_W) + s(BYc, by = ParEdu_E)`,
  paralleling the region-specific PGI smooths.
- **The gender mobility models** carry ParEdu at the matching cell flexibility:
  B4 a four-way `ParEdu × BY × Region × Gender`; B5 gender-specific ParEdu
  smooths; B6 Region × Gender ParEdu cell smooths.

**Plan.** Eq (4) (B1/B2) specifies the full ParEdu interaction set; Eq (5) (B3)
specifies ParEdu only as a flat linear covariate (`+ β₅·ParEdu_z`).

**Rationale.** Parental education has a strong nonlinear, cohort-varying
effect on mobility; entered only as a flat covariate it is unmodelled and
confounds the focal PGI cohort-smooth — the B3 focal nonlinearity test is
null under the flat covariate but significant once ParEdu is controlled at the
same flexibility (a West-specific PGI × cohort nonlinearity). Valid inference on
the focal PGI nonlinearity requires controlling ParEdu at matched flexibility;
the same logic extends to the gender mobility GAMs.

---

## 10. Gender analyses estimated within Sections A and B

**What.** The preregistered Section C gender analyses are reported as gender
moderators within Sections A and B (alongside the cohort and region moderators
of the same outcome) rather than as a standalone section, estimated on the
primary random-effects samples (§5):

- **A4 / B4** — linear four-way `PGI × BY × Region × Gender`.
- **A6 / B6** — Region × Gender cohort-varying GAM.
- **A0g / B0g** — overall Region × Gender PGI contrast (no BY × PGI term),
  full and pre-reunification.
- **A5 / B5** — per-gender cohort GAM, regions pooled (the gender analog of
  A3/B3).
- **A-spec** — the gender-specific logistic-transition model (plan Eqs. 8–10),
  reported only as a specification check: it does not converge / does not
  improve on the linear model, so its convergence status, AIC/BIC and boundary
  LRT are reported but it is not interpreted as a focal result.

**Plan.** Section C (Eqs. 6–10): gender-interaction models on attainment and
mobility.

**Rationale.** Reporting the gender analyses by phenotype rather than in a
separate section is a structural choice; the models are the preregistered
Section C analyses estimated on the §5 primary track. The logistic transition is
reported only as a specification check because it does not converge in practice.

---

## 11. Cross-section sensitivities (S1 / S2 / S4 / Method)

**What.** Each phenotype section owns its sensitivities
(`A_sensitivities_*`, `B_sensitivities_*`):

- **S1 — alternative PGIs.** PGI-Cog and PGI-NonCog substituted for PGI-Edu on
  the linear focal and the four-way gender focal (A1/A4, B1/B4), Benjamini–
  Hochberg per family (wf-PGI omitted, §8). Fit on the dedup-OLS samples with
  `lm` — a different estimator than the RE primary, but the relative comparison
  across PGIs (what S1 is for) is unaffected.
- **S2 — height negative control.** The height PGI substituted into the linear
  specification on education (S2-A) and on mobility (S2-B), plus a
  PGI-Height → height sanity check.
- **S4 — migration (A-only).** Two gates: outer (an A focal is significant) and
  inner (the within-region refit S4d runs only if all three migration
  preconditions hold — lower East mean PGI, lower East PGI variance, equal
  education variance).
- **Method comparison (A-only).** The A1 focal under dedup-OLS, cluster-robust
  CR2 SEs (`clubSandwich`, clustered by family), and RE (read from the A model
  cache).

**Plan.** S1/S2 (line 209), the S4 migration diagnostic, and the
family-clustering robustness checks.

**Rationale.** S1/S2 on the dedup samples match the relative-comparison purpose
and the conservative-independence cross-check; the CR2 row adds a third
family-clustering estimator so cross-estimator agreement is the substantive
read.

**Caveat.** The CR2 degrees of freedom for the high-order interaction can be low
(few effective family clusters carry it), making the CR2 test conservative.

---

## 12. Heteroscedasticity (S3) and per-focal dispersion robustness

**What.** Two related dispersion checks, both gated on focal significance (§1)
and both fit on the full-family samples (`dat_cluster` / `dat_cluster_mob`),
without the family random effect (the `gaulss` location-scale family cannot
carry it) — not on the dedup subsample:

- **S3-A / S3-B — heteroscedasticity.** Per region (East/West), Gaussian
  location-scale `gaulss()` models share a mean spec and vary the dispersion
  predictor (`σ ~ BYc`, `σ ~ PGI`, and a variance-aware `σ ~ BYc+PGI+BYc:PGI`
  when the focal smooth is significant). A-S3 uses `edu_z_kernel`; B-S3 uses
  mobility with the region-specific ParEdu smooth (§9). For B, the overall
  East–West PGI-mobility contrast (B0) is additionally tested OLS vs gaulss
  (full and pre-reunification).
- **Per-focal dispersion robustness.** `focal_dispersion_check()` re-estimates
  each focal that is significant in the main RE analysis with the residual
  variance free to follow the same structure as the focal (`scale ~ 1` vs
  `scale ~ <focal structure>`), reporting the focal estimate/p homoscedastic vs
  heteroscedastic and the scale LRT. Non-significant focals are listed but not
  fitted.

**Plan.** Sensitivity analysis for variance dispersion (run conditionally, §1).

**Rationale.** A significant genetic-dispersion slope or variance-aware GxE
dispersion term would indicate the focal GxE is partly a variance phenomenon,
not only a mean shift. On the current data the focal results are not
dispersion artifacts: B's findings survive (B3 West smooth, B0 contrast), and
every significant focal — including the gender effects in both sections —
remains significant when its variance is modelled (the gender effects, if
anything, strengthen, because the homoscedastic fit was mildly biased toward
zero by the variance heterogeneity). Strong genetic heteroscedasticity
(education variance rises with PGI) coexists with these mean effects without
explaining them.

**Caveat.** Because these fits carry no family random effect, the main RE
p-value is used only to gate (not refit); and the per-region S3 fits pool the
cohort smooth across cells rather than reproducing the full cell structure, so
they are approximations.

---

## 13. Per-study appendix (within-dataset)

**What.** Appendix-grade per-study analyses fit *within each dataset* (BASE-II /
SHIP / SOEP / TwinLife), in each section's `*_appendix_*`:

- **D / F (attainment, `A_appendix`)** — per-study region GAM and Region ×
  Gender GAM on `edu_z_kernel`.
- **E / E-RxG (mobility, `B_appendix`)** — per-study region GAM and Region ×
  Gender GAM on mobility (ParEdu-controlled).
- **F-linear** — per-study robust linear gender model
  (`outcome ~ PGI × Region × Gender + BY`), paper-comparable.

**Plan.** Not in the preregistered tests.

**Rationale.** Single-study complements to the pooled A/B analyses. Per-study
birth-year windows are narrow, so the GAMs are appendix-grade and some
study × cell combinations are skipped for lack of N (D/E: N ≥ 100 and both
regions ≥ 30, else single-region; F: each cell ≥ 50, else a single-region
2-cell fit, else skip; SHIP is East-only and absent from mobility). The pooled
A/B analyses plus the leave-one-cohort-out checks (§14) remain the
better-powered heterogeneity check.

**Note.** The SOEP F-linear row is the direct comparison to Spörlein, Zoch &
Schlueter (2025, *RSSM*), whose Gene-SOEP sample is the SOEP study analysed
here: men show a steeper PGI slope than women (gender × PGI ≈ −0.13, p ≈ .016 on
`edu_z_kernel`), gap descriptively larger in the East.

---

## 14. Gender findings — finding-matched leave-one-cohort-out (LOCO)

**What.** The gender results carry finding-matched leave-one-cohort-out checks
beyond the four-way `A4_LOO` / `B4_LOO` (which track only the
`PGI_Edu_z:BYc_z:east_west_c:gender_c` term). Cached alongside those LOO tables:

- **`A4_LOO_gender` / `B4_LOO_gender`** — the same primary A4/B4 models refit per
  fold (identical sample, family random intercept, covariates, centring,
  higher-order structure; four-way ParEdu control preserved for B4), extracting
  the pooled `PGI_Edu_z:gender_c` interaction with estimate, SE, 95% CI,
  exact *p*, supported birth-year range, change-from-full, and convergence.
- **`A5_LOO` / `B5_LOO`** — per fold, the per-gender GAM omnibus
  linear-vs-flexible LRT plus female- and male-specific curvature LRTs.
- **`B5_LOO_diffsmooth` / `B5_LOO_diff_intervals`** — per fold, the
  female-minus-male PGI-slope difference smooth over birth year with pointwise
  and simultaneous 95% bands (same `diff_smooth_pooled_gender_simul` procedure as
  the primary B5 export, `N_SIM` draws), trimmed to the fold's supported
  birth-year range, and the contiguous birth-year intervals where the
  simultaneous band excludes zero.

"LOCO" (leave-one-cohort-out) is used in code and docs to distinguish these
model checks from the leave-one-out (LOO) GWAS-weight construction upstream.

**Plan.** The preregistration specifies no leave-one-cohort-out analysis for any
finding (it is absent from `ANALYSIS_PLAN.md`).

**Rationale.** The four-way LOO does not speak to the reported pooled
`PGI × gender` interactions (§10, A4/B4) or to the nonlinear cohort-varying
gender pattern in mobility (B5). These finding-matched refits show whether each
gender result is carried or masked by any single cohort. Stability is read from
the direction, magnitude and uncertainty of the coefficient / test across folds,
not from each reduced model staying significant — omitting a cohort reduces
N and birth-year coverage.

**Status.** Additional sensitivity analysis (not preregistered).

**Caveats.** The fold that holds out TwinLife drops the family random intercept
(the only cohort carrying real family ids), matching the LOO engine used for the
four-way term. Folds that shrink birth-year coverage can widen or destabilise
the smooths; the supported birth-year range is reported per fold so trimmed
ranges are explicit.

---

## 15. Alternative-PGI substitution — finding-matched random-effects extension

**What.** The preregistered alternative-PGI substitution (§8; `ANALYSIS_PLAN.md`
line 209) is implemented in two layers. S1 (deduplicated OLS,
`S1_AltPGIs`) substitutes PGI-Cognition and PGI-Noncognitive into the A1/B1
three-way and A4/B4 four-way focal terms. Alongside it — not in place of it —
a finding-matched random-effects layer matches the actual focal findings, fit
with the primary `mgcv::bam(..., fREML) + s(fid_re, bs = "re")` specification on
a common complete-case sample across the three PGIs (PGI-Education,
PGI-Cognition, PGI-Noncognitive), each PGI re-residualised on `crossPC1–10` and
re-z-scored within that common sample, with PGI-Education refit on the same
sample:

- **`S1_finding1`** — attainment across cohorts: pooled `PGI × BY` (the reported
  increase) plus the contextual `PGI × BY × region` (A1-equivalent).
- **`S1_finding2`** — the East–West mobility contrast
  `PGI × region` with East/West standardised slopes (B0-equivalent).
- **`S1_finding3_att` / `S1_finding3_mob`** — the pooled `PGI × gender`
  interaction with gender-specific slopes (A4/B4-equivalent).
- **`S1_A5_nonlin` / `S1_B5_nonlin`** — the per-gender GAM omnibus + curvature
  tests; B5 additionally the female-minus-male difference smooth by the same
  procedure as the primary analysis.

**Multiplicity.** Benjamini–Hochberg is applied within each focal hypothesis
across the three PGIs (e.g. the three `PGI × region` mobility tests form one
family; the three `PGI × gender` attainment tests another). Unrelated
coefficients and model comparisons are not pooled into one correction
family. Both raw and adjusted *p* are retained; the families are recorded in the
output (`S1_FindingMatched_Status`).

**Plan.** "re-estimate Analyses A1–A3 substituting … PGI-Cognitive,
PGI-Non-Cognitive … Benjamini–Hochberg at q = .05" (`ANALYSIS_PLAN.md` line 209);
the substitution is prespecified, its BH family is three tests after §8.

**Rationale.** The prespecified S1 targets the A1/B1 three-way and A4/B4 four-way
terms, which are not the manuscript's focal findings (the pooled `PGI × BY`
increase, the B0 East–West contrast, the pooled `PGI × gender` interactions, the
B5 nonlinear gender pattern). Matching the alternative-PGI comparison to those
findings, on the primary RE specification and a common sample, is required for a
like-for-like read; the common sample and re-standardisation ensure differences
across PGIs are not confounded by sample composition or scaling. Differences
across PGIs are reported as differences in polygenic association, not as
changes in genetic effects.

**Status.** Prespecified substitution, with the finding-matched RE
implementation an extension of S1; it does not replace the dedup-OLS
`S1_AltPGIs`.

**Caveats.** The three PGIs are correlated, so BH across three is conservative.

The finding-matched refit re-residualises each raw PGI on `crossPC1–10` and
re-z-scores it within the analytic sample — the identical operation
`prepare_samples.R` applies per sample (`df`, `dat_cluster`,
`dat_cluster_mob`). It therefore reproduces the primary A1/B0/A4/B4 estimates by
construction: the finding-matched PGI-Education B0 `PGI × region` = 0.0727,
equal to the primary B0 (0.0727), and the finding-matched PGI-Education
`PGI × gender` equals the primary A4/B4 gender terms (−0.039 attainment,
−0.066 mobility). The LOCO full-sample rows (§14) provide the same cross-check
for the linear focal terms.

This equality holds only when the finding-matched and primary estimates are
built from the same PGI residualisation. Estimates on this scale are sensitive
to that step: a PGI that is not re-residualised and re-standardised within the
mobility sample is not unit-SD there, and the resulting focal estimate shifts
materially. Any refit or import must therefore share one residualisation with
the primary models.

---

## 16. One EA4 weight file used across all cohorts

**What.** BASE-II, SHIP-START, SHIP-TREND, SOEP-G, and TwinLife are all scored
with the SBayesR posterior derived from `EA4_excl_SHIP`. Thus, SHIP is excluded
from the discovery GWAS used to score it, whereas BASE-II is not separately
removed from the discovery GWAS. The same weight file and the same raw PLINK2
`BETA_SUM` scale are used in every cohort.

**Plan.** The plan specified discovery statistics excluding SHIP for SHIP and
discovery statistics excluding BASE-II for BASE-II, to avoid sample overlap in
both cohorts.

**Rationale.** BASE-II contributed approximately 0.07% of the EA4 discovery
sample, corresponding to an expected in-sample inflation of approximately
1.0007 in explained variance. This overlap effect is too small to distinguish
statistically in the present application. A common posterior also keeps all five
cohorts on one scoring scale, so cross-cohort raw-score differences are not
created by using different leave-one-out posteriors. A BASE-II-specific
`EA4_excl_BASEII` posterior is not used because its analytical calibration
factor does not map cleanly onto PLINK2's `BETA_SUM` scale, so a second
posterior would require an additional calibration step that the single-weight
approach avoids at no scientific cost at the precision relevant to this study.

**Caveat.** BASE-II is not strictly independent of its discovery GWAS. The
decision rests on the very small expected overlap bias and is therefore reported
transparently rather than described as complete leave-one-cohort-out scoring.

---

## 17. Monte Carlo power analysis for the PGI interaction models not conducted

**What.** The planned Monte Carlo power analysis for the PGI interaction models
was not conducted. This does not concern the separate, implemented power
calculation for the GREML analyses (§7a).

**Plan.** Before the main analyses, the plan proposed permuting educational
attainment and simulating a range of true effects in the observed design to
contextualize potential null interaction results (`ANALYSIS_PLAN.md`, lines
141–143).

**Rationale.** Once the final model specifications and analytic samples were
fixed, the fitted estimates, standard errors, and confidence intervals directly
characterized the precision of the interaction estimates. A retrospective power
calculation tied to the observed effect adds no independent information, because
it is a transformation of the same sampling uncertainty. Confidence intervals
and robustness analyses therefore describe what effect sizes remain compatible
with the data.

**Caveat.** A non-significant interaction is not evidence of equivalence. A
design-based simulation of minimum detectable effects under prespecified,
substantively meaningful interaction sizes remains a possible supplementary
precision analysis; it would not constitute post hoc support for a null
conclusion.

---

## 18. PGIs residualised on pooled cross-cohort ancestry PCs

**What.** Within each relevant pooled analytic sample, each PGI was residualised
by ordinary least squares on `crossPC1`–`crossPC10` and the residuals were then
standardised in that sample. The educational attainment and educational
mobility samples were processed separately. Genotype batch was not included
because no common batch variable was available in the harmonised data.

**Plan.** The plan specified residualising each PGI for the first 10
Germany-specific ancestry PCs and genotype batch using random effects within
cohorts.

**Rationale.** The cross-cohort PCs were the common ancestry coordinates
available for all contributing datasets. Pooled residualisation removes the
linear association between each PGI and these shared ancestry axes while
preserving a common cross-cohort score scale; processing each analytic sample
separately ensures that the residuals and subsequent one-*SD* standardisation
match the respondents entering that outcome's models.

**Caveat.** This implementation does not separately remove technical variation
associated with genotype batch that is not captured by the harmonised ancestry
coordinates or other upstream genotype processing.

---

## 19. Exploratory post-plan test of parental education as a moderator (S7)

**What.** Birth-year-standardised parental education is tested as a moderator of
the PGI-Education association with educational outcomes, and that moderation is
compared between East and West Germany. The estimates come from the Section B
(B1) right-hand side as published, which already carries
`parental_edu_z_kernel:PGI_Edu_z` and
`PGI_Edu_z:parental_edu_z_kernel:east_west_c`; no new focal model is
specified. A second specification adds `parental_edu_z_kernel:gender_c`, the
only covariate interaction absent from B1 under the Keller (2014)
recommendation that a G × E test carry G × C and E × C for every covariate C.
A third fit substitutes `edu_z_kernel` for `mobility` as the outcome to obtain
the intercept structure needed for the figure. Script:
`04_GxE_Analysis/01_PGI_Analysis/B_Mobility/B_paredu_moderation.R`; it exports
five CSVs — `S7_ParEdu_Keller.csv` (focal G × E terms, both specifications),
`S7_ParEdu_SimpleSlopes.csv` (PGI slope at −1/0/+1 SD parental education),
`S7_ParEdu_Lines.csv` (fitted prediction grid), `ParEdu_Gradient.csv` (the
regional social-origin gradient) and `ParEdu_Gradient_LOO.csv` (its
leave-one-study-out check).

**Plan.** The plan did not specify any test of socioeconomic moderation of the
PGI-Education association. Parental education entered Section B solely as a
control and as a component of the educational mobility outcome.

**Rationale.** Morris et al. (2026, *PNAS*) reported that PGI-Education
associations with educational attainment were disproportionately stronger among
participants from more advantaged origins. Because the terms needed to test the
same question were already present in B1, reporting them is more transparent
than leaving a directly comparable estimate unreported in a manuscript that
cites that study. The analysis is labelled exploratory throughout and is
reported in the Supplementary Results, not the main text.

**Caveat.** Estimated in the educational mobility analytic sample (*N* = 4,904,
no SHIP), in which East German participants are a minority rather than the
majority they form in the educational attainment sample. Neither region's
simple slope is individually distinguishable from zero, so the East–West
difference is a contrast between two undetected slopes. Because parental
education is a model term and educational mobility is the respondent-minus-
parental birth-year-standardised difference, all reported estimates are
identical for educational attainment and educational mobility; only the
parental-education main effect differs between the two parameterisations, by
exactly one.

---

## 20. Analytic sample capped at birth year 1996

**What.** Both analytic samples are restricted to `birth_year <= 1996`
(`04_GxE_Analysis/01_PGI_Analysis/00_setup/prepare_samples.R`). The restriction
is a single study-wide birth-year ceiling applied identically to every
contributing dataset, rather than a participant-specific age computed against
each study's own data-collection date. It is recorded as its own step in the
attrition tables ("Birth year <= 1996") for both the primary random-effects
track and the deduplicated track. No lower birth-year bound is imposed in code;
the lower edge of the sample is whatever the contributing datasets supply.

**Plan.** "We included all individuals born between 1920 and 2000, excluding
those younger than 25 years at data collection, who may not have completed their
education (OECD, 2024)" (`ANALYSIS_PLAN.md`, line 80).

**Rationale.** The 1996 ceiling operationalises the plan's age-25
education-completeness rule as a fixed birth-year cutoff. Birth year is
harmonised across all four studies, whereas an individual-level date of data
collection is not consistently available in the harmonised data, so a fixed
ceiling applies the completeness criterion uniformly and reproducibly across
studies and makes the resulting sample boundary explicit and auditable.

---

## 21. SHIP participants recorded in the West region are excluded

**What.** Rows satisfying `cohort == "SHIP" & east_west == "West"` are dropped
from both analytic samples
(`04_GxE_Analysis/01_PGI_Analysis/00_setup/prepare_samples.R`), applied to the
core sample and again when the full-family random-effects sample is built. SHIP
therefore enters every model as an East-only study, and the count of excluded
rows is reported with the other exclusion counts.

**Plan.** The plan defines the analytic sample by birth year and age at data
collection (`ANALYSIS_PLAN.md`, line 80) and treats region (East/West) as a
moderator throughout; it specifies no study-specific regional restriction.

**Rationale.** SHIP is a regional study recruited entirely in the East, so the
few records carrying a West region code cannot represent genuine West-recruited
participants and are treated as region-coding artefacts. Retaining them would
populate a SHIP × West cell of a handful of individuals that would nonetheless
carry SHIP's entire within-study East–West contrast, and would let study and
region be confounded in region-moderated models.

**Caveat.** SHIP contributes only to the East side of every East–West
comparison, so region remains partly confounded with study composition in the
attainment analyses, where SHIP is a large contributor. SHIP is absent from the
mobility analyses for a separate reason (no parental education is available), so
this exclusion affects the attainment sample only.

---

## 22. Published NonCog weights used in place of study-excluded NonCog

**What.** PGI-Non-Cognitive is computed from the published Malanchini weights
for all five studies. The study-excluded NonCog summary statistics were
constructed but are not used.

**Plan.** "we computed PGI-Non-Cognitive in SHIP and BASE-II from the GWAS
summary statistics in which these cohorts were excluded from the discovery
sample" (`ANALYSIS_PLAN.md`).

**Rationale.** The study-excluded construction did not work. Its
GWAS-by-subtraction over-subtracted cognition, giving an LDSC genetic
correlation with the cognitive GWAS of rg = −0.40 where the construction intends
rg ≈ 0, against rg = +0.16 for the published weights on the same cognitive GWAS.
The factor loadings matched the values implied by the genetic correlation
between the two inputs (rg = 0.60; loadings 0.604 on Cog and 0.797 on NonCog),
so the specification is not the cause; the imbalance between the inputs is
(educational attainment N ≈ 2.4M against N ≈ 224k for cognitive performance).
The published weights are the validated construct and put all five studies on
one posterior.

**Caveat.** The published weights derive from Lee et al. (2018) educational
attainment summary statistics, so overlap with SHIP and BASE-II is not excluded
by construction. PGI-Non-Cognitive is not a focal predictor; it enters the S1
alternative-PGI sensitivity, where it is null in both sections (attainment
*p* = .05 raw, .14 after Benjamini-Hochberg; mobility null). That null supports
the reported contrast attributing the birth-year amplification and the regional
mobility difference primarily to PGI-Cognition.
